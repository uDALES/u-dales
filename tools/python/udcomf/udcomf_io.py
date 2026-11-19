"""Streaming preparation of warm-start facet and statistics outputs for comfort postprocessing."""

from __future__ import annotations

from contextlib import ExitStack
from dataclasses import dataclass
import json
import os
from pathlib import Path
import tempfile
from typing import Sequence

import f90nml
from netCDF4 import Dataset
import numpy as np

from udconfig import parse_namoptions
from .heights import bracketing_levels, configured_receptor_heights
from .schedule import read_stats_schedule


@dataclass(frozen=True)
class FacetMergeResult:
    eb_path: Path
    temperature_path: Path
    source_cases: tuple[str, ...]
    record_count: int
    first_time: float
    last_time: float
    boundary_gaps: tuple[float, ...]


@dataclass(frozen=True)
class _Segment:
    case: str
    eb_path: Path
    temperature_path: Path
    times: np.ndarray


def _source_paths(case_path: Path) -> tuple[str, Path, Path]:
    path = Path(case_path).expanduser().resolve()
    if path.is_dir():
        case = path.name
        eb_path = path / f"facEB.{case}.nc"
    elif path.name.startswith("facEB.") and path.suffix == ".nc":
        case = path.name[len("facEB."):-len(".nc")]
        eb_path = path
    else:
        raise ValueError(f"Expected a case directory or facEB.<case>.nc path: {path}")
    temperature_path = eb_path.with_name(f"facT.{case}.nc")
    if not eb_path.is_file() or not temperature_path.is_file():
        raise FileNotFoundError(f"Paired facEB/facT files are required for {case} in {eb_path.parent}")
    return case, eb_path, temperature_path


def _schema(ds: Dataset) -> tuple:
    return (
        tuple((name, len(dim)) for name, dim in ds.dimensions.items() if name != "time"),
        tuple(
            (name, str(var.datatype), var.dimensions,
             tuple((attr, str(var.getncattr(attr))) for attr in var.ncattrs()))
            for name, var in ds.variables.items()
        ),
    )


def _inspect_segment(case: str, eb_path: Path, temperature_path: Path,
                     nfcts: int, nfaclyrs: int) -> _Segment:
    with Dataset(eb_path) as eb, Dataset(temperature_path) as temperature:
        if ("fct" not in eb.dimensions or "fct" not in temperature.dimensions
                or len(eb.dimensions["fct"]) != nfcts
                or len(temperature.dimensions["fct"]) != nfcts
                or "lyr" not in temperature.dimensions
                or len(temperature.dimensions["lyr"]) != nfaclyrs + 1):
            raise ValueError(f"Facet or layer dimensions do not match the destination case: {case}")
        required = ((eb, ("t", "LWout", "netsw")),
                    (temperature, ("t", "T", "dTdz")))
        for ds, names in required:
            if "time" not in ds.dimensions:
                raise ValueError(f"Missing time dimension in {ds.filepath()}")
            for name in names:
                if name not in ds.variables or ds.variables[name].dimensions[0] != "time":
                    raise ValueError(f"Missing or malformed {name} in {ds.filepath()}")
        times = np.asarray(eb.variables["t"][:], dtype=float)
        temperature_times = np.asarray(temperature.variables["t"][:], dtype=float)
        if (times.ndim != 1 or times.size == 0 or not np.isfinite(times).all()
                or np.any(np.diff(times) <= 0) or not np.array_equal(times, temperature_times)):
            raise ValueError(f"facEB and facT timestamps must match and increase in case {case}")
    return _Segment(case, eb_path, temperature_path, times)


def _create_output(source: Dataset, path: Path, source_cases: tuple[str, ...],
                   boundary_gaps: tuple[float, ...]) -> Dataset:
    output = Dataset(path, "w", format=source.data_model)
    try:
        for name, dim in source.dimensions.items():
            output.createDimension(name, None if name == "time" else len(dim))
        for name, var in source.variables.items():
            options = {}
            if "_FillValue" in var.ncattrs():
                options["fill_value"] = var.getncattr("_FillValue")
            chunks = var.chunking()
            if isinstance(chunks, list):
                options["chunksizes"] = tuple(chunks)
            filters = var.filters()
            for key in ("zlib", "complevel", "shuffle", "fletcher32"):
                if key in filters:
                    options[key] = filters[key]
            result = output.createVariable(name, var.datatype, var.dimensions, **options)
            result.setncatts({key: var.getncattr(key) for key in var.ncattrs()
                              if key != "_FillValue"})
            if "time" not in var.dimensions:
                result[:] = var[:]
        output.setncatts({key: source.getncattr(key) for key in source.ncattrs()})
        output.setncattr("udcomf_merged_cases", json.dumps(source_cases))
        output.setncattr("udcomf_boundary_gaps_seconds", json.dumps(boundary_gaps))
        return output
    except BaseException:
        output.close()
        raise


def merge_facet_outputs(
    destination_case: Path, case_paths: Sequence[Path], *, overwrite: bool = False,
) -> FacetMergeResult:
    """Merge paired warm-start facEB/facT files without filling time gaps.

    Each source is a case directory or a copied ``facEB.<case>.nc`` file. The
    output is ``facEB.<destination>.nc`` and ``facT.<destination>.nc`` in the
    destination directory. Records are ordered by their saved simulation time;
    overlapping or duplicate timestamps are rejected, and sources are never
    modified. The merge copies one time record at a time to bound memory use.
    """
    destination = Path(destination_case).expanduser().resolve()
    if isinstance(case_paths, (str, Path)) or not case_paths:
        raise ValueError("case_paths must be a nonempty sequence of paths")
    options = parse_namoptions(destination / f"namoptions.{destination.name}")
    nfcts = int(options["nfcts"])
    nfaclyrs = int(options["nfaclyrs"])
    segments = []
    seen = set()
    for path in case_paths:
        case, eb_path, temperature_path = _source_paths(Path(path))
        if eb_path in seen:
            raise ValueError(f"Duplicate warm-start source: {eb_path}")
        seen.add(eb_path)
        segments.append(_inspect_segment(case, eb_path, temperature_path, nfcts, nfaclyrs))
    segments.sort(key=lambda segment: segment.times[0])
    gaps = []
    for previous, current in zip(segments, segments[1:]):
        gap = float(current.times[0] - previous.times[-1])
        if gap <= 0:
            raise ValueError(f"Warm-start records overlap: {previous.case} and {current.case}")
        gaps.append(gap)

    for kind in ("eb_path", "temperature_path"):
        reference = None
        for segment in segments:
            with Dataset(getattr(segment, kind)) as ds:
                schema = _schema(ds)
            if reference is None:
                reference = schema
            elif schema != reference:
                raise ValueError(f"Incompatible {kind} NetCDF schema in case {segment.case}")

    eb_output = destination / f"facEB.{destination.name}.nc"
    temperature_output = destination / f"facT.{destination.name}.nc"
    if any(path in (eb_output, temperature_output) for segment in segments
           for path in (segment.eb_path, segment.temperature_path)):
        raise ValueError("Destination files cannot also be warm-start sources")
    if not overwrite and (eb_output.exists() or temperature_output.exists()):
        raise FileExistsError("Merged facet output already exists; use overwrite=True to replace it")

    temp_paths = []
    try:
        for output in (eb_output, temperature_output):
            fd, name = tempfile.mkstemp(prefix=f".{output.name}.", suffix=".tmp", dir=destination)
            os.close(fd)
            temp_paths.append(Path(name))
        cases = tuple(segment.case for segment in segments)
        with ExitStack() as stack:
            first_eb = stack.enter_context(Dataset(segments[0].eb_path))
            first_temperature = stack.enter_context(Dataset(segments[0].temperature_path))
            first_eb.set_auto_mask(False)
            first_temperature.set_auto_mask(False)
            out_eb = stack.enter_context(_create_output(first_eb, temp_paths[0], cases, tuple(gaps)))
            out_temperature = stack.enter_context(
                _create_output(first_temperature, temp_paths[1], cases, tuple(gaps))
            )
            offset = 0
            for segment in segments:
                with Dataset(segment.eb_path) as eb, Dataset(segment.temperature_path) as temperature:
                    eb.set_auto_mask(False)
                    temperature.set_auto_mask(False)
                    for source, target in ((eb, out_eb), (temperature, out_temperature)):
                        for name, variable in source.variables.items():
                            if "time" not in variable.dimensions:
                                continue
                            if variable.dimensions[0] != "time":
                                raise ValueError(f"Time must be the first dimension of {name}")
                            for i in range(segment.times.size):
                                target.variables[name][offset + i, ...] = variable[i, ...]
                offset += segment.times.size
        if not overwrite and (eb_output.exists() or temperature_output.exists()):
            raise FileExistsError("Destination appeared while writing merged facet output")
        os.replace(temp_paths[0], eb_output)
        os.replace(temp_paths[1], temperature_output)
    finally:
        for path in temp_paths:
            path.unlink(missing_ok=True)

    return FacetMergeResult(
        eb_path=eb_output,
        temperature_path=temperature_output,
        source_cases=cases,
        record_count=sum(len(segment.times) for segment in segments),
        first_time=float(segments[0].times[0]),
        last_time=float(segments[-1].times[-1]),
        boundary_gaps=tuple(gaps),
    )


# Statistics merge: only what the comfort pipeline reads, at the scalar levels
# bracketing the configured receptor heights.
_STATS_COORDS = ("time", "zt", "yt", "xt", "receptor_height")
_STATS_LEVEL_FIELDS = ("tha", "pabs", "qt", "ql")
_STATS_PLANE_FIELDS = ("ws_local", "ws_10")
_STATS_NAMES = ("stats_tkslice", "stats_t")


@dataclass(frozen=True)
class StatisticsMergeResult:
    path: Path
    source_cases: tuple[str, ...]
    source_files: tuple[Path, ...]
    record_count: int
    first_time: float
    last_time: float
    boundary_gaps: tuple[float, ...]
    selected_levels: tuple[int, ...]


@dataclass(frozen=True)
class _StatsSegment:
    case: str
    path: Path
    times: np.ndarray
    level_indices: tuple[int, ...]


def _statistics_source(case_path: Path) -> tuple[str, Path]:
    path = Path(case_path).expanduser().resolve()
    if path.is_dir():
        case = path.name
        for name in _STATS_NAMES:
            candidate = path / f"{name}.{case}.nc"
            if candidate.is_file():
                return case, candidate
        raise FileNotFoundError(f"No gathered stats_tkslice or stats_t file for case {case} in {path}")
    if path.suffix == ".nc":
        for name in _STATS_NAMES:
            prefix = f"{name}."
            if path.name.startswith(prefix):
                case = path.name[len(prefix):-len(".nc")]
                if not path.is_file():
                    raise FileNotFoundError(f"Statistics file is missing: {path}")
                return case, path
    raise ValueError(f"Expected a case directory or stats_t[kslice].<case>.nc path: {path}")


def _inspect_statistics_segment(case: str, path: Path, z_needed: np.ndarray,
                                tstatsdump: float) -> _StatsSegment:
    options_path = path.parent / f"namoptions.{case}"
    if options_path.is_file():
        segment_dump = float(f90nml.read(options_path).get("output", {}).get("tstatsdump", 10000.0))
        if not np.isclose(segment_dump, tstatsdump, rtol=0, atol=1e-6):
            raise ValueError(
                f"tstatsdump differs between the destination ({tstatsdump:g} s) and case {case} ({segment_dump:g} s)"
            )
    with Dataset(path) as ds:
        for name in ("time", "zt", "yt", "xt"):
            if name not in ds.dimensions:
                raise ValueError(f"Missing dimension {name} in {path.name}")
        for name in _STATS_COORDS + _STATS_LEVEL_FIELDS + _STATS_PLANE_FIELDS:
            if name not in ds.variables:
                raise ValueError(f"Missing variable {name} in {path.name}")
        for name in _STATS_LEVEL_FIELDS + _STATS_PLANE_FIELDS:
            dims = ds.variables[name].dimensions
            if not dims or dims[0] != "time":
                raise ValueError(f"Time must be the first dimension of {name} in {path.name}")
        if "zt" not in ds.variables["tha"].dimensions:
            raise ValueError(f"tha must carry the zt dimension in {path.name}")
        times = np.asarray(ds.variables["time"][:], dtype=float)
        if times.ndim != 1 or times.size == 0 or not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
            raise ValueError(f"Statistics times must be finite and increasing in {path.name}")
        z_saved = np.asarray(ds.variables["zt"][:], dtype=float)
        indices = []
        for z in z_needed:
            matches = np.flatnonzero(np.isclose(z_saved, z, rtol=1e-6, atol=1e-4))
            if len(matches) != 1:
                raise ValueError(f"{path.name} lacks the scalar level at z = {z:g} m")
            indices.append(int(matches[0]))
    return _StatsSegment(case, path, times, tuple(indices))


def _same_coordinates(segments: Sequence[_StatsSegment]) -> None:
    reference = None
    for segment in segments:
        with Dataset(segment.path) as ds:
            ds.set_auto_mask(False)
            coords = {
                name: np.asarray(ds.variables[name][:], dtype=float)
                for name in ("xt", "yt", "receptor_height")
            }
            coords["zt"] = np.asarray(ds.variables["zt"][:], dtype=float)[list(segment.level_indices)]
            layout = tuple((name, ds.variables[name].dimensions) for name in _STATS_LEVEL_FIELDS + _STATS_PLANE_FIELDS)
        if reference is None:
            reference = (coords, layout)
            continue
        if layout != reference[1]:
            raise ValueError(f"Variable layout differs in case {segment.case}")
        for name, value in coords.items():
            if value.shape != reference[0][name].shape or not np.allclose(value, reference[0][name], rtol=0, atol=1e-4):
                raise ValueError(f"Coordinate {name} differs in case {segment.case}")


def merge_statistics_outputs(
    destination_case: Path, case_paths: Sequence[Path], *, overwrite: bool = False,
) -> StatisticsMergeResult:
    """Merge the comfort subset of warm-start statistics into one stats_tkslice file.

    Each source is a case directory (its gathered ``stats_tkslice`` or
    ``stats_t`` file is used) or such a file. Only ``tha``, ``pabs``, ``qt``
    and ``ql`` at the scalar levels bracketing the destination's receptor
    heights, plus ``ws_local`` and ``ws_10``, are kept. Records are ordered by
    saved time; overlapping averaging windows are rejected and gaps are
    preserved. The output is ``stats_tkslice.<destination>.nc``.
    """
    destination = Path(destination_case).expanduser().resolve()
    if isinstance(case_paths, (str, Path)) or not case_paths:
        raise ValueError("case_paths must be a nonempty sequence of paths")
    heights = configured_receptor_heights(destination)
    schedule = read_stats_schedule(destination)
    profile = np.loadtxt(destination / f"prof.inp.{destination.name}", comments="#", ndmin=2)
    z_model = np.asarray(profile[:, 0], dtype=float)
    brackets = bracketing_levels(z_model, heights)
    levels = tuple(sorted({k for lo, hi, _ in brackets.values() for k in (lo, hi)}))
    z_needed = z_model[list(levels)]

    segments = []
    seen = set()
    for path in case_paths:
        case, stats_path = _statistics_source(Path(path))
        if stats_path in seen:
            raise ValueError(f"Duplicate warm-start source: {stats_path}")
        seen.add(stats_path)
        segments.append(_inspect_statistics_segment(case, stats_path, z_needed, schedule.tstatsdump))
    segments.sort(key=lambda segment: segment.times[0])
    gaps = []
    for previous, current in zip(segments, segments[1:]):
        gap = float(current.times[0] - schedule.tstatsdump - previous.times[-1])
        if gap < -schedule.tsample:
            raise ValueError(
                f"Averaging windows overlap between cases {previous.case} and {current.case}"
            )
        gaps.append(gap)
    _same_coordinates(segments)

    output = destination / f"stats_tkslice.{destination.name}.nc"
    if any(segment.path == output for segment in segments):
        raise ValueError("The destination stats_tkslice file cannot also be a warm-start source")
    if not overwrite and output.exists():
        raise FileExistsError("Merged statistics output already exists; use overwrite=True to replace it")

    fd, name = tempfile.mkstemp(prefix=f".{output.name}.", suffix=".tmp", dir=destination)
    os.close(fd)
    temp_path = Path(name)
    try:
        with Dataset(segments[0].path) as first:
            first.set_auto_mask(False)
            with Dataset(temp_path, "w", format="NETCDF4") as out:
                for dim_name, dim in first.dimensions.items():
                    if dim_name == "time":
                        out.createDimension(dim_name, None)
                    elif dim_name == "zt":
                        out.createDimension(dim_name, len(levels))
                    elif dim_name in ("xt", "yt", "receptor_height"):
                        out.createDimension(dim_name, len(dim))
                for var_name in _STATS_COORDS + _STATS_LEVEL_FIELDS + _STATS_PLANE_FIELDS:
                    var = first.variables[var_name]
                    options = {}
                    if "_FillValue" in var.ncattrs():
                        options["fill_value"] = var.getncattr("_FillValue")
                    if "time" in var.dimensions:
                        options.update(zlib=True, complevel=3)
                    result = out.createVariable(var_name, var.datatype, var.dimensions, **options)
                    result.setncatts({key: var.getncattr(key) for key in var.ncattrs() if key != "_FillValue"})
                    if "time" in var.dimensions:
                        continue
                    if var_name == "zt":
                        result[:] = np.asarray(var[:])[list(segments[0].level_indices)]
                    else:
                        result[...] = var[...]
                out.setncatts({key: first.getncattr(key) for key in first.ncattrs()})
                out.setncattr("udcomf_merged_cases", json.dumps([s.case for s in segments]))
                out.setncattr("udcomf_source_files", json.dumps([s.path.name for s in segments]))
                out.setncattr("udcomf_segment_first_times", json.dumps([float(s.times[0]) for s in segments]))
                out.setncattr("udcomf_segment_last_times", json.dumps([float(s.times[-1]) for s in segments]))
                out.setncattr("udcomf_boundary_gaps_seconds", json.dumps(gaps))
                out.setncattr("udcomf_selected_scalar_levels", json.dumps([k + 1 for k in levels]))
                out.setncattr("udcomf_averaging_window_s", float(schedule.tstatsdump))
                offset = 0
                for segment in segments:
                    with Dataset(segment.path) as source:
                        source.set_auto_mask(False)
                        for var_name in ("time",) + _STATS_LEVEL_FIELDS + _STATS_PLANE_FIELDS:
                            var = source.variables[var_name]
                            dims = var.dimensions
                            selector = [slice(None)] * len(dims)
                            if "zt" in dims:
                                selector[dims.index("zt")] = list(segment.level_indices)
                            for i in range(segment.times.size):
                                selector[0] = i
                                out.variables[var_name][offset + i, ...] = var[tuple(selector)]
                    offset += segment.times.size
        if not overwrite and output.exists():
            raise FileExistsError("Destination appeared while writing merged statistics")
        os.replace(temp_path, output)
    finally:
        temp_path.unlink(missing_ok=True)

    return StatisticsMergeResult(
        path=output,
        source_cases=tuple(segment.case for segment in segments),
        source_files=tuple(segment.path for segment in segments),
        record_count=sum(segment.times.size for segment in segments),
        first_time=float(segments[0].times[0]),
        last_time=float(segments[-1].times[-1]),
        boundary_gaps=tuple(gaps),
        selected_levels=levels,
    )
