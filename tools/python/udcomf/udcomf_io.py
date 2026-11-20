"""Streaming preparation of warm-start facet and statistics outputs for comfort postprocessing."""

from __future__ import annotations

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
from .sources import MERGED_FACETS, MERGED_STATISTICS


@dataclass(frozen=True)
class FacetMergeResult:
    path: Path
    source_cases: tuple[str, ...]
    source_files: tuple[Path, ...]
    record_count: int
    first_time: float
    last_time: float
    boundary_gaps: tuple[float, ...]


@dataclass(frozen=True)
class _Segment:
    case: str
    eb_path: Path
    times: np.ndarray


# Facet merge: only what the pedestrian longwave stage reads from facEB.
_FACET_FIELDS = ("t", "LWout")


def _source_paths(case_path: Path) -> tuple[str, Path]:
    path = Path(case_path).expanduser().resolve()
    if path.is_dir():
        case = path.name
        eb_path = path / f"facEB.{case}.nc"
    elif path.name.startswith("facEB.") and path.suffix == ".nc":
        case = path.name[len("facEB."):-len(".nc")]
        eb_path = path
    else:
        raise ValueError(f"Expected a case directory or facEB.<case>.nc path: {path}")
    if not eb_path.is_file():
        raise FileNotFoundError(f"facEB file is required for {case} in {eb_path.parent}")
    return case, eb_path


def _facet_schema(ds: Dataset) -> tuple:
    return tuple(
        (name, str(ds.variables[name].datatype), ds.variables[name].dimensions)
        for name in _FACET_FIELDS
    )


def _inspect_segment(case: str, eb_path: Path, nfcts: int) -> _Segment:
    with Dataset(eb_path) as eb:
        if "fct" not in eb.dimensions or len(eb.dimensions["fct"]) != nfcts:
            raise ValueError(f"Facet dimension does not match the destination case: {case}")
        if "time" not in eb.dimensions:
            raise ValueError(f"Missing time dimension in {eb_path}")
        for name in _FACET_FIELDS:
            if name not in eb.variables or eb.variables[name].dimensions[0] != "time":
                raise ValueError(f"Missing or malformed {name} in {eb_path}")
        if eb.variables["LWout"].dimensions != ("time", "fct"):
            raise ValueError(f"LWout must be dimensioned (time, fct) in {eb_path}")
        times = np.asarray(eb.variables["t"][:], dtype=float)
        if (times.ndim != 1 or times.size == 0 or not np.isfinite(times).all()
                or np.any(np.diff(times) <= 0)):
            raise ValueError(f"facEB timestamps must be finite and increasing in case {case}")
    return _Segment(case, eb_path, times)


def merge_facet_outputs(
    destination_case: Path, case_paths: Sequence[Path],
) -> FacetMergeResult:
    """Merge the longwave source of warm-start facEB files into one file.

    Only the facet timestamps ``t`` and the emitted longwave ``LWout`` are
    kept: these are the only facEB quantities the comfort indices need (MRT,
    and through it PET and UTCI, and the WBGT sensor balances). facT is not
    read by any index and is not merged. Each source is a case directory or a
    ``facEB.<case>.nc`` file; the destination's own facEB may be a source. The
    output is ``merged_facEB.<destination>.nc`` in the destination directory,
    replaced on every call because it is a cheap deterministic function of the
    sources, which are never modified. Records are ordered by saved simulation
    time, overlapping or duplicate timestamps are rejected, and gaps are
    preserved. One time record is copied at a time to bound memory use.
    """
    destination = Path(destination_case).expanduser().resolve()
    if isinstance(case_paths, (str, Path)) or not case_paths:
        raise ValueError("case_paths must be a nonempty sequence of paths")
    options = parse_namoptions(destination / f"namoptions.{destination.name}")
    nfcts = int(options["nfcts"])
    segments = []
    seen = set()
    for path in case_paths:
        case, eb_path = _source_paths(Path(path))
        if eb_path in seen:
            raise ValueError(f"Duplicate warm-start source: {eb_path}")
        seen.add(eb_path)
        segments.append(_inspect_segment(case, eb_path, nfcts))
    segments.sort(key=lambda segment: segment.times[0])
    gaps = []
    for previous, current in zip(segments, segments[1:]):
        gap = float(current.times[0] - previous.times[-1])
        if gap <= 0:
            raise ValueError(f"Warm-start records overlap: {previous.case} and {current.case}")
        gaps.append(gap)
    reference = None
    for segment in segments:
        with Dataset(segment.eb_path) as ds:
            schema = _facet_schema(ds)
        if reference is None:
            reference = schema
        elif schema != reference:
            raise ValueError(f"Incompatible facEB variable schema in case {segment.case}")

    output = destination / f"{MERGED_FACETS}.{destination.name}.nc"
    if any(segment.eb_path == output for segment in segments):
        raise ValueError("The merged facet file cannot also be a warm-start source")

    fd, name = tempfile.mkstemp(prefix=f".{output.name}.", suffix=".tmp", dir=destination)
    os.close(fd)
    temp_path = Path(name)
    try:
        cases = tuple(segment.case for segment in segments)
        with Dataset(segments[0].eb_path) as first, Dataset(temp_path, "w", format="NETCDF4") as out:
            first.set_auto_mask(False)
            out.createDimension("time", None)
            out.createDimension("fct", nfcts)
            for var_name in _FACET_FIELDS:
                var = first.variables[var_name]
                options = {}
                if "_FillValue" in var.ncattrs():
                    options["fill_value"] = var.getncattr("_FillValue")
                if var_name == "LWout":
                    options.update(zlib=True, complevel=3)
                result = out.createVariable(var_name, var.datatype, var.dimensions, **options)
                result.setncatts({key: var.getncattr(key) for key in var.ncattrs()
                                  if key != "_FillValue"})
            out.setncattr("udcomf_merged_cases", json.dumps(cases))
            out.setncattr("udcomf_source_files", json.dumps([s.eb_path.name for s in segments]))
            out.setncattr("udcomf_segment_first_times", json.dumps([float(s.times[0]) for s in segments]))
            out.setncattr("udcomf_segment_last_times", json.dumps([float(s.times[-1]) for s in segments]))
            out.setncattr("udcomf_boundary_gaps_seconds", json.dumps(gaps))
            out.setncattr("udcomf_contents", "facet timestamps t and emitted longwave LWout only")
            offset = 0
            for segment in segments:
                with Dataset(segment.eb_path) as eb:
                    eb.set_auto_mask(False)
                    for var_name in _FACET_FIELDS:
                        variable = eb.variables[var_name]
                        for i in range(segment.times.size):
                            out.variables[var_name][offset + i, ...] = variable[i, ...]
                offset += segment.times.size
        os.replace(temp_path, output)
    finally:
        temp_path.unlink(missing_ok=True)

    return FacetMergeResult(
        path=output,
        source_cases=cases,
        source_files=tuple(segment.eb_path for segment in segments),
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
_STATS_NAMES = ("stats_t", "stats_tkslice")


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


def _statistics_source(case_path: Path, source: str) -> tuple[str, Path]:
    """The ``source.<case>.nc`` file of a segment; no other file is tried."""
    path = Path(case_path).expanduser().resolve()
    prefix = f"{source}."
    if path.is_dir():
        case = path.name
        candidate = path / f"{source}.{case}.nc"
        if not candidate.is_file():
            raise FileNotFoundError(f"Gathered {source} file is missing for case {case}: {candidate}")
        return case, candidate
    if path.suffix == ".nc" and path.name.startswith(prefix):
        case = path.name[len(prefix):-len(".nc")]
        if not path.is_file():
            raise FileNotFoundError(f"Statistics file is missing: {path}")
        return case, path
    raise ValueError(f"Expected a case directory or {source}.<case>.nc path: {path}")


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
    destination_case: Path, case_paths: Sequence[Path], *, source: str = "stats_t",
) -> StatisticsMergeResult:
    """Merge the comfort subset of warm-start statistics into one file.

    ``source`` names the gathered file read from every segment, ``stats_t``
    (default) or ``stats_tkslice``; it must exist in each segment and no
    other file is tried. Each source path is a case directory or the
    ``source.<case>.nc`` file itself; the destination's own file may be a
    segment. Only ``tha``, ``pabs``, ``qt`` and ``ql`` at the scalar levels
    bracketing the destination's receptor heights, plus ``ws_local`` and
    ``ws_10``, are kept. Records are ordered by saved time; overlapping
    averaging windows are rejected and gaps are preserved. The output is
    always ``merged_stats_tkslice.<destination>.nc``, replaced on every call.
    """
    if source not in _STATS_NAMES:
        raise ValueError(f"source must be one of {list(_STATS_NAMES)}")
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
        case, stats_path = _statistics_source(Path(path), source)
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

    output = destination / f"{MERGED_STATISTICS}.{destination.name}.nc"
    if any(segment.path == output for segment in segments):
        raise ValueError("The merged statistics file cannot also be a warm-start source")

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
