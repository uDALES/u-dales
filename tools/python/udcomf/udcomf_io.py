"""Streaming preparation of facet outputs for comfort postprocessing."""

from __future__ import annotations

from contextlib import ExitStack
from dataclasses import dataclass
import json
import os
from pathlib import Path
import tempfile
from typing import Sequence

from netCDF4 import Dataset
import numpy as np

from udconfig import parse_namoptions


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
