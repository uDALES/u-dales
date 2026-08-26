#!/usr/bin/env python3
"""Strict output comparison helpers for CPU/GPU uDALES parity tests.

NetCDF outputs, and the Fortran unformatted restart files initd/inits."""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple

import netCDF4
import numpy as np


@dataclasses.dataclass(frozen=True)
class Tolerance:
    atol: float
    rtol: float


@dataclasses.dataclass
class ComparisonResult:
    passed: bool
    files_compared: int
    variables_compared: int
    failures: list[str]
    variable_results: list[Dict[str, Any]]

    def as_dict(self) -> Dict[str, Any]:
        return dataclasses.asdict(self)


def _tolerance_for(
    file_name: str,
    variable_name: str,
    dtype: np.dtype[Any],
    tolerance_config: Mapping[str, Any],
) -> Tolerance:
    if np.issubdtype(dtype, np.integer) or np.issubdtype(dtype, np.bool_):
        default = Tolerance(0.0, 0.0)
    else:
        raw_default = tolerance_config.get("default", {})
        default = Tolerance(
            float(raw_default.get("atol", 1.0e-6)),
            float(raw_default.get("rtol", 1.0e-6)),
        )

    overrides = tolerance_config.get("variables", {})
    raw_override = overrides.get(f"{file_name}:{variable_name}", overrides.get(variable_name))
    if raw_override is None:
        return default
    return Tolerance(
        float(raw_override.get("atol", default.atol)),
        float(raw_override.get("rtol", default.rtol)),
    )


# The restart files are Fortran sequential unformatted, not NetCDF, and their
# layout is fixed by modsave::writerestartfiles: one record per array, in this
# order, with (timee, dt) closing the main file and timee closing the scalar
# one. A record's element count follows from its byte length, so the domain
# size need not be known here; what must be known is the element type, and
# only wall is an integer. The record names are what the failure messages and
# the per-variable tolerance overrides use.
RESTART_PREFIXES = ("initd", "inits")
_RESTART_SCHEMA: Dict[str, Sequence[Tuple[str, np.dtype[Any]]]] = {
    "initd": (
        ("mindist", np.dtype("<f8")),
        ("wall", np.dtype("<i4")),
        ("u0", np.dtype("<f8")),
        ("v0", np.dtype("<f8")),
        ("w0", np.dtype("<f8")),
        ("pres0", np.dtype("<f8")),
        ("thl0", np.dtype("<f8")),
        ("e120", np.dtype("<f8")),
        ("ekm", np.dtype("<f8")),
        ("qt0", np.dtype("<f8")),
        ("ql0", np.dtype("<f8")),
        ("ql0h", np.dtype("<f8")),
        ("time", np.dtype("<f8")),
    ),
    "inits": (
        ("sv0", np.dtype("<f8")),
        ("time", np.dtype("<f8")),
    ),
}


def _read_fortran_records(path: Path) -> list[bytes]:
    """Split a sequential unformatted file into its records.

    Both compilers this project builds with, gfortran and nvfortran, frame each
    record with a 4-byte little-endian length before and after the payload.
    The trailing marker is checked against the leading one, so a file written
    with 8-byte markers, or truncated, is reported rather than misread.
    """
    data = path.read_bytes()
    records: list[bytes] = []
    offset = 0
    while offset < len(data):
        if offset + 4 > len(data):
            raise ValueError(f"{path.name}: truncated record marker at byte {offset}")
        length = int(np.frombuffer(data, dtype="<i4", count=1, offset=offset)[0])
        start = offset + 4
        end = start + length
        if length < 0 or end + 4 > len(data):
            raise ValueError(
                f"{path.name}: record at byte {offset} claims {length} bytes, "
                f"which runs past the end of the file"
            )
        trailer = int(np.frombuffer(data, dtype="<i4", count=1, offset=end)[0])
        if trailer != length:
            raise ValueError(
                f"{path.name}: record at byte {offset} has leading marker {length} "
                f"but trailing marker {trailer}"
            )
        records.append(data[start:end])
        offset = end + 4
    return records


def _restart_records(path: Path) -> Dict[str, np.ndarray[Any, Any]]:
    """Name and type each record of a restart file by its schema."""
    schema = _RESTART_SCHEMA[path.name[:5]]
    raw = _read_fortran_records(path)
    if len(raw) != len(schema):
        raise ValueError(
            f"{path.name}: {len(raw)} records, but the writer's layout has {len(schema)}"
        )
    decoded: Dict[str, np.ndarray[Any, Any]] = {}
    for (name, dtype), payload in zip(schema, raw):
        if len(payload) % dtype.itemsize != 0:
            raise ValueError(
                f"{path.name}:{name}: {len(payload)} bytes is not a whole number "
                f"of {dtype.itemsize}-byte elements"
            )
        decoded[name] = np.frombuffer(payload, dtype=dtype)
    return decoded


def _compare_restart_file(
    reference_path: Path,
    candidate_path: Path,
    tolerance_config: Mapping[str, Any],
) -> Tuple[list[str], list[Dict[str, Any]]]:
    """Compare two restart files record by record, under the NetCDF rules.

    Integer records must match exactly; real ones use the same tolerance
    lookup as a NetCDF variable of the same name, so the (timee, dt) record,
    named time, gets the tight time tolerance the matrix already sets.
    """
    failures: list[str] = []
    results: list[Dict[str, Any]] = []
    file_name = reference_path.name

    try:
        reference = _restart_records(reference_path)
        candidate = _restart_records(candidate_path)
    except ValueError as error:
        return [f"{file_name}: {error}"], results

    for name in reference:
        label = f"{file_name}:{name}"
        ref_values = reference[name]
        cand_values = candidate[name]
        if ref_values.shape != cand_values.shape:
            failures.append(
                f"{label}: record length differs: CPU={ref_values.shape}, GPU={cand_values.shape}"
            )
            continue
        valid = np.ones(ref_values.shape, dtype=bool)
        tolerance = _tolerance_for(file_name, name, ref_values.dtype, tolerance_config)
        passed, detail = _compare_numeric_variable(ref_values, cand_values, valid, tolerance)
        result = {"file": file_name, "variable": name, "passed": passed, **detail}
        results.append(result)
        if not passed:
            failures.append(f"{label}: {json.dumps(detail, sort_keys=True)}")

    return failures, results


def _output_files(directory: Path, prefix: str, case_id: int) -> Dict[str, Path]:
    """Files a required_outputs prefix names in one run directory.

    An exact ``<prefix>.<case>.nc`` wins over the glob when it exists. The glob
    is what lets ``fielddump`` find ``fielddump.000.000.103.nc``, but on its own
    it makes ``fac`` match ``facEB`` and ``facT`` as well, so a matrix entry
    listing all three would compare the two energy-balance files twice under
    two names. The output-file-set check would then also report a difference
    for ``fac`` whenever only one build wrote a ``facT`` file, blaming the
    wrong output.
    """
    if prefix in RESTART_PREFIXES:
        # initd00000123_000_000.103: the step count and the rank pair sit
        # between prefix and case number, and there is no .nc.
        return {
            path.name: path
            for path in sorted(directory.glob(f"{prefix}*.{case_id:03d}"))
            if path.is_file()
        }
    suffix = f".{case_id:03d}.nc"
    exact = directory / f"{prefix}{suffix}"
    if exact.is_file():
        return {exact.name: exact}
    return {
        path.name: path
        for path in sorted(directory.glob(f"{prefix}*{suffix}"))
        if path.is_file()
    }


def _masked_values(value: Any) -> Tuple[np.ndarray[Any, Any], np.ndarray[Any, Any]]:
    masked = np.ma.asarray(value)
    return np.asarray(masked.data), np.ma.getmaskarray(masked)


def _location(index: int, shape: Sequence[int]) -> Tuple[int, ...]:
    if not shape:
        return ()
    return tuple(int(item) for item in np.unravel_index(index, tuple(shape)))


def _compare_numeric_variable(
    reference: np.ndarray[Any, Any],
    candidate: np.ndarray[Any, Any],
    valid: np.ndarray[Any, Any],
    tolerance: Tolerance,
) -> Tuple[bool, Dict[str, Any]]:
    if not np.any(valid):
        return True, {
            "max_abs_error": 0.0,
            "max_rel_error": 0.0,
            "location": None,
        }

    ref_values = reference[valid]
    candidate_values = candidate[valid]
    finite = np.isfinite(ref_values) & np.isfinite(candidate_values)
    if not np.all(finite):
        invalid_flat = int(np.flatnonzero(valid)[int(np.flatnonzero(~finite)[0])])
        return False, {
            "reason": "unmasked NaN or infinity",
            "location": _location(invalid_flat, reference.shape),
        }

    difference = np.abs(candidate_values - ref_values)
    allowed = tolerance.atol + tolerance.rtol * np.abs(ref_values)
    passed_values = difference <= allowed
    max_position = int(np.argmax(difference))
    source_flat_indices = np.flatnonzero(valid)
    source_flat_index = int(source_flat_indices[max_position])
    denominator = np.maximum(np.abs(ref_values), np.finfo(np.float64).tiny)
    relative = difference / denominator

    return bool(np.all(passed_values)), {
        "max_abs_error": float(difference[max_position]),
        "max_rel_error": float(np.max(relative)),
        "location": _location(source_flat_index, reference.shape),
        "reference": float(ref_values[max_position]),
        "candidate": float(candidate_values[max_position]),
        "atol": tolerance.atol,
        "rtol": tolerance.rtol,
    }


def _compare_file(
    reference_path: Path,
    candidate_path: Path,
    tolerance_config: Mapping[str, Any],
) -> Tuple[list[str], list[Dict[str, Any]]]:
    failures: list[str] = []
    results: list[Dict[str, Any]] = []
    file_name = reference_path.name

    with netCDF4.Dataset(reference_path) as reference_ds, netCDF4.Dataset(candidate_path) as candidate_ds:
        reference_dimensions = {
            name: (len(dimension), dimension.isunlimited())
            for name, dimension in reference_ds.dimensions.items()
        }
        candidate_dimensions = {
            name: (len(dimension), dimension.isunlimited())
            for name, dimension in candidate_ds.dimensions.items()
        }
        if reference_dimensions != candidate_dimensions:
            failures.append(
                f"{file_name}: dimension mismatch: "
                f"CPU={reference_dimensions}, GPU={candidate_dimensions}"
            )

        reference_names = set(reference_ds.variables)
        candidate_names = set(candidate_ds.variables)
        if reference_names != candidate_names:
            missing_gpu = sorted(reference_names - candidate_names)
            missing_cpu = sorted(candidate_names - reference_names)
            if missing_gpu:
                failures.append(f"{file_name}: variables missing from GPU output: {missing_gpu}")
            if missing_cpu:
                failures.append(f"{file_name}: variables missing from CPU output: {missing_cpu}")

        for variable_name in sorted(reference_names & candidate_names):
            reference_variable = reference_ds.variables[variable_name]
            candidate_variable = candidate_ds.variables[variable_name]
            label = f"{file_name}:{variable_name}"

            if reference_variable.dimensions != candidate_variable.dimensions:
                failures.append(
                    f"{label}: dimension-name mismatch: "
                    f"CPU={reference_variable.dimensions}, GPU={candidate_variable.dimensions}"
                )
                continue
            if reference_variable.shape != candidate_variable.shape:
                failures.append(
                    f"{label}: shape mismatch: "
                    f"CPU={reference_variable.shape}, GPU={candidate_variable.shape}"
                )
                continue
            if reference_variable.dtype != candidate_variable.dtype:
                failures.append(
                    f"{label}: dtype mismatch: "
                    f"CPU={reference_variable.dtype}, GPU={candidate_variable.dtype}"
                )
                continue

            reference, reference_mask = _masked_values(reference_variable[:])
            candidate, candidate_mask = _masked_values(candidate_variable[:])
            if not np.array_equal(reference_mask, candidate_mask):
                failures.append(f"{label}: missing/fill-value mask differs")
                continue

            valid = ~reference_mask
            if np.issubdtype(reference.dtype, np.number):
                tolerance = _tolerance_for(
                    file_name,
                    variable_name,
                    reference.dtype,
                    tolerance_config,
                )
                passed, detail = _compare_numeric_variable(
                    reference,
                    candidate,
                    valid,
                    tolerance,
                )
            else:
                passed = bool(np.array_equal(reference[valid], candidate[valid]))
                detail = {} if passed else {"reason": "non-numeric values differ"}

            result = {
                "file": file_name,
                "variable": variable_name,
                "passed": passed,
                **detail,
            }
            results.append(result)
            if not passed:
                failures.append(f"{label}: {json.dumps(detail, sort_keys=True)}")

    return failures, results


def compare_output_directories(
    reference_dir: Path,
    candidate_dir: Path,
    case_id: int,
    required_outputs: Iterable[str],
    tolerance_config: Optional[Mapping[str, Any]] = None,
) -> ComparisonResult:
    """Compare every variable in the requested CPU and GPU outputs.

    NetCDF by default; a prefix in RESTART_PREFIXES selects the restart-file
    reader instead.
    """

    tolerance_config = tolerance_config or {}
    failures: list[str] = []
    variable_results: list[Dict[str, Any]] = []
    files_compared = 0

    for prefix in required_outputs:
        reference_files = _output_files(reference_dir, prefix, case_id)
        candidate_files = _output_files(candidate_dir, prefix, case_id)
        if not reference_files:
            failures.append(
                f"No CPU {prefix} output for case {case_id:03d} in {reference_dir}"
            )
        if not candidate_files:
            failures.append(
                f"No GPU {prefix} output for case {case_id:03d} in {candidate_dir}"
            )

        if set(reference_files) != set(candidate_files):
            failures.append(
                f"{prefix}: output-file set differs: "
                f"CPU={sorted(reference_files)}, GPU={sorted(candidate_files)}"
            )

        compare = _compare_restart_file if prefix in RESTART_PREFIXES else _compare_file
        for file_name in sorted(set(reference_files) & set(candidate_files)):
            file_failures, file_results = compare(
                reference_files[file_name],
                candidate_files[file_name],
                tolerance_config,
            )
            failures.extend(file_failures)
            variable_results.extend(file_results)
            files_compared += 1

    if files_compared == 0:
        failures.append("No output files were compared")
    if not variable_results:
        failures.append("No NetCDF variables were compared")

    return ComparisonResult(
        passed=not failures,
        files_compared=files_compared,
        variables_compared=len(variable_results),
        failures=failures,
        variable_results=variable_results,
    )
