#!/usr/bin/env python3
"""Run deterministic uDALES CPU/GPU parity cases from the GPU test matrix."""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Mapping, Optional, Sequence

if TYPE_CHECKING:
    from compare_outputs import ComparisonResult


TEST_DIR = Path(__file__).resolve().parent
REPO_ROOT = TEST_DIR.parents[2]
DEFAULT_MATRIX = TEST_DIR / "case_matrix.json"
CUDA_SELFTEST_PASS = "CUDA device self-tests passed."
CUDA_SELFTEST_FAIL = "CUDA device self-tests failed:"
CUDA_SELFTEST_ENV = "UDALES_RUN_CUDA_SELFTEST"
CUDA_SELFTEST_COVERAGE = "debug-device-selftests"


class ConfigurationError(RuntimeError):
    """The committed GPU test matrix is invalid."""


def _resolve_repo_path(raw_path: str) -> Path:
    path = Path(raw_path)
    return path if path.is_absolute() else REPO_ROOT / path


def _read_matrix(path: Path) -> Dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        matrix = json.load(handle)
    if not isinstance(matrix, dict):
        raise ConfigurationError(f"Top level of {path} must be an object")
    return matrix


def _read_namelist_integer(path: Path, key: str) -> Optional[int]:
    match = re.search(
        rf"(?mi)^\s*{re.escape(key)}\s*=\s*([0-9]+)\s*$",
        path.read_text(encoding="utf-8"),
    )
    return int(match.group(1)) if match else None


def _read_namelist_scalar(path: Path, key: str) -> Any:
    """Read one simple scalar assignment from a committed Fortran namelist."""

    match = re.search(
        rf"(?mi)^\s*{re.escape(key)}\s*=\s*([^!\n/]+)",
        path.read_text(encoding="utf-8"),
    )
    if match is None:
        raise KeyError(key)
    raw = match.group(1).strip().rstrip(",").strip()
    lowered = raw.lower()
    if lowered in {".true.", "true"}:
        return True
    if lowered in {".false.", "false"}:
        return False
    if (raw.startswith("'") and raw.endswith("'")) or (
        raw.startswith('"') and raw.endswith('"')
    ):
        return raw[1:-1]
    try:
        return int(raw)
    except ValueError:
        try:
            return float(raw.replace("d", "e").replace("D", "e"))
        except ValueError:
            return raw


def _discover_cuda_routines() -> set[str]:
    """Return production CUDA global/device routines declared under src/."""

    declaration = re.compile(
        r"(?im)^\s*attributes\s*\(\s*(?:global|device)\s*\)\s*"
        r"(?:(?:real|integer|logical|double\s+precision)(?:\s*\([^)]*\))?\s+)?"
        r"(?:subroutine|(?:real\s+)?function)\s+(\w+)"
    )
    discovered: set[str] = set()
    for source in sorted((REPO_ROOT / "src").glob("*.f90")):
        if source.name == "tests_cuda.f90":
            continue
        for match in declaration.finditer(source.read_text(encoding="utf-8", errors="replace")):
            discovered.add(f"src/{source.name}:{match.group(1)}")
    return discovered


def _discover_openacc_routines() -> set[str]:
    """Return production subroutines containing executable OpenACC regions."""

    discovered: set[str] = set()
    start = re.compile(
        r"^\s*(?!end\b)(?:[\w(),=*:+-]+\s+)*subroutine\s+(\w+)",
        re.IGNORECASE,
    )
    end = re.compile(r"^\s*end\s*subroutine", re.IGNORECASE)
    directive = re.compile(r"!\$acc\s+(?:parallel|kernels|routine)", re.IGNORECASE)
    for source in sorted((REPO_ROOT / "src").glob("*.f90")):
        stack: list[str] = []
        for line in source.read_text(encoding="utf-8", errors="replace").splitlines():
            match = start.match(line)
            if match:
                stack.append(match.group(1))
            if directive.search(line) and stack:
                discovered.add(f"src/{source.name}:{stack[-1]}")
            if end.match(line) and stack:
                stack.pop()
    return discovered


def validate_matrix(matrix: Mapping[str, Any]) -> None:
    errors: list[str] = []
    cases = matrix.get("cases")
    selections = matrix.get("selections")
    if not isinstance(cases, list) or not cases:
        errors.append("'cases' must be a non-empty list")
        cases = []
    if not isinstance(selections, dict) or not selections:
        errors.append("'selections' must be a non-empty object")
        selections = {}

    names: set[str] = set()
    required_keys = {
        "name",
        "case_id",
        "base_case",
        "namelist",
        "nprocs",
        "required_outputs",
        "coverage",
    }
    for index, case in enumerate(cases):
        label = f"case[{index}]"
        if not isinstance(case, dict):
            errors.append(f"{label} must be an object")
            continue
        missing = sorted(required_keys - set(case))
        if missing:
            errors.append(f"{label} is missing keys: {missing}")
            continue

        name = str(case["name"])
        label = f"case '{name}'"
        if name in names:
            errors.append(f"duplicate case name: {name}")
        names.add(name)

        base_case = _resolve_repo_path(str(case["base_case"]))
        namelist = _resolve_repo_path(str(case["namelist"]))
        nprocx: Optional[int] = None
        nprocy: Optional[int] = None
        if not base_case.is_dir():
            errors.append(f"{label}: base_case does not exist: {base_case}")
        if not namelist.is_file():
            errors.append(f"{label}: namelist does not exist: {namelist}")
        else:
            nprocx = _read_namelist_integer(namelist, "nprocx")
            nprocy = _read_namelist_integer(namelist, "nprocy")
            if nprocx is None or nprocy is None:
                errors.append(f"{label}: namelist must define nprocx and nprocy")
            elif nprocx * nprocy != int(case["nprocs"]):
                errors.append(
                    f"{label}: nprocs={case['nprocs']} but nprocx*nprocy={nprocx*nprocy}"
                )
            namelist_case_id = _read_namelist_integer(namelist, "iexpnr")
            if namelist_case_id != int(case["case_id"]):
                errors.append(
                    f"{label}: case_id={case['case_id']} but namelist iexpnr={namelist_case_id}"
                )

        for raw_override in case.get("overrides", []):
            override = _resolve_repo_path(str(raw_override))
            if not override.is_file():
                errors.append(f"{label}: override does not exist: {override}")

        outputs = case["required_outputs"]
        outputs_are_valid = isinstance(outputs, list) and bool(outputs)
        if not outputs_are_valid:
            errors.append(f"{label}: required_outputs must be a non-empty list")
        coverage = case["coverage"]
        if not isinstance(coverage, list) or not coverage:
            errors.append(f"{label}: coverage must be a non-empty list")
        if int(case["nprocs"]) < 1:
            errors.append(f"{label}: nprocs must be positive")
        if int(case.get("minimum_gpus", 1)) < 1:
            errors.append(f"{label}: minimum_gpus must be positive")
        if float(case.get("timeout_seconds", 300)) <= 0:
            errors.append(f"{label}: timeout_seconds must be positive")

        assertions = case.get("namelist_assertions", {})
        if not isinstance(assertions, dict):
            errors.append(f"{label}: namelist_assertions must be an object")
        elif namelist.is_file():
            for key, expected in assertions.items():
                try:
                    actual = _read_namelist_scalar(namelist, str(key))
                except KeyError:
                    errors.append(f"{label}: asserted namelist key is absent: {key}")
                    continue
                if actual != expected:
                    errors.append(
                        f"{label}: namelist assertion {key}={expected!r}, found {actual!r}"
                    )

        output_assertions = case.get("output_assertions", [])
        if not isinstance(output_assertions, list):
            errors.append(f"{label}: output_assertions must be a list")
        else:
            for assertion_index, assertion in enumerate(output_assertions):
                assertion_label = f"{label}: output_assertions[{assertion_index}]"
                if not isinstance(assertion, dict):
                    errors.append(f"{assertion_label} must be an object")
                    continue
                if assertion.get("type") != "global_peak":
                    errors.append(f"{assertion_label}: unsupported type {assertion.get('type')!r}")
                output = assertion.get("output")
                if not isinstance(output, str):
                    errors.append(f"{assertion_label}: output must be a string")
                elif outputs_are_valid and output not in outputs:
                    errors.append(
                        f"{assertion_label}: output {output!r} is not in required_outputs"
                    )
                if not isinstance(assertion.get("variable"), str):
                    errors.append(f"{assertion_label}: variable must be a string")
                expected_coordinates = assertion.get("expected_coordinates")
                if not isinstance(expected_coordinates, dict) or not expected_coordinates:
                    errors.append(
                        f"{assertion_label}: expected_coordinates must be a non-empty object"
                    )
                else:
                    for coordinate, expected in expected_coordinates.items():
                        try:
                            finite = math.isfinite(float(expected))
                        except (TypeError, ValueError):
                            finite = False
                        if not isinstance(coordinate, str) or not finite:
                            errors.append(
                                f"{assertion_label}: expected coordinate values must be finite numbers"
                            )
                try:
                    coordinate_atol = float(assertion.get("coordinate_atol", 0.0))
                except (TypeError, ValueError):
                    coordinate_atol = float("nan")
                if not math.isfinite(coordinate_atol) or coordinate_atol < 0.0:
                    errors.append(
                        f"{assertion_label}: coordinate_atol must be finite and non-negative"
                    )
                try:
                    minimum_value = float(assertion.get("minimum_value", 0.0))
                except (TypeError, ValueError):
                    minimum_value = float("nan")
                if not math.isfinite(minimum_value):
                    errors.append(f"{assertion_label}: minimum_value must be finite")
                record = assertion.get("record", -1)
                if not isinstance(record, int) or isinstance(record, bool):
                    errors.append(f"{assertion_label}: record must be an integer")
                expected_rank = assertion.get("expected_rank")
                if expected_rank is not None and (
                    not isinstance(expected_rank, list)
                    or len(expected_rank) != 2
                    or any(not isinstance(rank, int) or rank < 0 for rank in expected_rank)
                ):
                    errors.append(
                        f"{assertion_label}: expected_rank must be [x_rank, y_rank]"
                    )
                elif expected_rank is not None and nprocx is not None and nprocy is not None:
                    if expected_rank[0] >= nprocx or expected_rank[1] >= nprocy:
                        errors.append(
                            f"{assertion_label}: expected_rank {expected_rank} lies outside "
                            f"the {nprocx} x {nprocy} process grid"
                        )

    for selection, selected_names in selections.items():
        if not isinstance(selected_names, list) or not selected_names:
            errors.append(f"selection '{selection}' must be a non-empty list")
            continue
        unknown = sorted(set(selected_names) - names)
        if unknown:
            errors.append(f"selection '{selection}' references unknown cases: {unknown}")

    default_tolerance = matrix.get("tolerances", {}).get("default", {})
    for key in ("atol", "rtol"):
        value = float(default_tolerance.get(key, 0.0))
        if value < 0 or not value < float("inf"):
            errors.append(f"default tolerance {key} must be finite and non-negative")

    contract = matrix.get("ported_routines")
    if not isinstance(contract, dict):
        errors.append("'ported_routines' must be an object")
    else:
        known_targets = names | {CUDA_SELFTEST_COVERAGE}
        for kind, discovered in (
            ("cuda", _discover_cuda_routines()),
            ("openacc", _discover_openacc_routines()),
        ):
            declared = contract.get(kind)
            if not isinstance(declared, dict):
                errors.append(f"ported_routines.{kind} must be an object")
                continue
            declared_names = set(declared)
            missing = sorted(discovered - declared_names)
            stale = sorted(declared_names - discovered)
            if missing:
                errors.append(f"ported_routines.{kind} misses source routines: {missing}")
            if stale:
                errors.append(f"ported_routines.{kind} has stale routines: {stale}")
            for routine, targets in declared.items():
                if not isinstance(targets, list) or not targets:
                    errors.append(
                        f"ported_routines.{kind}.{routine} must name at least one test"
                    )
                    continue
                unknown_targets = sorted(set(targets) - known_targets)
                if unknown_targets:
                    errors.append(
                        f"ported_routines.{kind}.{routine} references unknown tests: "
                        f"{unknown_targets}"
                    )

    unsupported = matrix.get("unsupported_gpu_options")
    if not isinstance(unsupported, dict) or not unsupported:
        errors.append("'unsupported_gpu_options' must be a non-empty object")
    else:
        for option, values in unsupported.items():
            if not isinstance(values, list) or not values:
                errors.append(
                    f"unsupported_gpu_options.{option} must be a non-empty list"
                )
                continue
            for case in cases:
                if not isinstance(case, dict) or "namelist" not in case:
                    continue
                namelist = _resolve_repo_path(str(case["namelist"]))
                if not namelist.is_file():
                    continue
                try:
                    actual = _read_namelist_scalar(namelist, str(option))
                except KeyError:
                    continue
                if actual in values:
                    errors.append(
                        f"case '{case.get('name', '?')}': unsupported GPU option "
                        f"{option}={actual!r}"
                    )

    if errors:
        raise ConfigurationError("Invalid GPU test matrix:\n- " + "\n- ".join(errors))


def _case_lookup(matrix: Mapping[str, Any]) -> Dict[str, Dict[str, Any]]:
    return {str(case["name"]): dict(case) for case in matrix["cases"]}


def _selected_cases(matrix: Mapping[str, Any], selection: str) -> list[Dict[str, Any]]:
    if selection not in matrix["selections"]:
        choices = ", ".join(sorted(matrix["selections"]))
        raise ConfigurationError(f"Unknown selection '{selection}'. Available: {choices}")
    lookup = _case_lookup(matrix)
    return [lookup[name] for name in matrix["selections"][selection]]


def _mpi_command(executable: Path, namelist: Path, nprocs: int, gpu: bool) -> list[str]:
    launcher_variable = "UDALES_GPU_MPIEXEC" if gpu else "UDALES_CPU_MPIEXEC"
    extra_variable = "UDALES_GPU_MPI_ARGS" if gpu else "UDALES_CPU_MPI_ARGS"
    mpiexec = os.environ.get(launcher_variable, os.environ.get("MPIEXEC", "mpiexec"))
    extra = shlex.split(
        os.environ.get(extra_variable, os.environ.get("MPI_LAUNCH_EXTRA_ARGS", ""))
    )
    try:
        version = subprocess.run(
            [mpiexec, "--version"],
            check=False,
            capture_output=True,
            text=True,
            timeout=10,
        ).stdout
    except (OSError, subprocess.TimeoutExpired):
        version = ""
    if re.search(r"Open MPI|OpenRTE", version, flags=re.IGNORECASE) and "--oversubscribe" not in extra:
        extra.insert(0, "--oversubscribe")

    command = [mpiexec, *extra, "-n", str(nprocs)]
    use_bind = gpu and nprocs > 1 and os.environ.get("UDALES_GPU_BIND", "1") != "0"
    if use_bind:
        command.extend(["bash", str(REPO_ROOT / "tools" / "bind.sh")])
    command.extend([str(executable), namelist.name])
    return command


def _stage_case(case: Mapping[str, Any], destination: Path) -> Path:
    base_case = _resolve_repo_path(str(case["base_case"]))
    shutil.copytree(base_case, destination)
    case_id = int(case["case_id"])
    target_namelist = destination / f"namoptions.{case_id:03d}"
    shutil.copy2(_resolve_repo_path(str(case["namelist"])), target_namelist)
    for raw_override in case.get("overrides", []):
        override = _resolve_repo_path(str(raw_override))
        shutil.copy2(override, destination / override.name)
    return target_namelist


def _run_solver(
    executable: Path,
    run_dir: Path,
    namelist: Path,
    nprocs: int,
    timeout_seconds: float,
    gpu: bool,
    run_cuda_selftest: bool = False,
) -> Dict[str, Any]:
    log_path = run_dir / "run.log"
    command = _mpi_command(executable, namelist, nprocs, gpu)
    environment = os.environ.copy()
    environment.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    environment.setdefault("TMPDIR", "/tmp")
    # The runner, rather than every Debug executable invocation, owns the
    # policy for executing test-only CUDA assertions.
    environment.pop(CUDA_SELFTEST_ENV, None)
    if gpu and run_cuda_selftest:
        environment[CUDA_SELFTEST_ENV] = "1"
    started = time.monotonic()
    with log_path.open("w", encoding="utf-8") as log:
        log.write("Command: " + shlex.join(command) + "\n")
        log.flush()
        try:
            completed = subprocess.run(
                command,
                cwd=run_dir,
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
                timeout=timeout_seconds,
            )
            returncode = completed.returncode
            timed_out = False
        except subprocess.TimeoutExpired:
            returncode = 124
            timed_out = True
            log.write(f"\nGPU test timed out after {timeout_seconds} seconds\n")
    return {
        "command": command,
        "returncode": returncode,
        "timed_out": timed_out,
        "elapsed_seconds": time.monotonic() - started,
        "log": str(log_path),
    }


def _log_tail(path: Path, lines: int = 80) -> str:
    contents = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return "\n".join(contents[-lines:])


def _copy_artifacts(case_root: Path, artifact_root: Path, case_name: str) -> None:
    destination = artifact_root / case_name
    destination.mkdir(parents=True, exist_ok=True)
    for mode in ("cpu", "gpu"):
        source = case_root / mode
        if not source.is_dir():
            continue
        target = destination / mode
        target.mkdir(exist_ok=True)
        for path in source.iterdir():
            if path.name == "run.log" or path.suffix == ".nc":
                shutil.copy2(path, target / path.name)


def _write_junit(results: Sequence[Mapping[str, Any]], path: Path) -> None:
    testsuite = ET.Element(
        "testsuite",
        name="uDALES GPU parity",
        tests=str(len(results)),
        failures=str(sum(not bool(result["passed"]) for result in results)),
        time=f"{sum(float(result.get('elapsed_seconds', 0.0)) for result in results):.6f}",
    )
    for result in results:
        testcase = ET.SubElement(
            testsuite,
            "testcase",
            classname="gpu.cpu_parity",
            name=str(result["name"]),
            time=f"{float(result.get('elapsed_seconds', 0.0)):.6f}",
        )
        if not result["passed"]:
            failure = ET.SubElement(testcase, "failure", message="CPU/GPU parity failure")
            failure.text = "\n".join(str(item) for item in result.get("failures", []))
    path.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(testsuite).write(path, encoding="utf-8", xml_declaration=True)


def _check_executables(cpu_executable: Path, gpu_executable: Path) -> None:
    errors = []
    for label, path in (("CPU", cpu_executable), ("GPU", gpu_executable)):
        if not path.is_file():
            errors.append(f"{label} executable not found: {path}")
        elif not os.access(path, os.X_OK):
            errors.append(f"{label} executable is not executable: {path}")
    if errors:
        raise ConfigurationError("\n".join(errors))
    if cpu_executable.resolve() == gpu_executable.resolve():
        raise ConfigurationError("CPU and GPU executable paths must be different")

    linked_libraries: Dict[str, str] = {}
    for label, path in (("CPU", cpu_executable), ("GPU", gpu_executable)):
        completed = subprocess.run(
            ["ldd", str(path)],
            check=False,
            capture_output=True,
            text=True,
        )
        linked_libraries[label] = completed.stdout + completed.stderr
    cuda_pattern = re.compile(r"libcuda|libcudart|libcufft|libcudafor", re.IGNORECASE)
    if cuda_pattern.search(linked_libraries["CPU"]):
        raise ConfigurationError(f"CPU executable links CUDA libraries: {cpu_executable}")
    if not cuda_pattern.search(linked_libraries["GPU"]):
        raise ConfigurationError(f"GPU executable does not link CUDA libraries: {gpu_executable}")


def _check_cpu_executable(cpu_executable: Path) -> None:
    if not cpu_executable.is_file():
        raise ConfigurationError(f"CPU executable not found: {cpu_executable}")
    if not os.access(cpu_executable, os.X_OK):
        raise ConfigurationError(f"CPU executable is not executable: {cpu_executable}")


def _visible_gpu_count() -> Optional[int]:
    visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    if visible is not None:
        devices = [item.strip() for item in visible.split(",") if item.strip() and item.strip() != "-1"]
        return len(devices)
    try:
        completed = subprocess.run(
            ["nvidia-smi", "-L"],
            check=False,
            capture_output=True,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    if completed.returncode != 0:
        return None
    return sum(line.startswith("GPU ") for line in completed.stdout.splitlines())


def _require_cuda_selftest(log_path: Path, expected_ranks: int) -> Optional[str]:
    log = log_path.read_text(encoding="utf-8", errors="replace")
    if CUDA_SELFTEST_FAIL in log:
        return "CUDA device self-test suite reported failure"
    passed_ranks = [
        int(rank)
        for rank in re.findall(
            re.escape(CUDA_SELFTEST_PASS) + r"\s*rank=(\d+)",
            log,
        )
    ]
    if not passed_ranks:
        return (
            "CUDA device self-test suite did not run; "
            "use a Debug GPU executable or omit --require-debug-selftest"
        )
    expected = list(range(expected_ranks))
    if sorted(passed_ranks) != expected:
        return (
            "CUDA device self-test suite did not pass exactly once "
            f"on every MPI rank: expected={expected}, observed={sorted(passed_ranks)}"
        )
    return None


def _validate_output_assertions(
    case: Mapping[str, Any],
    run_dir: Path,
    mode: str,
) -> list[str]:
    """Validate semantic properties that CPU/GPU parity alone cannot prove."""

    import netCDF4
    import numpy as np

    failures: list[str] = []
    case_id = int(case["case_id"])
    for assertion in case.get("output_assertions", []):
        output = str(assertion["output"])
        variable_name = str(assertion["variable"])
        paths = sorted(run_dir.glob(f"{output}*.{case_id:03d}.nc"))
        label = f"{mode} {output}:{variable_name} global-peak assertion"
        if not paths:
            failures.append(f"{label}: no matching output files")
            continue

        best: Optional[Dict[str, Any]] = None
        for path in paths:
            with netCDF4.Dataset(path) as dataset:
                if variable_name not in dataset.variables:
                    failures.append(f"{label}: {variable_name} is absent from {path.name}")
                    continue
                variable = dataset.variables[variable_name]
                dimensions = list(variable.dimensions)
                values = np.ma.asarray(variable[:])
                if "time" in dimensions:
                    time_axis = dimensions.index("time")
                    record = int(assertion.get("record", -1))
                    record_count = values.shape[time_axis]
                    resolved_record = record if record >= 0 else record_count + record
                    if resolved_record < 0 or resolved_record >= record_count:
                        failures.append(
                            f"{label}: record {record} is outside {path.name}'s "
                            f"{record_count} time records"
                        )
                        continue
                    values = np.take(values, resolved_record, axis=time_axis)
                    dimensions.pop(time_axis)

                data = np.asarray(values.data)
                valid = ~np.ma.getmaskarray(values) & np.isfinite(data)
                if not np.any(valid):
                    failures.append(f"{label}: {path.name} has no finite unmasked values")
                    continue
                peak_indices = np.unravel_index(
                    int(np.argmax(np.where(valid, data, -np.inf))),
                    data.shape,
                )
                peak_value = float(data[peak_indices])
                coordinates: Dict[str, float] = {}
                coordinate_error = False
                for coordinate in assertion["expected_coordinates"]:
                    if coordinate not in dimensions or coordinate not in dataset.variables:
                        failures.append(
                            f"{label}: coordinate {coordinate!r} is unavailable in {path.name}"
                        )
                        coordinate_error = True
                        continue
                    coordinate_axis = dimensions.index(coordinate)
                    coordinate_values = np.asarray(dataset.variables[coordinate][:])
                    coordinates[coordinate] = float(coordinate_values[peak_indices[coordinate_axis]])
                if coordinate_error:
                    continue

                rank_match = re.fullmatch(
                    rf"{re.escape(output)}\.(\d+)\.(\d+)\.{case_id:03d}\.nc",
                    path.name,
                )
                rank = (
                    [int(rank_match.group(1)), int(rank_match.group(2))]
                    if rank_match is not None
                    else None
                )
                candidate = {
                    "value": peak_value,
                    "coordinates": coordinates,
                    "rank": rank,
                    "file": path.name,
                }
                if best is None or peak_value > float(best["value"]):
                    best = candidate

        if best is None:
            continue
        minimum_value = float(assertion.get("minimum_value", 0.0))
        if float(best["value"]) <= minimum_value:
            failures.append(
                f"{label}: peak {best['value']} does not exceed {minimum_value}; "
                "the source may not have been applied"
            )
        coordinate_atol = float(assertion.get("coordinate_atol", 0.0))
        for coordinate, expected in assertion["expected_coordinates"].items():
            actual = float(best["coordinates"][coordinate])
            if abs(actual - float(expected)) > coordinate_atol:
                failures.append(
                    f"{label}: peak {coordinate}={actual}, expected {expected} "
                    f"within {coordinate_atol}; file={best['file']}"
                )
        expected_rank = assertion.get("expected_rank")
        if expected_rank is not None and best["rank"] != expected_rank:
            failures.append(
                f"{label}: peak rank={best['rank']}, expected {expected_rank}; "
                f"file={best['file']}"
            )
    return failures


def run_cpu_fixture(
    case: Mapping[str, Any],
    cpu_executable: Path,
    run_root: Path,
    artifact_root: Optional[Path],
    tolerance_config: Mapping[str, Any],
) -> Dict[str, Any]:
    """Validate that a committed fixture runs and produces finite expected output."""

    from compare_outputs import compare_output_directories

    name = str(case["name"])
    print(f"\n==> CPU fixture validation: {name}")
    case_root = run_root / name
    cpu_dir = case_root / "cpu"
    namelist = _stage_case(case, cpu_dir)
    cpu_run = _run_solver(
        cpu_executable,
        cpu_dir,
        namelist,
        int(case["nprocs"]),
        float(case.get("timeout_seconds", 300)),
        gpu=False,
    )
    failures: list[str] = []
    comparison: Optional["ComparisonResult"] = None
    if cpu_run["returncode"] != 0:
        failures.append(
            f"CPU run exited with {cpu_run['returncode']}\n{_log_tail(cpu_dir / 'run.log')}"
        )
    else:
        # Self-comparison applies the complete file/variable/finite-value
        # contract without needing a GPU result.
        comparison = compare_output_directories(
            cpu_dir,
            cpu_dir,
            int(case["case_id"]),
            case["required_outputs"],
            tolerance_config,
        )
        failures.extend(comparison.failures)
        failures.extend(_validate_output_assertions(case, cpu_dir, "CPU"))
    result = {
        "name": name,
        "passed": not failures,
        "coverage": case["coverage"],
        "failures": failures,
        "cpu_run": cpu_run,
        "gpu_run": None,
        "elapsed_seconds": float(cpu_run["elapsed_seconds"]),
        "comparison": comparison.as_dict() if comparison else None,
    }
    if artifact_root is not None:
        _copy_artifacts(case_root, artifact_root, name)
    print(f"result: {'PASS' if result['passed'] else 'FAIL'}")
    for failure in failures:
        print(f"  {failure}", file=sys.stderr)
    return result


def run_case(
    case: Mapping[str, Any],
    cpu_executable: Path,
    gpu_executable: Path,
    run_root: Path,
    artifact_root: Optional[Path],
    tolerance_config: Mapping[str, Any],
    require_debug_selftest: bool,
) -> Dict[str, Any]:
    name = str(case["name"])
    print(f"\n==> GPU parity case: {name}")
    print("coverage: " + ", ".join(case["coverage"]))
    case_root = run_root / name
    cpu_dir = case_root / "cpu"
    gpu_dir = case_root / "gpu"
    cpu_namelist = _stage_case(case, cpu_dir)
    gpu_namelist = _stage_case(case, gpu_dir)
    timeout_seconds = float(case.get("timeout_seconds", 300))
    nprocs = int(case["nprocs"])

    cpu_run = _run_solver(
        cpu_executable,
        cpu_dir,
        cpu_namelist,
        nprocs,
        timeout_seconds,
        gpu=False,
    )
    gpu_run = _run_solver(
        gpu_executable,
        gpu_dir,
        gpu_namelist,
        nprocs,
        timeout_seconds,
        gpu=True,
        run_cuda_selftest=require_debug_selftest,
    )
    failures: list[str] = []
    if cpu_run["returncode"] != 0:
        failures.append(
            f"CPU run exited with {cpu_run['returncode']}\n{_log_tail(cpu_dir / 'run.log')}"
        )
    if gpu_run["returncode"] != 0:
        failures.append(
            f"GPU run exited with {gpu_run['returncode']}\n{_log_tail(gpu_dir / 'run.log')}"
        )
    if require_debug_selftest and gpu_run["returncode"] == 0:
        selftest_failure = _require_cuda_selftest(gpu_dir / "run.log", nprocs)
        if selftest_failure:
            failures.append(selftest_failure)

    comparison: Optional["ComparisonResult"] = None
    if cpu_run["returncode"] == 0 and gpu_run["returncode"] == 0:
        from compare_outputs import compare_output_directories

        comparison = compare_output_directories(
            cpu_dir,
            gpu_dir,
            int(case["case_id"]),
            case["required_outputs"],
            tolerance_config,
        )
        failures.extend(comparison.failures)
        failures.extend(_validate_output_assertions(case, cpu_dir, "CPU"))
        failures.extend(_validate_output_assertions(case, gpu_dir, "GPU"))

    elapsed = float(cpu_run["elapsed_seconds"]) + float(gpu_run["elapsed_seconds"])
    result = {
        "name": name,
        "passed": not failures,
        "coverage": case["coverage"],
        "failures": failures,
        "cpu_run": cpu_run,
        "gpu_run": gpu_run,
        "elapsed_seconds": elapsed,
        "comparison": comparison.as_dict() if comparison else None,
    }
    if artifact_root is not None:
        _copy_artifacts(case_root, artifact_root, name)
    status = "PASS" if result["passed"] else "FAIL"
    print(f"result: {status} ({elapsed:.2f} s)")
    for failure in failures:
        print(f"  {failure}", file=sys.stderr)
    return result


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("selection", nargs="?", default="smoke")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument(
        "--cpu-executable",
        type=Path,
        default=Path(os.environ.get("UDALES_CPU_BUILD", REPO_ROOT / "build/cpu/debug/u-dales")),
    )
    parser.add_argument(
        "--gpu-executable",
        type=Path,
        default=Path(os.environ.get("UDALES_GPU_BUILD", REPO_ROOT / "build/gpu/debug/u-dales")),
    )
    parser.add_argument("--work-root", type=Path, default=None)
    parser.add_argument(
        "--artifacts-dir",
        type=Path,
        default=Path(os.environ["UDALES_GPU_TEST_ARTIFACTS"])
        if "UDALES_GPU_TEST_ARTIFACTS" in os.environ
        else None,
    )
    parser.add_argument("--junit", type=Path, default=None)
    parser.add_argument("--summary", type=Path, default=None)
    parser.add_argument("--require-debug-selftest", action="store_true")
    parser.add_argument("--validate-config", action="store_true")
    parser.add_argument(
        "--cpu-only",
        action="store_true",
        help="Run and validate the CPU fixtures without requiring a GPU executable.",
    )
    parser.add_argument("--keep-work-dir", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    matrix = _read_matrix(args.matrix.resolve())
    validate_matrix(matrix)
    selected = _selected_cases(matrix, args.selection)
    if args.validate_config:
        print(f"PASSED: GPU test matrix configuration is valid ({len(matrix['cases'])} cases).")
        print(f"Selection '{args.selection}': {', '.join(case['name'] for case in selected)}")
        print("Scope: configuration only; no solver executable was built or run and no GPU was used.")
        print("Next: build and run the Debug GPU smoke tests with:")
        print("  python tests/run_tests.py gpu-smoke")
        return 0

    print(f"GPU test matrix is valid: {len(matrix['cases'])} cases")
    print(f"selection '{args.selection}': {', '.join(case['name'] for case in selected)}")

    cpu_executable = args.cpu_executable.resolve()
    gpu_executable = args.gpu_executable.resolve()
    if args.cpu_only:
        _check_cpu_executable(cpu_executable)
    else:
        _check_executables(cpu_executable, gpu_executable)
        bind_enabled = os.environ.get("UDALES_GPU_BIND", "1") != "0"
        required_gpus = (
            max(int(case.get("minimum_gpus", 1)) for case in selected)
            if bind_enabled
            else 1
        )
        visible_gpus = _visible_gpu_count()
        if visible_gpus is not None and visible_gpus < required_gpus:
            raise ConfigurationError(
                f"selection '{args.selection}' requires {required_gpus} GPU(s), "
                f"but only {visible_gpus} are visible"
            )
        if visible_gpus is None:
            print("warning: unable to determine visible GPU count; runtime launch will verify availability")

    work_parent = args.work_root.resolve() if args.work_root else None
    if work_parent:
        work_parent.mkdir(parents=True, exist_ok=True)
    temporary = tempfile.TemporaryDirectory(prefix="udales-gpu-parity-", dir=work_parent)
    run_root = Path(temporary.name)
    print(f"run root: {run_root}")

    artifact_root = args.artifacts_dir.resolve() if args.artifacts_dir else None
    if artifact_root:
        artifact_root.mkdir(parents=True, exist_ok=True)

    if args.cpu_only:
        results = [
            run_cpu_fixture(
                case,
                cpu_executable,
                run_root,
                artifact_root,
                matrix.get("tolerances", {}),
            )
            for case in selected
        ]
    else:
        results = [
            run_case(
                case,
                cpu_executable,
                gpu_executable,
                run_root,
                artifact_root,
                matrix.get("tolerances", {}),
                args.require_debug_selftest,
            )
            for case in selected
        ]
    summary = {
        "selection": args.selection,
        "passed": all(result["passed"] for result in results),
        "cpu_executable": str(cpu_executable),
        "gpu_executable": None if args.cpu_only else str(gpu_executable),
        "results": results,
    }

    summary_path = args.summary or (artifact_root / "summary.json" if artifact_root else None)
    if summary_path:
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    junit_path = args.junit or (artifact_root / "junit.xml" if artifact_root else None)
    if junit_path:
        _write_junit(results, junit_path)

    print("\nSummary")
    for result in results:
        print(f"- {result['name']}: {'PASS' if result['passed'] else 'FAIL'}")
    print(f"overall: {'PASS' if summary['passed'] else 'FAIL'}")

    if args.keep_work_dir:
        kept_path = Path(tempfile.mkdtemp(prefix="udales-gpu-kept-", dir=work_parent))
        shutil.copytree(run_root, kept_path / "runs")
        print(f"kept run directory: {kept_path / 'runs'}")
    temporary.cleanup()
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ConfigurationError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
