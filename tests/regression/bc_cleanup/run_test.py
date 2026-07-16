#!/usr/bin/env python3

"""Phase 0 regression gate for the BCcleanup refactor (BCcleanup_backlog.md).

Compares dumped fields between two refs, run under two MPI decompositions,
for two cases:

- 090: matched-anchor case, Vreman SGS (`lvreman = .true.`, the default).
  Bitwise across the whole branch (`default_atol = 0.0`) — the Task 4
  SGS re-point (`src/modsubgrid.f90:387,486`) lives in the `loneeqn`
  branch of `closure`, which 090 never executes, so that commit is a
  no-op here.
- 092: same fixture as 090 but with the one-equation TKE closure enabled
  (`lvreman = .false.`, `loneeqn = .true.`), so the Task 4 re-point
  actually runs. `default_atol = 5e-3` for this case (constant-thvs to
  evolving-thvf(k) is a deliberate physics change; see Task 4 and the
  correction note in `BCcleanup_backlog.md`).

Each `CaseSpec` carries its own `default_atol`, used when `--atol` is not
passed; an explicit `--atol` overrides the tolerance for every case.
Adapted from tests/regression/mpi_averaging_regression/run_test.py: the
worktree/build/run/stitch machinery is unchanged; the case matrix,
decomposition set, and field comparison were replaced for a tolerance-driven
gate that every BCcleanup commit is checked against.
"""

import argparse
import collections
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from shlex import quote as shlex_quote
from typing import Dict, List, Optional, Sequence, Tuple

import netCDF4 as nc
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REGRESSION_DIR = SCRIPT_DIR.parent
TESTS_DIR = REGRESSION_DIR.parent
REPO_ROOT = TESTS_DIR.parent

BUILD_SYSTEM = os.environ.get("UDALES_BUILD_SYSTEM", "icl")
GIT = os.environ.get("UDALES_GIT", "/usr/bin/git")
SUBMODULE_DIRS = ("2decomp-fft", "tools/View3D")
RUNTIME_MODULES = os.environ.get(
    "UDALES_RUNTIME_MODULES",
    "intel/2025a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a "
    "FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0",
)
DEFAULT_WORKDIR = Path(os.environ.get("UDALES_SCRATCH_ROOT", "/tmp")).resolve()
CACHED_FFTW_MODULE = REPO_ROOT / "build" / "debug" / "findFFTW-src" / "FindFFTW.cmake"

# Tolerance used only to match up time records between reference and current
# outputs (NOT the field comparison tolerance, which is --atol).
TIME_ATOL = 1.0e-12
TIME_RTOL = 1.0e-10

LOG_PATTERNS = {
    "courant": re.compile(r"Courant numbers \(x,y,z,tot\):\s*([Ee0-9+\-\.]+)"),
    "diffusion": re.compile(r"Diffusion number:\s*([Ee0-9+\-\.]+)"),
    "divergence": re.compile(r"divmax, divtot =\s*([Ee0-9+\-\.]+)\s+([Ee0-9+\-\.]+)"),
}

CaseSpec = collections.namedtuple("CaseSpec", "case nc_pattern fields abs_tol default_atol")

CASES = (
    CaseSpec(
        case="090",
        nc_pattern="tdump.*.090.nc",   # match the per-tile dump files case 090 produces
        fields=None,                    # None = compare every variable with ndim >= 3
        abs_tol=None,                   # filled from --atol or default_atol in main()
        default_atol=0.0,               # Vreman SGS: Task 4 re-point is dead code here (bitwise)
    ),
    CaseSpec(
        case="092",
        nc_pattern="tdump.*.092.nc",   # match the per-tile dump files case 092 produces
        fields=None,                    # None = compare every variable with ndim >= 3
        abs_tol=None,                   # filled from --atol or default_atol in main()
        default_atol=5.0e-3,            # one-equation SGS: Task 4 re-point is live here
    ),
)

CONFIGS = {
    "serial":   dict(nprocx=1, nprocy=1),
    "xy_split": dict(nprocx=2, nprocy=2),
}

# Cases available to --assert-only (single-branch smoke checks, not the
# reference-vs-current field comparison used by CASES/_run_case_matrix).
ASSERT_ONLY_CASES = {
    "090": CaseSpec(case="090", nc_pattern="tdump.*.090.nc", fields=None, abs_tol=None, default_atol=0.0),
    # Buried-slab fixture (backlog 6.6): floor-covering box, 2 grid cells
    # tall, so slabs kb and kb+1 are fully solid (IIcs(k) == 0). Exercises
    # the base-state fallback in diagfld (src/modthermodynamics.f90).
    "091": CaseSpec(case="091", nc_pattern="tdump.*.091.nc", fields=None, abs_tol=None, default_atol=0.0),
    "092": CaseSpec(case="092", nc_pattern="tdump.*.092.nc", fields=None, abs_tol=None, default_atol=5.0e-3),
}

BASE_STATE_LOG_MARKER = "Base state:"
SENTINEL_VALUE = -999.0
SENTINEL_ATOL = 1.0e-6


def _read_int_setting(namelist: Path, key: str) -> int:
    match = re.search(
        rf"(?m)^\s*{re.escape(key)}\s*=\s*(\d+)\s*$",
        namelist.read_text(encoding="utf-8"),
    )
    if not match:
        raise RuntimeError(f"Missing integer setting '{key}' in {namelist}")
    return int(match.group(1))


def _configure_namelist(run_dir: Path, spec: CaseSpec, nprocx: int, nprocy: int) -> None:
    target = run_dir / f"namoptions.{spec.case}"
    if not target.is_file():
        raise RuntimeError(f"Namelist not found after copying case: {target}")
    text = target.read_text(encoding="utf-8")
    text = re.sub(r"(?m)^(\s*nprocx\s*=\s*)\d+(\s*)$", rf"\g<1>{nprocx}\2", text)
    text = re.sub(r"(?m)^(\s*nprocy\s*=\s*)\d+(\s*)$", rf"\g<1>{nprocy}\2", text)
    target.write_text(text, encoding="utf-8")


def _tile_index_pattern(nc_pattern: str) -> "re.Pattern[str]":
    # nc_pattern is a shell glob such as "tdump.*.090.nc" where the single
    # "*" spans the two zero-padded tile indices, e.g. "000.000".
    escaped = re.escape(nc_pattern)
    escaped = escaped.replace(re.escape("*"), r"(\d{3})\.(\d{3})")
    return re.compile(f"^{escaped}$")


def _parse_tile_indices(path: Path, spec: CaseSpec) -> Tuple[int, int]:
    match = _tile_index_pattern(spec.nc_pattern).match(path.name)
    if not match:
        raise RuntimeError(f"Unexpected output file name: {path.name}")
    return int(match.group(1)), int(match.group(2))


def _dump_files(run_dir: Path, spec: CaseSpec) -> List[Path]:
    return sorted(run_dir.glob(spec.nc_pattern))


def _discover_fields(path: Path) -> List[str]:
    with nc.Dataset(path) as ds:
        return sorted(name for name, var in ds.variables.items() if var.ndim >= 3)


def _load_available_times(run_dir: Path, spec: CaseSpec) -> np.ndarray:
    files = _dump_files(run_dir, spec)
    if not files:
        raise RuntimeError(f"No files matching '{spec.nc_pattern}' found in {run_dir}")
    with nc.Dataset(files[0]) as ds:
        if "time" not in ds.variables:
            raise RuntimeError(f"Missing time variable in {files[0]}")
        times = np.asarray(ds.variables["time"][:], dtype=np.float64)
    if times.size == 0:
        raise RuntimeError(f"No time records found in {files[0]}")
    return times


def _select_common_time(reference_run_dir: Path, current_run_dir: Path, spec: CaseSpec) -> float:
    ref_times = _load_available_times(reference_run_dir, spec)
    cur_times = _load_available_times(current_run_dir, spec)

    common: List[float] = []
    for ref_time in ref_times:
        if np.any(np.isclose(cur_times, ref_time, atol=TIME_ATOL, rtol=TIME_RTOL)):
            common.append(float(ref_time))
    if not common:
        raise RuntimeError(
            f"No common output times between {reference_run_dir} ({ref_times}) and "
            f"{current_run_dir} ({cur_times})"
        )
    return max(common)


def _load_global_fields(run_dir: Path, spec: CaseSpec, target_time: float) -> Dict[str, np.ndarray]:
    namelist = run_dir / f"namoptions.{spec.case}"
    itot = _read_int_setting(namelist, "itot")
    jtot = _read_int_setting(namelist, "jtot")
    nprocx = _read_int_setting(namelist, "nprocx")
    nprocy = _read_int_setting(namelist, "nprocy")

    files = _dump_files(run_dir, spec)
    if not files:
        raise RuntimeError(f"No files matching '{spec.nc_pattern}' found in {run_dir}")

    field_names = list(spec.fields) if spec.fields is not None else _discover_fields(files[0])
    if not field_names:
        raise RuntimeError(f"No variables with ndim >= 3 found in {files[0]}")

    fields: Dict[str, np.ndarray] = {}
    for path in files:
        px, py = _parse_tile_indices(path, spec)
        with nc.Dataset(path) as ds:
            if "time" not in ds.variables:
                raise RuntimeError(f"Missing time variable in {path}")
            times = np.asarray(ds.variables["time"][:], dtype=np.float64)
            if times.size == 0:
                raise RuntimeError(f"No time records found in {path}")
            matches = np.where(np.isclose(times, target_time, atol=TIME_ATOL, rtol=TIME_RTOL))[0]
            if matches.size == 0:
                raise RuntimeError(f"No record at time {target_time} found in {path}; available times: {times}")
            record = int(matches[-1])
            istart = px * (itot // nprocx)
            jstart = py * (jtot // nprocy)
            for field in field_names:
                if field not in ds.variables:
                    raise RuntimeError(f"Field '{field}' not present in {path}")
                arr = np.asarray(ds.variables[field][record], dtype=np.float64)
                nz, ny, nx = arr.shape
                if field not in fields:
                    fields[field] = np.zeros((nz, jtot, itot), dtype=np.float64)
                fields[field][:nz, jstart : jstart + ny, istart : istart + nx] = arr
    return fields


def _mpi_exec_and_args() -> Tuple[str, str]:
    mpiexec = os.environ.get("MPIEXEC")
    if not mpiexec:
        mpiexec = "mpiexec"

    extra_args = os.environ.get("MPI_LAUNCH_EXTRA_ARGS", "").strip()
    try:
        version_proc = subprocess.run(
            [mpiexec, "--version"],
            check=False,
            capture_output=True,
            text=True,
        )
        version = (version_proc.stdout or "") + "\n" + (version_proc.stderr or "")
    except OSError:
        version = ""

    is_openmpi = re.search(r"Open MPI|OpenRTE", version, flags=re.IGNORECASE) is not None
    if is_openmpi and "--oversubscribe" not in extra_args:
        extra_args = f"--oversubscribe {extra_args}".strip()
    elif (not is_openmpi) and "--oversubscribe" in extra_args:
        extra_args = " ".join(arg for arg in extra_args.split() if arg != "--oversubscribe")
    return mpiexec, extra_args


def _run_case(path_to_exe: Path, run_dir: Path, spec: CaseSpec, nprocs: int) -> Path:
    mpiexec, extra_args = _mpi_exec_and_args()
    exe = path_to_exe / "u-dales"
    log_path = run_dir / "solver.log"
    diag_log = run_dir / "fortran_diagnostics.log"

    for stale in run_dir.glob("monitor*.txt"):
        stale.unlink()
    for stale in run_dir.glob(f"*.{spec.case}.nc"):
        stale.unlink()
    if diag_log.exists():
        diag_log.unlink()
    if log_path.exists():
        log_path.unlink()

    module_setup = ""
    if Path("/etc/profile.d/modules.sh").is_file():
        module_setup = "source /etc/profile.d/modules.sh >/dev/null 2>&1 || true; "
    shell_command = (
        f"{module_setup}"
        f"if command -v module >/dev/null 2>&1; then module load {RUNTIME_MODULES}; fi && "
        "export HDF5_USE_FILE_LOCKING=FALSE && "
        "export FOR_DISABLE_DIAGNOSTIC_DISPLAY=TRUE && "
        f"export FOR_DIAGNOSTIC_LOG_FILE={shlex_quote(str(diag_log))} && "
        f"cd {shlex_quote(str(run_dir))} && "
        f"{shlex_quote(mpiexec)} {extra_args} -n {nprocs} "
        f"{shlex_quote(str(exe))} namoptions.{spec.case}"
    )
    with log_path.open("w", encoding="utf-8") as log_handle:
        completed = subprocess.run(
            ["bash", "-lc", shell_command],
            cwd=REPO_ROOT,
            check=False,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if completed.returncode == 0:
        return log_path

    def _tail(label: str, text: str, limit: int = 40) -> str:
        stripped = text.strip()
        if not stripped:
            return f"{label}: <empty>"
        lines = stripped.splitlines()
        snippet = "\n".join(lines[-limit:])
        return f"{label} (last {min(len(lines), limit)} lines):\n{snippet}"

    details = [
        f"u-dales run failed for case {spec.case} in {run_dir}",
        f"command: {shell_command}",
        f"exit code: {completed.returncode}",
        _tail("solver log", log_path.read_text(encoding="utf-8", errors="replace")),
    ]
    if diag_log.exists():
        details.append(_tail("fortran diagnostics", diag_log.read_text(encoding="utf-8", errors="replace")))
    raise RuntimeError("\n\n".join(details))


def _extract_log_series(log_path: Path) -> Dict[str, List[Tuple[float, ...]]]:
    text = log_path.read_text(encoding="utf-8")
    series: Dict[str, List[Tuple[float, ...]]] = {}
    for key, pattern in LOG_PATTERNS.items():
        matches = []
        for match in pattern.finditer(text):
            matches.append(tuple(float(group) for group in match.groups()))
        series[key] = matches
    return series


def _safe_relative(diff_value: float, ref_value: float) -> float:
    if ref_value <= 0.0:
        return 0.0 if diff_value == 0.0 else float("inf")
    return diff_value / ref_value


def _compare_logs(reference_log: Path, current_log: Path, label: str) -> None:
    ref = _extract_log_series(reference_log)
    cur = _extract_log_series(current_log)
    failures: List[str] = []
    for key in sorted(LOG_PATTERNS):
        ref_series = ref[key]
        cur_series = cur[key]
        if len(ref_series) != len(cur_series):
            failures.append(f"{label}: {key} occurrence count mismatch {len(ref_series)} vs {len(cur_series)}")
            continue
        for idx, (ref_vals, cur_vals) in enumerate(zip(ref_series, cur_series)):
            for comp, (ref_val, cur_val) in enumerate(zip(ref_vals, cur_vals)):
                if np.isclose(cur_val, ref_val, atol=TIME_ATOL, rtol=TIME_RTOL):
                    continue
                abs_val = abs(cur_val - ref_val)
                rel_val = _safe_relative(abs_val, abs(ref_val))
                failures.append(
                    f"{label}: {key}[{idx}][{comp}] abs {abs_val:.6e} rel {rel_val:.6e} "
                    f"ref {ref_val:.6e} cur {cur_val:.6e}"
                )
    if failures:
        raise RuntimeError("Solver diagnostic drift exceeds tolerance:\n- " + "\n- ".join(failures))


def _compare_fields(
    spec: CaseSpec,
    reference: Dict[str, np.ndarray],
    current: Dict[str, np.ndarray],
    label: str,
) -> None:
    ref_fields = sorted(reference.keys())
    cur_fields = sorted(current.keys())
    if ref_fields != cur_fields:
        raise RuntimeError(
            f"{label}: field set mismatch between reference and current outputs\n"
            f"  reference: {ref_fields}\n  current:   {cur_fields}"
        )
    if not ref_fields:
        raise RuntimeError(f"{label}: no fields with ndim >= 3 found to compare")

    failures: List[str] = []
    for field in ref_fields:
        if reference[field].shape != current[field].shape:
            failures.append(
                f"{label}: {field} shape mismatch ref {reference[field].shape} cur {current[field].shape}"
            )
            continue
        # NaN comparisons are always False, so `max_diff > spec.abs_tol` below would
        # silently pass if either side contains NaN (np.abs(nan - x).max() is nan,
        # and `nan > anything` is False). Fail loudly instead.
        cur_nan_count = int(np.isnan(current[field]).sum())
        ref_nan_count = int(np.isnan(reference[field]).sum())
        if cur_nan_count or ref_nan_count:
            failures.append(
                f"{label}: {field} contains NaN (current: {cur_nan_count}, reference: {ref_nan_count})"
            )
            continue
        diff = np.abs(current[field] - reference[field])
        flat = int(np.argmax(diff))
        idx = tuple(int(v) for v in np.unravel_index(flat, diff.shape))
        max_diff = float(diff[idx])
        print(f"    {label}: {field} max|diff| = {max_diff:.6e}", flush=True)
        if max_diff > spec.abs_tol:
            ref_val = float(reference[field][idx])
            cur_val = float(current[field][idx])
            failures.append(
                f"{label}: {field} max abs {max_diff:.6e} at {idx} "
                f"ref {ref_val:.6e} cur {cur_val:.6e}"
            )
    if failures:
        raise RuntimeError("Field regression differences exceed tolerance:\n- " + "\n- ".join(failures))


def _create_worktree(worktree_root: Path, label: str, ref: str) -> Path:
    path = worktree_root / label
    try:
        subprocess.run([GIT, "worktree", "add", "--detach", str(path), ref], cwd=REPO_ROOT, check=True)
    except subprocess.CalledProcessError:
        subprocess.run([GIT, "clone", "--shared", "--no-checkout", str(REPO_ROOT), str(path)], cwd=worktree_root, check=True)
        subprocess.run([GIT, "checkout", "--detach", ref], cwd=path, check=True)
    _populate_local_submodules(path)
    return path


def _remove_worktree(path: Path) -> None:
    if not path.exists():
        return
    git_dir = path / ".git"
    if git_dir.is_dir():
        shutil.rmtree(path)
        return
    subprocess.run([GIT, "worktree", "remove", "--force", str(path)], cwd=REPO_ROOT, check=True)


def _populate_local_submodules(worktree: Path) -> None:
    for relative in SUBMODULE_DIRS:
        source = REPO_ROOT / relative
        target = worktree / relative
        if not source.exists():
            raise RuntimeError(f"Required local submodule directory is missing: {source}")
        if target.exists():
            if target.is_dir():
                shutil.rmtree(target)
            else:
                target.unlink()
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(source, target, symlinks=True)


def _prepare_offline_fftw(worktree: Path, build_dir: Path) -> None:
    if not CACHED_FFTW_MODULE.is_file():
        raise RuntimeError(
            f"Missing cached FindFFTW.cmake required for offline regression builds: {CACHED_FFTW_MODULE}"
        )
    cmake_lists = worktree / "CMakeLists.txt"
    text = cmake_lists.read_text(encoding="utf-8")

    def _copy_cached_fftw(dst: Path) -> None:
        dst.parent.mkdir(parents=True, exist_ok=True)
        if CACHED_FFTW_MODULE.resolve() == dst.resolve():
            return
        shutil.copy2(CACHED_FFTW_MODULE, dst)

    if 'set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)' in text:
        findfftw_dir = build_dir / "findFFTW-src"
        _copy_cached_fftw(findfftw_dir / "FindFFTW.cmake")
        return
    if 'set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")' in text:
        cmake_dir = worktree / "cmake"
        _copy_cached_fftw(cmake_dir / "FindFFTW.cmake")
        return
    old_block = """#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/findFFTW")
configure_file(downloadFindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
    else()
    message("CMake step for findFFTW completed (${result}).")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
endif()
set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")
"""
    new_block = """set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
"""
    if old_block not in text:
        raise RuntimeError(f"Unable to patch FFTW lookup in {cmake_lists}")
    cmake_dir = worktree / "cmake"
    _copy_cached_fftw(cmake_dir / "FindFFTW.cmake")
    cmake_lists.write_text(text.replace(old_block, new_block), encoding="utf-8")


def _build_in_worktree(worktree: Path, build_type: str) -> Path:
    build_dir = worktree / "build" / build_type.lower()
    build_script = worktree / "tools" / "build_executable.sh"
    if not build_script.is_file():
        raise RuntimeError(f"Missing build script in worktree: {build_script}")
    _prepare_offline_fftw(worktree, build_dir)
    subprocess.run(["bash", str(build_script), BUILD_SYSTEM, build_type.lower()], cwd=worktree, check=True)
    exe_path = build_dir / "u-dales"
    if not exe_path.is_file():
        raise RuntimeError(f"Build completed without producing executable: {exe_path}")
    return build_dir


def _build_current_workspace(build_type: str) -> Path:
    return _build_in_worktree(REPO_ROOT, build_type)


def _workspace_build_dir(build_type: str) -> Path:
    return REPO_ROOT / "build" / build_type.lower()


def _validate_build_dir(build_dir: Path, label: str) -> Path:
    resolved = build_dir.resolve()
    exe_path = resolved / "u-dales"
    if not exe_path.is_file():
        raise RuntimeError(f"{label} build dir does not contain executable: {exe_path}")
    return resolved


def _create_temp_root(workdir: Path, prefix: str, cleanup: bool):
    workdir.mkdir(parents=True, exist_ok=True)
    if cleanup:
        return tempfile.TemporaryDirectory(prefix=prefix, dir=workdir)

    path = Path(tempfile.mkdtemp(prefix=prefix, dir=workdir))

    class _KeepTempRoot:
        def __init__(self, keep_path: Path):
            self.name = str(keep_path)

        def __enter__(self) -> str:
            return self.name

        def __exit__(self, exc_type, exc, tb) -> bool:
            return False

    return _KeepTempRoot(path)


def _assert_base_state_logged(log_path: Path) -> None:
    text = log_path.read_text(encoding="utf-8", errors="replace")
    if BASE_STATE_LOG_MARKER not in text:
        raise RuntimeError(
            f"Startup log {log_path} does not contain '{BASE_STATE_LOG_MARKER}' "
            "(expected from modbasestate.f90's initbasestate)"
        )


def _assert_no_sentinel_or_nan(run_dir: Path, spec: CaseSpec) -> None:
    files = _dump_files(run_dir, spec)
    if not files:
        raise RuntimeError(f"No files matching '{spec.nc_pattern}' found in {run_dir}")

    failures: List[str] = []
    for path in files:
        with nc.Dataset(path) as ds:
            field_names = list(spec.fields) if spec.fields is not None else _discover_fields(path)
            for field in field_names:
                if field not in ds.variables:
                    failures.append(f"{path.name}: field '{field}' not present")
                    continue
                arr = np.asarray(ds.variables[field][:], dtype=np.float64)
                nan_mask = np.isnan(arr)
                if nan_mask.any():
                    idx = tuple(int(v) for v in np.unravel_index(int(np.argmax(nan_mask)), arr.shape))
                    failures.append(f"{path.name}: field '{field}' contains NaN at index {idx}")
                sentinel_mask = np.isclose(arr, SENTINEL_VALUE, atol=SENTINEL_ATOL, rtol=0.0)
                if sentinel_mask.any():
                    idx = tuple(int(v) for v in np.unravel_index(int(np.argmax(sentinel_mask)), arr.shape))
                    failures.append(
                        f"{path.name}: field '{field}' contains sentinel {SENTINEL_VALUE} at index {idx}"
                    )
    if failures:
        raise RuntimeError(
            "Dumped fields contain sentinel/NaN values (buried-slab base-state fallback not applied):\n- "
            + "\n- ".join(failures)
        )


def _run_assert_only_case(
    current_exe: Path,
    case_source: Path,
    tmp_root: Path,
    spec: CaseSpec,
) -> None:
    run_dir = tmp_root / "assert-only" / spec.case
    shutil.copytree(case_source, run_dir)
    namelist = run_dir / f"namoptions.{spec.case}"
    nprocx = _read_int_setting(namelist, "nprocx")
    nprocy = _read_int_setting(namelist, "nprocy")

    context = f"case {spec.case} [assert-only]"
    print(f"Running {context}", flush=True)
    print(f"  run dir: {run_dir}", flush=True)
    log_path = _run_case(current_exe, run_dir, spec, nprocx * nprocy)
    print(f"  run complete: {log_path}", flush=True)

    print(f"  checking startup log for '{BASE_STATE_LOG_MARKER}' ({context})", flush=True)
    _assert_base_state_logged(log_path)
    print(f"  base-state log check passed ({context})", flush=True)

    print(f"  scanning dumped fields for sentinel/NaN values ({context})", flush=True)
    _assert_no_sentinel_or_nan(run_dir, spec)
    print(f"  sentinel/NaN check passed ({context})", flush=True)


def _run_case_matrix(
    reference_exe: Path,
    current_exe: Path,
    case_source: Path,
    tmp_root: Path,
    spec: CaseSpec,
    compare_logs: bool,
    configs: Sequence[str],
) -> None:
    for label in configs:
        nprocx = CONFIGS[label]["nprocx"]
        nprocy = CONFIGS[label]["nprocy"]
        run_reference = tmp_root / "reference" / spec.case / label
        run_current = tmp_root / "current" / spec.case / label
        shutil.copytree(case_source, run_reference)
        shutil.copytree(case_source, run_current)
        _configure_namelist(run_reference, spec, nprocx, nprocy)
        _configure_namelist(run_current, spec, nprocx, nprocy)

        context = f"case {spec.case} [{label}]"
        print(f"Running {context}", flush=True)
        print(f"  reference run dir: {run_reference}", flush=True)
        reference_log = _run_case(reference_exe, run_reference, spec, nprocx * nprocy)
        print(f"  reference complete: {reference_log}", flush=True)
        print(f"  current run dir:   {run_current}", flush=True)
        current_log = _run_case(current_exe, run_current, spec, nprocx * nprocy)
        print(f"  current complete:   {current_log}", flush=True)
        print(f"  comparing fields for {context}", flush=True)
        target_time = _select_common_time(run_reference, run_current, spec)
        print(f"  comparing common time {target_time:.6f} for {context}", flush=True)
        _compare_fields(
            spec,
            _load_global_fields(run_reference, spec, target_time),
            _load_global_fields(run_current, spec, target_time),
            context,
        )
        print(f"  field comparison passed for {context}", flush=True)
        if compare_logs:
            print(f"  comparing logs for {context}", flush=True)
            _compare_logs(reference_log, current_log, context)
            print(f"  log comparison passed for {context}", flush=True)


def main(
    branch_a: str,
    branch_b: str,
    build_type: str,
    atol: Optional[float] = None,
    ci_mode: bool = False,
    reference_build_dir: Optional[Path] = None,
    reuse_current_build: bool = False,
    compare_logs: bool = False,
    workdir: Path = DEFAULT_WORKDIR,
    cleanup: bool = False,
    configs: Optional[Sequence[str]] = None,
) -> int:
    workdir = workdir.resolve()
    # An explicit --atol overrides every case's tolerance; otherwise each
    # case uses its own default_atol (see the CASES table and the module
    # docstring for why 090 and 092 differ).
    selected_specs = [
        spec._replace(abs_tol=(atol if atol is not None else spec.default_atol))
        for spec in CASES
    ]
    selected_configs = list(CONFIGS) if configs is None else list(configs)

    with _create_temp_root(workdir, "bc-cleanup-worktrees-", cleanup) as worktree_tmp, \
        _create_temp_root(workdir, "bc-cleanup-runs-", cleanup) as run_tmp:
        worktree_root = Path(worktree_tmp)
        run_root = Path(run_tmp)
        reference_worktree = worktree_root / "reference"
        current_worktree = worktree_root / "current"

        try:
            if reference_build_dir is not None:
                reference_exe = _validate_build_dir(reference_build_dir, "reference")
            else:
                reference_worktree = _create_worktree(worktree_root, "reference", branch_a)
                reference_exe = _build_in_worktree(reference_worktree, build_type)
            if ci_mode:
                current_worktree = _create_worktree(worktree_root, "current", branch_b)
                current_exe = _build_in_worktree(current_worktree, build_type)
            else:
                if reuse_current_build:
                    current_exe = _validate_build_dir(_workspace_build_dir(build_type), "current")
                else:
                    current_exe = _build_current_workspace(build_type)

            print("Running bc_cleanup regression", flush=True)
            print(f"reference branch:  {branch_a}", flush=True)
            print(f"current branch:    {branch_b}", flush=True)
            print(f"ci mode:           {ci_mode}", flush=True)
            print(f"build type:        {build_type}", flush=True)
            print(f"atol override:     {atol if atol is not None else '(none; using per-case default_atol)'}", flush=True)
            print(f"workdir:           {workdir}", flush=True)
            print(f"cleanup:           {cleanup}", flush=True)
            print(f"run root:          {run_root}", flush=True)
            print(f"cases:             {[(spec.case, spec.abs_tol) for spec in selected_specs]}", flush=True)
            print(f"configs:           {selected_configs}", flush=True)
            if reference_build_dir is not None:
                print(f"reference build dir:{reference_exe}", flush=True)
            else:
                print(f"reference worktree:{reference_worktree}", flush=True)
            if ci_mode:
                print(f"current worktree:  {current_worktree}", flush=True)
            else:
                print(f"current build dir: {current_exe}", flush=True)

            for spec in selected_specs:
                case_source = REPO_ROOT / "tests" / "cases" / spec.case
                if not case_source.is_dir():
                    raise RuntimeError(f"Missing case source: {case_source}")
                _run_case_matrix(reference_exe, current_exe, case_source, run_root, spec, compare_logs, selected_configs)
        finally:
            _remove_worktree(current_worktree)
            _remove_worktree(reference_worktree)
            subprocess.run([GIT, "worktree", "prune"], cwd=REPO_ROOT, check=True)

    print("bc_cleanup regression passed")
    return 0


def main_assert_only(
    case: str,
    build_type: str,
    reuse_current_build: bool = False,
    workdir: Path = DEFAULT_WORKDIR,
    cleanup: bool = False,
) -> int:
    """Single-branch smoke check: build the current workspace, run one case,
    and assert (a) the run completes, (b) the startup log contains
    'Base state:', and (c) no dumped field contains a -999. sentinel or NaN.

    Unlike main()/_run_case_matrix, this does not build or run a reference
    branch and does not compare fields against a reference — see
    .superpowers/sdd/task-7-brief.md step 2.
    """
    workdir = workdir.resolve()
    spec = ASSERT_ONLY_CASES.get(case)
    if spec is None:
        raise SystemExit(
            f"Unknown --assert-only case '{case}'; known cases: {sorted(ASSERT_ONLY_CASES)}"
        )

    with _create_temp_root(workdir, "bc-cleanup-assert-runs-", cleanup) as run_tmp:
        run_root = Path(run_tmp)

        if reuse_current_build:
            current_exe = _validate_build_dir(_workspace_build_dir(build_type), "current")
        else:
            current_exe = _build_current_workspace(build_type)

        print("Running bc_cleanup assert-only check", flush=True)
        print(f"case:              {case}", flush=True)
        print(f"build type:        {build_type}", flush=True)
        print(f"workdir:           {workdir}", flush=True)
        print(f"cleanup:           {cleanup}", flush=True)
        print(f"run root:          {run_root}", flush=True)
        print(f"current build dir: {current_exe}", flush=True)

        case_source = REPO_ROOT / "tests" / "cases" / spec.case
        if not case_source.is_dir():
            raise RuntimeError(f"Missing case source: {case_source}")
        _run_assert_only_case(current_exe, case_source, run_root, spec)

    print("bc_cleanup assert-only check passed")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("branch_a")
    parser.add_argument("branch_b")
    parser.add_argument("build_type", choices=["Debug", "Release"])
    parser.add_argument("--atol", type=float, default=None,
                        help="max allowed |current - reference| per field, applied to every "
                             "case. If omitted, each case uses its own default_atol from CASES "
                             "(090: 0.0 bitwise; 092: 5e-3).")
    parser.add_argument("--ci", action="store_true")
    parser.add_argument(
        "--reference-build-dir",
        type=Path,
        help="Reuse an existing reference build directory containing u-dales instead of rebuilding branch_a.",
    )
    parser.add_argument(
        "--reuse-current-build",
        action="store_true",
        help="Reuse the existing workspace build/<type>/u-dales instead of rebuilding the current workspace.",
    )
    parser.add_argument(
        "--compare-logs",
        action="store_true",
        help="Also compare solver log diagnostics (Courant, diffusion, divergence) in addition to field outputs.",
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=DEFAULT_WORKDIR,
        help="Directory under which reference worktrees and run outputs are created. Defaults to $UDALES_SCRATCH_ROOT or /tmp.",
    )
    parser.add_argument(
        "--cleanup",
        action="store_true",
        help="Remove temporary worktree and run directories at the end. Disabled by default so failed runs keep their artifacts.",
    )
    parser.add_argument(
        "--configs",
        type=str,
        help="Comma-separated decomposition labels to run, e.g. xy_split or serial,xy_split.",
    )
    parser.add_argument(
        "--assert-only",
        metavar="CASE",
        type=str,
        default=None,
        help=(
            "Run a single case (e.g. 091) against the current workspace build only, with no "
            "reference branch/build and no field comparison. Asserts the run completes, the "
            f"startup log contains '{BASE_STATE_LOG_MARKER}', and no dumped field contains a "
            f"{SENTINEL_VALUE} sentinel or NaN. branch_a/branch_b are accepted but unused. "
            f"Known cases: {sorted(ASSERT_ONLY_CASES)}."
        ),
    )
    args = parser.parse_args()
    configs = None
    if args.configs:
        configs = [item.strip() for item in args.configs.split(",") if item.strip()]
        unknown_configs = [label for label in configs if label not in CONFIGS]
        if unknown_configs:
            raise SystemExit(f"Unknown configs requested: {unknown_configs}")

    if args.assert_only:
        raise SystemExit(
            main_assert_only(
                args.assert_only,
                args.build_type,
                reuse_current_build=args.reuse_current_build,
                workdir=args.workdir,
                cleanup=args.cleanup,
            )
        )

    raise SystemExit(
        main(
            args.branch_a,
            args.branch_b,
            args.build_type,
            args.atol,
            ci_mode=args.ci,
            reference_build_dir=args.reference_build_dir,
            reuse_current_build=args.reuse_current_build,
            compare_logs=args.compare_logs,
            workdir=args.workdir,
            cleanup=args.cleanup,
            configs=configs,
        )
    )
