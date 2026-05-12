#!/usr/bin/env python3

"""Regression for MPI averaging/reduction behaviour against a reference branch."""

import argparse
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from shlex import quote as shlex_quote
from typing import Dict, List, Optional, Sequence, Tuple

import netCDF4 as nc
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REGRESSION_DIR = SCRIPT_DIR.parent
TESTS_DIR = REGRESSION_DIR.parent
REPO_ROOT = TESTS_DIR.parent
NAMELIST_DIR = SCRIPT_DIR / "namelists"

BUILD_SYSTEM = os.environ.get("UDALES_BUILD_SYSTEM", "icl")
GIT = os.environ.get("UDALES_GIT", "/usr/bin/git")
SUBMODULE_DIRS = ("2decomp-fft", "tools/View3D")
RUNTIME_MODULES = os.environ.get(
    "UDALES_RUNTIME_MODULES",
    "intel/2025a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a "
    "FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0",
)
DEFAULT_SCRATCH_ROOT = Path(os.environ.get("UDALES_SCRATCH_ROOT", "/tmp")).resolve()
CACHED_FFTW_MODULE = REPO_ROOT / "build" / "debug" / "findFFTW-src" / "FindFFTW.cmake"

ABS_TOL = 1.0e-9
NO_TREE_GLOBAL_ABS_TOL = 2.0e-8
LOG_ATOL = 1.0e-12
LOG_RTOL = 1.0e-10

CONFIGS = {
    "serial": (1, 1),
    "x_split": (2, 1),
    "y_split": (1, 2),
    "xy_split": (2, 2),
}

LOG_PATTERNS = {
    "courant": re.compile(r"Courant numbers \(x,y,z,tot\):\s*([Ee0-9+\-\.]+)"),
    "diffusion": re.compile(r"Diffusion number:\s*([Ee0-9+\-\.]+)"),
    "divergence": re.compile(r"divmax, divtot =\s*([Ee0-9+\-\.]+)\s+([Ee0-9+\-\.]+)"),
}


@dataclass(frozen=True)
class CaseSpec:
    case_id: int
    source_dir: Path
    prefix: str
    fields: Tuple[str, ...]
    vegetation_enabled: bool
    abs_tol: float

    @property
    def case_name(self) -> str:
        return f"{self.case_id:03d}"


CASES = (
    CaseSpec(
        case_id=100,
        source_dir=REPO_ROOT / "tests" / "cases" / "100",
        prefix="tdump",
        fields=("ut", "vt", "wt"),
        vegetation_enabled=False,
        abs_tol=NO_TREE_GLOBAL_ABS_TOL,
    ),
    CaseSpec(
        case_id=526,
        source_dir=REPO_ROOT / "tests" / "cases" / "526",
        prefix="treedump",
        fields=("tr_u", "tr_v", "tr_w"),
        vegetation_enabled=True,
        abs_tol=ABS_TOL,
    ),
)

CASE_BY_ID = {spec.case_id: spec for spec in CASES}


def _read_int_setting(namelist: Path, key: str) -> int:
    match = re.search(
        rf"(?m)^\s*{re.escape(key)}\s*=\s*(\d+)\s*$",
        namelist.read_text(encoding="utf-8"),
    )
    if not match:
        raise RuntimeError(f"Missing integer setting '{key}' in {namelist}")
    return int(match.group(1))


def _namelist_path(spec: CaseSpec, label: str, mode: str) -> Path:
    if spec.case_id == 100:
        return NAMELIST_DIR / "namoptions.100.serial"
    if spec.case_id == 526:
        return NAMELIST_DIR / f"namoptions.526.serial.{mode}"
    raise ValueError(f"Unsupported case id for namelist templates: {spec.case_id}")


def _copy_namelist(run_dir: Path, spec: CaseSpec, label: str, mode: str) -> None:
    namelist_source = _namelist_path(spec, label, mode)
    if not namelist_source.is_file():
        raise RuntimeError(f"Namelist file not found: {namelist_source}")
    target = run_dir / f"namoptions.{spec.case_id}"
    shutil.copy2(namelist_source, target)
    nprocx, nprocy = CONFIGS[label]
    text = target.read_text(encoding="utf-8")
    text = re.sub(r"(?m)^(\s*nprocx\s*=\s*)\d+(\s*)$", rf"\g<1>{nprocx}\2", text)
    text = re.sub(r"(?m)^(\s*nprocy\s*=\s*)\d+(\s*)$", rf"\g<1>{nprocy}\2", text)
    target.write_text(text, encoding="utf-8")


def _load_veg_mask(case_dir: Path, case_id: int) -> np.ndarray:
    namelist = case_dir / f"namoptions.{case_id}"
    itot = _read_int_setting(namelist, "itot")
    jtot = _read_int_setting(namelist, "jtot")
    ktot = _read_int_setting(namelist, "ktot")
    mask = np.zeros((ktot, jtot, itot), dtype=bool)
    with (case_dir / f"veg.inp.{case_id}").open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            i, j, k = (int(part) for part in stripped.split()[:3])
            mask[k - 1, j - 1, i - 1] = True
    return mask


def _support_mask(field: str, veg_mask: np.ndarray) -> np.ndarray:
    support = np.zeros_like(veg_mask, dtype=bool)
    if field == "tr_u":
        support[:, :, 1:] = veg_mask[:, :, :-1] | veg_mask[:, :, 1:]
    elif field == "tr_v":
        support[:, 1:, :] = veg_mask[:, :-1, :] | veg_mask[:, 1:, :]
    elif field == "tr_w":
        support[1:, :, :] = veg_mask[:-1, :, :] | veg_mask[1:, :, :]
    else:
        raise ValueError(f"Unknown field: {field}")
    return support


def _parse_tile_indices(path: Path, prefix: str, case_id: int) -> Tuple[int, int]:
    match = re.match(rf"{re.escape(prefix)}\.(\d{{3}})\.(\d{{3}})\.{case_id}\.nc$", path.name)
    if not match:
        raise RuntimeError(f"Unexpected output file name: {path.name}")
    return int(match.group(1)), int(match.group(2))


def _load_available_times(run_dir: Path, spec: CaseSpec) -> np.ndarray:
    files = sorted(run_dir.glob(f"{spec.prefix}.*.*.{spec.case_id}.nc"))
    if not files:
        raise RuntimeError(f"No {spec.prefix} outputs found in {run_dir}")
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
        if np.any(np.isclose(cur_times, ref_time, atol=LOG_ATOL, rtol=LOG_RTOL)):
            common.append(float(ref_time))
    if not common:
        raise RuntimeError(
            f"No common output times between {reference_run_dir} ({ref_times}) and "
            f"{current_run_dir} ({cur_times})"
        )
    return max(common)


def _load_global_fields(run_dir: Path, spec: CaseSpec, target_time: float) -> Dict[str, np.ndarray]:
    namelist = run_dir / f"namoptions.{spec.case_id}"
    itot = _read_int_setting(namelist, "itot")
    jtot = _read_int_setting(namelist, "jtot")
    ktot = _read_int_setting(namelist, "ktot")
    nprocx = _read_int_setting(namelist, "nprocx")
    nprocy = _read_int_setting(namelist, "nprocy")

    files = sorted(run_dir.glob(f"{spec.prefix}.*.*.{spec.case_id}.nc"))
    if not files:
        raise RuntimeError(f"No {spec.prefix} outputs found in {run_dir}")

    fields = {
        field: np.zeros((ktot, jtot, itot), dtype=np.float64) for field in spec.fields
    }
    for path in files:
        px, py = _parse_tile_indices(path, spec.prefix, spec.case_id)
        with nc.Dataset(path) as ds:
            if "time" not in ds.variables:
                raise RuntimeError(f"Missing time variable in {path}")
            times = np.asarray(ds.variables["time"][:], dtype=np.float64)
            if times.size == 0:
                raise RuntimeError(f"No time records found in {path}")
            matches = np.where(np.isclose(times, target_time, atol=LOG_ATOL, rtol=LOG_RTOL))[0]
            if matches.size == 0:
                raise RuntimeError(f"No record at time {target_time} found in {path}; available times: {times}")
            record = int(matches[-1])
            sample = np.asarray(ds.variables[spec.fields[0]][record], dtype=np.float64)
            nz, ny, nx = sample.shape
            istart = px * (itot // nprocx)
            jstart = py * (jtot // nprocy)
            for field in spec.fields:
                arr = np.asarray(ds.variables[field][record], dtype=np.float64)
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


def _run_case(path_to_exe: Path, run_dir: Path, case_id: int, nprocs: int) -> Path:
    mpiexec, extra_args = _mpi_exec_and_args()
    exe = path_to_exe / "u-dales"
    log_path = run_dir / "solver.log"
    diag_log = run_dir / "fortran_diagnostics.log"

    for stale in run_dir.glob("monitor*.txt"):
        stale.unlink()
    for stale in run_dir.glob(f"treedump.*.*.{case_id}.nc"):
        stale.unlink()
    for stale in run_dir.glob(f"tdump.*.*.{case_id}.nc"):
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
        f"{shlex_quote(str(exe))} namoptions.{case_id}"
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
        f"u-dales run failed for case {case_id} in {run_dir}",
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
                if np.isclose(cur_val, ref_val, atol=LOG_ATOL, rtol=LOG_RTOL):
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
    veg_mask: Optional[np.ndarray],
) -> None:
    failures: List[str] = []
    for field in spec.fields:
        if spec.vegetation_enabled:
            assert veg_mask is not None
            mask = _support_mask(field, veg_mask)
        else:
            mask = np.ones(reference[field].shape, dtype=bool)
        diff = np.abs(current[field] - reference[field])
        masked = np.where(mask, diff, -1.0)
        flat = int(np.argmax(masked))
        idx = tuple(int(v) for v in np.unravel_index(flat, diff.shape))
        max_diff = 0.0 if not np.any(mask) else float(diff[idx])
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


def _create_temp_root(scratch_root: Path, prefix: str, cleanup: bool):
    scratch_root.mkdir(parents=True, exist_ok=True)
    if cleanup:
        return tempfile.TemporaryDirectory(prefix=prefix, dir=scratch_root)

    path = Path(tempfile.mkdtemp(prefix=prefix, dir=scratch_root))

    class _KeepTempRoot:
        def __init__(self, keep_path: Path):
            self.name = str(keep_path)

        def __enter__(self) -> str:
            return self.name

        def __exit__(self, exc_type, exc, tb) -> bool:
            return False

    return _KeepTempRoot(path)


def _run_case_matrix(
    reference_exe: Path,
    current_exe: Path,
    case_source: Path,
    tmp_root: Path,
    spec: CaseSpec,
    compare_logs: bool,
    configs: Sequence[str],
) -> None:
    veg_mask = _load_veg_mask(case_source, spec.case_id) if spec.vegetation_enabled else None

    for label in configs:
        nprocx, nprocy = CONFIGS[label]
        run_reference = tmp_root / "reference" / spec.case_name / label
        run_current = tmp_root / "current" / spec.case_name / label
        shutil.copytree(case_source, run_reference)
        shutil.copytree(case_source, run_current)
        _copy_namelist(run_reference, spec, label, "reference")
        _copy_namelist(run_current, spec, label, "current")

        context = f"case {spec.case_name} [{label}]"
        print(f"Running {context}", flush=True)
        print(f"  reference run dir: {run_reference}", flush=True)
        reference_log = _run_case(reference_exe, run_reference, spec.case_id, nprocx * nprocy)
        print(f"  reference complete: {reference_log}", flush=True)
        print(f"  current run dir:   {run_current}", flush=True)
        current_log = _run_case(current_exe, run_current, spec.case_id, nprocx * nprocy)
        print(f"  current complete:   {current_log}", flush=True)
        print(f"  comparing fields for {context}", flush=True)
        target_time = _select_common_time(run_reference, run_current, spec)
        print(f"  comparing common time {target_time:.6f} for {context}", flush=True)
        _compare_fields(
            spec,
            _load_global_fields(run_reference, spec, target_time),
            _load_global_fields(run_current, spec, target_time),
            context,
            veg_mask,
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
    ci_mode: bool = False,
    reference_build_dir: Optional[Path] = None,
    reuse_current_build: bool = False,
    compare_logs: bool = False,
    scratch_root: Path = DEFAULT_SCRATCH_ROOT,
    cleanup: bool = False,
    case_ids: Optional[Sequence[int]] = None,
    configs: Optional[Sequence[str]] = None,
) -> int:
    scratch_root = scratch_root.resolve()
    selected_specs = list(CASES) if case_ids is None else [CASE_BY_ID[case_id] for case_id in case_ids]
    selected_configs = list(CONFIGS) if configs is None else list(configs)

    with _create_temp_root(scratch_root, "mpi-averaging-worktrees-", cleanup) as worktree_tmp, \
        _create_temp_root(scratch_root, "mpi-averaging-runs-", cleanup) as run_tmp:
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

            print("Running MPI averaging regression", flush=True)
            print(f"reference branch:  {branch_a}", flush=True)
            print(f"current branch:    {branch_b}", flush=True)
            print(f"ci mode:           {ci_mode}", flush=True)
            print(f"build type:        {build_type}", flush=True)
            print(f"scratch root:      {scratch_root}", flush=True)
            print(f"cleanup:           {cleanup}", flush=True)
            print(f"run root:          {run_root}", flush=True)
            print(f"cases:             {[spec.case_id for spec in selected_specs]}", flush=True)
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
                case_source = spec.source_dir
                if not case_source.is_dir():
                    raise RuntimeError(f"Missing case source: {case_source}")
                _run_case_matrix(reference_exe, current_exe, case_source, run_root, spec, compare_logs, selected_configs)
        finally:
            _remove_worktree(current_worktree)
            _remove_worktree(reference_worktree)
            subprocess.run([GIT, "worktree", "prune"], cwd=REPO_ROOT, check=True)

    print("MPI averaging regression passed")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("branch_a")
    parser.add_argument("branch_b")
    parser.add_argument("build_type", choices=["Debug", "Release"])
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
        "--scratch-root",
        type=Path,
        default=DEFAULT_SCRATCH_ROOT,
        help="Directory under which reference worktrees and run outputs are created. Defaults to $UDALES_SCRATCH_ROOT or /tmp.",
    )
    parser.add_argument(
        "--cleanup",
        action="store_true",
        help="Remove temporary worktree and run directories at the end. Disabled by default so failed runs keep their artifacts.",
    )
    parser.add_argument(
        "--cases",
        type=str,
        help="Comma-separated case ids to run, e.g. 100,526. Defaults to all configured cases.",
    )
    parser.add_argument(
        "--configs",
        type=str,
        help="Comma-separated decomposition labels to run, e.g. xy_split or serial,x_split.",
    )
    args = parser.parse_args()
    case_ids = None
    if args.cases:
        case_ids = [int(item.strip()) for item in args.cases.split(",") if item.strip()]
        unknown_case_ids = [case_id for case_id in case_ids if case_id not in CASE_BY_ID]
        if unknown_case_ids:
            raise SystemExit(f"Unknown case ids requested: {unknown_case_ids}")
    configs = None
    if args.configs:
        configs = [item.strip() for item in args.configs.split(",") if item.strip()]
        unknown_configs = [label for label in configs if label not in CONFIGS]
        if unknown_configs:
            raise SystemExit(f"Unknown configs requested: {unknown_configs}")

    raise SystemExit(
        main(
            args.branch_a,
            args.branch_b,
            args.build_type,
            ci_mode=args.ci,
            reference_build_dir=args.reference_build_dir,
            reuse_current_build=args.reuse_current_build,
            compare_logs=args.compare_logs,
            scratch_root=args.scratch_root,
            cleanup=args.cleanup,
            case_ids=case_ids,
            configs=configs,
        )
    )
