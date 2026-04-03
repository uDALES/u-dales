#!/usr/bin/env python3

from __future__ import annotations

"""Experimental legacy tree regression against the reference branch."""

import argparse
import os
import re
from shlex import quote as shlex_quote
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple

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
SCRATCH_ROOT = Path(
    os.environ.get("UDALES_SCRATCH_ROOT", "/tmp")
).resolve()

CASE_ID = 526
THERMO_FIELDS = ("tr_qt", "tr_thl")
MOMENTUM_FIELDS = ("tr_u", "tr_v", "tr_w")
INTERIOR_REL_TOL = 1.0e-3
NPROCS = 4
NPROCX = 2
NPROCY = 2
CACHED_FFTW_MODULE = REPO_ROOT / "build" / "debug" / "findFFTW-src" / "FindFFTW.cmake"


def _read_int_setting(namelist: Path, key: str) -> int:
    match = re.search(rf"(?m)^\s*{re.escape(key)}\s*=\s*(\d+)\s*$", namelist.read_text(encoding="utf-8"))
    if not match:
        raise RuntimeError(f"Missing integer setting '{key}' in {namelist}")
    return int(match.group(1))

def _patch_namelist(run_dir: Path, mode: str) -> None:
    namelist = run_dir / f"namoptions.{CASE_ID}"
    text = namelist.read_text(encoding="utf-8")
    text = re.sub(r"(?m)^(\s*nprocx\s*=\s*)\d+(\s*)$", rf"\g<1>{NPROCX}\2", text)
    text = re.sub(r"(?m)^(\s*nprocy\s*=\s*)\d+(\s*)$", rf"\g<1>{NPROCY}\2", text)

    if mode == "reference":
        text = re.sub(r"(?m)^\s*itree_mode\s*=.*\n?", "", text)
    elif mode == "current":
        if re.search(r"(?m)^\s*itree_mode\s*=", text):
            text = re.sub(r"(?m)^(\s*itree_mode\s*=\s*)\d+(\s*)$", r"\g<1>99\2", text)
        else:
            text = text.replace("&TREES\n", "&TREES\nitree_mode   = 99\n", 1)
    else:
        raise ValueError(f"Unknown mode: {mode}")

    namelist.write_text(text, encoding="utf-8")


def _run_case(path_to_exe: Path, run_dir: Path) -> None:
    launcher = os.environ.get("MPIEXEC", "mpiexec")
    exe = path_to_exe / "u-dales"
    diag_log = run_dir / "fortran_diagnostics.log"
    for stale in run_dir.glob("monitor*.txt"):
        stale.unlink()
    for stale in run_dir.glob(f"treedump.*.*.{CASE_ID}.nc"):
        stale.unlink()
    if diag_log.exists():
        diag_log.unlink()
    shell_command = (
        f"module load {RUNTIME_MODULES} && "
        f"export HDF5_USE_FILE_LOCKING=FALSE && "
        f"export FOR_DISABLE_DIAGNOSTIC_DISPLAY=TRUE && "
        f"export FOR_DIAGNOSTIC_LOG_FILE={shlex_quote(str(diag_log))} && "
        f"cd {shlex_quote(str(run_dir))} && "
        f"{shlex_quote(launcher)} -n {NPROCS} {shlex_quote(str(exe))} namoptions.{CASE_ID}"
    )
    print(shell_command)
    subprocess.run(["bash", "-lc", shell_command], cwd=REPO_ROOT, check=True)


def _treedump_files(run_dir: Path) -> List[Path]:
    files = sorted(run_dir.glob(f"treedump.*.*.{CASE_ID}.nc"))
    if not files:
        raise RuntimeError(f"No treedump outputs found in {run_dir}")
    return files


def _copy_for_local_read(path_a: Path, path_b: Path) -> Tuple[Path, Path]:
    scratch_dir = SCRATCH_ROOT / "udales-tree-regression-compare"
    scratch_dir.mkdir(parents=True, exist_ok=True)
    local_a = scratch_dir / "reference" / path_a.name
    local_b = scratch_dir / "current" / path_b.name
    local_a.parent.mkdir(parents=True, exist_ok=True)
    local_b.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path_a, local_a)
    shutil.copy2(path_b, local_b)
    return local_a, local_b


def _load_veg_mask(case_dir: Path) -> np.ndarray:
    namelist = case_dir / f"namoptions.{CASE_ID}"
    itot = _read_int_setting(namelist, "itot")
    jtot = _read_int_setting(namelist, "jtot")
    ktot = _read_int_setting(namelist, "ktot")
    mask = np.zeros((ktot, jtot, itot), dtype=bool)
    with (case_dir / f"veg.inp.{CASE_ID}").open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            i, j, k = (int(part) for part in stripped.split()[:3])
            if not (1 <= i <= itot and 1 <= j <= jtot and 1 <= k <= ktot):
                raise RuntimeError(
                    f"Vegetation point ({i}, {j}, {k}) is outside "
                    f"1..{itot}, 1..{jtot}, 1..{ktot}"
                )
            mask[k - 1, j - 1, i - 1] = True
    return mask


def _tile_mask(global_mask: np.ndarray, case_dir: Path, tile_name: str, tile_shape: Tuple[int, int, int]) -> np.ndarray:
    namelist = case_dir / f"namoptions.{CASE_ID}"
    nprocx = _read_int_setting(namelist, "nprocx")
    nprocy = _read_int_setting(namelist, "nprocy")
    itot = _read_int_setting(namelist, "itot")
    jtot = _read_int_setting(namelist, "jtot")
    match = re.match(rf"treedump\.(\d{{3}})\.(\d{{3}})\.{CASE_ID}\.nc$", tile_name)
    if not match:
        raise RuntimeError(f"Unexpected treedump tile name: {tile_name}")
    px = int(match.group(1))
    py = int(match.group(2))
    nx = tile_shape[2]
    ny = tile_shape[1]
    nz = tile_shape[0]
    if px >= nprocx or py >= nprocy:
        raise RuntimeError(f"Tile index out of range for {tile_name}")
    istart = px * (itot // nprocx)
    jstart = py * (jtot // nprocy)
    return global_mask[:nz, jstart : jstart + ny, istart : istart + nx]


def _fortran_local_idx(idx: Tuple[int, int, int]) -> Tuple[int, int, int]:
    z, y, x = idx
    return (x + 1, y + 1, z + 1)


def _format_loc(tile_name: str, idx: Tuple[int, int, int]) -> str:
    i, j, k = _fortran_local_idx(idx)
    return f"{tile_name} local ({i}, {j}, {k})"


def _format_optional_loc(tile_name: str | None, idx: Tuple[int, int, int] | None) -> str:
    if tile_name is None or idx is None:
        return "n/a"
    return _format_loc(tile_name, idx)


def _safe_relative(diff_value: float, ref_value: float) -> float:
    if ref_value <= 0.0:
        return 0.0 if diff_value == 0.0 else float("inf")
    return diff_value / ref_value


def _compare_outputs(
    reference_files: List[Path],
    current_files: List[Path],
    case_dir: Path,
    enforce_tolerances: bool = True,
) -> None:
    if len(reference_files) != len(current_files):
        raise RuntimeError(
            "Legacy tree regression produced a different number of treedump tiles: "
            f"{len(reference_files)} vs {len(current_files)}"
        )

    global_mask = _load_veg_mask(case_dir)
    thermo_max: Dict[str, Tuple[float, str | None, Tuple[int, int, int] | None]] = {
        field: (-1.0, None, None) for field in THERMO_FIELDS
    }
    thermo_ref_max = {field: 0.0 for field in THERMO_FIELDS}
    momentum_max: Dict[str, Tuple[float, str | None, Tuple[int, int, int] | None]] = {
        field: (-1.0, None, None) for field in MOMENTUM_FIELDS
    }
    momentum_ref_max = {field: 0.0 for field in MOMENTUM_FIELDS}
    interior_max: Dict[str, Tuple[float, str | None, Tuple[int, int, int] | None]] = {
        field: (-1.0, None, None) for field in MOMENTUM_FIELDS
    }
    interior_ref_max = {field: 0.0 for field in MOMENTUM_FIELDS}
    side_max: Dict[str, Tuple[float, str | None, Tuple[int, int, int] | None]] = {
        field: (-1.0, None, None) for field in MOMENTUM_FIELDS
    }
    side_ref_max = {field: 0.0 for field in MOMENTUM_FIELDS}
    momentum_counts = {field: 0 for field in MOMENTUM_FIELDS}
    interior_counts = {field: 0 for field in MOMENTUM_FIELDS}
    side_counts = {field: 0 for field in MOMENTUM_FIELDS}

    for path_reference, path_current in zip(reference_files, current_files):
        if path_reference.name != path_current.name:
            raise RuntimeError(
                f"Mismatched treedump tiles: {path_reference.name} vs {path_current.name}"
            )
        local_reference, local_current = _copy_for_local_read(path_reference, path_current)
        with nc.Dataset(local_reference) as ds_reference, nc.Dataset(local_current) as ds_current:
            sample = np.asarray(ds_reference.variables["tr_qt"][-1], dtype=np.float64)
            tile_mask = _tile_mask(global_mask, case_dir, path_reference.name, sample.shape)

            for field in THERMO_FIELDS:
                arr_reference = np.asarray(ds_reference.variables[field][-1], dtype=np.float64)
                arr_current = np.asarray(ds_current.variables[field][-1], dtype=np.float64)
                diff = np.abs(arr_current - arr_reference)
                thermo_ref_max[field] = max(thermo_ref_max[field], float(np.max(np.abs(arr_reference))))
                idx = tuple(int(v) for v in np.unravel_index(np.argmax(diff), diff.shape))
                val = float(diff[idx])
                if val > thermo_max[field][0]:
                    thermo_max[field] = (val, path_reference.name, idx)

            for field in MOMENTUM_FIELDS:
                arr_reference = np.asarray(ds_reference.variables[field][-1], dtype=np.float64)
                arr_current = np.asarray(ds_current.variables[field][-1], dtype=np.float64)
                diff = np.abs(arr_current - arr_reference)
                momentum_ref_max[field] = max(momentum_ref_max[field], float(np.max(np.abs(arr_reference))))
                idx = tuple(int(v) for v in np.unravel_index(np.argmax(diff), diff.shape))
                val = float(diff[idx])
                if val > momentum_max[field][0]:
                    momentum_max[field] = (val, path_reference.name, idx)
                momentum_counts[field] += int(np.count_nonzero(arr_current != arr_reference))

                if field == "tr_u":
                    interior = np.zeros_like(tile_mask, dtype=bool)
                    interior[:, :, 1:] = tile_mask[:, :, :-1] & tile_mask[:, :, 1:]
                    support = np.zeros_like(tile_mask, dtype=bool)
                    support[:, :, 1:] = tile_mask[:, :, :-1] | tile_mask[:, :, 1:]
                elif field == "tr_v":
                    interior = np.zeros_like(tile_mask, dtype=bool)
                    interior[:, 1:, :] = tile_mask[:, :-1, :] & tile_mask[:, 1:, :]
                    support = np.zeros_like(tile_mask, dtype=bool)
                    support[:, 1:, :] = tile_mask[:, :-1, :] | tile_mask[:, 1:, :]
                else:
                    interior = np.zeros_like(tile_mask, dtype=bool)
                    interior[1:, :, :] = tile_mask[:-1, :, :] & tile_mask[1:, :, :]
                    support = np.zeros_like(tile_mask, dtype=bool)
                    support[1:, :, :] = tile_mask[:-1, :, :] | tile_mask[1:, :, :]
                side = support & ~interior

                interior_counts[field] += int(np.count_nonzero(interior))
                if np.any(interior):
                    interior_ref_max[field] = max(
                        interior_ref_max[field],
                        float(np.max(np.where(interior, np.abs(arr_reference), 0.0))),
                    )
                    interior_diff = np.where(interior, diff, 0.0)
                    interior_idx = tuple(int(v) for v in np.unravel_index(np.argmax(interior_diff), interior_diff.shape))
                    interior_val = float(interior_diff[interior_idx])
                    if interior_val > interior_max[field][0]:
                        interior_max[field] = (interior_val, path_reference.name, interior_idx)
                side_counts[field] += int(np.count_nonzero(side))
                if np.any(side):
                    side_ref_max[field] = max(
                        side_ref_max[field],
                        float(np.max(np.where(side, np.abs(arr_reference), 0.0))),
                    )
                    side_diff = np.where(side, diff, 0.0)
                    side_idx = tuple(int(v) for v in np.unravel_index(np.argmax(side_diff), side_diff.shape))
                    side_val = float(side_diff[side_idx])
                    if side_val > side_max[field][0]:
                        side_max[field] = (side_val, path_reference.name, side_idx)

    for field in THERMO_FIELDS:
        value, tile_name, idx = thermo_max[field]
        rel = _safe_relative(value, thermo_ref_max[field])
        print(
            f"Relative diff {field}: {rel:.6e} "
            f"at {_format_optional_loc(tile_name, idx)}"
        )

    interior_failures = []
    for field in MOMENTUM_FIELDS:
        value, tile_name, idx = momentum_max[field]
        int_value, int_tile_name, int_idx = interior_max[field]
        side_value, side_tile_name, side_idx = side_max[field]
        rel = _safe_relative(value, momentum_ref_max[field])
        int_rel = _safe_relative(int_value, interior_ref_max[field])
        side_rel = _safe_relative(side_value, side_ref_max[field])
        print(
            f"Relative diff {field}: {rel:.6e} "
            f"at {_format_optional_loc(tile_name, idx)}"
        )
        print(f"Count diff {field}: {momentum_counts[field]}")
        print(
            f"Relative diff (interior) {field}: {int_rel:.6e} "
            f"at {_format_optional_loc(int_tile_name, int_idx)} count {interior_counts[field]}"
        )
        print(
            f"Relative diff (side) {field}: {side_rel:.6e} "
            f"at {_format_optional_loc(side_tile_name, side_idx)} count {side_counts[field]}"
        )
        if enforce_tolerances and int_rel > INTERIOR_REL_TOL:
            interior_failures.append(
                f"{field} interior relative diff {int_rel:.6e} at "
                f"{_format_optional_loc(int_tile_name, int_idx)}"
            )

    if interior_failures:
        raise RuntimeError(
            "Legacy tree regression interior relative differences exceed "
            f"{INTERIOR_REL_TOL:.1e}:\n- " + "\n- ".join(interior_failures)
        )


def analyze_existing_outputs(
    reference_run_dir: Path,
    current_run_dir: Path,
    case_dir: Path,
    enforce_tolerances: bool = True,
) -> None:
    _compare_outputs(
        _treedump_files(reference_run_dir),
        _treedump_files(current_run_dir),
        case_dir,
        enforce_tolerances=enforce_tolerances,
    )


def _sanitize_ref(ref: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", ref)


def _create_worktree(worktree_root: Path, label: str, ref: str) -> Path:
    path = worktree_root / label
    try:
        subprocess.run(
            [GIT, "worktree", "add", "--detach", str(path), ref],
            cwd=REPO_ROOT,
            check=True,
        )
    except subprocess.CalledProcessError:
        subprocess.run(
            [GIT, "clone", "--shared", "--no-checkout", str(REPO_ROOT), str(path)],
            cwd=worktree_root,
            check=True,
        )
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
            raise RuntimeError(
                f"Required local submodule directory is missing from the main checkout: {source}"
            )
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
    if 'set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)' in text:
        findfftw_dir = build_dir / "findFFTW-src"
        findfftw_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(CACHED_FFTW_MODULE, findfftw_dir / "FindFFTW.cmake")
        return

    if 'set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")' in text:
        cmake_dir = worktree / "cmake"
        cmake_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(CACHED_FFTW_MODULE, cmake_dir / "FindFFTW.cmake")
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
    cmake_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(CACHED_FFTW_MODULE, cmake_dir / "FindFFTW.cmake")
    cmake_lists.write_text(text.replace(old_block, new_block), encoding="utf-8")


def _build_in_worktree(worktree: Path, ref: str, build_type: str) -> Path:
    build_dir = worktree / "build" / build_type.lower()
    build_script = worktree / "tools" / "build_executable.sh"
    if not build_script.is_file():
        raise RuntimeError(f"Missing build script in worktree: {build_script}")
    _prepare_offline_fftw(worktree, build_dir)
    subprocess.run(
        ["bash", str(build_script), BUILD_SYSTEM, build_type.lower()],
        cwd=worktree,
        check=True,
    )
    exe_path = build_dir / "u-dales"
    if not exe_path.is_file():
        raise RuntimeError(f"Build completed without producing executable: {exe_path}")
    return build_dir


def _build_current_workspace(build_type: str) -> Path:
    return _build_in_worktree(REPO_ROOT, "workspace", build_type)


def main(branch_a: str, branch_b: str, build_type: str, ci_mode: bool = False) -> int:
    case_source = REPO_ROOT / "tests" / "cases" / str(CASE_ID)
    SCRATCH_ROOT.mkdir(parents=True, exist_ok=True)
    shutil.rmtree(SCRATCH_ROOT / "udales-tree-regression-compare", ignore_errors=True)
    with tempfile.TemporaryDirectory(prefix="tree-regression-worktrees-", dir=SCRATCH_ROOT) as worktree_tmp, tempfile.TemporaryDirectory(prefix="tree-reference-", dir=SCRATCH_ROOT) as reference_tmp, tempfile.TemporaryDirectory(prefix="tree-current-", dir=SCRATCH_ROOT) as current_tmp:
        worktree_root = Path(worktree_tmp)
        reference_worktree = worktree_root / "reference"
        current_worktree = worktree_root / "current"

        try:
            reference_worktree = _create_worktree(worktree_root, "reference", branch_a)
            reference_exe = _build_in_worktree(reference_worktree, branch_a, build_type)
            if ci_mode:
                current_worktree = _create_worktree(worktree_root, "current", branch_b)
                current_exe = _build_in_worktree(current_worktree, branch_b, build_type)
            else:
                current_exe = _build_current_workspace(build_type)

            run_reference = Path(reference_tmp)
            run_current = Path(current_tmp)
            shutil.copytree(case_source, run_reference, dirs_exist_ok=True)
            shutil.copytree(case_source, run_current, dirs_exist_ok=True)
            _patch_namelist(run_reference, "reference")
            _patch_namelist(run_current, "current")

            print("Running legacy tree regression against reference branch")
            print(f"reference branch:    {branch_a}")
            print(f"current branch:   {branch_b}")
            print(f"ci mode:          {ci_mode}")
            print(f"build type:       {build_type}")
            print(f"scratch root:        {SCRATCH_ROOT}")
            print(f"reference worktree:  {reference_worktree}")
            if ci_mode:
                print(f"current worktree:    {current_worktree}")
            else:
                print(f"current build dir:   {current_exe}")
            print(f"reference run dir:   {run_reference}")
            print(f"current run dir:     {run_current}")

            _run_case(reference_exe, run_reference)
            _run_case(current_exe, run_current)

            _compare_outputs(_treedump_files(run_reference), _treedump_files(run_current), case_source)
        finally:
            _remove_worktree(current_worktree)
            _remove_worktree(reference_worktree)
            subprocess.run([GIT, "worktree", "prune"], cwd=REPO_ROOT, check=True)

    print("Legacy tree regression against reference branch passed")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("branch_a", nargs="?")
    parser.add_argument("branch_b", nargs="?")
    parser.add_argument("build_type", nargs="?", choices=["Debug", "Release"])
    parser.add_argument("--analyze-only", action="store_true")
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--ci", action="store_true")
    parser.add_argument("--reference-run-dir", type=Path)
    parser.add_argument("--current-run-dir", type=Path)
    parser.add_argument("--case-dir", type=Path, default=REPO_ROOT / "tests" / "cases" / str(CASE_ID))
    args = parser.parse_args()

    if args.analyze_only:
        if args.reference_run_dir is None or args.current_run_dir is None:
            parser.error("--analyze-only requires --reference-run-dir and --current-run-dir")
        analyze_existing_outputs(
            args.reference_run_dir,
            args.current_run_dir,
            args.case_dir,
            enforce_tolerances=not args.report_only,
        )
        raise SystemExit(0)

    if args.branch_a is None or args.branch_b is None or args.build_type is None:
        parser.error("branch_a, branch_b, and build_type are required unless --analyze-only is used")

    raise SystemExit(main(args.branch_a, args.branch_b, args.build_type, ci_mode=args.ci))
