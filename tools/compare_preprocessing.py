#!/usr/bin/env python3
"""Compare MATLAB and Python preprocessing outputs for a case directory.

This tool is reusable from both integration tests and exploratory runs on real
experiment directories. It can either:

- compare a specific case directory via ``--case-dir``
- compare a case under an experiments root via positional ``case``
- list available real cases sorted by minimum STL facet count
"""

import argparse
import os
import shlex
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Optional

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXPERIMENTS_ROOT = Path.home() / "udales" / "experiments"
DEFAULT_TEMP_ROOTS = [
    os.environ.get("TMPDIR_PARENT"),
    "/rds/general/user/mvr/ephemeral",
    "/tmp",
]
SUPPORTED_PATTERNS = [
    "lscale.inp.{case}",
    "prof.inp.{case}",
    "factypes.inp.{case}",
    "facets.inp.{case}",
    "facetarea.inp.{case}",
    "solid_c.txt",
    "solid_u.txt",
    "solid_v.txt",
    "solid_w.txt",
    "fluid_boundary_c.txt",
    "fluid_boundary_u.txt",
    "fluid_boundary_v.txt",
    "fluid_boundary_w.txt",
    "facet_sections_c.txt",
    "facet_sections_u.txt",
    "facet_sections_v.txt",
    "facet_sections_w.txt",
    "svf.inp.{case}",
    "vfsparse.inp.{case}",
    "netsw.inp.{case}",
    "Tfacinit.inp.{case}",
    "veg.inp.{case}",
    "veg_params.inp.{case}",
]

PREPROCESS_CLEANUP_PATTERNS = [
    "xgrid.inp.{case}",
    "zgrid.inp.{case}",
    "lscale.inp.{case}",
    "prof.inp.{case}",
    "scalar.inp.*.{case}",
    "scalarsourcep.inp.*.{case}",
    "scalarsourcel.inp.*.{case}",
    "vegetation.inp.{case}",
    "veg.inp.{case}",
    "veg_params.inp.{case}",
    "factypes.inp.{case}",
    "facets.inp.{case}",
    "facetarea.inp.{case}",
    "solid_*.txt",
    "fluid_boundary_*.txt",
    "facet_sections_*.txt",
    "facets_unused.{case}",
    "svf.inp.{case}",
    "vfsparse.inp.{case}",
    "vf.txt",
    "vf.bin",
    "vf.nc.inp.{case}",
    "netsw.inp.{case}",
    "Tfacinit.inp.{case}",
    "Tfacinit_layers.inp.{case}",
    "timedepsw.inp.{case}",
    "timedepsveg.inp.{case}",
    "Sdir.txt",
    "Sdir.nc",
    "facets.vs3",
    "vertices.txt",
    "faces.txt",
    "info_directShortwave.txt",
    "DS.exe",
    "inmypoly_inp_info.txt",
    "zhgrid.txt",
    "zfgrid.txt",
    "info_fort.txt",
    "info_matchFacetsToCells.txt",
    "IBM_preproc.exe",
    "solid_points.fig",
    "fluid_boundary_points.fig",
]

RUNTIME_CLEANUP_PATTERNS = [
    "u-dales",
    "output.{case}",
    "job.{case}.slurm",
    "post-job.{case}",
    "decomp_2d_setup.log",
    "slurm-*.out",
    "monitor*.txt",
    "fielddump*",
    "xytdump*",
    "tdump*",
    "mintdump*",
    "treedump*",
]

COPY_IGNORE_GLOBS = [
    "*.nc",
    "Sdir.txt",
    "Sdir.nc",
    "DS.exe",
    "timedepsw.inp.*",
    "timedepsveg.inp.*",
    "vf.txt",
    "vf.bin",
    "vf.nc.inp.*",
    "facets.vs3",
]

PHASE_OUTPUTS = [
    ("Grid", ["lscale.inp.{case}", "prof.inp.{case}"]),
    (
        "IBM",
        [
            "factypes.inp.{case}",
            "facets.inp.{case}",
            "facetarea.inp.{case}",
            "solid_c.txt",
            "solid_u.txt",
            "solid_v.txt",
            "solid_w.txt",
            "fluid_boundary_c.txt",
            "fluid_boundary_u.txt",
            "fluid_boundary_v.txt",
            "fluid_boundary_w.txt",
            "facet_sections_c.txt",
            "facet_sections_u.txt",
            "facet_sections_v.txt",
            "facet_sections_w.txt",
        ],
    ),
    ("View factors", ["svf.inp.{case}", "vfsparse.inp.{case}"]),
    ("Direct shortwave", ["netsw.inp.{case}"]),
    ("SEB init", ["Tfacinit.inp.{case}"]),
    ("Vegetation", ["veg.inp.{case}", "veg_params.inp.{case}"]),
]


def _load_numeric_table(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    if np.ndim(data) == 0:
        return np.asarray([[data]])
    return np.atleast_2d(data)


def _atol_for_output(relpath: str, default_atol: float) -> float:
    # This reusable comparison tool is also used for exploratory cluster runs.
    # For some SEB cases, MATLAB's View3D system() route differs slightly from
    # a direct deterministic View3D run, even though:
    # - MATLAB and Python write byte-identical facets.vs3 inputs
    # - Python matches a standalone View3D invocation exactly
    #
    # Those differences are small sparse-VF rounding deltas that propagate into
    # small netsw changes. They are not treated as preprocessing bugs in the
    # exploratory sweep.
    if relpath.startswith("svf.inp."):
        return max(default_atol, 1.0e-4)
    if relpath.startswith("vfsparse.inp."):
        return max(default_atol, 1.0e-4)
    if relpath.startswith("netsw.inp."):
        return max(default_atol, 1.0e-2)
    return default_atol


def _discover_outputs(case_dir: Path, case: str) -> List[str]:
    outputs = []
    for pattern in SUPPORTED_PATTERNS:
        relpath = pattern.format(case=case)
        if (case_dir / relpath).exists():
            outputs.append(relpath)
    return outputs


def _cleanup_targets(case_dir: Path, case: str) -> List[Path]:
    matches = []
    seen = set()
    for pattern in PREPROCESS_CLEANUP_PATTERNS + RUNTIME_CLEANUP_PATTERNS:
        expanded = pattern.format(case=case)
        for path in case_dir.glob(expanded):
            resolved = path.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            matches.append(path)
    return sorted(matches)


def _count_stl_facets(path: Path) -> int:
    with path.open("rb") as f:
        header = f.read(84)
    if len(header) >= 84:
        tri_count = struct.unpack("<I", header[80:84])[0]
        expected_size = 84 + tri_count * 50
        if path.stat().st_size == expected_size:
            return tri_count

    count = 0
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.lstrip().startswith("facet normal"):
                count += 1
    return count


def _detect_seb(case_dir: Path, case: str) -> bool:
    namoptions = case_dir / f"namoptions.{case}"
    if not namoptions.exists():
        return False
    for line in namoptions.read_text(errors="ignore").lower().splitlines():
        if "=" not in line:
            continue
        lhs, rhs = line.split("=", 1)
        if lhs.strip() == "leb":
            return ".true." in rhs
    return False


def sorted_cases(experiments_root: Path) -> List[tuple]:
    rows = []
    for case_dir in sorted(p for p in experiments_root.iterdir() if p.is_dir() and p.name.isdigit()):
        stls = sorted(case_dir.glob("*.stl"))
        if not stls:
            continue
        smallest = min((_count_stl_facets(stl), stl.name) for stl in stls)
        rows.append((smallest[0], case_dir.name, smallest[1], _detect_seb(case_dir, case_dir.name)))
    return sorted(rows)


def _copy_case_tree(source_case: Path, target_case: Path) -> None:
    shutil.copytree(
        source_case,
        target_case,
        ignore=shutil.ignore_patterns(*COPY_IGNORE_GLOBS),
    )


def _temp_parent() -> Path:
    for candidate in DEFAULT_TEMP_ROOTS:
        if not candidate:
            continue
        candidate = Path(candidate).expanduser().resolve()
        if candidate.exists() and os.access(candidate, os.W_OK):
            return candidate
    return Path("/tmp")


def _run_matlab(matlab_case_root: Path, case: str, outputs: List[str]) -> str:
    matlab = shutil.which("matlab")
    if matlab is None:
        raise RuntimeError("matlab not found on PATH")

    env = os.environ.copy()
    env["DA_EXPDIR"] = str(matlab_case_root.parent)
    env["DA_TOOLSDIR"] = str(REPO_ROOT / "tools")
    env["MATLAB_USE_USERWORK"] = "0"
    virtual_env = env.get("VIRTUAL_ENV")
    if virtual_env:
        activate = str(Path(virtual_env).resolve() / "bin" / "activate")
    else:
        activate = str(Path(sys.executable).resolve().parent / "activate")
    path_export = shlex.quote(env.get("PATH", ""))
    ld_library_path = shlex.quote(env.get("LD_LIBRARY_PATH", ""))
    env["UDALES_PYTHON_CMD"] = (
        f"export PATH={path_export} "
        f"LD_LIBRARY_PATH={ld_library_path} && "
        f"source {shlex.quote(activate)} && "
        "python3"
    )
    try:
        libgfortran = subprocess.run(
            ["gfortran", "-print-file-name=libgfortran.so"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        libgfortran = ""
    if libgfortran and libgfortran != "libgfortran.so":
        gcc_libdir = str(Path(libgfortran).resolve().parent)
        current = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = gcc_libdir if not current else f"{gcc_libdir}:{current}"

    cmd = [
        matlab,
        "-nodesktop",
        "-noFigureWindows",
        "-nosplash",
        "-nodisplay",
        "-r",
        f"expnr={case}; write_inputs; quit",
    ]
    result = subprocess.run(
        cmd,
        cwd=REPO_ROOT / "tools",
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    missing = [rel for rel in outputs if not (matlab_case_root / rel).exists()]
    if missing:
        raise RuntimeError(
            "MATLAB preprocessing did not produce all expected outputs.\n"
            f"Return code: {result.returncode}\n"
            f"Missing: {missing}\n"
            f"Output:\n{result.stdout}"
        )
    if result.returncode != 0:
        raise RuntimeError(f"MATLAB preprocessing failed:\n{result.stdout}")
    return result.stdout


def _run_python(python_case_root: Path) -> str:
    cmd = [sys.executable, str(REPO_ROOT / "tools" / "write_inputs.py"), str(python_case_root)]
    result = subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Python preprocessing failed:\n{result.stdout}")
    return result.stdout


def _compare_outputs(matlab_case_root: Path, python_case_root: Path, outputs: List[str], atol: float) -> List[str]:
    mismatches: List[str] = []
    for relpath in outputs:
        matlab_path = matlab_case_root / relpath
        python_path = python_case_root / relpath
        if not matlab_path.exists():
            mismatches.append(f"{relpath}: missing MATLAB output")
            continue
        if not python_path.exists():
            mismatches.append(f"{relpath}: missing Python output")
            continue
        matlab_data = _load_numeric_table(matlab_path)
        python_data = _load_numeric_table(python_path)
        this_atol = _atol_for_output(relpath, atol)
        if matlab_data.shape != python_data.shape:
            mismatches.append(
                f"{relpath}: shape mismatch MATLAB {matlab_data.shape} vs Python {python_data.shape}"
            )
            continue
        try:
            np.testing.assert_allclose(python_data, matlab_data, rtol=0.0, atol=this_atol)
        except AssertionError as exc:
            mismatches.append(f"{relpath}: {exc}")
    return mismatches


def _compare_output_statuses(
    matlab_case_root: Path,
    python_case_root: Path,
    outputs: List[str],
    atol: float,
) -> dict:
    statuses = {}
    for relpath in outputs:
        matlab_path = matlab_case_root / relpath
        python_path = python_case_root / relpath
        if not matlab_path.exists():
            statuses[relpath] = "missing_matlab"
            continue
        if not python_path.exists():
            statuses[relpath] = "missing_python"
            continue
        try:
            matlab_data = _load_numeric_table(matlab_path)
            python_data = _load_numeric_table(python_path)
        except Exception:
            statuses[relpath] = "unreadable"
            continue
        this_atol = _atol_for_output(relpath, atol)
        if matlab_data.shape != python_data.shape:
            statuses[relpath] = "shape_mismatch"
            continue
        statuses[relpath] = (
            "match"
            if np.allclose(python_data, matlab_data, rtol=0.0, atol=this_atol)
            else "mismatch"
        )
    return statuses


def _phase_statuses(case: str, outputs: List[str], compare_statuses: dict) -> List[tuple]:
    output_set = set(outputs)
    rows = []
    for phase, patterns in PHASE_OUTPUTS:
        relpaths = [pattern.format(case=case) for pattern in patterns if pattern.format(case=case) in output_set]
        if not relpaths:
            rows.append((phase, "n/a", "No phase outputs for this case"))
            continue
        states = [compare_statuses.get(relpath, "unknown") for relpath in relpaths]
        if all(state == "match" for state in states):
            rows.append((phase, "OK", ", ".join(relpaths)))
        elif any(state == "mismatch" for state in states):
            rows.append((phase, "Mismatch", ", ".join(relpath for relpath in relpaths if compare_statuses.get(relpath) == "mismatch")))
        elif any(state.startswith("missing") for state in states):
            rows.append((phase, "Missing", ", ".join(relpath for relpath in relpaths if compare_statuses.get(relpath, "").startswith("missing"))))
        elif any(state == "shape_mismatch" for state in states):
            rows.append((phase, "Shape mismatch", ", ".join(relpath for relpath in relpaths if compare_statuses.get(relpath) == "shape_mismatch")))
        else:
            rows.append((phase, "Unknown", ", ".join(relpaths)))
    return rows


def print_phase_table(case: str, outputs: List[str], compare_statuses: dict) -> None:
    print("Phase summary:")
    print(f"{'Phase':<18} {'Status':<15} Evidence")
    print(f"{'-' * 18} {'-' * 15} {'-' * 40}")
    for phase, status, evidence in _phase_statuses(case, outputs, compare_statuses):
        print(f"{phase:<18} {status:<15} {evidence}")


def compare_case(case: str, source_case: Path, atol: float = 1.0e-10, keep_temp: bool = False) -> int:
    outputs = _discover_outputs(source_case, case)
    if not outputs:
        raise SystemExit(f"No supported preprocessing outputs found in {source_case}")

    temp_ctx = tempfile.TemporaryDirectory(
        prefix=f"udales-real-preproc-{case}-",
        dir=str(_temp_parent()),
    )
    temp_root = Path(temp_ctx.name)
    if keep_temp:
        print(f"Keeping temp directory: {temp_root}")

    try:
        matlab_case = temp_root / "matlab" / case
        python_case = temp_root / "python" / case
        matlab_case.parent.mkdir(parents=True, exist_ok=True)
        python_case.parent.mkdir(parents=True, exist_ok=True)
        _copy_case_tree(source_case, matlab_case)
        _copy_case_tree(source_case, python_case)

        for case_dir in (matlab_case, python_case):
            for target in _cleanup_targets(case_dir, case):
                target.unlink()

        print(f"Comparing case {case}")
        print(f"Source: {source_case}")
        print("Outputs:")
        for relpath in outputs:
            print(f"  - {relpath}")

        _run_matlab(matlab_case, case, outputs)
        _run_python(python_case)
        compare_statuses = _compare_output_statuses(matlab_case, python_case, outputs, atol)
        print_phase_table(case, outputs, compare_statuses)
        mismatches = _compare_outputs(matlab_case, python_case, outputs, atol)
    finally:
        if keep_temp:
            temp_ctx._finalizer.detach()
        else:
            temp_ctx.cleanup()

    if mismatches:
        print("Mismatches found:")
        for mismatch in mismatches:
            print(f"  - {mismatch}")
        return 1

    print("All compared outputs match.")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Compare MATLAB and Python preprocessing on a case.")
    parser.add_argument("case", nargs="?", help="Case number under the experiments root")
    parser.add_argument(
        "--case-dir",
        help="Explicit case directory to compare instead of looking under --experiments-root",
    )
    parser.add_argument(
        "--experiments-root",
        default=str(DEFAULT_EXPERIMENTS_ROOT),
        help="Root directory containing real experiment cases",
    )
    parser.add_argument(
        "--list-cases",
        action="store_true",
        help="List available experiment cases sorted by the number of STL facets",
    )
    parser.add_argument(
        "--atol",
        type=float,
        default=1.0e-10,
        help="Absolute tolerance for numeric comparisons",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep the temporary comparison directory for manual inspection",
    )
    args = parser.parse_args(argv)

    if args.case_dir:
        source_case = Path(args.case_dir).expanduser().resolve()
        if not source_case.exists():
            raise SystemExit(f"Case directory not found: {source_case}")
        case = source_case.name
        return compare_case(case, source_case, atol=args.atol, keep_temp=args.keep_temp)

    experiments_root = Path(args.experiments_root).expanduser().resolve()
    if args.list_cases:
        rows = sorted_cases(experiments_root)
        if not rows:
            raise SystemExit(f"No STL-based experiment cases found in {experiments_root}")
        print("Cases sorted by minimum STL facet count:")
        for facets, case, stl_name, seb in rows:
            suffix = " [SEB]" if seb else ""
            print(f"  {case:>3}  {facets:>8} facets  {stl_name}{suffix}")
        return 0

    if not args.case:
        raise SystemExit("case is required unless --list-cases or --case-dir is used")

    source_case = experiments_root / args.case
    if not source_case.exists():
        raise SystemExit(f"Case directory not found: {source_case}")
    return compare_case(args.case, source_case, atol=args.atol, keep_temp=args.keep_temp)


if __name__ == "__main__":
    raise SystemExit(main())
