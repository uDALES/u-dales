#!/usr/bin/env python3
"""Cross-solver check: ipoiss=0 on the CPU against ipoiss=3 on CPU and GPU.

The two Poisson solvers discretise the same operator but use different
libraries (FFTW for ipoiss=0, 2DECOMP&FFT with FFTW or cuFFT for ipoiss=3), so
the runs agree to round-off rather than bit for bit. One short simulation of
the dry case 103 is run with ipoiss=0 on one CPU rank as the reference; then
ipoiss=3 is run on the CPU build with 1, 2 (x), 2 (y) and 4 ranks and, when a
GPU executable is given, on the GPU build with the same decompositions. Every
requested output (fielddump, tdump, xytdump and the double-precision restart
files) of every candidate must match the reference within the tolerance.

Why it exists: the ipoiss=3 solver once indexed its local spectral-pencil
arrays with global bounds, which is invisible on one rank. Comparing every
multi-rank ipoiss=3 run against an independent solver is the check that the
decomposed solve gives the same physics, not merely a finite answer.

Environment (same names as the GPU parity harness): UDALES_CPU_MPIEXEC,
UDALES_GPU_MPIEXEC, MPIEXEC, UDALES_CPU_MPI_ARGS, UDALES_GPU_MPI_ARGS.
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
GPU_TEST_DIR = REPO_ROOT / "tests" / "integration" / "gpu"
sys.path.insert(0, str(GPU_TEST_DIR))

from compare_outputs import compare_output_directories  # noqa: E402
from run_gpu_tests import _run_solver  # noqa: E402

CASE_ID = 103
BASE_CASE = REPO_ROOT / "tests" / "regression" / "david_tests" / "cases" / "103"
BASE_NAMELIST = GPU_TEST_DIR / "namelists" / "namoptions.103.dry"
# The NetCDF dumps are single precision, so they hide anything below about
# 1e-7 relative; the restart files (initd) hold the double-precision state and
# are the sensitive part of the comparison. TRESTART makes every run write one
# at t = 1 and at t = 2.
REQUIRED_OUTPUTS = ("fielddump", "tdump", "xytdump", "initd")
TRESTART = 1.0
IRANDOM = 43
RANDU = 0.1
DECOMPOSITIONS = (("1x1", 1, 1), ("2x1", 2, 1), ("1x2", 1, 2), ("2x2", 2, 2))


def _namelist_text(ipoiss: int, nprocx: int, nprocy: int) -> str:
    text = BASE_NAMELIST.read_text(encoding="utf-8")
    for key, value in (("ipoiss", ipoiss), ("nprocx", nprocx), ("nprocy", nprocy), ("trestart", TRESTART)):
        text, count = re.subn(rf"^(\s*{key}\s*=\s*)\S+", rf"\g<1>{value}", text, flags=re.M)
        if count != 1:
            raise SystemExit(f"{BASE_NAMELIST}: expected exactly one '{key}' line, found {count}")
    # The dry case starts from a uniform flow: without perturbations the
    # divergence is zero, the pressure stays zero, and the solvers are never
    # exercised. randomize_field seeds by global cell index, so the perturbed
    # start is the same on every decomposition.
    text, count = re.subn(
        r"^(\s*lrandomize\s*=\s*)\S+",
        rf"\g<1>.true.\nirandom      = {IRANDOM}\nrandu        = {RANDU}",
        text,
        flags=re.M,
    )
    if count != 1:
        raise SystemExit(f"{BASE_NAMELIST}: expected exactly one 'lrandomize' line, found {count}")
    return text


def _stage(run_dir: Path, ipoiss: int, nprocx: int, nprocy: int) -> Path:
    shutil.copytree(BASE_CASE, run_dir)
    namelist = run_dir / f"namoptions.{CASE_ID:03d}"
    namelist.write_text(_namelist_text(ipoiss, nprocx, nprocy), encoding="utf-8")
    return namelist


def _run(
    label: str,
    executable: Path,
    ipoiss: int,
    nprocx: int,
    nprocy: int,
    gpu: bool,
    work_root: Path,
    timeout: float,
) -> Dict[str, Any]:
    run_dir = work_root / label
    namelist = _stage(run_dir, ipoiss, nprocx, nprocy)
    print(f"==> {label}: ipoiss={ipoiss} {nprocx}x{nprocy} on {executable}")
    result = _run_solver(executable, run_dir, namelist, nprocx * nprocy, timeout, gpu=gpu)
    result["dir"] = run_dir
    if result["returncode"] != 0:
        tail = "".join((run_dir / "run.log").read_text(errors="replace").splitlines(True)[-40:])
        print(f"    run exited with {result['returncode']}\n{tail}", file=sys.stderr)
    return result


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cpu-executable", type=Path, required=True)
    parser.add_argument("--gpu-executable", type=Path, default=None, help="add ipoiss=3 GPU candidates")
    parser.add_argument("--max-ranks", type=int, default=4, help="largest CPU decomposition to run")
    parser.add_argument("--max-gpu-ranks", type=int, default=None, help="largest GPU decomposition (default: --max-ranks)")
    parser.add_argument("--atol", type=float, default=1.0e-9)
    parser.add_argument("--rtol", type=float, default=1.0e-9)
    parser.add_argument("--timeout", type=float, default=300.0, help="seconds per run")
    parser.add_argument("--work-root", type=Path, default=None, help="keep run directories here")
    return parser.parse_args()


def _compare(reference_dir: Path, candidate_dir: Path, outputs, tolerance) -> tuple[bool, str]:
    comparison = compare_output_directories(reference_dir, candidate_dir, CASE_ID, outputs, tolerance)
    worst = max(
        (r for r in comparison.variable_results if "max_abs_error" in r),
        key=lambda r: r["max_abs_error"],
        default=None,
    )
    detail = ""
    if worst is not None:
        detail = (
            f"{comparison.variables_compared} variables in {comparison.files_compared} file(s); "
            f"largest difference {worst['max_abs_error']:.3e} in {worst['file']}:{worst['variable']} "
            f"(reference value {worst.get('reference', float('nan')):.6e})"
        )
    for failure in comparison.failures:
        print(f"    {failure}", file=sys.stderr)
    return comparison.passed, detail


def main() -> int:
    args = _parse_args()
    cpu = args.cpu_executable.resolve()
    gpu: Optional[Path] = args.gpu_executable.resolve() if args.gpu_executable else None
    for executable in (cpu, gpu):
        if executable is not None and not executable.is_file():
            print(f"ERROR: executable not found: {executable}", file=sys.stderr)
            return 2
    max_gpu_ranks = args.max_gpu_ranks if args.max_gpu_ranks is not None else args.max_ranks
    if gpu is None:
        max_gpu_ranks = 0
    tolerance = {
        "default": {"atol": args.atol, "rtol": args.rtol},
        "variables": {"time": {"atol": 1.0e-12, "rtol": 1.0e-12}},
    }

    if args.work_root:
        work_root = args.work_root.resolve()
        work_root.mkdir(parents=True, exist_ok=True)
        cleanup = None
    else:
        cleanup = tempfile.TemporaryDirectory(prefix="poisson-equivalence-")
        work_root = Path(cleanup.name)

    print(f"tolerance: atol={args.atol:g} rtol={args.rtol:g}   outputs: {', '.join(REQUIRED_OUTPUTS)}")
    print("reference for every decomposition: ipoiss=0 on the CPU build with the same nprocx x nprocy;")
    print("xytdump (one file for the whole domain) is also compared with the 1x1 reference.")

    results: list[tuple[str, bool, str]] = []
    reference_1x1: Optional[Dict[str, Any]] = None
    for name, nprocx, nprocy in DECOMPOSITIONS:
        nprocs = nprocx * nprocy
        if nprocs > max(args.max_ranks, max_gpu_ranks):
            print(f"==> {name}: skipped (needs {nprocs} ranks, limits cpu={args.max_ranks} gpu={max_gpu_ranks})")
            continue
        reference = _run(f"cpu-ipoiss0-{name}", cpu, 0, nprocx, nprocy, False, work_root, args.timeout)
        if reference["returncode"] != 0:
            results.append((f"cpu-ipoiss0-{name}", False, f"reference run exited with {reference['returncode']}"))
            continue
        if reference_1x1 is None:
            reference_1x1 = reference
        elif nprocs > 1:
            # The ipoiss=0 reference on this decomposition must itself agree
            # with the 1x1 reference in the domain-averaged statistics.
            passed, detail = _compare(reference_1x1["dir"], reference["dir"], ("xytdump",), tolerance)
            print(f"    xytdump vs 1x1 reference: {detail}\n    result: {'PASS' if passed else 'FAIL'}")
            results.append((f"cpu-ipoiss0-{name} xytdump vs 1x1", passed, ""))

        candidates = []
        if nprocs <= args.max_ranks:
            candidates.append(("cpu", cpu, False))
        if nprocs <= max_gpu_ranks:
            candidates.append(("gpu", gpu, True))
        for kind, executable, is_gpu in candidates:
            label = f"{kind}-ipoiss3-{name}"
            run = _run(label, executable, 3, nprocx, nprocy, is_gpu, work_root, args.timeout)
            if run["returncode"] != 0:
                results.append((label, False, f"run exited with {run['returncode']}"))
                print("    result: FAIL")
                continue
            passed, detail = _compare(reference["dir"], run["dir"], REQUIRED_OUTPUTS, tolerance)
            print(f"    all outputs vs ipoiss=0 {name}: {detail}")
            results.append((label, passed, ""))
            if nprocs > 1:
                passed_xy, detail_xy = _compare(reference_1x1["dir"], run["dir"], ("xytdump",), tolerance)
                print(f"    xytdump vs 1x1 reference: {detail_xy}")
                results.append((f"{label} xytdump vs 1x1", passed_xy, ""))
                passed = passed and passed_xy
            print(f"    result: {'PASS' if passed else 'FAIL'}")

    print("\nSummary (reference: ipoiss=0 on the CPU build)")
    for label, passed, note in results:
        print(f"- {label}: {'PASS' if passed else 'FAIL'}{'  ' + note if note else ''}")
    overall = bool(results) and all(passed for _, passed, _ in results)
    print(f"overall: {'PASS' if overall else 'FAIL'}")
    if cleanup is not None and overall:
        cleanup.cleanup()
    elif cleanup is not None:
        print(f"run directories kept in {work_root}")
        cleanup._finalizer.detach()  # type: ignore[attr-defined]
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())
