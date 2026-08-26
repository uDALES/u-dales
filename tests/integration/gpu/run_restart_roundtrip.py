#!/usr/bin/env python3
"""Restart round trip and warm-start parity, on both builds, from nothing.

Four solver runs of one small case:

    cold CPU  (0 -> T2, restart written at T1)     cold GPU  (same)
         |  its own restart                              |  its own restart
    warm CPU  (T1 -> T2)                            warm GPU  (T1 -> T2)

and four comparisons:

  - cold CPU vs cold GPU: the ordinary parity check, restart files included.
  - warm CPU vs warm GPU: parity of two runs that each restarted from their
    own build's file. This is what exercises the GPU's restart writer as part
    of a physical continuation: updateHostForRestart -> modsave -> readrestartfiles.
  - cold vs warm, per build: the round trip. At T2 the cold run's last field
    dump and the warm run's last field dump must agree, and so must the
    restart both runs write there. A restart that loses any piece of state
    changes the trajectory, and a lost piece is lost identically on both
    builds, so warm-versus-warm parity alone cannot see it.

    On the CPU the demand is bit-exactness: the case has no atomics, so the
    same binary is reproducible and any difference at all is a lost or
    mis-restored piece of state. On the GPU it is roundoff-exactness, 1e-11
    absolute, five orders tighter than the parity tolerance. The reason is
    not the restart: a warm start forms the slab averages (u0av, thvh and the
    rest) with the host reduction in readinitfiles, while the run it continues
    formed them with the device reduction in the loop, and the two sum in a
    different order. The columns differ in the last bit, and that carries.

Nothing is committed to the repository for this: the restart files are
produced by the cold runs inside a temporary directory and removed with it.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import netCDF4
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from compare_outputs import (  # noqa: E402
    Tolerance,
    _compare_numeric_variable,
    _compare_restart_file,
    _masked_values,
    compare_output_directories,
)
from run_gpu_tests import (  # noqa: E402
    REPO_ROOT,
    _check_executables,
    _log_tail,
    _read_matrix,
    _require_cuda_selftest,
    _run_solver,
    _stage_case,
    DEFAULT_MATRIX,
)

CASE: Dict[str, Any] = {
    "name": "restart-roundtrip",
    "case_id": 103,
    "base_case": "tests/regression/david_tests/cases/103",
    "namelist": "tests/integration/gpu/namelists/namoptions.103.restart-roundtrip",
    "overrides": [
        "tests/integration/gpu/fixtures/prof.inp.103",
        "tests/integration/gpu/fixtures/scalar.inp.103",
    ],
    "nprocs": 1,
    "timeout_seconds": 300,
}
PARITY_OUTPUTS_COLD = ["fielddump", "tdump", "xytdump", "initd", "inits"]
PARITY_OUTPUTS_WARM = ["fielddump", "tdump", "xytdump"]
EXACT = {"default": {"atol": 0.0, "rtol": 0.0}}
# Observed on the A4000: 1.7e-13 on thl0 (three ulps at 288 K), below 4e-15
# on everything else. 1e-11 leaves that sixty times over and is still 1e5
# tighter than the CPU/GPU parity tolerance.
ROUNDOFF = {"default": {"atol": 1.0e-11, "rtol": 0.0}}
ROUNDTRIP_TOLERANCE = {"cpu": EXACT, "gpu": ROUNDOFF}


def _namelist_value(text: str, key: str) -> str:
    match = re.search(rf"(?m)^\s*{key}\s*=\s*([^!\n]+)", text)
    if not match:
        raise ValueError(f"namelist has no {key}")
    return match.group(1).strip()


def _set_namelist(text: str, key: str, value: str) -> str:
    """Replace the value of an existing namelist entry, keeping its line."""
    pattern = rf"(?m)^(\s*{key}\s*=\s*)[^!\n]*"
    if not re.search(pattern, text):
        raise ValueError(f"namelist has no {key} to substitute")
    return re.sub(pattern, lambda m: f"{m.group(1)}{value}", text, count=1)


def warm_namelist(cold_text: str, startfile: str, t1: float) -> str:
    """The warm-start namelist derived from the cold one.

    runtime is relative to the restart time on a warm start (timeleft is set
    to runtime and btime to the restored timee), so the warm run gets
    T2 - T1 to reach the same end time as the cold run.
    """
    t2 = float(_namelist_value(cold_text, "runtime"))
    text = _set_namelist(cold_text, "lwarmstart", ".true.")
    text = _set_namelist(text, "startfile", f"'{startfile}'")
    text = _set_namelist(text, "runtime", f"{t2 - t1:g}")
    return text


def first_restart(run_dir: Path, case_id: int) -> Tuple[Path, List[Path]]:
    """The earliest initd file a run wrote, and every restart file with it."""
    initd = sorted(run_dir.glob(f"initd*_000_000.{case_id:03d}"))
    if not initd:
        raise FileNotFoundError(f"no initd file in {run_dir}")
    earliest = initd[0]
    step = earliest.name[5:13]
    companions = sorted(run_dir.glob(f"init?{step}_*.{case_id:03d}"))
    return earliest, companions


def compare_last_records(reference: Path, candidate: Path,
                         tolerance: Tolerance = Tolerance(0.0, 0.0)) -> Tuple[List[str], List[Dict[str, Any]]]:
    """Compare the final time record of every variable in two fielddump files."""
    failures: List[str] = []
    results: List[Dict[str, Any]] = []
    label_file = f"{reference.parent.name}/{reference.name} vs {candidate.parent.name}"
    with netCDF4.Dataset(reference) as ref_ds, netCDF4.Dataset(candidate) as cand_ds:
        if set(ref_ds.variables) != set(cand_ds.variables):
            failures.append(f"{label_file}: variable sets differ")
        for name in sorted(set(ref_ds.variables) & set(cand_ds.variables)):
            ref_var, cand_var = ref_ds.variables[name], cand_ds.variables[name]
            if "time" in ref_var.dimensions:
                ref_values, ref_mask = _masked_values(ref_var[-1, ...])
                cand_values, cand_mask = _masked_values(cand_var[-1, ...])
            else:
                ref_values, ref_mask = _masked_values(ref_var[:])
                cand_values, cand_mask = _masked_values(cand_var[:])
            if ref_values.shape != cand_values.shape:
                failures.append(f"{label_file}:{name}: shape {ref_values.shape} vs {cand_values.shape}")
                continue
            if not np.array_equal(ref_mask, cand_mask):
                failures.append(f"{label_file}:{name}: mask differs")
                continue
            if not np.issubdtype(ref_values.dtype, np.number):
                continue
            passed, detail = _compare_numeric_variable(ref_values, cand_values, ~ref_mask, tolerance)
            results.append({"file": label_file, "variable": name, "passed": passed, **detail})
            if not passed:
                failures.append(f"{label_file}:{name}: {json.dumps(detail, sort_keys=True)}")
    return failures, results


def _fielddump(run_dir: Path, case_id: int) -> Path:
    files = sorted(run_dir.glob(f"fielddump*.{case_id:03d}.nc"))
    if len(files) != 1:
        raise FileNotFoundError(f"expected one fielddump in {run_dir}, found {[f.name for f in files]}")
    return files[0]


def _last_restart(run_dir: Path, case_id: int) -> List[Path]:
    initd = sorted(run_dir.glob(f"initd*_000_000.{case_id:03d}"))
    if not initd:
        raise FileNotFoundError(f"no initd file in {run_dir}")
    step = initd[-1].name[5:13]
    return sorted(run_dir.glob(f"init?{step}_*.{case_id:03d}"))


def _run(label: str, executable: Path, run_dir: Path, namelist: Path, gpu: bool,
         selftest: bool, failures: List[str]) -> Dict[str, Any]:
    result = _run_solver(executable, run_dir, namelist, int(CASE["nprocs"]),
                         float(CASE["timeout_seconds"]), gpu, run_cuda_selftest=gpu and selftest)
    if result["returncode"] != 0:
        failures.append(f"{label} exited with {result['returncode']}\n{_log_tail(run_dir / 'run.log')}")
    elif gpu and selftest:
        problem = _require_cuda_selftest(run_dir / "run.log", int(CASE["nprocs"]))
        if problem:
            failures.append(f"{label}: {problem}")
    print(f"  {label}: {'ok' if result['returncode'] == 0 else 'FAILED'} ({result['elapsed_seconds']:.1f} s)")
    return result


def run_roundtrip(cpu_executable: Path, gpu_executable: Path, run_root: Path,
                  tolerances: Dict[str, Any], selftest: bool) -> Dict[str, Any]:
    case_id = int(CASE["case_id"])
    failures: List[str] = []
    checks: Dict[str, Any] = {}

    print("==> restart round trip")
    dirs = {name: run_root / name for name in ("cold_cpu", "cold_gpu", "warm_cpu", "warm_gpu")}
    cold_namelist = _stage_case(CASE, dirs["cold_cpu"])
    _stage_case(CASE, dirs["cold_gpu"])
    cold_text = cold_namelist.read_text()
    t1 = float(_namelist_value(cold_text, "trestart"))

    _run("cold CPU", cpu_executable, dirs["cold_cpu"], cold_namelist, False, selftest, failures)
    _run("cold GPU", gpu_executable, dirs["cold_gpu"], cold_namelist, True, selftest, failures)
    if failures:
        return {"passed": False, "failures": failures, "checks": checks}

    # Each warm run restarts from the file its own build wrote.
    for side in ("cpu", "gpu"):
        cold_dir, warm_dir = dirs[f"cold_{side}"], dirs[f"warm_{side}"]
        _stage_case(CASE, warm_dir)
        try:
            earliest, companions = first_restart(cold_dir, case_id)
        except FileNotFoundError as error:
            failures.append(str(error))
            continue
        for path in companions:
            shutil.copy2(path, warm_dir / path.name)
        (warm_dir / cold_namelist.name).write_text(warm_namelist(cold_text, earliest.name, t1))
        checks[f"warm_{side}_startfile"] = earliest.name
    if failures:
        return {"passed": False, "failures": failures, "checks": checks}

    _run("warm CPU", cpu_executable, dirs["warm_cpu"], dirs["warm_cpu"] / cold_namelist.name, False, selftest, failures)
    _run("warm GPU", gpu_executable, dirs["warm_gpu"], dirs["warm_gpu"] / cold_namelist.name, True, selftest, failures)
    if failures:
        return {"passed": False, "failures": failures, "checks": checks}

    # 1. cold parity, restart files included
    cold = compare_output_directories(dirs["cold_cpu"], dirs["cold_gpu"], case_id, PARITY_OUTPUTS_COLD, tolerances)
    checks["parity_cold"] = cold.as_dict()
    failures.extend(f"cold parity: {f}" for f in cold.failures)
    print(f"  cold CPU vs cold GPU: {'ok' if cold.passed else 'FAILED'} ({cold.variables_compared} variables)")

    # 2. warm parity
    warm = compare_output_directories(dirs["warm_cpu"], dirs["warm_gpu"], case_id, PARITY_OUTPUTS_WARM, tolerances)
    checks["parity_warm"] = warm.as_dict()
    failures.extend(f"warm parity: {f}" for f in warm.failures)
    print(f"  warm CPU vs warm GPU: {'ok' if warm.passed else 'FAILED'} ({warm.variables_compared} variables)")

    # 3. round trip per build: last field record, and the restart written at T2
    for side in ("cpu", "gpu"):
        cold_dir, warm_dir = dirs[f"cold_{side}"], dirs[f"warm_{side}"]
        side_failures: List[str] = []
        tol_config = ROUNDTRIP_TOLERANCE[side]
        tol = Tolerance(float(tol_config["default"]["atol"]), float(tol_config["default"]["rtol"]))
        try:
            f, r = compare_last_records(_fielddump(cold_dir, case_id), _fielddump(warm_dir, case_id), tol)
            side_failures.extend(f)
            n_fields = len(r)
            cold_restart = {p.name: p for p in _last_restart(cold_dir, case_id)}
            warm_restart = {p.name: p for p in _last_restart(warm_dir, case_id)}
            if set(cold_restart) != set(warm_restart):
                side_failures.append(
                    f"restart files at T2 differ: cold={sorted(cold_restart)}, warm={sorted(warm_restart)}"
                )
            n_records = 0
            for name in sorted(set(cold_restart) & set(warm_restart)):
                f, r = _compare_restart_file(cold_restart[name], warm_restart[name], tol_config)
                side_failures.extend(f)
                n_records += len(r)
        except FileNotFoundError as error:
            side_failures.append(str(error))
            n_fields = n_records = 0
        checks[f"roundtrip_{side}"] = {"passed": not side_failures, "failures": side_failures}
        failures.extend(f"round trip {side}: {x}" for x in side_failures)
        how = "bit-exact" if tol.atol == 0.0 else f"within {tol.atol:g}"
        print(f"  cold {side.upper()} vs warm {side.upper()} at T2: "
              f"{'ok' if not side_failures else 'FAILED'} "
              f"({how} over {n_fields} field variables and {n_records} restart records)")

    return {"passed": not failures, "failures": failures, "checks": checks}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cpu-executable", type=Path,
                        default=Path(os.environ.get("UDALES_CPU_BUILD", REPO_ROOT / "build/cpu/debug/u-dales")))
    parser.add_argument("--gpu-executable", type=Path,
                        default=Path(os.environ.get("UDALES_GPU_BUILD", REPO_ROOT / "build/gpu/debug/u-dales")))
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX,
                        help="only the tolerances block is read from it")
    parser.add_argument("--work-root", type=Path, default=None)
    parser.add_argument("--keep-work-dir", action="store_true")
    parser.add_argument("--require-debug-selftest", action="store_true")
    parser.add_argument("--summary", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    cpu_executable = args.cpu_executable.resolve()
    gpu_executable = args.gpu_executable.resolve()
    _check_executables(cpu_executable, gpu_executable)
    tolerances = _read_matrix(args.matrix.resolve()).get("tolerances", {})

    work_parent = args.work_root.resolve() if args.work_root else None
    if work_parent:
        work_parent.mkdir(parents=True, exist_ok=True)
    temporary = tempfile.TemporaryDirectory(prefix="udales-restart-roundtrip-", dir=work_parent)
    run_root = Path(temporary.name)
    print(f"run root: {run_root}")

    summary = run_roundtrip(cpu_executable, gpu_executable, run_root, tolerances, args.require_debug_selftest)
    summary.update({"cpu_executable": str(cpu_executable), "gpu_executable": str(gpu_executable)})
    if args.summary:
        args.summary.parent.mkdir(parents=True, exist_ok=True)
        args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True))

    if args.keep_work_dir:
        temporary._finalizer.detach()  # type: ignore[attr-defined]
        print(f"kept work dir: {run_root}")
    else:
        temporary.cleanup()

    print("\nSummary")
    for failure in summary["failures"]:
        print(f"- {failure}")
    print(f"overall: {'PASS' if summary['passed'] else 'FAIL'}")
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
