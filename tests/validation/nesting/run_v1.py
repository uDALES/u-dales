#!/usr/bin/env python3
"""End-to-end driver for the V1 Big Brother nesting validation.

Runs the whole pipeline at whichever preset is asked for:

    parent case -> spin-up -> production dumps
        -> child case (geometry + nesting input, cut from those dumps)
        -> nested child run
        -> analysis (JSON + CSV + plots)

The tiny preset exists so that this exact code path can be exercised in a few
minutes on a login node before the production run is submitted.  There is no
separate small-case code path: the only difference is the :class:`config.Preset`
that is passed in.

Usage
-----
    python run_v1.py <rundir> --preset tiny
    python run_v1.py <rundir> --preset production

Stages can be skipped so a failed run can be resumed:

    python run_v1.py <rundir> --preset production --start-at child
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from pathlib import Path

import make_child_case
import make_parent_case
import analyse
from caselib import run_solver
from config import get_preset

STAGES = ("parent-case", "spinup", "production", "child-case", "child", "analysis")


def _restart_file(casedir: Path, expnr: str) -> str:
    """Name of the restart file the spin-up wrote (rank 0,0).

    ``modsave`` names it ``initd<ntrun:08d>_<x>_<y>.<expnr>`` where ``ntrun`` is
    the timestep counter, which is not known in advance under adaptive time
    stepping -- so it is discovered rather than predicted.
    """
    pattern = re.compile(rf"^initd(\d{{8}})_000_000\.{re.escape(expnr)}$")
    matches = [(int(m.group(1)), p.name)
               for p in casedir.iterdir() if (m := pattern.match(p.name))]
    if not matches:
        raise FileNotFoundError(
            f"no initd????????_000_000.{expnr} restart file in {casedir}; "
            "did the spin-up phase finish?"
        )
    return max(matches)[1]


def _set_startfile(namoptions: Path, name: str) -> None:
    text = namoptions.read_text(encoding="ascii")
    new, n = re.subn(r"(?m)^(\s*startfile\s*=\s*).*$", lambda m: m.group(1) + f"'{name}'",
                     text)
    if n != 1:
        raise RuntimeError(f"could not set startfile in {namoptions}")
    namoptions.write_text(new, encoding="ascii")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("rundir", type=Path)
    parser.add_argument("--preset", default="tiny")
    parser.add_argument("--start-at", default="parent-case", choices=STAGES)
    parser.add_argument("--stop-after", default="analysis", choices=STAGES)
    parser.add_argument("--ibm-backend", default="auto")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--yes", action="store_true",
                        help="confirm a large run outside a batch job")
    args = parser.parse_args()

    preset = get_preset(args.preset)
    # Guard against the interactive footgun: the production preset launches a
    # 64-rank multi-hour run and writes tens of GB.  Inside a batch job (or with
    # an explicit --yes / UDALES_V1_CONFIRM) that is exactly what is wanted.
    nprocs = preset.nprocx * preset.nprocy
    if nprocs > 8 and not (args.yes or os.environ.get("PBS_JOBID")
                           or os.environ.get("SLURM_JOB_ID")
                           or os.environ.get("UDALES_V1_CONFIRM")):
        print(f"preset '{preset.name}' needs {nprocs} MPI ranks and hours of wall "
              "time.\nSubmit tests/validation/nesting/submit_cx3.pbs, or pass --yes "
              "if you really mean to run it here.", file=sys.stderr)
        return 2
    rundir = args.rundir
    rundir.mkdir(parents=True, exist_ok=True)
    parent_dir = rundir / preset.parent_expnr
    child_dir = rundir / preset.child_expnr
    outdir = rundir / "analysis"

    first = STAGES.index(args.start_at)
    last = STAGES.index(args.stop_after)

    def wanted(stage: str) -> bool:
        i = STAGES.index(stage)
        return first <= i <= last

    print(preset.summary())
    print()
    (rundir / "preset_summary.txt").write_text(preset.summary() + "\n", encoding="ascii")

    t_all = time.time()
    if wanted("parent-case"):
        t = time.time()
        make_parent_case.build(rundir, preset, ibm_backend=args.ibm_backend)
        print(f"[run_v1] parent case built in {time.time() - t:.1f} s")

    if wanted("spinup"):
        t = time.time()
        run_solver(parent_dir, f"namoptions_spinup.{preset.parent_expnr}",
                   preset.nprocx * preset.nprocy, parent_dir / "spinup.log")
        print(f"[run_v1] parent spin-up finished in {time.time() - t:.1f} s")

    if wanted("production"):
        _set_startfile(parent_dir / f"namoptions.{preset.parent_expnr}",
                       _restart_file(parent_dir, preset.parent_expnr))
        t = time.time()
        run_solver(parent_dir, f"namoptions.{preset.parent_expnr}",
                   preset.nprocx * preset.nprocy, parent_dir / "production.log")
        print(f"[run_v1] parent production finished in {time.time() - t:.1f} s")

    if wanted("child-case"):
        t = time.time()
        make_child_case.build(parent_dir, rundir, preset, ibm_backend=args.ibm_backend)
        manifest = json.loads((child_dir / "manifest.json").read_text())
        print(f"[run_v1] child case built in {time.time() - t:.1f} s: "
              f"{manifest['n_parent_levels']} parent levels, "
              f"|Phi|/A = {manifest['flux_residual_after_correction']['max_abs_normalised']:.2e}")

    if wanted("child"):
        t = time.time()
        run_solver(child_dir, f"namoptions.{preset.child_expnr}",
                   preset.child_nprocx * preset.child_nprocy, child_dir / "child.log")
        print(f"[run_v1] child run finished in {time.time() - t:.1f} s")

    if wanted("analysis"):
        t = time.time()
        metrics = analyse.run(parent_dir, child_dir, outdir, preset,
                              make_plots=not args.no_plots)
        print(f"[run_v1] analysis finished in {time.time() - t:.1f} s\n")
        print(analyse.summary(metrics))
        print(f"\nresults in {outdir}")

    print(f"\n[run_v1] total {time.time() - t_all:.1f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
