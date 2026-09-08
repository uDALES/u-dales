#!/usr/bin/env python3
"""End-to-end driver for V6 -- does mass drift over long nested runs?

(design section 10.4, row V6).  Not a physics question: with a boundary that
goes steady partway through and stays there, any later growth in ``divtot``,
``divmax`` or the nesting flux residual Phi is numerical drift and nothing
else.

    parent case -> spin-up -> production (a handful of dumped levels)
        -> child case (cut from those dumps, config.V6: nest_lendabort = .false.)
        -> patch RUN.runtime, NESTING.nest_statint, NAMCHECKSIM.tcheck
        -> long nested run, boundary frozen on the last stored level for most of it
        -> analysis (analyse_v6.py): a trend fit per diagnostic, pass/fail

``--runtime`` and ``--stat-interval`` are not on :class:`config.Preset`
because the right values depend on a measured step rate (see
``config.V6``'s docstring): a short **probe** run (small ``--runtime``) times
itself, and that measurement sizes the full run's ``--runtime`` and
``--stat-interval`` -- both patched into the identical case with
``caselib.set_namoption`` rather than by rebuilding it differently, so the
probe and the real run share the same code path up to that point.

Usage
-----
    # probe: measure the step rate on a login node, well past the 8-level record
    python run_v6.py $EPHEMERAL/v6-probe --runtime 300

    # the real run, sized from the probe (submit_cx3_v6.pbs does this)
    python run_v6.py $EPHEMERAL/v6 --runtime 45000 --stat-interval 500 --yes

    # re-run just the child after tuning --runtime, reusing the built parent
    python run_v6.py <rundir> --start-at child-case --runtime 45000 \\
        --stat-interval 500 --yes

Stages: parent-case, spinup, production, child-case, child, analysis.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import analyse_v6
import make_child_case
import make_parent_case
from caselib import read_namoption, run_solver, set_namoption
from config import get_preset
from run_v1 import _restart_file, _set_startfile

STAGES = ("parent-case", "spinup", "production", "child-case", "child", "analysis")

#: Above this, refuse to run interactively on a login node without --yes: a
#: multi-hour job belongs in a batch job (submit_cx3_v6.pbs), not a login-node
#: foreground process.  A probe run (a few hundred seconds of simulated time)
#: is well under it; the production runtime (tens of thousands of seconds) is
#: deliberately not.
_INTERACTIVE_RUNTIME_LIMIT_S = 1000.0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("rundir", type=Path)
    parser.add_argument("--preset", default="v6")
    parser.add_argument("--start-at", default="parent-case", choices=STAGES)
    parser.add_argument("--stop-after", default="analysis", choices=STAGES)
    parser.add_argument("--ibm-backend", default="auto")
    parser.add_argument("--runtime", type=float, default=None,
                        help="child RUN.runtime [s], patched in after child-case "
                             "is (re)built. Omit to keep make_child_case.build's "
                             "own default -- short, stops before the record ends, "
                             "so the freeze is never exercised.")
    parser.add_argument("--stat-interval", type=float, default=None,
                        help="NESTING.nest_statint AND NAMCHECKSIM.tcheck [s], "
                             "patched in after child-case is built. Omit to keep "
                             "the Fortran/formula defaults.")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--yes", action="store_true",
                        help="confirm a long run outside a batch job")
    args = parser.parse_args()

    preset = get_preset(args.preset)

    if args.runtime is not None and args.runtime > _INTERACTIVE_RUNTIME_LIMIT_S \
            and not (args.yes or os.environ.get("PBS_JOBID")
                     or os.environ.get("SLURM_JOB_ID")):
        print(f"--runtime {args.runtime:g} s is a long run (limit for an "
              f"unconfirmed interactive run is {_INTERACTIVE_RUNTIME_LIMIT_S:g} s).\n"
              "Submit tests/validation/nesting/submit_cx3_v6.pbs, or pass --yes "
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
        return first <= STAGES.index(stage) <= last

    print(preset.summary())
    print()

    t_all = time.time()
    if wanted("parent-case"):
        t = time.time()
        make_parent_case.build(rundir, preset, ibm_backend=args.ibm_backend)
        print(f"[run_v6] parent case built in {time.time() - t:.1f} s")

    if wanted("spinup"):
        t = time.time()
        run_solver(parent_dir, f"namoptions_spinup.{preset.parent_expnr}",
                   preset.nprocx * preset.nprocy, parent_dir / "spinup.log")
        print(f"[run_v6] parent spin-up finished in {time.time() - t:.1f} s")

    if wanted("production"):
        _set_startfile(parent_dir / f"namoptions.{preset.parent_expnr}",
                       _restart_file(parent_dir, preset.parent_expnr))
        t = time.time()
        run_solver(parent_dir, f"namoptions.{preset.parent_expnr}",
                   preset.nprocx * preset.nprocy, parent_dir / "production.log")
        print(f"[run_v6] parent production finished in {time.time() - t:.1f} s")

    if wanted("child-case"):
        t = time.time()
        make_child_case.build(parent_dir, rundir, preset, ibm_backend=args.ibm_backend)
        namoptions = child_dir / f"namoptions.{preset.child_expnr}"
        manifest = json.loads((child_dir / "manifest.json").read_text())
        record_runtime = manifest["runtime"]
        if args.runtime is not None:
            set_namoption(namoptions, "runtime", args.runtime)
        if args.stat_interval is not None:
            set_namoption(namoptions, "nest_statint", args.stat_interval)
            set_namoption(namoptions, "tcheck", args.stat_interval)
        actual_runtime = args.runtime if args.runtime is not None else record_runtime
        past_record = actual_runtime > record_runtime
        print(f"[run_v6] child case built in {time.time() - t:.1f} s: "
              f"{manifest['n_parent_levels']} parent levels, record margin ends "
              f"{record_runtime:.1f} s in; runtime set to {actual_runtime:.1f} s "
              f"({'PAST the record -- freeze will be exercised' if past_record else 'still inside the record -- nest_lendabort never triggers'})")
        if preset.nest_lendabort and past_record:
            print("[run_v6] WARNING: nest_lendabort = .true. on this preset and "
                  "runtime is past the record -- the run will ABORT, not freeze.",
                  file=sys.stderr)

    if wanted("child"):
        t = time.time()
        run_solver(child_dir, f"namoptions.{preset.child_expnr}",
                   preset.child_nprocx * preset.child_nprocy, child_dir / "child.log")
        print(f"[run_v6] child run finished in {time.time() - t:.1f} s")

    exit_code = 0
    if wanted("analysis"):
        t = time.time()
        # The child's own namoptions is the source of truth for how long it was
        # actually asked to run -- read it back rather than relying on
        # args.runtime, which is None (and the manifest's record_runtime does
        # not apply) whenever --start-at skips the child-case stage.
        namoptions = child_dir / f"namoptions.{preset.child_expnr}"
        expected_runtime_s = None
        if namoptions.exists():
            raw = read_namoption(namoptions, "runtime")
            if raw is not None:
                expected_runtime_s = float(raw)
        result = analyse_v6.run(child_dir, outdir, make_plots=not args.no_plots,
                                expected_runtime_s=expected_runtime_s)
        print(f"[run_v6] analysis finished in {time.time() - t:.1f} s\n")
        print(analyse_v6.summary(result))
        print(f"\nresults in {outdir}")
        if result["overall_verdict"] != "PASS":
            print(f"[run_v6] analysis verdict is {result['overall_verdict']}, not PASS "
                  "-- reporting failure", file=sys.stderr)
            exit_code = 1

    print(f"\n[run_v6] total {time.time() - t_all:.1f} s")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
