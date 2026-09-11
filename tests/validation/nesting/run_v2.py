#!/usr/bin/env python3
"""End-to-end driver for the V2 nesting sweeps.

V2 is not a survey of the parameter space; it is a **falsification test** of the
explanation section 10.5 of ``docs/udales-nesting-design.md`` gives for the V1
resolved-TKE deficit.  That explanation -- too little fetch for the child to
regenerate its own energy-containing turbulence, rather than a defect in the
boundary treatment -- makes two opposed predictions:

    P1  zone width should barely move the deficit;
    P2  child domain size should move it a lot.

So the driver runs two arms out of **one** parent: a zone-width ramp at a fixed
child size, and a child-size ramp at a fixed zone.  The parent is not re-run.
It is the V1 ``converged`` parent, already on disk, and every child is therefore
driven by numerically identical forcing -- which is what makes a few per cent of
difference between children mean anything at all.  ``config.Sweep.validate``
checks that each point really can share that parent (grid, geometry, forcing,
schedule and the plaza window) rather than trusting the presets to agree.

Usage
-----
    # the production sweep, against the V1 converged parent
    python run_v2.py $EPHEMERAL/nesting-v2 \\
        --sweep v2 --parent-dir $EPHEMERAL/nesting-v1-converged/903 --yes

    # the smoke test: build a tiny parent + reference child, then sweep it
    python run_v1.py $EPHEMERAL/v2-tiny --preset tiny
    python run_v2.py $EPHEMERAL/v2-tiny --sweep v2-tiny \\
        --parent-dir $EPHEMERAL/v2-tiny/903

    # C0a: cadence ladder off the same converged parent, same reused reference
    python run_v2.py $EPHEMERAL/nesting-c0a --sweep c0 \\
        --parent-dir $EPHEMERAL/nesting-v1-converged/903 \\
        --reuse-dir $EPHEMERAL/nesting-v1-converged/904 --yes

    # C0b: the sweep's parent does not exist yet -- build it warm-started from
    # the converged parent's restart files, run it, then the six children
    python run_v2.py $EPHEMERAL/nesting-c0b --sweep c0b \\
        --parent-restart-dir $EPHEMERAL/nesting-v1-converged/903 --yes

Stages are ``parent``, ``child-case``, ``child``, ``analysis`` and ``summary``;
``--only`` restricts the run to named sweep points, so a single failed child can
be redone without touching the rest.  The ``parent`` stage only does anything
with ``--parent-restart-dir``: it builds and runs the sweep's parent when
``--parent-dir`` has no field dumps yet, and is skipped when it has.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import analyse
import make_child_case
import make_parent_case
import sweep_summary
from caselib import run_solver
from config import Sweep, SweepPoint, get_sweep

STAGES = ("parent", "child-case", "child", "analysis", "summary")

#: Written next to each point's plots.  Deliberately *not* ``v1_metrics.json``:
#: a V2 point is one child of a sweep, and the sweep summary is the deliverable.
METRICS_NAME = "v2_metrics.json"


def _reuse_dir(point: SweepPoint, parent_dir: Path,
               override: Optional[Path]) -> Path:
    """Where an already-run child of this sweep lives.

    Defaults to the sibling of the parent case directory named after the
    child's experiment number -- i.e. exactly where ``run_v1.py`` left it, so
    the V1 rundir needs no rearranging and nothing is copied.
    """
    if override is not None:
        return Path(override)
    return Path(parent_dir).parent / point.preset.child_expnr


def _check_reused(point: SweepPoint, casedir: Path) -> Dict[str, object]:
    """Refuse to reuse a run that is not the run this point describes.

    Reuse is the whole economy of V2 -- it is what keeps the reference point
    from costing another 40 minutes -- but it is also the one place where the
    comparison could silently stop being paired.  So the manifest the previous
    run left behind is checked against this point's preset before its dumps are
    believed.
    """
    manifest_path = casedir / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"point '{point.key}' is marked for reuse but {manifest_path} does not "
            "exist; run it, or pass --reuse-dir at the directory that holds it"
        )
    manifest = json.loads(manifest_path.read_text())
    p = point.preset
    problems = []
    if manifest.get("child_expnr") != p.child_expnr:
        problems.append(f"child_expnr {manifest.get('child_expnr')} != {p.child_expnr}")
    for key, want in (("child_i0", p.child_i0), ("child_j0", p.child_j0),
                      ("ustar", p.ustar), ("dpdx", p.dpdx),
                      ("nest_timeinterp", p.timeinterp),
                      ("stats_start", p.child_spinup),
                      ("init_from_parent", bool(p.init_from_parent))):
        if manifest.get(key) != want:
            problems.append(f"{key} = {manifest.get(key)!r} != {want!r}")
    # The boundary cadence.  Manifests written before C0 carry no 'cadence'
    # block; for those the median parent dt is the cadence (every level was
    # used), and adaptive stepping puts it a fraction of a step above nominal
    # (3.004 s for the converged run), hence the tolerance.
    cad = manifest.get("cadence")
    have = (float(cad["seconds"]) if isinstance(cad, dict)
            else float(manifest.get("parent_dt_median", float("nan"))))
    if not (abs(have - p.cadence) <= 0.05 * p.cadence):
        problems.append(f"boundary cadence {have!r} s != {p.cadence!r} s")
    if not list(casedir.glob(f"fielddump.???.???.{p.child_expnr}.nc")):
        problems.append("no field dumps")
    if problems:
        raise RuntimeError(
            f"the run at {casedir} is not the child point '{point.key}' describes:\n  "
            + "\n  ".join(problems)
        )
    return manifest


def _disk_estimate(sweep: Sweep) -> str:
    """Rough per-point output volume, so a full disk is a prediction not a surprise."""
    lines, total = [], 0.0
    for pt in sweep.points:
        if pt.reuse:
            continue
        p = pt.preset
        # Levels stored in the nesting file follow the boundary cadence, the
        # child's own dumps its dump interval; neither is the parent's.
        nt_nest = p.production / p.cadence
        nt_dump = p.production / p.child_dtdump
        # 12 slabs, double precision; the two x-faces carry jtot x ktot x nzone
        # points and the two y-faces itot x ktot x nzone, three components each.
        slab = 3 * 2 * (p.child_itot + p.child_jtot) * p.child_ktot * p.nzone * 8
        dump = 3 * p.child_itot * p.child_jtot * p.child_ktot * 4
        gb = (nt_nest * slab + nt_dump * dump) / 1e9
        total += gb
        lines.append(f"    {pt.key:8s} nesting {nt_nest * slab / 1e9:6.1f} GB "
                     f"({nt_nest:.0f} levels, held in RAM while cutting) + dumps "
                     f"{nt_dump * dump / 1e9:6.1f} GB")
    lines.append(f"    {'total':8s} {total:6.1f} GB")
    return "\n".join(lines)


def _parent_has_dumps(parent_dir: Path, expnr: str) -> bool:
    return bool(list(Path(parent_dir).glob(f"fielddump.???.???.{expnr}.nc")))


def run_parent(rundir: Path, sweep: Sweep, parent_dir: Path, restart_dir: Path,
               *, ibm_backend: str, timings: Dict[str, float]) -> None:
    """Build the sweep's parent warm-started from ``restart_dir`` and run it.

    Only when ``parent_dir`` has no field dumps yet: a parent that has run is
    reused, exactly as V2 reuses V1's, and never re-run by accident.
    ``parent_dir`` must be ``rundir / <parent expnr>`` -- ``make_parent_case``
    names the directory itself.
    """
    p = sweep.parent
    nr = p.parent_expnr
    if _parent_has_dumps(parent_dir, nr):
        print(f"    parent {nr} already has field dumps in {parent_dir}; not re-run")
        return
    if Path(parent_dir).resolve() != (Path(rundir) / nr).resolve():
        raise SystemExit(
            f"--parent-restart-dir builds the parent at {Path(rundir) / nr}, but "
            f"--parent-dir is {parent_dir}; drop --parent-dir or make them agree")
    t = time.time()
    make_parent_case.build(rundir, p, ibm_backend=ibm_backend,
                           restart_dir=restart_dir)
    timings["parent-case"] = time.time() - t
    info = json.loads((Path(parent_dir) / "preset.json").read_text())
    print(f"    parent case built in {timings['parent-case']:.1f} s, warm start "
          f"from {info['startfile']} (t = {p.t_start:g} s), {p.production:g} s "
          f"dumping every {p.dtdump:g} s ({p.production / p.dtdump:.0f} levels)")
    t = time.time()
    run_solver(Path(parent_dir), f"namoptions.{nr}", p.nprocx * p.nprocy,
               Path(parent_dir) / "production.log")
    timings["parent-production"] = time.time() - t
    print(f"    parent production finished in {timings['parent-production']:.1f} s")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundir", type=Path)
    ap.add_argument("--sweep", default="v2")
    ap.add_argument("--parent-dir", type=Path, default=None,
                    help="the parent case directory (default <rundir>/<parent expnr>)")
    ap.add_argument("--reuse-dir", type=Path, default=None,
                    help="where the already-run reference child lives "
                         "(default: alongside --parent-dir)")
    ap.add_argument("--only", default="",
                    help="comma-separated sweep point keys to act on")
    ap.add_argument("--parent-restart-dir", type=Path, default=None,
                    help="build and run the sweep's parent first, warm-started "
                         "from the restart files in this directory, unless "
                         "--parent-dir already holds its field dumps (C0b)")
    ap.add_argument("--start-at", default="parent", choices=STAGES)
    ap.add_argument("--stop-after", default="summary", choices=STAGES)
    ap.add_argument("--ibm-backend", default="auto")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--prune-nesting", action="store_true",
                    help="delete each child's nesting.inp after it has run; the "
                         "analysis does not need it, but a re-run does")
    ap.add_argument("--yes", action="store_true",
                    help="confirm a large run outside a batch job")
    args = ap.parse_args()

    sweep = get_sweep(args.sweep)
    parent_dir = (args.parent_dir if args.parent_dir is not None
                  else args.rundir / sweep.parent.parent_expnr)
    nprocs = sweep.parent.child_nprocx * sweep.parent.child_nprocy
    if nprocs > 8 and not (args.yes or os.environ.get("PBS_JOBID")
                           or os.environ.get("SLURM_JOB_ID")
                           or os.environ.get("UDALES_V1_CONFIRM")):
        print(f"sweep '{sweep.name}' needs {nprocs} MPI ranks and hours of wall time.\n"
              "Submit tests/validation/nesting/submit_cx3_v2.pbs, or pass --yes if "
              "you really mean to run it here.", file=sys.stderr)
        return 2

    rundir = args.rundir
    rundir.mkdir(parents=True, exist_ok=True)
    only = [k for k in args.only.split(",") if k]
    points: List[SweepPoint] = ([sweep.point(k) for k in only] if only
                                else list(sweep.points))
    first, last = STAGES.index(args.start_at), STAGES.index(args.stop_after)

    def wanted(stage: str) -> bool:
        return first <= STAGES.index(stage) <= last

    print(sweep.summary())
    print(f"\nparent case          {parent_dir}")
    print(f"acting on            {', '.join(p.key for p in points)}")
    print("estimated new output:")
    print(_disk_estimate(sweep))
    print()
    (rundir / "sweep_summary.txt").write_text(sweep.summary() + "\n", encoding="ascii")

    t_all = time.time()
    timings: Dict[str, Dict[str, float]] = {}
    if wanted("parent") and args.parent_restart_dir is not None:
        print(f"\n=== parent {sweep.parent.parent_expnr} ({sweep.parent.name})")
        timings["parent"] = {}
        run_parent(rundir, sweep, parent_dir, args.parent_restart_dir,
                   ibm_backend=args.ibm_backend, timings=timings["parent"])

    if not (parent_dir / f"namoptions.{sweep.parent.parent_expnr}").exists():
        raise SystemExit(
            f"no parent case at {parent_dir}. V2 reuses the V1 parent rather than "
            "re-running it -- point --parent-dir at it (or, for a sweep whose "
            "parent warm-starts from another run, pass --parent-restart-dir)."
        )

    # One accumulated parent Bundle per child window, shared across the points
    # that use it.  Both halves the analysis I/O and guarantees the zone arm is
    # measured against literally the same parent statistics.
    parent_cache: Dict = {}

    for pt in points:
        p = pt.preset
        casedir = (_reuse_dir(pt, parent_dir, args.reuse_dir) if pt.reuse
                   else rundir / p.child_expnr)
        timings[pt.key] = {}
        print(f"\n=== {pt.key} ({p.name}, {p.child_expnr}): {pt.note}")
        print(f"    case {casedir}"
              + ("  [reused, not re-run]" if pt.reuse else ""))

        if pt.reuse:
            _check_reused(pt, casedir)
        else:
            if wanted("child-case"):
                t = time.time()
                make_child_case.build(parent_dir, rundir, p,
                                      ibm_backend=args.ibm_backend)
                manifest = json.loads((casedir / "manifest.json").read_text())
                timings[pt.key]["child-case"] = time.time() - t
                print(f"    child case built in {timings[pt.key]['child-case']:.1f} s: "
                      f"{manifest['n_parent_levels']} parent levels, |Phi|/A = "
                      f"{manifest['flux_residual_after_correction']['max_abs_normalised']:.2e}")
            if wanted("child"):
                t = time.time()
                run_solver(casedir, f"namoptions.{p.child_expnr}",
                           p.child_nprocx * p.child_nprocy, casedir / "child.log")
                timings[pt.key]["child"] = time.time() - t
                print(f"    child run finished in {timings[pt.key]['child']:.1f} s")
                if args.prune_nesting:
                    nest = casedir / f"nesting.inp.{p.child_expnr}.nc"
                    if nest.exists():
                        gb = nest.stat().st_size / 1e9
                        nest.unlink()
                        print(f"    pruned {nest.name} ({gb:.1f} GB)")

        if wanted("analysis"):
            t = time.time()
            outdir = rundir / "analysis" / pt.key
            metrics = analyse.run(parent_dir, casedir, outdir, p,
                                  make_plots=not args.no_plots,
                                  parent_cache=parent_cache,
                                  common_block_cells=sweep.common_block_cells,
                                  metrics_name=METRICS_NAME)
            timings[pt.key]["analysis"] = time.time() - t
            print(f"    analysis finished in {timings[pt.key]['analysis']:.1f} s")
            print(analyse.summary(metrics))

    if wanted("summary"):
        summary = sweep_summary.build(rundir, sweep, metrics_name=METRICS_NAME)
        summary["timings_s"] = timings
        sweep_summary.write(rundir / "analysis", sweep, summary,
                            make_plots=not args.no_plots)
        print()
        print(sweep_summary.arm_tables(sweep, summary["rows"]))
        print()
        print(json.dumps(summary["verdict"], indent=2))
        print(f"\nresults in {rundir / 'analysis'}")

    print(f"\n[run_v2] total {time.time() - t_all:.1f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
