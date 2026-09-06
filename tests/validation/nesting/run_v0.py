#!/usr/bin/env python3
"""End-to-end driver for V0 -- validation at refinement ratio > 1.

Design ``docs/udales-nesting-design.md`` section 10.4 row V0.  V1 and V2 both
run parent and child on the same grid, so the conservative interpolation of
section 1.3 has never carried a running simulation; V0 is the first end-to-end
test of refinement and, simultaneously, the first end-to-end test of that
interpolation.  The suite is arranged so that the two can be told apart.

The child is the same in every point -- and it is the **V1 child**: the same
128 x 128 x 64 cells at 2 m over the same box, the same buildings, the same
3 + 9 cell zone, the same forcing, the same window.  The truth is the same too:
the V1 converged parent's own sub-region.  Only the grid the boundary data
arrives on changes, and it changes in two independent ways:

``filtered`` arm
    the boundary data is the fine reference's own field, box-filtered onto the
    coarse grid.  A *perfect* coarse parent -- it knows exactly what the fine
    run did at the scales it can hold, which no real coarse LES would.  Paired
    with the truth (same realisation), so it is V1 with exactly one variable
    changed, and it isolates the prolongation and the parent's filter scale.

``coarse`` arm
    the boundary data comes from a genuinely coarse LES of the same domain,
    run here.  The real use case.  Unpaired -- a different realisation -- and it
    carries the coarse run's own biases, which ``analyse_v0.parent_deficit``
    measures separately so that the child's error can be split from its
    parent's.

Running both and reporting them apart is the point: the difference between them
at one ratio is how much the child suffers from its parent genuinely not knowing
the small scales, as opposed to from the interpolation.

Usage
-----
    # the production suite, against the V1 converged run already on disk
    python run_v0.py $EPHEMERAL/nesting-v0 --suite v0 \\
        --reference-dir $EPHEMERAL/nesting-v1-converged/903 --yes

    # the smoke test: builds and runs its own tiny fine reference too
    python run_v0.py $EPHEMERAL/v0-tiny --suite v0-tiny

``--only`` restricts the run to named points, so one failed child can be redone
without touching the rest, and ``--start-at``/``--stop-after`` restrict the
stages.
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
import analyse_v0
import make_child_case
import make_parent_case
import run_v1
from caselib import run_solver
from config import Preset, RefinedPoint, RefinementSuite, get_suite
from make_child_case import DrivingParent

STAGES = ("reference-case", "reference-spinup", "reference-production",
          "driver-case", "driver-spinup", "driver-production",
          "child-case", "child", "analysis", "summary")

#: Deliberately not ``v1_metrics.json``: a V0 point is one child of a suite and
#: the suite table is the deliverable.
METRICS_NAME = "v0_metrics.json"


def _run_parent(rundir: Path, preset: Preset, casedir: Path, *,
                stages: Dict[str, bool], ibm_backend: str,
                label: str, timings: Dict[str, float]) -> None:
    """Build, spin up and produce one parent case, whichever grid it is on."""
    if stages["case"]:
        t = time.time()
        make_parent_case.build(rundir, preset, ibm_backend=ibm_backend)
        timings[f"{label}-case"] = time.time() - t
        print(f"    {label} case built in {timings[f'{label}-case']:.1f} s")
    if stages["spinup"]:
        t = time.time()
        run_solver(casedir, f"namoptions_spinup.{preset.parent_expnr}",
                   preset.nprocx * preset.nprocy, casedir / "spinup.log")
        timings[f"{label}-spinup"] = time.time() - t
        print(f"    {label} spin-up finished in {timings[f'{label}-spinup']:.1f} s")
    if stages["production"]:
        run_v1._set_startfile(
            casedir / f"namoptions.{preset.parent_expnr}",
            run_v1._restart_file(casedir, preset.parent_expnr))
        t = time.time()
        run_solver(casedir, f"namoptions.{preset.parent_expnr}",
                   preset.nprocx * preset.nprocy, casedir / "production.log")
        timings[f"{label}-production"] = time.time() - t
        print(f"    {label} production finished in "
              f"{timings[f'{label}-production']:.1f} s")


def _disk_estimate(suite: RefinementSuite, points: List[RefinedPoint]) -> str:
    """Rough new output volume, so a full disk is a prediction and not a surprise."""
    lines, total = [], 0.0
    for d in suite.drivers_to_run:
        nt = d.production / d.dtdump
        gb = nt * 3 * d.itot * d.jtot * d.ktot * 4 / 1e9
        total += gb
        lines.append(f"    parent {d.parent_expnr} @ {d.dx:g} m   dumps {gb:6.1f} GB")
    for pt in points:
        p = pt.child
        nt = p.production / p.dtdump
        slab = 3 * 2 * (p.child_itot + p.child_jtot) * p.child_ktot * p.nzone * 8
        dump = 3 * p.child_itot * p.child_jtot * p.child_ktot * 4
        init = 3 * p.child_itot * p.child_jtot * p.child_ktot * 8
        gb = (nt * (slab + dump) + init) / 1e9
        total += gb
        lines.append(f"    child  {p.child_expnr} ({pt.key:12s}) nesting "
                     f"{(nt * slab + init) / 1e9:6.1f} GB + dumps "
                     f"{nt * dump / 1e9:6.1f} GB")
    lines.append(f"    {'total':38s} {total:6.1f} GB")
    return "\n".join(lines)


def _tallest_spectrum(metrics: Dict[str, object]) -> Optional[str]:
    spec = metrics.get("spectra", {})
    if not spec:
        return None
    return max(spec, key=lambda k: spec[k]["z_over_h"])


def _row(point: RefinedPoint, metrics: Dict[str, object],
         v0: Dict[str, object]) -> Dict[str, object]:
    d = metrics["v2"]["tke_deficit"]
    ca = metrics["v2"]["criterion_a"]
    pdf = v0.get("parent_deficit", {})
    rt = v0.get("runtime", {})
    height = _tallest_spectrum(metrics)
    bands = (v0["spectra_across_parent_nyquist"][height]["bands"]
             if height else {})

    def band(name: str) -> Optional[float]:
        return bands.get(name, {}).get("ratio_of_sums")

    return {
        "key": point.key,
        "arm": point.arm,
        "refine": point.refine,
        "parent_dx_m": point.driver.dx,
        "parent_nyquist_m": point.parent_nyquist_wavelength,
        "child_expnr": point.expnr,
        "samples_parent": metrics["samples"]["parent"],
        "samples_child": metrics["samples"]["child"],
        "tke_deficit_above_2h": d["above"]["mean_relative"],
        "tke_spread_above_2h_median": d["above"]["median_spread"],
        "tke_canopy": d["canopy"]["mean_relative"],
        "parent_tke_deficit_above_2h": pdf.get("tke_mean_relative_above_2h"),
        "parent_u_rms_over_ustar": pdf.get("u_rms_difference"),
        "spectra_height": height,
        "band_parent_resolved": band("parent_resolved"),
        "band_parent_marginal": band("parent_marginal"),
        "band_sub_parent_filter": band("sub_parent_filter"),
        "band_contrast": (v0["spectra_across_parent_nyquist"][height]["contrast"]
                          if height else None),
        "u_rms_difference_over_ustar":
            metrics["profile_metrics"]["u_rms_difference_over_ustar"],
        "u_noise_floor_over_ustar":
            metrics["profile_metrics"]["u_noise_floor_over_ustar"],
        "criterion_a": ca["max_interior_umean_error_over_ustar"],
        "criterion_a_passes": ca["passes"],
        "tke_error_at_zone_edge":
            metrics["v2"]["tke_error_vs_fetch"]["mean_error_at_zone_edge"],
        "tke_error_at_max_fetch":
            metrics["v2"]["tke_error_vs_fetch"]["mean_error_at_max_fetch"],
        "runtime_max_phi": rt.get("phi", {}).get("max_abs"),
        "runtime_max_divmax": rt.get("divmax", {}).get("max"),
        "runtime_gradp_ratio": rt.get("gradp_ratio", {}).get("median_abs"),
        "prolongation_parent_divmax":
            (v0.get("prolongation_offline") or {}).get("parent_before_prolongation"),
        "prolongation_child_divmax":
            (v0.get("prolongation_offline") or {}).get("before_projection"),
        "ramp_resolved_by_parent": point.resolves_the_ramp,
    }


_TABLE_COLUMNS = (
    ("key", "point", "{}"),
    ("refine", "r", "{}"),
    ("arm", "arm", "{}"),
    ("parent_dx_m", "dx_P [m]", "{:g}"),
    ("tke_deficit_above_2h", "TKE z/h>2", "{:+.2%}"),
    ("parent_tke_deficit_above_2h", "of which parent's", "{:+.2%}"),
    ("tke_canopy", "TKE canopy", "{:+.2%}"),
    ("band_parent_resolved", "E: resolved", "{:.3f}"),
    ("band_parent_marginal", "marginal", "{:.3f}"),
    ("band_sub_parent_filter", "sub-filter", "{:.3f}"),
    ("criterion_a", "crit A [u*]", "{:.4f}"),
    ("runtime_max_phi", "max |Phi|", "{:.1e}"),
    ("runtime_max_divmax", "max divmax", "{:.1e}"),
)


def table(rows: List[Dict[str, object]]) -> str:
    head = [c[1] for c in _TABLE_COLUMNS]
    body = []
    for r in rows:
        cells = []
        for key, _, fmt in _TABLE_COLUMNS:
            v = r.get(key)
            cells.append("n/a" if v is None else fmt.format(v))
        body.append(cells)
    widths = [max(len(head[i]), *(len(b[i]) for b in body)) if body else len(head[i])
              for i in range(len(head))]
    out = ["| " + " | ".join(h.ljust(w) for h, w in zip(head, widths)) + " |",
           "|" + "|".join("-" * (w + 2) for w in widths) + "|"]
    for b in body:
        out.append("| " + " | ".join(c.ljust(w) for c, w in zip(b, widths)) + " |")
    return "\n".join(out)


def write_summary(outdir: Path, suite: RefinementSuite,
                  rows: List[Dict[str, object]],
                  timings: Dict[str, Dict[str, float]]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    payload = {"suite": suite.name, "rows": rows, "timings_s": timings,
               "arms": sorted({r["arm"] for r in rows}),
               "refinements": sorted({r["refine"] for r in rows})}
    # The quantity V0 exists to produce: at one ratio, how much of the child's
    # error survives when the parent is *perfect* (filtered) and how much more
    # appears when it is real (coarse).
    paired = {}
    for refine in payload["refinements"]:
        f = next((r for r in rows if r["refine"] == refine and r["arm"] == "filtered"),
                 None)
        c = next((r for r in rows if r["refine"] == refine and r["arm"] == "coarse"),
                 None)
        if f is None or c is None:
            continue
        def diff(key):
            a, b = c.get(key), f.get(key)
            return None if a is None or b is None else a - b
        paired[f"r{refine}"] = {
            "filtered_tke_deficit": f["tke_deficit_above_2h"],
            "coarse_tke_deficit": c["tke_deficit_above_2h"],
            "cost_of_a_real_parent": diff("tke_deficit_above_2h"),
            "filtered_criterion_a": f["criterion_a"],
            "coarse_criterion_a": c["criterion_a"],
            "criterion_a_difference": diff("criterion_a"),
            "parent_own_tke_deficit_filtered": f["parent_tke_deficit_above_2h"],
            "parent_own_tke_deficit_coarse": c["parent_tke_deficit_above_2h"],
            "coarse_les_error_beyond_the_filter":
                diff("parent_tke_deficit_above_2h"),
        }
    payload["filtered_vs_coarse"] = paired
    (outdir / "v0_summary.json").write_text(json.dumps(payload, indent=2) + "\n",
                                            encoding="ascii")
    keys = list(rows[0].keys()) if rows else []
    with (outdir / "v0_summary.csv").open("w", encoding="ascii", newline="\n") as fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join("" if r[k] is None else
                              (f"{r[k]:.9g}" if isinstance(r[k], float) else str(r[k]))
                              for k in keys) + "\n")
    md = ["# V0 -- refinement validation, suite '%s'" % suite.name, "",
          suite.summary(), "", table(rows), "",
          "## filtered vs coarse, at the same ratio", "",
          "```", json.dumps(paired, indent=2), "```", ""]
    (outdir / "v0_summary.md").write_text("\n".join(md), encoding="ascii")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundir", type=Path)
    ap.add_argument("--suite", default="v0-tiny")
    ap.add_argument("--reference-dir", type=Path, default=None,
                    help="an existing fine reference case (the V1 converged parent); "
                         "when given, the reference is neither built nor run")
    ap.add_argument("--only", default="",
                    help="comma-separated point keys to act on")
    ap.add_argument("--start-at", default="reference-case", choices=STAGES)
    ap.add_argument("--stop-after", default="summary", choices=STAGES)
    ap.add_argument("--ibm-backend", default="auto")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--prune-nesting", action="store_true",
                    help="delete each child's nesting.inp once it has run; the "
                         "analysis does not need it, only a re-run does")
    ap.add_argument("--yes", action="store_true",
                    help="confirm a large run outside a batch job")
    args = ap.parse_args()

    suite = get_suite(args.suite)
    ref = suite.reference
    nprocs = ref.child_nprocx * ref.child_nprocy
    if nprocs > 8 and not (args.yes or os.environ.get("PBS_JOBID")
                           or os.environ.get("SLURM_JOB_ID")
                           or os.environ.get("UDALES_V1_CONFIRM")):
        print(f"suite '{suite.name}' needs {nprocs} MPI ranks and hours of wall "
              "time.\nSubmit tests/validation/nesting/submit_cx3_v0.pbs, or pass "
              "--yes if you really mean to run it here.", file=sys.stderr)
        return 2

    rundir = args.rundir
    rundir.mkdir(parents=True, exist_ok=True)
    only = [k for k in args.only.split(",") if k]
    points: List[RefinedPoint] = ([suite.point(k) for k in only] if only
                                  else list(suite.points))
    first, last = STAGES.index(args.start_at), STAGES.index(args.stop_after)

    def wanted(stage: str) -> bool:
        return first <= STAGES.index(stage) <= last

    ref_dir = (Path(args.reference_dir) if args.reference_dir is not None
               else rundir / ref.parent_expnr)

    print(suite.summary())
    print(f"\nreference (truth)    {ref_dir}"
          + ("  [existing, not re-run]" if args.reference_dir is not None else ""))
    print(f"acting on            {', '.join(p.key for p in points)}")
    print("estimated new output:")
    print(_disk_estimate(suite, points))
    print()
    (rundir / "v0_summary.txt").write_text(suite.summary() + "\n", encoding="ascii")

    timings: Dict[str, Dict[str, float]] = {}
    t_all = time.time()

    # ---- the fine reference ------------------------------------------------ #
    if args.reference_dir is None:
        timings["reference"] = {}
        print(f"=== fine reference {ref.name} ({ref.parent_expnr})")
        _run_parent(rundir, ref, ref_dir, ibm_backend=args.ibm_backend,
                    label="reference", timings=timings["reference"],
                    stages={"case": wanted("reference-case"),
                            "spinup": wanted("reference-spinup"),
                            "production": wanted("reference-production")})
    if not (ref_dir / f"namoptions.{ref.parent_expnr}").exists():
        raise SystemExit(
            f"no fine reference case at {ref_dir}.  V0 measures every child "
            "against a fine, unnested run -- point --reference-dir at one, or "
            "let this driver build it."
        )

    # ---- the coarse driving parents ---------------------------------------- #
    needed = {pt.driver.parent_expnr for pt in points if pt.runs_driver}
    for driver in suite.drivers_to_run:
        if driver.parent_expnr not in needed:
            continue
        casedir = rundir / driver.parent_expnr
        timings[driver.parent_expnr] = {}
        print(f"\n=== coarse driving parent {driver.parent_expnr} "
              f"({driver.itot}x{driver.jtot}x{driver.ktot} @ {driver.dx:g} m)")
        _run_parent(rundir, driver, casedir, ibm_backend=args.ibm_backend,
                    label="driver", timings=timings[driver.parent_expnr],
                    stages={"case": wanted("driver-case"),
                            "spinup": wanted("driver-spinup"),
                            "production": wanted("driver-production")})

    # ---- the children ------------------------------------------------------ #
    parent_cache: Dict = {}
    rows: List[Dict[str, object]] = []
    for pt in points:
        p = pt.child
        casedir = rundir / p.child_expnr
        driving_dir = ref_dir if pt.arm == "filtered" else rundir / pt.driver_expnr
        timings.setdefault(pt.key, {})
        print(f"\n=== {pt.key}: r = {pt.refine}, {pt.arm} arm, child {p.child_expnr}")
        print(f"    driven by {driving_dir}")

        if wanted("child-case"):
            t = time.time()
            driving = DrivingParent.refined(driving_dir, pt)
            try:
                make_child_case.build(driving_dir, rundir, p,
                                      ibm_backend=args.ibm_backend, driving=driving)
            finally:
                driving.dump.close()
            manifest = json.loads((casedir / "manifest.json").read_text())
            timings[pt.key]["child-case"] = time.time() - t
            print(f"    child case built in {timings[pt.key]['child-case']:.1f} s: "
                  f"{manifest['n_parent_levels']} parent levels from "
                  f"{manifest['refinement']['source']}, |Phi|/A = "
                  f"{manifest['flux_residual_after_correction']['max_abs_normalised']:.2e}")
            ic = manifest["initial_condition_divmax"]
            if ic.get("parent_before_prolongation") is not None:
                print(f"    prolongation: parent divmax "
                      f"{ic['parent_before_prolongation']:.3e} -> child "
                      f"{ic['before_projection']:.3e}")

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
            metrics = analyse.run(ref_dir, casedir, outdir, p,
                                  make_plots=not args.no_plots,
                                  parent_cache=parent_cache,
                                  metrics_name=METRICS_NAME)
            v0 = analyse_v0.augment(outdir, casedir, metrics, pt,
                                    metrics_name=METRICS_NAME)
            timings[pt.key]["analysis"] = time.time() - t
            print(f"    analysis finished in {timings[pt.key]['analysis']:.1f} s")
            print(analyse_v0.summary(v0, metrics))
            rows.append(_row(pt, metrics, v0))

    if wanted("summary"):
        if not rows:
            # Re-read whatever is on disk, so the table can be remade without
            # redoing the analysis.
            for pt in points:
                path = rundir / "analysis" / pt.key / METRICS_NAME
                if not path.exists():
                    continue
                metrics = json.loads(path.read_text())
                if "v0" in metrics:
                    rows.append(_row(pt, metrics, metrics["v0"]))
        if rows:
            write_summary(rundir / "analysis", suite, rows, timings)
            print()
            print(table(rows))
            print(f"\nresults in {rundir / 'analysis'}")
        else:
            print("\nno per-point metrics found; nothing to summarise")

    print(f"\n[run_v0] total {time.time() - t_all:.1f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
