#!/usr/bin/env python3
"""End-to-end driver for the V3 and V4 geometry-mismatch experiments.

    periodic runs (reference, parent) -> child cases -> child runs
        -> per-child analysis -> cross-child summary

Both experiments go through this one driver because they share every stage; what
differs is which periodic runs there are and what the per-child analysis reduces
to, and both of those live in :mod:`presets_geometry` and
:mod:`analyse_geometry`.

    # V3: the adjustment length and the standoff comparison
    python run_geometry.py $EPHEMERAL/nesting-v3 --experiment v3 --yes

    # V4: a mismatched parent layout, against the V1 converged child
    python run_geometry.py $EPHEMERAL/nesting-v4 --experiment v4 \\
        --baseline-dir $EPHEMERAL/nesting-v1-converged/904 --yes

    # the smoke tests
    python run_geometry.py $EPHEMERAL/v3-tiny --experiment v3-tiny
    python run_geometry.py $EPHEMERAL/v4-tiny --experiment v4-tiny

Stages are ``periodic-case``, ``periodic-spinup``, ``periodic-production``,
``periodic-stats``, ``child-case``, ``child``, ``analysis`` and ``summary``;
``--only`` restricts the children acted on, so one failed child can be redone.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np

import analyse
import analyse_geometry as ag
import make_child_case
import make_geometry_cases
import run_v1
from caselib import run_solver
from config import Preset, get_preset
from presets_geometry import ChildRun, Experiment, get_experiment

STAGES = ("periodic-case", "periodic-spinup", "periodic-production",
          "periodic-stats", "child-case", "child", "analysis", "summary")

#: Written next to each child's plots.  The V1/V2 child-versus-parent block is
#: written alongside it under a name that says what it is: for V3 and V4 the
#: parent sub-region is not a reference, so calling that file ``v1_metrics.json``
#: would invite exactly the misreading ``analyse_geometry``'s docstring warns
#: about.
LEGACY_NAME = "child_vs_parent_metrics.json"


# --------------------------------------------------------------------------- #
# Periodic runs
# --------------------------------------------------------------------------- #


def _periodic_stats_path(rundir: Path, key: str) -> Path:
    return rundir / "analysis" / f"periodic_{key}.json"


def run_periodic(rundir: Path, exp: Experiment, wanted, ibm_backend: str,
                 timings: Dict[str, Dict[str, float]]) -> None:
    """Build, spin up, run and reduce every periodic run, in order.

    Order matters: a run with ``bulk_from`` cannot be built until the run it
    names has been reduced, because its volume-flow-rate target *is* that
    reduction.  V3's flat parent is the case -- it is asked to carry the
    canopy's own bulk velocity, and that velocity is a measurement of the
    reference run, not a guess (see ``presets_geometry.V3_PARENT``).
    """
    for r in exp.periodic:
        p = r.preset
        casedir = rundir / p.parent_expnr
        t = timings.setdefault(f"periodic:{r.key}", {})
        print(f"\n=== periodic '{r.key}' ({p.name}, {p.parent_expnr}): {r.note}")

        rate: Optional[float] = p.uflowrate
        source = "preset"
        if r.bulk_from is not None:
            stats_path = _periodic_stats_path(rundir, r.bulk_from)
            if stats_path.exists():
                measured = json.loads(stats_path.read_text())["profiles"]["bulk_u"]
                rate, source = float(measured), f"measured from '{r.bulk_from}'"
            else:
                print(f"    [warn] {stats_path} does not exist; falling back to the "
                      f"preset's uflowrate = {rate}")
                source = "preset fallback -- the reference run has not been reduced"
            if rate is not None and p.uflowrate is not None:
                ratio = rate / p.uflowrate
                if not (0.5 < ratio < 2.0):
                    raise RuntimeError(
                        f"the measured bulk velocity {rate:.3f} m/s is {ratio:.2f} "
                        f"times the preset's {p.uflowrate:g} m/s.  That is too far "
                        "apart to be a calibration; something is wrong with the "
                        "reference run rather than with the guess.")
        if rate is not None:
            print(f"    uflowrate = {rate:.4f} m/s ({source})")

        if wanted("periodic-case"):
            t0 = time.time()
            make_geometry_cases.build_periodic(rundir, p, uflowrate=rate,
                                               ibm_backend=ibm_backend)
            t["case"] = time.time() - t0
            print(f"    case built in {t['case']:.1f} s")
        if wanted("periodic-spinup"):
            t0 = time.time()
            run_solver(casedir, f"namoptions_spinup.{p.parent_expnr}",
                       p.nprocx * p.nprocy, casedir / "spinup.log")
            t["spinup"] = time.time() - t0
            print(f"    spin-up finished in {t['spinup']:.1f} s")
        if wanted("periodic-production"):
            run_v1._set_startfile(
                casedir / f"namoptions.{p.parent_expnr}",
                run_v1._restart_file(casedir, p.parent_expnr))
            t0 = time.time()
            run_solver(casedir, f"namoptions.{p.parent_expnr}",
                       p.nprocx * p.nprocy, casedir / "production.log")
            t["production"] = time.time() - t0
            print(f"    production finished in {t['production']:.1f} s")
        if wanted("periodic-stats"):
            t0 = time.time()
            st, fluid, info = ag.accumulate_periodic(
                casedir, p, p.child_spinup, p.production, r.key)
            prof = ag.domain_profiles(st, fluid, p)
            out = {"key": r.key, "preset": p.name, "role": p.role,
                   "expnr": p.parent_expnr, "parent_layout": p.parent_layout,
                   "forcing": {"uflowrate": rate, "source": source,
                               "dpdx": 0.0 if rate is not None else p.dpdx},
                   "samples": info, "profiles": prof}
            path = _periodic_stats_path(rundir, r.key)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(out, indent=2) + "\n", encoding="ascii")
            t["stats"] = time.time() - t0
            print(f"    statistics in {t['stats']:.1f} s: bulk <u> = "
                  f"{prof['bulk_u']:.4f} m/s (half-window spread "
                  f"{prof['bulk_series_halfspread']:.4f}), canopy <u> = "
                  f"{prof['u_canopy']:.4f} m/s, <u'w'> at "
                  f"z = {prof['z_roof_m']:g} m = {prof['uw_at_roof']:+.5f} m2/s2")


# --------------------------------------------------------------------------- #
# Children
# --------------------------------------------------------------------------- #


def _fixed_station_residuals(blocks, equilibrium, preset,
                             stations_h=(5.0, 10.0, 20.0)) -> Dict[str, object]:
    """Relative canopy-velocity error at fixed distances from the inner zone edge.

    The cross-standoff comparison that needs no adjustment length to exist: at a
    given *physical station* -- the same distance into the domain for every
    child -- how far from equilibrium is the canopy?  Section 9.4's claim is
    that the 0-cell standoff should be closest at every station.  Stations are
    quoted in building heights and chosen to exist for every child in the sweep.
    """
    rows = [b for b in blocks if b["has_cube"]]
    out: Dict[str, object] = {"stations_h": list(stations_h)}
    eq = equilibrium.get("u_canopy")
    if not rows or eq in (None, 0.0):
        return out
    f = np.array([b["fetch_from_zone_h"] for b in rows], dtype=float)
    v = np.array([np.nan if b["u_canopy"] is None else b["u_canopy"] for b in rows])
    rel = np.abs(v - eq) / abs(eq)
    sp = np.array([np.nan if b.get("u_canopy_spread") is None
                   else b["u_canopy_spread"] for b in rows], dtype=float)
    out["at_station"] = {
        f"{s:g}h": (None if s < f[0] - 1.0e-9 or s > f[-1] + 1.0e-9
                    else float(np.interp(s, f, rel)))
        for s in stations_h}
    out["spread_at_station"] = {
        f"{s:g}h": (None if s < f[0] - 1.0e-9 or s > f[-1] + 1.0e-9
                    else float(np.interp(s, f, sp)))
        for s in stations_h}
    out["median_spread"] = float(np.nanmedian(sp)) if np.any(np.isfinite(sp)) else None
    return out


def analyse_child_v3(rundir: Path, exp: Experiment, c: ChildRun, casedir: Path,
                     parent_dir: Path, outdir: Path, legacy: bool) -> Dict:
    p = c.preset
    eq_path = _periodic_stats_path(rundir, exp.equilibrium_key)
    par_path = _periodic_stats_path(rundir, c.parent_key)
    for path in (eq_path, par_path):
        if not path.exists():
            raise FileNotFoundError(
                f"{path} is missing; run the periodic-stats stage first -- the "
                "adjustment length is measured against those profiles")
    equilibrium = json.loads(eq_path.read_text())["profiles"]
    imposed = json.loads(par_path.read_text())["profiles"]

    st, fluid, info = ag.accumulate_child(casedir, p, (), c.key)
    blocks = {scope: ag.streamwise_blocks(st, fluid, p, y_scope=scope)
              for scope in ("core", "canopy")}
    adjustment = ag.adjustment_length(blocks["core"], equilibrium, p)
    adjustment_full = ag.adjustment_length(blocks["canopy"], equilibrium, p)
    ii, jj = analyse.interior_indices(p)
    metrics = {
        "experiment": exp.name, "kind": "v3", "child_key": c.key,
        "preset": p.name, "note": c.note,
        "samples": {"child": info["n_samples"], "window_s": info["window_s"]},
        "configuration": {
            "standoff_cells": p.standoff_cells,
            "standoff_m": p.standoff_m,
            "first_row_fetch_m": p.first_row_fetch_m,
            "first_row_fetch_h": p.first_row_fetch_m / p.building_height,
            "n_rows": p.n_rows,
            "n_cubes": int(len(p.child_cube_centres())),
            "canopy_x_range_m": list(p.canopy_x_range),
            "canopy_y_range_m": list(p.canopy_y_range),
            "y_core_range_m": list(p.y_core_range),
            "zone_cells": p.zone_cells,
            "streamwise_interior_m": p.child_xlen - 2 * (p.guardwidth + p.zonewidth),
            "streamwise_interior_h": (p.child_xlen - 2 * (p.guardwidth + p.zonewidth))
                                     / p.building_height,
            "spanwise_interior_h": (p.child_ylen - 2 * (p.guardwidth + p.zonewidth))
                                   / p.building_height,
            "building_free_zone": p.building_free_zone,
            "nest_lparentgeom": not p.building_free_zone,
            "n_interior_cubes": int(len(p.cubes_in_analysis_interior())),
        },
        "time_series": {"times_s": st.times.tolist(),
                        "bulk_u": st.bulk_series.tolist(),
                        "note": "the child's own fluid-masked, depth-averaged "
                                "<u> at each sampled time; a flat trace over the "
                                "statistics window is what child_spinup was "
                                "chosen to leave"},
        "equilibrium": equilibrium,
        "imposed": imposed,
        "mismatch": {
            "bulk_relative": (imposed["bulk_u"] / equilibrium["bulk_u"] - 1.0
                              if equilibrium["bulk_u"] else None),
            "canopy_u_relative": (imposed["u_canopy"] / equilibrium["u_canopy"] - 1.0
                                  if equilibrium["u_canopy"] else None),
            "uw_at_roof_relative": (imposed["uw_at_roof"] / equilibrium["uw_at_roof"]
                                    - 1.0 if equilibrium["uw_at_roof"] else None),
            "note": "how far the imposed (parent) state is from the canopy's "
                    "equilibrium: the mismatch the child has to work off",
        },
        "blocks": blocks,
        "adjustment": adjustment,
        "adjustment_full_width": adjustment_full,
        "fixed_stations": _fixed_station_residuals(blocks["core"], equilibrium, p),
        "ibl": ag.ibl_depth(blocks["core"], imposed, p, p.ustar),
        "interior_profiles": ag.domain_profiles(st, fluid, p, ii, jj),
    }
    ag.write_v3(outdir, metrics)
    if legacy:
        analyse.run(parent_dir, casedir, outdir, p, make_plots=False,
                    metrics_name=LEGACY_NAME)
    return metrics


def analyse_child_v4(rundir: Path, exp: Experiment, c: ChildRun, casedir: Path,
                     parent_dir: Path, outdir: Path, legacy: bool,
                     baseline: Optional[Tuple[Path, Preset]],
                     cache: Dict) -> Dict:
    p = c.preset
    kk = [int(round(z / p.dz - 0.5)) for z in p.spectra_heights]
    kk = [k for k in kk if 0 <= k < p.child_ktot - 1]
    st, fluid, info = ag.accumulate_child(casedir, p, kk, c.key)

    metrics: Dict[str, object] = {
        "experiment": exp.name, "kind": "v4", "child_key": c.key,
        "preset": p.name, "note": c.note,
        "samples": {"child": info["n_samples"], "window_s": info["window_s"]},
        "configuration": {
            "parent_layout": p.parent_layout,
            "child_layout": p.child_layout,
            "n_cubes_parent": int(len(p.cube_centres())),
            "n_cubes_child": int(len(p.child_cube_centres())),
            "n_cubes_dropped_for_the_zone": int(len(p.child_cubes_removed())),
            "building_free_zone": p.building_free_zone,
            "nest_lparentgeom": not p.building_free_zone,
            "interior_cells": p.interior_cells,
            "interior_extent_h": p.interior_extent_h,
        },
        "interior_profiles": ag.domain_profiles(st, fluid, p,
                                                *analyse.interior_indices(p)),
        "time_series": {"times_s": st.times.tolist(),
                        "bulk_u": st.bulk_series.tolist()},
    }

    if baseline is not None:
        base_dir, base_preset = baseline
        cached = cache.get(str(base_dir))
        if cached is None:
            bst, bfluid, binfo = ag.accumulate_child(base_dir, base_preset, kk,
                                                     "baseline")
            cache[str(base_dir)] = (bst, bfluid, binfo)
        else:
            bst, bfluid, binfo = cached
            print("    [analysis] reusing the accumulated baseline child")
        if bfluid.shape != fluid.shape or not np.array_equal(bfluid, fluid):
            raise RuntimeError(
                "the baseline child's solid mask is not the mismatched child's.  "
                "V4 rests on the two children being geometrically identical -- "
                f"{int((bfluid ^ fluid).sum())} cells differ.  Check that the "
                "baseline directory really holds the V1 converged child.")
        metrics["samples"]["baseline"] = binfo["n_samples"]
        metrics["baseline"] = {"dir": str(Path(base_dir).resolve()),
                               "preset": base_preset.name,
                               "samples": binfo}
        metrics["comparison"] = ag.compare_children(
            st, bst, p, fluid, kk, c.key, "baseline")
    ag.write_v4(outdir, metrics)
    if legacy:
        analyse.run(parent_dir, casedir, outdir, p, make_plots=False,
                    metrics_name=LEGACY_NAME)
    return metrics


# --------------------------------------------------------------------------- #
# Cross-child summaries
# --------------------------------------------------------------------------- #


def summarise_v3(exp: Experiment, per_child: Dict[str, Dict]) -> Dict:
    """The standoff table, and the verdict on section 9.4's standoff claim.

    Section 9.4 says a building-free standoff is **actively counterproductive**:
    the flow adjusts once to the ground and then again to the canopy, so
    starting the buildings at the inner zone edge should give a *shorter* total
    adjustment.  Turned into something that can fail, that is two predictions:

    ``P-a``  the adjustment length measured **from the inner edge of the zone**
             does not decrease as the standoff grows;
    ``P-b``  at a fixed station measured from the zone edge, the canopy of the
             0-cell child is at least as close to equilibrium as any other's;
    ``P-c``  the adjustment length measured **from the first building face** does
             not decrease as the standoff grows.

    **P-c is the sharp one and P-a is nearly free**, which is worth being blunt
    about.  A standoff moves the canopy downstream, so the canopy's adjustment
    trivially finishes later measured from the zone edge -- by at least the
    standoff length -- whatever the physics.  The non-trivial content of "one
    adjustment instead of two" is that the *second* adjustment is slower than the
    single one would have been: the flow reaching the canopy has already
    equilibrated with the ground and has to be reworked.  That is exactly
    ``adjustment_from_first_row``.  If it comes out *shorter* behind a standoff
    -- if a decelerated approach flow makes the canopy adjustment quicker -- then
    the standoff is not counterproductive in the way section 9.4 claims; it
    merely costs its own length, which is a much weaker statement than the design
    makes.

    A refutation is a standoff that reaches equilibrium *sooner* from the zone
    edge, or that is closer to equilibrium at a fixed station, by more than the
    sampling spread.  Neither prediction is evaluated with a tolerance chosen
    afterwards: the spread is each child's own half-window spread of the same
    quantity, which is emitted per block.
    """
    rows = []
    for key, m in per_child.items():
        cfg, adj = m["configuration"], m["adjustment"]
        u = adj.get("u_canopy", {}) if isinstance(adj.get("u_canopy"), dict) else {}
        eq = u.get("vs_equilibrium", {}) if isinstance(u.get("vs_equilibrium"), dict) else {}
        last = u.get("vs_last_row", {}) if isinstance(u.get("vs_last_row"), dict) else {}
        w = adj.get("uw_at_roof", {}) if isinstance(adj.get("uw_at_roof"), dict) else {}
        weq = w.get("vs_equilibrium", {}) if isinstance(w.get("vs_equilibrium"), dict) else {}
        st = m["fixed_stations"].get("at_station", {})
        d = [x for x in m["ibl"]["delta_i_over_h"] if x is not None]
        rows.append({
            "key": key,
            "standoff_cells": cfg["standoff_cells"],
            "standoff_m": cfg["standoff_m"],
            "first_row_fetch_h": cfg["first_row_fetch_h"],
            "n_rows": cfg["n_rows"],
            "adj_u_from_zone_h": (None if eq.get("settled_only_at_the_last_row")
                                  else eq.get("adjustment_from_zone_h")),
            "adj_u_settled_only_at_last_row": eq.get("settled_only_at_the_last_row"),
            "adj_u_from_row1_h": eq.get("adjustment_from_first_row_h"),
            "adj_u_selfref_from_zone_h": last.get("adjustment_from_zone_h"),
            "adj_uw_from_zone_h": weq.get("adjustment_from_zone_h"),
            "err_row1": eq.get("error_at_first_row"),
            "err_last": eq.get("error_at_last_row"),
            "max_fetch_h": eq.get("max_fetch_from_zone_h"),
            "residual_at": {k: v for k, v in st.items()},
            "median_spread": m["fixed_stations"].get("median_spread"),
            "ibl_first_h": d[0] if d else None,
            "ibl_last_h": d[-1] if d else None,
        })
    rows.sort(key=lambda r: r["standoff_cells"])

    def series(name):
        return [(r["standoff_cells"], r[name]) for r in rows if r[name] is not None]

    verdict: Dict[str, object] = {}
    sa = series("adj_u_from_zone_h")
    if len(sa) >= 2:
        base_a = sa[0][1]
        verdict["P_a"] = {
            "adjustment_from_zone_h": {str(k): v for k, v in sa},
            "monotone_non_decreasing": all(b >= a - 1.0e-9
                                           for (_, a), (_, b) in zip(sa, sa[1:])),
            "any_standoff_shorter_than_zero": any(v < base_a - 1.0e-9
                                                  for _, v in sa[1:]),
        }
    else:
        verdict["P_a"] = {"note": "fewer than two children reached equilibrium; the "
                                  "adjustment length is not measurable here"}
    stations = sorted({k for r in rows for k in r["residual_at"]},
                      key=lambda t: float(t.rstrip("h")))
    pb: Dict[str, object] = {}
    for stn in stations:
        vals = [(r["standoff_cells"], r["residual_at"].get(stn), r["median_spread"])
                for r in rows if r["residual_at"].get(stn) is not None]
        if len(vals) < 2:
            continue
        zero = next((v for k, v, _ in vals if k == 0), None)
        spreads = [x for _, _, x in vals if x is not None]
        spread = np.nanmedian(spreads) if spreads else None
        better = [k for k, v, _ in vals[1:]
                  if zero is not None and spread is not None
                  and v < zero - spread]
        pb[stn] = {
            "residual": {str(k): v for k, v, _ in vals},
            "spread": None if spread is None else float(spread),
            "standoffs_beating_zero_by_more_than_the_spread": better,
        }
    verdict["P_b"] = pb
    sc = series("adj_u_from_row1_h")
    if len(sc) >= 2:
        base_c = sc[0][1]
        verdict["P_c"] = {
            "adjustment_from_first_row_h": {str(k): v for k, v in sc},
            "monotone_non_decreasing": all(b >= a - 1.0e-9
                                           for (_, a), (_, b) in zip(sc, sc[1:])),
            "any_standoff_shorter_than_zero": any(v < base_c - 1.0e-9
                                                  for _, v in sc[1:]),
        }
    else:
        verdict["P_c"] = {"note": "fewer than two children reached equilibrium"}
    refuted = (verdict["P_a"].get("any_standoff_shorter_than_zero") is True
               or verdict["P_c"].get("any_standoff_shorter_than_zero") is True
               or any(v.get("standoffs_beating_zero_by_more_than_the_spread")
                      for v in pb.values()))
    # Can the comparison tell the children apart at all?  Either two of them
    # reached equilibrium at measurably different fetches, or the residuals at
    # some fixed station spread further than the sampling noise.  If neither,
    # the answer is INCONCLUSIVE and must not be read as agreement.
    spread_out = []
    for stn, v in pb.items():
        vals = [x for x in v["residual"].values() if x is not None]
        if len(vals) >= 2 and v["spread"] is not None:
            spread_out.append((max(vals) - min(vals)) > v["spread"])
    discriminating = bool(spread_out and any(spread_out)) or len(sc) >= 2
    verdict["discriminating"] = discriminating
    verdict["design_9_4_standoff_claim"] = (
        "REFUTED" if refuted
        else ("SUPPORTED" if discriminating else "INCONCLUSIVE"))
    verdict["note"] = (
        "SUPPORTED means no standoff reached equilibrium sooner -- measured from "
        "the inner edge of the zone (P-a) or from the first building face (P-c, "
        "the sharp one) -- and none was closer to equilibrium at a fixed station "
        "by more than the sampling spread; and that the comparison could "
        "have shown otherwise.  INCONCLUSIVE means the children were not "
        "distinguishable at all: the residuals at every fixed station lay within "
        "the sampling spread of each other and fewer than two children reached "
        "equilibrium, so nothing was resolved and this is not agreement.")
    return {"rows": rows, "verdict": verdict}


def summarise_v4(exp: Experiment, per_child: Dict[str, Dict]) -> Dict:
    rows = []
    for key, m in per_child.items():
        c = m.get("comparison")
        if not c:
            rows.append({"key": key, "note": "no baseline comparison"})
            continue
        d, a, u = c["tke_difference"], c["criterion_a_prime"], c["umean_difference"]
        bands = {}
        for name, sp in c["spectra_child_over_baseline"].items():
            for b in ("band_16_64m", "band_8_16m"):
                bands[f"{name}:{b}"] = sp["bands"][b]["mean_of_ratios"]
        rows.append({
            "key": key,
            "parent_layout": m["configuration"]["parent_layout"],
            "child_layout": m["configuration"]["child_layout"],
            "dtke_above_pct": (None if d["above"]["mean_relative"] is None
                               else 100 * d["above"]["mean_relative"]),
            "spread_above_pct": (None if d["above"].get("mean_spread") is None
                                 else 100 * d["above"]["mean_spread"]),
            "sigma": d["above"].get("significance"),
            "sigma_median": d["above"].get("median_significance"),
            "dtke_canopy_pct": (None if d["canopy"]["mean_relative"] is None
                                else 100 * d["canopy"]["mean_relative"]),
            "umean_rms_over_ustar": u["rms_over_ustar"],
            "umean_floor_over_ustar": u["noise_floor_over_ustar"],
            "criterion_a_prime": a["max_interior_umean_error_over_ustar"],
            "criterion_a_prime_passes": a["passes"],
            "bands": bands,
        })
    verdict: Dict[str, object] = {}
    for r in rows:
        if r.get("dtke_above_pct") is None:
            continue
        sig = r.get("sigma_median")
        verdict[r["key"]] = {
            "tke_difference_pct": r["dtke_above_pct"],
            "significance_median": sig,
            "measurable": None if sig is None else bool(abs(sig) > 1.0),
            "criterion_a_prime_passes": r["criterion_a_prime_passes"],
            "reading": (
                "the mismatched parent layout leaves no interior signature above "
                "the sampling spread" if sig is not None and abs(sig) <= 1.0
                else "the mismatched parent layout leaves a measurable interior "
                     "signature; the sign and size are in tke_difference_pct"),
        }
    return {"rows": rows, "verdict": verdict}


def summary_table_v3(summary: Dict) -> str:
    head = (f"{'child':12s} {'standoff':>10s} {'row1 fetch':>11s} {'rows':>5s} "
            f"{'x_adj|zone':>11s} {'x_adj|row1':>11s} {'x_adj self':>11s} "
            f"{'uw x_adj':>9s} {'err row1':>9s} {'err last':>9s} "
            f"{'spread':>7s} {'IBL last':>9s}")
    lines = [head, "-" * len(head)]

    def f(x, fmt="{:.2f}"):
        return "  --  " if x is None else fmt.format(x)

    for r in summary["rows"]:
        lines.append(
            f"{r['key']:12s} {r['standoff_cells']:4d} c {r['standoff_m']:5.0f}m "
            f"{r['first_row_fetch_h']:10.2f}h {r['n_rows']:5d} "
            f"{f(r['adj_u_from_zone_h']):>10s}h {f(r['adj_u_from_row1_h']):>10s}h "
            f"{f(r['adj_u_selfref_from_zone_h']):>10s}h "
            f"{f(r['adj_uw_from_zone_h']):>8s}h "
            f"{f(r['err_row1'], '{:.1%}'):>9s} {f(r['err_last'], '{:.1%}'):>9s} "
            f"{f(r['median_spread'], '{:.1%}'):>7s} {f(r['ibl_last_h']):>8s}h")
    stations = sorted({k for r in summary["rows"] for k in r["residual_at"]},
                      key=lambda s: float(s.rstrip("h")))
    if stations:
        lines += ["", "residual canopy-velocity error at fixed stations from the "
                      "inner zone edge (the cross-standoff comparison):",
                  f"{'child':12s} " + " ".join(f"{s:>9s}" for s in stations)]
        for r in summary["rows"]:
            lines.append(f"{r['key']:12s} " + " ".join(
                f(r["residual_at"].get(s), "{:.1%}").rjust(9) for s in stations))
    return "\n".join(lines)


def summary_table_v4(summary: Dict) -> str:
    head = (f"{'child':12s} {'parent':10s} {'child geom':10s} {'dTKE>2h':>9s} "
            f"{'spread':>8s} {'sigma':>7s} {'sig med':>8s} {'dTKE cnp':>9s} "
            f"{'<u> rms':>8s} {'floor':>8s} {'crit A2':>8s} {'':>5s}")
    lines = [head, "-" * len(head)]
    for r in summary["rows"]:
        if "dtke_above_pct" not in r:
            lines.append(f"{r['key']:12s} {r.get('note', '')}")
            continue

        def f(x, fmt="{:8.3f}"):
            return "     -- " if x is None else fmt.format(x)
        lines.append(
            f"{r['key']:12s} {r['parent_layout']:10s} {r['child_layout']:10s} "
            f"{f(r['dtke_above_pct'], '{:8.2f}')}% {f(r['spread_above_pct'], '{:7.2f}')}% "
            f"{f(r['sigma'], '{:7.1f}')} {f(r['sigma_median'], '{:8.1f}')} "
            f"{f(r['dtke_canopy_pct'], '{:8.2f}')}% "
            f"{f(r['umean_rms_over_ustar'], '{:8.4f}')} "
            f"{f(r['umean_floor_over_ustar'], '{:8.4f}')} "
            f"{f(r['criterion_a_prime'], '{:8.4f}')} "
            f"{'PASS' if r['criterion_a_prime_passes'] else 'FAIL':>5s}")
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #


def _disk_estimate(exp: Experiment) -> str:
    """Rough disk budget, in the units of :func:`config.Preset.summary`.

    Periodic runs can write ``&OUTPUT`` (``fielddump_interval``, whole domain),
    ``&NESTDUMP`` (``dtdump``, the child's band + 1 cell only) or both -- see
    ``config.Preset.parent_output``.  For a child, the nesting file's levels
    are counted at its own boundary ``cadence`` and its field dump at
    ``child_dtdump``; both default to ``dtdump`` but V3/V4's parents decouple
    them (fine ``&NESTDUMP`` cadence, coarse everything else).
    """
    lines, total = [], 0.0
    for r in exp.periodic:
        p = r.preset
        gb = 0.0
        if p.writes_fielddump:
            nt = p.production / p.fielddump_interval
            gb += nt * 3 * p.itot * p.jtot * p.ktot * 4 / 1e9
        if p.writes_nestdump:
            nb = p.nestdump_nzone()
            ni, nj = p.child_itot, p.child_jtot
            band_cells = ni * nj - max(ni - 2 * nb, 0) * max(nj - 2 * nb, 0)
            nt = p.production / p.dtdump
            gb += nt * 3 * band_cells * p.ktot * 4 / 1e9
        total += gb
        lines.append(f"    periodic {r.key:12s} dumps {gb:7.1f} GB ({p.parent_output})")
    for c in exp.children:
        if not c.default:
            continue
        p = c.preset
        nt_nest = p.production / p.cadence
        nt_dump = p.production / p.child_dtdump
        slab = 3 * 2 * (p.child_itot + p.child_jtot) * p.child_ktot * p.nzone * 8
        dump = 3 * p.child_itot * p.child_jtot * p.child_ktot * 4
        gb = (nt_nest * slab + nt_dump * dump) / 1e9
        total += gb
        lines.append(f"    child    {c.key:12s} nesting {nt_nest * slab / 1e9:6.1f} GB + "
                     f"dumps {nt_dump * dump / 1e9:6.1f} GB")
    lines.append(f"    {'total':21s} {total:7.1f} GB")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundir", type=Path)
    ap.add_argument("--experiment", default="v3")
    ap.add_argument("--only", default="", help="comma-separated child keys")
    ap.add_argument("--start-at", default="periodic-case", choices=STAGES)
    ap.add_argument("--stop-after", default="summary", choices=STAGES)
    ap.add_argument("--ibm-backend", default="auto")
    ap.add_argument("--baseline-dir", type=Path, default=None,
                    help="V4: the matched-geometry child run to compare against "
                         "(the V1 'converged' child, .../904)")
    ap.add_argument("--baseline-preset", default="converged",
                    help="V4: the preset describing --baseline-dir")
    ap.add_argument("--no-legacy-metrics", action="store_true",
                    help="skip the V1/V2 child-versus-parent block, which for "
                         "these experiments is context rather than a criterion")
    ap.add_argument("--prune-nesting", action="store_true")
    ap.add_argument("--yes", action="store_true")
    args = ap.parse_args()

    exp = get_experiment(args.experiment)
    nprocs = max([r.preset.nprocx * r.preset.nprocy for r in exp.periodic]
                 + [c.preset.child_nprocx * c.preset.child_nprocy
                    for c in exp.children])
    if nprocs > 8 and not (args.yes or os.environ.get("PBS_JOBID")
                           or os.environ.get("SLURM_JOB_ID")
                           or os.environ.get("UDALES_V1_CONFIRM")):
        print(f"experiment '{exp.name}' needs {nprocs} MPI ranks and hours of wall "
              f"time.\nSubmit tests/validation/nesting/submit_cx3_{exp.kind}.pbs, or "
              "pass --yes if you really mean to run it here.", file=sys.stderr)
        return 2

    rundir = args.rundir
    rundir.mkdir(parents=True, exist_ok=True)
    (rundir / "analysis").mkdir(exist_ok=True)
    only = [k for k in args.only.split(",") if k]
    children = ([exp.child(k) for k in only] if only else exp.default_children)
    first, last = STAGES.index(args.start_at), STAGES.index(args.stop_after)

    def wanted(stage: str) -> bool:
        return first <= STAGES.index(stage) <= last

    print(exp.summary())
    print(f"\nrundir               {rundir}")
    print(f"children acted on    {', '.join(c.key for c in children)}")
    print("estimated new output:")
    print(_disk_estimate(exp))
    (rundir / "experiment_summary.txt").write_text(exp.summary() + "\n",
                                                   encoding="ascii")

    baseline: Optional[Tuple[Path, Preset]] = None
    if exp.kind == "v4":
        if exp.baseline_child_key is not None:
            b = exp.child(exp.baseline_child_key)
            baseline = (rundir / b.preset.child_expnr, b.preset)
            if b.key not in [c.key for c in children]:
                children = [b] + children
        elif args.baseline_dir is not None:
            baseline = (args.baseline_dir, get_preset(args.baseline_preset))
        elif exp.external_baseline is not None:
            raise SystemExit(
                f"experiment '{exp.name}' is measured as a difference from the V1 "
                f"'{exp.external_baseline.name}' child; pass --baseline-dir at that "
                "run's child case directory (.../904)")

    t_all = time.time()
    timings: Dict[str, Dict[str, float]] = {}
    run_periodic(rundir, exp, wanted, args.ibm_backend, timings)

    per_child: Dict[str, Dict] = {}
    baseline_cache: Dict = {}
    for c in children:
        p = c.preset
        casedir = rundir / p.child_expnr
        parent_dir = rundir / exp.periodic_run(c.parent_key).preset.parent_expnr
        outdir = rundir / "analysis" / c.key
        t = timings.setdefault(f"child:{c.key}", {})
        print(f"\n=== child '{c.key}' ({p.name}, {p.child_expnr}): {c.note}")
        print(f"    case {casedir}, driven by {parent_dir}")

        if wanted("child-case"):
            t0 = time.time()
            # The parent writes whatever p.parent_output says (config.Preset,
            # inherited from the periodic run this child was built as a
            # replace() of): 'both' for V3/V4 production and tiny presets, so
            # the child is driven from the fine &NESTDUMP band rather than the
            # coarse &OUTPUT dump that exists only for the parent's own
            # periodic-stats.  See presets_geometry's V3_PARENT/V4_PARENT for
            # why the two cadences differ.
            source = "nestdump" if p.writes_nestdump else "fielddump"
            driving = make_child_case.DrivingParent.matched(parent_dir, p, source=source)
            try:
                make_child_case.build(parent_dir, rundir, p,
                                      ibm_backend=args.ibm_backend, driving=driving)
            finally:
                driving.dump.close()
            manifest = json.loads((casedir / "manifest.json").read_text())
            t["case"] = time.time() - t0
            print(f"    child case built in {t['case']:.1f} s: "
                  f"{manifest['n_parent_levels']} parent levels, |Phi|/A = "
                  f"{manifest['flux_residual_after_correction']['max_abs_normalised']:.2e}"
                  f", nest_lparentgeom = "
                  f"{'.false.' if manifest['building_free_zone'] else '.true.'}")
        if wanted("child"):
            t0 = time.time()
            run_solver(casedir, f"namoptions.{p.child_expnr}",
                       p.child_nprocx * p.child_nprocy, casedir / "child.log")
            t["run"] = time.time() - t0
            print(f"    child run finished in {t['run']:.1f} s")
            if args.prune_nesting:
                nest = casedir / f"nesting.inp.{p.child_expnr}.nc"
                if nest.exists():
                    gb = nest.stat().st_size / 1e9
                    nest.unlink()
                    print(f"    pruned {nest.name} ({gb:.1f} GB)")
        if wanted("analysis"):
            t0 = time.time()
            if exp.kind == "v3":
                m = analyse_child_v3(rundir, exp, c, casedir, parent_dir, outdir,
                                     not args.no_legacy_metrics)
                print(ag.summary_v3(m))
            else:
                is_baseline = (baseline is not None
                               and Path(baseline[0]).resolve() == casedir.resolve())
                m = analyse_child_v4(rundir, exp, c, casedir, parent_dir, outdir,
                                     not args.no_legacy_metrics,
                                     None if is_baseline else baseline,
                                     baseline_cache)
                print(ag.summary_v4(m))
            per_child[c.key] = m
            t["analysis"] = time.time() - t0
            print(f"    analysis finished in {t['analysis']:.1f} s")

    if wanted("summary") and per_child:
        summary = (summarise_v3(exp, per_child) if exp.kind == "v3"
                   else summarise_v4(exp, per_child))
        summary["timings_s"] = timings
        summary["experiment"] = exp.name
        table = (summary_table_v3(summary) if exp.kind == "v3"
                 else summary_table_v4(summary))
        out = rundir / "analysis"
        (out / f"{exp.kind}_summary.json").write_text(
            json.dumps(summary, indent=2, default=float) + "\n", encoding="ascii")
        (out / f"{exp.kind}_summary.md").write_text(
            f"# {exp.name}: {exp.headline}\n\n```\n{table}\n```\n\n"
            f"## verdict\n\n```json\n"
            f"{json.dumps(summary['verdict'], indent=2, default=float)}\n```\n",
            encoding="ascii")
        print()
        print(table)
        print()
        print(json.dumps(summary["verdict"], indent=2, default=float))
        print(f"\nresults in {out}")

    print(f"\n[run_geometry] total {time.time() - t_all:.1f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
