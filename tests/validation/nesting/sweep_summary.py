#!/usr/bin/env python3
"""Reduce a V2 sweep to one table -- the deliverable of the experiment.

Each sweep point has already been analysed by :mod:`analyse`, which wrote a
``v2_metrics.json`` next to its plots.  This module reads those, puts the zone
arm and the size arm side by side in a single table, and writes it as CSV, as
JSON and as a Markdown block, plus two summary plots.

The table exists to answer two opposed predictions, stated in section 10.5 of
``docs/udales-nesting-design.md`` and in ``README.md`` here:

    P1  if the V1 resolved-TKE deficit is *fetch* limited, widening or
        narrowing the relaxation zone should barely move it;
    P2  shrinking the child domain -- and so the fetch -- should move it a lot.

So the columns that matter are the deficit and the spectral band ratios, read
**down** each arm.  A deficit that tracks the zone width falsifies the fetch
interpretation and implicates the boundary treatment; that is a result, not a
failure, and the table is laid out so it would be impossible to miss.

Usage
-----
    python sweep_summary.py <rundir> [--sweep v2]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from config import Sweep, get_sweep

#: Column key -> (header, format).  One place, so the CSV, the Markdown and the
#: plots cannot disagree about what a column means.
COLUMNS: Sequence[tuple] = (
    ("key", "point", "{}"),
    ("expnr", "nr", "{}"),
    ("n_rel_cells", "N_rel", "{}"),
    ("zone_cells", "zone", "{}"),
    ("child_cells", "child", "{}"),
    ("interior_cells", "int.cells", "{}"),
    ("interior_extent_h", "int./h", "{:.2f}"),
    ("zone_free", "zone clear", "{}"),
    ("samples", "samples", "{}"),
    ("tke_deficit_pct", "dTKE z/h>2 [%]", "{:+.2f}"),
    ("tke_spread_pct", "spread [%]", "{:.2f}"),
    ("tke_significance", "sigma", "{:+.1f}"),
    ("tke_significance_median", "sigma med", "{:+.1f}"),
    ("tke_deficit_common_pct", "dTKE common [%]", "{:+.2f}"),
    ("tke_canopy_pct", "dTKE z/h<1 [%]", "{:+.2f}"),
    ("band_16_64", "E ratio 16-64 m", "{:.3f}"),
    ("band_8_16", "E ratio 8-16 m", "{:.3f}"),
    ("band_large", "E ratio > L/4", "{:.3f}"),
    ("tke_err_zone_edge", "err@zone edge", "{:.3f}"),
    ("tke_err_1h", "err@1h fetch", "{:.3f}"),
    ("tke_err_2h", "err@2h fetch", "{:.3f}"),
    ("tke_err_max_fetch", "err@max fetch", "{:.3f}"),
    ("max_fetch_h", "max fetch/h", "{:.2f}"),
    ("crossing_faces", "faces at floor", "{}"),
    ("criterion_a", "crit. A [u*]", "{:.4f}"),
    ("criterion_a_pass", "crit. A", "{}"),
)


def _spectral_band(metrics: Dict, band: str, z_over_h: float = 2.0,
                   scope: str = "interior") -> Optional[float]:
    """``mean_of_ratios`` in one band, at the sampled height nearest ``z_over_h``.

    ``scope='common'`` reads the spectra taken over the sweep's common central
    block instead of over each child's own interior, which is the comparison to
    quote across the size arm: the two spans differ by a factor of 2.6 there, so
    the windowing is not the same even though the physical bands are.
    """
    spectra = (metrics["spectra"] if scope == "interior"
               else metrics["v2"]["common_block"].get("spectra", {}))
    if not spectra:
        return None
    name = min(spectra, key=lambda n: abs(spectra[n]["z_over_h"] - z_over_h))
    return spectra[name]["bands"].get(band, {}).get("mean_of_ratios")


def row_from_metrics(key: str, arms: Sequence[str], expnr: str,
                     metrics: Dict, reused: bool) -> Dict[str, object]:
    v2 = metrics["v2"]
    cfg, d, f = v2["configuration"], v2["tke_deficit"], v2["tke_error_vs_fetch"]
    cb = v2["common_block"]
    a = v2["criterion_a"]
    return {
        "key": key,
        "arms": "+".join(arms),
        "expnr": expnr,
        "reused": reused,
        "preset": metrics["preset"],
        "n_rel_cells": cfg["N_rel_cells"],
        "n_imp_cells": cfg["N_imp_cells"],
        "zone_cells": cfg["zone_cells"],
        "nzone": cfg["nzone"],
        "child_cells": f"{cfg['child_cells'][0]}x{cfg['child_cells'][1]}",
        "interior_cells": cfg["interior_cells"],
        "interior_extent_h": cfg["interior_extent_h"],
        "zone_free": "yes" if cfg["building_free_zone"] else "NO",
        "nest_lparentgeom": cfg["nest_lparentgeom"],
        "zone_fraction": cfg["zone_fraction"],
        "zone_fraction_warns": cfg["zone_fraction_warns"],
        "n_cubes_in_zone": cfg["n_cubes_in_zone"],
        "samples": metrics["samples"]["child"],
        "tke_deficit_pct": _pct(d["above"]["mean_relative"]),
        # ``.get`` rather than ``[]`` on the aggregates that were added after the
        # first production analyses were written: a metrics file from an older
        # analyse.py should degrade to a blank column, not crash the table.
        "tke_spread_pct": _pct(d["above"]["mean_spread"]),
        "tke_spread_median_pct": _pct(d["above"].get("median_spread")),
        "tke_significance": d["above"]["significance"],
        "tke_significance_median": d["above"].get("median_significance"),
        "tke_canopy_pct": _pct(d["canopy"]["mean_relative"]),
        "tke_deficit_common_pct": _pct(
            cb.get("tke_deficit", {}).get("above", {}).get("mean_relative")),
        "tke_deficit_common_spread_pct": _pct(
            cb.get("tke_deficit", {}).get("above", {}).get("mean_spread")),
        "band_16_64": _spectral_band(metrics, "band_16_64m"),
        "band_8_16": _spectral_band(metrics, "band_8_16m"),
        "band_large": _spectral_band(metrics, "lambda_gt_quarter_L"),
        "band_16_64_common": _spectral_band(metrics, "band_16_64m", scope="common"),
        "band_8_16_common": _spectral_band(metrics, "band_8_16m", scope="common"),
        "tke_err_zone_edge": f["mean_error_at_zone_edge"],
        "tke_err_1h": f["mean_error_at_fetch_h"].get("1h"),
        "tke_err_2h": f["mean_error_at_fetch_h"].get("2h"),
        "tke_err_max_fetch": f["mean_error_at_max_fetch"],
        "max_fetch_h": f["max_fetch_h"],
        "crossing_faces": f"{f['faces_crossing_the_floor']}/4",
        "mean_crossing_fetch_h": f["mean_crossing_fetch_h"],
        "criterion_a": a["max_interior_umean_error_over_ustar"],
        "criterion_a_pass": "PASS" if a["passes"] else "FAIL",
        "u_rms_over_ustar": metrics["profile_metrics"]["u_rms_difference_over_ustar"],
        "tke_rms_over_ustar2":
            metrics["profile_metrics"]["tke_rms_difference_over_ustar2"],
        "phi_after_correction": (metrics.get("flux_residual_after_correction") or {})
            .get("max_abs_normalised"),
    }


def _pct(x: Optional[float]) -> Optional[float]:
    return None if x is None else 100.0 * x


def _fmt(value, spec: str) -> str:
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "--"
    try:
        return spec.format(value)
    except (TypeError, ValueError):
        return str(value)


def markdown_table(rows: Sequence[Dict[str, object]],
                   columns: Sequence[tuple] = COLUMNS) -> str:
    head = [c[1] for c in columns]
    body = [[_fmt(r.get(c[0]), c[2]) for c in columns] for r in rows]
    width = [max(len(head[i]), *(len(b[i]) for b in body)) if body else len(head[i])
             for i in range(len(head))]
    out = ["| " + " | ".join(h.ljust(w) for h, w in zip(head, width)) + " |",
           "|" + "|".join("-" * (w + 2) for w in width) + "|"]
    for b in body:
        out.append("| " + " | ".join(v.ljust(w) for v, w in zip(b, width)) + " |")
    return "\n".join(out)


def arm_tables(sweep: Sweep, rows: Sequence[Dict[str, object]]) -> str:
    """The deliverable: the two arms, side by side, with the reference in both."""
    by_key = {r["key"]: r for r in rows}
    out: List[str] = []
    titles = {
        "zone": ("P1 -- zone width at a fixed child size (V2a)",
                 "read DOWN: if the deficit is fetch limited it should barely move"),
        "size": ("P2 -- child size at a fixed zone (V2b)",
                 "read DOWN: less fetch should mean a bigger deficit"),
    }
    for arm in ("zone", "size"):
        pts = [p for p in sweep.arm(arm) if p.key in by_key]
        if arm == "zone":
            pts.sort(key=lambda p: p.preset.zonewidth)
        else:
            pts.sort(key=lambda p: p.preset.child_itot)
        title, hint = titles[arm]
        out += [f"### {title}", "", f"*{hint}.*", "",
                markdown_table([by_key[p.key] for p in pts]), ""]
    missing = [p.key for p in sweep.points if p.key not in by_key]
    if missing:
        out += [f"**Incomplete:** no metrics for {', '.join(missing)}.", ""]
    return "\n".join(out)


def verdict(sweep: Sweep, rows: Sequence[Dict[str, object]]) -> Dict[str, object]:
    """State what the two arms did, without deciding for the reader.

    The numbers are ranges and slopes, not a pass or a fail: V2 is a
    falsification test, and both outcomes are informative.  What is reported is
    how much the deficit moved along each arm, in units of the parent's own
    sampling spread, so that "barely moved" is a measurement rather than an
    impression.
    """
    by_key = {r["key"]: r for r in rows}

    def arm_stats(arm: str, x_of) -> Dict[str, object]:
        pts = [p for p in sweep.arm(arm) if p.key in by_key]
        xs, ys, sp, spm = [], [], [], []
        for p in pts:
            r = by_key[p.key]
            if r["tke_deficit_pct"] is None:
                continue
            xs.append(float(x_of(p)))
            ys.append(float(r["tke_deficit_pct"]))
            sp.append(float(r["tke_spread_pct"] or np.nan))
            spm.append(float(r.get("tke_spread_median_pct") or np.nan))
        if len(ys) < 2:
            return {"n_points": len(ys)}
        order = list(np.argsort(xs))
        xs = [xs[i] for i in order]
        ys = [ys[i] for i in order]
        sp = [sp[i] for i in order]
        spm = [spm[i] for i in order]
        spread = float(np.nanmean(sp))
        spread_med = float(np.nanmean(spm))

        def in_spreads(s: float) -> Optional[float]:
            if not np.isfinite(s) or s <= 0:
                return None
            return float((max(ys) - min(ys)) / s)

        out = {
            "n_points": len(ys),
            "x": xs, "deficit_pct": ys,
            "range_pct": float(max(ys) - min(ys)),
            "mean_spread_pct": spread,
            "median_spread_pct": spread_med,
            "range_in_spreads": in_spreads(spread),
            "range_in_median_spreads": in_spreads(spread_med),
        }
        out["slope_pct_per_x"] = float(np.polyfit(xs, ys, 1)[0])
        return out

    zone = arm_stats("zone", lambda p: p.preset.zonewidth / p.preset.dx)
    size = arm_stats("size", lambda p: p.preset.interior_extent_h)
    note = ("Both arms measured. Compare `range_in_spreads` (and its more "
            "robust twin `range_in_median_spreads`, which uses the median "
            "per-height spread instead of the mean and so is not dragged up by "
            "the near-lid levels): the fetch "
            "interpretation of section 10.5 predicts a small number for the zone "
            "arm and a large one for the size arm. The opposite ordering "
            "falsifies it and implicates the boundary treatment.")
    return {"zone_arm": zone, "size_arm": size, "how_to_read": note,
            "confounds": [
                r["key"] for r in rows if r.get("zone_free") == "NO"],
            "confound_note": (
                "Points listed under 'confounds' have buildings inside the "
                "relaxation zone and run with nest_lparentgeom = .true.. That is "
                "forced by reusing the V1 parent: the widest street in the cube "
                "array is 16 m and the zone needs 26 m, so the only child whose "
                "zone sits over open ground is the one the plaza was carved for. "
                "Their deficits carry a geometry change as well as a fetch "
                "change and must not be read as pure fetch."),
            }


def build(rundir: Path, sweep: Sweep, metrics_name: str = "v2_metrics.json"
          ) -> Dict[str, object]:
    rundir = Path(rundir)
    rows: List[Dict[str, object]] = []
    for pt in sweep.points:
        path = rundir / "analysis" / pt.key / metrics_name
        if not path.exists():
            print(f"[sweep] no metrics for point '{pt.key}' at {path}")
            continue
        metrics = json.loads(path.read_text())
        rows.append(row_from_metrics(pt.key, pt.arms, pt.expnr, metrics, pt.reuse))
    return {"sweep": sweep.name,
            "parent_preset": sweep.parent.name,
            "common_block_cells": sweep.common_block_cells,
            "rows": rows,
            "verdict": verdict(sweep, rows)}


def write(outdir: Path, sweep: Sweep, summary: Dict[str, object],
          make_plots: bool = True) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = summary["rows"]
    (outdir / "sweep_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="ascii")

    keys = sorted({k for r in rows for k in r})
    with (outdir / "sweep_summary.csv").open("w", encoding="ascii",
                                             newline="\n") as fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(
                "" if r.get(k) is None else
                (f"{r[k]:.9g}" if isinstance(r[k], float) else str(r[k]))
                for k in keys) + "\n")

    text = "\n".join([
        f"# V2 sweep '{sweep.name}' -- the falsification table", "",
        sweep.summary(), "",
        arm_tables(sweep, rows),
        "## How the arms moved", "",
        "```json",
        json.dumps(summary["verdict"], indent=2),
        "```", ""])
    (outdir / "sweep_summary.md").write_text(text, encoding="ascii")

    if make_plots and rows:
        try:
            _plots(outdir, sweep, rows)
        except Exception as exc:  # pragma: no cover - plotting is a convenience
            print(f"[sweep] plotting skipped: {exc}")


def _plots(outdir: Path, sweep: Sweep, rows: Sequence[Dict[str, object]]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    by_key = {r["key"]: r for r in rows}
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    for a, arm, x_of, xlabel in (
            (ax[0], "zone", lambda p: p.preset.zonewidth / p.preset.dx,
             r"$N_{\rm rel}$ [cells]"),
            (ax[1], "size", lambda p: p.preset.interior_extent_h,
             r"interior extent [$h$]")):
        pts = [p for p in sweep.arm(arm)
               if p.key in by_key and by_key[p.key]["tke_deficit_pct"] is not None]
        pts.sort(key=x_of)
        x = [x_of(p) for p in pts]
        y = [by_key[p.key]["tke_deficit_pct"] for p in pts]
        e = [by_key[p.key]["tke_spread_pct"] or 0.0 for p in pts]
        clear = [by_key[p.key]["zone_free"] == "yes" for p in pts]
        a.errorbar(x, y, yerr=e, fmt="-", color="k", lw=1.0, capsize=3, zorder=1)
        for xi, yi, ok in zip(x, y, clear):
            a.plot([xi], [yi], "o" if ok else "s",
                   color="tab:blue" if ok else "tab:red", zorder=2,
                   label=("building-free zone" if ok else
                          "buildings in the zone"))
        a.axhline(0.0, color="0.6", lw=0.8)
        a.set_xlabel(xlabel)
        a.set_ylabel("resolved-TKE deficit above $z/h=2$ [%]")
        a.set_title({"zone": "P1: zone width", "size": "P2: child size"}[arm])
        a.grid(alpha=0.3)
        handles, labels = a.get_legend_handles_labels()
        seen: Dict[str, object] = {}
        for hh, ll in zip(handles, labels):
            seen.setdefault(ll, hh)
        if seen:
            a.legend(seen.values(), seen.keys(), fontsize=8)
    fig.suptitle(f"V2 '{sweep.name}': does the deficit track the zone, "
                 "or the fetch?")
    fig.savefig(outdir / "sweep_deficit.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    for a, arm, x_of, xlabel in (
            (ax[0], "zone", lambda p: p.preset.zonewidth / p.preset.dx,
             r"$N_{\rm rel}$ [cells]"),
            (ax[1], "size", lambda p: p.preset.interior_extent_h,
             r"interior extent [$h$]")):
        pts = [p for p in sweep.arm(arm) if p.key in by_key]
        pts.sort(key=x_of)
        x = [x_of(p) for p in pts]
        for band, style in (("band_16_64", "o-"), ("band_8_16", "s--"),
                            ("band_large", "^:")):
            xy = [(xi, by_key[p.key][band]) for xi, p in zip(x, pts)
                  if by_key[p.key][band] is not None]
            if xy:
                a.plot([q[0] for q in xy], [q[1] for q in xy], style,
                       label=band.replace("band_", "").replace("_", "-") + " m")
        a.axhline(1.0, color="0.6", lw=0.8)
        a.set_xlabel(xlabel)
        a.set_ylabel(r"child/parent spectral ratio at $z/h \approx 2$")
        a.grid(alpha=0.3)
        a.legend(fontsize=8)
    fig.suptitle("V2: the scale-selective part of the deficit")
    fig.savefig(outdir / "sweep_spectra.png", dpi=130)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundir", type=Path)
    ap.add_argument("--sweep", default="v2")
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()
    sweep = get_sweep(args.sweep)
    summary = build(args.rundir, sweep)
    write(Path(args.rundir) / "analysis", sweep, summary,
          make_plots=not args.no_plots)
    print(arm_tables(sweep, summary["rows"]))
    print(json.dumps(summary["verdict"], indent=2))
    print(f"\nwritten to {Path(args.rundir) / 'analysis'}")


if __name__ == "__main__":
    main()
