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
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from config import Sweep, SweepPoint, get_sweep

#: How each arm of a sweep is read: its title and reading hint, the abscissa
#: (a function of the :class:`config.SweepPoint`), and the axis label.  A
#: sweep declares which arms it has (``Sweep.arms``); this says what they mean.
ARMS: Dict[str, Dict[str, object]] = {
    "zone": {
        "title": "P1 -- zone width at a fixed child size (V2a)",
        "hint": "read DOWN: if the deficit is fetch limited it should barely move",
        "short": "P1: zone width",
        "x": lambda p: p.preset.zonewidth / p.preset.dx,
        "xlabel": r"$N_{\rm rel}$ [cells]",
        "xscale": "linear",
    },
    "size": {
        "title": "P2 -- child size at a fixed zone (V2b)",
        "hint": "read DOWN: less fetch should mean a bigger deficit",
        "short": "P2: child size",
        "x": lambda p: p.preset.interior_extent_h,
        "xlabel": r"interior extent [$h$]",
        "xscale": "linear",
    },
    "cadence": {
        "title": "C0 -- boundary-data cadence at the reference child",
        "hint": ("read DOWN: if the cadence causes the deficit, the 8-16 m and "
                 "16-64 m ratios at z/h = 2 rise towards 1 as the cadence falls "
                 "(the 0.5 s point reaching the z/h = 1 value, about 0.97) and the "
                 "Catmull-Rom point moves 16-64 m a little but 8-16 m not at all; "
                 "if the scheme causes it, every row reads the same"),
        "short": "C0: cadence",
        "x": lambda p: p.preset.cadence,
        "xlabel": "boundary cadence [s]",
        "xscale": "log",
    },
}


def _sort_key(arm: str):
    x = ARMS[arm]["x"]
    # ties (the same cadence, linear and Catmull-Rom) go linear first
    return lambda p: (float(x(p)), p.preset.timeinterp)


#: Column key -> (header, format).  One place, so the CSV, the Markdown and the
#: plots cannot disagree about what a column means.
COLUMNS: Sequence[tuple] = (
    ("key", "point", "{}"),
    ("expnr", "nr", "{}"),
    ("cadence_s", "cadence [s]", "{:g}"),
    ("timeinterp", "interp", "{}"),
    ("n_rel_cells", "N_rel", "{}"),
    ("zone_cells", "zone", "{}"),
    ("child_cells", "child", "{}"),
    ("interior_cells", "int.cells", "{}"),
    ("interior_extent_h", "int./h", "{:.2f}"),
    ("zone_free", "zone clear", "{}"),
    ("cubes_cleared", "cubes cut", "{}"),
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


#: The bands tabulated per sampled height for the cadence arm: the two section
#: 10.5 used, in the order the cadence hypothesis predicts them to recover.
HEIGHT_BANDS: Sequence[Tuple[str, str]] = (("band_8_16m", "8-16 m"),
                                           ("band_16_64m", "16-64 m"))


def bands_by_height(metrics: Dict) -> Dict[str, Dict[str, object]]:
    """``mean_of_ratios`` of every band at every sampled height, by height name.

    The whole spectral table of one point, not just the z/h = 2 row: C0's
    prediction is height-ordered (the lost band is ``2 U cadence``, and ``U``
    grows with height), so the deliverable is the ratio at each height.
    """
    out: Dict[str, Dict[str, object]] = {}
    for name, sp in sorted(metrics.get("spectra", {}).items(),
                           key=lambda kv: kv[1]["z_over_h"]):
        out[name] = {"z_over_h": sp["z_over_h"], "z_m": sp.get("z_m")}
        for band, b in sp["bands"].items():
            out[name][band] = b.get("mean_of_ratios")
            out[name][band + "_n_modes"] = b.get("n_modes")
    return out


def row_from_metrics(key: str, arms: Sequence[str], expnr: str,
                     metrics: Dict, reused: bool) -> Dict[str, object]:
    v2 = metrics["v2"]
    cfg, d, f = v2["configuration"], v2["tke_deficit"], v2["tke_error_vs_fetch"]
    cb = v2["common_block"]
    a = v2["criterion_a"]
    interp = cfg.get("nest_timeinterp", metrics.get("nest_timeinterp"))
    return {
        "key": key,
        "arms": "+".join(arms),
        "expnr": expnr,
        "reused": reused,
        "preset": metrics["preset"],
        # ``.get``: metrics written before C0 carry no cadence block; the
        # parent's dump interval was the cadence then.
        "cadence_s": cfg.get("cadence_s"),
        "cadence_stride": cfg.get("cadence_stride"),
        "c_dump_u0": cfg.get("C_dump_at_u0"),
        "timeinterp": ({1: "linear", 2: "CR"}.get(interp, interp)),
        "bands_by_height": bands_by_height(metrics),
        "n_rel_cells": cfg["N_rel_cells"],
        "n_imp_cells": cfg["N_imp_cells"],
        "zone_cells": cfg["zone_cells"],
        "nzone": cfg["nzone"],
        "child_cells": f"{cfg['child_cells'][0]}x{cfg['child_cells'][1]}",
        "interior_cells": cfg["interior_cells"],
        "interior_extent_h": cfg["interior_extent_h"],
        "zone_free": "yes" if cfg["building_free_zone"] else "NO",
        "cubes_cleared": cfg.get("n_cubes_cleared_from_child_zone"),
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


def _arm_points(sweep: Sweep, arm: str, by_key: Dict[str, Dict[str, object]]
                ) -> List[SweepPoint]:
    pts = [p for p in sweep.arm(arm) if p.key in by_key]
    pts.sort(key=_sort_key(arm))
    return pts


def cadence_band_table(pts: Sequence[SweepPoint],
                       by_key: Dict[str, Dict[str, object]]) -> str:
    """The C0 deliverable: the band ratios at **every** sampled height per
    point, next to the profile deficit above z/h = 2.

    Columns are generated from whatever heights the analysis sampled, so the
    tiny and production sweeps produce the same shape of table.
    """
    heights: Dict[str, float] = {}
    for p in pts:
        for name, hb in by_key[p.key].get("bands_by_height", {}).items():
            heights.setdefault(name, float(hb["z_over_h"]))
    names = sorted(heights, key=heights.get)
    columns: List[tuple] = [
        ("key", "point", "{}"), ("expnr", "nr", "{}"),
        ("cadence_s", "cadence [s]", "{:g}"), ("timeinterp", "interp", "{}"),
        ("c_dump_u0", "C_dump(u0)", "{:.2f}"),
    ]
    flat_rows: List[Dict[str, object]] = []
    for name in names:
        for band, label in HEIGHT_BANDS:
            columns.append((f"{band}@{name}", f"{label} @ z/h={heights[name]:.2f}",
                            "{:.3f}"))
    columns += [("tke_deficit_pct", "dTKE z/h>2 [%]", "{:+.2f}"),
                ("tke_spread_pct", "spread [%]", "{:.2f}"),
                ("samples", "samples", "{}")]
    for p in pts:
        r = dict(by_key[p.key])
        for name, hb in r.get("bands_by_height", {}).items():
            for band, _ in HEIGHT_BANDS:
                r[f"{band}@{name}"] = hb.get(band)
        flat_rows.append(r)
    return markdown_table(flat_rows, columns)


def arm_tables(sweep: Sweep, rows: Sequence[Dict[str, object]]) -> str:
    """The deliverable: every arm of the sweep, with the reference in each."""
    by_key = {r["key"]: r for r in rows}
    out: List[str] = []
    for arm in sweep.arms:
        spec = ARMS[arm]
        pts = _arm_points(sweep, arm, by_key)
        out += [f"### {spec['title']}", "", f"*{spec['hint']}.*", "",
                markdown_table([by_key[p.key] for p in pts]), ""]
        if arm == "cadence":
            out += ["#### Band ratios (child/parent) at every sampled height", "",
                    "*The pre-registered prediction is on 8-16 m and 16-64 m at "
                    "z/h = 2; the lower heights say where the deficit starts.*", "",
                    cadence_band_table(pts, by_key), ""]
    missing = [p.key for p in sweep.points if p.key not in by_key]
    if missing:
        out += [f"**Incomplete:** no metrics for {', '.join(missing)}.", ""]
    return "\n".join(out)


def verdict(sweep: Sweep, rows: Sequence[Dict[str, object]]) -> Dict[str, object]:
    """State what each arm did, without deciding for the reader.

    The numbers are ranges and slopes, not a pass or a fail: these are
    falsification tests, and both outcomes are informative.  What is reported
    is how much the deficit -- and, for the cadence arm, the band ratios --
    moved along each arm, in units of the parent's own sampling spread, so that
    "barely moved" is a measurement rather than an impression.
    """
    by_key = {r["key"]: r for r in rows}

    def arm_stats(arm: str) -> Dict[str, object]:
        x_of = ARMS[arm]["x"]
        pts = _arm_points(sweep, arm, by_key)
        xs, ys, sp, spm, keys = [], [], [], [], []
        b8, b16 = [], []
        for p in pts:
            r = by_key[p.key]
            if r["tke_deficit_pct"] is None:
                continue
            keys.append(p.key)
            xs.append(float(x_of(p)))
            ys.append(float(r["tke_deficit_pct"]))
            sp.append(float(r["tke_spread_pct"] or np.nan))
            spm.append(float(r.get("tke_spread_median_pct") or np.nan))
            b8.append(r.get("band_8_16"))
            b16.append(r.get("band_16_64"))
        if len(ys) < 2:
            return {"n_points": len(ys)}
        spread = float(np.nanmean(sp))
        spread_med = float(np.nanmean(spm))

        def in_spreads(s: float) -> Optional[float]:
            if not np.isfinite(s) or s <= 0:
                return None
            return float((max(ys) - min(ys)) / s)

        out = {
            "n_points": len(ys),
            "keys": keys,
            "x": xs, "deficit_pct": ys,
            "band_8_16_at_z2h": b8,
            "band_16_64_at_z2h": b16,
            "range_pct": float(max(ys) - min(ys)),
            "mean_spread_pct": spread,
            "median_spread_pct": spread_med,
            "range_in_spreads": in_spreads(spread),
            "range_in_median_spreads": in_spreads(spread_med),
        }
        if len(set(xs)) >= 2:
            out["slope_pct_per_x"] = float(np.polyfit(xs, ys, 1)[0])
        return out

    notes = {
        "zone": ("the fetch interpretation of section 10.5 predicts a small "
                 "`range_in_spreads` for the zone arm"),
        "size": "and a large one for the size arm; the opposite ordering "
                "falsifies it and implicates the boundary treatment",
        "cadence": ("the cadence hypothesis (plan of 2026-09-06, section 0) "
                    "predicts `band_8_16_at_z2h` and `band_16_64_at_z2h` rising "
                    "monotonically as `x` (the cadence) falls, reaching about 0.97 "
                    "at 0.5 s; a plateau above 0.82 and below 0.97 is the scheme's "
                    "own deficit; no movement refutes the hypothesis"),
    }
    out: Dict[str, object] = {f"{arm}_arm": arm_stats(arm) for arm in sweep.arms}
    out["how_to_read"] = (
        "Compare `range_in_spreads` (and its more robust twin "
        "`range_in_median_spreads`, which uses the median per-height spread "
        "instead of the mean and so is not dragged up by the near-lid levels): "
        + "; ".join(notes[a] for a in sweep.arms if a in notes) + ".")
    out.update({
            "confounds": [
                r["key"] for r in rows if r.get("zone_free") == "NO"],
            "confound_note": (
                "'confounds' lists points that carry buildings inside their "
                "relaxation zone and so run with nest_lparentgeom = .true.; their "
                "deficits would carry a geometry change as well as the swept one. "
                "It should be empty: every point of this sweep clears its own "
                "zone."),
            "child_geometry_note": (
                "Points with a nonzero 'cubes cut' do not carry the parent's "
                "cubes inside their guard + ramp band -- the child clears them, "
                "the parent keeps them and goes on imprinting them through the "
                "imposed velocity field. So those children are not exact "
                "sub-models of the parent, but only inside the band, where the "
                "solution is imposed and no criterion is applied: "
                "Preset.removed_cubes_reaching_the_interior is asserted empty, so "
                "no cleared cube reaches the region the statistics are taken "
                "over, and analyse.run additionally excludes any cell solid in "
                "either run from both averages."),
            })
    return out


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

    # Flat CSV: the per-height band table is spread into one column per
    # (band, height) so nothing nested reaches the file.
    flat: List[Dict[str, object]] = []
    for r in rows:
        f = {k: v for k, v in r.items() if k != "bands_by_height"}
        for name, hb in (r.get("bands_by_height") or {}).items():
            for band, _ in HEIGHT_BANDS:
                f[f"{band}@{name}"] = hb.get(band)
            f[f"z_over_h@{name}"] = hb.get("z_over_h")
        flat.append(f)
    keys = sorted({k for r in flat for k in r})
    with (outdir / "sweep_summary.csv").open("w", encoding="ascii",
                                             newline="\n") as fh:
        fh.write(",".join(keys) + "\n")
        for r in flat:
            fh.write(",".join(
                "" if r.get(k) is None else
                (f"{r[k]:.9g}" if isinstance(r[k], float) else str(r[k]))
                for k in keys) + "\n")

    text = "\n".join([
        f"# Sweep '{sweep.name}' ({', '.join(ARMS[a]['short'] for a in sweep.arms)})"
        " -- the falsification table", "",
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
    arms = list(sweep.arms)
    n = len(arms)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4.2), constrained_layout=True,
                             squeeze=False)
    for a, arm in zip(axes[0], arms):
        spec = ARMS[arm]
        x_of = spec["x"]
        pts = [p for p in _arm_points(sweep, arm, by_key)
               if by_key[p.key]["tke_deficit_pct"] is not None]
        x = [x_of(p) for p in pts]
        y = [by_key[p.key]["tke_deficit_pct"] for p in pts]
        e = [by_key[p.key]["tke_spread_pct"] or 0.0 for p in pts]
        clear = [by_key[p.key]["zone_free"] == "yes" for p in pts]
        a.errorbar(x, y, yerr=e, fmt="-", color="k", lw=1.0, capsize=3, zorder=1)
        for xi, yi, ok, p in zip(x, y, clear, pts):
            cr = p.preset.timeinterp == 2
            a.plot([xi], [yi], ("o" if ok else "s") if not cr else "D",
                   color="tab:blue" if ok else "tab:red", zorder=2,
                   label=("Catmull-Rom in time" if cr else
                          "building-free zone" if ok else "buildings in the zone"))
        a.axhline(0.0, color="0.6", lw=0.8)
        a.set_xscale(spec["xscale"])
        a.set_xlabel(spec["xlabel"])
        a.set_ylabel("resolved-TKE deficit above $z/h=2$ [%]")
        a.set_title(spec["short"])
        a.grid(alpha=0.3)
        handles, labels = a.get_legend_handles_labels()
        seen: Dict[str, object] = {}
        for hh, ll in zip(handles, labels):
            seen.setdefault(ll, hh)
        if seen:
            a.legend(seen.values(), seen.keys(), fontsize=8)
    fig.suptitle(f"'{sweep.name}': what does the deficit track?")
    fig.savefig(outdir / "sweep_deficit.png", dpi=130)
    plt.close(fig)

    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4.2), constrained_layout=True,
                             squeeze=False)
    for a, arm in zip(axes[0], arms):
        spec = ARMS[arm]
        x_of = spec["x"]
        pts = _arm_points(sweep, arm, by_key)
        x = [x_of(p) for p in pts]
        for band, style in (("band_16_64", "o-"), ("band_8_16", "s--"),
                            ("band_large", "^:")):
            xy = [(xi, by_key[p.key][band]) for xi, p in zip(x, pts)
                  if by_key[p.key][band] is not None and p.preset.timeinterp == 1]
            if xy:
                a.plot([q[0] for q in xy], [q[1] for q in xy], style,
                       label=band.replace("band_", "").replace("_", "-") + " m")
            cr = [(xi, by_key[p.key][band]) for xi, p in zip(x, pts)
                  if by_key[p.key][band] is not None and p.preset.timeinterp == 2]
            if cr:
                a.plot([q[0] for q in cr], [q[1] for q in cr], "D", mfc="none",
                       color="k", label=f"{band.replace('band_', '').replace('_', '-')} m, CR")
        a.axhline(1.0, color="0.6", lw=0.8)
        a.set_xscale(spec["xscale"])
        a.set_xlabel(spec["xlabel"])
        a.set_ylabel(r"child/parent spectral ratio at $z/h \approx 2$")
        a.set_title(spec["short"])
        a.grid(alpha=0.3)
        a.legend(fontsize=8)
    fig.suptitle(f"'{sweep.name}': the scale-selective part of the deficit")
    fig.savefig(outdir / "sweep_spectra.png", dpi=130)
    plt.close(fig)

    if "cadence" in arms:
        # The C0 deliverable as a picture: every band at every height against
        # the cadence, one panel per band.
        pts = _arm_points(sweep, "cadence", by_key)
        fig, axes = plt.subplots(1, len(HEIGHT_BANDS), figsize=(5 * len(HEIGHT_BANDS), 4.2),
                                 constrained_layout=True, squeeze=False)
        for a, (band, label) in zip(axes[0], HEIGHT_BANDS):
            heights: Dict[str, float] = {}
            for p in pts:
                for name, hb in by_key[p.key].get("bands_by_height", {}).items():
                    heights.setdefault(name, float(hb["z_over_h"]))
            for name in sorted(heights, key=heights.get):
                xy = [(ARMS["cadence"]["x"](p), by_key[p.key]["bands_by_height"][name].get(band))
                      for p in pts if p.preset.timeinterp == 1
                      and name in by_key[p.key].get("bands_by_height", {})]
                xy = [q for q in xy if q[1] is not None]
                if xy:
                    a.plot([q[0] for q in xy], [q[1] for q in xy], "o-",
                           label=f"z/h = {heights[name]:.2f}")
                cr = [(ARMS["cadence"]["x"](p), by_key[p.key]["bands_by_height"][name].get(band))
                      for p in pts if p.preset.timeinterp == 2
                      and name in by_key[p.key].get("bands_by_height", {})]
                cr = [q for q in cr if q[1] is not None]
                if cr:
                    a.plot([q[0] for q in cr], [q[1] for q in cr], "D", mfc="none",
                           color="k", label=f"z/h = {heights[name]:.2f}, CR")
            a.axhline(1.0, color="0.6", lw=0.8)
            a.set_xscale("log")
            a.set_xlabel(ARMS["cadence"]["xlabel"])
            a.set_ylabel(f"child/parent ratio, {label}")
            a.set_title(label)
            a.grid(alpha=0.3)
            a.legend(fontsize=8)
        fig.suptitle(f"C0 '{sweep.name}': band ratios against the boundary cadence")
        fig.savefig(outdir / "sweep_cadence_bands.png", dpi=130)
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
