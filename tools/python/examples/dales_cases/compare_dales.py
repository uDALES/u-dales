#!/usr/bin/env python3
"""Compare uDALES thermo cases (examples/301-303) with DALES output, single runs or seed ensembles.

Usage:
    python compare_dales.py <case> <out_dir> --udales U1 [U2 ...] [--dales D1 [D2 ...]]
                            [--t0 H --t1 H] [--fit0 H --fit1 H] [--block S]
    python compare_dales.py <case> <udales_run_dir> <dales_results_dir|none> <out_dir> [options]  (one pair)

case        301 (BOMEX), 302 (smoke) or 303 (DYCOMS-II RF01)
--udales    uDALES run directories holding thermostat_t.<case>.nc and thermostat_z.<case>.nc
--dales     DALES run directories holding tmser.001.nc and profiles.001.nc
out_dir     PNGs, summary.md and summary.json are written here

Outputs
  timeseries.png  lwp_bar, cfrac, zi (profile definition), ustar versus time; thin lines are members
  profiles.png    ensemble-mean thl, qt, ql, w'thl', w'qt' (resolved) over the window, +- one spread
  summary.md      per-run window means, ensemble means, spreads and the ratio |mean_U - mean_D| / spread
                  for every scalar and profile; the published-band table (see PUBLISHED); 302/303
                  entrainment rate we = dzi/dt + D zi
  summary.json    the same numbers, machine readable

Time alignment
  DALES profiles.001.nc holds block means over `timeav` (600 s in the shipped DALES cases) labelled
  at the block END: the record at t covers (t - timeav, t]. uDALES thermostat_z holds slab-mean
  snapshots every tthermostat (60 s). Both codes' profiles are put on the same blocks: uDALES
  snapshots with t in (tb - block, tb] are averaged into block tb, where block is the DALES record
  spacing (or --block, default 600 s without DALES). A block is in the averaging window when
  tb - block >= t0 and tb <= t1, so a 3-6 h window uses the 18 blocks ending at 3h10 ... 6h for
  both codes (the old tool also took the DALES block ending at 3 h, i.e. 2h50-3h).
  Time series (60 s samples in both codes) use t0 < t <= t1.

Boundary-layer height
  The tmser/thermostat zi diagnostics differ between the codes (DALES: thv gradient on a 4-column
  composite unless iblh_* are set; uDALES: column-wise max dthl/dz), so both are recomputed here
  from the block-mean profiles with one definition:
    zi       half level of the maximum dthl/dz of the block-mean thl profile (as uDALES zh(kmx));
    zi_qt8   (303 only) highest downward crossing of qt = 8 g/kg in the block-mean qt profile,
             linear between full levels (the Stevens et al. 2005 inversion height).
  The window mean of zi is the mean over the window blocks.

Ensemble statistics
  For each quantity: window mean per run; ensemble mean; spread = sample standard deviation of the
  members' window means (ddof = 1; reported as 0 for a single run). ratio (min) = |mean_U - mean_D| /
  min(spread_U, spread_D) and ratio (max) = |mean_U - mean_D| / max(spread_U, spread_D). "Very close"
  (the Task 14 bar: the means differ by less than the spread within either ensemble) is ratio (min) < 1
  in the strict reading, ratio (max) < 1 in the lenient one; both are reported. For profiles the same
  is done level by level; the tables give the maximum |difference|, the rms difference, the rms spreads,
  both rms ratios and the fraction of levels with |diff| < min(spreads), over z <= zmax (PROFILE_ZMAX).

Published bands
  PUBLISHED maps (case, variable) -> (t0_h, t1_h, lo, hi, source). A band must be a value read from
  a paper or table and the dict carries the citation. It is empty (no band could be cited); variables
  without a band are reported "no cited band" and are neither passed nor failed.

Variable names
  uDALES thermostat_t: time, lwp_bar, cfrac, zi, zi_qt, ustar; thermostat_z: time, thl, qt, ql (zt), wthl, wqt (zm)
  DALES tmser.001.nc: time, lwp_bar, cfrac, ustar; profiles.001.nc: time, zt, zm, thl, qt, ql, wthlr, wqtr
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

try:
    import netCDF4
except ImportError:  # pragma: no cover
    sys.exit("netCDF4 is required (pip install netCDF4)")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CASE_NAMES = {"301": "BOMEX", "302": "smoke", "303": "DYCOMS-II RF01"}
PUBLISHED: dict[tuple[str, str], tuple[float, float, float, float, str]] = {}
WINDOWS = {"301": (3.0, 6.0), "302": (2.0, 3.0), "303": (3.0, 4.0)}      # hours
FIT_WINDOWS = {"302": (1.0, 3.0), "303": (1.0, 4.0)}                      # hours
DIVERGENCE = {"302": 0.0, "303": 3.75e-6}                                 # D [1/s] for we
QT_THRES = {"303": 8.0e-3}                                                # zi_qt8 threshold [kg/kg]
PROFILE_ZMAX = {"301": 2500.0, "302": 1250.0, "303": 1200.0}              # m, below the sponge
TS_VARS = [("lwp_bar", "LWP", "g/m2", 1e3), ("cfrac", "cloud cover", "-", 1.0),
           ("ustar", "u*", "m/s", 1.0)]
ZI_VARS = [("zi", "zi (max dthl/dz)", "m"), ("zi_qt8", "zi_qt8 (qt = 8 g/kg)", "m")]
PROF_VARS = [("thl", "thl", "zt", 1.0, "K"), ("qt", "qt", "zt", 1e3, "g/kg"),
             ("ql", "ql", "zt", 1e3, "g/kg"), ("wthl", "w'thl' (resolved)", "zm", 1.0, "K m/s"),
             ("wqt", "w'qt' (resolved)", "zm", 1e3, "g/kg m/s")]
DALES_PROF = {"thl": "thl", "qt": "qt", "ql": "ql", "wthl": "wthlr", "wqt": "wqtr"}


def read_nc(path: Path) -> dict[str, np.ndarray]:
    out = {}
    with netCDF4.Dataset(path) as d:
        for k, v in d.variables.items():
            out[k] = np.ma.filled(v[:].astype(float), np.nan) if v.ndim else np.array(float(v[:]))
    return out


# ---------------------------------------------------------------- profile diagnostics
def zi_grad(thl, zt, zm):
    """Half level of the maximum dthl/dz (uDALES thermostat definition, on a mean profile)."""
    g = np.diff(thl) / np.diff(zt)                    # gradient at face k (between k-1 and k)
    return float(zm[1 + int(np.nanargmax(g))])


def zi_qt(qt, zt, thres):
    """Highest downward crossing of `thres` (linear between full levels); nan if none."""
    for k in range(qt.size - 2, -1, -1):
        if qt[k] > thres and qt[k + 1] <= thres:
            return float(zt[k] + (thres - qt[k]) * (zt[k + 1] - zt[k]) / (qt[k + 1] - qt[k]))
    return float("nan")


class Run:
    """One run, reduced to 60 s time series and block-mean profiles."""

    def __init__(self, code, path, case, block=None):
        self.code, self.path, self.case = code, Path(path), case
        p = self.path
        if code == "uDALES":
            ts = read_nc(p / f"thermostat_t.{case}.nc")
            pr = read_nc(p / f"thermostat_z.{case}.nc")
            self.prof = {v: pr[v] for v, *_ in PROF_VARS}
        else:
            ts = read_nc(p / "tmser.001.nc")
            pr = read_nc(p / "profiles.001.nc")
            self.prof = {v: pr[DALES_PROF[v]] for v, *_ in PROF_VARS}
        self.zt, self.zm = pr["zt"], pr["zm"]
        if self.zt.ndim > 1:
            self.zt, self.zm = self.zt[0], self.zm[0]
        self.t_ts = ts["time"]
        self.ts = {v: ts[v] for v, *_ in TS_VARS if v in ts}
        self.t_prof = pr["time"]
        sp = np.diff(self.t_prof)
        self.native_block = float(np.median(sp)) if sp.size else 60.0
        self.block = block

    def blocks(self, block):
        """Block-mean profiles on blocks (tb - block, tb], tb multiples of block. Returns tb, {var: [nb, nz]}."""
        t = self.t_prof
        if self.code == "DALES" and abs(self.native_block - block) < 1e-6:
            return t.copy(), self.prof                      # already block means ending at t
        tb_all = np.ceil(np.round(t / block, 6)) * block   # block end that contains each sample
        tbs = np.unique(tb_all)
        out = {v: np.array([np.nanmean(a[tb_all == tb], axis=0) for tb in tbs]) for v, a in self.prof.items()}
        return tbs, out

    def zi_series(self, block):
        tb, pr = self.blocks(block)
        res = {"zi": np.array([zi_grad(pr["thl"][i], self.zt, self.zm) for i in range(tb.size)])}
        if self.case in QT_THRES:
            res["zi_qt8"] = np.array([zi_qt(pr["qt"][i], self.zt, QT_THRES[self.case]) for i in range(tb.size)])
        return tb, res


def ts_window(t_s, y, t0, t1):
    m = (t_s > t0 * 3600 + 1e-6) & (t_s <= t1 * 3600 + 1e-6) & np.isfinite(y)
    return float(np.mean(y[m])) if m.any() else float("nan")


def block_mask(tb, block, t0, t1):
    return (tb - block >= t0 * 3600 - 1e-6) & (tb <= t1 * 3600 + 1e-6)


def entrainment(tb, zi, f0, f1, D):
    m = (tb >= f0 * 3600 - 1e-6) & (tb <= f1 * 3600 + 1e-6) & np.isfinite(zi)
    if m.sum() < 3:
        return float("nan")
    slope = np.polyfit(tb[m], zi[m], 1)[0]
    return float(slope + D * np.mean(zi[m]))


def stats(vals):
    v = np.array([x for x in vals if np.isfinite(x)])
    if v.size == 0:
        return float("nan"), float("nan")
    return float(v.mean()), float(v.std(ddof=1)) if v.size > 1 else 0.0


def ratio(mu, md, su, sd):
    s = min(su, sd)
    if s > 0:
        return abs(mu - md) / s
    return 0.0 if abs(mu - md) < 1e-12 * max(abs(mu), 1.0) else float("inf")


# ---------------------------------------------------------------- main
def parse(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("case", choices=sorted(CASE_NAMES))
    ap.add_argument("paths", nargs="+", help="<out_dir> with --udales, or <udales_dir> <dales_dir|none> <out_dir>")
    ap.add_argument("--udales", nargs="+", default=None)
    ap.add_argument("--dales", nargs="*", default=None)
    ap.add_argument("--t0", type=float); ap.add_argument("--t1", type=float)
    ap.add_argument("--fit0", type=float); ap.add_argument("--fit1", type=float)
    ap.add_argument("--block", type=float, help="profile block length [s] (default: DALES record spacing, else 600)")
    ap.add_argument("--label", default="", help="title suffix")
    a = ap.parse_args(argv)
    if a.udales is None:
        if len(a.paths) != 3:
            ap.error("give <udales_dir> <dales_dir|none> <out_dir>, or <out_dir> --udales ...")
        a.udales = [a.paths[0]]
        a.dales = [] if a.paths[1].lower() == "none" else [a.paths[1]]
        a.out_dir = Path(a.paths[2])
    else:
        if len(a.paths) != 1:
            ap.error("with --udales give only <out_dir>")
        a.out_dir = Path(a.paths[0])
        a.dales = a.dales or []
    return a


def main(argv=None) -> int:
    a = parse(argv)
    case = a.case
    a.out_dir.mkdir(parents=True, exist_ok=True)
    runs = {"uDALES": [Run("uDALES", p, case) for p in a.udales],
            "DALES": [Run("DALES", p, case) for p in a.dales]}
    codes = [c for c in ("uDALES", "DALES") if runs[c]]
    block = a.block or (runs["DALES"][0].native_block if runs["DALES"] else 600.0)
    w0, w1 = WINDOWS[case]
    w0 = a.t0 if a.t0 is not None else w0
    w1 = a.t1 if a.t1 is not None else w1
    tend = min(r.t_ts.max() for c in codes for r in runs[c]) / 3600.0
    if tend < w1 - 1e-6:
        print(f"note: shortest run covers {tend:.3f} h; window end {w1} h -> {tend:.3f} h")
        w1 = tend
    f0, f1 = FIT_WINDOWS.get(case, (w0, w1))
    f0 = a.fit0 if a.fit0 is not None else f0
    f1 = min(a.fit1 if a.fit1 is not None else f1, tend)
    zmax = PROFILE_ZMAX[case]
    zi_names = [z for z in ZI_VARS if z[0] == "zi" or case in QT_THRES]

    # ---- per-run reductions
    res = {c: [] for c in codes}
    for c in codes:
        for r in runs[c]:
            d = {"path": str(r.path)}
            for v, *_ in TS_VARS:
                d[v] = ts_window(r.t_ts, r.ts[v], w0, w1) if v in r.ts else float("nan")
            tb, zis = r.zi_series(block)
            m = block_mask(tb, block, w0, w1)
            for z, *_ in zi_names:
                d[z] = float(np.nanmean(zis[z][m])) if m.any() else float("nan")
            if case in DIVERGENCE:
                zv = "zi_qt8" if case in QT_THRES else "zi"
                d["we"] = entrainment(tb, zis[zv], f0, f1, DIVERGENCE[case])
            tbp, pr = r.blocks(block)
            mp = block_mask(tbp, block, w0, w1)
            d["nblocks"] = int(mp.sum())
            d["_prof"] = {v: np.nanmean(pr[v][mp], axis=0) for v, *_ in PROF_VARS}
            d["_zi"] = (tb, zis)
            res[c].append(d)

    scalars = [(v, f"{lab} [{u}]", s) for v, lab, u, s in TS_VARS] + [(z, f"{lab} [{u}]", 1.0) for z, lab, u in zi_names]
    if case in DIVERGENCE:
        scalars.append(("we", f"we = dzi/dt + D zi [mm/s] (fit {f0}-{f1} h)", 1e3))

    summary = {"case": case, "window_h": [w0, w1], "block_s": block, "runs": {}, "scalars": {}, "profiles": {}}
    for c in codes:
        summary["runs"][c] = [{k: v for k, v in d.items() if not k.startswith("_")} for d in res[c]]
    for v, lab, s in scalars:
        e = {}
        for c in codes:
            mu, sd = stats([d[v] * s for d in res[c]])
            e[c] = {"mean": mu, "spread": sd, "members": [d[v] * s for d in res[c]]}
        if len(codes) == 2:
            e["diff"] = e["uDALES"]["mean"] - e["DALES"]["mean"]
            e["ratio"] = ratio(e["uDALES"]["mean"], e["DALES"]["mean"], e["uDALES"]["spread"], e["DALES"]["spread"])
            e["ratio_max"] = ratio(e["uDALES"]["mean"], e["DALES"]["mean"],
                                   max(e["uDALES"]["spread"], e["DALES"]["spread"]),
                                   max(e["uDALES"]["spread"], e["DALES"]["spread"]))
        summary["scalars"][v] = e

    zt0 = (runs[codes[0]][0].zt, runs[codes[0]][0].zm)
    for v, lab, zn, s, u in PROF_VARS:
        e = {}
        for c in codes:
            arr = np.array([d["_prof"][v] * s for d in res[c]])
            e[c] = (arr.mean(axis=0), arr.std(axis=0, ddof=1) if len(arr) > 1 else np.zeros(arr.shape[1]))
        z = zt0[0] if zn == "zt" else zt0[1]
        sel = z <= zmax
        if len(codes) == 2 and e["uDALES"][0].size == e["DALES"][0].size:
            dlt = e["uDALES"][0] - e["DALES"][0]
            su, sd_ = e["uDALES"][1], e["DALES"][1]
            ok = np.abs(dlt) < np.minimum(su, sd_)
            rms = lambda x: float(np.sqrt(np.nanmean(x[sel] ** 2)))  # noqa: E731
            k = int(np.nanargmax(np.abs(dlt[sel])))
            summary["profiles"][v] = {"label": lab, "unit": u, "max_abs_diff": float(np.abs(dlt[sel])[k]),
                                      "z_max_diff": float(z[sel][k]), "rms_diff": rms(dlt),
                                      "rms_spread_u": rms(su), "rms_spread_d": rms(sd_),
                                      "ratio_rms": ratio(rms(dlt), 0.0, rms(su), rms(sd_)),
                                      "ratio_rms_max": ratio(rms(dlt), 0.0, max(rms(su), rms(sd_)), max(rms(su), rms(sd_))),
                                      "frac_levels_within": float(np.mean(ok[sel]))}
        summary[f"_p_{v}"] = (z, e)

    # ---- plots
    colors = {"uDALES": "C0", "DALES": "C3"}
    nts = len(TS_VARS) + len(zi_names)
    fig, axs = plt.subplots(2, (nts + 1) // 2, figsize=(5 * ((nts + 1) // 2), 7), sharex=True)
    axs = axs.flat
    i = 0
    for v, lab, u, s in TS_VARS:
        ax = axs[i]; i += 1
        for c in codes:
            for j, r in enumerate(runs[c]):
                if v in r.ts:
                    ax.plot(r.t_ts / 3600, r.ts[v] * s, color=colors[c], lw=0.7, alpha=0.6, label=c if j == 0 else None)
        ax.set_ylabel(f"{lab} [{u}]")
    for z, lab, u in zi_names:
        ax = axs[i]; i += 1
        for c in codes:
            for j, d in enumerate(res[c]):
                tb, zis = d["_zi"]
                ax.plot(tb / 3600, zis[z], color=colors[c], lw=0.7, alpha=0.6, label=c if j == 0 else None)
        ax.set_ylabel(f"{lab} [{u}]")
    for ax in list(axs)[:i]:
        ax.axvspan(w0, w1, color="C1", alpha=0.08); ax.grid(alpha=0.3); ax.legend(fontsize=7)
        ax.set_xlabel("time [h]")
    for ax in list(axs)[i:]:
        ax.set_visible(False)
    fig.suptitle(f"{case} {CASE_NAMES[case]} {a.label}: members (window {w0}-{w1} h)")
    fig.tight_layout(); fig.savefig(a.out_dir / "timeseries.png", dpi=130); plt.close(fig)

    fig, axs = plt.subplots(1, len(PROF_VARS), figsize=(3.2 * len(PROF_VARS), 5), sharey=True)
    for ax, (v, lab, zn, s, u) in zip(axs, PROF_VARS):
        z, e = summary.pop(f"_p_{v}")
        for c in codes:
            m, sd_ = e[c]
            ax.plot(m, z, color=colors[c], label=f"{c} (n={len(runs[c])})")
            ax.fill_betweenx(z, m - sd_, m + sd_, color=colors[c], alpha=0.2)
        ax.set_xlabel(f"{lab} [{u}]"); ax.grid(alpha=0.3); ax.set_ylim(0, zmax)
    axs[0].set_ylabel("z [m]"); axs[0].legend(fontsize=8)
    fig.suptitle(f"{case} {CASE_NAMES[case]} {a.label}: ensemble-mean profiles +- spread, {w0}-{w1} h")
    fig.tight_layout(); fig.savefig(a.out_dir / "profiles.png", dpi=130); plt.close(fig)

    # ---- summary
    with open(a.out_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=1, default=float)
    with open(a.out_dir / "summary.md", "w") as f:
        f.write(f"# {case} {CASE_NAMES[case]}: uDALES vs DALES {a.label}\n\n")
        f.write(f"Window {w0}-{w1} h; profile blocks {block:g} s (window = blocks ending in ({w0} h + block, {w1} h]).\n\n")
        for c in codes:
            f.write(f"{c} runs: " + ", ".join(f"`{r.path}`" for r in runs[c]) + "\n\n")
        hdr = "| quantity | " + " | ".join(f"{c} members | {c} mean | {c} spread" for c in codes)
        hdr += " | diff U-D | ratio (min) | ratio (max) |" if len(codes) == 2 else " |"
        f.write(hdr + "\n|" + "---|" * (hdr.count("|") - 1) + "\n")
        for v, lab, s in scalars:
            e = summary["scalars"][v]
            row = f"| {lab} | " + " | ".join(
                ", ".join(f"{x:.4g}" for x in e[c]["members"]) + f" | {e[c]['mean']:.4g} | {e[c]['spread']:.3g}" for c in codes)
            row += f" | {e['diff']:+.3g} | {e['ratio']:.2f} | {e['ratio_max']:.2f} |" if len(codes) == 2 else " |"
            f.write(row + "\n")
        if summary["profiles"]:
            f.write(f"\n## Profiles (window means, z <= {zmax:g} m)\n\n")
            f.write("| profile | max abs diff (z) | rms diff | rms spread U | rms spread D | rms diff / min rms spread | rms diff / max rms spread | levels within bar |\n|---|---|---|---|---|---|---|---|\n")
            for v, p in summary["profiles"].items():
                f.write(f"| {p['label']} [{p['unit']}] | {p['max_abs_diff']:.3g} ({p['z_max_diff']:.0f} m) | {p['rms_diff']:.3g} | "
                        f"{p['rms_spread_u']:.3g} | {p['rms_spread_d']:.3g} | {p['ratio_rms']:.2f} | {p['ratio_rms_max']:.2f} | "
                        f"{100 * p['frac_levels_within']:.0f} % |\n")
        f.write("\nratio (min) = |mean_U - mean_D| / min(spread_U, spread_D), ratio (max) uses max(spread_U, spread_D); "
                "spread = sample std of the members' window means. ratio (min) < 1: the means differ by less than the spread "
                "of both ensembles (the strict reading of the Task 14 bar); ratio (max) < 1: less than the spread of at least one. "
                "Levels within bar: |diff(z)| < min(spread_U(z), spread_D(z)).\n")
        nh = int(np.floor(tend + 1e-6))
        if nh >= 1:
            f.write("\n## Hourly LWP [g/m2] (ensemble mean, spread)\n\n| code | " + " | ".join(f"{h}-{h + 1} h" for h in range(nh)) + " |\n|---|" + "---|" * nh + "\n")
            for c in codes:
                cells = []
                for h in range(nh):
                    mu, sd_ = stats([ts_window(r.t_ts, r.ts["lwp_bar"] * 1e3, h, h + 1) for r in runs[c]])
                    cells.append(f"{mu:.2f} ({sd_:.2f})")
                f.write(f"| {c} | " + " | ".join(cells) + " |\n")
        if PUBLISHED:
            f.write("\n## Published bands\n\n| variable | band | uDALES mean | status |\n|---|---|---|---|\n")
            for (cs, v), (t0b, t1b, lo, hi, src) in PUBLISHED.items():
                if cs == case and v in summary["scalars"]:
                    mu = summary["scalars"][v]["uDALES"]["mean"]
                    f.write(f"| {v} | {lo:g}-{hi:g} ({src}) | {mu:.4g} | {'PASS' if lo <= mu <= hi else 'FAIL'} |\n")
        else:
            f.write("\nNo published band is cited (PUBLISHED is empty): no pass/fail against the literature.\n")
    print(f"wrote {a.out_dir}/timeseries.png, profiles.png, summary.md, summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
