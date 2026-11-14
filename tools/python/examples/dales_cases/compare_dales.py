#!/usr/bin/env python3
"""Compare a uDALES thermo case (examples/301-303) with DALES output.

Usage:
    python compare_dales.py <case> <udales_run_dir> <dales_results_dir|none> <out_dir>
                            [--t0 H --t1 H] [--fit0 H --fit1 H]

case              301 (BOMEX), 302 (smoke) or 303 (DYCOMS-II RF01)
udales_run_dir    run directory holding thermostat_t.<case>.nc and thermostat_z.<case>.nc
dales_results_dir directory with tmser.001.nc and profiles.001.nc (e.g. a DALES
                  case's results_3.2.0), or the word none
out_dir           PNGs and summary.md are written here

Outputs
  timeseries.png  lwp_bar, cfrac, zi, ustar versus time (uDALES, DALES, published bands)
  profiles.png    mean thl, qt, ql, w'thl' (resolved), w'qt' (resolved) over the averaging window
  summary.md      mean values over the window, uDALES vs DALES, and a pass/fail table against
                  the published bands (see PUBLISHED); for 302/303 also the entrainment rate
                  we = dzi/dt + D zi

Published bands
  PUBLISHED maps (case, variable) -> (t0_h, t1_h, lo, hi, source). A band must be a
  value read from a paper or table and the dict carries the citation. It is empty:
  the papers (Siebesma et al. 2003, Bretherton et al. 1999, Stevens et al. 2005) were
  not available when this tool was written, and the PyCLES and DALES files hold no
  published ranges. Variables without a band are reported with "no cited band" and
  are neither passed nor failed. Add entries only with a citation.

Averaging windows
  WINDOWS holds the default windows [hours] used for the profiles and the mean values;
  they are analysis choices (Siebesma et al. 2003 average BOMEX over hours 3-6 and
  Stevens et al. 2005 average RF01 profiles over the last hour, both recalled, not
  re-checked here). Override with --t0/--t1. The entrainment rate is a straight-line fit
  of zi over FIT_WINDOWS (--fit0/--fit1).

Variable names
  uDALES thermostat_t: time, lwp_bar, cfrac, zi, zi_qt, zi_sv, ustar, ...
  uDALES thermostat_z: time, thl, qt, ql (zt); wthl, wqt (zm)
  DALES tmser.001.nc:  time, lwp_bar, cfrac, zi, ustar, we, ...
  DALES profiles.001.nc: time, zt, zm, thl, qt, ql (zt); wtlr (resolved w'thl'), wqtr (zm)
Definitions of zi differ: uDALES zi is the height of the maximum dthl/dz (zi_qt and
zi_sv are threshold heights); DALES tmser zi is its own diagnostic, so compare
trends rather than absolute values.
"""
from __future__ import annotations

import argparse
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
# {(case, variable): (t0_h, t1_h, lo, hi, source)}; variable is a time-series name
# (lwp_bar [kg/m2], cfrac, zi [m], ustar [m/s], we [m/s]). Only cited values.
PUBLISHED: dict[tuple[str, str], tuple[float, float, float, float, str]] = {}
WINDOWS = {"301": (3.0, 6.0), "302": (2.0, 3.0), "303": (3.0, 4.0)}      # hours
FIT_WINDOWS = {"302": (1.0, 3.0), "303": (1.0, 4.0)}                      # hours
ZI_VAR = {"302": ("zi_sv", 0.0), "303": ("zi_qt", 3.75e-6)}           # (variable, D [1/s])
TS_VARS = [("lwp_bar", "LWP [g/m$^2$]", 1e3), ("cfrac", "cloud cover [-]", 1.0),
           ("zi", "zi [m]", 1.0), ("ustar", "u* [m/s]", 1.0)]
PROF_VARS = [("thl", "thl", "zt", 1.0, "K"), ("qt", "qt", "zt", 1e3, "g/kg"),
             ("ql", "ql", "zt", 1e3, "g/kg"), ("wthl", "resolved w'thl'", "zm", 1.0, "K m/s"),
             ("wqt", "resolved w'qt'", "zm", 1e3, "g/kg m/s")]
DALES_PROF = {"thl": "thl", "qt": "qt", "ql": "ql", "wthl": "wtlr", "wqt": "wqtr"}


def read_nc(path: Path) -> dict[str, np.ndarray]:
    out = {}
    with netCDF4.Dataset(path) as d:
        for k, v in d.variables.items():
            out[k] = np.ma.filled(v[:].astype(float), np.nan) if v.ndim else np.array(float(v[:]))
    return out


def window_mean(t_h, y, t0, t1):
    m = (t_h >= t0) & (t_h <= t1) & np.isfinite(y)
    return float(np.mean(y[m])) if m.any() else float("nan")


def band_status(case, var, mean_val):
    b = PUBLISHED.get((case, var))
    if b is None:
        return "no cited band", None
    t0, t1, lo, hi, src = b
    return ("PASS" if lo <= mean_val <= hi else "FAIL"), b


def entrainment(t_s, zi, t0_h, t1_h, D):
    """we = dzi/dt + D zi over the fit window (straight-line fit of zi)."""
    m = (t_s >= t0_h * 3600) & (t_s <= t1_h * 3600) & np.isfinite(zi)
    if m.sum() < 3:
        return float("nan"), float("nan")
    slope, icpt = np.polyfit(t_s[m], zi[m], 1)
    zmean = float(np.mean(zi[m]))
    return float(slope + D * zmean), float(slope)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("case", choices=sorted(CASE_NAMES))
    ap.add_argument("udales_run_dir", type=Path)
    ap.add_argument("dales_results_dir")
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--t0", type=float, help="averaging window start [h]")
    ap.add_argument("--t1", type=float, help="averaging window end [h]")
    ap.add_argument("--fit0", type=float, help="entrainment fit start [h]")
    ap.add_argument("--fit1", type=float, help="entrainment fit end [h]")
    a = ap.parse_args(argv)

    case = a.case
    a.out_dir.mkdir(parents=True, exist_ok=True)
    ut = read_nc(a.udales_run_dir / f"thermostat_t.{case}.nc")
    uz = read_nc(a.udales_run_dir / f"thermostat_z.{case}.nc")
    have_dales = a.dales_results_dir.lower() != "none"
    if have_dales:
        dd = Path(a.dales_results_dir)
        dt = read_nc(dd / "tmser.001.nc")
        dp = read_nc(dd / "profiles.001.nc")
    t_u = ut["time"] / 3600.0
    w0, w1 = WINDOWS[case]
    w0 = a.t0 if a.t0 is not None else w0
    w1 = a.t1 if a.t1 is not None else w1
    if t_u.size == 0:
        sys.exit("thermostat_t has no records")
    tend = t_u.max()
    if tend < w1:
        print(f"note: run covers {tend:.2f} h, shorter than the window end {w1} h; "
              f"averaging over {w0}-{min(w1, tend):.2f} h")
        w1 = tend
        w0 = min(w0, w1)

    # ---------- time series ----------
    fig, axs = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    rows = []
    for ax, (var, label, scale) in zip(axs.flat, TS_VARS):
        ax.plot(t_u, ut[var] * scale, label="uDALES")
        mu = window_mean(t_u, ut[var] * scale, w0, w1)
        md = float("nan")
        if have_dales and var in dt:
            t_d = dt["time"] / 3600.0
            ax.plot(t_d, dt[var] * scale, label="DALES", lw=1)
            md = window_mean(t_d, dt[var] * scale, w0, w1)
        b = PUBLISHED.get((case, var))
        if b is not None:
            t0b, t1b, lo, hi, src = b
            ax.fill_between([t0b, t1b], lo * scale, hi * scale, color="0.8", label=f"published: {src}")
        ax.axvspan(w0, w1, color="C1", alpha=0.08)
        ax.set_ylabel(label)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
        status, _ = band_status(case, var, mu / scale)
        rows.append((var, label, mu, md, status))
    for ax in axs[1]:
        ax.set_xlabel("time [h]")
    fig.suptitle(f"{case} {CASE_NAMES[case]}: time series (shaded: averaging window {w0}-{w1} h)")
    fig.tight_layout()
    fig.savefig(a.out_dir / "timeseries.png", dpi=130)
    plt.close(fig)

    # ---------- profiles ----------
    fig, axs = plt.subplots(1, len(PROF_VARS), figsize=(3.2 * len(PROF_VARS), 5), sharey=True)
    tz = uz["time"] / 3600.0
    mz = (tz >= w0) & (tz <= w1)
    prof_rows = []
    for ax, (var, label, zname, scale, unit) in zip(axs, PROF_VARS):
        if var in uz and mz.any():
            prof = np.nanmean(uz[var][mz], axis=0) * scale
            ax.plot(prof, uz[zname], label="uDALES")
        if have_dales and DALES_PROF[var] in dp:
            tdp = dp["time"] / 3600.0
            md_ = (tdp > w0) & (tdp <= w1 + 1e-9)
            if md_.any():
                dprof = np.nanmean(dp[DALES_PROF[var]][md_], axis=0) * scale
                ax.plot(dprof, dp[zname], label="DALES", lw=1)
                if var in uz and mz.any() and dprof.size == prof.size:
                    prof_rows.append((label, unit, float(np.nanmax(np.abs(prof - dprof)))))
        ax.set_xlabel(f"{label} [{unit}]")
        ax.grid(alpha=0.3)
    axs[0].set_ylabel("z [m]")
    axs[0].legend(fontsize=8)
    fig.suptitle(f"{case} {CASE_NAMES[case]}: mean profiles, {w0}-{w1} h")
    fig.tight_layout()
    fig.savefig(a.out_dir / "profiles.png", dpi=130)
    plt.close(fig)

    # ---------- entrainment ----------
    ent = None
    if case in ZI_VAR:
        zv, D = ZI_VAR[case]
        f0, f1 = FIT_WINDOWS[case]
        f0 = a.fit0 if a.fit0 is not None else f0
        f1 = a.fit1 if a.fit1 is not None else f1
        we_u, dz_dt = entrainment(ut["time"], ut[zv], f0, f1, D)
        we_d = float("nan")
        if have_dales and "we" in dt:
            we_d = window_mean(dt["time"] / 3600.0, dt["we"], f0, f1)
        status, _ = band_status(case, "we", we_u)
        ent = (zv, D, f0, f1, we_u, dz_dt, we_d, status)

    # ---------- summary ----------
    with open(a.out_dir / "summary.md", "w") as f:
        f.write(f"# {case} {CASE_NAMES[case]}: uDALES vs DALES\n\n")
        f.write(f"uDALES run: `{a.udales_run_dir}`  \nDALES: `{a.dales_results_dir}`  \n")
        f.write(f"Averaging window: {w0}-{w1} h ({t_u.size} thermostat records, last at {tend:.3f} h)\n\n")
        f.write("| variable | uDALES mean | DALES mean | published band | status |\n|---|---|---|---|---|\n")
        for var, label, mu, md, status in rows:
            b = PUBLISHED.get((case, var))
            bs = "none cited" if b is None else f"{b[2]:g} to {b[3]:g} ({b[4]})"
            f.write(f"| {label.replace('$', '')} | {mu:.5g} | {md:.5g} | {bs} | {status} |\n")
        if ent is not None:
            zv, D, f0, f1, we_u, dz_dt, we_d, status = ent
            b = PUBLISHED.get((case, "we"))
            bs = "none cited" if b is None else f"{b[2]:g} to {b[3]:g} ({b[4]})"
            f.write(f"| entrainment rate we [m/s] ({zv}, D = {D:g}, fit {f0}-{f1} h) | {we_u:.4g} "
                    f"| {we_d:.4g} | {bs} | {status} |\n")
        f.write("\nNo band is shown unless a value could be cited; `no cited band` is neither a pass nor a fail.\n")
        if prof_rows:
            f.write("\n## Profile differences (max |uDALES - DALES| over height)\n\n| profile | max abs diff |\n|---|---|\n")
            for label, unit, dmax in prof_rows:
                f.write(f"| {label} | {dmax:.4g} {unit} |\n")
    print(f"wrote {a.out_dir}/timeseries.png, profiles.png, summary.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
