#!/usr/bin/env python3
"""V6 analysis: does mass drift over a long nested run? (design section 10.4).

``run_v6.py``'s child prints nothing to netCDF worth reading here -- the whole
point of the construction is that the boundary goes steady partway through and
stays there, so every diagnostic that matters is already in the solver's own
stdout: ``modchecksim.chkdiv``'s ``divmax, divtot`` line (throttled by
``tcheck``) and ``modnesting.nesting_stats``'s block (throttled by
``nest_statint`` -- ``Phi (norm)``, the zone misfit rms, ``|grad p|`` zone /
interior / ratio).  This module parses ``child.log``, fits a trend to each
series over the whole run, and answers the row's question directly: a slope
indistinguishable from zero (small compared to its own standard error) *and*
too small to explain a meaningful fraction of the series' own excursion is a
pass; anything else is a monotone drift and is reported as a fail, not
smoothed over.

Usage
-----
    python analyse_v6.py <child_dir> <outdir>
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

# -- log line patterns -------------------------------------------------- #
# One line each, all printed by rank 0 only (modchecksim.f90, modnesting.f90).
_P_CHECKSIM = re.compile(r"Time of Simulation:\s*([-\d.Ee+]+)\s+dt:\s*([-\d.Ee+]+)")
# chkdiv's format is 2ES11.2 with NO literal separator between the two fields
# (modchecksim.f90: "write(6,'(A,2ES11.2)')'divmax, divtot = ', divmax, divtot"),
# so a negative value can butt straight up against the next field's sign with
# no whitespace at all ("1.23E-08-4.56E-07").  Matching each Fortran ES token
# by its own fixed shape (one leading digit, decimal point, exponent) rather
# than relying on \s+ to separate them handles that case too.
_ES_NUM = r"[+-]?\d\.\d+E[+-]\d+"
_P_DIVDIV = re.compile(rf"divmax, divtot =\s*({_ES_NUM})\s*({_ES_NUM})")
_P_NEST_T = re.compile(r"modnesting: t\s*=\s*([-\d.Ee+]+)")
_P_PHI = re.compile(r"modnesting: Phi \(norm\) =\s*([-\d.Ee+]+)")
_P_PHI_LID = re.compile(
    r"modnesting: Phi lid\s*=\s*([-\d.Ee+]+)\s*closed faces =\s*([-\d.Ee+]+)")
_P_MISFIT = re.compile(r"modnesting: zone misfit rms \[m/s\] =\s*([-\d.Ee+]+)")
_P_GRADP = re.compile(
    r"modnesting: \|grad p\| zone =\s*([-\d.Ee+]+)\s*interior =\s*([-\d.Ee+]+)\s*"
    r"ratio =\s*([-\d.Ee+]+)")
_P_FREEZE_WARN = re.compile(
    r"WARNING.*boundary (?:will freeze|now freezes) on (?:the last level|that level)")
_P_ABORT = re.compile(r"the run extends past the end of the parent record")


def parse_log(path: Path) -> Dict[str, object]:
    """Pull every scalar time series out of one ``child.log``.

    Two independent clocks are in the file: ``checksim_t``/``checksim_dt`` (the
    ``tcheck`` cadence, carrying ``divmax``/``divtot``) and ``nest_records``
    (the ``nest_statint`` cadence, one dict per report carrying whichever of
    ``phi``/``phi_lid``/``phi_closed``/``misfit_rms``/``gradp_zone``/
    ``gradp_interior``/``gradp_ratio`` that report printed).  Grouping the
    ``modnesting:`` block into one record per ``t`` line relies on
    ``nesting_stats`` writing its lines consecutively with nothing else
    interleaved, which holds because rank 0 is the only writer.
    """
    text = Path(path).read_text(encoding="utf-8", errors="replace")

    checksim_t: List[float] = []
    checksim_dt: List[float] = []
    divmax: List[float] = []
    divtot: List[float] = []
    nest_records: List[Dict[str, float]] = []

    current_t: Optional[float] = None
    current_dt: Optional[float] = None
    current_nest: Optional[Dict[str, float]] = None

    def flush_nest() -> None:
        nonlocal current_nest
        if current_nest is not None:
            nest_records.append(current_nest)
        current_nest = None

    for line in text.splitlines():
        m = _P_CHECKSIM.search(line)
        if m:
            current_t, current_dt = float(m.group(1)), float(m.group(2))
            continue
        m = _P_DIVDIV.search(line)
        if m and current_t is not None:
            checksim_t.append(current_t)
            checksim_dt.append(current_dt if current_dt is not None else float("nan"))
            divmax.append(float(m.group(1)))
            divtot.append(float(m.group(2)))
            continue
        m = _P_NEST_T.search(line)
        if m:
            flush_nest()
            current_nest = {"t": float(m.group(1))}
            continue
        if current_nest is not None:
            m = _P_PHI.search(line)
            if m:
                current_nest["phi"] = float(m.group(1))
                continue
            m = _P_PHI_LID.search(line)
            if m:
                current_nest["phi_lid"] = float(m.group(1))
                current_nest["phi_closed"] = float(m.group(2))
                continue
            m = _P_MISFIT.search(line)
            if m:
                current_nest["misfit_rms"] = float(m.group(1))
                continue
            m = _P_GRADP.search(line)
            if m:
                current_nest["gradp_zone"] = float(m.group(1))
                current_nest["gradp_interior"] = float(m.group(2))
                current_nest["gradp_ratio"] = float(m.group(3))
                continue
    flush_nest()

    return {
        "checksim_t": checksim_t,
        "checksim_dt": checksim_dt,
        "divmax": divmax,
        "divtot": divtot,
        "nest_records": nest_records,
        "n_freeze_warnings": len(_P_FREEZE_WARN.findall(text)),
        "aborted_past_record": bool(_P_ABORT.search(text)),
    }


def _nest_series(records: Sequence[Dict[str, float]], key: str):
    t = [r["t"] for r in records if key in r]
    y = [r[key] for r in records if key in r]
    return t, y


# -- regression ----------------------------------------------------------- #

def linfit(x: Sequence[float], y: Sequence[float]) -> Dict[str, float]:
    """Ordinary least squares ``y = a + b x``, with the slope's standard error.

    ``se_b = sqrt(s^2 / Sxx)``, ``s^2`` the residual variance (``n - 2`` degrees
    of freedom).  Degenerate cases (fewer than 3 points, or every ``x`` equal)
    report ``nan`` for ``se_b`` rather than raising, since a probe run may not
    have enough samples for this to be meaningful -- callers check ``n``.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = x.size
    out = {"n": int(n), "slope": float("nan"), "slope_se": float("nan"),
           "intercept": float("nan")}
    if n < 2:
        return out
    b, a = np.polyfit(x, y, 1)
    out["slope"] = float(b)
    out["intercept"] = float(a)
    if n < 3:
        return out
    resid = y - (a + b * x)
    sxx = float(np.sum((x - x.mean()) ** 2))
    if sxx <= 0.0:
        return out
    s2 = float(np.sum(resid ** 2)) / (n - 2)
    out["slope_se"] = float(np.sqrt(s2 / sxx))
    return out


def summarize_series(name: str, x: Sequence[float], y: Sequence[float],
                     sig_t: float = 3.0, rel_tol: float = 0.2) -> Dict[str, object]:
    """Fit, then answer "is this drifting?" the way the row's question needs.

    A series **fails** (drifts) only when the trend is both statistically
    distinguishable from zero (``|slope| > sig_t * slope_se``) AND practically
    large (the trend's total change over the run is more than ``rel_tol`` of
    the series' own peak-to-peak excursion).  Either alone is not enough: a
    long run with round-off-level noise gives a "significant" slope of an
    utterly negligible size, and a short, noisy run can show a large swing that
    is not a trend at all -- ``t_stat`` and ``drift_over_range`` are both
    reported so the numbers, not just the verdict, survive into the record.
    """
    x = list(x)
    y = list(y)
    n = len(y)
    out: Dict[str, object] = {"name": name, "n": n}
    if n == 0:
        out.update(verdict="NO DATA", reason="no samples parsed for this series")
        return out
    yarr = np.asarray(y, dtype=float)
    out["mean"] = float(yarr.mean())
    out["std"] = float(yarr.std())
    out["min"] = float(yarr.min())
    out["max"] = float(yarr.max())
    out["max_abs"] = float(np.max(np.abs(yarr)))
    out["range"] = float(yarr.max() - yarr.min())
    fit = linfit(x, y)
    out.update(fit)
    if n < 3:
        out["t_stat"] = float("nan")
        out["verdict"] = "INCONCLUSIVE"
        out["reason"] = f"only {n} samples (need >= 3 with spread in t to fit a trend)"
        return out
    duration = x[-1] - x[0]
    drift_over_run = fit["slope"] * duration
    scale = out["range"] if out["range"] > 0.0 else max(out["max_abs"], 1.0e-300)
    relevant = abs(drift_over_run) > rel_tol * scale
    if not np.isfinite(fit["slope_se"]):
        out["t_stat"] = float("nan")
        out["verdict"] = "INCONCLUSIVE"
        out["reason"] = "slope standard error is not finite (degenerate fit)"
        return out
    if fit["slope_se"] == 0.0:
        # Every point sits exactly on the fitted line (typically because the
        # series is exactly constant, e.g. phi_lid on a case with no lid
        # faces): the t-statistic is 0/0 for a flat series or unboundedly
        # "significant" for a noiseless slope, neither of which is a useful
        # number, so the verdict is read off the trend's size alone.
        out["t_stat"] = 0.0 if fit["slope"] == 0.0 else float("inf")
        out["drift_over_run"] = float(drift_over_run)
        out["drift_over_range"] = float(abs(drift_over_run) / scale)
        if fit["slope"] == 0.0:
            out["verdict"] = "PASS (bounded, no drift)"
            out["reason"] = "series is exactly constant over the run"
        elif relevant:
            out["verdict"] = "FAIL (drifting)"
            out["reason"] = (f"exact (zero-residual) slope {fit['slope']:.3e} would move "
                             f"the series by {abs(drift_over_run):.3e} over the run, "
                             f"{100*out['drift_over_range']:.0f}% of its own range")
        else:
            out["verdict"] = "PASS (bounded, no drift)"
            out["reason"] = (f"exact (zero-residual) slope {fit['slope']:.3e} moves the "
                             f"series by only {100*out['drift_over_range']:.0f}% of its "
                             "own range over the run")
        return out
    t_stat = fit["slope"] / fit["slope_se"]
    out["t_stat"] = float(t_stat)
    out["drift_over_run"] = float(drift_over_run)
    out["drift_over_range"] = float(abs(drift_over_run) / scale)
    significant = abs(t_stat) > sig_t
    if significant and relevant:
        out["verdict"] = "FAIL (drifting)"
        out["reason"] = (f"slope {fit['slope']:.3e} is {abs(t_stat):.1f} standard errors "
                          f"from zero and would move the series by "
                          f"{abs(drift_over_run):.3e} over the run, "
                          f"{100*out['drift_over_range']:.0f}% of its own range")
    else:
        out["verdict"] = "PASS (bounded, no drift)"
        out["reason"] = (f"slope {fit['slope']:.3e} ({abs(t_stat):.1f} SE from zero) "
                          f"would move the series by {abs(drift_over_run):.3e} over the "
                          f"run, {100*out['drift_over_range']:.0f}% of its own range")
    return out


SERIES = ("divmax", "divtot", "phi", "phi_lid", "phi_closed", "misfit_rms",
          "gradp_zone", "gradp_interior", "gradp_ratio")


def analyse(parsed: Dict[str, object]) -> Dict[str, object]:
    results: Dict[str, object] = {}
    for key in ("divmax", "divtot"):
        results[key] = summarize_series(key, parsed["checksim_t"], parsed[key])
    for key in ("phi", "phi_lid", "phi_closed", "misfit_rms",
                "gradp_zone", "gradp_interior", "gradp_ratio"):
        t, y = _nest_series(parsed["nest_records"], key)
        results[key] = summarize_series(key, t, y)

    dt = [d for d in parsed["checksim_dt"] if np.isfinite(d) and d > 0]
    checksim_t = parsed["checksim_t"]
    nest_t = [r["t"] for r in parsed["nest_records"]]
    total_time = 0.0
    if checksim_t:
        total_time = max(total_time, checksim_t[-1] - checksim_t[0])
    if nest_t:
        total_time = max(total_time, nest_t[-1] - nest_t[0])
    mean_dt = float(np.mean(dt)) if dt else float("nan")
    steps_estimate = int(round(total_time / mean_dt)) if mean_dt and mean_dt > 0 else None

    return {
        "series": results,
        "n_freeze_warnings": parsed["n_freeze_warnings"],
        "aborted_past_record": parsed["aborted_past_record"],
        "n_checksim_samples": len(checksim_t),
        "n_nesting_samples": len(nest_t),
        "total_simulated_time_s": total_time,
        "mean_dt_s": mean_dt,
        "steps_estimate": steps_estimate,
        "overall_verdict": ("FAIL" if any(
            str(r["verdict"]).startswith("FAIL") for r in results.values())
            else "PASS"),
    }


def run(child_dir: Path, outdir: Path, make_plots: bool = True) -> Dict[str, object]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = Path(child_dir) / "child.log"
    parsed = parse_log(log_path)
    result = analyse(parsed)
    (outdir / "v6_summary.json").write_text(json.dumps(result, indent=2) + "\n",
                                            encoding="ascii")
    (outdir / "v6_summary.md").write_text(summary(result) + "\n", encoding="ascii")
    if make_plots:
        try:
            _plot(parsed, outdir)
        except ImportError:
            pass
    return result


def _plot(parsed: Dict[str, object], outdir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=False)
    ax = axes[0]
    if parsed["checksim_t"]:
        ax.semilogy(parsed["checksim_t"], np.abs(np.asarray(parsed["divmax"])) + 1e-300,
                   label="|divmax|")
        ax.semilogy(parsed["checksim_t"], np.abs(np.asarray(parsed["divtot"])) + 1e-300,
                   label="|divtot|")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("divergence")
    ax.set_title("chkdiv, throttled by tcheck")
    ax.legend()

    ax = axes[1]
    t, phi = _nest_series(parsed["nest_records"], "phi")
    t2, misfit = _nest_series(parsed["nest_records"], "misfit_rms")
    if t:
        ax.semilogy(t, np.abs(np.asarray(phi)) + 1e-300, label="|Phi (norm)|")
    if t2:
        ax.semilogy(t2, np.abs(np.asarray(misfit)) + 1e-300, label="zone misfit rms")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("nesting diagnostics")
    ax.set_title("nesting_stats, throttled by nest_statint")
    ax.legend()

    fig.tight_layout()
    fig.savefig(outdir / "v6_series.png", dpi=110)
    plt.close(fig)


def summary(result: Dict[str, object]) -> str:
    lines = [
        "# V6 -- does mass drift over long nested runs?",
        "",
        f"Simulated time covered: {result['total_simulated_time_s']:.1f} s "
        f"(~{result['steps_estimate']} steps at mean dt = {result['mean_dt_s']:.4f} s)",
        f"checksim samples: {result['n_checksim_samples']}, "
        f"nesting_stats samples: {result['n_nesting_samples']}",
        f"freeze warnings seen: {result['n_freeze_warnings']} "
        f"({'exactly once, as designed' if result['n_freeze_warnings'] == 1 else 'NOT once -- check the construction'})",
        f"aborted past the record: {result['aborted_past_record']} "
        "(must be False -- nest_lendabort = .false. is the whole point)",
        "",
        "| series | n | slope | slope SE | t | drift/range | verdict |",
        "|---|---|---|---|---|---|---|",
    ]
    for key in SERIES:
        r = result["series"][key]
        if r.get("verdict") in ("NO DATA", "INCONCLUSIVE"):
            lines.append(f"| {key} | {r['n']} | - | - | - | - | {r['verdict']}: "
                        f"{r.get('reason', '')} |")
            continue
        lines.append(
            f"| {key} | {r['n']} | {r['slope']:.3e} | {r['slope_se']:.3e} | "
            f"{r['t_stat']:.2f} | {100*r['drift_over_range']:.1f}% | {r['verdict']} |")
    lines.append("")
    lines.append(f"**Overall: {result['overall_verdict']}**")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("child_dir", type=Path)
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    result = run(args.child_dir, args.outdir, make_plots=not args.no_plots)
    print(summary(result))
    return 0 if result["overall_verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
