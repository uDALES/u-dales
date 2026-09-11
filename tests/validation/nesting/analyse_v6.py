#!/usr/bin/env python3
"""V6 analysis: does mass drift over a long nested run? (design section 10.4).

``run_v6.py``'s child prints nothing to netCDF worth reading here -- the whole
point of the construction is that the boundary goes steady partway through and
stays there, so every diagnostic that matters is already in the solver's own
stdout: ``modchecksim.chkdiv``'s ``divmax, divtot`` line (throttled by
``tcheck``) and ``nesting.nesting_stats``'s block (throttled by
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
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# -- log line patterns -------------------------------------------------- #
# One line each, all printed by rank 0 only (modchecksim.f90, nesting_scheme.f90).
# Fortran prints a blown-up diagnostic as ``NaN``/``Infinity``, or -- when the
# value overflows an ES field -- as a run of asterisks.  Those tokens MUST be
# matched and carried through as non-finite floats rather than left unmatched.
# An unmatched line drops silently out of its series while the surrounding
# ``nesting: t =`` timestamps still parse, so a run that blew up keeps only
# the healthy samples that preceded it and then reads as bounded, trend-free
# and full-duration -- a blow-up certified as PASS.  ``_num`` converts.
_NONFIN = (r"(?:[-+]?(?:NaN|nan|NAN)"
           r"|[-+]?(?:Infinity|INFINITY|infinity|Inf|INF|inf)"
           r"|\*{2,})")
_FLOAT = r"(?:" + _NONFIN + r"|[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?)"
_ES_NUM = r"(?:" + _NONFIN + r"|[+-]?\d\.\d+E[+-]\d+)"


def _num(token: str) -> float:
    """Fortran numeric token -> float, mapping NaN/Inf/overflow onto floats."""
    t = token.strip()
    if not t or t.startswith("*"):
        return float("nan")
    low = t.lstrip("+-").lower()
    if low.startswith("nan"):
        return float("nan")
    if low in ("inf", "infinity"):
        return float("-inf") if t.startswith("-") else float("inf")
    return float(t.replace("D", "E").replace("d", "e"))

_P_CHECKSIM = re.compile(rf"Time of Simulation:\s*({_FLOAT})\s+dt:\s*({_FLOAT})")
# chkdiv's format is 2ES11.2 with NO literal separator between the two fields
# (modchecksim.f90: "write(6,'(A,2ES11.2)')'divmax, divtot = ', divmax, divtot"),
# so a negative value can butt straight up against the next field's sign with
# no whitespace at all ("1.23E-08-4.56E-07").  Matching each Fortran ES token
# by its own fixed shape (one leading digit, decimal point, exponent) rather
# than relying on \s+ to separate them handles that case too.
_P_DIVDIV = re.compile(rf"divmax, divtot =\s*({_ES_NUM})\s*({_ES_NUM})")
_P_NEST_T = re.compile(rf"nesting: t\s*=\s*({_FLOAT})")
_P_PHI = re.compile(rf"nesting: Phi \(norm\) =\s*({_FLOAT})")
_P_PHI_LID = re.compile(
    rf"nesting: Phi lid\s*=\s*({_FLOAT})\s*closed faces =\s*({_FLOAT})")
_P_MISFIT = re.compile(rf"nesting: zone misfit rms \[m/s\] =\s*({_FLOAT})")
_P_GRADP = re.compile(
    rf"nesting: \|grad p\| zone =\s*({_FLOAT})\s*interior =\s*({_FLOAT})\s*"
    rf"ratio =\s*({_FLOAT})")
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
    ``nesting:`` block into one record per ``t`` line relies on
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
            current_t, current_dt = _num(m.group(1)), _num(m.group(2))
            continue
        m = _P_DIVDIV.search(line)
        if m and current_t is not None:
            checksim_t.append(current_t)
            checksim_dt.append(current_dt if current_dt is not None else float("nan"))
            divmax.append(_num(m.group(1)))
            divtot.append(_num(m.group(2)))
            continue
        m = _P_NEST_T.search(line)
        if m:
            flush_nest()
            current_nest = {"t": _num(m.group(1))}
            continue
        if current_nest is not None:
            m = _P_PHI.search(line)
            if m:
                current_nest["phi"] = _num(m.group(1))
                continue
            m = _P_PHI_LID.search(line)
            if m:
                current_nest["phi_lid"] = _num(m.group(1))
                current_nest["phi_closed"] = _num(m.group(2))
                continue
            m = _P_MISFIT.search(line)
            if m:
                current_nest["misfit_rms"] = _num(m.group(1))
                continue
            m = _P_GRADP.search(line)
            if m:
                current_nest["gradp_zone"] = _num(m.group(1))
                current_nest["gradp_interior"] = _num(m.group(2))
                current_nest["gradp_ratio"] = _num(m.group(3))
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
    # A non-finite sample propagates through these aggregates; that is the
    # intended reading (they are not trustworthy), so silence the warning
    # rather than masking the values out.
    with np.errstate(invalid="ignore"):
        out["mean"] = float(yarr.mean())
        out["std"] = float(yarr.std())
        out["min"] = float(yarr.min())
        out["max"] = float(yarr.max())
        out["max_abs"] = float(np.max(np.abs(yarr)))
        out["range"] = float(yarr.max() - yarr.min())
    if x:
        out["t_first"] = float(x[0])
        out["t_last"] = float(x[-1])
        out["t_span"] = float(x[-1] - x[0])
    # A NaN/Inf diagnostic is the run blowing up, not evidence that is merely
    # incomplete: report it as a failure rather than letting it degrade into
    # INCONCLUSIVE (or, before the tokens were matched at all, vanish).
    nonfinite = ~np.isfinite(yarr)
    out["n_nonfinite"] = int(nonfinite.sum())
    if out["n_nonfinite"]:
        first_bad = next((xx for xx, bad in zip(x, nonfinite) if bad), None)
        out["verdict"] = "FAIL (non-finite)"
        out["reason"] = (
            f"{out['n_nonfinite']} of {n} samples are not finite"
            + ("" if first_bad is None else f", first at t = {first_bad:g}")
            + " -- the run produced NaN/Inf diagnostics")
        return out
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


def analyse(parsed: Dict[str, object], *, expected_runtime_s: Optional[float] = None,
            min_samples: int = 3, duration_tol: float = 0.95) -> Dict[str, object]:
    """Fit every series, then decide PASS/FAIL/NO-DATA/INCONCLUSIVE for the run.

    PASS is not merely "no FAIL was seen" -- that reads an empty or sparse log
    as a pass (an all-``/dev/null`` log has no FAIL verdicts because it has no
    verdicts of any useful kind). PASS instead requires positive evidence:

    * every series resolved to a real PASS/FAIL verdict (not NO DATA or
      INCONCLUSIVE), with finite statistics;
    * at least ``min_samples`` checksim and nesting_stats samples each;
    * the expected end-of-record freeze warning seen exactly once (never seen
      means the freeze was never exercised; more than once means something
      about the construction is off);
    * no abort past the record (``nest_lendabort = .false.`` is meant to
      freeze the boundary, not abort);
    * when ``expected_runtime_s`` is given (the caller's intended ``RUN.runtime``
      -- a tiny smoke run and the ~1e5-step production run pass different
      values here, so the acceptance bar scales with what was actually asked
      for rather than being weakened for the smoke test), the run actually
      covered at least ``duration_tol`` of it.

    Anything short of that is NO-DATA (nothing at all was parsed) or
    INCONCLUSIVE (something was parsed but the evidence is incomplete) --
    never PASS. An actual drift, or an unwanted abort, is FAIL. All of this is
    the aggregate's job, not each series': ``summarize_series`` above answers
    "is this one series drifting", not "does this run count as validated".
    """
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
    # Derived from sampled interval-mean timesteps (checksim_dt at the tcheck
    # cadence), NOT an exact step counter -- label it as approximate wherever
    # it is emitted (summary() does the same).
    steps_estimate = int(round(total_time / mean_dt)) if mean_dt and mean_dt > 0 else None

    n_checksim = len(checksim_t)
    n_nesting = len(nest_t)
    n_freeze = parsed["n_freeze_warnings"]
    aborted = parsed["aborted_past_record"]

    reasons: List[str] = []
    short_series: List[Tuple[str, float]] = []
    if expected_runtime_s is not None:
        for name, r in results.items():
            span = r.get("t_span")
            if span is None or float(span) < duration_tol * expected_runtime_s:
                short_series.append((name, 0.0 if span is None else float(span)))
    failing = [name for name, r in results.items()
              if str(r["verdict"]).startswith("FAIL")]
    blown_up = [name for name, r in results.items()
                if str(r["verdict"]).startswith("FAIL (non-finite")]
    drifting = [name for name in failing if name not in blown_up]
    if failing:
        overall = "FAIL"
        if blown_up:
            reasons.append(
                "non-finite diagnostics (the run blew up): " + ", ".join(sorted(blown_up)))
        if drifting:
            reasons.append(f"drifting series: {', '.join(sorted(drifting))}")
    elif aborted:
        overall = "FAIL"
        reasons.append(
            "run aborted past the parent record -- nest_lendabort = .false. "
            "should have frozen the boundary instead of aborting")
    elif n_checksim == 0 and n_nesting == 0:
        overall = "NO-DATA"
        reasons.append("no checksim or nesting_stats samples were parsed from the log")
    else:
        unresolved = [name for name, r in results.items()
                     if r["verdict"] in ("NO DATA", "INCONCLUSIVE")
                     or (r["n"] > 0 and not np.isfinite(r.get("max_abs", float("nan"))))]
        if unresolved:
            overall = "INCONCLUSIVE"
            reasons.append(
                f"series without a complete, finite verdict: {', '.join(sorted(unresolved))}")
        elif n_checksim < min_samples or n_nesting < min_samples:
            overall = "INCONCLUSIVE"
            reasons.append(
                f"too few samples (checksim={n_checksim}, nesting_stats={n_nesting}, "
                f"need >= {min_samples} of each)")
        elif n_freeze != 1:
            overall = "INCONCLUSIVE"
            reasons.append(
                "expected exactly one end-of-record freeze warning, log carries "
                f"{n_freeze}")
        elif expected_runtime_s is not None and total_time < duration_tol * expected_runtime_s:
            overall = "INCONCLUSIVE"
            reasons.append(
                f"run covered {total_time:.1f} s of the intended {expected_runtime_s:.1f} s "
                f"(< {100 * duration_tol:.0f}% reached)")
        elif expected_runtime_s is not None and short_series:
            # The record's timestamps and each diagnostic's own samples are
            # parsed independently, so a series can stop early -- or never
            # start -- while the run's apparent duration still looks complete.
            # Every series has to cover the interval it is claimed to certify.
            overall = "INCONCLUSIVE"
            reasons.append(
                "series not covering the intended interval: "
                + ", ".join(f"{nm} spans {sp:.1f} s" for nm, sp in sorted(short_series)))
        else:
            overall = "PASS"

    return {
        "series": results,
        "n_freeze_warnings": n_freeze,
        "aborted_past_record": aborted,
        "n_checksim_samples": n_checksim,
        "n_nesting_samples": n_nesting,
        "total_simulated_time_s": total_time,
        "expected_runtime_s": expected_runtime_s,
        "mean_dt_s": mean_dt,
        "steps_estimate": steps_estimate,
        "steps_estimate_is_approximate": True,
        "overall_verdict": overall,
        "overall_reasons": reasons,
    }


def run(child_dir: Path, outdir: Path, make_plots: bool = True, *,
        expected_runtime_s: Optional[float] = None, min_samples: int = 3,
        duration_tol: float = 0.95) -> Dict[str, object]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = Path(child_dir) / "child.log"
    parsed = parse_log(log_path)
    result = analyse(parsed, expected_runtime_s=expected_runtime_s,
                     min_samples=min_samples, duration_tol=duration_tol)
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
    expected = result.get("expected_runtime_s")
    duration_line = (
        f"Simulated time covered: {result['total_simulated_time_s']:.1f} s "
        f"(approx. {result['steps_estimate']} steps -- estimated from sampled "
        f"interval-mean dt = {result['mean_dt_s']:.4f} s, not an exact step counter)")
    if expected is not None:
        duration_line += f", intended runtime {expected:.1f} s"
    lines = [
        "# V6 -- does mass drift over long nested runs?",
        "",
        duration_line,
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
    for reason in result.get("overall_reasons", ()):
        lines.append(f"- {reason}")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("child_dir", type=Path)
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--expected-runtime-s", type=float, default=None,
                        help="the RUN.runtime [s] this child was actually asked to run "
                             "for -- PASS requires the log to show at least "
                             "--duration-tol of it was reached. Omit to skip that check "
                             "(e.g. exploring a partial or in-progress log).")
    parser.add_argument("--min-samples", type=int, default=3,
                        help="minimum checksim and nesting_stats sample count required "
                             "for PASS (default: 3)")
    parser.add_argument("--duration-tol", type=float, default=0.95,
                        help="fraction of --expected-runtime-s that must be covered for "
                             "PASS (default: 0.95)")
    args = parser.parse_args()
    result = run(args.child_dir, args.outdir, make_plots=not args.no_plots,
                expected_runtime_s=args.expected_runtime_s, min_samples=args.min_samples,
                duration_tol=args.duration_tol)
    print(summary(result))
    return 0 if result["overall_verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
