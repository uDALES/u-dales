#!/usr/bin/env python3
"""The reductions that are specific to refinement (design section 10.4 row V0).

``analyse.run`` already compares a child against the fine run it was cut from,
and for V0 that comparison needs no changes at all: every V0 child is on the
**same grid as the V1 child**, and the reference it is measured against is the
**same fine parent**, so profiles, resolved TKE, spectra, error-versus-distance
and criterion A are computed by exactly the code that produced section 10.5's
numbers.  That is deliberate -- it is what makes the V1 result the ``r = 1`` row
of V0's table rather than a separate experiment quoted alongside it.

What this module adds is the two things that only exist once ``r > 1``:

**Where the child's spectrum sits relative to the parent's filter scale.**
Section 10.5 located V1's deficit in the 8-64 m band.  At ``r = 2`` off a 4 m
parent the parent's Nyquist wavelength is 8 m and at ``r = 4`` off an 8 m parent
it is 16 m, so that band straddles the parent's cutoff at both ratios and the
single number quoted for it would average two physically different situations:
above the cutoff the child is *reproducing* structure the parent had, below it
the child must *generate* structure the parent never resolved.  Splitting the
ratio at the cutoff is the heart of the test, so the split is computed here from
the same spectra ``analyse`` already wrote, over the same wavenumbers -- the
child's grid and the reference's are identical, so the two spectra share their
``k`` axis bin for bin and no regridding is involved.

**Whether the divergence-preserving prolongation holds in the running solver.**
The writer's flux correction and the interpolation each have unit coverage
(P1-P17, U23-U39), but neither has ever carried a running refined case.  The
offline half is in the child's ``manifest.json`` -- the parent's own ``divmax``
next to the ``divmax`` of the prolonged field, which design section 1.3 says
must agree -- and the runtime half is in the child's log, which this module
parses: the boundary flux residual ``Phi``, the projection's ``divmax`` and
``divtot``, the zone misfit and the zone/interior pressure-gradient ratio.

**And one thing V0 needs that V1 did not: what the parent itself knew.**  A
child driven by a genuinely coarse parent inherits that parent's own biases, so
a mean-flow or TKE error measured against the fine truth is the parent's error
*plus* the nesting's.  ``make_child_case`` records the driving parent's own
interior profiles while it cuts the slabs, and :func:`parent_deficit` compares
them against the same fine truth.  Subtracting the ``filtered`` arm's value from
the ``coarse`` arm's at the same ratio separates "the filter threw the scales
away" from "the coarse LES got the remaining ones wrong".

Usage
-----
    python analyse_v0.py <point_analysis_dir> <child_case_dir>
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np


# --------------------------------------------------------------------------- #
# Spectra, split at the parent's filter scale
# --------------------------------------------------------------------------- #

#: Bands defined relative to the parent's Nyquist wavelength ``lambda_N = 2 dx_P``.
#: ``4 dx_P`` is the usual "actually resolved" limit of an LES -- between it and
#: the Nyquist a field is representable but badly damped by the discretisation --
#: so the middle band is kept separate rather than folded into either side.
NYQUIST_BANDS = (
    ("parent_resolved", 2.0, np.inf),      # lambda >= 4 dx_P
    ("parent_marginal", 1.0, 2.0),         # 2 dx_P <= lambda < 4 dx_P
    ("sub_parent_filter", 0.0, 1.0),       # lambda <  2 dx_P: the parent never had it
)


def nyquist_split(spectrum: Dict[str, object], parent_dx: float) -> Dict[str, object]:
    """Child/parent spectral ratio either side of the parent's Nyquist wavelength.

    ``spectrum`` is one entry of ``analyse``'s ``metrics["spectra"]``: the same
    ``k`` axis for both runs, ``E_parent`` from the fine reference and
    ``E_child`` from the nested child.  Bands are expressed in multiples of
    ``lambda_N = 2 dx_P`` so that ``r = 2`` and ``r = 4`` are read on the same
    axis even though their cutoffs differ by a factor of two.

    Both reductions of section 10.5 are kept: ``mean_of_ratios`` (what the
    design document quotes) and ``ratio_of_sums`` (energy conserving).
    """
    k = np.asarray(spectrum["wavenumber_rad_per_m"], dtype=float)
    ep = np.asarray(spectrum["E_parent"], dtype=float)
    ec = np.asarray(spectrum["E_child"], dtype=float)
    lam_n = 2.0 * float(parent_dx)
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = np.where(k > 0, 2.0 * np.pi / np.where(k > 0, k, 1.0), np.inf)
    out: Dict[str, object] = {
        "parent_dx_m": float(parent_dx),
        "parent_nyquist_wavelength_m": lam_n,
        "bands": {},
    }
    for name, lo, hi in NYQUIST_BANDS:
        sel = (k > 0) & (lam >= lo * lam_n) & (lam < hi * lam_n)
        n = int(sel.sum())
        entry: Dict[str, object] = {
            "n_modes": n,
            "wavelength_range_m": [lo * lam_n,
                                   None if np.isinf(hi) else hi * lam_n],
            "wavelength_range_over_nyquist": [lo, None if np.isinf(hi) else hi],
        }
        if n and np.any(ep[sel] > 0):
            ratio = ec[sel] / ep[sel]
            entry["mean_of_ratios"] = float(np.nanmean(ratio))
            entry["ratio_of_sums"] = float(ec[sel].sum() / ep[sel].sum())
            entry["E_parent_sum"] = float(ep[sel].sum())
            entry["E_child_sum"] = float(ec[sel].sum())
        else:
            entry["mean_of_ratios"] = None
            entry["ratio_of_sums"] = None
        out["bands"][name] = entry
    b = out["bands"]
    above = b["parent_resolved"]["ratio_of_sums"]
    below = b["sub_parent_filter"]["ratio_of_sums"]
    out["contrast"] = (None if above is None or below is None else below - above)
    out["note"] = (
        "'parent_resolved' is where the child is reproducing structure the parent "
        "had; 'sub_parent_filter' is where it must generate structure the parent "
        "never resolved.  'contrast' is the second minus the first: negative means "
        "the child is worse where it has to invent, positive means it is worse "
        "where it is being told."
    )
    return out


def spectra_split(metrics: Dict[str, object], parent_dx: float) -> Dict[str, object]:
    """:func:`nyquist_split` at every height ``analyse`` took a spectrum at."""
    return {name: nyquist_split(sp, parent_dx)
            for name, sp in metrics.get("spectra", {}).items()}


# --------------------------------------------------------------------------- #
# The runtime end-to-end check on the prolongation
# --------------------------------------------------------------------------- #

_LOG_PATTERNS = {
    "phi": r"Phi \(norm\) =\s*(\S+)",
    "phi_lid": r"Phi lid    =\s*(\S+)",
    "zone_misfit_rms": r"zone misfit rms \[m/s\] =\s*(\S+)",
    "gradp_ratio": r"\|grad p\| zone =\s*\S+\s+interior =\s*\S+\s+ratio =\s*(\S+)",
}


def _floats(text: str, pattern: str) -> List[float]:
    out = []
    for m in re.findall(pattern, text):
        try:
            out.append(float(m))
        except ValueError:
            pass
    return out


#: Design section 10.7 item 2 (plan section 7, R2(b)): the C1 pressure-response
#: diagnostic, ``|grad p|`` over the zone and over the trusted interior,
#: reported together with the ratio the same log line already carries -- the
#: number that says whether the linear reconstruction's extra divergence
#: source stays confined to the zone (ratio close to what the constant arm
#: shows) or leaks into the interior (ratio materially larger).
_GRADP_PATTERN = r"\|grad p\| zone =\s*(\S+)\s+interior =\s*(\S+)\s+ratio =\s*(\S+)"


def runtime_diagnostics(child_log: Path) -> Dict[str, object]:
    """Parse the solver's own nesting and divergence diagnostics from a run log.

    These are the end-to-end half of the prolongation check: the offline
    guarantee is that a solenoidal parent gives a solenoidal child target, and
    the guarantee is only worth what the running solver shows.  ``Phi`` is the
    normalised net volume flux through the boundary (``modnesting``), ``divmax``
    and ``divtot`` are what the projection left behind (``modpois``).

    Every reduction is reported as a **time-mean** (``mean_abs``/``mean``) in
    addition to the median/max already here: V0b's pressure-response question
    ("how far does it reach") is about the run's typical state, not its worst
    moment, so the summary table needs the mean of ``||Gp||`` zone, interior
    and their ratio -- plan section 7, R2(b).
    """
    text = Path(child_log).read_text(errors="replace")
    out: Dict[str, object] = {"log": str(child_log)}
    for key, pattern in _LOG_PATTERNS.items():
        vals = _floats(text, pattern)
        out[key] = {
            "n": len(vals),
            "mean_abs": (float(np.mean(np.abs(vals))) if vals else None),
            "max_abs": (float(np.max(np.abs(vals))) if vals else None),
            "median_abs": (float(np.median(np.abs(vals))) if vals else None),
            "last": (float(vals[-1]) if vals else None),
        }
    # The zone and interior magnitudes themselves, not only their ratio: the
    # generic loop above only captures the ratio's own capture group, and
    # "does the pressure response stay in the zone" needs both sides of it.
    gradp = re.findall(_GRADP_PATTERN, text)
    if gradp:
        gzone = np.asarray([float(a) for a, _, _ in gradp])
        gint = np.asarray([float(b) for _, b, _ in gradp])
        gratio = np.asarray([float(c) for _, _, c in gradp])
        out["gradp_zone"] = {"n": int(gzone.size), "mean": float(gzone.mean()),
                             "median": float(np.median(gzone)),
                             "max": float(gzone.max())}
        out["gradp_interior"] = {"n": int(gint.size), "mean": float(gint.mean()),
                                 "median": float(np.median(gint)),
                                 "max": float(gint.max())}
        out["gradp_ratio"]["mean"] = float(gratio.mean())
    else:
        out["gradp_zone"] = {"n": 0, "mean": None, "median": None, "max": None}
        out["gradp_interior"] = {"n": 0, "mean": None, "median": None, "max": None}
    div = re.findall(r"divmax, divtot =\s*(\S+)\s+(\S+)", text)
    dmax = [float(a) for a, _ in div]
    dtot = [float(b) for _, b in div]
    out["divmax"] = {"n": len(dmax),
                     "max": (float(np.max(dmax)) if dmax else None),
                     "median": (float(np.median(dmax)) if dmax else None)}
    out["divtot"] = {"n": len(dtot),
                     "max_abs": (float(np.max(np.abs(dtot))) if dtot else None),
                     "median_abs": (float(np.median(np.abs(dtot))) if dtot else None)}
    out["faces_forced"] = sorted(set(re.findall(r"face (\w+) is forced", text)))
    out["parent_io_warning"] = "parent I/O exceeds 1 % of runtime" in text
    return out


# --------------------------------------------------------------------------- #
# The staircase signature (plan section 7, R2; W8's original finding)
# --------------------------------------------------------------------------- #


def staircase_amplitude(profile: Dict[str, object], refine: int, ustar: float
                        ) -> Dict[str, object]:
    """Amplitude of the intra-parent-cell mean-flow error, at child resolution.

    W8 found that the piecewise-*constant* tangential reconstruction leaves a
    staircase in the child's mean wind: the target is identical across every
    child level inside one parent cell, so where the true profile is sheared
    the mismatch is a sawtooth of period ``refine`` child levels (period 2 at
    r = 2, period 4 at r = 4 -- design plan section 0, V0's "filtered" arm).
    The piecewise-*linear* reconstruction was written to remove exactly this,
    at the cost of a local divergence source (plan section 7, R2); reporting
    the amplitude for both arms is what makes the trade visible in the
    summary table rather than only in the mean-flow RMS, which mixes the
    staircase in with the smooth background error criterion A already scores.

    Isolated by removing, from the child-minus-truth mean-flow error at every
    level, the local mean of each group of ``refine`` consecutive child
    levels (one parent cell): what is left is exactly the intra-group
    (sawtooth) component, with the smooth background trend divided out group
    by group rather than assumed linear.  A matched-grid child (``refine ==
    1``) has one level per group and the residual is identically zero, as it
    must be -- there is no parent cell to be constant or linear across.
    """
    z = np.asarray(profile["z"], dtype=float)
    e = np.asarray(profile["u_child"], dtype=float) - np.asarray(profile["u_parent"], dtype=float)
    refine = int(refine)
    n = e.size
    ng = n // refine
    if ng < 1 or refine < 2:
        return {"available": False, "n_groups": 0, "refine": refine,
                "note": "refine < 2: no parent cell spans more than one child "
                        "level, so there is no staircase to measure"}
    e = e[: ng * refine].reshape(ng, refine)
    resid = e - e.mean(axis=1, keepdims=True)
    return {
        "available": True,
        "refine": refine,
        "n_groups": int(ng),
        "rms_over_ustar": float(np.sqrt(np.mean(resid ** 2)) / ustar),
        "max_abs_over_ustar": float(np.max(np.abs(resid)) / ustar),
        "note": ("RMS/max, over the whole profile, of the child-minus-truth "
                 "mean-flow error after subtracting each parent cell's own "
                 "group mean -- the intra-cell (sawtooth) component W8 "
                 "attributes to the tangential reconstruction, isolated from "
                 "the smooth background error criterion A already scores"),
    }


# --------------------------------------------------------------------------- #
# What the parent itself knew
# --------------------------------------------------------------------------- #


def parent_deficit(metrics: Dict[str, object], manifest: Dict[str, object],
                   building_height: float, ustar: float,
                   z_over_h_min: float = 2.0) -> Dict[str, object]:
    """The DRIVING parent's own interior statistics against the fine truth.

    Both are interior-only, fluid-only horizontal means; the fine truth comes
    from ``analyse``'s own parent accumulation and the driving parent's comes
    from ``make_child_case``, which accumulates it on the coarse grid while it
    is already reading every level.  The truth is interpolated onto the coarse
    levels rather than the other way round, so nothing is invented.

    This is the ceiling on what the child could possibly achieve.  A child
    cannot be more right than the boundary data it is given, and on the
    ``coarse`` arm the boundary data is wrong in two separate ways -- the filter
    threw scales away, and the coarse LES got the ones it kept slightly wrong.
    The ``filtered`` arm isolates the first, so the difference between the arms
    at one ratio is the second.
    """
    dp = manifest.get("driving_parent_profile")
    if not dp:
        return {"available": False,
                "note": "the child's manifest carries no driving-parent profile"}
    prof = metrics["profiles"]
    zf = np.asarray(prof["z"], dtype=float)
    zp = np.asarray(dp["z"], dtype=float)
    keep = (zp >= zf[0]) & (zp <= zf[-1])
    zp = zp[keep]
    # ``tke`` (and its dispersive/total siblings) is None for a band-only
    # (nestdump) driving source: no interior field and only one initial block
    # give no time series to take a temporal statistic over.  Compare the mean
    # flow regardless, mark the TKE comparison unavailable and say why, rather
    # than silently treating a missing number as zero or dropping the point.
    tke_available = dp.get("tke") is not None
    out: Dict[str, object] = {
        "available": True,
        "source": dp.get("source"),
        "parent_dx_m": dp.get("dx_m"),
        "z": zp.tolist(),
        "tke_available": tke_available,
        "note": ("the driving parent's own interior statistics against the fine "
                 "truth, on the driving grid; the ceiling on what the child "
                 "could reproduce.  'tke' is the per-cell temporal variance "
                 "(analyse.Bundle's definition, comparable to 'tke_parent' "
                 "below); the dispersive part (the spatial variance of the "
                 "time-mean field) is reported separately, not folded in"),
    }
    if not tke_available:
        out["tke_unavailable_reason"] = dp.get(
            "tke_unavailable_reason",
            "the driving source has no temporal statistic available")
    names = ["u"] + (["tke"] if tke_available else [])
    for name in names:
        key = f"{name}_parent"
        coarse = np.asarray(dp[name], dtype=float)[keep]
        fine = np.interp(zp, zf, np.asarray(prof[key], dtype=float))
        out[f"{name}_parent_coarse"] = coarse.tolist()
        out[f"{name}_truth_at_coarse_levels"] = fine.tolist()
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(np.abs(fine) > 0, (coarse - fine) / fine, np.nan)
        out[f"{name}_relative"] = rel.tolist()
        scale = ustar if name == "u" else ustar ** 2
        out[f"{name}_rms_difference"] = float(
            np.sqrt(np.nanmean((coarse - fine) ** 2)) / scale)
        sel = (zp / building_height) >= z_over_h_min
        out[f"{name}_mean_relative_above_{z_over_h_min:g}h"] = (
            float(np.nanmean(rel[sel])) if np.any(sel) else None)
        sel = (zp / building_height) < 1.0
        out[f"{name}_mean_relative_in_canopy"] = (
            float(np.nanmean(rel[sel])) if np.any(sel) else None)
    if tke_available:
        # The dispersive fraction is descriptive, not a deficit: there is no
        # truth-side dispersive number to compare it against (analyse.Bundle
        # does not compute one), only the driving parent's own split of its
        # own TKE into temporal and dispersive parts.  This is the number that
        # decides whether the (now consistent) TKE comparison above would have
        # looked materially different under the old, combined-mean estimator.
        disp = np.asarray(dp["tke_dispersive"], dtype=float)[keep]
        temp = np.asarray(dp["tke"], dtype=float)[keep]
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(np.abs(temp) > 0, disp / temp, np.nan)
        out["tke_dispersive_over_temporal"] = ratio.tolist()
        sel = (zp / building_height) >= z_over_h_min
        out[f"tke_dispersive_over_temporal_above_{z_over_h_min:g}h"] = (
            float(np.nanmean(ratio[sel])) if np.any(sel) else None)
        sel = (zp / building_height) < 1.0
        out["tke_dispersive_over_temporal_canopy"] = (
            float(np.nanmean(ratio[sel])) if np.any(sel) else None)
    return out


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #


def augment(outdir: Path, child_dir: Path, metrics: Dict[str, object],
            point, metrics_name: str = "v0_metrics.json") -> Dict[str, object]:
    """Add the V0 block to one point's metrics and rewrite the JSON in place."""
    outdir, child_dir = Path(outdir), Path(child_dir)
    manifest = json.loads((child_dir / "manifest.json").read_text())
    log = child_dir / "child.log"
    c = point.child
    v0: Dict[str, object] = {
        "point": point.key,
        "arm": point.arm,
        "refine": point.refine,
        "parent_dx_m": point.driver.dx,
        "parent_grid": [point.driver.itot, point.driver.jtot, point.driver.ktot],
        "parent_nyquist_wavelength_m": point.parent_nyquist_wavelength,
        "relaxation_ramp_resolved_by_parent": point.resolves_the_ramp,
        "L_rel_over_2dx_parent": c.zonewidth / (2.0 * point.driver.dx),
        "refinement": manifest.get("refinement"),
        "prolongation": manifest.get("prolongation"),
        "prolongation_requested": manifest.get("prolongation_requested"),
        "prolongation_offline": manifest.get("initial_condition_divmax"),
        "flux_residual_after_correction": manifest.get("flux_residual_after_correction"),
        "spectra_across_parent_nyquist": spectra_split(metrics, point.driver.dx),
        "parent_deficit": parent_deficit(metrics, manifest, c.building_height, c.ustar),
        "runtime": (runtime_diagnostics(log) if log.exists()
                    else {"log": str(log), "available": False}),
        "staircase": staircase_amplitude(metrics["profiles"], point.refine, c.ustar),
    }
    metrics["v0"] = v0
    (outdir / metrics_name).write_text(json.dumps(metrics, indent=2) + "\n",
                                       encoding="ascii")
    _write_csv(outdir, v0)
    return v0


def _write_csv(outdir: Path, v0: Dict[str, object]) -> None:
    rows = []
    for height, sp in v0["spectra_across_parent_nyquist"].items():
        for band, b in sp["bands"].items():
            rows.append((height, band, b["n_modes"],
                         b["wavelength_range_m"][0], b["wavelength_range_m"][1],
                         b["mean_of_ratios"], b["ratio_of_sums"]))
    if rows:
        with (Path(outdir) / "spectral_bands_across_parent_nyquist.csv").open(
                "w", encoding="ascii", newline="\n") as fh:
            fh.write("height,band,n_modes,lambda_min_m,lambda_max_m,"
                     "mean_of_ratios,ratio_of_sums\n")
            for r in rows:
                fh.write(",".join("" if x is None else
                                  (f"{x:.9g}" if isinstance(x, float) else str(x))
                                  for x in r) + "\n")
    pd = v0.get("parent_deficit", {})
    if pd.get("available"):
        cols = ["z", "u_parent_coarse", "u_truth_at_coarse_levels", "u_relative"]
        if pd.get("tke_available"):
            cols += ["tke_parent_coarse", "tke_truth_at_coarse_levels", "tke_relative",
                     "tke_dispersive_over_temporal"]
        with (Path(outdir) / "driving_parent_vs_truth.csv").open(
                "w", encoding="ascii", newline="\n") as fh:
            fh.write(",".join(cols) + "\n")
            for row in zip(*[pd[c] for c in cols]):
                fh.write(",".join("" if v is None or
                                  (isinstance(v, float) and np.isnan(v))
                                  else f"{v:.9g}" for v in row) + "\n")


def summary(v0: Dict[str, object], metrics: Dict[str, object]) -> str:
    """One point, in the shape section 10.5's V1 paragraph is written in."""
    d = metrics["v2"]["tke_deficit"]
    ca = metrics["v2"]["criterion_a"]
    lines = [
        f"V0 {v0['point']}: r = {v0['refine']} ({v0['arm']}), parent "
        f"{v0['parent_grid'][0]}x{v0['parent_grid'][1]}x{v0['parent_grid'][2]} @ "
        f"{v0['parent_dx_m']:g} m, lambda_Nyquist = "
        f"{v0['parent_nyquist_wavelength_m']:g} m",
    ]

    def pct(x):
        return "n/a" if x is None else f"{100 * x:+.2f}%"

    lines.append(f"  resolved-TKE deficit above z/h = 2:  "
                 f"{pct(d['above']['mean_relative'])}  "
                 f"(spread {pct(d['above']['median_spread'])} median)")
    lines.append(f"  inside the canopy (z/h < 1):         "
                 f"{pct(d['canopy']['mean_relative'])}")
    pdf = v0.get("parent_deficit", {})
    if pdf.get("available"):
        if pdf.get("tke_available"):
            lines.append(f"  the driving parent's own deficit:    "
                         f"{pct(pdf.get('tke_mean_relative_above_2h'))} TKE, "
                         f"{pdf.get('u_rms_difference', float('nan')):.4f} u* in the mean")
            lines.append(f"  dispersive/temporal TKE (parent):    "
                         f"z/h>2 {pct(pdf.get('tke_dispersive_over_temporal_above_2h'))}  "
                         f"canopy {pct(pdf.get('tke_dispersive_over_temporal_canopy'))}")
        else:
            lines.append(f"  the driving parent's own deficit:    "
                         f"TKE unavailable ({pdf.get('tke_unavailable_reason')}); "
                         f"{pdf.get('u_rms_difference', float('nan')):.4f} u* in the mean")
    for height, sp in v0["spectra_across_parent_nyquist"].items():
        b = sp["bands"]

        def r(name):
            v = b[name]["ratio_of_sums"]
            return "n/a" if v is None else f"{v:.3f}"

        lines.append(f"  spectra {height}: parent-resolved {r('parent_resolved')}  "
                     f"marginal {r('parent_marginal')}  "
                     f"sub-filter {r('sub_parent_filter')}")
    rt = v0.get("runtime", {})
    if rt.get("phi", {}).get("max_abs") is not None:
        lines.append(f"  runtime  max |Phi| = {rt['phi']['max_abs']:.2e}, "
                     f"max divmax = {rt['divmax']['max']:.2e}, "
                     f"zone misfit rms = {rt['zone_misfit_rms']['median_abs']:.3e} m/s, "
                     f"|grad p| zone/interior = {rt['gradp_ratio']['median_abs']:.2f}")
    if rt.get("gradp_zone", {}).get("mean") is not None:
        lines.append(f"  |grad p| (time-mean): zone = {rt['gradp_zone']['mean']:.3e}, "
                     f"interior = {rt['gradp_interior']['mean']:.3e}, "
                     f"ratio (of the means) = "
                     f"{rt['gradp_zone']['mean'] / rt['gradp_interior']['mean']:.2f}, "
                     f"mean of the per-report ratios = {rt['gradp_ratio']['mean']:.2f}")
    off = v0.get("prolongation_offline") or {}
    if off.get("parent_before_prolongation") is not None:
        lines.append(f"  prolongation  parent divmax "
                     f"{off['parent_before_prolongation']:.2e} -> child "
                     f"{off['before_projection']:.2e} -> projected "
                     f"{off['after_projection']:.2e}")
    if v0.get("prolongation") is not None:
        lines.append(f"  prolongation used: {v0['prolongation']!r} "
                     f"(requested {v0.get('prolongation_requested')!r})")
    st = v0.get("staircase") or {}
    if st.get("available"):
        lines.append(f"  staircase amplitude (intra-parent-cell mean-flow error): "
                     f"RMS {st['rms_over_ustar']:.4f} u*, "
                     f"max {st['max_abs_over_ustar']:.4f} u* "
                     f"over {st['n_groups']} parent-cell groups")
    lines.append(f"  criterion A (mean flow, interior): "
                 f"{ca['max_interior_umean_error_over_ustar']:.4f} u* against "
                 f"{ca['threshold']} -- {'PASS' if ca['passes'] else 'FAIL'}")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir", type=Path, help="the point's analysis directory")
    ap.add_argument("child_dir", type=Path)
    ap.add_argument("--suite", default="v0")
    ap.add_argument("--point", required=True)
    ap.add_argument("--metrics-name", default="v0_metrics.json")
    ns = ap.parse_args()

    from config import get_suite

    point = get_suite(ns.suite).point(ns.point)
    metrics = json.loads((ns.outdir / ns.metrics_name).read_text())
    v0 = augment(ns.outdir, ns.child_dir, metrics, point, ns.metrics_name)
    print(summary(v0, metrics))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
