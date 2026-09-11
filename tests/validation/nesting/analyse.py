#!/usr/bin/env python3
"""Compare the nested child against the parent sub-region it was cut from.

This is the measurement half of design section 10.4 row V1.  It answers four
questions, in increasing order of how much they are worth:

1. do the mean profiles ``<u>(z)``, ``<v>(z)`` agree?
2. does the resolved TKE profile agree?
3. do the streamwise spectra agree, below, at and above building height?
4. **how far from the lateral boundary does the error decay to the noise
   floor?**  Design section 0 names this the acceptance criterion, and it is the
   single most valuable number the experiment produces.

Everything is measured over the child's **interior only** -- outside the guard
strip and the relaxation ramp.  Comparing inside the zone would be circular:
the solution there is imposed, so agreement is arithmetic, not physics.  The
zone values are still reported in the error-versus-distance curve (they are the
left-hand end of it), but no pass criterion is applied to them.

Everything is also measured against a **noise floor**.  The parent's own
statistics carry sampling error, so a nonzero child-minus-parent difference
does not by itself mean the scheme is wrong.  The floor is estimated by
splitting each run's own window in half and comparing the halves with itself;
a difference between child and parent that sits at that level is
indistinguishable from sampling noise, and the decay length is quoted against
it as well as against a fixed threshold.

Usage
-----
    python analyse.py <parent_case_dir> <child_case_dir> <outdir> [--preset ...]
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from caselib import FieldDump, cell_centred, load_solid_mask
from config import Preset, get_preset


# --------------------------------------------------------------------------- #
# Accumulation
# --------------------------------------------------------------------------- #


@dataclass
class Bundle:
    """Time-accumulated statistics of one run on the child's cell centres.

    Everything is stored as a *pair* of half-window accumulations (``a`` and
    ``b``), so the full-window mean is their average and the half-to-half
    difference gives the sampling-noise floor for free.
    """

    label: str
    times: np.ndarray
    #: [half][component] -> mean field, shape (ni, nj, nk)
    mean: Dict[str, List[np.ndarray]] = field(default_factory=dict)
    #: [half][component] -> mean of the square
    msq: Dict[str, List[np.ndarray]] = field(default_factory=dict)
    counts: Tuple[int, int] = (0, 0)
    #: height index -> (nt, ni, nj) time series of u at cell centres
    planes: Dict[int, np.ndarray] = field(default_factory=dict)
    #: interior-mean resolved TKE at each sampled time (equilibration trace)
    tke_series: np.ndarray = field(default_factory=lambda: np.zeros(0))
    #: name of the preset whose interior ``tke_series`` was averaged over.  Only
    #: interesting for a parent bundle shared between sweep points (see
    #: :func:`run`); everything else in a Bundle is window-independent.
    tke_series_interior_from: str = ""

    def full_mean(self, comp: str) -> np.ndarray:
        na, nb = self.counts
        return (na * self.mean[comp][0] + nb * self.mean[comp][1]) / (na + nb)

    def full_msq(self, comp: str) -> np.ndarray:
        na, nb = self.counts
        return (na * self.msq[comp][0] + nb * self.msq[comp][1]) / (na + nb)

    def variance(self, comp: str, half: Optional[int] = None) -> np.ndarray:
        if half is None:
            return np.maximum(self.full_msq(comp) - self.full_mean(comp) ** 2, 0.0)
        return np.maximum(self.msq[comp][half] - self.mean[comp][half] ** 2, 0.0)

    def tke(self, half: Optional[int] = None) -> np.ndarray:
        return 0.5 * sum(self.variance(c, half) for c in ("u", "v", "w"))


def accumulate(dump: FieldDump, extract: Callable, mask: np.ndarray,
               levels: Sequence[int], plane_k: Sequence[int],
               interior: Tuple[np.ndarray, np.ndarray], label: str) -> Bundle:
    """One streaming pass over ``levels``, building a :class:`Bundle`.

    ``extract`` maps a global ``(u, v, w)`` level onto the child window in the
    *raw* (upper-face-missing) convention, so that :func:`caselib.cell_centred`
    reduces the parent and the child by exactly the same rule.
    """
    levels = list(levels)
    n = len(levels)
    half = n // 2
    shape = None
    sums: List[Dict[str, np.ndarray]] = []
    sqs: List[Dict[str, np.ndarray]] = []
    counts = [0, 0]
    planes: Dict[int, List[np.ndarray]] = {k: [] for k in plane_k}
    tke_series = np.zeros(n)
    ii, jj = interior

    for idx, lev in enumerate(levels):
        u, v, w = extract(*dump.read_level(lev))
        uc, vc, wc = cell_centred(u, v, w)
        if shape is None:
            shape = uc.shape
            for _ in range(2):
                sums.append({c: np.zeros(shape) for c in "uvw"})
                sqs.append({c: np.zeros(shape) for c in "uvw"})
        h = 0 if idx < half else 1
        counts[h] += 1
        for c, arr in zip("uvw", (uc, vc, wc)):
            sums[h][c] += arr
            sqs[h][c] += arr * arr
        for k in plane_k:
            planes[k].append(uc[:, :, k].copy())
        sel = mask[np.ix_(ii, jj)]
        block = np.stack([uc[np.ix_(ii, jj)], vc[np.ix_(ii, jj)], wc[np.ix_(ii, jj)]])
        tke_series[idx] = 0.5 * float(np.sum(block ** 2 * sel[None]) / max(sel.sum(), 1))

    bundle = Bundle(label=label, times=np.asarray([dump.times[l] for l in levels]))
    bundle.counts = (counts[0], counts[1])
    for c in "uvw":
        bundle.mean[c] = [sums[h][c] / max(counts[h], 1) for h in range(2)]
        bundle.msq[c] = [sqs[h][c] / max(counts[h], 1) for h in range(2)]
    bundle.planes = {k: np.stack(v) for k, v in planes.items()}
    bundle.tke_series = tke_series
    return bundle


# --------------------------------------------------------------------------- #
# Geometry of the comparison
# --------------------------------------------------------------------------- #


def interior_indices(preset: Preset) -> Tuple[np.ndarray, np.ndarray]:
    """Cell indices (into the reduced arrays) outside the guard + ramp.

    Index ``i`` of a reduced array is child cell ``i`` (0-based), whose centre
    sits ``(i + 1/2) dx`` from the west face and ``(itot - i - 1/2) dx`` from the
    east face.  A cell is interior when both distances are at least
    ``L_imp + L_rel``.
    """
    lz = preset.guardwidth + preset.zonewidth
    i = np.arange(preset.child_itot - 1)
    j = np.arange(preset.child_jtot - 1)
    di_w, di_e = (i + 0.5) * preset.dx, (preset.child_itot - i - 0.5) * preset.dx
    dj_s, dj_n = (j + 0.5) * preset.dy, (preset.child_jtot - j - 0.5) * preset.dy
    return (np.where((di_w >= lz) & (di_e >= lz))[0],
            np.where((dj_s >= lz) & (dj_n >= lz))[0])


def central_indices(preset: Preset, ncells: int) -> Tuple[np.ndarray, np.ndarray]:
    """The central ``ncells x ncells`` block of the reduced arrays.

    V2 needs a region that is the **same** for every child in the sweep.  The
    per-point interiors are not: the interior shrinks faster than the domain
    (two zone widths come off whatever the size), so a 64-cell child keeps 40
    interior cells where a 128-cell one keeps 104, and comparing each over its
    own interior confounds "shorter fetch" with "smaller measurement window".
    Comparing all of them over one central block separates the two.

    The block is centred on the full child domain; the reduced arrays are one
    cell shorter (``cell_centred`` drops the upper face), so the block is
    off-centre by at most half a cell, which is immaterial next to the zone
    widths involved.
    """
    out = []
    for ntot in (preset.child_itot, preset.child_jtot):
        if ncells > ntot - 1:
            raise ValueError(f"a {ncells}-cell block does not fit in {ntot} cells")
        lo = (ntot - ncells) // 2
        out.append(np.arange(lo, lo + ncells))
    return out[0], out[1]


def masked_profile(field3d: np.ndarray, mask: np.ndarray,
                   ii: np.ndarray, jj: np.ndarray) -> np.ndarray:
    """Interior, fluid-only horizontal average, as a function of z."""
    sel = mask[np.ix_(ii, jj)]
    blk = field3d[np.ix_(ii, jj)]
    n = sel.sum(axis=(0, 1)).astype(float)
    out = np.where(n > 0, (blk * sel).sum(axis=(0, 1)) / np.maximum(n, 1), np.nan)
    return out


def slab_rms_difference(a: np.ndarray, b: np.ndarray, mask: np.ndarray,
                        axis: int, other: np.ndarray, kmax: Optional[int] = None
                        ) -> np.ndarray:
    """RMS of ``a - b`` over the fluid cells of each slab normal to ``axis``.

    ``other`` selects which indices of the *other* horizontal direction take
    part, so that the error profile along ``x`` is measured over the interior
    span in ``y`` (and vice versa) rather than over cells that are themselves
    inside a zone.
    """
    d = a - b
    m = mask.astype(float)
    if kmax is not None:
        d = d[:, :, :kmax]
        m = m[:, :, :kmax]
    if axis == 0:
        d, m = d[:, other, :], m[:, other, :]
    else:
        d, m = d[other, :, :], m[other, :, :]
        d, m = np.swapaxes(d, 0, 1), np.swapaxes(m, 0, 1)
    num = (d ** 2 * m).sum(axis=(1, 2))
    den = m.sum(axis=(1, 2))
    return np.sqrt(np.where(den > 0, num / np.maximum(den, 1), np.nan))


# --------------------------------------------------------------------------- #
# Spectra
# --------------------------------------------------------------------------- #


def streamwise_spectrum(planes: np.ndarray, ii: np.ndarray, jj: np.ndarray,
                        mask2d: np.ndarray, dx: float
                        ) -> Tuple[np.ndarray, np.ndarray]:
    """Hann-windowed streamwise spectrum of ``u'`` on one horizontal plane.

    ``planes`` is ``(nt, ni, nj)``.  The time mean is removed pointwise, solid
    cells are set to zero (their velocity is a constant, so their fluctuation is
    zero anyway), a Hann window is applied over the interior ``x`` span, and the
    result is averaged over ``y`` and over time.

    The interior span is not periodic, so this is a *windowed* spectrum, not a
    Fourier series of a periodic signal.  That is fine for the purpose: the
    identical window and span are used for the child and for the parent
    sub-region, so the windowing bias cancels in the comparison.
    """
    blk = planes[:, ii, :][:, :, jj]
    blk = blk - blk.mean(axis=0, keepdims=True)
    m = mask2d[np.ix_(ii, jj)].astype(float)
    blk = blk * m[None, :, :]
    nx = ii.size
    win = np.hanning(nx)
    win = win / np.sqrt(np.mean(win ** 2))
    spec = np.fft.rfft(blk * win[None, :, None], axis=1)
    power = (np.abs(spec) ** 2).mean(axis=(0, 2)) * (2.0 * dx / nx)
    k = 2.0 * np.pi * np.fft.rfftfreq(nx, d=dx)
    return k, power


# --------------------------------------------------------------------------- #
# The decay length
# --------------------------------------------------------------------------- #


def _first_converged(d: np.ndarray, ok: np.ndarray) -> Optional[float]:
    """Distance from which every cell further in satisfies ``ok``."""
    idx = len(ok)
    while idx > 0 and ok[idx - 1]:
        idx -= 1
    return None if idx == len(ok) else float(d[idx])


def decay_length(distance: np.ndarray, error: np.ndarray, floor: np.ndarray,
                 threshold: float, zone: float) -> Dict[str, object]:
    """Smallest distance beyond the zone at which the error stays acceptable.

    ``error`` and ``floor`` are indexed by the same ``distance``.  The criterion
    is applied only where ``distance >= zone``: inside the zone the solution is
    imposed, so the error there is not evidence of anything.

    Three variants are reported rather than one, because "acceptable" has two
    defensible meanings and they can differ by a lot:

    ``vs_threshold``
        ``error <= threshold`` -- an absolute, run-independent bar.
    ``vs_noise_floor``
        ``error <= floor`` -- indistinguishable from the run's own sampling
        noise.  Note that the floor is estimated from two **independent** halves
        of the parent's window, whereas the child-versus-parent comparison is
        *paired* (same forcing, same geometry, and -- with
        ``nest_linitfromparent`` -- the same initial condition), so its noise is
        smaller than the floor suggests.  Read this variant as a lenient bound.
    ``vs_either``
        ``error <= max(threshold, floor)`` -- the headline number of design
        section 0, which does not ask for more accuracy than the reference
        itself carries.

    ``None`` means the criterion is never met within the available fetch.
    """
    order = np.argsort(distance)
    d, e, f = distance[order], error[order], floor[order]
    keep = d >= zone - 1.0e-9
    d, e, f = d[keep], e[keep], f[keep]
    if d.size == 0:
        return {"decay_length_m": None, "note": "no cells outside the zone"}
    out = {
        "vs_threshold_m": _first_converged(d, e <= threshold),
        "vs_noise_floor_m": _first_converged(d, e <= f),
        "decay_length_m": _first_converged(d, e <= np.maximum(threshold, f)),
        "error_at_zone_edge": float(e[0]),
        "max_error_outside_zone": float(np.max(e)),
        "median_noise_floor_outside_zone": float(np.median(f)),
        "innermost_distance_m": float(d[0]),
    }
    if out["decay_length_m"] is not None and out["decay_length_m"] <= d[0] + 1.0e-9:
        out["note"] = ("no measurable error anywhere outside the zone; the decay "
                       "length is at or below one cell past the ramp")
    elif out["decay_length_m"] is None:
        out["note"] = "never converges within the available fetch"
    return out


# --------------------------------------------------------------------------- #
# V2 reductions -- the falsification metrics
# --------------------------------------------------------------------------- #

#: Wavelength bands, in metres, that section 10.5 located the V1 deficit in.
#: Fixed **physical** bands, deliberately, so that children of different sizes
#: are compared over the same eddies.  ``lambda_gt_quarter_L`` is the one
#: exception: it is the "imposed large scales" band and is by definition
#: relative to the interior span, so it is computed per child and is *not*
#: comparable across the size arm.
SPECTRAL_BANDS: Tuple[Tuple[str, float, float], ...] = (
    ("band_16_64m", 16.0, 64.0),
    ("band_8_16m", 8.0, 16.0),
)


def _band_ratios(k: np.ndarray, e_parent: np.ndarray, e_child: np.ndarray,
                 dx: float, span_m: float) -> Dict[str, object]:
    """Child/parent spectral ratio, band by band.

    Two reductions, because they answer different questions and can differ by
    several per cent:

    ``mean_of_ratios``
        the mean of ``E_child/E_parent`` over the band's wavenumbers.  This is
        the reduction section 10.5 quotes (1.03 / 0.875 / 0.827 / 1.05), so it
        is the headline here too -- V2 has to be readable against V1.
    ``ratio_of_sums``
        the band-integrated energy ratio.  Energy conserving, and less swayed by
        a single noisy high-``k`` bin, but not what the design document quotes.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = np.where(k > 0, 2.0 * np.pi / np.where(k > 0, k, 1.0), np.inf)
    bands = list(SPECTRAL_BANDS) + [
        ("lambda_gt_quarter_L", 0.25 * span_m, np.inf),
        ("lambda_lt_4dx", 0.0, 4.0 * dx),
    ]
    out: Dict[str, object] = {}
    for name, lo, hi in bands:
        sel = (k > 0) & (lam >= lo) & (lam < hi)
        n = int(sel.sum())
        if n == 0 or not np.any(e_parent[sel] > 0):
            out[name] = {"n_modes": n, "mean_of_ratios": None,
                         "ratio_of_sums": None,
                         "wavelength_range_m": [lo, None if np.isinf(hi) else hi]}
            continue
        ratio = e_child[sel] / e_parent[sel]
        out[name] = {
            "n_modes": n,
            "mean_of_ratios": float(np.nanmean(ratio)),
            "ratio_of_sums": float(e_child[sel].sum() / e_parent[sel].sum()),
            "wavelength_range_m": [lo, None if np.isinf(hi) else hi],
        }
    return out


def tke_deficit(prof: Dict[str, List[float]], h: float,
                z_over_h_min: float = 2.0) -> Dict[str, object]:
    """The headline V1 number, as a profile and as one number above ``z/h``.

    ``relative`` is ``(TKE_child - TKE_parent) / TKE_parent`` at each height.
    ``spread`` is the parent's own **half-window spread**,
    ``|A - B| / TKE_parent``, where ``A`` and ``B`` are the two half-window
    estimates.  That is the quantity section 10.5 quotes as the per-height
    spread (1.1-1.5 % above the canopy for the converged V1 run), so
    ``relative / spread`` reproduces its "4-8 sigma".  It is conservative: the
    standard error of the *full*-window mean is about half of it, and the
    child-parent comparison is paired on top of that, so treat significance
    quoted this way as a lower bound.
    """
    z = np.asarray(prof["z"], dtype=float)
    tp = np.asarray(prof["tke_parent"], dtype=float)
    tc = np.asarray(prof["tke_child"], dtype=float)
    ta = np.asarray(prof["tke_parent_halfA"], dtype=float)
    tb = np.asarray(prof["tke_parent_halfB"], dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(tp > 0, (tc - tp) / tp, np.nan)
        spread = np.where(tp > 0, np.abs(ta - tb) / tp, np.nan)
    zh = z / h

    def band(sel: np.ndarray) -> Dict[str, object]:
        empty = {"n_levels": int(sel.sum()), "mean_relative": None,
                 "mean_spread": None, "median_spread": None,
                 "significance": None, "median_significance": None,
                 "spread_range": [None, None]}
        if not np.any(sel) or np.all(np.isnan(rel[sel])):
            return empty
        r = float(np.nanmean(rel[sel]))
        sp = float(np.nanmean(spread[sel]))
        med = float(np.nanmedian(spread[sel]))
        with np.errstate(divide="ignore", invalid="ignore"):
            per = np.where(spread[sel] > 0, rel[sel] / spread[sel], np.nan)
        return {
            "n_levels": int(sel.sum()),
            "mean_relative": r,
            # Two aggregates of the same per-height spread, because they differ
            # a lot and each can mislead on its own.  ``mean_spread`` averages
            # the spread over every level in the band, including the near-lid
            # ones where the resolved TKE is small and its relative sampling
            # error is large -- for the V1 converged child that is 10 % above
            # z/h = 6 against 1.2 % over z/h = 3-5, which drags the aggregate to
            # 4.5 % and the significance down to 2.2 sigma.  ``median_spread``
            # and ``median_significance`` are taken over the levels instead, and
            # give 4.0 sigma, closer to the 4-8 sigma design section 10.5 quotes
            # for z/h = 3-5.  The *deficit* is identical either way; only the
            # uncertainty aggregation differs, so both are reported and neither
            # is chosen after the fact.
            "mean_spread": sp,
            "median_spread": med,
            "spread_range": [float(np.nanmin(spread[sel])),
                             float(np.nanmax(spread[sel]))],
            "significance": (None if sp <= 0 else r / sp),
            "median_significance": (None if np.all(np.isnan(per))
                                    else float(np.nanmedian(per))),
        }

    return {
        "z_over_h": zh.tolist(),
        "absolute": (tc - tp).tolist(),
        "relative": rel.tolist(),
        "half_window_spread": spread.tolist(),
        "above": dict(band(zh >= z_over_h_min), z_over_h_min=z_over_h_min),
        "canopy": dict(band(zh < 1.0), z_over_h_max=1.0),
        "z_over_h_min": z_over_h_min,
    }


def _interp_at(x: np.ndarray, y: np.ndarray, xq: float) -> Optional[float]:
    """Linear interpolation of ``y(x)`` at ``xq``; ``None`` outside the range."""
    if x.size == 0 or xq < x[0] - 1.0e-9 or xq > x[-1] + 1.0e-9:
        return None
    return float(np.interp(xq, x, y))


def tke_error_vs_fetch(curves: Dict[str, Dict[str, List[float]]], preset: Preset,
                       fetches_h: Sequence[float] = (0.5, 1.0, 2.0)) -> Dict[str, object]:
    """The resolved-TKE error as a function of fetch beyond the inner zone edge.

    The abscissa is deliberately *fetch beyond the zone*, not distance from the
    face: the zone arm varies the zone width, so distance-from-the-face would
    compare a point that is 30 m into the interior of one child against a point
    still inside the zone of another.  Reported per face and averaged over the
    four, at the fixed fetches in ``fetches_h`` -- chosen small enough to exist
    for the smallest child in the sweep -- plus each face's own maximum fetch
    and the fetch at which the error first stays at or below the parent's
    sampling floor.
    """
    lz = preset.guardwidth + preset.zonewidth
    h = preset.building_height
    per_face: Dict[str, Dict[str, object]] = {}
    for label, ntot, spacing in (("x", preset.child_itot, preset.dx),
                                 ("y", preset.child_jtot, preset.dy)):
        c = curves[f"{label}_tke"]
        err = np.asarray(c["error"], dtype=float)
        floor = np.asarray(c["noise_floor"], dtype=float)
        for face, dist in (("low", np.asarray(c["distance_from_low_face_m"], float)),
                           ("high", np.asarray(c["distance_from_high_face_m"], float))):
            half = dist <= 0.5 * ntot * spacing + 1.0e-9
            order = np.argsort(dist[half])
            d, e, f = dist[half][order], err[half][order], floor[half][order]
            keep = d >= lz - 1.0e-9
            fetch, e, f = d[keep] - lz, e[keep], f[keep]
            ok = e <= f
            idx = len(ok)
            while idx > 0 and ok[idx - 1]:
                idx -= 1
            crossing = None if idx == len(ok) else float(fetch[idx])
            per_face[f"{label}_{face}"] = {
                "fetch_m": fetch.tolist(),
                "error": e.tolist(),
                "noise_floor": f.tolist(),
                "at_fetch_h": {f"{q:g}h": _interp_at(fetch, e, q * h)
                               for q in fetches_h},
                "max_fetch_h": float(fetch[-1] / h) if fetch.size else None,
                "error_at_max_fetch": float(e[-1]) if e.size else None,
                "error_at_zone_edge": float(e[0]) if e.size else None,
                "crossing_fetch_m": crossing,
                "crossing_fetch_h": None if crossing is None else crossing / h,
                "median_noise_floor": float(np.median(f)) if f.size else None,
            }

    def mean_over_faces(get) -> Optional[float]:
        vals = [get(v) for v in per_face.values()]
        vals = [v for v in vals if v is not None]
        return float(np.mean(vals)) if vals else None

    return {
        "fetches_h": list(fetches_h),
        "per_face": per_face,
        "mean_error_at_fetch_h": {
            f"{q:g}h": mean_over_faces(lambda v, q=q: v["at_fetch_h"][f"{q:g}h"])
            for q in fetches_h},
        "mean_error_at_zone_edge": mean_over_faces(lambda v: v["error_at_zone_edge"]),
        "mean_error_at_max_fetch": mean_over_faces(lambda v: v["error_at_max_fetch"]),
        "max_fetch_h": mean_over_faces(lambda v: v["max_fetch_h"]),
        "faces_crossing_the_floor": sum(
            1 for v in per_face.values() if v["crossing_fetch_h"] is not None),
        "mean_crossing_fetch_h": mean_over_faces(lambda v: v["crossing_fetch_h"]),
    }


def criterion_a(decay: Dict[str, Dict[str, object]], threshold: float = 0.05
                ) -> Dict[str, object]:
    """Design section 0 criterion A: a bound on the free interior of the mean flow.

    ``max_interior |<u>_child - <u>_parent| / u*`` over the four lateral faces,
    against the 0.05 bound.  Included in every V2 point so that a regression in
    the mean flow -- the thing the scheme *does* get right in V1 -- cannot pass
    unnoticed while attention is on the turbulence.
    """
    vals = {k: v.get("max_error_outside_zone")
            for k, v in decay.items() if k.startswith(("x_umean", "y_umean"))}
    finite = [v for v in vals.values() if v is not None and np.isfinite(v)]
    worst = max(finite) if finite else None
    return {
        "per_face": vals,
        "max_interior_umean_error_over_ustar": worst,
        "threshold": threshold,
        "passes": None if worst is None else bool(worst <= threshold),
    }


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #


def _levels_in_window(times: np.ndarray, offset: float, t0: float, t1: float,
                      stride: int) -> List[int]:
    shifted = times - offset
    idx = np.where((shifted >= t0 - 1.0e-9) & (shifted <= t1 + 1.0e-9))[0]
    return list(idx[::stride])


def run(parent_dir: Path, child_dir: Path, outdir: Path, preset: Preset,
        make_plots: bool = True,
        parent_cache: Optional[Dict[Tuple[int, int], Bundle]] = None,
        common_block_cells: int = 0,
        metrics_name: str = "v1_metrics.json") -> Dict[str, object]:
    """Compare one child against the parent sub-region it was cut from.

    ``parent_cache``, when given, is a caller-owned dict that this function
    fills with the accumulated parent :class:`Bundle`, keyed by the child window
    ``(child_itot, child_jtot)``.  A sweep that varies only the zone width drives
    several children out of the *same* parent sub-region, so passing the same
    dict to each call both halves the I/O and guarantees every one of them is
    measured against literally the same parent statistics -- which matters when
    the differences being compared are a few per cent.  The one thing that is
    not window-independent is ``Bundle.tke_series``, an equilibration trace
    averaged over the interior; the cached bundle keeps the interior of whichever
    point built it, recorded as ``tke_series_interior_from`` in the output.

    ``common_block_cells``, when nonzero, adds a second set of profiles over the
    central block of that many cells (:func:`central_indices`), so that children
    with different interiors can be compared over one common region.
    """
    parent_dir, child_dir, outdir = Path(parent_dir), Path(child_dir), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((child_dir / "manifest.json").read_text())
    t_offset = float(manifest["t_offset"])
    t0 = float(manifest["stats_start"])
    t1 = float(manifest["runtime"])

    ni, nj, nk = preset.child_itot - 1, preset.child_jtot - 1, preset.child_ktot - 1
    # Fluid in the child AND fluid in the parent over the same window.
    #
    # For V1 the two are the same mask -- the child is an exact sub-model.  For
    # the V2 size arm they are not: the child clears the cubes that would fall
    # in its guard + ramp band, so inside the band there are cells that are
    # fluid in the child and solid in the parent.  Averaging those with the
    # child's mask alone would fold the parent's near-zero in-building velocity
    # into the parent's statistics and bias the comparison.  Intersecting
    # removes that by construction, everywhere, and is a no-op wherever the
    # geometries agree.
    child_fluid = load_solid_mask(child_dir, (ni, nj, nk))
    parent_fluid_full = load_solid_mask(
        parent_dir, (preset.itot, preset.jtot, preset.ktot))
    parent_fluid = parent_fluid_full[preset.child_i0:preset.child_i0 + ni,
                                     preset.child_j0:preset.child_j0 + nj, :nk]
    mask = child_fluid & parent_fluid
    n_parent_only_solid = int((child_fluid & ~parent_fluid).sum())
    n_child_only_solid = int((parent_fluid & ~child_fluid).sum())
    ii, jj = interior_indices(preset)
    kk = [int(round(z / preset.dz - 0.5)) for z in preset.spectra_heights]
    kk = [k for k in kk if 0 <= k < nk]

    pdump = FieldDump(parent_dir, preset.parent_expnr, preset.dx)
    cdump = FieldDump(child_dir, preset.child_expnr, preset.dx)
    # The parent's stride is scaled by the ratio of the two dump intervals so
    # that both runs are sampled at the same interval (config.Preset
    # .analysis_parent_stride); identical to preset.stride whenever parent and
    # child dump at the same rate, which is every preset before C0b.
    plev = _levels_in_window(pdump.times, t_offset, t0, t1, preset.analysis_parent_stride)
    clev = _levels_in_window(cdump.times, 0.0, t0, t1, preset.stride)
    if len(plev) < 4 or len(clev) < 4:
        raise RuntimeError(
            f"too few samples in the window [{t0}, {t1}] s: "
            f"parent {len(plev)}, child {len(clev)}"
        )

    i0, j0 = preset.child_i0, preset.child_j0
    nif, njf = preset.child_itot, preset.child_jtot

    def cut(u, v, w):
        s = (slice(i0, i0 + nif), slice(j0, j0 + njf), slice(None))
        return u[s], v[s], w[s]

    def whole(u, v, w):
        return u, v, w

    key = (preset.child_itot, preset.child_jtot)
    cached = None if parent_cache is None else parent_cache.get(key)
    if cached is None:
        parent = accumulate(pdump, cut, mask, plev, kk, (ii, jj), "parent sub-region")
        parent.tke_series_interior_from = preset.name
        if parent_cache is not None:
            parent_cache[key] = parent
    else:
        print(f"[analyse] reusing the accumulated parent sub-region for a "
              f"{key[0]} x {key[1]} window (built for '{cached.tke_series_interior_from}')")
        parent = cached
    child = accumulate(cdump, whole, mask, clev, kk, (ii, jj), "nested child")

    metrics = _compare(parent, child, preset, mask, ii, jj, kk, manifest,
                       len(plev), len(clev), common_block_cells)
    metrics["parent_bundle_reused"] = cached is not None
    metrics["solid_mask"] = {
        "fluid_cells_compared": int(mask.sum()),
        "solid_in_the_parent_only": n_parent_only_solid,
        "solid_in_the_child_only": n_child_only_solid,
        "note": ("cells solid in either run are excluded from both; nonzero "
                 "'solid_in_the_parent_only' means the child cleared cubes out of "
                 "its relaxation zone, which is expected in the V2 size arm and "
                 "confined to the band"),
    }
    _write_outputs(outdir, parent, child, preset, mask, ii, jj, kk, metrics,
                   make_plots, metrics_name)
    return metrics


def _profiles_over(parent: Bundle, child: Bundle, mask: np.ndarray,
                   ii: np.ndarray, jj: np.ndarray, zf: np.ndarray
                   ) -> Dict[str, List[float]]:
    """Interior, fluid-only profiles of both runs over one horizontal block."""
    prof: Dict[str, List[float]] = {"z": zf.tolist()}
    for comp in ("u", "v", "w"):
        prof[f"{comp}_parent"] = masked_profile(parent.full_mean(comp), mask, ii, jj).tolist()
        prof[f"{comp}_child"] = masked_profile(child.full_mean(comp), mask, ii, jj).tolist()
    prof["tke_parent"] = masked_profile(parent.tke(), mask, ii, jj).tolist()
    prof["tke_child"] = masked_profile(child.tke(), mask, ii, jj).tolist()
    # Half-window split of the parent: the sampling-noise floor of a profile.
    prof["u_parent_halfA"] = masked_profile(parent.mean["u"][0], mask, ii, jj).tolist()
    prof["u_parent_halfB"] = masked_profile(parent.mean["u"][1], mask, ii, jj).tolist()
    prof["tke_parent_halfA"] = masked_profile(parent.tke(0), mask, ii, jj).tolist()
    prof["tke_parent_halfB"] = masked_profile(parent.tke(1), mask, ii, jj).tolist()
    return prof


def _compare(parent: Bundle, child: Bundle, preset: Preset, mask: np.ndarray,
             ii: np.ndarray, jj: np.ndarray, kk: Sequence[int],
             manifest: Dict, n_parent: int, n_child: int,
             common_block_cells: int = 0) -> Dict[str, object]:
    us = preset.ustar
    zf = (np.arange(preset.child_ktot - 1) + 0.5) * preset.dz
    lz = preset.guardwidth + preset.zonewidth

    prof = _profiles_over(parent, child, mask, ii, jj, zf)

    def prof_err(a: str, b: str, scale: float) -> float:
        x = np.asarray(prof[a]) - np.asarray(prof[b])
        return float(np.sqrt(np.nanmean(x ** 2)) / scale)

    profile_metrics = {
        "u_rms_difference_over_ustar": prof_err("u_child", "u_parent", us),
        "v_rms_difference_over_ustar": prof_err("v_child", "v_parent", us),
        "tke_rms_difference_over_ustar2": prof_err("tke_child", "tke_parent", us ** 2),
        "u_noise_floor_over_ustar": prof_err("u_parent_halfA", "u_parent_halfB", us)
                                     / np.sqrt(2.0),
        "tke_noise_floor_over_ustar2": prof_err("tke_parent_halfA", "tke_parent_halfB",
                                                us ** 2) / np.sqrt(2.0),
    }

    # -- error as a function of distance from each lateral boundary --------- #
    curves: Dict[str, Dict[str, List[float]]] = {}
    decay: Dict[str, Dict[str, object]] = {}
    threshold = 0.05  # of ustar (mean) / ustar^2 (TKE); reported alongside
    for axis, (label, ntot, spacing, other) in enumerate(
            [("x", preset.child_itot, preset.dx, jj),
             ("y", preset.child_jtot, preset.dy, ii)]):
        idx = np.arange(ntot - 1)
        dist_lo = (idx + 0.5) * spacing
        dist_hi = (ntot - idx - 0.5) * spacing
        for name, cscale, pfun in (
                ("umean", us, lambda b: b.full_mean("u")),
                ("tke", us ** 2, lambda b: b.tke())):
            err = slab_rms_difference(pfun(child), pfun(parent), mask, axis, other) / cscale
            floor = slab_rms_difference(
                {"umean": lambda b: b.mean["u"][0], "tke": lambda b: b.tke(0)}[name](parent),
                {"umean": lambda b: b.mean["u"][1], "tke": lambda b: b.tke(1)}[name](parent),
                mask, axis, other) / cscale / np.sqrt(2.0)
            key = f"{label}_{name}"
            curves[key] = {
                "index": idx.tolist(),
                "distance_from_low_face_m": dist_lo.tolist(),
                "distance_from_high_face_m": dist_hi.tolist(),
                "error": err.tolist(),
                "noise_floor": floor.tolist(),
            }
            half = dist_lo <= 0.5 * ntot * spacing
            decay[f"{key}_low_face"] = decay_length(
                dist_lo[half], err[half], floor[half], threshold, lz)
            decay[f"{key}_high_face"] = decay_length(
                dist_hi[~half], err[~half], floor[~half], threshold, lz)

    # -- spectra ------------------------------------------------------------ #
    def spectra_over(si: np.ndarray, sj: np.ndarray) -> Dict[str, Dict[str, object]]:
        out: Dict[str, Dict[str, object]] = {}
        span = si.size * preset.dx
        for k in kk:
            kp, pw = streamwise_spectrum(parent.planes[k], si, sj, mask[:, :, k], preset.dx)
            _, pc = streamwise_spectrum(child.planes[k], si, sj, mask[:, :, k], preset.dx)
            z = (k + 0.5) * preset.dz
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(pw > 0, pc / pw, np.nan)
            out[f"z_{z:g}m"] = {
                "z_m": z,
                "z_over_h": z / preset.building_height,
                "span_m": span,
                "wavenumber_rad_per_m": kp.tolist(),
                "E_parent": pw.tolist(),
                "E_child": pc.tolist(),
                "child_over_parent": ratio.tolist(),
                "band_mean_ratio_resolved": float(
                    np.nanmean(ratio[1:max(2, len(ratio) // 2)])),
                "bands": _band_ratios(kp, pw, pc, preset.dx, span),
            }
        return out

    spectra = spectra_over(ii, jj)

    # -- the V2 falsification block ----------------------------------------- #
    h = preset.building_height
    common: Dict[str, object] = {"cells": common_block_cells}
    if common_block_cells:
        ci, cj = central_indices(preset, common_block_cells)
        cprof = _profiles_over(parent, child, mask, ci, cj, zf)
        common.update({
            "extent_m": common_block_cells * preset.dx,
            "extent_h": common_block_cells * preset.dx / h,
            "profiles": cprof,
            "tke_deficit": tke_deficit(cprof, h),
            "spectra": spectra_over(ci, cj),
        })

    v2 = {
        "configuration": {
            "N_imp_cells": int(round(preset.guardwidth / preset.dx)),
            "N_rel_cells": int(round(preset.zonewidth / preset.dx)),
            "zone_cells": preset.zone_cells,
            "nzone": preset.nzone,
            "tau_s": preset.tau,
            "optical_depth": preset.optical_depth,
            "cadence_s": preset.cadence,
            "cadence_stride": preset.cadence_stride,
            "C_dump_at_u0": preset.c_dump_u0,
            "nest_timeinterp": preset.timeinterp,
            "parent_dtdump_s": preset.dtdump,
            "child_dtdump_s": preset.child_dtdump,
            "child_cells": [preset.child_itot, preset.child_jtot],
            "child_extent_m": [preset.child_xlen, preset.child_ylen],
            "interior_cells": preset.interior_cells,
            "interior_extent_m": preset.interior_extent_m,
            "interior_extent_h": preset.interior_extent_h,
            "zone_fraction": preset.zone_fraction,
            "zone_fraction_warns": preset.zone_fraction_warns,
            "building_free_zone": preset.building_free_zone,
            "nest_lparentgeom": not preset.building_free_zone,
            "clear_child_zone": bool(preset.clear_child_zone),
            "n_cubes_child": int(len(preset.child_cube_centres())),
            "n_cubes_cleared_from_child_zone": int(preset.n_child_cubes_removed),
            "n_cubes_in_zone": int(len(preset.cubes_in_zone())),
            "n_cubes_in_interior": int(len(preset.cubes_in_analysis_interior())),
            "building_clearance_available_m": preset.building_clearance_available,
            "zone_clearance_needed_m": preset.zone_clearance,
        },
        "tke_deficit": tke_deficit(prof, h),
        "tke_error_vs_fetch": tke_error_vs_fetch(curves, preset),
        "criterion_a": criterion_a(decay, threshold),
        "common_block": common,
    }

    return {
        "preset": preset.name,
        "samples": {"parent": n_parent, "child": n_child},
        "tke_series_interior_from": parent.tke_series_interior_from,
        "v2": v2,
        "window_s": [manifest["stats_start"], manifest["runtime"]],
        "ustar": us,
        "building_height_m": preset.building_height,
        "zone_width_m": lz,
        "threshold_used": threshold,
        "nest_timeinterp": manifest.get("nest_timeinterp"),
        "init_from_parent": manifest.get("init_from_parent"),
        "flux_residual_after_correction": manifest.get("flux_residual_after_correction"),
        "profiles": prof,
        "profile_metrics": profile_metrics,
        "error_vs_distance": curves,
        "decay_lengths": decay,
        "decay_lengths_in_building_heights": {
            k: (None if v.get("decay_length_m") is None
                else v["decay_length_m"] / preset.building_height)
            for k, v in decay.items()
        },
        "spectra": spectra,
        "tke_series": {
            "parent_times_s": (parent.times - parent.times[0]).tolist(),
            "parent": parent.tke_series.tolist(),
            "child_times_s": child.times.tolist(),
            "child": child.tke_series.tolist(),
        },
    }


def _write_outputs(outdir: Path, parent: Bundle, child: Bundle, preset: Preset,
                   mask: np.ndarray, ii: np.ndarray, jj: np.ndarray,
                   kk: Sequence[int], metrics: Dict, make_plots: bool,
                   metrics_name: str = "v1_metrics.json") -> None:
    (outdir / metrics_name).write_text(json.dumps(metrics, indent=2) + "\n",
                                       encoding="ascii")

    d = metrics["v2"]["tke_deficit"]
    _csv(outdir / "tke_deficit.csv",
         ["z_over_h", "tke_parent", "tke_child", "difference", "relative",
          "half_window_spread"],
         [d["z_over_h"], metrics["profiles"]["tke_parent"],
          metrics["profiles"]["tke_child"], d["absolute"], d["relative"],
          d["half_window_spread"]])

    for key, v in metrics["v2"]["tke_error_vs_fetch"]["per_face"].items():
        _csv(outdir / f"tke_error_vs_fetch_{key}.csv",
             ["fetch_m", "error", "noise_floor"],
             [v["fetch_m"], v["error"], v["noise_floor"]])

    rows = []
    for scope in ("interior", "common"):
        spec = (metrics["spectra"] if scope == "interior"
                else metrics["v2"]["common_block"].get("spectra", {}))
        for name, sp in spec.items():
            for band, b in sp["bands"].items():
                rows.append((scope, name, sp["z_over_h"], band, b["n_modes"],
                             b["mean_of_ratios"], b["ratio_of_sums"]))
    if rows:
        with (outdir / "spectral_bands.csv").open("w", encoding="ascii",
                                                  newline="\n") as fh:
            fh.write("scope,height,z_over_h,band,n_modes,mean_of_ratios,"
                     "ratio_of_sums\n")
            for r in rows:
                fh.write(",".join("" if x is None else
                                  (f"{x:.9g}" if isinstance(x, float) else str(x))
                                  for x in r) + "\n")

    prof = metrics["profiles"]
    cols = ["z", "u_parent", "u_child", "v_parent", "v_child",
            "tke_parent", "tke_child"]
    _csv(outdir / "profiles.csv", cols, [prof[c] for c in cols])

    for key, c in metrics["error_vs_distance"].items():
        _csv(outdir / f"error_vs_distance_{key}.csv",
             ["index", "distance_from_low_face_m", "distance_from_high_face_m",
              "error", "noise_floor"],
             [c["index"], c["distance_from_low_face_m"],
              c["distance_from_high_face_m"], c["error"], c["noise_floor"]])

    for key, s in metrics["spectra"].items():
        _csv(outdir / f"spectrum_{key}.csv",
             ["wavenumber_rad_per_m", "E_parent", "E_child"],
             [s["wavenumber_rad_per_m"], s["E_parent"], s["E_child"]])

    if make_plots:
        try:
            _plots(outdir, preset, metrics)
        except Exception as exc:  # pragma: no cover - plotting is a convenience
            print(f"[analyse] plotting skipped: {exc}")


def _csv(path: Path, header: Sequence[str], columns: Sequence[Sequence[float]]) -> None:
    rows = zip(*columns)
    with path.open("w", encoding="ascii", newline="\n") as fh:
        fh.write(",".join(header) + "\n")
        for row in rows:
            fh.write(",".join("" if v is None or (isinstance(v, float) and np.isnan(v))
                              else f"{v:.9g}" for v in row) + "\n")


def _plots(outdir: Path, preset: Preset, metrics: Dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    h = preset.building_height
    prof = metrics["profiles"]
    z = np.asarray(prof["z"])

    fig, ax = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    for a, keys, xlabel in (
            (ax[0], ("u_parent", "u_child"), r"$\langle u\rangle$ [m s$^{-1}$]"),
            (ax[1], ("v_parent", "v_child"), r"$\langle v\rangle$ [m s$^{-1}$]"),
            (ax[2], ("tke_parent", "tke_child"), r"resolved TKE [m$^2$ s$^{-2}$]")):
        a.plot(prof[keys[0]], z / h, "k-", label="parent sub-region")
        a.plot(prof[keys[1]], z / h, "r--", label="nested child")
        a.axhline(1.0, color="0.7", lw=0.8)
        a.set_xlabel(xlabel)
        a.set_ylabel(r"$z/h$")
        a.grid(alpha=0.3)
    ax[0].legend(fontsize=8)
    fig.suptitle(f"V1 Big Brother, preset '{preset.name}': interior profiles")
    fig.savefig(outdir / "profiles.png", dpi=130)
    plt.close(fig)

    keys = [k for k in metrics["error_vs_distance"] if k.startswith("x_")]
    fig, ax = plt.subplots(1, len(keys), figsize=(5 * len(keys), 4),
                           squeeze=False, constrained_layout=True)
    for a, key in zip(ax[0], keys):
        c = metrics["error_vs_distance"][key]
        d = np.asarray(c["distance_from_low_face_m"])
        a.semilogy(d / h, c["error"], "r-", label="child vs parent")
        a.semilogy(d / h, c["noise_floor"], "k:", label="sampling-noise floor")
        a.axvline(metrics["zone_width_m"] / h, color="b", ls="--",
                  label="inner edge of the zone")
        a.axvline((preset.child_xlen - metrics["zone_width_m"]) / h, color="b", ls="--")
        a.axhline(metrics["threshold_used"], color="0.6", lw=0.8)
        a.set_xlabel(r"distance from the west face, $x/h$")
        a.set_ylabel("normalised error")
        a.set_title(key)
        a.grid(alpha=0.3)
    ax[0][0].legend(fontsize=8)
    fig.suptitle("Error versus distance from the lateral boundary")
    fig.savefig(outdir / "error_vs_distance.png", dpi=130)
    plt.close(fig)

    spectra = metrics["spectra"]
    if spectra:
        fig, ax = plt.subplots(1, len(spectra), figsize=(5 * len(spectra), 4),
                               squeeze=False, constrained_layout=True)
        for a, (name, s) in zip(ax[0], spectra.items()):
            k = np.asarray(s["wavenumber_rad_per_m"])
            a.loglog(k[1:], np.asarray(s["E_parent"])[1:], "k-", label="parent")
            a.loglog(k[1:], np.asarray(s["E_child"])[1:], "r--", label="child")
            a.set_xlabel(r"$k_x$ [rad m$^{-1}$]")
            a.set_ylabel(r"$E_{uu}$")
            a.set_title(f"z = {s['z_m']:g} m  (z/h = {s['z_over_h']:.2f})")
            a.grid(alpha=0.3, which="both")
        ax[0][0].legend(fontsize=8)
        fig.suptitle("Streamwise spectra of $u'$ over the interior")
        fig.savefig(outdir / "spectra.png", dpi=130)
        plt.close(fig)

    v2 = metrics.get("v2")
    if v2:
        d = v2["tke_deficit"]
        zh = np.asarray(d["z_over_h"], dtype=float)
        rel = 100.0 * np.asarray(d["relative"], dtype=float)
        spread = 100.0 * np.asarray(d["half_window_spread"], dtype=float)
        fig, ax = plt.subplots(1, 2, figsize=(9, 4.5), constrained_layout=True)
        ax[0].fill_betweenx(zh, -spread, spread, color="0.85",
                            label="parent half-window spread")
        ax[0].plot(rel, zh, "r-", label="child - parent")
        ax[0].axvline(0.0, color="k", lw=0.8)
        ax[0].axhline(d["above"]["z_over_h_min"], color="b", ls="--", lw=0.8,
                      label=f"z/h = {d['above']['z_over_h_min']:g}")
        ax[0].axhline(1.0, color="0.7", lw=0.8)
        ax[0].set_xlabel("resolved-TKE difference [%]")
        ax[0].set_ylabel(r"$z/h$")
        ax[0].legend(fontsize=8)
        ax[0].grid(alpha=0.3)
        for key, v in v2["tke_error_vs_fetch"]["per_face"].items():
            ax[1].plot(np.asarray(v["fetch_m"]) / h, v["error"], lw=1.0, label=key)
        ax[1].plot(np.asarray(v["fetch_m"]) / h, v["noise_floor"], "k:",
                   label="sampling floor")
        ax[1].set_xlabel(r"fetch beyond the inner zone edge, $/h$")
        ax[1].set_ylabel(r"normalised TKE error [$u_\star^2$]")
        ax[1].legend(fontsize=7)
        ax[1].grid(alpha=0.3)
        cfg = v2["configuration"]
        fig.suptitle(f"{preset.name}: N_rel = {cfg['N_rel_cells']} cells, child "
                     f"{cfg['child_cells'][0]}$^2$, interior "
                     f"{cfg['interior_extent_h']:.2f}$h$")
        fig.savefig(outdir / "tke_deficit.png", dpi=130)
        plt.close(fig)

    ts = metrics["tke_series"]
    fig, a = plt.subplots(figsize=(6, 4), constrained_layout=True)
    a.plot(ts["parent_times_s"], ts["parent"], "k-", label="parent sub-region")
    a.plot(ts["child_times_s"], ts["child"], "r-", label="nested child")
    a.axvline(metrics["window_s"][0], color="b", ls="--", label="statistics start")
    a.set_xlabel("child time [s]")
    a.set_ylabel(r"interior mean $\frac{1}{2}\langle u_iu_i\rangle$ [m$^2$ s$^{-2}$]")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)
    fig.savefig(outdir / "tke_series.png", dpi=130)
    plt.close(fig)


def summary(metrics: Dict) -> str:
    m = metrics["profile_metrics"]
    lines = [
        f"V1 Big Brother, preset '{metrics['preset']}'",
        f"  samples: parent {metrics['samples']['parent']}, "
        f"child {metrics['samples']['child']}, window {metrics['window_s']} s",
        f"  nest_timeinterp = {metrics['nest_timeinterp']}, "
        f"init_from_parent = {metrics['init_from_parent']}",
        "",
        "  interior profiles (rms difference over the column):",
        f"    <u>   {m['u_rms_difference_over_ustar']:.4f} u*   "
        f"(sampling-noise floor {m['u_noise_floor_over_ustar']:.4f} u*)",
        f"    <v>   {m['v_rms_difference_over_ustar']:.4f} u*",
        f"    TKE   {m['tke_rms_difference_over_ustar2']:.4f} u*^2 "
        f"(floor {m['tke_noise_floor_over_ustar2']:.4f} u*^2)",
        "",
        "  decay lengths (distance at which the error falls to the floor "
        "and stays there):",
    ]
    for key, val in metrics["decay_lengths"].items():
        d = val.get("decay_length_m")
        inh = metrics["decay_lengths_in_building_heights"][key]
        head = "NOT REACHED" if d is None else f"{d:8.1f} m = {inh:5.2f} h"
        lines.append(
            f"    {key:26s} {head:>22s}   "
            f"err@zone-edge {val.get('error_at_zone_edge', float('nan')):.4f}, "
            f"max outside {val.get('max_error_outside_zone', float('nan')):.4f}, "
            f"floor {val.get('median_noise_floor_outside_zone', float('nan')):.4f}")
    lines.append("")
    lines.append("  spectra (child/parent power, mean over the resolved band):")
    for name, s in metrics["spectra"].items():
        b = s["bands"]
        def band(key):
            v = b[key]["mean_of_ratios"]
            return "   --  " if v is None else f"{v:6.3f}"
        lines.append(f"    {name:14s} z/h = {s['z_over_h']:.2f}   "
                     f"resolved {s['band_mean_ratio_resolved']:.3f}   "
                     f"16-64 m {band('band_16_64m')}   "
                     f"8-16 m {band('band_8_16m')}   "
                     f"> L/4 {band('lambda_gt_quarter_L')}   "
                     f"< 4dx {band('lambda_lt_4dx')}")
    v2 = metrics.get("v2")
    if v2:
        cfg, d = v2["configuration"], v2["tke_deficit"]
        f = v2["tke_error_vs_fetch"]
        a = v2["criterion_a"]

        def pct(x):
            return "  --  " if x is None else f"{100 * x:+6.2f}%"

        lines += [
            "",
            f"  V2: N_imp+N_rel = {cfg['N_imp_cells']}+{cfg['N_rel_cells']} cells, "
            f"child {cfg['child_cells'][0]}x{cfg['child_cells'][1]}, "
            f"interior {cfg['interior_cells']} cells = {cfg['interior_extent_h']:.2f} h, "
            f"zone {'building-free' if cfg['building_free_zone'] else 'CONTAINS BUILDINGS'} "
            f"(nest_lparentgeom = {'.false.' if cfg['building_free_zone'] else '.true.'})",
            f"    resolved-TKE deficit above z/h = {d['above']['z_over_h_min']:g}: "
            f"{pct(d['above']['mean_relative'])} against a half-window spread of "
            f"{pct(d['above']['mean_spread']).strip('+')} mean / "
            f"{pct(d['above']['median_spread']).strip('+')} median"
            + ("" if d['above']['significance'] is None
               else f"  ({d['above']['significance']:+.1f} sigma on the mean, "
                    f"{d['above']['median_significance']:+.1f} per-height median)"),
            f"    inside the canopy (z/h < 1):        {pct(d['canopy']['mean_relative'])} "
            f"against {pct(d['canopy']['mean_spread']).strip('+')}",
            "    TKE error vs fetch beyond the zone: "
            + ", ".join(f"{k} {('--' if v is None else format(v, '.3f'))}"
                        for k, v in f["mean_error_at_fetch_h"].items())
            + "; at max fetch ("
            + ("--" if f["max_fetch_h"] is None else format(f["max_fetch_h"], ".2f"))
            + " h) "
            + ("--" if f["mean_error_at_max_fetch"] is None
               else format(f["mean_error_at_max_fetch"], ".3f")),
            f"    faces whose TKE error reaches the sampling floor: "
            f"{f['faces_crossing_the_floor']}/4",
            "    criterion A (mean flow, interior): "
            + ("--" if a["max_interior_umean_error_over_ustar"] is None
               else format(a["max_interior_umean_error_over_ustar"], ".4f"))
            + f" u* against {a['threshold']:g} -- "
            + ("PASS" if a["passes"] else "FAIL"),
        ]
        cb = v2["common_block"]
        if cb.get("tke_deficit") and cb.get("extent_h") is not None:
            lines.append(
                f"    over the common {cb['cells']}-cell block ({cb['extent_h']:.2f} h): "
                f"deficit {pct(cb['tke_deficit']['above']['mean_relative'])}")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("parent_dir", type=Path)
    parser.add_argument("child_dir", type=Path)
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--preset", default="production")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--common-block", type=int, default=0,
                        help="also compare over this many central cells (V2)")
    parser.add_argument("--metrics-name", default="v1_metrics.json")
    args = parser.parse_args()

    preset = get_preset(args.preset)
    metrics = run(args.parent_dir, args.child_dir, args.outdir, preset,
                  make_plots=not args.no_plots,
                  common_block_cells=args.common_block,
                  metrics_name=args.metrics_name)
    print(summary(metrics))
    print(f"\nwritten to {args.outdir}")


if __name__ == "__main__":
    main()
