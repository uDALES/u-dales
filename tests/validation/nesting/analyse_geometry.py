#!/usr/bin/env python3
"""The V3 / V4 measurements: adjustment length, and the cost of a layout mismatch.

Why this is not part of ``analyse.py``
--------------------------------------

``analyse.py`` answers one question -- *does the child reproduce the parent
sub-region it was cut from?* -- and every number in it is a child-minus-parent
difference over cells that are fluid in **both**.  That is the right question
for V1 and V2, where the child is a sub-model of its parent.  It is the wrong
question here, and getting that wrong would be the easiest way to produce a
confident, meaningless answer:

* in **V3** the parent has no buildings at all, so ``<u>_child - <u>_parent``
  measures the canopy, not the nesting.  Criterion A will be far outside its
  bound and that is not a failure of anything;
* in **V4** the parent's buildings are in *different places*, so the
  fluid-in-both mask throws away roughly a quarter of the child's canopy layer,
  in a pattern set by the parent's array -- an averaging domain that has nothing
  to do with the one V1 averaged over.

So the child-versus-parent block is still computed and still reported (the brief
asks for the V1/V2 diagnostics for comparability, and for V4 it *is* interesting
as a measure of how far the interior is free to differ), but it carries no pass
criterion here.  The references are:

* **V3** -- a periodic run of the child's own canopy at the child's own forcing.
  That is what "in equilibrium with its own canopy" means, and without it an
  adjustment length could only ever be self-referential.
* **V4** -- the V1 ``converged`` child, which has *identical geometry, size,
  zone, forcing and schedule* to the V4 child.  Both are reduced here by the
  same code over each child's **own** fluid mask, so the two averaging domains
  are the same cells.

Everything below is measured over the child's interior -- outside the guard
strip and the relaxation ramp -- exactly as ``analyse.py`` defines it.

Usage
-----
    python analyse_geometry.py --help

is deliberately thin; ``run_geometry.py`` is the entry point.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

import analyse
from analyse import _band_ratios, interior_indices, masked_profile, slab_rms_difference
from caselib import FieldDump, cell_centred, load_solid_mask
from config import Preset

#: The bound design section 0 puts on the free interior of the mean flow, reused
#: here as the tolerance on "adjusted": a row counts as in equilibrium with the
#: canopy when its canopy-layer statistic is within 5 % of the equilibrium run's.
#: Fixed before any number was computed, and the same 0.05 criterion A uses, so
#: it is not a knob.
TOLERANCE = 0.05


# --------------------------------------------------------------------------- #
# Accumulation
# --------------------------------------------------------------------------- #


@dataclass
class GeoStats:
    """Time-accumulated statistics of one run on its own cell centres.

    The same half-window pair as :class:`analyse.Bundle`, plus the ``u'w'``
    cross-moment, which the adjustment measurement needs and ``analyse.py`` does
    not carry: the resolved shear stress at roof height is the sharpest
    streamwise signature of an internal boundary layer, and it is the direct
    proxy for the "visibly wrong facet stresses" section 9.4 expects on the
    first few building rows.
    """

    label: str
    times: np.ndarray
    counts: Tuple[int, int] = (0, 0)
    mean: Dict[str, List[np.ndarray]] = field(default_factory=dict)
    msq: Dict[str, List[np.ndarray]] = field(default_factory=dict)
    #: [half] -> mean of u*w
    cross_uw: List[np.ndarray] = field(default_factory=list)
    #: height index -> (nt, ni, nj) time series of u, for spectra
    planes: Dict[int, np.ndarray] = field(default_factory=dict)
    #: fluid-masked, depth-averaged <u> at each sampled time
    bulk_series: np.ndarray = field(default_factory=lambda: np.zeros(0))

    def _full(self, per_half: Sequence[np.ndarray]) -> np.ndarray:
        na, nb = self.counts
        return (na * per_half[0] + nb * per_half[1]) / (na + nb)

    def full_mean(self, comp: str) -> np.ndarray:
        return self._full(self.mean[comp])

    def variance(self, comp: str, half: Optional[int] = None) -> np.ndarray:
        if half is None:
            return np.maximum(self._full(self.msq[comp])
                              - self.full_mean(comp) ** 2, 0.0)
        return np.maximum(self.msq[comp][half] - self.mean[comp][half] ** 2, 0.0)

    def tke(self, half: Optional[int] = None) -> np.ndarray:
        return 0.5 * sum(self.variance(c, half) for c in ("u", "v", "w"))

    def uw(self, half: Optional[int] = None) -> np.ndarray:
        """Resolved kinematic shear stress ``<u'w'> = <uw> - <u><w>``."""
        if half is None:
            return self._full(self.cross_uw) - self.full_mean("u") * self.full_mean("w")
        return self.cross_uw[half] - self.mean["u"][half] * self.mean["w"][half]

    def mean_half(self, comp: str, half: int) -> np.ndarray:
        return self.mean[comp][half]


def accumulate_geo(dump: FieldDump, extract: Callable, levels: Sequence[int],
                   fluid: np.ndarray, plane_k: Sequence[int] = (),
                   label: str = "") -> GeoStats:
    """One streaming pass, building a :class:`GeoStats`.

    ``fluid`` is used only for the ``bulk_series`` trace; the 3-D fields are
    accumulated everywhere and masked at reduction time, exactly as
    ``analyse.accumulate`` does.
    """
    levels = list(levels)
    n = len(levels)
    half = n // 2
    sums: List[Dict[str, np.ndarray]] = []
    sqs: List[Dict[str, np.ndarray]] = []
    cross: List[np.ndarray] = []
    counts = [0, 0]
    planes: Dict[int, List[np.ndarray]] = {int(k): [] for k in plane_k}
    bulk = np.zeros(n)
    shape = None
    ncell_z = None

    for idx, lev in enumerate(levels):
        u, v, w = extract(*dump.read_level(lev))
        uc, vc, wc = cell_centred(u, v, w)
        if shape is None:
            shape = uc.shape
            for _ in range(2):
                sums.append({c: np.zeros(shape) for c in "uvw"})
                sqs.append({c: np.zeros(shape) for c in "uvw"})
                cross.append(np.zeros(shape))
            ncell_z = fluid.sum(axis=(0, 1)).astype(float)
        h = 0 if idx < half else 1
        counts[h] += 1
        for c, arr in zip("uvw", (uc, vc, wc)):
            sums[h][c] += arr
            sqs[h][c] += arr * arr
        cross[h] += uc * wc
        for k in planes:
            planes[k].append(uc[:, :, k].copy())
        # The volume-flow-rate controller averages u over the fluid cells of each
        # level and then over the depth (modforces.f90:404-410), so the bulk is
        # computed the same way here -- it is what V3's flat parent is asked to
        # hold, so it has to be the same number.
        plane = np.where(fluid, uc, 0.0).sum(axis=(0, 1))
        bulk[idx] = float(np.sum(np.where(ncell_z > 0, plane / np.maximum(ncell_z, 1),
                                          0.0)) / max(shape[2], 1))

    st = GeoStats(label=label, times=np.asarray([dump.times[l] for l in levels]))
    st.counts = (counts[0], counts[1])
    for c in "uvw":
        st.mean[c] = [sums[h][c] / max(counts[h], 1) for h in range(2)]
        st.msq[c] = [sqs[h][c] / max(counts[h], 1) for h in range(2)]
    st.cross_uw = [cross[h] / max(counts[h], 1) for h in range(2)]
    st.planes = {k: np.stack(v) for k, v in planes.items()}
    st.bulk_series = bulk
    return st


# --------------------------------------------------------------------------- #
# Reductions
# --------------------------------------------------------------------------- #


def _canopy_k(preset: Preset) -> Tuple[np.ndarray, int]:
    """Cell indices inside the canopy, and the first level above roof height."""
    nk = preset.child_ktot - 1
    zf = (np.arange(nk) + 0.5) * preset.dz
    inside = np.where(zf < preset.building_height)[0]
    above = int(min(int(preset.building_height / preset.dz), nk - 1))
    return inside, above


def domain_profiles(st: GeoStats, fluid: np.ndarray, preset: Preset,
                    ii: Optional[np.ndarray] = None,
                    jj: Optional[np.ndarray] = None) -> Dict[str, object]:
    """Horizontal, fluid-only profiles plus the canopy-layer reductions.

    With ``ii``/``jj`` omitted the whole domain is used, which is what a
    periodic run wants.
    """
    ni, nj = fluid.shape[0], fluid.shape[1]
    ii = np.arange(ni) if ii is None else ii
    jj = np.arange(nj) if jj is None else jj
    nk = fluid.shape[2]
    zf = (np.arange(nk) + 0.5) * preset.dz
    canopy, k_h = _canopy_k(preset)

    out: Dict[str, object] = {"z": zf.tolist(),
                              "z_over_h": (zf / preset.building_height).tolist()}
    for name, fld in (("u", st.full_mean("u")), ("v", st.full_mean("v")),
                      ("w", st.full_mean("w")), ("tke", st.tke()),
                      ("uw", st.uw())):
        out[name] = masked_profile(fld, fluid, ii, jj).tolist()
    for h in (0, 1):
        out[f"u_half{h}"] = masked_profile(st.mean_half("u", h), fluid, ii, jj).tolist()
        out[f"tke_half{h}"] = masked_profile(st.tke(h), fluid, ii, jj).tolist()

    u = np.asarray(out["u"], dtype=float)
    tke = np.asarray(out["tke"], dtype=float)
    uw = np.asarray(out["uw"], dtype=float)
    out["u_canopy"] = float(np.nanmean(u[canopy])) if canopy.size else None
    out["tke_canopy"] = float(np.nanmean(tke[canopy])) if canopy.size else None
    out["uw_at_roof"] = float(uw[k_h])
    out["z_roof_m"] = float(zf[k_h])
    # Depth-averaged, fluid-masked <u>: the quantity luvolflowr holds.
    out["bulk_u"] = float(np.nanmean(u))
    out["bulk_series_mean"] = float(np.mean(st.bulk_series))
    out["bulk_series_halfspread"] = float(abs(
        np.mean(st.bulk_series[:len(st.bulk_series) // 2])
        - np.mean(st.bulk_series[len(st.bulk_series) // 2:])))
    return out


def block_edges(preset) -> np.ndarray:
    """Streamwise analysis blocks, one cube period wide, in child metres.

    Blocks are phase-locked to the child's canopy lattice -- their edges are the
    canopy's building faces, continued *upstream* into the building-free
    standoff -- and clipped to the interior.  Two reasons, and both matter:

    * ``<u>(x)`` inside a cube array swings by tens of per cent within one
      period (in front of a cube, over it, in its wake).  The adjustment is the
      trend *through* that oscillation, so each block must contain exactly one
      cube and the same part of the pattern, or the streamwise signal is
      dominated by aliasing against the array.
    * continuing the same phase upstream means the standoff region is measured
      on the same abscissa as the canopy, which is what makes a 0-cell and a
      40-cell standoff comparable at a fixed distance from the zone.
    """
    p = preset.period
    lz = preset.guardwidth + preset.zonewidth
    span = preset.canopy_x_range
    if span is None:
        return np.zeros(0)
    x_first = span[0]
    k0 = int(np.ceil((lz - x_first) / p - 1.0e-9))
    k1 = int(np.floor((preset.child_xlen - lz - x_first) / p + 1.0e-9)) - 1
    if k1 < k0:
        return np.zeros(0)
    return x_first + p * np.arange(k0, k1 + 2)


def streamwise_blocks(st: GeoStats, fluid: np.ndarray, preset,
                      y_scope: str = "core") -> List[Dict[str, object]]:
    """Per-block canopy-layer statistics against fetch.

    ``y_scope`` is ``"core"`` (the central :attr:`y_core_fraction` of the
    canopy, the headline) or ``"canopy"`` (its whole spanwise extent).  Both are
    emitted so that the lateral internal boundary layers spreading in from the
    two spanwise zone edges can be seen rather than assumed away.
    """
    p = preset.period
    lz = preset.guardwidth + preset.zonewidth
    ni, nj = fluid.shape[0], fluid.shape[1]
    xc = (np.arange(ni) + 0.5) * preset.dx
    yc = (np.arange(nj) + 0.5) * preset.dy
    span = preset.y_core_range if y_scope == "core" else preset.canopy_y_range
    if span is None:
        return []
    jj = np.where((yc >= span[0] - 1.0e-9) & (yc <= span[1] + 1.0e-9))[0]
    canopy_k, k_h = _canopy_k(preset)
    edges = block_edges(preset)
    cubes = preset.child_cube_centres()
    cube_x = np.unique(np.round(cubes[:, 0], 6)) if cubes.size else np.zeros(0)

    umean = st.full_mean("u")
    tke = st.tke()
    uw = st.uw()
    umean_half = [st.mean_half("u", 0), st.mean_half("u", 1)]
    out: List[Dict[str, object]] = []
    for b in range(len(edges) - 1):
        lo, hi = edges[b], edges[b + 1]
        ii = np.where((xc >= lo - 1.0e-9) & (xc < hi - 1.0e-9))[0]
        if ii.size == 0 or jj.size == 0:
            continue
        u_p = masked_profile(umean, fluid, ii, jj)
        t_p = masked_profile(tke, fluid, ii, jj)
        w_p = masked_profile(uw, fluid, ii, jj)
        has_cube = bool(np.any((cube_x >= lo - 1.0e-9) & (cube_x < hi - 1.0e-9)))
        # Half-window estimates of the same block statistic: the sampling scale
        # every difference between standoffs has to be read against.
        uc_half = [float(np.nanmean(masked_profile(umean_half[h], fluid, ii, jj)
                                    [canopy_k])) if canopy_k.size else None
                   for h in (0, 1)]
        uc_full = float(np.nanmean(u_p[canopy_k])) if canopy_k.size else None
        out.append({
            "block": b,
            "x_lo_m": float(lo), "x_hi_m": float(hi),
            "fetch_from_zone_m": float(lo - lz),
            "fetch_from_zone_h": float((lo - lz) / preset.building_height),
            "fetch_from_first_row_m": float(lo - preset.canopy_x_range[0]),
            "fetch_from_first_row_h": float(
                (lo - preset.canopy_x_range[0]) / preset.building_height),
            "has_cube": has_cube,
            "n_fluid_cells": int(fluid[np.ix_(ii, jj)].sum()),
            "u": u_p.tolist(), "tke": t_p.tolist(), "uw": w_p.tolist(),
            "u_canopy": uc_full,
            "u_canopy_half0": uc_half[0], "u_canopy_half1": uc_half[1],
            "u_canopy_spread": (
                None if (uc_full in (None, 0.0) or uc_half[0] is None)
                else abs(uc_half[0] - uc_half[1]) / abs(uc_full)),
            "tke_canopy": float(np.nanmean(t_p[canopy_k])) if canopy_k.size else None,
            "uw_at_roof": float(w_p[k_h]),
        })
    return out


def _first_settled(fetch: np.ndarray, err: np.ndarray, tol: float) -> Optional[float]:
    """Smallest fetch beyond which **every** later value satisfies ``err <= tol``.

    The same "and stays there" rule ``analyse.decay_length`` uses, so an
    adjustment length here means what a decay length means there.  ``None`` --
    never settled within the available fetch -- is a result, not a failure.
    """
    ok = err <= tol
    idx = len(ok)
    while idx > 0 and ok[idx - 1]:
        idx -= 1
    return None if idx == len(ok) else float(fetch[idx])


def adjustment_length(blocks: Sequence[Dict[str, object]],
                      equilibrium: Dict[str, object], preset,
                      tolerance: float = TOLERANCE) -> Dict[str, object]:
    """**The V3 headline.**  How far in before the canopy is in equilibrium.

    Two quantities, both canopy-layer and both taken over the blocks that
    actually contain buildings:

    ``u_canopy``   the mean streamwise velocity below roof height;
    ``uw_at_roof`` the resolved shear stress at the first level above the roofs.

    Two references, reported side by side and neither chosen after the fact:

    ``vs_equilibrium``  against the periodic reference run -- the absolute
                        answer, and the one section 9.4 asks for;
    ``vs_last_row``     against this child's own last row -- the self-referential
                        answer, which is all that is available if the flow never
                        reaches equilibrium, and which is reported so that the
                        two can be told apart.

    Each is reported as a fetch from the **inner edge of the zone** (the domain
    a layout has to spend, which is what section 9.4's standoff argument is
    about) and from the **first building face** (whether the canopy adjustment
    itself is faster or slower behind a standoff).
    """
    rows = [b for b in blocks if b["has_cube"]]
    out: Dict[str, object] = {
        "tolerance": tolerance,
        "n_rows": len(rows),
        "n_blocks": len(blocks),
        "h_m": preset.building_height,
    }
    if not rows:
        out["note"] = "no building rows in the interior"
        return out
    fz = np.array([r["fetch_from_zone_m"] for r in rows], dtype=float)
    fr = np.array([r["fetch_from_first_row_m"] for r in rows], dtype=float)
    for q in ("u_canopy", "uw_at_roof"):
        vals = np.array([np.nan if r[q] is None else r[q] for r in rows], dtype=float)
        eq = equilibrium.get(q)
        entry: Dict[str, object] = {
            "values": vals.tolist(),
            "equilibrium": eq,
            "last_row": float(vals[-1]),
            "first_row": float(vals[0]),
        }
        for ref_name, ref in (("vs_equilibrium", eq), ("vs_last_row", float(vals[-1]))):
            if ref is None or not np.isfinite(ref) or abs(ref) < 1.0e-12:
                entry[ref_name] = {"note": "no usable reference value"}
                continue
            rel = np.abs(vals - ref) / abs(ref)
            zone = _first_settled(fz, rel, tolerance)
            row = _first_settled(fr, rel, tolerance)
            # "Settled at the last row" is not evidence that it stays settled --
            # there is no row after it to disagree.  Flagged rather than
            # silently reported as an adjustment length.
            trivial = zone is not None and len(fz) > 1 and abs(zone - fz[-1]) < 1e-9
            entry[ref_name] = {
                "settled_only_at_the_last_row": bool(trivial),
                "relative_error": rel.tolist(),
                "adjustment_from_zone_m": zone,
                "adjustment_from_zone_h": (None if zone is None
                                           else zone / preset.building_height),
                "adjustment_from_first_row_m": row,
                "adjustment_from_first_row_h": (None if row is None
                                                else row / preset.building_height),
                "error_at_first_row": float(rel[0]),
                "error_at_last_row": float(rel[-1]),
                "max_fetch_from_zone_h": float(fz[-1] / preset.building_height),
            }
        out[q] = entry
    out["fetch_from_zone_m"] = fz.tolist()
    out["fetch_from_first_row_m"] = fr.tolist()
    return out


def ibl_depth(blocks: Sequence[Dict[str, object]], imposed: Dict[str, object],
              preset, ustar: float, tolerance: float = TOLERANCE
              ) -> Dict[str, object]:
    """Internal boundary layer depth against fetch.

    ``delta_i(x)`` is the lowest height above which the block's ``<u>(z)`` is
    within ``tolerance`` of the **imposed** profile -- the parent's own, which is
    what the boundary supplies and what the flow above the IBL still remembers --
    for that height and every height above it.  Below it the flow has been
    reworked by the child's surface; above it, it has not.

    The tolerance is a fraction of the *imposed bulk velocity*, which is the
    conventional IBL criterion and the only scale here that is neither zero at
    the ground nor specific to one height.  Scaling by ``u*`` instead would put
    the bar at 0.02 m/s in a 2.5 m/s flow -- 0.8 %, far below the resolvable
    difference -- and would report no internal layer at all.

    ``None`` at a station means the profile never rejoins the imposed one within
    the domain depth, i.e. the IBL has reached the lid.
    """
    zf = np.asarray(imposed["z"], dtype=float)
    u_imp = np.asarray(imposed["u"], dtype=float)
    u_ref = float(imposed.get("bulk_u") or 0.0)
    if not np.isfinite(u_ref) or u_ref <= 0.0:
        u_ref = float(ustar)
    bar = tolerance * u_ref
    depths, fetches = [], []
    for b in blocks:
        u = np.asarray(b["u"], dtype=float)
        err = np.abs(u - u_imp)
        ok = np.where(np.isnan(err), True, err <= bar)
        idx = len(ok)
        while idx > 0 and ok[idx - 1]:
            idx -= 1
        depths.append(None if idx == len(ok) else float(zf[idx]))
        fetches.append(b["fetch_from_zone_m"])
    return {
        "definition": (f"lowest z above which |<u> - <u>_parent| <= "
                       f"{tolerance:g} U_bulk,imposed = {bar:.4f} m/s at every "
                       "height above"),
        "reference_velocity_m_s": u_ref,
        "fetch_from_zone_m": fetches,
        "delta_i_m": depths,
        "delta_i_over_h": [None if d is None else d / preset.building_height
                           for d in depths],
    }


# --------------------------------------------------------------------------- #
# V4 -- child against the matched-geometry baseline
# --------------------------------------------------------------------------- #


def compare_children(child: GeoStats, base: GeoStats, preset: Preset,
                     fluid: np.ndarray, plane_k: Sequence[int],
                     label_child: str, label_base: str) -> Dict[str, object]:
    """**The V4 headline.**  The mismatched child against the matched one.

    Both children have identical geometry, grid, zone, forcing and schedule --
    the only thing that differs is the layout their *parent* resolved -- so
    every difference below is attributable to the mismatch plus the fact that
    the two were driven by different realisations of the turbulence.  The
    second is bounded, not eliminated: each run's own half-window spread is
    reported alongside every difference, and the two are combined in quadrature
    into the ``spread`` against which the difference should be read.  A
    difference smaller than its spread means V4 has not measured anything, and
    that is exactly the outcome section 10.4 expects ("the interior is
    insensitive to the mismatch beyond the adjustment fetch").
    """
    us = preset.ustar
    ii, jj = interior_indices(preset)
    nk = fluid.shape[2]
    zf = (np.arange(nk) + 0.5) * preset.dz
    zh = zf / preset.building_height

    prof: Dict[str, List[float]] = {"z": zf.tolist(), "z_over_h": zh.tolist()}
    for tag, st in ((label_child, child), (label_base, base)):
        for comp in ("u", "v", "w"):
            prof[f"{comp}_{tag}"] = masked_profile(st.full_mean(comp), fluid,
                                                   ii, jj).tolist()
        prof[f"tke_{tag}"] = masked_profile(st.tke(), fluid, ii, jj).tolist()
        for h in (0, 1):
            prof[f"tke_{tag}_half{h}"] = masked_profile(st.tke(h), fluid,
                                                        ii, jj).tolist()
            prof[f"u_{tag}_half{h}"] = masked_profile(st.mean_half("u", h), fluid,
                                                      ii, jj).tolist()

    def arr(key: str) -> np.ndarray:
        return np.asarray(prof[key], dtype=float)

    tc, tb = arr(f"tke_{label_child}"), arr(f"tke_{label_base}")
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(tb > 0, (tc - tb) / tb, np.nan)
        sc = np.where(tb > 0, np.abs(arr(f"tke_{label_child}_half0")
                                     - arr(f"tke_{label_child}_half1")) / tb, np.nan)
        sb = np.where(tb > 0, np.abs(arr(f"tke_{label_base}_half0")
                                     - arr(f"tke_{label_base}_half1")) / tb, np.nan)
    spread = np.sqrt(sc ** 2 + sb ** 2)

    def band(sel: np.ndarray) -> Dict[str, object]:
        if not np.any(sel) or np.all(np.isnan(rel[sel])):
            return {"n_levels": int(sel.sum()), "mean_relative": None}
        r = float(np.nanmean(rel[sel]))
        sp = float(np.nanmean(spread[sel]))
        med = float(np.nanmedian(spread[sel]))
        with np.errstate(divide="ignore", invalid="ignore"):
            per = np.where(spread[sel] > 0, rel[sel] / spread[sel], np.nan)
        return {
            "n_levels": int(sel.sum()),
            "mean_relative": r,
            "mean_spread": sp,
            "median_spread": med,
            "significance": None if sp <= 0 else r / sp,
            "median_significance": (None if np.all(np.isnan(per))
                                    else float(np.nanmedian(per))),
        }

    du = arr(f"u_{label_child}") - arr(f"u_{label_base}")
    u_floor_c = arr(f"u_{label_child}_half0") - arr(f"u_{label_child}_half1")
    u_floor_b = arr(f"u_{label_base}_half0") - arr(f"u_{label_base}_half1")

    # Criterion A', the design section 0 interior bound taken against the right
    # reference: not the parent (whose buildings are elsewhere) but the child
    # V1 ran, which has the same buildings in the same places.
    crit: Dict[str, object] = {}
    for axis, (lbl, ntot, other) in enumerate([("x", preset.child_itot, jj),
                                               ("y", preset.child_jtot, ii)]):
        e = slab_rms_difference(child.full_mean("u"), base.full_mean("u"),
                                fluid, axis, other) / us
        idx = np.arange(ntot - 1)
        d_lo = (idx + 0.5) * preset.dx
        d_hi = (ntot - idx - 0.5) * preset.dx
        lz = preset.guardwidth + preset.zonewidth
        inside = (d_lo >= lz - 1.0e-9) & (d_hi >= lz - 1.0e-9)
        crit[f"{lbl}_max_interior"] = (float(np.nanmax(e[inside]))
                                       if np.any(inside) else None)
        crit[f"{lbl}_curve"] = {"distance_from_low_face_m": d_lo.tolist(),
                                "error": e.tolist()}
    worst = [v for k, v in crit.items()
             if k.endswith("max_interior") and v is not None]
    crit["max_interior_umean_error_over_ustar"] = max(worst) if worst else None
    crit["threshold"] = 0.05
    crit["passes"] = None if not worst else bool(max(worst) <= 0.05)

    # Resolved-TKE error against distance from each lateral face, between the
    # two children.  The direct analogue of V1's criterion B, with the right
    # reference: if a mismatched parent layout leaves a signature, it should be
    # largest near the boundary -- where the imposed field carries the parent's
    # wakes in the wrong places -- and decay inward.  A flat curve says the
    # difference is not coming from the boundary at all.
    tke_curves: Dict[str, object] = {}
    for axis, (lbl, ntot, other) in enumerate([("x", preset.child_itot, jj),
                                               ("y", preset.child_jtot, ii)]):
        e = slab_rms_difference(child.tke(), base.tke(), fluid, axis, other) / us ** 2
        f = slab_rms_difference(base.tke(0), base.tke(1), fluid, axis,
                                other) / us ** 2 / np.sqrt(2.0)
        idx = np.arange(ntot - 1)
        d_lo = (idx + 0.5) * preset.dx
        d_hi = (ntot - idx - 0.5) * preset.dx
        lz = preset.guardwidth + preset.zonewidth
        inside = (d_lo >= lz - 1.0e-9) & (d_hi >= lz - 1.0e-9)
        tke_curves[lbl] = {
            "distance_from_low_face_m": d_lo.tolist(),
            "distance_from_high_face_m": d_hi.tolist(),
            "error": e.tolist(),
            "baseline_noise_floor": f.tolist(),
            "max_interior": float(np.nanmax(e[inside])) if np.any(inside) else None,
            "error_at_inner_zone_edge": (float(e[inside][0]) if np.any(inside)
                                         else None),
            "error_at_mid_domain": (float(e[inside][e[inside].size // 2])
                                    if np.any(inside) else None),
            "median_noise_floor_interior": (float(np.nanmedian(f[inside]))
                                            if np.any(inside) else None),
        }

    spectra: Dict[str, object] = {}
    for k in plane_k:
        if k not in child.planes or k not in base.planes:
            continue
        kb, eb = analyse.streamwise_spectrum(base.planes[k], ii, jj,
                                             fluid[:, :, k], preset.dx)
        _, ec = analyse.streamwise_spectrum(child.planes[k], ii, jj,
                                            fluid[:, :, k], preset.dx)
        z = (k + 0.5) * preset.dz
        spectra[f"z_{z:g}m"] = {
            "z_m": z, "z_over_h": z / preset.building_height,
            "bands": _band_ratios(kb, eb, ec, preset.dx, ii.size * preset.dx),
        }

    return {
        "labels": {"child": label_child, "baseline": label_base},
        "profiles": prof,
        "tke_difference": {
            "relative": rel.tolist(),
            "combined_spread": spread.tolist(),
            "above": dict(band(zh >= 2.0), z_over_h_min=2.0),
            "canopy": dict(band(zh < 1.0), z_over_h_max=1.0),
        },
        "umean_difference": {
            "difference_over_ustar": (du / us).tolist(),
            "rms_over_ustar": float(np.sqrt(np.nanmean(du ** 2)) / us),
            "noise_floor_over_ustar": float(
                np.sqrt(np.nanmean(u_floor_c ** 2) + np.nanmean(u_floor_b ** 2))
                / us / np.sqrt(2.0)),
        },
        "criterion_a_prime": crit,
        "tke_error_vs_distance": tke_curves,
        "spectra_child_over_baseline": spectra,
    }


# --------------------------------------------------------------------------- #
# Drivers
# --------------------------------------------------------------------------- #


def _levels(times: np.ndarray, offset: float, t0: float, t1: float,
            stride: int) -> List[int]:
    shifted = np.asarray(times, dtype=float) - offset
    idx = np.where((shifted >= t0 - 1.0e-9) & (shifted <= t1 + 1.0e-9))[0]
    return list(idx[::stride])


def accumulate_periodic(casedir: Path, preset, t0: float, t1: float,
                        label: str) -> Tuple[GeoStats, np.ndarray, Dict]:
    """Accumulate a periodic run over its own whole domain."""
    casedir = Path(casedir)
    dump = FieldDump(casedir, preset.parent_expnr, preset.dx)
    shape = (preset.itot - 1, preset.jtot - 1, preset.ktot - 1)
    fluid = load_solid_mask(casedir, shape)
    levels = _levels(dump.times, float(dump.times[0]), t0, t1, preset.stride)
    if len(levels) < 4:
        raise RuntimeError(f"{label}: only {len(levels)} samples in [{t0}, {t1}] s")
    st = accumulate_geo(dump, lambda u, v, w: (u, v, w), levels, fluid, (), label)
    info = {"n_samples": len(levels),
            "window_s": [float(dump.times[levels[0]] - dump.times[0]),
                         float(dump.times[levels[-1]] - dump.times[0])],
            "n_solid_cells": int((~fluid).sum())}
    dump.close()
    return st, fluid, info


def accumulate_child(casedir: Path, preset, plane_k: Sequence[int] = (),
                     label: str = "") -> Tuple[GeoStats, np.ndarray, Dict]:
    """Accumulate a nested child over its **own** fluid mask.

    Deliberately not intersected with the parent's mask: for V3 and V4 the
    parent's buildings are elsewhere, and intersecting would silently reduce the
    child's averaging domain to a shape set by the parent's array.
    """
    casedir = Path(casedir)
    manifest = json.loads((casedir / "manifest.json").read_text())
    dump = FieldDump(casedir, preset.child_expnr, preset.dx)
    shape = (preset.child_itot - 1, preset.child_jtot - 1, preset.child_ktot - 1)
    fluid = load_solid_mask(casedir, shape)
    levels = _levels(dump.times, 0.0, float(manifest["stats_start"]),
                     float(manifest["runtime"]), preset.stride)
    if len(levels) < 4:
        raise RuntimeError(f"{label}: only {len(levels)} samples in the window")
    st = accumulate_geo(dump, lambda u, v, w: (u, v, w), levels, fluid,
                        plane_k, label)
    info = {"n_samples": len(levels),
            "window_s": [float(manifest["stats_start"]), float(manifest["runtime"])],
            "n_solid_cells": int((~fluid).sum())}
    dump.close()
    return st, fluid, info


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #


def write_csv(path: Path, header: Sequence[str],
              columns: Sequence[Sequence[object]]) -> None:
    with Path(path).open("w", encoding="ascii", newline="\n") as fh:
        fh.write(",".join(header) + "\n")
        for row in zip(*columns):
            fh.write(",".join(
                "" if v is None or (isinstance(v, float) and not np.isfinite(v))
                else (f"{v:.9g}" if isinstance(v, float) else str(v))
                for v in row) + "\n")


def write_v3(outdir: Path, metrics: Dict) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "v3_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n",
                                            encoding="ascii")
    for scope in ("core", "canopy"):
        blocks = metrics["blocks"][scope]
        if not blocks:
            continue
        write_csv(outdir / f"streamwise_blocks_{scope}.csv",
                  ["block", "x_lo_m", "fetch_from_zone_m", "fetch_from_zone_h",
                   "fetch_from_first_row_m", "has_cube", "u_canopy", "tke_canopy",
                   "uw_at_roof"],
                  [[b["block"] for b in blocks], [b["x_lo_m"] for b in blocks],
                   [b["fetch_from_zone_m"] for b in blocks],
                   [b["fetch_from_zone_h"] for b in blocks],
                   [b["fetch_from_first_row_m"] for b in blocks],
                   [int(b["has_cube"]) for b in blocks],
                   [b["u_canopy"] for b in blocks],
                   [b["tke_canopy"] for b in blocks],
                   [b["uw_at_roof"] for b in blocks]])
    ibl = metrics.get("ibl", {})
    if ibl.get("delta_i_m"):
        write_csv(outdir / "ibl_depth.csv",
                  ["fetch_from_zone_m", "delta_i_m", "delta_i_over_h"],
                  [ibl["fetch_from_zone_m"], ibl["delta_i_m"], ibl["delta_i_over_h"]])
    eq, imp = metrics.get("equilibrium"), metrics.get("imposed")
    if eq and imp:
        write_csv(outdir / "reference_profiles.csv",
                  ["z", "z_over_h", "u_equilibrium", "u_imposed",
                   "tke_equilibrium", "tke_imposed", "uw_equilibrium", "uw_imposed"],
                  [eq["z"], eq["z_over_h"], eq["u"], imp["u"], eq["tke"], imp["tke"],
                   eq["uw"], imp["uw"]])


def write_v4(outdir: Path, metrics: Dict) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "v4_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n",
                                            encoding="ascii")
    cmp_ = metrics.get("comparison")
    if not cmp_:
        return
    p = cmp_["profiles"]
    c, b = cmp_["labels"]["child"], cmp_["labels"]["baseline"]
    for key, v in cmp_.get("tke_error_vs_distance", {}).items():
        write_csv(outdir / f"tke_error_vs_distance_{key}.csv",
                  ["distance_from_low_face_m", "distance_from_high_face_m",
                   "error", "baseline_noise_floor"],
                  [v["distance_from_low_face_m"], v["distance_from_high_face_m"],
                   v["error"], v["baseline_noise_floor"]])
    write_csv(outdir / "profiles_vs_baseline.csv",
              ["z", "z_over_h", f"u_{c}", f"u_{b}", f"tke_{c}", f"tke_{b}",
               "tke_relative_difference", "combined_spread"],
              [p["z"], p["z_over_h"], p[f"u_{c}"], p[f"u_{b}"], p[f"tke_{c}"],
               p[f"tke_{b}"], cmp_["tke_difference"]["relative"],
               cmp_["tke_difference"]["combined_spread"]])


# --------------------------------------------------------------------------- #
# Human-readable summaries
# --------------------------------------------------------------------------- #


def summary_v3(m: Dict) -> str:
    cfg = m["configuration"]
    lines = [
        f"V3 {m['child_key']}: standoff {cfg['standoff_cells']} cells "
        f"({cfg['standoff_m']:g} m; {cfg['first_row_fetch_m']:g} m from the ramp), "
        f"{cfg['n_rows']} rows, interior fetch {cfg['streamwise_interior_h']:.1f} h",
        f"  samples {m['samples']['child']}, window {m['samples']['window_s']} s",
    ]
    eq, imp = m.get("equilibrium"), m.get("imposed")
    if eq and imp:
        lines.append(
            f"  imposed vs equilibrium: bulk {imp['bulk_u']:.3f} vs "
            f"{eq['bulk_u']:.3f} m/s ({100 * (imp['bulk_u'] / eq['bulk_u'] - 1):+.1f} %), "
            f"canopy <u> {imp['u_canopy']:.3f} vs {eq['u_canopy']:.3f} m/s "
            f"({100 * (imp['u_canopy'] / eq['u_canopy'] - 1):+.1f} %), "
            f"<u'w'> at roof {imp['uw_at_roof']:+.4f} vs {eq['uw_at_roof']:+.4f} m2/s2")
    adj = m.get("adjustment", {})
    for q in ("u_canopy", "uw_at_roof"):
        e = adj.get(q)
        if not isinstance(e, dict):
            continue
        for ref in ("vs_equilibrium", "vs_last_row"):
            r = e.get(ref)
            if not isinstance(r, dict) or "adjustment_from_zone_m" not in r:
                continue
            z = r["adjustment_from_zone_h"]
            f = r["adjustment_from_first_row_h"]
            lines.append(
                f"  {q:11s} {ref:15s} adjusted after "
                + ("NOT REACHED" if z is None else f"{z:5.2f} h from the zone")
                + ", "
                + ("NOT REACHED" if f is None else f"{f:5.2f} h from row 1")
                + f"; error row 1 {100 * r['error_at_first_row']:.1f} %, "
                f"last row {100 * r['error_at_last_row']:.1f} %, "
                f"max fetch {r['max_fetch_from_zone_h']:.1f} h")
    ibl = m.get("ibl", {})
    d = [x for x in ibl.get("delta_i_over_h", []) if x is not None]
    if d:
        lines.append(f"  IBL depth {d[0]:.2f} h at the first block to "
                     f"{d[-1]:.2f} h at the last")
    return "\n".join(lines)


def summary_v4(m: Dict) -> str:
    c = m.get("comparison")
    if not c:
        return "V4: no comparison (no baseline)"
    d = c["tke_difference"]
    a = c["criterion_a_prime"]
    u = c["umean_difference"]

    def pct(x):
        return "  --  " if x is None else f"{100 * x:+6.2f}%"

    lines = [
        f"V4 {m['child_key']}: {c['labels']['child']} against the matched-geometry "
        f"baseline {c['labels']['baseline']}",
        f"  samples child {m['samples']['child']}, baseline {m['samples']['baseline']}",
        f"  resolved-TKE difference above z/h = 2: {pct(d['above']['mean_relative'])} "
        f"against a combined half-window spread of "
        f"{pct(d['above'].get('mean_spread')).strip('+')} mean / "
        f"{pct(d['above'].get('median_spread')).strip('+')} median"
        + ("" if d['above'].get('significance') is None
           else f"  ({d['above']['significance']:+.1f} sigma, "
                f"{d['above']['median_significance']:+.1f} per-height median)"),
        f"  inside the canopy (z/h < 1):            "
        f"{pct(d['canopy']['mean_relative'])} against "
        f"{pct(d['canopy'].get('mean_spread')).strip('+')}",
        f"  <u> rms difference {u['rms_over_ustar']:.4f} u* "
        f"(sampling floor {u['noise_floor_over_ustar']:.4f} u*)",
        "  resolved-TKE error vs distance (u*^2): "
        + ", ".join(
            f"{k} edge {('--' if v['error_at_inner_zone_edge'] is None else format(v['error_at_inner_zone_edge'], '.3f'))}"
            f" mid {('--' if v['error_at_mid_domain'] is None else format(v['error_at_mid_domain'], '.3f'))}"
            f" max {('--' if v['max_interior'] is None else format(v['max_interior'], '.3f'))}"
            f" floor {('--' if v['median_noise_floor_interior'] is None else format(v['median_noise_floor_interior'], '.3f'))}"
            for k, v in c["tke_error_vs_distance"].items()),
        "  criterion A' (mean flow, interior, against the V1 child): "
        + ("--" if a["max_interior_umean_error_over_ustar"] is None
           else f"{a['max_interior_umean_error_over_ustar']:.4f}")
        + f" u* against {a['threshold']:g} -- "
        + ("PASS" if a["passes"] else "FAIL"),
    ]
    for name, sp in c["spectra_child_over_baseline"].items():
        b = sp["bands"]

        def g(key):
            v = b[key]["mean_of_ratios"]
            return "   --  " if v is None else f"{v:6.3f}"
        lines.append(f"  spectra {name:12s} z/h = {sp['z_over_h']:.2f}   "
                     f"16-64 m {g('band_16_64m')}   8-16 m {g('band_8_16m')}   "
                     f"> L/4 {g('lambda_gt_quarter_L')}   < 4dx {g('lambda_lt_4dx')}")
    return "\n".join(lines)


if __name__ == "__main__":  # pragma: no cover - run_geometry.py is the entry point
    raise SystemExit("this module is driven by run_geometry.py")
