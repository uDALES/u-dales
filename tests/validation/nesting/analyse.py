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
# Driver
# --------------------------------------------------------------------------- #


def _levels_in_window(times: np.ndarray, offset: float, t0: float, t1: float,
                      stride: int) -> List[int]:
    shifted = times - offset
    idx = np.where((shifted >= t0 - 1.0e-9) & (shifted <= t1 + 1.0e-9))[0]
    return list(idx[::stride])


def run(parent_dir: Path, child_dir: Path, outdir: Path, preset: Preset,
        make_plots: bool = True) -> Dict[str, object]:
    parent_dir, child_dir, outdir = Path(parent_dir), Path(child_dir), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((child_dir / "manifest.json").read_text())
    t_offset = float(manifest["t_offset"])
    t0 = float(manifest["stats_start"])
    t1 = float(manifest["runtime"])

    ni, nj, nk = preset.child_itot - 1, preset.child_jtot - 1, preset.child_ktot - 1
    mask = load_solid_mask(child_dir, (ni, nj, nk))
    ii, jj = interior_indices(preset)
    kk = [int(round(z / preset.dz - 0.5)) for z in preset.spectra_heights]
    kk = [k for k in kk if 0 <= k < nk]

    pdump = FieldDump(parent_dir, preset.parent_expnr, preset.dx)
    cdump = FieldDump(child_dir, preset.child_expnr, preset.dx)
    plev = _levels_in_window(pdump.times, t_offset, t0, t1, preset.stride)
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

    parent = accumulate(pdump, cut, mask, plev, kk, (ii, jj), "parent sub-region")
    child = accumulate(cdump, whole, mask, clev, kk, (ii, jj), "nested child")

    metrics = _compare(parent, child, preset, mask, ii, jj, kk, manifest,
                       len(plev), len(clev))
    _write_outputs(outdir, parent, child, preset, mask, ii, jj, kk, metrics,
                   make_plots)
    return metrics


def _compare(parent: Bundle, child: Bundle, preset: Preset, mask: np.ndarray,
             ii: np.ndarray, jj: np.ndarray, kk: Sequence[int],
             manifest: Dict, n_parent: int, n_child: int) -> Dict[str, object]:
    us = preset.ustar
    zf = (np.arange(preset.child_ktot - 1) + 0.5) * preset.dz
    lz = preset.guardwidth + preset.zonewidth

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
    spectra = {}
    for k in kk:
        kp, pw = streamwise_spectrum(parent.planes[k], ii, jj, mask[:, :, k], preset.dx)
        kc, pc = streamwise_spectrum(child.planes[k], ii, jj, mask[:, :, k], preset.dx)
        z = (k + 0.5) * preset.dz
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(pw > 0, pc / pw, np.nan)
        spectra[f"z_{z:g}m"] = {
            "z_m": z,
            "z_over_h": z / preset.building_height,
            "wavenumber_rad_per_m": kp.tolist(),
            "E_parent": pw.tolist(),
            "E_child": pc.tolist(),
            "child_over_parent": ratio.tolist(),
            "band_mean_ratio_resolved": float(np.nanmean(ratio[1:max(2, len(ratio) // 2)])),
        }

    return {
        "preset": preset.name,
        "samples": {"parent": n_parent, "child": n_child},
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
                   kk: Sequence[int], metrics: Dict, make_plots: bool) -> None:
    (outdir / "v1_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n",
                                            encoding="ascii")

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
        lines.append(f"    {name:14s} z/h = {s['z_over_h']:.2f}   "
                     f"ratio = {s['band_mean_ratio_resolved']:.3f}")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("parent_dir", type=Path)
    parser.add_argument("child_dir", type=Path)
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--preset", default="production")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()

    preset = get_preset(args.preset)
    metrics = run(args.parent_dir, args.child_dir, args.outdir, preset,
                  make_plots=not args.no_plots)
    print(summary(metrics))
    print(f"\nwritten to {args.outdir}")


if __name__ == "__main__":
    main()
