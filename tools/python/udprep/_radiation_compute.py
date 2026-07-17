"""Pure numeric helpers for the radiation section (no IO / case state).

Extracted so the shortwave/SEB compute can be tested and reused apart from
the IO/subprocess/netCDF orchestration in RadiationSection.
"""
from __future__ import annotations

import numpy as np


def interp_makima(x: np.ndarray, y: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x_new = np.asarray(x_new, dtype=float)

    if x.size < 2:
        raise ValueError("Need at least two points for interpolation")
    if np.any(np.diff(x) <= 0):
        raise ValueError("x must be strictly increasing for interpolation")

    if x.size == 2:
        # Modified Akima needs >=3 points for its slope estimates; with exactly
        # two it degenerates to linear interpolation (extrapolated on both sides).
        slope = (y[1] - y[0]) / (x[1] - x[0])
        return y[0] + slope * (x_new - x[0])

    m = np.diff(y) / np.diff(x)
    m1 = 2.0 * m[0] - m[1]
    m2 = 2.0 * m1 - m[0]
    m_end1 = 2.0 * m[-1] - m[-2]
    m_end2 = 2.0 * m_end1 - m[-1]
    m_ext = np.concatenate(([m2, m1], m, [m_end1, m_end2]))

    d = np.zeros_like(x)
    for i in range(x.size):
        w1 = abs(m_ext[i + 3] - m_ext[i + 2]) + abs(m_ext[i + 3] + m_ext[i + 2]) * 0.5
        w2 = abs(m_ext[i + 1] - m_ext[i]) + abs(m_ext[i + 1] + m_ext[i]) * 0.5
        if w1 + w2 > 0.0:
            d[i] = (w1 * m_ext[i + 1] + w2 * m_ext[i + 2]) / (w1 + w2)
        else:
            d[i] = 0.5 * (m_ext[i + 1] + m_ext[i + 2])

    idx = np.searchsorted(x, x_new, side="right") - 1
    idx = np.clip(idx, 0, x.size - 2)
    x0 = x[idx]
    x1 = x[idx + 1]
    y0 = y[idx]
    y1 = y[idx + 1]
    d0 = d[idx]
    d1 = d[idx + 1]
    h = x1 - x0
    t = (x_new - x0) / h

    t2 = t * t
    t3 = t2 * t
    h00 = 2.0 * t3 - 3.0 * t2 + 1.0
    h10 = t3 - 2.0 * t2 + t
    h01 = -2.0 * t3 + 3.0 * t2
    h11 = t3 - t2

    return h00 * y0 + h10 * h * d0 + h01 * y1 + h11 * h * d1

def net_shortwave_reflections(
    sdir: np.ndarray,
    dsky: float,
    vf,
    svf: np.ndarray,
    albedo: np.ndarray,
    tol: float = 0.01,
    max_iter: int = 1000,
) -> np.ndarray:
    """
    Compute net shortwave including reflections (tools/SEB/netShortwave.m); pure compute..

    Parameters
    ----------
    sdir : np.ndarray
        Direct shortwave on facets [W/m^2].
    dsky : float
        Diffuse sky irradiance [W/m^2].
    vf : array-like or sparse matrix
        View factor matrix between facets.
    svf : np.ndarray
        Sky view factor per facet.
    albedo : np.ndarray
        Facet albedo (0-1).
    tol : float
        Convergence threshold for additional absorbed energy.
    max_iter : int
        Safety cap on the iteration count.
    """
    sdir = np.asarray(sdir, dtype=float)
    svf = np.asarray(svf, dtype=float)
    albedo = np.asarray(albedo, dtype=float)
    if sdir.shape != svf.shape or sdir.shape != albedo.shape:
        raise ValueError("sdir, svf, and albedo must have matching shapes")
    if tol <= 0.0:
        raise ValueError("tol must be > 0")
    if max_iter <= 0:
        raise ValueError("max_iter must be > 0")

    kin0 = sdir + dsky * svf
    knet = (1.0 - albedo) * kin0
    kout = albedo * kin0

    for _ in range(max_iter):
        vf_kout = vf @ kout
        kadd = (1.0 - albedo) * vf_kout
        kout = albedo * vf_kout
        knet = knet + kadd

        denom = np.maximum(knet - kadd, 1.0e-12)
        if np.max(kadd / denom) < tol:
            break

    return knet


def net_shortwave_nonscattering(
    sdir: np.ndarray, dsky: float, fss: np.ndarray, albedo: np.ndarray
) -> np.ndarray:
    """Net shortwave without reflections: ``(1 - albedo) * (sdir + dsky*fss)``."""
    return (1.0 - np.asarray(albedo)) * (np.asarray(sdir) + dsky * np.asarray(fss))
