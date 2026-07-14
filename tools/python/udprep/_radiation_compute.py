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

