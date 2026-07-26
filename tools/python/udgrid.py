"""Pure grid-coordinate helpers for uDALES.

Array/scalar-in, array-out functions with no case state, extracted from UDBase so
the grid maths is separable and testable. UDBase does the file IO/fallback and
assigns the returned coordinate arrays.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np


def horizontal_grid(
    itot: int, dx: float, jtot: int, dy: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(xm, xt, ym, yt)``: cell edges and centres in x and y."""
    xm = np.arange(itot) * dx
    xt = xm + 0.5 * dx
    ym = np.arange(jtot) * dy
    yt = ym + 0.5 * dy
    return xm, xt, ym, yt


def z_grid_from_profile(zt: np.ndarray, zsize: float) -> Tuple[np.ndarray, np.ndarray]:
    """Derive ``(zm, dzt)`` (cell edges and spacings) from cell centres ``zt``."""
    zt = np.asarray(zt, dtype=float)
    zm_full = np.concatenate([[0], 0.5 * (zt[:-1] + zt[1:]), [zsize]])
    zm = zm_full[:-1]  # match dimensions
    dzt = np.diff(np.append(zm, zsize))
    return zm, dzt


def uniform_z_grid(
    zsize: float, ktot: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(zm, zt, dzt)`` for an equidistant z-grid of ``ktot`` cells."""
    dz = zsize / ktot
    zm = np.arange(ktot) * dz
    zt = zm + 0.5 * dz
    dzt = np.full(ktot, dz)
    return zm, zt, dzt
