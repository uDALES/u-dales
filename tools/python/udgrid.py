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
    """Derive ``(zm, dzt)`` (cell edges and spacings) from cell centres ``zt``.

    The solver rebuilds the faces from the ``prof.inp`` centres with the
    recursion ``zh(1) = 0``, ``zh(k+1) = 2*zf(k) - zh(k)`` (modglobal), which
    makes every centre the exact midpoint of its own faces. Midpoints of
    *adjacent centres* are not the faces on a stretched grid, so the same
    recursion is used here. ``zsize`` is only checked for consistency against
    the recovered domain top.
    """
    zt = np.asarray(zt, dtype=float)
    faces = np.empty(zt.size + 1, dtype=float)
    faces[0] = 0.0
    for k, centre in enumerate(zt):
        faces[k + 1] = 2.0 * centre - faces[k]
    if np.any(np.diff(faces) <= 0):
        import warnings

        warnings.warn(
            "profile centres do not describe a valid solver grid: the face "
            "recursion yields non-increasing faces. The solver would derive "
            "the same broken spacings from this prof.inp.",
            RuntimeWarning,
            stacklevel=2,
        )
    if zsize > 0 and abs(faces[-1] - zsize) > 1e-6 * zsize:
        import warnings

        warnings.warn(
            f"z grid from profile tops out at {faces[-1]:.6g} m but zsize is "
            f"{zsize:.6g} m; the profile centres and the namelist disagree.",
            RuntimeWarning,
            stacklevel=2,
        )
    zm = faces[:-1]
    dzt = np.diff(faces)
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
