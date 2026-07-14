"""Pure facet-data helpers for uDALES postprocessing.

Array-in/array-out functions with no case state, extracted from UDBase so facet
maths is separable and testable. UDBase keeps thin methods that gather case state
(grid sizes, spacings, facet sections) and delegate here.
"""
from __future__ import annotations

from typing import Dict, Tuple

import numpy as np


def facsec_to_field(
    var: np.ndarray,
    facsec: Dict[str, np.ndarray],
    dz: np.ndarray,
    shape: Tuple[int, int, int],
    dx: float,
    dy: float,
) -> np.ndarray:
    """Distribute a per-facet variable into a 3-D density field.

    Each facet section contributes ``var[facid] * area / (dx*dy*dz[k])`` to the
    grid cell it sits in, giving a density ``[var_units / m]``.

    Parameters
    ----------
    var : ndarray, shape (n_facets,)
    facsec : dict with keys ``facid`` (int per section), ``area`` (per section),
        and ``locs`` ((n_sections, 3) 0-based (i, j, k) cell indices).
    dz : ndarray, vertical grid spacing per k.
    shape : (itot, jtot, ktot).
    dx, dy : horizontal grid spacings.

    Returns
    -------
    ndarray, shape ``shape`` (float32).
    """
    itot, jtot, ktot = shape
    fld = np.zeros((itot, jtot, ktot), dtype=np.float32)

    facids = facsec["facid"]
    areas = facsec["area"]
    locs = facsec["locs"]  # (i, j, k) locations (0-based)

    i_idx = locs[:, 0].astype(int)
    j_idx = locs[:, 1].astype(int)
    k_idx = locs[:, 2].astype(int)

    for m in range(len(areas)):
        facid = facids[m]
        i, j, k = i_idx[m], j_idx[m], k_idx[m]
        cell_volume = dx * dy * dz[k]
        fld[i, j, k] += var[facid] * areas[m] / cell_volume

    return fld
