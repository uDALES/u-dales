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
    include_mask: np.ndarray | None = None,
) -> np.ndarray:
    """Distribute a per-facet variable into a 3-D density field.

    Each facet section contributes ``var[facid] * area / (dx*dy*dz[k])`` to the
    grid cell it sits in, giving a density ``[var_units / m]``. Vectorised with
    ``np.add.at`` and accumulated in float64 (returned as float32), replacing the
    former per-section Python loop.

    Parameters
    ----------
    var : ndarray, shape (n_facets,)
    facsec : dict with keys ``facid`` (int per section), ``area`` (per section),
        and ``locs`` ((n_sections, 3) 0-based (i, j, k) cell indices).
    dz : ndarray, vertical grid spacing per k.
    shape : (itot, jtot, ktot).
    dx, dy : horizontal grid spacings.
    include_mask : ndarray of bool, shape (n_facets,), optional
        Per-facet inclusion mask; sections whose ``facid`` is masked out are
        skipped (used for building-id filtering).

    Returns
    -------
    ndarray, shape ``shape`` (float32).
    """
    itot, jtot, ktot = shape

    facids = np.asarray(facsec["facid"]).astype(int)
    areas = np.asarray(facsec["area"], dtype=float)
    locs = np.asarray(facsec["locs"])
    i_idx = locs[:, 0].astype(int)
    j_idx = locs[:, 1].astype(int)
    k_idx = locs[:, 2].astype(int)

    contrib = np.asarray(var, dtype=float)[facids] * areas / (dx * dy * np.asarray(dz, dtype=float)[k_idx])

    if include_mask is not None:
        keep = np.asarray(include_mask, dtype=bool)[facids]
        i_idx, j_idx, k_idx, contrib = i_idx[keep], j_idx[keep], k_idx[keep], contrib[keep]

    fld = np.zeros((itot, jtot, ktot), dtype=np.float64)
    np.add.at(fld, (i_idx, j_idx, k_idx), contrib)
    return fld.astype(np.float32)
