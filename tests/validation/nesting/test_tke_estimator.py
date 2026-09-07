#!/usr/bin/env python3
"""Regression test for the driving-parent TKE estimator (nesting review R1).

``make_child_case.accumulate_profile`` used to form the driving parent's
resolved TKE from the combined space-and-time mean -- ``0.5 * sum(<q^2> -
<q>^2)`` over horizontally-averaged, then time-summed, moments.  That is a
real quantity (it includes the DISPERSIVE variance of the time-mean field, a
canopy's fixed spatial pattern), but it is not the same statistic
``analyse.Bundle.tke`` reports for the fine truth: ``Bundle`` subtracts each
cell's OWN temporal mean before averaging horizontally, so it is pure
temporal turbulence with the dispersive part excluded.  ``analyse_v0.py``
compared the two as if they were one quantity.

``make_child_case.decompose_tke`` is the fix: from the same accumulated
horizontal moments plus a new per-CELL temporal accumulation, it reports all
three quantities -- ``temporal`` (analyse.Bundle-consistent), ``dispersive``
(the part the old formula silently included) and ``total`` (their sum,
numerically identical to the old formula).  This test pins that function on
two independent problems:

1. the reviewer's own two-cell stationary example (1 and 3 m/s, no time
   variation at all): the old formula gives 0.5 m^2/s^2, the corrected
   temporal estimator gives 0.0.
2. a field with a KNOWN, separable dispersive part (a static horizontal
   pattern) and a KNOWN temporal part (a spatially uniform, time-varying
   signal), built directly from per-level fields the way
   ``accumulate_profile`` consumes them, so both the accumulation and the
   reduction are exercised, not just the reduction formula in isolation.

No solver, no case directory, no I/O: this is pure array arithmetic, so it
runs in a fraction of a second and needs nothing beyond ``numpy``.

Run directly, or through ``tests/run_tests.py nesting-validation``:

    python tests/validation/nesting/test_tke_estimator.py
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import numpy as np

from make_child_case import decompose_tke


def _accumulate(levels, fluid):
    """Reproduce ``accumulate_profile``'s reduction over a list of per-level
    ``{"u": arr2d, "v": arr2d, "w": arr2d}`` dicts, ``arr2d`` shape (ni, nj).

    Mirrors the real function's arithmetic exactly (mask, horizontal sum,
    divide by the fluid-cell count) without any of the parent-file, grid or
    case-directory machinery around it -- those are exercised by
    ``test_v0_tiny.py`` instead.
    """
    ni, nj = fluid.shape
    n_levels = len(levels)
    nkp = 1
    imean = {c: np.zeros(nkp) for c in "uvw"}
    imsq = {c: np.zeros(nkp) for c in "uvw"}
    cell_sum = {c: np.zeros((ni, nj, nkp)) for c in "uvw"}
    cell_sumsq = {c: np.zeros((ni, nj, nkp)) for c in "uvw"}
    icell_col = np.array([fluid.sum()], dtype=float)
    for level in levels:
        for c in "uvw":
            blk = level[c].reshape(ni, nj, 1)
            masked = np.where(fluid[:, :, None], blk, 0.0)
            imean[c] += masked.sum(axis=(0, 1)) / icell_col
            imsq[c] += np.where(fluid[:, :, None], blk * blk, 0.0).sum(axis=(0, 1)) / icell_col
            cell_sum[c] += masked
            cell_sumsq[c] += np.where(fluid[:, :, None], blk * blk, 0.0)
    return imean, imsq, cell_sum, cell_sumsq, icell_col, n_levels


class TestTkeEstimator(unittest.TestCase):
    """Pins ``decompose_tke`` against closed-form expected values."""

    def test_stationary_two_cell_field_matches_the_reviewers_example(self):
        """No time variation at all: old formula 0.5, corrected temporal 0.0."""
        fluid = np.ones((2, 1), dtype=bool)
        level = {"u": np.array([[1.0], [3.0]]), "v": np.zeros((2, 1)), "w": np.zeros((2, 1))}
        imean, imsq, cell_sum, cell_sumsq, icell, n = _accumulate([level], fluid)

        parts = decompose_tke(imean, imsq, cell_sum, cell_sumsq, float(n), icell)

        self.assertAlmostEqual(float(parts["total"][0]), 0.5, places=12)
        self.assertAlmostEqual(float(parts["temporal"][0]), 0.0, places=12)
        self.assertAlmostEqual(float(parts["dispersive"][0]), 0.5, places=12)
        # the identity the docstring promises: temporal + dispersive == total,
        # to floating-point round-off, by construction
        self.assertAlmostEqual(
            float(parts["temporal"][0] + parts["dispersive"][0]),
            float(parts["total"][0]), places=12)

    def test_separable_dispersive_and_temporal_parts_are_each_recovered(self):
        """A static spatial pattern plus a time-varying uniform signal.

        u(i, j, t) = pattern(i, j) + f(t), with ``pattern`` having a known
        spatial variance and ``f`` a known temporal variance, and both having
        zero mean (so their cross term in the combined space-time variance
        vanishes exactly, which is what makes ``total == dispersive +
        temporal`` an identity rather than an approximation).  v and w are
        held at zero throughout, so the whole TKE is carried by u alone and
        every expected number is a plain half-variance.
        """
        pattern = np.array([[-3.0, 1.0, 2.0],
                            [4.0, -1.0, -3.0]])  # 2 x 3, mean 0
        pattern -= pattern.mean()
        f = np.array([0.0, 2.0, -2.0, 4.0, -4.0])  # mean 0, 5 samples
        fluid = np.ones(pattern.shape, dtype=bool)

        var_spatial = float(np.mean(pattern ** 2))
        var_temporal = float(np.mean(f ** 2))

        levels = [{"u": pattern + ft, "v": np.zeros_like(pattern),
                   "w": np.zeros_like(pattern)} for ft in f]
        imean, imsq, cell_sum, cell_sumsq, icell, n = _accumulate(levels, fluid)

        parts = decompose_tke(imean, imsq, cell_sum, cell_sumsq, float(n), icell)

        expected_temporal = 0.5 * var_temporal
        expected_dispersive = 0.5 * var_spatial
        expected_total = expected_temporal + expected_dispersive

        self.assertAlmostEqual(float(parts["temporal"][0]), expected_temporal, places=12)
        self.assertAlmostEqual(float(parts["dispersive"][0]), expected_dispersive, places=12)
        self.assertAlmostEqual(float(parts["total"][0]), expected_total, places=12)
        # and the corrected estimator must differ from the old, inflated one
        # whenever there is any dispersive part at all
        self.assertGreater(abs(parts["total"][0] - parts["temporal"][0]), 1.0e-9)

    def test_zero_dispersive_part_recovers_the_old_formula_exactly(self):
        """A spatially uniform field: dispersive is 0, temporal == total."""
        f = np.array([1.0, -1.0, 3.0, -3.0])
        fluid = np.ones((2, 2), dtype=bool)
        levels = [{"u": np.full((2, 2), ft), "v": np.zeros((2, 2)),
                   "w": np.zeros((2, 2))} for ft in f]
        imean, imsq, cell_sum, cell_sumsq, icell, n = _accumulate(levels, fluid)

        parts = decompose_tke(imean, imsq, cell_sum, cell_sumsq, float(n), icell)

        self.assertAlmostEqual(float(parts["dispersive"][0]), 0.0, places=12)
        self.assertAlmostEqual(float(parts["temporal"][0]), float(parts["total"][0]),
                               places=12)
        self.assertAlmostEqual(float(parts["temporal"][0]), 0.5 * float(np.mean(f ** 2)),
                               places=12)


if __name__ == "__main__":
    unittest.main()
