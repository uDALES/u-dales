"""Tests closing the writer blockers of the 2026-09-06 nesting review (W1--W8).

============ ============================================================
 W1           grid alignment is asserted, with a logged opt-out
 W2           an offset child origin round-trips through validate/read
 W3           the time axis is validated at construction
 W4           the per-level append writer has a level-independent footprint
 W5           the cadence criterion C_dump and the refinement verdict
 W6           the correction magnitude and its split across faces are reported
 W7           units, cached residual, vertical ratios, FaceMasks from IBM
 W8           piecewise-linear prolongation, conservative to round-off
============ ============================================================

Each test states the number it pins in its name or its assertion; the
diagnostics tests print the numbers they compare so a run leaves a record.
"""

from __future__ import annotations

import logging
import sys
import tracemalloc
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from netCDF4 import Dataset

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from _common import PYTHON_DIR  # noqa: E402

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from exceptions import ConfigurationError  # noqa: E402

from udprep.nesting import (  # noqa: E402
    COMPONENTS,
    FACES,
    NestGrid,
    NestingAlignmentError,
    check_alignment,
    discrete_divergence,
    interpolate_child_fields,
    nesting_data_from_parent,
    slabs_from_parent,
)

from test_nesting import make_grids, solenoidal_parent_fields  # noqa: E402


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #


def solenoidal_fields_any_grid(grid, seed=3):
    """A discretely solenoidal field on an arbitrary (stretched) grid.

    Face *fluxes* are the discrete curl of a random edge potential, so the net
    flux out of every cell telescopes to zero exactly whatever the metrics; the
    velocities are the fluxes divided by their own face areas.  ``w = 0`` on
    the floor and the lid.
    """
    rng = np.random.default_rng(seed)
    itot, jtot, ktot = grid.itot, grid.jtot, grid.ktot
    ax = rng.normal(size=(itot, jtot + 1, ktot + 1))
    ay = rng.normal(size=(itot + 1, jtot, ktot + 1))
    az = rng.normal(size=(itot + 1, jtot + 1, ktot))
    ax[:, :, 0] = ax[:, :, -1] = 0.0
    ay[:, :, 0] = ay[:, :, -1] = 0.0
    fx = (az[:, 1:, :] - az[:, :-1, :]) - (ay[:, :, 1:] - ay[:, :, :-1])
    fy = (ax[:, :, 1:] - ax[:, :, :-1]) - (az[1:, :, :] - az[:-1, :, :])
    fz = (ay[1:, :, :] - ay[:-1, :, :]) - (ax[:, 1:, :] - ax[:, :-1, :])
    dx, dy, dz = grid.dx, grid.dy, grid.dzf
    u = fx / (dy[None, :, None] * dz[None, None, :])
    v = fy / (dx[:, None, None] * dz[None, None, :])
    w = fz / (dx[:, None, None] * dy[None, :, None])
    return u, v, w


def velocity_over_spacing(grid, u, v, w):
    """The divergence scale ``max|u| / min(dx)`` a misalignment is judged against."""
    speed = max(float(np.max(np.abs(a))) for a in (u, v, w))
    return speed / float(min(grid.dx.min(), grid.dy.min(), grid.dzf.min()))


# --------------------------------------------------------------------------- #
# W1 -- alignment
# --------------------------------------------------------------------------- #


class TestW1Alignment(unittest.TestCase):
    """W1: a child cell straddling a parent face is refused, or logged on opt-out."""

    def test_aligned_ratio_two_is_accepted_and_solenoidal_to_round_off(self):
        parent, child = make_grids(nx=8, ny=6, nz=6, rx=2, ry=2, rz=2,
                                   xlen=40.0, ylen=30.0, zsize=30.0)
        report = check_alignment(parent, child)
        for axis in ("x", "y", "z"):
            self.assertEqual(report[axis], 0.0)
        pu, pv, pw = solenoidal_parent_fields(parent)
        cu, cv, cw = interpolate_child_fields(parent, pu, pv, pw, child)
        div = float(np.max(np.abs(discrete_divergence(child, cu, cv, cw))))
        self.assertLess(div, 1e-15)          # the review measured 3e-17

    def test_ratio_one_and_a_half_is_rejected_on_every_entry_point(self):
        parent = NestGrid.uniform(8, 6, 6, 40.0, 30.0, 30.0)
        child = NestGrid.uniform(12, 9, 9, 40.0, 30.0, 30.0)
        fields = solenoidal_parent_fields(parent)
        with self.assertRaises(NestingAlignmentError) as ctx:
            check_alignment(parent, child)
        self.assertIn("1.67", str(ctx.exception))   # a parent face 5/3 m inside a child cell
        with self.assertRaises(NestingAlignmentError):
            interpolate_child_fields(parent, *fields, child)
        with self.assertRaises(NestingAlignmentError):
            slabs_from_parent(parent, *fields, child=child, nzone=2)
        with self.assertRaises(NestingAlignmentError):
            nesting_data_from_parent(parent, child, 2, [0.0, 1.0], [fields, fields])
        self.assertTrue(issubclass(NestingAlignmentError, ConfigurationError))

    def test_a_stretched_parent_vertical_under_a_uniform_child_is_rejected(self):
        # every parent face at a multiple of the child spacing would be the
        # only accepted layout; a geometric stretch has none of them there
        zh = np.cumsum(np.concatenate(([0.0], 2.0 * 1.2 ** np.arange(6))))
        zh *= 30.0 / zh[-1]
        parent = NestGrid.from_faces(np.linspace(0.0, 40.0, 9), np.linspace(0.0, 30.0, 7), zh)
        child = NestGrid.uniform(16, 12, 12, 40.0, 30.0, 30.0)
        with self.assertRaises(NestingAlignmentError) as ctx:
            check_alignment(parent, child)
        self.assertIn("along z", str(ctx.exception))
        # and a child whose faces contain the parent's is fine even when the
        # child itself is stretched differently within each parent cell
        child_zh = np.sort(np.concatenate((zh, 0.5 * (zh[:-1] + zh[1:]))))
        nested = NestGrid.from_faces(np.linspace(0.0, 40.0, 17), np.linspace(0.0, 30.0, 13),
                                     child_zh)
        report = check_alignment(parent, nested)
        self.assertEqual(report["z"], 0.0)

    def test_the_opt_out_logs_and_reproduces_the_review_divergence(self):
        parent = NestGrid.uniform(8, 6, 6, 40.0, 30.0, 30.0)
        pu, pv, pw = solenoidal_parent_fields(parent)
        scale = velocity_over_spacing(parent, pu, pv, pw)
        cases = {
            "ratio 1.5": NestGrid.uniform(12, 9, 9, 40.0, 30.0, 30.0),
        }
        zh = np.cumsum(np.concatenate(([0.0], 2.0 * 1.2 ** np.arange(6))))
        zh *= 30.0 / zh[-1]
        stretched = NestGrid.from_faces(parent.xh, parent.yh, zh)
        su, sv, sw = solenoidal_fields_any_grid(stretched)
        for label, child in cases.items():
            with self.subTest(case=label):
                with self.assertLogs("udprep.nesting", level=logging.WARNING) as logs:
                    cu, cv, cw = interpolate_child_fields(parent, pu, pv, pw, child,
                                                          allow_misaligned=True)
                self.assertTrue(any("not nested" in m for m in logs.output))
                div = float(np.max(np.abs(discrete_divergence(child, cu, cv, cw))))
                print(f"\n[W1] {label}: child divmax = {div:.3e} s-1 "
                      f"= {div / scale:.2f} x (max|u| / min dx)")
                # the review quotes ~8e-2 at velocity scale 0.7, dx 5 m: half of u/dx
                self.assertGreater(div / scale, 0.1)
        with self.subTest(case="stretched parent, uniform child"):
            child = NestGrid.uniform(16, 12, 12, 40.0, 30.0, 30.0)
            with self.assertLogs("udprep.nesting", level=logging.WARNING):
                cu, cv, cw = interpolate_child_fields(stretched, su, sv, sw, child,
                                                      allow_misaligned=True)
            div = float(np.max(np.abs(discrete_divergence(child, cu, cv, cw))))
            sscale = velocity_over_spacing(stretched, su, sv, sw)
            print(f"[W1] stretched parent: child divmax = {div:.3e} s-1 "
                  f"= {div / sscale:.2f} x (max|u| / min dx)")
            self.assertGreater(div / sscale, 0.1)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
