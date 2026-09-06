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
    PARENT_DT_RTOL,
    NestGrid,
    NestingAlignmentError,
    NestingData,
    NestingSchemaError,
    check_alignment,
    check_time_axis,
    discrete_divergence,
    interpolate_child_fields,
    nesting_data_from_parent,
    slabs_from_parent,
    validate_nesting_file,
    write_nesting_file,
)

from test_nesting import make_grids, random_nesting_data, solenoidal_parent_fields  # noqa: E402


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


# --------------------------------------------------------------------------- #
# W3 -- time axis
# --------------------------------------------------------------------------- #


class TestW3TimeAxis(unittest.TestCase):
    """W3: the stored time axis is the child's clock and parent_dt must match it."""

    def _kwargs(self, data):
        return dict(grid=data.grid, nzone=data.nzone, slabs=data.slabs,
                    rhobf=data.rhobf, rhobh=data.rhobh, parent_dx=data.parent_dx)

    def test_an_absolute_parent_time_axis_is_rejected(self):
        # the review's failure mode: a child frozen on the first level, silently
        data = random_nesting_data(seed=301, ntime=4)
        with self.assertRaises(ConfigurationError) as ctx:
            NestingData(times=3600.0 + data.times, parent_dt=60.0, **self._kwargs(data))
        self.assertIn("must start at exactly 0", str(ctx.exception))

    def test_a_non_monotone_time_axis_is_rejected(self):
        data = random_nesting_data(seed=302, ntime=4)
        for label, times in (("duplicate", [0.0, 60.0, 60.0, 180.0]),
                             ("unordered", [0.0, 120.0, 60.0, 180.0])):
            with self.subTest(case=label):
                with self.assertRaises(ConfigurationError) as ctx:
                    NestingData(times=times, parent_dt=60.0, **self._kwargs(data))
                self.assertIn("strictly increasing", str(ctx.exception))

    def test_a_parent_dt_that_is_not_the_cadence_is_rejected(self):
        data = random_nesting_data(seed=303, ntime=4)
        for label, dt in (("zero", 0.0), ("negative", -60.0), ("wrong", 30.0),
                          ("just outside tolerance", 60.0 * (1.0 + 2.0 * PARENT_DT_RTOL))):
            with self.subTest(case=label):
                with self.assertRaises(ConfigurationError) as ctx:
                    NestingData(times=data.times, parent_dt=dt, **self._kwargs(data))
                self.assertIn("parent_dt", str(ctx.exception))
        # inside the tolerance is fine, and the median is what is compared
        NestingData(times=data.times, parent_dt=60.0 * (1.0 + 0.5 * PARENT_DT_RTOL),
                    **self._kwargs(data))
        self.assertAlmostEqual(check_time_axis([0.0, 60.0, 120.0, 181.0], 60.0), 60.0)

    def test_a_single_level_needs_no_cadence(self):
        data = random_nesting_data(seed=304, ntime=1)
        NestingData(times=[0.0], parent_dt=0.0, **self._kwargs(data))
        with self.assertRaises(ConfigurationError):
            NestingData(times=[5.0], parent_dt=0.0, **self._kwargs(data))

    def test_validate_rejects_a_file_whose_axis_was_edited(self):
        data = random_nesting_data(seed=305, ntime=3)
        with TemporaryDirectory() as tmp:
            path = write_nesting_file(Path(tmp) / "t.nc", data)
            validate_nesting_file(path)
            with Dataset(path, "a") as ds:
                ds.variables["time"][0] = 1.0
            with self.assertRaises(NestingSchemaError) as ctx:
                validate_nesting_file(path)
            self.assertIn("start at exactly 0", str(ctx.exception))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
