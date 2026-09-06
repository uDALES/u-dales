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
    CORRECTION_WARN_FRACTION,
    FLUX_UNITS,
    PARENT_DT_RTOL,
    SPEC,
    FaceMasks,
    NestGrid,
    NestingAlignmentError,
    NestingData,
    NestingRefinementError,
    NestingSchemaError,
    apply_divergence_correction,
    boundary_faces,
    check_alignment,
    check_time_axis,
    correction_report,
    discrete_divergence,
    face_masks_from_ibm,
    fluid_face_area,
    fluid_lateral_area,
    interpolate_child_fields,
    nesting_data_from_parent,
    net_volume_flux,
    project_initial_condition,
    read_nesting_file,
    refinement_ratios,
    refinement_ratios_by_axis,
    slabs_from_parent,
    stored_coordinates,
    validate_nesting_file,
    verify_stored_residual,
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
# W2 -- offset child origin
# --------------------------------------------------------------------------- #


class TestW2OffsetOrigin(unittest.TestCase):
    """W2: a child cut at (128, 128) m of the parent writes a file the solver accepts."""

    def test_round_trip_of_an_offset_child(self):
        parent = NestGrid.uniform(16, 16, 4, 256.0, 256.0, 40.0)
        child = NestGrid.uniform(16, 8, 8, 64.0, 32.0, 40.0, x0=128.0, y0=128.0)
        fields = [solenoidal_parent_fields(parent, seed=s) for s in (1, 2, 3)]
        data = nesting_data_from_parent(parent, child, 3, [0.0, 10.0, 20.0], fields)
        self.assertEqual(data.child_origin_x, 128.0)
        self.assertEqual(data.child_origin_y, 128.0)
        self.assertEqual(float(data.grid.xh[0]), 128.0)     # built in parent coordinates
        with TemporaryDirectory() as tmp:
            path = write_nesting_file(Path(tmp) / "offset.nc", data)
            attrs = validate_nesting_file(path)
            back = read_nesting_file(path)
            with Dataset(path, "r") as ds:
                xh = np.asarray(ds.variables["xh"][:])
                yh = np.asarray(ds.variables["yh"][:])
            # the reader's rule: |file - run| <= 1e-10 max(|file|, |run|, xlen), run at 0
            for arr, n, length in ((xh, 16, 64.0), (yh, 8, 32.0)):
                run = np.linspace(0.0, length, n + 1)
                self.assertLessEqual(np.max(np.abs(arr - run)), 1e-10 * length)
            self.assertEqual(float(attrs["child_origin_x"]), 128.0)
            self.assertEqual(float(attrs["child_origin_y"]), 128.0)
            self.assertAlmostEqual(float(attrs["xlen"]), 64.0)
        self.assertEqual(float(back.grid.xh[0]), 0.0)
        self.assertEqual(float(back.grid.yf[0]), 2.0)
        np.testing.assert_array_equal(back.grid.zh, child.zh)   # z is not shifted
        self.assertEqual(back.child_origin_x, 128.0)
        for name in data.slabs:
            np.testing.assert_array_equal(back.slabs[name], data.slabs[name], name)
        # the raw back-end stores the same coordinates
        with TemporaryDirectory() as tmp:
            raw = read_nesting_file(write_nesting_file(Path(tmp) / "offset.dat", data))
        np.testing.assert_array_equal(raw.grid.xh, back.grid.xh)

    def test_a_grid_that_is_neither_relative_nor_at_the_origin_is_refused(self):
        data = random_nesting_data(seed=201, ntime=2)
        shifted = NestGrid.from_faces(data.grid.xh + 500.0, data.grid.yh, data.grid.zh)
        data.grid = shifted
        data.child_origin_x = 0.0
        with self.assertRaises(ConfigurationError) as ctx:
            stored_coordinates(data)
        self.assertIn("child_origin_x", str(ctx.exception))
        data.child_origin_x = 500.0
        self.assertEqual(float(stored_coordinates(data)["xh"][0]), 0.0)


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


# --------------------------------------------------------------------------- #
# W6 -- the correction magnitude and its split across faces are reported
# --------------------------------------------------------------------------- #


class TestW6CorrectionReport(unittest.TestCase):
    """W6: |delta| against the boundary velocity scale, and who carries the lid flux."""

    def _lid_flux_case(self, north_outflow, ntime=2):
        # through-flow west -> east at 1 m/s (balanced), plus an unbalanced
        # outflow through the north face: that is the net flux the parent
        # pushed through the child's lid, which the correction must absorb
        data = uniform_nesting_data(seed=61, ntime=ntime)   # 8 x 6 x 5; dx 10, dy 5, dz 5
        for arr in data.slabs.values():
            arr[...] = 0.0
        faces = boundary_faces(data)
        faces["west"][...] = 1.0
        faces["east"][...] = 1.0
        faces["north"][...] = north_outflow
        return data

    def test_the_report_is_analytic(self):
        data = self._lid_flux_case(0.2)
        area_we, area_sn = 30.0 * 25.0, 80.0 * 25.0            # 750 and 2000 m2
        total = 2 * area_we + 2 * area_sn                      # 5500 m2
        phi = 0.2 * area_sn                                    # +400 m3/s, outward
        delta = -phi / total
        scale = np.sqrt((2 * area_we * 1.0 + area_sn * 0.2 ** 2) / total)
        with self.assertLogs("udprep.nesting", level=logging.INFO) as logs:
            apply_divergence_correction(data)
        rep = data.correction
        self.assertAlmostEqual(rep["residual_max_abs"], phi, places=9)
        self.assertAlmostEqual(rep["delta_max_abs"], abs(delta), places=12)
        self.assertAlmostEqual(rep["delta_at_max"], delta, places=12)
        self.assertAlmostEqual(rep["velocity_scale_at_max"], scale, places=12)
        self.assertAlmostEqual(rep["delta_fraction_max"], abs(delta) / scale, places=12)
        for face, area in (("west", area_we), ("east", area_we),
                           ("south", area_sn), ("north", area_sn)):
            self.assertAlmostEqual(rep["faces"][face]["area_fraction"], area / total, places=12)
        self.assertEqual(rep["inflow_faces"], ["west"])
        self.assertAlmostEqual(rep["faces"]["west"]["mean_outward_normal_velocity"], -1.0)
        self.assertAlmostEqual(rep["faces"]["north"]["mean_outward_normal_velocity"], 0.2)
        # 13.6 % of the velocity scale: this one must be a WARNING that names the inflow face
        self.assertGreater(rep["delta_fraction_max"], CORRECTION_WARN_FRACTION)
        self.assertTrue(rep["exceeded"])
        warnings = [m for m in logs.output if m.startswith("WARNING")]
        self.assertEqual(len(warnings), 1)
        self.assertIn("inflow faces (west)", warnings[0])
        self.assertIn(f"{100 * abs(delta) / scale:.1f} %", warnings[0])
        print(f"\n[W6] delta = {delta:.4f} m/s = {100 * abs(delta) / scale:.1f} % of the "
              f"boundary velocity scale {scale:.3f} m/s; split "
              + ", ".join(f"{f} {100 * rep['faces'][f]['area_fraction']:.0f} %" for f in FACES))
        # the report is the same object correction_report returns, and survives copy()
        self.assertEqual(correction_report(data), rep)
        self.assertEqual(data.copy().correction, rep)

    def test_a_small_correction_is_reported_quietly(self):
        data = self._lid_flux_case(0.02)                      # 1.2 % of the scale
        with self.assertLogs("udprep.nesting", level=logging.INFO) as logs:
            apply_divergence_correction(data)
        self.assertFalse(data.correction["exceeded"])
        self.assertFalse(any(m.startswith("WARNING") for m in logs.output))
        self.assertTrue(any("spread over the fluid lateral faces" in m for m in logs.output))

    def test_the_report_before_correcting_predicts_the_correction(self):
        data = self._lid_flux_case(0.2)
        predicted = correction_report(data)
        self.assertIsNone(data.correction)
        apply_divergence_correction(data)
        self.assertEqual(predicted, data.correction)


# --------------------------------------------------------------------------- #
# W7 -- tidy: units, cached residual, vertical ratios, FaceMasks from the IBM
# --------------------------------------------------------------------------- #


def uniform_nesting_data(nzone=2, ntime=2, seed=5, **kwargs):
    """Random slabs on a uniform 8 x 6 x 5 grid with dx = 10, dy = 5, dz = 5 m."""
    rng = np.random.default_rng(seed)
    grid = NestGrid.uniform(8, 6, 5, 80.0, 30.0, 25.0)
    slabs = {}
    for face in FACES:
        for component in COMPONENTS:
            from udprep.nesting import slab_shape
            slabs[f"{component}_{face}"] = rng.normal(
                size=(ntime,) + slab_shape(grid, nzone, face, component))
    return NestingData(grid=grid, nzone=nzone, times=30.0 * np.arange(ntime), slabs=slabs,
                       parent_dx=20.0, parent_dt=30.0, **kwargs)


class TestW7Tidy(unittest.TestCase):

    def test_flux_units_are_volume_flux_everywhere(self):
        self.assertEqual(FLUX_UNITS, "m3 s-1")
        data = uniform_nesting_data()
        with TemporaryDirectory() as tmp:
            path = write_nesting_file(Path(tmp) / "u.nc", data)
            with Dataset(path, "r") as ds:
                self.assertEqual(ds.variables["net_volume_flux"].units, FLUX_UNITS)
                self.assertEqual(ds.variables["flux_residual"].units, FLUX_UNITS)
        # and the one error message that quotes a flux uses the same string
        grid = NestGrid.uniform(6, 6, 4, 12.0, 12.0, 8.0)
        u = np.zeros(grid.component_shape("u")); u[-1] = 1.0
        with self.assertRaises(ConfigurationError) as ctx:
            project_initial_condition(grid, {"u": u, "v": np.zeros(grid.component_shape("v")),
                                             "w": np.zeros(grid.component_shape("w"))})
        self.assertIn(FLUX_UNITS, str(ctx.exception))

    def test_a_stale_cached_residual_is_refused_not_trusted(self):
        data = uniform_nesting_data(seed=71)
        apply_divergence_correction(data)
        boundary_faces(data)["east"][1] += 0.5          # changed after the correction
        with self.assertRaises(ConfigurationError) as ctx:
            with TemporaryDirectory() as tmp:
                write_nesting_file(Path(tmp) / "stale.nc", data)
        self.assertIn("flux_residual[1]", str(ctx.exception))

    def test_masks_used_by_the_correction_travel_with_the_data(self):
        data = uniform_nesting_data(seed=72)
        west = np.ones((6, 5), dtype=bool); west[:2, :2] = False
        masks = FaceMasks(west=west)
        apply_divergence_correction(data, masks)
        self.assertIs(data.masks, masks)
        verify_stored_residual(data)                    # recomputes with data.masks: fine
        data.masks = None                               # lost: the area no longer matches
        with self.assertRaises(ConfigurationError) as ctx:
            verify_stored_residual(data)
        self.assertIn("fluid_lateral_area", str(ctx.exception))
        verify_stored_residual(data, masks)             # passing them explicitly works

    def test_refinement_ratio_includes_y_and_z(self):
        data = uniform_nesting_data(seed=73)            # child dx 10, dy 5, dz 5
        data.parent_dx = 10.0                           # x ratio 1, y ratio 2 (dy defaults to dx)
        self.assertEqual(refinement_ratios(data)[0], 2.0)
        data.parent_dy = 5.0
        data.parent_dz = 25.0                           # z ratio 5: over the limit
        by_axis = refinement_ratios_by_axis(data)
        self.assertEqual((by_axis["x"], by_axis["y"], by_axis["z"]), (1.0, 1.0, 5.0))
        self.assertEqual(refinement_ratios(data)[0], 5.0)
        with self.assertRaises(NestingRefinementError):
            with TemporaryDirectory() as tmp:
                write_nesting_file(Path(tmp) / "z.nc", data)
        data.parent_dz = 10.0
        with TemporaryDirectory() as tmp:
            back = read_nesting_file(write_nesting_file(Path(tmp) / "z.nc", data))
        self.assertEqual((back.parent_dy, back.parent_dz), (5.0, 10.0))   # optional attributes

    def test_face_masks_from_the_ibm_mask_pin_the_fluid_area_analytically(self):
        data = uniform_nesting_data(seed=74, ntime=2)    # 8 x 6 x 5, dx 10, dy 5, dz 5
        fluid = np.ones((8, 6, 5), dtype=bool)
        fluid[0:2, 1:3, 0:3] = False                     # a cube against the west face
        fluid[4:7, 4:6, 0:2] = False                     # a slab against the north face
        masks = face_masks_from_ibm(fluid)
        np.testing.assert_array_equal(masks.west, fluid[0])
        np.testing.assert_array_equal(masks.north, fluid[:, -1])
        self.assertTrue(np.all(masks.east) and np.all(masks.south))
        # west: 6 x 5 cells of 5 x 5 m minus 2 x 3 solid; north: 8 x 5 of 10 x 5 minus 3 x 2
        west = 30 * 25.0 - 6 * 25.0
        north = 40 * 50.0 - 6 * 50.0
        expected = west + 30 * 25.0 + 40 * 50.0 + north
        self.assertEqual(expected, 5050.0)
        self.assertAlmostEqual(fluid_lateral_area(data, masks), expected, places=9)
        self.assertAlmostEqual(fluid_face_area(data, masks), expected, places=9)  # rho = 1
        self.assertAlmostEqual(fluid_lateral_area(data), 5500.0, places=9)
        original = {f: v.copy() for f, v in boundary_faces(data).items()}
        apply_divergence_correction(data, masks)
        self.assertLess(np.max(np.abs(net_volume_flux(data, masks))), 1e-10)
        for face in ("west", "north"):
            solid = ~np.asarray(getattr(masks, face))
            np.testing.assert_array_equal(boundary_faces(data)[face][:, solid],
                                          original[face][:, solid])
        # the solid convention of solid_c.txt
        same = face_masks_from_ibm(~fluid, solid=True)
        np.testing.assert_array_equal(same.west, masks.west)

    def test_the_module_is_a_udprep_section(self):
        from udprep.udprep import UDPrep
        self.assertIn("nesting", [spec.name for spec in UDPrep.SECTION_SPECS])
        self.assertIn("lnesting", SPEC.fields)
        self.assertIs(SPEC.defaults["lnesting"], False)
        self.assertEqual(SPEC.defaults["nest_timeinterp"], 2)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
