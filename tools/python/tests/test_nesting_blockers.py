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
    CADENCE_COURANT_MAX,
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
    NestingWriter,
    analytic_slabs,
    apply_divergence_correction,
    boundary_faces,
    cadence_courant,
    check_alignment,
    check_refinement,
    check_time_axis,
    conservative_interpolate,
    correction_report,
    discrete_divergence,
    face_masks_from_ibm,
    fluid_face_area,
    initial_fields_from_fields,
    fluid_lateral_area,
    interpolate_child_fields,
    nesting_data_from_parent,
    nesting_diagnostics,
    net_volume_flux,
    project_initial_condition,
    read_nesting_file,
    refinement_ratios,
    refinement_ratios_by_axis,
    refinement_verdict,
    slabs_from_parent,
    stagger_masks_from_ibm,
    stored_coordinates,
    validate_nesting_file,
    verify_stored_residual,
    write_nesting_file,
)

from test_nesting import (  # noqa: E402
    closed_box_fields,
    face_flux_over_parent_cells,
    make_grids,
    random_nesting_data,
    random_parent_fields,
    solenoidal_parent_fields,
)


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
        cu, cv, cw = interpolate_child_fields(parent, pu, pv, pw, child, prolongation="constant")
        div = float(np.max(np.abs(discrete_divergence(child, cu, cv, cw))))
        self.assertLess(div, 1e-15)          # the review measured 3e-17
        # the linear scheme keeps the parent-cell integral, not the cell-wise value
        cu, cv, cw = interpolate_child_fields(parent, pu, pv, pw, child, prolongation="linear")
        vol = child.dx[:, None, None] * child.dy[None, :, None] * child.dzf[None, None, :]
        per_parent = (discrete_divergence(child, cu, cv, cw) * vol).reshape(
            8, 2, 6, 2, 6, 2).sum(axis=(1, 3, 5))
        self.assertLess(float(np.max(np.abs(per_parent))), 1e-12)

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
                                                          allow_misaligned=True,
                                                          prolongation="constant")
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
                                                      allow_misaligned=True,
                                                      prolongation="constant")
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
# W4 -- the per-level append writer
# --------------------------------------------------------------------------- #


def slab_bytes(data):
    return sum(arr.nbytes for arr in data.slabs.values())


class TestW4StreamingWriter(unittest.TestCase):
    """W4: NestingWriter writes level by level what write_nesting_file writes whole."""

    def setUp(self):
        self._tmp = TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def _metadata(self, data):
        return dict(rhobf=data.rhobf, rhobh=data.rhobh, parent_model="stream test",
                    parent_dx=data.parent_dx, parent_dt=data.parent_dt,
                    child_origin_x=data.child_origin_x, child_origin_y=data.child_origin_y,
                    child_dt=data.child_dt)

    def test_streamed_levels_equal_the_whole_record_route(self):
        data = random_nesting_data(seed=401, ntime=5, child_origin_x=64.0)
        west = np.ones((6, 5), dtype=bool); west[1:3, :2] = False
        masks = FaceMasks(west=west)
        raw = data.copy()
        untouched = {k: v.copy() for k, v in raw.slabs.items()}
        apply_divergence_correction(data, masks)
        whole = write_nesting_file(self.tmp / "whole.nc", data)
        with NestingWriter(self.tmp / "stream.nc", raw.grid, raw.nzone, masks=masks,
                           **self._metadata(raw)) as writer:
            for n in range(raw.ntime):
                info = writer.append_level(raw.times[n],
                                           {k: v[n] for k, v in raw.slabs.items()})
                self.assertEqual(info["level"], n)
                self.assertAlmostEqual(info["net_volume_flux"], data.net_volume_flux[n],
                                       places=9)
        validate_nesting_file(self.tmp / "stream.nc")
        a = read_nesting_file(whole)
        b = read_nesting_file(self.tmp / "stream.nc")
        scale = max(np.max(np.abs(v)) for v in data.slabs.values())
        for name in a.slabs:
            np.testing.assert_allclose(b.slabs[name], a.slabs[name], rtol=0, atol=1e-14 * scale,
                                       err_msg=name)
        np.testing.assert_array_equal(b.times, a.times)
        np.testing.assert_allclose(b.net_volume_flux, a.net_volume_flux, rtol=1e-12)
        self.assertLess(np.max(np.abs(b.flux_residual)), 1e-9)
        self.assertAlmostEqual(b.fluid_lateral_area, a.fluid_lateral_area, places=9)
        self.assertTrue(b.divergence_corrected)
        self.assertEqual(b.child_origin_x, 64.0)
        for name in ("xh", "yh", "zh"):
            np.testing.assert_array_equal(getattr(b.grid, name), getattr(a.grid, name))
        # the caller's arrays were not touched by the per-level correction
        for name in raw.slabs:
            np.testing.assert_array_equal(raw.slabs[name], untouched[name], name)
        diag = writer.diagnostics
        self.assertEqual(diag["ntime"], 5)
        self.assertAlmostEqual(diag["correction"]["delta_fraction_max"],
                               data.correction["delta_fraction_max"], places=12)
        self.assertEqual(diag["correction"]["inflow_faces"], data.correction["inflow_faces"])
        self.assertEqual(diag["cadence"]["C_dump"], cadence_courant(data, log=False)["C_dump"])
        self.assertTrue(diag["refinement"]["within_limits"])

    def test_the_initial_condition_is_synced_and_projected_on_the_stream(self):
        data = random_nesting_data(seed=402, ntime=3, rhobf=None, rhobh=None)
        data.initial_fields = initial_fields_from_fields(
            data.grid, *closed_box_fields(data.grid, seed=403))
        raw = data.copy()
        apply_divergence_correction(data)
        a = read_nesting_file(write_nesting_file(self.tmp / "whole.nc", data))
        with NestingWriter(self.tmp / "stream.nc", raw.grid, raw.nzone,
                           **self._metadata(raw)) as writer:
            for n in range(raw.ntime):
                writer.append_level(raw.times[n], {k: v[n] for k, v in raw.slabs.items()},
                                    initial_fields=raw.initial_fields if n == 0 else None)
        attrs = validate_nesting_file(self.tmp / "stream.nc")
        self.assertEqual(int(attrs["has_initial_condition"]), 1)
        b = read_nesting_file(self.tmp / "stream.nc")
        for c in COMPONENTS:
            np.testing.assert_allclose(b.initial_fields[c], a.initial_fields[c],
                                       rtol=0, atol=1e-12, err_msg=c)
        peak = float(np.max(np.abs(discrete_divergence(
            b.grid, b.initial_fields["u"], b.initial_fields["v"], b.initial_fields["w"]))))
        self.assertLess(peak, 1e-12)
        with self.assertRaises(ConfigurationError):
            with NestingWriter(self.tmp / "late.nc", raw.grid, raw.nzone,
                               **self._metadata(raw)) as writer:
                writer.append_level(0.0, {k: v[0] for k, v in raw.slabs.items()})
                writer.append_level(60.0, {k: v[1] for k, v in raw.slabs.items()},
                                    initial_fields=raw.initial_fields)
        self.assertFalse((self.tmp / "late.nc").exists())

    def test_a_failure_removes_the_partial_file(self):
        data = random_nesting_data(seed=404, ntime=4)
        level = {k: v[0] for k, v in data.slabs.items()}
        bad = {k: v.copy() for k, v in level.items()}
        bad["w_north"][0, 0, 0] = np.nan
        cases = {
            "NaN at level 2": ([(0.0, level), (60.0, level), (120.0, bad)], NestingSchemaError),
            "time going back": ([(0.0, level), (60.0, level), (30.0, level)], ConfigurationError),
            "first time not 0": ([(60.0, level)], ConfigurationError),
            "parent_dt off at close": ([(0.0, level), (30.0, level), (60.0, level)],
                                       ConfigurationError),
            "nothing appended": ([], ConfigurationError),
        }
        for label, (levels, exc) in cases.items():
            with self.subTest(case=label):
                path = self.tmp / f"{label.replace(' ', '_')}.nc"
                with self.assertRaises(exc):
                    with NestingWriter(path, data.grid, data.nzone,
                                       **self._metadata(data)) as writer:
                        for t, slabs in levels:
                            writer.append_level(t, slabs)
                self.assertFalse(path.exists(), label)

    def test_memory_high_water_mark_is_independent_of_the_number_of_levels(self):
        grid = NestGrid.uniform(16, 12, 10, 160.0, 120.0, 50.0)
        nzone = 3
        one_level = {k: v[0] for k, v in analytic_slabs(grid, nzone, [0.0]).items()}
        level_bytes = sum(a.nbytes for a in one_level.values())

        def stream(nlevels):
            tracemalloc.start()
            tracemalloc.reset_peak()
            with NestingWriter(self.tmp / f"n{nlevels}.nc", grid, nzone,
                               parent_dt=10.0, parent_dx=20.0) as writer:
                for n in range(nlevels):
                    slabs = analytic_slabs(grid, nzone, [10.0 * n])
                    writer.append_level(10.0 * n, {k: v[0] for k, v in slabs.items()})
            peak = tracemalloc.get_traced_memory()[1]
            tracemalloc.stop()
            return peak

        def whole(nlevels):
            tracemalloc.start()
            tracemalloc.reset_peak()
            times = 10.0 * np.arange(nlevels)
            data = NestingData(grid=grid, nzone=nzone, times=times,
                               slabs=analytic_slabs(grid, nzone, times),
                               parent_dt=10.0, parent_dx=20.0)
            write_nesting_file(self.tmp / f"w{nlevels}.nc", data)
            peak = tracemalloc.get_traced_memory()[1]
            tracemalloc.stop()
            return peak

        stream(5)                                        # warm up allocators
        s20, s200 = stream(20), stream(200)
        w200 = whole(200)
        print(f"\n[W4] one level = {level_bytes / 1e3:.0f} kB; streaming peak "
              f"{s20 / 1e6:.2f} MB at 20 levels, {s200 / 1e6:.2f} MB at 200; "
              f"whole-record route {w200 / 1e6:.2f} MB at 200")
        self.assertLess(s200, 1.25 * s20 + 2 * level_bytes)
        self.assertLess(s200, 20 * level_bytes)
        self.assertGreater(w200, 200 * level_bytes)      # the route being replaced
        # validation is level-wise too
        tracemalloc.start()
        tracemalloc.reset_peak()
        validate_nesting_file(self.tmp / "n200.nc")
        v200 = tracemalloc.get_traced_memory()[1]
        tracemalloc.stop()
        print(f"[W4] validate peak {v200 / 1e6:.2f} MB for a {200 * level_bytes / 1e6:.1f} MB record")
        self.assertLess(v200, 0.25 * 200 * level_bytes)


# --------------------------------------------------------------------------- #
# W5 -- the cadence criterion, and a refinement guard that is not inert
# --------------------------------------------------------------------------- #


class TestW5CadenceAndRefinement(unittest.TestCase):

    def _case(self, speed, seed=51):
        data = uniform_nesting_data(seed=seed, ntime=3)     # parent_dt 30 s, parent_dx 20 m
        for arr in data.slabs.values():
            arr[...] = 0.1
        boundary_faces(data)["east"][1, 2, 3] = speed        # one fast face at level 1
        return data

    def test_c_dump_is_max_normal_speed_times_dt_over_dx(self):
        self.assertEqual(CADENCE_COURANT_MAX, 2.0)
        data = self._case(3.0)
        rep = cadence_courant(data, log=False)
        self.assertAlmostEqual(rep["C_dump"], 3.0 * 30.0 / 20.0)      # 4.5
        self.assertEqual((rep["max_normal_speed"], rep["level"]), (3.0, 1))
        self.assertTrue(rep["exceeded"])
        # a solid face does not count
        east = np.ones((6, 5), dtype=bool); east[2, 3] = False
        self.assertAlmostEqual(cadence_courant(data, FaceMasks(east=east), log=False)["C_dump"],
                               0.1 * 30.0 / 20.0)
        with TemporaryDirectory() as tmp:
            with self.assertLogs("udprep.nesting", level=logging.WARNING) as logs:
                write_nesting_file(Path(tmp) / "fast.nc", data)
            self.assertTrue(any("C_dump" in m and "4.50" in m for m in logs.output))
            quiet = self._case(1.0)                                     # C_dump = 1.5
            with self.assertLogs("udprep.nesting", level=logging.INFO) as logs:
                write_nesting_file(Path(tmp) / "slow.nc", quiet)
            self.assertFalse(any(m.startswith("WARNING") for m in logs.output))
        data.parent_dx = 0.0
        self.assertIsNone(cadence_courant(data, log=False)["C_dump"])
        print(f"\n[W5] C_dump = {rep['C_dump']:.2f} at max|u_n| = {rep['max_normal_speed']:g} m/s, "
              f"dt_P = {rep['parent_dt']:g} s, dx_P = {rep['parent_dx']:g} m")

    def test_the_refinement_guard_yields_a_verdict_and_needs_a_stated_reason(self):
        data = uniform_nesting_data(seed=52)
        data.child_dt = 0.5                                             # temporal ratio 60
        verdict = refinement_verdict(data)
        self.assertFalse(verdict["within_limits"])
        self.assertEqual(verdict["temporal"], 60.0)
        self.assertEqual(len(verdict["violations"]), 1)
        with self.assertRaises(NestingRefinementError) as ctx:
            check_refinement(data)
        self.assertIn("allow_refinement_violation=True", str(ctx.exception))
        with self.assertLogs("udprep.nesting", level=logging.WARNING) as logs:
            allowed = check_refinement(data, allow_refinement_violation=True,
                                       reason="cadence study C0")
        self.assertTrue(allowed["allowed"])
        self.assertEqual(allowed["reason"], "cadence study C0")
        self.assertIn("cadence study C0", logs.output[0])
        with self.assertLogs("udprep.nesting", level=logging.WARNING) as logs:
            self.assertTrue(check_refinement(data, override=True)["allowed"])  # the old name
        self.assertIn("no reason given", logs.output[0])
        with TemporaryDirectory() as tmp:
            with self.assertRaises(NestingRefinementError):
                write_nesting_file(Path(tmp) / "no.nc", data)
            with self.assertLogs("udprep.nesting", level=logging.WARNING) as logs:
                write_nesting_file(Path(tmp) / "yes.nc", data, allow_refinement_violation=True,
                                   refinement_reason="test")
            self.assertTrue(any("allowed by the caller: test" in m for m in logs.output))
        data.child_dt = 1.0
        self.assertTrue(check_refinement(data)["within_limits"])

    def test_diagnostics_bundle_the_three_verdicts(self):
        data = self._case(3.0)
        diag = nesting_diagnostics(data)
        self.assertEqual(set(diag), {"cadence", "correction", "refinement"})
        self.assertTrue(diag["cadence"]["exceeded"])
        self.assertIn("delta_fraction_max", diag["correction"])
        self.assertTrue(diag["refinement"]["within_limits"])
        self.assertIsNone(data.correction)                              # nothing was changed


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


# --------------------------------------------------------------------------- #
# W8 -- piecewise-linear prolongation
# --------------------------------------------------------------------------- #


class TestW8LinearProlongation(unittest.TestCase):
    """W8: linear within a parent cell, conservative to round-off, no staircase."""

    def test_a_field_linear_in_z_prolongs_without_a_staircase(self):
        for r in (2, 4):
            with self.subTest(ratio=r):
                parent = NestGrid.uniform(4, 4, 8, 40.0, 40.0, 32.0)
                child = NestGrid.uniform(4 * r, 4 * r, 8 * r, 40.0, 40.0, 32.0)
                a, b = 1.5, 0.07
                pu = np.broadcast_to((a + b * parent.zf)[None, None, :],
                                     parent.component_shape("u")).copy()
                pv = np.broadcast_to((2.0 - 0.5 * b * parent.zf)[None, None, :],
                                     parent.component_shape("v")).copy()
                pw = np.zeros(parent.component_shape("w"))
                slabs = slabs_from_parent(parent, pu, pv, pw, child=child, nzone=2)
                # u on the x-faces of the west slab, (yf, zf, nzh); v on the same
                # slab, (yh, zf, nz): both must be the linear profile at the
                # child's own zf, to round-off
                want_u = (a + b * child.zf)[None, :, None]
                want_v = (2.0 - 0.5 * b * child.zf)[None, :, None]
                err_u = float(np.max(np.abs(slabs["u_west"] - want_u)))
                err_v = float(np.max(np.abs(slabs["v_west"] - want_v)))
                self.assertLess(err_u, 1e-13)
                self.assertLess(err_v, 1e-13)
                # and the staircase the constant scheme leaves, for the record
                const = slabs_from_parent(parent, pu, pv, pw, child=child, nzone=2,
                                          prolongation="constant")
                stair = float(np.max(np.abs(const["u_west"] - want_u)))
                print(f"\n[W8] r = {r}: linear max error {err_u:.1e}, constant staircase "
                      f"{stair:.4f} m/s = du/dz dz_c (r-1)/2 = {b * child.dzf[0] * (r - 1) / 2:.4f}")
                self.assertAlmostEqual(stair, b * child.dzf[0] * (r - 1) / 2, places=12)

    def test_fluxes_and_phi_match_the_constant_scheme_to_round_off(self):
        for ratios in ((2, 3, 2), (4, 1, 2)):
            with self.subTest(ratios=ratios):
                parent, child = make_grids(nx=6, ny=6, nz=4, rx=ratios[0], ry=ratios[1],
                                           rz=ratios[2])
                fields = [solenoidal_parent_fields(parent, seed=s) for s in (81, 82)]
                fields[1] = tuple(f + 0.3 * g for f, g in
                                  zip(fields[1], random_parent_fields(parent, seed=99)))
                out = {}
                for scheme in ("constant", "linear"):
                    data = nesting_data_from_parent(parent, child, 3, [0.0, 10.0], fields,
                                                    prolongation=scheme)
                    phi_before = net_volume_flux(data).copy()
                    apply_divergence_correction(data)
                    out[scheme] = (phi_before, data.flux_residual.copy(), data)
                scale = fluid_face_area(out["linear"][2]) * max(
                    np.max(np.abs(a)) for f in fields for a in f)
                np.testing.assert_allclose(out["linear"][0], out["constant"][0],
                                           rtol=0, atol=1e-13 * scale)
                np.testing.assert_allclose(out["linear"][1], out["constant"][1],
                                           rtol=0, atol=1e-13 * scale)
                # every parent-face integral is preserved by the linear scheme
                for component, pf in zip(COMPONENTS, fields[1]):
                    coords = [child.component_coords(component, ax) for ax in range(3)]
                    cf = conservative_interpolate(parent, pf, component, *coords,
                                                  prolongation="linear")
                    got = face_flux_over_parent_cells(child, cf, component, ratios)
                    area = {"u": parent.dy[None, :, None] * parent.dzf[None, None, :],
                            "v": parent.dx[:, None, None] * parent.dzf[None, None, :],
                            "w": parent.dx[:, None, None] * parent.dy[None, :, None]}[component]
                    np.testing.assert_allclose(got, pf * area, rtol=1e-13, atol=1e-13)

    def test_a_solid_neighbour_does_not_pollute_the_slope(self):
        parent = NestGrid.uniform(6, 6, 4, 60.0, 60.0, 40.0)
        child = NestGrid.uniform(12, 12, 8, 60.0, 60.0, 40.0)
        b = 0.1
        pu = np.broadcast_to((b * parent.yf)[None, :, None], parent.component_shape("u")).copy()
        fluid = np.ones((6, 6, 4), dtype=bool)
        fluid[:, 3, :] = False                          # a solid wall across y
        masks = stagger_masks_from_ibm(fluid)
        pu[~masks["u"]] = 0.0                           # the IBM parent has u = 0 in the wall
        coords = [child.component_coords("u", ax) for ax in range(3)]
        polluted = conservative_interpolate(parent, pu, "u", *coords)
        guarded = conservative_interpolate(parent, pu, "u", *coords, parent_mask=masks["u"])
        want = b * child.yf
        # child cells 4, 5 lie in parent cell 2, the fluid neighbour of the wall
        self.assertLess(float(np.max(np.abs(guarded[:, 4:6, :] - want[None, 4:6, None]))), 1e-13)
        self.assertGreater(float(np.max(np.abs(polluted[:, 4:6, :] - want[None, 4:6, None]))), 0.1)
        np.testing.assert_array_equal(guarded[:, 6:8, :], 0.0)     # the wall itself stays 0
        # still conservative
        got = face_flux_over_parent_cells(child, guarded, "u", (2, 2, 2))
        area = parent.dy[None, :, None] * parent.dzf[None, None, :]
        np.testing.assert_allclose(got, pu * area, rtol=1e-13, atol=1e-13)
        # the stagger masks follow the solver's IIu/IIv/IIw rule
        self.assertFalse(masks["v"][:, 3, :].any() or masks["v"][:, 4, :].any())
        self.assertTrue(masks["v"][:, 2, :].all() and masks["v"][:, 5, :].all())
        self.assertEqual(masks["u"].shape, parent.component_shape("u"))

    def test_log_profile_interior_mean_error_constant_vs_linear(self):
        ustar, kappa, z0 = 0.3, 0.4, 0.1
        parent = NestGrid.uniform(4, 4, 16, 40.0, 40.0, 32.0)
        child = NestGrid.uniform(8, 8, 32, 40.0, 40.0, 32.0)
        profile = lambda z: (ustar / kappa) * np.log(z / z0)       # noqa: E731
        pu = np.broadcast_to(profile(parent.zf)[None, None, :], parent.component_shape("u")).copy()
        pv = np.zeros(parent.component_shape("v"))
        pw = np.zeros(parent.component_shape("w"))
        truth = profile(child.zf)
        above = child.zf > 4.0                                     # away from the curved base
        rms = {}
        for scheme in ("constant", "linear"):
            cu, _, _ = interpolate_child_fields(parent, pu, pv, pw, child, prolongation=scheme)
            err = cu.mean(axis=(0, 1)) - truth
            rms[scheme] = (float(np.sqrt(np.mean(err ** 2))) / ustar,
                           float(np.sqrt(np.mean(err[above] ** 2))) / ustar,
                           float(np.max(np.abs(np.diff(err[above])))) / ustar)
        print(f"\n[W8] log profile, r = 2, interior mean error in u*: constant rms "
              f"{rms['constant'][0]:.3f} (above 4 m {rms['constant'][1]:.3f}, sawtooth "
              f"{rms['constant'][2]:.3f}); linear rms {rms['linear'][0]:.3f} (above 4 m "
              f"{rms['linear'][1]:.3f}, sawtooth {rms['linear'][2]:.3f})")
        self.assertLess(rms["linear"][0], 0.6 * rms["constant"][0])
        self.assertLess(rms["linear"][1], 0.1 * rms["constant"][1])
        self.assertLess(rms["linear"][2], 0.1 * rms["constant"][2])

    def test_an_unknown_prolongation_is_refused(self):
        parent, child = make_grids()
        pu = np.zeros(parent.component_shape("u"))
        with self.assertRaises(ConfigurationError):
            conservative_interpolate(parent, pu, "u",
                                     *[child.component_coords("u", ax) for ax in range(3)],
                                     prolongation="cubic")


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
