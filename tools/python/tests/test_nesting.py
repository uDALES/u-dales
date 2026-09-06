"""Unit tests P1--P17 of the nesting preprocessing (design doc §10.2, §10.6).

Each test isolates one mechanism of ``udprep.nesting`` and fails only if that
mechanism is wrong:

============ ============================================================
 P1           conservative interpolation, constants
 P2           conservative interpolation, flux identity
 P3           divergence preservation
 P4           local divergence bound and convergence
 P5           divergence correction
 P6           correction is minimal and shape-preserving
 P7           idempotence
 P8           masked correction
 P9           schema round-trip
 P10          refinement guard
 P11          container equivalence
 P12          schema 2: the stored flux_residual is the residual as stored
 P13          schema 2: the initial-condition block round-trips, both back-ends
 P14          the initial-condition projection makes the field solenoidal
 P15          the projection leaves every boundary-normal velocity untouched
 P16          the projection is a no-op on an already solenoidal field
 P17          schema 1 still writes, validates and reads (backwards compatibility)
============ ============================================================
"""

from __future__ import annotations

import sys
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
    INIT_VARIABLES,
    REQUIRED_GLOBAL_ATTRIBUTES,
    SCHEMA_VERSION,
    STAGGER,
    SUPPORTED_SCHEMA_VERSIONS,
    FaceMasks,
    NestGrid,
    NestingData,
    NestingRefinementError,
    NestingSchemaError,
    analytic_field,
    analytic_initial_fields,
    analytic_slabs,
    apply_divergence_correction,
    boundary_faces,
    conservative_interpolate,
    fluid_face_area,
    fluid_lateral_area,
    init_dimensions,
    initial_fields_from_fields,
    initial_fields_from_parent,
    project_initial_condition,
    interpolate_child_fields,
    nesting_data_from_parent,
    nesting_filename,
    net_volume_flux,
    read_nesting_file,
    refinement_ratios,
    slab_coordinates,
    slab_dimensions,
    slab_indices,
    slab_shape,
    slabs_from_fields,
    slabs_from_parent,
    sync_initial_condition,
    validate_nesting_file,
    write_analytic_nesting_file,
    write_nesting_file,
)

# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #


def make_grids(nx=4, ny=3, nz=5, rx=2, ry=2, rz=2, xlen=48.0, ylen=36.0, zsize=40.0):
    """A parent grid and an aligned child refined by ``(rx, ry, rz)``."""
    parent = NestGrid.uniform(nx, ny, nz, xlen, ylen, zsize)
    child = NestGrid.uniform(nx * rx, ny * ry, nz * rz, xlen, ylen, zsize)
    return parent, child


def random_parent_fields(parent, seed=0):
    rng = np.random.default_rng(seed)
    return tuple(rng.normal(size=parent.component_shape(c)) for c in COMPONENTS)


def solenoidal_parent_fields(parent, seed=1):
    """A parent field that is discretely divergence-free by construction.

    Built as the discrete curl of a random staggered vector potential, so the
    telescoping of the discrete divergence is exact to round-off.  The
    horizontal potential components are held constant on the bottom and top
    levels, which makes ``w = 0`` there -- an impermeable ground and lid, as in
    Case A -- so that the net flux through the *lateral* boundary vanishes too.
    Requires a uniform grid.
    """
    rng = np.random.default_rng(seed)
    itot, jtot, ktot = parent.itot, parent.jtot, parent.ktot
    dx, dy, dz = parent.dx[0], parent.dy[0], parent.dzf[0]
    ax = rng.normal(size=(itot, jtot + 1, ktot + 1))
    ay = rng.normal(size=(itot + 1, jtot, ktot + 1))
    az = rng.normal(size=(itot + 1, jtot + 1, ktot))
    ax[:, :, 0] = ax[:, :, -1] = 0.0
    ay[:, :, 0] = ay[:, :, -1] = 0.0
    u = (az[:, 1:, :] - az[:, :-1, :]) / dy - (ay[:, :, 1:] - ay[:, :, :-1]) / dz
    v = (ax[:, :, 1:] - ax[:, :, :-1]) / dz - (az[1:, :, :] - az[:-1, :, :]) / dx
    w = (ay[1:, :, :] - ay[:-1, :, :]) / dx - (ax[:, 1:, :] - ax[:, :-1, :]) / dy
    return u, v, w


def discrete_divergence(grid, u, v, w):
    """Cell-centred discrete divergence of a staggered field, shape (itot, jtot, ktot)."""
    return (
        (u[1:, :, :] - u[:-1, :, :]) / grid.dx[:, None, None]
        + (v[:, 1:, :] - v[:, :-1, :]) / grid.dy[None, :, None]
        + (w[:, :, 1:] - w[:, :, :-1]) / grid.dzf[None, None, :]
    )


def face_flux_over_parent_cells(child, field, component, ratios):
    """Sum child face fluxes over each coplanar parent face.

    Returns an array indexed by the parent face/cell indices, i.e.
    ``(itot_p+1, jtot_p, ktot_p)`` for ``u``.
    """
    rx, ry, rz = ratios
    steps = {"u": (rx, 1, 1), "v": (1, ry, 1), "w": (1, 1, rz)}[component]
    # keep only the child faces coplanar with a parent face, in the normal direction
    sel = [slice(None)] * 3
    axis = COMPONENTS.index(component)
    sel[axis] = slice(None, None, (rx, ry, rz)[axis])
    block = field[tuple(sel)]
    # multiply by the child face area, then sum over the tangential parent blocks
    dx, dy, dz = child.dx, child.dy, child.dzf
    area = {"u": (None, dy, dz), "v": (dx, None, dz), "w": (dx, dy, None)}[component]
    for ax, a in enumerate(area):
        if a is not None:
            shape = [1, 1, 1]
            shape[ax] = -1
            block = block * a.reshape(shape)
    shape = list(block.shape)
    new_shape = []
    for ax, n in enumerate(shape):
        r = (rx, ry, rz)[ax]
        if steps[ax] == 1 and r > 1:
            new_shape += [n // r, r]
        else:
            new_shape += [n]
    block = block.reshape(new_shape)
    # sum over the inserted refinement axes, from the back so indices stay valid
    axes = []
    pos = 0
    for ax, n in enumerate(shape):
        r = (rx, ry, rz)[ax]
        if steps[ax] == 1 and r > 1:
            axes.append(pos + 1)
            pos += 2
        else:
            pos += 1
    return block.sum(axis=tuple(axes))


def random_nesting_data(nzone=3, ntime=3, seed=2, stretched=True, **kwargs):
    """A NestingData with random slab values on a non-cubic, stretched grid."""
    rng = np.random.default_rng(seed)
    itot, jtot, ktot = 8, 6, 5
    xh = np.linspace(0.0, 80.0, itot + 1)
    yh = np.linspace(0.0, 30.0, jtot + 1)
    if stretched:
        zh = np.cumsum(np.concatenate(([0.0], 1.0 + np.arange(ktot) ** 1.5)))
    else:
        zh = np.linspace(0.0, 25.0, ktot + 1)
    grid = NestGrid.from_faces(xh, yh, zh)
    slabs = {}
    for face in FACES:
        for component in COMPONENTS:
            shape = (ntime,) + slab_shape(grid, nzone, face, component)
            slabs[f"{component}_{face}"] = rng.normal(size=shape)
    kwargs.setdefault("rhobf", 1.2 * np.exp(-grid.zf / 8000.0))
    kwargs.setdefault("rhobh", 1.2 * np.exp(-grid.zh / 8000.0))
    return NestingData(
        grid=grid,
        nzone=nzone,
        times=60.0 * np.arange(ntime),
        slabs=slabs,
        parent_dx=20.0,
        parent_dt=60.0,
        **kwargs,
    )


# --------------------------------------------------------------------------- #
# P1 -- conservative interpolation, constants
# --------------------------------------------------------------------------- #


class TestP1UniformField(unittest.TestCase):
    """P1: a uniform parent field maps to the identical uniform child field."""

    def test_uniform_parent_gives_uniform_child(self):
        for ratios in ((2, 2, 2), (4, 3, 2), (1, 4, 4)):
            with self.subTest(ratios=ratios):
                parent, child = make_grids(rx=ratios[0], ry=ratios[1], rz=ratios[2])
                values = {"u": -3.5, "v": 0.75, "w": 11.0}
                fields = [np.full(parent.component_shape(c), values[c]) for c in COMPONENTS]
                cu, cv, cw = interpolate_child_fields(parent, *fields, child)
                for component, arr in zip(COMPONENTS, (cu, cv, cw)):
                    self.assertEqual(arr.shape, child.component_shape(component))
                    np.testing.assert_allclose(arr, values[component], rtol=1e-15, atol=0.0)

    def test_uniform_field_survives_a_stretched_child_vertical(self):
        parent = NestGrid.uniform(4, 4, 4, 40.0, 40.0, 40.0)
        zh = np.array([0.0, 2.0, 6.0, 10.0, 16.0, 20.0, 26.0, 30.0, 40.0])
        child = NestGrid.from_faces(np.linspace(0.0, 40.0, 9),
                                    np.linspace(0.0, 40.0, 9), zh)
        field = np.full(parent.component_shape("u"), 2.25)
        coords = [child.component_coords("u", ax) for ax in range(3)]
        out = conservative_interpolate(parent, field, "u", *coords)
        np.testing.assert_allclose(out, 2.25, rtol=1e-15, atol=0.0)


# --------------------------------------------------------------------------- #
# P2 -- conservative interpolation, flux identity (the defining property)
# --------------------------------------------------------------------------- #


class TestP2FluxIdentity(unittest.TestCase):
    """P2: child face fluxes summed over a parent face equal the parent flux."""

    def test_child_fluxes_sum_to_the_parent_face_flux(self):
        for ratios in ((2, 2, 2), (4, 3, 2), (3, 1, 4)):
            with self.subTest(ratios=ratios):
                parent, child = make_grids(rx=ratios[0], ry=ratios[1], rz=ratios[2])
                fields = random_parent_fields(parent, seed=7)
                children = interpolate_child_fields(parent, *fields, child)
                for component, pf, cf in zip(COMPONENTS, fields, children):
                    got = face_flux_over_parent_cells(child, cf, component, ratios)
                    area = {"u": parent.dy[None, :, None] * parent.dzf[None, None, :],
                            "v": parent.dx[:, None, None] * parent.dzf[None, None, :],
                            "w": parent.dx[:, None, None] * parent.dy[None, :, None]}[component]
                    expected = pf * area
                    self.assertEqual(got.shape, expected.shape)
                    np.testing.assert_allclose(got, expected, rtol=1e-13, atol=1e-13)

    def test_coplanar_child_faces_reproduce_the_parent_value_exactly(self):
        parent, child = make_grids(rx=3, ry=2, rz=2)
        pu = random_parent_fields(parent, seed=11)[0]
        coords = [child.component_coords("u", ax) for ax in range(3)]
        cu = conservative_interpolate(parent, pu, "u", *coords)
        expected = np.repeat(np.repeat(pu, 2, axis=1), 2, axis=2)
        np.testing.assert_array_equal(cu[::3], expected)


# --------------------------------------------------------------------------- #
# P3 -- divergence preservation
# --------------------------------------------------------------------------- #


class TestP3DivergencePreservation(unittest.TestCase):
    """P3: a solenoidal parent gives a child target with zero parent-cell divergence."""

    def test_parent_cell_integrated_divergence_vanishes(self):
        ratios = (2, 3, 2)
        parent, child = make_grids(rx=ratios[0], ry=ratios[1], rz=ratios[2])
        pu, pv, pw = solenoidal_parent_fields(parent)
        np.testing.assert_allclose(discrete_divergence(parent, pu, pv, pw), 0.0, atol=1e-12)
        cu, cv, cw = interpolate_child_fields(parent, pu, pv, pw, child)
        fx = face_flux_over_parent_cells(child, cu, "u", ratios)
        fy = face_flux_over_parent_cells(child, cv, "v", ratios)
        fz = face_flux_over_parent_cells(child, cw, "w", ratios)
        net = ((fx[1:, :, :] - fx[:-1, :, :])
               + (fy[:, 1:, :] - fy[:, :-1, :])
               + (fz[:, :, 1:] - fz[:, :, :-1]))
        scale = max(np.max(np.abs(fx)), np.max(np.abs(fy)), np.max(np.abs(fz)))
        self.assertLess(np.max(np.abs(net)), 1e-12 * scale)

    def test_a_non_solenoidal_parent_is_not_silently_made_solenoidal(self):
        # guards the test above against trivially passing
        ratios = (2, 2, 2)
        parent, child = make_grids(rx=2, ry=2, rz=2)
        pu, pv, pw = random_parent_fields(parent, seed=13)
        cu, cv, cw = interpolate_child_fields(parent, pu, pv, pw, child)
        fx = face_flux_over_parent_cells(child, cu, "u", ratios)
        fy = face_flux_over_parent_cells(child, cv, "v", ratios)
        fz = face_flux_over_parent_cells(child, cw, "w", ratios)
        net = ((fx[1:, :, :] - fx[:-1, :, :])
               + (fy[:, 1:, :] - fy[:, :-1, :])
               + (fz[:, :, 1:] - fz[:, :, :-1]))
        self.assertGreater(np.max(np.abs(net)), 1.0)


# --------------------------------------------------------------------------- #
# P4 -- local divergence bound
# --------------------------------------------------------------------------- #


class TestP4LocalDivergence(unittest.TestCase):
    """P4: the interpolation creates no divergence inside a parent cell, and the
    residual converges at second order under parent refinement."""

    def test_child_cell_divergence_equals_its_parent_cell_divergence(self):
        rx, ry, rz = 2, 3, 2
        parent, child = make_grids(rx=rx, ry=ry, rz=rz)
        pu, pv, pw = random_parent_fields(parent, seed=17)
        cu, cv, cw = interpolate_child_fields(parent, pu, pv, pw, child)
        div_child = discrete_divergence(child, cu, cv, cw)
        div_parent = discrete_divergence(parent, pu, pv, pw)
        expanded = np.repeat(np.repeat(np.repeat(div_parent, rx, 0), ry, 1), rz, 2)
        np.testing.assert_allclose(div_child, expanded, rtol=1e-11, atol=1e-11)

    def test_residual_divergence_converges_at_second_order(self):
        # A continuum-solenoidal field with different wavenumbers in x and y, so
        # that its *discrete* divergence on the staggered grid is a genuine
        # O(h^2) sampling residual rather than an exact cancellation.
        length = 100.0
        kx = 2.0 * np.pi / length
        ky = 6.0 * np.pi / length

        def analytic(grid):
            ones = np.ones(grid.ktot)[None, None, :]
            u = (np.sin(kx * grid.xh)[:, None, None]
                 * np.cos(ky * grid.yf)[None, :, None] * ones)
            v = (-(kx / ky) * np.cos(kx * grid.xf)[:, None, None]
                 * np.sin(ky * grid.yh)[None, :, None] * ones)
            return u, v, np.zeros(grid.component_shape("w"))

        errors = []
        for n in (8, 16, 32):
            parent = NestGrid.uniform(n, n, 4, length, length, 40.0)
            child = NestGrid.uniform(2 * n, 2 * n, 8, length, length, 40.0)
            pu, pv, pw = analytic(parent)
            cu, cv, cw = interpolate_child_fields(parent, pu, pv, pw, child)
            def rms(values):
                return float(np.sqrt(np.mean(values ** 2)))

            child_residual = rms(discrete_divergence(child, cu, cv, cw))
            parent_residual = rms(discrete_divergence(parent, pu, pv, pw))
            # the interpolation adds nothing to what the parent sampling already has
            self.assertAlmostEqual(child_residual, parent_residual, delta=1e-12)
            errors.append(child_residual)
        for coarse, fine in zip(errors[:-1], errors[1:]):
            self.assertGreater(coarse / fine, 3.6)  # second order is a factor 4


# --------------------------------------------------------------------------- #
# P5 -- divergence correction
# --------------------------------------------------------------------------- #


class TestP5DivergenceCorrection(unittest.TestCase):
    """P5: Phi = 0 after correction, for random data, on every time level."""

    def test_flux_vanishes_on_every_time_level(self):
        data = random_nesting_data(seed=21, ntime=4)
        before = net_volume_flux(data)
        self.assertTrue(np.all(np.abs(before) > 1.0))
        residual = apply_divergence_correction(data)
        np.testing.assert_allclose(residual, before, rtol=0.0, atol=0.0)
        after = net_volume_flux(data)
        scale = fluid_face_area(data) * np.max(
            [np.max(np.abs(v)) for v in boundary_faces(data).values()]
        )
        self.assertLess(np.max(np.abs(after)), 1e-12 * scale)
        self.assertEqual(data.ntime, residual.size)

    def test_flux_matches_an_analytically_known_value(self):
        # pins the sign convention (outward positive) and the face areas
        # independently of the correction that uses them
        data = random_nesting_data(seed=23, ntime=1,
                                   rhobf=np.ones(5), rhobh=np.ones(6))
        for arr in data.slabs.values():
            arr[...] = 0.0
        faces = boundary_faces(data)
        faces["west"][...] = 1.0
        faces["east"][...] = 3.0
        faces["south"][...] = -2.0
        faces["north"][...] = 0.5
        grid = data.grid
        expected = ((3.0 - 1.0) * grid.ylen * grid.zsize
                    + (0.5 - (-2.0)) * grid.xlen * grid.zsize)
        np.testing.assert_allclose(net_volume_flux(data), [expected], rtol=1e-12)

    def test_pre_correction_residual_is_stored_and_the_flag_is_set(self):
        data = random_nesting_data(seed=22)
        before = net_volume_flux(data)
        apply_divergence_correction(data)
        np.testing.assert_array_equal(data.net_volume_flux, before)
        self.assertTrue(data.divergence_corrected)


# --------------------------------------------------------------------------- #
# P6 -- the correction is minimal and shape-preserving
# --------------------------------------------------------------------------- #


class TestP6CorrectionIsMinimal(unittest.TestCase):
    """P6: a constant increment per unit fluid area; face shapes unchanged."""

    def test_increment_is_a_single_constant_per_time_level(self):
        data = random_nesting_data(seed=31, ntime=3)
        original = {f: v.copy() for f, v in boundary_faces(data).items()}
        expected = -net_volume_flux(data) / fluid_face_area(data)
        apply_divergence_correction(data)
        sign = {"west": -1.0, "east": 1.0, "south": -1.0, "north": 1.0}
        for face, values in boundary_faces(data).items():
            increment = values - original[face]
            for n in range(data.ntime):
                span = np.ptp(increment[n])
                self.assertLess(span, 1e-12 * max(1.0, abs(expected[n])))
                np.testing.assert_allclose(
                    increment[n].mean(), sign[face] * expected[n], rtol=1e-12
                )

    def test_differences_within_a_face_are_unchanged(self):
        data = random_nesting_data(seed=32)
        original = {f: v.copy() for f, v in boundary_faces(data).items()}
        apply_divergence_correction(data)
        for face, values in boundary_faces(data).items():
            before = original[face] - original[face][:, :1, :1]
            after = values - values[:, :1, :1]
            np.testing.assert_allclose(after, before, rtol=1e-12, atol=1e-12)

    def test_interior_zone_columns_are_untouched(self):
        data = random_nesting_data(seed=33, nzone=3)
        before = {k: v.copy() for k, v in data.slabs.items()}
        apply_divergence_correction(data)
        for name, arr in data.slabs.items():
            component, face = name.split("_")
            if component != {"west": "u", "east": "u", "south": "v", "north": "v"}[face]:
                np.testing.assert_array_equal(arr, before[name])
                continue
            boundary = data.nzone if face in ("east", "north") else 0
            keep = [m for m in range(arr.shape[3]) if m != boundary]
            np.testing.assert_array_equal(arr[:, :, :, keep], before[name][:, :, :, keep])


# --------------------------------------------------------------------------- #
# P7 -- idempotence
# --------------------------------------------------------------------------- #


class TestP7Idempotence(unittest.TestCase):
    """P7: correcting an already-corrected file changes nothing to round-off."""

    def test_second_correction_is_a_no_op(self):
        data = random_nesting_data(seed=41, ntime=3)
        apply_divergence_correction(data)
        once = {k: v.copy() for k, v in data.slabs.items()}
        residual = apply_divergence_correction(data)
        scale = np.max([np.max(np.abs(v)) for v in once.values()])
        self.assertLess(np.max(np.abs(residual)), 1e-12 * fluid_face_area(data) * scale)
        for name, arr in data.slabs.items():
            np.testing.assert_allclose(arr, once[name], rtol=0.0, atol=1e-13 * scale)


# --------------------------------------------------------------------------- #
# P8 -- masked correction
# --------------------------------------------------------------------------- #


class TestP8MaskedCorrection(unittest.TestCase):
    """P8: with a partially solid boundary, the correction uses fluid area only."""

    def make_masks(self, grid):
        west = np.ones((grid.jtot, grid.ktot), dtype=bool)
        west[:2, :2] = False                      # a solid corner
        east = np.ones((grid.jtot, grid.ktot), dtype=bool)
        east[1, :] = False                        # a solid strip
        south = np.ones((grid.itot, grid.ktot), dtype=bool)
        south[:, 0] = False                       # a solid floor row
        north = np.ones((grid.itot, grid.ktot), dtype=bool)
        return FaceMasks(west=west, east=east, south=south, north=north)

    def test_masked_flux_vanishes_and_solid_faces_are_untouched(self):
        data = random_nesting_data(seed=51, ntime=2)
        masks = self.make_masks(data.grid)
        original = {f: v.copy() for f, v in boundary_faces(data).items()}
        apply_divergence_correction(data, masks)
        after = net_volume_flux(data, masks)
        scale = fluid_face_area(data, masks) * np.max(
            [np.max(np.abs(v)) for v in original.values()]
        )
        self.assertLess(np.max(np.abs(after)), 1e-12 * scale)
        for face, values in boundary_faces(data).items():
            solid = ~masks.get(data.grid, face).astype(bool)
            np.testing.assert_array_equal(values[:, solid], original[face][:, solid])

    def test_fluid_area_is_used_not_the_geometric_area(self):
        data = random_nesting_data(seed=52, ntime=2)
        masks = self.make_masks(data.grid)
        fluid = fluid_face_area(data, masks)
        geometric = fluid_face_area(data)
        self.assertLess(fluid, geometric)
        residual = net_volume_flux(data, masks)
        original = {f: v.copy() for f, v in boundary_faces(data).items()}
        apply_divergence_correction(data, masks)
        # the increment actually applied is -Phi / (fluid area), not -Phi / (geometric area)
        applied = (boundary_faces(data)["north"] - original["north"])[:, 0, 0]
        np.testing.assert_allclose(applied, -residual / fluid, rtol=1e-12)
        # correcting with the geometric area instead would leave this residual
        leftover = residual * (1.0 - fluid / geometric)
        self.assertGreater(np.max(np.abs(leftover)), 1e-6 * np.max(np.abs(residual)))
        np.testing.assert_allclose(net_volume_flux(data, masks), 0.0, atol=1e-9)


# --------------------------------------------------------------------------- #
# P9 -- schema round-trip
# --------------------------------------------------------------------------- #


class TestP9SchemaRoundTrip(unittest.TestCase):
    """P9: write -> read reproduces every field and attribute; a missing
    required attribute is rejected."""

    def setUp(self):
        self._tmp = TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def test_round_trip_reproduces_every_field(self):
        data = random_nesting_data(seed=61, ntime=3, parent_model="harmonie",
                                   child_origin_x=1200.0, child_origin_y=-450.0)
        apply_divergence_correction(data)
        path = write_nesting_file(self.tmp / nesting_filename(1), data)
        self.assertEqual(path.name, "nesting.inp.001.nc")
        attrs = validate_nesting_file(path)
        back = read_nesting_file(path)
        for name in data.slabs:
            np.testing.assert_array_equal(back.slabs[name], data.slabs[name])
        for name in ("xf", "xh", "yf", "yh", "zf", "zh"):
            np.testing.assert_array_equal(getattr(back.grid, name), getattr(data.grid, name))
        np.testing.assert_array_equal(back.times, data.times)
        np.testing.assert_array_equal(back.rhobf, data.rhobf)
        np.testing.assert_array_equal(back.rhobh, data.rhobh)
        np.testing.assert_array_equal(back.net_volume_flux, data.net_volume_flux)
        self.assertEqual(back.nzone, data.nzone)
        self.assertTrue(back.divergence_corrected)
        self.assertEqual(back.parent_model, "harmonie")
        self.assertEqual(back.parent_dx, data.parent_dx)
        self.assertEqual(back.parent_dt, data.parent_dt)
        self.assertEqual(back.child_origin_x, 1200.0)
        self.assertEqual(back.child_origin_y, -450.0)
        self.assertEqual(int(attrs["udales_nesting_schema"]), SCHEMA_VERSION)
        self.assertEqual(attrs["Conventions"], "CF-1.8")
        for key in REQUIRED_GLOBAL_ATTRIBUTES:
            self.assertIn(key, attrs)

    def test_dimensions_shapes_and_stagger_tags_match_the_contract(self):
        data = random_nesting_data(seed=62, ntime=2, nzone=3)
        path = write_nesting_file(self.tmp / "nesting.inp.002.nc", data)
        grid = data.grid
        with Dataset(path, "r") as ds:
            self.assertTrue(ds.dimensions["time"].isunlimited())
            for name, size in (("xf", grid.itot), ("xh", grid.itot + 1),
                               ("yf", grid.jtot), ("yh", grid.jtot + 1),
                               ("zf", grid.ktot), ("zh", grid.ktot + 1),
                               ("nz", data.nzone), ("nzh", data.nzone + 1)):
                self.assertEqual(len(ds.dimensions[name]), size, name)
            expected = {
                "u_west": ("time", "yf", "zf", "nzh"),
                "v_west": ("time", "yh", "zf", "nz"),
                "w_west": ("time", "yf", "zh", "nz"),
                "u_east": ("time", "yf", "zf", "nzh"),
                "v_east": ("time", "yh", "zf", "nz"),
                "w_east": ("time", "yf", "zh", "nz"),
                "u_south": ("time", "xh", "zf", "nz"),
                "v_south": ("time", "xf", "zf", "nzh"),
                "w_south": ("time", "xf", "zh", "nz"),
                "u_north": ("time", "xh", "zf", "nz"),
                "v_north": ("time", "xf", "zf", "nzh"),
                "w_north": ("time", "xf", "zh", "nz"),
            }
            for name, dims in expected.items():
                self.assertEqual(tuple(ds.variables[name].dimensions), dims, name)
                self.assertEqual(ds.variables[name].stagger, STAGGER[name[0]], name)
                self.assertEqual(ds.variables[name].dtype, np.dtype("f8"), name)

    def test_east_and_north_slab_index_convention(self):
        # slab index m maps to the *lowest* global index first
        grid = NestGrid.uniform(8, 6, 4, 80.0, 30.0, 20.0)
        nzone = 3
        for face, ntot in (("east", grid.itot), ("north", grid.jtot)):
            for component in COMPONENTS:
                idx = slab_indices(grid, nzone, face, component)
                self.assertEqual(idx[0], ntot - nzone)
                self.assertTrue(np.all(np.diff(idx) == 1))
        self.assertEqual(slab_indices(grid, nzone, "east", "u")[-1], grid.itot)
        self.assertEqual(slab_indices(grid, nzone, "east", "v")[-1], grid.itot - 1)
        self.assertEqual(slab_indices(grid, nzone, "north", "v")[-1], grid.jtot)
        for component in COMPONENTS:
            np.testing.assert_array_equal(
                slab_indices(grid, nzone, "west", component)[:1], [0]
            )

    def test_stored_values_sit_at_the_documented_stagger(self):
        # writes the analytic field and checks a handful of stored values
        # against f evaluated at the coordinate the contract implies
        grid = NestGrid.uniform(8, 6, 4, 80.0, 30.0, 20.0)
        nzone, times = 2, [0.0, 30.0]
        data = write_analytic_nesting_file(
            self.tmp / "nesting.inp.003.nc", grid, times, nzone
        )
        back = read_nesting_file(self.tmp / "nesting.inp.003.nc")
        for face in FACES:
            for component in COMPONENTS:
                xs, ys, zs = slab_coordinates(grid, nzone, face, component)
                for n, t in enumerate(times):
                    block = analytic_field(component, *np.ix_(xs, ys, zs), t)
                    perm = (1, 2, 0) if face in ("west", "east") else (0, 2, 1)
                    np.testing.assert_array_equal(
                        back.slabs[f"{component}_{face}"][n], block.transpose(perm)
                    )
        self.assertFalse(data.divergence_corrected)

    def test_missing_required_attribute_is_rejected(self):
        data = random_nesting_data(seed=63)
        for key in ("xlen", "itot", "rotation_deg", "udales_nesting_schema",
                    "divergence_corrected", "tool_version"):
            with self.subTest(attribute=key):
                path = self.tmp / f"missing_{key}.nc"
                write_nesting_file(path, data)
                with Dataset(path, "a") as ds:
                    ds.delncattr(key)
                with self.assertRaises(NestingSchemaError) as ctx:
                    validate_nesting_file(path)
                self.assertIn(key, str(ctx.exception))

    def test_corrupted_schema_is_rejected(self):
        data = random_nesting_data(seed=64)

        path = self.tmp / "bad_stagger.nc"
        write_nesting_file(path, data)
        with Dataset(path, "a") as ds:
            ds.variables["v_west"].stagger = "xh yf zf"
        with self.assertRaises(NestingSchemaError) as ctx:
            validate_nesting_file(path)
        self.assertIn("v_west", str(ctx.exception))

        path = self.tmp / "bad_itot.nc"
        write_nesting_file(path, data)
        with Dataset(path, "a") as ds:
            ds.itot = np.int32(data.grid.itot + 1)
        with self.assertRaises(NestingSchemaError):
            validate_nesting_file(path)

        path = self.tmp / "bad_rotation.nc"
        write_nesting_file(path, data)
        with Dataset(path, "a") as ds:
            ds.rotation_deg = np.float64(30.0)
        with self.assertRaises(NestingSchemaError):
            validate_nesting_file(path)

    def test_nan_values_are_an_error(self):
        data = random_nesting_data(seed=65)
        data.slabs["w_north"][0, 0, 0, 0] = np.nan
        with self.assertRaises(NestingSchemaError):
            write_nesting_file(self.tmp / "nan.nc", data)


# --------------------------------------------------------------------------- #
# P10 -- refinement guard
# --------------------------------------------------------------------------- #


class TestP10RefinementGuard(unittest.TestCase):
    """P10: spatial ratio > 4 or temporal > 30 is refused without an override."""

    def setUp(self):
        self._tmp = TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def test_spatial_ratio_beyond_four_is_refused(self):
        data = random_nesting_data(seed=71)     # child dx = 10 m, dy = 5 m
        data.parent_dx = 25.0                   # ratio 5 against dy
        spatial, _ = refinement_ratios(data)
        self.assertAlmostEqual(spatial, 5.0)
        with self.assertRaises(NestingRefinementError) as ctx:
            write_nesting_file(self.tmp / "too_coarse.nc", data)
        self.assertIn("spatial", str(ctx.exception))
        self.assertFalse((self.tmp / "too_coarse.nc").exists())

    def test_temporal_ratio_beyond_thirty_is_refused(self):
        data = random_nesting_data(seed=72)
        data.parent_dx = 10.0
        data.parent_dt = 60.0
        data.child_dt = 1.0                     # ratio 60
        _, temporal = refinement_ratios(data)
        self.assertAlmostEqual(temporal, 60.0)
        with self.assertRaises(NestingRefinementError) as ctx:
            write_nesting_file(self.tmp / "too_slow.nc", data)
        self.assertIn("temporal", str(ctx.exception))

    def test_override_allows_the_write(self):
        data = random_nesting_data(seed=73)
        data.parent_dx = 25.0
        data.parent_dt = 60.0
        data.child_dt = 1.0
        path = write_nesting_file(self.tmp / "override.nc", data, override=True)
        self.assertTrue(path.exists())
        validate_nesting_file(path)

    def test_ratios_at_the_limit_are_accepted(self):
        data = random_nesting_data(seed=74)           # times every 60 s, parent_dt 60
        data.parent_dx = 4.0 * min(data.grid.dx.min(), data.grid.dy.min())
        data.child_dt = 2.0                           # temporal ratio exactly 30
        path = write_nesting_file(self.tmp / "at_limit.nc", data)
        self.assertTrue(path.exists())

    def test_unknown_ratios_are_not_guessed(self):
        data = random_nesting_data(seed=75)
        data.parent_dx = 0.0                          # spatial unknown
        data.child_dt = None                          # temporal unknown
        self.assertEqual(refinement_ratios(data), (None, None))
        write_nesting_file(self.tmp / "unknown.nc", data)


# --------------------------------------------------------------------------- #
# P11 -- container equivalence
# --------------------------------------------------------------------------- #


class TestP11ContainerEquivalence(unittest.TestCase):
    """P11: the raw-binary and NetCDF back-ends give byte-identical buffers."""

    def setUp(self):
        self._tmp = TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def test_backends_agree_byte_for_byte(self):
        data = random_nesting_data(seed=81, ntime=3)
        apply_divergence_correction(data)
        write_nesting_file(self.tmp / "nesting.inp.001.nc", data)
        write_nesting_file(self.tmp / "nesting.inp.001.dat", data)
        from_nc = read_nesting_file(self.tmp / "nesting.inp.001.nc")
        from_raw = read_nesting_file(self.tmp / "nesting.inp.001.dat")
        for name in data.slabs:
            self.assertEqual(
                np.ascontiguousarray(from_nc.slabs[name]).tobytes(),
                np.ascontiguousarray(from_raw.slabs[name]).tobytes(),
                name,
            )
        for name in ("xf", "xh", "yf", "yh", "zf", "zh"):
            self.assertEqual(getattr(from_nc.grid, name).tobytes(),
                             getattr(from_raw.grid, name).tobytes(), name)
        self.assertEqual(from_nc.times.tobytes(), from_raw.times.tobytes())
        self.assertEqual(from_nc.rhobf.tobytes(), from_raw.rhobf.tobytes())
        self.assertEqual(from_nc.net_volume_flux.tobytes(),
                         from_raw.net_volume_flux.tobytes())
        self.assertEqual(from_nc.nzone, from_raw.nzone)
        self.assertEqual(from_nc.divergence_corrected, from_raw.divergence_corrected)
        self.assertEqual(from_nc.parent_dx, from_raw.parent_dx)


# --------------------------------------------------------------------------- #
# Supporting checks on the pieces the tests above rely on
# --------------------------------------------------------------------------- #


class TestSlabGeometry(unittest.TestCase):
    """The slab dimension names and shapes, which every other test assumes."""

    def test_slab_dimensions_follow_the_contract(self):
        self.assertEqual(slab_dimensions("west", "u"), ("yf", "zf", "nzh"))
        self.assertEqual(slab_dimensions("west", "v"), ("yh", "zf", "nz"))
        self.assertEqual(slab_dimensions("west", "w"), ("yf", "zh", "nz"))
        self.assertEqual(slab_dimensions("south", "u"), ("xh", "zf", "nz"))
        self.assertEqual(slab_dimensions("south", "v"), ("xf", "zf", "nzh"))
        self.assertEqual(slab_dimensions("south", "w"), ("xf", "zh", "nz"))
        self.assertEqual(slab_dimensions("east", "u"), slab_dimensions("west", "u"))
        self.assertEqual(slab_dimensions("north", "w"), slab_dimensions("south", "w"))

    def test_slabs_from_fields_matches_the_full_child_field(self):
        grid = NestGrid.uniform(8, 6, 5, 80.0, 30.0, 25.0)
        rng = np.random.default_rng(91)
        fields = {c: rng.normal(size=grid.component_shape(c)) for c in COMPONENTS}
        slabs = slabs_from_fields(grid, 2, fields["u"], fields["v"], fields["w"])
        # u_west[j, k, m] must be u[m, j, k]
        np.testing.assert_array_equal(
            slabs["u_west"], fields["u"][:3].transpose(1, 2, 0)
        )
        # v_north[i, k, m] must be v[i, jtot-nzone+m, k]
        np.testing.assert_array_equal(
            slabs["v_north"], fields["v"][:, grid.jtot - 2:, :].transpose(0, 2, 1)
        )

    def test_analytic_field_is_non_separable_and_component_specific(self):
        x, y, z = 12.0, 7.0, 3.0
        fu = analytic_field("u", x, y, z, 5.0)
        self.assertNotEqual(fu, analytic_field("v", x, y, z, 5.0))
        self.assertNotEqual(fu, analytic_field("u", y, x, z, 5.0))   # x/y transposed
        self.assertNotEqual(fu, analytic_field("u", x, y, z, 6.0))   # time dependence
        # non-separable: f(x1,y1) f(x2,y2) != f(x1,y2) f(x2,y1)
        a = analytic_field("u", 1.0, 2.0, z, 0.0) * analytic_field("u", 30.0, 40.0, z, 0.0)
        b = analytic_field("u", 1.0, 40.0, z, 0.0) * analytic_field("u", 30.0, 2.0, z, 0.0)
        self.assertNotAlmostEqual(a, b)

    def test_analytic_slabs_have_the_contract_shapes(self):
        grid = NestGrid.uniform(8, 6, 5, 80.0, 30.0, 25.0)
        slabs = analytic_slabs(grid, 2, [0.0, 10.0])
        for face in FACES:
            for component in COMPONENTS:
                self.assertEqual(
                    slabs[f"{component}_{face}"].shape,
                    (2,) + slab_shape(grid, 2, face, component),
                )

    def test_a_child_outside_the_parent_is_refused(self):
        parent = NestGrid.uniform(4, 4, 4, 40.0, 40.0, 40.0)
        child = NestGrid.uniform(8, 8, 8, 60.0, 40.0, 40.0)
        field = np.zeros(parent.component_shape("u"))
        coords = [child.component_coords("u", ax) for ax in range(3)]
        with self.assertRaises(ConfigurationError):
            conservative_interpolate(parent, field, "u", *coords)


class TestEndToEnd(unittest.TestCase):
    """The production path: interpolate a parent, correct, write, read back."""

    def setUp(self):
        self._tmp = TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def test_slab_only_interpolation_matches_the_full_field_path(self):
        parent, child = make_grids(nx=6, ny=6, nz=4, rx=2, ry=2, rz=2)
        fields = random_parent_fields(parent, seed=101)
        full = interpolate_child_fields(parent, *fields, child)
        direct = slabs_from_parent(parent, *fields, child=child, nzone=3)
        via_full = slabs_from_fields(child, 3, *full)
        for name in direct:
            np.testing.assert_array_equal(direct[name], via_full[name], name)

    def test_parent_to_file_gives_a_valid_flux_free_file(self):
        parent, child = make_grids(nx=6, ny=6, nz=4, rx=2, ry=2, rz=2)
        times = [0.0, 30.0, 60.0]
        fields = [solenoidal_parent_fields(parent, seed=200 + n) for n in range(len(times))]
        data = nesting_data_from_parent(parent, child, 3, times, fields,
                                        parent_model="udales")
        residual = apply_divergence_correction(data)
        # a solenoidal parent already satisfies Phi = 0 by construction (§1.3)
        scale = fluid_face_area(data) * np.max(
            [np.max(np.abs(v)) for v in data.slabs.values()]
        )
        self.assertLess(np.max(np.abs(residual)), 1e-12 * scale)
        path = write_nesting_file(self.tmp / nesting_filename(42), data)
        validate_nesting_file(path)
        back = read_nesting_file(path)
        self.assertLess(np.max(np.abs(net_volume_flux(back))), 1e-12 * scale)
        self.assertEqual(back.parent_model, "udales")
        self.assertAlmostEqual(back.parent_dt, 30.0)
        self.assertAlmostEqual(back.parent_dx, float(parent.dx[0]))


# --------------------------------------------------------------------------- #
# P12--P17 -- schema 2: the stored residual and the initial-condition block
# --------------------------------------------------------------------------- #


def closed_box_fields(grid, seed=7):
    """A random field with zero normal velocity on all six faces.

    Its net boundary flux is exactly zero, which is what the pure-Neumann
    projection needs, but its interior divergence is large -- so a projection
    has real work to do and cannot pass by doing nothing.
    """
    rng = np.random.default_rng(seed)
    u, v, w = (rng.normal(size=grid.component_shape(c)) for c in COMPONENTS)
    u[0, :, :] = u[-1, :, :] = 0.0
    v[:, 0, :] = v[:, -1, :] = 0.0
    w[:, :, 0] = w[:, :, -1] = 0.0
    return u, v, w


class TestP12StoredFluxResidual(unittest.TestCase):
    """P12: `flux_residual` describes the data **as stored**, not before correction.

    This is the whole point of the schema-2 addition: the solver validates
    against it instead of re-reading `4 x ntime` boundary slabs at
    initialisation (design section 10.6 item 3), so it has to be the residual of
    what is actually in the file.
    """

    def setUp(self):
        self._tmp = TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def test_it_is_the_post_correction_residual(self):
        data = random_nesting_data(seed=311)
        before = net_volume_flux(data)
        apply_divergence_correction(data)
        after = net_volume_flux(data)
        self.assertGreater(np.max(np.abs(before)), 1.0)      # not vacuous
        np.testing.assert_array_equal(data.net_volume_flux, before)
        np.testing.assert_allclose(data.flux_residual, after, rtol=0, atol=0)
        scale = fluid_face_area(data)
        self.assertLess(np.max(np.abs(data.flux_residual)), 1e-12 * scale)

    def test_an_uncorrected_file_stores_its_real_residual(self):
        data = random_nesting_data(seed=312)
        path = write_nesting_file(self.tmp / "raw.nc", data)
        back = read_nesting_file(path)
        np.testing.assert_allclose(back.flux_residual, net_volume_flux(data),
                                   rtol=1e-14, atol=0)
        self.assertGreater(np.max(np.abs(back.flux_residual)), 1.0)

    def test_the_fluid_lateral_area_is_written_and_matches(self):
        masks = FaceMasks(
            west=np.ones((6, 5), dtype=bool),
            east=np.ones((6, 5), dtype=bool),
            south=np.ones((8, 5), dtype=bool),
            north=np.ones((8, 5), dtype=bool),
        )
        masks.west[2, 1] = False
        masks.north[5, 3] = False
        data = random_nesting_data(seed=313)
        apply_divergence_correction(data, masks)
        path = write_nesting_file(self.tmp / "masked.nc", data)
        attrs = validate_nesting_file(path)
        self.assertAlmostEqual(float(attrs["fluid_lateral_area"]),
                               fluid_lateral_area(data, masks), places=9)
        # the geometric area is not the rho-weighted one this grid uses
        self.assertNotAlmostEqual(fluid_lateral_area(data, masks),
                                  fluid_face_area(data, masks), places=3)


class TestP13InitialConditionRoundTrip(unittest.TestCase):
    """P13: the full-domain block survives write -> read on both back-ends."""

    def setUp(self):
        self._tmp = TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def _data(self, seed=401):
        # rhobf == 1: the initial condition is projected with the solver's
        # density-free divergence operator (design F1), and the writer refuses
        # to carry a block on a file that claims any other density.
        data = random_nesting_data(seed=seed, rhobf=None, rhobh=None)
        data.initial_fields = initial_fields_from_fields(
            data.grid, *closed_box_fields(data.grid, seed=seed + 1)
        )
        return data

    def test_netcdf_round_trip_is_bitwise(self):
        data = self._data()
        path = write_nesting_file(self.tmp / "ic.nc", data)
        attrs = validate_nesting_file(path)
        self.assertEqual(int(attrs["has_initial_condition"]), 1)
        back = read_nesting_file(path)
        self.assertIsNotNone(back.initial_fields)
        for component in COMPONENTS:
            self.assertEqual(back.initial_fields[component].tobytes(),
                             data.initial_fields[component].tobytes(), component)

    def test_raw_backend_stores_the_same_bits(self):
        data = self._data(seed=402)
        write_nesting_file(self.tmp / "ic.nc", data)
        write_nesting_file(self.tmp / "ic.dat", data)
        a = read_nesting_file(self.tmp / "ic.nc")
        b = read_nesting_file(self.tmp / "ic.dat")
        for component in COMPONENTS:
            self.assertEqual(a.initial_fields[component].tobytes(),
                             b.initial_fields[component].tobytes(), component)
        self.assertEqual(a.flux_residual.tobytes(), b.flux_residual.tobytes())

    def test_the_dimensions_and_stagger_are_the_contract(self):
        data = self._data(seed=403)
        path = write_nesting_file(self.tmp / "ic.nc", data)
        with Dataset(path, "r") as ds:
            for component, name in zip(COMPONENTS, INIT_VARIABLES):
                var = ds.variables[name]
                self.assertEqual(tuple(var.dimensions), init_dimensions(component), name)
                self.assertEqual(var.stagger, STAGGER[component], name)

    def test_a_wrong_stagger_tag_is_rejected(self):
        data = self._data(seed=404)
        path = write_nesting_file(self.tmp / "ic.nc", data)
        with Dataset(path, "a") as ds:
            ds.variables["u_init"].stagger = "xf yf zf"
        with self.assertRaises(NestingSchemaError) as ctx:
            validate_nesting_file(path)
        self.assertIn("u_init", str(ctx.exception))

    def test_a_block_that_is_not_declared_is_rejected(self):
        data = self._data(seed=405)
        path = write_nesting_file(self.tmp / "ic.nc", data)
        with Dataset(path, "a") as ds:
            ds.setncattr("has_initial_condition", np.int32(0))
        with self.assertRaises(NestingSchemaError) as ctx:
            validate_nesting_file(path)
        self.assertIn("has_initial_condition", str(ctx.exception))

    def test_a_declared_block_that_is_missing_is_rejected(self):
        data = random_nesting_data(seed=406)
        path = write_nesting_file(self.tmp / "none.nc", data)
        with Dataset(path, "a") as ds:
            ds.setncattr("has_initial_condition", np.int32(1))
        with self.assertRaises(NestingSchemaError) as ctx:
            validate_nesting_file(path)
        self.assertIn("u_init", str(ctx.exception))

    def test_schema_1_refuses_to_carry_a_block(self):
        data = self._data(seed=407)
        with self.assertRaises(ConfigurationError):
            write_nesting_file(self.tmp / "v1.nc", data, schema=1)


class TestP14P16Projection(unittest.TestCase):
    """P14--P16: the projection of the initial condition (design 10.6 item 4)."""

    def _grid(self, stretched=False):
        if not stretched:
            return NestGrid.uniform(12, 10, 6, 24.0, 20.0, 12.0)
        zh = np.cumsum(np.concatenate(([0.0], np.linspace(1.0, 3.0, 6))))
        return NestGrid.from_faces(np.linspace(0.0, 24.0, 13),
                                   np.linspace(0.0, 20.0, 11), zh)

    def test_p14_the_projected_field_is_discretely_solenoidal(self):
        for stretched in (False, True):
            with self.subTest(stretched=stretched):
                grid = self._grid(stretched)
                fields = dict(zip(COMPONENTS, closed_box_fields(grid, seed=501)))
                before = np.max(np.abs(discrete_divergence(grid, *fields.values())))
                _, _, _, b, a = project_initial_condition(grid, fields)
                self.assertAlmostEqual(b, float(before), places=12)
                self.assertGreater(b, 1.0)          # not vacuous
                self.assertLess(a, 1e-12 * max(b, 1.0))

    def test_p15_boundary_normal_velocities_are_untouched(self):
        grid = self._grid()
        u, v, w = closed_box_fields(grid, seed=502)
        # give the lateral faces a non-zero but balanced normal velocity
        u[0, :, :] = 1.0
        u[-1, :, :] = 1.0
        fields = {"u": u, "v": v, "w": w}
        keep = {k: a.copy() for k, a in fields.items()}
        project_initial_condition(grid, fields)
        np.testing.assert_array_equal(fields["u"][0], keep["u"][0])
        np.testing.assert_array_equal(fields["u"][-1], keep["u"][-1])
        np.testing.assert_array_equal(fields["v"][:, 0], keep["v"][:, 0])
        np.testing.assert_array_equal(fields["v"][:, -1], keep["v"][:, -1])
        np.testing.assert_array_equal(fields["w"][:, :, 0], keep["w"][:, :, 0])
        np.testing.assert_array_equal(fields["w"][:, :, -1], keep["w"][:, :, -1])
        # and the interior did move, so the test is not vacuous
        self.assertGreater(float(np.max(np.abs(fields["v"] - keep["v"]))), 1e-3)

    def test_p16_it_is_a_no_op_on_a_solenoidal_field(self):
        grid = self._grid()
        fields = dict(zip(COMPONENTS, closed_box_fields(grid, seed=503)))
        project_initial_condition(grid, fields)
        keep = {k: a.copy() for k, a in fields.items()}
        _, _, _, b, a = project_initial_condition(grid, fields)
        moved = max(float(np.max(np.abs(fields[k] - keep[k]))) for k in fields)
        self.assertLess(moved, 1e-12)
        self.assertLess(a, 1e-12)

    def test_an_incompatible_field_is_refused_not_absorbed(self):
        grid = self._grid()
        u, v, w = closed_box_fields(grid, seed=504)
        u[-1, :, :] += 1.0          # net outflow with nothing to balance it
        with self.assertRaises(ConfigurationError) as ctx:
            project_initial_condition(grid, {"u": u, "v": v, "w": w})
        self.assertIn("net boundary flux", str(ctx.exception))

    def test_sync_takes_the_boundary_from_the_corrected_slabs(self):
        """The block's lateral faces must be the *corrected* boundary data.

        Otherwise the first substep sees a step change between the stored
        initial condition and the value `bcpup` imposes.
        """
        data = random_nesting_data(seed=505, stretched=False, rhobf=None, rhobh=None)
        data.initial_fields = initial_fields_from_fields(
            data.grid, *closed_box_fields(data.grid, seed=506)
        )
        apply_divergence_correction(data)
        faces = boundary_faces(data)
        block = data.initial_fields
        np.testing.assert_array_equal(block["u"][0, :, :], faces["west"][0])
        np.testing.assert_array_equal(block["u"][-1, :, :], faces["east"][0])
        np.testing.assert_array_equal(block["v"][:, 0, :], faces["south"][0])
        np.testing.assert_array_equal(block["v"][:, -1, :], faces["north"][0])
        peak = float(np.max(np.abs(discrete_divergence(data.grid, block["u"],
                                                       block["v"], block["w"]))))
        self.assertLess(peak, 1e-12)

    def test_a_non_unit_density_is_refused(self):
        """The one case where the two flux conventions genuinely disagree."""
        data = random_nesting_data(seed=508)      # rhobf = 1.2 exp(-z/8000)
        with self.assertRaises(ConfigurationError) as ctx:
            NestingData(
                grid=data.grid, nzone=data.nzone, times=data.times, slabs=data.slabs,
                rhobf=data.rhobf, rhobh=data.rhobh, parent_dt=data.parent_dt,
                initial_fields=initial_fields_from_fields(
                    data.grid, *closed_box_fields(data.grid, seed=509)),
            )
        self.assertIn("rhobf", str(ctx.exception))

    def test_the_analytic_helper_and_the_parent_path_agree(self):
        parent, child = make_grids(nx=6, ny=6, nz=4, rx=2, ry=2, rz=2)
        fields = random_parent_fields(parent, seed=507)
        block = initial_fields_from_parent(parent, *fields, child)
        full = interpolate_child_fields(parent, *fields, child)
        for component, arr in zip(COMPONENTS, full):
            np.testing.assert_array_equal(block[component], arr, component)
        grid = NestGrid.uniform(6, 5, 4, 12.0, 10.0, 8.0)
        analytic = analytic_initial_fields(grid, 3.0)
        for component in COMPONENTS:
            self.assertEqual(analytic[component].shape, grid.component_shape(component))
            xs, ys, zs = (grid.component_coords(component, ax) for ax in range(3))
            np.testing.assert_array_equal(
                analytic[component], analytic_field(component, *np.ix_(xs, ys, zs), 3.0)
            )


class TestP17SchemaOneCompatibility(unittest.TestCase):
    """P17: a schema 1 file must still write, validate and read.

    Backwards compatibility of the file format is a hard requirement: files
    written before schema 2 have no `flux_residual`, no `fluid_lateral_area` and
    no initial condition, and must load and run exactly as before.
    """

    def setUp(self):
        self._tmp = TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def test_both_versions_are_supported(self):
        self.assertEqual(SUPPORTED_SCHEMA_VERSIONS, (1, 2))
        self.assertEqual(SCHEMA_VERSION, 2)

    def test_a_schema_1_file_has_none_of_the_schema_2_items(self):
        data = random_nesting_data(seed=601)
        apply_divergence_correction(data)
        path = write_nesting_file(self.tmp / "v1.nc", data, schema=1)
        attrs = validate_nesting_file(path)
        self.assertEqual(int(attrs["udales_nesting_schema"]), 1)
        self.assertNotIn("has_initial_condition", attrs)
        self.assertNotIn("fluid_lateral_area", attrs)
        with Dataset(path, "r") as ds:
            self.assertNotIn("flux_residual", ds.variables)
            for name in INIT_VARIABLES:
                self.assertNotIn(name, ds.variables)

    def test_a_schema_1_file_reads_back_with_the_same_slabs(self):
        data = random_nesting_data(seed=602)
        apply_divergence_correction(data)
        path = write_nesting_file(self.tmp / "v1.nc", data, schema=1)
        back = read_nesting_file(path)
        for name in data.slabs:
            self.assertEqual(back.slabs[name].tobytes(), data.slabs[name].tobytes(), name)
        np.testing.assert_array_equal(back.net_volume_flux, data.net_volume_flux)
        self.assertIsNone(back.flux_residual)
        self.assertIsNone(back.fluid_lateral_area)
        self.assertIsNone(back.initial_fields)

    def test_a_schema_1_raw_file_round_trips(self):
        data = random_nesting_data(seed=603)
        apply_divergence_correction(data)
        write_nesting_file(self.tmp / "v1.dat", data, schema=1)
        back = read_nesting_file(self.tmp / "v1.dat")
        self.assertIsNone(back.flux_residual)
        self.assertIsNone(back.initial_fields)
        for name in data.slabs:
            self.assertEqual(back.slabs[name].tobytes(), data.slabs[name].tobytes(), name)

    def test_an_unknown_schema_is_still_rejected(self):
        data = random_nesting_data(seed=604)
        with self.assertRaises(ConfigurationError):
            write_nesting_file(self.tmp / "v9.nc", data, schema=9)
        path = write_nesting_file(self.tmp / "ok.nc", data)
        with Dataset(path, "a") as ds:
            ds.setncattr("udales_nesting_schema", np.int32(7))
        with self.assertRaises(NestingSchemaError) as ctx:
            validate_nesting_file(path)
        self.assertIn("udales_nesting_schema", str(ctx.exception))


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
