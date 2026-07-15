from __future__ import annotations

from types import SimpleNamespace
import unittest

import numpy as np
import trimesh

try:
    from _common import REPO_ROOT, RUN_SLOW_TESTS, requires_slow_tests
except ModuleNotFoundError:
    from tools.python.tests._common import REPO_ROOT, RUN_SLOW_TESTS, requires_slow_tests

# Gate BEFORE importing the solver: importing udprep.directshortwave triggers
# numba compilation, so even a skipped test class would pay minutes of import
# cost during ordinary discovery. Raising SkipTest at module level makes
# unittest discovery report the whole module as skipped without importing it.
if not RUN_SLOW_TESTS:
    raise unittest.SkipTest(
        "slow numba direct-shortwave tests; set UDALES_RUN_SLOW_TESTS=1 to run"
    )

from udgeom import UDGeom
from udprep.directshortwave import DirectShortwaveSolver

def _build_flat_terrain_fixture():
    """Create a minimal flat terrain with a single 1x1 ground cell."""
    vertices = np.array(
        [
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ],
        dtype=float,
    )
    faces = np.array([[0, 1, 2], [0, 2, 3]], dtype=int)
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    geom = UDGeom(stl=mesh)

    sim = SimpleNamespace(
        dx=1.0,
        dy=1.0,
        dzt=np.array([1.0, 1.0], dtype=float),
        itot=1,
        jtot=1,
        ktot=2,
        xlen=1.0,
        ylen=1.0,
        zsize=2.0,
        zm=np.array([0.0, 1.0], dtype=float),
        Sc=np.array([[[True, False]]], dtype=bool),
        facsec={
            "c": {
                "facid": np.array([0, 1], dtype=int),
                "area": np.array([0.5, 0.5], dtype=float),
                "locs": np.array([[0, 0, 1], [0, 0, 1]], dtype=int),
            }
        },
        facs={"area": np.array([0.5, 0.5], dtype=float)},
        geom=geom,
    )
    return sim, mesh

@requires_slow_tests
class TestDirectShortwaveFlatTerrain(unittest.TestCase):
    def setUp(self):
        self.sim, self.mesh = _build_flat_terrain_fixture()
        azimuth_deg = 20.0
        elevation_deg = 15.0
        az = np.deg2rad(azimuth_deg)
        el = np.deg2rad(elevation_deg)
        self.nsun = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)], dtype=float)
        self.irradiance = 800.0
        self.expected_flux = self.irradiance * np.sin(el)

    def test_facsec_matches_analytic_flux_on_flat_terrain(self):
        solver = DirectShortwaveSolver(
            self.sim,
            method="facsec",
            surface_mesh=self.mesh,
            ray_density=12.0,
            ray_jitter=0.0,
        )
        sdir, veg, bud = solver.compute(nsun=self.nsun, irradiance=self.irradiance, periodic_xy=False)

        self.assertEqual(veg.size, 0)
        self.assertTrue(np.allclose(sdir, self.expected_flux, rtol=0.0, atol=1.0e-9))
        self.assertAlmostEqual(bud["fac"], self.expected_flux, delta=1.0e-9)

    def test_moller_matches_analytic_total_flux_on_flat_terrain(self):
        solver = DirectShortwaveSolver(
            self.sim,
            method="moller",
            surface_mesh=self.mesh,
            ray_density=64.0,
            ray_jitter=0.0,
        )
        sdir, veg, bud = solver.compute(nsun=self.nsun, irradiance=self.irradiance, periodic_xy=False)

        self.assertEqual(veg.size, 0)
        self.assertTrue(np.all(np.isfinite(sdir)))
        self.assertAlmostEqual(bud["fac"], self.expected_flux, delta=1.0e-9)
        self.assertAlmostEqual(float(np.sum(sdir * self.mesh.area_faces)), self.expected_flux, delta=1.0e-9)

    def test_scanline_matches_analytic_total_flux_on_flat_terrain(self):
        solver = DirectShortwaveSolver(
            self.sim,
            method="scanline",
            surface_mesh=self.mesh,
            ray_density=64.0,
            ray_jitter=0.0,
        )
        sdir, veg, bud = solver.compute(nsun=self.nsun, irradiance=self.irradiance, periodic_xy=False)

        self.assertEqual(veg.size, 0)
        self.assertTrue(np.all(np.isfinite(sdir)))
        self.assertAlmostEqual(bud["fac"], self.expected_flux, delta=0.5)
        self.assertAlmostEqual(float(np.sum(sdir * self.mesh.area_faces)), self.expected_flux, delta=0.5)

def _build_veg_column_fixture():
    """Single 1x1 column, no solid; one vegetation cell at k=0.

    itot=jtot=1, two vertical cells (k=0 veg, k=1 empty above it). A vertical
    ray enters the empty top cell (r_in stays 1.0), then attenuates over the
    full dz of the vegetation cell, then exits the domain (no facet).
    """
    # Dummy ground quad only supplies face normals/areas; facsec ignores mesh
    # geometry for hits (it uses the solid mask), so it is never struck here.
    vertices = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        dtype=float,
    )
    faces = np.array([[0, 1, 2], [0, 2, 3]], dtype=int)
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    geom = UDGeom(stl=mesh)

    sim = SimpleNamespace(
        dx=1.0,
        dy=1.0,
        dzt=np.array([1.0, 1.0], dtype=float),
        itot=1,
        jtot=1,
        ktot=2,
        xlen=1.0,
        ylen=1.0,
        zsize=2.0,
        zm=np.array([0.0, 1.0], dtype=float),
        Sc=np.zeros((1, 1, 2), dtype=bool),
        facsec={
            "c": {
                "facid": np.empty((0,), dtype=int),
                "area": np.empty((0,), dtype=float),
                "locs": np.empty((0, 3), dtype=int),
            }
        },
        facs={"area": np.array([0.5, 0.5], dtype=float)},
        geom=geom,
    )
    return sim, mesh

def _build_veg_over_ground_block_fixture():
    """Single 1x1 column: solid ground block at k=0, vegetation cell at k=1
    directly above it, and a horizontal ground facet at z=1 (the top of the
    block, which shares grid cell k=1 with the vegetation).

    A vertical ray enters the domain at the top of the vegetation cell
    (r_in = 1), attenuates over the full dz of the vegetation, then strikes the
    ground facet at z=1 which lives in the SAME (i,j,k) cell as the vegetation.
    This is the moller double-counting scenario (facet + veg in one cell), and
    because the block is solid it is representable by facsec too, so the two
    kernels can be compared on an identical scene.
    """
    # Ground quad at z=1 (top face of the block), normal +z (up).
    vertices = np.array(
        [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]],
        dtype=float,
    )
    faces = np.array([[0, 1, 2], [0, 2, 3]], dtype=int)
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    geom = UDGeom(stl=mesh)

    solid = np.zeros((1, 1, 2), dtype=bool)
    solid[0, 0, 0] = True  # ground block occupies k=0

    sim = SimpleNamespace(
        dx=1.0,
        dy=1.0,
        dzt=np.array([1.0, 1.0], dtype=float),
        itot=1,
        jtot=1,
        ktot=2,
        xlen=1.0,
        ylen=1.0,
        zsize=2.0,
        zm=np.array([0.0, 1.0], dtype=float),
        Sc=solid,
        facsec={
            "c": {
                # Ground sections mapped to the fluid cell above the block.
                "facid": np.array([0, 1], dtype=int),
                "area": np.array([0.5, 0.5], dtype=float),
                "locs": np.array([[0, 0, 1], [0, 0, 1]], dtype=int),
            }
        },
        facs={"area": mesh.area_faces},
        geom=geom,
    )
    return sim, mesh

def _build_veg_against_wall_fixture():
    """2x1x1 domain: solid wall block at i=1, a vegetation cell at i=0 directly
    to its west, and a VERTICAL wall facet at x=1 (the west face of the block,
    normal -x) that shares grid cell (0,0,0) with the vegetation.

    A near-horizontal ray from a low western sun enters the vegetation cell
    (r_in = 1), attenuates over the segment up to the wall, then strikes the
    vertical facet, which the AABB overlap maps into the SAME (i,j,k) cell as
    the vegetation. This is the moller double-counting scenario for a hedge
    planted against a wall, and the horizontal-ray / vertical-facet counterpart
    of _build_veg_over_ground_block_fixture (vertical ray, horizontal facet).
    Because the sun faces the wall head-on the incidence is near-normal, so the
    scene is numerically well-conditioned despite the grazing angle to ground.
    """
    # Vertical quad at x=1 (west face of the block); winding gives normal -x
    # (pointing west, toward the sun) so the facet is lit.
    vertices = np.array(
        [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 1.0], [1.0, 0.0, 1.0]],
        dtype=float,
    )
    faces = np.array([[0, 2, 1], [0, 3, 2]], dtype=int)
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    geom = UDGeom(stl=mesh)

    solid = np.zeros((2, 1, 1), dtype=bool)
    solid[1, 0, 0] = True  # wall block occupies i=1

    sim = SimpleNamespace(
        dx=1.0,
        dy=1.0,
        dzt=np.array([1.0], dtype=float),
        itot=2,
        jtot=1,
        ktot=1,
        xlen=2.0,
        ylen=1.0,
        zsize=1.0,
        zm=np.array([0.0], dtype=float),
        Sc=solid,
        facsec={
            "c": {
                # Wall sections mapped to the fluid cell west of the block.
                "facid": np.array([0, 1], dtype=int),
                "area": np.array([0.5, 0.5], dtype=float),
                "locs": np.array([[0, 0, 0], [0, 0, 0]], dtype=int),
            }
        },
        facs={"area": mesh.area_faces},
        geom=geom,
    )
    return sim, mesh

@requires_slow_tests
class TestDirectShortwaveVegetation(unittest.TestCase):
    """Characterization + regression tests for vegetation radiation physics.

    Covers review findings P2a (sveg must be written in physical W/m3) and
    P2b (moller must attenuate a ray through vegetation before crediting the
    facet it strikes in the same cell).
    """

    LAD = 0.5
    DEC = 1.0  # extinction coefficient; tau = LAD * DEC * ds
    IRRADIANCE = 800.0

    def setUp(self):
        self.nsun_vertical = np.array([0.0, 0.0, 1.0], dtype=float)
        # tau over one full vegetation cell (dz = 1) for a vertical ray.
        self.tau = self.LAD * self.DEC * 1.0
        self.absorbed_fraction = 1.0 - np.exp(-self.tau)

    def _veg_data(self, i, j, k):
        return {
            "points": np.array([[i, j, k]], dtype=int),
            "params": {
                "lad": np.array([self.LAD], dtype=float),
                "dec": np.array([self.DEC], dtype=float),
            },
        }

    # -- P2a: sveg scaling contract ------------------------------------
    def test_facsec_sveg_is_physical_wm3(self):
        """s_veg returned by facsec must be absorbed power density [W/m3],
        i.e. veg_absorb * irradiance, consistent with bud['veg']."""
        sim, mesh = _build_veg_column_fixture()
        solver = DirectShortwaveSolver(
            sim,
            method="facsec",
            surface_mesh=mesh,
            ray_density=12.0,
            ray_jitter=0.0,
            veg_data=self._veg_data(0, 0, 0),
        )
        sdir, s_veg, bud = solver.compute(
            nsun=self.nsun_vertical, irradiance=self.IRRADIANCE, periodic_xy=False
        )

        self.assertEqual(s_veg.size, 1)
        cell_vol = sim.dx * sim.dy * sim.dzt[0]
        # Scaling contract: s_veg[k] * cell_vol == bud['veg'] (no double scale).
        self.assertAlmostEqual(
            float(s_veg[0]) * cell_vol, bud["veg"], delta=1.0e-6 * bud["veg"]
        )
        # Hand value: absorbed power per unit volume for a fully lit 1x1 column.
        expected = self.IRRADIANCE * self.absorbed_fraction / sim.dzt[0]
        self.assertAlmostEqual(float(s_veg[0]), expected, delta=0.02 * expected)
        # Budget closure: no facet here, so in == veg + out.
        self.assertAlmostEqual(
            bud["in"], bud["veg"] + bud["fac"] + bud["out"], delta=1.0e-3 * bud["in"]
        )

    def test_moller_sveg_is_physical_wm3(self):
        """s_veg returned by moller must be physical W/m3 as well."""
        sim, mesh = _build_veg_over_ground_block_fixture()
        solver = DirectShortwaveSolver(
            sim,
            method="moller",
            surface_mesh=mesh,
            ray_density=16.0,
            ray_jitter=0.0,
            veg_data=self._veg_data(0, 0, 1),
        )
        _sdir, s_veg, bud = solver.compute(
            nsun=self.nsun_vertical, irradiance=self.IRRADIANCE, periodic_xy=False
        )

        self.assertEqual(s_veg.size, 1)
        cell_vol = sim.dx * sim.dy * sim.dzt[1]
        self.assertAlmostEqual(
            float(s_veg[0]) * cell_vol, bud["veg"], delta=1.0e-6 * bud["veg"]
        )
        expected = self.IRRADIANCE * self.absorbed_fraction / sim.dzt[1]
        self.assertAlmostEqual(float(s_veg[0]), expected, delta=0.02 * expected)

    # -- P2b: moller must attenuate before crediting the facet ---------
    def test_moller_no_double_count_facet_and_veg(self):
        """With vegetation and a facet in the same cell, absorbed facet energy
        plus vegetation absorption plus outflow must equal the input.

        Before the fix the facet is credited with the pre-attenuation beam
        while the vegetation also absorbs part of it, so the sum EXCEEDS the
        input (double counting). This asserts the budget closes tightly."""
        sim, mesh = _build_veg_over_ground_block_fixture()
        solver = DirectShortwaveSolver(
            sim,
            method="moller",
            surface_mesh=mesh,
            ray_density=16.0,
            ray_jitter=0.0,
            veg_data=self._veg_data(0, 0, 1),
        )
        _sdir, _s_veg, bud = solver.compute(
            nsun=self.nsun_vertical, irradiance=self.IRRADIANCE, periodic_xy=False
        )

        total_out = bud["veg"] + bud["fac"] + bud["out"]
        self.assertAlmostEqual(bud["in"], total_out, delta=1.0e-3 * bud["in"])

    def test_moller_matches_facsec_on_veg_over_ground(self):
        """moller and facsec must agree on facet and vegetation totals for a
        tree cell sitting directly on a ground block (they disagree before the
        P2b fix, because moller credits the facet with the un-attenuated
        beam)."""
        veg = self._veg_data(0, 0, 1)

        sim_m, mesh_m = _build_veg_over_ground_block_fixture()
        moller = DirectShortwaveSolver(
            sim_m,
            method="moller",
            surface_mesh=mesh_m,
            ray_density=16.0,
            ray_jitter=0.0,
            veg_data=veg,
        )
        _sdm, sveg_m, bud_m = moller.compute(
            nsun=self.nsun_vertical, irradiance=self.IRRADIANCE, periodic_xy=False
        )

        sim_f, mesh_f = _build_veg_over_ground_block_fixture()
        facsec = DirectShortwaveSolver(
            sim_f,
            method="facsec",
            surface_mesh=mesh_f,
            ray_density=16.0,
            ray_jitter=0.0,
            veg_data=veg,
        )
        _sdf, sveg_f, bud_f = facsec.compute(
            nsun=self.nsun_vertical, irradiance=self.IRRADIANCE, periodic_xy=False
        )

        # Both budgets must close independently.
        self.assertAlmostEqual(
            bud_m["in"], bud_m["veg"] + bud_m["fac"] + bud_m["out"],
            delta=2.0e-3 * bud_m["in"],
        )
        self.assertAlmostEqual(
            bud_f["in"], bud_f["veg"] + bud_f["fac"] + bud_f["out"],
            delta=2.0e-3 * bud_f["in"],
        )
        # Facet and vegetation totals must agree between the two methods.
        self.assertAlmostEqual(bud_m["veg"], bud_f["veg"], delta=0.02 * bud_f["veg"])
        self.assertAlmostEqual(bud_m["fac"], bud_f["fac"], delta=0.02 * bud_f["fac"])
        self.assertAlmostEqual(
            float(sveg_m[0]), float(sveg_f[0]), delta=0.02 * float(sveg_f[0])
        )

    # -- P2b: near-horizontal ray hitting a VERTICAL wall through a hedge --
    def test_moller_matches_facsec_hedge_against_wall(self):
        """A hedge in the cell west of a wall, lit by a near-horizontal ray, is
        the vertical-facet counterpart of the veg-over-ground double-count test.
        moller must credit the wall with the POST-attenuation beam (after the
        hedge in the same cell has absorbed its share); before the P2b fix it
        credits the pre-attenuation beam and over-lights the wall. The fix is
        orientation-independent, so moller and facsec must agree here exactly as
        they do for the horizontal ground facet.

        Invariant asserted is kernel AGREEMENT (moller == facsec), not budget
        closure: this scene is deliberately open (a full-height wall whose top
        is exposed to the top-down ray seeding), so a fraction of the beam
        leaves the domain uncounted — but identically for both kernels, so their
        facet/vegetation totals still match. Agreement is the stronger check for
        the P2b regression, since the bug makes moller's wall credit exceed
        facsec's.
        """
        # Low western sun: near-horizontal ray (elevation ~5.7 deg, safely above
        # the |dir_z| < 1e-2 near-horizon guard) travelling east into the
        # west-facing wall at near-normal incidence.
        nsun = np.array([-1.0, 0.0, 0.1], dtype=float)
        veg = self._veg_data(0, 0, 0)

        sim_m, mesh_m = _build_veg_against_wall_fixture()
        moller = DirectShortwaveSolver(
            sim_m, method="moller", surface_mesh=mesh_m,
            ray_density=16.0, ray_jitter=0.0, veg_data=veg,
        )
        _sdm, sveg_m, bud_m = moller.compute(
            nsun=nsun, irradiance=self.IRRADIANCE, periodic_xy=False
        )

        sim_f, mesh_f = _build_veg_against_wall_fixture()
        facsec = DirectShortwaveSolver(
            sim_f, method="facsec", surface_mesh=mesh_f,
            ray_density=16.0, ray_jitter=0.0, veg_data=veg,
        )
        _sdf, sveg_f, bud_f = facsec.compute(
            nsun=nsun, irradiance=self.IRRADIANCE, periodic_xy=False
        )

        # Non-trivial scene: the vertical wall is lit through the hedge.
        self.assertGreater(bud_m["fac"], 0.0)
        self.assertGreater(bud_m["veg"], 0.0)
        # moller must not over-credit the wall vs facsec (the pre-fix symptom),
        # and the two kernels must agree on facet and vegetation absorption.
        self.assertAlmostEqual(bud_m["veg"], bud_f["veg"], delta=1.0e-3 * bud_f["veg"])
        self.assertAlmostEqual(bud_m["fac"], bud_f["fac"], delta=1.0e-3 * bud_f["fac"])
        self.assertAlmostEqual(
            float(sveg_m[0]), float(sveg_f[0]), delta=1.0e-3 * float(sveg_f[0])
        )

if __name__ == "__main__":
    unittest.main()
