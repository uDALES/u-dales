from __future__ import annotations

from types import SimpleNamespace
import unittest

import numpy as np
import trimesh

try:
    from _common import REPO_ROOT
except ModuleNotFoundError:
    from tools.python.tests._common import REPO_ROOT

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
        zf=np.array([0.0, 1.0], dtype=float),
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


if __name__ == "__main__":
    unittest.main()
