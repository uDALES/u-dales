from __future__ import annotations

import unittest

import numpy as np

from _common import REPO_ROOT
from udbase import UDBase
from udprep.directshortwave import DirectShortwaveSolver


class TestDirectShortwavePeriodic(unittest.TestCase):
    def setUp(self):
        self.sim = UDBase("101", REPO_ROOT / "examples" / "101", load_geometry=True)
        azimuth_deg = 20.0
        elevation_deg = 15.0
        az = np.deg2rad(azimuth_deg)
        el = np.deg2rad(elevation_deg)
        self.nsun = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)], dtype=float)

    def test_periodic_mode_runs_and_returns_finite_fluxes(self):
        solver = DirectShortwaveSolver(self.sim, method="moller", ray_density=2.0, ray_jitter=0.0)

        sdir_nonperiodic, _, bud_nonperiodic = solver.compute(
            nsun=self.nsun,
            irradiance=800.0,
            periodic_xy=False,
        )
        sdir_periodic, _, bud_periodic = solver.compute(
            nsun=self.nsun,
            irradiance=800.0,
            periodic_xy=True,
        )

        self.assertEqual(sdir_nonperiodic.shape, sdir_periodic.shape)
        self.assertTrue(np.all(np.isfinite(sdir_nonperiodic)))
        self.assertTrue(np.all(np.isfinite(sdir_periodic)))
        self.assertGreater(bud_nonperiodic["fac"], 0.0)
        self.assertGreater(bud_periodic["fac"], 0.0)
        self.assertGreater(np.count_nonzero(sdir_periodic > 0.0), 0)
        self.assertFalse(np.allclose(sdir_nonperiodic, sdir_periodic))


if __name__ == "__main__":
    unittest.main()
