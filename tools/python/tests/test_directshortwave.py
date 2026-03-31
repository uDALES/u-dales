from __future__ import annotations

import unittest

import numpy as np

from _common import REPO_ROOT
from udbase import UDBase
from udprep.directshortwave import DirectShortwaveSolver


class TestDirectShortwave(unittest.TestCase):
    def setUp(self):
        self.sim = UDBase("101", REPO_ROOT / "examples" / "101", load_geometry=True)
        azimuth_deg = 20.0
        elevation_deg = 15.0
        az = np.deg2rad(azimuth_deg)
        el = np.deg2rad(elevation_deg)
        self.nsun = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)], dtype=float)

    def test_scanline_matches_facsec_shape_and_scale(self):
        facsec = DirectShortwaveSolver(self.sim, method="facsec", ray_density=2.0, ray_jitter=0.0)
        scanline = DirectShortwaveSolver(self.sim, method="scanline", ray_density=2.0, ray_jitter=0.0)

        sdir_facsec, _, bud_facsec = facsec.compute(nsun=self.nsun, irradiance=800.0, periodic_xy=False)
        sdir_scanline, veg_scanline, bud_scanline = scanline.compute(
            nsun=self.nsun,
            irradiance=800.0,
            periodic_xy=False,
        )

        self.assertEqual(sdir_facsec.shape, sdir_scanline.shape)
        self.assertEqual(veg_scanline.size, 0)
        self.assertTrue(np.all(np.isfinite(sdir_facsec)))
        self.assertTrue(np.all(np.isfinite(sdir_scanline)))
        self.assertGreater(bud_facsec["fac"], 0.0)
        self.assertGreater(bud_scanline["fac"], 0.0)

        common = (sdir_facsec > 1.0e-6) & (sdir_scanline > 1.0e-6)
        self.assertTrue(np.any(common))
        rel = np.abs(sdir_scanline[common] - sdir_facsec[common]) / np.maximum(sdir_facsec[common], 1.0)
        corr = float(np.corrcoef(sdir_facsec[common], sdir_scanline[common])[0, 1])
        energy_ratio = float(bud_scanline["fac"] / bud_facsec["fac"])

        self.assertLess(float(np.median(rel)), 0.5)
        self.assertLess(float(np.quantile(rel, 0.9)), 1.0)
        self.assertGreater(corr, 0.85)
        self.assertGreater(energy_ratio, 0.65)
        self.assertLess(energy_ratio, 0.80)


if __name__ == "__main__":
    unittest.main()
