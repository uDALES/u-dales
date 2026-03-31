from __future__ import annotations

import unittest

from _common import REPO_ROOT, copy_case
from udprep import UDPrep


class TestView3D(unittest.TestCase):
    def setUp(self):
        self.temp_dir, self.case_dir = copy_case(REPO_ROOT / "examples" / "101")
        self.addCleanup(self.temp_dir.cleanup)

    def test_calc_view_factors(self):
        prep = UDPrep("101", self.case_dir, load_geometry=True)
        prep.radiation.view3d_out = 0

        vf, svf, paths = prep.radiation.calc_view_factors(maxD=100.0, force=True)

        self.assertEqual(vf.shape[0], prep.sim.geom.stl.faces.shape[0])
        self.assertEqual(vf.shape[0], vf.shape[1])
        self.assertEqual(svf.shape[0], vf.shape[0])
        self.assertTrue(paths["vs3"].exists())
        self.assertTrue(paths["vf"].exists())
        self.assertTrue(paths["svf"].exists())
        self.assertGreater(vf.nnz, 0)


if __name__ == "__main__":
    unittest.main()
