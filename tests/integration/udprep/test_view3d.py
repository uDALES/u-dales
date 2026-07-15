from __future__ import annotations

import unittest

from tools.python.tests._common import REPO_ROOT, copy_case
from udprep import UDPrep
from udgeom.view3d import count_sparse_entries, resolve_view3d_exe


class TestView3D(unittest.TestCase):
    def setUp(self):
        self.temp_dir, self.case_dir = copy_case(REPO_ROOT / "examples" / "101")
        self.addCleanup(self.temp_dir.cleanup)

        exe = resolve_view3d_exe()
        if not exe.exists():
            raise unittest.SkipTest(f"View3D executable not found at {exe}")

    def test_calc_view_factors(self):
        prep = UDPrep("101", self.case_dir, load_geometry=True)
        prep.radiation.view3d_out = 0
        prep.radiation.lvfsparse = False
        legacy_vs3 = self.case_dir / "facets.vs3"
        legacy_vs3.write_text("legacy view3d input\n", encoding="ascii")

        vf, svf, paths = prep.radiation.calc_view_factors(maxD=100.0, force=True)
        vf_cached, svf_cached, cached_paths = prep.radiation.calc_view_factors(maxD=100.0, force=False)

        self.assertEqual(vf.shape[0], prep.sim.geom.stl.faces.shape[0])
        self.assertEqual(vf.shape[0], vf.shape[1])
        self.assertEqual(svf.shape[0], vf.shape[0])
        self.assertEqual(paths["vs3"].name, "facets.101.vs3")
        self.assertTrue(paths["vs3"].exists())
        self.assertFalse(legacy_vs3.exists())
        self.assertFalse(paths["vf"].exists())
        self.assertTrue(paths["vf_nc"].exists())
        self.assertTrue(paths["svf"].exists())
        self.assertIsNone(paths["vfsparse"])
        self.assertGreater(vf.nnz, 0)
        self.assertNotIn("nnz", (self.case_dir / "namoptions.101").read_text(encoding="ascii").lower())
        self.assertEqual(cached_paths["vs3"], paths["vs3"])
        self.assertEqual(cached_paths["vf"], paths["vf"])
        self.assertEqual(cached_paths["vf_nc"], paths["vf_nc"])
        self.assertEqual(cached_paths["svf"], paths["svf"])
        self.assertTrue(cached_paths["vs3"].exists())
        self.assertFalse(cached_paths["vf"].exists())
        self.assertTrue(cached_paths["vf_nc"].exists())
        self.assertTrue(cached_paths["svf"].exists())
        self.assertIs(vf_cached, vf)
        self.assertIs(svf_cached, svf)

    def test_calc_view_factors_direct_sparse_writes_nnz(self):
        prep = UDPrep("101", self.case_dir, load_geometry=True)
        prep.radiation.view3d_out = 2
        prep.radiation.lvfsparse = True

        vf, svf, paths = prep.radiation.calc_view_factors(maxD=100.0, force=True)
        vf_cached, svf_cached, cached_paths = prep.radiation.calc_view_factors(maxD=100.0, force=False)

        self.assertEqual(vf.shape[0], prep.sim.geom.stl.faces.shape[0])
        self.assertEqual(vf.shape[0], vf.shape[1])
        self.assertEqual(svf.shape[0], vf.shape[0])
        self.assertEqual(paths["vs3"].name, "facets.101.vs3")
        self.assertEqual(paths["vf"].name, "vfsparse.inp.101")
        self.assertIsNone(paths["vf_nc"])
        self.assertIsNone(paths["vfsparse"])
        self.assertTrue(paths["vs3"].exists())
        self.assertTrue(paths["vf"].exists())
        self.assertTrue(paths["svf"].exists())
        self.assertGreater(vf.nnz, 0)
        self.assertEqual(count_sparse_entries(paths["vf"]), vf.nnz)
        self.assertIn(f"nnz = {vf.nnz}", (self.case_dir / "namoptions.101").read_text(encoding="ascii"))
        self.assertEqual(cached_paths["vs3"], paths["vs3"])
        self.assertEqual(cached_paths["vf"], paths["vf"])
        self.assertEqual(cached_paths["svf"], paths["svf"])
        self.assertIs(vf_cached, vf)
        self.assertIs(svf_cached, svf)


if __name__ == "__main__":
    unittest.main()
