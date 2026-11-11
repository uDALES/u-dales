from __future__ import annotations

import os
import subprocess
import unittest
from pathlib import Path

from tools.python.tests._common import REPO_ROOT, copy_case
from udprep import UDPrep
from udgeom.view3d import (
    count_sparse_entries,
    read_view3d_output,
    resolve_view3d_exe,
)


DATA_DIR = Path(__file__).resolve().parent / "data"


class TestView3D(unittest.TestCase):
    def setUp(self):
        self._old_preproc_ncpu = os.environ.get("PREPROC_NCPU")
        if self._old_preproc_ncpu is None:
            os.environ["PREPROC_NCPU"] = "8"
        self.addCleanup(self._restore_preproc_ncpu)

        self.temp_dir, self.case_dir = copy_case(REPO_ROOT / "examples" / "101")
        self.addCleanup(self.temp_dir.cleanup)

        exe = resolve_view3d_exe()
        if not exe.exists():
            raise unittest.SkipTest(f"View3D executable not found at {exe}")

    def _restore_preproc_ncpu(self):
        if self._old_preproc_ncpu is None:
            os.environ.pop("PREPROC_NCPU", None)
        else:
            os.environ["PREPROC_NCPU"] = self._old_preproc_ncpu

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

    def test_adjacent_obstructor_does_not_increase_view_factor(self):
        source = DATA_DIR / "view3d_adjacent_obstructor.vs3"
        obstructed_input = self.case_dir / source.name
        unobstructed_input = self.case_dir / "view3d_unobstructed.vs3"
        obstructed_input.write_bytes(source.read_bytes())
        unobstructed_input.write_text(
            "\n".join(
                line
                for line in source.read_text(encoding="ascii").splitlines()
                if not line.startswith("O 3 ")
            )
            + "\n",
            encoding="ascii",
        )

        exe = resolve_view3d_exe()
        obstructed_output = self.case_dir / "view3d_obstructed.out"
        unobstructed_output = self.case_dir / "view3d_unobstructed.out"
        env = os.environ.copy()
        env["VIEW3D_NUM_THREADS"] = "1"
        env["OMP_NUM_THREADS"] = "1"
        for input_path, output_path in (
            (obstructed_input, obstructed_output),
            (unobstructed_input, unobstructed_output),
        ):
            subprocess.run(
                [str(exe), str(input_path), str(output_path)],
                cwd=self.case_dir,
                env=env,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

        obstructed = read_view3d_output(
            obstructed_output, nfacets=2, outformat=2
        )
        unobstructed = read_view3d_output(
            unobstructed_output, nfacets=2, outformat=2
        )
        self.assertLessEqual(
            float(obstructed[0, 1]),
            float(unobstructed[0, 1]) + 1.0e-6,
        )
        self.assertLess(
            float(obstructed[0, 1]),
            float(unobstructed[0, 1]) - 1.0e-3,
        )


if __name__ == "__main__":
    unittest.main()
