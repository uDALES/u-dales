import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from _common import PYTHON_DIR

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udbase import UDBase  # noqa: E402
from udprep.udprep_vegetation import VegetationSection  # noqa: E402


class DummySim:
    save_veg = UDBase.save_veg

    def __init__(self, path: Path, expnr="777"):
        self.path = path
        self.expnr = expnr
        self.itot = 8
        self.jtot = 8
        self.ktot = 6
        self._lfgeom = False
        self.geom = None


class TestVegetationSection(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)
        self.sim = DummySim(self.workdir)

    def _make_section(self, **overrides):
        values = {
            "treesfile": f"trees.inp.{self.sim.expnr}",
            "lad": 1.0,
            "cd": 0.3,
            "ud": 0.0002,
            "dec": 0.0,
            "lsize": 0.15,
            "r_s": 0.0,
            "ltrees": True,
            "ltreesfile": True,
        }
        values.update(overrides)
        return VegetationSection("vegetation", values, sim=self.sim, defaults={})

    def test_load_block_populates_sparse_points_without_geometry(self):
        (self.workdir / "trees.inp.777").write_text("1 2 3 4 1 2\n", encoding="ascii")
        section = self._make_section()

        result = section.load_block()

        self.assertEqual(result, self.workdir / "trees.inp.777")
        self.assertEqual(section.ntrees, 8)
        np.testing.assert_array_equal(section.veg["points"][0], [0, 2, 0])
        self.assertEqual(section.veg["params"]["id"][0], 1)
        self.assertTrue(np.all(section.veg["params"]["lad"] == 1.0))

    def test_block_to_veg_writes_sparse_files(self):
        (self.workdir / "trees.inp.777").write_text("1 2 1 1 1 1\n", encoding="ascii")
        section = self._make_section()

        out = section.block_to_veg()

        veg_path = self.workdir / "veg.inp.777"
        params_path = self.workdir / "veg_params.inp.777"
        self.assertEqual(out["veg"], veg_path)
        self.assertEqual(out["params"], params_path)
        self.assertTrue(veg_path.exists())
        self.assertTrue(params_path.exists())
        veg_lines = veg_path.read_text(encoding="ascii").splitlines()
        params_lines = params_path.read_text(encoding="ascii").splitlines()
        self.assertEqual(veg_lines[0], "# position (i,j,k)")
        self.assertEqual(params_lines[0], "# id lad cd ud dec lsize r_s")
        self.assertEqual(len(veg_lines), 3)
        self.assertEqual(len(params_lines), 3)

    def test_load_block_rejects_out_of_bounds_block(self):
        (self.workdir / "trees.inp.777").write_text("1 9 1 1 1 1\n", encoding="ascii")
        section = self._make_section()
        with self.assertRaises(ValueError):
            section.load_block()

    def test_load_block_rejects_missing_tree_file(self):
        section = self._make_section()
        with self.assertRaises(FileNotFoundError):
            section.load_block()

    def test_load_block_accepts_legacy_seven_column_format(self):
        (self.workdir / "trees.inp.777").write_text("5 2 1 3 3 2 1\n", encoding="ascii")
        section = self._make_section()
        section.load_block()
        self.assertEqual(section.ntrees, 4)
        self.assertTrue(np.all(section.veg["params"]["id"] == 1))


if __name__ == "__main__":
    unittest.main()
