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

from udprep.udprep_scalars import ScalarsSection  # noqa: E402


class DummySim:
    def __init__(self, path: Path, expnr="100"):
        self.path = path
        self.expnr = expnr
        # Non-uniform grid: faces at [0, 1, 3, 6] → centres at [0.5, 2.0, 4.5], widths [1, 2, 3]
        self.zt = np.array([0.5, 2.0, 4.5])
        self.zf = self.zt  # legacy alias, set by GridSection
        self.dzt = np.array([1.0, 2.0, 3.0])
        self.ktot = 3
        self.zsize = 6.0


class TestScalarsSection(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)
        self.sim = DummySim(self.workdir, expnr="654")

    def test_generate_scalar_uses_requested_number_of_species(self):
        section = ScalarsSection(
            "scalars",
            {
                "nsv": 3,
                "sv10": 1.0,
                "sv20": 2.0,
                "sv30": 3.0,
                "sv40": 4.0,
                "sv50": 5.0,
                "lscasrc": 0,
                "lscasrcl": 0,
                "lscasrcr": 0,
            },
            sim=self.sim,
            defaults={},
        )
        section.generate_scalar()
        self.assertEqual(section.sc.shape, (3, 4))
        np.testing.assert_allclose(section.sc[:, 1:], [[1.0, 2.0, 3.0]] * 3)

    def test_write_scalar_creates_expected_file(self):
        section = ScalarsSection(
            "scalars",
            {
                "nsv": 2,
                "sv10": 1.0,
                "sv20": 2.0,
                "sv30": 0.0,
                "sv40": 0.0,
                "sv50": 0.0,
                "lscasrc": 0,
                "lscasrcl": 0,
                "lscasrcr": 0,
            },
            sim=self.sim,
            defaults={},
        )
        section.write_scalar()
        lines = (self.workdir / "scalar.inp.654").read_text(encoding="ascii").splitlines()
        self.assertEqual(lines[0], "# SDBL flow")
        self.assertTrue(lines[1].startswith("# z scaN"))
        self.assertEqual(len(lines), 5)

    def test_generate_and_write_scalar_sources(self):
        section = ScalarsSection(
            "scalars",
            {
                "nsv": 2,
                "lscasrc": 1,
                "lscasrcl": 1,
                "lscasrcr": 0,
                "nscasrc": 1,
                "nscasrcl": 1,
                "xS": 1.0,
                "yS": 2.0,
                "zS": 3.0,
                "SSp": 4.0,
                "sigSp": 5.0,
                "xSb": 1.0,
                "ySb": 2.0,
                "zSb": 3.0,
                "xSe": 4.0,
                "ySe": 5.0,
                "zSe": 6.0,
                "SSl": 7.0,
                "sigSl": 8.0,
                "sv10": 0.0,
                "sv20": 0.0,
                "sv30": 0.0,
                "sv40": 0.0,
                "sv50": 0.0,
            },
            sim=self.sim,
            defaults={},
        )
        section.generate_scalarsources()
        np.testing.assert_allclose(section.scasrcp[0], [1.0, 2.0, 3.0, 4.0, 5.0])
        np.testing.assert_allclose(section.scasrcl[0], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
        section.write_scalarsources()
        self.assertTrue((self.workdir / "scalarsourcep.inp.1.654").exists())
        self.assertTrue((self.workdir / "scalarsourcep.inp.2.654").exists())
        self.assertTrue((self.workdir / "scalarsourcel.inp.1.654").exists())
        self.assertTrue((self.workdir / "scalarsourcel.inp.2.654").exists())

    def test_network_scalar_sources_raise(self):
        section = ScalarsSection(
            "scalars",
            {"lscasrcr": 1, "lscasrc": 0, "lscasrcl": 0, "nsv": 0},
            sim=self.sim,
            defaults={},
        )
        with self.assertRaises(ValueError):
            section.generate_scalarsources()

    # --- generate_scalar edge cases ---

    def test_generate_scalar_nsv_zero_does_nothing(self):
        section = ScalarsSection(
            "scalars",
            {"nsv": 0, "lscasrc": 0, "lscasrcl": 0, "lscasrcr": 0},
            sim=self.sim,
            defaults={},
        )
        section.generate_scalar()
        self.assertFalse(hasattr(self.sim, "sc"))

    def test_generate_scalar_nsv_too_large_raises(self):
        section = ScalarsSection(
            "scalars",
            {
                "nsv": 6,
                "sv10": 0.0, "sv20": 0.0, "sv30": 0.0, "sv40": 0.0, "sv50": 0.0,
                "lscasrc": 0, "lscasrcl": 0, "lscasrcr": 0,
            },
            sim=self.sim,
            defaults={},
        )
        with self.assertRaises(ValueError):
            section.generate_scalar()

    # --- write_scalar edge cases ---

    def test_write_scalar_nsv_zero_writes_no_file(self):
        section = ScalarsSection(
            "scalars",
            {"nsv": 0, "lscasrc": 0, "lscasrcl": 0, "lscasrcr": 0},
            sim=self.sim,
            defaults={},
        )
        section.write_scalar()
        self.assertFalse((self.workdir / "scalar.inp.654").exists())

    # --- generate_scalarsources validation branches ---

    def test_generate_scalarsources_nsv_zero_with_lscasrc_raises(self):
        section = ScalarsSection(
            "scalars",
            {"nsv": 0, "lscasrc": 1, "lscasrcl": 0, "lscasrcr": 0, "nscasrc": 1},
            sim=self.sim,
            defaults={},
        )
        with self.assertRaises(ValueError):
            section.generate_scalarsources()

    def test_generate_scalarsources_zero_nscasrc_raises(self):
        section = ScalarsSection(
            "scalars",
            {"nsv": 1, "lscasrc": 1, "lscasrcl": 0, "lscasrcr": 0, "nscasrc": 0},
            sim=self.sim,
            defaults={},
        )
        with self.assertRaises(ValueError):
            section.generate_scalarsources()

    def test_generate_scalarsources_missing_point_coords_raises(self):
        # nscasrc==1 but coordinates left at sentinel -1.0
        section = ScalarsSection(
            "scalars",
            {
                "nsv": 1, "lscasrc": 1, "lscasrcl": 0, "lscasrcr": 0,
                "nscasrc": 1,
                "xS": -1.0, "yS": -1.0, "zS": -1.0, "SSp": -1.0, "sigSp": -1.0,
            },
            sim=self.sim,
            defaults={},
        )
        with self.assertRaises(ValueError):
            section.generate_scalarsources()

    def test_generate_scalarsources_zero_nscasrcl_raises(self):
        section = ScalarsSection(
            "scalars",
            {"nsv": 1, "lscasrc": 0, "lscasrcl": 1, "lscasrcr": 0, "nscasrcl": 0},
            sim=self.sim,
            defaults={},
        )
        with self.assertRaises(ValueError):
            section.generate_scalarsources()

    def test_generate_scalarsources_missing_line_coords_raises(self):
        # nscasrcl==1 but coordinates left at sentinel -1.0
        section = ScalarsSection(
            "scalars",
            {
                "nsv": 1, "lscasrc": 0, "lscasrcl": 1, "lscasrcr": 0,
                "nscasrcl": 1,
                "xSb": -1.0, "ySb": -1.0, "zSb": -1.0,
                "xSe": -1.0, "ySe": -1.0, "zSe": -1.0,
                "SSl": -1.0, "sigSl": -1.0,
            },
            sim=self.sim,
            defaults={},
        )
        with self.assertRaises(ValueError):
            section.generate_scalarsources()

    # --- write_scalarsources file-exists guard ---

    def test_write_scalarsources_existing_file_warns_and_skips(self):
        section = ScalarsSection(
            "scalars",
            {
                "nsv": 1, "lscasrc": 1, "lscasrcl": 0, "lscasrcr": 0,
                "nscasrc": 1,
                "xS": 1.0, "yS": 2.0, "zS": 3.0, "SSp": 4.0, "sigSp": 5.0,
                "sv10": 0.0, "sv20": 0.0, "sv30": 0.0, "sv40": 0.0, "sv50": 0.0,
            },
            sim=self.sim,
            defaults={},
        )
        section.generate_scalarsources()
        # Pre-create the output file so the guard triggers on the second write
        existing = self.workdir / "scalarsourcep.inp.1.654"
        existing.write_text("existing content", encoding="ascii")
        with self.assertWarns(UserWarning):
            section.write_scalarsources()
        # File should be unchanged (not overwritten)
        self.assertEqual(existing.read_text(encoding="ascii"), "existing content")


if __name__ == "__main__":
    unittest.main()
