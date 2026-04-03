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

from udprep.udprep_ic import ICSection  # noqa: E402
from udprep.udprep_scalars import ScalarsSection  # noqa: E402


class DummySim:
    def __init__(self, path: Path, expnr="100"):
        self.path = path
        self.expnr = expnr
        self.zf = np.array([1.0, 3.0, 5.0])
        self.ktot = 3
        self.zsize = 6.0


class TestICSection(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)
        self.sim = DummySim(self.workdir, expnr="321")

    def test_generate_prof_without_lapse(self):
        section = ICSection(
            "ic",
            {"thl0": 290.0, "qt0": 0.01, "u0": 4.0, "v0": 1.0, "tke": 0.5, "lapse": 0.0},
            sim=self.sim,
            defaults={},
        )
        section.generate_prof()
        np.testing.assert_allclose(section.pr[:, 0], self.sim.zf)
        np.testing.assert_allclose(section.pr[:, 1], [290.0, 290.0, 290.0])
        np.testing.assert_allclose(section.pr[:, 2:], [[0.01, 4.0, 1.0, 0.5]] * 3)

    def test_generate_prof_with_lapse(self):
        section = ICSection(
            "ic",
            {"thl0": 300.0, "qt0": 0.0, "u0": 2.0, "v0": 0.0, "tke": 0.1, "lapse": 0.5},
            sim=self.sim,
            defaults={},
        )
        section.generate_prof()
        np.testing.assert_allclose(section.pr[:, 1], [300.0, 301.0, 302.0])

    def test_write_prof_creates_expected_file(self):
        section = ICSection(
            "ic",
            {"thl0": 290.0, "qt0": 0.01, "u0": 4.0, "v0": 1.0, "tke": 0.5, "lapse": 0.0},
            sim=self.sim,
            defaults={},
        )
        section.write_prof()
        content = (self.workdir / "prof.inp.321").read_text(encoding="ascii").splitlines()
        self.assertEqual(content[0], "# SDBL flow ")
        self.assertTrue(content[1].startswith("# z thl qt u v tke"))
        self.assertEqual(len(content), 5)


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


if __name__ == "__main__":
    unittest.main()
