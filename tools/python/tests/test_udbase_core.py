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
from udvis import UDVis  # noqa: E402


class TestUDBaseCore(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)
        (self.workdir / "namoptions.001").write_text(
            "\n".join(
                [
                    "&DOMAIN",
                    " itot = 4",
                    " jtot = 3",
                    " ktot = 2",
                    " xlen = 40.0",
                    " ylen = 30.0",
                    " zsize = 20.0",
                    "/",
                ]
            )
            + "\n",
            encoding="ascii",
        )

    def test_constructor_initializes_visualization_facade(self):
        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

        self.assertIsInstance(sim.vis, UDVis)
        self.assertIs(sim.vis.sim, sim)

    def test_load_facsec_applies_one_based_to_zero_based_mapping(self):
        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)
        sim.ffacet_sections = "facet_sections"
        sim.ffluid_boundary = "fluid_boundary"
        (self.workdir / "facet_sections_c.txt").write_text(
            "\n".join(
                [
                    "# facid area fluid_idx distance",
                    "1 2.5 2 0.10",
                    "3 1.0 1 0.20",
                ]
            )
            + "\n",
            encoding="ascii",
        )
        (self.workdir / "fluid_boundary_c.txt").write_text(
            "\n".join(
                [
                    "# i j k",
                    "1 2 1",
                    "4 3 2",
                ]
            )
            + "\n",
            encoding="ascii",
        )

        facsec = sim.load_facsec("c")

        np.testing.assert_array_equal(facsec["facid"], np.array([0, 2]))
        np.testing.assert_allclose(facsec["area"], np.array([2.5, 1.0]))
        np.testing.assert_array_equal(facsec["locs"], np.array([[3, 2, 1], [0, 1, 0]]))
        np.testing.assert_allclose(facsec["distance"], np.array([0.10, 0.20]))

    def test_convert_fac_to_field_accumulates_cell_density(self):
        sim = UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)
        facsec = {
            "facid": np.array([0, 1, 0]),
            "area": np.array([2.0, 1.0, 3.0]),
            "locs": np.array([[0, 0, 0], [0, 0, 1], [0, 0, 0]]),
        }
        var = np.array([10.0, 20.0])

        field = sim.convert_fac_to_field(var, facsec=facsec)

        self.assertEqual(field.shape, (4, 3, 2))
        self.assertAlmostEqual(field[0, 0, 0], 0.05)
        self.assertAlmostEqual(field[0, 0, 1], 0.02)
        self.assertTrue(np.all(field[1:, :, :] == 0.0))
