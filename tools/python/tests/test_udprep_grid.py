import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from _common import PYTHON_DIR

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udprep.udprep_grid import GridSection  # noqa: E402


def make_grid_section(**overrides):
    values = {
        "xlen": 8.0,
        "ylen": 6.0,
        "zsize": 10.0,
        "itot": 4,
        "jtot": 3,
        "ktot": 5,
        "lzstretch": 0,
        "stretchconst": 2.0,
        "lstretchexp": 0,
        "lstretchexpcheck": 0,
        "lstretchtanh": 0,
        "lstretch2tanh": 0,
        "hlin": 2.0,
    }
    values.update(overrides)
    sim = SimpleNamespace()
    return GridSection("grid", values, sim=sim, defaults={}), sim


class TestGridSection(unittest.TestCase):
    def test_generate_xygrid_populates_staggered_and_face_grids(self):
        section, sim = make_grid_section()
        section._refresh_derived_grid_params()
        section.generate_xygrid()

        np.testing.assert_allclose(section.xt, [1.0, 3.0, 5.0, 7.0])
        np.testing.assert_allclose(section.yt, [1.0, 3.0, 5.0])
        np.testing.assert_allclose(section.xm, [0.0, 2.0, 4.0, 6.0, 8.0])
        np.testing.assert_allclose(section.ym, [0.0, 2.0, 4.0, 6.0])
        np.testing.assert_allclose(sim.xf, section.xt)
        np.testing.assert_allclose(sim.xh, section.xm)

    def test_generate_zgrid_uniform_matches_expected_faces_and_centers(self):
        section, sim = make_grid_section()
        section._refresh_derived_grid_params()
        section.generate_zgrid()

        np.testing.assert_allclose(section.zm, [0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
        np.testing.assert_allclose(section.zt, [1.0, 3.0, 5.0, 7.0, 9.0])
        np.testing.assert_allclose(section.dzt, [2.0, 2.0, 2.0, 2.0, 2.0])
        np.testing.assert_allclose(sim.zh, section.zm)
        np.testing.assert_allclose(sim.zf, section.zt)

    def test_refresh_derived_params_recovers_invalid_hlin(self):
        section, _ = make_grid_section(hlin="invalid")
        section._refresh_derived_grid_params()
        self.assertEqual(section.dx, 2.0)
        self.assertEqual(section.dy, 2.0)
        self.assertEqual(section.dz, 2.0)
        self.assertEqual(section.hlin, 1.0)

    def test_stretch_exp_check_builds_monotone_face_grid(self):
        section, _ = make_grid_section(
            lzstretch=1,
            lstretchexpcheck=1,
            zsize=12.0,
            ktot=6,
            hlin=2.0,
        )
        section._refresh_derived_grid_params()
        section.generate_zgrid()

        self.assertEqual(len(section.zm), 7)
        self.assertEqual(len(section.zt), 6)
        self.assertTrue(np.all(np.diff(section.zm) > 0.0))
        np.testing.assert_allclose(section.zt, 0.5 * (section.zm[:-1] + section.zm[1:]))
        np.testing.assert_allclose(section.dzt, np.diff(section.zm))

    def test_tanh_stretch_uses_computational_mapping_and_keeps_linear_prefix(self):
        section, _ = make_grid_section(
            lzstretch=1,
            lstretchtanh=1,
            zsize=20.0,
            ktot=10,
            hlin=2.0,
        )
        section._refresh_derived_grid_params()
        section.generate_zgrid()

        np.testing.assert_allclose(section.zm[:2], [0.0, 2.0])
        self.assertTrue(np.all(np.diff(section.zm) > 0.0))
        self.assertAlmostEqual(section.zm[-1], 20.0)

    def test_invalid_stretch_configuration_raises(self):
        section, _ = make_grid_section(lzstretch=1)
        section._refresh_derived_grid_params()
        with self.assertRaises(ValueError):
            section.generate_zgrid()


if __name__ == "__main__":
    unittest.main()
