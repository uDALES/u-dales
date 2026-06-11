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
        "dzlin": 1.0,
    }
    values.update(overrides)
    sim = SimpleNamespace()
    return GridSection("grid", values, sim=sim, defaults={}), sim


class TestGridSection(unittest.TestCase):
    def test_run_all_populates_complete_uniform_grid_state(self):
        section, sim = make_grid_section()

        section.run_all()

        self.assertEqual(section.dx, 2.0)
        self.assertEqual(section.dy, 2.0)
        self.assertEqual(section.dz, 2.0)
        np.testing.assert_allclose(sim.xt, [1.0, 3.0, 5.0, 7.0])
        np.testing.assert_allclose(sim.yt, [1.0, 3.0, 5.0])
        np.testing.assert_allclose(sim.xm, [0.0, 2.0, 4.0, 6.0])
        np.testing.assert_allclose(sim.ym, [0.0, 2.0, 4.0])
        np.testing.assert_allclose(sim.zm, [0.0, 2.0, 4.0, 6.0, 8.0])
        np.testing.assert_allclose(sim.zt, [1.0, 3.0, 5.0, 7.0, 9.0])
        np.testing.assert_allclose(sim.dzt, [2.0, 2.0, 2.0, 2.0, 2.0])
        np.testing.assert_allclose(sim.zf, sim.zt)
        np.testing.assert_allclose(sim.zh, sim.zm)

    def test_generate_xygrid_populates_staggered_and_face_grids(self):
        section, sim = make_grid_section()
        section._refresh_derived_grid_params()
        section.generate_xygrid()

        np.testing.assert_allclose(section.xt, [1.0, 3.0, 5.0, 7.0])
        np.testing.assert_allclose(section.yt, [1.0, 3.0, 5.0])
        np.testing.assert_allclose(section.xm, [0.0, 2.0, 4.0, 6.0])
        np.testing.assert_allclose(section.ym, [0.0, 2.0, 4.0])
        np.testing.assert_allclose(sim.xt, section.xt)
        np.testing.assert_allclose(sim.xm, section.xm)

    def test_generate_zgrid_uniform_matches_expected_faces_and_centers(self):
        section, sim = make_grid_section()
        section._refresh_derived_grid_params()
        section.generate_zgrid()

        np.testing.assert_allclose(section.zm, [0.0, 2.0, 4.0, 6.0, 8.0])
        np.testing.assert_allclose(section.zt, [1.0, 3.0, 5.0, 7.0, 9.0])
        np.testing.assert_allclose(section.dzt, [2.0, 2.0, 2.0, 2.0, 2.0])
        np.testing.assert_allclose(sim.zm, section.zm)
        np.testing.assert_allclose(sim.zt, section.zt)

    def test_refresh_derived_params_uses_computed_hlin_when_not_set(self):
        section, _ = make_grid_section(hlin=None, lzstretch=True, lstretchexp=True)
        section._refresh_derived_grid_params()
        self.assertEqual(section.dx, 2.0)
        self.assertEqual(section.dy, 2.0)
        self.assertEqual(section.dz, 2.0)
        self.assertEqual(section.hlin, 1.0)

    def test_stretch_exp_check_builds_monotone_face_grid(self):
        # ktot=20 / dzlin=0.5 gives alpha≈0.44 → cell stretch ratio≈1.028,
        # within the 0.95–1.05 quality band, so no RuntimeWarning fires.
        section, _ = make_grid_section(
            lzstretch=1,
            lstretchexpcheck=1,
            zsize=12.0,
            ktot=20,
            hlin=2.0,
            dzlin=0.5,
        )
        section._refresh_derived_grid_params()
        section.generate_zgrid()

        self.assertEqual(len(section.zm), 20)
        self.assertEqual(len(section.zt), 20)
        self.assertTrue(np.all(np.diff(section.zm) > 0.0))
        np.testing.assert_allclose(section.zt, section.zm + 0.5 * section.dzt)
        np.testing.assert_allclose(section.dzt, np.diff(np.append(section.zm, 12.0)))

    def test_stretch_exp_check_avoids_zero_alpha_root_for_case_993_grid(self):
        section, _ = make_grid_section(
            lzstretch=1,
            lstretchexpcheck=1,
            zsize=100.0,
            ktot=128,
            hlin=25.0,
            dzlin=0.4,
        )
        section._refresh_derived_grid_params()
        section.generate_zgrid()

        self.assertEqual(len(section.zm), 128)
        self.assertEqual(len(section.zt), 128)
        self.assertTrue(np.all(np.isfinite(section.zt)))
        self.assertTrue(np.all(np.diff(section.zm) > 0.0))
        np.testing.assert_allclose(section.zt[:63], 0.2 + 0.4 * np.arange(63))
        self.assertAlmostEqual(float(section.zt[63]), 25.402856293747867)
        self.assertAlmostEqual(float(section.zm[-1] + section.dzt[-1]), 100.0)

    def test_tanh_stretch_uses_computational_mapping_and_keeps_linear_prefix(self):
        # hlin=3.0 / dzlin=1.5 → loop converges to gf≈0.85,
        # last spacing ≈ 2.6 m ≤ 3×dzlin=4.5 m, so no RuntimeWarning fires.
        section, _ = make_grid_section(
            lzstretch=1,
            lstretchtanh=1,
            zsize=20.0,
            ktot=10,
            hlin=3.0,
            dzlin=1.5,
        )
        section._refresh_derived_grid_params()
        section.generate_zgrid()

        np.testing.assert_allclose(section.zm[:2], [0.0, 1.5])
        self.assertTrue(np.all(np.diff(section.zm) > 0.0))
        self.assertAlmostEqual(float(section.zm[-1] + section.dzt[-1]), 20.0)

    def test_invalid_stretch_configuration_raises(self):
        section, _ = make_grid_section(lzstretch=1)
        section._refresh_derived_grid_params()
        with self.assertRaises(ValueError):
            section.generate_zgrid()

    def test_multiple_stretch_methods_raise(self):
        section, _ = make_grid_section(
            lzstretch=1,
            lstretchexp=1,
            lstretchtanh=1,
        )
        section._refresh_derived_grid_params()
        with self.assertRaisesRegex(ValueError, "multiple stretch methods"):
            section.generate_zgrid()

    def test_stretch_exp_builds_monotone_face_grid_with_linear_prefix(self):
        # zsize=20 / ktot=10 / hlin=3 / dzlin=1.5 → loop converges to gf≈0.73,
        # last spacing ≈ 2.9 m ≤ 3×dzlin=4.5 m, so no RuntimeWarning fires.
        section, _ = make_grid_section(
            lzstretch=1,
            lstretchexp=1,
            zsize=20.0,
            ktot=10,
            hlin=3.0,
            dzlin=1.5,
            stretchconst=2.0,
        )
        section._refresh_derived_grid_params()
        section.generate_zgrid()

        self.assertEqual(len(section.zm), 10)
        self.assertEqual(len(section.zt), 10)
        # Linear prefix: first two faces should be 0.0 and dzlin
        np.testing.assert_allclose(section.zm[:2], [0.0, 1.5])
        # Strictly monotone
        self.assertTrue(np.all(np.diff(section.zm) > 0.0))
        # Centers lie at mid-points of their cells
        np.testing.assert_allclose(section.zt, section.zm + 0.5 * section.dzt)
        # dzt reconstructed from faces spans up to zsize
        np.testing.assert_allclose(section.dzt, np.diff(np.append(section.zm, 20.0)))

    def test_stretch_2tanh_builds_monotone_face_grid_with_linear_prefix(self):
        section, _ = make_grid_section(
            lzstretch=1,
            lstretch2tanh=1,
            zsize=20.0,
            ktot=10,
            hlin=2.0,
            dzlin=1.0,
            stretchconst=2.0,
        )
        section._refresh_derived_grid_params()
        section.generate_zgrid()

        self.assertEqual(len(section.zm), 10)
        self.assertEqual(len(section.zt), 10)
        # Linear prefix preserved
        np.testing.assert_allclose(section.zm[:2], [0.0, 1.0])
        # Strictly monotone
        self.assertTrue(np.all(np.diff(section.zm) > 0.0))
        # Domain top is respected
        self.assertAlmostEqual(float(section.zm[-1] + section.dzt[-1]), 20.0)

    def test_non_expcheck_stretch_methods_keep_hlin_as_transition_face(self):
        methods = ("lstretchexp", "lstretchtanh", "lstretch2tanh")

        for method in methods:
            with self.subTest(method=method):
                section, _ = make_grid_section(
                    lzstretch=1,
                    zsize=100.0,
                    ktot=128,
                    hlin=25.0,
                    dzlin=0.4,
                    stretchconst=2.0,
                    **{method: 1},
                )
                section._refresh_derived_grid_params()
                section.generate_zgrid()

                il = GridSection._matlab_round_int(25.0 / 0.4)
                linear_faces = np.append(section.zm[:il], 25.0)
                np.testing.assert_allclose(linear_faces, np.linspace(0.0, 25.0, il + 1))
                self.assertLess(section.zt[il - 1], 25.0)
                self.assertGreater(section.zt[il], 25.0)
                self.assertTrue(np.all(np.diff(section.zm) > 0.0))
                self.assertAlmostEqual(float(section.zm[-1] + section.dzt[-1]), 100.0)

    def test_warn_large_top_spacing_emits_runtime_warning(self):
        section, _ = make_grid_section(
            lzstretch=1,
            lstretchexp=1,
            # Very large domain with few cells forces a large top spacing
            zsize=100.0,
            ktot=6,
            hlin=2.0,
            dzlin=1.0,
            stretchconst=4.0,
        )
        section._refresh_derived_grid_params()
        with self.assertWarns(RuntimeWarning):
            section.generate_zgrid()


if __name__ == "__main__":
    unittest.main()
