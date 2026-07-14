"""Characterization tests for pure numeric helpers in the radiation section.

These pin the current behaviour of the IO-free compute helpers so they can later
be extracted from the IO/subprocess-heavy RadiationSection without changing
results (task 7 prerequisite).
"""
import sys
import unittest
from pathlib import Path

import numpy as np

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from _common import PYTHON_DIR  # noqa: E402

if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from udprep.udprep_radiation import RadiationSection  # noqa: E402


class TestInterpMakima(unittest.TestCase):
    """Modified-Akima interpolation (RadiationSection._interp_makima)."""

    interp = staticmethod(RadiationSection._interp_makima)

    def test_passes_through_data_points(self):
        x = np.array([0.0, 1, 2, 3, 4])
        y = np.array([0.0, 1, 4, 9, 16])
        np.testing.assert_allclose(self.interp(x, y, x), y)

    def test_reproduces_linear_data_exactly(self):
        x = np.array([0.0, 1, 2])
        y = np.array([0.0, 2, 4])
        np.testing.assert_allclose(self.interp(x, y, np.array([0.5, 1.5])), [1.0, 3.0])

    def test_known_values_on_quadratic_samples(self):
        x = np.array([0.0, 1, 2, 3, 4])
        y = np.array([0.0, 1, 4, 9, 16])
        np.testing.assert_allclose(
            self.interp(x, y, np.array([0.5, 2.5])), [0.3125, 6.239583], rtol=1e-6
        )

    def test_rejects_degenerate_input(self):
        with self.assertRaises(ValueError):
            self.interp(np.array([0.0]), np.array([0.0]), np.array([0.0]))
        with self.assertRaises(ValueError):
            self.interp(np.array([0.0, 0.0]), np.array([1.0, 2.0]), np.array([0.0]))


class TestNetShortwaveReflections(unittest.TestCase):
    """Iterative multi-bounce net-shortwave physics (calc_reflections_sw).

    Ground-truth values captured from the current implementation; also checks the
    physical invariants (no-reflection limit, energy split) so the numerics are
    pinned before the compute is separated from the section's IO (task 7).
    """

    @staticmethod
    def _reflect(*args, **kwargs):
        # self is unused by calc_reflections_sw; a bare instance is enough.
        section = object.__new__(RadiationSection)
        return section.calc_reflections_sw(*args, **kwargs)

    def test_no_view_factors_reduces_to_absorbed_direct_plus_diffuse(self):
        # With vf=0 there are no reflections: knet = (1-albedo)*(sdir + dsky*svf).
        knet = self._reflect(
            np.array([100.0, 200.0]), 0.0, np.zeros((2, 2)),
            np.array([1.0, 1.0]), np.array([0.5, 0.5]),
        )
        np.testing.assert_allclose(knet, [50.0, 100.0])

    def test_mutual_reflection_between_two_facets(self):
        vf = np.array([[0.0, 0.5], [0.5, 0.0]])
        knet = self._reflect(
            np.array([100.0, 0.0]), 0.0, vf,
            np.array([0.5, 0.5]), np.array([0.5, 0.5]),
        )
        np.testing.assert_allclose(knet, [53.320312, 13.28125], rtol=1e-6)

    def test_diffuse_sky_contribution(self):
        vf = np.array([[0.0, 0.5], [0.5, 0.0]])
        knet = self._reflect(
            np.array([50.0, 50.0]), 100.0, vf,
            np.array([0.3, 0.7]), np.array([0.2, 0.2]),
        )
        np.testing.assert_allclose(knet, [74.24, 103.36], rtol=1e-6)

    def test_zero_albedo_absorbs_all_incoming(self):
        # albedo=0 -> nothing reflected -> knet == sdir + dsky*svf, one iteration.
        knet = self._reflect(
            np.array([80.0, 120.0]), 50.0, np.array([[0.0, 1.0], [1.0, 0.0]]),
            np.array([0.5, 0.5]), np.array([0.0, 0.0]),
        )
        np.testing.assert_allclose(knet, [80.0 + 50.0 * 0.5, 120.0 + 50.0 * 0.5])

    def test_validation(self):
        with self.assertRaises(ValueError):
            self._reflect(np.array([1.0, 2.0]), 0.0, np.zeros((2, 2)),
                          np.array([1.0]), np.array([0.5, 0.5]))  # shape mismatch
        with self.assertRaises(ValueError):
            self._reflect(np.array([1.0]), 0.0, np.zeros((1, 1)),
                          np.array([1.0]), np.array([0.5]), tol=0.0)


if __name__ == "__main__":
    unittest.main()
