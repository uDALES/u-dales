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


if __name__ == "__main__":
    unittest.main()
