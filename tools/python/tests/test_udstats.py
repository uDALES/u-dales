"""Characterization tests for the extracted udstats helpers.

These pin the current numerical behaviour of merge_stat / time_average /
coarsegrain_field (values captured from the pre-extraction implementation), and
confirm the UDBase.* static-method wrappers still delegate to them.
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

import udstats  # noqa: E402
from udbase import UDBase  # noqa: E402


class TestMergeStat(unittest.TestCase):
    def test_mean_only_windows_the_last_axis(self):
        np.testing.assert_allclose(udstats.merge_stat(np.array([1.0, 2, 3, 4]), 2), [1.5, 3.5])

    def test_discards_oldest_incomplete_window(self):
        # 5 samples, window 2 -> 2 windows over the most recent 4 (drops the first).
        np.testing.assert_allclose(udstats.merge_stat(np.array([1.0, 2, 3, 4, 5]), 2), [2.5, 4.5])

    def test_variance_from_between_window_spread(self):
        mean, var = udstats.merge_stat(np.array([1.0, 2, 3, 4]), np.zeros(4), 2)
        np.testing.assert_allclose(mean, [1.5, 3.5])
        np.testing.assert_allclose(var, [0.25, 0.25])

    def test_covariance_two_variables(self):
        xm, ym, cov = udstats.merge_stat(
            np.array([1.0, 2, 3, 4]), np.array([2.0, 4, 6, 8]), np.zeros(4), 2
        )
        np.testing.assert_allclose(xm, [1.5, 3.5])
        np.testing.assert_allclose(ym, [3.0, 7.0])
        np.testing.assert_allclose(cov, [0.5, 0.5])

    def test_rejects_bad_args(self):
        with self.assertRaises(ValueError):
            udstats.merge_stat(np.array([1.0, 2, 3]), 0)  # n <= 0
        with self.assertRaises(ValueError):
            udstats.merge_stat(np.array([1.0]), 2)  # too few samples


class TestTimeAverage(unittest.TestCase):
    def test_mean_and_variance(self):
        mean, var = udstats.time_average(np.array([1.0, 2, 3, 4]))
        np.testing.assert_allclose(mean, [2.5])
        np.testing.assert_allclose(var, [1.25])

    def test_cross_covariance(self):
        xm, ym, cov = udstats.time_average(np.array([1.0, 2, 3, 4]), np.array([2.0, 4, 6, 8]))
        np.testing.assert_allclose((xm, ym, cov), ([2.5], [5.0], [2.5]))


class TestCoarsegrainField(unittest.TestCase):
    def test_shape_and_mass_conservation(self):
        var = np.arange(16.0).reshape(4, 4, 1)
        out = udstats.coarsegrain_field(var, [2.0], np.arange(4.0), np.arange(4.0))
        self.assertEqual(out.shape, (4, 4, 1, 1))
        # A normalised box filter conserves the field sum.
        self.assertAlmostEqual(float(out[..., 0].sum()), float(var.sum()), places=6)

    def test_rejects_non_3d(self):
        with self.assertRaises(ValueError):
            udstats.coarsegrain_field(np.zeros((4, 4)), [2.0], np.arange(4.0), np.arange(4.0))


class TestUDBaseWrappersDelegate(unittest.TestCase):
    def test_wrappers_match_module_functions(self):
        x = np.array([1.0, 2, 3, 4, 5, 6])
        np.testing.assert_allclose(UDBase.merge_stat(x, 3), udstats.merge_stat(x, 3))
        np.testing.assert_allclose(
            UDBase.time_average(x)[0], udstats.time_average(x)[0]
        )


if __name__ == "__main__":
    unittest.main()
