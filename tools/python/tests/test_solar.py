from __future__ import annotations

from datetime import datetime, timedelta
import unittest

import numpy as np

from _common import PYTHON_DIR
from udprep import solar

def _angle_diff_deg(a: float, b: float) -> float:
    diff = (a - b + 180.0) % 360.0 - 180.0
    return abs(diff)

class TestSolar(unittest.TestCase):
    def test_python_backend_tracks_pvlib(self):
        longitude = -0.13
        latitude = 51.5
        timezone = 0.0
        elevation = 0.0
        xazimuth = 90.0

        start_time = datetime(2011, 9, 30, 6, 0, 0)
        times = [start_time + n * timedelta(minutes=30) for n in range(25)]

        zen_diff = []
        az_diff = []
        i_diff = []
        dsky_diff = []

        for current_time in times:
            _, zen_py, az_py, i_py, dsky_py = solar.solar_state(
                current_time,
                longitude,
                latitude,
                timezone,
                elevation,
                xazimuth=xazimuth,
                backend="python",
            )
            _, zen_pv, az_pv, i_pv, dsky_pv = solar.solar_state(
                current_time,
                longitude,
                latitude,
                timezone,
                elevation,
                xazimuth=xazimuth,
                backend="pvlib",
            )
            zen_diff.append(_angle_diff_deg(float(zen_py), float(zen_pv)))
            az_diff.append(_angle_diff_deg(float(az_py), float(az_pv)))
            i_diff.append(abs(float(i_py) - float(i_pv)))
            dsky_diff.append(abs(float(dsky_py) - float(dsky_pv)))

        self.assertLess(max(zen_diff), 1.0)
        self.assertLess(max(az_diff), 2.0)
        self.assertLess(max(i_diff), 100.0)
        self.assertLess(max(dsky_diff), 100.0)

    def test_solar_position_entry_point_matches_backends_and_maps_isolarpos(self):
        when = datetime(2023, 8, 21, 12, 0, 0)
        zen_py, az_py = solar.solar_position(when, 2.33, 48.86, 0.0, 35.0, backend="python")
        spa = solar.solar_position_python(when, 2.33, 48.86, 0.0, 35.0)
        self.assertEqual((zen_py, az_py), (float(spa["zenith"]), float(spa["azimuth"])))
        zen_pv, az_pv = solar.solar_position(when, 2.33, 48.86, 0.0, 35.0, backend="pvlib")
        self.assertLess(_angle_diff_deg(zen_py, zen_pv), 0.1)
        self.assertLess(_angle_diff_deg(az_py, az_pv), 0.1)
        self.assertNotEqual((zen_py, az_py), (zen_pv, az_pv))
        self.assertEqual(solar.solar_position_backend(1), "python")
        self.assertEqual(solar.solar_position_backend(2), "pvlib")
        with self.assertRaisesRegex(ValueError, "isolarpos"):
            solar.solar_position_backend(3)
        with self.assertRaises(ValueError):
            solar.solar_position(when, 2.33, 48.86, 0.0, 35.0, backend="other")


if __name__ == "__main__":
    unittest.main()
