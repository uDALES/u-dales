"""Tests for the archived atmospheric shortwave forcing series."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from netCDF4 import Dataset
import numpy as np

from _common import PYTHON_DIR  # noqa: F401
from udprep.shortwave_forcing import write_shortwave_forcing


class TestShortwaveForcingWriter(unittest.TestCase):
    def test_preserves_source_values_and_metadata(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "shortwave_forcing.300.nc"
            times = np.array([0.0, 600.0])
            dni = np.array([0.0, 723.123456789])
            dsky = np.array([14.25, 92.987654321])
            zenith = np.array([95.0, 38.5])
            azimuth = np.array([-30.0, 172.0])
            ghi = np.array([14.25, 661.875])
            write_shortwave_forcing(
                path, times, dni, dsky, zenith, azimuth,
                source="HARMONIE ssrd; Erbs GHI split",
                start_time="2023-08-21T00:00:00", ghi=ghi,
            )

            with Dataset(path) as dataset:
                for name, expected in (
                    ("time", times), ("dni", dni), ("dsky", dsky),
                    ("solar_zenith", zenith), ("solar_azimuth_local", azimuth),
                    ("ghi", ghi),
                ):
                    np.testing.assert_array_equal(dataset.variables[name][:], expected)
                    self.assertEqual(dataset.variables[name].dtype, np.dtype("float64"))
                self.assertEqual(dataset.variables["dni"].units, "W m-2")
                self.assertEqual(dataset.variables["dsky"].units, "W m-2")
                self.assertEqual(dataset.variables["time"].units, "s")
                self.assertEqual(dataset.simulation_start, "2023-08-21T00:00:00")
                self.assertIn("Erbs", dataset.source)

            with self.assertRaises(FileExistsError):
                write_shortwave_forcing(
                    path, times, dni + 1, dsky, zenith, azimuth,
                    source="other", start_time="2023-08-21T00:00:00",
                )
            with Dataset(path) as dataset:
                np.testing.assert_array_equal(dataset.variables["dni"][:], dni)
            self.assertEqual(list(Path(tmp).glob(".*.tmp")), [])

    def test_without_ghi_and_invalid_inputs(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "shortwave_forcing.001.nc"
            values = np.array([1.0, 2.0])
            write_shortwave_forcing(
                path, np.array([0.0, 10.0]), values, values, values, values,
                source="uDALES isolar=1", start_time="2020-06-21T12:00:00",
            )
            with Dataset(path) as dataset:
                self.assertNotIn("ghi", dataset.variables)
            for times, dni in (
                (np.array([0.0, 0.0]), values),
                (np.array([0.0, 10.0]), np.array([1.0])),
                (np.array([0.0, 10.0]), np.array([np.nan, 2.0])),
            ):
                with self.subTest(times=times, dni=dni), self.assertRaises(ValueError):
                    write_shortwave_forcing(
                        path, times, dni, values, values, values,
                        source="test", start_time="2020-06-21T12:00:00",
                        overwrite=True,
                    )
            with Dataset(path) as dataset:
                np.testing.assert_array_equal(dataset.variables["dni"][:], values)
