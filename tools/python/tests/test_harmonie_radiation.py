"""Tests for HARMONIE radiation conversion helpers."""

from datetime import datetime
from pathlib import Path
from tempfile import TemporaryDirectory
import types
import unittest
from unittest import mock

import numpy as np

from _common import PYTHON_DIR  # noqa: F401

from udprep.harmonie_radiation import (
    ShortwaveAtmosphere,
    accumulated_flux_series,
    format_forecast_offset,
    generate_timedepsw_from_harmonie,
    interpolate_flux_to_times,
    make_longwave_times,
    make_model_times,
    split_global_horizontal_erbs,
    write_netsw,
    write_timedeplw,
)


class FakeAccumulationReader:
    field = "ssrd"

    def __init__(self, values):
        self.values = values

    def mean_accumulation(self, offset_seconds):
        return self.values[offset_seconds]


class TestHarmonieRadiationTimes(unittest.TestCase):
    def test_format_forecast_offset(self):
        self.assertEqual(format_forecast_offset(18 * 3600), "0018h00m00s")
        self.assertEqual(format_forecast_offset(18 * 3600 + 15 * 60), "0018h15m00s")
        self.assertEqual(format_forecast_offset(48 * 3600), "0048h00m00s")

    def test_make_model_times_uses_nt_when_available(self):
        np.testing.assert_allclose(make_model_times(108060.0, 600.0, 181)[[0, -1]], [0.0, 108000.0])
        self.assertEqual(make_model_times(108060.0, 600.0, 181).size, 181)

    def test_make_longwave_times_preserves_hourly_case_300_style(self):
        times = make_longwave_times(108060.0, 31)
        self.assertEqual(times.size, 31)
        np.testing.assert_allclose(times[[0, -1]], [0.0, 108000.0])


class TestAccumulatedFluxSeries(unittest.TestCase):
    def test_differences_accumulated_energy_to_flux(self):
        reader = FakeAccumulationReader(
            {
                900: 9000.0,
                1800: 27000.0,
                2700: 54000.0,
                3600: 90000.0,
            }
        )
        times, fluxes = accumulated_flux_series(
            reader,
            start_offset_seconds=1800,
            end_offset_seconds=3600,
            difference_interval_seconds=900,
        )
        np.testing.assert_allclose(times, [0.0, 900.0, 1800.0])
        np.testing.assert_allclose(fluxes, [20.0, 30.0, 40.0])

    def test_cycle_start_uses_zero_previous_accumulation(self):
        reader = FakeAccumulationReader(
            {
                0: 0.0,
                900: 9000.0,
                1800: 27000.0,
            }
        )
        times, fluxes = accumulated_flux_series(
            reader,
            start_offset_seconds=0,
            end_offset_seconds=1800,
            difference_interval_seconds=900,
        )
        np.testing.assert_allclose(times, [0.0, 900.0, 1800.0])
        np.testing.assert_allclose(fluxes, [0.0, 10.0, 20.0])

    def test_clamps_tiny_negative_accumulation_noise(self):
        reader = FakeAccumulationReader({900: 1000.0, 1800: 999.5})
        _, fluxes = accumulated_flux_series(
            reader,
            start_offset_seconds=1800,
            end_offset_seconds=1800,
            difference_interval_seconds=900,
        )
        np.testing.assert_allclose(fluxes, [0.0])

    def test_interpolates_native_flux_to_target_times(self):
        source_times = np.array([0.0, 900.0, 1800.0])
        source_fluxes = np.array([0.0, 90.0, 180.0])
        target_times = np.array([0.0, 600.0, 1200.0, 1800.0])
        np.testing.assert_allclose(
            interpolate_flux_to_times(source_times, source_fluxes, target_times),
            [0.0, 60.0, 120.0, 180.0],
        )


class TestErbsSplit(unittest.TestCase):
    def test_zero_or_nighttime_has_no_shortwave(self):
        when = datetime(2023, 8, 21, 0)
        self.assertEqual(split_global_horizontal_erbs(0.0, 45.0, when), (0.0, 0.0))
        self.assertEqual(split_global_horizontal_erbs(100.0, 95.0, when), (0.0, 0.0))

    def test_split_conserves_global_horizontal_irradiance(self):
        when = datetime(2023, 8, 21, 12)
        ghi = 650.0
        zenith = 35.0
        dni, dsky = split_global_horizontal_erbs(ghi, zenith, when)
        reconstructed = dni * np.cos(np.radians(zenith)) + dsky
        self.assertGreater(dni, 0.0)
        self.assertGreater(dsky, 0.0)
        self.assertLessEqual(dsky, ghi)
        self.assertAlmostEqual(reconstructed, ghi, places=10)

    def test_very_low_sun_is_treated_as_diffuse(self):
        when = datetime(2023, 8, 21, 20)
        self.assertEqual(
            split_global_horizontal_erbs(25.0, 89.9, when, min_direct_cos_zenith=0.01),
            (0.0, 25.0),
        )


class TestLongwaveWriter(unittest.TestCase):
    def test_write_timedeplw_format(self):
        from tempfile import TemporaryDirectory
        from pathlib import Path

        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "timedeplw.inp.300"
            write_timedeplw(
                path,
                np.array([0.0, 3600.0]),
                np.array([300.0, 301.25]),
                overwrite=False,
            )
            self.assertEqual(
                path.read_text(encoding="ascii").splitlines(),
                [
                    "Time-varying sky longwave flux from HARMONIE strd",
                    "time             LWsky",
                    "     0.000000   300.000000",
                    "  3600.000000   301.250000",
                ],
            )


class TestShortwaveWriters(unittest.TestCase):
    def test_write_netsw_format(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "netsw.inp.300"
            write_netsw(path, np.array([1.25, 2.5]), overwrite=False)

            lines = path.read_text(encoding="ascii").splitlines()
            self.assertIn("initial column", lines[0])
            np.testing.assert_allclose(np.loadtxt(path, skiprows=1), [1.25, 2.5])

    def test_write_netsw_requires_1d_values(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "netsw.inp.300"
            with self.assertRaisesRegex(ValueError, "one-dimensional"):
                write_netsw(path, np.ones((2, 2)), overwrite=False)

    def test_generate_timedepsw_writes_matching_netsw(self):
        with TemporaryDirectory() as tmp:
            case_dir = Path(tmp) / "300"
            case_dir.mkdir()
            prep = types.SimpleNamespace(
                sim=types.SimpleNamespace(path=case_dir, expnr="300"),
                radiation=types.SimpleNamespace(),
            )
            atmosphere = ShortwaveAtmosphere(
                times=np.array([0.0, 10.0]),
                ghi=np.array([100.0, 110.0]),
                dni=np.array([80.0, 90.0]),
                dsky=np.array([20.0, 20.0]),
                zenith=np.array([30.0, 31.0]),
                azimuth_local=np.array([180.0, 181.0]),
            )
            sdir = np.array([[1.0, 2.0], [3.0, 4.0]])
            knet = np.array([[10.0, 11.0], [20.0, 21.0]])

            with (
                mock.patch("udprep.udprep.UDPrep", return_value=prep),
                mock.patch(
                    "udprep.harmonie_radiation.prepare_harmonie_ssrd_atmosphere",
                    return_value=(atmosphere, types.SimpleNamespace(mask_points=4)),
                ),
                mock.patch(
                    "udprep.harmonie_radiation.map_atmosphere_to_facets",
                    return_value=(sdir, knet, None),
                ),
            ):
                result = generate_timedepsw_from_harmonie(
                    case_dir=case_dir,
                    nwp_root=Path(tmp),
                    write_sdir_nc=False,
                    verbose=False,
                )

            self.assertEqual(result.netsw_path, case_dir / "netsw.inp.300")
            np.testing.assert_allclose(
                np.loadtxt(result.netsw_path, skiprows=1), knet[:, 0]
            )


if __name__ == "__main__":
    unittest.main()
