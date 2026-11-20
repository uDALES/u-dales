"""Per-height uDALES thermal-comfort exchange tests."""

from datetime import datetime, timedelta, timezone
from pathlib import Path
from shutil import copyfile
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from netCDF4 import Dataset
import numpy as np

from udcomf.udcomf_export import (
    _ATMOSPHERE, _SHORTWAVE, _LONGWAVE, write_exchange_heights,
)
from udprep.solar import solar_position, solar_position_python


def _fixture(root: Path):
    case = root / "001"
    case.mkdir()
    (case / "namoptions.001").write_text(
        "&DOMAIN\n year=2023\n month=8\n day=21\n hour=0\n"
        " longitude=2.33\n latitude=48.86\n elevation=35.\n timezone=0.\n/\n"
        "&OUTPUT\n receptor_height=1.1\n tstatsdump=600.\n/\n", encoding="ascii",
    )
    sim = SimpleNamespace(
        path=case, expnr="001", xt=np.array([5., 15.]), yt=np.array([5., 15.]),
        dx=10., dy=10., dzt=np.array([1., 2.]), xazimuth=89.5,
    )
    times = np.array([1800., 5400., 9000.])
    for kind, fields in (("atmosphere", _ATMOSPHERE),
                         ("shortwave", _SHORTWAVE), ("longwave", _LONGWAVE)):
        path = case / f"pedestrian_{kind}.h1p1.001.nc"
        with Dataset(path, "w") as ds:
            for name, length in (("x", 2), ("y", 2), ("time", len(times)), ("bounds", 2)):
                ds.createDimension(name, length)
            ds.createVariable("x", "f8", ("x",))[:] = sim.xt
            ds.createVariable("y", "f8", ("y",))[:] = sim.yt
            ds.createVariable("time", "f8", ("time",))[:] = times
            ds.createVariable("time_bounds", "f8", ("time", "bounds"))[:] = np.column_stack((times - 600, times))
            ds.createVariable("receptor_height", "f8").assignValue(1.1)
            ds.createVariable("valid_mask", "i1", ("x", "y"))[:] = [[1, 1], [1, 0]]
            if kind != "atmosphere":
                ds.createVariable("time_complete", "i1", ("time",))[:] = 1
                ds.createVariable("processed_mask", "i1", ("x", "y"))[:] = [[1, 1], [1, 0]]
            for name in fields:
                units = ("K" if name == "ta" else "Pa" if name == "pabs" else
                         "kg kg-1" if name == "qv" else "m s-1" if name.startswith("ws_") else "W m-2")
                var = ds.createVariable(name, "f4", ("x", "y", "time"))
                var.units = units
                var[:] = (300. if name == "ta" else 100000. if name == "pabs" else
                          0.01 if name == "qv" else 5.)
    return sim


class TestExchange(unittest.TestCase):
    def test_writes_strict_schema_and_mask(self):
        with TemporaryDirectory() as tmp:
            sim = _fixture(Path(tmp))
            path = write_exchange_heights(sim)[1.1]
            with Dataset(path) as ds:
                ds.set_auto_mask(False)
                self.assertEqual(ds["ta"].dimensions, ("x", "y", "time"))
                self.assertEqual(ds["ta"].cell_methods, "time: mean")
                self.assertEqual(ds["time"].calendar, "gregorian")
                self.assertEqual(ds["time"].bounds, "time_bounds")
                self.assertEqual(float(ds["time"][0]), 1800.)
                self.assertEqual(ds["time"].units, "seconds since 2023-08-21 00:00:00 UTC")
                self.assertEqual(float(ds["time_bounds"][-1, 0]), 8400.)
                self.assertEqual(float(ds.averaging_window_s), 600.)
                self.assertEqual(ds.simulation_start_utc, "2023-08-21T00:00:00Z")
                self.assertEqual(float(ds["x"][0]), 5.)
                np.testing.assert_allclose(ds["x_bounds"][:], [[0., 10.], [10., 20.]])
                self.assertEqual(float(ds.site_longitude_degrees_east), 2.33)
                self.assertEqual(float(ds.model_axis_rotation_degrees_from_true_north), 89.5)
                for absent in ("longitude", "latitude", "z_ground", "crs"):
                    self.assertNotIn(absent, ds.variables)
                self.assertEqual(float(ds["qv"][0, 0, 0]), np.float32(.01))
                self.assertEqual(float(ds["ta"][1, 1, 0]), float(ds["ta"]._FillValue))
                self.assertEqual(int(ds["pedestrian_mask"][1, 1]), 0)
                self.assertEqual(len(ds.variables), 28)
                midpoint = datetime(2023, 8, 21, tzinfo=timezone.utc) + timedelta(seconds=1500)
                expected_azimuth = np.deg2rad(solar_position_python(
                    midpoint, 2.33, 48.86, 0., 35.
                )["azimuth"]) % (2 * np.pi)
                self.assertAlmostEqual(float(ds["solar_azimuth"][0]), expected_azimuth)

    def test_isolarpos_selects_the_preprocessing_solar_backend(self):
        with TemporaryDirectory() as tmp:
            sim = _fixture(Path(tmp))
            options = sim.path / "namoptions.001"
            options.write_text(options.read_text(encoding="ascii").replace(
                " tstatsdump=600.", " tstatsdump=600.\n isolarpos=2"
            ), encoding="ascii")
            path = write_exchange_heights(sim)[1.1]
            with Dataset(path) as ds:
                midpoint = datetime(2023, 8, 21, tzinfo=timezone.utc) + timedelta(seconds=1500)
                _, expected = solar_position(midpoint, 2.33, 48.86, 0., 35., backend="pvlib")
                self.assertAlmostEqual(float(ds["solar_azimuth"][0]), np.deg2rad(expected) % (2 * np.pi))
                self.assertIn("pvlib", ds.solar_geometry)
                self.assertIn("isolarpos=2", ds.solar_geometry)
            options.write_text(options.read_text(encoding="ascii").replace("isolarpos=2", "isolarpos=7"),
                               encoding="ascii")
            with self.assertRaisesRegex(ValueError, "isolarpos"):
                write_exchange_heights(sim, overwrite=True)

    def test_writes_separate_files_per_saved_height(self):
        with TemporaryDirectory() as tmp:
            sim = _fixture(Path(tmp))
            options = sim.path / "namoptions.001"
            options.write_text(options.read_text(encoding="ascii").replace(
                "receptor_height=1.1", "receptor_height=1.1\n nreceptor_heights=2\n receptor_heights=1.1,1.5"
            ), encoding="ascii")
            for kind in ("atmosphere", "shortwave", "longwave"):
                original = sim.path / f"pedestrian_{kind}.h1p1.001.nc"
                second = sim.path / f"pedestrian_{kind}.h1p5.001.nc"
                copyfile(original, second)
                with Dataset(second, "r+") as ds:
                    ds["receptor_height"].assignValue(1.5)
            outputs = write_exchange_heights(sim)
            self.assertEqual(set(outputs), {1.1, 1.5})
            with Dataset(outputs[1.5]) as ds:
                self.assertEqual(float(ds["receptor_height"][...]), 1.5)

    def test_rejects_incomplete_radiation_without_output(self):
        with TemporaryDirectory() as tmp:
            sim = _fixture(Path(tmp))
            with Dataset(sim.path / "pedestrian_longwave.h1p1.001.nc", "r+") as ds:
                ds["time_complete"][1] = 0
            with self.assertRaisesRegex(ValueError, "incomplete radiation"):
                write_exchange_heights(sim)
            self.assertFalse(list(sim.path.glob("thermal_comfort_inputs*")))

    def test_rejects_partial_mask_and_negative_humidity(self):
        with TemporaryDirectory() as tmp:
            sim = _fixture(Path(tmp))
            shortwave = sim.path / "pedestrian_shortwave.h1p1.001.nc"
            with Dataset(shortwave, "r+") as ds:
                ds["processed_mask"][0, 0] = 0
            with self.assertRaisesRegex(ValueError, "does not cover"):
                write_exchange_heights(sim)
            with Dataset(shortwave, "r+") as ds:
                ds["processed_mask"][0, 0] = 1
            with Dataset(sim.path / "pedestrian_atmosphere.h1p1.001.nc", "r+") as ds:
                ds["qv"][0, 0, 2] = -1e-9
            with self.assertRaisesRegex(ValueError, "nonnegative"):
                write_exchange_heights(sim)
            self.assertFalse(list(sim.path.glob("thermal_comfort_inputs*")))

    def test_rejects_radiation_records_that_differ_from_atmosphere(self):
        with TemporaryDirectory() as tmp:
            sim = _fixture(Path(tmp))
            with Dataset(sim.path / "pedestrian_shortwave.h1p1.001.nc", "r+") as ds:
                ds["time"][2] = 9002.
            with self.assertRaisesRegex(ValueError, "different time grid"):
                write_exchange_heights(sim)

    def test_rejects_bounds_inconsistent_with_tstatsdump(self):
        with TemporaryDirectory() as tmp:
            sim = _fixture(Path(tmp))
            with Dataset(sim.path / "pedestrian_atmosphere.h1p1.001.nc", "r+") as ds:
                ds["time_bounds"][1, 0] = 4500.
            with self.assertRaisesRegex(ValueError, "averaging window"):
                write_exchange_heights(sim)
