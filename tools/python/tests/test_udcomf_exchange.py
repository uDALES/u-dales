"""Strict per-height uDALES thermal-comfort exchange tests."""

from datetime import date, datetime, timezone
from pathlib import Path
from shutil import copyfile
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from netCDF4 import Dataset
import numpy as np

from udcomf.udcomf_export import (
    ExchangeMetadata, _ATMOSPHERE, _SHORTWAVE, _LONGWAVE, write_exchange_heights,
)
from udprep.solar import solar_position_python


def _fixture(root: Path):
    case = root / "001"
    case.mkdir()
    (case / "namoptions.001").write_text(
        "&DOMAIN\n year=2023\n month=8\n day=21\n hour=0\n"
        " longitude=2.33\n latitude=48.86\n elevation=35.\n timezone=0.\n/\n"
        "&OUTPUT\n receptor_height=1.1\n/\n", encoding="ascii",
    )
    sim = SimpleNamespace(
        path=case, expnr="001", xt=np.array([5., 15.]), yt=np.array([5., 15.]),
        dx=10., dy=10., dzt=np.array([1., 2.]), xazimuth=89.5,
    )
    meta = ExchangeMetadata(
        x=np.array([1005., 1015.]), y=np.array([2005., 2015.]),
        x_bounds=np.array([[1000., 1010.], [1010., 1020.]]),
        y_bounds=np.array([[2000., 2010.], [2010., 2020.]]),
        longitude=np.full((2, 2), 2.33), latitude=np.full((2, 2), 48.86),
        z_ground=np.full((2, 2), 77.),
        crs_attributes={"grid_mapping_name": "transverse_mercator", "crs_wkt": "TEST WKT"},
        vertical_datum="IGN69", model_version="test-revision", institution="test-institute",
        building_representation="IBM facets", terrain_convention="flat local ground",
        spinup_start_utc="2023-08-21T00:00:00Z", spinup_end_utc="2023-08-21T00:00:00Z",
    )
    times = np.arange(1, 25) * 3600.
    for kind, fields in (("atmosphere", _ATMOSPHERE),
                         ("shortwave", _SHORTWAVE), ("longwave", _LONGWAVE)):
        path = case / f"pedestrian_{kind}.h1p1.001.nc"
        with Dataset(path, "w") as ds:
            for name, length in (("x", 2), ("y", 2), ("time", 24), ("bounds", 2)):
                ds.createDimension(name, length)
            ds.createVariable("x", "f8", ("x",))[:] = sim.xt
            ds.createVariable("y", "f8", ("y",))[:] = sim.yt
            ds.createVariable("time", "f8", ("time",))[:] = times
            if kind == "atmosphere":
                ds.createVariable("source_time", "f8", ("time",))[:] = times
            ds.createVariable("time_bounds", "f8", ("time", "bounds"))[:] = np.column_stack((times - 900, times))
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
    return sim, meta


class TestExchange(unittest.TestCase):
    def test_writes_strict_schema_and_mask(self):
        with TemporaryDirectory() as tmp:
            sim, meta = _fixture(Path(tmp))
            path = write_exchange_heights(sim, date(2023, 8, 21), meta)[1.1]
            with Dataset(path) as ds:
                ds.set_auto_mask(False)
                self.assertEqual(ds["ta"].dimensions, ("x", "y", "time"))
                self.assertEqual(ds["ta"].cell_methods, "time: mean")
                self.assertEqual(ds["time"].calendar, "gregorian")
                self.assertEqual(ds["time"].bounds, "time_bounds")
                self.assertEqual(float(ds["time"][0]), 3600.)
                self.assertEqual(float(ds["time_bounds"][-1, 0]), 85500.)
                self.assertEqual(float(ds["x"][0]), 1005.)
                self.assertEqual(float(ds["qv"][0, 0, 0]), np.float32(.01))
                self.assertEqual(float(ds["ta"][1, 1, 0]), float(ds["ta"]._FillValue))
                self.assertEqual(int(ds["pedestrian_mask"][1, 1]), 0)
                self.assertEqual(len(ds.variables), 32)
                midpoint = datetime(2023, 8, 21, 0, 52, 30, tzinfo=timezone.utc)
                expected_azimuth = np.deg2rad(solar_position_python(
                    midpoint, 2.33, 48.86, 0., 35.
                )["azimuth"]) % (2 * np.pi)
                self.assertAlmostEqual(float(ds["solar_azimuth"][0]), expected_azimuth)

    def test_writes_separate_files_per_saved_height(self):
        with TemporaryDirectory() as tmp:
            sim, meta = _fixture(Path(tmp))
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
            outputs = write_exchange_heights(sim, date(2023, 8, 21), meta)
            self.assertEqual(set(outputs), {1.1, 1.5})
            with Dataset(outputs[1.5]) as ds:
                self.assertEqual(float(ds["receptor_height"][...]), 1.5)

    def test_rejects_incomplete_radiation_without_output(self):
        with TemporaryDirectory() as tmp:
            sim, meta = _fixture(Path(tmp))
            with Dataset(sim.path / "pedestrian_longwave.h1p1.001.nc", "r+") as ds:
                ds["time_complete"][4] = 0
            with self.assertRaisesRegex(ValueError, "incomplete radiation"):
                write_exchange_heights(sim, date(2023, 8, 21), meta)
            self.assertFalse(list(sim.path.glob("thermal_comfort_inputs*")))

    def test_rejects_partial_mask_and_negative_humidity(self):
        with TemporaryDirectory() as tmp:
            sim, meta = _fixture(Path(tmp))
            shortwave = sim.path / "pedestrian_shortwave.h1p1.001.nc"
            with Dataset(shortwave, "r+") as ds:
                ds["processed_mask"][0, 0] = 0
            with self.assertRaisesRegex(ValueError, "does not cover"):
                write_exchange_heights(sim, date(2023, 8, 21), meta)
            with Dataset(shortwave, "r+") as ds:
                ds["processed_mask"][0, 0] = 1
            with Dataset(sim.path / "pedestrian_atmosphere.h1p1.001.nc", "r+") as ds:
                ds["qv"][0, 0, 3] = -1e-9
            with self.assertRaisesRegex(ValueError, "nonnegative"):
                write_exchange_heights(sim, date(2023, 8, 21), meta)
            self.assertFalse(list(sim.path.glob("thermal_comfort_inputs*")))

    def test_rejects_statistics_record_at_wrong_hour(self):
        with TemporaryDirectory() as tmp:
            sim, meta = _fixture(Path(tmp))
            with Dataset(sim.path / "pedestrian_atmosphere.h1p1.001.nc", "r+") as ds:
                ds["source_time"][2] = 10802.
            with self.assertRaisesRegex(ValueError, "statistics records miss"):
                write_exchange_heights(sim, date(2023, 8, 21), meta)

    def test_rejects_mismatched_georeference(self):
        with TemporaryDirectory() as tmp:
            sim, meta = _fixture(Path(tmp))
            meta.x_bounds[0, 1] = 1009.
            with self.assertRaisesRegex(ValueError, "x coordinates/bounds"):
                write_exchange_heights(sim, date(2023, 8, 21), meta)

    def test_rejects_variable_ground_for_flat_model(self):
        with TemporaryDirectory() as tmp:
            sim, meta = _fixture(Path(tmp))
            meta.z_ground[0, 0] = 76.
            with self.assertRaisesRegex(ValueError, "flat model ground"):
                write_exchange_heights(sim, date(2023, 8, 21), meta)
