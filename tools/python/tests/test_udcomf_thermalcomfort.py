"""Model-neutral comfort schema, radiation physics, and index tests."""

from importlib.util import find_spec
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from netCDF4 import Dataset
import numpy as np

from udcomf.thermalcomfort import ComfortParameters, calculate_indices
from udcomf.thermalcomfort.physics import (
    EPSILON, SIGMA, direct_body_factor, mean_radiant_temperature,
    saturation_pressure_pa, sensor_temperatures, wet_bulb_globe_temperature,
)


def _exchange(path: Path, *, model: str = "PALM", local_wind: float = 1.,
              direct_sun: bool = True, rh: float = 50.) -> None:
    ta = 298.15
    pressure = 101325.
    e = rh / 100. * float(saturation_pressure_pa(np.array([ta]))[0])
    qv = EPSILON * e / (pressure - (1 - EPSILON) * e)
    mask = np.array([[1, 1], [1, 0]], dtype=np.int8)
    fill = np.float32(-9.96921e36)
    times = np.arange(1, 25, dtype=float) * 3600.
    with Dataset(path, "w") as ds:
        for name, size in (("x", 2), ("y", 2), ("time", 24), ("bounds", 2)):
            ds.createDimension(name, size)
        for name, values in (("x", [5., 15.]), ("y", [5., 15.])):
            ds.createVariable(name, "f8", (name,))[:] = values
            ds[name].units = "m"
            ds.createVariable(f"{name}_bounds", "f8", (name, "bounds"))[:] = [[0., 10.], [10., 20.]]
            ds[f"{name}_bounds"].units = "m"
        for name, values, units in (("longitude", 2.33, "degrees_east"),
                                    ("latitude", 48.86, "degrees_north"),
                                    ("z_ground", 77., "m")):
            ds.createVariable(name, "f8", ("x", "y"))[:] = values
            ds[name].units = units
        ds["z_ground"].vertical_datum = "IGN69"
        ds.createVariable("receptor_height", "f8").assignValue(1.1)
        ds["receptor_height"].units = "m above local terrain"
        ds.createVariable("pedestrian_mask", "i1", ("x", "y"))[:] = mask
        ds.createVariable("crs", "i4").assignValue(0)
        ds["crs"].grid_mapping_name = "transverse_mercator"
        ds["crs"].crs_wkt = "TEST WKT"
        ds.createVariable("time", "f8", ("time",))[:] = times
        ds["time"].units = "seconds since 2023-08-21 00:00:00 UTC"
        ds["time"].calendar = "gregorian"
        ds["time"].bounds = "time_bounds"
        ds.createVariable("time_bounds", "f8", ("time", "bounds"))[:] = np.column_stack(
            (times - 900, times)
        )
        ds["time_bounds"].units = ds["time"].units
        for name, value in (("solar_zenith", np.pi / 3), ("solar_azimuth", np.pi)):
            ds.createVariable(name, "f8", ("time",))[:] = value
            ds[name].units = "rad"
        for name, value in (
            ("ta", ta), ("pabs", pressure), ("qv", qv), ("ws_local", local_wind),
            ("ws_10", 1.), ("sw_direct_normal", 0.),
        ):
            units = ("K" if name == "ta" else "Pa" if name == "pabs" else
                     "kg kg-1" if name == "qv" else "m s-1" if name.startswith("ws_") else "W m-2")
            data = np.full((2, 2, 24), value, dtype=np.float32)
            if name == "sw_direct_normal" and direct_sun:
                data[0, 1, :] = 800.
            data[mask == 0, :] = fill
            var = ds.createVariable(name, "f4", ("x", "y", "time"), fill_value=fill)
            var[:] = data
            var.units = units
            var.cell_methods = "time: mean"
            var.grid_mapping = "crs"
            var.coordinates = "longitude latitude"
        for side in ("upface", "downface", "northface", "southface", "eastface", "westface"):
            for prefix, value in (("sw_nondirect_", 0.), ("lw_", SIGMA * ta**4)):
                name = f"{prefix}{side}"
                data = np.full((2, 2, 24), value, dtype=np.float32)
                data[mask == 0, :] = fill
                var = ds.createVariable(name, "f4", ("x", "y", "time"), fill_value=fill)
                var[:] = data
                var.units = "W m-2"
                var.cell_methods = "time: mean"
                var.grid_mapping = "crs"
                var.coordinates = "longitude latitude"
        ds.Conventions = "CF-1.10"
        for name, value in (
            ("model_name", model), ("model_version", "test"),
            ("institution", "test"), ("simulation_identifier", "fixture"),
            ("horizontal_grid_spacing_m", "x=10,y=10"),
            ("vertical_grid_spacing_m", "1,2"), ("terrain_convention", "flat"),
            ("building_representation", "none"),
            ("model_axis_rotation_degrees_from_true_north", 0.),
            ("vertical_interpolation_method", "linear"),
            ("spinup_start_utc", "2023-08-20T18:00:00Z"),
            ("spinup_end_utc", "2023-08-21T00:00:00Z"),
        ):
            ds.setncattr(name, value)


class TestPhysics(unittest.TestCase):
    def test_direct_factor_depends_on_zenith(self):
        self.assertAlmostEqual(direct_body_factor(0.), 0.06)
        self.assertGreater(direct_body_factor(np.pi / 3), direct_body_factor(0.))
        self.assertEqual(direct_body_factor(np.pi / 2), 0.)

    def test_uniform_enclosure_and_saturated_wick(self):
        ta = np.array([298.15])
        pressure = np.array([101325.])
        e = saturation_pressure_pa(ta)
        qv = EPSILON * e / (pressure - (1 - EPSILON) * e)
        fields = {"sw_direct_normal": np.zeros(1)}
        for side in ("upface", "downface", "northface", "southface", "eastface", "westface"):
            fields[f"sw_nondirect_{side}"] = np.zeros(1)
            fields[f"lw_{side}"] = SIGMA * ta**4
        np.testing.assert_allclose(mean_radiant_temperature(fields, np.pi / 3), ta, atol=1e-10)
        globe, wick, converged = sensor_temperatures(
            ta, pressure, qv, np.ones(1), fields, np.pi / 3
        )
        self.assertTrue(converged[0])
        self.assertAlmostEqual(globe[0], ta[0], places=4)
        self.assertAlmostEqual(wick[0], ta[0], places=3)

    def test_six_direction_mrt_fixed_reference(self):
        """Check the documented six-direction equation against a fixed result."""
        fields = {"sw_direct_normal": np.array([600.0])}
        for side, value in {
            "upface": 120.0, "downface": 30.0, "northface": 50.0,
            "southface": 80.0, "eastface": 70.0, "westface": 40.0,
        }.items():
            fields[f"sw_nondirect_{side}"] = np.array([value])
        for side, value in {
            "upface": 350.0, "downface": 450.0, "northface": 380.0,
            "southface": 410.0, "eastface": 390.0, "westface": 400.0,
        }.items():
            fields[f"lw_{side}"] = np.array([value])

        mrt = mean_radiant_temperature(fields, np.deg2rad(60.0))

        # Independent evaluation of the six-plane standing-person equation:
        # weighted SW=225.2922678358 and weighted LW=395.6 W m-2.
        self.assertAlmostEqual(float(mrt[0]), 314.985841050193, places=9)

    def test_iso_wbgt_reference_values(self):
        """Reproduce the ISO 7243 examples documented by pythermalcomfort."""
        ta = np.array([293.15, 293.15])
        globe = np.array([305.15, 305.15])
        wick = np.array([298.15, 298.15])
        direct = np.array([0.0, 1.0])

        result = wet_bulb_globe_temperature(ta, globe, wick, direct)

        np.testing.assert_allclose(result, [27.1, 25.9], rtol=0, atol=1e-12)

    def test_directional_wbgt_sensor_balances_match_independent_roots(self):
        """Compare fixed-point sensor temperatures with independently solved roots."""
        ta = np.array([303.15])
        pressure = np.array([100000.0])
        e = 0.5 * saturation_pressure_pa(ta)
        qv = EPSILON * e / (pressure - (1.0 - EPSILON) * e)
        fields = {"sw_direct_normal": np.array([700.0])}
        for side, value in {
            "upface": 120.0, "downface": 60.0, "northface": 80.0,
            "southface": 100.0, "eastface": 70.0, "westface": 90.0,
        }.items():
            fields[f"sw_nondirect_{side}"] = np.array([value])
        for side, value in {
            "upface": 360.0, "downface": 480.0, "northface": 440.0,
            "southface": 460.0, "eastface": 420.0, "westface": 430.0,
        }.items():
            fields[f"lw_{side}"] = np.array([value])

        globe, wick, converged = sensor_temperatures(
            ta, pressure, qv, np.array([1.2]), fields, np.deg2rad(60.0)
        )
        wbgt = wet_bulb_globe_temperature(
            ta, globe, wick, fields["sw_direct_normal"]
        )

        # These roots were evaluated independently with a bracketed scalar
        # solver applied to the Liljegren globe and wick balance equations.
        self.assertTrue(converged[0])
        self.assertAlmostEqual(float(globe[0]), 311.701170960989, delta=0.02)
        self.assertAlmostEqual(float(wick[0]), 296.051700718050, delta=0.02)
        self.assertAlmostEqual(float(wbgt[0]), 26.741424694833, delta=0.02)


@unittest.skipUnless(find_spec("pythermalcomfort"), "install tools/python[comfort]")
class TestComfortFiles(unittest.TestCase):
    def test_all_indices_native_grid_and_sunshade(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / "exchange.nc"
            _exchange(source)
            output = calculate_indices(source)
            with Dataset(output) as ds:
                ds.set_auto_mask(False)
                self.assertEqual(ds.model_name, "PALM")
                self.assertEqual(ds["mrt"].dimensions, ("x", "y", "time"))
                self.assertEqual(float(ds["mrt"][0, 0, 0]), 25.)
                self.assertGreater(float(ds["mrt"][0, 1, 0]), 25.)
                self.assertGreater(float(ds["wbgt"][0, 1, 0]), float(ds["wbgt"][0, 0, 0]))
                self.assertTrue(np.isfinite(float(ds["pet"][0, 0, 0])))
                self.assertTrue(np.isfinite(float(ds["utci"][0, 0, 0])))
                self.assertEqual(float(ds["mrt"][1, 1, 0]), float(ds["mrt"]._FillValue))
                self.assertEqual(int(ds["validity_flags"][1, 1, 0]), -1)
            self.assertEqual(source.name, "exchange.nc")

    def test_published_pet_and_utci_reference_values(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / "exchange.nc"
            _exchange(source, local_wind=0.1, direct_sun=False)
            output = calculate_indices(
                source, parameters=ComfortParameters(pet_met=1.2, pet_clo=0.5)
            )
            with Dataset(output) as ds:
                self.assertAlmostEqual(float(ds["pet"][0, 0, 0]), 24.67, delta=0.05)
                self.assertAlmostEqual(float(ds["utci"][0, 0, 0]), 24.6, delta=0.1)
                self.assertEqual(int(ds["validity_flags"][0, 0, 0]) & 4, 4)

    def test_same_physics_for_different_model_names(self):
        with TemporaryDirectory() as tmp:
            first = Path(tmp) / "palm.nc"
            second = Path(tmp) / "urbclim.nc"
            _exchange(first, model="PALM")
            _exchange(second, model="UrbClim")
            a = calculate_indices(first)
            b = calculate_indices(second)
            with Dataset(a) as left, Dataset(b) as right:
                for field in ("mrt", "pet", "utci", "wbgt"):
                    np.testing.assert_array_equal(left[field][:], right[field][:])

    def test_rejects_bad_schema_and_missing_inside_mask(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / "exchange.nc"
            _exchange(source)
            with Dataset(source, "r+") as ds:
                ds["qv"][0, 0, 0] = ds["qv"]._FillValue
            with self.assertRaisesRegex(ValueError, "qv has missing"):
                calculate_indices(source)
            self.assertFalse(list(Path(tmp).glob("thermal_comfort_indices*")))

    def test_invalid_humidity_is_flagged_without_clipping_or_losing_mrt(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / "exchange.nc"
            _exchange(source)
            with Dataset(source, "r+") as ds:
                ds["qv"][0, 0, :] = 0.04
            output = calculate_indices(source)
            with Dataset(output) as ds:
                ds.set_auto_mask(False)
                self.assertTrue(np.isfinite(ds["mrt"][0, 0, 0]))
                self.assertEqual(float(ds["pet"][0, 0, 0]), float(ds["pet"]._FillValue))
                self.assertEqual(float(ds["utci"][0, 0, 0]), float(ds["utci"]._FillValue))
                self.assertEqual(int(ds["validity_flags"][0, 0, 0]) & 1, 1)

    def test_utci_low_wind_is_missing_and_flagged(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / "exchange.nc"
            _exchange(source)
            with Dataset(source, "r+") as ds:
                ds["ws_10"][0, 0, :] = 0.1
            output = calculate_indices(source)
            with Dataset(output) as ds:
                ds.set_auto_mask(False)
                self.assertEqual(float(ds["utci"][0, 0, 0]), float(ds["utci"]._FillValue))
                self.assertEqual(int(ds["validity_flags"][0, 0, 0]) & 2, 2)
                self.assertTrue(np.isfinite(ds["pet"][0, 0, 0]))

    def test_rejects_direct_beam_below_horizon(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / "exchange.nc"
            _exchange(source)
            with Dataset(source, "r+") as ds:
                ds["solar_zenith"][0] = np.pi / 2
                ds["sw_direct_normal"][0, 1, 0] = 100.
            with self.assertRaisesRegex(ValueError, "Direct shortwave"):
                calculate_indices(source)
