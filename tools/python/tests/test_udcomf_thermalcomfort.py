"""Exchange schema, radiation physics, and comfort index tests."""

from dataclasses import fields, replace
from importlib.util import find_spec
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from netCDF4 import Dataset
import numpy as np

from udcomf.thermalcomfort import (
    DEFAULT_PARAMETERS_PATH, INDEX_NAMES, ComfortParameters, calculate_indices,
    default_output_path, required_fields, selected_indices,
)
from udcomf.thermalcomfort.physics import (
    EPSILON, SIGMA, direct_body_factor, mean_radiant_temperature,
    saturation_pressure_pa, sensor_temperatures, wet_bulb_globe_temperature,
)

P = ComfortParameters.defaults()
EXCHANGE_NAME = "thermal_comfort_inputs.h1p1.1.nc"
INDICES_NAME = "thermal_comfort_indices.h1p1.1.nc"


def _mrt(fields_, zenith):
    return mean_radiant_temperature(
        fields_, zenith, alpha_sw=P.human_shortwave_absorptivity,
        emissivity=P.human_longwave_emissivity, factor_horizontal=P.human_horizontal_factor,
        factor_vertical=P.human_vertical_factor, direct_side_factor=P.human_direct_side_factor,
    )


def _sensors(ta, pressure, qv, wind, fields_, zenith):
    return sensor_temperatures(
        ta, pressure, qv, wind, fields_, zenith,
        globe_diameter=P.wbgt_globe_diameter_m, globe_emissivity=P.wbgt_globe_emissivity,
        wick_diameter=P.wbgt_wick_diameter_m, wick_length=P.wbgt_wick_length_m,
        wick_emissivity=P.wbgt_wick_emissivity,
        wick_shortwave_absorptivity=P.wbgt_wick_shortwave_absorptivity,
        wind_min=P.wbgt_wind_min_m_s, iterations=P.wbgt_max_iterations,
        tolerance_k=P.wbgt_tolerance_k,
    )


def _wbgt(ta, globe, wick, direct):
    return wet_bulb_globe_temperature(
        ta, globe, wick, direct, sun_weights=P.wbgt_sun_weights,
        shade_weights=P.wbgt_shade_weights,
    )


def _exchange(path: Path, *, local_wind: float = 1.,
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
        ds.createVariable("receptor_height", "f8").assignValue(1.1)
        ds["receptor_height"].units = "m above model ground"
        ds.createVariable("pedestrian_mask", "i1", ("x", "y"))[:] = mask
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
        for side in ("upface", "downface", "northface", "southface", "eastface", "westface"):
            for prefix, value in (("sw_nondirect_", 0.), ("lw_", SIGMA * ta**4)):
                name = f"{prefix}{side}"
                data = np.full((2, 2, 24), value, dtype=np.float32)
                data[mask == 0, :] = fill
                var = ds.createVariable(name, "f4", ("x", "y", "time"), fill_value=fill)
                var[:] = data
                var.units = "W m-2"
                var.cell_methods = "time: mean"
        ds.Conventions = "CF-1.10"
        for name, value in (
            ("simulation_identifier", "fixture"),
            ("horizontal_grid_spacing_m", "x=10,y=10"),
            ("vertical_grid_spacing_m", "1,2"),
            ("model_axis_rotation_degrees_from_true_north", 0.),
            ("vertical_interpolation_method", "linear"),
        ):
            ds.setncattr(name, value)


class TestPhysics(unittest.TestCase):
    def test_direct_factor_depends_on_zenith(self):
        factor = dict(side_factor=P.human_direct_side_factor, top_factor=P.human_horizontal_factor)
        self.assertAlmostEqual(direct_body_factor(0., **factor), 0.06)
        self.assertGreater(direct_body_factor(np.pi / 3, **factor), direct_body_factor(0., **factor))
        self.assertEqual(direct_body_factor(np.pi / 2, **factor), 0.)

    def test_uniform_enclosure_and_saturated_wick(self):
        ta = np.array([298.15])
        pressure = np.array([101325.])
        e = saturation_pressure_pa(ta)
        qv = EPSILON * e / (pressure - (1 - EPSILON) * e)
        fields = {"sw_direct_normal": np.zeros(1)}
        for side in ("upface", "downface", "northface", "southface", "eastface", "westface"):
            fields[f"sw_nondirect_{side}"] = np.zeros(1)
            fields[f"lw_{side}"] = SIGMA * ta**4
        np.testing.assert_allclose(_mrt(fields, np.pi / 3), ta, atol=1e-10)
        globe, wick, converged = _sensors(ta, pressure, qv, np.ones(1), fields, np.pi / 3)
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

        mrt = _mrt(fields, np.deg2rad(60.0))

        # Independent evaluation of the six-plane standing-person equation:
        # weighted SW=225.2922678358 and weighted LW=395.6 W m-2.
        self.assertAlmostEqual(float(mrt[0]), 314.985841050193, places=9)

    def test_iso_wbgt_reference_values(self):
        """Reproduce the ISO 7243 examples documented by pythermalcomfort."""
        ta = np.array([293.15, 293.15])
        globe = np.array([305.15, 305.15])
        wick = np.array([298.15, 298.15])
        direct = np.array([0.0, 1.0])

        result = _wbgt(ta, globe, wick, direct)

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

        globe, wick, converged = _sensors(
            ta, pressure, qv, np.array([1.2]), fields, np.deg2rad(60.0)
        )
        wbgt = _wbgt(ta, globe, wick, fields["sw_direct_normal"])

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
            source = Path(tmp) / EXCHANGE_NAME
            _exchange(source)
            output = calculate_indices(source)
            self.assertEqual(output.name, INDICES_NAME)
            with Dataset(output) as ds:
                ds.set_auto_mask(False)
                self.assertEqual(ds.indices, "mrt pet utci wbgt")
                self.assertEqual(ds["mrt"].dimensions, ("x", "y", "time"))
                self.assertEqual(float(ds["mrt"][0, 0, 0]), 25.)
                self.assertGreater(float(ds["mrt"][0, 1, 0]), 25.)
                self.assertGreater(float(ds["wbgt"][0, 1, 0]), float(ds["wbgt"][0, 0, 0]))
                self.assertTrue(np.isfinite(float(ds["pet"][0, 0, 0])))
                self.assertTrue(np.isfinite(float(ds["utci"][0, 0, 0])))
                self.assertEqual(float(ds["mrt"][1, 1, 0]), float(ds["mrt"]._FillValue))
                self.assertEqual(int(ds["validity_flags"][1, 1, 0]), -1)
            self.assertEqual(source.name, EXCHANGE_NAME)

    def test_published_pet_and_utci_reference_values(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / EXCHANGE_NAME
            _exchange(source, local_wind=0.1, direct_sun=False)
            output = calculate_indices(
                source, parameters=replace(P, pet_met=1.2, pet_clo=0.5)
            )
            with Dataset(output) as ds:
                self.assertAlmostEqual(float(ds["pet"][0, 0, 0]), 24.67, delta=0.05)
                self.assertAlmostEqual(float(ds["utci"][0, 0, 0]), 24.6, delta=0.1)
                self.assertEqual(int(ds["validity_flags"][0, 0, 0]) & 4, 4)

    def test_same_physics_for_different_file_names(self):
        with TemporaryDirectory() as tmp:
            first = Path(tmp) / "model_a.nc"
            second = Path(tmp) / "model_b.nc"
            _exchange(first)
            _exchange(second)
            with self.assertRaisesRegex(ValueError, "thermal_comfort_inputs"):
                calculate_indices(first)
            a = calculate_indices(first, output_path=Path(tmp) / "a.nc")
            b = calculate_indices(second, output_path=Path(tmp) / "b.nc")
            with Dataset(a) as left, Dataset(b) as right:
                for field in ("mrt", "pet", "utci", "wbgt"):
                    np.testing.assert_array_equal(left[field][:], right[field][:])

    def test_default_output_name_mirrors_the_input_name(self):
        self.assertEqual(default_output_path(Path("x/thermal_comfort_inputs.h1p5.007.nc")).name,
                         "thermal_comfort_indices.h1p5.007.nc")
        with self.assertRaises(ValueError):
            default_output_path(Path("x/other.nc"))

    def test_missing_value_inside_mask_is_flagged_not_fatal(self):
        import json
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / EXCHANGE_NAME
            _exchange(source)
            with Dataset(source, "r+") as ds:
                ds["qv"][0, 0, 0] = ds["qv"]._FillValue
            output = calculate_indices(source)
            with Dataset(output) as ds:
                ds.set_auto_mask(False)
                self.assertEqual(int(ds["validity_flags"][0, 0, 0]) & 64, 64)
                self.assertEqual(float(ds["mrt"][0, 0, 0]), float(ds["mrt"]._FillValue))
                # only that record and cell is lost
                self.assertTrue(np.isfinite(ds["mrt"][0, 0, 1]))
                self.assertTrue(np.isfinite(ds["mrt"][0, 1, 0]))
                counts = json.loads(ds.validity_counts)
                self.assertEqual(counts["invalid_input"], 1)

    def test_invalid_humidity_is_flagged_without_clipping_or_losing_mrt(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / EXCHANGE_NAME
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

def _generic_exchange(path: Path, *, n_times: int = 3, rh: float = 50.) -> None:
    """Minimal file written by another tool: the required variables only, no attributes beyond units."""
    ta = 298.15
    pressure = 101325.
    e = rh / 100. * float(saturation_pressure_pa(np.array([ta]))[0])
    qv = EPSILON * e / (pressure - (1 - EPSILON) * e)
    with Dataset(path, "w") as ds:
        for name, size in (("x", 2), ("y", 2), ("time", n_times), ("bounds", 2)):
            ds.createDimension(name, size)
        for name in ("x", "y"):
            ds.createVariable(name, "f8", (name,))[:] = [5., 15.]
            ds[name].units = "m"
            ds.createVariable(f"{name}_bounds", "f8", (name, "bounds"))[:] = [[0., 10.], [10., 20.]]
        times = 600.0 * (1 + np.arange(n_times))
        ds.createVariable("time", "f8", ("time",))[:] = times
        ds["time"].units = "s since simulation start"
        ds.createVariable("time_bounds", "f8", ("time", "bounds"))[:] = np.column_stack((times - 600., times))
        ds.createVariable("receptor_height", "f8").assignValue(1.5)
        ds.createVariable("pedestrian_mask", "i1", ("x", "y"))[:] = 1
        for name, value in (("solar_zenith", np.pi / 3), ("solar_azimuth", np.pi)):
            ds.createVariable(name, "f8", ("time",))[:] = value
            ds[name].units = "rad"
        values = {"ta": ta, "pabs": pressure, "qv": qv, "ws_local": 1., "ws_10": 1.,
                  "sw_direct_normal": 0.}
        for side in ("upface", "downface", "northface", "southface", "eastface", "westface"):
            values[f"sw_nondirect_{side}"] = 0.
            values[f"lw_{side}"] = SIGMA * ta**4
        for name, value in values.items():
            units = ("K" if name == "ta" else "Pa" if name == "pabs" else
                     "kg kg-1" if name == "qv" else "m s-1" if name.startswith("ws_")
                     else "W m-2")
            var = ds.createVariable(name, "f4", ("x", "y", "time"))
            data = np.full((2, 2, n_times), value, dtype=np.float32)
            if name == "ta":
                data[1, 1, :] = np.nan  # one cell with unusable input
            var[:] = data
            var.units = units


    def test_slight_supersaturation_is_capped_and_flagged(self):
        import json
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / EXCHANGE_NAME
            _exchange(source, rh=102.)
            output = calculate_indices(source)
            with Dataset(output) as ds:
                ds.set_auto_mask(False)
                flags = int(ds["validity_flags"][0, 0, 0])
                self.assertEqual(flags & 32, 32)
                self.assertEqual(flags & 1, 0)
                for name in ("mrt", "pet", "utci", "wbgt"):
                    self.assertNotEqual(float(ds[name][0, 0, 0]),
                                        float(ds[name]._FillValue))
                counts = json.loads(ds.validity_counts)
                self.assertGreater(counts["humidity_capped"], 0)
                self.assertEqual(counts["invalid_humidity"], 0)
            saturated = Path(tmp) / "thermal_comfort_inputs.h1p1.2.nc"
            _exchange(saturated, rh=100.)
            sat_out = calculate_indices(saturated)
            with Dataset(sat_out) as ds:
                ds.set_auto_mask(False)
                # exact saturation: storage round-off must not raise bit 32
                self.assertEqual(int(ds["validity_flags"][0, 0, 0]), 0)
                self.assertNotEqual(float(ds["wbgt"][0, 0, 0]),
                                    float(ds["wbgt"]._FillValue))
            strict = calculate_indices(
                source, output_path=Path(tmp) / "strict.nc",
                parameters=replace(P, supersaturation_tolerance_percent=0.0),
            )
            with Dataset(strict) as ds:
                ds.set_auto_mask(False)
                self.assertEqual(int(ds["validity_flags"][0, 0, 0]) & 1, 1)
                self.assertEqual(float(ds["pet"][0, 0, 0]), float(ds["pet"]._FillValue))


    def test_accepts_minimal_files_with_arbitrary_record_counts(self):
        import json
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / "minimal.nc"
            _generic_exchange(source, n_times=3)
            output = calculate_indices(source, output_path=Path(tmp) / "indices.nc")
            with Dataset(output) as ds:
                ds.set_auto_mask(False)
                self.assertEqual(len(ds.dimensions["time"]), 3)
                np.testing.assert_array_equal(ds["pedestrian_mask"][:], [[1, 1], [1, 1]])
                self.assertNotEqual(float(ds["mrt"][0, 0, 0]),
                                    float(ds["mrt"]._FillValue))
                # the cell with nonfinite input inside the mask is flagged, not fatal
                self.assertEqual(float(ds["mrt"][1, 1, 0]),
                                 float(ds["mrt"]._FillValue))
                self.assertEqual(int(ds["validity_flags"][1, 1, 0]) & 64, 64)
                counts = json.loads(ds.validity_counts)
                self.assertEqual(counts["invalid_input"], 3)

    def test_rejects_files_without_mask_or_bounds(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / "minimal.nc"
            out = Path(tmp) / "indices.nc"
            _generic_exchange(source, n_times=1)
            with Dataset(source, "r+") as ds:
                ds.renameVariable("pedestrian_mask", "other")
            with self.assertRaisesRegex(ValueError, "pedestrian_mask"):
                calculate_indices(source, output_path=out)
            _generic_exchange(source, n_times=1)
            with Dataset(source, "r+") as ds:
                ds.renameVariable("time_bounds", "other")
            with self.assertRaisesRegex(ValueError, "time_bounds"):
                calculate_indices(source, output_path=out)

    def test_flags_nonphysical_cells_instead_of_refusing(self):
        import json
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / "minimal.nc"
            _generic_exchange(source, n_times=2)
            with Dataset(source, "r+") as ds:
                ds["pabs"][0, 1, :] = -5.0  # nonphysical pressure in one cell
            output = calculate_indices(source, output_path=Path(tmp) / "indices.nc")
            with Dataset(output) as ds:
                ds.set_auto_mask(False)
                self.assertEqual(int(ds["validity_flags"][0, 1, 0]) & 64, 64)
                self.assertEqual(float(ds["mrt"][0, 1, 0]), float(ds["mrt"]._FillValue))
                self.assertTrue(np.isfinite(ds["mrt"][0, 0, 0]))
                counts = json.loads(ds.validity_counts)
                self.assertGreater(counts["invalid_input"], 0)

    def test_output_path_can_be_chosen(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / EXCHANGE_NAME
            _exchange(source)
            output = calculate_indices(source, output_path=Path(tmp) / "out.nc")
            with Dataset(output) as ds:
                ds.set_auto_mask(False)
                self.assertTrue(np.isfinite(ds["mrt"][0, 0, 0]))

    def test_beyond_tolerance_supersaturation_stays_invalid(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / EXCHANGE_NAME
            _exchange(source, rh=110.)
            output = calculate_indices(source)
            with Dataset(output) as ds:
                ds.set_auto_mask(False)
                flags = int(ds["validity_flags"][0, 0, 0])
                self.assertEqual(flags & 1, 1)
                self.assertEqual(flags & 32, 0)
                self.assertEqual(float(ds["utci"][0, 0, 0]), float(ds["utci"]._FillValue))
                self.assertTrue(np.isfinite(ds["mrt"][0, 0, 0]))

    def test_utci_low_wind_is_missing_and_flagged(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / EXCHANGE_NAME
            _exchange(source)
            with Dataset(source, "r+") as ds:
                ds["ws_10"][0, 0, :] = 0.1
            output = calculate_indices(source)
            with Dataset(output) as ds:
                ds.set_auto_mask(False)
                self.assertEqual(float(ds["utci"][0, 0, 0]), float(ds["utci"]._FillValue))
                self.assertEqual(int(ds["validity_flags"][0, 0, 0]) & 2, 2)
                self.assertTrue(np.isfinite(ds["pet"][0, 0, 0]))

    def test_direct_beam_below_horizon_flags_the_cell(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / EXCHANGE_NAME
            _exchange(source)
            with Dataset(source, "r+") as ds:
                ds["solar_zenith"][0] = np.pi / 2
                ds["sw_direct_normal"][0, 1, 0] = 100.
            output = calculate_indices(source)
            with Dataset(output) as ds:
                ds.set_auto_mask(False)
                # inconsistent cell flagged and filled; the rest computed
                self.assertEqual(int(ds["validity_flags"][0, 1, 0]) & 64, 64)
                self.assertEqual(float(ds["mrt"][0, 1, 0]), float(ds["mrt"]._FillValue))
                self.assertTrue(np.isfinite(ds["mrt"][0, 0, 0]))


class TestParameters(unittest.TestCase):
    def test_defaults_file_is_complete_and_strict(self):
        with DEFAULT_PARAMETERS_PATH.open(encoding="utf-8") as stream:
            raw = json.load(stream)
        self.assertEqual(set(raw), {field.name for field in fields(ComfortParameters)})
        self.assertEqual(P.pet_met, 1.37)
        self.assertEqual(P.wbgt_globe_diameter_m, 0.0508)
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "parameters.json"
            partial = dict(raw)
            partial.pop("pet_clo")
            path.write_text(json.dumps(partial), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Missing comfort parameters.*pet_clo"):
                ComfortParameters.from_file(path)
            extra = dict(raw, unknown=1)
            path.write_text(json.dumps(extra), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Unknown comfort parameters.*unknown"):
                ComfortParameters.from_file(path)
            changed = dict(raw, pet_clo=0.5)
            path.write_text(json.dumps(changed), encoding="utf-8")
            self.assertEqual(ComfortParameters.from_file(path).pet_clo, 0.5)

    def test_invalid_values_are_rejected(self):
        with self.assertRaises(ValueError):
            replace(P, human_horizontal_factor=0.1)
        with self.assertRaises(ValueError):
            replace(P, wbgt_sun_weights=(0.5, 0.5, 0.5))
        with self.assertRaises(ValueError):
            replace(P, utci_ta_min_c=60.0)
        with self.assertRaises(ValueError):
            replace(P, pet_sex="other")

    def test_index_selection_and_field_requirements(self):
        self.assertEqual(selected_indices(None), INDEX_NAMES)
        self.assertEqual(selected_indices(["wbgt", "mrt"]), ("mrt", "wbgt"))
        with self.assertRaisesRegex(ValueError, "Unknown comfort indices"):
            selected_indices(["pmv"])
        with self.assertRaises(ValueError):
            selected_indices([])
        self.assertNotIn("ta", required_fields(("mrt",)))
        self.assertIn("ws_10", required_fields(("utci",)))
        self.assertNotIn("ws_local", required_fields(("utci",)))
        self.assertIn("ws_local", required_fields(("pet",)))
        self.assertNotIn("ws_10", required_fields(("pet", "wbgt")))


@unittest.skipUnless(find_spec("pythermalcomfort"), "install tools/python[comfort]")
class TestIndexSelection(unittest.TestCase):
    def test_mrt_only_needs_no_atmosphere_fields(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / EXCHANGE_NAME
            _exchange(source)
            full = calculate_indices(source, output_path=Path(tmp) / "full.nc")
            with Dataset(source, "r+") as ds:
                for name in ("ta", "pabs", "qv", "ws_local", "ws_10"):
                    ds.renameVariable(name, f"gone_{name}")
            with self.assertRaisesRegex(ValueError, "lacks a valid ta"):
                calculate_indices(source, output_path=Path(tmp) / "all.nc")
            output = calculate_indices(source, indices=["mrt"])
            self.assertEqual(output.name, INDICES_NAME)
            with Dataset(output) as ds, Dataset(full) as reference:
                ds.set_auto_mask(False)
                self.assertEqual(ds.indices, "mrt")
                self.assertEqual(ds.comfort_method, "standing six-direction MRT")
                for name in ("pet", "utci", "wbgt", "globe_temperature",
                             "natural_wet_bulb_temperature"):
                    self.assertNotIn(name, ds.variables)
                self.assertFalse(hasattr(ds, "pet_utci_library"))
                np.testing.assert_array_equal(ds["mrt"][:], reference["mrt"][:])
                np.testing.assert_array_equal(ds["validity_flags"][0, 0, :], 0)
                counts = json.loads(ds.validity_counts)
                self.assertEqual(counts["utci_outside_domain"], 0)
                self.assertEqual(counts["wbgt_outside_domain"], 0)

    def test_wbgt_only_skips_mrt_and_ws_10(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / EXCHANGE_NAME
            _exchange(source)
            full = calculate_indices(source, output_path=Path(tmp) / "full.nc")
            with Dataset(source, "r+") as ds:
                ds.renameVariable("ws_10", "gone")
            output = calculate_indices(source, indices=("wbgt",))
            with Dataset(output) as ds, Dataset(full) as reference:
                ds.set_auto_mask(False)
                self.assertEqual(sorted(set(ds.variables) & set(INDEX_NAMES)), ["wbgt"])
                self.assertIn("globe_temperature", ds.variables)
                np.testing.assert_array_equal(ds["wbgt"][:], reference["wbgt"][:])
                np.testing.assert_array_equal(ds["natural_wet_bulb_temperature"][:],
                                              reference["natural_wet_bulb_temperature"][:])

    def test_utci_only_skips_ws_local(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / EXCHANGE_NAME
            _exchange(source)
            full = calculate_indices(source, output_path=Path(tmp) / "full.nc")
            with Dataset(source, "r+") as ds:
                ds.renameVariable("ws_local", "gone")
            output = calculate_indices(source, indices=["utci"])
            with Dataset(output) as ds, Dataset(full) as reference:
                ds.set_auto_mask(False)
                self.assertEqual(sorted(set(ds.variables) & set(INDEX_NAMES)), ["utci"])
                np.testing.assert_array_equal(ds["utci"][:], reference["utci"][:])
                # WBGT and PET bits are never raised when they are not selected
                self.assertEqual(int(ds["validity_flags"][0, 0, 0]) & (4 | 8 | 16), 0)

    def test_sensor_parameters_change_wbgt(self):
        with TemporaryDirectory() as tmp:
            source = Path(tmp) / EXCHANGE_NAME
            _exchange(source)
            base = calculate_indices(source, indices=["wbgt"], output_path=Path(tmp) / "a.nc")
            other = calculate_indices(
                source, indices=["wbgt"], output_path=Path(tmp) / "b.nc",
                parameters=replace(P, wbgt_sun_weights=(0.6, 0.3, 0.1)),
            )
            with Dataset(base) as a, Dataset(other) as b:
                a.set_auto_mask(False)
                b.set_auto_mask(False)
                # sunlit cell changes with the weights, shaded cell does not
                self.assertNotEqual(float(a["wbgt"][0, 1, 0]), float(b["wbgt"][0, 1, 0]))
                self.assertEqual(float(a["wbgt"][0, 0, 0]), float(b["wbgt"][0, 0, 0]))
                self.assertIn('"wbgt_sun_weights": [0.6, 0.3, 0.1]', b.comfort_parameters)

