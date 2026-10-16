"""Multi-height configuration and atmospheric-plane extraction tests."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from netCDF4 import Dataset
import f90nml
import numpy as np

from udbase import UDBase
from udcomf.heights import configured_receptor_heights, height_tag


def _case(root: Path, *, omit_second_level: bool = False) -> UDBase:
    case = root / "001"
    case.mkdir()
    (case / "namoptions.001").write_text(
        "&DOMAIN\n itot=2\n jtot=2\n ktot=4\n xlen=20.\n ylen=20.\n zsize=12.\n/\n"
        "&OUTPUT\n receptor_height=1.1\n nreceptor_heights=2\n"
        " receptor_heights=1.1,1.5\n tstatsdump=900.\n/\n",
        encoding="ascii",
    )
    (case / "prof.inp.001").write_text(
        "# profile\n# z thl\n0.75 300\n2.25 300\n9.75 300\n11.25 300\n",
        encoding="ascii",
    )
    sim = UDBase(1, case, load_geometry=False, suppress_load_warnings=True)
    sim.Sc = np.zeros((2, 2, 4), dtype=bool)
    sim.Sc[1, 1, 0] = True

    levels = [2.25] if omit_second_level else [0.75, 2.25]
    with Dataset(case / "stats_kslice.001.nc", "w") as ds:
        for name, length in (("time", 2), ("zt", len(levels)), ("yt", 2), ("xt", 2),
                             ("receptor_height", 2)):
            ds.createDimension(name, length)
        ds.createVariable("time", "f8", ("time",))[:] = [3600., 7200.]
        ds.createVariable("zt", "f4", ("zt",))[:] = levels
        ds.createVariable("xt", "f4", ("xt",))[:] = sim.xt
        ds.createVariable("yt", "f4", ("yt",))[:] = sim.yt
        ds.createVariable("receptor_height", "f4", ("receptor_height",))[:] = [1.1, 1.5]
        for name, values in (
            ("tha", [300., 306.]), ("pabs", [100000., 99000.]),
            ("qt", [0.012, 0.014]), ("ql", [0.002, 0.002]),
        ):
            data = np.broadcast_to(np.asarray(values[:len(levels)])[None, :, None, None],
                                   (2, len(levels), 2, 2)).copy()
            ds.createVariable(name, "f4", ("time", "zt", "yt", "xt"))[:] = data
        wind = np.broadcast_to(np.array([1., 2.])[None, :, None, None], (2, 2, 2, 2)).copy()
        ds.createVariable("ws_local", "f4", ("time", "receptor_height", "yt", "xt"))[:] = wind
        ds.createVariable("ws_10", "f4", ("time", "yt", "xt"))[:] = 3.
    return sim


class TestMultiheightComfort(unittest.TestCase):
    def test_namelist_height_list_and_tags(self):
        with TemporaryDirectory() as tmp:
            sim = _case(Path(tmp))
            self.assertEqual(configured_receptor_heights(sim.path), (1.1, 1.5))
            self.assertNotEqual(height_tag(1.1), height_tag(1.5))
            self.assertEqual(float(sim.receptor_height), 1.1)
            with self.assertRaisesRegex(ValueError, "strictly increasing"):
                sim.comf.radiation.write_hourly_heights(np.array([3600.]), heights=(1.5, 1.1))

    def test_atmospheric_files_select_saved_speed_and_interpolate_scalars(self):
        with TemporaryDirectory() as tmp:
            sim = _case(Path(tmp))
            original = sim.receptor_height
            paths = sim.comf.atmosphere.write_hourly_heights(
                np.array([3600., 7200.]), source="stats_kslice"
            )
            self.assertEqual(set(paths), {1.1, 1.5})
            self.assertEqual(sim.receptor_height, original)
            for height, path in paths.items():
                with Dataset(path) as ds:
                    ds.set_auto_mask(False)
                    self.assertEqual(ds["ta"].dimensions, ("x", "y", "time"))
                    self.assertAlmostEqual(float(ds["receptor_height"][...]), height)
                    np.testing.assert_allclose(ds["time_bounds"][:],
                                               [[2700., 3600.], [6300., 7200.]])
                    expected_ta = 300. + (height - 0.75) / 1.5 * 6.
                    self.assertAlmostEqual(float(ds["ta"][0, 0, 0]), expected_ta, places=3)
                    self.assertAlmostEqual(float(ds["qv"][0, 0, 0]),
                                           0.010 + (height - 0.75) / 1.5 * 0.002, places=6)
                    self.assertEqual(float(ds["ws_local"][0, 0, 0]), 1. if height == 1.1 else 2.)
                    self.assertEqual(float(ds["ws_10"][0, 0, 0]), 3.)
                    self.assertEqual(int(ds["valid_mask"][1, 1]), 0)
                    self.assertTrue(np.isnan(ds["ta"][1, 1, 0]))
            with self.assertRaises(FileExistsError):
                sim.comf.atmosphere.write_hourly_heights(np.array([3600., 7200.]))

    def test_missing_brackets_and_times_fail_without_outputs(self):
        with TemporaryDirectory() as tmp:
            sim = _case(Path(tmp), omit_second_level=True)
            with self.assertRaisesRegex(ValueError, "bracketing scalar levels"):
                sim.comf.atmosphere.write_hourly_heights(np.array([3600.]))
            self.assertFalse(list(sim.path.glob("pedestrian_atmosphere.*.nc")))
        with TemporaryDirectory() as tmp:
            sim = _case(Path(tmp))
            with self.assertRaisesRegex(ValueError, "unique record"):
                sim.comf.atmosphere.write_hourly_heights(np.array([3900.]))

    def test_rejects_non_15_minute_statistics(self):
        with TemporaryDirectory() as tmp:
            sim = _case(Path(tmp))
            path = sim.path / "namoptions.001"
            namelist = f90nml.read(path)
            namelist["output"]["tstatsdump"] = 600.0
            namelist.write(path, force=True)
            with self.assertRaisesRegex(ValueError, "tstatsdump must be 900"):
                sim.comf.atmosphere.write_hourly_heights(np.array([3600.]))

    def test_does_not_clip_small_negative_qv(self):
        with TemporaryDirectory() as tmp:
            sim = _case(Path(tmp))
            with Dataset(sim.path / "stats_kslice.001.nc", "r+") as ds:
                ds["qt"][0, :, 0, 0] = 0.001
                ds["ql"][0, :, 0, 0] = 0.001000001
            with self.assertRaisesRegex(ValueError, "qv must be nonnegative"):
                sim.comf.atmosphere.write_hourly_heights(np.array([3600.]))
            self.assertFalse(list(sim.path.glob("pedestrian_atmosphere.*.nc")))

    def test_full_stats_selects_only_bracketing_levels(self):
        with TemporaryDirectory() as tmp:
            sim = _case(Path(tmp))
            source = sim.path / "stats_kslice.001.nc"
            target = sim.path / "stats_t.001.nc"
            with Dataset(source) as src, Dataset(target, "w") as dst:
                for name, dimension in src.dimensions.items():
                    dst.createDimension(name, len(sim.zt) if name == "zt" else len(dimension))
                for name in ("time", "xt", "yt", "receptor_height"):
                    original = src[name]
                    dst.createVariable(name, original.dtype, original.dimensions)[:] = original[:]
                dst.createVariable("zt", "f4", ("zt",))[:] = sim.zt
                for name in ("tha", "pabs", "qt", "ql"):
                    values = np.zeros((2, len(sim.zt), 2, 2), dtype=np.float32)
                    values[:, :2] = src[name][:]
                    values[:, 2:] = 9999.
                    dst.createVariable(name, "f4", ("time", "zt", "yt", "xt"))[:] = values
                for name in ("ws_local", "ws_10"):
                    original = src[name]
                    dst.createVariable(name, original.dtype, original.dimensions)[:] = original[:]

            selected = sim.load_stat_t("tha", time_index=0, vertical_indices=[0, 1])
            self.assertEqual(selected.shape, (2, 2, 2))
            paths = sim.comf.atmosphere.write_hourly_heights(
                np.array([3600.]), source="stats_t", heights=(1.5,)
            )
            with Dataset(paths[1.5]) as ds:
                self.assertAlmostEqual(float(ds["ta"][0, 0, 0]), 303., places=3)

    def test_legacy_single_height_stats_layout(self):
        with TemporaryDirectory() as tmp:
            sim = _case(Path(tmp))
            options_path = sim.path / "namoptions.001"
            namelist = f90nml.read(options_path)
            del namelist["output"]["nreceptor_heights"]
            del namelist["output"]["receptor_heights"]
            namelist.write(options_path, force=True)
            source_path = sim.path / "stats_kslice.001.nc"
            replacement = sim.path / "single.nc"
            with Dataset(source_path) as src, Dataset(replacement, "w") as dst:
                for name, dimension in src.dimensions.items():
                    if name != "receptor_height":
                        dst.createDimension(name, len(dimension))
                for name in ("time", "zt", "xt", "yt", "tha", "pabs", "qt", "ql", "ws_10"):
                    variable = src[name]
                    dst.createVariable(name, variable.dtype, variable.dimensions)[:] = variable[:]
                dst.createVariable("receptor_height", "f4").assignValue(1.1)
                dst.createVariable("ws_local", "f4", ("time", "yt", "xt"))[:] = src["ws_local"][:, 0]
            replacement.replace(source_path)
            paths = sim.comf.atmosphere.write_hourly_heights(np.array([3600.]))
            self.assertEqual(set(paths), {1.1})
            with Dataset(paths[1.1]) as ds:
                self.assertEqual(float(ds["ws_local"][0, 0, 0]), 1.0)
