"""Hourly integration and resumable pedestrian-radiation output tests."""

from pathlib import Path
from tempfile import TemporaryDirectory
import types
import unittest
from unittest import mock

from netCDF4 import Dataset
import numpy as np
import trimesh

from udcomf.checkpoints import hourly_windows
from udcomf.udcomf_radiation import UDComfRadiation


def _ground():
    return trimesh.Trimesh(
        vertices=[[-100, -100, 0], [100, -100, 0],
                  [100, 100, 0], [-100, 100, 0]],
        faces=[[0, 1, 2], [0, 2, 3]], process=False,
    )


def _case(path):
    mesh = _ground()
    times = np.arange(0, 7201, 300, dtype=float)
    solid = np.zeros((2, 2, 2), dtype=bool)
    solid[1, 1, 1] = True
    sim = types.SimpleNamespace(
        path=path, expnr=1, stl_file="ground.stl", receptor_height=1.1,
        xt=np.array([0.0, 1.0]), yt=np.array([0.0, 1.0]),
        zm=np.array([0.0, 1.0]), zsize=2.0, Sc=solid, xazimuth=0.0,
        geom=types.SimpleNamespace(stl=mesh),
        facs={"typeid": np.ones(2, dtype=int)},
        ltrees=False,
    )
    sw = {
        "time": times,
        "dni": np.zeros(len(times)),
        "dsky": 100.0 + 0.01 * times,
        "solar_zenith": np.full(len(times), 110.0),
        "solar_azimuth_local": np.zeros(len(times)),
    }
    sim.load_shortwave_forcing = lambda name: sw[name]
    sim.load_timedepsw = lambda: {
        "time": times, "netsw": np.full((2, len(times)), 50.0),
    }
    sim.assign_prop_to_fac = lambda name: np.full(2, 0.2)
    sim.load_fac_eb = lambda name: (
        times if name == "t" else np.full((2, len(times)), 400.0)
    )
    sim.load_timedeplw = lambda: {
        "time": np.array([0.0, 7200.0]),
        "LWsky": np.array([300.0, 300.0]),
    }
    for name in (
        "namoptions.1", "facets.inp.1", "factypes.inp.1",
        "solid_c.txt", "ground.stl", "shortwave_forcing.1.nc",
        "timedepsw.inp.1", "facEB.1.nc", "timedeplw.inp.1",
    ):
        (path / name).write_bytes(b"fixture")
    return sim


class TestHourlyWindows(unittest.TestCase):
    def test_linear_series_integrates_exact_preceding_quarter_hour(self):
        source = np.arange(0, 3901, 300, dtype=float)
        windows, gap = hourly_windows(source, np.array([3600.0]))
        self.assertEqual(gap, 330.0)
        self.assertTrue(windows[0].complete)
        self.assertAlmostEqual(windows[0].weights.sum(), 1.0)
        # f(t)=t: the [2700, 3600) mean is 3150.
        self.assertAlmostEqual(source[windows[0].indices] @ windows[0].weights, 3150.0)

    def test_gap_and_uncovered_endpoints_are_missing(self):
        source = np.array([0, 2700, 3000, 3600, 3900, 6300, 6600, 6900], dtype=float)
        windows, _ = hourly_windows(source, np.array([3600.0, 7200.0]), max_gap=450)
        self.assertFalse(windows[0].complete)
        self.assertFalse(windows[1].complete)
        with self.assertRaisesRegex(ValueError, "exact-hour"):
            hourly_windows(source, np.array([3500.0]))


class TestHourlyRadiation(unittest.TestCase):
    def test_multi_height_outputs_and_checkpoints_are_independent(self):
        with TemporaryDirectory() as tmp:
            sim = _case(Path(tmp))
            radiation = UDComfRadiation(sim)
            outputs = radiation.write_hourly_heights(
                np.array([3600.]), heights=(1.1, 1.5), kinds=("shortwave",),
                flat_indices=np.array([0]), tile_size=1, n_mu=2, n_azimuth=8,
            )
            self.assertEqual(sim.receptor_height, 1.1)
            self.assertNotEqual(outputs[1.1]["shortwave"], outputs[1.5]["shortwave"])
            for height, files in outputs.items():
                with Dataset(files["shortwave"]) as ds:
                    self.assertAlmostEqual(float(ds["receptor_height"][...]), height)
                tag = f"h{str(height).replace('.', 'p')}"
                self.assertTrue((Path(tmp) / "udcomf_radiation.checkpoints" / "shortwave"
                                 / tag / "time_0000_tile_00000.npz").is_file())
            with self.assertRaises(FileExistsError):
                radiation.write_hourly_heights(
                    np.array([3600.]), heights=(1.1, 1.5), kinds=("shortwave",),
                    flat_indices=np.array([0]), tile_size=1, n_mu=2, n_azimuth=8,
                )

    def test_native_grid_files_and_resume_do_not_retrace(self):
        with TemporaryDirectory() as tmp:
            sim = _case(Path(tmp))
            radiation = UDComfRadiation(sim)
            targets = np.array([3600.0, 7200.0])
            for kind, names, expected in (
                ("shortwave", ("sw_direct_normal", "sw_nondirect_upface"), 131.5),
                ("longwave", ("lw_upface", "lw_downface"), 400.0),
            ):
                output = radiation.write_hourly(
                    kind, targets, flat_indices=np.array([0, 1]),
                    tile_size=1, n_mu=2, n_azimuth=8,
                )
                with Dataset(output) as ds:
                    ds.set_auto_mask(False)
                    self.assertEqual(ds.variables[names[0]].dimensions, ("x", "y", "time"))
                    np.testing.assert_array_equal(ds.variables["time_complete"][:], [1, 1])
                    np.testing.assert_array_equal(
                        ds.variables["time_bounds"][:], [[2700, 3600], [6300, 7200]],
                    )
                    self.assertEqual(float(ds.variables["receptor_height"][...]), 1.1)
                    self.assertEqual(ds.variables["processed_mask"][0, 0], 1)
                    self.assertEqual(ds.variables["processed_mask"][1, 1], 0)
                    self.assertTrue(np.isnan(ds.variables[names[0]][1, 1, 0]))
                    self.assertAlmostEqual(float(ds.variables[names[-1]][0, 0, 0]), expected)
                checkpoint = (
                    Path(tmp) / "udcomf_radiation.checkpoints" / kind
                    / "time_0001_tile_00000.npz"
                )
                checkpoint.unlink()
                with mock.patch.object(radiation, "_source_series", wraps=radiation._source_series) as series:
                    radiation.write_hourly(
                        kind, targets, flat_indices=np.array([0, 1]),
                        tile_size=1, n_mu=2, n_azimuth=8,
                        resume=True, overwrite=True,
                    )
                    self.assertEqual(series.call_count, 1)
                with mock.patch.object(radiation, "_source_series", side_effect=AssertionError("recomputed")):
                    radiation.write_hourly(
                        kind, targets, flat_indices=np.array([0, 1]),
                        tile_size=1, n_mu=2, n_azimuth=8,
                        resume=True, overwrite=True,
                    )
                with self.assertRaisesRegex(ValueError, "settings changed"):
                    radiation.write_hourly(
                        kind, targets, flat_indices=np.array([0, 1]),
                        tile_size=1, n_mu=4, n_azimuth=8, overwrite=True,
                    )

    def test_missing_window_outputs_nan_without_tracing(self):
        with TemporaryDirectory() as tmp:
            sim = _case(Path(tmp))
            radiation = UDComfRadiation(sim)
            with mock.patch.object(radiation, "_source_series", side_effect=AssertionError("traced")):
                path = radiation.write_hourly(
                    "longwave", np.array([10800.0]),
                    flat_indices=np.array([0]), n_mu=2, n_azimuth=8,
                )
            with Dataset(path) as ds:
                ds.set_auto_mask(False)
                self.assertEqual(ds.variables["time_complete"][0], 0)
                self.assertTrue(np.isnan(ds.variables["lw_upface"][0, 0, 0]))
