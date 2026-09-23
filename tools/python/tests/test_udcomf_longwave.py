"""Longwave source and receptor tests using small synthetic scenes."""

import types
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock

import numpy as np
import trimesh
import xarray as xr

from udbase import UDBase
from udcomf.udcomf_radiation import LongwaveState, UDComfRadiation


def _case(mesh=None):
    return types.SimpleNamespace(
        xazimuth=90.0,
        receptor_height=1.1,
        xt=np.array([0.0, 1.0]),
        yt=np.array([0.0, 1.0]),
        zm=np.array([0.0, 1.0, 2.0]),
        zsize=3.0,
        Sc=np.zeros((2, 2, 3), dtype=bool),
        geom=types.SimpleNamespace(stl=mesh) if mesh is not None else None,
        facs={"typeid": np.ones(len(mesh.faces), dtype=int)} if mesh is not None else {},
    )


def _ground():
    return trimesh.Trimesh(
        vertices=[[-1000, -1000, 0], [1000, -1000, 0],
                  [1000, 1000, 0], [-1000, 1000, 0]],
        faces=[[0, 1, 2], [0, 2, 3]], process=False,
    )


def _state(sky=350.0, exitance=()):
    return LongwaveState(
        time=300.0, sky_irradiance=sky,
        facet_exitance=np.asarray(exitance, dtype=float),
    )


class TestLongwaveInputs(unittest.TestCase):
    def test_reads_one_emitted_flux_record_and_interpolates_sky(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp)
            (path / "namoptions.001").write_text(
                "&DOMAIN\n itot=2\n jtot=2\n ktot=2\n"
                " xlen=20\n ylen=20\n zsize=20\n nfcts=2\n/\n",
                encoding="ascii",
            )
            (path / "timedeplw.inp.001").write_text(
                "Longwave forcing\ntime LWsky\n0 300\n600 400\n",
                encoding="ascii",
            )
            xr.Dataset(
                {"t": ("time", np.array([100.0, 300.0])),
                 "LWout": (("time", "fct"), np.array([[310.0, 320.0], [330.0, 340.0]])),
                 "LWin": (("time", "fct"), np.array([[1.0, 2.0], [3.0, 4.0]]))},
                coords={"time": np.arange(2), "fct": np.arange(2)},
            ).to_netcdf(path / "facEB.001.nc")
            sim = UDBase(1, path, load_geometry=False, suppress_load_warnings=True)
            with mock.patch.object(sim, "load_fac_eb", wraps=sim.load_fac_eb) as loader:
                state = sim.comf.radiation.load_longwave_state(1)
            self.assertEqual(state.time, 300.0)
            self.assertEqual(state.sky_irradiance, 350.0)
            np.testing.assert_array_equal(state.facet_exitance, [330.0, 340.0])
            loader.assert_any_call("LWout", time_index=1)

    def test_rejects_uncovered_forcing_time_and_bad_emission(self):
        sim = _case()
        sim.load_fac_eb = mock.Mock(side_effect=lambda var, **kw: (
            np.array([300.0]) if var == "t" else np.array([320.0])
        ))
        sim.load_timedeplw = mock.Mock(return_value={
            "time": np.array([0.0, 200.0]), "LWsky": np.array([300.0, 400.0]),
        })
        with self.assertRaisesRegex(ValueError, "outside saved sky"):
            UDComfRadiation(sim).load_longwave_state(0)
        sim.load_timedeplw.return_value["time"] = np.array([0.0, 600.0])
        sim.load_fac_eb.side_effect = lambda var, **kw: (
            np.array([300.0]) if var == "t" else np.array([-1.0])
        )
        with self.assertRaisesRegex(ValueError, "exitance"):
            UDComfRadiation(sim).load_longwave_state(0)


class TestLongwaveReceptors(unittest.TestCase):
    def test_open_sky_over_emitting_ground(self):
        radiation = UDComfRadiation(_case(_ground()))
        result = radiation.longwave_at_receptor(
            np.array([0.0, 0.0, 1.1]), _state(exitance=[400.0, 400.0])
        )
        self.assertEqual(set(result), {
            "lw_upface", "lw_downface", "lw_northface", "lw_southface",
            "lw_eastface", "lw_westface",
        })
        self.assertAlmostEqual(result["lw_upface"], 350.0, places=10)
        self.assertAlmostEqual(result["lw_downface"], 400.0, places=10)
        for name in ("northface", "southface", "eastface", "westface"):
            self.assertAlmostEqual(result[f"lw_{name}"], 375.0, delta=2.0)

    def test_rejects_unmeshed_ground(self):
        with self.assertRaisesRegex(ValueError, "unmeshed ground"):
            UDComfRadiation(_case()).longwave_at_receptor(
                np.array([0.0, 0.0, 1.1]), _state()
            )

    def test_emitting_wall_blocks_sky_and_ray_map_reuses_hits(self):
        wall = trimesh.creation.box(extents=(0.2, 2.0, 2.0))
        wall.apply_translation((1.0, 0.0, 1.0))
        scene = trimesh.util.concatenate((wall, _ground()))
        radiation = UDComfRadiation(_case(scene))
        point = np.array([0.0, 0.0, 1.0])
        ray_map = radiation.trace_shortwave_rays(point, n_mu=4, n_azimuth=16)
        exitance = np.full(len(scene.faces), 420.0)
        with mock.patch.object(radiation, "_directional_source", side_effect=AssertionError("re-traced")):
            first = radiation.longwave_at_receptor(
                point, _state(sky=300.0, exitance=exitance),
                n_mu=4, n_azimuth=16, ray_map=ray_map,
            )
            second = radiation.longwave_at_receptor(
                point, _state(sky=300.0, exitance=exitance * 2),
                n_mu=4, n_azimuth=16, ray_map=ray_map,
            )
        self.assertGreater(first["lw_eastface"], first["lw_westface"])
        self.assertGreater(second["lw_eastface"], first["lw_eastface"])
        self.assertGreater(first["lw_upface"], 0.0)
        cold = radiation.longwave_at_receptor(
            point, _state(sky=300.0, exitance=np.zeros(len(scene.faces))),
            n_mu=4, n_azimuth=16, ray_map=ray_map,
        )
        self.assertLess(cold["lw_upface"], 300.0)
        self.assertGreater(first["lw_upface"], cold["lw_upface"])
        with self.assertRaisesRegex(ValueError, "does not match receptor"):
            radiation.longwave_at_receptor(
                np.array([0.0, 1.0, 1.0]), _state(exitance=exitance),
                n_mu=4, n_azimuth=16, ray_map=ray_map,
            )

    def test_sky_and_uniform_lower_emission_integrate_to_expected_planes(self):
        radiation = UDComfRadiation(_case())
        with (
            mock.patch.object(radiation, "_mesh", return_value=types.SimpleNamespace(n_cells=1)),
            mock.patch.object(
                radiation, "_directional_source",
                side_effect=lambda _point, direction: -1 if direction[2] > 0 else 0,
            ),
        ):
            result = radiation.longwave_at_receptor(
                np.array([0.0, 0.0, 1.1]), _state(sky=300.0, exitance=[400.0])
            )
        self.assertAlmostEqual(result["lw_upface"], 300.0, places=10)
        self.assertAlmostEqual(result["lw_downface"], 400.0, places=10)
        for name in ("northface", "southface", "eastface", "westface"):
            self.assertAlmostEqual(result[f"lw_{name}"], 350.0, delta=2.0)

    def test_uniform_temperature_enclosure_gives_emitted_irradiance(self):
        radiation = UDComfRadiation(_case())
        exitance = 0.9 * 5.670374419e-8 * 300.0**4
        with (
            mock.patch.object(radiation, "_mesh", return_value=types.SimpleNamespace(n_cells=1)),
            mock.patch.object(radiation, "_directional_source", return_value=0),
        ):
            result = radiation.longwave_at_receptor(
                np.array([0.0, 0.0, 1.1]), _state(sky=0.0, exitance=[exitance])
            )
        for value in result.values():
            self.assertAlmostEqual(value, exitance, delta=0.005 * exitance)

    def test_native_plane_marks_solid_and_unselected_cells_missing(self):
        sim = _case(_ground())
        sim.Sc[1, 1, 1] = True
        fields = UDComfRadiation(sim).longwave_plane(
            _state(exitance=[400.0, 400.0]), flat_indices=np.array([0, 2]),
            n_mu=2, n_azimuth=8,
        )
        self.assertEqual(len(fields), 6)
        for field in fields.values():
            self.assertEqual(field.shape, (2, 2))
            self.assertTrue(np.isfinite(field[0, 0]))
            self.assertTrue(np.isfinite(field[1, 0]))
            self.assertTrue(np.isnan(field[0, 1]))
            self.assertTrue(np.isnan(field[1, 1]))

    def test_vegetation_case_does_not_silently_omit_tree_longwave(self):
        sim = _case()
        sim.ltrees = True
        with self.assertRaisesRegex(NotImplementedError, "vegetation"):
            UDComfRadiation(sim).longwave_at_receptor(
                np.array([0.0, 0.0, 1.1]), _state()
            )
