"""Physics and data-contract tests for pedestrian shortwave irradiance."""

import types
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock

import numpy as np
import trimesh

from udbase import UDBase
from udprep.shortwave_forcing import write_shortwave_forcing
from udcomf.udcomf_radiation import (
    ShortwaveState,
    UDComfRadiation,
    _sphere_quadrature,
    facet_shortwave_exitance,
)


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


def _state(*, dni=0.0, dsky=0.0, zenith=60.0, azimuth=0.0, exitance=()):
    return ShortwaveState(
        time=0.0, dni=dni, dsky=dsky, zenith=zenith,
        azimuth_local=azimuth, facet_exitance=np.asarray(exitance, dtype=float),
    )


class TestShortwaveInputs(unittest.TestCase):
    def test_reads_real_archive_and_timedep_file_through_udbase(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp)
            (path / "namoptions.001").write_text(
                "&DOMAIN\n itot = 2\n jtot = 2\n ktot = 2\n"
                " xlen = 20.0\n ylen = 20.0\n zsize = 20.0\n nfcts = 1\n/\n",
                encoding="ascii",
            )
            (path / "timedepsw.inp.001").write_text(
                "# net shortwave\n0.00 600.00\n0.0000 80.0000\n", encoding="ascii"
            )
            write_shortwave_forcing(
                path / "shortwave_forcing.001.nc",
                np.array([0.0, 600.0]), np.array([0.0, 700.0]),
                np.array([0.0, 100.0]), np.array([95.0, 35.0]),
                np.array([0.0, 42.0]),
                source="test", start_time="2023-08-21T00:00:00",
            )
            sim = UDBase(1, path, load_geometry=False, suppress_load_warnings=True)
            sim._lffacets = sim._lffactypes = True
            sim.facs["typeid"] = np.array([2])
            sim.factypes["id"] = np.array([2])
            sim.factypes["al"] = np.array([0.2])

            state = sim.comf.radiation.load_shortwave_state(1)

            self.assertEqual((state.time, state.dni, state.dsky), (600.0, 700.0, 100.0))
            np.testing.assert_allclose(state.facet_exitance, [20.0])

    def test_reflected_exitance_uses_absorbed_flux_and_albedo(self):
        np.testing.assert_allclose(
            facet_shortwave_exitance(np.array([80.0, 50.0]), np.array([0.2, 0.0])),
            [20.0, 0.0],
        )
        for albedo in (-0.1, 1.0, np.nan):
            with self.subTest(albedo=albedo), self.assertRaisesRegex(ValueError, "albedo"):
                facet_shortwave_exitance(np.array([80.0]), np.array([albedo]))
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            facet_shortwave_exitance(np.array([-1.0]), np.array([0.2]))

    def test_state_reads_existing_udbase_loaders_and_checks_time(self):
        forcing = {
            "time": np.array([0.0, 600.0]),
            "dni": np.array([0.0, 700.0]),
            "dsky": np.array([0.0, 100.0]),
            "solar_zenith": np.array([95.0, 35.0]),
            "solar_azimuth_local": np.array([0.0, 42.0]),
        }
        sim = _case()
        sim.load_shortwave_forcing = mock.Mock(side_effect=lambda name: forcing[name])
        sim.load_timedepsw = mock.Mock(return_value={
            "time": 600.0, "netsw": np.array([80.0, 50.0])
        })
        sim.assign_prop_to_fac = mock.Mock(return_value=np.array([0.2, 0.0]))
        state = UDComfRadiation(sim).load_shortwave_state(1)
        self.assertEqual((state.time, state.dni, state.dsky), (600.0, 700.0, 100.0))
        np.testing.assert_allclose(state.facet_exitance, [20.0, 0.0])
        sim.load_timedepsw.assert_called_once_with(time_index=1)
        sim.assign_prop_to_fac.assert_called_once_with("al")
        sim.load_timedepsw.return_value["time"] = 300.0
        with self.assertRaisesRegex(ValueError, "timestamps do not match"):
            UDComfRadiation(sim).load_shortwave_state(1)

    def test_quadrature_integrates_full_sphere_solid_angle(self):
        directions, weights = _sphere_quadrature(8, 32)
        np.testing.assert_allclose(np.linalg.norm(directions, axis=1), 1.0, atol=1e-14)
        self.assertAlmostEqual(float(weights.sum()), 4.0 * np.pi, places=12)
        self.assertAlmostEqual(
            float(np.dot(np.maximum(directions[:, 2], 0), weights)), np.pi, places=12
        )


class TestShortwaveReceptors(unittest.TestCase):
    def test_open_sky_direct_is_separate_from_sky_diffuse(self):
        radiation = UDComfRadiation(_case())
        point = np.array([0.0, 0.0, 1.1])
        result = radiation.shortwave_at_receptor(point, _state(dni=750, dsky=100))
        self.assertEqual(result["sw_direct_normal"], 750.0)
        self.assertAlmostEqual(result["sw_nondirect_upface"], 100.0, places=10)
        self.assertEqual(result["sw_nondirect_downface"], 0.0)
        for name in ("northface", "southface", "eastface", "westface"):
            self.assertAlmostEqual(result[f"sw_nondirect_{name}"], 50.0, delta=1.0)
        no_direct = radiation.shortwave_at_receptor(point, _state(dsky=100))
        for name in result:
            if name != "sw_direct_normal":
                self.assertAlmostEqual(result[name], no_direct[name], places=12)
        direct_only = radiation.shortwave_at_receptor(point, _state(dni=750))
        self.assertEqual(direct_only["sw_direct_normal"], 750.0)
        for name, value in direct_only.items():
            if name != "sw_direct_normal":
                self.assertEqual(value, 0.0)
        self.assertEqual(
            radiation.shortwave_at_receptor(point, _state(dni=750, zenith=95))["sw_direct_normal"],
            0.0,
        )
        with self.assertRaisesRegex(ValueError, "Solar zenith"):
            radiation.shortwave_at_receptor(point, _state(zenith=200))

    def test_uniform_ground_and_sky_radiance_integrate_to_expected_planes(self):
        radiation = UDComfRadiation(_case())
        with (
            mock.patch.object(radiation, "_mesh", return_value=types.SimpleNamespace(n_cells=1)),
            mock.patch.object(
                radiation, "_directional_source",
                side_effect=lambda _point, direction: -1 if direction[2] > 0 else 0,
            ),
        ):
            result = radiation.shortwave_at_receptor(
                np.array([0.0, 0.0, 1.1]), _state(dsky=100, exitance=[40])
            )
        self.assertAlmostEqual(result["sw_nondirect_upface"], 100.0, places=10)
        self.assertAlmostEqual(result["sw_nondirect_downface"], 40.0, places=10)
        for name in ("northface", "southface", "eastface", "westface"):
            self.assertAlmostEqual(result[f"sw_nondirect_{name}"], 70.0, delta=1.0)

    def test_reflecting_wall_shades_direct_sun(self):
        wall = trimesh.creation.box(extents=(0.2, 2.0, 2.0))
        wall.apply_translation((1.0, 0.0, 1.0))
        radiation = UDComfRadiation(_case(wall))
        result = radiation.shortwave_at_receptor(
            np.array([0.0, 0.0, 1.0]),
            _state(dni=700, exitance=np.full(len(wall.faces), 80.0)),
            n_mu=8, n_azimuth=32,
        )
        self.assertEqual(result["sw_direct_normal"], 0.0)
        self.assertGreater(result["sw_nondirect_eastface"], 0.0)
        self.assertEqual(result["sw_nondirect_westface"], 0.0)
        self.assertTrue(all(value >= 0 for value in result.values()))

    def test_back_face_blocks_sky_without_emitting_its_outward_flux(self):
        roof = trimesh.Trimesh(
            vertices=[[-2, -2, 2], [2, -2, 2], [0, 2, 2]],
            faces=[[0, 1, 2]], process=False,
        )
        radiation = UDComfRadiation(_case(roof))
        self.assertEqual(
            radiation._directional_source(np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 1.0])),
            -2,
        )

    def test_vegetation_case_fails_instead_of_ignoring_attenuation(self):
        sim = _case()
        sim.ltrees = True
        with self.assertRaisesRegex(NotImplementedError, "vegetation attenuation"):
            UDComfRadiation(sim).shortwave_at_receptor(
                np.array([0.0, 0.0, 1.1]), _state(dsky=100)
            )

    def test_native_plane_marks_solid_and_unselected_cells_missing(self):
        sim = _case()
        sim.Sc[1, 1, 1] = True
        fields = UDComfRadiation(sim).shortwave_plane(
            _state(dni=500, dsky=100), flat_indices=np.array([0, 2]),
            n_mu=2, n_azimuth=8,
        )
        self.assertEqual(len(fields), 7)
        for field in fields.values():
            self.assertEqual(field.shape, (2, 2))
            self.assertTrue(np.isfinite(field[0, 0]))
            self.assertTrue(np.isfinite(field[1, 0]))
            self.assertTrue(np.isnan(field[0, 1]))
            self.assertTrue(np.isnan(field[1, 1]))

    def test_static_ray_map_reuses_visibility_across_timestamps(self):
        radiation = UDComfRadiation(_case())
        point = np.array([0.0, 0.0, 1.1])
        ray_map = radiation.trace_shortwave_rays(point, n_mu=2, n_azimuth=8)
        with mock.patch.object(radiation, "_directional_source", side_effect=AssertionError("re-traced")):
            first = radiation.shortwave_at_receptor(
                point, _state(dsky=100), n_mu=2, n_azimuth=8, ray_map=ray_map
            )
            second = radiation.shortwave_at_receptor(
                point, _state(dsky=200), n_mu=2, n_azimuth=8, ray_map=ray_map
            )
        self.assertAlmostEqual(second["sw_nondirect_upface"], 2 * first["sw_nondirect_upface"])
        with self.assertRaisesRegex(ValueError, "does not match receptor"):
            radiation.shortwave_at_receptor(
                np.array([1.0, 0.0, 1.1]), _state(dsky=100),
                n_mu=2, n_azimuth=8, ray_map=ray_map,
            )

        wall = trimesh.creation.box(extents=(0.2, 2.0, 2.0))
        wall.apply_translation((1.0, 0.0, 1.0))
        radiation = UDComfRadiation(_case(wall))
        wall_map = radiation.trace_shortwave_rays(point, n_mu=4, n_azimuth=16)
        first = radiation.shortwave_at_receptor(
            point, _state(exitance=np.full(len(wall.faces), 40.0)),
            n_mu=4, n_azimuth=16, ray_map=wall_map,
        )
        second = radiation.shortwave_at_receptor(
            point, _state(exitance=np.full(len(wall.faces), 80.0)),
            n_mu=4, n_azimuth=16, ray_map=wall_map,
        )
        self.assertGreater(first["sw_nondirect_eastface"], 0.0)
        self.assertAlmostEqual(
            second["sw_nondirect_eastface"], 2.0 * first["sw_nondirect_eastface"]
        )
