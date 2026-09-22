"""Small end-to-end uDALES thermal-comfort postprocessing test."""

from datetime import date
from importlib.util import find_spec
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from netCDF4 import Dataset
import numpy as np
import trimesh

from udcomf import UDComf
from udcomf.thermalcomfort import calculate_indices
from udcomf.thermalcomfort.physics import SIGMA
from udcomf.udcomf_export import ExchangeMetadata


def _small_case(root: Path):
    case = root / "001"
    case.mkdir()
    (case / "namoptions.001").write_text(
        "&DOMAIN\n"
        " itot=1\n jtot=1\n ktot=2\n xlen=10.\n ylen=10.\n zsize=2.\n"
        " nfcts=2\n nfctypes=1\n nfaclyrs=1\n"
        " year=2023\n month=8\n day=21\n hour=0\n timezone=0.\n"
        " longitude=2.33\n latitude=48.86\n elevation=35.\n/\n"
        "&OUTPUT\n receptor_height=1.1\n tstatsdump=900.\n/\n",
        encoding="ascii",
    )
    (case / "prof.inp.001").write_text(
        "# z thl\n0.5 298.15\n1.5 298.15\n", encoding="ascii"
    )
    for name in (
        "stats_kslice.001.nc", "facets.inp.001", "factypes.inp.001",
        "solid_c.txt", "ground.stl", "shortwave_forcing.001.nc",
        "timedepsw.inp.001", "facEB.001.nc", "timedeplw.inp.001",
    ):
        (case / name).write_bytes(b"validation fixture")

    mesh = trimesh.Trimesh(
        vertices=[[-100.0, -100.0, 0.0], [100.0, -100.0, 0.0],
                  [100.0, 100.0, 0.0], [-100.0, 100.0, 0.0]],
        faces=[[0, 1, 2], [0, 2, 3]], process=False,
    )
    targets = np.arange(1, 25, dtype=float) * 3600.0
    source_times = np.arange(0.0, targets[-1] + 300.0, 300.0)
    temperature = 298.15
    thermal_flux = SIGMA * temperature**4
    sim = SimpleNamespace(
        path=case, expnr="001", receptor_height=1.1,
        xt=np.array([5.0]), yt=np.array([5.0]),
        zt=np.array([0.5, 1.5]), zm=np.array([0.0, 1.0]),
        dzt=np.array([1.0, 1.0]), zsize=2.0,
        dx=10.0, dy=10.0, xazimuth=0.0,
        Sc=np.zeros((1, 1, 2), dtype=bool),
        stl_file="ground.stl", geom=SimpleNamespace(stl=mesh),
        facs={"typeid": np.ones(2, dtype=int), "normals": mesh.face_normals},
        ltrees=False,
    )

    profiles = {
        "tha": np.array([temperature, temperature]),
        "pabs": np.array([101325.0, 101325.0]),
        "qt": np.array([0.010, 0.010]),
        "ql": np.array([0.0, 0.0]),
    }

    def load_stats(name, *, time_index=None, vertical_indices=None):
        if name == "time":
            return targets
        if name == "zt":
            return sim.zt
        if name == "xt":
            return sim.xt
        if name == "yt":
            return sim.yt
        if name == "receptor_height":
            return np.array([1.1])
        if name == "ws_local":
            return np.full((1, 1), 1.0)
        if name == "ws_10":
            return np.full((1, 1), 2.0)
        levels = np.arange(2) if vertical_indices is None else np.asarray(vertical_indices)
        return profiles[name][levels].reshape(1, 1, -1)

    shortwave = {
        "time": source_times,
        "dni": np.zeros_like(source_times),
        "dsky": np.zeros_like(source_times),
        "solar_zenith": np.full_like(source_times, 100.0),
        "solar_azimuth_local": np.zeros_like(source_times),
    }
    sim.load_stat_kslice = load_stats
    sim.load_stat_t = load_stats
    sim.load_shortwave_forcing = lambda name: shortwave[name]
    sim.load_timedepsw = lambda: {
        "time": source_times,
        "netsw": np.zeros((2, len(source_times))),
    }
    sim.assign_prop_to_fac = lambda _name: np.full(2, 0.2)
    sim.load_fac_eb = lambda name, **_kwargs: (
        source_times if name == "t"
        else np.full((2, len(source_times)), thermal_flux)
    )
    sim.load_timedeplw = lambda: {
        "time": np.array([0.0, targets[-1]]),
        "LWsky": np.array([thermal_flux, thermal_flux]),
    }
    sim.comf = UDComf(sim)
    return sim, targets


@unittest.skipUnless(find_spec("pythermalcomfort"), "install tools/python[comfort]")
class TestEndToEndComfort(unittest.TestCase):
    def test_small_udales_case_reaches_model_neutral_indices(self):
        with TemporaryDirectory() as tmp:
            sim, targets = _small_case(Path(tmp))
            sim.comf.atmosphere.write_hourly_heights(targets, source="stats_kslice")
            sim.comf.radiation.write_hourly_heights(
                targets, tile_size=1, n_mu=8, n_azimuth=32, resume=False,
            )
            metadata = ExchangeMetadata(
                x=np.array([600005.0]), y=np.array([6860005.0]),
                x_bounds=np.array([[600000.0, 600010.0]]),
                y_bounds=np.array([[6860000.0, 6860010.0]]),
                longitude=np.array([[2.33]]), latitude=np.array([[48.86]]),
                z_ground=np.array([[77.0]]),
                crs_attributes={
                    "grid_mapping_name": "transverse_mercator",
                    "crs_wkt": "synthetic validation CRS",
                },
                vertical_datum="IGN69", model_version="validation",
                institution="uDALES test suite",
                building_representation="flat open-ground synthetic scene",
                terrain_convention="flat model ground",
                spinup_start_utc="2023-08-21T00:00:00Z",
                spinup_end_utc="2023-08-21T00:00:00Z",
            )
            exchange = sim.comf.export_exchange(date(2023, 8, 21), metadata)[1.1]
            output = calculate_indices(exchange)

            with Dataset(exchange) as source, Dataset(output) as result:
                self.assertEqual(source["ta"].dimensions, ("x", "y", "time"))
                self.assertEqual(result["mrt"].dimensions, ("x", "y", "time"))
                np.testing.assert_allclose(source["qv"][0, 0, :], 0.010, atol=1e-7)
                np.testing.assert_allclose(result["mrt"][0, 0, :], 25.0, atol=2e-4)
                for name in ("pet", "utci", "wbgt"):
                    self.assertTrue(np.isfinite(result[name][0, 0, :]).all())
                np.testing.assert_array_equal(result["validity_flags"][0, 0, :], 0)


if __name__ == "__main__":
    unittest.main()
