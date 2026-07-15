"""Characterization tests for the user-facing loader/analysis API of ``udbase``.

These lock in the CURRENT behaviour of the previously untested loaders so a
future refactor cannot silently change return shapes, coordinate transposition,
or the hand-checkable facet-averaging maths. They intentionally assert on the
present behaviour (a regression net), not on what the API "ought" to do.

Covered methods:
  load_field, load_stat_xyt, load_stat_t, load_stat_tree, load_slice,
  load_fac_momentum, load_fac_eb, load_fac_temperature, load_seb,
  area_average_fac, area_average_seb, convert_facvar_to_field,
  convert_facflx_to_field, calculate_frontal_properties.
"""

import io
import unittest
import warnings
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np
import xarray as xr

from udbase import UDBase  # noqa: E402


# Domain used by every case: itot=4, jtot=3, ktot=2 -> dx=dy=10, dzt=[10,10].
_NAMOPTIONS = "\n".join(
    [
        "&DOMAIN",
        " itot = 4",
        " jtot = 3",
        " ktot = 2",
        " xlen = 40.0",
        " ylen = 30.0",
        " zsize = 20.0",
        "/",
        "&WALLS",
        " nfaclyrs = 3",
        "/",
    ]
) + "\n"


def _factypes_row(nfaclyrs=3, wallid=1, lam=0.4):
    """One valid factypes data row: 6 header + d(k) + C(k) + l(k) + k(k+1)."""
    k = nfaclyrs
    fields = [str(wallid), "0", "0.1", "0.01", "0.2", "0.9"]
    fields += ["0.1"] * k          # d: layer thicknesses
    fields += ["1000"] * k         # C: heat capacity
    fields += ["1.0"] * k          # l: conductivity
    fields += [str(lam)] * (k + 1)  # k: conductivity (k+1 columns) -> factypes['lam']
    return " ".join(fields)


def _write_nc(path, data_vars, coords):
    ds = xr.Dataset(data_vars, coords=coords)
    ds.to_netcdf(path)
    ds.close()


class _CaseBase(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)
        (self.workdir / "namoptions.001").write_text(_NAMOPTIONS, encoding="ascii")

    def _sim(self):
        return UDBase("1", self.workdir, load_geometry=False, suppress_load_warnings=True)

    def _write_facet_files(self, n_facets=2, lam=0.4):
        """Write facetarea/facets/factypes so facet loaders + assign_prop work."""
        areas = "\n".join("1.0" for _ in range(n_facets))
        (self.workdir / "facetarea.inp.001").write_text("# area\n" + areas + "\n", encoding="ascii")
        rows = "\n".join("1 0.0 0.0 1.0" for _ in range(n_facets))
        (self.workdir / "facets.inp.001").write_text("# typeid nx ny nz\n" + rows + "\n", encoding="ascii")
        (self.workdir / "factypes.inp.001").write_text(
            "# h1\n# h2\n# h3\n" + _factypes_row(3, wallid=1, lam=lam) + "\n",
            encoding="ascii",
        )


class TestNetCDFLoaders(_CaseBase):
    """load_field / load_stat_* / load_slice / load_fac_* over _load_ncdata.

    Pins the MATLAB-style dimension reversal (C-order (time, z, y, x) becomes
    (x, y, z, time)) and the numpy-vs-Dataset return contract.
    """

    def test_load_field_returns_transposed_array(self):
        # NetCDF (time, zt, yt, xm) = (2, 2, 3, 4) -> reversed (xm, yt, zt, time).
        data = np.arange(2 * 2 * 3 * 4, dtype=float).reshape(2, 2, 3, 4)
        _write_nc(
            self.workdir / "fielddump.001.nc",
            {"u": (("time", "zt", "yt", "xm"), data)},
            {"time": [0, 1], "zt": [0, 1], "yt": [0, 1, 2], "xm": [0, 1, 2, 3]},
        )
        u = self._sim().load_field("u")
        self.assertIsInstance(u, np.ndarray)
        self.assertEqual(u.shape, (4, 3, 2, 2))
        # value at reversed index equals original at the mirrored index
        self.assertEqual(u[3, 2, 1, 0], data[0, 1, 2, 3])

    def test_load_field_browse_returns_dataset(self):
        data = np.zeros((2, 2, 3, 4))
        _write_nc(
            self.workdir / "fielddump.001.nc",
            {"u": (("time", "zt", "yt", "xm"), data)},
            {"time": [0, 1], "zt": [0, 1], "yt": [0, 1, 2], "xm": [0, 1, 2, 3]},
        )
        with redirect_stdout(io.StringIO()):
            ds = self._sim().load_field()  # var=None -> browse
        self.assertIsInstance(ds, xr.Dataset)
        self.assertIn("u", ds)
        self.assertEqual(ds["u"].shape, (4, 3, 2, 2))  # already transposed

    def test_load_field_missing_file_raises(self):
        with self.assertRaises(FileNotFoundError):
            self._sim().load_field("u")

    def test_load_stat_xyt_returns_transposed_array(self):
        # (time, zt) = (3, 2) -> (2, 3)
        data = np.arange(6.0).reshape(3, 2)
        _write_nc(
            self.workdir / "xytdump.001.nc",
            {"uxyt": (("time", "zt"), data)},
            {"time": [0, 1, 2], "zt": [0, 1]},
        )
        arr = self._sim().load_stat_xyt("uxyt")
        self.assertEqual(arr.shape, (2, 3))
        np.testing.assert_allclose(arr, data.T)

    def test_load_stat_t_returns_transposed_array(self):
        data = np.arange(2 * 2 * 3 * 4, dtype=float).reshape(2, 2, 3, 4)
        _write_nc(
            self.workdir / "tdump.001.nc",
            {"u": (("time", "zt", "yt", "xt"), data)},
            {"time": [0, 1], "zt": [0, 1], "yt": [0, 1, 2], "xt": [0, 1, 2, 3]},
        )
        arr = self._sim().load_stat_t("u")
        self.assertEqual(arr.shape, (4, 3, 2, 2))

    def test_load_stat_tree_returns_transposed_array(self):
        data = np.arange(6.0).reshape(3, 2)
        _write_nc(
            self.workdir / "treedump.001.nc",
            {"tree_drag_u": (("time", "zt"), data)},
            {"time": [0, 1, 2], "zt": [0, 1]},
        )
        arr = self._sim().load_stat_tree("tree_drag_u")
        self.assertEqual(arr.shape, (2, 3))

    def test_load_slice_k_returns_transposed_array(self):
        # horizontal slice: (time, yt, xm) = (2, 3, 4) -> (4, 3, 2)
        data = np.arange(2 * 3 * 4, dtype=float).reshape(2, 3, 4)
        _write_nc(
            self.workdir / "kslicedump.001.nc",
            {"u": (("time", "yt", "xm"), data)},
            {"time": [0, 1], "yt": [0, 1, 2], "xm": [0, 1, 2, 3]},
        )
        arr = self._sim().load_slice("k", "u")
        self.assertEqual(arr.shape, (4, 3, 2))

    def test_load_slice_invalid_plane_raises(self):
        with self.assertRaises(ValueError):
            self._sim().load_slice("z", "u")

    def test_load_fac_momentum_returns_transposed_array(self):
        # facet data (time, fct) = (2, 3) -> (fct, time) = (3, 2)
        data = np.arange(6.0).reshape(2, 3)
        _write_nc(
            self.workdir / "fac.001.nc",
            {"pres": (("time", "fct"), data)},
            {"time": [0, 1], "fct": [0, 1, 2]},
        )
        arr = self._sim().load_fac_momentum("pres")
        self.assertEqual(arr.shape, (3, 2))
        np.testing.assert_allclose(arr, data.T)

    def test_load_fac_eb_returns_transposed_array(self):
        data = np.arange(6.0).reshape(2, 3)
        _write_nc(
            self.workdir / "facEB.001.nc",
            {"hf": (("time", "fct"), data)},
            {"time": [0, 1], "fct": [0, 1, 2]},
        )
        arr = self._sim().load_fac_eb("hf")
        self.assertEqual(arr.shape, (3, 2))

    def test_load_fac_temperature_returns_transposed_array(self):
        # (time, lyr, fct) = (2, 3, 4) -> (fct, lyr, time) = (4, 3, 2)
        data = np.arange(2 * 3 * 4, dtype=float).reshape(2, 3, 4)
        _write_nc(
            self.workdir / "facT.001.nc",
            {"T": (("time", "lyr", "fct"), data)},
            {"time": [0, 1], "lyr": [0, 1, 2], "fct": [0, 1, 2, 3]},
        )
        arr = self._sim().load_fac_temperature("T")
        self.assertEqual(arr.shape, (4, 3, 2))


class TestLoadSeb(_CaseBase):
    """load_seb assembles facEB + facT + material conductivity into the SEB dict."""

    def _write_seb_inputs(self, n_facets=2, n_time=2, lam=0.4):
        fct = list(range(n_facets))
        time = list(range(n_time))
        lyr = [0, 1, 2]
        ones_ft = np.ones((n_time, n_facets))
        _write_nc(
            self.workdir / "facEB.001.nc",
            {
                "t": (("time",), np.asarray(time, dtype=float)),
                "netsw": (("time", "fct"), 5.0 * ones_ft),
                "LWin": (("time", "fct"), 4.0 * ones_ft),
                "LWout": (("time", "fct"), 1.0 * ones_ft),
                "hf": (("time", "fct"), 2.0 * ones_ft),
                "ef": (("time", "fct"), 3.0 * ones_ft),
            },
            {"time": time, "fct": fct},
        )
        T = 300.0 * np.ones((n_time, len(lyr), n_facets))
        dTdz = np.ones((n_time, len(lyr), n_facets))  # layer-0 gradient == 1
        _write_nc(
            self.workdir / "facT.001.nc",
            {"T": (("time", "lyr", "fct"), T), "dTdz": (("time", "lyr", "fct"), dTdz)},
            {"time": time, "lyr": lyr, "fct": fct},
        )
        self._write_facet_files(n_facets=n_facets, lam=lam)

    def test_load_seb_returns_expected_terms_and_shapes(self):
        self._write_seb_inputs(n_facets=2, n_time=2, lam=0.4)
        seb = self._sim().load_seb()

        self.assertEqual(
            set(seb), {"Kstar", "Lstar", "Lin", "Lout", "H", "E", "G", "Tsurf", "t"}
        )
        # facet-major (n_facets, n_time) = (2, 2); t is 1-D (n_time,)
        for key in ("Kstar", "Lstar", "Lin", "Lout", "H", "E", "G", "Tsurf"):
            self.assertEqual(seb[key].shape, (2, 2), msg=key)
        self.assertEqual(seb["t"].shape, (2,))

        np.testing.assert_allclose(seb["Kstar"], 5.0)
        np.testing.assert_allclose(seb["Lstar"], 4.0 - 1.0)  # Lin - Lout
        np.testing.assert_allclose(seb["H"], -2.0)           # sign convention: -hf
        np.testing.assert_allclose(seb["E"], -3.0)           # sign convention: -ef
        # G = -lam[:,0] * dTdz[:,0]; factypes['lam'] (conductivity) is the 'l'
        # block of the factypes row, i.e. 1.0 here, so G = -1.0 * 1 = -1.0.
        np.testing.assert_allclose(seb["G"], -1.0)
        np.testing.assert_allclose(seb["Tsurf"], 300.0)

    def test_load_seb_missing_file_raises(self):
        # No facEB file present.
        with self.assertRaises(FileNotFoundError):
            self._sim().load_seb()


class TestAreaAverage(_CaseBase):
    """Area-weighted facet averaging, hand-checkable on a 2-facet case."""

    def _sim_with_areas(self, areas):
        txt = "# area\n" + "\n".join(str(a) for a in areas) + "\n"
        (self.workdir / "facetarea.inp.001").write_text(txt, encoding="ascii")
        return self._sim()

    def test_area_average_fac_1d_hand_value(self):
        sim = self._sim_with_areas([1.0, 3.0])
        # (10*1 + 20*3) / (1 + 3) = 70/4 = 17.5
        self.assertAlmostEqual(float(sim.area_average_fac(np.array([10.0, 20.0]))), 17.5)

    def test_area_average_fac_2d_facets_first(self):
        sim = self._sim_with_areas([1.0, 3.0])
        var = np.array([[10.0, 10.0, 10.0], [20.0, 20.0, 20.0]])  # (n_facets, n_time)
        out = sim.area_average_fac(var)
        self.assertEqual(out.shape, (3,))
        np.testing.assert_allclose(out, 17.5)

    def test_area_average_fac_2d_facets_second(self):
        sim = self._sim_with_areas([1.0, 3.0])
        var = np.array([[10.0, 20.0], [10.0, 20.0], [10.0, 20.0]])  # (n_time, n_facets)
        out = sim.area_average_fac(var)
        self.assertEqual(out.shape, (3,))
        np.testing.assert_allclose(out, 17.5)

    def test_area_average_fac_with_selection(self):
        sim = self._sim_with_areas([1.0, 3.0])
        sel = np.array([False, True])
        # only facet 1 selected: (20*3)/3 = 20
        self.assertAlmostEqual(float(sim.area_average_fac(np.array([10.0, 20.0]), sel)), 20.0)

    def test_area_average_fac_requires_facetarea(self):
        sim = self._sim()  # no facetarea.inp -> _lffacetarea False
        with self.assertRaises(ValueError):
            sim.area_average_fac(np.array([1.0, 2.0]))

    def test_area_average_seb_hand_values(self):
        sim = self._sim_with_areas([1.0, 3.0])
        base = np.array([[10.0, 10.0, 10.0], [20.0, 20.0, 20.0]])  # -> avg 17.5
        seb = {
            "Kstar": base,
            "Lstar": 2.0 * np.ones((2, 3)),
            "Lin": 3.0 * np.ones((2, 3)),
            "Lout": 4.0 * np.ones((2, 3)),
            "H": 5.0 * np.ones((2, 3)),
            "E": 6.0 * np.ones((2, 3)),
            "G": 7.0 * np.ones((2, 3)),
            "t": np.array([0.0, 1.0, 2.0]),
        }
        out = sim.area_average_seb(seb)
        self.assertEqual(set(out), {"Kstar", "Lstar", "Lin", "Lout", "H", "E", "G", "t"})
        np.testing.assert_allclose(out["Kstar"], 17.5)
        np.testing.assert_allclose(out["Lstar"], 2.0)
        np.testing.assert_allclose(out["Lin"], 3.0)
        np.testing.assert_allclose(out["Lout"], 4.0)
        np.testing.assert_allclose(out["H"], 5.0)
        np.testing.assert_allclose(out["E"], 6.0)
        np.testing.assert_allclose(out["G"], 7.0)
        np.testing.assert_allclose(out["t"], [0.0, 1.0, 2.0])  # passed through unchanged


class TestConvertFacetToField(_CaseBase):
    """convert_facflx_to_field (density) and convert_facvar_to_field (weighted mean)."""

    def test_convert_facflx_to_field_hand_values(self):
        sim = self._sim()
        facsec = {
            "facid": np.array([0, 1]),
            "area": np.array([2.0, 4.0]),
            "locs": np.array([[0, 0, 0], [1, 1, 1]]),
        }
        var = np.array([10.0, 20.0])
        fld = sim.convert_facflx_to_field(var, facsec, sim.dzt)
        self.assertEqual(fld.shape, (4, 3, 2))
        self.assertEqual(fld.dtype, np.float32)
        # cell_volume = dx*dy*dz = 10*10*10 = 1000
        self.assertAlmostEqual(fld[0, 0, 0], 10.0 * 2.0 / 1000.0, places=6)  # 0.02
        self.assertAlmostEqual(fld[1, 1, 1], 20.0 * 4.0 / 1000.0, places=6)  # 0.08
        # everything else zero
        mask = np.ones_like(fld, dtype=bool)
        mask[0, 0, 0] = mask[1, 1, 1] = False
        self.assertTrue(np.all(fld[mask] == 0.0))

    def test_convert_facvar_to_field_recovers_facet_value(self):
        # Normalised transfer: a single facet per cell must recover its own value.
        sim = self._sim()
        facsec = {
            "facid": np.array([0, 1]),
            "area": np.array([2.0, 4.0]),
            "locs": np.array([[0, 0, 0], [1, 1, 1]]),
        }
        var = np.array([10.0, 20.0])
        fld = sim.convert_facvar_to_field(var, facsec)
        self.assertEqual(fld.shape, (4, 3, 2))
        self.assertAlmostEqual(fld[0, 0, 0], 10.0, places=4)
        self.assertAlmostEqual(fld[1, 1, 1], 20.0, places=4)
        # empty cells normalised to 0 (norm forced to 1 where flux is 0)
        self.assertAlmostEqual(fld[2, 2, 0], 0.0, places=6)


class TestFrontalProperties(_CaseBase):
    """calculate_frontal_properties over a single -x-facing facet."""

    def test_requires_geometry(self):
        sim = self._sim()  # load_geometry=False -> sim.geom is None
        with self.assertRaises(ValueError):
            sim.calculate_frontal_properties()

    def test_hand_values_single_frontal_facet(self):
        sim = self._sim()
        # Minimal stand-in geometry: one face whose outward normal points in -x
        # (so it is frontal to the x-direction flow). Only face_normals is read.
        sim.geom = SimpleNamespace(face_normals=np.array([[-1.0, 0.0, 0.0]]))
        sim.facsec["c"] = {
            "facid": np.array([0]),
            "area": np.array([100.0]),
            "locs": np.array([[0, 0, 0]]),
        }

        res = sim.calculate_frontal_properties()

        self.assertEqual(set(res), {"skylinex", "skyliney", "Afx", "Afy", "brx", "bry"})
        self.assertEqual(res["skylinex"].shape, (3, 2))  # (jtot, ktot)
        self.assertEqual(res["skyliney"].shape, (4, 2))  # (itot, ktot)
        # frontal area = phix * area = 1 * 100
        self.assertAlmostEqual(res["Afx"], 100.0, places=3)
        self.assertAlmostEqual(res["Afy"], 0.0, places=6)
        # blockage = (dy*dz of the one blocked column) / (ylen*zsize) = 100/600
        self.assertAlmostEqual(res["brx"], 100.0 / 600.0, places=5)
        self.assertAlmostEqual(res["bry"], 0.0, places=6)
        # only cell (j=0, k=0) is blocked in the x-skyline
        self.assertEqual(res["skylinex"][0, 0], 1.0)
        self.assertEqual(res["skylinex"].sum(), 1.0)
        self.assertEqual(res["skyliney"].sum(), 0.0)


if __name__ == "__main__":
    unittest.main()
