"""Characterization tests for udprep.udprep_seb (SEBSection).

These pin the CURRENT behaviour of the surface-energy-balance preprocessing
section as a regression net ahead of the planned SEB redesign. They must PASS
against the code as-is; where current behaviour looks suspect it is asserted
verbatim and flagged in the review report rather than "fixed".
"""
import sys
import types
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest import mock

import numpy as np

from _common import PYTHON_DIR  # noqa: F401  (kept for parity with sibling tests)

from udprep.udprep_seb import SEBSection, SPEC, DEFAULTS, FIELDS  # noqa: E402
from exceptions import DependencyError  # noqa: E402


def make_seb_section(workdir, *, expnr="777", nfcts=3, lEB=False,
                     iwallmom=1, iwalltemp=1, iwallmoist=1, **field_overrides):
    """Build a SEBSection over a minimal stub sim.

    The SEB fields come from DEFAULTS (overridable); the run-condition gates
    (lEB / iwallmom / iwalltemp / iwallmoist) and nfcts live on the stub sim,
    mirroring how other sections populate them.
    """
    sim = SimpleNamespace(
        path=Path(workdir),
        expnr=expnr,
        nfcts=nfcts,
        lEB=lEB,
        iwallmom=iwallmom,
        iwalltemp=iwalltemp,
        iwallmoist=iwallmoist,
    )
    values = {
        "facT": field_overrides.pop("facT", 288.0),
        "lfacTlyrs": field_overrides.pop("lfacTlyrs", False),
        "nfaclyrs": field_overrides.pop("nfaclyrs", 10),
        "facT_file": field_overrides.pop("facT_file", ""),
        "dtEB": field_overrides.pop("dtEB", 10.0),
    }
    values.update(field_overrides)
    section = SEBSection("seb", values, sim=sim, defaults=dict(values))
    return section, sim


def fake_netcdf_module(variables):
    """A stand-in netCDF4 module whose Dataset yields the given variables."""
    module = types.ModuleType("netCDF4")

    class _Dataset:
        def __init__(self, path, mode):
            self.variables = variables

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    module.Dataset = _Dataset
    return module


class TestSEBSpecAndDefaults(unittest.TestCase):
    def test_spec_metadata_matches_section(self):
        self.assertEqual(SPEC.name, "seb")
        self.assertIs(SPEC.section_cls, SEBSection)
        self.assertEqual(set(SPEC.fields), set(FIELDS))

    def test_defaults_expose_expected_fields(self):
        self.assertEqual(
            set(FIELDS),
            {"facT", "lfacTlyrs", "nfaclyrs", "facT_file", "dtEB"},
        )
        self.assertEqual(DEFAULTS["facT"], 288.0)
        self.assertEqual(DEFAULTS["lfacTlyrs"], False)
        self.assertEqual(DEFAULTS["nfaclyrs"], 10)
        self.assertEqual(DEFAULTS["facT_file"], "")
        self.assertEqual(DEFAULTS["dtEB"], 10.0)

    def test_construction_proxies_field_values_to_sim(self):
        with TemporaryDirectory() as tmp:
            section, sim = make_seb_section(tmp, facT=295.0, lfacTlyrs=True)
            self.assertEqual(section.facT, 295.0)
            self.assertEqual(sim.facT, 295.0)
            self.assertTrue(section.lfacTlyrs)
            self.assertEqual(section.nfaclyrs, 10)


class TestWriteTfacinit(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)

    def test_writes_uniform_facet_temperature_file(self):
        section, _ = make_seb_section(self.workdir, expnr="654", nfcts=4, facT=290.5)
        section.write_Tfacinit()

        path = self.workdir / "Tfacinit.inp.654"
        self.assertTrue(path.exists())
        lines = path.read_text(encoding="ascii").splitlines()
        self.assertTrue(lines[0].startswith("#"))
        # header + one line per facet
        self.assertEqual(len(lines), 1 + 4)

        values = np.loadtxt(path)
        self.assertEqual(values.shape, (4,))
        np.testing.assert_allclose(values, np.full(4, 290.5))

    def test_nfcts_zero_raises(self):
        section, _ = make_seb_section(self.workdir, nfcts=0)
        with self.assertRaisesRegex(ValueError, "nfcts"):
            section.write_Tfacinit()

    def test_sim_none_raises(self):
        section = SEBSection("seb", dict(DEFAULTS), sim=None, defaults=dict(DEFAULTS))
        with self.assertRaisesRegex(ValueError, "UDBase"):
            section.write_Tfacinit()


class TestRunAllGating(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)

    def _tfacinit(self, expnr="777"):
        return self.workdir / f"Tfacinit.inp.{expnr}"

    def test_skips_when_no_facet_temperature_needed(self):
        # lEB False and no iwall*==2 -> section is a no-op.
        section, _ = make_seb_section(self.workdir)
        section.run_all()
        self.assertFalse(self._tfacinit().exists())

    def test_runs_when_lEB_enabled(self):
        section, _ = make_seb_section(self.workdir, lEB=True)
        section.run_all()
        self.assertTrue(self._tfacinit().exists())

    def test_runs_when_iwallmom_is_2(self):
        # P7: SEBSection.run_all also gates on iwallmom/iwalltemp/iwallmoist==2,
        # unlike UDPrep.run_all which only gates SEB on (libm and lEB).
        section, _ = make_seb_section(self.workdir, lEB=False, iwallmom=2)
        section.run_all()
        self.assertTrue(self._tfacinit().exists())

    def test_runs_when_iwalltemp_is_2(self):
        section, _ = make_seb_section(self.workdir, lEB=False, iwalltemp=2)
        section.run_all()
        self.assertTrue(self._tfacinit().exists())

    def test_runs_when_iwallmoist_is_2(self):
        section, _ = make_seb_section(self.workdir, lEB=False, iwallmoist=2)
        section.run_all()
        self.assertTrue(self._tfacinit().exists())

    def test_lfacTlyrs_dispatches_to_layers_writer(self):
        section, _ = make_seb_section(self.workdir, lEB=True, lfacTlyrs=True)
        section.write_Tfacinit = mock.Mock()
        section.write_Tfacinit_layers = mock.Mock()
        section.run_all()
        section.write_Tfacinit_layers.assert_called_once_with()
        section.write_Tfacinit.assert_not_called()

    def test_default_dispatches_to_scalar_writer(self):
        section, _ = make_seb_section(self.workdir, lEB=True, lfacTlyrs=False)
        section.write_Tfacinit = mock.Mock()
        section.write_Tfacinit_layers = mock.Mock()
        section.run_all()
        section.write_Tfacinit.assert_called_once_with()
        section.write_Tfacinit_layers.assert_not_called()

    def test_run_all_sim_none_raises(self):
        section = SEBSection("seb", dict(DEFAULTS), sim=None, defaults=dict(DEFAULTS))
        with self.assertRaisesRegex(ValueError, "UDBase"):
            section.run_all()


class TestWriteTfacinitLayers(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)

    def test_writes_last_layer_temperatures(self):
        facT_path = self.workdir / "facT.nc"
        facT_path.write_bytes(b"placeholder")  # existence check only; Dataset mocked
        section, _ = make_seb_section(
            self.workdir, expnr="654", lfacTlyrs=True, facT_file="facT.nc"
        )
        # T shape (nfcts=2, ntimes=3, nlayers=4); writer takes the last layer.
        Tfac = np.arange(24, dtype=float).reshape(2, 3, 4)
        fake = fake_netcdf_module({"T": Tfac})
        with mock.patch.dict(sys.modules, {"netCDF4": fake}):
            section.write_Tfacinit_layers()

        out = self.workdir / "Tfacinit_layers.inp.654"
        self.assertTrue(out.exists())
        values = np.loadtxt(out)
        self.assertEqual(values.shape, (2, 3))
        np.testing.assert_allclose(values, Tfac[:, :, -1])

    def test_missing_facT_file_field_raises(self):
        section, _ = make_seb_section(self.workdir, lfacTlyrs=True, facT_file="")
        with self.assertRaisesRegex(ValueError, "facT_file"):
            section.write_Tfacinit_layers()

    def test_nonexistent_facT_file_raises(self):
        section, _ = make_seb_section(
            self.workdir, lfacTlyrs=True, facT_file="does_not_exist.nc"
        )
        with self.assertRaises(FileNotFoundError):
            section.write_Tfacinit_layers()

    def test_missing_T_variable_raises(self):
        facT_path = self.workdir / "facT.nc"
        facT_path.write_bytes(b"placeholder")
        section, _ = make_seb_section(
            self.workdir, lfacTlyrs=True, facT_file="facT.nc"
        )
        fake = fake_netcdf_module({"S": np.zeros((2, 3, 4))})
        with mock.patch.dict(sys.modules, {"netCDF4": fake}):
            with self.assertRaisesRegex(ValueError, "'T'"):
                section.write_Tfacinit_layers()

    def test_T_variable_not_3d_raises(self):
        facT_path = self.workdir / "facT.nc"
        facT_path.write_bytes(b"placeholder")
        section, _ = make_seb_section(
            self.workdir, lfacTlyrs=True, facT_file="facT.nc"
        )
        fake = fake_netcdf_module({"T": np.zeros((2, 3))})
        with mock.patch.dict(sys.modules, {"netCDF4": fake}):
            with self.assertRaisesRegex(ValueError, "3D"):
                section.write_Tfacinit_layers()

    def test_missing_netcdf4_dependency_raises_dependency_error(self):
        facT_path = self.workdir / "facT.nc"
        facT_path.write_bytes(b"placeholder")
        section, _ = make_seb_section(
            self.workdir, lfacTlyrs=True, facT_file="facT.nc"
        )
        # Setting the module to None in sys.modules makes `import netCDF4` raise
        # ImportError, which the writer wraps as DependencyError.
        with mock.patch.dict(sys.modules, {"netCDF4": None}):
            with self.assertRaises(DependencyError):
                section.write_Tfacinit_layers()

    def test_sim_none_raises(self):
        section = SEBSection("seb", dict(DEFAULTS), sim=None, defaults=dict(DEFAULTS))
        with self.assertRaisesRegex(ValueError, "UDBase"):
            section.write_Tfacinit_layers()


if __name__ == "__main__":
    unittest.main()
