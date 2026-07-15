import sys
import unittest
import warnings as _warnings
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from _common import PYTHON_DIR

from udprep.udprep_forcing import ForcingSection  # noqa: E402

class DummySim:
    def __init__(self, path: Path, expnr="100"):
        self.path = path
        self.expnr = expnr
        # Non-uniform grid: faces at [0, 1, 3, 6] → centres at [0.5, 2.0, 4.5], widths [1, 2, 3]
        self.zt = np.array([0.5, 2.0, 4.5])
        self.dzt = np.array([1.0, 2.0, 3.0])
        self.ktot = 3
        self.zsize = 6.0
        self.idriver = 0
        self.ltimedepnudge = False
        self.profsourcefile = ""

    def read_matrix(self, path, skiprows=0):
        return np.loadtxt(path, skiprows=skiprows)

    def load_prof(self):
        fpath = Path(self.path) / f"prof.inp.{self.expnr}"
        if not fpath.exists():
            raise FileNotFoundError(f"Profile file not found: {fpath}")
        return np.loadtxt(fpath, skiprows=2)

class TestForcingSection(unittest.TestCase):
    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.workdir = Path(self.temp_dir.name)
        self.sim = DummySim(self.workdir, expnr="321")

    # ---- generate_prof ---------------------------------------------------------

    def test_generate_prof_without_lapse(self):
        section = ForcingSection(
            "forcing",
            {"thl0": 290.0, "qt0": 0.01, "u0": 4.0, "v0": 1.0, "tke": 0.5, "lapse": 0.0},
            sim=self.sim,
            defaults={},
        )
        section.generate_prof()
        np.testing.assert_allclose(section.pr[:, 0], self.sim.zt)
        np.testing.assert_allclose(section.pr[:, 1], [290.0, 290.0, 290.0])
        np.testing.assert_allclose(section.pr[:, 2:], [[0.01, 4.0, 1.0, 0.5]] * 3)

    def test_generate_prof_with_lapse(self):
        section = ForcingSection(
            "forcing",
            {"thl0": 300.0, "qt0": 0.0, "u0": 2.0, "v0": 0.0, "tke": 0.1, "lapse": 0.5},
            sim=self.sim,
            defaults={},
        )
        section.generate_prof()
        # With non-uniform dzt=[1,2,3] and lapse=0.5:
        #   thl[1] = 300.0 + 0.5*(1+2)/2 = 300.75
        #   thl[2] = 300.75 + 0.5*(2+3)/2 = 302.0
        np.testing.assert_allclose(section.pr[:, 1], [300.0, 300.75, 302.0])

    # ---- write_prof ------------------------------------------------------------

    def test_write_prof_creates_expected_file(self):
        section = ForcingSection(
            "forcing",
            {"thl0": 290.0, "qt0": 0.01, "u0": 4.0, "v0": 1.0, "tke": 0.5, "lapse": 0.0},
            sim=self.sim,
            defaults={},
        )
        section.generate_prof()
        section.write_prof()
        content = (self.workdir / "prof.inp.321").read_text(encoding="ascii").splitlines()
        self.assertEqual(content[0], "# SDBL flow ")
        self.assertTrue(content[1].startswith("# z thl qt u v tke"))
        self.assertEqual(len(content), 5)

    def test_write_prof_force_overwrites_existing_file(self):
        section = ForcingSection(
            "forcing",
            {"thl0": 290.0, "qt0": 0.0, "u0": 0.0, "v0": 0.0, "tke": 0.0, "lapse": 0.0},
            sim=self.sim,
            defaults={},
        )
        section.generate_prof()
        section.write_prof()                         # creates file
        with self.assertWarns(UserWarning):          # without force: must warn
            section.write_prof()
        with _warnings.catch_warnings():             # with force: must NOT warn
            _warnings.simplefilter("error")
            section.write_prof(force=True)

    # ---- generate_lscale -------------------------------------------------------

    def _lscale_base(self, **overrides):
        """Return a minimal namoptions dict for generate_lscale with no active forcing."""
        base = {
            "w_s": 0.0, "dqtdxls": 0.0, "dqtdyls": 0.0, "dqtdtls": 0.0, "R": 0.0,
            "luoutflowr": 0, "lvoutflowr": 0, "luvolflowr": 0, "lvvolflowr": 0,
            "lprofforc": 0, "lcoriol": 0, "lnudge": 0,
        }
        base.update(overrides)
        return base

    def test_generate_lscale_ldp_sets_pressure_gradient_and_warns(self):
        section = ForcingSection(
            "forcing",
            {**self._lscale_base(), "dpdx": 0.001, "dpdy": 0.002},
            sim=self.sim,
            defaults={},
        )
        with self.assertWarns(UserWarning):
            section.generate_lscale()
        np.testing.assert_allclose(section.ls[:, 3], 0.001)
        np.testing.assert_allclose(section.ls[:, 4], 0.002)
        np.testing.assert_allclose(section.ls[:, 1:3], 0.0)  # geostrophic columns untouched

    def test_generate_lscale_profforc_sets_geostrophic_wind(self):
        section = ForcingSection(
            "forcing",
            {**self._lscale_base(lprofforc=1), "u0": 5.0, "v0": 2.0},
            sim=self.sim,
            defaults={},
        )
        section.generate_lscale()
        np.testing.assert_allclose(section.ls[:, 1], 5.0)
        np.testing.assert_allclose(section.ls[:, 2], 2.0)
        np.testing.assert_allclose(section.ls[:, 3:5], 0.0)  # pressure gradient untouched

    def test_generate_lscale_coriol_sets_geostrophic_wind(self):
        section = ForcingSection(
            "forcing",
            {**self._lscale_base(lcoriol=1), "u0": 3.0, "v0": 0.5},
            sim=self.sim,
            defaults={},
        )
        section.generate_lscale()
        np.testing.assert_allclose(section.ls[:, 1], 3.0)
        np.testing.assert_allclose(section.ls[:, 2], 0.5)
        np.testing.assert_allclose(section.ls[:, 3:5], 0.0)

    def test_generate_lscale_idriver2_no_forcing_not_ldp(self):
        section = ForcingSection(
            "forcing",
            {**self._lscale_base(), "idriver": 2},
            sim=self.sim,
            defaults={},
        )
        with _warnings.catch_warnings():
            _warnings.simplefilter("error")
            section.generate_lscale()          # must not warn about missing forcing
        np.testing.assert_allclose(section.ls[:, 1:5], 0.0)

    def test_generate_lscale_conflicting_forcing_raises(self):
        section = ForcingSection(
            "forcing",
            self._lscale_base(luoutflowr=1, lprofforc=1),
            sim=self.sim,
            defaults={},
        )
        with self.assertRaises(ValueError):
            section.generate_lscale()

    def test_generate_lscale_populates_vertical_profile_columns(self):
        section = ForcingSection(
            "forcing",
            {**self._lscale_base(), "w_s": 0.01, "dqtdxls": 1e-7, "dqtdyls": 2e-7,
             "dqtdtls": 3e-7, "R": 4e-7, "dpdx": 0.0, "dpdy": 0.0},
            sim=self.sim,
            defaults={},
        )
        with self.assertWarns(UserWarning):
            section.generate_lscale()
        np.testing.assert_allclose(section.ls[:, 0], self.sim.zt)
        np.testing.assert_allclose(section.ls[:, 5], 0.01)
        np.testing.assert_allclose(section.ls[:, 6], 1e-7)
        np.testing.assert_allclose(section.ls[:, 7], 2e-7)
        np.testing.assert_allclose(section.ls[:, 8], 3e-7)
        np.testing.assert_allclose(section.ls[:, 9], 4e-7)

    # ---- write_lscale ----------------------------------------------------------

    def test_write_lscale_creates_expected_file(self):
        section = ForcingSection(
            "forcing",
            {**self._lscale_base(), "dpdx": 0.0, "dpdy": 0.0},
            sim=self.sim,
            defaults={},
        )
        with self.assertWarns(UserWarning):
            section.generate_lscale()
        section.write_lscale()
        content = (self.workdir / "lscale.inp.321").read_text(encoding="ascii").splitlines()
        self.assertEqual(content[0], "# SDBL flow ")
        self.assertTrue(content[1].startswith("# z uq vq pqx pqy"))
        self.assertEqual(len(content), 5)  # 2 header lines + 3 data rows (ktot=3)

    def test_write_lscale_auto_generates_if_not_generated(self):
        section = ForcingSection(
            "forcing",
            {**self._lscale_base(), "dpdx": 0.0, "dpdy": 0.0},
            sim=self.sim,
            defaults={},
        )
        with self.assertWarns(UserWarning):   # ldp warning from auto-generated lscale
            section.write_lscale()
        self.assertTrue((self.workdir / "lscale.inp.321").exists())

    # ---- update_prof_from_nudge_data -------------------------------------------

    def test_update_prof_from_nudge_data_missing_profinp_raises(self):
        section = ForcingSection(
            "forcing", {}, sim=self.sim, defaults={},
        )
        with self.assertRaises(FileNotFoundError):
            section.update_prof_from_nudge_data("/nonexistent/source.dat")

    def test_update_prof_from_nudge_data_missing_sourcefile_warns_and_returns(self):
        section = ForcingSection(
            "forcing",
            {"thl0": 295.0, "qt0": 0.005, "u0": 3.0, "v0": 0.0, "tke": 0.2, "lapse": 0.0},
            sim=self.sim,
            defaults={},
        )
        section.generate_prof()
        section.write_prof()
        with self.assertWarns(UserWarning):
            section.update_prof_from_nudge_data(str(self.workdir / "missing_prdata.dat"))
        # prof.inp must remain unchanged — thl column still 295.0
        data = np.loadtxt(self.workdir / "prof.inp.321", skiprows=2)
        np.testing.assert_allclose(data[:, 1], 295.0, rtol=1e-5)

    def test_update_prof_from_nudge_data_invalid_sourcefile_warns_and_returns(self):
        section = ForcingSection(
            "forcing",
            {"thl0": 295.0, "qt0": 0.005, "u0": 3.0, "v0": 0.0, "tke": 0.2, "lapse": 0.0},
            sim=self.sim,
            defaults={},
        )
        section.generate_prof()
        section.write_prof()
        original = (self.workdir / "prof.inp.321").read_text(encoding="ascii")
        malformed = self.workdir / "malformed_prdata.dat"
        malformed.write_text("z u v\n1.0 2.0 3.0\n", encoding="ascii")

        with _warnings.catch_warnings(record=True) as caught:
            _warnings.simplefilter("always")
            section.update_prof_from_nudge_data(str(malformed))

        self.assertTrue(any("original prof.inp is kept" in str(w.message) for w in caught))
        self.assertEqual((self.workdir / "prof.inp.321").read_text(encoding="ascii"), original)

    def test_update_prof_from_nudge_data_empty_sourcefile_warns_and_returns(self):
        section = ForcingSection(
            "forcing",
            {"thl0": 295.0, "qt0": 0.005, "u0": 3.0, "v0": 0.0, "tke": 0.2, "lapse": 0.0},
            sim=self.sim,
            defaults={},
        )
        section.generate_prof()
        section.write_prof()
        original = (self.workdir / "prof.inp.321").read_text(encoding="ascii")
        empty = self.workdir / "empty_prdata.dat"
        empty.write_text("", encoding="ascii")

        with _warnings.catch_warnings(record=True) as caught:
            _warnings.simplefilter("always")
            section.update_prof_from_nudge_data(str(empty))

        self.assertTrue(any("original prof.inp is kept" in str(w.message) for w in caught))
        self.assertEqual((self.workdir / "prof.inp.321").read_text(encoding="ascii"), original)

    def test_update_prof_from_nudge_data_blank_sourcefile_warns_and_returns(self):
        section = ForcingSection(
            "forcing",
            {"thl0": 295.0, "qt0": 0.005, "u0": 3.0, "v0": 0.0, "tke": 0.2, "lapse": 0.0},
            sim=self.sim,
            defaults={},
        )
        section.generate_prof()
        section.write_prof()
        original = (self.workdir / "prof.inp.321").read_text(encoding="ascii")

        with _warnings.catch_warnings(record=True) as caught:
            _warnings.simplefilter("always")
            section.update_prof_from_nudge_data("")

        self.assertTrue(any("original prof.inp is kept" in str(w.message) for w in caught))
        self.assertEqual((self.workdir / "prof.inp.321").read_text(encoding="ascii"), original)

    def test_update_prof_from_nudge_data_applies_interpolation(self):
        from scipy.interpolate import CubicSpline
        section = ForcingSection(
            "forcing",
            {"thl0": 290.0, "qt0": 0.0, "u0": 0.0, "v0": 0.0, "tke": 0.0, "lapse": 0.0},
            sim=self.sim,
            defaults={},
        )
        section.generate_prof()
        section.write_prof()
        # Source file (z, u, v, thl, qt); z starts at 1.0 > 0 → surface point is prepended
        prdata_raw = np.array([
            [1.0, 2.0, 0.5, 300.0, 0.010],
            [3.0, 4.0, 1.0, 302.0, 0.009],
            [6.0, 6.0, 1.5, 306.0, 0.008],
        ])
        prdata_path = self.workdir / "prdata.dat"
        np.savetxt(prdata_path, prdata_raw, header="z u v thl qt")
        with self.assertWarns(UserWarning):
            section.update_prof_from_nudge_data(str(prdata_path))
        full_data = np.vstack(([0.0, 0.0, 0.0, 293.0, 0.0], prdata_raw))
        zt = self.sim.zt
        np.testing.assert_allclose(section.sim.pr[:, 1],
                                    CubicSpline(full_data[:, 0], full_data[:, 3])(zt), rtol=1e-5)
        np.testing.assert_allclose(section.sim.pr[:, 2],
                                    CubicSpline(full_data[:, 0], full_data[:, 4])(zt), rtol=1e-5)
        np.testing.assert_allclose(section.sim.pr[:, 3],
                                    CubicSpline(full_data[:, 0], full_data[:, 1])(zt), rtol=1e-5)
        np.testing.assert_allclose(section.sim.pr[:, 4],
                                    CubicSpline(full_data[:, 0], full_data[:, 2])(zt), rtol=1e-5)

    def test_update_prof_from_nudge_data_skips_surface_point_when_z_starts_at_zero(self):
        from scipy.interpolate import CubicSpline
        section = ForcingSection(
            "forcing",
            {"thl0": 290.0, "qt0": 0.0, "u0": 0.0, "v0": 0.0, "tke": 0.0, "lapse": 0.0},
            sim=self.sim,
            defaults={},
        )
        section.generate_prof()
        section.write_prof()
        # Source data starts at z=0 → no surface point prepended (would create duplicate otherwise)
        prdata_raw = np.array([
            [0.0, 0.0, 0.0, 293.0, 0.000],
            [3.0, 3.0, 0.5, 300.0, 0.010],
            [6.0, 5.0, 1.0, 304.0, 0.009],
        ])
        prdata_path = self.workdir / "prdata_zero.dat"
        np.savetxt(prdata_path, prdata_raw, header="z u v thl qt")
        with self.assertWarns(UserWarning):
            section.update_prof_from_nudge_data(str(prdata_path))
        expected_thl = CubicSpline(prdata_raw[:, 0], prdata_raw[:, 3])(self.sim.zt)
        np.testing.assert_allclose(section.sim.pr[:, 1], expected_thl, rtol=1e-5)

    # ---- update_prof_from_driver -----------------------------------------------

    def test_update_prof_from_driver_missing_driver_file_warns_and_returns(self):
        section = ForcingSection(
            "forcing",
            {"thl0": 295.0, "qt0": 0.0, "u0": 5.0, "v0": 0.0, "tke": 0.0, "lapse": 0.0},
            sim=self.sim,
            defaults={},
        )
        section.generate_prof()
        section.write_prof()
        with self.assertWarns(UserWarning):
            section.update_prof_from_driver("999", str(self.workdir / "nonexistent"), 1)
        # prof.inp must remain unchanged — u column still 5.0
        data = np.loadtxt(self.workdir / "prof.inp.321", skiprows=2)
        np.testing.assert_allclose(data[:, 3], 5.0, rtol=1e-5)

if __name__ == "__main__":
    unittest.main()
