"""Tier-1 tests for the moist thermodynamics / MO surface / radiation port."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

import thermo_case as tc


class _SolverTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not tc.UDALES_BUILD.is_file():
            raise unittest.SkipTest(f"no executable at {tc.UDALES_BUILD}")
        cls._tmp = tempfile.TemporaryDirectory(prefix="thermo_")
        cls.tmp = Path(cls._tmp.name)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def run_ok(self, run_dir, nprocs=1):
        res = tc.run(run_dir, nprocs)
        if res.returncode != 0:
            self.fail("solver failed:\n" + "\n".join((res.stdout + res.stderr).splitlines()[-40:]))
        return res


class TestInitChecks(_SolverTest):
    def test_sentinel_thls_aborts(self):
        d = tc.write_case(self.tmp / "sentinel", groups={"BC": {"thls": "-1."}})
        res = tc.run(d)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn("thls", res.stdout + res.stderr)

    def test_dry_case_without_qts_runs(self):
        d = tc.write_case(self.tmp / "dryqts",
                          groups={"PHYSICS": {"lmoist": ".false."}, "BC": {"qts": "-1."}})
        self.run_ok(d)


class TestCalthvGuard(_SolverTest):
    def test_uniform_saturated_column_does_not_divide_by_zero(self):
        # Default RF01 profile: uniform saturated mixed layer, zero vertical gradients,
        # so del_thv_dry - del_thv_sat == 0 in calthv (SIGFPE in the debug build).
        d = tc.write_case(self.tmp / "calthv", groups={"RUN": {"runtime": "0.5"}})
        self.run_ok(d)


class TestSurfaceSign(_SolverTest):
    """One step, laminar column with only a surface flux (no wind; neutral momentum wall function)."""

    def _first_level_change(self, name, wtsurf, wqsurf):
        d = tc.write_case(self.tmp / name, u=0.0, v=0.0, groups={
            "RUN": {"runtime": "0.5", "dtmax": "1.0"},
            "BC": {"BCbotm": "3", "wtsurf": str(wtsurf), "wqsurf": str(wqsurf)},
            "NAMSUBGRID": {"loneeqn": ".false.", "lvreman": ".true."},
            "OUTPUT": {"lbldump": ".true.", "tstatsdump": "0.5"}},
            thl=np.full(64, 289.0), qt=np.full(64, 5.0e-3))
        self.run_ok(d)
        with Dataset(d / "bldump_z.310.nc") as z:
            return float(z["thl"][-1, 0]) - 289.0, float(z["qt"][-1, 0]) - 5.0e-3

    def test_positive_flux_is_downward(self):
        # uDALES convention (D1): prescribed boundary fluxes are positive in -z,
        # so positive wtsurf/wqsurf cool/dry the first level by flux*dt/dz.
        dthl, dqt = self._first_level_change("sign", 1.0, 1.0e-3)
        dz, dt = 20.0, 1.0
        self.assertAlmostEqual(dthl, -1.0 * dt / dz, delta=2e-2 * 1.0 * dt / dz)
        self.assertAlmostEqual(dqt, -1.0e-3 * dt / dz, delta=2e-2 * 1.0e-3 * dt / dz)


class TestMO(_SolverTest):
    """Monin-Obukhov facets (iwallmom = 4) on a flat facet floor, laminar homogeneous column.

    Tolerances are 1e-3 relative: the records are float32 (1e-7) and the remaining
    differences are the drift of the air state within the last step (the drag slows the
    wind by ~5e-4 per step, and u*, L are computed from the state of the last RK substep,
    the recorded u1/v1/thl/qt from the end of the step; see test_stress_applied)."""

    Z0 = 2.e-4
    Z1 = 10.0

    def _mo_case(self, name, bctfz="-0.05", bcqfz="-5.0e-5", extra=None, runtime="19.5", u=7.0, v=-5.5, geo=True):
        # upward heat/moisture flux = negative in the uDALES (-z positive) convention.
        # 20 fixed 1-s steps so the laminar (tskin, L) pair of every section has converged.
        groups = {"RUN": {"runtime": runtime}, "BC": {"bctfz": bctfz, "bcqfz": bcqfz},
                  "OUTPUT": {"lbldump": ".true.", "tstatsdump": "1."}}
        for g, kv in (extra or {}).items():
            groups.setdefault(g, {}).update(kv)
        d = tc.write_case(self.tmp / name, u=u, v=v, floor={"z0": self.Z0, "z0h": self.Z0}, groups=groups,
                          ug=None if geo else 0.0, vg=None if geo else 0.0)
        res = self.run_ok(d)
        return d, res.stdout + res.stderr

    @staticmethod
    def _last(d):
        with Dataset(d / "bldump_t.310.nc") as t, Dataset(d / "bldump_z.310.nc") as z:
            r = {v: np.asarray(t[v][:], dtype=float) for v in ("ustar", "obl", "thls", "qts", "u1", "v1")}
            r["thl1"] = np.asarray(z["thl"][:, 0], dtype=float)
            r["qt1"] = np.asarray(z["qt"][:, 0], dtype=float)
        return r

    def test_unstable_ustar_obl_match_reference(self):
        d, out = self._mo_case("mo_unstable")
        self.assertIn("iwallmom = 4: 1 facets use Monin-Obukhov", out)
        r = self._last(d)
        ustar, obl, thls, qts = r["ustar"][-1], r["obl"][-1], r["thls"][-1], r["qts"][-1]
        thl1, qt1, u1, v1 = r["thl1"][-1], r["qt1"][-1], r["u1"][-1], r["v1"][-1]
        self.assertGreater(ustar, 0.01)
        self.assertLess(obl, 0.0)                         # upward buoyancy flux -> unstable
        self.assertGreater(thls, thl1)                    # surface warmer than air
        # near-converged: L follows the decelerating wind, about 2e-3 per step
        self.assertAlmostEqual(obl, r["obl"][-2], delta=5e-3 * abs(obl))
        # Reference: the same iteration on the recorded converged state (float32 records)
        L = tc.obl_iter(thl1, qt1, thls, qts, self.Z1, self.Z0, self.Z0, u1, v1, -0.01)
        self.assertAlmostEqual(obl, L, delta=1e-3 * abs(L))
        cm = tc.FKAR ** 2 / (np.log(self.Z1 / self.Z0) - tc.psim(self.Z1 / L) + tc.psim(self.Z0 / L)) ** 2
        self.assertAlmostEqual(ustar, np.sqrt(cm) * np.hypot(u1, v1), delta=1e-3 * ustar)
        # the flux definition of L closes the loop: L = -u*^3 thv_s / (kappa g w'thv'),
        # w'thv' = (1 + 0.61 qs) w'thl' + 0.61 thl_s w'qt' (upward-positive fluxes)
        wthv = (1 + (tc.RV / tc.RD - 1) * qts) * 0.05 + (tc.RV / tc.RD - 1) * thls * 5.0e-5
        thvs = thls * (1 + (tc.RV / tc.RD - 1) * qts)
        self.assertAlmostEqual(obl, -ustar ** 3 * thvs / (tc.FKAR * tc.GRAV * wthv), delta=1e-2 * abs(obl))

    def test_neutral_ustar_is_log_law(self):
        d, _ = self._mo_case("mo_neutral", bctfz="0.", bcqfz="0.")
        r = self._last(d)
        self.assertEqual(r["obl"][-1], 1.e6)              # Rib = 0: neutral cap
        expect = tc.FKAR * np.hypot(r["u1"][-1], r["v1"][-1]) / np.log(self.Z1 / self.Z0)
        self.assertAlmostEqual(r["ustar"][-1], expect, delta=1e-3 * expect)

    def test_stress_applied(self):
        """First step: the first-level wind loses u*^2 dt / dzf, along the wind (neutral and
        unstable). One RK step: the u*^2 of the three substeps differ by the drift only.
        ug = vg = 0, since uDALES applies the geostrophic pressure gradient of lscale.inp
        (f x ug, default xlat) even with lcoriol off."""
        for name, bctfz in (("mo_stress_neutral", "0."), ("mo_stress_unstable", "-0.05")):
            with self.subTest(case=name):
                d, _ = self._mo_case(name, bctfz=bctfz, bcqfz="0.", runtime="0.5", geo=False,
                                     extra={"NAMSUBGRID": {"loneeqn": ".false.", "lvreman": ".true."},
                                            "OUTPUT": {"tstatsdump": "0.5"}})
                r = self._last(d)
                du, dv = r["u1"][-1] - 7.0, r["v1"][-1] - (-5.5)
                tau = r["ustar"][-1] ** 2 * 1.0 / 20.0      # u*^2 dt / dzf
                self.assertAlmostEqual(np.hypot(du, dv), tau, delta=1e-3 * tau + 1e-6)
                self.assertAlmostEqual(du / dv, 7.0 / -5.5, delta=1e-3 * 7.0 / 5.5)

    def test_prescribed_ustar(self):
        d, _ = self._mo_case("mo_ustin", extra={"WALLS": {"ustin": "0.28"}})
        self.assertAlmostEqual(self._last(d)["ustar"][-1], 0.28, places=6)   # float32 record

    def test_moisture_requires_prescribed_facet_flux(self):
        d = tc.write_case(self.tmp / "mo_iwallmoist", u=7.0, v=-5.5, floor={"z0": 2e-4, "z0h": 2e-4},
                          groups={"RUN": {"runtime": "1."}, "WALLS": {"iwallmoist": "2"},
                                  "BC": {"bctfz": "-0.05", "bcqfz": "-5.0e-5"}})
        res = tc.run(d)
        self.assertNotEqual(res.returncode, 0)
        self.assertIn("iwallmoist", res.stdout + res.stderr)

    def test_slope_selects_mo_or_wall_function(self):
        """mo_maxslope = 30 deg: a facet tilted by 20 deg takes Monin-Obukhov, one tilted by
        60 deg the neutral wall function. The floor stays flat in the grid (only the normal,
        type and section geometry enter the wall functions); the west half is facet 1, the
        east half facet 2. The bldump u* is the mean over the Monin-Obukhov sections only,
        so it must be the neutral u* of facet 1 alone."""
        t20, t60 = np.radians(20.0), np.radians(60.0)
        normals = [(np.sin(t20), 0.0, np.cos(t20)), (np.sin(t60), 0.0, np.cos(t60))]
        d = tc.write_case(self.tmp / "mo_slope", u=7.0, v=0.0,
                          floor={"z0": self.Z0, "z0h": self.Z0, "normals": normals,
                                 "facet_of": lambda i, j: 1 if i <= 8 else 2},
                          groups={"RUN": {"runtime": "0.5"}, "BC": {"bctfz": "0.", "bcqfz": "0."},
                                  "OUTPUT": {"lbldump": ".true.", "tstatsdump": "0.5"}})
        res = self.run_ok(d)
        self.assertIn("iwallmom = 4: 1 facets use Monin-Obukhov (slope <=  30.0 deg), 1 use the neutral wall function",
                      res.stdout + res.stderr)
        with Dataset(d / "bldump_t.310.nc") as t:
            ustar = float(t["ustar"][-1])
        # tangential wind on facet 1: U minus its normal component; L = 1e6 (no flux).
        # The section reference point is z1 above the floor (the section distance).
        n1 = np.array(normals[0])
        uvec = np.array([7.0, 0.0, 0.0])
        utan = np.linalg.norm(uvec - uvec.dot(n1) * n1)
        expect = tc.FKAR * utan / np.log(self.Z1 / self.Z0)
        self.assertAlmostEqual(ustar, expect, delta=2e-3 * expect)   # u1 drifts within the first step


RV_RD = 461.5 / 287.04
RLV, CP = 2.53e6, 1004.0   # check modglobal: rlv, cp


def ql_udales(thl, qt, exn):
    """uDALES modthermodynamics:thermo (lqlnr = .false.): one linearised step, Tetens."""
    rd, rv, pref0 = 287.04, 461.5, 1.0e5
    tl = thl * exn
    p = pref0 * exn ** (CP / rd)
    es = 610.78 * np.exp(17.27 * (tl - 273.16) / (tl - 35.86))
    qsl = rd / rv * es / (p - (1 - rd / rv) * es)
    b1 = RLV ** 2 / (tl ** 2 * CP * rv)
    qs = qsl * (1 + b1 * qt) / (1 + b1 * qsl)
    return np.maximum(qt - qs, 0.0)


class TestCalthv(_SolverTest):
    def test_surface_level_is_zero(self):
        d = tc.write_case(self.tmp / "calthv_wf", groups={
            "RUN": {"runtime": "0.5"},
            "OUTPUT": {"lfielddump": ".true.", "tfielddump": "0.5", "fieldvars": "'dt'"}})
        self.run_ok(d)
        with Dataset(d / "fielddump.000.000.310.nc") as f:
            self.assertTrue(np.all(np.abs(np.asarray(f["dthvdz"][-1])[0]) <= 2e-10))  # float32(eps1)


class TestBldump(_SolverTest):
    def test_cloudy_column_diagnostics(self):
        # 3 steps, laminar RF01 column: horizontally homogeneous, so column = slab mean
        # A perfectly uniform saturated mixed layer divides by zero in calthv
        # (chi_sat, del_thv_dry = del_thv_sat = 0), which the debug build traps.
        # A weak slope below the inversion avoids that without changing the diagnostics.
        zf = (np.arange(64) + 0.5) * 20.0
        thl, qt = tc.rf01_profiles(zf)
        thl = thl + np.where(zf <= 840.0, 1.0e-3 * zf / 20.0, 0.0)
        qt = qt - np.where(zf <= 840.0, 1.0e-7 * zf / 20.0, 0.0)
        d = tc.write_case(self.tmp / "stat", thl=thl, qt=qt, groups={
            "RUN": {"runtime": "2.5"},
            "OUTPUT": {"lbldump": ".true.", "tstatsdump": "1.", "zi_qt_thres": "8.e-3"}})
        self.run_ok(d)
        with Dataset(d / "bldump_t.310.nc") as t, Dataset(d / "bldump_z.310.nc") as z:
            self.assertGreaterEqual(len(t["time"]), 2)
            ql = z["ql"][-1, :]
            rhof = z["rhof"][-1, :]
            dz = 20.0
            lwp_np = float(np.sum(rhof * ql * dz))
            self.assertAlmostEqual(float(t["lwp_bar"][-1]), lwp_np, delta=1e-6 * max(lwp_np, 1e-12))
            self.assertGreater(lwp_np, 0.0)             # RF01 mixed layer is saturated aloft
            self.assertEqual(float(t["cfrac"][-1]), 1.0)
            self.assertAlmostEqual(float(t["zi_qt"][-1]), 840.0, delta=20.0)
            self.assertAlmostEqual(float(t["zi"][-1]), 840.0, delta=20.0)
            # no radiation in this case: the radiation profiles stay zero
            for name in ("thlprad", "lwu", "lwd"):
                self.assertTrue(np.all(z[name][-1, :] == 0.0))

    def test_dry_column_has_no_moist_diagnostics(self):
        d = tc.write_case(self.tmp / "drystat", groups={
            "RUN": {"runtime": "1.5"},
            "PHYSICS": {"lmoist": ".false."}, "BC": {"qts": "-1."},
            "OUTPUT": {"lbldump": ".true.", "tstatsdump": "1."}})
        self.run_ok(d)
        with Dataset(d / "bldump_t.310.nc") as t:
            self.assertGreaterEqual(len(t["time"]), 1)
            self.assertNotIn("zi_qt", t.variables)      # zi_qt_thres not set
            self.assertEqual(float(t["lwp_bar"][-1]), 0.0)

    def _ql_record(self):
        """Record ql and ql recomputed from the same record's thl, qt, exnf with the uDALES
        one-step adjustment (Tetens, b1 = rlv^2/(Tl^2 cp rv)). thermo uses the pressure of
        the previous diagfld, so agreement is ~1e-5, not exact."""
        if not hasattr(self, "_qlrec"):
            zf = (np.arange(64) + 0.5) * 20.0
            thl, qt = tc.rf01_profiles(zf)
            thl = thl + np.where(zf <= 840.0, 1.0e-3 * zf / 20.0, 0.0)
            qt = qt - np.where(zf <= 840.0, 1.0e-7 * zf / 20.0, 0.0)
            d = tc.write_case(self.tmp / "stat_ql", thl=thl, qt=qt, groups={
                "RUN": {"runtime": "2.5"},
                "OUTPUT": {"lbldump": ".true.", "tstatsdump": "1."}})
            self.run_ok(d)
            with Dataset(d / "bldump_z.310.nc") as z:
                ql, thl_r, qt_r, exn_r = (np.asarray(z[v][-1, :], dtype=float)
                                          for v in ("ql", "thl", "qt", "exnf"))
            type(self)._qlrec = (ql, ql_udales(thl_r, qt_r, exn_r))
        return self._qlrec

    def test_record_ql_matches_record_thl_qt(self):
        # Checks rlv (a wrong rlv changes ql by ~12 %), that ql0(k) belongs to level k (#403:
        # it used to hold the ql of level k+1) and that the record's ql belongs to the same
        # state as its thl/qt (bldump after thermodynamics).
        ql, ref = self._ql_record()
        sat = ref > 0.1 * ref.max()
        self.assertGreater(int(sat.sum()), 5)
        dev = np.abs(ql[sat] - ref[sat]) / ref[sat]
        print(f"\n[ql] {int(sat.sum())} saturated levels, record ql vs recomputed ql: "
              f"max rel dev {dev.max():.3e}")
        np.testing.assert_allclose(ql[sat], ref[sat], rtol=1e-4)


def radpar_reference(absorber, rho, exn, dz, f0, f1, kappa):
    """Stevens et al. (2005) / DALES radpar longwave, one column. Returns thl tendency [K/s]."""
    n = absorber.size
    path = rho * absorber * dz
    lwpt = np.concatenate([np.cumsum(path[::-1])[::-1], [0.0]])     # above face k (incl. cell k)
    lwpb = np.concatenate([[0.0], np.cumsum(path)])                 # below face k
    flux = f0 * np.exp(-kappa * lwpt) + f1 * np.exp(-kappa * lwpb)  # net upward, faces 0..n
    return -(flux[1:] - flux[:-1]) / (rho * CP * exn * dz)


class TestRadiation(_SolverTest):
    RAD = {"iradiation": "2", "rad_longw": ".true.", "rka": "85.", "dlwtop": "70.", "dlwbot": "22."}

    def _run(self, name, rad, nprocs=1, **case):
        groups = {"RUN": {"runtime": "0.5"}, "NAMRADIATION": rad,
                  "OUTPUT": {"lbldump": ".true.", "tstatsdump": "0.5"}}
        groups.update(case.pop("groups", {}))
        d = tc.write_case(self.tmp / name, groups=groups, **case)
        self.run_ok(d, nprocs)
        return d

    def _report(self, label, got, ref):
        sel = np.abs(ref) > 1e-2 * np.abs(ref).max()
        dev = np.abs(got - ref)[sel] / np.abs(ref)[sel]
        print(f"\n[{label}] max|thlprad|={np.abs(got).max():.4e} max rel dev (|ref|>1% of max)={dev.max():.3e}")

    def test_liquid_longwave_matches_formula(self):
        d = self._run("rad_liq", self.RAD)
        with Dataset(d / "bldump_z.310.nc") as z:
            ql, rho, exn = (np.asarray(z[v][-1]) for v in ("ql", "rhof", "exnf"))
            got = np.asarray(z["thlprad"][-1])
        ref = radpar_reference(ql, rho, exn, 20.0, 70.0, 22.0, 85.0)
        self._report("liquid", got, ref)
        # thlprad uses start-of-step ql; the record's ql is end-of-step (~0.3 % change at
        # cloud top). A sign or face off-by-one error is O(1), so 1e-2 still discriminates.
        np.testing.assert_allclose(got, ref, rtol=1e-2, atol=1e-2 * np.abs(ref).max())
        self.assertLess(got.min(), -1e-4)          # cloud-top cooling present

    SMOKE = dict(RAD, rad_absorber="'sv'", rad_absorber_sv="1", rka="0.02", dlwtop="60.", dlwbot="0.")

    def _smoke_run(self, name, runtime="0.5", bldump=None):
        zf = (np.arange(64) + 0.5) * 20.0
        smoke = np.where(zf < 700, 1.0, 0.0)
        th = {"lbldump": ".true.", "tstatsdump": "0.5"}
        th.update(bldump or {})
        d = self._run(name, self.SMOKE, groups={"RUN": {"runtime": runtime}, "PHYSICS": {"lmoist": ".false."},
                                                "BC": {"qts": "-1."}, "SCALARS": {"nsv": "1"},
                                                "OUTPUT": th},
                      thl=np.where(zf < 700, 288.0, 295.0), qt=np.zeros(64), sv=smoke[:, None])
        return d, smoke

    def test_smoke_longwave_matches_formula(self):
        d, smoke = self._smoke_run("rad_smoke")
        with Dataset(d / "bldump_z.310.nc") as z:
            rho, exn = np.asarray(z["rhof"][-1]), np.asarray(z["exnf"][-1])
            got = np.asarray(z["thlprad"][-1])
        ref = radpar_reference(smoke, rho, exn, 20.0, 60.0, 0.0, 0.02)
        self._report("smoke", got, ref)
        np.testing.assert_allclose(got, ref, rtol=1e-2, atol=1e-2 * np.abs(ref).max())

    def test_smoke_identical_to_smoke_reference(self):
        # smoke_reference.txt: the same case recorded with an earlier smoke-specific longwave
        # implementation (Debug build). rad_absorber = 'sv' must reproduce it exactly.
        d, _ = self._smoke_run("rad_smoke_pre16", runtime="2.5",
                               bldump={"zi_sv_index": "1", "zi_sv_thres": "0.5"})
        ref = np.loadtxt(Path(__file__).with_name("smoke_reference.txt"), dtype=np.float32)
        with Dataset(d / "bldump_z.310.nc") as z, Dataset(d / "bldump_t.310.nc") as t:
            for col, v in enumerate(("thl", "thlprad", "lwd", "lwu")):
                with self.subTest(var=v):
                    np.testing.assert_array_equal(np.asarray(z[v][-1], dtype=np.float32), ref[:, col])
            np.testing.assert_array_equal(np.asarray(t["zi_sv"][:]), 700.0)

    def test_zi_sv_written_only_when_configured(self):
        d, _ = self._smoke_run("zi_sv_off")
        with Dataset(d / "bldump_t.310.nc") as t:
            self.assertNotIn("zi_sv", t.variables)
            self.assertNotIn("zi_qt", t.variables)

    # --- DYCOMS free-troposphere term (case hook dycoms_rf01)
    D = 3.75e-6
    DYCOMS = {"NAMCASES": {"case_hooks": "'dycoms_rf01'"},
              "NAMCASE_DYCOMS": {"dycoms_D": str(D), "dycoms_qt_thres": "8.e-3"}}

    @staticmethod
    def _ft_flux(z, zi, D):
        """F(z) / (rho_i cp) of Stevens et al. (2005), alpha_z = 1; 0 at and below zi."""
        s = np.maximum(z - zi, 0.0)
        return np.where(z > zi, D * (s ** (4 / 3) / 4 + zi * np.cbrt(s)), 0.0)

    def _ft_term(self):
        """FT tendency of the first step: hook run minus no-hook run (same initial state,
        so the radpar parts are identical), with zi from the initial qt profile."""
        if not hasattr(type(self), "_ft"):
            a = self._run("ft_hook", self.RAD, groups=self.DYCOMS)
            b = self._run("ft_none", self.RAD)
            with Dataset(a / "bldump_z.310.nc") as za, Dataset(b / "bldump_z.310.nc") as zb:
                ft = np.asarray(za["thlprad"][0], dtype=float) - np.asarray(zb["thlprad"][0], dtype=float)
            zf = (np.arange(64) + 0.5) * 20.0
            qt = tc.rf01_profiles(zf)[1]
            k = np.where(qt > 8e-3)[0].max()
            zi = zf[k] + (8e-3 - qt[k]) * (zf[k + 1] - zf[k]) / (qt[k + 1] - qt[k])   # 832.67 m
            type(self)._ft = (ft, zi)
        return self._ft

    def test_dycoms_flux_form_matches_pointwise_form(self):
        # The flux form is the cell mean of the pointwise form, so their difference is
        # dz^2/24 f''/f ~ (dz^2/24)(10/9)/(z-zi)^2 relative (the zi (z-zi)^(-2/3) term
        # dominates): 5e-3 at 3 dz above zi, below 1e-3 from 7 dz. A sign error or a
        # missing term is O(1).
        ft, zi = self._ft_term()
        dz = 20.0
        zf = (np.arange(64) + 0.5) * dz
        sel = zf - zi > 3 * dz
        pw = -self.D / 3 * (np.cbrt(zf[sel] - zi) + zi * (zf[sel] - zi) ** (-2 / 3))
        rel = np.abs(ft[sel] - pw) / np.abs(pw)
        bound = 1.2 * dz ** 2 / 24 * (10 / 9) / (zf[sel] - zi) ** 2 + 1e-5
        print(f"\n[dycoms] flux vs pointwise, z - zi > 3 dz: max rel {rel.max():.3e}; "
              f"z - zi > 7 dz: {rel[zf[sel] - zi > 7 * dz].max():.3e}")
        self.assertTrue(np.all(rel < bound), msg=f"rel {rel} bound {bound}")
        self.assertLess(rel[zf[sel] - zi > 7 * dz].max(), 1e-3)
        self.assertTrue(np.all(ft[zf + dz / 2 > zi] < 0))        # cooling in and above the zi cell
        np.testing.assert_array_equal(ft[zf + dz / 2 <= zi], 0.0)   # nothing below

    def test_dycoms_flux_form_conserves_column_integral(self):
        # sum(thlprad dz) = -(F(top) - F(zi)) / (rho_i cp); F(zi) = 0. Float32 records.
        ft, zi = self._ft_term()
        top = 64 * 20.0
        integral = (ft * 20.0).sum()
        ref = -self._ft_flux(np.array(top), zi, self.D)
        print(f"\n[dycoms] column integral {integral:.6e} vs -F(top) {float(ref):.6e}")
        np.testing.assert_allclose(integral, ref, rtol=1e-5)

    def test_dycoms_term_balances_subsidence(self):
        # 300 s so the free-troposphere subsidence signal (~2e-2 K) is far above
        # float32 resolution (3e-5 K at 297 K).
        zf = (np.arange(64) + 0.5) * 20.0
        D = self.D
        common = dict(wfls=-D * zf, groups={"RUN": {"runtime": "299.5"},
                                            "OUTPUT": {"lbldump": ".true.", "tstatsdump": "300."}})
        a = self._run("ft_both", self.RAD, groups=dict(common["groups"], **self.DYCOMS), wfls=common["wfls"])
        b = self._run("ft_subs", {"iradiation": "0"}, **common)
        with Dataset(a / "bldump_z.310.nc") as za, Dataset(b / "bldump_z.310.nc") as zb:
            dthl_both = np.asarray(za["thl"][-1], dtype=float) - tc.rf01_profiles(zf)[0]
            dthl_subs = np.asarray(zb["thl"][-1], dtype=float) - tc.rf01_profiles(zf)[0]
        # The top cell is excluded: there the subsidence stencil is one-sided/boundary-
        # dominated (subs run: -1.9e-3 K, cooling), while the FT term keeps cooling it
        # (-1.0e-2 K), so it measures the top boundary, not the balance.
        ft = (zf > 840.0 + 150.0) & (zf < zf[-1])
        ratio = np.abs(dthl_both[ft]).max() / np.abs(dthl_subs[ft]).max()
        print(f"\n[balance] max|dthl_both|/max|dthl_subs| in FT = {ratio:.4f}")
        # Residual sources: one-sided subsidence differences, O(dz/(z-zi)), and the
        # 832.7 vs 840 m zi offset. A sign error gives a ratio of ~2.
        self.assertLess(np.abs(dthl_both[ft]).max(), 0.25 * np.abs(dthl_subs[ft]).max())


class TestCaseHooks(_SolverTest):
    """&NAMCASES registry (modcases) and the checks of the data-driven options."""

    RAD = TestRadiation.RAD

    def _abort(self, name, groups, expect, **case):
        d = tc.write_case(self.tmp / name, groups=groups, **case)
        res = tc.run(d)
        out = res.stdout + res.stderr
        self.assertNotEqual(res.returncode, 0)
        self.assertIn(expect, out)
        return out

    def test_unknown_hook_aborts_and_lists_valid_names(self):
        out = self._abort("hook_unknown", {"NAMCASES": {"case_hooks": "'dycoms_rf02'"}},
                          "unknown case hook in case_hooks: dycoms_rf02")
        self.assertIn("Valid case hooks: dycoms_rf01", out)

    def test_dycoms_without_parameters_aborts(self):
        self._abort("hook_noparams", {"NAMRADIATION": self.RAD,
                                      "NAMCASES": {"case_hooks": "'dycoms_rf01'"}},
                    "requires dycoms_D > 0 and dycoms_qt_thres > 0")

    def test_dycoms_without_radiation_aborts(self):
        self._abort("hook_norad", dict(TestRadiation.DYCOMS), "requires iradiation = 2")

    def test_sv_absorber_without_index_aborts(self):
        self._abort("abs_noidx", {"NAMRADIATION": dict(self.RAD, rad_absorber="'sv'"),
                                  "SCALARS": {"nsv": "1"}},
                    "requires 1 <= rad_absorber_sv <= nsv", sv=np.ones((64, 1)))

    def test_dycoms_moist_to_top_column_is_skipped(self):
        # qt above the threshold up to the top level: no inversion inside the domain, so
        # no extrapolation; the column gets no FT term and init warns.
        zf = (np.arange(64) + 0.5) * 20.0
        d = tc.write_case(self.tmp / "hook_moisttop", thl=tc.rf01_profiles(zf)[0], qt=np.full(64, 9.0e-3),
                          groups=dict(TestRadiation.DYCOMS, NAMRADIATION=self.RAD,
                                      OUTPUT={"lbldump": ".true.", "tstatsdump": "0.5"}))
        res = self.run_ok(d)
        self.assertIn("the free-troposphere term is off", res.stdout + res.stderr)

    def test_duplicate_hook_aborts(self):
        groups = dict(TestRadiation.DYCOMS, NAMRADIATION=self.RAD)
        groups["NAMCASES"] = {"case_hooks": "'dycoms_rf01, dycoms_rf01'"}
        self._abort("hook_dup", groups,
                    "case hook listed twice in case_hooks: dycoms_rf01")

    def test_radiation_without_optics_aborts(self):
        self._abort("rad_noopt", {"NAMRADIATION": {"iradiation": "2"}},
                    "&NAMRADIATION: iradiation = 2 requires rka > 0 and dlwtop > 0 or dlwbot > 0")

    def test_invalid_absorber_aborts(self):
        # strict on the whole value: 'qlx' must not be read as 'ql'
        self._abort("abs_bad", {"NAMRADIATION": dict(self.RAD, rad_absorber="'qlx'")},
                    "rad_absorber must be 'ql' or 'sv'; got 'qlx'")

    def test_sv_absorber_index_above_nsv_aborts(self):
        self._abort("abs_big", {"NAMRADIATION": dict(self.RAD, rad_absorber="'sv'", rad_absorber_sv="2"),
                                "SCALARS": {"nsv": "1"}},
                    "requires 1 <= rad_absorber_sv <= nsv; got 2", sv=np.ones((64, 1)))

    def test_sv_index_with_ql_absorber_aborts(self):
        self._abort("abs_qlidx", {"NAMRADIATION": dict(self.RAD, rad_absorber_sv="1"),
                                  "SCALARS": {"nsv": "1"}},
                    "rad_absorber_sv is set but rad_absorber = 'ql'", sv=np.ones((64, 1)))

    def test_zi_sv_index_above_nsv_aborts(self):
        self._abort("zisv_big", {"SCALARS": {"nsv": "1"},
                                 "OUTPUT": {"lbldump": ".true.", "zi_sv_index": "2",
                                                   "zi_sv_thres": "0.5"}},
                    "zi_sv_index > 0 requires", sv=np.ones((64, 1)))

    def test_zi_qt_without_moisture_aborts(self):
        self._abort("ziqt_dry", {"PHYSICS": {"lmoist": ".false."}, "BC": {"qts": "-1."},
                                 "OUTPUT": {"lbldump": ".true.", "zi_qt_thres": "8.e-3"}},
                    "zi_qt_thres must be >= 0, and > 0 requires lmoist")

    def test_zi_sv_without_threshold_aborts(self):
        self._abort("zisv_nothres", {"SCALARS": {"nsv": "1"},
                                     "OUTPUT": {"lbldump": ".true.", "zi_sv_index": "1"}},
                    "zi_sv_index > 0 requires", sv=np.ones((64, 1)))

    def test_bldump_without_interval_aborts(self):
        self._abort("bld_noint", {"OUTPUT": {"lbldump": ".true.", "tstatsdump": "0."}},
                    "requires tstatsdump > 0")

    def test_negative_radiation_parameters_abort(self):
        self._abort("rad_negflux", {"NAMRADIATION": dict(self.RAD, dlwbot="-5.")},
                    "dlwtop and dlwbot must be >= 0")
        self._abort("rad_negtime", {"NAMRADIATION": dict(self.RAD, timerad="-1.")},
                    "timerad must be >= 0")

    def test_unreadable_namelist_aborts_once_on_all_ranks(self):
        # a read error in a new namelist group is collective: rank 0 reports once, all ranks abort
        d = tc.write_case(self.tmp / "nml_bad", groups={
            "RUN": {"nprocx": "2"}, "NAMRADIATION": {"iradiation": "banana"}})
        res = tc.run(d, 2)
        out = res.stdout + res.stderr
        self.assertNotEqual(res.returncode, 0)
        self.assertEqual(out.count("Problem in namoptions NAMRADIATION"), 1)

    def test_case_group_without_hook_warns(self):
        # &NAMCASE_DYCOMS is present but 'dycoms_rf01' is not in case_hooks: the run goes on, with a warning
        d = tc.write_case(self.tmp / "hook_unlisted", groups={
            "NAMCASE_DYCOMS": {"dycoms_D": "3.75e-6", "dycoms_qt_thres": "8.e-3"}})
        res = self.run_ok(d)
        self.assertIn("&NAMCASE_DYCOMS is in namoptions but 'dycoms_rf01' is not in case_hooks", res.stdout + res.stderr)

    def test_no_warning_when_hook_listed(self):
        d = tc.write_case(self.tmp / "hook_listed", groups=dict(TestRadiation.DYCOMS, NAMRADIATION=self.RAD))
        res = self.run_ok(d)
        self.assertNotIn("is not in case_hooks", res.stdout + res.stderr)

    def test_no_hooks_identical_to_run_without_namelist(self):
        # case_hooks = '' must be a no-op: every bldump output identical to a run
        # without &NAMCASES (laminar, one rank).
        groups = {"RUN": {"runtime": "2.5"}, "NAMRADIATION": self.RAD,
                  "OUTPUT": {"lbldump": ".true.", "tstatsdump": "0.5"}}
        runs = {}
        for label, extra in (("absent", {}), ("empty", {"NAMCASES": {"case_hooks": "''"}})):
            d = tc.write_case(self.tmp / f"nohook_{label}", groups=dict(groups, **extra))
            self.run_ok(d)
            runs[label] = d
        for f in ("bldump_t.310.nc", "bldump_z.310.nc"):
            with Dataset(runs["absent"] / f) as a, Dataset(runs["empty"] / f) as b:
                self.assertEqual(set(a.variables), set(b.variables))
                for v in a.variables:
                    with self.subTest(file=f, var=v):
                        np.testing.assert_array_equal(np.asarray(b[v][:]), np.asarray(a[v][:]))


class TestDecomposition(_SolverTest):
    """Radiation + Monin-Obukhov facet floor + moist SGS give the same answer on 1, 2 and 4 ranks."""

    def test_splits_agree(self):
        groups = {"RUN": {"runtime": "2.5", "lrandomize": ".true."},
                  "BC": {"bctfz": "-0.012", "bcqfz": "-4.7e-5"},
                  "NAMRADIATION": TestRadiation.RAD,
                  "OUTPUT": {"lbldump": ".true.", "tstatsdump": "1."},
                  **TestRadiation.DYCOMS}
        zf = (np.arange(64) + 0.5) * 20.0
        out = {}
        for label, npx, npy in (("serial", 1, 1), ("x", 2, 1), ("xy", 2, 2)):
            g = {k: dict(v) for k, v in groups.items()}
            g["RUN"].update(nprocx=str(npx), nprocy=str(npy))
            d = tc.write_case(self.tmp / f"dec_{label}", groups=g, wfls=-3.75e-6 * zf,
                              floor={"z0": 2.e-4, "z0h": 2.e-4})
            self.run_ok(d, npx * npy)
            with Dataset(d / "bldump_t.310.nc") as t, Dataset(d / "bldump_z.310.nc") as z:
                out[label] = {v: np.asarray(t[v][:]) for v in ("lwp_bar", "ustar", "obl", "thls")}
                out[label].update({v: np.asarray(z[v][:]) for v in ("thl", "qt", "ql", "thlprad", "wthl")})
        for label in ("x", "xy"):
            for v, ref in out["serial"].items():
                diff = np.abs(out[label][v] - ref)
                print(f"\n[decomp] {label}:{v} max abs diff={diff.max():.3e} "
                      f"max rel={(diff / np.maximum(np.abs(ref), 1e-30)).max():.3e}")
                # Reduction order (allreduce slab means) differs between splits; the round-off
                # is amplified by the turbulent spin-up, most at the inversion. A decomposition
                # or halo bug would give O(1) differences from the first step, so 1e-5 relative
                # (and 1e-6 of the field maximum absolute) is still orders of magnitude below it.
                # Resolved fluxes (wthl, wqt) are near-cancelling means of w*thl terms (~290|w|), so
                # their round-off scales with the terms, not the result; a halo bug would be O(1).
                atol_frac = 1e-4 if v in ("wthl", "wqt") else 1e-6
                with self.subTest(split=label, var=v):
                    np.testing.assert_allclose(out[label][v], ref, rtol=1e-5, atol=atol_frac * max(np.abs(ref).max(), 1e-30),
                                               err_msg=f"{label}:{v}")


class TestSGSHaloExchange(_SolverTest):
    """#405: the one-equation SGS TKE (e120/e12m) halos must be exchanged between ranks. Without it
    the e12 halo cells keep their initial value, so the diffusion of e12 at tile edges (and through
    ekm every field) depends on the decomposition."""

    VARS = ("u0", "v0", "w0", "th")

    def _fields(self, name, npx, npy):
        d = tc.write_case(self.tmp / name, u=3.0, v=-2.0, groups={
            "RUN": {"runtime": "3.5", "lrandomize": ".true.", "randu": "0.5", "randthl": "0.3",
                    "nprocx": str(npx), "nprocy": str(npy)},
            "OUTPUT": {"lfielddump": ".true.", "tfielddump": "3.5", "fieldvars": "'u0,v0,w0,th'"},
            "NAMSUBGRID": {"loneeqn": ".true.", "lvreman": ".false.", "ldelta": ".false."}})
        self.run_ok(d, npx * npy)
        tiles = {}
        for f in sorted(d.glob("fielddump.*.*.310.nc")):
            ix, iy = int(f.name.split(".")[1]), int(f.name.split(".")[2])
            with Dataset(f) as nc:
                names = {"u0": "u", "v0": "v", "w0": "w", "th": "thl"}
                tiles[ix, iy] = {v: np.asarray(nc[names[v]][-1]) for v in self.VARS}   # (z, y, x)
        nx, ny = 1 + max(k[0] for k in tiles), 1 + max(k[1] for k in tiles)
        return {v: np.concatenate([np.concatenate([tiles[ix, iy][v] for ix in range(nx)], axis=2)
                                   for iy in range(ny)], axis=1) for v in self.VARS}

    def test_decomposition_independent(self):
        ref = self._fields("sgshalo_1x1", 1, 1)
        for npx, npy in ((2, 1), (1, 2), (2, 2)):
            out = self._fields(f"sgshalo_{npx}x{npy}", npx, npy)
            for v in self.VARS:
                diff = np.abs(out[v] - ref[v]).max()
                print(f"\n[sgshalo] {npx}x{npy}:{v} max abs diff={diff:.3e}")
                with self.subTest(split=(npx, npy), var=v):
                    # bit-identical with the fix; a missing e12 halo exchange gives ~1e-5 differences
                    np.testing.assert_array_equal(out[v], ref[v])


class TestSGSOneEqn(_SolverTest):
    """In-solver test (runmode 1008, tests.f90:tests_sgs_oneeqn) of the one-equation
    SGS model aligned with DALES dev 075ac4e: Kh/Km = 1/Prandtl = 3 in unstable cells and
    ch1 + ch2*zlt/delta in stable ones (zlt with grav/thvf), dissipation without the old
    factor 2, and no TKE source at kb (the wall functions set the near-wall state)."""

    def test_closure_and_sources(self):
        d = tc.write_case(self.tmp / "sgs1eq", groups={
            "RUN": {"runmode": "1008"},
            # Prandtl is left at its default, which must be exactly 1/3 (checked in tests.f90)
            "NAMSUBGRID": {"loneeqn": ".true.", "lvreman": ".false.", "ldelta": ".false."}})
        res = self.run_ok(d)
        out = res.stdout + res.stderr
        print("\n" + "\n".join(l for l in out.splitlines() if "worst" in l or "FAIL" in l))
        self.assertIn("ALL TESTS PASSED: tests_sgs_oneeqn", out)


class TestLstend(_SolverTest):
    """In-solver test (runmode 1009, tests.f90:tests_lstend) of the local-field subsidence:
    every cell against the DALES dev 075ac4e lstend upwind formula, for thl, qt, sv, u and v
    (lmomsubs), with whls changing sign and fields varying in x."""

    def test_local_upwind_subsidence(self):
        d = tc.write_case(self.tmp / "lstend", groups={
            "RUN": {"runmode": "1009"}, "SCALARS": {"nsv": "1"}}, sv=np.ones((64, 1)))
        res = self.run_ok(d)
        out = res.stdout + res.stderr
        print("\n" + "\n".join(l for l in out.splitlines() if "lstend levels" in l or "FAIL" in l))
        self.assertIn("ALL TESTS PASSED: tests_lstend", out)


class TestThermoLevels(_SolverTest):
    """In-solver test (runmode 1010, tests.f90:tests_thermo_levels): after thermodynamics,
    ql0(k) is the saturation adjustment of thl0(k), qt0(k) and ql0h(k) that of thl0h(k),
    qt0h(k), cell by cell (#403)."""

    def test_full_and_half_level_ql(self):
        d = tc.write_case(self.tmp / "thermolev", groups={"RUN": {"runmode": "1010"}})
        res = self.run_ok(d)
        out = res.stdout + res.stderr
        print("\n" + "\n".join(l for l in out.splitlines() if "worst" in l or "FAIL" in l))
        self.assertIn("ALL TESTS PASSED: tests_thermo_levels", out)


class TestMOKernels(_SolverTest):
    """In-solver test (runmode 1012, tests.f90:tests_mo_kernels) of the pure Monin-Obukhov kernels:
    psim/psih against Businger-Dyer (unstable) and Beljaars-Holtslag (stable) reference values,
    phim/phih, the calc_obl_iter round trip from the bulk Richardson number of a known L, and the
    per-section surface update of iwallmom = 4 (mo_local_update) at its fixed point: the neutral log
    law, the unstable u* with the flux-defined Obukhov length, prescribed u* and the u* floor."""

    def test_stability_functions_and_obl_iteration(self):
        d = tc.write_case(self.tmp / "mokernels", groups={"RUN": {"runmode": "1012"}})
        res = self.run_ok(d)
        out = res.stdout + res.stderr
        print("\n" + "\n".join(l for l in out.splitlines() if "worst" in l or "FAIL" in l))
        self.assertIn("ALL TESTS PASSED: tests_mo_kernels", out)


class TestRandomPerturbation(_SolverTest):
    """Initial perturbations (modstartup:randomize_field, #395): counter-based hash of
    (irandom, field tag, iglob, jglob, k), so independent of the MPI decomposition, with
    uniform statistics, no spectral peak and independent noise per field."""

    AMPL = {"randu": 0.5, "randthl": 0.3, "randqt": 1.0e-4}
    NAMES = ("u", "v", "w", "thl", "qt")

    def _dump0(self, name, npx=1, npy=1, **amp):
        rnd = dict(self.AMPL)
        rnd.update(amp)
        groups = {"RUN": {"runmode": "1011", "runtime": "0.5", "lrandomize": ".true.", "nprocx": str(npx),
                          "nprocy": str(npy), "randu": str(rnd["randu"]),
                          "randthl": str(rnd["randthl"]), "randqt": str(rnd["randqt"])},
                  "DOMAIN": {"itot": "32", "jtot": "32", "xlen": "1120.", "ylen": "1120."},
                  "OUTPUT": {"lfielddump": ".true.", "tfielddump": "10.",
                             "fieldvars": "'u0,v0,w0,th,qt'"}}
        d = tc.write_case(self.tmp / name, groups=groups)
        self.run_ok(d, npx * npy)
        # runmode 1011 writes the initial field as record 0 and stops
        tiles = {}
        for f in sorted(d.glob("fielddump.*.*.310.nc")):
            ix, iy = int(f.name.split(".")[1]), int(f.name.split(".")[2])
            with Dataset(f) as nc:
                tiles[ix, iy] = {v: np.asarray(nc[v][0]) for v in self.NAMES}   # (z, y, x)
        nx = 1 + max(k[0] for k in tiles)
        ny = 1 + max(k[1] for k in tiles)
        return {v: np.concatenate([np.concatenate([tiles[ix, iy][v] for ix in range(nx)], axis=2)
                                   for iy in range(ny)], axis=1) for v in self.NAMES}

    def test_decomposition_independent(self):
        ser = self._dump0("rnd_ser")
        par = self._dump0("rnd_par", 2, 2)
        for v in self.NAMES:
            with self.subTest(var=v):
                np.testing.assert_array_equal(par[v], ser[v])

    def test_statistics(self):
        f = self._dump0("rnd_stat")
        a = {"u": self.AMPL["randu"], "v": self.AMPL["randu"], "w": self.AMPL["randu"],
             "thl": self.AMPL["randthl"], "qt": self.AMPL["randqt"]}
        noise = {}
        for v in self.NAMES:
            x = f[v] - f[v].mean(axis=(1, 2), keepdims=True)    # remove the level mean (profile)
            if v == "w":
                x = x[1:]                                      # skip the w = 0 wall level
            noise[v] = x / a[v]
            with self.subTest(var=v):
                self.assertAlmostEqual(float(noise[v].mean()), 0.0, delta=0.01)
                self.assertAlmostEqual(float(noise[v].var()), 1.0 / 3.0, delta=0.01)
                # spectrum: power at every (kx, ky) of each level, averaged over levels
                p = (np.abs(np.fft.fft2(noise[v], axes=(1, 2))) ** 2).mean(axis=0)
                p[0, 0] = 0.0
                self.assertLess(p.max(), 3.0 * p[p > 0].mean())
        for i, v1 in enumerate(self.NAMES):
            for v2 in self.NAMES[i + 1:]:
                n = min(noise[v1].shape[0], noise[v2].shape[0])
                c = np.corrcoef(noise[v1][:n].ravel(), noise[v2][:n].ravel())[0, 1]
                with self.subTest(pair=(v1, v2)):
                    self.assertLess(abs(c), 0.02)

    def test_thl_qt_noise_off_by_default(self):
        f = self._dump0("rnd_off", randthl=0.0, randqt=0.0)
        for v in ("thl", "qt"):
            with self.subTest(var=v):
                self.assertEqual(float(np.ptp(f[v], axis=(1, 2)).max()), 0.0)
        self.assertGreater(float(np.ptp(f["u"], axis=(1, 2)).max()), 0.1)

    def test_thl_qt_noise_present(self):
        f = self._dump0("rnd_on")
        self.assertGreater(float(np.ptp(f["thl"], axis=(1, 2)).min()), 0.1)
        self.assertGreater(float(np.ptp(f["qt"], axis=(1, 2)).min()), 1.0e-4)
        self.assertGreaterEqual(float(f["qt"].min()), 0.0)


if __name__ == "__main__":
    unittest.main()
