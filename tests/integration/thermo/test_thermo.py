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
            "NAMTHERMOSTAT": {"lthermostat": ".true.", "tthermostat": "0.5"}},
            thl=np.full(64, 289.0), qt=np.full(64, 5.0e-3))
        self.run_ok(d)
        with Dataset(d / "thermostat_z.310.nc") as z:
            return float(z["thl"][-1, 0]) - 289.0, float(z["qt"][-1, 0]) - 5.0e-3

    def test_positive_flux_is_downward(self):
        # uDALES convention (D1): prescribed boundary fluxes are positive in -z,
        # so positive wtsurf/wqsurf cool/dry the first level by flux*dt/dz.
        dthl, dqt = self._first_level_change("sign", 1.0, 1.0e-3)
        dz, dt = 20.0, 1.0
        self.assertAlmostEqual(dthl, -1.0 * dt / dz, delta=2e-2 * 1.0 * dt / dz)
        self.assertAlmostEqual(dqt, -1.0e-3 * dt / dz, delta=2e-2 * 1.0e-3 * dt / dz)


class TestMO(_SolverTest):
    def _mo_case(self, name, bcbotm, extra_bc=None):
        # upward heat/moisture flux = negative in the uDALES (-z positive) convention
        bc = {"BCbotm": str(bcbotm), "wtsurf": "-0.05", "wqsurf": "-5.0e-5",
              "z0": "2.e-4", "z0h": "2.e-4", "thls": "289.0", "qts": "9.0e-3"}
        bc.update(extra_bc or {})
        # 20 fixed 1-s steps so the laminar (thls, obl) pair has converged: surface_mo
        # runs every RK substep and re-diagnoses thls each time.
        d = tc.write_case(self.tmp / name, u=7.0, v=-5.5, groups={
            "RUN": {"runtime": "19.5"}, "BC": bc,
            "NAMTHERMOSTAT": {"lthermostat": ".true.", "tthermostat": "1."}})
        self.run_ok(d)
        return d

    def test_mo_ustar_obl_match_reference(self):
        d = self._mo_case("mo4", 4)
        with Dataset(d / "thermostat_t.310.nc") as t, Dataset(d / "thermostat_z.310.nc") as z:
            ustar, obl = float(t["ustar"][-1]), float(t["obl"][-1])
            thls, qts = float(t["thls"][-1]), float(t["qts"][-1])
            thl1, qt1 = float(z["thl"][-1, 0]), float(z["qt"][-1, 0])
            obl_prev = float(t["obl"][-2])
            u1, v1 = float(t["u1"][-1]), float(t["v1"][-1])   # drag slows the wind slightly
        self.assertGreater(ustar, 0.01)
        self.assertLess(obl, 0.0)                         # upward buoyancy flux -> unstable
        self.assertGreater(thls, thl1)                    # surface warmer than air
        self.assertAlmostEqual(obl, obl_prev, delta=0.02 * abs(obl))   # converged
        # Reference: same iteration on the recorded converged state (float32 records).
        z1 = 10.0
        L = tc.obl_iter(thl1, qt1, thls, qts, z1, 2e-4, 2e-4, u1, v1, -0.01)
        self.assertAlmostEqual(obl, L, delta=0.02 * abs(L))
        cm = tc.FKAR ** 2 / (np.log(z1 / 2e-4) - tc.psim(z1 / L) + tc.psim(2e-4 / L)) ** 2
        self.assertAlmostEqual(ustar, np.sqrt(cm) * np.hypot(u1, v1), delta=0.02 * ustar)

    def test_prescribed_ustar(self):
        d = self._mo_case("mo5", 5, {"ustin": "0.28"})
        with Dataset(d / "thermostat_t.310.nc") as t:
            self.assertAlmostEqual(float(t["ustar"][-1]), 0.28, places=6)   # float32 record


RV_RD = 461.5 / 287.04
RLV, CP = 2.26e6, 1004.0   # check modglobal: rlv, cp


def dthvdz_surface(thl, qt, ql, exn, dthldz, dqtdz):
    """DALES dev 075ac4e calthv, surface level."""
    if ql > 0:
        temp = thl * exn + RLV / CP * ql
        qs = qt - ql
        a = (1 - qt + RV_RD * qs * (1 + RLV / (461.5 * temp))) / (1 + RLV ** 2 * qs / (CP * 461.5 * temp ** 2))
        b = a * RLV / (temp * CP) - 1
    else:
        a = 1 + (RV_RD - 1) * qt
        b = RV_RD - 1
    return a * dthldz + b * thl * dqtdz


class TestCalthv(_SolverTest):
    def test_surface_level_unsaturated(self):
        d = tc.write_case(self.tmp / "calthv_sfc", u=7.0, v=-5.5, groups={
            "RUN": {"runtime": "2.5"},
            "OUTPUT": {"lfielddump": ".true.", "tfielddump": "1.", "fieldvars": "'th,qt,ql,dt'"},
            "BC": {"BCbotm": "4", "wtsurf": "-0.05", "wqsurf": "-5.0e-5", "z0": "2.e-4", "z0h": "2.e-4"},
            "NAMTHERMOSTAT": {"lthermostat": ".true.", "tthermostat": "1."}})
        self.run_ok(d)
        # Record alignment (approximate). Per RK substep the order is: bottom (-> surface_mo)
        # ... tstep_integrate -> fielddump/thermostat (write only at rk3step == 3 when due)
        # -> boundary -> thermodynamics; calthv runs on EVERY substep. Dump record n holds
        # end-of-step-n thl/qt. The dthvdz in dump record n+1 was computed by calthv at the
        # end of RK substep 2 of step n+1, from the substep-2 intermediate fields and the
        # ustar_mo/obl_mo of that substep. That is closer to thermostat record n+1 (substep 3
        # of step n+1) than to record n, so ustar/obl are taken from record n+1 and thl/qt
        # from dump record n (ql from n+1, diagnosed from them). The pairing is off by one
        # RK substep in both state and ustar/obl; it works because the laminar column is
        # near-steady and the 5e-3 tolerance absorbs it.
        n = 0
        with Dataset(d / "fielddump.000.000.310.nc") as f, \
             Dataset(d / "thermostat_t.310.nc") as t, Dataset(d / "thermostat_z.310.nc") as z:
            dthv = np.asarray(f["dthvdz"][n + 1])        # (z, y, x)
            thl1 = np.asarray(f["thl"][n])[0].mean()
            qt1 = np.asarray(f["qt"][n])[0].mean()
            ql1 = np.asarray(f["ql"][n + 1])[0].mean()
            exn1 = float(z["exnf"][n + 1, 0])
            ustar, obl = float(t["ustar"][n + 1]), float(t["obl"][n + 1])
        z1, kap = 10.0, tc.FKAR
        phih = (1 - 16 * z1 / obl) ** -0.5 if obl < 0 else min(1 + 5 * z1 / obl, 6.0)
        dthldz = -0.05 / ustar * phih / (kap * z1)        # upward fluxes: +0.05, +5e-5
        dqtdz = -5.0e-5 / ustar * phih / (kap * z1)
        ref = dthvdz_surface(thl1, qt1, ql1, exn1, dthldz, dqtdz)
        self.assertLess(ref, 0.0)                          # unstable surface layer
        self.assertAlmostEqual(float(dthv[0].mean()), ref, delta=5e-3 * abs(ref))

    def test_no_mo_keeps_surface_level(self):
        d = tc.write_case(self.tmp / "calthv_wf", groups={
            "RUN": {"runtime": "0.5"},
            "OUTPUT": {"lfielddump": ".true.", "tfielddump": "0.5", "fieldvars": "'dt'"}})
        self.run_ok(d)
        with Dataset(d / "fielddump.000.000.310.nc") as f:
            self.assertTrue(np.all(np.abs(np.asarray(f["dthvdz"][-1])[0]) <= 2e-10))  # float32(eps1)


class TestThermostat(_SolverTest):
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
            "NAMTHERMOSTAT": {"lthermostat": ".true.", "tthermostat": "1."}})
        self.run_ok(d)
        with Dataset(d / "thermostat_t.310.nc") as t, Dataset(d / "thermostat_z.310.nc") as z:
            self.assertGreaterEqual(len(t["time"]), 2)
            ql = z["ql"][-1, :]
            rhof = z["rhof"][-1, :]
            dz = 20.0
            lwp_np = float(np.sum(rhof * ql * dz))
            self.assertAlmostEqual(float(t["lwp_bar"][-1]), lwp_np, delta=1e-6 * max(lwp_np, 1e-12))
            self.assertGreater(lwp_np, 0.0)             # RF01 mixed layer is saturated aloft
            self.assertEqual(float(t["cfrac"][-1]), 1.0)
            self.assertAlmostEqual(float(t["zi_qt8"][-1]), 840.0, delta=20.0)
            self.assertAlmostEqual(float(t["zi"][-1]), 840.0, delta=20.0)
            # placeholder until the radiation task (Task 5) fills them
            for name in ("thlprad", "lwu", "lwd"):
                self.assertTrue(np.all(z[name][-1, :] == 0.0))

    def test_dry_column_has_no_moist_diagnostics(self):
        d = tc.write_case(self.tmp / "drystat", groups={
            "RUN": {"runtime": "1.5"},
            "PHYSICS": {"lmoist": ".false."}, "BC": {"qts": "-1."},
            "NAMTHERMOSTAT": {"lthermostat": ".true.", "tthermostat": "1."}})
        self.run_ok(d)
        with Dataset(d / "thermostat_t.310.nc") as t:
            self.assertGreaterEqual(len(t["time"]), 1)
            self.assertEqual(float(t["zi_qt8"][-1]), 0.0)
            self.assertEqual(float(t["lwp_bar"][-1]), 0.0)

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
                  "NAMTHERMOSTAT": {"lthermostat": ".true.", "tthermostat": "0.5"}}
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
        with Dataset(d / "thermostat_z.310.nc") as z:
            ql, rho, exn = (np.asarray(z[v][-1]) for v in ("ql", "rhof", "exnf"))
            got = np.asarray(z["thlprad"][-1])
        ref = radpar_reference(ql, rho, exn, 20.0, 70.0, 22.0, 85.0)
        self._report("liquid", got, ref)
        # thlprad uses start-of-step ql; the record's ql is end-of-step (~0.3 % change at
        # cloud top). A sign or face off-by-one error is O(1), so 1e-2 still discriminates.
        np.testing.assert_allclose(got, ref, rtol=1e-2, atol=1e-2 * np.abs(ref).max())
        self.assertLess(got.min(), -1e-4)          # cloud-top cooling present

    def test_smoke_longwave_matches_formula(self):
        rad = dict(self.RAD, rad_smoke=".true.", isvsmoke="1", rka="0.02", dlwtop="60.", dlwbot="0.")
        zf = (np.arange(64) + 0.5) * 20.0
        smoke = np.where(zf < 700, 1.0, 0.0)
        d = self._run("rad_smoke", rad, groups={"PHYSICS": {"lmoist": ".false."}, "BC": {"qts": "-1."},
                                                "SCALARS": {"nsv": "1"}},
                      thl=np.where(zf < 700, 288.0, 295.0), qt=np.zeros(64), sv=smoke[:, None])
        with Dataset(d / "thermostat_z.310.nc") as z:
            rho, exn = np.asarray(z["rhof"][-1]), np.asarray(z["exnf"][-1])
            got = np.asarray(z["thlprad"][-1])
        ref = radpar_reference(smoke, rho, exn, 20.0, 60.0, 0.0, 0.02)
        self._report("smoke", got, ref)
        np.testing.assert_allclose(got, ref, rtol=1e-2, atol=1e-2 * np.abs(ref).max())

    def test_dycoms_term_matches_formula(self):
        zf = (np.arange(64) + 0.5) * 20.0
        D = 3.75e-6
        rad = dict(self.RAD, lrad_dycoms=".true.", rad_dycoms_D=str(D))
        d = self._run("ft_formula", rad)
        with Dataset(d / "thermostat_z.310.nc") as z:
            ql, qt, rho, exn = (np.asarray(z[v][-1], dtype=float) for v in ("ql", "qt", "rhof", "exnf"))
            got = np.asarray(z["thlprad"][-1])
        # zi by the solver's rule: highest level with qt > 8 g/kg, linear to the next level.
        # (Laminar RF01: 832.7 m, between the 830/850 m cells, not 840 m.)
        k = np.where(qt > 8e-3)[0].max()
        zi = zf[k] + (8e-3 - qt[k]) * (zf[k + 1] - zf[k]) / (qt[k + 1] - qt[k])
        ref = radpar_reference(ql, rho, exn, 20.0, 70.0, 22.0, 85.0)
        ft = zf > zi
        ref[ft] -= D / 3 * (np.cbrt(zf[ft] - zi) + zi * (zf[ft] - zi) ** (-2 / 3))   # cooling
        self._report("dycoms", got, ref)
        np.testing.assert_allclose(got, ref, rtol=1e-2, atol=1e-2 * np.abs(ref).max())
        self.assertTrue(np.all(got[ft] < 0))

    def test_dycoms_term_balances_subsidence(self):
        # 300 s so the free-troposphere subsidence signal (~2e-2 K) is far above
        # float32 resolution (3e-5 K at 297 K).
        zf = (np.arange(64) + 0.5) * 20.0
        D = 3.75e-6
        rad = dict(self.RAD, lrad_dycoms=".true.", rad_dycoms_D=str(D))
        common = dict(wfls=-D * zf, groups={"RUN": {"runtime": "299.5"},
                                            "NAMTHERMOSTAT": {"lthermostat": ".true.", "tthermostat": "300."}})
        a = self._run("ft_both", rad, **common)
        b = self._run("ft_subs", {"iradiation": "0"}, **common)
        with Dataset(a / "thermostat_z.310.nc") as za, Dataset(b / "thermostat_z.310.nc") as zb:
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


class TestDecomposition(_SolverTest):
    """Radiation + MO surface + moist SGS give the same answer on 1, 2 and 4 ranks."""

    def test_splits_agree(self):
        groups = {"RUN": {"runtime": "2.5", "lrandomize": ".true."},
                  "BC": {"BCbotm": "4", "wtsurf": "-0.012", "wqsurf": "-4.7e-5",
                         "z0": "2.e-4", "z0h": "2.e-4"},
                  "NAMRADIATION": dict(TestRadiation.RAD, lrad_dycoms=".true."),
                  "NAMTHERMOSTAT": {"lthermostat": ".true.", "tthermostat": "1."}}
        zf = (np.arange(64) + 0.5) * 20.0
        out = {}
        for label, npx, npy in (("serial", 1, 1), ("x", 2, 1), ("xy", 2, 2)):
            g = {k: dict(v) for k, v in groups.items()}
            g["RUN"].update(nprocx=str(npx), nprocy=str(npy))
            d = tc.write_case(self.tmp / f"dec_{label}", groups=g, wfls=-3.75e-6 * zf)
            self.run_ok(d, npx * npy)
            with Dataset(d / "thermostat_t.310.nc") as t, Dataset(d / "thermostat_z.310.nc") as z:
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


if __name__ == "__main__":
    unittest.main()
