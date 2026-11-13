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


if __name__ == "__main__":
    unittest.main()
