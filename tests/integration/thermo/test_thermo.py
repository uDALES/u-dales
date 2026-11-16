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
        # neutral cap; its sign is that of a round-off Rib (FMA contraction at -O3 -march=native)
        self.assertEqual(abs(r["obl"][-1]), 1.e6)
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

    def test_horizontal_facet_cells_get_surface_layer_gradient(self):
        """The cells above a horizontal Monin-Obukhov facet take the surface-layer dthvdz
        (DALES calthv, first level) and the first-level one-equation TKE sources (whose
        formula tests.f90:tests_sgs_oneeqn (c) checks): dthvdz(kb) = a dthl/dz + b thl dqt/dz with
        dthl/dz = -w'thl' phih(z1/L) / (u* kappa z1). The c-grid state is computed from the
        cell-centre wind, the bldump u*, L from the u/v sections: identical in a homogeneous
        column up to the drift within the last step."""
        d, out = self._mo_case("mo_cells", extra={"OUTPUT": {"lfielddump": ".true.", "tfielddump": "19.5",
                                                               "fieldvars": "'dt,th,qt'"}})
        self.assertIn("iwallmom = 4: 256 cells above horizontal Monin-Obukhov facets", out)
        r = self._last(d)
        ustar, obl = r["ustar"][-1], r["obl"][-1]
        with Dataset(d / "fielddump.000.000.310.nc") as f:
            dthvdz = np.asarray(f["dthvdz"][-1], dtype=float)[0]
            thl = float(np.asarray(f["thl"][-1])[0].mean())
            qt = float(np.asarray(f["qt"][-1])[0].mean())
        zeta = self.Z1 / obl
        phih = (1 - 16 * zeta) ** -0.5
        scale = phih / (ustar * tc.FKAR * self.Z1)
        a, b = 1 + (tc.RV / tc.RD - 1) * qt, tc.RV / tc.RD - 1
        expect = a * (-0.05 * scale) + b * thl * (-5.0e-5 * scale)
        self.assertLess(expect, 0.)
        np.testing.assert_allclose(dthvdz, expect, rtol=1e-3)

    def test_sloped_facet_cells_keep_wall_treatment(self):
        """A 20 deg Monin-Obukhov facet (west half horizontal, east half 20 deg) gets no
        first-level treatment: dthvdz(kb) stays 0 there (and so do the kb TKE sources: the
        same cell list gates both), the horizontal half takes the surface-layer value."""
        t20 = np.radians(20.0)
        normals = [(0.0, 0.0, 1.0), (np.sin(t20), 0.0, np.cos(t20))]
        d = tc.write_case(self.tmp / "mo_cells_slope", u=7.0, v=-5.5,
                          floor={"z0": self.Z0, "z0h": self.Z0, "normals": normals,
                                 "facet_of": lambda i, j: 1 if i <= 8 else 2},
                          groups={"RUN": {"runtime": "1.5"}, "BC": {"bctfz": "-0.05", "bcqfz": "-5.0e-5"},
                                  "OUTPUT": {"lfielddump": ".true.", "tfielddump": "1.5", "fieldvars": "'dt'"}})
        res = self.run_ok(d)
        out = res.stdout + res.stderr
        self.assertIn("iwallmom = 4: 2 facets use Monin-Obukhov", out)
        self.assertIn("iwallmom = 4: 128 cells above horizontal Monin-Obukhov facets", out)
        with Dataset(d / "fielddump.000.000.310.nc") as f:
            dthvdz = np.asarray(f["dthvdz"][-1], dtype=float)[0]      # (y, x) at kb
        self.assertTrue(np.all(dthvdz[:, :8] < -1e-4))
        self.assertTrue(np.all(np.abs(dthvdz[:, 8:]) <= 2e-10))      # float32(eps1) floor of calthv



class TestMOFacetTemperature(_SolverTest):
    """Monin-Obukhov facets with a facet temperature (iwallmom = 4, iwalltemp = 2; the DALES
    isurf = 2 analogue): laminar homogeneous column over a flat facet floor whose temperature is
    fixed (lEB off), facet 2 K warmer than the air. Records every 1-s step.

    The reference air state is the mean of the last two records: the fluxes, u* and L recorded
    at the end of a step come from the state of its last RK3 substep, which starts at t + dt/2,
    so the record mean is that state to second order (the first-order drift is ~1e-3 per step
    here: heating Cs U dt / dzf of dT, deceleration u*^2 dt / dzf; against the end-of-step
    record the flux differs by 5.6e-4). L is re-iterated to convergence on that state in every
    update. Residuals: float32 records (3e-5 K ulp at 289 K on a 2 K dT: 1.5e-5 per record),
    the 1e-4 stopping tolerance of the L iteration (Cs and u* depend on L only through
    psi(z/L), |dCs/Cs| < 0.1 |dL/L| here), the second-order drift (~1e-6). Measured: 5e-6.
    Bound: 1e-4."""

    Z0, Z0H, Z1, TFAC = 2.e-4, 2.e-5, 10.0, 291.0
    RTOL = 1e-4

    def _case(self, name, extra=None, runtime="9.5", lgr=True, nprocs=1, tfac=None, vwind=-5.5, **kw):
        groups = {"RUN": {"runtime": runtime}, "BC": {"bctfz": "0.", "bcqfz": "0."},
                  "WALLS": {"iwalltemp": "2"},
                  "OUTPUT": {"lbldump": ".true.", "tstatsdump": "1."}}
        for g, kv in (extra or {}).items():
            groups.setdefault(g, {}).update(kv)
        floor = {"z0": self.Z0, "z0h": self.Z0H, "tfacinit": self.TFAC if tfac is None else tfac, "lgr": lgr}
        floor.update(kw.pop("floor", {}))
        d = tc.write_case(self.tmp / name, u=7.0, v=vwind, floor=floor, groups=groups, **kw)
        res = self.run_ok(d, nprocs)
        return d, res.stdout + res.stderr

    @staticmethod
    def _records(d):
        with Dataset(d / "bldump_t.310.nc") as t, Dataset(d / "bldump_z.310.nc") as z:
            r = {v: np.asarray(t[v][:], dtype=float) for v in ("ustar", "obl", "thls", "qts", "u1", "v1", "wtheta", "wq")}
            r["thl1"] = np.asarray(z["thl"][:, 0], dtype=float)
            r["qt1"] = np.asarray(z["qt"][:, 0], dtype=float)
        return r

    def _mid(self, r):
        """the outputs of the last step and the air state of its last substep (record mean)"""
        m = {k: v[-1] for k, v in r.items()}
        m.update({k: 0.5 * (r[k][-1] + r[k][-2]) for k in ("u1", "v1", "thl1", "qt1")})
        m["U"] = np.hypot(m["u1"], m["v1"])
        return m

    def test_heat_flux_is_cs_u_dT(self):
        d, out = self._case("mot_heat")
        self.assertIn("iwallmom = 4: 1 facets use Monin-Obukhov (slope <=  30.0 deg), 0 use the stability-corrected wall function", out)
        self.assertIn("iwallmom = 4: 256 cells above horizontal Monin-Obukhov facets", out)
        r = self._records(d)
        np.testing.assert_allclose(r["thls"], self.TFAC, rtol=1e-7)     # the facet temperature, fixed without lEB
        m = self._mid(r)
        # bcqfz = 0 on the green floor: no moisture flux, qskin = qt of the air
        self.assertEqual(np.max(np.abs(r["wq"])), 0.)
        L = tc.obl_iter(m["thl1"], m["qt1"], self.TFAC, m["qt1"], self.Z1, self.Z0, self.Z0H, m["u1"], m["v1"], -0.1)
        cm, cs = tc.drag_coefficients(self.Z1, L, self.Z0, self.Z0H)
        wtheta = cs * m["U"] * (self.TFAC - m["thl1"])
        print(f"\n[mo tskin] L {m['obl']:.5g} ref {L:.5g}; wtheta {m['wtheta']:.6g} ref {wtheta:.6g}; "
              f"u* {m['ustar']:.6g} ref {np.sqrt(cm) * m['U']:.6g}")
        self.assertLess(L, 0.)                                          # warm facet: unstable
        self.assertGreater(wtheta, 0.)                                  # upward
        np.testing.assert_allclose(m["obl"], L, rtol=self.RTOL)
        np.testing.assert_allclose(m["wtheta"], wtheta, rtol=self.RTOL)
        np.testing.assert_allclose(m["ustar"], np.sqrt(cm) * m["U"], rtol=self.RTOL)
        self.assertGreater(r["thl1"][-1], r["thl1"][0])                 # the air is heated

    def test_green_facet_resistance_law(self):
        """iwallmoist = 2 on a green Monin-Obukhov facet: w'qt' = -moist_flux(resa = 1/(Cs U)) with
        the facet saturation humidity qsat(Tfacet) and hurel = 1 (wsoil = wfc). The qskin of the
        Richardson number is diagnosed from the previous flux: the recorded qts (u/v sections)."""
        d, _ = self._case("mot_moist", extra={"WALLS": {"iwallmoist": "2"}, "ENERGYBALANCE": {"wsoil": "313.", "wfc": "313."}})
        r = self._records(d)
        m = self._mid(r)
        L = tc.obl_iter(m["thl1"], m["qt1"], self.TFAC, m["qts"], self.Z1, self.Z0, self.Z0H, m["u1"], m["v1"], -0.1)
        _, cs = tc.drag_coefficients(self.Z1, L, self.Z0, self.Z0H)
        resa = 1.0 / (cs * m["U"])
        wq = -tc.moist_flux(resa, m["qt1"], tc.qsat_facet(self.TFAC), 1.0)
        print(f"\n[mo tskin moist] wq {m['wq']:.6g} ref {wq:.6g}; qts {m['qts']:.6g} = qt + wq/(Cs U) "
              f"{m['qt1'] + wq * resa:.6g}")
        self.assertGreater(wq, 0.)                                      # evaporation
        np.testing.assert_allclose(m["wq"], wq, rtol=self.RTOL)
        np.testing.assert_allclose(m["qts"], m["qt1"] + wq * resa, rtol=self.RTOL)
        np.testing.assert_allclose(m["wtheta"], cs * m["U"] * (self.TFAC - m["thl1"]), rtol=self.RTOL)

    def test_steep_facet_uses_stability_corrected_wall_function(self):
        """A floor facet declared at 60 deg (steeper than mo_maxslope) under iwallmom = 4,
        iwalltemp = 2 runs the iwallmom = 2 wall function (wfuno with stability): the fields equal
        those of an iwallmom = 2 run bit for bit, and differ from the neutral iwallmom = 3 run."""
        t60 = np.radians(60.0)
        fields = {}
        for name, iwm in (("steep_mo", "4"), ("steep_wfuno", "2"), ("steep_neutral", "3")):
            d, out = self._case(name, runtime="2.5", floor={"normals": [(np.sin(t60), 0.0, np.cos(t60))]},
                                extra={"WALLS": {"iwallmom": iwm},
                                       "OUTPUT": {"lfielddump": ".true.", "tfielddump": "2.5", "fieldvars": "'u0,v0,w0,th'"}})
            if iwm == "4":
                self.assertIn("0 facets use Monin-Obukhov (slope <=  30.0 deg), 1 use the stability-corrected wall function", out)
            with Dataset(d / "fielddump.000.000.310.nc") as f:
                fields[name] = {v: np.asarray(f[v][-1], dtype=float) for v in ("u", "v", "w", "thl")}
        for v in ("u", "v", "w", "thl"):
            np.testing.assert_array_equal(fields["steep_mo"][v], fields["steep_wfuno"][v], err_msg=v)
        self.assertGreater(np.max(np.abs(fields["steep_mo"]["u"] - fields["steep_neutral"]["u"])), 1e-6)

    def test_slope_split(self):
        """West half 20 deg (Monin-Obukhov), east half 60 deg (stability-corrected wall function).
        The bldump u* is the mean over the Monin-Obukhov sections only: the MO u* of facet 1 at
        its tangential wind, with L from the facet temperature, against the initial state. The one
        recorded step takes its u* and flux from the last RK3 substep, which starts at t + dt/2, so the
        reference is off by half a step of drift: deceleration u*^2 dt / (dzf |U|) ~ 8e-4 and heating
        Cs U dt / dzf ~ 1e-3 of dT per step, so <= 5e-4 (measured 1.9e-4 for u*, 3.7e-4 for the flux).
        Bound 1e-3."""
        t20, t60 = np.radians(20.0), np.radians(60.0)
        normals = [(np.sin(t20), 0.0, np.cos(t20)), (np.sin(t60), 0.0, np.cos(t60))]
        d, out = self._case("mot_split", runtime="0.5", vwind=0.0,
                            floor={"normals": normals, "facet_of": lambda i, j: 1 if i <= 8 else 2},
                            extra={"OUTPUT": {"tstatsdump": "0.5"}})
        self.assertIn("iwallmom = 4: 1 facets use Monin-Obukhov (slope <=  30.0 deg), 1 use the stability-corrected wall function", out)
        with Dataset(d / "bldump_t.310.nc") as t:
            ustar, wtheta = float(t["ustar"][-1]), float(t["wtheta"][-1])
        n1 = np.array(normals[0])
        uvec = np.array([7.0, 0.0, 0.0])
        utan = np.linalg.norm(uvec - uvec.dot(n1) * n1)
        thl1 = 289.0
        L = tc.obl_iter(thl1, 9.0e-3, self.TFAC, 9.0e-3, self.Z1, self.Z0, self.Z0H, utan, 0.0, -0.1)
        cm, cs = tc.drag_coefficients(self.Z1, L, self.Z0, self.Z0H)
        print(f"\n[mo split] u* {ustar / (np.sqrt(cm) * utan) - 1:.2e} wtheta {wtheta / (cs * utan * (self.TFAC - thl1)) - 1:.2e}")
        self.assertAlmostEqual(ustar, np.sqrt(cm) * utan, delta=1e-3 * ustar)
        self.assertAlmostEqual(wtheta, cs * utan * (self.TFAC - thl1), delta=1e-3 * wtheta)

    def test_energy_balance_receives_the_flux(self):
        """lEB = .true. (dtEB = dt = 1 s, no radiation: netsw = 0, svf = 1, skyLW = 0): the facet
        sensible heat flux of the energy balance (facEB hf, W/m2, positive into the facet) is the
        Monin-Obukhov flux of the last substep, -rhoa cp wtheta, and the next step's Monin-Obukhov
        surface temperature is the facet temperature the energy balance computed (facT layer 1)."""
        d = self.tmp / "mot_eb"
        d.mkdir()
        for name, txt in (("svf", "1.0"), ("netsw", "0.0")):
            (d / f"{name}.inp.310").write_text(f"# {name}\n{txt}\n")
        (d / "vfsparse.inp.310").write_text("")
        d, _ = self._case("mot_eb", runtime="4.5", lgr=False,
                          extra={"ENERGYBALANCE": {"lEB": ".true.", "lwriteEBfiles": ".true.", "dtEB": "1.",
                                                   "lvfsparse": ".true.", "nnz": "0", "skyLW": "0."}})
        r = self._records(d)
        with Dataset(d / "facEB.310.nc") as f:
            hf = np.asarray(f["hf"][:], dtype=float).reshape(-1)
        with Dataset(d / "facT.310.nc") as f:
            T1 = np.asarray(f["T"][:], dtype=float).reshape(len(hf), -1)[:, 0]
        print(f"\n[mo EB] hf {hf}; -rhoa cp wtheta {-1.2 * 1004. * r['wtheta']}; facT1 {T1}; thls {r['thls']}")
        n = min(len(hf), len(r["wtheta"]))
        self.assertGreaterEqual(n, 3)
        np.testing.assert_allclose(hf[:n], -1.2 * CP * r["wtheta"][:n], rtol=1e-5)
        self.assertTrue(np.all(np.diff(T1) < 0.))                       # the facet loses heat (H up, LW out)
        np.testing.assert_allclose(r["thls"][1:n], T1[:n - 1], rtol=1e-7)   # MO uses the EB facet temperature

    def test_energy_balance_receives_the_latent_flux(self):
        """lEB with a green floor and iwallmoist = 2 (wsoil = wfc): the facet latent heat flux of the
        energy balance (facEB ef, W/m2, positive into the facet; modEB:updateGR converts with rhoa rlv)
        is the Monin-Obukhov resistance-law flux of the last substep, -rhoa rlv wq (float32: rtol 1e-5)."""
        d = self.tmp / "mot_eb_moist"
        d.mkdir()
        for name, txt in (("svf", "1.0"), ("netsw", "0.0")):
            (d / f"{name}.inp.310").write_text(f"# {name}\n{txt}\n")
        (d / "vfsparse.inp.310").write_text("")
        d, _ = self._case("mot_eb_moist", runtime="4.5", lgr=True,
                          extra={"WALLS": {"iwallmoist": "2"},
                                 "ENERGYBALANCE": {"lEB": ".true.", "lwriteEBfiles": ".true.", "dtEB": "1.",
                                                   "lvfsparse": ".true.", "nnz": "0", "skyLW": "0.",
                                                   "wsoil": "313.", "wfc": "313."}})
        r = self._records(d)
        with Dataset(d / "facEB.310.nc") as f:
            ef = np.asarray(f["ef"][:], dtype=float).reshape(-1)
        n = min(len(ef), len(r["wq"]))
        print(f"\n[mo EB moist] ef {ef[:n]}; -rhoa rlv wq {-1.2 * RLV * r['wq'][:n]}")
        self.assertGreaterEqual(n, 3)
        self.assertTrue(np.all(r["wq"][:n] > 0.))                      # evaporation
        np.testing.assert_allclose(ef[:n], -1.2 * RLV * r["wq"][:n], rtol=1e-5)

    def test_first_level_cells_get_surface_layer_gradient(self):
        """iwalltemp = 2 twin of TestMO.test_horizontal_facet_cells_get_surface_layer_gradient: the cells
        above the horizontal floor take dthvdz(kb) = (1 + 0.61 qt) dthl/dz with dthl/dz =
        -wtheta phih(z1/L) / (u* kappa z1) from the Monin-Obukhov state of the facet temperature (no
        moisture flux: bcqfz = 0). wtheta comes from the c-sections that also set the gradient; u* and
        L from the u/v sections: identical in a homogeneous column up to the drift within the last
        step and the float32 records, so the TestMO bound 1e-3 applies. Unstable and stable facets."""
        for name, tfac in (("mot_cells_unstable", self.TFAC), ("mot_cells_stable", 287.0)):
            with self.subTest(case=name):
                d, out = self._case(name, tfac=tfac, extra={"OUTPUT": {"lfielddump": ".true.", "tfielddump": "9.5",
                                                                        "fieldvars": "'dt,th,qt'"}})
                self.assertIn("iwallmom = 4: 256 cells above horizontal Monin-Obukhov facets", out)
                r = self._records(d)
                ustar, obl, wtheta = r["ustar"][-1], r["obl"][-1], r["wtheta"][-1]
                with Dataset(d / "fielddump.000.000.310.nc") as f:
                    dthvdz = np.asarray(f["dthvdz"][-1], dtype=float)[0]
                    qt = float(np.asarray(f["qt"][-1])[0].mean())
                zeta = self.Z1 / obl
                phih = (1 - 16 * zeta) ** -0.5 if zeta < 0 else (1 + 5 * zeta if zeta < 1 else 6.0)
                expect = (1 + (tc.RV / tc.RD - 1) * qt) * (-wtheta * phih / (ustar * tc.FKAR * self.Z1))
                print(f"\n[mo cells {name}] L {obl:.4g} wtheta {wtheta:.4g} dthvdz {dthvdz.mean():.6g} ref {expect:.6g}")
                self.assertEqual(np.sign(expect), -np.sign(tfac - 289.0))     # unstable: dthv/dz < 0
                np.testing.assert_allclose(dthvdz, expect, rtol=1e-3)

    def test_stable_facet(self):
        """Facet 2 K colder than the air (stable, L > 0, downward flux): L, wtheta = Cs U dT and
        u* = sqrt(Cm) U on the last-substep state, the bound of the unstable case (class docstring)."""
        tfac = 287.0
        d, _ = self._case("mot_stable", tfac=tfac)
        r = self._records(d)
        np.testing.assert_allclose(r["thls"], tfac, rtol=1e-7)
        m = self._mid(r)
        L = tc.obl_iter(m["thl1"], m["qt1"], tfac, m["qt1"], self.Z1, self.Z0, self.Z0H, m["u1"], m["v1"], -0.1)
        cm, cs = tc.drag_coefficients(self.Z1, L, self.Z0, self.Z0H)
        wtheta = cs * m["U"] * (tfac - m["thl1"])
        print(f"\n[mo tskin stable] L {m['obl']:.5g} ref {L:.5g}; wtheta {m['wtheta']:.6g} ref {wtheta:.6g}; "
              f"u* {m['ustar']:.6g} ref {np.sqrt(cm) * m['U']:.6g}")
        self.assertGreater(L, 0.)
        self.assertLess(m["wtheta"], 0.)
        np.testing.assert_allclose(m["obl"], L, rtol=self.RTOL)
        np.testing.assert_allclose(m["wtheta"], wtheta, rtol=self.RTOL)
        np.testing.assert_allclose(m["ustar"], np.sqrt(cm) * m["U"], rtol=self.RTOL)

    def test_decomposition_on_ridge(self):
        """1 vs 2 ranks on the sloped ridge (TestMOSlopedFacet geometry) with the facet temperature:
        the c-grid sections of sloped Monin-Obukhov facets carry their own state (mo_section_tskin in
        wallfunheat); only the order of the area-weighted bldump sums differs between the splits."""
        if not tc.ibm_preprocessor_available():
            self.skipTest("no gfortran for the udprep legacy IBM preprocessor")
        try:
            prep = tc.prepare_ridge(self.tmp / "prep", 20.0, self.Z0)
        except ImportError as exc:
            self.skipTest(f"udprep dependencies missing: {exc}")
        out = {}
        for npx in (1, 2):
            dd = tc.write_ridge_run(prep, self.tmp / f"ridge_t2_{npx}", u=7.0, v=0.0, thl=np.full(64, 289.0),
                                    qt=np.full(64, 5.0e-3), groups={
                "RUN": {"runtime": "3.0", "nprocx": str(npx)}, "WALLS": {"iwalltemp": "2"},
                "BC": {"bctfz": "0.", "bcqfz": "0."},
                "OUTPUT": {"lbldump": ".true.", "tstatsdump": "1.", "lfielddump": ".true.",
                           "tfielddump": "3.0", "fieldvars": "'u0,th'"}})
            nf = len((dd / "facets.inp.310").read_text().splitlines()) - 1
            (dd / "Tfacinit.inp.310").write_text("# initial facet temperature [K]\n" + f"{self.TFAC:.6f}\n" * nf)
            res = self.run_ok(dd, nprocs=npx)
            self.assertIn("iwallmom = 4: 4 facets use Monin-Obukhov", res.stdout + res.stderr)
            with Dataset(dd / "bldump_t.310.nc") as t:
                out[npx] = {v: np.asarray(t[v][:], dtype=float) for v in ("ustar", "obl", "wtheta")}
            if npx == 1:
                with Dataset(dd / "fielddump.000.000.310.nc") as f:
                    thl = np.asarray(f["thl"][-1], dtype=float)
        for v in ("ustar", "obl", "wtheta"):
            self.assertTrue(np.all(np.isfinite(out[1][v])))
            np.testing.assert_allclose(out[2][v], out[1][v], rtol=1e-6, err_msg=v)
        self.assertTrue(np.all(out[1]["wtheta"] > 0.))                  # warm facets heat the air
        self.assertGreater(thl.max(), 289.0)

    def test_offset_floor_reconstruction_point(self):
        """A horizontal STL floor that is not on a grid face (z = 9 m, first cell centre at 10 m) with
        z0 = 0.5 m: log(1 m / z0) <= 1, so its c-grid sections use reconstruction points (modibm
        initibmwallfun: lcomprec off). initibm_mocells leaves such sections out of the first-level
        treatment (it requires log(bnddst / max(z0, z0h)) > 1, which implies lcomprec), so they do
        not reuse the cell-centre state of mo_cells_update: wallfunheat updates them from the
        reconstruction-point air state and divides by the same |utan|, so fac_htc = Cs max(|utan|,
        0.1) / max(|utan|, 0.1) = fac_cth (Task 24 M5). Heat budget: the fluid gains, per step, dt A
        wtheta of the last substep (the RK3 step increment is dt times the last-substep tendency);
        residual: float32 thl (2e-5 summed over the first level) and the IBM corrections; measured
        4e-5, bound 5e-4. 1 and 2 ranks agree (order of the bldump area sums only)."""
        if not tc.ibm_preprocessor_available():
            self.skipTest("no gfortran for the udprep legacy IBM preprocessor")
        z0, z0h, zp, tfac = 0.5, 0.05, 9.0, 299.0
        try:
            prep = tc.prepare_stl(self.tmp / "prep_plane", lambda p: tc.write_plane_stl(p, zp, 560., 560.), z0, z0h)
        except ImportError as exc:
            self.skipTest(f"udprep dependencies missing: {exc}")
        dist = np.loadtxt(prep / "facet_sections_c.txt", skiprows=1)[:, 3]
        np.testing.assert_allclose(dist, 10.0 - zp)
        self.assertLessEqual(np.log((10.0 - zp) / z0), 1.0)              # reconstruction-point sections
        out = {}
        for npx in (1, 2):
            d = tc.write_ridge_run(prep, self.tmp / f"plane_{npx}", u=7.0, v=0.0, thl=np.full(64, 289.0),
                                   qt=np.full(64, 5.0e-3), groups={
                "RUN": {"runtime": "3.0", "nprocx": str(npx)},
                "WALLS": {"iwalltemp": "2", "lwritefac": ".true.", "dtfac": "1."},
                "BC": {"bctfz": "0.", "bcqfz": "0."},
                "OUTPUT": {"lbldump": ".true.", "tstatsdump": "1.", "lfielddump": ".true.",
                           "tfielddump": "3.0", "fieldvars": "'th'"}})
            nf = len((d / "facets.inp.310").read_text().splitlines()) - 1
            (d / "Tfacinit.inp.310").write_text("# initial facet temperature [K]\n" + f"{tfac:.6f}\n" * nf)
            res = self.run_ok(d, nprocs=npx)
            self.assertIn("iwallmom = 4: 0 cells above horizontal Monin-Obukhov facets", res.stdout + res.stderr)
            with Dataset(d / "bldump_t.310.nc") as t, Dataset(d / "fac.310.nc") as f:
                out[npx] = {"wtheta": np.asarray(t["wtheta"][:], dtype=float), "ustar": np.asarray(t["ustar"][:], dtype=float),
                            "htc": np.asarray(f["htc"][:], dtype=float), "cth": np.asarray(f["cth"][:], dtype=float)}
            if npx == 1:
                with Dataset(d / "fielddump.000.000.310.nc") as f:
                    thl = np.asarray(f["thl"][-1], dtype=float)
                area = np.sum(np.loadtxt(d / "facetarea.inp.310", skiprows=1))
        r = out[1]
        self.assertEqual(len(r["wtheta"]), 3)
        self.assertTrue(np.all(r["wtheta"] > 0.))                           # the warm floor heats the air
        self.assertTrue(np.all(r["cth"] > 0.))
        np.testing.assert_allclose(r["htc"], r["cth"], rtol=1e-6)
        gain = np.sum(thl - 289.0) * (560. / 16) ** 2 * 20.0               # no solid c-cells
        expect = np.sum(r["wtheta"]) * 1.0 * area
        print(f"\n[offset floor] gain/expect {gain / expect:.7f}; htc {r['htc'][-1]} cth {r['cth'][-1]}")
        self.assertAlmostEqual(gain / expect, 1.0, delta=5e-4)
        for v in r:
            np.testing.assert_allclose(out[2][v], r[v], rtol=1e-6, err_msg=v)

    def test_startup_checks(self):
        """Start-up checks: iwallmom = 4 accepts iwalltemp 1 and 2, rejects 3; with lmoist
        iwallmoist = 2 needs iwalltemp = 2, and iwallmoist must be 1 or 2."""
        for name, walls, ok, msg in (
                ("chk_t3", {"iwalltemp": "3"}, False, "iwalltemp = 1 or 2"),
                ("chk_m3", {"iwalltemp": "2", "iwallmoist": "3"}, False, "iwallmoist = 1 or 2"),
                ("chk_m2_t1", {"iwalltemp": "1", "iwallmoist": "2"}, False, "requires iwalltemp = 2"),
                ("chk_m2_t2", {"iwalltemp": "2", "iwallmoist": "2"}, True, None),
                ("chk_t1", {"iwalltemp": "1"}, True, None)):
            with self.subTest(case=name):
                dd = tc.write_case(self.tmp / name, floor={"z0": self.Z0, "z0h": self.Z0H, "tfacinit": self.TFAC},
                                   groups={"RUN": {"runtime": "1."}, "WALLS": walls,
                                           "BC": {"bctfz": "0.", "bcqfz": "0."}})
                res = tc.run(dd)
                if ok:
                    self.assertEqual(res.returncode, 0, (res.stdout + res.stderr)[-2000:])
                else:
                    self.assertNotEqual(res.returncode, 0)
                    self.assertIn(msg, res.stdout + res.stderr)


RV_RD = 461.5 / 287.04
RLV, CP = 2.53e6, 1004.0   # check modglobal: rlv, cp


class TestMOWarmStart(_SolverTest):
    """Warm start with Monin-Obukhov facets (iwallmom = 4, #415): a cold run of 2N fixed steps equals
    a warm start of N steps from its own checkpoint at step N, bit for bit (bldump and the restart
    files at 2N), on 1 and 2 ranks. The warm start restores the per-rank initm files (the
    per-section L, tskin, qskin, u*, Cs, csu, fluxes and the mocell surface-layer gradients) and
    re-diagnoses the thermodynamic fields (presf, exnf, rhof, thv0h, dthvdz) from the restored ones.
    The checkpoint at N is the one a cold run of N steps writes (the same binary and arithmetic).

    The SGS model is Vreman with lbuoycorr, which reads dthvdz, so the restored surface-layer dthvdz
    above the floor is exercised. The one-equation closure is tested to a bound instead
    (test_one_equation_closure): its e12 ghost plane kb - 1, which subgrid reads before
    modibm:bottom resets it, is not in the restart files. That loss is not in the Monin-Obukhov path;
    a flat lbottom case loses it the same way."""

    N = 5
    VREMAN = {"loneeqn": ".false.", "lvreman": ".true.", "lbuoycorr": ".true."}
    RAD = {"iradiation": "2", "rad_longw": ".true.", "rka": "85.", "dlwtop": "70.", "dlwbot": "22."}

    def _cold(self, name, case, nprocx, sgs=None):
        groups = {"RUN": {"runtime": f"{2 * self.N}.", "trestart": f"{self.N}.", "nprocx": str(nprocx)},
                  "OUTPUT": {"lbldump": ".true.", "tstatsdump": "1."},
                  "NAMSUBGRID": dict(self.VREMAN if sgs is None else sgs)}
        if case == "t1":
            # prescribed facet fluxes, moist, RF01 cloud deck with the longwave scheme
            groups.update({"BC": {"bctfz": "-0.05", "bcqfz": "-5.0e-5"}, "NAMRADIATION": self.RAD})
            floor = {"z0": 2.e-4, "z0h": 2.e-4}
        else:
            # facet temperature; green floor that evaporates (iwallmoist = 2), so qskin carries state
            groups.update({"BC": {"bctfz": "0.", "bcqfz": "0."}, "WALLS": {"iwalltemp": "2", "iwallmoist": "2"},
                           "ENERGYBALANCE": {"wsoil": "313.", "wfc": "313."}})
            floor = {"z0": 2.e-4, "z0h": 2.e-5, "tfacinit": 291.0, "lgr": True}
        d = tc.write_case(self.tmp / name, u=7.0, v=-5.5, floor=floor, groups=groups)
        self.run_ok(d, nprocx)
        return d

    def _warm(self, cold, name, nprocx, skip=()):
        w = tc.stage_warm_start(cold, self.tmp / name, self.N, float(self.N), skip=skip)
        return w, tc.run(w, nprocx)

    def _assert_continues(self, cold, warm):
        for f in ("bldump_t.310.nc", "bldump_z.310.nc"):
            with Dataset(cold / f) as c, Dataset(warm / f) as w:
                for v in w.variables:
                    if "time" not in w[v].dimensions:
                        continue
                    with self.subTest(file=f, var=v):
                        np.testing.assert_array_equal(np.asarray(w[v][:]), np.asarray(c[v][self.N:]))
        final = tc.restart_set(cold, 2 * self.N)
        self.assertTrue(any(p.name.startswith("initm") for p in final))
        for p in final:
            with self.subTest(restart=p.name):
                self.assertEqual((warm / p.name).read_bytes(), p.read_bytes())

    def _chain(self, case, nprocx):
        cold = self._cold(f"ws_{case}_{nprocx}", case, nprocx)
        warm, res = self._warm(cold, f"ws_{case}_{nprocx}_warm", nprocx)
        self.assertEqual(res.returncode, 0, (res.stdout + res.stderr)[-3000:])
        self.assertIn("iwallmom = 4: restored the Monin-Obukhov state", res.stdout)
        self._assert_continues(cold, warm)
        return cold

    def test_prescribed_flux_moist_1_rank(self):
        self._chain("t1", 1)

    def test_prescribed_flux_moist_2_ranks(self):
        self._chain("t1", 2)

    def test_facet_temperature_1_rank(self):
        cold = self._chain("t2", 1)
        with Dataset(cold / "bldump_t.310.nc") as t:
            self.assertGreater(float(t["wq"][-1]), 0.)          # evaporating: qskin carries state

    def test_facet_temperature_2_ranks(self):
        self._chain("t2", 2)

    def test_missing_state_files_start_cold_with_warning(self):
        """Without the initm files (a restart written with another iwallmom) the run continues with a
        rank-0 warning and the Monin-Obukhov state starts as at a cold start, so it differs."""
        cold = self._cold("ws_noinitm", "t1", 1)
        skip = {p.name for p in tc.restart_set(cold, self.N) if p.name.startswith("initm")}
        warm, res = self._warm(cold, "ws_noinitm_warm", 1, skip=skip)
        self.assertEqual(res.returncode, 0, (res.stdout + res.stderr)[-3000:])
        self.assertIn("WARNING: iwallmom = 4: no Monin-Obukhov restart files", res.stdout + res.stderr)
        with Dataset(cold / "bldump_t.310.nc") as c, Dataset(warm / "bldump_t.310.nc") as w:
            self.assertFalse(np.array_equal(np.asarray(w["ustar"][:]), np.asarray(c["ustar"][self.N:])))

    def test_incomplete_state_set_aborts(self):
        cold = self._cold("ws_partial", "t1", 2)
        warm, res = self._warm(cold, "ws_partial_warm", 2, skip={f"initm{self.N:08d}_001_000.310"})
        self.assertNotEqual(res.returncode, 0)
        self.assertIn("initm00000005_001_000.310 not found, but other ranks have theirs", res.stdout + res.stderr)

    def test_one_equation_closure(self):
        """One-equation SGS: the warm start misses only the e12 ghost plane kb - 1 (see the class
        docstring), a round-off difference in the TKE next to the floor (measured: 4e-13 relative in
        e120, 3e-16 in u). The float32 bldump Monin-Obukhov diagnostics are therefore identical; a
        warm start that lost the Monin-Obukhov state differs there by ~1e-3 (test above)."""
        cold = self._cold("ws_oneeqn", "t1", 1, sgs={"loneeqn": ".true.", "lvreman": ".false."})
        warm, res = self._warm(cold, "ws_oneeqn_warm", 1)
        self.assertEqual(res.returncode, 0, (res.stdout + res.stderr)[-3000:])
        with Dataset(cold / "bldump_t.310.nc") as c, Dataset(warm / "bldump_t.310.nc") as w:
            for v in ("ustar", "obl", "thls", "qts", "wtheta", "wq", "u1", "v1"):
                with self.subTest(var=v):
                    np.testing.assert_array_equal(np.asarray(w[v][:]), np.asarray(c[v][self.N:]))


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

    def test_samples_are_spaced_by_tstatsdump_after_tstatstart(self):
        # tstatstart > tstatsdump: the first sample is tstatsdump after tstatstart and the
        # records are tstatsdump apart (no burst of per-step samples when tstatstart is reached)
        d = tc.write_case(self.tmp / "tstatstart", groups={
            "RUN": {"runtime": "9.5", "ladaptive": ".false.", "dtmax": "0.5"},
            "OUTPUT": {"lbldump": ".true.", "tstatsdump": "2.", "tstatstart": "5."}})
        self.run_ok(d)
        with Dataset(d / "bldump_t.310.nc") as t:
            times = np.asarray(t["time"][:], dtype=float)
        print(f"\n[bldump] record times {times}")
        np.testing.assert_allclose(np.diff(times), 2.0, atol=1e-3)
        self.assertEqual(times[0], 7.0)       # tstatstart + tstatsdump
        self.assertGreaterEqual(len(times), 2)

    def test_time_step_longer_than_tstatsdump_writes_once_per_step(self):
        # dt = 1 s > tstatsdump = 0.4 s: one record per step, none duplicated. (With a constant
        # dt the catch-up of tnext does not change the records; it keeps tnext on the grid.)
        d = tc.write_case(self.tmp / "dtlong", groups={
            "RUN": {"runtime": "4.5", "ladaptive": ".false.", "dtmax": "1."},
            "OUTPUT": {"lbldump": ".true.", "tstatsdump": "0.4"}})
        self.run_ok(d)
        with Dataset(d / "bldump_t.310.nc") as t:
            times = np.asarray(t["time"][:], dtype=float)
        np.testing.assert_array_equal(times, [1.0, 2.0, 3.0, 4.0, 5.0])  # runtime 4.5 ends with the 5th step

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

    def test_warm_start(self):
        """A warm start re-diagnoses rhof and exnf before the first longwave update (they were zero,
        a floating-point exception), updates at its first step (thlprad is not in the restart files)
        and then rejoins the timerad grid of the run it continues: restart at 5 s, timerad = 4 s, so
        updates at 6 s (the first step) and 8 s, and thlprad is held over the step ending at 7 s
        (before the fix tnextrad started at 0 and caught up one interval per step: an update at 7 s)."""
        rad = dict(self.RAD, timerad="4.")
        d = self._run("rad_cold", rad, groups={"RUN": {"runtime": "10.", "trestart": "5."},
                                               "OUTPUT": {"lbldump": ".true.", "tstatsdump": "1."}})
        w = tc.stage_warm_start(d, self.tmp / "rad_warm", 5, 5.0)
        self.run_ok(w)
        with Dataset(w / "bldump_z.310.nc") as z:
            t = np.asarray(z["time"][:])
            q = np.asarray(z["thlprad"][:], dtype=float)
        np.testing.assert_array_equal(t, [6., 7., 8., 9., 10.])
        self.assertTrue(np.all(np.isfinite(q)))
        self.assertLess(q.min(), -1e-4)                      # cloud-top cooling present
        np.testing.assert_array_equal(q[1], q[0])            # held: no update at 7 s
        self.assertFalse(np.array_equal(q[2], q[1]))         # update at 8 s
        np.testing.assert_array_equal(q[3], q[2])            # held at 9 s

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
        # The top three cells are excluded: the top boundary reaches into them (thl of the
        # last cell: FT term -1.3e-2 K against subsidence -4e-3 K; second last: ratio 0.07),
        # so they measure the boundary, not the balance. Below that the two terms cancel to
        # a few 1e-3 of the subsidence signal, which is one float32 step of thl (3e-5 K at
        # 297 K) over the signal (1.2e-2 K).
        dz = zf[1] - zf[0]
        ft = (zf > 840.0 + 150.0) & (zf <= zf[-4])
        ratio = np.abs(dthl_both[ft]).max() / np.abs(dthl_subs[ft]).max()
        print(f"\n[balance] max|dthl_both|/max|dthl_subs| in FT = {ratio:.4f}")
        # The FT term is D/3 z (z - zi)^(-2/3) (flux form, cell mean); the subsidence term is the
        # upwind difference w (thl_{k+1} - thl_k)/dz of the initial profile
        # thl_i + (z - 840 m)^(1/3). Their difference, evaluated here for the window, is the
        # truncation error of the stencil plus the 7 m offset between the 840 m inversion and
        # zi = 832.7 m (about 1 % at the lowest level). Twice that, in addition to the float32
        # step, is the bound. A sign error gives a ratio of ~2, a missing term 1.
        zi = self._ft_term()[1]
        thl_i = tc.rf01_profiles(zf)[0]
        sub = np.zeros_like(zf)
        sub[:-1] = self.D * zf[:-1] * (thl_i[1:] - thl_i[:-1]) / dz
        ft_term = self.D / 3 * zf * np.maximum(zf - zi, 1e-9) ** (-2 / 3)
        model = np.abs(ft_term - sub)[ft].max() / np.abs(sub)[ft].max()
        print(f"[balance] stencil model {model:.4f}")
        self.assertLess(ratio, 2 * model)


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

    def test_leftover_namthermostat_aborts(self):
        # the group was replaced by lbldump in &OUTPUT; ignoring it would silently drop the output
        self._abort("namthermostat", {"NAMTHERMOSTAT": {"lthermostat": ".true.", "tthermostat": "60."}},
                    "lbldump")

    def test_bldump_without_interval_aborts(self):
        self._abort("bld_noint", {"OUTPUT": {"lbldump": ".true.", "tstatsdump": "0."}},
                    "requires tstatsdump > 0")

    def test_negative_radiation_parameters_abort(self):
        self._abort("rad_negflux", {"NAMRADIATION": dict(self.RAD, dlwbot="-5.")},
                    "dlwtop and dlwbot must be >= 0")
        self._abort("rad_negtime", {"NAMRADIATION": dict(self.RAD, timerad="-1.")},
                    "timerad must be >= 0")

    def test_unreadable_namelist_aborts_once_on_all_ranks(self):
        # a read error in a new namelist group is collective: rank 0 reports once, all ranks abort.
        # An unknown variable name is a read error (iostat > 0) under both gfortran and ifort; a
        # bad integer value ('banana') is not portable: gfortran 13-15 return iostat = -1 for it.
        d = tc.write_case(self.tmp / "nml_bad", groups={
            "RUN": {"nprocx": "2"}, "NAMRADIATION": {"iradiation": "2", "no_such_variable": "1"}})
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
                # Only the reduction order of the allreduce slab means differs between splits.
                # Measured over this run: identical to the last bit except the resolved flux
                # (a near-cancelling mean of w*thl terms), 1e-7 relative. A decomposition or
                # halo bug gives O(1) differences from the first step, so 1e-6 relative (and
                # 1e-6 of the field maximum absolute) leaves 10x over the round-off and 6
                # orders below a bug.
                atol_frac = 1e-6
                with self.subTest(split=label, var=v):
                    np.testing.assert_allclose(out[label][v], ref, rtol=1e-6, atol=atol_frac * max(np.abs(ref).max(), 1e-30),
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


class TestMOSlopedFacet(_SolverTest):
    """Monin-Obukhov on geometrically sloped facets: a periodic 20 degree ridge STL (below
    mo_maxslope = 30) through the Python IBM preprocessing (udprep, legacy backend), run with
    iwallmom = 4, iwalltemp = 1 on 1 and 2 ranks. Unlike the declared-normal tests above, the
    sections here are not grid-aligned: reconstruction point (trilinear thl, qt), section
    distance, rotated stress, w-grid sections and bctfz on sloped sections are all exercised."""

    ANGLE = 20.0
    Z0 = 2.e-4
    BCTFZ = -0.2          # K m/s, upward (negative = into the air)
    RUNTIME = 3.0         # three 1-s steps

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        if not tc.ibm_preprocessor_available():
            cls._tmp.cleanup()
            raise unittest.SkipTest("no gfortran for the udprep legacy IBM preprocessor")
        try:
            cls.prep = tc.prepare_ridge(cls.tmp / "prep", cls.ANGLE, cls.Z0)
        except ImportError as exc:
            cls._tmp.cleanup()
            raise unittest.SkipTest(f"udprep dependencies missing: {exc}")

    def _run(self, name, nprocx):
        # uniform, unsaturated thl and qt: the heat budget below needs a known initial integral
        d = tc.write_ridge_run(self.prep, self.tmp / name, u=7.0, v=0.0, thl=np.full(64, 289.0),
                               qt=np.full(64, 5.0e-3), groups={
            "RUN": {"runtime": str(self.RUNTIME), "nprocx": str(nprocx)},
            "BC": {"bctfz": str(self.BCTFZ), "bcqfz": "0."},
            "OUTPUT": {"lbldump": ".true.", "tstatsdump": "1.", "lfielddump": ".true.",
                       "tfielddump": str(self.RUNTIME), "fieldvars": "'th'"}})
        res = self.run_ok(d, nprocs=nprocx)
        return d, res.stdout + res.stderr

    def test_sloped_mo_facets(self):
        d1, out1 = self._run("ridge_1", 1)
        d2, out2 = self._run("ridge_2", 2)
        for out in (out1, out2):
            self.assertIn("iwallmom = 4: 4 facets use Monin-Obukhov (slope <=  30.0 deg), 0 use the neutral wall function", out)
            self.assertIn("iwallmom = 4: 0 cells above horizontal Monin-Obukhov facets", out)
        us = []
        for d in (d1, d2):
            with Dataset(d / "bldump_t.310.nc") as t:
                us.append(np.asarray(t["ustar"][:], dtype=float))
        for u in us:
            self.assertEqual(len(u), 3)
            self.assertTrue(np.all(np.isfinite(u)))
            self.assertTrue(np.all(u > 0.01))
        # the same sections and arithmetic: only the order of the area-weighted sum differs
        np.testing.assert_allclose(us[0], us[1], rtol=1e-6)
        # heat input: the air gains bctfz x facet area per unit time (thl was uniform, so the
        # integral over the fluid cells changes only through the facet flux; every section
        # applies flux * section area / (dx dy dzf)).
        dx = dy = 560.0 / 16
        dz = 20.0
        solid = np.zeros((16, 16, 64), dtype=bool)
        for l in (d1 / "solid_c.txt").read_text().splitlines():
            if l.strip() and not l.lstrip().startswith("#"):
                i, j, k = (int(x) for x in l.split()[:3])
                solid[i - 1, j - 1, k - 1] = True
        with Dataset(d1 / "fielddump.000.000.310.nc") as f:
            thl = np.asarray(f["thl"][-1], dtype=float)        # (k, j, i)
        gain = np.sum((thl.transpose(2, 1, 0) - 289.0)[~solid]) * dx * dy * dz
        area = np.sum(np.loadtxt(d1 / "facetarea.inp.310", skiprows=1))
        expect = -self.BCTFZ * area * self.RUNTIME
        # Tolerance 3e-3: the float32 thl record (ulp 3e-5 K at 289 K) gives a sum error of about
        # 1.3e-4 of the signal (random rounding over ~15000 cells), the rest is the conservation
        # error of the IBM advection and diffusion corrections near the ridge; measured 9e-4.
        self.assertAlmostEqual(gain / expect, 1.0, delta=3e-3)


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
