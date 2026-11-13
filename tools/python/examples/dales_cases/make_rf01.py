#!/usr/bin/env python3
"""Generate the uDALES inputs of the DYCOMS-II RF01 case (examples/303).

Usage: python make_rf01.py <out_dir> [--expnr 303 --dz 5 --ktot 300 --itot 96 --dx 35]

Specification: Stevens, B., et al., 2005: Evaluation of large-eddy simulations
via observations of nocturnal marine stratocumulus. Mon. Wea. Rev., 133,
1443-1462 (doi:10.1175/MWR2930.1). The numbers below were cross-checked against
the PyCLES case set-up (InitDYCOMS_RF01 in Initialization.pyx,
ForcingDyCOMS_RF01 in Forcing.pyx, SurfaceDYCOMS_RF01 in Surface.pyx,
RadiationDyCOMS_RF01 in Radiation.pyx); the paper itself was not available
here, so the PyCLES values are the cited source:

  profiles   thl = 289 K, qt = 9 g/kg for z <= 840 m;
             thl = 297.5 + (z - 840)^(1/3) K, qt = 1.5 g/kg above   (InitDYCOMS_RF01)
  winds      u = 7, v = -5.5 m/s, geostrophic wind equal to them
  subsidence w = -D z with D = 3.75e-6 1/s                          (ForcingDyCOMS_RF01)
  Coriolis   f = 2 Omega sin(31.5 deg) = 0.76e-4 1/s                (ForcingDyCOMS_RF01)
  surface    SHF = 15 W/m2, LHF = 115 W/m2, Cm = 0.0011             (SurfaceDYCOMS_RF01)
  radiation  F0 = 70, F1 = 22 W/m2, kappa = 85 m2/kg, D = 3.75e-6  (RadiationDyCOMS_RF01)

Files written: prof.inp.<expnr> (z thl qt u v e12) and lscale.inp.<expnr>
(z ug vg pgx pgy wfls dqtdxls dqtdyls dqtdtls thlpcar). The surface values
needed in namoptions are printed.

Choices that are not in the references:
 * e12 (= sqrt of the SGS TKE in uDALES, as in DALES) is 1 m/s below the initial
   inversion zi = 840 m and 0 above; uDALES raises 0 to e12min = 5e-5 itself, which
   is written directly here so that no start-up message appears. PyCLES starts from
   zero SGS TKE below 200 m.
 * Fluxes: uDALES prescribed fluxes are positive in -z (decision D1), so the upward
   RF01 fluxes are written with a negative sign.
 * Latent heat: the qt flux uses uDALES rlv = 2.26e6 J/kg (PyCLES divides by its own
   temperature-dependent L, about 2.5e6, and by rho = 1.22).
 * z0: the drag of RF01 is a constant Cm = 0.0011 (PyCLES u* = sqrt(Cm) U). With
   BCbotm = 4 uDALES computes Cm = kappa^2 / ln(z1/z0)^2 (neutral) at the first
   level z1 = dz/2, so z0 = z1 exp(-kappa / sqrt(Cm)), also used for z0h.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from convert_dales_case import write_table  # noqa: E402

# Constants taken from uDALES modglobal (cp, rlv) and PyCLES (rho, kappa_vk is DALES fkar)
CP = 1004.0       # J/kg/K, modglobal cp
RLV = 2.26e6      # J/kg, modglobal rlv
RHO = 1.22        # kg/m3, SurfaceDYCOMS_RF01 (qt_flux = fq / lv / 1.22)
FKAR = 0.41       # von Karman constant, modglobal fkar
SHF, LHF = 15.0, 115.0   # W/m2, upward
CM = 0.0011
ZI = 840.0
D = 3.75e-6
UG, VG = 7.0, -5.5


def profiles(zf: np.ndarray):
    below = zf <= ZI
    thl = np.where(below, 289.0, 297.5 + np.cbrt(np.maximum(zf - ZI, 0.0)))
    qt = np.where(below, 9.0e-3, 1.5e-3)
    e12 = np.where(below, 1.0, 5.0e-5)
    return thl, qt, e12


def surface(dz: float) -> dict:
    z1 = 0.5 * dz
    z0 = z1 * np.exp(-FKAR / np.sqrt(CM))
    return {
        "z0": z0,
        "wtsurf": -SHF / (RHO * CP),          # K m/s, positive downward (D1)
        "wqsurf": -LHF / (RHO * RLV),         # kg/kg m/s, positive downward (D1)
        "cm_check": FKAR ** 2 / np.log(z1 / z0) ** 2,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--expnr", default="303")
    ap.add_argument("--dz", type=float, default=5.0)
    ap.add_argument("--ktot", type=int, default=300)
    ap.add_argument("--itot", type=int, default=96)
    ap.add_argument("--dx", type=float, default=35.0)
    a = ap.parse_args(argv)

    a.out_dir.mkdir(parents=True, exist_ok=True)
    zf = (np.arange(a.ktot) + 0.5) * a.dz
    thl, qt, e12 = profiles(zf)
    n = a.ktot
    u = np.full(n, UG)
    v = np.full(n, VG)
    write_table(a.out_dir / f"prof.inp.{a.expnr}", "DYCOMS-II RF01 (Stevens et al. 2005; PyCLES InitDYCOMS_RF01)",
                "z thl qt u v e12", np.column_stack([zf, thl, qt, u, v, e12]))
    zero = np.zeros(n)
    wfls = -D * zf
    write_table(a.out_dir / f"lscale.inp.{a.expnr}", "DYCOMS-II RF01 (ForcingDyCOMS_RF01: ug, vg = u, v; wfls = -D z)",
                "z ug vg pgx pgy wfls dqtdxls dqtdyls dqtdtls thlpcar",
                np.column_stack([zf, u, v, zero, zero, wfls, zero, zero, zero, zero]))

    s = surface(a.dz)
    print(f"wrote prof.inp.{a.expnr}, lscale.inp.{a.expnr} in {a.out_dir}: ktot = {n}, "
          f"ztop = {n * a.dz:g} m, domain = {a.itot * a.dx:g} m")
    print(f"z0 = z0h = {s['z0']:.6e} m  (neutral Cm at z1 = {0.5 * a.dz:g} m: {s['cm_check']:.6f})")
    print(f"wtsurf = {s['wtsurf']:.8e} K m/s  (SHF {SHF} W/m2 upward, rho {RHO}, cp {CP:g}; negative = upward)")
    print(f"wqsurf = {s['wqsurf']:.8e} kg/kg m/s  (LHF {LHF} W/m2 upward, rho {RHO}, rlv {RLV:g})")
    print(f"f = {2 * 7.292e-5 * np.sin(np.radians(31.5)):.4e} 1/s at xlat = 31.5")
    return 0


if __name__ == "__main__":
    sys.exit(main())
