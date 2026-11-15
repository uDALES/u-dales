#!/usr/bin/env python3
"""Write a constant-density DALES baseprof.inp (ibas_prf = 5) that matches the uDALES base state.

Usage:
    python make_dales_baseprof.py <dales prof.inp> <ps [Pa]> <thls [K]> <qts [kg/kg]> <out baseprof.inp>
                                  [--rho 1.0]

uDALES is Boussinesq (rhobf = rhobh = 1, modfields) and computes the Exner function from the
hydrostatic slab-mean pressure (modthermodynamics fromztop). DALES with ibas_prf = 5 reads
rhobf, pb and tb from baseprof.inp.<iexpnr> (rows k = 1 .. kmax+1: height rhobf pb tb) and
sets exnf = (pb/pref0)^(rd/cp) once. This script gives DALES
  - a constant rhobf (default 1.0), so every anelastic factor in DALES (advection, Poisson,
    SGS) is 1, as in uDALES;
  - pb = the pressure that uDALES fromztop computes from the initial profile: the moist
    hydrostatic integral with thvh(kb) = thls (1 + (rv/rd - 1) qts), half-level thv from the
    level means, and ql from exact saturation adjustment (Tetens, rlv = 2.53e6, uDALES lqlnr);
  - tb = th0 * exnf (temperature, written for completeness; DALES does not use it here).
The level kmax+1 repeats thl, qt of level kmax (uDALES zero-gradient top ghost).
Constants are those of uDALES modglobal (identical to DALES dev 075ac4e except es).
"""
from __future__ import annotations

import argparse
import sys

import numpy as np

GRAV, RD, RV, CP, RLV, PREF0 = 9.81, 287.04, 461.5, 1004.0, 2.53e6, 1.0e5
ES0, AT, BT, TMELT = 610.78, 17.27, 35.86, 273.16
RDOCP = RD / CP


def read_prof(path):
    rows = [ln.split() for ln in open(path) if ln.strip() and not ln.lstrip().startswith("#")]
    a = np.array([[float(x) for x in r[:3]] for r in rows])
    return a[:, 0], a[:, 1], a[:, 2]


def ql_exact(thl, qt, p, exn):
    """Newton-Raphson saturation adjustment as uDALES thermo with lqlnr."""
    tl = thl * exn
    T = tl.copy()
    for _ in range(50):
        es = ES0 * np.exp(AT * (T - TMELT) / (T - BT))
        qs = RD / RV * es / (p - (1 - RD / RV) * es)
        T = T - (T + RLV / CP * qs - tl - RLV / CP * qt) / (1 + RLV**2 * qs / (RV * CP * T**2))
    es = ES0 * np.exp(AT * (T - TMELT) / (T - BT))
    qs = RD / RV * es / (p - (1 - RD / RV) * es)
    return np.maximum(qt - qs, 0.0)


def base_state(zf, thl, qt, ps, thls, qts):
    n = zf.size
    zh = np.zeros(n + 1)
    for k in range(n):
        zh[k + 1] = 2 * zf[k] - zh[k]
    dzf = np.diff(zh)
    dzh = np.empty(n)
    dzh[0] = zf[0]
    dzh[1:] = np.diff(zf)
    thvs = thls * (1 + (RV / RD - 1) * qts)
    exn = np.full(n, (ps / PREF0) ** RDOCP)
    ql = np.zeros(n)
    pres = np.full(n, ps)
    for _ in range(20):
        th = thl + RLV / (CP * exn) * ql
        thh = np.empty(n); qth = np.empty(n); qlh = np.empty(n)
        for k in range(1, n):
            w = 2 * dzh[k]
            thh[k] = (th[k] * dzf[k - 1] + th[k - 1] * dzf[k]) / w
            qth[k] = (qt[k] * dzf[k - 1] + qt[k - 1] * dzf[k]) / w
            qlh[k] = (ql[k] * dzf[k - 1] + ql[k - 1] * dzf[k]) / w
        pres[0] = (ps**RDOCP - GRAV * PREF0**RDOCP * zf[0] / (CP * thvs)) ** (1 / RDOCP)
        for k in range(1, n):
            thvh = thh[k] * (1 + (RV / RD - 1) * qth[k] - RV / RD * qlh[k])
            pres[k] = (pres[k - 1] ** RDOCP - GRAV * PREF0**RDOCP * dzh[k] / (CP * thvh)) ** (1 / RDOCP)
        exn = (pres / PREF0) ** RDOCP
        ql = ql_exact(thl, qt, pres, exn)
    th = thl + RLV / (CP * exn) * ql
    return pres, th * exn, ql


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("prof"); ap.add_argument("ps", type=float); ap.add_argument("thls", type=float)
    ap.add_argument("qts", type=float); ap.add_argument("out")
    ap.add_argument("--rho", type=float, default=1.0, help="constant rhobf [kg/m3] (uDALES: 1)")
    a = ap.parse_args(argv)
    zf, thl, qt = read_prof(a.prof)
    zf = np.append(zf, 2 * zf[-1] - zf[-2])        # level kmax+1
    thl = np.append(thl, thl[-1]); qt = np.append(qt, qt[-1])
    pb, tb, ql = base_state(zf, thl, qt, a.ps, a.thls, a.qts)
    with open(a.out, "w") as f:
        f.write(" #baseprofiles (constant rhobf; pb from the uDALES fromztop integral, make_dales_baseprof.py)\n")
        f.write(" #height rhobf pb tb\n")
        for k in range(zf.size):
            f.write(f"{zf[k]:9.2f} {a.rho:24.17E} {pb[k]:24.17E} {tb[k]:24.17E}\n")
    print(f"wrote {a.out}: {zf.size} rows, pb {pb[0]:.2f} .. {pb[-1]:.2f} Pa, max ql {ql.max()*1e3:.4f} g/kg")
    return 0


if __name__ == "__main__":
    sys.exit(main())
