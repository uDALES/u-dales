#!/usr/bin/env python3
"""Write the inputs of example 304: a dry convective boundary layer over a floor of fixed temperature.

Usage (from the repository root):
    python tools/python/examples/dales_cases/make_cbl304.py examples/304 \
        tools/python/examples/dales_cases/dales_inputs/cbl304

Writes, for uDALES, prof.inp.304, lscale.inp.304 and the flat-floor IBM files with Tfacinit.inp.304
(udprep.flat_floor, one horizontal facet at the fixed temperature TFAC, lGR = 0), and for DALES
dev 075ac4e prof.inp.001, lscale.inp.001 and the constant-density baseprof.inp.001
(make_dales_baseprof.py). The namoptions files of both codes are written by hand (examples/304,
dales_inputs/cbl304) from the constants printed at the end.

The case (an analysis choice, not a published case): dry, thl = 300 K up to 500 m, 6 K/km above,
geostrophic wind (5, 0) m/s at xlat = 45, floor at 302 K with z0 = 0.1 m, z0h = 0.01 m.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from udprep.flat_floor import write_flat_floor  # noqa: E402

EXPNR = "304"
ITOT = JTOT = 96
DX = DY = 50.0
KTOT, DZ = 80, 25.0
THL0, ZML, GAMMA = 300.0, 500.0, 6.0e-3     # K, m, K/m
UG, VG = 5.0, 0.0                           # m/s
E12 = 0.3                                   # m/s, initial sqrt(e) (both codes read e12)
TFAC = 302.0                                # K, facet (DALES thls) temperature
Z0, Z0H = 0.1, 0.01                         # m
PS = 1.0e5                                  # Pa


def profiles():
    zf = (np.arange(KTOT) + 0.5) * DZ
    thl = np.where(zf <= ZML, THL0, THL0 + GAMMA * (zf - ZML))
    return zf, thl


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    udir, ddir = Path(argv[0]), Path(argv[1])
    udir.mkdir(parents=True, exist_ok=True)
    ddir.mkdir(parents=True, exist_ok=True)
    zf, thl = profiles()

    with open(udir / f"prof.inp.{EXPNR}", "w") as f:
        f.write("# dry CBL over a fixed-temperature floor (make_cbl304.py)\n# z thl qt u v e12\n")
        for k in range(KTOT):
            f.write(f"{zf[k]:.9e} {thl[k]:.9e} {0.0:.9e} {UG:.9e} {VG:.9e} {E12:.9e}\n")
    with open(udir / f"lscale.inp.{EXPNR}", "w") as f:
        f.write("# dry CBL over a fixed-temperature floor (make_cbl304.py)\n"
                "# z ug vg pgx pgy wfls dqtdxls dqtdyls dqtdtls thlpcar\n")
        for k in range(KTOT):
            f.write(f"{zf[k]:.9e} {UG:.9e} {VG:.9e} 0. 0. 0. 0. 0. 0. 0.\n")
    counts = write_flat_floor(udir, EXPNR, ITOT, JTOT, DX, DY, z1=0.5 * DZ, z0=Z0, z0h=Z0H,
                              lgr=False, tfacinit=TFAC)

    with open(ddir / "prof.inp.001", "w") as f:
        f.write("# dry CBL over a fixed-temperature floor (uDALES examples/304, make_cbl304.py)\n"
                "# zf thl qt u v tke\n")
        for k in range(KTOT):
            f.write(f"{zf[k]:.9e} {thl[k]:.9e} {0.0:.9e} {UG:.9e} {VG:.9e} {E12:.9e}\n")
    with open(ddir / "lscale.inp.001", "w") as f:
        f.write("# dry CBL over a fixed-temperature floor (uDALES examples/304, make_cbl304.py)\n"
                "# height ug vg wfls dqtdxls dqtdyls dqtdtls dthlrad\n")
        for k in range(KTOT):
            f.write(f"{zf[k]:.9e} {UG:.9e} {VG:.9e} 0. 0. 0. 0. 0.\n")
    subprocess.run([sys.executable, str(HERE / "make_dales_baseprof.py"), str(ddir / "prof.inp.001"),
                    str(PS), str(TFAC), "0.0", str(ddir / "baseprof.inp.001")], check=True)

    print("&WALLS counts:")
    for k, v in counts.items():
        print(f"  {k:<12} = {v}")
    print(f"thls = {TFAC} (uDALES &BC, DALES &NAMSURFACE), z0 = {Z0}, z0h = {Z0H}, ps = {PS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
