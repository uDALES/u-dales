"""Programmatic flat-bottom test cases for the thermo integration tests."""
from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
UDALES_BUILD = Path(os.environ.get("UDALES_BUILD", REPO_ROOT / "build" / "debug" / "u-dales"))
RUNTIME_MODULES = os.environ.get(
    "UDALES_RUNTIME_MODULES",
    "intel/2021a netCDF/4.8.0-iimpi-2021a netCDF-Fortran/4.5.3-iimpi-2021a FFTW/3.3.9-intel-2021a",
)
sys.path.insert(0, str(REPO_ROOT / "tests" / "integration"))
sys.path.insert(0, str(REPO_ROOT / "tools" / "python"))
import mpi_env  # noqa: E402
from udprep.flat_floor import write_flat_floor  # noqa: E402

EXPNR = "310"

# Namelist groups as {group: {key: value-string}}. Values are Fortran literals.
BASE = {
    "RUN": {"iexpnr": EXPNR, "runtime": "0.5", "trestart": "1000.", "dtmax": "1.0",
            "ladaptive": ".false.", "libm": ".false.", "lles": ".true.",
            "nprocx": "1", "nprocy": "1", "lrandomize": ".false."},
    "OUTPUT": {"lfielddump": ".false.", "ltdump": ".false.", "lxytdump": ".false."},
    "DOMAIN": {"itot": "16", "jtot": "16", "ktot": "64", "xlen": "560.", "ylen": "560."},
    "PHYSICS": {"ltempeq": ".true.", "lmoist": ".true.", "lbuoyancy": ".true.",
                "lcoriol": ".false.", "ps": "101780."},
    "NAMSUBGRID": {"loneeqn": ".true.", "lvreman": ".false.", "ldelta": ".false."},
    "WALLS": {"nfcts": "0", "lbottom": ".true."},
    "BC": {"BCbotm": "3", "BCbotT": "1", "BCbotq": "1", "thls": "289.0", "qts": "9.0e-3",
           "wtsurf": "0.", "wqsurf": "0.", "z0": "1.e-4", "z0h": "1.e-4"},
}


def rf01_profiles(zf: np.ndarray, zi: float = 840.0):
    """DYCOMS-II RF01 initial thl [K] and qt [kg/kg] (Stevens et al. 2005)."""
    thl = np.where(zf <= zi, 289.0, 297.5 + np.cbrt(np.maximum(zf - zi, 0.0)))
    qt = np.where(zf <= zi, 9.0e-3, 1.5e-3)
    return thl, qt


def write_case(run_dir: Path, groups: dict | None = None, ktot: int = 64, dz: float = 20.0,
               thl=None, qt=None, u=7.0, v=-5.5, wfls=None, thlpcar=None, sv=None,
               floor: dict | None = None, ug=None, vg=None) -> Path:
    """Write namoptions/prof.inp/lscale.inp for a flat case. Returns run_dir.

    floor: if given, the bottom is a facet floor (udprep.flat_floor) with Monin-Obukhov
    facets (iwallmom = 4) instead of the flat bottom boundary (lbottom). Its keys are the
    keyword arguments of write_flat_floor (z0 and z0h are required); the facet fluxes
    bctfz/bcqfz go in groups["BC"]. groups are applied after the floor settings.
    ug, vg: geostrophic wind of lscale.inp (default u, v). Note that uDALES applies the
    geostrophic pressure gradient (modstartup: dpdxl = om23_gs vg) even with lcoriol off."""
    run_dir.mkdir(parents=True, exist_ok=True)
    nml = {g: dict(kv) for g, kv in BASE.items()}
    if floor is not None:
        itot, jtot = int(nml["DOMAIN"]["itot"]), int(nml["DOMAIN"]["jtot"])
        for g, kv in (groups or {}).items():          # a domain override must reach the floor too
            if g == "DOMAIN":
                itot, jtot = int(kv.get("itot", itot)), int(kv.get("jtot", jtot))
        xlen = float(((groups or {}).get("DOMAIN") or {}).get("xlen", nml["DOMAIN"]["xlen"]))
        ylen = float(((groups or {}).get("DOMAIN") or {}).get("ylen", nml["DOMAIN"]["ylen"]))
        counts = write_flat_floor(run_dir, EXPNR, itot, jtot, xlen / itot, ylen / jtot, z1=0.5 * dz, **floor)
        nml["RUN"]["libm"] = ".true."
        nml["WALLS"] = {"lbottom": ".false.", "iwallmom": "4", "iwalltemp": "1", "iwallmoist": "1",
                        **{k: str(v) for k, v in counts.items()}}
        nml["BC"] = {k: v for k, v in nml["BC"].items() if k in ("thls", "qts")}
    for g, kv in (groups or {}).items():
        nml.setdefault(g, {}).update(kv)
    nml["DOMAIN"]["ktot"] = str(ktot)
    zf = (np.arange(ktot) + 0.5) * dz
    if thl is None or qt is None:
        thl, qt = rf01_profiles(zf)
    with open(run_dir / f"namoptions.{EXPNR}", "w") as f:
        for g, kv in nml.items():
            f.write(f"&{g}\n")
            for k, val in kv.items():
                f.write(f"{k:<12} = {val}\n")
            f.write("/\n\n")
    with open(run_dir / f"prof.inp.{EXPNR}", "w") as f:
        f.write("# thermo test case\n# z thl qt u v e12\n")
        for k in range(ktot):
            f.write(f"{zf[k]:.6f} {thl[k]:.6f} {qt[k]:.8e} {u:.4f} {v:.4f} 0.1\n")
    ugv = u if ug is None else ug
    vgv = v if vg is None else vg
    wfls = np.zeros(ktot) if wfls is None else wfls
    thlpcar = np.zeros(ktot) if thlpcar is None else thlpcar
    with open(run_dir / f"lscale.inp.{EXPNR}", "w") as f:
        f.write("# thermo test case\n# z ug vg pgx pgy wfls dqtdxls dqtdyls dqtdtls thlpcar\n")
        for k in range(ktot):
            f.write(f"{zf[k]:.6f} {ugv:.4f} {vgv:.4f} 0. 0. {wfls[k]:.6e} 0. 0. 0. {thlpcar[k]:.6e}\n")
    if sv is not None:
        sv = np.asarray(sv, dtype=float).reshape(ktot, -1)
        with open(run_dir / f"scalar.inp.{EXPNR}", "w") as f:
            f.write("# thermo test case\n# z sv1..svn\n")
            for k in range(ktot):
                f.write(f"{zf[k]:.6f} " + " ".join(f"{x:.8e}" for x in sv[k]) + "\n")
    return run_dir


def _shell_prefix() -> str:
    mods = Path("/etc/profile.d/modules.sh")
    if mods.is_file() and RUNTIME_MODULES:
        return f"source {mods} && module load {RUNTIME_MODULES} && export HDF5_USE_FILE_LOCKING=FALSE && "
    return "export HDF5_USE_FILE_LOCKING=FALSE && "


def run(run_dir: Path, nprocs: int = 1) -> subprocess.CompletedProcess:
    mpi_env.preflight(nprocs)
    mpiexec, extra = mpi_env.mpi_exec_and_args()
    cmd = (f"{_shell_prefix()}cd '{run_dir}' && '{mpiexec}' {extra} -n {nprocs} "
           f"'{UDALES_BUILD}' namoptions.{EXPNR}")
    return subprocess.run(["bash", "-lc", cmd], cwd=REPO_ROOT, capture_output=True, text=True)


FKAR, GRAV, RD, RV = 0.40, 9.81, 287.04, 461.5


def psim(z):
    if z <= 0:
        x = (1 - 16 * z) ** 0.25
        return np.pi / 2 - 2 * np.arctan(x) + np.log((1 + x) ** 2 * (1 + x ** 2) / 8)
    return -2 / 3 * (z - 5 / 0.35) * np.exp(-0.35 * z) - z - (10 / 3) / 0.35


def psih(z):
    if z <= 0:
        x = (1 - 16 * z) ** 0.25
        return 2 * np.log((1 + x ** 2) / 2)
    return -2 / 3 * (z - 5 / 0.35) * np.exp(-0.35 * z) - (1 + 2 / 3 * z) ** 1.5 - (10 / 3) / 0.35 + 1


def rib_from_obl(z, L, z0h, z0m):
    return z / L * (np.log(z / z0h) - psih(z / L) + psih(z0h / L)) / \
        (np.log(z / z0m) - psim(z / L) + psim(z0m / L)) ** 2


def obl_iter(thl, qt, tskin, qskin, z, z0m, z0h, u, v, L):
    """Port of DALES dev 075ac4e modsurface.f90:calc_obl_iter."""
    thv = thl * (1 + (RV / RD - 1) * qt)
    thvs = tskin * (1 + (RV / RD - 1) * qskin)
    horv2 = max(u * u + v * v, 0.1 ** 2)
    rib = GRAV / thvs * z * (thv - thvs) / horv2
    if rib == 0:
        return 1e6
    if rib * L < 0 or abs(L) == 1e5:
        L = 0.01 if rib > 0 else -0.01
    Lold, it = 0.0, 0
    while abs((L - Lold) / L) > 1e-4 and it < 1000:
        it += 1
        Lold = L
        Ls, Le = L - 0.001 * L, L + 0.001 * L
        fx = rib - rib_from_obl(z, L, z0h, z0m)
        fxd = (-rib_from_obl(z, Ls, z0h, z0m) + rib_from_obl(z, Le, z0h, z0m)) / (Ls - Le)
        L = L - fx / fxd
        if rib * L < 0 or abs(L) == 1e5:
            L = 0.01 if rib > 0 else -0.01
        if abs(L) > 1e6:
            L = np.sign(L) * 1e6
    return L
