"""Programmatic flat-bottom test cases for the thermo integration tests."""
from __future__ import annotations

import os
import re
import shutil
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
from udprep.flat_floor import write_flat_floor, write_floor_factypes  # noqa: E402

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


def write_ridge_stl(path: Path, angle_deg: float, length: float, width: float) -> None:
    """ASCII STL of a periodic ridge z = h (1 - |x - length/2| / (length/2)), invariant in y.
    Both slopes face upward (normals have z > 0) and are tilted by angle_deg from horizontal."""
    h = np.tan(np.radians(angle_deg)) * length / 2
    tris = []
    for x0, z0, x1, z1 in ((0., 0., length / 2, h), (length / 2, h, length, 0.)):
        a, b, c, e = (x0, 0., z0), (x1, 0., z1), (x1, width, z1), (x0, width, z0)
        tris += [(a, b, c), (a, c, e)]
    with open(path, "w") as f:
        f.write("solid ridge\n")
        for t in tris:
            n = np.cross(np.subtract(t[1], t[0]), np.subtract(t[2], t[0]))
            n = n / np.linalg.norm(n)
            if n[2] < 0:
                t, n = (t[0], t[2], t[1]), -n
            f.write(f"facet normal {n[0]:.6f} {n[1]:.6f} {n[2]:.6f}\n outer loop\n")
            for v in t:
                f.write(f"  vertex {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
            f.write(" endloop\nendfacet\n")
        f.write("endsolid ridge\n")


def ibm_preprocessor_available() -> bool:
    """The legacy IBM preprocessor is compiled with gfortran on first use (udprep_ibm)."""
    return shutil.which("gfortran") is not None


def write_plane_stl(path: Path, z: float, length: float, width: float) -> None:
    """ASCII STL of a horizontal plane at height z covering the domain (normal +z)."""
    with open(path, "w") as f:
        f.write("solid plane\n")
        for t in (((0., 0., z), (length, 0., z), (length, width, z)), ((0., 0., z), (length, width, z), (0., width, z))):
            f.write("facet normal 0.000000 0.000000 1.000000\n outer loop\n")
            for v in t:
                f.write(f"  vertex {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
            f.write(" endloop\nendfacet\n")
        f.write("endsolid plane\n")


def prepare_ridge(parent: Path, angle_deg: float = 20.0, z0: float = 2.e-4, dz: float = 20.0,
                  ktot: int = 64) -> Path:
    """The ridge of write_ridge_stl through prepare_stl (z0 = z0h)."""
    xlen = float(BASE["DOMAIN"]["xlen"])
    ylen = float(BASE["DOMAIN"]["ylen"])
    return prepare_stl(parent, lambda p: write_ridge_stl(p, angle_deg, xlen, ylen), z0, z0, dz, ktot)


def prepare_stl(parent: Path, write_stl, z0: float, z0h: float, dz: float = 20.0, ktot: int = 64) -> Path:
    """Run the normal Python IBM preprocessing (udprep, legacy backend) on the STL that
    write_stl(path) writes.

    The case directory must be named after the experiment number (UDPrep checks it), so it is
    parent/EXPNR. Returns that directory, holding the IBM input files, the facet types
    (one floor type with z0, z0h and lGR = 1) and the STL; write_ridge_run copies it."""
    from udprep import UDPrep
    d = parent / EXPNR
    write_case(d, ktot=ktot, dz=dz, groups={
        "RUN": {"libm": ".true."}, "INPS": {"stl_file": f"geom.{EXPNR}.stl", "zsize": str(ktot * dz)},
        "WALLS": {"lbottom": ".false.", "iwallmom": "4", "iwalltemp": "1", "iwallmoist": "1"}})
    write_stl(d / f"geom.{EXPNR}.stl")
    write_floor_factypes(d, EXPNR, z0, z0h)
    prep = UDPrep(EXPNR, d)
    prep.ibm.run_ibm(backend="legacy")
    prep.ibm.write_facets()
    prep.ibm.write_facetarea()
    # the floor facet type of factypes.inp (negative = ground)
    fp = d / f"facets.inp.{EXPNR}"
    lines = fp.read_text().splitlines()
    fp.write_text("\n".join([lines[0]] + [f"{-1:5d} " + l.split(None, 1)[1] for l in lines[1:]]) + "\n")
    counts = {k: getattr(prep.sim, k) for k in (
        "nfcts", "nsolpts_u", "nsolpts_v", "nsolpts_w", "nsolpts_c", "nbndpts_u", "nbndpts_v",
        "nbndpts_w", "nbndpts_c", "nfctsecs_u", "nfctsecs_v", "nfctsecs_w", "nfctsecs_c")}
    (d / "counts.txt").write_text("\n".join(f"{k} {v}" for k, v in counts.items()) + "\n")
    return d


def write_ridge_run(prep_dir: Path, run_dir: Path, groups: dict | None = None, ktot: int = 64,
                    dz: float = 20.0, **kw) -> Path:
    """Copy a prepared ridge case to run_dir and write its namoptions (with the IBM counts)."""
    shutil.copytree(prep_dir, run_dir)
    counts = {k: v for k, v in (l.split() for l in (run_dir / "counts.txt").read_text().splitlines())}
    g = {"RUN": {"libm": ".true."},
         "WALLS": {"lbottom": ".false.", "iwallmom": "4", "iwalltemp": "1", "iwallmoist": "1", **counts}}
    for name, kv in (groups or {}).items():
        g.setdefault(name, {}).update(kv)
    return write_case(run_dir, groups=g, ktot=ktot, dz=dz, **kw)


def restart_set(run_dir: Path, ntrun: int) -> list:
    """The restart files a run wrote at step ntrun (initd/inits/initm per rank, initf)."""
    tag = f"{ntrun:08d}"
    return sorted(p for p in run_dir.glob(f"init?{tag}*.{EXPNR}") if p.name[4] in "dsmf")


def stage_warm_start(cold_dir: Path, warm_dir: Path, ntrun: int, runtime: float, skip=()) -> Path:
    """Copy cold_dir to warm_dir as a warm start from the checkpoint at step ntrun: the case inputs,
    the restart set of that step (minus the basenames in skip) and namoptions with lwarmstart,
    startfile and runtime; outputs and the other restart files are left out."""
    warm_dir.mkdir(parents=True)
    for p in cold_dir.iterdir():
        if p.is_file() and not p.name.startswith(("init", "bldump", "fielddump", "prof_restart")) \
                and not p.name.endswith(".nc"):
            shutil.copy(p, warm_dir)
    for p in restart_set(cold_dir, ntrun):
        if p.name not in skip:
            shutil.copy(p, warm_dir)
    nml = warm_dir / f"namoptions.{EXPNR}"
    text = nml.read_text()
    text = re.sub(r"^runtime\s*=.*$", f"runtime      = {runtime}", text, count=1, flags=re.M)
    text = text.replace("&RUN\n", f"&RUN\nlwarmstart   = .true.\nstartfile    = 'initd{ntrun:08d}_xxx_xxx.{EXPNR}'\n", 1)
    nml.write_text(text)
    return warm_dir


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


def drag_coefficients(z, L, z0m, z0h):
    """Cm, Cs of DALES dev 075ac4e calc_drag_coefficients (modsurface:mo_presc_tskin_update)."""
    lm = np.log(z / z0m) - psim(z / L) + psim(z0m / L)
    lh = np.log(z / z0h) - psih(z / L) + psih(z0h / L)
    return FKAR ** 2 / lm ** 2, FKAR ** 2 / (lm * lh)


def qsat_facet(T):
    """initfac:qsat, the facet saturation humidity (Bolton 1980) [kg/kg]."""
    gres = 611.00 * np.exp(17.27 * (T - 273.15) / (T - 35.85))
    return 0.62198 * 0.01 * gres / (1000 - 0.01 * gres)


def moist_flux(resa, qtair, qwall, hurel, resc=200.0, ress=50.0, cveg=0.8):
    """modibm:moist_flux, positive into the surface (evaporation only)."""
    return min(0.0, cveg * (qtair - qwall) / (resa + resc) + (1 - cveg) * (qtair - qwall * hurel) / (resa + ress))
