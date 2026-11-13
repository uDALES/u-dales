# Moist thermodynamics, MO surface and parametrised radiation — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port DALES's Monin–Obukhov surface layer, Deardorff moist buoyancy surface term and parametrised longwave radiation into uDALES, add a small diagnostics module, and validate on BOMEX, smoke and DYCOMS RF01.

**Architecture:** Each physics component is a separate, default-off Fortran module that follows DALES `dev` @ `075ac4e` line by line, with indexing translated to uDALES (`ib:ie, jb:je, kb:ke`, local indices, `ib = jb = kb = 1`). Validation has two tiers. Tier 1 is small Python-driven solver runs in `tests/integration/thermo/` that compare against independent NumPy implementations. Tier 2 is full-size cases in `examples/301`–`303`, run on CX3 and compared with DALES output and published ranges.

**Tech Stack:** Fortran 2008 + MPI + 2DECOMP&FFT + netCDF-Fortran (solver); Python 3.9 + NumPy + netCDF4 (tests, tooling); CMake; PBS on CX3.

**Spec:** `docs/superpowers/specs/2026-11-13-moist-thermo-port-design.md` (read it first; this plan argues from it).

## Global Constraints

- Worktree `~/udales/u-dales-thermo`, branch `thermo`, based on `origin/master` @ `752430c0`. Run every command from the worktree root.
- DALES reference source: `dev` @ `075ac4e`. A shallow clone is in the session scratchpad; if it is missing: `git clone --depth 1 -b dev https://github.com/dalesteam/dales.git $SCRATCH/dales`. Cite `DALES dev 075ac4e <file>:<routine>` in a comment above each ported routine.
- Everything new is off by default. With `iradiation = 0`, `lthermostat = .false.` and `BCbotm ∉ {4,5}`, results must be bit-for-bit identical to the parent commit. Exceptions: the #375 commit (Task 1), and the `wqsurf` sign commit (Task 3; no shipped case uses a non-zero `wqsurf`).
- **Flux sign convention (decision D1, decided by Maarten 2026-11-13):** every prescribed boundary flux in uDALES is **positive in −z** (downward): `wtsurf`, `wttop`, `wqtop`, facet `bctf*` already are; `wqsurf` was the exception and is flipped in Task 3. Physics code that needs DALES-style upward fluxes uses `-wtsurf`, `-wqsurf`. The long-term convention is uDALES issue #393 (milestone 3.0); reference it in the Task 3 sign commit message.
- **All netCDF output is single precision** (`NF90_FLOAT`, `modstat_nc.f90:251-255`). Test tolerances must allow float32 (≈6e-8 relative; 3e-5 K at 290 K).
- Review: Fable reviewed v1 of this plan (findings B1–B3, M1–M5 and minors are folded in below).
- Namelist variable names follow DALES where a DALES equivalent exists.
- Fortran style: match the surrounding uDALES code (lower-case keywords, `use ..., only :` lists, `write(0, *) 'ERROR: ...'` + `stop 1` for fatal errors, MPI via `comm3d`/`MY_REAL`/`mpierr` from `modmpi`).
- Sources are globbed by CMake (`CMakeLists.txt:103`), so new `src/*.f90` files need no CMake edit, only a re-configure.
- Stage files by explicit path; run `git status` before every commit; never `git add -A`. Do not push. The controller pushes, once, at the end.
- End every commit message with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` (or the implementing model's own line).

## Build and test commands (CX3)

These are used by every task. Load modules in the same command as the build (module state does not persist; never pipe `module load`).

```bash
# build (Debug for development, Release for tier-2 runs)
rm -f *.mod *__genmod.f90                 # stale in-source .mod files shadow the build
./tools/build_executable.sh icl debug     # -> build/debug/u-dales
./tools/build_executable.sh icl release   # -> build/release/u-dales

# Python env for tests
module purge; module load tools/prod Python/3.9.6-GCCcore-11.2.0; source ~/udales/.venv/bin/activate

# netCDF-writing solver tests must run with the BUILD stack, not the default runtime stack
export UDALES_RUNTIME_MODULES="intel/2021a netCDF/4.8.0-iimpi-2021a netCDF-Fortran/4.5.3-iimpi-2021a FFTW/3.3.9-intel-2021a"
export UDALES_BUILD=$PWD/build/debug/u-dales

python tests/integration/thermo/test_thermo.py -v          # this feature's tests
python tests/run_tests.py python-library                   # post-build contract (AGENTS.md)
python tests/run_tests.py supported --branch-a master --branch-b HEAD --build-type Debug
```

`netCDF4` must import in `~/udales/.venv`. Check with `python -c "import netCDF4"`; if it is missing, `pip install netCDF4` and note it in the commit message.

## File map

| File | Status | Responsibility |
|---|---|---|
| `src/modstartup.f90` | modify | #375 checks; `&BC` gains `ustin`; broadcast |
| `src/modsubgrid.f90` | modify | guard length-scale `sqrt` (#375) |
| `src/modglobal.f90` | modify | `BCbotm_mo = 4`, `BCbotm_mo_ustar = 5` |
| `src/modsurfdata.f90` | modify | `ustin` |
| `src/modsurface.f90` | **create** | Monin–Obukhov surface layer (DALES `isurf 3/4`) |
| `src/modibm.f90` | modify | `bottom`: dispatch `BCbotm 4/5`; flux sign (D1) |
| `src/modfields.f90` | modify | `tmp0` |
| `src/modthermodynamics.f90` | modify | `calthv` surface level + `tmp0` |
| `src/modfielddump.f90` | modify | `dt` → `dthvdz` dump code |
| `src/modraddata.f90` | **create** | radiation switches, parameters, fields |
| `src/modradiation.f90` | **create** | `initradiation`, `radiation`, `radpar`, `exitradiation` |
| `src/modthermostat.f90` | **create** | time series + slab profiles netCDF |
| `src/program.f90` | modify | init / step / exit calls |
| `tools/python/namelists.json` | modify | `NAMRADIATION`, `NAMTHERMOSTAT`, `ustin` |
| `docs/udales-namoptions-overview.md`, `docs/udales-boundary-conditions.md` | modify | document new options |
| `tests/integration/thermo/thermo_case.py` | **create** | programmatic case writer for tests |
| `tests/integration/thermo/test_thermo.py` | **create** | tier-1 tests |
| `tests/integration/thermo/README.md` | **create** | what each test proves |
| `tests/test_suites.yml` | modify | register suite in `supported` |
| `tests/cases/*`, `examples/*`, `tests/integration/**` namoptions | modify | set `thls`/`qts` (#375) |
| `examples/301`, `302`, `303` | **create** | BOMEX, smoke, RF01 full-size cases |
| `tools/python/examples/dales_cases/convert_dales_case.py` | **create** | DALES → uDALES input converter |
| `tools/python/examples/dales_cases/make_rf01.py` | **create** | RF01 input generator |
| `tools/python/examples/dales_cases/compare_dales.py` | **create** | tier-2 comparison plots + table |

---

### Task 0: Baseline

**Files:** none modified.

- [ ] **Step 1: Build and record baseline**

```bash
cd ~/udales/u-dales-thermo && git status && git log --oneline -3
rm -f *.mod *__genmod.f90 && ./tools/build_executable.sh icl debug 2>&1 | tail -5
ls -l build/debug/u-dales
```
Expected: clean tree at `b3fe6a0c` or later, and the executable exists.

- [ ] **Step 2: Run the supported selection once and save the result**

```bash
module purge; module load tools/prod Python/3.9.6-GCCcore-11.2.0; source ~/udales/.venv/bin/activate
python tests/run_tests.py supported --branch-a master --branch-b HEAD --build-type Debug 2>&1 | tee $SCRATCH/baseline_supported.log | tail -30
```
Expected: record which suites pass or fail **before** any change. A failure here is pre-existing; note it and do not fix it on this branch.

---

### Task 1: Fix #375 — validate `thls`/`qts`, guard SGS length scale

**Files:**
- Modify: `src/modstartup.f90` (`checkinitvalues`, ~697–944; dead `! call checkinitvalues` ~635)
- Modify: `src/modsubgrid.f90:391`
- Modify: every namoptions listed in Step 4
- Test: `tests/integration/thermo/test_thermo.py::TestInitChecks`, plus `thermo_case.py` (created here)

**Interfaces:**
- Produces: `thermo_case.write_case(run_dir, **overrides) -> Path` and `thermo_case.run(run_dir, nprocs=1) -> subprocess.CompletedProcess`; all later tests use them.

- [ ] **Step 1: Create the test helper `tests/integration/thermo/thermo_case.py`**

It writes a flat, IBM-free, periodic 16×16×64 case with expnr 310. All tests build on it.

```python
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
import mpi_env  # noqa: E402

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
               thl=None, qt=None, u=7.0, v=-5.5, wfls=None, thlpcar=None) -> Path:
    """Write namoptions/prof.inp/lscale.inp for a flat case. Returns run_dir."""
    run_dir.mkdir(parents=True, exist_ok=True)
    nml = {g: dict(kv) for g, kv in BASE.items()}
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
    wfls = np.zeros(ktot) if wfls is None else wfls
    thlpcar = np.zeros(ktot) if thlpcar is None else thlpcar
    with open(run_dir / f"lscale.inp.{EXPNR}", "w") as f:
        f.write("# thermo test case\n# z ug vg pgx pgy wfls dqtdxls dqtdyls dqtdtls thlpcar\n")
        for k in range(ktot):
            f.write(f"{zf[k]:.6f} {u:.4f} {v:.4f} 0. 0. {wfls[k]:.6e} 0. 0. 0. {thlpcar[k]:.6e}\n")
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
```

Before relying on `mpi_env.mpi_exec_and_args()`, open `tests/integration/mpi_env.py` and match its real return type. `tests/integration/tke_top_bc/test_tke_top_bc.py` shows the call; copy whatever it does.

- [ ] **Step 2: Write the failing test `tests/integration/thermo/test_thermo.py`**

```python
"""Tier-1 tests for the moist thermodynamics / MO surface / radiation port."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

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


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 3: Run it and confirm it fails**

```bash
python tests/integration/thermo/test_thermo.py -v TestInitChecks
```
Expected: `test_sentinel_thls_aborts` FAILS, because the run succeeds or traps with an FPE that does not mention `thls`. `test_dry_case_without_qts_runs` may fail with an FPE, which is #375 itself.

- [ ] **Step 4: Implement**

In `checkinitvalues` (`src/modstartup.f90`), add `thls, qts` to its `use modsurfdata, only :` list and `ltempeq, lbuoyancy, lmoist` to its `use modglobal` list. After the existing `ps < 0` check, add:

```fortran
      ! #375: thls/qts default to the sentinel -1; a negative thvs flips the sign
      ! of SGS buoyancy production and of the hydrostatic base state.
      if ((ltempeq .or. lbuoyancy) .and. (thls < 0.)) then
         write(0, *) 'ERROR: thls not set (sentinel -1) but ltempeq or lbuoyancy is on.'
         write(0, *) '       Set thls in &BC to the surface liquid water potential temperature [K].'
         stop 1
      end if
      if (qts < 0.) then
         if (lmoist) then
            write(0, *) 'ERROR: qts not set (sentinel -1) but lmoist is on. Set qts in &BC [kg/kg].'
            stop 1
         else
            qts = 0.   ! dry run: surface humidity is irrelevant but enters thvs
         end if
      end if
```

`checkinitvalues` runs on every rank after the broadcast, so the `qts = 0.` assignment is rank-consistent. Confirm this by reading where `program.f90:74` sits relative to the `MPI_BCAST(qts, ...)` in `readnamelists`. Then `thvs` must be recomputed: the existing `thvs = thls*(1.+(rv/rd - 1.)*qts)` at `modstartup.f90:525` runs **before** `checkinitvalues`. So add `thvs = thls*(1.+(rv/rd - 1.)*qts)` directly after the `qts = 0.` line, and add `thvs` to the `use modsurfdata` list.

Delete the dead commented-out `! call checkinitvalues` block at `modstartup.f90:~635`.

In `src/modsubgrid.f90:391`, replace the unguarded line with:

```fortran
                   zlt(i,j,k) = min(delta(i,k),cn*e120(i,j,k)/sqrt(grav/thvs*max(abs(dthvdz(i,j,k)),eps1)))
```
and add `eps1` to that routine's `use modglobal` list if it is not already there. Leave the commented line above it untouched.

Set `thls`/`qts` in every shipped namoptions that enables `ltempeq` or `lbuoyancy` without them:

```
examples/102/namoptions.102  examples/102/namoptions_old.102  examples/201/namoptions.201 (thls only)
tests/cases/064/namoptions.064  tests/cases/300/namoptions.300  tests/cases/525/namoptions.525
tests/cases/526/namoptions.526  tests/integration/processor_boundaries/namoptions.526.serial
tests/integration/sgs_statistics/namoptions.1007.{serial,xsplit,ysplit,xysplit}
tests/integration/warmstart_facets/namoptions.064.template
tests/regression/mpi_averaging_regression/namelists/namoptions.526.serial.{current,reference}
tests/regression/new_vegetation_module_against_v2.2/namoptions.526.{current,reference}
```
Rules:
- `thls` is the first-row `thl` of the case's `prof.inp`.
- `qts` is the first-row `qt` if `lmoist`, else leave it out (it defaults to 0 now).
- Put both in the existing `&BC` group, or add an `&BC` group if the case has none.

Re-run the grep from the spec's survey to confirm none remain:

```bash
for f in $(grep -rlE '^\s*(ltempeq|lbuoyancy)\s*=\s*\.true\.' examples tests --include='namoptions*'); do grep -qE '^\s*thls\s*=' $f || echo "MISSING thls: $f"; done
```
Expected: no output.

- [ ] **Step 5: Build, run the tests, run supported**

```bash
./tools/build_executable.sh icl debug 2>&1 | tail -3
python tests/integration/thermo/test_thermo.py -v TestInitChecks
python tests/run_tests.py supported --branch-a master --branch-b HEAD --build-type Debug 2>&1 | tail -30
```
Expected: both new tests PASS. Compare supported against the baseline. The branch-comparison regression will now differ on **essentially every `ltempeq` case**, which is expected for this commit. Before, a dry case with sentinel `qts = -1` had `thvs = 0.392·thls`; `qts = 0` multiplies `thvs` by 2.55. `thvs` feeds `thvh(kb)` in `diagfld` and the SGS length scale, and cases whose `thls` changed differ too. List the differing cases in the commit message with that explanation. No *other* suite may newly fail.

- [ ] **Step 6: Commit**

```bash
git add src/modstartup.f90 src/modsubgrid.f90 tests/integration/thermo/thermo_case.py tests/integration/thermo/test_thermo.py <each namoptions edited>
git status
git commit -m "Validate thls/qts sentinels and guard SGS length scale (fixes #375)

checkinitvalues aborts when thls is unset with ltempeq/lbuoyancy, or qts
is unset with lmoist; dry runs default qts to 0 and recompute thvs.
Shipped cases now set thls/qts. <list cases whose regression output changed>

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Diagnostics module `modthermostat`

**Files:**
- Create: `src/modthermostat.f90`
- Modify: `src/program.f90` (init after `initfielddump`; call after `statsdump` in loop; exit after `exitstatsdump`)
- Modify: `tools/python/namelists.json`
- Test: `test_thermo.py::TestThermostat`

**Interfaces:**
- Produces: files `thermostat_t.310.nc` (time series) and `thermostat_z.310.nc` (profiles), written by rank 0. Later tasks append variables by adding entries to the `nvt`/`nvz` tables and filling them in `thermostat`. Module routines: `initthermostat`, `thermostat`, `exitthermostat`.
- Time series variables (`time` dim): `lwp_bar` [kg m⁻²], `cfrac` [-], `zi` [m] (column-mean height of max `dthl/dz`), `zi_qt8` [m] (column-mean height of the 8 g/kg `qt` crossing; 0 if dry), `zi_smoke` [m] (column-mean 0.5 crossing of scalar `isvsmoke`; 0 if `nsv = 0`), `wmax` [m s⁻¹], `qlmax` [kg kg⁻¹], `ustar` [m s⁻¹], `obl` [m], `thls` [K], `qts` [kg kg⁻¹].
- Profile variables (`zt`,`time`; code `tt`): `thl`, `qt`, `ql`, `cfracz`, `rhof`, `exnf`, `thlprad`. Half-level (`zm`,`time`; code `mt`): `wthl`, `wqt`, `lwu`, `lwd`.
- Until Tasks 3 and 5 exist, `ustar`/`obl` write 0 and `thlprad`/`lwu`/`lwd` write 0. Those tasks replace the zeros.

- [ ] **Step 1: Write the failing test** (append to `test_thermo.py`)

```python
from netCDF4 import Dataset


class TestThermostat(_SolverTest):
    def test_cloudy_column_diagnostics(self):
        # 3 steps, laminar RF01 column: horizontally homogeneous, so column = slab mean
        d = tc.write_case(self.tmp / "stat", groups={
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
```

- [ ] **Step 2: Run it to verify it fails**

`python tests/integration/thermo/test_thermo.py -v TestThermostat`. Expected: FAIL, because `thermostat_t.310.nc` is missing. The namelist group `NAMTHERMOSTAT` is ignored if absent: the group-read pattern tolerates EOF.

- [ ] **Step 3: Implement `src/modthermostat.f90`**

```fortran
!> Cloud, surface and radiation diagnostics: domain-mean time series and slab
!> profiles, the uDALES counterpart of DALES tmser/profiles. Kept separate from
!> modstatsdump so it can move into the statsout framework as one unit.
module modthermostat
  implicit none
  save
  private
  public :: initthermostat, thermostat, exitthermostat, lthermostat, tthermostat

  logical :: lthermostat = .false.
  real    :: tthermostat = 60.    !< output interval [s]
  real    :: tnext = 0.
  integer :: ncidt, ncidz, nrect = 0, nrecz = 0
  integer, parameter :: nvt = 11, nvz = 11
  character(80) :: fnamet = 'thermostat_t.xxx.nc', fnamez = 'thermostat_z.xxx.nc'
  character(80) :: tnct(1,4), tncz(1,4), nct(nvt,4), ncz(nvz,4)

contains

  subroutine initthermostat
    use modglobal, only : ifnamopt, fname_options, cexpnr, kb, ke
    use modmpi,    only : myid, comm3d, mpierr, my_real, mpi_logical
    use modstat_nc, only : ncinfo, open_nc, define_nc, writestat_dims_nc
    implicit none
    integer :: ierr
    namelist/NAMTHERMOSTAT/ lthermostat, tthermostat

    if (myid == 0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read (ifnamopt, NAMTHERMOSTAT, iostat=ierr)
      if (ierr > 0) then
        write(0, *) 'ERROR: Problem in namoptions NAMTHERMOSTAT'
        write(0, *) 'iostat error: ', ierr
        stop 1
      end if
      close(ifnamopt)
    end if
    call MPI_BCAST(lthermostat, 1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(tthermostat, 1, MY_REAL,     0, comm3d, mpierr)
    if (.not. lthermostat) return

    tnext = tthermostat
    fnamet(14:16) = cexpnr
    fnamez(14:16) = cexpnr
    call ncinfo(tnct(1,:), 'time', 'Time', 's', 'time')
    tncz = tnct
    call ncinfo(nct( 1,:), 'lwp_bar',  'Liquid water path',             'kg/m^2', 'time')
    call ncinfo(nct( 2,:), 'cfrac',    'Cloud cover',                   '-',      'time')
    call ncinfo(nct( 3,:), 'zi',       'BL height (max dthl/dz)',       'm',      'time')
    call ncinfo(nct( 4,:), 'zi_qt8',   'BL height (qt = 8 g/kg)',       'm',      'time')
    call ncinfo(nct( 5,:), 'zi_smoke', 'BL height (smoke = 0.5)',       'm',      'time')
    call ncinfo(nct( 6,:), 'wmax',     'Maximum vertical velocity',     'm/s',    'time')
    call ncinfo(nct( 7,:), 'qlmax',    'Maximum liquid water',          'kg/kg',  'time')
    call ncinfo(nct( 8,:), 'ustar',    'Friction velocity',             'm/s',    'time')
    call ncinfo(nct( 9,:), 'obl',      'Obukhov length',                'm',      'time')
    call ncinfo(nct(10,:), 'thls',     'Surface thl',                   'K',      'time')
    call ncinfo(nct(11,:), 'qts',      'Surface qt',                    'kg/kg',  'time')
    call ncinfo(ncz( 1,:), 'thl',     'Liquid water potential temperature', 'K',      'tt')
    call ncinfo(ncz( 2,:), 'qt',      'Total water',                        'kg/kg',  'tt')
    call ncinfo(ncz( 3,:), 'ql',      'Liquid water',                       'kg/kg',  'tt')
    call ncinfo(ncz( 4,:), 'cfracz',  'Cloud fraction',                     '-',      'tt')
    call ncinfo(ncz( 5,:), 'rhof',    'Density',                            'kg/m^3', 'tt')
    call ncinfo(ncz( 6,:), 'exnf',    'Exner function',                     '-',      'tt')
    call ncinfo(ncz( 7,:), 'thlprad', 'Radiative thl tendency',             'K/s',    'tt')
    call ncinfo(ncz( 8,:), 'wthl',    'Resolved w thl flux',                'K m/s',  'mt')
    call ncinfo(ncz( 9,:), 'wqt',     'Resolved w qt flux',                 'm/s',    'mt')
    call ncinfo(ncz(10,:), 'lwu',     'Longwave flux, F1 term',             'W/m^2',  'mt')
    call ncinfo(ncz(11,:), 'lwd',     'Longwave flux, F0 term',             'W/m^2',  'mt')
    if (myid == 0) then
      call open_nc(fnamet, ncidt, nrect)
      if (nrect == 0) call define_nc(ncidt, 1, tnct)
      call define_nc(ncidt, nvt, nct)
      call open_nc(fnamez, ncidz, nrecz, n3=ke-kb+1)
      if (nrecz == 0) then
        call define_nc(ncidz, 1, tncz)
        call writestat_dims_nc(ncidz)      ! fill zt/zm, as modstatsdump.f90:219-222
      end if
      call define_nc(ncidz, nvz, ncz)
    end if
  end subroutine initthermostat

  subroutine thermostat
    use modglobal,  only : ib, ie, jb, je, kb, ke, itot, jtot, timee, rk3step, zf, zh, dzf, dzh, &
                           nsv, lmoist
    use modfields,  only : thl0, qt0, ql0, w0, sv0, rhof, exnf, IIc, IIcs, IIw, IIws
    use modsurfdata, only : thls, qts
    use modmpi,     only : myid, comm3d, mpierr, my_real, mpi_sum, mpi_max, avexy_ibm
    use modstat_nc, only : writestat_nc
    implicit none
    real :: vt(nvt), vz(kb:ke, nvz)
    real, allocatable :: tmp3(:,:,:)   ! allocatable: an automatic array is ~22 MB/rank on case 303
    real :: loc(7), glob(7), colwp, gmax, dthl
    integer :: i, j, k, kmx
    real :: ncol

    if (.not. lthermostat) return
    if (rk3step /= 3 .or. timee < tnext) return
    tnext = tnext + tthermostat
    ncol = real(itot) * real(jtot)
    allocate(tmp3(ib:ie,jb:je,kb:ke+1))

    ! --- column diagnostics: lwp, cloud cover, three BL heights
    loc = 0.
    do j = jb, je
      do i = ib, ie
        colwp = 0.
        do k = kb, ke
          colwp = colwp + rhof(k) * ql0(i,j,k) * dzf(k) * IIc(i,j,k)
        end do
        loc(1) = loc(1) + colwp
        if (any(ql0(i,j,kb:ke) * IIc(i,j,kb:ke) > 0.)) loc(2) = loc(2) + 1.
        gmax = -huge(1.); kmx = kb + 1
        do k = kb + 1, ke
          dthl = (thl0(i,j,k) - thl0(i,j,k-1)) / dzh(k)
          if (dthl > gmax) then
            gmax = dthl; kmx = k
          end if
        end do
        loc(3) = loc(3) + zh(kmx)
        if (lmoist) loc(4) = loc(4) + crossing(qt0(i,j,kb:ke), 8.e-3)
        if (nsv > 0) loc(5) = loc(5) + crossing(sv0(i,j,kb:ke,1), 0.5)   ! Task 5: index -> isvsmoke
      end do
    end do
    loc(6) = maxval(w0(ib:ie,jb:je,kb:ke))
    loc(7) = maxval(ql0(ib:ie,jb:je,kb:ke))
    call MPI_ALLREDUCE(loc(1:5), glob(1:5), 5, MY_REAL, MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(loc(6:7), glob(6:7), 2, MY_REAL, MPI_MAX, comm3d, mpierr)
    vt(1:5) = glob(1:5) / ncol
    vt(6:7) = glob(6:7)
    vt(8)  = 0.          ! ustar: Task 3
    vt(9)  = 0.          ! obl:   Task 3
    vt(10) = thls
    vt(11) = qts

    ! --- slab profiles (fluid cells only)
    call avexy_ibm(vz(kb:ke,1), thl0(ib:ie,jb:je,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIc(ib:ie,jb:je,kb:ke), IIcs(kb:ke), .false.)
    call avexy_ibm(vz(kb:ke,2), qt0 (ib:ie,jb:je,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIc(ib:ie,jb:je,kb:ke), IIcs(kb:ke), .false.)
    call avexy_ibm(vz(kb:ke,3), ql0 (ib:ie,jb:je,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIc(ib:ie,jb:je,kb:ke), IIcs(kb:ke), .false.)
    tmp3(:,:,kb:ke) = merge(1., 0., ql0(ib:ie,jb:je,kb:ke) > 0.)
    call avexy_ibm(vz(kb:ke,4), tmp3(:,:,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIc(ib:ie,jb:je,kb:ke), IIcs(kb:ke), .false.)
    vz(:,5) = rhof(kb:ke)
    vz(:,6) = exnf(kb:ke)
    vz(:,7) = 0.         ! thlprad: Task 5
    ! resolved fluxes at half levels k (face between k-1 and k); kb face is 0
    do k = kb + 1, ke
      tmp3(:,:,k) = w0(ib:ie,jb:je,k) * 0.5 * (thl0(ib:ie,jb:je,k) + thl0(ib:ie,jb:je,k-1))
    end do
    tmp3(:,:,kb) = 0.
    call avexy_ibm(vz(kb:ke,8), tmp3(:,:,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIw(ib:ie,jb:je,kb:ke), IIws(kb:ke), .false.)
    do k = kb + 1, ke
      tmp3(:,:,k) = w0(ib:ie,jb:je,k) * 0.5 * (qt0(ib:ie,jb:je,k) + qt0(ib:ie,jb:je,k-1))
    end do
    tmp3(:,:,kb) = 0.
    call avexy_ibm(vz(kb:ke,9), tmp3(:,:,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIw(ib:ie,jb:je,kb:ke), IIws(kb:ke), .false.)
    vz(:,10) = 0.        ! lwu: Task 5
    vz(:,11) = 0.        ! lwd: Task 5

    if (myid == 0) then
      call writestat_nc(ncidt, 1, tnct, (/timee/), nrect, .true.)
      call writestat_nc(ncidt, nvt, nct, vt, nrect, .false.)
      call writestat_nc(ncidz, 1, tncz, (/timee/), nrecz, .true.)
      call writestat_nc(ncidz, nvz, ncz, vz, nrecz, ke-kb+1)
    end if
    deallocate(tmp3)

  contains

    !> Height of the highest downward crossing of `thres` in a column; 0 if none.
    real function crossing(col, thres)
      real, intent(in) :: col(kb:ke), thres
      integer :: kk
      crossing = 0.
      do kk = ke - 1, kb, -1
        if (col(kk) > thres .and. col(kk+1) <= thres) then
          crossing = zf(kk) + (thres - col(kk)) * (zf(kk+1) - zf(kk)) / (col(kk+1) - col(kk))
          return
        end if
      end do
    end function crossing

  end subroutine thermostat

  subroutine exitthermostat
    use modmpi, only : myid
    use modstat_nc, only : exitstat_nc
    implicit none
    if (.not. lthermostat) return
    if (myid == 0) then
      call exitstat_nc(ncidt)
      call exitstat_nc(ncidz)
    end if
  end subroutine exitthermostat

end module modthermostat
```

Verification notes for the implementer. Check each against the actual code; do not assume.

- **`avexy_ibm` arguments.** It expects `aver(kb:ke+kh)` and `var(ib:ie,jb:je,kb:ke+kh)` with an explicit `kh`. The calls above pass `ke-1, 1` so that `ke+kh` resolves to `ke`, which keeps array sections exact. If you prefer, pass sections `kb:ke+kh` of full-size work arrays instead. Either is fine, as long as the bounds match the dummy arguments exactly.
- **Time-series writes.** Verified by review: a 1-D `vt(nvt)` resolves to `writestat_time_nc` and a 2-D `vz` to `writestat_1D_nc`.
- **Output timing.** `thermostat`, like `fielddump`, runs at `rk3step == 3` **after** `tstep_integrate` and **before** `thermodynamics`. Prognostic fields in a record are end-of-step. Diagnosed fields (`ql0`, `rhof`, `exnf`, `dthvdz`) are from the previous `thermodynamics` call, i.e. start-of-step. Tests below are written for this.
- **Profile writes.** `writestat_1D_nc` expects `vars(dim1,nvar)`, which `vz(kb:ke,nvz)` matches.
- **Filename slices.** `fnamet(14:16)` assumes `'thermostat_t.'` is 13 characters, which it is. Verify the slice anyway.

In `program.f90`:
- add `use modthermostat, only : initthermostat, thermostat, exitthermostat`;
- `call initthermostat` after `call initfielddump`;
- `call thermostat` after `call statsdump`;
- `call exitthermostat` after `call exitstatsdump`.

In `tools/python/namelists.json`, add `"NAMTHERMOSTAT": ["lthermostat", "tthermostat"]` under `namelists`, plus `"lthermostat": "NAMTHERMOSTAT"` and `"tthermostat": "NAMTHERMOSTAT"` under `variables`. Keep keys alphabetical, as the file does.

- [ ] **Step 4: Build and test**

```bash
./tools/build_executable.sh icl debug 2>&1 | tail -3
python tests/integration/thermo/test_thermo.py -v TestThermostat TestInitChecks
python tests/run_tests.py python-library 2>&1 | tail -5
```
Expected: PASS. If `zi` (max-gradient) is off by one cell, check the half-level indexing (`zh(k)` is the face below cell `k`). The RF01 jump sits between the cells straddling 840 m, so `zh` at that face is 840 m with `dz = 20 m`.

- [ ] **Step 5: Commit**

```bash
git add src/modthermostat.f90 src/program.f90 tools/python/namelists.json tests/integration/thermo/test_thermo.py
git status
git commit -m "Add modthermostat: cloud/surface/radiation diagnostics

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Surface flux sign (decision D1) and Monin–Obukhov surface layer

**Files:**
- Modify: `src/modibm.f90` (`bottom`, ~2063–2165)
- Modify: `src/modglobal.f90:163-176`, `src/modsurfdata.f90`, `src/modstartup.f90` (`&BC` list ~133–141; broadcast; `checkinitvalues`)
- Create: `src/modsurface.f90`
- Modify: `src/modthermostat.f90` (ustar/obl)
- Modify: `docs/udales-boundary-conditions.md`, `docs/udales-namoptions-overview.md`, `tools/python/namelists.json`
- Test: `test_thermo.py::TestSurfaceSign`, `::TestMO`

**Interfaces:**
- Consumes: `thermo_case`, thermostat variables `ustar`, `obl`, `thls`, `qts`, `thl`, `qt`.
- Produces, module `modsurface` (public): `real :: ustar_mo, obl_mo, dthldz_mo, dqtdz_mo`; `subroutine surface_mo` (called from `bottom`); `logical function lmo_surface()` returns `BCbotm == BCbotm_mo .or. BCbotm == BCbotm_mo_ustar`. Task 4 uses `dthldz_mo`, `dqtdz_mo`, `lmo_surface`.
- `modglobal`: `integer, parameter :: BCbotm_mo = 4, BCbotm_mo_ustar = 5`.
- `modsurfdata`: `real :: ustin = -1.` (`&BC`).

- [ ] **Step 1: Write the sign test** (append to `test_thermo.py`)

```python
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
```

The thermostat profile is written at the end of step 1. If the record timing differs (see the note in Task 2), use the last record and keep the rest of the test as is.

- [ ] **Step 2: Run it, record the outcome, apply decision D1**

`python tests/integration/thermo/test_thermo.py -v TestSurfaceSign`

Expected, from reading the code: `dthl` PASSES and `dqt` FAILS with the opposite sign. That confirms `wqsurf` is the one upward-positive flux.

D1 (decided: flip `wqsurf`). In the `BCbotq.eq.1` branch of `bottom` replace `+ wqsurf &` with `- wqsurf &`. Add one comment line above the `BCbotT` branch: `! wtsurf, wqsurf (like wttop, wqtop, bctf*): prescribed fluxes, positive in -z (downward)`. Then document the convention:
- `docs/udales-namoptions-overview.md:143-144`: "Kinematic heat flux at domain bottom, positive downward (into the surface)" and the same for moisture; same wording for `wttop`/`wqtop`.
- `docs/udales-boundary-conditions.md`: one sentence stating the global rule.
- `docs/udales-example-simulations.md:~452`: that `wtsurf = 0.01` has no effect in case 102 (no `lbottom`), so say so or remove the line.

If the test shows a different starting state than expected, stop and report to the controller. Do not improvise a sign change.

Re-run: both assertions PASS. Commit this alone:

```bash
git add src/modibm.f90 docs/udales-namoptions-overview.md docs/udales-boundary-conditions.md docs/udales-example-simulations.md tests/integration/thermo/test_thermo.py
git status
git commit -m "Make wqsurf downward-positive like all other prescribed fluxes

wtsurf, wttop, wqtop and facet fluxes are positive in -z; wqsurf alone
was upward-positive. No shipped case uses a non-zero wqsurf.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

- [ ] **Step 3: Write the failing MO tests**

The NumPy reference is a direct port of DALES `calc_obl_iter`, `psim`, `psih`. Put it in `thermo_case.py`:

```python
FKAR, GRAV, RD, RV = 0.41, 9.81, 287.04, 461.5


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
```

Verify `fkar` (von Kármán) in `modglobal`; use its value instead of 0.41 if it differs.

Tests:

```python
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
        self.assertGreater(ustar, 0.01)
        self.assertLess(obl, 0.0)                         # upward buoyancy flux -> unstable
        self.assertGreater(thls, thl1)                    # surface warmer than air
        self.assertAlmostEqual(obl, obl_prev, delta=0.02 * abs(obl))   # converged
        # Reference: same iteration on the recorded converged state (float32 records).
        z1 = 10.0
        L = tc.obl_iter(thl1, qt1, thls, qts, z1, 2e-4, 2e-4, 7.0, -5.5, -0.01)
        self.assertAlmostEqual(obl, L, delta=0.02 * abs(L))
        cm = tc.FKAR ** 2 / (np.log(z1 / 2e-4) - tc.psim(z1 / L) + tc.psim(2e-4 / L)) ** 2
        self.assertAlmostEqual(ustar, np.sqrt(cm) * np.hypot(7.0, -5.5), delta=0.02 * ustar)

    def test_prescribed_ustar(self):
        d = self._mo_case("mo5", 5, {"ustin": "0.28"})
        with Dataset(d / "thermostat_t.310.nc") as t:
            self.assertAlmostEqual(float(t["ustar"][-1]), 0.28, places=6)   # float32 record
```

The mean wind stays at the initial (7, −5.5) to within the 2 % tolerance over 20 s because the laminar column has no turbulence (`lrandomize = .false.`) and only surface drag acts. If the test shows otherwise, take `u0av`/`v0av` from a new thermostat variable rather than loosening the tolerance.

- [ ] **Step 4: Run to verify failure**

`python tests/integration/thermo/test_thermo.py -v TestMO`. Expected: FAIL, because `BCbotm = 4` stops with "bottom boundary type for momentum undefined".

- [ ] **Step 5: Implement `src/modsurface.f90`**

```fortran
!> Monin-Obukhov surface layer for the flat bottom (lbottom), ported from
!> DALES dev 075ac4e src/surface/modsurface.f90 (isurf = 3/4, lmostlocal = .false.):
!> domain-mean Obukhov length and friction velocity, locally applied stress.
!>   BCbotm = 4 (BCbotm_mo):       u* from MO similarity            (DALES isurf = 4)
!>   BCbotm = 5 (BCbotm_mo_ustar): u* = ustin                       (DALES isurf = 3)
!> Heat and moisture fluxes stay prescribed by wtsurf/wqsurf (BCbotT/BCbotq = 1).
module modsurface
  implicit none
  save
  private
  public :: surface_mo, lmo_surface, ustar_mo, obl_mo, dthldz_mo, dqtdz_mo

  real :: ustar_mo  = 0.      !< friction velocity [m/s]
  real :: obl_mo    = -0.1    !< Obukhov length [m] (also the iteration's first guess)
  real :: dthldz_mo = 0.      !< MO surface-layer thl gradient at zf(kb) [K/m]
  real :: dqtdz_mo  = 0.      !< MO surface-layer qt gradient at zf(kb) [1/m]
  real, parameter :: min_horv = 0.1   !< minimum wind speed [m/s] (DALES)

contains

  logical function lmo_surface()
    use modglobal, only : BCbotm, BCbotm_mo, BCbotm_mo_ustar
    lmo_surface = (BCbotm == BCbotm_mo) .or. (BCbotm == BCbotm_mo_ustar)
  end function lmo_surface

  subroutine surface_mo
    use modglobal,  only : ib, ie, jb, je, kb, ih, jh, kh, itot, jtot, zf, dzf, dzfi, dzhi, dzhiq, &
                           dxf, dxhi, fkar, rv, rd, lmoist, BCbotm, BCbotm_mo_ustar
    use modfields,  only : u0, v0, up, vp, thl0, qt0, u0av, v0av, thl0av, qt0av, momfluxb
    use modsurfdata, only : wtsurf, wqsurf, thls, qts, thvs, z0, z0h, ustin
    use modsubgriddata, only : ekm
    use modmpi,     only : comm3d, mpierr, my_real, mpi_sum
    implicit none
    integer :: i, j
    real :: z1, horvav, cm, cs, scaling, upcu, vpcv, horv, emom, utan2, bcmom
    real :: loc(2), glob(2), wthl_up, wqt_up

    z1 = zf(kb)
    ! DALES upward-positive kinematic fluxes; uDALES prescribes them positive in -z (D1)
    wthl_up = -wtsurf
    wqt_up  = -wqsurf

    ! drag coefficients (DALES calc_drag_coefficients, uniform here). DALES order for
    ! isurf 3/4: drag coefficients BEFORE get_obl, i.e. with the previous obl.
    cm = fkar**2 / (log(z1/z0) - psim(z1/obl_mo) + psim(z0/obl_mo))**2
    cs = fkar**2 / (log(z1/z0) - psim(z1/obl_mo) + psim(z0/obl_mo)) &
                 / (log(z1/z0h) - psih(z1/obl_mo) + psih(z0h/obl_mo))

    ! Obukhov length from the slab-mean first-level state and current surface values
    obl_mo = calc_obl_iter(thl0av(kb), qt0av(kb), thls, qts, z1, z0, z0h, u0av(kb), v0av(kb), obl_mo)

    ! friction velocity (DALES calc_friction_velocity / presc_friction_velocity)
    horvav = max(sqrt(u0av(kb)**2 + v0av(kb)**2), min_horv)
    if (BCbotm == BCbotm_mo_ustar) then
      ustar_mo = ustin
    else
      ustar_mo = sqrt(cm) * horvav
    end if
    ustar_mo = max(ustar_mo, 1.e-2)

    ! surface-layer gradients (DALES calc_surface_gradients; fluxes upward-positive)
    scaling   = 1. / (fkar * z1)
    dthldz_mo = - wthl_up / ustar_mo * phih(z1/obl_mo) * scaling
    if (lmoist) then
      dqtdz_mo = - wqt_up / ustar_mo * phih(z1/obl_mo) * scaling
    else
      dqtdz_mo = 0.
    end if

    ! surface scalars (DALES calc_surface_scalars): diagnose tskin/qskin per column, average
    loc = 0.
    do j = jb, je
      do i = ib, ie
        upcu = 0.5 * (u0(i,j,kb) + u0(i+1,j,kb))
        vpcv = 0.5 * (v0(i,j,kb) + v0(i,j+1,kb))
        horv = max(sqrt(upcu**2 + vpcv**2), min_horv)
        loc(1) = loc(1) + min(max(wthl_up / (cs * horv), -10.), 10.) + thl0(i,j,kb)
        loc(2) = loc(2) + min(max(wqt_up / (cs * horv), -5.e-2), 5.e-2) + qt0(i,j,kb)
      end do
    end do
    call MPI_ALLREDUCE(loc, glob, 2, MY_REAL, MPI_SUM, comm3d, mpierr)
    thls = glob(1) / (real(itot) * real(jtot))
    if (lmoist) qts = glob(2) / (real(itot) * real(jtot))
    thvs = thls * (1. + (rv/rd - 1.) * qts)

    ! momentum: replace the resolved diffusive flux through the bottom face by
    ! u'w' = -u*^2 u/|U| (pattern of modwallfunctions:wfmneutral, case 91)
    do j = jb, je
      do i = ib, ie
        utan2 = 0.25 * (v0(i,j,kb) + v0(i-1,j,kb) + v0(i,j+1,kb) + v0(i-1,j+1,kb))
        horv  = sqrt(max(min_horv**2, u0(i,j,kb)**2 + utan2**2))
        bcmom = ustar_mo**2 * u0(i,j,kb) / horv
        momfluxb(i,j,kb) = momfluxb(i,j,kb) + bcmom * dzfi(kb)
        emom = (dzf(kb-1)*(ekm(i,j,kb)*dxf(i-1) + ekm(i-1,j,kb)*dxf(i)) + &
                dzf(kb)*(ekm(i,j,kb-1)*dxf(i-1) + ekm(i-1,j,kb-1)*dxf(i)))*dxhi(i)*dzhiq(kb)
        up(i,j,kb) = up(i,j,kb) + (u0(i,j,kb) - u0(i,j,kb-1))*emom*dzhi(kb)*dzfi(kb) - bcmom*dzfi(kb)
      end do
    end do
    do j = jb, je
      do i = ib, ie
        utan2 = 0.25 * (u0(i,j,kb) + u0(i,j-1,kb) + u0(i+1,j-1,kb) + u0(i+1,j,kb))
        horv  = sqrt(max(min_horv**2, utan2**2 + v0(i,j,kb)**2))
        bcmom = ustar_mo**2 * v0(i,j,kb) / horv
        momfluxb(i,j,kb) = momfluxb(i,j,kb) + bcmom * dzfi(kb)
        emom = (dzf(kb-1)*(ekm(i,j,kb) + ekm(i,j-1,kb)) + dzf(kb)*(ekm(i,j,kb-1) + ekm(i,j-1,kb-1)))*dzhiq(kb)
        vp(i,j,kb) = vp(i,j,kb) + (v0(i,j,kb) - v0(i,j,kb-1))*emom*dzhi(kb)*dzfi(kb) - bcmom*dzfi(kb)
      end do
    end do
  end subroutine surface_mo

  !> Port of DALES dev 075ac4e modsurface.f90:calc_obl_iter (Newton iteration on Rib).
  real function calc_obl_iter(thl, qt, tskin, qskin, z, z0m, z0h, u, v, Lin) result(L)
    use modglobal, only : grav, rv, rd
    real, intent(in) :: thl, qt, tskin, qskin, z, z0m, z0h, u, v, Lin
    real :: horv2, rib, thv, thvsl, fx, fxdif, lold, lstart, lend
    integer :: iter
    L = Lin
    thv   = thl   * (1. + (rv/rd - 1.) * qt)
    thvsl = tskin * (1. + (rv/rd - 1.) * qskin)
    horv2 = max(u**2 + v**2, min_horv**2)
    rib   = grav / thvsl * z * (thv - thvsl) / horv2
    if (rib == 0.) then
      L = 1.e6
      return
    end if
    iter = 0
    lold = 0.
    if (rib * L < 0. .or. abs(L) == 1.e5) then
      if (rib > 0.) L =  0.01
      if (rib < 0.) L = -0.01
    end if
    do while (abs((L - lold) / L) > 1.e-4 .and. iter < 1000)
      iter   = iter + 1
      lold   = L
      lstart = L - 0.001 * L
      lend   = L + 0.001 * L
      fx     = rib - rib_from_obl(z, L, z0h, z0m)
      fxdif  = (- rib_from_obl(z, lstart, z0h, z0m) + rib_from_obl(z, lend, z0h, z0m)) / (lstart - lend)
      L      = L - fx / fxdif
      if (rib * L < 0. .or. abs(L) == 1.e5) then
        if (rib > 0.) L =  0.01
        if (rib < 0.) L = -0.01
      end if
      if (abs(L) > 1.e6) L = sign(1.e6, L)
    end do
  end function calc_obl_iter

  elemental real function rib_from_obl(z, L, z0h, z0m)
    real, intent(in) :: z, L, z0h, z0m
    rib_from_obl = z / L * (log(z/z0h) - psih(z/L) + psih(z0h/L)) &
                         / (log(z/z0m) - psim(z/L) + psim(z0m/L))**2
  end function rib_from_obl

  ! DALES dev 075ac4e modsurface.f90: psim, psih (Businger-Dyer / Beljaars-Holtslag), phim, phih
  elemental real function psim(zeta)
    real, intent(in) :: zeta
    real :: x
    if (zeta <= 0.) then
      x    = (1. - 16. * zeta) ** 0.25
      psim = 3.14159265 / 2. - 2. * atan(x) + log((1. + x)**2 * (1. + x**2) / 8.)
    else
      psim = -2./3. * (zeta - 5./0.35) * exp(-0.35 * zeta) - zeta - (10./3.) / 0.35
    end if
  end function psim

  elemental real function psih(zeta)
    real, intent(in) :: zeta
    real :: x
    if (zeta <= 0.) then
      x    = (1. - 16. * zeta) ** 0.25
      psih = 2. * log((1. + x**2) / 2.)
    else
      psih = -2./3. * (zeta - 5./0.35) * exp(-0.35 * zeta) - (1. + (2./3.) * zeta) ** 1.5 &
             - (10./3.) / 0.35 + 1.
    end if
  end function psih

  elemental real function phih(zeta)
    real, intent(in) :: zeta
    if (zeta < 0.) then
      phih = (1. - 16. * zeta) ** (-0.50)
    else if (zeta < 1.) then
      phih = 1. + 5. * zeta
    else
      phih = 6.
    end if
  end function phih

end module modsurface
```

Notes:
- **`dudz`/`dvdz`/`phim`.** DALES also computes these for its SGS surface treatment. uDALES's SGS does not consume them, so they are not ported (YAGNI).
- **`momfluxb` bounds.** Check them in `modfields` against `wfmneutral`'s `iomomflux(ib-hi:ie+hi, jb-hj:je+hj, kb-hk:ke+hk)`, and that `bottom` passes `momfluxb` the same way.
- **`u0av`, `v0av`, `thl0av`, `qt0av`.** These are the slab means from the last `diagfld`. Confirm they are filled before the first step: `readinitfiles` → `thermodynamics` at `timee == 0`. If they are not, call `surface_mo` only from the second step on, and set `ustar_mo` from neutral `cm` on the first.

Wiring:
- **`modglobal`:** add `integer, parameter :: BCbotm_mo = 4 ! Monin-Obukhov (DALES isurf 4)` and `BCbotm_mo_ustar = 5 ! Monin-Obukhov with prescribed ustar (DALES isurf 3)` after line 165.
- **`modsurfdata`:** `real :: ustin = -1. !< prescribed friction velocity for BCbotm = 5 [m/s]`.
- **`modstartup`:** add `ustin` to `namelist/BC/` and `call MPI_BCAST(ustin, 1, MY_REAL, 0, comm3d, mpierr)` next to the `z0` broadcast.
- **`checkinitvalues`:**

```fortran
      if ((BCbotm == BCbotm_mo) .or. (BCbotm == BCbotm_mo_ustar)) then
         if (.not. lbottom) then
            write(0, *) 'ERROR: BCbotm = 4/5 (Monin-Obukhov) requires lbottom = .true.'
            stop 1
         end if
         if ((.not. ltempeq) .or. BCbotT /= BCbotT_flux .or. (z0 <= 0.) .or. (z0h <= 0.)) then
            write(0, *) 'ERROR: BCbotm = 4/5 requires ltempeq, BCbotT = 1 and z0, z0h > 0'
            stop 1
         end if
         if (wtsurf == -1. .or. (lmoist .and. wqsurf == -1.)) then
            write(0, *) 'ERROR: BCbotm = 4/5: wtsurf (and wqsurf if lmoist) still at sentinel -1'
            stop 1
         end if
         if ((BCbotm == BCbotm_mo_ustar) .and. (ustin <= 0.)) then
            write(0, *) 'ERROR: BCbotm = 5 requires ustin > 0'
            stop 1
         end if
      end if
```
  (`lbottom` lives in `modibm`; check that `checkinitvalues` can `use modibm, only : lbottom` without a circular dependency. If it can't, put the check at the start of `bottom` on the first call.)
- **`bottom`** (`modibm.f90`): add `use modsurface, only : surface_mo` and in the momentum `if` chain:

```fortran
      elseif ((BCbotm.eq.BCbotm_mo) .or. (BCbotm.eq.BCbotm_mo_ustar)) then
      call surface_mo
```
  (`BCbotm_mo`, `BCbotm_mo_ustar` added to its `use modglobal` list).
- **`modthermostat`:** `use modsurface, only : ustar_mo, obl_mo` and set `vt(8) = ustar_mo; vt(9) = obl_mo`.
- **Docs:** add the two options to the `BCbotm` table in `docs/udales-boundary-conditions.md`, and `ustin` to `docs/udales-namoptions-overview.md`.
- **`namelists.json`:** `ustin` under `BC`.

- [ ] **Step 6: Build, test, regress**

```bash
./tools/build_executable.sh icl debug 2>&1 | tail -3
python tests/integration/thermo/test_thermo.py -v
python tests/run_tests.py supported --branch-a master --branch-b HEAD --build-type Debug 2>&1 | tail -30
```
Expected: all thermo tests PASS, and supported matches the Task 1 result.

- [ ] **Step 7: Commit**

```bash
git add src/modsurface.f90 src/modibm.f90 src/modglobal.f90 src/modsurfdata.f90 src/modstartup.f90 src/modthermostat.f90 docs/udales-boundary-conditions.md docs/udales-namoptions-overview.md tools/python/namelists.json tests/integration/thermo/thermo_case.py tests/integration/thermo/test_thermo.py
git status
git commit -m "Add Monin-Obukhov surface layer (BCbotm = 4/5, DALES isurf 4/3)

Port of DALES dev 075ac4e modsurface (lmostlocal = .false.).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Deardorff `calthv` surface level and `tmp0`

**Files:**
- Modify: `src/modthermodynamics.f90:130-235` (`calthv`)
- Modify: `src/modfielddump.f90` (code `dt` for `dthvdz`, both branches, ~134 and ~201; add `dthvdz` to the `use modfields` list at ~60)
- Test: `test_thermo.py::TestCalthv`

**Interfaces:**
- Consumes: `modsurface: lmo_surface, dthldz_mo, dqtdz_mo`.
- Produces: fielddump code `'dt'` → variable `dthvdz`. (`tmp0` from the spec is dropped: its only consumer is the deferred microphysics.)

- [ ] **Step 1: Write the failing test**

```python
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
        d = tc.write_case(self.tmp / "calthv", u=7.0, v=-5.5, groups={
            "RUN": {"runtime": "2.5"},
            "OUTPUT": {"lfielddump": ".true.", "tfielddump": "1.", "fieldvars": "'th,qt,ql,dt'"},
            "BC": {"BCbotm": "4", "wtsurf": "-0.05", "wqsurf": "-5.0e-5", "z0": "2.e-4", "z0h": "2.e-4"},
            "NAMTHERMOSTAT": {"lthermostat": ".true.", "tthermostat": "1."}})
        self.run_ok(d)
        # Record alignment: fielddump and thermostat run after tstep_integrate and BEFORE
        # thermodynamics. So dump record n holds end-of-step-n thl/qt, and the dthvdz that
        # calthv computes from them (with ustar/obl of step n, thermostat record n)
        # appears in dump record n+1. ql/exnf in record n are from the start of step n,
        # so take ql from record n+1 as well (it is diagnosed from record-n thl/qt).
        n = 0
        with Dataset(d / "fielddump.000.000.310.nc") as f, \
             Dataset(d / "thermostat_t.310.nc") as t, Dataset(d / "thermostat_z.310.nc") as z:
            dthv = np.asarray(f["dthvdz"][n + 1])        # (z, y, x)
            thl1 = np.asarray(f["thl"][n])[0].mean()
            qt1 = np.asarray(f["qt"][n])[0].mean()
            ql1 = np.asarray(f["ql"][n + 1])[0].mean()
            exn1 = float(z["exnf"][n + 1, 0])
            ustar, obl = float(t["ustar"][n]), float(t["obl"][n])
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
```

Check the fielddump variable names (`thl` for code `th`), the axis order against `modfielddump.f90` (`tttt` = x,y,z,time in Fortran order, i.e. `(time, z, y, x)` in netCDF4-python), the file name pattern (`fielddump.<myidx>.<myidy>.<expnr>.nc`, 3-digit fields), and that `tfielddump = 1` gives one record per step from t = 1 s. If the first record is at t = 0 instead, shift `n` by one and keep the alignment rule.

- [ ] **Step 2: Run to verify failure**

Expected: FAIL. Code `dt` is unknown and silently dumps `u0`, so the `dthvdz` variable is missing (KeyError).

- [ ] **Step 3: Implement**

`calthv` (`modthermodynamics.f90`):
- add `use modsurface, only : lmo_surface, dthldz_mo, dqtdz_mo`;
- declare `real :: a_surf, b_surf`.

Replace the moist-branch `dthvdz(i,j,kb)=0.` loop with:

```fortran
       ! surface level: DALES dev 075ac4e calthv. Needs MO surface gradients;
       ! with other bottom BCs dthvdz(kb) keeps its previous (zero) value.
       do j=jb,je
          do i=ib,ie
             if (ql0(i,j,kb) > 0) then
                temp = thl0(i,j,kb)*exnf(kb)+(rlv/cp)*ql0(i,j,kb)
                qs   = qt0(i,j,kb) - ql0(i,j,kb)
                a_surf = (1.-qt0(i,j,kb)+rv/rd*qs*(1.+rlv/(rv*temp))) &
                         /(1.+rlv**2*qs/(cp*rv*temp**2))
                b_surf = a_surf*rlv/(temp*cp)-1.
             else
                temp = thl0(i,j,kb)*exnf(kb)
                a_surf = 1.+(rv/rd-1.)*qt0(i,j,kb)
                b_surf = rv/rd-1.
             end if
             if (lmo_surface()) then
                dthvdz(i,j,kb) = a_surf*dthldz_mo + b_surf*thl0(i,j,kb)*dqtdz_mo
             else
                dthvdz(i,j,kb) = 0.
             end if
          end do
       end do
```

In the dry branch, replace `dthvdz(i,j,kb)=0.` with `dthvdz(i,j,kb) = merge(dthldz_mo, 0., lmo_surface())`.

The existing `eps1` floor loop afterwards stays, so the non-MO surface value becomes `±eps1` as before. Bit-for-bit is preserved.

`modfielddump.f90`, both `select case` branches:

```fortran
          case('dt')
            call ncinfo(ncname( n,:),'dthvdz','Vertical gradient of thv','K/m','tttt')
            pfields(n)%point => dthvdz(ib:ie,jb:je,kb:ke)
```
In the halo branch (`lhalos`, currently hard-coded off) the shapes do not match (`dthvdz` starts at `kb`, `vars` at `kb-kh`). Add the case there as `write(0, *) 'ERROR: fieldvars code dt is not available with halos'; stop 1`. Add a `dt` row to the `fieldvars` code list in `docs/udales-namoptions-overview.md`.

- [ ] **Step 4: Build, test, regress**

Same commands as Task 3, Step 6. Expected: all PASS, and supported unchanged.

- [ ] **Step 5: Commit**

```bash
git add src/modthermodynamics.f90 src/modfielddump.f90 docs/udales-namoptions-overview.md tests/integration/thermo/test_thermo.py
git status
git commit -m "Deardorff calthv surface level from MO gradients

Port of DALES dev 075ac4e modthermodynamics:calthv. Unchanged unless
BCbotm = 4/5. Adds fielddump code 'dt' (dthvdz).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Parametrised longwave radiation

**Files:**
- Create: `src/modraddata.f90`, `src/modradiation.f90`
- Modify: `src/program.f90` (`initradiation` after `initthermodynamics`; `radiation` after `lstend`; `exitradiation` before `exitmpi`)
- Modify: `src/modthermostat.f90` (`thlprad`, `lwu`, `lwd`)
- Modify: `tools/python/namelists.json`, `docs/udales-namoptions-overview.md`
- Test: `test_thermo.py::TestRadiation`

**Interfaces:**
- Produces: module `modraddata` (public data) `iradiation, irad_none=0, irad_par=2, timerad, rad_longw, rad_shortw, rad_smoke, rad_ls, lrad_dycoms, isvsmoke, rka, dlwtop, dlwbot, rad_dycoms_D, rad_dycoms_qtthres, thlprad(ib:ie,jb:je,kb:ke), lwd(ib:ie,jb:je,kb:ke+1), lwu(ib:ie,jb:je,kb:ke+1)`. Module `modradiation`: `initradiation`, `radiation`, `exitradiation`.

- [ ] **Step 1: Write the failing tests**

```python
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

    def test_liquid_longwave_matches_formula(self):
        d = self._run("rad_liq", self.RAD)
        with Dataset(d / "thermostat_z.310.nc") as z:
            ql, rho, exn = (np.asarray(z[v][-1]) for v in ("ql", "rhof", "exnf"))
            got = np.asarray(z["thlprad"][-1])
        ref = radpar_reference(ql, rho, exn, 20.0, 70.0, 22.0, 85.0)
        # thlprad uses start-of-step ql; the record's ql is end-of-step (~0.3 % change at
        # cloud top). A sign or face off-by-one error is O(1), so 1e-2 still discriminates.
        np.testing.assert_allclose(got, ref, rtol=1e-2, atol=1e-2 * np.abs(ref).max())
        self.assertLess(got.min(), -1e-4)          # cloud-top cooling present

    def test_smoke_longwave_matches_formula(self):
        rad = dict(self.RAD, rad_smoke=".true.", isvsmoke="1", rka="0.02", dlwtop="60.", dlwbot="0.")
        zf = (np.arange(64) + 0.5) * 20.0
        d = self._run("rad_smoke", rad, groups={"PHYSICS": {"lmoist": ".false."},
                                                "SCALARS": {"nsv": "1"}},
                      thl=np.where(zf < 700, 288.0, 295.0), qt=np.zeros(64))
        # scalar.inp: smoke = 1 below 700 m, 0 above (see note)
        with Dataset(d / "thermostat_z.310.nc") as z:
            rho, exn = np.asarray(z["rhof"][-1]), np.asarray(z["exnf"][-1])
            got = np.asarray(z["thlprad"][-1])
        smoke = np.where(zf < 700, 1.0, 0.0)
        ref = radpar_reference(smoke, rho, exn, 20.0, 60.0, 0.0, 0.02)
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
        ft = zf > 840.0 + 150.0
        # Residual sources: one-sided subsidence differences, O(dz/(z-zi)), and the
        # 832.7 vs 840 m zi offset. A sign error gives a ratio of ~2.
        self.assertLess(np.abs(dthl_both[ft]).max(), 0.25 * np.abs(dthl_subs[ft]).max())
```

Notes for the implementer:
- **Smoke test scalar input.** It needs a `scalar.inp.310`: 2 header lines, then `z sv1` per level, read when `nsv > 0` on a cold start (`modstartup.f90:~1545`). Add an optional `sv` argument to `write_case` that writes it, and check the column layout in the reader first. Also check which namelist group holds `nsv` (`&SCALARS`?) and fix the dict key if it differs.
- **Tolerances** (from the Fable review): records are float32, so every comparison uses ≥ 1e-2 relative, or a long enough run. Do not tighten them, and do not loosen them beyond the values shown without reporting back.

- [ ] **Step 2: Run to verify failure**

Expected: FAIL. `NAMRADIATION` is ignored and `thlprad` is all zeros.

- [ ] **Step 3: Implement `src/modraddata.f90`**

```fortran
!> Radiation switches, parameters and fields (subset of DALES dev 075ac4e modraddata).
module modraddata
  implicit none
  save
  integer, parameter :: irad_none = 0     !< no radiation
  integer, parameter :: irad_par  = 2     !< parametrised longwave (DALES radpar)
  integer :: iradiation = irad_none
  real    :: timerad    = 0.              !< update interval [s]; 0 = every step
  logical :: rad_longw  = .true.          !< longwave on
  logical :: rad_shortw = .false.         !< not supported: aborts
  logical :: rad_smoke  = .false.         !< absorber is scalar isvsmoke instead of ql0
  logical :: rad_ls     = .false.         !< DALES compatibility only (thlpcar is applied in forces)
  logical :: lrad_dycoms = .false.        !< DYCOMS free-troposphere term (Stevens et al. 2005)
  integer :: isvsmoke   = 1               !< scalar index of the smoke absorber
  real    :: rka        = 130.            !< extinction coefficient kappa [m2/kg]
  real    :: dlwtop     = 74.             !< F0 [W/m2]
  real    :: dlwbot     = 0.              !< F1 [W/m2]
  real    :: rad_dycoms_D       = 3.75e-6 !< large-scale divergence in the FT term [1/s]
  real    :: rad_dycoms_qtthres = 8.e-3   !< qt defining the local inversion height [kg/kg]
  real    :: tnextrad = 0.
  real, allocatable :: thlprad(:,:,:)     !< radiative thl tendency [K/s] (ib:ie,jb:je,kb:ke)
  real, allocatable :: lwd(:,:,:)         !< F0 term of net upward LW flux at faces [W/m2] (kb:ke+1)
  real, allocatable :: lwu(:,:,:)         !< F1 term of net upward LW flux at faces [W/m2] (kb:ke+1)
end module modraddata
```

- [ ] **Step 4: Implement `src/modradiation.f90`**

```fortran
!> Parametrised longwave radiation, ported from DALES dev 075ac4e
!> src/modradiation.f90 (radiation, radpar longwave branch) and the
!> DYCOMS rad_user of cases/dycoms_rf02/moduser.f90, which here is the
!> built-in option lrad_dycoms.
module modradiation
  use modraddata
  implicit none
  private
  public :: initradiation, radiation, exitradiation

contains

  subroutine initradiation
    use modglobal, only : ifnamopt, fname_options, ib, ie, jb, je, kb, ke, nsv, lmoist, ltempeq
    use modmpi,    only : myid, comm3d, mpierr, my_real, mpi_logical, mpi_integer
    implicit none
    integer :: ierr
    namelist/NAMRADIATION/ iradiation, timerad, rad_longw, rad_shortw, rad_smoke, rad_ls, &
         lrad_dycoms, isvsmoke, rka, dlwtop, dlwbot, rad_dycoms_D, rad_dycoms_qtthres

    if (myid == 0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read (ifnamopt, NAMRADIATION, iostat=ierr)
      if (ierr > 0) then
        write(0, *) 'ERROR: Problem in namoptions NAMRADIATION'
        write(0, *) 'iostat error: ', ierr
        stop 1
      end if
      close(ifnamopt)
    end if
    call MPI_BCAST(iradiation,  1, MPI_INTEGER, 0, comm3d, mpierr)
    call MPI_BCAST(timerad,     1, MY_REAL,     0, comm3d, mpierr)
    call MPI_BCAST(rad_longw,   1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(rad_shortw,  1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(rad_smoke,   1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(rad_ls,      1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(lrad_dycoms, 1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(isvsmoke,    1, MPI_INTEGER, 0, comm3d, mpierr)
    call MPI_BCAST(rka,         1, MY_REAL,     0, comm3d, mpierr)
    call MPI_BCAST(dlwtop,      1, MY_REAL,     0, comm3d, mpierr)
    call MPI_BCAST(dlwbot,      1, MY_REAL,     0, comm3d, mpierr)
    call MPI_BCAST(rad_dycoms_D,       1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(rad_dycoms_qtthres, 1, MY_REAL, 0, comm3d, mpierr)

    allocate(thlprad(ib:ie,jb:je,kb:ke)); thlprad = 0.
    allocate(lwd(ib:ie,jb:je,kb:ke+1));   lwd = 0.
    allocate(lwu(ib:ie,jb:je,kb:ke+1));   lwu = 0.
    if (iradiation == irad_none) return

    if (iradiation /= irad_par) then
      write(0, *) 'ERROR: iradiation must be 0 (none) or 2 (parametrised); got ', iradiation
      stop 1
    end if
    if (rad_shortw) then
      write(0, *) 'ERROR: rad_shortw is not supported in uDALES'
      stop 1
    end if
    if (.not. ltempeq) then
      write(0, *) 'ERROR: iradiation = 2 requires ltempeq = .true.'
      stop 1
    end if
    if (rad_smoke .and. (isvsmoke < 1 .or. isvsmoke > nsv)) then
      write(0, *) 'ERROR: rad_smoke requires 1 <= isvsmoke <= nsv'
      stop 1
    end if
    if ((.not. rad_smoke) .and. (.not. lmoist)) then
      write(0, *) 'ERROR: liquid-water radiation (rad_smoke = .false.) requires lmoist = .true.'
      stop 1
    end if
    if (lrad_dycoms .and. (.not. lmoist)) then
      write(0, *) 'ERROR: lrad_dycoms requires lmoist = .true.'
      stop 1
    end if
    if (rad_ls .and. myid == 0) then
      write(0, *) 'WARNING: rad_ls ignored; the thlpcar column of lscale.inp is always applied in forces'
    end if
  end subroutine initradiation

  !> Recompute thlprad at the first RK substep when due; add it to thlp every substep.
  subroutine radiation
    use modglobal, only : ib, ie, jb, je, kb, ke, timee, rk3step
    use modfields, only : thlp
    implicit none
    if (iradiation == irad_none) return
    if (rk3step == 1 .and. timee >= tnextrad) then
      tnextrad = tnextrad + timerad
      thlprad = 0.
      if (rad_longw) call radpar
    end if
    thlp(ib:ie,jb:je,kb:ke) = thlp(ib:ie,jb:je,kb:ke) + thlprad
  end subroutine radiation

  !> DALES dev 075ac4e modradiation.f90:radpar, longwave branch, plus the
  !> DYCOMS free-troposphere term (moduser rad_user of dycoms_rf02).
  subroutine radpar
    use modglobal, only : ib, ie, jb, je, kb, ke, cp, dzf, zf
    use modfields, only : ql0, qt0, sv0, rhof, exnf, IIc
    implicit none
    real :: lwpt(kb:ke+1), lwpb(kb:ke+1), absorber(kb:ke), zi
    integer :: i, j, k

    do j = jb, je
      do i = ib, ie
        if (rad_smoke) then
          absorber = sv0(i,j,kb:ke,isvsmoke) * IIc(i,j,kb:ke)
        else
          absorber = ql0(i,j,kb:ke) * IIc(i,j,kb:ke)
        end if
        lwpt(ke+1) = 0.
        do k = ke, kb, -1
          lwpt(k) = lwpt(k+1) + rhof(k) * absorber(k) * dzf(k)
        end do
        lwpb(kb) = 0.
        ! DALES writes rhof(k)*absorber(k-1)*dzf(k-1); rhof(k-1) is the consistent
        ! density of the cell whose path is added (deliberate, tiny difference).
        do k = kb + 1, ke + 1
          lwpb(k) = lwpb(k-1) + rhof(k-1) * absorber(k-1) * dzf(k-1)
        end do
        do k = kb, ke + 1
          lwd(i,j,k) = dlwtop * exp(-rka * lwpt(k))
          lwu(i,j,k) = dlwbot * exp(-rka * lwpb(k))
        end do
        do k = kb, ke
          thlprad(i,j,k) = - ((lwd(i,j,k+1) + lwu(i,j,k+1)) - (lwd(i,j,k) + lwu(i,j,k))) &
                           / (rhof(k) * cp * exnf(k) * dzf(k))
        end do

        if (lrad_dycoms) then
          ! Local inversion height: highest crossing of the qt threshold. Two fixes
          ! relative to DALES rad_user: zi is reset per column (DALES reuses the
          ! previous column's value) and the i loop covers the full range (DALES ran
          ! i = 2, j1). Sign: COOLING, -D/3[...]; DALES adds +D/3[...], which doubles
          ! subsidence warming instead of balancing it. With thl = 297.5 + (z-zi)^(1/3),
          ! -D/3[(z-zi)^(1/3) + zi (z-zi)^(-2/3)] = -(D/3) z (z-zi)^(-2/3)
          ! = w dthl/dz for w = -D z, i.e. exactly minus the subsidence warming.
          ! rho_i/rho ~ 1 and no 1/exner, as in DALES.
          zi = -1.
          do k = ke - 1, kb, -1
            if (qt0(i,j,k) > rad_dycoms_qtthres) then
              zi = zf(k) + (rad_dycoms_qtthres - qt0(i,j,k)) * (zf(k+1) - zf(k)) / (qt0(i,j,k+1) - qt0(i,j,k))
              exit
            end if
          end do
          if (zi > 0.) then
            do k = ke, kb, -1
              if (zf(k) <= zi) exit
              thlprad(i,j,k) = thlprad(i,j,k) - rad_dycoms_D / 3. * &
                               ((zf(k) - zi)**(1./3.) + zi * (zf(k) - zi)**(-2./3.))
            end do
          end if
        end if

        thlprad(i,j,kb:ke) = thlprad(i,j,kb:ke) * IIc(i,j,kb:ke)
      end do
    end do
  end subroutine radpar

  subroutine exitradiation
    implicit none
    if (allocated(thlprad)) deallocate(thlprad, lwd, lwu)
  end subroutine exitradiation

end module modradiation
```

Check these against the real code:
- **Initial state.** `rhof`, `exnf`, `ql0`, `qt0` must be valid at the first `radiation` call: the initial `thermodynamics` call happens during `readinitfiles`. Confirm.
- **Substep reset.** `thlp` must be reset each substep before the tendencies are added. Find where uDALES zeroes `thlp` (likely in `tstep_integrate`), so that adding `thlprad` once per substep is right, exactly as DALES does.
- **First-call timing.** `timee` at the first step is 0, so `tnextrad = 0` triggers on the first call.

`program.f90`:
- `use modradiation, only : initradiation, radiation, exitradiation`;
- `call initradiation` after `call initthermodynamics`;
- `call radiation` after `call lstend`;
- `call exitradiation` after `call exitthermostat`.

`modthermostat`: `use modraddata, only : thlprad, lwu, lwd, iradiation, irad_none, isvsmoke`; change the smoke `crossing` call to use `sv0(i,j,kb:ke,isvsmoke)` when `rad_smoke`, else keep scalar 1. When `iradiation /= irad_none`, average `thlprad` into `vz(:,7)` with `avexy_ibm` (IIc), and `lwu`/`lwd` at faces `kb:ke` into `vz(:,10)`/`vz(:,11)` (use a mask of ones and `IIws`; the top face `ke+1` is not written).

`namelists.json`: add a `NAMRADIATION` group with the 13 names; document the group in `docs/udales-namoptions-overview.md` with the parameter table from the spec.

- [ ] **Step 5: Build, test, regress**

Same commands as Task 3, Step 6. Expected: all PASS, and supported unchanged.

- [ ] **Step 6: Commit**

```bash
git add src/modraddata.f90 src/modradiation.f90 src/program.f90 src/modthermostat.f90 tools/python/namelists.json docs/udales-namoptions-overview.md tests/integration/thermo/thermo_case.py tests/integration/thermo/test_thermo.py
git status
git commit -m "Add parametrised longwave radiation (iradiation = 2, smoke, DYCOMS)

Port of DALES dev 075ac4e radpar (longwave). DYCOMS free-troposphere term
built in as lrad_dycoms, with sign and per-column zi fixed relative to
DALES rad_user (cooling balances subsidence warming; see comment).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Decomposition independence and CI registration

**Files:**
- Modify: `tests/integration/thermo/test_thermo.py`
- Create: `tests/integration/thermo/README.md`
- Modify: `tests/test_suites.yml` (`supported` group, after the `SGS TKE top boundary condition` entry ~150–161)

- [ ] **Step 1: Perturbations**

The review established that `randomize_field` uses global indices (`modstartup.f90:2504-2516`), so the perturbations are decomposition-independent. This test turns them on (`lrandomize = .true.`) to get a turbulent 10-step run.

- [ ] **Step 2: Write the test**

```python
class TestDecomposition(_SolverTest):
    """Radiation + MO surface + moist SGS give the same answer on 1, 2 and 4 ranks."""

    def test_splits_agree(self):
        groups = {"RUN": {"runtime": "9.5", "lrandomize": ".true."},
                  "BC": {"BCbotm": "4", "wtsurf": "-0.012", "wqsurf": "-4.7e-5",
                         "z0": "2.e-4", "z0h": "2.e-4"},
                  "NAMRADIATION": dict(TestRadiation.RAD, lrad_dycoms=".true."),
                  "NAMTHERMOSTAT": {"lthermostat": ".true.", "tthermostat": "3."}}
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
                # float32 records: 1e-6 is bitwise up to reduction-order round-off
                np.testing.assert_allclose(out[label][v], ref, rtol=1e-6, atol=1e-12, err_msg=f"{label}:{v}")
```

uDALES requires `jtot` and `ktot` to be divisible by `nprocy` (16 and 64: fine). Check `mpi_env.preflight(4)` on the login node. If 4 ranks are refused there, note it; the suite still runs in CI and on compute nodes.

- [ ] **Step 3: README** (`tests/integration/thermo/README.md`). One short section per test class: what it proves, and the independent reference it compares against. Also list the run requirements (`UDALES_BUILD`, `UDALES_RUNTIME_MODULES` with the build stack, `netCDF4`, up to 4 MPI ranks).

- [ ] **Step 4: Register in `tests/test_suites.yml`**, matching the existing entry format exactly (the parser is hand-written):

```yaml
    - label: "supported: moist thermodynamics, MO surface and radiation"
      class: supported
      kind: integration
      component: solver
      platform: hpc
      cost: fast
      env_UDALES_BUILD: "{repo_root}/build/{build_type_lower}/u-dales"
      command:
        - "{python}"
        - "{tests_dir}/integration/thermo/test_thermo.py"
```

Check `.github/workflows/ci.yml` for how CI provides `netCDF4` and how many MPI ranks it allows. If CI caps ranks at 2 or 3, skip the `xy` split in CI only: `if mpi_env` reports fewer than 4 slots, use `self.skipTest`-style handling for that sub-case. Look at how `supported-macos` documents its rank limit.

- [ ] **Step 5: Run the full selection**

```bash
python tests/run_tests.py supported --branch-a master --branch-b HEAD --build-type Debug 2>&1 | tail -40
```
Expected: the new suite PASSES, and every other suite matches the Task 1 result.

- [ ] **Step 6: Commit**

```bash
git add tests/integration/thermo/test_thermo.py tests/integration/thermo/README.md tests/test_suites.yml
git status
git commit -m "Register thermo integration tests in supported; add decomposition test

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Full-size cases `examples/301`–`303` and tooling

**Files:**
- Create: `tools/python/examples/dales_cases/convert_dales_case.py`, `make_rf01.py`, `compare_dales.py`, `README.md`
- Create: `examples/301/` (BOMEX), `examples/302/` (smoke), `examples/303/` (RF01): `namoptions.30x`, `prof.inp.30x`, `lscale.inp.30x`, (`scalar.inp.302`), `README.md`, `config.sh`

- [ ] **Step 1: `convert_dales_case.py`**

CLI: `python convert_dales_case.py <dales_case_dir> <expnr> <out_dir>`. It does the following:
- reads DALES `prof.inp.001` (`z thl qt u v tke`) and writes uDALES `prof.inp.<expnr>` (`z thl qt u v e12`), same columns;
- reads DALES `lscale.inp.001` (`z ug vg wfls dqtdxls dqtdyls dqtdtls dthlrad`) and writes the uDALES 10-column layout `z ug vg pgx=0 pgy=0 wfls dqtdxls dqtdyls dqtdtls thlpcar=dthlrad`;
- copies `scalar.inp.001` (if any) to `scalar.inp.<expnr>`, after checking the column layout against the uDALES reader;
- prints the DALES namoptions values the case needs, so they can be transcribed (`ps, thls, wtsurf, wqsurf, z0, ustin, isurf, itot, jtot, xsize, ysize, runtime, lcoriol`, Coriolis latitude, radiation block), and maps `isurf 3 → BCbotm 5`, `isurf 4 → BCbotm 4`;
- **prints `wtsurf`/`wqsurf` negated** (DALES fluxes are upward-positive; uDALES prescribed fluxes are positive in −z, decision D1).

Header lines are skipped by detecting lines that do not parse as floats, since DALES files have one or two header lines. Unit test: run it on `$SCRATCH/dales/cases/bomex` and check that row 1 of the output matches row 1 of the input column for column.

- [ ] **Step 2: `make_rf01.py`**

CLI: `python make_rf01.py <out_dir> [--dz 5 --ktot 300 --itot 96 --dx 35]`. Profiles:
- θl = 289 K and qt = 9 g/kg for z ≤ 840 m; θl = 297.5 + (z−840)^(1/3) and qt = 1.5 g/kg above;
- u = 7, v = −5.5, e12 = 1 below zi, 0 above (check the uDALES `e12` meaning: it is √e, not e; set accordingly);
- `lscale`: ug = 7, vg = −5.5, wfls = −3.75e−6·z.

Cite PyCLES `InitDYCOMS_RF01`/`ForcingDyCOMS_RF01`/`SurfaceDYCOMS_RF01`/`RadiationDyCOMS_RF01` and Stevens et al. (2005) in the module docstring. Surface:
- sensible heat 15 W m⁻² (upward) → `wtsurf = −15/(ρ cp)` with ρ = 1.22 kg m⁻³ (PyCLES); negative because uDALES fluxes are positive in −z;
- latent heat 115 W m⁻² (upward) → `wqsurf = −115/(ρ Lv)`, using uDALES `rlv = 2.26e6` (note that PyCLES uses 2.47e6 in its temperature formula; document the choice);
- `BCbotm = 4` with `z0` chosen so that neutral `Cm = κ²/ln²(zf1/z0) = 0.0011` at `zf1 = dz/2` (for dz = 5 m: `z0 = 2.5·exp(−0.41/√0.0011)`); print the value and use it for `z0` and `z0h`.

Radiation: `iradiation = 2, rka = 85, dlwtop = 70, dlwbot = 22, lrad_dycoms = .true.`. Coriolis on at `xlat = 31.5`. Domain 3360 × 3360 m, 4 h.

- [ ] **Step 3: Write the three case directories.**

Common settings:
- `libm = .false.`, `lbottom = .true.`, `nfcts = 0`;
- `ltempeq = lbuoyancy = .true.`;
- `loneeqn = .true.`;
- `ladaptive = .true.`, `courant = 0.5`;
- sponge as in the DALES case (check `igrw_damp`/`geodamptime` in uDALES `&PHYSICS`/`&DYNAMICS` and match DALES's damping depth);
- `lthermostat = .true.`, `tthermostat = 60.`;
- statsdump off;
- `config.sh` for `tools/hpc_execute.sh` with `NCPU` chosen so `jtot % nprocy == 0` and `ktot % nprocy == 0`.

Per case:
- 301 BOMEX: from `cases/bomex`, `BCbotm = 5`, `ustin = 0.28`, `z0 = 2e-4`, `wtsurf = -8e-3`, `wqsurf = -5.2e-5` (DALES +8e-3/+5.2e-5 upward), `thls = 299.1`, `lmoist`, `lcoriol` (`xlat` such that f = 0.376e-4 s⁻¹, i.e. 14.94°), 6 h.
- 302 smoke: from `cases/smoke`, `BCbotm = 4`, `z0 = 0.1`, zero fluxes, `lmoist = .false.`, `nsv = 1`, `iadv_sv = 7` (kappa, as DALES `iadv_sv = 7`; `&DYNAMICS`), `iradiation = 2`, `rad_smoke`, `rka = 0.02`, `dlwtop = 60`, 3 h.
- 303 RF01: from `make_rf01.py`.

All three cases use `timerad = 0` (radiation every step): `radiation` has no DALES-style `dt_lim` clamp, so a non-zero `timerad` with `ladaptive` would update late. Each case README lists every number with its source (DALES namoptions line or paper/PyCLES), plus the known deviations: uDALES has no 5th-order advection, its SGS details differ, and the RF01 drag uses an equivalent `z0`.

- [ ] **Step 4: Smoke-test each case for 60 s on a login node** with 4 ranks and `runtime = 60` (override via a temporary copy), to check that the inputs parse and nothing crashes. Expected: `thermostat_t.30x.nc` has records and no NaN.

- [ ] **Step 5: `compare_dales.py`**

CLI: `python compare_dales.py <case> <udales_run_dir> <dales_results_dir|none> <out_dir>`.
- It reads `thermostat_t/z.30x.nc`, and DALES `tmser.001.nc`/`profiles.001.nc` when given. DALES variable names include `lwp_bar`, `cfrac`, `zi`, `ustar`, `thl`, `qt`, `ql`; check the actual names in the files.
- It writes PNGs: time series of `lwp_bar`, `cfrac`, `zi`, `ustar`; mean profiles of `thl`, `qt`, `ql`, `wthl` over the paper's averaging window.
- Published ranges are shaded bands: a dict of `{(case, variable): (t0, t1, lo, hi, source)}`. Fill it only from values read in the papers or tables, and leave out any band whose number cannot be cited. A missing band is better than an invented one.
- It writes `summary.md` with a pass/fail table.
- Entrainment rate for 302/303: `we = dzi/dt + D·zi` (D = 0 for smoke, 3.75e−6 for RF01), using `zi_smoke`/`zi_qt8`.

- [ ] **Step 6: Commit**

```bash
git add tools/python/examples/dales_cases examples/301 examples/302 examples/303
git status
git commit -m "Add BOMEX, smoke and DYCOMS RF01 cases with DALES comparison tooling

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Tier-2 runs on CX3 and report

**Files:**
- Create: `examples/30x/RESULTS.md` (short, with plots referenced from a run directory outside the repo; do not commit large outputs)

- [ ] **Step 1: Build Release and submit**

```bash
./tools/build_executable.sh icl release 2>&1 | tail -3
for c in 301 302 303; do tools/hpc_execute.sh examples/$c; done
qstat -u $USER
```
Use `export I_MPI_OFI_PROVIDER=mlx` in the job environment (the CX3 fabric caveat in memory) for multi-node jobs. Check each job log for `libfabric provider:` with `I_MPI_DEBUG=5` before trusting timings. The performance of these runs is not the point; correctness is.

- [ ] **Step 2: Compare and write the report**

Run `compare_dales.py` for each case with the DALES `results_3.2.0` (301, 302) or `none` (303). Write `RESULTS.md` per case:
- pass/fail table;
- one paragraph on agreement;
- any deviation, with its most likely cause. Diagnose before blaming the port: surface `u*` first, then SGS, then advection.

- [ ] **Step 3: Commit**

```bash
git add examples/301/RESULTS.md examples/302/RESULTS.md examples/303/RESULTS.md
git status
git commit -m "Tier-2 validation results for BOMEX, smoke and RF01

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 9: Final verification and review

- [ ] Rebuild Debug and Release from clean (`rm -rf build/debug build/release`); both must compile without new warnings compared with the baseline build log.
- [ ] `python tests/run_tests.py python-library` and `python tests/run_tests.py supported --branch-a master --branch-b HEAD --build-type Debug`: all as in Task 6.
- [ ] Bit-for-bit check: run `tests/cases/526` and `tests/regression/david_tests` case 103 on the commit before Task 2 and on HEAD; compare restart files byte for byte. Expected: identical. Only Tasks 1 and 3-D1 may change results, and Task 2 onwards must not.
- [ ] Code review of the whole branch diff against `origin/master` by an independent reviewer (Fable), using the `requesting-code-review` skill. Fix confirmed findings and re-run the tests.
- [ ] Controller pushes `thermo` and opens a draft PR. Ask the user before pushing.

---

## Addendum (2026-11-13): align uDALES numerics with DALES for BOMEX and RF01

User decision: BOMEX and RF01 must match DALES closely. The audits are in `.superpowers/sdd/2026-11-13-moist-thermo-port/numerics-audit-{dynamics,sgs}.md`. The fixes below are applied **outright**: defaults change, so results change for all users. Each is recorded under **Breaking changes** in `CHANGELOG.md`, with the issue it closes. DALES reference: dev `075ac4e`.

Addendum global constraints:
- Every behaviour change gets a CHANGELOG "Breaking changes" entry: what changed, why, who is affected.
- Existing tests whose reference values change get updated, and each update is justified in the commit message. A test that starts failing for an unexplained reason is a stop.
- Do not touch `build/release` until the controller says the reference runs are finished; build and test with `build/debug`.
- Follow the surrounding code idiom. MPI is called through implicit externals (ifort #6580).

### Task 10: One-equation SGS model aligned with DALES (C1, C3, C4, C5)

**Files:** src/modsubgrid.f90, src/modsubgriddata.f90, tests (new or updated), CHANGELOG.md, docs (SGS section if one describes these constants).

- C1: the turbulent Prandtl convention. Set `ch = 1.0/Prandtl` (DALES modsubgrid.f90:74). Keep the default `Prandtl = 1/3`, so Kh = 3·Km in unstable/neutral cells, and drop the commented `Prandtl = 3.0` line. Grep every other use of `Prandtl`/`prandtli`/`prandtlturb` (Smagorinsky, Vreman, wall functions, statistics) and confirm each one already uses the Kh = Km/Pr convention. Report each use, and fix any that do not.
- C3: remove the factor `2.` in `sbdiss` ("added factor 2. for shear-driven flow"), to match DALES `- (ce1 + ce2*zlt/delta) * e120**2 /(2.*zlt)`. Keep the near-wall `damp`, which is uDALES-specific and 1 away from walls.
- C4:
  - Initialise `ekm, ekh, zlt, sbdiss, sbshr, sbbuo` to 0 after allocation, as DALES does at modsubgrid.f90:68.
  - When the MO surface is active (`lmo_surface()` from modsurface), compute the TKE sources at `kb` as DALES's "special treatment for lowest full level: k=1" block in `sources` (DALES modsubgrid.f90:588 onwards). It uses the MO gradients `dudz`/`dvdz` (uDALES equivalents: the MO module's surface gradients; add them to modsurface if they are not stored) and `dthvdz(kb)`.
  - Otherwise keep today's behaviour at `kb`: no sources, because the wall function sets the near-wall state. Apply this explicitly, without relying on uninitialised memory.
- C5: use `grav/thvf(k)` instead of `grav/thvs` in the stable length scale `zlt` and in `sbbuo`, as DALES does. Keep the `max(abs(dthvdz),eps1)` guard (#375).
- Tests:
  - A unit-style solver test that, for a known stratified state, checks the ekh/ekm ratio is 3 in unstable cells and `ch1 + ch2*zlt/delta` in stable ones.
  - A check that e12p(kb) is deterministic and 0-source without MO.
  - All existing supported suites pass. The tke_top_bc and fielddump tests may need reference updates; justify them.

### Task 11: Latent heat, local subsidence, thermostat timing (C2, C6, C8)

**Files:** src/modglobal.f90, src/modforces.f90, src/program.f90 (thermostat call position), tools/python/examples/dales_cases/make_rf01.py, tests/integration/thermo/test_thermo.py, examples/303 (regenerate inputs if the rlv change affects them), CHANGELOG.md.

- C2: set `rlv = 2.53e6` (the DALES value). Grep every literal `2.26e6` in src, tools, tests and docs, and every derived constant (`rlv/cp`, wall-moisture energy balance, modEB latent-heat flux), and update them consistently. Regenerate any committed `*.inp.303` that depends on it via make_rf01.py, and report which files changed.
- C6: subsidence on the local field, as DALES `lstend` (modforces.f90:298-414). Use the DALES upwind stencil on `thl0`, `qt0`, `sv0` (and `u0`, `v0` if `lmomsubs`), per cell, instead of the slab means `thl0av` etc. Keep uDALES's kb treatment, and say whether DALES does the same. Mask with IBM as the surrounding code does, if a mask is used there.
- C8: call `thermostat` after `thermodynamics`, so that ql and rhof are current. Check that the restart/tnext logic is unaffected.
- Tests: the thermo suite must pass, including the radiation balance test (which involves subsidence). Update the RLV constant in test_thermo.py.

### Task 12: Decomposition-independent random perturbations with thl/qt noise (C7, #395)

**Files:** src/modstartup.f90 (randomize_field and its callers), tests, CHANGELOG.md, docs (namoptions: randthl, randqt, krand).

- Replace the one-step LCG on `seed + linear_index`, which produces a sawtooth plane wave, with a counter-based hash. It must give a proper uniform deviate as a function of `(irandom, field_tag, iglob, jglob, k)`, for example a splitmix64- or murmur-style mix in integer(8) with masking to stay within the signed range. Keep the output independent of the MPI decomposition.
- Give each field its own `field_tag`, so that u, v, w, thl and qt noise are independent.
- Apply `randthl` to thl and `randqt` to qt (both the m and 0 fields) on cold start, up to `krand`, under the same guard as `randu`. Check how DALES `randomize_new` (modstartup.f90:1667) chooses levels and fields, and match it where the uDALES namelist allows.
- Tests:
  - decomposition independence (serial vs 2×2, bit-identical initial fields);
  - statistics of the deviates (mean about 0, variance of a uniform distribution, spectrum without a dominant wavenumber, correlation between field tags about 0), as a Python check on a dumped initial field;
  - thl/qt noise is present when randthl/randqt > 0.
- #395 is closed by this task.

### Task 13: DALES-compatible case settings and documentation

**Files:** examples/301, 302, 303 (namoptions, README, config.sh), tools/python/examples/dales_cases/README.md, docs/udales-namoptions-overview.md (defaults changed in Tasks 10-12).

- 301: `xlat = 15`; DALES BOMEX noise (`randthl = 0.1`, `randqt = 2.5e-5`, `randu = 0` (DALES applies no u/v/w noise), `krand = kmax`); `lqlnr = .true.`.
  - Time step: choose `courant`/`diffnr`/`dtmax` so that the uDALES dt matches the DALES run's dt (read it from /rds/general/user/mvr/ephemeral/dalesbomex tmser). Note that uDALES uses the L1 sum and DALES the L2 norm. Record the resulting dt in the README.
- 303 and 302: the equivalent of the above, for 303 against the DALES RF01 reference run (/rds/general/user/mvr/ephemeral/refruns).
- Fix the 301 README claim that 0.5 is the DALES courant default.
- Short tests: run each case for 60 s and check that it starts with no NaN.

### Task 14: Close-match validation

- Rebuild Release once the controller confirms it may.
- Run uDALES 301 and 303 (full length) with the Task 13 settings.
- Run DALES in a uDALES-compatible configuration: BOMEX and RF01 with a constant-density base profile (`ibas_prf = 5` with a constant-density baseprof.inp, giving Boussinesq-equivalent forcing) and BOMEX with `cu = 0`.
- Fix the comparison diagnostics:
  - the zi definitions: compute zi from the same definition (max dthl/dz on the mean profile, and the qt 8 g/kg threshold for RF01) for both codes, from profiles, in compare_dales.py;
  - the window alignment, i.e. the 10-minute block offset.
- Compare every pair against the seed-to-seed noise floor of each code from refruns-report.md. A difference counts as small when it is within the noise floor.
- Update RESULTS.md for 301 and 303 with the tables.
- Report any remaining difference above the noise floor, with a hypothesis.

### Task 15: Documented comparison of uDALES and DALES numerics

**Files:** new docs/udales-dales-numerics.md (linked from mkdocs.yml and from the tier-2 case READMEs).

A maintained reference describing, item by item, where uDALES numerics agree with and differ from classical DALES (dev 075ac4e, with notes where DALES 4.x releases differ). It should cover:
- governing equations (Boussinesq vs anelastic, base state, Exner);
- grid and staggering;
- advection schemes available in each code and the ones used by default and in the GCSS cases;
- time integration and time-step control (CFL norm, diffusion number vs Peclet);
- pressure solver;
- SGS models (one-equation TKE constants, length scale, wall damping, near-surface treatment, Smagorinsky and Vreman);
- surface schemes (MO vs wall functions);
- thermodynamics (saturation adjustment, es, constants);
- radiation;
- large-scale forcing and subsidence;
- the Galilean transform, sponge layer and boundary conditions;
- random perturbations;
- diagnostics definitions (zi, LWP, cloud cover).
For each item, give both formulas with file:line, the status after this branch (aligned / different by design / different, open issue #), and the measured effect on BOMEX/RF01 where Task 14 quantified it. Build it from numerics-audit-dynamics.md, numerics-audit-sgs.md and the Task 14 results, re-checking every row against the final code. Each claim of "identical" must be verified, not carried over from the audit.
