#!/usr/bin/env python3
"""The chemistry update of ``lchem`` must be carried into the next time step.

``modchem::chem`` applies the NO/NO2/O3 null-cycle reaction step to ``sv0`` on
the third Runge-Kutta substep. The result has to end up in ``svm``, which is
the starting point of the next step; otherwise the reaction is discarded and
only a small residual of it survives through the following advection
tendency (issue #427).

The test runs the energy-balance case 064 (immersed boundaries with ground
facets, 64^3 cells) with three passive scalars initialised to spatially
uniform values, no scalar sources, periodic lateral boundaries and a fixed
time step. The facets impose zero scalar flux and advection and diffusion of
a uniform field vanish, so the scalars stay uniform and the solver's chemistry
reduces to the implicit map of ``chem`` applied once per time step with
``dt = dtmax``. The test reproduces that map in Python and checks that

* with ``lchem = .true.`` the dumped scalars equal the map applied as many
  times as there were time steps up to the dump (the bug of #427 leaves them
  within a few per cent of the initial values instead);
* the molar sums NO + NO2 and O3 - NO are conserved by the reaction step;
* the fields stay spatially uniform over the fluid cells;
* with ``lchem = .false.`` the scalars are unchanged.

Only fluid cells (``mask_c = 1`` in the dump) are compared; the solver does
not evolve the scalar inside solid cells. Single rank; float32 output, so the
tolerances are float32-scale.

Environment:
  UDALES_BUILD            path to the u-dales executable
                          (default build/release/u-dales)
  UDALES_RUNTIME_MODULES  module stack loaded before the run
  MPIEXEC                 MPI launcher
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

REPO_ROOT = Path(__file__).resolve().parents[3]
CASE_ID = "064"
CASE_DIR = REPO_ROOT / "tests" / "cases" / CASE_ID
UDALES_BUILD = Path(os.environ.get("UDALES_BUILD", REPO_ROOT / "build" / "release" / "u-dales"))
RUNTIME_MODULES = os.environ.get(
    "UDALES_RUNTIME_MODULES",
    "intel/2025a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a "
    "FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0",
)

RUNTIME = 2.0
DT = 0.1
TFIELDDUMP = 1.0
K1 = 0.5
JNO2 = 0.8
#: initial NO, NO2, O3 in the solver's mass units; molar masses as in modchem.f90
INITIAL = np.array([10.0, 20.0, 30.0])
MOLAR_MASS = np.array([30.006, 46.005, 47.997])
SCALARS = ["sca1", "sca2", "sca3"]

sys.path.insert(0, str(REPO_ROOT / "tests" / "integration"))
import mpi_env  # noqa: E402


def chem_step(sv: np.ndarray, dt: float) -> np.ndarray:
    """One implicit null-cycle step of modchem::chem on uniform fields."""
    no, no2, o3 = sv / MOLAR_MASS
    rate = dt * (-K1 * no * o3 + JNO2 * no2) / (1.0 + ((no + o3) * K1 + JNO2) * dt)
    return sv + MOLAR_MASS * np.array([rate, -rate, rate])


def _shell_prefix() -> str:
    prefix = ""
    if Path("/etc/profile.d/modules.sh").is_file():
        prefix = "source /etc/profile.d/modules.sh >/dev/null 2>&1 || true; "
    return (
        f"{prefix}"
        f"if command -v module >/dev/null 2>&1; then module load {RUNTIME_MODULES}; fi && "
        f"export HDF5_USE_FILE_LOCKING=FALSE && "
    )


def _run(run_dir: Path) -> subprocess.CompletedProcess:
    mpi_env.preflight(1)
    mpiexec, extra_args = mpi_env.mpi_exec_and_args()
    command = (
        f"{_shell_prefix()}cd '{run_dir}' && "
        f"'{mpiexec}' {extra_args} -n 1 '{UDALES_BUILD}' namoptions.{CASE_ID}"
    )
    return subprocess.run(
        ["bash", "-lc", command], cwd=REPO_ROOT, check=False, capture_output=True, text=True
    )


def _set(text: str, key: str, value: str) -> str:
    text, n = re.subn(rf"(?m)^(\s*{key}\s*=\s*).*$", lambda m: m.group(1) + value, text)
    if n != 1:
        raise RuntimeError(f"expected exactly one '{key}' setting in namoptions.{CASE_ID}, found {n}")
    return text


def _write_case(run_dir: Path, lchem: bool) -> None:
    """Case 064 on one rank with three uniform scalars, fixed time step, chemistry on or off."""
    shutil.copytree(CASE_DIR, run_dir, dirs_exist_ok=True)
    path = run_dir / f"namoptions.{CASE_ID}"
    text = path.read_text(encoding="utf-8")
    text = _set(text, "runtime", f"{RUNTIME:g}")
    text = _set(text, "trestart", f"{RUNTIME:g}")
    text = _set(text, "dtmax", f"{DT:g}")
    text = _set(text, "ladaptive", ".false.")
    text = _set(text, "nprocx", "1")
    text = _set(text, "nprocy", "1")
    text = _set(text, "lfielddump", ".true.")
    text = _set(text, "tfielddump", f"{TFIELDDUMP:g}")
    text = _set(text, "fieldvars", "'s1,s2,s3,mc'")
    text = _set(text, "lxytdump", ".false.")
    for key in ("nsv", "lchem"):
        if re.search(rf"(?m)^\s*{key}\s*=", text):
            raise RuntimeError(f"namoptions.{CASE_ID} already sets '{key}'; the test expects to add it")
    text += (
        "\n&SCALARS\n"
        "nsv          = 3\n"
        "/\n\n"
        "&CHEMISTRY\n"
        f"lchem        = {'.true.' if lchem else '.false.'}\n"
        f"k1           = {K1:g}\n"
        f"JNO2         = {JNO2:g}\n"
        "/\n"
    )
    path.write_text(text, encoding="utf-8")

    levels = []
    with (run_dir / f"prof.inp.{CASE_ID}").open(encoding="utf-8") as stream:
        for line in stream:
            if line.strip() and not line.lstrip().startswith("#"):
                levels.append(float(line.split()[0]))
    with (run_dir / f"scalar.inp.{CASE_ID}").open("w", encoding="utf-8") as stream:
        stream.write("# uniform NO, NO2, O3 for the lchem regression test\n")
        stream.write("# z sv1 sv2 sv3\n")
        for z in levels:
            stream.write(f"{z:20.12f} {INITIAL[0]:12.6f} {INITIAL[1]:12.6f} {INITIAL[2]:12.6f}\n")


def _last_record(run_dir: Path) -> tuple[float, np.ndarray, np.ndarray]:
    """Time, fluid-cell means and standard deviations of the last dumped record."""
    files = sorted(run_dir.glob(f"fielddump.*.{CASE_ID}.nc"))
    if len(files) != 1:
        raise RuntimeError(f"expected one fielddump file in {run_dir}, found {len(files)}")
    with Dataset(files[0]) as nc:
        nc.set_auto_mask(False)
        time = float(np.asarray(nc.variables["time"][:], dtype=np.float64)[-1])
        if "mask_c" not in nc.variables:
            raise RuntimeError(f"{files[0].name} lacks mask_c")
        fluid = np.asarray(nc.variables["mask_c"][-1, ...], dtype=np.float64) == 1.0
        if not fluid.any() or fluid.all():
            raise RuntimeError("expected both fluid and solid cells in the dump")
        means, spreads = [], []
        for name in SCALARS:
            if name not in nc.variables:
                raise RuntimeError(f"{files[0].name} lacks {name}")
            values = np.asarray(nc.variables[name][-1, ...], dtype=np.float64)[fluid]
            means.append(values.mean())
            spreads.append(values.std())
    return time, np.array(means), np.array(spreads)


class ChemistryUpdateIsKept(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not UDALES_BUILD.is_file():
            raise RuntimeError(f"u-dales executable not found at {UDALES_BUILD}")
        cls._tmp = tempfile.TemporaryDirectory(prefix="udales-lchem-")
        root = Path(cls._tmp.name)
        cls.results = {}
        for lchem in (True, False):
            run_dir = root / ("lchem_on" if lchem else "lchem_off")
            _write_case(run_dir, lchem)
            done = _run(run_dir)
            if done.returncode != 0:
                tail = "\n".join((done.stdout or "").strip().splitlines()[-40:])
                raise RuntimeError(f"run in {run_dir.name} exited {done.returncode}\n{tail}\n{done.stderr[-2000:]}")
            cls.results[lchem] = _last_record(run_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def test_fields_stay_uniform(self) -> None:
        for lchem, (_time, means, spreads) in self.results.items():
            with self.subTest(lchem=lchem):
                np.testing.assert_array_less(spreads, 1.0e-5 * np.abs(means))

    def test_without_chemistry_the_scalars_are_unchanged(self) -> None:
        _time, means, _spreads = self.results[False]
        np.testing.assert_allclose(means, INITIAL, rtol=1.0e-6, atol=0.0)

    def test_chemistry_update_is_carried_to_the_next_step(self) -> None:
        time, means, _spreads = self.results[True]
        n_steps = int(round(time / DT))
        self.assertGreaterEqual(n_steps, 10, f"dump at t = {time} s covers too few steps")
        self.assertAlmostEqual(time, n_steps * DT, delta=1.0e-6 * time)
        expected = INITIAL.copy()
        for _ in range(n_steps):
            expected = chem_step(expected, DT)
        # the reaction moves NO by about 50 per cent over the run
        np.testing.assert_allclose(means, expected, rtol=1.0e-5, atol=0.0)
        self.assertGreater(abs(expected[0] - INITIAL[0]) / INITIAL[0], 0.1)

    def test_reaction_step_conserves_molar_sums(self) -> None:
        _time, means, _spreads = self.results[True]
        molar = means / MOLAR_MASS
        initial = INITIAL / MOLAR_MASS
        self.assertAlmostEqual(molar[0] + molar[1], initial[0] + initial[1], delta=1.0e-5 * (initial[0] + initial[1]))
        self.assertAlmostEqual(molar[2] - molar[0], initial[2] - initial[0], delta=1.0e-5 * abs(initial[2] - initial[0]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
