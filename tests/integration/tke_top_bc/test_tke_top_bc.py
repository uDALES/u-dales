#!/usr/bin/env python3
"""The SGS TKE condition at the top of the domain must match the momentum one.

`BCtopm_freeslip` and `BCtopm_pressure` impose zero momentum flux through the
lid, so the subgrid TKE must be zero flux there too.  In `diffe` the top face
flux is proportional to ``e120(i,j,ke+1) - e120(i,j,ke)``, so "zero flux" is
exactly "the ghost equals the top interior cell".  `BCtopm_noslip` imposes a
velocity at the lid, so the TKE gets a Dirichlet condition at the model floor
`e12min` instead.

Before the fix, `freeslip` and `pressure` pinned the ghost to `e12min` -- a
numerical floor, not a boundary value -- which drained subgrid energy through
the lid on every step, and `noslip` set no condition at all.

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
import struct
import subprocess
import tempfile
import unittest
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
TEST_DIR = Path(__file__).resolve().parent
CASE_ID = 526
CASE_DIR = REPO_ROOT / "tests" / "cases" / str(CASE_ID)

UDALES_BUILD = Path(os.environ.get("UDALES_BUILD", REPO_ROOT / "build" / "release" / "u-dales"))
RUNTIME_MODULES = os.environ.get(
    "UDALES_RUNTIME_MODULES",
    "intel/2025a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a "
    "FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0",
)

#: Restart record index of e120, and the model's TKE floor (modglobal.f90).
E120_RECORD = 7
E12MIN = 5.0e-5

#: BCtopm values (modglobal.f90).
BCTOPM = {"freeslip": 1, "noslip": 2, "pressure": 3}

NPROCX, NPROCY = 2, 2
RUNTIME = 40.0
DTMAX = 0.2
#: Initial e12 profile.  It must sit well above E12MIN or the two conditions
#: are numerically indistinguishable and the test cannot fail.
E12_INIT = 0.1


# --------------------------------------------------------------------------- #
# Running the solver
# --------------------------------------------------------------------------- #


def _mpi_exec_and_args() -> Tuple[str, str]:
    mpiexec = os.environ.get("MPIEXEC")
    if not mpiexec:
        mpiifort = shutil.which("mpiifort")
        mpiexec = str(Path(mpiifort).parent / "mpiexec") if mpiifort else "mpiexec"
    extra_args = os.environ.get("MPI_LAUNCH_EXTRA_ARGS", "").strip()
    try:
        version = subprocess.run(
            [mpiexec, "--version"], check=False, capture_output=True, text=True
        ).stdout
    except OSError:
        version = ""
    if re.search(r"Open MPI|OpenRTE", version, flags=re.IGNORECASE) and "--oversubscribe" not in extra_args:
        extra_args = f"--oversubscribe {extra_args}".strip()
    return mpiexec, extra_args


def _shell_prefix() -> str:
    prefix = ""
    if Path("/etc/profile.d/modules.sh").is_file():
        prefix = "source /etc/profile.d/modules.sh >/dev/null 2>&1 || true; "
    return (
        f"{prefix}"
        f"if command -v module >/dev/null 2>&1; then module load {RUNTIME_MODULES}; fi && "
        f"export HDF5_USE_FILE_LOCKING=FALSE && "
    )


def _run(run_dir: Path, nprocs: int) -> subprocess.CompletedProcess:
    mpiexec, extra_args = _mpi_exec_and_args()
    command = (
        f"{_shell_prefix()}cd '{run_dir}' && "
        f"'{mpiexec}' {extra_args} -n {nprocs} '{UDALES_BUILD}' namoptions.{CASE_ID}"
    )
    return subprocess.run(
        ["bash", "-lc", command], cwd=REPO_ROOT, check=False, capture_output=True, text=True
    )


# --------------------------------------------------------------------------- #
# Case construction
# --------------------------------------------------------------------------- #


def _write_case(run_dir: Path, bctopm: int) -> None:
    """Case 526 switched to the one-equation SGS model with a given lid.

    ``ldelta = .true.`` is required: with the stability-dependent length scale
    the one-equation closure divides by ``sqrt(grav/thvs*|dthvdz|)``
    (modsubgrid.f90) and raises a floating-point exception on this case.  That
    is a separate defect and not what this test is about.
    """
    shutil.copytree(CASE_DIR, run_dir, dirs_exist_ok=True)

    text = (run_dir / f"namoptions.{CASE_ID}").read_text(encoding="utf-8")
    text = text.replace(
        "lvreman      = .true.",
        "lvreman      = .false.\nloneeqn      = .true.\nldelta       = .true.",
    )
    edits = {
        "runtime": f"{RUNTIME:g}",
        "trestart": f"{RUNTIME:g}",
        "dtmax": f"{DTMAX:g}",
        "BCtopm": str(bctopm),
        "nprocx": str(NPROCX),
        "nprocy": str(NPROCY),
    }
    for key, value in edits.items():
        text, n = re.subn(rf"(?m)^(\s*{key}\s*=\s*).*$", lambda m: m.group(1) + value, text)
        if n == 0:
            text = text.replace("loneeqn      = .true.", f"loneeqn      = .true.\n{key} = {value}")
    # Output costs wall-clock and disk and nothing here reads it.
    text = re.sub(r"(?m)^(\s*l\w*dump\s*=\s*).*$", r"\g<1>.false.", text)
    (run_dir / f"namoptions.{CASE_ID}").write_text(text, encoding="utf-8")

    prof = run_dir / f"prof.inp.{CASE_ID}"
    lines = prof.read_text(encoding="utf-8").splitlines()
    body = [
        "  ".join(cols[:5] + [f"{E12_INIT:.6f}"])
        for cols in (line.split() for line in lines[2:] if line.strip())
    ]
    prof.write_text("\n".join(lines[:2] + body) + "\n", encoding="utf-8")


# --------------------------------------------------------------------------- #
# Reading e120 back out of the restart file
# --------------------------------------------------------------------------- #


def _records(path: Path) -> List[bytes]:
    """Split an ifort sequential-unformatted file into its records."""
    data = path.read_bytes()
    out: List[bytes] = []
    off = 0
    while off < len(data):
        (n,) = struct.unpack("<i", data[off:off + 4])
        off += 4
        out.append(data[off:off + n])
        off += n + 4
    return out


def _read_int(run_dir: Path, key: str) -> int:
    text = (run_dir / f"namoptions.{CASE_ID}").read_text(encoding="utf-8")
    match = re.search(rf"(?m)^\s*{key}\s*=\s*(\d+)", text)
    if not match:
        raise RuntimeError(f"no integer setting '{key}'")
    return int(match.group(1))


def _e120_top_two(run_dir: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Return (e120 at ke, e120 at the ke+1 ghost), stacked over all ranks.

    The restart holds each field on ``(ib-1:ie+1, jb-1:je+1, kb:ke+1)``, so the
    ghost is the last k plane and the top interior cell the one below it.
    """
    nx = _read_int(run_dir, "itot") // NPROCX + 2
    ny = _read_int(run_dir, "jtot") // NPROCY + 2
    nz = _read_int(run_dir, "ktot") + 1

    top, ghost = [], []
    files = sorted(run_dir.glob(f"initd*.{CASE_ID}"))
    if not files:
        raise RuntimeError(f"no restart file written in {run_dir}")
    for path in files:
        block = np.frombuffer(_records(path)[E120_RECORD], dtype="<f8").reshape((nz, ny, nx))
        inner = block[:, 1:-1, 1:-1]        # drop the x/y halos
        top.append(inner[-2].ravel())
        ghost.append(inner[-1].ravel())
    return np.concatenate(top), np.concatenate(ghost)


# --------------------------------------------------------------------------- #
# Tests
# --------------------------------------------------------------------------- #


class TkeTopBoundaryCondition(unittest.TestCase):
    """One run per BCtopm; the ghost must carry the right condition."""

    _cache: Dict[str, Path] = {}

    @classmethod
    def setUpClass(cls) -> None:
        if not UDALES_BUILD.is_file():
            raise RuntimeError(f"u-dales executable not found at {UDALES_BUILD}")
        cls._tmp = tempfile.TemporaryDirectory(prefix="udales-tke-bc-")
        cls._root = Path(cls._tmp.name)

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def _run_for(self, lid: str) -> Path:
        if lid not in self._cache:
            run_dir = self._root / lid
            _write_case(run_dir, BCTOPM[lid])
            done = _run(run_dir, NPROCX * NPROCY)
            if done.returncode != 0:
                tail = "\n".join((done.stdout or "").strip().splitlines()[-40:])
                self.fail(f"BCtopm = {lid} exited {done.returncode}\n{tail}")
            self._cache[lid] = run_dir
        return self._cache[lid]

    def _assert_zero_flux(self, lid: str) -> None:
        top, ghost = self._e120(lid)
        worst = float(np.max(np.abs(ghost - top)))
        print(f"[{lid}] max |e12(ke+1) - e12(ke)| = {worst:.3e}   "
              f"mean e12(ke) = {top.mean():.6e}", flush=True)
        self.assertLessEqual(
            worst, 1.0e-14,
            f"{lid} is a zero-flux lid, so the SGS TKE ghost must equal the top "
            f"interior cell; worst mismatch {worst:.3e}",
        )

    def _e120(self, lid: str) -> Tuple[np.ndarray, np.ndarray]:
        return _e120_top_two(self._run_for(lid))

    # -- the two zero-flux lids ------------------------------------------- #

    def test_freeslip_lid_is_zero_flux_for_tke(self) -> None:
        self._assert_zero_flux("freeslip")

    def test_pressure_lid_is_zero_flux_for_tke(self) -> None:
        self._assert_zero_flux("pressure")

    # -- the Dirichlet lid ------------------------------------------------ #

    def test_noslip_lid_pins_tke_to_the_floor(self) -> None:
        """A solid lid takes the SGS TKE to zero, i.e. to the model floor.

        This branch previously set no condition at all, so the ghost held
        whatever it happened to contain.
        """
        _, ghost = self._e120("noslip")
        worst = float(np.max(np.abs(ghost - E12MIN)))
        print(f"[noslip] max |e12(ke+1) - e12min| = {worst:.3e}", flush=True)
        self.assertLessEqual(worst, 1.0e-14, "a no-slip lid must pin the SGS TKE ghost to e12min")

    # -- the test must be able to fail ------------------------------------ #

    def test_the_lid_carries_real_turbulence(self) -> None:
        """Guard against a vacuous pass.

        If the top interior cell sat at e12min anyway, the zero-flux and
        Dirichlet conditions would be indistinguishable and the assertions
        above would hold no matter what the code did.
        """
        top, _ = self._e120("freeslip")
        margin = float(top.mean() / E12MIN)
        print(f"[guard] mean e12(ke) / e12min = {margin:.1f}", flush=True)
        self.assertGreater(
            margin, 10.0,
            "the top cell must carry SGS TKE well above the floor for this test to mean anything",
        )


if __name__ == "__main__":
    unittest.main()
