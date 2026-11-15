#!/usr/bin/env python3
"""The geostrophic pressure gradient must only act together with the Coriolis force.

`modstartup` builds the large-scale pressure gradient from the geostrophic wind
as ``dpdxl = om23_gs*vg - pgx - dpdx`` and ``dpdyl = -om23_gs*ug - pgy``.  The
``om23_gs`` terms are the geostrophic balance ``f k x u_g``: they only make
sense as the counterpart of the Coriolis acceleration on the resolved wind.
With ``lcoriol = .false.`` there is no such counterpart, so a nonzero ``ug``
accelerated the whole domain at ``f*|u_g|`` (issue #410).

The test runs a laminar column -- case 300 started from rest, no buildings, no
temperature or moisture gradient -- with ``ug = UG`` in ``lscale.inp`` and
reads the horizontal-mean wind back out of the restart file.

* ``lcoriol = .false.``: the wind must stay at rest.
* ``lcoriol = .true.``: unchanged behaviour.  From rest there is no Coriolis
  force yet, so the geostrophic gradient accelerates the flow at ``f*ug`` in
  the v direction, and ``v = om23_gs*ug*t`` to well within 1 %.

Environment:
  UDALES_BUILD            path to the u-dales executable
                          (default build/release/u-dales)
  UDALES_RUNTIME_MODULES  module stack loaded before the run
  MPIEXEC                 MPI launcher
"""

from __future__ import annotations

import math
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
CASE_ID = 300
CASE_DIR = REPO_ROOT / "tests" / "cases" / str(CASE_ID)

UDALES_BUILD = Path(os.environ.get("UDALES_BUILD", REPO_ROOT / "build" / "release" / "u-dales"))
RUNTIME_MODULES = os.environ.get(
    "UDALES_RUNTIME_MODULES",
    "intel/2021a netCDF/4.8.0-iimpi-2021a netCDF-Fortran/4.5.3-iimpi-2021a "
    "FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0",
)

#: Restart record indices (modsave.f90): u0, v0, then the (timee, dt) pair.
U0_RECORD, V0_RECORD, TIME_RECORD = 2, 3, 12

UG = 10.0          # geostrophic u-wind [m/s]
RUNTIME = 2.0
DTMAX = 0.2
#: om23_gs = 2*omega_gs*sin(xlat); omega_gs and the default xlat are in modglobal.f90.
OM23_GS = 2.0 * 7.292e-5 * math.sin(math.radians(52.0))

sys.path.insert(0, str(REPO_ROOT / "tests" / "integration"))
import mpi_env  # noqa: E402


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


def _write_case(run_dir: Path, lcoriol: bool) -> None:
    """Case 300 from rest, with ``ug = UG`` and the Coriolis switch set."""
    shutil.copytree(CASE_DIR, run_dir, dirs_exist_ok=True)

    path = run_dir / f"namoptions.{CASE_ID}"
    text = path.read_text(encoding="utf-8")
    edits = {
        "runtime": f"{RUNTIME:g}",
        "trestart": f"{RUNTIME:g}",
        "dtmax": f"{DTMAX:g}",
        "u0": "0.",
    }
    for key, value in edits.items():
        text, n = re.subn(rf"(?m)^(\s*{key}\s*=\s*).*$", lambda m: m.group(1) + value, text)
        if n != 1:
            raise RuntimeError(f"expected exactly one '{key}' setting in namoptions.{CASE_ID}, found {n}")
    text = re.sub(r"(?m)^(\s*l\w*dump\s*=\s*).*$", r"\g<1>.false.", text)
    flag = ".true." if lcoriol else ".false."
    text, n = re.subn(r"(?m)^(\s*ltempeq\s*=\s*).*$", rf"\g<1>.true.\nlcoriol      = {flag}", text)
    if n != 1:
        raise RuntimeError("expected exactly one ltempeq setting in the &PHYSICS group")
    path.write_text(text, encoding="utf-8")

    # Start from rest: u column of prof.inp (z thl qt u v e12) set to zero.
    prof = run_dir / f"prof.inp.{CASE_ID}"
    lines = prof.read_text(encoding="utf-8").splitlines()
    body = []
    for line in lines[2:]:
        if line.strip():
            cols = line.split()
            cols[3] = "0.000000"
            body.append("  ".join(cols))
    prof.write_text("\n".join(lines[:2] + body) + "\n", encoding="utf-8")

    # Geostrophic wind: ug = UG at every level, vg = 0, no other large-scale forcing.
    lscale = run_dir / f"lscale.inp.{CASE_ID}"
    lines = lscale.read_text(encoding="utf-8").splitlines()
    body = []
    for line in lines[2:]:
        if line.strip():
            cols = line.split()
            cols[1] = f"{UG:.6f}"
            body.append("  ".join(cols))
    lscale.write_text("\n".join(lines[:2] + body) + "\n", encoding="utf-8")


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


def _mean_wind(run_dir: Path) -> Tuple[float, float, float]:
    """Return (mean u, mean v, time) over the interior cells of the single rank."""
    nx = _read_int(run_dir, "itot") + 2
    ny = _read_int(run_dir, "jtot") + 2
    nz = _read_int(run_dir, "ktot") + 1
    files = sorted(run_dir.glob(f"initd*.{CASE_ID}"))
    if len(files) != 1:
        raise RuntimeError(f"expected one restart file in {run_dir}, found {len(files)}")
    recs = _records(files[0])
    u = np.frombuffer(recs[U0_RECORD], dtype="<f8").reshape((nz, ny, nx))[:-1, 1:-1, 1:-1]
    v = np.frombuffer(recs[V0_RECORD], dtype="<f8").reshape((nz, ny, nx))[:-1, 1:-1, 1:-1]
    t = float(np.frombuffer(recs[TIME_RECORD], dtype="<f8")[0])
    return float(u.mean()), float(v.mean()), t


class GeostrophicForcing(unittest.TestCase):
    _cache: Dict[bool, Tuple[float, float, float]] = {}

    @classmethod
    def setUpClass(cls) -> None:
        if not UDALES_BUILD.is_file():
            raise RuntimeError(f"u-dales executable not found at {UDALES_BUILD}")
        cls._tmp = tempfile.TemporaryDirectory(prefix="udales-geostrophic-")
        cls._root = Path(cls._tmp.name)

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def _wind(self, lcoriol: bool) -> Tuple[float, float, float]:
        if lcoriol not in self._cache:
            run_dir = self._root / ("coriol" if lcoriol else "nocoriol")
            _write_case(run_dir, lcoriol)
            done = _run(run_dir)
            if done.returncode != 0:
                tail = "\n".join((done.stdout or "").strip().splitlines()[-40:])
                self.fail(f"lcoriol = {lcoriol} exited {done.returncode}\n{tail}")
            self._cache[lcoriol] = _mean_wind(run_dir)
        return self._cache[lcoriol]

    def test_no_acceleration_without_coriolis(self) -> None:
        u, v, t = self._wind(False)
        expected_bug = OM23_GS * UG * t
        print(f"[lcoriol=F] t={t:g} s  mean u={u:.3e}  mean v={v:.3e}  "
              f"(the unbalanced gradient would give v={expected_bug:.3e})", flush=True)
        self.assertLessEqual(abs(u), 1.0e-4, "mean u accelerated with lcoriol = .false.")
        self.assertLessEqual(
            abs(v), 1.0e-4,
            f"mean v = {v:.3e} m/s after {t:g} s with lcoriol = .false. and ug = {UG:g}: "
            f"the geostrophic pressure gradient is acting without Coriolis "
            f"(f*ug*t = {expected_bug:.3e})",
        )

    def test_geostrophic_gradient_unchanged_with_coriolis(self) -> None:
        u, v, t = self._wind(True)
        expected = OM23_GS * UG * t
        print(f"[lcoriol=T] t={t:g} s  mean u={u:.3e}  mean v={v:.3e}  "
              f"expected v={expected:.3e}", flush=True)
        self.assertAlmostEqual(v / expected, 1.0, delta=0.01)


if __name__ == "__main__":
    unittest.main()
