#!/usr/bin/env python3
"""fielddump diagnostics that are computed at output time must be consistent.

The divergence of the corrected velocity (``di``) and its terms (``ux``,
``vy``, ``wz``), and the right-hand side of the Poisson equation (``rs``)
and its terms (``du``, ``dv``, ``dw``), are not stored by the model.  They
are rebuilt from ``um, vm, wm`` and from the predicted velocities
``pup, pvp, pwp`` when the instantaneous field output writes them
(``instant.compute_derived`` in out_instants.f90 and
``modpois.poisson_rhs_diagnostic``).  This test runs the small case 300
on a 2x2 decomposition with every such code requested and checks that

* each derivative field equals the finite difference of the field it is
  derived from, as written in the same record (``dudx`` against ``u``,
  ``dpupdx`` against ``pup``, and so on), which pins the stencil, the
  grid spacing and the staggering;
* each divergence equals the sum of its three written terms;
* every value is finite, including the first record;
* the solution does not depend on whether these fields are requested: the
  restart files of a run with the diagnostics and of a run without any
  fielddump are byte-identical.

The output file holds single precision, so all tolerances are float32-scale.

Environment:
  UDALES_BUILD            path to the u-dales executable
                          (default build/release/u-dales)
  UDALES_RUNTIME_MODULES  module stack loaded before the run
  MPIEXEC                 MPI launcher
"""

from __future__ import annotations

import filecmp
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Dict, List

import numpy as np
from netCDF4 import Dataset

REPO_ROOT = Path(__file__).resolve().parents[3]
CASE_ID = 300
CASE_DIR = REPO_ROOT / "tests" / "cases" / str(CASE_ID)

UDALES_BUILD = Path(os.environ.get("UDALES_BUILD", REPO_ROOT / "build" / "release" / "u-dales"))
RUNTIME_MODULES = os.environ.get(
    "UDALES_RUNTIME_MODULES",
    "intel/2025a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a "
    "FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0",
)

NPROCX, NPROCY = 2, 2
RUNTIME = 2.0
DTMAX = 0.2
TFIELDDUMP = 1.0
#: Every output-time diagnostic plus the fields they are derived from.
FIELDVARS = "u0,v0,w0,di,ux,vy,wz,rs,du,dv,dw,pd,pu,pv,pw"
VARIABLES = ["u", "v", "w", "div", "dudx", "dvdy", "dwdz",
             "rhs", "dpupdx", "dpvpdy", "dpwpdz", "p", "pup", "pvp", "pwp"]

#: float32 output: relative tolerance for the finite-difference checks and
#: for the "divergence equals the sum of its terms" identities.
RTOL = 1.0e-5

sys.path.insert(0, str(REPO_ROOT / "tests" / "integration"))
import mpi_env  # noqa: E402


# --------------------------------------------------------------------------- #
# Running the solver
# --------------------------------------------------------------------------- #


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
    mpi_env.preflight(nprocs)
    mpiexec, extra_args = mpi_env.mpi_exec_and_args()
    command = (
        f"{_shell_prefix()}cd '{run_dir}' && "
        f"'{mpiexec}' {extra_args} -n {nprocs} '{UDALES_BUILD}' namoptions.{CASE_ID}"
    )
    return subprocess.run(
        ["bash", "-lc", command], cwd=REPO_ROOT, check=False, capture_output=True, text=True
    )


def _write_case(run_dir: Path, with_dump: bool) -> None:
    """Case 300 on 2x2 ranks, with or without the fielddump of every diagnostic."""
    shutil.copytree(CASE_DIR, run_dir, dirs_exist_ok=True)
    path = run_dir / f"namoptions.{CASE_ID}"
    text = path.read_text(encoding="utf-8")
    edits = {
        "runtime": f"{RUNTIME:g}",
        "trestart": f"{RUNTIME:g}",
        "dtmax": f"{DTMAX:g}",
        "nprocx": str(NPROCX),
        "nprocy": str(NPROCY),
        "lfielddump": ".true." if with_dump else ".false.",
    }
    for key, value in edits.items():
        text, n = re.subn(rf"(?m)^(\s*{key}\s*=\s*).*$", lambda m: m.group(1) + value, text)
        if n != 1:
            raise RuntimeError(f"expected exactly one '{key}' setting in namoptions.{CASE_ID}, found {n}")
    for key in ("tinstantdump", "fieldvars"):
        if re.search(rf"(?m)^\s*{key}\s*=", text):
            raise RuntimeError(f"namoptions.{CASE_ID} already sets '{key}'; the test expects to add it")
    text, n = re.subn(
        r"(?m)^(\s*lfielddump\s*=.*)$",
        lambda m: m.group(1) + f"\ntinstantdump = {TFIELDDUMP:g}\nfieldvars    = '{FIELDVARS}'",
        text,
    )
    if n != 1:
        raise RuntimeError("could not add tinstantdump/fieldvars after lfielddump")
    path.write_text(text, encoding="utf-8")


def _read_real(run_dir: Path, key: str) -> float:
    text = (run_dir / f"namoptions.{CASE_ID}").read_text(encoding="utf-8")
    match = re.search(rf"(?m)^\s*{key}\s*=\s*([0-9.eE+-]+)", text)
    if not match:
        raise RuntimeError(f"no numeric setting '{key}'")
    return float(match.group(1))


# --------------------------------------------------------------------------- #
# Reading the dump
# --------------------------------------------------------------------------- #


class _Rank:
    """One x-processor column's ins_field file (y-merged), every variable as float64 (time, z, y, x)."""

    def __init__(self, path: Path) -> None:
        self.path = path
        with Dataset(path) as nc:
            missing = [v for v in VARIABLES if v not in nc.variables]
            if missing:
                raise RuntimeError(f"{path.name} lacks {missing}")
            self.zm = np.asarray(nc.variables["zm"][:], dtype=np.float64)
            self.fields: Dict[str, np.ndarray] = {}
            for name in VARIABLES:
                var = nc.variables[name]
                order = [next(i for i, d in enumerate(var.dimensions) if d.startswith(axis))
                         for axis in ("time", "z", "y", "x")]
                self.fields[name] = np.transpose(np.asarray(var[:], dtype=np.float64), order)
        self.nt = self.fields["u"].shape[0]


def _ranks(run_dir: Path) -> List[_Rank]:
    # out_instants writes one y-merged file per x-processor column
    files = sorted(run_dir.glob(f"ins_field.*.{CASE_ID}.nc"))
    if len(files) != NPROCX:
        raise RuntimeError(f"expected {NPROCX} ins_field files in {run_dir}, found {len(files)}")
    return [_Rank(f) for f in files]


# --------------------------------------------------------------------------- #
# Tests
# --------------------------------------------------------------------------- #


class FielddumpDiagnostics(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        if not UDALES_BUILD.is_file():
            raise RuntimeError(f"u-dales executable not found at {UDALES_BUILD}")
        cls._tmp = tempfile.TemporaryDirectory(prefix="udales-fielddump-diag-")
        root = Path(cls._tmp.name)
        cls.with_dump = root / "with_dump"
        cls.no_dump = root / "no_dump"
        for run_dir, flag in ((cls.with_dump, True), (cls.no_dump, False)):
            _write_case(run_dir, flag)
            done = _run(run_dir, NPROCX * NPROCY)
            if done.returncode != 0:
                tail = "\n".join((done.stdout or "").strip().splitlines()[-40:])
                raise RuntimeError(f"run in {run_dir.name} exited {done.returncode}\n{tail}\n{done.stderr[-2000:]}")
        cls.ranks = _ranks(cls.with_dump)
        cls.dx = _read_real(cls.with_dump, "xlen") / _read_real(cls.with_dump, "itot")
        cls.dy = _read_real(cls.with_dump, "ylen") / _read_real(cls.with_dump, "jtot")

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    # -- helpers ----------------------------------------------------------- #

    def _assert_close(self, what: str, got: np.ndarray, want: np.ndarray, *scale_from: np.ndarray) -> None:
        scale = max(float(np.max(np.abs(a))) for a in scale_from)
        self.assertGreater(scale, 0.0, f"{what}: the reference field is identically zero, the check is void")
        err = float(np.max(np.abs(got - want)))
        self.assertLessEqual(err, RTOL * scale, f"{what}: max abs error {err:.3e} exceeds {RTOL:g} x {scale:.3e}")

    def _x_derivative(self, rank: _Rank, field: str, derivative: str) -> None:
        f, d = rank.fields[field], rank.fields[derivative]
        want = (f[:, :, :, 1:] - f[:, :, :, :-1]) / self.dx
        self._assert_close(f"{derivative} vs {field} ({rank.path.name})", d[:, :, :, :-1], want, want)

    def _y_derivative(self, rank: _Rank, field: str, derivative: str) -> None:
        f, d = rank.fields[field], rank.fields[derivative]
        want = (f[:, :, 1:, :] - f[:, :, :-1, :]) / self.dy
        self._assert_close(f"{derivative} vs {field} ({rank.path.name})", d[:, :, :-1, :], want, want)

    def _z_derivative(self, rank: _Rank, field: str, derivative: str) -> None:
        # zm holds the half levels; dzf(k) = zm(k+1) - zm(k) spans cell k.
        dzf = (rank.zm[1:] - rank.zm[:-1])[None, :, None, None]
        f, d = rank.fields[field], rank.fields[derivative]
        want = (f[:, 1:, :, :] - f[:, :-1, :, :]) / dzf
        self._assert_close(f"{derivative} vs {field} ({rank.path.name})", d[:, :-1, :, :], want, want)

    # -- tests ------------------------------------------------------------- #

    def test_records_and_finite_values(self) -> None:
        for rank in self.ranks:
            self.assertGreaterEqual(rank.nt, 2, f"{rank.path.name}: fewer than two records")
            for name, values in rank.fields.items():
                self.assertTrue(np.all(np.isfinite(values)), f"{rank.path.name}: {name} has non-finite values")

    def test_velocity_derivatives_match_the_written_velocity(self) -> None:
        for rank in self.ranks:
            self._x_derivative(rank, "u", "dudx")
            self._y_derivative(rank, "v", "dvdy")
            self._z_derivative(rank, "w", "dwdz")

    def test_divergence_is_the_sum_of_its_terms(self) -> None:
        for rank in self.ranks:
            f = rank.fields
            self._assert_close(f"div vs dudx+dvdy+dwdz ({rank.path.name})",
                               f["div"], f["dudx"] + f["dvdy"] + f["dwdz"], f["dudx"], f["dvdy"], f["dwdz"])

    def test_poisson_terms_match_the_written_predicted_velocity(self) -> None:
        for rank in self.ranks:
            self._x_derivative(rank, "pup", "dpupdx")
            self._y_derivative(rank, "pvp", "dpvpdy")
            self._z_derivative(rank, "pwp", "dpwpdz")

    def test_poisson_rhs_is_the_sum_of_its_terms(self) -> None:
        for rank in self.ranks:
            f = rank.fields
            self._assert_close(f"rhs vs dpupdx+dpvpdy+dpwpdz ({rank.path.name})",
                               f["rhs"], f["dpupdx"] + f["dpvpdy"] + f["dpwpdz"],
                               f["dpupdx"], f["dpvpdy"], f["dpwpdz"])

    def test_diagnostics_do_not_change_the_solution(self) -> None:
        with_dump = sorted(p.name for p in self.with_dump.glob(f"initd*.{CASE_ID}"))
        no_dump = sorted(p.name for p in self.no_dump.glob(f"initd*.{CASE_ID}"))
        self.assertEqual(len(with_dump), NPROCX * NPROCY, "restart files missing from the run with fielddump")
        self.assertEqual(with_dump, no_dump, "the two runs wrote different restart files")
        for name in with_dump:
            self.assertTrue(filecmp.cmp(self.with_dump / name, self.no_dump / name, shallow=False),
                            f"{name} differs between the run with and without the fielddump diagnostics")


if __name__ == "__main__":
    unittest.main(verbosity=2)
