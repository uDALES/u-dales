"""Launching the solver from the nesting test drivers.

Shared by ``test_nesting.py`` (the unit runmodes) and ``test_nesting_cases.py``
(the integration matrix), which used to carry verbatim copies of this code.

Everything here is configured from the environment, with documented defaults:

``UDALES_RUNTIME_MODULES``
    Module stack loaded in the shell that runs the solver, when a ``module``
    command exists (a GitHub runner has none, so the value is ignored there).
    Default: the CX3 Intel 2021a stack that ``tools/build_executable.sh icl``
    builds with, see ``DEFAULT_RUNTIME_MODULES``.  Set it to the empty string
    to load nothing.  For a gfortran build on CX3 use the ``foss/2023a`` stack
    recorded in ``.github/skills/udales-exec/references/clusters.md``.
``UDALES_MPIEXEC``
    The MPI launcher.  Falls back to ``MPIEXEC`` (which
    ``.github/scripts/setup_mpi_env.sh`` exports in CI), then to the
    ``mpiexec`` next to ``mpiifort`` when the Intel stack is on ``PATH``, then
    to a bare ``mpiexec`` resolved in the run shell *after* the modules load.
``MPI_LAUNCH_EXTRA_ARGS``
    Extra launcher arguments.  ``--oversubscribe`` is added automatically when
    the launcher identifies itself as Open MPI, so that 2x2 runs work on a
    2-4 core runner or a login node.
``UDALES_REQUIRE_LAUNCHER``
    ``1`` makes an unusable launcher a test **failure** instead of a skip.  The
    manifest sets it for every suite in the merge gate: a suite that skips
    itself because ``module load`` failed is a silent pass, and the gate must
    not contain those.  Unset for ad-hoc runs, where a skip is the useful
    behaviour.
``TMPDIR``
    Where run directories are created (through ``tempfile.gettempdir()``).  On
    a shared cluster point it at scratch rather than ``/tmp``.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from typing import Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parents[3]

#: The CX3 Intel 2021a stack of ``tools/build_executable.sh icl``, applied only
#: when a ``module`` command exists.  Override with ``UDALES_RUNTIME_MODULES``.
DEFAULT_RUNTIME_MODULES = (
    "intel/2021a netCDF/4.8.0-iimpi-2021a netCDF-Fortran/4.5.3-iimpi-2021a "
    "FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0"
)

RUNTIME_MODULES = os.environ.get("UDALES_RUNTIME_MODULES", DEFAULT_RUNTIME_MODULES)


def scratch_dir() -> Path:
    """Where temporary run directories go: ``TMPDIR`` if set, else the system default."""
    return Path(tempfile.gettempdir())


def find_mpiexec() -> str:
    """The launcher, by the precedence documented in the module header."""
    for var in ("UDALES_MPIEXEC", "MPIEXEC"):
        value = os.environ.get(var, "").strip()
        if value:
            return value
    mpiifort = shutil.which("mpiifort")
    if mpiifort:
        return str(Path(mpiifort).parent / "mpiexec")
    return "mpiexec"


def mpi_exec_and_args() -> Tuple[str, str]:
    mpiexec = find_mpiexec()
    extra_args = os.environ.get("MPI_LAUNCH_EXTRA_ARGS", "").strip()
    try:
        version = subprocess.run(
            [mpiexec, "--version"], check=False, capture_output=True, text=True
        ).stdout
    except OSError:
        version = ""
    if (re.search(r"Open MPI|OpenRTE", version, flags=re.IGNORECASE)
            and "--oversubscribe" not in extra_args):
        extra_args = f"--oversubscribe {extra_args}".strip()
    return mpiexec, extra_args


def shell_prefix(quiet_intel_diagnostics: bool = False) -> str:
    """Shell preamble for a solver run: module stack and HDF5 locking.

    ``quiet_intel_diagnostics`` sets ``FOR_DISABLE_DIAGNOSTIC_DISPLAY=TRUE``,
    which silences the Intel runtime's ``forrtl: severe`` banner and traceback
    (a Debug build otherwise prints them on every trap).  Leave it off when the
    traceback is the thing being looked for.
    """
    prefix = ""
    if Path("/etc/profile.d/modules.sh").is_file():
        prefix = "source /etc/profile.d/modules.sh >/dev/null 2>&1 || true; "
    if RUNTIME_MODULES.strip():
        prefix += (f"if command -v module >/dev/null 2>&1; then "
                   f"module load {RUNTIME_MODULES}; fi && ")
    prefix += "export HDF5_USE_FILE_LOCKING=FALSE && "
    if quiet_intel_diagnostics:
        prefix += "export FOR_DISABLE_DIAGNOSTIC_DISPLAY=TRUE && "
    return prefix


def run_in_shell(command: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["bash", "-lc", command], cwd=REPO_ROOT, check=False, capture_output=True, text=True
    )


def run(run_dir: Path, namelist: str, nprocs: int, executable: Path,
        quiet_intel_diagnostics: bool = False) -> subprocess.CompletedProcess:
    """Run ``executable namelist`` on ``nprocs`` ranks inside ``run_dir``."""
    mpiexec, extra_args = mpi_exec_and_args()
    command = (
        f"{shell_prefix(quiet_intel_diagnostics)}cd '{run_dir}' && "
        f"'{mpiexec}' {extra_args} -n {nprocs} '{executable}' {namelist}"
    )
    return run_in_shell(command)


def launcher_unavailable_reason() -> Optional[str]:
    """``None`` when ``mpiexec -n 1 /bin/true`` works; otherwise why it did not."""
    mpiexec, extra_args = mpi_exec_and_args()
    probe = run_in_shell(f"{shell_prefix()}'{mpiexec}' {extra_args} -n 1 /bin/true")
    if probe.returncode == 0:
        return None
    return (probe.stderr or probe.stdout or f"exit code {probe.returncode}").strip()


def require_launcher() -> None:
    """Skip, or fail, when the launcher cannot start a one-rank job.

    A failed ``module load`` or a missing ``mpiexec`` used to turn a whole
    suite into ``SkipTest``, which ``run_tests.py`` reports as PASS.  Under
    ``UDALES_REQUIRE_LAUNCHER=1`` -- set by the manifest for the supported
    groups -- that is an error instead, so the merge gate cannot pass on a
    suite that never ran.
    """
    reason = launcher_unavailable_reason()
    if reason is None:
        return
    mpiexec, extra_args = mpi_exec_and_args()
    message = (f"MPI launcher is not usable here ({mpiexec} {extra_args}".strip()
               + f"): {reason}")
    if os.environ.get("UDALES_REQUIRE_LAUNCHER", "") == "1":
        raise RuntimeError(
            "UDALES_REQUIRE_LAUNCHER=1, so an unusable launcher is a failure, "
            "not a skip: " + message)
    raise unittest.SkipTest(message)


def tail(label: str, text: str, limit: int = 40) -> str:
    stripped = (text or "").strip()
    if not stripped:
        return f"{label}: <empty>"
    lines = stripped.splitlines()
    return f"{label} (last {min(len(lines), limit)} lines):\n" + "\n".join(lines[-limit:])
