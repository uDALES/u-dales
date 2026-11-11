"""Shared MPI launch settings for the Python test runners (mirror of tools/mpi_env.sh).

uDALES runs one MPI rank per CPU and never oversubscribes. ``preflight`` compares
the requested rank count with the CPUs granted to this process (the Slurm or PBS
allocation inside a batch job, otherwise what the operating system makes
available) and raises with both numbers when the request cannot be met.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Optional, Tuple

# Open MPI's session directory under TMPDIR: the default can exceed the socket
# path length on macOS. Other MPIs ignore these variables.
os.environ.setdefault("TMPDIR", "/tmp")
os.environ.setdefault("OMPI_MCA_prte_tmpdir_base", os.environ["TMPDIR"])
os.environ.setdefault("PRTE_MCA_prte_tmpdir_base", os.environ["TMPDIR"])


def mpi_exec_and_args() -> Tuple[str, str]:
    """Launcher and extra arguments: ``MPIEXEC`` and ``MPI_LAUNCH_EXTRA_ARGS``.

    Without ``MPIEXEC`` the Intel MPI launcher next to ``mpiifort`` is used when
    present, else ``mpiexec`` from ``PATH``. The extra arguments pass through
    unchanged.
    """
    mpiexec = os.environ.get("MPIEXEC")
    if not mpiexec:
        mpiifort = shutil.which("mpiifort")
        mpiexec = str(Path(mpiifort).parent / "mpiexec") if mpiifort else "mpiexec"
    return mpiexec, os.environ.get("MPI_LAUNCH_EXTRA_ARGS", "").strip()


def granted_cpus() -> Tuple[Optional[int], str]:
    """Number of CPUs granted to this process and where the number comes from."""
    env = os.environ
    if env.get("SLURM_JOB_ID"):
        if env.get("SLURM_NTASKS"):
            return int(env["SLURM_NTASKS"]) * int(env.get("SLURM_CPUS_PER_TASK", "1")), "Slurm allocation"
        if env.get("SLURM_CPUS_ON_NODE"):
            return int(env["SLURM_CPUS_ON_NODE"]) * int(env.get("SLURM_JOB_NUM_NODES", "1")), "Slurm allocation"
    if env.get("PBS_JOBID"):
        # PBS Pro: NCPUS is the ncpus of one chunk and PBS_NODEFILE lists one line
        # per mpiprocs (not per cpu), so the allocated cores are NCPUS x distinct hosts
        nodefile = env.get("PBS_NODEFILE")
        hosts = nodefile and Path(nodefile).is_file()
        if env.get("NCPUS") and hosts:
            nhosts = len({line.strip() for line in Path(nodefile).open() if line.strip()})
            return int(env["NCPUS"]) * nhosts, "PBS allocation"
        if env.get("NCPUS"):
            return int(env["NCPUS"]), "PBS allocation"
        if hosts:
            return sum(1 for line in Path(nodefile).open() if line.strip()), "PBS allocation"
    try:
        return len(os.sched_getaffinity(0)), "operating system"
    except AttributeError:  # macOS has no sched_getaffinity
        return os.cpu_count(), "operating system"


def preflight(nranks: int) -> None:
    """Raise ``RuntimeError`` when ``nranks`` exceeds the CPUs granted to this process."""
    n, source = granted_cpus()
    if n is None:
        return
    if nranks > n:
        hint = f"select=1:ncpus={nranks}:mpiprocs={nranks}" if source.startswith("PBS") else f"--ntasks={nranks}"
        raise RuntimeError(
            f"requested {nranks} MPI ranks but only {n} CPUs are available ({source}). "
            "uDALES runs one rank per CPU and never oversubscribes: use an allocation with "
            f"at least {nranks} CPUs (for example {hint}) or reduce nprocx*nprocy."
        )
