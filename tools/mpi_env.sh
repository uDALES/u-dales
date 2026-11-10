# Shared MPI launch settings for the uDALES test runners and tools/local_execute.sh.
# Source it (bash) after the script's own defaults; it never oversubscribes.
#
#   MPIEXEC               launcher; default: Intel MPI's mpiexec next to mpiifort if
#                         present, else "mpiexec" from PATH. Override from the environment.
#   MPI_LAUNCH_EXTRA_ARGS extra launcher arguments, passed through unchanged (may be empty).
#   mpi_granted_cpus      prints the number of CPUs granted to this process and where the
#                         number comes from: the Slurm or PBS allocation when running in a
#                         batch job, otherwise the CPUs the operating system makes available.
#   mpi_preflight NRANKS  exits with an error when NRANKS exceeds the granted CPUs. uDALES
#                         runs one rank per CPU: oversubscription is never requested, so a
#                         short allocation fails here with both numbers instead of running
#                         several ranks on one core.
#
# On nodes with hyperthreading request physical cores from the scheduler (Slurm:
# --hint=nomultithread): the launcher binds one rank per core and refuses an allocation of
# hardware threads that puts two ranks on one core. Outside a scheduler Open MPI counts
# one slot per core; on a hyperthreaded laptop with fewer cores than ranks set
# OMPI_MCA_hwloc_base_use_hwthreads_as_cpus=true to run one rank per hardware thread.
#
# Open MPI's session directory is placed under TMPDIR (default /tmp): the default location
# can exceed the socket path length on macOS. Other MPIs ignore these variables.

if [ -z "${MPIEXEC:-}" ] && command -v mpiifort >/dev/null 2>&1; then
    MPIEXEC="$(dirname "$(command -v mpiifort)")/mpiexec"
else
    MPIEXEC="${MPIEXEC:-mpiexec}"
fi
MPI_LAUNCH_EXTRA_ARGS="${MPI_LAUNCH_EXTRA_ARGS:-}"
export MPIEXEC MPI_LAUNCH_EXTRA_ARGS

export TMPDIR="${TMPDIR:-/tmp}"
export OMPI_MCA_prte_tmpdir_base="${OMPI_MCA_prte_tmpdir_base:-$TMPDIR}"
export PRTE_MCA_prte_tmpdir_base="${PRTE_MCA_prte_tmpdir_base:-$TMPDIR}"

mpi_granted_cpus() {
    local n="" src=""
    if [ -n "${SLURM_JOB_ID:-}" ]; then
        if [ -n "${SLURM_NTASKS:-}" ]; then
            n=$(( SLURM_NTASKS * ${SLURM_CPUS_PER_TASK:-1} ))
        elif [ -n "${SLURM_CPUS_ON_NODE:-}" ]; then
            n=$(( SLURM_CPUS_ON_NODE * ${SLURM_JOB_NUM_NODES:-1} ))
        fi
        src="Slurm allocation"
    elif [ -n "${PBS_JOBID:-}" ]; then
        # PBS Pro: NCPUS is the ncpus of one chunk and PBS_NODEFILE lists one line per
        # mpiprocs (not per cpu), so the allocated cores are NCPUS x distinct hosts
        if [ -n "${NCPUS:-}" ] && [ -n "${PBS_NODEFILE:-}" ] && [ -r "$PBS_NODEFILE" ]; then
            n=$(( NCPUS * $(sort -u "$PBS_NODEFILE" | wc -l) ))
        elif [ -n "${NCPUS:-}" ]; then
            n=$NCPUS
        elif [ -n "${PBS_NODEFILE:-}" ] && [ -r "$PBS_NODEFILE" ]; then
            n=$(wc -l < "$PBS_NODEFILE")
        fi
        src="PBS allocation"
    fi
    if [ -z "$n" ]; then
        # GNU nproc honours OMP_NUM_THREADS / OMP_THREAD_LIMIT; ask for the CPUs in the affinity mask
        if command -v nproc >/dev/null 2>&1; then n=$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc)
        elif n=$(getconf _NPROCESSORS_ONLN 2>/dev/null) && [ -n "$n" ]; then :
        elif command -v sysctl >/dev/null 2>&1; then n=$(sysctl -n hw.logicalcpu 2>/dev/null)
        fi
        src="operating system"
    fi
    [ -n "$n" ] || return 1
    echo "$n $src"
}

mpi_preflight() {
    local nranks=$1 n src
    read -r n src < <(mpi_granted_cpus) || { echo "WARNING: cannot determine the number of available CPUs; launching $nranks ranks unchecked" >&2; return 0; }
    if [ "$nranks" -gt "$n" ]; then
        local hint="--ntasks=$nranks"
        case "$src" in PBS*) hint="select=1:ncpus=$nranks:mpiprocs=$nranks";; esac
        echo "ERROR: requested $nranks MPI ranks but only $n CPUs are available ($src)." >&2
        echo "       uDALES runs one rank per CPU and never oversubscribes: use an allocation with at least" >&2
        echo "       $nranks CPUs (for example $hint) or reduce nprocx*nprocy." >&2
        return 1
    fi
    return 0
}
