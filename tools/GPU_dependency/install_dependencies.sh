#!/usr/bin/env bash
#
# Build every library the GPU build of uDALES needs, with one command.
#
# Usage:
#   tools/GPU_dependency/install_dependencies.sh [hx1|common] [options]
#
#   hx1      (default) submit the build as a PBS job on Imperial HX1
#   common   run the build in the background on this machine (nohup), so it
#            survives closing the terminal. Loads the NVIDIA HPC SDK with
#               module use /opt/nvidia/hpc_sdk/modulefiles; module load nvhpc/24.11
#            and takes cmake and the zlib headers from the system. Another
#            SDK location or version: edit NVHPC_MODULEFILES / NVHPC_MODULE
#            in tools/GPU_dependency/modules.sh (or export them). If the module
#            is missing the script stops before starting and says what exists
#            and what to change.
#
# Options:
#   --stack 1|2|3        MPI stack to build, see README.md (default 3)
#   --jobs N             cores for the compilers and, on hx1, the job's ncpus
#                        (default 8)
#   --walltime HH:MM:SS  hx1 job walltime (default 02:00:00 for stack 3, which
#                        builds OpenMPI as well, 01:00:00 for stacks 1 and 2;
#                        a stack 3 build took 46 min on 8 cores, so a short
#                        request that still has room keeps the queue wait short)
#   --mem SIZE           hx1 job memory (default 32gb)
#   --on-gpu-queue       hx1: request one A100 as well, so the job lands on the
#                        GPU queue (useful when the CPU queues are reserved)
#   --select SPEC        hx1: full PBS select string, replaces the one built
#                        from --jobs / --mem / --on-gpu-queue
#   --force              rebuild components that are already installed in
#                        build/ (by default they are skipped, so a failed run
#                        is resumed by running the same command again)
#   --dry-run            show what would be submitted or started, do nothing
#
# Everything is built under <stack directory>/build/ (git-ignored): the stage
# scripts are copied there and run in place, so the sources, the installed
# libraries and the logs all live in that one directory. Stack 3 builds its
# own OpenMPI first (build/openmpi/openmpi) and the libraries link against it.
#
# The log is build/logs/install_<time>.log; the script returns as soon as the
# job is submitted (or the background process started) and prints the command
# to follow it. The last lines of the log give the variables for the gpuhx1
# block of tools/build_executable.sh.
#
# All module names live in tools/GPU_dependency/modules.sh, selected by the
# hx1|common argument (exported to the stage scripts as UDALES_DEP_SYSTEM).
#
# Internal: --run executes the stages in the current shell (what the job and
# the background process call); --log PATH names its log file.
set -euo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

SYSTEM=hx1
STACK=3
JOBS=8
WALLTIME=""
MEM=32gb
SELECT=""
GPU_QUEUE=0
FORCE=0
DRY_RUN=0
RUN=0
LOG=""

usage() { sed -n '2,/^set -euo/p' "${BASH_SOURCE[0]}" | sed '$d' | cut -c3-; }

while [ $# -gt 0 ]; do
    case "$1" in
        hx1|common)     SYSTEM=$1 ;;
        --stack)        STACK=$2; shift ;;
        --stack=*)      STACK=${1#*=} ;;
        --jobs)         JOBS=$2; shift ;;
        --jobs=*)       JOBS=${1#*=} ;;
        --walltime)     WALLTIME=$2; shift ;;
        --walltime=*)   WALLTIME=${1#*=} ;;
        --mem)          MEM=$2; shift ;;
        --mem=*)        MEM=${1#*=} ;;
        --select)       SELECT=$2; shift ;;
        --select=*)     SELECT=${1#*=} ;;
        --on-gpu-queue) GPU_QUEUE=1 ;;
        --force)        FORCE=1 ;;
        --dry-run)      DRY_RUN=1 ;;
        --run)          RUN=1 ;;
        --log)          LOG=$2; shift ;;
        -h|--help)      usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
    shift
done

case "$STACK" in
    1) STACK_DIR=$here/stack1-NVHPC-OpenMPI-3.1.5 ;;
    2) STACK_DIR=$here/stack2-NVHPC-OpenMPI-4.1.5 ;;
    3) STACK_DIR=$here/stack3-OpenMPI-4.1.6-tm ;;
    *) echo "--stack must be 1, 2 or 3 (got '$STACK')" >&2; exit 2 ;;
esac
[[ "$JOBS" =~ ^[0-9]+$ && "$JOBS" -ge 1 ]] || { echo "--jobs must be a positive integer (got '$JOBS')" >&2; exit 2; }
# Measured on 8 cores: stack 3 takes 46 min (31 of them OpenMPI), the
# libraries alone 15 min. A short walltime is placed sooner; this still leaves
# more than twice the measured time.
if [ -z "$WALLTIME" ]; then
    if [ "$STACK" = 3 ]; then WALLTIME=02:00:00; else WALLTIME=01:00:00; fi
fi

BUILD=$STACK_DIR/build
LOGDIR=$BUILD/logs
STAMP=$(date +%Y%m%d_%H%M%S)
[ -n "$LOG" ] || LOG=$LOGDIR/install_$STAMP.log

LIBS="hdf5 netcdf_c netcdf_fortran fftw"

# Copy the stage scripts into build/. They locate everything relative to their
# own directory, so the sources and installs end up in build/ as well.
stage_scripts() {
    mkdir -p "$LOGDIR"
    cp -p "$STACK_DIR"/*.sh "$BUILD/"
    if [ "$STACK" = 3 ]; then
        mkdir -p "$BUILD/openmpi/logs"
        cp -p "$STACK_DIR"/openmpi/*.sh "$BUILD/openmpi/"
    fi
}

# Lmod's "module" is a shell function; a batch job or a nohup shell may not
# have it yet.
ensure_module_command() {
    if ! command -v module >/dev/null 2>&1; then
        local f
        for f in /etc/profile.d/modules.sh /etc/profile.d/lmod.sh /usr/share/lmod/lmod/init/bash; do
            if [ -r "$f" ]; then
                # shellcheck disable=SC1090
                source "$f" && break
            fi
        done
    fi
    command -v module >/dev/null 2>&1 || {
        echo "Error: no 'module' command. The stage scripts load the NVIDIA HPC SDK with it" >&2
        echo "(module use /opt/nvidia/hpc_sdk/modulefiles; module load nvhpc/24.11, see" >&2
        echo "$here/modules.sh). Install Environment Modules or Lmod, or make sure" >&2
        echo "/etc/profile.d/modules.sh exists, and run again." >&2
        return 1
    }
}

# Load the modules once here, before anything is submitted or started, so a
# missing module is reported in the terminal with the fix, not in a log.
check_modules() {
    local extra=""
    [ "$STACK" = 3 ] && extra=openmpi
    ensure_module_command || return 1
    ( export UDALES_DEP_SYSTEM=$SYSTEM; source "$here/modules.sh" $extra >/dev/null ) || {
        echo "Fix the above and run install_dependencies.sh again." >&2
        return 1
    }
}

# Download stages: no compilers needed, they skip what is already present.
download_all() {
    local lib
    if [ "$STACK" = 3 ]; then
        bash "$BUILD/openmpi/s0_download.sh" || return 1
    fi
    for lib in $LIBS; do
        bash "$BUILD/${lib}_s0_download.sh" || return 1
    done
}

# ---------------------------------------------------------------------------
# --run: the build itself, in this shell, appending everything to $LOG.
# ---------------------------------------------------------------------------
step() {
    local name=$1; shift
    local t0=$SECONDS
    echo
    echo "================ $name  ($(date '+%F %T')) ================"
    if "$@"; then
        echo "---------------- $name done in $(( SECONDS - t0 )) s"
    else
        echo "---------------- $name FAILED after $(( SECONDS - t0 )) s"
        echo
        echo "INSTALL FAILED at: $name"
        echo "Log: $LOG"
        exit 1
    fi
}

# installed NAME MARKER: true (and say so) when MARKER exists and --force is
# not given, so a rerun after a failure continues where it stopped.
installed() {
    if [ "$FORCE" = 0 ] && [ -e "$2" ]; then
        echo
        echo "================ $1: already installed ($2), skipping; --force rebuilds"
        return 0
    fi
    return 1
}

run_all() {
    stage_scripts
    echo "uDALES GPU dependencies, stack $STACK  ($(date '+%F %T'))"
    echo "host: $(hostname -s)   cores: $JOBS   build directory: $BUILD"
    echo "log: $LOG"
    ensure_module_command
    export BUILD_JOBS=$JOBS
    export UDALES_DEP_SYSTEM=$SYSTEM

    if [ "$STACK" = 3 ]; then
        # The libraries' s1 scripts take the MPI from OMPI_PREFIX.
        export OMPI_PREFIX=$BUILD/openmpi/openmpi
        if installed openmpi "$OMPI_PREFIX/bin/ompi_info"; then :; else
        step "openmpi download"            bash "$BUILD/openmpi/s0_download.sh"
        step "openmpi environment"         bash "$BUILD/openmpi/s1_set_env.sh"
        step "openmpi configure"           bash "$BUILD/openmpi/s2_configure.sh"
        step "openmpi build and install"   bash "$BUILD/openmpi/s3_build_install.sh"
        fi
    fi

    local lib s3 marker
    for lib in $LIBS; do
        case $lib in
            hdf5)           marker=$BUILD/hdf5-2.1.1/hdf5/bin/h5pcc ;;
            netcdf_c)       marker=$BUILD/netcdf-c-4.10.1/netcdfc/bin/nc-config ;;
            netcdf_fortran) marker=$BUILD/netcdf-fortran-4.6.4/netcdff/bin/nf-config ;;
            fftw)           marker=$BUILD/fftw-3.3.11/fftw3/lib/libfftw3_mpi.so ;;
        esac
        installed "$lib" "$marker" && continue
        s3=$(ls "$BUILD"/${lib}_s3_*.sh)
        step "$lib download"               bash "$BUILD/${lib}_s0_download.sh"
        step "$lib environment"            bash "$BUILD/${lib}_s1_set_env_for_install.sh"
        step "$lib configure"              bash "$BUILD/${lib}_s2_configure_install.sh"
        step "$lib build and install"      bash "$s3"
    done

    echo
    echo "INSTALL COMPLETE  ($(date '+%F %T'))"
    echo
    echo "Installed under $BUILD:"
    [ "$STACK" = 3 ] && echo "  OpenMPI         $OMPI_PREFIX"
    echo "  HDF5            $BUILD/hdf5-2.1.1/hdf5"
    echo "  netCDF-C        $BUILD/netcdf-c-4.10.1/netcdfc"
    echo "  netCDF-Fortran  $BUILD/netcdf-fortran-4.6.4/netcdff"
    echo "  FFTW            $BUILD/fftw-3.3.11/fftw3"
    echo
    if [ "$STACK" = 3 ]; then
        echo "Stack 3 in this directory is the default of tools/build_executable.sh (gpuhx1),"
        echo "tools/hpc_execute.sh (hx1:gpu) and tests/hpc_run_all_tests.sh: nothing to edit."
        echo "Build with:  tools/build_executable.sh gpuhx1 release"
    else
        local relbuild=tools/GPU_dependency/$(basename "$STACK_DIR")/build
        case $STACK in
            1) mpibin='$EBROOTNVHPC/Linux_x86_64/23.7/comm_libs/mpi/bin' ;;
            2) mpibin='$EBROOTNVHPC/Linux_x86_64/23.7/comm_libs/12.2/openmpi4/openmpi-4.1.5/bin' ;;
        esac
        echo "The uDALES scripts default to stack 3. To use stack $STACK, edit these three"
        echo "files (each keeps the stack $STACK lines as comments next to the stack 3 ones):"
        echo
        echo "1. tools/build_executable.sh, block 'gpuhx1': comment out the stack 3 lines and set"
        echo "     gpu_stack_dir=\"\${UDALES_GPU_STACK_DIR:-\$PWD/$relbuild}\""
        echo "     NVHPC_MPI_BIN=\"$mpibin\""
        echo "   (NETCDF_DIR, NETCDF_FORTRAN_DIR and FFTWDIR follow from gpu_stack_dir unchanged)"
        echo
        echo "2. tools/hpc_execute.sh, block 'hx1:gpu': use the stack $STACK job_modules, which puts"
        echo "     $mpibin"
        echo "   on PATH instead of the stack 3 OpenMPI."
        echo
        echo "3. tests/hpc_run_all_tests.sh: use the stack $STACK lines for"
        echo "     UDALES_GPU_FORTRAN_COMPILER, UDALES_GPU_MPIEXEC and PATH"
        echo "   (they take mpif90 and mpiexec from $mpibin)."
        echo
        echo "Stack $STACK runs on one node only; see tools/GPU_dependency/README.md."
    fi
}

if [ "$RUN" = 1 ]; then
    # A pipe rather than a process substitution, so that the last lines are in
    # the log before the script exits and the launcher sees them.
    mkdir -p "$(dirname "$LOG")"
    run_all 2>&1 | tee -a "$LOG"
    exit "${PIPESTATUS[0]}"
fi

# ---------------------------------------------------------------------------
# Launchers
# ---------------------------------------------------------------------------
FORCE_ARG=""; [ "$FORCE" = 1 ] && FORCE_ARG=" --force"
echo "uDALES GPU dependencies: stack $STACK on $SYSTEM, $JOBS cores"
echo "build directory: $BUILD"

if [ "$DRY_RUN" = 0 ]; then
    stage_scripts
fi
check_modules || exit 1

case "$SYSTEM" in
hx1)
    if [ -z "$SELECT" ]; then
        SELECT="select=1:ncpus=$JOBS:mem=$MEM"
        [ "$GPU_QUEUE" = 1 ] && SELECT="$SELECT:ngpus=1:gpu_type=A100"
    fi
    PBS_SCRIPT=$LOGDIR/install_$STAMP.pbs
    script=$(cat <<EOF
#!/bin/bash
#PBS -N gpudep-stack$STACK
#PBS -l $SELECT
#PBS -l walltime=$WALLTIME
#PBS -o $LOGDIR/install_$STAMP.pbs.out
#PBS -j oe
cd "$here"
exec bash "$here/install_dependencies.sh" --run $SYSTEM --stack $STACK --jobs $JOBS --log "$LOG"$FORCE_ARG
EOF
)
    if [ "$DRY_RUN" = 1 ]; then
        echo
        echo "Would download the sources on this node, then qsub this script:"
        echo "$script" | sed 's/^/    /'
        exit 0
    fi
    # Compute nodes may have no internet: fetch the tarballs here first (the
    # job runs the same stages again and finds them present).
    echo
    echo "Downloading the sources on this node..."
    download_all > "$LOGDIR/download_$STAMP.log" 2>&1 || {
        echo "Download failed, see $LOGDIR/download_$STAMP.log" >&2
        tail -5 "$LOGDIR/download_$STAMP.log" >&2
        exit 1
    }
    printf '%s\n' "$script" > "$PBS_SCRIPT"
    JOBID=$(qsub "$PBS_SCRIPT")
    echo "Submitted job $JOBID  ($SELECT, walltime $WALLTIME)"
    echo "PBS script: $PBS_SCRIPT"
    echo
    echo "The log appears once the job starts. Follow it with:"
    echo "    tail -f $LOG"
    echo "Job state:  qstat $JOBID"
    echo "It ends with INSTALL COMPLETE or INSTALL FAILED."
    ;;
common)
    cmd=(bash "$here/install_dependencies.sh" --run "$SYSTEM" --stack "$STACK" --jobs "$JOBS" --log "$LOG")
    [ "$FORCE" = 1 ] && cmd+=(--force)
    if [ "$DRY_RUN" = 1 ]; then
        echo
        echo "Would start in the background (nohup):"
        echo "    ${cmd[*]}"
        exit 0
    fi
    nohup "${cmd[@]}" > /dev/null 2>&1 &
    PID=$!
    echo "Started background process $PID (survives closing this terminal)"
    echo
    echo "Follow the log with:"
    echo "    tail -f $LOG"
    echo "It ends with INSTALL COMPLETE or INSTALL FAILED."
    ;;
esac
