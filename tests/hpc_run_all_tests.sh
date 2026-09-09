#!/usr/bin/env bash

# uDALES (https://github.com/uDALES/u-dales).

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Copyright (C) 2016-2026 the uDALES Team.

# Run the test suite on the Imperial HX1 cluster with one command.
#
# Usage, from the u-dales directory:
#   tests/hpc_run_all_tests.sh              fast mode (default): about 30 min of
#                                           run time, apart from queue waits
#   tests/hpc_run_all_tests.sh full         everything
#   tests/hpc_run_all_tests.sh [mode] --dry-run
#                                           write every job script, submit nothing
#
# The orchestrator is itself a PBS job, so it keeps going after the terminal
# is closed. Three steps:
#
#   1. Build all four executables from clean build directories, one by one:
#      CPU Debug, CPU Release, GPU Debug, GPU Release. Both modes.
#   2. Submit the jobs that can run at the same time and wait for them all:
#        python      python-library and the GPU matrix config check   (CPU)
#        gpu-debug   gpu-smoke, gpu-mpi                       (2 GPUs; Debug dirs)
#        gpu-nightly gpu-nightly                              (1 GPU; Release dirs)
#        gpu-full    gpu-full, chained after gpu-nightly      (4 GPUs)   full only
#        fixtures    the GPU fixture matrix on the CPU Debug binary,
#                    on a node of its own (place=excl)               full only
#   3. Submit, one after the other, the jobs that must have the tree to
#      themselves, because the supported regression checks out master and
#      HEAD in the working tree itself:
#        supported   supported, Debug   (master vs HEAD)
#        all         all, Release       (master vs HEAD)              full only
#
# Why the split: gpu-full needs four free A100s, which on the a100 queue can
# take hours to appear (88 minutes today) for five minutes of work; the CPU
# fixture matrix takes 40-90 minutes on a Debug binary depending on what else
# is on the node; and "all" repeats every supported suite. Fast keeps every
# build, the Python stream, both GPU build types with CPU/GPU parity and the
# restart round trip, and all supported solver suites.
#
# Everything is written under logdir/run_all_tests_<timestamp>/ in the
# repository, the same place tools/ud_compare_multiple_inputs.py logs to:
#   report.log            the final report: mode, every selection, every log
#   orchestrator.log      what the orchestrator itself printed
#   build_<target>.log    one per executable, all attempts
#   <job>.pbs             each generated job script
#   <job>.out             each job's PBS output, with its resource usage
#   <job>_<selection>.log each selection's full output
#   <job>_summary.txt     each job's START/END lines and summary
#
# The environment the jobs use is the one documented in tests/README.md,
# "On the Imperial HX1 cluster".

set -u

# ---------------------------------------------------------------------------
# Cluster environment and resources.
# ---------------------------------------------------------------------------
NVHPC_ROOT=/gpfs/easybuild/prod/software/NVHPC/23.7-CUDA-12.2.0
IMPI_BIN=/gpfs/easybuild/prod/software/impi/2021.9.0-intel-compilers-2023.1.0/mpi/2021.9.0/bin

# The a100 queue is one node per job and 12 GPUs per user.
RES_ORCH="select=1:ncpus=8:mpiprocs=8:mem=64gb"
WALL_ORCH="24:00:00"
RES_CPU="select=1:ncpus=8:mpiprocs=8:mem=64gb"
RES_GPU1="select=1:ncpus=8:mpiprocs=8:mem=128gb:ngpus=1:gpu_type=A100"
RES_GPU2="select=1:ncpus=8:mpiprocs=8:mem=128gb:ngpus=2:gpu_type=A100"
RES_GPU4="select=1:ncpus=8:mpiprocs=8:mem=128gb:ngpus=4:gpu_type=A100"

# The block every test job starts with. Single-quoted heredoc: nothing here is
# expanded when the job script is written, only when it runs on the node.
job_environment() {
    cat <<'EOF'
cd "$UDALES_ROOT"
# Modules for the Intel CPU build and its runtime, plus CMake for the harness.
module load intel/2023a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a FFTW/3.3.10-intel-compilers-2023.1.0 CMake/3.26.3-GCCcore-12.3.0
# Solver integration suites load their runtime themselves from this variable.
export UDALES_RUNTIME_MODULES="intel/2023a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a FFTW/3.3.10-intel-compilers-2023.1.0"
# The supported regression builds master and HEAD with plain cmake.
export UDALES_CMAKE_ARGS="-DCMAKE_Fortran_COMPILER=$(command -v mpiifort) -DNETCDF_DIR=$EBROOTNETCDF -DNETCDF_FORTRAN_DIR=$EBROOTNETCDFMINFORTRAN"
# The experimental regressions build older revisions with their own script;
# "common" is the one target every revision has.
export UDALES_BUILD_SYSTEM=common
export MPIEXEC="$(command -v mpiexec)"
export UDALES_CPU_MPIEXEC="$MPIEXEC"
# GPU parity: which targets and wrappers build_test_binaries.sh should use,
# and which launcher runs each binary. NVHPC's mpiexec finds orted via PATH.
export UDALES_CPU_SYSTEM=hx1
export UDALES_GPU_SYSTEM=gpuhx1
export UDALES_CPU_FORTRAN_COMPILER=$IMPI_BIN/mpiifort
export UDALES_GPU_FORTRAN_COMPILER=$NVHPC_ROOT/Linux_x86_64/23.7/comm_libs/mpi/bin/mpif90
export UDALES_GPU_MPIEXEC=$NVHPC_ROOT/Linux_x86_64/23.7/comm_libs/mpi/bin/mpiexec
export PATH="$PATH:$NVHPC_ROOT/Linux_x86_64/23.7/comm_libs/mpi/bin"
# Parity compares outputs to 1e-6; keep host threads out of the reductions.
export OMP_NUM_THREADS=1
# HX1 injects Lmod into every child bash through BASH_ENV, which defeats the
# module stub two contract tests rely on. The exported function survives.
unset BASH_ENV
source tools/python/.venv/bin/activate
export MPLCONFIGDIR="${TMPDIR:-/tmp}/u-dales-matplotlib-$USER"; mkdir -p "$MPLCONFIGDIR"

run_step() {
    name=$1; shift
    echo "[$(date +%H:%M:%S)] START $name" | tee -a "$SUMMARY"
    "$@" > "$RUN_DIR/${JOB_NAME}_$name.log" 2>&1
    rc=$?
    echo "[$(date +%H:%M:%S)] END   $name rc=$rc" | tee -a "$SUMMARY"
    grep -E "^overall:|^- " "$RUN_DIR/${JOB_NAME}_$name.log" | sed 's/^/    /' >> "$SUMMARY"
}
echo "host=$(hostname) HEAD=$(git rev-parse --short HEAD) start=$(date)" | tee "$SUMMARY"
EOF
}

# write_job <name> <resources> <walltime> <step>...
# Each <step> is "label|command". Produces $RUN_DIR/<name>.pbs.
write_job() {
    local name=$1 resources=$2 walltime=$3; shift 3
    local script="$RUN_DIR/$name.pbs"
    {
        echo "#!/bin/bash"
        echo "#PBS -N udt-$name"
        echo "#PBS -l $resources"
        echo "#PBS -l walltime=$walltime"
        [ -n "${EXTRA_PBS:-}" ] && echo "$EXTRA_PBS"
        echo "#PBS -o $RUN_DIR/$name.out"
        echo "#PBS -j oe"
        echo "set -u"
        echo "UDALES_ROOT=$UDALES_ROOT"
        echo "RUN_DIR=$RUN_DIR"
        echo "JOB_NAME=$name"
        echo "SUMMARY=$RUN_DIR/${name}_summary.txt"
        echo "NVHPC_ROOT=$NVHPC_ROOT"
        echo "IMPI_BIN=$IMPI_BIN"
        job_environment
        local step
        for step in "$@"; do
            echo "run_step ${step%%|*} ${step#*|}"
        done
        echo 'echo "JOB DONE $(date)" | tee -a "$SUMMARY"'
    } > "$script"
}

# ---------------------------------------------------------------------------
# The jobs of each mode.
#   CONCURRENT: submitted together in step 2. SEQUENTIAL: one at a time in
#   step 3. CHAIN_gpu_full: the job gpu-full waits for (PBS afterany): it
#   shares the Release build directories with gpu-nightly, and this way the
#   4-GPU wait never holds anything else.
# ---------------------------------------------------------------------------
BRANCHES="--branch-a master --branch-b HEAD"
define_jobs() {
    write_job python      "$RES_CPU"  "02:00:00" \
        "python-library|python tests/run_tests.py python-library" \
        "gpu-validate-config|python tests/integration/gpu/run_gpu_tests.py full --validate-config"
    write_job gpu-debug   "$RES_GPU2" "06:00:00" \
        "gpu-smoke|python tests/run_tests.py gpu-smoke" \
        "gpu-mpi|python tests/run_tests.py gpu-mpi"
    write_job gpu-nightly "$RES_GPU1" "06:00:00" \
        "gpu-nightly|python tests/run_tests.py gpu-nightly"
    write_job supported   "$RES_CPU"  "03:00:00" \
        "supported-debug|python tests/run_tests.py supported $BRANCHES --build-type Debug"
    CONCURRENT="python gpu-debug gpu-nightly"
    SEQUENTIAL="supported"
    SKIPPED=""
    CHAIN_gpu_full=""
    if [ "$MODE" = full ]; then
        write_job gpu-full  "$RES_GPU4" "06:00:00" \
            "gpu-full|python tests/run_tests.py gpu-full"
        # Reads only the CPU Debug binary. The GPU jobs re-run cmake in that
        # build directory as their first suite, but with nothing changed cmake
        # relinks nothing, so the binary this job runs is never rewritten.
        #
        # place=excl: the 128^3 ibm-reconstruction case is memory-bandwidth
        # bound at -O0 with all checks on, and a heavy neighbour on the node
        # slowed it 3.5x (45 s -> 157 s per step) and past its timeout. A node
        # of its own makes its run time predictable; the price is waiting for
        # an idle node in the small/medium pool.
        EXTRA_PBS="#PBS -l place=excl" write_job fixtures "$RES_CPU" "03:00:00" \
            "gpu-fixtures-cpu|python tests/integration/gpu/run_gpu_tests.py full --cpu-only --cpu-executable build/cpu/debug/u-dales"
        write_job all       "$RES_CPU"  "04:00:00" \
            "all-release|python tests/run_tests.py all $BRANCHES --build-type Release"
        CONCURRENT="python gpu-debug gpu-nightly gpu-full fixtures"
        SEQUENTIAL="supported all"
        CHAIN_gpu_full="gpu-nightly"
    else
        SKIPPED="gpu-full (4 GPUs), the CPU fixture matrix (40-90 min), all Release (repeats supported plus the experimental stream)"
    fi
}

# submit <name> [after-name]: qsub the generated script, record the id, and
# print it. With after-name, the job waits for that job to end first.
submit() {
    local name=$1 after=${2:-} id dep=""
    if [ -n "$after" ] && [ -s "$RUN_DIR/$after.jobid" ]; then
        dep="-W depend=afterany:$(cat "$RUN_DIR/$after.jobid")"
    fi
    id=$(qsub $dep "$RUN_DIR/$name.pbs") || { echo "qsub failed for $name" >> "$REPORT"; return 1; }
    echo "$id" > "$RUN_DIR/$name.jobid"
    echo "[$(date +%H:%M:%S)] submitted $name as $id${after:+ (after $after)}" >> "$ORCHLOG"
    echo "$id"
}

# wait_for <id>...: poll until none of the jobs is known to qstat any more.
wait_for() {
    local remaining id
    while true; do
        remaining=""
        for id in "$@"; do
            qstat "$id" >/dev/null 2>&1 && remaining="$remaining $id"
        done
        [ -z "$remaining" ] && return 0
        sleep 60
    done
}

# queue_wait <id>: "Xm YYs" from PBS's own qtime and stime, or "?".
queue_wait() {
    local info qt st w
    info=$(qstat -x -f "$1" 2>/dev/null) || { echo "?"; return; }
    qt=$(printf '%s\n' "$info" | awk -F' = ' '/^[[:space:]]*qtime = /{print $2}')
    st=$(printf '%s\n' "$info" | awk -F' = ' '/^[[:space:]]*stime = /{print $2}')
    if [ -n "$qt" ] && [ -n "$st" ]; then
        w=$(( $(date -d "$st" +%s) - $(date -d "$qt" +%s) ))
        printf '%dm %02ds' $((w / 60)) $((w % 60))
    else
        echo "?"
    fi
}

# report_job <name>: one line per selection into the report.
report_job() {
    local name=$1 id summary out log label rc passes fails wall status
    id=$(cat "$RUN_DIR/$name.jobid" 2>/dev/null || echo "?")
    summary="$RUN_DIR/${name}_summary.txt"
    out="$RUN_DIR/$name.out"
    wall=$(grep -oE "Walltime usage: [0-9:]+" "$out" 2>/dev/null | sed 's/Walltime usage: //')
    echo "  $name (job $id, queued ${wall:+$(queue_wait "$id"), }run ${wall:-?})" >> "$REPORT"
    if [ ! -f "$summary" ]; then
        if [ -f "$out" ]; then
            echo "      no summary: the job started but died before writing one; see $out" >> "$REPORT"
        else
            echo "      never started (no PBS output file): cancelled, or still queued when the orchestrator gave up" >> "$REPORT"
        fi
        return
    fi
    while read -r label rc; do
        log="$RUN_DIR/${name}_$label.log"
        passes=$(grep -cE "^- .*: PASS$" "$log" 2>/dev/null || true)
        fails=$(grep -cE "^- .*: FAIL$" "$log" 2>/dev/null || true)
        if [ "$rc" = 0 ]; then status=PASS; else status=FAIL; fi
        printf "      %-22s %-4s  rc=%s  suites %s pass / %s fail   %s\n" "$label" "$status" "$rc" "$passes" "$fails" "$log" >> "$REPORT"
        if [ "$fails" -gt 0 ]; then
            grep -E "^- .*: FAIL$" "$log" | sed 's/^- /          failed: /' >> "$REPORT"
        fi
    done < <(grep -E "^\[.*\] END " "$summary" | sed -E 's/^\[.*\] END +([^ ]+) rc=([0-9]+).*/\1 \2/')
}

# ---------------------------------------------------------------------------
# Orchestrator: what the submitted job executes.
# ---------------------------------------------------------------------------
orchestrate() {
    cd "$UDALES_ROOT"
    ORCHLOG="$RUN_DIR/orchestrator.log"
    REPORT="$RUN_DIR/report.log"
    : > "$ORCHLOG"
    define_jobs
    {
        echo "uDALES test suite on HX1  -  $(date)"
        echo "mode: $MODE  (argument given: ${MODE_ARG:-none, so the default})"
        echo "repository: $UDALES_ROOT  HEAD: $(git rev-parse --short HEAD) ($(git symbolic-ref --short -q HEAD || echo detached))"
        echo "run directory: $RUN_DIR"
        echo "selections run: concurrently [$CONCURRENT], then one at a time [$SEQUENTIAL]"
        [ -n "$SKIPPED" ] && echo "skipped in fast mode: $SKIPPED"
        if [ -n "$(git status --porcelain)" ]; then
            echo "NOTE: the working tree has uncommitted changes, so the supported"
            echo "      regression will skip itself (it checks out master and HEAD in place)."
        fi
        echo
    } > "$REPORT"

    # ---- Step 1: fresh builds -------------------------------------------
    # A clean build clones FindFFTW and 2decomp-fft from GitHub and then
    # fetches again for the update step, and any of those can fail on a
    # transient network error ("Empty reply from server"). Retry from clean
    # before calling a build failed; every attempt is kept in the log.
    echo "Step 1: builds" >> "$REPORT"
    local ok=1 target system type dir log t0 attempt
    for target in cpu:debug cpu:release gpu:debug gpu:release; do
        type=${target#*:}
        case ${target%%:*} in cpu) system=hx1 ;; gpu) system=gpuhx1 ;; esac
        dir="$UDALES_ROOT/build/${target%%:*}/$type"
        log="$RUN_DIR/build_${target%%:*}_$type.log"
        : > "$log"
        t0=$(date +%s)
        for attempt in 1 2 3; do
            echo "[$(date +%H:%M:%S)] building $target ($system $type) from clean, attempt $attempt" | tee -a "$ORCHLOG"
            echo "=================== attempt $attempt  $(date) ===================" >> "$log"
            rm -rf "$dir"
            ( ./tools/build_executable.sh "$system" "$type" ) >> "$log" 2>&1
            [ -x "$dir/u-dales" ] && break
            [ "$attempt" -lt 3 ] && sleep 60
        done
        if [ -x "$dir/u-dales" ]; then
            printf "  %-12s OK    %4ds  attempt %s  %s\n" "$target" $(( $(date +%s) - t0 )) "$attempt" "$log" >> "$REPORT"
        else
            printf "  %-12s FAILED after 3 attempts  %s\n" "$target" "$log" >> "$REPORT"
            ok=0
        fi
    done
    echo >> "$REPORT"
    if [ "$ok" -eq 0 ]; then
        echo "A build failed; no test jobs submitted." | tee -a "$REPORT" "$ORCHLOG"
        cat "$REPORT"
        return 1
    fi

    if [ "$DRY_RUN" = 1 ]; then
        echo "dry run: job scripts written to $RUN_DIR, nothing submitted" | tee -a "$REPORT"
        return 0
    fi

    # ---- Step 2: the jobs that can run at the same time -----------------
    local ids="" id name after
    echo "Step 2: concurrent jobs" >> "$REPORT"
    for name in $CONCURRENT; do
        after=""; [ "$name" = gpu-full ] && after="$CHAIN_gpu_full"
        id=$(submit "$name" "$after") || return 1
        ids="$ids $id"
    done
    echo "[$(date +%H:%M:%S)] waiting for:$ids" | tee -a "$ORCHLOG"
    wait_for $ids
    for name in $CONCURRENT; do report_job "$name"; done
    echo >> "$REPORT"

    # ---- Step 3: the jobs that must have the tree to themselves ---------
    echo "Step 3: sequential jobs" >> "$REPORT"
    for name in $SEQUENTIAL; do
        id=$(submit "$name") || return 1
        echo "[$(date +%H:%M:%S)] waiting for $id" | tee -a "$ORCHLOG"
        wait_for "$id"
        report_job "$name"
    done
    echo >> "$REPORT"

    # ---- Final report ----------------------------------------------------
    # Count selection status lines only, not the "failed:" detail lines.
    local failed
    failed=$(grep -cE "^      [^ ]+ +FAIL  rc=" "$REPORT" || true)
    if [ "$failed" -eq 0 ]; then
        echo "RESULT: every selection passed ($MODE mode)." >> "$REPORT"
    else
        echo "RESULT: $failed selection(s) failed ($MODE mode); see the logs listed above." >> "$REPORT"
    fi
    echo "finished: $(date)" >> "$REPORT"
    echo >> "$REPORT"
    echo "logs:" >> "$REPORT"
    ls -1 "$RUN_DIR" | sed "s#^#  $RUN_DIR/#" >> "$REPORT"
    cat "$REPORT"
}

# ---------------------------------------------------------------------------
# Entry point.
# ---------------------------------------------------------------------------
DRY_RUN=0
ACTION=submit
MODE=fast
MODE_ARG=""
for arg in "$@"; do
    case "$arg" in
        fast|full) MODE=$arg; MODE_ARG=$arg ;;
        --dry-run) DRY_RUN=1 ;;
        --orchestrate) ACTION=orchestrate ;;
        *) echo "Usage: tests/hpc_run_all_tests.sh [fast|full] [--dry-run]"; exit 2 ;;
    esac
done

if [ "$ACTION" = orchestrate ]; then
    # Inside the orchestrator job: UDALES_ROOT, RUN_DIR, MODE and MODE_ARG
    # arrived through qsub -v.
    orchestrate
    exit $?
fi

if [ ! -d src ] || [ ! -f tools/build_executable.sh ]; then
    echo "Run this from the u-dales directory."
    exit 1
fi
case "$(hostname -s)" in
    hx1*) ;;
    *) echo "This script is written for HX1; hostname is $(hostname -s)."; exit 1 ;;
esac
if [ ! -f tools/python/.venv/bin/activate ]; then
    echo "tools/python/.venv is missing; set it up with: bash tools/python/setup_venv.sh icl"
    exit 1
fi

UDALES_ROOT="$(pwd -P)"
RUN_DIR="$UDALES_ROOT/logdir/run_all_tests_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RUN_DIR"

if [ "$DRY_RUN" = 1 ]; then
    # Write only the job scripts; skip the builds.
    define_jobs
    echo "dry run ($MODE mode): job scripts written to $RUN_DIR"
    echo "  concurrent: $CONCURRENT"
    [ -n "$CHAIN_gpu_full" ] && echo "  gpu-full would be chained after: $CHAIN_gpu_full"
    echo "  sequential: $SEQUENTIAL"
    [ -n "$SKIPPED" ] && echo "  skipped:    $SKIPPED"
    exit 0
fi

# The orchestrator job: builds, then submits and waits for the test jobs.
ORCH="$RUN_DIR/orchestrator.pbs"
cat > "$ORCH" <<EOF
#!/bin/bash
#PBS -N udt-orchestrator
#PBS -l $RES_ORCH
#PBS -l walltime=$WALL_ORCH
#PBS -o $RUN_DIR/orchestrator.out
#PBS -j oe
#PBS -v UDALES_ROOT=$UDALES_ROOT,RUN_DIR=$RUN_DIR,MODE=$MODE,MODE_ARG=$MODE_ARG
cd "\$UDALES_ROOT"
bash tests/hpc_run_all_tests.sh --orchestrate
EOF
id=$(qsub "$ORCH") || { echo "qsub failed"; exit 1; }
echo "$id" > "$RUN_DIR/orchestrator.jobid"
echo "Submitted the test-suite orchestrator as $id ($MODE mode)."
echo "It builds the four executables, then submits and waits for the test jobs."
echo "Everything is written under:"
echo "  $RUN_DIR"
echo "Final report, once it finishes:  $RUN_DIR/report.log"
echo "Progress in the meantime:        $RUN_DIR/orchestrator.log  and  qstat -u \$USER"
if [ -n "$(git status --porcelain)" ]; then
    echo "Note: the working tree has uncommitted changes, so the supported regression"
    echo "will skip itself. Commit first if you want that harness to run."
fi
