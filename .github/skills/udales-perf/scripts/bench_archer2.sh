#!/usr/bin/env bash
# Benchmark harness for exp 900 on ARCHER2 (Slurm), adapted from bench.sh (PBS).
#
# Usage: bench_archer2.sh <tag> [nreps] [release|debug]
#   <tag>   short label for this measurement series, e.g. "base", "B4" (no spaces)
#   [nreps] number of repetitions to submit (default 3), each its own Slurm job
#   [build] build type to benchmark (default release). Debug runs are recorded
#           under tag "<tag>-dbg", use fewer steps (101 instead of 1001 -- debug
#           is several x slower) and the standard QoS (short's 20 min walltime
#           cap doesn't fit debug's slower steps at a comfortable margin).
#
# Env overrides:
#   BENCH_BIN=<path>      benchmark this binary instead of $REPO/build/<btype>/u-dales
#   BENCH_COMMIT=<sha>    the commit that binary was built from (required with BENCH_BIN)
#   BENCH_REP_START=<n>   first rep number (default 1)
#   BENCH_DEPEND=<jobid>  make the first submitted job depend (afterany) on jobid;
#                         subsequent reps chain serially. Last jobid printed as
#                         "LASTJOB <id>" for wrapper scripts.
#   BENCH_ACCOUNT=<acct>  Slurm account (default n02-assure)
#
# ARCHER2 allocates whole nodes per job (node-hour billing), so single-node
# runs here do not suffer the CX3 co-tenancy noise -- but still record the
# node id per rep and treat cross-node comparisons with the usual care.
#
# The binary is SNAPSHOTTED into the run dir at submit time so a rebuild while
# jobs are queued cannot race them. Slurm -o/-e go to logs/. Per rep: run,
# harvest one CSV line into results.csv, copy log home, delete run dir.
# Guards: refuses a dirty tree or a stale binary.

set -e

TAG=${1:?usage: bench_archer2.sh <tag> [nreps] [release|debug]}
NREPS=${2:-3}
BTYPE=${3:-release}
REP_START=${BENCH_REP_START:-1}
ACCOUNT=${BENCH_ACCOUNT:-n02-assure}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
PROJECT_DIR="$(dirname "$REPO")"
EXPDIR="$PROJECT_DIR/experiments/900"
RUNROOT="$PROJECT_DIR/outputs/bench900"
RESULTS=$EXPDIR/results.csv
LOGDIR=$EXPDIR/logs

case $BTYPE in
    release) WALLTIME=00:15:00; QOS=short;    RUNTIME_OVERRIDE="" ;;
    debug)   WALLTIME=00:30:00; QOS=standard; RUNTIME_OVERRIDE="5."; TAG=${TAG}-dbg ;;
    *) echo "ERROR: build type must be release or debug" >&2; exit 1 ;;
esac

# ---- provenance ------------------------------------------------------------
if [ -n "$BENCH_BIN" ]; then
    BUILD=$BENCH_BIN
    COMMIT=${BENCH_COMMIT:?BENCH_BIN requires BENCH_COMMIT=<sha of that binary>}
    [ -x "$BUILD" ] || { echo "ERROR: $BUILD not found/executable." >&2; exit 1; }
else
    BUILD=$REPO/build/$BTYPE/u-dales
    cd $REPO
    if [ -n "$(git status --porcelain)" ]; then
        echo "ERROR: $REPO working tree is dirty; commit first (benchmarks must map to a commit)." >&2
        git status --short >&2
        exit 1
    fi
    COMMIT=$(git rev-parse --short=10 HEAD)
    [ -x "$BUILD" ] || { echo "ERROR: $BTYPE binary $BUILD not found." >&2; exit 1; }
    newest_src=$(find src CMakeLists.txt cmake -name '*.f90' -o -name 'CMakeLists.txt' -o -name '*.cmake' 2>/dev/null | xargs stat -c '%Y' | sort -n | tail -1)
    bin_mtime=$(stat -c '%Y' $BUILD)
    if [ "$bin_mtime" -lt "$newest_src" ]; then
        echo "ERROR: $BTYPE binary is older than the newest compile input; rebuild first." >&2
        exit 1
    fi
fi

# ---- results file -----------------------------------------------------------
mkdir -p $LOGDIR $RUNROOT
if [ ! -f $RESULTS ]; then
    echo "tag,rep,jobid,commit,node,cpu_model,nsteps,total_cpu_time_s,s_per_step" > $RESULTS
fi

# ---- submit reps ------------------------------------------------------------
cd $LOGDIR   # sbatch cwd: any stray spool files land here, never in the repo
for rep in $(seq $REP_START $(( REP_START + NREPS - 1 ))); do
    rundir=$RUNROOT/${TAG}_rep${rep}
    jobfile=$LOGDIR/job.${TAG}_rep${rep}.slurm

    # snapshot inputs + binary NOW: later rebuilds/edits cannot affect this job
    mkdir -p $rundir
    cp $EXPDIR/namoptions.900 $EXPDIR/*.inp.900 $rundir/
    cp $EXPDIR/solid_*.txt $EXPDIR/fluid_boundary_*.txt $EXPDIR/facet_sections_*.txt $rundir/
    cp $BUILD $rundir/u-dales.bench
    if [ -n "$RUNTIME_OVERRIDE" ]; then
        sed -i "s/^runtime .*/runtime      = $RUNTIME_OVERRIDE/" $rundir/namoptions.900
    fi

    cat <<EOF > $jobfile
#!/bin/bash
#SBATCH --job-name=bench900_${TAG}_${rep}
#SBATCH --time=$WALLTIME
#SBATCH --nodes=1
#SBATCH --tasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --account=$ACCOUNT
#SBATCH --partition=standard
#SBATCH --qos=$QOS
#SBATCH -o $LOGDIR/slurm.${TAG}_rep${rep}.out
#SBATCH -e $LOGDIR/slurm.${TAG}_rep${rep}.err

module load epcc-job-env >/dev/null 2>&1 || true
export OMP_NUM_THREADS=1
set -e
cd $rundir

log=$rundir/output.900
{
  echo "# tag=${TAG} rep=${rep} commit=${COMMIT} build=${BTYPE}"
  echo "# node=\$(hostname)"
  echo "# cpu=\$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')"
  echo "# modules: \$(module -t list 2>&1 | tr '\n' ' ')"
} > \$log

srun --distribution=block:block --hint=nomultithread -n 64 $rundir/u-dales.bench $rundir/namoptions.900 >> \$log 2>&1

# ---- harvest ----
nsteps=\$(grep -c 'Time of Simulation' \$log || true)
cpu_s=\$(grep -i 'TOTAL CPU time' \$log | awk '{print \$NF}')
node=\$(hostname)
cpu_model=\$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //;s/,/;/g')
if [ -n "\$cpu_s" ] && [ "\$nsteps" -gt 0 ]; then
    sps=\$(awk -v c=\$cpu_s -v n=\$nsteps 'BEGIN{printf "%.5f", c/n}')
else
    sps=FAILED
fi
echo "${TAG},${rep},\${SLURM_JOB_ID},${COMMIT},\$node,\$cpu_model,\$nsteps,\$cpu_s,\$sps" >> $RESULTS

cp \$log $LOGDIR/output.${TAG}_rep${rep}.log
rm -rf $rundir
EOF

    if [ -n "$BENCH_DEPEND" ]; then
        jobid=$(sbatch --parsable --dependency=afterany:$BENCH_DEPEND $jobfile)
    else
        jobid=$(sbatch --parsable $jobfile)
    fi
    echo "submitted $TAG rep $rep ($BTYPE, $(basename $BUILD)): $jobid (commit $COMMIT)"
    BENCH_DEPEND=$jobid   # chain reps serially when dependencies are in use
done
[ -n "$jobid" ] && echo "LASTJOB $jobid"
