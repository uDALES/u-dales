#!/usr/bin/env bash
# Benchmark harness for exp 900 (uDALES issue #330 performance work).
#
# Usage: bench.sh <tag> [nreps] [release|debug]
#   <tag>   short label for this measurement series, e.g. "base", "B0", "B1" (no spaces)
#   [nreps] number of repetitions to submit (default 3), each as its own PBS job
#   [build] build type to benchmark (default release). Debug runs are recorded
#           under tag "<tag>-dbg", use a shorter simulation (101 steps instead of
#           1001 -- debug is several x slower; s/step is the metric either way)
#           and a 30-min walltime.
#
# Env overrides:
#   BENCH_BIN=<path>     benchmark this binary instead of $REPO/build/<btype>/u-dales
#                        (e.g. an archived bin/u-dales.<commit>.<btype>); requires
#   BENCH_COMMIT=<sha>   the commit that binary was built from (skips HEAD/guards)
#   BENCH_REP_START=<n>  first rep number (default 1); use to resubmit a lost rep
#   BENCH_DEPEND=<jobid> make the (first) submitted job depend afterany on jobid,
#                        and chain subsequent reps serially; the last jobid is
#                        printed as "LASTJOB <id>" for wrapper scripts
#
# The binary is SNAPSHOTTED into the run dir at submit time, so later rebuilds
# of build/<btype> can never race a queued job. PBS -o/-e files go to logs/.
# Per rep: run, harvest one CSV line into results.csv, copy log home, delete
# ephemeral run dir. Guards: refuses a dirty tree or a stale binary.

set -e

TAG=${1:?usage: bench.sh <tag> [nreps] [release|debug]}
NREPS=${2:-3}
BTYPE=${3:-release}
REP_START=${BENCH_REP_START:-1}

REPO=~/udales/u-dales
EXPDIR=~/udales/experiments/900
RESULTS=$EXPDIR/results.csv
LOGDIR=$EXPDIR/logs

case $BTYPE in
    release) WALLTIME=00:15:00; RUNTIME_OVERRIDE="" ;;
    debug)   WALLTIME=00:30:00; RUNTIME_OVERRIDE="5."; TAG=${TAG}-dbg ;;
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
mkdir -p $LOGDIR
if [ ! -f $RESULTS ]; then
    echo "tag,rep,jobid,commit,node,cpu_model,nsteps,total_cpu_time_s,s_per_step" > $RESULTS
fi

# ---- submit reps ------------------------------------------------------------
cd $LOGDIR   # qsub cwd: any PBS spool files land here, never in the repo
for rep in $(seq $REP_START $(( REP_START + NREPS - 1 ))); do
    rundir=$EPHEMERAL/900/${TAG}_rep${rep}
    jobfile=$LOGDIR/job.${TAG}_rep${rep}

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
#PBS -l walltime=$WALLTIME
#PBS -l select=1:ncpus=64:mpiprocs=64:mem=120gb:cpu_type=rome
#PBS -o $LOGDIR/
#PBS -e $LOGDIR/
module load intel/2021a netCDF/4.8.0-iimpi-2021a netCDF-Fortran/4.5.3-iimpi-2021a FFTW/3.3.9-intel-2021a
# some CX3 nodes firewall TCP-to-self, killing hydra bootstrap; loopback avoids it
export I_MPI_HYDRA_IFACE=lo
set -e
cd $rundir

log=$rundir/output.900
{
  echo "# tag=${TAG} rep=${rep} commit=${COMMIT} build=${BTYPE}"
  echo "# node=\$(hostname)"
  echo "# cpu=\$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')"
  echo "# modules: \$(module -t list 2>&1 | tr '\n' ' ')"
} > \$log

mpiexec -n 64 $rundir/u-dales.bench $rundir/namoptions.900 >> \$log 2>&1

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
echo "${TAG},${rep},\${PBS_JOBID%%.*},${COMMIT},\$node,\$cpu_model,\$nsteps,\$cpu_s,\$sps" >> $RESULTS

cp \$log $LOGDIR/output.${TAG}_rep${rep}.log
cd; rm -rf $rundir
EOF

    if [ -n "$BENCH_DEPEND" ]; then
        jobid=$(qsub -W depend=afterany:$BENCH_DEPEND $jobfile)
    else
        jobid=$(qsub $jobfile)
    fi
    echo "submitted $TAG rep $rep ($BTYPE, $(basename $BUILD)): $jobid (commit $COMMIT)"
    BENCH_DEPEND=$jobid   # chain reps serially when dependencies are in use
done
[ -n "$jobid" ] && echo "LASTJOB $jobid"
