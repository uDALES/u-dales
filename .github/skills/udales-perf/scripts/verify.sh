#!/usr/bin/env bash
# Verification runner for exp 900 (uDALES issue #330 performance work).
#
# Usage: verify.sh <tag> <binary> [nsteps]
#   <tag>     label, e.g. "base", "B0"; run dir = $EPHEMERAL/900/verify_<tag>
#   <binary>  path to the u-dales binary to run (archived bin/... recommended)
#   [nsteps]  number of fixed dt=0.05 steps (default 50; use 10 for Class-B
#             tolerance checks where chaotic growth must stay small)
#
# Runs the benchmark case for nsteps with restart output enabled (trestart hits
# just before the end), then LEAVES the run dir in place so the initd restart
# files (~1.2 GB/set, full double precision) can be compared:
#   Class A: for f in verify_A/initd*; do cmp $f verify_B/$(basename $f); done
#   Class B: compare_restarts.py verify_A verify_B
# Delete the run dirs promptly after comparison (disk rule: max two sets).

set -e

TAG=${1:?usage: verify.sh <tag> <binary> [nsteps]}
BIN=${2:?usage: verify.sh <tag> <binary> [nsteps]}
NSTEPS=${3:-50}

EXPDIR=~/udales/experiments/900
LOGDIR=$EXPDIR/logs
[ -x "$BIN" ] || { echo "ERROR: binary $BIN not found/executable." >&2; exit 1; }

RUNTIME=$(awk -v n=$NSTEPS 'BEGIN{printf "%.2f", n*0.05}')
TRESTART=$(awk -v n=$NSTEPS 'BEGIN{printf "%.2f", (n-1)*0.05}')

rundir=$EPHEMERAL/900/verify_${TAG}
rm -rf $rundir && mkdir -p $rundir
cp $EXPDIR/namoptions.900 $EXPDIR/*.inp.900 $rundir/
cp $EXPDIR/solid_*.txt $EXPDIR/fluid_boundary_*.txt $EXPDIR/facet_sections_*.txt $rundir/
cp $BIN $rundir/u-dales.verify
sed -i "s/^runtime .*/runtime      = $RUNTIME/;s/^trestart .*/trestart     = $TRESTART/" $rundir/namoptions.900

jobfile=$LOGDIR/job.verify_${TAG}
cat <<EOF > $jobfile
#!/bin/bash
#PBS -l walltime=00:20:00
#PBS -l select=1:ncpus=64:mpiprocs=64:mem=120gb:cpu_type=rome
#PBS -o $LOGDIR/
#PBS -e $LOGDIR/
module load intel/2021a netCDF/4.8.0-iimpi-2021a netCDF-Fortran/4.5.3-iimpi-2021a FFTW/3.3.9-intel-2021a
# some CX3 nodes firewall TCP-to-self, killing hydra bootstrap; loopback avoids it
export I_MPI_HYDRA_IFACE=lo
set -e
cd $rundir
mpiexec -n 64 $rundir/u-dales.verify $rundir/namoptions.900 > $rundir/output.900 2>&1
ls initd* > /dev/null   # fail loudly if no restart files were written
echo "verify_${TAG} done: \$(ls initd* | wc -l) restart files"
EOF

jobid=$(qsub $jobfile)
echo "submitted verify_${TAG} ($NSTEPS steps, $(basename $BIN)): $jobid"
