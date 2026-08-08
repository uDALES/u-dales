#!/bin/bash

set -e

if (( $# < 1 ))
then
 echo "The experiment directory must be set."
 echo "usage: FROM THE TOP LEVEL DIRECTORY run: bash ./u-dales/tools/archer_execute.sh <PATH_TO_CASE>"
 exit 1
fi

## go to experiment directory
pushd $1
inputdir=$(pwd)

## set experiment number via path
exp="${inputdir: -3}"

echo "experiment number: $exp"

## read in additional variables
if [ -f config.sh ]; then
    source config.sh
else
    echo "config.sh must be there inside $inputdir"
    exit 1
fi

## check if required variables are set
if [ -z $DA_WORKDIR ]; then
    echo "Output top-level directory DA_WORKDIR must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $DA_BUILD ]; then
    echo "Executable DA_BUILD must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $DA_TOOLSDIR ]; then
    echo "Script directory DA_TOOLSDIR must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $NNODE ]; then
    echo "Number of nodes NNODE must be set inside $inputdir/config.sh"
    echo "Product of NNODE and NCPU set in $inputdir/config.sh must be equal to the product of nprocx and nprocy set in $inputdir/namoptions.$exp"
    exit 1
fi;
if [ -z $NCPU ]; then
    echo "Number of CPU cores on each node NCPU must be set inside $inputdir/config.sh"
    echo "Product of NNODE and NCPU set in $inputdir/config.sh must be equal to the product of nprocx and nprocy set in $inputdir/namoptions.$exp"
    exit 1
fi;
if [ -z $WALLTIME ]; then
    echo "Wall clock time WALLTIME must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $QOS ]; then
    echo "Job quota QOS must be set inside $inputdir/config.sh"
    exit 1
fi;

## set the output directory
outdir=$DA_WORKDIR/$exp

## copy files to execution and output directory
mkdir -p $outdir
cp -r -P $inputdir/* $outdir
cp -r $DA_BUILD $outdir
pushd $outdir

echo "writing job.$exp.slurm"

## write new job.exp.slurm file
cat <<EOF > job.$exp.slurm
#!/bin/bash
#SBATCH --job-name=${exp}
#SBATCH --time=${WALLTIME}
#SBATCH --nodes=${NNODE}
#SBATCH --tasks-per-node=${NCPU}
#SBATCH --cpus-per-task=1
#SBATCH --account=n02-ASSURE
#SBATCH --partition=standard
#SBATCH --qos=${QOS}
module load epcc-job-env
export OMP_NUM_THREADS=1
srun --distribution=block:block --hint=nomultithread ./u-dales $outdir/namoptions.$exp > $outdir/output.$exp 2>&1
EOF

## submit job.exp file to queue
sbatch job.$exp.slurm

echo "job.$exp.slurm submitted."

popd
