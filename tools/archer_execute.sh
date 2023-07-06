#!/bin/bash

set -e

if (( $# < 1 ))
then
 echo "The experiment directory must be set."
 exit
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
fi

## set the output directory
outdir=$DA_WORKDIR/$exp

## copy files to execution and output directory
mkdir -p $outdir
cp -P $inputdir/* $outdir
cp $DA_BUILD $outdir
pushd $outdir

echo "writing job.$exp.slurm"

## write new job.exp.slurm file
echo "#!/bin/bash" > job.$exp.slurm

echo "#SBATCH --job-name=${exp}" >> job.$exp.slurm
echo "#SBATCH --time=${WALLTIME}" >> job.$exp.slurm
echo "#SBATCH --nodes=${NNODE}" >> job.$exp.slurm
echo "#SBATCH --tasks-per-node=${NCPU}" >> job.$exp.slurm
echo "#SBATCH --cpus-per-task=1" >> job.$exp.slurm
echo "#SBATCH --account=ecseae03" >> job.$exp.slurm
echo "#SBATCH --partition=standard" >> job.$exp.slurm
echo "#SBATCH --qos=${QOS}" >> job.$exp.slurm
echo "module load epcc-job-env" >> job.$exp.slurm
echo "export OMP_NUM_THREADS=1" >> job.$exp.slurm
echo "srun --distribution=block:block --hint=nomultithread ./u-dales $outdir/namoptions.$exp > $outdir/output.$exp 2>&1" >> job.$exp.slurm

## gather output files from cores in a single file
echo "$DA_TOOLSDIR/gather_outputs.sh $outdir " >> job.$exp.slurm

## submit job.exp file to queue
sbatch job.$exp.slurm

echo "job.$exp.slurm submitted."

popd
