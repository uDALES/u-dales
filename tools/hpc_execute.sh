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
outdir=$EPHEMERAL/$exp

echo "writing job.$exp."

## write new job.exp file for HPC
echo "#PBS -l walltime=${WALLTIME}" > job.$exp
echo "#PBS -l select=${NNODE}:ncpus=${NCPU}:mem=${MEM}" >> job.$exp

## load modules required for the execution
echo "module load intel-suite/2017.6 mpi/intel-2018 cmake/3.14.0 git/2.14.3" >> job.$exp

## copy files to execution and output directory
echo "mkdir -p $outdir" >> job.$exp
echo "cp $inputdir/* $outdir" >> job.$exp

## go to execution and output directory
echo "pushd $outdir" >> job.$exp

## execute program with mpi
echo "mpiexec -n $(( $NCPU * $NNODE )) $DA_BUILD $outdir/namoptions.$exp > $outdir/output.$exp 2>&1" >> job.$exp

## gather output files from cores in a single file
echo "$DA_TOOLSDIR/gather_outputs.sh $outdir " >> job.$exp

## submit job.exp file to queue
qsub job.$exp

echo "job.$exp submitted."
