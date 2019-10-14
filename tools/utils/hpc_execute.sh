#!/bin/bash

## automatically set the experiment number via path
cdir=$(pwd)
exp="${cdir: -3}"
echo "experiment number: ${exp}"

## some default parameters, changes can be made here
topdir=${DA_TOPDIR}
executable=${topdir}/u-dales/build/release/u-dales
utilspath=${topdir}/u-dales/tools/utils
ncpu=12

## HPC specific
nnode=5
walltime="24:00:00"
mem="20gb"

## directories
outdir=${EPHEMERAL}/${exp}
inputdir=$(pwd)

## write new job.exp file for HPC
echo "#PBS -l walltime=${walltime}" > job.${exp}
echo "#PBS -l select=${nnode}:ncpus=${ncpu}:mem=${mem}" >> job.${exp}

## load modules required for the execution
echo "module load intel-suite/2017.6 mpi/intel-2018" >> job.${exp}

## copy files to temporary directory
echo "mkdir $EPHEMERAL/${exp}" >> job.${exp}
echo "cp ${inputdir}/* $EPHEMERAL/${exp}" >> job.${exp}

## go to execution directory. This is important, otherwise warmstart inps cannot be read!
echo "cd $EPHEMERAL/${exp}" >> job.${exp}

## execute program with mpi
echo "mpiexec -n $(( ${ncpu} * ${nnode} )) ${executable} $EPHEMERAL/${exp}/namoptions.${exp} > $EPHEMERAL/${exp}/output.${exp} 2>&1" >> job.${exp}

## merge output files from cores to one file
echo "${utilspath}/mergehelper.sh ${exp} " >> job.${exp}

## submit job file to queue
qsub job.${exp}

echo "job.${exp} submitted."