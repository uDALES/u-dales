#!/bin/bash

## automatically set the experiment number via path
cdir=$(pwd)
exp="${cdir: -3}"
echo "experiment number: ${exp}"

# some default parameters, changes can be made here
# This should really be an arg and not harcoded!
executable=${DA_TOPDIR}/u-dales/build/release/u-dales
utilspath=${DA_TOPDIR}/u-dales/tools/utils
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
echo "module load git/2.14.3 cmake/3.14.0 netcdf/4.5.2-fortran" >> job.${exp}

## copy files to execution and output directory
echo "mkdir ${outdir}" >> job.${exp}
echo "cp ${inputdir}/* ${outdir}" >> job.${exp}

## go to execution and output directory
echo "cd ${outdir}" >> job.${exp}

## execute program with mpi
echo "mpiexec -n $(( ${ncpu} * ${nnode} )) ${executable} ${outdir}/namoptions.${exp} > ${outdir}/output.${exp} 2>&1" >> job.${exp}

## merge output files from several cpus to one file
echo "${utilspath}/mergehelper.sh ${exp} " >> job.${exp}

## submit job.exp file to queue
qsub job.${exp}

echo "job.${exp} submitted."
