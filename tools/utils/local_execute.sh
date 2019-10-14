#!/bin/bash

## automatically set the experiment number via path
cdir=$(pwd)
exp="${cdir: -3}"
echo "experiment number: ${exp}"

## some default parameters, changes can be made here
topdir=${DA_TOPDIR}
executable=${topdir}/u-dales/build/release/u-dales
utilspath=${topdir}/u-dales/tools/utils
ncpu=2

## directories
outdir=${DA_TOPDIR}/data/${exp}
inputdir=$(pwd)

echo "starting job.${exp}."

## copy files to output directory
mkdir $outdir
cp ${inputdir}/* $outdir

## go to execution and output directory
cd $outdir

## execute program with mpi
mpiexec -n ${ncpu} ${executable} ${inputdir}/namoptions.${exp} > ${outdir}/output.${exp} 2>&1

## merge output files from cores to one file
${utilspath}/mergehelper.sh ${exp}

echo "job.${exp} done."
