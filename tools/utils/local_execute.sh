#!/bin/bash

set -e


if (( $# < 1 ))
then
 echo "The experiment ID must be set."
 exit
fi

exp=$(printf "%03.0f" $1)    #pad argument 1 (target simulation) with zeros

if [ -z $NCPU ]; then
  NCPU=1
fi;
if [ -z $DA_WORKDIR ]; then
  echo "DA_WORKDIR must be set"
  exit
fi;
if [ -z $DA_EXPDIR ]; then
  echo "DA_EXPDIR must be set"
  exit
fi;
if [ -z $DA_UTILSDIR ]; then
  echo "DA_UTILSDIR must be set"
  exit
fi;
if [ -z $DA_BUILD ]; then
  echo "DA_BUILD must be set"
  exit
fi;

## automatically set the experiment number via path
inputdir=${DA_EXPDIR}/${exp}
echo "experiment number: ${exp}"
outdir=${DA_WORKDIR}/${exp}

echo "starting job.${exp}."

## copy files to output directory
mkdir -p $outdir
cp ${inputdir}/* $outdir

## go to execution and output directory
pushd $outdir

## execute program with mpi
mpiexec -n ${NCPU} ${DA_BUILD} namoptions.${exp} > output.${exp} 2>&1

## merge output files from cores to one file
${DA_UTILSDIR}/mergehelper.sh ${outdir}

popd

echo "job.${exp} done."
