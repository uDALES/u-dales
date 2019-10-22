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
if [ -z $DA_TOPDIR ]; then
  echo "DA_TOPDIR must be set"
  exit
fi;

path_to_exe=${DA_TOPDIR}/u-dales/build/u-dales
path_to_utils=${DA_TOPDIR}/u-dales/tools/utils


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
mpiexec -n ${NCPU} ${path_to_exe} namoptions.${exp} > output.${exp} 2>&1

## merge output files from cores to one file
${path_to_utils}/mergehelper.sh ${exp}

popd

echo "job.${exp} done."
