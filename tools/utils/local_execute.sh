#!/bin/bash

set -xe

this_dir=$(pwd)
exp_dir=$(dirname $this_dir)
proj_dir=$(dirname $exp_dir)

## Chnage the following paths based on your configuration
path_to_exe=${proj_dir}/build/u-dales
path_to_utils=${proj_dir}/tools/utils
NCPU=2

## automatically set the experiment number via path
inputdir=${this_dir}
exp="${inputdir: -3}"
echo "experiment number: ${exp}"
outdir=${proj_dir}/data/${exp}

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
