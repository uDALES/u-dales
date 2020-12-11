#!/usr/bin/env bash

# This is an example script to submit a udales job using PBS. Modify as required.

#PBS -lwalltime=00:03:00
#PBS -lselect=1:ncpus=32:mem=32gb

set -ex

CASE_ID=201
BUILD_TYPE=Release

cd $PBS_O_WORKDIR
$PBS_O_WORKDIR/tools/singularity/udales_build.sh 2 $BUILD_TYPE
$PBS_O_WORKDIR/tools/singularity/udales_run.sh 32 $BUILD_TYPE examples/$CASE_ID namoptions.$CASE_ID
