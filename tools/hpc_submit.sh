#!/usr/bin/env bash

#PBS -lwalltime=2:00:00
#PBS -lselect=1:ncpus=48:mem=124gb

set -ex

CASE_ID=201
BUILD_TYPE=Release

cd $PBS_O_WORKDIR
$PBS_O_WORKDIR/tools/build_udales.sh $BUILD_TYPE
$PBS_O_WORKDIR/tools/run_udales.sh 48 $BUILD_TYPE examples/$CASE_ID examples/$CASE_ID /namoptions.$CASE_ID
