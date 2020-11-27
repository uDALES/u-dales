#!/usr/bin/env bash
set -e

# Usage: run_udales.sh <NPROC> <BUILD_TYPE> <PATH_TO_CASE> <PATH_TO_NML>
# e.g. ./tools/run_udales.sh 1 Release examples/001 examples/001/namoptions.001

if [ ! -d src ]; then
    echo "Please run this script from the project folder"
    exit 1
fi

NPROC=$1
BUILD_TYPE=$2
PATH_TO_CASE=$3
PATH_TO_NML=$4

#ROOT_DIR=$(pwd)
ROOT_DIR=$PBS_O_WORKDIR
SIF_PATH=$ROOT_DIR/tools/singularity/image.sif
UDALES_EXE=$ROOT_DIR/build/$BUILD_TYPE/u-dales

singularity exec --containall \
    -B $ROOT_DIR:$ROOT_DIR \
    $SIF_PATH \
    bash -c "cd $ROOT_DIR/$PATH_TO_CASE && \
        mpiexec -n $NPROC $UDALES_EXE $ROOT_DIR/$PATH_TO_NML 2>&1 | tee -a run.log"
