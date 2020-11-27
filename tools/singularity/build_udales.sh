#!/usr/bin/env bash
set -e

# Usage: ./tools/build_udales.sh [Debug, Release]

if [ ! -d src ]; then
    echo "Please run this script from the project folder"
    exit 1
fi

NPROC=2 # TODO: make into a arg var.

BUILD_TYPE=$1

#ROOT_DIR=$(pwd)
ROOT_DIR=$PBS_O_WORKDIR
SIF_PATH=$ROOT_DIR/tools/singularity/image.sif
PATH_TO_BUILD_DIR=$ROOT_DIR/build/$BUILD_TYPE

mkdir -p $PATH_TO_BUILD_DIR

singularity exec --containall \
    -B $ROOT_DIR:$ROOT_DIR \
    $SIF_PATH \
    bash -c "cd $PATH_TO_BUILD_DIR && \
        cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE ../../ 2>&1 | tee -a config.log && \
        make -j$NPROC 2>&1 | tee -a $PATH_TO_BUILD_DIR/build.log"
