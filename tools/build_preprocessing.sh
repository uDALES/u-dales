#!/usr/bin/env bash

set -e

# Usage: ./tools/build_preprocessing.sh [common / icl]
if (( $# < 1 ))
then
 echo "The build type <common / icl> must be set."
 echo "usage: from being in u-dales directory run: tools/build_preprocessing.sh <build type>"
 exit 1
fi

if [ ! -d tools ]; then
    echo "Please run this script from being inside the u-dales folder"
    exit 1
fi

cd tools/View3D
mkdir build
cd build

system=$1
if [ $system == "icl" ]
then
    module load CMake/3.31.8-GCCcore-14.3.0
elif [ $system == "common" ]
then
    echo "Building View3D on local system."
else
    echo "This configuration is not avalable"
    exit 1
fi

cmake ..
echo "View3D configuration complete."

make
