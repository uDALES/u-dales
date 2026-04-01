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
mkdir -p build/src

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

if command -v cmake >/dev/null 2>&1
then
    cd build
    cmake ..
    echo "View3D configuration complete."
    make
    if [ -f src/view3d ]
    then
        echo "View3D executable available at tools/View3D/build/src/view3d"
        exit 0
    fi
    echo "CMake build completed but build/src/view3d was not created."
    exit 1
fi

echo "cmake not found; falling back to src/Makefile"
cd src
make clean
make
cp View3D ../build/src/view3d
make clean
echo "View3D executable available at tools/View3D/build/src/view3d"
