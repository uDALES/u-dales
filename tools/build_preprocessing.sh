#!/usr/bin/env bash

set -e

cd tools/View3D
mkdir build
cd build

cmake ..
echo "View3D configuration complete."

make
