#!/usr/bin/env bash

# uDALES (https://github.com/uDALES/u-dales).

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Copyright (C) 2016-2019 the uDALES Team.

set -e

# Usage: udales_build.sh <NPROC> [Debug, Release]
# e.g. ./tools/singularity/udales_build.sh 2 Release

NPROC=$1
BUILD_TYPE=$2

# https://stackoverflow.com/a/246128/8893833
THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SIF_PATH=$THIS_DIR/image.sif
ROOT_DIR=$THIS_DIR/../..
PATH_TO_BUILD_DIR=$ROOT_DIR/build/$BUILD_TYPE

# Always clean up
rm -rf $PATH_TO_BUILD_DIR
mkdir -p $PATH_TO_BUILD_DIR

singularity exec --containall \
    -B $ROOT_DIR:$ROOT_DIR \
    $SIF_PATH \
    bash -c "cd $PATH_TO_BUILD_DIR && \
        cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE ../../ 2>&1 | tee -a config.log && \
        make -j$NPROC 2>&1 | tee -a $PATH_TO_BUILD_DIR/build.log"
