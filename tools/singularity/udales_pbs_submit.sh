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

# This is an example script to submit a udales job using PBS. Modify as required.

#PBS -lwalltime=00:03:00
#PBS -lselect=1:ncpus=32:mem=32gb

set -ex

CASE_ID=201
BUILD_TYPE=Release

cd $PBS_O_WORKDIR
$PBS_O_WORKDIR/tools/singularity/udales_build.sh 2 $BUILD_TYPE
$PBS_O_WORKDIR/tools/singularity/udales_run.sh 32 $BUILD_TYPE examples/$CASE_ID namoptions.$CASE_ID
