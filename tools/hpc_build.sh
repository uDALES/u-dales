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

set -ex

# Usage: ./tools/hpc_build.sh [icl, archer, cca, common] [debug, release]

if [ ! -d src ]; then
    echo "Please run this script from the project folder"
    exit 1
fi

capitalize() {
    echo $* | sed -e "s/\b\(.\)/\u\1/g"
}

echo "--- Debug info ---"
echo "env: " `env`
echo "PATH: " ${PATH}

NPROC=2 # TODO: make into a arg var.
system=$1
build_type=$2


if [ $system == "icl" ]
then
    module load cmake/3.18.2 git/2.14.3 intel-suite/2020.2 mpi/intel-2019.8.254
    FC=mpiifort
    NETCDF_DIR=/apps/netcdf/4.4.1-c
    NETCDF_FORTRAN_DIR=/apps/netcdf/4.4.4-fortran

elif [ $system == "archer" ]
then
    module load cmake cray-hdf5 cray-netcdf cray-fftw
    FC=ftn
    NETCDF_DIR=$NETCDF_DIR
    NETCDF_FORTRAN_DIR=$NETCDF_FORTRAN_DIR
    #FFTW_DOUBLE_LIB=/opt/cray/pe/fftw/3.3.8.9/x86_rome/lib/libfftw3.so
    #FFTW_FLOAT_LIB=/opt/cray/pe/fftw/3.3.8.9/x86_rome/lib/libfftw3f.so

elif [ $system == "cca" ]
then
    module load git cmake
    FC=ftn
    NETCDF_DIR=$NETCDF_DIR
    NETCDF_FORTRAN_DIR=$NETCDF_DIR

elif [ $system == "common" ]
then
    FC=
    NETCDF_DIR=
    NETCDF_FORTRAN_DIR=

else
    echo "This configuration is not avalable"
    exit 1
fi


# Configure and Build
path_to_build_dir="$(pwd)/build/$build_type"
mkdir -p $path_to_build_dir
pushd $path_to_build_dir
cmake_build_type="$(capitalize $build_type)"
FC=$FC cmake -DNETCDF_DIR=$NETCDF_DIR \
             -DNETCDF_FORTRAN_DIR=$NETCDF_FORTRAN_DIR \
             -DCMAKE_BUILD_TYPE=cmake_build_type \
	     -DFFTW_DOUBLE_OPENMP_LIB=$FFTW_DOUBLE_LIB \
	     -DFFTW_FLOAT_OPENMP_LIB=$FFTW_FLOAT_LIB \
              ../../ 2>&1 | tee -a $path_to_build_dir/config.log
make -j$NPROC 2>&1 | tee -a $path_to_build_dir/build.log
popd
