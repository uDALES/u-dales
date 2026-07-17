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

set -euo pipefail

# Usage: ./tools/build_executable.sh [icl, archer, cca, common] [debug, release]

if [ ! -d src ]; then
    echo "Please run this script from being inside the u-dales folder"
    exit 1
fi

capitalize() {
    echo $* | sed -e "s/\b\(.\)/\u\1/g"
}

#echo "--- Debug info ---"
#echo "env: " `env`
#echo "PATH: " ${PATH}

NPROC=4 # TODO: make into a arg var.
system=$1
build_type=$2


if [ $system == "icl" ]
then
    # NOTE (2026-07): The line below no longer loads reliably on CX3. It mixes
    # three toolchain years (intel/2025a, netCDF built with iimpi-2023a, FFTW
    # built with intel-2021a). From a pristine shell Lmod silently collapses the
    # whole stack down to the 2021a toolchain (FFTW-intel-2021a pins it), so you
    # actually get ifort 2021.2 and intel/2025a is a no-op. But if the session
    # already has the 2025a toolchain loaded, Lmod cannot downgrade
    # intel-compilers 2025->2021 in one transaction and aborts with
    # "intel/2025a cannot be loaded", leaving nothing loaded. Replaced with a
    # self-consistent 2021a Intel stack (the version it resolved to anyway),
    # which loads regardless of prior session state.
    # NOTE (2026-07): CMake and git must use GCCcore-10.3.0 to match intel/2021a.
    # Using GCCcore-13.3.0 variants (CMake/3.29.3, git/2.45.1) causes Lmod to
    # swap GCCcore-10.3.0 -> 13.3.0, which then prevents cURL/7.76.0-GCCcore-10.3.0
    # and zlib/1.2.11-GCCcore-10.3.0 (deps of netCDF/HDF5) from loading, producing
    # "dependent module(s) are not currently loaded" warnings.
    # module load intel/2025a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0
    module load intel/2021a netCDF/4.8.0-iimpi-2021a netCDF-Fortran/4.5.3-iimpi-2021a FFTW/3.3.9-intel-2021a CMake/3.20.1-GCCcore-10.3.0 git/2.32.0-GCCcore-10.3.0-nodocs
    FC=mpiifort
    NETCDF_DIR=/sw-eb/software/netCDF/4.8.0-iimpi-2021a
    NETCDF_FORTRAN_DIR=/sw-eb/software/netCDF-Fortran/4.5.3-iimpi-2021a

elif [ $system == "archer" ]
then
    module load cmake cray-hdf5 cray-netcdf cray-fftw
    FC=ftn
    #NETCDF_DIR=$NETCDF_DIR
    #NETCDF_FORTRAN_DIR=$NETCDF_FORTRAN_DIR
    NETCDF_DIR=/opt/cray/pe/netcdf/4.9.0.7/crayclang/14.0/
    NETCDF_FORTRAN_DIR=/opt/cray/pe/netcdf/4.9.0.7/crayclang/14.0/
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
cmake_args=(
    -DNETCDF_DIR="$NETCDF_DIR"
    -DNETCDF_FORTRAN_DIR="$NETCDF_FORTRAN_DIR"
    -DCMAKE_BUILD_TYPE="$cmake_build_type"
)

if [ -n "${FFTW_DOUBLE_LIB:-}" ]; then
    cmake_args+=("-DFFTW_DOUBLE_OPENMP_LIB=$FFTW_DOUBLE_LIB")
fi

if [ -n "${FFTW_FLOAT_LIB:-}" ]; then
    cmake_args+=("-DFFTW_FLOAT_OPENMP_LIB=$FFTW_FLOAT_LIB")
fi

FC=$FC cmake "${cmake_args[@]}" ../../ 2>&1 | tee -a $path_to_build_dir/config.log
make -j$NPROC 2>&1 | tee -a $path_to_build_dir/build.log
popd
