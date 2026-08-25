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

# Usage: ./tools/build_executable.sh [icl, archer, cca, gpu, common] [debug, release]
#
# Optional environment overrides:
#   UDALES_BUILD_DIR          independent CMake build directory; overrides the
#                             default build/<cpu|gpu>/<debug|release> layout
#   UDALES_BUILD_JOBS         parallel build jobs (default: 8)
#   UDALES_FORTRAN_COMPILER   Fortran compiler/MPI wrapper

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

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <icl|archer|cca|gpu|common> <debug|release>"
    exit 2
fi

NPROC="${UDALES_BUILD_JOBS:-8}"
system="$1"
build_type="$2"

case "$build_type" in
    debug|release) ;;
    *)
        echo "Build type must be either 'debug' or 'release'"
        exit 2
        ;;
esac


if [ $system == "icl" ]
then
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

elif [ $system == "gpu" ]
then
    FC=mpif90
    NETCDF_DIR=/home/dipanjan/mygpu/netcdf-c-4.9.2/netcdfc
    NETCDF_FORTRAN_DIR=/home/dipanjan/mygpu/netcdf-fortran-4.6.1/netcdff

elif [ $system == "common" ]
then
    if [ -x /usr/bin/mpif90 ]; then
        FC=/usr/bin/mpif90
    else
        FC=
    fi
    NETCDF_DIR=
    NETCDF_FORTRAN_DIR=

else
    echo "This configuration is not avalable"
    exit 1
fi

FC="${UDALES_FORTRAN_COMPILER:-$FC}"


# Configure and Build
repo_root="$(pwd)"
if [ "$system" = "gpu" ]; then
    build_target="gpu"
else
    build_target="cpu"
fi
path_to_build_dir="${UDALES_BUILD_DIR:-$repo_root/build/$build_target/$build_type}"
if [[ "$path_to_build_dir" != /* ]]; then
    path_to_build_dir="$repo_root/$path_to_build_dir"
fi
echo "Build target:    $build_target"
echo "Build type:      $build_type"
echo "Build directory: $path_to_build_dir"
mkdir -p "$path_to_build_dir"
pushd "$path_to_build_dir"
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

if [ -n "$FC" ]; then
    cmake_args+=("-DCMAKE_Fortran_COMPILER=$FC")
fi

# A build directory remembers the compiler it was configured with. Hand it a
# different one and cmake deletes the cache and re-runs configure by itself -
# and CMAKE_BUILD_TYPE goes with it, because a command-line -D is a write into
# that cache rather than a standing instruction. The re-run does not re-read the
# command line, so it takes the "No build type selected, default to Release"
# branch in CMakeLists.txt, and a Release binary appears in a directory named
# debug with one STATUS line to say so.
#
# Notice the change here, where it can still be acted on, and clear the cache
# ourselves. cmake then configures once, from this command line, with the build
# type that was actually asked for. CMakeFiles goes too: its objects were built
# by the old compiler.
cmake_cache="$path_to_build_dir/CMakeCache.txt"
if [ -n "$FC" ] && [ -f "$cmake_cache" ]; then
    # The cache records this as STRING when it arrived by -D, FILEPATH when
    # cmake found it itself and UNINITIALIZED when -D carried no type, so the
    # pattern must not anchor on any one of them.
    cached_fc="$(sed -n 's/^CMAKE_Fortran_COMPILER:[^=]*=//p' "$cmake_cache")"
    if [ -n "$cached_fc" ] && [ "$cached_fc" != "$FC" ]; then
        echo "Fortran compiler changed for this build directory:"
        echo "  cached:    $cached_fc"
        echo "  requested: $FC"
        echo "Clearing the cmake cache so the build type survives the reconfigure."
        rm -rf "$cmake_cache" "$path_to_build_dir/CMakeFiles"
    fi
fi

cmake "${cmake_args[@]}" "$repo_root" 2>&1 | tee -a "$path_to_build_dir/config.log"
cmake --build . --parallel "$NPROC" 2>&1 | tee -a "$path_to_build_dir/build.log"
popd
