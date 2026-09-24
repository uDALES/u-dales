#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$workdir/hdf5_s1_set_env_for_install.sh"

cd "$workdir/hdf5-2.1.1"

cmake --fresh -S . -B build \
    -DCMAKE_INSTALL_PREFIX="$PWD/hdf5" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="$CC" \
    -DCMAKE_CXX_COMPILER="$CXX" \
    -DCMAKE_Fortran_COMPILER="$FC" \
    -DCMAKE_C_FLAGS="-fPIC" \
    -DCMAKE_CXX_FLAGS="-fPIC" \
    -DCMAKE_Fortran_FLAGS="-fPIC" \
    -DCMAKE_PREFIX_PATH="$EBROOTZLIB" \
    -DHDF5_ENABLE_PARALLEL=ON \
    -DHDF5_BUILD_FORTRAN=ON \
    -DHDF5_BUILD_CPP_LIB=OFF \
    -DHDF5_ENABLE_ZLIB_SUPPORT:BOOL=ON \
    -DZLIB_USE_EXTERNAL:BOOL=OFF \
    -DZLIB_INCLUDE_DIR:PATH="$EBROOTZLIB/include" \
    -DZLIB_LIBRARY:FILEPATH="$EBROOTZLIB/lib/libz.so" \
    -DBUILD_SHARED_LIBS=ON \
    -DBUILD_STATIC_LIBS=ON
