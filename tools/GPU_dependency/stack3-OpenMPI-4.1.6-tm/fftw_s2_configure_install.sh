#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$workdir/fftw_s1_set_env_for_install.sh"

cd "$workdir/fftw-3.3.11"

./configure \
    --prefix="$PWD/fftw3" \
    --enable-shared \
    --enable-static \
    --enable-openmp \
    --enable-mpi \
    CC="$CC" \
    CXX="$CXX" \
    FC="$FC" \
    F77="$F77" \
    MPICC="$MPICC" \
    CFLAGS="$CFLAGS" \
    CXXFLAGS="$CXXFLAGS" \
    FFLAGS="$FFLAGS" \
    FCFLAGS="$FCFLAGS"
