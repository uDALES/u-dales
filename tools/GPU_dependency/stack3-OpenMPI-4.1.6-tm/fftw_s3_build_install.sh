#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$workdir/fftw_s1_set_env_for_install.sh"

cd "$workdir/fftw-3.3.11"

make -j "${BUILD_JOBS:-8}"
make install

export FFTW_ROOT="$PWD/fftw3"
export PATH="$FFTW_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$FFTW_ROOT/lib:${LD_LIBRARY_PATH:-}"
export CPATH="$FFTW_ROOT/include:${CPATH:-}"
export LIBRARY_PATH="$FFTW_ROOT/lib:${LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="$FFTW_ROOT/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

for library in libfftw3 libfftw3_mpi libfftw3_omp; do
    if [[ ! -e "$FFTW_ROOT/lib/$library.so" ]]; then
        echo "Error: expected library was not installed: $FFTW_ROOT/lib/$library.so" >&2
        exit 1
    fi
done

echo "FFTW installed at $FFTW_ROOT"
ls -l "$FFTW_ROOT/lib"/libfftw3*.so*
