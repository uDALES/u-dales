#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$workdir/netcdf_fortran_s1_set_env_for_install.sh"

cd "$workdir/netcdf-fortran-4.6.4"

make -j "${BUILD_JOBS:-8}"
make install

export NETCDFF_ROOT="$PWD/netcdff"
export PATH="$NETCDFF_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$NETCDFF_ROOT/lib:${LD_LIBRARY_PATH:-}"

"$NETCDFF_ROOT/bin/nf-config" --version
"$NETCDFF_ROOT/bin/nf-config" --all
