#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$workdir/netcdf_c_s1_set_env_for_install.sh"

cd "$workdir/netcdf-c-4.10.1"

make -j "${BUILD_JOBS:-8}"
make install

export NETCDF_ROOT="$PWD/netcdfc"
export PATH="$NETCDF_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$NETCDF_ROOT/lib:${LD_LIBRARY_PATH:-}"

"$NETCDF_ROOT/bin/nc-config" --version
"$NETCDF_ROOT/bin/nc-config" --has-hdf5
"$NETCDF_ROOT/bin/nc-config" --has-nc4
"$NETCDF_ROOT/bin/nc-config" --has-parallel4
