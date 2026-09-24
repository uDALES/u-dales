#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$workdir/netcdf_c_s1_set_env_for_install.sh"

cd "$workdir/netcdf-c-4.10.1"

./configure \
    --prefix="$PWD/netcdfc" \
    --enable-netcdf-4 \
    --enable-parallel4 \
    --enable-shared \
    --enable-static \
    --disable-dap \
    --disable-byterange \
    --disable-libxml2 \
    CC="$CC" \
    CFLAGS="$CFLAGS" \
    CPPFLAGS="$CPPFLAGS" \
    LDFLAGS="$LDFLAGS"
