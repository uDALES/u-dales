#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
archive=netcdf-fortran-4.6.4.tar.gz
url=https://downloads.unidata.ucar.edu/netcdf-fortran/4.6.4/netcdf-fortran-4.6.4.tar.gz
source "$workdir/download_helper.sh"

cd "$workdir"

fetch_archive "$url" "$archive"

if [[ ! -d netcdf-fortran-4.6.4 ]]; then
    tar -xzf "$archive"
fi
