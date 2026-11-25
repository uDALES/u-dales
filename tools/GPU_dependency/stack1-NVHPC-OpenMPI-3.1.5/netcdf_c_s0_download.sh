#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
archive=netcdf-c-4.10.1.tar.gz
url=https://downloads.unidata.ucar.edu/netcdf-c/4.10.1/netcdf-c-4.10.1.tar.gz
source "$workdir/download_helper.sh"

cd "$workdir"

fetch_archive "$url" "$archive"

if [[ ! -d netcdf-c-4.10.1 ]]; then
    tar -xzf "$archive"
fi
