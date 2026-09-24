#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
archive=netcdf-fortran-4.6.4.tar.gz
url=https://downloads.unidata.ucar.edu/netcdf-fortran/4.6.4/netcdf-fortran-4.6.4.tar.gz

cd "$workdir"

if [[ ! -f "$archive" ]]; then
    if command -v curl >/dev/null 2>&1; then
        curl -fL -o "$archive" "$url"
    elif command -v wget >/dev/null 2>&1; then
        wget -O "$archive" "$url"
    else
        echo "Error: neither curl nor wget is available." >&2
        exit 1
    fi
fi

if [[ ! -d netcdf-fortran-4.6.4 ]]; then
    tar -xzvf "$archive"
fi
