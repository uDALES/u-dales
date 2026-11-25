#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
archive=hdf5-2.1.1.tar.gz
url=https://github.com/HDFGroup/hdf5/releases/download/2.1.1/hdf5-2.1.1.tar.gz
source "$workdir/download_helper.sh"

cd "$workdir"

fetch_archive "$url" "$archive"

if [[ ! -d hdf5-2.1.1 ]]; then
    tar -xzf "$archive"
fi
