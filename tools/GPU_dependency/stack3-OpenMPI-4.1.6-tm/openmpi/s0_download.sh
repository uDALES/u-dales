#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
archive=openmpi-4.1.6.tar.bz2
url=https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.6.tar.bz2
source "$workdir/../download_helper.sh"

cd "$workdir"

fetch_archive "$url" "$archive"

if [[ ! -d openmpi-4.1.6 ]]; then
    tar -xjf "$archive"
fi

mkdir -p logs
