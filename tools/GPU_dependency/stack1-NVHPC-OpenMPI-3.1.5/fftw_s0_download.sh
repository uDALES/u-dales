#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
archive=fftw-3.3.11.tar.gz
url=https://fftw.org/pub/fftw/fftw-3.3.11.tar.gz
source "$workdir/download_helper.sh"

cd "$workdir"

fetch_archive "$url" "$archive"

if [[ ! -d fftw-3.3.11 ]]; then
    tar -xzf "$archive"
fi
