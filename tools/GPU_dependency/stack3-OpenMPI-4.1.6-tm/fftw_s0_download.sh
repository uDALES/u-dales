#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
archive=fftw-3.3.11.tar.gz
url=https://fftw.org/pub/fftw/fftw-3.3.11.tar.gz

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

if [[ ! -d fftw-3.3.11 ]]; then
    tar -xzf "$archive"
fi
