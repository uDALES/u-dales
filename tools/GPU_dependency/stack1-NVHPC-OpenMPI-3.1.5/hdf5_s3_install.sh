#!/usr/bin/env bash
set -euo pipefail

workdir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$workdir/hdf5-2.1.1"

cmake --build build --parallel 8

cmake --install build


ls hdf5/bin
ls hdf5/lib
ls hdf5/include

# hdf5/bin/h5pcc -showconfig
# hdf5/bin/h5pfc -showconfig
./hdf5/bin/h5pcc -showconfig | grep -i parallel
./hdf5/bin/h5pfc -showconfig | grep -i parallel

grep -Ei "Parallel HDF5|I/O filters" hdf5/lib/libhdf5.settings
