#!/usr/bin/env bash

set -xe

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    sudo apt update
    sudo apt install gfortran libopenmpi-dev openmpi-bin libnetcdf-dev libnetcdff-dev graphviz libfftw3-dev
elif [[ "$OSTYPE" == "darwin"* ]]; then
    brew install netcdf netcdf-fortran open-mpi graphviz fftw
    brew link gcc
    brew link --overwrite pkgconf --dry-run
fi
