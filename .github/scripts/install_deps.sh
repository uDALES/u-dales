#!/usr/bin/env bash

set -xe

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    sudo apt update
    sudo apt install gfortran libopenmpi-dev openmpi-bin libnetcdf-dev libnetcdff-dev graphviz libfftw3-dev
elif [[ "$OSTYPE" == "darwin"* ]]; then
    brew unlink pkg-config@0.29.2
    brew install pkgconf
    brew link --overwrite pkgconf
    brew install netcdf netcdf-fortran open-mpi graphviz fftw
    brew link gcc
fi
