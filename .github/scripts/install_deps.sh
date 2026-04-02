#!/usr/bin/env bash

set -xe

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    sudo apt update
    sudo apt install gfortran libopenmpi-dev openmpi-bin libnetcdf-dev libnetcdff-dev graphviz libfftw3-dev
elif [[ "$OSTYPE" == "darwin"* ]]; then
    brew install pkgconf
    brew link --overwrite pkgconf
    brew install gcc netcdf netcdf-fortran open-mpi graphviz fftw
    if [[ -n "${GITHUB_PATH:-}" ]]; then
        echo "$(brew --prefix gcc)/bin" >> "$GITHUB_PATH"
        echo "$(brew --prefix open-mpi)/bin" >> "$GITHUB_PATH"
    fi
fi
