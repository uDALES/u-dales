#!/usr/bin/env bash

set -xe

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    sudo apt install gfortran libopenmpi-dev openmpi-bin libnetcdf-dev libnetcdff-dev
elif [[ "$OSTYPE" == "darwin"* ]]; then
    brew install netcdf open-mpi graphviz
    brew unlink gcc
    brew link gcc@9
fi
