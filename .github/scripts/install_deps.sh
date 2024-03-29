#!/usr/bin/env bash

set -xe

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    sudo apt update
    sudo apt install gfortran libopenmpi-dev openmpi-bin libnetcdf-dev libnetcdff-dev graphviz
elif [[ "$OSTYPE" == "darwin"* ]]; then
    brew install netcdf netcdf-fortran open-mpi graphviz
    brew link gcc
fi
