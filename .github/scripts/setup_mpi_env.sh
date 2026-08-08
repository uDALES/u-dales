#!/usr/bin/env bash

set -eu

if [[ "${1:-}" == "macos-latest" ]]; then
    export PATH="$(brew --prefix gcc)/bin:$(brew --prefix open-mpi)/bin:$PATH"
    MPI_PREFIX="$(brew --prefix open-mpi)"
    export FC="$(brew --prefix gcc)/bin/gfortran"
    export OMPI_FC="$FC"
    export UDALES_MPI_FORTRAN_COMPILER="$MPI_PREFIX/bin/mpif90"
    export UDALES_MPI_C_COMPILER="$MPI_PREFIX/bin/mpicc"
    export UDALES_MPIEXEC="$MPI_PREFIX/bin/mpiexec"
    # Temporary compatibility workaround for branch-switch regression builds
    # against master on macOS. master still carries an old
    # downloadFindFFTW.cmake.in minimum CMake version, which newer Homebrew
    # CMake rejects. Remove this flag once the findFFTW helper fix in this
    # branch has been merged into master.
    export UDALES_CMAKE_POLICY_MINIMUM="-DCMAKE_POLICY_VERSION_MINIMUM=3.5"
else
    export UDALES_MPI_FORTRAN_COMPILER=/usr/bin/mpif90
    export UDALES_MPI_C_COMPILER=/usr/bin/mpicc
    export UDALES_MPIEXEC=/usr/bin/mpiexec
    export UDALES_CMAKE_POLICY_MINIMUM=""
fi

export UDALES_CMAKE_ARGS="$UDALES_CMAKE_POLICY_MINIMUM -DCMAKE_Fortran_COMPILER=$UDALES_MPI_FORTRAN_COMPILER -DMPI_Fortran_COMPILER=$UDALES_MPI_FORTRAN_COMPILER -DMPI_C_COMPILER=$UDALES_MPI_C_COMPILER -DMPIEXEC_EXECUTABLE=$UDALES_MPIEXEC"
export MPIEXEC="$UDALES_MPIEXEC"
