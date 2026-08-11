#!/usr/bin/env bash

# Build isolated CPU and GPU executables for same-commit parity testing.

set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <Debug|Release>"
    exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

case "$1" in
    Debug|debug)
        BUILD_TYPE="debug"
        ;;
    Release|release)
        BUILD_TYPE="release"
        ;;
    *)
        echo "Build type must be Debug or Release"
        exit 2
        ;;
esac

CPU_SYSTEM="${UDALES_CPU_SYSTEM:-common}"
GPU_SYSTEM="${UDALES_GPU_SYSTEM:-gpu}"
CPU_BUILD_DIR="${UDALES_CPU_BUILD_DIR:-${REPO_ROOT}/build/cpu/${BUILD_TYPE}}"
GPU_BUILD_DIR="${UDALES_GPU_BUILD_DIR:-${REPO_ROOT}/build/gpu/${BUILD_TYPE}}"
CPU_COMPILER="${UDALES_CPU_FORTRAN_COMPILER:-/usr/bin/mpif90}"
GPU_COMPILER="${UDALES_GPU_FORTRAN_COMPILER:-$(command -v mpif90 || true)}"

if [ ! -x "$CPU_COMPILER" ]; then
    echo "ERROR: CPU MPI Fortran compiler is not executable: $CPU_COMPILER"
    echo "Set UDALES_CPU_FORTRAN_COMPILER to the required MPI compiler wrapper."
    exit 1
fi
if [ -z "$GPU_COMPILER" ] || [ ! -x "$GPU_COMPILER" ]; then
    echo "ERROR: GPU MPI Fortran compiler was not found."
    echo "Load NVHPC or set UDALES_GPU_FORTRAN_COMPILER explicitly."
    exit 1
fi
if [ "$(readlink -f "$CPU_COMPILER")" = "$(readlink -f "$GPU_COMPILER")" ]; then
    echo "ERROR: CPU and GPU compilers resolve to the same executable: $CPU_COMPILER"
    echo "Set UDALES_CPU_FORTRAN_COMPILER and UDALES_GPU_FORTRAN_COMPILER explicitly."
    exit 1
fi

echo "CPU compiler: $CPU_COMPILER"
echo "GPU compiler: $GPU_COMPILER"
echo "CPU build:    $CPU_BUILD_DIR"
echo "GPU build:    $GPU_BUILD_DIR"

(
    cd "$REPO_ROOT"
    UDALES_BUILD_DIR="$CPU_BUILD_DIR" \
    UDALES_FORTRAN_COMPILER="$CPU_COMPILER" \
        ./tools/build_executable.sh "$CPU_SYSTEM" "$BUILD_TYPE"

    UDALES_BUILD_DIR="$GPU_BUILD_DIR" \
    UDALES_FORTRAN_COMPILER="$GPU_COMPILER" \
        ./tools/build_executable.sh "$GPU_SYSTEM" "$BUILD_TYPE"
)

echo "Built CPU executable: $CPU_BUILD_DIR/u-dales"
echo "Built GPU executable: $GPU_BUILD_DIR/u-dales"
