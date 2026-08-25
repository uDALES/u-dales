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
        CMAKE_TYPE="Debug"
        ;;
    Release|release)
        BUILD_TYPE="release"
        CMAKE_TYPE="Release"
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

# Assert each directory really was configured as the type that was asked for.
#
# It is possible for it not to be, and the failure is silent. Hand cmake a
# compiler that differs from the one a directory was configured with and it
# deletes the cache to start over; CMAKE_BUILD_TYPE lives in that cache, so it
# goes too, and the reconfigure falls through to the "default to Release" branch
# in CMakeLists.txt. That is how a Release binary comes to sit in
# build/cpu/debug - and the parity matrix then compares it against a Debug GPU
# binary without complaint, because it compares outputs and never asks how they
# were built. --require-debug-selftest does not cover it either: that reads the
# GPU run's log, and nothing reads the CPU side at all.
#
# What is lost is the whole point of a Debug reference: -fcheck=all,
# -finit-real=nan, -ffpe-trap and UDALES_DEBUG. The comparison still passes,
# because -O3 against -O0 moves the CPU result by far less than the 1e-6
# tolerance. So it stays quiet until somebody wonders why bounds checking never
# catches anything.
#
# build_executable.sh now clears the cache itself when it sees the compiler
# change, so this should not fire. It is here because the check is one grep and
# the alternative is trusting that it never happens again.
check_build_type() {
    label="$1"
    dir="$2"
    cache="$dir/CMakeCache.txt"

    if [ ! -f "$cache" ]; then
        echo "ERROR: $label build has no CMakeCache.txt: $cache"
        exit 1
    fi

    actual="$(sed -n 's/^CMAKE_BUILD_TYPE:[^=]*=//p' "$cache")"
    if [ "$actual" != "$CMAKE_TYPE" ]; then
        echo "ERROR: the $label build directory is configured as '$actual', not '$CMAKE_TYPE'."
        echo "  $cache"
        echo "  Nothing downstream checks this: the parity matrix compares outputs,"
        echo "  not how the binaries producing them were built."
        echo "  Remove $dir and run this script again."
        exit 1
    fi
    echo "$label build type: $actual"
}

check_build_type CPU "$CPU_BUILD_DIR"
check_build_type GPU "$GPU_BUILD_DIR"

echo "Built CPU executable: $CPU_BUILD_DIR/u-dales"
echo "Built GPU executable: $GPU_BUILD_DIR/u-dales"
