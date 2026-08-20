#!/bin/bash
# Test script for TEST_PERIODIC_EBCORR - host periodic energy-balance correction
#
# Runs the runmode 1009 checks on modforces::periodicEBcorr, which removes over
# the volume above the canopy, and through the top level, the heat and moisture
# the facets put into the air over a periodic domain.
#
# The routine writes that out as three coupled scalars - R_theta_scaled,
# phi_theta_t and the level count M behind them - so checking those expressions
# would only restate the code. What runs here instead are the properties they
# exist to produce:
#   - the column-integrated tendency equals the domain-mean flux
#     tot_Tflux/(xlen*ylen), for any sinkbase and any fraction. This is the
#     closure the ke/M scaling preserves, and what breaks when M miscounts;
#   - of that total, (1-fraction) leaves through the top level and the rest
#     through the volume sink - the Grylls (2021) split;
#   - the flux entering is the sum over ranks, not the local contribution;
#   - nothing at or below sinkbase moves, and nothing outside ib:ie, jb:je
#     moves, halo cells included;
#   - the tendency accumulates, is linear in the flux, and is switched off by
#     lperiodicEBcorr, ltempeq and lmoist.
#
# Two passes, because the MPI_ALLREDUCE over totheatflux is invisible on a
# single rank: there, the local flux and the reduced one are the same number.
# The second pass gives each rank a different flux and checks every rank
# against the sum, so dropping the reduction fails on every rank but the first.
#
# This covers the host branch only. The OpenACC kernels are covered by
# tests_cuda.f90::test_periodic_ebcorr, which runs under
# UDALES_RUN_CUDA_SELFTEST on a Debug GPU build. Running this runmode against a
# GPU build fails with an explanatory message rather than passing vacuously.

set -eu

# Configuration
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_ROOT}/tests/integration/common/runtime_modules.sh"
load_udales_runtime_modules

UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/cpu/debug/u-dales}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/cases/101}"
NAMELIST_SOURCE="${NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1009}"
NAMELIST_SOURCE_MPI="${NAMELIST_SOURCE_MPI:-${SCRIPT_DIR}/namoptions.1009.tworank}"
NAMELIST="${NAMELIST:-namoptions.101}"
TMPDIR_PARENT="${TMPDIR_PARENT:-}"
if [ -z "${MPIEXEC:-}" ] && command -v mpiifort >/dev/null 2>&1; then
    MPIEXEC="$(dirname "$(command -v mpiifort)")/mpiexec"
else
    MPIEXEC="${MPIEXEC:-mpiexec}"
fi
MPI_LAUNCH_EXTRA_ARGS="${MPI_LAUNCH_EXTRA_ARGS:-}"
MPI_VERSION_OUTPUT="$("$MPIEXEC" --version 2>/dev/null || true)"

if printf '%s\n' "$MPI_VERSION_OUTPUT" | grep -Eqi "Open MPI|OpenRTE"; then
    MPI_LAUNCH_EXTRA_ARGS="--oversubscribe ${MPI_LAUNCH_EXTRA_ARGS}"
    TMPDIR="${TMPDIR:-/tmp}"
    export TMPDIR
    export OMPI_MCA_prte_tmpdir_base="${OMPI_MCA_prte_tmpdir_base:-$TMPDIR}"
    export PRTE_MCA_prte_tmpdir_base="${PRTE_MCA_prte_tmpdir_base:-$TMPDIR}"
    if [ -z "$TMPDIR_PARENT" ]; then
        TMPDIR_PARENT="$TMPDIR"
    fi
fi

# Check if executable exists
if [ ! -f "$UDALES_BUILD" ]; then
    echo "ERROR: u-dales executable not found at: $UDALES_BUILD"
    echo "Please build u-dales in debug mode first:"
    echo "  ${REPO_ROOT}/tools/build_executable.sh common debug"
    exit 1
fi

# Check if case source exists
if [ ! -d "$CASE_SOURCE" ]; then
    echo "ERROR: Case source not found: $CASE_SOURCE"
    exit 1
fi

# $1 label, $2 namelist source
run_one () {
    variant="$1"
    namelist_src="$2"

    if [ ! -f "$namelist_src" ]; then
        echo "ERROR: Namelist source not found: $namelist_src"
        return 1
    fi

    nprocx="$(awk -F= '/^[[:space:]]*nprocx[[:space:]]*=/ {gsub(/[[:space:]]/,"",$2); print $2; exit}' "$namelist_src")"
    nprocy="$(awk -F= '/^[[:space:]]*nprocy[[:space:]]*=/ {gsub(/[[:space:]]/,"",$2); print $2; exit}' "$namelist_src")"
    if [ -z "$nprocx" ] || [ -z "$nprocy" ]; then
        echo "ERROR: Failed to determine nprocx/nprocy from $namelist_src"
        return 1
    fi
    nprocs="${NPROCS:-$((nprocx * nprocy))}"

    if [ -n "$TMPDIR_PARENT" ]; then
        RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/periodic_ebcorr-XXXXXX")"
    else
        RUN_DIR="$(mktemp -d)"
    fi

    cp -r "$CASE_SOURCE"/. "$RUN_DIR"/
    cp "$namelist_src" "$RUN_DIR/$NAMELIST"

    echo "=========================================="
    echo "Running TEST_PERIODIC_EBCORR (${variant})"
    echo "=========================================="
    echo "MPI processes: $nprocs"
    echo "Executable:    $UDALES_BUILD"
    echo "Case source:   $(basename "$CASE_SOURCE")"
    echo "Run directory: $RUN_DIR"
    echo "Namelist:      $NAMELIST"
    echo "=========================================="
    echo ""

    ( cd "$RUN_DIR" && "$MPIEXEC" $MPI_LAUNCH_EXTRA_ARGS -n "$nprocs" "$UDALES_BUILD" "$NAMELIST" )
    rc=$?

    rm -rf "$RUN_DIR"
    return $rc
}

EXIT_CODE=0
run_one "single rank" "$NAMELIST_SOURCE" || EXIT_CODE=$?
if [ "$EXIT_CODE" -eq 0 ]; then
    run_one "two ranks, unequal fluxes" "$NAMELIST_SOURCE_MPI" || EXIT_CODE=$?
fi

echo ""
echo "=========================================="
if [ "$EXIT_CODE" -eq 0 ]; then
    echo "Test completed successfully"
else
    echo "Test failed with exit code: $EXIT_CODE"
fi
echo "=========================================="

exit $EXIT_CODE
