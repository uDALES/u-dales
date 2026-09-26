#!/bin/bash
# Warm-start facet-state continuity check (GitHub issue #269) on case 064.
#
# Before the fix for #269, a warm start silently re-seeded the facet temperatures
# from Tfacinit.inp.064 because no facet state was ever checkpointed. This suite
# runs a reference chain and a restarted chain and asserts that the restarted
# facT trajectory continues the reference one.

set -u

if ! command -v module >/dev/null 2>&1 && [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/debug/u-dales}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/cases/064}"
CASE_ID=064
NAMELIST="namoptions.${CASE_ID}"
TEMPLATE="${SCRIPT_DIR}/${NAMELIST}.template"
TMPDIR_PARENT="${TMPDIR_PARENT:-/tmp}"

# The solver modules must match the stack the executable was built with. See
# README.md: tools/build_executable.sh icl currently uses the 2021a stack while
# tools/hpc_execute.sh uses the 2025a stack, so this is deliberately overridable.
UDALES_RUNTIME_MODULES="${UDALES_RUNTIME_MODULES:-intel/2025a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0}"
# Modules providing a python3 with netCDF4 and numpy for the comparator. Leave
# empty to use whatever python3 is already on PATH.
UDALES_PYTHON_MODULES="${UDALES_PYTHON_MODULES:-}"
PYTHON="${PYTHON:-python3}"

# Simulation schedule. The reference runs for 2*T, checkpoints at T, and the
# restarted run continues from T for another T.
T_RESTART="${T_RESTART:-10.}"
T_TOTAL="${T_TOTAL:-20.}"
DTEB="${DTEB:-5}"
NFACLYRS=5

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
fi

if [ ! -f "$UDALES_BUILD" ]; then
    echo "ERROR: u-dales executable not found at: $UDALES_BUILD"
    exit 1
fi

if [ ! -d "$CASE_SOURCE" ]; then
    echo "ERROR: Case source not found: $CASE_SOURCE"
    exit 1
fi

if [ ! -f "$TEMPLATE" ]; then
    echo "ERROR: Namelist template not found: $TEMPLATE"
    exit 1
fi

ROOT_RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/warmstart-facets-XXXXXX")"
KEEP_RUN_DIR=0
cleanup() {
    if [ "$KEEP_RUN_DIR" -eq 0 ]; then
        rm -rf "$ROOT_RUN_DIR"
    else
        echo "Preserving run directory: $ROOT_RUN_DIR"
    fi
}
trap cleanup EXIT

failures=0

fail() {
    echo "FAIL: $*"
    failures=$((failures + 1))
    KEEP_RUN_DIR=1
}

# stage <label> <nprocx> <nprocy> <runmode> <runtime> <trestart> <lwarmstart> <startfile> [extra &RUN lines]
#
# The extra &RUN lines are left empty for every run except the lreadfacT
# negative control, so that the generated namelist stays readable by a uDALES
# build that predates lreadfacT. That keeps this suite able to run against
# unfixed code, where it must fail on the facT comparison rather than on a
# namelist error.
stage() {
    local label="$1" npx="$2" npy="$3" runmode="$4" runtime="$5" trestart="$6"
    local lwarmstart="$7" startfile="$8" extra_run="${9:-}"
    local run_dir="${ROOT_RUN_DIR}/${label}"

    mkdir -p "$run_dir"
    cp -r "$CASE_SOURCE"/. "$run_dir"/
    sed -e "s|@NPROCX@|${npx}|g" \
        -e "s|@NPROCY@|${npy}|g" \
        -e "s|@RUNMODE@|${runmode}|g" \
        -e "s|@RUNTIME@|${runtime}|g" \
        -e "s|@TRESTART@|${trestart}|g" \
        -e "s|@LWARMSTART@|${lwarmstart}|g" \
        -e "s|@STARTFILE@|${startfile}|g" \
        -e "s|@EXTRA_RUN@|${extra_run}|g" \
        -e "s|@DTEB@|${DTEB}|g" \
        -e "s|@NFACLYRS@|${NFACLYRS}|g" \
        "$TEMPLATE" > "${run_dir}/${NAMELIST}"
    echo "$run_dir"
}

# run_udales <run_dir> <nprocs>  -> returns solver exit code, log in run.log
run_udales() {
    local run_dir="$1" np="$2"
    local module_setup=""
    if [ -f /etc/profile.d/modules.sh ]; then
        module_setup="source /etc/profile.d/modules.sh >/dev/null 2>&1 || true; "
    fi
    bash -c "${module_setup}\
if command -v module >/dev/null 2>&1 && [ -n \"\$UDALES_RUNTIME_MODULES\" ]; then module load \$UDALES_RUNTIME_MODULES >/dev/null 2>&1; fi; \
export HDF5_USE_FILE_LOCKING=FALSE; \
export FOR_DISABLE_DIAGNOSTIC_DISPLAY=TRUE; \
cd '${run_dir}' && '${MPIEXEC}' ${MPI_LAUNCH_EXTRA_ARGS} -n ${np} '${UDALES_BUILD}' '${NAMELIST}' > run.log 2>&1"
}
export UDALES_RUNTIME_MODULES

# earliest_restart <run_dir> -> the initd basename (without rank/expnr) of the
# first checkpoint written, i.e. the one at t = T_RESTART
earliest_restart() {
    local run_dir="$1"
    ls "${run_dir}"/initd*_000_000."${CASE_ID}" 2>/dev/null \
        | xargs -r -n1 basename \
        | sed -e "s/_000_000\.${CASE_ID}\$//" \
        | sort \
        | head -1
}

# run_chain <tag> <nprocx> <nprocy>
# Reference cold start followed by a warm-started continuation.
run_chain() {
    local tag="$1" npx="$2" npy="$3"
    local np=$((npx * npy))

    local ref_dir
    ref_dir="$(stage "ref_${tag}" "$npx" "$npy" 1 "$T_TOTAL" "$T_RESTART" .false. '')"
    echo "--- reference chain [${tag}] (${npx}x${npy}), runtime=${T_TOTAL}, trestart=${T_RESTART}"
    if ! run_udales "$ref_dir" "$np"; then
        fail "reference run [${tag}] did not complete"
        tail -n 60 "${ref_dir}/run.log" || true
        return 1
    fi

    local stem
    stem="$(earliest_restart "$ref_dir")"
    if [ -z "$stem" ]; then
        fail "reference run [${tag}] wrote no initd restart files"
        return 1
    fi
    if [ ! -f "${ref_dir}/${stem/initd/initf}.${CASE_ID}" ]; then
        fail "reference run [${tag}] wrote no facet restart file ${stem/initd/initf}.${CASE_ID}"
    fi
    echo "--- checkpoint stem: ${stem}"

    local rst_dir
    rst_dir="$(stage "rst_${tag}" "$npx" "$npy" 2 "$T_RESTART" "$T_TOTAL" .true. \
        "${stem}_xxx_xxx.${CASE_ID}")"
    cp "${ref_dir}/${stem}"_*_*."${CASE_ID}" "$rst_dir"/
    cp "${ref_dir}/${stem/initd/initf}.${CASE_ID}" "$rst_dir"/ 2>/dev/null || true
    echo "--- restart chain [${tag}], runtime=${T_RESTART}"
    if ! run_udales "$rst_dir" "$np"; then
        fail "restart run [${tag}] did not complete"
        tail -n 60 "${rst_dir}/run.log" || true
        return 1
    fi
    if ! grep -q "Restored facet state from ${stem/initd/initf}" "${rst_dir}/run.log"; then
        fail "restart run [${tag}] did not report restoring the facet state"
        grep -i "facet" "${rst_dir}/run.log" || true
    fi

    printf '%s\n' "$ref_dir" "$rst_dir" "$stem" > "${ROOT_RUN_DIR}/chain_${tag}.paths"
    return 0
}

echo "=========================================="
echo "Warm-start facet continuity test (issue #269)"
echo "Executable   : $UDALES_BUILD"
echo "Case source  : $CASE_SOURCE"
echo "Run directory: $ROOT_RUN_DIR"
echo "=========================================="

run_chain 2x2 2 2 || true
run_chain 2x1 2 1 || true

# ---------------------------------------------------------------------------
# Negative control: lreadfacT = .false. reproduces the pre-fix behaviour, i.e.
# the facet temperatures are re-seeded from Tfacinit.inp. This keeps the
# comparator honest - if it passed here the assertions would be vacuous.
# ---------------------------------------------------------------------------
noread_dir=""
if [ -f "${ROOT_RUN_DIR}/chain_2x2.paths" ]; then
    ref_dir="$(sed -n 1p "${ROOT_RUN_DIR}/chain_2x2.paths")"
    stem="$(sed -n 3p "${ROOT_RUN_DIR}/chain_2x2.paths")"
    noread_dir="$(stage noread_2x2 2 2 2 "$T_RESTART" "$T_TOTAL" .true. \
        "${stem}_xxx_xxx.${CASE_ID}" 'lreadfacT    = .false.')"
    cp "${ref_dir}/${stem}"_*_*."${CASE_ID}" "$noread_dir"/
    cp "${ref_dir}/${stem/initd/initf}.${CASE_ID}" "$noread_dir"/ 2>/dev/null || true
    echo "--- negative control: lreadfacT = .false."
    if ! run_udales "$noread_dir" 4; then
        fail "lreadfacT = .false. run did not complete"
        tail -n 60 "${noread_dir}/run.log" || true
        noread_dir=""
    elif grep -q "Restored facet state from" "${noread_dir}/run.log"; then
        fail "lreadfacT = .false. still read the facet restart file"
    fi

    # -----------------------------------------------------------------------
    # Backward compatibility: a restart directory without an initf file (as
    # written by any older uDALES) must still run, and must warn loudly.
    # -----------------------------------------------------------------------
    compat_dir="$(stage compat_2x2 2 2 2 1. "$T_TOTAL" .true. \
        "${stem}_xxx_xxx.${CASE_ID}")"
    cp "${ref_dir}/${stem}"_*_*."${CASE_ID}" "$compat_dir"/
    rm -f "${compat_dir}/initf"*
    echo "--- backward compatibility: no initf file present"
    if ! run_udales "$compat_dir" 4; then
        fail "restart without a facet restart file did not complete"
        tail -n 60 "${compat_dir}/run.log" || true
    else
        if ! grep -q "facet restart file .* not found" "${compat_dir}/run.log"; then
            fail "restart without a facet restart file did not warn about the fallback"
        fi
        if ! grep -q "re-initialised from Tfacinit.inp" "${compat_dir}/run.log"; then
            fail "fallback warning does not name Tfacinit.inp"
        fi
        if grep -q "Restored facet state from" "${compat_dir}/run.log"; then
            fail "restart without a facet restart file claimed to restore facet state"
        fi
    fi
fi

# ---------------------------------------------------------------------------
# Numerical comparison of the facT trajectories.
# ---------------------------------------------------------------------------
# This runs even if an earlier check failed, so that a run against unfixed code
# reports the actual facT divergence rather than stopping at the log checks.
if [ -f "${ROOT_RUN_DIR}/chain_2x2.paths" ]; then
    module_setup=""
    if [ -f /etc/profile.d/modules.sh ]; then
        module_setup="source /etc/profile.d/modules.sh >/dev/null 2>&1 || true; "
    fi
    ref22="$(sed -n 1p "${ROOT_RUN_DIR}/chain_2x2.paths")"
    rst22="$(sed -n 2p "${ROOT_RUN_DIR}/chain_2x2.paths")"

    compare_args=(
        --reference "${ref22}/facT.${CASE_ID}.nc"
        --restart "${rst22}/facT.${CASE_ID}.nc"
    )
    if [ -f "${ROOT_RUN_DIR}/chain_2x1.paths" ]; then
        rst21="$(sed -n 2p "${ROOT_RUN_DIR}/chain_2x1.paths")"
        compare_args+=(--other-decomposition "${rst21}/facT.${CASE_ID}.nc")
    else
        fail "the 2x1 chain did not produce outputs, so the MPI shape check was skipped"
    fi
    if [ -n "$noread_dir" ]; then
        compare_args+=(--unfixed "${noread_dir}/facT.${CASE_ID}.nc")
    else
        fail "the lreadfacT negative control did not produce outputs"
    fi

    export COMPARE_ARGS="${compare_args[*]}"
    export UDALES_PYTHON_MODULES PYTHON
    bash -c "${module_setup}\
if command -v module >/dev/null 2>&1 && [ -n \"\$UDALES_PYTHON_MODULES\" ]; then module load \$UDALES_PYTHON_MODULES >/dev/null 2>&1; fi; \
\"\$PYTHON\" '${SCRIPT_DIR}/compare_facT.py' \$COMPARE_ARGS" || fail "facT comparison"
else
    fail "the 2x2 reference/restart chain did not produce outputs to compare"
fi

echo "=========================================="
if [ "$failures" -eq 0 ]; then
    echo "Warm-start facet continuity test PASSED"
    echo "=========================================="
    exit 0
else
    echo "Warm-start facet continuity test FAILED (${failures} failure(s))"
    echo "=========================================="
    exit 1
fi
