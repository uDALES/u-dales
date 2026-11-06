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

# Simulation schedule. The reference runs for T_TOTAL and checkpoints at
# T_RESTART; the restarted run continues from there to T_TOTAL.
# T_RESTART is a multiple of DTEB, so that checkpoint coincides with an
# energy-balance update. T_RESTART_MID deliberately is not: the checkpoint then
# falls inside an EB interval, which exercises the interval state (tEB, tnextEB
# and the accumulated fluxes) carried in the initf file.
# Without it the restarted run re-anchors the EB clock at the restart time and
# drops the flux accumulated since the last EB update.
T_RESTART="${T_RESTART:-10.}"
T_RESTART_MID="${T_RESTART_MID:-12.}"
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

# ---------------------------------------------------------------------------
# Optional GPU build. When UDALES_GPU_BUILD is set, the same chains are also run
# on that executable and compared with the CPU trajectories:
#   - a GPU warm start from the CPU checkpoint (the facet state restored on the
#     host has to reach the device before the first step),
#   - a GPU reference followed by a GPU warm start from the GPU checkpoint (the
#     device-side accumulators have to reach the restart file and come back),
#   - the same with the checkpoint inside an energy-balance interval.
# Each rank needs its own GPU, so the GPU decomposition follows
# UDALES_GPU_MAX_RANKS (>= 4: 2x2, 2-3: 2x1, 1: 1x1) and the launch goes
# through UDALES_GPU_LAUNCH_WRAPPER (tools/bind.sh by default), which maps the
# local rank to a device. UDALES_GPU_SELFTEST=1 switches the CUDA device
# self-tests on in the GPU runs (Debug builds only); among them is the check
# that the device facet state equals the host's after initialisation, which is
# what a warm start exercises.
# ---------------------------------------------------------------------------
UDALES_GPU_BUILD="${UDALES_GPU_BUILD:-}"
UDALES_GPU_MPIEXEC="${UDALES_GPU_MPIEXEC:-$MPIEXEC}"
UDALES_GPU_MAX_RANKS="${UDALES_GPU_MAX_RANKS:-4}"
UDALES_GPU_RUNTIME_MODULES="${UDALES_GPU_RUNTIME_MODULES:-}"
UDALES_GPU_SELFTEST="${UDALES_GPU_SELFTEST:-0}"
UDALES_GPU_LAUNCH_WRAPPER="${UDALES_GPU_LAUNCH_WRAPPER-${REPO_ROOT}/tools/bind.sh}"
GPU_MPI_LAUNCH_EXTRA_ARGS="${GPU_MPI_LAUNCH_EXTRA_ARGS:-}"
GPU_NPX=2
GPU_NPY=2
if [ -n "$UDALES_GPU_BUILD" ]; then
    if [ ! -f "$UDALES_GPU_BUILD" ]; then
        echo "ERROR: GPU u-dales executable not found at: $UDALES_GPU_BUILD"
        exit 1
    fi
    if [ -n "$UDALES_GPU_LAUNCH_WRAPPER" ] && [ ! -x "$UDALES_GPU_LAUNCH_WRAPPER" ]; then
        echo "ERROR: GPU launch wrapper not found or not executable: $UDALES_GPU_LAUNCH_WRAPPER"
        exit 1
    fi
    if "$UDALES_GPU_MPIEXEC" --version 2>/dev/null | grep -Eqi "Open MPI|OpenRTE"; then
        GPU_MPI_LAUNCH_EXTRA_ARGS="--oversubscribe ${GPU_MPI_LAUNCH_EXTRA_ARGS}"
    fi
    if [ "$UDALES_GPU_MAX_RANKS" -ge 4 ]; then
        GPU_NPX=2; GPU_NPY=2
    elif [ "$UDALES_GPU_MAX_RANKS" -ge 2 ]; then
        GPU_NPX=2; GPU_NPY=1
    else
        GPU_NPX=1; GPU_NPY=1
    fi
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
# The run directory is removed on success and kept on failure; KEEP_RUN_DIR=1
# in the environment keeps it either way, to inspect the solver logs.
KEEP_RUN_DIR="${KEEP_RUN_DIR:-0}"
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

# run_udales_gpu <run_dir> <nprocs> -> like run_udales, on the GPU executable
run_udales_gpu() {
    local run_dir="$1" np="$2"
    local module_setup="" selftest="" wrapper=""
    if [ -f /etc/profile.d/modules.sh ]; then
        module_setup="source /etc/profile.d/modules.sh >/dev/null 2>&1 || true; "
    fi
    if [ "$UDALES_GPU_SELFTEST" = 1 ]; then
        selftest="export UDALES_RUN_CUDA_SELFTEST=1; "
    fi
    if [ -n "$UDALES_GPU_LAUNCH_WRAPPER" ]; then
        wrapper="'${UDALES_GPU_LAUNCH_WRAPPER}' "
    fi
    bash -c "${module_setup}\
if command -v module >/dev/null 2>&1 && [ -n \"\$UDALES_GPU_RUNTIME_MODULES\" ]; then module load \$UDALES_GPU_RUNTIME_MODULES >/dev/null 2>&1; fi; \
export HDF5_USE_FILE_LOCKING=FALSE; ${selftest}\
cd '${run_dir}' && '${UDALES_GPU_MPIEXEC}' ${GPU_MPI_LAUNCH_EXTRA_ARGS} -n ${np} ${wrapper}'${UDALES_GPU_BUILD}' '${NAMELIST}' > run.log 2>&1"
}
export UDALES_GPU_RUNTIME_MODULES

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

# run_chain <tag> <nprocx> <nprocy> [trestart] [runner]
# Reference cold start followed by a warm-started continuation. The checkpoint
# time defaults to T_RESTART; the restarted run always runs on to T_TOTAL.
# runner is run_udales (default) or run_udales_gpu.
run_chain() {
    local tag="$1" npx="$2" npy="$3" trestart="${4:-$T_RESTART}" runner="${5:-run_udales}"
    local np=$((npx * npy))
    local rst_runtime
    rst_runtime="$(awk -v a="$T_TOTAL" -v b="$trestart" 'BEGIN{printf "%g", a-b}')"

    local ref_dir
    ref_dir="$(stage "ref_${tag}" "$npx" "$npy" 1 "$T_TOTAL" "$trestart" .false. '')"
    echo "--- reference chain [${tag}] (${npx}x${npy}), runtime=${T_TOTAL}, trestart=${trestart}"
    if ! "$runner" "$ref_dir" "$np"; then
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
    rst_dir="$(stage "rst_${tag}" "$npx" "$npy" 2 "$rst_runtime" "$T_TOTAL" .true. \
        "${stem}_xxx_xxx.${CASE_ID}")"
    cp "${ref_dir}/${stem}"_*_*."${CASE_ID}" "$rst_dir"/
    cp "${ref_dir}/${stem/initd/initf}.${CASE_ID}" "$rst_dir"/ 2>/dev/null || true
    echo "--- restart chain [${tag}], runtime=${rst_runtime}"
    if ! "$runner" "$rst_dir" "$np"; then
        fail "restart run [${tag}] did not complete"
        tail -n 60 "${rst_dir}/run.log" || true
        return 1
    fi
    if ! grep -q "Restored facet state from ${stem/initd/initf}" "${rst_dir}/run.log"; then
        fail "restart run [${tag}] did not report restoring the facet state"
        grep -i "facet" "${rst_dir}/run.log" || true
    fi
    if ! grep -q "Energy-balance interval restored from restart" "${rst_dir}/run.log"; then
        fail "restart run [${tag}] did not report restoring the energy-balance interval state"
        grep -i "energy-balance" "${rst_dir}/run.log" || true
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
# checkpoint inside an energy-balance interval (see T_RESTART_MID above)
run_chain mid_2x2 2 2 "$T_RESTART_MID" || true

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
    # No fallback: with lreadfacT (the default) readfacetfiles does not read
    # Tfacinit*.inp, so a restart directory without the initf file has no facet
    # state at all and the run must abort with an explicit message.
    # -----------------------------------------------------------------------
    missing_dir="$(stage missing_2x2 2 2 2 1. "$T_TOTAL" .true. \
        "${stem}_xxx_xxx.${CASE_ID}")"
    cp "${ref_dir}/${stem}"_*_*."${CASE_ID}" "$missing_dir"/
    rm -f "${missing_dir}/initf"*
    echo "--- missing initf file: the warm start must abort"
    if run_udales "$missing_dir" 4; then
        fail "restart without a facet restart file completed instead of aborting"
    else
        if ! grep -q "ERROR: facet restart file .* not found" "${missing_dir}/run.log"; then
            fail "restart without a facet restart file aborted without the expected error message"
            tail -n 40 "${missing_dir}/run.log" || true
        fi
        if ! grep -q "lreadfacT = .false." "${missing_dir}/run.log"; then
            fail "the missing-file error does not point at lreadfacT = .false."
        fi
    fi
    if grep -q "Restored facet state from" "${missing_dir}/run.log"; then
        fail "restart without a facet restart file claimed to restore facet state"
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

# ---------------------------------------------------------------------------
# Mid-interval checkpoint. The first facT record after the restart must be at
# the reference's next EB time (the restored tnextEB), not dtEB after the
# restart, and must match the reference there: the flux accumulated between the
# last EB update and the checkpoint has been carried across.
# ---------------------------------------------------------------------------
if [ -f "${ROOT_RUN_DIR}/chain_mid_2x2.paths" ]; then
    module_setup=""
    if [ -f /etc/profile.d/modules.sh ]; then
        module_setup="source /etc/profile.d/modules.sh >/dev/null 2>&1 || true; "
    fi
    refmid="$(sed -n 1p "${ROOT_RUN_DIR}/chain_mid_2x2.paths")"
    rstmid="$(sed -n 2p "${ROOT_RUN_DIR}/chain_mid_2x2.paths")"
    echo "--- mid-interval checkpoint at t = ${T_RESTART_MID} (dtEB = ${DTEB})"
    export COMPARE_ARGS="--reference ${refmid}/facT.${CASE_ID}.nc --restart ${rstmid}/facT.${CASE_ID}.nc"
    export UDALES_PYTHON_MODULES PYTHON
    bash -c "${module_setup}\
if command -v module >/dev/null 2>&1 && [ -n \"\$UDALES_PYTHON_MODULES\" ]; then module load \$UDALES_PYTHON_MODULES >/dev/null 2>&1; fi; \
\"\$PYTHON\" '${SCRIPT_DIR}/compare_facT.py' \$COMPARE_ARGS" || fail "mid-interval facT comparison"
else
    fail "the mid-interval reference/restart chain did not produce outputs to compare"
fi

# ---------------------------------------------------------------------------
# GPU build (optional, see the top of this file).
# ---------------------------------------------------------------------------
if [ -n "$UDALES_GPU_BUILD" ]; then
    gtag="${GPU_NPX}x${GPU_NPY}"
    gnp=$((GPU_NPX * GPU_NPY))
    echo "=========================================="
    echo "GPU build     : $UDALES_GPU_BUILD"
    echo "GPU launch    : $UDALES_GPU_MPIEXEC $GPU_MPI_LAUNCH_EXTRA_ARGS -n $gnp ${UDALES_GPU_LAUNCH_WRAPPER:-<no wrapper>} (${gtag})"
    echo "GPU self-tests: $UDALES_GPU_SELFTEST"
    echo "=========================================="

    # The CPU chain of the same decomposition is the parity reference; 2x2 and
    # 2x1 exist from above, anything else is run here.
    if [ ! -f "${ROOT_RUN_DIR}/chain_${gtag}.paths" ]; then
        run_chain "$gtag" "$GPU_NPX" "$GPU_NPY" || true
    fi
    if [ ! -f "${ROOT_RUN_DIR}/chain_mid_${gtag}.paths" ]; then
        run_chain "mid_${gtag}" "$GPU_NPX" "$GPU_NPY" "$T_RESTART_MID" || true
    fi

    # (a) GPU warm start from the CPU checkpoint: the facet state the host
    #     restores from initf must be what the device computes with.
    gpu_from_cpu_dir=""
    if [ -f "${ROOT_RUN_DIR}/chain_${gtag}.paths" ]; then
        cpu_ref="$(sed -n 1p "${ROOT_RUN_DIR}/chain_${gtag}.paths")"
        cpu_stem="$(sed -n 3p "${ROOT_RUN_DIR}/chain_${gtag}.paths")"
        rst_runtime="$(awk -v a="$T_TOTAL" -v b="$T_RESTART" 'BEGIN{printf "%g", a-b}')"
        gpu_from_cpu_dir="$(stage "gpu_from_cpu_${gtag}" "$GPU_NPX" "$GPU_NPY" 2 "$rst_runtime" "$T_TOTAL" .true. \
            "${cpu_stem}_xxx_xxx.${CASE_ID}")"
        cp "${cpu_ref}/${cpu_stem}"_*_*."${CASE_ID}" "$gpu_from_cpu_dir"/
        cp "${cpu_ref}/${cpu_stem/initd/initf}.${CASE_ID}" "$gpu_from_cpu_dir"/
        echo "--- GPU warm start [${gtag}] from the CPU checkpoint ${cpu_stem}"
        if ! run_udales_gpu "$gpu_from_cpu_dir" "$gnp"; then
            fail "GPU warm start from the CPU checkpoint [${gtag}] did not complete"
            tail -n 60 "${gpu_from_cpu_dir}/run.log" || true
            gpu_from_cpu_dir=""
        elif ! grep -q "Restored facet state from" "${gpu_from_cpu_dir}/run.log"; then
            fail "GPU warm start [${gtag}] did not report restoring the facet state"
        fi
    else
        fail "no CPU ${gtag} chain to warm-start the GPU build from"
    fi

    # (b) GPU reference and GPU restart, checkpoint on an energy-balance update.
    # (c) The same with the checkpoint inside an interval: the device-side flux
    #     accumulators have to reach the initf file and be picked up again.
    run_chain "gpu_${gtag}" "$GPU_NPX" "$GPU_NPY" "$T_RESTART" run_udales_gpu || true
    run_chain "gpu_mid_${gtag}" "$GPU_NPX" "$GPU_NPY" "$T_RESTART_MID" run_udales_gpu || true

    module_setup=""
    if [ -f /etc/profile.d/modules.sh ]; then
        module_setup="source /etc/profile.d/modules.sh >/dev/null 2>&1 || true; "
    fi
    export UDALES_PYTHON_MODULES PYTHON
    compare() {
        export COMPARE_ARGS="$*"
        bash -c "${module_setup}\
if command -v module >/dev/null 2>&1 && [ -n \"\$UDALES_PYTHON_MODULES\" ]; then module load \$UDALES_PYTHON_MODULES >/dev/null 2>&1; fi; \
\"\$PYTHON\" '${SCRIPT_DIR}/compare_facT.py' \$COMPARE_ARGS"
    }

    if [ -n "$gpu_from_cpu_dir" ] && [ -f "${ROOT_RUN_DIR}/chain_${gtag}.paths" ]; then
        cpu_ref="$(sed -n 1p "${ROOT_RUN_DIR}/chain_${gtag}.paths")"
        echo "--- GPU warm start from the CPU checkpoint vs the CPU reference (${gtag})"
        compare --reference "${cpu_ref}/facT.${CASE_ID}.nc" --restart "${gpu_from_cpu_dir}/facT.${CASE_ID}.nc" \
                --cross-build || fail "GPU-from-CPU-checkpoint facT comparison"
    fi
    if [ -f "${ROOT_RUN_DIR}/chain_gpu_${gtag}.paths" ]; then
        gref="$(sed -n 1p "${ROOT_RUN_DIR}/chain_gpu_${gtag}.paths")"
        grst="$(sed -n 2p "${ROOT_RUN_DIR}/chain_gpu_${gtag}.paths")"
        other=()
        if [ -f "${ROOT_RUN_DIR}/chain_${gtag}.paths" ]; then
            other=(--other-build "$(sed -n 1p "${ROOT_RUN_DIR}/chain_${gtag}.paths")/facT.${CASE_ID}.nc")
        fi
        echo "--- GPU restart vs GPU reference, and GPU reference vs CPU reference (${gtag})"
        compare --reference "${gref}/facT.${CASE_ID}.nc" --restart "${grst}/facT.${CASE_ID}.nc" \
                ${other[@]+"${other[@]}"} || fail "GPU facT comparison"
    else
        fail "the GPU reference/restart chain did not produce outputs to compare"
    fi
    if [ -f "${ROOT_RUN_DIR}/chain_gpu_mid_${gtag}.paths" ]; then
        gref="$(sed -n 1p "${ROOT_RUN_DIR}/chain_gpu_mid_${gtag}.paths")"
        grst="$(sed -n 2p "${ROOT_RUN_DIR}/chain_gpu_mid_${gtag}.paths")"
        other=()
        if [ -f "${ROOT_RUN_DIR}/chain_mid_${gtag}.paths" ]; then
            other=(--other-build "$(sed -n 1p "${ROOT_RUN_DIR}/chain_mid_${gtag}.paths")/facT.${CASE_ID}.nc")
        fi
        echo "--- GPU mid-interval checkpoint at t = ${T_RESTART_MID} (${gtag})"
        compare --reference "${gref}/facT.${CASE_ID}.nc" --restart "${grst}/facT.${CASE_ID}.nc" \
                ${other[@]+"${other[@]}"} || fail "GPU mid-interval facT comparison"
    else
        fail "the GPU mid-interval reference/restart chain did not produce outputs to compare"
    fi
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
