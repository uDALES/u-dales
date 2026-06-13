#!/usr/bin/env bash
# Build preprocessing tools for uDALES.
# Compiles the View3D executable and optionally the directshortwave and IBM
# f2py Python extension modules through the standalone CMake entry point at
# tools/preprocessing/.
#
# Usage (run from the repository root):
#   tools/build_preprocessing.sh <build_system> [build_target]
#
# Arguments:
#   build_system   Build environment to use (required).
#                  Allowed values : common | icl
#                    common - local Linux / WSL system (no module loading)
#                    icl    - Imperial College London HPC cluster
#                             (loads CMake/3.31.3-GCCcore-14.2.0 and
#                              Python/3.13.1-GCCcore-14.2.0 modules)
#
#   build_target   CMake target to build.
#                  Allowed values : view3d | preprocessing_tools
#                  Default        : view3d
#                    view3d               - View3D executable only
#                    preprocessing_tools  - View3D + directshortwave and
#                                          IBM f2py extension modules
#                                          (requires numpy and Python headers)
#
# Environment variables (optional overrides):
#   PREPROCESSING_PYTHON_EXECUTABLE
#                  Python interpreter used for the CMake build and f2py
#                  compilation. Use this when dependencies are provided by a
#                  managed environment, such as GitHub Actions' conda setup.
#                  If unset, preprocessing_tools requires tools/python/.venv.
#
# Output:
#   tools/preprocessing/build/bin/view3d
#     View3D executable (also symlinked to tools/View3D/build/src/view3d for
#     MATLAB compatibility).
#   tools/python/udprep/directshortwave_f2py*.so
#     Direct shortwave f2py module (preprocessing_tools target only).
#   tools/python/udprep/ibm_preproc_f2py*.so
#     IBM preprocessing f2py module (preprocessing_tools target only).
#
# Examples:
#   # Build View3D only on a local system
#   tools/build_preprocessing.sh common
#
#   # Build View3D + f2py modules on a local system
#   tools/build_preprocessing.sh common preprocessing_tools
#
#   # Full build on the ICL cluster
#   tools/build_preprocessing.sh icl preprocessing_tools
#
set -euo pipefail

# Usage: ./tools/build_preprocessing.sh [common / icl] [target]
if (( $# < 1 ))
then
 echo "The build type <common / icl> must be set."
 echo "usage: from being in u-dales directory run: tools/build_preprocessing.sh <build type> [target]"
 echo "default target: view3d"
 echo "available targets: view3d, preprocessing_tools (view3d + f2py modules)"
 exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UDALES_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${UDALES_ROOT}"

system=$1
target=${2:-view3d}
python_venv_dir="${UDALES_ROOT}/tools/python/.venv"
python_cmd=""

resolve_python_cmd() {
    local requested_python="${PREPROCESSING_PYTHON_EXECUTABLE:-}"

    if [ -n "${requested_python}" ]; then
        if [ -x "${requested_python}" ]; then
            python_cmd="${requested_python}"
            return 0
        fi

        if command -v "${requested_python}" >/dev/null 2>&1; then
            python_cmd="$(command -v "${requested_python}")"
            return 0
        fi

        echo "[ERROR] PREPROCESSING_PYTHON_EXECUTABLE is not executable: ${requested_python}"
        exit 1
    fi

    if [ -x "${python_venv_dir}/bin/python" ]; then
        python_cmd="${python_venv_dir}/bin/python"
        return 0
    fi

    if [ -x "${python_venv_dir}/bin/python3" ]; then
        python_cmd="${python_venv_dir}/bin/python3"
        return 0
    fi

    return 1
}

require_python_venv_for_preprocessing_tools() {
    if [ "$target" != "preprocessing_tools" ]; then
        return
    fi

    if [ -n "${python_cmd}" ]; then
        return
    fi

    echo "[ERROR] Python environment is required to build target 'preprocessing_tools'."
    echo "Expected virtual environment: ${UDALES_ROOT}/tools/python/.venv"
    echo "or explicit interpreter: PREPROCESSING_PYTHON_EXECUTABLE=/path/to/python"
    echo ""
    echo "Set it up from the repository root with:"
    echo "  bash ${UDALES_ROOT}/tools/python/setup_venv.sh $system preprocessing_tools"
    echo ""
    echo "The setup script creates the virtual environment and internally runs:"
    echo "  ${UDALES_ROOT}/tools/build_preprocessing.sh $system preprocessing_tools"
    echo "so you do not need to run build_preprocessing.sh explicitly after setup completes."
    exit 1
}

if [[ "$system" != "common" && "$system" != "icl" ]]; then
    echo "This configuration is not available"
    exit 1
fi

if [[ "$target" != "view3d" && "$target" != "preprocessing_tools" ]]; then
    echo "Invalid preprocessing target: $target"
    echo "available targets: view3d, preprocessing_tools (view3d + f2py modules)"
    exit 1
fi

if [ "$system" == "icl" ]
then
    module load CMake/3.31.3-GCCcore-14.2.0
    module load Python/3.13.1-GCCcore-14.2.0
elif [ "$system" == "common" ]
then
    echo "Building preprocessing target '${target}' on local system."
fi

if ! resolve_python_cmd; then
    python_cmd=""
fi
require_python_venv_for_preprocessing_tools

if [ -n "${python_cmd}" ]; then
    echo "Using preprocessing Python: ${python_cmd}"
    "${python_cmd}" --version 2>&1 || true
fi

if ! command -v cmake >/dev/null 2>&1; then
    echo "cmake not found"
    exit 1
fi

build_dir="tools/preprocessing/build"
cmake_args=(
    -S tools/preprocessing
    -B "${build_dir}"
    -DCMAKE_POLICY_VERSION_MINIMUM=3.5
)
if [ -n "${python_cmd}" ]; then
    cmake_args+=(-DPREPROCESSING_PYTHON_EXECUTABLE="${python_cmd}")
    missing_f2py_deps=0
    if ! "${python_cmd}" -c "import numpy" >/dev/null 2>&1; then
        echo "WARNING: numpy not available for ${python_cmd}; disabling f2py targets."
        missing_f2py_deps=1
    fi
    if ! "${python_cmd}" -c "import pathlib, sysconfig; header = pathlib.Path(sysconfig.get_paths().get('include','')) / 'Python.h'; raise SystemExit(0 if header.is_file() else 1)" >/dev/null 2>&1; then
        echo "WARNING: Python.h not found for ${python_cmd}; disabling f2py targets."
        missing_f2py_deps=1
    fi
    if [ "${missing_f2py_deps}" -eq 0 ] && [ "${target}" = "preprocessing_tools" ]; then
        cmake_args+=(-DBUILD_PREPROCESSING_DIRECTSHORTWAVE_F2PY=ON)
        cmake_args+=(-DBUILD_PREPROCESSING_IBM_F2PY=ON)
    fi
fi
cmake "${cmake_args[@]}"
cmake --build "${build_dir}" --target "${target}"

require_f2py=1

if [ -f "${build_dir}/bin/view3d" ]; then
    legacy_view3d_dir="tools/View3D/build/src"
    mkdir -p "${legacy_view3d_dir}"
    ln -sfn "../../../preprocessing/build/bin/view3d" "${legacy_view3d_dir}/view3d"
    echo "View3D executable available at ${UDALES_ROOT}/${build_dir}/bin/view3d"
    echo "MATLAB compatibility path available at ${UDALES_ROOT}/${legacy_view3d_dir}/view3d"
fi
if compgen -G "tools/python/udprep/directshortwave_f2py*.so" >/dev/null; then
    echo "directshortwave f2py module available at ${UDALES_ROOT}/tools/python/udprep/"
elif [ "${require_f2py}" -eq 1 ]; then
    echo "WARNING: directshortwave f2py module missing in ${UDALES_ROOT}/tools/python/udprep/"
    echo "         This usually means the preprocessing Python lacks numpy or Python headers."
fi
if compgen -G "tools/python/udprep/ibm_preproc_f2py*.so" >/dev/null; then
    echo "IBM preprocessing f2py module available at ${UDALES_ROOT}/tools/python/udprep/"
elif [ "${require_f2py}" -eq 1 ]; then
    echo "WARNING: ibm_preproc f2py module missing in ${UDALES_ROOT}/tools/python/udprep/"
    echo "         This usually means the preprocessing Python lacks numpy or Python headers."
fi
