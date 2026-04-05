#!/usr/bin/env bash

set -euo pipefail

# Usage: ./tools/build_preprocessing.sh [common / icl] [target]
if (( $# < 1 ))
then
 echo "The build type <common / icl> must be set."
 echo "usage: from being in u-dales directory run: tools/build_preprocessing.sh <build type> [target]"
 echo "default target: view3d"
 exit 1
fi

if [ ! -d tools ]; then
    echo "Please run this script from being inside the u-dales folder"
    exit 1
fi

system=$1
target=${2:-preprocessing_tools}

if [ $system == "icl" ]
then
    module load CMake/3.31.8-GCCcore-14.3.0
elif [ $system == "common" ]
then
    echo "Building preprocessing target '${target}' on local system."
else
    echo "This configuration is not avalable"
    exit 1
fi

python_cmd="${PREPROCESSING_PYTHON_EXECUTABLE:-}"
if [ -z "${python_cmd}" ]; then
    if command -v python >/dev/null 2>&1; then
        python_cmd="$(command -v python)"
    elif command -v python3 >/dev/null 2>&1; then
        python_cmd="$(command -v python3)"
    elif [ "${target}" != "view3d" ]; then
        echo "Neither python nor python3 was found; cannot configure preprocessing build."
        exit 1
    fi
fi

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
    if [ "${missing_f2py_deps}" -eq 1 ]; then
        cmake_args+=(-DBUILD_PREPROCESSING_DIRECTSHORTWAVE_F2PY=OFF)
        cmake_args+=(-DBUILD_PREPROCESSING_IBM_F2PY=OFF)
    fi
fi
cmake "${cmake_args[@]}"
cmake --build "${build_dir}" --target "${target}"

require_f2py=1

if [ -f "${build_dir}/bin/view3d" ]; then
    legacy_view3d_dir="tools/View3D/build/src"
    mkdir -p "${legacy_view3d_dir}"
    ln -sfn "../../../preprocessing/build/bin/view3d" "${legacy_view3d_dir}/view3d"
    echo "View3D executable available at ${build_dir}/bin/view3d"
    echo "MATLAB compatibility path available at ${legacy_view3d_dir}/view3d"
fi
if [ -f "${build_dir}/bin/IBM_preproc" ]; then
    echo "IBM preprocessing executable available at ${build_dir}/bin/IBM_preproc"
fi
if compgen -G "tools/python/udprep/directshortwave_f2py*.so" >/dev/null; then
    echo "directshortwave f2py module available at tools/python/udprep/"
elif [ "${require_f2py}" -eq 1 ]; then
    echo "WARNING: directshortwave f2py module missing in tools/python/udprep/"
    echo "         This usually means the preprocessing Python lacks numpy or Python headers."
fi
if compgen -G "tools/python/udprep/ibm_preproc_f2py*.so" >/dev/null; then
    echo "IBM preprocessing f2py module available at tools/python/udprep/"
elif [ "${require_f2py}" -eq 1 ]; then
    echo "WARNING: ibm_preproc f2py module missing in tools/python/udprep/"
    echo "         This usually means the preprocessing Python lacks numpy or Python headers."
fi
