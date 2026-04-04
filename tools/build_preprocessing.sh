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
target=${2:-view3d}

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
fi
cmake "${cmake_args[@]}"
cmake --build "${build_dir}" --target "${target}"

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
fi
if compgen -G "tools/python/udprep/ibm_preproc_f2py*.so" >/dev/null; then
    echo "IBM preprocessing f2py module available at tools/python/udprep/"
fi
