#!/usr/bin/env bash
# Compiler and tool modules for the GPU dependency stacks: the one place where
# machine-specific module names live. Sourced by every
# *_s1_set_env_for_install.sh and by stack3-*/openmpi/s1_set_env.sh.
#
#   source modules.sh            NVHPC (on HX1 also CMake and zlib)
#   source modules.sh openmpi    the above plus what the stack 3 OpenMPI links
#                                (on HX1: UCX-CUDA, PMIx, hwloc, libevent)
#
# UDALES_DEP_SYSTEM selects the machine, hx1 (default) or common;
# install_dependencies.sh exports it from its hx1|common argument.
#
# Exports used by the stage scripts:
#   DEP_NVHPC_HOME      versioned SDK root (.../Linux_x86_64/<version>)
#   DEP_CUDA_HOME       CUDA inside the SDK (OpenMPI --with-cuda)
#   DEP_MPI_DEFAULT     the OpenMPI the SDK calls comm_libs/mpi      (stack 1)
#   DEP_MPI_OPENMPI4    the OpenMPI 4 shipped in the SDK             (stack 2)
#   DEP_ZLIB_ROOT DEP_ZLIB_INCDIR DEP_ZLIB_LIBDIR DEP_ZLIB_LIBRARY
#   DEP_TM_ROOT DEP_UCX_ROOT DEP_PMIX_ROOT DEP_HWLOC_ROOT DEP_LIBEVENT_ROOT
#                       for the stack 3 OpenMPI configure; empty means "leave
#                       it out" (OpenMPI then builds its internal copy or does
#                       without)

# ---- common: a Linux machine with the NVIDIA HPC SDK ------------------------
# The SDK's modulefiles directory and the module to load. Change these two
# lines, or export the variables before running install_dependencies.sh, if
# your SDK lives elsewhere or is another version.
: "${NVHPC_MODULEFILES:=/opt/nvidia/hpc_sdk/modulefiles}"
: "${NVHPC_MODULE:=nvhpc/24.11}"

# ---- hx1: Imperial HX1 (EasyBuild modules) ----------------------------------
HX1_NVHPC_MODULE=NVHPC/23.7-CUDA-12.2.0
HX1_NVHPC_CUDA=12.2
HX1_TOOL_MODULES="CMake/3.26.3-GCCcore-12.3.0 zlib/1.2.13-GCCcore-12.3.0"
HX1_OPENMPI_MODULES="UCX-CUDA/1.14.1-GCCcore-12.3.0-CUDA-12.1.1 PMIx/4.2.4-GCCcore-12.3.0 hwloc/2.9.1-GCCcore-12.3.0 libevent/2.1.12-GCCcore-12.3.0"

_dep_self=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/$(basename -- "${BASH_SOURCE[0]}")
_dep_line() { command grep -n -F -m1 -- "$1" "$_dep_self" | cut -d: -f1; }

# $1: the module(s) that failed to load. Prints what exists instead and where
# to change the name, then exits.
_dep_module_missing() {
    local stem=${1%%/*} avail
    {
        echo
        echo "Error: module '$1' could not be loaded."
        avail=$(module -t avail "$stem" 2>&1 | grep -v -e ':$' -e '^$' -e '^No module' || true)
        if [ -n "$avail" ]; then
            echo "Modules with that name on this machine (module avail $stem):"
            echo "$avail" | sed 's/^/    /'
        else
            echo "No module named '$stem' exists on this machine (module avail $stem found nothing)."
        fi
        echo
        case "${UDALES_DEP_SYSTEM:-hx1}" in
        common)
            [ -d "$NVHPC_MODULEFILES" ] || echo "The modulefiles directory $NVHPC_MODULEFILES does not exist either."
            if [ -n "$avail" ]; then
                echo "A different version exists: put its name in NVHPC_MODULE at"
                echo "    $_dep_self  line $(_dep_line ': "${NVHPC_MODULE:=')"
                echo "(or export NVHPC_MODULE=<name> before running install_dependencies.sh) and run again."
            else
                echo "Install the NVIDIA HPC SDK from https://developer.nvidia.com/hpc-sdk-downloads"
                echo "(its default location /opt/nvidia/hpc_sdk has the modulefiles in"
                echo "/opt/nvidia/hpc_sdk/modulefiles) and run again. If the SDK is already installed"
                echo "somewhere else, set NVHPC_MODULEFILES and NVHPC_MODULE at"
                echo "    $_dep_self  lines $(_dep_line ': "${NVHPC_MODULEFILES:=') and $(_dep_line ': "${NVHPC_MODULE:=')"
                echo "(or export both before running install_dependencies.sh) and run again."
            fi
            ;;
        *)
            echo "Change the module name in"
            echo "    $_dep_self  lines $(_dep_line 'HX1_NVHPC_MODULE=')-$(_dep_line 'HX1_OPENMPI_MODULES=')"
            echo "to one that exists and run again. If the new NVHPC is another major version,"
            echo "set HX1_NVHPC_CUDA to the CUDA version it ships (ls \$EBROOTNVHPC/Linux_x86_64/*/cuda)."
            ;;
        esac
    } >&2
    exit 1
}

_dep_missing_tool() {   # $1 what, $2 how to get it
    echo "Error: $1 not found. $2" >&2
    exit 1
}

case "${UDALES_DEP_SYSTEM:-hx1}" in
hx1)
    module purge
    module load "$HX1_NVHPC_MODULE" || _dep_module_missing "$HX1_NVHPC_MODULE"
    # shellcheck disable=SC2086
    module load $HX1_TOOL_MODULES || _dep_module_missing "$HX1_TOOL_MODULES"
    if [ "${1:-}" = openmpi ]; then
        # shellcheck disable=SC2086
        module load $HX1_OPENMPI_MODULES || _dep_module_missing "$HX1_OPENMPI_MODULES"
        export DEP_UCX_ROOT=$EBROOTUCX DEP_PMIX_ROOT=$EBROOTPMIX
        export DEP_HWLOC_ROOT=$EBROOTHWLOC DEP_LIBEVENT_ROOT=$EBROOTLIBEVENT
    fi
    export DEP_NVHPC_HOME="$EBROOTNVHPC/Linux_x86_64/$EBVERSIONNVHPC"
    export DEP_CUDA_HOME="$DEP_NVHPC_HOME/cuda/$HX1_NVHPC_CUDA"
    export DEP_ZLIB_ROOT=$EBROOTZLIB
    export DEP_TM_ROOT=/opt/pbs
    ;;
common)
    if [ -d "$NVHPC_MODULEFILES" ]; then
        module use "$NVHPC_MODULEFILES"
    fi
    module load "$NVHPC_MODULE" || _dep_module_missing "$NVHPC_MODULE"
    # NVIDIA's modulefile sets NVHPC_ROOT; an EasyBuild one sets EBROOTNVHPC.
    if [ -n "${NVHPC_ROOT:-}" ]; then
        export DEP_NVHPC_HOME=$NVHPC_ROOT
    elif [ -n "${EBROOTNVHPC:-}" ]; then
        export DEP_NVHPC_HOME="$EBROOTNVHPC/Linux_x86_64/$EBVERSIONNVHPC"
    else
        _dep_missing_tool "the SDK root (module $NVHPC_MODULE set neither NVHPC_ROOT nor EBROOTNVHPC)" \
            "Export NVHPC_ROOT=<...>/Linux_x86_64/<version> and run again."
    fi
    # The SDK's cuda/ directory links to the default CUDA of that release.
    export DEP_CUDA_HOME="$DEP_NVHPC_HOME/cuda"
    command -v cmake >/dev/null 2>&1 || _dep_missing_tool "cmake" \
        "HDF5 is built with CMake (3.18 or newer): install it with the package manager and run again."
    # zlib: the system development package.
    _inc=""; for _d in /usr/include /usr/local/include; do [ -f "$_d/zlib.h" ] && { _inc=$_d; break; }; done
    _lib=$(ldconfig -p 2>/dev/null | awk '/libz\.so$/ {print $NF; exit}')
    [ -n "$_lib" ] || _lib=$(cc -print-file-name=libz.so 2>/dev/null || true)
    [ -n "$_inc" ] && [ -n "$_lib" ] && [ "${_lib#/}" != "$_lib" ] && [ -e "$_lib" ] || _dep_missing_tool "zlib headers and libz.so" \
        "Install the zlib development package (zlib1g-dev on Debian/Ubuntu, zlib-devel on RHEL/Fedora) and run again."
    export DEP_ZLIB_ROOT=$(dirname "$_inc") DEP_ZLIB_INCDIR=$_inc DEP_ZLIB_LIBDIR=$(dirname "$_lib") DEP_ZLIB_LIBRARY=$_lib
    unset _inc _lib _d
    # OpenMPI (stack 3) extras: PBS tm launch only if PBS is installed, UCX only
    # if UCX_ROOT points at one; PMIx, hwloc, libevent come from OpenMPI itself.
    export DEP_TM_ROOT=""; [ -f /opt/pbs/include/tm.h ] && export DEP_TM_ROOT=/opt/pbs
    export DEP_UCX_ROOT="${UCX_ROOT:-}" DEP_PMIX_ROOT="" DEP_HWLOC_ROOT="" DEP_LIBEVENT_ROOT=""
    ;;
*)
    echo "Error: UDALES_DEP_SYSTEM must be hx1 or common (got '${UDALES_DEP_SYSTEM}')" >&2
    exit 1
    ;;
esac

[ -n "${DEP_ZLIB_INCDIR:-}" ] || export DEP_ZLIB_INCDIR="$DEP_ZLIB_ROOT/include"
[ -n "${DEP_ZLIB_LIBDIR:-}" ] || export DEP_ZLIB_LIBDIR="$DEP_ZLIB_ROOT/lib"
[ -n "${DEP_ZLIB_LIBRARY:-}" ] || export DEP_ZLIB_LIBRARY="$DEP_ZLIB_LIBDIR/libz.so"

# The two OpenMPIs inside the SDK (stacks 1 and 2). comm_libs/mpi is the
# SDK's default (OpenMPI 3.1.5 in 23.7, HPC-X in newer releases). The
# OpenMPI 4 is installed once per CUDA version, comm_libs/<cuda>/openmpi4/
# openmpi-4.1.x; take the one for the SDK's default CUDA (what cuda/bin links
# to), else the newest. comm_libs/openmpi4 itself holds only wrapper scripts.
export DEP_MPI_DEFAULT="$DEP_NVHPC_HOME/comm_libs/mpi"
_cuda=$(basename "$(dirname "$(readlink -f "$DEP_NVHPC_HOME/cuda/bin")")")
DEP_MPI_OPENMPI4=""
for _d in "$DEP_NVHPC_HOME/comm_libs/$_cuda"/openmpi4/openmpi-* \
          $(ls -d "$DEP_NVHPC_HOME"/comm_libs/*/openmpi4/openmpi-* 2>/dev/null | sort -V -r); do
    if [ -x "$_d/bin/mpif90" ]; then DEP_MPI_OPENMPI4=$_d; break; fi
done
export DEP_MPI_OPENMPI4
unset _d _cuda
[ -d "$DEP_NVHPC_HOME/compilers/bin" ] || _dep_missing_tool "the NVHPC compilers under $DEP_NVHPC_HOME/compilers/bin" \
    "Check the module and DEP_NVHPC_HOME in $_dep_self."
