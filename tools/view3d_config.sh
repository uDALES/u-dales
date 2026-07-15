# View3D runtime configuration.
#
# This file is loaded automatically by the Python View3D launcher immediately
# before the external view3d executable is started. tools/write_inputs.sh also
# sources it while preparing local or PBS preprocessing runs.
#
# Usual workflow:
#   1. Edit values in this file.
#   2. Re-run write_inputs; no extra exports are needed.
#
# For experiment-specific overrides, set VIEW3D_CONFIG in the experiment
# config.sh to point to another file with the same syntax.
# export VIEW3D_CONFIG="/absolute/path/to/my_view3d_config.sh"

# Optional path to a specific View3D executable. Leave unset to use the normal
# build locations under tools/preprocessing/build or tools/View3D/build.
# export VIEW3D_EXE="/path/to/view3d"

# Number of OpenMP threads used by View3D. This is View3D-specific and takes
# priority over OMP_NUM_THREADS inside the modified View3D code.
#
# By default, use the preprocessing CPU request from tools/write_inputs.sh
# (PREPROC_NCPU). Set VIEW3D_NUM_THREADS explicitly here if you want View3D to
# use a different number of threads from the preprocessing job's requested CPU
# count.
if [ -z "${VIEW3D_NUM_THREADS:-}" ]; then
	if [ -n "${PREPROC_NCPU:-}" ]; then
		export VIEW3D_NUM_THREADS="$PREPROC_NCPU"
	fi
fi
# export VIEW3D_NUM_THREADS=128

# General OpenMP thread count. Keep this aligned with VIEW3D_NUM_THREADS unless
# you deliberately need a different OpenMP default for debugging.
if [ -n "${VIEW3D_NUM_THREADS:-}" ]; then
	export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$VIEW3D_NUM_THREADS}"
fi
# export OMP_NUM_THREADS=128

# Dense-matrix safety guard for the legacy dense View3D path, in GiB. The direct
# sparse path normally avoids this dense allocation, but keeping the guard active
# protects fallback/diagnostic runs.
#
# On Imperial HPC preprocessing jobs, derive the guard from PREPROC_MEM while
# leaving 16 GiB for Python, View3D metadata, OpenMP worker state, and system
# overhead. For example, PREPROC_MEM=128gb gives
# VIEW3D_MAX_DENSE_MATRIX_GIB=112. Set VIEW3D_MAX_DENSE_MATRIX_GIB explicitly
# before sourcing this file if you need a different limit.
_view3d_set_dense_limit_from_preproc_mem() {
	_mem="${PREPROC_MEM//[[:space:]]/}"
	case "$_mem" in
		*[gG][iI][bB])
			_mem_gib="${_mem%???}"
			;;
		*[gG][bB])
			_mem_gib="${_mem%??}"
			;;
		*[gG])
			_mem_gib="${_mem%?}"
			;;
		*[0-9])
			_mem_gib="$_mem"
			;;
		*)
			return 1
			;;
	esac
	case "$_mem_gib" in
		""|*[!0-9]*)
			return 1
			;;
	esac
	if [ "$_mem_gib" -gt 16 ]; then
		export VIEW3D_MAX_DENSE_MATRIX_GIB="$((_mem_gib - 16))"
	else
		export VIEW3D_MAX_DENSE_MATRIX_GIB=1
	fi
	return 0
}

if [ -z "${VIEW3D_MAX_DENSE_MATRIX_GIB:-}" ]; then
	if [ -n "${PREPROC_MEM:-}" ] && _view3d_set_dense_limit_from_preproc_mem; then
		:
	else
		export VIEW3D_MAX_DENSE_MATRIX_GIB=112
	fi
fi
unset -f _view3d_set_dense_limit_from_preproc_mem
unset _mem _mem_gib

# Debug/comparison controls. Leave these unset for normal optimized runs.
#
# Force serial View3D execution even when OpenMP is available.
# export VIEW3D_DISABLE_OPENMP=1
#
# Force the old dense output path instead of the direct sparse path.
# export VIEW3D_DISABLE_SPARSE_DIRECT=1
#
# Disable the dense-memory guard. This can allow legacy dense allocation to hit
# the scheduler memory limit, so use it only for controlled comparisons.
# export VIEW3D_DISABLE_DENSE_MEMORY_GUARD=1
