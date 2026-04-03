#!/usr/bin/env bash
set -euo pipefail

delete=0
preprocess_only=0
runtime_only=0
expnr=""

usage() {
  cat <<'EOF'
Usage: tools/cleansim.sh <case_dir> [--expnr NNN] [--delete] [--preprocess-only | --runtime-only]

Dry-run by default. Pass --delete to actually remove files.
EOF
}

if [[ $# -lt 1 ]]; then
  usage
  exit 1
fi

case_dir=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --delete)
      delete=1
      shift
      ;;
    --preprocess-only)
      preprocess_only=1
      shift
      ;;
    --runtime-only)
      runtime_only=1
      shift
      ;;
    --expnr)
      expnr="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      if [[ -z "$case_dir" ]]; then
        case_dir="$1"
      else
        echo "Unexpected argument: $1" >&2
        exit 1
      fi
      shift
      ;;
  esac
done

if [[ $preprocess_only -eq 1 && $runtime_only -eq 1 ]]; then
  echo "--preprocess-only and --runtime-only are mutually exclusive" >&2
  exit 1
fi

if [[ -z "$case_dir" ]]; then
  usage
  exit 1
fi

case_dir="$(cd "$case_dir" && pwd)"

if [[ -z "$expnr" ]]; then
  shopt -s nullglob
  namoptions=( "$case_dir"/namoptions.* )
  shopt -u nullglob
  numeric=()
  for path in "${namoptions[@]}"; do
    name="${path##*.}"
    if [[ "$name" =~ ^[0-9]+$ ]]; then
      numeric+=( "$name" )
    fi
  done
  if [[ ${#numeric[@]} -eq 1 ]]; then
    expnr="${numeric[0]}"
  fi
fi

preprocess_patterns=(
  "xgrid.inp.${expnr}"
  "zgrid.inp.${expnr}"
  "lscale.inp.${expnr}"
  "prof.inp.${expnr}"
  "scalar.inp.*.${expnr}"
  "scalarsourcep.inp.*.${expnr}"
  "scalarsourcel.inp.*.${expnr}"
  "trees.inp.${expnr}"
  "vegetation.inp.${expnr}"
  "veg.inp.${expnr}"
  "veg_params.inp.${expnr}"
  "factypes.inp.${expnr}"
  "facets.inp.${expnr}"
  "facetarea.inp.${expnr}"
  "solid_*.txt"
  "fluid_boundary_*.txt"
  "facet_sections_*.txt"
  "facets_unused.${expnr}"
  "svf.inp.${expnr}"
  "vfsparse.inp.${expnr}"
  "vf.txt"
  "vf.bin"
  "vf.nc.inp.${expnr}"
  "netsw.inp.${expnr}"
  "Tfacinit.inp.${expnr}"
  "Tfacinit_layers.inp.${expnr}"
  "timedepsw.inp.${expnr}"
  "timedepsveg.inp.${expnr}"
  "Sdir.txt"
  "Sdir.nc"
  "facets.vs3"
  "vertices.txt"
  "faces.txt"
  "info_directShortwave.txt"
  "DS.exe"
  "inmypoly_inp_info.txt"
  "zhgrid.txt"
  "zfgrid.txt"
  "info_fort.txt"
  "info_matchFacetsToCells.txt"
  "IBM_preproc.exe"
  "solid_points.fig"
  "fluid_boundary_points.fig"
)

runtime_patterns=(
  "u-dales"
  "output.${expnr}"
  "job.${expnr}.slurm"
  "post-job.${expnr}"
  "decomp_2d_setup.log"
  "slurm-*.out"
  "monitor*.txt"
  "fielddump*"
  "xytdump*"
  "tdump*"
  "mintdump*"
  "treedump*"
)

patterns=()
if [[ $runtime_only -eq 0 ]]; then
  patterns+=( "${preprocess_patterns[@]}" )
fi
if [[ $preprocess_only -eq 0 ]]; then
  patterns+=( "${runtime_patterns[@]}" )
fi

declare -A seen=()
matches=()
shopt -s nullglob
for pattern in "${patterns[@]}"; do
  for path in "$case_dir"/$pattern; do
    name="$(basename "$path")"
    if [[ -n "${seen[$name]:-}" ]]; then
      continue
    fi
    seen[$name]=1
    matches+=( "$path" )
  done
done
shopt -u nullglob

if [[ ${#matches[@]} -eq 0 ]]; then
  echo "No generated files matched in $case_dir"
  exit 0
fi

if [[ $delete -eq 1 ]]; then
  echo "Removing ${#matches[@]} path(s) from $case_dir"
else
  echo "Would remove ${#matches[@]} path(s) from $case_dir"
fi

printf '%s\n' "${matches[@]##*/}" | sort | sed 's/^/  - /'

if [[ $delete -eq 0 ]]; then
  echo "Dry run only. Re-run with --delete to remove these files."
  exit 0
fi

for path in "${matches[@]}"; do
  rm -rf "$path"
done
