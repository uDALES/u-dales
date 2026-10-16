#!/usr/bin/env bash

# uDALES (https://github.com/uDALES/u-dales).

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Copyright (C) 2016-2019 the uDALES Team.

# Runs cases in the examples folders.
# Data are saved under outputs.
# Usage: ./tools/examples/run_examples.sh [case ...]
#
# Without arguments, every case in examples/ is run (except those listed in
# EXCLUDED below), in ascending order. This guarantees that the driver
# (precursor) case 949 runs before the driven case 950.
#
# DA_BUILD, DA_WORKDIR and DA_TOOLSDIR default to the paths below and can be
# overridden from the environment. The number of MPI ranks for each case is
# taken from nprocx*nprocy in its namoptions file.

set -e

if [ ! -d src ]; then
  echo "Please run this script from the root folder"
  exit 1
fi

export DA_TOOLSDIR=${DA_TOOLSDIR:-$(pwd)/tools}
export DA_BUILD=${DA_BUILD:-$(pwd)/build/release/u-dales}
export DA_WORKDIR=${DA_WORKDIR:-$(pwd)/outputs}

# Cases that cannot sensibly be run by this script:
#   024: large-scale HPC case (1024^3 grid on nprocx*nprocy = 32*32 = 1024 cores).
EXCLUDED="024"

if (( $# > 0 )); then
    examples="$*"
else
    examples=""
    for dir in examples/*/; do
        example=$(basename "$dir")
        if [[ " $EXCLUDED " == *" $example "* ]]; then
            echo "Skipping excluded case $example"
            continue
        fi
        examples="$examples $example"
    done
fi

# Case inputs are staged in a temporary directory so that the config.sh
# shipped with each example (which contains machine-specific paths) does not
# override the settings above, and so that no files are added to examples/.
stagedir=$(mktemp -d)
trap 'rm -rf "$stagedir"' EXIT

for example in $examples
do
    namoptions=examples/$example/namoptions.$example
    if [ ! -f "$namoptions" ]; then
        echo "Case $example not found: $namoptions does not exist"
        exit 1
    fi

    # Always start from afresh
    rm -rf "${DA_WORKDIR:?}/$example"

    cp -r "examples/$example" "$stagedir/$example"
    rm -f "$stagedir/$example/config.sh"

    # The number of cores must equal nprocx*nprocy set in the namoptions file.
    nprocx=$(awk -F= '/^[[:space:]]*nprocx[[:space:]]*=/ {gsub(/[[:space:]]/,"",$2); print $2; exit}' "$namoptions")
    nprocy=$(awk -F= '/^[[:space:]]*nprocy[[:space:]]*=/ {gsub(/[[:space:]]/,"",$2); print $2; exit}' "$namoptions")
    export NCPU=$(( ${nprocx:-1} * ${nprocy:-1} ))

    if [[ $example == 102 ]]; then
        # Warmstart simulation: the restart files must be next to the other inputs.
        cp "examples/$example"/warmstart_files/init?00000267_*."$example" "$stagedir/$example/"
    fi

    if [[ $example == 950 ]]; then
        # Driven simulation: link the driver files written by the precursor
        # simulation 949 (which must have been run first).
        if [ ! -d "$DA_WORKDIR/949" ]; then
            echo "Case 950 requires the outputs of driver case 949 in $DA_WORKDIR/949; run 949 first"
            exit 1
        fi
        "$DA_TOOLSDIR/link_driver_files.sh" "$DA_WORKDIR/949" "$stagedir/$example"
    fi

    ./tools/local_execute.sh "$stagedir/$example"
done
