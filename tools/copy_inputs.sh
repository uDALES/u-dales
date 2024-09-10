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

# examples:
# new sim: "copy_inputs.sh 1 2" or "copy_inputs.sh 1 2 c"
# continue with new sim  "copy_inputs.sh 1 2 w"

set -e

if (( $# < 2 )) ; then
  echo "usage: `basename $0` sim#1_local_path sim#2 (start)"
  echo "prepares a new simulation based on an old one with uDALES"
  echo "   sim#1: local path of source case upon which the new one will be based"
  echo "   sim#2: three digit integer case number of the new simulation"
  echo "   start (optional): (c)old- or (w)arm-start, c is default"
  echo "... execution terminated"
  exit 1
fi

# go to experiment directory
pushd $1
	inputdir=$(pwd)

	## set source experiment number via path of the source directory
	src="${inputdir: -3}"

	## read in additional variables
	if [ -f config.sh ]; then
   	 source config.sh
	else
	 echo "config.sh must be set inside $inputdir"
     exit 1
	fi

	## check if required variables are set
	if [ -z $DA_EXPDIR ]; then
		echo "Experiment directory DA_EXPDIR must be set $inputdir/config.sh"
		exit 1
	fi;
  if [ -z $DA_WORKDIR ]; then
		echo "Experiment directory DA_WORKDIR must be set $inputdir/config.sh"
		exit 1
	fi;

popd

tar=$(printf "%03.0f" $2)    # pad argument 2 (target simulation) with zeros
start=${3:-c}                # default argument 3 to c (coldstart) if not provided
case=1                       # default value 1 (target folder does not exist or overwrite, 2=only copy files that don't exist, 3=interactive)

# check if target simulation already exists. If so, ask how to proceed.
if [ -d $DA_EXPDIR/$tar ]; then
  echo "target directory $tar already exists in $DA_EXPDIR"
  echo "continue with overwrite/copy nonexistent files/interactive/abort? (o/c/i/a)"
  read answer
  if [ "$answer" == "o" ]; then
    case=1
  elif [ "$answer" == "c" ]; then
    case=2
  elif [ "$answer" == "i" ]; then
    case=3
  else
    echo "abort"
    exit 1
  fi
else
  mkdir -p $DA_EXPDIR/$tar
fi

if [ -d $DA_WORKDIR/$tar ]; then
  echo "target directory $tar already exists in $DA_WORKDIR"
else
  mkdir -p $DA_WORKDIR/$tar
  echo "Creating target work directory, $tar, in $DA_WORKDIR"
fi

# test if original simulation exists
if [ ! -d $DA_EXPDIR/$src ]; then
  echo "original simulation $src does not exist in $DA_EXPDIR"
  echo "exit"
  exit 1
fi

# list of files to copy
declare -a tocopy=("/namoptions." "/lscale.inp." "/prof.inp." "/scalar.inp." "/purifs.inp." "/trees.inp." "/scals.inp." "/lad.inp." "/facetarea.inp." "/facets.inp." "/facets_unused." "/netsw.inp." "/svf.inp." "/Tfacinit.inp." "/vf.nc.inp." "/factypes.inp." "/timedepsw.inp." "/timedeplw.inp." "/timedepnudge.inp." "/vfsparse.inp.")

# copy and rename files
case $case in
  1)  #target folder didn't exist or just overwrite all files
    for i in "${tocopy[@]}"
      do
        cp $DA_EXPDIR/$src$i$src $DA_EXPDIR/$tar$i$tar 2>/dev/null || :
        #echo $i
      done
  ;;
  2)  #don't overwrite files in target folder #cp -n does not overwrite an existing file
    for i in "${tocopy[@]}"
      do
        cp -n $DA_EXPDIR/$src$i$src $DA_EXPDIR/$tar$i$tar 2>/dev/null || :
        #echo $i
      done
  ;;
  3)  #ask for every file if to overwrite or not
    for i in "${tocopy[@]}"
      do
        cp -i $DA_EXPDIR/$src$i$src $DA_EXPDIR/$tar$i$tar 2>/dev/null || :
        #echo $i
      done
  ;;
esac


####### modified function for sed -i
function sedi { if [[ "$OSTYPE" == "darwin"* ]]; then
        		sed -i '' "$1" "$2"
		elif [[ "$OSTYPE" == "linux"* ]]; then
        		sed -i "$1" "$2"
		fi;}

function update_namoptions() {
    local varname=$1
    local value=$2
    local sectionname=$3
    if grep -w -q "$varname" "$DA_EXPDIR/$tar/namoptions.$tar"; then
        sedi "/$varname /s/.*/$varname       = $value/g" "$DA_EXPDIR/$tar/namoptions.$tar"
    else
        sedi "/$sectionname/a\\
$varname = $value" "$DA_EXPDIR/$tar/namoptions.$tar"
    fi
}

# change simulation number in namoptions
update_namoptions 'iexpnr' $tar '&RUN'

# coldstart: set warmstart to false in namoptions
if [ $start == "c" ]; then
  sed -i.bak -e '/lwarmstart/s/.*/lwarmstart =  .false./g' $DA_EXPDIR/$tar"/namoptions."$tar
  rm $DA_EXPDIR/$tar"/namoptions."$tar".bak"

# warmstart: link newest warmstart files
elif [ $start == "w" ]; then
  # momentum warmstart files (minimum requirement for warmstarts)
  startfilen=$(ls -t $DA_WORKDIR/$src"/initd"* | head -1)
  if [ -z "$startfilen" ]; then
    echo "no restart files found in $DA_WORKDIR/$src"
    echo "exit"
    exit 1
  else
    # create links
    startfilen=${startfilen##*/}  # retain the part after the last slash
    startfilen=${startfilen%_*_*}   # retain the part before the underscore
    ln -s $DA_WORKDIR/$src"/"*$startfilen* $DA_EXPDIR/$tar/
  fi

  update_namoptions 'startfile' $startfilen\_xxx_xxx.$tar '&RUN'
  update_namoptions 'lwarmstart' '.true.' '&RUN'

  # scalar warmstart files
  scalarwarmstart=$DA_WORKDIR/$src"/inits"*
  if compgen -G "$scalarwarmstart" > /dev/null; then

    scalarfilen=$(ls -t $DA_WORKDIR/$src"/inits"* | head -1)
      
    # create links
    scalarfilen=${scalarfilen##*/}  # retain the part after the last slash
    scalarfilen=${scalarfilen%_*_*}   # retain the part before the underscore
    ln -s $DA_WORKDIR/$src"/"*$scalarfilen* $DA_EXPDIR/$tar/

    update_namoptions 'lreadscal' '.true.' '&SCALARS'

    else
      echo "Info: no scalar restart files found in $DA_WORKDIR/$src."
  fi

  # rename links
  for f in $DA_EXPDIR/$tar/*.$src; do
    mv $f "${f%.$src}.$tar"
  done
  echo "Creating links to warmstart files in $DA_WORKDIR/$src."

fi

# link any available driver files (without renaming)
driverfiles=$DA_WORKDIR/$src"/"?"driver_"*
for dfile in $driverfiles ; do
    if [ -f $dfile ]; then
      ln -s $dfile $DA_EXPDIR/$tar
    fi
  done
if [ -f $dfile ]; then
    echo "Added links to driverfiles."
fi

# copy config script for execution
cp $DA_EXPDIR/$src/config.sh $DA_EXPDIR/$tar

# copy STL for execution
stlfile=$(find $DA_EXPDIR/$src -iname "*.stl")
if [ -f "$stlfile" ] ; then
  cp "$stlfile" $DA_EXPDIR/$tar
fi

# copy *.txt files for IBM implementation
function copy_files() {
    local filename=$1
    filepath=$DA_EXPDIR/$src/$filename*".txt"
    for file in $filepath ; do
      if [ -f $file ]; then
        cp $file $DA_EXPDIR/$tar
      fi
    done
}
copy_files "solid_"
copy_files "facet_sections_"
copy_files "fluid_boundary_"
copy_files "info"

# copy any available scalar source files
scalarsourcefiles=$DA_EXPDIR/$src/"scalarsource"*
for sfile in $scalarsourcefiles ; do
  if [ -f $sfile ]; then
    sfilenew="${sfile%.$src}"
    sfilename="${sfilenew##*$src/}"
    cp $sfile $DA_EXPDIR/$tar/$sfilename"."$tar
  fi
done
