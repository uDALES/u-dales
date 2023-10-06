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
# new sim: "copy_inputs.sh 2 1" or "copy_inputs.sh 2 1 c"
# continue with old sim: "copy_inputs.sh 1 1 w"
# continue with new sim  "copy_inputs.sh 2 1 w"

set -e

if [ -z $DA_EXPDIR_SRC ]; then
  DA_EXPDIR_SRC=$DA_EXPDIR
fi

if [ -z $DA_WORKDIR_SRC ]; then
  DA_WORKDIR_SRC=$DA_WORKDIR
fi

if (( $# < 2 )) ; then
  echo "usage: `basename $0` sim#1 sim#2 (start)"
  echo "prepares a new simulation with DALES"
  echo "   sim#1: number of new simulation"
  echo "   sim#2: number of simulation upon which the newone is based"
  echo "   start (optional): (c)old- or (w)arm-start, c is default"
  echo "... execution terminated"
  exit 0
fi

tar=$(printf "%03.0f" $1)    # pad argument 1 (target simulation) with zeros
src=$(printf "%03.0f" $2)    # pad argument 2 (origin) with zeros
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
if [ ! -d $DA_EXPDIR_SRC/$src ]; then
  echo "original simulation $src does not exist in $DA_EXPDIR"
  echo "exit"
  exit 1
fi

# list of files to copy
declare -a tocopy=("/namoptions." "/lscale.inp." "/prof.inp." "/scalar.inp." "/purifs.inp." "/trees.inp." "/scals.inp." "/lad.inp." "/facetarea.inp." "/facets.inp." "/netsw.inp." "/svf.inp." "/Tfacinit.inp." "/vf.nc.inp." "/factypes.inp.")

# copy and rename files
case $case in
  1)  #target folder didn't exist or just overwrite all files
    for i in "${tocopy[@]}"
      do
        cp $DA_EXPDIR_SRC/$src$i$src $DA_EXPDIR/$tar$i$tar 2>/dev/null || :
        #echo $i
      done
  ;;
  2)  #don't overwrite files in target folder #cp -n does not overwrite an existing file
    for i in "${tocopy[@]}"
      do
        cp -n $DA_EXPDIR_SRC/$src$i$src $DA_EXPDIR/$tar$i$tar 2>/dev/null || :
        #echo $i
      done
  ;;
  3)  #ask for every file if to overwrite or not
    for i in "${tocopy[@]}"
      do
        cp -i $DA_EXPDIR_SRC/$src$i$src $DA_EXPDIR/$tar$i$tar 2>/dev/null || :
        #echo $i
      done
  ;;
esac

# change simulation number in namoptions
sed -i.bak -e "/iexpnr/s/.*/iexpnr       = $tar/g" $DA_EXPDIR/$tar"/namoptions."$tar
rm $DA_EXPDIR/$tar"/namoptions."$tar".bak"

# coldstart: set warmstart to false in namoptions
if [ $start == "c" ]; then
  sed -i.bak -e '/lwarmstart/s/.*/lwarmstart =  .false./g' $DA_EXPDIR/$tar"/namoptions."$tar
  rm $DA_EXPDIR/$tar"/namoptions."$tar".bak"
# warmstart: link newest warmstart files

elif [ $start == "w" ]; then
  # momentum warmstart files (minimum requirement for warmstarts)
  startfilen=$(ls -t $DA_WORKDIR_SRC/$src"/initd"* | head -1)
  if [ -z "$startfilen" ]; then
    echo "no restart files found in $DA_WORKDIR_SRC/$src"
    echo "exit"
    exit 1
  else
    # create links
    startfilen=${startfilen##*/}  # retain the part after the last slash
    startfilen=${startfilen%_*_*}   # retain the part before the underscore
    ln -s $DA_WORKDIR_SRC/$src"/"*$startfilen* $DA_EXPDIR/$tar/
  fi

  # scalar warmstart files
  scalarfilen=$(ls -t $DA_WORKDIR_SRC/$src"/inits"* | head -1)
  if [ -z "$scalarfilen" ]; then
    echo "Info: no scalar restart files found in $DA_WORKDIR_SRC/$src."
  else
    # create links
    scalarfilen=${scalarfilen##*/}  # retain the part after the last slash
    scalarfilen=${scalarfilen%_*_*}   # retain the part before the underscore
    ln -s $DA_WORKDIR_SRC/$src"/"*$scalarfilen* $DA_EXPDIR/$tar/
  fi

  # rename links
  for f in $DA_EXPDIR/$tar/*.$src; do
    mv $f "${f%.$src}.$tar"
  done
  echo "Creating links to warmstart files in $DA_WORKDIR_SRC/$src."

  sed -i.bak -e '/lwarmstart/s/.*/lwarmstart   = .true./g' $DA_EXPDIR/$tar"/namoptions."$tar # set warmstart to true in namoptions
  sed -i.bak -e "/startfile/s/.*/startfile    = '$startfilen\_xxx_xxx.$tar'/g" $DA_EXPDIR/$tar"/namoptions."$tar # change startfile to newest restartfiles
  rm $DA_EXPDIR/$tar"/namoptions."$tar".bak"
  # Note that this only works if lwarmstart and startfile are defined in namoptions!
  # TODO: If they are not defined under RUN section in namoptions, add them!
  echo "Warning! Check that lwarmstart = .true. and startfile = '${startfilen}_xxx_xxx.${tar}' in namoptions. If not, add them to the RUN section."

fi

# link any available driver files (without renaming)
driverfiles=$DA_WORKDIR_SRC/$src"/"?"driver_"*
for dfile in $driverfiles ; do
    if [ -f $dfile ]; then
      ln -s $dfile $DA_EXPDIR/$tar/
    fi
  done
if [ -f $dfile ]; then
    echo "Added links to driverfiles."
fi

# copy config script for execution
config_script=$DA_EXPDIR_SRC/$src/config.sh
if [ -f $config_script ] ; then
  cp $config_script $DA_EXPDIR/$tar
else
  echo 'Skipping config file (no config.sh file found).'
fi

# copy *.txt files for IBM implementation
  cp $DA_EXPDIR_SRC/$src/*.txt $DA_EXPDIR/$tar
