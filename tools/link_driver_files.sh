#!/usr/bin/env bash

set -e

if (( $# < 2 )) ; then
  echo "usage: `basename $0` sim#1 sim#2"
  echo "links the driver files for target uDALES simulations that read inlet from driver files."
  echo "   sim#1: path to output directory of driver simulation."
  echo "   sim#2: path to experiment directory of target simulation."
  echo "... execution terminated"
  exit 0
fi

## go to driver output directory
pushd $1 > /dev/null
srcdir=$(pwd)
popd > /dev/null

## go to target experiment directory
pushd $2 > /dev/null
tardir=$(pwd)
popd > /dev/null

src="${1: -3}"

# check if driver files already exist. If so, ask how to proceed.
if [ -f $tardir/"tdriver_000."$src ]; then
  echo "Few or all driver files already exist inside $tardir"
  echo "continue with overwrite/abort? (o/a)"
  read answer
  if [ "$answer" == "o" ]; then
    case=1
  else
    echo "abort"
    exit 1
  fi
else
  case=1
fi

if [ case==1 ]; then

  # link any available driver files (without renaming)
  driverfiles=$srcdir"/"?"driver_"*"."$src
  
  # check if any driver files exist in source directory by counting ls output
  file_count=$(ls $driverfiles 2>/dev/null | wc -l)
  
  if [ "$file_count" -eq 0 ]; then
      echo "Error: No driver files found in source directory $srcdir"
      echo "Looking for pattern: *driver_*.$src"
      exit 1
  fi
  
  # create symbolic links for all found driver files
  for dfile in $driverfiles ; do
      if [ -f $dfile ]; then
            ln -sf $dfile $tardir
      fi
  done
  
  echo "Added links to $file_count driverfiles."

fi;
