#!/usr/bin/env bash

set -e

if [ -z $DA_EXPDIR ]; then
	echo "source appropriate config file to set up DA_EXPDIR."
	echo "experiment directory DA_EXPDIR must be set inside the config.sh file inside the experiment case directory."
	echo "... execution terminated"
	exit 0
fi

if [ -z $DA_WORKDIR ]; then
	 echo "source appropriate config file to set up DA_WORKDIR."
	 echo "work directory DA_WORKDIR must be set inside the config.sh file inside the experiment case directory."
	 echo "... execution terminated"
	 exit 0
fi

if (( $# < 2 )) ; then
  echo "usage: `basename $0` sim#1 sim#2"
  echo "links the driver files for target uDALES simulations that read inlet from driver files."
  echo "   sim#1: number of new target simulation case."
  echo "   sim#2: number of simulation output directory which contains the driver files."
  echo "... execution terminated"
  exit 0
fi

tar=$(printf "%03.0f" $1)    # pad argument 1 (target simulation) with zeros
src=$(printf "%03.0f" $2)    # pad argument 2 (origin) with zeros

# check if driver files already exist. If so, ask how to proceed.
if [ -f $DA_EXPDIR/$tar/"tdriver_000."$src ]; then
  echo "Few or all driver files already exist inside $DA_EXPDIR/$tar"
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
  driverfiles=$DA_WORKDIR/$src"/"?"driver_"*"."$src
  for dfile in $driverfiles ; do
      if [ -f $dfile ]; then
            ln -s $dfile $DA_EXPDIR/$tar
      fi
  done
  if [ -f $dfile ]; then
      echo "Added links to driverfiles."
  fi

fi;