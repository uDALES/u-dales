#!/bin/bash

utildir='/projects/dales-u/utils'

echo $utildir

if (( $# == 1 )) ; then
    ## set experiment number from argument
    iexpnr=$1
    echo "experiment number: ${iexpnr}"

elif (( $# == 0 )) ; then
    ## set experiment number via path
    cdir=$(pwd)
    iexpnr="${cdir: -3}"
    echo "experiment number: ${iexpnr}"

else
	echo "error: call execute scipt as `basename $0` expnr or directly from exp directory."
	exit 0
    
fi

## get default variables
source ${utildir}/exe_setvar.sh

## specify parameters if you want to change from default

srcdir='/projects/dales-u/src/dales-urban'
workdir='/projects/dales-u/work'
exe=dalesurban_local
ncpu=16
# nnode=4
# walltime="48:00:00"
# mem=20
# queue=""
# queue="#PBS -q pqcdt"
# modules="intel-suite mpi/intel-5.1"

## execute simulation
source ${utildir}/exe_submit.sh
