#!/bin/bash

## get environment variables
topdir=${DA_TOPDIR}
srcdir=${DA_SRCDIR}
expdir=${DA_EXPDIR}
workdir=${DA_WORKDIR}
utldir=${DA_UTILSDIR}
location=${PLATFORM}

echo -e "we are at location=${location}"
echo -e "variables set to:\n topdir=${topdir}\n srcdir=${srcdir}"
echo -e "expdir=${expdir}\n workdir=${workdir}\n utldir=${utldir}"

## experiment number is automatically set via path
cdir=$(pwd)
iexpnr="${cdir: -3}"
echo "experiment number: ${iexpnr}"

## get platform dependent variables from execute-part1 script
. ${utldir}/exe_setvar.sh

## specify parameters now if you want to change from default
# exe=dalesurban_cx1
# ncpu=12
# nnode=4
# walltime="48:00:00"
# mem=20
# queue=""
# queue="#PBS -q pqcdt"
# modules="intel-suite mpi/intel-5.1"

## execute simulation
. ${utldir}/exe_submit.sh
