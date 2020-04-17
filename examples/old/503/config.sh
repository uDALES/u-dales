#!/bin/bash

# set environmental variables for uDALES
DA_UTILSDIR=$DA_TOPDIR/u-dales/tools/utils # Directory of utils scripts
DA_BUILD=$DA_TOPDIR/u-dales/build/test/u-dales # Executable
NCPU=1 # Number of CPUs to use for a simulation

# local only
DA_WORKDIR=$DA_TOPDIR/outputs # Output top-level directory

# HPC only
NNODE=1 # Number of nodes to use for a simulation
WALLTIME="00:30:00" # Maximum runtime for simulation in hours:minutes:seconds
MEM="2gb" # Memory request per node
