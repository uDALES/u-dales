#!/bin/bash

# set environmental variables for uDALES
DA_UTILSDIR=$DA_TOPDIR/u-dales/tools/utils # Directory of utils scripts
DA_BUILD=$DA_TOPDIR/u-dales/build/release/u-dales # Executable
NCPU=2 # Number of CPUs to use for a simulation

# local only
OUTPUT=$DA_TOPDIR/outputs # Output top-level directory

# HPC only
NNODE=5 # Number of nodes to use for a simulation
WALLTIME="24:00:00" # Maximum runtime for simulation in hours:minutes:seconds
MEM="20gb" # Memory request per node
