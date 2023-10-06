## To convert a text file filename.sh to an exicutable bash file filename.sh, use the following command
##	~$ chmod 777 filename.sh


## Usage: ./u-dales/tools/generate_synthetic_inflow.sh <PATH_TO_CASE>
## FROM THE PARENT DIRECTORY

#!/usr/bin/env bash

set -e

if (( $# < 1 )) 
then
    echo "The path to case/experiment folder must be set."
    exit
fi

## go to experiment directory
pushd $1
inputdir=$(pwd)

## set experiment number via path
exp="${inputdir: -3}"

echo "Creating driver files for uDALES target simulation case $exp..."

## read in additional variables
if [ -f config.sh ]; then
    source config.sh
fi;

## check if required variables are set
if [ -z $DA_TOOLSDIR ]; then
    echo "Script directory DA_TOOLSDIR must be set"
    exit
fi;

## compiling synthetic inflow generator source code
gfortran -O2 -fopenmp $DA_TOOLSDIR/syntheticInflow/modSyntheticInflow.f90 -o synInflow_executable

## copying syntheticInfow inputs to the experiment directory   (may be modified later)
cp syntheticInflow_inputs/*.txt .

## setting up OpenMP environmental variables
export OMP_NUM_THREADS=8    ## please increase this for large job based on availability of cores
export OMP_PLACES=cores

## run
./synInflow_executable $exp 2>&1 | tee -a log_synthetic_inflow_compute.$exp

## cleaning experiment dircetory
rm -f synInflow_executable
rm -f ut.txt upup.txt upvp.txt upwp.txt vpvp.txt vpwp.txt wpwp.txt
rm -f thlt.txt thlpthlpt.txt wpthlpt.txt
rm -f qtt.txt
rm -f length_time_scales_u.txt length_time_scales_v.txt length_time_scales_w.txt
rm -f length_time_scales_temp.txt length_time_scales_qt.txt

popd
