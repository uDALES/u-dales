#!/bin/bash

set -e

if (( $# < 1 ))
then
 echo "The output case directory must be set."
 echo "usage: FROM THE TOP LEVEL DIRECTORY run: bash ./u-dales/tools/archer_gather.sh <PATH_TO_output_CASE>"
 exit 1
fi

## go to output case directory
pushd $1
outdir=$(pwd)

## set experiment number via path
exp="${outdir: -3}"

## read in additional variables
if [ -f config.sh ]; then
    source config.sh
else
    echo "config.sh must be there inside $outdir"
    exit 1
fi

## check if required variables are set
if [ -z $DA_TOOLSDIR ]; then
    echo "Script directory DA_TOOLSDIR must be set inside $outdir/config.sh"
    exit 1
fi;
if [ -z $NNODE ]; then
    echo "Script directory NNODE must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $NCPU ]; then
    echo "Script directory NCPU must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $WALLTIME ]; then
    echo "Script directory WALLTIME must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $QOS ]; then
    echo "Script directory QOS must be set inside $inputdir/config.sh"
    exit 1
fi;


###### Job script
cat <<EOF > post-job.$exp.slurm
#!/bin/bash
#SBATCH --job-name=${exp}_gather
#SBATCH --time=${WALLTIME}
#SBATCH --nodes=${NNODE}
#SBATCH --tasks-per-node=${NCPU}
#SBATCH --cpus-per-task=1
#SBATCH --account=n02-ASSURE
#SBATCH --partition=standard
#SBATCH --qos=${QOS}

$DA_TOOLSDIR/gather_outputs.sh $outdir
EOF

## submit job file to queue
sbatch post-job.$exp.slurm

echo "post-job.$exp.slurm submitted."

popd
