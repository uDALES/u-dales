# Run with: tools/hpc_execute.sh examples/304   (from the repository root)
# Paths default to this repository layout on CX3; override by exporting the variable first.
export DA_EXPDIR=${DA_EXPDIR:-$HOME/udales/u-dales/examples}
export DA_TOOLSDIR=${DA_TOOLSDIR:-$HOME/udales/u-dales/tools}
export DA_BUILD=${DA_BUILD:-$HOME/udales/u-dales/build/release/u-dales}
export DA_WORKDIR=${DA_WORKDIR:-/rds/general/ephemeral/user/$USER/ephemeral/udales_outputs}
export NNODE=1
export NCPU=16         # = nprocx * nprocy in namoptions.304
export WALLTIME="02:00:00"
export MEM="16gb"
