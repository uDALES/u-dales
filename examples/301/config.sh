# Run with: tools/hpc_execute.sh examples/301   (from the repository root)
# Paths default to this repository layout on CX3; override by exporting the variable first.
# WALLTIME and MEM are provisional: sized from 60 s smoke tests on 4 ranks (see README.md).
export DA_EXPDIR=${DA_EXPDIR:-$HOME/udales/u-dales-thermo/examples}
export DA_TOOLSDIR=${DA_TOOLSDIR:-$HOME/udales/u-dales-thermo/tools}
export DA_BUILD=${DA_BUILD:-$HOME/udales/u-dales-thermo/build/release/u-dales}
export DA_WORKDIR=${DA_WORKDIR:-/rds/general/ephemeral/user/$USER/ephemeral/udales_outputs}
export NNODE=1
export NCPU=8          # = nprocx * nprocy in namoptions.301
export WALLTIME="00:30:00"
export MEM="16gb"
