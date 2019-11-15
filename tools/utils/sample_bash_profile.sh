# set environmental variables for uDALES
export DA_TOPDIR=$HOME/uDALES  # This is your top-level project directory.
export DA_EXPDIR=$DA_TOPDIR/experiments #  The top-level directory of the simulation setups.
export DA_WORKDIR=$DA_TOPDIR/outputs # Output directory
# export DA_WORKDIR=$EPHEMERAL # Output directory on HPC
export DA_UTILSDIR=$DA_TOPDIR/u-dales/tools/utils # Directory of utils scripts
export DA_BUILD=$DA_TOPDIR/u-dales/build/release/u-dales # Executable
export LOCAL_EXECUTE=1 # Do not set when executing on ICL HPC (used by `mergehelper.sh`).
export NCPU=2 # Number of CPUs to use for a simulation

export PATH=$PATH:$DA_UTILSDIR  # to call functions in utils from anywhere

# shortcuts and settings
alias ll='ls -ahl --color=auto'
alias ls='ls --color=auto'
alias grep='grep -In --color=auto'
alias vi='vim'
alias rm='rm -I'