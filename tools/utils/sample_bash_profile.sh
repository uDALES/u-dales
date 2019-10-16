export PS1='\[\e]2;\u@\h:$(pwd)\007\e]1;\h\007\]\u@\h:\w>'  # shows hostname

# set environmental variables for dales
export DA_TOPDIR=$HOME/uDALES
export DA_EXPDIR=$DA_TOPDIR/experiments
export DA_WORKDIR=$EPHEMERAL

export PATH=$PATH:$DA_UTILSDIR  # to call functions in utils from anywhere

# load modules
module load git
module load anaconda3/personal  # this is python
module load ncview

# export PATH=$PATH:~\bin    # not sure this is needed?

# shortcuts and settings
alias ll='ls -ahl --color=auto'
alias ls='ls --color=auto'
alias grep='grep -In --color=auto'
alias vi='vim'
alias rm='rm -I'
alias wd='cd $WORK'
alias xd='cd $DA_SRCDIR'