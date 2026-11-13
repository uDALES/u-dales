#!/bin/bash
# MPI viability probe (udales-perf skill, Phase 2).
#
# Submit this as a batch job ONCE PER NODE CLASS before committing any
# benchmark campaign to that class. Edit the scheduler directives and the
# module line to match the target machine/class (the resource request is what
# selects the class - probe with exactly the request the campaign will use).
#
# Interpreting results: if ssh/pbs_tmrsh to self succeed but every mpiexec
# variant fails with "error setting up the bootstrap proxies", the MPI stack
# is incompatible with this node class's OS image (observed on ICL CX3
# racks 12-15 with Intel MPI 2021.2) - use a different class or toolchain;
# no launcher flag will fix it.

#PBS -l walltime=00:15:00
#PBS -l select=1:ncpus=64:mpiprocs=64:mem=8gb
### Adapt the select line to the node class being probed, e.g. ncpus=128 for
### whole-node classes, and add site-specific selectors (cpu_type=...).

# module load <the exact stack the campaign will use>

echo "host: $(hostname)"
grep -m1 'model name' /proc/cpuinfo
echo "nodefile entries: $(wc -l < $PBS_NODEFILE 2>/dev/null || echo n/a)"
which mpiexec

echo "--- ssh to self ---"
timeout 30 ssh -o BatchMode=yes -o StrictHostKeyChecking=no $(hostname) true && echo SSH_OK || echo SSH_FAIL
echo "--- pbs_tmrsh to self ---"
timeout 30 pbs_tmrsh $(hostname) true 2>/dev/null && echo TMRSH_OK || echo TMRSH_FAIL

echo "--- T1 plain mpiexec ---"
timeout 90 mpiexec -n 4 hostname && echo T1_OK || echo T1_FAIL
echo "--- T2 -bootstrap ssh ---"
timeout 90 mpiexec -bootstrap ssh -n 4 hostname && echo T2_OK || echo T2_FAIL
echo "--- T3 I_MPI_HYDRA_BOOTSTRAP=fork ---"
I_MPI_HYDRA_BOOTSTRAP=fork timeout 90 mpiexec -n 4 hostname && echo T3_OK || echo T3_FAIL

echo "probe done"
