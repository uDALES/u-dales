#!/bin/bash
# MPI viability probe for ARCHER2 (Slurm), adapted from probe_mpi.sh (PBS).
#
# Submit this once per node class/QoS before committing a benchmark campaign.
# ARCHER2 allocates whole nodes per job, so "node class" here really means
# partition/QoS combination (standard vs short vs highmem) rather than a
# shared-node core count as on PBS clusters.
#
# Usage: sbatch probe_mpi_archer2.sh

#SBATCH --job-name=mpi-probe
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --account=n02-assure
#SBATCH --partition=standard
#SBATCH --qos=short

module load epcc-job-env >/dev/null 2>&1 || true

echo "host: $(hostname)"
grep -m1 'model name' /proc/cpuinfo
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
which srun

echo "--- T1 srun hostname (default) ---"
timeout 60 srun --distribution=block:block --hint=nomultithread -n 4 hostname && echo T1_OK || echo T1_FAIL

echo "--- T2 srun --mpi=pmi2 hostname ---"
timeout 60 srun --mpi=pmi2 -n 4 hostname && echo T2_OK || echo T2_FAIL

echo "probe done"
