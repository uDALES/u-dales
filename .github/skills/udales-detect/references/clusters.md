# Cluster Profiles (Detection)

Record environment fingerprints and the minimal instructions needed to run
uDALES tools on a given system. Append new entries as they are discovered.

## Imperial CX3 (ICL RDS/EPHEMERAL)

Fingerprint:
- `hostname -f` → `login-*.cx3.hpc.ic.ac.uk` (e.g. `login-b.cx3.hpc.ic.ac.uk`)
- Lmod modules (`module` is a function; `MODULEPATH` includes `/sw-eb/modules/all`);
  `module load tools/prod` exposes the EasyBuild software tree.
- `mpiexec` → `/opt/pbs/bin/mpiexec` (PBS); real build MPI comes from the Intel modules.
- Repo git dir is the `u-dales/` subfolder of the project dir.

Route execution to `udales-exec` (ICL/CX3 profile). Key facts it needs:
- Working build stack is the self-consistent Intel **2021a** set (see udales-exec
  clusters.md) — the older `intel/2025a` mixed line no longer co-loads.
- The `~/udales/.venv` is a Python 3.9.6 venv; it needs `Python/3.9.6-GCCcore-11.2.0`
  loaded (supplies `libpython3.9.so.1.0`). A `libpython3.9.so.1.0 missing` error means
  that module isn't loaded — the venv is fine, not broken.


## ECMWF Atos HPC2020

Fingerprint:
- `hostname` → `ac6-101.bullx` style names (`bullx` domain)
- Lmod modules with `prgenv/*` families; `MODULEPATH` under `/usr/local/apps/modulefiles/lmod`
- `sbatch` at `/usr/local/bin/sbatch`, no `qsub`; `$SCRATCH` = `/ec/res4/scratch/<user>`
- `account` command lists Slurm accounts

Route execution to `udales-exec` (ECMWF HPC2020 profile). Machine name for
`tools/build_executable.sh` and `bin/ud_run` is `ecmwf`.
