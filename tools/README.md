# uDALES workflow scripts

The shell wrappers in this directory automate the simulation lifecycle and
are the preferred interface — use them before hand-rolling
cmake/mpiexec/qsub commands. Fuller workflow context: `docs/udales-workflow.md`
and `docs/udales-simulation-setup.md`.

## The `config.sh` contract

Every case directory carries a `config.sh` that the wrappers source:

```bash
export DA_EXPDIR=...    # parent directory of the case directories
export DA_WORKDIR=...   # run directory (e.g. $EPHEMERAL on ICL)
export DA_TOOLSDIR=...  # this directory
export DA_BUILD=...     # path to the u-dales executable
# cluster runs additionally:
export NCPU=... NNODE=... WALLTIME="hh:mm:ss" MEM="..."
```

The case number (`iexpnr`) is derived from the LAST THREE characters of the
case directory name — directories may be purely numeric (`900`) or
descriptive with a numeric suffix (`benchmark-standard-900`).

## Lifecycle

| Step | Script | Notes |
|------|--------|-------|
| build solver | `build_executable.sh <system> <debug\|release>` | system = `icl`, `archer`, `cca`, `common`; sets modules + netCDF paths per platform |
| build preprocessing | `build_preprocessing.sh <system>` | View3D and facet-section tooling |
| new case from old | `copy_inputs.sh <src_case_path> <new_case_number> [c\|w]` | cold- or warm-start setup |
| preprocess geometry | `write_inputs.sh <case_path> [c\|l]` | runs the MATLAB preprocessing; `c` submits a compute-node batch job (use for anything big), `l`/default runs where you are; fills the `&WALLS` counts in namoptions |
| run locally | `local_execute.sh <case_path>` | workstation/login-node mpiexec run |
| run on ICL HPC (PBS) | `hpc_execute.sh <case_path>` | writes a PBS job from config.sh and submits it |
| run on ARCHER2 (Slurm) | `archer_execute.sh <case_path>` | Slurm equivalent |
| merge outputs | `gather_outputs.sh <exp_dir>` (+ `hpc_gather.sh`, `archer_gather.sh`) | merges per-rank NetCDF via NCO (`nco_concatenate_field*.sh`) |
| join two runs | `append_outputs.sh <exp1_dir> <exp2_dir>` | concatenate outputs across a restart |
| driver/inflow | `generate_synthetic_inflow.sh`, `link_driver_files.sh` | precursor/driver-simulation workflows |

## Porting to a new platform

The execute/gather scripts target the two production platforms (ICL CX3 =
PBS, ARCHER2 = Slurm) but are deliberately thin:

1. Copy the closest `*_execute.sh` / `*_gather.sh` pair and adapt the
   scheduler directives and module-load lines.
2. Add a system entry (modules, `FC`, netCDF prefixes) in
   `build_executable.sh`.
3. Record the new platform's stack and gotchas in
   `.github/skills/udales-exec/references/clusters.md`, and its performance
   profile in `.github/skills/udales-perf/references/machines.md`.
