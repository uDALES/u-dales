# Execute scripts and existing output directories

`run_test.sh` checks that the run wrappers never touch existing results:

- `tools/local_execute.sh` creates `$DA_WORKDIR/<expnr>` and runs a fresh case;
- when that directory already exists it exits with status 1 and an error
  message, leaves the directory untouched and does not start the solver;
- `hpc_execute.sh`, `archer_execute.sh` and `ecmwf_execute.sh` carry the same
  guard;
- `tools/copy_inputs.sh` prepares a new case without creating its output
  directory.

The solver and the MPI launcher are replaced by stubs (`DA_BUILD` and
`MPIEXEC`), so the test needs neither a build nor MPI and runs on any platform.
