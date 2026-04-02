# Legacy Tree Regression

This experimental regression test protects the legacy vegetation behaviour of
this branch against the `v2.2.0` release baseline.

It builds:

- `reference` from the `v2.2.0` release tag, using the legacy/block tree path
- the current branch/revision, also using the legacy/block tree path

By default, the current side uses the workspace build so the regression can be
used during local development before changes are committed. In `--ci` mode, the
current side is built from a clean detached worktree so the result is fully
reproducible from Git state.

Then it stages case `526` twice, forces `itree_mode = 99` on the current branch copy, runs both executables, and compares the
resulting `treedump.*.*.526.nc` outputs.

On Imperial HPC, it uses the existing `tools/build_executable.sh icl ...`
wrapper inside each temporary worktree so the regression build matches the
canonical local build settings.

The temporary worktrees reuse the already-populated local copies of
`2decomp-fft` and `tools/View3D` from the main checkout rather than fetching
submodules over HTTPS.

Current pass criteria:

- acceptance is based only on the interior relative errors of `tr_u`, `tr_v`,
  and `tr_w`
- the interior relative threshold is `1e-3`
- side/support-cell momentum differences are reported for diagnostics only
- `tr_qt` and `tr_thl` are also reported for diagnostics, but are not currently
  used as pass/fail criteria

This matches the intent of the original sparse-versus-block tree comparison:
interior canopy momentum should stay close, while side/support cells are more
sensitive to representation details in the full simulation path.

In particular, the side/support-cell values for `u`, `v`, and `w` are not
expected to match the original version exactly, because this branch applies the
staggering treatment more consistently than the old block-tree implementation.
Those side differences are therefore informative, but they are not used as the
acceptance criterion for this regression.

This test uses temporary `git worktree` checkouts while building, so it does
not need to switch the user’s active checkout or require a clean worktree.

By default, temporary worktrees, staged run directories, and scratch NetCDF
comparison copies are created under `$EPHEMERAL` when that environment
variable is available. You can override that location with
`UDALES_SCRATCH_ROOT`.

Typical entry points:

- direct regression run:
  `python tests/regression/new_vegetation_module_against_v2.2/run_test.py v2.2.0 HEAD Debug`
- CI-style reproducible run:
  `python tests/regression/new_vegetation_module_against_v2.2/run_test.py v2.2.0 HEAD Debug --ci`
- analysis only on existing outputs:
  `python tests/regression/new_vegetation_module_against_v2.2/run_test.py --analyze-only --report-only --reference-run-dir <reference_run_dir> --current-run-dir <current_run_dir> --case-dir tests/cases/526`
- curated experimental suite:
  `python tests/run_tests.py experimental --branch-b HEAD --build-type Debug`
