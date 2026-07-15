# bc_cleanup regression

Phase 0 regression gate for `BCcleanup_backlog.md`. Every commit in the
`BCcleanup` refactor is checked with this harness before moving to the
next step: it builds two refs (a reference and the current workspace),
runs case `090` under two MPI decompositions on each, and compares the
dumped fields.

This is a stripped-down adaptation of
`tests/regression/mpi_averaging_regression/run_test.py`: the
worktree/build/run/stitch machinery (git worktrees, offline FFTW
patching, `mpiexec` invocation, per-tile NetCDF stitching) is unchanged.
What changed is the case matrix (one case, not two), the decomposition
set (`serial` and `xy_split`, not four), and the field comparison
(every `ndim >= 3` variable in the dump, not a fixed 3-tuple), plus a
`--atol` knob so later, behaviour-changing commits can be checked
against a small non-zero tolerance instead of bitwise equality. The
tree-case (`526`) and `treedump`/vegetation-mask handling from the
parent script were deleted; case 090 has no vegetation.

## Case 090: matched-anchor construction

Case `090` (`tests/cases/090/`) is copied from `tests/cases/101/` (dry,
`ltempeq = .true.`, IBM buildings) and adjusted so that:

- `thls = 290.` in `namoptions.090` equals the bottom-row `thl` in
  `prof.inp.090` (z = 0.5 m).
- `qts = 0.0` in `namoptions.090` equals the bottom-row `qt` in
  `prof.inp.090`.

This "matched anchor" makes the old surface virtual temperature
`thvs` (computed from `thls`/`qts`) numerically equal to the new
`thv_b(kb)` (computed from the profile at the lowest level), so any
BCcleanup commit that is meant to be a pure refactor of how the
surface virtual temperature is derived can be checked at `--atol 0.0`:
if the two derivations agree at the matched anchor, a bitwise field
match confirms the refactor introduced no behavioural change on this
case.

Case 090 also enables `ltdump = .true.` with `tstatsdump = 5.` (and a
20 s `runtime`), so at least 3 `tdump.*.090.nc` time slices exist to
compare per run. See `tests/cases/090/README.md` for the full list of
changes from case 101.

## Tolerance regime per BCcleanup commit

The harness is generic; each task in `BCcleanup_backlog.md` decides its
own `--atol` depending on whether the commit is expected to be a
bitwise no-op or to change floating-point results at roundoff:

| Task   | `--atol`  | Rationale                                             |
| ------ | --------- | ------------------------------------------------------ |
| Task 2 | `0.0`     | Pure refactor, no behavioural change expected.          |
| Task 3 | `0.0`     | Pure refactor, no behavioural change expected.          |
| Task 4 | `5e-3`    | Reordering/reformulation expected to shift results at this scale. |
| Task 5 | `0.0`     | Pure refactor, no behavioural change expected.          |

Do not loosen `--atol` to make a failing comparison pass; a failure at
`0.0` for a task that expects `0.0` means the commit changed behaviour
and needs to be fixed, not the tolerance.

## Usage

```bash
python tests/regression/bc_cleanup/run_test.py BCcleanup BCcleanup Debug --atol 0.0
```

`branch_a` is the reference ref (built in a throwaway git worktree).
`branch_b` is only used when `--ci` is passed (it is otherwise a label;
without `--ci` the "current" side builds the current workspace as
checked out). Useful flags, inherited from the parent script:

- `--workdir DIR` — root directory for temporary worktrees and run
  outputs (defaults to `$UDALES_SCRATCH_ROOT` or `/tmp`).
- `--ci` — also build `branch_b` in its own worktree instead of using
  the current workspace.
- `--reuse-current-build` / `--reference-build-dir` — skip rebuilding
  when a suitable `u-dales` executable already exists.
- `--compare-logs` — additionally compare solver diagnostics (Courant,
  diffusion, divergence) between the two runs.
- `--configs serial,xy_split` — restrict which decompositions run.
- `--cleanup` — remove temporary worktree/run directories on exit
  (off by default so a failing run leaves its artifacts for
  inspection).

Exit code is `0` when every compared field is within `--atol` on both
decompositions, non-zero otherwise (including when zero matching
NetCDF dump files are found — that is always a failure, never a
silent pass).
