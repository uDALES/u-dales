# Nudge Test

Validates `modforces::nudge`, which relaxes the upper part of the domain towards
prescribed profiles.

## What nudge does

For every level from `kb+nnudge` upwards it applies

```
Xp(:,:,k) = Xp(:,:,k) - (X0av(k) - Xprof(k)) / tnudge
```

to the velocity tendencies (when `lnudgevel`), the scalars, the temperature
(when `ltempeq`) and the moisture (when `lmoist`).

## What this test pins

Two properties that are easy to get wrong and that the CPU and GPU branches must
agree on:

1. **The starting level.** The relaxation begins at `kb+nnudge`, so everything
   below it must be left untouched.
2. **The horizontal extent.** The CPU branch uses whole-array assignment, so the
   tendency is applied across each array's entire declared extent, halo cells
   included — `ib-ih..ie+ih` for the momentum and thermodynamic tendencies, and
   the wider `ib-ihc..ie+ihc` for the scalars. The GPU branch reproduces that
   with explicit loops. A port that covered only `ib..ie` would diverge, so the
   halo cells are asserted separately rather than relying on a whole-array
   comparison that an interior-only port would still pass.

It also checks that `lnudge = .false.` is a genuine no-op.

## Scope

This runmode covers the **host** branch. Running it against a GPU build fails
with an explanatory message rather than passing vacuously, because `nudge`
writes the device-resident tendencies there and the host arrays this test
inspects would never change.

The device kernels are covered in two other places:

- `tests_cuda.f90::test_nudge_profiles` — device self-test, runs under
  `UDALES_RUN_CUDA_SELFTEST=1` on a Debug GPU build.
- The `nudging` case in `tests/integration/gpu/case_matrix.json` — CPU/GPU
  parity on a real run.

## Requirements

- A CPU debug build: `tools/build_executable.sh common debug`
- MPI (`mpiexec`). Runs on a single rank.

## Running

```bash
tests/integration/nudge/run_test.sh
```

Overridable via environment: `UDALES_BUILD`, `CASE_SOURCE`, `NAMELIST_SOURCE`,
`NAMELIST`, `NPROCS`, `MPIEXEC`.

## Implementation

- Runmode `1007` (`TEST_NUDGE` in `src/modglobal.f90`)
- `tests_nudge` in `src/tests.f90`, dispatched from `execute_runmode_actions`
  in `src/program.f90`

## See also

`tests/integration/ibm_cell_lookup` follows the same runmode-test pattern.
