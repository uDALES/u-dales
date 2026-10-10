# SGS statistics integration test

Validates the subgrid-scale (SGS) fluxes written by `modstatsdump` against both
an analytic manufactured solution and the diffusion tendencies the solver
actually applies.

- Suite class: `experimental` (not part of the required merge gate)
- Kind: `integration`
- Component: `solver`
- Runmode: `1006` (`TEST_SGS_STATS`, dispatched from `program.f90` to
  `tests_sgs_statistics` in `src/tests.f90`)
- Fixture: `tests/cases/300` — a small (16x16x32) case with a **stretched**
  vertical grid and `libm = .false.`

## What it checks

The test calls the same public kernel the solver uses,
`modstatsdump::compute_sgs_fluxes`, so a change to the SGS stencils cannot drift
away from what is tested.

1. **Manufactured solution on a stretched grid.** With constant `ekm`/`ekh` and
   `u = a z^2`, `w = b x`, `thl = c z^2`, the SGS fluxes must equal
   `-nu (2 a z + b)` and `-kappa 2 c z`. Two assertions are made: an exact
   discrete identity (machine precision) and second-order consistency with the
   analytic profile, bounded by `max_k |dzf(k) - dzf(k-1)| / 2`. A stretched
   grid is essential here — on a uniform grid a `dzf`/`dzh` confusion is
   invisible.
2. **Consistency with the solver (anti-drift).** With pseudo-random `ekm`/`ekh`
   and velocity/temperature fields that vary in `z` only (and `w` in `x` only),
   the horizontal terms of `diffu`/`diffc` vanish identically, so the tendency
   they add to `up`/`thlp` must equal minus the vertical divergence of the
   statistics' SGS fluxes, to machine precision.
3. **Sign regression.** For a monotone positive shear the SGS momentum and heat
   fluxes must be negative (down-gradient). This is a one-line guard on the sign
   convention fixed in issue #306: the stored quantity is the physical flux
   `-K dphi/dn`, not `+K dphi/dn`.

All prescribed fields are pure functions of the **global** coordinates, so the
halos are filled analytically. The test therefore needs no halo exchange and is
decomposition-invariant; it is run on `1x1`, `2x1`, `1x2` and `2x2`
decompositions.

## Prerequisites

- A built solver. The path is taken from `UDALES_BUILD`, defaulting to
  `build/debug/u-dales`.
- An MPI launcher able to start 4 ranks (`--oversubscribe` is added
  automatically for Open MPI).

## Running

```bash
cd tests/integration/sgs_statistics
./run_test.sh
```

or via the dispatcher:

```bash
python tests/run_tests.py experimental
```

The script greps the run log for `ALL TESTS PASSED: tests_sgs_statistics` and
also checks the process exit code.
