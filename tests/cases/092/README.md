# Case 092

One-equation-SGS variant of case 090, for the `BCcleanup` regression harness
(`tests/regression/bc_cleanup/`). Copied from `tests/cases/090/` (itself
copied from `tests/cases/101/`, the IBM sparse-input fixture), renaming the
`*.090`-suffixed input files to `*.092`, with `iexpnr = 092` and one change
to `&NAMSUBGRID`:

- `lvreman = .false.` (was `.true.`)
- `loneeqn = .true.` (added)

This exists because case 090 leaves `lvreman = .true.` (the default), so it
only exercises the Vreman branch of `closure` in `src/modsubgrid.f90`; the
`zlt`/`sbbuo` lines re-pointed to the evolving `thvf(k)` (`src/modsubgrid.f90:387,486`)
live in the `loneeqn` branch and are never executed by 090 — see the
correction note in `BCcleanup_backlog.md` (§6, Phase 0) and
`docs/superpowers/plans/2026-07-16-bc-cleanup-phase0.md` (Task 4). Case 092
flips the SGS scheme to one-equation TKE closure so that branch actually
runs. `prof.inp.092`'s sixth column (`e12`) seeds the TKE field; any
initial value below `e12min` (5e-5, `src/modglobal.f90`) is clamped up on
the first step, so the all-zero column in the copied profile is fine.

All other namelist blocks, geometry (`stl_file = geom.101.STL`), `libm`,
and the grid are unchanged from case 090.

Current consumers:

- `tests/regression/bc_cleanup/run_test.py` (`tdump.*.092.nc`, default
  tolerance `5e-3`)
