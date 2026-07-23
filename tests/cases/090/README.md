# Case 090

Baseline case for the `BCcleanup` regression harness
(`tests/regression/bc_cleanup/`). Copied from `tests/cases/101/`
(the IBM sparse-input fixture) and adapted to a short, dry, IBM-buildings
run.

Historically (before Task 3 of the `BCcleanup` refactor), this case
carried matched anchors between the initial profile and the boundary
condition: `thls = 290.` in `namoptions.090` equaled the bottom-row
`thl` in `prof.inp.090` (z = 0.5 m), and `qts = 0.0` equaled the
bottom-row `qt`. That made the old `thvs` (computed from `thls`/`qts`)
equal the new `thv_b(kb)` (computed from the profile), so Phase 0
commits that were expected to be no-ops could be checked at
`--atol 0.0`. As of Task 3, `thls`/`qts` (and the other flat-surface
BC keys) have been pruned from `&BC` entirely, so `namoptions.090` no
longer sets them — the base state derives from `prof.inp` alone. The
bottom-row `thl`/`qt` in `prof.inp.090` are unchanged, so the anchor
property still holds automatically.

Differences from case 101:

- `iexpnr = 090` (was `101`).
- `runmode` removed (defaults to `RUN_COLDSTART`); case 101 uses the
  sparse-IJK test runmode (`1004`), which exits before any time
  integration and produces no dumps.
- `runtime = 20.` (was `101.`).
- `&OUTPUT`: `ltdump = .true.` added and `tstatsdump = 5.` (was `10.`)
  so at least 3 `tdump.*.090.nc` time slices are written within the
  20 s runtime.

Geometry (`stl_file = geom.101.STL`), `libm`, and the grid are
unchanged from case 101; the STL file keeps its original name because
`tests/cases/090/` was created by copying `101/` and renaming only the
`*.101`-suffixed input files to `*.090`.

Current consumers:

- `tests/regression/bc_cleanup/run_test.py`
