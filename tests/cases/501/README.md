# Cases 501 / 502 — precursor → driver inflow fixture

`501` (`idriver = 1`) writes inlet planes at `iplane`; `502` (`idriver = 2`)
consumes them. They are run as a pair by
`tests/integration/driver_inflow/run_test.sh`, which is registered as a
`supported` suite so CI runs it on every push.

## Why they exist

The driver lifecycle in `src/inflow.f90` (`initinflow`, `drivergen`,
`writedriverfile`, `readdriverfile`, and both branches of `exitinflow`) had no
runnable test. The existing driver cases are production-sized: `examples/949` is
256×128×128, `examples/950` consumes its planes, and `tests/cases/525` asks for
`runtime = 691200.` on 4096 ranks. None of them can run in CI, so the module was
only ever checked by reading it.

These two are 8×8×8 on 2 ranks and finish in seconds.

## Construction

Geometry (flat ground, 128 facets) is the 8³ fixture from
`tests/regression/david_tests/cases/103`; the grid is the smallest that already
had ground facets, so no preprocessing step is needed to build the pair.

Configuration follows what production actually uses rather than what was
convenient:

- **Vreman closure** (`lvreman = .true.`). Every case under `experiments/`
  (064, 065, 103, 105, 110, 525, 531) uses Vreman; none uses `loneeqn`.
- **`lfftwmeasure = .false.`**, so the FFT plans are reproducible. See
  `tests/regression/bc_cleanup/README.md` for why the default traps.
- **`ltempeq`/`lmoist` on, `BCxT = 3`, `BCxq = 3`.** This is a deliberate
  departure from `experiments/525`, which is neutral. Those switches are what
  set `lhdriver`/`lqdriver`, so the thl and qt driver files are written, read and
  deallocated as well as u/v/w. A neutral fixture would leave those branches —
  including half of `exitinflow` — untested. The test asserts the fallback
  warning ("not given by precursor") is absent, so a silent downgrade to a
  non-precursor inflow fails rather than passes.

## Timing contract

`501` records `driverstore = 5` planes from `tdriverstart = 1.` every
`dtdriver = 0.5`, i.e. at t = 1.0 … 3.0, and runs to `runtime = 4`. The solver
overrides `trestart` to `tdriverstart + (driverstore-1)*dtdriver` = 3.0 and says
so at startup.

Driver times are stored relative to `tdriverstart`, so `502` sees t = 0.0 … 2.0
and runs to `runtime = 1.5`, inside that span. Raising `502`'s `runtime` above
2.0 makes the solver warn that it will stop early, and beyond the last plane the
run ends. If you change any of `tdriverstart`, `dtdriver`, `driverstore` or
either `runtime`, keep `502`'s `runtime` below `(driverstore-1)*dtdriver`.

`jtot`, `ktot` and `nprocy` must match between the pair: the planes are per
y-rank slabs and their record length is fixed by the precursor's decomposition.
