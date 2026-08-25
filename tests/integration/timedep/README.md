# Time-Dependent Forcing Test

Drives `modtimedep::timedepsurf` and `modtimedep::timedepnudge`, the two
interpolations that move the prescribed surface fluxes and the nudging target
between the entries of a forcing table.

## Why this test exists

Every number `timedep` produces reaches the solution indirectly. The five
surface fluxes go into a facet wall function; the four nudging profiles go into
a tendency divided by `tnudge`; both are then integrated. A regression case can
tell that something moved. Only this can tell that it moved to the right value
at the right time.

There was no coverage of any kind before: no unit test, no integration test,
and no case in the repository setting any of the four `ltimedep*` switches.

## Tables, not files

A runmode test has no input file to read, so the tables are installed directly
in the module. That is also what makes the checks sharp. A table read from a
file is smooth in `k` and in `t`, and a smooth table hides a transposed index.
Here every entry is unique: linear in `k`, a decade apart between the four
profiles, and with a `t` dependence of a different shape from the `k` one.

## What the checks pin

| check | what it pins |
| --- | --- |
| the bracket | the search walks down from the last entry and takes the first whose time has passed, so an interval boundary belongs to the interval **above** it. Off by one and the whole forcing runs an interval early or late, which on a diurnal file is hours. |
| the endpoints | at a table time `fac` is zero and the result is the table entry itself, **bit for bit** - `a + 0*(b - a)` is `a` in IEEE arithmetic. These are checked with `/=`, so the obvious tidy-up to `(1-fac)*a + fac*b`, which is not the same expression, fails them. |
| the interior | two intervals of different widths, at two different fractions, so a `fac` divided by the wrong width is caught. |
| the freeze | it starts **at** the last table time, not after it. The search returns the last index there and the guard writes nothing, so the last column is never installed and what stands from then on is the interpolated value from the step before. In a run that is the last column to within one step, which is why it has never mattered; it is still not the same thing, and a port that clamped to the last column instead would pass every other check here. Poisoned first, so a routine that ran anyway and landed on the last column is not mistaken for one that correctly did nothing. |
| one column each | each of the five fluxes and each of the four profiles carries its own column of the table, at its own level. Five scalars of the same magnitude interpolated in one routine is the shape a copy-paste error survives in. |
| the switches | each routine returns before touching anything when its own switch is off. Setting one of the four is the normal way to use the module. |

## Coverage

Host code, and therefore worth running against either build - both routines are
compiled into the GPU build unchanged, and `timedep_device` calls `timedepsurf`
directly. What that leaves uncovered here is the OpenACC branch of the profile
interpolation, which is:

- `tests_cuda.f90::test_timedep_nudge_device`, under `UDALES_RUN_CUDA_SELFTEST`
  on a Debug GPU build, and
- the `timedep-nudging` case in `tests/integration/gpu/case_matrix.json`.

`ltimedepsurf` has no end-to-end coverage on either device. Its five fluxes are
read by the facet wall function under `iwalltemp = 1`, and no case in the
repository combines that with `libm`. This test is the only thing that exercises
the interpolation at all; the gap in the wall function predates it.

`timedeplw` and `timedepsw` are not driven here. Both need an initialised facet
set and both feed the energy balance, which is host code in either build and is
covered by `tests_eb`.

## Running

```bash
tests/integration/timedep/run_test.sh
```

Environment overrides: `UDALES_BUILD`, `CASE_SOURCE`, `NAMELIST_SOURCE`,
`MPI_NAMELIST_SOURCE`, `NPROCS`, `MPIEXEC`.

## Passes

| pass | what it adds |
| --- | --- |
| single rank | the interpolation itself |
| two ranks in y | that every rank brackets the same way. The profiles are not decomposed and each rank interpolates its own copy from its own broadcast table, so a rank that brackets differently nudges towards a different target on each side of the decomposition - and nothing downstream would say so. |
