# Driver Plane Test

Validates `moddriver::drivergen` with `idriver = 2`: that one call leaves all
twelve driver inlet planes current, in each of the four time branches it can
take, and that the `m` planes are written on `rk3step` 0 and 3 only.

## What drivergen does

A driven simulation reads a precursor's recycle plane from disk at startup
(`readdriverfile`) into `store*driver(j, k, record)`, and then on every stage
interpolates the record set to the current `timee` into twelve planes:

```
u0driver   umdriver
v0driver   vmdriver
w0driver   wmdriver
thl0driver thlmdriver     (ltempeq .and. lhdriver)
qt0driver  qtmdriver      (lmoist  .and. lqdriver)
sv0driver  svmdriver      (nsv > 0 .and. lsdriver)
```

`modboundary::xmi_driver` and its temperature, moisture and scalar siblings
read all twelve to set the inlet ghost cells. Four branches pick the records:

| condition | records used |
| --- | --- |
| `abs(elapsrec) < 1e-4` | the nearest record, exactly |
| `elapsrec > 0 .and. x == 1` | record 1, clamped - before the first record |
| `elapsrec < 0` | `x` and `x+1` |
| `elapsrec > 0` | `x-1` and `x` |

`u` and `v` are then rotated by the inflow angle `iangle`, and only after that
are the `m` planes taken from the `0` planes - and only on `rk3step` 0 or 3.

## Why this is tested on the host

This is the host half of the invariant the GPU driver inlet rests on.
`initCUDA` mirrors the twelve planes onto the device and fills them there and
then, rather than leaving the fill to the time loop. That is only correct if
`drivergen` really does write all twelve in the one call the initialisation
sequence makes - which is what the passes here pin down.

It has to be that way round because `boundary_device` uploads the planes from
`updateDriverPlanesDevice` only on the stage `drivergen` itself runs on, the
last of the three, so nothing uploads a complete set until the end of the
first step. Allocating the device planes without filling them left the first
two stages of the first step reading whatever the allocation happened to land
on. `updateDevice` uploads `u0driver_d` alone, for `bcpup`, so `u0` was the
one inlet field that still came out right - which is what made the failure
hard to read. Case 452 came out of its first step with a divergence of 4.33
against 2.7E-13 for the same case on the host, a 40 K temperature anomaly at
the inlet, and never recovered.

The final pass is that fact stated directly: at `rk3step` 1 the `0` planes
move and the `m` planes do not.

## What the test asserts

Per pass, three separate claims:

1. every `0` plane holds the interpolated record,
2. every `m` plane holds exactly what its `0` plane holds,
3. no plane was left at the sentinel it was filled with beforehand - which is
   what lets a plane nobody thought to check individually still fail.

The store is built in the test rather than read from a driver file. Every
record is linear in the record index and every time used is a dyadic fraction,
so the interpolated answer is exact and follows from the interpolation
arithmetic rather than from a fixture that would have to be regenerated with
the code it checks. The `m` store is filled with a poison value, so sourcing
the `m` planes from the file instead of from the `0` planes fails claim 2.

The rotated pass runs with `iangle` non-zero. `u` and `v` are rotated after
the interpolation and before the `m` planes are taken, so a copy made ahead of
the rotation leaves `um` and `u0` disagreeing.

## Coverage

This covers the host branch only. The device side of the same invariant - that
the mirror is live when `initCUDA` returns, and that
`allocDriverPlanesDevice` fills what it allocates - is
`tests_cuda.f90::test_driver_plane_seed`, which runs under
`UDALES_RUN_CUDA_SELFTEST` on a Debug GPU build.

`drivergen` has no device branch, so unlike the checksim runmode this one runs
against a GPU build as well.

## Running

```bash
tests/integration/driver_planes/run_test.sh
```

Environment overrides: `UDALES_BUILD`, `CASE_SOURCE`, `NAMELIST_SOURCE`,
`MPI_NAMELIST_SOURCE`, `NPROCS`, `MPIEXEC`.

## Passes

| pass | what it adds |
| --- | --- |
| one rank | the four time branches, the rotation, and the `rk3step` gate |
| two ranks in y | that the planes are indexed by the rank's own `jb..je`; the expected values carry `zstart(2)` |
