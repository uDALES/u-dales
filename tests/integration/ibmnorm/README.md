# IBM Normal Velocity Test

Validates the host routines `modibm::ibmnorm` is built from: `solid`, which
pins a field inside the solid, and the two second-order advection corrections
that remove the flux the scheme sent across solid faces.

## What ibmnorm does

Three things, once per Runge–Kutta stage:

```
solid(solid_info_u, um, up, 0., …)          ! velocity and its tendency → 0
solid(solid_info_c, thlm, thlp, interior, …, mask_c)
                                            ! scalars → interior value, or the
                                            !   mean of the fluid neighbours
advecc2nd_corr_{conservative,liberal}(thl0, thlp)
                                            ! undo advection across solid faces
```

The scalar branch is the interesting one. A solid cell with no fluid neighbour
is pinned to a fixed interior value; one with fluid neighbours takes their mean
instead, which makes the field flat across the boundary — a zero-flux
condition.

## What this test pins

1. **Velocity.** A solid point ends at exactly zero, in both the field and the
   tendency, and no other cell moves.
2. **The scalar boundary condition.** A walled-in point takes the interior
   value exactly. A point with fluid neighbours takes their mean, worked out
   here from the mask rather than read back — so a constant field stays
   constant across the boundary, and running twice with different interior
   values leaves those points unchanged.
3. **The liberal correction** vanishes on a constant field, for any velocity,
   and otherwise equals the velocity across each solid face times the jump the
   scheme saw there — one term per face, where the routine forms two.
4. **The conservative correction** on a constant field equals that constant
   times the velocity divergence taken over the solid faces alone.
5. Both corrections vanish at zero velocity, are linear in the scalar, and move
   only cell-centred boundary points.

## Two invariants the parallel port rests on

The device runs one thread per point, so the test also checks what a sequential
loop does not care about:

- **No cell appears twice in a solid list.** A repeat is harmless in a
  sequential loop and a race in a parallel one.
- **The mask agrees with the list**: a cell is masked solid exactly when it is
  listed. `createmasks` builds it that way, by running `solid` over the same
  list. It matters because the masked branch reads its neighbours in place: if
  a listed cell were another's neighbour, the host would read it after it had
  been updated and the device before. The invariant makes that impossible, and
  it is why the device kernel may accumulate in registers.

## Negative controls

Applied to the host routines and confirmed to fail the test:

| Injected fault | Caught by |
| --- | --- |
| `solid` sets the field but not the tendency | solid velocity point not exactly zero |
| `var = var/count` (interior value not removed) | constant field not preserved, interior value leaked, neighbour mean |
| liberal correction substitutes `var(i+1)` where it should use `var(i)` | the liberal closed form |
| conservative correction drops the `k-1` face | the conservative closed form |

Note the third: it survives every *vanishing* check — a constant field makes it
zero either way — and is caught only by the check that pins the value.

Applied to the device kernels (see the scope note below on where they run):

| Injected fault | Caught by |
| --- | --- |
| `solid_device` sets the field but not the tendency | `solid u rhs` |
| `solid_masked_device` keeps the interior value in the mean | `solid masked var` |
| `advecc2nd_corr_liberal_device` drops the `j+1` face | `advecc2nd_corr_liberal` |
| `um_d`/`vm_d`/`wm_d` refresh dropped from `updateDevice` | the `ibm-scalar-wall` parity case, **but only without the self-tests** — `test_masscorr` and `test_ibmnorm` both restore the mirrors on the way out, which hides a missing refresh in a Debug run. The Release nightly, which runs no self-tests, is what catches it. |

The device test needed a second mask to become able to fail. `solid_masked` has
to run against the mask the solid list implies, but the advection corrections
do not, and the case geometries are uniform in `y` — so with that same mask the
`j` faces are never solid and two of the six branches are never taken. The
corrections are compared on a checkerboard instead, which is safe for them
because they only read the mask.

## Scope

This runmode covers the **host** branch. Running it against a GPU build fails
with an explanatory message rather than passing vacuously.

The device kernels are covered by `tests_cuda.f90::test_ibmnorm`, which runs
under `UDALES_RUN_CUDA_SELFTEST=1` on a Debug GPU build. Unlike the wall
function and mass correction self-tests it is a direct comparison against the
host rather than an independent derivation, which is what the other `modibm`
kernels do — the host side is property-tested here, so what is left for the
port is that it reproduces it.

It needs a case with `libm = .true.`, and **no case in the `smoke` selection has
one**, so it first bites in `nightly`/`full` — on `ibm-scalar-wall`,
`ibm-reconstruction` and `surface-energy-balance`.

## A precompute that did not pay

Every one of the six neighbour tests in `solid`'s masked branch and in both
advection corrections reads only the masks, and the masks are finished by the
end of `initibm` and never written again — so each test gives the same answer on
every Runge–Kutta stage of the run. The obvious move is to derive them once into
a flag word per point, plus the reciprocal of the neighbour count, and have both
the host and the device read that instead.

Built and measured, on case 442 (128³, 345 kernel launches):

| kernel | reading the masks | reading a cache |
| --- | --- | --- |
| `solid_masked_device` | 27.6 µs | 28.8 µs |
| `advecc2nd_corr_liberal_device` | 11.8 µs | 13.0 µs |

Slower both times, and by more than run-to-run scatter for the second (median
11.79 vs 12.93 µs, standard deviation 0.3 µs). The reason is that the six mask
reads were never the problem: the points in a list are spatially clustered and
the six neighbours of one cell share cache lines, so those loads almost always
hit. The cache replaces them with two loads that have no reuse at all, one per
point, and costs an extra pass over memory to do it.

The change was reverted. Recorded here so the idea is not re-tried blind: it is
a real precompute, it just loses to the cache hierarchy.

## One thing that is not covered

`updateHost` used to hand `um`, `vm`, `wm`, `thlm` and `qtm` back for `ltrees`
runs, because `vegetation_forcing` read them at the neighbours of every tree
point and one of those can be a solid cell that `ibmnorm` has just pinned on
the device. Removing that handover did **not** fail the `vegetation` case, even
though it has both `ltrees` and `libm`: its trees evidently do not sit against
a solid cell — so the guard rested on inspection rather than on a test.

That is moot now. `vegetation_forcing` runs on the device and reads the
mirrors directly, so there is no handover left to get wrong; the reordering it
compensated for is handled by `ibmnorm` and the reader being on the same side
of the bus.

## Requirements

- A CPU debug build: `tools/build_executable.sh common debug`
- MPI (`mpiexec`), able to run two ranks.

## Running

```bash
tests/integration/ibmnorm/run_test.sh
```

Overridable via environment: `UDALES_BUILD`, `CASE_SOURCE`, `NAMELIST_SOURCE`,
`NAMELIST_SOURCE_MPI`, `NAMELIST`, `NPROCS`, `MPIEXEC`.
