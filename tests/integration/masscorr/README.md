# Mass Flow Rate Correction Test

Validates `modforces::masscorr`, which corrects the predicted velocity so the
domain meets a prescribed flow rate.

## What masscorr does

Four branches, one per controller: volume-averaged and outlet-plane, for `u`
and for `v`. Each reduces the current and predicted velocity over the fluid
cells, forms the shortfall against the target, and adds it uniformly to the
tendency. For the volume branch in `u`:

```
uvol(k)      = <up>_xy at level k, fluid cells only     (avexy_ibm)
uvolold(k)   = <um>_xy at level k, fluid cells only
uoutflow     = rk3coef * sum_k uvol(k)*dzf(k) / zh(ke+1)
uflowrateold =           sum_k uvolold(k)*dzf(k) / zh(ke+1)
udef         = uflowrate - (uoutflow + uflowrateold)

up(i,j,k)   += udef * rk3coefi                          ! ib:ie, jb:je, kb:ke
```

The outlet branch does the same over the plane `i = ie`, divided by the outlet
area rather than the domain height.

## What this test pins

Checking `uoutflow`, `uflowrateold` and `udef` against their own expressions
would restate the code. The test instead measures, with its own reduction, the
rate the routine claims to set:

1. **The target is met.** After the call, the fluid-volume mean of
   `um + rk3coef*up` equals `uflowrate` — and the outlet-plane mean does too,
   for the outflow branch. Both directions, all four branches.
2. **A second call is a no-op.** The target is already met, so `udef` must come
   out zero and the tendency must not move. A defect that is right once but
   inconsistent with the correction actually applied fails this even though it
   passes check 1.
3. **The rate is the global one.** Each rank is given a different slab of one
   global field, and every rank is checked against the same answer, so a
   reduction that stayed rank-local fails on every rank.
4. **The two directions are independent.** The `u` branches must not move `vp`,
   and the `v` branches must not move `up`.
5. **Shape and support.** The correction is one and the same constant over
   `ib:ie, jb:je, kb:ke`, is non-zero, and reaches nothing outside that — no
   halo cell, no level above `ke`.
6. **Switches.** All four controllers off, or `linoutflow` on with all four
   set, is a genuine no-op.

## The two mask phases

The volume branches average over the fluid cells, so the first phase runs
against the case's real IBM masks. The outflow branches divide by an outlet
area that `calcfluidvolumes` derives from `IIc` while the branch itself masks
with `IIu` or `IIv`; those agree only when the outlet plane is clear, so the
second phase sets every mask to one and recomputes the areas. That is the
configuration those branches are written for — *"Assumes ie=itot"*, with no
mention of blocks at the outlet.

The device self-test punches its own mask instead of using the case's. Every
GPU smoke case runs with `libm = .false.`, where every mask is one and a device
reduction that ignored the mask would agree with a correct one; the check has
to make its own solid cells to have anything to say. This was found by a
negative control that failed to fire.

## The two passes

| Pass | Ranks | Why |
| --- | --- | --- |
| single rank | 1 | The ordinary configuration; the only one where the `v` outlet row is whole. |
| two ranks, split domain | 2 (`nprocy = 2`) | On one rank the local reduction and the global one are the same number, so a dropped `MPI_ALLREDUCE` is invisible. |

The `v` outflow branch is skipped at two ranks. Its outlet is a row of constant
`j`, which lives on one rank once the domain is split in `y`, while the branch
all-reduces over every rank and would add up one row per rank. That is a
property of the scheme, not of this port, so the check runs only where the row
is whole. The `u` outflow branch has no such problem under a `y` split: its
outlet is a plane of constant `i` that every rank holds part of, and the
all-reduce assembles it.

## Negative controls

Applied to the host branch and confirmed to fail the test:

| Injected fault | Caught by |
| --- | --- |
| `total = local` in `sum_ibm_reduce` (no reduction) | outflow rate, `uouttot`, **two-rank pass only** |
| `udef = uflowrate - uoutflow` (drop the current rate) | volume and outflow rates |
| `up += udef` instead of `udef*rk3coefi` | rates, and the second-call no-op |
| correction loop runs to `ke+kh` | support outside the interior |
| `averl(k) = sum(var)` with no mask in `avexy_ibm` | volume rates, second-call no-op |
| the original `sumy_ibm(vout, vp(ib:je,…))` v outlet call | aborts: `munmap_chunk(): invalid pointer` |

Applied to the device data path:

| Injected fault | Caught by |
| --- | --- |
| `um_d = um` / `vm_d = vm` dropped from `updateDevice` | the `surface-energy-balance` and `ibm-scalar-wall` parity cases |

Applied to the device kernels:

| Injected fault | Caught by |
| --- | --- |
| `add_uniform_device` loops over `jb-jh:je+jh` | support outside the interior |
| `modcuda::avexy_ibm_device` drops `*II` | masked volume rate |
| `modcuda::sumx_ibm_column_device` scaled by 0.5 | `v` outflow rate |
| `avexy_ibm_device` declares `var` from `kb`, not `kzb` | volume rate (reads `um_d` shifted by `kh`) |
| `modcuda::sumy_ibm_column_device` scaled by 1.001 | `u` outflow rate |
| `col_stage = col_d` dropped (stale staging) | masked volume rate |

## Why `updateDevice` refreshes `um_d` and `vm_d`

Running on the device moved `masscorr` ahead of `updateHost`, and there it
reduces `um_d` and `vm_d` — which no update routine had written by that point
in the stage. `updateDevicePriorPoiss` refreshes them, but only later in the
same stage, so on the first stage of a run they held whatever the allocation
left behind. Case 064 diverged on step one with a Courant number of 2×10⁷
against a divergence of 2×10⁻⁸: a uniform, divergence-free velocity, which is
exactly what a wrong `udef` adds.

Refreshing them once at startup would not be enough. The host is the last
writer before the stage: `ibmnorm` zeroes `um` at solid points, and `boundary`
rewrites it afterwards — for a profile or driver inlet at the interior plane
`i = ib`, not only in the halo. So `updateDevice` copies them, guarded by the
controllers that actually read them, and a run without flow forcing pays
nothing.

## A side effect the self-test has to undo

`masscorr` leaves `udef`, `vdef` and `uouttot` behind in `modfields`, and
`uouttot` is read by the outflow boundary conditions during the first step —
`boundary` only recomputes it at the *end* of that step, after the damage would
be done. The device self-test runs between `initCUDA` and the time loop, so a
value it left in `uouttot` perturbed the first step of every inlet/outlet case.
It showed up as `boundary-x-profile` failing parity while `masscorr` was inert
in that case. The self-test now saves and restores all three.

## Scope

This runmode covers the **host** branch. Running it against a GPU build fails
with an explanatory message rather than passing vacuously, because `masscorr`
writes the device-resident tendencies there and the host arrays this test
inspects would never change.

The device kernels are covered in two other places:

- `tests_cuda.f90::test_masscorr` — device self-test, runs under
  `UDALES_RUN_CUDA_SELFTEST=1` on a Debug GPU build. It repeats the same rate
  argument rather than diffing against the host, so a mistake shared by both
  branches cannot pass. It allocates the mask mirrors and the staging column
  that `initCUDA` leaves out when a controller is off, so every case exercises
  all four branches whatever its namelist asks for, and gives back exactly what
  it took.
- The `surface-energy-balance` and `ibm-scalar-wall` cases in
  `tests/integration/gpu/case_matrix.json` — CPU/GPU parity on real runs with
  `luvolflowr = .true.`.

## Requirements

- A CPU debug build: `tools/build_executable.sh common debug`
- MPI (`mpiexec`), able to run two ranks.

## Running

```bash
tests/integration/masscorr/run_test.sh
```

Overridable via environment: `UDALES_BUILD`, `CASE_SOURCE`, `NAMELIST_SOURCE`,
`NAMELIST_SOURCE_MPI`, `NAMELIST`, `NPROCS`, `MPIEXEC`.
