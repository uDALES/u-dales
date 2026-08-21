# Facet Energy Balance Test

Validates `modEB`: the energy balance solve itself, the flux time integration
in `intqH`, and the longwave exchange in `calclw`.

## Why this one has no device branch

Every other routine in this stretch of the time loop has been ported to the
GPU. `EB` deliberately has not, and the reasons are recorded at the top of
`modEB::EB`:

- It is not field work. Nothing in it is indexed by `(i,j,k)`; the whole
  routine walks facet arrays, and `nfcts` is in the hundreds where the grid is
  in the millions.
- The body runs on rank 0 alone, between two MPI collectives — the all-reduce
  in `intqH` that gathers a facet split across ranks, and the broadcast that
  hands the new temperatures back. MPI here is not CUDA-aware, so a device-side
  solve would have to come down for the reduction and go back up after the
  broadcast, on the one step in `dtEB` where there is any work at all.
- The largest input is the view-factor matrix `vf`, `nfcts` by `nfcts`. On a
  case big enough for the dense `calclw` loop to cost anything, mirroring it
  would be hundreds of megabytes of device memory for a matrix–vector product
  evaluated once per `dtEB`.

Measured on case 064 (140 facets, 5 layers, `dtEB = 10 s`, CPU Release, one
rank, 200 s of simulation): the whole of `EB` is **2.87 ms of a 24.4 s time
loop, 0.012%**. Of that, 2.52 ms is the per-facet matrix solve, 8 µs is
`calclw` on the sparse path, and 0.33 ms is `intqH` over all 600 stages.

So what the "port" of this routine consists of is not moving arithmetic to the
device but making the facet arrays cross the bus only when they carry something
new. That is `updateFacFluxHost` and `updateFacetPropsDevice` in `modcuda`, and
`lfacetprops_dirty` in `initfac`.

## Where EB sits in the time loop

`EB` runs above `updateHost`, immediately after `ibmnorm`:

```
call ibmnorm
  cudaDeviceSynchronize
  integrateFacFluxDevice                        ! every stage, no traffic
  if (eb_will_run()) updateFacIntegralsHost     ! one step in dtEB/dt
call EB
  updateHost                                    ! the field tendencies
call vegetation_forcing
```

**Nothing crosses the bus per step.** `fachf_d` and `facef_d` are integrated
into `fachfi_d`/`facefi_d` on the device and cleared each stage; the totals
come down only when the balance actually fires. That works because the time
integral is per-rank:

```
sum_steps dt * (sum_ranks fachf) == sum_ranks (sum_steps dt * fachf)
```

The two sums commute and `dt` is identical on every rank, so each rank
integrates its own partial flux with no communication at all, and the single
`MPI_ALLREDUCE` moves from twice per step to once per `dtEB`.

Measured with nsys on the parity case over 31 steps at `dtEB = 10 s`:
device-to-host memcpy operations **3,146 → 3,090**, which is exactly the 62
per-step drains (2 arrays × 31 steps) replaced by 6 event drains (2 arrays × 3
energy balances).

Two traps live here, and both are recorded because both are silent rather than
loud:

- **The drain used to sit inside `updateHost`.** Moving `EB` above it without
  moving the drain does not crash — `intqH` reduces a `fachf` that is still
  zero, so the facet heat fluxes are simply absent from the energy balance on
  the GPU and present everywhere else. Measured, `facEB:hf` came out as 0.0
  against a CPU reference of −103.75 W/m², carrying a 103.7 error into `dTdz`.
- **The clear is not optional even though the copy is.** The wall functions
  add into the accumulators every stage; let them run on and the third stage
  integrates the sum of all three, which is the facet heat flux three times
  too large. That bug shipped once.

`EB` now depends on `updateHost` not at all, so this arrangement survives
unchanged when `vegetation_forcing` is ported and `updateHost` goes away.

## What this test pins

1. **`intqH`'s stage contract, and that it communicates nothing.** The facet
   fluxes are integrated on the third Runge–Kutta stage only and cleared on
   every stage. The integral is *local*, so the probe is rank-dependent and
   every rank checks its own number — a routine that reduced here would leave
   the same total on all of them, which the two-rank pass catches on both.
2. **The reduction, once, when the balance fires.** Every rank is given a
   distinct local integral and the result is checked on both sides: normalised
   to a flux on rank 0, where the division by `tEB` and `faca` happens, and as
   the raw reduced total elsewhere — which also pins that the reduction is an
   all-reduce and that the normalisation is rank 0's alone. The local
   integrals must come back reset, or the next interval double-counts.
3. **When the energy balance fires, and what it tells the GPU mirror.** Not on
   the first two stages however late it is; not on the third until `tnextEB`
   has arrived; and when it does fire it sets `lfacetprops_dirty`. That flag is
   the other half of the contract `tests_cuda.f90::test_facet_props_refresh`
   checks from the device side.
4. **Equilibrium is a fixed point.** Uniform layer temperatures with the net
   surface flux exactly balancing the emitted longwave must leave the facet
   where it was. Nothing conducts and nothing radiates on net, so the scheme
   has to return its input — and almost any slip in the matrix assembly, the
   inverse or the two matmuls destroys that.
5. **The surface energy balance closes.** Row one of the system says

   ```
   faclam(n,1)*facTdash(n,1) = boltz*facem(n)*Told^3*Tnew - (SW + LWin + H + E)
   ```

   which the test checks against `faclam`, the returned `facTdash`, and the
   four flux terms it set itself. That is the physics rather than a second copy
   of the algebra: it never builds A, B, C, D or E.
6. **The innermost layer is held.** The last row of the system exists to pin
   it, so `facT(n,nfaclyrs+1)` must not move.
7. **`calclw` agrees with itself.** The sparse and dense view-factor paths are
   two separate loops that must produce the same incoming longwave. The test
   builds the dense matrix from the sparse triplets and runs both.

## Negative controls

Every check has one, and each control asserts that it fires — a control that
cannot fail says nothing about the check beside it.

- The stage-count and reduction probes are checked to be distinguishable from
  what a broken routine would produce, so a silently zero probe cannot make
  phase 1 pass for free.
- A 1% flux imbalance must move a facet that the balanced case left still.
- The surface closure is re-evaluated with the **new** surface temperature in
  the linearisation instead of the old one, and again with the sensible heat
  term dropped. Both must break it. The first one also proves the surface moved
  far enough for `Told` and `Tnew` to be told apart at all — without enough
  heating that control quietly passes for the wrong reason.
- The largest view factor is removed from the dense matrix, which must change
  some facet's incoming longwave. Otherwise the sparse/dense comparison is
  between two arrays that never depended on the view factors.

## Mutation testing

The checks were run against fourteen deliberate faults in `modEB` when the
balance was first written, and against twelve more when the flux integration
moved to the device. All twenty-six are caught. The original fourteen:

| # | Fault | Caught by |
|---|---|---|
| 1 | `intqH` clears the fluxes only on the consumed stage | flux left behind on stage 1 |
| 2 | `intqH` integrates on every stage | `fachfi` after stage 1 |
| 3 | `intqH` integrates the rank-local flux | `fachfi` after stage 3, two ranks |
| 4 | `EB` never sets `lfacetprops_dirty` | dirty flag not set |
| 5 | `EB` sets it on every call | dirty flag set on stage 1 |
| 6 | Linearisation uses `faclam(n,2)` | fixed point drifts |
| 7 | Sensible heat dropped from the net flux | surface closure |
| 8 | Innermost layer no longer pinned | innermost layer moved |
| 9 | Sparse longwave reads `facT(i,1)` for `facT(j,1)` | sparse vs dense |
| 10 | `BM` assembled from `facd(n,1)` for every layer | conduction relation |
| 11 | `CM` uses `faccp(n,1)` for every layer | energy conservation |
| 12 | `DM` gradient correction halved | energy conservation |
| 13 | `EM` uses `faclam(n,m)` on both faces | energy conservation |
| 14 | Forcing loses half its `tEB` scaling | fixed point drifts |

The twelve added with the per-rank integral — eight against the host contract
and four against the device one:

| # | Fault | Caught by |
|---|---|---|
| 15 | `intqH` integrates on every stage | `fachfi` after stage 1 |
| 16 | `intqH` stops clearing the accumulators | flux left behind on stage 1 |
| 17 | `intqH` integrates without `dt` | `fachfi` after stage 3 |
| 18 | `EB` does not reset the per-rank integrals | integrals not reset |
| 19 | `EB` reduces the wrong array | reduced heat flux wrong |
| 20 | `EB` normalises the raw integral, not the reduced one | reduced heat flux wrong |
| 21 | `eb_will_run` ignores `tnextEB` | ran before its time |
| 22 | `eb_will_run` ignores the stage | moved the surface on stage 1 |
| 23 | Device kernel integrates on every stage | `facEB` parity |
| 24 | Device kernel stops clearing the accumulators | `facEB` parity |
| 25 | Download leaves the device integral in place | `facEB` parity |
| 26 | Download runs every step, not just on EB events | `facEB` parity |

Faults 10 to 13 are why the test imposes a layered wall profile before it
starts. **Every wall type in the shipped case has identical layers** — same
thickness, same heat capacity, same conductivity all the way through — so a
routine that reads the wrong layer produces exactly the right answer. Fault 6
survived the whole file until the profile was added, and 10 to 13 survived
until the conduction relation and the energy identity were added on top of it.
The test now asserts the profile is not flat, so if that ever stops holding it
says so rather than passing on a degenerate wall.

## Running

```bash
bash tests/integration/eb/run_test.sh
```

Two passes, one rank and then two. The two-rank pass is what separates the
reduced facet flux from the rank-local one; on one rank they coincide and the
test says so rather than pretending to have checked.

CPU build by default. Nothing in the routine is device-side, so it runs on
either.

## A note on `nnz`

This work turned up a silent defect in case 064 itself, now fixed.

`nnz` is the number of entries read from `vfsparse.inp.xxx`, and nothing else
bounds that loop. It matters only under `lvfsparse` — but case 064 is the one
case in the tree that sets `lvfsparse = .true.`, and it omitted `nnz`. The
variable was an uninitialised module integer, so in practice it was zero: the
file was opened, none of its 992 entries were read, and the facets could not
see each other. Only the sky term reached `facLWin`.

That is not a rounding detail. Sky view factors in this case go down to 0.53,
and measured on a 31 s run the worst facet's incoming longwave was **179.5
W/m² instead of 336.3 — short by 157 W/m², 47% of its hemisphere**. Surface
temperatures differ by only 0.01 K over so short a run, because the wall has
real thermal inertia, which is exactly why nothing noticed.

Three things changed:

- `nnz = 992` is now set in `namoptions.064`, `namoptions.064.seb` and
  `namoptions.1012`. The `sparse-view-factors` coverage the parity case
  declares is finally exercised.
- `nnz` defaults to 0 in `modglobal`, so the value is defined rather than
  relying on static zero-initialisation.
- `checkinitvalues` rejects `lEB .and. lvfsparse .and. nnz <= 0` at startup
  with a message naming the file, rather than letting the physics go missing.

Both preprocessors already write `nnz` alongside the sparse file
(`tools/write_inputs.m`, `tools/python/udprep/udprep_radiation.py`), so a
namelist without it is stale rather than deliberate.
