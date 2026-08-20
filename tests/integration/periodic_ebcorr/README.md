# Periodic Energy Balance Correction Test

Validates `modforces::periodicEBcorr`, which removes from the air the heat and
moisture the facets put into it, so a periodic domain does not accumulate
energy without bound.

## What periodicEBcorr does

The facet wall functions sum their flux into `totheatflux` and `totqflux`. Once
per Runge–Kutta stage the correction reduces those across ranks, converts them
to a domain-mean flux, and takes that back out in two pieces:

```
H_proj         = tot_Tflux / (xlen*ylen)          ! domain-mean flux
abl_height     = zh(ke+1) / fraction              ! reverse-engineered ABL depth
R_theta        = H_proj / abl_height
M              = ke - sinkbase                    ! levels above the canopy
R_theta_scaled = R_theta * ke / M

thlp(i,j,k)  += R_theta_scaled                    ! k = sinkbase+1 .. ke
thlp(i,j,ke) += (1-fraction)*tot_Tflux/(xlen*ylen*dzf(ke))
```

with the identical algebra for moisture. Following Grylls (2021), a share
`fraction` of the flux is removed by the volume sink above the canopy and the
remaining `1-fraction` leaves through the top level.

## What this test pins

Checking those three scalars against themselves would only restate the code.
The test asserts the properties they exist to produce:

1. **Column closure.** The depth-weighted column integral of the tendency
   equals the domain-mean flux `tot_Tflux/(xlen*ylen)`, for every `sinkbase`
   and every `fraction` tried. This is what the `ke/M` scaling exists to
   preserve, and it is what breaks when `M` miscounts the levels.
2. **The Grylls split.** Exactly `(1-fraction)` of the total leaves through the
   top level. Both `ke` and `ke-1` carry the volume sink, so their difference
   isolates the top-level term — a quantity the routine never forms.
3. **The reduction.** Each rank contributes a different flux, and every rank is
   checked against the sum. Using the rank-local flux instead of the reduced
   one fails on every rank but the first.
4. **Support.** Nothing at or below `sinkbase` moves; every level above it
   does; nothing outside `ib:ie, jb:je` moves, halo cells included; nothing
   above `ke` moves.
5. **Shape and switches.** The tendency accumulates rather than assigns, is
   horizontally uniform, is linear in the flux (and exactly zero at zero flux),
   and is silenced independently by `lperiodicEBcorr`, `ltempeq` and `lmoist`.

Heat and moisture are also compared against each other at equal fluxes, which
catches a copy-paste slip between the two branches.

## The two passes

| Pass | Ranks | Why |
| --- | --- | --- |
| single rank | 1 | The ordinary configuration. |
| two ranks, unequal fluxes | 2 (`nprocy = 2`) | On one rank the local flux and the reduced flux are the same number, so a dropped `MPI_ALLREDUCE` is invisible. Here the ranks contribute 1 and 2 units and every rank must see 3. |

## A limit of the scheme, not of this test

The column integral is the domain-mean flux only on a **uniform** vertical
grid: the `ke/M` scaling weights levels by count, not by depth. That is a
property of the scheme rather than of the port, so the test checks the grid is
uniform and says so if it is not, instead of quietly assuming it.

## Negative controls

Each of these was applied to the host branch and confirmed to fail the test:

| Injected fault | Caught by |
| --- | --- |
| `tot_Tflux = totheatflux` (no reduction) | column closure, **two-rank pass only** |
| `M = ke - sinkbase + 1` | column closure |
| top-level term dropped | column closure, Grylls split |
| top-level term applied at `ke-1` | Grylls split |
| loops run over `ib-1:ie+1, jb-1:je+1` | halo check |
| `thlp = R_theta_scaled` (assign, not accumulate) | column closure |
| sink loop starts at `sinkbase` | support below `sinkbase`, column closure |

## Scope

This runmode covers the **host** branch. Running it against a GPU build fails
with an explanatory message rather than passing vacuously, because
`periodicEBcorr` writes the device-resident tendencies there and the host
arrays this test inspects would never change.

The device kernels are covered in two other places:

- `tests_cuda.f90::test_periodic_ebcorr` — device self-test, runs under
  `UDALES_RUN_CUDA_SELFTEST=1` on a Debug GPU build. It repeats the same
  closure argument rather than diffing against the host, so a mistake shared by
  both branches cannot pass.
- The `surface-energy-balance` case in `tests/integration/gpu/case_matrix.json`
  — CPU/GPU parity on a real run with `lperiodicEBcorr = .true.`.

## Requirements

- A CPU debug build: `tools/build_executable.sh common debug`
- MPI (`mpiexec`), able to run two ranks.

## Running

```bash
tests/integration/periodic_ebcorr/run_test.sh
```

Overridable via environment: `UDALES_BUILD`, `CASE_SOURCE`, `NAMELIST_SOURCE`,
`NAMELIST_SOURCE_MPI`, `NAMELIST`, `NPROCS`, `MPIEXEC`.
