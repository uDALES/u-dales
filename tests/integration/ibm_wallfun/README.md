# TEST_IBM_WALLFUN (runmode 1008)

Host-side validation of the IBM wall functions in `src/modibm.f90`:
`diffc_corr`, `wallfunmom` and `wallfunheat`.

## Why this exists separately from the device self-tests

`src/tests_cuda.f90` compares the OpenACC port against these same host
routines. That comparison is valuable but has two limits:

- it is symmetric, so a mistake present on both sides passes it;
- it needs a GPU, so on a CPU-only machine these routines have no coverage.

This runmode instead checks properties the host routines must hold on their
own, derived without reference to how they are written.

## What is checked

| # | Property |
|---|---|
| 1 | With no solid neighbours there is nothing to cancel, so no correction is applied. Fails if a mask test is written the wrong way round. |
| 2 | A constant field has zero differences, so the correction is zero whatever the coefficients are. |
| 3 | One isolated solid face, uniform diffusivity and a unit jump gives `-ekh*dx2i`, worked out by hand. Pins the sign and the coefficient. |
| 4 | Running the correction twice doubles it: it accumulates into the tendency rather than replacing it. |
| 5 | Only cells in the boundary point list are written. |
| 6 | A fluid at rest produces no stress and reports none. |
| 7 | The wall stress opposes the flow, and the momentum removed from the field equals the momentum reported on the facets. The second is a cross-check between two outputs written under different indexing. |
| 8 | The heat removed from `thlp` equals `totheatflux`, and `fac_pres` is recomputed here directly from the section list. |
| 9 | The moisture wall function: zero flux at equilibrium, clamped when the air is wetter than the surface, moisture only ever added when it is drier, and nothing at all without a green roof. |
| 10 | `dxdydzfi` and `dxdydzhi` are each the reciprocal of their own cell volume. Asserted as a definition rather than by recomputing the same product, so the check and the code cannot share a mistake. |
| 11 | `local_coords` returns an orthonormal frame, orthogonal to the facet normal, unchanged when the velocity is scaled, and invalid for a wall-normal velocity. |

## Two vertical grids

The runmode is executed twice: once on the case's own uniform vertical grid,
and once with `prof.inp.stretched`, a geometrically stretched grid of the same
depth.

The stretched pass is not redundant. On a uniform grid `dzf` and `dzh` are equal
at every level, so `dxdydzfi` and `dxdydzhi` hold identical values - a swap
between them, or either built from the wrong spacing, changes nothing. Every
other case in the repository is uniform as well, the CPU-GPU parity matrix
included, so on a uniform grid alone nothing anywhere would catch it. Checks 7,
8 and 10 can only fail on the stretched pass.

## Running

```bash
bash tests/integration/ibm_wallfun/run_test.sh
```

Runs as part of `python tests/run_tests.py solver-suites`, so it is covered by
the `supported` stream on both Linux and macOS. Needs a case with
`libm = .true.` and facet sections on the rank; the bundled namelist uses
`tests/cases/101`.

## Verified against

Each check was confirmed to fail when the corresponding behaviour is broken:

| Break | Caught by |
|---|---|
| mask polarity inverted (`< eps1` to `> eps1`) | 1 and 3 |
| correction assigns instead of accumulating | 4 |
| `rhs - momvol` becomes `rhs + momvol` | 7, both halves |
| heat scatter scaled by 1.001 | 8 |
| facet pressure accumulated without the area weight | 8 |
| `dxdydzfi` built from `dzh` | 7 and 10, stretched grid only |
| `local_coords` normalisation removed | 11 |

The two entries marked "stretched grid only" pass on the uniform grid. That is
the reason the second pass exists.
