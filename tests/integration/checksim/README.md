# Checksim Test

Validates the three rank-local reductions behind `modchecksim::checksim` —
`courant_local`, `diffnr_local` and `div_local` — and the diffusion-number
geometry cache `diffnrgeom` that `initchecksim` builds for them.

## What checksim does

Once every `tcheck` seconds (every step when `tcheck = 0`) it prints three
numbers, each a reduction over the cells the rank owns followed by an
all-reduce:

```
courant   max over ib..ie, jb..je, kb..ke of  (um*dxhi(i) + vm*dyi + wm*dzhi(k)) * dtmn
diffnr    max over the same box of  max(ekm, ekh) * (1/dzh(k)^2 + dxh2i(i) + dy2i) * dtmn
div       max |div| and the volume integral of div, where
          div = (u0(i+1)-u0(i))*dxi + (v0(j+1)-v0(j))*dyi + (w0(k+1)-w0(k))*dzfi(k)
```

Two things in there depend only on the grid, so they are built once and both
the host loop and the device kernel read them from there. `initchecksim` builds
`diffnrgeom(ib:ie, kb:ke)` for the geometry factor in the middle line, since
nothing else wants it. The cell volume `dx*dy*dzf(k)` in the last line is a
generic grid quantity, so it lives in `modglobal` as `dvcell`, beside the
reciprocal `dxdydzfi` that the tendency routines divide by, and is mirrored to
`dvcell_d` in `modcuda` alongside `dxdydzfi_d`.

`dtmn` multiplies the maximum rather than every cell in the first two lines,
which is exact: rounding is monotonic and `dtmn` is a mean time step, so
`max(fl(a*dtmn), fl(b*dtmn))` is `fl(max(a,b)*dtmn)`.

## What this test pins

Each reduction is driven with a state whose answer follows from the grid rather
than from re-running the loop, so a transcription error in the port cannot be
reproduced by the expectation.

1. **The caches are built from the right quantities.** `diffnrgeom` is
   asserted equal to the expression it replaced with `/=`, not a tolerance:
   the point of precomputing it was that the printed diffusion number must not
   move, and only exact equality says that. Rebuilding it as `dzhi(k)*dzhi(k)`
   instead of `1/dzh(k)**2` — one rounding instead of two, a difference well
   below print precision — fails this. `dvcell` is checked the same way against
   `dx*dy*dzf(k)`; it is *not* claimed bit-identical to the `div*dx*dy*dzf(k)`
   it replaces, but it must be the product of exactly those three, so a cache
   built from `dzh` or from `dzfi` is caught. Only the stretched pass can see
   the `dzh` version, since on a uniform grid `dzf` and `dzh` are equal. The
   device mirror is checked in `tests_cuda.f90::test_cell_volume_reciprocals`,
   with `dxdydzfi_d` and `dxdydzhi_d`, rather than here.
2. **One term at a time.** A single spike in one field, with everything else
   zero, isolates one term of a sum, so the maximum is that term alone and can
   be written down. `ekm` and `ekh` are checked separately, because a run whose
   Prandtl number is below one is limited by `ekh` and that is the term a port
   is most likely to lose.
3. **The reduction box.** The same spike is planted at each bound of
   `ib..ie, jb..je, kb..ke` in turn, and a large poison value just outside each
   of those bounds must be ignored. A port that reduces over the halo still
   reproduces every interior value and would pass an interior-only comparison.
4. **The forward reach.** `div_local` reads `i+1`, `j+1` and `k+1`, so its
   spikes go at `ie+1`, `je+1` and `ke+1`: the cell that sees them is the last
   interior one. A port written `u0(i) - u0(i-1)`, or one stopping at `ie-1`,
   reads zero there and returns nothing.
5. **The volume integral covers the whole box.** A ramp in x gives every cell
   the same divergence, so `divtot` is a cell count times a known value. That
   is the check that fails when the sum is one plane short, which no
   single-spike test can see.
6. **Every grid factor and every index.** After the cache-exactness check, the
   test substitutes strictly increasing values for `dxhi`, `dzhi`, `dzf`,
   `dzfi`, `diffnrgeom` and `dvcell`, and distinct values for `dxi`, `dyi`,
   `dx` and `dy`, restoring them before it returns. Every case in the repository has
   `dx = dy` and a uniform grid in x, so without this `dxhi(i)` equals `dyi`
   and equals `dxhi(ib)`, and a loop that pairs the `v` term with the x spacing
   or reads the array at a fixed index returns exactly the right answer. `dzfi`
   is deliberately not the reciprocal of `dzf`, so a kernel that builds the
   gradient from the wrong one of the pair is caught. `dvcell` is substituted
   too, and `dx` and `dy` are given values unrelated to it, so a `div_local`
   that reached for `dx*dy*dzf(k)` directly instead of the cache would not
   agree with the expectation.

   A cyclic substitution is not enough: with a period of seven on a 64-wide
   domain, `dxhi(ib)` and `dxhi(ie)` come out equal and a kernel pinned to `ib`
   passes. The values increase strictly instead, so no two indices collide.

## Passes

The runner makes three:

| pass | what it adds |
| --- | --- |
| uniform vertical grid | the case as it ships |
| stretched vertical grid | `prof.inp.stretched`, so `dzh`, `dzf` and their reciprocals differ at every level and `initchecksim` builds a cache that actually varies |
| two ranks in y | `jb..je` is no longer the whole domain, so a loop written over the global extent is distinguishable from one written over the rank's own |

`prof.inp.stretched` is a copy of the file of the same name in
`tests/integration/ibm_wallfun`; it is a fixed geometric grid at ratio 1.02
over the 64 levels case 101 uses.

## Scope

This runmode covers the **host** branch. Running it against a GPU build fails
with an explanatory message rather than passing vacuously, because the three
reductions read the device fields there.

The device kernels are covered by `tests_cuda.f90::test_checksim_reductions`,
which runs under `UDALES_RUN_CUDA_SELFTEST=1` on a Debug GPU build. It builds
its reference with an explicit host loop over the same seed — not by calling
the host branch, which is not compiled in that build at all — and it makes the
same grid substitutions, since no GPU parity case has a stretched or
anisotropic grid either.

Nothing else covers these routines: `checksim` only prints, so a CPU/GPU parity
run compares none of its output.

## Requirements

- A CPU debug build: `tools/build_executable.sh common debug`
- MPI (`mpiexec`). The third pass needs two ranks.

## Running

```bash
tests/integration/checksim/run_test.sh
```

Overridable via environment: `UDALES_BUILD`, `CASE_SOURCE`, `NAMELIST_SOURCE`,
`MPI_NAMELIST_SOURCE`, `PROF_STRETCHED`, `NPROCS`, `MPIEXEC`.

## Implementation

- Runmode `1014` (`TEST_CHECKSIM` in `src/modglobal.f90`)
- `tests_checksim` in `src/tests.f90`, dispatched from
  `execute_runmode_actions` in `src/program.f90`

## See also

`tests/integration/nudge` follows the same runmode-test pattern.
