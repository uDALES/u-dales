# Adaptive Time Step Limits Test

Drives `modtstep::tstep_limits`, the fused reduction that produces the Courant
and diffusion numbers `tstep_update` chooses `dt` from, and the geometry cache
`inittstep` builds for it.

## Why this test exists

`tstep_update` is the one routine in the time loop whose output is a single
number. Everything else in it is bookkeeping over clocks and counters; the part
worth testing is the reduction over the interior of the domain, and through
`tstep_update` that reduction is visible only as a shift in `dt`, after an
all-reduce, a division by a namelist Courant number and a `min` against two
other limits. Splitting `tstep_limits` out is what makes it observable, and it
is also the split the GPU port needed: the loop is the only part of the module
with a host form and a device form.

## Not modchecksim's pair

`courant_local` and `diffnr_local` in `modchecksim` look like the same two
numbers and are not. Those are diagnostics - computed, printed and discarded.
These choose the time step. They are also a different expression:

| | `modchecksim` | `modtstep` |
| --- | --- | --- |
| velocities | signed | absolute |
| horizontal spacing | `dxhi(i)`, `dy2i` | `dxi`, `dy2i` |
| vertical spacing | `dzhi(k)`, `1/dzh(k)**2` | `1/dzh(k)`, `dzh2i(k)` |
| geometry cache | `diffnrgeom(i,k)` | `diffgeom(k)` |

The two caches agree only on a grid uniform in x, and only to the last bit even
there, which is why `modtstep` has its own rather than borrowing one.

## What the checks pin

| check | what it pins |
| --- | --- |
| the cache | `diffgeom(k)` is exactly `dzh2i(k) + dx2i + dy2i`, at exact equality. The whole justification for precomputing it is that it is the same arithmetic in the same order, so that `dt` does not move by a bit. |
| one term at a time | each of `um`, `vm` and `wm` carries the grid factor it is supposed to carry, with the others zeroed. |
| the absolute values | a strictly negative velocity produces the same number as its positive twin. Without `abs` the reduction, which starts at zero, never leaves it - and a flow with a strong downdraught gets a Courant number several times too small and a `dt` several times too large. |
| the reduction box | a spike at each of its bounds is seen; the same spike one cell outside any of them, in all five fields at once, is not. |
| `ekm` and `ekh` | one at a time with the other zeroed. A run whose Prandtl number is below one is limited by `ekh`, and that is the term a port is most likely to lose. |
| independence | one fused loop still produces two numbers: a velocity must not move the diffusion number, and a diffusivity must not move the Courant number. |
| the time step | it multiplies each maximum once rather than every cell, at **exact** equality, because that is what the hoist claims. |

## The grid substitution

Every case in the repository has `dx = dy`, and all but one a uniform grid in
z. There `dxi` equals `dyi` and `dzh` is constant in `k`, so a loop that pairs
the `v` term with the `x` spacing, or reads `dzh` at a fixed index, returns
exactly the right answer. The test substitutes strictly increasing, exactly
representable values for `dzh`, `diffgeom`, `dxi` and `dyi`, and restores them
before returning.

`diffgeom` is substituted rather than rebuilt from the new `dzh`, so that it
varies in `k` independently: a diffusion term that reached for the Courant
number's spacing instead would otherwise still agree.

The exactness check on the cache runs before the substitution, against the grid
`inittstep` actually saw.

## Coverage

Host branch only. The device kernel is `tests_cuda.f90::test_tstep_limits`,
which runs under `UDALES_RUN_CUDA_SELFTEST` on a Debug GPU build, and every
adaptive parity case in `tests/integration/gpu/case_matrix.json`.

## Running

```bash
tests/integration/tstep/run_test.sh
```

Environment overrides: `UDALES_BUILD`, `CASE_SOURCE`, `NAMELIST_SOURCE`,
`MPI_NAMELIST_SOURCE`, `PROF_STRETCHED`, `NPROCS`, `MPIEXEC`.

## Passes

| pass | what it adds |
| --- | --- |
| uniform vertical grid | the checks, on the grid the cases use |
| stretched vertical grid | a `diffgeom` that varies with `k`, so the exact-equality check on the cache is able to fail |
| two ranks in y | that the reduction runs over the rank's own `jb..je`, not the global extent |
