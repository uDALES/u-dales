# SGS TKE top boundary condition

Regression test for the fix to
[#373](https://github.com/uDALES/u-dales/issues/373): the subgrid TKE (`e12`)
condition at the top of the domain must match the momentum condition in the
same `BCtopm` branch.

| `BCtopm` | momentum | SGS TKE |
|---|---|---|
| `freeslip` (1) | zero flux | zero flux — `e120(:,:,ke+1) = e120(:,:,ke)` |
| `noslip` (2) | Dirichlet, `Uinf` | Dirichlet, `e12min` |
| `pressure` (3) | zero flux | zero flux |

In `diffe` (`modsubgrid.f90`) the flux through the top face is proportional to
`e120(i,j,ke+1) - e120(i,j,ke)`, so "zero flux at the lid" is exactly "the
ghost equals the top interior cell". That is what the test asserts, by reading
`e120` back out of the restart file.

## Why it was wrong

`freeslip` and `pressure` used to pin the ghost to `e12min = 5e-5`, which is
the numerical floor the one-equation model is kept above, not a physical
boundary value. Since the interior TKE is always well above the floor, the lid
was a permanent subgrid-energy sink. Measured on the tree case `526` when the
fix was made, the top interior cell was depressed by **35 %** (`1.35e-2`
against `2.09e-2`), with the influence still visible about six cells down.

`noslip` set no TKE condition at all, so the ghost held whatever it last
contained.

## Running it

```bash
UDALES_BUILD=$PWD/build/release/u-dales \
  python tests/integration/tke_top_bc/test_tke_top_bc.py
```

Three runs of case `300` (one per `BCtopm`), 2 s of simulated time each on 4
ranks, about 12 s in total with a Release build.

Case `300` is the 16x16x32 stretched-grid case without buildings. The lid
condition is a plane operation on `e120(:,:,ke+1)` in `modboundary.f90` and
nothing about it depends on the domain size or on immersed boundaries, so the
smallest case that runs the one-equation model is the right one; the stretched
grid additionally exercises `dzf(ke+1)` and `dzh2i(ke+1)` in the top-face term
of `diffe`. The condition is applied on every step, so a handful of steps is
enough: the assertions hold from the first one. The test adds the `&BC` group
to the case itself, because `namoptions.300` carries none.

Against the pre-fix binary three of the four tests fail: `freeslip` and
`pressure` by `9.5e-2` — thirteen orders above the tolerance and comparable to
the TKE itself — and `noslip` by exactly `e12min`, the ghost having been left
at zero.

## Two things the test deliberately handles

**`ldelta = .true.` is set.** With the stability-dependent length scale the
one-equation closure evaluates `cn*e120/sqrt(grav/thvs*dthvdz)`
(`modsubgrid.f90`), which raised a floating-point exception before the first
timestep on the tree case this test was first written against. That is a
separate defect; `ldelta` takes the other branch and keeps the length-scale
choice out of the test.

**The initial `e12` profile is written as 0.1.** The test rewrites the TKE
column of `prof.inp` itself rather than relying on what the case ships with. A
zero column is clamped to `e12min` by `modstartup`, and left there the
zero-flux and Dirichlet conditions would be numerically indistinguishable and
every assertion would pass regardless of what the code did.
`test_the_lid_carries_real_turbulence` guards that: it requires the top cell to
carry TKE at least ten times the floor, and reports the actual ratio (~1800).
