# Thermodynamics Contract Test

Pins the facts about `modthermodynamics::thermodynamics` that the GPU port
reproduces rather than corrects.

## Why this test exists

The GPU port of thermodynamics has its own test:
`tests_cuda.f90::test_thermodynamics_device` seeds one set of fields, runs the
host routine and the device routine over it, and requires everything they
produce to agree - six 3D fields and thirteen profiles. That catches any
disagreement between the two implementations, and nothing that both get wrong
together.

This test is that blind spot. Each check below is a fact about the host
routine that is odd enough that someone would reasonably change it, and where
the device kernel is written to match rather than to be right.

## The saturation level shift

`thermo` writes `ql0` one k level low.

```fortran
real, intent(in)  :: qt (ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh), &
                     thl(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh), ...
real, intent(out) :: ql (ib-ih:ie+ih, jb-jh:je+jh, kb   :ke+kh)
```

The `ql` dummy is a k level shorter than the other two - the commented-out
declarations above the routine show `qt` and `thl` being widened to `kb-kh`
and `ql` being left behind. The first of the two calls passes the whole of
`ql0`, whose own lower bound is `kb-kh`. Explicit shape means sequence
association, so dummy level `kb` lands on actual level `kb-kh`: the saturation
computed from level `k` is stored at level `k-kh`, and `ql0(:,:,ke+kh)` is
never written at all.

`ql0h` is not affected. Its lower bound is `kb`, so the second call is
aligned, and the test asserts that too - both halves of the same fact.

The shift reaches `ql0av`, `calthv` and the `ql` field in a fielddump, so it
is not cosmetic. Correcting it is a physics change to the CPU solver and
belongs in its own commit. Until then `thermo_device` declares its dummies
character for character the same way, so the two shift together and stay
comparable. **If this test starts failing because the host was corrected, the
device declaration has to be corrected in the same commit.**

## The other four

| check | what it pins |
| --- | --- |
| `dim(x,0.)` is `max(x,0.)` | `thermo_device` spells the clamp as `max`, because no device runtime is guaranteed to implement `DIM`. Written down rather than assumed. |
| half levels | `calc_halflev` interpolates by cell height and then overrides level `kb` with `thls` and `qts`, after the loop that already wrote it. |
| the `dthvdz` clamp | `calthv` writes a hard zero at `kb` and the clamp then turns it into `+eps1`, because `sign(eps1, 0.)` is `+eps1`. The clamp runs over `kb..ke`, so `ke+kh` keeps the whole-array zero. subgrid divides by `dthvdz`, so a level left at zero would be a division by zero rather than a small number. |
| the two `thv0h` halos | The saturated branch writes `ib:ie, jb:je` and leaves the halo columns alone. The unsaturated branch assigns `thv0h = thl0h(:,:,kb:ke+kh)` whole, halo included. The device has two kernels for the same reason. |

The unsaturated branch is unreachable from a moist namelist, so the test flips
`lmoist` and calls `thermodynamics` a second time to get at it.

## How the reference is built

The seeded profile has a real vertical gradient, so the finite differences
`calthv` forms are nowhere near zero - which matters, because its `chi_sat`
denominator divides by one of them. `qt0` alternates by level between one
value no saturation specific humidity reaches and one every plausible value
sits below, so both arms of the saturated test are taken and neither sits near
its threshold.

`timee` is set away from zero so the leading `diagfld` does not run. That is
what makes the saturation reference computable at all: the first `thermo` then
reads exactly the hydrostatic column the test seeded, where after a `diagfld`
it would read one the test cannot reproduce without reimplementing `fromztop`.

## Coverage

Host branch only. The device side is
`tests_cuda.f90::test_thermodynamics_device`, which runs under
`UDALES_RUN_CUDA_SELFTEST` on a Debug GPU build, and the parity cases
`thermo-scalar-wall`, `nudging`, `temperature-kappa` and `stability-wall` in
`tests/integration/gpu/case_matrix.json`.

`thermodynamics` keeps its host branch in a GPU build - `readinitfiles` calls
it before `initCUDA`, and that call is what seeds every device array the loop
path reads before it writes - so this runmode runs against a GPU build too.

## Running

```bash
tests/integration/thermodynamics/run_test.sh
```

Environment overrides: `UDALES_BUILD`, `CASE_SOURCE`, `NAMELIST_SOURCE`,
`MPI_NAMELIST_SOURCE`, `NPROCS`, `MPIEXEC`.

## Passes

| pass | what it adds |
| --- | --- |
| one rank | the five checks, over the whole y extent |
| two ranks in y | that the checks run over the rank's own `jb..je`, and that the slab reduction inside `diagfld` still produces a usable column when it does |
