# Fortran Intrinsic Support Probe

Checks that the Fortran intrinsics `src/` depends on are actually implemented by
the compiler's runtime, not merely accepted at compile time.

## Why this exists

NVHPC 24.11 compiles `findloc()` on a logical array without complaint but has no
runtime implementation for it. A GPU build succeeds, and the program aborts only
when the call is first reached:

```
0: FINDLOC: unimplemented for data type
```

In `modibm` that call sat inside the IBM wall-function reconstruction branch,
which is only entered by geometries with non-grid-aligned or near-wall facets.
Every grid-aligned case built and ran fine, so the gap surfaced on a user's
simulation rather than in CI. gfortran implements the intrinsic, so no CPU build
ever showed a problem.

The general shape of the problem is a compiler-specific runtime gap in a rarely
executed code path. This test closes it at build time instead.

## What it does

1. **Capability matrix.** Runs each intrinsic in its own process, so a runtime
   abort in one probe does not hide the others, and records supported or
   unsupported per compiler.
2. **Source attribution.** For each unsupported intrinsic, greps `src/` for code
   that uses it, ignoring comments. If `src/` uses it, the test fails and prints
   the offending lines. If nothing uses it, the gap is reported but not fatal.

Step 2 is what makes the test self-maintaining: it does not encode a list of
banned constructs, it fails when the solver actually depends on something the
target compiler cannot run. Reintroducing a logical `findloc` into `src/` fails
the test on the NVHPC toolchain.

## Running

```bash
tests/unit/fortran_intrinsics/run_test.sh
```

By default every compiler in `gfortran nvfortran` that is on `PATH` is tested.
Set `FC` to test specific compilers:

```bash
FC="gfortran nvfortran" tests/unit/fortran_intrinsics/run_test.sh
```

Needs no uDALES build, and runs in a couple of seconds.

## Expected output

On a machine with both compilers and the current sources, `findloc-logical` is
reported unsupported for nvfortran and not fatal, because `src/` no longer uses
it:

```
-- gfortran (GNU Fortran ...)
     findloc-logical    supported
     ...
-- nvfortran (nvfortran 24.11-0 ...)
     findloc-logical    UNSUPPORTED (not used by src/, not fatal): 0: FINDLOC: unimplemented for data type
     findloc-real       supported
     findloc-integer    supported
     count-logical      supported
```

## Adding a check

1. Add a `probe_*` subroutine to `intrinsic_probe.f90` and a `case` for it.
2. Add its name to `CHECKS` in `run_test.sh`.
3. Add a pattern to `pattern_for()` matching how it appears in Fortran source.
   Leave the pattern empty to report the capability without ever failing.

## See also

`tests/integration/ibm_cell_lookup` validates the replacement that `modibm` now
uses in place of `findloc`.
