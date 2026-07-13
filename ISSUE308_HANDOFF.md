# Handoff prompt — Issue #308 "Fix compiler warnings"

> Paste everything below into Claude Code on the machine that has the Fortran/MPI
> toolchain. It documents work already done on branch `fixcompilerwarnings` and
> what remains. Delete this file before opening the final PR.

---

You are continuing work on **uDALES issue #308 — "Fix compiler warnings"**, on the
branch **`fixcompilerwarnings`** (already pushed to `origin`). The goal is a clean
debug build with the temporary warning suppression removed.

## Context: the three warning classes (Intel diagnostic codes from the issue)

When Intel debug builds are run with all diagnostics on, the code emits ~466
warnings and ~1,681 remarks in three classes:

1. **`#6717`** — implicit typing (no `IMPLICIT NONE`).
2. **`#8889`** — missing explicit `EXTERNAL`/interface, mostly around MPI calls.
3. **`#7712`** — unused variables and unused `use ... only:` imports.

The Intel debug build currently hides classes 2 and 3 via
`-diag-disable=7712,8889` in `CMakeLists.txt` (see line ~57). **"Done" = fix the
underlying code and remove that suppression, with a clean debug build on both
gfortran and Intel.**

Key build facts:
- `CMakeLists.txt` GNU debug flags (~line 43): `-Wall -Wextra -Wuninitialized -Warray-bounds -Wconversion` (NOT suppressed — CI already surfaces these).
- `CMakeLists.txt` Intel debug flags (~line 57): `-warn all -diag-disable=7712,8889`.
- CI (`.github/workflows/ci.yml`) builds **gfortran only** (Ubuntu + macOS, Debug + Release) and runs the regression suite. It does **not** exercise Intel, so `#8889`/`#7712` Intel codes can only be verified with a local Intel build.

## What has ALREADY been done on this branch (2 commits, pushed)

**Commit `4e6ef86` — `#6717` (complete):**
- `src/wfmneutral.f90` was the *only* source file lacking `IMPLICIT NONE`.
  Added `IMPLICIT NONE`, declared the previously-implicit real `utangInt`, and
  removed 7 unused integer locals (`kl, ku, im, jm, ip, jp, kp`).
- This class is believed fully resolved (single file).

**Commit `29e8fa5` — `#8889` (partial — the unambiguous cases):**
- Five files called `MPI_*` with no `use mpi` anywhere in scope → no explicit
  interface. Added `use mpi`:
  - `src/modforces.f90` (module level)
  - `src/modchecksim.f90` (module level)
  - `src/modglobal.f90` (inside `initglobal` only — avoid polluting this core module)
  - `src/scalsource.f90` (inside `createscals`; file has bare subroutines, no module)
  - `src/tstep.f90` (inside `tstep_update`; bare subroutines)
- **Why this is safe (verified by reading `src/modmpi.f90`):** `modmpi` does
  `use mpi` with no `private`, so it already re-exports the real MPI entities
  (`MPI_SUM`, `MPI_MAX`, `MPI_INTEGER`, ...). `my_real` is `modmpi`'s own integer,
  not an MPI name. So a scope that has both `use mpi` and
  `use modmpi, only: mpi_max` references the *same* entity — no name clash.
- The other 12 MPI-using files already have **module-level `use mpi`**
  (`heatpump, initfac, modEB, modboundary, modfielddump, modibm, modinlet,
  modmpi, modpurifiers, modsave, modstartup, modstatsdump, modsubgrid,
  modtimedep, readinput, tests, vegetation`), so their contained subroutines get
  MPI interfaces via **host association** — expected to satisfy `#8889` without
  edits. **This assumption is unverified** (needs Intel). See task 2 below.

## What REMAINS

**Task 1 — Confirm the pushed changes build clean on gfortran.**
- Check CI for commit `29e8fa5` passed (Debug + Release, both OSes, regression
  tests green), OR build locally:
  ```
  source .github/scripts/setup_mpi_env.sh <os>   # sets UDALES_MPI_FORTRAN_COMPILER etc.
  cmake -S . -B build/Debug -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_Fortran_COMPILER="$UDALES_MPI_FORTRAN_COMPILER" \
        -DMPI_Fortran_COMPILER="$UDALES_MPI_FORTRAN_COMPILER" \
        -DMPI_C_COMPILER="$UDALES_MPI_C_COMPILER"
  cmake --build build/Debug -- -j 4 2>&1 | tee build_debug.log
  ```

**Task 2 — Verify/finish `#8889` on Intel.**
- Do an Intel debug build with the suppression temporarily lifted:
  edit `CMakeLists.txt` line ~57 to drop `-diag-disable=7712,8889` (or add
  `-diag-enable=8889`), configure with `ifx`/`ifort`, build, capture the log.
- Confirm the 12 host-association files emit **no** `#8889`. If any still do,
  the host-association assumption is wrong → add `use mpi` per-subroutine to
  each scope that calls `MPI_*` in those files (surgical, same pattern as the
  5 files already fixed).

**Task 3 — `#7712` unused variables / imports (the bulk — ~1,681 remarks, ~40 files).**
- Drive this from the **actual compiler output**, not by eye:
  - gfortran: `-Wall -Wextra` prints `Unused variable` / `Unused dummy argument`.
  - Intel: `-warn all` (suppression removed) prints remark `#7712`.
  - Grep the build log: `grep -E "Unused|#7712|7712" build_debug.log`
- Per the issue, the highest-volume files are:
  `src/modstartup.f90, src/modibm.f90, src/modsave.f90, src/modforces.f90,
   src/moddriver.f90, src/modsubgrid.f90, src/modpois.f90, src/modboundary.f90`.
- Fix file-by-file: remove genuinely-unused locals and trim over-broad
  `use ..., only:` lists. **Caution — these are physics kernels.** Before deleting
  a "unused" variable confirm it is not used via `equivalence`, namelist I/O, an
  `#ifdef`/commented block that may be re-enabled, or passed by host association.
  Rebuild after each file and keep the regression suite green.
- Rebuild after each batch and re-grep to watch the count fall.

**Task 4 — Remove the suppression.**
- Delete `-diag-disable=7712,8889` from `CMakeLists.txt` (~line 57) and the
  associated `# FIXME:` comment block above it (~lines 50-54).

**Task 5 — Verify + PR.**
- Clean debug build on **both** gfortran and Intel (zero warnings in the three
  classes). Run the regression suite: `python3 tests/run_tests.py <selection>
  --build-type Debug` (see `.github/workflows/ci.yml` for how CI invokes it, and
  `tests/test_suites.yml` for selections). Results must be unchanged.
- Scope decision the user already made: **full fix of all three classes**, and
  **do NOT** add `-Werror`/`-diag-error` to CI in this PR.
- Open a PR referencing #308. In the description, note that Intel verification
  was done locally (CI is gfortran-only). Delete `ISSUE308_HANDOFF.md`.

## Useful facts already established
- Only `wfmneutral.f90` ever lacked `IMPLICIT NONE`.
- `modmpi` re-exports MPI names (no `private`); `my_real` is its own integer.
- Files that call `MPI_*`: see the two lists above (5 fixed, 12 rely on host
  association, plus `modpois`/`moddriver` call no live MPI).
- Branch is based on `master` as of 2026-07-13 (submodules `2decomp-fft` and
  `tools/View3D` initialized).
