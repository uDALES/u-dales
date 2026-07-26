# Extending uDALES

This page collects practical recipes for common developer tasks: adding a namelist option, adding an output field, adding a statistic, adding an example case, and adding a test. Each recipe is derived by tracing a real, working example through the code, with file and function references, so you can use that example as a template rather than starting from scratch. See the [Fortran API reference](udales-software-docs.md) for generated source documentation and the [namoptions overview](udales-namoptions-overview.md) for the full input parameter list.

## Adding a namelist option

Template: the `lydump` switch (`&OUTPUT`, y-averaged statistics output) is a simple, working logical option you can copy directly.

1. **Declare the variable and its default** in `src/modglobal.f90`, e.g. `logical :: lydump = .false.` (declared alongside the other output switches around line 200). Pick a sensible default that preserves current behaviour when the option is absent from `namoptions`.
2. **Add it to the relevant namelist group** in `src/modstartup.f90`, subroutine `readnamelists`. There are two places to touch:
   - the `use modglobal, only : ...` import list at the top of `readnamelists` (`lydump` is imported there alongside the other `OUTPUT`-group switches)
   - the `namelist/OUTPUT/ &` declaration itself (`lydump` is listed there together with `ltdump, lytdump, lxydump, lxytdump, lmintdump, ...`)
   If you are adding a wholly new namelist group instead of extending an existing one, also add a `read (ifnamopt, YOURGROUP, iostat=ierr)` block (with the same `iostat` error-handling pattern used for `OUTPUT`, `BC`, etc.) and a `rewind (ifnamopt)` after it.
3. **Broadcast it from rank 0 to all MPI ranks**, still in `readnamelists`: `call MPI_BCAST(lydump, 1, MPI_LOGICAL, 0, comm3d, mpierr)`. Match the MPI type to the Fortran type (`MPI_LOGICAL`, `MPI_INTEGER`, `MY_REAL`, or `MPI_CHARACTER` with `len(...)` for strings).
4. **Use the value** wherever it is needed. `lydump` is consumed in two ways that are both worth knowing about as patterns:
   - a sanity check in `checknamelistvalues` (also in `modstartup.f90`): the value is pulled in via `use modglobal, only : ..., lydump, lytdump, ...` and checked, e.g. `if ((lydump .or. lytdump) .and. (nprocx > 1)) then ... stop 1`.
   - the actual behaviour, in `src/modstatsdump.f90`, which imports it with `use modglobal, only : dt, lydump, lytdump, ...` at module scope and gates both the NetCDF setup in `initstatsdump` and the sampling/writing in `statsdump` behind `if (lydump) then ... end if`.
5. **Document the new option** in `docs/udales-namoptions-overview.md`, adding a row to the table for the namelist group it belongs to (name, default, possible values, description, unit).

Checklist:

- [ ] Declare in `src/modglobal.f90` with a safe default
- [ ] Add to the `use modglobal, only : ...` list and the `namelist/GROUP/` statement in `readnamelists` (`src/modstartup.f90`)
- [ ] Broadcast with `MPI_BCAST` in the same subroutine
- [ ] Consume the value at its point of use, importing it with `use modglobal, only : ...` in that module
- [ ] Add a row to `docs/udales-namoptions-overview.md`

**Do not rename existing namelist keys or groups.** Live, user-facing namelist names (e.g. `idriver`, `&DRIVER`) must stay stable across releases — case files in the wild reference them by name. Renaming an existing option is a separate, explicitly-approved task, not something to bundle into feature work.

## Adding a variable to fielddump

Template: `src/modfielddump.f90` already handles this generically — adding a new instantaneous 3D output field only means adding one `case` branch, using an existing one (e.g. `case('mc')` for the IBM cell mask, or `case('di')` for divergence, which is computed in-place) as a model.

How it works: the namelist option `fieldvars` (declared in `src/modglobal.f90`, read in the `&OUTPUT` group in `modstartup.f90`) is a comma-separated string of two-character field codes, e.g. `'u0,v0,w0,th,s1'`. In `initfielddump`, `nvar = (LEN(trim(fieldvars))+1)/3` counts the fields (each code plus its separator is 3 characters), then a loop over `n=1,nvar` does `select case(fieldvars(3*n-2:3*n-1))` to pick out each two-character code in turn. Each `case` branch does two things: calls `ncinfo(ncname(n,:), 'ncvarname', 'long description', 'units', 'gridtag')` to describe the NetCDF variable, and points a module-level pointer at the data with `pfields(n)%point => somefield(ib:ie,jb:je,kb:ke)`. The generic `fielddump` subroutine later copies every `pfields(n)%point` into the output buffer and writes it via `writestat_nc` — it has no per-variable logic, so you never need to touch it.

There are two nearly-identical `select case` blocks in `initfielddump`: one for `lhalos` (ghost cells included) and one for the normal case. `lhalos` is hard-coded `.false.` (see the assignment right above the `MPI_BCAST` calls), so in practice only the non-halo branch (starting `else` around the `do n=1,nvar` loop that includes `case('u0')`, ..., `case('mc')`, `case('di')`) is exercised — new fields normally only need adding there, and to the `lhalos` branch too if you want the option kept consistent (some existing fields, like `tx`/`ty`/`tz`/`hf`/`mu`/`mv`/`mw`/`mc`/`di`, are only implemented in the non-halo branch).

The grid stagger tag (last `ncinfo` argument, e.g. `'tttt'`, `'mttt'`, `'ttmt'`) records which grid (cell-centre `t` vs. staggered `m`) each of the x, y, z, time dimensions uses — copy the tag from a field with the same staggering as yours (a cell-centred scalar like pressure uses `'tttt'`; a u-velocity-staggered flux like `tau_x` uses `'mttt'`).

Checklist:

- [ ] Confirm the field you want to output already exists as a module-level array (e.g. in `modfields.f90` or wherever it's computed) and is reachable via a `use` statement in `modfielddump.f90`
- [ ] Add a `use modsomething, only : yourfield` import in `initfielddump` if it isn't already imported
- [ ] Add a `case('xx')` branch (pick an unused two-character code) in the non-halo `select case(fieldvars(3*n-2:3*n-1))` block, with `ncinfo(...)` and `pfields(n)%point => yourfield(ib:ie,jb:je,kb:ke)`
- [ ] Mirror the branch in the `lhalos` block only if halo output for this field is wanted
- [ ] Add the new code to the `fieldvars` row's list of labels in `docs/udales-namoptions-overview.md`
- [ ] No changes needed in `fielddump` or `exitfielddump` — they iterate over `pfields` generically

## Adding a statistic to xytdump

Template: `pxyt` (domain-mean pressure) in `src/modstatsdump.f90`. It is the simplest xy/time-averaged statistic because it is a plain spatial-then-time average with no correlation/variance term to subtract — a good starting point before tackling flux or variance statistics.

`xytdump.xxx.nc` holds statistics that are averaged over x, y, *and* time (as opposed to `xydump`, which is only x/y-averaged, or `tdump`, which is only time-averaged). The pipeline for `pxyt` runs through three stages in `statsdump`:

1. **Instantaneous xy-average** (`if (lxydump .or. lxytdump) then` block, module-level variable `pxy`, declared `real, dimension(kb:ke+kh) :: pxy` inside `modstatsdump.f90`): `call avexy_ibm(pxy(kb:ke+kh), pres0(ib:ie,jb:je,kb:ke+kh), ib,ie,jb,je,kb,ke,kh, IIc(...), IIcs(kb:ke+kh), .false.)` averages the instantaneous `pres0` field over x and y at each sample.
2. **Running time-average** (`if (lxytdump) then` block further down): `pxyt(kb:ke+kh) = (pxyt(kb:ke+kh)*(tstatsdumpp-tsamplep) + pxy(kb:ke+kh)*tsamplep)*tstatsdumppi` accumulates `pxy` into the persistent, module-level array `pxyt`. `pxyt` is declared in `src/modfields.f90` (`real, allocatable :: pxyt(:)`) and allocated/zeroed alongside the other `xyt` accumulators (`allocate(pxyt(kb:ke+kh))`, `pxyt=0.`) in the same allocation block as `uxyt`, `vxyt`, etc.
3. **NetCDF definition and write**, both in `initstatsdump`/`statsdump`:
   - in `initstatsdump`, inside `if (lxytdump) then`, `call ncinfo(ncstatxyt(6,:), 'pxyt', 'Pressure', 'm^2/s^2', 'tt')` registers the variable at index 6 of the `nstatxyt`-sized `ncstatxyt` array (`ncstatxyt` is declared in `modfields.f90` and allocated `allocate(ncstatxyt(nstatxyt,4))` in `initstatsdump`);
   - in `statsdump`, inside the final `if (lxytdump) then ... if (myid == 0) then` block, `varsxyt(:,6) = pxyt(kb:ke)` copies it into the write buffer, which is then written in one call: `call writestat_1D_nc(ncidxyt, nstatxyt, ncstatxyt, varsxyt, nrecxyt, khigh-klow+1)`.

Note the index (`6`) must be consistent between the `ncinfo` call and the `varsxyt(:,6) = ...` assignment, and `nstatxyt` (currently `23`, set in the module header alongside `nstaty`, `nstatt`, etc.) must be bumped if you add a new slot rather than reusing one.

Checklist (for `xytdump`; the same pattern applies to `tdump`/`ncstatt`/`varst`, `xydump`/`ncstatxy`/`varsxy`, etc., with `t`/`xy` suffixes swapped accordingly):

- [ ] If accumulating in time, add a persistent array (e.g. `real, allocatable :: myvarxyt(:)`) in `src/modfields.f90`, and allocate + zero it in the same block as `uxyt`, `vxyt`, ...
- [ ] Compute the instantaneous spatial average with `avexy_ibm(...)` (or reuse an existing `xy` intermediate) inside the `if (lxydump .or. lxytdump)` block in `statsdump`
- [ ] Accumulate it into the time-average with the same `(old*(tstatsdumpp-tsamplep) + new*tsamplep)*tstatsdumppi` running-mean pattern, inside `if (lxytdump)`
- [ ] Bump `nstatxyt` in the module header and `call ncinfo(ncstatxyt(N,:), 'name', 'description', 'units', 'gridtag')` in `initstatsdump`
- [ ] Add `varsxyt(:,N) = myvarxyt(kb:ke)` next to the write call in `statsdump`
- [ ] Document the new output variable in `docs/udales-output-files.md` if that page enumerates `xytdump` fields

## Adding an example case

Template: compare `examples/001` (minimal, no buildings, pressure-gradient forcing) with `examples/102` (buildings, warmstart, scalar point source) — together they show the full range of what an example directory can contain.

Every example under `examples/<NNN>/` needs, at minimum:

- `namoptions.<NNN>` — the namelist file; the last three digits of every filename in the case must match `iexpnr` in `&RUN`.
- `config.sh` — sets `DA_EXPDIR`, `DA_TOOLSDIR`, `DA_BUILD`, `DA_WORKDIR`, `NCPU` (as absolute, user-specific paths — copy and edit one from an existing example, e.g. `examples/001/config.sh` or `examples/102/config.sh`, don't reuse it as-is).
- Geometry/IBM inputs produced by pre-processing: `facets.inp.<NNN>`, `factypes.inp.<NNN>`, `facetarea.inp.<NNN>`, `facet_sections_{c,u,v,w}.txt`, `fluid_boundary_{c,u,v,w}.txt`, `solid_{c,u,v,w}.txt`, and either an `.stl` file (e.g. `examples/001/flat_ground.stl`, `examples/102/geom.102.STL`) or a MATLAB/Python-generated set from `stl_file` in `&INPS`. The corresponding facet/point counts (`nfcts`, `nsolpts_*`, `nbndpts_*`, `nfctsecs_*`) go in `&WALLS` in `namoptions.<NNN>` and are produced by the pre-processing tools, not hand-written.
- `prof.inp.<NNN>` and `lscale.inp.<NNN>` — initial/forcing profiles, also pre-processing outputs.
- `info.txt` — a short human-readable description of the case (see `examples/001/info.txt`, `examples/102/info.txt`).
- Optional, case-specific inputs: `scalar.inp.<NNN>` and `scalarsourcep.inp.1.<NNN>` for scalar sources (as in `102`), `warmstart_files/` and restart files (`initd...`, matched by `startfile` in `&RUN`) for a warmstart (as in `102`), driver/`*driver*` files for driver-driven cases (as in `949`/`950`).

Getting the geometry/IBM inputs, `prof.inp`, and `lscale.inp` right normally means running the pre-processing pipeline (see `docs/udales-pre-processing.md`) rather than writing them by hand.

To exercise a new example locally:

```sh
./u-dales/tools/local_execute.sh examples/<NNN>
```

`tools/local_execute.sh` reads `config.sh` from the case directory, checks `DA_WORKDIR`/`DA_BUILD`/`DA_TOOLSDIR`/`NCPU` are set, and runs the case. `tools/examples/run_examples.sh` is the CI/local sweep that runs several examples back to back (currently listing `001 002 101 102 201 501 502`, cross-check against the current contents of `examples/` before relying on this list — it also handles downloading the extra warmstart/driver assets for `102`/`502` via `curl`+`unzip`). Add your new case number to that loop if it should be part of the routine sweep.

Checklist:

- [ ] Create `examples/<NNN>/` with `namoptions.<NNN>` and `config.sh` (edit paths, don't commit real absolute paths pointing outside the repo if avoidable)
- [ ] Generate geometry/IBM inputs via pre-processing (or copy+adapt from a similar example) and make sure `&WALLS` counts match
- [ ] Add `prof.inp.<NNN>`, `lscale.inp.<NNN>`, and any case-specific `.inp` files
- [ ] Add `info.txt` describing the case
- [ ] Test with `./u-dales/tools/local_execute.sh examples/<NNN>`
- [ ] Add the case to `tools/examples/run_examples.sh` if it should run in the routine sweep
- [ ] Document the case in `docs/udales-example-simulations.md` (setup table row plus a walkthrough section, following the `001`/`002`/`101`/`102` pattern)

## Adding a test

Template: `tests/integration/ibm_sparse_input/run_test.sh` (an MPI/solver-driven integration test) for a case-based suite, or any file in `tools/python/tests/` for a pure-Python unit test.

Tests live in layered locations by scope, documented in `tests/README.md`:

- `tools/python/tests/` — unit tests for Python modules (`udbase`, `udprep`, `udgeom`, ...); picked up automatically by `unittest discover`, no manifest entry needed.
- `tests/unit/` — isolated Fortran/solver routine tests with minimal setup.
- `tests/integration/` — end-to-end, multi-component checks, often against committed fixtures in `tests/cases/` (e.g. `tests/cases/101` is shared by `integration/ibm_sparse_input/`).
- `tests/system/` — heavier whole-code, solver-driven validation.
- `tests/regression/` — branch-to-branch or reference-output comparisons (e.g. `tests/regression/david_tests/`).

New automated test *runs* are wired up through `tests/test_suites.yml`, a hand-parsed YAML manifest read by `tests/run_tests.py`. Its structure (see the schema comment at the top of the file):

- `groups.<name>.suites` — concrete suite entries belonging to that selection.
- `groups.<name>.includes` — other group names folded in (e.g. `supported` includes `python-library`).
- Each suite has `label`, `class` (`supported` = required merge-gate, `experimental` = not yet gating), `kind` (`unit`/`integration`/`reference`/`system`/`regression`), `component`, `platform` (`linux`/`macos`/`hpc`/`any`), `cost` (`fast`/`medium`/`slow`), and `command` — an argv list passed straight to `subprocess.run`, with `{python}`, `{repo_root}`, `{tests_dir}`, `{branch_a}`, `{branch_b}`, `{build_type}`, `{build_type_lower}` available as format placeholders. Optional `env_<NAME>` keys inject environment variables into the child process (e.g. `env_UDALES_BUILD: "{repo_root}/build/{build_type_lower}/u-dales"` used by the `ibm_sparse_input`, `mpi_operators`, and `processor_boundaries` suites).

`tests/integration/ibm_sparse_input/run_test.sh` is a good template for a solver-driven suite: it resolves `UDALES_BUILD`/`CASE_SOURCE`/`NAMELIST_SOURCE` from environment variables with sane defaults, copies a fixture case (`tests/cases/101`) plus a test-local namelist into a scratch `mktemp -d` directory, derives `NPROCS` from `nprocx`/`nprocy` in the namelist, runs `mpiexec ... "$UDALES_BUILD" "$NAMELIST"`, and propagates the exit code — the same pattern `tests/run_tests.py` expects from any suite command.

Checklist:

- [ ] Decide the scope (`unit`/`integration`/`system`/`regression`) and put the test file(s) under the matching `tests/<layer>/` subdirectory (or `tools/python/tests/` for a pure-Python unit test)
- [ ] If it needs a case fixture, add it under `tests/cases/` (reuse an existing one, e.g. `tests/cases/100` or `101`, if it fits) rather than duplicating a full case
- [ ] Write the test/script so it exits non-zero on failure (Python: raise/assert; shell: propagate `$?`, as in `run_test.sh`)
- [ ] Add a suite entry to `tests/test_suites.yml` under the appropriate group (usually `supported` if it should merge-gate, `experimental` otherwise), filling in `label`, `class`, `kind`, `component`, `platform`, `cost`, `command`, and any `env_<NAME>` needed
- [ ] Run it via the dispatcher to confirm wiring: `python tests/run_tests.py <group>` (e.g. `python tests/run_tests.py supported`)
- [ ] Update `tests/README.md` if you introduced a new fixture under `tests/cases/` or a new suite category worth documenting there
