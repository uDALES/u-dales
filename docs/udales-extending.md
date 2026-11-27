# Extending uDALES

This page collects practical recipes for common developer tasks: adding a namelist option, adding an output field, adding a statistic, adding an example case, and adding a test. Each recipe is derived by tracing a real, working example through the code, with file and function references, so you can use that example as a template rather than starting from scratch. See the [Fortran API reference](udales-software-docs.md) for generated source documentation and the [namoptions overview](udales-namoptions-overview.md) for the full input parameter list.

## Adding a namelist option

Template: the `lydump` switch (`&OUTPUT`, y-averaged statistics output) is a simple, working logical option you can copy directly.

1. **Declare the variable and its default** in `src/modglobal.f90`, e.g. `logical :: lydump = .false.` (declared alongside the other output switches around line 200). Pick a sensible default that preserves current behaviour when the option is absent from `namoptions`.
2. **Add it to the relevant namelist group** in `src/in_readnamelists.f90`, subroutine `read_namelist_inputs`. There are two places to touch:
   - the `use modglobal, only : ...` import list in module `readnamelists` (`lydump` is imported there alongside the other `OUTPUT`-group switches)
   - the `namelist/OUTPUT/ &` declaration itself (`lydump` is listed there together with `ltdump, lytdump, lxydump, lxytdump, lmintdump, ...`)
   If you are adding a wholly new namelist group instead of extending an existing one, also add a `read (ifnamopt, YOURGROUP, iostat=ierr)` block (with the same `iostat` error-handling pattern used for `OUTPUT`, `BC`, etc.) and a `rewind (ifnamopt)` after it.
3. **Broadcast it from rank 0 to all MPI ranks** in `broadcast_config_parameters`, also in `src/in_readnamelists.f90`: `call MPI_BCAST(lydump, 1, MPI_LOGICAL, 0, comm3d, mpierr)`. Match the MPI type to the Fortran type (`MPI_LOGICAL`, `MPI_INTEGER`, `MY_REAL`, or `MPI_CHARACTER` with `len(...)` for strings).
4. **Use the value** wherever it is needed. `lydump` is consumed in two ways that are both worth knowing about as patterns:
   - a sanity check right after the `read (ifnamopt, OUTPUT, iostat=ierr)` in `read_namelist_inputs` (see the `receptor_height` checks there) or in `checkinitvalues` (`src/modstartup.f90`) for checks that need other groups, each ending in `stop 1` with an explicit message.
   - the actual behaviour, in `src/out_stats.f90`, which imports it with `use modglobal, only : ..., lydump, lytdump, ...` at module scope and gates the allocation in `stats_prepare`, the NetCDF setup (`stats_createnc_yavg`) and the sampling/writing in `stats_main` behind `if (lydump) then ... end if`.
5. **Document the new option** in `docs/udales-namoptions-overview.md`, adding a row to the table for the namelist group it belongs to (name, default, possible values, description, unit), and in `docs/schemas/udales_input_schema.json` (type, default, description). `tests/integration/input/udales_variable_indexer.py` indexes namelist, broadcast and schema entries so that the three can be compared.

Checklist:

- [ ] Declare in `src/modglobal.f90` with a safe default
- [ ] Add to the `use modglobal, only : ...` list and the `namelist/GROUP/` statement in `read_namelist_inputs` (`src/in_readnamelists.f90`)
- [ ] Broadcast with `MPI_BCAST` in `broadcast_config_parameters`
- [ ] Consume the value at its point of use, importing it with `use modglobal, only : ...` in that module
- [ ] Add a row to `docs/udales-namoptions-overview.md` and an entry to `docs/schemas/udales_input_schema.json`

**Do not rename existing namelist keys or groups.** Live, user-facing namelist names (e.g. `idriver`, `&DRIVER`) must stay stable across releases — case files in the wild reference them by name. Renaming an existing option is a separate, explicitly-approved task, not something to bundle into feature work.

## Adding a variable to the instantaneous output

Template: `src/out_instants.f90` (module `instant`) handles the 3-D fields generically — adding a new instantaneous 3-D output field only means adding one `case` branch in `instant_field_init`, using an existing one (e.g. `case('mc')` for the IBM cell mask, or `case('di')` for the divergence, which is computed at output time) as a model.

How it works: the namelist option `fieldvars` (declared in `src/modglobal.f90`, read in the `&OUTPUT` group in `src/in_readnamelists.f90`) is a comma-separated string of two-character labels, e.g. `'u0,v0,w0,th,s1'`. In `instant_field_init`, `nfieldvars` counts the labels (each label plus its separator is 3 characters), then a loop over `n = 1, nfieldvars` does `select case(fieldvars(3*n-2:3*n-1))` to pick out each label in turn. A `case` branch does two things: it calls `ncinfo(fldVars(n,:), 'ncvarname', 'long description', 'units', 'gridtag')` to describe the NetCDF variable, and it points a module-level pointer at the data with `pfields(n)%point => somefield(ib:ie,jb:je,kb:ke)`. For a quantity that has to be computed at output time instead of read from an existing array, the branch calls `register_derived(n, DERIV_XXX)` and `compute_derived` fills a buffer for that `DERIV_XXX` code (add a code and a `case` there). The generic `instant_field_main` then writes every variable with `writeoffset` — it has no per-variable logic beyond the derived fields `tha`, `rh` and `pabs`, so you normally never need to touch it. `instant_main` runs the field, slice and probe writers at the same instants, `tinstantstart + n*tinstantdump`, on substep 3 of the Runge–Kutta step.

There are two nearly identical `select case` blocks in `instant_field_init`: one for `lhalos_out` (ghost cells included) and one for the normal case. `lhalos_out` is hard-coded `.false.`, so in practice only the second block (the one that also knows `tx`, `ty`, `tz`, `hf`, `mu`, `mv`, `mw`, `mc`, `pd`–`pw`, `rs`–`dw`, `di`–`wz`) is exercised; new fields only need adding there. Update the error message of its `case default` as well, since it lists the accepted labels.

The grid stagger tag (last `ncinfo` argument, e.g. `'tttt'`, `'mttt'`, `'ttmt'`) records which grid (cell centre `t` or staggered `m`) each of the x, y, z, time dimensions uses — copy the tag from a field with the same staggering as yours (a cell-centred scalar like pressure uses `'tttt'`; a u-staggered quantity like `tau_x` uses `'mttt'`).

Slices (`slicevars`) and probes (`probevars`) have their own, shorter label lists: `instant_ncdescription_islice`/`_jslice`/`_kslice` for the slices and `instant_gather_probe_var` plus `instant_add_var_atts` for the probes. A new label must be added to each of those if it should be available there too.

Checklist:

- [ ] Confirm the field you want to output already exists as a module-level array (e.g. in `modfields.f90` or wherever it is computed) and is reachable via a `use` statement in `out_instants.f90`; otherwise add a `DERIV_*` code and compute it in `compute_derived`
- [ ] Add a `case('xx')` branch (pick an unused two-character label) in the non-halo `select case(fieldvars(3*n-2:3*n-1))` block of `instant_field_init`, with `ncinfo(...)` and `pfields(n)%point => yourfield(ib:ie,jb:je,kb:ke)` (or `register_derived`)
- [ ] Extend the label list in the `case default` error message
- [ ] Add the label to the `fieldvars` row in `docs/udales-namoptions-overview.md` and to the `fieldvars` description in `docs/schemas/udales_input_schema.json`, and the variable to the `ins_field` row of `docs/udales-output-files.md`
- [ ] Add it to the slice and probe label lists only if it is wanted there

## Adding a statistic to stats_xyt

Template: the air temperature `tha` in `src/out_stats.f90` (module `stats`). It is the simplest time-averaged statistic because it is a plain time average of a 3-D field followed by a spatial average, with no correlation or variance term to subtract — a good starting point before tackling flux or variance statistics.

`stats_xyt.<expnr>.nc` holds statistics averaged over x, y *and* time (as opposed to `stats_xy`, which is only x/y-averaged, or `stats_t`, which is only time-averaged). All time-averaged files share the same 3-D accumulators; the module is organised by physical group (`vel`, `temp`, `moist`, `relative_humidity`, `scalar`, `tree`) and by stage (`allocate`, `ncdescription`/`init`, `reset`, `compute`, `write`), and `stats_main` calls the stages in order every `tsample` and at the end of every `tstatsdump` window. For `tha` the pipeline is:

1. **3-D accumulator.** `stats_allocate_tavg_temp` allocates and zeroes `tempt(ib:ie,jb:je,kb:ke+kh)`; `stats_reset_tavg_temp` zeroes it again at the start of every averaging window. (The accumulators are allocated in `stats_prepare`, at the first sample, not at startup.)
2. **Time average.** `stats_compute_tavg_temp` calls `stats_compute_tavg(tempt, tempnow)` every `tsample`, the running mean `(old*(tstatsdumpp-tsamplep) + new*tsamplep)*tstatsdumppi` that every accumulator uses. `tempnow` is the instantaneous field, computed in `stats_calculate_air_temperature` just before.
3. **Output in `stats_t`.** `stats_ncdescription_tavg_temp` registers the variable with `ncinfo(tVars(ctrt+7,:), 'tha', 'Air temperature', 'K', 'tttt')` and bumps the counter `ctrt`; `stats_tavg_field` maps the name back to the accumulator (`case ('tha'); fld = tempt(:,:,kb:ke)`), and `stats_write_tavg` loops over all registered names, so no per-variable write code is needed. The same table serves the plane files `stats_tislice`, `stats_tjslice` and `stats_tkslice`.
4. **Output in `stats_xyt`.** `stats_init_xytavg_temp` allocates the profile `tempxyt(kb:ke+kh)` and registers `ncinfo(xytVars(ctrxyt+6,:), 'tha', 'Air temperature', 'K', 'tt')`; `stats_compute_xytavg_temp` reduces the 3-D accumulator with `spatial_avg(tempxyt, tempt(...), kb, ke, IIc(...), IIcs(...), .false.)` (use `IIw`/`IIws` and `.true.` for w-staggered quantities so that masked levels become `-999.`); `stats_write_xytavg_temp` writes it with `writestat_nc(ncidxyt, 'tha', tempxyt(kb:ke), nrecxyt, zdim)`.

The index in the `ncinfo` call (`ctrxyt+6`) only has to be unique within the group and the group counter must be bumped by the number of variables the group registers; the `writestat_nc` call looks the variable up by name.

Checklist (for `stats_xyt`; `stats_xy`, `stats_yt` and `stats_y` follow the same pattern with `xy`, `yt`, `y` in the subroutine names):

- [ ] If the quantity is a time average, add an accumulator in the matching `stats_allocate_tavg_<group>`, zero it in `stats_reset_tavg_<group>`, and accumulate it in `stats_compute_tavg_<group>` with `stats_compute_tavg`
- [ ] Register it for `stats_t` in `stats_ncdescription_tavg_<group>` (bump the counter) and add its `case` to `stats_tavg_field`
- [ ] Allocate the profile and register the variable in `stats_init_xytavg_<group>` (bump `ctrxyt`), reduce it in `stats_compute_xytavg_<group>` with `spatial_avg`, and write it in `stats_write_xytavg_<group>`
- [ ] Guard a new physical group behind its switch (`ltempeq`, `lmoist`, `nsv > 0`, ...) in `stats_prepare` and `stats_main`, as the existing groups are
- [ ] Use the same variable name in every file; the file name carries the averaging
- [ ] If `gather_outputs.sh` must concatenate the variable along x (u-staggered) or y (v-staggered), add it to the `proposed_xmparam`/`proposed_ymparam` lists there
- [ ] Document the variable in the relevant rows of `docs/udales-output-files.md`

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

`tools/local_execute.sh` reads `config.sh` from the case directory, checks `DA_WORKDIR`/`DA_BUILD`/`DA_TOOLSDIR`/`NCPU` are set, and runs the case. `tools/examples/run_examples.sh` is the local sweep that runs the examples back to back: by default every case in `examples/` in ascending order (so the driver case `949` runs before the driven case `950`, whose driver files are linked in with `tools/link_driver_files.sh`), except those listed in its `EXCLUDED` variable (currently `024`, a 1024-core HPC case). It ignores each case's `config.sh`, taking `DA_BUILD`/`DA_WORKDIR`/`DA_TOOLSDIR` from the environment (defaulting to `build/release/u-dales`, `outputs/` and `tools/`) and `NCPU` from `nprocx*nprocy`. uDALES runs one MPI rank per CPU and never oversubscribes, so the full sweep needs 8 CPUs (cases `001` and `999` use 8 ranks, the others 4); on a 4-CPU machine pass the 4-rank cases explicitly, e.g. `./tools/examples/run_examples.sh 002 101 102 201 949 950`. A new case is therefore picked up automatically; add it to `EXCLUDED` (with a comment) if it cannot run in the routine sweep, or add any extra setup it needs to the loop, as done for `102` (warmstart files) and `950` (driver files).

Checklist:

- [ ] Create `examples/<NNN>/` with `namoptions.<NNN>` and `config.sh` (edit paths, don't commit real absolute paths pointing outside the repo if avoidable)
- [ ] Generate geometry/IBM inputs via pre-processing (or copy+adapt from a similar example) and make sure `&WALLS` counts match
- [ ] Add `prof.inp.<NNN>`, `lscale.inp.<NNN>`, and any case-specific `.inp` files
- [ ] Add `info.txt` describing the case
- [ ] Test with `./u-dales/tools/local_execute.sh examples/<NNN>`
- [ ] Check the case runs with `tools/examples/run_examples.sh <NNN>` (it is included in the routine sweep automatically; add it to `EXCLUDED` there if it should not be)
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

`tests/integration/ibm_sparse_input/run_test.sh` is a good template for a solver-driven suite: it resolves `UDALES_BUILD`/`CASE_SOURCE`/`NAMELIST_SOURCE` from environment variables with sane defaults, copies a fixture case (`tests/cases/101`) plus a test-local namelist into a scratch `mktemp -d` directory, derives `NPROCS` from `nprocx`/`nprocy` in the namelist, sources `tools/mpi_env.sh` for the launcher and the one-rank-per-CPU pre-flight check (`mpi_preflight "$NPROCS"`; uDALES never oversubscribes), runs `"$MPIEXEC" $MPI_LAUNCH_EXTRA_ARGS -n "$NPROCS" "$UDALES_BUILD" "$NAMELIST"`, and propagates the exit code — the same pattern `tests/run_tests.py` expects from any suite command.

Checklist:

- [ ] Decide the scope (`unit`/`integration`/`system`/`regression`) and put the test file(s) under the matching `tests/<layer>/` subdirectory (or `tools/python/tests/` for a pure-Python unit test)
- [ ] If it needs a case fixture, add it under `tests/cases/` (reuse an existing one, e.g. `tests/cases/100` or `101`, if it fits) rather than duplicating a full case
- [ ] Write the test/script so it exits non-zero on failure (Python: raise/assert; shell: propagate `$?`, as in `run_test.sh`)
- [ ] Add a suite entry to `tests/test_suites.yml` under the appropriate group (usually `supported` if it should merge-gate, `experimental` otherwise), filling in `label`, `class`, `kind`, `component`, `platform`, `cost`, `command`, and any `env_<NAME>` needed
- [ ] Run it via the dispatcher to confirm wiring: `python tests/run_tests.py <group>` (e.g. `python tests/run_tests.py supported`)
- [ ] Update `tests/README.md` if you introduced a new fixture under `tests/cases/` or a new suite category worth documenting there
