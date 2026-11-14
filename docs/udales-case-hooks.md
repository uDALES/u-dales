# Adding case-specific physics with hooks

Physics that belongs to one test case, such as the free-troposphere radiation term of DYCOMS-II RF01, does not go into the core modules as a namelist switch. It goes into its own module in `src/cases/` and is switched on at run time by name:

```fortran
&NAMCASES
case_hooks = 'dycoms_rf01'
/
```

The core modules call the registered hooks at fixed points. Without `&NAMCASES`, or with `case_hooks = ''`, no hook is registered and every hook point does nothing, so the results are identical to a build without the hook system.

## How it works

`src/modcases.f90` is the registry:

- `initcases` reads `&NAMCASES` and splits `case_hooks` at commas. For each name, `register_case` calls the matching `register_<name>` routine (a `select case`). An unknown name aborts and lists the valid names; a name listed twice aborts. After registering, `warn_unlisted_groups` tries to read the namelist group of every known case that is not listed and warns on rank 0 if it is present, since the group would be ignored (each case module provides a `<name>_group_present` function for this).
- A `register_<name>` routine fills one `hook_set`, a derived type with one procedure pointer per hook point. Fortran has no arrays of procedure pointers, so the registry keeps an array of these derived types instead.
- The `register_<name>` routine also reads and checks the case's namelist (through a `<name>_config` routine of the case module), so a misconfigured case aborts at startup, before the facet and initial-field setup.
- The core code calls `call_hooks_<point>`, which runs the hook of every registered case that set that point, in `case_hooks` order.
- Errors are reported once, from rank 0, and every rank then stops through `modmpi:abort_all_ranks`.

| Hook point | Called from | When | Typical use |
| ---------- | ----------- | ---- | ----------- |
| `init` | `program.f90`, after `readinitfiles` (cold and warm starts) | once | checks that need the initial fields |
| `radiation` | `modradiation:radiation`, after the generic `radpar` tendency is computed and before `thlprad` is added to `thlp` (only with `iradiation /= 0`) | at each radiation update (`rk3step == 1`, every `timerad`); the contribution is held in `thlprad` until the next update | add to `thlprad` |
| `forcing` | the time loop, after `lstend` | every Runge-Kutta substep | add to the tendencies (`up`, `thlp`, ...) |
| `stats` | the time loop, after `thermostat` | every Runge-Kutta substep: the hook must gate itself on `rk3step == 3` and its own output time, as `thermostat` does | write case-specific diagnostics |

Hooks take no arguments: they use the model state through `use` statements, like any other module.

## Worked example: DYCOMS-II RF01

`src/cases/case_dycoms_rf01.f90` adds the free-troposphere term of the Stevens et al. (2005) radiation to the generic longwave scheme.

1. **Write the case module** in `src/cases/` (CMake compiles every `src/**/*.f90`; re-run the configure step after adding a file). Keep it private by default, with an explicit public list of the hook procedures:

    ```fortran
    module case_dycoms_rf01
      implicit none
      private
      public :: dycoms_rf01_config, dycoms_rf01_init, dycoms_rf01_radiation
      real :: dycoms_D        = 0.   !< required: no nonzero defaults
      real :: dycoms_qt_thres = 0.
    contains
      subroutine dycoms_rf01_config      ! reads &NAMCASE_DYCOMS, aborts if a parameter is unset
      subroutine dycoms_rf01_init        ! warns if no column of the initial qt has an inversion
      subroutine dycoms_rf01_radiation   ! thlprad(i,j,k) -= IIc * (G(zh(k+1)) - G(zh(k))) / dzf(k)
    end module
    ```

    Here `G = F / (rho_i cp) = dycoms_D [(z-zi)^(4/3)/4 + zi (z-zi)^(1/3)]`: the flux is normalised by `rho_i cp` as in DALES, not by the local `rhof cp exnf`.

    The case has its own namelist, `&NAMCASE_DYCOMS`, read from the same `namoptions` file at registration (rank 0 reads, then `MPI_BCAST`), with the usual error handling. Parameters that are not obviously nonzero default to 0, and the run aborts with a clear message when they are not set; the case's `namoptions` sets them explicitly. **Case namelist variables carry the case prefix** (`dycoms_D`, not `D`): `tools/python/namelists.json` maps variable names to namelist groups globally, so a generic name would capture the same name in every other group.

2. **Register it** in `src/modcases.f90`: add the name to `valid_names`, a branch to the `select case` in `register_case`, a check of its group in `warn_unlisted_groups`, and a `register_<name>` routine that sets the pointers of the hook points the case uses:

    ```fortran
    case ('dycoms_rf01')
      call register_dycoms_rf01(hooks(nhooks))
    ...
    subroutine register_dycoms_rf01(h)
      use case_dycoms_rf01, only : dycoms_rf01_config, dycoms_rf01_init, dycoms_rf01_radiation
      type(hook_set), intent(inout) :: h
      call dycoms_rf01_config
      h%init      => dycoms_rf01_init
      h%radiation => dycoms_rf01_radiation
    end subroutine
    ```

    The `register_<name>` routines live in `modcases` rather than in the case modules, so the case modules depend only on the core modules and there is no circular `use`.

3. **Mask solid cells.** A hook that writes a tendency multiplies its contribution by `IIc` (or the matching `IIu`/`IIv`/`IIw`), as the core modules do.

4. **Switch it on** in the case's `namoptions` (see `examples/303/namoptions.303`):

    ```fortran
    &NAMCASES
    case_hooks   = 'dycoms_rf01'
    /

    &NAMCASE_DYCOMS
    dycoms_D        = 3.75e-6
    dycoms_qt_thres = 8.0e-3
    /
    ```

5. **Test and document it.** Add tests next to the feature's integration tests (here `tests/integration/thermo/test_thermo.py`, `TestRadiation` and `TestCaseHooks`), add the case namelist to `docs/udales-namoptions-overview.md` and to `tools/python/namelists.json`.

## When to add a hook point

Keep the set of hook points small. Add a new point only when a case needs it, call it from one place in the core code, document it in the table above and in the header of `modcases.f90`, and keep it a no-op when no hook is registered. A generic, data-driven option (for example `rad_absorber = 'sv'` for a smoke absorber) is preferable to a hook when the physics is the same and only the input differs.
