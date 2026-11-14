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

- `initcases` reads `&NAMCASES` and splits `case_hooks` at commas. For each name, `register_case` calls the matching `register_<name>` routine (a `select case`). An unknown name aborts and lists the valid names; a name listed twice aborts.
- A `register_<name>` routine fills one `hook_set`, a derived type with one procedure pointer per hook point. Fortran has no arrays of procedure pointers, so the registry keeps an array of these derived types instead.
- The core code calls `call_hooks_<point>`, which runs the hook of every registered case that set that point.

| Hook point | Called from | Typical use |
| ---------- | ----------- | ----------- |
| `init` | `program.f90`, once after startup (after `readinitfiles`, so for cold and warm starts) | read the case namelist, check the configuration |
| `radiation` | `modradiation:radiation`, after the generic `radpar` tendency is computed and before `thlprad` is added to `thlp` (only with `iradiation /= 0`, at the radiation update) | add to `thlprad` |
| `forcing` | the time loop, after `lstend` | add to the tendencies (`up`, `thlp`, ...) |
| `stats` | the time loop, after `thermostat` | write case-specific diagnostics |

Hooks take no arguments: they use the model state through `use` statements, like any other module.

## Worked example: DYCOMS-II RF01

`src/cases/case_dycoms_rf01.f90` adds the free-troposphere term of the Stevens et al. (2005) radiation to the generic longwave scheme.

1. **Write the case module** in `src/cases/` (CMake compiles every `src/**/*.f90`; re-run the configure step after adding a file). Keep it private by default, with an explicit public list of the hook procedures:

    ```fortran
    module case_dycoms_rf01
      implicit none
      private
      public :: dycoms_rf01_init, dycoms_rf01_radiation
      real :: D        = 0.   !< required: no nonzero defaults
      real :: qt_thres = 0.
    contains
      subroutine dycoms_rf01_init        ! reads &NAMCASE_DYCOMS, aborts if D or qt_thres unset
      subroutine dycoms_rf01_radiation   ! thlprad(i,j,k) -= IIc * (F(zh(k+1)) - F(zh(k))) / dzf(k)
    end module
    ```

    The case has its own namelist, `&NAMCASE_DYCOMS`, read from the same `namoptions` file in the `init` hook (rank 0 reads, then `MPI_BCAST`), with the usual error handling. Parameters that are not obviously nonzero default to 0 and the `init` hook aborts with a clear message when they are not set; the case's `namoptions` sets them explicitly.

2. **Register it** in `src/modcases.f90`: add the name to `valid_names`, a branch to the `select case` in `register_case`, and a `register_<name>` routine that sets the pointers of the hook points the case uses:

    ```fortran
    case ('dycoms_rf01')
      call register_dycoms_rf01(hooks(nhooks))
    ...
    subroutine register_dycoms_rf01(h)
      use case_dycoms_rf01, only : dycoms_rf01_init, dycoms_rf01_radiation
      type(hook_set), intent(inout) :: h
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
    D            = 3.75e-6
    qt_thres     = 8.0e-3
    /
    ```

5. **Test and document it.** Add tests next to the feature's integration tests (here `tests/integration/thermo/test_thermo.py`, `TestRadiation` and `TestCaseHooks`), add the case namelist to `docs/udales-namoptions-overview.md` and to `tools/python/namelists.json`.

## When to add a hook point

Keep the set of hook points small. Add a new point only when a case needs it, call it from one place in the core code, document it in the table above and in the header of `modcases.f90`, and keep it a no-op when no hook is registered. A generic, data-driven option (for example `rad_absorber = 'sv'` for a smoke absorber) is preferable to a hook when the physics is the same and only the input differs.
