# Changelog

Notable changes to uDALES, newest first. This project uses
[semantic versioning](https://semver.org/); released versions are also tagged on
[GitHub](https://github.com/uDALES/u-dales/releases).

Entries under **Unreleased** are on `master` but not yet in a tagged release.
Changes that alter the meaning, sign, units or location of an output variable,
or that remove or disable a feature, are listed under **Breaking changes** so
that they are easy to find when upgrading.

## Unreleased

### Breaking changes

- **The SGS TKE halos are exchanged between MPI ranks** ([issue #405](https://github.com/uDALES/u-dales/issues/405)).
  `modboundary:halos` never called `exchange_halo_z` for `e120` and `e12m`, so with the
  one-equation closure (`loneeqn = .true.`) and `nprocx > 1` or `nprocy > 1` the e12 halo
  cells on inter-rank boundaries kept their initial values for the whole run. The TKE
  diffusion (`diffe`) and the horizontal averages of the closure used them, giving a
  permanent spurious TKE source at the tile edges. Both fields are now exchanged when
  `loneeqn` is set. Affected: every multi-rank run with `loneeqn = .true.` (results change).
  Not affected: Smagorinsky and Vreman runs, and 1-rank runs, which are bit-identical.
- **Initial perturbations use a decomposition-independent hash, and thl/qt noise is
  available** ([issue #395](https://github.com/uDALES/u-dales/issues/395)).
  `modstartup:randomize_field` added one step of a linear congruential generator, seeded
  with the cell's linear index, so the noise was a sawtooth plane wave (58 % of its
  variance at a single wavenumber), and u, v and w received the same number at every cell.
  It is replaced by a splitmix64-style counter-based hash of
  `(irandom, field, iglob, jglob, k)`: uniform, white, independent between fields, and still
  independent of the MPI decomposition. `randthl` and `randqt` (previously ignored,
  default 0 and still 0) now perturb `thl` and `qt` (the m and 0 fields alike, up to `krand`,
  qt capped at its local value so that it stays non-negative). Affected: every case with
  `lrandomize = .true.` (the default) gets different initial velocity noise of the same
  amplitude, so results are no longer bit-identical to earlier versions (statistically
  equivalent). Cases with `lrandomize = .false.` are unchanged. New test runmode 1011 writes
  the initial fielddump record and stops.

- **Cloud liquid water `ql0` is now stored at the right level**
  ([issue #403](https://github.com/uDALES/u-dales/issues/403)). `modthermodynamics:thermo`
  declared its `ql` argument from `kb`, while `ql0` is allocated from `kb-1`, so the
  full-level saturation adjustment of level `k+1` was stored in `ql0(k)` (since v2.0).
  The cloud therefore sat one level too low with respect to `thl` and `qt`. The
  half-level `ql0h` (and so the buoyancy term through `thv0h`) was computed at the
  right level, and changes only indirectly, through `presh`/`exnh` (about 1e-4
  relative). Now
  fixed: `thermo` takes assumed-shape arguments with the actual lower bounds. Every
  run with `lmoist = .true.` and any saturated cell changes, through these routes:
  - `ql` itself in field dumps, restart files and the thermostat (`ql`, `cfracz`,
    `lwp_bar`, `cfrac`, `qlmax`);
  - the slab means `ql0av`, hence `th0av`, `thvf`, `thvh(kb)`, `presf`/`presh`,
    `exnf`/`exnh` and `rhof`;
  - the moist `dthvdz` of `calthv` (SGS stability and the one-equation buoyancy source);
  - the cloud-water longwave radiation absorber (`iradiation = 2`).
  Dry runs and moist runs that never saturate are unchanged. A warm start from a
  restart file written before this fix carries the shifted `ql0` into the first
  substep's radiation, which is a one-substep transient.
- **Monin-Obukhov surface state is computed before the SGS model.** With
  `BCbotm = 4/5`, the Obukhov length, `u*`, the surface-layer gradients and the
  re-diagnosed `thls`/`qts`/`thvs` are now computed at the start of each RK substep
  (`modsurface:surface_mo_state`, called before `advection`/`subgrid`), as in DALES
  dev 075ac4e (`surface` before `advection` and `subgrid`). Before, they were computed
  in `bottom`, after `subgrid`, so the one-equation SGS TKE sources at `kb` used the
  surface-layer shear `dudz`/`dvdz` of the previous substep (and zero on the first
  substep). `dthvdz(kb)` keeps its timing: `calthv` sets it at the end of the
  previous substep, as in DALES. `bottom` still applies the stress. On a cold start
  the MO state is now also computed during startup, before a second
  `thermodynamics` call, so the first substep already has the MO `dthvdz(kb)`
  (DALES `readinitfiles` does the same). Who is affected: runs with `BCbotm = 4/5`;
  with `loneeqn = .true.` every substep changes; with Vreman or Smagorinsky the
  change comes only from the cold-start initial state (the time-loop reorder alone
  was verified bit-identical for them). All other runs are unchanged bit for bit.
- **Thermostat diagnostics are sampled after `thermodynamics`.** `thermostat` now
  runs at the end of the time step, after `boundary` and `thermodynamics`, so `ql`,
  `rhof`, `exnf` and the slab means it writes belong to the same state as `thl`,
  `qt` and `w` (DALES samples its time series after thermodynamics too). Before,
  `ql0` and `rhof` were one RK substep old. The output times are unchanged. Who is
  affected: only runs with `lthermostat = .true.`; the solution itself is unchanged.
- **Large-scale subsidence acts on the local fields.** `lstend` now applies the
  upwind subsidence term `-whls dphi/dz` to `thl0`, `qt0`, `sv0` (and `u0`, `v0` with
  `lmomsubs`) cell by cell, as DALES dev 075ac4e does, instead of to the slab means
  `thl0av`, `qt0av`, `sv0av`. The mean tendency is the same; the term now also
  advects the fluctuations, which matters at a sharp, wavy inversion (DYCOMS RF01).
  The first level keeps the downward-only treatment (as in DALES). Who is affected:
  runs with a non-zero `wfls` in `lscale.inp`; urban cases (`wfls = 0`) are unchanged
  bit for bit.
- **Latent heat of vaporisation `rlv` is now 2.53e6 J/kg** (was 2.26e6), the DALES
  dev 075ac4e value. It enters the saturation adjustment, `thv`, `dthvdz`, `th0av`
  and the derived constants `rlvi` and `rlvocp`, and through them the latent heat
  flux of green-roof facets in the energy balance (`modEB:updateGR`; the soil-water
  budget uses `rlv*rlvi` and is unchanged) and the psychrometric constant and
  moisture sources of the tree model (`vegetation`, `lmoist` and `ltempeq`). For
  the same `thl` and `qt`, cloud liquid water is about 10-15 % lower than before
  (more in warmer clouds).
  Who is affected: runs with `lmoist = .true.` that form cloud, moist runs with
  trees, and energy-balance runs with green roofs. Dry runs are unchanged. Case 303
  (DYCOMS RF01) now prescribes `wqsurf = -3.72578241e-5` (115 W/m2 / (1.22 x 2.53e6)).
- **One-equation SGS model (`loneeqn`) aligned with DALES.** Four changes to
  `modsubgrid.f90`, following DALES dev 075ac4e (#401). Who is affected: only runs with
  `loneeqn = .true.` (Vreman and Smagorinsky runs change only through the new
  `Prandtl` default, see the next entry).
  - The turbulent Prandtl number now means `Kh = Km/Prandtl`, as in DALES and as
    the Smagorinsky and Vreman models already used it (`ch = 1/Prandtl`). With
    the default `Prandtl = 1/3`, `Kh = 3 Km` in unstable and neutral cells,
    going to `Km` in strongly stable ones. The model previously used
    `ch = Prandtl`, which gave `Kh = Km/3`: 9 times less SGS heat and moisture
    mixing than intended. Cases that set `Prandtl` for the one-equation model
    must invert it.
  - The SGS TKE dissipation no longer carries the extra factor 2 ("added factor
    2. for shear-driven flow"); it is DALES's `-(ce1 + ce2 l/Δ) e^2/(2l)`. The
    near-wall `damp` factor (1 away from walls) is kept.
  - The stable length scale and the buoyancy source use `grav/thvf(k)` (slab
    mean virtual potential temperature) instead of `grav/thvs`, as in DALES.
  - The TKE sources at the first level `kb`: with the Monin-Obukhov surface
    (`BCbotm = 4/5`) they are now computed as in DALES's k = 1 block, using the
    MO surface-layer shears `dudz`, `dvdz` and `dthvdz(kb)`. With the other
    bottom conditions there are still no sources at `kb`, but they are now set
    to zero explicitly; before, uninitialised memory was added to `e12p(kb)`.
- **Default turbulent Prandtl number is exactly 1/3** (was 0.333), the DALES value.
  `Kh = Km/Prandtl` was `3.003 Km`; it is now `3 Km`, about 0.1 % less SGS heat,
  moisture and scalar diffusivity. Who is affected: every Smagorinsky, Vreman and
  one-equation run that does not set `Prandtl` in `&NAMSUBGRID`.

- **SGS flux sign convention.** All subgrid-scale fluxes in the statistics
  output (`usgs*`, `vsgs*`, `wsgs*`, `thlsgs*`, `qtsgs*`, `sv1sgs`–`sv4sgs`,
  `sca1sgs*`–`sca3sgs*`) were previously stored as `+K dphi/dn`, i.e. **minus**
  the physical flux, and were therefore opposite in sign to the resolved fluxes
  written into the same files with the same units. They are now stored as the
  flux itself, `-K dphi/dn`. Any post-processing that compensated for the old
  convention — by negating the SGS term, or by subtracting rather than adding it
  when forming a total flux — must be updated. See
  [Output files](udales-output-files.md) for the explicit definitions.
- **`wsgs` is now a stress, not a tendency.** `wsgs` (and `wsgsyt` in `ytdump`)
  previously held the vertical part of the `w`-diffusion tendency, in `m/s^2`,
  stored at the w-level. It now holds the normal stress
  `tau_33 = -2 K_m dw/dz` in `m^2/s^2`, stored at the **cell centre**. The
  NetCDF grid code changed from `t0mt` to `t0tt` accordingly, so the variable's
  vertical coordinate changed from `zm` to `zt`.
- **`ltkedump` is disabled.** The TKE-budget path never accumulated its budget
  terms, so `tkedump.<expnr>.nc` was written full of zeros. Setting
  `ltkedump = .true.` now aborts the run with an explanatory error rather than
  producing a plausible-looking but meaningless file. Repairing the path is
  tracked in [issue #352](https://github.com/uDALES/u-dales/issues/352).
- **Masked bottom level in SGS profiles.** The w-located xy-averaged profiles
  (`wxy*`, `usgsxy*`, `vsgsxy*`, `thlsgsxy*`) now report `-999.` at the bottom
  level `kb` in IBM runs, consistent with every other fully masked level.
  Previously that level held an unmasked sum over solid points divided by the
  fluid-point count at the top of the domain.
- **Units and metadata corrections** in the statistics NetCDF headers:
  `vsgsxyt` and `vpwpxy` are now `m^2/s^2` (were `K m/s`/`Km/s`); `sv1sgs`–
  `sv4sgs` are now `g/m^2 s` (were `gm/s`); `wqtyt` is now `kg/kg m/s` and
  `wsca1yt`–`wsca3yt` are now `kg/m^2 s` (all were `K m/s`). Placeholder
  `long_name` strings in the TKE-budget header were replaced with the actual
  term names.

### Fixed

- One-equation SGS model: in a fully solid slab above `kb` the slab-mean `thvf(k)`
  holds the `-999` fill, so `grav/thvf(k)` gave a negative buoyancy parameter (a NaN
  length scale where `dthvdz > 0`, a wrong-sign buoyancy source). Those levels now use
  the neutral length scale and no buoyancy source. The values change in the solid
  cells, but the fluid cells directly above or below such a slab see the change
  through the SGS diffusion of `e12` (`diffe` has no IBM mask); before, a NaN could
  reach them the same way.
- The SGS statistics no longer read uninitialised data at the top level
  `ke+kh`, and all SGS output arrays are zeroed before use.
- The 3-D `sv1sgs`–`sv4sgs` output in `tdump` is now masked with `IIw`, so solid
  and wall-adjacent faces report `-999.` instead of an interior-formula value
  the solver never applied.
- `upwpxy`, `vpwpxy` and `wpthlpxy` report `-999.` at masked levels instead of
  the product of two sentinels (`-999000`).
- Several latent bugs inside the (disabled) TKE-budget path were repaired so
  that it is correct if revived: the mean rate-of-strain was a scalar reused
  across the whole domain, the vertical SGS transport term paired the `w`-based
  fluxes with the mean `v`, a shear-production term was double-counted, and a
  four-point average varied the wrong index.

### Added

- **Case hooks** (`src/modcases.f90`, `&NAMCASES case_hooks`): case-specific physics
  lives in `src/cases/` and is switched on by name at run time instead of by switches in
  the core modules; see `docs/udales-case-hooks.md`. Hook points: init, radiation,
  forcing, stats; without hooks the results are bit-identical. First hook:
  `dycoms_rf01` (`&NAMCASE_DYCOMS dycoms_D, dycoms_qt_thres`), the DYCOMS-II RF01 free-troposphere
  radiation term, now computed as the divergence of the analytic flux between cell
  faces (finite at the inversion, conserves the column integral). The generic radiation
  absorber `rad_absorber = 'ql' | 'sv'` (`rad_absorber_sv`) and the thermostat
  diagnostics `zi_sv` (`zi_sv_index`, `zi_sv_thres`) and `zi_qt` (`zi_qt_thres`), each
  written only when set, replace the smoke switches and the hard-coded `zi_qt8`. The
  radiation optics `rka` and `dlwtop` have no defaults (0) and `iradiation = 2` aborts
  without them. The removed `lrad_dycoms`, `rad_dycoms_D`, `rad_dycoms_qtthres`,
  `rad_smoke`, `isvsmoke`, the `zi_smoke`/`zi_qt8` outputs and the former
  `rka = 130`/`dlwtop = 74` defaults were never released.

- **Engineering hygiene of the modules added on this branch** (issue #402 applied to the new
  code only): `modsurface`, `modradiation`, `modraddata`, `modthermostat`, `modcases` and
  `case_dycoms_rf01` are private by default with explicit public lists; configuration errors
  and namelist read errors in them abort all ranks through `modmpi:abort_all_ranks`, after
  one message from rank 0 (this includes the new start-up checks in `checkinitvalues` and the
  `fieldvars = 'dt'` check); `tthermostat` has no default (0) and `lthermostat = .true.`
  aborts without it (it was 60 s, never released); `iradiation = 2` also rejects negative
  `timerad`, `dlwtop`, `dlwbot`; a case namelist group whose hook is not in `case_hooks`
  gives a warning; `randomize_field` no longer aliases its inout field with the qt cap
  (results unchanged to the bit). New test runmode 1012 (`tests_mo_kernels`) checks the
  Monin-Obukhov kernels against reference values. Bit-identical results otherwise.

- In-solver test `runmode` 1008 (`tests_sgs_oneeqn`, run from
  `tests/integration/thermo/test_thermo.py:TestSGSOneEqn`) checking the
  one-equation SGS closure (`Kh/Km`, stable length scale) and its TKE sources.

- `docs/udales-output-files.md` now documents the SGS flux definitions, their
  sign convention, and the fact that the SGS statistics carry **no** surface or
  wall-model contribution, so the total-stress profile does not close at the
  wall ([issue #353](https://github.com/uDALES/u-dales/issues/353)).
- New experimental test suite `tests/integration/sgs_statistics` (`runmode`
  1006) validating the SGS fluxes against a manufactured solution on a
  stretched vertical grid, against the solver's own diffusion tendencies, and
  against the sign convention, on a new fixture `tests/cases/300`.
