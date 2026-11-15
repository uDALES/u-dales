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

- **`&NAMTHERMOSTAT` is removed; the diagnostics are `lbldump` in `&OUTPUT`.** The cloud,
  boundary-layer and radiation time series and profiles are now part of the statistics
  machinery (`modstatsdump`, helper `modbldiag`) and written to `bldump_t.<expnr>.nc` and
  `bldump_z.<expnr>.nc` (were `thermostat_t`/`thermostat_z`); variables, units and
  definitions are unchanged. `lthermostat` becomes `lbldump` (default `.false.`) and
  `tthermostat` is replaced by `tstatsdump` (samples every `tstatsdump` from `tstatstart` on,
  at the same point of the time loop); `zi_qt_thres`, `zi_sv_index` and `zi_sv_thres`
  move to `&OUTPUT` with the same defaults and checks. An old `&NAMTHERMOSTAT` group in
  `namoptions` is silently ignored (no diagnostics are written), so set `lbldump = .true.`
  and `tstatsdump` in `&OUTPUT`.
- **The von Karman constant default `fkar` is 0.40 instead of 0.41.** It is now the value used by DALES
  (`fkar = 0.4`), PALM (`kappa = 0.4`) and WRF; the 2.3 % u* gap between uDALES and DALES in the RF01 comparison
  was entirely this constant. `fkar` remains a `&WALLS` input. The inlet generator (`modinlet`) keeps its own
  0.41, which is calibrated together with its Coles wake constants. Affected: every run that does not set `fkar`,
  through the wall functions (`modwallfunctions`, `modibm`), including the Monin-Obukhov facets
  (`iwallmom = 4`, new); for the same roughness length the wall stresses fall by about 5 % (kappa^2 is 0.160 instead of
  0.168). To keep the old behaviour set `fkar = 0.41` in `&WALLS`. Example 303 and the DALES RF01 inputs use the
  recomputed `z0 = z0h = 1.445944e-5 m` (same neutral Cm = 0.0011).
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
  equivalent). Cases with `lrandomize = .false.` are unchanged.

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
  - `ql` itself in field dumps, restart files and the `bldump` files (`ql`, `cfracz`,
    `lwp_bar`, `cfrac`, `qlmax`);
  - the slab means `ql0av`, hence `th0av`, `thvf`, `thvh(kb)`, `presf`/`presh`,
    `exnf`/`exnh` and `rhof`;
  - the moist `dthvdz` of `calthv` (SGS stability and the one-equation buoyancy source);
  - the cloud-water longwave radiation absorber (`iradiation = 2`).
  Dry runs and moist runs that never saturate are unchanged. A warm start from a
  restart file written before this fix carries the shifted `ql0` into the first
  substep's radiation, which is a one-substep transient.
- **Large-scale subsidence acts on the local fields.** `lstend` now applies the
  upwind subsidence term `-whls dphi/dz` to `thl0`, `qt0`, `sv0` (and `u0`, `v0` with
  `lmomsubs`) cell by cell, as DALES dev 075ac4e does, instead of to the slab means
  `thl0av`, `qt0av`, `sv0av`. The mean tendency is the same; the term now also
  advects the fluctuations, which matters at a sharp, wavy inversion (DYCOMS RF01).
  The first level keeps the downward-only treatment (as in DALES). Who is affected:
  runs with a non-zero `wfls` in `lscale.inp`; urban cases (`wfls = 0`) are unchanged
  bit for bit.
- **Latent heat of vaporisation `rlv` is now 2.53e6 J/kg**, the DALES dev 075ac4e
  value. The old 2.26e6 is the value at 100 °C; the atmospheric value is 2.50e6 at
  0 °C and 2.45e6 at 20 °C. 2.53e6 was chosen so that uDALES and DALES use identical
  constants; it is about 3 % above the 20 °C value. It enters the saturation adjustment, `thv`, `dthvdz`, `th0av`
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
    mixing than intended. No user action is needed, whether or not a case sets
    `Prandtl`: a value set in `&NAMSUBGRID` now means `Kh = Km/Prandtl` in all three
    SGS models. Every one-equation run with `Prandtl = 1/3` (the default; no shipped
    case sets another value) gets 9 times more SGS heat and moisture mixing than
    before.
  - The SGS TKE dissipation no longer carries the extra factor 2 ("added factor
    2. for shear-driven flow"); it is DALES's `-(ce1 + ce2 l/Δ) e^2/(2l)`. The
    near-wall `damp` factor (1 away from walls) is kept.
  - The stable length scale and the buoyancy source use `grav/thvf(k)` (slab
    mean virtual potential temperature) instead of `grav/thvs`, as in DALES. In a
    fully solid slab above `kb` (IBM) `thvf(k)` holds the `-999` fill, so those
    levels use the neutral length scale and no buoyancy source; the fluid cells next
    to such a slab see this through the SGS diffusion of `e12`.
  - The TKE sources at the first level `kb`: there are still none next to the wall
    functions (which set the near-wall state), but they are now set to zero
    explicitly; before, uninitialised memory was added to `e12p(kb)`. Above
    horizontal Monin-Obukhov facets (`iwallmom = 4`, new) they follow DALES, from
    the surface-layer shears; see the [numerics reference](udales-dales-numerics.md).
- **Default turbulent Prandtl number is exactly 1/3** (was 0.333), the DALES value.
  `Kh = Km/Prandtl` was `3.003 Km`; it is now `3 Km`, about 0.1 % less SGS heat,
  moisture and scalar diffusivity. Who is affected: every Smagorinsky, Vreman and
  one-equation run that does not set `Prandtl` in `&NAMSUBGRID`.
- **`wqsurf` is positive downward**, like `wtsurf`, `wttop`, `wqtop` and the facet
  fluxes (`modibm:bottom`; the long-term convention is
  [issue #393](https://github.com/uDALES/u-dales/issues/393)). Before, `wqsurf` alone
  was positive upward. Who is affected: runs that prescribe a non-zero `wqsurf`
  (no shipped case did); change the sign of the input to keep the old flux.
- **`thls` and `qts` must be set** ([issue #375](https://github.com/uDALES/u-dales/issues/375)).
  `checkinitvalues` aborts when `thls` is unset with `ltempeq` or `lbuoyancy`, or `qts`
  is unset with `lmoist`; dry runs set `qts = 0`. Before, the -1 sentinels entered
  `thvs`, so dry runs had `thvs = 0.392 thls` instead of `thls`. The one-equation
  stable length scale also floors `|dthvdz|` at `eps1`. Who is affected: runs that
  relied on the sentinels (they now stop at start-up), and dry `ltempeq` runs through
  `thvs`.

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

- `modthermodynamics:thermo` no longer reads the iteration counter `niter` uninitialised
  in the Newton-Raphson saturation adjustment (`lqlnr = .true.`)
  ([issue #400](https://github.com/uDALES/u-dales/issues/400)). `niter` is a diagnostic
  counter only (`nitert` is not used), so results are unchanged.
- `modthermodynamics:calthv`: in a saturated cell where the dry and saturated gradient
  coefficients coincide the moist `dthvdz` divided by zero; it now uses the saturated
  value. Only affects `lmoist = .true.` runs with exactly uniform saturated cells.
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

- **Monin-Obukhov surface layer on facets** (`iwallmom = 4` in `&WALLS`, with the new inputs
  `mo_maxslope` and `ustin`; `modsurface:mo_local_update`, `modibm:wallfunmom`), ported from DALES
  dev 075ac4e (`isurf = 4`, or `isurf = 3` with a prescribed friction velocity `ustin > 0`) with the
  Obukhov length computed locally (DALES default `lmostlocal = .true.`), per facet section. Facets whose
  normal is within `mo_maxslope` (default 30 degrees) of the vertical use it, from the air state at the
  section's reference point (the interpolation of the other wall functions), the facet `z0`/`z0h` and the
  prescribed facet fluxes `bctfz` and `bcqfz` (the latter on facets whose type has `lGR = 1`, i.e. that
  exchange moisture); such a facet takes `bctfz`/`bcqfz` whatever its slope. Steeper facets and walls use the
  neutral log-law wall function (as `iwallmom = 2` does under `iwalltemp = 1`). Requires `libm`, `ltempeq`,
  `iwalltemp = 1`, `iwallmoist = 1` with `lmoist`, and `z0, z0h > 0` on the Monin-Obukhov facets; a start-up
  line counts the Monin-Obukhov and wall-function facets. `bldump` writes the area means over the
  Monin-Obukhov sections (`ustar`, `obl` as `1/mean(1/L)`, the diagnosed surface `thls`, `qts`). The
  per-section state is not in the restart files (cold-started again on a warm start). In the cells directly
  above horizontal Monin-Obukhov facets (normal +z) the one-equation SGS TKE sources and `dthvdz` use the
  Monin-Obukhov surface-layer gradients, the DALES first-level treatment (`modibm:mo_cells_update`, evaluated
  before `subgrid` as in DALES); sloped Monin-Obukhov facets and all other walls keep no TKE source at `kb`.
  Off by default (`iwallmom = 2`); other runs are bit-identical. `tools/python/udprep/flat_floor.py` writes the IBM
  inputs of a flat floor made of facets.
- **Parametrised longwave radiation** (`&NAMRADIATION`, `iradiation = 2`; `modradiation`,
  `modraddata`), DALES `radpar`: cloud-water absorber `rad_absorber = 'ql' | 'sv'`
  (`rad_absorber_sv`), optics `rka`, `dlwtop`, `dlwbot` and `timerad`. `rka` and `dlwtop` (or `dlwbot`)
  have no defaults (0) and `iradiation = 2` aborts without them. Off by default (`iradiation = 0`).
- **Case hooks** (`src/modcases.f90`, `&NAMCASES case_hooks`): case-specific physics
  lives in `src/cases/` and is switched on by name at run time instead of by switches in
  the core modules; see `docs/udales-case-hooks.md`. Hook points: init, radiation,
  forcing, stats; without hooks the results are bit-identical. First hook:
  `dycoms_rf01` (`&NAMCASE_DYCOMS dycoms_D, dycoms_qt_thres`), the DYCOMS-II RF01
  free-troposphere radiation term, computed as the divergence of the analytic flux between cell
  faces (finite at the inversion, conserves the column integral).
- **`bldump` diagnostics** (`lbldump`, `tstatsdump`, `zi_qt_thres`, `zi_sv_index`, `zi_sv_thres` in
  `&OUTPUT`; `modbldiag`, called from `modstatsdump`): `bldump_t.<expnr>.nc` (liquid water path, cloud
  cover, boundary-layer heights, `wmax`, `qlmax`, `ustar`, `obl`, `thls`, `qts`, first-level `u`, `v`) and
  `bldump_z.<expnr>.nc` (slab-mean profiles). Instantaneous samples every `tstatsdump` from `tstatstart`
  on, at the end of the time step, after
  `thermodynamics`, so `ql`, `rhof`, `exnf` and the slab means belong to the same state as
  `thl`, `qt` and `w`. The tracer boundary-layer heights `zi_sv` (`zi_sv_index`, `zi_sv_thres`) and
  `zi_qt` (`zi_qt_thres`) are written only when set. Off by default; the solution itself is unaffected.
- `fieldvars = 'dt'` writes `dthvdz` to the field dump.
- **Initial `thl` and `qt` noise** (`randthl`, `randqt`; see Breaking changes, issue #395).
- **Examples 301 (BOMEX), 302 (smoke cloud, scalar radiation absorber) and 303 (DYCOMS-II RF01)**, on a
  floor made of one horizontal facet with the Monin-Obukhov facet surface (`iwallmom = 4`), with
  the DALES comparison tooling in `tools/python/examples/dales_cases` (input converter,
  `compare_dales.py`, DALES input sets and patches); results in each example's `RESULTS.md`.
- **Tests**: in-solver runmodes 1008 (`tests_sgs_oneeqn`: one-equation SGS closure and TKE
  sources), 1009 (`tests_lstend`), 1010 (`tests_thermo_levels`: `ql0` level, #403), 1011 (cold
  start, write the initial fielddump record and stop), 1012 (`tests_mo_kernels`:
  Monin-Obukhov kernels and the per-section surface update), driven by `tests/integration/thermo/test_thermo.py` (includes the
  1-rank against multi-rank check of the e12 halos, and the Monin-Obukhov facet floor: u* and L
  against the analytic values in a neutral and an unstable case, the applied stress, prescribed u*, and
  the slope criterion).
- **Module hygiene for the new modules**: `modsurface`, `modradiation`, `modraddata`,
  `modbldiag`, `modcases` and `case_dycoms_rf01` are private by default with explicit
  public lists; configuration and namelist read errors in them abort all ranks through
  `modmpi:abort_all_ranks`, after one message from rank 0 (this includes the new start-up
  checks in `checkinitvalues`); `iradiation = 2` also rejects negative `timerad`, `dlwtop`,
  `dlwbot`; a case namelist group whose hook is not in `case_hooks` gives a warning.

- `docs/udales-output-files.md` now documents the SGS flux definitions, their
  sign convention, and the fact that the SGS statistics carry **no** surface or
  wall-model contribution, so the total-stress profile does not close at the
  wall ([issue #353](https://github.com/uDALES/u-dales/issues/353)).
- New experimental test suite `tests/integration/sgs_statistics` (`runmode`
  1006) validating the SGS fluxes against a manufactured solution on a
  stretched vertical grid, against the solver's own diffusion tendencies, and
  against the sign convention, on a new fixture `tests/cases/300`.
