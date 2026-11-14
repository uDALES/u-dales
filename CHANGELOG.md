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

- **Monin-Obukhov surface state is computed before the SGS model.** With
  `BCbotm = 4/5`, the Obukhov length, `u*`, the surface-layer gradients and the
  re-diagnosed `thls`/`qts`/`thvs` are now computed at the start of each RK substep
  (`modsurface:surface_mo_state`, called before `advection`/`subgrid`), as in DALES
  dev 075ac4e (`surface` before `advection` and `subgrid`). Before, they were computed
  in `bottom`, after `subgrid`, so the one-equation SGS sources at `kb` and
  `dthvdz(kb)` used the gradients of the previous substep. `bottom` still applies
  the stress. Who is affected: only runs with `BCbotm = 4/5`; all other runs are
  unchanged bit for bit.
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
  the same `thl` and `qt`, cloud liquid water is about 12-14 % lower than before.
  Who is affected: runs with `lmoist = .true.` that form cloud, moist runs with
  trees, and energy-balance runs with green roofs. Dry runs are unchanged. Case 303
  (DYCOMS RF01) now prescribes `wqsurf = -3.72578241e-5` (115 W/m2 / (1.22 x 2.53e6)).
- **One-equation SGS model (`loneeqn`) aligned with DALES.** Four changes to
  `modsubgrid.f90`, following DALES dev 075ac4e (#401). Who is affected: only runs with
  `loneeqn = .true.`; Vreman and Smagorinsky runs are unchanged (bit for bit).
  - The turbulent Prandtl number now means `Kh = Km/Prandtl`, as in DALES and as
    the Smagorinsky and Vreman models already used it (`ch = 1/Prandtl`). With
    the default `Prandtl = 0.333`, `Kh = 3 Km` in unstable and neutral cells,
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
  the neutral length scale and no buoyancy source. Only solid cells are affected.
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
