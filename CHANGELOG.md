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

- **Execute scripts refuse an existing output directory.** `local_execute.sh`,
  `hpc_execute.sh`, `archer_execute.sh` and `ecmwf_execute.sh` (and so
  `ud_run <machine> sim`) now exit with an error when `$DA_WORKDIR/<expnr>`
  already exists, instead of copying the inputs and the executable over it and
  writing the new run next to the old results. To rerun a case, move or delete
  its output directory yourself first. `copy_inputs.sh` no longer creates the
  (empty) output directory of the new case; cases prepared with the previous
  version have an empty output directory that must be removed before the first
  run.

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

- **Prescribed facet heat and moisture flux (`iwalltemp = 1`, `iwallmoist = 1`).**
  `wallfunheat` only assigned the prescribed flux (`bctf*`, `bcqf*`) to facets
  facing exactly `+x`, `-x`, `+y`, `-y` or `+z`; any other facet (sloped, or
  facing down) was left with the flux of the previously processed section, or
  an undefined value, which also entered the total heat and moisture flux and
  the facet energy balance. Such facets now get the flux of their facet type
  when the new optional `factypeflux.inp.<expnr>` lists it (`wallid`, heat flux
  [K m/s], moisture flux [kg/kg m/s]; `factypes.inp` is unchanged and still read
  as before), otherwise `bctf*`/`bcqf*` by direction as before, otherwise zero,
  and a start-up line reports how many facets have no prescribed flux when a
  prescribed flux (including `ltimedepsurf`) is nonzero. The heat and moisture
  fluxes of a section now use separate variables, and the remaining unset paths
  are defined: with `iwallmoist = 2` the moisture flux is zero where the heat
  transfer coefficient times the tangential velocity is zero.
  `iwalltemp` or `iwallmoist` other than 1 or 2, and `iwallmoist = 2` without
  `ltempeq = .true.` and `iwalltemp = 2`, now stop at start-up. Python
  pre-processing writes the new file with `IBMSection.write_factypeflux`.
  Who is affected: (1) runs with non-grid-aligned or downward-facing facets and
  a nonzero `bctf*`/`bcqf*` change (those facets previously got an arbitrary
  flux); (2) runs with `iwalltemp = 2`, `lmoist` and the default `iwallmoist = 1`
  change on green facets that are sloped or face down: they previously received
  the section's wall-function heat flux as a moisture flux, now zero (or the
  `bcqf*`/`factypeflux.inp` value); (3) `iwalltemp = 1` with `iwallmoist = 2`
  no longer lets the moisture flux leak into the heat flux of unmatched facets,
  and now stops at start-up because `iwallmoist = 2` needs `iwalltemp = 2`;
  (4) Debug builds that trap uninitialised reals aborted on these facets.
  Runs with `iwalltemp = 1`, `iwallmoist = 1` (or no moisture) and all
  prescribed fluxes zero are unchanged. See #409.
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

- `docs/udales-output-files.md` now documents the SGS flux definitions, their
  sign convention, and the fact that the SGS statistics carry **no** surface or
  wall-model contribution, so the total-stress profile does not close at the
  wall ([issue #353](https://github.com/uDALES/u-dales/issues/353)).
- New experimental test suite `tests/integration/sgs_statistics` (`runmode`
  1006) validating the SGS fluxes against a manufactured solution on a
  stretched vertical grid, against the solver's own diffusion tendencies, and
  against the sign convention, on a new fixture `tests/cases/300`.
