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

- **Output files and variables renamed; output modules rewritten.** The
  statistics and instantaneous outputs are now written by `out_stats.f90`
  (time- and space-averaged statistics), `out_instants.f90` (fields, slices,
  probes) and `out_write_nc.f90`; `modstatsdump.f90` keeps only the legacy
  `mintdump` and the disabled `tkedump` paths. Files are written **one per
  x-rank** (`<name>.<px>.<expnr>.nc`, covering the rank's whole y-extent)
  instead of one per 2-D pencil; xz-plane output is written per y-rank and
  fully reduced output once. `gather_outputs.sh` was rewritten for the new
  layout. File names:

  | before | now |
  |---|---|
  | `tdump.<px>.<py>.<expnr>.nc` | `stats_t.<px>.<expnr>.nc` |
  | `xytdump.<expnr>.nc` | `stats_xyt.<expnr>.nc` |
  | `xydump.<expnr>.nc` | `stats_xy.<expnr>.nc` |
  | `ytdump.<expnr>.nc` | `stats_yt.<px>.<expnr>.nc` |
  | `ydump.<expnr>.nc` | `stats_y.<px>.<expnr>.nc` |
  | `treedump.<px>.<py>.<expnr>.nc` | `stats_tree.<px>.<expnr>.nc` |
  | `fielddump.<px>.<py>.<expnr>.nc` | `ins_field.<px>.<expnr>.nc` |
  | `islicedump`/`jslicedump`/`kslicedump.<px>.<py>.<expnr>.nc` | `ins_islice.<px>`, `ins_jslice.<py>`, `ins_kslice.<px>.<expnr>.nc` |

  Variable names no longer encode the averaging; the file name does, and the
  same quantity has the same name in every file: `ut`/`uxyt`/`uyt`/`uy` →
  `u`; `upwpt`/`upwpxyt` → `upwp`; `uwxyt` → `uw`; `upuptc` → `upup`;
  `tketc` → `tke`; `thlt` → `thl`; `wpthlpt`/`wthlxyt` → `wpthlp`/`wthl`;
  `usgsxyt` → `usgs`; scalars `sca1t` → `s1`, `sca1psca1pt` → `s1ps1p`,
  `wpsca1pt` → `wps1p`, `sv1sgs` → `s1sgs`; instantaneous `sca1` → `s1`;
  slices `u_kslice` → `u`. Earlier entries in this changelog refer to the
  variables by their former names.
- **`&OUTPUT` namelist.** `tfielddump` is replaced by `tinstantdump`, the
  interval of *all* instantaneous output (fields, slices and probes), and
  `tinstantstart` sets the start of that schedule; a namelist that still
  contains `tfielddump` fails to read. `kslice`, `islice` and `jslice` are
  now lists of `nkslice`, `nislice`, `njslice` planes (several planes per
  direction) and `slicevars` selects the slice variables. `tcheck` moved from
  `&NAMCHECKSIM` to `&OUTPUT`. All example and test cases are updated.
- **Time-averaged statistics on an absolute schedule.** A warm start no
  longer shortens an averaging window: the accumulators restart from zero and
  the next record is written `tstatsdump` seconds later. `stats_xy` and
  `stats_y` are instantaneous spatial averages written every `tsample`, as
  before.
- **Fatal output errors abort all ranks** (`MPI_ABORT`) instead of stopping
  rank 0 and leaving the others waiting.
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

- Cell-centre velocities used by the statistics (`wc`) are computed
  correctly, and the statistics flag gating and window resets are consistent
  after every dump and across restarts.
- The chemistry update of `lchem` is kept: since the fusion of the scalar
  update loop the `svm = sv0` copy ran before `chem`, discarding the
  NO/NO2/O3 update every step (#427). A regression test covers the path.
- **The geostrophic pressure gradient is now applied only with `lcoriol = .true.`**
  (#410). `modstartup` added `om23_gs*vg` and `-om23_gs*ug` to the large-scale
  pressure gradient whether or not the Coriolis force was on, so a run with
  nonzero `ug`/`vg` in `lscale.inp` and `lcoriol = .false.` (and `lprofforc =
  .false.`) was accelerated at a constant `f*|u_g|` (about 1.15e-4 s^-1 times
  `|u_g|` at the default `xlat = 52`). Such runs change; they no longer get
  this spurious body force. Runs with `lcoriol = .true.`, with `ug = vg = 0`
  (the default), or with `lprofforc = .true.` are bit-identical.
- **The f2py preprocessing modules run in double precision.** The IBM
  (`ibm_preproc_f2py`) and direct-shortwave (`directshortwave_f2py`) extension
  modules were compiled with default 4-byte reals while the solver uses
  `-fdefault-real-8`. Geometry, grid coordinates and facet sections were
  therefore computed with about seven significant digits. Facet sections
  lying exactly on a cell face could be assigned to either neighbouring cell
  depending on compiler rounding, small sections near cell faces were lost
  to the 1e-5 m^2 area threshold, and, more seriously, the first cell edge
  (`xt(1) - dx/2`) evaluated in single precision exceeded the 1e-8 bounding-box
  tolerance, so facets touching the `x = 0` or `y = 0` domain boundary were
  rejected as out of bounds and received no facet sections at all (in the
  `997`/`988` test geometry, 105 ground facets covering 5% of the surface).
  Both modules are now built with `-fdefault-real-8 -fdefault-double-8`, their
  signatures declare `real(kind=8)`, and the Python callers pass float64 arrays.
  The flux-point selection in `matchFacetsCells` now treats candidate scores
  that agree to within a relative 1e-6 as tied and takes the first, as the
  MATLAB implementation's `max()` does; in double precision, mirror-image
  candidates (the cells above and below a section on a horizontal cell face)
  otherwise differ in the last bit and the choice was decided by round-off.
  Rebuild the preprocessing tools
  (`tools/build_preprocessing.sh <system> preprocessing_tools`) after updating.
- **Preprocessing x/y grids are built by index, like the solver.** `udprep`
  built `xm`, `xt`, `ym`, `yt` with `np.arange` and a float step, which returns
  `itot+1` points for some `xlen/itot` pairs (for example 120 cells over 5.5 m)
  and so put an extra face at the domain top. They now come from
  `udgrid.horizontal_grid`, the same `(i-1)*dx`, `+dx/2` construction as
  `modglobal.f90`, as the z-grid already did since #388.
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

- **Time-averaged statistics on selected planes.** `ltislicedump`,
  `ltjslicedump`, `ltkslicedump` write the full `stats_t` variable set only on
  the planes `islice(1:nislice)`, `jslice(1:njslice)`, `kslice(1:nkslice)`
  (`stats_tislice`, `stats_tjslice`, `stats_tkslice`), independently of
  `ltdump`.
- **Point probes.** `lprobedump` writes the `probevars` variables at the
  `nprobe` grid points listed in `probe.inp.<expnr>` (one header line, then
  `i j k` per line) to `ins_probe.<expnr>.nc`; `copy_inputs.sh` copies
  `probe.inp`.
- **Gapped averaging schedule.** `tstatsgap > 0` writes the time-averaged
  statistics at `tstatstart + n*(tstatsdump + tstatsgap)`, each record the
  mean of the `tstatsdump` seconds preceding it, with no sampling in the gap.
  The schedule is kept across warm starts.
- **New output quantities.** Air temperature `tha` (`ltempeq`), relative
  humidity `rh` (`ltempeq` and `lmoist`), hydrostatic absolute pressure
  `pabs` and liquid water `ql` in the statistics files and, through the
  `ta`, `rh`, `pa`, `ql` labels of `fieldvars`/`slicevars`/`probevars`, in the
  instantaneous files. `fieldvars` also accepts the divergence of the
  corrected velocity and its terms (`di,ux,vy,wz`), the pressure correction
  and predicted velocities of the last substep (`pd,pu,pv,pw`) and the
  Poisson right-hand side and its terms (`rs,du,dv,dw`), all computed at
  output time.
- **Pedestrian-level wind speed.** `ws_local`, the time-mean horizontal wind
  speed at `receptor_height` (or at each of `receptor_heights`), and `ws_10`
  at z = 10 m, in `stats_t` and as full xy-planes in `stats_tkslice`.
- **Water-body facets.** Facets with wall types in (-30, -20] are treated as
  open water with a fixed anchor temperature `waterT` (`&ENERGYBALANCE`;
  negative falls back to `flrT`); see [Water bodies](udales-water-bodies.md).
- **Namelist input refactor.** Namelists are read in `in_readnamelists.f90`
  (before the domain decomposition is set up) and sparse `(i, j, k)` input
  files in `in_readsparsefiles.f90`. A round-trip test mode writes the
  namelists back for `tests/integration/input`, which also checks that every
  namelist variable is broadcast and described in
  `docs/schemas/udales_input_schema.json`. New entry points in `bin/`:
  `ud_run` (run or gather a case on any supported machine), `ud_nml2json`,
  `ud_nml2v3`, `ud_testenvironment`.
- **Thermal-comfort post-processing** (`tools/python/udcomf`): one command,
  `tools/ud_run_thermal_comfort.sh <case> [l|f|c] [common|icl|ecmwf] -- --heights ...`
  (or `ud_run_thermal_comfort.py`), takes a finished case from `stats_t` (or
  `stats_tkslice`) and `facEB` to MRT, UTCI, PET and WBGT at the requested
  receptor heights, optionally merging the outputs of several restarted
  segments first. Physical and sensor constants live in
  `comfort_defaults.json`. See `tools/python/examples/thermal_comfort/`.
- `isolarpos` (`&INPS`) selects the solar-position algorithm of the
  preprocessing and post-processing (1 = built-in NREL SPA port, 2 = pvlib).
- Tests reorganised into `tests/integration`, `tests/regression`,
  `tests/system` and `tests/unit`, driven by `tests/test_suites.yml`.
- **Preprocessing on ECMWF HPC2020, and a test suite that survives the
  terminal.** `tools/write_inputs.sh` takes a fourth argument, the system
  (`common`, `icl` or `ecmwf`), and a new start mode `f` that runs the
  preprocessing in the foreground and returns its exit code. With `c ecmwf` it
  submits a Slurm job (`PREPROC_ACCOUNT`, `PREPROC_QOS`, `PREPROC_WALLTIME`,
  `PREPROC_MEM`), next to the existing PBS job for `icl`.
  `tests/system/ud_test_preprocess.sh` accepts `--system ecmwf` (or the bare
  system word): on `common`/`icl` it now detaches from the terminal and writes
  its console output to `logdir/test_inputs.console.log` (`--foreground` to
  opt out); on `ecmwf` it submits the whole suite as one Slurm job. Each case is
  preprocessed in the foreground, so failures are detected from the exit code.
- **Facet-size uniformisation in the Python preprocessing.** With
  `lfacetuniform = .true.` (namelist `&INPS`), `udprep` first rebuilds the STL
  so that facet areas fall within a band around a target area: facets above the
  band are split and coplanar facets of the same type below it are merged, with
  the geometry, the creases, the facet types and the mesh conformity preserved.
  The result is written to `<stem>_uniform.stl` (plus facet types, a face map
  and a JSON sidecar), `stl_file`/`types_path` are pointed at it, and the
  original STL is left untouched. Options `facet_target_edge`,
  `facet_size_band`, `facet_allow_split`, `facet_allow_merge` and
  `facet_count_max`; the default target keeps the facet count about the same.
  The per-patch work runs in up to `nompthreads` worker processes.
  Also available as `udgeom.uniformize_facets` (`workers=`). See
  [Pre-processing](docs/udales-pre-processing.md#facet-size-uniformisation).
- `udprep` now refuses an existing `facets.inp` whose facet count does not
  match the geometry instead of silently reusing it, and `Section.save_param`
  updates a preprocessing-only variable (such as `stl_file`) in the namelist
  block where it already appears instead of appending a new `&INP` block.
- `docs/udales-output-files.md` now documents the SGS flux definitions, their
  sign convention, and the fact that the SGS statistics carry **no** surface or
  wall-model contribution, so the total-stress profile does not close at the
  wall ([issue #353](https://github.com/uDALES/u-dales/issues/353)).
- New experimental test suite `tests/integration/sgs_statistics` (`runmode`
  1006) validating the SGS fluxes against a manufactured solution on a
  stretched vertical grid, against the solver's own diffusion tendencies, and
  against the sign convention, on a new fixture `tests/cases/300`.
