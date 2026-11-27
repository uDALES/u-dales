# Output files

uDALES writes its outputs as NetCDF files. Output that keeps the full horizontal extent is written as one file per x-rank, covering the whole y-extent of that rank's column: `<name>.<px>.<expnr>.nc`, where `<px>` is the three-digit x-rank. Output that lives on an xz-plane (`ins_jslice`, `stats_tjslice`) is written per y-rank instead, `<name>.<py>.<expnr>.nc`. Output that is already reduced over the whole domain (`stats_xy`, `stats_xyt`, `ins_probe`, the facet files) is written once, by rank 0, as `<name>.<expnr>.nc`. Use `gather_outputs.sh` (or `ud_run <machine> gather`) to merge the per-rank files into one `<name>.<expnr>.nc` per experiment — see [Post-processing](udales-post-processing.md).

The same quantity has the same variable name in every file; the averaging is encoded in the file name, not in the variable name (`u` is the mean streamwise velocity in `stats_t`, `stats_xyt`, `stats_yt`, ... and the instantaneous one in `ins_field` and the slices). Scalars are `s1`–`s5`, their variances `s1ps1p`, turbulent fluxes `wps1p` (and `ups1p`, `vps1p` in the 3-D statistics), SGS fluxes `s1sgs`.

All switches and frequency parameters below live in the `&OUTPUT`, `&TREES`, `&WALLS` and `&ENERGYBALANCE` namelists of `namoptions`; see the [Configuration](udales-namoptions-overview.md) page for their defaults and constraints. Two schedules exist:

- **Time-averaged statistics** (`stats_*`, `stats_tree`) are sampled every `tsample` from `tstatstart` and one record is written at the end of every `tstatsdump` window; each record is the mean over the preceding `tstatsdump` seconds. With `tstatsgap > 0` the windows are placed on the absolute schedule `tstatstart + n*(tstatsdump + tstatsgap)` and nothing is sampled in the gaps. A warm start never shortens a window: the accumulators restart from zero and the next record is written `tstatsdump` seconds later. `stats_xy` and `stats_y` are the exceptions: they hold an instantaneous spatial average written every `tsample`.
- **Instantaneous output** (`ins_*`) is written at `tinstantstart + n*tinstantdump`; fields, slices and probes are written at the same instants.

## Output files

| File | Enabled by | Frequency | Contents |
| ---- | ---------- | --------- | -------- |
| **Time-averaged statistics** | | | |
| `stats_t.<px>.<expnr>.nc` | `ltdump` | `tstatsdump` | Time-averaged 3-D statistics: `u`, `v`, `w`, `p`, `pabs`; turbulent fluxes `upwp`, `vpwp`, `upvp`; variances `upup`, `vpvp`, `wpwp` and `tke`; SGS stresses `usgs`, `vsgs`, `wsgs`; with `ltempeq`: `thl`, `tha`, `upthlp`, `vpthlp`, `wpthlp`, `thlpthlp`, `thlsgs`; with `lmoist`: `qt`, `ql`, `wpqtp`, `qtpqtp`, `qtsgs`, and `rh` when also `ltempeq`; per scalar n: `sn`, `upsnp`, `vpsnp`, `wpsnp`, `snpsnp`, `snsgs`; `PSS` (photostationary-state defect, `lchem` with 3 scalars). Also the horizontal wind speed at the receptor height(s) (`ws_local`) and at z = 10 m (`ws_10`), masked where an interpolation level is solid. |
| `stats_tislice.<px>.<expnr>.nc`, `stats_tjslice.<py>.<expnr>.nc`, `stats_tkslice.<px>.<expnr>.nc` | `ltislicedump`, `ltjslicedump`, `ltkslicedump` | `tstatsdump` | The `stats_t` variables on the planes `islice(1:nislice)`, `jslice(1:njslice)`, `kslice(1:nkslice)` only. Independent of `ltdump`. `stats_tkslice` also carries the full xy-planes of `ws_local` and `ws_10`. |
| `stats_xyt.<expnr>.nc` | `lxytdump` | `tstatsdump` | x-, y- and time-averaged profiles: `u`, `v`, `w`, `p`, `pabs`; turbulent fluxes `upwp`, `vpwp`, `upvp`; dispersive fluxes `uw`, `vw`, `uv`; `uu`, `vv`, `ww`; variances `upup`, `vpvp`, `wpwp`, `tke`; `usgs`, `vsgs`, `wsgs`; `thl`, `tha`, `wpthlp`, `wthl`, `thlpthlp`, `thlsgs`; `qt`, `ql`, `wpqtp`, `wqt`, `qtpqtp`, `qtsgs`; `rh`. |
| `stats_yt.<px>.<expnr>.nc` | `lytdump` | `tstatsdump` | y- and time-averaged (x, z) statistics: `u`, `v`, `w`, `p`, `pabs`, `upwp`, `uw`, `upup`, `vpvp`, `wpwp`, `usgs`, `wsgs`; `thl`, `tha`, `wpthlp`, `wthl`, `thlpthlp`, `thlsgs`; `qt`, `ql`, `wpqtp`, `wqt`, `qtpqtp`, `qtsgs`; `rh`; per scalar n: `sn`, `wpsnp`, `wsn`, `snpsnp`, `snsgs`. |
| `stats_tree.<px>.<expnr>.nc` | `ltreedump` (`&TREES`) | `tstatsdump` | Time-averaged vegetation terms: drag `tr_u`, `tr_v`, `tr_w`; `tr_thl`; `tr_qt`, `tr_qtR`, `tr_qtA`, `tr_omega`; scalar source/sink terms. |
| `mintdump.<px>.<py>.<expnr>.nc` | `lmintdump` | `tstatsdump` | Legacy reduced time-averaged statistics per pencil: mean `u`, `v`, `w`, `thl`, `qt` and `p` only. |
| **Instantaneous spatial averages** | | | |
| `stats_xy.<expnr>.nc` | `lxydump` | `tsample` | x- and y-averaged instantaneous profiles: the `stats_xyt` variables except the variances. |
| `stats_y.<px>.<expnr>.nc` | `lydump` | `tsample` | y-averaged instantaneous (x, z) fields: `u`, `v`, `w`, `p`, `pabs`, `upwp`, `uw`, `usgs`, `wsgs`; `thl`, `tha`, `wpthlp`, `wthl`, `thlsgs`; `qt`, `ql`, `wpqtp`, `wqt`, `qtsgs`; `rh`; per scalar n: `sn`, `wpsnp`, `wsn`, `snsgs`. |
| **Instantaneous fields, slices and probes** | | | |
| `ins_field.<px>.<expnr>.nc` | `lfielddump` | `tinstantdump` | Instantaneous 3-D snapshots of the variables selected by `fieldvars`: `u`, `v`, `w`, `p`, `pabs`, `thl`, `tha`, `rh`, `qt`, `ql`, `s1`–`s5`; wall stresses `tau_x`, `tau_y`, `tau_z` and heat flux `thl_flux`; IBM masks `mask_u`, `mask_v`, `mask_w`, `mask_c`; the pressure correction and predicted velocities of the last substep; the Poisson right-hand side and its terms; the divergence of the corrected velocity and its terms. |
| `ins_islice.<px>.<expnr>.nc`, `ins_jslice.<py>.<expnr>.nc`, `ins_kslice.<px>.<expnr>.nc` | `lislicedump`, `ljslicedump`, `lkslicedump` | `tinstantdump` | Instantaneous yz-, xz- and xy-slices at `islice(1:nislice)`, `jslice(1:njslice)`, `kslice(1:nkslice)` of the variables selected by `slicevars` (`u`, `v`, `w`, `p`, `pabs`, `thl`, `tha`, `rh`, `qt`, `ql`, `s1`–`s4`). Each file has a slice dimension holding the plane positions; written only by the ranks that own a plane. |
| `ins_probe.<expnr>.nc` | `lprobedump` | `tinstantdump` | Time series of the `probevars` variables at the `nprobe` points of `probe.inp.<expnr>`, with the point coordinates `xt`, `yt`, `zt`. Dimensions `(time, point)`. |
| **Facets** | | | |
| `fac.<expnr>.nc` | `lwritefac` (`&WALLS`) | `dtfac` | Facet momentum and heat-transfer data per facet, averaged over the `dtfac` window: surface shear stresses `tau_x`, `tau_y`, `tau_z`, pressure and pressure fluctuation, heat transfer coefficients. |
| `facT.<expnr>.nc` | `lwriteEBfiles` (requires `lEB`) | `dtEB` | Facet layer temperature (`T`) and temperature gradient (`dTdz`) per facet and facet layer. |
| `facEB.<expnr>.nc` | `lwriteEBfiles` (requires `lEB`) | `dtEB` | Facet surface energy-balance terms per facet: net shortwave, incoming longwave, outgoing longwave, sensible heat flux, latent heat flux, soil water content. |

`tkedump.<expnr>.nc` (`ltkedump`) is **disabled**: the TKE-budget terms are never accumulated, so setting `ltkedump = .true.` aborts the run with an error instead of writing a file of zeros. Tracked in [issue #352](https://github.com/uDALES/u-dales/issues/352).

`tha` (air temperature) requires `ltempeq`, `rh` (relative humidity) requires `ltempeq` and `lmoist`, `pabs` is the hydrostatic absolute pressure profile, `ql` the liquid-water specific humidity. In the instantaneous files they are selected with the labels `ta`, `rh`, `pa`, `ql`.

Files written before this release used different names: `tdump` → `stats_t`, `xytdump` → `stats_xyt`, `xydump` → `stats_xy`, `ytdump` → `stats_yt`, `ydump` → `stats_y`, `treedump` → `stats_tree`, `fielddump` → `ins_field`, `islicedump`/`jslicedump`/`kslicedump` → `ins_islice`/`ins_jslice`/`ins_kslice`, and variable names carried the averaging as a suffix (`ut`, `uxyt`, `upwpt`, `sca1t`, `sv1sgs`, `u_kslice`, ...). See the [changelog](udales-changelog.md) for the full mapping.

For `ws_local` and `ws_10`, u and v are first destaggered to scalar-cell centres and
linearly interpolated in height at each `tsample`. The instantaneous horizontal
speed `sqrt(u^2 + v^2)` is then accumulated, so the result is the time mean of
wind-speed magnitude rather than the magnitude of time-mean u and v. Heights are
relative to model ground (`z=0`); a point is masked if either vertical interpolation
level is solid or the requested height lies outside the scalar-level range. Each
`stats_t` and `stats_tkslice` store `receptor_height` as a scalar coordinate for one
height, preserving the existing `(time, yt, xt)` `ws_local` layout. Set
`nreceptor_heights` and `receptor_heights` to request several heights; then
`ws_local` has dimensions `(time, receptor_height, yt, xt)` and the height
coordinate is one-dimensional. The first list entry must match the scalar
`receptor_height`, which remains the single-height setting used by preprocessing.
`ws_10` is always a `(time, yt, xt)` field fixed at 10 m.

## Subgrid-scale (SGS) flux variables

The statistics files listed above contain SGS counterparts of the resolved
fluxes, named `usgs`, `vsgs`, `wsgs`, `thlsgs`, `qtsgs` and `s1sgs`–`s4sgs`.

### Sign convention

**These are fluxes, and they carry the same sign as the resolved fluxes written
alongside them in the same file with the same units.** Explicitly, with `K_m`
the eddy viscosity (`ekm`) and `K_h` the eddy diffusivity (`ekh`):

| Variable | Definition |
| -------- | ---------- |
| `usgs` | `<u'w'>_sgs = -K_m (du/dz + dw/dx)`, at the uw corner |
| `vsgs` | `<v'w'>_sgs = -K_m (dv/dz + dw/dy)`, at the vw corner |
| `wsgs` | `tau_33 = -2 K_m dw/dz`, at the cell centre |
| `thlsgs` | `<w'thl'>_sgs = -K_h dthl/dz`, at the w-level |
| `qtsgs` | `<w'qt'>_sgs = -K_h dqt/dz`, at the w-level |
| `s1sgs`–`s4sgs` | `<w's_n'>_sgs = -K_h ds_n/dz`, at the w-level |

The total flux is therefore the sum of the resolved and SGS parts, for example
`uw + upwp + usgs` in `stats_xyt` for the total vertical flux of streamwise momentum.

!!! warning "Breaking change"
    Before the fix of issue #306 (PR #358, November 2026) these variables were stored **without** the leading minus
    sign, i.e. as `+K dphi/dn`, and were therefore opposite in sign to the
    resolved fluxes in the same file. Any analysis that compensated for that
    (by negating the SGS term, or by subtracting rather than adding it) must be
    updated. Note also that `wsgs` was previously not a flux at all but the
    vertical part of the `w`-diffusion tendency, in `m/s^2` rather than
    `m^2/s^2`, and located at the w-level rather than the cell centre.

### No surface contribution

The SGS statistics are computed from the interior eddy-diffusivity formulas
only. The fluxes that the wall model and the immersed-boundary method actually
impose at surfaces (`momfluxb`/`tfluxb`, the facet stresses `tau_x`, `tau_y`,
`tau_z`, and `thl_flux`) are **not** accumulated into them. Consequences:

- The total-stress profile does not close at the wall: near solid surfaces the
  resolved plus SGS flux is smaller in magnitude than the true total flux by the
  missing surface contribution.
- Levels that the masks exclude — in particular the bottom level `kb` of the
  w-located profiles, which the IBM masks out entirely — are reported as `-999.`
  rather than as a number.
- In the 3-D `stats_t` output, `s1sgs`–`s4sgs` are likewise set to `-999.` at
  solid and wall-adjacent faces, where the interior formula does not describe
  what the solver applied.

Adding the surface contribution is tracked in
[issue #353](https://github.com/uDALES/u-dales/issues/353).

## Restart files

If `trestart` is smaller than `runtime`, uDALES periodically writes unformatted (non-NetCDF) restart files per CPU: `initd<ntrun>_<px>_<py>.<expnr>` holds the flow fields (`mindist`, `wall`, `u0`, `v0`, `w0`, `pres0`, `thl0`, `e120`, `ekm`, `qt0`, `ql0`, `ql0h`, plus the current time and timestep), and, if `nsv > 0`, a companion `inits<ntrun>_<px>_<py>.<expnr>` holds the scalar fields `sv0`. `ntrun` is an 8-digit run counter embedded in the filename. These are written by `writerestartfiles` in `modsave.f90`, at the interval set by `trestart`, and additionally whenever an `exit_now.<expnr>` file is detected or (for driver simulations) when the inlet file store is exhausted.

A third restart file, `initf<ntrun>.<expnr>`, holds the facet state. It is written whenever `nfcts > 0` and at least one of `lEB`, `iwalltemp = 2`, `iwallmom = 2` or `iwallmoist = 2` is active. Unlike the flow fields, the facet arrays are global rather than domain-decomposed, so this file is written by rank 0 only and carries no `_<px>_<py>` rank fields in its name. It is self-describing: record 1 is an integer header (format version, `nfcts`, `nfaclyrs`), followed by the wall temperature profile `facT` (`nfcts` facets by `nfaclyrs + 1` layer interfaces), the soil moisture `facwsoil` (only when `lEB`), the saturation humidity `facqsat`, the relative humidity `fachurel`, the green-roof resistance factors `facf`, and the current time. When `lEB` is active the file continues with the state of the energy-balance interval in progress: the time of the last and of the next energy-balance update (`tEB`, `tnextEB`) and the facet heat and moisture fluxes accumulated since the last update (`fachfi`, `facefi`). On a warm start these let the interval continue where the checkpoint interrupted it. The file is the only source of the facet state on a warm start: with `lreadfacT = .true.` (the default) `Tfacinit.inp.<expnr>` is not read, and a missing, truncated or differently versioned file aborts the run.

To resume from a restart file, set `lwarmstart = .true.` (or `lstratstart = .true.`) and point `startfile` at the `initd...` filename of the desired restart point (the trailing digits of `startfile` must match the experiment's `iexpnr`); the matching `inits...` files (when `nsv > 0`) and the `initf...` file are read automatically. If `initf...` is missing, or its format version, `nfcts` or `nfaclyrs` do not match the current build and `namoptions`, the run aborts with an explicit message. Set `lreadfacT = .false.` to deliberately re-seed the facet temperatures from `Tfacinit.inp.<expnr>` (or `Tfacinit_layers.inp.<expnr>` with `lfacTlyrs`) on a warm start; that is the only case in which those files are read on a warm start. See the `RUN` and `SCALARS` namelists on the [Configuration](udales-namoptions-overview.md) page for the exact parameter definitions.
