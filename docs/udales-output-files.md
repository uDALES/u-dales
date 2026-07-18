# Output files

uDALES writes its outputs as NetCDF files. Simulations run on more than one CPU write one file per CPU for most output types; files that are already reduced across the whole domain (e.g. the y-, xy- and xy-time-averaged statistics, the TKE budget and the facet energy balance files) are written once, by rank 0. Per-CPU files use the naming pattern `<name>.<px>.<py>.<expnr>.nc`, where `<px>` and `<py>` are the three-digit x- and y-pencil ranks; reduced files use `<name>.<expnr>.nc`. Use `gather_outputs.sh` (or its HPC/ARCHER2 wrappers) to merge the per-CPU files into single files per experiment — see [Post-processing](udales-post-processing.md).

All switches and frequency parameters below live in the `&OUTPUT`, `&TREES`, `&WALLS` and `&ENERGYBALANCE` namelists of `namoptions`; see the [Configuration](udales-namoptions-overview.md) page for their defaults and constraints.

## Output files

| File | Enabled by | Frequency | Contents |
| ---- | ---------- | --------- | -------- |
| **Statistics (time- and/or space-averaged)** | | | |
| `tdump.<px>.<py>.<expnr>.nc` | `ltdump` | `tstatsdump` | Time-averaged 3D statistics per CPU: mean velocity, temperature, moisture, pressure and scalar fields; PSS defect; turbulent momentum/heat/scalar fluxes; variances and TKE; SGS scalar fluxes. |
| `mintdump.<px>.<py>.<expnr>.nc` | `lmintdump` | `tstatsdump` | Reduced time-averaged 3D statistics per CPU: mean u, v, w, temperature and moisture, and pressure only (lighter-weight alternative to `tdump`). |
| `xytdump.<expnr>.nc` | `lxytdump` | `tstatsdump` | x-, y- and time-averaged 1D (height) statistics: mean velocity/temperature/moisture/pressure profiles, turbulent/kinematic/SGS fluxes, temperature and momentum variances, TKE. |
| `xydump.<expnr>.nc` | `lxydump` | `tsample` | x- and y-averaged instantaneous 1D (height) profiles: velocity, temperature, moisture, pressure, turbulent/SGS momentum and heat fluxes, advective fluxes. |
| `ytdump.<expnr>.nc` | `lytdump` | `tstatsdump` | y- and time-averaged 2D (x, height) statistics: mean velocity/temperature/moisture/scalar profiles, turbulent/kinematic/SGS fluxes, variances. Only functional if the x-direction is not parallelised. |
| `ydump.<expnr>.nc` | `lydump` | `tsample` | y-averaged instantaneous 2D (x, height) profiles: velocity, temperature, moisture, scalars, turbulent/SGS momentum and heat fluxes, advective fluxes. Only functional if the x-direction is not parallelised. |
| `tkedump.<expnr>.nc` | `ltkedump` | `tstatsdump` | xy-averaged TKE budget terms vs. height: buoyancy production, total transport, advection, turbulent transport, SGS transport, shear production, viscous and SGS dissipation. Currently marked unsupported in the namoptions documentation — treat as experimental. |
| `treedump.<px>.<py>.<expnr>.nc` | `ltreedump` (`&TREES`) | `tstatsdump` | Time-averaged 3D vegetation/tree terms per CPU: drag in x, y, z; temperature and moisture source/sink terms; scalar source/sink terms; decoupling factor. |
| **Instantaneous fields** | | | |
| `fielddump.<px>.<py>.<expnr>.nc` | `lfielddump` | `tfielddump` | Instantaneous 3D snapshots of the variables listed in `fieldvars` (e.g. velocity components, temperature, moisture, pressure, scalars, wall stresses, heat flux, IBM masks, divergence), per CPU. |
| **Slices** | | | |
| `kslicedump.<px>.<py>.<expnr>.nc` | `lkslicedump` | `tsample` | Instantaneous horizontal (x, y) slice at `k = kslice` of u, v, w, temperature and moisture, per CPU. |
| `islicedump.<px>.<py>.<expnr>.nc` | `lislicedump` | `tsample` | Instantaneous vertical (y, z) slice at `i = islice` of u, v, w, temperature and moisture; written only by the CPU(s) owning that i-index. |
| `jslicedump.<px>.<py>.<expnr>.nc` | `ljslicedump` | `tsample` | Instantaneous vertical (x, z) slice at `j = jslice` of u, v, w, temperature and moisture; written only by the CPU(s) owning that j-index. |
| **Facets** | | | |
| `fac.<expnr>.nc` | `lwritefac` (`&WALLS`) | `dtfac` | Facet momentum and heat-transfer data per facet, averaged over the `dtfac` window: surface shear stresses `tau_x`, `tau_y`, `tau_z`, pressure and pressure fluctuation, heat transfer coefficients. |
| `facT.<expnr>.nc` | `lwriteEBfiles` (requires `lEB`) | `dtEB` | Facet layer temperature (`T`) and temperature gradient (`dTdz`) per facet and facet layer. |
| `facEB.<expnr>.nc` | `lwriteEBfiles` (requires `lEB`) | `dtEB` | Facet surface energy-balance terms per facet: net shortwave, incoming longwave, outgoing longwave, sensible heat flux, latent heat flux, soil water content. |

Statistics that are sampled every `tsample` and only dumped every `tstatsdump` (`tdump`, `mintdump`, `xytdump`, `ytdump`, `treedump`, `tkedump`) report the average over the preceding `tstatsdump` window; the switches whose frequency is `tsample` (`xydump`, `ydump`, the slice dumps) instead write an instantaneous sample every `tsample`. Sampling/averaging only starts once the simulation time passes `tstatstart`.

## Restart files

If `trestart` is smaller than `runtime`, uDALES periodically writes unformatted (non-NetCDF) restart files per CPU: `initd<ntrun>_<px>_<py>.<expnr>` holds the flow fields (`mindist`, `wall`, `u0`, `v0`, `w0`, `pres0`, `thl0`, `e120`, `ekm`, `qt0`, `ql0`, `ql0h`, plus the current time and timestep), and, if `nsv > 0`, a companion `inits<ntrun>_<px>_<py>.<expnr>` holds the scalar fields `sv0`. `ntrun` is an 8-digit run counter embedded in the filename. These are written by `writerestartfiles` in `modsave.f90`, at the interval set by `trestart`, and additionally whenever an `exit_now.<expnr>` file is detected or (for driver simulations) when the inlet file store is exhausted.

To resume from a restart file, set `lwarmstart = .true.` (or `lstratstart = .true.`) and point `startfile` at the `initd...` filename of the desired restart point (the trailing digits of `startfile` must match the experiment's `iexpnr`); the matching `inits...` file, if present, is read automatically. See the `RUN` and `SCALARS` namelists on the [Configuration](udales-namoptions-overview.md) page for the exact parameter definitions.
