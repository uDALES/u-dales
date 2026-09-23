# Thermal-comfort postprocessing

Runnable command-line examples for the complete uDALES-to-exchange workflow
and the model-neutral MRT/PET/UTCI/WBGT calculation are in
[`tools/python/examples/thermal_comfort`](../examples/thermal_comfort/README.md).

This document is also the technical handoff for the implementation. It records
the package boundaries, solver and preprocessing changes, physical equations,
file contracts, setup, validation, and unresolved limits. A new developer or
AI agent should read this file and [`VALIDATION.md`](VALIDATION.md) before
changing the pipeline.

## Purpose and status

The objective is to calculate hourly two-dimensional pedestrian-level maps of:

- mean radiant temperature (`mrt`), degree Celsius;
- physiological equivalent temperature (`pet`), degree Celsius;
- Universal Thermal Climate Index (`utci`), degree Celsius;
- wet-bulb globe temperature (`wbgt`), degree Celsius.

The indices are postprocessed after the LES. They are not calculated inside
uDALES. Each map is calculated from atmospheric and radiative quantities
averaged over the preceding 15 minutes at an exact hourly endpoint. Therefore,
the result is an index evaluated from mean inputs, not the time average of an
index evaluated at every LES timestep.

The implementation is complete through the initial validated index engine in
commit `42abbc9a`. Its main development sequence after commit `75835290` was:

| Commit | Main change |
|---|---|
| `34497934` | Opt-in air-temperature output (`ta` selector, NetCDF `tha`). |
| `5527ba52` | Relative-humidity diagnostics (`rh`). |
| `aa8528b5` | Hydrostatic absolute pressure (`pabs`) and liquid water (`ql`). |
| `d193a136` | Initial pedestrian and 10 m mean wind-speed statistics. |
| `84c4ba68` | Configurable scalar `receptor_height`. |
| `f257a4a2` | Generalized scalar-height output as `ws_local`. |
| `5484eef0` | Current statistics and saved-radiation UDBase loaders. |
| `d0847a28` | `udcomf` facade and shared-case package structure. |
| `010407af` | Preservation of original atmospheric shortwave forcing. |
| `43482cee` | Multi-height `ws_local` and `stats_kslice` support. |
| `f01deb53` | Full uDALES-to-exchange input pipeline. |
| `42abbc9a` | MRT, PET, UTCI, WBGT, reference tests, and quadrature correction. |

Use the actual Git history as the source of truth if these commits have since
been rebased or merged.

## Non-negotiable design decisions

1. **No repeated LES or surface-energy-balance physics.** The postprocessor
   reads quantities already saved by preprocessing or the simulation. The only
   new physical calculations before the index stage are receptor visibility,
   directional irradiance, exact-window integration, vertical interpolation of
   atmospheric scalars, and midpoint solar position for exchange metadata.
2. **Reuse existing loaders.** `UDComf` owns the same `UDBase` instance exposed
   as `sim.comf`. It uses `UDBase.load_stat_t`, `load_stat_kslice`,
   `load_fac_eb`, `load_fac_temperature`, `load_timedepsw`, `load_timedeplw`,
   `load_shortwave_forcing`, and existing facet-property assignment. Do not
   create model-specific duplicate NetCDF readers for these files.
3. **Separate model-specific preparation from common physics.** uDALES-specific
   code ends at the strict exchange NetCDF. `calculate_indices` reads that same
   format for uDALES, PALM, UrbClim, or another model.
4. **Keep each model on its native horizontal grid.** There is no common-grid
   remapping or horizontal interpolation. Every final field has dimensions
   `(x,y,time)` on the source model's native two-dimensional plane.
5. **One file per receptor height.** Radiation and atmospheric interpolation
   depend on height. Multiple configured heights are handled sequentially in
   one postprocessing run; the LES is not repeated.
6. **No hidden fallback, clipping, or gap filling.** Invalid units, missing
   fields, incomplete 15-minute windows, negative humidity/irradiance, schema
   disagreement, supersaturation, or values outside an index's validated
   domain are rejected or explicitly flagged and written as missing.
7. **Use one comfort protocol across models.** Human, clothing, activity, and
   sensor parameters must be identical for every participant in an
   intercomparison.

## Data flow and ownership

```text
Python preprocessing
  shortwave_forcing.<expnr>.nc  (DNI, diffuse horizontal, solar angles)
  timedepsw.inp.<expnr>         (absorbed facet shortwave)
  timedeplw.inp.<expnr>         (downward sky longwave)
  geometry/facet/grid inputs
                |
                v
uDALES simulation
  stats_t.<expnr>.nc or stats_kslice.<expnr>.nc
  facEB.<expnr>.nc and facT.<expnr>.nc
                |
                v
uDALES-specific udcomf stage
  pedestrian_atmosphere.h<height>.<expnr>.nc
  pedestrian_shortwave.h<height>.<expnr>.nc
  pedestrian_longwave.h<height>.<expnr>.nc
                |
                v
strict model-neutral exchange stage
  thermal_comfort_inputs_udales.h<height>.<expnr>.nc
                |
                v
model-neutral thermalcomfort stage
  thermal_comfort_indices_<model>.h<height>.nc
```

The exchange boundary is deliberate. A PALM or UrbClim contributor should
produce the same strict exchange variables from that model's output, then run
the identical index routine. Model name or file name does not alter the
physics.

## Package design

| File/module | Responsibility |
|---|---|
| `udbase.py` | Loads the case once and exposes it through `UDBase`; attaches `sim.comf`. |
| `udnetcdf.py` | Shared selective NetCDF read used by UDBase; reverses dimensions to the established Python/MATLAB convention and can select one time and selected `zt` levels before loading. |
| `udcomf/udcomf.py` | Thin facade holding the shared `UDBase` state and entry points for merge and exchange export. |
| `udcomf/heights.py` | Reads and validates scalar or multi-height namelist configuration and creates stable file tags such as `h1p1`. |
| `udcomf/udcomf_io.py` | Streams and merges warm-start `facEB`/`facT` segments, preserving real gaps and rejecting overlaps/schema changes. |
| `udcomf/udcomf_atmosphere.py` | Extracts `ta`, `pabs`, `qv`, `ws_local`, and `ws_10` from gathered `stats_t` or `stats_kslice`. |
| `udcomf/udcomf_radiation.py` | Builds receptor grids, traces first-hit rays, and calculates six shortwave and six longwave receiving-plane irradiances plus local DNI. |
| `udcomf/checkpoints.py` | Exact preceding-window integration, atomic per-hour/per-tile radiation checkpoints, manifest validation, and NetCDF consolidation. |
| `udcomf/udcomf_export.py` | Preflight audit and strict 24-hour, per-height, model-neutral exchange writer. |
| `udcomf/thermalcomfort/physics.py` | MRT radiation equation and globe/wick heat- and mass-balance physics. |
| `udcomf/thermalcomfort/indices.py` | Strict exchange validation, tiled PET/UTCI/WBGT/MRT calculation, flags, metadata, and atomic final output. |
| `udcomf/VALIDATION.md` | Verification matrix, fixed references, tolerances, scientific scope, and remaining physical validation. |
| `examples/thermal_comfort/` | Runnable uDALES-specific and model-neutral command-line workflows. |

`udnetcdf.py` is not a second uDALES case abstraction. It is the low-level
load-and-close helper that `UDBase` itself uses. Selective `time_index` and
`vertical_indices` reads avoid loading a complete multi-terabyte statistics
record when only a few horizontal planes are required.

## Changes to the uDALES solver

The thermal-comfort pipeline required new diagnostics, but existing `thl` and
pressure behavior was preserved. The main files changed were
`src/modthermodynamics.f90`, `src/out_instants.f90`, `src/out_stats.f90`,
`src/out_write_nc.f90`, `src/modglobal.f90`, and
`src/in_readnamelists.f90`.

### Thermodynamic diagnostics

`modthermodynamics` now provides shared `pure elemental` functions so every
output family uses the same pointwise physics:

```text
ta = thl * exner + (Lv/cp) * ql
qv = qt - ql
epsilon = Rd/Rv
e = qv * pabs / (epsilon + (1 - epsilon) * qv)
es = es0 * exp[at * (ta - Tmelt) / (ta - bt)]
rh = 100 * e/es
```

Here `ta` is actual air temperature in K, `thl` is liquid-water potential
temperature in K, `ql` and `qt` are specific humidities in kg/kg, and `pabs` is
the hydrostatic base-state pressure `presf` in Pa. The `p0` output remains the
dynamic kinematic pressure in `m2 s-2`; it is not suitable for humidity or
temperature conversion.

The solver-side `rh` uses the model's existing Tetens saturation expression.
The common comfort engine intentionally does not consume this `rh`: it derives
humidity from `ta`, `qv`, and `pabs` using one common Buck/Liljegren convention
for every model. Solver `rh` remains useful for direct diagnostics and
observation comparisons.

### New output variables

| Quantity | NetCDF name | Units | Availability and meaning |
|---|---|---|---|
| Air temperature | `tha` | K | Statistical output when `ltempeq`; instantaneous selector is `ta`. Calculated at every sample before averaging. |
| Relative humidity | `rh` | `%` | Statistical output when both `ltempeq` and `lmoist`; instantaneous selector is `rh`. Uncapped, so supersaturation remains visible. |
| Hydrostatic absolute pressure | `pabs` | Pa | Base-state `presf`; instantaneous selector is `pa`. No full-domain time accumulator is needed because it is horizontally uniform and time independent in this formulation. |
| Liquid-water specific humidity | `ql` | kg/kg | Output when `lmoist`; instantaneous selector is `ql`. Accumulated directly rather than reconstructed from means. |
| Local horizontal wind speed | `ws_local` | m/s | Time mean of instantaneous horizontal speed at each configured receptor height. Statistics only. |
| UTCI reference wind speed | `ws_10` | m/s | Time mean of instantaneous horizontal speed at 10 m above model ground. Statistics only. |

`th` still selects NetCDF `thl`; adding `ta` did not change the value or normal
output behavior of `thl`. Instantaneous `ta` requires `ltempeq`, `ql` requires
`lmoist`, and `rh` requires both. These checks apply to fields, i/j/k slices,
and probes. `pa` is independent of moisture. The statistics families include
the corresponding variables according to enabled physics.

### Pedestrian wind sampling

At each statistics sample, u and v are first destaggered to scalar-cell
centres. Both components are linearly interpolated in physical height between
bracketing scalar levels, then the instantaneous magnitude is formed:

```text
ws(z) = sqrt(u(z)^2 + v(z)^2)
```

The statistics accumulator therefore stores `mean(sqrt(u^2+v^2))`, not
`sqrt(mean(u)^2+mean(v)^2)`. Vertical velocity is not included. A receptor is
valid only if all scalar levels used for interpolation are fluid. Out-of-range
heights produce fill values and a warning. Restarts reject an existing output
file if its saved receptor-height coordinate differs from the current
namelist.

For one configured height, Python sees `ws_local(x,y,time)` and a scalar
`receptor_height`. For several heights, it sees
`ws_local(x,y,receptor_height,time)`. `ws_10(x,y,time)` remains two-dimensional
per time. These variables are written to both `stats_t` and `stats_kslice`,
including when only k-slice statistics are enabled.

### Receptor-height namelist contract

The `OUTPUT` namelist accepts:

```fortran
receptor_height = 1.1
nreceptor_heights = 3
receptor_heights = 1.1, 1.5, 2.0
```

`receptor_height` is the backward-compatible scalar and defaults to 1.1 m.
When `nreceptor_heights > 0`, the list must contain exactly that many finite,
positive, strictly increasing values; its first value must equal
`receptor_height`. The solver validates on rank 0 and broadcasts the settings.
The maximum supported count is 1000.

`tools/python/udprep/defaults.json` intentionally retains only the scalar
`receptor_height=1.1`. This is not a bug: preprocessing inputs are not
height-specific. `udcomf` reads the full list directly from `namoptions` with
`f90nml`.

### Gathering output

No special variable list is needed in `gather_outputs.sh`. The diagnostics are
ordinary variables inside the existing NetCDF output families, so the normal
gather route carries them with those files. The multi-height coordinate and
`ws_local` dimension were covered by NCO/MPI gathering tests. Always use the
gathered `stats_t.<expnr>.nc` or `stats_kslice.<expnr>.nc`; `udcomf` does not
read rank-split files.

## Changes to the Python suite

1. UDBase statistics loading now uses `stats_t` and `stats_xyt`, replacing the
   obsolete `tdump`/`xytdump` names. `load_stat_kslice` was added for gathered
   horizontal time-average slices.
2. `udnetcdf.load_ncdata` can select a record and selected `zt` levels before
   materializing an array. It preserves the established reversed dimension
   convention used by UDBase and MATLAB-side workflows.
3. Existing facet loaders remain authoritative. `load_fac_eb` gained an
   optional time selector; `load_fac_temperature` remains the `facT` reader.
4. UDBase gained validated readers for `Sdir.nc`,
   `shortwave_forcing.<expnr>.nc`, `timedepsw.inp.<expnr>`, and
   `timedeplw.inp.<expnr>`. Large `timedepsw` text files may be scanned for a
   selected timestamp/facet set without returning the complete matrix.
5. Python shortwave preprocessing and the HARMONIE SSRD converter now write
   `shortwave_forcing.<expnr>.nc` before atmospheric radiation is mapped to
   facets. It stores float64 `time`, DNI, diffuse horizontal irradiance, solar
   zenith, local solar azimuth, source, simulation start, and GHI where
   available. Cached preprocessing can backfill this archive without repeating
   facet tracing after validating cached timestamps.
6. `UDBase` attaches `UDComf(self)` as `sim.comf`; geometry, masks, facet
   properties, forcing, and output state are not loaded a second time.
7. The optional `comfort` package extra pins `pythermalcomfort==4.6.0`; PET and
   UTCI use that external scientific implementation instead of duplicating the
   iterative MEMI model or the operational UTCI sixth-order polynomial.

## Required saved inputs

The uDALES-specific postprocessor expects the following in the case directory
unless an explicit input/output directory is passed:

| Input | Why it is needed |
|---|---|
| `namoptions.<expnr>` and `prof.inp.<expnr>` | Clock, grid, physics flags, output cadence, receptor heights, and scalar-level heights. |
| STL named by the case plus `facets.inp.<expnr>` and `factypes.inp.<expnr>` | First-hit visibility, face ordering/orientation, albedo, and facet radiation mapping. |
| `solid_c.txt` | Excludes receptors inside IBM solids. |
| Gathered `stats_t.<expnr>.nc` or `stats_kslice.<expnr>.nc` | `tha`, `pabs`, `qt`, `ql`, `ws_local`, `ws_10`, coordinates, and saved times. |
| `shortwave_forcing.<expnr>.nc` | Atmospheric DNI, diffuse horizontal irradiance, and solar angles before facet mapping. |
| `timedepsw.inp.<expnr>` | Facet absorbed net shortwave after uDALES preprocessing shading/reflection. |
| `facEB.<expnr>.nc` | Facet timestamps and emitted `LWout`; also carries archived SEB output. |
| `timedeplw.inp.<expnr>` | Downward sky-longwave forcing interpolated to facet-output times. |
| `facT.<expnr>.nc` | Paired and preserved during warm-start merging for schema/time consistency, although receptor longwave uses archived `facEB.LWout`. |

Do not rerun preprocessing or the LES merely to duplicate a quantity already
present in these files. Receptor radiation is calculated after the LES because
it is not one of the saved facet or atmospheric outputs.

## Atmospheric extraction

For each configured height, `UDComfAtmosphere.write_hourly_heights`:

- opens gathered `stats_t` or `stats_kslice` through UDBase;
- reads only the scalar levels bracketing that height;
- linearly interpolates saved `tha`, `pabs`, `qt`, and `ql` in physical height;
- writes `ta=tha`, `qv=qt-ql`, and `pabs`;
- copies the solver-saved `ws_local` plane for the matching height and the
  solver-saved `ws_10` plane; it never derives speed from mean u/v;
- requires fluid cells at the receptor and both interpolation levels;
- preserves the actual statistics timestamp as `source_time` and describes the
  requested `[hour-900 s,hour)` bounds separately.

Negative `qv`, missing bracketing levels, duplicate/missing times, mismatched
height coordinates, or grid disagreement are hard errors.

## Receptor radiation physics

Radiation is evaluated at points `(xt,yt,receptor_height)` over fluid cells.
The height is above uDALES model ground `z=0`; geodetic ground elevation is a
separate exchange coordinate. Horizontal receiving-plane normals are rotated
with the same `azimuth - xazimuth` convention as preprocessing.

### Angular integration and visibility

The sphere uses Gauss-Legendre sampling in vertical cosine and uniform azimuth,
with defaults `n_mu=8`, `n_azimuth=32`. Every ray is assigned to its first
front-facing STL facet, unobstructed sky, or opaque unmeshed model ground.
Projected weights for each receiving plane are normalized so constant radiance
integrates exactly. This removed a 0.1255 K warm MRT bias in the default
uniform-enclosure test; checkpoint manifest version 2 prevents mixing results
from before and after that correction.

Static first-hit directions are traced once per receptor and reused over all
source times. The changing direct-sun ray is tested separately at each source
time. This is directional hemispherical quadrature, not a stored
receptor-to-facet view-factor matrix.

### Shortwave

`timedepsw.netsw` is absorbed facet shortwave irradiance. For facet albedo `a`,
the reflected exitance recovered for a Lambertian source is:

```text
M_reflected = a * netsw/(1-a)
L_reflected = M_reflected/pi
```

Albedo must satisfy `0 <= a < 1`. Archived `dsky` is diffuse horizontal
irradiance and is represented as isotropic upper-sky radiance `dsky/pi`.
Each receptor output contains six non-direct receiving-plane irradiances
(visible sky diffuse plus visible facet reflection) and a separate
`sw_direct_normal` equal to DNI only when the local sun ray is clear. Direct
beam is never folded into the six non-direct fields. Unmeshed downward ground
is opaque and contributes zero reflected shortwave.

### Longwave

`facEB.LWout = emissivity * sigma * T_surface^4` is treated as Lambertian facet
exitance, so facet radiance is `LWout/pi`. `timedeplw.LWsky` is represented as
isotropic sky radiance `LWsky/pi` and linearly interpolated to each facet-output
timestamp. `facEB.LWin` is emissivity-weighted absorbed longwave and is not a
source radiance.

The uDALES facet exchange does not include reflected longwave, so the
postprocessor does not invent it. A longwave ray that reaches unmeshed ground
or a back-facing facet is a hard error because its source emission is unknown;
a complete, consistently oriented ground/facet mesh is required.

### Time integration and checkpoints

Hourly radiation values are exact integrals of a piecewise-linear interpolation
of saved source samples over `[hour-900 s,hour)`. They are not claimed to be
exact physical sub-cadence means. A window is incomplete when its endpoints
are not bracketed or an internal source interval exceeds `max_gap` (default
1.1 times the median source cadence). No extrapolation or gap filling occurs.

Checkpoints are atomic NPZ files per hour and receptor tile. Their manifest
contains the case, height, grid, receptor selection hash, quadrature, target
times, file paths/sizes/modification times, and algorithm version. Resume is
allowed only when the manifest is identical. `tile_size` controls checkpoint
granularity, not parallel execution.

## Strict model-neutral exchange contract

Each exchange file contains exactly 24 records at `01:00,...,24:00`, with time
bounds `[time-900 s,time)`, on one native `(x,y)` grid and one scalar receptor
height. Every physical field has dimensions `(x,y,time)`, `cell_methods =
"time: mean"`, a CF grid mapping, longitude/latitude coordinates, and an
explicit fill value outside one fixed pedestrian mask.

Required atmospheric and radiative fields are:

| Group | Variables | Units |
|---|---|---|
| Atmosphere | `ta`, `pabs`, `qv`, `ws_local`, `ws_10` | K, Pa, kg kg-1, m s-1, m s-1 |
| Shortwave | `sw_direct_normal`; `sw_nondirect_{up,down,north,south,east,west}face` | W m-2 |
| Longwave | `lw_{up,down,north,south,east,west}face` | W m-2 |

It also requires increasing projected `x/y` centres and contiguous cell
bounds, WGS84 longitude/latitude on `(x,y)`, geodetic `z_ground` and its named
vertical datum, a scalar CF CRS with complete WKT, solar zenith and true-north
clockwise azimuth at each window midpoint, model/grid/provenance attributes,
and exact UTC spin-up endpoints.

The uDALES exchange writer requires constant `z_ground` because the current
model ground is flat. It does not infer projection, origin, model rotation, or
vertical datum from solar coordinates or the STL. These must be supplied from
a verified geospatial transform. Assigning EPSG:2154 directly to rotated or
offset local LES axes is not valid.

## Thermal-comfort physics

### Mean radiant temperature

The standing-person six-direction factors are 0.06 for up and down and 0.22
for each vertical side, which sum to one. Human shortwave absorptivity defaults
to 0.70 and longwave emissivity to 0.97. The orientation-averaged direct-beam
factor above the horizon is:

```text
f_direct = 0.28*sin(zenith) + 0.06*cos(zenith)
```

The body-weighted non-direct shortwave, separate direct beam, and body-weighted
longwave are converted to MRT through:

```text
K_body = sum(f_i * K_nondirect_i) + f_direct * DNI
L_body = sum(f_i * L_i)
Tmrt = [(alpha_sw*K_body + epsilon_p*L_body)/(epsilon_p*sigma)]^(1/4)
```

MRT is evaluated internally in K and written in degree Celsius. Human azimuth
is averaged, so saved solar azimuth is validated but does not alter this
standing-person factor. This is a six-plane approximation, not angularly
resolved human geometry.

### Humidity shared by PET, UTCI, and WBGT

The common engine derives vapour pressure from `qv` and `pabs` and relative
humidity from the Buck liquid-water saturation equation in
`thermalcomfort/physics.py`. Supersaturation is not capped to 100 percent; it
is flagged invalid. This guarantees the same humidity convention for every
model supplying the exchange format.

```text
epsilon = 18.015/28.97
e = qv*pabs/[epsilon + (1-epsilon)*qv]
es = 100*1.004*6.1121*exp[17.502*(ta-273.15)/(ta-32.18)]
RH = 100*e/es
```

### PET

PET uses `pythermalcomfort.models.pet_steady`, the steady MEMI heat-balance
implementation. It receives air temperature, MRT, local receptor wind,
relative humidity, and actual pressure in hPa. The default standing reference
person is 35 years old, male, 1.75 m, 75 kg, 1.37 met (79.7 W m-2), and 0.9 clo.
The `pet_met` argument is in met; 1.37 met is equivalent to 79.7 W m-2 of body
surface area, not 80 W total-body power. A vector calculation is retried
cell-by-cell only to
isolate numerical failures, which are flagged rather than filled with a
fallback value.

### UTCI

UTCI uses `pythermalcomfort.models.utci` with air temperature, MRT, relative
humidity, and the saved 10 m horizontal wind. Output rounding is disabled.
Cells are calculated only inside the operational approximation domain:

```text
-50 <= Ta <= 50 degree C
-30 <= Tmrt-Ta <= 70 K
0.5 <= ws_10 <= 17 m/s
vapour pressure <= 5000 Pa
```

The wind is not clipped to 0.5 or 17 m/s.

### WBGT

WBGT solves local Liljegren-style black-globe and naturally ventilated wetted-
wick balances, replacing the reference program's flat-ground radiation
reconstruction with the saved six-direction urban radiation:

- globe: 50.8 mm diameter, emissivity/absorptivity 0.95, six-plane mean
  non-direct irradiance, direct DNI over one quarter of sphere area;
- wick: 7 mm diameter by 25.4 mm length, longwave emissivity 0.95, shortwave
  absorptivity 0.60, lateral/end-face area weighting, and direct side
  projection based on `sin(zenith)/pi` plus the end face;
- convection, viscosity, conductivity, vapour diffusion, latent heat, and
  saturation pressure follow the implemented Liljegren correlations;
- fixed-point iteration uses up to 80 iterations, 0.02 K convergence tolerance,
  and under-relaxation for unconverged cells.

The forced-flow domain requires `ws_local >= 0.13 m/s`; the implemented wick
property range requires `283.15 <= Ta <= 313.15 K`. No wind floor is imposed.
ISO 7243 weighting is selected per cell: with local direct sun,
`0.7*Tnwb + 0.2*Tg + 0.1*Ta`; without direct sun,
`0.7*Tnwb + 0.3*Tg`. The solved globe and natural wet-bulb temperatures are
retained in the output for audit.

This is an urban directional adaptation of Liljegren, not a claim of bitwise
equivalence to or instrument validation of the original flat-ground program.

### Validity and output

The final file preserves coordinates, masks, time bounds, CRS, and provenance.
It writes `mrt`, `pet`, `utci`, `wbgt`, `globe_temperature`, and
`natural_wet_bulb_temperature` in degree Celsius. The `validity_flags` bits are:

| Bit | Meaning |
|---:|---|
| 1 | Invalid relative humidity or vapour pressure. |
| 2 | UTCI outside its operational domain or non-finite. |
| 4 | WBGT outside its implemented sensor-property domain. |
| 8 | WBGT globe or wick balance did not converge. |
| 16 | PET solver failed. |

Global `validity_counts` records counts for all flags. MRT remains available
when it is physically finite even if another index is invalid.

## Virtual environment setup

Use the one canonical repository environment, `tools/python/.venv`. Do not add
notebook `sys.path` modifications or maintain a second local environment. From
the repository root on Linux/WSL:

```bash
bash tools/python/setup_venv.sh common preprocessing_tools
tools/python/.venv/bin/python -m pip install -e 'tools/python[comfort]'
tools/python/.venv/bin/python -m pip check
```

On the ICL system use:

```bash
bash tools/python/setup_venv.sh icl preprocessing_tools
tools/python/.venv/bin/python -m pip install -e 'tools/python[comfort]'
```

The setup script installs the uDALES tools editable and builds mandatory
View3D, direct-shortwave f2py, and IBM f2py preprocessing artifacts. If the
chosen Python lacks `Python.h`, set `PYTHON_BIN` to an interpreter with
development headers; `/opt/pbs/python/bin/python3` has previously been a valid
cluster choice.

Core tools declare Python 3.9 or newer, but `pythermalcomfort` 4.x requires
Python 3.10 or newer. The comfort extra pins `pythermalcomfort==4.6.0`, whose
metadata requires `numpy>=1.21,<2.3`. The environment validated during this
work used Python 3.13.9 and NumPy 2.2.6. Do not force NumPy 2.3 or newer while
using this pythermalcomfort version.

For a fresh clone on another machine, initialize the repository submodules
before setup:

```bash
git submodule update --init --recursive
```

Case directories, gathered LES outputs, forcing files, radiation checkpoints,
and generated exchange/index NetCDF files are runtime data and are normally not
stored in Git. Transfer the required case files listed above separately. Also
verify that the checked-out branch contains the implementation commits listed
under "Purpose and status" plus the runnable example scripts before relying on
this workflow.

## Solver configuration before an expensive run

For the agreed hourly protocol, configure the LES before running it. For a
simulation starting at 18:00 UTC on the preceding day and an analysis day
starting at model time 21600 s, the essential settings are:

```fortran
! In the appropriate existing namelist groups:
ltempeq = .true.
lmoist = .true.
leb = .true.
lwriteebfiles = .true.
ltimedepsw = .true.
ltimedeplw = .true.

! OUTPUT settings:
receptor_height = 1.1
! For several heights only:
nreceptor_heights = 3
receptor_heights = 1.1, 1.5, 2.0

tstatstart = 21600.
tstatsdump = 900.
tstatsgap = 2700.
tsample = 5.

! Choose full statistics or the disk-saving horizontal slices:
ltdump = .false.
ltkslicedump = .true.
nkslice = <number of selected scalar levels>
kslice = <one-based levels bracketing every receptor height>
```

The case must run through the final requested endpoint. `tsample` must be
positive and no longer than 900 s; preflight cautions above 5 s because rapidly
varying conditions may be undersampled. Every receptor height and 10 m must lie
within the scalar-centre height range. For k-slice output, include both scalar
levels bracketing every receptor height (only one if the height exactly equals
a scalar level). `ws_local` and `ws_10` themselves are saved directly and do
not require reconstruction from the selected k slices.

`stats_kslice` greatly reduces output volume, but the current solver still
allocates full three-dimensional statistics accumulators. It is a disk-I/O
optimization, not an equivalent memory optimization.

Run the audit before submitting the LES:

```python
from datetime import date
from pathlib import Path
from udcomf.udcomf_export import check_comfort_preflight

report = check_comfort_preflight(Path("/path/to/case"), date(2023, 8, 21))
print(report.blockers)
print(report.cautions)
```

Do not submit while `report.blockers` is nonempty. Preflight cannot prove that
a warm start did not interrupt one of the final 15 minutes before a statistics
dump, so restart times still require a separate audit.

## End-to-end usage

The case directory name must be its experiment number because the current
preflight and merge interfaces derive `namoptions.<expnr>` from that name.

1. Run Python preprocessing and confirm that
   `shortwave_forcing.<expnr>.nc`, `timedepsw.inp.<expnr>`, and
   `timedeplw.inp.<expnr>` exist and are populated.
2. Run uDALES with the output settings above, then gather `stats_t` or
   `stats_kslice` and all required facet output.
3. If the LES used warm starts under several case numbers, merge paired
   `facEB`/`facT` files with `sim.comf.prepare_facet_outputs(...)` or the CLI's
   repeated `--facet-source` option. Never fill restart gaps.
4. Prepare verified georeferencing arrays and metadata as documented in the
   runnable example README. These facts are intentionally not guessed.
5. Generate atmospheric planes, both radiation sets, and exchange files:

```bash
PY=tools/python/.venv/bin/python
CASE=/path/to/experiments/305

$PY tools/python/examples/thermal_comfort/udales_to_comfort.py \
  --case-dir "$CASE" \
  --analysis-day 2023-08-21 \
  --metadata "$CASE/thermal_comfort_metadata.json" \
  --statistics-source stats_kslice \
  --tile-size 256 \
  --n-mu 8 \
  --n-azimuth 32
```

6. Apply the common index implementation to one or several exchange files:

```bash
$PY tools/python/examples/thermal_comfort/calculate_indices.py \
  "$CASE/thermal_comfort_inputs_udales.h1p1.305.nc"
```

Use `--skip-atmosphere` or `--skip-radiation` only when the corresponding final
plane files already exist and match the requested case, heights, hours, and
settings. Existing outputs are not replaced without explicit `--overwrite`.
Radiation checkpoint resume is enabled by default. The full command set,
metadata JSON/NPZ format, multi-model batch calculation, custom common-person
parameters, warm-start example, and `nohup` command are in the runnable example
README linked at the top of this document.

## Performance and operational limits

- Receptor radiation is currently a single-process Python/PyVista ray sweep.
  `tile_size` bounds checkpoint work and consolidation memory; it is not a
  worker count. No GPU implementation exists.
- Heights are processed sequentially to limit memory, but visibility must be
  recalculated because it changes with height.
- The full facet source series is loaded once per radiation kind. Peak memory
  includes geometry plus roughly `nfcts * nsource_times` facet values.
- A complete Paris full-plane run has not yet been performance-qualified.
  Test a receptor subset and quadrature convergence before reserving production
  resources.
- Vegetation attenuation and vegetation longwave emission are not implemented.
  Cases with loaded trees/vegetation raise instead of silently ignoring them.
- Diffuse sky is isotropic; facet reflection/emission is Lambertian; longwave
  reflection is absent because the uDALES facet model does not provide it.
- The radiation output contains irradiance in W m-2. It is never multiplied by
  receptor or facet area during the comfort calculation.
- A complete, correctly ordered and consistently oriented triangular STL is
  essential. Face count and normals are checked against `facets.inp`.
- Warm-start merge preserves gaps. An hourly window crossing a gap remains
  incomplete, and strict exchange export rejects incomplete radiation.

## Verification and remaining scientific evidence

Run the focused and complete checks from the repository root:

```bash
MPLCONFIGDIR=/tmp/matplotlib tools/python/.venv/bin/python -m unittest -v \
  tools/python/tests/test_udcomf_examples.py \
  tools/python/tests/test_udcomf_end_to_end.py

MPLCONFIGDIR=/tmp/matplotlib tools/python/.venv/bin/python \
  tests/run_tests.py python-library
```

At the handoff represented by this document, the complete stream passed 516
Python unit tests with 11 expected skips, followed by all configured
direct-shortwave, preprocessing, UDBase/MATLAB-reference, and View3D integration
groups. The tests cover multi-height wind, masks, selective statistics reads,
warm starts, exact windows, radiation limits, strict schema, model-label
invariance, fixed MRT/PET/UTCI/WBGT references, and a synthetic full chain.

This is software verification, not complete field validation. In particular,
the six-direction urban WBGT adaptation still needs comparison with measured
globe and natural-wet-bulb temperatures. See [`VALIDATION.md`](VALIDATION.md)
for numerical reference values, citations, tolerances, and the exact boundary
between verified equations and pending measurement evidence.

## Detailed component notes

The sections below retain lower-level API behavior, Paris warm-start findings,
file names, and restart details that are useful during real case processing.

## Warm-start facet outputs

When a simulation was restarted under several case numbers, merge its paired
`facEB` and `facT` outputs before loading longwave states. Pass case directories
to `sim.comf.prepare_facet_outputs(case_paths=[...])`. If the files were copied
into one directory, pass the `facEB.<case>.nc` paths instead; matching `facT`
files are found beside them. The same operation is available without loading
the full `UDBase` case:

```python
from pathlib import Path
from udcomf.udcomf_io import merge_facet_outputs

case = Path("/home/dipanjan/simulation/udtest/experiments/305")
sources = [case / f"facEB.{number}.nc" for number in range(312, 318)]
result = merge_facet_outputs(case, sources)
```

This writes `facEB.305.nc` and `facT.305.nc` without changing the sources.
The merge sorts segments by saved simulation time, rejects overlaps and schema
mismatches, and preserves restart gaps exactly; it does not interpolate or
claim that an averaging window crossing a gap is complete. Existing merged
files are not replaced unless `overwrite=True` is explicit.

## Shortwave

`sim.comf.radiation.load_shortwave_state(time_index)` reads the existing
`shortwave_forcing.<expnr>.nc` and `timedepsw.inp.<expnr>` files through
`UDBase`. It checks the selected timestamps and maps the existing facet albedo.
No solar forcing, facet direct beam, or facet-to-facet reflection is rerun.

For each opaque facet, uDALES `netsw` is absorbed shortwave irradiance. With
albedo `a < 1`, reflected exitance is `a * netsw / (1 - a)` in W m-2. The
receptor calculation assumes that this reflection is Lambertian, so its
radiance is reflected exitance divided by pi. An albedo of exactly one cannot
be inverted from absorbed flux and is rejected.

The HARMONIE-derived `dsky` is diffuse **horizontal** irradiance. This stage
assumes an isotropic upper sky with radiance `dsky / pi`. A ray from the receptor
sees the first front-facing facet, unobstructed sky, or opaque unmeshed ground.
The six `sw_nondirect_*face` fields integrate sky diffuse and reflected facet
radiance over each receiving hemisphere. `sw_direct_normal` is DNI when the
sun ray is clear and zero when shaded or below the horizon. Direct beam is
never added to a non-direct field. Outputs are irradiance in W m-2, not power
multiplied by facet area.

`shortwave_at_receptor` calculates one point; `shortwave_plane` returns native
(x, y) arrays with NaN at solid or unselected cells. Hourly means and files
are written by `write_hourly`. Angular quadrature is
configurable via `n_mu` and `n_azimuth`. Results should be checked for
quadrature convergence in obstructed scenes.

## Several receptor heights

The solver keeps `receptor_height` as the legacy scalar and uses
`nreceptor_heights` plus `receptor_heights` for multi-height `ws_local` output.
`udcomf` reads the list with `f90nml` rather than relying on UDBase's scalar
namelist parser. The first list entry must equal `receptor_height`. The
`udprep/defaults.json` scalar stays unchanged: preprocessing inputs are not
receptor-specific.

After gathering `stats_kslice` (or `stats_t`) and merging warm-start facet
outputs, process all saved heights without changing the case configuration:

```python
import numpy as np
from udbase import UDBase

sim = UDBase(305, "/path/to/305", load_geometry=True)
hours = np.arange(25200.0, 108001.0, 3600.0)
atmosphere = sim.comf.atmosphere.write_hourly_heights(
    hours, source="stats_kslice"
)
radiation = sim.comf.radiation.write_hourly_heights(
    hours, tile_size=256, n_mu=8, n_azimuth=32, resume=True
)
```

Each height has its own `pedestrian_atmosphere.h<height>.<expnr>.nc`,
`pedestrian_shortwave.h<height>.<expnr>.nc`, and
`pedestrian_longwave.h<height>.<expnr>.nc`. Fields remain `(x,y,time)`;
`receptor_height` is a scalar coordinate in each file. Radiation checkpoints
live in separate `<kind>/h<height>/` directories. The old `write_hourly`
single-height call retains its existing filenames and checkpoint path.

The atmospheric exporter reads only the saved scalar levels bracketing the
requested heights. It writes `ta` from saved `tha`, `pabs`, `qv` from saved
`qt-ql`, and the saved `ws_local` and `ws_10` means. It never calculates wind
speed from time-mean velocity components. Gathered statistics must contain
each requested hourly record, the complete `ws_local` height coordinate, and
all bracketing scalar levels; missing data are rejected. Valid cells must be
fluid at the receptor and both interpolation levels. The output time axis is
seconds since simulation start, with the actual source record time retained.
The 900 s bounds describe the requested averaging interval; whether a warm
start interrupted that interval must still be checked separately.

Radiation is computed afresh at every height because visibility changes with
height. The STL mesh and existing forcing/facet files are reused; neither
preprocessing nor the LES is rerun. Heights are processed sequentially to
limit peak memory. A Paris full-plane ray sweep has not yet been performance
qualified. These files are uDALES-side inputs to the common-model thermal-index
engine described below; this radiation stage does not itself calculate MRT,
PET, UTCI, or WBGT.
If a multi-height radiation run stops after finishing some heights, rerun with
`resume=True, overwrite=True`: matching completed tile checkpoints are reused
and existing final files are reconsolidated. Changed input files or settings
are rejected by the checkpoint manifest.

## Model-neutral exchange files

After all three per-height files exist, `sim.comf.export_exchange(day, metadata)`
assembles `thermal_comfort_inputs_udales.h<height>.<expnr>.nc` for each configured
height. It streams the saved `(x,y,time)` fields without re-running the LES,
surface energy balance, or receptor radiation. It requires the exact 24 hourly
windows on the analysis day, complete shortwave and longwave coverage of every
atmospheric-valid receptor, matching native grids/heights/times, finite physical
fields, and actual statistics records at the requested hours. Incomplete or
subset radiation files are rejected; outside the fixed pedestrian mask every
field is written with `_FillValue`. Negative humidity is rejected, never clipped.
The only new physical evaluation is solar geometry at each exact 15-minute
window midpoint using the existing `udprep.solar` SPA routine, with the
configured solar longitude, latitude, and elevation. The output azimuth is
clockwise from true north in radians.

Construct an `ExchangeMetadata` object from **verified case georeferencing**:
projected model-axis `x/y` centres and `(n,2)` bounds, WGS84 `longitude` and
`latitude` on the `(x,y)` grid, constant geodetic `z_ground` and its named
vertical datum, and a CF grid-mapping attribute dictionary including
`grid_mapping_name` plus complete `crs_wkt` or `spatial_ref`. Also provide the
actual solver version/commit, institution, terrain and building conventions,
and exact UTC spin-up endpoints. A rotated model grid needs a CRS that genuinely
describes its projected model axes; assigning EPSG:2154 directly to rotated or
offset local `x/y` is incorrect. The writer checks bounds against LES spacing,
but it cannot independently verify the geographic transform you supply.

```python
from datetime import date
from udcomf.udcomf_export import ExchangeMetadata

metadata = ExchangeMetadata(
    x=projected_x, y=projected_y,
    x_bounds=projected_x_bounds, y_bounds=projected_y_bounds,
    longitude=longitude_grid, latitude=latitude_grid,
    z_ground=constant_ground_elevation,
    crs_attributes=verified_cf_crs,
    vertical_datum="<verified vertical datum>",
    model_version="<solver commit used for this run>",
    institution="<contributing institution>",
    building_representation="<actual IBM/facet representation>",
    terrain_convention="<actual flat-ground convention>",
    spinup_start_utc="2023-08-20T18:00:00Z",
    spinup_end_utc="2023-08-21T00:00:00Z",
)
files = sim.comf.export_exchange(date(2023, 8, 21), metadata)
```

These geographic and provenance facts are not inferred from the solar site
coordinates or from the STL. For the Paris case, verify the STL's projected
origin/rotation and geodetic ground reference before constructing metadata.
The assembler cannot prove that a warm start did not interrupt a 15-minute
statistics window; audit restart boundaries separately. It also does not
compute thermal-comfort indices, which belong to the next pipeline step.

## Model-neutral MRT, PET, UTCI and WBGT

Install the optional scientific implementation on Python 3.10 or newer:

```bash
tools/python/.venv/bin/python -m pip install -e 'tools/python[comfort]'
```

The same function reads a strict exchange file from **any** model; it never
opens uDALES case inputs or horizontally remaps another model's grid:

```python
from pathlib import Path
from udcomf.thermalcomfort import calculate_indices, ComfortParameters

output = calculate_indices(
    Path("thermal_comfort_inputs_udales.h1p1.305.nc"),
    parameters=ComfortParameters(),
)
```

The default output is `thermal_comfort_indices_<model>.h<height>.nc`, with
native coordinates, mask, time bounds, and hourly `(x,y,time)` planes for
`mrt`, `pet`, `utci` and `wbgt` in degrees Celsius. `globe_temperature` and
`natural_wet_bulb_temperature` are retained for WBGT auditing. All values are
indices **of the preceding 15-minute mean inputs**, not time means of the
instantaneous indices. The input file must contain every required variable,
units, grid, fixed mask, 24 time bounds, solar angles, georeference and
provenance; missing or negative source values inside the mask are rejected.
Output outside the mask, and indices outside their physical or validated
domains, use `_FillValue`. No meteorological input is clipped. A per-cell
`validity_flags` bit mask and global `validity_counts` record humidity, UTCI,
WBGT sensor-domain/convergence and PET failures.

MRT uses the documented six-direction standing-person weights (0.06 up/down,
0.22 each side), human shortwave absorptivity 0.70 and emissivity 0.97. The
direct beam is separate: its orientation-averaged standing-body factor varies
with solar zenith as `0.28 sin(zenith) + 0.06 cos(zenith)` above the horizon.
Azimuth is checked but not used because person orientation is averaged. This
is a stated six-plane approximation, not angularly resolved human geometry.
The factors and PET person can be changed through `ComfortParameters` but
**must be kept identical across models** in one intercomparison.

PET is `pythermalcomfort`'s steady MEMI implementation for a standing
35-year-old man, 1.75 m, 75 kg, 1.37 met (79.7 W m-2) and 0.9 clo by default.
The `pet_met` input unit is **met**; the W m-2 value is its body-surface-area
equivalent, not 80 W total-body power. UTCI uses that
library's unrounded operational calculation with `ws_10`; outside its
temperature, radiant-temperature, wind or vapour-pressure domain it is missing
and flagged, not extrapolated. PET uses `ws_local` and actual pressure.
Both receive RH derived from `qv` and `pabs`; supersaturation is flagged, not
silently capped at 100%.

WBGT solves a Liljegren-type black-globe and natural-wet-bulb heat/mass balance
with the saved **local** six-plane radiation, air temperature, pressure,
humidity and `ws_local`. A 50.8 mm globe receives direct DNI over one-quarter
of its area and the six-plane mean non-direct irradiance. A 7 mm by 25.4 mm
vertical wick receives lateral and end-face irradiance by surface-area
weighting; its direct side projection is `sin(zenith)/pi`. This replaces the
published Liljegren program's flat-ground shortwave/sky/surface estimate with
an urban directional estimate: it is **not** claimed to be the unmodified,
instrument-validated Liljegren model. Cells below its 0.13 m s-1 forced-flow
domain, outside the wick property fit's 283.15--313.15 K temperature range, or
with nonconvergent sensor balances are missing and flagged; wind is not raised
to the reference code's floor. The ISO 7243 solar-load expression
is selected where local direct beam is positive, and the no-direct-solar
expression elsewhere. This per-cell sun/shade rule is part of the protocol.

References: [SOLWEIG standing-person factors](https://umep-dev.github.io/solweig/physics/tmrt/),
[pythermalcomfort PET/UTCI](https://pythermalcomfort.readthedocs.io/en/stable/documentation/models.html),
[Liljegren reference source](https://github.com/mdljts/wbgt/blob/master/src/wbgt.c),
[ISO 7243](https://www.iso.org/standard/67188.html).
The numerical verification cases and their limits are recorded in
[`VALIDATION.md`](VALIDATION.md). The directional WBGT sensor balances have an
independent root-solver check, but the urban radiation adaptation still needs
comparison with physical globe and wick measurements before it can be called
measurement-validated.

`trace_shortwave_rays` returns a receptor-specific first-hit map. Pass it back
as `ray_map` to `shortwave_at_receptor` at later times to reuse static geometry;
the changing direct-sun ray is still tested for each timestamp.

Current limits: the sky distribution is isotropic, facet reflection is
Lambertian, unmeshed ground is non-emitting, and vegetation attenuation is not
implemented (vegetation cases raise). A full Paris-plane ray sweep has not
been performance-qualified; process small receptor batches before scaling up.

## Longwave

`sim.comf.radiation.load_longwave_state(time_index)` reads one `facEB` output
record and the saved `timedeplw.inp.<expnr>` forcing through `UDBase`. The
forcing is linearly interpolated to the `facEB` timestamp, matching uDALES's
time-dependent sky-longwave interpolation within the forcing range. Outside
that range, the loader rejects the timestamp instead of guessing a value.

The facet source is archived `LWout = emissivity * sigma * T_surface**4` in
W m-2. It already incorporates the facet surface temperature and emissivity;
the comfort calculation does not rerun the surface energy balance. In contrast,
`facEB.LWin` is emissivity-weighted **absorbed** longwave at the facet, not a
source radiance or pedestrian irradiance. uDALES's facet longwave exchange
includes sky and facet emission, but not longwave reflection; this calculation
does not add reflected longwave that the simulation did not represent.

`longwave_at_receptor` integrates isotropic sky radiance `LWsky / pi` and
Lambertian facet radiance `LWout / pi` over six receiving hemispheres. A ray
uses the first visible, front-facing facet or the unobstructed sky. The six
`lw_*face` outputs are irradiance in W m-2, never facet-area-weighted fluxes.
The same `trace_shortwave_rays` first-hit map can be passed as `ray_map` for
multiple longwave times. `longwave_plane` returns native `(x, y)` arrays with
NaN at solid or unselected cells. Hourly means and files are written by
`write_hourly`.

An unmeshed-ground or back-face hit is rejected: its emitted longwave is
unknown, so zero would silently bias the comfort calculation. A complete,
consistently oriented ground facet mesh is required. Vegetation emission and
attenuation are not implemented (vegetation cases raise). The sky is isotropic,
facet emission Lambertian, and the quadrature should be checked for convergence
in obstructed scenes. No full Paris-plane performance claim is made here.

## Hourly output and restart

The write_hourly method takes exact-hour endpoints in simulation seconds and
writes means over the preceding 900 seconds. Each saved source timestamp is a
point sample. The integral is exact under piecewise-linear interpolation
between those samples; it is not an exact physical sub-cadence mean. The
default maximum accepted source interval is 1.1 times the median cadence;
set max_gap explicitly when that assumption is unsuitable. Windows outside
the source range or crossing a longer interval are incomplete and contain
NaN, never gap-filled.

Both files use native (x,y,time) fields in W m-2, time_bounds,
time_complete, valid_mask, processed_mask, coordinates, receptor_height,
and provenance. processed_mask distinguishes a subset check from a full
plane. Unselected and invalid cells are NaN.

For example, use this Python code to request the same hourly endpoints for
both radiation kinds:

    import numpy as np
    from udbase import UDBase

    sim = UDBase(305, "/path/to/305", load_geometry=True)
    hours = np.arange(21600.0, 108001.0, 3600.0)
    for kind in ("shortwave", "longwave"):
        sim.comf.radiation.write_hourly(
            kind, hours, tile_size=256, n_mu=8, n_azimuth=32,
            resume=True,
        )

By default, the files are pedestrian_shortwave.<expnr>.nc and
pedestrian_longwave.<expnr>.nc in the case directory, with separate
udcomf_radiation.checkpoints/<kind>/ directories. Pass flat_indices for a
small receptor subset, or output_path and checkpoint_dir to write elsewhere.
Checkpoint identity includes input path, size, modification time, selected
receptors, quadrature, and target times. Resume skips finished hour/tile
checkpoints; changed inputs or settings require a new checkpoint directory.
Use overwrite=True to replace an existing final file.

The full facet source series is loaded once through UDBase for each kind.
Static first-hit rays are traced once per receptor tile and reused across
needed source times; direct-sun rays remain time-dependent. Peak memory still
includes the full geometry and one full facet-source matrix. A Paris full-plane
sweep has not been timed or qualified. Start with a subset and check angular
convergence before allocating production resources.

For the merged 305 warm starts, among hours 6 through 30 inclusive, all
shortwave windows have source coverage. Longwave hours 8, 25, and 29 cross
restart gaps; hour 30 is also incomplete because the final saved facet
record is at 107999.016 s rather than 108000 s. These hours remain present
with time_complete=0 and NaN longwave fields.

## Before the LES run

Run the configuration check before submitting a Paris simulation:

    from datetime import date
    from pathlib import Path
    from udcomf.udcomf_export import check_comfort_preflight

    report = check_comfort_preflight(Path("/path/to/case"), date(2023, 8, 21))
    print(report.blockers)
    print(report.cautions)

For a run starting at 2023-08-20
18:00 UTC, the 24 required comparison planes end at model times 25200,
28800, ..., 108000 s. The needed statistics settings are tstatstart=21600,
tstatsdump=900, and tstatsgap=2700, with ltempeq and lmoist enabled. Enable
either ltdump or ltkslicedump; for k-slice-only output, configure kslice to
include both scalar levels bracketing every receptor height. The case must run
through the last dump. Every receptor height must lie between scalar-centre
heights, as must 10 m for ws_10. Preflight checks these requirements without
loading the large geometry or running the solver.

The check does not prove that a warm-started averaging window is complete:
restarts inside its last 15 minutes can yield a partial mean at the correct
timestamp. Nor does it validate surface radiation or geospatial metadata.
Check actual restart boundaries and output records after the run.

Full stats_t includes at least 28 three-dimensional float32 fields at each
dump when temperature and moisture are enabled. For 24 dumps, the raw data
payload is at least about 1.03 TB for the 305 grid and 3.72 TB for the 350
grid, before slices, probes, restarts and other outputs. stats_kslice is the
disk-saving alternative, but the solver still allocates full 3D time-average
accumulators; k-slice-only does not deliver an equivalent memory reduction.
