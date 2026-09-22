# Thermal-comfort input preparation

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
qualified. These files are uDALES-side inputs to the future common-model
thermal-index engine; they do not yet calculate MRT, PET, UTCI, or WBGT.
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
