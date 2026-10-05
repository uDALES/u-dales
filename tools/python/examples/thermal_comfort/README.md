# Thermal-comfort workflow examples

These scripts implement the two postprocessing boundaries deliberately kept
separate in `udcomf`:

1. `udales_to_comfort.py` reads saved uDALES results and writes one strict,
   model-neutral exchange file per receptor height.
2. `calculate_indices.py` reads exchange files from uDALES or another model and
   calculates MRT, PET, UTCI, and WBGT with the same protocol.

Neither script runs the LES. The first script reuses saved uDALES statistics,
facet energy-balance output, shortwave forcing, longwave forcing, geometry,
and preprocessing files. Pedestrian radiation is the expensive postprocessing
step because visibility must be evaluated at every requested receptor.

## Install

Use the repository virtual environment and install the optional comfort
dependency:

```bash
tools/python/.venv/bin/python -m pip install -e 'tools/python[comfort]'
```

## Verified georeferencing

The exchange writer does not infer geographic coordinates from the solar
longitude/latitude or from the STL. Prepare an NPZ containing these numeric
arrays:

```python
import numpy as np

np.savez_compressed(
    "paris_georeference.npz",
    x=projected_model_x_centres,          # (nx,)
    y=projected_model_y_centres,          # (ny,)
    x_bounds=projected_model_x_bounds,    # (nx, 2)
    y_bounds=projected_model_y_bounds,    # (ny, 2)
    longitude=longitude_on_les_grid,      # (nx, ny), WGS84 degrees east
    latitude=latitude_on_les_grid,        # (nx, ny), WGS84 degrees north
    z_ground=geodetic_ground_elevation,   # (nx, ny), constant for flat uDALES ground
)
```

Create a JSON file beside it. Every value must describe the actual simulation;
the placeholders below are not defaults:

```json
{
  "georeference_npz": "paris_georeference.npz",
  "crs_attributes": {
    "grid_mapping_name": "<verified CF grid mapping name>",
    "crs_wkt": "<complete WKT describing the projected model axes>"
  },
  "vertical_datum": "<vertical datum of z_ground>",
  "model_version": "<uDALES commit used for the simulation>",
  "institution": "<contributing institution>",
  "building_representation": "<actual IBM and facet representation>",
  "terrain_convention": "<actual flat-ground convention>",
  "spinup_start_utc": "2023-08-20T18:00:00Z",
  "spinup_end_utc": "2023-08-21T00:00:00Z"
}
```

For a rotated LES, the CRS must genuinely describe the projected model axes.
Labelling rotated or offset local coordinates directly as EPSG:2154 is wrong.
The script checks shape, spacing, bounds, coordinate ranges, and constant
ground elevation, but it cannot prove that the supplied transform is correct.

## uDALES outputs to exchange files

The case directory name must be its experiment number and contain gathered
`stats_kslice.<expnr>.nc` or `stats_t.<expnr>.nc`, complete facet output, and
the inputs documented in `tools/python/udcomf/README.md`.

```bash
PY=tools/python/.venv/bin/python
CASE=/home/dipanjan/simulation/udtest/experiments/305

$PY tools/python/examples/thermal_comfort/udales_to_comfort.py \
  --case-dir "$CASE" \
  --analysis-day 2023-08-21 \
  --metadata "$CASE/thermal_comfort_metadata.json" \
  --statistics-source stats_kslice \
  --tile-size 256 \
  --n-mu 8 \
  --n-azimuth 32
```

All `receptor_heights` configured in `namoptions` are processed by default.
Use `--heights 1.1 1.5 2.0` to select a configured subset. Existing final files
are never replaced unless `--overwrite` is supplied. Matching radiation tile
checkpoints are reused by default; `--no-resume` disables that behavior.

For warm-started facet output, repeat `--facet-source` for each source case or
`facEB` file before generating longwave radiation:

```bash
$PY tools/python/examples/thermal_comfort/udales_to_comfort.py \
  --case-dir "$CASE" \
  --analysis-day 2023-08-21 \
  --metadata "$CASE/thermal_comfort_metadata.json" \
  --facet-source "$CASE/facEB.312.nc" \
  --facet-source "$CASE/facEB.313.nc" \
  --facet-source "$CASE/facEB.314.nc"
```

When the pedestrian atmosphere and both radiation files already exist and
have passed their earlier checks, assemble only the exchange files with:

```bash
$PY tools/python/examples/thermal_comfort/udales_to_comfort.py \
  --case-dir "$CASE" \
  --analysis-day 2023-08-21 \
  --metadata "$CASE/thermal_comfort_metadata.json" \
  --skip-atmosphere \
  --skip-radiation
```

For a long headless radiation run:

```bash
LOG="$CASE/udcomf_postprocess.log"
nohup "$PY" -u tools/python/examples/thermal_comfort/udales_to_comfort.py \
  --case-dir "$CASE" \
  --analysis-day 2023-08-21 \
  --metadata "$CASE/thermal_comfort_metadata.json" \
  --statistics-source stats_kslice \
  > "$LOG" 2>&1 < /dev/null &
printf 'PID %s, log %s\n' "$!" "$LOG"
```

The outputs are:

```text
pedestrian_atmosphere.h<height>.<expnr>.nc
pedestrian_shortwave.h<height>.<expnr>.nc
pedestrian_longwave.h<height>.<expnr>.nc
thermal_comfort_inputs_udales.h<height>.<expnr>.nc
```

## Exchange files to comfort indices

Run the same model-neutral command on exchange files from every participating
model. With the documented default person and sensor parameters:

```bash
$PY tools/python/examples/thermal_comfort/calculate_indices.py \
  "$CASE/thermal_comfort_inputs_udales.h1p1.305.nc"
```

Several model or height files can be processed together. Each output is placed
beside its input:

```bash
$PY tools/python/examples/thermal_comfort/calculate_indices.py \
  outputs/thermal_comfort_inputs_udales.h1p1.305.nc \
  outputs/thermal_comfort_inputs_palm.h1p1.nc \
  outputs/thermal_comfort_inputs_urbclim.h1p1.nc
```

The output name is `thermal_comfort_indices_<model>.h<height>.nc`. It contains
hourly `(x,y,time)` fields for `mrt`, `pet`, `utci`, and `wbgt`, plus globe and
natural-wet-bulb temperatures and validity flags.

To use a different common intercomparison protocol, put only changed
`ComfortParameters` fields in one JSON file and pass the same file for every
model. For example:

```json
{
  "human_shortwave_absorptivity": 0.7,
  "human_longwave_emissivity": 0.97,
  "pet_met": 1.37,
  "pet_clo": 0.9
}
```

```bash
$PY tools/python/examples/thermal_comfort/calculate_indices.py \
  --parameters common_comfort_parameters.json \
  thermal_comfort_inputs_*.nc
```

Do not use different person or sensor parameters for different models in one
comparison. Inputs outside an index's validated domain are flagged and written
as missing; the calculation does not clip them into range.

## Scientific scope

The equations, fixed reference cases, and current validation limits are in
`tools/python/udcomf/VALIDATION.md`. The full input contract, restart rules,
quadrature assumptions, and current vegetation limitation are documented in
`tools/python/udcomf/README.md`.
