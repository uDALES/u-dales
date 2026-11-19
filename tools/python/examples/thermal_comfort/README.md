# Thermal-comfort workflow examples

These scripts implement the two postprocessing boundaries deliberately kept
separate in `udcomf`:

1. `udales_to_comfort.py` reads saved uDALES results and writes one exchange
   input file per receptor height.
2. `calculate_indices.py` reads exchange files and calculates MRT, PET, UTCI,
   and WBGT with one set of person and sensor parameters.

Neither script runs the LES. The first script reuses saved uDALES statistics,
facet energy-balance output, shortwave forcing, longwave forcing, geometry,
and preprocessing files. Pedestrian radiation is the expensive postprocessing
step because visibility must be evaluated at every requested receptor.

Every saved statistics record becomes one output record. Its value is the
solver's mean over the `tstatsdump` seconds preceding the record time, as set
in the `OUTPUT` namelist, and the solar angles stored with it are evaluated at
the midpoint of that window. The number of records is therefore whatever
`tstatstart`, `tstatsdump`, `tstatsgap` and `runtime` produced.

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
    "city_georeference.npz",
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
  "georeference_npz": "city_georeference.npz",
  "crs_attributes": {
    "grid_mapping_name": "<verified CF grid mapping name>",
    "crs_wkt": "<complete WKT describing the projected model axes>"
  },
  "vertical_datum": "<vertical datum of z_ground>",
  "model_version": "<uDALES commit used for the simulation>",
  "institution": "<contributing institution>",
  "building_representation": "<actual IBM and facet representation>",
  "terrain_convention": "<actual flat-ground convention>"
}
```

For a rotated LES, the CRS must genuinely describe the projected model axes.
Labelling rotated or offset local coordinates directly with the unrotated national projection is wrong.
The script checks shape, spacing, bounds, coordinate ranges, and constant
ground elevation, but it cannot prove that the supplied transform is correct.

## uDALES outputs to exchange files

The case directory name must be its experiment number and contain gathered
`stats_tkslice.<expnr>.nc` or `stats_t.<expnr>.nc`, complete facet output, and
the inputs documented in `tools/python/udcomf/README.md`.

```bash
PY=tools/python/.venv/bin/python
CASE=/path/to/experiments/1

$PY tools/python/examples/thermal_comfort/udales_to_comfort.py \
  --case-dir "$CASE" \
  --metadata "$CASE/thermal_comfort_metadata.json" \
  --statistics-source stats_tkslice \
  --tile-size 256 \
  --n-mu 8 \
  --n-azimuth 32
```

All `receptor_heights` configured in `namoptions` are processed by default.
Use `--heights 1.1 1.5 2.0` to select a configured subset. Existing final files
are never replaced unless `--overwrite` is supplied. Matching radiation tile
checkpoints are reused by default; `--no-resume` disables that behavior.

For a simulation that was warm-started under several case numbers, repeat
`--facet-source` for each source case or `facEB` file, and `--stats-source`
for each source case or gathered `stats_tkslice`/`stats_t` file. The facet
merge writes complete `facEB`/`facT` files; the statistics merge keeps only
the comfort subset (`tha`, `pabs`, `qt`, `ql` at the scalar levels bracketing
the configured receptor heights, `ws_local`, `ws_10`) and writes it as
`stats_tkslice.<expnr>.nc` of the case:

```bash
$PY tools/python/examples/thermal_comfort/udales_to_comfort.py \
  --case-dir "$CASE" \
  --metadata "$CASE/thermal_comfort_metadata.json" \
  --facet-source "$CASE/facEB.2.nc" \
  --facet-source "$CASE/facEB.3.nc" \
  --facet-source "$CASE/facEB.4.nc" \
  --stats-source /path/to/experiments/2 \
  --stats-source /path/to/experiments/3 \
  --stats-source /path/to/experiments/4
```

The merged files must not already exist in the case unless
`--overwrite-merge` is given, and the case's own files cannot be sources.

When the pedestrian atmosphere and both radiation files already exist and
have passed their earlier checks, assemble only the exchange files with:

```bash
$PY tools/python/examples/thermal_comfort/udales_to_comfort.py \
  --case-dir "$CASE" \
  --metadata "$CASE/thermal_comfort_metadata.json" \
  --skip-atmosphere \
  --skip-radiation
```

For a long headless radiation run:

```bash
LOG="$CASE/udcomf_postprocess.log"
nohup "$PY" -u tools/python/examples/thermal_comfort/udales_to_comfort.py \
  --case-dir "$CASE" \
  --metadata "$CASE/thermal_comfort_metadata.json" \
  --statistics-source stats_tkslice \
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

Run the index calculation on one or several exchange files. With the
documented default person and sensor parameters:

```bash
$PY tools/python/examples/thermal_comfort/calculate_indices.py \
  "$CASE/thermal_comfort_inputs_udales.h1p1.1.nc"
```

Several height files can be processed together. Each output is placed beside
its input:

```bash
$PY tools/python/examples/thermal_comfort/calculate_indices.py \
  outputs/thermal_comfort_inputs_udales.h1p1.1.nc \
  outputs/thermal_comfort_inputs_udales.h1p5.1.nc
```

The output name is `thermal_comfort_indices_<model>.h<height>.nc`. It contains
one `(x,y,time)` record per exchange record for `mrt`, `pet`, `utci`, and
`wbgt`, plus globe and natural-wet-bulb temperatures and validity flags.

To use a different person or sensor configuration, put only the changed
`ComfortParameters` fields in one JSON file. For example:

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

Use the same parameter file for every exchange file that is to be compared.
Inputs outside an index's validated domain are flagged and written as missing;
the calculation does not clip them into range.

## Scientific scope

The equations, fixed reference cases, and current validation limits are in
`tools/python/udcomf/VALIDATION.md`. The full input contract, restart rules,
quadrature assumptions, and current vegetation limitation are documented in
`tools/python/udcomf/README.md`.
