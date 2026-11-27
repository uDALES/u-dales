# Thermal-comfort command-line workflow

## One command for the whole chain

After the simulation has run and its outputs are gathered, one command gives
the index maps:

```bash
# from the top-level project directory
u-dales/tools/ud_run_thermal_comfort.sh <PATH_TO_CASE> [l|f|c] [common|icl|ecmwf] -- --heights 1.1 1.5 [options]
```

The shell script runs `tools/python/udcomf/ud_run_thermal_comfort.py` in the
background, in the foreground or as a scheduler job, so that it survives a
closed terminal; everything after `--` is passed to the Python wrapper, whose
options are the union of the two step scripts described below
(`--heights` required; `--indices`, `--parameters`, `--statistics-source`,
`--data-source`, the radiation settings, `--skip-*`, `--overwrite`). The
wrapper writes one log of the whole chain, `<case>/thermal_comfort.<expnr>.log`,
with a header line per stage and the details indented under it, and finishes
with `thermal_comfort_indices.h<height>.<expnr>.nc` in the case directory.
A short user guide with the three steps, the common options and a tutorial
notebook is in `tools/python/examples/thermal_comfort/`.

## The two steps

The two step scripts in this directory implement the two postprocessing
boundaries deliberately kept separate in `udcomf`; the wrapper calls them
in-process, and they can be run on their own, for example to recompute the
indices with another parameter file without repeating the radiation step:

1. `udales_to_comfort.py` reads saved uDALES results and writes one exchange
   input file per receptor height.
2. `calculate_indices.py` reads exchange files and calculates MRT, PET, UTCI
   and WBGT, or any subset of them, with one set of person, sensor and domain
   parameters read from a JSON file.

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

## Site information

Everything geographic the pipeline needs is already in `namoptions`: the site
`longitude`, `latitude` and `elevation` and the UTC start date for the solar
position, and `xazimuth` for the rotation of the model x-axis from true
north. The exchange file carries the model grid in metres and records these
values as global attributes. No other metadata file is read.

## uDALES outputs to exchange files

The case directory name must be its experiment number and contain gathered
`stats_t.<expnr>.nc` (the default `--statistics-source`) or
`stats_tkslice.<expnr>.nc` (`--statistics-source stats_tkslice`), complete
facet output, and the inputs documented in `tools/python/udcomf/README.md`.

```bash
PY=tools/python/.venv/bin/python
CASE=/path/to/experiments/1

$PY tools/python/udcomf/udales_to_comfort.py \
  --case-dir "$CASE" \
  --heights 1.1 1.5 \
  --statistics-source stats_tkslice \
  --tile-size 256 \
  --n-mu 8 \
  --n-azimuth 32
```

`--heights` is required and lists the receptor heights to process; each must
be one of the `receptor_heights` saved by the simulation. All outputs are
written into the case directory. Existing final files are never replaced
unless `--overwrite` is supplied. Matching radiation tile checkpoints are
reused by default; `--no-resume` disables that behavior.

For a simulation that was warm-started under several case numbers, repeat
`--data-source` with the case directory of each segment, including the case
itself if its own records belong to the chain. From every directory the files
`facEB.<case>.nc` and `<statistics-source>.<case>.nc` are read following the
uDALES naming rule; the statistics file is the one named by
`--statistics-source` (default `stats_t`) and must exist in every segment, no
other file is tried. Only what the indices need is merged:
the facet timestamps `t` and emitted longwave `LWout` into
`merged_facEB.<expnr>.nc`, and the comfort subset of the statistics (`tha`,
`pabs`, `qt`, `ql` at the scalar levels bracketing the requested receptor
heights, `ws_local`, `ws_10`) into `merged_stats_tkslice.<expnr>.nc`. The
gathered files are never modified and `facT` is not needed. The rest of the
run reads the two merged files:

```bash
$PY tools/python/udcomf/udales_to_comfort.py \
  --case-dir "$CASE" \
  --heights 1.1 \
  --data-source /path/to/experiments/1 \
  --data-source /path/to/experiments/2 \
  --data-source "$CASE"
```

The merged files are rebuilt on every run with `--data-source`; they are
cheap and a deterministic function of the segments.

When the pedestrian atmosphere and both radiation files already exist and
have passed their earlier checks, assemble only the exchange files with:

```bash
$PY tools/python/udcomf/udales_to_comfort.py \
  --case-dir "$CASE" \
  --heights 1.1 1.5 \
  --skip-atmosphere \
  --skip-radiation
```

For a long headless radiation run:

```bash
LOG="$CASE/udcomf_postprocess.log"
nohup "$PY" -u tools/python/udcomf/udales_to_comfort.py \
  --case-dir "$CASE" \
  --heights 1.1 1.5 \
  --statistics-source stats_tkslice \
  > "$LOG" 2>&1 < /dev/null &
printf 'PID %s, log %s\n' "$!" "$LOG"
```

The outputs are:

```text
pedestrian_atmosphere.h<height>.<expnr>.nc
pedestrian_shortwave.h<height>.<expnr>.nc
pedestrian_longwave.h<height>.<expnr>.nc
thermal_comfort_inputs.h<height>.<expnr>.nc
```

## Exchange files to comfort indices

Run the index calculation on one or several exchange files. With the
parameters of `tools/python/udcomf/thermalcomfort/comfort_defaults.json`:

```bash
$PY tools/python/udcomf/calculate_indices.py \
  "$CASE/thermal_comfort_inputs.h1p1.1.nc"
```

Several height files can be processed together. Each output is placed beside
its input, with `thermal_comfort_inputs` in the name replaced by
`thermal_comfort_indices`:

```bash
$PY tools/python/udcomf/calculate_indices.py \
  outputs/thermal_comfort_inputs.h1p1.1.nc \
  outputs/thermal_comfort_inputs.h1p5.1.nc
```

The output `thermal_comfort_indices.h<height>.<expnr>.nc` contains one
`(x,y,time)` record per exchange record for each selected index in degrees
Celsius, globe and natural-wet-bulb temperatures when WBGT is selected, and
validity flags. It is always written in the directory of its exchange file,
whose name must start with `thermal_comfort_inputs`; the script prints the
full path of every result when it finishes.

`--indices` takes a comma-separated subset. Only the exchange fields those
indices need must be present:

| index | exchange fields read |
|---|---|
| `mrt` | the thirteen radiation planes and `solar_zenith` |
| `utci` | MRT inputs, `ta`, `pabs`, `qv`, `ws_10` |
| `pet` | MRT inputs, `ta`, `pabs`, `qv`, `ws_local` |
| `wbgt` | `ta`, `pabs`, `qv`, `ws_local`, the radiation planes, `solar_zenith` |

```bash
$PY tools/python/udcomf/calculate_indices.py \
  --indices mrt,wbgt "$CASE"/thermal_comfort_inputs.h*.1.nc
```

Every person, sensor and domain value of the calculation (six-direction body
factors and optical properties, the PET reference person, the UTCI
operational domain, the WBGT globe and wick geometry and properties, the
sensor solver settings and the ISO 7243 weights) is read from
`comfort_defaults.json`; nothing is fixed in the code. To use other values,
copy that file, edit the entries and pass the copy:

```bash
cp tools/python/udcomf/thermalcomfort/comfort_defaults.json my_comfort_parameters.json
# edit my_comfort_parameters.json
$PY tools/python/udcomf/calculate_indices.py \
  --parameters my_comfort_parameters.json \
  "$CASE"/thermal_comfort_inputs.h*.1.nc
```

The parameter file must contain every key of the defaults file and no others;
the values used are stored in the output attribute `comfort_parameters`. Use
the same parameter file for every exchange file that is to be compared.
Inputs outside an index's validated domain are flagged and written as missing;
the calculation does not clip them into range.

## Scientific scope

The equations, fixed reference cases, and current validation limits are in
`tools/python/udcomf/VALIDATION.md`. The full input contract, restart rules,
quadrature assumptions, and current vegetation limitation are documented in
`tools/python/udcomf/README.md`.
