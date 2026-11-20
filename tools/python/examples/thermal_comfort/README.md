# Thermal comfort from a uDALES simulation

Three steps give pedestrian thermal-comfort maps (MRT, PET, UTCI, WBGT) from a
uDALES case:

1. **Run the simulation** with the output settings described in
   [`tools/python/udcomf/README.md`](../../udcomf/README.md): time-averaged
   statistics (`ltdump` or `ltkslicedump`), temperature and moisture
   equations, the facet energy balance with its output files, time-dependent
   shortwave and longwave forcing, and the pedestrian receptor heights
   (`receptor_heights` in `&OUTPUT`).
2. **Gather the outputs** with `tools/gather_outputs.sh <case>` as for any run.
3. **Run the thermal-comfort chain** with one command:

```bash
# from the top-level project directory
u-dales/tools/ud_run_thermal_comfort.sh <PATH_TO_CASE> [start] [system] -- --heights 1.1 1.5
```

`start` is `l` (background on this node, default), `f` (foreground) or `c`
(scheduler job); `system` is `common` (default), `icl` or `ecmwf`. Everything
after `--` goes to the Python wrapper `tools/python/udcomf/ud_run_thermal_comfort.py`;
`--heights` lists the receptor heights to process and is required.

The chain runs the preflight check of `namoptions` and the gathered files,
merges warm-start segments if asked, extracts the receptor-level means,
traces the receptor shortwave and longwave radiation (the slow step, with
checkpoints that let a rerun resume), assembles one exchange file per height
and calculates the indices. The results are

```text
<case>/thermal_comfort_indices.h<height>.<expnr>.nc
```

with `(x, y, time)` maps of `mrt`, `pet`, `utci` and `wbgt` in degrees Celsius,
the WBGT sensor temperatures and a `validity_flags` bit field. One log of the
whole chain is `<case>/thermal_comfort.<expnr>.log`.

## Common options

```bash
# a subset of the indices
... -- --heights 1.1 --indices mrt,utci

# your own person, sensor and domain parameters
... -- --heights 1.1 --parameters u-dales/tools/python/examples/thermal_comfort/comfort_parameters.json

# a simulation that was warm-started under several case numbers
... -- --heights 1.1 --data-source /cases/001 --data-source /cases/002 --data-source /cases/003

# statistics saved as receptor slices instead of the full 3-D stats_t
... -- --heights 1.1 --statistics-source stats_tkslice

# compute-node job on ECMWF with 32 cores (resources via environment or the case's config.sh)
COMFORT_NCPU=32 COMFORT_MEM=64gb u-dales/tools/ud_run_thermal_comfort.sh <case> c ecmwf -- --heights 1.1
```

## Comfort parameters

Every person, sensor and domain value of the index calculation lives in a
JSON file. The packaged defaults, `tools/python/udcomf/thermalcomfort/comfort_defaults.json`,
are what the chain uses when `--parameters` is not given; do not edit that
file. `comfort_parameters.json` in this directory is a copy meant for editing:
change the values you need and pass the file with `--parameters`. A parameter
file must contain every key of the defaults and no others; the values used are
stored in the output attribute `comfort_parameters`.

## Tutorial notebook

`thermal_comfort_tutorial.ipynb` walks through the chain on a gathered case:
the preflight report, running the wrapper, and reading and plotting the
exchange and index files. It needs a completed and gathered simulation; set
the case path in its first code cell.

## Running the steps separately

The two scripts `tools/python/udcomf/udales_to_comfort.py` (uDALES outputs to
exchange files) and `tools/python/udcomf/calculate_indices.py` (exchange files
to indices) can be run on their own, for example to recompute the indices with
another parameter file without repeating the radiation step. Their options and
the file formats are described in
[`tools/python/udcomf/WORKFLOW.md`](../../udcomf/WORKFLOW.md).
