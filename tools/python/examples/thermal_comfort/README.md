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

## Pipeline overview

```text
                 USER RUNS                        FILES PRODUCED (case directory, <expnr> = case number)
                 =========                        ==============================================

 namoptions.<expnr>  ─────────────────────────►  settings: receptor_heights, tstatsdump/tstatsgap/tstatstart,
 (edited by hand)                                 ltdump|ltkslicedump, ltempeq, lmoist, lEB + lwriteEBfiles,
                                                  ltimedepsw/ltimedeplw, year..second, longitude, latitude,
                                                  elevation, xazimuth, isolarpos
        │
        ▼
 1. tools/write_inputs.sh  (preprocessing) ───►  geometry STL, facets.inp, factypes.inp, solid_c.txt, prof.inp,
                                                  shortwave_forcing.<expnr>.nc, timedepsw.inp, timedeplw.inp
        │
        ▼
 2. run uDALES ───────────────────────────────►  per-rank stats_t / stats_tkslice, facEB, ...
        │
        ▼
 3. tools/gather_outputs.sh ──────────────────►  stats_t.<expnr>.nc  or  stats_tkslice.<expnr>.nc
                                                  facEB.<expnr>.nc
        │
        ▼
 4. tools/ud_run_thermal_comfort.sh <case> [l|f|c] [system] -- --heights ...
      │  (shell launcher: background / foreground / PBS or Slurm job)
      ▼
   udcomf/ud_run_thermal_comfort.py          log: thermal_comfort.<expnr>.log
      │
      ├─ preflight ─── namoptions, prof.inp ─────────────────────────► report; stops on blockers
      │
      ├─ merge (only --data-source) ── facEB + statistics of each
      │                                 warm-start segment ───────────► merged_facEB.<expnr>.nc
      │                                                                 merged_stats_tkslice.<expnr>.nc
      │
      ├─ atmosphere ── stats_t | stats_tkslice | merged_stats_tkslice ► pedestrian_atmosphere.h<h>.<expnr>.nc
      │
      ├─ radiation ─── STL, facets/factypes, solid_c, shortwave_forcing,
      │                timedepsw, timedeplw, facEB | merged_facEB ────► pedestrian_shortwave.h<h>.<expnr>.nc
      │                (ray tracing; checkpointed, resumes on rerun)    pedestrian_longwave.h<h>.<expnr>.nc
      │
      ├─ exchange ──── the three per-height files
      │                + namoptions site and start date ──────────────► thermal_comfort_inputs.h<h>.<expnr>.nc
      │
      └─ indices ───── exchange file + comfort parameters JSON ───────► thermal_comfort_indices.h<h>.<expnr>.nc
```

Steps 1 to 3 are the normal uDALES workflow; step 4 is the only new action.
The two step scripts `udcomf/udales_to_comfort.py` (preflight to exchange)
and `udcomf/calculate_indices.py` (indices) can also be run alone, for
example to recompute the indices with another parameter file without
repeating the radiation tracing.

### Files from the user's perspective

| file | role | user edits it? |
|---|---|---|
| `namoptions.<expnr>` | all settings, including the comfort-relevant output switches and the site | yes |
| `tools/ud_run_thermal_comfort.sh` | launches the chain so that it survives a closed terminal | no, run it |
| `examples/thermal_comfort/comfort_parameters.json` | copy of the defaults; edit and pass with `--parameters` | optional |
| `<case>/thermal_comfort.<expnr>.log` | one log of the chain: stage headers with indented details | read |
| `<case>/thermal_comfort_indices.h<h>.<expnr>.nc` | the result, one file per receptor height `h` | read |

### Stages: inputs, outputs and their content

| stage | reads | writes | content of the output |
|---|---|---|---|
| preflight | `namoptions`, `prof.inp` | nothing | blockers and cautions in the log |
| merge (optional) | from every segment: `facEB.<case>.nc`, `<statistics-source>.<case>.nc` | `merged_facEB`, `merged_stats_tkslice` | facet times `t` and `LWout` only; `tha`, `pabs`, `qt`, `ql` at the scalar levels bracketing the heights, `ws_local`, `ws_10`; gaps recorded in attributes |
| atmosphere | statistics file: `time`, `tha`, `pabs`, `qt`, `ql`, `ws_local`, `ws_10`; `tstatsdump` | `pedestrian_atmosphere` | per record and cell at the receptor height: `ta`, `pabs`, `qv` (interpolated between the bracketing levels), `ws_local`, `ws_10`; time bounds; pedestrian mask |
| radiation | STL + `facets.inp` + `factypes.inp` (albedo), `solid_c.txt` (fluid cells), `shortwave_forcing` (dni, dsky, sun angles), `timedepsw` (netsw per facet), `timedeplw` (LWsky), `facEB` (`t`, `LWout`) | `pedestrian_shortwave`, `pedestrian_longwave` | per record and cell: `sw_direct_normal` and six `sw_nondirect_*face` planes; six `lw_*face` planes; each the mean over the window `[t - tstatsdump, t)` |
| exchange | the three per-height files; `namoptions` longitude, latitude, elevation, start date, `xazimuth`, `isolarpos` | `thermal_comfort_inputs` | the 18 fields above on the model grid with cell bounds, time bounds, `solar_zenith` and `solar_azimuth` at the window midpoint, `pedestrian_mask`, `receptor_height`, site attributes |
| indices | exchange file; comfort parameters JSON | `thermal_comfort_indices` | `mrt`, `pet`, `utci`, `wbgt` in °C, `globe_temperature`, `natural_wet_bulb_temperature`, `validity_flags` bit field, `validity_counts`, the parameters used |

Every record in every file is a mean over the `tstatsdump` seconds preceding
its time stamp; the indices are computed from those means, not averaged
afterwards.

### Which input each index needs

| index | exchange fields | upstream files they come from |
|---|---|---|
| MRT | 13 radiation planes, `solar_zenith` | `shortwave_forcing`, `timedepsw` + albedo, `timedeplw`, `facEB` `LWout`; geometry for visibility |
| UTCI | MRT, `ta`, `pabs`, `qv`, `ws_10` | the above plus statistics `tha`, `pabs`, `qt`, `ql`, `ws_10` |
| PET | MRT, `ta`, `pabs`, `qv`, `ws_local` | the above plus statistics `ws_local` |
| WBGT | `ta`, `pabs`, `qv`, `ws_local`, 13 radiation planes, `solar_zenith`; not MRT | same sources as PET |

With `--indices mrt` the index step requires only the radiation planes of the
exchange file. `facT` is read by no stage.

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
