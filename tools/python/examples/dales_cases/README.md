# DALES validation cases for the moist-thermodynamics port

Tooling for the tier-2 validation cases `examples/301` (BOMEX), `examples/302` (smoke) and
`examples/303` (DYCOMS-II RF01). Python environment: `module purge; module load tools/prod
Python/3.9.6-GCCcore-11.2.0; source ~/udales/.venv/bin/activate` (numpy, matplotlib, netCDF4).

| script | purpose |
|---|---|
| `convert_dales_case.py <dales_case_dir> <expnr> <out_dir>` | converts DALES `prof.inp.001`, `lscale.inp.001` (8 columns to the 10-column uDALES layout, `pgx = pgy = 0`, `dthlrad` to `thlpcar`) and `scalar.inp.001`; prints the namoptions values to transcribe, with `wtsurf`/`wqsurf` negated (decision D1: uDALES fluxes are positive downward) and `isurf 3/4` mapped to `BCbotm 5/4`; re-reads its output and checks row 1 and all rows against the input |
| `make_rf01.py <out_dir> [--dz 5 --ktot 300 --itot 96 --dx 35]` | writes the RF01 `prof.inp.303`/`lscale.inp.303` and prints the surface values (`z0`, `wtsurf`, `wqsurf`) |
| `compare_dales.py <case> <udales_run_dir> <dales_results_dir\|none> <out_dir>` | time series, profiles, `summary.md` |

Case inputs were produced with
`convert_dales_case.py $SCRATCH/dales/cases/bomex 301 examples/301`, the same for `smoke`/302, and
`make_rf01.py examples/303`. The namoptions files were written by hand from the printed values.

## Published ranges

`PUBLISHED` in `compare_dales.py` maps `(case, variable)` to `(t0_h, t1_h, lo, hi, source)`. It is empty.
The papers (Siebesma et al. 2003; Bretherton et al. 1999; Stevens et al. 2005) were not available when
the tool was written, and the PyCLES and DALES files hold no published ranges, so no number could be cited.
Add a band only with a value read from a paper or table and its citation. Until then the summary table reports
`no cited band` for every variable, and a pass/fail verdict needs the bands filled in. The default averaging
windows (BOMEX 3-6 h, smoke 2-3 h, RF01 3-4 h) and entrainment fit windows (smoke 1-3 h, RF01 1-4 h) are
analysis choices, not values from the papers; the BOMEX and RF01 ones follow the paper conventions as recalled,
not re-checked. Override with `--t0/--t1/--fit0/--fit1`.

## DALES versus uDALES variables

DALES `tmser.001.nc`: `lwp_bar` [kg/m2], `cfrac`, `zi`, `ustar`, `we`. DALES `profiles.001.nc`: `thl`, `qt`,
`ql`, `wtlr`/`wqtr` (resolved fluxes). uDALES `thermostat_t`/`thermostat_z`: see `NAMTHERMOSTAT`. The
definitions of `zi` differ (uDALES: height of the maximum `dthl/dz`), so compare trends.
The uDALES `thermostat_z` fluxes `wthl`, `wqt` are resolved fluxes only.

## Running the cases

```bash
tools/build_executable.sh icl release
tools/hpc_execute.sh examples/301    # likewise 302, 303; NCPU in config.sh = nprocx * nprocy
python tools/python/examples/dales_cases/compare_dales.py 301 <outdir>/301 \
    $SCRATCH/dales/cases/bomex/results_3.2.0 plots301
```

Caveat for BOMEX: the DALES `results_3.2.0` reference ran without the prescribed -2 K/day radiative
cooling (`tltendrad` is zero), while `examples/301` applies it through `thlpcar`. For a like-for-like
check against 3.2.0, set the `thlpcar` column of `lscale.inp.301` to 0. See `examples/301/RESULTS.md`.

Smoke-test timings (60 s of simulated time, login node, release build; main-loop CPU time per simulated
minute): see the Task 7 report.
