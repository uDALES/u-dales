# DALES validation cases for the moist-thermodynamics port

Tooling for the tier-2 validation cases `examples/301` (BOMEX), `examples/302` (smoke) and
`examples/303` (DYCOMS-II RF01). Python environment: `module purge; module load tools/prod
Python/3.9.6-GCCcore-11.2.0; source ~/udales/.venv/bin/activate` (numpy, matplotlib, netCDF4).

| script | purpose |
|---|---|
| `convert_dales_case.py <dales_case_dir> <expnr> <out_dir>` | converts DALES `prof.inp.001`, `lscale.inp.001` (8 columns to the 10-column uDALES layout, `pgx = pgy = 0`, `dthlrad` to `thlpcar`) and `scalar.inp.001`; prints the namoptions values to transcribe, with `wtsurf`/`wqsurf` negated (decision D1: uDALES fluxes are positive downward) and `isurf 3/4` mapped to `BCbotm 5/4`; re-reads its output and checks row 1 and all rows against the input |
| `make_rf01.py <out_dir> [--dz 5 --ktot 300 --itot 96 --dx 35]` | writes the RF01 `prof.inp.303`/`lscale.inp.303` and prints the surface values (`z0`, `wtsurf`, `wqsurf`) |
| `compare_dales.py <case> <out_dir> --udales U1 [U2 ...] [--dales D1 ...]` | seed ensembles of both codes: per-run window means, ensemble means, spreads and the ratio \|mean diff\| / spread for LWP, cloud cover, zi, u* and the thl/qt/ql/w'thl'/w'qt' profiles; `timeseries.png`, `profiles.png`, `summary.md`, `summary.json`. The one-pair form `compare_dales.py <case> <udales_run_dir> <dales_dir\|none> <out_dir>` still works |
| `make_dales_baseprof.py <prof.inp.001> <ps> <thls> <qts> <baseprof.inp.001>` | constant-density DALES `baseprof.inp` (`ibas_prf = 5`) with `pb` from the uDALES hydrostatic integral, for Boussinesq-equivalent DALES runs |

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
`ql`, `wthlr`/`wqtr` (resolved fluxes). uDALES `bldump_t`/`bldump_z`: see `lbldump` in the OUTPUT namelist overview.
The `zi` diagnostics of the two codes differ (DALES tmser: thv gradient on a 4-column composite by default;
uDALES bldump: column-wise max `dthl/dz`), so `compare_dales.py` recomputes zi for both codes from the
block-mean profiles with one definition: the half level of the maximum `dthl/dz`, and for 303 also `zi_qt8`,
the highest crossing of qt = 8 g/kg. DALES profiles are block means over `timeav` (600 s) labelled at the
block end; the tool averages the 60 s uDALES snapshots into the same blocks and uses the blocks that lie
inside the window (a 3-6 h window starts with the block 3h00-3h10 in both codes; the earlier version also
took the DALES block 2h50-3h00, and read the resolved w'thl' under a name DALES does not write, `wtlr`).
The uDALES `bldump_z` fluxes `wthl`, `wqt` are resolved fluxes only.

## Running the cases

```bash
tools/build_executable.sh icl release
tools/hpc_execute.sh examples/301    # likewise 302, 303; NCPU in config.sh = nprocx * nprocy
python tools/python/examples/dales_cases/compare_dales.py 301 <outdir>/301 \
    $SCRATCH/dales/cases/bomex/results_3.2.0 plots301
```

Caveat for BOMEX: the DALES `results_3.2.0` reference ran without the prescribed -2 K/day radiative
cooling (`tltendrad` is zero; reported as dalesteam/dales#348), while `examples/301` applies it through
`thlpcar`. Compare against a run with current DALES `dev` (which applies it), or set the `thlpcar` column of
`lscale.inp.301` to 0 for a check against 3.2.0. See `examples/301/RESULTS.md`.

## Numerics settings of the cases (301-303)

The converter prints DALES `courant`, `randthl`, `randqt`, `randu` and `xlat` for transcription, but most of
them do not carry over unchanged. Every numerics setting and its reason is in the table of each case README
(`examples/301`, `302`, `303`). In short:

- Noise: DALES perturbs thl (`randthl`) and qt (`randqt`) at all levels and applies no u, v, w noise, so the
  cases set `randu = 0`, `randthl`, `randqt` and `krand = ktot` explicitly (uDALES defaults: `randu = 0.01`, `randthl = randqt = 0`).
- `xlat` is copied from the DALES case file (15 for BOMEX); `lqlnr = .true.` as DALES (which always iterates).
- Time step: uDALES `courant` is an L1 sum, DALES's an L2 norm with default 1.0, and uDALES has no Galilean
  shift (`cu`, `cv`); so `courant` is tuned until the mean dt matches a DALES run (301: 0.55; see `examples/301/README.md`),
  and `diffnr` = DALES `peclet` (default 0.15) x (`min(dx, dy, dz)^2 (1/dx^2 + 1/dy^2 + 1/dz^2)`). Do not transcribe DALES's `courant`.
- Fixed step for an identical-dt comparison: `ladaptive = .false.` with `dtmax` as the step (commented block in
  `namoptions.301` and `namoptions.303`); DALES needs the same `dtmax` and `ladaptive = .false.`.

## DALES inputs for the close-match comparison (`dales_inputs/`)

DALES dev 075ac4e input files that run BOMEX and RF01 with the numerics of `examples/301` and `examples/303`
(settings table in each case's `RESULTS.md`): `dales_inputs/bomex` and `dales_inputs/rf01` hold `namoptions.001`,
`prof.inp.001`, `lscale.inp.001` and `baseprof.inp.001`. They set a fixed time step (`ladaptive = .false.`,
`dtmax` = 2.0 s and 0.3 s, the same as the fixed-dt blocks of the uDALES namoptions), no Galilean shift
(`cu = cv = 0`), `ibas_prf = 5` with a constant-density base profile (rhobf = 1, `pb` from the uDALES hydrostatic
integral, written by `make_dales_baseprof.py`), cd2 advection for every variable, thl/qt noise only
(`krandumax = 0`) and, for RF01, `wqsurf` with rlv = 2.53e6. Only `irandom` changes between ensemble members.
For an adaptive run set `ladaptive = .true.` and `dtmax = 20.` (the fixed `dtmax` otherwise caps the step).

The DALES build needs `dales_inputs/dales-075ac4e-user-es.patch`, applied to a copy of the DALES tree after copying
`cases/dycoms_rf02/moduser.f90` to `src/moduser.f90`:
- `src/moduser.f90` `rad_user`: the RF01 free-troposphere term in the flux-divergence form at half levels with zi
  reset per column, identical to uDALES `src/cases/case_dycoms_rf01.f90`. The dycoms_rf02 version uses the
  pointwise `(z - zi)^(-2/3)` form (singular next to zi), has the sign error of dalesteam/dales#347, loops
  `i = 2, j1` and does not reset zi per column.
- `src/modthermodynamics.f90`: the liquid saturation-vapour-pressure table uses Tetens with the uDALES constants
  instead of Murphy-Koop (about 0.1 % in es), so both codes use the same es.
Then apply `dales_inputs/dales-075ac4e-iradiation10.patch` on top:
- `src/modradiation.f90`: DALES dev 075ac4e switches `rad_longw`, `rad_shortw` and `rad_smoke` off for
  `iradiation = 10` (`irad_user`), so the `radpar` call in the `irad_user` branch is unreachable and only
  `rad_user` runs. Without this hunk the RF01 reference has **no cloud-top longwave cooling** (F0 = 70, F1 = 22
  W/m2), only the free-troposphere term. A DALES RF01 reference run without this hunk has this defect.
- von Karman constant: DALES uses `fkar = 0.4`; uDALES used 0.41 until its default `fkar` became 0.40 (as DALES and
  PALM), so no edit is needed. To compare with a uDALES run that sets `fkar = 0.41`, change `src/modglobal.f90`
  line 101 of DALES dev 075ac4e (`real,parameter :: fkar = 0.4`) to 0.41; with 0.4 against 0.41 the DALES u* is about
  2.5 % lower at the same wind.
- The `iradiation = 10` reset has been reported upstream as a comment on [dalesteam/dales#347](https://github.com/dalesteam/dales/issues/347#issuecomment-6878639996); until DALES changes it, this patch is needed.
Build as in the DALES README; with Intel 2021a on CX3, compile `modheterostats.f90` at `-O1` (ifort ICE).

Do not set `krandumax > 0` with the default `krandumin = 1` in DALES: that perturbs w at k = 1 (the surface),
which DALES never resets, so a fixed random vertical velocity (|w| up to `randu`) stays at the ground for the
whole run and `divmax` = rhobf randu / dzf is constant in the log. The first DALES RF01 reference
(`refruns/d303`) did this.
