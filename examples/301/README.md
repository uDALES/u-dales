# 301: BOMEX shallow cumulus (Siebesma et al. 2003)

Full-size validation case for the moist thermodynamics, Monin-Obukhov (MO) surface and
thermostat diagnostics. Inputs are converted from the DALES case `cases/bomex`
(`tools/python/examples/dales_cases/convert_dales_case.py`), grid and forcing unchanged.
Run on CX3 with `tools/hpc_execute.sh examples/301` (8 ranks, 4 x 2); compare with
`tools/python/examples/dales_cases/compare_dales.py` (see the README there).
Outputs of interest: `thermostat_t.301.nc`, `thermostat_z.301.nc`.

## Settings and sources

| setting | value | source |
|---|---|---|
| grid | 64 x 64 x 80, dx = dy = 100 m, dz = 40 m (z up to 3200 m) | DALES `itot, jtot, kmax, xsize, ysize`; z from `prof.inp.001` |
| duration | 21600 s (6 h) | task spec; the DALES file has `runtime = 28800` |
| time step | `ladaptive`, `courant = 0.5`, `dtmax = 10` | DALES `ladaptive`; `courant` is the DALES default; `dtmax` is an upper bound chosen here |
| `ps` | 101500 Pa | DALES `ps` |
| `thls` | 299.1 K | DALES `thls` (initial value; MO diagnoses it every step) |
| `qts` | 1.7e-2 kg/kg | first guess equal to `qt` at the first level; MO diagnoses it every step |
| `wtsurf` | -8.0e-3 K m/s | DALES `wtsurf = +8e-3` (upward), negated: uDALES prescribed fluxes are positive downward (D1) |
| `wqsurf` | -5.2e-5 kg/kg m/s | DALES `wqsurf = +5.2e-5`, negated |
| `BCbotm` | 5 (MO, u* = ustin) | DALES `isurf = 3` |
| `ustin` | 0.28 m/s | DALES `ustin` |
| `z0`, `z0h` | 2e-4 m | DALES `z0` sets both |
| Coriolis | `lcoriol`, `xlat = 14.94` (f = 0.376e-4 1/s) | Siebesma et al. 2003 f; the DALES file has `xlat = 15` (f = 0.3775e-4) |
| forcing | `lscale.inp.301` | DALES `lscale.inp.001`: `ug, vg, wfls, dqtdtls`; `dthlrad` becomes `thlpcar` (-2.315e-5 K/s) |
| radiation | `iradiation = 0`; the cooling is the `thlpcar` column | DALES (no interactive radiation in BOMEX) |
| SGS | `loneeqn = .true.` (1.5-order TKE), `ldelta = .false.` | DALES `ldelta = .false.`; one-equation closure as DALES |
| advection | `iadv_mom = iadv_tke = iadv_thl = iadv_qt = 2` | DALES `iadv_* = 2` |
| sponge | `igrw_damp = 2`, `ksp` default (top quarter, `rnu0 = 2.75e-3`) | DALES default `igrw_damp = 2`; the `ksp` default is the same formula |
| surface | `libm = .false.`, `lbottom = .true.`, `nfcts = 0`, `BCbotT = BCbotq = 1` | flat bottom without buildings |
| thermostat | `lthermostat`, `tthermostat = 60` s | task spec |
| `timerad` | 0 | task spec (radiation every step; unused with `iradiation = 0`) |
| initial noise | `lrandomize`, `randu = 0.1` m/s on u, v, w, levels 1-20 (z < 800 m), `irandom = 43` | chosen here, see deviations |
| decomposition | `nprocx = 4`, `nprocy = 2`; `jtot % nprocy = 0` and `ktot % nprocy = 0` | |

`prof.inp.301` row 1 equals the DALES row 1 column for column (checked by the converter,
which re-reads its output). The `lscale.inp.001` row at z = 3220 m (above the model top) is dropped.

## Deviations from DALES

- MO surface: `BCbotm = 5` follows DALES `isurf = 3` (u* = ustin, prescribed fluxes), so the
  spec's earlier wall-function caveat does not apply to this case.
- No Galilean shift: DALES runs with `cu = -8` (u is shifted by 8 m/s so that the domain-mean
  wind is small). uDALES has no `cu`; `prof.inp` and `lscale.inp` carry the physical winds, so
  the adaptive time step is shorter (u about 9 m/s with dx = 100 m gives dt about 5 s).
- No `thl`/`qt` noise: DALES perturbs thl (0.1 K) and qt (2.5e-5 kg/kg); uDALES reads `randthl` and
  `randqt` but does not apply them. The case is started with uniform noise on u, v, w instead.
- Advection and SGS details differ from DALES in the places uDALES differs (no 5th-order scheme, SGS
  details); BOMEX uses second-order central differences in DALES too.
- Duration 6 h instead of 8 h; `xlat` 14.94 instead of 15 (task spec).
- Output: only the thermostat files (statistics, field dumps and cross sections are off).

## Concerns

- No `thl`/`qt` noise, so spin-up differs from DALES; check the first hour.
- The `qts` and `thls` namelist values are first guesses only; the comparison should use the
  diagnosed `thls`, `qts`, `ustar` and `obl` in `thermostat_t`.
