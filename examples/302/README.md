# 302: smoke cloud (Bretherton et al. 1999)

Full-size validation case for the scalar radiation absorber (`rad_absorber = 'sv'`), the dry MO surface and the
`zi_sv` thermostat diagnostic (`zi_sv_index = 1`, `zi_sv_thres = 0.5`). Inputs are converted from the DALES case `cases/smoke`
(`convert_dales_case.py`), grid and profiles unchanged. Run with `tools/hpc_execute.sh examples/302`
(16 ranks, 8 x 2); compare with `compare_dales.py` (entrainment rate from `zi_sv`, D = 0).

## Settings and sources

| setting | value | source |
|---|---|---|
| grid | 64 x 64 x 250, dx = dy = 50 m, dz = 5 m (z up to 1250 m) | DALES `itot, jtot, kmax, xsize = ysize = 3200`; z from `prof.inp.001` |
| duration | 10800 s (3 h) | DALES `runtime` |
| time step | `ladaptive`, `courant = 0.5`, `dtmax = 10` | DALES `ladaptive, courant`; `dtmax` chosen here |
| `ps`, `thls` | 100000 Pa, 288 K | DALES |
| moisture | `lmoist = .false.` | DALES |
| surface | `BCbotm = 4` (MO, u* from MO), `z0 = z0h = 0.1` m, `wtsurf = wqsurf = 0` | DALES `isurf = 4`, `z0`, `wtsurf`, `wqsurf`; `ustin = 0.2` is unused with `BCbotm = 4`; `z0h` set equal to `z0` |
| scalar | `nsv = 1`, `iadv_sv = 7` (kappa), `scalar.inp.302` | DALES `nsv`, `iadv_sv = 7`, `scalar.inp.001` (smoke 1 up to 700 m, then a ramp to 0); `wsvsurf = 0` is the default |
| advection | `iadv_mom = iadv_tke = 2`, `iadv_thl = 7` | DALES `iadv_*`; DALES `iadv_qt = 7` becomes 2: uDALES supports only 2 for qt, and qt is inert in this dry case |
| radiation | `iradiation = 2`, `rad_longw`, `rad_absorber = 'sv'`, `rad_absorber_sv = 1`, `rka = 0.02`, `dlwtop = 60`, `dlwbot = 0` | DALES `namoptions.001` (`rad_smoke`, `isvsmoke = 1`; `dlwbot` is the default 0) |
| `timerad` | 0 | task spec; the DALES file has 600 s |
| Coriolis, sponge | `lcoriol = .false.`, `igrw_damp = 0` | DALES `lcoriol = .false.`, `igrw_damp = 0` |
| SGS | `loneeqn = .true.`, `ldelta = .false.` | as DALES (default) |
| thermostat | `lthermostat`, `tthermostat = 60` s | task spec |
| initial noise | `lrandomize`, `randu = 0.1` m/s on u, v, w, levels 1-140 (z < 700 m), `irandom = 43` | chosen here, see deviations |
| decomposition | `nprocx = 8`, `nprocy = 2` (`ktot = 250` is divisible by 2 only) | |

## Deviations from DALES

- No `thl` noise: DALES perturbs thl (`randthl = 0.1`); uDALES does not apply `randthl`. Uniform noise on
  u, v, w is used instead.
- `iadv_qt = 2` (see the table); no 5th-order scheme in uDALES; SGS details differ.
- Radiation updated every step rather than every 600 s.
- The `lscale.inp` of the case is all zero except the grid; `thlpcar = 0`, the cooling is the smoke radiation.
- Output: only the thermostat files.

## Concerns

- The smoke profile in `scalar.inp.302` has a sharp edge at 700 m; kappa advection keeps it monotone, but
  the thermostat `zi_sv` (smoke = 0.5 crossing) is quantised by the 5 m grid.
