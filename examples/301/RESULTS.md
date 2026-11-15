# Tier-2 validation results: 301 BOMEX

Note on the von Karman constant: the flat-bottom runs (the sections from "Earlier close-match validation" on) used kappa = 0.41 in uDALES (its default then) and
kappa = 0.4 in DALES (built in). The uDALES default `fkar` is now 0.40, the DALES value. BOMEX has a prescribed u*
(`BCbotm = 5`), so kappa enters only through the surface-layer gradients (see the settings table); the runs were not
repeated with `fkar = 0.40`. `namoptions.301` prescribes z0 = 2e-4 m independently of kappa, so it needs no change.

## Facet floor with the Monin-Obukhov facet surface and the first-level SGS treatment (shipped configuration)

The shipped `namoptions.301`: a floor made of one horizontal facet with `iwallmom = 4` (local Monin-Obukhov per facet
section, u* prescribed by `ustin = 0.28`, fluxes `bctfz`, `bcqfz`), `fkar = 0.40`, and in the 4096 first-level cells
above the facet the DALES first-level treatment of the one-equation SGS model (TKE shear and buoyancy sources and
`dthvdz` from the Monin-Obukhov surface-layer gradients, evaluated before `subgrid`). Runs: Release build, FFTW planner
`FFTW_ESTIMATE`, fixed dt = 2.0 s, 8 ranks (4 x 2), seeds 43-45. DALES: the runs of the sections below (DALES dev
075ac4e, fixed dt = 2.0 s, `lmostlocal = .true.`, kappa 0.4). Window 3-6 h, 600 s blocks.

**Verdict.** The near-surface difference of the section below is gone: the resolved fluxes at 40 m are within 2 % of
DALES (w'thl' 1.79e-3 against 1.82e-3 K m/s, w'qt' 1.98e-5 against 2.01e-5), and all five profiles now meet the bar
(rms difference 0.38-0.52 of the smaller rms spread; 57-90 % of the levels inside both spreads; resolved w'thl' and
w'qt' were 1.27 and 2.40 without the treatment). zi and cloud cover are within the larger spread. LWP is 0.14 g/m2
(2 %) below DALES, as without the treatment (the treatment changes it by -0.01 g/m2, 0.26 of the smaller spread); it is
within the uDALES spread (0.29, one low member) but 3.5 DALES spreads, so the strict bar is still not met for LWP.

| quantity | uDALES members | uDALES mean (spread) | DALES members | DALES mean (spread) | diff | ratio (min) | ratio (max) |
|---|---|---|---|---|---|---|---|
| LWP [g/m2] | 6.457, 6.857, 7.028 | 6.781 (0.293) | 6.926, 6.952, 6.875 | 6.918 (0.039) | -0.137 | 3.50 | 0.47 |
| cloud cover | 0.1962, 0.1938, 0.2070 | 0.1990 (0.0070) | 0.2067, 0.2005, 0.2040 | 0.2037 (0.0031) | -0.0047 | 1.53 | 0.67 |
| zi, max dthl/dz [m] | 1653, 1691, 1667 | 1670 (19) | 1647, 1669, 1651 | 1656 (12) | +15 | 1.26 | 0.77 |
| u* [m/s] | 0.28 (prescribed) | 0.28 | 0.28 | 0.28 | 0 | - | - |

| profile (z <= 2500 m) | max abs diff (z) | rms diff | rms spread U / D | rms diff / min rms spread | levels with diff < both spreads |
|---|---|---|---|---|---|
| thl [K] | 0.021 (1740 m) | 0.0068 | 0.028 / 0.018 | 0.38 | 65 % |
| qt [g/kg] | 0.026 (1740 m) | 0.011 | 0.036 / 0.023 | 0.47 | 73 % |
| ql [g/kg] | 0.0002 (620 m) | 0.00010 | 0.00023 / 0.00018 | 0.52 | 57 % |
| w'thl' resolved [K m/s] | 0.0008 (1560 m) | 0.00026 | 0.0008 / 0.0006 | 0.47 | 90 % |
| w'qt' resolved [g/kg m/s] | 0.0011 (1440 m) | 0.00043 | 0.0015 / 0.0009 | 0.48 | 89 % |

Near-surface resolved fluxes, window 3-6 h, ensemble mean (spread), zm = 40, 80, 120 m:

| run | w'thl' [K m/s] | w'qt' [kg/kg m/s] |
|---|---|---|
| facet floor with the first-level treatment | 1.79e-3 (1.5e-5), 2.84e-3 (2.8e-5), 2.45e-3 (5.2e-5) | 1.98e-5 (1.7e-7), 3.55e-5 (2.6e-7), 4.26e-5 (3.9e-7) |
| facet floor without it (section below) | 3.47e-3 (3.4e-5), 3.86e-3 (4.3e-5), 2.84e-3 (2.7e-5) | 3.18e-5 (1.2e-7), 4.31e-5 (2.9e-7), 4.56e-5 (3.4e-7) |
| DALES | 1.82e-3 (2.8e-5), 2.88e-3 (4.3e-5), 2.49e-3 (4.2e-5) | 2.01e-5 (4e-8), 3.59e-5 (9e-8), 4.28e-5 (1.1e-7) |

The near-surface resolved fluxes remain 1.5-2 % below DALES at all three levels (1-2 spreads); not investigated
further.

Hourly LWP [g/m2], ensemble mean (spread):

| code | 0-1 h | 1-2 h | 2-3 h | 3-4 h | 4-5 h | 5-6 h |
|---|---|---|---|---|---|---|
| uDALES, first-level treatment | 3.18 (0.19) | 5.21 (0.20) | 6.28 (0.67) | 6.54 (0.28) | 7.03 (0.62) | 6.78 (1.07) |
| DALES | 3.08 (0.10) | 5.14 (0.41) | 6.33 (0.84) | 6.22 (0.53) | 6.98 (0.65) | 7.55 (0.33) |

Against the facet-floor runs without the treatment (section below; 60 s blocks): LWP -0.01 g/m2 (0.26 spreads), cloud
cover -0.010 (1.4 / 1.2 spreads), zi +13 m, thl and qt profiles 1.2 / 1.1 of the smaller rms spread (0.6 / 0.7 of the
larger), resolved w'qt' 2.4 (the near-surface change).

## Facet floor without the first-level SGS treatment (earlier runs)

These runs predate the first-level SGS treatment above horizontal Monin-Obukhov facets; otherwise the configuration
is the shipped one.

The floor is one facet with `iwallmom = 4` (local Monin-Obukhov per facet section, u* prescribed by `ustin = 0.28`,
fluxes `bctfz`, `bcqfz`) and `fkar = 0.40` (default), without TKE sources at `kb`. Runs: Release build,
FFTW planner `FFTW_ESTIMATE`, fixed dt = 2.0 s, 8 ranks (4 x 2), seeds 43-45. DALES: the same DALES dev 075ac4e runs as
in the section below (fixed dt = 2.0 s, `lmostlocal = .true.` by default, kappa 0.4). Window 3-6 h, 600 s blocks.

**Verdict.** The cloud layer and the mean profiles agree within the seed spread (thl, qt, ql: rms difference 0.81-0.85
of the smaller rms spread; zi within 0.3 spreads). LWP is 0.13 g/m2 (2 %) lower than DALES, which is 3.2 spreads with
both spreads near 0.04 g/m2, so the strict bar is not met for LWP; cloud cover is within the larger spread. The resolved
fluxes do not meet the bar below about 100 m: at 40 m the resolved w'thl' is 3.47e-3 against 1.82e-3 K m/s in DALES
(w'qt' 3.2e-5 against 2.0e-5). This is the one difference left between the codes near the surface: DALES computes the
SGS TKE sources at its first level from the Monin-Obukhov surface-layer shears (and `dthvdz(1)`), uDALES has no TKE
source at `kb` next to any wall. With less SGS TKE at the first level, more of the (prescribed) surface flux is carried
by the resolved motion at 40 m. An A/B test on one seed confirms it: the earlier flat-bottom Monin-Obukhov route with
these k = 1 sources gives 1.83e-3 K m/s at 40 m (DALES 1.79e-3, same seed), the same route without them 3.47e-3, and the
facet floor 3.44e-3.

| quantity | uDALES members | uDALES mean (spread) | DALES members | DALES mean (spread) | diff | ratio (min) | ratio (max) |
|---|---|---|---|---|---|---|---|
| LWP [g/m2] | 6.754, 6.787, 6.832 | 6.791 (0.039) | 6.926, 6.952, 6.875 | 6.918 (0.039) | -0.127 | 3.25 | 3.24 |
| cloud cover | 0.2024, 0.2058, 0.2175 | 0.2086 (0.0080) | 0.2067, 0.2005, 0.2040 | 0.2037 (0.0031) | +0.0049 | 1.57 | 0.61 |
| zi, max dthl/dz [m] | 1664, 1649, 1660 | 1658 (8) | 1647, 1669, 1651 | 1656 (12) | +2 | 0.28 | 0.19 |
| u* [m/s] | 0.28 (prescribed) | 0.28 | 0.28 | 0.28 | 0 | - | - |

| profile (z <= 2500 m) | max abs diff (z) | rms diff | rms spread U / D | rms diff / min rms spread | levels with diff < both spreads |
|---|---|---|---|---|---|
| thl [K] | 0.030 (1700 m) | 0.011 | 0.013 / 0.018 | 0.85 | 35 % |
| qt [g/kg] | 0.038 (1620 m) | 0.017 | 0.021 / 0.023 | 0.81 | 38 % |
| ql [g/kg] | 0.0003 (1500 m) | 0.00013 | 0.00015 / 0.00018 | 0.84 | 24 % |
| w'thl' resolved [K m/s] | 0.0016 (40 m) | 0.0006 | 0.0005 / 0.0006 | 1.27 | 33 % |
| w'qt' resolved [g/kg m/s] | 0.0117 (40 m) | 0.0020 | 0.0008 / 0.0009 | 2.40 | 33 % |

Hourly LWP [g/m2], ensemble mean (spread):

| code | 0-1 h | 1-2 h | 2-3 h | 3-4 h | 4-5 h | 5-6 h |
|---|---|---|---|---|---|---|
| uDALES facet floor | 3.34 (0.01) | 4.40 (0.21) | 6.26 (0.39) | 7.03 (0.93) | 6.46 (0.62) | 6.88 (0.70) |
| DALES | 3.08 (0.10) | 5.14 (0.41) | 6.33 (0.84) | 6.22 (0.53) | 6.98 (0.65) | 7.55 (0.33) |

Secondary reference, the earlier flat-bottom runs below (`BCbotm = 5`, removed since; kappa 0.41): against them the
facet floor has LWP -0.16 g/m2 (0.86 of their spread), cloud cover +0.003 (0.4 spreads), zi -19 m and the same
near-surface flux difference (w'thl' +1.6e-3 K m/s at 40 m; rms over height 1.5 times the smaller rms spread).

Near-surface resolved fluxes, window 3-6 h, ensemble means (zm = 40, 80, 120 m):

| run | w'thl' [K m/s] | w'qt' [kg/kg m/s] |
|---|---|---|
| facet floor (3 seeds) | 3.47e-3, 3.86e-3, 2.84e-3 | 3.18e-5, 4.31e-5, 4.56e-5 |
| flat-bottom MO with k = 1 sources (3 seeds, kappa 0.41) | 1.83e-3, 2.88e-3, 2.48e-3 | 2.02e-5, 3.61e-5, 4.30e-5 |
| DALES (3 seeds) | 1.82e-3, 2.88e-3, 2.49e-3 | 2.01e-5, 3.59e-5, 4.28e-5 |

## Earlier close-match validation: flat-bottom Monin-Obukhov route (`BCbotm = 5`, since removed)

**Verdict.** Before the #405 fix uDALES fails the close-match bar by a wide margin (LWP 2.3 vs 6.9 g/m2, cloud
cover 0.03 vs 0.20). The cause is a uDALES bug, found here: **the SGS TKE `e120`/`e12m` halos are never exchanged
between MPI ranks** (issue #405, fixed in this change). With that fix, uDALES and DALES agree within the
seed-to-seed spread for LWP and cloud cover, and for all five profiles the rms difference is below the smaller rms
spread of the two ensembles (ratio 0.38-0.66; 49-83 % of the levels are inside both spreads). u* is prescribed in
BOMEX (`BCbotm = 5` / `isurf = 3`), so its match is not a test. zi is 19 m (half a 40 m grid cell) higher in uDALES,
1.2-1.6 spreads, not significant with 3 seeds (t = 1.7).

Runs: uDALES Release build with the e12 halo exchange (and, for comparison, the same source without it), 8 ranks
(4 x 2). DALES dev 075ac4e with `tools/python/examples/dales_cases/dales_inputs/dales-075ac4e-user-es.patch`, 16 ranks
(4 x 4), inputs `tools/python/examples/dales_cases/dales_inputs/bomex`. Three seeds per code (`irandom` 43, 44, 45).
Analysis: `compare_dales.py 301 <out> --udales ... --dales ...`. BOMEX was not rerun after the uDALES default `fkar` changed to 0.40.

### Identical settings, side by side

| setting | uDALES (examples/301) | DALES (dales_inputs/bomex) | matched? |
|---|---|---|---|
| grid, domain | 64 x 64 x 80, 100 m x 100 m x 40 m | same | yes |
| time step | fixed 2.0 s (`ladaptive = .false.`, `dtmax = 2.0`) | fixed 2.0 s (`ladaptive = .false.`, `dtmax = 2.0`; 10800 steps, checked in the log) | yes |
| time scheme | RK3 | RK3 | yes |
| advection | cd2 for u, v, w, e12, thl, qt (`iadv_* = 2`) | cd2 (`iadv_* = 2`) | yes |
| Galilean shift | none | `cu = cv = 0` | yes |
| base state | Boussinesq (rhobf = rhobh = 1) | `ibas_prf = 5`, rhobf = 1 (`baseprof.inp.001`) | yes |
| Exner function | from the hydrostatic slab-mean pressure, every step | fixed from `pb` = the uDALES hydrostatic pressure of the initial profile (`make_dales_baseprof.py`); `presf` for qsat recomputed every step in both | at t = 0; drift of order 1e-5 in exner over the run (few Pa change of the hydrostatic pressure). Not matchable without a code change |
| constants | grav 9.81, rd 287.04, rv 461.5, cp 1004, rlv 2.53e6, pref0 1e5, tmelt 273.16 | same | yes |
| saturation vapour pressure | Tetens (es0 610.78, at 17.27, bt 35.86) | Tetens, same constants (patched `esatltab`, 0.2 K table, linear interpolation) | yes (table error < 1e-4 relative) |
| saturation adjustment | exact (Newton-Raphson, `lqlnr`) | two-step linearised (DALES icethermo) | no; within 0.2 % of exact (an audit of the saturation adjustment) |
| SGS | 1.5-order TKE, `ldelta = .false.`, cf 2.5, cn 0.76, Rigc 0.25, Prandtl 1/3 (Kh/Km = 1 + 2 l/Delta), no factor 2 in dissipation, MO TKE sources at kb, grav/thvf | same model and constants (DALES log), `sgs_surface_fix = .false.` | yes |
| von Karman constant | `fkar = 0.41` | `fkar = 0.4` (built in) | no; with u* prescribed it enters only the MO surface-layer gradients (k = 1 TKE sources, dthvdz(1)); not rerun. With the uDALES default now 0.40 both codes use 0.40 |
| SGS floors | numol added to ekm/ekh; e12min 5e-5 | ekmin 1e-6; e12min 5e-5 | no; negligible (1e-5 m2/s) |
| surface | `BCbotm = 5` (MO, u* = 0.28), wtsurf 8e-3 K m/s, wqsurf 5.2e-5 up, z0 = z0h = 2e-4; one domain-mean Obukhov length | `isurf = 3`, same values; per-column Obukhov length (`lmostlocal = .true.`, the DALES dev default, not set in the inputs) | values yes; local vs domain-mean MO no (u* is prescribed, so only the surface-layer gradients see it) |
| large-scale forcing | ug, vg, wfls, dqtdtls, thlpcar (-2 K/day) from DALES `lscale.inp.001` | same file | yes |
| subsidence | upwind on the local fields (DALES lstend) | lstend | yes |
| Coriolis | xlat 15 | xlat 15 | yes |
| sponge | `igrw_damp = 2`, ksp = min(3kmax/4, kmax-15), rnu0 2.75e-3 | default `igrw_damp = 2`, same ksp and rnu0 | yes; DALES also pins thl/qt at kmax to the slab mean (inside the sponge) |
| top BC | zero gradient for thl, qt, e12 | thl, qt extrapolated with the top gradient; e12 = e12min | no; only in the top sponge cell |
| initial noise | thl +-0.1 K, qt +-2.5e-5 on all 80 levels, no u, v, w noise | same amplitudes and levels | amplitudes yes; generators differ (uDALES splitmix64 with one draw for thlm and thl0, DALES LCG with separate draws), a different realisation by design |
| initial e12 | `prof.inp` e12 column | same | yes |
| decomposition | 4 x 2 | 4 x 4 | irrelevant once the e12 halo bug is fixed |
| diagnostics | 60 s slab-mean snapshots (bldump) | 60 s samples, 600 s block means (genstat) | aligned in compare_dales.py: same 600 s blocks, same zi definition from the block-mean profiles |

### Results, fixed dt = 2 s, window 3-6 h (3 seeds each)

uDALES + e12 halo fix vs DALES. Spread = sample std over the 3 members; ratio (min) =
|mean_U - mean_D| / min(spread_U, spread_D) (bar met if < 1), ratio (max) uses the larger spread.

| quantity | uDALES members | uDALES mean (spread) | DALES members | DALES mean (spread) | diff | ratio (min) | ratio (max) |
|---|---|---|---|---|---|---|---|
| LWP [g/m2] | 7.01, 6.74, 7.11 | 6.95 (0.19) | 6.93, 6.95, 6.88 | 6.92 (0.04) | +0.03 | 0.87 | 0.18 |
| cloud cover | 0.195, 0.212, 0.208 | 0.205 (0.009) | 0.207, 0.201, 0.204 | 0.204 (0.003) | +0.002 | 0.52 | 0.18 |
| zi, max dthl/dz [m] | 1673, 1691, 1660 | 1675 (16) | 1647, 1669, 1651 | 1656 (12) | +19 | 1.64 | 1.23 |
| u* [m/s] | 0.28 (prescribed) | 0.28 | 0.28 | 0.28 | 0 | - | - |

Profile bar: the rms over height of the ensemble-mean difference, divided by the smaller of the two rms spreads
(rms over height of the per-level sample std); below 1 means the profile difference is smaller than the noise of
either ensemble on average. The last column is the stricter per-level count.

| profile (z <= 2500 m) | max abs diff (z) | rms diff | rms spread U / D | rms diff / min rms spread | levels with diff < both spreads |
|---|---|---|---|---|---|
| thl [K] | 0.035 (1660 m) | 0.011 | 0.017 / 0.018 | 0.64 | 57 % |
| qt [g/kg] | 0.050 (1660 m) | 0.015 | 0.027 / 0.023 | 0.66 | 60 % |
| ql [g/kg] | 0.0002 (860 m) | 0.00007 | 0.00025 / 0.00018 | 0.39 | 49 % |
| w'thl' resolved [K m/s] | 0.0005 (1280 m) | 0.0002 | 0.0006 / 0.0006 | 0.38 | 83 % |
| w'qt' resolved [g/kg m/s] | 0.0010 (520 m) | 0.0004 | 0.0012 / 0.0009 | 0.49 | 65 % |

Hourly LWP [g/m2], ensemble mean (spread):

| code | 0-1 h | 1-2 h | 2-3 h | 3-4 h | 4-5 h | 5-6 h |
|---|---|---|---|---|---|---|
| uDALES + fix | 3.18 (0.19) | 5.32 (0.25) | 5.84 (1.13) | 7.09 (0.67) | 6.72 (0.47) | 7.05 (0.67) |
| DALES | 3.08 (0.10) | 5.14 (0.41) | 6.33 (0.84) | 6.22 (0.53) | 6.98 (0.65) | 7.55 (0.33) |
| uDALES without the e12 fix | 0.00 | 0.00 | 0.05 | 0.46 | 1.54 | 4.95 |

uDALES without the e12 halo fix: LWP 2.32 (0.06), cloud cover 0.029 (0.003), zi 1621 (29) m;
thl 0.30 K warmer and qt 0.60 g/kg drier at 300 m; ratio 117 (LWP), 58 (cloud cover), 13-22 (profiles).

The remaining zi difference: the profiles agree within the spread at every level, and thl/qt differ most at the
inversion (1660 m, 0.035 K, 0.05 g/kg), so the uDALES inversion sits slightly higher. A parabolic sub-grid
refinement of the gradient maximum gives the same offset (1677 vs 1661 m, spreads 13 and 10 m). Welch t = 1.7 with
3 + 3 members (p about 0.2), so it is not resolved from noise. Next step if it matters: 6 or more seeds per code.

### Secondary comparison: adaptive dt (3 seeds each)

uDALES `courant = 0.55`, `diffnr = 0.2` (case file); DALES its defaults (`courant` 1.0 L2, `peclet` 0.15, `dtmax = 20`).
The steps differ: with the e12 fix uDALES is Courant-limited at hourly means 4.0, 2.7, 2.6, 2.3-2.5 s
(minimum 1.66 s), DALES runs 6.1, 4.4, 4.1, 3.7 s. The `courant = 0.55` was tuned against uDALES with
the e12 bug, whose weak convection allowed longer steps.

| quantity | uDALES + fix | DALES | ratio (min / max) | uDALES without the e12 fix |
|---|---|---|---|---|
| LWP [g/m2] | 7.08 (0.35) | 6.72 (0.16) | 2.21 / 1.04 | 1.75 (0.63) |
| cloud cover | 0.202 (0.012) | 0.199 (0.004) | 0.52 / 0.20 | 0.024 |
| zi [m] | 1681 (21) | 1667 (8) | 1.85 / 0.70 | 1623 |
| thl, qt, ql, w'thl', w'qt' profiles: rms diff / min rms spread | 0.57, 0.54, 2.05, 1.30, 1.12 | | | 23, 30, 29, 14, 14 |

With different time steps the adaptive pair is close but not within the strict bar for LWP and ql; the fixed-dt pair is.

### The e12 halo bug (finding)

`src/modboundary.f90:halos` exchanges the halos of u, v, w, thl, qt and the scalars with `exchange_halo_z`, but not
those of `e120` and `e12m`. `xm_periodic`/`ym_periodic` copy the e12 halos only when a direction is not split over
ranks. With `nprocx > 1` or `nprocy > 1` the e12 halos at every subdomain edge keep their initial value (the
`prof.inp` e12 profile, 0.2-1 m/s) for the whole run. The same code is on origin/master (61d81c83). Evidence:
- e12 budget at z = 2420 m (free troposphere), slab mean on rank 0: shear + buoyancy + dissipation = -4.4e-4 m/s2,
  balanced by +4.4e-4 from e12 advection (`advecc_2nd`), which for a quiet homogeneous layer should be about 0; e12
  stays at 0.15-0.16 (initial 0.19) instead of decaying to e12min as in DALES (1e-4).
- 6 h restart fields, slab means (fixed dt, seed 43): e12 0.61 vs DALES 0.24 at 300 m, 0.15 vs 0.002 at 2420 m;
  ekm 5.3 vs 2.0 m2/s at 300 m. Resolved fluxes are weak (cloud-layer w'qt' 0.007 vs 0.045 g/kg m/s in hours 3-6).
- The same binary without the fix on 1 rank (1 x 1, so the periodic copy is used): thl/qt at 300 m and 1 h = 298.717 K /
  16.723 g/kg, equal to the fixed binary on 4 x 2 (298.717 / 16.723) and to DALES (298.715 / 16.719); on 4 x 2 without
  the fix 298.814 / 16.602.
- Bisection (6 h BOMEX): the commit that aligned the one-equation SGS model with DALES (Kh = Km/Prandtl and no factor 2 in the
  dissipation) and every later commit was bad. It did not introduce the bug; it removed two compensating errors
  (Kh = Km/3 and doubled dissipation) that had been hiding it. `Prandtl = 3.` (Kh = Km/3) with the bug gives LWP 6.8 (compensation, not a fix).
- The fix (issue #405, with the regression test `TestSGSHaloExchange` and a CHANGELOG entry): in
  `halos`, `if (loneeqn) then; call exchange_halo_z(e120); call exchange_halo_z(e12m); end if`, with `e120, e12m`
  from modfields and `loneeqn` from modsubgriddata. It changes every one-equation run with more than one rank in x
  or y. Every earlier uDALES BOMEX number in this file, and the time-step tuning (`courant = 0.55`), were affected.

## Earlier results (superseded; all affected by the e12 halo bug)

Early run (Release, Intel/CX3, single node, 8 ranks), main loop 521 s, no NaN.

**Published-range checks are pending citable values**: `PUBLISHED` in compare_dales.py is empty by design, so there is no pass/fail table; numbers below are reported without verdicts. Interpretation is tentative.

### Reference: DALES dev, run like-for-like

The reference is DALES `dev` at 075ac4e, built on CX3 (Intel 2021a, Release) and run with `cases/bomex` for 6 h on 16 ranks. Only the runtime and the decomposition were changed, and `lsamptend` was switched on for budget output. That run applies the prescribed -2 K/day radiative cooling (`thltendradall` = -2.315e-5 K/s up to ~1500 m), as `examples/301` does via `thlpcar`.

Window 3-6 h:

| quantity | uDALES | DALES dev |
|---|---|---|
| LWP [g/m2] | 11.37 | 7.53 |
| cloud cover | 0.226 | 0.164 |
| zi [m] | 1679 | 1684 |
| u* [m/s] | 0.28 | 0.28 |

Max abs profile differences (3-6 h mean): thl 0.133 K, qt 0.240 g/kg, ql 0.0054 g/kg, resolved w'qt' 0.028 g/kg m/s.

Hourly LWP [g/m2], hours 0-1 to 5-6:
- uDALES: 1.24, 7.24, 12.63, 9.19, 11.64, 13.27
- DALES dev: 3.85, 4.48, 5.97, 6.84, 8.98, 6.78

Observations (tentative):
- The mean state agrees: zi matches to 5 m and thl to 0.13 K. u* is prescribed (0.28) in both.
- uDALES has about 1.5x the LWP and 1.4x the cloud cover of DALES dev. uDALES LWP is noisy hour to hour, so the 3-h means carry a large sampling error.
- DALES has cloud within the first hour; uDALES needs about two hours. uDALES applies no initial thl/qt noise because `randthl`/`randqt` are never applied (#395).
- Candidates for the residual (the A/B below tests (a) and (d)):
  - (a) Doubled SGS TKE dissipation: `src/modsubgrid.f90:495` has `sbdiss = - 2. * (ce1 + ce2*zlt/delta) ...` ("added factor 2 for shear-driven flow"); DALES dev has no factor 2.
  - (b) No TKE sources at the first level: the `sources` loop starts at `kb+1` (`src/modsubgrid.f90:454`), so `e12p(kb)` gets no shear, buoyancy or dissipation source (DALES has a dedicated k = 1 block).
  - (c) `grav/thvs` instead of DALES `grav/thvf(k)` in `zlt` (`src/modsubgrid.f90:392`) and `sbbuo` (`src/modsubgrid.f90:492`).
  - (d) No initial thl/qt noise (#395). Hours 3-6 are still 1.4-1.9x, which is not a spin-up signature, so this alone is unlikely to explain it.
  - (e) The one-step linearised saturation adjustment (`thermo` in `src/modthermodynamics.f90`) overestimates `ql` by about 1-3 %. It has the same sign but is too small, so it is a quantified non-candidate.

### A/B test of (a) and (d)

Throwaway builds of the code at that time, same case and ranks, one 6 h realisation each. Variant A removes the factor `2.` in `sbdiss`. Variant B applies thl/qt initial noise (`randthl = 0.1`, `randqt = 2.5e-5`, the DALES BOMEX values, up to `krand`). Window 3-6 h. The last two columns are the max abs difference from DALES dev over height.

| run | LWP [g/m2] | cloud cover | zi [m] | thl [K] | qt [g/kg] |
|---|---|---|---|---|---|
| uDALES as is | 11.37 | 0.226 | 1679 | 0.13 | 0.24 |
| A: no factor 2 | 7.40 | 0.153 | 1717 | 0.24 | 0.36 |
| B: thl/qt noise | 12.52 | 0.212 | 1696 | 0.13 | 0.20 |
| A+B | 8.50 | 0.172 | 1693 | 0.25 | 0.39 |
| DALES dev | 7.53 | 0.164 | 1684 | - | - |

The doubled dissipation (a) accounts for the cloud excess. Noise (d) makes no significant difference: hourly LWP varies by up to 2x, so about 1 g/m2 is within sampling noise. With A the mean profiles move further from DALES and zi rises by about 35 m, so other SGS differences remain. Reported as #401. These runs predate the SGS alignment (ch = 1/Prandtl, no factor 2, DALES kb sources); results with the aligned model replace them in the close-match validation.

### The DALES 3.2.0 archive is not a like-for-like reference

`cases/bomex/results_3.2.0` in the DALES repository has `tltendrad` = 0, and its column thl budget closes (residual 3e-5 K m/s) only without the -2 K/day cooling. The uDALES budget closes only with it (residual 4e-4). The qt budgets close in both, so the surface-flux sign and magnitude are right. This was reported upstream as dalesteam/dales#348.

For context, the archive compared with uDALES run with and without the cooling (`thlpcar = 0` in `lscale.inp.301`), 3-6 h:

| quantity | uDALES with cooling | uDALES `thlpcar = 0` | DALES 3.2.0 |
|---|---|---|---|
| LWP [g/m2] | 11.37 | 7.28 | 4.34 |
| cloud cover | 0.226 | 0.162 | 0.138 |
| zi [m] | 1679 | 1616 | 1620 |

The no-cooling pair shows the same pattern as the dev pair: matching zi, thl within 0.11 K, and about 1.7x the LWP.
