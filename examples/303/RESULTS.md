# Tier-2 validation results: 303 DYCOMS-II RF01

## Facet floor with the Monin-Obukhov facet surface and the first-level SGS treatment (shipped configuration)

The shipped `namoptions.303`: a floor made of one horizontal facet with `iwallmom = 4` (local Monin-Obukhov per facet
section, u* from MO, fluxes `bctfz`, `bcqfz`), kappa = 0.40, `z0 = z0h = 1.445944e-5 m`, and in the 9216 first-level
cells above the facet the DALES first-level treatment of the one-equation SGS model (TKE shear and buoyancy sources and
`dthvdz` from the Monin-Obukhov surface-layer gradients, evaluated before `subgrid`). Runs: Release build, FFTW planner
`FFTW_ESTIMATE`, fixed dt = 0.3 s, 32 ranks (8 x 4), seeds 43-45, 4 h. DALES: the runs of the sections below. Window
3-4 h, 600 s blocks; entrainment fit 1-4 h.

**Verdict.** The near-floor difference of the section below is gone: at 15 m the resolved w'thl' is 5.34e-3 against
5.36e-3 K m/s in DALES and w'qt' 1.78e-5 against 1.80e-5 kg/kg m/s (within 1 %; without the treatment 7.86e-3 and
2.56e-5). All five profiles meet the bar (rms difference 0.19-0.69 of the smaller rms spread; resolved fluxes were 2.64
and 1.58). LWP, zi and u* agree within the spread (u* 0.2466 against 0.2469 m/s, within the larger spread). zi_qt8 (+0.8 m)
and the entrainment rate (4.82 against 4.75 mm/s, +1.5 %) are within the larger spread, the uDALES spread being the
larger. Cloud cover is now 1.2 % above DALES (0.818 against 0.808, 2.7 of the smaller and 1.0 of the larger spread);
without the treatment it was 1.2 % below (0.799), so the treatment moves it by +0.019, more than either spread.

| quantity | uDALES members | uDALES mean (spread) | DALES members | DALES mean (spread) | diff | ratio (min) | ratio (max) |
|---|---|---|---|---|---|---|---|
| LWP [g/m2] | 15.75, 16.28, 14.71 | 15.58 (0.80) | 17.19, 12.99, 15.44 | 15.20 (2.11) | +0.38 | 0.47 | 0.18 |
| cloud cover | 0.8099, 0.8288, 0.8164 | 0.8184 (0.0096) | 0.8116, 0.8043, 0.8093 | 0.8084 (0.0038) | +0.0100 | 2.66 | 1.04 |
| u* [m/s] | 0.2463, 0.2470, 0.2464 | 0.2466 (0.0004) | 0.2467, 0.2472, 0.2468 | 0.2469 (0.0003) | -0.0003 | 1.15 | 0.78 |
| zi, max dthl/dz [m] | 875.0, 875.0, 874.2 | 874.7 (0.5) | 875.8, 874.2, 874.2 | 874.7 (1.0) | 0.0 | 0.00 | 0.00 |
| zi_qt8 [m] | 859.2, 859.8, 858.2 | 859.1 (0.8) | 858.7, 857.8, 858.4 | 858.3 (0.5) | +0.8 | 1.61 | 0.94 |
| we [mm/s] (fit 1-4 h) | 4.858, 4.923, 4.680 | 4.821 (0.126) | 4.761, 4.726, 4.760 | 4.749 (0.020) | +0.072 | 3.57 | 0.57 |

| profile (z <= 1200 m) | max abs diff (z) | rms diff | rms spread U / D | rms diff / min rms spread | levels with diff < both spreads |
|---|---|---|---|---|---|
| thl [K] | 0.146 (868 m) | 0.015 | 0.034 / 0.031 | 0.50 | 82 % |
| qt [g/kg] | 0.102 (868 m) | 0.013 | 0.025 / 0.026 | 0.52 | 78 % |
| ql [g/kg] | 0.0054 (848 m) | 0.0011 | 0.0021 / 0.0042 | 0.54 | 35 % |
| w'thl' resolved [K m/s] | 0.0014 (865 m) | 0.0003 | 0.0004 / 0.0005 | 0.69 | 73 % |
| w'qt' resolved [g/kg m/s] | 0.0009 (865 m) | 0.0003 | 0.0017 / 0.0017 | 0.19 | 89 % |

Near-floor resolved fluxes, window 3-4 h, ensemble mean (spread), zm = 5, 10, 15 m:

| run | w'thl' [K m/s] | w'qt' [kg/kg m/s] |
|---|---|---|
| facet floor with the first-level treatment | 1.74e-3 (1e-5), 3.54e-3 (2e-5), 5.34e-3 (2e-5) | 5.81e-6 (4e-8), 1.18e-5 (7e-8), 1.78e-5 (8e-8) |
| facet floor without it (section below) | 3.12e-3 (1e-5), 5.93e-3 (2e-5), 7.86e-3 (2e-5) | 1.01e-5 (7e-8), 1.92e-5 (1.3e-7), 2.56e-5 (1.5e-7) |
| DALES | 1.75e-3 (1e-5), 3.57e-3 (2e-5), 5.36e-3 (3e-5) | 5.88e-6 (2e-8), 1.20e-5 (4e-8), 1.80e-5 (4e-8) |

Hourly LWP [g/m2], ensemble mean (spread):

| code | 0-1 h | 1-2 h | 2-3 h | 3-4 h |
|---|---|---|---|---|
| uDALES, first-level treatment | 27.38 (0.19) | 14.02 (0.26) | 15.72 (0.55) | 15.58 (0.80) |
| DALES | 26.85 (0.08) | 14.60 (0.23) | 15.66 (0.29) | 15.20 (2.11) |

Against the facet-floor runs without the treatment (section below; 60 s blocks): LWP -0.19 g/m2 (0.24 spreads), u*
+0.0001 m/s (0.49), cloud cover +0.019 (3.5 / 2.0 spreads), zi +0.6 m, zi_qt8 +1.0 m (2.0 / 1.3), we +0.09 mm/s
(1.8 / 0.7); thl, qt, ql profiles 0.42-0.76 of the smaller rms spread; resolved w'thl', w'qt' 2.17 and 1.57 (the
near-floor change).

## Facet floor without the first-level SGS treatment (earlier runs)

These runs predate the first-level SGS treatment above horizontal Monin-Obukhov facets; otherwise the configuration
is the shipped one.

The floor is one facet with `iwallmom = 4` (local Monin-Obukhov per facet section, u* from MO, fluxes `bctfz`,
`bcqfz`), without TKE sources at `kb`, kappa = 0.40, `z0 = z0h = 1.445944e-5 m` in the floor facet type. Runs: Release
build, FFTW planner `FFTW_ESTIMATE`, fixed dt = 0.3 s, 32 ranks (8 x 4), seeds 43-45, 4 h. DALES: the runs of the next
section (DALES dev 075ac4e with both patches, `lmostlocal = .true.` by default). Window 3-4 h, 600 s blocks; entrainment
fit 1-4 h.

**Verdict.** LWP, zi, zi_qt8 and the thl, qt, ql profiles agree within the seed spread (profile rms ratios 0.47-0.79).
u* is 0.2464 against 0.2469 m/s (0.2 %, 1.7 of the tiny spreads, as with the earlier flat-bottom runs) and the entrainment
rate 4.71 against 4.75 mm/s (0.8 %, within the larger spread). Cloud cover is 0.799 against 0.808 (-1.2 %, 2.5 spreads).
The resolved fluxes differ at the lowest levels (largest difference at 15 m), the same signature as in BOMEX
(examples/301/RESULTS.md): DALES has SGS TKE sources at its first level from the MO surface-layer shears, uDALES has none
at `kb`, so more of the surface flux is resolved near the floor.

| quantity | uDALES members | uDALES mean (spread) | DALES members | DALES mean (spread) | diff | ratio (min) | ratio (max) |
|---|---|---|---|---|---|---|---|
| LWP [g/m2] | 15.36, 15.28, 16.66 | 15.77 (0.78) | 17.19, 12.99, 15.44 | 15.20 (2.11) | +0.56 | 0.73 | 0.27 |
| cloud cover | 0.7928, 0.8017, 0.8028 | 0.7991 (0.0055) | 0.8116, 0.8043, 0.8093 | 0.8084 (0.0038) | -0.0093 | 2.49 | 1.70 |
| u* [m/s] | 0.2468, 0.2462, 0.2464 | 0.2464 (0.0003) | 0.2467, 0.2472, 0.2468 | 0.2469 (0.0003) | -0.0005 | 1.66 | 1.59 |
| zi, max dthl/dz [m] | 874.2, 874.2, 874.2 | 874.2 (0) | 875.8, 874.2, 874.2 | 874.7 (1.0) | -0.6 | - | 0.58 |
| zi_qt8 [m] | 858.4, 857.4, 858.0 | 858.0 (0.5) | 858.7, 857.8, 858.4 | 858.3 (0.5) | -0.3 | 0.72 | 0.67 |
| we [mm/s] (fit 1-4 h) | 4.750, 4.669, 4.707 | 4.709 (0.041) | 4.761, 4.726, 4.760 | 4.749 (0.020) | -0.040 | 1.98 | 0.98 |

| profile (z <= 1200 m) | max abs diff (z) | rms diff | rms spread U / D | rms diff / min rms spread | levels with diff < both spreads |
|---|---|---|---|---|---|
| thl [K] | 0.109 (872 m) | 0.015 | 0.034 / 0.031 | 0.49 | 34 % |
| qt [g/kg] | 0.058 (872 m) | 0.011 | 0.025 / 0.026 | 0.47 | 62 % |
| ql [g/kg] | 0.0056 (822 m) | 0.0012 | 0.0015 / 0.0042 | 0.79 | 19 % |
| w'thl' resolved [K m/s] | 0.0025 (15 m) | 0.0009 | 0.0003 / 0.0005 | 2.64 | 33 % |
| w'qt' resolved [g/kg m/s] | 0.0077 (15 m) | 0.0016 | 0.0010 / 0.0017 | 1.58 | 49 % |

Hourly LWP [g/m2], ensemble mean (spread):

| code | 0-1 h | 1-2 h | 2-3 h | 3-4 h |
|---|---|---|---|---|
| uDALES facet floor | 27.28 (0.02) | 13.62 (0.20) | 14.36 (0.95) | 15.77 (0.77) |
| DALES | 26.85 (0.08) | 14.60 (0.23) | 15.66 (0.29) | 15.20 (2.11) |

Secondary reference, the earlier flat-bottom runs of the next section (`BCbotm = 4`, removed since; kappa 0.40): LWP
+0.32 g/m2 (0.42 spreads), u* -0.0001 m/s (0.35), zi_qt8 -0.2 m (0.47), we -0.02 mm/s (0.39), cloud cover -0.012 (2.3, 1.1
of the larger spread); thl, qt, ql profiles within the spread (0.49-0.97), resolved fluxes again differ near the floor.

## Earlier close-match validation: flat-bottom Monin-Obukhov route (`BCbotm = 4`, since removed)

**Verdict: very close.** With the configuration that shipped then (uDALES Release build of this change, `examples/303/namoptions.303`
with the fixed-dt block enabled, kappa = 0.40, `z0 = z0h = 1.445944e-5 m`; DALES dev 075ac4e with the two patches below and
its native kappa = 0.4, inputs `tools/python/examples/dales_cases/dales_inputs/rf01`), three seeds each over 4 h,
uDALES and DALES agree within the seed-to-seed spread for LWP, cloud cover, zi, zi_qt8 and (within 0.3 %) u* and the
entrainment rate. For all five profiles the rms difference is below the smaller rms spread (ratio 0.13-0.51). Two items sit
just above the strict bar (the difference exceeds the smaller of the two spreads) with tiny absolute differences: u*
(ratio 1.29: 0.2465 vs 0.2469 m/s, 0.14 %, because the DALES spread is only 0.0003) and the entrainment rate (ratio 1.27:
4.72 vs 4.75 mm/s, 0.5 %). Level by level the fraction of levels with |diff| below both spreads is 44 % (thl), 72 % (qt),
34 % (ql), 83 % (w'thl'), 86 % (w'qt').

Settings (identical in both codes unless stated): 96 x 96 x 300 cells (35 m x 35 m x 5 m), fixed dt = 0.3 s, cd2 advection
for all fields, one-equation SGS model with the DALES constants, `radpar` longwave radiation every step plus the
free-troposphere term (uDALES `iradiation = 2` with the `dycoms_rf01` case hook; DALES `iradiation = 10` with a patched
`rad_user`), Monin-Obukhov surface with prescribed fluxes, `igrw_damp = 2` sponge from 1122.5 m, thl and qt noise
(+-0.1 K, +-2.5e-5 kg/kg, all levels, no velocity noise; the random generators differ), 8 x 4 ranks, seeds (`irandom`) 43, 44, 45,
window 3-4 h. DALES needs `dales_inputs/dales-075ac4e-user-es.patch` (flux-form free-troposphere term, Tetens saturation
vapour pressure with the uDALES constants) and `dales_inputs/dales-075ac4e-iradiation10.patch` (without it DALES switches
`radpar` off for `iradiation = 10`, see below); no `fkar` edit is needed. Not matched: DALES uses a per-column Obukhov length
(`lmostlocal = .true.`, its dev default), uDALES one domain-mean value; u* agrees to 0.0004 m/s. Analysis:
`compare_dales.py 303 <out> --udales ... --dales ...`.

| quantity | uDALES members | uDALES mean (spread) | DALES members | DALES mean (spread) | diff | ratio (min / max spread) |
|---|---|---|---|---|---|---|
| LWP [g/m2] | 13.9, 17.8, 14.6 | 15.4 (2.1) | 17.2, 13.0, 15.4 | 15.2 (2.1) | +0.24 | 0.12 / 0.11 |
| cloud cover | 0.812, 0.800, 0.823 | 0.8115 (0.011) | 0.812, 0.804, 0.809 | 0.8084 (0.004) | +0.003 | 0.82 / 0.28 |
| u* [m/s] | 0.2479, 0.2455, 0.2463 | 0.2465 (0.0013) | 0.2467, 0.2472, 0.2468 | 0.2469 (0.0003) | -0.0004 | 1.29 / 0.28 |
| zi, max dthl/dz [m] | 874.2, 874.2, 875.0 | 874.4 (0.5) | 875.8, 874.2, 874.2 | 874.7 (1.0) | -0.3 | 0.58 / 0.29 |
| zi_qt8 [m] | 857.8, 857.9, 858.7 | 858.2 (0.5) | 858.7, 857.8, 858.4 | 858.3 (0.5) | -0.15 | 0.31 / 0.31 |
| we = dzi_qt8/dt + D zi, fit 1-4 h [mm/s] | 4.67, 4.74, 4.76 | 4.72 (0.05) | 4.76, 4.73, 4.76 | 4.75 (0.02) | -0.03 | 1.27 / 0.54 |

| profile (z <= 1200 m) | max abs diff (z) | rms diff | rms spread U / D | rms diff / min rms spread | levels with diff < both spreads |
|---|---|---|---|---|---|
| thl [K] | 0.042 (862 m) | 0.0059 | 0.016 / 0.031 | 0.38 | 44 % |
| qt [g/kg] | 0.020 (862 m) | 0.0055 | 0.014 / 0.026 | 0.39 | 72 % |
| ql [g/kg] | 0.0021 (828 m) | 0.0005 | 0.0041 / 0.0042 | 0.13 | 34 % |
| w'thl' resolved [K m/s] | 0.0012 (865 m) | 0.0003 | 0.0008 / 0.0005 | 0.51 | 83 % |
| w'qt' resolved [g/kg m/s] | 0.0007 (420 m) | 0.0004 | 0.0021 / 0.0017 | 0.25 | 86 % |

Hourly LWP [g/m2], ensemble mean (spread): uDALES 27.4 (0.2), 14.0 (0.5), 15.5 (0.3), 15.4 (2.1); DALES 26.9 (0.1),
14.6 (0.2), 15.7 (0.3), 15.2 (2.1). The 3 members per code give a noisy spread estimate; more seeds would tighten the
two ratios just above 1.

### Earlier record: kappa = 0.41 in both codes, `z0 = 1.07e-5 m`

An earlier three-seed pair of the same comparison used kappa = 0.41 in both codes (DALES built with `fkar = 0.41`) and
`z0 = z0h = 1.07e-5 m` (the same neutral drag coefficient Cm = 0.0011 as the shipped values). It gave the same picture:
LWP 15.8 vs 16.2 g/m2, cloud cover 0.813 vs 0.815, u* 0.2466 vs 0.2466 m/s, zi_qt8 858.8 vs 858.8 m, we 4.80 vs 4.81 mm/s.
The 20-40 min window was above the strict bar (LWP 17.0 vs 15.8 g/m2, we 9.6 vs 10.0 mm/s; the spreads are tiny that
early) and closed by 1 h; it was not re-measured for the shipped configuration.

### The DALES reference defect (iradiation = 10)

An earlier mismatch (uDALES entraining about twice as fast, LWP 15.8 vs 9.4 g/m2) was a defect in the DALES reference, not in
uDALES: DALES dev 075ac4e `modradiation.f90:initradiation` sets `rad_longw = rad_shortw = rad_smoke = .false.` for
`iradiation = 10`, so the `radpar` call in its `irad_user` branch never runs and those DALES runs had no cloud-top longwave
cooling (F0 = 70, F1 = 22 W/m2), only the free-troposphere term. Found with a horizontally homogeneous single-column run
of both codes (all noise off): every process tendency agreed to 3-4 digits except radiation, which was exactly zero in the
DALES cloud. Reported upstream as a comment on dalesteam/dales#347. The 2.3 % u* gap of that comparison was the von Karman
constant (DALES 0.4, uDALES 0.41), which is why the uDALES default `fkar` is now 0.40. A 2 x 2 against 8 x 4 uDALES
check showed no decomposition dependence.

## Earlier close-match comparison: superseded (DALES reference without cloud-top cooling)

The first close-match comparison (fixed dt = 0.3 s, three seeds, window 3-4 h, uDALES with the e12 halo fix, DALES with
the Tetens patch but without the `iradiation = 10` patch and with kappa = 0.4) found uDALES entraining about twice as fast:
LWP 15.8 vs 9.4 g/m2, cloud cover 0.81 vs 0.63, zi_qt8 859 vs 834 m, we 4.8 vs 2.6 mm/s, u* +2.3 %. The cause was the DALES
reference (no cloud-top longwave cooling) and the kappa difference, see above; the conclusions drawn then about DALES
behaviour, cloud loss and the role of cd2 are withdrawn, and the adaptive-dt notes and the cd2 against 5th-order advection
comparison were made with that defective reference and have not been repeated (open item, see
`docs/udales-dales-numerics.md`).

Two DALES input problems found on the way are fixed in `dales_inputs/rf01`: w noise at the surface level (`krandumax > 0`
with the default `krandumin = 1` leaves a fixed random vertical velocity, `divmax = rhobf(1) randu / dzf`, at the ground for the whole
run), and the sponge (`igrw_damp = 2` is the DALES default and is now set explicitly).

## Earlier results (superseded; e12 halo bug and older SGS)

Early run (Release, Intel/CX3, single node, 32 ranks), main loop 5704 s, no NaN.

**Published-range checks are pending citable values**: `PUBLISHED` in compare_dales.py is empty by design, so there is no pass/fail table; numbers below are reported without verdicts. Interpretation is tentative.

No DALES reference. Time series (interpolated; window 3-4 h means from compare_dales):

| quantity | 0.5 h | 1 h | 2 h | 3 h | 4 h | 3-4 h mean |
|---|---|---|---|---|---|---|
| LWP [g/m2] | 26.6 | 42.7 | 54.8 | ~60 | ~63 | 59.8 |
| cloud cover | 1.000 | 0.9967 | 0.9978 | 0.9986 | 0.9997 | 0.9995 |
| zi (max dthl/dz) [m] | 850.6 | 858.0 | 869.4 | 879.7 | 889.5 | 884.4 |
| zi_qt8 [m] (now `zi_qt`, `zi_qt_thres = 8e-3`) | 845.4 | 845.5 | 856.9 | 868.5 | 879.3 | - |
| u* [m/s] | 0.255 | 0.245 | 0.253 | 0.250 | 0.248 | 0.248 |

Entrainment rate (zi_qt8, large-scale divergence D = 3.75e-6 1/s, fit 1-4 h): we = 6.26e-3 m/s. (LWP at 0.5-2 h computed from bldump_t; 3 h and 4 h are read from the plot; the 3-4 h mean is from summary.md.)

Observations (tentative):
- No cloud collapse: cloud cover stays above 0.994 throughout and LWP, after an initial spin-up dip to ~19 g/m2 at 0.8 h, recovers and rises slowly. Stratocumulus-top ql peaks at ~0.43 g/kg near 850 m, inversion near 870-890 m.
- zi rises steadily, about 39 m from 0.5 to 4 h; the net `we` is the sum of that and the subsidence term D*zi (about 3.3e-3 m/s), so about half of the fitted `we` is the prescribed divergence, not a sign of an anomaly. No published-range comparison made (no citable values).
