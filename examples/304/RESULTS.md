# Tier-2 validation results: 304 dry CBL over a floor of fixed temperature

Validated binary: commit f19067f8 (Release, Intel 2021a, FFTW planner `FFTW_ESTIMATE`). The later commits on the
branch change only docs, examples and tooling. DALES: dev 075ac4e unpatched (Intel 2021a build), inputs
`tools/python/examples/dales_cases/dales_inputs/cbl304` (`isurf = 2`, `lmostlocal` default). Both codes: the shipped
namoptions, fixed dt = 1.0 s, 16 ranks (4 x 4), seeds 43-45, 4 h. Window 2-4 h, 600 s profile blocks
(`compare_dales.py 304`). Run directories (CX3): `/rds/general/user/mvr/ephemeral/t24/{u304,d304}/f4{3,4,5}`,
analysis `t24/an304`.

## Verdict

The surface scheme reproduces DALES `isurf = 2`: per unit first-level wind speed, u* agrees within 0.03 % and the
surface heat flux within 0.1 % (window means: u*/|U1| = 0.09166 against 0.09163, wtheta/|U1| = 7.330e-3 against
7.321e-3 m/s K per m/s). The absolute u* and surface heat flux are 0.35 % and 0.28 % below DALES, 1.0 and 1.1 of the
larger ensemble spread (1.25 and 1.21 of the smaller): the first-level wind speed |U1| in the window is 0.4 % lower in
uDALES (3.060 against 3.072 m/s; u1 lower, v1 higher, a slightly different evolution of the mean wind under Coriolis),
and u* and the flux scale with it. zi agrees to within a grid level (+1.4 m, dz = 25 m; the DALES members are
identical, so its spread is 0). The mean thl profile is within the spread (rms difference 0.74 of the larger and 1.06
of the smaller rms spread; largest difference 0.003 K). The resolved w'thl' profile is 0.6-1 % below DALES through the
mixed layer (the surface flux is 0.3 % lower) and 3.3 % and 1.8 % lower at 25 and 50 m, where the resolved/SGS split
next to the floor matters (DALES SGS flux at 25 m 1.46e-2 of a total 2.16e-2 K m/s); the ensemble spreads are only
0.3-0.7 % of the flux, so this is 2.0-2.3 rms spreads. That near-floor signature is the one of 301/303 (resolved w'thl'
1.9 % low at 40 m in BOMEX), not specific to the facet-temperature path. qt, ql and w'qt' are zero in uDALES (dry); DALES
has 1e-10 qt from its dry-surface setting (see README): ignore those rows.

| quantity | uDALES members | uDALES mean (spread) | DALES members | DALES mean (spread) | diff | ratio (min) | ratio (max) |
|---|---|---|---|---|---|---|---|
| u* [m/s] | 0.2814, 0.2794, 0.2805 | 0.2805 (0.0010) | 0.2822, 0.2816, 0.2806 | 0.2815 (0.0008) | -0.0010 | 1.25 | 0.99 |
| surface heat flux wtheta [K m/s] | 0.02248, 0.02238, 0.02241 | 0.02243 (5.1e-5) | 0.02254, 0.02251, 0.02242 | 0.02249 (5.9e-5) | -6.2e-5 | 1.21 | 1.06 |
| zi, max dthl/dz [m] | 737.5, 735.4, 737.5 | 736.8 (1.2) | 735.4, 735.4, 735.4 | 735.4 (0) | +1.4 | inf | 1.15 |
| first-level |U1| [m/s] | | 3.0602 | | 3.0721 | -0.0119 (-0.4 %) | | |
| u* / |U1| | | 0.09166 | | 0.09163 | +0.03 % | | |
| wtheta / |U1| [K] | | 7.330e-3 | | 7.321e-3 | +0.1 % | | |

| profile (z <= 1450 m) | max abs diff (z) | rms diff | rms spread U / D | rms diff / min rms spread | rms diff / max rms spread | levels with diff < both spreads |
|---|---|---|---|---|---|---|
| thl [K] | 0.0030 (162 m) | 0.0013 | 0.0018 / 0.0013 | 1.06 | 0.74 | 47 % |
| w'thl' resolved [K m/s] | 2.35e-4 (50 m) | 8.7e-5 | 3.9e-5 / 4.5e-5 | 2.25 | 1.96 | 29 % |

Near-floor resolved w'thl' [K m/s], window 2-4 h, ensemble mean (spread):

| zm [m] | uDALES | DALES | difference |
|---|---|---|---|
| 25 | 6.752e-3 (2.7e-5) | 6.981e-3 (4.7e-5) | -3.3 % |
| 50 | 1.316e-2 (1.6e-5) | 1.339e-2 (6.8e-5) | -1.8 % |
| 75 | 1.592e-2 (2.6e-5) | 1.607e-2 (7.3e-5) | -0.9 % |
| 100 | 1.643e-2 (4.0e-5) | 1.654e-2 (7.4e-5) | -0.7 % |
| 200 | 1.396e-2 (5.5e-5) | 1.406e-2 (5.8e-5) | -0.7 % |

Evolution of the differences (ensemble means, uDALES / DALES - 1): surface heat flux +0.16 % in 0-10 min, -0.50 % in
10-30 min, -0.08 % in 30-60 min, -0.51 % in 1-2 h, -0.28 % in 2-4 h; u* +0.11 %, -0.29 %, -0.06 %, -0.61 %, -0.35 %. The
first-level thl agrees to 0.0005 K in the window, so the flux difference is not a temperature difference.

The bldump `obl` (`1/mean(1/L)` over the facet sections, -71.9 m) is not comparable with the DALES tmser `obukh`
(-74.9 m), which DALES computes from the slab-mean state (`oblav`), so it is not in the bar.

## Bit identity of the iwalltemp = 1 cases

Examples 301 and 303 (both `iwallmom = 4`, `iwalltemp = 1`, with first-level cells), Release builds with
`FFTW_ESTIMATE` of 3a428b09 (the parent) and f19067f8, 8 ranks (4 x 2), fixed dt (301: 100 steps of 1 s; 303: 100
steps of 0.3 s): every common variable of `bldump_t` and `bldump_z` identical (numpy `array_equal`, 10 records) and
all 8 restart files byte-identical in both cases. `bldump_t` gains `wtheta` and `wq`.
