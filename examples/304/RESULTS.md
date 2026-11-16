# Tier-2 validation results: 304 dry CBL over a floor of fixed temperature

Validated binary: commit f19067f8 (Release, Intel 2021a, FFTW planner `FFTW_ESTIMATE`). The later source changes
(review fixes) do not change the solution: one `Tair` interpolation is no longer repeated (same value), and `fac_htc`
(output only) has a new wind floor. DALES: dev 075ac4e unpatched (Intel 2021a build), inputs
`tools/python/examples/dales_cases/dales_inputs/cbl304` (`isurf = 2`, `lmostlocal` default). Both codes: the shipped
namoptions, fixed dt = 1.0 s, 16 ranks (4 x 4), 4 h, seeds 43-48 (first 43-45, then 46-48). Window 2-4 h, 600 s
profile blocks (`compare_dales.py 304`). Run directories (CX3): `/rds/general/user/mvr/ephemeral/t24/{u304,d304}/f4{3..8}`;
analyses `t24/an304` (3 seeds), `t24/an304_6` (6 seeds), per-run normalised values `t24/norm304.py`.

Ratios are |mean difference| / spread, with the smaller (min) and the larger (max) of the two ensemble spreads
(sample standard deviation of the members' window means). The bar is ratio (min) < 1.

## Verdict (6 seeds)

| quantity | uDALES mean (spread) | DALES mean (spread) | diff | ratio min / max | bar |
|---|---|---|---|---|---|
| u* [m/s] | 0.2807 (0.0009) | 0.2815 (0.0005) | -0.27 % | 1.42 / 0.85 | missed (within the larger spread) |
| surface heat flux wtheta [K m/s] | 0.02245 (6.1e-5) | 0.02249 (4.0e-5) | -0.16 % | 0.92 / 0.60 | met |
| zi, max dthl/dz [m] | 736.1 (1.1) | 736.5 (1.7) | -0.3 m | 0.32 / 0.20 | met |
| thl profile (rms diff / rms spread) | | | | 0.81 / 0.57 | met |
| resolved w'thl' profile (rms diff / rms spread) | | | | 0.54 / 0.44 | met |
| resolved w'thl' at 25 m | 6.839e-3 (1.0e-4) | 6.970e-3 (6.2e-5) | -1.9 % | 2.11 / 1.30 | missed |
| resolved w'thl' at 50 m | 1.327e-2 (1.2e-4) | 1.339e-2 (7.3e-5) | -0.9 % | 1.66 / 0.98 | missed (within the larger spread) |
| first-level wind |U1| [m/s] | 3.0639 (0.0100) | 3.0731 (0.0076) | -0.30 % | 1.21 / 0.88 | (diagnostic) |
| u* / |U1| | 0.091628 (5.1e-5) | 0.091603 (7.2e-5) | +0.03 % | 0.50 / 0.35 | met |
| wtheta / |U1| [K] | 7.3281e-3 (1.2e-5) | 7.3180e-3 (6.4e-6) | +0.14 % | 1.57 / 0.81 | missed (within the larger spread) |

With 6 seeds the surface heat flux, zi, the thl profile and the resolved w'thl' profile meet the bar. u* misses it
(1.42 of the smaller spread, 0.85 of the larger). The resolved w'thl' at the two lowest half levels misses it. The 3-seed
numbers (below) were further off, and the spreads themselves changed by up to a factor 2 between 3 and 6 seeds.

What the runs show about the cause:
- The surface transfer per unit first-level wind agrees: u*/|U1| within 0.03 % (0.50 / 0.35 spreads), wtheta/|U1| within
  0.14 % (within the larger spread). The first-level wind is 0.30 % lower in uDALES, and u* and the flux follow it.
- **Control with a prescribed flux** (below: uDALES `iwalltemp = 1`, DALES `isurf = 4`, same flux in both): the first-level
  wind is 0.48 % lower in uDALES and u* 0.43 % lower, with u*/|U1| again within 0.05 %. So the lower wind and u* do not
  come from the facet-temperature path. Their cause is not identified.
- One candidate is the u* stress stencil: uDALES applies `Cm(L_u) max(|U_u|, 0.1)^2` at the u point
  (src/modibm.f90:1886-1891, and the iwalltemp = 1 branch below it), DALES the square of the centre u* averaged to the u
  point (DALES modsubgrid.f90:942-951). Both stencils were applied to the same saved uDALES first-level fields (seed 43,
  8 snapshots 0.5-4 h; `t24/stencil304.py`, run `t24/u304fd`). The uDALES stencil gives a 0.04-0.06 % larger mean stress
  than the DALES stencil. That has the sign of the lower wind, but it is about an order of magnitude smaller than the
  wind difference, so it does not account for it by itself on this evidence.
- The near-floor deficit of the resolved w'thl' (-1.9 % at 25 m, -0.9 % at 50 m) does **not** appear in the control (-0.4 %
  and -0.1 %, within the spreads). It is not explained by the surface heat flux (-0.16 %). It is therefore associated with
  the facet-temperature runs. Its cause is not identified. One difference that exists only there is the local
  heterogeneity of the surface flux, which follows the local first-level wind and temperature in both codes.

### 3 seeds (43-45), the first analysis

| quantity | uDALES mean (spread) | DALES mean (spread) | diff | ratio min / max |
|---|---|---|---|---|
| u* [m/s] | 0.2805 (0.0010) | 0.2815 (0.0008) | -0.35 % | 1.25 / 0.99 |
| surface heat flux wtheta [K m/s] | 0.02243 (5.1e-5) | 0.02249 (5.9e-5) | -0.28 % | 1.21 / 1.06 |
| zi [m] | 736.8 (1.2) | 735.4 (0) | +1.4 m | inf / 1.15 |
| thl profile (rms) | | | | 1.06 / 0.74 |
| resolved w'thl' profile (rms) | | | | 2.25 / 1.96 |
| resolved w'thl' at 25 m | 6.752e-3 (2.7e-5) | 6.981e-3 (4.7e-5) | -3.3 % | 8.4 / 4.9 |
| first-level wind |U1| [m/s] | 3.0602 (0.012) | 3.0722 (0.011) | -0.39 % | 1.10 / 0.97 |
| u* / |U1| | 0.091650 (4.3e-5) | 0.091616 (6.7e-5) | +0.04 % | 0.78 / 0.51 |
| wtheta / |U1| [K] | 7.3284e-3 (1.4e-5) | 7.3201e-3 (7.4e-6) | +0.11 % | 1.12 / 0.60 |

### Per run (window 2-4 h)

| seed | uDALES u* | uDALES wtheta | uDALES |U1| | uDALES u*/|U1| | DALES u* | DALES wtheta | DALES |U1| | DALES u*/|U1| |
|---|---|---|---|---|---|---|---|---|
| 43 | 0.28143 | 0.022482 | 3.0723 | 0.091603 | 0.28219 | 0.022536 | 3.0821 | 0.091559 |
| 44 | 0.27943 | 0.022382 | 3.0476 | 0.091687 | 0.28157 | 0.022506 | 3.0739 | 0.091601 |
| 45 | 0.28054 | 0.022414 | 3.0607 | 0.091660 | 0.28061 | 0.022423 | 3.0605 | 0.091689 |
| 46 | 0.28004 | 0.022428 | 3.0575 | 0.091593 | 0.28119 | 0.022463 | 3.0699 | 0.091596 |
| 47 | 0.28184 | 0.022556 | 3.0747 | 0.091665 | 0.28170 | 0.022497 | 3.0729 | 0.091673 |
| 48 | 0.28111 | 0.022450 | 3.0703 | 0.091557 | 0.28174 | 0.022507 | 3.0792 | 0.091497 |

|U1| = |(window mean u1, window mean v1)| at the first level: uDALES bldump `u1`, `v1` (60 s samples), DALES `u`, `v` of
`profiles.001.nc` at k = 1 (600 s blocks in the window); u* and wtheta are window means of the time series.

### Near-floor resolved w'thl' [K m/s], 6 seeds

| zm [m] | uDALES | DALES | difference | ratio min / max |
|---|---|---|---|---|
| 25 | 6.839e-3 (1.0e-4) | 6.970e-3 (6.2e-5) | -1.9 % | 2.11 / 1.30 |
| 50 | 1.327e-2 (1.2e-4) | 1.339e-2 (7.3e-5) | -0.9 % | 1.66 / 0.98 |
| 75 | 1.599e-2 (9.0e-5) | 1.606e-2 (6.2e-5) | -0.4 % | 1.01 / 0.70 |
| 100 | 1.648e-2 (8.1e-5) | 1.652e-2 (6.7e-5) | -0.2 % | 0.54 / 0.45 |
| 200 | 1.401e-2 (8.0e-5) | 1.402e-2 (8.3e-5) | -0.1 % | 0.16 / 0.16 |

DALES SGS w'thl' at 25 m is 1.46e-2 of a total 2.16e-2 K m/s, so the resolved part is a third of the flux there.

## Control: prescribed surface flux (uDALES iwalltemp = 1, DALES isurf = 4)

Case 304 with the surface heat flux prescribed at 0.02246 K m/s upward in both codes. That is the 2-4 h mean of the
facet-temperature runs, constant in time, so the boundary layer is shallower (zi 709 m against 736 m). uDALES
`iwalltemp = 1`, `bctfz = -0.02246`, the same release binary; DALES `isurf = 4`, `wtsurf = 0.02246`, `wqsurf = 0`.
Otherwise everything as above, seeds 43-48. Runs `t24/{u304c,d304c}/f4{3..8}`, analysis `t24/an304c`, set-up `t24/setup_ctl.sh`.

| quantity, 6 seeds | uDALES mean (spread) | DALES mean (spread) | diff | ratio min / max |
|---|---|---|---|---|
| u* [m/s] | 0.2812 (0.0011) | 0.2825 (0.0010) | -0.43 % | 1.17 / 1.13 |
| first-level wind |U1| [m/s] | 3.0701 (0.013) | 3.0850 (0.014) | -0.48 % | 1.14 / 1.09 |
| u* / |U1| | 0.091605 (4.9e-5) | 0.091561 (6.9e-5) | +0.05 % | 0.90 / 0.64 |
| zi [m] | 709.0 (1.7) | 709.7 (2.2) | -0.7 m | 0.41 / 0.32 |
| thl profile (rms) | | | | 0.49 / 0.45 |
| resolved w'thl' profile (rms) | | | | 1.41 / 1.13 (largest difference at 425 m) |
| resolved w'thl' at 25 m | 7.199e-3 (7.2e-5) | 7.225e-3 (4.4e-5) | -0.4 % | 0.58 / 0.36 |
| resolved w'thl' at 50 m | 1.337e-2 (9.1e-5) | 1.338e-2 (5.5e-5) | -0.1 % | 0.16 / 0.10 |

With 3 seeds (43-45) the control showed only -0.16 % in u* and -0.19 % in |U1| (0.62 / 0.27 and 0.58 / 0.28 spreads),
and +0.8 % in the resolved w'thl' at 25 m.

The bldump `obl` (`1/mean(1/L)` over the facet sections) is not comparable with the DALES tmser `obukh`, which DALES
computes from the slab-mean state (`oblav`), so it is not in the bar.

## Bit identity of the iwalltemp = 1 cases

Examples 301 and 303 (both `iwallmom = 4`, `iwalltemp = 1`, with first-level cells), Release builds with
`FFTW_ESTIMATE` of 3a428b09 (the parent) and f19067f8, 8 ranks (4 x 2), fixed dt (301: 100 steps of 1 s; 303: 100
steps of 0.3 s): every common variable of `bldump_t` and `bldump_z` identical (numpy `array_equal`, 10 records) and
all 8 restart files byte-identical in both cases. `bldump_t` gains `wtheta` and `wq`.
