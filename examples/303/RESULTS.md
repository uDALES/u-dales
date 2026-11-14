# Tier-2 validation results: 303 DYCOMS-II RF01

## Close-match validation against DALES dev (Task 19, 2026-11-14): corrected DALES reference

**Verdict: very close.** With identical cd2 numerics and a fixed dt of 0.3 s, uDALES and DALES agree within the
seed-to-seed spread for LWP, cloud cover, u*, zi_qt8 and the entrainment rate, and for all five profiles the rms
difference is below the smaller rms spread (0.18-0.65). Level by level the bar is stricter: the fraction of levels
with |diff| below both spreads is 42 % (thl), 72 % (qt), 36 % (ql), 62 % (w'thl'), 86 % (w'qt'). At 20-40 min the
LWP and entrainment differences are about 7 times the (then very small) spread; see "Early times". The one item above the strict bar is zi from the maximum
dthl/dz (0.8 m, a sixth of the 5 m grid spacing; DALES members 875.8, 875.8, 874.2 m against 874.2, 874.2, 875.0;
ratio 1.7 with the smaller spread, 0.9 with the larger).

The Task 14 mismatch below (uDALES entraining twice as fast) was a **defect in the DALES reference**, not in
uDALES: DALES dev 075ac4e `modradiation.f90:initradiation` sets `rad_longw = rad_shortw = rad_smoke = .false.`
for `iradiation = 10`, so the `radpar` call in its `irad_user` branch never runs. The Task 14 DALES runs had
therefore no cloud-top longwave cooling (F0 = 70, F1 = 22 W/m2), only the free-troposphere term of `rad_user`.
Found with a horizontally homogeneous single-column run of both codes (all noise off, so advection and pressure
vanish): every process tendency agreed to 3-4 digits except radiation, which was exactly zero in the DALES cloud.
The fix is the one-line hunk in `tools/python/examples/dales_cases/dales_inputs/dales-075ac4e-t19.patch`.

The 2.3 % u* gap of Task 14 was the von Karman constant (DALES 0.4, uDALES 0.41): the MO drag coefficient scales
with kappa^2. **The runs reported here used kappa = 0.41 in both codes** (DALES built with `fkar = 0.41` in
`modglobal.f90`). From the change of the uDALES default `fkar` to 0.40 (the DALES and PALM value) onward, both
codes use 0.40 and DALES needs no `fkar` edit.

Note (Task 20): the reported runs used kappa = 0.41 in both codes; the committed `namoptions.303` and
`dales_inputs/rf01/namoptions.001` now use kappa = 0.40 with the recomputed `z0 = z0h = 1.445944e-5 m` (same neutral
Cm = 0.0011). Not rerun.

Runs: uDALES `/rds/general/user/mvr/ephemeral/t14/x303/f4{3,4,5}` (binary code-equivalent to 674395fe, 8 x 4 ranks;
the same runs as Task 14). DALES `/rds/general/user/mvr/ephemeral/t19/f303/dfix_4{3,4,5}` (dev 075ac4e + the t14
and t19 patches, `fkar = 0.41`, inputs `dales_inputs/rf01`, 8 x 4 ranks), seeds 43, 44, 45; window 3-4 h;
`compare_dales.py 303 ...`, output `/rds/general/user/mvr/ephemeral/t19/an/303fixed_t19`. Settings as in the
Task 14 table below, with radiation now matched (radpar longwave every step in both codes) and kappa matched.

| quantity | uDALES members | uDALES mean (spread) | DALES members | DALES mean (spread) | diff | ratio (min / max spread) |
|---|---|---|---|---|---|---|
| LWP [g/m2] | 14.6, 16.3, 16.6 | 15.8 (1.1) | 16.5, 17.1, 14.9 | 16.2 (1.2) | -0.35 | 0.33 / 0.30 |
| cloud cover | 0.775, 0.839, 0.827 | 0.813 (0.034) | 0.841, 0.801, 0.804 | 0.815 (0.023) | -0.002 | 0.09 / 0.06 |
| u* [m/s] | 0.2470, 0.2470, 0.2458 | 0.2466 (0.0007) | 0.2464, 0.2466, 0.2469 | 0.2466 (0.0002) | -0.0001 | 0.29 / 0.09 |
| zi, max dthl/dz [m] | 874.2, 874.2, 875.0 | 874.4 (0.5) | 875.8, 875.8, 874.2 | 875.3 (1.0) | -0.8 | 1.73 / 0.87 |
| zi_qt8 [m] | 858.0, 858.6, 859.7 | 858.8 (0.9) | 859.3, 858.7, 858.5 | 858.8 (0.4) | -0.05 | 0.10 / 0.05 |
| we = dzi_qt8/dt + D zi, fit 1-4 h [mm/s] | 4.68, 4.81, 4.91 | 4.80 (0.11) | 4.87, 4.82, 4.75 | 4.81 (0.06) | -0.01 | 0.24 / 0.13 |

| profile (z <= 1200 m) | max abs diff (z) | rms diff | rms spread U / D | rms diff / min rms spread | levels with diff < both spreads |
|---|---|---|---|---|---|
| thl [K] | 0.083 (872 m) | 0.0095 | 0.033 / 0.016 | 0.59 | 42 % |
| qt [g/kg] | 0.048 (878 m) | 0.0082 | 0.025 / 0.016 | 0.52 | 72 % |
| ql [g/kg] | 0.0021 (778 m) | 0.0007 | 0.0022 / 0.0023 | 0.31 | 36 % |
| w'thl' resolved [K m/s] | 0.0005 (760 m) | 0.0002 | 0.0003 / 0.0004 | 0.65 | 62 % |
| w'qt' resolved [g/kg m/s] | 0.0006 (15 m) | 0.0002 | 0.0015 / 0.0012 | 0.18 | 86 % |

Hourly LWP [g/m2], ensemble mean (spread): uDALES 27.4 (0.2), 14.0 (0.4), 15.3 (0.7), 15.8 (1.1); DALES 26.9 (0.1),
14.7 (0.2), 15.2 (0.6), 16.2 (1.2).

Early times (20-40 min, `t19/an/303s_fix`; the members barely differ yet, so the spreads are tiny): LWP 17.0 vs 15.8
g/m2 (+7 %, ratio 7.1 with the smaller spread), we 9.6 vs 10.0 mm/s (-4 %, ratio 7.5), so above the bar at that
time (it closes by 1 h, hourly LWP above); zi_qt8 845.2 vs 845.9 m, u* 0.2716 vs 0.2717 m/s. With the old reference the same
window gave LWP 17.0 vs 22.0 and we 9.6 vs 5.9. A DALES run with the radiation hunk but kappa = 0.4 (seed 43,
`t19/an/303s_rad`) gives the same entrainment and LWP as the full fix and the old 2 % u* deficit, so the radiation
defect explains the entrainment/LWP gap and kappa the u* gap.

## Close-match validation against DALES dev (Task 14, 2026-11-14): superseded, DALES reference without cloud-top cooling

The DALES runs of this section (`t14/d303/*`) were made with `iradiation = 10` without the t19 patch, so DALES had
no radpar longwave cooling (see the Task 19 section above), and with kappa = 0.4. The uDALES runs are the ones used
above. Kept for the record.

### Task 14 text (as written then)


**Caveat (Task 19):** every DALES run in this section, including the "standard numerics" run, had no cloud-top
longwave cooling. Conclusions drawn here about DALES behaviour, cloud loss and the role of cd2 do not hold. With the
corrected reference both codes keep LWP of about 16 g/m2 at 3-4 h with cd2 (Task 19 section above).

**Verdict (then): not close.** With identical cd2 numerics and a fixed dt of 0.3 s, LWP falls from 60 to 10-15 g/m2
in the first hour in both runs. But uDALES keeps more cloud and entrains faster than DALES, far outside the seed spread: LWP 15.8 vs
9.4 g/m2, cloud cover 0.81 vs 0.63, zi_qt8 859 vs 834 m, we 4.8 vs 2.6 mm/s, u* +2.3 %. Diagnosed in Task 19 (see
"Resolution" below). All uDALES numbers here include the e12 halo fix (issue #405, commit 674395fe; the Task 14
binary `u-dales.c220463a+e12halo` is code-equivalent to 674395fe, see `examples/301/RESULTS.md`); c220463a has the bug
(adaptive run shown for reference).

Runs: `/rds/general/user/mvr/ephemeral/t14`: uDALES `x303/f4{3,4,5}` (c220463a + e12 halo fix, 32 ranks 8 x 4),
DALES `d303/f4{3,4,5}` (dev 075ac4e + `dales_inputs/dales-075ac4e-t14.patch`, inputs `dales_inputs/rf01`, 32 ranks
8 x 4), seeds 43, 44, 45; window 3-4 h; `compare_dales.py 303 ...`, outputs in `t14/an/303*`.

### Identical settings, side by side

As in `examples/301/RESULTS.md` (advection, time scheme, Galilean shift, base state, exner, constants, es, saturation
adjustment, SGS, subsidence, sponge, top BC, noise), with these case-specific rows:

| setting | uDALES (examples/303) | DALES (dales_inputs/rf01) | matched? |
|---|---|---|---|
| grid | 96 x 96 x 300, 35 m x 35 m x 5 m | same | yes |
| time step | fixed 0.3 s | fixed 0.3 s (48000 steps) | yes |
| surface | `BCbotm = 4` (MO, prescribed fluxes, u* diagnosed), z0 = z0h = 1.07e-5 m, wtsurf 1.2246e-2 K m/s, wqsurf 3.7258e-5 up | `isurf = 4`, z0mav = z0hav, same fluxes (rlv 2.53e6, rho 1.22) | settings yes; u* differs by 2.3 % from the first minute (finding 2). DALES ran per-column MO (`lmostlocal = .true.`, its dev default), uDALES one domain-mean Obukhov length and u*: not matched (with kappa aligned, Task 19, u* agrees to 0.0001 m/s) |
| radiation | `iradiation = 2` (radpar longwave: rka 85, F0 70, F1 22), every step | `iradiation = 10` (radpar + rad_user), same values, every step | yes (same formula, both computed at the first RK substep) |
| free-troposphere term | flux-divergence form at half levels, D = 3.75e-6, qt = 8 g/kg zi per column (case hook) | patched rad_user, identical | yes |
| noise | thl +-0.1 K, qt +-2.5e-5, all 300 levels, no velocity noise | same (`krandumax = 0`) | amplitudes yes (generators differ) |
| sponge | `igrw_damp = 2`, ksp = 225 (zf = 1122.5 m) | default `igrw_damp = 2`, same | yes |
| Coriolis, geostrophic wind, subsidence | xlat 31.5; ug = (7, -5.5); wfls = -3.75e-6 z | same `lscale.inp` | yes |

Fixed dt: DALES adaptive (diffusion-limited, peclet 0.15) hourly means 0.43, 0.40, 0.42, 0.42 s; uDALES + fix
(`diffnr = 0.16`) 0.47-0.50 s, minimum 0.43 s. 0.3 s is below both and stable in all runs.

### Results, fixed dt = 0.3 s, window 3-4 h

| quantity | uDALES members | uDALES mean (spread) | DALES members | DALES mean (spread) | diff | ratio (min / max spread) |
|---|---|---|---|---|---|---|
| LWP [g/m2] | 14.6, 16.3, 16.6 | 15.8 (1.1) | 7.9, 9.9, 10.3 | 9.4 (1.3) | +6.4 | 6.0 / 5.1 |
| cloud cover | 0.775, 0.839, 0.827 | 0.813 (0.034) | 0.613, 0.648, 0.641 | 0.634 (0.019) | +0.18 | 9.6 / 5.3 |
| u* [m/s] | 0.2470, 0.2470, 0.2458 | 0.2466 (0.0007) | 0.2406, 0.2412, 0.2410 | 0.2409 (0.0003) | +0.0056 | 18 / 8 |
| zi, max dthl/dz [m] | 874.2, 874.2, 875.0 | 874.4 (0.5) | 850, 850, 850 | 850 (0) | +24 | - / 51 |
| zi_qt8 [m] | 858.0, 858.6, 859.7 | 858.8 (0.9) | 833.8, 833.8, 833.3 | 833.6 (0.3) | +25 | 85 / 29 |
| we = dzi_qt8/dt + D zi, fit 1-4 h [mm/s] | 4.68, 4.81, 4.91 | 4.80 (0.11) | 2.59, 2.59, 2.52 | 2.57 (0.04) | +2.2 | 55 / 20 |

| profile (z <= 1200 m) | max abs diff (z) | rms diff | rms spread U / D | rms diff / min rms spread |
|---|---|---|---|---|
| thl [K] | 9.5 (862 m) | 1.28 | 0.033 / 0.010 | 123 |
| qt [g/kg] | 6.3 (862 m) | 0.89 | 0.025 / 0.010 | 89 |
| ql [g/kg] | 0.12 (848 m) | 0.020 | 0.002 / 0.002 | 9 |
| w'thl' resolved [K m/s] | 0.030 (855 m) | 0.004 | 0.0003 / 0.0012 | 13 |
| w'qt' resolved [g/kg m/s] | 0.026 (855 m) | 0.009 | 0.0015 / 0.0019 | 6 |

The maxima sit at the inversion and mostly measure its 25 m offset. Inside the boundary layer (3-4 h means, uDALES /
DALES): thl 289.47 / 289.42 K at 100 m, 289.78 / 290.11 at 800 m; qt 9.24 / 9.54 g/kg at 100 m, 8.61 / 8.57 at 800 m;
ql 0.102 / 0.072 g/kg at 800 m; resolved w'qt' at 400 m 0.038 / 0.028 g/kg m/s. uDALES is better mixed (thl rises
0.29 K from 400 to 800 m against 0.67 K in DALES), i.e. DALES decouples more.

Hourly LWP [g/m2], ensemble mean (spread):

| run | 0-1 h | 1-2 h | 2-3 h | 3-4 h |
|---|---|---|---|---|
| uDALES + fix, dt 0.3 s | 27.4 (0.2) | 14.0 (0.4) | 15.3 (0.7) | 15.8 (1.1) |
| DALES, dt 0.3 s | 30.1 (0.1) | 10.1 (0.3) | 9.7 (1.0) | 9.4 (1.3) |
| uDALES + fix, adaptive, seed 43 | 27.9 | 13.5 | 14.5 | 14.2 |
| DALES, adaptive, seed 43 | 30.2 | 9.6 | 10.0 | 8.7 |
| uDALES c220463a (e12 halo bug), adaptive, seed 43 | 22.4 | 14.6 | 17.3 | 24.5 |
| DALES standard numerics (iadv 5, cu 7, cv -5.5, anelastic ibas_prf 3), adaptive, seed 43 | 45.8 | 23.4 | 15.1 | 11.7 |

Adaptive pair (secondary, seed 43): LWP 14.2 vs 8.7, cloud cover 0.79 vs 0.63, zi_qt8 857 vs 834 m, we 4.6 vs 2.6 mm/s,
u* 0.247 vs 0.241: the same gap as at fixed dt. The DALES standard-numerics run (iadv 5, `t14/d303/std43`) used the
same `iradiation = 10` inputs, so it also had no cloud-top cooling; the cd2 vs iadv 5 comparison has to be redone
with the corrected DALES (a Task 15 / advection-PR item).

### The DALES RF01 reference

The first DALES RF01 run (`refruns/d303`) had a constant `divmax = 2.44e-2`: w noise at the surface level
(`krandumax = 168` with the default `krandumin = 1`; divmax = rhobf(1) randu / dzf exactly). DALES also has a
sponge (`igrw_damp = 2` is the default). Task 14 removed the surface noise and set the sponge explicitly. Its
conclusion that "DALES cd2 loses the cloud" in every variant (`t14/diag/kr2`, `t14/diag/ftnew`, `t14/d303/*`) came
from runs without cloud-top longwave cooling and is withdrawn. The fixed DALES inputs are in
`tools/python/examples/dales_cases/dales_inputs/rf01` and need both DALES patches (t14 and t19).

### Resolution (Task 19)

The Task 14 hypotheses (SGS at the inversion, MO surface, radiation, decomposition) were settled by one
homogeneous single-column comparison: the DALES radpar tendency was exactly zero (`iradiation = 10` resets
`rad_longw` in DALES dev 075ac4e), and kappa (0.4 vs 0.41) explained the u* gap. A 2 x 2 vs 8 x 4 uDALES check showed
no decomposition dependence. The "better mixed / DALES decouples more" difference was the missing DALES cloud-top
cooling. See the Task 19 section at the top and `.superpowers/sdd/2026-11-13-moist-thermo-port/task-19-report.md`.

## Earlier results (superseded; e12 halo bug and older SGS)

Run at commit e2b1a0a5 (Release, Intel/CX3), single node, PBS job 4640941 on cx3 (32 ranks), main loop 5704 s, walltime used 1 h 36 min (limit 3 h), no NaN. Run directory: `/rds/general/user/mvr/ephemeral/t2/303` (outside repo).
Plots: `/rds/general/user/mvr/ephemeral/t2/plots303/timeseries.png`, `/rds/general/user/mvr/ephemeral/t2/plots303/profiles.png`, `summary.md`.

**Published-range checks are pending citable values**: `PUBLISHED` in compare_dales.py is empty by design, so there is no pass/fail table; numbers below are reported without verdicts. Interpretation is tentative.

No DALES reference. Time series (interpolated; window 3-4 h means from compare_dales):

| quantity | 0.5 h | 1 h | 2 h | 3 h | 4 h | 3-4 h mean |
|---|---|---|---|---|---|---|
| LWP [g/m2] | 26.6 | 42.7 | 54.8 | ~60 | ~63 | 59.8 |
| cloud cover | 1.000 | 0.9967 | 0.9978 | 0.9986 | 0.9997 | 0.9995 |
| zi (max dthl/dz) [m] | 850.6 | 858.0 | 869.4 | 879.7 | 889.5 | 884.4 |
| zi_qt8 [m] (now `zi_qt`, `zi_qt_thres = 8e-3`) | 845.4 | 845.5 | 856.9 | 868.5 | 879.3 | - |
| u* [m/s] | 0.255 | 0.245 | 0.253 | 0.250 | 0.248 | 0.248 |

Entrainment rate (zi_qt8, large-scale divergence D = 3.75e-6 1/s, fit 1-4 h): we = 6.26e-3 m/s. (LWP at 0.5-2 h computed from thermostat_t; 3 h and 4 h are read from the plot; the 3-4 h mean is from summary.md.)

Observations (tentative):
- No cloud collapse: cloud cover stays above 0.994 throughout and LWP, after an initial spin-up dip to ~19 g/m2 at 0.8 h, recovers and rises slowly. Stratocumulus-top ql peaks at ~0.43 g/kg near 850 m, inversion near 870-890 m.
- zi rises steadily, about 39 m from 0.5 to 4 h; the net `we` is the sum of that and the subsidence term D*zi (about 3.3e-3 m/s), so about half of the fitted `we` is the prescribed divergence, not a sign of an anomaly. No published-range comparison made (no citable values).
