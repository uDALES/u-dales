# Tier-2 validation results: 301 BOMEX

Note on the von Karman constant: the runs reported here used kappa = 0.41 in uDALES (its default then) and
kappa = 0.4 in DALES (built in). The uDALES default `fkar` is now 0.40, the DALES value. BOMEX has a prescribed u*
(`BCbotm = 5`), so kappa enters only through the surface-layer gradients (see the settings table); the runs were not
repeated with `fkar = 0.40`. `namoptions.301` prescribes z0 = 2e-4 m independently of kappa, so it needs no change.

## Close-match validation against DALES dev

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
| diagnostics | 60 s slab-mean snapshots (thermostat) | 60 s samples, 600 s block means (genstat) | aligned in compare_dales.py: same 600 s blocks, same zi definition from the block-mean profiles |

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
