# Tier-2 validation results: 301 BOMEX

Run at commit e2b1a0a5 (Release, Intel/CX3), single node, PBS job 4640939 on cx3 (8 ranks), main loop 521 s, walltime used ~9 min (limit 30 min), no NaN. Run directory: `/rds/general/user/mvr/ephemeral/t2/301` (outside repo).

**Published-range checks are pending citable values**: `PUBLISHED` in compare_dales.py is empty by design, so there is no pass/fail table; numbers below are reported without verdicts. Interpretation is tentative.

## Reference: DALES dev, run like-for-like

The reference is DALES `dev` at 075ac4e, built on CX3 (Intel 2021a, Release) and run with `cases/bomex` for 6 h on 16 ranks. Only the runtime and the decomposition were changed, and `lsamptend` was switched on for budget output. That run applies the prescribed -2 K/day radiative cooling (`thltendradall` = -2.315e-5 K/s up to ~1500 m), as `examples/301` does via `thlpcar`. Run directory: `/rds/general/user/mvr/ephemeral/dalesbomex`. Plots: `/rds/general/user/mvr/ephemeral/dalesbomex/cmp_cool/{timeseries,profiles}.png`.

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

## A/B test of (a) and (d)

Throwaway builds of this branch (e5634881), same case and ranks, one 6 h realisation each; report and plots in `/rds/general/user/mvr/ephemeral/sgsab`. Variant A removes the factor `2.` in `sbdiss`. Variant B applies thl/qt initial noise (`randthl = 0.1`, `randqt = 2.5e-5`, the DALES BOMEX values, up to `krand`). Window 3-6 h. The last two columns are the max abs difference from DALES dev over height.

| run | LWP [g/m2] | cloud cover | zi [m] | thl [K] | qt [g/kg] |
|---|---|---|---|---|---|
| uDALES as is | 11.37 | 0.226 | 1679 | 0.13 | 0.24 |
| A: no factor 2 | 7.40 | 0.153 | 1717 | 0.24 | 0.36 |
| B: thl/qt noise | 12.52 | 0.212 | 1696 | 0.13 | 0.20 |
| A+B | 8.50 | 0.172 | 1693 | 0.25 | 0.39 |
| DALES dev | 7.53 | 0.164 | 1684 | - | - |

The doubled dissipation (a) accounts for the cloud excess. Noise (d) makes no significant difference: hourly LWP varies by up to 2x, so about 1 g/m2 is within sampling noise. With A the mean profiles move further from DALES and zi rises by about 35 m, so other SGS differences remain. Reported as #401. These runs predate the SGS alignment (Task 10: ch = 1/Prandtl, no factor 2, DALES kb sources); results with the aligned model replace them in the close-match validation.

## The DALES 3.2.0 archive is not a like-for-like reference

`cases/bomex/results_3.2.0` in the DALES repository has `tltendrad` = 0, and its column thl budget closes (residual 3e-5 K m/s) only without the -2 K/day cooling. The uDALES budget closes only with it (residual 4e-4). The qt budgets close in both, so the surface-flux sign and magnitude are right. This was reported upstream as dalesteam/dales#348.

For context, the archive compared with uDALES run with and without the cooling (`thlpcar = 0` in `lscale.inp.301`, run dir `/rds/general/user/mvr/ephemeral/bdiag/301nr`), 3-6 h:

| quantity | uDALES with cooling | uDALES `thlpcar = 0` | DALES 3.2.0 |
|---|---|---|---|
| LWP [g/m2] | 11.37 | 7.28 | 4.34 |
| cloud cover | 0.226 | 0.162 | 0.138 |
| zi [m] | 1679 | 1616 | 1620 |

The no-cooling pair shows the same pattern as the dev pair: matching zi, thl within 0.11 K, and about 1.7x the LWP.
