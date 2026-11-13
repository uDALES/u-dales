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
- DALES has cloud within the first hour; uDALES needs about two hours. uDALES applies no initial thl/qt noise because `randthl`/`randqt` are never applied (#395). This is the leading candidate for the residual. Others are not using a Galilean shift and different SGS constants. None of these has been tested.

## The DALES 3.2.0 archive is not a like-for-like reference

`cases/bomex/results_3.2.0` in the DALES repository has `tltendrad` = 0, and its column thl budget closes (residual 3e-5 K m/s) only without the -2 K/day cooling. The uDALES budget closes only with it (residual 4e-4). The qt budgets close in both, so the surface-flux sign and magnitude are right. This was reported upstream as dalesteam/dales#348.

For context, the archive compared with uDALES run with and without the cooling (`thlpcar = 0` in `lscale.inp.301`, run dir `/rds/general/user/mvr/ephemeral/bdiag/301nr`), 3-6 h:

| quantity | uDALES with cooling | uDALES `thlpcar = 0` | DALES 3.2.0 |
|---|---|---|---|
| LWP [g/m2] | 11.37 | 7.28 | 4.34 |
| cloud cover | 0.226 | 0.162 | 0.138 |
| zi [m] | 1679 | 1616 | 1620 |

The no-cooling pair shows the same pattern as the dev pair: matching zi, thl within 0.11 K, and about 1.7x the LWP.
