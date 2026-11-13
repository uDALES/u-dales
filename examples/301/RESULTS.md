# Tier-2 validation results: 301 BOMEX

Run at commit e2b1a0a5 (Release, Intel/CX3), single node, PBS job 4640939 on cx3 (8 ranks), main loop 521 s, walltime used ~9 min (limit 30 min), no NaN. Run directory: `/rds/general/user/mvr/ephemeral/t2/301` (outside repo).
Plots: `/rds/general/user/mvr/ephemeral/t2/plots301/timeseries.png`, `/rds/general/user/mvr/ephemeral/t2/plots301/profiles.png`, `summary.md`.

**Published-range checks are pending citable values**: `PUBLISHED` in compare_dales.py is empty by design, so there is no pass/fail table; numbers below are reported without verdicts. Interpretation is tentative.

**The DALES 3.2.0 reference is not like-for-like.** It ran without the prescribed -2 K/day radiative cooling: its `tltendrad` is exactly 0, and its column thl budget closes (residual 3e-5) only without the cooling. uDALES applies the cooling via `thlpcar` in `lscale.inp.301`, as the BOMEX specification requires, and its budget closes only with it (residual 4e-4). The qt budgets close in both, so the surface flux sign and size are right. See the like-for-like table at the end.

uDALES vs DALES 3.2.0 (`results_3.2.0`), window 3-6 h (uDALES run is 6 h, DALES 16 h):

| quantity | uDALES | DALES |
|---|---|---|
| LWP [g/m2] | 11.37 | 4.34 |
| cloud cover | 0.226 | 0.138 |
| zi [m] | 1679 | 1620 |
| u* [m/s] | 0.28 | 0.28 |

Profile differences (max abs over height, 3-6 h mean): thl 0.76 K, qt 0.55 g/kg, ql 0.0077 g/kg, resolved w'thl' 0.0079 K m/s, resolved w'qt' 0.028 g/kg m/s.

Observations (tentative):
- The gaps in the table above are mostly the radiative-cooling mismatch, not a uDALES defect: see the like-for-like comparison below.
- u* matches (prescribed 0.28 in both), so a surface-stress mismatch is not the cause.
- uDALES has about 2.6x the LWP and 1.6x the cloud cover of DALES. ql profile peaks at about 0.011 g/kg (uDALES) vs 0.0065 g/kg (DALES), with cloud base lower and cloud top higher in uDALES. uDALES LWP is very intermittent (5-20 g/m2).
- The uDALES mixed layer is about 0.5 K colder in thl and slightly moister than DALES near cloud base; zi in uDALES drifts up slightly (1630 to 1700 m) while DALES drifts down.
- Resolved w'thl' is about 30 % larger in magnitude in the cloud layer in uDALES. The deviation is in the cloud layer, so SGS/advection scheme differences (not checked) are candidate causes alongside the thermo port; not diagnosed here.

Like-for-like rerun (same binary and namoptions, `thlpcar = 0` in `lscale.inp.301`, run dir `/rds/general/user/mvr/ephemeral/bdiag/301nr`), window 3-6 h:

| quantity | uDALES with cooling | uDALES `thlpcar = 0` | DALES 3.2.0 |
|---|---|---|---|
| LWP [g/m2] | 11.4 | 7.3 | 4.3 |
| cloud cover | 0.226 | 0.162 | 0.138 |
| zi [m] | 1679 | 1616 | 1620 |

Without the cooling the mixed-layer thl matches DALES to within 0.01-0.03 K and zi agrees. A residual LWP excess of about 1.7x remains. It is not diagnosed. uDALES LWP is very noisy hour to hour, and the candidates are no initial thl/qt noise (randthl/randqt are not applied in modstartup), no Galilean shift, and different SGS constants.
