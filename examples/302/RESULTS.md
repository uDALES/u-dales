# Tier-2 validation results: 302 smoke

Early run (Release, Intel/CX3, single node, 16 ranks), main loop 1250 s, no NaN.

**Published-range checks are pending citable values**: `PUBLISHED` in compare_dales.py is empty by design, so there is no pass/fail table; numbers below are reported without verdicts. Interpretation is tentative.

uDALES vs DALES 3.2.0, window 2-3 h:

| quantity | uDALES | DALES |
|---|---|---|
| LWP, cloud cover | 0 | 0 (cfrac nan) |
| zi (smoke, diagnostic) [m] | 731.5 | 396.6 |
| u* [m/s] | 0.01357 | 0.01278 |
| entrainment rate we, fit 1-3 h [m/s] | 0.00339 | -0.00073 |

Profile differences (max abs): thl 1.06 K, qt 5e-4 g/kg (DALES noise, uDALES qt = 0), ql 0, resolved w'thl' 0.0039 K m/s, resolved w'qt' 2.7e-4 g/kg m/s.

Observations (tentative):
- The mean thl profiles overlap visually (mixed layer ~287.6 K, inversion near 720-730 m, 295 K above); the 1.06 K maximum difference sits at the inversion step. The resolved w'thl' profile has the same shape and a ~5 % smaller maximum in uDALES.
- The `zi` series differ (731 vs 397 m) most likely because of a diagnostic-definition difference rather than a physical one (not verified): the DALES `zi` series (about 400 m after 1 h) disagrees with the inversion height visible in the DALES thl profile (~720 m). The entrainment-rate row inherits this and should not be read as a physical discrepancy.
- u* is within 6 %; both show the same intermittent spikes (~1.4 h, ~2.4 h).
