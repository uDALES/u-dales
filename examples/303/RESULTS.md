# Tier-2 validation results: 303 DYCOMS-II RF01

Run at commit e2b1a0a5 (Release, Intel/CX3), single node, PBS job 4640941 on cx3 (32 ranks), main loop 5704 s, walltime used 1 h 36 min (limit 3 h), no NaN. Run directory: `/rds/general/user/mvr/ephemeral/t2/303` (outside repo).
Plots: `/rds/general/user/mvr/ephemeral/t2/plots303/timeseries.png`, `/rds/general/user/mvr/ephemeral/t2/plots303/profiles.png`, `summary.md`.

**Published-range checks are pending citable values**: `PUBLISHED` in compare_dales.py is empty by design, so there is no pass/fail table; numbers below are reported without verdicts. Interpretation is tentative.

No DALES reference. Time series (interpolated; window 3-4 h means from compare_dales):

| quantity | 0.5 h | 1 h | 2 h | 3 h | 4 h | 3-4 h mean |
|---|---|---|---|---|---|---|
| LWP [g/m2] | 26.6 | 42.7 | 54.8 | ~60 | ~63 | 59.8 |
| cloud cover | 1.000 | 0.9967 | 0.9978 | 0.9986 | 0.9997 | 0.9995 |
| zi (max dthl/dz) [m] | 850.6 | 858.0 | 869.4 | 879.7 | 889.5 | 884.4 |
| zi_qt8 [m] | 845.4 | 845.5 | 856.9 | 868.5 | 879.3 | - |
| u* [m/s] | 0.255 | 0.245 | 0.253 | 0.250 | 0.248 | 0.248 |

Entrainment rate (zi_qt8, large-scale divergence D = 3.75e-6 1/s, fit 1-4 h): we = 6.26e-3 m/s. (LWP at 0.5-2 h computed from thermostat_t; 3 h and 4 h are read from the plot; the 3-4 h mean is from summary.md.)

Observations (tentative):
- No cloud collapse: cloud cover stays above 0.994 throughout and LWP, after an initial spin-up dip to ~19 g/m2 at 0.8 h, recovers and rises slowly. Stratocumulus-top ql peaks at ~0.43 g/kg near 850 m, inversion near 870-890 m.
- zi rises steadily, about 39 m from 0.5 to 4 h; the net `we` is the sum of that and the subsidence term D*zi (about 3.3e-3 m/s), so about half of the fitted `we` is the prescribed divergence, not a sign of an anomaly. No published-range comparison made (no citable values).
