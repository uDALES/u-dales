# Geostrophic forcing requires `lcoriol`

Regression test for [#410](https://github.com/uDALES/u-dales/issues/410).

`modstartup` builds the large-scale pressure gradient from the geostrophic wind
(`dpdxl = om23_gs*vg - pgx - dpdx`, `dpdyl = -om23_gs*ug - pgy`). The `om23_gs`
terms balance the Coriolis force on the resolved wind, so they must only be
applied with `lcoriol = .true.`. Before the fix a run with `lcoriol = .false.`
and nonzero `ug`/`vg` was accelerated at `f*|u_g|`.

The test runs case 300 (laminar, no buildings) from rest with `ug = 10 m/s` and
reads the horizontal-mean wind from the restart file:

| `lcoriol` | expected |
|---|---|
| `.false.` | wind stays at rest (`|mean u|, |mean v| < 1e-4 m/s`; the bug gives `v = f*ug*t = 2.5e-3 m/s` after 2.2 s) |
| `.true.` | unchanged: `v = om23_gs*ug*t` to within 1 % |

Environment: `UDALES_BUILD`, `UDALES_RUNTIME_MODULES`, `MPIEXEC` (as for the other
solver integration tests).
