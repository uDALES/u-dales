# 301: BOMEX shallow cumulus (Siebesma et al. 2003)

Full-size validation case for the moist thermodynamics, the Monin-Obukhov (MO) facet surface and
bldump diagnostics. Inputs are converted from the DALES case `cases/bomex`
(`tools/python/examples/dales_cases/convert_dales_case.py`), grid and forcing unchanged.
Run on CX3 with `tools/hpc_execute.sh examples/301` (8 ranks, 4 x 2); compare with
`tools/python/examples/dales_cases/compare_dales.py` (see the README there).
Outputs of interest: `bldump_t.301.nc`, `bldump_z.301.nc`.

## Floor

The floor is one horizontal facet covering the domain (`libm = .true.`, `lbottom = .false.`, `nfcts = 1`),
as in every uDALES case: the flat bottom boundary (`lbottom`) is not used. Its IBM files (`facets.inp`,
`factypes.inp`, `facetarea.inp`, `solid_*.txt`, `fluid_boundary_*.txt`, `facet_sections_*.txt`) were written by
`python -m udprep.flat_floor examples/301 301 --itot 64 --jtot 64 --dx 100 --dy 100 --z1 20 --z0 2e-4 --z0h 2e-4`
(from `tools/python`), which prints the `&WALLS` counts. The surface layer is the Monin-Obukhov facet scheme
(`iwallmom = 4`, the DALES surface with a local Obukhov length per facet section). The floor's facet type (-1)
has `lGR = 1` in `factypes.inp.301`: in uDALES this green flag only means that the facet exchanges
moisture (`iwallmoist = 1` applies `bcqfz` to it)).

## Settings and sources

| setting | value | source |
|---|---|---|
| grid | 64 x 64 x 80, dx = dy = 100 m, dz = 40 m (z up to 3200 m) | DALES `itot, jtot, kmax, xsize, ysize`; z from `prof.inp.001` |
| duration | 21600 s (6 h) | case setup; the DALES file has `runtime = 28800` |
| time step | `ladaptive`, `courant = 0.55`, `diffnr = 0.2`, `dtmax = 20` | DALES `ladaptive`; `courant` tuned so that the mean dt matches the DALES run (see "Time step"); `diffnr` = DALES `peclet = 0.15` x 1.32 (not limiting); `dtmax` = DALES default |
| `ps` | 101500 Pa | DALES `ps` |
| `thls` | 299.1 K | DALES `thls`; the reference surface value (`thvs`, base state); the floor facet diagnoses its own surface thl every step |
| `qts` | 1.7e-2 kg/kg | `qt` at the first level; reference value only, as `thls` |
| `bctfz` | -8.0e-3 K m/s | DALES `wtsurf = +8e-3` (upward), negated: uDALES prescribed fluxes are positive downward (D1) |
| `bcqfz` | -5.2e-5 kg/kg m/s | DALES `wqsurf = +5.2e-5`, negated |
| `iwallmom`, `ustin` | 4 (MO facets), 0.28 m/s (prescribed u*) | DALES `isurf = 3`, `ustin` |
| `z0`, `z0h` | 2e-4 m (floor facet type, `factypes.inp.301`) | DALES `z0` sets both |
| Coriolis | `lcoriol`, `xlat = 15` (f = 0.3775e-4 1/s) | DALES `xlat = 15` (Siebesma et al. 2003 give f = 0.376e-4; the case file is followed) |
| forcing | `lscale.inp.301` | DALES `lscale.inp.001`: `ug, vg, wfls, dqtdtls`; `dthlrad` becomes `thlpcar` (-2.315e-5 K/s) |
| radiation | `iradiation = 0`; the cooling is the `thlpcar` column | DALES (no interactive radiation in BOMEX) |
| SGS | `loneeqn = .true.` (1.5-order TKE), `ldelta = .false.` | DALES `ldelta = .false.`; one-equation closure as DALES |
| advection | `iadv_mom = iadv_tke = iadv_thl = iadv_qt = 2` | DALES `iadv_* = 2` |
| sponge | `igrw_damp = 2`, `ksp` default (top quarter, `rnu0 = 2.75e-3`) | DALES default `igrw_damp = 2`; the `ksp` default is the same formula |
| `lqlnr` | `.true.` | exact liquid water by Newton-Raphson; DALES dev uses two linearised steps (within 0.2 % of exact) and ignores its `lqlnr` |
| surface | `libm = .true.`, `lbottom = .false.`, one floor facet, `iwalltemp = iwallmoist = 1` | a facet floor without buildings (see "Floor") |
| bldump | `lbldump`, `tstatsdump = 60` s | case setup |
| `timerad` | 0 | case setup (radiation every step; unused with `iradiation = 0`) |
| initial noise | `lrandomize`, `randthl = 0.1` K, `randqt = 2.5e-5` kg/kg, `randu = 0`, `krand = 80` (all levels), `irandom = 43` | DALES `randthl`, `randqt`, `irandom`; DALES BOMEX applies no u, v, w noise (its momentum noise is off) and perturbs all levels |
| `tcheck` | 60 s | mean dt to the log every minute (used for the time-step comparison) |
| decomposition | `nprocx = 4`, `nprocy = 2`; `jtot % nprocy = 0` and `ktot % nprocy = 0` | |

`prof.inp.301` row 1 equals the DALES row 1 column for column (checked by the converter,
which re-reads its output). The `lscale.inp.001` row at z = 3220 m (above the model top) is dropped.

## Time step

DALES BOMEX (dev 075ac4e, `courant` default 1.0, `peclet` 0.15, seed 43, 64 x 64 x 80, `cu = -8`) is limited by
diffusion (5528 of its 5895 steps); its Courant number (L2) stays at about 0.6. uDALES without the Galilean shift
is limited by the Courant number of the 9 m/s mean wind and, later, of the cloud updrafts; the diffusion number is
never the limiting criterion for `diffnr >= 0.2`. Mean dt in s (step-weighted, so the average dt per step; uDALES from
the 60 s `tcheck` means), 6 h runs on CX3 with `build/release`:

| window | DALES | uDALES `courant` 0.55 (final) | 0.5 | 0.6 | 1.0 (`diffnr` 0.25) |
|---|---|---|---|---|---|
| 0-1 h | 5.59 | 5.54 | 5.03 | 6.05 | 10.2 |
| 1-3 h | 3.82 | 4.81 | 4.32 | 5.27 | 9.5 |
| 3-6 h | 3.21 | 3.07 | 2.88 | 3.35 | 5.45 |
| 0-6 h | 3.66 | 3.81 | 3.52 | 4.17 | 6.98 |

`courant = 0.55` matches the DALES mean within 5 % over 0-1 h, 3-6 h and 0-6 h; hours 1-3 are 26 % longer, because DALES's
diffusion limit and uDALES's Courant limit do not evolve alike and one parameter pair cannot track both. `diffnr = 0.2`
and 0.25 gave identical series; `diffnr = 0.1` made diffusion limit the step (4.6, 4.3, 2.8 s in the three windows with `courant = 1`).

For a **fixed** step (needed for an identical-dt comparison with DALES, where `ladaptive = .false.` and `dtmax`
is the step), uncomment the two commented lines in `&RUN` of `namoptions.301`: `ladaptive = .false.`, `dtmax = 2.0`;
`courant` and `diffnr` are then ignored and `runtime` should be a multiple of `dtmax`. 2.0 s is below the adaptive step
of both codes in the close-match runs (RESULTS.md); DALES needs `ladaptive = .false.` and `dtmax = 2.0`, too
(`tools/python/examples/dales_cases/dales_inputs/bomex`). The time-step table above was measured with the e12 halo
bug of RESULTS.md present; with that bug fixed, uDALES's adaptive step is shorter (hourly means 2.3-4.0 s).

## Numerics compared with DALES

Short version; the item-by-item comparison and the reasons are in [uDALES and DALES numerics](../../docs/udales-dales-numerics.md), the full settings of the
reported runs in [RESULTS.md](RESULTS.md).

| item | uDALES 301 | DALES reference (`dales_inputs/bomex`) | matched? |
|---|---|---|---|
| advection | cd2 for u, v, w, e12, thl, qt | cd2 (the DALES BOMEX case too) | yes |
| time step | fixed 2.0 s for the comparison (adaptive L1 `courant = 0.55` by default) | fixed 2.0 s (adaptive: L2 `courant = 1.0`, `peclet = 0.15`) | fixed yes; adaptive no |
| base state, Exner | Boussinesq; Exner updated every step | `ibas_prf = 5`, constant density; Exner frozen at the uDALES initial profile | by construction; Exner drift about 1e-5 |
| Galilean shift | none | `cu = 0` (the DALES case uses `cu = -8`) | yes; the shift itself is not tested |
| SGS | one-equation TKE, DALES constants | the same | yes |
| saturation adjustment, es | exact (`lqlnr`), Tetens | two-step (within 0.2 %), Tetens (patched) | nearly |
| surface | MO facet floor, u* = 0.28 prescribed, L per facet section | `isurf = 3`, per-column L (`lmostlocal` default) | yes |
| sponge, noise | `igrw_damp = 2`, `ksp = 60`; thl 0.1 K, qt 2.5e-5 on all levels | the same | yes (different generators) |

## Deviations from DALES

- MO surface: `iwallmom = 4` with `ustin` follows DALES `isurf = 3` (u* = ustin, prescribed fluxes) on the
  floor facet; the surface state lives on the u and v points of the facet sections, not on the cell centres.
- No Galilean shift: DALES runs with `cu = -8` (u is shifted by 8 m/s so that the domain-mean
  wind is small). uDALES has no `cu`; `prof.inp` and `lscale.inp` carry the physical winds
  (u about 9 m/s), so the uDALES Courant number is dominated by the advection of the mean wind
  (see "Time step").
- Courant norm: uDALES limits the L1 sum `|u|/dx + |v|/dy + |w|/dz`, DALES the L2 norm, and DALES's
  default `courant` is 1.0 (not 0.5); `courant` is therefore a tuned value here, not a transcribed one.
- uDALES has no 5th-order advection scheme; BOMEX uses second-order central differences in DALES too. The
  remaining differences are listed in [uDALES and DALES numerics](../../docs/udales-dales-numerics.md).
- Duration 6 h instead of 8 h (case setup).
- Output: only the bldump files (statistics, field dumps and cross sections are off).

## Concerns

- The dt match is a 6 h-mean match (see "Time step"); hours 1-3 run about 25 % longer steps than DALES.
- The `qts` and `thls` namelist values are reference values only; the comparison should use the
  diagnosed `thls`, `qts`, `ustar` and `obl` in `bldump_t` (area means over the floor sections).
