# uDALES and DALES numerics compared

This page lists, item by item, where the numerics of uDALES agree with classical DALES and where they differ,
and which settings the DALES validation cases (BOMEX, DYCOMS-II RF01, smoke) use. It is meant to be kept up to
date: when either code changes one of the items below, change the row.

**Versions.** uDALES: the moist thermodynamics port, cited as `src/<file>:<line>` as of the pull request that
introduced this page (line numbers drift; the routine name is the stable key). DALES: `dev` at 075ac4e, cited as `DALES <file>:<line>` (paths under `src/`). Where DALES 4.4.2 (the
last release checked) differs from dev, the row says so.

**Status labels.**

| label | meaning |
|---|---|
| aligned | the same formula in both codes, checked in the code at the versions above |
| aligned in the validation | the defaults differ, but the validation runs set both codes to the same choice |
| different by design | a deliberate difference, left in place; its effect is stated |
| different, open | a difference that should be removed or tested; the follow-up is named |

**Validation.** The DALES comparison runs both codes with the same numerics: cd2 advection for every variable,
no Galilean shift, a constant-density DALES base state and an identical fixed time step. The DALES inputs and
patches are in `tools/python/examples/dales_cases/dales_inputs`; the runs, settings and results are in
[301 results](https://github.com/uDALES/u-dales/blob/master/examples/301/RESULTS.md) and
[303 results](https://github.com/uDALES/u-dales/blob/master/examples/303/RESULTS.md).

## Summary: what still differs

| # | item | status | effect on BOMEX / RF01 |
|---|---|---|---|
| 1 | Boussinesq (uDALES) vs anelastic (DALES default `ibas_prf = 3`) | different by design | not measured: the validation used a constant-density DALES. An offline estimate gives about 8-10 % more mass-weighted surface-flux moistening of the BOMEX cloud layer in an anelastic DALES (about 0.05-0.1 g/kg qt over 6 h) |
| 2 | Exner function: updated every step (uDALES) vs frozen at the base state (DALES) | different by design | in the validation DALES's frozen Exner is uDALES's initial one; the drift over a run is of order 1e-5 |
| 3 | Saturation adjustment: exact Newton-Raphson (uDALES `lqlnr`) vs two linearised steps (DALES dev) | different by design | DALES two-step is within about 0.2 % of the exact ql for BOMEX and RF01 states |
| 4 | Galilean shift: none in uDALES; DALES BOMEX uses `cu = -8` | different by design | not measured (the validation used `cu = 0` in DALES); changes the cd2 dispersion error and the dt |
| 5 | Courant number: L1 sum (uDALES) vs L2 norm (DALES); diffusion number vs Peclet number | different by design | the adaptive steps differ between the codes (measured steps in the case RESULTS); the comparison uses a fixed dt |
| 6 | Advection schemes available | different, open (advection PR) | the validation used cd2 for every variable in both codes, while DALES's GCSS stratocumulus case uses 5th order; cd2 for qt and thl at the 5 m RF01 inversion is the setting most likely to matter. The cd2 vs 5th-order comparison with DALES has not been done |
| 7 | MO surface: per facet section on a facet floor (uDALES `iwallmom = 4`) vs per column (DALES dev default `lmostlocal = .true.`) | aligned (local in both) | the u, v points of uDALES carry their own L and u* for the stress, DALES averages the u* of two columns at a u point; second-order equivalent. The first-level SGS treatment uses a per-cell state from the cell-centre wind, as DALES per column (identical to the u/v state in a homogeneous column; RESULTS of 301, 303) |
| 8 | Saturation vapour pressure: Tetens (uDALES) vs Murphy-Koop table (DALES) | aligned in the validation (DALES patched) | Tetens is 0.10-0.13 % lower, about +1 % ql at 0.5 g/kg (offline estimate) |
| 9 | Top boundary and sponge details | different by design | only in the top sponge cell |
| 10 | SGS floors: `numol` added (uDALES) vs `ekmin = 1e-6` floor (DALES) | different by design | negligible (1e-5 m2/s) |
| 11 | Random perturbation generators | different by design | different realisations of the same amplitudes; the ensembles absorb it |
| 12 | Smagorinsky: no stability correction in uDALES; no Vreman in DALES | different by design | not used by the validation cases |

## Governing equations, base state, Exner

| item | DALES dev 075ac4e | uDALES | status | effect |
|---|---|---|---|---|
| density | anelastic: `rhobf(k)`, `rhobh(k)` from the base state; advection `-(1/rhobf) d(rhobf w a)/dz` (DALES advec_2nd.f90:101, 127), Poisson RHS `rhobf (du/dx + dv/dy) + d(rhobh w)/dz` (DALES modpois.f90:318-320), SGS fluxes weighted by `rhobh/rhobf` (DALES modsubgrid.f90:712-716) | Boussinesq: `rhobf = rhobh = 1` (src/modfields.f90:556-557); no density in advection (src/modadvection.f90:103-146), Poisson (src/modpois.f90:962) or diffusion | different by design (see summary #1) | validation: DALES `ibas_prf = 5` with constant `rhobf` |
| base state | `ibas_prf = 3` default: standard lapse rate from `thls` (DALES modglobal.f90:133, modstartup.f90:1764-1774); 4.4.2 the same | none | different by design | DALES base `pb` set to the uDALES hydrostatic pressure of the initial profile (`make_dales_baseprof.py`) |
| Exner | fixed at start-up: `exnf = (pb/pref0)^(rd/cp)` (DALES modstartup.f90:1832); `presf` is recomputed hydrostatically every step (DALES modthermodynamics.f90:583-639) and used for qsat | recomputed every step from the hydrostatic slab-mean `presf` (src/modthermodynamics.f90:350-357, 363-367) | different by design (summary #2) | drift of order 1e-5 over the validation runs |
| buoyancy | `wp += grav (thv0h - thvh)/thvh`, `thvh` the slab mean of `thv0h`, first level from `th0av` (DALES modforces.f90:135, modthermodynamics.f90:288) | the same (src/modforces.f90:83, src/modthermodynamics.f90:90, 100) | aligned | - |
| hydrostatic pressure start | `thvh(1)` from `th0av(1)`, `qt0av(1)`, `ql0av(1)` (DALES modthermodynamics.f90:611) | `thvh(kb) = thvs` (src/modthermodynamics.f90:434) | different by design | less than 1 Pa at the first level |

## Grid and staggering

| item | DALES | uDALES | status |
|---|---|---|---|
| staggering | Arakawa C, periodic laterally | the same | aligned |
| vertical grid | `zf` from `prof.inp`, `zh(1) = 0`, `zh(k+1) = zh(k) + 2 (zf(k) - zh(k))`, `dzh(1) = 2 zf(1)` (DALES modglobal.f90:445-460) | the same (src/modglobal.f90:768-783) | aligned |
| filter width | `delta = (dx dy dzf)^(1/3)` (DALES modglobal.f90:468) | `(dxf dy dzf)^(1/3)` (src/modglobal.f90:816) | aligned for uniform dx |
| halos | 1 (cd2), 2 (kappa), 3 (5th/6th) (DALES modglobal.f90:330-366) | `ih = 1` for the prognostic fields; scalars and kappa-advected `thl` use a 2-cell copy (`ihc = 2`, `thl0c`) (src/modglobal.f90:606-630) | - |

## Advection

| item | DALES | uDALES | status |
|---|---|---|---|
| momentum | cd2, 5th order (5), 5th horizontal / 2nd vertical (52), cd6 (6), 6th / 2nd (62), hybrid 5th-order/WENO (55, 555) (DALES modadvection.f90:67-140) | cd2 only (src/modadvection.f90:46-54) | different, open |
| e12 | as momentum, plus kappa and upwind | cd2 only (src/modadvection.f90:56-64) | different, open |
| thl | cd2, 5, 52, cd6, 62, kappa (7), upwind (1), hybrid (DALES modadvection.f90:182-218) | cd2 or kappa (src/modadvection.f90:66-76) | different, open |
| qt | as thl (DALES modadvection.f90:223-) | cd2 only (src/modadvection.f90:78-86) | different, open |
| passive scalars | one scheme for all scalars (`iadv_sv`, a scalar in dev; an array in 4.4.2) | cd2, kappa or upwind in the dispatch (src/modadvection.f90:87-99), but `initglobal` overwrites `iadv_sv` with kappa for every scalar, whatever the namelist says (src/modglobal.f90:577-580) | different, open (#407) |
| default | `iadv_mom = 5`; the others default to `iadv_mom` (DALES modglobal.f90:142, 373-381) | `iadv_mom = 2`; the others default to `iadv_mom` (src/modglobal.f90:421, 565-567) | - |
| cd2 stencil | flux form `u (a_i + a_{i-1})/2`, vertical with `rhobf` weights (DALES advec_2nd.f90:36-139) | the same without density, general `dzf` weights (src/modadvection.f90:103-146) | aligned for constant density and uniform dz |
| kappa limiter | `phi(r) = max(0, min(2r, min(1/3 + 2r/3, 2)))` (DALES advec_kappa.f90:409-421; the advection loops use an inlined form without the division and `eps1`, advec_kappa.f90:26) | the same function, with `eps1` (src/modadvection.f90:376-387) | aligned (up to the `eps1` regularisation) |
| half levels for buoyancy | `thl0h`, `qt0h` from the kappa reconstruction when that variable is kappa-advected, linear otherwise (DALES modthermodynamics.f90:225-226, 932-973) | always linear (src/modthermodynamics.f90:551-581) | different, open (advection PR) |

**What the GCSS DALES cases use** (DALES `cases/*/namoptions.001` at 075ac4e):

| DALES case | mom | tke | thl | qt | sv | Galilean shift |
|---|---|---|---|---|---|---|
| bomex | 2 | 2 | 2 | 2 | 2 | `cu = -8` |
| atex | 2 | 2 | 2 | 2 | 7 | `cu = -6.5, cv = -1` |
| rico | 2 | 2 | 2 | 2 | 7 | `cu = -5, cv = -4` |
| dycoms_rf02 | 5 | 5 | 5 | 5 | 7 | `cu = 5, cv = -5.5` |
| smoke | 2 | 2 | 7 | 7 | 7 | none |

DALES has no RF01 case; its stratocumulus case (RF02) uses the 5th-order scheme for everything and kappa for
passive scalars. The cumulus cases use cd2 for momentum, thl and qt, with a Galilean shift that keeps the mean
wind relative to the grid small.

**What the validation used:** cd2 for u, v, w, e12, thl and qt in both codes (examples 301 and 303;
DALES inputs in `tools/python/examples/dales_cases/dales_inputs`). That is the DALES BOMEX choice, but not the
DALES stratocumulus choice. Where it matters:

- **RF01 inversion (5 m grid).** thl jumps by about 8.5 K and qt by about 7.5 g/kg across the inversion. cd2 is
  not monotone, so thl and qt overshoot and undershoot there, which feeds spurious mixing and cloud-top
  evaporation. How much this changes RF01 relative to DALES's 5th-order scheme has not been measured; it is an
  open item of the advection pull request.
- **BOMEX** (dz = 40 m, a smooth inversion) uses cd2 in the DALES case as well, so the validation is like for
  like with the DALES standard set-up, except for the Galilean shift.
- **Smoke (302)** uses kappa for thl and the smoke scalar in both codes; qt is inert there.

uDALES has kappa only for thl and the passive scalars, and only for thl through a dedicated copy (`thl0c`).
Extending kappa (or a 5th-order scheme) to qt and e12 is the planned advection pull request.

## Time integration and time-step control

| item | DALES | uDALES | status | effect |
|---|---|---|---|---|
| scheme | low-storage RK3 (Wicker and Skamarock), `rk3coef = dt/(4 - rk3step)` (DALES tstep.f90:259) | the same (src/modtstep.f90:191) | aligned | - |
| e12 floor | `e120 = max(e12min, ...)` every substep (DALES tstep.f90:272) | the same (src/modtstep.f90:215, 226) | aligned | - |
| CFL number | L2 norm `sqrt((u dt/dx)^2 + (v dt/dy)^2 + (w dt/dzh)^2)`, in the Galilean frame (DALES tstep.f90:190-206) | L1 sum `(abs(u)/dx + abs(v)/dy + abs(w)/dzh) dt`, physical winds (src/modtstep.f90:119, 139) | different by design (summary #5) | for the same `courant` uDALES takes shorter steps (up to sqrt(3) for equal components) |
| default `courant` | 1.0 for cd2 and 5th order, 0.7 with kappa, cd6 or 62 (DALES tstep.f90:63-98) | 1.5 for cd2 momentum, 1.4 otherwise, at most 1.1 with kappa or upwind (src/modglobal.f90:584-598) | - | - |
| diffusion limit | Peclet `max(ekm, ekh) dt / min(dzh, dx, dy)^2 <= peclet`, default 0.15 (DALES tstep.f90:194-198, modglobal.f90:184) | `max(ekm, ekh)(1/dzh^2 + 1/dx^2 + 1/dy^2) dt <= diffnr`, default 0.25 (src/modtstep.f90:120-121, src/modglobal.f90:482) | different by design | the cases convert `peclet` to `diffnr` by the ratio of the two formulas (case READMEs) |
| output times | adaptive mode: `dt` is clipped to hit radiation and output times (`dt_lim`, DALES tstep.f90:206, modradiation.f90:243-245) | not clipped | different by design | none for statistics sampled at multiples of a fixed dt |
| spin-up | `dt` doubles from a small value until the CFL and Peclet numbers settle, even with `ladaptive = .false.` | no spin-up by default (src/modtstep.f90:64) | different by design | first steps only (about 1 % e12 offset in the first minute of an RF01 column test) |
| logged Courant number | from absolute velocities, L2 per level (DALES modchecksim.f90:297-357) | `tcheck` logs `max((um/dx + vm/dy + wm/dz) dtmn)` with signed velocities and the interval-mean dt (src/modchecksim.f90:115), not the controlling L1 sum | different, open (#408) | diagnostic only; with a large mean wind of either sign the logged number can be far from the controlling one |

Because of these differences the same case does not take the same adaptive steps in both codes; the measured
steps are in the case RESULTS files. The DALES comparison therefore uses an identical fixed dt.

## Pressure solver

| item | DALES | uDALES | status |
|---|---|---|---|
| method | FFT in x and y, tridiagonal in z (`solver_id = 0` default; FFTW, cuFFT and HYPRE optional) (DALES modpois.f90:85-136, 411-) | FFT in x and y (2DECOMP&FFT), tridiagonal in z (`ipoiss = 0` default) (src/modpois.f90:64-415) | aligned for periodic domains |
| RHS | `rhobf (du/dx + dv/dy) + d(rhobh w)/dz` of `u + dt up` (DALES modpois.f90:318-320) | plain divergence (src/modpois.f90:962) | aligned for constant density |
| form | full pressure each substep | incremental: `pres0` gradient added in the advection routines, the solve gives the increment, `pres0 += p` (src/modadvection.f90:180, 224, 274; src/modpois.f90:1146) | aligned (the same projected velocity) |
| boundaries | w = 0 at bottom and top, periodic laterally | the same for the cases here | aligned |

## Subgrid-scale models

The three validation cases use the one-equation (Deardorff) TKE model with `ldelta = .false.`, as DALES.

| item | DALES | uDALES | status | effect |
|---|---|---|---|---|
| constants | `cm = cf/(2 pi) (1.5 alpha)^-1.5` (0.118), `ch = 1/Prandtl`, `ch2 = ch - ch1`, `ceps`, `ce1 = cn^2 (cm/Rigc - ch1 cm)`, `ce2 = ceps - ce1`; cf 2.5, cn 0.76, Rigc 0.25, Prandtl 1/3, ch1 1 (DALES modsubgrid.f90:71-79, modsubgriddata.f90:42-54) | the same (src/modsubgrid.f90:75-83, src/modsubgriddata.f90:43-55) | aligned | - |
| length scale | `zlt = delta`, or `min(delta, cn e12 / N)` where `dthvdz > 0` (DALES modsubgrid.f90:462-465) | the same, with `max(abs(dthvdz), eps1)` under the root (src/modsubgrid.f90:398-406) | aligned | - |
| Km, Kh | `ekm = cm zlt e12`, `ekh = (ch1 + ch2 zlt/delta) ekm` (DALES modsubgrid.f90:469-470) | the same (src/modsubgrid.f90:400-401, 407-408) | aligned | - |
| floors | `max(ekm, ekmin)`, `max(ekh, ekmin)`, `ekmin = 1e-6` (DALES modsubgrid.f90:472-473) | `+ numol` (1.5e-5) and `+ numol/0.71`, no floor (src/modsubgrid.f90:402-403, 409-410) | different by design | negligible (1e-5 m2/s in the free troposphere) |
| TKE sources | shear `Km S^2/(2 e12)`, buoyancy `-Kh (g/thvf) dthvdz/(2 e12)`, dissipation `-(ce1 + ce2 zlt/delta) e12^2/(2 zlt)` (DALES modsubgrid.f90:579-581) | the same, molecular parts subtracted, near-wall `damp` (1 away from walls) (src/modsubgrid.f90:502-509) | aligned | - |
| first level | MO shears `dudz`, `dvdz` and `dthvdz(1)` in the k = 1 block (`sgs_surface_fix = .false.`) (DALES modsubgrid.f90:648-678) | the same in the fluid cells that own a c-grid section of a horizontal Monin-Obukhov facet (at a roof edge this can be a cell beside the roof) (`iwallmom = 4`, facet normal +z to within 1e-6; a start-up line counts the cells): `du/dz`, `dv/dz` from a per-cell Monin-Obukhov state on the cell-centre `thl`, `qt` and wind (the DALES per-column surface; not mapped from the u/v section states, which would need section data from the neighbouring rank at tile edges), `dthvdz` from the surface-layer `dthl/dz`, `dqt/dz` (src/modibm.f90:mo_cells_update, src/modsubgrid.f90:523-556, src/modthermodynamics.f90:calthv). Next to sloped Monin-Obukhov facets and every other wall: no sources at `kb`, `dthvdz(kb) = 0` | aligned (horizontal Monin-Obukhov facets); different by design elsewhere | BOMEX: resolved w'thl' at 40 m 1.79e-3 against DALES 1.82e-3 K m/s (3.47e-3 without the treatment); RF01: at 15 m 5.34e-3 against 5.36e-3 K m/s (7.86e-3 without); all mean and flux profiles within the spread in both cases; RF01 cloud cover moves from 1.2 % below to 1.2 % above DALES (examples/301, 303 RESULTS.md) |
| SGS fluxes | `rhobh/rhobf`-weighted diffusion; surface flux enters as `flux/dzf` (DALES modsubgrid.f90:701-737) | the same stencil without density; the facet flux replaces the diffusion across the floor (`modibm:diffu_corr`, `diffc_corr`, `wallfunmom`, `wallfunheat`) | aligned for constant density | - |
| e12 halos | exchanged with the other prognostic fields | exchanged when `loneeqn` (src/modboundary.f90:91-94; #405) | aligned | - |
| wall damping | none (optional Mason length scale, `lmason`) | `damp = 1` (Piomelli damping commented out) (src/modsubgrid.f90:385-390) | aligned in effect | - |
| Smagorinsky | `ekm = (cs delta)^2 sqrt(2 S^2) sqrt(1 - Ri/Pr)`, MO `dudz` at k = 1 (DALES modsubgrid.f90:262-363) | neutral `(cs delta damp)^2 sqrt(2 S^2)`, no stability correction (src/modsubgrid.f90:219-279) | different by design | not used by the validation cases |
| Vreman | not available | Vreman (2004), optional buoyancy correction (`lbuoycorr`) (src/modsubgrid.f90:281-371) | uDALES only | - |
| top | `ekm(k1) = ekm(kmax)` | the same | aligned | - |

## Surface schemes

| item | DALES | uDALES | status | effect |
|---|---|---|---|---|
| where | the flat bottom, per column (`isurf`) | facets: `iwallmom = 4` on facets within `mo_maxslope` (30 degrees) of horizontal, per facet section of the u, v and w grids; the floor of 301-303 is one horizontal facet (src/modibm.f90:289-317, 1821-1884) | different by design | - |
| MO with prescribed u* | `isurf = 3` (DALES surface/modsurface.f90:478-485) | `iwallmom = 4` with `ustin > 0` in `&WALLS` (src/modsurface.f90:39-59) | aligned | BOMEX |
| MO with diagnosed u* | `isurf = 4` (DALES surface/modsurface.f90:486-493) | `iwallmom = 4` (src/modsurface.f90:39-59) | aligned | RF01, smoke |
| MO with prescribed surface temperature | `isurf = 2`: `calc_mean_wind`, `get_obl` (L from Rib between the first level and `tskin = thls`, `qskin` of the previous step), `calc_drag_coefficients` with that new L, `ra = 1/(Cs horv)`, `tskin = thls`, `qtsurf`, `u* = sqrt(Cm) horv`, `thlflux = (tskin - thl)/ra` (DALES surface/modsurface.f90:468-477, 507-555, 558-573, 724-797, 905-958) | `iwallmom = 4` with `iwalltemp = 2`: the same sequence per facet section with `tskin` = the facet temperature `facT(fac,1)`; `w'thl' = Cs U (Tfacet - thl)` (src/modsurface.f90:72-91, src/modibm.f90:1615-1662); applied by `wallfunheat` on the c-grid sections (src/modibm.f90:2065-2085), the stress by `wallfunmom` (src/modibm.f90:1860-1865) | aligned | example 304 (RESULTS.md) |
| surface humidity with a prescribed surface temperature | `qtsurf`: `qskin = surfwet qsat(tskin) + (1 - surfwet) qt`, `surfwet = ra/(ra + rsisurf2)` (DALES surface/modsurface.f90:867-902) | facets with `lGR = 1`: `bcqfz` (`iwallmoist = 1`) or the uDALES resistance law `moist_flux` with `resa = 1/(Cs U)` (vegetation fraction 0.8, canopy and soil resistances `facf`, `fachurel`, Bolton `qsat`) (`iwallmoist = 2`); `qskin = qt + w'qt'/(Cs U)` from the previous flux, as DALES uses the `qskin` of its previous `qtsurf`; other facets `qskin = qt`; the diagnosed `qskin` is bounded to `qt +- 5e-2` (the `calc_surface_scalars` bound of isurf 3/4, which DALES isurf 2 does not apply) (src/modibm.f90:1615-1662) | different by design (uDALES facet moisture laws) | 304 is dry |
| `ustin` with a prescribed surface temperature | not available (`isurf = 2` diagnoses u*) | `ustin > 0` prescribes u* also under `iwalltemp = 2` (src/modsurface.f90:85-89) | uDALES extension | - |
| cold start, prescribed surface temperature | `tskin = thls`, `obl = -0.1`, `qskin` from `qtsurf` with `ra = 50` (DALES modstartup.f90:833-834, 855-862) | per section at its first visit: `L = -0.1`, `qskin = qt` (zero flux), one extra update (src/modibm.f90:1626-1633) | aligned (except the DALES `ra = 50` start of `qskin`) | - |
| steep facets with a facet temperature | - | `iwallmom = 4`, `iwalltemp = 2`: facets steeper than `mo_maxslope` use `wfuno` with its stability correction (the `iwallmom = 2` path; src/modibm.f90:1797); neutral under `iwalltemp = 1` | uDALES only | - |
| stability functions, Obukhov length iteration | unstable Businger-Dyer, stable Beljaars-Holtslag; `calc_obl_iter` Newton iteration on Rib (DALES surface/modsurface.f90) | the same formulas (src/modsurface.f90:94-130, 140-) | aligned | - |
| local or domain-mean MO | dev default `lmostlocal = .true.`: per-column L, Cm and (isurf 4) u* = sqrt(Cm) times the local wind speed (DALES surface/modsurfdata.f90:240, surface/modsurface.f90:651-690, 905-958); 4.4.2 default `.false.` | local: per facet section, L, Cm, Cs, u*, tskin and qskin from the air state interpolated to the section's point (the u, v points on a flat floor) (src/modsurface.f90:39-59) | aligned (summary #7) | RF01 u* 0.2464 vs 0.2469 m/s (0.2 %); BOMEX u* prescribed |
| surface fluxes | `wtsurf`, `wqsurf`, positive upward | facet fluxes `bctfz`, `bcqfz` (on facets of a type with `lGR = 1`), positive downward | aligned (sign converted) | - |
| cold start | `tskin = thlprof(1)`, `qskin = qtprof(1)`, `obl = -0.1`, one `surface` call in `readinitfiles` (DALES modstartup.f90:836-889) | the same per section at its first visit (src/modibm.f90:1873-1880); the first-level cells at a cold start in `readinitfiles` (src/modibm.f90:1680-1742) | aligned | - |
| warm start | the restart file holds `ustar`, `thlflux`, `qtflux`, `dthldz`, `dqtdz`, `obl`, `tskin`, `qskin`, and `dthvdz`, `ekh`, `presf`, `presh` (DALES modstartup.f90:1413-1433) | per-rank `initm` files with the per-section L, `tskin`, `qskin`, u*, Cs, csu, fluxes and the first-level-cell gradients (src/modibm.f90:468-534, read in src/modstartup.f90:2472-2494); `presf`, `presh`, `exnf`, `rhof`, `thv0h`, `dthvdz` re-diagnosed from the restored fields instead of stored, and the cold-start values of the `thl0`, `qt0` ghost level `kb - 1` restored (src/modthermodynamics.f90:63, src/modstartup.f90:1877-1896); without `initm` the cold start | aligned (bit-identical continuation; with the one-equation closure up to round-off: the e12 ghost level `kb - 1` is not stored, `ekh` is not stored either but is recomputed before use) | the re-diagnosis runs only with `iwallmom = 4` or `iradiation /= 0`; other warm starts keep their approximate `thv` profiles |
| von Karman constant | 0.4 (parameter, DALES modglobal.f90:101) | 0.40 (`&WALLS fkar`, src/modglobal.f90:338); the inlet generator keeps 0.41 with its Coles wake constants | aligned | 0.41 against 0.4 gives a 2.3 % u* difference in RF01 |
| order | `surface` before advection and subgrid (DALES program.f90:232-238) | the stress in the wall function (`ibmwallfun`, after subgrid), from the same `u0`, `thl0`, `qt0` of the substep; the first-level cells above horizontal Monin-Obukhov facets (`mo_cells_update`) before subgrid, and at a cold start once before the start-up thermodynamics, as DALES `readinitfiles` (src/program.f90:166, src/modstartup.f90:1689-1696) | aligned | - |
| stress stencil | cell-centred `u + cu`, 4-point v, u* averaged over two cells, wind floor 1e-10 (DALES modsubgrid.f90:942-951) | u* and L of the section at the u (v) point, 4-point v (u), `min_horv = 0.1` in Cm (src/modibm.f90:1860-1884) | different by design | second-order equivalent; on the saved first-level fields of example 304 the uDALES stencil gives a 0.04-0.06 % larger mean stress (examples/304/RESULTS.md) |

### Wall functions

uDALES's flat-bottom conditions and every facet outside `iwallmom = 4` use wall functions instead of Monin-Obukhov
similarity. DALES has no wall-function bottom conditions; its grid-conforming immersed-boundary method
(`lapply_ibm`) applies a neutral log-law drag on vertical walls and a Monin-Obukhov drag on horizontal surfaces of
obstacles (DALES modibm.f90:183-193, 640-649), without a facet energy balance. The comparison cases use neither.

| item | DALES | uDALES | status |
|---|---|---|---|
| flat bottom | MO only (`isurf`) | `BCbotm = 2`: `wfuno`, a log-law stress and heat flux with a Richardson-number stability correction (Uno et al. 1995, Cai 2012; extension of Louis 1979 to rough walls); `BCbotm = 3`: `wfmneutral`, neutral log law; for MO use a facet floor with `iwallmom = 4` (above) (src/modibm.f90:bottom, src/modwallfunctions.f90:24, 262) | different by design |
| building facets | `lapply_ibm`: constant neutral law-of-the-wall coefficients for vertical walls, MO drag with stability for horizontal faces; prescribed wall temperatures (`lwallheat`) | the same wall functions at every IBM facet, `iwallmom = 2` (stability) or 3 (neutral), with `iwalltemp`/`iwallmoist` for scalars; with `iwallmom = 4` MO on the near-horizontal facets and on the others the stability-corrected wall function (`iwalltemp = 2`) or the neutral law (`iwalltemp = 1`) (src/modibm.f90:wallfunmom, wallfunheat) | different by design |
| TKE at the first level | MO sources at k = 1 | MO sources above horizontal Monin-Obukhov facets (`iwallmom = 4`, see Subgrid model); no sources at `kb` next to the wall functions and sloped Monin-Obukhov facets (src/modsubgrid.f90:515-556) | aligned on horizontal Monin-Obukhov facets; different by design elsewhere |

## Thermodynamics

| item | DALES | uDALES | status | effect |
|---|---|---|---|---|
| constants | grav 9.81, rd 287.04, rv 461.5, cp 1004, rlv 2.53e6, pref0 1e5, tmelt 273.16 (DALES modglobal.f90:78-95) | the same (src/modglobal.f90:287-331) | aligned | `rlv = 2.53e6` is the DALES value, chosen for identical constants; the real value is 2.50e6 at 0 °C and 2.45e6 at 20 °C (2.53e6 is about 3 % above it), and the old uDALES 2.26e6 is the 100 °C value |
| es | Murphy and Koop (2005) table, 0.2 K spacing, linear interpolation; liquid-ice blend between 253 and 268 K (DALES modthermodynamics.f90:136-156, 701-716) | Tetens, `610.78 exp(17.27 (T - 273.16)/(T - 35.86))` (src/modthermodynamics.f90:507, 530); no ice | aligned in the validation (DALES patched to Tetens) | Tetens is 0.10-0.13 % lower at 285-298 K, about +1 % ql at 0.5 g/kg (offline estimate) |
| qsat | `(rd/rv) es / (p - (1 - rd/rv) es)` | the same (src/modthermodynamics.f90:508) | aligned | - |
| saturation adjustment | two linearised steps, the second re-linearised at the first estimate (DALES modthermodynamics.f90:766-797); its `lqlnr` is ignored. 4.4.2: Newton-Raphson iteration (`icethermo0`) by default | `lqlnr = .true.`: Newton-Raphson to 1e-5 relative (src/modthermodynamics.f90:491-520), used by 301-303; default `.false.`: one linearised step (src/modthermodynamics.f90:521-544) | different by design (summary #3) | dev two-step within about 0.2 % of exact; the one-step default overestimates ql by 0.6-9 % in BOMEX cores (offline estimate) |
| half levels | linear (or kappa, see Advection) with `thls`, `qts` at the surface; saturation adjustment repeated at half levels | linear with `thls`, `qts` (src/modthermodynamics.f90:551-581); repeated at half levels (src/modthermodynamics.f90:81-83); with `iwallmom = 4` the namelist `thls`, `qts` are kept (DALES re-diagnoses them every step from `tskin`, `qskin`; `bldump` reports the diagnosed means) | aligned for cd2 | - |
| `calthv` | Deardorff moist/dry `dthvdz` with `chi_half = 0.5`, MO gradient at k = 1 (DALES modthermodynamics.f90:333-479) | the same, plus a guard against a zero denominator (src/modthermodynamics.f90:140-274) | aligned | - |
| density | `rhof = presf/(rd Tv)` (DALES modthermodynamics.f90:575) | the same (src/modthermodynamics.f90:386) | aligned | - |

## Radiation

| item | DALES | uDALES | status | effect |
|---|---|---|---|---|
| parametrised longwave | `radpar`: `lwd = dlwtop exp(-rka LWP_above)`, `lwu = dlwbot exp(-rka LWP_below)`, `thlprad = -dF/dz /(rhof cp exnf)` (DALES modradiation.f90:339-452) | `iradiation = 2`, the same (src/modradiation.f90:135-171); absorber ql or a scalar (`rad_absorber`) | aligned | - |
| update | at `rk3step = 1` every `timerad`, held for all substeps (DALES modradiation.f90:243-248) | the same (src/modradiation.f90:122-130) | aligned | - |
| warm start | `tnext_radiation` and `thlprad` in the restart file (DALES modstartup.f90:1436-1437) | not stored: updated at the first step of a warm start, then on the `timerad` grid of the continued run (src/modradiation.f90:122-126) | different by design | with `timerad > 0` and a checkpoint between updates, one interval differs |
| DYCOMS free-troposphere term | `rad_user` in `cases/dycoms_rf02/moduser.f90`: pointwise derivative, sign that warms (dalesteam/dales#347) | case hook `dycoms_rf01`: divergence of the analytic flux between faces, cooling (src/cases/case_dycoms_rf01.f90:112-149) | different by design; the DALES reference was patched to the uDALES form | - |
| `iradiation = 10` | `initradiation` resets `rad_longw = rad_shortw = rad_smoke = .false.` for `iradiation = 10`, so `radpar` never runs in the `irad_user` branch (DALES modradiation.f90:187; also in 4.4.2) | not applicable | DALES defect (dalesteam/dales#347); the DALES comparison inputs patch it | without the patch DALES RF01 has no cloud-top longwave cooling |
| prescribed cooling | `thlpcar` added through `rad_ls` (DALES modradiation.f90:570) | added in `forces` (src/modforces.f90:107) | aligned | BOMEX -2 K/day |

## Large-scale forcing and subsidence

| item | DALES | uDALES | status |
|---|---|---|---|
| subsidence | upwind `-whls dphi/dz` on the local fields, ascent term not at k = 1 (DALES modforces.f90:339-345) | the same (src/modforces.f90:766-783) | aligned |
| momentum subsidence | `lmomsubs`, default off | the same | aligned |
| prescribed tendencies | `dqtdtls`, advective `dthldxls` etc. (DALES modforces.f90:381) | the same (src/modforces.f90:782-785) | aligned |
| geostrophic forcing | `dpdxl = om23_gs vg`, `dpdyl = -om23_gs ug` | the same (src/modstartup.f90:2243-2244) | aligned |
| Coriolis | 4-point averages with the `2 Omega cos(lat)` terms (`om22`) in u and w, `omega = 7.292e-5`, plus `cu`, `cv` terms (DALES modforces.f90:200-286) | the same without `cu`, `cv` (src/modforces.f90:631-675) | aligned for `cu = cv = 0` |

## Galilean transform, sponge layer and boundary conditions

| item | DALES | uDALES | status | effect |
|---|---|---|---|---|
| Galilean transform | `cu`, `cv`: fields stored relative to the moving frame, correction terms in Coriolis, damping and the surface wind (DALES modglobal.f90:162-163) | none | different by design (summary #4) | not measured; DALES BOMEX uses `cu = -8` |
| sponge | `igrw_damp = 2` default: u, v relax to `ug - cu`, `vg - cv`, w to 0, thl, qt to slab means; `ksp = min(3 kmax/4, kmax - 15)`, `tsc = rnu0 sin^2(...)`, `rnu0 = 2.75e-3` (DALES modglobal.f90:117, modboundary.f90:81-92, 261-268) | `igrw_damp = 0` default; option 2 is the same with `cu = 0` (src/modboundary.f90:47-59, 1485-1492) | aligned in the validation (both 2) | BOMEX from zf = 2380 m, RF01 from 1122.5 m |
| kmax pinning | `thl0(kmax)`, `qt0(kmax)` set to the slab mean every substep (DALES modboundary.f90:305-308) | none | different by design | inside the sponge only |
| top boundary | `e12(k1) = e12min`; thl, qt extrapolated with the top gradient (DALES modboundary.f90:340-420) | free slip, zero gradient for thl, qt and e12 (src/modboundary.f90:175-190) | different by design | only in the top cell (inside the sponge); RF01 top cell cools by about 0.2 K over 4 h |
| lateral | periodic | periodic for these cases; inflow, outflow and driver options | aligned for these cases | - |

## Random perturbations

| item | DALES | uDALES | status |
|---|---|---|---|
| generator | sequential LCG over the global domain, carried across levels and fields (DALES modstartup.f90:1623-1662) | counter-based splitmix64 hash of `(irandom, field, i, j, k)` (src/modstartup.f90:2653-2701) | different by design (other realisations) |
| fields and levels | thl and qt on levels `1..krand` (default all), m and 0 fields drawn separately; u, v, w on `krandumin..krandumax` (default none) (DALES modstartup.f90:763-786) | u, v, w on `kb..krand` with `randu`; thl, qt with `randthl`, `randqt` (one draw for the m and 0 fields) (src/modstartup.f90:1295-1319) | aligned in the validation (no velocity noise, thl 0.1 K and qt 2.5e-5 on all levels) |
| qt cap | amplitude capped at the local qt | the same | aligned |
| defaults | `randthl = 0.1`, `randqt = 1e-5`, `randu = 0.5`, `irandom = 0` (DALES modstartup.f90:44-47) | `randu = 0.01`, `randthl = randqt = 0`, `irandom = 43` (src/modstartup.f90:42-44) | different |

## Diagnostics

| item | DALES | uDALES bldump | status |
|---|---|---|---|
| LWP | mean over columns of `sum rhof ql dzf` (DALES modtimestat.f90:424-435, 617) | the same (src/modbldiag.f90:115-119) | aligned |
| cloud cover | fraction of columns with column ql > 0 (DALES modtimestat.f90:432, 616) | fraction of columns with any ql > 0 (src/modbldiag.f90:120) | aligned |
| zi | max gradient of thv on 4-column composites, refined by the zero of the gradient change (DALES modtimestat.f90:60-61, 1191-1221) | per column, half level of the max dthl/dz, then averaged (src/modbldiag.f90:121-129) | different; the validation computes zi identically for both codes from the 600 s block-mean thl profile (`compare_dales.py`) |
| zi from qt | no timestat output (the DYCOMS `rad_user` uses a qt = 8 g/kg crossing per column for the radiation only) | `zi_qt` with `zi_qt_thres` (src/modbldiag.f90:130) | the validation computes `zi_qt8` for both codes from the block-mean qt profile (`compare_dales.py`) |
| sampling | after `thermodynamics`; profiles are 600 s block means of 60 s samples | after `thermodynamics` (src/program.f90:240-243); 60 s snapshots, block-averaged in `compare_dales.py` | aligned |

## Settings of the validation runs

The full side-by-side settings of the comparison runs are in
[examples/301/RESULTS.md](https://github.com/uDALES/u-dales/blob/master/examples/301/RESULTS.md) and
[examples/303/RESULTS.md](https://github.com/uDALES/u-dales/blob/master/examples/303/RESULTS.md); the input
sources are in each case README. The tables below give the reason for each numerics setting and flag the ones
known to matter. "DALES case" means the DALES `cases/<name>/namoptions.001`; "match" means chosen so that both
codes run the same numerics.

### 301 BOMEX

| setting | uDALES 301 | DALES reference (`dales_inputs/bomex`) | reason |
|---|---|---|---|
| grid | 64 x 64 x 80, 100 m x 100 m x 40 m | the same | DALES case |
| advection | cd2 for u, v, w, e12, thl, qt | cd2 for all | DALES case |
| time step | fixed 2.0 s (commented block in `namoptions.301`); adaptive by default (`courant = 0.55`, `diffnr = 0.2`, `dtmax = 20`) | fixed 2.0 s (10800 steps); adaptive: `courant` 1.0 (L2), `peclet` 0.15 | match; 2.0 s is below both adaptive steps |
| SGS | one-equation TKE, `ldelta = .false.`, DALES constants | the same | DALES case |
| surface | facet floor, `iwallmom = 4`, `ustin = 0.28`, `bctfz`/`bcqfz` = 8e-3 and 5.2e-5 upward, z0 = z0h = 2e-4 | `isurf = 3`, the same values; `lmostlocal` default (`.true.`) | DALES case; local MO in both |
| von Karman | 0.40 (default) | 0.4 | DALES value; u* is prescribed, so only the surface-layer gradients see it (the value used in the comparison runs is recorded in RESULTS) |
| saturation | `lqlnr = .true.` (exact) | two-step | closest uDALES option |
| radiation | none (`iradiation = 0`, `timerad` unused); -2 K/day `thlpcar` | `rad_ls` | DALES case |
| forcing, subsidence | DALES `lscale.inp`, local upwind subsidence | the same | DALES case |
| Coriolis | `xlat = 15` | 15 | DALES case |
| sponge | `igrw_damp = 2`, `ksp = 60` (zf = 2380 m), `rnu0 = 2.75e-3` | the same | DALES default |
| Galilean shift | none | `cu = cv = 0` | uDALES constraint (DALES case: `cu = -8`) |
| noise | thl 0.1 K, qt 2.5e-5 on all 80 levels, no u, v, w noise, seeds 43-45 | the same | DALES case |
| base state | Boussinesq | `ibas_prf = 5`, constant `rhobf`, `pb` from the uDALES hydrostatic profile | uDALES constraint |
| run, statistics | 6 h; window 3-6 h, 600 s blocks | the same | 6 h is enough for a quasi-steady 3-6 h window (the DALES file runs 8 h) |

Known to matter at this resolution: the Galilean shift (DALES runs BOMEX at `cu = -8`; not tested) and the
anelastic base state (not tested). The adaptive steps differ (see [Time integration](#time-integration-and-time-step-control)).

### 303 DYCOMS-II RF01

| setting | uDALES 303 | DALES reference (`dales_inputs/rf01`) | reason |
|---|---|---|---|
| grid | 96 x 96 x 300, 35 m x 35 m x 5 m | the same | Stevens et al. (2005) |
| advection | cd2 for u, v, w, e12, thl, qt | cd2 for all | match; DALES's own stratocumulus case (RF02) uses 5th order |
| time step | fixed 0.3 s (commented block in `namoptions.303`); adaptive by default (`courant = 0.5`, `diffnr = 0.16`, `dtmax = 20`) | fixed 0.3 s (48000 steps); `dtmax = 0.3` | match; below both adaptive steps |
| SGS | one-equation TKE, DALES constants | the same | DALES default |
| surface | facet floor, `iwallmom = 4`, fixed SHF 15 and LHF 115 W/m2 (`bctfz`, `bcqfz`), z0 = z0h giving neutral Cm = 0.0011 | `isurf = 4`, the same fluxes and z0; `lmostlocal` default (`.true.`) | Stevens et al. (2005), cross-checked against PyCLES; local MO in both |
| von Karman | 0.40 (default), `z0 = z0h = 1.445944e-5` m for neutral Cm = 0.0011 | 0.4, the same z0 | match (the value used in the comparison runs is recorded in RESULTS) |
| saturation | `lqlnr = .true.` | two-step | closest uDALES option |
| radiation | `iradiation = 2`, rka 85, F0 70, F1 22 W/m2, `timerad = 0` (every step), case hook `dycoms_rf01` (D = 3.75e-6) | `iradiation = 10` with radpar made reachable by the patch, patched `rad_user`, `timerad = 0` | Stevens et al. (2005) via PyCLES; the DALES patches give the same formula |
| subsidence, Coriolis | `wfls = -3.75e-6 z`, local upwind; `xlat = 31.5`, ug = (7, -5.5) | the same | Stevens et al. (2005) via PyCLES |
| sponge | `igrw_damp = 2`, `ksp = 225` (zf = 1122.5 m), `rnu0 = 2.75e-3` (time scale 6 min at the top) | the same | DALES default |
| Galilean shift | none | none | uDALES constraint (DALES RF02 case uses `cu = 5, cv = -5.5`) |
| noise | thl 0.1 K, qt 2.5e-5 on all 300 levels, no velocity noise, seeds 43-45 | the same (`krandumax = 0`) | DALES RF02 case |
| base state | Boussinesq | `ibas_prf = 5`, constant `rhobf` | uDALES constraint |
| run, statistics | 4 h; window 3-4 h, 600 s blocks; entrainment fit 1-4 h | the same | Stevens et al. (2005) run 4 h |

Known to matter at this resolution: **cd2 for qt and thl at the 5 m inversion** (see
[Advection](#advection)); the cd2 vs 5th-order comparison with DALES is open.

### 304 dry CBL over a floor of fixed temperature

uDALES 304 is the counterpart of DALES `isurf = 2` (not a DALES case; settings table in `examples/304/README.md`):
96 x 96 x 80 at 50 m x 50 m x 25 m, cd2 everywhere, fixed dt 1.0 s, one-equation SGS, `iwallmom = 4` with
`iwalltemp = 2` and the facet at 302 K (`lEB` off) against DALES `isurf = 2` with `thls = 302`, z0 = 0.1 m, z0h = 0.01 m,
geostrophic wind (5, 0) m/s at 45 degrees, thl noise 0.1 K, seeds 43-45, constant-density DALES base state. uDALES
is dry; DALES runs `lmoist` with qt = 0 and `rsisurf2 = 1e10` (dry `isurf = 2` leaves `qskin` uninitialised). Result
with 6 seeds (examples/304/RESULTS.md):
- The surface heat flux, zi and the thl and resolved w'thl' profiles are within the spread.
- u* is 0.27 % low (1.4 of the smaller, 0.85 of the larger spread). The transfer per unit first-level wind agrees (u*/|U1|
  within 0.03 %). The first-level wind is 0.3 % lower.
- A prescribed-flux control (`iwalltemp = 1` against `isurf = 4`) shows the same lower wind and u* (0.48 %, 0.43 %), so
  that difference does not come from the facet-temperature path. Its cause is not identified. The u-point stress stencil
  (row "stress stencil") accounts for about 0.05 % of the stress on the saved fields.
- The resolved w'thl' at 25 m is 1.9 % low (2.1 spreads). This does not appear in the control, and its cause is not
  identified.

### 302 smoke cloud

uDALES 302 follows DALES `cases/smoke`: 64 x 64 x 250 at 50 m x 50 m x 5 m, cd2 for momentum and e12, kappa for
thl and the smoke scalar (qt is cd2 in uDALES and kappa in DALES, but inert: `lmoist = .false.`), adaptive
`courant = 0.5` (L2 in DALES, L1 in uDALES), `diffnr = 0.16` for `peclet = 0.15`, `dtmax = 20`, MO with diagnosed u*
(z0 = 0.1 m), smoke longwave radiation (rka 0.02, F0 60 W/m2) every step instead of every 600 s, no Coriolis, no
sponge (`igrw_damp = 0`), thl noise 0.1 K on all levels. `xlat` is not set (uDALES default 52 degrees), which has no effect
because `lcoriol = .false.`. No DALES smoke reference run exists.

## Changes that alter results

The changes that brought uDALES closer to DALES (SGS model, constants, surface timing, subsidence, perturbations,
the e12 halo exchange) alter the results of existing runs. Who is affected, and by how much, is recorded in the
[CHANGELOG](udales-changelog.md) under "Breaking changes", which is the single source for that information.

## Open items

1. **Advection** (planned pull request): kappa or 5th-order advection for qt and e12 as well as thl, kappa half
   levels for the buoyancy, and the RF01 comparison of cd2 against DALES's 5th-order scheme (no result yet).
2. **Local MO in DALES dev** (`lmostlocal = .true.`): rerun the DALES references with `lmostlocal = .false.`, or
   port the local variant, to remove the last unmatched surface item.
3. **Galilean shift and anelastic base state**: their effect on BOMEX has not been measured (DALES runs with
   `cu = -8` and with `ibas_prf = 3` against the present references would show it).
4. **Adaptive time step**: retune `courant`/`diffnr` in 301 and 303 if a matched adaptive pair is wanted.
5. **BOMEX zi**: half a grid cell above DALES, not significant with three seeds; more seeds would resolve it.
6. **Scalar advection override** (#407) and the **logged Courant number** (#408).
