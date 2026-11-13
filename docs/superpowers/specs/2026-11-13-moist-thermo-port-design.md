# Moist thermodynamics and parametrised radiation from DALES — design

Date: 2026-11-13
Branch: `thermo` (from `origin/master` @ `752430c0`)
Source: DALES `dev` @ `075ac4e` (2026-11-12), https://github.com/dalesteam/dales

## Goal

Bring the DALES Deardorff moist-buoyancy scheme and the parametrised
(longwave) radiation into uDALES, and validate them on three standard
DALES/GCSS cases: BOMEX, the smoke cloud, and DYCOMS-II RF01.

## Scope

In:

1. Fix for issue #375 (`thls`/`qts` sentinels give a negative `thvs`).
1b. Monin–Obukhov surface layer for the flat bottom (DALES `isurf = 3/4`).
2. Deardorff all-or-nothing moist buoyancy (`calthv`) brought in line with DALES `dev`.
3. Parametrised longwave radiation (`iradiation = 2`), including the smoke
   absorber option and the DYCOMS free-troposphere term.
4. A diagnostic statistics module for cloud/radiation quantities.
5. CI regression tests and three full-size validation cases.

Out (deferred):

- Bulk microphysics (`modbulkmicro`, KK00/SB). DYCOMS RF02 is therefore
  replaced by the non-drizzling RF01.
- Shortwave radiation (`rad_shortw`): no target case uses it, and it needs
  solar geometry, surface albedo and skin temperature that uDALES does not
  provide in this form. Setting it aborts with an error.
- Full/RRTMG/land-surface radiation (`iradiation = 1, 3, 4, ...`).
- Integration with the statsout (v3.0) statistics rewrite.

## Design principles

- Port faithfully from DALES: same algorithms, namelist group and variable
  names, so a line-by-line diff against DALES stays readable. Translate only
  the indexing (`2:i1, 1:kmax` → `ib:ie, kb:ke`) and the module wiring.
- Everything new is off by default. With `iradiation = 0` and
  `lthermostat = .false.`, results are bit-for-bit identical to the parent
  commit (checked by the existing `david_tests` regression after each commit
  from step 2 onward).
- Correct in the presence of immersed boundaries even though the target
  cases are flat: tendencies and averages respect `IIc`.

## Commit sequence

Each commit builds, passes the supported test selection, and can be reviewed alone.

1. **#375 fix.** Abort in `checkinitvalues` when `thls`/`qts` are still at
   their sentinel while `ltempeq` or `lbuoyancy` is on; delete the dead
   commented call in `modstartup`; set `thls`/`qts` in `tests/cases/526` and
   any other shipped case that omits them; guard the length-scale `sqrt` in
   `modsubgrid`. This commit intentionally changes results for misconfigured
   cases; it may be split into its own PR if preferred.
2. **Diagnostics** (`modthermostat`) first, because later tests read it.
3. **Monin–Obukhov surface** (`modsurface.f90`), see Component 0.
4. **Deardorff `calthv`.** See below.
5. **Radiation.** New `modraddata.f90`, `modradiation.f90`.
6. **Tests and cases.** CI tests, `examples/301`–`303`, converter and
   comparison scripts, docs.

## Component 0: Monin–Obukhov surface layer (`modsurface.f90`)

Requested so that uDALES and DALES runs share the same surface treatment.
New `BCbotm` options on the flat bottom (`lbottom = .true.`):

- `BCbotm = 4` (`BCbotm_mo`, DALES `isurf = 4`): fluxes prescribed by
  `wtsurf`/`wqsurf`; `u*` from the Obukhov length, iterated from the bulk
  Richardson number (`calc_obl_iter`), with Businger–Dyer (unstable) and
  Beljaars–Holtslag (stable) `psim`/`psih`, roughness `z0`/`z0h`.
- `BCbotm = 5` (`BCbotm_mo_ustar`, DALES `isurf = 3`): as 4, but
  `u* = ustin` (new `&BC` parameter).

Both use DALES's domain-mean formulation (`lmostlocal = .false.`):
domain-mean `obl` and `u*`, and stress applied locally as
`u'w' = -u*² u/|U|` with `|U| ≥ 0.1 m/s`. The momentum update follows the
existing `wfmneutral` pattern (add back the resolved diffusive flux through
the bottom face, subtract the wall flux). Surface gradients
`dthldz = -wθ/u* · φh/(κ z1)` and `dqtdz` likewise feed Component 1. Surface
values `tskin`/`qskin` are diagnosed per column and averaged into
`thls`/`qts`, and `thvs` is updated, as in DALES.

**Flux sign convention (decision D1).** In `bottom`, `wtsurf` enters with
a minus sign and `wqsurf` with a plus sign, so a positive `wtsurf` is a
downward (cooling) flux and a positive `wqsurf` an upward (moistening) one.
This is undocumented, and every shipped case uses 0. Before building on it, a
one-step test establishes the actual sign. If confirmed, `wtsurf` is changed to
upward-positive (DALES convention) in its own commit, with docs updated.
The documentation example using `wtsurf = 0.01` is reviewed.

## Component 1: Deardorff moist buoyancy (`modthermodynamics.f90`)

uDALES `calthv` already matches DALES for interior levels. Changes:

- **Surface level `kb`.** Replace `dthvdz(kb) = 0` with DALES's
  saturated/unsaturated surface coefficients:
  - saturated (`ql0 > 0`): `a_surf = (1 - qt + rv/rd·qs·(1 + Lv/(Rv T))) / (1 + Lv² qs/(cp Rv T²))`,
    `b_surf = a_surf·Lv/(T cp) - 1`
  - unsaturated: `a_surf = 1 + (rv/rd - 1)·qt`, `b_surf = rv/rd - 1`
  - `dthvdz(kb) = a_surf·dthldz + b_surf·thl0(kb)·dqtdz`

  `dthldz`/`dqtdz` come from the Monin–Obukhov surface (Component 0). With
  any other bottom BC, `dthvdz(kb)` keeps its current value, so existing
  cases are unchanged.
- **`tmp0`.** New field in `modfields`, the absolute temperature
  `thl0·exnf + (Lv/cp)·ql0`, filled in `calthv` as in DALES. It is used by
  statistics and is needed by microphysics later.
- The half-level virtual potential temperature is unchanged (the inline
  uDALES formula is identical to DALES's `calc_virt_pot_temp`).

The scheme acts only through the one-equation SGS model (`loneeqn`), so
the moist test cases run with `loneeqn = .true.`.

## Component 2: parametrised radiation (`modraddata.f90`, `modradiation.f90`)

Namelist `&NAMRADIATION` (same names as DALES):

| name | default | meaning |
|---|---|---|
| `iradiation` | 0 | 0 none, 2 parametrised; anything else aborts |
| `itimerad` | 0 | update interval [s]; 0 = every step |
| `rad_longw` | `.true.` | longwave on |
| `rad_shortw` | `.false.` | must stay `.false.` (aborts otherwise) |
| `rad_smoke` | `.false.` | absorber is scalar `isvsmoke` instead of `ql0` |
| `rad_ls` | `.false.` | accepted for DALES compatibility; no-op with a warning, because uDALES already applies `thlpcar` in `forces` |
| `isvsmoke` | 1 | scalar index of the smoke absorber |
| `rka` | 130. | extinction coefficient κ [m² kg⁻¹] |
| `dlwtop` | 74. | F0 [W m⁻²] |
| `dlwbot` | 0. | F1 [W m⁻²] |
| `lrad_dycoms` | `.false.` | add the DYCOMS free-troposphere term |
| `rad_dycoms_D` | 3.75e-6 | divergence D [s⁻¹] in that term |
| `rad_dycoms_qtthres` | 8e-3 | qt threshold defining the local `zi` [kg kg⁻¹] |

The `lrad_dycoms` namelist switch replaces DALES's `iradiation = 10` +
case `moduser` (uDALES has no `moduser`).

Algorithm (`radpar`, per column, per DALES):

- `lwpt(k) = Σ_{k'≥k} ρ a dz`, `lwpb(k) = Σ_{k'<k} ρ a dz`, with `a` = `ql0` or smoke
- `lwd = F0 exp(-κ lwpt)`, `lwu = F1 exp(-κ lwpb)`
- `thlprad(k) = -[(lwd+lwu)(k+1) - (lwd+lwu)(k)] / (ρ cp exnf dzf)`
- with `lrad_dycoms`: find the highest level with `qt0 > qtthres`, interpolate
  `zi`, and for `zf > zi` add `D/3·[(z-zi)^{1/3} + zi (z-zi)^{-2/3}]`
  (DALES form, α_z = 1, ρ_i/ρ ≈ 1). **Sign:** this term is a **cooling**,
  `thlprad -= D/3·[...]`. DALES adds it with a plus sign. In Stevens et al.
  (2005), the flux `ρ_i cp D α_z [(z-zi)^{4/3}/4 + zi (z-zi)^{1/3}]` grows
  with height. Its tendency `-(D/3)·z·(z-zi)^{-2/3}` exactly cancels
  subsidence warming `-w ∂θl/∂z` of the RF01 free-troposphere profile
  `θl = 297.5 + (z-zi)^{1/3}`. The deviation from DALES is commented in the
  source, and a tier-1 test checks that the two cancel.

DALES `rad_user` bugs fixed in the port (noted in a comment): the `i` loop
runs to `j1` instead of `i1`, and `zi` is not reset per column (a column
without a crossing reuses the previous column's `zi`). Columns without a
crossing get no free-troposphere term.

Wiring:

- `initradiation` after `initthermodynamics`.
- `radiation` at the start of the RK step, before `advection`, as in DALES.
  It recomputes `thlprad` on `rk3step == 1` when due, then adds it to `thlp`
  on every substep.
- `exitradiation` in finalisation.
- `thlprad`, `lwd`, `lwu` are zeroed in solid cells; absorber inside solids is 0.
- Consistency checks: `iradiation = 2` requires `ltempeq`; `rad_smoke`
  requires `nsv ≥ isvsmoke`; without `rad_smoke`, it requires `lmoist`.
- Namelist read and broadcast in `modstartup`; entries in
  `tools/python/namelists.json`; docs page section.

## Component 3: diagnostics (`modthermostat.f90`)

A self-contained counterpart of DALES `tmser`/`profiles`, kept separate from
`modstatsdump` so it does not collide with the statsout rewrite and can later
move into it as one unit. Namelist `&NAMTHERMOSTAT`: `lthermostat`
(default `.false.`), `tthermostat` [s].

- `thermostat_t.<expnr>.nc` (time series): `lwp_bar`, `cfrac`, `zi`
  (maximum-`thl`-gradient definition as in DALES), `wmax`, `qlmax`; for smoke
  runs also the smoke-based `zi`.
- `thermostat_z.<expnr>.nc` (slab profiles): `thl`, `qt`, `ql`, `cfrac(z)`,
  resolved `wthl`/`wqt`, `lwu`, `lwd`, `thlprad`.
- Variable names and units match DALES netCDF output so comparison scripts
  use one set of keys. Averages are over fluid cells (`IIc`). Reductions use
  the existing uDALES `avexy_ibm`/MPI helpers.

## Testing

### Tier 1 — CI (`tests/integration/thermo/`)

Registered in `tests/test_suites.yml` as `class: supported`, `cost: fast`,
`component: solver`. Small runs (16×16×~64, ≤4 ranks, ≤1 min each) with
independent NumPy reference calculations:

1. **calthv:** one step on a prescribed profile with saturated and
   unsaturated layers; `dthvdz` (including `kb`) against a NumPy version.
   Output via a new optional `fielddump` variable `dthvdz` (off unless
   listed in `fieldvars`).
2. **Radiation, liquid:** one step from the RF01 initial state; `thlprad`
   against the analytic Stevens et al. (2005) formula evaluated in NumPy on
   the dumped `ql`, `rhof`, `exnf`.
3. **Radiation, smoke:** same for the smoke absorber.
4. **Decomposition independence:** test 2 serial vs x-split vs xy-split must
   agree to round-off (pattern of `tests/integration/sgs_statistics`).
5. **Off by default:** existing `david_tests` regression unchanged.
6. **#375:** a case with sentinel `thls`/`qts` and `ltempeq` aborts with the
   new message.
7. **Surface sign (D1):** one step with only a surface heat/moisture flux.
   Positive `wtsurf`/`wqsurf` must warm/moisten the first level by
   `flux·dt/dzf`.
8. **Monin–Obukhov:** `u*`, `obl`, `thls` from the first output, against a
   NumPy port of the same iteration on the logged slab-mean inputs; for
   `BCbotm = 5`, `u* = ustin` exactly.
9. **Free-troposphere balance:** a laminar RF01 column (no perturbations)
   with `lrad_dycoms`: above `zi + 50 m`, radiative plus subsidence `thl`
   tendency is ≈0 (relative to either term, ≤1e-3).

### Tier 2 — full-size validation on CX3 (`examples/301`–`303`, heavy)

| case | grid / domain | duration | reference |
|---|---|---|---|
| 301 BOMEX | DALES case inputs (64×64, 6.4 km, kmax 80) | 6 h | DALES `results_3.2.0`; Siebesma et al. 2003 ranges (cloud cover, LWP, mean profiles) |
| 302 Smoke | DALES case inputs (64×64×250, 3.2 km) | 3 h | DALES `results_3.2.0`; Bretherton et al. 1999 (entrainment rate, `zi` growth) |
| 303 DYCOMS RF01 | Stevens et al. 2005 specification, cross-checked against PyCLES `InitDYCOMS_RF01` (96×96, 35 m, dz 5 m, 3.36 km; `BCbotm = 4` with `z0` chosen for Cm ≈ 0.0011) | 4 h | Stevens et al. 2005 ensemble (LWP, `zi`, entrainment rate, cloud cover) |

- BOMEX and smoke inputs are converted from DALES by
  `examples/30x/convert_dales_case.py` (8-column DALES `lscale.inp` to the
  10-column uDALES layout with `pgx = pgy = 0`; `dthlrad` → `thlpcar`).
  Surface: fixed fluxes via `BCbotT/BCbotq = flux`; momentum via the wall
  function. This differs from DALES `isurf = 3/4`, where `u*` comes from
  Monin–Obukhov similarity; the difference is a known caveat in the report.
- RF01 inputs are generated from the paper's specification. Every number
  (profiles, fluxes, D, F0/F1/κ, f) is cross-checked against the paper
  during implementation and cited in the case README.
- `compare_dales.py` per case writes overlay plots (uDALES vs DALES vs
  published ranges) and a pass/fail table. Pass = inside the published
  ensemble range over the averaging window the paper uses.

## Risks

- Surface-layer differences (wall function vs MO) may shift BOMEX cloud
  statistics. Mitigation: compare `u*` and surface fluxes first.
- RF01 is sensitive to numerics: many LES in Stevens et al. 2005 lose the
  cloud. A low-LWP result needs diagnosing (SGS, advection scheme) before it
  counts as a port bug.
- `thermostat` reductions cost a little per output step; negligible at the
  planned `tthermostat`.
