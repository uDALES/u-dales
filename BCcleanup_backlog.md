# Boundary-condition cleanup: retiring the DALES flat-surface scheme + inflow-subsystem tidy

> Backlog + rationale for retiring the legacy flat-domain surface scheme in favour of the IBM
> facet path, and folding in the overlapping parts of the long-standing BC-subsystem cleanup —
> removing the recycling/rescaling inlet generator and reshaping `modinlet` into a general
> inflow-generation module (`modinflow`).
> Covers issues **#299**, **#302**, and **#68** (in part — see §5).
> Untracked scratch file — not part of the build. Source refs are `file:line` into `src/`,
> verified against the current tree.

---

## Summary

In uDALES every solid surface (ground, walls, roofs) is an **immersed-boundary facet**. The facet
path is fully self-contained: per-facet surface temperature `facT`, per-facet roughness
`facz0/facz0h`, energy-balance coupling — it **never reads `modsurfdata`** (the module's single
`use` site in `modibm.f90` is inside the legacy `subroutine bottom` itself). The old DALES scheme
(one domain-wide `thls`, `z0`, `wtsurf`, pinned to the flat plane `k=kb`) is redundant machinery
running alongside the real one.

Two distinct problems, one root variable (`thls`/`qts`) — keep them separate:

- **Path A (#299), latent:** `calc_halflev` plants the skin value into the air half-level field at
  `kb`, unguarded, every timestep. Verified **inert for all prognostic fields** in the current
  code (§1.3) — a landmine, not an active leak. Remove it for hygiene and to unblock the rest.
- **Path B (#302), live:** the thermodynamic base state (Exner, hydrostatic pressure, `thvs`) is
  anchored on `thls`, whose default is the **unvalidated sentinel `-1.`** — an unset `thls`
  poisons buoyancy at *every* level on *every* run with active thermodynamics. This is the real
  every-run hazard and the reason Phase 0 introduces an explicit reference.

The flat-surface scheme's guarded body is a no-op wherever `lbottom=.false.` (the default), but
`subroutine bottom` also does **unconditional** work every run (TKE ghost cells, fielddump
diagnostics — §3, Phase 1), so retiring it is a relocation + delete, not a bare delete. The two
cases that set `lbottom=.true.` (§4) are migrated to ground facets as part of the same change.

Phases:

- **Phase 0 — fix #299/#302.** Remove the `kb` overwrite; decouple the base state from "surface
  temperature" by **deriving it from the initial profiles + `ps`** — no new namelist input
  (decision resolved, see §1.5 and the Phase 0 checklist). Fold in three adjacent consistency
  fixes: `ps`-anchored Exner, SGS buoyancy re-point to `thvf(k)`, and a base-profile fallback for
  fully-solid slabs.
- **Phase 1 — retire `lbottom`.** Migrate the two cases to ground facets; relocate the
  unconditional code; delete the `bottom` / `wfuno` CASE(91)/(92) scheme.
- **Phase 2 — dissolve `modsurfdata`.** Remove the legacy scalars and dead members; move survivors
  to a base-state module.
- **Phase 3 — remove the recycling inlet generator; reshape into `modinflow` (#68, overlapping
  part).** The Lund machinery is not merely dormant — it is **unreachable** (§3 Phase 3), so there
  is no behaviour-preservation obligation. The wider #68 restructuring is out of scope — see §5.

All phases are backed by **new tests in the simulation suite** (`tests/`) that lock the
near-surface behaviour (§6).

---

## 1. Why this is needed

### 1.1 `thls`/`qts` are overloaded — three roles in one variable

1. **Flat-ground surface BC** — `thl0h(kb)=thls`, the `bottom` wall function. *Redundant; IBM owns
   all surfaces.*
2. **Base-state anchor** — `thvs = thls·(1+(rv/rd−1)·qts)`, hydrostatic `presf`, Exner `exnf/exnh`.
   *Real physics; must be re-sourced from an explicit reference.*
3. **Inlet-generator & stats reference** — turbulent-inflow scaling (`modinlet.f90`), buoyancy-flux
   diagnostic (`modstatsdump.f90`). *`modinlet` is removed outright (Phase 3); `modstatsdump` is
   re-pointed at the explicit reference.*

Exhaustive consumer table (every live `thls`/`qts` read in `src/`), so "nothing left" is checkable:

| Consumer | Where | Disposition |
|---|---|---|
| `thl0h(kb)=thls`, `qt0h(kb)=qts` | modthermodynamics.f90:526,536 | delete (Phase 0) |
| Exner first guess `exnf/exnh` | modthermodynamics.f90:291-292 | re-point to `thlref` (Phase 0) |
| `thvs = thls·(…qts)` | modstartup.f90:523 | re-point to `thlref`/`qref` (Phase 0) |
| `subroutine bottom` → `wfuno` | modibm.f90:2005,2023,2046 | deleted with `bottom` (Phase 1) |
| Buoyancy-flux diagnostic | modstatsdump.f90:2126 | re-point to `thlref` (Phase 2) |
| Inlet generator | modinlet.f90 (throughout); inlet branches modstartup.f90:1383,1826,2348 | deleted (Phase 3) |
| — comments only, no action | modsubgrid.f90:388,488; modforces.f90 (commented) | none |

`thvs` (derived from `thls` at modstartup.f90:523) has its own live readers, not in the table:
modsubgrid.f90:388,488 (SGS length scale / TKE buoyancy — re-pointed to `thvf(k)` in Phase 0) and
modthermodynamics.f90:395 (`thvh(kb)` hydrostatic anchor — re-pointed to `thv_b(kb)`). The
resolved buoyancy stopped reading it in the `ILS13` edit (modforces.f90:75,83; the `thvsi` fossil
survives only in comments at lines 75-76 — its unused declaration was removed by the #323 warning
cleanup) — see Appendix B.

### 1.2 The physical subtlety: *which* temperature, and *where*?

**(a) Surface (skin) temperature ≠ air temperature at the surface.** In DALES, `thls` is the
**skin** temperature; it drives a surface **flux** via Monin–Obukhov using the **air** temperature
at the first level and a transfer coefficient — there is deliberately a jump `dT = T_air − T_surf`
across the surface layer. But `calc_halflev` does `thl0h(i,j,kb) = thls`
(modthermodynamics.f90:526): it plants the **skin** value straight into the **air** half-level
field, as if air at the bottom face equals skin temperature. The facet path never conflates them:
`heat_transfer_coef_flux(... Tair, facT(fac,1) ...)` (modibm.f90:1538) keeps air and skin as
separate inputs bridged by a transfer coefficient.

**(b) Where is `kb`, physically?** `thl0h(kb)=thls` assumes `kb` is "just above the ground." With
the IBM that is false: over a **building or elevated terrain**, `kb` is deep inside solid and the
lowest **fluid** cell sits at a column-dependent `k_air(i,j) ≫ kb`. A single scalar at a single
flat level cannot represent a rough, built-up surface; any consistent "surface value" is
intrinsically per-column / per-facet — exactly the geometry the facet `bound_info` + wall distance
`bnddst` already encode (modibm.f90:70-105, 1345-1362). The flat `kb` override is a relic of a
world with one ground plane at one height.

**Implication for #299.** "A consistent approach to link `thl0h` to the ground surface temperature"
cannot be a flat `thl0h(kb)=thls`. Consistency means letting the IBM own the bottom interface and
giving `calc_halflev` a purely internal, geometry-agnostic value at `kb`.

### 1.3 Path A (#299) is latent — verified consumption chain

The `kb` overwrite (modthermodynamics.f90:526,536) is unguarded and runs every timestep
(`thermodynamics` called unconditionally from program.f90:215; `calc_halflev` unconditionally at
modthermodynamics.f90:70). But tracing consumption shows it currently reaches no prognostic field:

- `thl0h`/`qt0h` are read **only** by `calthv` (→ `thv0h`, modthermodynamics.f90:147) and `thermo`
  (→ `ql0h`). No advection, diffusion, or statistics route reads them.
- `thv0h` has exactly **one** consumer: the buoyancy term modforces.f90:83, whose loop runs
  `k = kb+1, ke` — `thv0h(kb)` is never read. `w(kb)` is additionally pinned to zero
  (modboundary.f90:165-166; `mask_w(:,:,kb)=0` modibm.f90:156).
- The hydrostatic pressure (`fromztop`, modthermodynamics.f90:364-406) is built from slab averages
  `th0av/qt0av` and anchored on `thvs` (line 395) — it never touches `thl0h(kb)`.

So the harm today is indirect: with `thls` at its `-1.` sentinel the overwrite drives the
saturation iteration in `thermo` unphysical, which is why the `tl<100 → tl=100` clamp exists
(modthermodynamics.f90:486-491, "X. Long: fix to tackle incorrect thls input"). The clamp protects
a value nobody consumes; fixing the root cause removes it. The removal must nonetheless keep
`thl0h(kb)`/`qt0h(kb)` sane (dumps, FP-trap builds) — see the ghost-cell note in Phase 0.

### 1.4 Path B (#302) is live — base state rides on an unvalidated sentinel

`thls` defaults to `-1.` (modsurfdata.f90:61; likewise `qts`, `thvs`) and **no namelist validation
exists**: it is read, broadcast (modstartup.f90:525) and used as-is. The chain
`thls → thvs` (modstartup.f90:523) → Exner first guess (modthermodynamics.f90:291-292) →
hydrostatic `presf(kb)` anchor (modthermodynamics.f90:395-398) → `exnh` → `thv0h` at **all**
levels → buoyancy (modforces.f90:83) means an unset `thls` corrupts buoyancy everywhere,
independent of Path A. Phase 0 replaces this with a base state **derived** from the initial
profiles and a validated `ps` (§1.5).

### 1.5 Sensitivity analysis — the reference is a linearisation constant (resolves Phase 0's open decision)

Where the `thls`-derived reference actually enters the dynamics today (post-`ILS13`, Appendix B):

| Consumer | Formula | Sensitivity to a 10 K error |
|---|---|---|
| Exner first guess | `exnf = 1−g·zf/(cp·thls)` (modthermodynamics.f90:291-292); refined from `ps`-anchored pressure before any use (312, 329-334) | ≈ zero (enters only via `ql0av/exnf` in the first-guess `th0av`, then twice refined; *exactly* zero for dry runs) |
| Hydrostatic anchor | `thvh(kb)=thvs` → `presf(kb)` over half a cell (modthermodynamics.f90:395-398) | parts in 10⁵ of pressure |
| SGS length scale / TKE buoyancy | `grav/thvs` (modsubgrid.f90:388,488) | ~3% of a modelled term |
| Resolved buoyancy | none — anomaly vs evolving `thvh(k)` (modforces.f90:83); loop `kb+1..ke` never reads `thvh(kb)` | zero |

Any value within tens of K of the near-surface air temperature is dynamically indistinguishable;
only the unset `-1.` sentinel is catastrophic: it reaches `thermo` through the Path A overwrite
(`thl0h(kb)=−1` → the X. Long clamp), corrupts the `presf(kb)` anchor by ~4–5% via `thvs<0`, and
flips/NaNs the SGS buoyancy via `grav/thvs` (note `sqrt(grav/thvs·|dthvdz|)` in `zlt`). Worse, the
sentinel is **partially live in the current test suite**: dry cases that set `thls` but leave
`qts` unset (e.g. tests/cases/101) run with `thvs = thls·(1+0.61·(−1)) ≈ 0.39·thls ≈ 116 K`, so
their SGS buoyancy terms are ~2.5× too strong — a real bug the derived `thv_b` fixes. An input whose *value* carries no information but whose
*absence* is fatal should not be an input. A startup-time constant is also sufficient: all
first-order time dependence rides on the evolving slab means (`th0av/qt0av/ql0av`, `thvh/thvf`);
the reference only linearises. Even a 15 K canopy-air drift over an EB run leaves the
linearisation error at the 0.1 K / few-% level.

**Decision (2026-07-15): derive the base state, don't ask for it.** At cold start, from
`prof.inp` + validated `ps`: `thl_b(k)=thlprof(k)`, `qt_b(k)=qtprof(k)` → `thv_b(k)`; hydrostatic
`p_b(k)` integrated upward from `p_b(z=0)=ps`; `exn_b(k)=(p_b/pref0)**(rd/cp)`. Fixed for the
run and logged at startup. No restart-format change is needed: `prof.inp` is read on warm starts
too (modstartup.f90:1103 cold, 1624 warm, 1006 stratstart), so the base state is recomputed
*identically* at every start — never from the evolving fields, which would break restart
continuity. No `thlref` namelist key exists to forget or to set inconsistently with the profiles.
`ps` stays a required input — it genuinely sets absolute pressure/saturation and cannot be
derived.

**IBM base height, and the profiles below the ground surface.** The base state is a property of
the full domain column `kb..ke+kh`, not of the fluid: `prof.inp` defines values at every level
(including levels inside terrain), and the hydrostatic integration runs from z = 0 regardless of
where the lowest fluid cell sits. Below the ground surface, `p_b/exn_b/thv_b` are the **analytic
continuation of the reference column** — the terrain replaced by reference air. Under variable
terrain height this is the only consistent 1-D choice (there is no single ground level to
re-anchor on), and it is harmless: solid-level values are never consumed by the dynamics (IBM
masks); they matter only through (a) the integration constant carried to the first fluid level,
where replacing terrain by reference air errs at second order, ~(g·Δz_terrain/cp)·(Δθ/θ), and
(b) dumps/statistics, where a smooth reference beats a sentinel. `ps` keeps a sharp meaning:
pressure of the reference column at the domain bottom z = 0, not at the (ill-defined) terrain
top.

**The current code gets the below-ground march wrong wherever a slab is fully solid.**
`avexy_ibm` returns `-999.` for slabs with no fluid cells, papering over `kb` only with a hack
its own comment disowns ("value at kb is used in modthermo … potentially should account for in
modthermo", modmpi.f90:644-659). `fromztop` then marches `presf` upward through the
`-999`-poisoned `thvh` at buried levels: the march mis-weights the solid layers (a
garbage-dependent bias in `presf` at *all* levels above, of order the buried layers' air weight,
~0.12 hPa per metre — a real if small qsat/T bias for moist runs) and the buried levels dump
garbage. Phase 0 fixes this with the base-profile fallback (`IIcs(k)==0 → thl_b/qt_b`, `ql=0`)
and deletes the `kb` hack.

**Correction (2026-07-16): there is no Exner anchoring inconsistency.** An earlier note here
claimed `thermo` mixes a `pref0`-anchored π with `ps`-anchored pressure. That was wrong:
`diagfld` refines the linear first guess through two `fromztop` passes and **recomputes**
`exnf/exnh` from the `ps`-anchored pressure (modthermodynamics.f90:312,318,329-334 — DALES-4.0
machinery) before `thermo` consumes them. The remaining Phase 0 action is hygiene only: seed the
iteration with `exn_b` instead of the `thls` formula, which removes the last `thls` read in
`modthermodynamics` and is bitwise-neutral for dry runs (the first guess enters `th0av` only
through the `ql0av/exnf` term).

---

## 2. Inventory

### (A) Legacy flat-ground BC — remove (with relocations)

| What | Where | Note |
|---|---|---|
| `thl0h(kb)=thls`, `qt0h(kb)=qts` | modthermodynamics.f90:526,536 | unguarded; latent (§1.3) — Phase 0 |
| `subroutine bottom` guarded body | modibm.f90:2020-2091 | guarded by `lbottom`; consumes `thls,z0,z0h,wtsurf,wqsurf` (Phase 1) |
| `subroutine bottom` **unconditional** code | modibm.f90:2010-2011 (e120/e12m `kb-1` ghost), 2014-2017 + 2093-2096 (`tau_*`/`thl_flux` bracket) | ⚠ runs on every run — **relocate/retire, not delete** (Phase 1) |
| `call bottom` + import | program.f90:153, program.f90:38 | remove with the subroutine (Phase 1) |
| `wfuno`/`wfmneutral` `CASE(91)/(92)` | modwallfunctions.f90:81-164, 307-350 | in-code TODO "should be moved out"; only caller is `bottom` (Phase 1) |
| `BCbotT/BCbotq/BCbots/BCbotm` | modglobal.f90:160-174; branches in `bottom`; namelist read/broadcast modstartup.f90:136,422-425; `BCbotm` write modstartup.f90:816 | value only ever *read* inside `bottom` (Phase 1) |
| `z0,z0h` (the `modsurfdata` scalars) | modsurfdata.f90:72-73 | ⚠ per-facet `facz0/z0` locals in modibm/modwallfunctions are DIFFERENT — do not touch |
| `wtsurf,wqsurf` | modsurfdata.f90:80,84 | flat-plane prescribed fluxes (Phase 1) |
| `wsvsurf` | modsurfdata.f90:86 | allocated+broadcast but **never consumed** (kb scalar flux hard-coded `+0.`, modibm.f90:2080) — dead |
| `tau_x/tau_y/tau_z/thl_flux` fielddump fields | modfields.f90:82,476-479; modfielddump.f90:230-240 | only ever record `bottom`'s contribution — retire or re-point at IBM wall functions (Phase 1 decision) |

### (B) Keep — genuine top-of-domain BC or base state (not "the ground")

- Top BCs via `fluxtop`/`valuetop` at `ke+1` in `modboundary.f90`:
  `thl_top, qt_top, wttop, wqtop, sv_top, wsvtop`. Real physical top boundary, unrelated to IBM.
- `ps` (reference pressure): keep as a required, **validated** input — it sets absolute
  pressure/saturation and cannot be derived. `thvs` (base-state virtual temp): keep the
  *quantity*, but as the derived base-state profile `thv_b(k)` (§1.5), no longer an input-coupled
  scalar.

### (C) Dead already — safe to delete on sight

`tskin, qskin, thlflux, qtflux, svflux, obl, oblav, Cm, Cs, Cmav, Csav, horvel,
ustar (the modsurfdata one), svs, dqtdz, lmostlocal`, and `dudz/dvdz/dthldz` (comments only).
Zero live references (no `use modsurfdata, only:` list anywhere imports them). Deleting these
shrinks `modsurfdata` to almost nothing — matching the module's own header wish: *"whole module
should be removed and variables moved"* (modsurfdata.f90:30).

---

## 3. Backlog

Single work stream. Sequence Phase 0 → 1 → 3a → 2 (Phase 2's final deletion needs Phase 3a's
`thls` removal); 3b can trail. Tests in §6 land with the phase whose behaviour they lock in.

### Phase 0 — fix #299/#302 (thermodynamics + base-state decoupling)

- [x] Remove `thl0h(ib:ie,jb:je,kb)=thls` and `qt0h(ib:ie,jb:je,kb)=qts`
      (modthermodynamics.f90:526,536).
- [x] Decide the replacement value at `kb`. Bare deletion leaves `calc_halflev` interpolating
      against the ghost `thl0(kb-1)`, which is set **once at startup** (modstartup.f90:1209) and
      never refreshed — that staleness is why the override exists. Options: one-sided
      `thl0h(kb)=thl0(kb)`, or refresh the scalar ghost at `kb-1` in `modboundary` alongside the
      existing `ekm/ekh` ghosts (modboundary.f90:452-453). Either is fine for prognostics (§1.3);
      pick the one that keeps dumps/FP-trap builds clean.
- [x] **Derive the base state** (decision: §1.5 — no new namelist input): at every start (cold,
      warm, stratstart — `prof.inp` is read on all three, modstartup.f90:1103/1624/1006) build
      `thl_b/qt_b/thv_b/p_b/exn_b(kb:ke+kh)` from the profiles + `ps`; no restart-format change;
      log at startup. Re-point the consumers: `exnf/exnh → exn_b`
      (modthermodynamics.f90:291-292), `thvh(kb) → thv_b(kb)` (modthermodynamics.f90:395); delete
      `thvs = thls·(…)` (modstartup.f90:523).
- [x] **Re-seed the Exner first guess with `exn_b`** — hygiene, not a bug fix (see the §1.5
      correction: the final Exner is already `ps`-anchored by the two-pass recompute,
      modthermodynamics.f90:312,329-334); removes the last `thls` read in `modthermodynamics`;
      bitwise-neutral when dry.
- [x] **Re-point the SGS buoyancy terms** `grav/thvs → grav/thvf(k)`
      (modsubgrid.f90:388,488), completing the `ILS13` modernisation (Appendix B) and matching
      DALES ≥ 4.0.
- [x] **Base-profile fallback for solid slabs**: in the thermo-facing slab averages, where
      `IIcs(k)==0` use `thl_b/qt_b` (and `ql=0`) instead of `avexy_ibm`'s `-999.`; delete the
      `kb` hack (modmpi.f90:644-650). Fixes the hydrostatic march below the ground surface
      (§1.5).
- [x] **Validate `ps`** at startup (required whenever thermodynamics is active; cannot be
      derived). `thls`/`qts` stay in the namelist for now — they are still consumed by the
      `lbottom` wall function (modibm.f90:2023,2046) — but are removed from the base-state path,
      and `lbottom=.true.` gains validation that `thls` is set. Full namelist removal (with the
      loud unknown-key migration signal) moves to **Phase 1**, when `lbottom` is retired.
- [x] Remove the `tl<100` clamp (modthermodynamics.f90:486-491) once the root cause is gone and
      the new tests pass.
- [ ] **Staging for regression:** land as two commits — (1) the `kb`-overwrite deletion, compared
      **bitwise**; (2) the base-state re-source, which is *not* bitwise even with matched values
      (discrete hydrostatic `exn_b` differs from the linear formula at O(10⁻⁵)) — compare under a
      tight numerical tolerance bounded by the §1.5 sensitivities.
- [ ] Tests: §6.1, §6.3, §6.6.

### Phase 1 — retire `lbottom` and the flat-wall scheme

- [x] **Retire the `lbottom` namelist switch** (modibm.f90:49; modstartup.f90:153,445) — done via
      wholesale deletion of the flat-surface scheme (Task 2) rather than migrate-then-retire.
- [x] **Migrate the cases that use `lbottom`** (§4) to ground facets. `examples/999/namoptions.999`
      and `tests/regression/david_tests/cases/103/namoptions.103` set `lbottom=.true.` and were
      deliberately left broken-parse (unknown namelist key) — Tasks 4-5 rewrite them wholesale as
      facet migrations. **Case 103 done (Task 4, 2026-07-16):** flat 8x8 ground generated via
      `udgeom.create_flat_surface` + IBM f2py preprocessing with `stl_ground=True`; 128 floor
      facets cover the full footprint, `solid_c/u/v` are empty everywhere (no solid volume
      introduced — `kb` stays fluid), `factypes.inp.103` carries `z0=0.01`/`z0h=0.000067`,
      `Tfacinit.inp.103` is 288 K on all facets, `namoptions.103` sets `libm=.true.`,
      `iwallmom=2`/`iwalltemp=2`. **Case 999 done (Task 5, 2026-07-16):** flat 64x64 ground
      generated as a *minimal 2-triangle* STL (`udgeom.create_flat_surface(64, 64, edgelength=64)`)
      rather than per-cell facets — this case's 128x128 floor would otherwise mean 32768 facets;
      the IBM f2py backend maps all 128x128 boundary cells to the 2 facets regardless (confirmed:
      `nbndpts_c`/`nsolpts_w` = 16384, full floor coverage, `solid_c/u/v` empty). This case is
      neutral (`ltempeq`/`lmoist` both off, unset in the file), so `namoptions.999` sets
      `iwallmom=3` (the neutral wall function the solver would force anyway per
      `modstartup.f90:803-805` when `ltempeq=.false.`) rather than `iwallmom=2`; no `Tfacinit`
      file is needed since `initfac.f90:297`'s read guard
      (`lEB .or. iwalltemp==2 .or. iwallmom==2 .or. iwallmoist==2`) is false throughout.
      `factypes.inp.999`'s existing wallid=1 row already carried the old `&BC` roughness
      (`z0=0.05`, `z0h=0.00035` — the generic-concrete default), so no factypes edit was needed.
      Both migrations fixed `&DOMAIN xsize/ysize` where present (not valid Fortran namelist
      keys — replaced with `xlen`/`ylen`; `examples/999` already used `xlen`/`ylen`, no fix
      needed there) since the files wouldn't have parsed otherwise; grid *values* are unchanged
      throughout. Still open: `examples/024/namoptions.024:48`
      (which set the now-deleted `BCbotT = 2`, dead even pre-refactor since it defaults
      `lbottom=.false.`) has been cleaned as part of Task 3's sweep, along with every other
      namoptions file under `tests/` and `examples/` that set a pruned key (except the two above).
- [x] **Relocate the unconditional code first** (it runs on every run regardless of `lbottom`):
  - [x] `e120/e12m` ghost at `kb-1` (modibm.f90:2010-2011) → `modboundary`, next to the `ekm/ekh`
        ghosts. Nothing else sets it; the subgrid model reads it at `kb`.
  - [x] Decide fate of `tau_x/tau_y/tau_z/thl_flux`. **Corrected picture (2026-07-16):** the
        original claim that they only captured `bottom`'s contribution was wrong — `ibmwallfun`
        (guarded by `libm`, not `lbottom`) accumulated real IBM wall-function tendencies into
        them, with `bottom`'s overwrite doubling as the per-step reset; they were a live but
        undocumented, in-repo-unexercised per-cell diagnostic whose reset semantics were
        entangled with the scheme being deleted. **Decision: retire** (fields, allocations,
        `ibmwallfun` writes, fielddump 'tx'/'ty'/'tz'/'hf' cases). The per-facet `fac_tau_*`
        statistics remain; reinstating a per-cell version cleanly needs its own reset in
        `ibmwallfun` (one-commit revert + small fix if ever wanted).
- [x] Delete `subroutine bottom` (modibm.f90:1997-2099), its call site (program.f90:153) and
      import (program.f90:38).
- [x] Delete `CASE(91)/(92)` surface blocks in `wfuno`/`wfmneutral`
      (modwallfunctions.f90:81-164, 307-350). **Went further:** post-deletion both
      subroutines had zero remaining cases and zero call sites anywhere in src/ (the
      IBM path uses `wallfunmom`/`wallfunheat` in modibm — the scouted facet-code
      callers never existed), so `wfuno`, `wfmneutral`, and their orphaned Uno1995
      helpers `unoh`/`unom` were deleted wholesale. `modwallfunctions` is now an
      empty module shell (tombstone comment points here); the file can be dropped
      entirely in a follow-up — the build globs `src/*.f90`, so no CMake edit needed.
- [x] Remove `BCbot*` end-to-end: declarations (modglobal.f90:160-174), namelist read/broadcast
      (modstartup.f90:136,422-425), and the `BCbotm` write (modstartup.f90:816).
- [x] Update the pre-processing surfaces: `tools/preprocessing.m:224` (`addvar 'lbottom'`) and
      `tools/python/namelists.json:282,393`.
- [ ] Tests: §6.2 (plus §6.5 bitwise check for `lbottom=.false.` runs).

### Phase 2 — dissolve `modsurfdata`

- [x] Remove the (A) `modsurfdata` scalar **declarations**: `thls, qts, z0, z0h, wtsurf, wqsurf,
      wsvsurf`. Task 2 already dropped their namelist membership, reads, and broadcasts
      (modstartup.f90 &BC + `use modsurfdata` lists) — the declarations were kept deliberately
      because the unreachable `iinletgen` branches in `readinitfiles`/`readrestartfiles`
      (modstartup.f90 ~1380,1827,2349) still dereference `thls`. Precondition: the consumer
      table in §1.1 is empty apart from the re-point targets — i.e. Phases 0, 1 and 3a have
      landed.
- [x] Delete the (C) dead `modsurfdata` members.
- [x] Re-point `modstatsdump`'s buoyancy-flux diagnostic (`tkestatsdump`, modstatsdump.f90:2126)
      from `grav/thls` to `grav/thv_b(kb)` — pulled forward into Task 2 (Step 1) to unblock the
      `thls` namelist prune below; `modbasestate` already existed pre-Phase-2.
- [x] Move survivors `ps` and the derived base-state profiles (`thl_b/qt_b/thv_b/p_b/exn_b`,
      §1.5) into a new `modbasestate`; then remove `modsurfdata` entirely.
- [x] Update docs: `docs/udales-namoptions-overview.md` (rows 139-152, 200) and
      `docs/udales-example-simulations.md` (lbottom section, line 198). Done in Task 3, ahead
      of the rest of Phase 2 landing.

### Phase 3 — remove the recycling/rescaling inlet generator (#68, overlapping part)

Scope definition: **"the inlet generator"** = the Lund-1998 recycling/rescaling machinery
(`modinlet.f90` in full, plus every parameter that exists *solely* to serve it), gated by
`iinletgen` (default 0, modglobal.f90:179).

**It is unreachable, not just dormant.** `initinlet` and `inletgen` are never called — the only
call sites are commented out (program.f90:85; modstartup.f90:629). `iinletgen=1` is a silent
no-op; `iinletgen=2` is half-wired (`readinletfile` is called at modstartup.f90:1447,1877, but
`initinlet` never runs and `modboundary` never applies inlet planes). No case in `examples/` or
`tests/` sets `iinletgen ≠ 0`. Consequence: **no behaviour-preservation obligation for any
`iinletgen` mode** — this is dead-code removal, and §6.4 reduces to a build + suite check.

The trap is a set of **false-friend** parameters that carry "inlet" in their name or merely *live*
in `modinletdata` but are consumed by kept subsystems — relocate, don't delete (§3c).

#### 3a. Remove (pure recycling/rescaling — nothing else uses it)
- [x] `modinlet.f90` in full (all of `inletgen`, `blthickness`, `dispthickness`,
      `momentumthickness`, `enthalpythickness`, `writeinletfile`/`readinletfile`, `wallawinlet`,
      z-interpolators), plus the remaining live call sites: `readinletfile`
      (modstartup.f90:1447,1877), `exitinlet` (modstartup.f90:2368), `zinterpolate*` uses
      (modstartup.f90:2176 ff.), and the commented `initinlet` stubs. **Done (Task 1,
      2026-07-16):** `git rm src/modinlet.f90`; all call sites/imports removed from
      modstartup.f90 (build globs `src/*.f90`, no CMake edit needed).
- [x] The recycling half of `modinletdata.f90` (lines 27-129, 132-190 except the shared vars in
      §3c): recycle/inlet-station scaling arrays (`Uinl/Urec`, `zir*/zii*/zor*/zoi*`, `heavi*`,
      `loc*`), thicknesses (`di/dr/dti/dtr/theta*`), friction scales (`utaui/utaur/ttaui/ttaur`),
      `irecy`, and the y-interpolation block for reading inlet files. **Done:** module trimmed to
      `ubulk, vbulk, iangle, iangledeg` + the driver block (`storeu0driver…chunkread_e`);
      every surviving symbol has a live consumer in modboundary/moddriver/modtstep (verified by
      grepping every `use modinletdata` site in `src/`). `irecy` itself was initially kept (its
      one write, `irecy = ib + iplane` in `modboundary.f90`, looked live) but turned out to be
      write-only — no read anywhere in `src/`; deleted in Task 3 along with the write and the
      `use modinletdata, only : irecy` import (`irecydriver`, a different variable consumed by
      the driver path, is unaffected).
- [x] Switches owned only by the generator: `iinletgen`, `lfixinlet`, `lfixutauin`, `linletRA`
      (modglobal.f90:179,195-197), and namelist params `di`, `dti` (read modstartup.f90:142,
      broadcast 544-545); `totinletav` reads used only by inlet-gen. **Done**, plus two
      scout-drift extras discovered during the cut (not in the original itemised list, but
      unreachable once the above landed): `lreadminl`'s modglobal.f90 declaration (its namelist
      entry/broadcast/`use` sites were already slated for removal, but the declaration itself was
      missed), and `nstore` (modglobal.f90) + the `nstepread==nstore+1` disjunct in
      modsave.f90's restart-trigger condition — `nstore` was "number of rk steps in inletfile",
      referenced live only there; the disjunct was always `.false.` (nstepread stuck at its old
      default 1, nstore defaulting to 1002, never namelist-configurable), so dropping it is
      bitwise-neutral. `lstoreplane`/`lwallfunc` were considered and initially **kept** in
      `&INLET` — they were equally unreachable post-deletion (only consumer was `modinlet.f90`)
      but outside that task's itemised scope; flagged as a follow-up. **Done (Task 3):** both
      removed end-to-end — modglobal.f90 declarations, `&INLET` namelist membership, and their
      `MPI_BCAST` calls in modstartup.f90; zero references left in `src/` bar a dead comment
      (`moddriver.f90:1194`, inside an already-commented-out block) and a stale-but-harmless
      mention in another broadcast's inline comment (fixed). `&INLET` now carries only
      `Uinf, Vinf, inletav` — matching this plan's original intent.
- [x] The inlet-gen `thls` branches in modstartup.f90 (1383, 1826, 2348). **Done** — removed along
      with the arms that contained them; `use modsurfdata, only:thls` dropped from both
      `readinitfiles` and `readrestartfiles` (no longer referenced in either).
- [ ] Tests: §6.4.

#### 3b. Repurpose the module as *the inflow-generation module* (data + functionality together)

**Carry-overs from the Phase-2 final review (2026-07-16, cosmetic — deferred here to keep the
frozen gate spans clean):** orphaned `INFO` group in tools/python/namelists.json (:100 + reverse
maps — described modinlet's inlet-metadata namelist); write-only locals
`uaverager/uaveragei/taverager/taveragei/waverage` in modstartup.f90 (~:950, ~:1644); stale
only-imports `Uinf, jgb, jge` in `readinitfiles`' modglobal list (verify `jgb/jge` strat-start use
first); duplicate `timee` in modsave.f90's only-list (pre-existing); comment fossils referencing
modinlet/lstoreplane in moddriver.f90:7,1194.

Rather than delete the file, turn it into the single, self-contained home for inflow generation:
the precursor/driver method moves in now, and a planned **synthetic turbulence generator** lands
here later. This applies the "data + functionality together" principle (the same reason
`modinletdata`/`modsurfdata` — data-only modules — are retired).

**Naming (proposed): `modinflow`.** `inlet` is narrow and legacy once the recycling method is
gone; `modinflow` is method-agnostic (driver + synthetic + whatever next) and pairs with an
`iinflow` selector. Keep the `mod` prefix — dropping it is a repo-wide convention change, out of
scope here. (Fallback: keep `modinlet` to avoid renaming call sites.)

**Naming, actual (Task 1, 2026-07-16):** `src/inflow.f90` / module `inflow` — no `mod` prefix,
per explicit user direction overriding the proposal above. `modinlet.f90` was already deleted in
Phase 2 (Task 1), so there was no surviving file to repurpose in place; `moddriver.f90` and
`modinletdata.f90` were merged into this new file instead, with `moddriver`'s subroutine bodies
carried over verbatim (only `initdriver`→`initinflow`, `exitdriver`→`exitinflow`, dropped
now-redundant internal `use modinletdata`/`use moddriver` only-lists, and a dropped dead `use
modsave, only: writerestartfiles` in `drivergen`).

- [ ] Strip the Lund recycling/rescaling code (§3a), leaving the module as the inflow-generation
      shell. **Superseded:** §3a already removed `modinlet.f90`/the recycling half of
      `modinletdata.f90` in Phase 2 — nothing left to strip here.
- [x] Fold `moddriver`'s functionality into it, and move the driver data
      (`storeu0driver…storesv0driver`, modinletdata.f90:194-200; `iangle`/`iangledeg`, consumed at
      moddriver.f90:465-466) in with it. Retire `moddriver` as a separate module. The three driver
      cases (examples/949, examples/950, tests/cases/525 — `idriver` 1/2) are the regression
      anchors for this move. **Done (Task 1, 2026-07-16):** all 7 subroutines and surviving data
      moved into `src/inflow.f90`; five never-read members (`tdriver`, `storee120driver`,
      `e120driver`, `storee12mdriver`, `e12mdriver` — zero live references beyond declaration,
      verified by repo grep) dropped along with the move; `storeumdriver` kept (allocated, though
      never subsequently read — pre-existing behaviour, not this task's to fix). Driver
      cases/regression run deferred per this task's gate (static verification only, no Fortran
      toolchain in this environment).
- [x] Relocate `ubulk`/`vbulk` (modinletdata.f90:130-131; assigned across modstartup.f90
      1342-1988, consumed only at modboundary.f90:159-160; *not* inflow-generation state) to the
      BC-owning module, **not** into the inflow module. **Done (Task 1, 2026-07-16):** now
      module data in `modboundary` (public, defaults `0.`, comments preserved).
- [x] Delete `modinletdata.f90` once emptied. **Done (Task 1, 2026-07-16):** `git rm`, alongside
      `moddriver.f90`.
- [ ] Consolidate the method selector: replace `iinletgen` (0/1/2) + `idriver` (0/1/2) with one
      switch `iinflow` (0=none, 1=driver, 2=synthetic). Deferred to Task 2 (selector rename).
- Sequencing: §3a is well-bounded and low-risk (dead code) and lands first; the
  `moddriver`→`modinflow` consolidation in §3b is a larger refactor — separate commit/PR.

#### 3c. Do NOT touch — genuinely shared, already outside `modinletdata`
- `Uinf`, `Vinf`, `ifixuinf` (declared in `modglobal`:333,334,252; free-stream forcing in
  `modforces`/`modtstep`; top BC value in `modboundary` `valuetop`). Named "used in inlet
  generator" but general.
- `inletav` (running-average window reused by the `ifixuinf==2` free-stream forcing,
  modforces.f90:190). Rename later if desired, but keep.
- `lzerogradtop`, `lzerogradtopscal` (top-BC switches).
- `uouttot`/`vouttot` (declared modfields.f90:389-390; general convective-outflow velocity for
  all iolet BCs) — see §5.

---

## 4. Cases in the repo that flag `lbottom` (must be migrated in Phase 1)

**Both cases migrated (Tasks 4-5, 2026-07-16)** — `lbottom = .true.` no longer appears anywhere
in the repo:

- ~~examples/999/namoptions.999:43~~ → `namoptions.999` now sets `libm=.true.` with a 2-facet
  ground STL (`geom.999.STL`) and generated `&WALLS` counts (Task 5).
- ~~tests/regression/david_tests/cases/103/namoptions.103:68~~ → `namoptions.103` now sets
  `libm=.true.` with a 128-facet ground STL (`geom.103.STL`) and generated `&WALLS` counts
  (Task 4).
- Documented in docs/udales-example-simulations.md:198 and docs/udales-namoptions-overview.md
  (which states `lbottom` is *"Used only if no ground facets"*, row 200) — these docs describe
  the retired scheme; not updated as part of this migration (tracked separately if still needed).
- Pre-processing surfaces: tools/preprocessing.m:224, tools/python/namelists.json:282,393 — already
  cleaned (Phase 1, `[x]` above).

These represented flat-floor setups with no ground facets. Migration = generate floor facets for
them (the conceptual replacement the docs already imply). Solver-run confirmation that results are
statistically equivalent to the pre-migration `lbottom` behaviour (§6.2 — the facet wall function
is not the CASE(91/92) one, so bitwise identity is not the criterion here) is **deferred** to a
Linux session with a full build toolchain; only static input-generation verification has been done
in this environment (see §6.2 below).

---

## 5. Scope of #68: what is in vs out

**§3a done (Task 1, 2026-07-16):** the recycling/rescaling inlet generator is removed —
`src/modinlet.f90` deleted; `modinletdata.f90`, `modstartup.f90`, `modglobal.f90`, `modsave.f90`,
`modtstep.f90` trimmed of every generator-only symbol (`iinletgen`, `lfixinlet`, `lfixutauin`,
`linletRA`, `lreadminl`, `di`, `dti`, `totinletav`, `nstore`, and the full recycling half of
`modinletdata`); `&INLET` namelist originally kept `Uinf, Vinf, inletav, lstoreplane, lwallfunc` —
**Task 3 (2026-07-16)** finished the sweep by removing `irecy` (write-only after Task 1) and
`lstoreplane`/`lwallfunc` (dead post-Task-1, deferred there as out-of-scope), so `&INLET` now
carries only `Uinf, Vinf, inletav`. Zero matches for the deleted names across `src/` bar deliberate
comment survivors (see `.superpowers/sdd/phase2-task-1-report.md` and
`.superpowers/sdd/phase2-task-3-report.md`). §3b (repurpose into `modinflow`, fold in `moddriver`)
is not started.

Issue #68 ("Restructure and retest modinlet and moddriver. Improved modboundary.") is a broad
BC-subsystem cleanup. Only the part that directly couples to this surface-BC work is pulled in
here.

**In scope (Phase 3):**
- Remove the recycling/rescaling inlet generator (`modinlet` + generator-only params, §3a). Done.
  Unreachable dead code (§3 Phase 3), recommended for removal in #68, and the last `thls` consumer
  outside the base state once Phases 0-1 land — removing it makes the `modsurfdata` dissolution
  clean.
- Retire the `modinletdata` data-only module by relocating its survivors (§3b).

**Out of scope (track separately, e.g. a follow-on "modboundary/driver consistency" issue):**
- `moddriver`'s internal restructuring beyond absorbing the relocated data.
- `uouttot`/`vouttot` *forcing-consistency* across the momentum-forcing modes (ties to #54) — the
  general convective-outflow velocity stays; only its `= ubulk` wiring is incidentally tidied when
  `ubulk` is relocated.
- `linoutflow` / `lper2inout` / `ltempinout` / `lmoistinout` switch revision.
- Consistent BC notation across directions (`BCxm`, `BCyq`, …) and symmetric x/y inflow-outflow.
- Lateral (y) inflow/outflow for scalars.
- Scalar top BC set to `ke+1` should be `ke+khc` (and `svprof` allocation); `sv_top = svprof(ke)`
  overwrite in `modstartup` — touches the same `sv_top`/`svprof` kept in Inventory (B), so
  coordinate if both land close together.
- ~~`scalSIRANE` rename/removal; `iosi`/`iohi`/`ioqi`~~ — already absent from `src/`; report on
  #68 as done rather than tracking here.
- `wqtop` is broadcast (modstartup.f90:512) but absent from the `namelist/BC/` list
  (modstartup.f90:131-138, which carries `wttop` but not `wqtop`) — pre-existing (predates this
  branch), found during Phase 2 review (Task 3). The moist top-flux BC has therefore never been
  settable from `namoptions`; it silently rides the `real :: wqtop = 0.` default
  (modboundary.f90:40) on every run. Not fixed here — out of scope for this cleanup pass; track as
  a follow-up (add `wqtop` to `namelist/BC/`, decide/validate a sensible default).

Rationale for the split: the out-of-scope items are about *lateral / inflow-outflow* boundary
conditions and momentum forcing — orthogonal to the *surface* (bottom) and *base-state* concerns
driving #299/#302. The line for Phase 3 is "grows-a-BL-by-recycling" machinery in, general
inflow/outflow/forcing infrastructure out.

---

## 6. Testing (simulation suite, `tests/`)

The suite is driven by `tests/test_suites.yml` (groups → suites; `kind` ∈ unit/integration/
reference/**regression**; `component` includes `ibm`) and `tests/run_tests.py`, with cases under
`tests/regression/`. Two distinct comparison regimes — say which one each test uses:

- **Bitwise** (identical prognostic fields over a short run): for pure dead-code removals where
  nothing consumed should change — the Phase 0 overwrite deletion with `thlref` set equal to the
  old `thls`, and Phase 1 on `lbottom=.false.` cases.
- **Statistical equivalence** (time/slab-averaged profiles over an averaging window, stated
  tolerance): for the `lbottom` → floor-facet migration, where the wall-function formulation
  legitimately changes.

New coverage to add:

- [ ] **§6.1 Near-surface correctness with elevated ground (Phase 0).** A small case with
      buildings / raised terrain and active thermodynamics (`ltempeq`) + moisture (`lmoist`),
      asserting that `thl0h`, `thv0h` and the near-surface buoyancy are physically sane and
      **independent of the removed `thls`** — the prognostic fields must be bitwise-unchanged by
      the overwrite deletion. This pins the skin-vs-air handling and turns the §1.3 inertness
      analysis into an executable guarantee.
- [ ] **§6.2 `lbottom` retirement regression (Phase 1).** Convert the migrated cases
      (examples/999, regression 103) to floor-facet equivalents; assert statistical equivalence
      to the pre-migration `lbottom` result (regime and tolerance stated in the suite entry).
      **Deferred (2026-07-16, Task 4):** case 103's inputs were generated and statically verified
      (facet/solid/fluid file parsing only — see `.superpowers/sdd/phase1-task-4-report.md`); no
      Fortran/MPI/NetCDF-Fortran toolchain is available in this environment to actually run the
      old-`lbottom` vs new-facet solver comparison. The statistical-equivalence run (old-vs-new
      103, slab-averaged profiles over an averaging window) is deferred to a Linux session with a
      full build — the facet wall-function formulation legitimately differs from the old
      `BCbotT=2` flat-surface one, so the acceptance criterion is statistical, not bitwise.
      **Deferred (2026-07-16, Task 5):** case 999's inputs were likewise generated and statically
      verified only (2-facet minimal ground STL; full 128x128 floor coverage confirmed on
      `solid_w`/`fluid_boundary_c` — see `.superpowers/sdd/phase1-task-5-report.md`); same
      toolchain gap. This case is neutral (`ltempeq=.false.`), so the old-vs-new comparison is
      against the old `z0=0.05`/`z0h=0.00035` flat-surface neutral wall function rather than a
      `BCbotT`-driven one — still a statistical (not bitwise) criterion, since the new facet
      wall function (`iwallmom=3`, `mom_transfer_coef_neutral`) is a different formulation from
      the old flat-surface one even in the neutral limit.
- [ ] **§6.3 Base-state invariance + validation (Phase 0).** Assert that Exner / hydrostatic
      pressure / buoyancy depend only on the derived base state (initial profiles + `ps`); assert
      startup **fails loudly** when `ps` is unset with active thermodynamics, and when
      `lbottom=.true.` without `thls` — a guard against the #302 coupling and the `-1`-sentinel
      failure mode re-appearing. (Rejection of `thls`/`qts` as unknown namelist keys lands with
      their Phase 1 removal.)
- [ ] **§6.4 modinlet removal (Phase 3a).** Build + full supported suite unchanged with
      `modinlet` gone. Cheap by construction: the generator is unreachable and no case sets
      `iinletgen ≠ 0`.
- [ ] **§6.5 `bottom` relocation check (Phase 1).** Bitwise regression on an `lbottom=.false.`
      case across the `e120/e12m`-ghost relocation, plus a decision-backed check on the
      `tau_*`/`thl_flux` fielddump fields (removed or re-pointed, not silently zeroed).
- [ ] **§6.6 Base state below the ground surface (Phase 0).** A terrain case whose lowest slabs
      above `kb` are fully solid (`IIcs(k)==0`): assert no `-999` reaches the thermo-facing slab
      means, that `presf` at the fluid levels equals the reference-column continuation (weight of
      the buried layers included), and that dumps at buried levels contain the base profiles
      rather than sentinels.

      **Progress (2026-07-16): fixture landed, run unexecuted.** `tests/cases/091/` (copy of
      `tests/cases/090/` with `geom.091.STL` — a programmatic, watertight, floor-covering box
      2 grid cells tall) was generated and regenerated end-to-end with the Python `udprep` IBM
      preprocessing (`ibm.run_all(backend="f2py")`, plus `seb`/`scalars` reruns for the
      facet-count- and geometry-dependent side files). Statically verified: `solid_c.txt`,
      `solid_u.txt`, `solid_v.txt` have exactly `itot*jtot=4096` unique `(i,j)` at both `k=1`
      and `k=2` and nothing at `k>=3` — i.e. slabs `kb`/`kb+1` are fully solid on the relevant
      grids. `tests/regression/bc_cleanup/run_test.py --assert-only <case>` was added (checks:
      run completes; startup log contains `Base state:`; no dumped field is `-999.`/NaN),
      wired to both `090` and `091`, and verified with `python -m py_compile` plus standalone
      unit checks of the new assertion helpers against synthetic fixtures. **Not done:** the
      solver was never actually run against case 091 in this environment (no Fortran/MPI/
      NetCDF-Fortran toolchain for building `u-dales` itself, only for the Python preprocessing
      f2py extensions) — the `presf`/reference-column-continuation and no-sentinel-in-dumps
      assertions above are implemented but unexercised. Whoever picks this up next should run
      `python tests/regression/bc_cleanup/run_test.py BCcleanup BCcleanup Debug --assert-only 091`
      on a machine with a full build toolchain before relying on this fixture.
- [ ] Register these under an `ibm` / `reference` group in `tests/test_suites.yml` so they run in
      the merge gate.

      **Correction (2026-07-16): Task 4's `5e-3` tolerance claim on case 090 was vacuous.**
      Case 090 only sets `lvreman = .true.` (the default), so it exercises solely the Vreman
      branch of `closure` in `src/modsubgrid.f90`. Task 4's re-point of `zlt`/`sbbuo` to the
      evolving `thvf(k)` (`src/modsubgrid.f90:387,486`) lives in the `loneeqn` branch, which
      090 never enters — so that commit is dead code from 090's point of view and 090 is
      expected **bitwise** (`--atol 0.0`) across the whole branch, including through Task 4.
      `tests/cases/092/` (090 + `lvreman = .false.` / `loneeqn = .true.`) was added so a case
      actually runs the re-pointed lines; the `5e-3` tolerance expectation from Task 4 applies
      to 092, not 090. See `tests/regression/bc_cleanup/run_test.py`'s per-case `default_atol`
      (090: `0.0`, 092: `5e-3`) and `docs/superpowers/plans/2026-07-16-bc-cleanup-phase0.md`
      (Task 4) for the matching correction.

---

## Appendix: how the IBM facet path applies surface BCs (for contrast)

- **Representation.** `nfcts` facets; per-facet arrays loaded in `initfac.f90:readfacetfiles`
  (facets.inp, factypes.inp, Tfacinit.inp, …). Surface temperature `facT(0:nfcts, nfaclyrs+1)`,
  column 1 = outdoor skin (initfac.f90:64,129); roughness `facz0/facz0h` (initfac.f90:215-216);
  `facqsat`, `fachurel`, resistances `facf`.
- **Fluid–facet geometry.** `bound_info_type` per staggered var (modibm.f90:70-105), loaded from
  `fluid_boundary_*.txt` / `facet_sections_*.txt`; each section ties a fluid boundary point, a
  facet id, an area and a **wall distance** `bnddst`.
- **Momentum.** `wallfunmom` (modibm.f90:1286): tangential velocity at reconstruction point,
  transfer coeff (`mom_transfer_coef_stability/neutral`), stress rotated from local
  (strm,span,norm) to global and projected on `dir` — any facet orientation, no floor/wall
  special-casing — applied to RHS.
- **Heat/moisture.** `wallfunheat` (modibm.f90:1436): `heat_transfer_coef_flux(Tair, facT(fac,1))`
  → flux into `thlp`; latent via `facqsat/fachurel`; accumulates `fachf/facef` for the EB.
- **Bottom of domain.** *Not* an IBM facet — special-cased in `subroutine bottom`
  (modibm.f90:1997, active only when `lbottom=.true.`, plus the unconditional code noted in §2A).
  `modboundary` at `kb` sets `w=0` and the no-slip `ekm/ekh` ghosts (modboundary.f90:165-166,
  452-453).
- **EB coupling loop.** `facT(fac,1)` → wall-function flux into `thlp/qtp` + `fachf/facef` →
  `modEB.f90` conduction solve → new `facT` → refreshed `facqsat`. Genuinely per-facet, closed
  loop.

---

## Appendix B: buoyancy-formulation provenance (the `ILS13` edit)

Recovered from DALES releases — uDALES's own history is squashed at the fork ("Initial
restructure"):

| Code | Resolved buoyancy | SGS reference |
|---|---|---|
| DALES 3.2 (fork-era base) | `wp += g/thvs·thv0h` — full thv, constant reference; the slab-mean part left to the pressure projection (exact only in periodic domains) | `g/thvs` |
| DALES ≥ 4.0 (incl. main) | `wp += g·(thv0h−thvh(k))/thvh(k)` — anomaly vs evolving slab mean | `g/thvf(k)` |
| uDALES today | DALES-4.0 form via the `ILS13` edit (modforces.f90:75,83; `thvsi` fossil comments at 75-76, declaration removed by #323) | still `g/thvs` (modsubgrid.f90:388,488) — Phase 0 re-points |

The explicit anomaly form is **required** in uDALES, not cosmetic: with inflow–outflow BCs and
IBM solid fractions, the k=0 pressure mode no longer exactly absorbs a horizontally uniform
forcing on `w`, so the old reliance on the projection breaks. `thvh(k)` is rebuilt every
timestep from IBM-masked slab averages (`avexy_ibm`, `IIcs` masks); the buoyancy loop runs
`kb+1..ke`, so `thvh(kb)` — the base-state anchor — is never read by it. uDALES extras vs DALES:
the `lbuoyancy` gate (off → θ fully passive, no buoyancy term at all), and no rain-water loading
(no microphysics).
