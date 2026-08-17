# BC Cleanup Phase 2 Implementation Plan — Phase 3a precondition + dissolve `modsurfdata`

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the provably-unreachable recycling/rescaling inlet generator (backlog §3a — the documented precondition for Phase 2), then dissolve `modsurfdata` (backlog Phase 2), plus the ledgered dead-array cleanup.

**Architecture:** Pure deletion + relocation, no physics changes: `modinlet.f90` and every generator-only parameter go; `modsurfdata`'s survivors move to their owning modules (`ps` → `modbasestate`; the top-BC cluster `thl_top/qt_top/wttop/wqtop/sv_top/wsvtop/wsvtopdum` → `modboundary`); everything else in `modsurfdata` is deleted with the module file. All gates bitwise.

**Tech Stack:** Fortran 90; bc_cleanup harness (090/092 bitwise, 091 assert-only) — all runs deferred to the Linux session.

## Global Constraints

- No Fortran toolchain here: static verification only (symbol-by-symbol use-list audits, greps); every skipped run recorded. All three tasks expect **bitwise** gates (dead-code deletion / relocation only).
- Conventional commits ending `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`; tick matching `BCcleanup_backlog.md` boxes (§3a list, Phase 2 list) in the same commits.
- Scout facts (verified 2026-07-16 on `BCcleanup` HEAD 34e9747d — line refs below are from that state):
  - `iinletgen` (modglobal.f90:161) is in NO namelist — always 0, generator statically unreachable. `lreadminl` (&INLET, modstartup.f90:142, bcast :423) gates a second dead class (restart inlet-mean reads, modstartup.f90:2290-2351). No namoptions in the repo sets either.
  - `inletgen/inletgennotemp` have zero call sites (not even commented). Live-but-inert call sites to delete: `readinletfile` (modstartup.f90:1444,1878 inside `iinletgen==2`), `exitinlet` (modstartup.f90:2371, unconditional no-op — its body deallocates only under iinletgen guards), zinterpolate1d/2d/w1d/t1d (modstartup.f90:2302-2319 inside the `lreadminl` block), dead import `initinlet` (modstartup.f90:88).
  - Keep (live, misleading names): `Uinf`, `Vinf`, `ifixuinf` (modglobal), `inletav` (modglobal:317, consumed modforces.f90:188), `irecy` (modinletdata:181, consumed modboundary.f90:40,60), `ubulk/vbulk` (modinletdata:130-131, written unconditionally modstartup.f90:1339,1346, read modboundary.f90:159-160), ALL driver data (modinletdata:143-144,194-234; consumed by moddriver/modboundary/modtstep). `modinletdata` itself SURVIVES this phase (its recycling half is deleted; driver merge is Phase 3b, out of scope).
  - `modsurfdata` import sites (complete): modbasestate.f90:20 (ps), modboundary.f90:132,399 (top-BC cluster), modinlet.f90 ×5 (dies in Task 1), modstartup.f90:81,695,949,2166, modthermodynamics.f90:263,397 (ps). NO other file imports it.
  - `thvs` has zero readers (three bridge writes modstartup.f90:1030,1151,1678). `thls` dereferences survive only in modinlet + the dead modstartup branches (1380,1827 in `iinletgen==1`; 2349 in `lreadminl`). Everything else in modsurfdata except `ps` + top-BC cluster + `wsvtopdum` (live at modstartup.f90:507) is import-free dead weight.
  - Dead extras (zero references beyond declaration/alloc): `momfluxb/tfluxb/qfluxb/cth` (modfields.f90:87-90, 849-854, 872 — the modibm `cth` hits are unrelated locals/dummies), `bcTfluxA` (modibmdata.f90:53), and the four dead `rhs` stores in `ibmwallfun` (modibm.f90:1172,1177,1180,1183,1186,1210,1237 — `rhs` is never read; wallfunmom's `rhs` dummy is a different variable).

---

### Task 1: Phase 3a — delete the inlet generator

**Files:** Delete: `src/modinlet.f90`. Modify: `src/modinletdata.f90`, `src/modstartup.f90`, `src/modglobal.f90`, `src/modsave.f90`, `src/modtstep.f90`, `BCcleanup_backlog.md`.

- [ ] **Step 1: `git rm src/modinlet.f90`** (build globs src/*.f90 — no CMake edit).
- [ ] **Step 2: modstartup.f90 surgery** (locate by content, not stale line numbers):
  - Delete imports: `use modinlet, only : initinlet` (:88), `use modinlet, only:readinletfile` (:959), the `use modinlet, only:zinterpolate…` block (:2177), `use modinlet, only:exitinlet` (:2362) and `call exitinlet` (:2371).
  - Delete the `iinletgen==1 / elseif iinletgen==2` startup branches in `readinitfiles` (both copies: ~:1349-1460 and ~:1783-1920) — KEEP the unconditional `ubulk`/`vbulk` assignments that precede them (~:1339,1346) and whatever follows the closing `end if ! iinletgen/idriver`; the `idriver` branches inside the same construct STAY (driver is live). Read the whole construct carefully first; only the `iinletgen` arms go.
  - Delete the `if (lreadminl) then …` block in `readrestartfiles` (~:2290-2351) including its `totinletav` reads and zinterpolate calls.
  - Delete the `use modinletdata` inlet-array imports (:955 list; from :89 remove `di, dr, di_test, dti` — keep `iangledeg, iangle`).
  - &INLET namelist (:140-143): remove `di, dti, linletRA, lfixinlet, lfixutauin, lreadminl` and their broadcasts (:394,395,396,423,527,528); KEEP `Uinf, Vinf, inletav` (+ their broadcasts) — the group survives with only those members.
- [ ] **Step 3: modglobal.f90:** delete declarations `iinletgen` (:161), `lfixinlet, lfixutauin` (:178-179), `linletRA` (:177), `totinletav` (:318). KEEP `inletav` (:317), `Uinf/Vinf/ifixuinf`.
- [ ] **Step 4: modsave.f90:** remove `iinletgen` from the modglobal import (:45), the `use modinletdata, only : nstepread` (:48), and the dead `if ((iinletgen==2) .and. (nstepread==nstore))` block (:58 — read it: delete the whole always-false branch, keep any else/fallthrough).
- [ ] **Step 5: modtstep.f90:** in the monitor-logging block (:289), remove the dead `if (iinletgen == 1)` arm keeping the else body's write as the unconditional path; drop `iinletgen`-adjacent imports (:175) but KEEP `nstepreaddriver, irecydriver` (:181 — live driver reads, :293).
- [ ] **Step 6: modinletdata.f90:** delete the recycling-only inventory — lines 27-129 except nothing (all recycling) minus keeping `di`? NO: `di/dti` leave with the namelist; delete 27-129 wholesale, delete 132-142 (`totalu…totalreadu`) and 146-177 (y-interpolation block) and 183-190 (`nfile…lzinzsim`) EXCEPT keep `irecy` (:181) and `ubulk, vbulk` (:130-131) and ALL driver data (:143-144 `iangle/iangledeg`, :194-234). Re-read the file after editing: every surviving symbol must have a live consumer (irecy→modboundary, ubulk/vbulk→modboundary, driver block→moddriver/modboundary/modtstep, iangle→moddriver).
- [ ] **Step 7: static checks + commit.** `grep -rn "modinlet\b\|iinletgen\|lfixinlet\|lfixutauin\|linletRA\|totinletav\|lreadminl\|readinletfile\|exitinlet\|zinterpolate\|inletgen" src/` → zero (list any deliberate comment survivors). Symbol audit of every touched use list. Commit `refactor(inlet): delete unreachable recycling/rescaling inlet generator (backlog 3a, #68)`; tick §3a boxes; note in §5 that the 3a scope is done.

**Deferred gate:** 090/092 bitwise, 091 assert-only (pure dead-code removal; the deleted `exitinlet` call and modsave/modtstep dead branches had no observable effect).

---

### Task 2: dissolve `modsurfdata`

**Files:** Delete: `src/modsurfdata.f90`. Modify: `src/modbasestate.f90`, `src/modboundary.f90`, `src/modstartup.f90`, `src/modthermodynamics.f90`, `BCcleanup_backlog.md`.

- [ ] **Step 1: move `ps` into `modbasestate`** — declaration `real :: ps = 101325.` with its doc comment; `initbasestate` then uses its own module variable (drop the `use modsurfdata` at :20). Re-point the other `ps` importers: modthermodynamics.f90:263,397 and modstartup.f90:695 (checkinitvalues) and the `ps` entry at modstartup.f90:81 → `use modbasestate, only : ps`. The &PHYSICS namelist read + bcast (:121,:509) keep working — namelist members just need `ps` accessible in `readnamelists`' scope via the new import.
- [ ] **Step 2: move the top-BC cluster into `modboundary`** — declarations `thl_top, qt_top, wttop, wqtop` (with defaults/comments from modsurfdata.f90:62,65,81,83) and allocatables `sv_top, wsvtop` + `wsvtopdum(1:99) = 0.` (:70,87,89). They are consumed in modboundary (:132→ delete that use line; :399 likewise) and set/broadcast in modstartup (:81 import → `use modboundary, only : …`; allocation of wsvtop/sv_top ~:506-508, `wsvtop = wsvtopdum(1:nsv)` :507; sv_top writes :1573-1574,2016-2017 + bcasts). Check modboundary's existing module structure (public/save) and place the declarations with the module-level data. Audit for circular use: modboundary must not `use modstartup` (it doesn't); modstartup already uses modboundary? Check — if not, add the only-list import.
- [ ] **Step 3: delete the rest** — `thvs` (+ the 3 bridge lines modstartup.f90:1030,1151,1678 and `thvs` in the :949 import), `thls` (+ from :949, :2166 imports — the dereferencing dead branches died in Task 1), `qts, z0, z0h, wtsurf, wqsurf, wsvsurf`, and the whole dead block (`tskin, qskin, lmostlocal, obl, oblav, Cm, Cs, ustar, thlflux, qtflux, svflux, dudz, dvdz, dqtdz, dthldz, svs, Cmav, Csav, horvel`). Then `git rm src/modsurfdata.f90`.
- [ ] **Step 4: static checks + commit.** `grep -rn "modsurfdata" src/ tools/ docs/` → zero in src/ (docs mentions of the historical module are fine — list them); every moved symbol resolves at every consumer (walk all former import sites); &BC and &PHYSICS namelist members all in scope in `readnamelists`. Commit `refactor(surface): dissolve modsurfdata; ps to modbasestate, top BCs to modboundary (backlog phase 2)`; tick the Phase 2 boxes.

**Deferred gate:** 090/092 bitwise, 091 assert-only (relocation only).

---

### Task 3: dead-array extras + input-surface sweep

**Files:** Modify: `src/modfields.f90`, `src/modibmdata.f90`, `src/modibm.f90`, `tools/python/namelists.json`, `docs/udales-namoptions-overview.md`, `BCcleanup_backlog.md`; possibly namoptions under tests/examples.

- [x] **Step 1:** delete `momfluxb/tfluxb/qfluxb/cth` (modfields.f90:87-90 decls, 849-852 allocs, 854 zero-init, 872 deallocate — adjust the deallocate list, don't orphan it) and `bcTfluxA` (modibmdata.f90:53). Verify with grep that the only `cth` survivors are modibm's unrelated locals/dummy args.
- [x] **Step 2:** delete the dead `rhs` in `ibmwallfun` (modibm.f90: declaration :1172, allocate :1177, the four `rhs = up/vp/wp/thlp` stores, deallocate :1237). Do NOT touch `wallfunmom`'s `rhs` dummy argument — different variable.
- [x] **Step 3: input surfaces.** `grep -rn "^\s*\(di\|dti\|linletRA\|lfixinlet\|lfixutauin\|lreadminl\)\s*=" tests examples` and strip hits (expect zero — verify). tools/python/namelists.json: remove those keys from the "INLET" list + reverse mappings (keep uinf/vinf/inletav); validate with `python -m json.tool`. docs/udales-namoptions-overview.md: remove/mark the &INLET rows for the deleted keys; note the inlet generator's removal (points at driver/synthetic inflow as the supported mechanisms). **Plus three review-mandated additions:** (a) deleted write-only `irecy` (modinletdata.f90 decl, modboundary.f90 write + import); (b) deleted dead `lstoreplane`/`lwallfunc` (modglobal.f90 decls, &INLET namelist, MPI_BCASTs) — `&INLET` now carries only `Uinf, Vinf, inletav`; (c) documented the pre-existing `wqtop`-absent-from-&BC gap in `BCcleanup_backlog.md` §5.
- [x] **Step 4: commit** `refactor(fields): drop dead wall-flux arrays and rhs stores; sweep retired &INLET keys`; tick the remaining ledgered-minor boxes in the backlog.

**Deferred gate:** 090/092 bitwise, 091 assert-only.

---

### Task 4: final Phase-2 branch review

Whole-phase review (base = Phase 1 HEAD 34e9747d) on the strongest model; one fix wave; re-verify; close the ledger with the updated deferred-gate spans (all Phase 2 commits expected bitwise on 090 AND 092).
