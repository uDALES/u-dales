# BC Cleanup Phase 1 Implementation Plan — retire `lbottom` and the flat-wall scheme

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Delete the legacy flat-surface scheme (`subroutine bottom`, `lbottom`, `BCbot*`, wall-function `CASE(91)/(92)`), relocate its unconditional code, retire its dead namelist keys and output fields, and migrate the two `lbottom` cases to ground facets (Phase 1 of `BCcleanup_backlog.md`).

**Architecture:** The unconditional e12 ghost moves into `closurebc` (next to the existing `ekm/ekh` ghosts — deliberately fixing a one-step staleness). The always-zero `tau_x/tau_y/tau_z/thl_flux` fielddump fields are retired. The last live `thls` consumer (`tkestatsdump`) is re-pointed to `thv_b(kb)`, after which `thls, qts, z0, z0h, wtsurf, wqsurf, wsvsurfdum, BCbot*, lbottom` all leave the namelists (variables stay declared in `modsurfdata` until Phase 2). Cases 103 and 999 are regenerated with ground facets via the Python preprocessing.

**Tech Stack:** Fortran 90, existing bc_cleanup harness (cases 090/092), Python preprocessing (venv at `../.venv`, as used for case 091).

## Global Constraints

- No Fortran toolchain on this machine: builds and harness gates are **deferred to the Linux session** and recorded in `.superpowers/sdd/progress.md`. Static verification (symbol-by-symbol use-list audit, greps) is mandatory per task.
- Comparison regimes per task, stated below; never claim bitwise without a passing `--atol 0.0` run (deferred).
- Fortran style: lowercase, `use <module>, only : <list>`; comments state constraints, not narration.
- Every commit ends with `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
- Tick the matching `BCcleanup_backlog.md` Phase 1 checkboxes in the same commit as each item.
- Scout facts (verified 2026-07-16 on `BCcleanup` HEAD): `bottom` = modibm.f90:1997-2099 with unconditional ghost lines 2010-2011 and tau bracket 2014-2017/2093-2096; `call bottom` program.f90:154, import :39; CASE(91)=modwallfunctions.f90:81-130 & 307-349, CASE(92)=133-162, called only from `bottom` (modibm.f90:2023,2025,2046); BCbot* decls modglobal.f90:158-174, namelist modstartup.f90:136, bcasts :422-425, stray write :828; `lbottom` decl modibm.f90:49, namelist &WALLS modstartup.f90:153, bcast :445, validation :774-779; e120(kb-1) sole consumer = `diffe` modsubgrid.f90:659 (loneeqn only); ekm/ekh ghosts = `closurebc` modboundary.f90:447-462; tau fields modfields.f90:82,476-479, read only by modfielddump.f90:229-240 (fieldvars codes 'tx','ty','tz','hf'), always zero when `lbottom=.false.`; live `thls` consumer = modstatsdump.f90:2126 (`ltkedump`); `qts`/`wsvsurfdum` already dead.

---

### Task 1: Relocate the e12 ghost; retire the tau/thl_flux dump fields

**Files:**
- Modify: `src/modboundary.f90` (`closurebc`), `src/modibm.f90` (`bottom` unconditional lines), `src/modfields.f90`, `src/modfielddump.f90`, `docs/udales-namoptions-overview.md` (fieldvars codes), `BCcleanup_backlog.md`

**Interfaces:**
- Produces: `bottom` contains ONLY the `lbottom`-guarded body (Task 2 deletes it whole). `closurebc` owns the e12 kb-1 ghost.

- [ ] **Step 1: Move the ghost into `closurebc`**

In `src/modboundary.f90`, `subroutine closurebc`, after the top/bottom `ekm/ekh` if-block (after line ~470, i.e. once, not per-branch), add:

```fortran
      ! e12 ghost at kb-1: zero-gradient mirror for the TKE diffusion stencil at kb.
      ! Relocated from the retired flat-surface scheme; setting it here makes the
      ! value current within the same subgrid evaluation (previously one step stale,
      ! set in subroutine bottom after diffe had already consumed it).
      e120(:, :, kb - 1) = e120(:, :, kb)
      e12m(:, :, kb - 1) = e12m(:, :, kb)
```

Add `e120, e12m` to `closurebc`'s `use modfields` list (check the existing list first). Delete modibm.f90:2010-2011 (the two ghost lines in `bottom`).

- [ ] **Step 2: Retire `tau_x/tau_y/tau_z/thl_flux`**

They only ever record `bottom`'s tendency contribution and are identically zero for `lbottom=.false.` (modibm.f90:2093-2096 computes `up-up` etc.). Delete:
- modibm.f90:2014-2017 and 2093-2096 (the bracket in `bottom`) and the four names from `bottom`'s `use modfields` list (modibm.f90:2004);
- declarations modfields.f90:82 and allocations modfields.f90:476-479 (and any matching deallocate);
- the fielddump cases modfielddump.f90:229-240 (the `'tx'/'ty'/'tz'/'hf'` cases — read the surrounding `select`/`if` chain and remove the whole case blocks) and the four names from modfielddump.f90:64;
- any mention of `tx/ty/tz/hf` fieldvars codes in `docs/udales-namoptions-overview.md`.
Check how modfielddump handles an unrecognized fieldvars code (read the chain's else/default); if it errors loudly, that is the desired behaviour for old namoptions requesting 'tx' — say so in the commit message.

- [ ] **Step 3: Static checks + commit**

`grep -n "tau_x\|tau_y\|tau_z\|thl_flux" src/` → zero matches. `grep -n "e120(:, :, kb - 1)\|e12m(:, :, kb - 1)" src/` → only modboundary.f90. Commit `refactor(ibm): relocate e12 kb-1 ghost to closurebc; retire tau/thl_flux dump fields`. Tick the backlog "Relocate the unconditional code" boxes.

**Deferred gates:** 090 bitwise (Vreman: `diffe` never runs, tau fields not dumped); 092 ≤ 5e-3 (loneeqn: ghost staleness fix is a real, deliberate behaviour change — record max diffs).

---

### Task 2: Delete `bottom`, CASE(91)/(92), `BCbot*`, `lbottom`; re-point `tkestatsdump`; prune the namelists

**Files:**
- Modify: `src/modibm.f90`, `src/program.f90`, `src/modwallfunctions.f90`, `src/modglobal.f90`, `src/modstartup.f90`, `src/modstatsdump.f90`, `BCcleanup_backlog.md`

- [ ] **Step 1: Re-point the last live `thls` consumer**

modstatsdump.f90:2126: replace `(grav/thls)` with `(grav/thv_b(kb))`; add `use modbasestate, only : thv_b` to `tkestatsdump` (check the routine's existing imports; remove `thls` from its `use modsurfdata` list, deleting the line if empty). This pulls forward the backlog's Phase-2 re-point; it changes `ltkedump` output values by design (better reference) — no registered gate case uses `ltkedump`.

- [ ] **Step 2: Delete the scheme**

- `subroutine bottom` in full (modibm.f90:1997-2099 as it stands after Task 1), `call bottom` (program.f90:154) and `bottom` from the import at program.f90:39; `lbottom` declaration + `public` export (modibm.f90:49, :30).
- CASE(91) and CASE(92) blocks in `wfuno` (modwallfunctions.f90:81-130, 133-162) and CASE(91) in `wfmneutral` (:307-349), including the "should be moved out" comments.
- `BCbot*`: declarations + named constants (modglobal.f90:158-174), namelist membership (modstartup.f90:136), broadcasts (:422-425), and the stray no-op write `BCbotm = BCbotm_wfneutral` at modstartup.f90:828 (KEEP the `iwallmom = 3` flip on the previous line — it is live IBM behaviour).
- `lbottom` namelist membership (&WALLS, modstartup.f90:153), import (:95-98 list), broadcast (:445), and the Phase-0 validation block (modstartup.f90:774-779).

- [ ] **Step 3: Prune the &BC namelist keys whose consumers are now gone**

Remove from `namelist/BC/` (modstartup.f90:132-140) and their reads/broadcasts: `thls`, `qts`, `z0`, `z0h`, `wtsurf`, `wqsurf`, `wsvsurfdum` (plus the `wsvsurf = wsvsurfdum(1:nsv)` assignment, its allocation/broadcast at modstartup.f90:516-518). KEEP `ds`, `BCzp`, `wttop/thl_top/qt_top/wsvtopdum` and all other top-BC keys (live). The modsurfdata declarations stay (Phase 2); unreachable inlet-generator reads of `thls`/`z0` still compile against the declared defaults. Remove now-unused imports of the pruned names from modstartup's `use modsurfdata` list at :81 and :966 — but verify each name is truly unreferenced in the file first (e.g. `ps` stays).

- [ ] **Step 4: Static checks + commit**

`grep -n "lbottom\|BCbot" src/` → zero matches. `grep -n "wforient" src/modwallfunctions.f90` → no 91/92 cases remain. `grep -n "thls\|qts\b\|wtsurf\|wqsurf\|wsvsurf" src/modstartup.f90` → only unreachable inlet branches (list them in the report). Compile-risk audit of every touched `use` list. Commit `refactor(bc): delete flat-surface scheme (bottom, BCbot*, lbottom, wf CASE 91/92); prune dead BC keys`. Tick the backlog boxes (delete bottom / CASE(91,92) / BCbot* / migrate-lbottom precursor).

**Deferred gates:** 090 AND 092 bitwise (everything deleted was `lbottom`-guarded or dead; the statsdump re-point is `ltkedump`-gated and off in both cases).

---

### Task 3: Sweep the case files, preprocessing surfaces, and docs

**Files:**
- Modify: every namoptions under `tests/` and `examples/` that sets a pruned key; `tools/preprocessing.m`; `tools/python/namelists.json`; `docs/udales-namoptions-overview.md`; `docs/udales-example-simulations.md`; `tests/regression/bc_cleanup/README.md`; `tests/cases/090|091|092/README.md` if they mention thls/qts

- [ ] **Step 1: Strip pruned keys from all namoptions**

`grep -rln "^\s*\(thls\|qts\|z0\|z0h\|wtsurf\|wqsurf\|BCbot[mTqs]\|lbottom\)\s*=" tests examples` and delete those lines (whole-word; do NOT touch `z0` inside factypes/facet files, only namoptions `&BC`/`&WALLS` lines; do NOT touch cases 103/999 — Tasks 4-5 rewrite them wholesale). Scout counts for cross-check: thls 7 files, qts 5, z0/z0h 10, wtsurf/wqsurf 12. List every file in the commit body.

- [ ] **Step 2: Preprocessing + docs surfaces**

- `tools/preprocessing.m:224`: delete the `addvar(obj, 'lbottom', 0)` line and the dead commented block at ~:1214.
- `tools/python/namelists.json`: remove `"bcbotm","bcbotq","bcbots","bcbott"` (lines 4-7), `"thls"` (:35), `"qts"` (:33), and the z0/z0h/wtsurf/wqsurf/wsvsurfdum entries from the "BC" list; remove `"lbottom"` from the "WALLS" list (:282) and the `"lbottom": "WALLS"` (:393), `"qts": "BC"` (:502), `"thls": "BC"` (:518) mappings (and the matching mappings for the other pruned keys). Validate the file with `python -m json.tool`.
- `docs/udales-namoptions-overview.md`: delete the thls/qts/z0/z0h/wtsurf/wqsurf/BCbot*/lbottom rows (they were marked deprecated in Phase 0); note in the &BC section that legacy flat-surface keys were removed and old namoptions fail loudly on them.
- `docs/udales-example-simulations.md` (~line 198): replace the lbottom description with the ground-facet approach.
- Update the bc_cleanup READMEs: cases 090/092 no longer set thls/qts (the matched-anchor construction is now historical — the base state derives from prof.inp alone); keep the tolerance table.

- [ ] **Step 3: Commit**

`test/docs: strip retired flat-surface keys from cases, preprocessing and docs`. Tick the backlog preprocessing/docs boxes.

**Deferred gates:** 090/092 bitwise (namoptions key removal only — but this REQUIRES Task 2's code to be in the same build; the harness runs each ref's own case files, so cross-ref comparisons stay valid).

---

### Task 4: Migrate `tests/regression/david_tests/cases/103` to ground facets

**Files:**
- Rewrite: `tests/regression/david_tests/cases/103/` inputs

- [ ] **Step 1: Generate ground-facet inputs**

Case 103 is 8×8×8, flat, `libm=.false.`, `lbottom`, `BCbotT=2` (wall-function T with `thls=288`), `z0=0.01`, `z0h=0.000067`, `ltempeq/lmoist=.true.`. Using the Python preprocessing (same venv/driver as case 091 — see `.superpowers/sdd/task-7-report.md` for the exact invocation): generate a **flat ground-facet** floor at z=0 (no volume blocks), producing facets/solid/fluid/facet-section files. Facet type: roughness `z0=0.01`, `z0h=0.000067` via `factypes.inp.103`; initial facet temperature 288 K via `Tfacinit.inp.103` (file already exists — regenerate consistently).

- [ ] **Step 2: Rewrite `namoptions.103`**

`libm=.true.`, `&WALLS` counts from the generated files, `iwallmom=2`, `iwalltemp=2` (facet-temperature wall function — nearest equivalent of the old `BCbotT=2`); delete `lbottom`, `BCbotT`, `thls`, `qts`, `z0`, `z0h`, `wtsurf`, `wqsurf`. Keep grid/runtime/physics unchanged.

- [ ] **Step 3: Static verification + commit**

Parse the generated solid/fluid files with python: floor facets must cover the full 8×8 floor; NO fully-solid slab may exist (ground facets at z=0 leave kb fluid). Commit `test: migrate case 103 from lbottom to ground facets (backlog phase 1)`. Record in the backlog §6.2 note: statistical-equivalence run (old-vs-new 103, slab profiles over an averaging window) deferred to the Linux session; the wall-function formulation legitimately changes, so the criterion is statistical, not bitwise.

---

### Task 5: Migrate `examples/999` to ground facets

Same procedure as Task 4 for `examples/999` (128³, neutral — `ltempeq` off, only `z0=0.05`, `z0h=0.00035`): flat floor facets with that roughness, `libm=.true.`, `nfcts` from generation, delete `lbottom`/`z0`/`z0h`; keep everything else. Static check as Task 4 (floor coverage 128×128, no solid slab). Commit `docs(examples): migrate example 999 from lbottom to ground facets`. Tick the backlog §4 migration boxes.

---

### Task 6: Final Phase-1 branch review

Whole-branch review of the Phase 1 commits (base = Phase 0 HEAD `e5585afe`) on the strongest model, findings fixed in one wave, ledger updated. Then update the deferred-gates list in `.superpowers/sdd/progress.md`: 090 bitwise across all Phase-1 commits; 092 ≤5e-3 for Task 1 then bitwise for Tasks 2-3; §6.2 statistical for 103/999.
