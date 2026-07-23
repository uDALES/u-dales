# BC Cleanup Phase 3b Implementation Plan — the `inflow` module

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Merge `moddriver` + the surviving `modinletdata` into a single self-contained module **`inflow`** (file `src/inflow.f90` — user-directed name, no `mod` prefix), relocate `ubulk/vbulk` to `modboundary`, rename the selector `idriver` → `iinflow` (namelist `&DRIVER` → `&INFLOW`), and clear the Phase-2 carry-overs (backlog §3b).

**Architecture:** `inflow` = modinletdata's driver data (module-level) + moddriver's subroutines, with module lifecycle entry points `initinflow`/`exitinflow` (dispatching on `iinflow`) and the method-specific routines (`drivergen`, `readdriverfile`, `readdriverfile_chunk`, `driverchunkread`, `writedriverfile`) keeping their names — the driver *method* is still legitimately called "driver"; `iinflow` selects the method (0 = none, 1 = write precursor driver planes, 2 = driver inflow from files; 3 reserved for the future synthetic generator). `ubulk/vbulk` go to `modboundary` (their only reader), NOT into `inflow` — they are outflow-BC state, not inflow-generation state (backlog §3b). Pure relocation + rename: all gates bitwise.

**Tech Stack:** Fortran 90; bc_cleanup harness — runs deferred to the Linux session; the three driver cases (examples/949, examples/950, tests/cases/525) are the regression anchors for the merge and get their namoptions migrated.

## Global Constraints

- No Fortran toolchain: static verification only; record skipped runs. Both code tasks expect **bitwise** gates on 090/091/092 (no driver active there) — the only behavioural addition is `call exitinflow` wired into `exitmodules` (teardown-only deallocation, no observable effect).
- Conventional commits + `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`; tick backlog §3b boxes in the same commits.
- DO NOT TOUCH (name look-alikes, scout-verified): `RUN_DRIVER` runmode constant (modglobal.f90:75, program.f90:29,242); the `BCxm_driver`/`BCxT_driver`/`BCxq_driver`/`BCxs_driver` BC-type enums; `driverid/cdriverid`; `linoutflow/luoutflowr/luvolflowr` family; `lsdriver/lhdriver/lqdriver` and the other `&DRIVER`-adjacent modglobal config (`tdriverstart`, `dtdriver`, `driverstore`, `driverjobnr`, `lchunkread`, `chunkread_size`, `iplane`) — that config STAYS in modglobal per repo convention; only its namelist group is renamed.
- Scout facts (2026-07-16, HEAD 27056eff): moddriver.f90:32-38 has NO module data (blanket `use modinletdata` at :33), public list :36 (`initdriver,exitdriver,readdriverfile,drivergen,readdriverfile_chunk,driverchunkread`); subroutines initdriver :39-172, drivergen :174-513, writedriverfile :515-748, readdriverfile :750-931, readdriverfile_chunk :933-1165, driverchunkread :1167-1188, exitdriver :1190-1219 (NEVER called — wire `exitinflow` into exitmodules). Call sites: program.f90:88 (`call initdriver`, unconditional, self-guarded), modboundary.f90:256,270-271 (drivergen/driverchunkread), modstartup.f90:1331-1335,1388,1600,1653-1657 (drivergen/readdriverfile[_chunk]). `idriver` decl modglobal.f90:160; auto-set `idriver = 2` + re-broadcast in checkinitvalues modstartup.f90:828-832 under `case(BCxm_driver)`; &DRIVER namelist modstartup.f90:141-144 (9 members), read :238, bcasts :396-402,455,518. modinletdata survivors: ubulk/vbulk (:27-28), iangle/iangledeg (:29-30), driver arrays (:34-66), irecydriver/nstepreaddriver/chunkread counters (:68-73). Suspected-dead members (grep-verify each before deleting; keep any with a live ref): `tdriver` (:48, never allocated), `storee120driver` (:38), `e120driver` (:47), `storee12mdriver` (:59), `e12mdriver` (:60), `storeumdriver` (:53). ubulk/vbulk writers modstartup.f90:1318,1325,1363,1739; readers modboundary.f90:166-167. modinletdata importers: moddriver (blanket + 4 only-lists), modboundary :142,736,812,843,903,1214, modstartup :88,940, modtstep :181.

---

### Task 1: create `inflow`; retire `moddriver` + `modinletdata`

**Files:** Create: `src/inflow.f90`. Delete: `src/moddriver.f90`, `src/modinletdata.f90`. Modify: `src/modboundary.f90`, `src/modstartup.f90`, `src/modtstep.f90`, `src/program.f90`, `BCcleanup_backlog.md`.

- [ ] **Step 1: build `src/inflow.f90`** — module `inflow`, `implicit none`, `save`, explicit public list. Contents:
  - Module data: everything from modinletdata EXCEPT `ubulk/vbulk` (→ modboundary) and any suspected-dead member you PROVE dead by repo grep (list the evidence per deleted name in your report; keep anything with a single live reference). Keep declarations/comments/defaults byte-identical for everything moved.
  - Subroutines: all seven from moddriver, verbatim bodies, with two renames: `initdriver` → `initinflow`, `exitdriver` → `exitinflow` (update their self-guards' comments accordingly; internal `idriver` reads stay — Task 2 renames the selector). Drop moddriver's blanket `use modinletdata` and the now-internal redundant only-lists (the data is host-associated from the module itself); PRESERVE each subroutine's other use statements exactly (note: drivergen's `use modsave, only: writerestartfiles` is a dead import — drop it and say so).
  - Public: `initinflow, exitinflow, drivergen, readdriverfile, readdriverfile_chunk, driverchunkread` + every datum other modules import (`u0driver…svmdriver` set for modboundary, `iangle, iangledeg` for modstartup, `nstepreaddriver, irecydriver, chunkread…` for modtstep — enumerate from the importer lists and make them public explicitly).
- [ ] **Step 2: `ubulk`/`vbulk` → modboundary** module data (defaults `0.`, comments preserved), added to its public list; delete modboundary's `use modinletdata, only : ubulk, vbulk` (:142) — module-local now; re-point modstartup:940 to `use modboundary, only : ubulk, vbulk`.
- [ ] **Step 3: re-point all importers** — modboundary :736,812,843,903,1214 and modtstep :181 → `use inflow, only : …` (same symbol lists); modstartup :88 → `use inflow, only : iangledeg, iangle`; modboundary's `use moddriver, only : drivergen, driverchunkread`-style imports (find them near :256,270) → `use inflow, …`; modstartup's driver-call imports near :1331/:1653 likewise; program.f90:88 `call initdriver` → `call initinflow` (+ its import). Wire `call exitinflow` into `exitmodules` (modstartup) alongside the other exit calls, with the module import.
- [ ] **Step 4: delete the two old files** (`git rm`); CMake globs — no build edit.
- [ ] **Step 5: static checks + commit.** `grep -rn "moddriver\|modinletdata" src/` → zero (list comment survivors — the moddriver.f90:7/:1194 fossils die with the file); symbol-by-symbol audit: every `use inflow` only-list member is declared public in the module; every former importer compiles by inspection. Commit `refactor(inflow): merge moddriver + modinletdata into new inflow module (backlog 3b)`; tick the §3b merge/relocate/delete boxes.

**Deferred gate:** 090/091/092 bitwise; driver anchors 949/950/525 build-and-run deferred (statistical identity expected — pure relocation).

---

### Task 2: internal carry-overs only — REVISED 2026-07-16 per user decision

> **User decision:** NO renaming of `idriver` or `&DRIVER` — "no changes for the user at this
> stage"; the selector consolidation (`iinflow`) is DEFERRED to a future dedicated task. Steps
> 1-4 below are superseded; only Step 5's carry-overs plus the orphaned-INFO-group cleanup in
> namelists.json (tooling-internal, describes a namelist that no longer exists) remain in scope.

### ~~Task 2 (superseded): `iinflow` selector + `&INFLOW` group + carry-overs~~

**Files:** Modify: `src/modglobal.f90`, `src/modstartup.f90`, `src/inflow.f90`, `src/modboundary.f90`, `src/modtstep.f90`, `src/modsave.f90`, 9 namoptions files, `tools/python/namelists.json`, `tools/preprocessing.m`, `docs/udales-driver-simulations.md`, `docs/udales-namoptions-overview.md`, `BCcleanup_backlog.md`.

- [ ] **Step 1: rename the selector.** modglobal.f90:160 → `integer :: iinflow = 0 !< inflow generation: 0 = none, 1 = write precursor driver planes, 2 = driver inflow from files (3 reserved: synthetic)`. Rename every `idriver` read/write repo-wide in src/ (scout table: modboundary :135,:256; inflow.f90 [ex-moddriver] ×10; modstartup :76,:142,:396,:698,:830-832,:933,:1327,:1369,:1559,:1649 + comments; modtstep :175,:289). The checkinitvalues auto-set becomes `iinflow = 2` under `case(BCxm_driver)` with its re-broadcast. Verify with `grep -rn "idriver" src/` → zero (excluding `driverid/cdriverid/lsdriver/lhdriver/lqdriver` — use word-boundary grep `\bidriver\b`).
- [ ] **Step 2: rename the group.** modstartup.f90:141-144 → `namelist/INFLOW/ iinflow, tdriverstart, driverjobnr, dtdriver, driverstore, iplane, iangledeg, lchunkread, chunkread_size`; update the read statement's group + error message (:238-243).
- [ ] **Step 3: sweep the inputs.** All 9 namoptions with `&DRIVER` (examples/949, examples/950, tests/cases/525 active; tests/cases/526, tests/integration/processor_boundaries/namoptions.526.serial, 2× mpi_averaging_regression namelists, 2× new_vegetation_module namelists inert with idriver=0): `&DRIVER` → `&INFLOW`, `idriver` → `iinflow`, values unchanged. namelists.json: group `"DRIVER"` → `"INFLOW"` with `idriver` → `iinflow` in the list + all 9 reverse-map entries re-pointed; DELETE the orphaned `"INFO"` group (:100-107) and its 6 reverse-map entries (dtin, jgtotinl, kmaxin, nprocsinl, totalreadu, wtop); validate with `python -m json.tool`. preprocessing.m :232 (`addvar 'idriver'` → `'iinflow'`) and :234 (`obj.idriver~=2` → `obj.iinflow~=2`).
- [ ] **Step 4: docs.** `docs/udales-driver-simulations.md`: retitle to inflow simulations framing; `&DRIVER`→`&INFLOW`, `idriver = 1/2` → `iinflow = 1/2` throughout (keep the driver-method terminology in prose — the files are still `*driver_*`); note old namoptions fail loudly. `docs/udales-namoptions-overview.md` :237-247 (section header `Namelist DRIVER` → `Namelist INFLOW`, the idriver row → iinflow) and :359 (`idriver` mention).
- [ ] **Step 5: carry-overs** (backlog §3b list): delete the write-only locals `uaverager/uaveragei/taverager/taveragei/waverage` (modstartup.f90:950-954 decls + :1644-1648 zero-writes); remove the dead `Uinf` from readinitfiles' modglobal only-import (:928) — KEEP `jgb, jge` (live at :1734, scout-verified); fix modsave.f90:43-45 duplicate `timee`. (INFO group handled in Step 3; moddriver fossils died in Task 1.)
- [ ] **Step 6: static checks + commit.** `\bidriver\b` and `&DRIVER`/`namelist/DRIVER` greps → zero repo-wide (code + inputs + tools + docs; list historical-doc survivors); json valid; commit `refactor(inflow)!: rename idriver to iinflow, &DRIVER to &INFLOW; clear 3b carry-overs` (breaking-change marker for the input rename); tick the remaining §3b boxes.

**Deferred gate:** 090/091/092 bitwise (no &DRIVER/&INFLOW group set in them — verify; if any sets the inert block, migrate it in the sweep); driver cases parse-checkable only on Linux.

---

### Task 3: final Phase-3b review

Whole-phase review (base = 27056eff) on the strongest model; fix wave if needed; ledger closed with gate spans; backlog §3/§5 updated (Phase 3a+3b complete — remaining §5 out-of-scope items stand alone); memory update.
