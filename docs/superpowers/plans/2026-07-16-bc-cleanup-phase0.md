# BC Cleanup Phase 0 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the `thls`/`qts`-anchored thermodynamic base state with one derived from `prof.inp` + `ps`, delete the flat-surface `kb` overwrite, and lock the behaviour with a branch-comparison regression harness (Phase 0 of `BCcleanup_backlog.md`, issues #299/#302).

**Architecture:** A new `src/modbasestate.f90` builds fixed base-state profiles (`thl_b/qt_b/thv_b/pf_b/ph_b/exnf_b/exnh_b`) at every start from the initial profiles and `ps`. `modthermodynamics` consumers are re-pointed to it; `modsubgrid` moves to the evolving `thvf(k)` (DALES ≥ 4.0 parity). Each commit is gated by a NetCDF field comparison between the pre-change and post-change branch states, with a per-commit tolerance regime (bitwise where the backlog proves neutrality, tolerance otherwise).

**Tech Stack:** Fortran 90 (existing module style), CMake via `tools/build_executable.sh`, Python 3 + numpy + netCDF4 for the harness (pattern: `tests/regression/mpi_averaging_regression/run_test.py`).

## Global Constraints

- Build/run requires Linux or WSL (`docs/udales-installation.md:11`); native Windows is unsupported. Build: `./tools/build_executable.sh common debug` (or `release`) from the repo root. Run: `mpiexec -n <N> <repo>/build/<type>/u-dales namoptions.<case>` with cwd = case directory.
- Branch for all work: `BCcleanup`. Reference for comparisons: the commit tagged at the start of each task (the harness takes `branch_a`/`branch_b` git refs).
- Comparison regimes (backlog §6): **bitwise** = `--atol 0.0`; **tolerance** = stated per task. Never claim bitwise without a passing `--atol 0.0` run.
- Fortran style: lowercase keywords, `use <module>, only : <list>`, validation via the `checkinitvalues` patterns (`src/modstartup.f90:764-767, 792-802`). Match surrounding comment density; comments state constraints, not narration.
- Every commit message ends with `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
- Executors: use a cheaper model (sonnet) for these tasks; the orchestrator reviews between tasks.
- Update the matching checkbox in `BCcleanup_backlog.md` (Phase 0 section) in the same commit as each completed item.

---

### Task 1: Baseline case 090 + bc_cleanup regression harness

**Files:**
- Create: `tests/cases/090/` (copied from `tests/cases/101/`, suffixes renamed)
- Create: `tests/regression/bc_cleanup/run_test.py`
- Create: `tests/regression/bc_cleanup/README.md`

**Interfaces:**
- Produces: CLI `python tests/regression/bc_cleanup/run_test.py <branch_a> <branch_b> <build_type> [--atol 0.0] [--workdir DIR]`, exit 0 on pass, non-zero on any field mismatch above atol. Later tasks call exactly this.
- Produces: case `090` — dry, `ltempeq=.true.`, IBM buildings, **matched anchors**: `thls = 290.` (equal to `prof.inp` bottom `thl`), `qts = 0.0` (equal to `prof.inp` `qt`), `ps` set explicitly. All Phase 0 comparisons run on this case.

- [ ] **Step 1: Copy case 101 → 090 and rename suffixes**

```bash
cd tests/cases && cp -r 101 090 && cd 090
for f in *.101; do mv "$f" "${f%.101}.090"; done
# STL filename referenced in namoptions keeps its own name; check below.
```

- [ ] **Step 2: Edit `namoptions.090`**

Read `namoptions.090` first. Apply these settings (add keys if absent, in the groups shown):
- `&RUN`: `iexpnr = 090`, `runtime = 20.` (short; case 101 default is 101 s)
- `&BC`: `thls = 290.` (was 295 — must equal `prof.inp.090` bottom thl), `qts = 0.0`
- `&PHYSICS`: ensure `ps` is set; if absent add `ps = 101500.`
- `&OUTPUT` (or the fielddump group used by 101 — check `docs/udales-namoptions-overview.md` for the exact keys): enable the instantaneous/`tdump` field output with a sample interval ≤ 5 s so at least 3 dump slices exist within 20 s.
- Keep `libm`, geometry (`stl_file`), and grid unchanged.

Verify `prof.inp.090` bottom row has `thl = 290.0`, `qt = 0.0` (case 101's does: `290.0`, `0.0` at z = 0.5).

- [ ] **Step 3: Write the harness**

Copy `tests/regression/mpi_averaging_regression/run_test.py` to `tests/regression/bc_cleanup/run_test.py`, then apply these changes (keep the worktree/build/run/stitch machinery intact):

Replace the `CASES` tuple with:

```python
CaseSpec = collections.namedtuple("CaseSpec", "case nc_pattern fields abs_tol")

CASES = (
    CaseSpec(
        case="090",
        nc_pattern="tdump.*.090.nc",   # match the per-tile dump files case 090 produces
        fields=None,                    # None = compare every variable with ndim >= 3
        abs_tol=None,                   # filled from --atol
    ),
)
```

Replace `CONFIGS` with just two decompositions:

```python
CONFIGS = {
    "serial":   dict(nprocx=1, nprocy=1),
    "xy_split": dict(nprocx=2, nprocy=2),
}
```

In `_compare_fields`, when `spec.fields is None`, iterate `sorted(ds.data_vars)` and compare every variable with `ds[name].ndim >= 3`; report per-field `max|diff|` and fail if any exceeds `abs_tol`. Add to the argparser:

```python
parser.add_argument("--atol", type=float, default=0.0,
                    help="max allowed |current - reference| per field (0.0 = bitwise)")
```

and thread it into each `CaseSpec` before comparison. Point the case directory to `tests/cases/090` and the run command to `mpiexec -n <nprocs> <exe> namoptions.090`. Delete the tree-case (`526`) handling and the `treedump` logic.

- [ ] **Step 4: Sanity-run the harness A == B**

```bash
python tests/regression/bc_cleanup/run_test.py BCcleanup BCcleanup Debug --atol 0.0
```

Expected: builds one worktree per ref (identical), runs 090 under `serial` and `xy_split`, all fields identical, exit 0. If the executable fails to run case 090, fix the namoptions (most likely the dump group keys) before proceeding — do not loosen atol.

- [ ] **Step 5: Write `tests/regression/bc_cleanup/README.md`**

State: purpose (Phase 0 gate for BCcleanup_backlog.md), the matched-anchor construction of case 090 (thls == prof bottom, qts == 0 so old `thvs` equals new `thv_b(kb)`), the per-commit tolerance regimes table (Task 2: 0.0; Task 3: 0.0; Task 4: 5e-3; Task 5: 0.0), and the usage line from Step 4.

- [ ] **Step 6: Commit**

```bash
git add tests/cases/090 tests/regression/bc_cleanup
git commit -m "test: add bc_cleanup regression harness + matched-anchor case 090"
```

---

### Task 2: Delete the `kb` overwrite (backlog Phase 0, Path A — bitwise commit)

**Files:**
- Modify: `src/modthermodynamics.f90` (`calc_halflev`, lines ~508-540)

**Interfaces:**
- Consumes: Task 1 harness.
- Produces: `thl0h/qt0h(kb)` become one-sided copies of the full-level `kb` value; no other field changes.

- [ ] **Step 1: Re-verify the consumer chain on the current tree**

The backlog §1.3 claims `thl0h/qt0h(kb)` reach no prognostic field. Re-verify on HEAD (the #323 merge reshuffled code):

```bash
grep -n "thl0h\|qt0h\|thv0h\|ql0h" src/*.f90 | grep -v "!"
```

Confirm: `thv0h` is consumed only by the buoyancy loop (`src/modforces.f90:83`, loop `kb+1..ke`) and produced by `calthv`; `ql0h` is produced by `thermo` and only written to restart files (`src/modsave.f90`); no advection/diffusion/statistics path reads `thl0h/qt0h/thv0h/ql0h` at `kb`. If anything new consumes them at `kb`, STOP and report — do not proceed.

- [ ] **Step 2: Replace the overwrite**

In `calc_halflev` (`src/modthermodynamics.f90`), replace:

```fortran
        thl0h(ib:ie,jb:je,kb) = thls
```
and (a few lines below)
```fortran
          qt0h(ib:ie,jb:je,kb)  = qts
```
with:
```fortran
        ! one-sided: the IBM owns the bottom interface; kb may lie inside solid (#299)
        thl0h(ib:ie,jb:je,kb) = thl0(ib:ie,jb:je,kb)
```
```fortran
          qt0h(ib:ie,jb:je,kb)  = qt0(ib:ie,jb:je,kb)
```

Remove `thls`/`qts` from the `use modsurfdata` list of `calc_halflev` (line ~514); delete the whole `use` line if empty.

- [ ] **Step 3: Build**

```bash
./tools/build_executable.sh common debug
```
Expected: clean compile, no new warnings.

- [ ] **Step 4: Bitwise gate**

```bash
git stash list  # ensure clean; commit locally first if the harness needs a ref
git add -A && git commit -m "wip" # temporary if needed, or run harness with workdir on dirty tree per its docs
python tests/regression/bc_cleanup/run_test.py <pre-task-2-sha> BCcleanup Debug --atol 0.0
```
Expected: PASS (prognostic dump fields identical). `thl0h/qt0h/ql0h(kb)` are not in the tdump output, so bitwise must hold exactly. If it fails, the Step 1 trace missed a consumer — STOP, report which field differs and where.

- [ ] **Step 5: Commit (squash the wip if used)**

```bash
git add src/modthermodynamics.f90 BCcleanup_backlog.md
git commit -m "fix(thermo): remove flat-surface thls/qts overwrite at kb (#299)"
```
Tick the two overwrite checkboxes in `BCcleanup_backlog.md` in this commit.

---

### Task 3: Derived base state — `modbasestate` (bitwise-on-matched-case commit)

**Files:**
- Create: `src/modbasestate.f90`
- Modify: `src/modstartup.f90` (three `initbasestate` call sites; delete `thvs=` line; bridge assignment; `checkinitvalues` validation)
- Modify: `src/modthermodynamics.f90` (`diagfld` seed + solid-slab fallback; `fromztop` anchor)
- Modify: `src/program.f90` (exit call)
- Modify: `CMakeLists.txt` only if sources are listed explicitly (check for a glob first)

**Interfaces:**
- Consumes: profile arrays `thlprof, qtprof` as read in `readinitfiles` (verify their declared bounds — assumed `(kb:ke)`; adapt the call if they are `1:kmax`).
- Produces: module `modbasestate` with `thl_b, qt_b, thv_b, pf_b, ph_b, exnf_b, exnh_b` (all `kb:ke+kh`) and `subroutine initbasestate(thlprof, qtprof)`, `subroutine exitbasestate`. Tasks 4-7 and Phases 1-2 rely on these names.

- [ ] **Step 1: Create `src/modbasestate.f90`**

```fortran
!> \file modbasestate.f90
!! Hydrostatic base state derived from the initial profiles and ps (issue #302).
!! Fixed per start; recomputed identically on cold, warm and strat starts.
module modbasestate
   implicit none
   save
   public :: initbasestate, exitbasestate
   real, allocatable :: thl_b(:)  !< base-state liquid water potential temperature [K]
   real, allocatable :: qt_b(:)   !< base-state total specific humidity [kg/kg]
   real, allocatable :: thv_b(:)  !< base-state virtual potential temperature [K]
   real, allocatable :: pf_b(:)   !< base-state hydrostatic pressure, full levels [Pa]
   real, allocatable :: ph_b(:)   !< base-state hydrostatic pressure, half levels [Pa]
   real, allocatable :: exnf_b(:) !< base-state Exner function, full levels
   real, allocatable :: exnh_b(:) !< base-state Exner function, half levels

contains

   subroutine initbasestate(thlprof, qtprof)
      use modglobal,   only : kb, ke, kh, zf, dzf, dzh, grav, cp, rd, rv, pref0
      use modsurfdata, only : ps
      use modmpi,      only : myid
      real, intent(in) :: thlprof(kb:ke), qtprof(kb:ke)
      real :: rdocp, thvh_b
      integer :: k

      rdocp = rd/cp

      if (.not. allocated(thl_b)) then
         allocate (thl_b(kb:ke+kh), qt_b(kb:ke+kh), thv_b(kb:ke+kh))
         allocate (pf_b(kb:ke+kh), ph_b(kb:ke+kh))
         allocate (exnf_b(kb:ke+kh), exnh_b(kb:ke+kh))
      end if

      thl_b(kb:ke) = thlprof
      qt_b(kb:ke)  = qtprof
      thl_b(ke+kh) = thlprof(ke)
      qt_b(ke+kh)  = qtprof(ke)
      thv_b = thl_b*(1.+(rv/rd - 1.)*qt_b) ! ql = 0 in the base state

      ! hydrostatic integration from ps at the domain bottom (z = 0), same
      ! discrete scheme as fromztop; defined over the full column, including
      ! levels inside terrain (reference-column continuation, backlog section 1.5)
      ph_b(kb) = ps
      pf_b(kb) = (ps**rdocp - grav*(pref0**rdocp)*zf(kb)/(cp*thv_b(kb)))**(1./rdocp)
      do k = kb + 1, ke + kh
         thvh_b  = (thv_b(k)*dzf(k-1) + thv_b(k-1)*dzf(k))/(2.*dzh(k))
         pf_b(k) = (pf_b(k-1)**rdocp - grav*(pref0**rdocp)*dzh(k)/(cp*thvh_b))**(1./rdocp)
         ph_b(k) = (ph_b(k-1)**rdocp - grav*(pref0**rdocp)*dzf(k-1)/(cp*thv_b(k-1)))**(1./rdocp)
      end do

      exnf_b = (pf_b/pref0)**rdocp
      exnh_b = (ph_b/pref0)**rdocp

      if (myid == 0) then
         write (*, '(A,F8.3,A,F10.1,A)') ' Base state: thv_b(kb) = ', thv_b(kb), &
            ' K, ps = ', ps, ' Pa (derived from prof.inp, #302)'
      end if
   end subroutine initbasestate

   subroutine exitbasestate
      if (allocated(thl_b)) deallocate (thl_b, qt_b, thv_b, pf_b, ph_b, exnf_b, exnh_b)
   end subroutine exitbasestate
end module modbasestate
```

Check `CMakeLists.txt`: if `src` sources are globbed, nothing to do; if listed, add `src/modbasestate.f90`.

- [ ] **Step 2: Call it from all three start paths in `readinitfiles`**

Add `use modbasestate, only : initbasestate` to `readinitfiles`. Insert
```fortran
      call initbasestate(thlprof(kb:ke), qtprof(kb:ke))
```
immediately after the profile MPI broadcasts in each branch:
1. cold start — after the `MPI_BCAST(e12prof...)` group that follows the `prof.inp` read at ~line 1108;
2. warm start — after the broadcasts at ~lines 1674-1678;
3. `lstratstart` — after its `prof.inp` read/broadcast (~line 1006 onwards).

The call must precede the first `call thermodynamics` inside the cold branch (~line 1607). Adapt the slice if `thlprof` is dimensioned `1:kmax` rather than `kb:ke`.

- [ ] **Step 3: Delete the `thvs` computation; bridge remaining readers**

At `src/modstartup.f90:522-525`, replace
```fortran
      thvs = thls*(1.+(rv/rd - 1.)*qts)
      call MPI_BCAST(thvs, 1, MY_REAL, 0, comm3d, mpierr)
```
with nothing here; instead, immediately after each `call initbasestate(...)` from Step 2, add:
```fortran
      thvs = thv_b(kb) ! bridge for remaining modsurfdata readers; removed in Phase 2
```
(add `thv_b` to the `use modbasestate` list). This keeps `lbottom`, `modsave`, `modfielddump` consumers defined — with a *correct* value (it also fixes the live qts-sentinel bug, backlog §1.5, in those paths).

- [ ] **Step 4: Re-point `modthermodynamics`**

In `diagfld`, replace (lines ~291-292):
```fortran
    exnf   = 1-grav*zf/(cp*thls)
    exnh  = 1-grav*zh/(cp*thls)
```
with:
```fortran
    ! base-state seed; refined below from the ps-anchored pressure
    exnf = exnf_b
    exnh = exnh_b
```
Add the solid-slab fallback immediately after the `ql0av` `avexy_ibm` call (~line 289):
```fortran
    ! slabs without fluid cells: avexy_ibm returns -999.; use the base state
    ! (reference-column continuation) so the hydrostatic march stays sane
    do k=kb,ke+kh
       if (IIcs(k) == 0) then
          thl0av(k) = thl_b(k)
          qt0av(k)  = qt_b(k)
          ql0av(k)  = 0.
       end if
    end do
```
In `fromztop`, replace `thvh(kb) = thvs` (line ~395) with `thvh(kb) = thv_b(kb)` and change its `use modsurfdata, only : ps, thvs` to `use modsurfdata, only : ps`. Add `use modbasestate, only : exnf_b, exnh_b, thl_b, qt_b, thv_b` where needed (`diagfld`, `fromztop`). Remove `thls` from `diagfld`'s `use modsurfdata` list (keep `ps`).

- [ ] **Step 5: Validation in `checkinitvalues` + exit call**

In `checkinitvalues` (`src/modstartup.f90`, after the existing runtime check ~line 767), add:
```fortran
      if ((ltempeq .or. lmoist) .and. ps < 0.) then
         if (myid == 0) then
            write (0, *) 'ERROR: ps must be set in &PHYSICS when ltempeq or lmoist is enabled.'
            write (0, *) 'The base state is derived from prof.inp and ps (issue #302). ps = ', ps
         end if
         stop 1
      end if
      if (lbottom .and. thls < 0.) then
         if (myid == 0) then
            write (0, *) 'ERROR: lbottom=.true. requires thls in &BC (legacy scheme, see #302).'
         end if
         stop 1
      end if
```
Verify the needed symbols are importable there (`ltempeq`, `lmoist` from `modglobal`; `lbottom` from `modibm`; `ps`, `thls` from `modsurfdata`) and extend the `use` lists accordingly. In `src/program.f90`, add `call exitbasestate` next to `call exitthermodynamics` in the shutdown sequence (with the module `use` added).

- [ ] **Step 6: Case audit**

```bash
grep -L "ps *=" tests/cases/*/namoptions.* examples/*/namoptions.*
grep -l "ltempeq *= *.true." tests/cases/*/namoptions.* examples/*/namoptions.*
```
Any case with `ltempeq`/`lmoist` but no `ps` now fails validation: add `ps = 101500.` to those namoptions (list them in the commit message). Do not touch anything else in the cases.

- [ ] **Step 7: Build + bitwise gate on 090**

```bash
./tools/build_executable.sh common debug
python tests/regression/bc_cleanup/run_test.py <pre-task-3-sha> BCcleanup Debug --atol 0.0
```
Expected: PASS bitwise. Reasoning (must hold, else investigate): case 090 is dry ⇒ the Exner seed is inert (`ql0av = 0` in the first-guess `th0av`); matched anchors ⇒ old `thvs = 290·(1+0.61·0) = 290 = thv_b(kb)`; no fully-solid slabs in the 090 geometry ⇒ fallback is a no-op. Startup log must show the new `Base state:` line.

- [ ] **Step 8: Commit**

```bash
git add src/modbasestate.f90 src/modstartup.f90 src/modthermodynamics.f90 src/program.f90 CMakeLists.txt tests/cases examples BCcleanup_backlog.md
git commit -m "feat(thermo): derive base state from prof.inp + ps; retire thls anchor (#302)"
```
Tick the derive/re-seed/fallback/validate checkboxes in `BCcleanup_backlog.md`.

---

### Task 4: SGS buoyancy re-point to `thvf(k)` (tolerance commit)

**Files:**
- Modify: `src/modsubgrid.f90` (lines ~388, ~488 and their `use` lists)
- Modify: `src/modforces.f90` (drop the dead `use modsurfdata, only : thvs`, line ~66-70; delete the `thvsi` fossil comments at 75-76)

**Interfaces:**
- Consumes: `thvf(kb:ke+kh)` from `modfields` (computed every step in `diagfld`).

- [ ] **Step 1: Re-point the two SGS terms**

`src/modsubgrid.f90:388`:
```fortran
                   zlt(i,j,k) = min(delta(i,k),cn*e120(i,j,k)/sqrt(grav/thvf(k)*abs(dthvdz(i,j,k))))
```
`src/modsubgrid.f90:488`:
```fortran
             sbbuo(i,j,k)  = -(ekh(i,j,k)-numol*prandtlmoli)*grav/thvf(k)*dthvdz(i,j,k)/ ( 2*e120(i,j,k))     ! subtract molecular diffusivity
```
Remove `thvs` from both routines' `use modsurfdata` lists; add `thvf` to their `use modfields` lists. Delete the stale trailing comments referring to `thls`/`thvs` on those lines. In `src/modforces.f90` delete the unused `use modsurfdata, only : thvs`, the `real thvsi` remnant comments at lines 75-76.

- [ ] **Step 2: Build + tolerance gate**

```bash
./tools/build_executable.sh common debug
python tests/regression/bc_cleanup/run_test.py <pre-task-4-sha> BCcleanup Debug --atol 5e-3
```
Expected: PASS at `5e-3` but FAIL at `0.0` (this is a deliberate physics change: constant 290 K → evolving `thvf(k)`; DALES ≥ 4.0 parity). Record the reported max diffs in the commit message. If diffs exceed 5e-3, inspect whether `thvf` is defined at all `k` used (it is set in `diagfld` lines ~336-345) before loosening anything.

- [ ] **Step 3: Commit**

```bash
git add src/modsubgrid.f90 src/modforces.f90 BCcleanup_backlog.md
git commit -m "fix(subgrid): use evolving thvf(k) in SGS buoyancy terms (DALES 4.0 parity)"
```
Tick the SGS checkbox in the backlog.

**Correction (2026-07-16):** the `--atol 5e-3` expectation above does not hold on case 090.
090 only sets `lvreman = .true.` (the default), so `closure` in `src/modsubgrid.f90` takes the
Vreman branch, never the `loneeqn` branch that contains the two re-pointed lines
(`zlt` at line 387, `sbbuo` at line 486). On 090 this task is dead code and the gate passes
bitwise (`--atol 0.0`), same as Tasks 2/3/5. `tests/cases/092/` (090 with
`lvreman = .false.` / `loneeqn = .true.`) was added to actually exercise the `loneeqn` branch;
the `5e-3` tolerance applies there. See `tests/regression/bc_cleanup/run_test.py`'s per-case
`default_atol` and the matching correction in `BCcleanup_backlog.md` §6.6.

---

### Task 5: Remove the `tl<100` clamp

**Files:**
- Modify: `src/modthermodynamics.f90` (`thermo`, lines ~484-494)

- [ ] **Step 1: Delete the clamp and the stale import**

In `thermo`, delete the X. Long comment block and the clamp:
```fortran
                !! X. Long: This is a fix to tackle incorrect thls input. ...
                if (tl<100.0) then 
                    tl=100.0
                end if
```
Delete `use modsurfdata, only : thls` from `thermo` (line ~442) — verify with `grep -n thls src/modthermodynamics.f90` that no live reference remains in the file.

- [ ] **Step 2: Build + bitwise gate**

```bash
./tools/build_executable.sh common debug
python tests/regression/bc_cleanup/run_test.py <pre-task-5-sha> BCcleanup Debug --atol 0.0
```
Expected: PASS bitwise — after Task 2, `tl` at `kb` derives from real air values and the clamp is dead code on any valid case. (Case 090 is dry so `thermo` only runs if `lmoist`; the gate still proves no accidental change.)

- [ ] **Step 3: Commit**

```bash
git add src/modthermodynamics.f90 BCcleanup_backlog.md
git commit -m "refactor(thermo): remove tl<100 clamp; root cause fixed (#299/#302)"
```

---

### Task 6: Register the suite + docs

**Files:**
- Modify: `tests/test_suites.yml` (experimental group)
- Modify: `docs/udales-namoptions-overview.md` (`thls`, `qts`, `ps` rows)

- [ ] **Step 1: Suite entry**

Add under the experimental group, following the schema of the entry at `tests/test_suites.yml:479-490`:
```yaml
      - label: "experimental: bc-cleanup phase-0 regression"
        class: experimental
        kind: regression
        component: solver
        platform: hpc
        cost: slow
        command:
          - "{python}"
          - "{tests_dir}/regression/bc_cleanup/run_test.py"
          - "{branch_a}"
          - "{branch_b}"
          - "{build_type}"
```

- [ ] **Step 2: Docs**

In `docs/udales-namoptions-overview.md`: mark `thls`/`qts` as *deprecated — used only by the legacy `lbottom` scheme; the thermodynamic base state is now derived from `prof.inp` and `ps` (removal planned with `lbottom`)*; on the `ps` row note it is *required when `ltempeq` or `lmoist`*.

- [ ] **Step 3: Run the harness once more (current HEAD vs Task 1 tag) and commit**

```bash
python tests/regression/bc_cleanup/run_test.py <task-1-sha> BCcleanup Debug --atol 5e-3
git add tests/test_suites.yml docs/udales-namoptions-overview.md
git commit -m "test/docs: register bc_cleanup suite; document derived base state"
```
Expected: PASS at 5e-3 (cumulative diff is Task 4's only).

---

### Task 7: Buried-slab fixture (backlog §6.6) — best effort

**Files:**
- Create: `tests/cases/091/` (copy of 090 with a full-floor solid slab)
- Modify: `tests/regression/bc_cleanup/run_test.py` (add case 091, assertion-only mode)

- [ ] **Step 1: Geometry**

Copy case 090 → 091 (rename suffixes, `iexpnr = 091`). Replace the geometry with a flat box covering the entire floor to a height of 2 grid cells (read `xsize/ysize/ktot/zsize` from `namoptions.091` for the exact extents), regenerate the IBM inputs with the Python preprocessing (`docs/udales-preprocessing-windows.md` — the repo venv at `tools/python` works on Windows; geometry how-to in `docs/udales-geometry-tutorial.md`). The result must make slabs `kb..kb+1` fully solid (`IIcs = 0` there).

- [ ] **Step 2: Assertions (single-branch mode)**

Add to the harness an `--assert-only <case>` mode that runs one branch only and checks: (a) the run completes; (b) the startup log contains `Base state:`; (c) no dumped field contains `-999.` or NaN anywhere, including the buried levels. Wire case 091 into it.

```bash
python tests/regression/bc_cleanup/run_test.py BCcleanup BCcleanup Debug --assert-only 091
```

- [ ] **Step 3: Commit — or defer with a written note**

If preprocessing cannot produce the fixture in reasonable time, STOP: commit nothing, add a dated note to `BCcleanup_backlog.md` §6.6 that the fixture is deferred to the Phase 1 PR, and report why.

```bash
git add tests/cases/091 tests/regression/bc_cleanup BCcleanup_backlog.md
git commit -m "test: buried-slab base-state fixture (backlog 6.6)"
```
