# Test design for #332 (retire flat-surface scheme)

Each test below is tied to a specific claim the PR makes, and is written so that
it fails if that claim is false. The ordering is by value, not by effort.

The guiding rule: a regression comparison shows a number has not *changed*; it
cannot show the number was ever *right*. Everything the PR asserts as a
deliberate improvement needs an independent expectation, not a self-comparison.

Status: (done) landed and verified · (built) written, verification pending ·
(todo) designed here only.

---

## 1. Base state arithmetic — (built) `runmode 1006`, `tests/integration/basestate`

**Claim.** The base state (Exner, hydrostatic pressure, `thv` reference) is
correctly derived from `prof.inp` + `ps`, replacing the `thls` sentinel.

**Why the existing tests do not cover it.** 090/092 compare a build against
another build of the same branch; 091 asserts only that the string `Base state:`
was printed. A base state that was wrong from the first commit passes all three.
That is precisely the shape of the bug #302 fixes (`thvs ≈ 0.39·thls` → ~116 K).

**Design.** Dispatch at the runmode point (after `initglobal`, before
`readinitfiles`), supply profiles in-test, call `initbasestate` directly. For a
uniform-`thv` column the discrete integration telescopes, so there is a closed
form to compare against:

    p(z)**(Rd/cp) = ps**(Rd/cp) - g*pref0**(Rd/cp)*z/(cp*thv)

at `zf` for full levels, `zh` for half levels.

**Assertions.** `thv_b(kb) == thl(kb)` (the #302 signature); `ph_b(kb) == ps`
bitwise; `pf_b`/`ph_b` vs closed form (rtol 1e-10); `exnf_b == (pf_b/pref0)**rdocp`;
pressure strictly decreasing; moist `thv = thl*(1+(rv/rd-1)*qt)`.

**Falsifiable.** Verified by mutation: `grav -> 1.001*grav` in the integration
must make it fail. A test that has only ever been green proves nothing about its
own sensitivity.

---

## 2. Retired-key contract — (todo) `tests/integration/retired_keys`

**Claim.** "Old namoptions setting any retired key abort at the namelist read."
This is the PR's main user-facing promise; it is what stops six of the eight
cases under `experiments/` at startup.

**Why it matters.** Nothing tests it. If a key were silently accepted and
ignored, a user's case would run with a setting that no longer does anything —
the exact failure mode the loud-abort design exists to prevent, and it would look
like success.

**Design.** One minimal namelist per retired key, each a copy of a valid one plus
that single key. Run each; require a non-zero exit *and* the key's name in the
output. No MPI, no geometry, no timesteps — it fails at the namelist read, so the
whole set costs seconds.

**Keys.** `thls`, `qts`, `z0`, `z0h`, `wtsurf`, `wqsurf`, `wsvsurf`, `BCbotm`,
`BCbotT`, `BCbotq`, `BCbots`, `lbottom`, and the inlet-generator set (`di`, `dti`,
`linletRA`, `lfixinlet`, `lfixutauin`, `lreadminl`, `lstoreplane`, `lwallfunc`).

**Assertion detail.** Require the *message* to name the offending key, not merely
a non-zero exit: a case that aborts for an unrelated reason would otherwise pass
and the test would rot silently.

---

## 3. Buried-slab fallback values — (todo) extend `tests/cases/091`

**Claim (spec §6.6).** Where a slab is fully solid (`IIcs(k) == 0`), the
thermo-facing slab means fall back to the base profiles instead of `avexy_ibm`'s
`-999.`, and `presf` at the fluid levels equals the reference-column continuation
with the buried layers' weight included.

**Why the existing 091 assert is not enough.** It checks only that no dumped
field is `-999`/NaN. "Not `-999`" is a much weaker claim than "equals `thl_b`": a
fallback writing zeros, or the wrong profile, or the right profile shifted by one
level, passes today. The spec asks for three assertions; only the weakest exists.
(The backlog describes the `presf` check as "implemented but unexercised" — it
was never written.)

**Design.** Two parts.

*Values.* In the assert-only path, for every `k` where the slab is fully solid
(known from `solid_c.txt`), require the dumped `thl`/`qt` to equal the `prof.inp`
profile at that `k` to roundoff, not merely differ from `-999`.

*Continuation.* This is the check that measures the PR's own justification — the
backlog argues the poisoned march biases `presf` at **all levels above** by
~0.12 hPa per metre of buried layer. Reconstruct the reference-column
continuation from `prof.inp` + `ps` (same closed form as test 1, since 090/091's
profile is uniform-`thv`) and compare `presf` at the fluid levels. If `presf` is
not dumped, add it to the assert-only path via a diagnostic write rather than
inferring it.

**Falsifiable.** Mutation: force the fallback to write `thl_b(k)+0.1` and confirm
failure; force the march to skip a buried layer's weight and confirm the
continuation check fails.

---

## 4. Driver inflow lifecycle — (done) `tests/integration/driver_inflow`

**Claim.** `moddriver` + `modinletdata` merged into `inflow.f90` is a rename-only
relocation, and `exitinflow` is now wired into teardown.

**Design.** Cases 501 (`idriver = 1`) → 502 (`idriver = 2`), 8³ on 2 ranks. The
precursor writes planes; the driver consumes them; both must complete. Covers
`initinflow`, `drivergen`, `writedriverfile`, `readdriverfile`, and both branches
of `exitinflow`.

**Assertion detail.** Requires the "not given by precursor" warning to be absent.
Without that, the solver would silently fall back to a non-precursor inflow and
the test would pass while covering half of what it claims.

**Note.** `ltempeq`/`lmoist` are deliberately on, unlike the neutral
`experiments/525`: they are what make `BCxT = 3`/`BCxq = 3` engage
`lhdriver`/`lqdriver`, so the thl/qt planes are exercised too.

---

## 5. Retired dump fields — (todo) fold into `retired_keys`

**Claim (spec §6.5).** The per-cell `tau_x`/`tau_y`/`tau_z`/`thl_flux` dump
fields are removed, "not silently zeroed", and an unknown `fieldvars` code is now
a loud startup error rather than a silent mislabelled `u0` dump.

**Design.** `fieldvars = 'tau_x'` must abort at startup naming the field. Same
harness as test 2. This covers the removal and the new `fieldvars` validation in
one cheap check — and the old behaviour (silently dumping `u0` under another
name) is exactly the kind of failure that no output comparison would catch,
because the file would look perfectly well-formed.

---

## 6. `ps` validation — (todo) fold into `retired_keys`

**Claim.** "`ps` is validated at startup when thermodynamics is active."

**Design.** `ltempeq = .true.` with no `ps` must abort. One namelist, same
harness. Guards the replacement for the sentinel the PR removed: having deleted
`thls`, an unset `ps` must not become the new silent-garbage input.

---

## 7. Migrated cases 103 / 999 — (todo) needs a decision first

**Claim.** Both `lbottom` cases migrated to ground facets: 103 to a fixed 288 K
ground facet (the equivalent of its old `BCbotT=2`/`thls=288`), 999 to a neutral
2-facet floor. The PR asks for statistical-equivalence runs because the
wall-function formulation legitimately changed.

**Problem.** Neither case has *any* automated coverage: 103 lives under
`david_tests`, whose `run_and_compare` is commented out *and* points at
directories that do not exist; 999 is referenced by no suite at all. With
`david_tests` being removed, 103 needs a new home regardless.

**Design (pragmatic).** Statistical equivalence against the old scheme needs a
pre-PR build and a criterion nobody has written down, which is a poor merge gate.
Cheaper and honest:

- move 103 into `tests/cases/`, register a fast smoke suite: runs clean, no NaN,
  and the ground facet holds 288 K (its migration's whole point);
- same for 999 (neutral floor, 2 facets).

That converts "no coverage" into "the migration's stated intent is asserted",
without pretending a statistical claim has been checked. If statistical
equivalence is genuinely wanted, it should be a one-off study attached to the PR,
not a permanent suite member.

---

## Not worth building

**Bitwise comparison against `f5de904c`.** The PR body prescribes it. It cannot
work, for three independent reasons, the third fatal: (a) the reference silently
runs with `thls = -1` because Task 3 pruned the key the old code requires;
(b) `lfftwmeasure` is unknown to the old namelist reader; (c) **a nondeterministic
executable cannot produce a bitwise baseline** — the old code has no planner
switch, so its own output varies run to run. Any historical comparison must be
tolerance-based above the ~1e-7 noise floor, which demonstrates "no large change",
not "no behavioural change".

**A `loneeqn` fixture as a merge gate.** Case 092 is the repo's only `loneeqn`
case and every case under `experiments/` uses Vreman. Task 4 (`thvs` → evolving
`thvf(k)`) only executes in the `loneeqn` branch, so the PR's one
behaviour-changing commit does not touch production usage. Keep 092 as a
correctness check on a path that is real but rarely used; do not treat it as
load-bearing.

---

## Infrastructure notes for whoever runs these

- **Reproducibility is a precondition, not a detail.** `FFTW_MEASURE` picks FFT
  plans by runtime benchmarking, so two identical runs differ in ~15/32 dumped
  fields (~1 ULP of the single-precision dumps). Any bitwise test must set
  `lfftwmeasure = .false.`; see `tests/regression/bc_cleanup/README.md`.
- **Every test must state why it failed.** Two runs in this work were killed by a
  compute node where Hydra could not connect to itself; the symptom was
  indistinguishable from "the code broke". Grep the failure text before believing
  a verdict, and prefer `I_MPI_HYDRA_BOOTSTRAP=fork` for single-node jobs.
- **The stale `intel/2025a` module stack** is duplicated in at least four
  harnesses (`bc_cleanup` fixed; `mpi_averaging_regression`,
  `new_vegetation_module_against_v2.2`, `integration/mpi_operators` not). It
  mixes three toolchain years and does not co-load on CX3.
