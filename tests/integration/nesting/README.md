# Nesting tests

Two layers live in this directory:

* **`test_nesting.py`** -- the in-solver unit runmodes 1006-1010 (U1-U34 of
  `docs/udales-nesting-design.md` section 10.1), which exercise one mechanism
  at a time through the public test hooks.
* **`test_nesting_cases.py`** -- the integration matrix I1-I8 of section 10.3,
  which runs the whole solver on nested cases built by
  `make_case_fixtures.py`.  See "Integration tests" below.

---

## Unit tests (runmodes 1006-1010)

In-solver unit tests for the one-way nesting feature, implementing the U1-U34
matrix of `docs/udales-nesting-design.md` section 10.1. They call the
production routines of `src/modnesting.f90` and `src/modnestingio.f90` through
the public test hooks of `docs/udales-nesting-spec.md` section 7 -- nothing
under test is reimplemented here.

Each runmode is a `logical function` in `src/tests.f90`, dispatched from
`execute_runmode_actions` in `src/program.f90`, and exits 0 (all subtests
passed) or 1.

| Runmode | Function | Covers |
|---|---|---|
| 1006 | `tests_nesting_weights` | U1-U7: shape function value/range/monotonicity, C1 and C2 joints with zero slope at both ends, bounded-union bounds/symmetry/degeneracy, the integral identity `int W ds = L_imp + L_rel/2` |
| 1007 | `tests_nesting_geometry` | U8-U14: stagger coordinates against an independent global formula, zone membership and `sum W` against a serial global reference, IBM masking, wall erosion, the building-free rule, width reporting |
| 1008 | `tests_nesting_io` | U15-U22: spatial read against the analytic field, per-rank hyperslab equals the full read bitwise, time-interpolation exactness, Hermite C1, buffer roll vs full reload, restart repositioning, header validation |
| 1009 | `tests_nesting_flux` | U23-U28: `nest_flux_residual` against an analytic net flux, decomposition invariance, fluid-face masking, a corrected file giving `Phi = 0`, the assertion firing on an uncorrected one, linearity of `Phi` in time |
| 1010 | `tests_nesting_update` | U29-U34: bitwise no-op off zone, Dirichlet limit, linear limit, full-step composition `exp(-W dt/tau)`, stability over `dt/tau` in `[1e-3, 1e6]`, solid points untouched |

## Running

```bash
./tools/build_executable.sh icl release
python tests/integration/nesting/test_nesting.py
```

The Python driver builds the fixtures, runs every runmode on 1x1, 2x1, 1x2 and
2x2 ranks and asserts the exit code, then runs the abort cases and the
message cases. To run one runmode by hand:

```bash
python tests/integration/nesting/make_fixtures.py <rundir>
cp tests/integration/nesting/{namoptions.*,prof.inp.901} <rundir>
cd <rundir> && mpiexec -n 4 .../build/release/u-dales namoptions.1008
```

(set `nprocx`/`nprocy` in the namoptions to match `-n`).

Also run the tests against a **Debug** build: `-check bounds` catches things a
release build reads straight past (see "Known failures" below).

## Fixtures

`make_fixtures.py` writes everything through the production writer in
`tools/python/udprep/nesting.py`, so the writer and the Fortran reader cannot
drift apart. Nothing is committed; the files are regenerated on every run.

| File | Content |
|---|---|
| `nesting_analytic.901.nc` | exactly `analytic_field()` with `ANALYTIC_COEFFS`, **not** divergence corrected, so every stored value equals the analytic function to the last bit. Used by U15-U21 and U23-U34, with `nest_fluxtol` raised so the init-time flux check does not reject it. |
| `nesting_corrected.901.nc` | the same field with the offline divergence correction (U26). |
| `assertfire.901.nc` | a copy of the uncorrected file, run at the production tolerance so the flux assertion must fire (U27). |
| `nesting_nonlinear.901.nc` | analytic in space, **non**-linear in time: stored level `n` holds the field at pseudo-time `s_n = t_n^2/t_max`. U18 and U19 need this -- a field linear in `t` is reproduced exactly by every interpolant and so cannot separate linear from Hermite. |
| `bad_schema/itot/xlen/zf/stagger.901.nc` | one corrupted header item each (U22). |

The case is 32 x 32 x 16 cells over 32 x 32 x 16 m (`dx = dy = dz = 1 m`),
`nzone = 12`, six parent levels at `t = 0, 10, ..., 50 s`, and
`L_imp = 3 m`, `L_rel = 8 m`. It carries no IBM input files: the solid points
the IBM subtests need are written straight into `IIu/IIv/IIw` from a
global-index rule, so the same physical cells are marked on every
decomposition.

## How the private state is reached

The zone weights and the time-interpolated target live in private variables of
`modnesting` and there is no accessor. Rather than recompute them, the tests
read them back out through `nesting_apply`, which is public:

* **weight probe** -- `nesting_apply` leaves
  `up = (target - um)(1 - exp(-W dt_s/tau))/dt_s`. Running it twice from
  `um = 0` and `um = 1` cancels the unknown target and leaves
  `1 - exp(-W dt_s/tau)`; with `rk3step = 0` and `tau = 1` that inverts to `W`
  exactly. Points outside the zone are untouched and come back as `W = 0`.
* **target probe** -- `tau <= 0` makes the update Dirichlet, so with `um = 0`
  and `dt_s = 1`, `up` *is* the time-interpolated parent target.

Both probes are exact, use only production code, and are independent of each
other's subject.

## Abort cases

U13, U22 and U27 assert that a routine *stops*, so each runs as its own
process. They are selected from the namelist, and the driver checks both the
non-zero exit and the message:

| Case | Namelist selector | Expected message |
|---|---|---|
| U13 | `nest_lparentgeom = .false.` (runmode 1007) | `solid points found inside the relaxation zone` |
| U22 | `nestfile = 'bad_*.901.nc'` (runmode 1008) | `mismatch in <field>` |
| U27 | `nestfile = 'assertfire.901.nc'` (runmode 1009) | `not flux balanced` |

## Known failures and gaps

**U12 / `build_eroded_mask` -- FIXED, this note is kept for the history.**
`build_eroded_mask` used to read `prev(i+di, j+dj, k+dk)` one plane past the
upper bound at `k = ke+kh`, which aborted a Debug build and silently wiped the
`k = ke` zone layer in a Release build once `nest_nwall >= 2`. It now clamps
`k+dk` into `kb-kh:ke+kh` (`src/modnesting.f90:1219-1225`). Runmodes 1006-1010
pass against **both** the Release and the Debug build, on 1x1, 2x1, 1x2 and
2x2, as of this writing.

**U22 stagger tag.** Design section 10.1 asks for a mismatched `stagger`
attribute to abort. The normative reader contract
(`docs/udales-nesting-spec.md` section 6) lists the header items
`nestio_validate` compares and the per-variable `stagger` attribute is not
among them, and `modnestingio` never reads it. The runmode reports the
discrepancy as an `INFO` line rather than asserting a requirement the contract
does not make; the driver checks that the line is still printed. Either the
spec or the reader needs a decision here.

**U21 prefetch.** `modnesting` has no prefetch on/off switch; `set_interval`
has an incremental roll-with-read-ahead branch and a full-reload branch. The
test compares the two branches against each other and against a fresh
initialisation, which is the invariance the design asks for.

**U9/U10, U24 decomposition independence.** These are checked against a serial
reference computed from global indices only, so agreeing with it on any layout
means the layouts agree with each other. The driver additionally runs every
runmode on all four layouts.

---

# Integration tests (I1-I8)

`test_nesting_cases.py` implements the I1-I8 matrix of
`docs/udales-nesting-design.md` section 10.3. Where the runmodes above call one
routine at a time, these run the **whole solver** on a nested case and check
the properties the composed scheme is supposed to have.

```bash
./tools/build_executable.sh icl release      # and 'debug' -- run both
module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
source /rds/general/user/mvr/home/udales/.venv/bin/activate
UDALES_BUILD=$PWD/build/release/u-dales \
  python tests/integration/nesting/test_nesting_cases.py
```

Useful environment variables: `UDALES_BUILD`, `UDALES_BASELINE` (the pre-branch
binary I1 compares against, default `build/u-dales.baseline` -- **not in
version control**; `build/` is gitignored, so build it yourself with
`git stash && ./tools/build_executable.sh icl release &&
cp build/release/u-dales build/u-dales.baseline && git stash pop`, or point
`UDALES_BASELINE` at one. I1 skips itself if it is missing),
`UDALES_RUNTIME_MODULES`, `MPIEXEC`, and `UDALES_NESTING_KEEP=1` to leave the
run directories behind for inspection.

| Test | Class | Verdict |
|---|---|---|
| I1 | `TestI1NoOpSmallCase`, `TestI1NoOpExistingCase` | passes |
| I2 | `TestI2FaceImposition` | passes |
| I3 | `TestI3UniformFlow` | passes |
| I4 | `TestI4ManufacturedSolenoidal` | passes |
| I5 | `TestI5DecompositionParity` | passes |
| I6 | `TestI6RestartParity` | passes (after the `program.f90` fix below) |
| I7 | `TestI7ZoneIsolation` | passes |
| I8 | `TestI8IbmInteraction` | passes |

## Cases and fixtures

`make_case_fixtures.py` builds a complete run directory -- `nesting.inp`,
`prof.inp` (which also *defines* the vertical grid), `lscale.inp` and a
`namoptions` filled in from the `namoptions.902` template in this directory.
Every nesting file goes through the production writer
`tools/python/udprep/nesting.py`, as for the unit fixtures.

The parent is defined **on the child grid** (refinement ratio 1). That is
deliberate: it removes the interpolation from the loop, so an error seen here
is the solver's. The interpolation is covered separately by P1-P11 in
`tools/python/tests/test_nesting.py`.

Three geometries:

| Name | Grid | Zone | Used by |
|---|---|---|---|
| `BASE` | 32 x 32 x 16 over 32 x 32 x 16 m | `L_imp` = 16 m, `L_rel` = 0, `tau` = 0: `W == 1` **everywhere** | I2, I3, I5 |
| `CONVERGENCE` | the same box at 16/32/64 cells per side | as `BASE` | I4 |
| `ZONED` | 32 x 32 x 16 | `L_imp` = 3 m, `L_rel` = 8 m, `tau` = 4 s | I1 (small case), I6, I7 |

`W == 1` everywhere is obtained by making the guard strip half the domain, so
every point is within `L_imp` of a lateral face and the bounded union of design
section 1.1 returns 1. `nzone` is then half the horizontal size, the smallest
value for which the west/east (south/north) slabs together cover every child
index. Such a case trips the "zone occupies 100 % of the domain" warning from
`nesting_init`; that is expected.

Four manufactured fields, all with `Phi = 0` exactly:

| Field | Definition | Property |
|---|---|---|
| `uniform` | `u = U`, `v = w = 0` | a fixed point of every term on a free-slip box |
| `face_forced` | `u = U0 + A cos(2 pi x/Lx) cos(2 pi y/Ly) sin(2 pi z/Lz)` | identical on both `x` faces, so `Phi = 0`, but `du/dx != 0` so the projection has real work |
| `taylor_green` | `u = sin(2 pi x/L) cos(2 pi y/L)`, `v = -cos(2 pi x/L) sin(2 pi y/L)` | design section 10.3's field for I4 |
| `mixed_mode` | `u = dpsi/dy`, `v = -dpsi/dx` for `psi = sin(2 pi x/L) sin(4 pi y/L)` | solenoidal analytically, discrete divergence genuinely `O(h^2)` |

`python tests/integration/nesting/make_case_fixtures.py` prints each field's
peak discrete divergence and net boundary flux on all three grids, without
running the solver.

## Reading the solver's fields

Every NetCDF output in uDALES is written as `NF90_FLOAT`
(`src/modstat_nc.f90`), which cannot resolve a 1e-9 parity claim let alone a
round-off one. All field assertions therefore read the `initd` **restart**
file, which is unformatted double precision; `read_restart_fields` stitches the
per-rank files into global arrays laid out like the solver's own
`(kb:ke+kh, jb-1:je+1, ib-1:ie+1)` block, so the domain ghost planes -- which
is where the east face `u(ie+1)` and the north face `v(je+1)` live -- come
along. Global reduced diagnostics (`Phi`, the zone misfit, `|grad p|` zone and
interior, the guard/relaxation energy split, `divmax`/`divtot`) are parsed from
stdout.

## What each test establishes

**I1 -- no-op guarantee.** Split in two, because "bitwise identical to the
pre-branch binary" is achievable on a small case and not on a large one, for a
reason that has nothing to do with this branch: uDALES plans all of its FFTs
with `FFTW_MEASURE` (`src/modpois.f90:110-191`,
`2decomp-fft/src/fft_fftw3.f90:26`), which chooses the transform algorithm by
timing it at run time. Two runs of *the same* binary on `tests/cases/526` differ
at the 1e-16..1e-15 relative level, with a fixed timestep and on one rank, with
or without this branch.

* `TestI1NoOpSmallCase` runs a 32 x 32 x 16 periodic case with
  `lnesting = .false.`, first verifying that the baseline reproduces *itself*
  bitwise, and then requiring the branch to match it bitwise. It does: every
  restart record is byte-identical and stdout is line-identical.
* `TestI1NoOpExistingCase` runs `tests/cases/526` unchanged (IBM, trees,
  statistics, 4 ranks, `ladaptive` forced off so a round-off difference cannot
  feed back into a different timestep). It measures the baseline's own
  run-to-run spread and requires the branch to be no further away than that.
  Measured: baseline vs baseline 7.3e-12 worst relative field difference
  (dominated by `pres0`, whose values sit near zero); baseline vs branch
  6.4e-12, i.e. *smaller*.

`divmax`/`divtot` are compared with a tolerance rather than bitwise on case 526,
and the test prints the baseline-vs-baseline spread that justifies it. On the
small case, which is bitwise reproducible, they are required to match exactly.

**I2 -- face imposition survives projection (design finding F2).** Measured:
`max |u - u_parent|` on all four imposed faces is 2.2e-16 (one ulp) after
`poisson` + `tstep_integrate`, on a parent whose peak interior divergence is
3.8e-2 and whose `|grad p|` is far above round-off. A companion test asserts
both of those, so the result cannot pass vacuously.

**I3 -- uniform flow.** Measured: `max|u - U| = 0`, `max|v| = max|w| = 0`,
`max|pres0| = 0`, `|Phi| <= 1.4e-17`, zone misfit 1.1e-17, `|grad p|_zone = 0`,
`divmax = divtot = 0`. Exactly preserved, not approximately.

**I4 -- manufactured solenoidal field.** Two halves, because the field design
section 10.3 names turns out to be a *stronger* case than the design expected.

On a staggered grid the Taylor-Green field's discrete divergence cancels
identically rather than to `O(h^2)`: for `u` at `(xh_i, yf_j)` and `v` at
`(xf_i, yh_j)`,

```
du/dx = 2 sin(pi h/L)/h * cos(2 pi xf_i/L) cos(2 pi yf_j/L)
dv/dy = -2 sin(pi h/L)/h * cos(2 pi xf_i/L) cos(2 pi yf_j/L)
```

which sum to zero for any `h` when `dx = dy`. Measured peak `|div|` is 4.2e-16,
8.3e-16 and 1.8e-15 on `h` = 2, 1 and 0.5 m -- round-off, and *growing* with
resolution, which is the signature of accumulation and not of truncation. So
`||grad p||` sits at round-off (3.1e-16, 3.1e-16, 3.0e-16) on all three grids
and there is no convergence rate to measure; asserting second order would be
asserting a property of the round-off. `test_taylor_green_needs_no_correction`
therefore requires `|grad p|`, `divmax` and `Phi` to be at round-off on every
grid, and `test_taylor_green_divergence_is_exact` records why.

The convergence half runs `mixed_mode` instead, built the same way from a
stream function but with different wavenumbers in `x` and `y`, so the
cancellation above does not occur. Measured discrete divergence 1.33e-3,
3.62e-4, 9.23e-5 (orders 1.88, 1.97) and `||grad p||` 4.04e-3, 9.91e-4,
2.46e-4 (**orders 2.03, 2.01**). The timestep is held fixed across the three
grids, which it must be: `grad p ~ div * L / dt_s`.

**I5 -- decomposition parity.** `uniform` and `taylor_green` on 1x1, 2x1, 1x2
and 2x2. Measured: `u0`, `v0`, `w0` and `pres0` agree **exactly** (0.0, not
1e-9) on every layout. Since both fields are imposed every substep, what this
really exercises is the slab hyperslab read, the zone point lists and the
`facval` indexing, all of which do change with decomposition.

**I6 -- restart parity.** This test found the one real defect of the M5 pass,
in `src/program.f90`; the fix has landed and I6 is now fully bitwise on both
restart points. See "One defect, two symptoms" below for the write-up, kept
because the failure mode is worth not reintroducing.

**I7 -- zone isolation.** Constructed so that the two runs differ *only* in the
interior: a spin-up writes a restart file, the file is copied, and in the copy a
solenoidal perturbation built from `psi = A Wx(x) Wy(y)` is added, where `Wx`
and `Wy` are raised-cosine windows that are identically zero within
`L_imp + L_rel` of a lateral face. The test asserts that the perturbation is
exactly zero everywhere in the zone before either run starts. Both files are
then warm-started for 1 s, over which the mean flow advances one cell -- so
anything that reaches the zone got there through the elliptic pressure solve,
which is concern C1 of design section 7.

Measured, with a 3 m guard and an 8 m ramp:

```
d =  0..2 cells (guard):   6.4e-06  7.0e-06  8.5e-06
d =  3..10     (ramp):     1.5e-05 ... 1.28e-02
d = 11..15     (interior): 3.6e-02  5.7e-02  2.8e-02  2.6e-03  1.5e-03
```

The interior difference of 5.7e-2 m/s is attenuated to 8.5e-6 in the guard
strip: a factor of 6.7e3, e-folding over roughly 1.3 cells inside the ramp. On
the imposed boundary faces themselves the difference is **exactly zero** for
all four faces -- the pressure Neumann condition of `bcp` leaves the
boundary-normal velocity untouched, so C1's global pressure response does not
move the imposed values at all. Reported alongside: `|grad p|` zone/interior
ratio 1.06-1.14 and the guard/relaxation energy split.

**I8 -- IBM interaction.** `tests/cases/064` is a single 6 m cube in a
64 x 64 x 64 m box with its windward face at `x = 24 m`, so a `3 + 20 = 23 m`
zone puts the inner edge exactly one cell upstream of the building --
"buildings adjacent to the zone edge", and the tightest arrangement the
building-free rule of design section 5 allows. `nest_lparentgeom` is left at
`.false.`, so `nesting_init` aborts if any of the geometry falls inside the
zone; it does not. The reference is the same case and the same initial
condition, driven periodically at the same volume flow rate. Energy balance,
temperature and moisture are off so the comparison is purely mechanical.

Measured on the two windward facets (`fac.064.nc`, `lwritefac = .true.`):
`tau_y` 3.1e-2, `tau_z` 3.4e-2, `pres` 4.9e-2 relative difference, against a
stated tolerance of 10 %; the leeward facets differ by 1.0e-2, 8.6e-3 and
5.1e-3. This is a **modelling** comparison, not an exactness one -- the two
runs have genuinely different boundary conditions -- so the number is printed
as well as asserted. The nested run's `divmax` is 1.9e-15, `Phi` is 0, and the
`|grad p|` zone/interior ratio is 0.16: with buildings, the projection works
*less* hard in the zone than in the interior.

## One defect, two symptoms: `nesting_init` ran before `timee` existed

**Fixed** -- `src/program.f90` now calls `nesting_init` after `readinitfiles`.
Recorded here because the failure mode is easy to reintroduce and neither
symptom points at the cause.

As found: `TestI6RestartParity` failed, and *every* nested case crashed
outright against a Debug build. Both were the same single line, and the test
suite was deliberately not weakened around either.

### The defect

`program.f90:79` dispatches `execute_runmode_actions`, so the unit runmodes
1006-1010 return before `nesting_init` is ever called -- which is why the
34 unit tests pass against both builds and could never have caught this.

`nesting_init` was called from `program.f90:103`. `readinitfiles`, which is
what assigns `timee` (`src/modstartup.f90:1203` on a cold start, and
`readrestartfiles` on a warm one), was called from `program.f90:105` -- *after*
it. So `nesting_init` ran

```fortran
      it_lo = 0
      ...
         call set_interval(timee)
      end if
      call reload_all
      ...
         call eval_target(timee)
```

(`src/modnesting.f90:238-250`) on a `timee` that has not been assigned yet.
`real :: timee` in `src/modglobal.f90:441` carries no initialiser.

### Symptom 1 -- every nested run dies on a Debug build

The Intel Debug flags include `-init=snan` and `-fpe0`
(`CMakeLists.txt:62-66`), so the undefined `timee` is a signalling NaN and the
first arithmetic that touches it traps:

```
forrtl: error (75): floating point exception
  modnesting_mp_eval_target_   1627  modnesting.f90     ! th = (t - tlo)/h2
  modnesting_mp_nesting_init_   249  modnesting.f90
  MAIN__                        103  program.f90
```

The trap location follows the *first use of `t`*, which pins the operand: with
two stored parent levels the `do while` loops in `set_interval` do not execute
and the trap lands in `eval_target` at line 1627; with six levels it moves up
to `set_interval` line 1580, `if (t < nestio_hdr%time(il + 1)) exit`. `tlo`,
`thi` and `h2` are all computed before line 1627 without trapping, so the NaN
is `t` and nothing else.

In a Release build `timee` reads as 0 from static storage, which is why cold
starts appear to work.

### Symptom 2 -- nested runs do not restart bitwise

On a warm start `timee` is still 0 at `nesting_init`, so the four-slot parent
buffer and the interpolated target are positioned at `t = 0` instead of at the
restart time. `program.f90:124` then calls `boundary`, whose nesting branch
(`src/modboundary.f90:384` -> `nesting_boundary`) fills the lateral ghost planes
**and** `u0`/`um` on the boundary faces from that `t = 0` target; the first
substep after the restart advects against them.
`nesting_update_target` only repairs the buffer at `program.f90:144`, one
`boundary` call too late.

Measured, 100 steps versus 50 + restart + 50:

```
mid-interval  u0 max abs diff 8.7e-08   v0 4.7e-08   w0 4.6e-08   pres0 1.7e-09
on-boundary   u0 max abs diff 1.7e-07   v0 9.9e-08   w0 1.2e-07   pres0 4.8e-09
```

Three tests localise it rather than merely reporting it:

* `test_control_restart_without_nesting` runs the identical case with periodic
  laterals and `lnesting = .false.` and is **bitwise** on every restart record.
  uDALES's restart machinery is sound; the difference is nesting's.
* `test_restart_with_a_time_constant_parent` replaces the unsteady parent with
  a constant one and **passes**: with the target the same at `t = 0` as at the
  restart time, mispositioning the buffer cannot matter and the restart comes
  back at round-off.
* `test_nesting_init_positions_the_buffer_at_the_restart_time` reads the defect
  straight out of the solver's own output. Restarting at `t = 4.625 s` from a
  file whose parent levels are `0, 1.25, 2.5, 3.75, 5.0, ...`, `nesting_init`
  prints `buffer at interval 1` where it should print `4`.

### The fix (applied)

Run `nesting_init` after the thing that defines what it reads:

```diff
   call calcfluidvolumes

-  call nesting_init
-
   call readinitfiles

+  ! nesting_init reads `timee`, and readinitfiles is what assigns it -- 0 on a
+  ! cold start, the restart time on a warm one. Initialising nesting before it
+  ! positioned the parent buffer at an undefined time (a signalling NaN under
+  ! the Debug -init=snan) and, on a warm start, at t = 0 instead of the restart
+  ! time (docs/udales-nesting-design.md section 9.5).
+  call nesting_init
+
   call createscals
```

Nothing between the two calls needs nesting: `nesting_boundary`,
`nesting_apply`, `nesting_bcpup` and `nesting_stats` all return immediately
while `linit` is false, and `readinitfiles` does not call `boundary`.
`nesting_init` needs `createmasks` and `calcfluidvolumes`, which still precede
it.

With the fix in and both builds rebuilt:

* the Debug build's floating-point trap is gone -- the whole matrix runs there
  (22 tests, 2 I1 skips), with values matching the Release build;
* `TestI6RestartParity` is fully green -- all ten restart records bitwise
  identical for both the mid-interval and the on-a-parent-level restart, and
  `nesting_init reports interval 4, expected 4`.

**Design section 9.5 is deliberately not wired in.** `modnesting` provides
`nesting_restart_write` and `nesting_restart_read`; `grep -rn nesting_restart
src/` finds call sites only in `src/tests.f90`. They are unit-tested (U20) and
unused by the solver, and the doc comments now say so. Once `nesting_init` runs
after `readinitfiles`, `set_interval(timee)` reconstructs the buffer state
*exactly* -- I6 proves it bitwise -- so storing it would pin state that is
already recoverable, at the cost of changing the `initd` record layout.


## Running against a Debug build

Every test except I1 is self-contained and should be run against both builds:

```bash
UDALES_BUILD=$PWD/build/debug/u-dales \
  python tests/integration/nesting/test_nesting_cases.py
```

I1 skips itself when `UDALES_BUILD` looks like a Debug build. The committed
baseline `build/u-dales.baseline` is a Release binary (see above -- you build
it yourself, it is not in the repository), so comparing a Debug
build against it would measure the optimisation level rather than the branch.

Measured: **22 tests, 22 pass, 2 skips** (the two I1-vs-baseline tests), in
1030 s. Before the `program.f90` fix in the section above, every nested case
(I2-I8) trapped here on `-init=snan`, which is precisely what makes the Debug
run worth keeping in the loop.

## Current status

Against `build/release/u-dales` on this branch: **28 tests, 28 pass** (561 s).
Against `build/debug/u-dales`: **22 tests, 22 pass, 2 skips** (1030 s) -- the
skips are the two I1-vs-baseline comparisons, which need a Release binary.
The unit runmodes (1006-1010) pass against both builds.
