# V1 -- the "Big Brother" nesting validation

This directory implements test **V1** of `docs/udales-nesting-design.md`
section 10.4:

> Does matched LES-to-LES nesting reproduce the parent?

It is the test that decides whether the nesting *scheme* works, as opposed to
whether it is correctly *implemented* -- which the U1-U34 runmodes, P1-P11 and
the I1-I8 integration matrix in `tests/integration/nesting/` already establish.
Nothing here re-tests those; if they are red, fix them first, because V1 tells
you nothing when a part is broken.

---

## The experiment

A perfect-model, or "Big Brother", protocol:

1. **Parent** ("big brother") -- a periodic, neutral urban LES over an array of
   cubes, driven by a fixed mean pressure gradient, spun up to a statistically
   steady state and then run through a production window during which it dumps
   `u, v, w`.
2. **Reference truth** -- the sub-region of the parent that the child covers,
   extracted from the parent's own dumps.  Not a separate run: literally the
   same field.
3. **Child** -- that same sub-region, at the **same resolution** and with the
   **same buildings**, run with `lnesting = .true.` and driven by boundary
   slabs cut from the parent's dumps.

If the scheme works, the child's statistics match the parent sub-region's, and
whatever error there is stays near the boundary.  Parent and child share the
grid (refinement ratio exactly 1) and the geometry, so nothing here tests the
interpolation -- `slabs_from_fields` *cuts* the slabs rather than interpolating
them.  The interpolation is covered separately by P1-P11 in
`tools/python/tests/test_nesting.py`.

### The building-free zone, and how it is enforced

The design's rule (section 5) is **no buildings anywhere `W > 0`**, plus a
`nest_nwall`-cell erosion around any solid point.  A uniform 32 m cube period
admits no building-free 24 m zone, so the geometry gives way rather than the
zone: the **parent** carries an aligned cube array with the cubes over the
child's guard + ramp band removed -- a plaza.  Every cube is required to be
either entirely inside the child's interior or entirely outside the child, with
a `nest_nwall * dx` margin, so the child's zone sits over open ground.  The
parent is still periodic and still statistically steady; it simply has a square
of open ground where the child's boundary will land, which is the configuration
the scheme was designed for.

The point of doing it this way is that the child can then set

```fortran
nest_lparentgeom = .false.
```

so that `nesting_init` **asserts** the rule and aborts if any solid point has
`W > 0`.  The constraint becomes a checked invariant of every run rather than a
comment in a namelist.  `test_v1_tiny` checks it from the other side too: zero
solid cells in the guard + ramp band, and more than zero in the interior.

`config.py` also offers `geometry = "uniform"` -- the unbroken array, which puts
buildings inside the zone and needs `nest_lparentgeom = .true.`.  That is legal
only because this is self-nesting (the parent resolves the same buildings), it
is the harder case, and it is worth having; it is not the V1 default.

### Why a fixed `dpdx` rather than `luvolflowr`

The two runs have to share the momentum source **exactly**, or the comparison
measures the difference in forcing instead of the difference the boundaries
make.  A volume-flow-rate controller adjusts a body force in time to hold a
target flux, and the child's controller would compensate for precisely the
boundary error being measured.  With a fixed pressure gradient the source is
identical by construction, and the friction velocity
`u* = sqrt(dpdx * zsize)` is known a priori, which gives every error metric a
run-independent scale.

`lscale.inp` therefore carries `pgx = 0` and the forcing lives in `dpdx` in
`&PHYSICS`.  Setting both would double it -- see "Finding N3" below, and the
docstring of `caselib.write_lscale`, which is where the reason is recorded so
that the hand-written `lscale.inp` does not get "simplified" away later.

### Why the child starts from the parent's own field

The nesting file carries the optional schema-2 full-domain initial condition
(`u_init`/`v_init`/`w_init`) and the child sets `nest_linitfromparent = .true.`,
so it cold-starts from the parent's instantaneous field at the first stored
time, made consistent with the corrected boundary data and projected
divergence-free by the writer.  Without it the child's interior turbulence has
to grow in from the boundaries and most of the run is spin-up.

The consequence for the statistics is worth stating: the comparison is
**paired**.  Child and parent share the forcing, the geometry, the grid and the
initial condition, and the boundary keeps them phase-locked at the domain edge
while the interior decorrelates over an eddy-turnover time.  So the
child-minus-parent difference is *smaller* than two independent runs of the same
length would give, and it grows with distance from the boundary purely from
sampling.  Two things follow:

* the sampling-noise floor quoted below (estimated from two independent halves
  of the parent's own window) is an **upper bound** on the noise of the paired
  comparison, so treat a decay length measured against it as lenient;
* the error-versus-distance curve on a **short** window measures decorrelation
  rather than statistical bias.  It only becomes a statement about the scheme
  once the averaging window is long compared with the integral time scale.  The
  production window is 1500 s, about 30 flow-throughs of the 256 m child; the
  `tke_series.png` trace and the sample counts in `v1_metrics.json` are there so
  this can be checked rather than assumed.

---

## Layout

| File | What it is |
|---|---|
| `config.py` | **every** parameter of the experiment, as two `Preset` objects: `tiny` and `production`, including the cube layout.  Run it to print both. |
| `caselib.py` | shared machinery: namelist rendering, the cube-array + ground mesh, the 1-D input profiles, the field-dump reader, the solid mask, the solver launcher |
| `make_parent_case.py` | builds the parent case directory (two namelists: spin-up and production) |
| `make_child_case.py` | cuts the child case -- geometry, `namoptions`, `prof.inp` and `nesting.inp.<nr>.nc` -- out of the parent's dumps |
| `analyse.py` | the comparison: profiles, TKE, spectra, error vs distance, JSON + CSV + PNG |
| `run_v1.py` | end-to-end driver, stage by stage |
| `test_v1_tiny.py` | the `tiny` preset as a unittest -- the harness smoke test |
| `submit_cx3.pbs` | the production job for CX3, 64 cores / 8 h.  **Review before submitting.** |

Nothing is committed as data: the STL, the IBM sparse inputs, the profiles, the
namelists and the nesting file are all regenerated from `config.py` on every
run.  The only committed things are code and this README.

`caselib.cube_mesh` calls two private helpers of
`udgeom.geometry_generation` (`_operate_unit_cube` and
`_generate_ground_matlab_style`) because the public `create_cubes` can only
express the unbroken regular array, not the plaza.  It makes the same two calls
in the same order that `create_cubes(..., 'AC')` makes, and
`test_v1_tiny.test_uniform_layout_reproduces_create_cubes` requires the mesh
built from the *full* layout to match `create_cubes` vertex for vertex, so the
coupling cannot drift silently.  If `create_cubes` ever grows a layout argument,
delete `cube_mesh` and call it.

---

## How to regenerate it

```bash
module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
source /rds/general/user/mvr/home/udales/.venv/bin/activate
./tools/build_executable.sh icl release      # if the binary is not already there
```

Smoke test (about 70 s on a login node, 4 ranks):

```bash
python tests/validation/nesting/run_v1.py $EPHEMERAL/v1-tiny --preset tiny
```

or, as the suite entry does,

```bash
python tests/run_tests.py nesting-validation      # runs BOTH entries, see below
python tests/validation/nesting/test_v1_tiny.py   # just the smoke test
```

Production (**do not run this on a login node**):

```bash
qsub tests/validation/nesting/submit_cx3.pbs      # qsub is in /opt/pbs/bin
```

`run_v1.py` refuses the `production` preset outside a batch job unless given
`--yes`, because it wants 64 ranks and tens of gigabytes.

Individual stages, for a re-run after a failure:

```bash
python tests/validation/nesting/run_v1.py <rundir> --preset production \
    --start-at child-case --yes
```

Stages are `parent-case`, `spinup`, `production`, `child-case`, `child`,
`analysis`.

### Suite registration

`tests/test_suites.yml` gains one group, `nesting-validation`, with two entries:
the tiny smoke test (`cost: fast`) and the production campaign (`cost: slow`),
both `class: experimental`, `kind: system`.  The group is **not** included by
`all`, `supported` or `nesting` -- it has to be asked for by name.  Verify with
`run_tests.py`'s own group expansion; at the time of writing `all` expands to 30
suites and none of them is a `nesting-validation` one.

---

## The parameters, and why

Production preset (`config.py`):

```
parent    256 x 256 x 64 cells, 512 x 512 x 128 m, dx = dy = dz = 2 m
child     128 x 128 x 64 cells, origin (128, 128) m -- the central sub-region
buildings 16 m cubes, h = 16 m (8 cells), 16 m streets -> 32 m period
geometry  'plaza': 228 of 256 cubes; the 28 over the child's zone are removed
zone      L_imp = 6 m (3 cells), L_rel = 18 m (9 cells), tau = 1 s, nzone = 12
interior  104 x 104 cells = 208 x 208 m = 13 h x 13 h, holding 6 x 6 cubes
forcing   dpdx = 1.25e-3 m/s^2 -> u* = 0.4 m/s
schedule  spin-up 3600 s, production [3600, 5400] s, dumps every 3 s
          child spin-up 300 s, statistics over the remaining 1500 s
neutral   ltempeq = .false., lmoist = .false., lbuoyancy = .false.
```

* **The zone width** follows design section 1.4: `N_imp = 3` (one cell for the
  advection stencil, one of margin, one for `ihc = 2`), `N_rel = 9`, total 12
  cells -- the recommended default.  `L_rel = 18 m` clears
  `max(8 dx, h) = 16 m`.
* **`tau = 1 s`** gives an absorption optical depth
  `D = (L_imp + L_rel/2) / (U tau) = 5.0` at the bulk velocity, comfortably past
  the `D >~ 2.3` that design section 1.4(b) asks for; a reflection makes the
  round trip and comes back at `exp(-2D) = 5e-5`.  `config.py` prints `D`.
* **The child window is a whole number of cube periods** and its faces sit
  mid-street.  Two things follow, and both matter: the child's cube layout is
  *identical* to the parent's restriction -- there is only one description of
  the geometry, `Preset.cube_centres()`, and the child's is a window on it, and
  the smoke test still checks the resulting IBM masks cell by cell -- and the
  four lateral boundary planes are entirely fluid, so the offline flux
  correction, which the Python writer applies over all faces, matches the
  solver's `IIu`/`IIv`-masked sum exactly.
* **The zone is square and symmetric.**  104 interior cells in both `x` and `y`,
  so the `x_*` and `y_*` decay lengths carry equal weight; there is no direction
  to discount.
* **`nest_timeinterp = 1` (linear), not 2.**  See "Finding N1" below.  This is
  not a preference.

### Two things not to chase

* **`nesting_init` prints no zone warnings at this size, and that is correct.**
  Its two checks are `L_imp + L_rel < 6 cells` (here 12) and
  `L_imp + L_rel > 15 % of the domain`.  The second compares **metres against
  metres**: 24 m against the child's `xlen = ylen = 256 m` is 9.4 %, so it does
  not fire.  Comparing 24 against the child's 128 *cells* gives 18.75 % and the
  false impression that it should -- it is a cells-versus-metres slip, and worth
  not repeating.  Should a future preset make it fire, note that the warning is
  about the *fraction* of the domain the zone eats, not about the absolute size
  of the interior; 104 interior cells would still be ample.
* **The parent is 16 x 16 cube periods with a plaza in the middle**, so it is
  not horizontally homogeneous, and the plaza is a real feature of the flow.
  That is fine and intended: parent and child see exactly the same geometry, so
  the comparison is unaffected.  It does mean the parent's own statistics are
  not those of an infinite uniform array, and should not be quoted as such.

The `tiny` preset (96 x 96 x 32 parent, 64 x 64 x 32 child, 120 s + 120 s, 24 of
36 cubes, 2 x 2 cubes in the child interior) is a smoke test only.  It runs the
identical code path; it makes no physical claim.

---

## What comes out, and what each number means

`<rundir>/analysis/` after a run:

| File | Content |
|---|---|
| `v1_metrics.json` | everything below, machine readable, including the full curves |
| `profiles.csv` | `z`, `<u>`, `<v>`, resolved TKE, for parent and child |
| `error_vs_distance_{x,y}_{umean,tke}.csv` | the decay-length data |
| `spectrum_z_*.csv` | streamwise spectra at each height |
| `profiles.png`, `error_vs_distance.png`, `spectra.png`, `tke_series.png` | the same, plotted |

Everything is computed on **cell centres**, over **fluid cells only**
(`solid_c.txt`), and over the child's **interior** -- cells whose centre is at
least `L_imp + L_rel` from every lateral face.  Comparing inside the zone would
be circular: the solution there is imposed, so agreement is arithmetic.  The
zone values still appear as the left-hand end of the error-versus-distance
curve, but no criterion is applied to them.

Parent and child are reduced by exactly the same rule.  `modfielddump` writes
`u0(ib:ie, ...)`, so the upper face of each direction is missing from the file;
`caselib.cell_centred` therefore drops the last cell in each direction for both
datasets rather than closing the child's with a boundary value it does not have.

### `profile_metrics`

* `u_rms_difference_over_ustar` -- rms over the column of
  `<u>_child(z) - <u>_parent(z)`, in units of `u*`.
* `tke_rms_difference_over_ustar2` -- the same for resolved TKE
  `k = 1/2 (var u + var v + var w)`, in units of `u*^2`, where the variances are
  temporal at each point and then averaged over the interior fluid cells.
* `*_noise_floor_*` -- the same quantity computed between the **first and second
  halves** of the parent's own window, divided by sqrt(2) (two half-length
  samples differ by sqrt(2) more than a full-length sample does).  This is what
  pure sampling error looks like.  A difference at or below the floor is not
  evidence of anything.

### `error_vs_distance` and `decay_lengths` -- the headline

For each slab normal to `x` (and to `y`), the rms over the fluid cells of that
slab of `child - parent`, normalised by `u*` (mean) or `u*^2` (TKE).  The span
in the other horizontal direction is restricted to that direction's interior, so
a cell inside one zone never contributes to the other's error profile.

`decay_lengths` reports, for each face, the smallest distance beyond the zone
from which **every** cell further in satisfies the criterion, in three variants:

* `vs_threshold_m` -- `error <= 0.05` (of `u*` or `u*^2`): an absolute bar;
* `vs_noise_floor_m` -- `error <= floor`: indistinguishable from sampling noise
  (lenient, see the "paired comparison" note above);
* `decay_length_m` -- `error <= max(threshold, floor)`: the headline number of
  design section 0, which does not demand more accuracy than the reference
  itself carries.

Alongside them: `error_at_zone_edge`, `max_error_outside_zone` and
`median_noise_floor_outside_zone`, so the decay length can be read in context
rather than as a bare number.  `decay_lengths_in_building_heights` divides by
`h = 16 m`.  `None` means the criterion is never met within the available fetch
-- that is a result, not a crash.

A decay length "at or below one cell past the ramp" means the error was already
acceptable everywhere outside the zone, i.e. the measurement is limited by the
zone width, not by the scheme.

### `spectra`

Streamwise spectra of `u'` at three heights (below, at and above `h`), on the
interior `x` span, Hann-windowed, averaged over `y` and over time.  The interior
span is not periodic, so these are windowed spectra, not Fourier series -- but
the identical window and span are used for both datasets, so the bias cancels in
the comparison, which is the only thing being asked of them.  Inside the canopy
the transect passes through buildings; solid cells are set to zero (their
fluctuation is zero anyway) identically in both datasets.
`band_mean_ratio_resolved` is the mean of `E_child / E_parent` over the lower
half of the resolved wavenumber band -- 1.0 is perfect, below 1 means the child
has lost energy.

### `tke_series`

Interior-mean `1/2 <u_i u_i>` at every sampled time, for both runs, with the
statistics start marked.  This is how the choice of `child_spinup` is checked
rather than assumed: if the child's trace has not joined the parent's by the
time the statistics window opens, lengthen `child_spinup` in `config.py` and
re-run the analysis stage.

---

## Pass criterion

V1 passes when, over the child's interior and with the production preset:

1. **Mean profiles.** `u_rms_difference_over_ustar` and
   `v_rms_difference_over_ustar` are at or below the corresponding sampling
   noise floor, or below 0.05 `u*`, whichever is larger.
2. **Resolved TKE.** `tke_rms_difference_over_ustar2` at or below its floor, or
   below 0.05 `u*^2`.
3. **Spectra.** `band_mean_ratio_resolved` within 10 % of 1 at all three
   heights; in particular no systematic loss at high wavenumber, which would
   mean the zone is damping the child's small scales.
4. **Decay length -- the acceptance criterion of design section 0.**
   `decay_length_m` for `x_umean`, `x_tke`, `y_umean` and `y_tke`, measured from
   both faces of each direction, is at most **2 h = 32 m** beyond the inner edge
   of the zone.  Report it in metres and in `h`; it is the number V2
   (zone-width sensitivity) and V5 (parent coarsening) are measured against.
5. **Bookkeeping.**  The child's own diagnostics stay clean: normalised `Phi` at
   round-off throughout, `divmax` at round-off, and the reported
   `|grad p|` zone/interior ratio order 1.

Criteria 1-3 are about the *interior* and say whether the scheme reproduces the
parent at all; criterion 4 is what the design actually asks for and is the
number to quote.

**No tuning.**  Do not narrow the analysis window, move the interior boundary,
or change `tau` to make a number look better.  If V1 fails, that is the result;
V2 exists to vary the zone width deliberately and see what happens.

---

## Findings

### Finding N1 -- the Fritsch-Carlson limiter breaks the flux compatibility

**Found while bringing this harness up, in code this directory does not own
(`src/modnesting.f90`).  Reported, not worked around: the harness selects
`nest_timeinterp = 1` and says so.  A fix is in progress by the owner of
`src/`; this section should be re-pointed at it once it lands.**

Design section 3.1(2) rests on `Phi` being a *linear* functional of the boundary
data, so that an interpolant which is linear **in the data** carries `Phi = 0`
from the stored levels to every intermediate time.  `nest_timeinterp = 1`
(linear) has that property.  `nest_timeinterp = 2` does not -- but the defect is
the **slope limiter**, not cubic Hermite as such.  `hermite()` in
`src/modnesting.f90` uses Fritsch-Carlson limited slopes,

```fortran
if (s1*s2 <= 0.) then
   d2 = 0.
else
   d2 = (wa + wb)/(wa/s1 + wb/s2)
end if
```

-- a weighted harmonic mean with a sign-change branch, both nonlinear in
`y1..y4`.  Each boundary face therefore acquires its own effective weighting and
the cancellation that makes `Phi` vanish is destroyed.  An *unlimited* cubic
Hermite (Catmull-Rom) is linear in the data and preserves `Phi` identically;
monotonicity buys nothing for a sign-unconstrained quantity like a velocity
component, so specifying a monotone limiter was a design error rather than an
implementation one.  Independently confirmed by the owner of `src/`: with
boundary data corrected to `|Phi| ~ 1e-14` at every stored level, across an
interval linear gives 1.1e-14, unlimited Catmull-Rom 1.4e-14, and
Fritsch-Carlson 1.0e+00.

Measured here on the `tiny` preset, with a real turbulent parent and a nesting
file whose stored levels are corrected to `|Phi|/A = 7.8e-15`:

| `nest_timeinterp` | normalised `Phi` at runtime | `divmax` | `divtot` |
|---|---|---|---|
| 1 (linear) | 0 to 8.9e-16 | 4.0e-16 | -1.6e-14 |
| 2 (monotone Hermite) | 5.2e-5 | 2.1e-4 | 1.7 |

With `nest_lfluxassert = .true.` (the default) mode 2 **aborts the run** at the
first substep: `5.2e-5` against `nest_fluxtol = 1e-10`.  With the assertion off
it runs, and the pure-Neumann pressure problem is left incompatible, which is
exactly the Case-A symptom design section 3.2 predicts.

Freezing the parent data in time (all levels equal, so no interpolation happens)
gives `Phi = 2.7e-14` with mode 2, which localises the defect to the
interpolation and rules out the masks, the slab indexing and the stagger --
those were checked independently by imposing known constants on each of the four
boundary planes and recovering each face's fluid area to five digits.

Why the existing tests do not catch it: U28 checks that `Phi` at an interpolated
time equals the interpolation of the endpoint values, which is the *linear*
claim; and I3/I4/I5 use parent data that is constant or nearly constant in time,
for which the limiter is inactive.  A regression test would need a parent whose
boundary flux varies non-monotonically between stored levels -- which is what
any real turbulent parent does.

Until the fix lands, V1 runs with linear time interpolation, which means it does
**not** exercise the Hermite path and cannot measure the once-per-`dt_P`
pressure transient that design section 1.3 predicts for the linear interpolant.
That measurement should be added to V1 once an unlimited Hermite mode is
available.

### Finding N2 -- the offline flux correction is doing real work

The parent's field dumps are `NF90_FLOAT`.  Truncating a discretely solenoidal
field to single precision leaves a net boundary flux of `|Phi|/A ~ 1e-9` m/s on
the tiny preset -- ten times the solver's `nest_fluxtol = 1e-10`.  So the
offline correction is not a formality here: without it every V1 run would abort
at initialisation.  After correction, `|Phi|/A = 1.8e-14`.  The same holds for
the initial-condition block: `divmax` goes from `~1e-7` (the single-precision
truncation) to `~1e-16` after `sync_initial_condition`'s projection.

This is worth keeping in mind if the parent's output is ever promoted to double
precision: the numbers above will change, but the correction should stay.

### Finding N3 -- `UDPrep.forcing` doubles a mean pressure-gradient forcing

Not a nesting issue, but it bit this harness, and it has been confirmed
independently as a production bug in the preprocessing.
`generate_lscale` writes `ls[:, 3] = dpdx` into the `pgx` column of
`lscale.inp` when no forcing switch is set
(`tools/python/udprep/udprep_forcing.py`), while `modstartup` computes
`dpdxl(k) = om23_gs*vg(k) - pgx(k) - dpdx` from `pgx` **and** the `&PHYSICS`
`dpdx` (`src/modstartup.f90:2231`).  A case that sets `dpdx` in `namoptions` and
then runs the standard forcing preprocessing is forced twice.

It is unrelated to nesting and is deliberately **not** being fixed on this
branch.  This harness therefore writes its own `prof.inp` and `lscale.inp` with
`pgx = 0` and uses `UDPrep` only for the IBM/geometry section.  The reason is
recorded in the docstring of `caselib.write_lscale` and at both call sites, so
that the hand-written profiles are not later "cleaned up" into a
`UDPrep.forcing.run_all()` that would silently double the forcing.

---

## Status

The whole pipeline has been driven end to end at the `tiny` preset on CX3
(4 ranks, login node, ~70 s wall time): parent spin-up and production, slab cut,
nested child run, analysis, numbers out.  `test_v1_tiny.py` is 10 tests, all
passing (64 s).  The production preset has been **prepared but not submitted**.

The production parent case has separately been built and checked at full size:
143 s to preprocess, 228 cubes, `nfcts = 12992`, **zero solid cells in the
child's guard + ramp band** and 6 x 6 cubes in its interior; 30 s of simulated
time on 4 ranks took 80 steps in 62.0 s, i.e. 5.4e6 cell-steps/s with a mean
`dt` of 0.379 s and `divmax = 8.3e-16`.  `submit_cx3.pbs` sizes against those
numbers.

What the tiny run showed, for orientation only -- 24 samples over 71 s cannot
support a physical claim:

* the child ran clean with `nest_lparentgeom = .false.`, i.e. `nesting_init`
  verified the zone is building-free rather than warning about it: `Phi` at
  round-off (`-3.7e-16` at the first substep), `divmax` at round-off;
* interior profiles of `<u>`, `<v>` and TKE lie close to the parent's, with rms
  differences of 0.018 `u*`, 0.008 `u*` and 0.052 `u*^2` against noise floors of
  0.092 `u*` and 0.029 `u*^2` -- note the TKE difference is *above* its floor at
  this window length, which is exactly the kind of statement only the production
  run can settle;
* the error-versus-distance curve is smallest at the boundary and largest in the
  middle -- the signature of a paired comparison too short to average out the
  decorrelation, exactly as the "paired" note above predicts, and the reason the
  production window is 1500 s rather than 80 s.
