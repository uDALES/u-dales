# V1 and V2 -- the "Big Brother" nesting validation

This directory implements tests **V1** and **V2** of
`docs/udales-nesting-design.md` section 10.4.

**V1** asks:

> Does matched LES-to-LES nesting reproduce the parent?

It is answered, in section 10.5: yes for the mean flow and yes inside the
canopy, but with a real ~10 % resolved-TKE deficit above the canopy, which that
section attributes to insufficient **fetch** rather than to the boundary
treatment.

**V2** puts that attribution on trial.  It is a *falsification test*, not a
parameter survey, and it is written up in "V2 -- the falsification test" below.
Read that section, not this one, if V2 is what you came for; the V1 material
here is its foundation and its vocabulary.

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
| `analyse.py` | the comparison for **one** child: profiles, TKE, spectra, error vs distance, the V2 falsification metrics; JSON + CSV + PNG |
| `sweep_summary.py` | reduces a whole V2 sweep to one table; CSV + JSON + Markdown + two plots |
| `run_v1.py` | end-to-end V1 driver, stage by stage |
| `run_v2.py` | end-to-end V2 sweep driver: builds, runs and analyses every point out of one parent |
| `test_v1_tiny.py` | the `tiny` preset as a unittest -- the V1 harness smoke test |
| `test_v2_tiny.py` | the `v2-tiny` sweep as a unittest -- the V2 harness smoke test, plus the checks on the production sweep's configuration that need no run |
| `submit_cx3.pbs` | the V1 production job for CX3, 64 cores / 8 h.  **Review before submitting.** |
| `submit_cx3_v2.pbs` | the V2 sweep job for CX3, 64 cores / 8 h, reusing the V1 parent.  **Review before submitting.** |

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

`tests/test_suites.yml` gains one group, `nesting-validation`, with four
entries -- the V1 tiny smoke test (`cost: fast`), the V1 production campaign
(`cost: slow`), the V2 tiny sweep (`cost: medium`) and the V2 production sweep
(`cost: slow`) -- all `class: experimental`, `kind: system`.  The group is
**not** included by `all`, `supported` or `nesting`: it has to be asked for by
name.  Verify with `run_tests.py`'s own group expansion rather than by reading
the YAML:

```python
import sys; sys.path.insert(0, "tests")
import run_tests
m = run_tests._load_manifest()
for grp in sorted(m["groups"]):
    suites = run_tests._expand_groups(m, grp)
    print(grp, len(suites),
          [s["label"] for s in suites if s["label"].startswith("nesting-validation")])
```

At the time of writing `all` expands to 33 suites and none of them is a
`nesting-validation` one; `nesting-validation` expands to exactly the four.

The V2 production entry reuses the V1 parent, so the V1 production entry has to
have run first -- it points `--parent-dir` at `build/nesting-validation-v1/903`.
On CX3, submit `submit_cx3_v2.pbs` instead: it defaults to the `converged` V1 run
in `$EPHEMERAL` and is sized against it.

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
| `tke_deficit.csv` | the resolved-TKE difference profile and the parent's half-window spread |
| `tke_error_vs_fetch_{x,y}_{low,high}.csv` | the TKE error against fetch beyond the zone, per face |
| `spectral_bands.csv` | the child/parent spectral ratio in each named band |
| `profiles.png`, `error_vs_distance.png`, `spectra.png`, `tke_series.png`, `tke_deficit.png` | the same, plotted |

The last four were added for V2 and are emitted by *every* run, V1 included:
`analyse.py` has one code path, and `v1_metrics.json` therefore now carries a
`v2` block as well.  The name is historical -- it is the per-child metrics file
-- and a V2 sweep writes the identical content as `v2_metrics.json` in a
per-point subdirectory.

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

**Both production runs are done.**  `production` (PBS job 3990472, 53 min) and
`converged` (3991175, 4 h 08) have run and are written up in
`docs/udales-nesting-design.md` section 10.5.  The `converged` output is on CX3
at `$EPHEMERAL/nesting-v1-converged` -- 170 GB of parent dumps in `903`, the
child in `904` -- and V2 reuses it rather than re-running a parent.  This
paragraph used to say "prepared but not submitted"; it is kept only so the
reader knows which of the two runs the numbers below come from.

The whole pipeline has also been driven end to end at the `tiny` preset on CX3
(4 ranks, login node, ~90 s wall time): parent spin-up and production, slab cut,
nested child run, analysis, numbers out.  `test_v1_tiny.py` is 10 tests, all
passing (100 s).

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

---

# V2 -- the falsification test

## What is on trial

V1 measured, on the `converged` run, a resolved-TKE deficit of **-9.9 %** above
`z/h = 2`, against a per-height parent half-window spread of 1.1-1.5 %.  It is
scale-selective: at `z/h = 2.06` the child/parent spectral ratio is 1.03 for
`lambda > L/4`, 0.875 for 16-64 m, 0.827 for 8-16 m and 1.05 below `4 dx`.  In
physical space the same thing shows up as a TKE error that decays monotonically
inward and is *still falling* where the far zone truncates it.

Section 10.5 of the design document reads that as a **fetch limitation**: the
imposed large scales come through strong, and what is missing is the
energy-containing middle, which the child has to regenerate for itself and does
not have room to.  The alternative reading is that the **boundary treatment** is
damping those scales -- that the relaxation zone is eating them rather than the
domain being too short to remake them.

Those two readings are not a matter of taste.  They make opposite predictions:

> **P1 -- zone width should barely move the deficit.**  If the boundary
> treatment were at fault, widening or narrowing the relaxation zone would
> change how much energy it removes, and the deficit would track it.
>
> **P2 -- child domain size should move it a lot.**  Less fetch, bigger deficit.

**V2 is built to fail.**  If P1 fails -- if the deficit tracks the zone width --
then the fetch interpretation is wrong, the boundary treatment is implicated,
and section 10.5 needs rewriting.  That is a more valuable outcome than a
confirmation, and nothing in this harness is arranged to avoid it: the analysis
window, the height range (`z/h >= 2`), the spectral bands (16-64 m, 8-16 m) and
the metric (relative resolved-TKE difference against the parent's own
half-window spread) are all fixed to what section 10.5 already used, before any
V2 number existed.  Do not adjust them, and in particular do not change `tau`,
the statistics window or `child_spinup` to make a curve flatter.

## Why the parent is reused, and what that costs

Both arms are driven by the V1 `converged` parent that is already on disk:
`$EPHEMERAL/nesting-v1-converged/903`, 170 GB, 3600 dumped levels, spin-up
10800 s, production window [10800, 21600] s.  It is not re-run.  Two reasons,
and the second is the important one:

* it would cost 3.2 h of the job's 8;
* it would be a **different realisation** of the turbulence.  The differences
  V2 is trying to resolve are a few per cent, comparable with the parent's own
  slow oscillation (section 10.5, "one correction to an earlier reading").
  Driving every child from numerically identical forcing is what makes those
  few per cent attributable to the child rather than to the weather.

`config.Sweep.validate()` turns "these points can share that parent" into a
checked invariant rather than a comment: grid, geometry, forcing, schedule,
experiment numbers and -- the subtle one -- the **plaza window** must all match
the parent preset, and the regenerated cube layout must be cube-for-cube the one
the parent actually ran.  `run_v2.py` additionally checks the reused child's
`manifest.json` against its preset before believing its dumps.

### What the fixed parent allows

The parent carries the V1 plaza: the 28 cubes over the V1 child's guard + ramp
band were removed, leaving **40 m of open ground inside each lateral face** of
the 128-cell child.  A child's zone needs `L_imp + L_rel + nest_nwall*dx` of
that.  Two consequences, both measured by `Preset.building_clearance_available`
rather than assumed:

* **The zone ramp stops at `N_rel = 16`, not the 20 of the design table.**
  `N_rel = 16` needs `6 + 32 + 2 = 40` m, which is exactly what is there: the
  weights reach zero at 38 m and the nearest solid cell centre is at 41 m, so
  `nesting_init`'s "no solid point where `W > 0`" assertion passes with a
  cell of margin and `nest_lparentgeom` stays `.false.`.  `N_rel = 20` needs
  48 m.  Running it would mean either re-running the parent -- which
  destroys the identical-forcing property above -- or putting buildings inside
  the relaxation zone, which changes the boundary treatment and so confounds
  precisely the variable P1 is about.  4 -> 16 cells is still a factor of four
  in `L_rel` and 7 -> 19 cells of total zone thickness, which is ample lever.
  **This is a deliberate substitution and it is the one deviation from the
  brief.**

* **The size arm cannot be kept building-free at all**, and that is a property
  of the parent, not a choice.  The widest street in the cube array is 16 m and
  the zone needs 26 m, so the *only* child whose zone sits over open ground is
  the one the plaza was carved for.  Shifting a smaller child off centre does
  not help -- the free frame is 40 m wide but it is a frame, so no smaller
  square has all four faces inside it.  The 96- and 64-cell children therefore
  run with `nest_lparentgeom = .true.`, which is legal here because this is
  self-nesting: the parent resolves the same buildings, so the imposed field is
  the right field, wakes included.  **Their deficits carry a geometry change as
  well as a fetch change**, they are listed under `confounds` in
  `sweep_summary.json`, they are marked `NO` in the `zone clear` column, and
  they must not be read as pure fetch.  The zone arm is clean and carries no
  such caveat.

Because the size arm is confounded and the zone arm is not, an asymmetry in the
evidence is built in from the start, and it is worth stating plainly: **P1 is
the sharp test and P2 is the corroboration.** A clean result on P1 stands on its
own; a large P2 effect is consistent with fetch but does not by itself prove it.

## The points

```
$ python tests/validation/nesting/config.py --sweep v2
sweep 'v2': 6 points (5 to run, 1 reused), one parent 'converged' (903)
common comparison block 40 cells = 5.00h

key        arms        nr   child     N_imp+N_rel  nzone  interior         zone
ref        zone+size   904  128x128     3+9        12     104 cells 13.00h  clear      9.4%  reuse
nrel4      zone        905  128x128     3+4         7     114 cells 14.25h  clear      5.5%  run
nrel12     zone        906  128x128     3+12       15      98 cells 12.25h  clear     11.7%  run
nrel16     zone        907  128x128     3+16       19      90 cells 11.25h  clear     14.8%  run
size64     size        908   64x64      3+9        12      40 cells  5.00h  BUILDINGS 18.8%  run
size96     size        909   96x96      3+9        12      72 cells  9.00h  BUILDINGS 12.5%  run
```

`ref` is the V1 `converged` child.  It sits in **both** arms -- it is the
`N_rel = 9` point of the zone arm and the `128^2` point of the size arm -- and
it is reused, not repeated, so five children run rather than six.

Everything except the zone width and the child size is held fixed: `L_imp = 3`
cells, `tau = 1 s`, `nest_timeinterp = 1`, `nest_nwall = 1`,
`nest_linitfromparent = .true.`, the schedule, the forcing and `child_spinup =
600 s`.  In particular `tau` is *not* rescaled with the zone width, so the
absorption optical depth `D = (L_imp + L_rel/2)/(U tau)` varies along the ramp:
3.3 at `N_rel = 4`, 5.0 at 9, 6.0 at 12, 7.3 at 16.  All four clear the
`D >~ 2.3` that design section 1.4(b) asks for (the weakest, 3.3, still returns
a reflection at `exp(-2D) = 1.3e-3`), so a trend along the zone arm cannot be
blamed on having lost absorption at the narrow end.  `config.py` prints `D` for
every preset.  Rescaling `tau` to hold `D` fixed would have been the other
defensible choice; it was rejected because it changes two things at once, and
P1 is about the zone width.

`nzone` is the one derived quantity that has to follow the zone width, because
the nesting file must store at least as many cells as the weights are nonzero
over; `config._sweep_child` computes it from `Preset.zone_cells`, so it cannot
be forgotten, and `test_v2_tiny` asserts it at every point of both sweeps.

**The interior shrinks faster than the domain.**  Two zone widths come off
whatever the child's size, so the 64-cell child keeps 40 interior cells (5 h)
where the 128-cell one keeps 104 (13 h).  Every table reports the interior in
cells *and* in building heights next to every result, and the `zone fraction`
column records that `nesting_init` warns above 15 % -- which it does at `size64`
(18.8 %) and at `nrel16` (14.8 % -- just under, so it does not).  **That warning
is expected at the small sizes and is not a fault**: it is about how much of the
domain the zone eats, not about the zone being wrong.  `test_v2_tiny` pins it
down by asserting the warning fires exactly where `Preset.zone_fraction_warns`
says it will.

## What is measured

Per child, against the same parent sub-region, all in `v2_metrics.json` under
`v2`:

| Quantity | Where | What it is |
|---|---|---|
| resolved-TKE difference profile | `v2.tke_deficit` | `(TKE_child - TKE_parent)/TKE_parent` at each height, its mean above `z/h = 2`, and the parent's own **half-window spread** `\|A - B\|/TKE_parent` at each height |
| the same over a common region | `v2.common_block` | the central 40 x 40 cells, the same *physical* block for every point |
| spectral band ratios | `spectra.*.bands` | `mean_of_ratios` and `ratio_of_sums` in 16-64 m, 8-16 m, `lambda > L/4` and `lambda < 4 dx`, at each of the three heights |
| TKE error vs fetch | `v2.tke_error_vs_fetch` | per face, on an abscissa of **fetch beyond the inner zone edge**; the error at fixed fetches, at each face's own maximum fetch, and the fetch at which it first stays at or below the parent's sampling floor |
| mean-flow interior fidelity | `v2.criterion_a` | design section 0 criterion A: `max_interior \|<u>_child - <u>_parent\|/u*` over the four faces, against 0.05 |

Three choices in there are worth defending, because each could have been made
to flatter a hypothesis and was not:

* **The spread is the half-window spread, not a standard error.**  `|A - B|`
  between the parent's two half-window estimates, divided by the full-window
  mean.  That is what section 10.5 quotes (1.1-1.5 % at `z/h = 3-5`, giving the
  V1 deficit 4-8 sigma), so V2's significances are readable against V1's.  It
  is conservative by about a factor of two: the standard error of the
  full-window mean is roughly half of it, and the child-parent comparison is
  paired on top of that.  Quoted significance is therefore a lower bound.

  **Two aggregates of it are reported and neither is chosen after the fact.**
  The band `z/h >= 2` runs all the way to the lid at `z/h = 7.8`, and up there
  the resolved TKE is small and its relative sampling error is large: for the V1
  reference child the per-height spread is 1.2 % over `z/h = 3-5` but 10.5 %
  above `z/h = 6`.  So `mean_spread` (the mean over the band) comes out at 4.5 %
  and gives the V1 deficit 2.2 sigma, while the median spread is 1.7 % and the
  median of the *per-height* significances is 4.0 sigma, closer to section
  10.5's 4-8 for `z/h = 3-5`.  The *deficit* is -9.91 % either
  way -- only the uncertainty aggregation moves -- so both go in the table
  (`sigma` and `sigma med`) and the per-height profile is in `tke_deficit.csv`
  for anyone who wants to look at a specific height.
* **The headline band reduction is `mean_of_ratios`, not `ratio_of_sums`.**
  Also because it is what section 10.5 quotes.  Recomputing the V1 `converged`
  spectra with `mean_of_ratios` gives 1.033 / 0.869 / 0.833 / 1.059 against the
  published 1.03 / 0.875 / 0.827 / 1.05; `ratio_of_sums` gives 1.034 / 0.906 /
  0.815 / 0.928, which is defensible but is not the same number.  Both are
  emitted; the tables show the first.
* **Fetch, not distance from the face, is the abscissa for the TKE error.**
  The zone arm varies the zone width, so at a fixed distance from the face one
  child is in its interior while another is still inside its zone.  Fixed
  fetches of `0.5h`, `1h` and `2h` were chosen to exist for the smallest child
  in the sweep (`size64` has 2.5 h of fetch from each face); anything deeper is
  reported as each face's own maximum.

The 16-64 m and 8-16 m bands are **fixed physical bands**, so children of
different sizes are compared over the same eddies.  The `lambda > L/4` band is
by definition relative to the interior span and is therefore *not* comparable
across the size arm; `sweep_summary` also reports the bands taken over the
common 40-cell block, where the spans do match, as `band_16_64_common` and
`band_8_16_common` in the CSV and JSON.

## How to read the summary table

`<rundir>/analysis/sweep_summary.md` is the deliverable.  It is two tables --
one per arm, with `ref` appearing in both -- ordered along the swept variable,
so **each is read downwards**.

| Column | Read it as |
|---|---|
| `N_rel`, `zone`, `child`, `int.cells`, `int./h` | the configuration; `int./h` is the free fetch in building heights and is the abscissa of P2 |
| `zone clear` | `yes` = building-free zone, `nest_lparentgeom = .false.`; `NO` = the confounded points |
| `dTKE z/h>2 [%]` | **the headline.**  Mean resolved-TKE difference above `z/h = 2`.  V1 measured `-9.9` here |
| `spread [%]`, `sigma`, `sigma med` | the parent's own half-window spread over the same heights, and the deficit in units of it -- `sigma` from the mean spread over the band, `sigma med` from the median.  `sigma med` is the one comparable with section 10.5; `sigma` is the conservative one.  See "the spread is the half-window spread" above |
| `dTKE common [%]` | the same deficit over the **common central block**, which is the same physical region for every point.  Compare this column across the size arm: it separates "shorter fetch" from "smaller measurement window" |
| `dTKE z/h<1 [%]` | inside the canopy.  V1 measured `+1.5` -- the canopy is clean, and a V2 point that spoils it is telling you something |
| `E ratio 16-64 m`, `E ratio 8-16 m` | the scale-selective part, at `z/h ~ 2`.  V1: 0.869 and 0.833 |
| `E ratio > L/4` | the imposed large scales.  V1: 1.033.  Not comparable across the size arm |
| `err@zone edge`, `err@1h`, `err@2h`, `err@max fetch`, `max fetch/h` | the TKE error decay, in `u*^2`, averaged over the four faces |
| `faces at floor` | how many of the four faces' errors reach the parent's sampling floor and stay there.  V1: 2 of 4 -- both `x` faces, neither `y` face |
| `crit. A [u*]`, `crit. A` | the mean flow, worst of the four faces.  V1: 0.039 (design section 0 quotes 0.028, which is the west face alone) against the 0.05 bound.  A `FAIL` here is a regression and outranks everything else in the table |

Then read the `verdict` block:

```json
"zone_arm": { "range_pct": ..., "mean_spread_pct": ..., "range_in_spreads": ..., "slope_pct_per_x": ... },
"size_arm": { ... }
```

`range_in_spreads` is how far the deficit moved along an arm, in units of the
parent's own sampling spread; `range_in_median_spreads` is the same thing using
the median per-height spread, which is not dragged up by the near-lid levels and
is the more discriminating of the two.  The fetch interpretation predicts

* zone arm: `range_in_spreads` small -- of order 1 or less -- and no monotone
  trend in `slope_pct_per_x`;
* size arm: `range_in_spreads` large, with the deficit growing as `int./h`
  falls.

**The falsifying outcome is a zone arm whose `range_in_spreads` is large and
monotone in `N_rel`.**  If that is what comes out, say so: the deficit is being
produced by the relaxation zone, not by the fetch, and section 10.5's
explanation is wrong.  A middle outcome -- both arms moving comparably -- means
the two effects are entangled at this domain size and V2 has not separated
them; report that too, rather than picking whichever arm looks cleaner.

Note also what would make the whole table uninterpretable: a `dTKE z/h>2`
smaller than its `spread`, at every point.  Then the sweep is measuring sampling
noise and neither prediction has been tested.  V1's `converged` window was sized
so that this does not happen (deficit 9.9 % against a 1.1-1.5 % spread); V2
inherits that window unchanged, which is why `child_spinup` and the statistics
window are not swept.

## Running it

Smoke test (about 2 min on a login node, 4 ranks -- runs a tiny parent, a tiny
reference child, then the tiny sweep):

```bash
module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
source /rds/general/user/mvr/home/udales/.venv/bin/activate
python tests/validation/nesting/test_v2_tiny.py
```

The `TestSweepConfiguration` half of that file needs no solver and no run: it
checks the **production** sweep's consistency -- that every point regenerates
the V1 parent's cube layout, that the reference point is shared and reused, that
`nzone` covers every zone, that the common block lands in the same physical
place for every point, that the "cannot share this parent" guard actually
fires, and that `building_free_zone` is measured rather than declared.  Run it
alone with

```bash
python tests/validation/nesting/test_v2_tiny.py TestSweepConfiguration
```

Production (**do not run this on a login node**):

```bash
qsub tests/validation/nesting/submit_cx3_v2.pbs
```

Individual points, for a re-run after a failure:

```bash
python tests/validation/nesting/run_v2.py $EPHEMERAL/nesting-v2 --sweep v2 \
    --parent-dir $EPHEMERAL/nesting-v1-converged/903 \
    --reuse-dir $EPHEMERAL/nesting-v1-converged/904 \
    --only nrel16 --start-at child --yes
```

and re-make the table alone, from the per-point JSON that is already there:

```bash
python tests/validation/nesting/sweep_summary.py $EPHEMERAL/nesting-v2 --sweep v2
```

Stages are `child-case`, `child`, `analysis`, `summary`.

## Cost, from measurement

The V1 `converged` job (3991175, 64 ranks, same machine, same case) measured:

| stage | time |
|---|---|
| slab cut, 128^2 child, `nzone = 12`, 3600 parent levels | 608 s |
| child run, 128^2, 10190 s of simulated time | 1621 s |
| analysis, 3400 parent + 3400 child levels | 1040 s |

The slab cut reads all 170 GB of parent dumps whatever the child's size, so it
splits into ~400 s of fixed read and ~170 MB/s of nesting-file write.  The child
run scales with cell-steps at 1.7e7 cell-steps/s on 64 ranks.  The analysis
splits into ~800 s of parent accumulation and ~240 s per 128^2 child -- and the
parent accumulation is **shared** between points that use the same child window
(`analyse.run`'s `parent_cache`), so it is paid three times, not six.  That
sharing is not only cheaper: it guarantees the whole zone arm is measured
against literally the same parent statistics.

|  | nrel4 | nrel12 | nrel16 | size64 | size96 | sum |
|---|---|---|---|---|---|---|
| slab cut [s] | 525 | 661 | 728 | 505 | 557 | 2976 |
| child run [s] | 1621 | 1621 | 1621 | 650 | 1200 | 6713 |
| child analysis [s] | 240 | 240 | 240 | 60 | 135 | 915 |
| parent accumulation [s] | | | | | | 2400 |
| | | | | | **total** | **~13 000 s = 3.7 h** |

So the brief's "about 3 h" is close; 3.7 h is the measurement-based figure, and
the 8 h walltime leaves 2.2x headroom -- the right amount, since none of the
three sizes below 128^2 has been timed at production volume.

**Memory** peaks in the slab cut, which holds the whole nesting file in memory
before writing it (`udprep.nesting.write_nesting_file` takes arrays, not a
stream).  The widest point, `nrel16` at `nzone = 19`, is 55.2 GB of slabs; the
V1 job peaked at 39 GB for a 35.2 GB file, so budget ~60 GB and `mem=128gb` is
comfortable.  This is the second reason not to widen the ramp past `N_rel = 16`
without thought: `N_rel = 20` would be 66.6 GB.

**Disk**, all new, in `$RUNDIR`: 164 GB of nesting files (20.9 + 43.8 + 55.2 +
17.6 + 26.4) and 173 GB of child field dumps, ~337 GB total.  The V1 run it
reads (246 GB) is not touched.  `--prune-nesting` deletes each child's
`nesting.inp` once it has run, if that ever gets tight -- the analysis does not
need it, only a re-run does.

## The reference row, measured

The `ref` point has been run through the production analysis already -- it needs
no solver, only the V1 dumps -- so the table's anchor is known before the job is
submitted, and the new code is verified at production volume (3398 parent and
3397 child levels, 991 s, `$EPHEMERAL/nesting-v2/analysis/ref/`):

```
V2: N_imp+N_rel = 3+9 cells, child 128x128, interior 104 cells = 13.00 h,
    zone building-free (nest_lparentgeom = .false.)
  resolved-TKE deficit above z/h = 2:  -9.91%  (spread 4.53% mean / 1.67% median;
                                                -2.2 sigma on the mean spread,
                                                -4.0 per-height median)
  inside the canopy (z/h < 1):         +1.50%  against 0.67%
  spectra at z/h = 2.06:   16-64 m 0.869   8-16 m 0.833   > L/4 1.033   < 4dx 1.059
  TKE error vs fetch:  0.5h 0.170,  1h 0.161,  2h 0.155;  at 6.44 h 0.138
  faces whose TKE error reaches the sampling floor: 2/4
  criterion A (mean flow, interior): 0.0386 u* against 0.05 -- PASS
  over the common 40-cell block (5.00 h): deficit -8.54%
```

Every one of those that section 10.5 also quotes reproduces it: `-9.9 %` above
`z/h = 2`, `+1.5 %` in the canopy, `1.03 / 0.875 / 0.827 / 1.05` for the four
spectral bands against `1.033 / 0.869 / 0.833 / 1.059` here.  The differences
are in the fourth decimal of the band edges, not in the physics, and they are
the strongest available evidence that `analyse.py`'s new V2 block measures what
section 10.5 measured rather than something adjacent to it.

The one number that is *new* is the last: over the central 40 x 40 cells the
deficit is **-8.54 %** rather than -9.91 %.  So restricting the measurement to
the middle of a 128-cell child recovers about 1.4 points of the deficit --
which is the size of the "smaller measurement window" effect that the size arm
has to be read against.  If the 64-cell child's deficit over that same block is
much worse than -8.54 %, the extra is fetch; if it is close to it, it was the
window.

## Status

The whole V2 pipeline has been driven end to end at the `v2-tiny` sweep on CX3
(4 ranks, login node): tiny parent, tiny reference child, three swept children
including one with buildings in its zone, per-point analysis, cached parent
accumulation, summary table.  `test_v2_tiny.py` is 21 tests, all passing
(141 s); `test_v1_tiny.py` is unchanged at 10 tests, all passing (100 s), so the
shared `analyse.py` did not regress.  The production sweep has been **prepared
but not submitted**; `submit_cx3_v2.pbs` is sized from the numbers above.

What the tiny sweep showed, for orientation only -- 23 samples over 71 s cannot
support any physical claim, and its half-window spreads (20-55 %) are larger
than its deficits:

* the zone arm moved the deficit by 0.40 percentage points over `N_rel = 2, 4,
  8` (`range_in_spreads = 0.015`), the size arm by 3.58 points between the
  32- and 64-cell children (`range_in_spreads = 0.089`);
* every child ran with `Phi` and `divmax` at round-off, all four faces forced;
* the 32-cell child ran with `nest_lparentgeom = .true.`, `nesting_init`
  reporting 3936 solid points in the zone and the 15 % zone-fraction warning,
  both as predicted rather than as surprises.

Those numbers are the *shape* the production table will have, not a preview of
its content.
