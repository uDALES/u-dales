# V0, V1 and V2 -- the "Big Brother" nesting validation

This directory implements tests **V0**, **V1** and **V2** of
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

**C0** asks whether that deficit is the **boundary cadence**'s doing at all: the
3 s dumps cannot carry the band the child was short of, and the review of
2026-09-06 says so quantitatively.  It is written up in "C0 -- the cadence
discriminator" below, with its predictions registered before the runs.

**V0** asks the question nesting exists for:

> Does a child at higher resolution than its parent reproduce it?

V1 and V2 both run at refinement ratio exactly 1, so neither validates
refinement and neither should be cited as if it did.  V0 is written up in
"V0 -- refinement" below.  Its child is the **V1 child**, unchanged, measured
against the **V1 parent**, unchanged; only the grid the boundary data arrives on
changes.  So the V1 result is the `r = 1` row of V0's table rather than a
separate experiment quoted next to it, and everything in this section applies to
V0 as written.

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
| `config.py` | **every** parameter of the experiments: `Preset` (one case), `Sweep` (V2), `RefinedPoint`/`RefinementSuite` (V0), including the cube layout.  Run it to print them. |
| `caselib.py` | shared machinery: namelist rendering, the cube-array + ground mesh, the 1-D input profiles, the field-dump reader and its flux-conservative coarse view, the `&NESTPARENT` band reader (`NestParent`), the solid mask, the solver launcher |
| `make_parent_case.py` | builds a parent case directory (two namelists: spin-up and production), at whatever resolution the preset asks for; `Preset.parent_output` selects full field dumps, the `&NESTPARENT` band, or both |
| `make_child_case.py` | cuts (`r = 1`) or interpolates (`r > 1`) the child case -- geometry, `namoptions`, `prof.inp` and `nesting.inp.<nr>.nc` -- out of a parent's dumps, full (`--source fielddump`) or band-only (`--source nestparent`) |
| `analyse.py` | the comparison for **one** child: profiles, TKE, spectra, error vs distance, the V2 falsification metrics; JSON + CSV + PNG |
| `analyse_v0.py` | the reductions that only exist at `r > 1`: the spectral ratio split at the parent's Nyquist wavelength, the solver's runtime `Phi`/`divmax`, and the driving parent's own deficit |
| `sweep_summary.py` | reduces a whole V2 sweep to one table; CSV + JSON + Markdown + two plots |
| `run_v1.py` | end-to-end V1 driver, stage by stage |
| `run_v2.py` | end-to-end V2 sweep driver: builds, runs and analyses every point out of one parent |
| `run_v0.py` | end-to-end V0 driver: the coarse parents, the four refined children, the per-point analysis and the suite table |
| `test_v1_tiny.py` | the `tiny` preset as a unittest -- the V1 harness smoke test |
| `test_v2_tiny.py` | the `v2-tiny` sweep as a unittest -- the V2 harness smoke test, plus the checks on the production sweep's configuration that need no run |
| `test_c0_tiny.py` | the `c0-tiny` and `c0b-tiny` sweeps as a unittest -- the C0 harness smoke test, plus the production C0 sweeps' configuration and the experiment-number register |
| `test_v0_tiny.py` | the `v0-tiny` suite as a unittest -- the V0 harness smoke test, plus the coarsening/prolongation invariants and the production suite's configuration |
| `test_nestparent_tiny.py` | the parent-side zone dump (D1) on the `tiny` parent: both outputs at one cadence, bit-identical child builds from each, the storage ratio, the parent's I/O accounting, and a refined child built from a band-only driver |
| `submit_cx3.pbs` | the V1 production job for CX3, 64 cores / 8 h.  **Review before submitting.** |
| `submit_cx3_v2.pbs` | the V2 sweep job for CX3, 64 cores / 8 h, reusing the V1 parent.  **Review before submitting.** |
| `submit_cx3_c0b.pbs` | the C0b job for CX3, 64 cores / 4 h / 128 GB: the fine-cadence parent warm-started from the V1 restart, six children, the table.  **Review before submitting.**  C0a goes through `submit_cx3_v2.pbs` with `UDALES_V2_SWEEP=c0`. |
| `submit_cx3_v0.pbs` | the V0 refinement job for CX3, 64 cores / 6 h, reusing the V1 parent as both reference and filtered-arm source.  **Review before submitting.** |

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

### The parent-side zone dump (D1)

A parent that only drives children need not write full field dumps: with
`Preset.parent_output = "nestparent"` its production namelist carries a
`&NESTPARENT` block instead of `lfielddump` (`"both"` writes the two side by side),
and `src/nesting_parent.f90` writes, per rank, only the band of `nestparent_nzone`
parent cells inside each lateral face of the child box plus one whole-box
initial block (`docs/udales-nesting-spec.md`, section 9).  The band is
`Preset.nestparent_nzone(child)`: the child's guard + ramp on the parent grid,
rounded up, **plus one cell** for the prolongation's slope stencil -- so a
driver sized for a refined child (`make_parent_case.build(...,
nestparent_child=<fine preset>)`) can be prolonged without reaching into the
interior it did not write.  `make_child_case --source nestparent` builds the
same `nesting.inp` from those files, through the same writer; the interior is
`NaN` on the way in, so a slab that reached past the band fails the build rather
than storing zeros.  Because the band has no interior, the `prof.inp` seed and
the manifest's `driving_parent_profile` come from the single initial block.

Measured on `tiny` (`test_nestparent_tiny.py`, about 3 min on a login node): the
two sources give bit-identical nesting files; the band files are 4.0x smaller
than the field dumps (5.1x by cell count for this geometry -- 96^2 parent, 64^2
box, 8-cell band; the production geometry gives ~15x); the parent wrote 35 MB
in 0.06 s of write calls over 40 dumps on 4 ranks.  The V0/V2/C0 drivers keep
`parent_output = "fielddump"`: those parents are also the reference the child is
compared against, and the filtered arm box-filters full dumps.

```bash
python tests/validation/nesting/test_nestparent_tiny.py
python tests/validation/nesting/make_child_case.py <parent_dir> <outdir> --preset tiny --source nestparent
```

### Suite registration

`tests/test_suites.yml` gains one group, `nesting-validation`, holding a tiny
smoke-test entry and a production entry for each experiment in this directory --
V1 (`cost: fast` / `slow`), V2 (`medium` / `slow`) and V0 (`medium` / `slow`) --
all `class: experimental`, `kind: system`.  The group is **not** included by
`all`, `supported` or `nesting`: it has to be asked for by name.  Verify with `run_tests.py`'s own group expansion rather than by reading
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

No `nesting-validation` entry appears in `all`, `supported`, `nesting`,
`experimental`, `python-library`, `supported-macos` or `lint`, and it is that
*absence* rather than any particular count that matters, since the group grows
as experiments are added.  `test_v0_tiny.TestSuiteRegistration` asserts exactly
it for the V0 entries -- both present in `nesting-validation`, neither reachable
from any other group -- so it cannot rot.

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
(`src/nesting_scheme.f90`).  Reported, not worked around: the harness selects
`nest_timeinterp = 1` and says so.  A fix is in progress by the owner of
`src/`; this section should be re-pointed at it once it lands.**

Design section 3.1(2) rests on `Phi` being a *linear* functional of the boundary
data, so that an interpolant which is linear **in the data** carries `Phi = 0`
from the stored levels to every intermediate time.  `nest_timeinterp = 1`
(linear) has that property.  `nest_timeinterp = 2` does not -- but the defect is
the **slope limiter**, not cubic Hermite as such.  `hermite()` in
`src/nesting_scheme.f90` uses Fritsch-Carlson limited slopes,

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

Both arms are clean tests: every point runs the same boundary treatment
(`nest_lparentgeom = .false.`, a building-free zone) so that in each arm exactly
one thing moves.  See "the size arm keeps its zone building-free by clearing it
in the child" below for how that is arranged, and what it costs.

**V2 is built to fail.**  If P1 fails -- if the deficit tracks the zone width --
then the fetch interpretation is wrong, the boundary treatment is implicated,
and section 10.5 needs rewriting.  That is a more valuable outcome than a
confirmation, and nothing in this harness is arranged to avoid it: the analysis
window, the height range (`z/h >= 2`), the spectral bands (16-64 m, 8-16 m) and
the metric (relative resolved-TKE difference against the parent's own
half-window spread) are all fixed to what section 10.5 already used, before any
V2 number existed.  Do not adjust them, and in particular do not change `tau`,
the statistics window or `child_spinup` to make a curve flatter.

## What V1 and V2 do not test: refinement

**Both V1 and V2 run at refinement ratio exactly 1.**  Parent and child share
the grid; `slabs_from_fields` *cuts* the boundary slabs rather than
interpolating them, and the child's `dx` equals the parent's.

That is deliberate -- ratio 1 is what isolates the nesting scheme from the
interpolation, so that when V1 finds a 10 % TKE deficit there is no question of
it being an interpolation artefact -- and the interpolation itself is covered
separately by P1-P11 in `tools/python/tests/test_nesting.py`.

But **the main purpose of nesting is running the child at higher resolution than
the parent**, and that end-to-end case is what **V0** tests -- see
"V0 -- refinement" below, which is where the interpolation is finally exercised
by a running solver rather than by P1-P17 alone.  Nothing in V1 or V2 validates
refinement, and neither should be cited as if it did.  V5 (parent coarsening)
approaches the same question from the other side but is not the same experiment.

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

* **The size arm keeps its zone building-free by clearing it in the child.**
  No child smaller than 128 cells can *inherit* a clear zone from this parent --
  the widest street in the cube array is 16 m and the zone needs 26 m, so the
  only child whose zone sits over open ground by inheritance is the one the
  plaza was carved for, and shifting a smaller child off centre does not help
  (the free region is a 40 m frame, so no smaller square has all four faces
  inside it).  But it does not have to inherit one.

  **Parent and child geometry are not required to match.**  Design section 9.4
  says so explicitly, and V3 and V4 exist precisely to vary them.  So the 96-
  and 64-cell children simply do not carry the cubes that fall in their guard +
  ramp band: `Preset.clear_child_zone` drops them from the *child's* layout
  while the parent keeps them.  The parent's buildings still reach the child --
  their wakes are in the velocity field imposed on the boundary -- so the zone
  gets physically meaningful forcing without containing a single solid cell.
  Every point of the sweep therefore runs `nest_lparentgeom = .false.`, the
  solver *asserts* the design section 5 rule at all six, and **child size is the
  only variable moving along the arm.  P2 is a sharp test, not corroboration.**
  `size64` clears 12 cubes and `size96` 20; the zone arm clears none, because
  its zones fit in the plaza already.

### What clearing the child's zone does and does not confound

Worth being exact about, because it is the one place where the child stops
being a perfect sub-model of the parent.

**What differs.**  Inside the child's guard + ramp band, and only there, the
child is open ground where the parent has cubes.  Two children of the size arm
therefore have a different near-boundary geometry from the parent sub-region
they are compared against.

**Why that is acceptable here.**  The band is where the solution is *imposed*.
The child does not compute the flow there in any meaningful sense -- it is
relaxed onto the parent's -- and no criterion is applied to it: `analyse.py`
measures over the interior, and the error-versus-fetch curves start at the inner
zone edge.  The comparison never asks the child to reproduce the parent
somewhere the two are built differently.

**Two things make that a checked claim rather than a hopeful one.**

1. `Preset.removed_cubes_reaching_the_interior()` must be empty, and
   `validate()` refuses a preset where it is not.  A cleared cube that also
   poked into the analysis interior would put the statistics over two different
   geometries; here none does, because in this array a cube always occupies
   8-24 m in from a child face (centres sit at 16 mod 32, faces at 0 mod 32) and
   the zone is 24 m deep, so a cleared cube stops exactly where the interior
   begins.  This is not a coincidence to be relied on quietly -- it is why the
   tiny sweep carries the production zone rather than a tiny one, and why the
   guard has a test that makes it fire.
2. `analyse.run` compares over cells that are fluid in **both** runs.  A cell
   that is fluid in the child and solid in the parent would otherwise fold the
   parent's near-zero in-building velocity into the parent's statistics.  The
   intersection is a no-op wherever the geometries agree, and
   `v2_metrics.json`'s `solid_mask` block reports how many cells it removed.

**What a reader should still keep in mind: the cleared band is a soft
obstacle, not open ground.**  Inside the band the child is relaxed towards the
parent's velocity field, and that field contains the parent's cubes -- as
near-zero velocity where a cube stands, and as wakes downstream of it.  So the
child's ramp carries a low-velocity imprint of a building that it does not
itself resolve: no IBM enforcing it, no wall stress, no ongoing production.  It
is neither a building nor a plaza.

That is the mechanism by which "the parent's buildings still imprint on the
child", stated precisely, and it has three consequences worth naming:

* **It does not touch the mass budget.**  The compatibility condition is
  evaluated on the domain's boundary faces, and the guard strip -- the first
  3 cells, where `W = 1` -- is over open ground in the parent too, because a
  cube in this array never comes closer than 8 m to a child face.  The offline
  correction still drives the stored `Phi` to round-off and the runtime
  diagnostics stay there; the tiny sweep's cleared child shows exactly that.
* **It does not change the boundary treatment.**  `nest_lparentgeom`, the
  weights, `tau`, the shape function and the guard width are identical at all
  six points.  What differs is *what the imposed field describes*, not how it is
  imposed -- which is precisely the difference between this and the earlier
  arrangement, where the smaller children would have had solid cells inside
  `W > 0` and a different `nest_lparentgeom`.
* **It biases P2 towards confirming the fetch interpretation, not away.**  The
  imposed field injects the parent's wake turbulence at the boundary, but the
  child has no body there to sustain it, so that turbulence decays inward
  instead of being regenerated.  The smaller children therefore get, if
  anything, *less* self-sustaining turbulence than a child with buildings all
  the way to its edge would -- which deepens the deficit at small size.  So a
  **null result on P2 would be the surprising and the more trustworthy
  outcome**, and a large P2 effect should be read with this in mind rather than
  as pure fetch.

Design section 9.4 covers exactly this situation, and V3 exists to measure the
adjustment length it implies; V2 does not measure it, and does not need to,
because the invariant above keeps it out of the region being compared.

The interior also still shrinks faster than the domain: a 64-cell child has 40
interior cells where a 128-cell one has 104, so the smaller points have less
fetch *and* a smaller measurement window.  That is what the common central block
separates, and it is now the only thing left for it to separate.

## The points

```
$ python tests/validation/nesting/config.py --sweep v2
sweep 'v2': 6 points (5 to run, 1 reused), one parent 'converged' (903)
common comparison block 40 cells = 5.00h

key        arms        nr   child     N_imp+N_rel  nzone  interior         zone                   run
ref        zone+size   904  128x128     3+9        12     104 cells 13.00h  clear      9.4%  cut 0    reuse
nrel4      zone        905  128x128     3+4         7     114 cells 14.25h  clear      5.5%  cut 0    run
nrel12     zone        906  128x128     3+12       15      98 cells 12.25h  clear     11.7%  cut 0    run
nrel16     zone        907  128x128     3+16       19      90 cells 11.25h  clear     14.8%  cut 0    run
size64     size        908   64x64      3+9        12      40 cells  5.00h  clear     18.8%  cut 12   run
size96     size        909   96x96      3+9        12      72 cells  9.00h  clear     12.5%  cut 20   run
```

`cut` is how many of the parent's cubes the child does **not** carry, because
they would have fallen in its guard + ramp band.  `clear` in every row is the
consequence: every point runs `nest_lparentgeom = .false.`.

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
| what was averaged over | `solid_mask` (top level) | how many cells were compared, and how many were excluded for being solid in one run but not the other -- nonzero only where the child cleared its zone, and confined to the band |

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
| `zone clear`, `cubes cut` | `yes` everywhere -- every point runs `nest_lparentgeom = .false.` and the solver asserts it.  `cubes cut` is how many of the parent's cubes the child dropped to get there; nonzero means that child is not an exact sub-model of the parent *inside its band*, and nothing else |
| `dTKE z/h>2 [%]` | **the headline.**  Mean resolved-TKE difference above `z/h = 2`.  V1 measured `-9.9` here |
| `spread [%]`, `sigma`, `sigma med` | the parent's own half-window spread over the same heights, and the deficit in units of it -- `sigma` from the mean spread over the band, `sigma med` from the median.  `sigma med` is the one comparable with section 10.5; `sigma` is the conservative one.  See "the spread is the half-window spread" above |
| `dTKE common [%]` | the same deficit over the **common central block** -- the same physical region, 40 x 40 cells at `x, y = [216, 296] m`, for every point.  Compare this column across the size arm: with the geometry confound gone it is the one remaining thing that separates "shorter fetch" from "smaller measurement window", and the V1 reference gives -8.54 % over it against -9.91 % over its own 104-cell interior, so about 1.4 points of the headline number is window rather than fetch |
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

Smoke test (about 3 min on a login node, 4 ranks -- runs a tiny parent, a tiny
reference child, then the tiny sweep):

```bash
module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
source /rds/general/user/mvr/home/udales/.venv/bin/activate
python tests/validation/nesting/test_v2_tiny.py
```

The tiny sweep is **not** built on the `tiny` preset.  It has its own parent,
`tiny-sweep`: a 128 x 128 x 32 parent with a 96-cell reference child, carrying
the *production* zone (3 + 9 cells) on a tiny domain.  It has to.  For a cleared
cube to stay out of the analysis interior the zone must be at least as deep as a
cube sits in from a child face, which in this 32 m array is always 24 m -- the
production zone exactly, and nearly twice `tiny`'s 14 m.  A `tiny`-based sweep
would have had its smallest child's *interior* geometry changed by the clearing,
which `Preset.validate` refuses.  `tiny` itself is untouched, so `test_v1_tiny`
is unaffected.

The `TestSweepConfiguration` half of that file needs no solver and no run: it
checks the **production** sweep's consistency -- that every point regenerates
the V1 parent's cube layout, that the reference point is shared and reused, that
`nzone` covers every zone, that the common block lands in the same physical
place for every point, that every point's zone is building-free, that no cleared
cube reaches the region compared, that the "cannot share this parent" and
"cleared cube reaches the interior" guards actually fire, and that
`building_free_zone` is measured rather than declared.  Run it alone with

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

**Memory** no longer scales with the record.  The slab cut used to hold the
whole nesting file before writing it (`udprep.nesting.write_nesting_file`
takes arrays, not a stream): 55.2 GB of slabs at the widest point, `nrel16` at
`nzone = 19`, and a measured 39 GB peak for V1's 35.2 GB file.  `make_child_case`
now streams each parent level through `udprep.nesting.NestingWriter` as it is
cut -- corrected, residuals stored, initial condition synced and projected
against the corrected level 0 -- and holds nothing but the current level.
Measured on the `tiny` preset with `tracemalloc` (the slab-cut stage alone,
preprocessing excluded): the old path peaked at 161.5 MB for 40 levels and
123.3 MB for 14 (`c0-tiny/cad9`), 1.47 MB per level; the new path peaks at
103.1 MB and 103.0 MB -- flat in the level count.  What remains is not the
record: ~74 MB is retained by the geometry generation and the preprocessing
imports (`cube_geometry` 46 MB, `run_preprocessing` 28 MB, measured per stage)
and ~25 MB is the per-level working set -- one parent level (7 MB at 96^2 x 32)
and the writer's transient copy of it -- which scales with the parent's level
size, not with how many there are.  The output is unchanged: every
slab, the initial condition, `time` and `net_volume_flux` are bit-identical to
the old path's; only the stored post-correction `flux_residual` -- a
cancellation to ~6e-10 m3/s summed in a different BLAS order -- moves by
1.2e-11 m3/s, 1e-16 of the gross boundary flux.  The production sweep therefore
needs the parent level (~100 MB at 256^2 x 64) plus the initial block, not
60 GB; `mem=128gb` is now headroom rather than a requirement, and
`N_rel = 20` is no longer a memory question.

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

Re-verified after the size arm was changed to clear the child's zone: the
reference row is identical to every digit, and its `solid_mask` block reports
zero cells solid in one run but not the other, which is the direct check that
the mask intersection is a no-op where the two geometries agree.

The one number that is *new* is the last: over the central 40 x 40 cells the
deficit is **-8.54 %** rather than -9.91 %.  So restricting the measurement to
the middle of a 128-cell child recovers about 1.4 points of the deficit --
which is the size of the "smaller measurement window" effect that the size arm
has to be read against.  If the 64-cell child's deficit over that same block is
much worse than -8.54 %, the extra is fetch; if it is close to it, it was the
window.

## Status

The whole V2 pipeline has been driven end to end at the `v2-tiny` sweep on CX3
(4 ranks, login node): parent, reference child, three swept children -- one of
which clears 12 cubes out of its own zone -- per-point analysis, cached parent
accumulation, summary table.  `test_v2_tiny.py` is 26 tests, all passing
(202 s); `test_v1_tiny.py` is unchanged at 10 tests, all passing, so the shared
`analyse.py` and `caselib.py` did not regress.  The production sweep has been
**prepared but not submitted**; `submit_cx3_v2.pbs` is sized from the numbers
above.

What the tiny sweep showed, for orientation only -- 23 samples over 71 s cannot
support any physical claim, and its half-window spreads (12-43 %) are larger
than its deficits:

* the zone arm moved the deficit by 3.4 percentage points over `N_rel = 4, 9,
  12` (`range_in_median_spreads = 0.30`), the size arm by 4.3 points between the
  64- and 96-cell children (`0.17`);
* every child ran `nest_lparentgeom = .false.` -- `nesting_init` *verified* the
  zone was clear rather than warning about it -- with `Phi` and `divmax` at
  round-off and all four faces forced, the 64-cell one included, which is the
  point: clearing the child's zone is what lets the smallest child assert the
  rule instead of excusing itself from it;
* the child's solid cells matched the parent's exactly outside the band and were
  zero inside it, at every point;
* the 15 % zone-fraction warning fired at `nrel12` (15.6 %) and `size64`
  (18.8 %) and nowhere else, as `Preset.zone_fraction_warns` predicted.

Those numbers are the *shape* the production table will have, not a preview of
its content -- and note that at this window length the tiny zone arm moves as
much as the size arm does, which is exactly the situation the production
window's 3400 samples exist to resolve.

---
---

# C0 -- the cadence discriminator

> Is the V1 TKE deficit above the canopy caused by the 3 s boundary cadence?

`nesting-plan-2026-09-06.md` section 0 and 1 (the repo's parent directory) is
the specification; this section is the operating manual and the pre-registered
prediction.  V2's zone arm came in with the 8-16 m ratio at `z/h = 2` at
0.825 / 0.815 / 0.810 for `N_rel = 4 / 9 / 12` -- P1 holds -- but both readings
of the deficit predicted that, so it discriminates nothing.  What the review of
2026-09-06 added is a **quantified third reading**: the stored boundary data is
sampled every 3 s and interpolated linearly, and by Taylor's hypothesis that
removes every wavelength below `2 U dt` from the imposed field and attenuates
the octave above it by sinc^4.  Computed from the converged parent's own spectra
and mean wind, the fraction of each band that survives the boundary at 3 s is

| z/h | U (m/s) | 2 U dt at 3 s | 8-16 m: boundary keeps / child has | 16-32 m | 32-64 m |
|---|---|---|---|---|---|
| 0.56 | 1.09 | 6.6 m | 0.62 / **0.99** | 0.87 / 1.00 | 0.96 / 1.02 |
| 1.06 | 2.16 | 12.9 m | 0.14 / **0.97** | 0.59 / 0.98 | 0.86 / 1.01 |
| 2.06 | 3.58 | 21.5 m | 0.00 / **0.82** | 0.21 / 0.85 | 0.66 / 0.95 |

The child regenerates most of what the boundary lost, and the residual deficit
orders band by band and height by height with what was lost.  So the reading
is: **the loss is at the boundary, the fetch is the recovery**, and V1's
"fetch, not the boundary treatment" named the recovery and missed the loss.

The dimensionless number is the **dump Courant number** `C_dump = U dt_P / dx_P`
(`Preset.dump_courant`; `Preset.describe()` prints it at `u0`).  Nothing the
parent resolved is lost at the boundary when `2 U dt_P <= 4 dx_P`, i.e.
`C_dump <= 2` at the largest wind in the zone.  V1 ran at 1.6 (z/h = 0.6), 3.2
(z/h = 1), 5.4 (z/h = 2), about 7.5 at the lid.

## What C0 does

`Preset.cadence` (seconds, default `dtdump`, must be a whole multiple of it)
is the interval of the boundary data handed to the child.  `make_child_case`
reads every `cadence / dtdump`-th parent dump level **and never opens the
rest** -- the levels are read by index from the per-rank dump files, so a 6 s
child costs half the I/O of a 3 s one -- and everything downstream (`n_use`,
`runtime`, `t_offset`, `parent_dt`) is taken from the subsampled axis.  The
manifest records it under `cadence`: `seconds`, `stride`, `parent_dt`,
`n_levels_dumped`, `n_levels_used`, `C_dump_at_u0`.

Two arms, both `Sweep`s with a single `"cadence"` arm, through `run_v2.py`:

**C0a -- `c0`, coarser cadences from the existing dumps.**  Parent = the
converged 903 (170 GB of 3 s dumps, reused exactly as V2 reuses it); reference
= the converged child 904, reused.  Three children: `cad6` (961, 6 s, linear),
`cad9` (962, 9 s, linear) and `cr3` (963, 3 s, `nest_timeinterp = 2`).  Mode 2
is now the **unlimited** Catmull-Rom cubic Hermite -- the Fritsch-Carlson
monotone limiter that broke `Phi = 0` in Finding N1 is gone, and
`test_c0_tiny` checks that the CR child's `Phi` and `divmax` sit at round-off
like the linear ones'.  A smoother interpolant cannot restore a band the samples
do not contain, so CR should move 16-32 m a little and 8-16 m not at all.

**C0b -- `c0b`, the fine-cadence ladder.**  The converged parent left its
end-of-spin-up restart (`initd00031204_*.903`, t = 10800 s).  `config.C0_FINE`
(960) **continues** it: `make_parent_case.build(..., restart_dir=...)` symlinks
the 64 restart files into the new case under the new experiment number
(`readrestartfiles` builds each rank's file name from `startfile` by overwriting
the rank fields, so the extension follows the namelist), writes a warm-start
namelist with `startfile` set, and no spin-up phase; the run is 2400 s with
`tfielddump = 0.5` (dt is about 0.38 s, so every 1-2 steps; 4800 levels,
240 GB).  Six children (964-969) are sliced from those dumps at 0.5, 1, 1.5, 3,
6 and 9 s, linear, each with the 600 s discard and an 1800 s statistics
window.  One parent realisation drives the whole ladder, so the comparison is
paired, and the 3 s point cross-checks V1 and C0a.  The children dump every
3 s as V1's did (`Preset.child_dtdump`); the analysis samples the parent every
6th level to match (`Preset.analysis_parent_stride`).  `run_v2.py
--parent-restart-dir` builds and runs the parent first when the run directory
has no dumps for it, and reuses it when it has.

**Primary metric: the band ratios**, 8-16 m and 16-64 m, at every sampled
height (`z/h = 0.56, 1.06, 2.06`).  Between the converged run's 1491 s and
10191 s windows they reproduced to 0.004, so 1800 s is enough for them.  It is
*not* enough for the profile deficit above `z/h = 2` (the parent's own
half-window spread at that window length is of the deficit's size), which is
reported but secondary.  The analysis's bands are section 10.5's; the plan's
16-32 m prediction is read against the 16-64 m column, where the 32-64 m half
(0.66 kept at 3 s) dilutes it.

## Pre-registered predictions

8-16 m and 16-32 m band ratios at `z/h = 2`, plus the height at which the
deficit first exceeds 3 %.  Written before any C0 number existed.

| | 9 s | 6 s | 3 s | 1.5 s | 0.5 s | CR at 3 s |
|---|---|---|---|---|---|---|
| cadence causes it | < 0.75 | about 0.78 | 0.82 (V1) | about 0.9 | **>= 0.97** (the z/h = 1 value) | 16-32 m up a little, 8-16 m unchanged |
| something in the scheme causes it | 0.82 | 0.82 | 0.82 | 0.82 | 0.82 | 0.82 |
| both | falls with dt | | | | plateau above 0.82, below 0.97 | |

**Outcome (2026-09-07, jobs 3993705 C0a and 3993706 C0b, both exit 0).** The
cadence causes it. 8-16 m / 16-64 m ratios at z/h = 2, one fine parent
realisation (960) sliced to six cadences: 0.935 / 0.971 at 0.5 s, 0.920 / 0.960
at 1 s, 0.914 / 0.950 at 1.5 s, 0.831 / 0.873 at 3 s, 0.713 / 0.725 at 6 s,
0.663 / 0.650 at 9 s; deficits -2.1, -3.3, -4.9, -11.3, -20.5, -24.5 %. C0a's
6 s and 9 s subsamples of the old dumps reproduce the C0b points to 0.002, and
the 3 s point reproduces V1. The 16-64 m prediction for the 0.5 s point is met
(0.971 >= 0.97); 8-16 m falls just short (0.935), which is the linear
interpolant's own sinc^4 attenuation of an 8 m eddy at 0.5 s (about 0.84 at
the boundary), not the scheme. The prediction that Catmull-Rom would leave
8-16 m unchanged was wrong: it moved 0.833 -> 0.892 and halved the deficit
(-9.9 -> -5.7 %), so `nest_timeinterp = 2` is now the default and sweep `c0c`
(expnrs 970-975, the same six cadences with the cubic, off the same 960 dumps)
gives its operating curve. The criterion-A flags on the C0b rows are the
600-sample mean-flow floor (spread 8.2 %), not a finding. Full tables:
`$EPHEMERAL/nesting-c0a/analysis/sweep_summary.md` and
`$EPHEMERAL/nesting-c0b/analysis/sweep_summary.md`; design section 10.5.

The 1 s point of C0b sits between the 1.5 s and 0.5 s columns (about 0.93 under
the first row).  `C_dump` at `u0 = 3 m/s`: 13.5, 9, 4.5, 2.25, 1.5, 0.75 for
9, 6, 3, 1.5, 1, 0.5 s; at the z/h = 2 wind of 3.58 m/s the 0.5 s point is
0.9 and the 1 s point 1.8, both inside the criterion.

**Decision rule.**  If the 0.5 s point reaches >= 0.97: the cause is cadence;
design sections 0, 1.3, 6.2, 9.4, 10.5, this README's V1 status and the PR
text are rewritten, `C_dump <= 2` becomes a stated requirement, the writer
warns when a file violates it (plan item W5), and the fine parent becomes the
V1 reference for everything after.  If a plateau remains above 0.82 and below
0.97: that plateau is the scheme's own deficit, and the next arm is `tau` in
{0.5, 1, 4} s and the guard width.  If nothing moves, the review's hypothesis
is refuted and section 10.5 stands as written, minus the word "not".

## Layout

| File | What it adds |
|---|---|
| `config.py` | `Preset.cadence`, `child_dtdump`, `cadence_stride`, `dump_courant`, `analysis_parent_stride`; `Sweep.arms`; the `c0`, `c0b`, `c0-tiny`, `c0b-tiny` sweeps and `C0_FINE`; the experiment-number register (C0 owns 960-969) |
| `make_child_case.py` | subsamples the parent levels before reading; the manifest's `cadence` block |
| `make_parent_case.py` | `--restart-dir`: a parent warm-started from another run's restart set |
| `run_v2.py` | the `parent` stage and `--parent-restart-dir`; the reuse guard checks the cadence |
| `analyse.py` | samples the parent at `analysis_parent_stride` |
| `sweep_summary.py` | the `ARMS` registry (`zone`, `size`, `cadence`); the per-height band table and `sweep_cadence_bands.png` for the cadence arm |
| `test_c0_tiny.py` | both tiny sweeps end to end, plus the production sweeps' configuration and the register |
| `submit_cx3_c0b.pbs` | the C0b job: parent + six children + table, 64 cores / 4 h / 128 GB |

## Running it

Smoke test (about 10 min on a login node, 4 ranks: the tiny V1 parent and
child, three C0a children, the warm-started 0.5 s parent, three C0b children):

```bash
module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
source /rds/general/user/mvr/home/udales/.venv/bin/activate
python tests/validation/nesting/test_c0_tiny.py
```

`TestC0Configuration` alone needs no solver: the four sweeps' consistency, that
C0a is the V1 child at other cadences off the V1 parent, that `C0_FINE` is the
converged parent with only the schedule changed and `t_start = 10800`, that
960-969 collide with nothing, that a cadence that is not a multiple of `dtdump`
is refused, and the registration.

Production.  C0a goes through the V2 job with the sweep and run directory
overridden (the script's own 8 h walltime is overridden on the command line;
`qsub -v` values may not contain commas):

```bash
qsub -N udales-nesting-c0a -l walltime=03:00:00 \
     -v UDALES_V2_SWEEP=c0,UDALES_V2_RUNDIR=$EPHEMERAL/nesting-c0a \
     tests/validation/nesting/submit_cx3_v2.pbs
qsub tests/validation/nesting/submit_cx3_c0b.pbs
```

Results: `$EPHEMERAL/nesting-c0a/analysis/sweep_summary.{md,csv,json}` and
`$EPHEMERAL/nesting-c0b/analysis/...`; the per-height band table is the
"Band ratios (child/parent) at every sampled height" block of the Markdown and
the `band_8_16m@...` / `band_16_64m@...` columns of the CSV, with
`sweep_cadence_bands.png` next to them.  A single point is redone with

```bash
python tests/validation/nesting/run_v2.py $EPHEMERAL/nesting-c0b --sweep c0b \
    --only cad1.5 --start-at child-case --yes
```

## Cost

C0a: the slab cut reads 1/2, 1/3 and all of the 170 GB (about 5, 4 and 10 min),
the three children run 1621 s each as V1's did, the four analyses share one
800 s parent accumulation; about 2.3 h against 3 h.  C0b: `submit_cx3_c0b.pbs`
carries the sizing -- parent 45 min including 240 GB of dump I/O, slab cuts
53 min, children 37 min, analysis 12 min, about 2.6 h against 4 h.  The 0.5 s
point's 47 GB nesting file is streamed level by level, not held (see V2's
"Cost, from measurement"), so `mem=128gb` is far more than it needs.

## Status

Both tiny sweeps have been driven end to end on a CX3 login node (4 ranks)
through `run_v2.py`: the tiny V1 parent and child, then `c0-tiny` (ref reused,
`cad6` cut from every 2nd of the 40 dumped levels, `cad9` from every 3rd,
`cr3` at `nest_timeinterp = 2`), then `c0b-tiny` (the 0.5 s parent 960
warm-started from 903's spin-up restart -- 4 symlinked restart
files, `lwarmstart = .true.`, dumps starting at t = 120 s every 0.437 s median
-- and three children at 0.5 / 1.5 / 3 s cut from 240 / 80 / 40 of its 240
levels).  `test_c0_tiny.py` is 18 tests, all passing (8 configuration tests in
0.2 s; 10 end-to-end in 212 s); every child, the Catmull-Rom one included, kept
`Phi` below 1e-9 and `divmax` below 1e-10, and the parent/child sample counts
of the C0b analysis agree to within 2 with the parent sampled every 6th level.
`test_v1_tiny.py` and `test_v2_tiny.py` are unchanged and pass.

What the tiny sweeps showed, for orientation only -- 20-26 samples over 60-80 s
with a 20-25 % half-window spread support no physical claim -- the 8-16 m and
16-64 m ratios at `z/h = 2.06` were: C0a `ref` 0.520 / 0.647, `cr3` 0.640 /
0.777, `cad6` 0.282 / 0.352, `cad9` 0.186 / 0.225; C0b `cad0.5` 0.799 / 0.937,
`cad1.5` 0.686 / 0.854, `cad3` 0.487 / 0.644.  Monotone in the cadence at every
height, and the 3 s points of the two sweeps agree with each other -- which is
the *shape* the production table is predicted to have, at a window far too
short to say anything about its content.

Production: **submitted 2026-09-06, not yet run.**

| job | what | walltime | results |
|---|---|---|---|
| 3993705 | C0a, `c0`: 961 (6 s), 962 (9 s), 963 (CR 3 s) off the converged 903 dumps, against the reused 904 | 03:00:00 | `$EPHEMERAL/nesting-c0a/analysis/sweep_summary.{md,csv,json}`, `sweep_cadence_bands.png` |
| 3993706 | C0b, `c0b`: parent 960 warm-started from `nesting-v1-converged/903/initd00031204_*.903`, 2400 s at 0.5 s; children 964-969 at 0.5 / 1 / 1.5 / 3 / 6 / 9 s | 04:00:00 | `$EPHEMERAL/nesting-c0b/analysis/sweep_summary.{md,csv,json}`, `sweep_cadence_bands.png` |

The two jobs are independent and use the V1 binary (`build/release/u-dales` of
2026-09-06), so the comparison with V1 and V2 stays paired; the merge-blocker
work of the plan's section 3 does not touch them.

---
---

# V0 -- refinement

> Does a child at higher resolution than its parent reproduce it?

This is the use case nesting exists for, and until now nothing had tested it end
to end.  Design section 10.4 lists V0 first for that reason, not because it is
next in sequence.

Two things are being tested at once and the experiment is arranged so they can
be told apart:

* **refinement** -- a child on a finer mesh than the data driving it;
* **the interpolation** -- `tools/python/udprep/nesting.py`'s conservative
  prolongation (design section 1.3), which at `r = 1` is the identity and has
  therefore never carried a running simulation.  P1-P17 cover it as a function;
  nothing has covered it as part of a solver.

## What is held fixed, and why that matters

**The child is the V1 child.**  Not a child like it: the same 128 x 128 x 64
cells at `dx = 2 m` over the same 256 x 256 x 128 m box, the same 36 cubes, the
same `L_imp = 6 m` + `L_rel = 18 m` zone, `tau = 1 s`, `nest_timeinterp = 1`,
the same `dpdx`, the same cold start from the parent block, the same 10 800 s
window.  `config.RefinementSuite.validate` checks every one of those against
`CONVERGED` rather than trusting the presets to agree, and
`test_v0_tiny.test_every_child_is_the_v1_child` checks it again.

**The reference truth is the V1 parent**, `$EPHEMERAL/nesting-v1-converged/903`,
read and never written -- the same 2 m unnested run, the same sub-region, the
same accumulation code, the same window.

So the only thing that changes is the mesh the boundary data lives on, and
**the V1 result is the `r = 1` row of V0's table**: `-9.91 %` resolved TKE above
`z/h = 2`, `+1.50 %` in the canopy, `0.0386 u*` on criterion A, spectral band
ratios `1.033 / 0.869 / 0.833 / 1.059`.  Every V0 row is read against those
numbers directly, produced by the same `analyse.py` over the same window.

A consequence worth stating: because parent and child are on different grids but
child and *reference* are on the same one, **no regridding enters the
comparison**.  The child's spectrum and the reference's share their `k` axis bin
for bin.  Nothing is coarse-grained to make the two comparable, so nothing about
the measurement can be blamed for what it measures.

## The reference truth: why it is the fine run, in both arms

The obvious-looking choice -- compare the child against the coarse parent that
drove it -- is wrong, and wrong in a way that would flatter or damn the scheme
arbitrarily.  A 2 m child *should* carry turbulence a 4 m or 8 m parent cannot
represent; scoring it against the parent would count that as error.  The only
defensible truth for a refined child is a run at the child's own resolution, and
one already exists.

## The two arms, and why both

The user's question was which reference configuration to use, and the answer is
that the two candidates measure different things and the difference between them
is itself the interesting quantity.  Both are run and both are reported.

### `filtered` -- a perfect coarse parent

The boundary data is the fine reference's own field, box-filtered onto the
coarse grid by `caselib.coarsen_staggered`, then prolonged back onto the child.

* **It is an idealisation, and is labelled as one.**  A box-filtered fine field
  is a coarse field that knows *exactly* what the fine run was doing at every
  scale it can hold -- filtered fine-scale information a genuinely coarse LES
  would never have had, because a coarse LES has to model those scales rather
  than filter them.
* **What it buys is that it is paired.**  Same realisation, same eddies, same
  slow modes.  It is V1 with exactly one variable changed, so the difference
  between its numbers and V1's is attributable to the prolongation and the
  parent's filter scale, and to nothing else.
* **What it can claim:** whether the interpolation and the parent's cutoff cost
  the child anything.  **What it cannot claim:** anything about running off a
  real coarse model.

### `coarse` -- a genuine coarse LES

The boundary data comes from an actual run at 4 m (`911`, 128 x 128 x 32) or 8 m
(`912`, 64 x 64 x 16) over the same 512 x 512 x 128 m domain, with the same 228
cubes, the same `dpdx` and the same schedule.

* **This is the real use case.**
* **It is not paired.**  A separate run is a separate realisation, so its eddies
  are not the reference's and only statistics can be compared.  Read the
  error-versus-distance curves accordingly: V1's were paired, and part of why
  its interior error was so small is that its child and its truth shared an
  initial condition and a history.  A coarse-arm child cannot share those.
* **It carries its parent's biases.**  A 16 m cube is 4 cells wide at 4 m and 2
  cells at 8 m, and its drag will not be the 8-cell version's.  A mean-flow
  error measured against the fine truth is therefore *the parent's error plus
  the nesting's*, and quoting it as if it were the nesting's would be wrong.

### Why both, together

`make_child_case` records the **driving parent's own** interior profiles while
it is already reading every level, and `analyse_v0.parent_deficit` compares them
against the same fine truth.  That gives the ceiling: a child cannot be more
right than the data it is given.  Then, at one ratio,

    filtered arm's parent deficit  =  what the filter threw away
    coarse arm's parent deficit    =  that, plus what the coarse LES got wrong
    coarse minus filtered          =  the coarse LES's own error
    child deficit minus parent's   =  what the nesting cost on top

and `v0_summary.json`'s `filtered_vs_coarse` block reports exactly that
decomposition, per ratio.  Neither arm alone can produce it.

## The grids

| | reference / child | `r = 2` parent | `r = 4` parent |
|---|---|---|---|
| cells | 256 x 256 x 64 / 128 x 128 x 64 | 128 x 128 x 32 | 64 x 64 x 16 |
| `dx` | 2 m | 4 m | 8 m |
| cube width | 8 cells | 4 cells | 2 cells |
| Nyquist wavelength `2 dx_P` | -- | 8 m | 16 m |
| `L_rel / 2 dx_P` | -- | 2.25 | 1.13 |
| ranks | 64 / 64 | 64 | 16 |

`r = 4` is the writer's validated maximum (`MAX_SPATIAL_REFINEMENT`), so the
suite spans the whole supported range.

**The zone does not change between V1 and V0**, and it does not have to:
design section 1.4(c) asks for `L_rel >= max(8 dx_child, h/2..h, 2 dx_P)`, and
`L_rel = 18 m` clears `2 dx_P` at both ratios.
`test_v0_tiny.test_the_production_ramp_is_resolved_by_both_parents` asserts it
rather than leaving it to be noticed.  Keeping the zone fixed is also what makes
the V1 row comparable, so it would have been worth some cost; it happens to cost
nothing.

**Time is not coarsened.**  `dtdump = 3 s` at every ratio, as in V1.  V0 varies
space only; temporal coarsening is V5's question and mixing the two would make
neither answerable.

**The coarse parents run at the same `dtmax = 0.5 s` as the fine one**, so they
integrate at a smaller Courant number than they need to.  That is deliberate --
the mesh is meant to be the only difference between the parents -- and it costs
about 1 400 s of the job.  `RefinedPoint.validate` enforces the equality.

## What makes the geometry comparison legitimate

Buildings are 16 m cubes on 16 m streets, so every face is aligned to all three
grids and the same cube array is representable exactly at 2, 4 and 8 m.  That is
checked, not assumed, in two places:

* `RefinedPoint.validate` requires the coarse parent's cube layout to be the
  fine reference's, cube for cube, and its plaza window to be the same object --
  the plaza is carved in metres, and a coarse preset left to compute its own
  `nest_nwall` margin would carve a slightly *larger* one and the child's
  buildings would stop being the parent's;
* `test_v0_tiny.test_the_coarse_parent_is_the_fine_one_block_averaged` requires
  the coarse run's IBM solid mask to be the fine run's, block-ANDed down, **cell
  for cell**, after both have been through the real preprocessing.  At the tiny
  scale that is 1 536 of 1 536 solid cells at `r = 2` and 192 of 192 at `r = 4`.

## The two refinement-specific measurements

### Where the child's spectrum sits relative to the parent's filter scale

Section 10.5 located V1's deficit in the 8-64 m band.  At `r = 2` the parent's
Nyquist wavelength is 8 m and at `r = 4` it is 16 m, so that band **straddles
the cutoff at both ratios** and a single number for it would average two
physically different situations.  `analyse_v0.nyquist_split` therefore reports
the child/parent ratio in three bands defined relative to `lambda_N = 2 dx_P`:

| band | wavelengths | what the child is doing there |
|---|---|---|
| `parent_resolved` | `>= 4 dx_P` | reproducing structure the parent had |
| `parent_marginal` | `2 dx_P` to `4 dx_P` | representable by the parent but badly damped |
| `sub_parent_filter` | `< 2 dx_P` | **generating** structure the parent never resolved |

plus `contrast = sub_parent_filter - parent_resolved`, negative when the child
is worse where it has to invent and positive when it is worse where it is being
told.  The three bands are asserted to partition every resolved mode, so nothing
falls between them.  The fixed physical bands of section 10.5 (16-64 m, 8-16 m)
are still reported alongside, unchanged, so the V1 row stays readable.

This is the measurement whose outcome is genuinely not predictable in advance:
the child has more capacity to build its own inertial range than the parent had,
but it also has a larger gap to bridge.

### Whether the divergence-preserving prolongation holds in the running solver

Design section 1.3 claims something stronger than conservation: because the
prolongation is piecewise constant tangentially and *linear* normally, each of
`du/dx`, `dv/dy`, `dw/dz` is constant inside a parent cell and equal to the
parent's, so the child target carries the parent's discrete divergence cell by
cell.  Three checks, offline and online:

* **The filter.**  `coarsen_staggered` takes the mean of the `r x r` fine faces
  co-planar with each coarse face, so a coarse face flux is exactly the sum of
  the fine ones and a coarse cell's net flux is the sum of the `r^3` fine cells'.
  Verified on a field built from a vector potential: fine `divmax` `5.0e-16`,
  coarse `1.1e-16` at `r = 2` and `1.7e-17` at `r = 4`.  It is also the exact
  left inverse of the prolongation's tangential half, so `coarsen(prolong(x))
  == x` to round-off -- checked, and it matters because otherwise a round trip
  would move the field and V0 would be measuring the round trip.
* **The prolongation, on a real field.**  `make_child_case` records the parent's
  own `divmax` over the child window next to the `divmax` of the field prolonged
  from it.  On the tiny suite they agree to every digit printed -- e.g.
  `4.9429218051955104e-08` from both sides at `r = 2` -- and the smoke test
  compares their *ratio* against 1, not merely bounds them.
* **The running solver.**  `analyse_v0.runtime_diagnostics` parses the child's
  log for `nesting`'s `Phi` and `modpois`'s `divmax`/`divtot`, plus the zone
  misfit and the zone/interior pressure-gradient ratio that design section 7
  names as the C1 diagnostic.  On the tiny suite, `max |Phi| <= 9.4e-14` and
  `max divmax <= 9.1e-16` at both ratios and on both arms.

## Cost, and the walltime

Sized from the V1 converged job's measured stage times plus two rates measured
directly for this on the real 903 dumps.  The number that mattered most was the
one nobody had: **refinement costs almost nothing in the slab cut**.  Coarsening
plus interpolating twelve slabs adds 0.056 s/level at `r = 2` and 0.025 s/level
at `r = 4`, against a V1 cut of 0.012 s/level of compute inside 0.17 s/level of
I/O -- so `+160 s` and `+46 s` over the whole 3 600-level record.  The slab cut
stays I/O bound at `r > 1`.

| stage | s |
|---|---|
| coarse parent 911, 4 m, 22.6 G cell-steps, 64 ranks | 2 500 |
| coarse parent 912, 8 m, 2.8 G cell-steps, 16 ranks | 1 500 |
| slab cut `r2-filtered` (608 read/write + 160 refinement) | 770 |
| slab cut `r4-filtered` (608 + 46) | 655 |
| slab cut `r2-coarse` (reads 22.6 GB, not 170) | 500 |
| slab cut `r4-coarse` (reads 2.8 GB) | 400 |
| 4 child runs at 1 621 s | 6 484 |
| reference accumulation, 3 400 levels, **shared** by all four | 800 |
| 4 child analyses at 240 s | 960 |
| case preprocessing | 360 |
| **total** | **14 929 s = 4.15 h** |

`submit_cx3_v0.pbs` asks for **6 h**, 1.45x that -- less headroom than the V2
job's 2.2x, deliberately.  The queue is deep and every hour of walltime asked
for is paid in queue time; the only genuinely new items are the two coarse
parent runs and the coarse arm's slab cuts, 4 900 s together, so being 2x wrong
about all of them still fits.  The header says how to split it into per-arm jobs
of 2.5 h and 3.5 h with `-v UDALES_V0_POINTS=...`, which is the right thing to
do when the queue is deeper still: the points are independent once the coarse
parents exist.

**Memory** no longer peaks in the slab cut (the 34.0 GB of slabs is streamed,
see V2's "Cost, from measurement"); `mem=128gb`.  **Disk**, all new: 25 GB of coarse
parent dumps, 136 GB of nesting files, 181 GB of child dumps, ~343 GB.
`--prune-nesting` drops the nesting files once their children have run.

One efficiency deliberately not taken: the two `filtered` points each read all
170 GB of 903 independently.  Building both in one pass would save ~400 s and
would mean a second code path through `make_child_case`; it is not worth it.

## How to run it

```bash
module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
source /rds/general/user/mvr/home/udales/.venv/bin/activate
```

Smoke test (about 160 s on a login node, 4 ranks -- it builds and runs its own
tiny fine reference, two tiny coarse parents and four tiny children):

```bash
python tests/validation/nesting/run_v0.py $EPHEMERAL/v0-tiny --suite v0-tiny
python tests/validation/nesting/test_v0_tiny.py     # the same, as 26 assertions
```

Production (**do not run this on a login node**):

```bash
qsub tests/validation/nesting/submit_cx3_v0.pbs
```

Individual points and stages, for a re-run after a failure:

```bash
python tests/validation/nesting/run_v0.py $EPHEMERAL/nesting-v0 --suite v0 \
    --reference-dir $EPHEMERAL/nesting-v1-converged/903 \
    --only r2-coarse --start-at child-case --yes
```

Stages are `reference-case`, `reference-spinup`, `reference-production`,
`driver-case`, `driver-spinup`, `driver-production`, `child-case`, `child`,
`analysis`, `summary`.  With `--reference-dir` the three `reference-*` stages do
nothing.  `--start-at summary` rebuilds the table from the per-point JSON
without redoing any analysis.

## What comes out

Per point, in `<rundir>/analysis/<key>/`: everything `analyse.py` writes for a
V1/V2 child (`profiles.csv`, `error_vs_distance_*.csv`, `spectrum_*.csv`,
`tke_deficit.csv`, `spectral_bands.csv`, plots) with the metrics in
`v0_metrics.json`, whose `v2` block holds the V1/V2 diagnostics and whose `v0`
block holds the refinement-specific ones, plus
`spectral_bands_across_parent_nyquist.csv` and `driving_parent_vs_truth.csv`.

For the suite, in `<rundir>/analysis/`: `v0_summary.{md,csv,json}`.  The
Markdown carries the table below and the `filtered_vs_coarse` decomposition.

## Status

**Prepared, validated at tiny scale, and not submitted.**  `submit_cx3_v0.pbs`
is written and sized; submitting it is a human decision.

The whole pipeline has been driven end to end at the `v0-tiny` suite on CX3
(4 ranks, login node, 174 s): a tiny fine reference, two tiny coarse parents
(48^2 x 16 at 4 m and 24^2 x 8 at 8 m), four refined children, per-point
analysis with a cached reference accumulation, and the suite table.
`test_v0_tiny.py` is **26 tests, all passing (156 s)**.  `test_v1_tiny.py`
(10 tests, 68 s) and `test_v2_tiny.py` (26 tests, 189 s) were re-run afterwards
and are unchanged, so the shared `config.py`, `caselib.py`, `make_child_case.py`
and `analyse.py` did not regress -- the `r = 1` path through the child builder
is untouched by construction, and that is what those two runs confirm.

The tiny suite deliberately violates one of the design's own conditions:
`tiny`'s `L_rel = 8 m` is below `2 dx_P = 16 m` at `r = 4`, so the relaxation
ramp is not resolved in that parent's terms.  `config.py` prints a note saying
so and `test_v0_tiny.test_the_tiny_suite_reports_its_own_violation` asserts the
note's condition, so the violation is visible rather than quietly absorbed.  It
does not matter for a smoke test -- what is being exercised is the code path --
and the production suite clears the condition at both ratios.

**What the tiny suite showed is structural only.**  24 samples over 71 s of a
24 x 24 x 8 parent cannot support any physical claim, and its half-window
spreads (16 %) are comparable with its deficits; the numbers below are the
*shape* of the production table, not a preview of its content:

| point | r | arm | TKE `z/h>2` | of which the parent's | E resolved | marginal | sub-filter | crit A | max \|Phi\| | max divmax |
|---|---|---|---|---|---|---|---|---|---|---|
| r2-filtered | 2 | filtered | -34.93% | +5.81% | 0.636 | 0.268 | 0.563 | 0.3161 | 4.1e-15 | 8.0e-16 |
| r2-coarse | 2 | coarse | -5.40% | +34.71% | 1.361 | 1.120 | 5.509 | 1.2699 | 8.8e-14 | 9.1e-16 |
| r4-filtered | 4 | filtered | -50.35% | -14.95% | 0.623 | 0.224 | 0.154 | 0.5527 | 9.2e-15 | 7.6e-16 |
| r4-coarse | 4 | coarse | -63.43% | -52.77% | 0.640 | 0.044 | 0.023 | 1.4083 | 9.4e-14 | 8.3e-16 |

What is worth reading in it, and only this:

* **every point ran.**  `nesting_init` verified the zone was building-free
  (`nest_lparentgeom = .false.`, an assertion not a warning), all four faces
  were forced, the stored levels passed the flux check, and the cold start came
  from the interpolated parent block -- at `r = 2` and at `r = 4`, on both arms;
* **`Phi` and `divmax` stayed at round-off**, `<= 9.4e-14` and `<= 9.1e-16`, so
  refinement does not disturb the flux compatibility of design section 3.1 or
  give the projection anything extra to clean up;
* **the parent's own deficit column is doing its job.**  At `r = 2` the filtered
  parent is nearly the truth (`+5.8 %`) while the genuine 4 m LES is `+34.7 %`
  off it, and at `r = 4` they are `-14.9 %` and `-52.8 %`.  That is the ceiling
  the child is working under, separated from the child's own error, which is
  exactly what running both arms is for;
* the criterion A column fails everywhere, as it does for `test_v1_tiny` at this
  size, and means nothing at 71 s.

Nothing physical about refinement is claimed here.  That is what the production
job is for, and it has not been run.

---

# V0b -- which tangential prolongation ships

`tools/python/udprep/nesting.py`'s conservative interpolation reconstructs each
parent cell's velocity tangentially in one of two ways (`Preset.prolongation`,
`udprep.nesting.PROLONGATIONS`): `"constant"` (piecewise-constant, divergence-preserving; was the shipped
default) is exactly divergence-preserving -- design section 1.3's guarantee --
but leaves a staircase in the child's mean wind wherever the true profile is
sheared (W8's finding from V0's own filtered arm: `±0.1 u*`, alternating in
sign between adjacent child levels inside one parent cell).  `"linear"`
(conservative piecewise-linear) removes that staircase but gives up the
divergence guarantee: the reconstructed target carries a local divergence the
pressure solver then has to correct for, which V0 never measured because it
ran the old constant reconstruction throughout (nesting review 2026-09-07,
finding 2; plan section 7, R2).

**What is on trial.**  Both reconstructions, at the ratios V0 already
validated (`r = 2`, `r = 4`), driven by the C0b fine parent (960) at the FINAL
configuration -- cadence 0.5 s, `nest_timeinterp = 2` (Catmull-Rom) -- rather
than V0's own 3 s/linear-time setup.  Unlike V0 there is no `coarse` arm and no
new parent run: 960 is a genuine LES, and box-filtering its own dumps onto the
`r = 2`/`r = 4` grids (exactly V0's `filtered` arm) isolates the reconstruction
as the only free variable.  Four children: `r2-constant`, `r2-linear`,
`r4-constant`, `r4-linear` (`config.V0B`).

**What is measured**, beyond what V0 already reports for a filtered-arm point:

* pre-projection target divergence, both the driving parent's own
  (`manifest.json`: `initial_condition_divmax.parent_before_prolongation`) and
  the interpolated child's before the writer's own projection
  (`.before_projection`) -- the design section 1.3 identity holds for
  `constant` and is expected NOT to for `linear`;
* the pressure response and how far it reaches: `analyse_v0.runtime_diagnostics`
  now reports the time-mean of `|grad p|` over the zone and the interior
  separately, and of their ratio (design section 10.7 item 2 / concern C1),
  not only the per-report ratio V0 already parsed;
* the staircase amplitude (`analyse_v0.staircase_amplitude`): the
  child-minus-truth mean-flow error with each parent cell's own group mean
  subtracted, isolating the intra-cell sawtooth from the smooth background
  error criterion A already scores;
* the existing V0 machinery unchanged: criterion A, the spectral band ratios
  across the parent Nyquist, TKE deficit.

`run_v0.write_summary`'s `constant_vs_linear` block pairs the two
reconstructions at each ratio, next to the existing (empty, for this suite)
`filtered_vs_coarse` block.

**Cost.**  No parent build, spin-up or production -- 960 is read-only. Costed
from V0's own measured filtered-arm rate (box-filtering included) scaled to
960's 4800 levels (against V0's 3600): ~2.0 h total, `walltime=06:00:00`
requested (3.0x headroom). Disk ~222 GB new output before pruning
(4 x (45.3 GB nesting + 10.1 GB dumps)); `--prune-nesting` is used, so the
nesting files do not accumulate. Full derivation and a real-data smoke-test
cross-check in `submit_cx3_v0b.pbs`'s header.

## Status

**Prepared, validated at tiny scale, submitted (see job id below).**
`test_v0b_tiny.py` is **19 tests, all passing (157 s on a login node)**:
both prolongations build and run at both ratios; the divergence identity holds
for `constant` and is confirmed NOT to for `linear` (recorded, not asserted
away); the pressure-response and staircase reductions are finite and reported
for every point; the `constant_vs_linear` summary pairing is populated at both
ratios; and `prolongation=None` reproduces the writer's own default (`"linear"` since V0b settled it; the test tracks `DEFAULT_PROLONGATION` rather than naming a scheme)
exactly -- every stored variable in the nesting file, checked automatically,
not only by the one-off git-stash comparison in the implementing report.  On
real data: a single real-data smoke build of `r2-constant` against the actual
960 parent measured 926 s against the header's 1030 s estimate, and the
divergence identity held exactly (`parent divmax 9.470e-08 -> child
9.470e-08`).  `python -m unittest discover -s tools/python/tests -p
"test_nesting*.py"` (96 tests) and `python tests/validation/nesting/test_v0_tiny.py`
(26 tests) were re-run afterwards and are unchanged.

| job | experiment | submitted | job id |
|---|---|---|---|
| V0b | `v0b` (4 children: r2/r4 x constant/linear) | 2026-09-07 | `4000813.pbs-7` |

Queued in `v1_medium24` at submission (`qstat -u $USER`); check
`$EPHEMERAL/nesting-v0b/analysis/v0_summary.md` once it finishes.

---

# V0c -- is criterion A genuinely violated, or is the window too short?

The review of V0b (`nesting-review-2026-09-08-codex.md`, finding 1) found that
its headline criterion-A failures cannot be told apart from a
**matched-resolution control at the same window**: C0c's own `cr0.5` point
(expnr 970, matched 2 m resolution, same fine parent 960, same 1800 s window)
also fails (0.083 against the 0.05 bound), so the refined children's cost
relative to that control is ~17 % (`r = 2`) and ~45 % (`r = 4`), not a
pass-to-fail collapse. V1 converged, at a 10 191 s window, passes (0.0386).
The open question is whether refinement genuinely violates the bound at a
window long enough to resolve it, with its OWN control at that SAME window.

**The construction.** Three children off one fresh fine-truth parent (982),
all `arm = "filtered"`, all with `prolongation` left unset (V0b already
settled that; V0c does not sweep the reconstruction again):

* `r1` -- the matched-resolution **control**. No coarsening: `RefinedPoint
  .coarsen == 1` degenerates `CoarsenedFieldDump` to the identity, so the
  boundary data is a plain slab cut of 982's own dumps and the writer never
  reaches the prolongation. Built through the identical `RefinedPoint`/
  `DrivingParent.refined` code path as `r2`/`r4`
  (`RefinementSuite.allow_matched_control = True` lifts the normal "a
  `refine = 1` point is V1, use that instead" refusal) -- the point being
  that every number downstream comes from one analysis code path across all
  three arms, not a bespoke one for the control.
* `r2`, `r4` -- box-filtered from 982's own dumps, exactly V0b's `filtered`
  arm, at the shipped cadence (0.5 s) and interpolant (Catmull-Rom).

**The window.** The sampling floor on a temporal-mean statistic falls as
1/sqrt(N); V0b's 1800 s window (600 samples) would need ~6x the samples to
roughly halve it -- 0.098 -> ~0.049 (`r2`, right at the bound), 0.121 ->
~0.061 (`r4`, still failing but much closer), 0.083 -> ~0.042 (control,
passes). Rather than pick a new number, V0c reuses `CONVERGED.production =
10800` and `child_spinup = 600` unchanged -- a 10200 s statistics window
(~10191 s once the trailing margin is dropped), the same ~5.7x V0b's window
**and the same window V1 converged itself used**, so every V0c number is
directly comparable to the V1 converged row already in section 10.5.

**The parent (982).** Warm-started from the SAME 903 spin-up restart 960
used -- **not** a continuation of 960's own end state, which does not exist:
960's production namelist sets `trestart = 1e9`, so `writerestartfiles`
never fires during production and 960 wrote no restart of its own. 982 runs
the FULL 10800 s production from 903's restart, which reproduces 960's own
first 2400 s deterministically before continuing another 8400 s beyond where
960 stopped, and is still far cheaper than starting over: it skips the
10800 s spin-up entirely.

**The pressure diagnostic.** `nest_statint = 30 s` (not the inherited
`tstatsdump` default V0b left it at) gives ~360 `nesting_stats` reports over
the run, ~340 after the 600 s discard -- `analyse_v0.runtime_diagnostics` now
pairs each report with its own `nesting: t = ...` timestamp and reports
the zone/interior `|grad p|` norms **with the startup discarded and a spread
(mean +/- std)**, not V0b's two compulsory endpoints (review finding 2).

**Two statistics, kept separate (review finding 1).** `run_v0`'s summary
table now carries a `profile RMS [u*]` column
(`profile_metrics.u_rms_difference_over_ustar`, the RMS over height of the
area-averaged mean-flow error) right next to `crit A [u*]`
(`criterion_a.max_interior_umean_error_over_ustar`, the MAXIMUM over interior
fetch stations and faces of a single-slab RMS) -- two different reductions of
the same field, never one standing in for the other.

**Uncertainty.** The paired child-minus-truth difference's uncertainty should
come from temporal blocks of the paired differences, keeping the spatial
structure -- **not** from the parent's own half-window spread, which is not a
valid uncertainty for a paired comparison (the same correction V0b's own
prose needed, design section 10.5). This is a post-hoc analysis step on the
saved per-level differences and is not automated by this harness; do it
against the child/reference dumps once the run has produced them, before
drawing a conclusion from the raw criterion-A numbers alone.

**Disk and memory.** The dominant cost is 982's own full-domain 0.5 s field
dump (~1.08 TB): box-filtering the `r2`/`r4` boundary data needs
`CoarsenedFieldDump` over the WHOLE domain, and there is no cheap way to get
that from the much smaller `&NESTPARENT` band instead -- `caselib.NestParent`'s
box is NaN outside a thin strip sized for a *ratio-1* margin, and
`DrivingParent`'s addressing assumes `CoarsenedFieldDump` wraps a
full-domain, globally-indexed array, which a `NestParent` box is not (its
footprint already equals the child window, addressed from its own local
origin). Reworking both was judged not cheap enough to do under this task's
time budget, so V0c falls back to full `FieldDump` dumps and budgets for
them. Each child's nesting file (~204 GB, `--prune-nesting` keeps it
transient) needs ~225-230 GB of peak RSS to build (`udprep.nesting
.write_nesting_file` holds the whole array in memory) -- `submit_cx3_v0c.pbs`
requests `mem=300gb` for that. Net new/persistent disk ~1.22 TB. Full
derivation in `submit_cx3_v0c.pbs`'s header and `config.py`'s V0c comment
block.

**That memory estimate is stale -- see the "V0c16" section below.**
`make_child_case.build` (what `run_v0.py` actually calls) stopped calling
`write_nesting_file` before V0c was even written (`71011e96`, the commit
immediately before this suite's own `42367112`); it streams the slab cut
through `NestingWriter.append_level` instead, one parent level at a time, so
peak RSS does not grow with the 21600-level record. Measured directly (V0c16
section): ~650 MB, flat across a 10x range in level count. `submit_cx3_v0c
.pbs`'s `mem=300gb` still holds with a lot of room to spare -- 4004496 is not
wrong to have asked for it, just not sized from the real number -- so it is
left as submitted rather than edited retroactively.

## How to run it

```bash
qsub tests/validation/nesting/submit_cx3_v0c.pbs
```

## Status

**Prepared, validated at tiny scale, submitted (see job id below).**
`test_v0c_tiny.py` is **21 tests, all passing (~151-171 s on a login node)**:
all three points (`r1`, `r2`, `r4`) build and run; `r1`'s boundary data is a
plain slab cut and `r2`/`r4`'s are interpolated; the offline divergence
identity is correctly left unmeasured for `r1` (there is nothing to compare a
single, non-interpolated field against) and holds for `r2`/`r4`; the pressure
response is a resolved, startup-discarded series with a finite mean and
spread for every point; and profile RMS and criterion A are reported
separately and are both finite for every point. `python -m unittest discover
-s tools/python/tests -p "test_nesting*.py"` (96 tests),
`python tests/validation/nesting/test_v0_tiny.py` (26 tests) and
`python tests/validation/nesting/test_v0b_tiny.py` (19 tests) were re-run
afterwards and are unchanged (the `V0(?!b)` label-matching regex in
`test_v0_tiny.py` was widened to `V0(?![bc])` so it does not also match V0c's
own registration).

| job | experiment | submitted | job id | outcome |
|---|---|---|---|---|
| V0c | `v0c` (982 fine-truth parent, 3 children: r1/r2/r4) | 2026-09-08 | `4004496.pbs-7` | **DONE**, exit 0, 08:37:45 on `cx3-6-15` |

**Result (design section 10.5).** Criterion A is **0.0392** (r = 1 control, against V1 converged's
independently measured 0.0386), **0.0479** at r = 2 -- inside the 0.05 bound -- and 0.0746 at r = 4,
which is not. Every value roughly halves against the same case at 1800 s, as 1/sqrt(N) predicts, so
the earlier "refinement fails the mean-flow standard" reading was a short-window artefact. Resolved
TKE above z/h = 2 is -1.32 / -2.18 / -9.21 % against driving-parent deficits of 0.00 / -9.20 /
-22.82 %. The interior pressure norm is identical to 0.5 % across pre-projection divergences
spanning six orders, so the linear prolongation's extra divergence carries no measurable pressure
cost. `nest_statint = 30 s` gave 361 pressure samples per arm instead of V0b's two.

Three limits are recorded in design section 10.5 under "what V0c does not settle": all three arms
are **box-filtered** (the suite ran no genuine-coarse-parent arm), the paired differences carry
**no uncertainty estimate** (the harness's two spreads are unpaired and ~100x too large to serve),
and this reference is **much noisier aloft than V1's** (13.1 % half-window TKE spread at
z/h = 3-5 against V1's 1.1-1.5 %). It also **failed V7**: the children read boundary data for
48.2 / 45.7 / 49.2 % of their runtime against a < 1 % criterion.

Full table: `$EPHEMERAL/nesting-v0c/analysis/v0_summary.md` (dumps retained, ~1.1 TB, so the
paired-uncertainty re-analysis needs no new run -- but `$EPHEMERAL` is purged periodically).

---

# V0c16 -- V0c's own suite, decomposed for CX3's 16-core queues

V0c (job 4004496, `submit_cx3_v0c.pbs`, `ncpus=64:mem=300gb`) sat queued in
`v1_medium24` for six hours: `qstat` reported "Insufficient amount of
resource: ncpus" against 281 queued / 18 running in that queue. PBS routes a
job by `ncpus` and walltime band only -- a shorter walltime at `ncpus=64`
keeps a job in the same `v1_medium24/72` pair, it does not move it to a
smaller, faster-clearing queue. The 16-core queues turn over far better
(`qstat -Qf` snapshot: `v1_small72` ~211 queued / ~118 running,
`v1_small24` ~328 queued / ~71 running -- `v1_small72`'s running:queued ratio
is ~2.6x `v1_small24`'s). V0c16 is that 16-core twin: the SAME three-point
suite (matched-resolution control `r1` plus box-filtered `r2`/`r4`, same
128 x 128 x 64 child, same window), run through a SEPARATE preset/suite
(`config.V0C16_FINE`/`V0C16`, fresh experiment numbers 992-997) and a
SEPARATE `$EPHEMERAL` directory, so 4004496 is untouched and keeps working if
it ever starts. Whichever job starts first wins; once one is producing
output, `qdel` the other rather than let both run to completion.

**The memory claim on V0c is stale.** `config.py`'s V0c comment block and
`submit_cx3_v0c.pbs` both size `mem=300gb` from "`udprep.nesting
.write_nesting_file` holds the whole nesting file in RAM" -- ~204 GB for the
full 21600-level record, times a 1.11x file-to-RSS ratio measured on a much
smaller, pre-streaming build (964, 35 GB). That description no longer matches
the code: `make_child_case.build` -- what `run_v0.py` actually calls for
every V0c/V0c16 point -- stopped calling `write_nesting_file` when the slab
cut was streamed through `NestingWriter.append_level` (`71011e96`, "stream
the slab cut through NestingWriter (W4 adoption)", the commit immediately
before V0c's own `42367112` in history). That commit's own docstring says so
directly: "the slab cut holds one parent level rather than the whole record:
its peak memory does not grow with the number of levels" -- and its own
`tracemalloc` numbers already showed this flat (103.1 MB at 40 levels,
103.0 MB at 14).

Verified directly on this worktree rather than trusted from the commit
message: wrapped `make_child_case.build` for `v0c-tiny`'s `r1` point in
`/usr/bin/time -v` (whole-process RSS) plus `tracemalloc` + `resource
.getrusage` (in-process), against the real 240-level `v0c-fine-tiny` parent
-- **654.7 MB peak RSS / 111.2 MB tracemalloc peak**. Re-ran the identical
build against a copy of that parent with its `fielddump` netCDF time axis
artificially extended 10x (2400 levels, `u`/`v`/`w` duplicated) -- **657.5 MB
/ 111.3 MB**, flat to within 0.4%. `caselib.FieldDump.read_level` also only
ever holds one global `(itot, jtot, ktot)` parent level at a time (freed
every loop iteration, ~96 MB of `float64` for the three components at the
full 256 x 256 x 64 parent domain), so nothing in this pipeline scales with
the 21600-level production record; the true floor is imports/geometry/IBM
arrays (the tiny case's measured ~650 MB) plus one such slab, of order 1 GB
at production scale, not 225-230 GB. `submit_cx3_v0c.pbs`'s `mem=300gb` is
not *wrong* for that reason -- it holds the real number with room to spare --
just not sized from it, so it is left as submitted rather than edited
retroactively.

`mem=120gb` below is therefore not sized from the python build step at all
(it needs a tiny fraction of that): it is sized for headroom against the
SOLVER's own memory, unmeasured directly here. V0b shipped `mem=128gb` for
the identical 128 x 128 x 64 child geometry, and the 64-rank `converged`/V0c
parent (256 x 256 x 64, `nfcts = 12992`) already runs at `mem=128gb` in the
V1 converged job (3991175) -- 16 ranks holding the same global domain should
not need materially more. Since queue selection is by `ncpus`/walltime band
only, not `mem` (above), asking for the `v1_small*` queues' near-full
headroom costs nothing.

**The decomposition.** 4x4 = 16 ranks for both the parent (`itot=jtot=256`,
64 cells/rank) and the child (`itot=jtot=128`, 32 cells/rank), not V0c's 8x8.
`Preset.validate` checks divisibility; both domains clear it. The r2/r4
driving-parent placeholders are bookkeeping only (`RefinedPoint.runs_driver`
is `False` for every 'filtered' point) so their own `nprocx`/`nprocy` is set
to 4x4 too, uniformly, rather than mirroring V0c's own asymmetric
`(8,8)`/`(4,4)` split.

**Walltime -- measured, and chosen for the queue, not only for margin.**
Measured on this worktree (login node, `build/release/u-dales`), reusing
already-built production-scale case files (`$EPHEMERAL/nesting-v1-
converged/903` and `.../904`) with the decomposition patched to 4x4 and a
short cold-start `runtime`, the SAME domains/`nfcts` this job actually runs:

| case | domain, nfcts | 30 s sim, 3 (parent) / 2 (child) probes | sim/wall ratio |
|---|---|---|---|
| parent | 256x256x64, 12992 | main-loop wall 27.72 / 18.32 / 26.37 s | 1.08 / 1.65 / 1.15 (mean 1.30, median 1.15) |
| child (nested) | 128x128x64, 2240 | main-loop wall 10.78 / 10.17 s | 2.81 / 2.97 (mean 2.89) |

The spread reflects login-node contention (load average ~8.5 against 64
cores at the time), not a real swing in solver cost -- the same case gave
bit-identical `dt` sequences at every rank count tried, since the
decomposition does not change the physics. Consistency check: composing the
measured 4->16-rank speedup with the measured 16->64 speedup reproduces
`clusters.md`'s independently-known 4->64 speedup (3.65x) to within
rounding, so these probes are not an outlier. Using the SLOWER ratio from
each pair: parent production `10800 / 1.15 = 9391 s = 2.61 h`; three child
runs `3 x (10800 / 2.81 = 3844 s) = 192 min`; children slab-cut/box-filter
reads and analysis stay at V0c's own rank-independent python-side figures
(~246 min and ~65 min). Total **~11.4 h**.

11.4 h is under 24 h -- a walltime sized purely for margin over that number
would stay in `v1_small24`, not move to `v1_small72`. `walltime=30:00:00` is
requested instead specifically so this job lands in the better-clearing
72-hour queue band, with ~2.6x headroom over the 11.4 h estimate folded in
anyway for the measured run-to-run noise and because this is the first real
run at this exact decomposition. Full derivation in
`submit_cx3_v0c16.pbs`'s header.

**Disk.** Same domains as V0c, so the same volumes: ~1.08 TB parent dumps +
3 x ~45 GB child field dumps = ~1.22 TB persistent, plus one ~204 GB nesting
file transient per child (`--prune-nesting`). A SEPARATE run directory
(`$EPHEMERAL/nesting-v0c16`) keeps this from colliding with 4004496's
(`nesting-v0c`) if both happen to run at once.

## How to run it

```bash
qsub tests/validation/nesting/submit_cx3_v0c16.pbs
```

## Status

**Prepared, validated at tiny scale (real 4x4/16-rank decomposition),
submitted (see job id below).** `test_v0c16_tiny.py` is **14 tests, all
passing (~115 s on a login node)**: all three points (`r1`, `r2`, `r4`) build
and run at 4x4; `r1`'s boundary data is a plain slab cut and `r2`/`r4`'s are
interpolated; the nesting file is flux-balanced at every ratio; profile RMS
and criterion A are both finite; the summary table is written.
`python -m unittest discover -s tools/python/tests -p "test_nesting*.py"`
(96 tests) and `python tests/validation/nesting/test_v0c_tiny.py` (21 tests,
its registration check widened to `V0c(?!\d)` so it does not also match
V0c16's own label) were re-run afterwards and are unchanged.

| job | experiment | submitted | job id | outcome |
|---|---|---|---|---|
| V0c16 | `v0c16` (992 fine-truth parent, 4x4/16 ranks, 3 children: r1/r2/r4) | 2026-09-08 | `4008364.pbs-7` | **superseded -- never ran** |

**Outcome: V0c won the race.** 4004496 started in `v1_medium24` at 03:19 on 2026-09-09 and finished
at 11:57 with exit 0, so its 16-core twin was cancelled at 06:49 without ever starting
(`qstat -x` records "Not Running: Insufficient amount of resource: ncpus and terminated"), exactly
as the twin-job rule above prescribes. No `$EPHEMERAL/nesting-v0c16` directory was created and no
V0c16 analysis exists. The preset, suite and submission script are kept: they are the working
16-core decomposition of the same three-point suite, and the next time a 64-core job will not clear
this is the one to submit.

---

# V3 and V4 -- when parent and child geometry differ

`docs/udales-nesting-design.md` section 9.4 is about the configuration the
scheme actually exists to support: a parent that resolves *different* buildings
from the child, or none at all.  Section 10.4 turns that into two rows.

> **V3 -- parent without buildings.**  The parent resolves no geometry; the
> child has buildings starting **at** the inner zone edge.  Deliverable: the
> adjustment length, and confirmation that a standoff lengthens rather than
> shortens it.
>
> **V4 -- different parent geometry.**  A parent with a different building
> layout.  Deliverable: interior statistics, confirming that the interior is
> insensitive to the mismatch beyond the adjustment fetch.

Both are one-way and at refinement ratio 1, like V1 and V2, so nothing here
tests the interpolation either (V0 does).

**Two wordings in that table are worth flagging rather than inheriting**, and
this directory does not own the design document, so they are recorded here.  V3's
deliverable as written -- "*confirmation that* a standoff lengthens rather than
shortens it" -- presupposes its own answer; what is built below is a test that
can return `REFUTED`, and if it does, the design document is what has to change.
V4's deliverable says "interior statistics vs S3", and no S3 appears anywhere
else in the document; it is read here as V1, which is the only matched-geometry
reference that exists.

## What section 9.4 claims, and how these experiments can refute it

Section 9.4 makes two falsifiable claims.

1. **An internal boundary layer must develop** between the inner edge of the
   zone and the first building row, because the imposed near-surface profile is
   in equilibrium with the parent's roughness, not the child's canopy.  The
   adjustment length is what V3 measures.
2. **The buildings should start immediately at the inner zone edge, not after a
   standoff.**  A building-free standoff is claimed to be *actively
   counterproductive*: the flow there adjusts only to the ground roughness
   through weak shear-driven mixing and then has to adjust a second time on
   reaching the canopy, so `W = 0` gives one adjustment instead of two.

Claim 2 is the one that can fail, and V3 is arranged so that it can.  Turned
into predictions that the harness evaluates without any tolerance chosen after
the fact:

> **P-a**  The adjustment length measured **from the inner edge of the zone**
> does not *decrease* as the standoff grows.  (That is the domain a layout has
> to spend, which is what the claim is about.)
>
> **P-b**  At a fixed station measured from the zone edge, the canopy of the
> 0-cell child is at least as close to equilibrium as every other child's.
>
> **P-c**  The adjustment length measured **from the first building face** does
> not *decrease* as the standoff grows.

**P-c is the sharp one, and P-a is nearly free -- worth being blunt about.**  A
standoff moves the canopy downstream, so the canopy's adjustment trivially
finishes later measured from the zone edge, by at least the standoff length,
whatever the physics.  The non-trivial content of "one adjustment instead of
two" is that the *second* adjustment is slower than the single one would have
been, because the flow reaching the canopy has already equilibrated with the
ground and has to be reworked.  That is `adjustment_from_first_row`.  If it comes
out *shorter* behind a standoff -- if a decelerated approach flow makes the
canopy adjustment quicker -- then a standoff is not counterproductive in the way
section 9.4 claims; it merely costs its own length, which is a much weaker
statement than the design makes, and this README should then say so.

A refutation is a standoff that reaches equilibrium sooner from the zone edge or
from the first building face, or that is closer to equilibrium at a fixed station
by more than the sampling spread.  `run_geometry.py` reports `design_9_4_standoff_claim` as `SUPPORTED`,
`REFUTED` or `INCONCLUSIVE`.  `INCONCLUSIVE` is a real outcome and is emitted
whenever the children were not distinguishable at all -- the residuals at every
fixed station lying within each other's sampling spread, and fewer than two
children reaching equilibrium.  It is **not** agreement, and `SUPPORTED` is
withheld unless the comparison could have shown otherwise.

**The standoffs are 0, 5, 15 and 40 cells.**  The first three are the design
table's own comparison.  40 cells = 80 m = 5 h was added because 0, 5 and 15
cells are 0, 0.6 and 1.9 building heights -- all short compared with any
plausible adjustment scale -- so a sweep of only those could fail to separate
the hypotheses for want of lever rather than because the claim is right.  Giving
the claimed mechanism room to act makes the test sharper, not kinder.

## The parent sub-region is **not** a reference here, and that is the whole design problem

`analyse.py` answers one question -- does the child reproduce the parent
sub-region it was cut from? -- and every number in it, criterion A included, is a
child-minus-parent difference over cells that are fluid in **both** runs.  That
is right for V1 and V2, where the child is a sub-model of its parent.  It is
wrong here, and using it anyway would be the easiest way to produce a confident
and meaningless answer:

* in **V3** the parent has no buildings at all, so `<u>_child - <u>_parent`
  measures the canopy rather than the nesting.  Criterion A will be far outside
  its bound and that is not a failure of anything;
* in **V4** the parent's buildings are in *different places*, so the
  fluid-in-both mask discards roughly a quarter of the child's canopy layer in a
  pattern set by the parent's array -- an averaging domain with nothing to do
  with the one V1 averaged over.

So the child-versus-parent block is still computed for every child and is still
written out -- it is asked for as a comparable diagnostic, and for V4 it is
genuinely interesting as a measure of how far the interior is *free* to differ
-- but it is written as `child_vs_parent_metrics.json`, not `v1_metrics.json`,
and **no pass criterion is applied to it.**  The references are elsewhere:

| | reference for the interior | why |
|---|---|---|
| **V3** | a periodic run of the child's own canopy at the child's own forcing (`920`) | that is what "in equilibrium with its own canopy" *means*; without it an adjustment length can only be self-referential, which cannot tell an adjusted canopy from one that stopped adjusting short of equilibrium |
| **V4** | the **V1 `converged` child** | identical geometry, size, zone, forcing and schedule; the only difference is the layout its parent resolved |

Both references are reduced by the same code, over each child's **own** fluid
mask, in `analyse_geometry.py`.  For V4 that means the mismatched child and the
V1 child are averaged over the same cells, which they are not if either is
intersected with its parent's mask.  `run_geometry.py` aborts if the two
children's solid masks are not identical.

**Where the V1/V2 diagnostics are.**  `analysis/<child>/child_vs_parent_metrics.json`
carries exactly what `v1_metrics.json` carries for a V1 or V2 child -- interior
mean and resolved-TKE profiles with the parent's half-window noise floor,
error-versus-distance curves and decay lengths for `<u>` and TKE from all four
faces, streamwise spectra and band ratios at three heights, the `tke_series`
equilibration trace, the `v2` block (TKE deficit, TKE error versus fetch,
criterion A, the common-block reduction) and the `solid_mask` census -- together
with the CSVs of all of it.  It is written for every child of both experiments,
so the two campaigns are comparable with V1 and V2 line for line.  Read it as
*context*: for V3 and V4 the parent is the forcing, not the truth, and the
criterion A line in it will say FAIL for reasons that have nothing to do with
the scheme.  `--no-legacy-metrics` skips it.  V4's own criterion A -- against the
V1 child -- is `criterion_a_prime` in `v4_metrics.json`, and *that* one is a
criterion.

## V3 -- parent without buildings

```
reference     920   96 x  96 x 64,  192 x 192 x 128 m, periodic 6 x 6 cube array
                    fixed dpdx = 1.25e-3 m/s^2 -> u* = 0.4 m/s (V1's forcing)
                    spin-up 3600 s, production 3600 s, dumps every 3 s
parent        921  320 x 160 x 64,  640 x 320 x 128 m, FLAT -- no buildings
                    volume-flow forced at the reference's measured bulk velocity
                    spin-up 9000 s, production 3600 s
parent-cubes  926  320 x 160 x 64, V1-style ALIGNED canopy where the child's
                    zone sits, fixed dpdx (same forcing as everything else)
                    spin-up 3600 s, production 3600 s
children   922-925  256 x 128 x 64, 512 x 256 x 128 m, origin (64, 32) m, on 'parent'
                 canopy of 16 m cubes on a 32 m period, starting 0 / 5 / 15 / 40
                 cells past the clear box; 14 / 14 / 13 / 11 streamwise rows
                 zone L_imp = 6 m (3 cells), L_rel = 18 m (9 cells), tau = 1 s
                 streamwise interior 464 m = 29 h; spanwise 208 m = 13 h
child        927  same window and zone as 922 (standoff 0), on 'parent-cubes':
                 the parent resolved buildings where the child's zone now sits
                 and the child clears them (nest_lparentgeom = .false.) --
                 'cleared-parent-cubes', compared against 922
```

Everything the child's canopy is made of -- 16 m cubes, 16 m streets, `h = 16`
m, `dx = 2` m, `u* = 0.4` m/s, a 128 m rigid lid -- is V1's, so the canopy V3
measures is the canopy V1 and V2 measured.

**Boundary cadence, interpolant and source.**  Both parents write their
child's boundary at 0.5 s (`&NESTPARENT`, `parent_output = "both"`,
`config.Preset.parent_output` via `presets_geometry._NESTPARENT_CADENCE`): with
the mean wind at the domain top ~5.5 m/s (V1's converged run, `z/h = 2`) and
`dx = 2` m, the operating rule `C_dump <= 2` (design section 10.5) needs
`dtdump <= 0.73` s, so `C_dump = 1.4` at 5.5 m/s (0.75 at the preset's own
`u0 = 3` m/s -- the smaller, less conservative number `summary()` prints).
Every child interpolates with the unlimited Catmull-Rom cubic
(`nest_timeinterp = 2`), C0c's recommended default, rather than V1's linear.
Each parent's own `&OUTPUT` full-domain field dump stays at the old 3 s
cadence -- `periodic-stats` reads it back for the bulk velocity and canopy
statistics `run_geometry.py` needs, and 3 s was always enough for those --
so the fine cadence costs only the `&NESTPARENT` band, not a 6x-larger
full-domain dump.  `make_child_case`'s manifest records
`driving_source = "nestparent"` for every V3 child, and `Preset.summary()`
prints the `C_dump` and interpolant lines for every preset, parent or child.

**The cleared-parent-cubes arm** (nesting-plan-2026-09-06.md section 0, "New
from V2").  V2 found that a child which clears cubes out of its *own* zone
loses 0.05-0.07 `u*` of canopy-layer mean flow, because the wakes those cubes
would have shed inside the zone are missing from the inflow -- V3's own
question, seen from the other side of the boundary.  `standoff0` (parent has
no buildings anywhere) and `cleared-parent-cubes` (parent has V1's own
aligned buildings filling the same domain, including where the child's zone
sits) share the **identical** child construction: both use
`child_phase = "standoff"` at `standoff_cells = 0`
(`presets_geometry._v3_child`), which places the child's own lattice from the
clear box outward and never reads the parent's layout at all, so the two
children's cube centres are asserted equal cube for cube
(`test_v34_tiny.py`).  Only the periodic run driving the boundary changes --
flat (`V3_PARENT`) versus V1's own aligned canopy (`V3_PARENT_CUBES`) -- so
the comparison isolates what a genuinely absent boundary condition costs the
adjustment length from what "parent had cubes there" costs.  (An earlier
version of this arm used `child_phase = "parent"`, which regenerates the
*global* lattice and clips it to the clear box; that lattice is anchored to
the parent's coordinate origin rather than to the clear box the way
`"standoff"` is, so it placed an entirely different, disjoint set of cubes --
caught by comparing the two child layouts directly rather than assuming the
mechanism did what its name suggested.)  Both parents are driven by the same
fixed `dpdx`: a real canopy carries its own drag, unlike the flat parent, so
`parent-cubes` needs no volume-flow calibration against the reference.

### Why the flat parent is flow-rate forced, and what that costs

This is the one design decision in V3 that is not forced by the brief, so it is
worth stating in full.

A flat surface at the cases' ground roughness (`factypes` type 1, `z0 = 0.05` m
-- the default the preprocessing writes, and the same one the cubes carry)
driven by the canopy's own `dpdx` runs at roughly **twice** the canopy's
equilibrium bulk velocity, because the canopy drag it is missing is most of the
drag.  That is a *bulk momentum* mismatch, and it is not the one section 9.4 is
about: the child's interior would then be decelerating everywhere, at a rate
that puts full adjustment hundreds of building heights downstream, and there
would be no near-surface equilibrium for anything to adjust *to*.  V3 would
measure a bulk imbalance and report "not reached".

So the flat parent is driven by `luvolflowr` at the reference run's **measured**
bulk velocity instead.  The child keeps V1's fixed `dpdx`.  The child's interior
is then in global momentum balance -- its own forcing matches its own canopy's
drag at that bulk speed -- and the only thing out of equilibrium is the *shape*
of the near-surface profile: a log layer down to the ground where the child has a
canopy.  That is precisely the imposed profile section 9.4 says will be out of
equilibrium, and the adjustment length becomes a property of the internal
boundary layer rather than of a bulk imbalance.

The bulk velocity is a **measurement, not a guess**: `run_geometry.py` runs the
reference first, reduces it, and writes its fluid-masked depth-averaged `<u>`
into the parent's `uflowrate` -- computed the same way `modforces.f90:404-410`
computes it, or the parent would be asked to hold a different number from the one
it reports.  The preset carries a nominal value only as a sanity bound: a
measured bulk more than a factor of two from it aborts, because that would mean
the reference run is broken rather than that the guess was poor.  The tiny
experiment measured 2.5493 m/s in the reference and the flat parent then held
2.5434 m/s, 0.2 % low.

**What this costs, stated plainly.**  A flat wall at `z0 = 0.05` m carrying the
canopy's bulk velocity carries a *weaker* surface stress than the canopy does --
`u*` about 0.18 m/s against 0.4 -- because a smooth surface needs less stress to
hold the same wind.  So the V3 parent also delivers less turbulence to the
boundary than a matched parent would, and the child has to regenerate the
difference as well as reshape the profile.  Both mismatches are measured and
reported (`mismatch` in `v3_metrics.json`: the bulk, canopy-velocity and
roof-level stress differences between the imposed and the equilibrium states),
so the adjustment length is quoted *for a stated mismatch* rather than as a
universal number.

There is no third option available at this resolution.  Matching both the bulk
and the stress would need a flat wall with a city-scale roughness length,
`z0 ~ 2 m` by Macdonald's morphometry for this array -- and the IBM wall function
needs `log(dist/z0) > 1` at the first cell centre (`modibm.f90:1383`), which on a
2 m grid caps `z0` at about 0.2 m.  A genuinely coarse parent could carry it;
that is V0's territory, not V3's.

### What is measured

Everything is reduced onto **streamwise blocks one cube period wide**, phase
locked to the child's canopy lattice and continued upstream into the
building-free standoff, clipped to the interior.  Two reasons and both matter:
`<u>(x)` inside a cube array swings by tens of per cent within one period, so a
block must contain exactly one cube and the same part of the pattern or the
streamwise signal is aliasing against the array; and continuing the same phase
upstream puts the standoff region on the same abscissa as the canopy, which is
what makes a 0-cell and a 40-cell standoff comparable at a fixed station.

| quantity | where | what it is |
|---|---|---|
| `u_canopy` | `blocks.*[].u_canopy` | mean `<u>` below roof height, over the fluid cells of the block |
| `uw_at_roof` | `blocks.*[].uw_at_roof` | resolved `<u'w'>` at the first level above the roofs -- the sharpest streamwise signature of an internal boundary layer, and the direct proxy for the "visibly wrong facet stresses" section 9.4 expects on the first few rows |
| **adjustment length** | `adjustment` | the smallest fetch beyond which **every** later row is within 5 % of the reference, for both quantities, against **two** references and reported from **two** origins |
| IBL depth | `ibl` | the lowest height above which the block's `<u>(z)` is within 5 % of the imposed bulk velocity of the parent's own profile, at that height and every height above |
| residual at fixed stations | `fixed_stations` | the relative canopy-velocity error at 5 h, 10 h and 20 h from the inner zone edge -- the cross-standoff comparison that needs no adjustment length to exist |
| the imposed-vs-equilibrium mismatch | `mismatch` | how far the parent's state is from the canopy's equilibrium: the size of what the child has to work off |

The 5 % tolerance is design section 0's own interior bound, reused; it was fixed
before any number was computed and is not a knob.  "Adjusted" uses the same "and
stays there" rule `analyse.decay_length` uses, so an adjustment length here means
what a decay length means there.  A settle that happens only at the very last
row is flagged (`settled_only_at_the_last_row`) and not reported as an
adjustment length -- there is no row after it to disagree.

**Two references, side by side, neither chosen after the fact.**
`vs_equilibrium` is against the periodic reference: the absolute answer, and the
one section 9.4 asks for.  `vs_last_row` is against the child's own last row: the
self-referential answer, which is all that exists if the flow never reaches
equilibrium within the domain, and which is reported *so that the two can be
told apart*.  `None` -- never reached within the available fetch -- is a result.

**Two origins.**  From the **inner edge of the zone**, which is the domain a
layout has to spend and the quantity P-a is about; and from the **first building
face**, which says whether the canopy adjustment itself is faster or slower
behind a standoff.  A standoff can perfectly well shorten the second while
lengthening the first, and section 9.4's claim is about the first.

**The spanwise scope.**  The child's spanwise faces are nested too, so lateral
internal boundary layers spread inward from both of them, further with fetch.
The headline statistics are taken over the central half of the canopy's spanwise
extent (`y_core_fraction = 0.5`); the full width is emitted alongside as
`blocks.canopy` and `adjustment_full_width`, so the contamination can be seen
rather than assumed away.

### V3's zone contains no soft obstacles at all

The child's zone is building-free, and the field it is relaxed towards -- a flat
parent's -- contains no buildings either.  So unlike V2 (and unlike V4 below),
V3's ramp carries no low-velocity imprint of a building the child does not
resolve.  It is the cleanest configuration in the whole campaign: the only thing
crossing the boundary is a horizontally homogeneous flat-wall flow, which is
exactly what makes the adjustment behind it attributable to the canopy.

## V4 -- different parent geometry

```
parent    930  256 x 256 x 64 -- V1's parent, cube for cube and second for
               second, except that its array is STAGGERED rather than aligned
               (udgeom.create_cubes 'SC'; same 16 m cubes, same 32 m period,
               same plan area density, so the same solid cell count)
child     931  128 x 128 x 64 at origin (128, 128) m carrying V1's ALIGNED array
               V1's zone (3 + 9 cells), tau, nest_nwall, nest_linitfromparent,
               forcing and schedule, to the digit -- but NOT V1's boundary
               cadence or interpolant, see below
baseline       the V1 'converged' child on disk, read and reduced, not re-run
```

**Boundary cadence, interpolant and source -- and the one thing that is no
longer "to the digit".**  The parent writes its child's boundary at 0.5 s
(`&NESTPARENT`, `parent_output = "both"`, `presets_geometry._NESTPARENT_CADENCE`):
`C_dump <= 2` (design section 10.5) needs `dtdump <= 0.73` s at the domain
top's ~5.5 m/s, and 0.5 s gives `C_dump = 1.4` there.  The child interpolates
with the unlimited Catmull-Rom cubic (`nest_timeinterp = 2`), C0c's
recommended default, rather than V1's linear.  The parent's own `&OUTPUT`
full-domain dump stays at the old 3 s cadence for `periodic-stats`, so the
fine boundary costs only the `&NESTPARENT` band.  `manifest.json` records
`driving_source = "nestparent"`.  **This is the one respect in which the child
is deliberately not V1's, to the digit**: the V1 `converged` baseline this
child is diffed against ran at V1's original 3 s cadence and linear
interpolant, so `criterion_a_prime` and every other number in `v4_metrics.json`
now compares a mismatched-parent child with a *better* boundary treatment
against a matched-parent baseline with the *old* one.  Per C0's own numbers
(design section 10.5) that is not a small effect at these heights -- the
resolved-TKE deficit above the canopy falls from about 11 % at 3 s/linear to
2 % at 0.5 s/cubic on a matched parent -- so a measurable difference from V1
in this comparison is not evidence of a geometry-mismatch cost until the
cadence/interpolant confound is ruled out, and this README says so rather than
letting the number speak for itself.  Re-running the V1 baseline at the same
cadence and interpolant would remove the confound; it has not been done here
because the campaign's verdict (nesting-plan-2026-09-06.md, C0) is that the
fine cadence is the correct boundary treatment going forward and the old V1
number is itself superseded, not a fixed reference to be matched.

**The child is V1's child.**  Not merely the same size -- the same cubes in the
same places.  V1 got its building-free zone by carving a plaza out of the parent;
V4 gets it by clearing the child, which section 9.4 (and the `clear_child_zone`
work in `config.py`) says is the right way round.  Those are different mechanisms
and they have to produce the same layout, so `test_v34_tiny.py` asserts it cube
for cube against `config.CONVERGED`, and asserts that both drop the same 28
cubes.  Every other parameter is asserted equal too.  So the only thing that
differs between this child and V1's is **the layout its parent resolved**, and
the answer to "does a mismatched parent layout cost anything in the interior" is
a difference from V1's numbers.

### What is measured

| quantity | what it is |
|---|---|
| `tke_difference.above` / `.canopy` | `(TKE_mismatch - TKE_V1)/TKE_V1` above `z/h = 2` and inside the canopy, against the two runs' half-window spreads combined in quadrature |
| `umean_difference` | the `<u>(z)` difference, rms over the column in `u*`, with its own sampling floor |
| **criterion A'** | design section 0's interior bound taken against the *right* reference: `max_interior rms\|<u>_mismatch - <u>_V1\|/u*` against 0.05.  Not against the parent, whose buildings are elsewhere |
| `tke_error_vs_distance` | the resolved-TKE difference between the two children, per slab, against distance from each lateral face, with the baseline's own half-window floor.  V1's criterion B with the right reference: if a mismatched parent layout leaves a signature it should be largest near the boundary, where the imposed field carries the parent's wakes in the wrong places, and decay inward.  A flat curve says the difference is not coming from the boundary at all |
| `spectra_child_over_baseline` | the child/child spectral ratio in the 16-64 m, 8-16 m, `> L/4` and `< 4 dx` bands at each of V1's three heights |

The half-window spreads are the point.  V4's child and V1's were driven by
*different realisations* of the turbulence, so the difference between them
contains weather as well as mismatch.  That is bounded, not eliminated: every
difference is reported next to the combined spread, and a difference smaller
than its spread means **V4 has not measured a cost** -- which is exactly what
section 10.4 expects ("the interior is insensitive to the mismatch beyond the
adjustment fetch").  `verdict.<child>.measurable` is that comparison, made on the
median per-height significance, as V2's is.

### The soft obstacle, for V4 specifically

The V2 README sets this out in general: inside the child's cleared band the
child is relaxed towards the parent's velocity field, and that field contains the
parent's cubes -- as near-zero velocity where a cube stands, and as wakes
downstream of one.  So the ramp carries a low-velocity imprint of a building the
child does not itself resolve: no IBM enforcing it, no wall stress, no ongoing
production.  It is neither a building nor a plaza.  Three statements specific to
V4:

* **It does not touch the mass budget.**  `Phi` is evaluated on the child's
  boundary faces, and both sides of the comparison are unmasked there -- the
  child has no solid cells in its band, and the Python writer corrects whatever
  field is stored over all four faces.  The correction drives the stored `Phi` to
  round-off and the runtime diagnostic stays there, cube in the plane or not.
  `test_v34_tiny.py` checks `Phi` and `divmax` at round-off on every child.
* **It does not change the boundary treatment.**  `nest_lparentgeom`, the
  weights, `tau`, the shape function, the guard width and `nest_nwall` are V1's,
  unchanged.  What differs is *what the imposed field describes*, not how it is
  imposed.
* **It is not separable from the effect being measured, and V4 does not pretend
  otherwise.**  Here V4 differs from V2 in a way worth spelling out.  In V2 the
  parent's array was aligned with the child window and no cube came closer than
  8 m to a child face, so the guard strip was over open ground in the parent too.
  A **staggered** array's displaced columns put cube centres on the child's
  spanwise faces: at four streamwise stations, the parent has a cube straddling
  the child's south face and four more straddling its north face, reaching 8 m
  into a 24 m band.  The imposed inflow there is near zero over 8 x 16 m patches
  below roof height.  Those patches sit inside the band, their wakes advect
  downstream along the boundary rather than into the interior, and the analysis
  interior begins 24 m in -- so they do not enter the region compared.  But they
  *are* part of what "a mismatched parent layout" does, and V4's answer
  therefore bundles two things: the parent's interior wakes arriving in the wrong
  places, and the parent's buildings intersecting the child's boundary planes.
  If V4 finds a measurable interior cost, separating those two would need a
  further experiment -- a child window offset so that no parent cube straddles a
  face -- and this README should be the place that says so rather than the place
  that quietly implies it was already done.

### Why there is no matched-geometry control on this parent

The obvious control is the same child carrying the parent's *staggered* array --
V1 repeated on this realisation.  **It cannot have a building-free zone at any
size, and that is a theorem rather than a budget.**  Write `c` for a cube centre
in child metres and `L` for the child extent.  A cube is dropped from the clear
box when `c < 34` or `c > L - 34` (26 m of clearance plus 8 m of half width), and
it also reaches the analysis interior when `16 < c < L - 16`.  So the residues
modulo the 32 m period that a dropped cube may not occupy span two 18 m windows
-- 36 m of a 32 m period -- and a staggered array's two column families sit
exactly half a period apart, so one of them always lands in a blocked window
whatever the child's origin or size.  Clearing the zone would remove buildings
from the region the statistics are taken over, which
`removed_cubes_reaching_the_interior` refuses and should refuse.

The alternatives are to keep the parent's cubes inside the zone and run
`nest_lparentgeom = .true.` -- legal for self-nesting, but then the control no
longer shares V1's boundary treatment and stops being a control -- or to widen
the zone, which changes the variable under test.  `test_v34_tiny.py` turns the
argument into a checked invariant by asserting that `validate()` refuses such a
preset.  If it ever stops firing, the argument needs revisiting, not deleting.

## Layout

| File | What it is |
|---|---|
| `presets_geometry.py` | **every** parameter of V3 and V4, as `GeoPreset` objects whose parent and child layouts are independent, grouped into `Experiment` objects.  Imports the shared machinery from `config.py` and registers its presets there; adds nothing to `config.py` itself.  Run it to print both experiments |
| `make_geometry_cases.py` | the periodic case builder: flat ground, a staggered array, and volume-flow-rate forcing.  The namelist itself is still `make_parent_case.parent_sections` -- three keys are patched, not restated.  The nested children need nothing new: `make_child_case.build` already takes the child's layout from `Preset.child_cube_centres` |
| `analyse_geometry.py` | the V3/V4 measurements: streamwise blocks, adjustment length, IBL depth, and the child-against-child comparison.  Carries the `u'w'` cross-moment, which `analyse.py` does not |
| `run_geometry.py` | end-to-end driver for both, stage by stage, plus the cross-child summary and the verdicts |
| `test_v34_tiny.py` | both experiments as a unittest at login-node size, plus the production configuration checks that need no run |
| `submit_cx3_v3.pbs`, `submit_cx3_v4.pbs` | the two production jobs.  **Review before submitting.** |

## How to run it

Smoke test (about 3 minutes on a login node, 4 ranks -- two tiny experiments,
nine solver runs):

```bash
module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
source /rds/general/user/mvr/home/udales/.venv/bin/activate
export TMPDIR=$EPHEMERAL
python tests/validation/nesting/test_v34_tiny.py
```

or one experiment on its own, which prints the full tables:

```bash
python tests/validation/nesting/run_geometry.py $EPHEMERAL/v3-tiny --experiment v3-tiny
python tests/validation/nesting/run_geometry.py $EPHEMERAL/v4-tiny --experiment v4-tiny
```

Production (**do not run this on a login node**):

```bash
qsub tests/validation/nesting/submit_cx3_v3.pbs
qsub tests/validation/nesting/submit_cx3_v4.pbs
```

V3 needs nothing from V1 -- it runs its own reference and its own parent.  V4
reads the V1 `converged` child at `$EPHEMERAL/nesting-v1-converged/904` as its
baseline and refuses to start without it.

Stages are `periodic-case`, `periodic-spinup`, `periodic-production`,
`periodic-stats`, `child-case`, `child`, `analysis`, `summary`; `--only`
restricts the children acted on, so a single failed child can be redone:

```bash
python tests/validation/nesting/run_geometry.py $EPHEMERAL/nesting-v3 \
    --experiment v3 --only standoff40 --start-at child-case --yes
```

### What comes out

```
<rundir>/analysis/
  periodic_<key>.json          the equilibrium reference / the imposed state:
                               profiles, canopy-layer reductions, bulk velocity
  <child>/v3_metrics.json      V3: blocks, adjustment, ibl, fixed_stations,
                               mismatch, the child's own bulk trace
  <child>/v4_metrics.json      V4: the comparison against the baseline child
  <child>/child_vs_parent_metrics.json   the V1/V2 diagnostics, as context
  <child>/*.csv                streamwise_blocks_{core,canopy}, ibl_depth,
                               reference_profiles, profiles_vs_baseline,
                               tke_error_vs_distance_{x,y}, plus the V1/V2 CSVs
  v3_summary.{md,json}         the standoff table and the P-a/P-b/P-c verdict
  v4_summary.{md,json}         the difference-from-V1 table and its verdict
```

`v3_summary.md` and `v4_summary.md` are the deliverables; read them first, then
the `verdict` block, then a per-child JSON if a row needs explaining.

`periodic-stats` writes `analysis/periodic_<key>.json`, which is both the
equilibrium reference the adjustment length is measured against *and* the source
of the flat parent's `uflowrate`.  A V3 re-run that starts after that stage reads
the cached file; one that skips it and rebuilds the parent case falls back to the
preset's nominal value and says so in the log.

## Cost

Both are sized from the V1 `converged` job's measured stage times (PBS 3991175,
recorded in `.github/skills/udales-exec/references/clusters.md`), not from a
blanket allowance.  The two PBS headers carry the full derivation.

| | V3 | V4 |
|---|---|---|
| periodic runs | reference ~15 min, flat parent ~71 min, cubes parent ~84 min | staggered parent ~3 h 15 (+ nestparent I/O) |
| slab cuts | 5 x ~11.8 min (nestparent read, 0.5 s nesting.inp write) | ~23.9 min |
| child runs | 5 x ~15.5 min | ~27 min |
| analysis | ~31 min | ~25 min |
| **estimate** | **~5 h 53** | **~4 h 32** |
| **walltime requested** | **8 h** (1.35x) | **6 h** (1.32x) |
| new disk | ~866 GB | ~529 GB |
| `mem=` | 64 GB | 96 GB |

Revised from the pre-nestparent figures (V3 5 h / 245 GB, V4 6 h / 249 GB) once
every parent moved to the 0.5 s `&NESTPARENT` cadence and V3 gained its fifth
child; the two PBS headers carry the full derivation, including which lines
are still V1-measured rates and which are scaled estimates.  Almost all of the
extra disk is `nesting.inp` itself, which now stores the boundary at 6x the
temporal resolution (0.5 s against 3 s) -- not `&NESTPARENT`, whose own write
volume is a couple of times the old single-cadence figure at most, which is
the saving the parent-side zone dump exists to deliver (nesting-plan-
2026-09-06.md section 2: a full 0.5 s field dump of these parents would be
~10-75x larger again).  `mem=` is unchanged: the slab cut has been a streaming
writer since the D1 zone-dump work (commit 71011e96), so its peak no longer
scales with the nesting file's size the way the pre-nestparent estimate assumed.

Tighter than V1's and V2's blanket 8 h where the numbers support it (V4 still
is); V3's second periodic run is the least-measured figure in either estimate
(scaled from the reference's own small-case rate, not separately timed), so
its walltime carries more headroom than a tight backfill number would.  Note
what a short request does and does not buy: on CX3 a 64-cpu job is routed by
`ncpus` and walltime band only -- `<= 64` cpus and `<= 24` h both go to
`v1_medium24` -- so any walltime up to 24 h lands in the *same* queue, and the
whole benefit of a tighter number is backfill
eligibility (`backfill_depth = 1` on that queue).  If it stays deep, the other
lever is `ncpus`: `<= 16` routes to `v1_small24`/`v1_small72`, which were running
thousands of jobs when this was written.  Whether that trades well is not settled
by the existing measurements -- `clusters.md` records 5.4e6 cell-steps/s on 4
ranks and 1.7e7 on 64, but for *different case sizes*, so "3.1x for 16x the
ranks" is a warning that scaling is far from perfect rather than a measured
efficiency.  Sizing a 16-rank run would need its own timing, and
`nprocx`/`nprocy` changed in `presets_geometry.py`.

## Status

**Prepared, exercised end to end at tiny size on the Release binary, not yet
submitted / submitted (see job ids below).**  Revised 2026-09-07 for the fine
`&NESTPARENT` cadence and V3's new cleared-parent-cubes arm.  Both tiny
experiments run the identical production code path on a login node in a few
minutes: an equilibrium reference, a flat parent calibrated from it, three
standoffs plus the cleared-parent-cubes child (V3, now four periodic runs and
five children); an aligned parent, a staggered parent, a matched child and a
mismatched child (V4).  `test_v34_tiny.py` covers both plus the production
configuration.

What the tiny runs showed, **for orientation only** -- a handful of samples
over a 51 s window on a 40 s spin-up cannot support any physical claim:

* the volume-flow calibration works: the reference measured a bulk of
  2.5493 m/s and the flat parent then held 2.5434 m/s, 0.2 % low;
* the premise of V3 holds at tiny size -- the imposed canopy-layer velocity was
  2.3135 m/s against the equilibrium's 0.8129, and the roof-level resolved
  stress had the wrong sign in the flat parent, so there is a large mismatch
  for the child to work off;
* the streamwise blocks resolve it: with a 15-cell standoff, the first canopy
  row's residual mean-flow error is 122.4 % against the flat-parent
  `standoff0`'s 182.0 % and the cleared-parent-cubes arm's 5.0 % -- the child
  that inherits a real (if displaced) canopy's turbulence from its boundary
  starts far closer to equilibrium than either standoff arm, which is the
  qualitative direction the cleared-parent-cubes arm exists to check;
* every child is driven from the parent's `&NESTPARENT` band, not a full-domain
  dump (`manifest.json`: `driving_source = "nestparent"`), at `C_dump = 0.75`
  (`u0 = 3` m/s) and the Catmull-Rom cubic interpolant;
* `test_v34_tiny.py` is 37 tests, all passing in 388 s on a fresh Release
  build (was 29 in 184 s before the fine cadence and the new arm);
* every child ran `nest_lparentgeom = .false.` with all four faces forced,
  `Phi` at round-off (`~1e-14` in the stored file) and `divmax` at round-off,
  which is the point: a canopy generated with no reference to the parent's
  geometry still lets the solver *assert* the design section 5 rule;
* the cleared-parent-cubes child's solid mask is identical, cell for cell, to
  `standoff0`'s -- the two children really are the same canopy, and only the
  driving parent (flat vs. V1's own aligned array) differs;
* the V4 comparison ran against a baseline child whose solid mask was checked
  identical to the mismatched child's, and refused a baseline of the wrong shape.

No production number exists yet.  Nothing in this section should be read as one.

### Job status

| job | experiment | submitted | job id |
|---|---|---|---|
| V3 | `v3` (5 children: 4 standoffs + cleared-parent-cubes) | 2026-09-07 12:35 UTC | `3996513.pbs-7` |
| V4 | `v4` (1 child: mismatch, against the V1 `converged` baseline) | 2026-09-07 12:35 UTC | `3996514.pbs-7` |

Both queued in `v1_medium24` at submission (`qstat -u $USER`); check
`$EPHEMERAL/nesting-v3/analysis/v3_summary.md` and
`$EPHEMERAL/nesting-v4/analysis/v4_summary.md` once they finish.

---

# V6 -- does mass drift over long nested runs?

Design section 10.4 row V6: `10^5`-step run -> `divtot` bounded, not drifting.
Every other row asks whether nesting reproduces some physical quantity; this
one asks nothing about physics at all -- only whether the scheme's own
bookkeeping (the pressure projection, the flux-corrected boundary, the
relaxation forcing) stays bounded when it is asked to run for a very long
time. It was also the only row in the table that had never been run.

## The construction

The cleanest way to isolate numerical drift from everything else is to remove
everything else: drive the child from a boundary that is **constant in
time**, then run far beyond it, and read any growth in the solver's own
conservation diagnostics as drift and nothing else.

`src/nesting_scheme.f90`'s `check_record_end` already has the mechanism.  Past the
last stored parent time level, `read_level` and `eval_target` clamp their
level index to `ntime` (confirmed by reading the code -- `it = min(max(ilev,
1), ntime)` and the two clamped calls to `nestio_hdr%time` in `eval_target`),
so `tlo = thi` and the interpolation weight `th` degenerates to `0`: the
"interpolated" target is just the last stored level, held.  This is fatal by
default (`nest_lendabort = .true.`) -- an aborted run means the user's setup
outlived its own boundary data, which is a real mistake worth catching -- but
`nest_lendabort = .false.` reports it once (`nendwarn`) and lets the run
continue, frozen.  `config.V6` sets exactly that.

The parent stores only **8 dumped levels** (`dtdump = 3 s`, `production =
24 s`, a `spinup` of 15 s just long enough for the cold random initial
condition to be pressure-projected before the first dump -- the geometry is
TINY's own 96x96x32 aligned-cube-array parent / 64x64x32 child, reused as-is
because this row is about numerics, not turbulence, so the cost belongs in
the step count rather than the grid).  `run_v6.py` then patches `RUN.runtime`
out to roughly 100,000 steps' worth of simulated time -- about 3.7x the
stored record's own ~21 s span -- so the boundary is time-invariant for the
overwhelming majority of the run.  `nest_timeinterp = 2` (Catmull-Rom cubic
Hermite) rather than TINY's linear default: the freeze exercises
`eval_target`'s `h2 <= 0` fallback, which only exists on the cubic path, so
this is the interpolant that actually tests the frozen-boundary code.

`test_v6_tiny.py` checks the construction itself before trusting it for
anything: the freeze warning fires **exactly once** (`nendwarn`'s guard) and
the run does not abort, at a scale that runs on a login node in about a
minute past the record.

## What is measured

Everything comes from the solver's own stdout, at two independent, patched
throttles (`caselib.set_namoption`, since the right cadence depends on a
measured step rate `config.py` cannot know in advance):

* `NAMCHECKSIM.tcheck` -- `chkdiv`'s `divmax`, `divtot`
  ([modchecksim.f90:161](../../../src/modchecksim.f90#L161));
* `NESTING.nest_statint` -- `nesting_stats`'s flux residual `Phi` (norm and
  lid/closed-face split), the zone misfit rms, and `|grad p|` zone / interior
  / ratio (design section 6.4).

Both are set to the same value so the two clocks line up.  `child_dtdump` is
set to 300000 s -- far beyond any runtime this experiment uses -- so the
child never writes a single field dump; nothing here needs 3-D output.

`analyse_v6.py` fits an ordinary-least-squares trend to each series over the
whole run and calls it a **fail only when the slope is both statistically
distinguishable from zero** (more than 3 standard errors from 0) **and large
enough to move the series by a real fraction of its own range** (more than
20 % of peak-to-peak). Either alone is not enough: round-off noise over a
long run gives a "significant" slope of a negligible size, and a short noisy
window can show a large swing that is not a trend at all. A flat, bounded
series is a pass; a monotone trend is a fail, reported as one, not smoothed
over.

## Cost, from measurement

Two probes on a login node (`UDALES_BUILD=build/release/u-dales`, the child
on its own single rank):

| probe runtime | steps | main-loop CPU time | rate |
|---|---|---|---|
| 300 s | 793 | 61.67 s | 12.86 steps/s |
| 900 s | 2383 | 186.81 s | 12.76 steps/s |

("`TOTAL CPU time by main time loop`", printed at the end of `child.log`;
this excludes the run's own fixed ~13-20 s start-up, which does not grow with
the run length.)  Mean `dt` was 0.3775-0.3777 s in both probes.  The two
measurements agree to within 1 %, so the mean rate (12.81 steps/s) is used
directly rather than extrapolated from one sample.

For 101,000 steps (a deliberate ~1 % overshoot of the 1e5 target, since the
actual count tracks the adaptive `dt` actually taken): `--runtime 38000`
(101000 x 0.3775 s, rounded up), main-loop wall time ~0.82 h, everything else
(parent build + spin-up + production + child-case build + analysis) ~75 s
measured directly above. `submit_cx3_v6.pbs` requests `walltime=06:00:00`
(~7.1x that estimate -- more headroom than this campaign's other jobs,
deliberately: this is the first real run of the construction, only exercised
at probe scale so far) on `ncpus=4` (the parent's own 2x2; the child runs on
1x1 for almost the whole walltime -- "few cores is fine and preferable" for a
long, thin job). `--stat-interval 300` gives ~127 samples over the run: a
trend fit wants `>= 3` and this clears that by two orders of magnitude,
without the log growing past a few thousand lines.

## Layout

| File | What it is |
|---|---|
| `config.V6` | the preset: TINY's geometry, `nest_lendabort = False`, `nest_timeinterp = 2`, an 8-level parent record |
| `run_v6.py` | end-to-end driver; `--runtime` and `--stat-interval` patch the built case rather than being preset fields, because the right values depend on a measured step rate |
| `analyse_v6.py` | parses `child.log`, fits a trend per series, reports pass/fail |
| `test_v6_tiny.py` | harness smoke test (the construction, past the record, in about a minute) plus a pure-Python check of the trend/verdict logic against synthetic series and a synthetic log fragment |
| `submit_cx3_v6.pbs` | the CX3 production job, sized from measurement in its own header. **Review before submitting.** |

## How to run it

```bash
module purge && module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
source /rds/general/user/mvr/home/udales/.venv/bin/activate
export UDALES_BUILD=$PWD/build/release/u-dales

python tests/run_tests.py nesting-validation      # includes the V6 tiny entry
python tests/validation/nesting/test_v6_tiny.py   # just the smoke test

# a probe, to re-measure the step rate on different hardware
python tests/validation/nesting/run_v6.py $EPHEMERAL/v6-probe --runtime 300

qsub tests/validation/nesting/submit_cx3_v6.pbs   # the real thing
```

## Status

**Prepared, validated at tiny scale (15 tests, ~75-90 s on a login node), a
step rate measured from two probes (300 s and 900 s simulated, agreeing to
within 1 %), submitted.**  The full `tools/python/tests/test_nesting*.py`
suite (96 tests) was re-run afterwards and is unchanged.

| job | experiment | submitted | job id |
|---|---|---|---|
| V6 | `v6` (101,000-step frozen-boundary run) | 2026-09-08 00:16 UTC | `4001085.pbs-7` |

Queued in `v1_small24` at submission (`qstat -u $USER`), behind V0b's
`4000813.pbs-7`; check `$EPHEMERAL/nesting-v6/analysis/v6_summary.md` and
`$EPHEMERAL/nesting-v6/991/child.log` once it finishes.
