# Vegetation forcing (runmode 1013)

`vegetation_forcing` was the last host routine left inside the GPU time loop's
transfer window. This directory holds the host-side test that came with its
port; the device side is covered by `tests_cuda.f90::test_vegetation_forcing`
and by three GPU parity cases.

Run it with:

```
bash tests/integration/vegetation/run_test.sh
```

It runs twice, on one rank and on two.

## What the port changed

**Residency.** All five loops - the drag on the u, v and w face lists, the
canopy energy balance, and the scalar deposition - now run as OpenACC kernels
over the same sparse lists, against the field mirrors `modcuda` already holds.
The tendencies are added straight into `up_d`/`vp_d`/`wp_d`/`thlp_d`/`qtp_d`/
`svp_d`.

**Call site.** `call vegetation_forcing` moved from below `updateHost` to above
it. It was the one routine in that window that read the previous-step fields on
the host, and it read them at the *neighbours* of every tree point - which can
be a solid cell `ibmnorm` has just pinned on the device. That is why
`updateHost` carried an `ltrees` block copying `um`, `vm`, `wm`, `thlm` and
`qtm` back down every step. The block is gone; a tree run now costs the same as
any other. On a CPU build nothing moved: `updateHost` does not exist there and
the order against `EB` is unchanged.

**Cache.** Everything in the canopy loop that depends only on the vegetation
properties is now computed once, in `init_vegetation_cache`: both Beer-Lambert
exponentials, the leaf-size and stomatal-resistance clamps, and the divisions
by `rhoa*rlv`, `rhoa*cp`, `gam*rs` and `130*sqrt(lsize)`. Six divisions and two
exponentials per vegetation point per step become three divisions. Writing `s`
for the saturation-curve slope and `g` for the psychrometric constant:

- `s = 4098*e_sat/(T - 35.85)**2`, and the loop only needs the ratio
  `f = s/(s + 2g) = 4098*e_sat/(4098*e_sat + 2g*(T - 35.85)**2)`, which removes
  the slope and its division outright;
- `g/(s + 2g) = (1 - f)/2`, so `omega = 1/(1 + (1 - f)*rs/r_a)`;
- `r_a = 130*sqrt(lsize/sqrt(wind2)) = 130*sqrt(lsize)*wind2**(-1/4)`, so
  `rs/r_a = veg_rs_ra * wind2**0.25` with `veg_rs_ra` cached.

The cache also collapses `vegetation_forcing_legacy` and
`vegetation_forcing_sveg` into one routine. The only thing that separated them
was how `q_av_leaf` is formed, and that is now a table lookup times a
mode-dependent scalar - `Qstar` in legacy mode, `1` in sveg mode. `Qstar` is
deliberately left outside the cache: it is a namelist constant today, but
keeping it in the loop means the cache stays valid if `dQdt` is ever wired up.

The cache is used by the host branch too, so the two branches share one source
of truth and the CPU gets the same saving.

**Diagnostics.** `veg_up`/`veg_vp`/`veg_wp` and the `vegp%*` arrays are read by
`statsdump` under `ltreedump` and by nothing else. They stay on the device and
come down through `updateVegDiagHost`, gated on
`modstatsdump::statsdump_will_sample` - so the cadence follows `tsample`, not
`dt`, and a run without `ltreedump` moves nothing at all. This is the same
shape as the facet-integral drain above `EB`.

## Measured effect

`nsys`, case 526 on one rank, 20 timesteps (60 Runge-Kutta stages), GPU Release
build, before and after the port:

| | before | after |
|---|---|---|
| device-to-host bytes | 9,469.1 MB | 8,112.5 MB |
| device-to-host operations | 2,237 | 2,097 |
| device-to-host time | 865.6 ms | 703.4 ms |
| host-to-device bytes | 9,491.7 MB | 9,492.0 MB |

The 1,356.6 MB that stopped moving is 299.5 copies of a 4.53 MB field, which is
exactly the 5 arrays x 60 stage-calls the `ltrees` block in `updateHost` was
doing. The operation count falls by only 140 rather than 300 because the
diagnostic drain adds about 160 sub-megabyte copies in their place - 8 sparse
arrays on each of the 20 sampling steps.

Host-to-device traffic was untouched by this change: `updateDevicePriorPoiss`
still ran, and with `vegetation_forcing` gone from the window it sat directly
after `updateHost` with no host work between them. That is what made it
removable, which was the next change rather than this one - it took host-to-
device down from 9,492.0 MB to 6,522.6 MB on the same case, and its guard is
`tests_cuda.f90::test_bc_profile_upload`.

## What this test pins

The reference is not the new code in another form. `canopy_reference` inside
`tests_vegetation` writes out the *original* expressions from the raw
vegetation properties, keeping `slope_sat`, `r_a` and `gam` explicit and
sharing nothing with `vegetation.f90`. If the cache or the folded divisions
changed anything, that is where it shows.

On top of that:

- the drag on all three face lists, and that nothing outside those lists moves
  - checked over the whole declared extent, halos included;
- `qt = qtR + qtA`, so the diagnostic split stays exhaustive;
- `qe + qh = q_av_leaf`, recovered from the two tendencies through their own
  scalings, so the absorbed radiation has to go somewhere;
- the mode dispatch, by running legacy and sveg over the same fields and
  requiring each to match its own reference *and* to differ from the other;
- drag-only mode leaving heat and moisture untouched while still applying drag;
- scalar deposition, per component.

### Non-degeneracy

The shipped vegetation cases give every point the same `lad`, `rs`, `lsize`,
`ud` and `dec`, taken from one set of `&TREES` namelist values. A per-point
index mix-up would therefore reproduce the right answer everywhere. The test
imposes a bounded, non-monotone per-point profile on all of them plus `sveg`,
and asserts the result is not flat before trusting any comparison. The drag
coefficients `dcoef_u/v/w` are deliberately left alone - they were built from
the original `lad` and are read straight out by both the routine and the
checks, so the two paths stay independent.

### Two ranks

With the domain split in x the single tree of case 526 lands on one rank and
the other has no vegetation points at all. `init_vegetation` broadcasts and
barriers, and the test calls it once per radiation mode, so no rank may return
early: the point count is reduced across ranks instead, and a rank with nothing
local runs every collective and skips only the per-point checks.

## Mutation results

Every mutation below was applied to the shipped source, built, and run against
this test (host, `m*`) or against the GPU parity cases with
`--require-debug-selftest` (device, `d*`; sampling predicate, `s*`). All were
caught; the unmutated baseline passes.

### Host branch - `tests_vegetation`

| # | Mutation | Caught by |
|---|---|---|
| m0 | none (baseline) | - (passes) |
| m1 | `veg_rs_ra` drops the `sqrt` on `lsize` | decoupling factor |
| m2 | `omega` uses `1 + f` instead of `1 - f` | decoupling factor |
| m3 | `f` doubles the psychrometric term | decoupling factor |
| m4 | legacy `veg_qleaf` swaps the two exponentials | moisture and heat tendency |
| m5 | `qh = qe - q_av_leaf` | heat tendency, energy closure |
| m6 | `veg_lad_rlv` divides by `cp` instead of `rlv` | moisture tendency, both shares |
| m7 | `veg_aero` uses `rlv` instead of `cp` | moisture and heat tendency |
| m8 | u-drag reads `wm(i,j,k)` where it wants `wm(i,j,k+1)` | u drag |
| m9 | `up` gets `- veg_up` instead of `+ veg_up` | u drag scattered into up |
| m10 | deposition adds instead of subtracts | scalar deposition |
| m11 | sveg mode falls through to the legacy formula | sveg tendencies |
| m12 | `wind2**0.5` instead of `wind2**0.25` in `rs/r_a` | decoupling factor |
| m13 | u face list indexed at `i + 1` | u drag |
| m14 | `veg_lad_ud` drops `ud` | scalar deposition |
| m15 | sveg `veg_qleaf` drops the division by `lad` | sveg tendencies |

### Device branch - `test_vegetation_forcing` and the parity cases

| # | Mutation | Caught by |
|---|---|---|
| d0 | none (baseline) | - (passes) |
| d1 | u drag never scattered into `up_d` | self-test: `up` |
| d2 | `wind2**0.5` instead of `wind2**0.25` | self-test: canopy qt |
| d3 | `vegp_qt_d` accumulates instead of assigning | self-test: tendencies accumulate across calls |
| d4 | deposition sign flipped | self-test: scalar deposition, and the deposition parity case |
| d5 | u drag reads `wm_d(i,j,k)` for `wm_d(i,j,k+1)` | self-test: u drag |
| d6 | `updateVegDiagHost` skips `vegp%thl` | self-test: drained thl |
| d7 | canopy kernel drops the `Qstar` factor | self-test: canopy qt |

### Sampling predicate

| # | Mutation | Caught by |
|---|---|---|
| s1 | `statsdump_will_sample` always false | all three parity cases |
| s2 | fires one step late (`tsamplep >= tsample + dt`) | the `vegetation-sveg` cadence case only |

`s2` is why `vegetation-sveg` runs five steps with `tsample = 2*dtmax` instead
of the single step the other two use: a one-step case samples on its only step,
so an off-by-one in the predicate is invisible to it.

## Known limitation, unchanged by this port

`vegetation_forcing` returns immediately when the local rank has no vegetation
points of its own. The staggered face lists are built from `dcoef_3d` *after* a
halo exchange, so a rank with no vegetation can still have `npts_u > 0` at a
face it shares with a vegetated neighbour - and that drag is dropped. This is
pre-existing behaviour, reproduced exactly by the port and documented at the
early return in `vegetation_forcing`; it is not something the port introduced,
and fixing it would change results at rank boundaries.
