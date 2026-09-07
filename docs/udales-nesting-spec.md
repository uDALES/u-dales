# uDALES nesting — implementation contract

Companion to `docs/udales-nesting-design.md`. **This file is normative**: every work package codes
against the names, shapes and semantics below. Do not change them unilaterally — if something here
is wrong, say so rather than silently diverging.

Rationale for every design choice is in the design document; this file is the interface only.

---

## 1. Build and verify

```bash
cd /rds/general/user/mvr/home/udales/u-dales
./tools/build_executable.sh icl release      # ifort 2021.2; ~3 min; must end "Built target u-dales"
```

A pre-change baseline binary is at `build/u-dales.baseline`.

Python (separate module stack — `module purge` first):
```bash
module load tools/prod && module load Python/3.9.6-GCCcore-11.2.0
source /rds/general/user/mvr/home/udales/.venv/bin/activate
python -m unittest discover -s tools/python/tests -p 'test_nesting*.py'
```

## 2. Fortran conventions

* Free-form F90, `implicit none`, `save`, `private` + explicit `public` list, as in `modibm.f90`.
* `use modx, only : a, b` — always `only`.
* Two-space indent; match the surrounding file's style exactly.
* Real kind: plain `real` (the build is `-r8` equivalent via 2DECOMP `DOUBLE_PREC`).
* No trailing whitespace. No tabs.
* Every public routine gets a `!>` doc comment naming its call site.
* Guard every routine with `if (.not. lnesting) return` as the first executable statement.

## 3. Constants (owned by `src/modglobal.f90`)

```fortran
integer, parameter :: BCxm_nesting = 4
integer, parameter :: BCym_nesting = 3

integer, parameter :: TEST_NESTING_WEIGHTS  = 1006
integer, parameter :: TEST_NESTING_GEOMETRY = 1007
integer, parameter :: TEST_NESTING_IO       = 1008
integer, parameter :: TEST_NESTING_FLUX     = 1009
integer, parameter :: TEST_NESTING_UPDATE   = 1010
integer, parameter :: TEST_NESTING_INIT     = 1011
```

## 4. Namelist `&NESTING` (read in `modstartup`, broadcast to all ranks)

| Name | Type | Default | Meaning |
|---|---|---|---|
| `lnesting` | logical | `.false.` | master switch |
| `nestfile` | char(256) | `''` | empty ⇒ `nesting.inp.<expnr>.nc` |
| `nest_guardwidth` | real | `0.` | $L_{\rm imp}$ [m] |
| `nest_zonewidth` | real | `0.` | $L_{\rm rel}$ [m] |
| `nest_tau` | real | `0.` | [s]. The relaxation rate is `W/tau`, so `tau<=0` means an infinite rate **wherever `W>0`**, i.e. the whole zone becomes Dirichlet and the ramp is defeated. It is therefore only legal when `nest_zonewidth == 0` (pure guard strip, no ramp); `checkinitvalues` enforces this. |
| `nest_shape` | integer | `1` | 1 raised cosine, 2 quintic |
| `nest_lateral(4)` | logical | `.true.` | W, E, S, N. A nested direction must impose **both** of its faces: with `BCxm = BCxm_nesting`, `nest_lateral(1:2)` must both be `.true.`, and likewise `nest_lateral(3:4)` under `BCym_nesting`. The nesting BC replaces the convective outflow as well as the inflow, so a face left out would get neither and its ghost plane would never be set. `checkinitvalues` and `nesting_init` both reject it. |
| `nest_top` | logical | `.false.` | Case C — **not implemented in v1**, must error if `.true.` |
| `nest_timeinterp` | integer | `2` | Default 2: cubic Hermite with **unlimited Catmull-Rom** slopes; 1 = linear. The C0 cadence experiment (2026-09) showed the Hermite halves the interior TKE deficit at the same parent cadence (−9.9 % → −5.7 % at 3 s). The slopes MUST stay unlimited: the interpolant has to be linear in the data or it breaks the flux compatibility of design §3.1 (U44). |
| `nest_nwall` | integer | `1` | wall erosion, cells |
| `nest_lparentgeom` | logical | `.false.` | parent resolves the child geometry |
| `nest_fluxtol` | real | `1.e-10` | abort threshold on \|Φ\| (normalised, §7) |
| `nest_lfluxassert` | logical | `.true.` | **on by default** |
| `nest_lfluxcheckall` | logical | `.false.` | Recompute the flux residual of **every** stored time level from the boundary slabs at init (4 × `ntime` reads on every perimeter rank) instead of validating the residual the writer stored. Off by default; the reader falls back to the recompute on its own for a schema 1 file, or when the file's `fluid_lateral_area` is not this run's. |
| `nest_statint` | real | `-1.` | Interval [s] between `nesting_stats` reports. `< 0` means "use `tstatsdump`"; `0` reports every timestep. A report costs three `MPI_ALLREDUCE`s, three full-domain sweeps and seven lines of stdout, so it is throttled like `statsdump` rather than run every substep. Whatever the interval, the first and the last timestep always report; the energy-injection accumulators sum over every substep since the previous report, and the Φ lines carry the largest-magnitude residuals (six faces, lid, closed faces) seen since the previous report rather than one substep's instantaneous value. With `nest_lfluxassert = .false.` the flux residual is only evaluated on the substep that reports, so those maxima reduce to it. |
| `nest_lendabort` | logical | `.true.` | What happens when the simulation time passes the last stored parent level (`time(ntime)`), beyond which the buffer clamps and the boundary freezes on that level. `.true.`: abort — at initialisation already if `timee + runtime` lies beyond the record, at the crossing otherwise. `.false.`: warn once and run on frozen. One timestep of grace is allowed at the end. A single-level file is a steady parent and never triggers this. |
| `nest_linitfromparent` | logical | `.false.` | **Cold start only.** Fill `u0`/`um`, `v0`/`vm`, `w0`/`wm` from the file's full-domain initial-condition block instead of from `prof.inp`. Requires a schema 2 file with `has_initial_condition = 1`; errors otherwise. Ignored, with a message, on a warm start — the restart file already holds a state consistent with the parent and overwriting it would break restart parity. |

`checkinitvalues` must, when `lnesting`:
1. require `ipoiss == POISS_FFT2D`, else stop;
2. require `BCxm == BCxm_nesting .or. BCym == BCym_nesting`, else stop;
3. **not** force `BCtopm = BCtopm_pressure` (unlike the `BCxm_profile`/`BCxm_driver` branches);
4. stop if `nest_top`;
5. stop if `nest_guardwidth <= 0.` or `nest_zonewidth < 0.`;
6. stop if `nest_tau <= 0.` while `nest_zonewidth > 0.` (see the `nest_tau` row above);
7. stop if a nested direction does not impose both of its faces (see the `nest_lateral` row).

## 5. File format `nesting.inp.<expnr>.nc`

**Two schema versions are normative.** Version 1 is the original file. Version 2 adds

* `flux_residual(time)` — the net boundary flux of the data **as stored**, i.e. *after* any
  divergence correction — together with the `fluid_lateral_area` attribute giving the area it was
  summed over, so the solver can validate every stored level at initialisation without re-reading
  the boundary slabs;
* an **optional** full-domain initial condition `u_init`/`v_init`/`w_init`, flagged by the
  `has_initial_condition` attribute, for `nest_linitfromparent`.

**Both versions must load and run.** A version 1 file has neither addition and behaves exactly as
before: the solver recomputes the residual from the boundary slabs, with a warning, and
`nest_linitfromparent` is an error against it. Everything version 2 adds is therefore optional on
read and required only of a file that declares `udales_nesting_schema = 2`.

CDL dimension order below; **Fortran sees the reverse**. The decomposed index is deliberately the
outermost spatial dimension so a rank's hyperslab is one contiguous run.

```
dimensions:
  time = UNLIMITED ;
  zf = ktot ;  zh = ktot+1 ;
  xf = itot ;  xh = itot+1 ;  yf = jtot ;  yh = jtot+1 ;
  nz  = nzone ;  nzh = nzone+1 ;
variables:
  double time(time) ;
  double xf(xf), xh(xh), yf(yf), yh(yh), zf(zf), zh(zh) ;
  double rhobf(zf), rhobh(zh) ;
  double net_volume_flux(time) ;   // net boundary flux BEFORE correction (provenance).
                                   // Post-correction it is ~1e-16 and carries no information,
                                   // so readers must NOT assume this variable is zero.

  // SCHEMA 2 ONLY, required:
  double flux_residual(time) ;     // net boundary flux of the data AS STORED, same functional
                                   // and same sign convention as net_volume_flux. This is what
                                   // the solver validates at init, normalised by its own total
                                   // fluid boundary area, instead of re-reading the slabs.

  // SCHEMA 2 ONLY, present iff has_initial_condition = 1:
  double u_init(xh, yf, zf) ;      u_init:stagger = "xh yf zf" ;
  double v_init(xf, yh, zf) ;      v_init:stagger = "xf yh zf" ;
  double w_init(xf, yf, zh) ;      w_init:stagger = "xf yf zh" ;
                                   // The full child-grid velocity at time(1). Discretely
                                   // solenoidal on the child grid, with the boundary-normal
                                   // velocities equal to the (corrected) slab values at that
                                   // time and w = 0 on the floor and the lid, so a cold start
                                   // from it begins divergence free and consistent with the
                                   // imposed boundary. Fortran sees (z, y, x): a rank's (i,j)
                                   // block is contiguous in z.

  // west / east slabs: decomposed index is y
  double u_west (time, yf, zf, nzh) ;   u_west:stagger  = "xh yf zf" ;
  double v_west (time, yh, zf, nz ) ;   v_west:stagger  = "xf yh zf" ;
  double w_west (time, yf, zh, nz ) ;   w_west:stagger  = "xf yf zh" ;
  //  ... u_east, v_east, w_east: identical shapes

  // south / north slabs: decomposed index is x
  double u_south(time, xh, zf, nz ) ;   u_south:stagger = "xh yf zf" ;
  double v_south(time, xf, zf, nzh) ;   v_south:stagger = "xf yh zf" ;
  double w_south(time, xf, zh, nz ) ;   w_south:stagger = "xf yf zh" ;
  //  ... u_north, v_north, w_north

// global attributes (all required):
  :Conventions = "CF-1.8" ;  :udales_nesting_schema = 1 or 2 ;
  :divergence_corrected = 1 ;
  :itot = ; :jtot = ; :ktot = ; :nzone = ; :xlen = ; :ylen = ;
  :parent_model = ; :parent_dx = ; :parent_dt = ;
  :child_origin_x = ; :child_origin_y = ; :rotation_deg = 0. ;
  :created = ; :creator = ; :tool_version = ;
  :child_dt = ;        // OPTIONAL: child timestep; without it the temporal-refinement guard is skipped
  :parent_dy = ; :parent_dz = ;   // OPTIONAL: smallest parent spacings in y and z; the writer's
                                  // spatial-refinement guard is the largest per-axis ratio it knows

// SCHEMA 2 ONLY, both required:
  :has_initial_condition = 0 or 1 ;   // whether u_init/v_init/w_init are present
  :fluid_lateral_area = ;             // GEOMETRIC (no density) fluid area of the four lateral
                                      // boundary faces that flux_residual was summed over. The
                                      // solver compares it against its own IIu/IIv area and
                                      // falls back to the full recompute if they disagree, so a
                                      // mask mismatch between writer and solver is caught rather
                                      // than trusted.
```

**Coordinates are child-relative.** `xf`, `xh`, `yf`, `yh` are the child's own grid, starting at
0 -- what `nestio_validate` compares against the run's `xh(1:itot+1)` to `nestio_tol = 1e-10` of
`xlen`. The child's position in the parent is carried **only** by `child_origin_x`/`child_origin_y`
(metres, parent coordinates); a writer working in parent coordinates subtracts them before
writing. `zf`/`zh` are the run's vertical and are not shifted.

**Index conventions.** The west slab covers child cells `i = 1..nzone` (centres) and faces
`i = 1..nzone+1`. The east slab covers centres `i = itot-nzone+1..itot` and faces
`i = itot-nzone+1..itot+1`; **its first slab index corresponds to the lowest global index**, i.e.
slab index `m` ↔ global `i = itot-nzone+m` for centres and `i = itot-nzone+m` for faces
(so `m = nzone+1` ↔ `i = itot+1`). South/north follow the same rule in `j`.

**Units** `m s-1` for the slabs; `net_volume_flux` and `flux_residual` carry `m3 s-1` (they are
`sum(rho u_n dA)` with `rhobf == 1`, see the density convention below). Missing/NaN values are an
error, not a sentinel: the reader checks every slab and every initial-condition block it reads
(`nestio_check_values`) for non-finite elements and for the variable's `_FillValue` (netCDF's
default fill for doubles when none is declared) and aborts, naming the variable, the time level
and the first offending element. A NaN would otherwise pass the flux assertion silently, since
`abs(NaN) > tol` is false.

**Density convention.** Fluxes are weighted by `rhobf(k)` on both the writer and the solver side,
matching DALES's `openboundary_divcorr`. uDALES's own Poisson RHS carries no density (design F1),
so the two agree exactly while `rhobf == 1`, which is always. Both sides must use the same
convention or the runtime flux assertion can fire spuriously.

The one place the two conventions cannot be reconciled is the initial-condition block, which is
projected with the solver's own **density-free** divergence operator: a file carrying an initial
condition must therefore have `rhobf == rhobh == 1`, and the writer refuses anything else rather
than storing a field whose boundary flux does not close.

## 6. `src/modnestingio.f90` — input only, no scheme knowledge

Must not `use modnesting`. Standalone and separately compilable.

```fortran
module modnestingio
  implicit none;  save;  private
  public :: nestio_open, nestio_validate, nestio_read, nestio_read_block, &
            nestio_close, nestio_hdr, nestio_header_type

  !> Schema versions this reader understands. A schema 1 file must keep loading
  !! and running exactly as before, so everything schema 2 adds is OPTIONAL here.
  integer, parameter :: NESTIO_SCHEMA_MIN = 1
  integer, parameter :: NESTIO_SCHEMA_MAX = 2

  type nestio_header_type
    integer :: schema = 0, itot = 0, jtot = 0, ktot = 0, nzone = 0, ntime = 0
    real    :: xlen = 0., ylen = 0., rotation_deg = 0.
    logical :: divergence_corrected = .false.
    logical :: has_flux_residual = .false.      ! schema 2
    real    :: fluid_lateral_area = 0.          ! schema 2
    logical :: has_initial_condition = .false.  ! schema 2
    real, allocatable :: time(:), xf(:), xh(:), yf(:), yh(:), zf(:), zh(:)
    real, allocatable :: rhobf(:), rhobh(:), net_volume_flux(:), flux_residual(:)
  end type
  type(nestio_header_type) :: nestio_hdr

  !> Open read-only on every rank and populate nestio_hdr. ierr/=0 on failure.
  subroutine nestio_open(fname, ierr)
    character(len=*), intent(in)  :: fname
    integer,          intent(out) :: ierr

  !> Compare the header against modglobal (itot,jtot,ktot,xlen,ylen,xf,xh,yf,yh,zf,zh),
  !! the schema version, AND the per-variable `stagger` attribute of every slab
  !! variable present. On mismatch, write a message naming the offending field
  !! and its two values, then stop 1. The stagger check matters because a file
  !! can have a correct grid and still lay the data out at the wrong staggered
  !! location, which would otherwise be read silently and wrongly.
  !! When has_initial_condition is set, u_init/v_init/w_init are additionally
  !! checked for presence, shape and stagger, on the same terms.
  subroutine nestio_validate()

  !> Read one time level of one variable, for this rank's range of the
  !! decomposed index. varname is e.g. 'u_west'. start2/count2 are 1-based in
  !! the decomposed (outermost spatial) index. buf is (n_zone_dim, n_z_dim, count2).
  subroutine nestio_read(varname, it, start2, count2, buf, ierr)
    character(len=*), intent(in)  :: varname
    integer,          intent(in)  :: it, start2, count2
    real,             intent(out) :: buf(:,:,:)
    integer,          intent(out) :: ierr

  !> Read a rectangular block of one full-domain, time-independent variable
  !! (u_init/v_init/w_init, schema 2). start2/count2 index y and start3/count3
  !! index x, 1-based; the whole vertical is read. buf is (nz, count2, count3),
  !! matching the file's Fortran dimension order (z, y, x).
  subroutine nestio_read_block(varname, start2, count2, start3, count3, buf, ierr)
    character(len=*), intent(in)  :: varname
    integer,          intent(in)  :: start2, count2, start3, count3
    real,             intent(out) :: buf(:,:,:)
    integer,          intent(out) :: ierr

  subroutine nestio_close()
end module
```

`nestio_read` must report the elapsed wall time cumulatively in a public
`real :: nestio_tread = 0.` so the I/O fraction can be reported (design §6.2).

## 7. `src/modnesting.f90` — the scheme

**Call order (normative).** `nesting_init` must be called *after* `readinitfiles`, because it
positions the parent time buffer on `timee` and `readinitfiles` is what assigns `timee`
(`modglobal.f90:441` declares it with no initialiser). Its own prerequisites `createmasks` and
`calcfluidvolumes` must still precede it. Getting this wrong traps on a Debug build and silently
breaks restart parity on a Release one — see design §9.5.

Public surface exactly as in the design document §9.2/§9.3, plus these **test hooks**, which must
be the same code the solver uses (no duplicated logic):

```fortran
  !> Shape function W(s). ishape: 1 raised cosine, 2 quintic.
  !! W=1 for s<=limp; 0 for s>=limp+lrel; C1 (ishape=1) or C2 (ishape=2) between,
  !! with zero slope at BOTH ends.
  real function nest_shape_fn(s, limp, lrel, ishape)

  !> Bounded union W = 1 - product(1-w(1:n)).
  real function nest_union(w, n)

  !> Physical coordinate of staggered variable ivar (1=u,2=v,3=w) at LOCAL indices i,j,k.
  !! u -> (xh(ig), yf(jg), zf(k));  v -> (xf(ig), yh(jg), zf(k));  w -> (xf(ig), yf(jg), zh(k))
  !! with ig = i + zstart(1) - 1, jg = j + zstart(2) - 1.
  subroutine nest_stagger_coord(ivar, i, j, k, x, y, z)

  !> Net volume flux of the predicted velocity through the domain boundary,
  !! summed over FLUID faces only (IIu/IIv/IIw), MPI-reduced over comm3d,
  !! normalised by the total fluid boundary area so the tolerance is dimensionless.
  !! This is Phi over all six faces (design section 3.1).
  real function nest_flux_residual(pup, pvp, pwp, rk3coef)

  !> The same quantity, split into the total and the part the LID carries.
  subroutine nest_flux_split(pup, pvp, pwp, rk3coef, phi_all, phi_lid)
```

**What the flux assertion asserts.** `nesting_bcpup` tests `phi_all` under a rigid lid
(`BCtopm_freeslip`/`BCtopm_noslip`, design case A) and `phi_all - phi_lid` under a leaky one
(`BCtopm_pressure`, case B). Under a rigid lid `bcpup` forces `w* = 0` at `ke+1`, so `phi_lid` is
identically zero and the two are the same number. Under a leaky lid the top face is not a datum:
`bcpup` sets `w*` there from the accumulated pressure and `tderive` adds the matching increment,
which is exactly the Dirichlet-in-the-mean-mode row the solver pins (design F3). The projection is
therefore complete for any net flux, `Phi_total = 0` is **not** a solvability requirement, and the
lid flux is the child breathing against its reservoir rather than an error. What must still vanish
is the flux through the faces the scheme controls. `nest_lfluxassert` therefore stays on by default
in both cases.

**Faces inside solids.** A boundary face or ghost value whose first interior point is solid at the
component's own stagger (`IIu`/`IIv`/`IIw`) is set to **zero**, not to the parent value, in both
`nesting_bcpup` and `nesting_boundary` (design §4). This is what makes the fluid-face $\Phi$ the
whole flux the Poisson right-hand side sees: an unmasked parent value on a solid face would push its
flow through the building, and since the assertion sums fluid faces only it could not see the
divergence source that creates. The writer must balance the fluid faces only (`FaceMasks`) for the
same reason. `nesting_bcpup` also sets `um`/`u0` at the imposed faces to the target it imposes, so
the integrated field is consistent with the projection at every substep rather than only after the
next `boundary` call. Unit test U45 pins both.

**Target time.** `nesting_update_target` evaluates the target at `timee`, which `tstep_update`
sets to $t^{n+1}$ at the first substep, so all three RK3 substeps of a step relax towards, and
`nesting_bcpup` imposes, the same target $\tilde q(t^{n+1})$ (design §1.2). This matches
`timedep` and the driver inflow.

**Scalars are not nested in v1.** Only `u`, `v`, `w` are imposed; temperature, humidity and passive
scalars keep the profile inlet plus convective outflow treatment under the nesting BCs.

**I/O, as implemented.** Every rank opens the file read-only and reads its own hyperslab; the one
new level at a parent-interval crossing is read synchronously at that crossing. There is no rank-0
scatter and no read-ahead (design §6.3, as revised).

Zone storage uses `zone_type` (design §9.2). The relaxation update is design §1.2, verbatim:

```fortran
qstar = qm(i,j,k) + rk3coef*qp(i,j,k)
qnew  = tgt + (qstar - tgt) * exp(-w * rk3coef / tau)      ! tau==0 => qnew = tgt
qp(i,j,k) = (qnew - qm(i,j,k)) / rk3coef
```

## 8. Ownership — do not edit files you do not own

| Package | Owns |
|---|---|
| **W1-IO** | `src/modnestingio.f90` |
| **W1-PY** | `tools/python/udprep/nesting.py`, `tools/python/tests/test_nesting.py` |
| **W2-CORE** | `src/modnesting.f90` |
| **W2-INT** | `src/modglobal.f90`, `src/modstartup.f90`, `src/program.f90`, `src/modboundary.f90` |
| **W3-TEST** | `src/tests.f90`, `tests/test_suites.yml`, `tests/integration/nesting/**` |
| **D1-DUMP** | `src/modnestdump.f90`, `tests/validation/nesting/test_nestdump_tiny.py` (the parent-side zone dump, section 9) |

`CMakeLists.txt` needs no change (`GLOB_RECURSE` + `CONFIGURE_DEPENDS`).

## 9. nestdump files `nestdump.<ipx>.<ipy>.<expnr>.nc`, `nestdump_init.<ipx>.<ipy>.<expnr>.nc`

**Written by the parent** (`src/modnestdump.f90`, namelist `&NESTDUMP`, design section 6.2), read by
whatever builds a child's `nesting.inp` -- `tests/validation/nesting/caselib.NestDump` today, other
tools tomorrow, which is why the layout is specified here.  These files are **raw parent velocity**:
no flux correction, no projection; section 5's file is produced from them by the same writer that
produces it from full field dumps.

**Namelist.**

| Name | Default | Meaning |
|---|---|---|
| `lnestdump` | `.false.` | switch |
| `tnestdump` | `1.` | dump interval [s]; tested at `rk3step == 3` like `tfielddump`, first dump at `btime + tnestdump`; at or below the timestep every step is written |
| `nestdump_x0`, `nestdump_y0` | `0.` | child box origin, parent coordinates [m]; **must coincide with parent cell faces** (`xh`, `yh`) or the run aborts |
| `nestdump_xsize`, `nestdump_ysize` | `-1.` | child box size [m]; the far faces must be parent faces too |
| `nestdump_nzone` | `0` | band thickness in **parent** cells inside each lateral face of the box: the child's guard + ramp (or its stored `nzone`, whichever is wider) on the parent grid, rounded up, **plus one cell** for the tangential-slope stencil of `udprep.nesting.conservative_interpolate`; `2 * nestdump_nzone` must be below the box size |
| `nestdump_linit` | `.true.` | write the whole box once, at the first dump, for the child's initial condition |

**Index convention.** Global 1-based parent cell indices, as `itot`/`jtot` count them. `u(i)` is the
x-face `xh(i)` west of cell `i`, `v(j)` the y-face `yh(j)` south of cell `j`, `w(k)` the z-face
`zh(k)` below cell `k`.  A block covering cells `i1..i2`, `j1..j2` carries `u` on faces `i1..i2+1`,
`v` on `j1..j2+1` and `w` on `1..ktot+1` (the lid value included), i.e. the complete staggered set
of those cells.  Blocks written by different ranks overlap on their shared faces (and the four strips
of the band overlap at the box corners); overlapping values are identical, the upper face coming from
the halo exchanged immediately before the dump.

**Band file** `nestdump.<ipx>.<ipy>.<expnr>.nc`, one per rank whose subdomain meets the band; a rank
that misses it writes nothing.  For each of the four strips the rank intersects -- west and east are
`nestdump_nzone` cells wide over the box's full `j` range, south and north `nestdump_nzone` cells
deep over the full `i` range -- one block, named by the strip:

```
dimensions:
  time = UNLIMITED ;  zt = ktot ;  zm = ktot+1 ;
  xt_west = i2-i1+1 ;  xm_west = i2-i1+2 ;  yt_west = j2-j1+1 ;  ym_west = j2-j1+2 ;
  // ... _east, _south, _north for the strips this rank writes; absent otherwise
variables:
  float time(time) ;  float zt(zt), zm(zm), xt_west(xt_west), xm_west(xm_west), ... ;   // m
  float u_west(time, zt, yt_west, xm_west) ;   u_west:stagger = "xh yf zf" ;
  float v_west(time, zt, ym_west, xt_west) ;   v_west:stagger = "xf yh zf" ;
  float w_west(time, zm, yt_west, xt_west) ;   w_west:stagger = "xf yf zh" ;
  //  each carries i_start, i_end, j_start, j_end (the block's cells, global 1-based),
  //  units = "m/s", _FillValue = -999.f ;  CDL order shown, Fortran sees the reverse
// global attributes:
  :udales_nestdump_schema = 1 ;  :itot, :jtot, :ktot, :dx, :dy ;
  :box_x0, :box_y0, :box_xsize, :box_ysize ;                 // the namelist box
  :box_i_start, :box_i_end, :box_j_start, :box_j_end, :box_ni, :box_nj ;   // the box in cells
  :nzone = nestdump_nzone ;  :tnestdump ;
  :myidx, :myidy, :nprocx, :nprocy, :rank_i_start, :rank_i_end, :rank_j_start, :rank_j_end ;
  :index_convention = <the paragraph above, in one line> ;
```

**Initial-block file** `nestdump_init.<ipx>.<ipy>.<expnr>.nc`, one per rank whose subdomain meets the
box, written once at the first dump time (the instant of the band's first record, from the same
field).  Same global attributes; no time dimension, a scalar `time`; one block `u(zt, yt, xm)`,
`v(zt, ym, xt)`, `w(zm, yt, xt)` with dimensions `xt, xm, yt, ym` over the rank's part of the box.

**Assembling the box.** Allocate `u(ni+1, nj, ktot)`, `v(ni, nj+1, ktot)`, `w(ni, nj, ktot+1)` over
the box, place every block by its `i_start`/`j_start` offset from `box_i_start`/`box_j_start`, and
treat what no block covered as missing (the reader uses `NaN`, so a slab cut or a prolongation stencil
that reaches past the band fails loudly instead of taking zeros).  `time` is single precision, as is
`fielddump`'s, and the velocities are single precision converted by the netCDF library from the
solver's doubles exactly as `fielddump`'s are -- a child built from either source is bit-identical.

**Cost accounting.** The parent prints, for the first dump and every 100th, the bytes written over
all ranks and the wall time of the write calls on the slowest rank, and a run total at exit
(`nestdump: <n> dumps, <MB> MB in total (all ranks), <s> s in the write calls (slowest rank)`).

