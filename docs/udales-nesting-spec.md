# uDALES nesting — implementation contract (v1)

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
| `nest_lateral(4)` | logical | `.true.` | W, E, S, N |
| `nest_top` | logical | `.false.` | Case C — **not implemented in v1**, must error if `.true.` |
| `nest_timeinterp` | integer | `2` | 1 linear, 2 monotone cubic Hermite |
| `nest_nwall` | integer | `1` | wall erosion, cells |
| `nest_lparentgeom` | logical | `.false.` | parent resolves the child geometry |
| `nest_fluxtol` | real | `1.e-10` | abort threshold on \|Φ\| (normalised, §7) |
| `nest_lfluxassert` | logical | `.true.` | **on by default** |

`checkinitvalues` must, when `lnesting`:
1. require `ipoiss == POISS_FFT2D`, else stop;
2. require `BCxm == BCxm_nesting .or. BCym == BCym_nesting`, else stop;
3. **not** force `BCtopm = BCtopm_pressure` (unlike the `BCxm_profile`/`BCxm_driver` branches);
4. stop if `nest_top`;
5. stop if `nest_guardwidth <= 0.` or `nest_zonewidth < 0.`;
6. stop if `nest_tau <= 0.` while `nest_zonewidth > 0.` (see the `nest_tau` row above).

## 5. File format `nesting.inp.<expnr>.nc`

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
  :Conventions = "CF-1.8" ;  :udales_nesting_schema = 1 ;
  :divergence_corrected = 1 ;
  :itot = ; :jtot = ; :ktot = ; :nzone = ; :xlen = ; :ylen = ;
  :parent_model = ; :parent_dx = ; :parent_dt = ;
  :child_origin_x = ; :child_origin_y = ; :rotation_deg = 0. ;
  :created = ; :creator = ; :tool_version = ;
  :child_dt = ;        // OPTIONAL: child timestep; without it the temporal-refinement guard is skipped
```

**Index conventions.** The west slab covers child cells `i = 1..nzone` (centres) and faces
`i = 1..nzone+1`. The east slab covers centres `i = itot-nzone+1..itot` and faces
`i = itot-nzone+1..itot+1`; **its first slab index corresponds to the lowest global index**, i.e.
slab index `m` ↔ global `i = itot-nzone+m` for centres and `i = itot-nzone+m` for faces
(so `m = nzone+1` ↔ `i = itot+1`). South/north follow the same rule in `j`.

**Units** `m s-1`. Missing/NaN values are an error, not a sentinel.

**Density convention.** Fluxes are weighted by `rhobf(k)` on both the writer and the solver side,
matching DALES's `openboundary_divcorr`. uDALES's own Poisson RHS carries no density (design F1),
so the two agree exactly while `rhobf == 1`, which is always. Both sides must use the same
convention or the runtime flux assertion can fire spuriously.

## 6. `src/modnestingio.f90` — input only, no scheme knowledge

Must not `use modnesting`. Standalone and separately compilable.

```fortran
module modnestingio
  implicit none;  save;  private
  public :: nestio_open, nestio_validate, nestio_read, nestio_close, nestio_hdr, nestio_header_type

  type nestio_header_type
    integer :: schema = 0, itot = 0, jtot = 0, ktot = 0, nzone = 0, ntime = 0
    real    :: xlen = 0., ylen = 0., rotation_deg = 0.
    logical :: divergence_corrected = .false.
    real, allocatable :: time(:), xf(:), xh(:), yf(:), yh(:), zf(:), zh(:)
    real, allocatable :: rhobf(:), rhobh(:), net_volume_flux(:)
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
  subroutine nestio_validate()

  !> Read one time level of one variable, for this rank's range of the
  !! decomposed index. varname is e.g. 'u_west'. start2/count2 are 1-based in
  !! the decomposed (outermost spatial) index. buf is (n_zone_dim, n_z_dim, count2).
  subroutine nestio_read(varname, it, start2, count2, buf, ierr)
    character(len=*), intent(in)  :: varname
    integer,          intent(in)  :: it, start2, count2
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
  real function nest_flux_residual(pup, pvp, pwp, rk3coef)
```

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

`CMakeLists.txt` needs no change (`GLOB_RECURSE` + `CONFIGURE_DEPENDS`).
