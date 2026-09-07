!> \file modnesting.f90
!!  One-way nesting by velocity imposition and relaxation zones.
!!
!!  See docs/udales-nesting-design.md for the scheme and its derivation, and
!!  docs/udales-nesting-spec.md for the normative interface contract.
!!
!!  The parent data are read by modnestingio as four lateral slabs per velocity
!!  component. Each slab is stored flat, per time slot, and the zone is stored
!!  as a point list with weights (design section 9.2) rather than as three more
!!  3-D arrays.
!!
!!  \author Maarten van Reeuwijk, Imperial College London
!
!  This file is part of uDALES.
!
! uDALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! uDALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1993-2009 Delft University of Technology, Wageningen University,
! Utrecht University, KNMI
!
module modnesting

   use modnestingio, only : nestio_open, nestio_validate, nestio_read, &
                            nestio_read_block, nestio_close, nestio_hdr, nestio_tread

   implicit none
   save
   private

   public :: nesting_init, nesting_update_target, nesting_apply, &
             nesting_boundary, nesting_bcpup, nesting_stats,     &
             nesting_finalize
   ! Test hooks: exercised directly by src/tests.f90 (runmodes TEST_NESTING_*).
   public :: nest_shape_fn, nest_union, nest_stagger_coord, nest_flux_residual, &
             nest_flux_split, nest_time_interp, nest_record_end_warnings,        &
             nest_injection, nest_nsolid_zone
   ! Namelist variables: read and broadcast by modstartup.
   public :: lnesting, nestfile, nest_guardwidth, nest_zonewidth, nest_tau,   &
             nest_shape, nest_lateral, nest_top, nest_timeinterp, nest_nwall, &
             nest_lparentgeom, nest_fluxtol, nest_lfluxassert,                &
             nest_lfluxcheckall, nest_linitfromparent, nest_statint, nest_lendabort

   logical            :: lnesting         = .false.
   character(len=256) :: nestfile         = ''
   real               :: nest_guardwidth  = 0.
   real               :: nest_zonewidth   = 0.
   real               :: nest_tau         = 0.
   integer            :: nest_shape       = 1
   logical            :: nest_lateral(4)  = .true.
   logical            :: nest_top         = .false.
   !> Time interpolation of the parent data: 1 linear, 2 cubic Hermite with
   !! unlimited Catmull-Rom slopes (default). The C0 cadence experiment
   !! (2026-09) showed 2 halves the interior TKE deficit at the same parent
   !! cadence; it is flux safe because the interpolant is linear in the data
   !! (U44, design section 3.1).
   integer            :: nest_timeinterp  = 2
   integer            :: nest_nwall       = 1
   logical            :: nest_lparentgeom = .false.
   real               :: nest_fluxtol     = 1.e-10
   logical            :: nest_lfluxassert = .true.
   !> Recompute the flux residual of every stored time level from the boundary
   !! slabs at init, instead of validating the residual the writer stored
   !! (schema 2). Off by default: the recompute costs 4 x ntime reads on every
   !! perimeter rank (design section 10.6 item 3).
   logical            :: nest_lfluxcheckall = .false.
   !> Cold start only: initialise u0/um, v0/vm and w0/wm from the full-domain
   !! block of a schema 2 file instead of from prof.inp.
   logical            :: nest_linitfromparent = .false.
   !> Interval [s] between nesting_stats reports. Each report costs three
   !! MPI_ALLREDUCEs, three full-domain sweeps and seven lines of stdout, so it
   !! is throttled like statsdump rather than run every substep. < 0 (the
   !! default) means "use tstatsdump"; 0 reports every timestep. Whatever the
   !! interval, the first and the last timestep always report, and the
   !! injection accumulators einj_* sum over every substep in between.
   real               :: nest_statint     = -1.
   !> What to do when the simulation time passes the last stored parent level.
   !! .true. (default): abort -- at initialisation already, if the run's end
   !! time timee + runtime lies beyond the record, and at the crossing
   !! otherwise. .false.: warn once and freeze the boundary on the last level.
   !! One timestep of grace is allowed at the end, so a record that ends
   !! exactly at the run's end time is not an error.
   logical            :: nest_lendabort   = .true.

   !----------------------------------------------------------------- internals

   integer, parameter :: NCOMP = 3   !< 1 = u, 2 = v, 3 = w
   integer, parameter :: NFACE = 4   !< 1 = west, 2 = east, 3 = south, 4 = north
   integer, parameter :: NSLOT = 4   !< buffered parent levels: lo-1, lo, lo+1, lo+2

   character(len=1), parameter :: cmp_name(NCOMP) = (/ 'u', 'v', 'w' /)
   character(len=5), parameter :: fac_name(NFACE) = (/ 'west ', 'east ', 'south', 'north' /)

   !> One staggered component's zone: local indices, weight, flat slab offset.
   type zone_type
      integer              :: npts = 0
      integer, allocatable :: ijk(:,:)      !< (npts,3) LOCAL indices
      real,    allocatable :: w(:)          !< (npts) weight, IBM-masked and eroded
      integer, allocatable :: src(:)        !< (npts) offset into the flat slab buffer
   end type zone_type

   type(zone_type) :: zone_u, zone_v, zone_w

   !> Flat per-component slab storage: all active faces concatenated.
   type slab_type
      integer           :: n = 0
      real, allocatable :: b(:,:)           !< (n, NSLOT) buffered parent levels
      real, allocatable :: cur(:)           !< (n) time-interpolated target
   end type slab_type

   type(slab_type) :: sbuf(NCOMP)

   ! Slab geometry, per component and face. Fortran dimension order of the file
   ! is (zone, z, decomposed), see modnestingio.
   logical :: sl_on(NCOMP, NFACE) = .false. !< slab present on this rank
   integer :: sn1(NCOMP, NFACE)   = 0       !< zone-direction length
   integer :: sn2(NCOMP, NFACE)   = 0       !< vertical length
   integer :: sn3(NCOMP, NFACE)   = 0       !< decomposed-direction count on this rank
   integer :: soff(NCOMP, NFACE)  = 0       !< offset of this slab in the flat buffer
   integer :: sds(NCOMP, NFACE)   = 1       !< first decomposed index read (global, 1-based)
   integer :: sdt(NCOMP, NFACE)   = 0       !< 1 = yf, 2 = yh, 3 = xf, 4 = xh

   logical :: linit    = .false.            !< nesting_init has completed
   logical :: lface(NFACE) = .false.        !< face active (namelist and BCxm/BCym)
   integer :: nzone    = 0                  !< zone thickness in cells, from the file
   integer :: ntime    = 0                  !< number of parent time levels
   integer :: it_lo    = 0                  !< parent level bracketing the current time
   real    :: phi_last = 0.            !< last normalised flux residual, all six faces
   real    :: phi_lid_last = 0.        !< of which the lid contributed this much
   ! Kinetic energy per unit density injected by the zone forcing since the
   ! last nesting_stats REPORT, sum over zone points of q_new * dq * dV with
   ! dq the velocity change the forcing made over the substep and dV the cell
   ! volume dx*dy*dzf(k) -- units m^5 s^-2 (multiply by rho for J). Split by
   ! where it was injected: guard strip (W >= 1) and relaxation ramp
   ! (0 < W < 1). Accumulated over every substep in nesting_apply, reported
   ! and reset in nesting_stats (design section 6.4).
   real    :: einj_guard = 0., einj_relax = 0.
   real    :: tnextstat = 0.                !< time of the next nesting_stats report
   integer :: nendwarn  = 0                 !< end-of-record warnings issued (0 or 1)
   real    :: area_bnd = 0.                 !< total FLUID domain-boundary area
   real    :: area_lat = 0.                 !< FLUID area of the four LATERAL faces only
   real    :: twall0   = 0.                 !< wall clock at the end of nesting_init
   real    :: tread0   = 0.                 !< nestio_tread at the end of nesting_init
   integer :: nsolid_zone = 0               !< solid points found inside the zone

contains

   !> Called from program.f90 AFTER readinitfiles (and after its own
   !! prerequisites createmasks and calcfluidvolumes). The order is load
   !! bearing: this routine positions the parent time buffer on timee, and
   !! readinitfiles is what assigns timee -- 0 on a cold start, the restart time
   !! on a warm one (design section 9.5). There is no separate restart record:
   !! the buffer state is reconstructed from timee alone.
   !! Opens and validates the parent file, builds the weights and the three zone
   !! point lists, enforces the building-free rule, checks the stored flux
   !! residuals and loads the first parent time levels.
   subroutine nesting_init
      use mpi,       only : MPI_Wtime
      use modglobal, only : cexpnr, timee, runtime, dx, dy, xlen, ylen, tstatsdump, &
                            BCxm, BCym, BCxm_nesting, BCym_nesting
      use modmpi,    only : myid

      integer :: ierr, f
      real    :: ltot, tflux

      if (.not. lnesting) return

      if (len_trim(nestfile) == 0) nestfile = 'nesting.inp.'//cexpnr//'.nc'

      if (nest_top) call nest_abort('nest_top (case C) is not implemented in v1')
      if (nest_shape /= 1 .and. nest_shape /= 2) &
         call nest_abort('nest_shape must be 1 (raised cosine) or 2 (quintic)')
      if (nest_timeinterp /= 1 .and. nest_timeinterp /= 2) &
         call nest_abort('nest_timeinterp must be 1 (linear) or 2 (Hermite)')
      if (nest_nwall < 0) call nest_abort('nest_nwall must be >= 0')
      if (nest_guardwidth <= 0.) call nest_abort('nest_guardwidth must be > 0')
      if (nest_zonewidth < 0.) call nest_abort('nest_zonewidth must be >= 0')
      if (nest_statint < 0.) nest_statint = tstatsdump
      tnextstat = 0.
      nendwarn  = 0
      call reset_injection

      call nestio_open(trim(nestfile), ierr)
      if (ierr /= 0) call nest_abort('cannot open '//trim(nestfile))
      call nestio_validate()

      nzone = nestio_hdr%nzone
      ntime = nestio_hdr%ntime
      if (nzone < 1) call nest_abort('nzone < 1 in '//trim(nestfile))
      if (ntime < 1) call nest_abort('no time levels in '//trim(nestfile))

      ! A face is forced only if the namelist selects it AND the corresponding
      ! momentum BC is the nesting one. A nested direction must impose BOTH of
      ! its faces: the nesting BC replaces the convective outflow as well as the
      ! inflow, so a face dropped from nest_lateral would get neither and its
      ! ghost plane would never be set. checkinitvalues rejects this from the
      ! namelist; the check is repeated here for callers that bypass it.
      if ((BCxm == BCxm_nesting) .and. .not. (nest_lateral(1) .and. nest_lateral(2))) &
         call nest_abort('BCxm = BCxm_nesting needs both x faces: nest_lateral(1:2) must be .true.')
      if ((BCym == BCym_nesting) .and. .not. (nest_lateral(3) .and. nest_lateral(4))) &
         call nest_abort('BCym = BCym_nesting needs both y faces: nest_lateral(3:4) must be .true.')
      lface(1) = nest_lateral(1) .and. (BCxm == BCxm_nesting)
      lface(2) = nest_lateral(2) .and. (BCxm == BCxm_nesting)
      lface(3) = nest_lateral(3) .and. (BCym == BCym_nesting)
      lface(4) = nest_lateral(4) .and. (BCym == BCym_nesting)
      if (.not. any(lface)) call nest_abort('lnesting is on but no lateral face is active')

      ! ---- width bookkeeping, in metres and in cells (design section 1.4) ----
      ltot = nest_guardwidth + nest_zonewidth
      if (myid == 0) then
         write(*,'(a)') ' modnesting: relaxation zone'
         write(*,'(a,es12.5,a,es12.5,a)') '   L_imp = ', nest_guardwidth, &
                                          ' m, L_rel = ', nest_zonewidth, ' m'
         if (lface(1) .or. lface(2)) &
            write(*,'(a,f8.2,a,f8.2,a,f8.2,a)') '   x: ', nest_guardwidth/dx, ' + ', &
               nest_zonewidth/dx, ' = ', ltot/dx, ' cells'
         if (lface(3) .or. lface(4)) &
            write(*,'(a,f8.2,a,f8.2,a,f8.2,a)') '   y: ', nest_guardwidth/dy, ' + ', &
               nest_zonewidth/dy, ' = ', ltot/dy, ' cells'
         write(*,'(a,i0,a)') '   file zone thickness = ', nzone, ' cells'
         if (lface(1) .or. lface(2)) then
            if (ltot/dx < 6.) write(*,'(a,f8.2,a)') &
               ' modnesting: WARNING zone is only ', ltot/dx, ' cells wide in x (< 6)'
            if (ltot > 0.15*xlen) write(*,'(a,f6.2,a)') &
               ' modnesting: WARNING zone occupies ', 100.*ltot/xlen, ' % of the domain in x (> 15 %)'
         end if
         if (lface(3) .or. lface(4)) then
            if (ltot/dy < 6.) write(*,'(a,f8.2,a)') &
               ' modnesting: WARNING zone is only ', ltot/dy, ' cells wide in y (< 6)'
            if (ltot > 0.15*ylen) write(*,'(a,f6.2,a)') &
               ' modnesting: WARNING zone occupies ', 100.*ltot/ylen, ' % of the domain in y (> 15 %)'
         end if
      end if

      if (lface(1) .or. lface(2)) then
         if (ltot > nzone*dx*(1. + 1.e-12)) call nest_abort( &
            'L_imp + L_rel exceeds the zone thickness stored in the file (x)')
      end if
      if (lface(3) .or. lface(4)) then
         if (ltot > nzone*dy*(1. + 1.e-12)) call nest_abort( &
            'L_imp + L_rel exceeds the zone thickness stored in the file (y)')
      end if

      ! ---- slab layout and zone point lists ----
      call setup_slabs

      nsolid_zone = 0
      call build_zone(1, zone_u)
      call build_zone(2, zone_v)
      call build_zone(3, zone_w)
      call report_zone

      ! ---- building-free rule (design section 9.4 item 3) ----
      if (nsolid_zone > 0) then
         if (.not. nest_lparentgeom) then
            call nest_abort('solid points found inside the relaxation zone and '// &
                            'nest_lparentgeom = .false. (the zone must be building-free)')
         else if (myid == 0) then
            write(*,'(a,i0,a)') ' modnesting: WARNING ', nsolid_zone, &
               ' solid points inside the relaxation zone (allowed by nest_lparentgeom)'
         end if
      end if

      ! ---- the stored input must be divergence corrected (design section 3.3) ----
      area_bnd = fluid_boundary_area()
      area_lat = fluid_lateral_boundary_area()
      if (area_bnd <= 0.) call nest_abort('no fluid domain-boundary area found')

      if (.not. nestio_hdr%divergence_corrected .and. myid == 0) then
         write(*,'(a)') ' modnesting: WARNING the input file is not marked divergence_corrected'
      end if

      tflux = MPI_Wtime()
      call check_stored_flux
      tflux = MPI_Wtime() - tflux
      if (myid == 0) write(*,'(a,es12.4,a)') ' modnesting: flux check took ', tflux, ' s'

      ! ---- the record must cover the run (design section 1.3) ----
      call check_record_end(timee + runtime, .true.)
      call check_record_end(timee, .false.)

      ! ---- load the first parent time levels ----
      it_lo = 0
      call set_interval(timee)
      call reload_all
      call eval_target(timee)

      ! ---- optional cold-start initialisation from the parent ----
      call init_from_parent

      linit  = .true.
      twall0 = MPI_Wtime()
      tread0 = nestio_tread

      if (myid == 0) then
         write(*,'(a,i0,a,i0,a)') ' modnesting: initialised, ', ntime, &
            ' parent time levels, buffer at interval ', it_lo, &
            ' (see docs/udales-nesting-design.md)'
      end if

      do f = 1, NFACE
         if (myid == 0 .and. lface(f)) &
            write(*,'(a,a,a)') ' modnesting: face ', trim(fac_name(f)), ' is forced'
      end do

   end subroutine nesting_init


   !> Called from program.f90 after timedep, once per RK substep. Rolls the time
   !! buffer when a parent interval is crossed (prefetching one interval ahead)
   !! and evaluates the time-interpolated target (design section 1.3).
   subroutine nesting_update_target
      use modglobal, only : timee

      if (.not. lnesting) return
      if (.not. linit) return

      call check_record_end(timee, .false.)
      call set_interval(timee)
      call eval_target(timee)

   end subroutine nesting_update_target


   !> Called from program.f90 between grwdamp and poisson. Overwrites up/vp/wp
   !! in the zone; nothing may be inserted between this call and poisson.
   !! Implements the substep-implicit update of design section 1.2.
   subroutine nesting_apply
      use modglobal, only : rk3step, dt
      use modfields, only : um, up, vm, vp, wm, wp

      integer :: n, i, j, k
      real    :: rk3coef, rk3coefi, ww, tgt, qstar, qnew, fac

      if (.not. lnesting) return
      if (.not. linit) return

      if (rk3step == 0) then
         rk3coef = 1.
      else
         rk3coef = dt/(4. - real(rk3step))
      end if
      rk3coefi = 1./rk3coef

      do n = 1, zone_u%npts
         ww = zone_u%w(n)
         if (ww <= 0.) cycle
         i = zone_u%ijk(n,1); j = zone_u%ijk(n,2); k = zone_u%ijk(n,3)
         tgt   = sbuf(1)%cur(zone_u%src(n))
         qstar = um(i,j,k) + rk3coef*up(i,j,k)
         fac   = relax_factor(ww, rk3coef)
         qnew  = tgt + (qstar - tgt)*fac
         up(i,j,k) = (qnew - um(i,j,k))*rk3coefi
         call accum_injection(ww, qnew, qnew - qstar, k)
      end do

      do n = 1, zone_v%npts
         ww = zone_v%w(n)
         if (ww <= 0.) cycle
         i = zone_v%ijk(n,1); j = zone_v%ijk(n,2); k = zone_v%ijk(n,3)
         tgt   = sbuf(2)%cur(zone_v%src(n))
         qstar = vm(i,j,k) + rk3coef*vp(i,j,k)
         fac   = relax_factor(ww, rk3coef)
         qnew  = tgt + (qstar - tgt)*fac
         vp(i,j,k) = (qnew - vm(i,j,k))*rk3coefi
         call accum_injection(ww, qnew, qnew - qstar, k)
      end do

      do n = 1, zone_w%npts
         ww = zone_w%w(n)
         if (ww <= 0.) cycle
         i = zone_w%ijk(n,1); j = zone_w%ijk(n,2); k = zone_w%ijk(n,3)
         tgt   = sbuf(3)%cur(zone_w%src(n))
         qstar = wm(i,j,k) + rk3coef*wp(i,j,k)
         fac   = relax_factor(ww, rk3coef)
         qnew  = tgt + (qstar - tgt)*fac
         wp(i,j,k) = (qnew - wm(i,j,k))*rk3coefi
         call accum_injection(ww, qnew, qnew - qstar, k)
      end do

   end subroutine nesting_apply


   !> Zero the injection accumulators. Called at init and after every
   !! nesting_stats report, so a report covers every substep since the last one.
   subroutine reset_injection
      einj_guard = 0.
      einj_relax = 0.
   end subroutine reset_injection


   !> Accumulate the kinetic energy (per unit density) the zone forcing put
   !! in at one point over this substep, q_new * dq * dV, split by guard strip
   !! versus relaxation ramp. The cell volume is what makes the sum a volume
   !! integral on a stretched grid; without it the diagnostic weighted every
   !! level equally.
   subroutine accum_injection(ww, qnew, dq, k)
      use modglobal, only : dx, dy, dzf

      real,    intent(in) :: ww, qnew, dq
      integer, intent(in) :: k

      if (ww >= 1.) then
         einj_guard = einj_guard + qnew*dq*dx*dy*dzf(k)
      else
         einj_relax = einj_relax + qnew*dq*dx*dy*dzf(k)
      end if

   end subroutine accum_injection


   !> Test hook: this rank's injection accumulators since the last report
   !! (or the last call with lreset), and optionally zero them.
   subroutine nest_injection(eguard, erelax, lreset)
      real,    intent(out) :: eguard, erelax
      logical, intent(in)  :: lreset

      eguard = einj_guard
      erelax = einj_relax
      if (lreset) call reset_injection

   end subroutine nest_injection


   !> Called from modboundary::boundary. Fills the ghost planes of u0/um, v0/vm
   !! and w0/wm from the parent, in the xmi_driver pattern
   !! (src/modboundary.f90:720).
   !!
   !! Boundary faces and ghost values inside solids are masked to zero rather
   !! than imposed (design section 4: "mask explicitly regardless"). A face
   !! value or a ghost plane serves the first interior cell; where that cell
   !! is solid at the component's own stagger the value is 0, which is what
   !! ibmnorm holds at every other solid point. Without the mask a parent that
   !! does not resolve the child's buildings (or resolves them differently)
   !! would push its flow through the solid, and nest_flux_split -- which sums
   !! fluid faces only -- could not see the divergence source it creates.
   subroutine nesting_boundary
      use modglobal, only : ib, ie, jb, je, kb, ke, ibrank, ierank, jbrank, jerank
      use modfields, only : u0, um, v0, vm, w0, wm

      integer :: i, j, k
      real    :: uu, vv, ww

      if (.not. lnesting) return
      if (.not. linit) return

      if (lface(1) .and. ibrank) then
         do j = jb - 1, je + 1
            do k = kb, ke
               uu = facval(1, 1, 1, k, j)*bnd_mask(1, ib, j, k)
               u0(ib, j, k)     = uu
               um(ib, j, k)     = uu
               u0(ib - 1, j, k) = uu
               um(ib - 1, j, k) = uu
               vv = facval(2, 1, 1, k, j)*bnd_mask(2, ib, j, k)
               v0(ib - 1, j, k) = vv
               vm(ib - 1, j, k) = vv
            end do
            do k = kb, ke + 1
               ww = facval(3, 1, 1, k, j)*bnd_mask(3, ib, j, k)
               w0(ib - 1, j, k) = ww
               wm(ib - 1, j, k) = ww
            end do
         end do
      end if

      if (lface(2) .and. ierank) then
         do j = jb - 1, je + 1
            do k = kb, ke
               uu = facval(1, 2, nzone + 1, k, j)*bnd_mask(1, ie + 1, j, k)
               u0(ie + 1, j, k) = uu
               um(ie + 1, j, k) = uu
               vv = facval(2, 2, nzone, k, j)*bnd_mask(2, ie, j, k)
               v0(ie + 1, j, k) = vv
               vm(ie + 1, j, k) = vv
            end do
            do k = kb, ke + 1
               ww = facval(3, 2, nzone, k, j)*bnd_mask(3, ie, j, k)
               w0(ie + 1, j, k) = ww
               wm(ie + 1, j, k) = ww
            end do
         end do
      end if

      if (lface(3) .and. jbrank) then
         do i = ib - 1, ie + 1
            do k = kb, ke
               vv = facval(2, 3, 1, k, i)*bnd_mask(2, i, jb, k)
               v0(i, jb, k)     = vv
               vm(i, jb, k)     = vv
               v0(i, jb - 1, k) = vv
               vm(i, jb - 1, k) = vv
               uu = facval(1, 3, 1, k, i)*bnd_mask(1, i, jb, k)
               u0(i, jb - 1, k) = uu
               um(i, jb - 1, k) = uu
            end do
            do k = kb, ke + 1
               ww = facval(3, 3, 1, k, i)*bnd_mask(3, i, jb, k)
               w0(i, jb - 1, k) = ww
               wm(i, jb - 1, k) = ww
            end do
         end do
      end if

      if (lface(4) .and. jerank) then
         do i = ib - 1, ie + 1
            do k = kb, ke
               vv = facval(2, 4, nzone + 1, k, i)*bnd_mask(2, i, je + 1, k)
               v0(i, je + 1, k) = vv
               vm(i, je + 1, k) = vv
               uu = facval(1, 4, nzone, k, i)*bnd_mask(1, i, je, k)
               u0(i, je + 1, k) = uu
               um(i, je + 1, k) = uu
            end do
            do k = kb, ke + 1
               ww = facval(3, 4, nzone, k, i)*bnd_mask(3, i, je, k)
               w0(i, je + 1, k) = ww
               wm(i, je + 1, k) = ww
            end do
         end do
      end if

   end subroutine nesting_boundary


   !> Called from modboundary::bcpup. Imposes the parent boundary-normal velocity
   !! on the predicted velocity and asserts the flux residual. Follows the
   !! BCxm_profile / BCxm_driver pattern of src/modboundary.f90:1247-1302: the
   !! face value is set in the predicted field and the tendency is zeroed, so
   !! the face only evolves through the pressure correction.
   !!
   !! Faces inside solids are masked to zero, see nesting_boundary. The start-
   !! of-step value um (and u0) at the face is set to the same target: the
   !! projection this substep uses the target evaluated at t^{n+1}, and
   !! tstep_integrate leaves the face at um + rk3coef*up = um, so without this
   !! u0 at the face would carry the previous step's target until boundary
   !! resets it -- a boundary-cell divergence that chkdiv reported and boundary
   !! then removed, polluting divtot, the case-A symptom design section 6.4
   !! says to watch. The far faces (ie+1, je+1) are outside tstep_integrate's
   !! loop, so u0/v0 are set here as well.
   subroutine nesting_bcpup(pup, pvp, pwp, rk3coef)
      use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh, &
                            ibrank, ierank, jbrank, jerank,     &
                            BCxm, BCym, BCxm_nesting, BCym_nesting, &
                            BCtopm, BCtopm_pressure
      use modfields, only : up, vp, um, vm, u0, v0
      use modmpi,    only : myid

      real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh), intent(inout) :: pup, pvp, pwp
      real, intent(in) :: rk3coef

      integer :: i, j, k
      real    :: rk3coefi, phi_closed, uu, vv

      if (.not. lnesting) return
      if (.not. linit) return

      rk3coefi = 1./rk3coef

      if (BCxm == BCxm_nesting) then
         if (lface(1) .and. ibrank) then
            do k = kb, ke
               do j = jb - 1, je + 1
                  uu = facval(1, 1, 1, k, j)*bnd_mask(1, ib, j, k)
                  pup(ib, j, k) = uu*rk3coefi
                  up(ib, j, k)  = 0. ! u(ib) only evolves according to pressure correction
                  um(ib, j, k)  = uu
                  u0(ib, j, k)  = uu
               end do
            end do
         end if

         if (lface(2) .and. ierank) then
            do k = kb, ke
               do j = jb - 1, je + 1
                  uu = facval(1, 2, nzone + 1, k, j)*bnd_mask(1, ie + 1, j, k)
                  pup(ie + 1, j, k) = uu*rk3coefi
                  up(ie + 1, j, k)  = 0.
                  um(ie + 1, j, k)  = uu
                  u0(ie + 1, j, k)  = uu
               end do
            end do
         end if
      end if

      if (BCym == BCym_nesting) then
         if (lface(3) .and. jbrank) then
            do k = kb, ke
               do i = ib - 1, ie + 1
                  vv = facval(2, 3, 1, k, i)*bnd_mask(2, i, jb, k)
                  pvp(i, jb, k) = vv*rk3coefi
                  vp(i, jb, k)  = 0.
                  vm(i, jb, k)  = vv
                  v0(i, jb, k)  = vv
               end do
            end do
         end if

         if (lface(4) .and. jerank) then
            do k = kb, ke
               do i = ib - 1, ie + 1
                  vv = facval(2, 4, nzone + 1, k, i)*bnd_mask(2, i, je + 1, k)
                  pvp(i, je + 1, k) = vv*rk3coefi
                  vp(i, je + 1, k)  = 0.
                  vm(i, je + 1, k)  = vv
                  v0(i, je + 1, k)  = vv
               end do
            end do
         end if
      end if

      ! The compatibility condition is a property of the boundary faces only
      ! (design section 3.1) and the solver will not complain on its own (F3).
      ! Under a leaky lid the top face is a free response rather than an imposed
      ! datum, so only the closed faces are asserted on - see nest_flux_split.
      ! The split is one MPI_ALLREDUCE per substep, so with the assertion off
      ! it is only evaluated on the substep whose report is about to print it.
      if (nest_lfluxassert .or. stats_due()) &
         call nest_flux_split(pup, pvp, pwp, rk3coef, phi_last, phi_lid_last)

      if (.not. nest_lfluxassert) return

      if (BCtopm == BCtopm_pressure) then
         phi_closed = phi_last - phi_lid_last
      else
         phi_closed = phi_last
      end if

      if (abs(phi_closed) > nest_fluxtol) then
         if (myid == 0) then
            write(*,'(a,es12.5,a,es12.5)') ' modnesting: flux residual ', phi_closed, &
               ' exceeds nest_fluxtol = ', nest_fluxtol
            if (BCtopm == BCtopm_pressure) write(*,'(a,es12.5,a,es12.5)') &
               '   (total over all six faces ', phi_last, ', of which the lid carries ', &
               phi_lid_last
         end if
         call nest_abort('normalised boundary flux residual out of tolerance')
      end if

   end subroutine nesting_bcpup


   !> Called from program.f90 alongside statsdump, i.e. every RK3 substep;
   !! returns at once unless a report is due (stats_due: third substep, and
   !! either nest_statint has elapsed or this is the last timestep). Reports
   !! the flux residual, the zone misfit, the pressure-gradient ratio, the
   !! energy injected since the previous report and the parent-file read cost
   !! (design section 6.4).
   subroutine nesting_stats(p)
      use mpi,       only : MPI_Wtime
      use modglobal, only : timee, ib, ie, ih, jb, je, jh, kb, ke, kh
      use modfields, only : u0, v0, w0
      use modmpi,    only : myid, comm3d, mpierr, my_real, mpi_sum

      !> Pressure increment from the last projection. Passed in from program.f90
      !! rather than taken with a use statement: modnesting cannot use modpois
      !! without closing the cycle modnesting -> modpois -> modboundary ->
      !! modnesting.
      real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh), intent(in) :: p

      real    :: sl(2), sg(2), rmsmis, twall, frac
      real    :: gl(4), gg(4), gzone, gint, gratio, el(2), eg(2)

      if (.not. lnesting) return
      if (.not. linit) return
      if (.not. stats_due()) return

      sl = 0.
      call accum_misfit(zone_u, 1, u0, sl)
      call accum_misfit(zone_v, 2, v0, sl)
      call accum_misfit(zone_w, 3, w0, sl)

      call MPI_ALLREDUCE(sl, sg, 2, MY_REAL, MPI_SUM, comm3d, mpierr)

      if (sg(2) > 0.) then
         rmsmis = sqrt(sg(1)/sg(2))
      else
         rmsmis = 0.
      end if

      ! Pressure-gradient magnitude in the zone versus the trusted interior.
      ! This is the headline diagnostic for concern C1 (design section 7): a
      ! large ratio means the projection is working hard against the imposed
      ! field and the zone is contaminating the interior through the pressure.
      gl = 0.
      call accum_gradp(zone_u, 1, p, gl)
      call accum_gradp(zone_v, 2, p, gl)
      call accum_gradp(zone_w, 3, p, gl)
      call MPI_ALLREDUCE(gl, gg, 4, MY_REAL, MPI_SUM, comm3d, mpierr)

      gzone = 0.; gint = 0.
      if (gg(2) > 0.) gzone = sqrt(gg(1)/gg(2))
      if (gg(4) > 0.) gint  = sqrt(gg(3)/gg(4))
      if (gint > 0.) then
         gratio = gzone/gint
      else
         gratio = 0.
      end if

      el(1) = einj_guard
      el(2) = einj_relax
      call MPI_ALLREDUCE(el, eg, 2, MY_REAL, MPI_SUM, comm3d, mpierr)

      twall = MPI_Wtime() - twall0
      if (twall > 0.) then
         frac = 100.*(nestio_tread - tread0)/twall
      else
         frac = 0.
      end if

      if (myid == 0) then
         write(*,'(a,f12.3)')  ' modnesting: t          = ', timee
         write(*,'(a,es12.4)') ' modnesting: Phi (norm) = ', phi_last
         write(*,'(a,es12.4,a,es12.4)') ' modnesting: Phi lid    = ', phi_lid_last, &
            '  closed faces = ', phi_last - phi_lid_last
         write(*,'(a,es12.4)') ' modnesting: zone misfit rms [m/s] = ', rmsmis
         write(*,'(a,es12.4,a,es12.4,a,f8.3)') ' modnesting: |grad p| zone = ', gzone, &
            '  interior = ', gint, '  ratio = ', gratio
         write(*,'(a,es12.4,a,es12.4,a)') ' modnesting: energy injected guard = ', eg(1), &
            '  relaxation = ', eg(2), '  [m5 s-2 per unit density, since the last report]'
         write(*,'(a,es12.4,a,f6.2,a)') ' modnesting: read time = ', nestio_tread - tread0, &
            ' s (', frac, ' % of run)'
         if (frac > 1.) write(*,'(a)') ' modnesting: WARNING parent I/O exceeds 1 % of runtime'
      end if

      call reset_injection
      if (nest_statint > 0.) then
         ! next multiple of the interval, so a warm start does not replay
         ! every report it "missed" and the cadence is the same on every run
         tnextstat = (aint(timee/nest_statint) + 1.)*nest_statint
      else
         tnextstat = timee
      end if

   end subroutine nesting_stats


   !> .true. on the substep whose nesting_stats call will report: the third
   !! RK3 substep, when the report interval has elapsed or the run is on its
   !! last timestep (timeleft is decremented at the first substep, so it is
   !! already <= 0 on every substep of that step). nesting_bcpup uses the same
   !! test to decide whether the flux split is needed for the report.
   logical function stats_due()
      use modglobal, only : rk3step, timee, timeleft

      stats_due = (rk3step == 3) .and. ((timee >= tnextstat) .or. (timeleft <= 0.))

   end function stats_due


   !> Accumulate sum(|grad p|^2) and the point count, separately over this
   !! component's zone points (acc 1,2) and over its trusted-interior points
   !! (acc 3,4), each evaluated at the component's own staggered location.
   subroutine accum_gradp(zn, ivar, p, acc)
      use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh

      type(zone_type), intent(in)    :: zn
      integer,         intent(in)    :: ivar
      real,            intent(in)    :: p(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh)
      real,            intent(inout) :: acc(4)

      integer :: n, i, j, k
      real    :: g
      logical, allocatable :: lz(:,:,:)

      ! Zone points.
      do n = 1, zn%npts
         if (zn%w(n) <= 0.) cycle
         i = zn%ijk(n,1); j = zn%ijk(n,2); k = zn%ijk(n,3)
         g = gradp_at(ivar, i, j, k, p)
         acc(1) = acc(1) + g*g
         acc(2) = acc(2) + 1.
      end do

      ! Trusted interior: every point of this rank the zone does not touch.
      ! Mark the zone points in a scratch mask first -- a per-point search of
      ! the point list would be O(npts) per interior point.
      allocate(lz(ib:ie, jb:je, kb:ke))
      lz = .false.
      do n = 1, zn%npts
         if (zn%w(n) <= 0.) cycle
         i = zn%ijk(n,1); j = zn%ijk(n,2); k = zn%ijk(n,3)
         if (i >= ib .and. i <= ie .and. j >= jb .and. j <= je .and. &
             k >= kb .and. k <= ke) lz(i,j,k) = .true.
      end do

      do k = kb, ke
         do j = jb, je
            do i = ib, ie
               if (lz(i,j,k)) cycle
               g = gradp_at(ivar, i, j, k, p)
               acc(3) = acc(3) + g*g
               acc(4) = acc(4) + 1.
            end do
         end do
      end do
      deallocate(lz)

   end subroutine accum_gradp


   !> Pressure gradient component along ivar's normal direction, at ivar's
   !! own staggered location.
   real function gradp_at(ivar, i, j, k, p)
      use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh, dxi, dyi, dzhi

      integer, intent(in) :: ivar, i, j, k
      real,    intent(in) :: p(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh)

      select case (ivar)
      case (1)
         gradp_at = (p(i,j,k) - p(i-1,j,k))*dxi
      case (2)
         gradp_at = (p(i,j,k) - p(i,j-1,k))*dyi
      case default
         gradp_at = (p(i,j,k) - p(i,j,k-1))*dzhi(max(k, kb+1))
      end select

   end function gradp_at



   !> Called from program.f90 at the end of the run (and by the unit tests
   !! between re-initialisations). Closes the parent file and releases the
   !! zone lists and slab buffers.
   subroutine nesting_finalize
      integer :: c

      if (.not. lnesting) return

      call nestio_close()

      if (allocated(zone_u%ijk)) deallocate(zone_u%ijk, zone_u%w, zone_u%src)
      if (allocated(zone_v%ijk)) deallocate(zone_v%ijk, zone_v%w, zone_v%src)
      if (allocated(zone_w%ijk)) deallocate(zone_w%ijk, zone_w%w, zone_w%src)
      zone_u%npts = 0; zone_v%npts = 0; zone_w%npts = 0

      do c = 1, NCOMP
         if (allocated(sbuf(c)%b))   deallocate(sbuf(c)%b)
         if (allocated(sbuf(c)%cur)) deallocate(sbuf(c)%cur)
         sbuf(c)%n = 0
      end do

      linit = .false.

   end subroutine nesting_finalize


   !> Shape function W(s): 1 for s<=limp, 0 for s>=limp+lrel, smooth between
   !! with zero slope at both ends. ishape: 1 raised cosine (C1), 2 quintic (C2).
   real function nest_shape_fn(s, limp, lrel, ishape)
      use modglobal, only : pi

      real,    intent(in) :: s, limp, lrel
      integer, intent(in) :: ishape

      real :: xi

      if (s <= limp) then
         nest_shape_fn = 1.
      else if (lrel <= 0.) then
         nest_shape_fn = 0.
      else if (s >= limp + lrel) then
         nest_shape_fn = 0.
      else
         xi = (s - limp)/lrel
         select case (ishape)
         case (2)
            ! quintic: C2, zero first and second derivative at both ends
            nest_shape_fn = 1. - (6.*xi**5 - 15.*xi**4 + 10.*xi**3)
         case default
            ! raised cosine: C1, zero slope at both ends
            nest_shape_fn = 0.5*(1. + cos(pi*xi))
         end select
      end if

   end function nest_shape_fn


   !> Bounded union of face weights: W = 1 - product(1 - w(1:n)).
   real function nest_union(w, n)
      integer, intent(in) :: n
      real,    intent(in) :: w(n)

      integer :: i
      real    :: prod

      prod = 1.
      do i = 1, n
         prod = prod*(1. - w(i))
      end do

      nest_union = 1. - prod

   end function nest_union


   !> Physical coordinate of staggered variable ivar (1=u, 2=v, 3=w) at LOCAL indices i,j,k.
   !! u -> (xh(ig), yf(jg), zf(k)); v -> (xf(ig), yh(jg), zf(k)); w -> (xf(ig), yf(jg), zh(k)),
   !! with ig = i + zstart(1) - 1 and jg = j + zstart(2) - 1.
   subroutine nest_stagger_coord(ivar, i, j, k, x, y, z)
      use modglobal, only : xf, xh, yf, yh, zf, zh
      use decomp_2d, only : zstart

      integer, intent(in)  :: ivar, i, j, k
      real,    intent(out) :: x, y, z

      integer :: ig, jg

      ig = i + zstart(1) - 1
      jg = j + zstart(2) - 1

      select case (ivar)
      case (1)
         x = xh(ig); y = yf(jg); z = zf(k)
      case (2)
         x = xf(ig); y = yh(jg); z = zf(k)
      case (3)
         x = xf(ig); y = yf(jg); z = zh(k)
      case default
         x = 0.; y = 0.; z = 0.
      end select

   end subroutine nest_stagger_coord


   !> Net volume flux of the predicted velocity through the domain boundary over
   !! fluid faces only, MPI-reduced, normalised by total fluid boundary area.
   !! The predicted velocity is rk3coef*p?p (see fillps, modpois.f90:942).
   !! This is Phi of design section 3.1, over all six faces; nest_flux_split
   !! additionally reports how much of it the lid carries.
   real function nest_flux_residual(pup, pvp, pwp, rk3coef)
      use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh

      real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh), intent(in) :: pup, pvp, pwp
      real, intent(in) :: rk3coef

      real :: phi_lid

      call nest_flux_split(pup, pvp, pwp, rk3coef, nest_flux_residual, phi_lid)

   end function nest_flux_residual


   !> Phi (design section 3.1) split into the total over all six domain-boundary
   !! faces and the part the LID carries, both normalised by the total fluid
   !! boundary area so the tolerance is dimensionless.
   !!
   !! The split is what makes the flux assertion correct under a leaky lid
   !! (BCtopm_pressure, design case B). There the top face is not a datum: bcpup
   !! sets w* at ke+1 from the accumulated pressure and tderive adds the matching
   !! increment, which is exactly the Dirichlet-in-the-mean-mode row the solver
   !! pins (design F3), so the projection is complete and Phi_total = 0 is NOT a
   !! requirement -- the lid flux is the child breathing against its reservoir.
   !! What must still vanish is the flux through the faces the scheme controls,
   !! phi_all - phi_lid, and that is what nesting_bcpup asserts on. Under a rigid
   !! lid (freeslip/noslip) bcpup forces w* = 0 there, phi_lid is identically
   !! zero, and the two are the same number.
   subroutine nest_flux_split(pup, pvp, pwp, rk3coef, phi_all, phi_lid)
      use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh, dx, dy, dzf, &
                            ibrank, ierank, jbrank, jerank
      use modfields, only : IIu, IIv, IIw, rhobf, rhobh
      use modmpi,    only : comm3d, mpierr, my_real, mpi_sum

      real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh), intent(in) :: pup, pvp, pwp
      real, intent(in)  :: rk3coef
      real, intent(out) :: phi_all, phi_lid

      integer :: i, j, k
      real    :: sl(3), sg(3), af

      sl = 0.

      if (ibrank) then
         do k = kb, ke
            af = dy*dzf(k)
            do j = jb, je
               if (IIu(ib, j, k) == 1) then
                  sl(1) = sl(1) - rhobf(k)*pup(ib, j, k)*rk3coef*af
                  sl(2) = sl(2) + af
               end if
            end do
         end do
      end if

      if (ierank) then
         do k = kb, ke
            af = dy*dzf(k)
            do j = jb, je
               if (IIu(ie + 1, j, k) == 1) then
                  sl(1) = sl(1) + rhobf(k)*pup(ie + 1, j, k)*rk3coef*af
                  sl(2) = sl(2) + af
               end if
            end do
         end do
      end if

      if (jbrank) then
         do k = kb, ke
            af = dx*dzf(k)
            do i = ib, ie
               if (IIv(i, jb, k) == 1) then
                  sl(1) = sl(1) - rhobf(k)*pvp(i, jb, k)*rk3coef*af
                  sl(2) = sl(2) + af
               end if
            end do
         end do
      end if

      if (jerank) then
         do k = kb, ke
            af = dx*dzf(k)
            do i = ib, ie
               if (IIv(i, je + 1, k) == 1) then
                  sl(1) = sl(1) + rhobf(k)*pvp(i, je + 1, k)*rk3coef*af
                  sl(2) = sl(2) + af
               end if
            end do
         end do
      end if

      af = dx*dy
      do j = jb, je
         do i = ib, ie
            if (IIw(i, j, kb) == 1) then
               sl(1) = sl(1) - rhobh(kb)*pwp(i, j, kb)*rk3coef*af
               sl(2) = sl(2) + af
            end if
            if (IIw(i, j, ke + 1) == 1) then
               sl(1) = sl(1) + rhobh(ke + 1)*pwp(i, j, ke + 1)*rk3coef*af
               sl(2) = sl(2) + af
               sl(3) = sl(3) + rhobh(ke + 1)*pwp(i, j, ke + 1)*rk3coef*af
            end if
         end do
      end do

      call MPI_ALLREDUCE(sl, sg, 3, MY_REAL, MPI_SUM, comm3d, mpierr)

      if (sg(2) > 0.) then
         phi_all = sg(1)/sg(2)
         phi_lid = sg(3)/sg(2)
      else
         phi_all = 0.
         phi_lid = 0.
      end if

   end subroutine nest_flux_split


   ! ------------------------------------------------------------------ private

   !> Report and abort the whole job. Rank 0 prints the message; the abort
   !! itself goes through MPI_Abort rather than `stop`, because most callers
   !! are collective but not all are: a rank-local read failure (read_level,
   !! facval) would otherwise leave the other ranks blocked in their next
   !! reduction with the job still occupying its allocation. MPI_Abort tears
   !! every rank down, whichever one called it.
   subroutine nest_abort(message)
      use mpi,    only : MPI_COMM_WORLD, MPI_Abort
      use modmpi, only : myid

      character(len=*), intent(in) :: message

      integer :: ierr

      if (myid == 0) write(*,'(a,a)') ' modnesting: ERROR ', trim(message)
      flush(6)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      stop 1   ! not reached; keeps the compiler's flow analysis honest

   end subroutine nest_abort


   !> Past the last stored parent level read_level and eval_target clamp, so
   !! the boundary silently freezes on that level -- a different problem from
   !! the one the user set up. Under nest_lendabort (the default) that is
   !! fatal; otherwise it is reported once and the run goes on frozen. One
   !! timestep (dtmax: dt itself is undefined until the first tstep_update)
   !! of grace is allowed, so a record ending exactly at the run's end time,
   !! which the last step may overshoot by up to dt, is not an error.
   !! lend: t is the run's END time, checked at init so a run that will
   !! outlast its record fails at once rather than hours in.
   subroutine check_record_end(t, lend)
      use modglobal, only : dtmax
      use modmpi,    only : myid

      real,    intent(in) :: t
      logical, intent(in) :: lend

      real :: tlast, grace

      if (ntime < 2) return   ! a single level is a steady parent, valid for all t

      tlast = nestio_hdr%time(ntime)
      grace = max(dtmax, 1.e-9*max(abs(tlast), 1.))
      if (t <= tlast + grace) return

      if (nest_lendabort) then
         if (myid == 0) then
            if (lend) then
               write(*,'(a,es12.5,a,es12.5)') ' modnesting: the run ends at t = ', t, &
                  ' but the parent record ends at t = ', tlast
            else
               write(*,'(a,es12.5,a,es12.5)') ' modnesting: t = ', t, &
                  ' is past the last parent time level at t = ', tlast
            end if
         end if
         call nest_abort('the run extends past the end of the parent record'// &
            ' (nest_lendabort = .false. freezes the boundary on the last level instead)')
      end if

      if (nendwarn == 0) then
         nendwarn = 1
         if (myid == 0) then
            if (lend) then
               write(*,'(a,es12.5,a,es12.5,a)') ' modnesting: WARNING the run ends at t = ', t, &
                  ' but the parent record ends at t = ', tlast, &
                  '; the boundary will freeze on the last level (nest_lendabort = .false.)'
            else
               write(*,'(a,es12.5,a,es12.5,a)') ' modnesting: WARNING t = ', t, &
                  ' is past the last parent time level at t = ', tlast, &
                  '; the boundary now freezes on that level (nest_lendabort = .false.)'
            end if
         end if
      end if

   end subroutine check_record_end


   !> Test hook: the number of solid points (all three staggers, MPI-reduced)
   !! nesting_init found inside the zone -- what the building-free rule
   !! judges and what the nest_lparentgeom warning reports.
   integer function nest_nsolid_zone()

      nest_nsolid_zone = nsolid_zone

   end function nest_nsolid_zone


   !> Test hook: how many end-of-record warnings have been issued since
   !! nesting_init (0 or 1 -- the warning is issued once).
   integer function nest_record_end_warnings()

      nest_record_end_warnings = nendwarn

   end function nest_record_end_warnings


   !> exp(-W*dt_s/tau), with tau <= 0 meaning Dirichlet (design section 1.2).
   real function relax_factor(ww, rk3coef)
      real, intent(in) :: ww, rk3coef

      real :: arg

      if (nest_tau <= 0.) then
         relax_factor = 0.
      else
         arg = ww*rk3coef/nest_tau
         if (arg > 700.) then
            relax_factor = 0.
         else
            relax_factor = exp(-arg)
         end if
      end if

   end function relax_factor


   !> Total FLUID area of the domain boundary, summed over the same faces as
   !! nest_flux_residual. Static, so computed once at init.
   real function fluid_boundary_area()
      use modglobal, only : ib, ie, jb, je, kb, ke, dx, dy, dzf, &
                            ibrank, ierank, jbrank, jerank
      use modfields, only : IIu, IIv, IIw
      use modmpi,    only : comm3d, mpierr, my_real, mpi_sum

      integer :: i, j, k
      ! rank-1 of size 1, not scalar: see check_stored_flux -- gfortran's
      ! `use mpi` allows only one rank per MPI_ALLREDUCE argument per file.
      real    :: al(1), ag(1), af

      al(1) = 0.

      if (ibrank) then
         do k = kb, ke
            af = dy*dzf(k)
            do j = jb, je
               if (IIu(ib, j, k) == 1) al(1) = al(1) + af
            end do
         end do
      end if
      if (ierank) then
         do k = kb, ke
            af = dy*dzf(k)
            do j = jb, je
               if (IIu(ie + 1, j, k) == 1) al(1) = al(1) + af
            end do
         end do
      end if
      if (jbrank) then
         do k = kb, ke
            af = dx*dzf(k)
            do i = ib, ie
               if (IIv(i, jb, k) == 1) al(1) = al(1) + af
            end do
         end do
      end if
      if (jerank) then
         do k = kb, ke
            af = dx*dzf(k)
            do i = ib, ie
               if (IIv(i, je + 1, k) == 1) al(1) = al(1) + af
            end do
         end do
      end if
      af = dx*dy
      do j = jb, je
         do i = ib, ie
            if (IIw(i, j, kb) == 1)     al(1) = al(1) + af
            if (IIw(i, j, ke + 1) == 1) al(1) = al(1) + af
         end do
      end do

      call MPI_ALLREDUCE(al, ag, 1, MY_REAL, MPI_SUM, comm3d, mpierr)

      fluid_boundary_area = ag(1)

   end function fluid_boundary_area


   !> FLUID area of the four LATERAL domain-boundary faces, geometric (no
   !! density). Compared against the fluid_lateral_area the writer stored, so a
   !! mask mismatch between writer and solver is caught rather than trusted.
   real function fluid_lateral_boundary_area()
      use modglobal, only : ib, ie, jb, je, kb, ke, dx, dy, dzf, &
                            ibrank, ierank, jbrank, jerank
      use modfields, only : IIu, IIv
      use modmpi,    only : comm3d, mpierr, my_real, mpi_sum

      integer :: i, j, k
      ! rank-1 of size 1, not scalar: see check_stored_flux -- gfortran's
      ! `use mpi` allows only one rank per MPI_ALLREDUCE argument per file.
      real    :: al(1), ag(1), af

      al(1) = 0.

      if (ibrank) then
         do k = kb, ke
            af = dy*dzf(k)
            do j = jb, je
               if (IIu(ib, j, k) == 1) al(1) = al(1) + af
            end do
         end do
      end if
      if (ierank) then
         do k = kb, ke
            af = dy*dzf(k)
            do j = jb, je
               if (IIu(ie + 1, j, k) == 1) al(1) = al(1) + af
            end do
         end do
      end if
      if (jbrank) then
         do k = kb, ke
            af = dx*dzf(k)
            do i = ib, ie
               if (IIv(i, jb, k) == 1) al(1) = al(1) + af
            end do
         end do
      end if
      if (jerank) then
         do k = kb, ke
            af = dx*dzf(k)
            do i = ib, ie
               if (IIv(i, je + 1, k) == 1) al(1) = al(1) + af
            end do
         end do
      end if

      call MPI_ALLREDUCE(al, ag, 1, MY_REAL, MPI_SUM, comm3d, mpierr)

      fluid_lateral_boundary_area = ag(1)

   end function fluid_lateral_boundary_area


   !> Slab dimensions and this rank's hyperslab of the decomposed index.
   !! The decomposed range is extended by one on each side (clamped to the file
   !! dimension) so the ghost planes at ib-1/ie+1/jb-1/je+1 can be served.
   subroutine setup_slabs
      use modglobal, only : itot, jtot, ktot
      use decomp_2d, only : zstart, zend

      integer :: c, f, n1, n2, n3, idt, glo, ghi, dmax, ds, de, ntot(NCOMP)

      sl_on = .false.
      sn1 = 0; sn2 = 0; sn3 = 0; soff = 0; sds = 1; sdt = 0
      ntot = 0

      do f = 1, NFACE
         if (.not. lface(f)) cycle

         ! does this rank overlap the slab's global index range?
         select case (f)
         case (1)
            if (zstart(1) > nzone + 1) cycle
         case (2)
            if (zend(1) < itot - nzone) cycle
         case (3)
            if (zstart(2) > nzone + 1) cycle
         case (4)
            if (zend(2) < jtot - nzone) cycle
         end select

         do c = 1, NCOMP
            if (f <= 2) then
               select case (c)
               case (1)
                  n1 = nzone + 1; n2 = ktot;     idt = 1
               case (2)
                  n1 = nzone;     n2 = ktot;     idt = 2
               case default
                  n1 = nzone;     n2 = ktot + 1; idt = 1
               end select
            else
               select case (c)
               case (1)
                  n1 = nzone;     n2 = ktot;     idt = 4
               case (2)
                  n1 = nzone + 1; n2 = ktot;     idt = 3
               case default
                  n1 = nzone;     n2 = ktot + 1; idt = 3
               end select
            end if

            select case (idt)
            case (1)
               glo = zstart(2); ghi = zend(2);     dmax = jtot
            case (2)
               glo = zstart(2); ghi = zend(2) + 1; dmax = jtot + 1
            case (3)
               glo = zstart(1); ghi = zend(1);     dmax = itot
            case default
               glo = zstart(1); ghi = zend(1) + 1; dmax = itot + 1
            end select

            ds = max(1, glo - 1)
            de = min(dmax, ghi + 1)
            n3 = de - ds + 1
            if (n3 < 1) cycle

            sl_on(c,f) = .true.
            sn1(c,f)   = n1
            sn2(c,f)   = n2
            sn3(c,f)   = n3
            sds(c,f)   = ds
            sdt(c,f)   = idt
            soff(c,f)  = ntot(c)
            ntot(c)    = ntot(c) + n1*n2*n3
         end do
      end do

      do c = 1, NCOMP
         sbuf(c)%n = ntot(c)
         if (allocated(sbuf(c)%b))   deallocate(sbuf(c)%b)
         if (allocated(sbuf(c)%cur)) deallocate(sbuf(c)%cur)
         allocate(sbuf(c)%b(max(ntot(c),1), NSLOT)); sbuf(c)%b = 0.
         allocate(sbuf(c)%cur(max(ntot(c),1)));      sbuf(c)%cur = 0.
      end do

   end subroutine setup_slabs


   !> Flat index of slab element (m, kk, d) of component c, face f.
   integer function sidx(c, f, m, kk, d)
      integer, intent(in) :: c, f, m, kk, d

      sidx = soff(c,f) + m + (kk - 1)*sn1(c,f) + (d - 1)*sn1(c,f)*sn2(c,f)

   end function sidx


   !> Local index (i or j, depending on the slab's decomposed direction) mapped
   !! to this rank's slab index, clamped at the domain edges.
   integer function dloc(c, f, l)
      use decomp_2d, only : zstart

      integer, intent(in) :: c, f, l

      integer :: g

      select case (sdt(c,f))
      case (1, 2)
         g = l + zstart(2) - 1
      case default
         g = l + zstart(1) - 1
      end select

      dloc = min(max(g - sds(c,f) + 1, 1), sn3(c,f))

   end function dloc


   !> Time-interpolated parent value on a boundary plane: component c, face f,
   !! zone index m, vertical index kk, local index l along the face.
   real function facval(c, f, m, kk, l)
      integer, intent(in) :: c, f, m, kk, l

      integer :: kc

      ! Every rank that imposes a face has that face's slab (setup_slabs marks
      ! it for every rank overlapping the slab's index range, and the face
      ! ranks always do). Reaching here without it is an inconsistency between
      ! the slab layout and the face bookkeeping, not a value to return 0 for.
      if (.not. sl_on(c,f)) call nest_abort('facval: component '//cmp_name(c)// &
         ' of face '//trim(fac_name(f))//' is imposed on a rank that holds no slab for it')

      kc = min(max(kk, 1), sn2(c,f))
      facval = sbuf(c)%cur(sidx(c, f, m, kc, dloc(c, f, l)))

   end function facval


   !> IBM mask at the stagger of ivar (1 = fluid, 0 = solid).
   integer function iimask(ivar, i, j, k)
      use modfields, only : IIu, IIv, IIw

      integer, intent(in) :: ivar, i, j, k

      select case (ivar)
      case (1)
         iimask = IIu(i,j,k)
      case (2)
         iimask = IIv(i,j,k)
      case default
         iimask = IIw(i,j,k)
      end select

   end function iimask


   !> Fluid mask (1. fluid, 0. solid) of component c at the stagger point
   !! (i,j,k), for masking a domain-boundary face or ghost value: the caller
   !! passes the index of the FIRST INTERIOR point the value serves. IIu/IIv/IIw
   !! carry no valid halos, so the indices are clamped to the range that is
   !! defined for that stagger (u: i up to ie+1, v: j up to je+1, w: k up to
   !! ke+1, everything else interior). The clamp only ever acts on the ghost
   !! rows and columns at the domain corners, which enter neither the Poisson
   !! RHS nor any interior cell's stencil.
   real function bnd_mask(c, i, j, k)
      use modglobal, only : ib, ie, jb, je, kb, ke

      integer, intent(in) :: c, i, j, k

      integer :: ic, jc, kc

      ic = min(max(i, ib), ie)
      jc = min(max(j, jb), je)
      kc = min(max(k, kb), ke)
      select case (c)
      case (1)
         ic = min(max(i, ib), ie + 1)
      case (2)
         jc = min(max(j, jb), je + 1)
      case default
         kc = min(max(k, kb), ke + 1)
      end select

      bnd_mask = real(iimask(c, ic, jc, kc))

   end function bnd_mask


   !> Build the fluid mask (1 fluid, 0 solid) at the stagger of ivar, eroded by
   !! nest_nwall cells away from any solid point (design section 5 item 2).
   !! IIu/IIv/IIw have no valid halos, so a real copy is made and exchanged; the
   !! erosion is then applied one cell at a time with a halo exchange after each
   !! pass, which keeps it exact and decomposition independent for any
   !! nest_nwall. The k halo planes are left fluid: the ground is a boundary
   !! condition, not a solid point, and IIw(:,:,kb) = 0 is a convention of
   !! createmasks rather than geometry.
   subroutine build_eroded_mask(ivar, msk)
      use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh, libm
      use decomp_2d, only : exchange_halo_z

      integer,           intent(in)  :: ivar
      real, allocatable, intent(out) :: msk(:,:,:)

      real, allocatable :: prev(:,:,:)
      integer :: i, j, k, n, di, dj, dk, kk
      logical :: lsolid

      allocate(msk(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh))
      msk = 1.

      do k = kb, ke + kh
         do j = jb, je
            do i = ib, ie
               msk(i,j,k) = real(iimask(ivar, i, j, k))
            end do
         end do
      end do
      if (ivar == 3) msk(:,:,kb) = 1.

      call exchange_halo_z(msk)

      if (.not. libm) return
      if (nest_nwall < 1) return

      allocate(prev(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh))

      do n = 1, nest_nwall
         prev = msk
         do k = kb, ke + kh
            do j = jb, je
               do i = ib, ie
                  if (prev(i,j,k) < 0.5) cycle
                  lsolid = .false.
                  do dk = -1, 1
                     ! Clamp in k: the array stops at ke+kh, and there is no cell
                     ! above the lid to erode from. Reading k+dk unclamped ran one
                     ! plane past the upper bound, which silently wiped the whole
                     ! k = ke zone layer for nest_nwall >= 2 and made the zone
                     ! decomposition-dependent.
                     kk = min(max(k + dk, kb - kh), ke + kh)
                     do dj = -1, 1
                        do di = -1, 1
                           if (prev(i + di, j + dj, kk) < 0.5) lsolid = .true.
                        end do
                     end do
                  end do
                  if (lsolid) msk(i,j,k) = 0.
               end do
            end do
         end do
         call exchange_halo_z(msk)
      end do

      deallocate(prev)

   end subroutine build_eroded_mask


   !> Unmasked weight of variable ivar at LOCAL indices i,j,k, and the per-face
   !! weights used to pick the slab a zone point is read from.
   subroutine zone_weight(ivar, i, j, k, wf, ww)
      use modglobal, only : xlen, ylen

      integer, intent(in)  :: ivar, i, j, k
      real,    intent(out) :: wf(NFACE), ww

      real :: x, y, z

      call nest_stagger_coord(ivar, i, j, k, x, y, z)

      wf = 0.
      if (lface(1)) wf(1) = nest_shape_fn(x,        nest_guardwidth, nest_zonewidth, nest_shape)
      if (lface(2)) wf(2) = nest_shape_fn(xlen - x, nest_guardwidth, nest_zonewidth, nest_shape)
      if (lface(3)) wf(3) = nest_shape_fn(y,        nest_guardwidth, nest_zonewidth, nest_shape)
      if (lface(4)) wf(4) = nest_shape_fn(ylen - y, nest_guardwidth, nest_zonewidth, nest_shape)

      ww = nest_union(wf, NFACE)

   end subroutine zone_weight


   !> Zone index m of variable ivar at global indices (ig,jg) within slab face f,
   !! or -1 when the point is outside that slab.
   integer function zone_m(ivar, f, ig, jg)
      use modglobal, only : itot, jtot

      integer, intent(in) :: ivar, f, ig, jg

      integer :: g, nlen, base

      select case (f)
      case (1, 2)
         g = ig
         if (ivar == 1) then
            nlen = nzone + 1
         else
            nlen = nzone
         end if
         if (f == 1) then
            base = 0
         else
            base = itot - nzone
         end if
      case default
         g = jg
         if (ivar == 2) then
            nlen = nzone + 1
         else
            nlen = nzone
         end if
         if (f == 3) then
            base = 0
         else
            base = jtot - nzone
         end if
      end select

      zone_m = g - base
      if (zone_m < 1 .or. zone_m > nlen) zone_m = -1

   end function zone_m


   !> Build the point list, weights and slab offsets for one staggered component.
   !! Points are dropped where the IBM mask is solid or where a solid point lies
   !! within nest_nwall cells, so the trusted interior and the solid points are
   !! left bit-identical by nesting_apply.
   subroutine build_zone(ivar, zn)
      use modglobal, only : ib, ie, jb, je, kb, ke
      use modmpi,    only : comm3d, mpierr, my_real, mpi_sum
      use decomp_2d, only : zstart

      integer,         intent(in)  :: ivar
      type(zone_type), intent(out) :: zn

      integer :: i, j, k, ipass, n, f, fbest, m, ig, jg, l, nbad
      real    :: wf(NFACE), ww, wbest
      ! rank-1 of size 1, not scalar: see check_stored_flux -- gfortran's
      ! `use mpi` allows only one rank per MPI_ALLREDUCE argument per file.
      real    :: rl(1), rg(1)
      real, allocatable :: msk(:,:,:)

      nbad = 0
      call build_eroded_mask(ivar, msk)

      do ipass = 1, 2
         n = 0
         do k = kb, ke
            do j = jb, je
               do i = ib, ie
                  ! w on the ground plane is set by the bottom BC, not by nesting
                  if (ivar == 3 .and. k == kb) cycle

                  call zone_weight(ivar, i, j, k, wf, ww)
                  if (ww <= 0.) cycle

                  if (iimask(ivar, i, j, k) == 0) then
                     if (ipass == 1) nsolid_zone = nsolid_zone + 1
                     cycle
                  end if
                  if (msk(i,j,k) < 0.5) cycle

                  ! pick the strongest face whose slab actually contains the point
                  ig = i + zstart(1) - 1
                  jg = j + zstart(2) - 1
                  fbest = 0
                  wbest = -1.
                  do f = 1, NFACE
                     if (.not. sl_on(ivar,f)) cycle
                     if (zone_m(ivar, f, ig, jg) < 0) cycle
                     if (wf(f) > wbest) then
                        wbest = wf(f)
                        fbest = f
                     end if
                  end do

                  if (fbest == 0) then
                     if (ipass == 1) nbad = nbad + 1
                     cycle
                  end if

                  n = n + 1
                  if (ipass == 2) then
                     m = zone_m(ivar, fbest, ig, jg)
                     if (fbest <= 2) then
                        l = j
                     else
                        l = i
                     end if
                     zn%ijk(n,1) = i
                     zn%ijk(n,2) = j
                     zn%ijk(n,3) = k
                     zn%w(n)     = ww
                     zn%src(n)   = sidx(ivar, fbest, m, k, dloc(ivar, fbest, l))
                  end if
               end do
            end do
         end do

         if (ipass == 1) then
            zn%npts = n
            allocate(zn%ijk(max(n,1),3)); zn%ijk = 0
            allocate(zn%w(max(n,1)));     zn%w   = 0.
            allocate(zn%src(max(n,1)));   zn%src = 1
         end if
      end do

      deallocate(msk)

      rl(1) = real(nbad)
      call MPI_ALLREDUCE(rl, rg, 1, MY_REAL, MPI_SUM, comm3d, mpierr)
      if (rg(1) > 0.) call nest_abort('zone points fall outside the slabs stored in the '// &
                                   'nesting file - increase nzone or reduce the zone width')

   end subroutine build_zone


   !> Report the zone size, MPI-reduced.
   subroutine report_zone
      use modmpi, only : myid, comm3d, mpierr, my_real, mpi_sum

      real :: rl(4), rg(4)
      integer :: n

      rl(1) = real(zone_u%npts)
      rl(2) = real(zone_v%npts)
      rl(3) = real(zone_w%npts)
      rl(4) = real(nsolid_zone)

      call MPI_ALLREDUCE(rl, rg, 4, MY_REAL, MPI_SUM, comm3d, mpierr)

      nsolid_zone = nint(rg(4))

      if (myid == 0) then
         n = nint(rg(1)) + nint(rg(2)) + nint(rg(3))
         write(*,'(a,i0,a,i0,a,i0,a,i0,a)') ' modnesting: zone points u/v/w = ', &
            nint(rg(1)), '/', nint(rg(2)), '/', nint(rg(3)), ' (total ', n, ')'
      end if

   end subroutine report_zone


   !> Verify that every stored parent time level is flux balanced.
   !!
   !! Phi is a linear functional of the boundary data (design section 3.1
   !! item 2), so checking every stored level is sufficient for every
   !! time-interpolated target.
   !!
   !! A schema 2 file carries flux_residual(time), the residual of the data AS
   !! STORED, and fluid_lateral_area, the area it was summed over. When the two
   !! sides agree on that area -- i.e. the writer masked the same boundary faces
   !! the solver calls solid -- the stored residual is validated directly and no
   !! slab is read at all. Otherwise, and for a schema 1 file, which stores only
   !! the PRE-correction net_volume_flux and hence nothing about the corrected
   !! data, the residual is recomputed from the boundary slabs: 4 x ntime reads
   !! on every perimeter rank. nest_lfluxcheckall forces the recompute
   !! unconditionally (design section 10.6 item 3).
   subroutine check_stored_flux
      use modglobal, only : ib, ie, jb, je, kb, ke, dx, dy, dzf, &
                            ibrank, ierank, jbrank, jerank
      use modfields, only : IIu, IIv, rhobf
      use modmpi,    only : myid, comm3d, mpierr, my_real, mpi_sum

      integer :: it, i, j, k, ierr
      ! sl/sg are rank-1 of size 1 rather than scalars: gfortran's `use mpi`
      ! builds a single implicit interface for MPI_ALLREDUCE per file, and this
      ! file also reduces rank-1 buffers elsewhere. Mixing the two ranks is a
      ! hard error there (it compiles under ifort). Keep every MPI_ALLREDUCE
      ! actual argument in this file rank-1.
      real    :: sl(1), sg(1), phi, af, aerr
      logical :: lw, le, ls, ln, lfull
      real, allocatable :: bw(:,:,:), be(:,:,:), bs(:,:,:), bn(:,:,:)

      ! ---- decide between the cheap check and the full recompute ----
      lfull = nest_lfluxcheckall

      if (.not. nestio_hdr%has_flux_residual) then
         if (myid == 0 .and. .not. lfull) write(*,'(a)') &
            ' modnesting: WARNING the input file predates schema 2 and stores no'// &
            ' post-correction flux residual; recomputing it from the boundary slabs'
         lfull = .true.
      else
         aerr = abs(nestio_hdr%fluid_lateral_area - area_lat)
         if (aerr > 1.e-8*max(area_lat, abs(nestio_hdr%fluid_lateral_area), 1.e-30)) then
            if (myid == 0) then
               write(*,'(a)') ' modnesting: WARNING the fluid lateral boundary area of the'// &
                  ' input file does not match this run'
               write(*,'(a,es22.14,a,es22.14)') '   file = ', nestio_hdr%fluid_lateral_area, &
                  ', run = ', area_lat
               write(*,'(a)') '   the stored flux residual cannot be trusted;'// &
                  ' recomputing it from the boundary slabs'
            end if
            lfull = .true.
         end if
      end if

      if (.not. lfull) then
         do it = 1, ntime
            phi = nestio_hdr%flux_residual(it)/area_bnd
            if (abs(phi) > nest_fluxtol) then
               if (myid == 0) then
                  write(*,'(a,i0,a,es12.5,a,es12.5)') &
                     ' modnesting: stored flux residual of time level ', it, ' is ', phi, &
                     ' (normalised), tolerance ', nest_fluxtol
               end if
               call nest_abort('the parent file is not flux balanced - rerun the offline correction')
            end if
         end do
         if (myid == 0) write(*,'(a,i0,a)') ' modnesting: ', ntime, &
            ' stored time levels are flux balanced (from the stored residual, no slab read)'
         return
      end if

      lw = lface(1) .and. ibrank .and. sl_on(1,1)
      le = lface(2) .and. ierank .and. sl_on(1,2)
      ls = lface(3) .and. jbrank .and. sl_on(2,3)
      ln = lface(4) .and. jerank .and. sl_on(2,4)

      if (lw) allocate(bw(sn1(1,1), sn2(1,1), sn3(1,1)))
      if (le) allocate(be(sn1(1,2), sn2(1,2), sn3(1,2)))
      if (ls) allocate(bs(sn1(2,3), sn2(2,3), sn3(2,3)))
      if (ln) allocate(bn(sn1(2,4), sn2(2,4), sn3(2,4)))

      do it = 1, ntime
         sl(1) = 0.

         if (lw) then
            call nestio_read('u_west', it, sds(1,1), sn3(1,1), bw, ierr)
            if (ierr /= 0) call nest_abort('failed to read u_west from '//trim(nestfile))
            do k = kb, ke
               af = dy*dzf(k)
               do j = jb, je
                  if (IIu(ib,j,k) == 1) &
                     sl(1) = sl(1) - rhobf(k)*bw(1, k, dloc(1,1,j))*af
               end do
            end do
         end if

         if (le) then
            call nestio_read('u_east', it, sds(1,2), sn3(1,2), be, ierr)
            if (ierr /= 0) call nest_abort('failed to read u_east from '//trim(nestfile))
            do k = kb, ke
               af = dy*dzf(k)
               do j = jb, je
                  if (IIu(ie+1,j,k) == 1) &
                     sl(1) = sl(1) + rhobf(k)*be(nzone + 1, k, dloc(1,2,j))*af
               end do
            end do
         end if

         if (ls) then
            call nestio_read('v_south', it, sds(2,3), sn3(2,3), bs, ierr)
            if (ierr /= 0) call nest_abort('failed to read v_south from '//trim(nestfile))
            do k = kb, ke
               af = dx*dzf(k)
               do i = ib, ie
                  if (IIv(i,jb,k) == 1) &
                     sl(1) = sl(1) - rhobf(k)*bs(1, k, dloc(2,3,i))*af
               end do
            end do
         end if

         if (ln) then
            call nestio_read('v_north', it, sds(2,4), sn3(2,4), bn, ierr)
            if (ierr /= 0) call nest_abort('failed to read v_north from '//trim(nestfile))
            do k = kb, ke
               af = dx*dzf(k)
               do i = ib, ie
                  if (IIv(i,je+1,k) == 1) &
                     sl(1) = sl(1) + rhobf(k)*bn(nzone + 1, k, dloc(2,4,i))*af
               end do
            end do
         end if

         call MPI_ALLREDUCE(sl, sg, 1, MY_REAL, MPI_SUM, comm3d, mpierr)
         phi = sg(1)/area_bnd

         if (abs(phi) > nest_fluxtol) then
            if (myid == 0) then
               write(*,'(a,i0,a,es12.5,a,es12.5)') &
                  ' modnesting: flux residual of stored time level ', it, ' is ', phi, &
                  ' (normalised), tolerance ', nest_fluxtol
            end if
            call nest_abort('the parent file is not flux balanced - rerun the offline correction')
         end if
      end do

      if (allocated(bw)) deallocate(bw)
      if (allocated(be)) deallocate(be)
      if (allocated(bs)) deallocate(bs)
      if (allocated(bn)) deallocate(bn)

      if (myid == 0) write(*,'(a,i0,a)') ' modnesting: ', ntime, &
         ' stored time levels are flux balanced (recomputed from the boundary slabs)'

   end subroutine check_stored_flux


   !> Cold-start initialisation of the interior from the parent (design section
   !! 10.6 item 4). Called at the end of nesting_init, i.e. after readinitfiles
   !! has filled u0/um from prof.inp, which this then overwrites.
   !!
   !! Only a cold start is touched: on a warm start the restart file already
   !! holds a state consistent with the parent, and overwriting it would break
   !! restart parity (test I6). The stored block is the field at the FIRST
   !! stored time and has been made discretely solenoidal on the child grid by
   !! the writer, with the boundary-normal velocities equal to the slab values
   !! at that time, so the first projection has nothing to clean up.
   subroutine init_from_parent
      use modglobal, only : ib, ie, jb, je, kb, ke, kh, ktot, &
                            ierank, jerank, timee, lwarmstart, lstratstart
      use modfields, only : u0, um, v0, vm, w0, wm, u0av, v0av, &
                            IIu, IIus, IIv, IIvs
      use modmpi,    only : myid, avexy_ibm
      use decomp_2d, only : zstart, exchange_halo_z

      integer :: ierr, i, j, k, ni, nj, i0, j0
      real, allocatable :: buf(:,:,:)

      if (.not. nest_linitfromparent) return

      if (lwarmstart .or. lstratstart) then
         if (myid == 0) write(*,'(a)') ' modnesting: nest_linitfromparent is set but this'// &
            ' is a warm start; the restart file wins and the parent block is not read'
         return
      end if

      if (.not. nestio_hdr%has_initial_condition) call nest_abort( &
         'nest_linitfromparent is set but '//trim(nestfile)//' carries no initial-condition'// &
         ' block (it needs schema 2 with has_initial_condition = 1)')

      if (abs(timee - nestio_hdr%time(1)) > 1.e-8*max(abs(timee), 1.) .and. myid == 0) then
         write(*,'(a,es12.5,a,es12.5)') ' modnesting: WARNING the initial-condition block is'// &
            ' stored at t = ', nestio_hdr%time(1), ' but the run starts at t = ', timee
      end if

      i0 = zstart(1)
      j0 = zstart(2)

      ! u lives on xh: the east-most rank additionally owns the face at itot+1.
      ni = ie - ib + 1
      if (ierank) ni = ni + 1
      nj = je - jb + 1
      allocate(buf(ktot, nj, ni))
      call nestio_read_block('u_init', j0, nj, i0, ni, buf, ierr)
      if (ierr /= 0) call nest_abort('failed to read, or found non-finite or fill values in,'// &
         ' u_init of '//trim(nestfile))
      do i = 1, ni
         do j = 1, nj
            do k = 1, ktot
               u0(ib + i - 1, jb + j - 1, kb + k - 1) = buf(k, j, i)
            end do
         end do
      end do
      deallocate(buf)

      ! v lives on yh: the north-most rank additionally owns the face at jtot+1.
      ni = ie - ib + 1
      nj = je - jb + 1
      if (jerank) nj = nj + 1
      allocate(buf(ktot, nj, ni))
      call nestio_read_block('v_init', j0, nj, i0, ni, buf, ierr)
      if (ierr /= 0) call nest_abort('failed to read, or found non-finite or fill values in,'// &
         ' v_init of '//trim(nestfile))
      do i = 1, ni
         do j = 1, nj
            do k = 1, ktot
               v0(ib + i - 1, jb + j - 1, kb + k - 1) = buf(k, j, i)
            end do
         end do
      end do
      deallocate(buf)

      ! w lives on zh: every rank owns ktot+1 levels, kb .. ke+1.
      ni = ie - ib + 1
      nj = je - jb + 1
      allocate(buf(ktot + 1, nj, ni))
      call nestio_read_block('w_init', j0, nj, i0, ni, buf, ierr)
      if (ierr /= 0) call nest_abort('failed to read, or found non-finite or fill values in,'// &
         ' w_init of '//trim(nestfile))
      do i = 1, ni
         do j = 1, nj
            do k = 1, ktot + 1
               w0(ib + i - 1, jb + j - 1, kb + k - 1) = buf(k, j, i)
            end do
         end do
      end do
      deallocate(buf)

      call exchange_halo_z(u0)
      call exchange_halo_z(v0)
      call exchange_halo_z(w0)

      um = u0
      vm = v0
      wm = w0

      ! readinitfiles computed these from the prof.inp fields we have just
      ! replaced, so they would otherwise describe a state the run no longer has.
      call avexy_ibm(u0av(kb:ke+kh), u0(ib:ie,jb:je,kb:ke+kh), ib, ie, jb, je, kb, ke, kh, &
                     IIu(ib:ie,jb:je,kb:ke+kh), IIus(kb:ke+kh), .false.)
      call avexy_ibm(v0av(kb:ke+kh), v0(ib:ie,jb:je,kb:ke+kh), ib, ie, jb, je, kb, ke, kh, &
                     IIv(ib:ie,jb:je,kb:ke+kh), IIvs(kb:ke+kh), .false.)

      if (myid == 0) write(*,'(a,es12.5)') ' modnesting: cold start initialised from the'// &
         ' parent initial-condition block at t = ', nestio_hdr%time(1)

   end subroutine init_from_parent


   !> Read one parent time level into a buffer slot. Levels outside [1,ntime]
   !! are clamped, so the ends of the record replicate rather than extrapolate.
   subroutine read_level(ilev, islot)
      integer, intent(in) :: ilev, islot

      integer :: c, f, it, ierr, nn
      real, allocatable :: tmp(:,:,:)

      it = min(max(ilev, 1), ntime)

      do c = 1, NCOMP
         do f = 1, NFACE
            if (.not. sl_on(c,f)) cycle
            allocate(tmp(sn1(c,f), sn2(c,f), sn3(c,f)))
            call nestio_read(cmp_name(c)//'_'//trim(fac_name(f)), it, &
                             sds(c,f), sn3(c,f), tmp, ierr)
            if (ierr /= 0) then
               deallocate(tmp)
               call nest_abort('failed to read, or found non-finite or fill values in, slab '// &
                  cmp_name(c)//'_'//trim(fac_name(f))//' of '//trim(nestfile)// &
                  ' (see the modnestingio line above)')
            end if
            nn = sn1(c,f)*sn2(c,f)*sn3(c,f)
            sbuf(c)%b(soff(c,f) + 1:soff(c,f) + nn, islot) = reshape(tmp, (/ nn /))
            deallocate(tmp)
         end do
      end do

   end subroutine read_level


   !> Fill all NSLOT buffer slots for the current interval it_lo.
   subroutine reload_all
      integer :: s

      do s = 1, NSLOT
         call read_level(it_lo - 2 + s, s)
      end do

   end subroutine reload_all


   !> Position it_lo so that time(it_lo) <= t < time(it_lo+1) and keep the
   !! buffer in step. Crossing one interval rolls the slots and reads the level
   !! one interval ahead (design section 6.2 item 3); a larger jump reloads.
   subroutine set_interval(t)
      real, intent(in) :: t

      integer :: il, ilold, d, s

      ilold = it_lo
      il = max(min(it_lo, max(ntime - 1, 1)), 1)

      if (ntime > 1) then
         do while (il < ntime - 1)
            if (t < nestio_hdr%time(il + 1)) exit
            il = il + 1
         end do
         do while (il > 1)
            if (t >= nestio_hdr%time(il)) exit
            il = il - 1
         end do
      else
         il = 1
      end if

      it_lo = il
      if (ilold == 0) return   ! nesting_init loads the slots itself

      d = it_lo - ilold
      if (d == 0) return

      if (d > 0 .and. d < NSLOT) then
         do s = 1, NSLOT - d
            sbuf(1)%b(:,s) = sbuf(1)%b(:,s + d)
            sbuf(2)%b(:,s) = sbuf(2)%b(:,s + d)
            sbuf(3)%b(:,s) = sbuf(3)%b(:,s + d)
         end do
         do s = NSLOT - d + 1, NSLOT
            call read_level(it_lo - 2 + s, s)
         end do
      else
         call reload_all
      end if

   end subroutine set_interval


   !> Evaluate the time-interpolated target for every buffered slab value.
   !! nest_timeinterp: 1 linear, 2 cubic Hermite (Catmull-Rom, unlimited).
   subroutine eval_target(t)
      real, intent(in) :: t

      integer :: c, m
      real    :: tlo, thi, h1, h2, h3, th
      real    :: y1, y2, y3, y4

      tlo = nestio_hdr%time(min(it_lo, ntime))
      thi = nestio_hdr%time(min(it_lo + 1, ntime))

      h2 = thi - tlo
      if (h2 > 0.) then
         th = (t - tlo)/h2
      else
         th = 0.
      end if
      th = min(max(th, 0.), 1.)

      h1 = tlo - nestio_hdr%time(max(it_lo - 1, 1))
      h3 = nestio_hdr%time(min(it_lo + 2, ntime)) - thi

      do c = 1, NCOMP
         if (sbuf(c)%n <= 0) cycle
         if (nest_timeinterp == 1 .or. h2 <= 0.) then
            do m = 1, sbuf(c)%n
               y2 = sbuf(c)%b(m,2)
               y3 = sbuf(c)%b(m,3)
               sbuf(c)%cur(m) = y2 + th*(y3 - y2)
            end do
         else
            do m = 1, sbuf(c)%n
               y1 = sbuf(c)%b(m,1)
               y2 = sbuf(c)%b(m,2)
               y3 = sbuf(c)%b(m,3)
               y4 = sbuf(c)%b(m,4)
               sbuf(c)%cur(m) = hermite(y1, y2, y3, y4, h1, h2, h3, th)
            end do
         end if
      end do

   end subroutine eval_target


   !> Test hook: the time interpolant itself, for one scalar sample.
   !! mode 1 linear on [t2,t3], mode 2 the Hermite of `hermite` below. This is
   !! the SAME code the solver uses (eval_target calls the same routines), so
   !! U44 cannot pass while the solver does something else.
   real function nest_time_interp(y1, y2, y3, y4, h1, h2, h3, th, mode)
      real,    intent(in) :: y1, y2, y3, y4, h1, h2, h3, th
      integer, intent(in) :: mode

      if (mode == 1) then
         nest_time_interp = y2 + th*(y3 - y2)   ! exactly as in eval_target
      else
         nest_time_interp = hermite(y1, y2, y3, y4, h1, h2, h3, th)
      end if

   end function nest_time_interp


   !> Cubic Hermite on [t2,t3] (Catmull-Rom, generalised to non-uniform
   !! spacing). C1 across interval crossings and exact for data linear in time.
   !!
   !! The slopes are deliberately UNLIMITED. The value is a fixed linear
   !! combination of y1..y4 whose coefficients depend only on the spacings and
   !! on th -- never on the data -- so a boundary field with zero net flux at
   !! every stored level keeps zero net flux everywhere in between, which is
   !! what design section 3.1 needs. A monotone limiter (Fritsch-Carlson, used
   !! here originally) makes the coefficients data dependent: each boundary
   !! face is then weighted differently and the cancellation collapses.
   !! Measured, with every stored level corrected to |Phi| ~ 1e-14: linear
   !! 1.1e-14, this routine 1.4e-14, Fritsch-Carlson 1.0e+00. Monotonicity
   !! buys nothing for a sign-unconstrained quantity like velocity. Do not
   !! reintroduce a limiter here; U44 will fail if you do.
   real function hermite(y1, y2, y3, y4, h1, h2, h3, th)
      real, intent(in) :: y1, y2, y3, y4, h1, h2, h3, th

      real :: s1, s2, s3, d2, d3, t2, t3c, h00, h10, h01, h11

      s2 = (y3 - y2)/h2

      ! Slopes are the spacing-weighted ARITHMETIC mean of the one-sided
      ! slopes (Catmull-Rom generalised to non-uniform spacing), which is a
      ! fixed linear combination of y1..y4. See the comment on the routine:
      ! this linearity is what makes the interpolant flux compatible, and a
      ! monotone limiter would destroy it.
      if (h1 > 0.) then
         s1 = (y2 - y1)/h1
         d2 = (s1*h2 + s2*h1)/(h1 + h2)
      else
         d2 = s2
      end if

      if (h3 > 0.) then
         s3 = (y4 - y3)/h3
         d3 = (s2*h3 + s3*h2)/(h2 + h3)
      else
         d3 = s2
      end if

      t2  = th*th
      t3c = t2*th
      h00 =  2.*t3c - 3.*t2 + 1.
      h10 =     t3c - 2.*t2 + th
      h01 = -2.*t3c + 3.*t2
      h11 =     t3c -    t2

      hermite = h00*y2 + h10*h2*d2 + h01*y3 + h11*h2*d3

   end function hermite


   !> Accumulate sum((tgt-q)^2) and the point count for the zone misfit diagnostic.
   subroutine accum_misfit(zn, ic, q, s)
      use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh

      type(zone_type), intent(in)    :: zn
      integer,         intent(in)    :: ic
      real,            intent(in)    :: q(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh)
      real,            intent(inout) :: s(2)

      integer :: n, i, j, k
      real    :: d

      do n = 1, zn%npts
         if (zn%w(n) <= 0.) cycle
         i = zn%ijk(n,1); j = zn%ijk(n,2); k = zn%ijk(n,3)
         d = sbuf(ic)%cur(zn%src(n)) - q(i,j,k)
         s(1) = s(1) + d*d
         s(2) = s(2) + 1.
      end do

   end subroutine accum_misfit

end module modnesting
