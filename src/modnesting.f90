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
                            nestio_close, nestio_hdr, nestio_tread

   implicit none
   save
   private

   public :: nesting_init, nesting_update_target, nesting_apply, &
             nesting_boundary, nesting_bcpup, nesting_stats,     &
             nesting_restart_write, nesting_restart_read, nesting_finalize
   ! Test hooks: exercised directly by src/tests.f90 (runmodes TEST_NESTING_*).
   public :: nest_shape_fn, nest_union, nest_stagger_coord, nest_flux_residual
   ! Namelist variables: read and broadcast by modstartup.
   public :: lnesting, nestfile, nest_guardwidth, nest_zonewidth, nest_tau,   &
             nest_shape, nest_lateral, nest_top, nest_timeinterp, nest_nwall, &
             nest_lparentgeom, nest_fluxtol, nest_lfluxassert

   logical            :: lnesting         = .false.
   character(len=256) :: nestfile         = ''
   real               :: nest_guardwidth  = 0.
   real               :: nest_zonewidth   = 0.
   real               :: nest_tau         = 0.
   integer            :: nest_shape       = 1
   logical            :: nest_lateral(4)  = .true.
   logical            :: nest_top         = .false.
   integer            :: nest_timeinterp  = 2
   integer            :: nest_nwall       = 1
   logical            :: nest_lparentgeom = .false.
   real               :: nest_fluxtol     = 1.e-10
   logical            :: nest_lfluxassert = .true.

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
   real    :: ttarget  = -1.                !< time the buffer cur(:) was evaluated at
   real    :: phi_last = 0.
   ! Energy injected by the zone forcing since the last nesting_stats call,
   ! split by where it was injected: guard strip (W >= 1) and relaxation ramp
   ! (0 < W < 1). Accumulated in nesting_apply, reported and reset in
   ! nesting_stats (design section 6.4).
   real    :: einj_guard = 0., einj_relax = 0.                 !< last normalised flux residual
   real    :: area_bnd = 0.                 !< total FLUID domain-boundary area
   real    :: twall0   = 0.                 !< wall clock at the end of nesting_init
   real    :: tread0   = 0.                 !< nestio_tread at the end of nesting_init
   integer :: nsolid_zone = 0               !< solid points found inside the zone

   ! restart state, applied by nesting_init when set by nesting_restart_read
   logical :: lrestart_pending = .false.
   integer :: it_lo_restart = 0
   real    :: t_restart = 0.

contains

   !> Called from program.f90 after createmasks/calcfluidvolumes, before readinitfiles.
   !! Opens and validates the parent file, builds the weights and the three zone
   !! point lists, enforces the building-free rule, checks the stored flux
   !! residuals and loads the first parent time levels.
   subroutine nesting_init
      use mpi,       only : MPI_Wtime
      use modglobal, only : cexpnr, timee, dx, dy, xlen, ylen, &
                            BCxm, BCym, BCxm_nesting, BCym_nesting
      use modmpi,    only : myid

      integer :: ierr, f
      real    :: ltot

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

      call nestio_open(trim(nestfile), ierr)
      if (ierr /= 0) call nest_abort('cannot open '//trim(nestfile))
      call nestio_validate()

      nzone = nestio_hdr%nzone
      ntime = nestio_hdr%ntime
      if (nzone < 1) call nest_abort('nzone < 1 in '//trim(nestfile))
      if (ntime < 1) call nest_abort('no time levels in '//trim(nestfile))

      ! A face is forced only if the namelist selects it AND the corresponding
      ! momentum BC is the nesting one.
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
      if (area_bnd <= 0.) call nest_abort('no fluid domain-boundary area found')

      if (.not. nestio_hdr%divergence_corrected .and. myid == 0) then
         write(*,'(a)') ' modnesting: WARNING the input file is not marked divergence_corrected'
      end if

      call check_stored_flux

      ! ---- load the first parent time levels ----
      it_lo = 0
      if (lrestart_pending) then
         call set_interval(t_restart)
         if (it_lo_restart /= it_lo) call nest_abort( &
            'restart time does not fall in the parent interval stored in the restart file')
      else
         call set_interval(timee)
      end if
      call reload_all
      if (lrestart_pending) then
         call eval_target(t_restart)
      else
         call eval_target(timee)
      end if
      lrestart_pending = .false.

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

      call reset_injection

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
         call accum_injection(ww, qnew, qnew - qstar)
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
         call accum_injection(ww, qnew, qnew - qstar)
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
         call accum_injection(ww, qnew, qnew - qstar)
      end do

   end subroutine nesting_apply


   !> Zero the per-substep injection accumulators (nesting_apply overwrites the
   !! tendencies each substep, so the injection is a per-substep quantity).
   subroutine reset_injection
      einj_guard = 0.
      einj_relax = 0.
   end subroutine reset_injection


   !> Accumulate the kinetic energy the zone forcing put in at one point,
   !! q * dq, split by guard strip versus relaxation ramp.
   subroutine accum_injection(ww, qnew, dq)
      real, intent(in) :: ww, qnew, dq

      if (ww >= 1.) then
         einj_guard = einj_guard + qnew*dq
      else
         einj_relax = einj_relax + qnew*dq
      end if

   end subroutine accum_injection


   !> Called from modboundary::boundary. Fills the ghost planes of u0/um, v0/vm
   !! and w0/wm from the parent, in the xmi_driver pattern
   !! (src/modboundary.f90:720).
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
               uu = facval(1, 1, 1, k, j)
               u0(ib, j, k)     = uu
               um(ib, j, k)     = uu
               u0(ib - 1, j, k) = uu
               um(ib - 1, j, k) = uu
               vv = facval(2, 1, 1, k, j)
               v0(ib - 1, j, k) = vv
               vm(ib - 1, j, k) = vv
            end do
            do k = kb, ke + 1
               ww = facval(3, 1, 1, k, j)
               w0(ib - 1, j, k) = ww
               wm(ib - 1, j, k) = ww
            end do
         end do
      end if

      if (lface(2) .and. ierank) then
         do j = jb - 1, je + 1
            do k = kb, ke
               uu = facval(1, 2, nzone + 1, k, j)
               u0(ie + 1, j, k) = uu
               um(ie + 1, j, k) = uu
               vv = facval(2, 2, nzone, k, j)
               v0(ie + 1, j, k) = vv
               vm(ie + 1, j, k) = vv
            end do
            do k = kb, ke + 1
               ww = facval(3, 2, nzone, k, j)
               w0(ie + 1, j, k) = ww
               wm(ie + 1, j, k) = ww
            end do
         end do
      end if

      if (lface(3) .and. jbrank) then
         do i = ib - 1, ie + 1
            do k = kb, ke
               vv = facval(2, 3, 1, k, i)
               v0(i, jb, k)     = vv
               vm(i, jb, k)     = vv
               v0(i, jb - 1, k) = vv
               vm(i, jb - 1, k) = vv
               uu = facval(1, 3, 1, k, i)
               u0(i, jb - 1, k) = uu
               um(i, jb - 1, k) = uu
            end do
            do k = kb, ke + 1
               ww = facval(3, 3, 1, k, i)
               w0(i, jb - 1, k) = ww
               wm(i, jb - 1, k) = ww
            end do
         end do
      end if

      if (lface(4) .and. jerank) then
         do i = ib - 1, ie + 1
            do k = kb, ke
               vv = facval(2, 4, nzone + 1, k, i)
               v0(i, je + 1, k) = vv
               vm(i, je + 1, k) = vv
               uu = facval(1, 4, nzone, k, i)
               u0(i, je + 1, k) = uu
               um(i, je + 1, k) = uu
            end do
            do k = kb, ke + 1
               ww = facval(3, 4, nzone, k, i)
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
   subroutine nesting_bcpup(pup, pvp, pwp, rk3coef)
      use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh, &
                            ibrank, ierank, jbrank, jerank,     &
                            BCxm, BCym, BCxm_nesting, BCym_nesting
      use modfields, only : up, vp
      use modmpi,    only : myid

      real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh), intent(inout) :: pup, pvp, pwp
      real, intent(in) :: rk3coef

      integer :: i, j, k
      real    :: rk3coefi

      if (.not. lnesting) return
      if (.not. linit) return

      rk3coefi = 1./rk3coef

      if (BCxm == BCxm_nesting) then
         if (lface(1) .and. ibrank) then
            do k = kb, ke
               do j = jb - 1, je + 1
                  pup(ib, j, k) = facval(1, 1, 1, k, j)*rk3coefi
                  up(ib, j, k)  = 0. ! u(ib) only evolves according to pressure correction
               end do
            end do
         end if

         if (lface(2) .and. ierank) then
            do k = kb, ke
               do j = jb - 1, je + 1
                  pup(ie + 1, j, k) = facval(1, 2, nzone + 1, k, j)*rk3coefi
                  up(ie + 1, j, k)  = 0.
               end do
            end do
         end if
      end if

      if (BCym == BCym_nesting) then
         if (lface(3) .and. jbrank) then
            do k = kb, ke
               do i = ib - 1, ie + 1
                  pvp(i, jb, k) = facval(2, 3, 1, k, i)*rk3coefi
                  vp(i, jb, k)  = 0.
               end do
            end do
         end if

         if (lface(4) .and. jerank) then
            do k = kb, ke
               do i = ib - 1, ie + 1
                  pvp(i, je + 1, k) = facval(2, 4, nzone + 1, k, i)*rk3coefi
                  vp(i, je + 1, k)  = 0.
               end do
            end do
         end if
      end if

      ! The compatibility condition is a property of the boundary faces only
      ! (design section 3.1) and the solver will not complain on its own (F3).
      phi_last = nest_flux_residual(pup, pvp, pwp, rk3coef)

      if (nest_lfluxassert .and. abs(phi_last) > nest_fluxtol) then
         if (myid == 0) then
            write(*,'(a,es12.5,a,es12.5)') ' modnesting: flux residual ', phi_last, &
               ' exceeds nest_fluxtol = ', nest_fluxtol
         end if
         call nest_abort('normalised boundary flux residual out of tolerance')
      end if

   end subroutine nesting_bcpup


   !> Called from program.f90 alongside statsdump. Reports the flux residual,
   !! the zone misfit and the parent-file read cost (design section 6.4).
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
         write(*,'(a,es12.4)') ' modnesting: zone misfit rms [m/s] = ', rmsmis
         write(*,'(a,es12.4,a,es12.4,a,f8.3)') ' modnesting: |grad p| zone = ', gzone, &
            '  interior = ', gint, '  ratio = ', gratio
         write(*,'(a,es12.4,a,es12.4)') ' modnesting: energy injected guard = ', eg(1), &
            '  relaxation = ', eg(2)
         write(*,'(a,es12.4,a,f6.2,a)') ' modnesting: read time = ', nestio_tread - tread0, &
            ' s (', frac, ' % of run)'
         if (frac > 1.) write(*,'(a)') ' modnesting: WARNING parent I/O exceeds 1 % of runtime'
      end if

   end subroutine nesting_stats


   !> Accumulate sum(|grad p|^2) and the point count, separately over this
   !! component's zone points (acc 1,2) and over its trusted-interior points
   !! (acc 3,4), each evaluated at the component's own staggered location.
   subroutine accum_gradp(zn, ivar, p, acc)
      use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh, dxi, dyi, dzhi

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



   !> Writes the parent interval index and the buffered times so the buffer can
   !! be repositioned on restart (design section 9.5).
   !!
   !! Not currently wired into modsave: nesting_init reconstructs the buffer
   !! state exactly from `timee`, which readinitfiles sets to the restart time,
   !! so a restart needs no stored nesting record. Kept, and exercised by the
   !! unit tests, so that wiring it later is a change of call site only.
   subroutine nesting_restart_write(unit)
      integer, intent(in) :: unit

      if (.not. lnesting) return

      write(unit) it_lo, ttarget, ntime, nzone

   end subroutine nesting_restart_write


   !> Not currently wired into modstartup -- see nesting_restart_write. Reads
   !> the state written by nesting_restart_write.
   !! When nesting_init has already run the buffer is repositioned immediately,
   !! otherwise the state is applied at the end of nesting_init.
   subroutine nesting_restart_read(unit)
      integer, intent(in) :: unit

      integer :: itl, nt, nz
      real    :: tt

      if (.not. lnesting) return

      read(unit) itl, tt, nt, nz

      if (linit) then
         if (nt /= ntime .or. nz /= nzone) call nest_abort( &
            'restart file was written against a different nesting input file')
         it_lo = 0
         call set_interval(tt)
         call reload_all
         call eval_target(tt)
         if (itl /= it_lo) call nest_abort( &
            'restart time does not fall in the parent interval stored in the restart file')
      else
         it_lo_restart    = itl
         t_restart        = tt
         lrestart_pending = .true.
      end if

   end subroutine nesting_restart_read


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
   real function nest_flux_residual(pup, pvp, pwp, rk3coef)
      use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh, dx, dy, dzf, &
                            ibrank, ierank, jbrank, jerank
      use modfields, only : IIu, IIv, IIw, rhobf, rhobh
      use modmpi,    only : comm3d, mpierr, my_real, mpi_sum

      real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh), intent(in) :: pup, pvp, pwp
      real, intent(in) :: rk3coef

      integer :: i, j, k
      real    :: sl(2), sg(2), af

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
            end if
         end do
      end do

      call MPI_ALLREDUCE(sl, sg, 2, MY_REAL, MPI_SUM, comm3d, mpierr)

      if (sg(2) > 0.) then
         nest_flux_residual = sg(1)/sg(2)
      else
         nest_flux_residual = 0.
      end if

   end function nest_flux_residual


   ! ------------------------------------------------------------------ private

   !> Report and stop. Collective: every rank must reach the same condition.
   subroutine nest_abort(message)
      use modmpi, only : myid

      character(len=*), intent(in) :: message

      if (myid == 0) write(*,'(a,a)') ' modnesting: ERROR ', trim(message)
      stop 1

   end subroutine nest_abort


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
      real    :: al, ag, af

      al = 0.

      if (ibrank) then
         do k = kb, ke
            af = dy*dzf(k)
            do j = jb, je
               if (IIu(ib, j, k) == 1) al = al + af
            end do
         end do
      end if
      if (ierank) then
         do k = kb, ke
            af = dy*dzf(k)
            do j = jb, je
               if (IIu(ie + 1, j, k) == 1) al = al + af
            end do
         end do
      end if
      if (jbrank) then
         do k = kb, ke
            af = dx*dzf(k)
            do i = ib, ie
               if (IIv(i, jb, k) == 1) al = al + af
            end do
         end do
      end if
      if (jerank) then
         do k = kb, ke
            af = dx*dzf(k)
            do i = ib, ie
               if (IIv(i, je + 1, k) == 1) al = al + af
            end do
         end do
      end if
      af = dx*dy
      do j = jb, je
         do i = ib, ie
            if (IIw(i, j, kb) == 1)     al = al + af
            if (IIw(i, j, ke + 1) == 1) al = al + af
         end do
      end do

      call MPI_ALLREDUCE(al, ag, 1, MY_REAL, MPI_SUM, comm3d, mpierr)

      fluid_boundary_area = ag

   end function fluid_boundary_area


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

      if (.not. sl_on(c,f)) then
         facval = 0.
         return
      end if

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
      real    :: rl, rg
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

      rl = real(nbad)
      call MPI_ALLREDUCE(rl, rg, 1, MY_REAL, MPI_SUM, comm3d, mpierr)
      if (rg > 0.) call nest_abort('zone points fall outside the slabs stored in the '// &
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
   !! Phi is a linear functional of the boundary data (design section 3.1
   !! item 2), so checking every stored level is sufficient for every
   !! time-interpolated target. Only the four boundary-normal slabs are read,
   !! and only on the ranks that own a domain-boundary face, so the cost is a
   !! perimeter read rather than a full read of the file.
   !! Note that net_volume_flux in the file is the PRE-correction flux
   !! (spec section 5) and carries no information about the corrected data.
   subroutine check_stored_flux
      use modglobal, only : ib, ie, jb, je, kb, ke, dx, dy, dzf, &
                            ibrank, ierank, jbrank, jerank
      use modfields, only : IIu, IIv, rhobf
      use modmpi,    only : myid, comm3d, mpierr, my_real, mpi_sum

      integer :: it, i, j, k, ierr
      real    :: sl, sg, phi, af
      logical :: lw, le, ls, ln
      real, allocatable :: bw(:,:,:), be(:,:,:), bs(:,:,:), bn(:,:,:)

      lw = lface(1) .and. ibrank .and. sl_on(1,1)
      le = lface(2) .and. ierank .and. sl_on(1,2)
      ls = lface(3) .and. jbrank .and. sl_on(2,3)
      ln = lface(4) .and. jerank .and. sl_on(2,4)

      if (lw) allocate(bw(sn1(1,1), sn2(1,1), sn3(1,1)))
      if (le) allocate(be(sn1(1,2), sn2(1,2), sn3(1,2)))
      if (ls) allocate(bs(sn1(2,3), sn2(2,3), sn3(2,3)))
      if (ln) allocate(bn(sn1(2,4), sn2(2,4), sn3(2,4)))

      do it = 1, ntime
         sl = 0.

         if (lw) then
            call nestio_read('u_west', it, sds(1,1), sn3(1,1), bw, ierr)
            if (ierr /= 0) call nest_abort('failed to read u_west from '//trim(nestfile))
            do k = kb, ke
               af = dy*dzf(k)
               do j = jb, je
                  if (IIu(ib,j,k) == 1) &
                     sl = sl - rhobf(k)*bw(1, k, dloc(1,1,j))*af
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
                     sl = sl + rhobf(k)*be(nzone + 1, k, dloc(1,2,j))*af
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
                     sl = sl - rhobf(k)*bs(1, k, dloc(2,3,i))*af
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
                     sl = sl + rhobf(k)*bn(nzone + 1, k, dloc(2,4,i))*af
               end do
            end do
         end if

         call MPI_ALLREDUCE(sl, sg, 1, MY_REAL, MPI_SUM, comm3d, mpierr)
         phi = sg/area_bnd

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
         ' stored time levels are flux balanced'

   end subroutine check_stored_flux


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
               call nest_abort('failed to read a slab from '//trim(nestfile))
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
   !! nest_timeinterp: 1 linear, 2 monotone cubic Hermite (Fritsch-Carlson).
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

      ttarget = t

   end subroutine eval_target


   !> Monotone cubic Hermite on [t2,t3], with Fritsch-Carlson limited slopes
   !! from the neighbouring intervals. Exact for data linear in time.
   real function hermite(y1, y2, y3, y4, h1, h2, h3, th)
      real, intent(in) :: y1, y2, y3, y4, h1, h2, h3, th

      real :: s1, s2, s3, d2, d3, wa, wb, t2, t3c, h00, h10, h01, h11

      s2 = (y3 - y2)/h2

      if (h1 > 0.) then
         s1 = (y2 - y1)/h1
         if (s1*s2 <= 0.) then
            d2 = 0.
         else
            wa = 2.*h2 + h1
            wb = h2 + 2.*h1
            d2 = (wa + wb)/(wa/s1 + wb/s2)
         end if
      else
         d2 = s2
      end if

      if (h3 > 0.) then
         s3 = (y4 - y3)/h3
         if (s2*s3 <= 0.) then
            d3 = 0.
         else
            wa = 2.*h3 + h2
            wb = h3 + 2.*h2
            d3 = (wa + wb)/(wa/s2 + wb/s3)
         end if
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
