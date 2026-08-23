!> \file modboundary.f90
!! All boundary conditions are in this file, except for immersed boundaries.
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!
!! \par Revision list
!! \par Authors
!!
module modboundary
   use mpi

   implicit none
   save
   private
   public :: initboundary, boundary_conditions, boundary, grwdamp, ksp, tqaver, exchange_halos, halos, bcp, bcpup, closurebc
#if defined(_GPU)
   public :: halos_device, boundary_device
#if defined(UDALES_DEBUG)
   public :: xm_periodic_device, xT_periodic_device, xq_periodic_device, xs_periodic_device, &
             ym_periodic_device, yT_periodic_device, yq_periodic_device, ys_periodic_device
#endif
#endif
   integer :: ksp = -1 !<    lowest level of sponge layer
   real :: rnu0 = 2.75e-3
contains
   !>
   !! Initializing Boundary; specifically the sponge layer
   !>
   subroutine initboundary
      use modglobal,    only : ib, kb, ke, kh, kmax, pi, zf, iplane
      use modinletdata, only : irecy
      use modfields,    only : tsc
      implicit none

      real    :: zspb, zspt
      integer :: k
      allocate (tsc(kb:ke + kh))
      ! Sponge layer
      if (ksp == -1) then
         !      ksp  = min(3*kmax/4,kmax - 15)
         ksp = (kb - 1) + max(min(3*kmax/4, kmax - 15),1)
      end if

      zspb = zf(ksp)
      zspt = zf(ke)

      tsc(kb:ksp - 1) = 0.0
      do k = ksp, ke
         tsc(k) = rnu0*sin(0.5*pi*(zf(k) - zspb)/(zspt - zspb))**2
      end do
      tsc(ke + 1) = tsc(ke)
      irecy = ib + iplane

   end subroutine initboundary

   !>
   !! Fill halo cells, including ghost cells outside domain
   ! Needs to be called before divergence is calculated
   subroutine exchange_halos
    implicit none
#if defined(_GPU)
    call halos_device
#else
    call halos
#endif
   end subroutine exchange_halos
   subroutine halos
      use modglobal, only : ihc, jhc, khc, nsv, &
                            BCxm, BCym, BCxT, BCyT, BCxq, BCyq, BCxs, BCys, &
                            BCxm_periodic, BCxT_periodic, BCxq_periodic, BCxs_periodic, &
                            BCym_periodic, BCyT_periodic, BCyq_periodic, BCys_periodic, &
                            ibrank, ierank, jbrank, jerank
      use modfields, only : u0, v0, w0, um, vm, wm, thl0, thlm, qt0, qtm, sv0, svm, thl0c
      use m_halo, only : halo_exchange
      implicit none
      integer n

#if defined(_GPU)
      ! The GPU build of 2DECOMP accepts device arrays. During startup, before
      ! initCUDA allocates the persistent CUDA arrays, use temporary OpenACC
      ! mappings for the host fields and copy the exchanged halos back.
      !$acc data create(u0, v0, w0, um, vm, wm, thl0, thlm, thl0c, qt0, qtm, sv0, svm)
      !$acc update device(u0, v0, w0, um, vm, wm, thl0, thlm, thl0c, qt0, qtm, sv0, svm)
      !$acc host_data use_device(u0, v0, w0, um, vm, wm, thl0, thlm, thl0c, qt0, qtm, sv0, svm)
#endif
      call halo_exchange(u0, 3)
      call halo_exchange(v0, 3)
      call halo_exchange(w0, 3)
      call halo_exchange(um, 3)
      call halo_exchange(vm, 3)
      call halo_exchange(wm, 3)
      call halo_exchange(thl0, 3)
      call halo_exchange(thlm, 3)
      call halo_exchange(thl0c, 3, opt_levels=(/ihc,jhc,khc/))
      call halo_exchange(qt0, 3)
      call halo_exchange(qtm, 3)
      do n = 1, nsv
         call halo_exchange(sv0(:, :, :, n), 3, opt_levels=(/ihc,jhc,khc/))
         call halo_exchange(svm(:, :, :, n), 3, opt_levels=(/ihc,jhc,khc/))
      enddo
#if defined(_GPU)
      !$acc end host_data
      !$acc update host(u0, v0, w0, um, vm, wm, thl0, thlm, thl0c, qt0, qtm, sv0, svm)
      !$acc end data
#endif

      if (ibrank .and. ierank) then ! not parallelized in x
        if (BCxm == BCxm_periodic) call xm_periodic
        if (BCxT == BCxT_periodic) call xT_periodic
        if (BCxq == BCxq_periodic) call xq_periodic
        if (BCxs == BCxs_periodic) call xs_periodic
      end if

      if (jbrank .and. jerank) then ! not parallelized in y
        if (BCym == BCym_periodic) call ym_periodic
        if (BCyT == BCyT_periodic) call yT_periodic
        if (BCyq == BCyq_periodic) call yq_periodic
        if (BCys == BCys_periodic) call ys_periodic
      end if
   end subroutine halos
#if defined(_GPU)
   !> Fill device-resident halo and periodic ghost cells
   subroutine halos_device
      use modglobal, only : ihc, jhc, khc, nsv, &
                            BCxm, BCym, BCxT, BCyT, BCxq, BCyq, BCxs, BCys, &
                            BCxm_periodic, BCxT_periodic, BCxq_periodic, BCxs_periodic, &
                            BCym_periodic, BCyT_periodic, BCyq_periodic, BCys_periodic, &
                            ibrank, ierank, jbrank, jerank, &
                            ltempeq, lmoist, iadv_thl, iadv_kappa
      use m_halo,    only : halo_exchange
      use modcuda,   only : u0_d, v0_d, w0_d, um_d, vm_d, wm_d, &
                            thl0_d, thlm_d, qt0_d, qtm_d, sv0_d, svm_d, thl0c_d
      implicit none
      integer n

      call halo_exchange(u0_d, 3)
      call halo_exchange(v0_d, 3)
      call halo_exchange(w0_d, 3)
      call halo_exchange(um_d, 3)
      call halo_exchange(vm_d, 3)
      call halo_exchange(wm_d, 3)
      if (ltempeq) then
        call halo_exchange(thl0_d, 3)
        call halo_exchange(thlm_d, 3)
        if (iadv_thl == iadv_kappa) then
          call halo_exchange(thl0c_d, 3, opt_levels=(/ihc,jhc,khc/))
        end if
      end if
      if (lmoist) then
        call halo_exchange(qt0_d, 3)
        call halo_exchange(qtm_d, 3)
      end if
      do n = 1, nsv
         call halo_exchange(sv0_d(:, :, :, n), 3, opt_levels=(/ihc,jhc,khc/))
         call halo_exchange(svm_d(:, :, :, n), 3, opt_levels=(/ihc,jhc,khc/))
      enddo

      if (ibrank .and. ierank) then ! not parallelized in x
        if (BCxm == BCxm_periodic) call xm_periodic_device
        if (BCxT == BCxT_periodic) call xT_periodic_device
        if (BCxq == BCxq_periodic) call xq_periodic_device
        if (BCxs == BCxs_periodic) call xs_periodic_device
      end if

      if (jbrank .and. jerank) then ! not parallelized in y
        if (BCym == BCym_periodic) call ym_periodic_device
        if (BCyT == BCyT_periodic) call yT_periodic_device
        if (BCyq == BCyq_periodic) call yq_periodic_device
        if (BCys == BCys_periodic) call ys_periodic_device
      end if
   end subroutine halos_device
#endif


   !>
   !! Set boundary conditions for the next timestep
   ! Will result in velocity field being not divergence-free

   !> Apply the boundary conditions, wherever the fields live.
   !!
   !! Only for the call inside the time loop. The one in the initialisation
   !! sequence has to name boundary directly, because it runs before initCUDA
   !! and there is nothing on the device for boundary_device to write.
   subroutine boundary_conditions
    implicit none
#if defined(_GPU)
    call boundary_device
#else
    call boundary
#endif
   end subroutine boundary_conditions

   subroutine boundary
      use modglobal,      only : kb, ke, khc, dzf, zh, nsv, &
                                 ltempeq, lmoist, luvolflowr, luoutflowr, &
                                 BCxm, BCym, BCxT, BCyT, BCxq, BCyq, BCxs, BCys, BCtopm, BCtopT, BCtopq, BCtops, &
                                 BCtopm_freeslip, BCtopm_noslip, BCtopm_pressure, &
                                 BCtopT_flux, BCtopT_value, BCtopq_flux, BCtopq_value, BCtops_flux, BCtops_value, &
                                 BCxm_periodic, BCxm_profile, BCxm_driver, &
                                 BCxT_periodic, BCxT_profile, BCxT_driver, &
                                 BCxq_periodic, BCxq_profile, BCxq_driver, &
                                 BCxs_periodic, BCxs_profile, BCxs_driver, BCxs_custom, &
                                 BCym_periodic, BCym_profile, BCyT_periodic, BCyT_profile, &
                                 BCyq_periodic, BCyq_profile, BCys_periodic, &
                                 ibrank, ierank, jbrank, jerank, e12min, idriver, &
                                 Uinf, Vinf, &
                                 rk3step, lchunkread
      use modfields,      only : u0, v0, w0, um, vm, wm, thl0, thlm, qt0, qtm, e120, e12m, u0av, v0av, uouttot, vouttot, thl0c
      use modsubgriddata, only : ekh, ekm, loneeqn
      use modsurfdata,    only : thl_top, qt_top, sv_top, wttop, wqtop, wsvtop
      use modmpi,         only : slabsum, avey_ibm
      use moddriver,      only : drivergen, driverchunkread
      use modinletdata,   only : ubulk, vbulk !, iangle
      implicit none
      real, dimension(kb:ke) :: uaverage, vaverage
      integer k, n

     ! if not using massflowrate need to set outflow velocity
     if (luoutflowr) then
        ! do nothing - calculated in modforces
     elseif (.not. luvolflowr) then
        !ubulk = sum(u0av)/(ke-kb+1)
        do k = kb, ke
           uaverage(k) = u0av(k)*dzf(k)
        end do

        do k = kb, ke
           vaverage(k) = v0av(k)*dzf(k)
        end do
        ! need a method to know if we have all blocks at lowest cell kb
        ! assuming this for now (hence kb+1)
        uouttot = sum(uaverage(kb:ke))/(zh(ke + 1) - zh(kb+1))
        vouttot = sum(vaverage(kb:ke))/(zh(ke + 1) - zh(kb+1))
     else
        uouttot = ubulk
        vouttot = vbulk
     end if

     ! Bottom BC - many ways of enforcing this but this is simplest
     ! Other variables handled by bottom
     wm(:, :, kb) = 0.
     w0(:, :, kb) = 0.

     !! Top
     ! Momentum
     select case(BCtopm)
     case(BCtopm_freeslip)
        !free-slip = zero-flux
        call fluxtop(um, ekm, 0.0)
        call fluxtop(u0, ekm, 0.0)
        call fluxtop(vm, ekm, 0.0)
        call fluxtop(v0, ekm, 0.0)
        w0(:, :, ke + 1) = 0.0
        wm(:, :, ke + 1) = 0.0
        if (loneeqn) then
          e120(:, :, ke + 1) = e12min
          e12m(:, :, ke + 1) = e12min
        end if
     case(BCtopm_noslip)
        !no-slip = fixed velocity at wall
        call valuetop(um, Uinf)
        call valuetop(u0, Uinf)
        call valuetop(vm, Vinf)
        call valuetop(v0, Vinf)
        w0(:, :, ke + 1) = 0.0
        wm(:, :, ke + 1) = 0.0
      case(BCtopm_pressure)
         call fluxtop(um, ekm, 0.0)
         call fluxtop(u0, ekm, 0.0)
         call fluxtop(vm, ekm, 0.0)
         call fluxtop(v0, ekm, 0.0)
         if (loneeqn) then
           e120(:, :, ke + 1) = e12min
           e12m(:, :, ke + 1) = e12min
         end if
         ! w considered in modpois
      case default
        write(0, *) "ERROR: top boundary type for velocity undefined"
        stop 1
     end select

     ! Temperature
     select case(BCtopT)
     case(BCtopT_flux)
        call fluxtop(thlm, ekh, wttop)
        call fluxtop(thl0, ekh, wttop)
        do n=1,khc
           thl0c(:,:,ke+n) = thl0c(:,:,ke+n-1)
        end do
     case(BCtopT_value)
        call valuetop(thlm, thl_top)
        call valuetop(thl0, thl_top)
     case default
        write(0, *) "ERROR: top boundary type for temperature undefined"
        stop 1

     end select

     ! Moisture
     select case(BCtopq)
     case(BCtopq_flux)
        call fluxtop(qtm, ekh, wqtop)
        call fluxtop(qt0, ekh, wqtop)
     case(BCtopq_value)
        call valuetop(qtm, qt_top)
        call valuetop(qt0, qt_top)
     case default
        write(0, *) "ERROR: top boundary type for moisture undefined"
        stop 1
     end select

     ! Scalars
     select case(BCtops)
     case(BCtops_flux)
        call fluxtopscal(wsvtop)
        call fluxtopscal(wsvtop)
     case(BCtops_value)
        call valuetopscal(sv_top)
        call valuetopscal(sv_top)
     case default
        write(0, *) "ERROR: top boundary type for scalars undefined"
        stop 1
     end select

     if (idriver == 1) call drivergen ! Should be moved elsewhere, as not related to boundary conditions.

     ! x inlet
     if (ibrank) then ! set inlet
       ! Momentum
       select case(BCxm)
       case(BCxm_periodic)
         ! Handled in halos
       case(BCxm_profile)
         !uouttot = cos(iangle)*ubulk
         call xmi_profile
       case(BCxm_driver)
         !uouttot = ubulk ! does this hold for all forcings of precursor simulations? tg3315
         if(rk3step==0 .or. rk3step==3) then
          if (lchunkread) call driverchunkread
          call drivergen ! think this should be done at the start of an rk3 loop?
         end if
         call xmi_driver
       case default
         write(0, *) "ERROR: lateral boundary type for veloctiy in x-direction undefined"
         stop 1
       end select

       ! Temperature
       if (ltempeq) then
         select case(BCxT)
         case(BCxT_periodic) ! periodic
           ! Handled in halos
         case(BCxT_profile) ! profile
           call xTi_profile
         case(BCxT_driver)
           call xTi_driver
         case default
           write(0, *) "ERROR: lateral boundary type for temperature in x-direction undefined"
           stop 1
         end select
       end if

       ! Moisture
       if (lmoist) then
         select case(BCxq)
         case(BCxq_periodic)
           ! Handled in halos
         case(BCxq_profile)
           call xqi_profile
         case(BCxq_driver)
           call xqi_driver
         case default
           write(0, *) "ERROR: lateral boundary type for humidity in x-direction undefined"
           stop 1
         end select
       end if

       ! Scalars
       if (nsv > 0) then
         select case(BCxs)
         case(BCxs_periodic)
           ! Handled in halos
         case(BCxs_profile)
           call xsi_profile
         case(BCxs_driver)
           call xsi_driver
        case(BCxs_custom)
           call xsi_custom
         case default
           write(0, *) "ERROR: lateral boundary type for scalars in x-direction undefined"
           stop 1
         end select
       end if

     end if !ibrank

     if (jbrank) then ! set y inlet
       ! Momentum
       select case(BCym)
       case(BCym_periodic)
         ! Handled in halos
       case(BCym_profile)
         call ymi_profile
       case default
         write(0, *) "ERROR: lateral boundary type for veloctiy in y-direction undefined"
         stop 1
       end select

       ! Temperature
       if (ltempeq) then
         select case(BCyT)
         case(BCyT_periodic)
           ! Handled in halos
         case(BCyT_profile)
           call yTi_profile
         case default
           write(0, *) "ERROR: lateral boundary type for temperature in y-direction undefined"
           stop 1
         end select
       end if

       ! Moisture
       if (lmoist) then
         select case(BCyq)
         case(BCyq_periodic)
           ! Handled in halos
         case(BCyq_profile)
           call yqi_profile
         case default
           write(0, *) "ERROR: lateral boundary type for humidity in y-direction undefined"
           stop 1
         end select
       end if

       if (nsv > 0) then !scalars
         select case(BCys)
         case(1)
           ! Handled in halos
         case(2)
           call ysi_profile
         case default
           write(0, *) "ERROR: lateral boundary type for scalars in y-direction undefined"
           stop 1
         end select
       end if

     end if !jbrank

     !> Outlet
     ! Currently only outflow boundary conditions are convective
     if (ierank) then
       if (BCxm .ne. BCxm_periodic) call xmo_convective
       if ((BCxT .ne. BCxT_periodic) .and. ltempeq) call xTo_convective
       if ((BCxq .ne. BCxq_periodic) .and. lmoist ) call xqo_convective
       if ((BCxs .ne. BCxs_periodic) .and. nsv > 0) call xso_convective
     end if

     if (jerank) then
       if (BCym .ne. BCym_periodic) call ymo_convective
       if ((BCyT .ne. BCyT_periodic) .and. ltempeq) call yTo_convective
       if ((BCyq .ne. BCyq_periodic) .and. lmoist ) call yqo_convective
       if ((BCys .ne. BCys_periodic) .and. nsv > 0) call yso_convective
     end if

   end subroutine boundary


#if defined(_GPU)
   !> Boundary conditions, applied to the device copies.
   !!
   !! Mirrors boundary above branch for branch, and every leaf it calls sits
   !! beside the host leaf it was written from. What is not mirrored is where
   !! it sits in the step. The host version had to run after the post-Poisson
   !! handover, because it wrote host memory that the next updateDevice
   !! carried back up; this one runs before the handover, so the boundary
   !! planes are on the device first and the host copies are fetched with them
   !! already in place.
   !!
   !! That reordering is what makes the invalidateHostFields call at the end
   !! load-bearing: fielddump and statsdump have pulled their fields down
   !! before this point, and those copies stop matching the device the moment
   !! the kernels below run.
   !!
   !! uouttot and vouttot stay on the host. They reduce u0av and v0av, which
   !! thermodynamics produces on the host anyway, and they are two scalars
   !! over ktot levels - there is nothing here for a kernel to do.
   subroutine boundary_device
      use modglobal,      only : kb, ke, dzf, zh, nsv, &
                                 ltempeq, lmoist, luvolflowr, luoutflowr, &
                                 BCxm, BCym, BCxT, BCyT, BCxq, BCyq, BCxs, BCys, BCtopm, BCtopT, BCtopq, BCtops, &
                                 BCtopm_freeslip, BCtopm_noslip, BCtopm_pressure, &
                                 BCtopT_flux, BCtopT_value, BCtopq_flux, BCtopq_value, BCtops_flux, BCtops_value, &
                                 BCxm_periodic, BCxm_profile, BCxm_driver, &
                                 BCxT_periodic, BCxT_profile, BCxT_driver, &
                                 BCxq_periodic, BCxq_profile, BCxq_driver, &
                                 BCxs_periodic, BCxs_profile, BCxs_driver, BCxs_custom, &
                                 BCym_periodic, BCym_profile, BCyT_periodic, BCyT_profile, &
                                 BCyq_periodic, BCyq_profile, BCys_periodic, &
                                 ibrank, ierank, jbrank, jerank, idriver, &
                                 rk3step, lchunkread
      use modfields,      only : u0av, v0av, uouttot, vouttot
      use modsubgriddata, only : loneeqn
      use modcuda,        only : invalidateHostFields, updateDriverPlanesDevice
      use moddriver,      only : drivergen, driverchunkread
      use modinletdata,   only : ubulk, vbulk
      implicit none
      real, dimension(kb:ke) :: uaverage, vaverage
      integer k

     ! if not using massflowrate need to set outflow velocity
     if (luoutflowr) then
        ! do nothing - calculated in modforces
     elseif (.not. luvolflowr) then
        do k = kb, ke
           uaverage(k) = u0av(k)*dzf(k)
        end do

        do k = kb, ke
           vaverage(k) = v0av(k)*dzf(k)
        end do
        uouttot = sum(uaverage(kb:ke))/(zh(ke + 1) - zh(kb+1))
        vouttot = sum(vaverage(kb:ke))/(zh(ke + 1) - zh(kb+1))
     else
        uouttot = ubulk
        vouttot = vbulk
     end if

     ! Bottom BC
     call bottom_w_device

     !! Top
     ! Momentum
     select case(BCtopm)
     case(BCtopm_freeslip)
        !free-slip = zero-flux
        call fluxtop_uv_device
        call topw_zero_device
        if (loneeqn) call tope12_device
     case(BCtopm_noslip)
        !no-slip = fixed velocity at wall
        call valuetop_uv_device
        call topw_zero_device
     case(BCtopm_pressure)
        call fluxtop_uv_device
        if (loneeqn) call tope12_device
        ! w considered in modpois
     case default
        write(0, *) "ERROR: top boundary type for velocity undefined"
        stop 1
     end select

     ! Temperature
     select case(BCtopT)
     case(BCtopT_flux)
        call fluxtop_thl_device
        call topthl0c_device
     case(BCtopT_value)
        call valuetop_thl_device
     case default
        write(0, *) "ERROR: top boundary type for temperature undefined"
        stop 1
     end select

     ! Moisture
     select case(BCtopq)
     case(BCtopq_flux)
        call fluxtop_qt_device
     case(BCtopq_value)
        call valuetop_qt_device
     case default
        write(0, *) "ERROR: top boundary type for moisture undefined"
        stop 1
     end select

     ! Scalars. The host branch calls each of these twice; both are writes to
     ! the ghost levels from level ke, which they do not touch, so the second
     ! call cannot change anything the first did not already do.
     select case(BCtops)
     case(BCtops_flux)
        call fluxtop_sv_device
     case(BCtops_value)
        call valuetop_sv_device
     case default
        write(0, *) "ERROR: top boundary type for scalars undefined"
        stop 1
     end select

     ! Recording a driver plane reads the prognostic fields on the host, and
     ! this is the point in the sequence the host path read them at. See
     ! updateHostForDriverDump, which writedriverfile calls for that.
     if (idriver == 1) call drivergen ! Should be moved elsewhere, as not related to boundary conditions.

     ! x inlet
     if (ibrank) then ! set inlet
       ! Momentum
       select case(BCxm)
       case(BCxm_periodic)
         ! Handled in halos
       case(BCxm_profile)
         call xmi_profile_device
       case(BCxm_driver)
         if(rk3step==0 .or. rk3step==3) then
          if (lchunkread) call driverchunkread
          call drivergen ! think this should be done at the start of an rk3 loop?
          ! drivergen has just rewritten every driver plane on the host, and
          ! the four driver leaves below read them on the device.
          call updateDriverPlanesDevice
         end if
         call xmi_driver_device
       case default
         write(0, *) "ERROR: lateral boundary type for veloctiy in x-direction undefined"
         stop 1
       end select

       ! Temperature
       if (ltempeq) then
         select case(BCxT)
         case(BCxT_periodic) ! periodic
           ! Handled in halos
         case(BCxT_profile) ! profile
           call xTi_profile_device
         case(BCxT_driver)
           call xTi_driver_device
         case default
           write(0, *) "ERROR: lateral boundary type for temperature in x-direction undefined"
           stop 1
         end select
       end if

       ! Moisture
       if (lmoist) then
         select case(BCxq)
         case(BCxq_periodic)
           ! Handled in halos
         case(BCxq_profile)
           call xqi_profile_device
         case(BCxq_driver)
           call xqi_driver_device
         case default
           write(0, *) "ERROR: lateral boundary type for humidity in x-direction undefined"
           stop 1
         end select
       end if

       ! Scalars
       if (nsv > 0) then
         select case(BCxs)
         case(BCxs_periodic)
           ! Handled in halos
         case(BCxs_profile)
           call xsi_profile_device
         case(BCxs_driver)
           call xsi_driver_device
        case(BCxs_custom)
           call xsi_custom_device
         case default
           write(0, *) "ERROR: lateral boundary type for scalars in x-direction undefined"
           stop 1
         end select
       end if

     end if !ibrank

     if (jbrank) then ! set y inlet
       ! Momentum
       select case(BCym)
       case(BCym_periodic)
         ! Handled in halos
       case(BCym_profile)
         call ymi_profile_device
       case default
         write(0, *) "ERROR: lateral boundary type for veloctiy in y-direction undefined"
         stop 1
       end select

       ! Temperature
       if (ltempeq) then
         select case(BCyT)
         case(BCyT_periodic)
           ! Handled in halos
         case(BCyT_profile)
           call yTi_profile_device
         case default
           write(0, *) "ERROR: lateral boundary type for temperature in y-direction undefined"
           stop 1
         end select
       end if

       ! Moisture
       if (lmoist) then
         select case(BCyq)
         case(BCyq_periodic)
           ! Handled in halos
         case(BCyq_profile)
           call yqi_profile_device
         case default
           write(0, *) "ERROR: lateral boundary type for humidity in y-direction undefined"
           stop 1
         end select
       end if

       if (nsv > 0) then !scalars
         select case(BCys)
         case(1)
           ! Handled in halos
         case(2)
           call ysi_profile_device
         case default
           write(0, *) "ERROR: lateral boundary type for scalars in y-direction undefined"
           stop 1
         end select
       end if

     end if !jbrank

     !> Outlet
     ! Currently only outflow boundary conditions are convective
     if (ierank) then
       if (BCxm .ne. BCxm_periodic) call xmo_convective_device
       if ((BCxT .ne. BCxT_periodic) .and. ltempeq) call xTo_convective_device
       if ((BCxq .ne. BCxq_periodic) .and. lmoist ) call xqo_convective_device
       if ((BCxs .ne. BCxs_periodic) .and. nsv > 0) call xso_convective_device
     end if

     if (jerank) then
       if (BCym .ne. BCym_periodic) call ymo_convective_device
       if ((BCyT .ne. BCyT_periodic) .and. ltempeq) call yTo_convective_device
       if ((BCyq .ne. BCyq_periodic) .and. lmoist ) call yqo_convective_device
       if ((BCys .ne. BCys_periodic) .and. nsv > 0) call yso_convective_device
     end if

     ! Everything the host holds of these fields is now one boundary condition
     ! behind the device. Nothing downstream may skip a transfer on the grounds
     ! that it already has the field.
     call invalidateHostFields

   end subroutine boundary_device
#endif


   !> Re-apply the flux-type top boundary conditions after the diffusivity moved.
   !!
   !! boundary sets the top ghost level from the eddy diffusivity of the step
   !! before; closurebc then recomputes ekm and ekh, which changes what the
   !! same flux implies for that level. So the flux branches - and only those,
   !! the value branches do not involve a diffusivity - are applied again here.
   !!
   !! The device path calls the same routines boundary_device calls, rather
   !! than a second copy of them. It used to be a second copy, with the
   !! division of fluxtop's expression arranged differently on the two sides.
   subroutine reassure_fluxtop_boundary
    use modglobal,      only : ltempeq, lmoist, nsv, &
                               BCtopm, BCtopT, BCtopq, BCtops, &
                               BCtopm_freeslip, BCtopm_pressure, &
                               BCtopT_flux, BCtopq_flux, BCtops_flux
#if !defined(_GPU)
    use modfields,      only : u0, v0, um, vm, thl0, thlm, qt0, qtm
    use modsubgriddata, only : ekh, ekm
    use modsurfdata,    only : wttop, wqtop, wsvtop
#endif
    implicit none

#if defined(_GPU)
    select case(BCtopm)
      case(BCtopm_freeslip, BCtopm_pressure)
        call fluxtop_uv_device
      case default
    end select

    if (ltempeq .and. (BCtopT .eq. BCtopT_flux)) call fluxtop_thl_device
    if (lmoist  .and. (BCtopq .eq. BCtopq_flux)) call fluxtop_qt_device
    if (nsv > 0 .and. (BCtops .eq. BCtops_flux)) call fluxtop_sv_device
#else
    select case(BCtopm)
      case(BCtopm_freeslip)
        !free-slip = zero-flux
        call fluxtop(um, ekm, 0.0)
        call fluxtop(u0, ekm, 0.0)
        call fluxtop(vm, ekm, 0.0)
        call fluxtop(v0, ekm, 0.0)
      case(BCtopm_pressure)
        call fluxtop(um, ekm, 0.0)
        call fluxtop(u0, ekm, 0.0)
        call fluxtop(vm, ekm, 0.0)
        call fluxtop(v0, ekm, 0.0)
      case default
    end select

    if (ltempeq .and. (BCtopT .eq. BCtopT_flux)) then
      call fluxtop(thlm, ekh, wttop)
      call fluxtop(thl0, ekh, wttop)
    end if

    if (lmoist .and. (BCtopq .eq. BCtopq_flux)) then
      call fluxtop(qtm, ekh, wqtop)
      call fluxtop(qt0, ekh, wqtop)
    end if

    if (nsv > 0 .and. (BCtops .eq. BCtops_flux)) then
      call fluxtopscal(wsvtop)
    end if
#endif
   end subroutine reassure_fluxtop_boundary


   subroutine closurebc
     use modglobal,      only : ib, ie, jb, je, kb, ke, numol, prandtlmoli, &
                                ibrank, ierank, jbrank, jerank, BCtopm, BCxm, BCym, &
                                BCtopm_freeslip, BCtopm_noslip, BCtopm_pressure, &
                                BCxm_periodic, BCym_periodic
     use m_halo,         only : halo_exchange
#if defined(_GPU)
     use modcuda,        only : ekm_d, ekh_d
#else
     use modsubgriddata, only : ekm, ekh
#endif
     integer :: i, j

#if defined(_GPU)
     call halo_exchange(ekm_d, 3)
     call halo_exchange(ekh_d, 3)
#else
     call halo_exchange(ekm, 3)
     call halo_exchange(ekh, 3)
#endif

     ! Top and bottom
     if ((BCtopm .eq. BCtopm_freeslip) .or. (BCtopm .eq. BCtopm_pressure)) then
       !$acc kernels default(present)
       do j = jb - 1, je + 1
         do i = ib - 1, ie + 1
#if defined(_GPU)
           ekm_d(i, j, ke + 1) = ekm_d(i, j, ke)
           ekh_d(i, j, ke + 1) = ekh_d(i, j, ke)
           ekm_d(i, j, kb - 1) = 2.*numol - ekm_d(i, j, kb)
           ekh_d(i, j, kb - 1) = (2.*numol*prandtlmoli) - ekh_d(i, j, kb)
#else
           ekm(i, j, ke + 1) = ekm(i, j, ke) ! zero-gradient top wall
           ekh(i, j, ke + 1) = ekh(i, j, ke) ! zero-gradient top wall
           ekm(i, j, kb - 1) = 2.*numol - ekm(i, j, kb) ! no-slip lower wall
           ekh(i, j, kb - 1) = (2.*numol*prandtlmoli) - ekh(i, j, kb) ! no-slip lower wall
#endif
         end do
       end do
       !$acc end kernels
     else if (BCtopm .eq. BCtopm_noslip) then
       !$acc kernels default(present)
       do j = jb - 1, je + 1
         do i = ib - 1, ie + 1
#if defined(_GPU)
           ekm_d(i, j, ke + 1) = 2.*numol - ekm_d(i, j, ke)
           ekh_d(i, j, ke + 1) = (2.*numol*prandtlmoli) - ekh_d(i, j, ke)
           ekm_d(i, j, kb - 1) = 2.*numol - ekm_d(i, j, kb)
           ekh_d(i, j, kb - 1) = (2.*numol*prandtlmoli) - ekh_d(i, j, kb)
#else
           ekm(i, j, ke + 1) = 2.*numol - ekm(i, j, ke) ! no-slip top wall
           ekh(i, j, ke + 1) = (2.*numol*prandtlmoli) - ekh(i, j, ke) ! no-slip top wall
           ekm(i, j, kb - 1) = 2.*numol - ekm(i, j, kb) ! no-slip lower wall
           ekh(i, j, kb - 1) = (2.*numol*prandtlmoli) - ekh(i, j, kb) ! no-slip lower wall
#endif
         end do
       end do
       !$acc end kernels
     end if

     if (BCxm .ne. BCxm_periodic) then ! inflow/outflow
       if (ibrank) then
#if defined (_GPU)
         !$acc kernels default(present)
         ekm_d(ib - 1, :, :) = ekm_d(ib, :, :)
         ekh_d(ib - 1, :, :) = ekh_d(ib, :, :)
         !$acc end kernels
#else
         ekm(ib - 1, :, :) = ekm(ib, :, :)
         ekh(ib - 1, :, :) = ekh(ib, :, :)
#endif
       end if
       if (ierank) then
#if defined (_GPU)
         !$acc kernels default(present)
         ekm_d(ie + 1, :, :) = ekm_d(ie, :, :)
         ekh_d(ie + 1, :, :) = ekh_d(ie, :, :)
         !$acc end kernels
#else
         ekm(ie + 1, :, :) = ekm(ie, :, :)
         ekh(ie + 1, :, :) = ekh(ie, :, :)
#endif
       end if
     else ! periodic
       if (ibrank .and. ierank) then
#if defined (_GPU)
         !$acc kernels default(present)
         ekm_d(ib - 1, :, :) = ekm_d(ie, :, :)
         ekm_d(ie + 1, :, :) = ekm_d(ib, :, :)
         ekh_d(ib - 1, :, :) = ekh_d(ie, :, :)
         ekh_d(ie + 1, :, :) = ekh_d(ib, :, :)
         !$acc end kernels
#else
         ekm(ib - 1, :, :) = ekm(ie, :, :)
         ekm(ie + 1, :, :) = ekm(ib, :, :)
         ekh(ib - 1, :, :) = ekh(ie, :, :)
         ekh(ie + 1, :, :) = ekh(ib, :, :)
#endif
       end if
     end if

     if (BCym .ne. BCym_periodic) then ! inflow/outflow
       if (jbrank) then
#if defined (_GPU)
         !$acc kernels default(present)
         ekm_d(:,jb-1,:) = ekm_d(:,jb,:)
         ekh_d(:,jb-1,:) = ekh_d(:,jb,:)
         !$acc end kernels
#else
         ekm(:,jb-1,:) = ekm(:,jb,:)
         ekh(:,jb-1,:) = ekh(:,jb,:)
#endif
       end if
       if (jerank) then
#if defined (_GPU)
         !$acc kernels default(present)
         ekm_d(:,je+1,:) = ekm_d(:,je,:)
         ekh_d(:,je+1,:) = ekh_d(:,je,:)
         !$acc end kernels
#else
         ekm(:,je+1,:) = ekm(:,je,:)
         ekh(:,je+1,:) = ekh(:,je,:)
#endif
       end if
     else ! periodic
       if (jbrank .and. jerank) then
#if defined (_GPU)
         !$acc kernels default(present)
         ekm_d(:, jb - 1, :) = ekm_d(:, je, :)
         ekm_d(:, je + 1, :) = ekm_d(:, jb, :)
         ekh_d(:, jb - 1, :) = ekh_d(:, je, :)
         ekh_d(:, je + 1, :) = ekh_d(:, jb, :)
         !$acc end kernels
#else
         ekm(:, jb - 1, :) = ekm(:, je, :)
         ekm(:, je + 1, :) = ekm(:, jb, :)
         ekh(:, jb - 1, :) = ekh(:, je, :)
         ekh(:, je + 1, :) = ekh(:, jb, :)
#endif
       end if
     end if

     call reassure_fluxtop_boundary

   end subroutine closurebc


   !> Set lateral periodic boundary conditions for momentum in x/i direction on the host
   subroutine xm_periodic
      use modglobal,      only : ib, ie, ih
      use modfields,      only : u0, um, v0, vm, w0, wm, e120, e12m
      use modsubgriddata, only : loneeqn
      implicit none
      integer :: m

      do m = 1, ih
         u0(ib - m, :, :) = u0(ie + 1 - m, :, :)
         u0(ie + m, :, :) = u0(ib - 1 + m, :, :)
         v0(ib - m, :, :) = v0(ie + 1 - m, :, :)
         v0(ie + m, :, :) = v0(ib - 1 + m, :, :)
         w0(ib - m, :, :) = w0(ie + 1 - m, :, :)
         w0(ie + m, :, :) = w0(ib - 1 + m, :, :)
         um(ib - m, :, :) = um(ie + 1 - m, :, :)
         um(ie + m, :, :) = um(ib - 1 + m, :, :)
         vm(ib - m, :, :) = vm(ie + 1 - m, :, :)
         vm(ie + m, :, :) = vm(ib - 1 + m, :, :)
         wm(ib - m, :, :) = wm(ie + 1 - m, :, :)
         wm(ie + m, :, :) = wm(ib - 1 + m, :, :)
      end do

      if (loneeqn) then
         do m = 1, ih
            e120(ib - m, :, :) = e120(ie + 1 - m, :, :)
            e120(ie + m, :, :) = e120(ib - 1 + m, :, :)
            e12m(ib - m, :, :) = e12m(ie + 1 - m, :, :)
            e12m(ie + m, :, :) = e12m(ib - 1 + m, :, :)
         end do
      end if
   end subroutine xm_periodic
#if defined(_GPU)
   subroutine xm_periodic_device
      use modglobal,      only : ib, ie, ih, jb, je, jh, kb, ke, kh
      use modsubgriddata, only : loneeqn
      use modcuda,        only : u0_d, um_d, v0_d, vm_d, w0_d, wm_d, e120_d, e12m_d
      implicit none
      integer :: j, k, m

      !$acc parallel loop collapse(2) default(present) private(m)
      do k = kb - kh, ke + kh
         do j = jb - jh, je + jh
            do m = 1, ih
               u0_d(ib - m, j, k) = u0_d(ie + 1 - m, j, k)
               u0_d(ie + m, j, k) = u0_d(ib - 1 + m, j, k)
               v0_d(ib - m, j, k) = v0_d(ie + 1 - m, j, k)
               v0_d(ie + m, j, k) = v0_d(ib - 1 + m, j, k)
               w0_d(ib - m, j, k) = w0_d(ie + 1 - m, j, k)
               w0_d(ie + m, j, k) = w0_d(ib - 1 + m, j, k)
               um_d(ib - m, j, k) = um_d(ie + 1 - m, j, k)
               um_d(ie + m, j, k) = um_d(ib - 1 + m, j, k)
               vm_d(ib - m, j, k) = vm_d(ie + 1 - m, j, k)
               vm_d(ie + m, j, k) = vm_d(ib - 1 + m, j, k)
               wm_d(ib - m, j, k) = wm_d(ie + 1 - m, j, k)
               wm_d(ie + m, j, k) = wm_d(ib - 1 + m, j, k)
            end do
         end do
      end do
      !$acc end parallel loop

      if (loneeqn) then
         !$acc parallel loop collapse(2) default(present) private(m)
         do k = kb - kh, ke + kh
            do j = jb - jh, je + jh
               do m = 1, ih
                  e120_d(ib - m, j, k) = e120_d(ie + 1 - m, j, k)
                  e120_d(ie + m, j, k) = e120_d(ib - 1 + m, j, k)
                  e12m_d(ib - m, j, k) = e12m_d(ie + 1 - m, j, k)
                  e12m_d(ie + m, j, k) = e12m_d(ib - 1 + m, j, k)
               end do
            end do
         end do
         !$acc end parallel loop
      end if
   end subroutine xm_periodic_device
#endif

   !> Set x/i periodic boundary conditions for temperature on the host
   subroutine xT_periodic
      use modglobal, only : ib, ie, ih, ihc
      use modfields, only : thl0, thlm, thl0c
      implicit none
      integer :: m

      do m = 1, ih
         thl0(ib - m, :, :) = thl0(ie + 1 - m, :, :)
         thl0(ie + m, :, :) = thl0(ib - 1 + m, :, :)
         thlm(ib - m, :, :) = thlm(ie + 1 - m, :, :)
         thlm(ie + m, :, :) = thlm(ib - 1 + m, :, :)
      end do

      do m = 1, ihc
         thl0c(ib - m, :, :) = thl0c(ie + 1 - m, :, :)
         thl0c(ie + m, :, :) = thl0c(ib - 1 + m, :, :)
      end do
   end subroutine xT_periodic
#if defined(_GPU)
   subroutine xT_periodic_device
      use modglobal, only : jb, je, jh, jhc, kb, ke, kh, khc, &
                            ib, ie, ih, ihc, ltempeq, iadv_thl, iadv_kappa
      use modcuda,   only : thl0_d, thlm_d, thl0c_d
      implicit none
      integer :: j, k, m

      if (ltempeq) then
        !$acc parallel loop collapse(2) default(present) private(m)
        do k = kb - kh, ke + kh
          do j = jb - jh, je + jh
            do m = 1, ih
              thl0_d(ib - m, j, k) = thl0_d(ie + 1 - m, j, k)
              thl0_d(ie + m, j, k) = thl0_d(ib - 1 + m, j, k)
              thlm_d(ib - m, j, k) = thlm_d(ie + 1 - m, j, k)
              thlm_d(ie + m, j, k) = thlm_d(ib - 1 + m, j, k)
            end do
          end do
        end do
        !$acc end parallel loop

        if (iadv_thl == iadv_kappa) then
          !$acc parallel loop collapse(2) default(present) private(m)
          do k = kb - khc, ke + khc
            do j = jb - jhc, je + jhc
              do m = 1, ihc
                thl0c_d(ib - m, j, k) = thl0c_d(ie + 1 - m, j, k)
                thl0c_d(ie + m, j, k) = thl0c_d(ib - 1 + m, j, k)
              end do
            end do
          end do
          !$acc end parallel loop
        end if
      end if
   end subroutine xT_periodic_device
#endif

   !> Set x/i periodic boundary conditions for humidity on the host
   subroutine xq_periodic
      use modglobal, only : ib, ie, ih
      use modfields, only : qt0, qtm
      implicit none
      integer :: m

      do m = 1, ih
         qt0(ib - m, :, :) = qt0(ie + 1 - m, :, :)
         qt0(ie + m, :, :) = qt0(ib - 1 + m, :, :)
         qtm(ib - m, :, :) = qtm(ie + 1 - m, :, :)
         qtm(ie + m, :, :) = qtm(ib - 1 + m, :, :)
      end do
   end subroutine xq_periodic
#if defined(_GPU)
   subroutine xq_periodic_device
      use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh, lmoist
      use modcuda,   only : qt0_d, qtm_d
      implicit none
      integer :: j, k, m

      if (lmoist) then
        !$acc parallel loop collapse(2) default(present) private(m)
        do k = kb - kh, ke + kh
          do j = jb - jh, je + jh
            do m = 1, ih
              qt0_d(ib - m, j, k) = qt0_d(ie + 1 - m, j, k)
              qt0_d(ie + m, j, k) = qt0_d(ib - 1 + m, j, k)
              qtm_d(ib - m, j, k) = qtm_d(ie + 1 - m, j, k)
              qtm_d(ie + m, j, k) = qtm_d(ib - 1 + m, j, k)
            end do
          end do
        end do
        !$acc end parallel loop
      end if
   end subroutine xq_periodic_device
#endif

   !> Set x/i periodic boundary conditions for scalars on the host
   subroutine xs_periodic
      use modglobal, only : ib, ie, ihc
      use modfields, only : sv0, svm
      implicit none
      integer :: m

      do m = 1, ihc
         sv0(ib - m, :, :, :) = sv0(ie + 1 - m, :, :, :)
         sv0(ie + m, :, :, :) = sv0(ib - 1 + m, :, :, :)
         svm(ib - m, :, :, :) = svm(ie + 1 - m, :, :, :)
         svm(ie + m, :, :, :) = svm(ib - 1 + m, :, :, :)
      end do
   end subroutine xs_periodic
#if defined(_GPU)
   subroutine xs_periodic_device
      use modglobal, only : ib, ie, ihc, jb, je, jhc, kb, ke, khc, nsv
      use modcuda,   only : sv0_d, svm_d
      implicit none
      integer :: j, k, m, n

      if (nsv>0) then
        !$acc parallel loop collapse(3) default(present) private(m)
        do n = 1, nsv
          do k = kb - khc, ke + khc
            do j = jb - jhc, je + jhc
              do m = 1, ihc
                sv0_d(ib - m, j, k, n) = sv0_d(ie + 1 - m, j, k, n)
                sv0_d(ie + m, j, k, n) = sv0_d(ib - 1 + m, j, k, n)
                svm_d(ib - m, j, k, n) = svm_d(ie + 1 - m, j, k, n)
                svm_d(ie + m, j, k, n) = svm_d(ib - 1 + m, j, k, n)
              end do
            end do
          end do
        end do
        !$acc end parallel loop
      end if
   end subroutine xs_periodic_device
#endif

   !> Set lateral periodic boundary conditions for momentum in y/j direction on the host
   subroutine ym_periodic
      use modglobal,      only : jb, je, jh
      use modfields,      only : u0, um, v0, vm, w0, wm, e120, e12m
      use modsubgriddata, only : loneeqn
      implicit none
      integer :: m

      do m = 1, jh
         u0(:, jb - m, :) = u0(:, je + 1 - m, :)
         u0(:, je + m, :) = u0(:, jb - 1 + m, :)
         v0(:, jb - m, :) = v0(:, je + 1 - m, :)
         v0(:, je + m, :) = v0(:, jb - 1 + m, :)
         w0(:, jb - m, :) = w0(:, je + 1 - m, :)
         w0(:, je + m, :) = w0(:, jb - 1 + m, :)
         um(:, jb - m, :) = um(:, je + 1 - m, :)
         um(:, je + m, :) = um(:, jb - 1 + m, :)
         vm(:, jb - m, :) = vm(:, je + 1 - m, :)
         vm(:, je + m, :) = vm(:, jb - 1 + m, :)
         wm(:, jb - m, :) = wm(:, je + 1 - m, :)
         wm(:, je + m, :) = wm(:, jb - 1 + m, :)
      end do

      if (loneeqn) then
         do m = 1, jh
            e120(:, jb - m, :) = e120(:, je + 1 - m, :)
            e120(:, je + m, :) = e120(:, jb - 1 + m, :)
            e12m(:, jb - m, :) = e12m(:, je + 1 - m, :)
            e12m(:, je + m, :) = e12m(:, jb - 1 + m, :)
         end do
      end if
   end subroutine ym_periodic
#if defined(_GPU)
   subroutine ym_periodic_device
      use modglobal,      only : ib, ie, ih, jb, je, jh, kb, ke, kh
      use modsubgriddata, only : loneeqn
      use modcuda,        only : u0_d, um_d, v0_d, vm_d, w0_d, wm_d, e120_d, e12m_d
      implicit none
      integer :: i, k, m

      !$acc parallel loop collapse(2) default(present) private(m)
      do k = kb - kh, ke + kh
         do i = ib - ih, ie + ih
            do m = 1, jh
               u0_d(i, jb - m, k) = u0_d(i, je + 1 - m, k)
               u0_d(i, je + m, k) = u0_d(i, jb - 1 + m, k)
               v0_d(i, jb - m, k) = v0_d(i, je + 1 - m, k)
               v0_d(i, je + m, k) = v0_d(i, jb - 1 + m, k)
               w0_d(i, jb - m, k) = w0_d(i, je + 1 - m, k)
               w0_d(i, je + m, k) = w0_d(i, jb - 1 + m, k)
               um_d(i, jb - m, k) = um_d(i, je + 1 - m, k)
               um_d(i, je + m, k) = um_d(i, jb - 1 + m, k)
               vm_d(i, jb - m, k) = vm_d(i, je + 1 - m, k)
               vm_d(i, je + m, k) = vm_d(i, jb - 1 + m, k)
               wm_d(i, jb - m, k) = wm_d(i, je + 1 - m, k)
               wm_d(i, je + m, k) = wm_d(i, jb - 1 + m, k)
            end do
         end do
      end do
      !$acc end parallel loop

      if (loneeqn) then
         !$acc parallel loop collapse(2) default(present) private(m)
         do k = kb - kh, ke + kh
            do i = ib - ih, ie + ih
               do m = 1, jh
                  e120_d(i, jb - m, k) = e120_d(i, je + 1 - m, k)
                  e120_d(i, je + m, k) = e120_d(i, jb - 1 + m, k)
                  e12m_d(i, jb - m, k) = e12m_d(i, je + 1 - m, k)
                  e12m_d(i, je + m, k) = e12m_d(i, jb - 1 + m, k)
               end do
            end do
         end do
         !$acc end parallel loop
      end if
   end subroutine ym_periodic_device
#endif

   !> Set y/j periodic boundary conditions for temperature on the host
   subroutine yT_periodic
      use modglobal, only : jb, je, jh, jhc
      use modfields, only : thl0, thlm, thl0c
      implicit none
      integer :: m

      do m = 1, jh
         thl0(:, jb - m, :) = thl0(:, je + 1 - m, :)
         thl0(:, je + m, :) = thl0(:, jb - 1 + m, :)
         thlm(:, jb - m, :) = thlm(:, je + 1 - m, :)
         thlm(:, je + m, :) = thlm(:, jb - 1 + m, :)
      end do

      do m = 1, jhc
         thl0c(:, jb - m, :) = thl0c(:, je + 1 - m, :)
         thl0c(:, je + m, :) = thl0c(:, jb - 1 + m, :)
      end do
   end subroutine yT_periodic
#if defined(_GPU)
   subroutine yT_periodic_device
      use modglobal, only : ib, ie, ih, ihc, kb, ke, kh, khc, &
                            jb, je, jh, jhc, ltempeq, iadv_thl, iadv_kappa
      use modcuda,   only : thl0_d, thlm_d, thl0c_d
      implicit none
      integer :: i, k, m

      if (ltempeq) then
         !$acc parallel loop collapse(2) default(present) private(m)
         do k = kb - kh, ke + kh
            do i = ib - ih, ie + ih
               do m = 1, jh
                  thl0_d(i, jb - m, k) = thl0_d(i, je + 1 - m, k)
                  thl0_d(i, je + m, k) = thl0_d(i, jb - 1 + m, k)
                  thlm_d(i, jb - m, k) = thlm_d(i, je + 1 - m, k)
                  thlm_d(i, je + m, k) = thlm_d(i, jb - 1 + m, k)
               end do
            end do
         end do
         !$acc end parallel loop

         if (iadv_thl == iadv_kappa) then
          !$acc parallel loop collapse(2) default(present) private(m)
          do k = kb - khc, ke + khc
              do i = ib - ihc, ie + ihc
                do m = 1, jhc
                    thl0c_d(i, jb - m, k) = thl0c_d(i, je + 1 - m, k)
                    thl0c_d(i, je + m, k) = thl0c_d(i, jb - 1 + m, k)
                end do
              end do
          end do
          !$acc end parallel loop
         end if
      end if
   end subroutine yT_periodic_device
#endif

   !> Set y/j periodic boundary conditions for humidity
   subroutine yq_periodic
      use modglobal, only : jb, je, jh
      use modfields, only : qt0, qtm
      implicit none
      integer :: m

      do m = 1, jh
         qt0(:, jb - m, :) = qt0(:, je + 1 - m, :)
         qt0(:, je + m, :) = qt0(:, jb - 1 + m, :)
         qtm(:, jb - m, :) = qtm(:, je + 1 - m, :)
         qtm(:, je + m, :) = qtm(:, jb - 1 + m, :)
      end do
   end subroutine yq_periodic
#if defined(_GPU)
   subroutine yq_periodic_device
      use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh, lmoist
      use modcuda,   only : qt0_d, qtm_d
      implicit none
      integer :: i, k, m

      if (lmoist) then
         !$acc parallel loop collapse(2) default(present) private(m)
         do k = kb - kh, ke + kh
            do i = ib - ih, ie + ih
               do m = 1, jh
                  qt0_d(i, jb - m, k) = qt0_d(i, je + 1 - m, k)
                  qt0_d(i, je + m, k) = qt0_d(i, jb - 1 + m, k)
                  qtm_d(i, jb - m, k) = qtm_d(i, je + 1 - m, k)
                  qtm_d(i, je + m, k) = qtm_d(i, jb - 1 + m, k)
               end do
            end do
         end do
         !$acc end parallel loop
      end if
   end subroutine yq_periodic_device
#endif

   !> Set y/j periodic boundary conditions for scalars
   subroutine ys_periodic
      use modglobal, only : jb, je, jhc, nsv
      use modfields, only : sv0, svm
      implicit none
      integer :: n, m

      do n = 1, nsv
         do m = 1, jhc
            sv0(:, jb - m, :, :) = sv0(:, je + 1 - m, :, :)
            sv0(:, je + m, :, :) = sv0(:, jb - 1 + m, :, :)
            svm(:, jb - m, :, :) = svm(:, je + 1 - m, :, :)
            svm(:, je + m, :, :) = svm(:, jb - 1 + m, :, :)
         end do
      end do
   end subroutine ys_periodic
#if defined(_GPU)
   subroutine ys_periodic_device
      use modglobal, only : ib, ie, ihc, jb, je, jhc, kb, ke, khc, nsv
      use modcuda,   only : sv0_d, svm_d
      implicit none
      integer :: i, k, m, n

      if (nsv>0) then
         !$acc parallel loop collapse(3) default(present) private(m)
         do n = 1, nsv
            do k = kb - khc, ke + khc
               do i = ib - ihc, ie + ihc
                  do m = 1, jhc
                     sv0_d(i, jb - m, k, n) = sv0_d(i, je + 1 - m, k, n)
                     sv0_d(i, je + m, k, n) = sv0_d(i, jb - 1 + m, k, n)
                     svm_d(i, jb - m, k, n) = svm_d(i, je + 1 - m, k, n)
                     svm_d(i, je + m, k, n) = svm_d(i, jb - 1 + m, k, n)
                  end do
               end do
            end do
         end do
         !$acc end parallel loop
      end if
   end subroutine ys_periodic_device
#endif


     subroutine xmi_profile
       use modglobal,      only : ib, jb, je, kb, ke
       use modfields,      only : u0, um, v0, vm, w0, wm, e120, e12m, uprof, vprof, e12prof
       use modsubgriddata, only : loneeqn

       integer j, k

       do j = jb - 1, je + 1
         do k = kb, ke + 1
           u0(ib, j, k) = uprof(k)
           um(ib, j, k) = uprof(k)
           u0(ib - 1, j, k) = 2*u0(ib, j, k) - u0(ib + 1, j, k) ! (u(ib+1)+u(ib-1))/2 = u(ib)
           um(ib - 1, j, k) = 2*um(ib, j, k) - um(ib + 1, j, k) ! (u(ib+1)+u(ib-1))/2 = u(ib)
           v0(ib - 1, j, k) = 2*vprof(k) - v0(ib, j, k) ! (v(ib)+v(ib-1))/2 = vprof
           vm(ib - 1, j, k) = 2*vprof(k) - vm(ib, j, k) ! (v(ib)+v(ib-1))/2 = vprof
           w0(ib - 1, j, k) = -w0(ib, j, k)
           wm(ib - 1, j, k) = -wm(ib, j, k)
         end do
       end do

       if (loneeqn) then
         do j = jb - 1, je + 1
           do k = kb, ke + 1
             e120(ib - 1, j, k) = 2*e12prof(k) - e120(ib, j, k) ! (e12(ib)+e12(ib-1))/2=e12prof
             e12m(ib - 1, j, k) = 2*e12prof(k) - e12m(ib, j, k) ! (e12(ib)+e12(ib-1))/2=e12prof
           end do
         end do
       end if

     end subroutine xmi_profile
#if defined(_GPU)
     !> The loops are nested the other way round from the host routine.
     !!
     !! Every lateral condition writes a plane at one fixed i or j of an array
     !! whose i is contiguous, so which index the threads of a warp walk
     !! decides how far apart their addresses are. For an x plane that index
     !! has to be j - a stride of one i-row - rather than k, which is a whole
     !! i-j slab away. Getting this backwards is not a small penalty: the
     !! scalar inlet below cost 23 ms a stage with k innermost and a fraction
     !! of a millisecond with j innermost.
     subroutine xmi_profile_device
       use modglobal,      only : ib, jb, je, kb, ke
       use modsubgriddata, only : loneeqn
       use modcuda,        only : u0_d, um_d, v0_d, vm_d, w0_d, wm_d, e120_d, e12m_d, &
                                  uprof_d, vprof_d, e12prof_d
       implicit none
       integer j, k

       !$acc parallel loop collapse(2) default(present)
       do k = kb, ke + 1
         do j = jb - 1, je + 1
           u0_d(ib, j, k) = uprof_d(k)
           um_d(ib, j, k) = uprof_d(k)
           u0_d(ib - 1, j, k) = 2*u0_d(ib, j, k) - u0_d(ib + 1, j, k) ! (u(ib+1)+u(ib-1))/2 = u(ib)
           um_d(ib - 1, j, k) = 2*um_d(ib, j, k) - um_d(ib + 1, j, k) ! (u(ib+1)+u(ib-1))/2 = u(ib)
           v0_d(ib - 1, j, k) = 2*vprof_d(k) - v0_d(ib, j, k) ! (v(ib)+v(ib-1))/2 = vprof
           vm_d(ib - 1, j, k) = 2*vprof_d(k) - vm_d(ib, j, k) ! (v(ib)+v(ib-1))/2 = vprof
           w0_d(ib - 1, j, k) = -w0_d(ib, j, k)
           wm_d(ib - 1, j, k) = -wm_d(ib, j, k)
         end do
       end do
       !$acc end parallel loop

       if (loneeqn) then
         !$acc parallel loop collapse(2) default(present)
         do k = kb, ke + 1
           do j = jb - 1, je + 1
             e120_d(ib - 1, j, k) = 2*e12prof_d(k) - e120_d(ib, j, k) ! (e12(ib)+e12(ib-1))/2=e12prof
             e12m_d(ib - 1, j, k) = 2*e12prof_d(k) - e12m_d(ib, j, k) ! (e12(ib)+e12(ib-1))/2=e12prof
           end do
         end do
         !$acc end parallel loop
       end if

     end subroutine xmi_profile_device
#endif


     subroutine xmi_driver
       use modglobal,      only : ib, jb, je, kb, ke
       use modinletdata,   only : u0driver, umdriver, v0driver, vmdriver, w0driver, wmdriver
       use modfields,      only : u0, um, v0, vm, w0, wm, e120, e12m
       use modsubgriddata, only : loneeqn

       integer j, k

       do j = jb - 1, je + 1
         do k = kb, ke !tg3315 removed +1 following above...
           u0(ib,j,k) = u0driver(j,k) !max(0.,u0driver(j,k))
           um(ib,j,k) = umdriver(j,k) !max(0.,umdriver(j,k))
           u0(ib-1,j,k) = u0driver(j,k) !max(0.,u0driver(j,k))
           um(ib-1,j,k) = umdriver(j,k) !max(0.,umdriver(j,k))
           ! u0(ib-1,j,k)= 2*u0(ib, j, k) - u0(ib + 1, j, k) ! (u(ib+1)+u(ib-1))/2 = u(ib)
           ! um(ib-1,j,k)= 2*um(ib, j, k) - um(ib + 1, j, k) ! (u(ib+1)+u(ib-1))/2 = u(ib)

           !v0(ib,j,k)   = v0driver(j,k) !max(0.,v0driver(j,k))
           !vm(ib,j,k)   = vmdriver(j,k) !max(0.,vmdriver(j,k))
           v0(ib-1,j,k)   = v0driver(j,k) !max(0.,v0driver(j,k))
           vm(ib-1,j,k)   = vmdriver(j,k) !max(0.,vmdriver(j,k))
         end do

         do k=kb,ke+1
           !w0(ib,j,k)   = w0driver(j,k) !max(0.,w0driver(j,k))
           !wm(ib,j,k)   = wmdriver(j,k) !max(0.,wmdriver(j,k))
           w0(ib-1,j,k) = w0driver(j,k) !max(0.,w0driver(j,k))
           wm(ib-1,j,k) = wmdriver(j,k) !max(0.,wmdriver(j,k))
         end do
       end do

       if (loneeqn) then
         do j = jb - 1, je + 1
           do k = kb, ke + 1
             ! to be changed in the future: e12 should be taken from recycle plane!
             !e120(ib-1,j,k) = e120driver(j,k)      ! extrapolate e12 from interior
             !e12m(ib-1,j,k) = e12mdriver(j,k)      ! extrapolate e12 from interior
             e120(ib - 1, j, k) = e120(ib, j, k) ! (e12(ib)+e12(ib-1))/2=e12prof
             e12m(ib - 1, j, k) = e12m(ib, j, k) ! (e12(ib)+e12(ib-1))/2=e12prof
           end do
         end do
       end if

     end subroutine xmi_driver
#if defined(_GPU)
     subroutine xmi_driver_device
       use modglobal,      only : ib, jb, je, kb, ke
       use modsubgriddata, only : loneeqn
       use modcuda,        only : u0_d, um_d, v0_d, vm_d, w0_d, wm_d, e120_d, e12m_d, &
                                  u0driver_d, umdriver_d, v0driver_d, vmdriver_d, &
                                  w0driver_d, wmdriver_d
       implicit none
       integer j, k

       !$acc parallel loop collapse(2) default(present)
       do k = kb, ke !tg3315 removed +1 following above...
         do j = jb - 1, je + 1
           u0_d(ib,j,k) = u0driver_d(j,k)
           um_d(ib,j,k) = umdriver_d(j,k)
           u0_d(ib-1,j,k) = u0driver_d(j,k)
           um_d(ib-1,j,k) = umdriver_d(j,k)

           v0_d(ib-1,j,k) = v0driver_d(j,k)
           vm_d(ib-1,j,k) = vmdriver_d(j,k)
         end do
       end do
       !$acc end parallel loop

       !$acc parallel loop collapse(2) default(present)
       do k = kb, ke + 1
         do j = jb - 1, je + 1
           w0_d(ib-1,j,k) = w0driver_d(j,k)
           wm_d(ib-1,j,k) = wmdriver_d(j,k)
         end do
       end do
       !$acc end parallel loop

       if (loneeqn) then
         !$acc parallel loop collapse(2) default(present)
         do k = kb, ke + 1
           do j = jb - 1, je + 1
             ! to be changed in the future: e12 should be taken from recycle plane!
             e120_d(ib - 1, j, k) = e120_d(ib, j, k)
             e12m_d(ib - 1, j, k) = e12m_d(ib, j, k)
           end do
         end do
         !$acc end parallel loop
       end if

     end subroutine xmi_driver_device
#endif


     subroutine xTi_profile
       use modglobal, only : ib, jb, je, kb, ke
       use modfields, only : thl0, thlm, thlprof
       integer j, k

       ! set ghost cell
       ! do j = jb - 1, je + 1
       !   do k = kb, ke + 1
       !     thl0(ib - 1, j, k) = 2*thlprof(k) - thl0(ib, j, k)
       !     thlm(ib - 1, j, k) = 2*thlprof(k) - thlm(ib, j, k)
       !   end do
       ! end do
       do j = jb - 1, je + 1
         do k = kb, ke + 1
           thl0(ib - 1, j, k) = thlprof(k)
           thlm(ib - 1, j, k) = thlprof(k)
         end do
       end do

       ! set first internal cell as well
       do j = jb - 1, je + 1
        do k = kb, ke
           thl0(ib, j, k) = thlprof(k)
           thlm(ib, j, k) = thlprof(k)
        end do
       end do

     end subroutine xTi_profile
#if defined(_GPU)
     subroutine xTi_profile_device
       use modglobal, only : ib, jb, je, kb, ke
       use modcuda,   only : thl0_d, thlm_d, thlprof_d
       implicit none
       integer j, k

       ! set ghost cell
       !$acc parallel loop collapse(2) default(present)
       do k = kb, ke + 1
         do j = jb - 1, je + 1
           thl0_d(ib - 1, j, k) = thlprof_d(k)
           thlm_d(ib - 1, j, k) = thlprof_d(k)
         end do
       end do
       !$acc end parallel loop

       ! set first internal cell as well
       !$acc parallel loop collapse(2) default(present)
       do k = kb, ke
         do j = jb - 1, je + 1
           thl0_d(ib, j, k) = thlprof_d(k)
           thlm_d(ib, j, k) = thlprof_d(k)
         end do
       end do
       !$acc end parallel loop

     end subroutine xTi_profile_device
#endif


     subroutine xTi_driver
       use modglobal,    only : ib, jb, je, kb, ke
       use modinletdata, only : thl0driver, thlmdriver
       use modfields,    only : thl0, thlm
       integer j, k

       do j = jb - 1, je + 1
         do k = kb, ke + 1
           thl0(ib - 1, j, k) = thl0driver(j, k)
           thlm(ib - 1, j, k) = thlmdriver(j, k)
         end do
       end do

     end subroutine xTi_driver
#if defined(_GPU)
     subroutine xTi_driver_device
       use modglobal, only : ib, jb, je, kb, ke
       use modcuda,   only : thl0_d, thlm_d, thl0driver_d, thlmdriver_d
       implicit none
       integer j, k

       !$acc parallel loop collapse(2) default(present)
       do k = kb, ke + 1
         do j = jb - 1, je + 1
           thl0_d(ib - 1, j, k) = thl0driver_d(j, k)
           thlm_d(ib - 1, j, k) = thlmdriver_d(j, k)
         end do
       end do
       !$acc end parallel loop

     end subroutine xTi_driver_device
#endif


     subroutine xqi_profile
       use modglobal,    only : ib, jb, je, kb, ke
       use modfields,    only : qt0, qtm, qtprof
       integer j, k

       do j = jb - 1, je + 1
         do k = kb, ke + 1
           qt0(ib - 1, j, k) = 2*qtprof(k) - qt0(ib, j, k)
           qtm(ib - 1, j, k) = 2*qtprof(k) - qtm(ib, j, k)
         end do
       end do

   end subroutine xqi_profile
#if defined(_GPU)
   subroutine xqi_profile_device
     use modglobal, only : ib, jb, je, kb, ke
     use modcuda,   only : qt0_d, qtm_d, qtprof_d
     implicit none
     integer j, k

     !$acc parallel loop collapse(2) default(present)
     do k = kb, ke + 1
       do j = jb - 1, je + 1
         qt0_d(ib - 1, j, k) = 2*qtprof_d(k) - qt0_d(ib, j, k)
         qtm_d(ib - 1, j, k) = 2*qtprof_d(k) - qtm_d(ib, j, k)
       end do
     end do
     !$acc end parallel loop

   end subroutine xqi_profile_device
#endif


   subroutine xqi_driver
     use modglobal,    only : ib, jb, je, kb, ke
     use modinletdata, only : qt0driver, qtmdriver
     use modfields,    only : qt0, qtm

     integer j, k

     do j = jb - 1, je + 1
       do k = kb, ke + 1
         qt0(ib - 1, j, k) = qt0driver(j, k)
         qtm(ib - 1, j, k) = qtmdriver(j, k)
       end do
     end do

   end subroutine xqi_driver
#if defined(_GPU)
   subroutine xqi_driver_device
     use modglobal, only : ib, jb, je, kb, ke
     use modcuda,   only : qt0_d, qtm_d, qt0driver_d, qtmdriver_d
     implicit none
     integer j, k

     !$acc parallel loop collapse(2) default(present)
     do k = kb, ke + 1
       do j = jb - 1, je + 1
         qt0_d(ib - 1, j, k) = qt0driver_d(j, k)
         qtm_d(ib - 1, j, k) = qtmdriver_d(j, k)
       end do
     end do
     !$acc end parallel loop

   end subroutine xqi_driver_device
#endif


   subroutine xsi_profile
     use modglobal,    only : ib, jb, je, kb, ke, nsv, ihc
     use modfields,    only : sv0, svm, svprof

     integer j, k, n, m

     do j = jb, je
       do k = kb, ke + 1
         do n = 1, nsv
           do m = 1, ihc
             sv0(ib - m, j, k, n) = 2*svprof(k, n) - sv0(ib - m + 1, j, k, n)
             svm(ib - m, j, k, n) = 2*svprof(k, n) - svm(ib - m + 1, j, k, n)
           end do
         end do
       end do
     end do

   end subroutine xsi_profile
#if defined(_GPU)
   !> Scalar inlet profile: ghost columns from the first interior column.
   !!
   !! The recurrence over m is carried in a register rather than through
   !! sv0_d. The host routine reads back the column it wrote on the previous
   !! pass, which on the device is a dependent round trip to global memory at
   !! a fixed i - a stride of one i-row between neighbouring threads, so every
   !! lane of a warp waits on its own transaction. Reading the interior column
   !! once and stepping the recurrence locally is the same arithmetic in the
   !! same order, and it took this routine from 23 ms a stage to under one.
   subroutine xsi_profile_device
     use modglobal, only : ib, jb, je, kb, ke, nsv, ihc
     use modcuda,   only : sv0_d, svm_d, svprof_d
     implicit none
     integer j, k, n, m
     real    :: c0, cm, twoprof

     !$acc parallel loop collapse(3) default(present) private(m, c0, cm, twoprof)
     do n = 1, nsv
       do k = kb, ke + 1
         do j = jb, je
           twoprof = 2*svprof_d(k, n)
           c0 = sv0_d(ib, j, k, n)
           cm = svm_d(ib, j, k, n)
           do m = 1, ihc
             c0 = twoprof - c0
             cm = twoprof - cm
             sv0_d(ib - m, j, k, n) = c0
             svm_d(ib - m, j, k, n) = cm
           end do
         end do
       end do
     end do
     !$acc end parallel loop

   end subroutine xsi_profile_device
#endif


     subroutine xsi_custom
       use modglobal,    only : ib, jb, je, jtot, kb, ke, nsv, ihc
       use modfields,    only : sv0, svm, svprof
       use decomp_2d,    only : zstart

       integer j, k, n, m

       do j = jb, je
          if (j + zstart(2) - 1 == jtot/2) then
             do k = kb, ke + 1
                do n = 1, nsv
                   do m = 1, ihc
                      sv0(ib - m, j-1:j+1, k, n) = 2*svprof(k, n) - sv0(ib - m + 1, j-1:j+1, k, n)
                      svm(ib - m, j-1:j+1, k, n) = 2*svprof(k, n) - svm(ib - m + 1, j-1:j+1, k, n)
                   end do
                end do
             end do
          end if
       end do

   end subroutine xsi_custom
#if defined(_GPU)
     subroutine xsi_custom_device
       use modglobal, only : ib, jb, je, jtot, kb, ke, nsv, ihc
       use modcuda,   only : sv0_d, svm_d, svprof_d
       use decomp_2d, only : zstart
       implicit none
       integer j, k, n, m, jmid
       real    :: c0, cm, twoprof

       ! The host branch walks every j to find the one that carries the middle
       ! of the domain; the rank knows which that is without looking.
       jmid = jtot/2 - zstart(2) + 1
       if (jmid < jb .or. jmid > je) return

       ! The recurrence is carried in a register, as in xsi_profile_device.
       !$acc parallel loop collapse(3) default(present) private(m, c0, cm, twoprof)
       do n = 1, nsv
         do k = kb, ke + 1
           do j = jmid - 1, jmid + 1
             twoprof = 2*svprof_d(k, n)
             c0 = sv0_d(ib, j, k, n)
             cm = svm_d(ib, j, k, n)
             do m = 1, ihc
               c0 = twoprof - c0
               cm = twoprof - cm
               sv0_d(ib - m, j, k, n) = c0
               svm_d(ib - m, j, k, n) = cm
             end do
           end do
         end do
       end do
       !$acc end parallel loop

   end subroutine xsi_custom_device
#endif


   subroutine xsi_driver
     use modglobal,    only : ib, ihc, jb, je, kb, ke, nsv
     use modinletdata, only : sv0driver, svmdriver
     use modfields,    only : sv0, svm

     integer j, k, n, m

     do j = jb - 1, je + 1
       do k = kb, ke + 1
         do n = 1, nsv
           do m = 1, ihc
             sv0(ib - m, j, k, n) = sv0driver(j, k, n)
             svm(ib - m, j, k, n) = svmdriver(j, k, n)
           end do
         end do
       end do
     end do

   end subroutine xsi_driver
#if defined(_GPU)
   subroutine xsi_driver_device
     use modglobal, only : ib, ihc, jb, je, kb, ke, nsv
     use modcuda,   only : sv0_d, svm_d, sv0driver_d, svmdriver_d
     implicit none
     integer j, k, n, m

     !$acc parallel loop collapse(3) default(present) private(m)
     do n = 1, nsv
       do k = kb, ke + 1
         do j = jb - 1, je + 1
           do m = 1, ihc
             sv0_d(ib - m, j, k, n) = sv0driver_d(j, k, n)
             svm_d(ib - m, j, k, n) = svmdriver_d(j, k, n)
           end do
         end do
       end do
     end do
     !$acc end parallel loop

   end subroutine xsi_driver_device
#endif


   subroutine xmo_convective
     use modglobal,      only : ie, dxi, rk3coef
     use modfields,      only : v0, vm, w0, wm, e120, e12m, uouttot
     use modsubgriddata, only : loneeqn

     v0(ie + 1, :, :) = v0(ie+1, :, :) - (v0(ie+1, :, :) - v0(ie, :, :))*dxi*rk3coef*uouttot
     w0(ie + 1, :, :) = w0(ie+1, :, :) - (w0(ie+1, :, :) - w0(ie, :, :))*dxi*rk3coef*uouttot
     vm(ie + 1, :, :) = vm(ie+1, :, :) - (vm(ie+1, :, :) - vm(ie, :, :))*dxi*rk3coef*uouttot
     wm(ie + 1, :, :) = wm(ie+1, :, :) - (wm(ie+1, :, :) - wm(ie, :, :))*dxi*rk3coef*uouttot

     if (loneeqn) then
       e120(ie + 1, :, :) = e120(ie, :, :) - (e120(ie + 1, :, :) - e120(ie, :, :))*dxi*rk3coef*uouttot
       e12m(ie + 1, :, :) = e12m(ie, :, :) - (e12m(ie + 1, :, :) - e12m(ie, :, :))*dxi*rk3coef*uouttot
     end if

   end subroutine xmo_convective
#if defined(_GPU)
   subroutine xmo_convective_device
     use modglobal,      only : ie, jb, je, jh, kb, ke, kh, dxi, rk3coef
     use modfields,      only : uouttot
     use modsubgriddata, only : loneeqn
     use modcuda,        only : v0_d, vm_d, w0_d, wm_d, e120_d, e12m_d
     implicit none
     integer j, k
     real    :: l_dxi, l_rk3coef, l_uouttot

     l_dxi = dxi
     l_rk3coef = rk3coef
     l_uouttot = uouttot

     !$acc parallel loop collapse(2) default(present)
     do k = kb - kh, ke + kh
       do j = jb - jh, je + jh
         v0_d(ie + 1, j, k) = v0_d(ie+1, j, k) - (v0_d(ie+1, j, k) - v0_d(ie, j, k))*l_dxi*l_rk3coef*l_uouttot
         w0_d(ie + 1, j, k) = w0_d(ie+1, j, k) - (w0_d(ie+1, j, k) - w0_d(ie, j, k))*l_dxi*l_rk3coef*l_uouttot
         vm_d(ie + 1, j, k) = vm_d(ie+1, j, k) - (vm_d(ie+1, j, k) - vm_d(ie, j, k))*l_dxi*l_rk3coef*l_uouttot
         wm_d(ie + 1, j, k) = wm_d(ie+1, j, k) - (wm_d(ie+1, j, k) - wm_d(ie, j, k))*l_dxi*l_rk3coef*l_uouttot
       end do
     end do
     !$acc end parallel loop

     if (loneeqn) then
       !$acc parallel loop collapse(2) default(present)
       do k = kb - kh, ke + kh
         do j = jb - jh, je + jh
           e120_d(ie + 1, j, k) = e120_d(ie, j, k) - (e120_d(ie + 1, j, k) - e120_d(ie, j, k))*l_dxi*l_rk3coef*l_uouttot
           e12m_d(ie + 1, j, k) = e12m_d(ie, j, k) - (e12m_d(ie + 1, j, k) - e12m_d(ie, j, k))*l_dxi*l_rk3coef*l_uouttot
         end do
       end do
       !$acc end parallel loop
     end if

   end subroutine xmo_convective_device
#endif


  !  subroutine xmo_Neumann
  !    use modglobal,      only : ie
  !    use modfields,      only : v0, vm, w0, wm, e120, e12m
  !    use modsubgriddata, only : loneeqn

  !    v0(ie + 1, :, :) = v0(ie, :, :)
  !    w0(ie + 1, :, :) = w0(ie, :, :)
  !    vm(ie + 1, :, :) = vm(ie, :, :)
  !    wm(ie + 1, :, :) = wm(ie, :, :)

  !    if (loneeqn) then
  !      e120(ie + 1, :, :) = e120(ie, :, :)
  !      e12m(ie + 1, :, :) = e12m(ie, :, :)
  !    end if

  !  end subroutine xmo_Neumann


   subroutine xTo_convective
     use modglobal, only : ie, dxi, rk3coef
     use modfields, only : thl0, thlm, uouttot

     thl0(ie + 1, :, :) = thl0(ie+1, :, :) - (thl0(ie + 1, :, :) - thl0(ie, :, :))*dxi*rk3coef*uouttot
     thlm(ie + 1, :, :) = thlm(ie+1, :, :) - (thlm(ie + 1, :, :) - thlm(ie, :, :))*dxi*rk3coef*uouttot

   end subroutine xTo_convective
#if defined(_GPU)
   subroutine xTo_convective_device
     use modglobal, only : ie, jb, je, jh, kb, ke, kh, dxi, rk3coef
     use modfields, only : uouttot
     use modcuda,   only : thl0_d, thlm_d
     implicit none
     integer j, k
     real    :: l_dxi, l_rk3coef, l_uouttot

     l_dxi = dxi
     l_rk3coef = rk3coef
     l_uouttot = uouttot

     !$acc parallel loop collapse(2) default(present)
     do k = kb - kh, ke + kh
       do j = jb - jh, je + jh
         thl0_d(ie + 1, j, k) = thl0_d(ie+1, j, k) - (thl0_d(ie + 1, j, k) - thl0_d(ie, j, k))*l_dxi*l_rk3coef*l_uouttot
         thlm_d(ie + 1, j, k) = thlm_d(ie+1, j, k) - (thlm_d(ie + 1, j, k) - thlm_d(ie, j, k))*l_dxi*l_rk3coef*l_uouttot
       end do
     end do
     !$acc end parallel loop

   end subroutine xTo_convective_device
#endif


  !  subroutine xTo_Neumann
  !    use modglobal, only : ie
  !    use modfields, only : thl0, thlm

  !    thl0(ie + 1, :, :) = thl0(ie, :, :)
  !    thlm(ie + 1, :, :) = thlm(ie, :, :)

  !  end subroutine xTo_Neumann


   subroutine xqo_convective
     use modglobal, only : ie, dxi, rk3coef
     use modfields, only : qt0, qtm, uouttot

     qt0(ie + 1, :, :) = qt0(ie, :, :) - (qt0(ie + 1, :, :) - qt0(ie, :, :))*dxi*rk3coef*uouttot
     qtm(ie + 1, :, :) = qtm(ie, :, :) - (qtm(ie + 1, :, :) - qtm(ie, :, :))*dxi*rk3coef*uouttot

   end subroutine xqo_convective
#if defined(_GPU)
   subroutine xqo_convective_device
     use modglobal, only : ie, jb, je, jh, kb, ke, kh, dxi, rk3coef
     use modfields, only : uouttot
     use modcuda,   only : qt0_d, qtm_d
     implicit none
     integer j, k
     real    :: l_dxi, l_rk3coef, l_uouttot

     l_dxi = dxi
     l_rk3coef = rk3coef
     l_uouttot = uouttot

     !$acc parallel loop collapse(2) default(present)
     do k = kb - kh, ke + kh
       do j = jb - jh, je + jh
         qt0_d(ie + 1, j, k) = qt0_d(ie, j, k) - (qt0_d(ie + 1, j, k) - qt0_d(ie, j, k))*l_dxi*l_rk3coef*l_uouttot
         qtm_d(ie + 1, j, k) = qtm_d(ie, j, k) - (qtm_d(ie + 1, j, k) - qtm_d(ie, j, k))*l_dxi*l_rk3coef*l_uouttot
       end do
     end do
     !$acc end parallel loop

   end subroutine xqo_convective_device
#endif


   subroutine xso_convective
     use modglobal, only : ie, rk3coef, dxi, nsv
     use modfields, only :sv0, svm, uouttot
     integer n

     do n = 1, nsv
       sv0(ie + 1, :, :, n) = sv0(ie + 1, :, :, n) - (sv0(ie + 1, :, :, n) - sv0(ie, :, :, n))*dxi*rk3coef*uouttot
       svm(ie + 1, :, :, n) = svm(ie + 1, :, :, n) - (svm(ie + 1, :, :, n) - svm(ie, :, :, n))*dxi*rk3coef*uouttot
     end do

   end subroutine xso_convective
#if defined(_GPU)
   subroutine xso_convective_device
     use modglobal, only : ie, jb, je, jhc, kb, ke, khc, rk3coef, dxi, nsv
     use modfields, only : uouttot
     use modcuda,   only : sv0_d, svm_d
     implicit none
     integer j, k, n
     real    :: l_dxi, l_rk3coef, l_uouttot

     l_dxi = dxi
     l_rk3coef = rk3coef
     l_uouttot = uouttot

     !$acc parallel loop collapse(3) default(present)
     do n = 1, nsv
       do k = kb - khc, ke + khc
         do j = jb - jhc, je + jhc
           sv0_d(ie + 1, j, k, n) = sv0_d(ie + 1, j, k, n) - (sv0_d(ie + 1, j, k, n) - sv0_d(ie, j, k, n))*l_dxi*l_rk3coef*l_uouttot
           svm_d(ie + 1, j, k, n) = svm_d(ie + 1, j, k, n) - (svm_d(ie + 1, j, k, n) - svm_d(ie, j, k, n))*l_dxi*l_rk3coef*l_uouttot
         end do
       end do
     end do
     !$acc end parallel loop

   end subroutine xso_convective_device
#endif


  !  subroutine xso_Neumann
  !    use modglobal, only : ie, ihc, nsv
  !    use modfields, only :sv0, svm
  !    integer n, m

  !    do n = 1, nsv
  !      do m = 1, ihc
  !        sv0(ie + m, :, :, n) = sv0(ie, :, :, n)
  !        svm(ie + m, :, :, n) = svm(ie, :, :, n)
  !      end do
  !    end do

  !  end subroutine xso_Neumann


   subroutine ymi_profile
     use modglobal,      only : ib, ie, jb, kb, ke
     use modfields,      only : u0, um, v0, vm, w0, wm, e120, e12m, uprof, vprof, e12prof
     use modsubgriddata, only : loneeqn
     integer i, k

     do i = ib - 1, ie + 1
       do k = kb, ke + 1
         v0(i, jb, k) = vprof(k)
         vm(i, jb, k) = vprof(k)
         v0(i, jb - 1, k) = 2*v0(i, jb, k) - v0(i, jb + 1, k)
         vm(i, jb - 1, k) = 2*vm(i, jb, k) - vm(i, jb + 1, k)
         u0(i, jb - 1, k) = 2*uprof(k) - u0(i, jb, k)
         um(i, jb - 1, k) = 2*uprof(k) - um(i, jb, k)
         w0(i, jb - 1, k) = -w0(i, jb, k)
         wm(i, jb - 1, k) = -wm(i, jb, k)
       end do
     end do

     if (loneeqn) then
       do i = ib - 1, ie + 1
         do k = kb, ke + 1
           e120(i, jb - 1, k) = 2*e12prof(k) - e120(i, jb - 1, k)
           e12m(i, jb - 1, k) = 2*e12prof(k) - e12m(i, jb - 1, k)
         end do
       end do
     end if

   end subroutine ymi_profile
#if defined(_GPU)
   !> i innermost: a y plane is a fixed j, so consecutive threads walking i are
   !! consecutive in memory. See the note on xmi_profile_device.
   subroutine ymi_profile_device
     use modglobal,      only : ib, ie, jb, kb, ke
     use modsubgriddata, only : loneeqn
     use modcuda,        only : u0_d, um_d, v0_d, vm_d, w0_d, wm_d, e120_d, e12m_d, &
                                uprof_d, vprof_d, e12prof_d
     implicit none
     integer i, k

     !$acc parallel loop collapse(2) default(present)
     do k = kb, ke + 1
       do i = ib - 1, ie + 1
         v0_d(i, jb, k) = vprof_d(k)
         vm_d(i, jb, k) = vprof_d(k)
         v0_d(i, jb - 1, k) = 2*v0_d(i, jb, k) - v0_d(i, jb + 1, k)
         vm_d(i, jb - 1, k) = 2*vm_d(i, jb, k) - vm_d(i, jb + 1, k)
         u0_d(i, jb - 1, k) = 2*uprof_d(k) - u0_d(i, jb, k)
         um_d(i, jb - 1, k) = 2*uprof_d(k) - um_d(i, jb, k)
         w0_d(i, jb - 1, k) = -w0_d(i, jb, k)
         wm_d(i, jb - 1, k) = -wm_d(i, jb, k)
       end do
     end do
     !$acc end parallel loop

     if (loneeqn) then
       !$acc parallel loop collapse(2) default(present)
       do k = kb, ke + 1
         do i = ib - 1, ie + 1
           e120_d(i, jb - 1, k) = 2*e12prof_d(k) - e120_d(i, jb - 1, k)
           e12m_d(i, jb - 1, k) = 2*e12prof_d(k) - e12m_d(i, jb - 1, k)
         end do
       end do
       !$acc end parallel loop
     end if

   end subroutine ymi_profile_device
#endif


   subroutine yTi_profile
     use modglobal, only : ib, ie, jb, kb, ke
     use modfields, only : thl0, thlm, thlprof

     integer i, k

     do i = ib - 1, ie + 1
       do k = kb, ke + 1
         thl0(i, jb - 1, k) = 2*thlprof(k) - thl0(i, jb, k)
         thlm(i, jb - 1, k) = 2*thlprof(k) - thlm(i, jb, k)
       end do
     end do

   end subroutine yTi_profile
#if defined(_GPU)
   subroutine yTi_profile_device
     use modglobal, only : ib, ie, jb, kb, ke
     use modcuda,   only : thl0_d, thlm_d, thlprof_d
     implicit none
     integer i, k

     !$acc parallel loop collapse(2) default(present)
     do k = kb, ke + 1
       do i = ib - 1, ie + 1
         thl0_d(i, jb - 1, k) = 2*thlprof_d(k) - thl0_d(i, jb, k)
         thlm_d(i, jb - 1, k) = 2*thlprof_d(k) - thlm_d(i, jb, k)
       end do
     end do
     !$acc end parallel loop

   end subroutine yTi_profile_device
#endif


   subroutine yqi_profile
     use modglobal, only : ie, jb, kb, ke
     use modfields, only : qt0, qtm, qtprof

     integer i, k

     do i = jb - 1, ie + 1
       do k = kb, ke + 1
         qt0(i, jb - 1, k) = 2*qtprof(k) - qt0(i, jb, k)
         qtm(i, jb - 1, k) = 2*qtprof(k) - qtm(i, jb, k)
       end do
     end do

   end subroutine yqi_profile
#if defined(_GPU)
   subroutine yqi_profile_device
     use modglobal, only : ie, jb, kb, ke
     use modcuda,   only : qt0_d, qtm_d, qtprof_d
     implicit none
     integer i, k

     ! The i loop starts from jb - 1, as it does in the host routine. On every
     ! configuration in the test suite jb equals ib, so the two agree; it is
     ! copied rather than corrected so that the device path cannot be the one
     ! that changes an answer.
     !$acc parallel loop collapse(2) default(present)
     do k = kb, ke + 1
       do i = jb - 1, ie + 1
         qt0_d(i, jb - 1, k) = 2*qtprof_d(k) - qt0_d(i, jb, k)
         qtm_d(i, jb - 1, k) = 2*qtprof_d(k) - qtm_d(i, jb, k)
       end do
     end do
     !$acc end parallel loop

   end subroutine yqi_profile_device
#endif


   subroutine ysi_profile
     use modglobal, only : ib, ie, jb, kb, ke, nsv, ihc
     use modfields, only : sv0, svm, svprof

     integer i, k, n, m

     do i = ib - 1, ie + 1
       do k = kb, ke + 1
         do n = 1, nsv
           do m = 1, ihc
             sv0(i, jb - m, k, n) = 2*svprof(k, n) - sv0(i, jb - m + 1, k, n)
             svm(i, jb - m, k, n) = 2*svprof(k, n) - svm(i, jb - m + 1, k, n)
           end do
         end do
       end do
     end do

   end subroutine ysi_profile
#if defined(_GPU)
   subroutine ysi_profile_device
     use modglobal, only : ib, ie, jb, kb, ke, nsv, ihc
     use modcuda,   only : sv0_d, svm_d, svprof_d
     implicit none
     integer i, k, n, m
     real    :: c0, cm, twoprof

     ! The recurrence is carried in a register, as in xsi_profile_device.
     !$acc parallel loop collapse(3) default(present) private(m, c0, cm, twoprof)
     do n = 1, nsv
       do k = kb, ke + 1
         do i = ib - 1, ie + 1
           twoprof = 2*svprof_d(k, n)
           c0 = sv0_d(i, jb, k, n)
           cm = svm_d(i, jb, k, n)
           do m = 1, ihc
             c0 = twoprof - c0
             cm = twoprof - cm
             sv0_d(i, jb - m, k, n) = c0
             svm_d(i, jb - m, k, n) = cm
           end do
         end do
       end do
     end do
     !$acc end parallel loop

   end subroutine ysi_profile_device
#endif


   subroutine ymo_convective
     use modglobal,      only : je, dyi, rk3coef
     use modfields,      only : u0, um, w0, wm, e120, e12m, vouttot
     use modsubgriddata, only : loneeqn

     ! change to vouttot
     u0(:, je + 1, :) = u0(:, je + 1, :) - (u0(:, je + 1, :) - u0(:, je, :))*dyi*rk3coef*vouttot
     um(:, je + 1, :) = um(:, je + 1, :) - (um(:, je + 1, :) - um(:, je, :))*dyi*rk3coef*vouttot
     w0(:, je + 1, :) = w0(:, je + 1, :) - (w0(:, je + 1, :) - w0(:, je, :))*dyi*rk3coef*vouttot
     wm(:, je + 1, :) = wm(:, je + 1, :) - (wm(:, je + 1, :) - wm(:, je, :))*dyi*rk3coef*vouttot

     if (loneeqn) then
       e120(:, je + 1, :) = e120(:, je + 1, :) - (e120(:, je + 1, :) - e120(:, je, :))*dyi*rk3coef*vouttot
       e12m(:, je + 1, :) = e12m(:, je + 1, :) - (e12m(:, je + 1, :) - e12m(:, je, :))*dyi*rk3coef*vouttot
     end if

   end subroutine ymo_convective
#if defined(_GPU)
   subroutine ymo_convective_device
     use modglobal,      only : ib, ie, ih, je, kb, ke, kh, dyi, rk3coef
     use modfields,      only : vouttot
     use modsubgriddata, only : loneeqn
     use modcuda,        only : u0_d, um_d, w0_d, wm_d, e120_d, e12m_d
     implicit none
     integer i, k
     real    :: l_dyi, l_rk3coef, l_vouttot

     l_dyi = dyi
     l_rk3coef = rk3coef
     l_vouttot = vouttot

     !$acc parallel loop collapse(2) default(present)
     do k = kb - kh, ke + kh
       do i = ib - ih, ie + ih
         u0_d(i, je + 1, k) = u0_d(i, je + 1, k) - (u0_d(i, je + 1, k) - u0_d(i, je, k))*l_dyi*l_rk3coef*l_vouttot
         um_d(i, je + 1, k) = um_d(i, je + 1, k) - (um_d(i, je + 1, k) - um_d(i, je, k))*l_dyi*l_rk3coef*l_vouttot
         w0_d(i, je + 1, k) = w0_d(i, je + 1, k) - (w0_d(i, je + 1, k) - w0_d(i, je, k))*l_dyi*l_rk3coef*l_vouttot
         wm_d(i, je + 1, k) = wm_d(i, je + 1, k) - (wm_d(i, je + 1, k) - wm_d(i, je, k))*l_dyi*l_rk3coef*l_vouttot
       end do
     end do
     !$acc end parallel loop

     if (loneeqn) then
       !$acc parallel loop collapse(2) default(present)
       do k = kb - kh, ke + kh
         do i = ib - ih, ie + ih
           e120_d(i, je + 1, k) = e120_d(i, je + 1, k) - (e120_d(i, je + 1, k) - e120_d(i, je, k))*l_dyi*l_rk3coef*l_vouttot
           e12m_d(i, je + 1, k) = e12m_d(i, je + 1, k) - (e12m_d(i, je + 1, k) - e12m_d(i, je, k))*l_dyi*l_rk3coef*l_vouttot
         end do
       end do
       !$acc end parallel loop
     end if

   end subroutine ymo_convective_device
#endif


   subroutine yTo_convective

     use modglobal, only : je, dyi, rk3coef
     use modfields, only : thl0, thlm, vouttot

     thl0(:, je + 1, :) = thl0(:, je + 1, :) - (thl0(:, je + 1, :) - thl0(:, je, :))*dyi*rk3coef*vouttot
     thlm(:, je + 1, :) = thlm(:, je + 1, :) - (thlm(:, je + 1, :) - thlm(:, je, :))*dyi*rk3coef*vouttot

   end subroutine yTo_convective
#if defined(_GPU)
   subroutine yTo_convective_device
     use modglobal, only : ib, ie, ih, je, kb, ke, kh, dyi, rk3coef
     use modfields, only : vouttot
     use modcuda,   only : thl0_d, thlm_d
     implicit none
     integer i, k
     real    :: l_dyi, l_rk3coef, l_vouttot

     l_dyi = dyi
     l_rk3coef = rk3coef
     l_vouttot = vouttot

     !$acc parallel loop collapse(2) default(present)
     do k = kb - kh, ke + kh
       do i = ib - ih, ie + ih
         thl0_d(i, je + 1, k) = thl0_d(i, je + 1, k) - (thl0_d(i, je + 1, k) - thl0_d(i, je, k))*l_dyi*l_rk3coef*l_vouttot
         thlm_d(i, je + 1, k) = thlm_d(i, je + 1, k) - (thlm_d(i, je + 1, k) - thlm_d(i, je, k))*l_dyi*l_rk3coef*l_vouttot
       end do
     end do
     !$acc end parallel loop

   end subroutine yTo_convective_device
#endif


   subroutine yqo_convective

     use modglobal, only : je, dyi, rk3coef
     use modfields, only : qt0, qtm, vouttot

     qt0(:, je + 1, :) = qt0(:, je + 1, :) - (qt0(:, je + 1, :) - qt0(:, je, :))*dyi*rk3coef*vouttot
     qtm(:, je + 1, :) = qtm(:, je + 1, :) - (qtm(:, je + 1, :) - qtm(:, je, :))*dyi*rk3coef*vouttot

   end subroutine yqo_convective
#if defined(_GPU)
   subroutine yqo_convective_device
     use modglobal, only : ib, ie, ih, je, kb, ke, kh, dyi, rk3coef
     use modfields, only : vouttot
     use modcuda,   only : qt0_d, qtm_d
     implicit none
     integer i, k
     real    :: l_dyi, l_rk3coef, l_vouttot

     l_dyi = dyi
     l_rk3coef = rk3coef
     l_vouttot = vouttot

     !$acc parallel loop collapse(2) default(present)
     do k = kb - kh, ke + kh
       do i = ib - ih, ie + ih
         qt0_d(i, je + 1, k) = qt0_d(i, je + 1, k) - (qt0_d(i, je + 1, k) - qt0_d(i, je, k))*l_dyi*l_rk3coef*l_vouttot
         qtm_d(i, je + 1, k) = qtm_d(i, je + 1, k) - (qtm_d(i, je + 1, k) - qtm_d(i, je, k))*l_dyi*l_rk3coef*l_vouttot
       end do
     end do
     !$acc end parallel loop

   end subroutine yqo_convective_device
#endif


   subroutine yso_convective

     use modglobal, only : je, rk3coef, dyi, nsv
     use modfields, only :sv0, svm, vouttot
     integer n

     do n = 1, nsv
       sv0(:, je + 1, :, n) = sv0(:, je + 1, :, n) - (sv0(:, je + 1, :, n) - sv0(:, je, :, n))*dyi*rk3coef*vouttot
       svm(:, je + 1, :, n) = svm(:, je + 1, :, n) - (svm(:, je + 1, :, n) - svm(:, je, :, n))*dyi*rk3coef*vouttot
     end do

   end subroutine yso_convective
#if defined(_GPU)
   subroutine yso_convective_device
     use modglobal, only : ib, ie, ihc, je, kb, ke, khc, rk3coef, dyi, nsv
     use modfields, only : vouttot
     use modcuda,   only : sv0_d, svm_d
     implicit none
     integer i, k, n
     real    :: l_dyi, l_rk3coef, l_vouttot

     l_dyi = dyi
     l_rk3coef = rk3coef
     l_vouttot = vouttot

     !$acc parallel loop collapse(3) default(present)
     do n = 1, nsv
       do k = kb - khc, ke + khc
         do i = ib - ihc, ie + ihc
           sv0_d(i, je + 1, k, n) = sv0_d(i, je + 1, k, n) - (sv0_d(i, je + 1, k, n) - sv0_d(i, je, k, n))*l_dyi*l_rk3coef*l_vouttot
           svm_d(i, je + 1, k, n) = svm_d(i, je + 1, k, n) - (svm_d(i, je + 1, k, n) - svm_d(i, je, k, n))*l_dyi*l_rk3coef*l_vouttot
         end do
       end do
     end do
     !$acc end parallel loop

   end subroutine yso_convective_device
#endif


  !  subroutine yso_Neumann

  !    use modglobal, only : je, jhc, nsv
  !    use modfields, only : sv0, svm
  !    integer n, m

  !    do n = 1, nsv
  !      do m = 1, jhc
  !        sv0(:, je + m, :, n) = sv0(:, je, :, n)
  !        svm(:, je + m, :, n) = svm(:, je, :, n)
  !      end do
  !    end do

  !  end subroutine yso_Neumann


#if defined(_GPU)
   attributes(global) subroutine bcpup_pwp_BCtopm_xslip_cuda
     use modcuda, only: pwp_d, ie_d, je_d, kb_d, ke_d, kh_d, tidandstride
     implicit none
     integer :: i, j, tidx, tidy, tidz, stridex, stridey, stridez
     call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
     if (tidz==kb_d) then
       do j = tidy, je_d, stridey
         do i = tidx, ie_d, stridex
           pwp_d(i, j, kb_d) = 0.
         end do
       end do
     end if
     if (tidz==stridez) then
       do j = tidy, je_d, stridey
         do i = tidx, ie_d, stridex
           pwp_d(i, j, ke_d + kh_d) = 0.
         end do
       end do
     end if
   end subroutine bcpup_pwp_BCtopm_xslip_cuda

   attributes(global) subroutine bcpup_pwp_BCtopm_pressure_cuda(rk3coefi, pres0ij_ke)
     use modcuda, only: pwp_d, wm_d, wp_d, ie_d, je_d, kb_d, ke_d, dzhi_d, tidandstride
     implicit none
     real, value, intent(in) :: rk3coefi, pres0ij_ke
     integer :: i, j, tidx, tidy, tidz, stridex, stridey, stridez
     call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
     if (tidz==kb_d) then
       do j = tidy, je_d, stridey
         do i = tidx, ie_d, stridex
           pwp_d(i, j, kb_d)  = 0.
         end do
       end do
     end if
     if (tidz==stridez) then
       do j = tidy, je_d, stridey
         do i = tidx, ie_d, stridex
           pwp_d(i, j, ke_d + 1) = wm_d(i, j, ke_d+1) * rk3coefi + 2 * pres0ij_ke*dzhi_d(ke_d+1)
           wp_d(i, j, ke_d + 1) = pwp_d(i, j, ke_d+1) - wm_d(i,j,ke_d+1) * rk3coefi
         end do
       end do
     end if
   end subroutine bcpup_pwp_BCtopm_pressure_cuda

   attributes(global) subroutine bcpup_pup_BCxm_periodic_cuda
     use modcuda, only: pup_d, ib_d, ie_d, je_d, ke_d, tidandstride
     implicit none
     integer :: j, k, tidx, tidy, tidz, stridex, stridey, stridez
     call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
     if (tidx == ib_d) then
       do k = tidz, ke_d, stridez
         do j = tidy, je_d, stridey
           pup_d(ie_d+1, j, k) = pup_d(ib_d, j, k)
         end do
       end do
     end if
   end subroutine bcpup_pup_BCxm_periodic_cuda

   attributes(global) subroutine bcpup_pup_BCxm_profile_cuda(ibrank, ierank, rk3coefi, dxi, uouttot)
     use modcuda, only: pup_d, up_d, um_d, u0_d, uprof_d, ib_d, ie_d, je_d, kb_d, ke_d, tidandstride
     implicit none
     logical, value, intent(in) :: ibrank, ierank
     real   , value, intent(in) :: rk3coefi, dxi, uouttot
     integer :: j, k, tidx, tidy, tidz, stridex, stridey, stridez
     call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
     if (ibrank) then
       if (tidx == ib_d) then
         do k = tidz, ke_d, stridez
           do j = tidy-1, je_d+1, stridey
             pup_d(ib_d, j, k) = uprof_d(k) * rk3coefi
             up_d(ib_d, j, k) = 0.
           end do
         end do
       end if
     end if
     if (ierank) then
       if (tidx == stridex) then
         do k = tidz, ke_d, stridez
           do j = tidy-1, je_d+1, stridey
             if (k==kb_d) then
               pup_d(ie_d+1, j, k) = pup_d(ie_d, j, k)
             else
               pup_d(ie_d+1, j, k) = um_d(ie_d+1, j, k) * rk3coefi - (u0_d(ie_d+1, j, k) - u0_d(ie_d,j,k)) * dxi * uouttot
             end if
             up_d(ie_d+1, j, k)  = pup_d(ie_d+1, j, k) - um_d(ie_d+1, j, k) * rk3coefi
           end do
         end do
       end if
     end if
   end subroutine bcpup_pup_BCxm_profile_cuda 
   
   attributes(global) subroutine bcpup_pup_BCxm_driver_cuda(ibrank, ierank, rk3coefi, dxi, uouttot)
     use modcuda, only: pup_d, up_d, um_d, u0_d, u0driver_d, ib_d, ie_d, je_d, ke_d, tidandstride
     implicit none
     logical, value, intent(in) :: ibrank, ierank
     real   , value, intent(in) :: rk3coefi, dxi, uouttot
     integer :: j, k, tidx, tidy, tidz, stridex, stridey, stridez
     call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
     if (ibrank) then
       if (tidx == ib_d) then
         do k = tidz, ke_d, stridez
           do j = tidy-1, je_d+1, stridey
             pup_d(ib_d, j, k) = u0driver_d(j, k) * rk3coefi
             up_d(ib_d, j, k) = 0.
           end do
         end do
       end if
     end if
     if (ierank) then
       if (tidx == stridex) then
         do k = tidz, ke_d, stridez
           do j = tidy-1, je_d+1, stridey
             pup_d(ie_d+1, j, k) = um_d(ie_d+1, j, k) * rk3coefi - ( u0_d(ie_d+1, j, k) - u0_d(ie_d, j, k) ) * dxi * uouttot
             up_d(ie_d+1, j, k)  = pup_d(ie_d+1, j, k) - um_d(ie_d+1, j, k)*rk3coefi
           end do
         end do
       end if
     end if
   end subroutine bcpup_pup_BCxm_driver_cuda
 
   attributes(global) subroutine bcpup_pvp_BCym_periodic_cuda
     use modcuda, only: pvp_d, ie_d, jb_d, je_d, ke_d, tidandstride
     implicit none
     integer :: i, k, tidx, tidy, tidz, stridex, stridey, stridez
     call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
     if (tidy == jb_d) then
       do k = tidz, ke_d, stridez
         do i = tidx, ie_d, stridex
           pvp_d(i, je_d+1, k) = pvp_d(i, jb_d, k)
         end do
       end do
     end if
   end subroutine bcpup_pvp_BCym_periodic_cuda
   
   attributes(global) subroutine bcpup_pvp_BCym_profile_cuda(jbrank, jerank, rk3coefi, dyi, vouttot)
     use modcuda, only: pvp_d, vp_d, vm_d, v0_d, vprof_d, ie_d, jb_d, je_d, kb_d, ke_d, tidandstride
     implicit none
     logical, value, intent(in) :: jbrank, jerank
     real   , value, intent(in) :: rk3coefi, dyi, vouttot
     integer :: i, k, tidx, tidy, tidz, stridex, stridey, stridez
     call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
     if (jbrank) then
       if (tidy == jb_d) then
         do k = tidz, ke_d, stridez
           do i = tidx-1, ie_d+1, stridex
             pvp_d(i, jb_d, k) = vprof_d(k) * rk3coefi
             vp_d(i, jb_d, k) = 0.
           end do
         end do
       end if
     end if
     if (jerank) then
       if (tidy == stridey) then
         do k = tidz, ke_d, stridez
           do i = tidx-1, ie_d+1, stridex
             if (k==kb_d) then
               pvp_d(i, je_d+1, k) = pvp_d(i, je_d, k)
             else
               pvp_d(i, je_d+1, k) = vm_d(i, je_d+1, k) * rk3coefi - ( v0_d(i, je_d+1, k) - v0_d(i, je_d, k) ) * dyi * vouttot
             end if
             vp_d(i, je_d+1, k)  = pvp_d(i, je_d+1, k) - vm_d(i, je_d+1, k) * rk3coefi
           end do
         end do
       end if
     end if
   end subroutine bcpup_pvp_BCym_profile_cuda
#endif

   !>set boundary conditions pup,pvp,pwp in subroutine fillps in modpois.f90
   subroutine bcpup
     use modglobal,    only : ib, ie, ih, jb, je, jh, kb, ke, kh, dxi, dyi, dzhi, rk3coefi, &
                              ibrank, ierank, jbrank, jerank, BCxm, BCym, BCtopm, &
                              BCtopm_freeslip, BCtopm_noslip, BCtopm_pressure, &
                              BCxm_periodic, BCxm_profile, BCxm_driver, &
                              BCym_periodic, BCym_profile
     use modfields,    only : uouttot, vouttot, pres0, pres0, IIc, IIcs
     use modmpi,       only : avexy_ibm
     use m_halo,       only : halo_exchange
#if defined(_GPU)
     use cudafor
     use modcuda,      only : pup_d, pvp_d, pwp_d, griddim, blockdim, checkCUDA
#else
     use modfields,    only : pup, pvp, pwp, up, vp, wp, um, vm, wm, u0, v0, uprof, vprof
     use modinletdata, only : u0driver
#endif
     implicit none
     real, dimension(kb:ke+kh) :: pres0ij
     integer :: i, j, k

#if defined(_GPU)
     call halo_exchange(pup_d, 3, opt_levels=(/ih,jh,0/))
     call halo_exchange(pvp_d, 3, opt_levels=(/ih,jh,0/))
     call halo_exchange(pwp_d, 3, opt_levels=(/ih,jh,0/))
#else
     call halo_exchange(pup, 3, opt_levels=(/ih,jh,0/))
     call halo_exchange(pvp, 3, opt_levels=(/ih,jh,0/))
     call halo_exchange(pwp, 3, opt_levels=(/ih,jh,0/))
#endif

     select case(BCtopm)
       case(BCtopm_freeslip, BCtopm_noslip)
#if defined(_GPU)
         call bcpup_pwp_BCtopm_xslip_cuda<<<griddim,blockdim>>>
         call checkCUDA( cudaGetLastError(), 'bcpup_pwp_BCtopm_xslip_cuda' )
#else
         do j = jb, je
           do i = ib, ie
             pwp(i, j, kb) = 0.
             pwp(i, j, ke + kh) = 0.
           end do
         end do
#endif

       case(BCtopm_pressure)
         call avexy_ibm(pres0ij(kb:ke+kh),pres0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)

#if defined(_GPU)
         call bcpup_pwp_BCtopm_pressure_cuda<<<griddim,blockdim>>>(rk3coefi, pres0ij(ke))
         call checkCUDA( cudaGetLastError(), 'bcpup_pwp_BCtopm_pressure_cuda' )
#else
         do j = jb, je
           do i = ib, ie
             pwp(i, j, kb)  = 0.
             !pwp(i, j, ke + 1) = wm(i, j, ke+1) * rk3coefi - (-pres0ij(ke) - pres0(i,j,ke)) * dzhi(ke+1) ! Doesn't work
             pwp(i, j, ke + 1) = wm(i, j, ke+1) * rk3coefi + 2 * pres0ij(ke)*dzhi(ke+1)
             wp(i, j, ke + 1) = pwp(i, j, ke+1) - wm(i,j,ke+1) * rk3coefi
           end do
         end do
#endif
     end select !BCtopm

     select case(BCxm)
       case(BCxm_periodic)
         if (ibrank .and. ierank) then ! not parallelised in x
#if defined(_GPU)
           call bcpup_pup_BCxm_periodic_cuda<<<griddim,blockdim>>>
           call checkCUDA( cudaGetLastError(), 'bcpup_pup_BCxm_periodic_cuda' )
#else
           do k = kb, ke
              do j = jb, je
                 pup(ie+1, j, k) = pup(ib, j, k) ! cyclic
              end do
           end do
#endif
         end if

       case(BCxm_profile)
#if defined(_GPU)
         call bcpup_pup_BCxm_profile_cuda<<<griddim,blockdim>>>(ibrank, ierank, rk3coefi, dxi, uouttot)
         call checkCUDA( cudaGetLastError(), 'bcpup_pup_BCxm_profile_cuda' )
#else
         if (ibrank) then
           do k=kb,ke
             do j=jb-1,je+1
               pup(ib, j, k) = uprof(k) * rk3coefi
               up(ib, j, k) = 0.
             end do
           end do
         end if

         if (ierank) then
           do k = kb+1, ke
             do j = jb-1, je+1
               ! convective
               pup(ie+1, j, k) = um(ie+1, j, k) * rk3coefi - (u0(ie+1, j, k) - u0(ie,j,k)) * dxi * uouttot !u0(ie,j,k) ! du/dt +u*du/dx=0 -> pup(i)=um(i)/rk3coef -um(i)*(um(i)-um(i-1))/dxf(i-1)
               ! Neumann
               !pup(ie+1,j,k) = pup(ie,j,k)
               up(ie+1, j, k) = pup(ie+1, j, k) - um(ie+1,j,k)*rk3coefi
             end do
           end do
           ! Neumann at bottom - performs better with no slip
           pup(ie+1, :, kb) = pup(ie, :, kb)
           up(ie+1, :, kb) = pup(ie+1,: , kb) - um(ie+1, :, kb) * rk3coefi
         end if
#endif

       case(BCxm_driver)
#if defined(_GPU)
         call bcpup_pup_BCxm_driver_cuda<<<griddim,blockdim>>>(ibrank, ierank, rk3coefi, dxi, uouttot)
         call checkCUDA( cudaGetLastError(), 'bcpup_pup_BCxm_driver_cuda' )
#else
         if (ibrank) then
           do k = kb, ke
             do j = jb - 1, je + 1
               pup(ib, j, k) = u0driver(j, k) * rk3coefi
               up(ib, j, k) = 0. ! u(ib) only evolves according to pressure correction
             end do
           end do
         end if

         if (ierank) then
           do k = kb, ke
             do j = jb-1, je+1
               pup(ie+1, j, k) = um(ie+1, j, k) * rk3coefi - (u0(ie+1, j, k) - u0(ie, j, k)) * dxi * uouttot    ! du/dt +u*du/dx=0 -> pup(i)=um(i)/rk3coef -um(i)*(um(i)-um(i-1))/dxf(i-1)
               ! !Neumann
               !pup(ie+1,j,k) = pup(ie,j,k)
               up(ie+1, j, k) = pup(ie+1, j, k) - um(ie+1, j, k) * rk3coefi
             end do
           end do
           ! Neumann at bottom - performs better with no slip
           ! pup(ie+1, :, kb) = pup(ie, :, kb)
           ! up(ie+1, :, kb) = pup(ie+1, :, kb) - um(ie+1, :, kb) * rk3coefi
         end if
#endif
     end select ! BCxm
    
     select case(BCym)
       case(BCym_periodic)
         if (jbrank .and. jerank) then ! not parallelised in y
#if defined(_GPU)
           call bcpup_pvp_BCym_periodic_cuda<<<griddim,blockdim>>>
           call checkCUDA( cudaGetLastError(), 'bcpup_pvp_BCym_periodic_cuda' )
#else
           do k = kb, ke
             do i = ib, ie
               pvp(i, je+1, k) = pvp(i, jb, k) ! cyclic
             end do
           end do
#endif
         end if

       case(BCym_profile)
#if defined(_GPU)
         call bcpup_pvp_BCym_profile_cuda<<<griddim,blockdim>>>(jbrank, jerank, rk3coefi, dyi, vouttot)
         call checkCUDA( cudaGetLastError(), 'bcpup_pvp_BCym_profile_cuda' )
#else
         if (jbrank) then
           do k = kb, ke
             do i = ib-1, ie+1
               pvp(i,jb,k) = vprof(k)*rk3coefi
               vp(i,jb,k) = 0.
             end do
           end do
         end if

         if (jerank) then
           do k = kb, ke
             do i = ib-1, ie+1
               ! change to vouttot
               pvp(i, je+1, k) = vm(i, je+1, k) * rk3coefi - (v0(i, je+1, k) - v0(i, je, k)) * dyi * vouttot
               vp(i, je+1, k) = pvp(i, je+1, k) - vm(i,je+1,k)*rk3coefi
             end do
           end do
           pvp(:, je+1, kb) = pvp(:, je, kb)
           vp(:, je+1, kb) = pvp(:, je+1, kb) - vm(:, je+1, kb)*rk3coefi
         end if
#endif
     end select ! BCym

   end subroutine bcpup


   !>set pressure boundary conditions
   subroutine bcp
     use modglobal, only : ib, ie, jb, je, kb, ke, &
                           ibrank, ierank, jbrank, jerank, BCxm, BCym, BCxm_periodic, BCym_periodic
     use m_halo,    only : halo_exchange
#if defined(_GPU)
     use cudafor
     use modcuda,   only : p_d, pres0_d
#else
     use modfields, only : p, pres0
#endif
     implicit none
     integer :: i, j, k

#if defined(_GPU)
     call halo_exchange(p_d, 3)
     call halo_exchange(pres0_d, 3)
#else
     call halo_exchange(p, 3)
     call halo_exchange(pres0, 3)
#endif

     if (BCxm .eq. BCxm_periodic) then
     
       if (ibrank .and. ierank) then
         !$acc kernels default(present)
         do j = jb, je
           do k = kb, ke
#if defined(_GPU)
             p_d(ib-1, j, k) = p_d(ie, j, k)
             p_d(ie+1, j, k) = p_d(ib, j, k)
#else
             p(ib-1, j, k) = p(ie, j, k)
             p(ie+1, j, k) = p(ib, j, k)
             !pres0(ib - 1, j, k) = pres0(ie, j, k)
             !pres0(ie + 1, j, k) = pres0(ib, j, k)
#endif
           end do
         end do
         !$acc end kernels
       end if

     else
     
       if (ibrank) then
         !$acc kernels default(present)
         do k = kb, ke
           do j = jb-1, je+1
#if defined(_GPU)
             p_d(ib-1, j, k) = p_d(ib, j, k)
             pres0_d(ib-1, j, k) = pres0_d(ib, j, k)
#else
             p(ib-1, j, k) = p(ib, j, k)
             pres0(ib-1, j, k) = pres0(ib, j, k)
#endif
           end do
         end do
         !$acc end kernels
       end if

       if (ierank) then
         !$acc kernels default(present)
         do k = kb, ke
           do j = jb-1, je+1
#if defined(_GPU)
             p_d(ie+1, j, k) = p_d(ie, j, k)
             pres0_d(ie+1, j, k) = pres0_d(ie, j, k)
#else
             p(ie+1, j, k) = p(ie, j, k)
             pres0(ie+1, j, k) = pres0(ie, j, k)
#endif
           end do
         end do
         !$acc end kernels
       end if

     end if ! BCxm

     if (BCym .eq. BCym_periodic) then

       if (jbrank .and. jerank) then
         !$acc kernels default(present)
         do i = ib, ie
           do k = kb, ke
#if defined(_GPU)
             p_d(i, jb-1, k) = p_d(i, je, k)
             p_d(i, je+1, k) = p_d(i, jb, k)
#else
             p(i, jb-1, k) = p(i, je, k)
             p(i, je+1, k) = p(i, jb, k)
             !pres0(ib - 1, j, k) = pres0(ie, j, k)
             !pres0(ie + 1, j, k) = pres0(ib, j, k)
#endif
           end do
         end do
         !$acc end kernels
       end if

     else

       if (jbrank) then
         !$acc kernels default(present)
         do k = kb, ke
           do i = ib-1, ie+1
#if defined(_GPU)
             p_d(i,jb-1,k) = p_d(i,jb,k)
             pres0_d(i,jb-1,k) = pres0_d(i,jb,k)
#else
             p(i,jb-1,k) = p(i,jb,k)
             pres0(i,jb-1,k) = pres0(i,jb,k)
#endif
           end do
         end do
         !$acc end kernels
       end if

       if (jerank) then
         !$acc kernels default(present)
         do k = kb, ke
           do i = ib-1, ie+1
#if defined(_GPU)
             p_d(i, je+1, k) = p_d(i,je,k)
             pres0_d(i, je+1, k) = pres0_d(i,je,k)
#else
             p(i, je+1, k) = p(i,je,k)
             pres0(i, je+1, k) = pres0(i,je,k)
#endif
           end do
         end do
         !$acc end kernels
       end if

     end if !BCym

   end subroutine bcp

   !>
   !! grwdamp damps gravity waves in the upper part of the domain.
   !>
   !! The lower limit of the damping region is set by ksp
   !! Horizontal fluctuations at the top of the domain (for instance gravity waves)
   !! are damped out by a sponge layer through an additional forcing/source term.
   !! \latexonly
   !! \begin{eqnarray}
   !! \force{i}{sp}(z) &=& -\frac{1}{t^{\mr{sp}}}\left(\xav{\fav{u_i}}-\fav{u_i}\right), \\\\
   !!  \source{\varphi}{sp}(z) &=& -\frac{1}{t^{\mr{sp}}}\left(\xav{\varphi}-\varphi\right),
   !! \end{eqnarray}
   !! with $t^{\mr{sp}}$ a relaxation time scale that goes from
   !! $t^{\mr{sp}}_0=1/(2.75\times10^{-3})\mr{s}\approx 6$min at the top of the domain
   !! to infinity at the bottom of the sponge layer.
   !! \endlatexonly
   subroutine grwdamp
      use modglobal, only: ke, ltempeq, lmoist, lcoriol, igrw_damp, geodamptime
#if defined(_GPU)
      use modcuda,   only: up_d, vp_d, wp_d, thlp_d, qtp_d, u0_d, v0_d, w0_d, thl0_d, qt0_d, ug_d, vg_d, u0av_d, v0av_d, thl0av_d, qt0av_d, tsc_d
#else
      use modfields, only: up, vp, wp, thlp, qtp, u0, v0, w0, thl0, qt0, ug, vg, u0av, v0av, thl0av, qt0av, tsc
#endif
      implicit none
      integer :: k
      select case (igrw_damp)
      case (0) !do nothing
      case (1)
         !$acc kernels default(present)
         do k = ksp, ke
#if defined(_GPU)
            up_d(:, :, k) = up_d(:, :, k) - (u0_d(:, :, k) - u0av_d(k))*tsc_d(k)
            vp_d(:, :, k) = vp_d(:, :, k) - (v0_d(:, :, k) - v0av_d(k))*tsc_d(k)
            wp_d(:, :, k) = wp_d(:, :, k) - w0_d(:, :, k)*tsc_d(k)
#else
            up(:, :, k) = up(:, :, k) - (u0(:, :, k) - u0av(k))*tsc(k)
            vp(:, :, k) = vp(:, :, k) - (v0(:, :, k) - v0av(k))*tsc(k)
            wp(:, :, k) = wp(:, :, k) - w0(:, :, k)*tsc(k)
#endif
         end do
         !$acc end kernels

         if (ltempeq) then
           !$acc kernels default(present)
           do k = ksp, ke
#if defined(_GPU)
             thlp_d(:, :, k) = thlp_d(:, :, k) - (thl0_d(:, :, k) - thl0av_d(k))*tsc_d(k)
#else
             thlp(:, :, k) = thlp(:, :, k) - (thl0(:, :, k) - thl0av(k))*tsc(k)
#endif
           end do
           !$acc end kernels
         end if

         if (lmoist) then
           !$acc kernels default(present)
           do k = ksp, ke
#if defined(_GPU)
              qtp_d(:, :, k) = qtp_d(:, :, k) - (qt0_d(:, :, k) - qt0av_d(k))*tsc_d(k)
#else
              qtp(:, :, k) = qtp(:, :, k) - (qt0(:, :, k) - qt0av(k))*tsc(k)
#endif
           end do
           !$acc end kernels
         end if

         if (lcoriol) then
            !$acc kernels default(present)
            do k = ksp, ke
#if defined(_GPU)
               up_d(:, :, k) = up_d(:, :, k) - (u0_d(:, :, k) - ug_d(k))*((1./(geodamptime*rnu0))*tsc_d(k))
               vp_d(:, :, k) = vp_d(:, :, k) - (v0_d(:, :, k) - vg_d(k))*((1./(geodamptime*rnu0))*tsc_d(k))
#else
               up(:, :, k) = up(:, :, k) - (u0(:, :, k) - ug(k))*((1./(geodamptime*rnu0))*tsc(k))
               vp(:, :, k) = vp(:, :, k) - (v0(:, :, k) - vg(k))*((1./(geodamptime*rnu0))*tsc(k))
#endif
            end do
            !$acc end kernels
         end if
      case (2)
         !$acc kernels default(present)
         do k = ksp, ke
#if defined(_GPU)
            up_d(:, :, k) = up_d(:, :, k) - (u0_d(:, :, k) - ug_d(k))*tsc_d(k)
            vp_d(:, :, k) = vp_d(:, :, k) - (v0_d(:, :, k) - vg_d(k))*tsc_d(k)
            wp_d(:, :, k) = wp_d(:, :, k) - w0_d(:, :, k)*tsc_d(k)
#else
            up(:, :, k) = up(:, :, k) - (u0(:, :, k) - ug(k))*tsc(k)
            vp(:, :, k) = vp(:, :, k) - (v0(:, :, k) - vg(k))*tsc(k)
            wp(:, :, k) = wp(:, :, k) - w0(:, :, k)*tsc(k)
#endif
         end do
         !$acc end kernels

         if (ltempeq) then
           !$acc kernels default(present)
           do k = ksp, ke
#if defined(_GPU)
             thlp_d(:, :, k) = thlp_d(:, :, k) - (thl0_d(:, :, k) - thl0av_d(k))*tsc_d(k)
#else
             thlp(:, :, k) = thlp(:, :, k) - (thl0(:, :, k) - thl0av(k))*tsc(k)
#endif
           end do
           !$acc end kernels
         end if

         if (lmoist) then
           !$acc kernels default(present)
           do k = ksp, ke
#if defined(_GPU)
              qtp_d(:, :, k) = qtp_d(:, :, k) - (qt0_d(:, :, k) - qt0av_d(k))*tsc_d(k)
#else
              qtp(:, :, k) = qtp(:, :, k) - (qt0(:, :, k) - qt0av(k))*tsc(k)
#endif
           end do
           !$acc end kernels
         end if
      case (3)
         !$acc kernels default(present)
         do k = ksp, ke
#if defined(_GPU)
            up_d(:, :, k) = up_d(:, :, k) - (u0_d(:, :, k) - u0av_d(k))*tsc_d(k)
            vp_d(:, :, k) = vp_d(:, :, k) - (v0_d(:, :, k) - v0av_d(k))*tsc_d(k)
            wp_d(:, :, k) = wp_d(:, :, k) - w0_d(:, :, k)*tsc_d(k)
#else
            up(:, :, k) = up(:, :, k) - (u0(:, :, k) - u0av(k))*tsc(k)
            vp(:, :, k) = vp(:, :, k) - (v0(:, :, k) - v0av(k))*tsc(k)
            wp(:, :, k) = wp(:, :, k) - w0(:, :, k)*tsc(k)
#endif
         end do
         !$acc end kernels

         if (ltempeq) then
           !$acc kernels default(present)
           do k = ksp, ke
#if defined(_GPU)
             thlp_d(:, :, k) = thlp_d(:, :, k) - (thl0_d(:, :, k) - thl0av_d(k))*tsc_d(k)
#else
             thlp(:, :, k) = thlp(:, :, k) - (thl0(:, :, k) - thl0av(k))*tsc(k)
#endif
           end do
           !$acc end kernels
         end if

         if (lmoist) then
           !$acc kernels default(present)
           do k = ksp, ke
#if defined(_GPU)
              qtp_d(:, :, k) = qtp_d(:, :, k) - (qt0_d(:, :, k) - qt0av_d(k))*tsc_d(k)
#else
              qtp(:, :, k) = qtp(:, :, k) - (qt0(:, :, k) - qt0av(k))*tsc(k)
#endif
           end do
           !$acc end kernels
         end if
      case default
         write(0, *) "ERROR: no gravity wave damping option selected"
         stop 1
      end select
   end subroutine grwdamp


   subroutine fluxtop(field, ek, flux)
      use modglobal, only:ib, ie, ih, jb, je, jh, kb, ke, kh, dzf, dzh, dzhi, eps1

      real, intent(inout) :: field(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh)
      real, intent(in)    ::    ek(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh)
      real, intent(in)    :: flux
      !
      if (abs(flux) .le. eps1) then !it's zero-flux, we don't need to do the calculation
         field(:, :, ke + 1) = field(:, :, ke)
      else
         field(:, :, ke + 1) = field(:, :, ke) + dzh(ke + 1)*flux/(dzhi(ke + 1)*(0.5*(dzf(ke)*ek(:, :, ke + 1) + dzf(ke + 1)*ek(:, :, ke))))
      end if
      !
   end subroutine fluxtop

   subroutine valuetop(field, val)
      use modglobal, only:ib, ie, ih, jb, je, jh, kb, ke, kh
      real, intent(inout) :: field(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh)
      real, intent(in)    :: val

      ! (field(i, j, kp)*dzf(k) + field(i, j, k)*dzf(kp))*dzhi(kp)*0.5 = val
      !field(:,:,ke+1) = (2.*val*dzh(ke+1) - field(:,:,ke)*dzf(ke+1)) * dzfi(ke)
      field(:, :, ke + 1) = 2*val - field(:, :, ke)
      !if (myid == 0) write(*,*) (field(40, 1, ke+1)*dzf(ke) + field(40, 1, ke)*dzf(ke+1))*dzhi(ke+1)*0.5

   end subroutine valuetop

   subroutine fluxtopscal(flux)
      use modglobal, only:ib, ie, ih, jb, je, jh, ke, dzf, dzh, dzhi, nsv, khc
      use modfields, only:sv0, svm
      use modsubgriddata, only:ekh

      real, intent(in)    :: flux(:)
      integer :: m, n
      !
      !all the ghost cells have the same value?
      do m = 1, khc
      do n = 1, nsv
  sv0(ib-ih:ie+ih,jb-jh:je+jh,ke+m,n) = sv0(ib-ih:ie+ih,jb-jh:je+jh,ke,n) + dzh(ke+1) * flux(n) / ( dzhi(ke+1) * (0.5*(dzf(ke)*ekh(ib-ih:ie+ih,jb-jh:je+jh,ke+1)+dzf(ke+1)*ekh(ib-ih:ie+ih,jb-jh:je+jh,ke))))
  svm(ib-ih:ie+ih,jb-jh:je+jh,ke+m,n) = svm(ib-ih:ie+ih,jb-jh:je+jh,ke,n) + dzh(ke+1) * flux(n) / ( dzhi(ke+1) * (0.5*(dzf(ke)*ekh(ib-ih:ie+ih,jb-jh:je+jh,ke+1)+dzf(ke+1)*ekh(ib-ih:ie+ih,jb-jh:je+jh,ke))))
      end do
      end do
      !
   end subroutine fluxtopscal

   subroutine valuetopscal(val)
      use modglobal, only:ke, nsv, khc
      use modfields, only:sv0, svm
      real, intent(in)    :: val(1:nsv)
      integer :: m, n
      !
      ! all the ghost cells have the same vlaue?
      do m = 1, khc
      do n = 1, nsv
         sv0(: , : , ke + m, n) = 2*val(n) - sv0(: , : , ke, n)
         svm(: , : , ke + m, n) = 2*val(n) - svm(: , : , ke, n)
      end do
      end do
      !
   end subroutine valuetopscal

#if defined(_GPU)
   ! Device twins of the four top-boundary helpers above, plus the plane
   ! assignments boundary makes around them.
   !
   ! They take no arguments, unlike fluxtop and valuetop. An assumed-shape or
   ! explicit-shape device dummy would put a descriptor between the kernel and
   ! the array for no gain: each of these is wanted for one fixed set of
   ! fields, and naming them directly is what lets the compiler see the
   ! addresses it saw before this existed.
   !
   ! fluxtop_uv_device, fluxtop_thl_device, fluxtop_qt_device and
   ! fluxtop_sv_device are called from two places - boundary_device here and
   ! reassure_fluxtop_boundary, which applies them again once closurebc has
   ! moved the eddy diffusivity. They were duplicated between the two before,
   ! with the division arranged differently on each side.

   !> Zero-flux top boundary for the horizontal velocities.
   !!
   !! fluxtop(field, ek, 0.) is a copy of the top fluid level into the ghost
   !! level - the diffusivity never enters - so the four fields are one kernel.
   subroutine fluxtop_uv_device
      use modglobal, only : ib, ie, ih, jb, je, jh, ke
      use modcuda,   only : u0_d, um_d, v0_d, vm_d
      implicit none
      integer :: i, j

      !$acc parallel loop collapse(2) default(present)
      do j = jb - jh, je + jh
         do i = ib - ih, ie + ih
            um_d(i, j, ke + 1) = um_d(i, j, ke)
            u0_d(i, j, ke + 1) = u0_d(i, j, ke)
            vm_d(i, j, ke + 1) = vm_d(i, j, ke)
            v0_d(i, j, ke + 1) = v0_d(i, j, ke)
         end do
      end do
      !$acc end parallel loop
   end subroutine fluxtop_uv_device

   !> Fixed-velocity top boundary for the horizontal velocities.
   subroutine valuetop_uv_device
      use modglobal, only : ib, ie, ih, jb, je, jh, ke, Uinf, Vinf
      use modcuda,   only : u0_d, um_d, v0_d, vm_d
      implicit none
      integer :: i, j

      !$acc parallel loop collapse(2) default(present)
      do j = jb - jh, je + jh
         do i = ib - ih, ie + ih
            um_d(i, j, ke + 1) = 2*Uinf - um_d(i, j, ke)
            u0_d(i, j, ke + 1) = 2*Uinf - u0_d(i, j, ke)
            vm_d(i, j, ke + 1) = 2*Vinf - vm_d(i, j, ke)
            v0_d(i, j, ke + 1) = 2*Vinf - v0_d(i, j, ke)
         end do
      end do
      !$acc end parallel loop
   end subroutine valuetop_uv_device

   !> Flux top boundary for temperature.
   !!
   !! The two scalars hoisted out of the kernel are the parts of fluxtop's
   !! expression that do not depend on i or j. What is left inside is the same
   !! division by the same product the host forms, so the two agree to the bit.
   subroutine fluxtop_thl_device
      use modglobal,   only : ib, ie, ih, jb, je, jh, ke, ltempeq, dzf, dzh, dzhi, eps1
      use modcuda,     only : thl0_d, thlm_d, ekh_d
      use modsurfdata, only : wttop
      implicit none
      integer :: i, j
      real :: flux_factor, dzhi_top, ek_lower, ek_upper

      if (.not. ltempeq) return

      if (abs(wttop) .le. eps1) then
         !$acc parallel loop collapse(2) default(present)
         do j = jb - jh, je + jh
            do i = ib - ih, ie + ih
               thlm_d(i, j, ke + 1) = thlm_d(i, j, ke)
               thl0_d(i, j, ke + 1) = thl0_d(i, j, ke)
            end do
         end do
         !$acc end parallel loop
      else
         flux_factor = dzh(ke + 1)*wttop
         dzhi_top = dzhi(ke + 1)
         ek_lower = 0.5*dzf(ke)
         ek_upper = 0.5*dzf(ke + 1)
         !$acc parallel loop collapse(2) default(present)
         do j = jb - jh, je + jh
            do i = ib - ih, ie + ih
               thlm_d(i, j, ke + 1) = thlm_d(i, j, ke) + flux_factor / &
                  (dzhi_top*(ek_lower*ekh_d(i, j, ke + 1) + ek_upper*ekh_d(i, j, ke)))
               thl0_d(i, j, ke + 1) = thl0_d(i, j, ke) + flux_factor / &
                  (dzhi_top*(ek_lower*ekh_d(i, j, ke + 1) + ek_upper*ekh_d(i, j, ke)))
            end do
         end do
         !$acc end parallel loop
      end if
   end subroutine fluxtop_thl_device

   !> Fixed-value top boundary for temperature.
   subroutine valuetop_thl_device
      use modglobal,   only : ib, ie, ih, jb, je, jh, ke, ltempeq
      use modcuda,     only : thl0_d, thlm_d
      use modsurfdata, only : thl_top
      implicit none
      integer :: i, j

      if (.not. ltempeq) return

      !$acc parallel loop collapse(2) default(present)
      do j = jb - jh, je + jh
         do i = ib - ih, ie + ih
            thlm_d(i, j, ke + 1) = 2*thl_top - thlm_d(i, j, ke)
            thl0_d(i, j, ke + 1) = 2*thl_top - thl0_d(i, j, ke)
         end do
      end do
      !$acc end parallel loop
   end subroutine valuetop_thl_device

   !> Zero-gradient ghost levels for the kappa-scheme copy of temperature.
   subroutine topthl0c_device
      use modglobal, only : ib, ie, ihc, jb, je, jhc, ke, khc, &
                            ltempeq, iadv_thl, iadv_kappa
      use modcuda,   only : thl0c_d
      implicit none
      integer :: i, j, n

      if (.not. ltempeq) return
      if (iadv_thl /= iadv_kappa) return

      ! n runs outside the kernel: each level is filled from the one below it,
      ! which the previous pass wrote.
      do n = 1, khc
         !$acc parallel loop collapse(2) default(present)
         do j = jb - jhc, je + jhc
            do i = ib - ihc, ie + ihc
               thl0c_d(i, j, ke + n) = thl0c_d(i, j, ke + n - 1)
            end do
         end do
         !$acc end parallel loop
      end do
   end subroutine topthl0c_device

   !> Flux top boundary for moisture.
   subroutine fluxtop_qt_device
      use modglobal,   only : ib, ie, ih, jb, je, jh, ke, lmoist, dzf, dzh, dzhi, eps1
      use modcuda,     only : qt0_d, qtm_d, ekh_d
      use modsurfdata, only : wqtop
      implicit none
      integer :: i, j
      real :: flux_factor, dzhi_top, ek_lower, ek_upper

      if (.not. lmoist) return

      if (abs(wqtop) .le. eps1) then
         !$acc parallel loop collapse(2) default(present)
         do j = jb - jh, je + jh
            do i = ib - ih, ie + ih
               qtm_d(i, j, ke + 1) = qtm_d(i, j, ke)
               qt0_d(i, j, ke + 1) = qt0_d(i, j, ke)
            end do
         end do
         !$acc end parallel loop
      else
         flux_factor = dzh(ke + 1)*wqtop
         dzhi_top = dzhi(ke + 1)
         ek_lower = 0.5*dzf(ke)
         ek_upper = 0.5*dzf(ke + 1)
         !$acc parallel loop collapse(2) default(present)
         do j = jb - jh, je + jh
            do i = ib - ih, ie + ih
               qtm_d(i, j, ke + 1) = qtm_d(i, j, ke) + flux_factor / &
                  (dzhi_top*(ek_lower*ekh_d(i, j, ke + 1) + ek_upper*ekh_d(i, j, ke)))
               qt0_d(i, j, ke + 1) = qt0_d(i, j, ke) + flux_factor / &
                  (dzhi_top*(ek_lower*ekh_d(i, j, ke + 1) + ek_upper*ekh_d(i, j, ke)))
            end do
         end do
         !$acc end parallel loop
      end if
   end subroutine fluxtop_qt_device

   !> Fixed-value top boundary for moisture.
   subroutine valuetop_qt_device
      use modglobal,   only : ib, ie, ih, jb, je, jh, ke, lmoist
      use modcuda,     only : qt0_d, qtm_d
      use modsurfdata, only : qt_top
      implicit none
      integer :: i, j

      if (.not. lmoist) return

      !$acc parallel loop collapse(2) default(present)
      do j = jb - jh, je + jh
         do i = ib - ih, ie + ih
            qtm_d(i, j, ke + 1) = 2*qt_top - qtm_d(i, j, ke)
            qt0_d(i, j, ke + 1) = 2*qt_top - qt0_d(i, j, ke)
         end do
      end do
      !$acc end parallel loop
   end subroutine valuetop_qt_device

   !> Flux top boundary for the scalars.
   !!
   !! The i and j extents are the momentum halo, not the scalar one: that is
   !! what fluxtopscal writes, and the wider ghost columns of sv0 are left to
   !! the halo exchange.
   subroutine fluxtop_sv_device
      use modglobal,   only : ib, ie, ih, jb, je, jh, ke, khc, nsv, dzf, dzh, dzhi
      use modcuda,     only : sv0_d, svm_d, ekh_d
      use modsurfdata, only : wsvtop
      implicit none
      integer :: i, j, m, n
      real :: flux_factor, dzhi_top, ek_lower, ek_upper

      if (nsv < 1) return

      dzhi_top = dzhi(ke + 1)
      ek_lower = 0.5*dzf(ke)
      ek_upper = 0.5*dzf(ke + 1)

      do n = 1, nsv
         flux_factor = dzh(ke + 1)*wsvtop(n)
         !$acc parallel loop collapse(3) default(present)
         do m = 1, khc
            do j = jb - jh, je + jh
               do i = ib - ih, ie + ih
                  sv0_d(i, j, ke + m, n) = sv0_d(i, j, ke, n) + flux_factor / &
                     (dzhi_top*(ek_lower*ekh_d(i, j, ke + 1) + ek_upper*ekh_d(i, j, ke)))
                  svm_d(i, j, ke + m, n) = svm_d(i, j, ke, n) + flux_factor / &
                     (dzhi_top*(ek_lower*ekh_d(i, j, ke + 1) + ek_upper*ekh_d(i, j, ke)))
               end do
            end do
         end do
         !$acc end parallel loop
      end do
   end subroutine fluxtop_sv_device

   !> Fixed-value top boundary for the scalars.
   !!
   !! valuetopscal writes the full scalar halo, unlike fluxtopscal above.
   subroutine valuetop_sv_device
      use modglobal,   only : ib, ie, ihc, jb, je, jhc, ke, khc, nsv
      use modcuda,     only : sv0_d, svm_d
      use modsurfdata, only : sv_top
      implicit none
      integer :: i, j, m, n
      real :: val

      if (nsv < 1) return

      do n = 1, nsv
         val = sv_top(n)
         !$acc parallel loop collapse(3) default(present)
         do m = 1, khc
            do j = jb - jhc, je + jhc
               do i = ib - ihc, ie + ihc
                  sv0_d(i, j, ke + m, n) = 2*val - sv0_d(i, j, ke, n)
                  svm_d(i, j, ke + m, n) = 2*val - svm_d(i, j, ke, n)
               end do
            end do
         end do
         !$acc end parallel loop
      end do
   end subroutine valuetop_sv_device

   !> No flow through the bottom of the domain.
   subroutine bottom_w_device
      use modglobal, only : ib, ie, ih, jb, je, jh, kb
      use modcuda,   only : w0_d, wm_d
      implicit none
      integer :: i, j

      !$acc parallel loop collapse(2) default(present)
      do j = jb - jh, je + jh
         do i = ib - ih, ie + ih
            wm_d(i, j, kb) = 0.
            w0_d(i, j, kb) = 0.
         end do
      end do
      !$acc end parallel loop
   end subroutine bottom_w_device

   !> No flow through the top of the domain.
   subroutine topw_zero_device
      use modglobal, only : ib, ie, ih, jb, je, jh, ke
      use modcuda,   only : w0_d, wm_d
      implicit none
      integer :: i, j

      !$acc parallel loop collapse(2) default(present)
      do j = jb - jh, je + jh
         do i = ib - ih, ie + ih
            w0_d(i, j, ke + 1) = 0.0
            wm_d(i, j, ke + 1) = 0.0
         end do
      end do
      !$acc end parallel loop
   end subroutine topw_zero_device

   !> Background turbulence level in the ghost level above the domain.
   subroutine tope12_device
      use modglobal, only : ib, ie, ih, jb, je, jh, ke, e12min
      use modcuda,   only : e120_d, e12m_d
      implicit none
      integer :: i, j

      !$acc parallel loop collapse(2) default(present)
      do j = jb - jh, je + jh
         do i = ib - ih, ie + ih
            e120_d(i, j, ke + 1) = e12min
            e12m_d(i, j, ke + 1) = e12min
         end do
      end do
      !$acc end parallel loop
   end subroutine tope12_device
#endif


!>Set thl, qt and sv(n) equal to slab average at level kmax
! Think this can be removed
   subroutine tqaver

      use modmpi, only:comm3d, mpierr, my_real, mpi_sum
      use modglobal, only:ib, ie, jb, je, ke, nsv, rslabs
      use modfields, only:thl0, qt0, sv0
      implicit none

      real thl0a, qt0a
      real thl0al, qt0al
      integer n
      real, allocatable, dimension(:) :: sv0al, sv0a
      allocate (sv0al(nsv), sv0a(nsv))

      thl0al = sum(thl0(ib:ie, jb:je, ke))
      qt0al = sum(qt0(ib:ie, jb:je, ke))

      do n = 1, nsv
         sv0al(n) = sum(sv0(ib:ie, jb:je, ke, n))
      enddo

      call MPI_ALLREDUCE(thl0al, thl0a, 1, MY_REAL, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(qt0al, qt0a, 1, MY_REAL, &
                         MPI_SUM, comm3d, mpierr)
      if (nsv > 0) then
         call MPI_ALLREDUCE(sv0al, sv0a, nsv, MY_REAL, &
                            MPI_SUM, comm3d, mpierr)
      end if

      thl0a = thl0a/rslabs
      qt0a = qt0a/rslabs
      sv0a = sv0a/rslabs

      thl0(ib:ie, jb:je, ke) = thl0a
      qt0(ib:ie, jb:je, ke) = qt0a
      do n = 1, nsv
         sv0(ib:ie, jb:je, ke, n) = sv0a(n)
      enddo
      deallocate (sv0al, sv0a)

      return
   end subroutine tqaver

end module modboundary
