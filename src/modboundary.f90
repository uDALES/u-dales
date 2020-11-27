!> \file modboundary.f90
!>
!>
!! \par Revision list
!! All boundary conditions are in this file, except for immersed boundaries.
!! \par Authors
!!
module modboundary

   implicit none
   save
   private
   public :: initboundary, boundary, grwdamp, ksp, tqaver, &
             bcp, bcpup, closurebc
   integer :: ksp = -1 !<    lowest level of sponge layer
   real, allocatable :: tsc(:) !<   damping coefficients to be used in grwdamp.
   real :: rnu0 = 2.75e-3
contains
   !>
   !! Initializing Boundary; specifically the sponge layer
   !>
   subroutine initboundary
      use modglobal, only:ib, kb, ke, kh, kmax, pi, zf, iplane
      use modinletdata, only:irecy,irecydriver
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
      irecydriver = iplane! + ib

   end subroutine initboundary

   !>
   !! Execute boundary conditions
   subroutine boundary

      use modglobal, only:ib, ie, ih, jb, je, jgb, jge, jh, kb, ke, kh, linoutflow, dzf, zh, dy, &
         timee, ltempeq, lmoist, BCxm, BCym, BCxT, BCyT, BCxq, BCyq, BCxs, BCys, BCtopm, BCtopT,&
         BCtopq, BCtops, e12min, idriver, luvolflowr, luoutflowr
      use modfields, only:u0, v0, w0, um, vm, wm, thl0, thlm, qt0, qtm, uout, uouttot, e120, e12m,&
                          u0av
      use modsubgriddata, only:ekh, ekm
      use modsurfdata, only:thl_top, qt_top, sv_top, wttop, wqtop, wsvtop
      use modmpi, only:myid, slabsum
      use modinlet, only:inletgen, inletgennotemp
      use moddriver, only : drivergen
      use modinletdata, only:irecy, ubulk, iangle
!    use modsurface, only : getobl
      implicit none
      real, dimension(kb:ke) :: uaverage
      integer i, k

     ! if not using massflowrate need to set outflow velocity
     if (luoutflowr) then
        ! do nothing - calculated in modforces
     elseif (.not. luvolflowr) then
        !ubulk = sum(u0av)/(ke-kb+1)
        do k = kb, ke
           uaverage(k) = u0av(k)*dzf(k)                                                        
        end do
        ! need a method to know if we have all blocks at lowest cell kb
        ! assuming this for now (hence kb+1)
        uouttot = sum(uaverage(kb:ke))/(zh(ke + 1) - zh(kb+1)) 
     else
        uouttot = ubulk
     end if

     !BCxm!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !periodic or inflow/outflow conditions for momentum
      if (BCxm .eq. 1) then  !periodic
         call cyclicmi

         if (idriver == 1) then ! write driver files
           call drivergen
         end if

      else if (BCxm .eq. 2) then !previously iinletgen 1
         uouttot = cos(iangle)*ubulk
         if (ltempeq) then
            call inletgen
         else
            call inletgennotemp
         end if

         ! iolet - called due to BCtopm = 3

      else if (BCxm .eq. 3) then ! previously iinletgen 2
         uouttot = cos(iangle)*ubulk
         if (ltempeq) then
            call inletgen
         else
            call inletgennotemp
         end if

         ! iolet - called due to BCtopm = 3

      else if (BCxm .eq. 4) then !previously (inoutflow without iinlet)
         uouttot = cos(iangle)*ubulk
         if (ltempeq) then
            call inletgen
         else
            call inletgennotemp
         end if

         ! iolet - called due to BCtopm = 3

      else if (BCxm .eq. 5) then ! driver from drivergen (idriver == 2)

         uouttot = ubulk ! does this hold for all forcings of precursor simulations? tg3315

         call drivergen

         ! iolet - called due to BCtopm = 3

      else
         write(0, *) "ERROR: lateral boundary type for veloctiy in x-direciton undefined"
         stop 1
      end if

      !BCym!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !currently BC in y is always periodic for momentum
      if (BCym .eq. 1) then
         call cyclicmj
      else
         write(0, *) "ERROR: lateral boundary type for velocity in y-direction undefined"
         stop 1
      end if

      !BCxT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (BCxT .eq. 1) then
         call cyclichi

      else if (BCxT .eq. 2) then !inoutflow - will be overwritten unless BCxm == 1
         call iohi    ! make sure uouttot is known and realistic
      else if (BCxT .eq. 3) then
         !do nothing, temperature is considered in iolet
      else
         write(0, *) "ERROR: lateral boundary type for temperature in x-direction undefined"
         stop 1
      end if

      !BCyT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (BCyT .eq. 1) then
         call cyclichj
      else
         write(0, *) "ERROR: lateral boundary type for temperature in y-direction undefined"
         stop 1
      end if

      !BCxq!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (BCxq .eq. 1) then
         call cyclicqi
      else if (BCxq .eq. 2) then !inoutflow  - will be overwritten unless BCxm == 1
        call ioqi ! tg3315 - make sure uouttot is known and realistic
      elseif (BCxq .eq. 3) then 
        !do nothing, temperature is considered in iolet
      else
         write(0, *) "ERROR: lateral boundary type for humidity in x-direction undefined"
         stop 1
      end if

      !BCyq!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (BCyq .eq. 1) then
         call cyclicqj
      else
         write(0, *) "ERROR: lateral boundary type for humidity in y-direction undefined"
         stop 1
      end if

      !BCys!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (BCys .eq. 1) then
         call cyclicsj
      elseif (BCys .eq. 5) then
         ! done in scalSIRANE
      else
         write(0, *) "ERROR: lateral boundary type for scalars in y-direction undefined"
         stop 1
      end if

      !BCxs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (BCxs .eq. 1) then
         call cyclicsi
      else if (BCxs .eq. 2) then !inoutflow  - will be overwritten unless BCxm == 1
         call iosi ! make sure uouttot is known and correct for the running set-up
      else if (BCxs .eq. 3) then
         ! do nothing - considered in iolet

      else if (BCxs .eq. 4) then !scalrec - will be overwritten unless BCxm == 1
         call scalrec

      else if (BCxs .eq. 5) then !previously SIRANE - will be overwritten unless BCxm == 1
         call scalSIRANE !  make sure uouttot/ vouttot is known and realistic

      else
         write(0, *) "ERROR: lateral boundary type for scalars in x-direction undefined"
         stop 1
      end if

      !BCtopm!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (BCtopm .eq. 1) then
         !free-slip = zero-flux
         call fluxtop(um, ekm, 0.0)
         call fluxtop(u0, ekm, 0.0)
         call fluxtop(vm, ekm, 0.0)
         call fluxtop(v0, ekm, 0.0)
         e120(:, :, ke + 1) = e12min ! free slip top wall
         e12m(:, :, ke + 1) = e12min
         w0(:, :, ke + 1) = 0.0
         wm(:, :, ke + 1) = 0.0
      else if (BCtopm .eq. 2) then
         !no-slip = zero velocity at wall
         call valuetop(um, 0.0)
         call valuetop(u0, 0.0)
         call valuetop(vm, 0.0)
         call valuetop(v0, 0.0)
         w0(:, :, ke + 1) = 0.0
         wm(:, :, ke + 1) = 0.0
      else if (BCtopm .eq. 3) then
         call fluxtop(um, ekm, 0.0)
         call fluxtop(u0, ekm, 0.0)
         call fluxtop(vm, ekm, 0.0)
         call fluxtop(v0, ekm, 0.0)
         e120(:, :, ke + 1) = e12min ! free slip top wall
         e12m(:, :, ke + 1) = e12min
         if (idriver==2) then ! does not use ddispdx, Uinf etc.
           w0(:, :, ke + 1) = 0.0
           wm(:, :, ke + 1) = 0.0
         else
           call inlettop ! for iinletgen...
         end if
         call iolet  !ils13, 13.8.18: iolet also deals with lateral boundaries!!
      else
         write(0, *) "ERROR: top boundary type for velocity undefined"
         stop 1
      end if

      !BCtopT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (BCtopT .eq. 1) then
         call fluxtop(thlm, ekh, wttop)
         call fluxtop(thl0, ekh, wttop)
      else if (BCtopT .eq. 2) then
         call valuetop(thlm, thl_top)
         call valuetop(thl0, thl_top)
      else
         write(0, *) "ERROR: top boundary type for temperature undefined"
         stop 1
      end if

      !BCtopq!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (BCtopq .eq. 1) then
         call fluxtop(qtm, ekh, wqtop)
         call fluxtop(qt0, ekh, wqtop)
      else if (BCtopq .eq. 2) then
         call valuetop(qtm, qt_top)
         call valuetop(qt0, qt_top)
      else
         write(0, *) "ERROR: top boundary type for humidity undefined"
         stop 1
      end if

      !BCtops!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (BCtops .eq. 1) then
         call fluxtopscal(wsvtop)
         call fluxtopscal(wsvtop)
      else if (BCtops .eq. 2) then
         call valuetopscal(sv_top)
         call valuetopscal(sv_top)
      else
         write(0, *) "ERROR: top boundary type for scalars undefined"
         stop 1
      end if

   end subroutine boundary

   subroutine closurebc
      use modsubgriddata, only:ekm, ekh
      use modglobal, only:ib, ie, jb, je, kb, ke, ih, jh, kh, numol, prandtlmoli, linoutflow, BCtopm
      use modmpi, only:excjs
      integer i, j

      ! Top and bottom
      ! ils13, 13.8.18: what should it be for slip or mixed BCs ?
      if ((BCtopm.eq.1) .or. (BCtopm.eq.3)) then !free-slip
         do j = jb - 1, je + 1
            do i = ib - 1, ie + 1
               ekm(i, j, ke + 1) = ekm(i, j, ke) ! zero-gradient top wall
               ekh(i, j, ke + 1) = ekh(i, j, ke) ! zero-gradient top wall
               ekm(i, j, kb - 1) = 2.*numol - ekm(i, j, kb) ! no-slip lower wall
               ekh(i, j, kb - 1) = (2.*numol*prandtlmoli) - ekh(i, j, kb) ! no-slip lower wall
            end do
         end do
      else if (BCtopm.eq.2)  then !no-slip
         do j = jb - 1, je + 1
            do i = ib - 1, ie + 1
               ekm(i, j, ke + 1) = 2.*numol - ekm(i, j, ke) ! no-slip top wall
               ekh(i, j, ke + 1) = (2.*numol*prandtlmoli) - ekh(i, j, ke) ! no-slip top wall
               ekm(i, j, kb - 1) = 2.*numol - ekm(i, j, kb) ! no-slip lower wall
               ekh(i, j, kb - 1) = (2.*numol*prandtlmoli) - ekh(i, j, kb) ! no-slip lower wall
            end do
         end do
      end if

      ! horizontal BC's
      if (linoutflow) then ! inflow/outflow
         ekm(ib - 1, :, :) = ekm(ib, :, :)
         ekm(ie + 1, :, :) = ekm(ie, :, :)
         ekh(ib - 1, :, :) = ekh(ib, :, :)
         ekh(ie + 1, :, :) = ekh(ie, :, :)
      else
         ekm(ib - 1, :, :) = ekm(ie, :, :) ! periodic
         ekm(ie + 1, :, :) = ekm(ib, :, :)
         ekh(ib - 1, :, :) = ekh(ie, :, :)
         ekh(ie + 1, :, :) = ekh(ib, :, :)
      end if

      call excjs(ekm, ib, ie, jb, je, kb - kh, ke + kh, ih, jh)
      call excjs(ekh, ib, ie, jb, je, kb - kh, ke + kh, ih, jh)

   end subroutine closurebc


   !> Sets lateral periodic boundary conditions for the scalars
   subroutine iosi
      use modglobal, only:ib, ie, jb, je, ih, jh, kb, ke, kh, nsv, dt, rk3step, dxhi, ltempeq, &
         ihc, jhc, khc, dy
      use modfields, only:sv0, svm, svprof, uouttot
      use modinletdata, only:ubulk
      real rk3coef
      integer k, n, m

      rk3coef = dt/(4.-dble(rk3step))
      do n = 1, nsv
         do k = kb, ke + 1
            sv0(ib - 1, :, k, n) = 2*svprof(k, n) - sv0(ib, :, k, n)
            svm(ib - 1, :, k, n) = 2*svprof(k, n) - svm(ib, :, k, n)
         end do
         sv0(ie + 1, :, :, n) = sv0(ie, :, :, n) - (sv0(ie + 1, :, :, n) - sv0(ie, :, :, n))*dxhi(ie + 1)*rk3coef*uouttot  ! tg3315 should be uouttot and will have to change depending on forcing
         svm(ie + 1, :, :, n) = svm(ie, :, :, n) - (svm(ie + 1, :, :, n) - svm(ie, :, :, n))*dxhi(ie + 1)*rk3coef*uouttot
      enddo

      return
   end subroutine iosi

   subroutine scalrec
      use modglobal, only:ib, ie, jb, je, ih, jh, kb, ke, kh, nsv, dt, rk3step, dxhi, ltempeq, &
         ihc, jhc, khc, dy
      use modfields, only:sv0, svm, svprof, uouttot, um, u0, vm, v0
      use modinletdata, only:ubulk
      real rk3coef
      integer k, n, m
      ! recycling method for scalar fields following Matheou and Bowman (2015)

      if (nsv > 0) then
         rk3coef = dt/(4.-dble(rk3step))
         do m = 1, ihc ! loop over virtual cells
            do n = 1, nsv - 1
               sv0(ib - m, :, :, n + 1) = sv0(ie + 1 - m, :, :, n)
               sv0(ie + m, :, :, n) = sv0(ib - 1 + m, :, :, n + 1)
               svm(ib - m, :, :, n + 1) = svm(ie + 1 - m, :, :, n)
               svm(ie + m, :, :, n) = svm(ib - 1 + m, :, :, n + 1)
            end do

            ! zero conc. on scalar 1 !tg3315 should be changed to as above in
            sv0(ib - m, :, :, 1) = 0.
            svm(ib - m, :, :, 1) = 0.

            ! DIY outflow BC (advection step as linout) tg3315
          sv0(ie+m,:,:,nsv)=sv0(ie+1-m,:,:,nsv)-(sv0(ie+m,:,:,nsv)-sv0(ie+1-m,:,:,nsv))*dxhi(ie+m)*rk3coef*uouttot
          svm(ie+m,:,:,nsv)=svm(ie+1-m,:,:,nsv)-(svm(ie+m,:,:,nsv)-svm(ie+1-m,:,:,nsv))*dxhi(ie+m)*rk3coef*uouttot
         end do
      end if
      return
   end subroutine scalrec

   subroutine scalSIRANE
    use modglobal, only : ib,ie,jb,je,ih,jh,kb,ke,kh,nsv,dt,lscalrec,lmoistinout,ltempinout,rk3step,dxhi,ltempeq,&
         ihc, jhc, khc, lSIRANEinout, dy
      use modfields, only:sv0, svm, svprof
      use modinletdata, only:ubulk
      use modmpi, only:myid, nprocs
      real rk3coef
      integer k, n, m
      if (nsv > 0) then
         !rk3coef = dt / (4. - dble(rk3step))
         do n = 1, nsv
            do m = 1, ihc
               do k = kb, ke + 1
                  sv0(ib - m, :, k, n) = 2*svprof(k, n) - sv0(ib - m + 1, :, k, n) !scalars have two ghost cells...???
                  svm(ib - m, :, k, n) = 2*svprof(k, n) - svm(ib - m + 1, :, k, n)
               end do
!              sv0(ie+m,:,:,n)= sv0(ie+m-1,:,:,n) - (sv0(ie+m,:,:,n)-sv0(ie+m-1,:,:,n))*dxhi(ie+m)*rk3coef*ubulk
!              svm(ie+m,:,:,n)= svm(ie+m-1,:,:,n) - (svm(ie+m,:,:,n)-svm(ie+m-1,:,:,n))*dxhi(ie+m)*rk3coef*ubulk !changed from uouttot to ubulk here !tg3315 08/11/2017
!              sv0(ie+m,:,:,n)= sv0(ie+m-1,:,:,n) - (sv0(ie+m,:,:,n)-sv0(ie+m-1,:,:,n))*dxhi(ie+m)*rk3coef*u0(ie+m,:,:)
!              svm(ie+m,:,:,n)= svm(ie+m-1,:,:,n) - (svm(ie+m,:,:,n)-svm(ie+m-1,:,:,n))*dxhi(ie+m)*rk3coef*um(ie+m,:,:) !changed from uouttot to ubulk here !tg3315 08/11/2017
               svm(ie + m, :, :, n) = svm(ie + m - 1, :, :, n)
               sv0(ie + m, :, :, n) = sv0(ie + m - 1, :, :, n)
            end do !m, ihc
         end do !n, nsv

         do m = 1, jhc
            do n = 1, nsv
            if (myid == 0) then
               do k = kb, ke + 1
                  sv0(:, jb - m, k, n) = 2*svprof(k, n) - sv0(:, jb - m + 1, k, n)
                  svm(:, jb - m, k, n) = 2*svprof(k, n) - svm(:, jb - m + 1, k, n)
               end do
            end if
            if (myid == nprocs - 1) then
               !sv0(:,je+m,:,n)= sv0(:,je+m-1,:,n) - (sv0(:,je+m,:,n)-sv0(:,je+m-1,:,n))*dy*rk3coef*ubulk
               !svm(:,je+m,:,n)= svm(:,je+m-1,:,n) - (svm(:,je+m,:,n)-svm(:,je+m-1,:,n))*dy*rk3coef*ubulk !changed from uouttot to ubulk here !tg3315 08/11/2017
               !      sv0(:,je+m,:,n)= sv0(:,je+m-1,:,n) - (sv0(:,je+m,:,n)-sv0(:,je+m-1,:,n))*dy*rk3coef*v0(:,je+m,:)
               !      svm(:,je+m,:,n)= svm(:,je+m-1,:,n) - (svm(:,je+m,:,n)-svm(:,je+m-1,:,n))*dy*rk3coef*vm(:,je+m,:) !changed from uouttot to ubulk here !tg3315 08/11/2017
               svm(:, je + m, :, n) = svm(:, je + m - 1, :, n)
               sv0(:, je + m, :, n) = sv0(:, je + m - 1, :, n)
            end if
            end do !n, nsv
         end do !m, jhc

      end if !nsv>0
      return
   end subroutine scalSIRANE

   !!!!!!!!!!! x/i periodic BC for scalars!!!!!!!!!!!
   !> Sets x/i periodic boundary conditions for the temperature
   subroutine cyclichi
      use modglobal, only:ib, ie, jb, je, ih, jh, kb, ke, kh, nsv, dt, rk3step, dxhi, ihc, jhc, khc, dy
      use modfields, only:thl0, thlm
      integer m

      do m = 1, ih
         thl0(ib - m, :, :) = thl0(ie + 1 - m, :, :)
         thl0(ie + m, :, :) = thl0(ib - 1 + m, :, :)
         thlm(ib - m, :, :) = thlm(ie + 1 - m, :, :)
         thlm(ie + m, :, :) = thlm(ib - 1 + m, :, :)
      end do

      return
   end subroutine cyclichi

   !> Sets x/i periodic boundary conditions for the humidity
   subroutine cyclicqi
      use modglobal, only:ib, ie, jb, je, ih, jh, kb, ke, kh, nsv, dt, rk3step, dxhi, ihc, jhc, khc, dy
      use modfields, only:qt0, qtm
      integer m

      do m = 1, ih
         qt0(ib - m, :, :) = qt0(ie + 1 - m, :, :)
         qt0(ie + m, :, :) = qt0(ib - 1 + m, :, :)
         qtm(ib - m, :, :) = qtm(ie + 1 - m, :, :)
         qtm(ie + m, :, :) = qtm(ib - 1 + m, :, :)
      end do

      return
   end subroutine cyclicqi

   !> Sets x/iperiodic boundary conditions for the scalars
   subroutine cyclicsi
      use modglobal, only:ib, ie, jb, je, ih, jh, kb, ke, kh, nsv, dt, rk3step, dxhi, ihc, jhc, khc, dy
      use modfields, only:sv0, svm
      integer m

      do m = 1, ihc
         sv0(ib - m, :, :, :) = sv0(ie + 1 - m, :, :, :)
         sv0(ie + m, :, :, :) = sv0(ib - 1 + m, :, :, :)
         svm(ib - m, :, :, :) = svm(ie + 1 - m, :, :, :)
         svm(ie + m, :, :, :) = svm(ib - 1 + m, :, :, :)
      end do

      return
   end subroutine cyclicsi
   !!!!!!!!!!!end x/i periodic BC for scalars!!!!!!!!

   !> Sets x/inlet-outlet boundary conditions for moisture
   subroutine ioqi
     use modglobal, only: ib, ie, jb, je, ih, jh, kb, ke, kh, dxhi, rk3step, dt
     use modfields, only: qt0, qtm, qtprof, uouttot
     use modinletdata, only: ubulk
     integer k,j
     real rk3coef                                                                                   

     rk3coef = dt/(4.-dble(rk3step))

     do k = kb, ke
       do j = jb, je
         qt0(ib - 1, j, k) = 2*qtprof(k) - qt0(ib, j, k) !watch!
         qtm(ib - 1, j, k) = 2*qtprof(k) - qtm(ib, j, k)
       end do
    end do
    
    !uouttot is zero unless lmassflowr 
    qt0(ie + 1, :, :) = qt0(ie, :, :) - (qt0(ie + 1, :, :) - qt0(ie, :, :))*dxhi(ie + 1)*rk3coef*uouttot ! tg3315 should be uouttot and will have to change depending on forcing
    qtm(ie + 1, :, :) = qtm(ie, :, :) - (qtm(ie + 1, :, :) - qtm(ie, :, :))*dxhi(ie + 1)*rk3coef*uouttot

   end subroutine ioqi
   
   !> Sets x/in;et-outlet boundary conditions for temperature
   subroutine iohi
     use modglobal, only: ib, ie, jb, je, ih, jh, kb, ke, kh, dxhi, rk3step, dt
     use modfields, only: thl0, thlm, thlprof, uouttot
     use modinletdata, only: ubulk
     integer k,j   
     real rk3coef                                                                                  

     rk3coef = dt/(4.-dble(rk3step))

     do k = kb, ke
       do j = jb, je
         thl0(ib - 1, j, k) = 2*thlprof(k) - thl0(ib, j, k) !watch!
         thlm(ib - 1, j, k) = 2*thlprof(k) - thlm(ib, j, k)
       end do
    end do

    thl0(ie + 1, :, :) = thl0(ie, :, :) - (thl0(ie + 1, :, :) - thl0(ie, :, :))*dxhi(ie + 1)*rk3coef*uouttot ! tg3315 should be uouttot and will have to change depending on forcing
    thlm(ie + 1, :, :) = thlm(ie, :, :) - (thlm(ie + 1, :, :) - thlm(ie, :, :))*dxhi(ie + 1)*rk3coef*uouttot

   end subroutine iohi

   !!!!!!!!!!! y/j periodic BC for scalars!!!!!!!!!!!
   !> Sets y/j periodic boundary conditions for the temperature
   subroutine cyclichj
      use modglobal, only:ib, ie, jb, je, ih, jh, kb, ke, kh, nsv, dt, rk3step, dxhi, ihc, jhc, khc, dy
      use modfields, only:thl0, thlm
      use modmpi, only:excjs, myid, nprocs

      call excjs(thl0, ib, ie, jb, je, kb, ke + kh, ih, jh)
      call excjs(thlm, ib, ie, jb, je, kb, ke + kh, ih, jh)

      return
   end subroutine cyclichj

   !> Sets y/j periodic boundary conditions for the humidity
   subroutine cyclicqj
      use modglobal, only:ib, ie, jb, je, ih, jh, kb, ke, kh, nsv, dt, rk3step, dxhi, ihc, jhc, khc, dy
      use modfields, only:qt0, qtm
      use modmpi, only:excjs, myid, nprocs

      call excjs(qt0, ib, ie, jb, je, kb, ke + kh, ih, jh)
      call excjs(qtm, ib, ie, jb, je, kb, ke + kh, ih, jh)

      return
   end subroutine cyclicqj

   !> Sets y/j periodic boundary conditions for the scalars
   subroutine cyclicsj
      use modglobal, only:ib, ie, jb, je, ih, jh, kb, ke, kh, nsv, dt, rk3step, dxhi, ihc, jhc, khc, dy
      use modfields, only:sv0, svm
      use modmpi, only:excjs, myid, nprocs
      integer n

      do n = 1, nsv
         call excjs(sv0(:, :, :, n), ib, ie, jb, je, kb - khc, ke + khc, ihc, jhc)
         call excjs(svm(:, :, :, n), ib, ie, jb, je, kb - khc, ke + khc, ihc, jhc)
      enddo

      return
   end subroutine cyclicsj
   !!!!!!!!!!!end y/j periodic BC for scalars!!!!!!!!

   !>set lateral periodic boundary conditions for momentum in x/i direction
   subroutine cyclicmi

      use modglobal, only:ib, ie, jb, je, ih, jh, kb, ke, kh, jmax
      use modfields, only:u0, um, v0, vm, w0, wm, e120, e12m
      use modsubgriddata, only:loneeqn, lsmagorinsky

      integer n, m

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

         e120(ib - m, :, :) = e120(ie + 1 - m, :, :)
         e120(ie + m, :, :) = e120(ib - 1 + m, :, :)
         e12m(ib - m, :, :) = e12m(ie + 1 - m, :, :)
         e12m(ie + m, :, :) = e12m(ib - 1 + m, :, :)

      end do

      if (loneeqn) then
         e120(ib - m, :, :) = e120(ie + 1 - m, :, :)
         e120(ie + m, :, :) = e120(ib - 1 + m, :, :)
         e12m(ib - m, :, :) = e12m(ie + 1 - m, :, :)
         e12m(ie + m, :, :) = e12m(ib - 1 + m, :, :)
      end if

      return
   end subroutine cyclicmi

   !>set lateral periodic boundary conditions for momentum in y/j direction
   subroutine cyclicmj

      use modglobal, only:ib, ie, jb, je, ih, jh, kb, ke, kh, jmax
      use modfields, only:u0, um, v0, vm, w0, wm, e120, e12m, shear
      use modsubgriddata, only:loneeqn, lsmagorinsky
      use modmpi, only:excjs

      integer n, m

      call excjs(u0, ib, ie, jb, je, kb, ke + kh, ih, jh)
      call excjs(v0, ib, ie, jb, je, kb, ke + kh, ih, jh)
      call excjs(w0, ib, ie, jb, je, kb, ke + kh, ih, jh)
      call excjs(um, ib, ie, jb, je, kb, ke + kh, ih, jh)
      call excjs(vm, ib, ie, jb, je, kb, ke + kh, ih, jh)
      call excjs(wm, ib, ie, jb, je, kb, ke + kh, ih, jh)

      if (loneeqn) then
         call excjs(e120, ib, ie, jb, je, kb, ke + kh, ih, jh)
         call excjs(e12m, ib, ie, jb, je, kb, ke + kh, ih, jh)
         ! exchange shear components between processors
         do n = 1, 12 ! for all 12 components
            call excjs(shear(:, :, :, n), ib, ie, jb, je, kb, ke, 0, 1)
         end do
      end if

      if (lsmagorinsky) then
         ! exchange shear components between processors
         do n = 1, 12 ! for all 12 components
            call excjs(shear(:, :, :, n), ib, ie, jb, je, kb, ke, 0, 1)
         end do
      end if

      return
   end subroutine cyclicmj

   !>set inlet and outlet boundary conditions in i-direction
   subroutine iolet

      use modglobal, only:dxhi, dxhci, xh, zh, ib, ie, jb, je, ih, jh, kb, ke, kh, nsv, rk3step, dt, iinletgen, ltempeq, lmoist, ihc, idriver, dy, dzf, jtot, zh
      use modfields, only:u0, um, v0, vm, w0, wm, e120, e12m, thl0, thlm, qt0, qtm, sv0, svm, uprof, vprof, e12prof, thlprof, &
         qtprof, svprof, uouttot, wouttot
      use modmpi, only:excjs, myid, slabsum
      use modinletdata, only:u0inletbcold, v0inletbcold, w0inletbcold, uminletbc, vminletbc, wminletbc, totaluold, &
         t0inletbcold, tminletbc, u0driver, v0driver, w0driver, e120driver, thl0driver, qt0driver, umdriver, vmdriver, wmdriver,& 
         e12mdriver, thlmdriver, qtmdriver, sv0driver, svmdriver

      real rk3coef
      real, dimension(kb:ke) :: uin
      integer n, i, j, k, m

      rk3coef = dt/(4.-dble(rk3step))

      ! Inlet boundary is located at ib (not ib-1)!

      ! Inlet
      if ((iinletgen == 1) .or. (iinletgen == 2)) then
         do j = jb, je
            do k = kb, ke
               u0(ib, j, k) = u0inletbcold(j, k)
               um(ib, j, k) = uminletbc(j, k)
               u0(ib - 1, j, k) = 2*u0(ib, j, k) - u0(ib + 1, j, k)
               um(ib - 1, j, k) = 2*um(ib, j, k) - um(ib + 1, j, k)

               v0(ib - 1, j, k) = v0inletbcold(j, k)
               vm(ib - 1, j, k) = vminletbc(j, k)

               ! to be changed in the future: e12 should be taken from recycle plane!
               e120(ib - 1, j, k) = e120(ib, j, k) ! extrapolate e12 from interior
               e12m(ib - 1, j, k) = e12m(ib, j, k) ! extrapolate e12 from interior

               do n = 1, nsv
                  do m = 1, ihc
                     sv0(ib - m, j, k, n) = 2*svprof(k, n) - sv0(ib + (m - 1), j, k, n)
                     svm(ib - m, j, k, n) = 2*svprof(k, n) - svm(ib + (m - 1), j, k, n)
                  enddo
               enddo
            end do
            do k = kb, ke + 1
               w0(ib - 1, j, k) = w0inletbcold(j, k)
               wm(ib - 1, j, k) = wminletbc(j, k)
            end do
         end do

         ! Heat
         if (ltempeq) then
            do k = kb, ke
               do j = jb, je
                  thl0(ib - 1, j, k) = t0inletbcold(j, k)
                  thlm(ib - 1, j, k) = tminletbc(j, k)
               end do
            end do
         end if

         if (lmoist) then
            do k = kb, ke
               do j = jb, je
                  qt0(ib - 1, j, k) = 2*qtprof(k) - qt0(ib, j, k) !watch!
                  qtm(ib - 1, j, k) = 2*qtprof(k) - qtm(ib, j, k)
               end do
            end do
         end if

    ! Driver inlet
    elseif (idriver == 2) then
       do j=jb-1,je+1
         do k=kb,ke !tg3315 removed +1 following above...
           u0(ib,j,k)=u0driver(j,k) !max(0.,u0driver(j,k))
           um(ib,j,k)=umdriver(j,k) !max(0.,umdriver(j,k))
           u0(ib-1,j,k)= u0driver(j,k) !max(0.,2.*u0(ib,j,k)-u0(ib+1,j,k))
           um(ib-1,j,k)= umdriver(j,k)  !max(0.,2.*um(ib,j,k)-um(ib+1,j,k))

           v0(ib,j,k)   = v0driver(j,k) !max(0.,v0driver(j,k))
           vm(ib,j,k)   = vmdriver(j,k) !max(0.,vmdriver(j,k))
           v0(ib-1,j,k)   = v0driver(j,k) !max(0.,v0driver(j,k))
           vm(ib-1,j,k)   = vmdriver(j,k) !max(0.,vmdriver(j,k))

           ! to be changed in the future: e12 should be taken from recycle plane!
           !e120(ib-1,j,k) = e120driver(j,k)      ! extrapolate e12 from interior
           !e12m(ib-1,j,k) = e12mdriver(j,k)      ! extrapolate e12 from interior
             
           do n=1,nsv
              do m = 1,ihc
                 sv0(ib-m,j,k,n) = sv0driver(j,k,n)
                 svm(ib-m,j,k,n) = svmdriver(j,k,n)
                 !sv0(ib-m,j,k,n) = 2*svprof(k,n) - sv0(ib+(m-1),j,k,n)
                 !svm(ib-m,j,k,n) = 2*svprof(k,n) - svm(ib+(m-1),j,k,n)
              enddo
              sv0(ib,j,k,n) = sv0driver(j,k,n)
              svm(ib,j,k,n) = svmdriver(j,k,n)
           enddo
         end do
         do k=kb,ke+1
           w0(ib-1,j,k)   = w0driver(j,k) !max(0.,w0driver(j,k))
           wm(ib-1,j,k)   = wmdriver(j,k) !max(0.,wmdriver(j,k))
           w0(ib,j,k)   = w0driver(j,k) !max(0.,w0driver(j,k))
           wm(ib,j,k)   = wmdriver(j,k) !max(0.,wmdriver(j,k))
         end do
       end do

       ! Heat
       if (ltempeq ) then
          do j=jb-1,je+1
             do k=kb,ke+1
                thl0(ib,j,k) = thl0driver(j,k)
                thlm(ib,j,k) = thlmdriver(j,k)
                thl0(ib-1,j,k) = thl0driver(j,k)
                thlm(ib-1,j,k) = thlmdriver(j,k)
                !thlm(ib-1,j,k) = 2*thlm(ib,j,k) - thlm(ib+1,j,k)
                !thl0(ib-1,j,k) = 2*thl0(ib,j,k) - thl0(ib+1,j,k)
             end do
          end do
       end if

       if (lmoist ) then
          do j=jb-1,je+1
             do k=kb,ke+1
                qt0(ib,j,k) = qt0driver(j,k)
                ! qt0(ib-1,j,k) = 2*qtprof(k) - qt0(ib,j,k)
                qtm(ib,j,k) = qtmdriver(j,k)
                ! qtm(ib-1,j,k) = 2*qtprof(k) - qtm(ib,j,k)
                qt0(ib-1,j,k) = qt0driver(j,k)
                ! qt0(ib-1,j,k) = 2*qtprof(k) - qt0(ib,j,k)  !watch!
                qtm(ib-1,j,k) = qtmdriver(j,k)
                ! qtm(ib-1,j,k) = 2*qtprof(k) - qtm(ib,j,k)
                ! qt0(ib-1,j,k) = 2*qt0(ib,j,k) - qt0(ib+1,j,k)
                ! qtm(ib-1,j,k) = 2*qtm(ib,j,k) - qtm(ib+1,j,k)
             end do
          end do
       end if

       else ! (if iinetgen==0)

         do j = jb - 1, je + 1
            do k = kb, ke + 1
               ! Momentum
               u0(ib, j, k) = uprof(k)
               um(ib, j, k) = uprof(k)

               v0(ib - 1, j, k) = 2*vprof(k) - v0(ib, j, k) ! (v(ib)+v(ib-1))/2 = vprof
               vm(ib - 1, j, k) = 2*vprof(k) - vm(ib, j, k) ! (v(ib)+v(ib-1))/2 = vprof

               e120(ib - 1, j, k) = 2*e12prof(k) - e120(ib, j, k) ! (e12(ib)+e12(ib-1))/2=e12prof
               e12m(ib - 1, j, k) = 2*e12prof(k) - e12m(ib, j, k) ! (e12(ib)+e12(ib-1))/2=e12prof

               do n = 1, nsv
                  do m = 1, ihc
                     sv0(ib - m, j, k, n) = 2*svprof(k, n) - sv0(ib + (m - 1), j, k, n)
                     svm(ib - m, j, k, n) = 2*svprof(k, n) - svm(ib + (m - 1), j, k, n)
                  end do
               end do
               !end if

            enddo
         enddo

         ! Heat
         if (ltempeq) then
            do j = jb - 1, je + 1
               do k = kb, ke + 1
                  thl0(ib - 1, j, k) = 2*thlprof(k) - thl0(ib, j, k)
                  thlm(ib - 1, j, k) = 2*thlprof(k) - thlm(ib, j, k)
               end do
            end do
         end if

         if (lmoist) then
            do j = jb - 1, je + 1
               do k = kb, ke + 1
                  qt0(ib - 1, j, k) = 2*qtprof(k) - qt0(ib, j, k)
                  qtm(ib - 1, j, k) = 2*qtprof(k) - qtm(ib, j, k)
               end do
            end do
         end if

         u0(ib - 1, :, :) = 2*u0(ib, :, :) - u0(ib + 1, :, :) ! (u(ib+1)+u(ib-1))/2 = u(ib)
         um(ib - 1, :, :) = 2*um(ib, :, :) - um(ib + 1, :, :) ! (u(ib+1)+u(ib-1))/2 = u(ib)
         w0(ib - 1, :, :) = -w0(ib, :, :) ! (w(ib)+w(ib-1))/2 = 0
         wm(ib - 1, :, :) = -wm(ib, :, :)
      end if ! iinletgen==1 .or. iinletgen==2

      ! tg3315 added to ensure that uouttot matches driven inflow regardless of forcing
      ! set up assuming we have a block at lowest cell kb
      if (idriver==2) then
         uin = 0.
         uouttot = 0.
         call slabsum(uin,kb,ke,um,ib-ih,ie+ih,jb-jh,je+jh,kb-kh,ke+kh,ie,ie,jb,je,kb,ke) ! determine horizontal (j) average outflow velocity old
         do k=kb,ke
            uin(k) = uin(k)*dzf(k)*dy  ! flow rate through each slab at ib
         end do
         uouttot = sum(uin(kb+1:ke))/((zh(ke+1)-zh(kb+1))*jtot*dy)     ! convective outflow velocity
      end if

      ! Outlet
      ! Momentum
      v0(ie + 1, :, :) = v0(ie, :, :) - (v0(ie + 1, :, :) - v0(ie, :, :))*dxhi(ie + 1)*rk3coef*uouttot
      w0(ie + 1, :, :) = w0(ie, :, :) - (w0(ie + 1, :, :) - w0(ie, :, :))*dxhi(ie + 1)*rk3coef*uouttot
      vm(ie + 1, :, :) = vm(ie, :, :) - (vm(ie + 1, :, :) - vm(ie, :, :))*dxhi(ie + 1)*rk3coef*uouttot
      wm(ie + 1, :, :) = wm(ie, :, :) - (wm(ie + 1, :, :) - wm(ie, :, :))*dxhi(ie + 1)*rk3coef*uouttot
      e120(ie + 1, :, :) = e120(ie, :, :) - (e120(ie + 1, :, :) - e120(ie, :, :))*dxhi(ie + 1)*rk3coef*uouttot
      e12m(ie + 1, :, :) = e12m(ie, :, :) - (e12m(ie + 1, :, :) - e12m(ie, :, :))*dxhi(ie + 1)*rk3coef*uouttot

      ! Heat
      if (ltempeq) then
         thl0(ie + 1, :, :) = thl0(ie, :, :) - (thl0(ie + 1, :, :) - thl0(ie, :, :))*dxhi(ie + 1)*rk3coef*uouttot
         thlm(ie + 1, :, :) = thlm(ie, :, :) - (thlm(ie + 1, :, :) - thlm(ie, :, :))*dxhi(ie + 1)*rk3coef*uouttot
      end if

      if (lmoist) then
         qt0(ie + 1, :, :) = qt0(ie, :, :) - (qt0(ie + 1, :, :) - qt0(ie, :, :))*dxhi(ie + 1)*rk3coef*uouttot
         qtm(ie + 1, :, :) = qtm(ie, :, :) - (qtm(ie + 1, :, :) - qtm(ie, :, :))*dxhi(ie + 1)*rk3coef*uouttot
      end if

      ! tg3315 !changed dxhi to dxhci!?
      do n = 1, nsv
         sv0(ie + 1, :, :, n) = sv0(ie, :, :, n) - (sv0(ie + 1, :, :, n) - sv0(ie, :, :, n))*dxhci(ie + 1)*rk3coef*uouttot
         svm(ie + 1, :, :, n) = svm(ie, :, :, n) - (svm(ie + 1, :, :, n) - svm(ie, :, :, n))*dxhci(ie + 1)*rk3coef*uouttot
      end do

      return

   end subroutine iolet

   !>set boundary conditions pup,pvp,pwp in subroutine fillps in modpois.f90
   subroutine bcpup(pup, pvp, pwp, rk3coef)

      use modglobal, only:ib, ie, jb, je, ih, jh, kb, ke, kh, linoutflow, dxfi, iinletgen, &
         Uinf, libm, jmax, idriver
      use modfields, only:pres0, up, vp, wp, um, w0, u0, uouttot
      use modmpi, only:excjs, myid
      use modinletdata, only:irecy, u0inletbc, ddispdx, u0driver 

      real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh), intent(inout) :: pup
      real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh), intent(inout) :: pvp
      real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh), intent(inout) :: pwp

      real, intent(in) :: rk3coef
      real rk3coefi

      integer i, j, k

      rk3coefi = 1./rk3coef
      if (linoutflow) then
         if ((iinletgen == 1) .or. (iinletgen == 2)) then
            do j = jb, je
               do i = ib, ie
                  pwp(i, j, kb) = 0.
                  pwp(i, j, ke + kh) = (Uinf*ddispdx)*rk3coefi
               end do
            end do
            do k = kb, ke
               do j = jb, je
                  pup(ie + 1, j, k) = -(u0(ie + 1, j, k) - u0(ie, j, k))*dxfi(ie)*uouttot + um(ie + 1, j, k)*rk3coefi ! du/dt +u*du/dx=0 -> pup(i)=um(i)/rk3coef -um(i)*(um(i)-um(i-1))/dxf(i-1)
                  pup(ib, j, k) = u0inletbc(j, k)*rk3coefi
               end do
            end do
       elseif (idriver == 2) then
          do j=jb,je
             do i=ib,ie
                pwp(i,j,kb)  = 0.
                pwp(i,j,ke+kh)= 0. !(Uinf*ddispdx ) *rk3coefi ! tg3315 - idriver does not use Uinf ddisp etc.
             end do
          end do
          do k=kb,ke
             do j=jb,je
                pup(ie+1,j,k) = - (u0(ie+1,j,k)-u0(ie,j,k))*dxfi(ie)*uouttot + um(ie+1,j,k)*rk3coefi   ! du/dt +u*du/dx=0 -> pup(i)=um(i)/rk3coef -um(i)*(um(i)-um(i-1))/dxf(i-1)
                pup(ib,j,k) = u0driver(j,k)*rk3coefi
             end do
          end do 
         else ! if not iinletgen
            do j = jb, je
               do i = ib, ie
                  pwp(i, j, kb) = 0.
                  pwp(i, j, ke + kh) = 0.
               end do
            end do
            do k = kb, ke
               do j = jb, je
                  pup(ie + 1, j, k) = -(u0(ie + 1, j, k) - u0(ie, j, k))*dxfi(ie)*uouttot + um(ie + 1, j, k)*rk3coefi ! du/dt +u*du/dx=0 -> pup(i)=um(i)/rk3coef -um(i)*(um(i)-um(i-1))/dxf(i-1)
                  pup(ib, j, k) = pup(ib, j, k) - up(ib, j, k) ! pup(ib)= up(ib) + um(ib)/rk3coef, where up should be zero!
               end do
            end do
         end if ! inletgen
      else ! if not linoutflow
         do j = jb, je
            do i = ib, ie
               pwp(i, j, kb) = 0.
               pwp(i, j, ke + kh) = 0.
            end do
         end do
         do k = kb, ke
            do j = jb, je
               pup(ie + 1, j, k) = pup(ib, j, k) ! cyclic
               !pup(ib - 1, j, k) = pup(ie, j, k) ! tg3315 is this condition not needed? Was not here before but I think exists in Dales4.0 in modpois...!?
            end do
         end do
      endif

      call excjs(pup, ib, ie, jb, je, kb, ke + kh, ih, jh) ! cyclic
      call excjs(pvp, ib, ie, jb, je, kb, ke + kh, ih, jh) ! cyclic
      call excjs(pwp, ib, ie, jb, je, kb, ke + kh, ih, jh) ! cyclic

   end subroutine bcpup

   !>set pressure boundary conditions
   subroutine bcp(p)

      use modglobal, only:ib, ie, jb, je, ih, jh, kb, ke, kh, linoutflow, dxfi
      use modfields, only:pres0, up, u0, um, uouttot
      use modmpi, only:excj

      real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh), intent(inout) :: p !< pressure
      integer i, j, k

      if (linoutflow) then
         do k = kb, ke
            do j = jb, je
               p(ib - 1, j, k) = p(ib, j, k) ! inflow:  dp/dn=0
               pres0(ib - 1, j, k) = pres0(ib, j, k) ! inflow:  dp/dn=0
               p(ie + 1, j, k) = -p(ie, j, k) ! outflow: p=0
               pres0(ie + 1, j, k) = -pres0(ie, j, k) ! outflow: p=0
               up(ie + 1, j, k) = -(u0(ie + 1, j, k) - u0(ie, j, k))*dxfi(ie)*uouttot
            enddo
         enddo
      else
         do k = kb, ke
            do j = jb, je
               p(ib - 1, j, k) = p(ie, j, k)
               p(ie + 1, j, k) = p(ib, j, k)
            enddo
         enddo
      endif

      call excj(p, ib - 1, ie + 1, jb - 1, je + 1, kb - 1, ke + 1) ! cyclic
      call excj(pres0, ib - 1, ie + 1, jb - 1, je + 1, kb - 1, ke + 1) ! cyclic

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
      use modglobal, only:ke, kmax, lcoriol, igrw_damp, geodamptime
      use modfields, only:up, vp, wp, thlp, qtp, u0, v0, w0, thl0, qt0, ug, vg, thl0av, qt0av, u0av, v0av
      use modmpi, only:myid
      implicit none

      integer k
      select case (igrw_damp)
      case (0) !do nothing
      case (1)
         do k = ksp, ke
            up(:, :, k) = up(:, :, k) - (u0(:, :, k) - u0av(k))*tsc(k)
            vp(:, :, k) = vp(:, :, k) - (v0(:, :, k) - v0av(k))*tsc(k)
            wp(:, :, k) = wp(:, :, k) - w0(:, :, k)*tsc(k)
            thlp(:, :, k) = thlp(:, :, k) - (thl0(:, :, k) - thl0av(k))*tsc(k)
            qtp(:, :, k) = qtp(:, :, k) - (qt0(:, :, k) - qt0av(k))*tsc(k)
         end do
         if (lcoriol) then
            do k = ksp, ke
               up(:, :, k) = up(:, :, k) - (u0(:, :, k) - ug(k))*((1./(geodamptime*rnu0))*tsc(k))
               vp(:, :, k) = vp(:, :, k) - (v0(:, :, k) - vg(k))*((1./(geodamptime*rnu0))*tsc(k))
            end do
         end if
      case (2)
         do k = ksp, ke
            up(:, :, k) = up(:, :, k) - (u0(:, :, k) - ug(k))*tsc(k)
            vp(:, :, k) = vp(:, :, k) - (v0(:, :, k) - vg(k))*tsc(k)
            wp(:, :, k) = wp(:, :, k) - w0(:, :, k)*tsc(k)
            thlp(:, :, k) = thlp(:, :, k) - (thl0(:, :, k) - thl0av(k))*tsc(k)
            qtp(:, :, k) = qtp(:, :, k) - (qt0(:, :, k) - qt0av(k))*tsc(k)
         end do
      case (3)
         do k = ksp, ke
            up(:, :, k) = up(:, :, k) - (u0(:, :, k) - u0av(k))*tsc(k)
            vp(:, :, k) = vp(:, :, k) - (v0(:, :, k) - v0av(k))*tsc(k)
            wp(:, :, k) = wp(:, :, k) - w0(:, :, k)*tsc(k)
            thlp(:, :, k) = thlp(:, :, k) - (thl0(:, :, k) - thl0av(k))*tsc(k)
            qtp(:, :, k) = qtp(:, :, k) - (qt0(:, :, k) - qt0av(k))*tsc(k)
         end do
      case default
         write(0, *) "ERROR: no gravity wave damping option selected"
         stop 1
      end select

      return
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
      use modglobal, only:ib, ie, ih, jb, je, jh, kb, ke, kh, eps1
      real, intent(inout) :: field(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh)
      real, intent(in)    :: val
      !
      field(:, :, ke + 1) = 2*val - field(:, :, ke)
      !
   end subroutine valuetop

   subroutine fluxtopscal(flux)
      use modglobal, only:ib, ie, ih, jb, je, jh, kb, ke, kh, dzf, dzh, dzhi, nsv, khc
      use modfields, only:sv0, svm
      use modsubgriddata, only:ekh

      real, intent(in)    :: flux(1:nsv)
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
      use modglobal, only:ib, ie, ih, jb, je, jh, kb, ke, kh, eps1, nsv, khc
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

   !> Sets top boundary conditions for momentum
   subroutine inlettop

      use modglobal, only:ib, ie, jb, je, kb, ke, ih, jh, kh, dzh, dzf, &
         e12min, dxfi, dxf, dxhi, xh, jgb, jge, Uinf, dzfi
      use modfields, only:w0, wm, wout, wouttot
      use modinletdata, only:Uinl, ddispdxold
      use modmpi, only:slabsumi, myid
      implicit none
      integer :: i
      real    :: nji

      do i = ib, ie
         w0(i, :, ke + 1) = Uinf*ddispdxold
         wm(i, :, ke + 1) = Uinf*ddispdxold
      end do
      call slabsumi(wout, ib, ie, w0, ib - ih, ie + ih, jb - jh, je + jh, kb - kh, ke + kh, ib, ie, jb, je, ke + 1, ke + 1) ! determine vertical (j) average outflow velocity
      nji = 1./(jge - jgb + 1)
      do i = ib, ie
         wout(i) = wout(i)*dxf(i)*nji
      end do
      wouttot = sum(wout(ib:ie))/(xh(ie + 1) - xh(ib)) ! Area-averaged outflow velocity

      return
   end subroutine inlettop

!>Set thl, qt and sv(n) equal to slab average at level kmax
   subroutine tqaver

      use modmpi, only:comm3d, mpierr, my_real, mpi_sum
      use modglobal, only:ib, ie, jb, je, ih, jh, kb, ke, nsv, rslabs
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
