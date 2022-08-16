!> \file modpois.f90
!!  Solves the Poisson equation for the pressure fluctuations
!>
!!  \author Jasper Tomas, TU Delft, 31 March 2014
!!  \author Harm Jonker, TU Delft
!!  \author Hans Cuijpers, IMAU
!!  \todo documentation
!!  \par Revision list
!! Jasper Tomas: Now a different Poisson solver is implemented that can be used for inflow/outflow boundary conditions in the i-direction.
!! Periodic BC's can also still be used.
!
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

module modpois
use decomp_2d
use decomp_2d_fft

implicit none

include "fftw3.f"

private
public :: initpois,poisson,exitpois,p,pup,pvp,pwp,rhs,dpupdx,dpvpdy,dpwpdz,xyzrt
save

  real, allocatable, target :: p(:,:,:)   ! difference between p at previous and new step (p = P_new - P_old)
  real, allocatable, target :: pup(:,:,:), pvp(:,:,:), pwp(:,:,:)
  real, allocatable, target :: rhs(:,:,:)   ! rhs of pressure solver
  real, allocatable, target :: dpupdx(:,:,:)
  real, allocatable, target :: dpvpdy(:,:,:)
  real, allocatable, target :: dpwpdz(:,:,:)
  real, allocatable :: xrt(:), yrt(:)
  real, allocatable, target :: xyzrt(:,:,:)
  real, allocatable :: a(:), b(:), c(:), ax(:), bx(:), cx(:) ! coefficients for tridiagonal matrix
  integer :: ibc1, ibc2, kbc1, kbc2

  !integer*8 :: plan_r2fc_x, plan_r2fc_y, plan_fc2r_x, plan_fc2r_y
  integer*8 :: plan_r2fr_x, plan_r2fr_y, plan_fr2r_x, plan_fr2r_y
  real, allocatable :: Sxr(:), Sxfr(:)
  real, allocatable :: Syr(:), Syfr(:)
  ! complex, dimension(0:itot/2):: Sxfc
  ! complex, dimension(0:jtot/2):: Syfc

contains
  subroutine initpois
    use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh,imax,jmax,itot,jtot,ktot,dyi,dxfi,ipoiss,POISS_FFT2D,POISS_FFT3D,POISS_CYC,pi,dy,linoutflow,dzh,dzf,dxh,dxf
    use modmpi,    only : myidx, myidy
    use modfields, only : rhobf, rhobh
    !use decomp_2d
    !use decomp_2d_fft
    implicit none
    real dxi, fac
    integer i,j,k,iv,jv

    !allocate(p(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(pup(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(pvp(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(pwp(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    call alloc_z(p)
    call alloc_z(rhs, opt_zlevel=(/0,0,0/))
    call alloc_z(dpupdx, opt_zlevel=(/0,0,0/))
    call alloc_z(dpvpdy, opt_zlevel=(/0,0,0/))
    call alloc_z(dpwpdz, opt_zlevel=(/0,0,0/))


    if (ipoiss == POISS_FFT2D) then
      dxi = dxfi(1) ! Assumes equidistant in x
      allocate(xrt(itot))
      allocate(yrt(jtot))
      allocate(xyzrt(imax,jmax,ktot))
      allocate(a(ktot), b(ktot), c(ktot))

      allocate(Sxr(0:itot-1), Sxfr(0:itot-1))
      allocate(Syr(0:jtot-1), Syfr(0:jtot-1))

      kbc1 = 1 ! Neumann

      ! generate Eigenvalues xrt and yrt
      if (.not. linoutflow) then
        ! streamwise - periodic
        fac = 1./(2.*itot)
        do i=3,itot,2
          xrt(i-1)=-4.*dxi*dxi*(sin(float((i-1))*pi*fac))**2
          xrt(i)  = xrt(i-1)
        end do
        xrt(1    ) = 0.
        xrt(itot ) = -4.*dxi*dxi

        fac = 1./(2.*jtot)
        do j=3,jtot,2
          yrt(j-1)=-4.*dyi*dyi*(sin(float((j-1))*pi*fac))**2
          yrt(j  )= yrt(j-1)
        end do
        yrt(1    ) = 0.
        yrt(jtot ) = -4.*dyi*dyi

        ! top - Neumann
        kbc2 = 1

      else
        ! streamise - Neumann
        fac = 1./(2.*itot)
        do i=1,itot
          xrt(i)=-4.*dxi*dxi*(sin(float((i-1))*pi*fac))**2
          !xrt(i)=-4.*dxi*dxi*(sin(float((i))*pi*fac))**2
        end do

        fac = 1./(2.*jtot)
        do j=1,jtot
          yrt(j)=-4.*dyi*dyi*(sin(float((j-1))*pi*fac))**2
          !yrt(j)=-4.*dyi*dyi*(sin(float((j))*pi*fac))**2
        end do

        ! top - Dirichlet
        !kbc2 = 3
        kbc2 = 1
        ! Could also set a Dirichlet BC on the inlet at the top, like in SPARKLE?

        call dfftw_plan_r2r_1d(plan_r2fr_x,itot,Sxr,Sxfr,FFTW_REDFT10,FFTW_MEASURE)
        call dfftw_plan_r2r_1d(plan_r2fr_y,jtot,Syr,Syfr,FFTW_REDFT10,FFTW_MEASURE)
        call dfftw_plan_r2r_1d(plan_fr2r_x,itot,Sxfr,Sxr,FFTW_REDFT01,FFTW_MEASURE)
        call dfftw_plan_r2r_1d(plan_fr2r_y,jtot,Syfr,Syr,FFTW_REDFT01,FFTW_MEASURE)

      end if

      ! ! spanwise - periodic
      ! fac = 1./(2.*jtot)
      ! do j=3,jtot,2
      !   yrt(j-1)=-4.*dyi*dyi*(sin(float((j-1))*pi*fac))**2
      !   yrt(j  )= yrt(j-1)
      ! end do
      ! yrt(1    ) = 0.
      ! yrt(jtot ) = -4.*dyi*dyi

      !xyrt = 0.
      ! SO: rhobf is not here in DALES 4.0 - fine because it is set to 1 currently but check.

      do k=1,zsize(3)
        do j=1,zsize(2)
          jv = j + myidy*zsize(2)!jmax
          do i=1,zsize(1)
            iv = i + myidx*zsize(1)!imax
            xyzrt(i,j,k) = rhobf(k)*(xrt(iv)+yrt(jv))
          end do
        end do
      end do

      ! generate matrix coefficients
      do k=1,ktot
        a(k) = rhobh(k)  /(dzf(k)*dzh(k))
        c(k) = rhobh(k+1)/(dzf(k)*dzh(k+1))
        b(k) = -(a(k)+c(k))
      end do

      ! bottom - Neumann
      !b(1) = b(1) + a(1)
      b(1) = -c(1)

      ! top
      if (kbc2 .eq. 1) then ! Neumann
        !b(ktot) = b(ktot) + c(ktot)
        b(ktot) = -a(ktot)
      elseif (kbc2 .eq. 3) then ! Dirichlet
        b(ktot) = b(ktot) - c(ktot) ! not convinced this is right, but performs better than below
        !b(ktot) = a(ktot)
      end if

      a(1) = 0.
      c(ktot) = 0.

    elseif (ipoiss == POISS_FFT3D) then
      !call decomp_2d_fft_init
      ! Define 3D FFT coefficients

    elseif (ipoiss == POISS_CYC) then

      allocate(yrt(jtot))
      allocate(ax(itot), bx(itot), cx(itot))
      allocate(a(ktot), b(ktot), c(ktot))

      ! spanwise - periodic
      yrt(1) = 0.
      yrt(jtot)=-4./(dy*dy)
      do j=3,jtot,2
        yrt(j-1)=(-4./(dy*dy))*(sin(float((j-1))*pi/(2.*jtot)))**2
        yrt(j)= yrt(j-1)
      end do

      if (linoutflow ) then
        ibc1 = 1 ! Neumann
        ibc2 = 3 ! Dirichlet
      else
        ibc1 = 2
        ibc2 = 2
      endif
      kbc1 = 1
      kbc2 = 1

      do i=1,itot
        ax(i) =  1./(dxf(i)*dxh(i))
        cx(i) =  1./(dxf(i)*dxh(i+1))
        bx(i) =  - (ax(i) + cx(i))
      end do

      if (ibc1 .eq.1) then
        !bx(1) = bx(1) + ax(1)
        bx(1) = -cx(1)
      !elseif (ibc1 .eq. 2) then
        !bx(1) = bx(1)
      elseif (ibc1 .eq. 3) then
        bx(1) = bx(1) - ax(1) ! not convinced this is right
        !bx(1) = cx(1)
      end if

      if (ibc2 .eq. 1) then
        !bx(itot) = bx(itot) + cx(itot)
        bx(itot) = -ax(1)
      !elseif (ibc2 .eq. 2) then
        !bx(itot) = bx(itot)
      elseif (ibc2 .eq. 3) then
        !bx(itot) = bx(itot) - cx(itot)
        bx(itot) = ax(itot)
      end if

      if (ibc1 .ne. 2) then
        ax(1)    = 0.
      end if

      if (ibc2 .ne. 2) then
        cx(itot) = 0.
      end if

    end if ! ipoiss

  end subroutine initpois

  subroutine poisson
    use modglobal, only : ib,ie,ih,kb,ke,kh,ktot,dxh,dxf,dy,dzf,dzh,linoutflow,iinletgen,ipoiss,POISS_FFT2D,POISS_FFT3D,POISS_CYC
    use modmpi, only : myid,nprocs,barrou
    implicit none
    integer ibc1,ibc2,kbc1,kbc2,ksen

    call fillps

!  ibc?=1: neumann
!  ibc?=2: periodic
!  ibc?=3: dirichlet

    select case (ipoiss)
    case (POISS_FFT2D)
      call solmpj
    case (POISS_FFT3D)

    case (POISS_CYC)
      ! if (linoutflow ) then
      !   ibc1 = 1      ! inlet
      !   ibc2 = 3      ! outlet
      ! else
      !   ibc1 = 2
      !   ibc2 = 2
      ! endif
      ! kbc1 = 1
      ! kbc2 = 1
      !ksen = ktot/nprocs ! This needs changing - nprocs no longer defined
      !call poisr(p,dxf,dxh,dy,dzf,dzh, &
                 !ibc1,ibc2,kbc1,kbc2,ksen)
      call poisr
    case default
      write(0, *) "Invalid choice for Poisson solver"
       stop 1
    end select

    call tderive

  end subroutine poisson

  subroutine exitpois
    implicit none
    deallocate(p, xrt, yrt, xyzrt, a, b, c)
  end subroutine exitpois
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fillps

  ! Chiel van Heerwaarden,  19 June 2007
  ! Adapted fillps for RK3 time loop


    use modfields, only : up, vp, wp, um, vm, wm,u0,v0,w0
    use modglobal, only : rk3step, ib,ie,jb,je,kb,ke,ih,jh,kh, dxfi,dyi,dzfi,dt,&
                          linoutflow,libm,dtmax,ierank,jerank,pi,dy,imax,jmax,ylen,xf,zf
    use modmpi,    only : excjs, myidx, myidy
    use modboundary, only: bcpup
!    use modibm,    only : ibmnorm

    implicit none
    !real,allocatable :: pup(:,:,:), pvp(:,:,:), pwp(:,:,:)
    integer i,j,k
    real rk3coef,rk3coefi

    ! allocate(pup(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    ! allocate(pvp(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    ! allocate(pwp(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))

    if (rk3step == 0) then ! dt not defined yet
      rk3coef = 1.
    else
      rk3coef = dt / (4. - dble(rk3step))
    end if
    rk3coefi = 1. / rk3coef

    do k=kb,ke
      do j=jb,je
        do i=ib,ie
          pup(i,j,k) = up(i,j,k) + um(i,j,k) * rk3coefi         ! see equation 5.81
          pvp(i,j,k) = vp(i,j,k) + vm(i,j,k) * rk3coefi
          pwp(i,j,k) = wp(i,j,k) + wm(i,j,k) * rk3coefi
        end do
      end do
    end do

  !****************************************************************

  !     Fill the right hand for the poisson solver.
  !     The values for up(i2,j,k) and vp(i,j2,k) are still
  !     unknown and have to be set cyclic.
  !     Also we take wp(i,j,1) and wp(i,j,k1) equal to zero.

  !     NOTE:

  !     The poisson-solver only accepts values for i from 2 to i1,
  !     for j from 1 to jmax and for k from 1 to kmax.
  !     The right-hand p is therefore filled in this partical way.

  !**************************************************************

   call bcpup(pup,pvp,pwp,rk3coef)   ! boundary conditions for pup,pvp,pwp

    do k=kb,ke
      do j=jb,je
        do i=ib,ie
          dpupdx(i,j,k) = (pup(i+1,j,k)-pup(i,j,k)) * dxfi(i)
          dpvpdy(i,j,k) = (pvp(i,j+1,k)-pvp(i,j,k)) * dyi
          dpwpdz(i,j,k) = (pwp(i,j,k+1)-pwp(i,j,k)) * dzfi(k)

          p(i,j,k)  =  ( pup(i+1,j,k)-pup(i,j,k) ) * dxfi(i) &         ! see equation 5.72
                          +( pvp(i,j+1,k)-pvp(i,j,k) ) * dyi &
                          +( pwp(i,j,k+1)-pwp(i,j,k) ) * dzfi(k)
        end do
      end do
    end do

    ! ! MOMS
    ! do k=kb,ke
    !   do j=jb,je
    !     do i=ib,ie
    !       p(i,j,k) = -12 * pi**2 / ylen**2 * cos(2. * pi * xf(i+myidx*imax) / ylen) &
    !                                       * cos(2. * pi * (dy*((0.5+(j-1))+myidy*jmax)) / ylen) &
    !                                       * cos(2. * pi * (zf(k)) / ylen)
    !     end do
    !   end do
    ! end do

    !rhs = p(ib:ie,jb:je,kb:ke)

    !deallocate( pup,pvp,pwp )

  end subroutine fillps


  subroutine tderive

    !-----------------------------------------------------------------|
    !                                                                 |
    !*** *tderive*  read input fields for initialisation              |
    !                                                                 |
    !      Hans Cuijpers   I.M.A.U.     06/01/1995                    |
    !                                                                 |
    !     purpose.                                                    |
    !     --------                                                    |
    !                                                                 |
    !     Refill array p with pressure values. The poisson-solver     |
    !     produced a pressure array p in which the i-index varied     |
    !     between 2 and i1, the j- and k-index between 1 and resp.    |
    !     jmax and kmax. For our further calculations we'll change    |
    !     the range for the j-index to vary between 2 and j1.         |
    !                                                                 |
    !     Further we set cyclic boundary conditions for the pressure- |
    !     fluctuations in the x-y plane.                              |
    !                                                                 |
    !**   interface.                                                  |
    !     ----------                                                  |
    !                                                                 |
    !             *tderive* is called from *program*.                 |
    !                                                                 |
    !-----------------------------------------------------------------|

    use modfields, only : u0, v0, w0, up, vp, wp, pres0, IIc, IIcs, uouttot
    use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh,dxhi,dyi,dzhi,linoutflow,rslabs,ibrank,ierank,jbrank,jerank,dxfi
    use modmpi,    only : myid,excj,slabsum,avexy_ibm
    use modboundary,only : bcp
    implicit none
    integer i,j,k
    real, dimension(kb-kh:ke+kh) :: pij
    !    logical, dimension(ib:ie, jb:je, kb:ke) :: pnan
    real :: pijk
    !    integer :: ipnan

    ! Mathieu ATTTT: CHANGED!!! Loop removed!!!

    ! **  Boundary conditions **************

    call bcp(p) ! boundary conditions for p.

    ! write(*,*) "uouttot", uouttot
    !
    ! call exchange_halo_z(p)
    ! call exchange_halo_z(pres0)
    !
    ! if (linoutflow) then
    !   if (ibrank) then
    !     do k = kb, ke
    !       do j = jb-1, je+1
    !         p(ib - 1, j, k) = p(ib, j, k) ! inflow:  dp/dn=0
    !         pres0(ib - 1, j, k) = pres0(ib, j, k) ! inflow:  dp/dn=0
    !       enddo
    !     enddo
    !   end if
    !
    !   if (ierank) then
    !     do k = kb, ke
    !       do j = jb-1, je+1
    !         ! p(ie + 1, j, k) = -p(ie, j, k) ! outflow: p=0
    !         ! pres0(ie + 1, j, k) = -pres0(ie, j, k) ! outflow: p=0
    !         p(ie + 1, j, k) = p(ie, j, k) ! outflow: dp/dn=0
    !         pres0(ie + 1, j, k) = pres0(ie, j, k) ! outflow: dp/dn=0
    !         up(ie + 1, j, k) = -(u0(ie + 1, j, k) - u0(ie, j, k))*dxfi(ie)*uouttot
    !       enddo
    !     enddo
    !   end if
    ! endif

    !*****************************************************************
    ! **  Calculate time-derivative for the velocities with known ****
    ! **  pressure gradients.  ***************************************
    !*****************************************************************

    !   if (myid == 0) then
    !     write(*,*) "net ts", p(ib, jb, :)
    !   end if

    !pnan = isnan(p(ib:ie, jb:je, kb:ke))
    !ipnan = 0
    !do k=kb,ke
    !do j=jb,je
    !do i=ib,ie
    !  if (pnan(i,j,k)) then
    !    ipnan = ipnan + 1
    !  end if
    !end do
    !end do
    !end do
    !   write(*,*) "NaN", myid, ipnan

    !    write(*,*) "NaN", myid, dble(isnan(p)) !sum(real(isnan(p)))

    ! do k=kb,ke
    !   do j=jb,je
    !     do i=ib,ie
    !       vp(i,j,k) = vp(i,j,k)-(p(i,j,k)-p(i,j-1,k))*dyi
    !     end do
    !   end do
    ! end do
    !
    ! if (linoutflow .and. ierank) then
    !   do k=kb,ke
    !     do j=jb,je
    !       do i=ib,ie+1
    !         up(i,j,k) = up(i,j,k)-(p(i,j,k)-p(i-1,j,k))*dxhi(i)           ! see equation 5.82 (u is computed from the mass conservation)
    !       end do
    !     end do
    !   end do
    ! else
    !   do k=kb,ke
    !     do j=jb,je
    !       do i=ib,ie
    !         up(i,j,k) = up(i,j,k)-(p(i,j,k)-p(i-1,j,k))*dxhi(i)           ! see equation 5.82 (u is computed from the mass conservation)
    !       end do
    !     end do
    !   end do
    ! endif
    !
    ! do k=kb+1,ke
    !   do j=jb,je
    !     do i=ib,ie
    !       wp(i,j,k) = wp(i,j,k)-(p(i,j,k)-p(i,j,k-1))*dzhi(k)
    !     end do
    !   end do
    ! end do

    do i=ib,ie
      do j=jb,je
        up(i,j,kb) = up(i,j,kb)-(p(i,j,kb)-p(i-1,j,kb))*dxhi(i)
        vp(i,j,kb) = vp(i,j,kb)-(p(i,j,kb)-p(i,j-1,kb))*dyi
        do k=kb+1,ke
          up(i,j,k) = up(i,j,k)-(p(i,j,k)-p(i-1,j,k))*dxhi(i)
          vp(i,j,k) = vp(i,j,k)-(p(i,j,k)-p(i,j-1,k))*dyi
          wp(i,j,k) = wp(i,j,k)-(p(i,j,k)-p(i,j,k-1))*dzhi(k)
        end do
      end do
    end do

    if (linoutflow) then
      if (ierank) then
        do j=jb,je
          do k=kb,ke
            up(ie+1,j,k) = up(ie+1,j,k)-(p(ie+1,j,k)-p(ie,j,k))*dxhi(i)
          end do
        end do
      end if

      if (jerank) then
        do i=ib,ie
          do k=kb,ke
            vp(i,je+1,k) = vp(i,je+1,k)-(p(i,je+1,k)-p(i,je,k))*dyi
          end do
        end do
      end if
    end if


    ! tg3315 02/02/2019
    ! account for pressure offset that results from ill-defined problem in pressure
    ! correction method when periodic horizontal BCs are applied. Arises within cyclic
    ! reduction scheme (called by BLKTRI) and due to the numerics in PRODP in
    ! cycred.f which define the BC in the periodic case. A linear offset existed in
    ! the pressure correction term (p) and this can be controlled by subtracting
    ! the volume averaged modified pressure from this value at all time steps.
    ! Periodic: p - <p>_ijk
    ! Makes no change on physical effect of modified pressure in code.

    ! tg3315 - update 24/06/19 -- there is a missing term in the application of the periodic BCs for pup, could this be part of the problem? Test with this to see if can avoid use of pijk below.
    ! refer to mvr for necessity of this

    ! useful refs:
    ! https://opensky.ucar.edu/islandora/object/technotes%3A98/datastream/PDF/download/citation.pdf
    ! https://epubs.siam.org/doi/pdf/10.1137/0711042

    pij =0.; pijk=0.;

    if (.not. linoutflow) then
      call slabsum(pij(kb:ke),kb,ke,p(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,ib,ie,jb,je,kb,ke)
      pij = pij/rslabs
      pijk = sum(pij(kb:ke))/(ke-kb)
    end if

    !write(*,*) pijk

    do k=kb-1,ke+1
      do j=jb-1,je+1
        do i=ib-1,ie+1
          pres0(i,j,k)=pres0(i,j,k)+p(i,j,k)-pijk ! update of the pressure: P_new = P_old + p
        enddo
      enddo
    enddo


    ! if (linoutflow) then
    !   if (ibrank) then
    !     do k = kb, ke
    !       do j = jb-1, je+1
    !         pres0(ib - 1, j, k) = pres0(ib, j, k) ! inflow:  dp/dn=0
    !       enddo
    !     enddo
    !   end if
    !
    !   if (ierank) then
    !     do k = kb, ke
    !       do j = jb-1, je+1
    !         pres0(ie + 1, j, k) = pres0(ie, j, k) ! outflow: dp/dn=0
    !       enddo
    !     enddo
    !   end if
    !
    !   if (jbrank) then
    !     do k = kb, ke
    !       do i = ib-1, ie+1
    !         pres0(i,jb-1,k) = pres0(i,jb,k) ! inflow:  dp/dn=0
    !       enddo
    !     enddo
    !   end if
    !
    !   if (jerank) then
    !     do k = kb, ke
    !       do i = ib-1, ie+1
    !         pres0(i,je+1,k) = pres0(i,je,k) ! outflow: dp/dn=0
    !       enddo
    !     enddo
    !   end if
    ! end if

    return
  end subroutine tderive

!ils13, 13.08.18: currently unused, not callled
! SO: replaced by 2DECOMP routines
  subroutine ALL_ALL_j(p,ptrans,iaction)
! purpose: do all-to-all communication
! data are only distributed over the j-direction for p
! data are only distributed over the k-direction for ptrans
! NOTE: p     (0:imax+1  etc
!       ptrans(1:imax    etc


  use modglobal, only : imax,isen,jmax,jsen,jtot,kmax
  use modmpi,    only : comm3d,mpierr,my_real,nprocs, barrou

  implicit none


  integer  iaction
  real     p(0:imax+1,0:jmax+1,0:kmax+1)
  real     ptrans(1:isen,1:jtot,1:kmax)


! help arrays for sending and receiving

  real,allocatable,dimension(:) :: bufin, bufout


! help variables

  integer ii, jbegin, jend, proc
  integer     ibegin, iend
  integer     i, j, k

  allocate(bufin(imax*jmax*kmax),bufout(imax*jmax*kmax))
  if(iaction==0)then
    ii = 0
    do proc=0,nprocs-1
      ibegin =  (proc)*isen+1
      iend   =  (proc+1)*isen
      do i=ibegin,iend
      do j=1,jmax
      do k=1,kmax
        ii = ii + 1
        bufin(ii) = p(i,j,k)
      enddo
      enddo
      enddo
    enddo

    ii = 0
!     call barrou()
    call MPI_ALLTOALL(bufin,   (isen*jsen*kmax),MY_REAL, &
                          bufout,(isen*jsen*kmax),MY_REAL, &
                          comm3d,mpierr)
!     call barrou()
    ii = 0
    do proc = 0,nprocs-1
      jbegin =  proc   *jsen + 1
      jend   = (proc+1)*jsen
      do i=1,isen
      do j=jbegin,jend
      do k=1,kmax
        ii = ii + 1
        ptrans(i,j,k) = bufout(ii)
      enddo
      enddo
      enddo
    enddo
!     call barrou()

  elseif(iaction==1)then
    ii = 0
    do  proc = 0,nprocs-1
      jbegin =  proc   *jsen + 1
      jend   = (proc+1)*jsen
      do i=1,isen
      do j=jbegin,jend
      do k=1,kmax
        ii = ii + 1
        bufin(ii)  = ptrans(i,j,k)
      enddo
      enddo
      enddo
    enddo

!     call barrou()
    call MPI_ALLTOALL(bufin,   (isen*jsen*kmax),MY_REAL, &
                          bufout,(isen*jsen*kmax),MY_REAL, &
                          comm3d,mpierr)
!     call barrou()

    ii = 0
    do proc=0,nprocs-1
      ibegin =    (proc)*isen+1
      iend   =  (proc+1)*isen
      do i=ibegin,iend
      do j=1,jmax
      do k=1,kmax
        ii = ii + 1
        p(i,j,k) = bufout(ii)
      enddo
      enddo
      enddo
    enddo
!     call barrou()
  endif

  deallocate(bufin,bufout)

  return
  end subroutine ALL_ALL_j
!
! Mathieu added ALL_ALL_j2 which transposes between j and k instead of i and k
! This was to be able to use solver poisr, maybe we should tidy this up later
! on! ALL_ALL_j2 uses ksen which is NOT in the modules... it is provided in the
! header for this reason
!

      subroutine ALL_ALL_j2(p,ptrans,iaction,ksen)
!
  use modglobal, only : imax,isen,jmax,jsen,jtot,kmax
  use modmpi,    only : comm3d,mpierr,my_real,nprocs, barrou

      implicit none
!
!     include 'param.txt'
!
!     include 'mpif.h'
!     include 'mpi_cons.txt'
!
      integer  iaction, ksen
      real     p(0:imax+1,0:jmax+1,0:kmax+1)
      real     ptrans(1:imax,1:jtot,1:ksen)
      integer i,j,k
!
!
! purpose: do all-to-all communication
!
! data are only distributed over the j-direction for p
! data are only distributed over the k-direction for ptrans
!
! help arrays for sending and receiving
!
       real bufin ((imax)*jmax*(kmax))
       real bufout((imax)*jmax*(kmax))
!
! help variables
!
       integer ii, jstart, jend, proc
       integer     kstart, kend, jvalue, kvalue
!
!
       if(iaction.eq.0)then
!
       ii = 0
       do proc=0,nprocs-1
          kstart =    (proc)*ksen+1
          kend   =  (proc+1)*ksen
       do k=kstart,kend
       do j=1,jmax
       do i=1,imax
          ii = ii + 1
          bufin(ii) = p(i,j,k)
       enddo
       enddo
       enddo
       enddo

!
!
       ii = 0
!
!
       call barrou()
!
       call MPI_ALLTOALL(bufin,   (imax*jsen*ksen),MY_REAL, &
                           bufout,(imax*jsen*ksen),MY_REAL, &
                           comm3d,mpierr)
!      bufout = bufin
!
       call barrou()
!
       ii = 0
!
       do  proc = 0,nprocs-1
           jstart =  proc   *jsen + 1
           jend   = (proc+1)*jsen
       do k=1,ksen
       do j=jstart,jend
       do i=1,imax
          ii = ii + 1
          ptrans(i,j,k) = bufout(ii)
       enddo
       enddo
       enddo
!
       enddo
!
!
       call barrou()
!
!
       elseif(iaction.eq.1)then
!
       ii = 0
!
       do  proc = 0,nprocs-1
           jstart =  proc   *jsen + 1
           jend   = (proc+1)*jsen
       do k=1,ksen
       do j=jstart,jend
       do i=1,imax
          ii = ii + 1
          bufin(ii)  = ptrans(i,j,k)
       enddo
       enddo
       enddo
!
       enddo
!
       call barrou()
!
       call MPI_ALLTOALL(bufin,   (imax*jsen*ksen),MY_REAL, &
                           bufout,(imax*jsen*ksen),MY_REAL, &
                           comm3d,mpierr)
!      bufout = bufin
!
       call barrou()
!
!
       ii = 0
!
       do proc=0,nprocs-1
          kstart =    (proc)*ksen+1
          kend   =  (proc+1)*ksen
       do k=kstart,kend
       do j=1,jmax
       do i=1,imax
          ii = ii + 1
          p(i,j,k) = bufout(ii)
       enddo
       enddo
       enddo
       enddo
!
!
       call barrou()
       endif
!
       return
       end subroutine ALL_ALL_j2

  subroutine poisr
  !subroutine poisr(rhs,dx,dxh,dy,dz,dzh, &
                       !ibc1,ibc2,kbc1,kbc2,ksen)
  !
  ! CHANGES:
  ! includes jtot and ksen (calculated in poisson.f) in header
  !          help array rhst with dimensions (imax,jtot,ksen)
  !          ALL_ALL copy rhs to rhst and back for FFT's
  !
  !  ibc?=1: neumann
  !  ibc?=2: periodic
  !  ibc?=3: dirichlet
  !
  ! only FFT in j-direction, cyclic reduction in the others

  !      include'param.txt'
  !     include 'mpif.h'
  !     include 'mpi_cons.txt'
    use modglobal, only : imax,itot,jmax,jtot,kmax,ktot,poisrcheck,imax2,jmax2,kmax2,ib,ie,jb,je,kb,ke
    use modfields, only : worksave
    use modmpi,    only : myid,comm3d,mpierr,my_real,nprocs, barrou,MPI_SUM
    implicit none

    ! authors: m.j.b.m. pourquie, b.j. boersma

    integer i, j, k, ksen
    real rhs(0:imax+1,0:jmax+1,0:kmax+1)!,dx(0:IMAX+1),dxh(1:IMAX+1),dy
    real, allocatable, dimension(:,:,:) :: py, pz
    real, allocatable, dimension(:,:,:) ::  rhs2
    real, allocatable, dimension(:) ::  work
    integer  ier, iperio, kperio
    integer  ibc1,ibc2,kbc1,kbc2
    !real dz(0:kmax+1),dzh(1:kmax+1),pi       !   dz: kb-1:ke+1,  dzh: kb:ke+1
    !real a(imax),b(imax),c(imax),bin(imax)
    !real ax(itot),bx(itot),cx(itot)
    real bin(itot)
    !real az(kmax),bz(kmax),cz(kmax)
    !real az(ktot),bz(ktot),cz(ktot)
    !real yrt(jtot)
    real, allocatable, dimension(:,:) ::  vfftj
    real, allocatable, dimension(:,:) ::  FFTJ
    real, allocatable, dimension(:,:) ::  y
    real wj(jtot+15)
    real  angle, tst
    real suml, sum
    integer ipos, jv

    !allocate(rhs2(imax,jtot,ksen))
    !allocate(work(2*imax*jmax*kmax))
    !allocate(vfftj(imax*ksen,jtot))
    allocate(FFTJ(imax2*kmax2,jtot))
    !allocate(y(imax,kmax))

    call alloc_y(py, opt_ylevel=(/0,0,0/))
    call alloc_z(pz, opt_zlevel=(/0,0,0/))
    pz = p(ib:ie,jb:je,kb:ke)

    !pi=4.*atan(1.)

    ! !do i=1,imax
    ! do i=1,itot
    !   ax(i) =  1./(dx(i)*dxh(i))
    !   cx(i) =  1./(dx(i)*dxh(i+1))
    !   bx(i) =  - (ax(i) + cx(i))
    ! end do
    !
    ! if ((ibc1) .eq.1) then
    !   bx(1) = bx(1) + ax(1)
    ! elseif (ibc1 .eq. 2) then
    !   bx(1) = bx(1)
    ! elseif (ibc1 .eq. 3) then
    !   bx(1) = bx(1) - ax(1)
    ! end if
    !
    ! if ((ibc2).eq.1) then
    !   bx(itot) = bx(itot) + cx(itot)
    ! elseif ((ibc2).eq.2) then
    !   bx(itot) = bx(itot)
    ! elseif ((ibc2).eq.3) then
    !   bx(itot) = bx(itot) - cx(itot)
    ! end if
    !
    ! if (ibc1 .ne. 2) then
    !   cx(itot) = 0.
    !   ax(1)    = 0.
    ! end if

    ! fill coefficients in k-direction

    ! do k=1,kmax                          ! Mathieu's version
    !   az(k) =  1./(dz(k)*dzh(k-1))
    !   cz(k) =  1./(dz(k)*dzh(k))
    !   bz(k) =  - (az(k) + cz(k))
    ! enddo

    !   !do k=1,kmax
    ! do k=1,ktot
    !   az(k) = 1./(dz(k) * dzh(k))
    !   cz(k) = 1./(dz(k) * dzh(k+1))
    !   bz(k) = -(az(k) + cz(k))
    ! end do
    !
    ! if((kbc1).eq.1)then
    !   bz(1) = bz(1) + az(1)
    ! elseif(kbc1.eq.2)then
    !   bz(1) = bz(1)
    ! end if
    !
    ! if ((kbc2) .eq. 1) then
    !   bz(ktot) = bz(ktot) + cz(ktot)
    ! elseif ((kbc2) .eq. 2) then
    !   bz(ktot) = bz(ktot)
    ! elseif ((kbc2) .eq. 3) then
    !   bz(ktot) = bz(ktot)  - cz(ktot)
    ! end if
    !
    ! if (kbc1.ne.2) then
    !   cz(ktot) = 0.
    !   az(1)    = 0.
    ! end if

    ! initialise for FFT
    call vrffti(jtot,wj)

    ! yrt(1) = 0.
    ! yrt(jtot)=-4./(dy*dy)
    ! do j=3,jtot,2
    !   yrt(j-1)=(-4./(dy*dy))*(sin(float((j-1))*pi/(2.*jtot)))**2
    !   yrt(j)= yrt(j-1)
    ! end do

    ! call barrou()
    ! call ALL_ALL_j2(rhs,rhs2,0,ksen)
    ! call barrou()
    call transpose_z_to_y(pz,py)

    ! now in y-pencil, do FFT in y-direction
    !do k=1,ksen
    do k=1,ysize(3)!kmax2
      !do i=1,imax
      do i=1,ysize(1)!imax2
        !ipos = (k-1)*imax+i
        ipos = (k-1)*ysize(1) + i
        !do j=1,jtot
        do j=1,ysize(2)!jmax2=jtot
          !vfftj(ipos,j)=rhs2(i,j,k)
          FFTJ(ipos,j) = py(i,j,k)
        end do
      end do
    end do

    !call vrfftf(imax*ksen,jtot,vfftj,rhs2,imax*ksen,wj)
    call vrfftf(ysize(1)*ysize(3),ysize(2),FFTJ,py,ysize(1)*ysize(3),wj)

    !do k=1,ksen
    do k=1,ysize(3)!kmax2
      !do i=1,imax
      do i=1,ysize(1)!imax2
        !ipos=(k-1)*imax+i
        ipos = (k-1)*ysize(1) + i
        !do j=1,jtot
        do j=1,ysize(2)!jmax2=jtot
          !rhs2(i,j,k) = vfftj(ipos,j)
          py(i,j,k) = vfftj(ipos,j)
        end do
      end do
    end do

    ! call barrou()
    ! call ALL_ALL_j2(rhs,rhs2,1,ksen)
    ! call barrou()

    do j=1,jmax
      jv = j + myid*jmax
      !do i=1,imax
      do i=1,itot
        !bin(i) = b(i) + yrt(jv)
        bin(i) = bx(i) + yrt(jv)
      end do

    !do k=1,kmax
    do k=1,ktot
      do i=1,imax
         ipos=(k-1)*imax + i
         !y(i,k) = rhs(i,j,k)
         y(i,k) = pz(i,j,k)
      end do
    end do

    iperio = 1
    kperio = 1
    if (ibc1.eq.2) iperio = 0
    if (kbc1.eq.2) kperio = 0

    if (poisrcheck .eq. 0) then
      poisrcheck = 1
      !CALL BLKTRI(0,kperio,kmax,az,bz,cz,iperio,imax,a,bin,c,imax,y,ier,work) !0 for periodic BC
      CALL BLKTRI(0,kperio,kmax,a,b,c,iperio,imax,ax,bin,cx,imax,y,ier,work) !0 for periodic BC
      worksave = work
      write(6,*) 'First time step in POISR, poisrcheck=', poisrcheck
    end if

    !CALL BLKTRI(1,kperio,kmax,az,bz,cz,iperio,imax,a,bin,c,imax,y,ier,worksave)
    CALL BLKTRI(1,kperio,kmax,a,b,c,iperio,imax,ax,bin,cx,imax,y,ier,worksave)

      do k=1,kmax
      do i=1,imax
         ipos=(k-1)*imax+i
!PAR CHECK
!        vfftj(ipos,j) = y(i,k)
         rhs(i,j,k) = y(i,k)
      enddo
      enddo

      enddo         ! end loop over angles

1234  continue
!     call sumchk3(rhs,imax,jmax,kmax,9,1)
!
      call barrou()
      call ALL_ALL_j2(rhs,rhs2,0,ksen)
      call barrou()
!
      do k=1,ksen
      do i=1,imax
      ipos=(k-1)*imax+i
      do j=1,jtot
      vfftj(ipos,j)=rhs2(i,j,k)
      enddo
      enddo
      enddo


      call vrfftb(imax*ksen,jtot,vfftj,rhs2,imax*ksen,wj)
!!ATTTT      dok=1,kmax
      do k=1,ksen
      do i=1,imax
      ipos=(k-1)*imax+i
      do j=1,jtot
      rhs2(i,j,k)=vfftj(ipos,j)
!     if(abs(rhs(i,j,k)-help(i,j,k)).gt.1.e-13)then
!       write(6,*)'ERRR', i, j, k, rhs(i,j,k), help(i,j,k)
!     endif
      enddo
      enddo
      enddo
      call barrou()
      call ALL_ALL_j2(rhs,rhs2,1,ksen)
      call barrou()


!     call sumchk3(rhs,imax,jmax,kmax,2,1)


      deallocate(rhs2)
      deallocate(work)
      deallocate(vfftj)
      deallocate(y)
      return
  end subroutine poisr

  subroutine solmpj
  ! version: working version, barrou's removed,
  !          correct timing fft's
  !          AAPC with MPI-provided routines
  !          uses only 2 AAPC, using MPI-rovided routines,
  !          to pre-distribute the arrays s.t. complete
  !          2-D planes are present on each processor
  !          uses ALLTOALL instead of ALLTOALLV
  !          ONLY distribution in j-direction allowed

  ! NOTE: input array p1 is supposed to have the ip1ray distribution,
  !       i.e. the entire range of the first index must be present on
  !       each processor

  !******************************************************************
  !********************  FAST POISSON SOLVER ************************
  !*****                                                        *****
  !***               P_xx + P_yy + P_zz  =f(x,y,z)                ***
  !****                                                         *****
  !******************************************************************
  !   FOURIER TRANSFORMS IN X AND Y DIRECTION   GIVE:
  !   a^2 P + b^2 P + P_zz =  F(x,y,z) = FFT_i [ FTT_j (f(x,y,z))]

  !   where a and b are the KNOWN eigenvalues, and P_zz is

  !   P_zz =[ P_{i,j,k+1} - 2 P_{i,j,k} +P_{i,j,k-1} ] / (dz * dz)

  !   a^2 P + b^2 +P_zz =
  !   [P_{i,j,k+1}-(2+a^2+ b^2) P_{i,j,k}+P_{i,j,k-1}]/(dz*dz)=F( x,y,z)

  !   The equation above results in a tridiagonal system in k which
  !   can be solved with Gaussian elemination --> P
  !   The P we have found with the Gaussian elemination is still in
  !   the Fourier Space and 2 backward FFTS are necessary to compute
  !   the physical P
  !******************************************************************
  !******************************************************************
  !******************************************************************
  !****   Programmer: Bendiks Jan Boersma                      ******
  !****               email : b.j.boersma@wbmt.tudelft.nl      ******
  !****                                                        ******
  !****   USES      :  VFFTPACK   (netlib)                     ******
  !****             :  FFTPACK    (netlib)                     ******
  !****                (B.J. Boersma & L.J.P. Timmermans)      ******
  !****                                                        ******
  !******************************************************************
  !******************************************************************

  ! mpi-version, no master region for timing
  !              copy times all included
    use modmpi,    only : myid,comm3d,mpierr,nprocs, barrou
    use modglobal, only : imax,jmax,ktot,isen,itot,jtot,pi,dyi,dzf,dzh,dxfi, kb, ke, kh,kmax, ib, ie, jb, je, kb, ke, linoutflow, ierank, jerank, ibrank, jbrank
    use modfields, only : rhobf, rhobh

    implicit none

    real, allocatable, dimension(:,:,:) :: px, py, pz
    real, allocatable, dimension(:,:,:) :: d
    real, allocatable, dimension(:) :: FFTI, FFTJ, winew, wjnew
    real, dimension(1:ktot):: vout
    real    z,ak,bk,bbk,fac
    integer jv
    integer i, j, k

    call alloc_x(px, opt_xlevel=(/0,0,0/))
    call alloc_y(py, opt_ylevel=(/0,0,0/))
    call alloc_z(pz, opt_zlevel=(/0,0,0/))
    pz = p(ib:ie,jb:je,kb:ke)

    rhs = pz

    allocate(d(imax,jmax,kmax), FFTI(itot), FFTJ(jtot))

    if (.not. linoutflow) then
      allocate(winew(2*itot+15))
      call rffti(itot,winew)
      allocate(wjnew(2*jtot+15))
      call rffti(jtot,wjnew)
    else
      allocate(winew(3*itot+15))
      call costi(itot,winew)
      allocate(wjnew(3*jtot+15))
      call costi(jtot,wjnew)
    end if

    ! allocate(wjnew(2*jtot+15))
    ! call rffti(jtot,wjnew)

    ! ! Generate tridiagonal matrix
    !   do k=1,ktot
    !     ! SB fixed the coefficients
    !     a(k) = rhobh(k)  /(dzf(k)*dzh(k))
    !     c(k) = rhobh(k+1)/(dzf(k)*dzh(k+1))
    !     b(k) = -(a(k)+c(k))
    !   end do
    !   b(1) = b(1)+a(1)
    !   a(1) = 0.
    !   b(ktot) = b(ktot)+c(ktot)
    !   c(ktot) = 0.

    call transpose_z_to_y(pz, py)
    call transpose_y_to_x(py, px)

    !if (jbrank .and. ibrank) write(*,*) "nrank, p(:,1,1): ", nrank, px(:,1,1)

    ! Now in x-pencil, do FFT in x-direction
    if (.not. linoutflow) then
      fac = 1./sqrt(itot*1.)
      do k=1,xsize(3)
        do j=1,xsize(2)
          do i=1,xsize(1)!itot
            FFTI(i) = px(i,j,k)
          end do
          call rfftf(itot,FFTI,winew)
          do i=1,xsize(1)
            px(i,j,k) = FFTI(i)*fac
          end do
        end do
      end do

    else
      fac = 1./sqrt(2.*itot)
      !fac = 1./(4.*itot)
      !fac = 1.
      do k=1,xsize(3)
        do j=1,xsize(2)
          ! do i=1,xsize(1)!itot
          !   FFTI(i) = px(i,j,k)
          ! end do
          Sxr = px(:,j,k)
          !call cost(itot,FFTI,winew)
          call dfftw_execute(plan_r2fr_x)
          ! do i=1,xsize(1)
          !   px(i,j,k) = FFTI(i)*fac
          ! end do
          px(:,j,k) = Sxfr*fac
        end do
      end do

    end if

    call transpose_x_to_y(px, py)

    ! Now in y-pencil, do FFT in y-direction
    if (.not. linoutflow) then
      fac = 1./sqrt(jtot*1.)
      do i=1,ysize(1)
        do k=1,ysize(3)
          do j=1,ysize(2)!jtot
            FFTJ(j) = py(i,j,k)
          end do
          call rfftf(jtot,FFTJ,wjnew)
          do j=1,ysize(2)
            py(i,j,k)=FFTJ(j)*fac
          end do
        end do
      end do

    else
      fac = 1./sqrt(2.*jtot)
      !fac = 1./(4.*jtot)
      !fac = 1.
      do i=1,ysize(1)
        do k=1,ysize(3)
          ! do j=1,ysize(2)!jtot
          !   FFTJ(j) = py(i,j,k)
          ! end do
          Syr = py(i,:,k)
          call dfftw_execute(plan_r2fr_y)
          !call cost(jtot,FFTJ,wjnew)
          ! do j=1,ysize(2)
          !   py(i,j,k)=FFTJ(j)*fac
          ! end do
          py(i,:,k) = Syfr*fac
        end do
      end do

      !if (ibrank .and. jbrank) write(*,*) "nrank, Fy{p(1,:,1)}: ", nrank, py(1,:,1)

    end if

    call transpose_y_to_z(py,pz)

    ! Now in z-pencil, solve system using Gaussian elimination
    do j=1,zsize(2)
      !jv = j + myidy*jmax
      do i=1,zsize(1)
        z         = 1./(b(1)+xyzrt(i,j,1))
        d(i,j,1)  = c(1)*z
        pz(i,j,1) = pz(i,j,1)*z
      end do
    end do

    do k=2,zsize(3)-1
      do j=1,zsize(2)
        !jv = j + myid*jmax
        do i=1,zsize(1)
          bbk       = b(k)+xyzrt(i,j,k)
          z         = 1./(bbk-a(k)*d(i,j,k-1))
          d(i,j,k)  = c(k)*z
          pz(i,j,k) = (pz(i,j,k)-a(k)*pz(i,j,k-1))*z
        end do
      end do
    end do

    ak = a(ktot)
    bk = b(ktot)
    do j=1,zsize(2)
      !jv = j + myid*jmax
      do i=1,zsize(1)
        bbk = bk + xyzrt(i,j,ktot)
        z        = bbk-ak*d(i,j,ktot-1)
        if(z/=0.) then
          pz(i,j,ktot) = (pz(i,j,ktot)-ak*pz(i,j,ktot-1))/z
        else
          pz(i,j,ktot) =0.
        end if
      end do
    end do

    do k=zsize(3)-1,1,-1
      do j=1,zsize(2)
        do i=1,zsize(1)
          pz(i,j,k) = pz(i,j,k)-d(i,j,k)*pz(i,j,k+1)
        end do
      end do
    end do

    ! do j = 1,zsize(2)
    !   do i = 1,zsize(1)
    !     call tridiagonal(a, (b + xyzrt(i,j,:)), c, pz(i,j,:), vout)
    !     pz(i,j,:) = vout
    !   end do
    ! end do

    call transpose_z_to_y(pz,py)

    call transpose_y_to_x(py,px)
    !if (jbrank .and. ibrank) write(*,*) "p(:,1,1): ", px(:,1,1)
    call transpose_x_to_y(px,py)


    ! Now in y-pencil, do backward FFT in y direction
    if (.not. linoutflow) then
      fac = 1./sqrt(jtot*1.)
      do i=1,ysize(1)
        do k=1,ysize(3)
          do j=1,ysize(2)!jtot
            FFTJ(j) = py(i,j,k)
          end do
          call rfftb(jtot,FFTJ,wjnew)
          do j=1,ysize(2)
            py(i,j,k) = FFTJ(j)*fac
          end do
        end do
      end do

    else
      fac = 1./sqrt(2.*jtot)
      !fac = 1./(4.*jtot)
      !fac = 1./(2.*jtot)
      !fac = 1.
      do i=1,ysize(1)
        do k=1,ysize(3)
          ! do j=1,ysize(2)!jtot
          !   FFTJ(j) = py(i,j,k)
          ! end do
          Syfr = py(i,:,k)
          !call cost(jtot,FFTJ,wjnew)
          call dfftw_execute(plan_fr2r_y)
          ! do j=1,ysize(2)
          !   py(i,j,k) = FFTJ(j)*fac
          ! end do
          py(i,:,k) = Syr*fac
        end do
      end do
    end if

    call transpose_y_to_x(py,px)

    ! Now in x-pencil, do backward FFT in x-direction
    if (.not. linoutflow) then
      fac = 1./sqrt(itot*1.)
      do k=1,xsize(3)
        do j=1,xsize(2)
          do i=1,xsize(1)
            FFTI(i) = px(i,j,k)
          end do
          call rfftb(itot,FFTI,winew)
          do i=1,xsize(1)
            px(i,j,k) = FFTI(i)*fac
          end do
        end do
      end do

    else
      fac = 1./sqrt(2.*itot)
      !fac = 1./(4.*itot)
      !fac = 1./(2.*itot)
      !fac = 1.
      do k=1,xsize(3)
        do j=1,xsize(2)
          ! do i=1,xsize(1)
          !   FFTI(i) = px(i,j,k)
          ! end do
          Sxfr = px(:,j,k)
          !call cost(itot,FFTI,winew)
          call dfftw_execute(plan_fr2r_x)
          ! do i=1,xsize(1)
          !   px(i,j,k) = FFTI(i)*fac
          ! end do
          px(:,j,k) = Sxr*fac
        end do
      end do

    end if

    !if (jbrank .and. ibrank) write(*,*) "p(:,1,1): ", px(:,1,1)

    call transpose_x_to_y(px, py)
    call transpose_y_to_z(py, pz)

    p(ib:ie,jb:je,kb:ke) = pz

    deallocate(px,py,pz,d,FFTI,FFTJ,winew,wjnew)

    return
  end subroutine solmpj

  SUBROUTINE tridiagonal(a, b, c, r, u)
  real, dimension(:), intent(in)  :: a,b,c,r
  real, dimension(:), intent(out) :: u
  real, dimension(size(b)) :: gam
  integer :: n, j, k
  real    :: bet

  n = size(b)

  ! Forward substitution
  bet = b(1)
  u(1) = r(1) / bet
  do j = 2, n
    gam(j) = c(j-1) / bet
    bet = b(j) - a(j-1) * gam(j)
    u(j) = (r(j) - a(j-1) * u(j-1)) / bet
  end do

  ! Backward substitution
  do j = n-1, 1, -1
    u(j) = u(j) - gam(j+1) * u(j+1)
  end do
END SUBROUTINE tridiagonal


end module modpois
