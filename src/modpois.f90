!> \file modpois.f90
!!  Solves the Poisson equation for the pressure fluctuations
!>
!!  \author Jasper Tomas, TU Delft, 31 March 2014
!!  \author Harm Jonker, TU Delft
!!  \author Hans Cuijpers, IMAU
!!  \todo documentation
!!  \par Revision list
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
public :: initpois,poisson,exitpois,p,pup,pvp,pwp,rhs,dpupdx,dpvpdy,dpwpdz,xyzrt,sp,Fxy,Fxyz,dpdztop,pij
save

  real, allocatable, target :: p(:,:,:)   ! difference between p at previous and new step (p = P_new - P_old)
  real, allocatable, target :: pup(:,:,:), pvp(:,:,:), pwp(:,:,:)
  real, allocatable, target :: rhs(:,:,:)   ! rhs of pressure solver
  real, allocatable, target :: dpupdx(:,:,:)
  real, allocatable, target :: dpvpdy(:,:,:)
  real, allocatable, target :: dpwpdz(:,:,:)
  real, allocatable, target :: Fxy(:,:,:), Fxyz(:,:,:)
  real, allocatable :: xrt(:), yrt(:), zrt(:)
  real, allocatable, target :: xyzrt(:,:,:), bxyzrt(:,:,:)
  real, allocatable :: a(:), b(:), c(:), ax(:), bx(:), cx(:) ! coefficients for tridiagonal matrix
  real, allocatable :: pz_top(:,:,:), py_top(:,:,:), px_top(:,:,:)
  integer :: ibc1, ibc2, kbc1, kbc2

  integer :: plan_r2fc_x, plan_r2fc_y, plan_fc2r_x, plan_fc2r_y
  integer :: plan_r2fr_x, plan_r2fr_y, plan_fr2r_x, plan_fr2r_y
  integer :: plan_r2fr_z, plan_fr2r_z
  real, allocatable :: Sxr(:), Sxfr(:), Syr(:), Syfr(:), Szr(:), Szfr(:)
  complex, allocatable :: Sxfc(:), Syfc(:)
  type(DECOMP_INFO) :: sp
  type(DECOMP_INFO) :: decomp_top
  real, allocatable :: dpdztop(:,:)
  real, allocatable :: pij(:)

contains
  subroutine initpois
    use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh,imax,jmax,itot,jtot,ktot, &
                          dxi,dzh,dzf,dy,dyi,dxfi,dzfi,ipoiss,pi,linoutflow,&
                          POISS_FFT2D,POISS_FFT3D,POISS_CYC,POISS_FFT2D_2DECOMP,&
                          BCxm,BCym,BCzp,BCtopm,BCtopm_pressure
    use modmpi,    only : myidx, myidy, myid
    use modfields, only : rhobf, rhobh

    implicit none
    real dzi, fac, b_top_D, b_top_N
    integer i,j,k,iv,jv,kv

    allocate(dpdztop(ib:ie,jb:je), pij(kb:ke+kh))

    !allocate(p(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(pup(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(pvp(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(pwp(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    call alloc_z(p)
    call alloc_z(rhs, opt_zlevel=(/0,0,0/))
    call alloc_z(dpupdx)
    call alloc_z(dpvpdy)
    call alloc_z(dpwpdz)
    call alloc_z(Fxy, opt_zlevel=(/0,0,0/))
    call alloc_z(Fxyz, opt_zlevel=(/0,0,0/))

    if (ipoiss == POISS_FFT2D) then
      allocate(xrt(itot))
      allocate(yrt(jtot))
      allocate(zrt(ktot))
      allocate(xyzrt(imax,jmax,ktot))
      allocate(a(ktot), b(ktot), c(ktot))
      allocate(bxyzrt(imax,jmax,ktot))

      ! generate Eigenvalues xrt and yrt
      if (BCxm == 1) then ! periodic
        fac = 1./(2.*itot)
        do i=3,itot,2
          xrt(i-1)=-4.*dxi*dxi*(sin(float((i-1))*pi*fac))**2
          xrt(i)  = xrt(i-1)
        end do
        xrt(1    ) = 0.
        xrt(itot ) = -4.*dxi*dxi

        allocate(Sxr(0:itot-1), Sxfc(0:itot/2))
        call dfftw_plan_dft_r2c_1d(plan_r2fc_x,itot,Sxr,Sxfc,FFTW_MEASURE)
        call dfftw_plan_dft_c2r_1d(plan_fc2r_x,itot,Sxfc,Sxr,FFTW_MEASURE)

      else ! Neumann-Neumann
        fac = 1./(2.*itot)
        do i=1,itot
          xrt(i)=-4.*dxi*dxi*(sin(float((i-1))*pi*fac))**2
        end do

        allocate(Sxr(0:itot-1), Sxfr(0:itot-1))
        call dfftw_plan_r2r_1d(plan_r2fr_x,itot,Sxr,Sxfr,FFTW_REDFT10,FFTW_MEASURE)
        call dfftw_plan_r2r_1d(plan_fr2r_x,itot,Sxfr,Sxr,FFTW_REDFT01,FFTW_MEASURE)
      end if

      if (BCym == 1) then ! periodic
        fac = 1./(2.*jtot)
        do j=3,jtot,2
          yrt(j-1)=-4.*dyi*dyi*(sin(float((j-1))*pi*fac))**2
          yrt(j  )= yrt(j-1)
        end do
        yrt(1    ) = 0.
        yrt(jtot ) = -4.*dyi*dyi

        allocate(Syr(0:jtot-1), Syfc(0:jtot/2))
        call dfftw_plan_dft_r2c_1d(plan_r2fc_y,jtot,Syr,Syfc,FFTW_MEASURE)
        call dfftw_plan_dft_c2r_1d(plan_fc2r_y,jtot,Syfc,Syr,FFTW_MEASURE)

      else ! Neumann-Neumann
        fac = 1./(2.*jtot)
        do j=1,jtot
          yrt(j)=-4.*dyi*dyi*(sin(float((j-1))*pi*fac))**2
        end do

        allocate(Syr(0:jtot-1), Syfr(0:jtot-1))
        call dfftw_plan_r2r_1d(plan_r2fr_y,jtot,Syr,Syfr,FFTW_REDFT10,FFTW_MEASURE)
        call dfftw_plan_r2r_1d(plan_fr2r_y,jtot,Syfr,Syr,FFTW_REDFT01,FFTW_MEASURE)
      end if

      if (BCzp == 1) then ! solve using GE
        kbc1 = 1 ! Neumann
        ! top - Neumann
        kbc2 = 1

        ! generate matrix coefficients
        do k=1,ktot
          a(k) = rhobh(k)  /(dzf(k)*dzh(k))
          c(k) = rhobh(k+1)/(dzf(k)*dzh(k+1))
          b(k) = -(a(k)+c(k))
        end do

        ! bottom - Neumann
        b(1) = b(1) + a(1)
        !b(1) = -c(1)

        ! top
        b_top_N = b(ktot) + c(ktot)
        b_top_D = b(ktot) - c(ktot)
        if (kbc2 .eq. 1) then ! Neumann
          b(ktot) = b_top_N
          !b(ktot) = -a(ktot)
        elseif (kbc2 .eq. 3) then ! Dirichlet
          b(ktot) = b_top_D
          !b(ktot) = a(ktot)
        end if

        a(1) = 0.
        c(ktot) = 0.
        zrt = 0.

      else ! cosine transform
        dzi = dzfi(1) ! Assumes equidistant in z

        fac = 1./(2.*ktot)

        do k=2,ktot
          zrt(k)=-4.*dzi*dzi*(sin(float((k-1))*pi*fac))**2
        end do
        zrt(1) = 0.

        allocate(Szr(0:ktot-1), Szfr(0:ktot-1))
        call dfftw_plan_r2r_1d(plan_r2fr_z,ktot,Szr,Szfr,FFTW_REDFT10,FFTW_MEASURE)
        call dfftw_plan_r2r_1d(plan_fr2r_z,ktot,Szfr,Szr,FFTW_REDFT01,FFTW_MEASURE)

      end if

      ! In z-pencil in spectral space - same dims as physical space in this case
      do k=1,zsize(3)
        do j=1,zsize(2)
          !jv = j + myidy*zsize(2)!jmax
          jv = j-1 + zstart(2)
          do i=1,zsize(1)
            !iv = i + myidx*zsize(1)!imax ! is this correct in all cases? what if zsize is not constant on each rank?
            iv = i-1 + zstart(1)
            xyzrt(i,j,k) = rhobf(k)*(xrt(iv)+yrt(jv)+zrt(k))
          end do
        end do
      end do

      !if (BCtopm .eq. BCtop_pressure) then ! Set a Dirichlet condition for the zero mode
      do k=1,zsize(3)
        do j=1,zsize(2)
          do i=1,zsize(1)
            if ((xyzrt(i,j,k) .eq. 0.) .and. (k .eq. ktot)) then! .and. (BCtopm .eq. BCtopm_pressure)) then
              ! Set a Dirichlet BC ACROSS the top cell
              bxyzrt(i,j,k) = b_top_D !+ xyzrt(i,j,k) !(zero anyway)
            else
              bxyzrt(i,j,k) = b(k) + xyzrt(i,j,k)
            end if
          end do
        end do
      end do

      ! call decomp_info_init(itot, jtot, 1, decomp_top)
      ! allocate(px_top(1:decomp_top%xsz(1),1:decomp_top%xsz(2),1:decomp_top%xsz(3)))
      ! allocate(py_top(1:decomp_top%ysz(1),1:decomp_top%ysz(2),1:decomp_top%ysz(3)))
      ! allocate(pz_top(1:decomp_top%zsz(1),1:decomp_top%zsz(2),1:decomp_top%zsz(3)))


    elseif (ipoiss == POISS_FFT2D_2DECOMP) then
      call decomp_2d_fft_init(1) ! 1 means x pencil
      call decomp_info_init(itot/2+1, jtot, ktot, sp)
      allocate(xrt(itot/2+1))
      allocate(yrt(jtot))
      allocate(zrt(ktot))
      allocate(xyzrt(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
      allocate(a(ktot), b(ktot), c(ktot))

      fac = 1./itot
      do i=2,itot/2
        xrt(i)=-4.*dxi*dxi*(sin(float((i-1))*pi*fac))**2
      end do
      xrt(1) = 0.
      xrt(itot/2+1) = -4.*dxi*dxi

      ! This is not correct - needs to be complex
      ! fac = 1./(2.*jtot)
      ! do j=3,jtot,2
      !   yrt(j-1)=-4.*dyi*dyi*(sin(float((j-1))*pi*fac))**2
      !   yrt(j  )= yrt(j-1)
      ! end do
      ! yrt(1    ) = 0.
      ! yrt(jtot ) = -4.*dyi*dyi

      fac = 1./jtot ! Periodic
      do j=1,jtot
        yrt(j)=-4.*dyi*dyi*(sin(float((j-1))*pi*fac))**2
      end do

      kbc1 = 1 ! Neumann
      ! top - Neumann
      kbc2 = 1
      ! Could also set a Dirichlet BC on the inlet at the top, like in SPARKLE?

      ! generate matrix coefficients
      do k=1,ktot
        a(k) = rhobh(k)  /(dzf(k)*dzh(k))
        c(k) = rhobh(k+1)/(dzf(k)*dzh(k+1))
        b(k) = -(a(k)+c(k))
      end do

      ! bottom - Neumann
      b(1) = b(1) + a(1)
      !b(1) = -c(1)

      ! top
      if (kbc2 .eq. 1) then ! Neumann
        b(ktot) = b(ktot) + c(ktot)
        !b(ktot) = -a(ktot)
      elseif (kbc2 .eq. 3) then ! Dirichlet
        b(ktot) = b(ktot) - c(ktot)
        !b(ktot) = a(ktot)
      end if

      a(1) = 0.
      c(ktot) = 0.

      zrt = 0.

      ! In z-pencil in spectral space
      do k=1,sp%zsz(3)
        do j=1,sp%zsz(2)
          jv = j-1 + sp%zst(2)
          do i=1,sp%zsz(1)
            iv = i-1 + sp%zst(1)
            xyzrt(i,j,k) = rhobf(k)*(xrt(iv)+yrt(jv)+zrt(k))
          end do
        end do
      end do

    elseif (ipoiss == POISS_FFT3D) then ! periodic in all 3 dimension (probably not useful)
      call decomp_2d_fft_init(3) ! 3 means z pencil
      call decomp_info_init(itot, jtot, ktot/2+1, sp) ! have to do this because sp is not public
      ! Generate wavenumbers assuming FFTW implementation

      dzi = dzfi(1) ! Assumes equidistant in z
      allocate(xrt(itot))
      allocate(yrt(jtot))
      allocate(zrt(ktot/2+1))

      fac = 1./itot
      do i=2,itot
        xrt(i)=-4.*dxi*dxi*(sin(float((i-1))*pi*fac))**2
      end do
      xrt(1    ) = 0.
      xrt(itot ) = -4.*dxi*dxi

      fac = 1./jtot
      do j=2,jtot
        yrt(j)=-4.*dyi*dyi*(sin(float((j-1))*pi*fac))**2
      end do
      yrt(1    ) = 0.
      yrt(jtot ) = -4.*dyi*dyi

      fac = 1./ktot
      do k=2,ktot/2
        zrt(k)=-4.*dzi*dzi*(sin(float((k-1))*pi*fac))**2
      end do
      zrt(1) = 0.
      zrt(ktot/2+1) = -4.*dzi*dzi

      ! ! In z-pencil in spectral space
      ! allocate(xyzrt(imax,jmax,ktot)) ! should be ktot/2+1?
      ! do k=1,ktot/2+1
      !   do j=1,jmax
      !     jv = j + myidy*jmax
      !     do i=1,imax
      !       iv = i + myidx*imax
      !       xyzrt(i,j,k) = rhobf(k)*(xrt(iv)+yrt(jv)+zrt(k))
      !     end do
      !   end do
      ! end do

      ! In x-pencil in spectral space: sp%xsz(1)=itot
      allocate(xyzrt(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
      do k=1,sp%xsz(3)
        kv = k-1 + sp%xst(3)
        do j=1,sp%xsz(2)
          jv = j-1 + sp%xst(2)
          do i=1,sp%xsz(1)
            iv = i-1 + sp%xst(1)
            xyzrt(i,j,k) = rhobf(k)*(xrt(iv)+yrt(jv)+zrt(kv))
          end do
        end do
      end do


    elseif (ipoiss == POISS_CYC) then
      ! Not supported

      ! allocate(yrt(jtot))
      ! allocate(ax(itot), bx(itot), cx(itot))
      ! allocate(a(ktot), b(ktot), c(ktot))
      !
      ! ! spanwise - periodic
      ! yrt(1) = 0.
      ! yrt(jtot)=-4./(dy*dy)
      ! do j=3,jtot,2
      !   yrt(j-1)=(-4./(dy*dy))*(sin(float((j-1))*pi/(2.*jtot)))**2
      !   yrt(j)= yrt(j-1)
      ! end do
      !
      ! if (linoutflow ) then
      !   ibc1 = 1 ! Neumann
      !   ibc2 = 3 ! Dirichlet
      ! else
      !   ibc1 = 2
      !   ibc2 = 2
      ! endif
      ! kbc1 = 1
      ! kbc2 = 1
      !
      ! do i=1,itot
      !   ax(i) =  1./(dxf(i)*dxh(i))
      !   cx(i) =  1./(dxf(i)*dxh(i+1))
      !   bx(i) =  - (ax(i) + cx(i))
      ! end do
      !
      ! if (ibc1 .eq.1) then
      !   !bx(1) = bx(1) + ax(1)
      !   bx(1) = -cx(1)
      ! !elseif (ibc1 .eq. 2) then
      !   !bx(1) = bx(1)
      ! elseif (ibc1 .eq. 3) then
      !   bx(1) = bx(1) - ax(1) ! not convinced this is right
      !   !bx(1) = cx(1)
      ! end if
      !
      ! if (ibc2 .eq. 1) then
      !   !bx(itot) = bx(itot) + cx(itot)
      !   bx(itot) = -ax(1)
      ! !elseif (ibc2 .eq. 2) then
      !   !bx(itot) = bx(itot)
      ! elseif (ibc2 .eq. 3) then
      !   !bx(itot) = bx(itot) - cx(itot)
      !   bx(itot) = ax(itot)
      ! end if
      !
      ! if (ibc1 .ne. 2) then
      !   ax(1)    = 0.
      ! end if
      !
      ! if (ibc2 .ne. 2) then
      !   cx(itot) = 0.
      ! end if

    end if ! ipoiss

  end subroutine initpois

  subroutine poisson
    use modglobal, only : ib,ie,ih,jb,je,kb,ke,kh,itot,jtot,ktot,dy,dzf,dzh,linoutflow,iinletgen,ipoiss,POISS_FFT2D,POISS_FFT3D,POISS_CYC,POISS_FFT2D_2DECOMP,imax,jmax,eps1,BCxm,BCym,BCzp,ibrank,jbrank
    use modmpi, only : myid,nprocs,barrou
    implicit none
    integer ibc1,ibc2,kbc1,kbc2,ksen
    complex, allocatable, dimension(:,:,:) :: Fx, Fy, Fz
    real, allocatable, dimension(:,:,:) :: px, py, pz
    !real, allocatable, dimension(:) :: FFTI, FFTJ, winew, wjnew
    real, allocatable, dimension(:,:,:) :: Fzr, d
    real, dimension(1:ktot) :: vout
    real    z,ak,bk,bbk
    integer i, j, k
    real fac

    call fillps

    rhs = p(ib:ie,jb:je,kb:ke)

!  ibc?=1: neumann
!  ibc?=2: periodic
!  ibc?=3: dirichlet

    select case (ipoiss)
    case (POISS_FFT2D)
      call alloc_x(px, opt_xlevel=(/0,0,0/))
      call alloc_y(py, opt_ylevel=(/0,0,0/))
      call alloc_z(pz, opt_zlevel=(/0,0,0/))

      pz = p(ib:ie,jb:je,kb:ke)

      ! if (BCxm == 1) then
      !   allocate(FFTI(itot))
      !   allocate(winew(2*itot+15))
      !   call rffti(itot,winew)
      ! end if

      ! if (BCym == 1) then
      !   allocate(FFTJ(jtot))
      !   allocate(wjnew(2*jtot+15))
      !   call rffti(jtot,wjnew)
      ! end if

      call transpose_z_to_y(pz, py)
      call transpose_y_to_x(py, px)

      ! In x-pencil, do FFT in x-direction
      if (BCxm == 1) then
        fac = 1./sqrt(itot*1.)

        ! do k=1,xsize(3)
        !   do j=1,xsize(2)
        !     do i=1,xsize(1)!itot
        !       FFTI(i) = px(i,j,k)
        !     end do
        !     call rfftf(itot,FFTI,winew)
        !     do i=1,xsize(1)
        !       px(i,j,k) = FFTI(i)*fac
        !     end do
        !   end do
        ! end do

        do k=1,xsize(3)
          do j=1,xsize(2)
            Sxr = px(:,j,k)
            call dfftw_execute(plan_r2fc_x)
            px(1,j,k) = REAL(Sxfc(0))
            do i=1,itot/2-1
              px(2*i,j,k) = REAL(Sxfc(i))
              px(2*i+1,j,k) = AIMAG(Sxfc(i))
            end do
            px(itot,j,k) = REAL(Sxfc(itot/2))
          end do
        end do
        px = px*fac

      else
        fac = 1./sqrt(2.*itot)
        do k=1,xsize(3)
          do j=1,xsize(2)
            Sxr = px(:,j,k)
            call dfftw_execute(plan_r2fr_x)
            px(:,j,k) = Sxfr*fac
          end do
        end do

      end if

      call transpose_x_to_y(px, py)

      ! In y-pencil, do FFT in y-direction
      if (BCym == 1) then
        fac = 1./sqrt(jtot*1.)

        ! do i=1,ysize(1)
        !   do k=1,ysize(3)
        !     do j=1,ysize(2)!jtot
        !       FFTJ(j) = py(i,j,k)
        !     end do
        !     call rfftf(jtot,FFTJ,wjnew)
        !     do j=1,ysize(2)
        !       py(i,j,k)=FFTJ(j)*fac
        !     end do
        !   end do
        ! end do

        do i=1,ysize(1)
          do k=1,ysize(3)
            Syr = py(i,:,k)
            call dfftw_execute(plan_r2fc_y)
            py(i,1,k) = REAL(Syfc(0))
            do j=1,jtot/2-1
              py(i,2*j,k) = REAL(Syfc(j))
              py(i,2*j+1,k) = AIMAG(Syfc(j))
            end do
            py(i,jtot,k) = REAL(Syfc(jtot/2))
          end do
        end do
        py = py*fac

      else
        fac = 1./sqrt(2.*jtot)
        do i=1,ysize(1)
          do k=1,ysize(3)
            Syr = py(i,:,k)
            call dfftw_execute(plan_r2fr_y)
            py(i,:,k) = Syfr*fac
          end do
        end do

      end if

      call transpose_y_to_z(py,pz)
      Fxy = pz

      ! In z-pencil
      if (BCzp == 1) then
        call solmpj(pz)

      ! pz_top(:,:,1) = pz(:,:,ktot)
      ! if (ibrank .and. jbrank) pz_top(1,1,1) =-pz(1,1,ktot) ! Dirichlet for zero mode

      else
        ! Cosine transform in z
        fac = 1./sqrt(2.*ktot)
        do i=1,zsize(1)
          do j=1,zsize(2)
            Szr = pz(i,j,:)
            call dfftw_execute(plan_r2fr_z)
            pz(i,j,:) = Szfr*fac
          end do
        end do

        ! Divide by wavenumbers
        do i=1,imax
          do j=1,jmax
            do k=1,ktot
              if (xyzrt(i,j,k) .ne. 0.) then
                pz(i,j,k) = pz(i,j,k) / xyzrt(i,j,k)
              else
                pz(i,j,k) = 0.
              end if
            end do
          end do
        end do

        ! invert FFT
        fac = 1./sqrt(2.*ktot)
        do i=1,zsize(1)
          do j=1,zsize(2)
            Szfr = pz(i,j,:)
            call dfftw_execute(plan_fr2r_z)
            pz(i,j,:) = Szr*fac
          end do
        end do

      end if ! BCzp

      Fxyz = pz

      call transpose_z_to_y(pz,py)
      ! call transpose_z_to_y(pz_top, py_top, decomp_top)

      ! In y-pencil, do backward FFT in y direction
      if (BCym == 1) then
        fac = 1./sqrt(jtot*1.)

        ! do i=1,ysize(1)
        !   do k=1,ysize(3)
        !     do j=1,ysize(2)!jtot
        !       FFTJ(j) = py(i,j,k)
        !     end do
        !     call rfftb(jtot,FFTJ,wjnew)
        !     do j=1,ysize(2)
        !       py(i,j,k) = FFTJ(j)*fac
        !     end do
        !   end do
        ! end do

        do i=1,ysize(1)
          do k=1,ysize(3)
            Syfc(0) = CMPLX(py(i,1,k), 0.)
            do j=1,jtot/2-1
              Syfc(j) = CMPLX(py(i,2*j,k), py(i,2*j+1,k))
            end do
            Syfc(jtot/2) = CMPLX(py(i,jtot,k), 0.)
            call dfftw_execute(plan_fc2r_y)
            py(i,:,k) = Syr*fac
          end do
        end do

        ! if (size(py_top,3) > 0) then
        !   do i=1,decomp_top%ysz(1)
        !     Syfc(0) = CMPLX(py_top(i,1,1), 0.)
        !     do j=1,jtot/2-1
        !       Syfc(j) = CMPLX(py_top(i,2*j,1), py_top(i,2*j+1,1))
        !     end do
        !     Syfc(jtot/2) = CMPLX(py_top(i,jtot,1), 0.)
        !     call dfftw_execute(plan_fc2r_y)
        !     py_top(i,:,1) = Syr*fac
        !   end do
        ! end if

      else
        fac = 1./sqrt(2.*jtot)
        do i=1,ysize(1)
          do k=1,ysize(3)
            Syfr = py(i,:,k)
            call dfftw_execute(plan_fr2r_y)
            py(i,:,k) = Syr*fac
          end do
        end do
      end if

      call transpose_y_to_x(py,px)
      ! call transpose_y_to_x(py_top,px_top,decomp_top)

      ! In x-pencil, do backward FFT in x-direction
      if (BCxm == 1) then
        fac = 1./sqrt(itot*1.)

        ! do k=1,xsize(3)
        !   do j=1,xsize(2)
        !     do i=1,xsize(1)
        !       FFTI(i) = px(i,j,k)
        !     end do
        !     call rfftb(itot,FFTI,winew)
        !     do i=1,xsize(1)
        !       px(i,j,k) = FFTI(i)*fac
        !     end do
        !   end do
        ! end do

        do k=1,xsize(3)
          do j=1,xsize(2)
            Sxfc(0) = CMPLX(px(1,j,k), 0.)
            do i=1,itot/2-1
              Sxfc(i) = CMPLX(px(2*i,j,k), px(2*i+1,j,k))
            end do
            Sxfc(itot/2) = CMPLX(px(itot,j,k), 0.)
            call dfftw_execute(plan_fc2r_x)
            px(:,j,k) = Sxr*fac
          end do
        end do

      else
        fac = 1./sqrt(2.*itot)
        do k=1,xsize(3)
          do j=1,xsize(2)
            Sxfr = px(:,j,k)
            call dfftw_execute(plan_fr2r_x)
            px(:,j,k) = Sxr*fac
          end do
        end do

        ! if (size(px_top,3) > 0) then
        !   do j=1,decomp_top%xsz(2)
        !     Sxfr = px_top(:,j,1)
        !     call dfftw_execute(plan_fr2r_x)
        !     px_top(:,j,1) = Sxr*fac
        !   end do
        ! end if

      end if

      call transpose_x_to_y(px, py)
      call transpose_y_to_z(py, pz)

      ! call transpose_x_to_y(px_top,py_top,decomp_top)
      ! call transpose_y_to_z(py_top,pz_top,decomp_top)

      p(ib:ie,jb:je,kb:ke) = pz
      ! p(ib:ie,jb:je, ke+1) = pz_top(:,:,1)

      deallocate(px,py,pz)
      ! if (BCxm == 1) deallocate(FFTI, winew)
      ! if (BCym == 1) deallocate(FFTJ, wjnew)

    case(POISS_FFT2D_2DECOMP)
      write(0, *) 'ERROR: POISS_FFT2D_2DECOMP cannot be used.'
      stop 1
      ! DMajumdar: this POISS_FFT2D_2DECOMP has been commented as it's creating compilation error when
      ! using pbartholomew/2decomp-fft. If one wishes to use POISS_FFT2D_2DECOMP, needs to be rectified.
      
      ! call alloc_x(px, opt_xlevel=(/0,0,0/))
      ! call alloc_y(py, opt_ylevel=(/0,0,0/))
      ! call alloc_z(pz, opt_zlevel=(/0,0,0/))
      ! allocate(Fx(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
      ! allocate(Fy(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
      ! allocate(Fz(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
      ! allocate(d (sp%zsz(1),sp%zsz(2),sp%zsz(3)))

      ! pz = p(ib:ie,jb:je,kb:ke)

      ! ! Starting in z-pencil, transpose to x-pencil
      ! call transpose_z_to_y(pz, py)
      ! call transpose_y_to_x(py, px)

      ! ! Do forward FFT in x direction
      ! call r2c_1m_x(px, Fx)
      ! Fx = Fx/sqrt(1.*itot)

      ! ! Transpose to y-pencil
      ! call transpose_x_to_y(Fx, Fy, sp)

      ! ! Do forward FFT in y direction
      ! call c2c_1m_y(Fy, -1, plan(0,2))
      ! Fy = Fy/sqrt(1.*jtot)

      ! ! Transpose to z-pencil
      ! call transpose_y_to_z(Fy, Fz, sp)

      ! ! Solve system using Gaussian elimination
      ! do j=1,sp%zsz(2)
      !   do i=1,sp%zsz(1)
      !     z         = 1./(b(1)+xyzrt(i,j,1))
      !     d(i,j,1)  = c(1)*z
      !     Fz(i,j,1) = Fz(i,j,1)*z
      !   end do
      ! end do

      ! do k=2,sp%zsz(3)-1
      !   do j=1,sp%zsz(2)
      !     do i=1,sp%zsz(1)
      !       bbk       = b(k)+xyzrt(i,j,k)
      !       z         = 1./(bbk-a(k)*d(i,j,k-1))
      !       d(i,j,k)  = c(k)*z
      !       Fz(i,j,k) = (Fz(i,j,k)-a(k)*Fz(i,j,k-1))*z
      !     end do
      !   end do
      ! end do

      ! ak = a(ktot)
      ! bk = b(ktot)
      ! do j=1,sp%zsz(2)
      !   do i=1,sp%zsz(1)
      !     bbk = bk + xyzrt(i,j,ktot)
      !     z        = bbk-ak*d(i,j,ktot-1)
      !     if(z/=0.) then
      !       Fz(i,j,ktot) = (Fz(i,j,ktot)-ak*Fz(i,j,ktot-1))/z
      !     else
      !       Fz(i,j,ktot) =0.
      !     end if
      !   end do
      ! end do

      ! do k=sp%zsz(3)-1,1,-1
      !   do j=1,sp%zsz(2)
      !     do i=1,sp%zsz(1)
      !       Fz(i,j,k) = Fz(i,j,k)-d(i,j,k)*Fz(i,j,k+1)
      !     end do
      !   end do
      ! end do

      ! ! Tranpose to y-pencil
      ! call transpose_z_to_y(Fz, Fy, sp)

      ! ! Do backward FFT in y direction
      ! call c2c_1m_y(Fy, 1, plan(2,2))
      ! Fy = Fy/sqrt(1.*jtot)

      ! ! Transpose to x-pencil
      ! call transpose_y_to_x(Fy, Fx, sp)

      ! ! Do backward FFT in x direction
      ! call c2r_1m_x(Fx, px)
      ! px = px/sqrt(1.*itot)

      ! ! Tranpose to z-pencil
      ! call transpose_x_to_y(px, py)
      ! call transpose_y_to_z(py, pz)

      ! p(ib:ie,jb:je,kb:ke) = pz

      ! deallocate(px,py,pz,Fx,Fy,Fz,d)


    case (POISS_FFT3D)
      allocate(Fx(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
      allocate(Fy(sp%ysz(1),sp%ysz(2),sp%ysz(3)))
      allocate(Fz(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
      allocate(Fzr(imax,jmax,ktot))

      call alloc_z(pz, opt_zlevel=(/0,0,0/))
      pz = p(ib:ie,jb:je,kb:ke)

      call decomp_2d_fft_3d(pz,Fx) ! start in z-pencil in physical space, end in x-pencil in Fourier space

      ! call transpose_x_to_y(Fx,Fy,sp)
      ! call transpose_y_to_z(Fy,Fz,sp)
      !
      ! ! !Convert to real
      ! ! do i=ib,ie
      ! !   do j=jb,je
      ! !     Fzr(i,j,1) = REAL(Fz(i,j,1))
      ! !     do k=2,ktot/2
      ! !       Fzr(i,j,2*(k-1))   = REAL(Fz(i,j,k))
      ! !       Fzr(i,j,2*(k-1)+1) = AIMAG(Fz(i,j,k))
      ! !     end do
      ! !     Fzr(i,j,ktot) = REAL(Fz(i,j,ktot/2+1))
      ! !   end do
      ! ! end do
      ! ! Divide by wavenumbers in Fourier space
      ! do i=1,sp%zsz(1)
      !   do j=1,sp%zsz(2)
      !     do k=1,sp%zsz(3)
      !     !do k=1,ktot
      !       if (xyzrt(i,j,k) .ne. 0.) then
      !         Fz(i,j,k) = Fz(i,j,k) / CMPLX(xyzrt(i,j,k) * (itot*jtot*ktot))
      !         !Fzr(i,j,k) = Fzr(i,j,k) / (xyzrt(i,j,k) * (itot*jtot*ktot*2))
      !       else
      !         Fz(i,j,k) = CMPLX(0.)
      !       end if
      !
      !     end do
      !   end do
      ! end do
      !
      ! ! !convert back to complex
      ! ! do i=ib,ie
      ! !   do j=jb,je
      ! !     Fz(i,j,1) = CMPLX(Fzr(i,j,1), 0.)
      ! !     do k=2,ktot/2
      ! !       Fz(i,j,k) = CMPLX(Fzr(i,j,2*(k-1)), Fzr(i,j,2*(k-1)+1))
      ! !     end do
      ! !     Fz(i,j,ktot/2+1) = CMPLX(Fzr(i,j,ktot), 0.)
      ! !   end do
      ! ! end do
      !
      ! call transpose_z_to_y(Fz,Fy,sp)
      ! call transpose_y_to_x(Fy,Fx,sp)

      ! Divide by wavenumbers in Fourier space
      do i=1,sp%xsz(1)
        do j=1,sp%xsz(2)
          do k=1,sp%xsz(3)
            if (xyzrt(i,j,k) .ne. 0.) then
              Fx(i,j,k) = Fx(i,j,k) / CMPLX(xyzrt(i,j,k) * (itot*jtot*ktot))
            else
              Fx(i,j,k) = CMPLX(0.)
            end if
          end do
        end do
      end do

      call decomp_2d_fft_3d(Fx,pz)

      p(ib:ie,jb:je,kb:ke) = pz

      deallocate(Fx,Fy,Fz,pz)

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
      !call poisr
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


    use modfields, only : up, vp, wp, um, vm, wm,u0,v0,w0,IIw,IIws
    use modglobal, only : rk3step, ib,ie,jb,je,kb,ke,ih,jh,kh,dxi,dxfi,dyi,dzfi,dt,&
                          linoutflow,libm,dtmax,ierank,jerank,pi,dy,imax,jmax,ylen,xf,zf,jbrank
    use modmpi,    only : excjs, myidx, myidy, avexy_ibm, myid
    use modboundary, only: bcpup
!    use modibm,    only : ibmnorm

    implicit none
    !real,allocatable :: pup(:,:,:), pvp(:,:,:), pwp(:,:,:)
    integer i,j,k
    real rk3coef,rk3coefi
    real, dimension(kb:ke+kh) :: wpxy

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
          p(i,j,k)  =  ( pup(i+1,j,k)-pup(i,j,k) ) * dxi &         ! see equation 5.72
                          +( pvp(i,j+1,k)-pvp(i,j,k) ) * dyi &
                          +( pwp(i,j,k+1)-pwp(i,j,k) ) * dzfi(k)
        end do
      end do
    end do

    do k=kb,ke
      do j=jb,je
        do i=ib,ie
          dpupdx(i,j,k) = (pup(i+1,j,k)-pup(i,j,k)) * dxi
          dpvpdy(i,j,k) = (pvp(i,j+1,k)-pvp(i,j,k)) * dyi
          dpwpdz(i,j,k) = (pwp(i,j,k+1)-pwp(i,j,k)) * dzfi(k)
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
    use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh,dxi,dxhi,dyi,dzhi,linoutflow,rslabs,ibrank,ierank,jbrank,jerank,dxfi,BCtopm,BCtopm_pressure
    use modmpi,    only : myid,slabsum,avexy_ibm
    use modboundary,only : bcp
    implicit none
    integer i,j,k
    !real, dimension(kb:ke+kh) :: pij
    !real :: pijk

    ! **  Boundary conditions **************

    call bcp(p) ! boundary conditions for p.

    !*****************************************************************
    ! **  Calculate time-derivative for the velocities with known ****
    ! **  pressure gradients.  ***************************************
    !*****************************************************************

    do i=ib,ie
      do j=jb,je
        up(i,j,kb) = up(i,j,kb)-(p(i,j,kb)-p(i-1,j,kb))*dxi
        vp(i,j,kb) = vp(i,j,kb)-(p(i,j,kb)-p(i,j-1,kb))*dyi
        do k=kb+1,ke
          up(i,j,k) = up(i,j,k)-(p(i,j,k)-p(i-1,j,k))*dxi
          vp(i,j,k) = vp(i,j,k)-(p(i,j,k)-p(i,j-1,k))*dyi
          wp(i,j,k) = wp(i,j,k)-(p(i,j,k)-p(i,j,k-1))*dzhi(k)
        end do
      end do
    end do

    if (BCtopm .eq. BCtopm_pressure) then
      ! Get out the slab averaged dp/dz = <rhw>
      call avexy_ibm(pij(kb:ke+kh),p(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)

      do i=ib,ie
        do j=jb,je
          !wp(i,j,ke+1) = -(-pij(ke)-p(i,j,ke))*dzhi(ke+1) ! Doesn't work
          wp(i,j,ke+1) = wp(i,j,ke+1) + 2*pij(ke)*dzhi(ke+1)
        end do
      end do

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

    !pij =0.; pijk=0.;

    ! if (.not. linoutflow) then
    !   call slabsum(pij(kb:ke),kb,ke,p(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,ib,ie,jb,je,kb,ke)
    !   pij = pij/rslabs
    !   pijk = sum(pij(kb:ke))/(ke-kb)
    ! end if

    do k=kb-1,ke+1
      do j=jb-1,je+1
        do i=ib-1,ie+1
          pres0(i,j,k)=pres0(i,j,k)+p(i,j,k)!-pijk ! update of the pressure: P_new = P_old + p
        enddo
      enddo
    enddo

    return
  end subroutine tderive

  subroutine solmpj(x)
    use modmpi,    only : myid,comm3d,mpierr,nprocs, barrou
    use modglobal, only : imax,jmax,ktot,isen,itot,jtot,pi,dyi,dzf,dzh,dxfi, kb, ke, kh,kmax, ib, ie, jb, je, kb, ke, linoutflow, ierank, jerank, ibrank, jbrank, BCxm, BCym
    use modfields, only : rhobf, rhobh

    implicit none

    !real, allocatable, dimension(:,:,:) :: px, py, pz
    real, dimension(imax,jmax,ktot) :: d
    real, dimension(imax,jmax,ktot), intent(INOUT) :: x
    !real, allocatable, dimension(:) :: FFTI, FFTJ, winew, wjnew
    real    z,ak,bk,bbk,fac
    integer i, j, k

    do j=1,zsize(2)
      do i=1,zsize(1)
        !z         = 1./(b(1)+xyzrt(i,j,1))
        z         = 1./(bxyzrt(i,j,1))
        d(i,j,1)  = c(1)*z
        x(i,j,1) = x(i,j,1)*z
      end do
    end do

    do k=2,zsize(3)-1
      do j=1,zsize(2)
        do i=1,zsize(1)
          !bbk       = b(k)+xyzrt(i,j,k)
          bbk       = bxyzrt(i,j,k)
          z         = 1./(bbk-a(k)*d(i,j,k-1))
          d(i,j,k)  = c(k)*z
          x(i,j,k) = (x(i,j,k)-a(k)*x(i,j,k-1))*z
        end do
      end do
    end do

    ak = a(ktot)
    bk = b(ktot)
    do j=1,zsize(2)
      do i=1,zsize(1)
        !bbk = bk + xyzrt(i,j,ktot)
        bbk = bxyzrt(i,j,ktot)
        z        = bbk-ak*d(i,j,ktot-1)
        !if(z/=0.) then
          x(i,j,ktot) = (x(i,j,ktot)-ak*x(i,j,ktot-1))/z
        !else
        !  write(*,*) myid, i, j
        !  x(i,j,ktot) =0.
        !end if
      end do
    end do

    do k=zsize(3)-1,1,-1
      do j=1,zsize(2)
        do i=1,zsize(1)
          x(i,j,k) = x(i,j,k)-d(i,j,k)*x(i,j,k+1)
        end do
      end do
    end do

    return
  end subroutine solmpj

end module modpois
