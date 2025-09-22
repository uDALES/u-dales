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
use decomp_2d_constants

use modglobal, only : ipoiss, POISS_FFT2D, POISS_FFT2D_2DECOMP, POISS_FFT3D, &
                      BCxm, BCym, BCzp, BCxm_periodic, BCym_periodic
#if defined(_GPU)
use modcuda,   only : pup_d, pvp_d, pwp_d, p_d
#else
use modfields, only : pup, pvp, pwp
#endif
use modfields, only : p
implicit none

include "fftw3.f"

private
public :: initpois,poisson,exitpois,p,rhs,dpupdx,dpvpdy,dpwpdz,xyzrt,sp,Fxy,Fxyz,pij
save
  real, allocatable, target :: rhs(:,:,:)   ! rhs of pressure solver
  real, allocatable, target :: dpupdx(:,:,:)
  real, allocatable, target :: dpvpdy(:,:,:)
  real, allocatable, target :: dpwpdz(:,:,:)

  real,    allocatable :: px(:,:,:), py(:,:,:), pz(:,:,:)
  complex, allocatable :: Fx(:,:,:), Fz(:,:,:)
  real, allocatable, target :: Fxy(:,:,:), Fxyz(:,:,:)

  real, allocatable :: xrt(:), yrt(:), zrt(:)
  real, allocatable, target :: xyzrt(:,:,:), bxyzrt(:,:,:)

  !! coefficients for tridiagonal matrix
  real, allocatable :: b(:), c(:)
#if defined(_GPU)
  real, allocatable, pinned :: a(:)
  real, allocatable, pinned :: zz(:,:,:), dd(:,:,:)
  real, allocatable, device :: a_d(:)
  real, allocatable, device :: zz_d(:,:,:), dd_d(:,:,:)
#else
  real, allocatable :: a(:)
  real, allocatable :: zz(:,:,:), dd(:,:,:)
#endif

  integer(kind=selected_int_kind(18)) :: plan_r2fc_x, plan_r2fc_y, plan_fc2r_x, plan_fc2r_y
  integer(kind=selected_int_kind(18)) :: plan_r2fr_x, plan_r2fr_y, plan_fr2r_x, plan_fr2r_y
  integer(kind=selected_int_kind(18)) :: plan_r2fr_z, plan_fr2r_z
  real, allocatable :: Sxr(:), Sxfr(:), Syr(:), Syfr(:), Szr(:), Szfr(:)
  complex, allocatable :: Sxfc(:), Syfc(:)
  type(DECOMP_INFO) :: sp
  type(DECOMP_INFO) :: decomp_top
  real, allocatable :: pij(:)

  integer :: sp_zst3, sp_zst2, sp_zst1, sp_zen3, sp_zen2, sp_zen1

contains
  subroutine initpois
    use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh,imax,jmax,itot,jtot,ktot, &
                          dxi,dzh,dzf,dy,dyi,dxfi,dzfi,pi,linoutflow,&
                          BCtopm,BCtopm_pressure
    use modfields, only : rhobf, rhobh
    implicit none
    integer :: kbc1, kbc2
    integer :: i, j, k, iv, jv, kv
    real    :: dzi, fac, b_top_D, b_top_N
    logical, dimension(3) :: skip_c2c = [.false., .false., .true.]

#if defined(_GPU)
    allocate(pup_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(pvp_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(pwp_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(p_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
#else
    allocate(pup(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(pvp(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(pwp(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
#endif
    call alloc_z(p)
    ! call alloc_z(rhs, opt_levels=(/0,0,0/))
    ! call alloc_z(dpupdx)
    ! call alloc_z(dpvpdy)
    ! call alloc_z(dpwpdz)

    allocate(pij(kb:ke+kh))

    if (ipoiss == POISS_FFT2D) then

      call alloc_x(px, opt_levels=(/0,0,0/))
      call alloc_y(py, opt_levels=(/0,0,0/))
      call alloc_z(pz, opt_levels=(/0,0,0/))

      ! call alloc_z(Fxy, opt_levels=(/0,0,0/))
      ! call alloc_z(Fxyz, opt_levels=(/0,0,0/))

      allocate(xrt(itot))
      allocate(yrt(jtot))
      allocate(zrt(ktot))
      allocate(xyzrt(imax,jmax,ktot))
      allocate(a(ktot), b(ktot), c(ktot))
      allocate(bxyzrt(imax,jmax,ktot))

      ! generate Eigenvalues xrt and yrt
      if (BCxm == BCxm_periodic) then ! periodic
        fac = 1./(2.*itot)
        do i=3,itot,2
          xrt(i-1)=-4.*dxi*dxi*(sin(real((i-1))*pi*fac))**2
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
          xrt(i)=-4.*dxi*dxi*(sin(real((i-1))*pi*fac))**2
        end do

        allocate(Sxr(0:itot-1), Sxfr(0:itot-1))
        call dfftw_plan_r2r_1d(plan_r2fr_x,itot,Sxr,Sxfr,FFTW_REDFT10,FFTW_MEASURE)
        call dfftw_plan_r2r_1d(plan_fr2r_x,itot,Sxfr,Sxr,FFTW_REDFT01,FFTW_MEASURE)
      end if

      if (BCym == BCym_periodic) then ! periodic
        fac = 1./(2.*jtot)
        do j=3,jtot,2
          yrt(j-1)=-4.*dyi*dyi*(sin(real((j-1))*pi*fac))**2
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
          yrt(j)=-4.*dyi*dyi*(sin(real((j-1))*pi*fac))**2
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
          zrt(k)=-4.*dzi*dzi*(sin(real((k-1))*pi*fac))**2
        end do
        zrt(1) = 0.

        allocate(Szr(0:ktot-1), Szfr(0:ktot-1))
        call dfftw_plan_r2r_1d(plan_r2fr_z,ktot,Szr,Szfr,FFTW_REDFT10,FFTW_MEASURE)
        call dfftw_plan_r2r_1d(plan_fr2r_z,ktot,Szfr,Szr,FFTW_REDFT01,FFTW_MEASURE)
      end if

      ! In z-pencil in spectral space - same dims as physical space in this case
      do k=1,zsize(3)
        do j=1,zsize(2)
          jv = j-1 + zstart(2)
          do i=1,zsize(1)
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

    elseif (ipoiss == POISS_FFT2D_2DECOMP) then
      
      if ((BCxm .ne. BCxm_periodic) .or. (BCym .ne. BCym_periodic)) then
        write(*,*) "Error: Invalid choice for ipois in namoptions. Use ipoiss = 0"
        stop 1
      end if

      call decomp_2d_fft_init(PHYSICAL_IN_X, opt_skip_XYZ_c2c=skip_c2c)
      call decomp_info_init(itot/2+1, jtot, ktot, sp)

      call alloc_x(px, opt_levels=(/0,0,0/))
      call alloc_y(py, opt_levels=(/0,0,0/))
      call alloc_z(pz, opt_levels=(/0,0,0/))
      allocate(Fz(sp%zsz(1),sp%zsz(2),sp%zsz(3)))

      sp_zst1=sp%zst(1);sp_zst2=sp%zst(2);sp_zst3=sp%zst(3);
      sp_zen1=sp%zen(1);sp_zen2=sp%zen(2);sp_zen3=sp%zen(3);

      allocate(xrt(itot/2+1))
      allocate(yrt(jtot))
      allocate(zrt(ktot))
      allocate(xyzrt(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
      allocate(a(ktot), b(ktot), c(ktot))
      allocate(zz(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
      allocate(dd(sp%zsz(1),sp%zsz(2),sp%zsz(3)))

      fac = 1./itot
      do i=2,itot/2
        xrt(i)=-4.*dxi*dxi*(sin(real((i-1))*pi*fac))**2
      end do
      xrt(1) = 0.
      xrt(itot/2+1) = -4.*dxi*dxi

      fac = 1./jtot ! Periodic
      do j=1,jtot
        yrt(j)=-4.*dyi*dyi*(sin(real((j-1))*pi*fac))**2
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
      
      !! compute the coefficients necessary for the Gaussian elimination along z-direction
      do k=1,ktot
        do j=sp_zst2,sp_zen2
          do i=sp_zst1,sp_zen1
            if (k==1) then
              zz(i,j,k) = 1./( b(k) + xyzrt(i,j,k)                    )
              dd(i,j,k) = c(k)*zz(i,j,k)
            elseif (k==ktot) then
              zz(i,j,k) =      b(k) + xyzrt(i,j,k) - a(k)*dd(i,j,k-1)
            else
              zz(i,j,k) = 1./( b(k) + xyzrt(i,j,k) - a(k)*dd(i,j,k-1) )
              dd(i,j,k) = c(k)*zz(i,j,k)
            end if
          end do
        end do
      end do
      
#if defined(_GPU)
      allocate(a_d(ktot))
      allocate(zz_d(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
      allocate(dd_d(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
      a_d  = a
      zz_d = zz
      dd_d = dd
#endif

    elseif (ipoiss == POISS_FFT3D) then ! periodic in all 3 dimension (may not work correctly, needs proper inspection)
      
      write(*,*) "Error: Invalid choice for ipois in namoptions. Use eiter 0 or 3"
      stop 1

      call decomp_2d_fft_init(PHYSICAL_IN_Z) ! means z pencil
      call decomp_info_init(itot, jtot, ktot/2+1, sp) ! have to do this because sp is not public

      call alloc_z(pz, opt_levels=(/0,0,0/))
      allocate(Fx(sp%xsz(1),sp%xsz(2),sp%xsz(3)))

      !!! Generate wavenumbers assuming FFTW implementation

      dzi = dzfi(1) ! Assumes equidistant in z
      allocate(xrt(itot))
      allocate(yrt(jtot))
      allocate(zrt(ktot/2+1))

      fac = 1./itot
      do i=2,itot
        xrt(i)=-4.*dxi*dxi*(sin(real((i-1))*pi*fac))**2
      end do
      xrt(1    ) = 0.
      xrt(itot ) = -4.*dxi*dxi

      fac = 1./jtot
      do j=2,jtot
        yrt(j)=-4.*dyi*dyi*(sin(real((j-1))*pi*fac))**2
      end do
      yrt(1    ) = 0.
      yrt(jtot ) = -4.*dyi*dyi

      fac = 1./ktot
      do k=2,ktot/2
        zrt(k)=-4.*dzi*dzi*(sin(real((k-1))*pi*fac))**2
      end do
      zrt(1) = 0.
      zrt(ktot/2+1) = -4.*dzi*dzi

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

    end if ! ipoiss

  end subroutine initpois

  subroutine poisson
    use modglobal, only: ib,ie,jb,je,kb,ke,itot,jtot,ktot
    implicit none
    real    :: fac
    integer :: i, j, k

    call fillps
    ! rhs = p(ib:ie,jb:je,kb:ke)   ! needs to be uncommented for fielddump to write correct value of rhs

    select case (ipoiss)
    case (POISS_FFT2D)

      !$acc data create(px, py, pz)
#if defined(_GPU)
      !$acc kernels default(present)
      pz = p_d(ib:ie,jb:je,kb:ke)
      !$acc end kernels
#else
      pz = p(ib:ie,jb:je,kb:ke)
#endif
      call transpose_z_to_y(pz, py)
      call transpose_y_to_x(py, px)
      !$acc update host(px)
      !$acc end data

      ! In x-pencil, do FFT in x-direction
      if (BCxm == BCxm_periodic) then
        fac = 1./sqrt(itot*1.)
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

      !$acc data create(px, py)
      !$acc update device(px)
      call transpose_x_to_y(px, py)
      !$acc update host(py)
      !$acc end data

      ! In y-pencil, do FFT in y-direction
      if (BCym == BCym_periodic) then
        fac = 1./sqrt(jtot*1.)
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

      !$acc data create(py, pz)
      !$acc update device(py)
      call transpose_y_to_z(py,pz)
      !$acc update host(pz)
      !$acc end data

      ! Fxy = pz

      ! In z-pencil
      if (BCzp == 1) then
        call solmpj(pz)
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
        do k=1,zsize(3)
          do j=1,zsize(2)
           do i=1,zsize(1)
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
      ! Fxyz = pz

      !$acc data create(py, pz)
      !$acc update device(pz)
      call transpose_z_to_y(pz,py)
      !$acc update host(py)
      !$acc end data

      ! In y-pencil, do backward FFT in y direction
      if (BCym == BCym_periodic) then
        fac = 1./sqrt(jtot*1.)
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

      !$acc data create(px, py)
      !$acc update device(py)
      call transpose_y_to_x(py,px)
      !$acc update host(px)
      !$acc end data

      ! In x-pencil, do backward FFT in x-direction
      if (BCxm == BCxm_periodic) then
        fac = 1./sqrt(itot*1.)
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
      end if

      !$acc data create(px, py, pz)
      !$acc update device(px)
      call transpose_x_to_y(px, py)
      call transpose_y_to_z(py, pz)
#if defined (_GPU)
      !$acc kernels default(present)
      p_d(ib:ie,jb:je,kb:ke) = pz
      !$acc end kernels
#else
      p(ib:ie,jb:je,kb:ke) = pz
#endif
      !$acc end data

    case(POISS_FFT2D_2DECOMP)

#if defined(_GPU)
      !$acc data create(px, py, pz, Fz)

      !$acc kernels default(present)
      pz = p_d(ib:ie,jb:je,kb:ke)
      !$acc end kernels
#else
      pz = p(ib:ie,jb:je,kb:ke)
#endif

      !!!! Starting in z-pencil, transpose to x-pencil
      call transpose_z_to_y(pz, py)
      call transpose_y_to_x(py, px)

      !!!! Forward FFT
      call decomp_2d_fft_3d(px, Fz)

      !$acc kernels default(present)
      Fz = Fz/sqrt(1.*itot*jtot)
      !$acc end kernels

      !!!! Solve system using Gaussian elimination
      !$acc kernels default(present)
      do j=1,sp_zen2
        do i=1,sp_zen1
#if defined(_GPU)
          Fz(i,j,1) = Fz(i,j,1)*zz_d(i,j,1)
#else
          Fz(i,j,1) = Fz(i,j,1)*zz(i,j,1)
#endif
        end do
      end do
      !$acc end kernels

      !$acc kernels default(present)
      do k=2,sp_zen3-1
        do j=1,sp_zen2
          do i=1,sp_zen1
#if defined(_GPU)
            Fz(i,j,k) = (Fz(i,j,k)-a_d(k)*Fz(i,j,k-1))*zz_d(i,j,k)
#else
            Fz(i,j,k) = (Fz(i,j,k)-a(k)*Fz(i,j,k-1))*zz(i,j,k)
#endif
          end do
        end do
      end do
      !$acc end kernels

      !$acc kernels default(present)
      do j=1,sp_zen2
        do i=1,sp_zen1
#if defined(_GPU)
          if(zz_d(i,j,ktot)/=0.) then
            Fz(i,j,ktot) = (Fz(i,j,ktot)-a_d(ktot)*Fz(i,j,ktot-1))/zz_d(i,j,ktot)
          else
            Fz(i,j,ktot) =0.
          end if
#else
          if(zz(i,j,ktot)/=0.) then
            Fz(i,j,ktot) = (Fz(i,j,ktot)-a(ktot)*Fz(i,j,ktot-1))/zz(i,j,ktot)
          else
            Fz(i,j,ktot) =0.
          end if
#endif
        end do
      end do
      !$acc end kernels

      !$acc kernels default(present)
      do k=sp_zen3-1,1,-1
        do j=1,sp_zen2
          do i=1,sp_zen1
#if defined(_GPU)
            Fz(i,j,k) = Fz(i,j,k)-dd_d(i,j,k)*Fz(i,j,k+1)
#else
            Fz(i,j,k) = Fz(i,j,k)-dd(i,j,k)*Fz(i,j,k+1)
#endif
          end do
        end do
      end do
      !$acc end kernels

      !!!! Inverse FFT
      call decomp_2d_fft_3d(Fz, px)

      !$acc kernels default(present)
      px = px/sqrt(1.*itot*jtot)
      !$acc end kernels

      !!!! From x-pencil transpose back to z-pencil
      call transpose_x_to_y(px, py)
      call transpose_y_to_z(py, pz)

#if defined (_GPU)
      !$acc kernels default(present)
      p_d(ib:ie,jb:je,kb:ke) = pz
      !$acc end kernels

      !$acc end data
#else
      p(ib:ie,jb:je,kb:ke) = pz
#endif

    case (POISS_FFT3D)
      ! This block needs careful rectification for GPU run

      pz = p(ib:ie,jb:je,kb:ke)

      call decomp_2d_fft_3d(pz,Fx) ! start in z-pencil in physical space, end in x-pencil in Fourier space

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

    case default
      write(0, *) "Invalid choice for Poisson solver"
      stop 1
    end select

    call tderive

  end subroutine poisson

  subroutine exitpois
    implicit none

#if defined(_GPU)
    deallocate(p_d, pup_d, pvp_d, pwp_d)
#else
    deallocate(pup, pvp, pwp)
#endif

    deallocate(p, pij)
    ! deallocate(rhs, dpupdx, dpvpdy, dpwpdz)

    select case (ipoiss)
      case (POISS_FFT2D)
        deallocate(px, py, pz)
        ! deallocate(Fxy, Fxyz)
        deallocate(xrt, yrt, zrt, xyzrt, a, b, c, bxyzrt)
      case (POISS_FFT2D_2DECOMP)
        deallocate(px, py, pz, Fz)
        deallocate(xrt, yrt, zrt, xyzrt, a, b, c, zz, dd)
#if defined(_GPU)
        deallocate(a_d, zz_d, dd_d)
#endif
      case (POISS_FFT3D)
        deallocate(pz, Fx)
        deallocate(xrt, yrt, zrt, xyzrt)
      case default
        write(0, *) "Invalid choice for Poisson solver"
        stop 1
    end select
  end subroutine exitpois
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined(_GPU)
  attributes(global) subroutine compute_pup_cuda(rk3coefi)
    use modcuda, only    : ie_d, je_d, ke_d, &
                           up_d, vp_d, wp_d, um_d, vm_d, wm_d, &
                           tidandstride
    implicit none
    real   , value, intent(in) :: rk3coefi
    integer :: tidx, tidy, tidz, stridex, stridey, stridez
    integer :: i, j, k
    call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
    do k = tidz, ke_d, stridez
      do j = tidy, je_d, stridey
        do i = tidx, ie_d, stridex
          pup_d(i,j,k) = up_d(i,j,k) + um_d(i,j,k) * rk3coefi
          pvp_d(i,j,k) = vp_d(i,j,k) + vm_d(i,j,k) * rk3coefi
          pwp_d(i,j,k) = wp_d(i,j,k) + wm_d(i,j,k) * rk3coefi
        end do
      end do
    end do
  end subroutine compute_pup_cuda

  attributes(global) subroutine compute_p_cuda(dxi, dyi)
    use modcuda, only : ie_d, je_d, ke_d, dzfi_d, &
                        p_d, pup_d, pvp_d, pwp_d, tidandstride
    implicit none
    real   , value, intent(in) :: dxi, dyi
    integer :: tidx, tidy, tidz, stridex, stridey, stridez
    integer :: i, j, k
    call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
    do k = tidz, ke_d, stridez
      do j = tidy, je_d, stridey
        do i = tidx, ie_d, stridex
          p_d(i,j,k) = ( pup_d(i+1,j,k)-pup_d(i,j,k) ) * dxi &
                     + ( pvp_d(i,j+1,k)-pvp_d(i,j,k) ) * dyi &
                     + ( pwp_d(i,j,k+1)-pwp_d(i,j,k) ) * dzfi_d(k)
        end do
      end do
    end do
  end subroutine compute_p_cuda
#endif
  subroutine fillps
  ! Chiel van Heerwaarden,  19 June 2007
  ! Adapted fillps for RK3 time loop
    use modfields,   only : up, vp, wp, um, vm, wm
    use modglobal,   only : rk3step, ib, ie, jb, je, kb, ke, dxi, dyi, dzfi, dt, &
                            pi, dy, imax, jmax, ylen, xf, zf
    use modmpi,      only : myidx, myidy
    use modboundary, only : bcpup
#if defined(_GPU)
    use cudafor
    use modcuda,      only: griddim, blockdim, checkCUDA
#endif
    implicit none

    integer :: i, j, k
    real    :: rk3coef, rk3coefi

    if (rk3step == 0) then ! dt not defined yet
      rk3coef = 1.
    else
      rk3coef = dt / (4. - dble(rk3step))
    end if
    rk3coefi = 1. / rk3coef

#if defined(_GPU)
    call compute_pup_cuda<<<griddim,blockdim>>>(rk3coefi)
    call checkCUDA( cudaGetLastError(), 'compute_pup_cuda' )
#else
    do k=kb,ke
      do j=jb,je
        do i=ib,ie
          pup(i,j,k) = up(i,j,k) + um(i,j,k) * rk3coefi         ! see equation 5.81
          pvp(i,j,k) = vp(i,j,k) + vm(i,j,k) * rk3coefi
          pwp(i,j,k) = wp(i,j,k) + wm(i,j,k) * rk3coefi
        end do
      end do
    end do
#endif

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
    call bcpup(rk3coefi)   ! boundary conditions for pup,pvp,pwp

#if defined(_GPU)
    call compute_p_cuda<<<griddim,blockdim>>>(dxi, dyi)
    call checkCUDA( cudaGetLastError(), 'compute_p_cuda' )
#else
    do k=kb,ke
      do j=jb,je
        do i=ib,ie
          p(i,j,k) = ( pup(i+1,j,k)-pup(i,j,k) ) * dxi &         ! see equation 5.72
                   + ( pvp(i,j+1,k)-pvp(i,j,k) ) * dyi &
                   + ( pwp(i,j,k+1)-pwp(i,j,k) ) * dzfi(k)
        end do
      end do
    end do
#endif
    ! This block needs to be uncommented to write correct value of dpupdx, dpvpdy, dpwpdz through fielddump
    !do k=kb,ke
    !  do j=jb,je
    !    do i=ib,ie
    !      dpupdx(i,j,k) = (pup(i+1,j,k)-pup(i,j,k)) * dxi
    !      dpvpdy(i,j,k) = (pvp(i,j+1,k)-pvp(i,j,k)) * dyi
    !      dpwpdz(i,j,k) = (pwp(i,j,k+1)-pwp(i,j,k)) * dzfi(k)
    !    end do
    !  end do
    !end do

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
  end subroutine fillps

#if defined(_GPU)
  attributes(global) subroutine tderive_update_upvpwp_cuda(dxi, dyi)
    use modcuda, only: ie_d, je_d, kb_d, ke_d, up_d, vp_d, wp_d, p_d, dzhi_d, tidandstride
    implicit none
    real, value, intent(in) :: dxi, dyi
    integer :: tidx, tidy, tidz, stridex, stridey, stridez, i, j, k

    call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

    do i = tidx, ie_d, stridex
      do j = tidy, je_d, stridey
        do k = tidz, ke_d, stridez
          up_d(i,j,k) = up_d(i,j,k)-(p_d(i,j,k)-p_d(i-1,j,k))*dxi
          vp_d(i,j,k) = vp_d(i,j,k)-(p_d(i,j,k)-p_d(i,j-1,k))*dyi
          if (k .ne. kb_d) then
            wp_d(i,j,k) = wp_d(i,j,k)-(p_d(i,j,k)-p_d(i,j,k-1))*dzhi_d(k)
          end if
        end do
      end do
    end do
  end subroutine tderive_update_upvpwp_cuda
   
  attributes(global) subroutine tderive_update_wptop_cuda(pijke)
    use modcuda, only: ie_d, je_d, kb_d, ke_d, wp_d, dzhi_d, tidandstride
    implicit none
    real, value, intent(in) :: pijke
    integer :: tidx, tidy, tidz, stridex, stridey, stridez, i, j

    call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

    if (tidz == kb_d) then
      do i = tidx, ie_d, stridex
        do j = tidy, je_d, stridey
          wp_d(i,j,ke_d+1) = wp_d(i,j,ke_d+1) + 2*pijke*dzhi_d(ke_d+1)
        end do
      end do
    end if
  end subroutine tderive_update_wptop_cuda
#endif
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
    !-----------------------------------------------------------------|

    use modfields,   only: up, vp, wp, pres0, IIc, IIcs
    use modglobal,   only: ib,ie,ih,jb,je,jh,kb,ke,kh,dxi,dyi,dzhi,linoutflow,rslabs,BCtopm,BCtopm_pressure
    use modmpi,      only: slabsum, avexy_ibm
    use modboundary, only: bcp
#if defined(_GPU)
     use cudafor
     use modcuda,    only: griddim, blockdim, checkCUDA, pres0_d
#endif
    implicit none
    integer :: i, j, k

    ! **  Boundary conditions **************
    call bcp ! boundary conditions for p.

    !*****************************************************************
    ! **  Calculate time-derivative for the velocities with known ****
    ! **  pressure gradients.  ***************************************
    !*****************************************************************
#if defined(_GPU)
    call tderive_update_upvpwp_cuda<<<griddim,blockdim>>>(dxi, dyi)
    call checkCUDA( cudaGetLastError(), 'tderive_update_upvpwp_cuda' )
#else
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
#endif

    if (BCtopm .eq. BCtopm_pressure) then
      ! Get out the slab averaged dp/dz = <rhw>
#if defined(_GPU)
      p = p_d
#endif
      call avexy_ibm(pij(kb:ke+kh),p(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)

#if defined(_GPU)
      call tderive_update_wptop_cuda<<<griddim,blockdim>>>(pij(ke))
      call checkCUDA( cudaGetLastError(), 'tderive_update_wptop_cuda' )
#else
      do i=ib,ie
        do j=jb,je
          !wp(i,j,ke+1) = -(-pij(ke)-p(i,j,ke))*dzhi(ke+1) ! Doesn't work
          wp(i,j,ke+1) = wp(i,j,ke+1) + 2*pij(ke)*dzhi(ke+1)
        end do
      end do
#endif
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

    !$acc kernels default(present)
    do k=kb-1,ke+1
      do j=jb-1,je+1
        do i=ib-1,ie+1
#if defined(_GPU)
          pres0_d(i,j,k)=pres0_d(i,j,k)+p_d(i,j,k)
#else
          pres0(i,j,k)=pres0(i,j,k)+p(i,j,k)!-pijk ! update of the pressure: P_new = P_old + p
#endif
        end do
      end do
    end do
    !$acc end kernels

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
