!!> \file modsubgrid.f90
!!!  Calculates and applies the Sub Filter Scale diffusion
!
!>
!!  Calculates and applies the Sub Filter Scale diffusion
!>
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  \todo Documentation
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

module modsubgrid
implicit none
save
! private
public :: subgrid, initsubgrid,exitsubgrid
public :: ldelta, lmason,lsmagorinsky,cf, Rigc,prandtl, cm, cn, ch1, ch2, ce1, ce2, ekm,ekh, sbdiss,sbshr,sbbuo

  !cstep: set default values
  !cstep: user has the option to change ldelta, cf, cn, and Rigc in namoptions

  logical :: ldelta   = .false. !<  switch for subgrid length formulation (on/off)
  logical :: lmason   = .false. !<  switch for decreased length scale near the surface
  logical :: lsmagorinsky= .false. !<  switch for smagorinsky subgrid scheme
  logical :: ldynsub     = .false. !<  switch for dynamic subgrid scheme
  real :: cf      = 2.5  !< filter constant
  real :: Rigc    = 0.25 !< critical Richardson number
  real :: Prandtl = 3
  real :: cm      = 0.12
  real :: cn      = 0.76
  real :: ch1     = 1.
  real :: ch2     = 2.
  real :: ce1     = 0.19
  real :: ce2     = 0.51
  real :: cs      = -1
  real :: nmason  = 2.   !< exponent in Mason correction function
  real :: alpha_kolm   = 1.5     !< factor in Kolmogorov expression for spectral energy
  real :: beta_kolm    = 1.      !< factor in Kolmogorov relation for temperature spectrum

  real, allocatable :: ekm(:,:,:)  !<   k-coefficient for momentum
  real, allocatable :: ekh(:,:,:)  !<   k-coefficient for heat and q_tot
  real, allocatable :: sbdiss(:,:,:)!< dissiation
  real, allocatable :: sbshr(:,:,:) !< shear production
  real, allocatable :: sbbuo(:,:,:) !< buoyancy production / destruction
  real, allocatable :: zlt(:,:,:)  !<   filter width

  !CvH Allocate
  !CvH Dynamic subgrid model variables
  real, allocatable :: u_bar(:,:), v_bar(:,:), w_bar(:,:)
  real, allocatable :: u_hat(:,:), v_hat(:,:), w_hat(:,:)
  real, allocatable :: S11(:,:), S12(:,:), S13(:,:), S22(:,:), S23(:,:), S33(:,:)
  real, allocatable :: S11_bar(:,:), S12_bar(:,:), S13_bar(:,:), S22_bar(:,:), S23_bar(:,:), S33_bar(:,:)
  real, allocatable :: S11_hat(:,:), S12_hat(:,:), S13_hat(:,:), S22_hat(:,:), S23_hat(:,:), S33_hat(:,:)
  real, allocatable :: S_S11_bar(:,:), S_S12_bar(:,:), S_S13_bar(:,:), S_S22_bar(:,:), S_S23_bar(:,:), S_S33_bar(:,:)
  real, allocatable :: S_S11_hat(:,:), S_S12_hat(:,:), S_S13_hat(:,:), S_S22_hat(:,:), S_S23_hat(:,:), S_S33_hat(:,:)
  real, allocatable :: S(:,:), S_bar(:,:), S_hat(:,:)
  real, allocatable :: L11(:,:), L12(:,:), L13(:,:), L22(:,:), L23(:,:), L33(:,:)
  real, allocatable :: Q11(:,:), Q12(:,:), Q13(:,:), Q22(:,:), Q23(:,:), Q33(:,:)
  real, allocatable :: M11(:,:), M12(:,:), M13(:,:), M22(:,:), M23(:,:), M33(:,:)
  real, allocatable :: N11(:,:), N12(:,:), N13(:,:), N22(:,:), N23(:,:), N33(:,:)
  real, allocatable :: LM(:,:), MM(:,:), QN(:,:), NN(:,:)
  real, allocatable :: csz(:)

  real              :: const, beta, cs4_2, cs2_2
  real              :: LMav, LMavl, MMav, MMavl, QNav, QNavl, NNav, NNavl


contains
  subroutine initsubgrid
    use modglobal, only : ih,i1,jh,j1,k1,delta,zf,fkar, &
                          pi,ifnamopt,fname_options
    use modmpi,    only : myid, nprocs, comm3d, mpierr, my_real, mpi_logical, mpi_integer

    implicit none

    integer   :: k, ierr

    real :: ceps, ch
    real :: mlen

    namelist/NAMSUBGRID/ &
        ldelta,lmason, cf,cn,Rigc,Prandtl,lsmagorinsky,ldynsub,cs,nmason

    allocate(ekm(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(ekh(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(zlt(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(sbdiss(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(sbshr(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(sbbuo(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(csz(k1))

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMSUBGRID,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMSUBGRID'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMSUBGRID'
      endif
      write(6 ,NAMSUBGRID)
      close(ifnamopt)
    end if
    call MPI_BCAST(ldelta     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lmason     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(nmason     ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(lsmagorinsky,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(ldynsub     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(cs         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(cf         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(cn         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(Rigc       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(Prandtl    ,1,MY_REAL   ,0,comm3d,mpierr)

    cm = cf / (2. * pi) * (1.5*alpha_kolm)**(-1.5)

!     ch   = 2. * alpha_kolm / beta_kolm
    ch   = prandtl
    ch2  = ch-ch1

    ceps = 2. * pi / cf * (1.5*alpha_kolm)**(-1.5)
    ce1  = (cn**2)* (cm/Rigc - ch1*cm )
    ce2  = ceps - ce1

    if(cs == -1.) then
      csz(:)  = (cm**3/ceps)**0.25   !< Smagorinsky constant
    else
      csz(:)  = cs
    end if

    if(lmason .and. lsmagorinsky) then
      do k = 1,k1
        mlen   = (1. / (csz(k) * delta(k))**nmason + 1. / (fkar * zf(k))**nmason)**(-1./nmason)
        csz(k) = mlen / delta(k)
      end do
    end if

    if (myid==0) then
      write (6,*) 'cf    = ',cf
      write (6,*) 'cm    = ',cm
      write (6,*) 'ch    = ',ch
      write (6,*) 'ch1   = ',ch1
      write (6,*) 'ch2   = ',ch2
      write (6,*) 'ceps  = ',ceps
      write (6,*) 'ceps1 = ',ce1
      write (6,*) 'ceps2 = ',ce2
      write (6,*) 'cs    = ',cs
      write (6,*) 'Rigc  = ',Rigc
    endif

    !CvH init dynamic subgrid model
    if(ldynsub) then
      
      allocate(S(2-ih:i1+ih,2-jh:j1+jh)) 
      
      allocate(S11(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S12(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S13(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S22(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S23(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S33(2-ih:i1+ih,2-jh:j1+jh)) 

      allocate(u_bar(2-ih:i1+ih,2-jh:j1+jh))
      allocate(v_bar(2-ih:i1+ih,2-jh:j1+jh))
      allocate(w_bar(2-ih:i1+ih,2-jh:j1+jh))
      allocate(u_hat(2-ih:i1+ih,2-jh:j1+jh))
      allocate(v_hat(2-ih:i1+ih,2-jh:j1+jh))
      allocate(w_hat(2-ih:i1+ih,2-jh:j1+jh))
      
      allocate(S11_bar(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S12_bar(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S13_bar(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S22_bar(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S23_bar(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S33_bar(2-ih:i1+ih,2-jh:j1+jh)) 
      
      allocate(S11_hat(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S12_hat(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S13_hat(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S22_hat(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S23_hat(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S33_hat(2-ih:i1+ih,2-jh:j1+jh)) 
       
      allocate(S_S11_bar(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S_S12_bar(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S_S13_bar(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S_S22_bar(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S_S23_bar(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S_S33_bar(2-ih:i1+ih,2-jh:j1+jh)) 
      
      allocate(S_S11_hat(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S_S12_hat(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S_S13_hat(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S_S22_hat(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S_S23_hat(2-ih:i1+ih,2-jh:j1+jh))
      allocate(S_S33_hat(2-ih:i1+ih,2-jh:j1+jh)) 
  
      allocate(S_bar(2-ih:i1+ih,2-jh:j1+jh)) 
      allocate(S_hat(2-ih:i1+ih,2-jh:j1+jh)) 
      
      allocate(L11(2-ih:i1+ih,2-jh:j1+jh))
      allocate(L12(2-ih:i1+ih,2-jh:j1+jh))
      allocate(L13(2-ih:i1+ih,2-jh:j1+jh))
      allocate(L22(2-ih:i1+ih,2-jh:j1+jh))
      allocate(L23(2-ih:i1+ih,2-jh:j1+jh))
      allocate(L33(2-ih:i1+ih,2-jh:j1+jh))
  
      allocate(Q11(2-ih:i1+ih,2-jh:j1+jh))
      allocate(Q12(2-ih:i1+ih,2-jh:j1+jh))
      allocate(Q13(2-ih:i1+ih,2-jh:j1+jh))
      allocate(Q22(2-ih:i1+ih,2-jh:j1+jh))
      allocate(Q23(2-ih:i1+ih,2-jh:j1+jh))
      allocate(Q33(2-ih:i1+ih,2-jh:j1+jh)) 
  
      allocate(M11(2-ih:i1+ih,2-jh:j1+jh))
      allocate(M12(2-ih:i1+ih,2-jh:j1+jh))
      allocate(M13(2-ih:i1+ih,2-jh:j1+jh))
      allocate(M22(2-ih:i1+ih,2-jh:j1+jh))
      allocate(M23(2-ih:i1+ih,2-jh:j1+jh))
      allocate(M33(2-ih:i1+ih,2-jh:j1+jh)) 
  
      allocate(N11(2-ih:i1+ih,2-jh:j1+jh))
      allocate(N12(2-ih:i1+ih,2-jh:j1+jh))
      allocate(N13(2-ih:i1+ih,2-jh:j1+jh))
      allocate(N22(2-ih:i1+ih,2-jh:j1+jh))
      allocate(N23(2-ih:i1+ih,2-jh:j1+jh))
      allocate(N33(2-ih:i1+ih,2-jh:j1+jh)) 
  
      allocate(LM(2-ih:i1+ih,2-jh:j1+jh))
      allocate(MM(2-ih:i1+ih,2-jh:j1+jh))
      allocate(QN(2-ih:i1+ih,2-jh:j1+jh))
      allocate(NN(2-ih:i1+ih,2-jh:j1+jh))
    
    end if
   
!
  end subroutine initsubgrid
  subroutine subgrid

 ! Diffusion subroutines
! Thijs Heus, Chiel van Heerwaarden, 15 June 2007

    use modglobal, only : i1,ih,i2,j1,jh,j2,k1,nsv, lmoist
    use modfields, only : up,vp,wp,e12p,thl0,thlp,qt0,qtp,sv0,svp
    use modsurfdata,only : ustar,thlflux,qtflux,svflux
    implicit none
    integer n

    call closure
    call diffu(up)
    call diffv(vp)
    call diffw(wp)
    if (.not. lsmagorinsky) call diffe(e12p)
    call diffc(thl0,thlp,thlflux)
    if (lmoist) call diffc( qt0, qtp, qtflux)
    do n=1,nsv
      call diffc(sv0(:,:,:,n),svp(:,:,:,n),svflux(:,:,n))
    end do
    if ((.not. lsmagorinsky) .and. (.not. ldynsub)) call sources
  end subroutine

  subroutine exitsubgrid
    implicit none
    deallocate(ekm,ekh,zlt,sbdiss,sbbuo,sbshr,csz)
  end subroutine exitsubgrid

  subroutine filter(v2f, ndx)
    use modglobal, only : i1,ih,i2,j1,jh,j2,kmax,k1
    ! top hat filter
    ! for now, only filter in horizontal
    implicit none

    integer, intent(in)    :: ndx
    real,    intent(inout) :: v2f(2-ih:i1+ih,2-jh:j1+jh)
    real                   :: v2fin(2-ih:i1+ih,2-jh:j1+jh)
    integer                :: i,j

    v2fin(:,:) = v2f(:,:)
   if(ndx == 2) then
      do j = 2,j1
        do i = 2,i1
          !v2f(i,j) = 0.5 * (0.25 * v2fin(i-1,j) + 0.5 * v2fin(i,j) + 0.25 * v2fin(i+1,j)) &
          ! + 0.5 * (0.25 * v2fin(i,j-1) + 0.5 * v2fin(i,j) + 0.25 * v2fin(i,j+1))
          v2f(i,j) = 0.25 * ( &
              0.25 * v2fin(i-1,j-1) + 0.5 * v2fin(i,j-1) + 0.25 * v2fin(i+1,j-1) &
            + 0.5  * v2fin(i-1,j  ) +       v2fin(i,j  ) + 0.5  * v2fin(i+1,j  ) &
            + 0.25 * v2fin(i-1,j+1) + 0.5 * v2fin(i,j+1) + 0.25 * v2fin(i+1,j+1) )
        end do
      end do
    elseif(ndx == 4) then
      do j = 2,j1
        do i = 2,i1
          !v2f(i,j) = 0.5 * (0.125 * v2fin(i-2,j) + 0.25 * v2fin(i-1,j) + 0.25 * v2fin(i,j) + 0.25 * v2fin(i+1,j) + 0.125 * v2fin(i+2,j)) &
          ! + 0.5 * (0.125 * v2fin(i,j-2) + 0.25 * v2fin(i,j-1) + 0.25 * v2fin(i,j) + 0.25 * v2fin(i,j+1) + 0.125 * v2fin(i,j+2))
          v2f(i,j) = 0.0625 * ( &
              0.25 * v2fin(i-2,j-2) + 0.5 * v2fin(i-1,j-2) + 0.5 * v2fin(i  ,j-2) + 0.5 * v2fin(i+1,j-2) + 0.25 * v2fin(i+2,j-2) &
            + 0.5  * v2fin(i-2,j-1) + 1.0 * v2fin(i-1,j-1) + 1.0 * v2fin(i  ,j-1) + 1.0 * v2fin(i+1,j-1) + 0.5  * v2fin(i+2,j-1) &
            + 0.5  * v2fin(i-2,j  ) + 1.0 * v2fin(i-1,j  ) + 1.0 * v2fin(i  ,j  ) + 1.0 * v2fin(i+1,j  ) + 0.5  * v2fin(i+2,j  ) &
            + 0.5  * v2fin(i-2,j+1) + 1.0 * v2fin(i-1,j+1) + 1.0 * v2fin(i  ,j+1) + 1.0 * v2fin(i+1,j+1) + 0.5  * v2fin(i+2,j+1) &
            + 0.25 * v2fin(i-2,j+2) + 0.5 * v2fin(i-1,j+2) + 0.5 * v2fin(i  ,j+2) + 0.5 * v2fin(i+1,j+2) + 0.25 * v2fin(i+2,j+2) )
        end do
      end do
    end if

    return
  end subroutine filter

  subroutine closure

!-----------------------------------------------------------------|
!                                                                 |
!*** *closure*  calculates K-coefficients                         |
!                                                                 |
!      Hans Cuijpers   I.M.A.U.   06/01/1995                      |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!     All the K-closure factors are calculated.                   |
!                                                                 |
!     ekm(i,j,k) = k sub m : for velocity-closure                 |
!     ekh(i,j,k) = k sub h : for temperture-closure               |
!     ekh(i,j,k) = k sub h = k sub c : for concentration-closure  |
!                                                                 |
!     We will use the next model for these factors:               |
!                                                                 |
!     k sub m = 0.12 * l * sqrt(E)                                |
!                                                                 |
!     k sub h = k sub c = ( 1 + (2*l)/D ) * k sub m               |
!                                                                 |
!           where : l = mixing length  ( in model = z2 )          |
!                   E = subgrid energy                            |
!                   D = grid-size distance                        |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!             *closure* is called from *program*.                 |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal,   only : i1, j1,kmax,k1,ih,jh,i2,j2,delta,ekmin,grav, zf, fkar, &
                         dxi,dyi,dzf,dzh,rk3step,rslabs
  use modfields,   only : dthvdz,e120,u0,v0,w0
  use modsurfdata, only : dudz,dvdz,thvs,z0m

  use modmpi,    only : excjs, myid, nprocs, comm3d, mpierr, my_real, mpi_sum
  implicit none

  real    :: strain,mlen
  integer :: i,j,k,kp,km,jp,jm

!********************************************************************
!*********************************************************************
!  if (lsmagorinsky) then
!    do j=2,j1
!    do i=2,i1
!      k = 1
!      kp=k+1
!      km=k-1
!      jp=j+1
!      jm=j-1
!
!      strain =  ( &
!              ((u0(i+1,j,k)-u0(i,j,k))   *dxi        )**2    + &
!              ((v0(i,jp,k)-v0(i,j,k))    *dyi         )**2    + &
!              ((w0(i,j,kp)-w0(i,j,k))    /dzf(k)     )**2    )
!
!      strain = strain + 0.5 * ( &
!                ((w0(i,j,kp)-w0(i-1,j,kp))  *dxi     + &
!                (u0(i,j,kp)-u0(i,j,k))     / dzh(kp)  )**2    + &
!                ((w0(i,j,k)-w0(i-1,j,k))    *dxi     + &
!                dudz(i,j)   )**2    + &
!                ((w0(i+1,j,k)-w0(i,j,k))    *dxi     + &
!                dudz(i+1,j)   )**2    + &
!                ((w0(i+1,j,kp)-w0(i,j,kp))  *dxi     + &
!                (u0(i+1,j,kp)-u0(i+1,j,k)) / dzh(kp)  )**2    )
!
!      strain = strain + 0.5 * ( &
!                ((u0(i,jp,k)-u0(i,j,k))     *dyi     + &
!                (v0(i,jp,k)-v0(i-1,jp,k))  *dxi        )**2    + &
!                ((u0(i,j,k)-u0(i,jm,k))     *dyi     + &
!                (v0(i,j,k)-v0(i-1,j,k))    *dxi        )**2    + &
!                ((u0(i+1,j,k)-u0(i+1,jm,k)) *dyi     + &
!                (v0(i+1,j,k)-v0(i,j,k))    *dxi        )**2    + &
!                ((u0(i+1,jp,k)-u0(i+1,j,k)) *dyi     + &
!                (v0(i+1,jp,k)-v0(i,jp,k))  *dxi        )**2    )
!      strain = strain + 0.5 * ( &
!                ((v0(i,j,kp)-v0(i,j,k))     / dzh(kp) + &
!                (w0(i,j,kp)-w0(i,jm,kp))   *dyi        )**2    + &
!                (dvdz(i,j)+ &
!                (w0(i,j,k)-w0(i,jm,k))     *dyi        )**2    + &
!                (dvdz(i,j+1)+ &
!                (w0(i,jp,k)-w0(i,j,k))     *dyi        )**2    + &
!                ((v0(i,jp,kp)-v0(i,jp,k))   / dzh(kp) + &
!                (w0(i,jp,kp)-w0(i,j,kp))   *dyi        )**2    )
!      if(lmason) then
!        mlen        = (1. / (csz(k) * delta(k))**nmason + 1. / (fkar * (zf(k) + z0m(i,j)))**nmason)**(-1./nmason)
!      else
!        mlen        = csz(k) * delta(k)
!      end if
!      ekm(i,j,k)  = mlen**2.*sqrt(0.5*strain)
!      !if(i == 2 .and. j == 2) write(6,*)"CvH mlen", k, mlen, mlen / delta(k), delta(k)
!      ekh(i,j,k)  = ekm(i,j,k)/prandtl
!    end do
!    end do
!
!    do k=2,kmax
!    do j=2,j1
!    do i=2,i1
!      kp=k+1
!      km=k-1
!      jp=j+1
!      jm=j-1
!
!      strain =  ( &
!              ((u0(i+1,j,k)-u0(i,j,k))   *dxi        )**2    + &
!              ((v0(i,jp,k)-v0(i,j,k))    *dyi         )**2    + &
!              ((w0(i,j,kp)-w0(i,j,k))    /dzf(k)     )**2    )
!
!      strain = strain + 0.5 * ( &
!                ((w0(i,j,kp)-w0(i-1,j,kp))  *dxi     + &
!                (u0(i,j,kp)-u0(i,j,k))     / dzh(kp)  )**2    + &
!                ((w0(i,j,k)-w0(i-1,j,k))    *dxi     + &
!                (u0(i,j,k)-u0(i,j,km))     / dzh(k)   )**2    + &
!                ((w0(i+1,j,k)-w0(i,j,k))    *dxi     + &
!                (u0(i+1,j,k)-u0(i+1,j,km)) / dzh(k)   )**2    + &
!                ((w0(i+1,j,kp)-w0(i,j,kp))  *dxi     + &
!                (u0(i+1,j,kp)-u0(i+1,j,k)) / dzh(kp)  )**2    )
!
!      strain = strain + 0.5 * ( &
!                ((u0(i,jp,k)-u0(i,j,k))     *dyi     + &
!                (v0(i,jp,k)-v0(i-1,jp,k))  *dxi        )**2    + &
!                ((u0(i,j,k)-u0(i,jm,k))     *dyi     + &
!                (v0(i,j,k)-v0(i-1,j,k))    *dxi        )**2    + &
!                ((u0(i+1,j,k)-u0(i+1,jm,k)) *dyi     + &
!                (v0(i+1,j,k)-v0(i,j,k))    *dxi        )**2    + &
!                ((u0(i+1,jp,k)-u0(i+1,j,k)) *dyi     + &
!                (v0(i+1,jp,k)-v0(i,jp,k))  *dxi        )**2    )
!      strain = strain + 0.5 * ( &
!                ((v0(i,j,kp)-v0(i,j,k))     / dzh(kp) + &
!                (w0(i,j,kp)-w0(i,jm,kp))   *dyi        )**2    + &
!                ((v0(i,j,k)-v0(i,j,km))     / dzh(k)+ &
!                (w0(i,j,k)-w0(i,jm,k))     *dyi        )**2    + &
!                ((v0(i,jp,k)-v0(i,jp,km))   / dzh(k)+ &
!                (w0(i,jp,k)-w0(i,j,k))     *dyi        )**2    + &
!                ((v0(i,jp,kp)-v0(i,jp,k))   / dzh(kp) + &
!                (w0(i,jp,kp)-w0(i,j,kp))   *dyi        )**2    )
!      if(lmason) then
!        mlen        = (1. / (csz(k) * delta(k))**nmason + 1. / (fkar * (zf(k) + z0m(i,j)))**nmason)**(-1./nmason)
!      else
!        mlen        = csz(k) * delta(k)
!      end if
!      ekm(i,j,k)  = mlen**2.*sqrt(0.5*strain)
!      !if(i == 2 .and. j == 2) write(6,*)"CvH mlen", k, mlen, mlen / delta(k), delta(k)
!
!      ekh(i,j,k)  = ekm(i,j,k)/prandtl
!    end do
!    end do
!    end do

  if(lsmagorinsky) then
    if(ldynsub .and. rk3step == 1) then
      ! CvH dynamic subgrid model, compute cs only in first RK3 iteration step
      ! go through the model layer by layer
      do k = 1,kmax
        S(:,:)   = 0.
        
        S11(:,:) = 0.
        S12(:,:) = 0.
        S13(:,:) = 0.
        S22(:,:) = 0.
        S23(:,:) = 0.
        S33(:,:) = 0. 

        u_bar(:,:) = 0.
        v_bar(:,:) = 0.
        w_bar(:,:) = 0.
        u_hat(:,:) = 0.
        v_hat(:,:) = 0.
        w_hat(:,:) = 0.
     
        S11_bar(:,:) = 0.
        S12_bar(:,:) = 0.
        S13_bar(:,:) = 0.
        S22_bar(:,:) = 0.
        S23_bar(:,:) = 0.
        S33_bar(:,:) = 0. 
        
        S11_hat(:,:) = 0.
        S12_hat(:,:) = 0.
        S13_hat(:,:) = 0.
        S22_hat(:,:) = 0.
        S23_hat(:,:) = 0.
        S33_hat(:,:) = 0.
        
        S_S11_bar(:,:) = 0.
        S_S12_bar(:,:) = 0.
        S_S13_bar(:,:) = 0.
        S_S22_bar(:,:) = 0.
        S_S23_bar(:,:) = 0.
        S_S33_bar(:,:) = 0.
        
        S_S11_hat(:,:) = 0.
        S_S12_hat(:,:) = 0.
        S_S13_hat(:,:) = 0.
        S_S22_hat(:,:) = 0.
        S_S23_hat(:,:) = 0.
        S_S33_hat(:,:) = 0.
  
        S_bar(:,:) = 0.
        S_hat(:,:) = 0.
        
        L11(:,:) = 0.
        L12(:,:) = 0.
        L13(:,:) = 0.
        L22(:,:) = 0.
        L23(:,:) = 0.
        L33(:,:) = 0.
        
        Q11(:,:) = 0.
        Q12(:,:) = 0.
        Q13(:,:) = 0.
        Q22(:,:) = 0.
        Q23(:,:) = 0.
        Q33(:,:) = 0. 
        
        M11(:,:) = 0.
        M12(:,:) = 0.
        M13(:,:) = 0.
        M22(:,:) = 0.
        M23(:,:) = 0.
        M33(:,:) = 0. 
        
        N11(:,:) = 0.
        N12(:,:) = 0.
        N13(:,:) = 0.
        N22(:,:) = 0.
        N23(:,:) = 0.
        N33(:,:) = 0. 
  
        LM(:,:) = 0.
        MM(:,:) = 0.
        QN(:,:) = 0.
        NN(:,:) = 0.

        ! CHECK FOR NONEQUIDISTANT GRID!!!
        if(k == 1) then
          do j = 2 - jh + 1,j1 + jh - 1
            do i = 2 - ih + 1,i1 + ih - 1
              S11(i,j) =  0.5 * ( (u0(i+1,j,k) - u0(i,j,k)) * dxi &
                + (u0(i+1,j,k) - u0(i,j,k)) * dxi )          ! dudx + dudx
  
              S12(i,j) = 0.5 * ( 0.25*(u0(i,j+1,k)+u0(i+1,j+1,k) - (u0(i,j-1,k)+u0(i+1,j-1,k))) * dyi &
                + 0.25*(v0(i+1,j,k)+v0(i+1,j+1,k) - (v0(i-1,j,k)+v0(i-1,j+1,k))) * dxi )         ! dudy + dvdx
  
              S13(i,j) = 0.5 * ( dudz(i,j) &
                + 0.25*(w0(i+1,j,k)+w0(i+1,j,k+1) - (w0(i-1,j,k)+w0(i-1,j,k+1))) * dxi )         ! dudz + dwdx
  
              S22(i,j) = 0.5 * ( (v0(i,j+1,k) - v0(i,j,k)) * dyi &
                + (v0(i,j+1,k) - v0(i,j,k)) * dyi )         ! dvdy + dvdy
  
              S23(i,j) = 0.5 * ( dvdz(i,j) &
                + 0.25*(w0(i,j+1,k)+w0(i,j+1,k+1) - (w0(i,j-1,k)+w0(i,j-1,k+1))) * dyi )         ! dvdz + dwdy
  
              S33(i,j) = 0.5 * ( (w0(i,j,k+1) - w0(i,j,k)) / dzf(k) &
                +(w0(i,j,k+1) - w0(i,j,k)) / dzf(k) )       ! dwdz + dwdz
            end do
          end do

        else
          do j = 2 - jh + 1,j1 + jh - 1
            do i = 2 - ih + 1,i1 + ih - 1
              S11(i,j) =  0.5 * ( (u0(i+1,j,k) - u0(i,j,k)) * dxi &
                + (u0(i+1,j,k) - u0(i,j,k)) * dxi )          ! dudx + dudx
  
              S12(i,j) = 0.5 * ( 0.25*(u0(i,j+1,k)+u0(i+1,j+1,k) - (u0(i,j-1,k)+u0(i+1,j-1,k))) * dyi &
                + 0.25*(v0(i+1,j,k)+v0(i+1,j+1,k) - (v0(i-1,j,k)+v0(i-1,j+1,k))) * dxi )         ! dudy + dvdx
  
              S13(i,j) = 0.5 * ( 0.25*(u0(i,j,k+1)+u0(i+1,j,k+1) - (u0(i,j,k-1)+u0(i+1,j,k-1))) / dzf(k) &
                + 0.25*(w0(i+1,j,k)+w0(i+1,j,k+1) - (w0(i-1,j,k)+w0(i-1,j,k+1))) * dxi )         ! dudz + dwdx
  
              S22(i,j) = 0.5 * ( (v0(i,j+1,k) - v0(i,j,k)) * dyi &
                + (v0(i,j+1,k) - v0(i,j,k)) * dyi )         ! dvdy + dvdy
  
              S23(i,j) = 0.5 * ( 0.25*(v0(i,j,k+1)+v0(i,j+1,k+1) - (v0(i,j,k-1)+v0(i,j+1,k-1))) / dzf(k) &
                + 0.25*(w0(i,j+1,k)+w0(i,j+1,k+1) - (w0(i,j-1,k)+w0(i,j-1,k+1))) * dyi )         ! dvdz + dwdy
  
              S33(i,j) = 0.5 * ( (w0(i,j,k+1) - w0(i,j,k)) / dzf(k) &
                +(w0(i,j,k+1) - w0(i,j,k)) / dzf(k) )       ! dwdz + dwdz
            end do
          end do
        end if
 
        S(:,:) = sqrt(2.*(S11(:,:)**2. + S22(:,:)**2. + S33(:,:)**2. + &
          2. * (S12(:,:)**2. + S13(:,:)**2. + S23(:,:)**2.)))
    
        ! Unstagger the grid
        do j = 2 - jh,j1 + jh - 1
          do i = 2 - ih,i1 + ih - 1
            u_bar(i,j) = 0.5 * (u0(i,j,k) + u0(i+1,j,k))
            v_bar(i,j) = 0.5 * (v0(i,j,k) + v0(i,j+1,k))
            w_bar(i,j) = 0.5 * (w0(i,j,k) + w0(i,j,k+1))
          end do
        end do
     
        u_hat(:,:) = u_bar(:,:)
        v_hat(:,:) = v_bar(:,:)
        w_hat(:,:) = w_bar(:,:)
  
        L11(:,:) = u_bar(:,:) * u_bar(:,:)
        L12(:,:) = u_bar(:,:) * v_bar(:,:)
        L13(:,:) = u_bar(:,:) * w_bar(:,:)
        L22(:,:) = v_bar(:,:) * v_bar(:,:)
        L23(:,:) = v_bar(:,:) * w_bar(:,:)
        L33(:,:) = w_bar(:,:) * w_bar(:,:)
  
        Q11(:,:) = u_bar(:,:) * u_bar(:,:)
        Q12(:,:) = u_bar(:,:) * v_bar(:,:)
        Q13(:,:) = u_bar(:,:) * w_bar(:,:)
        Q22(:,:) = v_bar(:,:) * v_bar(:,:)
        Q23(:,:) = v_bar(:,:) * w_bar(:,:)
        Q33(:,:) = w_bar(:,:) * w_bar(:,:)
  
        ! filter the variable at twice the grid size
        call filter(u_bar,2)
        call filter(v_bar,2)
        call filter(w_bar,2)
  
        ! follow Bou-Zeid, 2005
        call filter(L11,2)
        L11(:,:) = L11(:,:) - u_bar(:,:)*u_bar(:,:)
        call filter(L12,2)
        L12(:,:) = L12(:,:) - u_bar(:,:)*v_bar(:,:)
        call filter(L13,2)
        L13(:,:) = L13(:,:) - u_bar(:,:)*w_bar(:,:)
        call filter(L22,2)
        L22(:,:) = L22(:,:) - v_bar(:,:)*v_bar(:,:)
        call filter(L23,2)
        L23(:,:) = L23(:,:) - v_bar(:,:)*w_bar(:,:)
        call filter(L33,2)
        L33(:,:) = L33(:,:) - w_bar(:,:)*w_bar(:,:)
  
        ! filter the variable at four times the grid size
        call filter(u_hat,4)
        call filter(v_hat,4)
        call filter(w_hat,4)

        ! follow Bou-Zeid, 2005
        call filter(Q11,4)
        Q11(:,:) = Q11(:,:) - u_hat(:,:)*u_hat(:,:)
        call filter(Q12,4)
        Q12(:,:) = Q12(:,:) - u_hat(:,:)*v_hat(:,:)
        call filter(Q13,4)
        Q13(:,:) = Q13(:,:) - u_hat(:,:)*w_hat(:,:)
        call filter(Q22,4)
        Q22(:,:) = Q22(:,:) - v_hat(:,:)*v_hat(:,:)
        call filter(Q23,4)
        Q23(:,:) = Q23(:,:) - v_hat(:,:)*w_hat(:,:)
        call filter(Q33,4)
        Q33(:,:) = Q33(:,:) - w_hat(:,:)*w_hat(:,:)
  
        S11_bar(:,:) = S11(:,:)
        S12_bar(:,:) = S12(:,:)
        S13_bar(:,:) = S13(:,:)
        S22_bar(:,:) = S22(:,:)
        S23_bar(:,:) = S23(:,:)
        S33_bar(:,:) = S33(:,:)
  
        S11_hat(:,:) = S11_bar(:,:)
        S12_hat(:,:) = S12_bar(:,:)
        S13_hat(:,:) = S13_bar(:,:)
        S22_hat(:,:) = S22_bar(:,:)
        S23_hat(:,:) = S23_bar(:,:)
        S33_hat(:,:) = S33_bar(:,:)
  
        call filter(S11_bar,2)
        call filter(S12_bar,2)
        call filter(S13_bar,2)
        call filter(S22_bar,2)
        call filter(S23_bar,2)
        call filter(S33_bar,2)
  
        call filter(S11_hat,4)
        call filter(S12_hat,4)
        call filter(S13_hat,4)
        call filter(S22_hat,4)
        call filter(S23_hat,4)
        call filter(S33_hat,4)
  
        S_bar(:,:) = sqrt(2.*(S11_bar(:,:)**2. + S22_bar(:,:)**2. + S33_bar(:,:)**2. + &
          2. *(S12_bar(:,:)**2. + S13_bar(:,:)**2. + S23_bar(:,:)**2.)))
  
        S_hat(:,:) = sqrt(2.*(S11_hat(:,:)**2. + S22_hat(:,:)**2. + S33_hat(:,:)**2. + &
          2. *(S12_hat(:,:)**2. + S13_hat(:,:)**2. + S23_hat(:,:)**2.)))
  
        S_S11_bar(:,:) = S(:,:) * S11(:,:)
        S_S12_bar(:,:) = S(:,:) * S12(:,:)
        S_S13_bar(:,:) = S(:,:) * S13(:,:)
        S_S22_bar(:,:) = S(:,:) * S22(:,:)
        S_S23_bar(:,:) = S(:,:) * S23(:,:)
        S_S33_bar(:,:) = S(:,:) * S33(:,:)
  
        S_S11_hat(:,:) = S_S11_bar(:,:)
        S_S12_hat(:,:) = S_S12_bar(:,:)
        S_S13_hat(:,:) = S_S13_bar(:,:)
        S_S22_hat(:,:) = S_S22_bar(:,:)
        S_S23_hat(:,:) = S_S23_bar(:,:)
        S_S33_hat(:,:) = S_S33_bar(:,:)
  
        call filter(S_S11_bar,2)
        call filter(S_S12_bar,2)
        call filter(S_S13_bar,2)
        call filter(S_S22_bar,2)
        call filter(S_S23_bar,2)
        call filter(S_S33_bar,2)
  
        call filter(S_S11_hat,4)
        call filter(S_S12_hat,4)
        call filter(S_S13_hat,4)
        call filter(S_S22_hat,4)
        call filter(S_S23_hat,4)
        call filter(S_S33_hat,4)
  
        !Compute Mij and Nij as in Bou-Zeid, 2005
        const = 2. * delta(k) ** 2.
  
        M11 = const*(S_S11_bar - 4.  * S_bar * S11_bar)
        M12 = const*(S_S12_bar - 4.  * S_bar * S12_bar)
        M13 = const*(S_S13_bar - 4.  * S_bar * S13_bar)
        M22 = const*(S_S22_bar - 4.  * S_bar * S22_bar)
        M23 = const*(S_S23_bar - 4.  * S_bar * S23_bar)
        M33 = const*(S_S33_bar - 4.  * S_bar * S33_bar)
  
        N11 = const*(S_S11_hat - 16. * S_hat * S11_hat)
        N12 = const*(S_S12_hat - 16. * S_hat * S12_hat)
        N13 = const*(S_S13_hat - 16. * S_hat * S13_hat)
        N22 = const*(S_S22_hat - 16. * S_hat * S22_hat)
        N23 = const*(S_S23_hat - 16. * S_hat * S23_hat)
        N33 = const*(S_S33_hat - 16. * S_hat * S33_hat)
  
        ! Contracting the tensors
        LM  = L11*M11 + L22*M22 + L33*M33 + 2. * (L12*M12 + L13*M13 + L23*M23)
        MM  = M11*M11 + M22*M22 + M33*M33 + 2. * (M12*M12 + M13*M13 + M23*M23)
        
        QN  = Q11*N11 + Q22*N22 + Q33*N33 + 2. * (Q12*N12 + Q13*N13 + Q23*N23)
        NN  = N11*N11 + N22*N22 + N33*N33 + 2. * (N12*N12 + N13*N13 + N23*N23)

        LMavl = 0.
        MMavl = 0.
        QNavl = 0.
        NNavl = 0.

        LMav = 0.
        MMav = 0.
        QNav = 0.
        NNav = 0.

        do j = 2,j1
          do i = 2,i1
            LMavl = LMavl + LM(i,j)
            MMavl = MMavl + MM(i,j)
            QNavl = QNavl + QN(i,j)
            NNavl = NNavl + NN(i,j)
          end do
        end do

        call MPI_ALLREDUCE(LMavl, LMav, 1, MY_REAL, MPI_SUM, comm3d,mpierr)
        call MPI_ALLREDUCE(MMavl, MMav, 1, MY_REAL, MPI_SUM, comm3d,mpierr)
        call MPI_ALLREDUCE(QNavl, QNav, 1, MY_REAL, MPI_SUM, comm3d,mpierr)
        call MPI_ALLREDUCE(NNavl, NNav, 1, MY_REAL, MPI_SUM, comm3d,mpierr)

        LMav = LMav / rslabs
        MMav = MMav / rslabs
        QNav = QNav / rslabs
        NNav = NNav / rslabs

        LMav  = max(1.e-24, LMav)
        MMav  = max(1.e-24, MMav)
        QNav  = max(1.e-24, QNav)
        NNav  = max(1.e-24, NNav)

        cs2_2 = LMav / MMav
        cs4_2 = QNav / NNav
        
        cs2_2 = max(1.e-24, cs2_2)
        cs4_2 = max(1.e-24, cs4_2)

        beta = cs4_2 / cs2_2
        if(beta < 0.125) then
          write(6,*) "WARNING! Beta clip at k = ", k
          beta = max(cs4_2 / cs2_2, 0.125)
        end if

        !csz(k) = sqrt( (LMav / MMav) / ( (QNav * MMav) / (NNav * LMav) ))
        csz(k) = sqrt(cs2_2 / beta)

      end do
    end if
    ! End dynamic part of smagorinksy computation

    ! For ekm and ekh computation, revert to highest possible accuracy in strain
    ! discretization
    do k = 1,kmax
      mlen        = csz(k) * delta(k)

      do i = 2,i1
        do j = 2,j1

          kp=k+1
          km=k-1
          jp=j+1
          jm=j-1

          if(k == 1) then
            strain =  ( &
              ((u0(i+1,j,k)-u0(i,j,k))   *dxi        )**2    + &
              ((v0(i,jp,k)-v0(i,j,k))    *dyi        )**2    + &
              ((w0(i,j,kp)-w0(i,j,k))    /dzf(k)     )**2    )

            strain = strain + 0.5 * ( &
              ( 0.25*(w0(i+1,j,kp)-w0(i-1,j,kp))*dxi + &
              dudz(i,j)   )**2 )

            strain = strain + 0.125 * ( &
              ((u0(i,jp,k)-u0(i,j,k))     *dyi     + &
              (v0(i,jp,k)-v0(i-1,jp,k))  *dxi        )**2    + &
              ((u0(i,j,k)-u0(i,jm,k))     *dyi     + &
              (v0(i,j,k)-v0(i-1,j,k))    *dxi        )**2    + &
              ((u0(i+1,j,k)-u0(i+1,jm,k)) *dyi     + &
              (v0(i+1,j,k)-v0(i,j,k))    *dxi        )**2    + &
              ((u0(i+1,jp,k)-u0(i+1,j,k)) *dyi     + &
              (v0(i+1,jp,k)-v0(i,jp,k))  *dxi        )**2    )

            strain = strain + 0.125 * ( &
              ( 0.25*(w0(i,jp,kp)-w0(i,jm,kp))*dyi + &
              dvdz(i,j)   )**2 )
      
          else

            strain =  ( &
              ((u0(i+1,j,k)-u0(i,j,k))   *dxi        )**2    + &
              ((v0(i,jp,k)-v0(i,j,k))    *dyi        )**2    + &
              ((w0(i,j,kp)-w0(i,j,k))    /dzf(k)     )**2    )

            strain = strain + 0.125 * ( &
              ((w0(i,j,kp)-w0(i-1,j,kp))  *dxi     + &
              (u0(i,j,kp)-u0(i,j,k))     / dzh(kp)  )**2    + &
              ((w0(i,j,k)-w0(i-1,j,k))    *dxi     + &
              (u0(i,j,k)-u0(i,j,km))     / dzh(k)   )**2    + &
              ((w0(i+1,j,k)-w0(i,j,k))    *dxi     + &
              (u0(i+1,j,k)-u0(i+1,j,km)) / dzh(k)   )**2    + &
              ((w0(i+1,j,kp)-w0(i,j,kp))  *dxi     + &
              (u0(i+1,j,kp)-u0(i+1,j,k)) / dzh(kp)  )**2    )

            strain = strain + 0.125 * ( &
              ((u0(i,jp,k)-u0(i,j,k))     *dyi     + &
              (v0(i,jp,k)-v0(i-1,jp,k))  *dxi        )**2    + &
              ((u0(i,j,k)-u0(i,jm,k))     *dyi     + &
              (v0(i,j,k)-v0(i-1,j,k))    *dxi        )**2    + &
              ((u0(i+1,j,k)-u0(i+1,jm,k)) *dyi     + &
              (v0(i+1,j,k)-v0(i,j,k))    *dxi        )**2    + &
              ((u0(i+1,jp,k)-u0(i+1,j,k)) *dyi     + &
              (v0(i+1,jp,k)-v0(i,jp,k))  *dxi        )**2    )

            strain = strain + 0.125 * ( &
              ((v0(i,j,kp)-v0(i,j,k))     / dzh(kp) + &
              (w0(i,j,kp)-w0(i,jm,kp))   *dyi        )**2    + &
              ((v0(i,j,k)-v0(i,j,km))     / dzh(k)+ &
              (w0(i,j,k)-w0(i,jm,k))     *dyi        )**2    + &
              ((v0(i,jp,k)-v0(i,jp,km))   / dzh(k)+ &
              (w0(i,jp,k)-w0(i,j,k))     *dyi        )**2    + &
              ((v0(i,jp,kp)-v0(i,jp,k))   / dzh(kp) + &
              (w0(i,jp,kp)-w0(i,j,kp))   *dyi        )**2    )
          end if

          !if(i == 10 .and. j == 10) write(6,*) "strain:", k, "old", sqrt(0.5*strain), "new", S(10,10)
          
          ekm(i,j,k)  = mlen ** 2. * sqrt(2. * strain)
          ekh(i,j,k)  = ekm(i,j,k) / prandtl
        end do
      end do
    end do

  else
    do k=1,kmax
      do j=2,j1
        do i=2,i1
          if (ldelta .or. (dthvdz(i,j,k)<=0)) then
            zlt(i,j,k) = delta(k)
            if (lmason) zlt(i,j,k) = (1. / zlt(i,j,k) ** nmason + 1. / ( fkar * (zf(k) + z0m(i,j)))**nmason) ** (-1./nmason)
            ekm(i,j,k) = cm * zlt(i,j,k) * e120(i,j,k)
            ekh(i,j,k) = (ch1 + ch2) * ekm(i,j,k)
          else
            zlt(i,j,k) = min(delta(k),cn*e120(i,j,k)/sqrt(grav/thvs*abs(dthvdz(i,j,k))))
            if (lmason) zlt(i,j,k) = (1. / zlt(i,j,k) ** nmason + 1. / ( fkar * (zf(k) + z0m(i,j)))**nmason) ** (-1./nmason)

            ekm(i,j,k) = cm * zlt(i,j,k) * e120(i,j,k)
            ekh(i,j,k) = (ch1 + ch2 * zlt(i,j,k)/delta(k)) * ekm(i,j,k)
          endif
        end do
      end do
    end do
  end if

  ekm(:,:,:) = max(ekm(:,:,:),ekmin)
  ekh(:,:,:) = max(ekh(:,:,:),ekmin)

!*************************************************************
!     Set cyclic boundary condition for K-closure factors.
!*************************************************************

  ekm(1, :,:) = ekm(i1,:,:)
  ekm(i2,:,:) = ekm(2, :,:)
  ekh(1, :,:) = ekh(i1,:,:)
  ekh(i2,:,:) = ekh(2, :,:)

  call excjs( ekm           , 2,i1,2,j1,1,k1,ih,jh)
  call excjs( ekh           , 2,i1,2,j1,1,k1,ih,jh)

  do j=1,j2
    do i=1,i2
      ekm(i,j,k1)  = ekm(i,j,kmax)
      ekh(i,j,k1)  = ekh(i,j,kmax)
    end do
  end do

  return
  end subroutine closure
  subroutine sources


!-----------------------------------------------------------------|
!                                                                 |
!*** *sources*                                                    |
!      calculates various terms from the subgrid TKE equation     |
!                                                                 |
!     Hans Cuijpers   I.M.A.U.     06/01/1995                     |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!      Subroutine sources calculates all other terms in the       |
!      subgrid energy equation, except for the diffusion terms.   |
!      These terms are calculated in subroutine diff.             |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!     *sources* is called from *program*.                         |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal,   only : i1,j1,kmax,delta,dx,dy,dxi,dyi,dzf,zf,dzh,grav
  use modfields,   only : u0,v0,w0,e120,e12p,dthvdz
  use modsurfdata,  only : dudz,dvdz,thvs
  implicit none

  real    tdef2
  integer i,j,k,jm,jp,km,kp

  do k=2,kmax
  do j=2,j1
  do i=2,i1
    kp=k+1
    km=k-1
    jp=j+1
    jm=j-1

    tdef2 = 2. * ( &
             ((u0(i+1,j,k)-u0(i,j,k))   /dx         )**2    + &
             ((v0(i,jp,k)-v0(i,j,k))    /dy         )**2    + &
             ((w0(i,j,kp)-w0(i,j,k))    /dzf(k)     )**2    )

    tdef2 = tdef2 + 0.25 * ( &
              ((w0(i,j,kp)-w0(i-1,j,kp))  / dx     + &
               (u0(i,j,kp)-u0(i,j,k))     / dzh(kp)  )**2    + &
              ((w0(i,j,k)-w0(i-1,j,k))    / dx     + &
               (u0(i,j,k)-u0(i,j,km))     / dzh(k)   )**2    + &
              ((w0(i+1,j,k)-w0(i,j,k))    / dx     + &
               (u0(i+1,j,k)-u0(i+1,j,km)) / dzh(k)   )**2    + &
              ((w0(i+1,j,kp)-w0(i,j,kp))  / dx     + &
               (u0(i+1,j,kp)-u0(i+1,j,k)) / dzh(kp)  )**2    )

    tdef2 = tdef2 + 0.25 * ( &
              ((u0(i,jp,k)-u0(i,j,k))     / dy     + &
               (v0(i,jp,k)-v0(i-1,jp,k))  / dx        )**2    + &
              ((u0(i,j,k)-u0(i,jm,k))     / dy     + &
               (v0(i,j,k)-v0(i-1,j,k))    / dx        )**2    + &
              ((u0(i+1,j,k)-u0(i+1,jm,k)) / dy     + &
               (v0(i+1,j,k)-v0(i,j,k))    / dx        )**2    + &
              ((u0(i+1,jp,k)-u0(i+1,j,k)) / dy     + &
               (v0(i+1,jp,k)-v0(i,jp,k))  / dx        )**2    )

    tdef2 = tdef2 + 0.25 * ( &
              ((v0(i,j,kp)-v0(i,j,k))     / dzh(kp) + &
               (w0(i,j,kp)-w0(i,jm,kp))   / dy        )**2    + &
              ((v0(i,j,k)-v0(i,j,km))     / dzh(k)+ &
               (w0(i,j,k)-w0(i,jm,k))     / dy        )**2    + &
              ((v0(i,jp,k)-v0(i,jp,km))   / dzh(k)+ &
               (w0(i,jp,k)-w0(i,j,k))     / dy        )**2    + &
              ((v0(i,jp,kp)-v0(i,jp,k))   / dzh(kp) + &
               (w0(i,jp,kp)-w0(i,j,kp))   / dy        )**2    )


    sbshr(i,j,k)  = ekm(i,j,k)*tdef2/ ( 2*e120(i,j,k))
    sbbuo(i,j,k)  = -ekh(i,j,k)*grav/thvs*dthvdz(i,j,k)/ ( 2*e120(i,j,k))
    sbdiss(i,j,k) = - (ce1 + ce2*zlt(i,j,k)/delta(k)) * e120(i,j,k)**2 /(2.*zlt(i,j,k))
  end do
  end do
  end do
!     ----------------------------------------------end i,j,k-loop

!     --------------------------------------------
!     special treatment for lowest full level: k=1
!     --------------------------------------------


  do j=2,j1
    jp=j+1
    jm=j-1
  do i=2,i1



! **  Calculate "shear" production term: tdef2  ****************

    tdef2 =  2. * ( &
            ((u0(i+1,j,1)-u0(i,j,1))*dxi)**2 &
          + ((v0(i,jp,1)-v0(i,j,1))*dyi)**2 &
          + ((w0(i,j,2)-w0(i,j,1))/dzf(1))**2   )

    tdef2 = tdef2 + ( 0.25*(w0(i+1,j,2)-w0(i-1,j,2))*dxi + &
                          dudz(i,j)   )**2

    tdef2 = tdef2 +   0.25 *( &
          ((u0(i,jp,1)-u0(i,j,1))*dyi+(v0(i,jp,1)-v0(i-1,jp,1))*dxi)**2 &
         +((u0(i,j,1)-u0(i,jm,1))*dyi+(v0(i,j,1)-v0(i-1,j,1))*dxi)**2 &
         +((u0(i+1,j,1)-u0(i+1,jm,1))*dyi+(v0(i+1,j,1)-v0(i,j,1))*dxi)**2 &
         +((u0(i+1,jp,1)-u0(i+1,j,1))*dyi+ &
                                 (v0(i+1,jp,1)-v0(i,jp,1))*dxi)**2   )

    tdef2 = tdef2 + ( 0.25*(w0(i,jp,2)-w0(i,jm,2))*dyi + &
                          dvdz(i,j)   )**2

! **  Include shear and buoyancy production terms and dissipation **

    sbshr(i,j,1)  = ekm(i,j,1)*tdef2/ ( 2*e120(i,j,1))
    sbbuo(i,j,1)  = -ekh(i,j,1)*grav/thvs*dthvdz(i,j,1)/ ( 2*e120(i,j,1))
    sbdiss(i,j,1) = - (ce1 + ce2*zlt(i,j,1)/delta(1)) * e120(i,j,1)**2 /(2.*zlt(i,j,1))
  end do
  end do

  e12p(2:i1,2:j1,1:kmax) = e12p(2:i1,2:j1,1:kmax)+ &
            sbshr(2:i1,2:j1,1:kmax)+sbbuo(2:i1,2:j1,1:kmax)+sbdiss(2:i1,2:j1,1:kmax)

  return
  end subroutine sources

  subroutine diffc (putin,putout,flux)

    use modglobal, only : i1,ih,i2,j1,jh,j2,k1,kmax,dx2i,dzf,dy2i,dzh
    implicit none

    real, intent(in)    :: putin(2-ih:i1+ih,2-jh:j1+jh,k1)
    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real, intent(in)    :: flux (i2,j2)

    integer i,j,k,jm,jp,km,kp

    do k=2,kmax
      kp=k+1
      km=k-1

      do j=2,j1
        jp=j+1
        jm=j-1

        do i=2,i1
          putout(i,j,k) = putout(i,j,k) &
                    +  0.5 *  ( &
                  ( (ekh(i+1,j,k)+ekh(i,j,k))*(putin(i+1,j,k)-putin(i,j,k)) &
                    -(ekh(i,j,k)+ekh(i-1,j,k))*(putin(i,j,k)-putin(i-1,j,k)))*dx2i &
                    + &
                  ( (ekh(i,jp,k)+ekh(i,j,k)) *(putin(i,jp,k)-putin(i,j,k)) &
                    -(ekh(i,j,k)+ekh(i,jm,k)) *(putin(i,j,k)-putin(i,jm,k)) )*dy2i &
                    + &
                  ( (dzf(kp)*ekh(i,j,k) + dzf(k)*ekh(i,j,kp)) &
                    *  (putin(i,j,kp)-putin(i,j,k)) / dzh(kp)**2 &
                    - &
                    (dzf(km)*ekh(i,j,k) + dzf(k)*ekh(i,j,km)) &
                    *  (putin(i,j,k)-putin(i,j,km)) / dzh(k)**2           )/dzf(k) &
                            )

        end do
      end do
    end do

    do j=2,j1
      do i=2,i1

        putout(i,j,1) = putout(i,j,1) &
                  + 0.5 * ( &
                ( (ekh(i+1,j,1)+ekh(i,j,1))*(putin(i+1,j,1)-putin(i,j,1)) &
                  -(ekh(i,j,1)+ekh(i-1,j,1))*(putin(i,j,1)-putin(i-1,j,1)) )*dx2i &
                  + &
                ( (ekh(i,j+1,1)+ekh(i,j,1))*(putin(i,j+1,1)-putin(i,j,1)) &
                  -(ekh(i,j,1)+ekh(i,j-1,1))*(putin(i,j,1)-putin(i,j-1,1)) )*dy2i &
                  + &
                ( (dzf(2)*ekh(i,j,1) + dzf(1)*ekh(i,j,2)) &
                  *  (putin(i,j,2)-putin(i,j,1)) / dzh(2)**2 &
                  + flux(i,j) *2.                        )/dzf(1) &
                          )

      end do
    end do

  end subroutine diffc



  subroutine diffe(putout)

    use modglobal, only : i1,ih,i2,j1,jh,j2,k1,kmax,dx2i,dzf,dy2i,dzh
    use modfields, only : e120
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    integer             :: i,j,k,jm,jp,km,kp

    do k=2,kmax
      kp=k+1
      km=k-1

      do j=2,j1
        jp=j+1
        jm=j-1

        do i=2,i1

          putout(i,j,k) = putout(i,j,k) &
                  +  1.0 *  ( &
              ((ekm(i+1,j,k)+ekm(i,j,k))*(e120(i+1,j,k)-e120(i,j,k)) &
              -(ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k)))*dx2i &
                  + &
              ((ekm(i,jp,k)+ekm(i,j,k)) *(e120(i,jp,k)-e120(i,j,k)) &
              -(ekm(i,j,k)+ekm(i,jm,k)) *(e120(i,j,k)-e120(i,jm,k)) )*dy2i &
                  + &
              ((dzf(kp)*ekm(i,j,k) + dzf(k)*ekm(i,j,kp)) &
              *(e120(i,j,kp)-e120(i,j,k)) / dzh(kp)**2 &
              -(dzf(km)*ekm(i,j,k) + dzf(k)*ekm(i,j,km)) &
              *(e120(i,j,k)-e120(i,j,km)) / dzh(k)**2        )/dzf(k) &
                            )

        end do
      end do
    end do

  !     --------------------------------------------
  !     special treatment for lowest full level: k=1
  !     --------------------------------------------

    do j=2,j1
      do i=2,i1

        putout(i,j,1) = putout(i,j,1) + &
            ( (ekm(i+1,j,1)+ekm(i,j,1))*(e120(i+1,j,1)-e120(i,j,1)) &
              -(ekm(i,j,1)+ekm(i-1,j,1))*(e120(i,j,1)-e120(i-1,j,1)) )*dx2i &
            + &
            ( (ekm(i,j+1,1)+ekm(i,j,1))*(e120(i,j+1,1)-e120(i,j,1)) &
              -(ekm(i,j,1)+ekm(i,j-1,1))*(e120(i,j,1)-e120(i,j-1,1)) )*dy2i &
            + &
              ( (dzf(2)*ekm(i,j,1) + dzf(1)*ekm(i,j,2)) &
              *  (e120(i,j,2)-e120(i,j,1)) / dzh(2)**2              )/dzf(1)

      end do
    end do

  end subroutine diffe


  subroutine diffu (putout)

    use modglobal, only : i1,ih,i2,j1,jh,j2,k1,kmax,dxi,dx2i,dzf,dy,dyi,dy2i,dzh, cu,cv
    use modfields, only : u0,v0,w0
    use modsurfdata,only : ustar
    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: emmo,emom,emop,empo
    real                :: fu
    real                :: ucu, upcu
    integer             :: i,j,k,jm,jp,km,kp

    do k=2,kmax
      kp=k+1
      km=k-1

      do j=2,j1
        jp=j+1
        jm=j-1

        do i=2,i1

          emom = ( dzf(km) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i-1,j,km) ) ) / &
                    ( 4.   * dzh(k) )

          emop = ( dzf(kp) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,kp) + ekm(i-1,j,kp) ) ) / &
                    ( 4.   * dzh(kp) )

          empo = 0.25 * ( &
                  ekm(i,j,k)+ekm(i,jp,k)+ekm(i-1,jp,k)+ekm(i-1,j,k)  )

          emmo = 0.25 * ( &
                  ekm(i,j,k)+ekm(i,jm,k)+ekm(i-1,jm,k)+ekm(i-1,j,k)  )


          putout(i,j,k) = putout(i,j,k) &
                  + &
                  ( ekm(i,j,k)  * (u0(i+1,j,k)-u0(i,j,k)) &
                    -ekm(i-1,j,k)* (u0(i,j,k)-u0(i-1,j,k)) ) * 2. * dx2i &
                  + &
                  ( empo * ( (u0(i,jp,k)-u0(i,j,k))   *dyi &
                            +(v0(i,jp,k)-v0(i-1,jp,k))*dxi) &
                    -emmo * ( (u0(i,j,k)-u0(i,jm,k))   *dyi &
                            +(v0(i,j,k)-v0(i-1,j,k))  *dxi)   ) / dy &
                  + &
                  ( emop * ( (u0(i,j,kp)-u0(i,j,k))   /dzh(kp) &
                            +(w0(i,j,kp)-w0(i-1,j,kp))*dxi) &
                    -emom * ( (u0(i,j,k)-u0(i,j,km))   /dzh(k) &
                            +(w0(i,j,k)-w0(i-1,j,k))  *dxi)   ) /dzf(k)

        end do
      end do
    end do

  !     --------------------------------------------
  !     special treatment for lowest full level: k=1
  !     --------------------------------------------

    do j=2,j1
      jp = j+1
      jm = j-1

      do i=2,i1

        empo = 0.25 * ( &
              ekm(i,j,1)+ekm(i,jp,1)+ekm(i-1,jp,1)+ekm(i-1,j,1)  )

        emmo = 0.25 * ( &
              ekm(i,j,1)+ekm(i,jm,1)+ekm(i-1,jm,1)+ekm(i-1,j,1)  )

        emop = ( dzf(2) * ( ekm(i,j,1) + ekm(i-1,j,1) )  + &
                    dzf(1) * ( ekm(i,j,2) + ekm(i-1,j,2) ) ) / &
                  ( 4.   * dzh(2) )


        ucu   = 0.5*(u0(i,j,1)+u0(i+1,j,1))+cu

        if(ucu >= 0.) then
          upcu  = max(ucu,1.e-10)
        else
          upcu  = min(ucu,-1.e-10)
        end if


        fu = ( 0.5*( ustar(i,j)+ustar(i-1,j) ) )**2  * &
                upcu/sqrt(upcu**2  + &
                ((v0(i,j,1)+v0(i-1,j,1)+v0(i,jp,1)+v0(i-1,jp,1))/4.+cv)**2)

        putout(i,j,1) = putout(i,j,1) &
                + &
              ( ekm(i,j,1)  * (u0(i+1,j,1)-u0(i,j,1)) &
              -ekm(i-1,j,1)* (u0(i,j,1)-u0(i-1,j,1)) ) * 2. * dx2i &
                + &
              ( empo * ( (u0(i,jp,1)-u0(i,j,1))   *dyi &
                        +(v0(i,jp,1)-v0(i-1,jp,1))*dxi) &
              -emmo * ( (u0(i,j,1)-u0(i,jm,1))   *dyi &
                        +(v0(i,j,1)-v0(i-1,j,1))  *dxi)   ) / dy &
                + &
              ( emop * ( (u0(i,j,2)-u0(i,j,1))    /dzh(2) &
                        +(w0(i,j,2)-w0(i-1,j,2))  *dxi) &
                -fu   ) / dzf(1)

      end do
    end do

  end subroutine diffu


  subroutine diffv (putout)

    use modglobal, only : i1,ih,i2,j1,jh,j2,k1,kmax,dx,dxi,dx2i,dzf,dyi,dy2i,dzh, cu,cv
    use modfields, only : u0,v0,w0
    use modsurfdata,only : ustar

    implicit none

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: emmo, eomm,eomp,epmo
    real                :: fv, vcv,vpcv
    integer             :: i,j,k,jm,jp,km,kp

    do k=2,kmax
      kp=k+1
      km=k-1

      do j=2,j1
        jp=j+1
        jm=j-1

        do i=2,i1

          eomm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i,jm,km) ) ) / &
                    ( 4.   * dzh(k) )

          eomp = ( dzf(kp) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                      dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) / &
                    ( 4.   * dzh(kp) )

          emmo = 0.25  * ( &
                ekm(i,j,k)+ekm(i,jm,k)+ekm(i-1,jm,k)+ekm(i-1,j,k)  )

          epmo = 0.25  * ( &
                ekm(i,j,k)+ekm(i,jm,k)+ekm(i+1,jm,k)+ekm(i+1,j,k)  )


        putout(i,j,k) = putout(i,j,k) &
                + &
              ( epmo * ( (v0(i+1,j,k)-v0(i,j,k))   *dxi &
                        +(u0(i+1,j,k)-u0(i+1,jm,k))*dyi) &
                -emmo * ( (v0(i,j,k)-v0(i-1,j,k))   *dxi &
                        +(u0(i,j,k)-u0(i,jm,k))    *dyi)   ) / dx &
                + &
              (ekm(i,j,k) * (v0(i,jp,k)-v0(i,j,k)) &
              -ekm(i,jm,k)* (v0(i,j,k)-v0(i,jm,k))  ) * 2. * dy2i &
                + &
              ( eomp * ( (v0(i,j,kp)-v0(i,j,k))    /dzh(kp) &
                        +(w0(i,j,kp)-w0(i,jm,kp))  *dyi) &
                -eomm * ( (v0(i,j,k)-v0(i,j,km))    /dzh(k) &
                        +(w0(i,j,k)-w0(i,jm,k))    *dyi)   ) / dzf(k)

        end do
      end do
    end do

  !     --------------------------------------------
  !     special treatment for lowest full level: k=1
  !     --------------------------------------------

    do j=2,j1
      jp = j+1
      jm = j-1
      do i=2,i1

        emmo = 0.25 * ( &
              ekm(i,j,1)+ekm(i,jm,1)+ekm(i-1,jm,1)+ekm(i-1,j,1)  )

        epmo = 0.25  * ( &
              ekm(i,j,1)+ekm(i,jm,1)+ekm(i+1,jm,1)+ekm(i+1,j,1)  )

        eomp = ( dzf(2) * ( ekm(i,j,1) + ekm(i,jm,1)  )  + &
                    dzf(1) * ( ekm(i,j,2) + ekm(i,jm,2) ) ) / &
                  ( 4.   * dzh(2) )

        vcv   = 0.5*(v0(i,j,1)+v0(i,j+1,1))+cv
        if(vcv >= 0.) then
          vpcv  = max(vcv,1.e-10)
        else
          vpcv  = min(vcv,-1.e-10)
        end if


        fv    = ( 0.5*( ustar(i,j)+ustar(i,j-1) ) )**2  * &
                    vpcv/sqrt(vpcv**2  + &
                ((u0(i,j,1)+u0(i+1,j,1)+u0(i,jm,1)+u0(i+1,jm,1))/4.+cu)**2)

        putout(i,j,1) = putout(i,j,1) &
                  + &
                  ( epmo * ( (v0(i+1,j,1)-v0(i,j,1))   *dxi &
                            +(u0(i+1,j,1)-u0(i+1,jm,1))*dyi) &
                    -emmo * ( (v0(i,j,1)-v0(i-1,j,1))   *dxi &
                            +(u0(i,j,1)-u0(i,jm,1))    *dyi)   ) / dx &
                  + &
                ( ekm(i,j,1) * (v0(i,jp,1)-v0(i,j,1)) &
                  -ekm(i,jm,1)* (v0(i,j,1)-v0(i,jm,1))  ) * 2. * dy2i &
                  + &
                ( eomp * ( (v0(i,j,2)-v0(i,j,1))     /dzh(2) &
                          +(w0(i,j,2)-w0(i,jm,2))    *dyi) &
                  -fv   ) / dzf(1)

      end do
    end do

  end subroutine diffv



  subroutine diffw(putout)

    use modglobal, only : i1,ih,i2,j1,jh,j2,k1,kmax,dx,dxi,dx2i,dy,dyi,dy2i,dzf,dzh
    use modfields, only : u0,v0,w0
    implicit none

  !*****************************************************************

    real, intent(inout) :: putout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real                :: emom, eomm, eopm, epom
    integer             :: i,j,k,jm,jp,km,kp

    do k=2,kmax
      kp=k+1
      km=k-1
      do j=2,j1
        jp=j+1
        jm=j-1
        do i=2,i1

          emom = ( dzf(km) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i-1,j,km) ) ) / &
                    ( 4.   * dzh(k) )

          eomm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i,jm,km) ) ) / &
                    ( 4.   * dzh(k) )

          eopm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jp,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i,jp,km) ) ) / &
                    ( 4.   * dzh(k) )

          epom = ( dzf(km) * ( ekm(i,j,k)  + ekm(i+1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i+1,j,km) ) ) / &
                    ( 4.   * dzh(k) )


          putout(i,j,k) = putout(i,j,k) &
                + &
                  ( epom * ( (w0(i+1,j,k)-w0(i,j,k))    *dxi &
                            +(u0(i+1,j,k)-u0(i+1,j,km)) /dzh(k) ) &
                    -emom * ( (w0(i,j,k)-w0(i-1,j,k))    *dxi &
                            +(u0(i,j,k)-u0(i,j,km))     /dzh(k) ))/dx &
                + &
                  ( eopm * ( (w0(i,jp,k)-w0(i,j,k))     *dyi &
                            +(v0(i,jp,k)-v0(i,jp,km))   /dzh(k) ) &
                    -eomm * ( (w0(i,j,k)-w0(i,jm,k))     *dyi &
                            +(v0(i,j,k)-v0(i,j,km))     /dzh(k) ))/dy &
                + &
                  ( ekm(i,j,k) * (w0(i,j,kp)-w0(i,j,k)) /dzf(k) &
                  -ekm(i,j,km)* (w0(i,j,k)-w0(i,j,km)) /dzf(km) ) * 2. &
                                                              / dzh(k)

        end do
      end do
    end do

  end subroutine diffw

end module
