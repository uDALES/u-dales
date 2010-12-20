!> \file modpois.f90
!!  Solves the Poisson equation for the pressure fluctuations

!>
!!  Solves the Poisson equation for the pressure fluctuations
!>
!!  \author Harm Jonker, TU Delft
!!  \author Hans Cuijpers, IMAU
!!  \todo documentation
!!  \par Revision list
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

implicit none
private
public :: initpois,poisson,exitpois,p
save

  real,allocatable :: p(:,:,:)                            ! pressure fluctuations

contains
  subroutine initpois
    use modglobal, only : i2,j2,k1
    implicit none

    allocate(p(i2,j2,0:k1))

  end subroutine initpois

  subroutine poisson
    implicit none
    call fillps
    call solmpj(p)
    call tderive

  end subroutine poisson

  subroutine exitpois
    implicit none
    deallocate(p)
  end subroutine exitpois
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fillps

  ! Chiel van Heerwaarden,  19 June 2007
  ! Adapted fillps for RK3 time loop


    use modfields, only : up, vp, wp, um, vm, wm
    use modglobal, only : rk3step, i1,i2,j1,kmax,k1,ih,jh, dx,dy,dzf,dzh,rdt
    use modmpi,    only : excjs
    implicit none
    real,allocatable :: pup(:,:,:), pvp(:,:,:), pwp(:,:,:)
    integer i,j,k
    real rk3coef

    allocate(pup(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(pvp(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(pwp(2-ih:i1+ih,2-jh:j1+jh,k1))

    rk3coef = rdt / (4. - dble(rk3step))

    do k=1,kmax
      do j=2,j1
        do i=2,i1
          pup(i,j,k) = up(i,j,k) + um(i,j,k) / rk3coef
          pvp(i,j,k) = vp(i,j,k) + vm(i,j,k) / rk3coef
          pwp(i,j,k) = wp(i,j,k) + wm(i,j,k) / rk3coef
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

    do j=2,j1
      do i=2,i1
        pwp(i,j,1)  = 0.
        pwp(i,j,k1) = 0.
      end do
    end do

    do k=1,kmax
      do j=2,j1
        pup(i2,j,k) = pup(2,j,k)
      end do
    end do

    call excjs( pup          , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( pvp          , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( pwp          , 2,i1,2,j1,1,k1,ih,jh)

    do k=1,kmax
      do j=2,j1
        do i=2,i1
          p(i,j,k)  =  ( pup(i+1,j,k)-pup(i,j,k) ) / dx &
                          +( pvp(i,j+1,k)-pvp(i,j,k) ) / dy &
                          +( pwp(i,j,k+1)-pwp(i,j,k) ) / dzf(k)
        end do
      end do
    end do

    deallocate( pup,pvp,pwp )

  end subroutine fillps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine solmpj(p1)
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
    use modglobal, only : imax,jmax,kmax,i1,j1,k1,kmax,isen,jtot,pi,dxi,dyi,dzi,dzf,dzh
    implicit none

    real, intent (inout) :: p1(0:i1,0:j1,0:k1)
    real, allocatable, dimension(:,:,:) :: d,p2
    real, allocatable, dimension(:,:) :: xyrt
    real, allocatable, dimension(:) :: xrt,yrt,a,b,c,FFTI,FFTJ,winew,wjnew
    real    z,ak,bk,bbk,fac
    integer jv
    integer i, j, k
  ! p and d distributed equally:
    allocate(d(imax,jmax,kmax))

  ! re-distributed p:

    allocate(p2(isen,jtot,kmax))

  ! re-distributed p1:

    allocate(xyrt(0:i1,0:jtot+1),xrt(0:i1),yrt(0:jtot+1))
    allocate(a(0:kmax+1),b(0:kmax+1),c(0:kmax+1))
    allocate(FFTI(imax),FFTJ(jtot),winew(2*imax+15),wjnew(2*jtot+15))

    call MPI_COMM_RANK( comm3d, myid, mpierr )
    call MPI_COMM_SIZE( comm3d, nprocs, mpierr )
    call rffti(imax,winew)
    call rffti(jtot,wjnew)
!     call barrou()
  !FFT  ---> I direction
    fac = 1./sqrt(imax*1.)
    do k=1,kmax
      do j=1,jmax
          do i=1,imax
            FFTI(i) =p1(i,j,k)
          end do
          call rfftf(imax,FFTI,winew)
          do i=1,imax
  ! ATT: First back to p1, then re-distribution!!!
            p1(i,j,k)=FFTI(i)*fac
          end do
      end do
    end do
    call ALL_ALL_j(p1,p2,0)
  !FFT  ---> J direction

    fac = 1./sqrt(jtot*1.)
    do i=1,isen
      do k=1,kmax
          do j=1,jtot
            FFTJ(j) =p2(i,j,k)
          end do
          call rfftf(jtot,FFTJ,wjnew)
          do j=1,jtot
  ! ATTT back to pl
            p2(i,j,k)=FFTJ(j)*fac
          end do
      end do
    end do
!     call barrou()
    call ALL_ALL_j(p1,p2,1)
!     call barrou()


  ! Generate Eigenvalues  (xrt and yrt )

  !  I --> direction


    fac = 1./(2.*imax)
    do i=3,imax,2
      xrt(i-1)=-4.*dxi*dxi*(sin(float((i-1))*pi*fac))**2
      xrt(i)  = xrt(i-1)
    end do
    xrt(1    ) = 0.
    xrt(imax ) = -4.*dxi*dxi

  !  J --> direction

    fac = 1./(2.*jtot)
    do j=3,jtot,2
      yrt(j-1)=-4.*dyi*dyi*(sin(float((j-1))*pi*fac))**2
      yrt(j  )= yrt(j-1)
    end do
    yrt(1    ) = 0.
    yrt(jtot ) = -4.*dyi*dyi

  ! Generate tridiagonal matrix
    do k=1,kmax
      a(k)=1./(dzf(k)*dzh(k))
      c(k)=1./(dzf(k)*dzh(k+1))
      b(k)=-(a(k)+c(k))
    end do
    b(1   )=b(1   )+a(1   )
    a(1   )=0.
    b(kmax)=b(kmax)+c(kmax)
    c(kmax)=0.


  ! SOLVE TRIDIAGONAL SYSTEMS WITH GAUSSIAN ELEMINATION
    do j=1,jmax
      jv = j + myid*jmax
      do i=1,imax
        xyrt(i,j)= xrt(i)+yrt(jv) !!! changed!!!
  !         xyrt(i,j)= xrt(i)+yrt(j) !!! changed!!!
        z        = 1./(b(1)+xyrt(i,j))
        d(i,j,1) = c(1)*z
        p1(i,j,1) = p1(i,j,1)*z
      end do
    end do

    do k=2,kmax-1
      do  j=1,jmax
        do  i=1,imax
          bbk      = b(k)+xyrt(i,j)
          z        = 1./(bbk-a(k)*d(i,j,k-1))
          d(i,j,k) = c(k)*z
          p1(i,j,k) = (p1(i,j,k)-a(k)*p1(i,j,k-1))*z
        end do
      end do
    end do



    ak =a(kmax)
    bk =b(kmax)
    do j=1,jmax
      do i=1,imax
  !         bbk = bk +xrt(i)+yrt(j)
        bbk = bk +xyrt(i,j)
        z        = bbk-ak*d(i,j,kmax-1)
        if(z/=0.) then
          p1(i,j,kmax) = (p1(i,j,kmax)-ak*p1(i,j,kmax-1))/z
        else
          p1(i,j,kmax) =0.
        end if
      end do
    end do
    do k=kmax-1,1,-1
      do j=1,jmax
        do i=1,imax
          p1(i,j,k) = p1(i,j,k)-d(i,j,k)*p1(i,j,k+1)
        end do
      end do
    end do
!     call barrou()


  ! MPI_ALL CALL!!!



    call ALL_ALL_j(p1,p2,0)
!     call barrou()


  ! BACKWARD FFT ---> I direction


  ! BACKWARD FFT ---> J direction

    fac = 1./sqrt(jtot*1.)
    do i=1,isen
      do k=1,kmax
        do j=1,jtot
  ! ATT, ADAPTED!!!
          FFTJ(j) =p2(i,j,k)
        end do
        call rfftb(jtot,FFTJ,wjnew)
        do j=1,jtot
  ! ATT back to p2!!!
          p2(i,j,k)=FFTJ(j)*fac
        end do
      end do
    end do
!     call barrou()
    call ALL_ALL_j(p1,p2,1)

    fac = 1./sqrt(imax*1.)
    do k=1,kmax
      do j=1,jmax
        do i=1,imax
          FFTI(i) =p1(i,j,k)
        end do
        call rfftb(imax,FFTI,winew)
        do i=1,imax
  ! ATT back to p1 !!!
          p1(i,j,k)=FFTI(i)*fac
        end do
      end do
    end do
    deallocate(d,p2,xyrt,xrt,yrt,a,b,c,FFTI,FFTJ,winew,wjnew)
!     call barrou()
    return
  end subroutine solmpj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

    use modfields, only : up, vp, wp
    use modglobal, only : i1,j1,i2,j2,kmax,k1,dx,dy,dzh
    use modmpi,    only : excj
    implicit none
    integer i,j,k


  ! Mathieu ATTTT: CHANGED!!! Loop removed!!!

  ! **  Cyclic boundary conditions **************

    do k=1,kmax
    do j=2,j1
      p(1,j,k)  = p(i1,j,k)
      p(i2,j,k) = p(2,j,k)
    end do

    end do
    call excj( p           , 1, i2, 1, j2, 0,k1)

  !*****************************************************************
  ! **  Calculate time-derivative for the velocities with known ****
  ! **  pressure gradients.  ***************************************
  !*****************************************************************

    do k=1,kmax
    do j=2,j1
    do i=2,i1

      up(i,j,k) = up(i,j,k)-(p(i,j,k)-p(i-1,j,k))/dx
      vp(i,j,k) = vp(i,j,k)-(p(i,j,k)-p(i,j-1,k))/dy

  !      up(i,j,k) = up(i,j,k)-(p(i,j,k)-p(i-1,j,k))/dx
  !     1            - 0.5 * um(i,j,k)/dt
  !      vp(i,j,k) = vp(i,j,k)-(p(i,j,k)-p(i,j-1,k))/dy
  !     1            - 0.5 * vm(i,j,k)/dt

    end do
    end do
    end do

    do k=2,kmax
    do j=2,j1
    do i=2,i1
  !      wp(i,j,k) = wp(i,j,k)+presgrad(i,j,k)
  !     1            - 0.5 * wm(i,j,k)/dt
      wp(i,j,k) = wp(i,j,k)-(p(i,j,k)-p(i,j,k-1))/dzh(k)
    end do
    end do
    end do


    return
  end subroutine tderive


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


end module modpois



