!----------------------------------------------------------------------------
! This file is part of DALES.
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
! Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!----------------------------------------------------------------------------
!
!

  subroutine advecc_kappa(var,varp)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Modified by Hans Cuijpers (19 June 1995)
!      Subroutine originates from:
!                              kappa1.0a.f
!          A program for the solution of advection problems.
!               (F77-real Version 1.0a, January 1995)
!                               by
!                           Cees Beets
!                               of
!      The Institute for Marine and Atmospheric Research Utrecht
!                    e-mail: beets@fys.ruu.nl
!_______________________________________________________________________
!____________________________PURPOSE____________________________________
!     Subroutine kappa solves the advection equation:
!       c  +  (u c) + (v c) + (w c)  = 0
!        t         x       y       z
!     on structered grids.
!     For the time discretization a positive second order
!     Runge-Kutta scheme is used (ref. CWI REPORT NM-R930,
!     W.H. Hundsdorfer, B. Koren, M. van Loon and J.G. Verwer).
!     For the space discretization the positive limited
!     k=1/3 scheme is used (ref. CWI REPORT NM-R9308, B. Koren)
!_________________________ON INPUT______________________________________
!     var(i2,j2,k1):         Scalar (theta_l, q_tot or concentration)
!_________________________ON OUTPUT____________________________________
!     varp(i2,j2,k1):        Advection part of scalar
!_____________________PROGRAM SPECIFICATIONS____________________________
!     uses
!       rlim:      limiter
!     other specs
!       i/o:       none
!       language:  F77 (compilable under F90)
!       bc:        assumes periodic boundary conditions in x,y
!                  and zero flux boundary condition at the surface
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  use modglobal, only : i1,i2,ih,j1,j2,jh,k1,kmax,dxi,dyi,dzi
  use modfields, only : u0, v0, w0
  implicit none

  real,external :: rlim
  real      var(2-ih:i1+ih,2-jh:j1+jh,k1),varp(2-ih:i1+ih,2-jh:j1+jh,k1)
  real      d1,d2,cf
  integer   i,j,k
  do k=1,kmax
    do j=2,j1
      do i=2,i2
        if (u0(i,j,k)>0) then
          d1 = var(i-1,j,k)-var(i-2,j,k)
          d2 = var(i  ,j,k)-var(i-1,j,k)
          cf = var(i-1,j,k)
        else
          d1 = var(i  ,j,k)-var(i+1,j,k)
          d2 = var(i-1,j,k)-var(i  ,j,k)
          cf = var(i  ,j,k)
        end if
        cf = cf + rlim(d1,d2)
        varp(i,j,k)   = varp(i,j,k)   + cf * u0(i,j,k) * dxi
        varp(i-1,j,k) = varp(i-1,j,k) - cf * u0(i,j,k) * dxi
      end do
    end do
  end do

  do k=1,kmax
    do j=2,j2
      do i=2,i1
        if (v0(i,j,k)>0) then
          d1 = var(i,j-1,k)-var(i,j-2,k)
          d2 = var(i,j  ,k)-var(i,j-1,k)
          cf = var(i,j-1,k)
        else
          d1 = var(i,j  ,k)-var(i,j+1,k)
          d2 = var(i,j-1,k)-var(i,j  ,k)
          cf = var(i,j  ,k)
        end if
        cf = cf + rlim(d1,d2)
        varp(i,j-1,k) = varp(i,j-1,k) - cf * v0(i,j,k) * dyi
        varp(i,j,k)   = varp(i,j,k)   + cf * v0(i,j,k) * dyi
      end do
    end do
  end do

  do k=3,kmax
    do j=2,j1
      do i=2,i1
        if (w0(i,j,k)>0) then
          d1 = var(i,j,k-1)-var(i,j,k-2)
          d2 = var(i,j,k  )-var(i,j,k-1)
          cf = var(i,j,k-1)
        else
          d1 = var(i,j,k  )-var(i,j,k+1)
          d2 = var(i,j,k-1)-var(i,j,k  )
          cf = var(i,j,k  )
        end if
        cf = cf + rlim(d1,d2)
        varp(i,j,k-1) = varp(i,j,k-1) - cf * w0(i,j,k) * dzi
        varp(i,j,k)   = varp(i,j,k)   + cf * w0(i,j,k) * dzi
      end do
    end do
  end do

  do j=2,j1
    do i=2,i1
      if (w0(i,j,2)>0) then
        d1 = 0
        d2 = var(i,j,1)-var(i,j,2)
        cf = var(i,j,1)
      else
        d1 = var(i,j,2)-var(i,j,3)
        d2 = var(i-1,j,1)-var(i  ,j,2)
        cf = var(i  ,j,2)
      end if
      cf = cf + rlim(d1,d2)
      varp(i,j,1) = varp(i,j,1) - cf * w0(i,j,2) * dzi
      varp(i,j,2) = varp(i,j,2) + cf * w0(i,j,2) * dzi
    end do
  end do


  return
  end subroutine advecc_kappa

  subroutine  halflev_kappa(var,varout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax
    use modfields, only : u0, v0, w0
    implicit none
    real,external :: rlim
    real      var(2-ih:i1+ih,2-jh:j1+jh,k1),varout(2-ih:i1+ih,2-jh:j1+jh,k1)
    real      d1,d2,cf
    integer   i,j,k

    do k=3,k1
      do j=2,j1
        do i=2,i1
          if (w0(i,j,k)>=0) then
            d1 = var(i,j,k-1)-var(i,j,k-2)
            d2 = var(i,j,k  )-var(i,j,k-1)
            cf = var(i,j,k-1)
          else
            d1 = var(i,j,k  )-var(i,j,k+1)
            d2 = var(i,j,k-1)-var(i,j,k  )
            cf = var(i,j,k  )
          end if
          varout(i,j,k) = cf + rlim(d1,d2)

        end do
      end do
    end do

    do j=2,j1
      do i=2,i1
        if (w0(i,j,2)>=0) then
          d1 = 0
          d2 = var(i,j,2)-var(i,j,1)
          cf = var(i,j,1)
        else
          d1 = var(i,j,2)-var(i,j,3)
          d2 = var(i,j,1)-var(i,j,2)
          cf = var(i,j,2)
        end if
        varout(i,j,2) = cf + rlim(d1,d2)
      end do
    end do
  end subroutine halflev_kappa
!***********************************************************************
  real function rlim(d1,d2)
    use modglobal, only : eps1
    implicit none

    real ri,phir,d1,d2

    ri    = (d2+eps1)/(d1+eps1)
    phir  = max(0.,min(2.*ri,min(1./3.+2./3.*ri,2.)))
    rlim  = 0.5*phir*d1
    end function rlim



