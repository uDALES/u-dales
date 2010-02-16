!> \file modforces.f90
!!  Calculates the other forces and sources in the equations.

!>
!!  Calculates the other forces and sources in the equations.
!>
!!  This includes the large scale forcings, the coriolis and the subsidence
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



module modforces
!Calculates additional forces and large scale tendencies
implicit none
save
private
public :: forces, coriolis, lstend,lforce_user
logical :: lforce_user = .false.
contains
  subroutine forces

!-----------------------------------------------------------------|
!                                                                 |
!      Hans Cuijpers   I.M.A.U.                                   |
!      Pier Siebesma   K.N.M.I.     06/01/1995                    |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!      Calculates all other terms in the N-S equation,            |
!      except for the diffusion and advection terms.              |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!     *forces* is called from *program*.                          |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal, only : i1,j1,kmax,dzh,dzf,grav
  use modfields, only : u0,v0,w0,up,vp,wp,thv0h,dpdxl,dpdyl
  use modsurfdata,only : thvs
  use moduser,   only : force_user
  implicit none

  integer i, j, k, jm, jp, km, kp

  if (lforce_user) call force_user

  do k=2,kmax
    kp=k+1
    km=k-1
  do j=2,j1
    jp=j+1
    jm=j-1
  do i=2,i1

    up(i,j,k) = up(i,j,k) - dpdxl(k)
    vp(i,j,k) = vp(i,j,k) - dpdyl(k)

    wp(i,j,k) = wp(i,j,k)+ grav/thvs * thv0h(i,j,k)

  end do
  end do
!     -------------------------------------------end i&j-loop
  end do
!     -------------------------------------------end k-loop

!     --------------------------------------------
!     special treatment for lowest full level: k=1
!     --------------------------------------------

  do j=2,j1
    jp = j+1
    jm = j-1
  do i=2,i1

    up(i,j,1) = up(i,j,1) - dpdxl(1)

    vp(i,j,1) = vp(i,j,1) - dpdyl(1)

    wp(i,j,1) = 0.0

  end do
  end do
!     ----------------------------------------------end i,j-loop


  return
  end subroutine forces
  subroutine coriolis

!-----------------------------------------------------------------|
!                                                                 |
!      Thijs Heus TU Delft                                        |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!      Calculates the Coriolis force.                             |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!     *coriolis* is called from *program*.                          |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal, only : i1,j1,kmax,dzh,dzf,cu,cv,om22,om23
  use modfields, only : u0,v0,w0,up,vp,wp
  implicit none

  integer i, j, k, jm, jp, km, kp


  do k=2,kmax
    kp=k+1
    km=k-1
  do j=2,j1
    jp=j+1
    jm=j-1
  do i=2,i1

    up(i,j,k) = up(i,j,k)+ cv*om23 &
          +(v0(i,j,k)+v0(i,jp,k)+v0(i-1,j,k)+v0(i-1,jp,k))*om23*0.25 &
          -(w0(i,j,k)+w0(i,j,kp)+w0(i-1,j,kp)+w0(i-1,j,k))*om22*0.25

    vp(i,j,k) = vp(i,j,k)  - cu*om23 &
          -(u0(i,j,k)+u0(i,jm,k)+u0(i+1,jm,k)+u0(i+1,j,k))*om23*0.25


    wp(i,j,k) = wp(i,j,k) + cu*om22 +( (dzf(km) * (u0(i,j,k)  + u0(i+1,j,k) )    &
                +    dzf(k)  * (u0(i,j,km) + u0(i+1,j,km))  ) / dzh(k) ) &
                * om22*0.25
  end do
  end do
!     -------------------------------------------end i&j-loop
  end do
!     -------------------------------------------end k-loop

!     --------------------------------------------
!     special treatment for lowest full level: k=1
!     --------------------------------------------

  do j=2,j1
    jp = j+1
    jm = j-1
  do i=2,i1

    up(i,j,1) = up(i,j,1)  + cv*om23 &
          +(v0(i,j,1)+v0(i,jp,1)+v0(i-1,j,1)+v0(i-1,jp,1))*om23*0.25 &
          -(w0(i,j,1)+w0(i,j ,2)+w0(i-1,j,2)+w0(i-1,j ,1))*om22*0.25

    vp(i,j,1) = vp(i,j,1) - cu*om23 &
          -(u0(i,j,1)+u0(i,jm,1)+u0(i+1,jm,1)+u0(i+1,j,1))*om23*0.25

    wp(i,j,1) = 0.0

  end do
  end do
!     ----------------------------------------------end i,j-loop


  return
  end subroutine coriolis

  subroutine lstend

!-----------------------------------------------------------------|
!                                                                 |
!*** *lstend*  calculates large-scale tendencies                  |
!                                                                 |
!      Pier Siebesma   K.N.M.I.     06/01/1995                    |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!     calculates and adds large-scale tendencies due to           |
!     large scale advection and subsidence.                       |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!             *lstend* is called from *program*.                  |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal, only : i1,j1,k1,kmax,dzh,nsv,lmomsubs
  use modfields, only : up,vp,thlp,qtp,svp,&
                        whls, u0av,v0av,thl0av,qt0av,sv0av,&
                        dudxls,dudyls,dvdxls,dvdyls,dthldxls,dthldyls,dqtdxls,dqtdyls,dqtdtls
  implicit none

  integer k,n
  real subsplus, subsmin, subs
  real,allocatable,dimension(:) :: subsu, subsv
  allocate(subsu(k1))
  allocate(subsv(k1))
  subsu=0.0
  subsv=0.0
!   if (ltimedep) then
! !     call ls
!   end if


!     1. DETERMINE LARGE SCALE TENDENCIES
!        --------------------------------

!     1.1 lowest model level above surface : only downward component

  k = 1
  subs        = 0.5*whls(k+1)  *(thl0av(k+1)-thl0av(k)  )/dzh(k+1)
  thlp(2:i1,2:j1,1) = thlp(2:i1,2:j1,1) -u0av(k)*dthldxls(k)-v0av(k)*dthldyls(k)-subs

  subs        = 0.5*whls(k+1)  *(qt0av(k+1)-qt0av(k)  )/dzh(k+1)
  qtp(2:i1,2:j1,1)  = qtp(2:i1,2:j1,1)-u0av(k)*dqtdxls(k)-v0av(k)*dqtdyls(k)-subs+dqtdtls(k)

  if(lmomsubs) subsu(k)    = 0.5*whls(k+1)  *(u0av(k+1)-u0av(k)  )/dzh(k+1)
  up(2:i1,2:j1,1)   = up(2:i1,2:j1,1) -u0av(k)*dudxls(k)-v0av(k)*dudyls(k)-subsu(k)

  if(lmomsubs) subsv(k)    = 0.5*whls(k+1)  *(v0av(k+1)-v0av(k)  )/dzh(k+1)
  vp(2:i1,2:j1,1)   = vp(2:i1,2:j1,1) -u0av(k)*dvdxls(k)-v0av(k)*dvdyls(k)-subsv(k)

  do n=1,nsv
    subs =  0.5*whls(k+1)  *(sv0av(k+1,n)-sv0av(k,n)  )/dzh(k+1)
    svp(2:i1,2:j1,1,n) = svp(2:i1,2:j1,1,n)-subs
  enddo

!     1.2 other model levels twostream

  do k=2,kmax
    subs    = whls(k+1)  *(thl0av(k+1)-thl0av(k)  )/dzh(k+1)
    thlp(2:i1,2:j1,k) = thlp(2:i1,2:j1,k)-u0av(k)*dthldxls(k)-v0av(k)*dthldyls(k)-subs

    subs    = whls(k+1)  *(qt0av(k+1) - qt0av(k)  ) /dzh(k+1)
    qtp(2:i1,2:j1,k) = qtp(2:i1,2:j1,k) -u0av(k)*dqtdxls(k)-v0av(k)*dqtdyls(k)-subs+dqtdtls(k)

    subsplus    = whls(k+1) *(u0av(k+1) - u0av(k)  )/dzh(k+1)
    subsmin     = whls(k)   *(u0av(k)   - u0av(k-1))/dzh(k)
    if(lmomsubs) subsu(k)    = (dzh(k)*subsplus + dzh(k+1)*subsmin) &
                     /(dzh(k)+dzh(k+1))
    up(2:i1,2:j1,k)   = up(2:i1,2:j1,k)-u0av(k)*dudxls(k)-v0av(k)*dudyls(k)-subsu(k)

    subsplus    = whls(k+1) *(v0av(k+1) - v0av(k)  )/dzh(k+1)
    subsmin     = whls(k)   *(v0av(k)   - v0av(k-1))/dzh(k)
    if(lmomsubs) subsv(k)    = (dzh(k)*subsplus + dzh(k+1)*subsmin) &
                     /(dzh(k)+dzh(k+1))
    vp(2:i1,2:j1,k)   = vp(2:i1,2:j1,k)-u0av(k)*dvdxls(k)-v0av(k)*dvdyls(k)-subsv(k)
    do n=1,nsv
      subs  = whls(k+1)  *(sv0av(k+1,n) - sv0av(k,n)  ) /dzh(k+1)
      svp(2:i1,2:j1,k,n) = svp(2:i1,2:j1,k,n)-subs
    enddo

  enddo

  deallocate(subsu,subsv)

  return
  end subroutine lstend

end module modforces
