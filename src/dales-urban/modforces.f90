!> \file modforces.f90
!!  Calculates the other forces and sources in the equations.

!>
!!  Calculates the other forces and sources in the equations.
!>
!!  This includes the large scale forcings, the coriolis and the subsidence
!!  \author Jasper Tomas, TU Delft  March 31 2014
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  Only the routine 'forces' is used
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
public :: forces, coriolis, lstend,fixuinf1,fixuinf2,fixthetainf,&
          detfreestream,detfreestrtmp,initforces
contains
  subroutine initforces
  use modglobal, only : freestreamav,freestrtmpav,ifixuinf,ltempeq
  implicit none
  real freestream,freestrtmp
    
    if (ifixuinf==2) then
      call detfreestream(freestream)
      freestreamav = freestream
      if (ltempeq==.true.) then
        call detfreestrtmp(freestrtmp)
        freestrtmpav = freestrtmp
      end if
    end if
  
  end subroutine initforces

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

    !  use modglobal, only : i1,j1,kmax,dzh,dzf,grav
    use modglobal, only : ib,ie,jb,je,kb,ke,kh,dzhi,dzf,grav,lbuoyancy
    use modfields, only : u0,v0,w0,up,vp,wp,thv0h,dpdxl,dpdyl,thlp,thlpcar,thvh
    use modibmdata, only : nxwallsnorm, xwallsnorm
    use modsurfdata,only : thvs
    implicit none

    real thvsi
    integer i, j, k, n, jm, jp, km, kp


    if (lbuoyancy == .true.) then
!ILS13 replace thvsi by thvh      
! thvsi = 1./thvsi
       do k=kb+1,ke
          do j=jb,je
             do i=ib,ie
                up(i,j,k) = up(i,j,k) - dpdxl(k)
                vp(i,j,k) = vp(i,j,k) - dpdyl(k)
                !write(*,*), 'k',k
                !write(*,*), 'thvh',thvh(k)
                wp(i,j,k) = wp(i,j,k) + grav * (thv0h(i,j,k)-thvh(k))/thvh(k)
                !IS+HJ      wp(i,j,k) = wp(i,j,k)+ grav*thlsi * (thl0h(i,j,k)-thls)
             end do
          end do
       end do
    else
       do k=kb+1,ke
          do j=jb,je
             do i=ib,ie
                up(i,j,k) = up(i,j,k) - dpdxl(k)
                vp(i,j,k) = vp(i,j,k) - dpdyl(k)
                ! IS+HJ      wp(i,j,k) = wp(i,j,k)
             end do
          end do
       end do
    end if

    !     ----------------------------------------------
    !     add radiative heating to potential temperature
    !     ----------------------------------------------
    do k=kb,ke
       do j=jb,je
          do i=ib,ie
             thlp(i,j,k) = thlp(i,j,k)+thlpcar(k)
          end do
       end do
    end do

    !     --------------------------------------------
    !     special treatment for lowest full level: k=1
    !     --------------------------------------------

    do j=jb,je
       jp = j+1
       jm = j-1
       do i=ib,ie

          up(i,j,kb) = up(i,j,kb) - dpdxl(kb)

          vp(i,j,kb) = vp(i,j,kb) - dpdyl(kb)

          wp(i,j,kb) = 0.0

       end do
    end do
    !     ----------------------------------------------end i,j-loop


    return
  end subroutine forces


  subroutine detfreestream(freestream)
  use modglobal, only : ib,ie,jb,je,kb,ke,kh,dxf,xh,dt,&
                        Uinf
  use modfields, only : u0,dpdxl,dgdt,dpdx
  use modmpi, only    : myid,comm3d,mpierr,mpi_sum,my_real,nprocs
  implicit none

  real, intent(out) :: freestream

  real  utop
  integer i,j

    utop = 0.

    do j =jb,je
      do i =ib,ie
        utop = utop + 0.5*(u0(i,j,ke)+u0(i+1,j,ke))*dxf(i)
      end do
    end do
    utop = utop / ( (je-jb+1)*(xh(ie+1)-xh(ib) ) )
    call MPI_ALLREDUCE(utop,    freestream,1,MY_REAL,MPI_SUM,comm3d,mpierr)
    freestream = freestream / nprocs

  end subroutine detfreestream



  subroutine detfreestrtmp(freestrtmp)
  use modglobal, only : ib,ie,jb,je,kb,ke,kh,dxf,xh,dt,&
                        Uinf
  use modfields, only : thl0,dpdxl,dgdt,dpdx
  use modmpi, only    : myid,comm3d,mpierr,mpi_sum,my_real,nprocs
  implicit none

  real, intent(out) :: freestrtmp

  real  ttop
  integer i,j

    ttop = 0.

    do j =jb,je
      do i =ib,ie
        ttop = ttop + thl0(i,j,ke)*dxf(i)
      end do
    end do
    ttop = ttop / ( (je-jb+1)*(xh(ie+1)-xh(ib) ) )
    call MPI_ALLREDUCE(ttop,    freestrtmp,1,MY_REAL,MPI_SUM,comm3d,mpierr)
    freestrtmp = freestrtmp / nprocs

  end subroutine detfreestrtmp

  subroutine fixuinf2
  use modglobal, only : ib,ie,jb,je,kb,ke,kh,dxf,xh,dt,&
                        Uinf,ifixuinf,tscale,timee,rk3step,inletav,&
                        freestreamav,freestrtmpav,ltempeq
  use modsurfdata,only: thl_top
  use modfields, only : u0,thl0,dpdxl,dgdt,dpdx,thlsrcdt
  use modmpi, only    : myid,comm3d,mpierr,mpi_sum,my_real,nprocs
  implicit none

  real  utop,freestream,freestrtmp,rk3coef
  integer i,j

  utop = 0.

  if (ifixuinf==2 .and. rk3step==3) then
  
    call detfreestream(freestream)

    freestreamav=  freestream*dt/inletav + (1.-dt/inletav)*freestreamav

!  Write some statistics to monitoring file 
      if (myid==0) then
        open(unit=11,file='freestr.txt',position='append')
        write(11,3002) timee,freestream,freestreamav
3002    format (13(6e14.6))
        close(11)
      endif


!    dgdt =  (1./tscale) * (freestream - Uinf)
    dgdt =  (1./tscale) * (freestreamav - Uinf)            ! plus sign because dpdx is SUBTRACTED from Navier-Stokes eqs
!    dgdt =  (1./inletav) * (freestreamav - Uinf)

    if (ltempeq==.true.) then
      call detfreestrtmp(freestrtmp)
      freestrtmpav=  freestrtmp*dt/inletav + (1.-dt/inletav)*freestrtmpav
      thlsrcdt = -(1./tscale) * (freestrtmpav - thl_top)   ! minus sign because thlsr is ADDED to Navier-Stokes eqs.

      if (myid==0) then
        open(unit=11,file='theta_top.txt',position='append')
        write(11,3009) timee,freestrtmp,freestrtmpav
3009    format (13(6e20.12))
        close(11)
      endif
    end if

  end if
  end subroutine fixuinf2


  subroutine fixuinf1
  use modglobal, only : ib,ie,jb,je,kb,ke,kh,dxf,xh,dt,&
                        Uinf,ifixuinf,tscale,timee,rk3step,inletav,&
                        freestreamav
  use modfields, only : u0,dpdxl,dgdt,dpdx
  use modmpi, only    : myid,comm3d,mpierr,mpi_sum,my_real,nprocs
  implicit none

  real  utop,freestream,rk3coef
  integer i,j

  utop = 0.


  if (ifixuinf==1 .and. rk3step==3) then
    
!    rk3coef = dt / (4. - dble(rk3step))
 
    do j =jb,je
      do i =ib,ie
        utop = utop + 0.5*(u0(i,j,ke)+u0(i+1,j,ke))*dxf(i)
      end do
    end do
    utop = utop / ( (je-jb+1)*(xh(ie+1)-xh(ib) ) )
    call MPI_ALLREDUCE(utop,    freestream,1,MY_REAL,MPI_SUM,comm3d,mpierr)
    freestream = freestream / nprocs

!  Write some statistics to monitoring file 
!      if (myid==0 .and. rk3step==3) then


!    dpdxl(:) = dpdx + (1./rk3coef) * (freestream - Uinf)  
    dpdxl(:) = dpdx + (1./dt) * (freestream - Uinf)  

      if (myid==0) then
        open(unit=11,file='freestr.txt',position='append')
        write(11,3003) timee,freestream
3003    format (13(6e20.12))
        close(11)

        open(unit=11,file='dpdx___.txt',position='append')
        write(11,3002) timee,dpdxl(kb),dpdxl(kb)-dpdx
3002    format (13(6e20.12))
        close(11)
      endif
  
  end if

  end subroutine fixuinf1



  subroutine fixthetainf
  use modglobal, only : ib,ie,jb,je,kb,ke,kh,dxf,xh,dt,&
                        Uinf,ifixuinf,tscale,timee,rk3step,inletav,&
                        freestreamav,thlsrc,ltempeq
  use modfields, only : thl0
  use modmpi, only    : myid,comm3d,mpierr,mpi_sum,my_real,nprocs
  use modsurfdata, only: thl_top
  implicit none

  real  ttop,freestreamtheta,rk3coef
  integer i,j

  ttop = 0.


  if (ifixuinf==1 .and. rk3step==3 .and. ltempeq==.true.) then

!    rk3coef = dt / (4. - dble(rk3step))

    do j =jb,je
      do i =ib,ie
        ttop = ttop + thl0(i,j,ke)*dxf(i)
      end do
    end do
    ttop = ttop / ( (je-jb+1)*(xh(ie+1)-xh(ib) ) )
    call MPI_ALLREDUCE(ttop,    freestreamtheta,1,MY_REAL,MPI_SUM,comm3d,mpierr)
    freestreamtheta = freestreamtheta / nprocs


    thlsrc = -(1./dt) * (freestreamtheta - thl_top)
      if (myid==0) then
        open(unit=11,file='theta_top.txt',position='append')
        write(11,3003) timee,freestreamtheta
3003    format (13(6e20.12))
        close(11)

        open(unit=11,file='thlsrc.txt',position='append')
        write(11,3002) timee,thlsrc
3002    format (13(6e20.12))
        close(11)
      endif

  end if

  end subroutine fixthetainf




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

!  use modglobal, only : i1,j1,kmax,dzh,dzf,cu,cv,om22,om23
  use modglobal, only : ib,ie,jb,je,kb,ke,kh,dzh,dzf,cu,cv,om22,om23,lcoriol
  use modfields, only : u0,v0,w0,up,vp,wp
  implicit none

  integer i, j, k, jm, jp, km, kp

  if (lcoriol == .true.) then
  
  do k=kb+1,ke
    kp=k+1
    km=k-1
  do j=jb,je
    jp=j+1
    jm=j-1
  do i=ib,ie

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

  do j=jb,je
    jp = j+1
    jm = j-1
  do i=ib,ie

    up(i,j,kb) = up(i,j,kb)  + cv*om23 &
          +(v0(i,j,kb)+v0(i,jp,kb)+v0(i-1,j,kb)+v0(i-1,jp,kb))*om23*0.25 &
          -(w0(i,j,kb)+w0(i,j ,kb+1)+w0(i-1,j,kb+1)+w0(i-1,j ,kb))*om22*0.25

    vp(i,j,kb) = vp(i,j,kb) - cu*om23 &
          -(u0(i,j,kb)+u0(i,jm,kb)+u0(i+1,jm,kb)+u0(i+1,j,kb))*om23*0.25

    wp(i,j,kb) = 0.0

  end do
  end do
!     ----------------------------------------------end i,j-loop
  end if

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

  use modglobal, only : ib,ie,jb,je,kb,ke,kh,dzh,nsv,lmomsubs
  use modfields, only : up,vp,thlp,qtp,svp,&
                        whls, u0av,v0av,thl0av,qt0av,sv0av,&
                        dudxls,dudyls,dvdxls,dvdyls,dthldxls,dthldyls,dqtdxls,dqtdyls,dqtdtls
  implicit none

  integer k,n
  real subs_thl,subs_qt,subs_u,subs_v,subs_sv
!   if (ltimedep) then
! !     call ls
!   end if


!     1. DETERMINE LARGE SCALE TENDENCIES
!        --------------------------------

!     1.1 lowest model level above surface : only downward component

  subs_u   = 0.
  subs_v   = 0.
  subs_thl = 0.
  subs_qt  = 0.
  subs_sv  = 0.

  k = kb
  if (whls(k+1).lt.0) then !neglect effect of mean ascending on tendencies at the lowest full level
  subs_thl     = 0.5*whls(k+1)  *(thl0av(k+1)-thl0av(k))/dzh(k+1)
  subs_qt      = 0.5*whls(k+1)  *(qt0av (k+1)-qt0av(k) )/dzh(k+1)
    if(lmomsubs) then
       subs_u  = 0.5*whls(k+1)  *(u0av  (k+1)-u0av(k)  )/dzh(k+1)
       subs_v  = 0.5*whls(k+1)  *(v0av  (k+1)-v0av(k)  )/dzh(k+1)
    endif
    do n=1,nsv
      subs_sv =  0.5*whls(k+1)  *(sv0av(k+1,n)-sv0av(k,n)  )/dzh(k+1)
!      svp(2:i1,2:j1,1,n) = svp(2:i1,2:j1,1,n)-subs_sv
      svp(ib:ie,jb:je,kb,n) = svp(ib:ie,jb:je,kb,n)-subs_sv
    enddo
  endif

  thlp(ib:ie,jb:je,k) = thlp(ib:ie,jb:je,k) -u0av(k)*dthldxls(k)-v0av(k)*dthldyls(k)-subs_thl
  qtp(ib:ie,jb:je,k)  = qtp (ib:ie,jb:je,k) -u0av(k)*dqtdxls (k)-v0av(k)*dqtdyls (k)-subs_qt +dqtdtls(k)
  up(ib:ie,jb:je,k)   = up  (ib:ie,jb:je,k) -u0av(k)*dudxls(k)  -v0av(k)*dudyls  (k)-subs_u
  vp(ib:ie,jb:je,k)   = vp  (ib:ie,jb:je,k) -u0av(k)*dvdxls(k)  -v0av(k)*dvdyls  (k)-subs_v


!     1.2 other model levels twostream

  do k=kb+1,ke

    if (whls(k+1).lt.0) then   !downwind scheme for subsidence
      subs_thl    = whls(k+1) * (thl0av(k+1) - thl0av(k))/dzh(k+1)
      subs_qt     = whls(k+1) * (qt0av (k+1) - qt0av (k))/dzh(k+1)
      do n=1,nsv
        subs_sv   = whls(k+1)  *(sv0av(k+1,n) - sv0av(k,n))/dzh(k+1)
        svp(ib:ie,jb:je,k,n) = svp(ib:ie,jb:je,k,n)-subs_sv
      enddo
      if(lmomsubs) then
         subs_u   = whls(k+1) * (u0av  (k+1) - u0av  (k))/dzh(k+1)
         subs_v   = whls(k+1) * (v0av  (k+1) - v0av  (k))/dzh(k+1)
      endif
    else !downwind scheme for mean upward motions
      subs_thl    = whls(k) * (thl0av(k) - thl0av(k-1))/dzh(k)
      subs_qt     = whls(k) * (qt0av (k) - qt0av (k-1))/dzh(k)
      do n=1,nsv
        subs_sv   = whls(k) * (sv0av(k,n) - sv0av(k-1,n))/dzh(k)
        svp(ib:ie,jb:je,k,n) = svp(ib:ie,jb:je,k,n)-subs_sv
      enddo
      if(lmomsubs) then
         subs_u   = whls(k) * (u0av  (k) - u0av  (k-1))/dzh(k)
         subs_v   = whls(k) * (v0av  (k) - v0av  (k-1))/dzh(k)
      endif
    endif

    thlp(ib:ie,jb:je,k) = thlp(ib:ie,jb:je,k)-u0av(k)*dthldxls(k)-v0av(k)*dthldyls(k)-subs_thl
    qtp (ib:ie,jb:je,k) = qtp (ib:ie,jb:je,k)-u0av(k)*dqtdxls (k)-v0av(k)*dqtdyls (k)-subs_qt+dqtdtls(k)
    up  (ib:ie,jb:je,k) = up  (ib:ie,jb:je,k)-u0av(k)*dudxls  (k)-v0av(k)*dudyls  (k)-subs_u
    vp  (ib:ie,jb:je,k) = vp  (ib:ie,jb:je,k)-u0av(k)*dvdxls  (k)-v0av(k)*dvdyls  (k)-subs_v

  enddo

  return
  end subroutine lstend

end module modforces
