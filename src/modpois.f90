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

implicit none
private
public :: initpois,poisson,exitpois,p
save

  real,allocatable :: p(:,:,:)   ! difference between p at previous and new step (p = P_new - P_old)

contains
  subroutine initpois
    use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh
    implicit none

    allocate(p(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))

  end subroutine initpois

  subroutine poisson
    use modglobal, only : ib,ie,ih,kb,ke,kh,kmax,dxh,dxf,dy,dzf,dzh,linoutflow,iinletgen,ipoiss,POISS_FFT,POISS_CYC
    use modmpi, only : myid,nprocs,barrou
    implicit none
    integer ibc1,ibc2,kbc1,kbc2,ksen

    call fillps

!  ibc?=1: neumann
!  ibc?=2: periodic
!  ibc?=3: dirichlet

    select case (ipoiss)
    case (POISS_FFT)
      call solmpj(p)
    case (POISS_CYC)
      if (linoutflow ) then
        ibc1 = 1      ! inlet
        ibc2 = 3      ! outlet
      else
        ibc1 = 2
        ibc2 = 2
      endif
      kbc1 = 1
      kbc2 = 1
      ksen = kmax/nprocs
      call poisr(p,dxf,dxh,dy,dzf,dzh, &
                 ibc1,ibc2,kbc1,kbc2,ksen)
    case default
      write(0, *) "Invalid choice for Poisson solver"
       stop 1
    end select

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


    use modfields, only : up, vp, wp, um, vm, wm,u0,v0,w0
    use modglobal, only : rk3step, ib,ie,jb,je,kb,ke,ih,jh,kh, dxfi,dyi,dzfi,dt,&
                          linoutflow,libm
    use modmpi,    only : excjs
    use modboundary, only: bcpup
!    use modibm,    only : ibmnorm

    implicit none
    real,allocatable :: pup(:,:,:), pvp(:,:,:), pwp(:,:,:)
    integer i,j,k
    real rk3coef,rk3coefi

    allocate(pup(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(pvp(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(pwp(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))

    rk3coef = dt / (4. - dble(rk3step))
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
          p(i,j,k)  =  ( pup(i+1,j,k)-pup(i,j,k) ) * dxfi(i) &         ! see equation 5.72
                          +( pvp(i,j+1,k)-pvp(i,j,k) ) * dyi &
                          +( pwp(i,j,k+1)-pwp(i,j,k) ) * dzfi(k)
        end do
      end do
    end do

    deallocate( pup,pvp,pwp )

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

    use modfields, only : up, vp, wp, pres0, IIc, IIcs
    use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh,dxhi,dyi,dzhi,linoutflow,rslabs
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

    do k=kb,ke
    do j=jb,je
    do i=ib,ie
      vp(i,j,k) = vp(i,j,k)-(p(i,j,k)-p(i,j-1,k))*dyi
    end do
    end do
    end do

    if (linoutflow ) then
      do k=kb,ke
      do j=jb,je
      do i=ib,ie+1
        up(i,j,k) = up(i,j,k)-(p(i,j,k)-p(i-1,j,k))*dxhi(i)           ! see equation 5.82 (u is computed from the mass conservation)
      end do
      end do
      end do
    else
      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        up(i,j,k) = up(i,j,k)-(p(i,j,k)-p(i-1,j,k))*dxhi(i)           ! see equation 5.82 (u is computed from the mass conservation)
      end do
      end do
      end do
    endif

    do k=kb+1,ke
    do j=jb,je
    do i=ib,ie
      wp(i,j,k) = wp(i,j,k)-(p(i,j,k)-p(i,j,k-1))*dzhi(k)
    end do
    end do
    end do

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

    do k=kb-1,ke+1
    do j=jb-1,je+1
    do i=ib-1,ie+1
      pres0(i,j,k)=pres0(i,j,k)+p(i,j,k)-pijk ! update of the pressure: P_new = P_old + p
    enddo
    enddo
    enddo

    return
  end subroutine tderive

!ils13, 13.08.18: currently unused, not callled
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

! subroutine poisr

      SUBROUTINE poisr(rhs,dx,dxh,dy,dz,dzh, &
                       ibc1,ibc2,kbc1,kbc2,ksen)
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
  use modglobal, only : imax,isen,jmax,jsen,jtot,kmax,poisrcheck
  use modfields, only : worksave
  use modmpi,    only : myid,comm3d,mpierr,my_real,nprocs, barrou,MPI_SUM
      implicit none

! authors: m.j.b.m. pourquie, b.j. boersma

      integer i, j, k, ksen
      real rhs(0:imax+1,0:jmax+1,0:kmax+1),dx(0:IMAX+1),dxh(1:IMAX+1),dy
      real, allocatable, dimension(:,:,:) ::  rhs2
      real, allocatable, dimension(:) ::  work
      integer  ier, iperio, kperio
      integer  ibc1,ibc2,kbc1,kbc2
      real dz(0:kmax+1),dzh(1:kmax+1),pi       !   dz: kb-1:ke+1,  dzh: kb:ke+1
      real a(imax),b(imax),c(imax),bin(imax)
      real az(kmax),bz(kmax),cz(kmax)
      real yrt(jtot)
      real, allocatable, dimension(:,:) ::  vfftj
      real, allocatable, dimension(:,:) ::  y
      real wj(jtot+15)
      real   angle, tst
      integer ipos, jv
      real  suml, sum
!
! test write
!     write(6,*)'POISR, kmax, ksen, nprocs,jmax,jtot ',kmax, ksen, nprocs ,jmax,jtot
      allocate(rhs2(imax,jtot,ksen))
      allocate(work(2*imax*jmax*kmax))
      allocate(vfftj(imax*ksen,jtot))
      allocate(y(imax,kmax))

      pi=4.*atan(1.)
      do i=1,imax
         a(i) =  1./(dx(i)*dxh(i))
         c(i) =  1./(dx(i)*dxh(i+1))
         b(i) =  - (a(i) + c(i))
      enddo

      if((ibc1).eq.1)then
! Neumann
         b(1)    = b(1)    + a(1)
      elseif(ibc1.eq.2)then
! periodic
         b(1)    = b(1)    
      elseif(ibc1.eq.3)then
! Dir
         b(1)    = b(1)    - a(1)
      endif

      if((ibc2).eq.1)then
! Neumann
         b(imax) = b(imax) + c(imax)
      elseif((ibc2).eq.2)then
         b(imax) = b(imax) 
      elseif((ibc2).eq.3)then
!         b(imax) = b(imax) - c(kmax)   ! Jasper T. : bug? Should be b(imax) - c(imax)?
         b(imax) = b(imax) - c(imax) 
      endif
      

      if(ibc1.ne.2)then
      c(imax) = 0.
      a(1)    = 0.
      endif
!     do i=1,imax
!        write(6,*)'i, abc(i)',i, a(i),b(i),c(i)
!     enddo
      
!
! fill coefficients in k-direction

!      do k=1,kmax                          ! Mathieu's version
!         az(k) =  1./(dz(k)*dzh(k-1))
!         cz(k) =  1./(dz(k)*dzh(k))
!         bz(k) =  - (az(k) + cz(k))
!      enddo

      do k=1,kmax
         az(k) =  1./(dz(k)*dzh(k))
         cz(k) =  1./(dz(k)*dzh(k+1))
         bz(k) =  - (az(k) + cz(k))
      enddo

!
! BC:
!
!periodic     bz(1)    = bz(1) - az(1)
!periodic      az(1)    = 0.
!periodic     bz(kmax) = bz(kmax) - az(kmax)
!periodic      az(kmax) = 0.
      if((kbc1).eq.1)then
! Neumann
         bz(1)    = bz(1)    + az(1)
      elseif(kbc1.eq.2)then
! periodic
         bz(1)    = bz(1)    
      endif

      if((kbc2).eq.1)then
! Neumann
         bz(kmax) = bz(kmax) + cz(kmax)
      elseif((kbc2).eq.2)then
         bz(kmax) = bz(kmax) 
      elseif((kbc2).eq.3)then
!
! p = 0.

         bz(kmax) = bz(kmax)  - cz(kmax)
      endif
      if(kbc1.ne.2)then
      cz(kmax) = 0.
      az(1)    = 0.
      endif
!     do k=1,kmax
!        write(6,*)'k, abc(k)',k, az(k),bz(k),cz(k)
!     enddo
!
! initialise for FFT

! PAR
!     call vrffti(jmax,wj)
      call vrffti(jtot,wj)
!         

      yrt(1)=0.
! PAR
!     yrt(jmax)=-4./(dy*dy)
      yrt(jtot)=-4./(dy*dy)
      do j=3,jtot,2
!
! 2.*(cos2(alpha) -1) = 2.*(-2.*sin(alpha)**2

      yrt(j-1)=(-4./(dy*dy))*(sin(float((j-1))*pi/(2.*jtot)))**2
      yrt(j)= yrt(j-1)
      enddo 

!     help = rhs
!
! sum check (comment out if not deeded)
!
!!      suml = 0.
!!      do k=1,kmax
!!         do j=1,jmax
!!         do i=1,imax
!!         suml = suml + rhs(i,j,k)*dx(i)*dy*dz(k)
!!         enddo
!!         enddo
!!      enddo
!!      call MPI_ALLREDUCE(suml, sum, 1, MY_REAL, MPI_SUM, comm3d,mpierr)
!!      if(myid.eq.0) write(6,*)'solver sum = ', sum
!
! end sum check
!

      call barrou()
      call ALL_ALL_j2(rhs,rhs2,0,ksen)
      call barrou()

      do k=1,ksen
      do i=1,imax
      ipos=(k-1)*imax+i
      do j=1,jtot
      vfftj(ipos,j)=rhs2(i,j,k)
!     suml = suml + rhs2(i,j,k)*dy*dz(k)*dx(i)
      enddo
      enddo
      enddo
!     call MPI_ALLREDUCE(suml, sum, 1, MY_REAL, MPI_SUM, comm3d,mpierr)

!     if(myid.eq.0) write(6,*)'solver sum = ', sum
!
      call vrfftf(imax*ksen,jtot,vfftj,rhs2,imax*ksen,wj)
!
!      goto 1234
      do k=1,ksen
      do i=1,imax
      ipos=(k-1)*imax+i
      do j=1,jtot
      rhs2(i,j,k) = vfftj(ipos,j)
      enddo
      enddo
      enddo
      call barrou()
      call ALL_ALL_j2(rhs,rhs2,1,ksen)
      call barrou()
!     call sumchk3(rhs,imax,jmax,kmax,8,1)

! PAR CHECK!!!
      do j=1,jmax       ! begin loop over angles
!     do j=1,jtot       ! begin loop over angles

!
! add proper part to diagonal

      jv = j + myid*jmax
      do i=1,imax
! PAR CHECK!!
      bin(i)=b(i) + yrt(jv)!!!!!!!!!/(Rp(i)**2)
      enddo

!ATTT      do k=1,ksen
      do k=1,kmax
      do i=1,imax
         ipos=(k-1)*imax+i
!PAR CHECK
         y(i,k) = rhs(i,j,k)!vfftj(ipos,j)
      enddo
      enddo
      iperio=1
      kperio=1
      if(ibc1.eq.2)iperio=0 ! i periodic
      if(kbc1.eq.2)kperio=0 ! k periodic
      if(poisrcheck.eq.0) then
      poisrcheck=1
      CALL BLKTRI(0,kperio,kmax,az,bz,cz,iperio,imax,a,bin,c,imax,y &
!
!                   ^ 0 for periodic BC
           ,ier,work)
!     write(6,*)'ier = ', ier, iperio, kperio
!     if(ier.ne.0)stop 'IER'
      worksave = work
      write(6,*) 'First time step in POISR, poisrcheck=', poisrcheck 
      end if
      CALL BLKTRI(1,kperio,kmax,az,bz,cz,iperio,imax,a,bin,c,imax,y &
           ,ier,worksave)
!
!     write(6,*)'myid, ier = ', myid, ier, j
!      if(ier.ne.0)stop 'IER'

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
    use modglobal, only : imax,jmax,kmax,isen,jtot,pi,dyi,dzf,dzh, dxfi, kb, ke, kh
!, rhoa
!    use modfields, only : rhobf, rhobh
    implicit none
    real :: dxi
    integer :: i1, j1, k1
!   real, intent(inout), dimension(:,:,:) :: p1 
    real p1(0:imax+1,0:jmax+1,0:kmax+1)

    real, allocatable, dimension(:,:,:) :: d,p2
    real, allocatable, dimension(:,:,:) :: xyzrt
    real, allocatable, dimension(:) :: xrt,yrt,a,b,c,FFTI,FFTJ,winew,wjnew, rhobf, rhobh
    real    z,ak,bk,bbk,fac
    integer jv
    integer i, j, k
    !real dzl(ke+kh-(kb-kh)),dzhl(ke+kh-(kb-kh))
    real dzl(0:kmax+1),dzhl(1:kmax+1)

    dxi = dxfi(1)
    dzl(0:kmax+1) = dzf(kb-kh:ke+kh)
    dzhl(1:kmax+1) = dzh(kb:ke+kh)

    i1 = imax+1
    j1 = jmax+1
    k1 = kmax+1
 !   allocate(p1(0:i1,0:j1,0:k1))
  ! p and d distributed equally:
    allocate(d(imax,jmax,kmax))

  ! re-distributed p:

    allocate(p2(isen,jtot,kmax))

  ! re-distributed p1:

    allocate(rhobf(1:kmax), rhobh(1:kmax+1))
    allocate(xyzrt(0:i1,0:j1,0:k1),xrt(0:i1),yrt(0:jtot+1))
    allocate(a(0:kmax+1),b(0:kmax+1),c(0:kmax+1))
    allocate(FFTI(imax),FFTJ(jtot),winew(2*imax+15),wjnew(2*jtot+15))

    rhobf=1; rhobh = 1;

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

    ! NOTE -- dzfm dzh are defined from 0 -- hence ..-1
    do k=1,kmax
      ! SB fixed the coefficients
      a(k)=rhobh(k)/(dzl(k)*dzhl(k))
      c(k)=rhobh(k+1)/(dzl(k)*dzhl(k+1))
      b(k)=-(a(k)+c(k))
    end do
   b(1   )=b(1)+a(1)
    a(1   )=0.
    b(kmax)=b(kmax)+c(kmax)
    c(kmax)=0.

    do k=1,kmax
    do j=1,jmax
    jv = j + myid*jmax
    do i=1,imax
      xyzrt(i,j,k)= rhobf(k)*(xrt(i)+yrt(jv)) !!! LH
    end do
    end do
    end do

  ! SOLVE TRIDIAGONAL SYSTEMS WITH GAUSSIAN ELEMINATION
    do j=1,jmax
      jv = j + myid*jmax
      do i=1,imax
        z        = 1./(b(1)+xyzrt(i,j,1))
        d(i,j,1) = c(1)*z
        p1(i,j,1) = p1(i,j,1)*z
      end do
    end do

    do k=2,kmax-1
      do  j=1,jmax
      jv = j + myid*jmax
        do  i=1,imax
          bbk      = b(k)+xyzrt(i,j,k)
          z        = 1./(bbk-a(k)*d(i,j,k-1))
          d(i,j,k) = c(k)*z
          p1(i,j,k) = (p1(i,j,k)-a(k)*p1(i,j,k-1))*z
        end do
      end do
    end do



    ak =a(kmax)
    bk =b(kmax)
    do j=1,jmax
      jv = j + myid*jmax
      do i=1,imax
        bbk = bk +xyzrt(i,j,kmax)
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
   deallocate(d,p2,xyzrt,xrt,yrt,a,b,c,FFTI,FFTJ,winew,wjnew)
!     call barrou()
    return
  end subroutine solmpj


end module modpois
