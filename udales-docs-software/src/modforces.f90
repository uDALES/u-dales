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
            detfreestream,detfreestrtmp,nudge,&
            masscorr,uoutletarea,voutletarea,fluidvolume,calcfluidvolumes
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

    !  use modglobal, only : i1,j1,kmax,dzh,dzf,grav
    use modglobal, only : ib,ie,jb,je,kb,ke,kh,dzhi,dzf,grav,lbuoyancy
    use modfields, only : u0,v0,w0,up,vp,wp,thv0h,dpdxl,dpdyl,thlp,thlpcar,thvh
    use modibmdata, only : nxwallsnorm, xwallsnorm
    use modsurfdata,only : thvs
    use modmpi, only : myid
    implicit none

    real thvsi
    integer i, j, k, n, jm, jp, km, kp


    if (lbuoyancy ) then
    !ILS13 replace thvsi by thvh
    ! thvsi = 1./thvsi

 
       !write(*,*) 'thvh',thvh
       do k=kb+1,ke
          do j=jb,je
             do i=ib,ie
                up(i,j,k) = up(i,j,k) - dpdxl(k)
                vp(i,j,k) = vp(i,j,k) - dpdyl(k)
                wp(i,j,k) = wp(i,j,k) + grav * (thv0h(i,j,k)-thvh(k))/thvh(k)
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
                          Uinf,lvinf,dy
    use modfields, only : u0,dpdxl,dgdt,dpdx,v0
    use modmpi, only    : myid,comm3d,mpierr,mpi_sum,my_real,nprocs
    implicit none
    real, intent(out) :: freestream

    real  utop,vtop,dum
    integer i,j

    if (lvinf) then
        vtop = 0.
        do j =jb,je
          do i =ib,ie
            vtop = vtop + 0.5*(v0(i,j,ke)+u0(i,j+1,ke))*dy
          end do
        end do
        vtop = vtop / ( (je-jb+1)*(xh(ie+1)-xh(ib) ) )
        call MPI_ALLREDUCE(vtop,   dum,1,MY_REAL,MPI_SUM,comm3d,mpierr)
        freestream = dum / nprocs
    else
        utop = 0.
        do j =jb,je
          do i =ib,ie
            !dum=0.5*(u0(i,j,ke)+u0(i+1,j,ke))*dxf(i)
            !utop = utop + dum
            utop = utop + 0.5*(u0(i,j,ke)+u0(i+1,j,ke))*dxf(i)
          end do
        end do
        utop = utop / ( (je-jb+1)*(xh(ie+1)-xh(ib) ) )
        call MPI_ALLREDUCE(utop,   dum,1,MY_REAL,MPI_SUM,comm3d,mpierr)
        freestream = dum / nprocs
    !write(*,*) "myid,utop,dum,freestream",myid,utop,dum,freestream
    end if

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

    if ((ifixuinf==2) .and. (rk3step==3)) then
      call detfreestream(freestream)

      freestreamav=  freestream*dt/inletav + (1.-dt/inletav)*freestreamav

      ! Write some statistics to monitoring file
      ! if (myid==0) then
      !   open(unit=11,file='freestr.txt',position='append')
      !   write(11,3002) timee,freestream,freestreamav
      !   3002      format (13(6e14.6))
      !   close(11)
      ! endif


      !    dgdt =  (1./tscale) * (freestream - Uinf)
      !    dgdt =  (1./dt) * (freestreamav - Uinf)
          dgdt =  (1./tscale) * (freestreamav - Uinf)            ! plus sign because dpdx is SUBTRACTED from Navier-Stokes eqs
      !    dgdt =  (1./inletav) * (freestreamav - Uinf)

      !    if (ltempeq) then  !tg3315 commented
      !      call detfreestrtmp(freestrtmp)
      !      freestrtmpav=  freestrtmp*dt/inletav + (1.-dt/inletav)*freestrtmpav
      !      thlsrcdt = -(1./tscale) * (freestrtmpav - thl_top)   ! minus sign because thlsr is ADDED to Navier-Stokes eqs.

      !      if (myid==0) then
      !        open(unit=11,file='theta_top.txt',position='append')
      !        write(11,3009) timee,freestrtmp,freestrtmpav
      !3009    format (13(6e20.12))
      !        close(11)
      !      endif
      !    end if

    end if
  end subroutine fixuinf2

  subroutine fixuinf1
    use modglobal, only : ib,ie,jb,je,kb,ke,kh,dxf,xh,dt,&
                          Uinf,Vinf,ifixuinf,tscale,timee,rk3step,inletav,&
                          freestreamav,lvinf
    use modfields, only : u0,dpdxl,dgdt,dpdx,up,vp
    use modmpi, only    : myid,comm3d,mpierr,mpi_sum,my_real,nprocs
    implicit none

    real  utop,freestream,rk3coef
    integer i,j,k

    utop = 0.


    if ((ifixuinf==1) .and. (rk3step==3)) then

      ! rk3coef = dt / (4. - dble(rk3step))

      ! do j =jb,je
      !   do i =ib,ie
      !     utop = utop + 0.5*(u0(i,j,ke)+u0(i+1,j,ke))*dxf(i)
      !   end do
      ! end do
      ! utop = utop / ( (je-jb+1)*(xh(ie+1)-xh(ib) ) )
      ! call MPI_ALLREDUCE(utop,    freestream,1,MY_REAL,MPI_SUM,comm3d,mpierr)
      ! freestream = freestream / nprocs

      ! Write some statistics to monitoring file
      ! if (myid==0 .and. rk3step==3) then


      ! ! dpdxl(:) = dpdx + (1./rk3coef) * (freestream - Uinf)
      ! dpdxl(:) = dpdx + (1./dt) * (freestream - Uinf)
      call detfreestream(freestream)
      ! write(*,*) "freestream",freestream
      if (lvinf) then
        do k=kb,ke
          do i=ib,ie
            do j=jb,je
              vp(i,j,k) = vp(i,j,k) - (1./dt) * (freestream - Vinf)
            enddo
          enddo
        enddo
      else
        do k=kb,ke
          do j=jb,je
            do i=ib,ie
              up(i,j,k) = up(i,j,k) - (1./dt) * (freestream - Uinf)
            enddo
          enddo
        enddo
      endif
      ! if (myid==0) then
      !   write(*,*), "freestream", freestream
      !   write(*,*), "Uinf", Uinf
      !   open(unit=11,file='freestr.txt',position='append')
      !   write(11,3003) timee,freestream
      !   3003    format (13(6e20.12))
      !   close(11)

      !   open(unit=11,file='dpdx___.txt',position='append')
      !   write(11,3002) timee,dpdxl(kb),dpdxl(kb)-dpdx
      !   3002    format (13(6e20.12))
      !   close(11)
      ! endif
 
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

    ! if (ifixuinf==1 .and. rk3step==3 .and. ltempeq) then !tg3315 commented

    !   rk3coef = dt / (4. - dble(rk3step))

    !   do j =jb,je
    !     do i =ib,ie
    !       ttop = ttop + thl0(i,j,ke)*dxf(i)
    !     end do
    !   end do
    !   ttop = ttop / ( (je-jb+1)*(xh(ie+1)-xh(ib) ) )
    !   call MPI_ALLREDUCE(ttop,    freestreamtheta,1,MY_REAL,MPI_SUM,comm3d,mpierr)
    !   freestreamtheta = freestreamtheta / nprocs


    !   thlsrc = -(1./dt) * (freestreamtheta - thl_top)
    !     if (myid==0) then
    !       open(unit=11,file='theta_top.txt',position='append')
    !       write(11,3003) timee,freestreamtheta
    !       3003    format (13(6e20.12))
    !       close(11)

    !       open(unit=11,file='thlsrc.txt',position='append')
    !       write(11,3002) timee,thlsrc
    !       3002    format (13(6e20.12))
    !       close(11)
    !     endif

    ! end if

  end subroutine fixthetainf

  subroutine masscorr
    !> correct the velocities to get prescribed flow rate

    use modglobal, only : ib,ie,jb,je,ih,jh,kb,ke,dzf,dxf,dy,dt,rk3step,&
                          uflowrate,vflowrate,linoutflow,&
                          luoutflowr,lvoutflowr,luvolflowr,lvvolflowr
    use modfields, only : um,up,vm,vp,uout,uouttot,udef,vout,vouttot,vdef,&
                          uoutarea,voutarea,fluidvol,IIu,IIv
    use modmpi,    only : myid,comm3d,mpierr,nprocs,MY_REAL,sumy_ibm

    real, dimension(ib:ie, kb:ke) :: uvol
    real, dimension(ib:ie, kb:ke) :: uvolold
    real, dimension(ib:ie, kb:ke) :: vvol
    real, dimension(ib:ie, kb:ke) :: vvolold
    real, dimension(kb:ke)        :: uoutold
    real, dimension(kb:ke)        :: voutold
    real                          rk3coef,rk3coefi,&
                                  uoutflow,voutflow,&
                                  uflowrateold,vflowrateold
    integer                       i,j,k

    if ((.not.linoutflow) .and. (luoutflowr)) then
      rk3coef = dt / (4. - dble(rk3step))
      rk3coefi = 1 / rk3coef

      udef = 0.
      uout = 0.
      uoutflow = 0.
      uoutold = 0.

      ! integrate u fixed at outlet ie along y
      call sumy_ibm(uout,up(ie,jb:je,kb:ke)*dy,ie,ie,jb,je,kb,ke,IIu(ie,jb:je,kb:ke))  ! u tendency at previous time step
      call sumy_ibm(uoutold,um(ie,jb:je,kb:ke)*dy,ie,ie,jb,je,kb,ke,IIu(ie,jb:je,kb:ke))  ! u at previous time step

      ! integrate u in z
      do k=kb,ke
        uout(k) = rk3coef*uout(k)*dzf(k)
        uoutold(k) = uoutold(k)*dzf(k)
      end do
      uoutflow = sum(uout(kb:ke))
      uflowrateold = sum(uoutold(kb:ke))

      ! average over outflow area
      uoutflow = uoutflow/uoutarea
      uflowrateold = uflowrateold/uoutarea

      ! flow correction to match outflow rate
      udef = uflowrate - (uoutflow + uflowrateold)
      do k = kb,ke
        do j = jb,je
            do i = ib,ie
              up(i,j,k) = up(i,j,k) + udef*rk3coefi
            end do
        end do
      end do

      ! bss116 calculate uouttot which is used in modboundary.
      ! this really should be in the routine directly!
      uouttot = sum(uout(kb:ke))  ! mass flow rate at outlet

    elseif ((.not.linoutflow) .and. (luvolflowr)) then
      rk3coef = dt / (4. - dble(rk3step))
      rk3coefi = 1 / rk3coef

      udef = 0.
      uout = 0.
      uoutflow = 0.
      uoutold = 0.
      uvol = 0.
      uvolold = 0.

      ! integrate u in y
      call sumy_ibm(uvol,up(ib:ie,jb:je,kb:ke)*dy,ib,ie,jb,je,kb,ke,IIu(ib:ie,jb:je,kb:ke))  ! u tendency at previous time step
      call sumy_ibm(uvolold,um(ib:ie,jb:je,kb:ke)*dy,ib,ie,jb,je,kb,ke,IIu(ib:ie,jb:je,kb:ke))  ! u at previous time step

      ! integrate u in x
      do k=kb,ke
        uout(k) = sum(uvol(ib:ie,k)*dxf(ib:ie))
        uoutold(k) = sum(uvolold(ib:ie,k)*dxf(ib:ie))
      end do

      ! integrate u in z
      do k=kb,ke
        uout(k) = rk3coef*uout(k)*dzf(k)
        uoutold(k) = uoutold(k)*dzf(k)
      end do
      uoutflow = sum(uout(kb:ke))
      uflowrateold = sum(uoutold(kb:ke))

      ! average over fluid volume
      uoutflow = uoutflow/fluidvol
      uflowrateold = uflowrateold/fluidvol

      ! flow correction to match outflow rate
      udef = uflowrate - (uoutflow + uflowrateold)

      do k = kb,ke
        do j = jb,je
            do i = ib,ie
              up(i,j,k) = up(i,j,k) + udef*rk3coefi
            end do
        end do
      end do

    end if

    if ((.not.linoutflow) .and. (lvoutflowr)) then
      rk3coef = dt / (4. - dble(rk3step))
      rk3coefi = 1 / rk3coef

      vdef = 0.
      vout = 0.
      voutflow = 0.
      voutold = 0.

      ! integrate v fixed at outlet je along x
      if (myid==nprocs-1) then
        do k=kb,ke
          vout(k) = sum(vp(ib:ie,je,k)*IIv(ib:ie,je,k)*dxf(ib:ie))  ! v tendency at previous time step
          voutold(k) = sum(vm(ib:ie,je,k)*IIv(ib:ie,je,k)*dxf(ib:ie))  ! v at previous time step
        end do
      end if

      call MPI_BCAST(vout,ke-kb+1,MY_REAL ,nprocs-1,comm3d,mpierr)
      call MPI_BCAST(voutold,ke-kb+1,MY_REAL ,nprocs-1,comm3d,mpierr)

      ! integrate v in z
      do k=kb,ke
        vout(k) = rk3coef*vout(k)*dzf(k)
        voutold(k) = voutold(k)*dzf(k)
      end do
      voutflow = sum(vout(kb:ke))
      vflowrateold = sum(voutold(kb:ke))

      ! average over outflow area
      voutflow = voutflow/voutarea
      vflowrateold = vflowrateold/voutarea

      ! flow correction to match outflow rate
      vdef = vflowrate - (voutflow + vflowrateold)
      do k = kb,ke
        do j = jb,je
            do i = ib,ie
              vp(i,j,k) = vp(i,j,k) + vdef*rk3coefi
            end do
        end do
      end do

    elseif ((.not.linoutflow) .and. (lvvolflowr)) then
      rk3coef = dt / (4. - dble(rk3step))
      rk3coefi = 1 / rk3coef

      vdef = 0.
      vout = 0.
      voutflow = 0.
      voutold = 0.
      vvol = 0.
      vvolold = 0.

      ! integrate v in y
      call sumy_ibm(vvol,vp(ib:ie,jb:je,kb:ke)*dy,ib,ie,jb,je,kb,ke,IIv(ib:ie,jb:je,kb:ke))  ! v tendency at previous time step
      call sumy_ibm(vvolold,vm(ib:ie,jb:je,kb:ke)*dy,ib,ie,jb,je,kb,ke,IIv(ib:ie,jb:je,kb:ke))  ! v at previous time step

      ! integrate v in x
      do k=kb,ke
        vout(k) = sum(vvol(ib:ie,k)*dxf(ib:ie))
        voutold(k) = sum(vvolold(ib:ie,k)*dxf(ib:ie))
      end do

      ! integrate v in z
      do k=kb,ke
        vout(k) = rk3coef*vout(k)*dzf(k)
        voutold(k) = voutold(k)*dzf(k)
      end do
      voutflow = sum(vout(kb:ke))
      vflowrateold = sum(voutold(kb:ke))

      ! average over fluid volume
      voutflow = voutflow/fluidvol
      vflowrateold = vflowrateold/fluidvol

      ! flow correction to match outflow rate
      vdef = vflowrate - (voutflow + vflowrateold)
      do k = kb,ke
        do j = jb,je
            do i = ib,ie
              vp(i,j,k) = vp(i,j,k) + vdef*rk3coefi
            end do
        end do
      end do

    end if

  end subroutine masscorr

  subroutine uoutletarea(area)
    ! calculates outlet area of domain for u-velocity excluding blocks

    use modglobal, only   : ib,ie,jb,je,kb,ke,dy,dzf
    use modfields, only   : IIc
    use modmpi, only      : sumy_ibm

    implicit none
    real, intent(out)       :: area
    real, dimension(kb:ke)  :: sumy
    integer                    k

    sumy = 0.
    ! integrate fluid area at outflow plane in y
    call sumy_ibm(sumy,IIc(ie,jb:je,kb:ke)*dy,ie,ie,jb,je,kb,ke,IIc(ie,jb:je,kb:ke))

    ! integrate fluid area at outflow plane in z
    do k=kb,ke
      sumy(k) = sumy(k)*dzf(k)
    end do
    area = sum(sumy(kb:ke))

  end subroutine uoutletarea

  subroutine voutletarea(area)
    ! calculates outlet area of domain for v-velocity excluding blocks

    use modglobal, only : ib,ie,jb,je,kb,ke,dxf,dzf
    use modfields, only : IIc
    use modmpi,    only : myid,comm3d,mpierr,nprocs,MY_REAL

    implicit none
    real, intent(out)       :: area
    real, dimension(kb:ke)  :: sumx
    integer                    k

    sumx = 0.
    ! integrate fluid area at outflow plane in x
    if (myid==nprocs-1) then
      do k=kb,ke
        sumx(k) = sum(IIc(ib:ie,je,k)*dxf(ib:ie))
      end do
    end if

    call MPI_BCAST(sumx,ke-kb+1,MY_REAL,nprocs-1,comm3d,mpierr)

    ! integrate fluid area at outflow plane in z
    do k=kb,ke
      sumx(k) = sumx(k)*dzf(k)
    end do
    area = sum(sumx(kb:ke))

  end subroutine voutletarea

  subroutine fluidvolume(volume)
    ! calculates fluid volume of domain excluding blocks

    use modglobal, only   : ib,ie,jb,je,kb,ke,dy,dxf,dzf
    use modfields, only   : IIc
    use modmpi, only      : sumy_ibm

    implicit none
    real, intent(out)             :: volume
    real, dimension(ib:ie,kb:ke)  :: sumy
    real, dimension(kb:ke)        :: sumxy
    integer                          k

    sumy = 0.
    sumxy = 0.

    ! integrate fluid volume in y
    call sumy_ibm(sumy,IIc(ib:ie,jb:je,kb:ke)*dy,ib,ie,jb,je,kb,ke,IIc(ib:ie,jb:je,kb:ke))

    ! integrate fluid area in x
    do k=kb,ke
      sumxy(k) = sum(sumy(ib:ie,k)*dxf(ib:ie))
    end do

    ! integrate fluid area in z
    volume = sum(sumxy(kb:ke)*dzf(kb:ke))

  end subroutine fluidvolume

  subroutine calcfluidvolumes
    !> calculates fluid volume and outlet areas, excluding blocks
    !> and saves it to variables from modfields

    use modfields, only : uoutarea, voutarea, fluidvol
    implicit none
    real :: volume

    ! calculate outlet area
    call uoutletarea(volume)
    uoutarea = volume
    ! calculate outlet area
    call voutletarea(volume)
    voutarea = volume
    ! calculate fluid volume
    call fluidvolume(volume)
    fluidvol = volume

  end subroutine calcfluidvolumes

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
    !     *coriolis* is called from *program*.                        |
    !                                                                 |
    !-----------------------------------------------------------------|

    ! use modglobal, only : i1,j1,kmax,dzh,dzf,om22,om23
    use modglobal, only : ib,ie,jb,je,kb,ke,kh,dzh,dzf,om22,om23,lcoriol,lprofforc,timee
    use modfields, only : u0,v0,w0,up,vp,wp,ug,vg
    use modmpi, only : myid
    implicit none

    integer i, j, k, jm, jp, km, kp
    real, dimension(kb:ke+kh) :: ugg
    real om23g

    if (lcoriol ) then
      ! if (myid==0) then
      !   write(*,*) "up before coriol",up(3,3,ke)
      ! end if
      do k=kb+1,ke
        kp=k+1
        km=k-1
        do j=jb,je
          jp=j+1
          jm=j-1
          do i=ib,ie

            up(i,j,k) = up(i,j,k)  &
                  +((v0(i,j,k)+v0(i,jp,k)+v0(i-1,j,k)+v0(i-1,jp,k))*om23*0.25) &
                  -((w0(i,j,k)+w0(i,j,kp)+w0(i-1,j,kp)+w0(i-1,j,k))*om22*0.25)

            vp(i,j,k) = vp(i,j,k)  &
                  -((u0(i,j,k)+u0(i,jm,k)+u0(i+1,jm,k)+u0(i+1,j,k))*om23*0.25)


            wp(i,j,k) = wp(i,j,k) +(( (dzf(km) * (u0(i,j,k)  + u0(i+1,j,k) )    &
                        +    dzf(k)  * (u0(i,j,km) + u0(i+1,j,km))  ) / dzh(k) ) &
                        * om22*0.25)

          end do
        end do
        ! -------------------------------------------end i&j-loop
      end do
      ! -------------------------------------------end k-loop

      ! --------------------------------------------
      ! special treatment for lowest full level: k=1
      ! --------------------------------------------

      do j=jb,je
        jp = j+1
        jm = j-1
        do i=ib,ie

          up(i,j,kb) = up(i,j,kb)  &
                +(v0(i,j,kb)+v0(i,jp,kb)+v0(i-1,j,kb)+v0(i-1,jp,kb))*om23*0.25 &
                -(w0(i,j,kb)+w0(i,j ,kb+1)+w0(i-1,j,kb+1)+w0(i-1,j ,kb))*om22*0.25

          vp(i,j,kb) = vp(i,j,kb) &
                -(u0(i,j,kb)+u0(i,jm,kb)+u0(i+1,jm,kb)+u0(i+1,j,kb))*om23*0.25

          wp(i,j,kb) = 0.0

        end do
      end do
      ! ----------------------------------------------end i,j-loop
      ! if (myid==0) then
      !   write(*,*) "up after coriol",up(3,3,ke)
      ! end if

    elseif (lprofforc) then

      ugg(:) = ug(:)
      om23g = om23

      do k=kb+1,ke
        do j=jb,je
          do i=ib,ie

            up(i,j,k) = up(i,j,k) + om23g*(ugg(k) - u0(i,j,k))

          enddo
        enddo
      enddo

      ! --------------------------------------------
      ! special treatment for lowest full level: k=1
      ! --------------------------------------------

      do j=jb,je
        jp = j+1
        jm = j-1
        do i=ib,ie

          up(i,j,kb) = up(i,j,kb) + om23g*(ugg(kb) - u0(i,j,kb))

        enddo
      enddo
      ! if (myid==0) then
      !   write(*,*) "up after profforc",up(3,3,ke)
      ! end if

    endif !lcoriol and lprofforc

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
    use modmpi, only: myid
    implicit none

    integer k,n
    real subs_thl,subs_qt,subs_u,subs_v,subs_sv

    ! if (ltimedep) then
    !   ! call ls
    ! end if
    ! if (myid==0) then
    !   write(*,*) "up before lstend",up(3,3,ke)
    ! end if


    ! 1. DETERMINE LARGE SCALE TENDENCIES
    !    --------------------------------

    ! 1.1 lowest model level above surface : only downward component

    subs_u   = 0.
    subs_v   = 0.
    subs_thl = 0.
    subs_qt  = 0.
    subs_sv  = 0.

    k = kb
    if (whls(k+1).lt.0) then !neglect effect of mean ascending on tendencies at the lowest full level
      subs_thl     = whls(k+1)  *(thl0av(k+1)-thl0av(k))/dzh(k+1) ! tg3315 ils13 bss116 31/07/18 Dales 4.0 multiplies these by 0.5. To reduce subsidence towards the ground? Have removed
      subs_qt      = whls(k+1)  *(qt0av (k+1)-qt0av(k) )/dzh(k+1)
      if(lmomsubs) then
        subs_u  = whls(k+1)  *(u0av  (k+1)-u0av(k)  )/dzh(k+1)
        subs_v  = whls(k+1)  *(v0av  (k+1)-v0av(k)  )/dzh(k+1)
      endif
      do n=1,nsv
        subs_sv =  whls(k+1)  *(sv0av(k+1,n)-sv0av(k,n)  )/dzh(k+1)
        ! svp(2:i1,2:j1,1,n) = svp(2:i1,2:j1,1,n)-subs_sv
        svp(ib:ie,jb:je,kb,n) = svp(ib:ie,jb:je,kb,n)-subs_sv
      enddo
    endif

    thlp(ib:ie,jb:je,k) = thlp(ib:ie,jb:je,k) -u0av(k)*dthldxls(k)-v0av(k)*dthldyls(k)-subs_thl
    qtp(ib:ie,jb:je,k)  = qtp (ib:ie,jb:je,k) -u0av(k)*dqtdxls (k)-v0av(k)*dqtdyls (k)-subs_qt +dqtdtls(k)
    up(ib:ie,jb:je,k)   = up  (ib:ie,jb:je,k) -u0av(k)*dudxls(k)  -v0av(k)*dudyls  (k)-subs_u
    vp(ib:ie,jb:je,k)   = vp  (ib:ie,jb:je,k) -u0av(k)*dvdxls(k)  -v0av(k)*dvdyls  (k)-subs_v


    ! 1.2 other model levels twostream
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

  subroutine nudge
    use modglobal,  only : kb,ke,lmoist,ltempeq,lnudge,tnudge,nnudge,numol,nsv
    use modfields,  only : thlp,qtp,svp,sv0av,thl0av,qt0av
    use modmpi,     only : myid
    implicit none
    integer :: k
    real :: numoli

    numoli = 1/numol

    if (lnudge .eqv. .false.) return

    if (nsv>0) then
      do k=ke-nnudge,ke
        svp(:,:,k,1) = svp(:,:,k,1) - ( sv0av(k,1) - 0. ) / (tnudge/2 + (ke-k)*tnudge/nnudge)
      end do
    end if

    if (ltempeq) then
      do k=ke-nnudge,ke
        thlp(:,:,k) = thlp(:,:,k) - ( thl0av(k) - 288 ) / (tnudge/2 + (ke-k)*tnudge/nnudge)
      end do
    end if !ltempeq

    if (lmoist) then
      do k=ke-nnudge,ke
        qtp(:,:,k) = qtp(:,:,k) - ( qt0av(k) - 0. ) / (tnudge/2 +(ke-k)*tnudge/nnudge)
      end do
    end if !lmoist

  end subroutine nudge

end module modforces
