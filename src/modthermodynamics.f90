!>\file modthermodynamics.f90
!! Do the thermodynamics

!>
!! Do the thermodynamics
!>
!! Timeseries of the most relevant parameters. Written to tmser1.expnr and tmsurf.expnr
!! If netcdf is true, this module leads the tmser.expnr.nc output
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Thijs Heus,MPI-M
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

module modthermodynamics

  implicit none
  !   private
  public :: thermodynamics,calc_halflev
  public :: lqlnr
  logical :: lqlnr    = .false. !< switch for ql calc. with Newton-Raphson (on/off)
  real, allocatable :: th0av(:)
  real :: chi_half=0.5  !< set wet, dry or intermediate (default) mixing over the cloud edge
  real, allocatable :: thv0(:,:,:)
contains

  !> Allocate and initialize arrays
  subroutine initthermodynamics
    use modglobal, only : kb,ke,kh,ib,ie,jb,je

    implicit none

    allocate(th0av(kb:ke+kh))
    allocate(thv0(ib:ie,jb:je,kb:ke+kh))
    th0av = 0.

  end subroutine initthermodynamics

  !> Do moist thermodynamics.
  !! Calculate the liquid water content, do the microphysics, calculate the mean hydrostatic pressure, calculate the fields at the half levels, and finally calculate the virtual potential temperature.
  subroutine thermodynamics
    use modglobal, only : lmoist, timee, kb, ke, kh, ib, ih, ie, jb, jh, je,rlv, cp, rslabs, rd, rv, libm, eps1
    use modfields, only : thl0,thl0h,qt0,qt0h,ql0,ql0h,presf,presh,exnf,exnh,thvh,thv0h,qt0av,ql0av,thvf,rhof,IIc,IIw,IIcs,IIws
    use modmpi,    only : slabsum,avexy_ibm,myid
!ILS13 added variables behind "exnh"
    implicit none
    integer :: k
    if (timee==0) call diagfld
    if (lmoist) then
       call thermo(thl0,qt0,ql0,presf,exnf)
    end if

    call diagfld
    call calc_halflev !calculate halflevel values of qt0 and thl0
    if (lmoist) then
       call thermo(thl0h,qt0h,ql0h,presh,exnh)
    end if
    call calthv

!ILS13 introduced from DALES4.0   13.05.2015
    thvh=0.
!    call slabsum(thvh,kb,ke+kh,thv0h(:,:,kb:ke+kh),ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh) !redefine halflevel thv using calculated thv
!    thvh = thvh/rslabs
    call avexy_ibm(thvh(kb:ke+kh),thv0h(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)

!    if (libm) then
!      call avexy_ibm(thvh(kb:ke),thv0h(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIw(ib:ie,jb:je,kb:ke),IIws(kb:ke))    
!    else
!      call slabsum(thvh,kb,ke+kh,thv0h(:,:,kb:ke+kh),ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
!     !redefine halflevel thv using calculated thv
!     thvh = thvh/rslabs
!    end if

    thvh(kb) = th0av(kb)*(1+(rv/rd-1)*qt0av(kb)-rv/rd*ql0av(kb)) ! override first level
    if (abs(thvh(kb+1))<eps1) then
      thvh(kb+1) = th0av(kb+1)*(1+(rv/rd-1)*qt0av(kb+1)-rv/rd*ql0av(kb+1)) ! override second level if all blocks at kb
    end if
!    where (thvh==0) !override slabs completely covered by blocks
!      thvh = th0av(kb)*(1+(rv/rd-1)*qt0av(kb)-rv/rd*ql0av(kb))
!    endwhere

    do k=kb,ke+kh
!    thv0(ib+ih:ie,jb+jh:je,k) = (thl0(ib+ih:ie,jb+ih:je,k)+rlv*ql0(ib+ih:ie,jb+ih:je,k)/(cp*exnf(k)))*(1+(rv/rd-1)*qt0(ib+ih:ie,jb+ih:je,k)-rv/rd*ql0(ib+ih:ie,jb+ih:je,k))
    thv0(ib:ie,jb:je,k) = (thl0(ib:ie,jb:je,k)+rlv*ql0(ib:ie,jb:je,k)/(cp*exnf(k)))*(1+(rv/rd-1)*qt0(ib:ie,jb:je,k)-rv/rd*ql0(ib:ie,jb:je,k))
    enddo
    thvf = 0.0

    !write(*,*) "thv0",thv0
!    call slabsum(thvf,kb,ke+kh,thv0,ib,ie+ih,jb,je+jh,kb,ke+kh,ib+ih,ie,jb+ih,je,kb,ke+kh)
!    call slabsum(thvf,kb,ke+kh,thv0,ib,ie,jb,je,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
    call avexy_ibm(thvf(kb:ke+kh),thv0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
!    write(*,*) 'IIc(2,2,:), myid' , IIc(12,2,:), myid

!    where (thvf==0) !override slabs completely covered by blocks
!      thvf = th0av(kb)*(1+(rv/rd-1)*qt0av(kb)-rv/rd*ql0av(kb))
!    endwhere

!    thvf = thvf/rslabs
    !write(*,*) "thvf",thvf
    !write(*,*) "exnf",exnf
   
!    do k=1,k1
!      rhof(k) = presf(k)/(rd*thvf(k)*exnf(k))
!    end do

  end subroutine thermodynamics
  !> Cleans up after the run
  subroutine exitthermodynamics
    implicit none
    deallocate(th0av)
  end subroutine exitthermodynamics

  !> Calculate thetav and dthvdz
  subroutine calthv
    use modglobal, only : lmoist,ib,ie,jb,je,kb,ke,kh,zf,zh,dzh,rlv,rd,rv,cp,eps1
    use modfields, only : thl0,thl0h,ql0,ql0h,qt0,qt0h,sv0,exnf,exnh,thv0h,dthvdz
    use modsurfdata,only : dthldz,dqtdz

    implicit none

    integer i, j, k
    real    qs
    real    c1,c2,dq,dth,dthv,temp
    real    a_dry, b_dry, a_moist, b_moist, c_liquid, epsilon, eps_I, chi_sat, chi
    real    del_thv_sat, del_thv_dry
    dthvdz = 0.0
    if (lmoist) then

       do  k=kb,ke+kh
          do  j=jb,je
             do  i=ib,ie
                thv0h(i,j,k) = (thl0h(i,j,k)+rlv*ql0h(i,j,k)/(cp*exnh(k))) &
                     *(1+(rv/rd-1)*qt0h(i,j,k)-rv/rd*ql0h(i,j,k))
             end do
          end do
       end do

       do k=kb+1,ke
          do j=jb,je
             do i=ib,ie
                !
                !         default thv jump computed unsaturated
                !
                epsilon = rd/rv
                eps_I = 1/epsilon - 1.  !cstep approx 0.608

                a_dry = 1. + eps_I * qt0(i,j,k)
                b_dry = eps_I * thl0(i,j,k)

                dth = thl0(i,j,k+1)-thl0(i,j,k-1)
                dq  = qt0(i,j,k+1)-qt0(i,j,k-1)

                del_thv_dry = a_dry   * dth + b_dry * dq

                dthv = del_thv_dry

                if  (ql0(i,j,k)> 0) then  !include moist thermodynamics
                   temp = thl0(i,j,k)*exnf(k)+(rlv/cp)*ql0(i,j,k)
                   qs   = qt0(i,j,k) - ql0(i,j,k)

                   a_moist = (1.-qt0(i,j,k)+qs/epsilon*(1.+rlv/(rv*temp))) &
                        /(1.+rlv**2*qs/(cp*rv*temp**2))
                   b_moist = a_moist*rlv/cp-temp
                   c_liquid = a_dry * rlv / cp - thl0(i,j,k) / epsilon

                   del_thv_sat = a_moist * dth + b_moist * dq

                   chi     = 2*chi_half*(zf(k) - zf(k-1))/(dzh(k)+dzh(k+1))
                   chi_sat = c_liquid * ql0(i,j,k) / (del_thv_dry - del_thv_sat)

                   if (chi < chi_sat) then  !mixed parcel is saturated
                      dthv = del_thv_sat
                   end if
                end if

                dthvdz(i,j,k) = dthv/(dzh(k+1)+dzh(k))
             end do
          end do
       end do

       do j=jb,je
          do i=ib,ie
            dthvdz(i,j,kb)=0.
          end do
       end do

    else
       !      thv0h = thl0h
       thv0h = thl0h(:,:,kb:ke+kh)
       do k=kb+1,ke
          do j=jb,je
             do i=ib,ie
                dthvdz(i,j,k) = (thl0(i,j,k+1)-thl0(i,j,k-1))/(dzh(k+1)+dzh(k))
             end do
          end do
       end do
       do  j=jb,je
          do  i=ib,ie
            dthvdz(i,j,kb)=0.
          end do
       end do
    end if

    !CvH remove WHERE
    !where (abs(dthvdz)<eps1)
    !  dthvdz = sign(eps1,dthvdz)
    !end where
    do k=kb,ke
       do j=jb,je
          do i=ib,ie
             if(abs(dthvdz(i,j,k)) < eps1) then
                dthvdz(i,j,k) = sign(eps1, dthvdz(i,j,k))
             end if
          end do
       end do
    end do



  end subroutine calthv
  !> Calculate diagnostic slab averaged fields.
  !!     Calculates slab averaged fields assuming
  !!     hydrostatic equilibrium for: u,v,theta_l,theta_v,
  !!     qt,ql,exner,pressure and the density
  !! \author      Pier Siebesma   K.N.M.I.     06/01/1995
  subroutine diagfld

    use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh,khc,nsv,zh,zf,rslabs,grav,rlv,cp,rd,rv,pref0
    use modfields, only : u0,v0,thl0,qt0,ql0,sv0,u0av,v0av,thl0av,qt0av,ql0av,sv0av, &
         presf,presh,exnf,exnh,rhof,thvf,IIc,IIcs,IIu,IIus,IIv,IIvs
    use modsurfdata,only : thls,ps
    use modmpi,    only : slabsum,myid,avexy_ibm
    implicit none

    real     tv
    integer  k, n


    !!*********************************************************
    !!  1.0   calculate average profiles of u,v,thl,qt and ql *
    !!        assuming hydrostatic equilibrium                *
    !!*********************************************************

    ! initialise local MPI arrays
    u0av = 0.0
    v0av = 0.0
    thl0av = 0.0
    th0av  = 0.0
    qt0av  = 0.0
    ql0av  = 0.0
    sv0av = 0.


    !CvH changed momentum array dimensions to same value as scalars!
    !  call slabsum(u0av  ,kb,ke+kh,u0  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
!    call slabsum(u0av  ,kb,ke+kh,u0(:,:,kb:ke+kh)  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
    call avexy_ibm(u0av(kb:ke+kh),u0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIu(ib:ie,jb:je,kb:ke+kh),IIus(kb:ke+kh),.false.)
    !  call slabsum(v0av  ,kb,ke+kh,v0  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
!    call slabsum(v0av  ,kb,ke+kh,v0(:,:,kb:ke+kh)  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
    call avexy_ibm(v0av(kb:ke+kh),v0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIv(ib:ie,jb:je,kb:ke+kh),IIvs(kb:ke+kh),.false.)
    !  call slabsum(thl0av,kb,ke+kh,thl0,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
!    call slabsum(thl0av,kb,ke+kh,thl0(:,:,kb:ke+kh),ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
    call avexy_ibm(thl0av(kb:ke+kh),thl0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)

    !write(*,*) 'thl0av(kb), thl0av(kb+1)', thl0av(kb), thl0av(kb+1)
   
    !if (IIbl == 0) then ! as lEB applies blocks to kb and masking matrices average this to zero
    !  thl0av(kb) = thl0av(kb+1)
    !end if

    !  call slabsum(qt0av ,kb,ke+kh,qt0 ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
!    call slabsum(qt0av ,kb,ke+kh,qt0(:,:,kb:ke+kh) ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
    call avexy_ibm(qt0av(kb:ke+kh),qt0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
!    call slabsum(ql0av ,kb,ke+kh,ql0 ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
    call avexy_ibm(ql0av(kb:ke+kh),ql0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)

    exnf   = 1-grav*zf/(cp*thls)
    exnh  = 1-grav*zh/(cp*thls)
    th0av  = thl0av + (rlv/cp)*ql0av/exnf

    !write(*,*) 'thl0av',thl0av
    !write(*,*) 'th0av',th0av
    

    do n=1,nsv
!       call slabsum(sv0av(kb,n),kb,ke+kh,sv0(ib-ih,jb-jh,kb,n),ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
    call avexy_ibm(sv0av(kb:ke+khc,n),sv0(ib:ie,jb:je,kb:ke+khc,n),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+khc),IIcs(kb:ke+khc),.false.)
    end do
!    sv0av = sv0av/rslabs

    !***********************************************************
    !  2.0   calculate average profile of pressure at full and *
    !        half levels, assuming hydrostatic equilibrium.    *
    !***********************************************************

    !    2.1 Use first guess of theta, then recalculate theta
    call fromztop
    exnf = (presf/pref0)**(rd/cp)
    th0av = thl0av + (rlv/cp)*ql0av/exnf


    !    2.2 Use new updated value of theta for determination of pressure

    call fromztop



    !***********************************************************
    !  3.0   Construct density profiles and exner function     *
    !       for further use in the program                     *
    !***********************************************************

    !    3.1 determine exner

    exnh(kb) = (ps/pref0)**(rd/cp)
    exnf(kb) = (presf(kb)/pref0)**(rd/cp)
    do k=kb+1,ke+kh
       exnf(k) = (presf(k)/pref0)**(rd/cp)
       exnh(k) = (presh(k)/pref0)**(rd/cp)
    end do

    thvf(kb)      = th0av(kb)*exnf(kb)*(1+(rv/rd-1)*qt0av(kb)-rv/rd*ql0av(kb))
    rhof(kb) = presf(kb)/(rd*thvf(kb))

    !    3.2 determine rho

    do k=kb+1,ke     !+kh    ?
    !   write(*,*) "exnf(k)",exnf(k)
    !   write(*,*) "th0av(k)",th0av(k)
    !   write(*,*) "qt0av(k)",qt0av(k)
    !   write(*,*) "ql0av(k)",ql0av(k)
       thvf(k) = th0av(k)*exnf(k)*(1.+(rv/rd-1)*qt0av(k)-rv/rd*ql0av(k))
       rhof(k) = presf(k)/(rd*thvf(k))
    end do
    return
  end subroutine diagfld


  !> Calculates slab averaged pressure
  !!      Input :  zf,zh,theta and qt profile
  !!      Output:  pressure profile at full and
  !!               half levels
  !!
  !!      Method: Using hydrostatic equilibrium
  !!
  !!                              -g*pref0**(rd/cp)
  !! =====>       dp**(rd/cp)/dz = --------------
  !!                                 cp*thetav
  !! \author Pier Siebesma   K.N.M.I.     06/01/1995
  subroutine fromztop

    use modglobal, only : kmax,kb,ke,kh,dzf,dzh,rv,rd,cp,tmelt,zf,grav,pref0,lEB
    use modfields, only : qt0av,ql0av,presf,presh,thvh,thvf,IIcs
    use modsurfdata,only : ps,thvs
    implicit none

    integer   k
    real      rdocp
    real,allocatable,dimension (:) :: thetah, qth, qlh

    allocate(thetah(kb:ke+kh), qth(kb:ke+kh), qlh(kb:ke+kh))
    rdocp = rd/cp

    !**************************************************
    !    1.0 Determine theta and qt at half levels    *
    !**************************************************

    do k=kb+1,ke+kh
       thetah(k) = (th0av(k)*dzf(k-1) + th0av(k-1)*dzf(k))/(2*dzh(k))
       qth   (k) = (qt0av(k)*dzf(k-1) + qt0av(k-1)*dzf(k))/(2*dzh(k))
       qlh   (k) = (ql0av(k)*dzf(k-1) + ql0av(k-1)*dzf(k))/(2*dzh(k))
    end do

    !**************************************************
    !     2.1  calculate pressures at full levels     *
    !          assuming hydrostatic equilibrium       *
    !**************************************************

    !     1: lowest level: use thvs

    thvh(kb) = thvs
    presf(kb) = ps**rdocp - &
         grav*(pref0**rdocp)*zf(kb) /(cp*thvh(kb))
    presf(kb) = presf(kb)**(1./rdocp)

    !     2: higher levels

    do k=kb+1,ke+kh 
       thvh(k)  = thetah(k)*(1+(rv/rd-1)*qth(k)-rv/rd*qlh(k))
       presf(k) = presf(k-1)**rdocp - grav*(pref0**rdocp)*dzh(k) /(cp*thvh(k))
       presf(k) = presf(k)**(1./rdocp)
    end do

    !**************************************************
    !     2.2   calculate pressures at half levels    *
    !           assuming hydrostatic equilibrium      *
    !**************************************************

    presh(kb) = ps
    thvf(kb) = th0av(kb)*(1+(rv/rd-1)*qt0av(kb)-rv/rd*ql0av(kb))
    do k=kb+1,ke+kh
       thvf(k)  = th0av(k)*(1+(rv/rd-1)*qt0av(k)-rv/rd*ql0av(k))
       presh(k) = presh(k-1)**rdocp - grav*(pref0**rdocp)*dzf(k-1) / (cp*thvf(k-1))
       presh(k) = presh(k)**(1./rdocp)
    end do

    deallocate(thetah, qth, qlh)

    return
  end subroutine fromztop

  !> Calculates liquid water content.
  !!     Given theta_l and q_tot the liquid water content
  !!     is calculated, making an "all-or-nothing" assumption.
  !!     if lfull=true   ==> ql at full levels on output
  !!     if lfull=false  ==> ql at half levels on output
  !!
  !! \author Hans Cuijpers   I.M.A.U.
  !! \author Pier Siebesma   K.N.M.I.     06/01/1995
  subroutine thermo (thl,qt,ql,pressure,exner)



    !  use modglobal, only : ih,jh,i1,j1,k1,es0,at,bt,rd,rv,rlv,cp,tmelt
    use modglobal, only : ih,jh,ib,ie,jb,je,kb,ke,kh,es0,at,bt,rd,rv,rlv,cp,tmelt
    use modsurfdata, only : thls
    implicit none

    integer i, j, k
    real tl, es, qs, qsl, b1
    !  real, intent(in)  :: qt(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh),thl(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh),exner(kb:ke+kh),pressure(kb:ke+kh)
    real, intent(in)  :: qt(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh),thl(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh),exner(kb:ke+kh),pressure(kb:ke+kh)
    real, intent(out) :: ql(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
    real :: Tnr,qsatur,Tnr_old
    integer :: niter,nitert

    if (lqlnr) then
       !mc      calculation of T with Newton-Raphson method
       !mc      first guess is Tnr=tl
       !mc
       nitert = 0
       do k=kb,ke+kh
          do j=jb,je
             do i=ib,ie

                tl  = thl(i,j,k)*exner(k)
                Tnr=tl
                Tnr_old=0.
                do while (abs(Tnr-Tnr_old)/Tnr>1e-5)
                   niter = niter+1
                   Tnr_old = Tnr
                   es    = es0*exp(at*(Tnr-tmelt)/(Tnr-bt))
                   qsatur= rd/rv*es/(pressure(k)-(1-rd/rv)*es)
                   Tnr = Tnr - (Tnr+(rlv/cp)*qsatur-tl- &
                        (rlv/cp)*qt(i,j,k))/(1+(rlv**2*qsatur)/ &
                        (rv*cp*Tnr**2))
                end do
                nitert =max(nitert,niter)
                niter = 0

                ql(i,j,k) = dim(qt(i,j,k)-qsatur,0.)

             end do
          end do
       end do
    else


       do k=kb,ke+kh
          do j=jb,je
             do i=ib,ie
                tl  = thl(i,j,k)*exner(k)
                es  = es0*exp(at*(tl-tmelt)/(tl-bt))
                qsl = rd/rv*es/(pressure(k)-(1-rd/rv)*es)
                b1  = rlv**2/(tl**2*cp*rv)
                qs  = qsl*(1.+b1*qt(i,j,k))/(1.+b1*qsl)
                ql(i,j,k) = dim(qt(i,j,k)-qs,0.)
             end do
          end do
       end do
    end if

    return
  end subroutine thermo

  !> Calculates the scalars at half levels.
  !! If the kappa advection scheme is active, interpolation needs to be done consistently.
  subroutine calc_halflev
    use modglobal, only : ib,ie,jb,je,kb,ke,kh,dzf,dzh,iadv_thl, iadv_qt, iadv_kappa
    use modfields, only : thl0,thl0h,qt0,qt0h
    use modsurfdata,only: qts,thls
    implicit none

    integer :: i,j,k


    !      do  k=kb+1,ke+kh
    do  k=kb,ke+kh
       do  j=jb,je
          do  i=ib,ie
             thl0h(i,j,k) = (thl0(i,j,k)*dzf(k-1)+thl0(i,j,k-1)*dzf(k))/(2*dzh(k))
          end do
       end do

    end do
        thl0h(ib:ie,jb:je,kb) = thls

    !      do  k=kb+1,ke+kh
    do  k=kb,ke+kh
       do  j=jb,je
          do  i=ib,ie
             qt0h(i,j,k)  = (qt0 (i,j,k)*dzf(k-1)+qt0 (i,j,k-1)*dzf(k))/(2*dzh(k))
          end do
       end do
    end do
          qt0h(ib:ie,jb:je,kb)  = qts

  end subroutine calc_halflev

end module modthermodynamics
