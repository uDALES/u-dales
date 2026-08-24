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
#if defined(_GPU)
  use cudafor
#endif

  implicit none
  !   private
  public :: thermodynamics,calc_halflev,thermodynamics_step
  public :: lqlnr
#if defined(_GPU) && defined(UDALES_DEBUG)
  public :: thermodynamics_device, diagfld_device, calc_halflev_device, &
            calthv_device, thermo_device, avexy_thv0_device
#endif
  logical :: lqlnr    = .false. !< switch for ql calc. with Newton-Raphson (on/off)
  real, allocatable :: th0av(:)
  real :: chi_half=0.5  !< set wet, dry or intermediate (default) mixing over the cloud edge
  real, allocatable :: thv0(:,:,:)

  !> Grid-only combinations the per-cell loops used to re-derive.
  !!
  !! None of the three has any time dependence - each is the vertical grid and
  !! at most one constant - and each was being formed inside a triple loop
  !! over the whole domain, on every stage. They are derived once here for the
  !! host and the device both, which is the only way the two can be relied on
  !! to hold the same number.
  !!
  !! Each is the host's original expression, evaluated in the host's order, so
  !! every value is what the loop used to compute bit for bit. That rules out
  !! the obvious further step: folding twodzh and dzhpair into reciprocals
  !! would turn a division into a multiplication and change the last bit of
  !! every half level and every buoyancy gradient in the model.
  real, allocatable :: twodzh(:)   !< 2*dzh(k), the half-level interpolation divisor
  real, allocatable :: dzhpair(:)  !< dzh(k+1)+dzh(k), the dthvdz divisor
  real, allocatable :: chi_thv(:)  !< the mixing fraction calthv's saturated branch tests
#if defined(_GPU)
  real, device, allocatable :: twodzh_d(:), dzhpair_d(:), chi_thv_d(:)
#endif
contains

  !> Allocate and initialize arrays
  subroutine initthermodynamics
    use modglobal, only : kb,ke,kh,ib,ie,jb,je,zf,dzh
#if defined(_GPU)
    use modglobal, only : lmoist
#endif

    implicit none
    integer :: k

    allocate(th0av(kb:ke+kh))
    allocate(thv0(ib:ie,jb:je,kb:ke+kh))
    th0av = 0.

#if defined(_GPU)
    ! The Newton-Raphson saturation solve iterates until convergence per cell,
    ! so every thread of a warp would run to the slowest cell's iteration
    ! count. It defaults off and no case in the repository turns it on, so
    ! rather than write a kernel nothing exercises, say so here - at
    ! initialisation, where the namelist has just been read and the run has
    ! produced no output yet.
    if (lqlnr) then
       write(*,*) 'lqlnr has no GPU implementation: the Newton-Raphson &
                  &saturation solve is host-only. Run with lqlnr = .false.'
       error stop 1
    end if
#endif

    ! The grid-only divisors. See the declarations for why they are formed
    ! from the original expressions rather than as reciprocals.
    allocate(twodzh(kb:ke+kh))
    do k = kb, ke+kh
       twodzh(k) = 2*dzh(k)
    end do

    allocate(dzhpair(kb+1:ke))
    allocate(chi_thv(kb+1:ke))
    do k = kb+1, ke
       dzhpair(k) = dzh(k+1)+dzh(k)
       chi_thv(k) = 2*chi_half*(zf(k) - zf(k-1))/(dzh(k)+dzh(k+1))
    end do

#if defined(_GPU)
    allocate(twodzh_d(kb:ke+kh));  twodzh_d  = twodzh
    allocate(dzhpair_d(kb+1:ke));  dzhpair_d = dzhpair
    ! chi_thv has one reader, calthv's saturated branch, which a run without
    ! moisture never reaches.
    if (lmoist) then
       allocate(chi_thv_d(kb+1:ke)); chi_thv_d = chi_thv
    end if
#endif

  end subroutine initthermodynamics

  !> Liquid water content and the diagnostic fields, wherever the fields live.
  !!
  !! Only for the call inside the time loop. readinitfiles names
  !! thermodynamics directly, because it runs before initCUDA and there is
  !! nothing on the device to derive from - and that host call is also what
  !! seeds every device array the loop path reads before it first writes.
  subroutine thermodynamics_step
    implicit none
#if defined(_GPU)
    call thermodynamics_device
#else
    call thermodynamics
#endif
  end subroutine thermodynamics_step

  !> Do moist thermodynamics.
  !! Calculate the liquid water content, do the microphysics, calculate the mean hydrostatic pressure, calculate the fields at the half levels, and finally calculate the virtual potential temperature.
  subroutine thermodynamics
    use modglobal, only : lmoist, timee, kb, ke, kh, ib, ie, jb, je,rlv, cp, rd, rv, eps1
    use modfields, only : thl0,thl0h,qt0,qt0h,ql0,ql0h,presf,presh,exnf,exnh,thvh,thv0h,qt0av,ql0av,thvf,IIc,IIw,IIcs,IIws
    use modmpi,    only : slabsum,avexy_ibm
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
    call avexy_ibm(thvh(kb:ke+kh),thv0h(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)

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
    call avexy_ibm(thvf(kb:ke+kh),thv0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
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
    deallocate(twodzh, dzhpair, chi_thv)
#if defined(_GPU)
    deallocate(twodzh_d, dzhpair_d)
    if (allocated(chi_thv_d)) deallocate(chi_thv_d)
#endif
  end subroutine exitthermodynamics

  !> Calculate thetav and dthvdz
  subroutine calthv
    use modglobal, only : lmoist,ib,ie,jb,je,kb,ke,kh,rlv,rd,rv,cp,eps1
    use modfields, only : thl0,thl0h,ql0,ql0h,qt0,qt0h,exnf,exnh,thv0h,dthvdz

    implicit none

    integer i, j, k
    real    qs
    real    dq,dth,dthv,temp
    real    a_dry, b_dry, a_moist, b_moist, c_liquid, epsilon, eps_I, chi_sat, chi
    real    del_thv_sat, del_thv_dry
    real    rvord, rvordm1

    ! Lifted out of the loops below, where they were two divisions per cell of
    ! quantities that depend on nothing. Same expressions in the same order,
    ! so the same values.
    epsilon = rd/rv
    eps_I   = 1/epsilon - 1.   !cstep approx 0.608
    rvord   = rv/rd
    rvordm1 = rv/rd - 1.

    dthvdz = 0.0
    if (lmoist) then

       do  k=kb,ke+kh
          do  j=jb,je
             do  i=ib,ie
                thv0h(i,j,k) = (thl0h(i,j,k)+rlv*ql0h(i,j,k)/(cp*exnh(k))) &
                     *(1+rvordm1*qt0h(i,j,k)-rvord*ql0h(i,j,k))
             end do
          end do
       end do

       do k=kb+1,ke
          do j=jb,je
             do i=ib,ie
                !
                !         default thv jump computed unsaturated
                !
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

                   chi     = chi_thv(k)
                   chi_sat = c_liquid * ql0(i,j,k) / (del_thv_dry - del_thv_sat)

                   if (chi < chi_sat) then  !mixed parcel is saturated
                      dthv = del_thv_sat
                   end if
                end if

                dthvdz(i,j,k) = dthv/dzhpair(k)
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
                dthvdz(i,j,k) = (thl0(i,j,k+1)-thl0(i,j,k-1))/dzhpair(k)
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

    use modglobal, only : ib,ie,jb,je,kb,ke,kh,khc,nsv,zh,zf,grav,rlv,cp,rd,rv,pref0
    use modfields, only : u0,v0,thl0,qt0,ql0,sv0,u0av,v0av,thl0av,qt0av,ql0av,sv0av, &
         presf,presh,exnf,exnh,rhof,thvf,IIc,IIcs,IIu,IIus,IIv,IIvs
    use modsurfdata,only : thls,ps
    use modmpi,    only : slabsum,avexy_ibm
    implicit none

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
    call avexy_ibm(u0av(kb:ke+kh),u0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,kh,IIu(ib:ie,jb:je,kb:ke+kh),IIus(kb:ke+kh),.false.)
    !  call slabsum(v0av  ,kb,ke+kh,v0  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
!    call slabsum(v0av  ,kb,ke+kh,v0(:,:,kb:ke+kh)  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
    call avexy_ibm(v0av(kb:ke+kh),v0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,kh,IIv(ib:ie,jb:je,kb:ke+kh),IIvs(kb:ke+kh),.false.)
    !  call slabsum(thl0av,kb,ke+kh,thl0,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
!    call slabsum(thl0av,kb,ke+kh,thl0(:,:,kb:ke+kh),ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
    call avexy_ibm(thl0av(kb:ke+kh),thl0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)

    !write(*,*) 'thl0av(kb), thl0av(kb+1)', thl0av(kb), thl0av(kb+1)
   
    !if (IIbl == 0) then ! as lEB applies blocks to kb and masking matrices average this to zero
    !  thl0av(kb) = thl0av(kb+1)
    !end if

    !  call slabsum(qt0av ,kb,ke+kh,qt0 ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
!    call slabsum(qt0av ,kb,ke+kh,qt0(:,:,kb:ke+kh) ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
    call avexy_ibm(qt0av(kb:ke+kh),qt0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
!    call slabsum(ql0av ,kb,ke+kh,ql0 ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
    call avexy_ibm(ql0av(kb:ke+kh),ql0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)

    exnf   = 1-grav*zf/(cp*thls)
    exnh  = 1-grav*zh/(cp*thls)
    th0av  = thl0av + (rlv/cp)*ql0av/exnf

    !write(*,*) 'thl0av',thl0av
    !write(*,*) 'th0av',th0av
    

    do n=1,nsv
!       call slabsum(sv0av(kb,n),kb,ke+kh,sv0(ib-ih,jb-jh,kb,n),ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
    call avexy_ibm(sv0av(kb:ke+khc,n),sv0(ib:ie,jb:je,kb:ke+khc,n),ib,ie,jb,je,kb,ke,kh,IIc(ib:ie,jb:je,kb:ke+khc),IIcs(kb:ke+khc),.false.)
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

    use modglobal, only : kb,ke,kh,dzf,dzh,rv,rd,cp,zf,grav,pref0
    use modfields, only : qt0av,ql0av,presf,presh,thvh,thvf
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
                
                !! X. Long: This is a fix to tackle incorrect thls input. The reason why tl is going less than 100K (unphysical) 
                !! is that it is calculated from thl which dose not change over time here.Probably this is happening at
                !'calc_halflev and call thermo(thl0h,qt0h,ql0h,presh,exnh)'. This problem should be fixed later.  
                if (tl<100.0) then 
                    tl=100.0
                end if
                
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
    use modglobal, only : ib,ie,jb,je,kb,ke,kh,dzf
    use modfields, only : thl0,thl0h,qt0,qt0h
    use modsurfdata,only: qts,thls
    implicit none

    integer :: i,j,k


    !      do  k=kb+1,ke+kh
    do  k=kb,ke+kh
       do  j=jb,je
          do  i=ib,ie
             thl0h(i,j,k) = (thl0(i,j,k)*dzf(k-1)+thl0(i,j,k-1)*dzf(k))/twodzh(k)
          end do
       end do

    end do
        thl0h(ib:ie,jb:je,kb) = thls

    !      do  k=kb+1,ke+kh
    do  k=kb,ke+kh
       do  j=jb,je
          do  i=ib,ie
             qt0h(i,j,k)  = (qt0 (i,j,k)*dzf(k-1)+qt0 (i,j,k-1)*dzf(k))/twodzh(k)
          end do
       end do
    end do
          qt0h(ib:ie,jb:je,kb)  = qts

  end subroutine calc_halflev

#if defined(_GPU)
  !> Do moist thermodynamics on the device.
  !!
  !! The same six steps as thermodynamics above, in the same order, and the
  !! split between what moves and what stays is one line: anything a column
  !! deep stays on the host. fromztop's hydrostatic recurrence is sequential
  !! in k, and the exner, density and slab-average profiles beside it are
  !! ke-kb+2 elements each - there is nothing there for a kernel to do, and
  !! reusing the host arithmetic verbatim is also what keeps the two paths
  !! from drifting. What moves is the 3D work, ql and the half levels and
  !! thv0h and dthvdz, plus the six slab reductions that feed the column,
  !! because those read fields the device owns.
  !!
  !! Everything downstream of the first diagfld is skipped when ltempeq is
  !! off, rather than run. Not an optimisation: with no temperature equation
  !! nothing writes thl0 at all - tstep_integrate, advection and every branch
  !! of the boundary conditions that touches it are themselves under ltempeq -
  !! so thl0h, thv0h, dthvdz, thvh and thvf cannot move from the values the
  !! host thermodynamics in readinitfiles produced, and initCUDA carries those
  !! onto the device. Running the kernels would recompute the same numbers
  !! from device arrays that a dry run does not allocate.
  subroutine thermodynamics_device
    use modglobal, only : lmoist, ltempeq, timee, kb, ke, kh, rd, rv, eps1
    use modfields, only : thvh, thvf, qt0av, ql0av, IIws
    use modcuda,   only : thl0_d, qt0_d, ql0_d, thl0h_d, qt0h_d, ql0h_d, &
                          thv0h_d, presf_d, presh_d, exnf_d, exnh_d, IIw_d, &
                          avexy_ibm_device
    implicit none

    real, dimension(kb:ke+kh) :: thvh_keep, thvf_keep

    ! diagfld's hydrostatic pass writes thvh and thvf on its way to presf, and
    ! the tail of this routine then replaces both with slab averages of thv0h
    ! and thv0 - which are the kernels a run without a temperature equation
    ! skips. So the two profiles have to be carried across the call rather
    ! than left at the intermediate the pass happens to leave in them.
    ! Neither can have moved: nothing writes thl0, qt0 or ql0 in that
    ! configuration. thvh in particular is read - it is the buoyancy reference
    ! modforces divides by - so this is what makes the skip exact rather than
    ! nearly right.
    if (.not. ltempeq) then
       thvh_keep = thvh
       thvf_keep = thvf
    end if

    if (timee==0) call diagfld_device
    if (lmoist) then
       call thermo_device(thl0_d, qt0_d, ql0_d, presf_d, exnf_d)
    end if

    call diagfld_device

    if (.not. ltempeq) then
       thvh = thvh_keep
       thvf = thvf_keep
       return
    end if

    call calc_halflev_device
    if (lmoist) then
       call thermo_device(thl0h_d, qt0h_d, ql0h_d, presh_d, exnh_d)
    end if
    call calthv_device

    thvh = 0.
    call avexy_ibm_device(thvh(kb:ke+kh), thv0h_d, lbound(thv0h_d,3), &
                          IIw_d, IIws(kb:ke+kh), .false.)

    thvh(kb) = th0av(kb)*(1+(rv/rd-1)*qt0av(kb)-rv/rd*ql0av(kb)) ! override first level
    if (abs(thvh(kb+1))<eps1) then
      thvh(kb+1) = th0av(kb+1)*(1+(rv/rd-1)*qt0av(kb+1)-rv/rd*ql0av(kb+1)) ! override second level if all blocks at kb
    end if

    thvf = 0.0
    call avexy_thv0_device(thvf(kb:ke+kh))

  end subroutine thermodynamics_device

  !> Slab averages on the device, hydrostatic column on the host.
  !!
  !! The host routine zeroes all seven profiles before filling them, and
  !! avexy_ibm assigns rather than accumulates, so six of those zeros change
  !! nothing. The seventh does: sv0av's tail level, ke+khc, is past what the
  !! reduction writes and keeps the zero. Zeroing the other six here would be
  !! worse than useless - the temperature profiles are not recomputed when
  !! ltempeq is off, and zeroing them would throw away what initialisation
  !! left in favour of nothing.
  subroutine diagfld_device
    use modglobal, only : ib,ie,jb,je,kb,ke,kh,khc,nsv,zh,zf,grav,rlv,cp,rd,rv,pref0, &
                          ltempeq,lmoist
    use modfields, only : u0av,v0av,thl0av,qt0av,ql0av,sv0av, &
                          presf,presh,exnf,exnh,rhof,thvf,IIcs,IIus,IIvs
    use modsurfdata,only : thls,ps
    use modcuda,   only : u0_d,v0_d,thl0_d,qt0_d,ql0_d,sv0_d,IIc_d,IIu_d,IIv_d, &
                          avexy_ibm_device,avexy_ibm_c_device,updateThermoProfilesDevice
    implicit none

    integer  k, n

    sv0av = 0.

    call avexy_ibm_device(u0av(kb:ke+kh), u0_d, lbound(u0_d,3), IIu_d, IIus(kb:ke+kh), .false.)
    call avexy_ibm_device(v0av(kb:ke+kh), v0_d, lbound(v0_d,3), IIv_d, IIvs(kb:ke+kh), .false.)

    ! thl0, qt0 and ql0 have no writer without a temperature equation, so
    ! their averages are the ones initialisation left. See thermodynamics_device.
    if (ltempeq) then
      call avexy_ibm_device(thl0av(kb:ke+kh), thl0_d, lbound(thl0_d,3), IIc_d, IIcs(kb:ke+kh), .false.)
      call avexy_ibm_device(qt0av(kb:ke+kh),  qt0_d,  lbound(qt0_d,3),  IIc_d, IIcs(kb:ke+kh), .false.)
      call avexy_ibm_device(ql0av(kb:ke+kh),  ql0_d,  lbound(ql0_d,3),  IIc_d, IIcs(kb:ke+kh), .false.)
    end if

    do n=1,nsv
      call avexy_ibm_c_device(sv0av(kb:ke+khc,n), sv0_d(:,:,:,n), lbound(sv0_d,3), &
                              IIc_d, IIcs(kb:ke+khc), .false.)
    end do

    exnf   = 1-grav*zf/(cp*thls)
    exnh  = 1-grav*zh/(cp*thls)
    th0av  = thl0av + (rlv/cp)*ql0av/exnf

    !    2.1 Use first guess of theta, then recalculate theta
    call fromztop
    exnf = (presf/pref0)**(rd/cp)
    th0av = thl0av + (rlv/cp)*ql0av/exnf

    !    2.2 Use new updated value of theta for determination of pressure
    call fromztop

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
    do k=kb+1,ke
       thvf(k) = th0av(k)*exnf(k)*(1.+(rv/rd-1)*qt0av(k)-rv/rd*ql0av(k))
       rhof(k) = presf(k)/(rd*thvf(k))
    end do

    ! The four columns the saturation and virtual-temperature kernels read.
    call updateThermoProfilesDevice

  end subroutine diagfld_device

  !> Liquid water content on the device.
  !!
  !! The dummy declarations are copied from thermo above character for
  !! character, and that is load-bearing rather than tidy. ql is declared one
  !! k level shorter than qt and thl, so the whole-array ql0 the first call
  !! passes lands on the dummy one level low and the routine writes
  !! ql0(:,:,kb-kh:ke) where the loop reads kb:ke+kh. That shift is what the
  !! host does - the commented-out declarations above thermo show qt and thl
  !! being widened to kb-kh and ql being left behind - and it reaches ql0av,
  !! calthv and the fielddump output. Declaring the device dummy the obvious
  !! way would silently disagree with the CPU build on every moist case; this
  !! way the two are wrong together and stay comparable. ql0h, whose lower
  !! bound is kb, is not affected either way.
  !!
  !! lqlnr's Newton-Raphson branch is absent; initthermodynamics refuses the
  !! run rather than let it reach here.
  subroutine thermo_device(thl, qt, ql, pressure, exner)
    use modglobal, only : ih,jh,ib,ie,jb,je,kb,ke,kh,es0,at,bt,rd,rv,rlv,cp,tmelt
    implicit none

    real, device, intent(in)  :: qt(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh), &
                                 thl(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh), &
                                 exner(kb:ke+kh), pressure(kb:ke+kh)
    real, device, intent(out) :: ql(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)

    integer :: i, j, k
    real    :: tl, es, qs, qsl, b1
    real    :: rdorv, omrdorv, rlv2

    ! Formed once and reaching the kernel as scalars. collapse(3) gives every
    ! cell its own thread, so an expression written inside the loop body is
    ! evaluated once per cell however invariant it is - and these three
    ! carried two divisions between them.
    rdorv   = rd/rv
    omrdorv = 1-rd/rv
    rlv2    = rlv**2

    ! firstprivate spelled out for the three hoisted constants. They are the
    ! default for a scalar read but not written inside the region, so this
    ! changes nothing - it is there so that the next person to add a hoisted
    ! scalar sees where it belongs. Putting one in private instead hands every
    ! thread an uninitialised copy, silently and only at -O2.
    !$acc parallel loop collapse(3) default(present) &
    !$acc   private(tl,es,qs,qsl,b1) firstprivate(rdorv,omrdorv,rlv2)
    do k=kb,ke+kh
       do j=jb,je
          do i=ib,ie
             tl  = thl(i,j,k)*exner(k)

             !! X. Long: This is a fix to tackle incorrect thls input. The reason why tl is going less than 100K (unphysical)
             !! is that it is calculated from thl which dose not change over time here.
             if (tl<100.0) then
                 tl=100.0
             end if

             es  = es0*exp(at*(tl-tmelt)/(tl-bt))
             qsl = rdorv*es/(pressure(k)-omrdorv*es)
             b1  = rlv2/(tl**2*cp*rv)
             qs  = qsl*(1.+b1*qt(i,j,k))/(1.+b1*qsl)
             ! dim(x,0.) with x = qt-qs, spelled as the max the standard
             ! defines it to be, so no device runtime has to implement DIM.
             ql(i,j,k) = max(qt(i,j,k)-qs, 0.)
          end do
       end do
    end do
    !$acc end parallel loop

  end subroutine thermo_device

  !> Scalars at the half levels, on the device.
  subroutine calc_halflev_device
    use modglobal, only : ib,ie,jb,je,kb,ke,kh
    use modcuda,   only : thl0_d,thl0h_d,qt0_d,qt0h_d,dzf_d
    use modsurfdata,only: qts,thls
    implicit none

    integer :: i,j,k
    real    :: thls_l, qts_l, dzfm, dzfp, dzi

    ! The surface values reach the kernel as scalars rather than as module
    ! variables, which is what every other device routine here does.
    thls_l = thls
    qts_l  = qts

    ! One kernel for both fields, where the host has two loops. Same traffic
    ! in the prognostic fields, but the three grid values per level are
    ! fetched once instead of twice, and twodzh_d is the divisor formed at
    ! initialisation rather than 2*dzh(k) re-derived by every thread.
    !$acc parallel loop collapse(3) default(present) &
    !$acc   private(dzfm,dzfp,dzi) firstprivate(thls_l,qts_l)
    do  k=kb,ke+kh
       do  j=jb,je
          do  i=ib,ie
             dzfm = dzf_d(k-1)
             dzfp = dzf_d(k)
             dzi  = twodzh_d(k)
             thl0h_d(i,j,k) = (thl0_d(i,j,k)*dzfm+thl0_d(i,j,k-1)*dzfp)/dzi
             qt0h_d(i,j,k)  = (qt0_d(i,j,k)*dzfm+qt0_d(i,j,k-1)*dzfp)/dzi
          end do
       end do
    end do
    !$acc end parallel loop

    !$acc parallel loop collapse(2) default(present)
    do j=jb,je
       do i=ib,ie
          thl0h_d(i,j,kb) = thls_l
          qt0h_d(i,j,kb)  = qts_l
       end do
    end do
    !$acc end parallel loop

  end subroutine calc_halflev_device

  !> thetav at the half levels and its vertical gradient, on the device.
  subroutine calthv_device
    use modglobal, only : lmoist,ib,ie,jb,je,kb,ke,kh,ih,jh,rlv,rd,rv,cp,eps1
    use modcuda,   only : thl0_d,thl0h_d,ql0_d,ql0h_d,qt0_d,qt0h_d, &
                          exnf_d,exnh_d,thv0h_d,dthvdz_d
    implicit none

    integer :: i, j, k
    real    :: qs
    real    :: dq,dth,dthv,temp
    real    :: a_dry, b_dry, a_moist, b_moist, c_liquid, epsilon, eps_I, chi_sat, chi
    real    :: del_thv_sat, del_thv_dry
    real    :: rvord, rvordm1, rlv2, cprv, eps1_l

    ! Formed once, as the host now does too. epsilon and eps_I in particular
    ! were two divisions per cell of quantities that depend on nothing at all,
    ! and under collapse(3) there is no outer loop for the compiler to hoist
    ! them into - every thread is one cell.
    epsilon = rd/rv
    eps_I   = 1/epsilon - 1.   !cstep approx 0.608
    rvord   = rv/rd
    rvordm1 = rv/rd - 1.
    rlv2    = rlv**2
    cprv    = cp*rv
    eps1_l  = eps1

    dthvdz_d = 0.0

    if (lmoist) then

       !$acc parallel loop collapse(3) default(present) firstprivate(rvord,rvordm1)
       do  k=kb,ke+kh
          do  j=jb,je
             do  i=ib,ie
                thv0h_d(i,j,k) = (thl0h_d(i,j,k)+rlv*ql0h_d(i,j,k)/(cp*exnh_d(k))) &
                     *(1+rvordm1*qt0h_d(i,j,k)-rvord*ql0h_d(i,j,k))
             end do
          end do
       end do
       !$acc end parallel loop

       ! epsilon and eps_I are deliberately absent from the private list.
       ! They are assigned above the region and only read inside it, so they
       ! have to arrive as firstprivate; naming them private hands every
       ! thread an uninitialised copy instead. They were private while they
       ! were still being assigned in the loop body, and leaving them there
       ! after hoisting them out cost case 502 a 0.2% error by its second
       ! step - with every parity case, Debug and Release, still passing.
       !$acc parallel loop collapse(3) default(present) &
       !$acc   private(qs,dq,dth,dthv,temp,a_dry,b_dry,a_moist,b_moist) &
       !$acc   private(c_liquid,chi_sat,chi,del_thv_sat,del_thv_dry) &
       !$acc   firstprivate(epsilon,eps_I,rlv2,cprv)
       do k=kb+1,ke
          do j=jb,je
             do i=ib,ie
                !
                !         default thv jump computed unsaturated
                !
                a_dry = 1. + eps_I * qt0_d(i,j,k)
                b_dry = eps_I * thl0_d(i,j,k)

                dth = thl0_d(i,j,k+1)-thl0_d(i,j,k-1)
                dq  = qt0_d(i,j,k+1)-qt0_d(i,j,k-1)

                del_thv_dry = a_dry   * dth + b_dry * dq

                dthv = del_thv_dry

                if  (ql0_d(i,j,k)> 0) then  !include moist thermodynamics
                   temp = thl0_d(i,j,k)*exnf_d(k)+(rlv/cp)*ql0_d(i,j,k)
                   qs   = qt0_d(i,j,k) - ql0_d(i,j,k)

                   a_moist = (1.-qt0_d(i,j,k)+qs/epsilon*(1.+rlv/(rv*temp))) &
                        /(1.+rlv2*qs/(cprv*temp**2))
                   b_moist = a_moist*rlv/cp-temp
                   c_liquid = a_dry * rlv / cp - thl0_d(i,j,k) / epsilon

                   del_thv_sat = a_moist * dth + b_moist * dq

                   chi     = chi_thv_d(k)
                   chi_sat = c_liquid * ql0_d(i,j,k) / (del_thv_dry - del_thv_sat)

                   if (chi < chi_sat) then  !mixed parcel is saturated
                      dthv = del_thv_sat
                   end if
                end if

                dthvdz_d(i,j,k) = dthv/dzhpair_d(k)
             end do
          end do
       end do
       !$acc end parallel loop

    else

       ! The host assigns thv0h = thl0h(:,:,kb:ke+kh) whole, so the halo
       ! columns come across too - unlike the moist branch above, which writes
       ! the interior only.
       !$acc parallel loop collapse(3) default(present)
       do k=kb,ke+kh
          do j=jb-jh,je+jh
             do i=ib-ih,ie+ih
                thv0h_d(i,j,k) = thl0h_d(i,j,k)
             end do
          end do
       end do
       !$acc end parallel loop

       !$acc parallel loop collapse(3) default(present)
       do k=kb+1,ke
          do j=jb,je
             do i=ib,ie
                dthvdz_d(i,j,k) = (thl0_d(i,j,k+1)-thl0_d(i,j,k-1))/dzhpair_d(k)
             end do
          end do
       end do
       !$acc end parallel loop

    end if

    ! dthvdz(:,:,kb) is left at the zero above, which is what the host writes
    ! there in both branches.

    !$acc parallel loop collapse(3) default(present) firstprivate(eps1_l)
    do k=kb,ke
       do j=jb,je
          do i=ib,ie
             if(abs(dthvdz_d(i,j,k)) < eps1_l) then
                dthvdz_d(i,j,k) = sign(eps1_l, dthvdz_d(i,j,k))
             end if
          end do
       end do
    end do
    !$acc end parallel loop

  end subroutine calthv_device

  !> Fluid-cell slab mean of thv0, formed without materialising thv0.
  !!
  !! thv0 is a full 3D field on the host and its only reader is the slab
  !! average on the next line, so the device forms the column sums straight
  !! out of thl0, ql0 and qt0. One 3D array not allocated at all - 137 MB at
  !! 256^3 - for the same expression evaluated in the same order.
  !!
  !! The all-solid lowest level is handled the way avexy_ibm_device handles
  !! it, by forming the unmasked sum for that one level and handing both to
  !! the shared avexy_ibm_finish, so the two cannot disagree about it.
  subroutine avexy_thv0_device(aver)
    use modglobal, only : ib,ie,jb,je,kb,ke,kh,rlv,cp,rd,rv
    use modfields, only : IIcs
    use modmpi,    only : avexy_ibm_finish
    use modcuda,   only : thl0_d,ql0_d,qt0_d,exnf_d,IIc_d,col_d,col_stage
    implicit none

    real, intent(out) :: aver(kb:ke+kh)

    real    :: s, averl_kb_nomask, rvord, rvordm1
    integer :: i, j, k

    rvord   = rv/rd
    rvordm1 = rv/rd - 1.

    !$acc parallel loop gang default(present) private(s) firstprivate(rvord,rvordm1)
    do k = kb, ke+kh
      s = 0.
      !$acc loop vector collapse(2) reduction(+:s)
      do j = jb, je
        do i = ib, ie
          s = s + ((thl0_d(i,j,k)+rlv*ql0_d(i,j,k)/(cp*exnf_d(k))) &
                  *(1+rvordm1*qt0_d(i,j,k)-rvord*ql0_d(i,j,k)))*IIc_d(i,j,k)
        end do
      end do
      col_d(k) = s
    end do
    !$acc end parallel loop

    averl_kb_nomask = 0.
    if (IIcs(kb) == 0) then
      s = 0.
      !$acc parallel loop collapse(2) default(present) reduction(+:s) &
      !$acc   firstprivate(rvord,rvordm1)
      do j = jb, je
        do i = ib, ie
          s = s + (thl0_d(i,j,kb)+rlv*ql0_d(i,j,kb)/(cp*exnf_d(kb))) &
                 *(1+rvordm1*qt0_d(i,j,kb)-rvord*ql0_d(i,j,kb))
        end do
      end do
      !$acc end parallel loop
      averl_kb_nomask = s
    end if

    col_stage = col_d
    call avexy_ibm_finish(aver, col_stage, averl_kb_nomask, kb, ke, kh, IIcs(kb:ke+kh), .false.)

  end subroutine avexy_thv0_device
#endif

end module modthermodynamics
