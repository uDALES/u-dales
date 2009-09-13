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
  logical :: lqlnr    = .true. !< switch for ql calc. with Newton-Raphson (on/off)
  real, allocatable :: th0av(:)
  real :: chi_half=0.5  !< set wet, dry or intermediate (default) mixing over the cloud edge


contains

!> Allocate and initialize arrays
  subroutine initthermodynamics
    use modglobal, only : k1

    implicit none

    allocate(th0av(k1))
    th0av = 0.

  end subroutine initthermodynamics

!> Do moist thermodynamics.
!! Calculate the liquid water content, do the microphysics, calculate the mean hydrostatic pressure, calculate the fields at the half levels, and finally calculate the virtual potential temperature.
  subroutine thermodynamics
    use modglobal, only : lmoist, timee
    use modfields, only : thl0,thl0h,qt0,qt0h,ql0,ql0h,presf,presh,exnf,exnh
    use modmicrophysics, only : microphysics
    use modmicrophysics, only : imicro, imicro_none, imicro_drizzle
    implicit none
    if (timee==0) call diagfld
    if (lmoist) then
!       if (imicro==imicro_none .or. imicro == imicro_drizzle)
      call thermo(thl0,qt0,ql0,presf,exnf)
      if (imicro/=imicro_none) call microphysics
    end if

    call diagfld
    call calc_halflev !calculate halflevel values of qt0 and thl0
    if (lmoist) then
!       if (imicro==imicro_none .or. imicro == imicro_drizzle)
      call thermo(thl0h,qt0h,ql0h,presh,exnh)
    end if
    call calthv

  end subroutine thermodynamics
!> Cleans up after the run
  subroutine exitthermodynamics
  implicit none
    deallocate(th0av)
  end subroutine exitthermodynamics

!> Calculate thetav and dthvdz
  subroutine calthv
    use modglobal, only : lmoist,i1,j1,k1,kmax,zf,zh,dzh,rlv,rd,rv,cp,eps1
    use modfields, only : thl0,thl0h,ql0,ql0h,qt0,qt0h,sv0,exnf,exnh,thv0h,dthvdz
    use modsurface,only : dthldz,dqtdz

    implicit none

    integer i, j, k
    real    qs
    real    c1,c2,dq,dth,dthv,temp
    real    a_dry, b_dry, a_moist, b_moist, c_liquid, epsilon, eps_I, chi_sat, chi
    real    del_thv_sat, del_thv_dry
    dthvdz = 0
    if (lmoist) then

      do  j=2,j1
      do  i=2,i1
      do  k=2,k1
        thv0h(i,j,k) = (thl0h(i,j,k)+rlv*ql0h(i,j,k)/(cp*exnh(k))) &
                      *(1+(rv/rd-1)*qt0h(i,j,k)-rv/rd*ql0h(i,j,k))
      end do
      end do
      end do

      do k=2,kmax
      do j=2,j1
      do i=2,i1
!
!     default thv jump computed unsaturated
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
      do j=2,j1
      do i=2,i1

        if(ql0(i,j,1)>0) then
          temp = thl0(i,j,1)*exnf(1)+(rlv/cp)*ql0(i,j,1)
          qs   = qt0(i,j,1) - ql0(i,j,1)
          c1   = (1.-qt0(i,j,1)+rv/rd*qs*(1.+rlv**2/(rv*temp))) &
                    /(1.+rlv**2*qs/(cp*rv*temp**2))
          c2   = c1*rlv/(temp*cp)-1.

        else
          c1 = 1.+(rv/rd-1)*qt0(i,j,1)
          c2 = rv/rd-1

        end if
        dthvdz(i,j,1) = c1*dthldz(i,j) + c2*thl0(i,j,1)*dqtdz(i,j)

      end do
      end do

    else
      thv0h = thl0h
      do k=2,kmax
      do j=2,j1
      do i=2,i1
        dthvdz(i,j,k) = (thl0(i,j,k+1)-thl0(i,j,k-1))/(dzh(k+1)+dzh(k))
      end do
      end do
      end do
      do  j=2,j1
      do  i=2,i1
        dthvdz(i,j,1) = dthldz(i,j)
      end do
      end do
    end if
    where (abs(dthvdz)<eps1)
      dthvdz = sign(eps1,dthvdz)
    end where

  end subroutine calthv
!> Calculate diagnostic slab averaged fields.
!!     Calculates slab averaged fields assuming
!!     hydrostatic equilibrium for: u,v,theta_l,theta_v,
!!     qt,ql,exner,pressure and the density
!! \author      Pier Siebesma   K.N.M.I.     06/01/1995
  subroutine diagfld

  use modglobal, only : i1,ih,j1,jh,k1,nsv,zh,zf,cu,cv,rslabs,grav,rlv,cp,rd,rv,pref0
  use modfields, only : u0,v0,w0,thl0,qt0,ql0,sv0,u0av,v0av,thl0av,qt0av,ql0av,sv0av, &
                        presf,presh,exnf,exnh, rhof
  use modsurface,only : thls,ps
  use modmpi,    only : slabsum
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
  call slabsum(u0av  ,1,k1,u0  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
 call slabsum(v0av  ,1,k1,v0  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
  call slabsum(thl0av,1,k1,thl0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
  call slabsum(qt0av ,1,k1,qt0 ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
  call slabsum(ql0av ,1,k1,ql0 ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)

    u0av  = u0av  /rslabs + cu
    v0av  = v0av  /rslabs + cv
    thl0av = thl0av/rslabs
    qt0av = qt0av /rslabs
    ql0av = ql0av /rslabs
   exnf   = 1-grav*zf/(cp*thls)
    exnh  = 1-grav*zh/(cp*thls)
   th0av  = thl0av + (rlv/cp)*ql0av/exnf

  do n=1,nsv
    call slabsum(sv0av(1,n),1,k1,sv0(1,1,1,n),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
  end do
  sv0av = sv0av/rslabs

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

  exnh(1) = (ps/pref0)**(rd/cp)
  exnf(1) = (presf(1)/pref0)**(rd/cp)
  do k=2,k1
    exnf(k) = (presf(k)/pref0)**(rd/cp)
    exnh(k) = (presh(k)/pref0)**(rd/cp)
  end do

  tv      = th0av(1)*exnf(1)*(1+(rv/rd-1)*qt0av(1)-rv/rd*ql0av(1))
  rhof(1) = presf(1)/(rd*tv)

!    3.2 determine rho

  do k=2,k1
    tv      = th0av(k)*exnf(k)*(1.+(rv/rd-1)*qt0av(k)-rv/rd*ql0av(k))
    rhof(k) = presf(k)/(rd*tv)
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

  use modglobal, only : k1,dzf,dzh,rv,rd,cp,tmelt,zf,grav,pref0
  use modfields, only : qt0av,ql0av,presf,presh
  use modsurface,only : ps,thls
  implicit none

  integer   k
  real      thetav, rdocp
  real,allocatable,dimension (:) :: thetah, qth, qlh

  allocate(thetah(k1), qth(k1), qlh(k1))
  rdocp = rd/cp

!**************************************************
!    1.0 Determine theta and qt at half levels    *
!**************************************************

  do k=2,k1
    thetah(k) = (th0av(k)*dzf(k-1) + th0av(k-1)*dzf(k))/(2*dzh(k))
    qth   (k) = (qt0av(k)*dzf(k-1) + qt0av(k-1)*dzf(k))/(2*dzh(k))
    qlh   (k) = (ql0av(k)*dzf(k-1) + ql0av(k-1)*dzf(k))/(2*dzh(k))
  end do

!**************************************************
!     2.1  calculate pressures at full levels     *
!          assuming hydrostatic equilibrium       *
!**************************************************

!     1: lowest level

  thetav   = th0av(1)*(1+(rv/rd-1)*qt0av(1)-rv/rd*ql0av(1))
  presf(1) = ps**rdocp - &
                 grav*(pref0**rdocp)*zf(1) /(cp*thetav)
  presf(1) = presf(1)**(1./rdocp)

!     2: higher levels

  do k=2,k1

    thetav   = thetah(k)*(1+(rv/rd-1)*qth(k)-rv/rd*qlh(k))
    presf(k) = presf(k-1)**rdocp - &
                   grav*(pref0**rdocp)*dzh(k) /(cp*thetav)
    presf(k) = presf(k)**(1./rdocp)
  end do

!**************************************************
!     2.2   calculate pressures at half levels    *
!           assuming hydrostatic equilibrium      *
!**************************************************

  presh(1) = ps
  do k=2,k1
    thetav   = th0av(k-1)*(1+(rv/rd-1)*qt0av(k-1)-rv/rd*ql0av(k-1))
    presh(k) = presh(k-1)**rdocp - &
                   grav*(pref0**rdocp)*dzf(k-1) / (cp*thetav)
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



  use modglobal, only : ih,jh,i1,j1,k1,es0,at,bt,rd,rv,rlv,cp,tmelt
  use modsurface, only : thls
  implicit none

  integer i, j, k
  real tl, es, qs, qsl, b1
  real, intent(in)  :: qt(2-ih:i1+ih,2-jh:j1+jh,k1),thl(2-ih:i1+ih,2-jh:j1+jh,k1),exner(k1),pressure(k1)
  real, intent(out) :: ql(2-ih:i1+ih,2-jh:j1+jh,k1)
  real :: Tnr,qsatur,Tnr_old
  integer :: niter,nitert
    if (lqlnr) then

!mc      calculation of T with Newton-Raphson method
!mc      first guess is Tnr=tl
!mc
      nitert = 0
      do j=2,j1
        do i=2,i1
          do k=1,k1

            tl  = thl(i,j,k)*exner(k)
            Tnr=tl
            Tnr_old=0
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
            niter = 0.0

            ql(i,j,k) = dim(qt(i,j,k)-qsatur,0.)

          end do
        end do
      end do
    else


      do j=2,j1
        do i=2,i1
          do k=1,k1
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
    use modglobal, only : i1,j1,k1,dzf,dzh,iadv_thl, iadv_qt, iadv_kappa
    use modfields, only : thl0,thl0h,qt0,qt0h
    use modsurface,only: qts,thls
    implicit none

    integer :: i,j,k


    if (iadv_thl==iadv_kappa) then
      call halflev_kappa(thl0,thl0h)
    else
      do  j=2,j1
      do  i=2,i1
      do  k=2,k1
        thl0h(i,j,k) = (thl0(i,j,k)*dzf(k-1)+thl0(i,j,k-1)*dzf(k))/(2*dzh(k))
      end do
      end do
      end do
    end if
    thl0h(2:i1,2:j1,1) = thls

    if (iadv_qt==iadv_kappa) then
        call halflev_kappa(qt0,qt0h)
    else
      do  j=2,j1
      do  i=2,i1
      do  k=2,k1
        qt0h(i,j,k)  = (qt0 (i,j,k)*dzf(k-1)+qt0 (i,j,k-1)*dzf(k))/(2*dzh(k))
      end do
      end do
      end do
      qt0h(2:i1,2:j1,1)  = qts
    end if
  end subroutine calc_halflev

end module modthermodynamics
