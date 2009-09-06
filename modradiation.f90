!> \file modradiation.f90
!!  Calculates the radiative sources

!>
!!  Calculates the radiative sources
!>
!!  \author Stephan de Roode,TU Delft
!!  \todo Full Radiation
!!  \todo Documentation
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



module modradiation

implicit none

SAVE

  integer, parameter :: irad_none  = 0   !< 0=no radiation
  integer, parameter :: irad_full  = 1   !< 1=full radiation
  integer, parameter :: irad_par   = 2   !< 2=parameterized radiation
  integer, parameter :: irad_user  = 10  !< 10=user specified radiation

  logical :: rad_ls      = .true.   !< prescribed radiative forcing
  logical :: rad_longw   = .true.   !< parameterized longwave radiative forcing
  logical :: rad_shortw  = .true.   !< parameterized shortwave radiative forcing
  logical :: rad_smoke   = .false.  !< longwave divergence for smoke cloud

  real              :: timerad = 0 !<  timescale of the radiation scheme
  real              :: tnext   = 0
  real :: rho_air_mn = 1.1436 !< mean air density used in radiation computation
  real :: rka        = 130.   !< extinction coefficient in radpar scheme
  real :: dlwtop     = 74.    !< longwave radiative flux divergence at top of domain
  real :: dlwbot     = 0.     !< longwave radiative flux divergence near the surface
  real :: sw0        = 1100.0 !< direct component at top of the cloud (W/m^2), diffuse not possible
  real :: gc         = 0.85   !< asymmetry factor of droplet scattering angle distribution
  real :: sfc_albedo = 0.05   !< ground surface albedo
  real :: reff       = 1.e-5  !< cloud droplet effective radius (m)
  integer :: isvsmoke = 1     !< number of passive scalar to be used for optical depth calculation
  integer :: iradiation = irad_par
  integer :: irad    = -1


  real mu                    !< cosine of the solar zenith angle

  real, allocatable :: thlprad(:,:,:)                      !<   the radiative tendencies
  real, allocatable :: swn(:,:,:),lwd(:,:,:),lwu(:,:,:)    !<   shortwave radiative flux


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initradiation
    use modglobal, only : kmax,i1,j1,k1,nsv,dtmax,ih,jh,btime
    use modmpi,    only : myid
    implicit none

    allocate(thlprad(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(swn(2:i1,2:j1,k1))
    allocate(lwd(2:i1,2:j1,k1))
    allocate(lwu(2:i1,2:j1,k1))
    thlprad = 0.
    swn = 0.
    lwd = 0.
    lwu = 0.
    if (irad/=-1) then
      if (myid==0) write (*,*) 'WARNING: The use of irad is deprecated. Please use the iradiation switch'
      select case (irad)
      case (0)
        iradiation = 0
      case (1)
        iradiation = 2
        rad_ls     = .true.
        rad_longw  = .false.
        rad_shortw = .false.
        rad_smoke  = .false.
      case (2)
        iradiation = 2
        rad_ls     = .false.
        rad_longw  = .true.
        rad_shortw = .false.
        rad_smoke  = .false.
      case (3)
        iradiation = 1
      case (4)
        iradiation = 2
        rad_ls     = .false.
        rad_longw  = .true.
        rad_shortw = .true.
        rad_smoke  = .false.
      case (10)
        iradiation = 2
        rad_ls     = .false.
        rad_longw  = .false.
        rad_shortw = .false.
        rad_smoke  = .true.
      end select
    end if
    tnext = -1e-3+btime

    if (rad_smoke.and.isvsmoke>nsv) then
      if (rad_shortw) then
         stop 'you want to compute solar radiative transfer through a smoke cloud'
      endif
      stop 'Smoke radiation with wrong (non-existent?) scalar field'
    endif

    if (iradiation==1) then
       stop 'full scheme not yet implemented'
    endif

  end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine radiation
    use modglobal, only : timee, dt_lim,rk3step
    use modfields, only : thlp
    implicit none

    if(timee<tnext .and. rk3step==3) then
      dt_lim = min(dt_lim,tnext-timee)
    end if

    if(timerad==0 .or. (timee>tnext .and. rk3step==1)) then
      tnext = tnext+timerad
      thlprad = 0.0
      select case (iradiation)
          case (irad_none)
          case (irad_full)
            stop 'full scheme not yet implemented'
            call rad_full
          case (irad_par)
            if (rad_ls) then
                call radprof
            endif

            if(rad_longw.or.rad_shortw) then
              call radpar
            endif
          case (irad_user)
            call rad_user
            if (rad_ls) then
                call radprof
            endif

            if(rad_longw.or.rad_shortw) then
              call radpar
            endif
      end select
    end if
    thlp = thlp + thlprad


  end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine exitradiation
    implicit none
    deallocate(thlprad,swn,lwd,lwu)
  end subroutine exitradiation
!
!
  subroutine radfull
  use radiation, only : d4stream
  use modglobal, only : i1,j1,k1
  use modsurface, only : thls
    implicit none
!       integer,parameter :: nzp=k1, nxp=1, nyp=1
      real    :: cntlat, time_in, sst,  CCN
      real, dimension (k1)                 :: dn0, pi0, pi1, dzt
      real, dimension (k1,i1,j1)          :: a_pexnr, a_scr1, a_rv, a_rc
      real, dimension (k1,i1,j1)        :: a_tt, a_rflx, a_sflx, albedo, a_rpp
!convert our variables into bjorn variables, including switching column order

!run radiation
!         if (present(time_in) .and. present(cntlat) .and. present(thls)) then
!           if (level == 3) then
!                       d4stream(n1, n2, n3, alat, time, sknt, sfc_albedo, CCN, dn0, &
!          pi0, pi1, dzm, pip, tk, rv, rc, tt, rflx, sflx, albedo, rr)
                 call d4stream(k1,i1,j1, cntlat, time_in, sst, 0.05, CCN,&
                  dn0, pi0, pi1, dzt, a_pexnr, a_scr1, a_rv, a_rc, a_tt,  &
                  a_rflx, a_sflx, albedo, rr=a_rpp)
!           else
!              call d4stream(nzp, nxp, nyp, cntlat, time_in, sst, 0.05, CCN,&
!                   dn0, pi0, pi1, dzt, a_pexnr, a_scr1, a_rv, a_rc, a_tt,  &
!                   a_rflx, a_sflx, albedo)
!           end if
!         end if
!convert back

end subroutine radfull


 subroutine radpar
!-----------------------------------------------------------------|
!                                                                 |
!*** *radpar*  calculates tendency due to parameterized radiation |
!                                                                 |
!      Hans Cuijpers   I.M.A.U.  06/01/1995                       |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!      In this subroutine the radiative flux divergence of the    |
!      cloud is calculated                                        |

!      Margreet van Zanten IMAU juli 2000
!      SW flux calculation added
!      short version of radiat.f
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!-----------------------------------------------------------------|




  use modglobal,    only : i1,j1,kmax, k1,ih,jh,dzf,timee,pi,cp,xtime,xday,xlon,xlat,rslabs
  use modfields,    only : ql0, sv0
  implicit none
  real, allocatable :: lwpt(:),lwpb(:)
  real, allocatable :: tau(:)
  real, allocatable :: absorber(:,:,:)  !  liquid water content or smoke

  real thlpld,thlplu,thlpsw
  real obliq,xlam,phi,el,declin,hora,rtime,tauc
  integer :: i=0, j=0, k

  real :: rho_l= 1000.  !liquid water density (kg/m3)

  allocate(lwpt(k1),lwpb(k1))
  allocate(tau(k1))
  allocate(absorber(2-ih:i1+ih,2-jh:j1+jh,k1))
  absorber = 0.0
  tau  = 0.0
  lwpt = 0.0
  lwpb = 0.0


! MPI
! initialise local variables

   if (rad_longw) then
     if (rad_smoke) then
        absorber(2-ih:i1+ih,2-jh:j1+jh,1:k1) = sv0(2-ih:i1+ih,2-jh:j1+jh,1:k1,isvsmoke)
     else
        absorber(2-ih:i1+ih,2-jh:j1+jh,1:k1) = ql0(2-ih:i1+ih,2-jh:j1+jh,1:k1)
     endif

     do j=2,j1
     do i=2,i1
       lwpt = 0.
       lwpb = 0.

! **   Downward LWP

       do k=kmax,1,-1
        lwpt(k) = lwpt(k+1) + rho_air_mn*absorber(i,j,k)*dzf(k)
       end do

! **   Upward LWP

       do k=2,k1
         lwpb(k) = lwpb(k-1) + rho_air_mn*absorber(i,j,k-1)*dzf(k-1)
       end do

       do k=1,k1
        lwd(i,j,k) = dlwtop*exp(-rka*lwpt(k))
        lwu(i,j,k) = dlwbot*exp(-rka*lwpb(k))
      end do

       do k=1,kmax
         thlpld         = -(lwd(i,j,k+1)-lwd(i,j,k))/(rho_air_mn*cp*dzf(k))
         thlplu         = -(lwu(i,j,k+1)-lwu(i,j,k))/(rho_air_mn*cp*dzf(k))
         thlprad(i,j,k) =   thlprad(i,j,k) + thlpld+thlplu
       end do

    end do
    end do  ! end i,j loop

  endif  !end longwave loop

!----------------------------------------------------------------------

  if (rad_shortw) then

    !compute solar zenith angle

    rtime   = xtime + timee/3600
    phi    = xlat * pi/180.
    el     = xlon * pi/180.
    obliq  = 23.45 * pi/180.
    xlam   = 4.88 + 0.0172 * xday
    declin = asin(sin(obliq)*sin(xlam))
    hora   = el-pi + 2.*pi*(rtime/24.)
    mu   = max(0.,sin(declin)*sin(phi)+cos(declin)*cos(phi)* &
                                                         cos(hora))
    swn = 0.0

    do j=2,j1
    do i=2,i1

      if (mu > 0.035) then  !factor 0.035 needed for security
        tauc = 0.           ! tau cloud
        do k = 1,kmax
          tau(k) = 0.      ! tau laagje dz
          if (ql0(i,j,k) > 1e-5) then
            tau(k)=1.5*ql0(i,j,k)*rho_air_mn*dzf(k)/reff/rho_l
            tauc=tauc+tau(k)
          end if
        end do
        call sunray(tau,tauc,i,j)
      end if

      do k=1,kmax
        thlpsw          = (swn(i,j,k+1)-swn(i,j,k))/(rho_air_mn*cp*dzf(k))
        thlprad(i,j,k)  = thlprad(i,j,k) + thlpsw

      end do

    end do
    end do

  endif  !end shortwave loop

  deallocate(lwpt,lwpb,tau,absorber)



  return
  end subroutine radpar



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
  subroutine sunray(tau,tauc,i,j)

!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!   the sunray model is described by fouquart and  bonnel                  c
!   (1980, contr. atmos. phys.).                                                               c

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  use modglobal, only :  k1
  implicit none

  real, intent(inout), dimension (k1) :: tau
  integer, intent(in) :: i,j
  real,allocatable, dimension (:) :: taude
  real gcde,tauc,taucde &
           ,taupath,t1,t2,t3,c1,c2 &
           ,omega,omegade,ff,x1,x2,x3,rk,mu2,rp,alpha,beta,rtt &
           ,exmu0,expk,exmk,xp23p,xm23p,ap23b
  integer k
  allocate(taude(k1))

  taucde = 0.         ! tau' cloud
  taupath = 0.
  do k=1,k1
     taude(k) =0.     ! tau' laagje dz
  end do

! omega=0.9989-4.e-3*exp(-0.15*tauc)  !   fouquart and bonnel (1980)
  omega=1.-1.e-3*(0.9+2.75*(mu+1.)*exp(-0.09*tauc)) !fouquart

! the equations for the delta-eddington approximation are equal to those
! for the eddington approximation with transformed parameters g, omega
! and tau (joseph, wiscomb and weinman, 1976, j.a.s.).
! parameternames: x -> xde (delta-eddington)

  ff=gc*gc
  gcde=gc/(1.+gc)
  taucde=(1.0-omega*ff)*tauc

  do k =1,k1
    taude(k)=(1.e0-omega*ff)*tau(k)
  end do

  omegade=(1.0-ff)*omega/(1.e0-omega*ff)

! the solution of the eddington equations are given by shettle and weinman
! (1970, j.a.s.).

  x1=1.0-omegade*gcde
  x2=1.0-omegade
  rk=sqrt(3.0*x2*x1)
  mu2=mu*mu
  x3=4.0*(1.0-rk*rk*mu2)
  rp=sqrt(3.0*x2/x1)
  alpha=3.e0*omegade*mu2*(1.0+gcde*x2)/x3
  beta=3.*omegade*mu*(1.0+3.0*gcde*mu2*x2)/x3

  rtt=2.0/3.0
  exmu0= exp(-taucde/mu)
  expk=  exp(rk*taucde)
  exmk=1.0/expk
  xp23p=1.0+rtt*rp
  xm23p=1.0-rtt*rp
  ap23b=alpha+rtt*beta

  t1=1-sfc_albedo-rtt*(1.+sfc_albedo)*rp
  t2=1-sfc_albedo+rtt*(1.+sfc_albedo)*rp
  t3=(1-sfc_albedo)*alpha-rtt*(1+sfc_albedo)*beta+sfc_albedo*mu
  c2=(xp23p*t3*exmu0-t1*ap23b*exmk)/(xp23p*t2*expk-xm23p*t1*exmk)
  c1=(ap23b-c2*xm23p)/xp23p

  do k = k1,1,-1
      taupath = taupath + taude(k)
      swn(i,j,k)=sw0*(4./3.)*(rp*(c1*exp(-rk*taupath) &
                 -c2*exp(rk*taupath)) &
                 -beta*exp(-taupath/mu)) &
                 +mu*sw0*exp(-taupath/mu)
  end do

  deallocate(taude)

  return
  end subroutine sunray


!***********************************************************************
!***  In this subroutine the a precribed radiative tendency     ********
!***  is taken into account.                                    ********
!***********************************************************************

  subroutine radprof
  use modglobal,    only : i1,j1,kmax
  use modfields,    only : thlpcar
  implicit none
  integer k

  do k=1,kmax
    thlprad(2:i1,2:j1,k) = thlprad(2:i1,2:j1,k) + thlpcar(k)
  end do

  return
  end subroutine radprof


end module
