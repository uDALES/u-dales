!> \file modradiation.f90
!!  Calculates the radiative sources

!>
!!  Calculates the radiative sources
!>
!!  \author Stephan de Roode,TU Delft
!!  \author Thijs Heus, MPI-M
!!  \todo Full Radiation : Namelists and init; get rid of datafiles
!!  \todo Full Radiation : Test and debug
!!  \todo Full Radiation : Clean up the code
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
use modraddata
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initradiation
    use modglobal, only : kmax,i1,ih,j1,jh,k1,nsv,dtmax,ih,jh,btime
    use modmpi,    only : myid
    implicit none

    allocate(thlprad(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(swd(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(swu(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(lwd(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(lwu(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(albedo(2-ih:i1+ih,2-jh:j1+jh))
    thlprad = 0.
    swd = 0.
    swu = 0.
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
            call radfull
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
    deallocate(thlprad,swd,swu,lwd,lwu,albedo)
  end subroutine exitradiation
!
!
  subroutine radfull
  use radiation, only : d4stream
  use modglobal, only : imax,i1,ih,jmax,j1,jh,kmax,k1,cp,dzf,rlv
  use modfields, only : rhof, exnf,exnh, thl0,qt0,ql0,sv0
  use modsurface, only : thls,qts
  use modmicrodata, only : imicro, imicro_bulk, Nc_0,iqr
    implicit none
  real :: thlpld,thlplu,thlpsd,thlpsu
  real, dimension(k1)  :: rhof_b, exnf_b
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: temp_b, qv_b, ql_b,rr_b
  integer :: i,j,k
!take care of UCLALES z-shift for thermo variables.
      do k=1,kmax
        rhof_b(k+1)     = rhof(k)
        exnf_b(k+1)     = exnf(k)
        qv_b(:,:,k+1)   = qt0(:,:,k) - ql0(:,:,k)
        ql_b(:,:,k+1)   = ql0(:,:,k)
        temp_b(:,:,k+1) = thl0(:,:,k)*exnf(k)+(rlv/cp)*ql0(:,:,k)
        if (imicro==imicro_bulk) rr_b(:,:,k+1) = sv0(:,:,k,iqr)

      end do
      rhof_b(1)     = rhof(1)
      exnf_b(1)     = exnh(1) + 0.5*dzf(1)*(exnh(1)-exnf(1))
      ql_b(:,:,1)   = ql0(:,:,1)
      qv_b(:,:,1)   = qts + 0.5*dzf(1)*(qts-qt0(:,:,1))-ql_b(:,:,1)
      temp_b(:,:,1) = (thls + 0.5*dzf(1)*(thls-thl0(:,:,1)))*exnf_b(1)+ &
             (rlv/cp)*ql_b(:,:,1)
      if (imicro==imicro_bulk) rr_b(:,:,1) = rr_b(:,:,2)
!run radiation
             call d4stream(i1,ih,j1,jh,k1,   &  !done
             thls,          &       !done
             sfc_albedo,          &   !done
             Nc_0,           &     !done
                  rhof_b,      &     !done
                  exnf_b*cp,      &    !done
                  temp_b,   &  !done
                  qv_b,    &   !done
                  ql_b,   &    !done
                  swd,     &       !done     -OUT
                  swu,   &       !done -OUT
                  lwd,     &       !done     -OUT
                  lwu,   &       !done -OUT
                  albedo,   &      !done       -OUT
                  rr=rr_b)

!Add up thl tendency
        do k=1,kmax
        do j=2,j1
        do i=2,i1
          thlpld          = -(lwd(i,j,k+1)-lwd(i,j,k))/(rho_air_mn*cp*dzf(k))
          thlplu          = -(lwu(i,j,k+1)-lwu(i,j,k))/(rho_air_mn*cp*dzf(k))
          thlpsd          = -(swd(i,j,k+1)-swd(i,j,k))/(rho_air_mn*cp*dzf(k))
          thlpsu          = -(swu(i,j,k+1)-swu(i,j,k))/(rho_air_mn*cp*dzf(k))

          thlprad(i,j,k)  = thlprad(i,j,k) + thlpld+thlplu+thlpsu+thlpsd
        end do
        end do
        end do


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




  use modglobal,    only : i1,j1,kmax, k1,ih,jh,dzf,cp,rslabs,xtime,timee,xday,xlat,xlon
  use modfields,    only : ql0, sv0
  implicit none
  real, allocatable :: lwpt(:),lwpb(:)
  real, allocatable :: tau(:)
  real, allocatable :: absorber(:,:,:)  !  liquid water content or smoke

  real thlpld,thlplu,thlpsw
  real tauc
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


    swd = 0.0
    mu=zenith(xtime + timee/3600,xday,xlat,xlon)
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
        thlpsw          = (swd(i,j,k+1)-swd(i,j,k))/(rho_air_mn*cp*dzf(k))
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
      swd(i,j,k)=sw0*(4./3.)*(rp*(c1*exp(-rk*taupath) &
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
