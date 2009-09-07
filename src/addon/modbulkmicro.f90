!> \file modbulkmicro.f90

!>
!!  Bulk microphysics.
!>
!! Calculates bulk microphysics using a two moment scheme.
!! \see   Seifert and Beheng (Atm. Res., 2001)
!! \see  Seifert and Beheng (Met Atm Phys, 2006)
!! \see  Stevens and Seifert (J. Meteorol. Soc. Japan, 2008)  (rain sedim, mur param)
!! \see  Seifert (J. Atm Sc., 2008) (rain evap)
!! \see  Khairoutdinov and Kogan (2000) (drizzle param : auto, accr, sedim, evap)
!!  \author Olivier Geoffroy, K.N.M.I.
!!  \author Margreet van Zanten, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \par Revision list
!! \todo documentation
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


module modbulkmicro

!
!
!
!   Amount of liquid water is splitted into cloud water and precipitable
!   water (as such it is a two moment scheme). Cloud droplet number conc. is
!   fixed in place and time.
!
!   rhoz is used to convert mixing ratio (kg/kg) in water content (kg/m^3)
!   rhoz chosen to be rhof(k) (in modthermodynamics,modstartup), is rhof(1) better?
!   same rhoz value used for some diagnostics calculation (in modbulkmicrostat, modtimestat)
!
!   Cond. sampled timeav averaged profiles are weighted with fraction of condition,
!   similarly as is done in sampling.f90
!
!   bulkmicro is called from *modmicrophysics*
!
!*********************************************************************
  use modbulkmicrodata

  implicit none
  save
  contains

!> Initializes and allocates the arrays
  subroutine initbulkmicro
    use modglobal, only : ih,i1,jh,j1,k1,dzf,dtmax,rk3step
    use modmpi,    only : myid
    implicit none

    allocate (Nr(2-ih:i1+ih,2-jh:j1+jh,k1)        &
             ,Nrp(2-ih:i1+ih,2-jh:j1+jh,k1)       &
             ,qltot(2-ih:i1+ih,2-jh:j1+jh,k1)     &
             ,qr(2-ih:i1+ih,2-jh:j1+jh,k1)        &
             ,qrp(2-ih:i1+ih,2-jh:j1+jh,k1)       &
             ,qc(2-ih:i1+ih,2-jh:j1+jh,k1)        &
             ,Nc(2-ih:i1+ih,2-jh:j1+jh,k1)        &
             ,nuc(2-ih:i1+ih,2-jh:j1+jh,k1)        &
             ,thlpmcr(2-ih:i1+ih,2-jh:j1+jh,k1)   &
             ,qtpmcr(2-ih:i1+ih,2-jh:j1+jh,k1)    &
             ,sedc(2-ih:i1+ih,2-jh:j1+jh,k1)      &
             ,sed_qr(2-ih:i1+ih,2-jh:j1+jh,k1)    &
             ,sed_Nr(2-ih:i1+ih,2-jh:j1+jh,k1)    &
             ,exnz(2-ih:i1+ih,2-jh:j1+jh,k1)      &
             ,presz(2-ih:i1+ih,2-jh:j1+jh,k1)     &
             ,Dvc(2-ih:i1+ih,2-jh:j1+jh,k1)       &
             ,xc(2-ih:i1+ih,2-jh:j1+jh,k1)        &
             ,Dvr(2-ih:i1+ih,2-jh:j1+jh,k1)       &
             ,xr(2-ih:i1+ih,2-jh:j1+jh,k1)        &
             ,mur(2-ih:i1+ih,2-jh:j1+jh,k1)       &
             ,lbdr(2-ih:i1+ih,2-jh:j1+jh,k1)      &
             ,au(2-ih:i1+ih,2-jh:j1+jh,k1)        &
             ,phi(2-ih:i1+ih,2-jh:j1+jh,k1)       &
             ,tau(2-ih:i1+ih,2-jh:j1+jh,k1)       &
             ,rhoz(2-ih:i1+ih,2-jh:j1+jh,k1)      &
             ,ac(2-ih:i1+ih,2-jh:j1+jh,k1)        &
             ,sc(2-ih:i1+ih,2-jh:j1+jh,k1)        &
             ,br(2-ih:i1+ih,2-jh:j1+jh,k1)        &
             ,evap(2-ih:i1+ih,2-jh:j1+jh,k1)      &
             ,Nevap(2-ih:i1+ih,2-jh:j1+jh,k1)     &
             ,qr_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    &
             ,Nr_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    &
             ,wfall_qr(2-ih:i1+ih,2-jh:j1+jh,k1)  &
             ,wfall_Nr(2-ih:i1+ih,2-jh:j1+jh,k1))

    allocate(precep(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(qrmask(2-ih:i1+ih,2-jh:j1+jh,k1),qcmask(2-ih:i1+ih,2-jh:j1+jh,k1))

!


  end subroutine initbulkmicro

!> Cleaning up after the run
  subroutine exitbulkmicro
  !*********************************************************************
  ! subroutine exitbulkmicro
  !*********************************************************************
    implicit none

    deallocate(Nr,Nrp,qltot,qr,qrp,thlpmcr,qtpmcr  )

    deallocate(nuc,sedc,sed_qr,sed_Nr,exnz,presz,Dvc,xc,Dvr,mur,lbdr, &
               xr,au,phi,tau,ac,sc,br,evap,Nevap,qr_spl,Nr_spl,wfall_qr,wfall_Nr)

    deallocate(precep)

  end subroutine exitbulkmicro

!> Calculates the microphysical source term.
  subroutine bulkmicro
    use modglobal, only : dt,rk3step,timee,kmax,rlv,cp
    use modfields, only : sv0,svm,svp,qtp,thlp,qt0,ql0,presf, exnf,rhof
    use modbulkmicrostat, only : bulkmicrotend
    use modmpi,    only : myid
    implicit none
    integer :: k

    Nr     = sv0(2-ih:i1+ih,2-jh:j1+jh,1:k1,inr)
    qr    = sv0(2-ih:i1+ih,2-jh:j1+jh,1:k1,iqr)
    Nrp    = 0.0
    qrp   = 0.0
    thlpmcr = 0.0
    qtpmcr  = 0.0
    Nc    = 0.0

    delt = dt/ (4. - dble(rk3step))
    do k = 1, k1
      exnz(2:i1,2:j1,k) = exnf(k)
      presz(2:i1,2:j1,k) = presf(k)
      rhoz(2:i1,2:j1,k) = rhof(k)
    enddo

    if ( timee .eq. 0. .and. rk3step .eq. 1 .and. myid .eq. 0) then
      write(*,*) 'l_lognormal',l_lognormal
      write(*,*) 'rhoz(1)', rhoz(2,2,1),' rhoz(10)', rhoz(2,2,10)
      write(*,*) 'l_mur_cst',l_mur_cst,' mur_cst',mur_cst
      write(*,*) 'nuc = param'
    endif

  !*********************************************************************
  ! remove neg. values of Nr and qr
  !*********************************************************************
    if (l_rain) then
    if (sum(qr, qr<0.) > 0.000001*sum(qr)) then
      write(*,*)'amount of neg. qr and Nr thrown away is too high  ',timee,' sec'
    end if
    if (sum(Nr, Nr<0.) > 0.000001*sum(Nr)) then
      write(*,*)'amount of neg. qr and Nr thrown away is too high  ',timee,' sec'
    end if

      where (Nr  < 0. ) Nr  = 0.
      where (qr < 0. ) qr = 0.
    endif   ! l_rain

!    where (Nr <= 0. .or. qr<=qrmin )
    where (qr<=qrmin )
!     Nr  = 0.
!     qr  = 0.
     qrmask = .false.
    elsewhere
     qrmask = .true.
    end where

  !*********************************************************************
  ! calculate qltot and initialize cloud droplet number Nc
  !*********************************************************************
    qc = ql0
    where (qc > qcmin)
      Nc = Nc_0
      qcmask = .true.
    elsewhere
      qcmask = .false.
    end where
    qltot = qc + qr

  !*********************************************************************
  ! calculate Rain DSD integral properties & parameters xr, Dvr, lbdr, mur
  !*********************************************************************
    if (l_rain) then
      xr(2:i1,2:j1,1:k1) = 0.
      Dvr(2:i1,2:j1,1:k1) = 0.
      mur(2:i1,2:j1,1:k1) = 30.
      lbdr(2:i1,2:j1,1:k1) = 0.
      if (l_sb ) then
        where (qrmask(2:i1,2:j1,1:k1))
          xr(2:i1,2:j1,1:k1) = rhoz(2:i1,2:j1,1:k1)*qr(2:i1,2:j1,1:k1)/(Nr(2:i1,2:j1,1:k1))
          xr(2:i1,2:j1,1:k1) = min(max(xr(2:i1,2:j1,1:k1),xrmin),xrmax) ! to ensure xr is within borders
          Dvr(2:i1,2:j1,1:k1) = (xr(2:i1,2:j1,1:k1)/pirhow)**(1./3.)
        endwhere
        if (l_mur_cst) then
        ! mur = cst
          where (qrmask(2:i1,2:j1,1:k1))
            mur(2:i1,2:j1,1:k1) = mur_cst
            lbdr(2:i1,2:j1,1:k1) = ((mur(2:i1,2:j1,1:k1)+3.)*(mur(2:i1,2:j1,1:k1)+2.)*(mur(2:i1,2:j1,1:k1)+1.))**(1./3.)/Dvr(2:i1,2:j1,1:k1)
          endwhere
        else
        ! mur = f(Dv)
          where (qrmask(2:i1,2:j1,1:k1))
!            mur(2:i1,2:j1,1:k1) = 10. * (1+tanh(1200.*(Dvr(2:i1,2:j1,1:k1)-0.0014))) ! Stevens & Seifert (2008) param
            mur(2:i1,2:j1,1:k1) = min(30.,- 1. + 0.008/ (qr(2:i1,2:j1,1:k1)*rhoz(2:i1,2:j1,1:k1))**0.6)  ! G09b

            lbdr(2:i1,2:j1,1:k1) = ((mur(2:i1,2:j1,1:k1)+3.)*(mur(2:i1,2:j1,1:k1)+2.)*(mur(2:i1,2:j1,1:k1)+1.))**(1./3.)/Dvr(2:i1,2:j1,1:k1)
          endwhere
        endif

      else

        where (qrmask(2:i1,2:j1,1:k1))
          xr(2:i1,2:j1,1:k1) = rhoz(2:i1,2:j1,1:k1)*qr(2:i1,2:j1,1:k1)/(Nr(2:i1,2:j1,1:k1))
          xr(2:i1,2:j1,1:k1) = min(xr(2:i1,2:j1,1:k1),xrmaxkk) ! to ensure x_pw is within borders
          Dvr(2:i1,2:j1,1:k1) = (xr(2:i1,2:j1,1:k1)/pirhow)**(1./3.)
        endwhere
      endif ! l_sb
    endif   ! l_rain

  !*********************************************************************
  ! call microphysical processes subroutines
  !*********************************************************************
    if (l_sedc)  call sedimentation_cloud

    if (l_rain) then
      call bulkmicrotend
      call autoconversion
      call bulkmicrotend
      call accretion
      call bulkmicrotend
      call evaporation
      call bulkmicrotend
      call sedimentation_rain
      call bulkmicrotend
    endif

    sv0(2:i1,2:j1,1:k1,inr)=Nr(2:i1,2:j1,1:k1)
    sv0(2:i1,2:j1,1:k1,iqr)=qr(2:i1,2:j1,1:k1)

    svp(2:i1,2:j1,1:k1,inr)=svp(2:i1,2:j1,1:k1,inr)+Nrp(2:i1,2:j1,1:k1)
    svp(2:i1,2:j1,1:k1,iqr)=svp(2:i1,2:j1,1:k1,iqr)+qrp(2:i1,2:j1,1:k1)

    thlp(2:i1,2:j1,1:k1)=thlp(2:i1,2:j1,1:k1)+thlpmcr(2:i1,2:j1,1:k1)
    qtp(2:i1,2:j1,1:k1)=qtp(2:i1,2:j1,1:k1)+qtpmcr(2:i1,2:j1,1:k1)

  !*********************************************************************
  ! remove negative values and non physical low values
  !*********************************************************************
    where (svp(2:i1,2:j1,1:k1,iqr)+svm(2:i1,2:j1,1:k1,iqr)/delt .lt. qrmin)
      svp(2:i1,2:j1,1:k1,iqr) = - svm(2:i1,2:j1,1:k1,iqr)/delt
      svp(2:i1,2:j1,1:k1,inr)  = - svm(2:i1,2:j1,1:k1,inr)/delt
      qtp(2:i1,2:j1,1:k1) = qtp(2:i1,2:j1,1:k1) + svm(2:i1,2:j1,1:k1,iqr)/delt
      thlp(2:i1,2:j1,1:k1) = thlp(2:i1,2:j1,1:k1) - (rlv/(cp*exnz(2:i1,2:j1,1:k1)))*svm(2:i1,2:j1,1:k1,iqr)/delt
    endwhere
    where (svp(2:i1,2:j1,1:k1,inr)+svm(2:i1,2:j1,1:k1,inr)/delt .lt. 0.)
      svp(2:i1,2:j1,1:k1,iqr) = - svm(2:i1,2:j1,1:k1,iqr)/delt
      svp(2:i1,2:j1,1:k1,inr)  = - svm(2:i1,2:j1,1:k1,inr)/delt
      qtp(2:i1,2:j1,1:k1) = qtp(2:i1,2:j1,1:k1) + svm(2:i1,2:j1,1:k1,iqr)/delt
      thlp(2:i1,2:j1,1:k1) = thlp(2:i1,2:j1,1:k1) - (rlv/(cp*exnz(2:i1,2:j1,1:k1)))*svm(2:i1,2:j1,1:k1,iqr)/delt
    endwhere

  end subroutine bulkmicro
  !> Determine autoconversion rate and adjust qrp and Nrp accordingly
  !!
  !!   The autoconversion rate is formulated for f(x)=A*x**(nuc)*exp(-Bx),
  !!   decaying exponentially for droplet mass x.
  !!   It can easily be reformulated for f(x)=A*x**(nuc)*exp(-Bx**(mu)) and
  !!   by chosing mu=1/3 one would get a gamma distribution in drop diameter
  !!   -> faster rain formation. (Seifert)
  subroutine autoconversion
    use modglobal, only : i1,j1,kmax,eps1,rlv,cp
    use modmpi,    only : myid
    implicit none
    au = 0.

    if (l_sb ) then
    !
    ! SB autoconversion
    !
      tau(2:i1,2:j1,1:k1) = 0.
      phi(2:i1,2:j1,1:k1) = 0.
      nuc(2:i1,2:j1,1:k1) = 0.
      k_au = k_c/(20*x_s)

      where (qcmask(2:i1,2:j1,1:k1))

        nuc(2:i1,2:j1,1:k1) = 1.58*(rhoz(2:i1,2:j1,1:k1)*qc(2:i1,2:j1,1:k1)*1000.) +0.72-1. !G09a
!        nuc(2:i1,2:j1,1:k1) = 0. !
        xc(2:i1,2:j1,1:k1) = rhoz(2:i1,2:j1,1:k1)*qc(2:i1,2:j1,1:k1)/Nc(2:i1,2:j1,1:k1)
        au(2:i1,2:j1,1:k1) = k_au* (nuc(2:i1,2:j1,1:k1)+2.)*(nuc(2:i1,2:j1,1:k1)+4.)/(nuc(2:i1,2:j1,1:k1)+1.)**2.    &
                    * (qc(2:i1,2:j1,1:k1) * xc(2:i1,2:j1,1:k1))**2. * 1.225 ! *rho**2/rho/rho (= 1)
        tau(2:i1,2:j1,1:k1) = 1.0 - qc(2:i1,2:j1,1:k1)/qltot(2:i1,2:j1,1:k1)
        phi(2:i1,2:j1,1:k1) = k_1 * tau(2:i1,2:j1,1:k1)**k_2 * (1.0 -tau(2:i1,2:j1,1:k1)**k_2)**3
        au(2:i1,2:j1,1:k1) = au(2:i1,2:j1,1:k1) * (1.0 + phi(2:i1,2:j1,1:k1)/(1.0 -tau(2:i1,2:j1,1:k1))**2)

        qrp(2:i1,2:j1,1:k1) = qrp(2:i1,2:j1,1:k1) + au(2:i1,2:j1,1:k1)
        Nrp(2:i1,2:j1,1:k1) = Nrp(2:i1,2:j1,1:k1) + au(2:i1,2:j1,1:k1)/x_s
        qtpmcr(2:i1,2:j1,1:k1) = qtpmcr(2:i1,2:j1,1:k1) - au(2:i1,2:j1,1:k1)
        thlpmcr(2:i1,2:j1,1:k1) = thlpmcr(2:i1,2:j1,1:k1) + (rlv/(cp*exnz(2:i1,2:j1,1:k1)))*au(2:i1,2:j1,1:k1)
      endwhere

    else
    !
    ! KK00 autoconversion
    !
      where (qcmask(2:i1,2:j1,1:k1))
        au(2:i1,2:j1,1:k1) = 1350.0 * qc(2:i1,2:j1,1:k1)**(2.47) * (Nc(2:i1,2:j1,1:k1)/1.0E6)**(-1.79)

        qrp(2:i1,2:j1,1:k1) = qrp(2:i1,2:j1,1:k1) + au(2:i1,2:j1,1:k1)
        Nrp(2:i1,2:j1,1:k1) = Nrp(2:i1,2:j1,1:k1) + au(2:i1,2:j1,1:k1) * rhoz(2:i1,2:j1,1:k1)/(pirhow*D0_kk**3.)
        qtpmcr(2:i1,2:j1,1:k1) = qtpmcr(2:i1,2:j1,1:k1) - au(2:i1,2:j1,1:k1)
        thlpmcr(2:i1,2:j1,1:k1) = thlpmcr(2:i1,2:j1,1:k1) + (rlv/(cp*exnz(2:i1,2:j1,1:k1)))*au(2:i1,2:j1,1:k1)
      endwhere

    end if !l_sb


    if (any(qc(2:i1,2:j1,1:kmax)/delt - au(2:i1,2:j1,1:kmax) .lt. 0.)) then
      write(6,*)'au too large', count(qc(2:i1,2:j1,1:kmax)/delt - au(2:i1,2:j1,1:kmax) .lt. 0.),myid
    end if

  end subroutine autoconversion

  subroutine accretion
  !*********************************************************************
  ! determine accr. + self coll. + br-up rate and adjust qrp and Nrp
  ! accordingly. Break-up : Seifert (2007)
  !*********************************************************************
    use modglobal, only : i1,j1,kmax,eps1,rlv,cp,dt,dzf
    use modfields, only : rhof
    use modmpi,    only : myid
    implicit none
    real :: phi_br(2-ih:i1+ih,2-jh:j1+jh,k1)

    ac(2:i1,2:j1,1:k1)=0

    if (l_sb ) then
    !
    ! SB accretion
    !
     where (qrmask(2:i1,2:j1,1:k1) .and. qcmask(2:i1,2:j1,1:k1))
       tau(2:i1,2:j1,1:k1) = 1.0 - qc(2:i1,2:j1,1:k1)/(qltot(2:i1,2:j1,1:k1))
       phi(2:i1,2:j1,1:k1) = (tau(2:i1,2:j1,1:k1)/(tau(2:i1,2:j1,1:k1) + k_l))**4.
       ac(2:i1,2:j1,1:k1) = k_r *rhoz(2:i1,2:j1,1:k1)*qc(2:i1,2:j1,1:k1) * qr(2:i1,2:j1,1:k1) * phi(2:i1,2:j1,1:k1) * (1.225/rhoz(2:i1,2:j1,1:k1))**0.5

       qrp(2:i1,2:j1,1:k1) = qrp(2:i1,2:j1,1:k1) + ac(2:i1,2:j1,1:k1)
       qtpmcr(2:i1,2:j1,1:k1) = qtpmcr(2:i1,2:j1,1:k1) - ac(2:i1,2:j1,1:k1)
       thlpmcr(2:i1,2:j1,1:k1) = thlpmcr(2:i1,2:j1,1:k1) + (rlv/(cp*exnz(2:i1,2:j1,1:k1)))*ac(2:i1,2:j1,1:k1)
     endwhere
    !
    ! SB self-collection & Break-up
    !
      sc(2:i1,2:j1,1:k1)=0
      br(2:i1,2:j1,1:k1)=0
      where (qrmask(2:i1,2:j1,1:k1))
        sc(2:i1,2:j1,1:k1) = k_rr *rhoz(2:i1,2:j1,1:k1)* qr(2:i1,2:j1,1:k1) * Nr(2:i1,2:j1,1:k1)  &
                       * (1 + kappa_r/lbdr(2:i1,2:j1,1:k1))**(-9.)* (1.225/rhoz(2:i1,2:j1,1:k1))**0.5

      endwhere
      where (  Dvr(2:i1,2:j1,1:k1) .gt. 0.30E-3 .and. qrmask(2:i1,2:j1,1:k1))
        phi_br(2:i1,2:j1,1:k1) = k_br * (Dvr(2:i1,2:j1,1:k1)-D_eq)
        br(2:i1,2:j1,1:k1) = (phi_br(2:i1,2:j1,1:k1) + 1.) * sc(2:i1,2:j1,1:k1)
      elsewhere
        br(2:i1,2:j1,1:k1) = 0. ! (phi_br = -1)
      endwhere

      Nrp(2:i1,2:j1,1:k1) = Nrp(2:i1,2:j1,1:k1) - sc(2:i1,2:j1,1:k1) + br(2:i1,2:j1,1:k1)

    else
    !
    ! KK00 accretion
    !
      where (qrmask(2:i1,2:j1,1:k1) .and. qcmask(2:i1,2:j1,1:k1))
        ac(2:i1,2:j1,1:k1) = 67.0 * ( qc(2:i1,2:j1,1:k1) * qr(2:i1,2:j1,1:k1) )**1.15

        qrp(2:i1,2:j1,1:k1) = qrp(2:i1,2:j1,1:k1) + ac(2:i1,2:j1,1:k1)
        qtpmcr(2:i1,2:j1,1:k1) = qtpmcr(2:i1,2:j1,1:k1) - ac(2:i1,2:j1,1:k1)
        thlpmcr(2:i1,2:j1,1:k1) = thlpmcr(2:i1,2:j1,1:k1) + (rlv/(cp*exnz(2:i1,2:j1,1:k1)))*ac(2:i1,2:j1,1:k1)
     endwhere

    end if !l_sb


   if (any(qc(2:i1,2:j1,1:kmax)/delt - ac(2:i1,2:j1,1:kmax) .lt. 0.)) then
     write(6,*)'ac too large', count(qc(2:i1,2:j1,1:kmax)/delt - ac(2:i1,2:j1,1:kmax) .lt. 0.),myid
   end if

  end subroutine accretion

!> Sedimentation of cloud water
!!
!!   The sedimentation of cloud droplets assumes a lognormal DSD in which the
!!   geometric std dev. is assumed to be fixed at 1.3.
!! sedimentation of cloud droplets
!! lognormal CDSD is assumed (1 free parameter : sig_g)
!! terminal velocity : Stokes velocity is assumed (v(D) ~ D^2)
!! flux is calc. anal.
  subroutine sedimentation_cloud
    use modglobal, only : i1,j1,k1,kmax,eps1,rlv,cp,dzf,pi
    use modmpi,    only : myid

    implicit none
    integer :: k

!    real    :: qc_spl(2-ih:i1+ih,2-jh:j1+jh,k1)       &! work variable
!              ,Nc_spl(2-ih:i1+ih,2-jh:j1+jh,k1)
!    real,save :: dt_spl,wfallmax
!
!    qc_spl(2:i1,2:j1,1:k1) = qc(2:i1,2:j1,1:k1)
!    Nc_spl(2:i1,2:j1,1:k1)  = Nc(2:i1,2:j1,1:k1)
!
!    wfallmax = 9.9
!    n_spl = ceiling(wfallmax*delt/(minval(dzf)))
!    dt_spl = delt/real(n_spl)
!
!    do jn = 1 , n_spl  ! time splitting loop

    sedc(2:i1,2:j1,1:k1) = 0.
    csed = c_St*(3./(4.*pi*rhow))**(2./3.)*exp(5.*log(sig_g)**2.)
    where (qcmask(2:i1,2:j1,1:k1))
      sedc(2:i1,2:j1,1:k1) = csed*(Nc(2:i1,2:j1,1:k1))**(-2./3.)*(qc(2:i1,2:j1,1:k1)*rhoz(2:i1,2:j1,1:k1))**(5./3.)
    endwhere
    do k = 1,kmax
      qtpmcr(2:i1,2:j1,k) = qtpmcr(2:i1,2:j1,k) + (sedc(2:i1,2:j1,k+1)-sedc(2:i1,2:j1,k))/(dzf(k)*rhoz(2:i1,2:j1,k))
     thlpmcr(2:i1,2:j1,k) = thlpmcr(2:i1,2:j1,k) - (rlv/(cp*exnz(2:i1,2:j1,k))) &
                                    *(sedc(2:i1,2:j1,k+1)-sedc(2:i1,2:j1,k))/(dzf(k)*rhoz(2:i1,2:j1,k))
    enddo

  end subroutine sedimentation_cloud


!> Sedimentaion of rain
!! sedimentation of drizzle water
!! - gen. gamma distr is assumed. Terminal velocities param according to
!!   Stevens & Seifert. Flux are calc. anal.
!! - l_lognormal =T : lognormal DSD is assumed with D_g and N known and
!!   sig_g assumed. Flux are calc. numerically with help of a
!!   polynomial function
  subroutine sedimentation_rain
    use modglobal, only : i1,j1,k1,kmax,eps1,dzf,pi,dt
    use modfields, only : rhof
    use modmpi,    only : myid,mpi_max,mpi_integer,mpierr,comm3d
    implicit none
    integer :: i,j,k,jn
    integer :: n_spl      !<  sedimentation time splitting loop
    real    :: pwcont
    real    :: wvar(2-ih:i1+ih,2-jh:j1+jh,k1)       &!<  work variable
              ,xr_spl(2-ih:i1+ih,2-jh:j1+jh,k1)     &!<  for time splitting
              ,Dvr_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    &!<     -
              ,mur_spl(2-ih:i1+ih,2-jh:j1+jh,k1)    &!<     -
              ,lbdr_spl(2-ih:i1+ih,2-jh:j1+jh,k1)   &!<     -
              ,Dgr(2-ih:i1+ih,2-jh:j1+jh,k1)         !<  lognormal geometric diameter
    real,save :: dt_spl,wfallmax

    qr_spl(2:i1,2:j1,1:k1) = qr(2:i1,2:j1,1:k1)
    Nr_spl(2:i1,2:j1,1:k1)  = Nr(2:i1,2:j1,1:k1)

    wfallmax = 9.9
    n_spl = ceiling(wfallmax*delt/(minval(dzf)))
    dt_spl = delt/real(n_spl)

    do jn = 1 , n_spl ! time splitting loop

      sed_qr(2:i1,2:j1,1:k1) = 0.
      sed_Nr(2:i1,2:j1,1:k1) = 0.

    if (l_sb ) then

      where (qr_spl(2:i1,2:j1,1:k1) > qrmin)
        xr_spl(2:i1,2:j1,1:k1) = rhoz(2:i1,2:j1,1:k1)*qr_spl(2:i1,2:j1,1:k1)/(Nr_spl(2:i1,2:j1,1:k1))
        xr_spl(2:i1,2:j1,1:k1) = min(max(xr_spl(2:i1,2:j1,1:k1),xrmin),xrmax) ! to ensure xr is within borders
        Dvr_spl(2:i1,2:j1,1:k1) = (xr_spl(2:i1,2:j1,1:k1)/pirhow)**(1./3.)
      endwhere

      if (l_lognormal) then
        do j = 2,j1
        do i = 2,i1
        do k = 1,kmax
          if (qr_spl(i,j,k) > qrmin) then
            Dgr(i,j,k) = (exp(4.5*(log(sig_gr))**2))**(-1./3.)*Dvr_spl(i,j,k) ! correction for width of DSD
            sed_qr(i,j,k) = 1.*sed_flux(Nr_spl(i,j,k),Dgr(i,j,k),log(sig_gr)**2,D_s,3)
            sed_Nr(i,j,k) = 1./pirhow*sed_flux(Nr_spl(i,j,k),Dgr(i,j,k) ,log(sig_gr)**2,D_s,0)
!        correction for the fact that pwcont .ne. qr_spl
!        actually in this way for every grid box a fall velocity is determined
            pwcont = liq_cont(Nr_spl(i,j,k),Dgr(i,j,k),log(sig_gr)**2,D_s,3)         ! note : kg m-3
            if (pwcont > eps1) then
              sed_qr(i,j,k) = (qr_spl(i,j,k)*rhoz(i,j,k)/pwcont)*sed_qr(i,j,k)  ! or qr_spl*(sed_qr/pwcont) = qr_spl*fallvel.
            end if
          end if ! qr_spl threshold statement
        end do
        end do
        end do

      else
    !
    ! SB rain sedimentation
    !
        if (l_mur_cst) then
          mur_spl(2:i1,2:j1,1:k1) = mur_cst
        else
          where (qr_spl(2:i1,2:j1,1:k1) > qrmin)
!          mur_spl(2:i1,2:j1,1:k1) = 10. * (1+tanh(1200.*(Dvr_spl(2:i1,2:j1,1:k1)-0.0014))) ! SS08
          mur_spl(2:i1,2:j1,1:k1) = min(30.,- 1. + 0.008/ (qr_spl(2:i1,2:j1,1:k1)*rhoz(2:i1,2:j1,1:k1))**0.6)  ! G09b
          endwhere
        endif

        where (qr_spl(2:i1,2:j1,1:k1) > qrmin)
          lbdr_spl(2:i1,2:j1,1:k1) = ((mur_spl(2:i1,2:j1,1:k1)+3.)*(mur_spl(2:i1,2:j1,1:k1)+2.)*(mur_spl(2:i1,2:j1,1:k1)+1.))**(1./3.)/Dvr_spl(2:i1,2:j1,1:k1)
          wfall_qr(2:i1,2:j1,1:k1)  = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr_spl(2:i1,2:j1,1:k1))**(-1.*(mur_spl(2:i1,2:j1,1:k1)+4.))))
          wfall_Nr(2:i1,2:j1,1:k1) = max(0.,(a_tvsb-b_tvsb*(1.+c_tvsb/lbdr_spl(2:i1,2:j1,1:k1))**(-1.*(mur_spl(2:i1,2:j1,1:k1)+1.))))
          sed_qr(2:i1,2:j1,1:k1) = wfall_qr(2:i1,2:j1,1:k1)*qr_spl(2:i1,2:j1,1:k1)*rhoz(2:i1,2:j1,1:k1)
          sed_Nr(2:i1,2:j1,1:k1) = wfall_Nr(2:i1,2:j1,1:k1)*Nr_spl(2:i1,2:j1,1:k1)
        endwhere
       endif !l_lognormal

    else
    !
    ! KK00 rain sedimentation
    !
      where (qr_spl(2:i1,2:j1,1:k1) > qrmin)
        xr_spl(2:i1,2:j1,1:k1) = rhoz(2:i1,2:j1,1:k1)*qr_spl(2:i1,2:j1,1:k1)/(Nr_spl(2:i1,2:j1,1:k1))
        xr_spl(2:i1,2:j1,1:k1) = min(xr_spl(2:i1,2:j1,1:k1),xrmaxkk) ! to ensure xr is within borders
        Dvr_spl(2:i1,2:j1,1:k1) = (xr_spl(2:i1,2:j1,1:k1)/pirhow)**(1./3.)
        sed_qr(2:i1,2:j1,1:k1) = max(0., 0.006*1.0E6*Dvr_spl(2:i1,2:j1,1:k1)- 0.2) * qr_spl(2:i1,2:j1,1:k1)*rhoz(2:i1,2:j1,1:k1)
        sed_Nr(2:i1,2:j1,1:k1) = max(0.,0.0035*1.0E6*Dvr_spl(2:i1,2:j1,1:k1)- 0.1) * Nr_spl(2:i1,2:j1,1:k1)
      endwhere

    end if !l_sb
!
    do k = 1,kmax

      wvar(2:i1,2:j1,k) = qr_spl(2:i1,2:j1,k) + (sed_qr(2:i1,2:j1,k+1) - sed_qr(2:i1,2:j1,k))*dt_spl/(dzf(k)*rhoz(2:i1,2:j1,k))
      if (any(wvar(2:i1,2:j1,k) .lt. 0.)) then
        write(6,*)'sed qr too large', count(wvar(2:i1,2:j1,k) .lt. 0.),myid, minval(wvar), minloc(wvar)
      end if

      Nr_spl(2:i1,2:j1,k) = Nr_spl(2:i1,2:j1,k) + (sed_Nr(2:i1,2:j1,k+1) - sed_Nr(2:i1,2:j1,k))*dt_spl/dzf(k)
      qr_spl(2:i1,2:j1,k) = qr_spl(2:i1,2:j1,k) + (sed_qr(2:i1,2:j1,k+1) - sed_qr(2:i1,2:j1,k))*dt_spl/(dzf(k)*rhoz(2:i1,2:j1,k))
      if ( jn == 1. ) then
        precep(2:i1,2:j1,k) =  sed_qr(2:i1,2:j1,k)/rhoz(2:i1,2:j1,k)   ! kg kg-1 m s-1
      endif

    end do  ! second k loop
!
    enddo ! time splitting loop

    Nrp(2:i1,2:j1,1:k1)= Nrp(2:i1,2:j1,1:k1) + (Nr_spl(2:i1,2:j1,1:k1) - Nr(2:i1,2:j1,1:k1))/delt
    qrp(2:i1,2:j1,1:k1)= qrp(2:i1,2:j1,1:k1) + (qr_spl(2:i1,2:j1,1:k1) - qr(2:i1,2:j1,1:k1))/delt

  end subroutine sedimentation_rain

  !*********************************************************************
  !*********************************************************************

  subroutine evaporation
  !*********************************************************************
  ! Evaporation of prec. : Seifert (2008)
  ! Cond. (S>0.) neglected (all water is condensed on cloud droplets)
  !*********************************************************************
    use modglobal, only : i1,j1,kmax,eps1,es0,rd,rv,tmelt,rlv,cp,at,bt,pi,ep
    use modfields, only : exnf,thl0,qt0,svm
    use modmpi,    only : myid
    implicit none
    real    :: F(2-ih:i1+ih,2-jh:j1+jh,k1)     & ! ventilation factor
              ,S(2-ih:i1+ih,2-jh:j1+jh,k1)     & ! super or undersaturation
              ,G(2-ih:i1+ih,2-jh:j1+jh,k1)     & ! cond/evap rate of a drop
              ,T(2-ih:i1+ih,2-jh:j1+jh,k1)     & ! temperature
              ,esat(2-ih:i1+ih,2-jh:j1+jh,k1)  & ! sat. pressure
              ,qsat(2-ih:i1+ih,2-jh:j1+jh,k1)  & ! sat. water vapor
              ,tl(2-ih:i1+ih,2-jh:j1+jh,k1)    & ! T_l
              ,b1(2-ih:i1+ih,2-jh:j1+jh,k1)    &
              ,qsl(2-ih:i1+ih,2-jh:j1+jh,k1)

    evap(2:i1,2:j1,1:k1) = 0.
    Nevap(2:i1,2:j1,1:k1) = 0.

    where ( qrmask(2:i1,2:j1,1:k1))
      tl(2:i1,2:j1,1:k1) = thl0(2:i1,2:j1,1:k1)*exnz(2:i1,2:j1,1:k1)
      esat(2:i1,2:j1,1:k1)  = es0*exp(at*(tl(2:i1,2:j1,1:k1)-tmelt)/(tl(2:i1,2:j1,1:k1)-bt))
      qsl(2:i1,2:j1,1:k1) = ep*esat(2:i1,2:j1,1:k1)/(presz(2:i1,2:j1,1:k1)-(1-ep)*esat(2:i1,2:j1,1:k1))
      b1(2:i1,2:j1,1:k1)  = rlv**2/(cp*rv*tl(2:i1,2:j1,1:k1)**2)
      qsat(2:i1,2:j1,1:k1)  = qsl(2:i1,2:j1,1:k1)*(1.+b1(2:i1,2:j1,1:k1)*qt0(2:i1,2:j1,1:k1))/(1.+b1(2:i1,2:j1,1:k1)*qsl(2:i1,2:j1,1:k1))
      S(2:i1,2:j1,1:k1) =  min(0.,(qt0(2:i1,2:j1,1:k1) - qc(2:i1,2:j1,1:k1))/qsat(2:i1,2:j1,1:k1) - 1.)
      T(2:i1,2:j1,1:k1) = tl(2:i1,2:j1,1:k1) + (rlv/cp)*qc(2:i1,2:j1,1:k1)
      G(2:i1,2:j1,1:k1) = (Rv*T(2:i1,2:j1,1:k1))/(Dv*esat(2:i1,2:j1,1:k1)) + rlv/(Kt*T(2:i1,2:j1,1:k1))*(rlv/(Rv*T(2:i1,2:j1,1:k1)) -1.)
      G(2:i1,2:j1,1:k1) = 1./G(2:i1,2:j1,1:k1)
    end where
    if (l_sb ) then
      where ( qrmask(2:i1,2:j1,1:k1))
        F(2:i1,2:j1,1:k1) = avf*f_gamma(mur(2:i1,2:j1,1:k1)+2.)/lbdr(2:i1,2:j1,1:k1)**(mur(2:i1,2:j1,1:k1)+2.)    &
                  +                                                                  &
                   bvf*Sc_num**(1./3.)*(a_tvsb/nu_a)**0.5*f_gamma(mur(2:i1,2:j1,1:k1)+2.5)/lbdr(2:i1,2:j1,1:k1)**(mur(2:i1,2:j1,1:k1)+2.5) * &
                   (1.-(1./2.)  *(b_tvsb/a_tvsb)    *(lbdr(2:i1,2:j1,1:k1)/(   c_tvsb+lbdr(2:i1,2:j1,1:k1)))**(mur(2:i1,2:j1,1:k1)+2.5)  &
                      -(1./8.)  *(b_tvsb/a_tvsb)**2.*(lbdr(2:i1,2:j1,1:k1)/(2.*c_tvsb+lbdr(2:i1,2:j1,1:k1)))**(mur(2:i1,2:j1,1:k1)+2.5)  &
                      -(1./16.) *(b_tvsb/a_tvsb)**3.*(lbdr(2:i1,2:j1,1:k1)/(3.*c_tvsb+lbdr(2:i1,2:j1,1:k1)))**(mur(2:i1,2:j1,1:k1)+2.5)  &
                      -(5./128.)*(b_tvsb/a_tvsb)**4.*(lbdr(2:i1,2:j1,1:k1)/(4.*c_tvsb+lbdr(2:i1,2:j1,1:k1)))**(mur(2:i1,2:j1,1:k1)+2.5)  )
        evap(2:i1,2:j1,1:k1) = 2*pi*(Nr(2:i1,2:j1,1:k1)*lbdr(2:i1,2:j1,1:k1)**(mur(2:i1,2:j1,1:k1)+1.)/f_gamma(mur(2:i1,2:j1,1:k1)+1.)) &
                          *G(2:i1,2:j1,1:k1)*F(2:i1,2:j1,1:k1)*S(2:i1,2:j1,1:k1)/rhoz(2:i1,2:j1,1:k1)
        Nevap(2:i1,2:j1,1:k1) = c_Nevap*evap(2:i1,2:j1,1:k1)*rhoz(2:i1,2:j1,1:k1)/xr(2:i1,2:j1,1:k1)
      endwhere
    else
      where ( qrmask(2:i1,2:j1,1:k1))
        evap(2:i1,2:j1,1:k1) = c_evapkk*2*pi*Dvr(2:i1,2:j1,1:k1)*G(2:i1,2:j1,1:k1)*S(2:i1,2:j1,1:k1)*Nr(2:i1,2:j1,1:k1)/rhoz(2:i1,2:j1,1:k1)
        Nevap(2:i1,2:j1,1:k1) = evap(2:i1,2:j1,1:k1)*rhoz(2:i1,2:j1,1:k1)/xr(2:i1,2:j1,1:k1)
      end where
    endif

    where (evap(2:i1,2:j1,1:k1) < -svm(2:i1,2:j1,1:k1,2)/delt .and. qrmask(2:i1,2:j1,1:k1))   ! total evap of drops and qr
      Nevap(2:i1,2:j1,1:k1) = - svm(2:i1,2:j1,1:k1,1)/delt
      evap(2:i1,2:j1,1:k1) = - svm(2:i1,2:j1,1:k1,2)/delt
    endwhere
    qrp(2:i1,2:j1,1:k1) = qrp(2:i1,2:j1,1:k1) + evap(2:i1,2:j1,1:k1)
    Nrp(2:i1,2:j1,1:k1) = Nrp(2:i1,2:j1,1:k1) + Nevap(2:i1,2:j1,1:k1)
    qtpmcr(2:i1,2:j1,1:k1) = qtpmcr(2:i1,2:j1,1:k1) -evap(2:i1,2:j1,1:k1)
    thlpmcr(2:i1,2:j1,1:k1) = thlpmcr(2:i1,2:j1,1:k1) + (rlv/cp)*evap(2:i1,2:j1,1:k1)

  end subroutine evaporation

  !*********************************************************************
  !*********************************************************************

  function f_gamma(xx) result(gam)

  !*********************************************************************
  ! Gamma function
  !
  ! http://pagesperso-orange.fr/jean-pierre.moreau/Fortran/gamma_f90.txt
  ! Reference:
  ! "Numerical Recipes, by W.H. Press, B.P. Flannery, S.A. Teukolsky
  !  and T. Vetterling, Cambridge University Press, 1986" [BIBLI 08].
  ! --------------------------------------------------------------- *
  ! Returns the value of Gamma(x) for X>0.
  !
  ! O. Geoffroy february 2008
  !*********************************************************************
    use modglobal, only : ih,i1,jh,j1,k1
    implicit none
    real :: xx(2:i1,2:j1,k1),  &
            x(2:i1,2:j1,k1),  &
            tmp(2:i1,2:j1,k1),  &
            ser(2:i1,2:j1,k1),  &
            gam(2:i1,2:j1,k1)


    real cof(6),stp,half,one,fpf

  integer j
  data cof,stp /76.18009173,-86.50532033,24.01409822,  &
       -1.231739516,0.120858003e-2,-0.536382e-5,2.50662827465/
  data half,one,fpf /0.5,1.0,5.5/

  x(2:i1,2:j1,1:k1)=xx(2:i1,2:j1,1:k1)-one
  tmp(2:i1,2:j1,1:k1)=x(2:i1,2:j1,1:k1)+fpf
  tmp(2:i1,2:j1,1:k1)=(x(2:i1,2:j1,1:k1)+half)*log(tmp(2:i1,2:j1,1:k1))-tmp(2:i1,2:j1,1:k1)
  ser(2:i1,2:j1,1:k1)=one
  do j=1,6
    where (qrmask(2:i1,2:j1,1:k1))
      x(2:i1,2:j1,1:k1)=x(2:i1,2:j1,1:k1)+one
      ser(2:i1,2:j1,1:k1)=ser(2:i1,2:j1,1:k1)+cof(j)/x(2:i1,2:j1,1:k1)
    endwhere
  end do

  gam(2:i1,2:j1,1:k1) = exp(tmp(2:i1,2:j1,1:k1)+log(stp*ser(2:i1,2:j1,1:k1)))


  end function f_gamma

  !*********************************************************************
  !*********************************************************************

  real function sed_flux(Nin,Din,sig2,Ddiv,nnn)

  !*********************************************************************
  ! Function to calculate numerically the analytical solution of the
  ! sedimentation flux between Dmin and Dmax based on
  ! Feingold et al 1986 eq 17 -20.
  ! fall velocity is determined by alfa* D^beta with alfa+ beta taken as
  ! specified in Rogers and Yau 1989 Note here we work in D and in SI
  ! (in Roger+Yau in cm units + radius)
  ! flux is multiplied outside sed_flux with 1/rho_air to get proper
  ! kg/kg m/s units
  !
  ! M.C. van Zanten    August 2005
  !*********************************************************************
    use modglobal, only : pi,rhow
    implicit none

    real, intent(in) :: Nin, Din, sig2, Ddiv
    integer, intent(in) :: nnn
    !para. def. lognormal DSD (sig2 = ln^2 sigma_g), D sep. droplets from drops
    !,power of of D in integral

    real, parameter ::   C = rhow*pi/6.     &
                        ,D_intmin = 1e-6    &
                        ,D_intmax = 4.3e-3

    real ::  alfa         & ! constant in fall velocity relation
            ,beta         & ! power in fall vel. rel.
            ,D_min        & ! min integration limit
            ,D_max        & ! max integration limit
            ,flux           ![kg m^-2 s^-1]

    integer :: k

    flux = 0.0

    if (Din < Ddiv) then
      alfa = 3.e5*100  ![1/ms]
      beta = 2
      D_min = D_intmin
      D_max = Ddiv
      flux = C*Nin*alfa*erfint(beta,Din,D_min,D_max,sig2,nnn)
    else
      do k = 1,3
        select case(k)
        case(1)        ! fall speed ~ D^2
          alfa = 3.e5*100 ![1/m 1/s]
          beta = 2
          D_min = Ddiv
          D_max = 133e-6
        case(2)        ! fall speed ~ D
          alfa = 4e3     ![1/s]
          beta = 1
          D_min = 133e-6
          D_max = 1.25e-3
        case default         ! fall speed ~ sqrt(D)
          alfa = 1.4e3 *0.1  ![m^.5 1/s]
          beta = .5
          D_min = 1.25e-3
          D_max = D_intmax
        end select
        flux = flux + C*Nin*alfa*erfint(beta,Din,D_min,D_max,sig2,nnn)
      end do
    end if
      sed_flux = flux
  end function sed_flux

  !*********************************************************************
  !*********************************************************************

  real function liq_cont(Nin,Din,sig2,Ddiv,nnn)

  !*********************************************************************
  ! Function to calculate numerically the analytical solution of the
  ! liq. water content between Dmin and Dmax based on
  ! Feingold et al 1986 eq 17 -20.
  !
  ! M.C. van Zanten    September 2005
  !*********************************************************************
    use modglobal, only : pi,rhow
    implicit none

    real, intent(in) :: Nin, Din, sig2, Ddiv
    integer, intent(in) :: nnn
    !para. def. lognormal DSD (sig2 = ln^2 sigma_g), D sep. droplets from drops
    !,power of of D in integral

    real, parameter :: beta = 0           &
                      ,C = pi/6.*rhow     &
                      ,D_intmin = 80e-6    &   ! value of start of rain D
                      ,D_intmax = 3e-3         !4.3e-3    !  value is now max value for sqrt fall speed rel.

    real ::  D_min        & ! min integration limit
            ,D_max          ! max integration limit

    if (Din < Ddiv) then
    D_min = D_intmin
    D_max = Ddiv
    else
    D_min = Ddiv
    D_max = D_intmax
    end if

    liq_cont = C*Nin*erfint(beta,Din,D_min,D_max,sig2,nnn)

  end function liq_cont

  !*********************************************************************
  !*********************************************************************

  real function erfint(beta, D, D_min, D_max, sig2,nnn )

  !*********************************************************************
  ! Function to calculate erf(x) approximated by a polynomial as
  ! specified in 7.1.27 in Abramowitz and Stegun
  ! NB phi(x) = 0.5(erf(0.707107*x)+1) but 1 disappears by substraction
  !
  !*********************************************************************
    implicit none
    real, intent(in) :: beta, D, D_min, D_max, sig2
    integer, intent(in) :: nnn

    real, parameter :: eps = 1e-10       &
                      ,a1 = 0.278393    & !a1 till a4 constants in polynomial fit to the error
                      ,a2 = 0.230389    & !function 7.1.27 in Abramowitz and Stegun
                      ,a3 = 0.000972    &
                      ,a4 = 0.078108
    real :: nn, ymin, ymax, erfymin, erfymax, D_inv

    D_inv = 1./(eps + D)
    nn = beta + nnn

    ymin = 0.707107*(log(D_min*D_inv) - nn*sig2)/(sqrt(sig2))
    ymax = 0.707107*(log(D_max*D_inv) - nn*sig2)/(sqrt(sig2))

    erfymin = 1.-1./((1.+a1*abs(ymin) + a2*abs(ymin)**2 + a3*abs(ymin)**3 +a4*abs(ymin)**4)**4)
    erfymax = 1.-1./((1.+a1*abs(ymax) + a2*abs(ymax)**2 + a3*abs(ymax)**3 +a4*abs(ymax)**4)**4)
    if (ymin < 0.) then
      erfymin = -1.*erfymin
    end if
    if (ymax < 0.) then
      erfymax = -1.*erfymax
    end if
    erfint = D**nn*exp(0.5*nn**2*sig2)*0.5*(erfymax-erfymin)
  !  if (erfint < 0.) write(*,*)'erfint neg'
    if (erfint < 0.) erfint = 0.
  end function erfint



end module modbulkmicro


