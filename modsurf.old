!----------------------------------------------------------------------------
! This file is part of DALES.
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
! Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!----------------------------------------------------------------------------
!
!

!       ----------------------------------------------------------------
!*    module *modsurf* relevant variables for the surface layer
!       ----------------------------------------------------------------

module modsurf
implicit none
save

contains
  subroutine initsurf
    use modglobal, only : i2,j2,nsv
    use modsurfdata, only : ustar,dudz,dvdz,tstar,qstar,dqtdz,dthldz,svstar,svs
    implicit none

    allocate(ustar (i2,j2))
    allocate(dudz  (i2,j2))
    allocate(dvdz  (i2,j2))
    allocate(tstar (i2,j2))
    allocate(qstar (i2,j2))
    allocate(dqtdz (i2,j2))
    allocate(dthldz(i2,j2))
    allocate(svstar(i2,j2,nsv))
    allocate(svs(nsv))

  end subroutine initsurf

  subroutine surf
    use modsurfdata, only : isurf
    implicit none

    select case (isurf)
    case(1:2)
      call surface
    case(3:4)
      call surflux
    end select
  end subroutine surf

  subroutine exitsurf
    use modsurfdata, only : ustar,dudz,dvdz,tstar,qstar,dqtdz,dthldz,svstar,svs
    implicit none
    deallocate(ustar,dudz,dvdz,tstar,qstar,dqtdz,dthldz,svstar,svs)
  end subroutine exitsurf

  subroutine surface

!-----------------------------------------------------------------|
!                                                                 |
!*** *surface*  calculates surface layer parameters               |
!                                                                 |
!      Hans Cuijpers   I.M.A.U.     06/01/1995                    |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!    calculation of ustar, tstar, qstar,                          |
!    the obukhov length and the gradients of the                  |
!    prognostic fields. The results are used to                   |
!    calculate the tendencies at the lowest full level.           |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!     *surface* is called from *program*.                         |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal, only : i1,i2,j1,j2,nsv,timee,zf,&
                        fkar,rv,rd,grav,cu, cv,rslabs
  use modfields, only : u0,v0,thl0,qt0,sv0,u0av,v0av,qt0av,thl0av
  use modmpi,    only : myid,comm3d,my_real, mpierr, mpi_sum, excj
  use modsurfdata, only : ustar,tstar,qstar,dudz,dvdz,dqtdz,dthldz,&
                          svs,svstar,z0,qts,thls,thvs,isurf,lneutraldrag, obl
  implicit none


  integer :: i, j, n, itel
  real ::  chku, chkt, chkq
  real ::  psihz0, psiqz0, psimz0
  real ::  phihz0, phiqz0, phimz0
  real ::  phihzf, phimzf
  real ::  psihzf, psimzf
  real ::  z0m,    z0h,    z0q
  real ::  usold, tsold, qsold
  real ::  ust,   tst,   qst
  real ::  upcu, vpcv
  real ::  horv, stab
  real ::  am,    ah,   aq
  real ::  alpha, rnu, eps
  real ::  dthz1, dqz1, dsvz1
  real :: obli

!**** CONSTANTS ****

  alpha = 0.018
  am    = 0.11
  ah    = 0.40
  aq    = 0.62
  rnu   = 1.5e-5
  eps   = 1.e-2
  itel  = 1

  chku  = 10*eps
  chkt  = 10*eps
  chkq  = 10*eps

  if (timee==0.0) then

    dthz1 = thl0av(1)-thls
    dqz1  = qt0av(1)-qts
    horv  = sqrt(u0av(1)**2 + v0av(1)**2)

    ust = horv *fkar/log(zf(1)/z0) + eps
    tst = dthz1*fkar/log(zf(1)/z0)
    qst = dqz1 *fkar/log(zf(1)/z0)

  else

! Mathieu MPI added

! ASSUMPTION: 2,2 is on PE 0!!!!
     if(myid==0)then
       ust = ustar(2,2)
       tst = tstar(2,2)
       qst = qstar(2,2)
     endif

    call MPI_BCAST(ust  ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(tst  ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(qst  ,1,MY_REAL   ,0,comm3d,mpierr)

  endif

!***********************************************************************
!***  Calculate ust, tst, qst and obukhov-length iteratively   *********
!***********************************************************************

  dthz1 = thl0av(1)-thls
  dqz1  = qt0av(1)-qts
  horv  = sqrt(u0av(1)**2 + v0av(1)**2)

  z0m   = z0
  z0h   = z0
  z0q   = z0
  psimz0 = 0.
  psihz0 = 0.
  psiqz0 = 0.
  psimzf = 0.
  psihzf = 0.
  phimzf = 0.
  phihzf = 0.

  if (lneutraldrag) then

    if (isurf==1) then
      z0m = am*rnu/ust + alpha*ust**2/grav
      z0h = ah*rnu/ust
      z0q = aq*rnu/ust
    endif

    ust = fkar * horv / (log(zf(1)/z0m))
    tst = fkar * dthz1/ (log(zf(1)/z0h))
    qst = fkar * dqz1 / (log(zf(1)/z0q))

    psimz0 = 0.
    psihz0 = 0.
    psiqz0 = 0.
    psimzf = 0.
    psihzf = 0.

    stab   = 0.

  else
    !   obl   = ust**2/(fkar*(grav/thvs)*(tst+(rv/rd-1)*thvs*qst))
!     horv  = sqrt(u0av(1)**2 + v0av(1)**2)
!     dthz1 = thl0av(1)-thls
!     dqz1  = qt0av(1)-qts
    !   stab  = dthz1+(rv/rd-1)*thvs*dqz1
    if(abs(tst)>eps)then
      obl   = ust**2/(fkar*(grav/thvs)*(tst+(rv/rd-1)*thvs*qst))
    else
      obl   = 0.
    end if
    stab = obl

    !CvH loop only over myid = 0
    if(myid==0) then
      do while (chku>eps .or. chkt>eps .or. chkq>eps)
        usold = ust
        tsold = tst
        qsold = qst

        if (isurf==1) then
          z0m = am*rnu/ust + alpha*ust**2/grav
          z0h = ah*rnu/ust
          z0q = aq*rnu/ust
        endif

      !     Calculate Psi-functions for an unstable surface layer
        if (stab < 0.) then

          phimz0 = (1.-16.*z0m/obl)**0.25
          phihz0 = (1.-16.*z0h/obl)**0.50
          phiqz0 = (1.-16.*z0q/obl)**0.50
          phimzf = (1.-16.*zf(1)/obl)**0.25
          phihzf = (1.-16.*zf(1)/obl)**0.50

          psimz0 = 2.*log(0.5+phimz0/2.)+log(0.5+phimz0**2/2.) &
                      -2.*atan(phimz0)+2.*atan(1.)
          psihz0 = 2.*log(0.5+phihz0/2.)
          psiqz0 = 2.*log(0.5+phiqz0/2.)
          psimzf = 2.*log(0.5+phimzf/2.)+log(0.5+phimzf**2/2.) &
                      -2.*atan(phimzf)+2.*atan(1.)
          psihzf = 2.*log(0.5+phihzf/2.)

        endif

      !     End calculations for an unstable surface layer
      !     ------------------------------------------------------
      !     Calculate Psi-functions for a neutral boundary layer

        if (stab == 0.) then

          psimz0 = 0.
          psihz0 = 0.
          psiqz0 = 0.
          psimzf = 0.
          psihzf = 0.

        endif

      !     End calculations for a neutral boundary layer
      !     ------------------------------------------------------
      !     Calculate Psi-functions for a stable surface layer

        if (stab > 0.) then

          psimz0 = -5.*z0m/obl
          psihz0 = -5.*z0h/obl
          psiqz0 = -5.*z0q/obl
          psimzf = -5.*zf(1)/obl
          psihzf = -5.*zf(1)/obl

        endif

      !     End calculations for a stable surface layer
      !     -------------------------------------------------------

        ust = fkar* horv/(log(zf(1)/z0m)-psimzf+psimz0)
        tst = fkar*dthz1/(log(zf(1)/z0h)-psihzf+psihz0)
        qst = fkar* dqz1/(log(zf(1)/z0q)-psihzf+psiqz0)
    !   obl = ust**2/(fkar*(grav/thvs)*(tst+(rv/rd-1)*thvs*qst))
        obli = (fkar*(grav/thvs)*(tst+(rv/rd-1)*thvs*qst))
        if(abs(obli)<1e-3*eps) then
          obli = sign(1e-3*eps,obli)
        end if
        obl = ust**2/obli
        if((ust>eps).and.(tst>eps).and.(qst>eps))then
          chku = abs((ust-usold)/ust)
          chkt = abs((tst-tsold)/tst)
          chkq = abs((qst-qsold)/qst)
        else
          chku = 0
          chkt = 0
          chkq = 0
        end if

        if (itel>10) then
          call stopsurf(obl,ust,tst,qst,chku,chkt,chkq,z0m,z0h,z0q)
        endif
        itel  = itel + 1
      end do

    end if

    ! CvH distribute ust,tst and qst
    call MPI_BCAST(ust  ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(tst  ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(qst  ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(obl  ,1,MY_REAL   ,0,comm3d,mpierr)
    if (abs(obl)<1e-10) obl =0

    stab = obl

  end if

  !***********************************************************************
  !***  Calculation of ustar(i,j),tstar(i,j), qstar(i,j),  ***************
  !***      dudz(i,j),dvdz(i,j), dthldz(i,j), dqtdz(i,j)   ***************
  !***********************************************************************

  if (stab < 0.) then
    phimzf = (1.-16.*zf(1)/obl)**(-0.25)
    phihzf = (1.-16.*zf(1)/obl)**(-0.50)
  endif

  if (stab == 0.) then
   phimzf = 1.
   phihzf = 1.
 endif

  if (stab > 0.) then
    phimzf = (1.+5.*zf(1)/obl)
    phihzf = (1.+5.*zf(1)/obl)
  endif

  do j=2,j1
    do i=2,i1

      dthz1 = thl0(i,j,1)-thls
      dqz1  = qt0(i,j,1)-qts
      upcu  = 0.5*(u0(i,j,1)+u0(i+1,j,1))+cu
      vpcv  = 0.5*(v0(i,j,1)+v0(i,j+1,1))+cv
      horv  = sqrt(upcu**2 + vpcv**2)

      ustar(i,j) = fkar * horv  / (log(zf(1) / z0m) - psimzf + psimz0)
      tstar(i,j) = fkar * dthz1 / (log(zf(1) / z0h) - psihzf + psihz0)
      qstar(i,j) = fkar * dqz1  / (log(zf(1) / z0q) - psihzf + psiqz0)

      dudz(i,j)   = ustar(i,j)*(phimzf/(fkar*zf(1)))*(upcu/horv)
      dvdz(i,j)   = ustar(i,j)*(phimzf/(fkar*zf(1)))*(vpcv/horv)
      dthldz(i,j) = tstar(i,j)*(phihzf/(fkar*zf(1)))
      dqtdz(i,j)  = qstar(i,j)*(phihzf/(fkar*zf(1)))

    end do
  end do

  do j=1,j2
    ustar(1,j)=ustar(i1,j)
    ustar(i2,j)=ustar(2,j)
  end do

  call excj( ustar  , 1, i2, 1, j2, 1,1)

  do n=1,nsv
    do j=2,j1
    do i=2,i1
      dsvz1 = sv0(i,j,1,n) - svs(n)
      svstar(i,j,n) = fkar*dsvz1/(log(zf(1)/z0h)-psihzf+psihz0)
    enddo
    enddo
  enddo

  return

  end subroutine surface



  subroutine surflux

!-----------------------------------------------------------------|
!                                                                 |
!*** *surflux*  calculates surface layer parameters               |
!                                                                 |
!      Hans Cuijpers   I.M.A.U.     06/01/1995                    |
!      Pier Siebesma   K.N.M.I.     11/09/1998                    |
!      Roel Neggers        KNMI       01/09/2000                  |
!      Margreet van Zanten IMAU       01/09/2000                  |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!     calculation of ustar, tstar, qstar,                         |
!     the obukhov length and the gradients of the                 |
!     prognostic fields. The results are used to                  |
!     calculate the tendencies at the lowest full level.          |
!                                                                 |
!     Free convective limit and neutral limit added 01/09/2000    |
!                                                                 |
!     Loops have been modified to make them faster (Stephan)      |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!     if isurf>=3 then *surflux* is called from *program*.        |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal, only : i1,i2,jmax,j2,j1,jtot,nsv,zf,cu,cv,fkar,rv,rd,grav,eps1
  use modfields, only : u0,v0,thl0av,qt0av,sv0av, u0av, v0av
  use modmpi,    only : myid,excj
  use modsurfdata, only : svs,dthldz,dqtdz,dudz,dvdz,ustar,tstar,qstar,svstar,&
                          wtsurf,thls,qts,ps,wqsurf,z0,isurf,ustin,wsvsurf,thvs,&
                          lneutraldrag, obl
  implicit none

  real :: horv,ust,tst,qst,wthvsurf
  real :: c1,c2
  real :: dthldz_all, dqtdz_all
  real :: phihz0inv, phihzfinv ,phimzf,phihzf
  real :: psihz0,psihzf,upcu,vpcv,Cm !, obl
  real, allocatable, dimension (:) :: svst,dsvdz_all
  integer :: i,j,n,iii

  allocate(svst(nsv), dsvdz_all(nsv))

!       1 Determine ustar
!       -----------------
  if (isurf==3) then
    ust = max(ustin,eps1)
  endif
  if (isurf==4) then
    ust = 0.0
    call calust(ust)
    ust = max(ust, eps1)
  endif

!       2 Determine tstar, qstar and sstar
!       ----------------------------------

  tst = -wtsurf/ust
  qst = -wqsurf/ust

  do n=1,nsv
    svst(n)= -wsvsurf(n)/ust
  enddo

!       3 Determine du/dz, dthl/dz, dqt/dz at lowest level
!       --------------------------------------------------

  c1       = 1.+(rv/rd-1)*qts
  c2       = (rv/rd-1)
  thvs     = thls*c1
  wthvsurf = c1*wtsurf+c2*thls*wqsurf

!      3.1. FREE CONVECTION LIMIT
!       --------------------------------------------------

  if (ust == 0.0) then
    if (wthvsurf <= 0.0) then
       stop 'stable & neutral BL with ustar=0 not coded!'
    endif

    !step  flux-gradient relationships according to Garratt, text book 'the ABL'
    !step  note 1: the factor -0.7 is discussed on page 52 for thv
    !step  note 2: the code below assumes that -0.7 is appropriate for q
    !step          and for thl

    obl        = 0.0
    dthldz_all = -0.7*(wtsurf)**(.666) &
                          *(grav/thvs)**(-0.3333)*zf(1)**(-1.3333)
    dqtdz_all  = -0.7*(wqsurf)**(.666) &
                          *(grav/thvs)**(-0.3333)*zf(1)**(-1.3333)
    do n=1,nsv
      dsvdz_all(n)  = -0.7*(wsvsurf(n))**(.666) &
                          *(grav/thvs)**(-0.3333)*zf(1)**(-1.3333)
    enddo

    thls = thl0av(1) - zf(1)*dthldz_all
    qts  = qt0av(1) - zf(1)*dqtdz_all
    do n=1,nsv
      svs(n)  = sv0av(1,n) - zf(1)*dsvdz_all(n)
    enddo
    c1 = 1.+(rv/rd-1)*qts
    c2 = (rv/rd-1)
    wthvsurf =c1*wtsurf+c2*thls*wqsurf
    thvs = thls*c1

    do j=1,j2
    do i=1,i2
      dthldz(i,j) = dthldz_all
      dqtdz(i,j)  = dqtdz_all
      dudz(i,j) = 0.0
      dvdz(i,j) = 0.0
    enddo
    enddo

  else  !stable, neutral or unstable

!      3.2. NEUTRAL SURFACE LAYER
!       --------------------------------------------------

    if (wthvsurf == 0.0) then

      obl    = -1.e20 ! hj: like infinity ...
      phimzf = 1.
      phihzf = 1.

    else   !stable or unstable

!       **** First guess for surface properties
      do iii=1,2

        obl = -ust**3/(fkar*(grav/thvs)*wthvsurf)

!        3.3  UNSTABLE SURFACE LAYER
!         -------------------------------------------------------------

        if (obl < 0.)  then
!                                                      !see Garratt:
          phimzf    = (1.-16.*zf(1)/obl) ** (-0.25)  ! (3.33a)
          phihz0inv = (1.-16.*z0/obl) ** 0.50        ! (3.33b)
          phihzf    = (1.-16.*zf(1)/obl) ** (-0.50)
          phihzfinv = 1./phihzf
          psihz0    = 2.*log(0.5+phihz0inv/2.)     ! (3.39)
          psihzf    = 2.*log(0.5+phihzfinv/2.)

        else   ! (obl > 0) since ustin/=0 thus obl/=0

!        3.4 STABLE SURFACE LAYER
!         -------------------------------------------------------------

          phimzf = (1.+5.*zf(1)/obl)               ! below (3.33)
          phihzf = (1.+5.*zf(1)/obl)               !
          psihz0 = -5.*z0/obl                      ! (3.40)
          psihzf = -5.*zf(1)/obl

        endif

!          thls and qts are computed at z0:

        thls = thl0av(1)-tst/fkar*(log(zf(1)/z0)-psihzf+psihz0)
        qts  = qt0av(1) -qst/fkar*(log(zf(1)/z0)-psihzf+psihz0)

        do n=1,nsv
          svs(n)  = sv0av(1,n) -svst(n)/fkar*(log(zf(1)/z0)-psihzf+psihz0)
        enddo

        c1 = 1.+(rv/rd-1)*qts
        c2 = (rv/rd-1)
        wthvsurf =c1*wtsurf+c2*thls*wqsurf
        thvs = thls*c1

      enddo

      ! CvH - calculate drag coefficients for local ustar
      ! ustar**2  = Cm * U **2
      horv = max(sqrt(u0av(1) **2. + v0av(1) **2.), eps1)
      Cm  =  ust**2 / (horv**2)

      ! CvH - do loop to j2 and i2 instead of j1 and i1
      do j=2,j1
        do i=2,i1
          upcu  = 0.5*(u0(i,j,1)+u0(i+1,j,1))+cu
          vpcv  = 0.5*(v0(i,j,1)+v0(i,j+1,1))+cv
          horv  = max(sqrt(upcu**2 + vpcv**2),eps1)

          ! CvH - Get local values for ustar(i,j)
          ustar (i,j) = Cm ** 0.5 * horv
          tstar (i,j) = -wtsurf / ustar(i,j)
          qstar (i,j) = -wqsurf / ustar(i,j)

          do n=1,nsv
            svstar(i,j,n) = -wsvsurf(n) / ustar(i,j)
          enddo

          dudz  (i,j) = ustar(i,j)*phimzf/(fkar*zf(1))*(upcu/horv)
          dvdz  (i,j) = ustar(i,j)*phimzf/(fkar*zf(1))*(vpcv/horv)
          dthldz(i,j) = tstar(i,j)*phihzf/(fkar*zf(1))
          dqtdz (i,j) = qstar(i,j)*phihzf/(fkar*zf(1))

        enddo
      enddo

      do j=1,j2
        ustar(1,j)=ustar(i1,j)
        ustar(i2,j)=ustar(2,j)
      end do

      call excj( ustar  , 1, i2, 1, j2, 1,1)

    endif  !checks wthvsurf == 0

  endif   !end of if statement that checks ustin==0

  deallocate(svst, dsvdz_all)

  return
  end subroutine surflux

  subroutine calust(ust)

  use modglobal,   only : fkar, grav, rv,rd,pi,zf, eps1
  use modfields,   only : v0av,u0av,thl0av
  use modmpi,      only : my_real,mpierr,comm3d,mpi_sum,myid
  use modsurfdata, only : z0,wqsurf,wtsurf,qts,lneutraldrag
  implicit none

  !CvH Include z0 term in correction
  real obl
  real phimzf, psimzf, phimz0, psimz0
  real c1, c2
  real ust, u1
  real conv, eps, ustold
  real wthvsurf

  psimz0 = 0.0
  u1 = (u0av(1)**2 + v0av(1)**2)**0.5
  c1 =  1.+(rv/rd-1)*qts
  c2 =  (rv/rd-1)
  wthvsurf = c1*wtsurf + c2*thl0av(1)*wqsurf
  wthvsurf = max(1e-20, wthvsurf)
  ust = fkar*u1/log(zf(1)/z0)
  conv = 1.
  eps  = 0.001

  if (lneutraldrag) then
    ust = fkar*u1 / (log(zf(1)/z0))
  else
    if(myid==0) then
      do
        obl    = -ust**3*thl0av(1)/(fkar*grav*wthvsurf)
        !--- unstable case ---
        if (obl<-eps) then
          phimzf = (1.-16.*zf(1)/obl)**(-0.25)
          psimzf = 2.*log(0.5*(1+1/phimzf))+log(0.5*(1.+1/phimzf**2)) &
                      -2.*atan(1/phimzf)+pi*0.5
          phimz0 = (1.-16.*z0/obl)**(-0.25)
          psimz0 = 2.*log(0.5*(1+1/phimz0))+log(0.5*(1.+1/phimz0**2)) &
                      -2.*atan(1/phimz0)+pi*0.5


        !--- stable case ---
        elseif (obl>eps) then
    !         phimzf = (1.+5.*zf(1)/obl)
          psimzf = -5.*zf(1)/obl

        !--- neutral case ---
        else
          psimzf = 0.
        endif

        if (abs(psimzf) >= log(zf(1)/z0)) then
          psimzf = 0.
        end if

        ustold = ust
        ust = fkar*u1/(log(zf(1)/z0)-psimzf + psimz0)

        conv = abs(ust-ustold)
        if (conv<=eps) exit
      end do
    end if

    call MPI_BCAST(ust  ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(obl  ,1,MY_REAL   ,0,comm3d,mpierr)

  end if

  return
  end subroutine calust

  subroutine stopsurf(obl,ust,tst,qst,chku,chkt,chkq,z0m,z0h,z0q)

!**************************************************************
!     Extra subroutine used only in case of bad convergence   *
!**************************************************************

  use modglobal, only : i1,j1,rslabs,timee
  use modmpi,    only : my_real,mpierr,comm3d,mpi_sum,myid
  use modsurfdata,only : ustar,tstar,qstar
  implicit none

  real ustav ,  tstav ,  qstav ,  sigu  ,  sigt  ,  sigq
  real ustavl,  tstavl,  qstavl,  sigul ,  sigtl ,  sigql
  real ust, tst, qst
  real chku, chkt, chkq
  real z0q, z0h, z0m
  real obl
  integer i, j

  ustavl= 0.
  tstavl= 0.
  qstavl= 0.
  sigul = 0.
  sigtl = 0.
  sigql = 0.

  ustav = 0.
  tstav = 0.
  qstav = 0.
  sigu  = 0.
  sigt  = 0.
  sigq  = 0.

  do j=2,j1
  do i=2,i1
    ustavl= ustavl+ ustar(i,j)
    tstavl= tstavl+ tstar(i,j)
    qstavl= qstavl+ qstar(i,j)
  end do
  end do

! Mathieu MPI communication

  call MPI_ALLREDUCE(ustavl, ustav, 1,  MY_REAL, &
                         MPI_SUM, comm3d,mpierr)
  call MPI_ALLREDUCE(tstavl, tstav, 1,  MY_REAL, &
                         MPI_SUM, comm3d,mpierr)
  call MPI_ALLREDUCE(qstavl, qstav, 1,  MY_REAL, &
                         MPI_SUM, comm3d,mpierr)

  ustav = ustav/rslabs
  tstav = tstav/rslabs
  qstav = qstav/rslabs

  do j=2,j1
  do i=2,i1
    sigul= sigul+ (ustar(i,j)-ustav)**2
    sigtl= sigtl+ (tstar(i,j)-tstav)**2
    sigql= sigql+ (qstar(i,j)-qstav)**2
  end do
  end do


! Mathieu MPI communication

  call MPI_ALLREDUCE(sigul, sigu, 1,  MY_REAL, &
                         MPI_SUM, comm3d,mpierr)
  call MPI_ALLREDUCE(sigtl, sigt, 1,  MY_REAL, &
                         MPI_SUM, comm3d,mpierr)
  call MPI_ALLREDUCE(sigql, sigq, 1,  MY_REAL, &
                         MPI_SUM, comm3d,mpierr)

  sigu = sqrt(sigu/rslabs)
  sigt = sqrt(sigt/rslabs)
  sigq = sqrt(sigq/rslabs)


  if(myid==0)then

    write(6,'(/A,F6.0)') &
          'RESULTS OF SUBROUTINE SURFACE FOR TIME =',timee

    write(6,'(/A//A,/3E11.3)') &
          'Iteration procedure statistics', &
          '    chku       chkt       chkq  ', &
          chku,chkt,chkq

    write(6,'(/A,/4E15.7)') &
          '      ustar          tstar          qstar           L  ', &
          ust,tst,qst,obl

    write(6,'(/A,/3E15.7)') &
          '     z0m(m)         z0h(m)         z0q(m)  ', &
          z0m,z0h,z0q

    write(6,'(/A//A/A,3E15.7/A,3E15.7)') &
          'Slab average and sigma of local values', &
          '              ustar          tstar          qstar', &
          'average:',ustav,tstav,qstav, &
          'sigma:  ',sigu,sigt,sigq

  endif ! end if(myid==0)

  stop 'Bad convergence in iteration procedure'


  end subroutine stopsurf

  !----------------------------------------------------------------
  subroutine qtsurf

  use modglobal, only : tmelt,bt,at,rd,rv,cp,es0, pref0
  use modsurfdata, only : qts,thls,ps
  implicit none
  real exner, tsurf, es
  exner = (ps/pref0)**(rd/cp)
  tsurf = thls*exner
  es    = es0*exp(at*(tsurf-tmelt)/(tsurf-bt))
  qts   = rd/rv*es/(ps-(1-rd/rv)*es)

  return
  end subroutine qtsurf
!----------------------------------------------------------------

end module modsurf
