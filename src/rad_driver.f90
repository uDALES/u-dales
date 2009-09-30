!> \file rad_driver.f90
!!  Manages the full radiation McICA scheme

!>
!!  Manages the full radiation McICA scheme
!>
!!  \author Robert Pincus
!!  \author Bjorn Stevens
!!  \author Thijs Heus
!!  \todo Documentation
!!  \par Revision list
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
! Copyright 1999-2008, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module radiation

  use fuliou,      only : rad
  use modglobal,   only : cexpnr,cp,rcp,cpr,rhow,pref0,pi,xlat,xlon,xday,xtime,timee
  use rad_solver,  only : nv1, nv, SolarConstant
  use modraddata,  only : zenith,useMcICA
  implicit none


  logical, save     :: first_time = .True.
  real, allocatable, save ::  pp(:), pt(:), ph(:), po(:), pre(:), pde(:), &
       plwc(:), piwc(:), prwc(:), pgwc(:), fds(:), fus(:), fdir(:), fuir(:)

  integer :: i,j,k, npts
  real    :: ee, u0, day, time, alat, zz

  contains

    subroutine d4stream(i1,ih,j1,jh,k1, tskin, albedo, CCN, dn0, &
         pi0,  tk, rv, rc, fds3D,fus3D,fdir3D,fuir3D, rr)
      implicit none

      integer, intent (in) :: i1,ih,j1,jh,k1
      real, intent (in)    :: CCN
      ! CvH real, intent (in)    :: sknt, sfc_albedo, CCN
      real, dimension (k1), intent (in)                 :: dn0, pi0
      real, dimension (2-ih:i1+ih,2-jh:j1+jh,k1), intent (in)  :: tk, rv, rc
      real, optional, dimension (2-ih:i1+ih,2-jh:j1+jh,k1), intent (in) :: rr
      real, dimension (2-ih:i1+ih,2-jh:j1+jh,k1), intent (out) :: fus3D,fds3D,fuir3D,fdir3D
      real, dimension (i1+1,i1+1), intent (in) :: tskin, albedo

      integer :: kk
      real    :: prw, p0(k1), exner(k1), pres(k1)
      character (len=19) :: background

      if (first_time) then
         p0(k1) = (pref0*(pi0(k1)/cp)**cpr) / 100.
         p0(k1-1) = (pref0*(pi0(k1-1)/cp)**cpr) / 100.
         background  = 'backrad.inp.'//cexpnr
         call setup(background,k1,npts,nv1,nv,p0)
         first_time = .False.
         if (allocated(pre))   pre(:) = 0.
         if (allocated(pde))   pde(:) = 0.
         if (allocated(piwc)) piwc(:) = 0.
         if (allocated(prwc)) prwc(:) = 0.
         if (allocated(plwc)) plwc(:) = 0.
         if (allocated(pgwc)) pgwc(:) = 0.
      end if
      !
      ! initialize surface albedo, emissivity and skin temperature.
      !
      ee = 1.0
      !
      ! determine the solar geometery, as measured by u0, the cosine of the
      ! solar zenith angle
      !
      u0 = zenith(xtime + timee/3600,xday,xlat,xlon)
      !
      ! call the radiation
      !
      prw = (4./3.)*pi*rhow
      do j=2,j1
         do i=2,i1
            do k=1,k1
               exner(k)= (pi0(k))/cp
               pres(k) = pref0 * (exner(k))**cpr
            end do
            pp(nv1) = 0.5*(pres(1)+pres(2)) / 100.
            do k=2,k1
               kk = nv-(k-2)
               pt(kk) = tk(i,j,k)
               ph(kk) = rv(i,j,k)
               if (present(rr)) then
                  plwc(kk) = 1000.*dn0(k)*max(0.,(rc(i,j,k)-rr(i,j,k)))
                  prwc(kk) = 1000.*dn0(k)*rr(i,j,k)
               else
                  plwc(kk) = 1000.*dn0(k)*rc(i,j,k)
                  prwc(kk) = 0.
               end if
               pre(kk)  = 1.e6*(plwc(kk)/(1000.*prw*CCN*dn0(k)))**(1./3.)
               if (plwc(kk).le.0.) pre(kk) = 0.
               if (k < k1) pp(kk) = 0.5*(pres(k)+pres(k+1)) / 100.
            end do
            pp(nv-k1+2) = pres(k1)/100. - 0.5*(pres(k1-1)-pres(k1)) / 100.

            call rad( albedo(i,j), u0, SolarConstant, tskin(i,j), ee, pp, pt, ph, po,&
                 fds, fus, fdir, fuir, plwc=plwc, pre=pre, useMcICA=useMcICA)

            do k=1,k1
               kk = nv1 - (k-1)
               fus3d(i,j,k) =  fus(kk)
               fds3d(i,j,k) =  fds(kk)
               fuir3d(i,j,k) = fuir(kk)
               fdir3d(i,j,k) = fdir(kk)
!                rflx(i,j,k) = sflx(i,j,k) + fuir(kk) - fdir(kk)
            end do
            
!           CvH TO BE REMOVED
!           if (u0 > 0.) then
!              albedo(i,j) = fus(1)/fds(1)
!           else
!              albedo(i,j) = -999.
!           end if

!             do k=2,k1-3
!                xfact  = exner(k)*dzm(k)/(cp*dn0(k))
!                tp(i,j,k) = - (rflx(i,j,k) - rflx(k-1,i,j))*xfact
!             end do

         end do
      end do

    end subroutine d4stream

  !>
  !! sets up the input data to extend through an atmopshere of appreiciable
  !! depth using a background souding specified as a paramter, match this to
  !! the original sounding using p0 as this does not depend on time and thus
  !! allows us to recompute the same background matching after a history start
  !!
  subroutine setup(filenm,k1,npts,nv1,nv,zp)
  implicit none

    character (len=19), intent (in) :: filenm
    integer, intent (in) :: k1
    integer, intent (out):: npts,nv1,nv
    real, intent (in)    :: zp(k1)

    real, allocatable  :: sp(:), st(:), sh(:), so(:), sl(:)

    integer :: k, ns, norig, index
    logical :: blend
    real    :: pa, pb, ptop, ptest, test, dp1, dp2, dp3, Tsurf

    norig = 0
    open ( unit = 08, file = filenm, status = 'old' )
    print *, 'Reading Background Sounding: ',filenm
    read (08,*) Tsurf, ns
    allocate ( sp(ns), st(ns), sh(ns), so(ns), sl(ns))
    do k=1,ns
       read ( 08, *) sp(k), st(k), sh(k), so(k), sl(k)
    enddo
    close (08)

    !
    ! identify what part, if any, of background sounding to use
    !
    ptop = zp(k1)
    if (sp(2) < ptop) then
       pa = sp(1)
       pb = sp(2)
       k = 3
       do while (sp(k) < ptop)
          pa = pb
          pb = sp(k)
          k  = k+1
       end do
       k=k-1           ! identify first level above top of input
       blend = .True.
    else
       blend = .False.
    end if
    !
    ! if blend is true then the free atmosphere above the sounding will be
    ! specified based on the specified background climatology, here the
    ! pressure levels for this part of the sounding are determined
    !
    if (blend) then
       dp1 = pb-pa
       dp2 = ptop - pb
       dp3 = zp(k1-1) - zp(k1)
       if (dp1 > 2.*dp2) k = k-1 ! first level is too close, blend from prev
       npts  = k
       norig = k
       ptest = sp(k)
       test = ptop-ptest
       do while (test > 2*dp3)
          ptest = (ptest+ptop)*0.5
          test  = ptop-ptest
          npts  = npts + 1
       end do
       nv1 = npts + k1
    else
       nv1 = k1
    end if
    nv = nv1-1
    !
    ! allocate the arrays for the sounding data to be used in the radiation
    ! profile and then fill them first with the sounding data, by afill, then
    ! by interpolating the background profile at pressures less than the
    ! pressure at the top fo the sounding
    !
    allocate (pp(nv1),fds(nv1),fus(nv1),fdir(nv1),fuir(nv1))
    allocate (pt(nv),ph(nv),po(nv),pre(nv),pde(nv),plwc(nv),prwc(nv))

    if (blend) then
       pp(1:norig) = sp(1:norig)
       pt(1:norig) = st(1:norig)
       ph(1:norig) = sh(1:norig)
       po(1:norig) = so(1:norig)

       do k=norig+1,npts
          pp(k) = (ptop + pp(k-1))*0.5
          index = getindex(sp,ns,pp(k))
          pt(k) =  intrpl(sp(index),st(index),sp(index+1),st(index+1),pp(k))
          ph(k) =  intrpl(sp(index),sh(index),sp(index+1),sh(index+1),pp(k))
          po(k) =  intrpl(sp(index),so(index),sp(index+1),so(index+1),pp(k))
       end do
       !
       ! set the ozone constant below the reference profile
       !
       do k=npts+1,nv
          po(k) =  po(npts)
       end do
    end if

  end subroutine setup
  !>
  !!  locate the index closest to a value
  !!
  integer function getindex(x,n,xval)
  implicit none

    integer, intent (in)  :: n
    real,    intent (in)  :: x(n),xval

    integer :: ia, ib

    ia=1
    ib=n
    if (xval < x(1)) then
       getindex = 1
    elseif (xval > x(n)) then
       getindex = n-1
    else
       getindex = (ia+ib)/2
       do while (getindex /= ia .or. ib /= getindex+1)
          getindex = (ia+ib)/2
          if ((xval-x(getindex)) >= 0.0) then
             ia = getindex
          else
             ib = getindex
          end if
       end do
    endif

  end function getindex

  !> linear interpolation between two points,
  !
  real function intrpl(x1,y1,x2,y2,x)
  implicit none

    real, intent (in)  :: x1,y1,x2,y2,x

    real :: slope

    slope  = (y2-y1)/(x2 - x1 + epsilon(1.))
    intrpl = y1+slope*(x-x1)

  end function intrpl

!   !> Return the cosine of the solar zenith angle give the decimal day and
!   !> the latitude
!   !
!   real function zenith()
!     use modglobal, only : xtime, timee, pi, xday,xlat,xlon
!     implicit none
!     real :: rtime,phi,el,obliq,xlam,declin,hora
!
!     rtime   = xtime + timee/3600
!     phi    = xlat * pi/180.
!     el     = xlon * pi/180.
!     obliq  = 23.45 * pi/180.
!     xlam   = 4.88 + 0.0172 * xday
!     declin = asin(sin(obliq)*sin(xlam))
!     hora   = el-pi + 2.*pi*(rtime/24.)
!     zenith = max(0.,sin(declin)*sin(phi)+cos(declin)*cos(phi)* &
!                                                          cos(hora))
!   end function zenith
end module radiation
