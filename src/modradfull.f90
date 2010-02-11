!> \file modradfull.f90
!!  Full radiative transfer

!>
!!  Calculates the radiative sources
!>
!!  \author Bjorn Stevens
!!  \author Robert Pincus
!!  \author Thijs Heus, MPI-M
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

module modradfull

  use RandomNumbers

  implicit none
  private
  public :: radfull

  logical, save     :: d4stream_initialized = .False.
  real, allocatable, save ::  pp(:), pt(:), ph(:), po(:), pre(:), pde(:), &
       plwc(:), piwc(:), prwc(:), pgwc(:), fds(:), fus(:), fdir(:), fuir(:)

  integer :: i,j,k, npts
  real    :: ee, u0, day, time, alat, zz


 !
  ! u denotes the Double-Gauss quadratures and weights (Sykes, 1951).
  ! the p0-p3 denote the legendre polynomials for x=u(1) to u(4)
  ! respectively.  p11d(4,4), p22d(4,4), and p33d(4,4) are defined as
  ! 0.5*p1d(i)*p1d(j), 0.5*p2d(i)*p2d(j), and 0.5*p3d(i)*p3d(j) resp.
  !

  integer :: nv,nv1,mb
  real :: totalpower
  real,parameter :: SolarConstant      = 1.365d+3  !< The Solar radiation constant

  real, parameter :: u(4) = (/-0.7886752,-0.2113247,0.2113247,0.7886752/)
  real, parameter :: p0d(4)= (/ 1., 1., 1., 1./)
  real, parameter :: p1d(4)= (/-0.788675, -.211325,  .211325, .788675/)
  real, parameter :: p2d(4)= (/ 0.433013, -.433013, -.433013, .433013/)
  real, parameter :: p3d(4)= (/-.0433940,  .293394, -.293394, .043394/)
  real, parameter, dimension(4, 4) ::                       &
    p11d = reshape(                                         &
    (/ .311004E+00, .833334E-01,-.833334E-01,-.311004E+00,  &
       0.833334E-01, .223291E-01,-.223291E-01,-.833334E-01, &
       -.833334E-01,-.223291E-01, .223291E-01, .833334E-01, &
       -.311004E+00,-.833334E-01, .833334E-01, .311004E+00  /), (/ 4, 4 /))
  real, parameter, dimension(4, 4) ::                       &
   p22d = reshape(                                          &
    (/ .937501E-01,-.937501E-01,-.937501E-01, .937501E-01,  &
       -.937501E-01, .937501E-01, .937501E-01,-.937501E-01, &
       -.937501E-01, .937501E-01, .937501E-01,-.937501E-01, &
       0.937501E-01,-.937501E-01,-.937501E-01, .937501E-01  /), (/ 4, 4 /))
  real, parameter, dimension(4, 4) ::                       &
    p33d = reshape(                                         &
    (/ .941520E-03,-.636577E-02, .636577E-02,-.941520E-03,  &
       -.636577E-02, .430400E-01,-.430400E-01, .636577E-02, &
       0.636577E-02,-.430400E-01, .430400E-01,-.636577E-02, &
       -.941520E-03, .636577E-02,-.636577E-02, .941520E-03  /), (/ 4, 4 /))

  logical, save :: fuliou_Initialized = .False.
  type(randomNumberSequence), save :: randoms
  real, parameter :: minSolarZenithCosForVis = 1.e-4

    logical, save      :: ckd_Initialized = .False.
  real, save         :: bllmx, brlmn
  integer, save      :: ngases

  TYPE ckd_properties
     character (len=5) :: name
     integer           :: ng, nt, np, noverlap, iband
     real              :: mweight, default_conc, tbase
     real, allocatable :: hk(:),sp(:),xk(:,:,:,:)
  END TYPE ckd_properties

  TYPE band_properties
     private
     integer              :: kg, ngases
     real                 :: llimit, rlimit, center, power = 0.
     integer, allocatable :: gas_id(:)
     real, allocatable    :: hk(:)
  END TYPE band_properties

  TYPE (ckd_properties),  allocatable :: gas(:)
  TYPE (band_properties), allocatable :: band(:), solar_bands(:), ir_bands(:)

  integer, save :: nsizes = 8
  logical, save :: cldwtr_Initialized = .False.

  real, allocatable    :: re(:), fl(:), bz(:,:), wz(:,:), gz(:,:)

contains
    subroutine radfull
  !   use radiation,    only : d4stream
    use modglobal,    only : imax,i1,ih,jmax,j1,jh,kmax,k1,cp,dzf,dzh,rlv,rd,zf,pref0
    use modfields,    only : rhof, exnf,exnh, thl0,qt0,ql0,sv0
    use modsurfdata,  only : albedo, tskin, qskin, thvs, qts, ps
    use modmicrodata, only : imicro, imicro_bulk, Nc_0,iqr
    use modraddata,   only : thlprad, lwd,lwu,swd,swu,rho_air_mn
use modmpi, only : myid
      implicit none
    real :: thlpld,thlplu,thlpsd,thlpsu
    real, dimension(k1)  :: rhof_b, exnf_b
    real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: temp_b, qv_b, ql_b,rr_b
    real, dimension(1:i1+1,1:j1+1) :: tempskin
    integer :: i,j,k

    real :: exnersurf

!take care of UCLALES z-shift for thermo variables.
      do k=1,kmax
        rhof_b(k+1)     = rhof(k)
        exnf_b(k+1)     = exnf(k)
        do j=2,j1
          do i=2,i1
            qv_b(i,j,k+1)   = qt0(i,j,k) - ql0(i,j,k)
            ql_b(i,j,k+1)   = ql0(i,j,k)
            temp_b(i,j,k+1) = thl0(i,j,k)*exnf(k)+(rlv/cp)*ql0(i,j,k)
            if (imicro==imicro_bulk) rr_b(i,j,k+1) = sv0(i,j,k,iqr)
          end do
        end do
      end do

      !take care of the surface boundary conditions
      !CvH edit, extrapolation creates instability in surface scheme
      exnersurf = (ps/pref0) ** (rd/cp)
      rhof_b(1) = rhof(1) + dzh(1)/dzh(2)*(rhof(1)-rhof(2))
      exnf_b(1) = exnh(1) + dzh(1)/dzh(2)*(rhof(1)-rhof(2))
    

      do j=2,j1
        do i=2,i1
          ql_b(i,j,1)   = 0
          qv_b(i,j,1)   = qv_b(i,j,2) +dzh(1)/dzh(2)*(qv_b(i,j,2)-qv_b(i,j,3))
          temp_b(i,j,1) = temp_b(i,j,2) +dzh(1)/dzh(2)*(temp_b(i,j,2)-temp_b(i,j,3))
        end do
      end do
      tempskin = 0.5*(temp_b(1:i1+1,1:j1+1,1)+temp_b(1:i1+1,1:j1+1,2))
      ! tempskin = tskin*exnh(1)
      !CvH end edit

      if (imicro==imicro_bulk) then
        rr_b(:,:,1) = 0.
        call d4stream(i1,ih,j1,jh,k1,tempskin,albedo,Nc_0,rhof_b,exnf_b*cp,temp_b,qv_b,ql_b,swd,swu,lwd,lwu,rr=rr_b)
      else
        call d4stream(i1,ih,j1,jh,k1,tempskin,albedo,Nc_0,rhof_b,exnf_b*cp,temp_b,qv_b,ql_b,swd,swu,lwd,lwu)
      end if

!Downward radiation fluxes are pointing downward in UCLALES, pointing upward in DALES
      lwd = -lwd
      swd = -swd
!      lwd(:,:,1) = lwd(:,:,1)+0.3333333*(lwd(:,:,1)-lwd(:,:,2))
!      swd(:,:,1) = swd(:,:,1)+0.3333333*(swd(:,:,1)-swd(:,:,2))
!      lwu(:,:,1) = lwu(:,:,1)+0.3333333*(lwu(:,:,1)-lwu(:,:,2))
!      swu(:,:,1) = swu(:,:,1)+0.3333333*(swu(:,:,1)-swu(:,:,2))

!Add up thl tendency
      do k=1,kmax
        do j=2,j1
          do i=2,i1
            thlpld          = -(lwd(i,j,k+1)-lwd(i,j,k))/(rhof(k)*cp*dzf(k))
            thlplu          = -(lwu(i,j,k+1)-lwu(i,j,k))/(rhof(k)*cp*dzf(k))
            thlpsd          = -(swd(i,j,k+1)-swd(i,j,k))/(rhof(k)*cp*dzf(k))
            thlpsu          = -(swu(i,j,k+1)-swu(i,j,k))/(rhof(k)*cp*dzf(k))

            thlprad(i,j,k)  = thlprad(i,j,k) + thlpld+thlplu+thlpsu+thlpsd
          end do
        end do
      end do

    end subroutine radfull

    subroutine d4stream(i1,ih,j1,jh,k1, tskin, albedo, CCN, dn0, &
         pi0,  tk, rv, rc, fds3D,fus3D,fdir3D,fuir3D, rr)
      use modglobal, only : cexpnr,cp,cpr,pi,pref0,timee,xday,xlat,xlon,xtime,rhow
      use modraddata,only : useMcICA,zenith
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

      if (.not. d4stream_initialized) then
         p0(k1) = (pref0*(pi0(k1)/cp)**cpr) / 100.
         p0(k1-1) = (pref0*(pi0(k1-1)/cp)**cpr) / 100.
         background  = 'backrad.inp.'//cexpnr
         call d4stream_setup(background,k1,npts,nv1,nv,p0)
         d4stream_initialized = .True.
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
            end do

         end do
      end do

    end subroutine d4stream

  !>
  !! sets up the input data to extend through an atmopshere of appreiciable
  !! depth using a background souding specified as a paramter, match this to
  !! the original sounding using p0 as this does not depend on time and thus
  !! allows us to recompute the same background matching after a history start
  !!
  subroutine d4stream_setup(filenm,k1,npts,nv1,nv,zp)
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

  end subroutine d4stream_setup
  !> coefficient calculations for four first-order differential equations.
  !>
  !> See the paper by Liou, Fu and Ackerman (1988) for the formulation of
  !> the delta-four-stream approximation in a homogeneous layer.
  subroutine coefft(solar,w,w1,w2,w3,t0,t1,u0,f0,aa,zz,a1,z1,fk1,fk2)

    logical, intent (in) :: solar
    real, intent(in)     :: w, w1, w2, w3, t0, t1, u0, f0
    real, intent (out)   :: aa(4,4,2), zz(4,2), a1(4,4), z1(4), fk1, fk2

    integer :: i, j
    real :: a(2,2,2), b(4,3), c(4,5), d(4), z(4), b1, c1, fw, fw1, fw2 &
         ,fw3, fw4, w0w, w1w, w2w, w3w, q1, q2, q3, zx, x, y, a2, b2   &
         ,fq0, fq1, fq, dt

    x = 0.5 * w
    w0w = x
    w1w = x * w1
    w2w = x * w2
    w3w = x * w3
    if ( solar ) then
       fw = u0 * u0
       q1 = - w1w * u0
       q2 =   w2w * ( 1.5 * fw - 0.5 )
       q3 = - w3w * ( 2.5 * fw - 1.5 ) * u0
       do i = 1, 4
          c(i,5) = (w0w + q1*p1d(i) + q2*p2d(i) + q3*p3d(i))/u(i)
       end do
    else
       do i = 1, 4
          c(i,5) = 1.0/ u(i)
       end do
    endif

    fq = 0.5 * w0w
    do i = 3, 4
       do j = 1, 4
          c(i,j) = fq + w1w*p11d(i,j) + w2w*p22d(i,j) + w3w*p33d(i,j)
          if ( i .eq. j ) then
             c(i,j) = ( c(i,j) - 1.0 ) / u(i)
          else
             c(i,j) = c(i,j) / u(i)
          endif
       end do
    end do

    b(1,1) = c(4,4) - c(4,1)
    b(1,2) = c(4,4) + c(4,1)
    b(2,1) = c(4,3) - c(4,2)
    b(2,2) = c(4,3) + c(4,2)
    b(3,1) = c(3,4) - c(3,1)
    b(3,2) = c(3,4) + c(3,1)
    b(4,1) = c(3,3) - c(3,2)
    b(4,2) = c(3,3) + c(3,2)
    b(1,3) = c(4,5) - c(1,5)
    b(2,3) = c(3,5) - c(2,5)
    b(3,3) = c(3,5) + c(2,5)
    b(4,3) = c(4,5) + c(1,5)
    ! **********************************************************************
    ! coefficient calculations for second order differential equations.
    ! **********************************************************************
    fw1 = b(1,1) * b(1,2)
    fw2 = b(2,1) * b(3,2)
    fw3 = b(3,1) * b(2,2)
    fw4 = b(4,1) * b(4,2)
    a(2,2,1) = fw1 + fw2
    a(2,1,1) = b(1,1) * b(2,2) + b(2,1) * b(4,2)
    a(1,2,1) = b(3,1) * b(1,2) + b(4,1) * b(3,2)
    a(1,1,1) = fw3 + fw4
    a(2,2,2) = fw1 + fw3
    a(2,1,2) = b(1,2) * b(2,1) + b(2,2) * b(4,1)
    a(1,2,2) = b(3,2) * b(1,1) + b(4,2) * b(3,1)
    a(1,1,2) = fw2 + fw4
    d(1) = b(3,2) * b(4,3) + b(4,2) * b(3,3) + b(2,3) / u0
    d(2) = b(1,2) * b(4,3) + b(2,2) * b(3,3) + b(1,3) / u0
    d(3) = b(3,1) * b(1,3) + b(4,1) * b(2,3) + b(3,3) / u0
    d(4) = b(1,1) * b(1,3) + b(2,1) * b(2,3) + b(4,3) / u0

    ! **********************************************************************
    ! coefficient calculations for fourth-order differential equations.
    ! **********************************************************************
    x = u0 * u0
    b1 = a(2,2,1) + a(1,1,1)
    c1 = a(2,1,1) * a(1,2,1) - a(1,1,1) * a(2,2,1)
    z(1) = a(2,1,1) * d(3) + d(4) / x - a(1,1,1) * d(4)
    z(2) = a(1,2,1) * d(4) - a(2,2,1) *d(3) + d(3) / x
    z(3) = a(2,1,2) * d(1) + d(2) / x - a(1,1,2) * d(2)
    z(4) = a(1,2,2) * d(2) - a(2,2,2) * d(1) + d(1) / x

    ! **********************************************************************
    ! fk1 and fk2 are the eigenvalues.
    ! **********************************************************************
    dt = t1 - t0
    x = sqrt ( b1 * b1 + 4.0 * c1 )
    fk1 = sqrt ( ( b1 + x ) * 0.5 )
    fk2 = sqrt ( ( b1 - x ) * 0.5 )
    fw = u0 * u0
    x = 1.0 / ( fw * fw ) - b1 / fw - c1
    fw = 0.5 * f0 / x
    z(1) = fw * z(1)
    z(2) = fw * z(2)
    z(3) = fw * z(3)
    z(4) = fw * z(4)
    z1(1) = 0.5 * ( z(1) + z(3) )
    z1(2) = 0.5 * ( z(2) + z(4) )
    z1(3) = 0.5 * ( z(2) - z(4) )
    z1(4) = 0.5 * ( z(1) - z(3) )
    a2 = ( fk1 * fk1 - a(2,2,1) ) / a(2,1,1)
    b2 = ( fk2 * fk2 - a(2,2,1) ) / a(2,1,1)
    x = b(1,1) * b(4,1) - b(3,1) * b(2,1)
    fw1 = fk1 / x
    fw2 = fk2 / x
    y = fw2 * ( b2 * b(2,1) - b(4,1) )
    zx = fw1 * ( a2 * b(2,1) - b(4,1) )
    a1(1,1) = 0.5 * ( 1 - y )
    a1(1,2) = 0.5 * ( 1 - zx )
    a1(1,3) = 0.5 * ( 1 + zx )
    a1(1,4) = 0.5 * ( 1 + y )
    y = fw2 * ( b(3,1) - b2 * b(1,1) )
    zx = fw1 * ( b(3,1) - a2 * b(1,1) )
    a1(2,1) = 0.5 * ( b2 - y )
    a1(2,2) = 0.5 * ( a2 - zx )
    a1(2,3) = 0.5 * ( a2 + zx )
    a1(2,4) = 0.5 * ( b2 + y )
    a1(3,1) = a1(2,4)
    a1(3,2) = a1(2,3)
    a1(3,3) = a1(2,2)
    a1(3,4) = a1(2,1)
    a1(4,1) = a1(1,4)
    a1(4,2) = a1(1,3)
    a1(4,3) = a1(1,2)
    a1(4,4) = a1(1,1)
    if ( solar ) then
       fq0 = exp ( - t0 / u0 )
       fq1 = exp ( - t1 / u0 )
    else
       fq0 = 1.0
       fq1 = exp ( - dt / u0 )
    endif
    x = exp ( - fk1 * dt )
    y = exp ( - fk2 * dt )
    do i = 1, 4
       zz(i,1) = z1(i) * fq0
       zz(i,2) = z1(i) * fq1
       aa(i,1,1) = a1(i,1)
       aa(i,2,1) = a1(i,2)
       aa(i,3,1) = a1(i,3) * x
       aa(i,4,1) = a1(i,4) * y
       aa(i,3,2) = a1(i,3)
       aa(i,4,2) = a1(i,4)
       aa(i,1,2) = a1(i,1) * y
       aa(i,2,2) = a1(i,2) * x
    end do

  end subroutine coefft

  !> In the limits of no scattering ( Fu, 1991 ), fk1 = 1.0 / u(3) and
  !> fk2 = 1.0 / u(4).
  subroutine coefft0( solar,t0,t1,u0,f0,aa,zz,a1,z1,fk1,fk2)

    logical, intent (in) :: solar
    real, intent(in)     :: t0, t1, u0, f0
    real, intent (out)   :: aa(4,4,2), zz(4,2), a1(4,4), z1(4), fk1, fk2

    integer :: i, j, k, jj
    real    :: x, y, dt, fw

    fk1 = 4.7320545
    fk2 = 1.2679491
    y = exp ( - ( t1 - t0 ) / u0 )
    fw = 0.5 * f0
    do i = 1, 4
       if ( solar ) then
          z1(i) = 0.0
          zz(i,1) = 0.0
          zz(i,2) = 0.0
       else
          jj = 5 - i
          z1(i) = fw / ( 1.0 + u(jj) / u0 )
          zz(i,1) = z1(i)
          zz(i,2) = z1(i) * y
       endif
       do j = 1, 4
          a1(i,j) = 0.0
          do k = 1, 2
             aa(i,j,k) = 0.0
          end do
       end do
    end do
    do  i = 1, 4
       j = 5 - i
       a1(i,j) = 1.0
    end do
    dt = t1 - t0
    x = exp ( - fk1 * dt )
    y = exp ( - fk2 * dt )
    aa(1,4,1) = y
    aa(2,3,1) = x
    aa(3,2,1) = 1.0
    aa(4,1,1) = 1.0
    aa(1,4,2) = 1.0
    aa(2,3,2) = 1.0
    aa(3,2,2) = x
    aa(4,1,2) = y

  end subroutine coefft0

  !> **********************************************************************
  !> In the solar band  asbs is the surface albedo, while in the infrared
  !> band asbs is  blackbody intensity emitted at the surface temperature
  !> times surface emissivity.  In this subroutine, the delta-four-stream
  !> is applied to nonhomogeneous atmospheres. See comments in subroutine
  !> 'qcfel' for array AB(13,4*n).
  !> **********************************************************************
  subroutine qccfe (solar,asbs,ee,t,w,w1,w2,w3,u0,f0,fk1,fk2,a4,g4,z4 )

    logical, intent (in)              :: solar
    real, intent (in)                 :: asbs, ee
    real, dimension (nv), intent (in) :: t, w, w1, w2, w3, u0, f0
    real, dimension (nv), intent (out):: fk1, fk2
    real, intent(out)                 :: a4(4,4,nv), g4(4,nv), z4(4,nv)

    integer :: i, j, k, m18, m28, n4, kf, j1, j2, j3, i1
    integer :: i2, i3, i8, m2, m1
    real    :: fu(4,4), wu(4), ab(13,4*nv), bx(4*nv), xx(4*nv)
    real    :: aa(4,4,2), zz(4,2), a1(4,4), z1(4), fw1, fw2, wn, v1, v2, v3

    n4 = 4*nv
    do i = 1, n4
       do j = 1, 13
          ab(j,i) = 0.0
       end do
    end do

    wn = w(1)
    if ( wn <= 1.0e-4 ) then
       call coefft0(solar,0.0,t(1),u0(1),f0(1),aa,zz,a1,z1,fk1(1),fk2(1))
    else
       if ( wn >= 0.999999 ) wn = 0.999999
       call coefft(solar,wn,w1(1),w2(1),w3(1),0.0,t(1),u0(1),f0(1), &
            aa,zz,a1,z1,fk1(1),fk2(1))
    endif
    do  i = 1, 4
       z4(i,1) = z1(i)
       do  j = 1, 4
          a4(i,j,1) = a1(i,j)
       end do
    end do
    do  i = 1, 2
       bx(i) = - zz(i+2,1)
       i8 = i + 8
       do j = 1, 4
          ab(i8-j,j) = aa(i+2,j,1)
       end do
    end do
    do i = 1, 4
       wu(i) = zz(i,2)
       do j = 1, 4
          fu(i,j) = aa(i,j,2)
       end do
    end do
    do k = 2, nv
       wn = w(k)
       if ( w(k) .le. 1.0e-4 ) then
          call coefft0(solar,t(k-1),t(k),u0(k),f0(k),aa,zz,a1,z1,fk1(k),fk2(k))
       else
          if ( wn .ge. 0.999999 ) wn = 0.999999
          call coefft(solar,wn,w1(k),w2(k),w3(k),t(k-1),t(k),u0(k),f0(k),  &
               aa,zz,a1,z1,fk1(k),fk2(k))
       endif
       do  i = 1, 4
          z4(i,k) = z1(i)
          do  j = 1, 4
             a4(i,j,k) = a1(i,j)
          end do
       end do
       kf = k + k + k + k
       i1 = kf - 5
       i2 = i1 + 3
       j1 = kf - 7
       j2 = j1 + 3
       i3 = 0
       do  i = i1, i2
          i3 = i3 + 1
          bx(i) = - wu(i3) + zz(i3,1)
          j3 = 0
          i8 = i + 8
          do  j = j1, j2
             j3 = j3 + 1
             ab(i8-j,j) = fu(i3,j3)
          end do
          j3 = 0
          do j = j2 + 1, j2 + 4
             j3 = j3 + 1
             ab(i8-j,j) = - aa(i3,j3,1)
          end do
       end do
       do  i = 1, 4
          wu(i) = zz(i,2)
          do j = 1, 4
             fu(i,j) = aa(i,j,2)
          end do
       end do
    end do
    if ( solar ) then
       v1 = 0.2113247 * asbs
       v2 = 0.7886753 * asbs
       v3 = asbs * u0(1) * f0(1) * exp ( - t(nv) / u0(1) )
    else
       v1 = 0.2113247 * ( 1.0 - ee )
       v2 = 0.7886753 * ( 1.0 - ee )
       v3 = asbs
    end if
    m1 = n4 - 1
    m2 = n4
    m18 = m1 + 8
    m28 = m2 + 8
    fw1 = v1 * wu(3)
    fw2 = v2 * wu(4)
    bx(m1) = - ( wu(1) - fw1 - fw2 - v3 )
    bx(m2) = - ( wu(2) - fw1 - fw2 - v3 )
    do j = 1, 4
       j1 = n4 - 4 + j
       fw1 = v1 * fu(3,j)
       fw2 = v2 * fu(4,j)
       ab(m18-j1,j1) = fu(1,j) - fw1 - fw2
       ab(m28-j1,j1) = fu(2,j) - fw1 - fw2
    end do
    call qcfel (ab, bx, xx)
    do k = 1, nv
       j = k + k + k + k - 4
       do i = 1, 4
          j = j + 1
          g4(i,k) = xx(j)
       end do
    end do

  end subroutine qccfe

  !> **********************************************************************
  !> 1. `qcfel' is the abbreviation of ` qiu constants for each layer'.
  !> 2. The inhomogeneous atmosphere is divided into n adjacent homogeneous
  !>    layers where the  single scattering properties are constant in each
  !>    layer and allowed to vary from one to another. Delta-four-stream is
  !>    employed for each homogeneous layer. The boundary conditions at the
  !>    top and bottom of the atmosphere,  together with  continuity condi-
  !>    tions  at  layer interfaces lead to a system of algebraic equations
  !>    from which 4*n unknown constants in the problom can be solved.
  !> 3. This subroutine is used for solving the 4*n unknowns of A *X = B by
  !>    considering the fact that the coefficient matrix is a sparse matrix
  !>    with the precise pattern in this special problom.
  !> 4. The method is not different in principle from the general scheme of
  !>    Gaussian elimination with backsubstitution, but carefully optimized
  !>    so as to minimize arithmetic operations.  Partial  pivoting is used
  !>    to quarantee  method's numerical stability,  which will  not change
  !>    the basic pattern of sparsity of the matrix.
  !> 5. Scaling special problems so as to make  its nonzero matrix elements
  !>    have comparable magnitudes, which will ameliorate the stability.
  !> 6. a, b and x present A, B and X in A*X=B, respectively. and n4=4*n.
  !> 7. AB(13,4*n) is the matrix A in band storage, in rows 3 to 13; rows 1
  !>    and 2 and other unset elements should be set to zero on entry.
  !> 8. The jth column of A is stored in the jth column of the array AB  as
  !>    follows:
  !>            AB(8+i-j,j) = A(i,j) for max(1,j-5) <= i <= min(4*n,j+5).
  !>    Reversedly, we have
  !>            A(ii+jj-8,jj) = AB(ii,jj).
  !> **********************************************************************
  subroutine qcfel(ab, b, x)

    real, intent (inout)  :: ab(13,4*nv),b(4*nv)
    real, intent (out)    :: x(4*nv)

    integer :: i, j, k, l, m, n1, n2, n3, n4, m1, m2, m3, m4
    integer :: k44, n44, m18, m28, m38, m48, m1f, im1, i0m1, i0, i0f, ifq
    real    :: xx, yy, t, p

    i0=0
    n4 = 4*nv
    do  k = 1, nv - 1
       k44 = 4 * k - 4
       do  l= 1, 4
          m1 = k44 + l
          p = 0.0
          do  i = 8, 14 - l
             if ( abs ( ab(i,m1) ) .gt. abs ( p ) ) then
                p = ab(i,m1)
                i0 = i
             endif
          end do
          i0m1 = i0 + m1
          m18 = m1 + 8
          if ( i0 /= 8 ) then
             do  j = m1, m1 + 8 - l
                i0f = i0m1 - j
                m1f = m18 - j
                t = ab(i0f,j)
                ab(i0f,j) = ab(m1f,j)
                ab(m1f,j) = t
             end do
             i0f = i0m1 - 8
             t = b(i0f)
             b(i0f) = b(m1)
             b(m1) = t
          end if
          yy = ab(8,m1)
          ab(8,m1) = 1.0
          do  j = m1 + 1, m1 + 8 - l
             m1f = m18 - j
             ab(m1f,j) = ab(m1f,j) / yy
          end do
          b(m1) = b(m1) / yy
          do i = 9, 14 - l
             xx = ab(i,m1)
             ab(i,m1) = 0.0
             im1 = i + m1
             do  j = m1 + 1, m1 + 8 - l
                ifq = im1 - j
                m1f = m18 - j
                ab(ifq,j) = ab(ifq,j) - ab(m1f,j) * xx
             end do
             ifq = im1 - 8
             b(ifq) = b(ifq) - b(m1) * xx
          end do
       end do
    end do
    n44 = n4 - 4
    do l = 1, 3
       m1 = n44 + l
       p = 0.0
       do i = 8, 12 - l
          if ( abs ( ab(i,m1) ) .gt. abs ( p ) ) then
             p = ab(i,m1)
             i0 = i
          endif
       end do
       i0m1 = i0 + m1
       m18 = m1 + 8
       if( i0 /= 8 ) then
          do  j = m1, m1 + 4 - l
             i0f = i0m1 - j
             m1f = m18 - j
             t = ab(i0f,j)
             ab(i0f,j) = ab(m1f,j)
             ab(m1f,j) = t
          end do
          i0f = i0m1 - 8
          t = b(i0f)
          b(i0f) = b(m1)
          b(m1) = t
       end if
       yy = ab(8,m1)
       ab(8,m1) = 1.0
       do  j = m1 + 1, m1 + 4 - l
          m1f = m18 - j
          ab(m1f,j) = ab(m1f,j) / yy
       end do
       b(m1) = b(m1) / yy
       do  i = 9, 12 - l
          xx = ab(i,m1)
          ab(i,m1) = 0.0
          im1 = i + m1
          do  j = m1 + 1, m1 + 4 - l
             ifq = im1 - j
             m1f = m18 - j
             ab(ifq,j) = ab(ifq,j) - ab(m1f,j) * xx
          end do
          ifq = im1 - 8
          b(ifq) = b(ifq) - b(m1) * xx
       end do
    end do
    yy = ab(8,n4)
    ab(8,n4) = 1.0
    b(n4) = b(n4) / yy
    n3 = n4 - 1
    n2 = n3 - 1
    n1 = n2 - 1
    x(n4) = b(n4)
    x(n3) = b(n3) - ab(7,n4) * x(n4)
    x(n2) = b(n2) - ab(7,n3) * x(n3) - ab(6,n4) * x(n4)
    x(n1) = b(n1) - ab(7,n2) * x(n2) - ab(6,n3) * x(n3) - ab(5,n4) * x(n4)
    do  k = 1, nv - 1
       m4 = 4 * ( nv - k )
       m3 = m4 - 1
       m2 = m3 - 1
       m1 = m2 - 1
       m48 = m4 + 8
       m38 = m3 + 8
       m28 = m2 + 8
       m18 = m1 + 8
       x(m4) = b(m4)
       do  m = m4 + 1, m4 + 4
          x(m4) = x(m4) - ab(m48-m,m) * x(m)
       end do
       x(m3) = b(m3)
       do  m = m3 + 1, m3 + 5
          x(m3) = x(m3) - ab(m38-m,m) * x(m)
       end do
       x(m2) = b(m2)
       do  m = m2 + 1, m2 + 6
          x(m2) = x(m2) - ab(m28-m,m) * x(m)
       end do
       x(m1) = b(m1)
       do  m = m1 + 1, m1 + 7
          x(m1) = x(m1) - ab(m18-m,m) * x(m)
       end do
    end do

  end subroutine qcfel

  !> **********************************************************************
  !> In this subroutine, we incorporate a delta-function adjustment to
  !> account for the  forward  diffraction  peak in the context of the
  !> four-stream approximation ( Liou, Fu and Ackerman, 1988 ). w1(n),
  !> w2(n), w3(n), w(n), and t(n) are the adjusted parameters.
  !> **********************************************************************
  subroutine adjust (tt,ww,ww1,ww2,ww3,ww4,t,w,w1,w2,w3)

    real, dimension (nv), intent (in) :: tt,ww,ww1,ww2,ww3,ww4
    real, dimension (nv), intent (out):: t,w,w1,w2,w3

    integer :: k
    real    :: tt0, f, fw, dt(nv)

    tt0 = 0.0
    do  k = 1, nv
       f = ww4(k) / 9.0
       fw = 1.0 - f * ww(k)
       w1(k) = ( ww1(k) - 3.0 * f ) / ( 1.0 - f )
       w2(k) = ( ww2(k) - 5.0 * f ) / ( 1.0 - f )
       w3(k) = ( ww3(k) - 7.0 * f ) / ( 1.0 - f )
       w(k) = ( 1.0 - f ) * ww(k) / fw
       dt(k) = (tt(k) - tt0) * fw
       tt0 = tt(k)
    end do
    t(1) = dt(1)
    do k = 2, nv
       t(k) = dt(k) + t(k-1)
    end do

  end subroutine adjust
  !> --------------------------------------------------------------------------
  !> Subroutine qft: Delta 4-stream solver for fluxes
  !>
  subroutine qft (solar, ee, as, u0, bf, tt, ww, ww1, ww2, ww3, ww4, ffu, ffd)
    use modglobal, only : pi
    logical, intent (in) :: solar
    real, intent (in)    :: ee, as, u0
    real, dimension (nv), intent (in)   :: tt,ww,ww1,ww2,ww3,ww4
    real, dimension (nv1), intent (in)  :: bf
    real, dimension (nv1), intent (out) :: ffu, ffd

    real, dimension (nv) :: t,w,w1,w2,w3,u0a,f0a,fk1,fk2
    integer :: k, kk, ii, jj
    real    :: x(4), fi(4), a4(4,4,nv), z4(4,nv), g4(4,nv)
    real    :: tkm1, fw3, fw4, y1, xy, xas, xee
    real, parameter :: fw1 = 0.6638960, fw2 = 2.4776962

    call adjust(tt,ww,ww1,ww2,ww3,ww4,t,w,w1,w2,w3)

    if (solar) then
       fw3 = u0
       xee = 0.0
       xas = as
       do k = 1, nv
          u0a(k) = u0
          f0a(k) = 1./pi
       end do
    else
       fw3 = 0.
       xas = bf(nv1) * ee
       xee = ee
       tkm1 = 0.0
       do k = 1, nv
          f0a(k) = 2.0 * ( 1.0 - w(k) ) * bf(k)
          u0a(k) = -(t(k)-tkm1) / ( alog( bf(k+1)/bf(k) ) + epsilon(1.))
          tkm1 = t(k)
       end do
    end if
    call qccfe (solar,xas,xee,t,w,w1,w2,w3,u0a,f0a,fk1,fk2,a4,g4,z4 )

    tkm1 = 0.
    do k = 1, nv1
       if ( k == 1 ) then
          x(1) = 1.0
          x(2) = 1.0
          x(3) = exp ( - fk1(1) * t(1) )
          x(4) = exp ( - fk2(1) * t(1) )
          kk = 1
          xy = 1.0
       else
          kk = k - 1
          y1 = t(kk) - tkm1
          x(1) = exp ( - fk2(kk) * y1 )
          x(2) = exp ( - fk1(kk) * y1 )
          x(3) = 1.0
          x(4) = 1.0
          if (solar) y1 = t(kk)
          xy =  exp ( - y1 / u0a(kk) )
       endif
       if (kk > 1) tkm1 = t(kk)

       do  jj = 1, 4
          fi(jj) = z4(jj,kk) * xy
       end do
       do ii = 1, 4
          fw4 = g4(ii,kk) * x(ii)
          do jj = 1, 4
             fi(jj) = fi(jj) + a4(jj,ii,kk) * fw4
          end do
       end do

       ffu(k)= fw1 * fi(2) + fw2 * fi(1)
       ffd(k)= fw1 * fi(3) + fw2 * fi(4)  + fw3 * xy
    end do

  end subroutine qft
  !> Subroutine rad_init initialize data arrays for gases, ice model and water
  !> model on first call
  subroutine rad_init

    if (.not.fuliou_Initialized) then
       call init_ckd
       randoms = new_RandomNumberSequence(1)
       call init_cldwtr
       fuliou_Initialized = .True.
    end if

    fuliou_Initialized = .True.
  end subroutine rad_init
  !> Subroutine rad: Computes radiative fluxes using a band structure
  !> defined by input ckd file
  subroutine rad (as, u0, ss, pts, ee, pp, pt, ph, po, fds, fus, fdir, fuir, &
       plwc, pre, useMcICA )

    real, intent (in)  :: pp (nv1) ! pressure at interfaces

    real, dimension(nv), intent (in)  :: &
         pt,   & ! temperature [K] at mid points
         ph,   & ! humidity mixing ratio in kg/kg
         po      ! ozone mixing ratio

    real, optional, dimension(nv), intent (in)  :: &
         plwc, & ! cloud liquid water content [g/m^3]
         pre!,  & ! effective radius of cloud droplets [microns]
!          piwc, & ! cloud ice water content [g/m^3]
!          pde,  & ! effective diameter of ice particles [microns]
!          prwc, & ! rain water content [g/m^3]
!          pgwc    ! graupel water content

    real, intent (in) :: &
         as, & ! broadband albedo (all visible bands given this value)
         ee, & ! broadband surface emissivity (all IR bands given this value)
         u0, & ! cosine of solar zenith angle
         ss, & ! Solar constant
         pts   ! Surface skin temperature

    logical, optional, intent (in ) :: useMcICA

    real, dimension(nv1), intent (out)::  &
         fds, fus,  & ! downward and upward solar flux
         fdir, fuir   ! downward and upward ir flux

    logical            :: McICA = .False.

    if(present(useMcICA)) McICA = useMcICA

    call rad_ir(pts, ee, pp, pt, ph, po, fdir, fuir, &
                 plwc, pre, McICA  )
    call rad_vis(as, u0, ss, pp, pt, ph, po, fds, fus,  &
                 plwc, pre, McICA  )

  end subroutine rad

  !> Subroutine rad_ir
  !> Computes IR radiative fluxes using a band structure
  !> defined by input ckd file
  !>
  subroutine rad_ir (pts, ee, pp, pt, ph, po, fdir, fuir, &
       plwc, pre, useMcICA  )
    use modglobal, only : pi

    real, intent (in)  :: pp (nv1) ! pressure at interfaces

    real, dimension(nv), intent (in)  :: &
         pt,   & ! temperature [K] at mid points
         ph,   & ! humidity mixing ratio in kg/kg
         po      ! ozone mixing ratio

    real, optional, dimension(nv), intent (in)  :: &
         plwc, & ! cloud liquid water content [g/m^3]
         pre!,  & ! effective radius of cloud droplets [microns]
!          piwc, & ! cloud ice water content [g/m^3]
!          pde,  & ! effective diameter of ice particles [microns]
!          prwc, & ! rain water content [g/m^3]
!          pgwc    ! graupel water content

    real, intent (in) :: &
         ee, & ! broadband surface emissivity (all IR bands given this value)
         pts   ! Surface skin temperature

    logical, optional, intent (in ) :: useMcICA

    real, dimension(nv1), intent (out)::  &
         fdir, fuir   ! downward and upward ir flux

    ! ----------------------------------------
    logical            :: McICA = .False.
    logical, parameter :: irWeighted = .False.

    real, dimension (nv)   :: tw,ww,tg,dz,tauNoGas, wNoGas, Tau, w
    real, dimension (nv1)  :: fu1, fd1, bf
    real, dimension (nv,4) :: www, pfNoGas, pf

    integer :: ib, ig, k, ig1, ig2, ibandloop, iblimit
    real :: fuq2, xir_norm
    real, dimension(:), allocatable, save :: bandWeights
    real :: randomNumber
    ! ----------------------------------------

    ib=0;ig1=0;ig2=0
    if (.not.fuliou_Initialized) then
       call init_ckd
       randoms = new_RandomNumberSequence(1)
       call init_cldwtr
       fuliou_Initialized = .True.
    end if

    if(.not. allocated(bandweights)) then
      allocate(bandweights(size(ir_bands)))
      call computeIRBandWeights(ir_bands, irWeighted, bandWeights)
    end if
    if(present(useMcICA)) McICA = useMcICA

    fdir(:) = 0.0; fuir(:) = 0.0

    call thicks(pp, pt, ph, dz)

    if (McICA) then
       !
       ! Select a single band and g-point (ib, ig1) and use these as the limits
       !   in the loop through the spectrum below.
       !
       randomNumber = getRandomReal(randoms)
       call select_bandg(ir_bands, bandweights, randomNumber, ib, ig1)
       ig2 = ig1
       iblimit = 1
    else
       iblimit = size(ir_bands)
    end if

    bandLoop: do ibandloop = 1, iblimit
      if (.not. McICA) then
         ib  = ibandloop
         ig1 = 1
         ig2 = kg(ir_bands(ib))
      end if
      !
      ! Water vapor continuum optical depth
      !
      call gascon ( center(ir_bands(ib)), pp, pt, ph, TauNoGas )
      wNoGas = 0.; pfNoGas  = 0.
      if (present(plwc)) then
        call cloud_water(ib + size(solar_bands), pre, plwc, dz, tw, ww, www)
        call combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, tw, ww, www)
      end if

      call planck(pt, pts, llimit(ir_bands(ib)), rlimit(ir_bands(ib)), bf)

      gPointLoop: do ig = ig1, ig2
         tau = TauNoGas; w = wNoGas; pf = pfNoGas
         call gases (ir_bands(ib), ig, pp, pt, ph, po, tg )
         call combineOpticalProperties(tau, w, pf, tg)

         !
         ! Solver expects cumulative optical depth
         !
         do k = 2, nv
           tau(k) = tau(k) + tau(k - 1)
         end do
         call qft (.False., ee, 0., 0., bf, tau, w, pf(:, 1), pf(:, 2),      &
              pf(:, 3), pf(:, 4), fu1, fd1)

         if (McICA) then
            xir_norm = 1./bandweights(ib)
         else
            xir_norm = gPointWeight(ir_bands(ib), ig)
         end if

         fdir(:) = fdir(:) + fd1(:) * xir_norm
         fuir(:) = fuir(:) + fu1(:) * xir_norm
      end do gPointLoop
    end do bandLoop
    !
    ! fuq2 is the surface emitted flux in the band 0 - 280 cm**-1 with a
    ! hk of 0.03.
    !
    fuq2 = bf(nv1) * 0.03 * pi * ee
    fuir(:) = fuir(:) + fuq2
  end subroutine rad_ir
  !> Subroutine rad_vis: Computes radiative fluxes using a band structure
  !> defined by input ckd file
  subroutine rad_vis (as, u0, ss, pp, pt, ph, po, fds, fus,  &
       plwc, pre, useMcICA  )

    real, intent (in)  :: pp (nv1) ! pressure at interfaces

    real, dimension(nv), intent (in)  :: &
         pt,   & ! temperature [K] at mid points
         ph,   & ! humidity mixing ratio in kg/kg
         po      ! ozone mixing ratio

    real, optional, dimension(nv), intent (in)  :: &
         plwc, & ! cloud liquid water content [g/m^3]
         pre!,  & ! effective radius of cloud droplets [microns]
!          piwc, & ! cloud ice water content [g/m^3]
!          pde,  & ! effective diameter of ice particles [microns]
!          prwc, & ! rain water content [g/m^3]
!          pgwc    ! graupel water content

    real, intent (in) :: &
         as, & ! broadband albedo (all visible bands given this value)
         u0, & ! cosine of solar zenith angle
         ss    ! Solar constant

    logical, optional, intent (in ) :: useMcICA

    real, dimension(nv1), intent (out)::  &
         fds, fus    ! downward and upward solar flux

    ! ----------------------------------------
    logical            :: McICA = .False.
    logical, parameter :: solarWeighted = .false.

    real, dimension (nv)   :: tw,ww,tg,tgm,dz, tauNoGas, wNoGas, tau, w
    real, dimension (nv1)  :: fu1, fd1, bf
    real, dimension (nv,4) :: www, pfNoGas, pf

    real, dimension(:), allocatable, save :: bandWeights

    integer :: ib, ig, k, ig1, ig2, ibandloop, iblimit
    real    :: fuq1, xs_norm
    real    :: randomNumber
    ! ----------------------------------------
    ig2=0
    if (.not.fuliou_Initialized) call rad_init

    if (.not. allocated(bandweights)) then
      allocate(bandweights(size(solar_bands)))
      call computeSolarBandWeights(solar_bands, solarWeighted, bandWeights)
    end if
    if(present(useMcICA)) McICA = useMcICA

    fds(:)  = 0.0
    fus(:)  = 0.0
    bf(:)   = 0.0

    if(u0 > minSolarZenithCosForVis) then
      call thicks(pp, pt, ph, dz)

      if (McICA) then
         randomNumber = getRandomReal(randoms)
         !
         ! Select a single band and g-point (ib, ig1) and use these as the
         ! limits in the loop through the spectrum below.
         !
         call select_bandg(solar_bands, bandweights, randomNumber, ib, ig1)
         ig2 = ig1
         iblimit = 1
      else
         iblimit = size(solar_bands)
      end if

      bandLoop: do ibandloop =  1, iblimit
         !
         ! select g points either all, or one depending on McICA
         !
         if (.not. McICA) then
           ib  = ibandloop
           ig1 = 1
           ig2 = kg(solar_bands(ib))
         end if

         !
         ! Rayleigh scattering
         !
         call rayle ( ib, u0, power(solar_bands(ib)), pp, pt, dz, tauNoGas, &
              wNoGas, pfNoGas)
         !
         ! Water vapor continuum
         !
         call gascon ( center(solar_bands(ib)), pp, pt, ph, tgm )
         if(any(tgm > 0.)) &
           call combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, tgm)
         !
         ! Cloud water
         !
         if (present(plwc)) then
           call cloud_water(ib, pre, plwc, dz, tw, ww, www)
           call combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, tw,ww,www)
         end if

         gPointLoop: do ig =  ig1, ig2
           tau = tauNoGas; w = wNoGas; pf = pfNoGas
           call gases (solar_bands(ib), ig, pp, pt, ph, po, tg )
           call combineOpticalProperties(tau, w, pf, tg)

           !
           ! Solver expects cumulative optical depth
           !
           do k = 2, nv
             tau(k) = tau(k) + tau(k - 1)
           end do
           call qft (.true., 0., as, u0, bf, tau, w, pf(:, 1), pf(:, 2),    &
                pf(:, 3), pf(:, 4), fu1, fd1)
           if (McICA) then
              xs_norm =power(solar_bands(ib))/ bandweights(ib)
           else
              xs_norm =gPointWeight(solar_bands(ib), ig)*power(solar_bands(ib))
           end if
           fds(:) = fds(:) + fd1(:) * xs_norm
           fus(:) = fus(:) + fu1(:) * xs_norm
         end do gPointLoop
      end do bandLoop
      !
      ! In this model, we used the solar spectral irradiance determined by
      ! Thekaekara (1973), and 1340.0 W/m**2 is the solar energy contained
      ! in the spectral region 0.2 - 4.0 um., thus scale solar fluxes by
      ! fuq1
      !
      fuq1 = ss / totalpower
      fds(:)  = fds(:)*fuq1
      fus(:)  = fus(:)*fuq1
    end if
  end subroutine rad_vis
  !> Subroutine select_bandg
  !>
  !> selects the band (i) and the g point (j) based on the probability a
  !> photon would be found in the wavelengths covered by the band and g
  !> point, which is given by g_prob.  Note g_prob sums to unity for both
  !> the solar bands (power > 0.) and the infrared bands respectively.
  !>
  subroutine select_bandg(bands, bandweights, randomNumber, i, j)
    type(band_properties), &
          dimension(:),    &
             intent ( in) :: bands
    real, dimension(:), &
             intent ( in) :: bandweights
    real,    intent ( in) :: randomNumber
    integer, intent (out) :: i, j

    real :: cumulative

    i=1; j=1
    ! The probability contained in the first g point of the first band
    cumulative = gPointWeight(bands(i), j) * bandweights(i)

    do while (randomNumber > cumulative .and. cumulative < 1.0)
       j = j+1
       if (j > kg(bands(i)) ) then
          i=i+1
          j=1
       end if
       cumulative = cumulative + gPointWeight(bands(i), j) * bandweights(i)
    end do
  end subroutine select_bandg

  !> subroutine thicks: Integrates the hydrostatic equation to provide
  !> layer thicknesses
  !>
  subroutine thicks(pp, pt, ph, dz)
    use modglobal, only : ep2,grav,rd
    real, intent (in) :: pp(nv1), pt(nv), ph(nv)
    real, intent (out):: dz(nv)

    integer :: i
    real    :: tv

    do  i = 1, nv
       tv = pt(i)*(1+0. + ep2*ph(i) )
       dz(i) = (rd/grav) * tv * alog( pp(i+1) / pp(i) )
    end do

  end subroutine thicks

    !> Adds optical properties to running sum
    !>   If ssa and/or w[1-4] are not present we assume the new medium is
    !> strictly absorbring
    !>
  subroutine combineOpticalProperties(tau,      ssa,      pF, &
                                      tauToAdd, ssaToAdd, pFtoAdd)
    use modglobal, only : eps1
    real, dimension(:),    intent(inout) :: tau, ssa
    real, dimension(:, :), intent(inout) :: pF   ! Phs function (level, moment)
    real, dimension(:),    intent(in)    :: tautoAdd
    real, dimension(:),    optional, intent(in) :: ssaToAdd
    real, dimension(:, :), optional, intent(in) :: pFToAdd ! Phs function

    integer :: j
    if(present(ssaToAdd) .and. present(pfToAdd)) then
       do j = 1, size(pF, 2)
          where (ssa(:) * tau(:) + ssaToAdd(:) * tauToAdd(:) > 0.)
             pf(:, j) = (ssa(:)*tau(:)*pf(:, j) + ssaToAdd(:)*tauToAdd(:)     &
                  * pfToAdd(:, j))/(ssa(:)*tau(:) + ssaToAdd(:) * tauToAdd(:))
          elsewhere
             pf(:, j) = 0.
          end where
       end do
       where (tau(:) + tauToAdd(:) > 0.)
          ssa(:) = (ssa(:) * tau(:) + ssaToAdd(:) * tauToAdd(:)) /            &
               (tau(:) + tauToAdd(:))
       elsewhere
          ssa(:) = 0.
       end where
       tau(:) = tau(:) + tauToAdd(:)
    else
      !
      ! New medium is absorbing - phase function doesn't change
      !
       ssa(:) = (ssa(:) * tau(:)) / (tau(:) + tauToAdd(:))
       tau(:) = tau(:) + tauToAdd(:)
    end if

  end subroutine combineOpticalProperties
  !> Subroutine rayle:  computes optical properties associated with rayleigh
  !> scattering
  !>
  !> ri is the coefficient in Eq.(4.8) of Fu (1991) to compute the optical
  !> depth due to Rayleigh scattering in the solar bands.
  !>
  !> tr, wr, and wwr are the optical depth, single scattering albedo,
  !> and expansion coefficients of the phase function ( 1, 2, 3, and
  !> 4 ) due to the Rayleigh scattering for a given layer.
  !>
  subroutine rayle ( ib, u0, power, pp, pt, dz, tr, wr, wwr)
    integer, intent (in) :: ib
    real, intent (in)    :: u0, power, pp(nv1), pt(nv), dz(nv)
    real, intent (out)   :: tr(nv), wr(nv), wwr(nv,4)

    real, parameter :: ri(6)=(/ 0.9022e-5, 0.5282e-6, 0.5722e-7, &
         0.1433e-7, 0.4526e-8, 0.1529e-8 /)

    integer :: i
    real    :: x

    if ( ib == 1 ) then
       x = -3.902860e-6*u0*u0 + 6.120070e-6*u0 + 4.177440e-6
    else
       x = ri(ib)
    endif

    if(power > 0.) then
      do  i = 1, nv
        tr(i) = x * ( pp(i) + pp(i+1) ) * dz(i) * 0.5 / pt(i)
      end do
      wr(:) = 1.0
      wwr(:, :) = 0.
      wwr(:, 2) = 0.5
    else
      tr(:)     = 0.
      wr(:)     = 0.
      wwr(:, :) = 0.
    end if
  end subroutine rayle

  !> tgm(nv) are the optical depthes due to water vapor continuum absorp-
  !> tion in nv layers for a given band ib. We include continuum absorp-
  !> tion in the 280 to 1250 cm**-1 region. vv(11)-vv(17) are the central
  !> wavenumbers of each band in this region.
  subroutine gascon ( center, pp, pt, ph, tgm)
    real,  intent (in) :: center, pp(nv1), pt(nv), ph(nv)
    real, intent (out) :: tgm(nv)

    integer :: k
    real    :: ff, pe, s, pmid

    if ( center >= 280 .and. center <= 1250.) then
       s = ( 4.18 +  5577.8 * exp ( - 0.00787 * center ) ) / 1013.25
       do k = 1, nv
          pmid   = (pp(k) + pp(k+1))*0.5
          pe     = pmid * ph(k) / ( 0.622 + 0.378 * ph(k) )
          ff     = s*(pe +  0.002*pmid ) *exp (1800.0/pt(k) -6.08108)
          tgm(k) = ff * ph(k) * ( pp(k+1) - pp(k) ) * 1.019767
       end do
    else
       tgm(:) = 0.
    end if
  end subroutine gascon
  !> Subroutine planck:  Integrates planck function over band in xk sub
  !> intervals.  The temperatures at the interfaces are taken as the
  !> average of the mid-point temperatures, the temperature at the lowest
  !> pressure interface is set to the temperature at the first mid point,
  !> and surface temperatures are taken as the skin temperature.
  !>
  subroutine planck ( pt, tskin, llimit, rlimit, bf)
    real, intent (in)    :: pt(nv), tskin, llimit, rlimit
    real, intent (out)   :: bf(nv1) ! intensity [W/m^2/Sr]

    real, parameter :: xk = 10.

    integer :: k
    real    :: v1, v2, vmid, fq1, fq2, tk

    do k = 1, nv1
       bf(k) = 0.0
    end do

    v1 = llimit
    do while (v1 > rlimit+epsilon(rlimit))
       v2 = max(v1 - xk, rlimit)
       vmid = ( v1 + v2 ) * 0.5
       fq1 = 1.19107e-8 * vmid * vmid * vmid
       fq2 = 1.43884 * vmid
       do k = 2, nv
          tk = (pt(k)+pt(k-1))*0.5
          bf(k) = bf(k) + (fq1/(exp(fq2/tk) - 1.0))*(v1-v2)
       end do
       bf(1) = bf(1) + (fq1/(exp(fq2/pt(1)) - 1.0))*(v1-v2)
       bf(nv1) = bf(nv1) + (fq1/(exp(fq2/tskin) - 1.0))*(v1-v2)
       v1 = v2
    end do

  end subroutine planck
  !
    !>
    !> find the weighting for band points so that the probability of a photon
    !> existing in the g-point range of a band can be calculated, and used for
    !> McICA calculations.  This is the relative band width for the IR bands
    !>
  subroutine computeIRBandWeights(bands, weighted, bandweights)
    type(band_properties), &
          dimension(:), intent( in) :: bands
    logical,            intent( in) :: weighted
    real, dimension(:), intent(out) :: bandweights

    integer :: ib

    if(size(bands) /= size(bandweights)) &
      stop "Didn't provide the right amount of storage for band weights"
    if(any(isSolar(bands))) stop "Can't compute IR band weights for solar bands."

    if (weighted) then
       do ib = 1, size(bands)
         bandweights(ib) = (llimit(bands(ib)) - rlimit(bands(ib)))/(bllmx-brlmn)
       end do
    else
       bandweights(:) = 1./(real(size(bands)))
    end if
  end subroutine computeIRBandWeights
    !>
    !> find the weighting for band points so that the probability of a photon
    !> existing in the g-point range of a band can be calculated, and used for
    !> McICA calculations.  This is the solar power for the solar bands
    !>
  subroutine computeSolarBandWeights(bands, weighted, bandweights)
    type(band_properties), &
          dimension(:), intent( in) :: bands
    logical,            intent( in) :: weighted
    real, dimension(:), intent(out) :: bandweights

    integer :: i
    i=1
        if(size(bands) /= size(bandweights)) &
         stop "Didn't provide the right amount of storage for band weights"
    if(any(.not. isSolar(bands))) stop "Can't compute solar band weights in IR"

    if(weighted) then
       bandweights(:) = (/ (power(bands(i))/totalpower, i = 1, size(bands)) /)
    else
       bandweights(:) = 1./(real(size(bands)))
    end if
  end subroutine computeSolarBandWeights
  !
  ! ---------------------------------------------------------------------------
  !> Subroutine ckd_init:  Reads the correlated K distribution data an assures
  !> that it conforms to expected properties
  !>
  subroutine init_ckd
    use modglobal, only : cexpnr
    implicit none

    integer :: i, j, k, l, n, ib, ii, mbs, mbir
    logical :: check
    real    :: bllmx, brlmn
    real, dimension(2000) :: realVars

    real, allocatable :: gasesinband(:)
    character (len=20) :: gasfile
    gasfile = 'ckd.inp.'//cexpnr
    OPEN ( unit = 66, file = trim(gasfile), status = 'old' )
    read (66, '(2I5)') mb, ngases
    allocate (band(mb), gas(ngases))
    read (66, '(300(6E12.4,/))') realVars(1:mb+1)
    do ib = 1, mb
       band(ib)%llimit = realVars(ib)
    end do
    band(mb)%rlimit =  realVars(mb+1)
    read (66, '(300(6E12.4,/))') realVars(1:mb)
    do ib = 1, mb
       band(ib)%power = realVars(ib)
    end do

    mbs = 0.
    do ib=1,mb-1
       band(ib)%rlimit = band(ib+1)%llimit
       band(ib)%center =(band(ib)%rlimit+band(ib)%llimit)*0.5
       if (band(ib)%power > 0.) mbs = mbs + 1
    end do
    band(mb)%center =(band(mb)%rlimit+band(mb)%llimit)*0.5

    if (band(mb)%power > 0.) mbs = mbs + 1
    mbir = mb - mbs
    print 600, trim(gasfile), ngases, mb, mbs, mbir, sum(band%power)

    do n=1,ngases
       read (66,'(A5,I4)') gas(n)%name,gas(n)%iband
       read (66,'(4I4)') gas(n)%noverlap,gas(n)%ng,gas(n)%np,gas(n)%nt
       read (66,'(3E13.5)') gas(n)%mweight, gas(n)%default_conc, gas(n)%tbase
       allocate (gas(n)%hk(gas(n)%ng))
       allocate (gas(n)%sp(gas(n)%np))
       allocate (gas(n)%xk(gas(n)%nt,gas(n)%np,gas(n)%ng,gas(n)%noverlap))
       read (66, '(300(6E12.4,/))') gas(n)%hk(1:gas(n)%ng)
       read (66,'(300(6E12.4,/))') gas(n)%sp(1:gas(n)%np)
       read (66,'(300(6E12.4,/))')                                            &
            realVars(:(gas(n)%ng * gas(n)%np * gas(n)%nt * gas(n)%noverlap))
       ib = 1
       do l = 1, gas(n)%noverlap
         do k = 1, gas(n)%nt
           do j = 1, gas(n)%np
             do i = 1, gas(n)%ng
               gas(n)%xk(k,j,i,l) = realVars(ib)
               ib = ib + 1
             end do
           end do
         end do
       end do

       if (abs(sum(gas(n)%hk) - 1.) <= 1.1 * spacing(1.) ) then
          print 601, gas(n)%name, gas(n)%iband, gas(n)%noverlap,              &
               gas(n)%ng, gas(n)%np, gas(n)%nt
       else
          print *, gas(n)%hk, sum(gas(n)%hk(:))
          stop 'TERMINATING: gas did not occur with probability one in band'
       end if
    end do

    bllmx = tiny(1.)
    brlmn = huge(1.)
    allocate (gasesinband(ngases))
    totalpower = 0.
    do ib=1,mb
       check = .False.
       ii = 0
       do n=1,ngases
          if (gas(n)%iband == ib) then
             check = .True.
             ii = ii+1
             gasesinband(ii) = n
             if (gas(n)%ng > 1 .and. .not. allocated(band(ib)%hk)) then
                allocate(band(ib)%hk(gas(n)%ng))
                band(ib)%hk = gas(n)%hk
                band(ib)%kg = gas(n)%ng
             end if
          end if
       end do
       if (.not.check) stop 'TERMINATING: Gases do not span bands'
       band(ib)%ngases = ii
       allocate(band(ib)%gas_id(ii))
       band(ib)%gas_id(:) = gasesinband(1:ii)
       if (band(ib)%power < epsilon(1.)) then
          bllmx = max(band(ib)%llimit, bllmx)
          brlmn = min(band(ib)%rlimit, brlmn)
       end if
       totalpower = totalpower+band(ib)%power
    end do
    deallocate (gasesinband)

    !
    ! Make separate solar and IR spectra
    !
    allocate(solar_bands(mbs), ir_bands(mbir))
    i = 1; j = 1
    do ib = 1, mb
       if(isSolar(band(ib))) then
          if(i > size(solar_bands)) stop 'TERMINATING: mismatch in solar bands'
          solar_bands(i) = copy_band_properties(band(ib))
          i = i + 1
       else
          if(j > size(ir_bands))    stop 'TERMINATING: mismatch in solar bands'
          ir_bands(j) = copy_band_properties(band(ib))
          j = j + 1
       end if
    end do

    do ib=1,mb
       print 602, ib, band(ib)%power, band(ib)%llimit, band(ib)%rlimit,    &
            band(ib)%ngases, band(ib)%kg
    end do
    print 604

    ckd_Initialized = .True.

600 format ('-----------------------------------------------------------', &
         /3x,'Computing with file: ',A20,' containing ',I3,' gases',   &
         /3x,'Number of Bands (Total/Solar/IR): ',I3,'/',I3,'/',I3     &
         /3x,'Total Power in Solar Bands: ',F8.2)
601 format ('                                                -----------', &
         /3x,'Reading Gas: ',A5,', Band ',I3,', Overlap Type ',I1,     &
         /3x,'# of g-points/ P Levels/ T Coeffs :',3I3)
602 format ('-----------                                                ', &
         /3x,'Band: ',I3,': ',F8.2,' Wm^-2',', between ',F6.0,' and ',F6.0,&
         ' cm^-1',/3x,I3,' gase(s): and ',I3,' g-points')
604 format ('---------------------------------------- Finished band init ')

  end subroutine init_ckd
  !
  ! ---------------------------------------------------------------------------
  !> Subroutine gases:  given an atmospheric state, a band and g point number
  !> this routroutine calculates the optical depth at that g-point within that
  !> band for all the gases (nongray gaseous absorbers) in the band.
  !>
  subroutine gases ( this_band, ig, pp, pt, ph, po, tg )

    implicit none

    type(band_properties), &
             intent (in) :: this_band
    integer, intent (in) :: ig
    real, intent(in)     :: pp(nv1), pt(nv), ph(nv), po(nv)
    real, intent (out)   :: tg(nv)

    real, dimension (nv) :: fkg, fkga, fkgb, pq

    integer :: k, n, igg, nn
    real    :: xfct

    if (.not.ckd_Initialized) stop 'TERMINATING:  ckd_gases not initialized'
    do k = 1, nv
       tg(k) = 0.
    end do
    !
    ! loop through gases finding the cumulative probability for the band and g
    ! point, weight this by the power if in the solar bands, for which power>0.
    !
    do nn=1,this_band%ngases

       n = this_band%gas_id(nn)
       igg=min(ig,gas(n)%ng)
       select case(gas(n)%noverlap)
       case (1)
          igg=min(ig,gas(n)%ng)
          call qk (gas(n)%nt, gas(n)%np, gas(n)%sp, gas(n)%tbase,            &
               gas(n)%xk(1,1,igg,1), pp, pt, fkg )
          call select_gas(gas(n)%name, gas(n)%default_conc, gas(n)%mweight,  &
               pp, ph, po, pq)
          xfct = (2.24e4/gas(n)%mweight) * 10./9.81
          do k = 1, nv
             tg(k) = tg(k) + fkg(k)*pq(k)*(pp(k+1)-pp(k))*xfct
          end do

       case (2)
          call qk (gas(n)%nt, gas(n)%np, gas(n)%sp, gas(n)%tbase,            &
               gas(n)%xk(1,1,igg,1), pp, pt, fkga )
          call qk (gas(n)%nt, gas(n)%np, gas(n)%sp, gas(n)%tbase,            &
               gas(n)%xk(1,1,igg,2), pp, pt, fkgb )
          call select_gas(gas(n)%name, gas(n)%default_conc, gas(n)%mweight,  &
               pp, ph, po, pq)
          do k = 1, nv
             fkg(k) = fkga(k) + pq(k) * fkgb(k)
          end do
          call select_gas('  CO2', gas(n)%default_conc, gas(n)%mweight,       &
               pp, ph, po, pq)
          xfct = (2.24e4/gas(n)%mweight) * 10./9.81
          do k = 1, nv
             tg(k) = tg(k) + fkg(k)*pq(k)*(pp(k+1)-pp(k))*xfct
          end do

       case default
          stop 'TERMINATING: overlap type not supported'
       end select
    end do

  end subroutine gases
  !> Subroutine select_gas determines the mixing ratio to use in the optical
  !> depth calculation.  For simple overlaps, this amounts to selecting the
  !> correct input array, ph for water vapor, po for ozone. For fixed gases
  !> this converts a concentration to a mixing ratio given the molecular weight
  !> of the gas and some specified background concentration.  For C02 the
  !> mixing ratio is chosen so that if conc = 330.e-6, the pq*xfct in subroutine
  !> po = 0.5
  !>
  subroutine select_gas (name, conc_x, mx, pp, ph, po, pq)
    use modglobal, only : mair
    character (len=5), intent (in) :: name
    real, intent (in) :: conc_x,mx
    real, intent (in) :: pp(nv1), ph(nv), po(nv)
    real, intent (out):: pq(nv)

    integer :: k
    real    :: xx

    select case(name)
    case ('  H2O')
       pq(:) = ph(:)
    case ('   O3')
       pq(:) = po(:)
    case ('  CO2')
       xx = conc_x/330.e-6 / ((2.24e4/mx) * 10./9.81)
       pq(:) = xx
    case ('OVRLP')
       do k=1,nv
          if ( pp(k) >= 63.1 ) then
             pq(k) = ph(k)
          else
             pq(k) = 0.0
          endif
       end do
    case default
       xx = conc_x*(mx/mair)
       pq(:) = xx
    end select

    return
  end subroutine select_gas
  !
  ! ---------------------------------------------------------------------------
  !> Subroutine qk: interpolates the gasesous absorption coefficients in units
  !> of (cm-atm)**-1 to the given temperature and pressure in each layer
  !> following: ln k = a + b * ( t - tbase ) + c * ( t - tbase ) ** 2 in
  !> temperature and  linear interpolation in pressure.
  !>
  subroutine qk (nt, np, stanp, tbase, coefki, pp, pt, fkg )

    implicit none

    integer, intent (in) :: np, nt
    real, intent(in)  :: pp(nv1), pt(nv), coefki(nt,np), stanp(np), tbase
    real, intent(out) :: fkg(nv)

    integer :: i1, k
    real    :: x1, x2, y1, y2, xp, pmid

    if (nt < 3) then
       do k = 1, nv
          fkg(k) = coefki(1,1)
       end do
    else
       i1 = 1
       x1 = 0.
       xp = 0.
       do k = 1, nv
          pmid = 0.5*(pp(k)+pp(k+1))
          do while ( pmid >= stanp(i1) .and. i1 < np)
             i1 = i1 + 1
          end do
          y1 = (pt(k)-tbase)
          y2 = y1 * y1
          x2 = exp (coefki(1,i1) + coefki(2,i1)*y1 + coefki(3,i1)*y2)
          if (i1 > 1) then
             x1 = exp (coefki(1,i1-1) + coefki(2,i1-1)*y1 + coefki(3,i1-1)*y2)
             xp = stanp(i1-1)
          endif
          fkg(k) = x1 + (x2-x1)*(pmid-xp)/(stanp(i1)-xp)
       end do
    end if
    return
  end subroutine qk
  !
  ! ----------------------------------------------------------------------
  !
  elemental integer function kg(thisBand)
    type(band_properties), intent(in) :: thisBand

    kg = thisBand%kg
  end function kg
  !
  ! ----------------------------------------------------------------------
  !
  elemental real function power(thisBand)
    type(band_properties), intent(in) :: thisBand

    power = thisBand%power
  end function power
  !
  ! ----------------------------------------------------------------------
  !
  elemental real function llimit(thisBand)
    type(band_properties), intent(in) :: thisBand

    llimit = thisBand%llimit
  end function llimit
  !
  ! ----------------------------------------------------------------------
  !
  elemental real function rlimit(thisBand)
    type(band_properties), intent(in) :: thisBand

    rlimit = thisBand%rlimit
  end function rlimit
  !
  ! ----------------------------------------------------------------------
  !
  elemental real function center(thisBand)
    type(band_properties), intent(in) :: thisBand

    center = (thisBand%llimit + thisBand%rlimit) / 2.
  end function center
  !
  ! ----------------------------------------------------------------------
  !
  elemental real function gPointWeight(thisBand, ig)
    type(band_properties), intent(in) :: thisBand
    integer,               intent(in) :: ig

    gPointWeight = thisBand%hk(ig)
  end function gPointWeight
  !> function isSolar:  Returns True if the band is in the Solar region
  elemental logical function isSolar(thisBand)
    type(band_properties), intent(in) :: thisBand

    isSolar = thisBand%power > 0.
  end function isSolar
  function copy_band_properties(original) result(copy)

    type(band_properties), intent(in) :: original
    type(band_properties)             :: copy

    copy%kg     = original%kg
    copy%ngases = original%ngases
    copy%llimit = original%llimit
    copy%rlimit = original%rlimit
    copy%center = original%center
    copy%power  = original%power

    if(allocated(original%gas_id)) then
       allocate(copy%gas_id(size(original%gas_id)))
       copy%gas_id(:) = original%gas_id(:)
    end if

    if(allocated(original%hk)) then
       allocate(copy%hk(size(original%hk)))
       copy%hk(:) = original%hk(:)
    end if
  end function copy_band_properties
  !> Surbourine cloud_init initialize data arrays for the cloud model,
  !> checking for consistency between band structure of cloud model and CKD
  !>
  subroutine init_cldwtr
    use modglobal, only : cexpnr
    integer, parameter  :: nrec = 21600

    real, dimension(mb) :: cntrs

    integer             :: ib, i, nbands=18, ierr
    character (len=12)  :: frmt
    character (len=20)  :: filenm

    filenm = 'cldwtr.inp.'//cexpnr
    open ( unit = 71, file = filenm, status = 'old', recl=nrec,iostat=ierr)
    if (ierr==0) read (71,'(2I3)') nsizes, nbands
    if (nbands /= mb .or. nsizes*nbands*15 > nrec) &
         stop 'TERMINATING: incompatible cldwtr.dat file'

    allocate (re(nsizes),fl(nsizes),bz(nsizes,mb),wz(nsizes,mb),gz(nsizes,mb))
    write(frmt,'(A1,I2.2,A8)') '(',mb,'E15.7)    '
    if (ierr/=0) call initvar_cldwtr(cntrs,re,fl,bz,wz,gz)
    if (ierr==0) read (71,frmt) (cntrs(i), i=1,mb)
    do i=1,mb
       if (spacing(1.) < abs(cntrs(i)- center(band(i))) ) &
            stop 'TERMINATING: cloud properties not matched to band structure'
    end do

    write(frmt,'(A1,I2.2,A9)') '(',nsizes,'E15.7)   '
    if (ierr==0) read (71,frmt) (re(i), i=1,nsizes)
    if (ierr==0) read (71,frmt) (fl(i), i=1,nsizes)

    write(frmt,'(A1,I4.4,A7)') '(',nsizes*mb,'E15.7) '
    if (ierr==0) read (71,frmt) ((bz(i,ib), i=1,nsizes), ib=1,mb)
    if (ierr==0) read (71,frmt) ((wz(i,ib), i=1,nsizes), ib=1,mb)
    if (ierr==0) read (71,frmt) ((gz(i,ib), i=1,nsizes), ib=1,mb)
    if (ierr==0) close (71)

    if (minval((/bz,wz,gz/)) < 0.) &
         stop 'TERMINATING: cloud properties out of bounds'

    cldwtr_Initialized = .True.

  end subroutine init_cldwtr

  !> Subroutine cloud_water:  calculates the optical depth (tw), single
  !> scattering albedo (ww), and phase function (www(4)) given the cloud
  !> water [g/m^3] and effective radius [microns] by interpolating based on
  !> known optical properties at predefined sizes
  !>
  subroutine cloud_water ( ib, pre, pcw, dz, tw, ww, www )

    implicit none

    integer, intent (in) :: ib
    real, dimension (nv), intent (in) :: pre, pcw, dz
    real, intent (out) :: tw(nv), ww(nv), www(nv,4)

    integer :: k, j, j0, j1
    real    :: gg, wght, cwmks

    if (.not.cldwtr_Initialized) stop 'TERMINATING: Cloud not cldwtr_Initialized'

    do k = 1, nv
       cwmks = pcw(k)*1.e-3
       if ( cwmks .ge. 1.e-8) then
          j = 0
          do while (j<nsizes)
             if (pre(k) > re(j+1)) then
               j = j + 1
             else
               exit
             end if
          end do
          if (j >= 1 .and. j < nsizes) then
             j1 = j+1
             wght = (pre(k)-re(j))/(re(j1)-re(j))
             tw(k) = dz(k) * cwmks * ( bz(j,ib) / fl(j) +   &
                  ( bz(j1,ib) / fl(j1) - bz(j,ib) / fl(j) ) /    &
                  ( 1.0 / re(j1) - 1.0 / re(j) ) * ( 1.0 / pre(k) &
                  - 1.0 / re(j) ) )
             ww(k) = wz(j,ib) + (wz(j1,ib) - wz(j,ib) ) * wght
             gg    = gz(j,ib) + (gz(j1,ib) - gz(j,ib) ) * wght
          else
             j0 = max(j,1)
             tw(k) = dz(k) * cwmks * (bz(j0,ib)/fl(j0))
             ww(k) = wz(j0,ib)
             gg    = gz(j0,ib)
          end if
          www(k,1) = 3.0 * gg
          do j=2,4
             wght = real(2*j+1)/real(2*j-1)
             www(k,j) = www(k,j-1) * gg * wght
          end do
       else
          www(k,:) = 0.0
          tw(k) = 0.0
          ww(k) = 0.0
          gg    = 0.
       end if
    end do

    return
  end subroutine cloud_water

  !> linear interpolation between two points, returns indicies of the
  !> interpolation points and weights
  !>
  subroutine interpolate(x,ny,y,i1,i2,alpha)

    integer, intent (in) :: ny
    real, intent (in)    :: x, y(ny)

    integer, intent (out) :: i1, i2
    real, intent (out)    :: alpha

    if (y(1) < y(2)) stop 'TERMINATING: band centers increasing'

    i2 = 1
    do while (x < y(i2) .and. i2 < ny)
       i2 = i2+1
    end do
    i1 = max(1,i2-1)
    alpha = 1.

    if(i2.ne.i1) alpha = (x-y(i1))/(y(i2)-y(i1))
    if (alpha <0 .or. alpha >1) print 600, x, y(1), y(ny), alpha

    return

600 format(/'CLOUD_INIT WARNING:  Extrapolating because data out of range', &
         /1x,'x = ',F8.1,', ymax = ',F7.0,', ymin =',F7.0,', alpha = ',F6.3)
  end subroutine interpolate
  subroutine initvar_cldwtr(cntrs,re,fl,bz,wz,gz)
  real, dimension(:),intent(out) :: cntrs,re,fl
  real, dimension(:,:),intent(out) :: bz,wz,gz
    cntrs = (/0.3225000E+05, 0.1110000E+05, 0.6475000E+04, 0.4625000E+04, 0.3425000E+04, 0.2675000E+04, 0.2200000E+04, 0.1800000E+04, 0.1550000E+04, 0.1325000E+04, 0.1175000E+04, 0.1040000E+04, 0.8900000E+03, 0.7350000E+03, 0.6050000E+03, 0.4700000E+03, 0.3400000E+03, 0.1400000E+03/)
    re = (/0.5000000E-01, 0.1400000E+00, 0.2200000E+00, 0.2800000E+00, 0.5000000E+00, 0.4700000E+00, 0.1000000E+01, 0.2500000E+01/)
    fl = (/0.5000000E-01, 0.1400000E+00, 0.2200000E+00, 0.2800000E+00, 0.5000000E+00, 0.4700000E+00, 0.1000000E+01, 0.2500000E+01/)
    bz = reshape((/0.1511000E+02, 0.4025000E+02, 0.5981000E+02, 0.7243000E+02, 0.8369000E+02, 0.7399000E+02, 0.1281700E+03, 0.1209100E+03, 0.1574000E+02, 0.4170000E+02, 0.6152000E+02, 0.7447000E+02, 0.8578000E+02, 0.7559000E+02, 0.1304600E+03, 0.1218400E+03, 0.1638000E+02, 0.4352000E+02, 0.6484000E+02, 0.7797000E+02, 0.8731000E+02, 0.7736000E+02, 0.1343000E+03, 0.1240600E+03, 0.1757000E+02, 0.4578000E+02, 0.6644000E+02, 0.8015000E+02, 0.9049000E+02, 0.7990000E+02, 0.1375600E+03, 0.1259200E+03, 0.1819000E+02, 0.4663000E+02, 0.6939000E+02, 0.8220000E+02, 0.9146000E+02, 0.7999000E+02, 0.1382100E+03, 0.1260800E+03, 0.2130000E+02, 0.5188000E+02, 0.7777000E+02, 0.8702000E+02, 0.9491000E+02, 0.8355000E+02, 0.1434600E+03, 0.1284500E+03, 0.2244000E+02, 0.5735000E+02, 0.8441000E+02, 0.1035000E+03, 0.1034900E+03, 0.8417000E+02, 0.1527700E+03, 0.1320700E+03, 0.1832000E+02, 0.5269000E+02, 0.7667000E+02, 0.1003100E+03, 0.1054600E+03, 0.9286000E+02, 0.1578200E+03, 0.1330300E+03, 0.1727000E+02, 0.5044000E+02, 0.7418000E+02, 0.9676000E+02, 0.1053200E+03, 0.9525000E+02, 0.1580700E+03, 0.1344800E+03, 0.1373000E+02, 0.4490000E+02, 0.6770000E+02, 0.9085000E+02, 0.1091600E+03, 0.1054800E+03, 0.1631100E+03, 0.1362100E+03, 0.1030000E+02, 0.3628000E+02, 0.5723000E+02, 0.7643000E+02, 0.1064500E+03, 0.1049000E+03, 0.1617300E+03, 0.1366200E+03, 0.7160000E+01, 0.2640000E+02, 0.4351000E+02, 0.5724000E+02, 0.9255000E+02, 0.9055000E+02, 0.1491000E+03, 0.1351300E+03, 0.6390000E+01, 0.2100000E+02, 0.3381000E+02, 0.4336000E+02, 0.6690000E+02, 0.6358000E+02, 0.1138300E+03, 0.1256500E+03, 0.1033000E+02, 0.3087000E+02, 0.4763000E+02, 0.6033000E+02, 0.7954000E+02, 0.7392000E+02, 0.1274600E+03, 0.1282100E+03, 0.1186000E+02, 0.3564000E+02, 0.5481000E+02, 0.6985000E+02, 0.9039000E+02, 0.8416000E+02, 0.1424900E+03, 0.1352500E+03, 0.1027000E+02, 0.3308000E+02, 0.5181000E+02, 0.6726000E+02, 0.9324000E+02, 0.8860000E+02, 0.1487100E+03, 0.1404200E+03, 0.6720000E+01, 0.2409000E+02, 0.3942000E+02, 0.5168000E+02, 0.8334000E+02, 0.8072000E+02, 0.1401400E+03, 0.1435700E+03, 0.3920000E+01, 0.1476000E+02, 0.2532000E+02, 0.3263000E+02, 0.6085000E+02, 0.5881000E+02, 0.1123000E+03, 0.1456200E+03/),(/nsizes,mb/))
    wz = reshape((/0.9999990E+00, 0.9999990E+00, 0.9999990E+00, 0.9999990E+00, 0.9999980E+00, 0.9999990E+00, 0.9999980E+00, 0.9999970E+00, 0.9997530E+00, 0.9997000E+00, 0.9996670E+00, 0.9996460E+00, 0.9994920E+00, 0.9994700E+00, 0.9993440E+00, 0.9986670E+00, 0.9959140E+00, 0.9949670E+00, 0.9943790E+00, 0.9938420E+00, 0.9913850E+00, 0.9907530E+00, 0.9889080E+00, 0.9748310E+00, 0.9837610E+00, 0.9789810E+00, 0.9765680E+00, 0.9747000E+00, 0.9634660E+00, 0.9599340E+00, 0.9538650E+00, 0.8976900E+00, 0.7029490E+00, 0.6832410E+00, 0.6797230E+00, 0.6690450E+00, 0.6426160E+00, 0.6329960E+00, 0.6297760E+00, 0.5888200E+00, 0.9473430E+00, 0.9296190E+00, 0.9248060E+00, 0.9145570E+00, 0.8771690E+00, 0.8670470E+00, 0.8536610E+00, 0.7374260E+00, 0.9193560E+00, 0.8962740E+00, 0.8859240E+00, 0.8810970E+00, 0.8127720E+00, 0.7816370E+00, 0.7754180E+00, 0.6373410E+00, 0.8747170E+00, 0.8611220E+00, 0.8478500E+00, 0.8516770E+00, 0.7871710E+00, 0.7729520E+00, 0.7531430E+00, 0.6186560E+00, 0.7647500E+00, 0.7524100E+00, 0.7365290E+00, 0.7434350E+00, 0.6712720E+00, 0.6593920E+00, 0.6394920E+00, 0.5499410E+00, 0.8075360E+00, 0.8087000E+00, 0.7959940E+00, 0.8054890E+00, 0.7505770E+00, 0.7555240E+00, 0.7094720E+00, 0.5719890E+00, 0.7533460E+00, 0.7720260E+00, 0.7672730E+00, 0.7770790E+00, 0.7512640E+00, 0.7609730E+00, 0.7125360E+00, 0.5682860E+00, 0.6327220E+00, 0.6763320E+00, 0.6846310E+00, 0.6935520E+00, 0.7079860E+00, 0.7177240E+00, 0.6824300E+00, 0.5528670E+00, 0.2888850E+00, 0.3484890E+00, 0.3716530E+00, 0.3803670E+00, 0.4545400E+00, 0.4657690E+00, 0.4754090E+00, 0.4938810E+00, 0.2618270E+00, 0.3062830E+00, 0.3213400E+00, 0.3330510E+00, 0.3929170E+00, 0.4068760E+00, 0.4174500E+00, 0.4845930E+00, 0.2958040E+00, 0.3399290E+00, 0.3524940E+00, 0.3655020E+00, 0.4162290E+00, 0.4303690E+00, 0.4352670E+00, 0.4913560E+00, 0.3012140E+00, 0.3547460E+00, 0.3693460E+00, 0.3819060E+00, 0.4336020E+00, 0.4473970E+00, 0.4474060E+00, 0.4869680E+00, 0.2437140E+00, 0.3187610E+00, 0.3446420E+00, 0.3527700E+00, 0.4279060E+00, 0.4389790E+00, 0.4459720E+00, 0.4772640E+00, 0.1090120E+00, 0.1872300E+00, 0.2268490E+00, 0.2249760E+00, 0.3313820E+00, 0.3359170E+00, 0.3748820E+00, 0.4570670E+00/),(/nsizes,mb/))
    gz = reshape((/0.8380000E+00, 0.8390000E+00, 0.8440000E+00, 0.8470000E+00, 0.8490000E+00, 0.8600000E+00, 0.8530000E+00, 0.8590000E+00, 0.8090000E+00, 0.8100000E+00, 0.8190000E+00, 0.8230000E+00, 0.8230000E+00, 0.8490000E+00, 0.8330000E+00, 0.8430000E+00, 0.7740000E+00, 0.7870000E+00, 0.7810000E+00, 0.7920000E+00, 0.8120000E+00, 0.8360000E+00, 0.8150000E+00, 0.8330000E+00, 0.8010000E+00, 0.8020000E+00, 0.7930000E+00, 0.7930000E+00, 0.8140000E+00, 0.8290000E+00, 0.8180000E+00, 0.8320000E+00, 0.8770000E+00, 0.8730000E+00, 0.8790000E+00, 0.8800000E+00, 0.8850000E+00, 0.8990000E+00, 0.8910000E+00, 0.9080000E+00, 0.7830000E+00, 0.7690000E+00, 0.7770000E+00, 0.7560000E+00, 0.7640000E+00, 0.7760000E+00, 0.7700000E+00, 0.7970000E+00, 0.8180000E+00, 0.8050000E+00, 0.8240000E+00, 0.8300000E+00, 0.8150000E+00, 0.8010000E+00, 0.8200000E+00, 0.8450000E+00, 0.8100000E+00, 0.8020000E+00, 0.8260000E+00, 0.8400000E+00, 0.8290000E+00, 0.8530000E+00, 0.8400000E+00, 0.8680000E+00, 0.7740000E+00, 0.7660000E+00, 0.7990000E+00, 0.8180000E+00, 0.8150000E+00, 0.8690000E+00, 0.8340000E+00, 0.8690000E+00, 0.7340000E+00, 0.7280000E+00, 0.7670000E+00, 0.7970000E+00, 0.7960000E+00, 0.8710000E+00, 0.8180000E+00, 0.8540000E+00, 0.6930000E+00, 0.6880000E+00, 0.7360000E+00, 0.7720000E+00, 0.7800000E+00, 0.8800000E+00, 0.8080000E+00, 0.8460000E+00, 0.6430000E+00, 0.6460000E+00, 0.6980000E+00, 0.7410000E+00, 0.7590000E+00, 0.8820000E+00, 0.7930000E+00, 0.8390000E+00, 0.5640000E+00, 0.5820000E+00, 0.6370000E+00, 0.6900000E+00, 0.7190000E+00, 0.8710000E+00, 0.7640000E+00, 0.8190000E+00, 0.4660000E+00, 0.4940000E+00, 0.5460000E+00, 0.6090000E+00, 0.6510000E+00, 0.8230000E+00, 0.7010000E+00, 0.7660000E+00, 0.3750000E+00, 0.4100000E+00, 0.4550000E+00, 0.5250000E+00, 0.5830000E+00, 0.7730000E+00, 0.6370000E+00, 0.7100000E+00, 0.2620000E+00, 0.3010000E+00, 0.3340000E+00, 0.4060000E+00, 0.4850000E+00, 0.6950000E+00, 0.5450000E+00, 0.6310000E+00, 0.1440000E+00, 0.1810000E+00, 0.2000000E+00, 0.2560000E+00, 0.3520000E+00, 0.5620000E+00, 0.4130000E+00, 0.5170000E+00, 0.6000000E-01, 0.7700000E-01, 0.8800000E-01, 0.1120000E+00, 0.1810000E+00, 0.3100000E+00, 0.2220000E+00, 0.3270000E+00/),(/nsizes,mb/))

  end subroutine initvar_cldwtr


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



end module modradfull
