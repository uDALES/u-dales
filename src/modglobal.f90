!> \file modglobal.f90
!!  Declares the global constants

!>
!!  Declares the global constants
!>
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
module modglobal

implicit none
save

      ! Simulation dimensions (parconst.f90)
      integer :: imax = 64
      integer :: jtot = 64
      integer :: jmax
      integer :: jsen
      integer :: kmax = 96
      integer :: isen
      integer ::  i1
      integer ::  j1
      integer ::  k1
      integer ::  k2
      integer ::  i2
      integer ::  j2
      integer ::  nsv = 0       !< Number of additional scalar fields
      integer ::  ncosv = 0

      integer ::  ih=3
      integer ::  jh=3
      integer ::  kh=1

      character(50) :: fname_options = 'namoptions'
      integer, parameter :: longint=8
      logical :: lwarmstart = .false.!<   flag for "cold" or "warm" start
      real    :: trestart  = 3600. !<     * each trestart sec. a restart file is written to disk
      integer(kind=longint) :: itrestart !<     * each trestart sec. a restart file is written to disk
      integer(kind=longint)    :: tnextrestart    !<     * each trestart sec. a restart file is written to disk
      character(50) :: startfile    !<    * name of the restart file

      logical :: llsadv   = .false. !<  switch for large scale forcings

      !<  Global constants modconst.f90
      !< File numbers

      integer, parameter :: ifinput    = 1
      integer, parameter :: ifoutput   = 2
      integer, parameter :: ifnamopt   = 3

      real,parameter :: pi       = 3.141592653589793116
      real,parameter :: grav     = 9.81             !<    *gravity acceleration.
      real,parameter :: rd       = 287.04           !<    *gas constant for dry air.
      real,parameter :: rv       = 461.5            !<    *gas constant for water vapor.
      real,parameter :: cp       = 1004.            !<    *specific heat at constant pressure (dry air).
      real,parameter :: rlv      = 2.5e6            !<    *latent heat for vaporisation.
      real,parameter :: ep       = rd/rv            !<    0.622
      real,parameter :: ep2      = rv/rd - 1.       !<    0.61
      !< real,parameter :: cv       = cp-rd            !<    716.96
      real,parameter :: rcp      = rd/cp            !<    0.286
      real,parameter :: cpr      = cp/rd            !<    3.50
      real,parameter :: rlvocp   = rlv/cp           !<    2.49
      real, parameter :: mair    = 28.967          !< Molar mass of air

      real,parameter :: rhow     = 0.998e3          !<    * Density of water
      real,parameter :: pref0    = 1.e5             !<    *standard pressure used in exner function.
      real,parameter :: tmelt    = 273.16           !<    *temperature of melting of ice.
      real,parameter :: es0      = 610.78           !<    * constants used for computation
      real,parameter :: at       = 17.27            !<    * of saturation mixing ratio
      real,parameter :: bt       = 35.86            !<    * using Tetens Formula.
      real,parameter :: ekmin    = 1.e-6            !<    *minimum value for k-coefficient.
      real,parameter :: e12min   = 5.e-5            !<    *minimum value for TKE.
      real,parameter :: fkar     = 0.4              !<    *Von Karman constant
      real,parameter :: eps1     = 1.e-10           !<    *very small number*
      real,parameter :: epscloud = 1.e-5            !<    *limit for cloud calculation 0.01 g/kg
      real,parameter :: boltz    = 5.67e-8          !<    *Stefan-Boltzmann constant

      logical :: lcoriol  = .true.  !<  switch for coriolis force
      integer :: igrw_damp = 1 !< switch to enable gravity wave damping 
      real    :: geodamptime = 7200. !< time scale for nudging to geowind in sponge layer, prevents oscillations
      real    :: om22                       !<    *2.*omega_earth*cos(lat)
      real    :: om23                       !<    *2.*omega_earth*sin(lat)
      real    :: om22_gs                       !<    *2.*omega_earth*cos(lat)
      real    :: om23_gs                       !<    *2.*omega_earth*sin(lat)
      real    :: xlat    = 52.              !<    *latitude  in degrees.
      real    :: xlon    = 0.               !<    *longitude in degrees.


      !Advection scheme

      integer :: iadv_mom = 5, iadv_tke = -1, iadv_thl = -1,iadv_qt = -1,iadv_sv(100) = -1
      integer, parameter :: iadv_upw    = 1
      integer, parameter :: iadv_cd2    = 2
      integer, parameter :: iadv_5th    = 5
      integer, parameter :: iadv_cd6    = 6
      integer, parameter :: iadv_kappa  = 7

      logical :: lmoist   = .true.  !<   switch to calculate moisture fields
      logical :: lsgbucorr= .false.  !<   switch to enable subgrid buoyancy flux


      ! Global variables (modvar.f90)
      real :: xday      = 1.    !<     * day number
      real :: xtime     = 0.    !<     * GMT time
      real :: cu        = 0.    !<     * translation velocity in x-direction
      real :: cv        = 0.    !<     * translation velocity in y-direction
      real :: runtime   = 300.  !<     * simulation time in secs
      real :: dtmax     = 20.    !<     * maximum time integration interval
      integer(kind=longint) :: idtmax        !<     * maximum time integration interval
      real :: dtav_glob   = 60.
      real :: timeav_glob = 3600.
      real :: tres     = 0.001
      real :: thres     = 5.e-3 !<     * threshold value for inversion height calculations
      real :: dqt               !<     * applied gradient of qt at top of model
      real :: dtheta            !<     * applied gradient of theta at top of model
      real,allocatable :: dsv(:)          !<     * applied gradient of sv(n) at top of model
    !<     real :: dsv(nsv)          !<     * applied gradient of sv(n) at top of model

      integer(kind=longint) :: dt                !<     * time integration interval
      real :: rdt                !<     * time integration interval
      integer(kind=longint) :: timee             !<     * elapsed time since the "cold" start
      real :: rtimee             !<     * elapsed time since the "cold" start
      integer(kind=longint) :: btime             !<     * time of (re)start
      integer :: ntimee         !<     * number of timesteps since the cold start
      integer :: ntrun          !<     * number of timesteps since the start of the run
      integer(kind=longint) :: timeleft
      
      logical :: ladaptive   = .false.    !<    * adaptive timestepping on or off

      real    :: courant = -1
      real    :: peclet  = 0.15
      integer(kind=longint) :: dt_lim


      integer :: rk3step = 0

      integer :: iexpnr = 0     !<     * number of the experiment

      character(3) cexpnr



      ! modphsgrd.f90

      real :: dx              !<  grid spacing in x-direction
      real :: dy              !<  grid spacing in y-direction
      real :: dz              !<  grid spacing in z-direction
      real :: dxi             !<  1/dx
      real :: dyi             !<  1/dy
      real :: dzi             !<  1/dz
      real :: dxiq            !<  1/(dx*4)
      real :: dyiq            !<  1/(dy*4)
      real :: dziq            !<  1/(dz*4)
      real :: dxi5            !<  1/(dx*2)
      real :: dyi5            !<  1/(dy*2)
      real :: dzi5            !<  1/(dz*2)
      real :: dx2i            !<  (1/dx)**2
      real :: dy2i            !<  (1/dy)**2


      real :: rslabs
      real, allocatable :: dzf(:)         !<  thickness of full level
      real, allocatable :: dzh(:)         !<  thickness of half level
      real, allocatable :: zh(:)          !<  height of half level [m]
      real, allocatable :: zf(:)          !<  height of full level [m]
      real :: xsize    = -1 !<  domain size in x-direction
      real :: ysize    = -1 !<  domain size in y-direction
      real, allocatable :: delta(:)       !<  (dx*dy*dz)**(1/3)

      logical :: leq      = .true.  !<  switch for (non)-equidistant mode.
      logical :: lmomsubs = .false.  !<  switch to apply subsidence on the momentum or not
      character(80) :: author='', version='DALES 3.2'
contains

!> Initialize global settings.
!!
!! Set courant number, calculate the grid sizes (both computational and physical), and set the coriolis parameter
  subroutine initglobal
    use modmpi, only: nprocs, myid,comm3d, my_real, mpierr
    implicit none

    integer :: advarr(4)
    real phi, colat, silat, omega, omega_gs
    integer :: k, n
    character(80) chmess

    !timestepping
    if (courant<0) then
      select case(iadv_mom)
      case(iadv_cd2)
        courant = 3
      case(iadv_cd6)
        courant = 1.4
      case(iadv_5th)
        courant = 1.4
      case default
        courant = 1.4
      end select
    end if


    ! phsgrid

    jmax = jtot/nprocs
    isen = imax/nprocs
    jsen = jmax
    i1=imax+1
    j1=jmax+1
    k1=kmax+1
    k2=kmax+2
    i2=imax+2
    j2=jmax+2
    !set the number of ghost cells. NB: This switch has to run in order of required ghost cells
    advarr = (/iadv_mom,iadv_tke,iadv_thl,iadv_qt/)
    if     (any(advarr==iadv_cd6).or.any(iadv_sv(1:nsv)==iadv_cd6)) then
      ih = 3
      jh = 3
      kh = 1
    elseif (any(advarr==iadv_5th).or.any(iadv_sv(1:nsv)==iadv_5th)) then
      ih = 3
      jh = 3
      kh = 1
    elseif (any(advarr==iadv_kappa).or.any(iadv_sv(1:nsv)==iadv_kappa)) then
      ih = 2
      jh = 2
      kh = 1
    elseif (any(advarr==iadv_cd2).or.any(iadv_sv(1:nsv)==iadv_cd2)) then
      ih = 1
      jh = 1
      kh = 1
    end if
    ncosv = max(2*nsv-3,0)


    ! Global constants

    ! Select advection scheme for scalars. If not set in the options file, the momentum scheme is used
    if (iadv_tke<0) iadv_tke = iadv_mom
    if (iadv_thl<0) iadv_thl = iadv_mom
    if (iadv_qt<0)  iadv_qt  = iadv_mom

    !CvH remove where
    !where (iadv_sv<0)  iadv_sv  = iadv_mom
    do n = 1, nsv
      if(iadv_sv(n) < 0) then
        iadv_sv(n) = iadv_mom
      end if
    end do

    phi    = xlat*pi/180.
    colat  = cos(phi)
    silat  = sin(phi)
    if (lcoriol) then
      omega = 7.292e-5
      omega_gs = 7.292e-5
    else
      omega = 0.
      omega_gs = 7.292e-5
    end if
    om22   = 2.*omega*colat
    om23   = 2.*omega*silat
    om22_gs   = 2.*omega_gs*colat
    om23_gs   = 2.*omega_gs*silat

    ! Variables
    allocate(dsv(nsv))
    write(cexpnr,'(i3.3)') iexpnr


    ! Create the physical grid variables
    allocate(dzf(k1))
    allocate(dzh(k1))
    allocate(zh(k1))
    allocate(zf(k1))
    allocate(delta(k1))


    rslabs = real(imax*jtot)

    dx = xsize / float(imax)
    dy = ysize / float(jtot)

    ! MPI

    ! Note, that the loop for reading zf and calculating zh
    ! has been split so that reading is only done on PE 1


    if(myid==0)then
      open (ifinput,file='prof.inp.'//cexpnr)
      read(ifinput,'(a72)') chmess
      read(ifinput,'(a72)') chmess

      do k=1,kmax
        read(ifinput,*) zf(k)
      end do
      close(ifinput)

    end if ! end if myid==0

  ! MPI broadcast kmax elements from zf

    call MPI_BCAST(zf,kmax,MY_REAL   ,0,comm3d,mpierr)

    zh(1) = 0.0
    do k=1,kmax
      zh(k+1) = zh(k) + 2.0*(zf(k)-zh(k))
    end do
    zf(k1)  = zf(kmax)+ 2.0*(zh(k1)-zf(kmax))


    do  k=1,kmax
      dzf(k) = zh(k+1) - zh(k)
    end do
    dzf(k1) = dzf(kmax)

    dzh(1) = 2*zf(1)
    do k=2,k1
      dzh(k) = zf(k) - zf(k-1)
    end do

    do k=1,k1

      delta(k) = (dx*dy*dzf(k))**(1./3.)
    end do

  !--------------------------------------------------
  ! *** Check whether the grid is equidistant *****
  !--------------------------------------------------

    leq=.true.
    dz = dzf(1)
    do k=1,k1
      if (dzf(k)/=dz) then
        leq = .false.
      end if
    end do

  ! MPI

    if(myid==0)then
      if (.not.leq) then
        write(6,*) &
            'WARNING, You are working with a non-equidistant grid!!!!'
      end if
    end if ! end if myid==0

    dxi     = 1./dx
    dyi     = 1./dy
    dzi     = 1./dz

    dxiq    = 0.25*dxi
    dyiq    = 0.25*dyi
    dziq    = 0.25*dzi

    dx2i    = dxi*dxi
    dy2i    = dyi*dyi

    dxi5    = 0.5*dxi
    dyi5    = 0.5*dyi
    dzi5    = 0.5*dzi

    if(myid==0)then
      write (6,*) 'lev    dz     zf      zh       dzh    delta'
      do k=k1,1,-1
        write(6,'(i4,5f8.2)') k,dzf(k),zf(k),zh(k),dzh(k),delta(k)
      end do
    end if
    tnextrestart = trestart/tres
    timeleft     = btime+runtime

  end subroutine initglobal
!> Clean up when leaving the run
  subroutine exitglobal
    deallocate(dsv,dzf,dzh,zh,zf,delta)
  end subroutine exitglobal

end module modglobal
