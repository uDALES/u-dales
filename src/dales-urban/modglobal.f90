!> \file modglobal.f90
!!  Declares the global constants

!> 
!! \author Jasper Tomas, TU Delft 31 March 2014

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

  integer :: poisrcheck = 0   ! switch to check if it is the first (RK) time step
  ! Simulation dimensions (parconst.f90)
  integer :: imax = 64
  integer :: jtot = 64
  integer :: jmax
  integer :: jsen
  integer :: kmax = 96
  integer :: isen
  integer ::  ib
  integer ::  ie
  integer ::  jb
  integer ::  je
  integer ::  jgb           ! global j range
  integer ::  jge           ! global j range
  integer ::  offset
  integer ::  kb 
  integer ::  ke
  integer ::  nsv = 0       !< Number of additional scalar fields
  integer ::  ncosv = 0
  integer ::  nsvl = 0
  integer ::  nsvp = 0
  integer ::  nvar=10
  integer ::  ih=3
  integer ::  jh=3
  integer ::  kh=1
  integer ::  ihc=2         ! used in k-scheme
  integer ::  jhc=2         ! used in k-scheme
  integer ::  khc=2         ! used in k-scheme

  integer :: nxwalls=0        ! no. of x-walls in IBM
  integer :: nywalls=0        ! no. of y-walls in IBM
  integer :: nzwalls=0        ! no. of z-walls in IBM
  integer :: nblocks=0        ! no. of blocks in IBM 

  integer ::  iplane        ! ib+iplane is the plane that is stored when lstoreplane=.true.
  integer ::  linletgen  = 0 !<  0: no inletgen, 1: turb. inlet generator (Lund (1998)), 2: read inlet from file
  integer ::  nstore = 1002  ! number of rk steps in inletfile. This should be a multiple of three!
  character(50) :: fname_options = 'namoptions'
  integer, parameter :: longint=8
  logical :: lwarmstart = .false.!<   flag for "cold" or "warm" start
  logical :: lreadscal  = .false.!<   flag for reading scalar pollutant field (warm start)
  real    :: trestart           !<     * each trestart sec. a restart file is written to disk
  real    :: tfielddump         !<
  real    :: trestart2          !<     * each trestart2 sec. a data file is written to disk
  real    :: tnextrestart       !<     * each trestart sec. a restart file is written to disk
  real    :: tnextrestart2      !<     * each trestart2 sec. a data file is written to disk
  real    :: tscale             !       timescale: domain height*Uinf/utau**2
  real    :: tnextfielddump     !<
  character(50) :: startfile    !<    * name of the restart file

  logical :: llsadv     = .false. !<  switch for large scale forcings
  logical :: linoutflow = .false. !<  switch for periodic BC in both horizontal directions (false) or inflow/outflow in i and periodic in j.
  logical :: lzerogradtop = .false. !<  switch for zero gradient BC's at top wall (linletgen 1 and 2 are seperate).
  logical :: lzerogradtopscal = .false.  !
  logical :: ltfluxtop  = .false. !fixed flux or fixed value top BC
  logical :: ltfluxbot  = .false. !fixed flux or fixed values bottom BC
  logical :: lbuoyancy  = .true.  !<  switch for buoyancy force in modforces
  logical :: lneutral  = .false.
  logical :: ltempeq    = .true.  !<  switch for solving temperature equation (either with or without buoyancy term)
  logical :: lscalinout = .false. !<  seperate switch for inflow/outflow BC for scalar (only necessary when linoutflow==.false.).
  logical :: lscalrec   = .false. !<  
  logical :: ltempinout = .false. !<  seperate switch for inflow/outflow BC for temperature (only necessary when linoutflow==.false.).
  logical :: lmoistinout= .false. !<  seperate switch for inflow/outflow BC for moisture (only necessary when linoutflow==.false.).
  logical :: lper2inout = .false. !<  switch that determines type of restart: .true. means switching from periodic to in/outflow: inlet profile is read from prof.inp
  logical :: libm       = .false. !<  switch that determines wether the Immersed Boundary Method is turned on
  logical :: lwalldist  = .true.  !<  switch that determines wether the wall distances should be computed
  logical :: lles       = .true.  !<  switch that determines wether the subgrid model is turned on or constant ekm and ekh are used (DNS)
  logical :: ltec2d     = .false. !   2D tec files written
  logical :: ltec3d     = .false. !<  switch that determines wether tec3d files are written
  logical :: linletRA   = .false. !<  switch that determines wether a Running Average should be used (.true.) in inlet generator
  logical :: lfixinlet  = .false. !<  switch that determines wether the average inlet profiles can evolve or not (only used when linletgen=1,2)
  logical :: lfixutauin = .false. !<  switch that determines wether the utau is kept fixed at the inlet (only used when linletgen=1,2)
  logical :: lscasrc   = .false.  !
  logical :: lreadminl = .false. !<  switch for reading mean inlet/recycle plane profiles (used in inletgenerator)
  logical :: lwallfunc  = .false. !<  switch that determines wether wall functions are used to compute the wall-shear stress 
  logical :: lwallfuncs  = .false. !<  switch that determines wether wall functions are used to compute the wall-shear stress 
  logical :: lMOST      = .false. !< MOST or Werner-Wengle wall-function
  logical :: lmassflowr = .false. !<  switch that determines wether u-velocity is corrected to get a fixed mass flow rate
  logical :: lstoreplane =.false. !<  switch that determines wether i-plane data is stored.
  logical :: lstore3d    =.false. !<  switch that determines wether 3d fields are stored in subdir's each trestart2.
  logical :: lstorexz    =.false. !<  switch that determines wether xz fields are stored in subdir's each trestart2.
  logical :: lstorexy   =  .false. !xy files stored
  logical :: lreadmean  = .false. !<  switch that determines wether mean variables should be read from means#myid#.#expnr#
  logical :: lnetcdf    = .true. !<  write NETCDFs
  logical :: ifixuinf   = .true. !dpdxl relaxed to have Uinf 1. dpdx = (1/dt)*(Uh-Uinf)2. d/dt(dpdx) = 1/tau*(Uh-Uinf)                                                                              
  
  real    :: freestreamav = 0.  !
  real    :: freestrtmpav = 0.   !
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
  real, parameter :: mair    = 28.967           !< Molar mass of air


  real             :: numol                     !< kinematic viscosity for couette flow Re=5000 (Re=Uinf*H/(2*nu)) H=1, Uinf=1
  real             :: numoli                    !< 1/numol

  real             :: prandtlmol                !< Prandtl number (for air at 300K). Fluid property! 
  real             :: prandtlmoli               !< Inverse of Prandtl number

  real,parameter :: rhow     = 0.998e3          !<    * Density of water
  real,parameter :: pref0    = 1.e5             !<    *standard pressure used in exner function.
  real,parameter :: tmelt    = 273.16           !<    *temperature of melting of ice.
  real,parameter :: es0      = 610.78           !<    * constants used for computation
  real,parameter :: at       = 17.27            !<    * of saturation mixing ratio
  real,parameter :: bt       = 35.86            !<    * using Tetens Formula.
  !      real,parameter :: ekmin    = 1.e-6            !<    *minimum value for k-coefficient.
  real,parameter :: ekmin    = 1.e-12           !<    *minimum value for k-coefficient.
  real,parameter :: e12min   = 5.e-5            !<    *minimum value for TKE.
  real,parameter :: fkar     = 0.4              !<    *Von Karman constant
  real,parameter :: eps1     = 1.e-10           !<    *very small number*
  real,parameter :: epscloud = 1.e-5            !<    *limit for cloud calculation 0.01 g/kg
  real,parameter :: boltz    = 5.67e-8          !<    *Stefan-Boltzmann constant

  logical :: lcoriol  = .true.                  !<  switch for coriolis force
  integer :: igrw_damp = 2                      !< switch to enable gravity wave damping 
  real    :: geodamptime = 7200.                !< time scale for nudging to geowind in sponge layer, prevents oscillations
  real    :: massflowrate=    1.                !< fixed mass flow rate used for u-velocity correction
  real    :: Uinf=    0.                        !< fixed U_inf (used in inlet generator)
  real    :: inletav=         0.                !< averaging interval for inlet generator
  real    :: totinletav=      0.                !< averaging interval for inlet generator (used in Running Average)
  real    :: om22                               !<    *2.*omega_earth*cos(lat)
  real    :: om23                               !<    *2.*omega_earth*sin(lat)
  real    :: om22_gs                            !<    *2.*omega_earth*cos(lat)
  real    :: om23_gs                            !<    *2.*omega_earth*sin(lat)
  real    :: xlat    = 52.                      !<    *latitude  in degrees.
  real    :: xlon    = 0.                       !<    *longitude in degrees.


  !scalar source in fluid domain
  integer :: xS = 0, yS = 0, zS=0
  real    :: SS = 0.
  real    :: sigS = 0.


  !Advection scheme

  integer :: iadv_mom = 5, iadv_tke = -1, iadv_thl = -1,iadv_qt = -1,iadv_sv(100) = -1
  integer, parameter :: iadv_upw    = 1
  integer, parameter :: iadv_cd2    = 2
  integer, parameter :: iadv_5th    = 5
  integer, parameter :: iadv_cd6    = 6
  integer, parameter :: iadv_62     = 62
  integer, parameter :: iadv_52     = 52
  integer, parameter :: iadv_kappa  = 7

  logical :: lmoist   = .false.  !<   switch to calculate moisture fields
  logical :: lsgbucorr= .false.  !<   switch to enable subgrid buoyancy flux


  ! Global variables (modvar.f90)
  real :: xday      = 1.    !<     * day number
  real :: xtime     = 0.    !<     * GMT time
  real :: cu        = 0.    !<     * translation velocity in x-direction
  real :: cv        = 0.    !<     * translation velocity in y-direction
  real :: runtime   = 300.  !<     * simulation time in secs
  real :: dtmax     = 20.    !<     * maximum time integration interval
  !      integer(kind=longint) :: idtmax        !<     * maximum time integration interval
  real :: dtav_glob   = 60.
  real :: timeav_glob = 3600.
  real :: totavtime   = 0    !<    * the total time over which the values are averaged in meansXXX.XXX
  real :: thres     = 5.e-3 !<     * threshold value for inversion height calculations
  real :: dqt               !<     * applied gradient of qt at top of model
  real :: dtheta            !<     * applied gradient of theta at top of model
  real,allocatable :: dsv(:)          !<     * applied gradient of sv(n) at top of model
  !<     real :: dsv(nsv)          !<     * applied gradient of sv(n) at top of model

  !     integer(kind=longint) :: dt                !<     * time integration interval
  real ::  dt                !<     * time integration interval
  !      integer(kind=longint) :: timee             !<     * elapsed time since the "cold" start
  real :: timee             !<     * elapsed time since the "cold" start
  !      integer(kind=longint) :: btime             !<     * time of (re)start
  real :: btime            !<     * time of (re)start
  real :: startmean      !
  real :: runavtime        !<     * time of starting running average
  integer :: ntimee         !<     * number of timesteps since the cold start
  integer :: ntrun          !<     * number of timesteps since the start of the run
  real    :: timeleft

  logical :: ladaptive   = .false.    !<    * adaptive timestepping on or off

  real    :: courant = -1
  !      real    :: peclet  = 0.15
  real    :: peclet  = 0.25
  !      real    :: peclet  = 0.2
  !      real    :: peclet  = 0.4
  !integer(kind=longint) :: dt_lim
  real    :: dt_lim


  integer :: rk3step = 0

  integer :: iexpnr = 0     !<     * number of the experiment

  character(3) cexpnr

   real :: thlsrc  = 0.  

   integer :: kplane(100)             ! k-index of planes that are stored in time
 integer :: nkplane = 0             ! number of kplanes being stored  

  ! modphsgrd.f90

  real :: dy              !<  grid spacing in y-direction
  real :: dy2             !<  grid spacing in y-direction squared
  real :: dz              !<  grid spacing in z-direction
  real :: dyi             !<  1/dy
  real :: dyiq            !<  1/(dy*4)
  real :: dyi5            !<  1/(dy*2)
  real :: dy2i            !<  (1/dy)**2


  real :: rslabs
  real, allocatable :: dzf(:)         !<  thickness of full level
  real, allocatable :: dzfc(:)        !<  thickness of full level (extra ghost nodes (used in k-scheme)
  real, allocatable :: dzfci(:)       !<  1/dzfc
  real, allocatable :: dzf2(:)        !<  thickness of full level squared
  real, allocatable :: dzh(:)         !<  thickness of half level
  real, allocatable :: zh(:)          !<  height of half level [m]
  real, allocatable :: zf(:)          !<  height of full level [m]
  real, allocatable :: dzfi(:)         !<  1/dzf
  real, allocatable :: dzfiq(:)        !<  0.25*(1/dzf)
  real, allocatable :: dzfi5(:)        !<  0.5*(1/dzf)
  real, allocatable :: dzhi(:)         !<  1/dzh
  real, allocatable :: dzhci(:)        !<  1/dzh (extra ghost nodes (used in k-scheme)
  real, allocatable :: dzhiq(:)        !<  0.25*(1/dzh)
  real, allocatable :: dzh2i(:)        !<  1/dzh^2
  real, allocatable :: zhi(:)          !<  1/zh
  real, allocatable :: zfi(:)          !<  1/zf
  real, allocatable :: dxf(:)         !<  thickness of full level
  real, allocatable :: dxfc(:)        !<  thickness of full level (extra ghost nodes (used in k-scheme)
  real, allocatable :: dxfci(:)       !<  1/dxfc
  real, allocatable :: dxf2(:)        !<  thickness of full level squared
  real, allocatable :: dxfi(:)        !<  = 1/dxf
  real, allocatable :: dxfiq(:)       !<  = 0.25*(1/dxf)
  real, allocatable :: dxfi5(:)       !<  = 0.5*(1/dxf)
  real, allocatable :: dxh(:)         !<  thickness of half level
  real, allocatable :: dxhi(:)        !<  = 1/dxh
  real, allocatable :: dxhci(:)       !<  = 1/dxh (with extra ghost nodes (used in k-scheme))
  real, allocatable :: dxhiq(:)       !<  = 0.25*(1/dxh)
  real, allocatable :: dxh2i(:)       !<  = 1/dxh^2
  real, allocatable :: xh(:)          !<  height of half level [m]
  real, allocatable :: xf(:)          !<  height of full level [m]
  real :: xsize    = -1 !<  domain size in x-direction
  real :: ysize    = -1 !<  domain size in y-direction
  real, allocatable :: delta(:,:)       !<  (dx*dy*dz)**(1/3)

  real, allocatable :: scalsrc(:,:,:)   !<   i,k plane with scalar sources
  real, allocatable :: scalloc(:,:)     !<   x/z locations of line source + sigma_x and sigma_z
real, allocatable :: scalptsrc(:,:,:,:) !<   i,j,k distribution from scalar point sources
real, allocatable :: scalptloc(:,:)     !<   x/y/z locations of point source+sigma


  logical :: leq      = .true.  !<  switch for (non)-equidistant mode.
  logical :: lmomsubs = .false.  !<  switch to apply subsidence on the momentum or not
  character(80) :: author='', version='DALES 3.2'
contains

  !> Initialize global settings.
  !!
  !! Set courant number, calculate the grid sizes (both computational and physical), and set the coriolis parameter
  subroutine initglobal
    use modmpi, only : nprocs, myid,comm3d, my_real, mpierr
    implicit none

    integer :: advarr(4)
    real phi, colat, silat, omega, omega_gs
    integer :: i, k, n
    character(80) chmess

    !timestepping
    if (courant<0) then
       select case(iadv_mom)
       case(iadv_cd2)
          courant = 1.5
       case(iadv_cd6)
          courant = 1.1
       case(iadv_62)
          courant = 1.1
       case(iadv_5th)
          courant = 1.4
       case(iadv_52)
          courant = 1.4
       case default
          courant = 1.4
       end select
       if (any(iadv_sv(1:nsv)==iadv_cd6) .or. any((/iadv_thl,iadv_qt,iadv_tke/)==iadv_cd6)) then
          courant = min(courant, 1.1)
       elseif (any(iadv_sv(1:nsv)==iadv_62) .or. any((/iadv_thl,iadv_qt,iadv_tke/)==iadv_62)) then
          courant = min(courant, 1.1)
       elseif (any(iadv_sv(1:nsv)==iadv_kappa) .or. any((/iadv_thl,iadv_qt,iadv_tke/)==iadv_kappa)) then
          courant = min(courant, 1.1)
       elseif (any(iadv_sv(1:nsv)==iadv_upw) .or. any((/iadv_thl,iadv_qt,iadv_tke/)==iadv_upw)) then
          courant = min(courant, 1.1)
       elseif (any(iadv_sv(1:nsv)==iadv_5th) .or. any((/iadv_thl,iadv_qt,iadv_tke/)==iadv_5th)) then
          courant = min(courant, 1.4)
       elseif (any(iadv_sv(1:nsv)==iadv_52 ).or. any((/iadv_thl,iadv_qt,iadv_tke/)==iadv_52)) then
          courant = min(courant, 1.4)
       elseif (any(iadv_sv(1:nsv)==iadv_cd2) .or. any((/iadv_thl,iadv_qt,iadv_tke/)==iadv_cd2)) then
          courant = min(courant, 1.5)
       end if
    end if

    ! phsgrid

    jmax = jtot/nprocs
    isen = imax/nprocs
    jsen = jmax
    !set the number of ghost cells. NB: This switch has to run in order of required ghost cells
    advarr = (/iadv_mom,iadv_tke,iadv_thl,iadv_qt/)
    if     (any(advarr==iadv_cd6).or.any(iadv_sv(1:nsv)==iadv_cd6)) then
       ih = 3
       jh = 3
       kh = 1
    elseif (any(advarr==iadv_62).or.any(iadv_sv(1:nsv)==iadv_62)) then
       ih = 3
       jh = 3
       kh = 1
    elseif (any(advarr==iadv_5th).or.any(iadv_sv(1:nsv)==iadv_5th)) then
       ih = 3
       jh = 3
       kh = 1
    elseif (any(advarr==iadv_52).or.any(iadv_sv(1:nsv)==iadv_52)) then
       ih = 3
       jh = 3
       kh = 1
       !    elseif (any(advarr==iadv_kappa).or.any(iadv_sv(1:nsv)==iadv_kappa)) then
    elseif (any(advarr==iadv_kappa)) then
       ih = 2
       jh = 2
       kh = 1
    elseif (any(advarr==iadv_cd2).or.any(iadv_sv(1:nsv)==iadv_cd2)) then
       ih = 1
       jh = 1
       kh = 1
       ihc = 1
       jhc = 1
       khc = 1
    end if

    ! J. Tomas added this for using only kappa scheme for sv(:)
    if (any(iadv_sv(1:nsv)==iadv_kappa)) then
       ihc = 2
       jhc = 2
       khc = 2
    end if


    offset = 1         
    ib = 2-offset
    ie = imax+1-offset
    jb = 2-offset
    je = jmax+1-offset
    jgb= jb               ! global j range (starting at the same as j as the processor j range) 
    jge= jtot+1-offset    ! global j range
    kb = 1-offset
    ke = kmax-offset



    ncosv = max(2*nsv-3,0)

    ! Global constants

    ! Select advection scheme for scalars. If not set in the options file, the momentum scheme is used
    if (iadv_tke<0) iadv_tke = iadv_mom
    if (iadv_thl<0) iadv_thl = iadv_mom
    if (iadv_qt<0)  iadv_qt  = iadv_mom

    !CvH remove where
    !where (iadv_sv<0)  iadv_sv  = iadv_mom
    
    ! tg3315 commented out
    !do n = 1, nsv
    !   if(iadv_sv(n) < 0) then
    !      iadv_sv(n) = iadv_mom
    !   end if
    !end do
    
    !tg3315 added
    do n = 1, nsv
        iadv_sv(n) = iadv_kappa
    end do
    !ends here

    phi    = xlat*pi/180.
    colat  = cos(phi)
    silat  = sin(phi)
    if (lcoriol) then
       omega = 7.292e-5
       omega_gs = 7.292e-5
    else
       omega = 0.
       omega_gs = 0.
    end if
    om22   = 2.*omega*colat
    om23   = 2.*omega*silat
    om22_gs   = 2.*omega_gs*colat
    om23_gs   = 2.*omega_gs*silat

    ! Variables
    allocate(dsv(nsv))
    write(cexpnr,'(i3.3)') iexpnr


    ! Create the physical grid variables
    allocate(dzf(kb-kh:ke+kh))
    allocate(dzf2(kb-kh:ke+kh))
    allocate(dzfi(kb-kh:ke+kh))
    allocate(dzfiq(kb-kh:ke+kh))
    allocate(dzfi5(kb-kh:ke+kh))
    allocate(dzh(kb:ke+kh))
    allocate(dzhi(kb:ke+kh))
    allocate(dzhiq(kb:ke+kh))
    allocate(dzh2i(kb:ke+kh))
    allocate(zh(kb:ke+kh))
    allocate(zf(kb:ke+kh))

    allocate(dxf(ib-ih:ie+ih))
    allocate(dxf2(ib-ih:ie+ih))
    allocate(dxfi(ib-ih:ie+ih))
    allocate(dxfiq(ib-ih:ie+ih))
    allocate(dxfi5(ib-ih:ie+ih))
    allocate(dxh(ib:ie+ih))
    allocate(dxhi(ib:ie+ih))
    allocate(dxhiq(ib:ie+ih))
    allocate(dxh2i(ib:ie+ih))
    allocate(xh(ib:ie+ih))
    allocate(xf(ib:ie+ih))
    allocate(delta(ib-ih:ie+ih,kb:ke+kh))

    rslabs = real(imax*jtot)

    dy = ysize / float(jtot)

    ! MPI

    ! Note, that the loop for reading zf and calculating zh
    ! has been split so that reading is only done on PE 1


    if(myid==0)then
       open (ifinput,file='prof.inp.'//cexpnr)
       read(ifinput,'(a72)') chmess
       read(ifinput,'(a72)') chmess

       do k=kb,ke
          read(ifinput,*) zf(k)
       end do
       close(ifinput)

       ! J. Tomas: Read the x-coordinates of the cell centers from xgrid.inp.XXX
       open (ifinput,file='xgrid.inp.'//cexpnr)
       read(ifinput,'(a72)') chmess
       read(ifinput,'(a72)') chmess

       do i=ib,ie
          read(ifinput,*) xf(i)
       end do
       close(ifinput)

    end if ! end if myid==0

    ! MPI broadcast kmax elements from zf
    call MPI_BCAST(zf,kmax,MY_REAL   ,0,comm3d,mpierr)
    ! MPI broadcast imax elements from xf
    call MPI_BCAST(xf,imax,MY_REAL   ,0,comm3d,mpierr)

    zh(kb) = 0.0
    do k=kb,ke
       zh(k+1) = zh(k) + 2.0*(zf(k)-zh(k))
    end do
    zf(ke+kh)  = zf(ke)+ 2.0*(zh(ke+kh)-zf(ke))


    do  k=kb,ke
       dzf(k) = zh(k+1) - zh(k)
    end do
    dzf(ke+1) = dzf(ke)
    dzf(kb-1) = dzf(kb)

    dzh(kb) = 2*zf(kb)
    do k=kb+1,ke+kh
       dzh(k) = zf(k) - zf(k-1)
    end do


    ! j. tomas: same trick for x-direction...
    xh(ib) = 0.0
    do i=ib,ie
       xh(i+1) = xh(i) + 2.0*(xf(i)-xh(i))
    end do
    xf(ie+ih)  = xf(ie)+ 2.0*(xh(ie+ih)-xf(ie))


    do  i=ib,ie
       dxf(i) = xh(i+1) - xh(i)
    end do
    dxf(ie+1) = dxf(ie)
    dxf(ib-1) = dxf(ib)

    dxh(ib) = 2*xf(ib)
    do i=ib+1,ie+ih
       dxh(i) = xf(i) - xf(i-1)
    end do

    do k=kb,ke+kh
       do i=ib-ih,ie+ih
          delta(i,k) = (dxf(i)*dy*dzf(k))**(1./3.)
       end do
    end do

    !--------------------------------------------------
    ! *** Check whether the grid is equidistant *****
    !--------------------------------------------------

    !    leq=.true.
    leq=.false.          ! grid is now always non-equidistant
    dz = dzf(kb)
    do k=kb,ke+kh
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

    dzhi    = 1./dzh
    dzfi    = 1./dzf
    dzf2    = dzf*dzf
    dxhi    = 1./dxh
    dxfi    = 1./dxf
    dxf2    = dxf*dxf
    dyi     = 1./dy
    dy2     = dy*dy

    dzhiq   = 0.25*dzhi
    dzfiq   = 0.25*dzfi
    dxhiq   = 0.25*dxhi
    dxfiq   = 0.25*dxfi
    dyiq    = 0.25*dyi

    dzh2i   = dzhi*dzhi
    dxh2i   = dxhi*dxhi
    dy2i    = dyi*dyi

    dzfi5   = 0.5*dzfi
    dxfi5   = 0.5*dxfi
    dyi5    = 0.5*dyi

    ! Grid used in kappa scheme advection (extra ghost nodes)
    if (any(iadv_sv(1:nsv)==iadv_kappa)) then  
       allocate(dzfc(kb-khc:ke+khc))
       allocate(dxfc(ib-ihc:ie+ihc))
       allocate(dzfci(kb-khc:ke+khc))
       allocate(dxfci(ib-ihc:ie+ihc))
       allocate(dzhci(kb-1:ke+khc))
       allocate(dxhci(ib-1:ie+ihc))

       dzfc(kb-kh:ke+kh) = dzf(kb-kh:ke+kh)
       dzfc(kb-khc) = dzfc(kb-kh)
       dzfc(ke+khc) = dzfc(ke+kh)

       dxfc(ib-ih:ie+ih) = dxf(ib-ih:ie+ih)
       dxfc(ib-ihc) = dxfc(ib-ih)
       dxfc(ie+ihc) = dxfc(ie+ih)

       dzhci(kb:ke+kh)   = dzhi(kb:ke+kh)
       dzhci(kb-1)       = dzhci(kb)
       dzhci(ke+khc)     = dzhci(ke+kh)

       dxhci(ib:ie+ih)   = dxhi(ib:ie+ih)
       dxhci(ib-1)       = dxhci(ib)
       dxhci(ie+ihc)     = dxhci(ie+ih) 

       dzfci = 1./dzfc
       dxfci = 1./dxfc
    end if


    if(myid==0)then

       write (6,*) 'lev    dz     zf      zh       dzh    delta(ib,k)'
       do k=ke+1,kb,-1
          write(6,'(i4,5f8.5)') k,dzf(k),zf(k),zh(k),dzh(k),delta(ib,k)
       end do
       ! same for x:
       write (6,*) 'lev    dxf     xf      xh       dxh    delta(i,kb)'
       do i=ie+1,ib,-1
          write(6,'(i4,5f9.5)') i,dxf(i),xf(i),xh(i),dxh(i),delta(i,kb)
       end do
    end if
    tnextrestart = trestart
    tnextrestart2 = trestart2
    tnextfielddump = tfielddump
    timeleft     = btime+runtime


  end subroutine initglobal

  !> Clean up when leaving the run
  subroutine exitglobal
    deallocate(dsv,dzf,dzh,zh,zf,delta)
  end subroutine exitglobal

end module modglobal
