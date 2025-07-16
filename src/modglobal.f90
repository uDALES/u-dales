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

   integer :: poisrcheck = 0 ! switch to check if it is the first (RK) time step
   ! Simulation dimensions (parconst.f90)
   integer :: itot = 96 ! Used to be called imax
   integer :: jtot = 96
   integer :: ktot = 96 ! Rename to ktot?
   integer :: imax
   integer :: imax1
   integer :: imax2
   integer :: isen
   integer :: jmax
   integer :: jmax1
   integer :: jmax2
   integer :: jsen
   integer :: kmax
   integer :: kmax1
   integer :: kmax2
   integer ::  ib
   integer ::  ie
   integer ::  jb
   integer ::  je
   integer ::  jgb ! global j range - remove eventually
   integer ::  jge ! global j range - remove eventually
   integer ::  offset
   integer ::  kb
   integer ::  ke
   integer ::  nsv = 0 !< Number of additional scalar fields
   integer ::  nvar = 0
   character(50) :: fieldvars = ''

   integer ::  ih = 3
   integer ::  jh = 3
   integer ::  kh = 1
   integer ::  ihc = 2 ! used in k-scheme
   integer ::  jhc = 2 ! used in k-scheme
   integer ::  khc = 2 ! used in k-scheme

   integer :: nblocks = 0 ! no. of blocks in IBM
   integer, allocatable :: block(:,:)
   integer :: nfcts = -1 ! no. of wall facets
   integer ::  iplane ! ib+iplane is the plane that is stored when lstoreplane=.true.
   integer ::  nstore = 1002 ! number of rk steps in inletfile. This should be a multiple of three!
   character(90) :: fname_options = 'namoptions'
   integer, parameter :: longint = 8
   logical :: lwarmstart = .false. !<   flag for "cold" or "warm" start
   logical :: lstratstart = .false.
   logical :: lfielddump = .false. !< switch to enable the fielddump
   logical :: lreadscal = .false. !<   flag for reading scalar pollutant field (warm start)

   !Switches for boundary conditions
   !momentum (m), temperature (T), humidity (q) and scalars (s)
   !lateral in x/i direction (x), in y/j direction (y) at the top (top) and at the bottom (bot)

   !> x direction
   ! momentum
   integer, parameter :: BCxm_periodic = 1
   integer, parameter :: BCxm_profile = 2
   integer, parameter :: BCxm_driver = 3
   ! temperature
   integer, parameter :: BCxT_periodic = 1
   integer, parameter :: BCxT_profile = 2
   integer, parameter :: BCxT_driver = 3
   ! moisture
   integer, parameter :: BCxq_periodic = 1
   integer, parameter :: BCxq_profile = 2
   integer, parameter :: BCxq_driver = 3
   ! scalars
   integer, parameter :: BCxs_periodic = 1
   integer, parameter :: BCxs_profile = 2
   integer, parameter :: BCxs_driver = 3
   integer, parameter :: BCxs_custom = 4 ! Used to demonstrate flow
   ! set defaults
   integer :: BCxm = BCxm_periodic
   integer :: BCxT = BCxT_periodic
   integer :: BCxq = BCxq_periodic
   integer :: BCxs = BCxs_periodic

   !> y direction
   ! momentum
   integer, parameter :: BCym_periodic = 1
   integer, parameter :: BCym_profile = 2
   ! temperature
   integer, parameter :: BCyT_periodic = 1
   integer, parameter :: BCyT_profile = 2
   ! moisture
   integer, parameter :: BCyq_periodic = 1
   integer, parameter :: BCyq_profile = 2
   ! scalars
   integer, parameter :: BCys_periodic = 1
   integer, parameter :: BCys_profile = 2
   ! set defaults
   integer :: BCym = BCym_periodic
   integer :: BCyT = BCyT_periodic
   integer :: BCyq = BCyq_periodic
   integer :: BCys = BCys_periodic

   !> top
   ! momentum
   integer, parameter :: BCtopm_freeslip = 1 ! zero flux
   integer, parameter :: BCtopm_noslip = 2   ! fixed velocity
   integer, parameter :: BCtopm_pressure = 3 ! vertical velocity can vary according to pressure gradient
   ! temperature
   integer, parameter :: BCtopT_flux = 1  ! determined by flux wttop
   integer, parameter :: BCtopT_value = 2 ! determined by value thl_top
   ! moisture
   integer, parameter :: BCtopq_flux = 1  ! determined by flux wqtop
   integer, parameter :: BCtopq_value = 2 ! determined by value qt_top
   ! scalars
   integer, parameter :: BCtops_flux = 1  ! determined by flux wstop
   integer, parameter :: BCtops_value = 2 ! determined by value sv_top
   ! set defaults
   integer :: BCtopm = BCtopm_freeslip
   integer :: BCtopT = BCtopT_flux
   integer :: BCtopq = BCtopq_flux
   integer :: BCtops = BCtops_flux

   !> bottom
   ! tangential velocities (vertical is always impermeable)
   integer, parameter :: BCbotm_freeslip = 1  ! do nothing
   integer, parameter :: BCbotm_wf = 2        ! wall function
   integer, parameter :: BCbotm_wfneutral = 3 ! neutral wall function
   ! temperature
   integer, parameter :: BCbotT_flux = 1       ! determined by wtsurf
   integer, parameter :: BCbotT_wf = 2         ! wall function
   ! moisture
   integer, parameter :: BCbotq_flux = 1       ! determined by wtsurf
   ! scalars
   integer, parameter :: BCbots_flux = 1       ! zero flux
   ! set defaults
   integer :: BCbotm = BCbotm_wf
   integer :: BCbotT = BCbotT_flux
   integer :: BCbotq = BCbotq_flux
   integer :: BCbots = BCbots_flux

   integer :: BCzp = 1 ! 1: solve poisson equation using GE. 2: solve using cosine transform
   real :: ds = 0 ! Shifted boundary conditions

   integer :: iinletgen = 0 !<  0: no inletgen, 1: turb. inlet generator (Lund (1998)), 2: read inlet from file
   integer :: idriver = 0 !<  0: no inlet driver store, 1: Save inlet driver data, 2: read inlet driver data from file
   logical :: linoutflow = .false. !<  switch for periodic BC in both horizontal directions (false) or inflow/outflow in i and periodic in j.
   logical :: lzerogradtop = .false. !<  switch for zero gradient BC's at top wall (iinletgen 1 and 2 are seperate).
   logical :: lzerogradtopscal = .false. !
   logical :: lbuoyancy = .false. !<  switch for buoyancy force in modforces
   logical :: ltempeq = .false. !<  switch for solving temperature equation (either with or without buoyancy term)
   logical :: lscalrec = .false. !<
   logical :: lSIRANEinout = .false. !<
   logical :: ltempinout = .false. !<  seperate switch for inflow/outflow BC for temperature (only necessary when linoutflow.eqv..false.).
   logical :: lmoistinout = .false. !<  seperate switch for inflow/outflow BC for moisture (only necessary when linoutflow.eqv..false.).
   logical :: lper2inout = .false. !<  switch that determines type of restart: .true. means switching from periodic to in/outflow: inlet profile is read from prof.inp
   logical :: libm = .true. !<  switch that determines whether the Immersed Boundary Method is turned on
   logical :: lwalldist = .false. !<  switch that determines whether the wall distances should be computed
   logical :: lles = .true. !<  switch that determines whether the subgrid model is turned on or constant ekm and ekh are used (DNS)
   logical :: linletRA = .false. !<  switch that determines whether a Running Average should be used (.true.) in inlet generator
   logical :: lfixinlet = .false. !<  switch that determines whether the average inlet profiles can evolve or not (only used when iinletgen=1,2)
   logical :: lfixutauin = .false. !<  switch that determines whether the utau is kept fixed at the inlet (only used when iinletgen=1,2)
   logical :: lscasrc = .false. !
   logical :: lscasrcl = .false. !tg3315
   logical :: lydump = .false.  !<  switch to output y-averaged statistics every tsample
   logical :: lytdump = .false.  !<  switch to output y- and time- averaged statistics every tstatsdump
   logical :: lxydump    = .false.  !<  switch to output x- and y-avewraged statistics every tsample
   logical :: lxytdump   = .false.  !<  switch to output x-, y- and time-averaged statistics every tstatsdump
   logical :: lscasrcr  = .false.  !<  switch for network of point sources at lowest level
   logical :: ltkedump = .false. !tg3315
   logical :: lkslicedump= .false.  !<  switch to output slices in the xy-plane every tsample
   logical :: lislicedump= .false.  !<  switch to output slices in the yz-plane every tsample
   logical :: ljslicedump= .false.  !<  switch to output slices in the xz-plane every tsample
   integer :: kslice    = 1! k at which to output slice in xy-plane
   integer :: islice    = 1! i at which to output slice in yz-plane
   integer :: jslice    = 1! j at which to output slice in xz-plane
   integer :: isliceloc    ! local islice on core
   logical :: islicerank    ! cpu that islice is on
   integer :: jsliceloc    ! local jslice on core
   logical :: jslicerank    ! cpu that jslice is on
   logical :: ltdump    = .false.      !<  switch to output time-averaged statistics every tstatsdump
   logical :: lmintdump    = .false.      !<  switch to output prognostic statistics every tstatsdump

   logical :: ltrees = .false.         !<  switch to turn on trees module
   logical :: lpurif = .false.         !<  switch to turn on purifiers module
   logical :: ltreedump = .false.   !<  switch to output tree results time-averaged statistics every tstatsdump

   logical :: lreadminl = .false. !<  switch for reading mean inlet/recycle plane profiles (used in inletgenerator)
   logical :: lwallfunc = .true. !<  switch that determines whether wall functions are used to compute the wall-shear stress
   logical :: luoutflowr = .false. !<  switch that determines whether u-velocity is corrected to get a fixed outflow rate
   logical :: lvoutflowr = .false. !<  switch that determines whether u-velocity is corrected to get a fixed outflow rate
   logical :: luvolflowr = .false. !<  switch that determines whether u-velocity is corrected to get a fixed volume flow rate
   logical :: lvvolflowr = .false. !<  switch that determines whether u-velocity is corrected to get a fixed volume flow rate
   logical :: lstoreplane = .false. !<  switch that determines whether i-plane data is stored.
   logical :: lstorexy = .false. !xy files stored
   logical :: lreadmean = .false. !<  switch that determines whether mean variables should be read from means#myid#.#expnr#
   logical :: lstat = .false.
   logical :: lEB = .false.
   logical :: lwriteEBfiles = .false.
   logical :: lwritefac = .false.
   real    :: dtfac = 10.
   real    :: tfac = 0. !time of last calculation of facet quantites
   real    :: tnextfac = 0. !time for next calculation of facet energy balance
   logical :: lperiodicEBcorr = .false. ! Switch used to correct periodic heat build up.
   integer :: sinkbase = 0 ! This is the z index above which a sink is applied in periodicEBcorr scheme
   real    :: fraction = 1 ! Fraction of excess heat removed by volume sink in periodic energy balance correction.

   logical :: lvfsparse = .false. !< whether to read in view factors in sparse format
   integer :: nnz !< number of non-zero view factors
   logical :: lconstW = .false.  ! The evaporated water can be removed from the soil (lconstW=false) or the soil moisture can be assumed as constant in time (lconstW=true)
   logical :: lfacTlyrs = .false.
!  logical :: ifixuinf   = .true. !dpdxl relaxed to have Uinf 1. dpdx = (1/dt)*(Uh-Uinf)2. d/dt(dpdx) = 1/tau*(Uh-Uinf)
   integer :: ifixuinf = 0
   logical :: lvinf = .false. !use Vinf instead of Uinf for the fixed velocity at infinity
   logical :: lrandomize = .true.

   logical :: ibrank
   logical :: ierank
   logical :: jbrank
   logical :: jerank

   real    :: freestreamav = 0. !
   real    :: freestrtmpav = 0. !
   !<  Global constants modconst.f90
   !< File numbers

   integer, parameter :: ifinput = 1
   integer, parameter :: ifoutput = 2
   integer, parameter :: ifnamopt = 3

   real, parameter :: pi = 3.141592653589793116
   real, parameter :: grav = 9.81 !<    *gravity acceleration.
   real, parameter :: rd = 287.04 !<    *gas constant for dry air.
   real, parameter :: rv = 461.5 !<    *gas constant for water vapor.
   real, parameter :: cp = 1004. !<    *specific heat at constant pressure (dry air).
   real, parameter :: rlv = 2.26e6 !<    *latent heat for vaporisation.
   real, parameter :: rlvi = 1/rlv !inverse
   real, parameter :: ep = rd/rv !<    0.622
   real, parameter :: ep2 = rv/rd - 1. !<    0.61
   !< real,parameter :: cv       = cp-rd            !<    716.96
   real, parameter :: rcp = rd/cp !<    0.286
   real, parameter :: cpr = cp/rd !<    3.50
   real, parameter :: rlvocp = rlv/cp !<    2.49
   real, parameter :: mair = 28.967 !< Molar mass of air
   real, parameter :: rhoa = 1.2 !density of air used in some calculations

   real :: wfc = 313. !water content at field capacity (kg/m3)
   real :: wwilt = 171. !water ocntent at wilting point (kg/m3)
   real :: wgrmax = 450. !maximum water content (kg/m3)
   real :: rsmin = 110. !minimum resistance of soil/plant
   real :: rsmax = 5000. !maximum resistance of soil/plant
   real :: GRLAI = 2. !Leave area index of green roof
   real :: wsoil = 0. !water content of soil (kg/m3)
   real :: bldT = 0. !building internal temperature, currently also ground temperature at a depth equal to floor facet thickness
   real :: flrT = 0. !ground internal temperature
   real :: skyLW = 0. !longwave radiation from the sky
   real :: gres = 0. !saturation vapour pressure of green roof
   real :: grqs = 0. !saturation humidity of green roof
   real :: grdqdt = 0. !gradient of saturation humidity for green roof

   real, parameter :: numol = 1.5e-5 !< kinematic viscosity for couette flow Re=5000 (Re=Uinf*H/(2*nu)) H=1, Uinf=1
   real, parameter :: numoli = 1./numol !< 1/numol
   real, parameter :: prandtlmol = 0.71 !< Prandtl number (for air at 300K). Fluid property!
   real, parameter :: prandtlmoli = 1./prandtlmol !< Inverse of Prandtl number
   real :: prandtlturb = prandtlmol

   integer         :: iwallmom = 2, iwalltemp = 1, iwallmoist = 1, iwallscal = 1

   real, parameter :: rhow = 0.998e3 !<    * Density of water
   real, parameter :: pref0 = 1.e5 !<    *standard pressure used in exner function.
   real, parameter :: tmelt = 273.16 !<    *temperature of melting of ice.
   real, parameter :: es0 = 610.78 !<    * constants used for computation
   real, parameter :: at = 17.27 !<    * of saturation mixing ratio
   real, parameter :: bt = 35.86 !<    * using Tetens Formula.
   !      real,parameter :: ekmin    = 1.e-6            !<    *minimum value for k-coefficient.
   real, parameter :: ekmin = 1.e-12 !<    *minimum value for k-coefficient.
   real, parameter :: e12min = 5.e-5 !<    *minimum value for TKE.
   real :: fkar = 0.41 !<   *Von Karman constant
   real, parameter :: eps1 = 1.e-10 !<    *very small number*
   real, parameter :: epscloud = 1.e-5 !<    *limit for cloud calculation 0.01 g/kg
   real, parameter :: boltz = 5.67e-8 !<    *Stefan-Boltzmann constant

   real, parameter, dimension(3) :: xhat = (/1.,0.,0./)
   real, parameter, dimension(3) :: yhat = (/0.,1.,0./)
   real, parameter, dimension(3) :: zhat = (/0.,0.,1./)
   real, parameter, dimension(3) :: vec0 = (/0.,0.,0./) ! zero vector

   logical :: lprofforc = .false. !<  nudge flow to a profile !
   logical :: lcoriol = .false. !<  switch for coriolis force
   integer :: igrw_damp = 0 !< switch to enable gravity wave damping
   real    :: geodamptime = 7200. !< time scale for nudging to geowind in sponge layer, prevents oscillations
   real    :: uflowrate = 1. !< fixed flow rate used for u-velocity correction
   real    :: vflowrate = 1. !< fixed flow rate used for v-velocity correction
   real    :: Uinf = 0. !< fixed U_inf (used in inlet generator), also in conjunction with ifixuinf
   real    :: Vinf = 0. !fixed V_inf
   real    :: inletav = 0. !< averaging interval for inlet generator
   real    :: totinletav = 0. !< averaging interval for inlet generator (used in Running Average)
   real    :: om22 !<    *2.*omega_earth*cos(lat)
   real    :: om23 !<    *2.*omega_earth*sin(lat)
   real    :: om22_gs !<    *2.*omega_earth*cos(lat)
   real    :: om23_gs !<    *2.*omega_earth*sin(lat)
   real    :: xlat = 52. !<    *latitude  in degrees.
   real    :: xlon = 0. !<    *longitude in degrees.

   !scalar source in fluid domain
   real, allocatable :: xSa(:)
   real, allocatable :: ySa(:)
   real, allocatable :: zSa(:)
   real    :: xS = 0., yS = 0., zS = 0.
   real    :: xSb = 0., ySb = 0., zSb = 0.
   real    :: xSe = 0., ySe = 0., zSe = 0.
   real    :: SS = 0.
   real    :: sigS = 0.
   integer :: nscasrc = 0              !< number of scalar point sources
   integer :: nscasrcl = 0              !< number of scalar line sources
   real, allocatable :: scasrcp(:,:,:)    !< field with data from scalarsourcep.inp.xxx containing coordinates of the source points, strength and standard deviation
   real, allocatable :: scasrcl(:,:,:)    !< field with data from scalarsourcel.inp.xxx containing coordinates of the end points of line sources, strength per unit length and standard deviation

   !trees
   integer, allocatable :: tree(:,:)             !< field with data from tree.inp.xxx
   integer :: ntree_max = 0
   integer :: ntrees = 0
   !real, allocatable :: ladz(:)                  !< field with leaf area density data
   real    :: cd = 0., ud = 0., Qstar = 0., dQdt = 0., dec = 0., lad = 0., lsize = 0., r_s = 0.  !< current set of tree parameters
            ! volumetric drag coefficient, deposition velocity, net radiation, dQ*/dt , extinction coefficient, leaf area density, characteristic leaf size, stomatal resistance, respectively
   real    :: tr_A = 0.

   logical :: lnudge = .false.                   !< switch for applying nudging at the top of the domain
   logical :: lnudgevel = .true.                 !< switch for nudging velocities
   real    :: tnudge = 60.                       !< time scale for nudging
   integer :: nnudge = 0                         !< number of points from kb to start nudging

   !chemistry
   logical :: lchem = .false.    ! switch for basic chemistry
   real    :: k1 = 0., JNO2 = 0.   ! k1 = rate constant (O3 + NO -> NO2 + 02 ), JNO2 = NO2 photolysis rate

   !purifiers
   integer, allocatable :: purif(:,:)            !< field with data from purif.inp.xxx
   integer :: npurif = 0
   real    :: Qpu = 0., epu = 0.                 !< flowrate and efficiency of purifiers

   ! Heat pump
   logical :: lheatpump = .false.   !< switch for heat pump
   logical :: lfan_hp = .true.      !< switch for heat pump fan on/off
   integer :: nhppoints = 0         !< number of heat pump grid points
   real    :: QH_dot_hp = 0.        !< Total rate of heat extracted from the ambient air by the heat pump (W)
   real    :: Q_dot_hp = 0.         !< Total volume flux going out from the heat pump in vertical direction(m^3/s)

   ! Poisson solver
   integer, parameter :: POISS_FFT2D = 0, &
                         POISS_CYC   = 1, &
                         POISS_FFT3D = 2, &
                         POISS_FFT2D_2DECOMP = 3

   integer :: ipoiss   = POISS_FFT2D

   !Advection scheme
   integer, parameter :: iadv_upw = 1  !< first order upwind scheme
   integer, parameter :: iadv_cd2 = 2  !< second order central difference scheme
   integer, parameter :: iadv_kappa = 7  !< Kappa scheme
   integer :: iadv_mom = 2, iadv_tke = -1, iadv_thl = -1, iadv_qt = -1, iadv_sv(100) = -1

   logical :: lmoist = .false. !<   switch to calculate moisture fields
   ! Global variables (modvar.f90)
   real :: xday = 1. !<     * day number
   real :: xtime = 0. !<     * GMT time
   real :: runtime = 300. !<     * simulation time in secs
   real :: dtmax = 20. !<     * maximum time integration interval

   real    :: trestart = 10000. !<     * each trestart sec. a restart file is written to disk. bss116: per default do not write restart files
   real    :: tfielddump = 10000. !< Time step for field outputs
   real    :: tsample = 5. !<    Sample time steps for statistics
   real    :: tstatsdump = 10000. !< Time step for statistics outputs tg3315
   real    :: tstatstart = 0.      !< Starting time of statistic computation
   real    :: tnextrestart !<     * each trestart sec. a restart file is written to disk
   real    :: tscale !       timescale: domain height*Uinf/utau**2
   real    :: tnextfielddump !<
   character(90) :: startfile = '' !<    * name of the restart file

   real :: totavtime = 0. !<    * the total time over which the values are averaged in meansXXX.XXX
   real :: dtEB = 10. !time interval between calculations of facet energy balance
   real :: tEB = 0. !time of last calculation of facet energy balance
   real :: tnextEB = 0. !time for next calculation of facet energy balance
   real :: totheatflux = 0. ! Total sensible heat flux from facs into air in one timestep
   real :: totqflux  = 0. ! Total latent heat flux from facs into air in one timestep


   real :: thres = 5.e-3 !<     * threshold value for inversion height calculations
   real :: dqt !<     * applied gradient of qt at top of model
   real :: dtheta !<     * applied gradient of theta at top of model
   real, allocatable :: dsv(:) !<     * applied gradient of sv(n) at top of model

   real ::  dt !<     * time integration interval
   !      integer(kind=longint) :: timee             !<     * elapsed time since the "cold" start
   real :: timee !<     * elapsed time since the "cold" start
   !      integer(kind=longint) :: btime             !<     * time of (re)start
   real :: btime !<     * time of (re)start
   real :: runavtime !<     * time of starting running average
   integer :: ntimee !<     * number of timesteps since the cold start
   integer :: ntrun !<     * number of timesteps since the start of the run
   real    :: timeleft
   logical :: ladaptive = .false. !<    * adaptive timestepping on or off

   real    :: tdriverstart = 0.   !<     * time at which to start recording inlet driver file (only necessary if idriver == 1)
   real    :: tdriverstart_cold = 0.   !< to store tdriverstart of cold started simulation while doing warmstart
   real    :: tdriverdump         !<     * time in inlet driver simulation at which data dumps are made (idriver == 1)
   real    :: dtdriver = 0.1      !<     * time frequency at which inlet driver data dumps are made (idriver == 1)
   integer :: driverstore         !<     * number of stored driver steps for inlet (automatically calculated)
   integer :: driverjobnr         !<     * Job number of the driver inlet generation run (idriver == 2)
   character(3) :: cdriverjobnr
   logical :: lhdriver = .false.    !<     * switch for reading temperature driver files
   logical :: lqdriver = .false.    !<     * switch for reading temperature driver files
   logical :: lsdriver = .false.   !<     * switch for reading scalar driver files
   logical :: iplanerank = .false.
   integer :: driverid
   character(3) :: cdriverid

   logical :: lchunkread = .false.     !< * logical switch for chunkwise reading of driver files
   integer :: chunkread_size = 100      !< * chunk size of each reading

   real    :: courant = -1.
   real    :: diffnr = 0.25
   real    :: dt_lim

   integer :: rk3step = 0

   integer :: iexpnr = 0 !<     * number of the experiment

   character(3) cexpnr

   real :: thlsrc = 0.

   ! modphsgrd.f90

   real :: dx !<  grid spacing in x-direction
   real :: dx2 !<  grid spacing in x-direction squared
   real :: dxi !<  1/dx
   real :: dxiq !<  1/(dx*4)
   real :: dxi5 !<  1/(dx*2)
   real :: dx2i !<  (1/dx)**2
   real :: dy !<  grid spacing in y-direction
   real :: dy2 !<  grid spacing in y-direction squared
   real :: dz !<  grid spacing in z-direction
   real :: dyi !<  1/dy
   real :: dyiq !<  1/(dy*4)
   real :: dyi5 !<  1/(dy*2)
   real :: dy2i !<  (1/dy)**2

   integer :: nfaclyrs = 3
   real, allocatable   :: AM(:,:), BM(:,:), CM(:,:), DM(:,:), EM(:,:), FM(:,:), GM(:,:), HM(:,:), inAM(:,:),IDM(:,:) !matrices for the facet energy balance
   real, allocatable   :: bb(:),w(:),dumv(:),Tdash(:) !vector for the facet energy balance

   real :: rslabs
   real, allocatable :: dzf(:) !<  thickness of full level
   real, allocatable :: dzfc(:) !<  thickness of full level (extra ghost nodes (used in k-scheme)
   real, allocatable :: dzfci(:) !<  1/dzfc
   real, allocatable :: dzf2(:) !<  thickness of full level squared
   real, allocatable :: dzh(:) !<  thickness of half level
   real, allocatable :: zh(:) !<  height of half level [m]
   real, allocatable :: zf(:) !<  height of full level [m]
   real, allocatable :: dzfi(:) !<  1/dzf
   real, allocatable :: dzfiq(:) !<  0.25*(1/dzf)
   real, allocatable :: dzfi5(:) !<  0.5*(1/dzf)
   real, allocatable :: dzhi(:) !<  1/dzh
   real, allocatable :: dzhci(:) !<  1/dzh (extra ghost nodes (used in k-scheme)
   real, allocatable :: dzhiq(:) !<  0.25*(1/dzh)
   real, allocatable :: dzh2i(:) !<  1/dzh^2
   real, allocatable :: zhi(:) !<  1/zh
   real, allocatable :: zfi(:) !<  1/zf
   real, allocatable :: dxf(:) !<  thickness of full level
   real, allocatable :: dxfc(:) !<  thickness of full level (extra ghost nodes (used in k-scheme)
   real, allocatable :: dxfci(:) !<  1/dxfc
   real, allocatable :: dxf2(:) !<  thickness of full level squared
   real, allocatable :: dxfi(:) !<  = 1/dxf
   real, allocatable :: dxfiq(:) !<  = 0.25*(1/dxf)
   real, allocatable :: dxfi5(:) !<  = 0.5*(1/dxf)
   real, allocatable :: dxh(:) !<  thickness of half level
   real, allocatable :: dxhi(:) !<  = 1/dxh
   real, allocatable :: dxhci(:) !<  = 1/dxh (with extra ghost nodes (used in k-scheme))
   real, allocatable :: dxhiq(:) !<  = 0.25*(1/dxh)
   real, allocatable :: dxh2i(:) !<  = 1/dxh^2
   real, allocatable :: xh(:) !<  height of half level [m]
   real, allocatable :: xf(:) !<  height of full level [m]
   real, allocatable :: yh(:) !<  height of half level [m]
   real, allocatable :: yf(:) !<  height of full level [m]
   real :: xlen = -1. !<  domain size in x-direction
   real :: ylen = -1. !<  domain size in y-direction
   real, allocatable :: delta(:, :) !<  (dx*dy*dz)**(1/3)

   logical :: lmomsubs = .false. !<  switch to apply subsidence on the momentum or not
   character(80) :: author = '', version = 'DALES U'
contains

   !> Initialize global settings.
   !!
   !! Set courant number, calculate the grid sizes (both computational and physical), and set the coriolis parameter
   subroutine initglobal
      use modmpi,   only : myid, comm3d, my_real, mpierr
      use decomp_2d
      implicit none

      integer :: advarr(4)
      real phi, colat, silat, omega, omega_gs
      integer :: i, j, k, n
      character(80) chmess

      ! Global constants

      ! Select advection scheme for scalars. If not set in the options file, the momentum scheme is used
      if (iadv_tke < 0) iadv_tke = iadv_mom
      if (iadv_thl < 0) iadv_thl = iadv_mom
      if (iadv_qt < 0) iadv_qt = iadv_mom

      !CvH remove where
      !where (iadv_sv<0)  iadv_sv  = iadv_mom

      !tg3315 added - only uses kappa advection scheme...
      do n = 1, nsv
         iadv_sv(n) = iadv_kappa
      end do
      !ends here

      !timestepping
      if (courant < 0) then
         select case (iadv_mom)
         case (iadv_cd2)
            courant = 1.5
         case default
            courant = 1.4
         end select
         if (any(iadv_sv(1:nsv) == iadv_kappa) .or. any((/iadv_thl, iadv_qt, iadv_tke/) == iadv_kappa)) then
            courant = min(courant, 1.1)
         elseif (any(iadv_sv(1:nsv) == iadv_upw) .or. any((/iadv_thl, iadv_qt, iadv_tke/) == iadv_upw)) then
            courant = min(courant, 1.1)
         elseif (any(iadv_sv(1:nsv) == iadv_cd2) .or. any((/iadv_thl, iadv_qt, iadv_tke/) == iadv_cd2)) then
            courant = min(courant, 1.5)
         end if
      end if

      ! phsgrid

      !jmax = jtot/nprocy ! Only in z-pencil and not true generally - uneven no.
      !imax = itot/nprocx ! Only in z-pencil
      isen = imax ! Only in z-pencil - replace eventually so it is pencil-independent (in poisson)
      jsen = jmax ! Only in z-pencil - replace eventually so it is pencil-independent (in poisson)
      !set the number of ghost cells. NB: This switch has to run in order of required ghost cells
      advarr = (/iadv_mom, iadv_tke, iadv_thl, iadv_qt/)
      if (any(advarr == iadv_kappa)) then
         ih = 2
         jh = 2
         kh = 1
         ! SO: think this is inconsistent
      elseif (any(advarr == iadv_cd2) .or. any(iadv_sv(1:nsv) == iadv_cd2)) then
         ih = 1
         jh = 1
         kh = 1
         ihc = 1
         jhc = 1
         khc = 1
      end if

      ! J. Tomas added this for using only kappa scheme for sv(:)
      if (any(iadv_sv(1:nsv) == iadv_kappa) .or. (iadv_thl == iadv_kappa)) then
         ih = 1
         jh = 1
         kh = 1
         ihc = 2
         jhc = 2
         khc = 2
      end if

      ! Eventually ib etc should be completely replaced.
      ! All arrays start at 1, like in 2DECOMP, and end at e.g. zsize(1) = imax in old terminology
      ib = 1    ! Remove eventually
      jb = 1
      jgb = jb ! global j range (starting at the same as j as the processor j range)
      jge = jtot ! global j range
      kb = 1 ! Make redundant

      !kmax = ktot

      ! Define indices in terms of 2DECOMP's. Subject to change! z=pencil 'special' for now, but could rename e.g. imax -> imax3
      imax1 = xsize(1) !=itot
      imax2 = ysize(1)
      imax = zsize(1)

      jmax1 = xsize(2)
      jmax2 = ysize(2) !=jtot
      jmax = zsize(2)

      kmax1 = xsize(3)
      kmax2 = ysize(3)
      kmax = zsize(3)

      ie = imax
      je = jmax
      ke = kmax

      decomp_main%zlevel = (/ih, jh, kh/)

      if (zstart(1) == 1) then
        ibrank = .true.
      else
        ibrank = .false.
      end if

      if (zend(1) == itot) then
        ierank = .true.
      else
        ierank = .false.
      end if

      if (zstart(2) == 1) then
        jbrank = .true.
      else
        jbrank = .false.
      end if

      if (zend(2) == jtot) then
        jerank = .true.
      else
        jerank = .false.
      end if

      !write(*,*) "myid, ibrank, ierank", myid, ibrank, ierank

      phi = xlat*pi/180.
      colat = cos(phi)
      silat = sin(phi)

      omega = 7.292e-5
      omega_gs = 7.292e-5
      om22 = 2.*omega*colat
      om23 = 2.*omega*silat
      om22_gs = 2.*omega_gs*colat
      om23_gs = 2.*omega_gs*silat

      ! Variables
      allocate (dsv(nsv))

      ! Create the physical grid variables
      allocate (dzf(kb - kh:ke + kh))
      allocate (dzf2(kb - kh:ke + kh))
      allocate (dzfi(kb - kh:ke + kh))
      allocate (dzfiq(kb - kh:ke + kh))
      allocate (dzfi5(kb - kh:ke + kh))
      allocate (dzh(kb:ke + kh))
      allocate (dzhi(kb:ke + kh))
      allocate (dzhiq(kb:ke + kh))
      allocate (dzh2i(kb:ke + kh))
      allocate (zh(kb:ke + kh))
      allocate (zf(kb:ke + kh))

      allocate (dxf(ib-ih:itot+ih))
      allocate (dxf2(ib-ih:itot+ih))
      allocate (dxfi(ib-ih:itot+ih))
      allocate (dxfiq(ib-ih:itot+ih))
      allocate (dxfi5(ib-ih:itot+ih))
      allocate (dxh(ib:itot+ih))
      allocate (dxhi(ib:itot+ih))
      allocate (dxhiq(ib:itot+ih))
      allocate (dxh2i(ib:itot+ih))
      allocate (xh(ib:itot+ih))
      allocate (xf(ib:itot+ih))
      allocate (yh(jb:jtot+jh))
      allocate (yf(jb:jtot+jh))
      allocate (delta(ib-ih:itot+ih, kb:ke + kh))

      rslabs = real(itot*jtot)

      dx = xlen/float(itot)
      dy = ylen/float(jtot)

      ! MPI

      ! Note, that the loop for reading zf and calculating zh
      ! has been split so that reading is only done on PE 1

      write (cexpnr, '(i3.3)') iexpnr
      if (nrank == 0) then
         open (ifinput, file='prof.inp.'//cexpnr)
         read (ifinput, '(a72)') chmess
         read (ifinput, '(a72)') chmess

         do k = kb, ke
            read (ifinput, *) zf(k)
         end do
         close (ifinput)

         ! ! J. Tomas: Read the x-coordinates of the cell centers from xgrid.inp.XXX
         ! ! SO: still reads for now, but need to remove any reference to xf, xh, etc eventually
         ! open (ifinput, file='xgrid.inp.'//cexpnr)
         ! read (ifinput, '(a72)') chmess
         ! read (ifinput, '(a72)') chmess
         !
         ! do i = ib, itot
         !    read (ifinput, *) xf(i)
         ! end do
         ! close (ifinput)

      end if ! end if nrank==0

      ! MPI broadcast ktot elements from zf
      call MPI_BCAST(zf, ktot, MY_REAL, 0, comm3d, mpierr)
      ! MPI broadcast itot elements from xf
      ! call MPI_BCAST(xf, itot, MY_REAL, 0, comm3d, mpierr)

      zh(kb) = 0.0
      do k = kb, ke
         zh(k + 1) = zh(k) + 2.0*(zf(k) - zh(k))
      end do
      zf(ke + kh) = zf(ke) + 2.0*(zh(ke + kh) - zf(ke))

      do k = kb, ke
         dzf(k) = zh(k + 1) - zh(k)
      end do
      dzf(ke + 1) = dzf(ke)
      dzf(kb - 1) = dzf(kb)

      dzh(kb) = 2*zf(kb)
      do k = kb + 1, ke + kh
         dzh(k) = zf(k) - zf(k - 1)
      end do

      ! j. tomas: same trick for x-direction...
      ! xh(ib) = 0.0
      ! do i = ib, itot
      !    xh(i + 1) = xh(i) + 2.0*(xf(i) - xh(i))
      ! end do
      ! xf(itot + ih) = xf(itot) + 2.0*(xh(itot + ih) - xf(itot))

      do i=ib,itot+ih
        xh(i) = (i-1) * dx
        xf(i) = xh(i) + dx/2
      end do

      do j=jb,jtot+jh
        yh(j) = (j-1) * dy
        yf(j) = yh(j) + dy/2
      end do

      ! These should be removed eventually
      do i = ib, itot
         dxf(i) = xh(i + 1) - xh(i)
      end do
      dxf(itot + 1) = dxf(itot)
      dxf(ib - 1) = dxf(ib)

      dxh(ib) = 2*xf(ib)
      do i = 2, itot + ih
         dxh(i) = xf(i) - xf(i - 1)
      end do

      do k = kb, ke + kh
         do i = ib - ih, itot + ih
            delta(i, k) = (dxf(i)*dy*dzf(k))**(1./3.)
         end do
      end do

      !--------------------------------------------------
      ! *** Check whether the grid is equidistant *****
      !--------------------------------------------------

      !if (myid == 0) then
      !do k=kb,ke+kh
      !if (.not.(dzf(k).eq.dzf(1)))
      !      write (6, *) &
      !      'WARNING, You are working with a non-equidistant grid!!!!'
      !end if
      !end do
      !end if ! end if myid==0

      dzhi = 1./dzh
      dzfi = 1./dzf
      dzf2 = dzf*dzf
      dxhi = 1./dxh
      dxfi = 1./dxf
      dxf2 = dxf*dxf
      dxi = 1./dx
      dx2 = dx*dx
      dyi = 1./dy
      dy2 = dy*dy

      dzhiq = 0.25*dzhi
      dzfiq = 0.25*dzfi
      dxhiq = 0.25*dxhi
      dxfiq = 0.25*dxfi
      dyiq = 0.25*dyi
      dxiq = 0.25*dxi

      dzh2i = dzhi*dzhi
      dxh2i = dxhi*dxhi
      dy2i = dyi*dyi
      dx2i = dxi*dxi

      dzfi5 = 0.5*dzfi
      dxfi5 = 0.5*dxfi
      dyi5 = 0.5*dyi
      dxi5 = 0.5*dxi

      ! Grid used in kappa scheme advection (extra ghost nodes)
      if (any(iadv_sv(1:nsv) == iadv_kappa) .or. (iadv_thl == iadv_kappa)) then
         allocate (dzfc(kb - khc:ke + khc))
         allocate (dxfc(ib - ihc:itot + ihc))
         allocate (dzfci(kb - khc:ke + khc))
         allocate (dxfci(ib - ihc:itot + ihc))
         allocate (dzhci(kb - 1:ke + khc))
         allocate (dxhci(ib - 1:itot + ihc))

         dzfc(kb - kh:ke + kh) = dzf(kb - kh:ke + kh)
         dzfc(kb - khc) = dzfc(kb - kh)
         dzfc(ke + khc) = dzfc(ke + kh)

         dxfc(ib - ih:itot + ih) = dxf(ib - ih:itot + ih)
         dxfc(ib - ihc) = dxfc(ib - ih)
         dxfc(itot + ihc) = dxfc(itot + ih)

         dzhci(kb:ke + kh) = dzhi(kb:ke + kh)
         dzhci(kb - 1) = dzhci(kb)
         dzhci(ke + khc) = dzhci(ke + kh)

         dxhci(ib:itot + ih) = dxhi(ib:itot + ih)
         dxhci(ib - 1) = dxhci(ib)
         dxhci(itot + ihc) = dxhci(itot + ih)

         dzfci = 1./dzfc
         dxfci = 1./dxfc
      end if

      tnextrestart = trestart
      tnextfielddump = tfielddump
!    tnextstatsdump = tstatsdump
      timeleft = runtime ! tg3315 previously btime + runtime

   end subroutine initglobal

   !> Clean up when leaving the run
   subroutine exitglobal
      deallocate (dsv, dzf, dzh, zh, zf, delta)
   end subroutine exitglobal

end module modglobal
