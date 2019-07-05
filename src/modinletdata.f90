!!
!!  \author Jasper Tomas,TU Delft, 31 March 2014
!!  \par Revision list
!!  \todo Documentation
!

module modinletdata
  implicit none
  save

    real, allocatable :: storeu0inletbc(:,:,:)
    real, allocatable :: storev0inletbc(:,:,:)
    real, allocatable :: storew0inletbc(:,:,:)
    real, allocatable :: storet0inletbc(:,:,:)
    real, allocatable :: u0rot(:,:,:)    ! dummy variable to eventually rotate u0 in horizontal plane
    real, allocatable :: v0rot(:,:,:)    ! dummy variable to eventually rotate v0 in horizontal plane

    real, allocatable :: Utav(:,:)       !< j-averaged and time averaged mean velocity
    real, allocatable :: QLtav(:,:)       !< j-averaged and time averaged mean velocity
    real, allocatable :: QTtav(:,:)       !< j-averaged and time averaged mean velocity
    real, allocatable :: Ttav(:,:)       !< j-averaged and time averaged mean temperature
    real, allocatable :: uaver(:,:)     !< j-averaged  mean velocity
    real, allocatable :: taver(:,:)     !< j-averaged  mean temperature
    real, allocatable :: u0inletbc(:,:)  !< final computed inlet bc for u
    real, allocatable :: v0inletbc(:,:)  !< final computed inlet bc for v
    real, allocatable :: w0inletbc(:,:)  !< final computed inlet bc for t
    real, allocatable :: t0inletbc(:,:)  !< final computed inlet bc for w
    real, allocatable :: u0inletbcold(:,:)  !< old final computed inlet bc for u
    real, allocatable :: v0inletbcold(:,:)  !< old final computed inlet bc for v
    real, allocatable :: w0inletbcold(:,:)  !< old final computed inlet bc for w
    real, allocatable :: t0inletbcold(:,:)  !< old final computed inlet bc for t
    real, allocatable :: uminletbc(:,:)  !< final computed inlet bc for u
    real, allocatable :: vminletbc(:,:)  !< final computed inlet bc for v
    real, allocatable :: wminletbc(:,:)  !< final computed inlet bc for w
    real, allocatable :: tminletbc(:,:)  !< final computed inlet bc for t
    real, allocatable :: Uinl(:)   !< mean inlet u-velocity (time-averaged and y-averaged)
    real, allocatable :: QLinl(:)   !< mean inlet u-velocity (time-averaged and y-averaged)
    real, allocatable :: QTinl(:)   !< mean inlet u-velocity (time-averaged and y-averaged)
    real, allocatable :: Winl(:)   !< mean inlet w-velocity  (time-averaged and y-averaged)
    real, allocatable :: Tinl(:)   !< mean inlet temperature  (time-averaged and y-averaged)
    real, allocatable :: Urec(:)   !< mean recycle u-velocity (time-averaged and y-averaged)
    real, allocatable :: QLrec(:)   !< mean recycle u-velocity (time-averaged and y-averaged)
    real, allocatable :: QTrec(:)   !< mean recycle u-velocity (time-averaged and y-averaged)
    real, allocatable :: Wrec(:)   !< mean recycle w-velocity  (time-averaged and y-averaged)
    real, allocatable :: Trec(:)   !< mean recycle temperature  (time-averaged and y-averaged)
    real, allocatable :: zirf(:)      !< zf-coordinate at recycle station (inner scaling)
    real, allocatable :: ziif(:)      !< zf-coordinate at inlet station (inner scaling)
    real, allocatable :: zirh(:)      !< zh-coordinate at recycle station (inner scaling)
    real, allocatable :: ziih(:)      !< zh-coordinate at inlet station (inner scaling)
    real, allocatable :: zorf(:)      !< zf-coordinate at recycle station (outer scaling)
    real, allocatable :: zoif(:)      !< zf-coordinate at inlet station (outer scaling)
    real, allocatable :: zorh(:)      !< zh-coordinate at recycle station (outer scaling)
    real, allocatable :: zoih(:)      !< zh-coordinate at inlet station (outer scaling)
    real, allocatable :: zotr(:)      !< zf-coordinate at recycle station (temperature outer scaling)
    real, allocatable :: zoti(:)      !< zf-coordinate at inlet station (temperature outer scaling)
    real, allocatable :: displ(:)     !< displacement thickness
    real, allocatable :: displold(:)  !< old displacement thickness
    real, allocatable :: upupavinl(:)   !< j-averaged time-averaged u'u' at the inlet
    real, allocatable :: vpvpavinl(:)   !< j-averaged time-averaged v'v' at the inlet
    real, allocatable :: wpwpavinl(:)   !< j-averaged time-averaged w'w' at the inlet
    real, allocatable :: upwpavinl(:)   !< j-averaged time-averaged u'w' at the inlet
    real, allocatable :: thlpthlpavinl(:) !< j-averaged time-averaged thl'thl' at the inlet
    real, allocatable :: thlpupavinl(:)   !< j-averaged time-averaged thl'u' at the inlet
    real, allocatable :: thlpwpavinl(:)   !< j-averaged time-averaged thl'w' at the inlet
    real, allocatable :: qlpqlpavinl(:) !< j-averaged time-averaged thl'thl' at the inlet
    real, allocatable :: qlpupavinl(:)   !< j-averaged time-averaged thl'u' at the inlet
    real, allocatable :: qlpwpavinl(:)   !< j-averaged time-averaged thl'w' at the inlet
    real, allocatable :: qtpqtpavinl(:) !< j-averaged time-averaged thl'thl' at the inlet
    real, allocatable :: qtpupavinl(:)   !< j-averaged time-averaged thl'u' at the inlet
    real, allocatable :: qtpwpavinl(:)   !< j-averaged time-averaged thl'w' at the inlet
    real, allocatable :: zfin(:)                       ! zf from inlet simulation 
    real, allocatable :: zhin(:)                       ! zh from inlet simulation 
    real, allocatable :: dzfin(:)                      ! dzf from inlet simulation 
    real, allocatable :: dzhin(:)                      ! dzh from inlet simulation 
    real, allocatable :: heavif(:)          !< Heaviside function for u,v 
    real, allocatable :: heavih(:)          !< Heaviside function for w 
    real, allocatable :: heavit(:)          !< Heaviside function for t 
    integer, allocatable :: loclowif(:)   !< index of lower zir at full level
    integer, allocatable :: locupif(:)    !< index of upper zir at full level
    integer, allocatable :: loclowih(:)   !< index of lower zir at half level
    integer, allocatable :: locupih(:)    !< index of upper zir at half level
    integer, allocatable :: loclowof(:)   !< index of lower zor at full level
    integer, allocatable :: locupof(:)    !< index of upper zor at full level
    integer, allocatable :: loclowoh(:)   !< index of lower zor at half level
    integer, allocatable :: locupoh(:)    !< index of upper zor at half level
    integer, allocatable :: loclowot(:)   !< index of lower zot at full level
    integer, allocatable :: locupot(:)    !< index of upper zot at full level
    integer, allocatable :: linlf(:)      !< index of lower zfin 
    integer, allocatable :: linuf(:)      !< index of upper zfin 
    integer, allocatable :: linlh(:)      !< index of lower zhin 
    integer, allocatable :: linuh(:)      !< index of upper zhin 
    
    real :: di=0.09        !< delta at inlet  (should be prescribed in namoptions!)  corresponds to critertion 0.99
!    real :: di=0.12        !< delta at inlet  (should be prescribed in namoptions!)  corresponds to critertion 1.0
!    real :: di=0.0645        !< delta at inlet  (should be prescribed in namoptions!) corresponds to criterion 0.95
!    real :: di=0.100       !< delta at inlet  (should be prescribed in namoptions!)
    real :: di_test          !< measured delta at the inlet
    real :: dti_test         !< measured deltat at the inlet
    real :: dr        !< delta at recycle station
    real :: dti       !< delta_t at inlet
    real :: dtr       !< delta_t at recycle
    real :: thetai    !< momentum thickness at inlet
    real :: thetar    !< momentum thickness at recycle
    real :: thetati   !< enthalpy thickness at inlet
    real :: thetatr   !< enthalpy thickness at recycle
    real :: utaui     !< u_tau at inlet
    real :: utaur     !< u_tau at recycle
    real :: ttaui     !< t_tau at inlet
    real :: ttaur     !< t_tau at recycle
    real :: lmoi      !< Obukhov length at inlet
    real :: lmor      !< Obukhov length at recycle
    real :: q0        !< wall heat flux at recycle
    real :: deltat=0.  !< full time step (set to zero at start of sim)
    real :: ubulk=0.   !< Bulk velocity (to be determined at first time step)   
    real :: totalu=0.  !< Bulk velocity inlet
    real :: totaluold=0.  !< old bulk velocity inlet
    real :: ddispdx=0.    !< spatial variation of displacement thickness (d/dx(delta*))
    real :: ddispdxold=0. !< old ddispdx
    real :: wtop=0.       !< mean vertical velocity at top read from zgrid.inf
    real :: xfm       !< mean (xf)
    real :: xf2m      !< mean (xf^2)
!    real :: dtin=0.0055  !< time step used in inletgenerator
    real :: dtin      !< time step used in inletgenerator
    real :: elapstep=0.  !< elapsed time in this time step. (used in time interpolation whe reading inletfiles)
    real :: totalreadu !< bulk velocity of inlet data (computed once at first time step)
    real :: iangle    !< inflow angle in radians (change with respect to inlet velocity that is read in)
    real :: iangledeg=0. !< inflow angle in degrees (change with respect to inlet velocity that is read in)

! Needed for interpolation in y-direction
    integer :: jgbin
    integer :: jgein
    integer :: jgtotinl      !< total number of cells in y-direction of inlet files (all procs together)   
    integer :: jbin
    integer :: jein
    integer :: jtotin
    integer :: jbdum
    integer :: jedum
    integer :: jtotdum
    integer :: filenumstart
    integer :: filenumend
    integer :: filestoread
    integer :: procinlo
    integer :: procinup
    integer :: jend
    integer :: jgend
    integer :: jbeg
    integer :: jgbeg

    real, allocatable :: yh(:)
    real, allocatable :: yf(:)
    real, allocatable :: yhin(:)
    real, allocatable :: yfin(:)
    real, allocatable :: yhdum(:)
    real, allocatable :: yfdum(:)
    integer, allocatable :: ylocupf(:)
    integer, allocatable :: yloclowf(:)
    integer, allocatable :: ylocuph(:)
    integer, allocatable :: yloclowh(:)
    
    real :: dyin


 
    integer :: irecy  !< ib + irecy is the i-index of recycle station
    
    integer :: nfile=0   !< file number to be read or written
    integer :: nstepread=1   !< time step number in file containing inlet plane
    integer :: rk3stepin=1   !< rk3step in inlet plane data
    integer :: kbin
    integer :: kein    
    integer :: nprocsinl     !< number of procs at used in inletdata files (only used for inletgen==2)    
    integer :: inlfactor     !< ratio of number of processors in this sim and in inlet data files (only used for inletgen==2)    
    logical :: lzinzsim = .true.      ! lzinzsim is .true. when inlet zgrid equals sim zgrid
    
! Inlet driver simulation variables - idriver - ae1212

    real, allocatable :: storeu0driver(:,:,:)
    real, allocatable :: storev0driver(:,:,:)
    real, allocatable :: storew0driver(:,:,:)
    real, allocatable :: storethl0driver(:,:,:)
    real, allocatable :: storee120driver(:,:,:)
    real, allocatable :: storeqt0driver(:,:,:)
    real, allocatable :: storesv0driver(:,:,:,:)
    real, allocatable :: storetdriver(:)
    real, allocatable :: u0driver(:,:)
    real, allocatable :: v0driver(:,:)
    real, allocatable :: w0driver(:,:)
    real, allocatable :: e120driver(:,:)
    real, allocatable :: tdriver(:)
    real, allocatable :: thl0driver(:,:)
    real, allocatable :: qt0driver(:,:)
    real, allocatable :: sv0driver(:,:,:)

    real, allocatable :: storeumdriver(:,:,:)
    real, allocatable :: umdriver(:,:)
    real, allocatable :: storevmdriver(:,:,:)
    real, allocatable :: vmdriver(:,:)
    real, allocatable :: storewmdriver(:,:,:)
    real, allocatable :: wmdriver(:,:)
    real, allocatable :: storee12mdriver(:,:,:)
    real, allocatable :: e12mdriver(:,:)
    real, allocatable :: storethlmdriver(:,:,:)
    real, allocatable :: thlmdriver(:,:)
    real, allocatable :: storeqtmdriver(:,:,:)
    real, allocatable :: qtmdriver(:,:)
    real, allocatable :: storesvmdriver(:,:,:,:)
    real, allocatable :: svmdriver(:,:,:)
    
    integer :: irecydriver
    integer :: nstepreaddriver=0
end module
