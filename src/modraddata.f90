!> Radiation switches, parameters and fields (subset of DALES dev 075ac4e modraddata).
module modraddata
  implicit none
  save
  integer, parameter :: irad_none = 0     !< no radiation
  integer, parameter :: irad_par  = 2     !< parametrised longwave (DALES radpar)
  integer :: iradiation = irad_none
  real    :: timerad    = 0.              !< update interval [s]; 0 = every step
  logical :: rad_longw  = .true.          !< longwave on
  logical :: rad_shortw = .false.         !< not supported: aborts
  logical :: rad_smoke  = .false.         !< absorber is scalar isvsmoke instead of ql0
  logical :: rad_ls     = .false.         !< DALES compatibility only (thlpcar is applied in forces)
  logical :: lrad_dycoms = .false.        !< DYCOMS free-troposphere term (Stevens et al. 2005)
  integer :: isvsmoke   = 1               !< scalar index of the smoke absorber
  real    :: rka        = 130.            !< extinction coefficient kappa [m2/kg]
  real    :: dlwtop     = 74.             !< F0 [W/m2]
  real    :: dlwbot     = 0.              !< F1 [W/m2]
  real    :: rad_dycoms_D       = 3.75e-6 !< large-scale divergence in the FT term [1/s]
  real    :: rad_dycoms_qtthres = 8.e-3   !< qt defining the local inversion height [kg/kg]
  real    :: tnextrad = 0.
  real, allocatable :: thlprad(:,:,:)     !< radiative thl tendency [K/s] (ib:ie,jb:je,kb:ke)
  real, allocatable :: lwd(:,:,:)         !< F0 term of net upward LW flux at faces [W/m2] (kb:ke+1)
  real, allocatable :: lwu(:,:,:)         !< F1 term of net upward LW flux at faces [W/m2] (kb:ke+1)
end module modraddata
