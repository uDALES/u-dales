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
  logical :: rad_ls     = .false.         !< DALES compatibility only (thlpcar is applied in forces)
  character(16) :: rad_absorber = 'ql'    !< absorber: 'ql' (liquid water) or 'sv' (scalar rad_absorber_sv)
  integer :: rad_absorber_sv = 0          !< scalar index of the absorber if rad_absorber = 'sv'; required
  real    :: rka        = 0.              !< extinction coefficient kappa [m2/kg]; required (> 0)
  real    :: dlwtop     = 0.              !< F0 [W/m2]; dlwtop or dlwbot required (> 0)
  real    :: dlwbot     = 0.              !< F1 [W/m2]
  real    :: tnextrad = 0.
  real, allocatable :: thlprad(:,:,:)     !< radiative thl tendency [K/s] (ib:ie,jb:je,kb:ke)
  real, allocatable :: lwd(:,:,:)         !< F0 term of net upward LW flux at faces [W/m2] (kb:ke+1)
  real, allocatable :: lwu(:,:,:)         !< F1 term of net upward LW flux at faces [W/m2] (kb:ke+1)
end module modraddata
