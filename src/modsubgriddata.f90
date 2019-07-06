!!> \file modsubdata.f90
!!!  Provides variable definitions for Calculates and applies the Sub Filter Scale diffusion
!
!>
!!  Calculates and applies the Sub Filter Scale diffusion
!>
!!  \author Jasper Tomas, TU Delft
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  \todo Documentation

module modsubgriddata
implicit none
save
! private
  logical :: ldelta       = .false. !<  switch for subgrid length formulation (on/off)
  logical :: lmason       = .false. !<  switch for decreased length scale near the surface
  logical :: lsmagorinsky = .false. !<  switch for smagorinsky subgrid scheme
  logical :: lvreman      = .false. !<  switch for Vreman (2004) subgrid scheme
  logical :: lbuoycorr    = .false. !<  switch for buoyancy correction in Vreman (2004) subgrid scheme
  logical :: loneeqn      = .false. !<  switch for one-eqn subgrid scheme
  real :: cf      = 2.5  !< filter constant
  real :: Rigc    = 0.25 !< critical Richardson number
  real :: Prandtl = 0.333
!  real :: Prandtl = 3.0
!  real :: prandtli= 1./3.
  real :: prandtli
  real :: cm      = 0.12
  real :: cn      = 0.76
  real :: ch1     = 1.
  real :: ch2     = 2.
  real :: ce1     = 0.19
  real :: ce2     = 0.51
  real :: cs      = -1.
  real :: nmason  = 2.   !< exponent in Mason correction function
  real :: alpha_kolm  = 1.5     !< factor in Kolmogorov expression for spectral energy
  real :: beta_kolm   = 1.      !< factor in Kolmogorov relation for temperature spectrum
!  real :: damp   = 1.      !< used in van driest damping function
  real :: dampmin   = 1e-10     !< maximum damping used in van driest/Piomelli damping function
  real :: c_vreman  = 0.07      !< model constant for subgrid-scale model by Vreman (2004)        
!  real :: c_vreman  = 0.025      !< model constant for subgrid-scale model by Vreman (2004) corresponds with smag_const=0.1       

  real, allocatable :: ekm(:,:,:)   !< k-coefficient for momentum
  real, allocatable :: ekh(:,:,:)   !< k-coefficient for heat and q_tot

  real, allocatable :: sbdiss(:,:,:)!< dissiation
  real, allocatable :: sbshr(:,:,:) !< shear production
  real, allocatable :: sbbuo(:,:,:) !< buoyancy production / destruction
  real, allocatable :: zlt(:,:,:)   !< filter width
  
  real, allocatable :: csz(:,:)       !< Smagorinsky constant
  real, allocatable :: damp(:,:,:)      !< used in van Driest/Piomelli damping function
end module modsubgriddata

