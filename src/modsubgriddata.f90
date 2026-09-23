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
!
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

module modsubgriddata
implicit none
save
! private
  logical :: ldelta       = .false. !<  switch for subgrid length formulation (on/off)
  logical :: lmason       = .false. !<  switch for decreased length scale near the surface
  logical :: lsmagorinsky = .false. !<  switch for smagorinsky subgrid scheme
  logical :: lvreman      = .true. !<  switch for Vreman (2004) subgrid scheme
  logical :: lbuoycorr    = .false. !<  switch for buoyancy correction in Vreman (2004) subgrid scheme
  logical :: loneeqn      = .false. !<  switch for one-eqn subgrid scheme
  logical :: lnarrowdamp  = .false. !<  switch for local Smagorinsky-type eddy-viscosity floor in geometrically
                                    !<  under-resolved (narrow) fluid cells listed in narrow_c.txt; the base SGS
                                    !<  model (e.g. Vreman) is unchanged elsewhere
  real    :: cs_narrow    = 0.2     !<  Smagorinsky constant used for the narrow-cell viscosity floor
  real    :: cs_floor     = 0.      !<  if > 0: global weak Smagorinsky-type lower bound on the Vreman eddy
                                    !<  viscosity, ekm >= (cs_floor*delta)**2 * |S|. In resolved turbulence the
                                    !<  Vreman viscosity exceeds this bound; it only engages where Vreman is
                                    !<  anomalously small relative to strain (grid-scale unidirectional jets).
  integer :: nnarrowpts   = 0       !<  total number of narrow cells (data lines in narrow_c.txt)
  integer :: nnarrowptsrank = 0     !<  number of narrow cells on this rank
  integer, allocatable :: narrowpts_loc(:,:) !< local (i,j,k) of narrow cells on this rank
  logical :: lnarrowcolfull = .false. !< extend the narrow-cell floor over the FULL column height (kb..ke) of every
                                      !< column containing narrow cells, so jets cannot escape the damped region
                                      !< by climbing above the masked levels
  integer :: nnarrowcolsrank = 0      !< number of such columns on this rank
  integer, allocatable :: narrowcols_loc(:,:) !< local (i,j) of those columns
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

