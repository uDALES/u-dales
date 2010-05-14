!!> \file modsubdata.f90
!!!  Provides variable definitions for Calculates and applies the Sub Filter Scale diffusion
!
!>
!!  Calculates and applies the Sub Filter Scale diffusion
!>
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  \todo Documentation
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
! public :: ldelta, lmason,lsmagorinsky,cf, Rigc,prandtl, cm, cn, ch1, ch2, ce1, ce2, ekm,ekh, sbdiss,sbshr,sbbuo

  logical :: ldelta   = .false. !<  switch for subgrid length formulation (on/off)
  logical :: lmason   = .false. !<  switch for decreased length scale near the surface
  logical :: lsmagorinsky= .false. !<  switch for smagorinsky subgrid scheme
  logical :: ldynsub     = .false. !<  switch for dynamic subgrid scheme
  real :: cf      = 2.5  !< filter constant
  real :: Rigc    = 0.25 !< critical Richardson number
  real :: Prandtl = 3
  real :: cm      = 0.12
  real :: cn      = 0.76
  real :: ch1     = 1.
  real :: ch2     = 2.
  real :: ce1     = 0.19
  real :: ce2     = 0.51
  real :: cs      = -1
  real :: nmason  = 2.   !< exponent in Mason correction function
  real :: alpha_kolm   = 1.5     !< factor in Kolmogorov expression for spectral energy
  real :: beta_kolm    = 1.      !< factor in Kolmogorov relation for temperature spectrum

  real, allocatable :: ekm(:,:,:)  !<   k-coefficient for momentum
  real, allocatable :: ekh(:,:,:)  !<   k-coefficient for heat and q_tot
  real, allocatable :: sbdiss(:,:,:)!< dissiation
  real, allocatable :: sbshr(:,:,:) !< shear production
  real, allocatable :: sbbuo(:,:,:) !< buoyancy production / destruction
  real, allocatable :: zlt(:,:,:)  !<   filter width

  !CvH Allocate
  !CvH Dynamic subgrid model variables
  real, allocatable :: u_bar(:,:), v_bar(:,:), w_bar(:,:)
  real, allocatable :: u_hat(:,:), v_hat(:,:), w_hat(:,:)
  real, allocatable :: S11(:,:), S12(:,:), S13(:,:), S22(:,:), S23(:,:), S33(:,:)
  real, allocatable :: S11_bar(:,:), S12_bar(:,:), S13_bar(:,:), S22_bar(:,:), S23_bar(:,:), S33_bar(:,:)
  real, allocatable :: S11_hat(:,:), S12_hat(:,:), S13_hat(:,:), S22_hat(:,:), S23_hat(:,:), S33_hat(:,:)
  real, allocatable :: S_S11_bar(:,:), S_S12_bar(:,:), S_S13_bar(:,:), S_S22_bar(:,:), S_S23_bar(:,:), S_S33_bar(:,:)
  real, allocatable :: S_S11_hat(:,:), S_S12_hat(:,:), S_S13_hat(:,:), S_S22_hat(:,:), S_S23_hat(:,:), S_S33_hat(:,:)
  real, allocatable :: S(:,:), S_bar(:,:), S_hat(:,:)
  real, allocatable :: L11(:,:), L12(:,:), L13(:,:), L22(:,:), L23(:,:), L33(:,:)
  real, allocatable :: Q11(:,:), Q12(:,:), Q13(:,:), Q22(:,:), Q23(:,:), Q33(:,:)
  real, allocatable :: M11(:,:), M12(:,:), M13(:,:), M22(:,:), M23(:,:), M33(:,:)
  real, allocatable :: N11(:,:), N12(:,:), N13(:,:), N22(:,:), N23(:,:), N33(:,:)
  real, allocatable :: LM(:,:), MM(:,:), QN(:,:), NN(:,:)
  real, allocatable :: csz(:)
  real, allocatable :: weighttf1(:,:)
  real, allocatable :: weighttf2(:,:)

  integer           :: tf1 = 2
  integer           :: tf2 = 4

  real              :: const, beta, cs2_tf1, cs2_tf2
  real              :: LMav, LMavl, MMav, MMavl, QNav, QNavl, NNav, NNavl

end module modsubgriddata

