!> \file modsurfdata.f90
!! Variable definitions and auxilary routines for the surface model

!>
!! Variable definitions and auxilary routines for surface model
!>
!! This routine should have no dependency on any other routine, save perhaps modglobal or modfields.
!!  \author Thijs Heus, MPI-M
!!  \todo Documentation
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


!!whole module should be removed and variables moved


module modsurfdata

! implicit none

SAVE
  ! Surface properties
  real, allocatable :: tskin      (:,:) !<  Skin temperature [K]
  real, allocatable :: qskin      (:,:) !<  Skin specific humidity [kg/kg]
  real              :: ps       = -1    !<  Surface pressure [Pa]

  ! Surface energy balance
  logical           :: lmostlocal  = .false.  !<  Switch to apply MOST locally to get local Obukhov length
  real, allocatable :: obl   (:,:)      !<  Obukhov length [m]
  real              :: oblav =   0.001           !<  Spatially averaged obukhov length [m]
  real, allocatable :: Cm    (:,:)      !<  Drag coefficient for momentum [-]
  real, allocatable :: Cs    (:,:)      !<  Drag coefficient for scalars [-]
  real, allocatable :: ustar (:,:)      !<  Friction velocity [m/s]
  real, allocatable :: thlflux (:,:)    !<  Kinematic temperature flux [K m/s]
  real, allocatable :: qtflux  (:,:)    !<  Kinematic specific humidity flux [kg/kg m/s]
  real, allocatable :: svflux  (:,:,:)  !<  Kinematic scalar flux [- m/s]

  ! Surface gradients of prognostic variables
  real, allocatable :: dudz  (:,:)      !<  U-wind gradient in surface layer [1/s]
  real, allocatable :: dvdz  (:,:)      !<  V-wind gradient in surface layer [1/s]
  real, allocatable :: dqtdz (:,:)      !<  Specific humidity gradient in surface layer [kg/kg/m]
  real, allocatable :: dthldz(:,:)      !<  Liquid water potential temperature gradient in surface layer [K/m]

  ! Surface properties in case of prescribed conditions (previous isurf 2, 3 and 4)
  real              :: thls     = -1.    !<  Surface liquid water potential temperature [K]
  real              :: thl_top  = -1.    !<  Surface liquid water potential temperature [K] at top wall

  real              :: qts     = -1.         !<  Surface specific humidity [kg/kg]
  real              :: qt_top  = -1.     !<  Top value of specific humidity [kg/kg]

  real              :: thvs    = -1.         !<  Surface virtual temperature [K]

  real, allocatable :: svs   (:)        !<  Surface scalar concentration [-]
  real, allocatable :: sv_top (:)       ! top scalar concentration concentrations

  real              :: z0    = -1.       !<  Surface roughness length [m]
  real              :: z0h   = -1.      !<  Surface roughness for heat [m]
  ! prescribed surface fluxes
  real              :: Cmav             !<  Average drag coefficient for momentum [-]
  real              :: Csav             !<  Average drag coefficient for scalars [-]
  real              :: horvel           !<  Average horizontal velocity at first level


  real              :: wtsurf = -1.      !<  Prescribed kinematic temperature flux [K m/s]
  real              :: wttop  = 0.

  real              :: wqtop  = 0.
  real              :: wqsurf = -1.      !<  Prescribed kinematic moisture flux [kg/kg m/s]

  real, allocatable :: wsvsurf(:) !<  Prescribed surface scalar(n) flux [- m/s]
  real, allocatable :: wsvtop(:)
  real :: wsvsurfdum(1:99) = 0. !<  Dummy variables as nsv allocated variable
  real :: wsvtopdum(1:99)  = 0.


end module modsurfdata
