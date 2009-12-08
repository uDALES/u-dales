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



module modsurfdata

! implicit none

SAVE
  integer :: isurf        = -1            !<   Flag for surface parametrization

  ! Soil properties

  ! Domain-uniform properties
  integer, parameter  :: ksoilmax = 4       !<  Number of soil layers [-]

  real              :: lambdasat          !<  heat conductivity saturated soil [W/m/K]
  real              :: Ke                 !<  Kersten number [-]

  real, allocatable :: zsoil  (:)         !<  Height of bottom soil layer from surface [m]
  real, allocatable :: dzsoil (:)         !<  Depth of soil layer [m]
  real, allocatable :: dzsoilh(:)         !<  Depth of soil layer between center of layers [m]

  ! Spatially varying properties
  real, allocatable :: lambda  (:,:,:)    !<  Heat conductivity soil layer [W/m/K]
  real, allocatable :: lambdah (:,:,:)    !<  Heat conductivity soil layer half levels [W/m/K]
  real, allocatable :: Dh      (:,:,:)    !<  Heat diffusivity
  real, allocatable :: phiw    (:,:,:)    !<  Water content soil matrix [-]
  real              :: phiwav  (ksoilmax)
  real, allocatable :: phitot  (:,:)      !<  Total soil water content [-]
  real, allocatable :: pCs     (:,:,:)    !<  Volumetric heat capacity [J/m3/K]
  real, allocatable :: rootf   (:,:,:)    !<  Root fraction per soil layer [-]
  real              :: rootfav (ksoilmax)
  real, allocatable :: tsoil   (:,:,:)    !<  Soil temperature [K]
  real              :: tsoilav (ksoilmax)
  real, allocatable :: tsoildeep (:,:)    !<  Soil temperature [K]
  real              :: tsoildeepav

  real, allocatable :: swdavn  (:,:,:)
  real, allocatable :: swuavn  (:,:,:)
  real, allocatable :: lwdavn  (:,:,:)
  real, allocatable :: lwuavn  (:,:,:)

  integer           :: nradtime  = 60

  ! Soil related constants [adapted from ECMWF]
  real, parameter   :: phi       = 0.472  !<  volumetric soil porosity [-]
  real, parameter   :: phifc     = 0.323  !<  volumetric moisture at field capacity [-]
  real, parameter   :: phiwp     = 0.171  !<  volumetric moisture at wilting point [-]

  real, parameter   :: pCm       = 2.19e6 !<  Volumetric soil heat capacity [J/m3/K]
  real, parameter   :: pCw       = 4.2e6  !<  Volumetric water heat capacity [J/m3/K]

  real, parameter   :: lambdadry = 0.190  !<  Heat conductivity dry soil [W/m/K]
  real, parameter   :: lambdasm  = 3.11   !<  Heat conductivity soil matrix [W/m/K]
  real, parameter   :: lambdaw   = 0.57   !<  Heat conductivity water [W/m/K]


  ! Land surface properties

  ! Surface properties
  real, allocatable :: z0m        (:,:) !<  Roughness length for momentum [m]
  real              :: z0mav = -1
  real, allocatable :: z0h        (:,:) !<  Roughness length for heat [m]
  real              :: z0hav = -1
  real, allocatable :: tskin      (:,:) !<  Skin temperature [K]
  real, allocatable :: tskinm     (:,:) !<  Skin temperature previous timestep [K]
  real, allocatable :: qskin      (:,:) !<  Skin specific humidity [kg/kg]
  real, allocatable :: albedo     (:,:) !<  Surface albedo [-]
  real              :: albedoav = -1
  real, allocatable :: LAI        (:,:) !<  Leaf area index vegetation [-]
  real              :: LAIav    = -1
  real, allocatable :: Cskin      (:,:) !<  Heat capacity skin layer [J]
  real              :: Cskinav  = -1
  real, allocatable :: lambdaskin (:,:) !<  Heat conductivity skin layer [W/m/K]
  real              :: lambdaskinav
  real              :: ps       = -1    !<  Surface pressure [Pa]

  ! Surface energy balance
  real, allocatable :: Qnet     (:,:)   !<  Net radiation [W/m2]
  real              :: Qnetav   = -1
  real, allocatable :: LE       (:,:)   !<  Latent heat flux [W/m2]
  real, allocatable :: H        (:,:)   !<  Sensible heat flux [W/m2]
  real, allocatable :: G0       (:,:)   !<  Ground heat flux [W/m2]
  real, allocatable :: ra       (:,:)   !<  Aerodynamic resistance [s/m]
  real, allocatable :: rs       (:,:)   !<  Vegetation resistance [s/m]
  real              :: rsisurf2 = 0.    !<  Vegetation resistance [s/m] if isurf2 is used
  real, allocatable :: rsmin    (:,:)   !<  Minimum vegetation resistance [s/m]
  real              :: rsminav = -1
  real, allocatable :: tendskin (:,:)   !<  Tendency of skin [W/m2]
  real, allocatable :: gD       (:,:)   !<  Response factor vegetation to vapor pressure deficit [-]
  real              :: gDav

  ! Turbulent exchange variables
  logical           :: lmostlocal  = .false.  !<  Switch to apply MOST locally to get local Obukhov length
  logical           :: lsmoothflux = .false.  !<  Create uniform sensible and latent heat flux over domain
  real, allocatable :: obl   (:,:)      !<  Obukhov length [m]
  real              :: oblav            !<  Spatially averaged obukhov length [m]
  real, allocatable :: Cm    (:,:)      !<  Drag coefficient for momentum [-]
  real, allocatable :: Cs    (:,:)      !<  Drag coefficient for scalars [-]
  real, allocatable :: ustar (:,:)      !<  Friction velocity [m/s]
  real, allocatable :: thlflux (:,:)    !<  Kinematic temperature flux [K m/s]
  real, allocatable :: qtflux  (:,:)    !<  Kinematic specific humidity flux [kg/kg m/s]
  real, allocatable :: svstar  (:,:,:)  !<  Turbulent scalar scale [- m/s]

  ! Surface gradients of prognostic variables
  real, allocatable :: dudz  (:,:)      !<  U-wind gradient in surface layer [1/s]
  real, allocatable :: dvdz  (:,:)      !<  V-wind gradient in surface layer [1/s]
  real, allocatable :: dqtdz (:,:)      !<  Specific humidity gradient in surface layer [kg/kg/m]
  real, allocatable :: dthldz(:,:)      !<  Liquid water potential temperature gradient in surface layer [K/m]

  ! Surface properties in case of prescribed conditions (previous isurf 2, 3 and 4)
  real              :: thls  = -1       !<  Surface liquid water potential temperature [K]
  real              :: qts              !<  Surface specific humidity [kg/kg]
  real              :: thvs             !<  Surface virtual temperature [K]
  real, allocatable :: svs   (:)        !<  Surface scalar concentration [-]
  real              :: z0    = -1       !<  Surface roughness length [m]

  ! prescribed surface fluxes
  real              :: ustin  = -1      !<  Prescribed friction velocity [m/s]
  real              :: wtsurf = -1      !<  Prescribed kinematic temperature flux [K m/s]
  real              :: wqsurf = -1      !<  Prescribed kinematic moisture flux [kg/kg m/s]
  real              :: wsvsurf(100)     !<  Prescribed surface scalar(n) flux [- m/s]

end module modsurfdata
