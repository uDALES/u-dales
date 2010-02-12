!> \file modraddata.f90
!! Variable definitions and auxilary routines for radiation

!>
!! Variable definitions and auxilary routines for radiation.
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



module modraddata

! implicit none

SAVE

  integer, parameter :: irad_none  = 0   !< 0=no radiation
  integer, parameter :: irad_full  = 1   !< 1=full radiation
  integer, parameter :: irad_par   = 2   !< 2=parameterized radiation
  integer, parameter :: irad_lsm   = 3   !< 3=simple surface radiation for land surface model
  integer, parameter :: irad_user  = 10  !< 10=user specified radiation

  logical :: rad_ls      = .true.   !< prescribed radiative forcing
  logical :: rad_longw   = .true.   !< parameterized longwave radiative forcing
  logical :: rad_shortw  = .true.   !< parameterized shortwave radiative forcing
  logical :: rad_smoke   = .false.  !< longwave divergence for smoke cloud
  logical :: useMcICA    = .true.   !< Use the Monte Carlo Independent Column Approach

  real              :: timerad = 0 !<  timescale of the radiation scheme
  real              :: tnext   = 0 !<  time of the first upcoming call of the radiation scheme
  real :: rho_air_mn = 1.1436 !< mean air density used in radiation computation
  real :: rka        = 130.   !< extinction coefficient in radpar scheme
  real :: dlwtop     = 74.    !< longwave radiative flux divergence at top of domain
  real :: dlwbot     = 0.     !< longwave radiative flux divergence near the surface
  real :: sw0        = 1100.0 !< direct component at top of the cloud (W/m^2), diffuse not possible
  real :: gc         = 0.85   !< asymmetry factor of droplet scattering angle distribution
  !CvH real :: sfc_albedo = 0.05   !< ground surface albedo
  real :: reff       = 1.e-5  !< cloud droplet effective radius (m)
  integer :: isvsmoke = 1     !< number of passive scalar to be used for optical depth calculation
  integer :: iradiation = irad_none !< Selection parameter for type of radiation scheme
  integer :: irad    = -1  !< Deprecated selection parameter for the type of radiation scheme


  real mu                    !< cosine of the solar zenith angle

  real, allocatable :: thlprad(:,:,:)!<   the radiative tendencies
  real, allocatable :: swd(:,:,:)    !<   shortwave downward radiative flux
  real, allocatable :: swu(:,:,:)    !<   shortwave upward radiative flux
  real, allocatable :: lwd(:,:,:)    !<   longwave downward radiative flux
  real, allocatable :: lwu(:,:,:)    !<   longwave upward radiative flux
  !CvH albedo real, allocatable :: albedo(:,:)   !<   Albedo


contains
!< Calculation of the cosine of the zenith angle
!< \param time UTC Time of the simulation
!< \param xday Day at the start of the simulation
!< \param xlat Latitude of the domain
!< \param xlon Longitude of the domain
  real function zenith(time, xday, xlat,xlon)
    use modglobal, only : pi
!     implicit none
    real, intent(in) :: time, xday, xlat, xlon
    real :: phi,el,obliq,xlam,declin,hora
    real :: day,daytime
    day    = xday + floor(time/86400.)
    daytime= mod(time,86400.)
    phi    = xlat * pi/180.
    el     = xlon * pi/180.
    obliq  = 23.45 * pi/180.
    xlam   = 4.88 + 0.0172 * day
    declin = asin(sin(obliq)*sin(xlam))
    hora   = el-pi + 2.*pi*(daytime/24.)
    zenith = max(0.,sin(declin)*sin(phi)+cos(declin)*cos(phi)* &
                                                         cos(hora))
  end function zenith

end module modraddata
