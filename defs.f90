!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1999-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
!
module defs

  integer :: nv=1, nv1=2, mb=1
  real    :: totalpower
  real, parameter :: R      = 287.04
  real, parameter :: Rm     = 461.5
!   real, parameter :: ep     = R/Rm
  real, parameter :: ep2    = Rm/R - 1.
  real, parameter :: cp     = 1005.
!   real, parameter :: cv     = cp-R
  real, parameter :: rcp    = R/cp
  real, parameter :: cpr    = cp/R
  real, parameter :: g      = 9.8
  real, parameter :: p00    = 1.e+05
  real, parameter :: p00i   = 1./p00
!   real, parameter :: omega  = 7.292e-05
!   real, parameter :: alvl   = 2.5e+06
!   real, parameter :: alvi   = 2.834e+06
  real, parameter :: rowt   = 1.e+3
!   real, parameter :: vonk   = 0.40
!   real, parameter :: stefan = 5.6696e-8
  real, parameter :: SolarConstant  = 1.365d+3
  real, parameter :: mair   = 28.967
  real, parameter :: pi     = 3.14159265358979323846264338327

end module defs
