!----------------------------------------------------------------------------
! This file is part of DALES.
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
! Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!----------------------------------------------------------------------------
!
!
module modsurfdata
save
  integer :: isurf        = -1       !  flag for surface parametrization
  logical :: lneutraldrag = .false.  !  do not apply stability correction

  real, allocatable :: ustar (:,:)    !       friction velocity
  real, allocatable :: dudz  (:,:)    !       veloc. grad. in x-direction in surf. layer
  real, allocatable :: dvdz  (:,:)    !       veloc. grad. in y-direction in surf. layer
  real, allocatable :: tstar (:,:)    !       turb. temperature scale
  real, allocatable :: qstar (:,:)    !       turb. specific humidity scale
  real, allocatable :: dqtdz (:,:)    !       qt-gradient in surface layer
  real, allocatable :: dthldz(:,:)    !       thl-gradient in surface layer
  real :: thls   =-1          !       surface liq. water pot. temperature
  real :: qts              !       surface total water specific humidity
  real :: ps     =-1         !       surface pressure
  real :: thvs                       !   *theta_v constant used in buoyancy equation.
  
  real :: obl                         !       Obukhov length
  real :: z0     = -1   !       surface roughness
  real :: ustin  = -1     !       prescribed friction velocity
  real :: wtsurf = -1      !       prescribed surface thl-flux
  real :: wqsurf = -1      !       prescribed surface qt-flux

  real, allocatable :: svstar(:,:,:)!      concentration scale
  real, allocatable :: svs(:)         !      surface scalar variable concentration
  real              :: wsvsurf(100)     !      prescribed surface sv(n)-flux
end module modsurfdata