!!
!!  \author Jasper Tomas,TU Delft, 31 March 2014
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
module modinletdata
  implicit none
  save

    real :: ubulk=0.   !< Bulk velocity (to be determined at first time step)
    real :: vbulk=0.   !< Bulk velocity (to be determined at first time step)
    real :: iangle    !< inflow angle in radians (change with respect to inlet velocity that is read in)
    real :: iangledeg=0. !< inflow angle in degrees (change with respect to inlet velocity that is read in)

    integer :: irecy  !< ib + irecy is the i-index of recycle station

! Inlet driver simulation variables - idriver - ae1212

    real, allocatable :: storeu0driver(:,:,:)
    real, allocatable :: storev0driver(:,:,:)
    real, allocatable :: storew0driver(:,:,:)
    real, allocatable :: storethl0driver(:,:,:)
    real, allocatable :: storee120driver(:,:,:)
    real, allocatable :: storeqt0driver(:,:,:)
    real, allocatable :: storesv0driver(:,:,:,:)
    real, allocatable :: storetdriver(:)
    real, allocatable :: u0driver(:,:)
    real, allocatable :: v0driver(:,:)
    real, allocatable :: u0driverrot(:,:)
    real, allocatable :: v0driverrot(:,:)
    real, allocatable :: w0driver(:,:)
    real, allocatable :: e120driver(:,:)
    real, allocatable :: tdriver(:)
    real, allocatable :: thl0driver(:,:)
    real, allocatable :: qt0driver(:,:)
    real, allocatable :: sv0driver(:,:,:)

    real, allocatable :: storeumdriver(:,:,:)
    real, allocatable :: umdriver(:,:)
    real, allocatable :: storevmdriver(:,:,:)
    real, allocatable :: vmdriver(:,:)
    real, allocatable :: storewmdriver(:,:,:)
    real, allocatable :: wmdriver(:,:)
    real, allocatable :: storee12mdriver(:,:,:)
    real, allocatable :: e12mdriver(:,:)
    real, allocatable :: storethlmdriver(:,:,:)
    real, allocatable :: thlmdriver(:,:)
    real, allocatable :: storeqtmdriver(:,:,:)
    real, allocatable :: qtmdriver(:,:)
    real, allocatable :: storesvmdriver(:,:,:,:)
    real, allocatable :: svmdriver(:,:,:)

    integer :: irecydriver
    integer :: nstepreaddriver=0

    integer :: chunkreadctr = 1             ! chunk reading counter
    integer :: chunkread_s = 0              ! chunk reading loop start
    integer :: chunkread_e = 0              ! chunk reading loop end

end module
