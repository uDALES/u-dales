!> \file moduser.f90
!! A dummy file for cases where one wants additional forcings
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
module moduser

contains
subroutine initsurf_user
  implicit none
end subroutine initsurf_user

subroutine force_user
  implicit none
end subroutine force_user

subroutine rad_user
! DYCOMS version of rad_user adds a radiative forcing above the inversion layer
use modfields, only : qt0
use modglobal, only : i1,j1,kmax,dzh,zf
use modraddata, only : thlprad
implicit none
integer :: i,j,k
real    :: thres = 8e-3,a=1,D=3.75e-6
real    :: zi=0.0

  do i=2,j1
  do j=2,j1
  !Determine local BL-height following the specifications.
    do k=kmax,1,-1
      if (qt0(i,j,k) > thres) then
        zi =zf(k) + (thres-qt0(i,j,k)) &
                    *dzh(k+1)/(qt0(i,j,k+1)-qt0(i,j,k))
        exit
      endif
    enddo
    ! Apply the radiative flux divergence
    do k=kmax,1,-1
      if (zf(k)< zi) exit
      thlprad(i,j,k) = thlprad(i,j,k) + a*D/3.*((zf(k)-zi)**(1./3.)+zi*(zf(k)-zi)**(-2./3.))
    end do
  enddo
  enddo
end subroutine rad_user

subroutine micro_user
  implicit none
end subroutine micro_user

subroutine surf_user
  implicit none
end subroutine surf_user

end module moduser
