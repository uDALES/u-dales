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
  use modglobal, only : kmax,dzh,zf,zh, timee
  use modfields, only : qt0av,thlp,qtp,whls
  implicit none
  integer :: k,kinv,ktrot
  real    :: ztrot, zinv
  if (timee<5400) return
  do k=1,kmax
    if (qt0av(k)<6.5e-3) exit
  end do
  kinv = k-1
  zinv = zf(kinv) + (qt0av(kinv)-6.5e-3)/(qt0av(kinv)-qt0av(kinv+1))*dzh(k)
  do k=kinv,kmax
    if (zf(k) > zinv+300) exit
  end do
  ktrot = k-1
  ztrot = zf(ktrot) + (zf(ktrot)-(zinv+300))/(zf(ktrot)-zf(ktrot+1))*dzh(k)

  whls = 0.
  do k=1,kinv
    whls(k)     = -zh(k)/zinv*6.5e-3
    thlp(:,:,k) = thlp(:,:,k) -1.1575e-5*(3-zf(k)/zinv)
    qtp(:,:,k)  = qtp(:,:,k)  -1.58e-8*(3-zf(k)/zinv)
  end do
  do k=kinv+1,ktrot
    whls(k)     = 6.5e-3*(zh(k)-ztrot)/300
    thlp(:,:,k) = thlp(:,:,k) -1.1575e-5*(2.)
    qtp(:,:,k)  = qtp(:,:,k)  -1.58e-8*(2.)
  end do
end subroutine force_user

subroutine rad_user
  implicit none
end subroutine rad_user

subroutine micro_user
  implicit none
end subroutine micro_user

subroutine surf_user
  implicit none
end subroutine surf_user

end module moduser
