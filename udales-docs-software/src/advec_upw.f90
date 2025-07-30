!> \file advec_upw.f90
!!  Does advection with a 1st order upwind scheme.
!! \par Revision list
!! \par Authors
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

!> Advection at cell center
subroutine advecc_upw(hi, hj, hk, putin, putout)

   use modglobal, only:ib, ie, ih, jb, je, jh, kb, ke, dyi, dxfci, dzfci
   use modfields, only:u0, v0, w0
   implicit none

   integer, intent(in) :: hi !< size of halo in i
   integer, intent(in) :: hj !< size of halo in j
   integer, intent(in) :: hk !< size of halo in k
   real, dimension(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk), intent(in)  :: putin !< Input: the cell centered field
   real, dimension(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk), intent(inout) :: putout !< Output: the tendency

   real, allocatable, dimension(:, :, :) :: put
   integer :: i, j, k

   allocate (put(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk))

   do k = kb, ke
      do j = jb, je
         do i = ib, ie + 1
            if (u0(i, j, k) > 0) then
               put(i, j, k) = putin(i - 1, j, k)
            else
               put(i, j, k) = putin(i, j, k)
            endif
         enddo
      enddo
   enddo

   do k = kb, ke
      do j = jb, je
         do i = ib, ie
            putout(i, j, k) = putout(i, j, k) - &
                              (u0(i + 1, j, k)*put(i + 1, j, k) - u0(i, j, k)*put(i, j, k))*dxfci(i)
         enddo
      enddo
   enddo

   do k = kb, ke
      do j = jb, je + 1
         do i = ib, ie
            if (v0(i, j, k) > 0) then
               put(i, j, k) = putin(i, j - 1, k)
            else
               put(i, j, k) = putin(i, j, k)
            endif
         enddo
      enddo
   enddo
   do k = kb, ke
      do j = jb, je
         do i = ib, ie
            putout(i, j, k) = putout(i, j, k) - &
                              (v0(i, j + 1, k)*put(i, j + 1, k) - v0(i, j, k)*put(i, j, k))*dyi
         enddo
      enddo
   enddo

   do k = kb, ke + 1
      do j = jb, je
         do i = ib, ie
            if (w0(i, j, k) > 0) then
               put(i, j, k) = putin(i, j, k - 1)
            else
               put(i, j, k) = putin(i, j, k)
            endif
         enddo
      enddo
   enddo
   do k = kb, ke
      do j = jb, je
         do i = ib, ie
            putout(i, j, k) = putout(i, j, k) - &
                              (w0(i, j, k + 1)*put(i, j, k + 1) - w0(i, j, k)*put(i, j, k))*dzfci(k)
         enddo
      enddo
   enddo

   deallocate (put)

end subroutine advecc_upw

