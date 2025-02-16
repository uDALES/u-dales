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

module advection_upw
   implicit none
   save

   contains

      !> Advection at cell center through upwind scheme
#if defined(_GPU)
      attributes(global) subroutine advecc_upw_cuda(hi, hj, hk, putin, putout)
         use modcuda, only: ib_d, ie_d, jb_d, je_d, kb_d, ke_d, dxfci_d, dyi_d, dzfci_d, &
                            u0_d, v0_d, w0_d, &
                            tidandstride
         implicit none
         integer, value, intent(in) :: hi, hj, hk
         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d - hk:ke_d + hk), intent(in)    :: putin
         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d     :ke_d + hk), intent(inout) :: putout
         integer :: i, j, k, tidx, tidy, tidz, stridex, stridey, stridez
         real    :: fluxr, fluxl, fluxb, fluxf, fluxu, fluxd
         
         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
         
         do k = tidz, ke_d, stridez
            do j = tidy, je_d, stridey
               do i = tidx, ie_d, stridex
                   if (u0_d(i+1, j, k) > 0) then   
                      fluxr = putin(i, j, k)
                   else
                      fluxr = putin(i + 1, j, k)
                   end if
 
                   if (u0_d(i, j, k) > 0) then
                      fluxl = putin(i - 1, j, k)
                   else
                      fluxl = putin(i, j, k)
                   end if
 
                   if (v0_d(i, j+1, k) > 0) then
                      fluxb = putin(i, j, k)
                   else
                      fluxb = putin(i, j + 1, k)
                   end if
 
                   if (v0_d(i, j, k) > 0) then
                      fluxf = putin(i, j - 1, k)
                   else
                      fluxf = putin(i, j, k)
                   end if
 
                   if (w0_d(i, j, k+1) > 0) then
                      fluxu = putin(i, j, k)
                   else
                      fluxu = putin(i, j, k + 1)
                   end if
 
                   if (w0_d(i, j, k) > 0) then
                      fluxd = putin(i, j, k - 1)
                   else
                      fluxd = putin(i, j, k)
                   end if
 
                   putout(i, j, k) = putout(i, j, k) &
                                     - (u0_d(i + 1, j, k)*fluxr - u0_d(i, j, k)*fluxl)*dxfci_d(i) & ! -d(uc)/dx (stretched grid)
                                     - (v0_d(i, j + 1, k)*fluxb - v0_d(i, j, k)*fluxf)*dyi_d &      ! -d(vc)/dy (no stretched grid)
                                     - (w0_d(i, j, k + 1)*fluxu - w0_d(i, j, k)*fluxd)*dzfci_d(k)   ! -d(wc)/dz (stretched grid)
               end do
            end do
         end do
      end subroutine advecc_upw_cuda
#else
      subroutine advecc_upw(hi, hj, hk, putin, putout)
         use modglobal, only: ib, ie, ih, jb, je, jh, kb, ke, dyi, dxfci, dzfci
         use modfields, only: u0, v0, w0
         implicit none

         integer, intent(in) :: hi !< size of halo in i
         integer, intent(in) :: hj !< size of halo in j
         integer, intent(in) :: hk !< size of halo in k
         real, dimension(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk), intent(in)    :: putin  !< Input: the cell centered field
         real, dimension(ib - hi:ie + hi, jb - hj:je + hj, kb     :ke + hk), intent(inout) :: putout !< Output: the tendency

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
                  end if
               end do
            end do
         end do

        do k = kb, ke
            do j = jb, je
               do i = ib, ie
                  putout(i, j, k) = putout(i, j, k) - &
                                    (u0(i + 1, j, k)*put(i + 1, j, k) - u0(i, j, k)*put(i, j, k))*dxfci(i)
               end do
            end do
        end do

        do k = kb, ke
           do j = jb, je + 1
              do i = ib, ie
                 if (v0(i, j, k) > 0) then
                    put(i, j, k) = putin(i, j - 1, k)
                 else
                    put(i, j, k) = putin(i, j, k)
                 end if
              end do
           end do
        end do
        do k = kb, ke
           do j = jb, je
              do i = ib, ie
                 putout(i, j, k) = putout(i, j, k) - &
                                   (v0(i, j + 1, k)*put(i, j + 1, k) - v0(i, j, k)*put(i, j, k))*dyi
              end do
           end do
        end do

        do k = kb, ke + 1
           do j = jb, je
              do i = ib, ie
                 if (w0(i, j, k) > 0) then
                    put(i, j, k) = putin(i, j, k - 1)
                 else
                    put(i, j, k) = putin(i, j, k)
                 end if
              end do
           end do
        end do
        do k = kb, ke
           do j = jb, je
              do i = ib, ie
                 putout(i, j, k) = putout(i, j, k) - &
                                  (w0(i, j, k + 1)*put(i, j, k + 1) - w0(i, j, k)*put(i, j, k))*dzfci(k)
              end do
           end do
        end do

        deallocate (put)

     end subroutine advecc_upw
#endif

end module advection_upw
