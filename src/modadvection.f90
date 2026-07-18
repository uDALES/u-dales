
!> \file advection.f90
!!  Advection management

!>
!!  Advection management
!! \par Revision list
!! variable x-grid now possible
!! Thijs Heus, Chiel van Heerwaarden, 15 June 2007
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

!> Advection redirection function
module modadvection

   implicit none

contains

subroutine advection

   use modglobal, only:lmoist, nsv, iadv_mom, iadv_tke, iadv_thl, iadv_qt, iadv_sv, &
      iadv_cd2, iadv_kappa, iadv_upw, &
      ltempeq, ih, jh, kh, ihc, jhc, khc, kb, ke, ib, ie, jb, je
   use modfields, only:u0, up, v0, vp, w0, wp, e120, e12p, thl0, thl0c, thlp, thlpc, qt0, qtp, sv0, svp
   use modsubgriddata, only:loneeqn
   implicit none
   integer :: n

   select case (iadv_mom)
   case (iadv_cd2)
     call advecu_2nd(u0,up)
     call advecv_2nd(v0,vp)
     call advecw_2nd(w0,wp)
   case default
      write(0, *) "ERROR: Unknown advection scheme"
      stop 1
   end select

   if (loneeqn) then
      select case (iadv_tke)
      case (iadv_cd2)
         call advecc_2nd(ih, jh, kh, e120, e12p)
      case default
         write(0, *) "ERROR: Unknown advection scheme"
         stop 1
      end select
   end if

   select case (iadv_thl)
   case (iadv_cd2)
      if (ltempeq) call advecc_2nd(ih, jh, kh, thl0, thlp)
   case (iadv_kappa)
      thlpc(ib:ie,jb:je,kb:ke) = thlp(ib:ie,jb:je,kb:ke)
      if (ltempeq) call advecc_kappa(ihc, jhc, khc, thl0c, thlpc)
      thlp(ib:ie,jb:je,kb:ke) = thlpc(ib:ie,jb:je,kb:ke)
   case default
      write(0, *) "ERROR: Unknown advection scheme"
      stop 1
   end select

   if (lmoist) then
      select case (iadv_qt)
      case (iadv_cd2)
         call advecc_2nd(ih, jh, kh, qt0, qtp)
      case default
         write(0, *) "ERROR: Unknown advection scheme"
         stop 1
      end select
   end if
   do n = 1, nsv
      select case (iadv_sv (n))
      case (iadv_cd2)
         call advecc_2nd(ihc, jhc, khc, sv0(:, :, :, n), svp(:, :, :, n))
      case (iadv_kappa)
         call advecc_kappa(ihc, jhc, khc, sv0(:, :, :, n), svp(:, :, :, n))
      case (iadv_upw)
         call advecc_upw(ihc, jhc, khc, sv0(:, :, :, n), svp(:, :, :, n))
      case default
         write(0, *) "ERROR: Unknown advection scheme"
         stop 1
      end select
   end do

end subroutine advection

subroutine advecc_2nd(hi, hj, hk, putin, putout)

   use modglobal, only:kb, ke, ib, ie, jb, je, dxi5, dyi5, dzf, dzhi, dzfi5
   use modfields, only:u0, v0, w0
   implicit none

   integer, intent(in) :: hi !< size of halo in i
   integer, intent(in) :: hj !< size of halo in j
   integer, intent(in) :: hk !< size of halo in k
   real, dimension(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk), intent(in)  :: putin !< Input: the cell centered field
   real, dimension(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk), intent(inout) :: putout !< Output: the tendency

   integer :: i, j, k, ip, im, jp, jm, kp, km

   ! Single fused pass, i innermost (column-major; see #330). The horizontal and
   ! vertical updates are kept as two statements so the rounding sequence is
   ! identical to the former two-pass version.
   do k = kb, ke
      km = k - 1
      kp = k + 1
      do j = jb, je
         jm = j - 1
         jp = j + 1
         do i = ib, ie
            im = i - 1
            ip = i + 1
            putout(i, j, k) = putout(i, j, k) - ( &
                              ( &
                              u0(ip, j, k)*(putin(ip, j, k) + putin(i, j, k)) &
                            - u0(i, j, k)*(putin(im, j, k) + putin(i, j, k)) & ! d(uc)/dx
                              )*dxi5 &
                            + ( & !
                              v0(i, jp, k)*(putin(i, jp, k) + putin(i, j, k)) &
                            - v0(i, j, k)*(putin(i, jm, k) + putin(i, j, k)) & ! d(vc)/dy
                              )*dyi5)
            putout(i, j, k) = putout(i, j, k) - ( &
                              w0(i, j, kp)*(putin(i, j, kp)*dzf(k) + putin(i, j, k)*dzf(kp))*dzhi(kp) &
                            - w0(i, j, k)*(putin(i, j, km)*dzf(k) + putin(i, j, k)*dzf(km))*dzhi(k) &
                              )*dzfi5(k)
         end do
      end do
   end do

end subroutine advecc_2nd

!> Advection at the u point.
subroutine advecu_2nd(putin, putout)

   use modglobal, only:ih, ib, ie, jb, je, jh, kb, ke, kh, dxi, dxiq, dyiq, dzf, dzfi5, dzhi
   use modfields, only:u0, v0, w0, pres0
   implicit none

   real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh), intent(in)  :: putin !< Input: the u-field
   real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh), intent(inout) :: putout !< Output: the tendency

   integer :: i, j, k, ip, im, jp, jm, kp, km

   ! Single fused pass, i innermost; two statements preserve the two-pass
   ! rounding sequence (see advecc_2nd / #330).
   do k = kb, ke
      km = k - 1
      kp = k + 1
      do j = jb, je
         jm = j - 1
         jp = j + 1
         do i = ib, ie
            im = i - 1
            ip = i + 1
            putout(i, j, k) = putout(i, j, k) - ( &
                              ( &
                              (putin(i, j, k) + putin(ip, j, k))*(u0(i, j, k) + u0(ip, j, k)) &
                            - (putin(i, j, k) + putin(im, j, k))*(u0(i, j, k) + u0(im, j, k)) & ! d(uu)/dx
                              )*dxiq &
                            + ( &
                              (putin(i, j, k) + putin(i, jp, k))*(v0(i, jp, k) + v0(im, jp, k)) &
                            - (putin(i, j, k) + putin(i, jm, k))*(v0(i, j, k) +  v0(im, j, k)) & ! d(vu)/dy
                              )*dyiq) &
                            - ((pres0(i, j, k) - pres0(i - 1, j, k))*dxi) ! - dp/dx
            putout(i, j, k) = putout(i, j, k) - ( &
                              (putin(i, j, kp)*dzf(k) + putin(i, j, k)*dzf(kp))*dzhi(kp) &
                            * (w0(i, j, kp) + w0(im, j, kp)) &
                            - (putin(i, j, k)*dzf(km) + putin(i, j, km)*dzf(k))*dzhi(k) &
                            * (w0(i, j, k) + w0(im, j, k)) &
                              )*0.5*dzfi5(k)
         end do
      end do
   end do

end subroutine advecu_2nd

!> Advection at the v point.
subroutine advecv_2nd(putin, putout)

   use modglobal, only:ih, ib, ie, jh, jb, je, kb, ke, kh, dxiq, dyiq, dzf, dzfi5, dzhi, dyi
   use modfields, only:u0, v0, w0, pres0
   implicit none

   real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh), intent(in)  :: putin !< Input: the v-field
   real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh), intent(inout) :: putout !< Output: the tendency

   integer :: i, j, k, ip, im, jp, jm, kp, km
   do k = kb, ke
      km = k - 1
      kp = k + 1
      do j = jb, je
         jm = j - 1
         jp = j + 1
         do i = ib, ie
            im = i - 1
            ip = i + 1

            putout(i, j, k) = putout(i, j, k) - ( &
                              ( &
                              (u0(ip, j, k) + u0(ip, jm, k))*(putin(i, j, k) + putin(ip, j, k)) &
                            - (u0(i, j, k)  + u0(i, jm, k)) *(putin(i, j, k) + putin(im, j, k)) & ! d(uv)/dx
                              )*dxiq &
                            + ( &
                              ( v0(i, jp, k) + v0(i, j, k))*(putin(i, j, k) + putin(i, jp, k)) &
                            - (v0(i, jm, k) + v0(i, j, k))*(putin(i, j, k) + putin(i, jm, k)) & ! d(vv)/dy
                              )*dyiq &
                              ) &
                            - ((pres0(i, j, k) - pres0(i, jm, k))*dyi) ! - dp/dy
            putout(i, j, k) = putout(i, j, k) - ( &
                              (w0(i, j, kp) + w0(i, jm, kp)) &
                            * (putin(i, j, kp)*dzf(k) + putin(i, j, k)*dzf(kp))*dzhi(kp) &
                            - (w0(i, j, k) + w0(i, jm, k)) &
                            * (putin(i, j, km)*dzf(k) + putin(i, j, k)*dzf(km))*dzhi(k) &
                              )*0.5*dzfi5(k)
         end do
      end do
   end do

end subroutine advecv_2nd

!> Advection at the w point.
subroutine advecw_2nd(putin, putout)

   use modglobal, only:ih, ib, ie, jh, jb, je, kb, ke, kh, dxiq, dyiq, dzf, dzhi, dzhiq
   use modfields, only:u0, v0, w0, pres0
   ! use modmpi, only : myid
   implicit none

   real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh), intent(in)  :: putin !< Input: the w-field
   real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh), intent(inout) :: putout !< Output: the tendency

   integer :: i, j, k, ip, im, jp, jm, kp, km

   do k = kb + 1, ke
      km = k - 1
      kp = k + 1
      do j = jb, je
         jm = j - 1
         jp = j + 1
         do i = ib, ie
            im = i - 1
            ip = i + 1

            putout(i, j, k) = putout(i, j, k) - ( &
                              ( &
                              (putin(ip, j, k) + putin(i, j, k))*(dzf(km)*u0(ip, j, k) + dzf(k)*u0(ip, j, km)) &
                            - (putin(i, j, k)  + putin(im, j, k))*(dzf(km)*u0(i, j, k) + dzf(k)*u0(i, j, km)) &
                              )*dxiq*dzhi(k) & ! d(uw)/dx
                            + ( &
                              (putin(i, jp, k) + putin(i, j, k))*(dzf(km)*v0(i, jp, k) + dzf(k)*v0(i, jp, km)) &
                            - (putin(i, j, k) + putin(i, jm, k))*(dzf(km)*v0(i, j, k) + dzf(k)*v0(i, j, km)) &
                              )*dyiq*dzhi(k) & ! d(vw)/dy
                            + ( &
                              (putin(i, j, k) + putin(i, j, kp))*(w0(i, j, k) + w0(i, j, kp)) &
                            - (putin(i, j, k) + putin(i, j, km))*(w0(i, j, k) + w0(i, j, km)) &
                              )*dzhiq(k) & ! d(ww)/dz
                              ) &
                            - ((pres0(i, j, k) - pres0(i, j, km))*dzhi(k)) ! - dp/dz

         end do
      end do
   end do
end subroutine advecw_2nd

  subroutine advecc_kappa(hi, hj, hk, var, varp)

!  use modglobal, only : i1,i2,ih,j1,j2,jh,k1,kmax,dxi,dyi,dzi
     use modglobal, only:ib, ie, jb, je, kb, ke, dxhci, dyi, dzhci, dxfc, dzfc, dxfci, dzfci
     use modfields, only:u0, v0, w0
     implicit none
     integer, intent(in) :: hi !< size of halo in i
     integer, intent(in) :: hj !< size of halo in j
     integer, intent(in) :: hk !< size of halo in k
     real, dimension(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk), intent(in)  :: var !< Input: the cell centered field
     real, dimension(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk), intent(inout) :: varp !< Output: the tendency
     real, dimension(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk)      ::  duml ! 3d dummy variable: lower cell side
     real, dimension(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk)      ::  dumu ! 3d dummy variable: upper cell side

     integer i, j, k
     real :: cf, d1, d2

     dumu(:, :, :) = 0.
     duml(:, :, :) = 0.
! -d(uc)/dx (stretched grid)
     do k = kb, ke
        do j = jb, je
           do i = ib, ie + 1
              if (u0(i, j, k) > 0) then
                 d1 = (var(i - 1, j, k) - var(i - 2, j, k))*dxhci(i - 1)
                 d2 = (var(i, j, k) - var(i - 1, j, k))*dxhci(i)
                 cf = var(i - 1, j, k)
              else
                 d1 = (var(i, j, k) - var(i + 1, j, k))*dxhci(i + 1)
                 d2 = (var(i - 1, j, k) - var(i, j, k))*dxhci(i)
                 cf = var(i, j, k)
              end if
              cf = cf + dxfc(i)*rlim(d1, d2)
              dumu(i - 1, j, k) = -cf*u0(i, j, k)*dxfci(i - 1) !swapped the -1s here !tg3315 !now also swapped the signs...
              duml(i, j, k) = cf*u0(i, j, k)*dxfci(i)
           end do
        end do
     end do

  varp(:,:,:) = varp(:,:,:) + dumu(:,:,:)+duml(:,:,:)

  dumu(:,:,:) = 0.
  duml(:,:,:) = 0.
! -d(vc)/dy (no stretched grid)
     do k = kb, ke
        do j = jb, je + 1
           do i = ib, ie
              if (v0(i, j, k) > 0) then
                 d1 = var(i, j - 1, k) - var(i, j - 2, k)
                 d2 = var(i, j, k) - var(i, j - 1, k)
                 cf = var(i, j - 1, k)
              else
                 d1 = var(i, j, k) - var(i, j + 1, k)
                 d2 = var(i, j - 1, k) - var(i, j, k)
                 cf = var(i, j, k)
              end if
              cf = cf + rlim(d1, d2)
              duml(i, j, k) = cf*v0(i, j, k)*dyi !tg3315
              dumu(i, j - 1, k) = -cf*v0(i, j, k)*dyi
           end do
        end do
     end do

  varp(:,:,:) = varp(:,:,:) + dumu(:,:,:)+duml(:,:,:)

  dumu(:,:,:) = 0.
  duml(:,:,:) = 0.
! -d(wc)/dz (stretched grid)
!  do k=kb,ke+1
     do k = kb + 1, ke + 1
        do j = jb, je
           do i = ib, ie
              if (w0(i, j, k) > 0) then
                 d1 = (var(i, j, k - 1) - var(i, j, k - 2))*dzhci(k - 1)
                 d2 = (var(i, j, k) - var(i, j, k - 1))*dzhci(k)
                 cf = var(i, j, k - 1)
              else
                 d1 = (var(i, j, k) - var(i, j, k + 1))*dzhci(k + 1)
                 d2 = (var(i, j, k - 1) - var(i, j, k))*dzhci(k)
                 cf = var(i, j, k)
              end if
              cf = cf + dzfc(k)*rlim(d1, d2)
              duml(i, j, k) = cf*w0(i, j, k)*dzfci(k) !tg3315 swapped
              dumu(i, j, k - 1) = -cf*w0(i, j, k)*dzfci(k - 1)
           end do
        end do
     end do

     varp(:,:,:) = varp(:,:,:) + dumu(:,:,:)+duml(:,:,:)

     return
  end subroutine advecc_kappa

!> Determination of the limiter function
  real function rlim(d1, d2)
     use modglobal, only:eps1
     implicit none
     real, intent(in) :: d1 !< Scalar flux at 1.5 cells upwind
     real, intent(in) :: d2 !< Scalar flux at 0.5 cells upwind

     real ri, phir

     ri = (d2 + eps1)/(d1 + eps1)
     phir = max(0., min(2.*ri, min(1./3.+2./3.*ri, 2.)))
     rlim = 0.5*phir*d1
  end function rlim

subroutine advecc_upw(hi, hj, hk, putin, putout)

   use modglobal, only:ib, ie, jb, je, kb, ke, dyi, dxfci, dzfci
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

end module modadvection
