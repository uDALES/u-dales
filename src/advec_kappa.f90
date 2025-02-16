!> \file advec_kappa.f90
!!  Does advection with a kappa limiter scheme.
!! \par Revision list
!! \par Authors
!! \see Hundsdorfer et al 1995
!!
!! For advection of scalars that need to be strictly monotone (for example chemically reacting species)
!! the kappa scheme has been implemented:
!! \latexonly
!! \begin{eqnarray}
!!  F_{i-\frac{1}{2}}^{\kappa} &=& \fav{u}_{i-\frac{1}{2}}
!!  \left[\phi_{i-1}+\frac{1}{2}\kappa_{i-\frac{1}{2}}\left(\phi_{i-1}-\phi_{i-2}\right)\right],
!! \end{eqnarray}
!! in case $\fav{u}>0$. $\kappa_{i-\smfrac{1}{2}}$ serves as a switch between higher order advection and
!! first order upwind in case of strong upwind gradients of $\phi$.
!! \endlatexonly
!! This makes the scheme monotone, but also rather dissipative.
!!
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

module advection_kappa
   implicit none
   save

   contains
   
      !> Advection at cell center through kappa scheme
#if defined(_GPU)
      attributes(global) subroutine advecc_kappa_reset_cuda(hi, hj, hk)
         use modcuda, only: ie_d, je_d, ke_d, dumu_d, duml_d, tidandstride
         implicit none
         integer, value, intent(in) :: hi, hj, hk
         integer :: i, j, k, tidx, tidy, tidz, stridex, stridey, stridez
         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
         do k = tidz, ke_d + hk, stridez
            do j = tidy - hj, je_d + hj, stridey
               do i = tidx - hi, ie_d + hi, stridex
                  dumu_d(i,j,k) = 0.
                  duml_d(i,j,k) = 0.
               end do
            end do
         end do
      end subroutine advecc_kappa_reset_cuda

      ! -d(uc)/dx (stretched grid)
      attributes(global) subroutine advecc_kappa_ducdx_cuda(hi, hj, hk, var)
         use modcuda, only: ib_d, ie_d, jb_d, je_d, kb_d, ke_d, dxfc_d, dxfci_d, dxhci_d, &
                            u0_d, dumu_d, duml_d, &
                            tidandstride
         implicit none
         integer, value, intent(in) :: hi, hj, hk
         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d - hk:ke_d + hk), intent(in) :: var
         real    :: cf, d1, d2
         integer :: i, j, k, tidx, tidy, tidz, stridex, stridey, stridez

         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

         do k = tidz, ke_d, stridez
            do j = tidy, je_d, stridey
               do i = tidx, ie_d + 1, stridex
                  if (u0_d(i, j, k) > 0) then
                     d1 = (var(i - 1, j, k) - var(i - 2, j, k))*dxhci_d(i - 1)
                     d2 = (var(i, j, k) - var(i - 1, j, k))*dxhci_d(i)
                     cf = var(i - 1, j, k)
                  else
                     d1 = (var(i, j, k) - var(i + 1, j, k))*dxhci_d(i + 1)
                     d2 = (var(i - 1, j, k) - var(i, j, k))*dxhci_d(i)
                     cf = var(i, j, k)
                  end if
                  cf = cf + dxfc_d(i)*rlim_cuda(d1, d2)
                  dumu_d(i - 1, j, k) = -cf*u0_d(i, j, k)*dxfci_d(i - 1)
                  duml_d(i, j, k) = cf*u0_d(i, j, k)*dxfci_d(i)
               end do
            end do
         end do
      end subroutine advecc_kappa_ducdx_cuda

      ! -d(vc)/dy (no stretched grid)
      attributes(global) subroutine advecc_kappa_dvcdy_cuda(hi, hj, hk, var)
         use modcuda, only: ib_d, ie_d, jb_d, je_d, kb_d, ke_d, dyi_d, &
                            v0_d, dumu_d, duml_d, &
                            tidandstride
         implicit none
         integer, value, intent(in) :: hi, hj, hk
         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d - hk:ke_d + hk), intent(in) :: var
         real    :: cf, d1, d2
         integer :: i, j, k, tidx, tidy, tidz, stridex, stridey, stridez

         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

         do k = tidz, ke_d, stridez
            do j = tidy, je_d + 1, stridey
               do i = tidx, ie_d, stridex
                  if (v0_d(i, j, k) > 0) then
                     d1 = var(i, j - 1, k) - var(i, j - 2, k)
                     d2 = var(i, j, k) - var(i, j - 1, k)
                     cf = var(i, j - 1, k)
                  else
                     d1 = var(i, j, k) - var(i, j + 1, k)
                     d2 = var(i, j - 1, k) - var(i, j, k)
                     cf = var(i, j, k)
                  end if
                  cf = cf + rlim_cuda(d1, d2)
                  duml_d(i, j, k) = cf*v0_d(i, j, k)*dyi_d
                  dumu_d(i, j - 1, k) = -cf*v0_d(i, j, k)*dyi_d
               end do
            end do
         end do
      end subroutine advecc_kappa_dvcdy_cuda

      ! -d(wc)/dz (stretched grid)
      attributes(global) subroutine advecc_kappa_dwcdz_cuda(hi, hj, hk, var)
         use modcuda, only: ib_d, ie_d, jb_d, je_d, kb_d, ke_d, dzfc_d, dzfci_d, dzhci_d, &
                            w0_d, dumu_d, duml_d, &
                            tidandstride
         implicit none
         integer, value, intent(in) :: hi, hj, hk
         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d - hk:ke_d + hk), intent(in) :: var
         real    :: cf, d1, d2
         integer :: i, j, k, tidx, tidy, tidz, stridex, stridey, stridez

         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

         do k = tidz + 1, ke_d + 1, stridez
            do j = tidy, je_d, stridey
               do i = tidx, ie_d, stridex
                  if (w0_d(i, j, k) > 0) then
                     d1 = (var(i, j, k - 1) - var(i, j, k - 2))*dzhci_d(k - 1)
                     d2 = (var(i, j, k) - var(i, j, k - 1))*dzhci_d(k)
                     cf = var(i, j, k - 1)
                  else
                     d1 = (var(i, j, k) - var(i, j, k + 1))*dzhci_d(k + 1)
                     d2 = (var(i, j, k - 1) - var(i, j, k))*dzhci_d(k)
                     cf = var(i, j, k)
                  end if
                  cf = cf + dzfc_d(k)*rlim_cuda(d1, d2)
                  duml_d(i, j, k) = cf*w0_d(i, j, k)*dzfci_d(k)
                  dumu_d(i, j, k - 1) = -cf*w0_d(i, j, k)*dzfci_d(k - 1)
               end do
            end do
         end do
      end subroutine advecc_kappa_dwcdz_cuda

      attributes(global) subroutine advecc_kappa_add_cuda(hi, hj, hk, varp)
         use modcuda, only: ib_d, ie_d, jb_d, je_d, kb_d, ke_d, dumu_d, duml_d, tidandstride
         implicit none
         integer, value, intent(in) :: hi, hj, hk
         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d:ke_d + hk), intent(inout) :: varp
         integer :: i, j, k, tidx, tidy, tidz, stridex, stridey, stridez
         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
         do k = tidz, ke_d + hk, stridez
            do j = tidy - hj, je_d + hj, stridey
               do i = tidx - hi, ie_d + hi, stridex
                  varp(i,j,k) = varp(i,j,k) + dumu_d(i,j,k) + duml_d(i,j,k)
               end do
            end do
         end do
      end subroutine advecc_kappa_add_cuda

      !> Determination of the limiter function
      attributes(device) real function rlim_cuda(d1, d2)
         use modcuda, only: eps1_d
         implicit none
         real, intent(in) :: d1 !< Scalar flux at 1.5 cells upwind
         real, intent(in) :: d2 !< Scalar flux at 0.5 cells upwind
         real :: ri, phir
         ri = (d2 + eps1_d)/(d1 + eps1_d)
         phir = max(0., min(2.*ri, min(1./3.+2./3.*ri, 2.)))
         rlim_cuda = 0.5*phir*d1
      end function rlim_cuda

#else
      subroutine advecc_kappa(hi, hj, hk, var, varp)
         use modglobal, only: ib, ie, ihc, jb, je, jhc, kb, ke, khc, dxhci, dyi, dzhci, dxfc, dzfc, dxfci, dzfci
         use modfields, only: u0, v0, w0
         implicit none

         integer, intent(in) :: hi !< size of halo in i
         integer, intent(in) :: hj !< size of halo in j
         integer, intent(in) :: hk !< size of halo in k
         real, dimension(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk), intent(in)    :: var  !< Input: the cell centered field
         real, dimension(ib - hi:ie + hi, jb - hj:je + hj, kb     :ke + hk), intent(inout) :: varp !< Output: the tendency
         
         real, dimension(ib - hi:ie + hi, jb - hj:je + hj, kb     :ke + hk)                :: duml ! 3d dummy variable: lower cell side
         real, dimension(ib - hi:ie + hi, jb - hj:je + hj, kb     :ke + hk)                :: dumu ! 3d dummy variable: upper cell side

         integer :: i, j, k, il, iu, jl, ju, kl, ku, n
         real    :: cf, d1, d2

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
#endif

end module advection_kappa
