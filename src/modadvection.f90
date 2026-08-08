
!> \file modadvection.f90
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
         use modglobal,       only: ib, ie, jb, je, kb, ke, ih, jh, kh, ihc, jhc, khc, &
                                    iadv_mom, iadv_tke, iadv_thl, iadv_qt, iadv_sv, &
                                    iadv_cd2, iadv_kappa, iadv_upw, &
                                    ltempeq, lmoist, nsv
         use modsubgriddata,  only: loneeqn

#if defined(_GPU)
         use cudafor
         use modcuda,         only: griddim, blockdim, checkCUDA, thlptothlpc_cuda, thlpctothlp_cuda, &
                                    u0_d, up_d, v0_d, vp_d, w0_d, wp_d, e120_d, e12p_d, &
                                    thl0_d, thl0c_d, thlp_d, thlpc_d, qt0_d, qtp_d, sv0_d, svp_d
#else
         use modfields,       only: u0, up, v0, vp, w0, wp, e120, e12p, thl0, thl0c, thlp, thlpc, qt0, qtp, sv0, svp
#endif
         implicit none
         integer :: n

         select case (iadv_mom)
            case (iadv_cd2)
#if defined(_GPU)
               call advecu_2nd_cuda<<<griddim,blockdim>>>(u0_d, up_d)
               call checkCUDA( cudaGetLastError(), 'advecu_2nd_cuda' )

               call advecv_2nd_cuda<<<griddim,blockdim>>>(v0_d, vp_d)
               call checkCUDA( cudaGetLastError(), 'advecv_2nd_cuda' )

               call advecw_2nd_cuda<<<griddim,blockdim>>>(w0_d, wp_d)
               call checkCUDA( cudaGetLastError(), 'advecw_2nd_cuda' )
#else
               call advecu_2nd(u0,up)
               call advecv_2nd(v0,vp)
               call advecw_2nd(w0,wp)
#endif
            case default
               write(0, *) "ERROR: Unknown advection scheme"
               stop 1
         end select

         if (loneeqn) then
            select case (iadv_tke)
               case (iadv_cd2)
#if defined(_GPU)
                  call advecc_2nd_cuda<<<griddim,blockdim>>>(ih, jh, kh, e120_d, e12p_d)
                  call checkCUDA( cudaGetLastError(), 'advecc_2nd_cuda for e12p' )
#else
                  call advecc_2nd(ih, jh, kh, e120, e12p)
#endif
               case default
                  write(0, *) "ERROR: Unknown advection scheme"
                  stop 1
            end select
         end if

         if (ltempeq) then
            select case (iadv_thl)
               case (iadv_cd2)
#if defined(_GPU)
                  call advecc_2nd_cuda<<<griddim,blockdim>>>(ih, jh, kh, thl0_d, thlp_d)
                  call checkCUDA( cudaGetLastError(), 'advecc_2nd_cuda for thlp' )
#else
                  if (ltempeq) call advecc_2nd(ih, jh, kh, thl0, thlp)
#endif
               case (iadv_kappa)
#if defined(_GPU)
                  call thlptothlpc_cuda<<<griddim,blockdim>>>
                  call checkCUDA( cudaGetLastError(), 'thlptothlpc_cuda' )

                  ! -d(u tlh)/dx
                  call advecc_kappa_reset_cuda<<<griddim,blockdim>>>(ihc, jhc, khc)
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_reset_cuda 1st call in temp' )
                  call advecc_kappa_ducdx_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, thl0c_d)
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_ducdx_cuda in temp' )
                  call advecc_kappa_add_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, thlpc_d)
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_add_cuda 1st call in temp' )

                  ! -d(v thl)/dy
                  call advecc_kappa_reset_cuda<<<griddim,blockdim>>>(ihc, jhc, khc)
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_reset_cuda 2nd call in temp' )
                  call advecc_kappa_dvcdy_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, thl0c_d)
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_dvcdy_cuda in temp' )
                  call advecc_kappa_add_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, thlpc_d)
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_add_cuda 2nd call in temp' )

                  ! -d(w thl)/dz
                  call advecc_kappa_reset_cuda<<<griddim,blockdim>>>(ihc, jhc, khc)
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_reset_cuda 3rd call in temp' )
                  call advecc_kappa_dwcdz_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, thl0c_d)
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_dwcdz_cuda in temp' )
                  call advecc_kappa_add_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, thlpc_d)
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_add_cuda 3rd call in temp' )

                  call thlpctothlp_cuda<<<griddim,blockdim>>>
                  call checkCUDA( cudaGetLastError(), 'thlpctothlp_cuda' )
#else
                  thlpc(ib:ie,jb:je,kb:ke) = thlp(ib:ie,jb:je,kb:ke)
                  if (ltempeq) call advecc_kappa(ihc, jhc, khc, thl0c, thlpc)
                  thlp(ib:ie,jb:je,kb:ke) = thlpc(ib:ie,jb:je,kb:ke)
#endif
               case default
                  write(0, *) "ERROR: Unknown advection scheme"
                  stop 1
            end select
         end if

         if (lmoist) then
            select case (iadv_qt)
               case (iadv_cd2)
#if defined(_GPU)
                  call advecc_2nd_cuda<<<griddim,blockdim>>>(ih, jh, kh, qt0_d, qtp_d)
                  call checkCUDA( cudaGetLastError(), 'advecc_2nd_cuda for qtp' )
#else
                  call advecc_2nd(ih, jh, kh, qt0, qtp)
#endif
               case default
                  write(0, *) "ERROR: Unknown advection scheme"
                  stop 1
            end select
         end if

         do n = 1, nsv
            select case (iadv_sv (n))
               case (iadv_cd2)
#if defined(_GPU)
                  call advecc_2nd_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, sv0_d(:, :, :, n), svp_d(:, :, :, n))
                  call checkCUDA( cudaGetLastError(), 'advecc_2nd_cuda for svp' )
#else
                  call advecc_2nd(ihc, jhc, khc, sv0(:, :, :, n), svp(:, :, :, n))
#endif
               case (iadv_kappa)
#if defined(_GPU)
                  ! -d(uc)/dx
                  call advecc_kappa_reset_cuda<<<griddim,blockdim>>>(ihc, jhc, khc)
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_reset_cuda 1st call in scalar' )
                  call advecc_kappa_ducdx_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, sv0_d(:, :, :, n))
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_ducdx_cuda in scalar' )
                  call advecc_kappa_add_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, svp_d(:, :, :, n))
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_add_cuda 1st call in scalar' )

                  ! -d(vc)/dy
                  call advecc_kappa_reset_cuda<<<griddim,blockdim>>>(ihc, jhc, khc)
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_reset_cuda 2nd call in scalar' )
                  call advecc_kappa_dvcdy_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, sv0_d(:, :, :, n))
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_dvcdy_cuda in scalar' )
                  call advecc_kappa_add_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, svp_d(:, :, :, n))
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_add_cuda 2nd call in scalar' )

                  ! -d(wc)/dz
                  call advecc_kappa_reset_cuda<<<griddim,blockdim>>>(ihc, jhc, khc)
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_reset_cuda 3rd call in scalar' )
                  call advecc_kappa_dwcdz_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, sv0_d(:, :, :, n))
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_dwcdz_cuda in scalar' )
                  call advecc_kappa_add_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, svp_d(:, :, :, n))
                  call checkCUDA( cudaGetLastError(), 'advecc_kappa_add_cuda 3rd call in scalar' )
#else
                  call advecc_kappa(ihc, jhc, khc, sv0(:, :, :, n), svp(:, :, :, n))
#endif
               case (iadv_upw)
#if defined(_GPU)
                  call advecc_upw_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, sv0_d(:, :, :, n), svp_d(:, :, :, n))
                  call checkCUDA( cudaGetLastError(), 'advecc_upw_cuda for svp' )
#else
                  call advecc_upw(ihc, jhc, khc, sv0(:, :, :, n), svp(:, :, :, n))
#endif
               case default
                  write(0, *) "ERROR: Unknown advection scheme"
                  stop 1
            end select
         end do
      end subroutine advection

#if !defined(_GPU)

      !> Advection at cell center central difference.
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
               end do
            end do
         end do

         do j = jb, je
            jm = j - 1
            jp = j + 1
            do i = ib, ie
               im = i - 1
               ip = i + 1
               do k = kb, ke
                  km = k - 1
                  kp = k + 1
                  putout(i, j, k) = putout(i, j, k) - ( &
                                    w0(i, j, kp)*(putin(i, j, kp)*dzf(k) + putin(i, j, k)*dzf(kp))*dzhi(kp) &
                                 - w0(i, j, k)*(putin(i, j, km)*dzf(k) + putin(i, j, k)*dzf(km))*dzhi(k) &
                                    )*dzfi5(k)
               end do
            end do
         end do
      end subroutine advecc_2nd

      !> Advection at the u point central difference.
      subroutine advecu_2nd(putin, putout)
         use modglobal, only:ih, ib, ie, jb, je, jh, kb, ke, kh, dxi, dxiq, dyiq, dzf, dzfi5, dzhi
         use modfields, only:u0, v0, w0, pres0
         implicit none

         real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh), intent(in)  :: putin !< Input: the u-field
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
                                    (putin(i, j, k) + putin(ip, j, k))*(u0(i, j, k) + u0(ip, j, k)) &
                                 - (putin(i, j, k) + putin(im, j, k))*(u0(i, j, k) + u0(im, j, k)) & ! d(uu)/dx
                                    )*dxiq &
                                 + ( &
                                    (putin(i, j, k) + putin(i, jp, k))*(v0(i, jp, k) + v0(im, jp, k)) &
                                 - (putin(i, j, k) + putin(i, jm, k))*(v0(i, j, k) +  v0(im, j, k)) & ! d(vu)/dy
                                    )*dyiq) &
                                 - ((pres0(i, j, k) - pres0(i - 1, j, k))*dxi) ! - dp/dx

               end do
            end do
         end do

         do j = jb, je
            jm = j - 1
            jp = j + 1
            do i = ib, ie
               im = i - 1
               ip = i + 1
               do k = kb, ke
                  km = k - 1
                  kp = k + 1
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

      !> Advection at the v point central difference.
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

               end do
            end do
         end do

         do j = jb, je
            jm = j - 1
            jp = j + 1
            do i = ib, ie
               im = i - 1
               ip = i + 1
               do k = kb, ke
                  km = k - 1
                  kp = k + 1
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

      !> Advection at the w point central difference.
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
         use modglobal, only: ib, ie, jb, je, kb, ke, dxhci, dyi, dzhci, dxfc, dzfc, dxfci, dzfci
         use modfields, only: u0, v0, w0
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

#else

      ! CUDA implementations for advec_2nd
      attributes(global) subroutine advecc_2nd_cuda(hi, hj, hk, putin, putout)
         use modcuda, only: ib_d, ie_d, jb_d, je_d, kb_d, ke_d, dxi5_d, dyi5_d, dzf_d, dzfi5_d, dzhi_d, &
                            u0_d, v0_d, w0_d, &
                            tidandstride
         implicit none

         integer, value, intent(in) :: hi, hj, hk
         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d - hk:ke_d + hk), intent(in)    :: putin
         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d     :ke_d + hk), intent(inout) :: putout

         integer :: tidx, tidy, tidz, stridex, stridey, stridez
         integer :: i, j, k, im, ip, jm, jp, km, kp

         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

         do i = tidx, ie_d, stridex
            im = i - 1
            ip = i + 1
            do j = tidy, je_d, stridey
               jm = j - 1
               jp = j + 1
               do k = tidz, ke_d, stridez
                  putout(i, j, k) = putout(i, j, k) - ( &
                                    ( &
                                    u0_d(ip, j, k)*(putin(ip, j, k) + putin(i, j, k)) &
                                  - u0_d(i, j, k)*(putin(im, j, k) + putin(i, j, k)) & ! d(uc)/dx
                                    )*dxi5_d &
                                  + ( & !
                                    v0_d(i, jp, k)*(putin(i, jp, k) + putin(i, j, k)) &
                                  - v0_d(i, j, k)*(putin(i, jm, k) + putin(i, j, k)) & ! d(vc)/dy
                                    )*dyi5_d)
               end do
            end do
         end do

         do i = tidx, ie_d, stridex
            im = i - 1
            ip = i + 1
            do j = tidy, je_d, stridey
               jm = j - 1
               jp = j + 1
               do k = tidz, ke_d, stridez
                  km = k - 1
                  kp = k + 1
                  putout(i, j, k) = putout(i, j, k) - ( &
                                    w0_d(i, j, kp)*(putin(i, j, kp)*dzf_d(k) + putin(i, j, k)*dzf_d(kp))*dzhi_d(kp) &
                                  - w0_d(i, j, k)*(putin(i, j, km)*dzf_d(k) + putin(i, j, k)*dzf_d(km))*dzhi_d(k) &
                                    )*dzfi5_d(k)
               end do
            end do
         end do
      end subroutine advecc_2nd_cuda
      attributes(global) subroutine advecu_2nd_cuda(putin, putout)
         use modcuda, only: ib_d, ie_d, ih_d, jb_d, je_d, jh_d, kb_d, ke_d, kh_d, dxi_d, dxiq_d, dyiq_d, dzf_d, dzfi5_d, dzhi_d, &
                            u0_d, v0_d, w0_d, pres0_d, &
                            tidandstride
         implicit none

         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d - kh_d:ke_d + kh_d), intent(in)    :: putin
         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d       :ke_d + kh_d), intent(inout) :: putout

         integer :: tidx, tidy, tidz, stridex, stridey, stridez
         integer :: i, j, k, im, ip, jm, jp, km, kp

         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

         do i = tidx, ie_d, stridex
            im = i - 1
            ip = i + 1
            do j = tidy, je_d, stridey
               jm = j - 1
               jp = j + 1
               do k = tidz, ke_d, stridez
                  km = k - 1
                  kp = k + 1
                  putout(i, j, k) = putout(i, j, k) - ( &
                                    ( &
                                    (putin(i, j, k) + putin(ip, j, k))*(u0_d(i, j, k) + u0_d(ip, j, k)) &
                                  - (putin(i, j, k) + putin(im, j, k))*(u0_d(i, j, k) + u0_d(im, j, k)) & ! d(uu)/dx
                                    )*dxiq_d &
                                  + ( &
                                    (putin(i, j, k) + putin(i, jp, k))*(v0_d(i, jp, k) + v0_d(im, jp, k)) &
                                  - (putin(i, j, k) + putin(i, jm, k))*(v0_d(i, j, k) + v0_d(im, j, k)) & ! d(vu)/dy
                                    )*dyiq_d) &
                                  - ((pres0_d(i, j, k) - pres0_d(im, j, k))*dxi_d) ! - dp/dx
               end do
            end do
         end do

         do k = tidz, ke_d, stridez
            km = k - 1
            kp = k + 1
            do j = tidy, je_d, stridey
               jm = j - 1
               jp = j + 1
               do i = tidx, ie_d, stridex
                  im = i - 1
                  ip = i + 1
                  putout(i, j, k) = putout(i, j, k) - ( &
                                    (putin(i, j, kp)*dzf_d(k) + putin(i, j, k)*dzf_d(kp))*dzhi_d(kp) &
                                  * (w0_d(i, j, kp) + w0_d(im, j, kp)) &
                                  - (putin(i, j, k)*dzf_d(km) + putin(i, j, km)*dzf_d(k))*dzhi_d(k) &
                                  * (w0_d(i, j, k) + w0_d(im, j, k)) &
                                    )*0.5*dzfi5_d(k)
               end do
            end do
         end do
      end subroutine advecu_2nd_cuda
      attributes(global) subroutine advecv_2nd_cuda(putin, putout)
         use modcuda, only: ib_d, ie_d, ih_d, jb_d, je_d, jh_d, kb_d, ke_d, kh_d, dxiq_d, dyi_d, dyiq_d, dzf_d, dzfi5_d, dzhi_d, &
                            u0_d, v0_d, w0_d, pres0_d, &
                            tidandstride
         implicit none

         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d - kh_d:ke_d + kh_d), intent(in)    :: putin
         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d       :ke_d + kh_d), intent(inout) :: putout

         integer :: tidx, tidy, tidz, stridex, stridey, stridez
         integer :: i, j, k, im, ip, jm, jp, km, kp

         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

         do i = tidx, ie_d, stridex
            im = i - 1
            ip = i + 1
            do j = tidy, je_d, stridey
               jm = j - 1
               jp = j + 1
               do k = tidz, ke_d, stridez
                  km = k - 1
                  kp = k + 1
                  putout(i, j, k) = putout(i, j, k) - ( &
                                    ( &
                                    (u0_d(ip, j, k) + u0_d(ip, jm, k))*(putin(i, j, k) + putin(ip, j, k)) &
                                  - (u0_d(i, j, k)  + u0_d(i, jm, k)) *(putin(i, j, k) + putin(im, j, k)) & ! d(uv)/dx
                                    )*dxiq_d &
                                  + ( &
                                    (v0_d(i, jp, k) + v0_d(i, j, k))*(putin(i, j, k) + putin(i, jp, k)) &
                                  - (v0_d(i, jm, k) + v0_d(i, j, k))*(putin(i, j, k) + putin(i, jm, k)) & ! d(vv)/dy
                                    )*dyiq_d &
                                    ) &
                                  - ((pres0_d(i, j, k) - pres0_d(i, jm, k))*dyi_d) ! - dp/dy
               end do
            end do
         end do

         do k = tidz, ke_d, stridez
            km = k - 1
            kp = k + 1
            do j = tidy, je_d, stridey
               jm = j - 1
               jp = j + 1
               do i = tidx, ie_d, stridex
                  im = i - 1
                  ip = i + 1
                  putout(i, j, k) = putout(i, j, k) - ( &
                                    (w0_d(i, j, kp) + w0_d(i, jm, kp)) &
                                  * (putin(i, j, kp)*dzf_d(k) + putin(i, j, k)*dzf_d(kp))*dzhi_d(kp) &
                                  - (w0_d(i, j, k) + w0_d(i, jm, k)) &
                                  * (putin(i, j, km)*dzf_d(k) + putin(i, j, k)*dzf_d(km))*dzhi_d(k) &
                                    )*0.5*dzfi5_d(k)
               end do
            end do
         end do
      end subroutine advecv_2nd_cuda
      attributes(global) subroutine advecw_2nd_cuda(putin, putout)
         use modcuda, only: ib_d, ie_d, ih_d, jb_d, je_d, jh_d, kb_d, ke_d, kh_d, dxiq_d, dyiq_d, dzf_d, dzhi_d, dzhiq_d, &
                            u0_d, v0_d, w0_d, pres0_d, &
                            tidandstride
         implicit none

         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d - kh_d:ke_d + kh_d), intent(in)    :: putin
         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d       :ke_d + kh_d), intent(inout) :: putout

         integer :: tidx, tidy, tidz, stridex, stridey, stridez
         integer :: i, j, k, im, ip, jm, jp, km, kp

         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

         do i = tidx, ie_d, stridex
            im = i - 1
            ip = i + 1
            do j = tidy, je_d, stridey
               jm = j - 1
               jp = j + 1
               do k = tidz, ke_d, stridez
                  km = k - 1
                  kp = k + 1
                  putout(i, j, k) = putout(i, j, k) - ( &
                                    ( &
                                    (putin(ip, j, k) + putin(i, j, k))*(dzf_d(km)*u0_d(ip, j, k) + dzf_d(k)*u0_d(ip, j, km)) &
                                  - (putin(i, j, k)  + putin(im, j, k))*(dzf_d(km)*u0_d(i, j, k) + dzf_d(k)*u0_d(i, j, km)) &
                                    )*dxiq_d*dzhi_d(k) & ! d(uw)/dx
                                  + ( &
                                    (putin(i, jp, k) + putin(i, j, k))*(dzf_d(km)*v0_d(i, jp, k) + dzf_d(k)*v0_d(i, jp, km)) &
                                  - (putin(i, j, k) + putin(i, jm, k))*(dzf_d(km)*v0_d(i, j, k) + dzf_d(k)*v0_d(i, j, km)) &
                                    )*dyiq_d*dzhi_d(k) & ! d(vw)/dy
                                  + ( &
                                    (putin(i, j, k) + putin(i, j, kp))*(w0_d(i, j, k) + w0_d(i, j, kp)) &
                                  - (putin(i, j, k) + putin(i, j, km))*(w0_d(i, j, k) + w0_d(i, j, km)) &
                                    )*dzhiq_d(k) & ! d(ww)/dz
                                    ) &
                                  - ((pres0_d(i, j, k) - pres0_d(i, j, km))*dzhi_d(k)) ! - dp/dz
               end do
            end do
         end do
      end subroutine advecw_2nd_cuda

      ! CUDA implementations for advec_kappa
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

      ! CUDA implementation for advec_upw
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

#endif

end module modadvection
