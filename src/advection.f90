
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
subroutine advection
   use modglobal, only:lmoist, nsv, iadv_mom, iadv_tke, iadv_thl, iadv_qt, iadv_sv, &
      iadv_cd2, iadv_kappa, iadv_upw, &
      ltempeq, ih, jh, kh, ihc, jhc, khc, kb, ke, ib, ie, jb, je
   use modfields, only:u0, up, v0, vp, w0, wp, e120, e12p, thl0, thl0c, thlp, thlpc, qt0, qtp, sv0, svp, pres0, uh, vh, wh, pres0h
   use modsubgriddata, only:loneeqn
   use decomp_2d
#if defined(_GPU)
   use cudafor
   use cudamodule, only : griddim, blockdim, checkCUDA, &
                          u0_d, v0_d, w0_d, e120_d, thl0_d, thl0c_d, qt0_d, sv0_d, up_d, vp_d, wp_d, e12p_d, thlp_d, thlpc_d, qtp_d, svp_d, &
                          advecc_2nd_cuda, advecu_2nd_cuda, advecv_2nd_cuda, advecw_2nd_cuda, advecc_upw_cuda, &
                          advecc_kappa_reset_cuda, advecc_kappa_ducdx_cuda, advecc_kappa_dvcdy_cuda, advecc_kappa_dwcdz_cuda, advecc_kappa_add_cuda, thlptothlpc_cuda, thlpctothlp_cuda
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
            call advecc_2nd(ih, jh, kh, thl0, thlp)
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
            call advecc_kappa(ihc, jhc, khc, thl0c, thlpc)
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
            write(*,*) "Inside kappa scheme.", ih, jh, kh, ihc, jhc, khc
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
            write(*,*) "Inside upwind scheme.", ih, jh, kh, ihc, jhc, khc
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
