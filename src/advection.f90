
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
   use mpi
   use modglobal, only:lmoist, nsv, iadv_mom, iadv_tke, iadv_thl, iadv_qt, iadv_sv, &
      iadv_cd2, iadv_kappa, iadv_upw, &
      ltempeq, ih, jh, kh, ihc, jhc, khc, kb, ke, ib, ie, jb, je
   use modfields, only:u0, up, v0, vp, w0, wp, e120, e12p, thl0, thl0c, thlp, thlpc, qt0, qtp, sv0, svp, pres0, uh, vh, wh, pres0h
   use modsubgriddata, only:loneeqn
   use decomp_2d
#if defined(_GPU)
   use cudafor
   use cudamodule, only : griddim, blockdim, checkCUDA, &
                          u0_d, v0_d, w0_d, e120_d, thl0_d, qt0_d, sv0_d, up_d, vp_d, wp_d, e12p_d, thlp_d, qtp_d, svp_d, &
                          advecc_2nd_cuda, advecu_2nd_cuda, advecv_2nd_cuda, advecw_2nd_cuda
#endif
   implicit none
   integer :: n
   real    :: stime

   stime = MPI_Wtime()

  select case (iadv_mom)
     case (iadv_cd2)
#if defined(_GPU)

 !        stime = MPI_Wtime()
         call advecu_2nd_cuda<<<griddim,blockdim>>>(u0_d, up_d)
         call checkCUDA( cudaGetLastError(), 'advecu_2nd_cuda' )
  !       write(6,*)'advecu_2nd time = ', MPI_Wtime() - stime

 !        stime = MPI_Wtime()
         call advecv_2nd_cuda<<<griddim,blockdim>>>(v0_d, vp_d)
         call checkCUDA( cudaGetLastError(), 'advecv_2nd_cuda' )
 !       write(6,*)'advecv_2nd time = ', MPI_Wtime() - stime

 !        stime = MPI_Wtime()
         call advecw_2nd_cuda<<<griddim,blockdim>>>(w0_d, wp_d)
         call checkCUDA( cudaGetLastError(), 'advecw_2nd_cuda' )
 !       write(6,*)'advecw_2nd time = ', MPI_Wtime() - stime

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
   !         stime = MPI_Wtime()
            call advecc_2nd_cuda<<<griddim,blockdim>>>(ih, jh, kh, e120_d, e12p_d)
            call checkCUDA( cudaGetLastError(), 'advecc_2nd_cuda for e12p' )
   !         write(6,*)'advecc_2nd time for e12p = ', MPI_Wtime() - stime
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
    !        stime = MPI_Wtime()
            call advecc_2nd_cuda<<<griddim,blockdim>>>(ih, jh, kh, thl0_d, thlp_d)
            call checkCUDA( cudaGetLastError(), 'advecc_2nd_cuda for thlp' )
     !       write(6,*)'advecc_2nd GPU time for thlp = ', MPI_Wtime() - stime
#else
            call advecc_2nd(ih, jh, kh, thl0, thlp)
#endif
         case (iadv_kappa)
            thlpc(ib:ie,jb:je,kb:ke) = thlp(ib:ie,jb:je,kb:ke)
            call advecc_kappa(ihc, jhc, khc, thl0c, thlpc)
            thlp(ib:ie,jb:je,kb:ke) = thlpc(ib:ie,jb:je,kb:ke)
         case default
            write(0, *) "ERROR: Unknown advection scheme"
            stop 1
      end select
   end if

   if (lmoist) then
      select case (iadv_qt)
         case (iadv_cd2)
#if defined(_GPU)
     !       stime = MPI_Wtime()
            call advecc_2nd_cuda<<<griddim,blockdim>>>(ih, jh, kh, qt0_d, qtp_d)
            call checkCUDA( cudaGetLastError(), 'advecc_2nd_cuda for qtp' )
     !       write(6,*)'advecc_2nd GPU time for qtp = ', MPI_Wtime() - stime
#else
            call advecc_2nd(ih, jh, kh, qt0, qtp)
#endif
         case default
            write(0, *) "ERROR: Unknown advection scheme"
            stop 1
      end select
   end if

   write(6,*)'advecc_2nd time = ', MPI_Wtime() - stime

   do n = 1, nsv
      select case (iadv_sv (n))
         case (iadv_cd2)
#if defined(_GPU)
      !      stime = MPI_Wtime()
            call advecc_2nd_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, sv0_d(:, :, :, n), svp_d(:, :, :, n))
            call checkCUDA( cudaGetLastError(), 'advecc_2nd_cuda for svp' )
      !      write(6,*)'advecc_2nd GPU time for svp = ', MPI_Wtime() - stime
#else
            call advecc_2nd(ihc, jhc, khc, sv0(:, :, :, n), svp(:, :, :, n))
#endif
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
