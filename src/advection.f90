!> \file advection.f90
!!  Advection management

!>
!!  Advection management
!! \par Revision list
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

  use modglobal, only : lmoist, nsv, iadv_mom,iadv_tke,iadv_thl,iadv_qt,iadv_sv, &
                        iadv_cd2,iadv_5th,iadv_cd6,iadv_kappa,iadv_upw
  use modfields, only : u0,up,v0,vp,w0,wp,e120,e12p,thl0,thlp,qt0,qtp,sv0,svp
  use modsubgrid, only: lsmagorinsky
  implicit none
  integer :: n

  select case(iadv_mom)
    case(iadv_cd2)
      call advecu_2nd(u0,up)
      call advecv_2nd(v0,vp)
      call advecw_2nd(w0,wp)
    case(iadv_5th)
      call advecu_5th(u0,up)
      call advecv_5th(v0,vp)
      call advecw_5th(w0,wp)
    case(iadv_cd6)
      call advecu_6th(u0,up)
      call advecv_6th(v0,vp)
      call advecw_6th(w0,wp)
    case default
      stop "Unknown advection scheme "
  end select

  if (.not. lsmagorinsky) then
    select case(iadv_tke)
      case(iadv_cd2)
        call advecc_2nd(e120,e12p)
      case(iadv_5th)
        call advecc_5th(e120,e12p)
      case(iadv_cd6)
        call advecc_6th(e120,e12p)
      case(iadv_kappa)
        call advecc_kappa(e120,e12p)
      case default
        stop "Unknown advection scheme "
    end select
  end if

  select case(iadv_thl)
    case(iadv_cd2)
      call advecc_2nd(thl0,thlp)
    case(iadv_5th)
      call advecc_5th(thl0,thlp)
    case(iadv_cd6)
      call advecc_6th(thl0,thlp)
    case(iadv_kappa)
      call advecc_kappa(thl0,thlp)
    case(iadv_upw)
      call advecc_upw(thl0,thlp)
    case default
      stop "Unknown advection scheme "
  end select
  if (lmoist) then
    select case(iadv_qt)
      case(iadv_cd2)
        call advecc_2nd(qt0,qtp)
      case(iadv_5th)
        call advecc_5th(qt0,qtp)
      case(iadv_cd6)
        call advecc_6th(qt0,qtp)
      case(iadv_kappa)
        call advecc_kappa(qt0,qtp)
      case(iadv_upw)
        call advecc_upw(qt0,qtp)
      case default
        stop "Unknown advection scheme "
    end select
  end if
  do n=1,nsv
    select case(iadv_sv(n))
    case(iadv_cd2)
      call advecc_2nd(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_5th)
      call advecc_5th(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_cd6)
      call advecc_6th(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_kappa)
      call advecc_kappa(sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_upw)
      call advecc_upw(sv0(:,:,:,n),svp(:,:,:,n))
    case default
      stop "Unknown advection scheme "
    end select
  end do

end subroutine advection
