!> \file advection.f90
!!  Advection management
!>
!!  Advection management
!! \par Revision list
!! variable x-grid now possible 
!! \par Authors
!! Jasper Tomas, 31 March 2014
!! Thijs Heus, Chiel van Heerwaarden, 15 June 2007

!> Advection redirection function
subroutine advection

  use modglobal,  only : lmoist, nsv, iadv_mom,iadv_tke,iadv_thl,iadv_qt,iadv_sv,&
                         iadv_cd2,iadv_5th,iadv_52,iadv_cd6,iadv_62,iadv_kappa,iadv_upw,&
                         ltempeq,ih,jh,kh,ihc,jhc,khc
  use modfields,  only : u0,up,v0,vp,w0,wp,e120,e12p,thl0,thlp,qt0,qtp,sv0,svp
  use modsubgriddata, only : lsmagorinsky,lvreman,loneeqn
  implicit none
  integer :: n

  select case(iadv_mom)
    case(iadv_cd2)
      call advecu_2nd(u0,up)
      call advecv_2nd(v0,vp)
      call advecw_2nd(w0,wp)
    case default
      stop "Unknown advection scheme "
  end select

  if (loneeqn == .true.) then
    select case(iadv_tke)
      case(iadv_cd2)
        call advecc_2nd(ih,jh,kh,e120,e12p)
      case default
        stop "Unknown advection scheme "
    end select
  end if

  select case(iadv_thl)
    case(iadv_cd2)
      if (ltempeq) call advecc_2nd(ih,jh,kh,thl0,thlp)
    case default
      stop "Unknown advection scheme "
  end select
  if (lmoist) then
    select case(iadv_qt)
      case(iadv_cd2)
        call advecc_2nd(ih,jh,kh,qt0,qtp)
      case default
        stop "Unknown advection scheme "
    end select
  end if
  do n=1,nsv
    select case(iadv_sv(n))
    case(iadv_cd2)
      call advecc_2nd(ihc,jhc,khc,sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_kappa)
!      call advecc_kappa(sv0(:,:,:,n),svp(:,:,:,n))
      call advecc_kappa(ihc,jhc,khc,sv0(:,:,:,n),svp(:,:,:,n))
    case(iadv_upw)
!      call advecc_kappa(sv0(:,:,:,n),svp(:,:,:,n))
      call advecc_upw(ihc,jhc,khc,sv0(:,:,:,n),svp(:,:,:,n))
    case default
      stop "Unknown advection scheme "
    end select
  end do

end subroutine advection
