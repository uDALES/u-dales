!> \file modmicrophysics.f90
!!  Microphysics abstraction layer.

!>
!!  Microphysics abstraction layer.
!>
!!  Also provides the drizzle routine
!!  \author Hans Cuijpers, IMAU
!!  \author Thijs Heus,MPI-M
!!  \todo Documentation
!!  \par Revision list
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



module modmicrophysics


implicit none
private
public :: initmicrophysics, microphysics, microsources,exitmicrophysics
public :: imicro,imicro_none,imicro_drizzle,imicro_bulk,imicro_bin

integer :: imicro = 0

integer, parameter :: imicro_none    = 0
integer, parameter :: imicro_drizzle = 1
integer, parameter :: imicro_bulk    = 2
integer, parameter :: imicro_bin     = 3
integer, parameter :: imicro_user    = 10

contains
  subroutine initmicrophysics
    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_integer,mpi_real,mpi_logical
    use modglobal,only :ifnamopt,fname_options
    use modbulkmicro, only : initbulkmicro
    use modbulkmicrodata, only : l_sb,l_rain,l_sedc,l_mur_cst,mur_cst,& ! OG
                          Nc_0, sig_g, sig_gr                           ! SdeR
!     use modbinmicro,  only : initbinmicro
    implicit none
    integer :: ierr
    namelist/NAMMICROPHYSICS/ &
    imicro,l_sb,l_rain,l_sedc,l_mur_cst,mur_cst, &     ! OG
    Nc_0, sig_g, sig_gr                                ! SdeR

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMMICROPHYSICS,iostat=ierr)
      write(6 ,NAMMICROPHYSICS)
      close(ifnamopt)
    end if

    call MPI_BCAST(imicro,   1, MPI_INTEGER ,0,comm3d,ierr)
    call MPI_BCAST(l_sb,     1, MPI_LOGICAL ,0,comm3d,ierr)
    call MPI_BCAST(l_rain,   1, MPI_LOGICAL ,0,comm3d,ierr)
    call MPI_BCAST(l_sedc,   1, MPI_LOGICAL ,0,comm3d,ierr)
    call MPI_BCAST(l_mur_cst,1, MPI_INTEGER ,0,comm3d,ierr)
    call MPI_BCAST(mur_cst,  1, MPI_REAL    ,0,comm3d,ierr)
    call MPI_BCAST(Nc_0  ,   1, MPI_LOGICAL ,0,comm3d,ierr)
    call MPI_BCAST(sig_g,    1, MPI_INTEGER ,0,comm3d,ierr)
    call MPI_BCAST(sig_gr,   1, MPI_REAL    ,0,comm3d,ierr)

    select case (imicro)
    case(imicro_none)
    case(imicro_drizzle)
    case(imicro_bulk)
      call initbulkmicro
    case(imicro_bin)
!       call initbinmicro
    case(imicro_user)
    end select
  end subroutine initmicrophysics


  subroutine microphysics
!     use modbulkmicro, only : bulkmicro
!     use modbinmicro,  only : binmicro
    implicit none
    select case (imicro)
    case(imicro_none)
    case(imicro_drizzle)
    case(imicro_bulk)
!       call bulkmicro
    case(imicro_bin)
!       call binmicro
    case(imicro_user)
    end select
  end subroutine microphysics

  subroutine microsources
   use modbulkmicro, only : bulkmicro
!     use modbinmicro,  only : binmicrosources
    implicit none

    select case (imicro)
    case(imicro_none)
    case(imicro_drizzle)
      call drizzle
    case(imicro_bulk)
      call bulkmicro
    case(imicro_bin)
!       call binmicrosources
    case(imicro_user)
      call micro_user
    end select

  end subroutine microsources
  subroutine exitmicrophysics
    use modbulkmicro, only : exitbulkmicro
 !     use modbinmicro,  only : exitbinmicro
    implicit none

     select case (imicro)
     case(imicro_none)
     case(imicro_drizzle)
     case(imicro_bulk)
       call exitbulkmicro
     case(imicro_bin)
 !       call exitbinmicro
      case(imicro_user)
     end select
  end subroutine exitmicrophysics
 subroutine drizzle

!-----------------------------------------------------------------|
!                                                                 |
!      Hans Cuijpers   I.M.A.U.  23 May 1995                      |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!      Calculates gravitational settling (or rainfall rate)       |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!     *drizzle* is called from *program*.                         |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal, only : i1,j1,kmax,rlv,cp,dzf,pi
  use modfields, only : qtp,ql0,thlp,rhof,exnf
  use modbulkmicrodata, only : rhow, sig_g, Nc_0,c_St
  implicit none
  real :: sedc,csed
  integer :: i, j, k
    csed = c_St*(3./(4.*pi*rhow))**(2./3.)*exp(5.*log(sig_g)**2.)*Nc_0**(-2./3.)
  sedc = 0.

  do k = 1,kmax
  do i=2,i1
  do j=2,j1
  if (ql0(i,j,k)>0.0) then
    sedc= csed*((ql0(i,j,k+1)*rhof(k+1))**(5./3.)-(ql0(i,j,k)*rhof(k))**(5./3.))/(dzf(k)*rhof(k))
    qtp(i,j,k) = qtp(i,j,k) + sedc
    thlp(i,j,k) = thlp(i,j,k) - (rlv/(cp*exnf(k)))*sedc
  endif
  enddo
  enddo
  enddo
  return
  end subroutine drizzle

end module modmicrophysics
