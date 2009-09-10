!> \file stat_nc_dummy.f90
!!  Replacement file for the netcdf routine.

!>
!!  Replacement file for the netcdf routine.
!>
!! If netcdf is not desired, copy this file to the source dir, and remove modstat_nc.f90
!!  \author Thijs Heus,MPI-M
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
module modstat_nc
    implicit none
    logical :: lnetcdf
    integer, save :: timeID=0, ztID=0, zmID=0, xtID=0, xmID=0, ytID=0, ymID=0
!> The only interface necessary to write data to netcdf, regardless of the dimensions.
    interface writestat_nc
      module procedure writestat_time_nc
      module procedure writestat_1D_nc
      module procedure writestat_2D_nc
      module procedure writestat_3D_nc
      module procedure writestat_3D_short_nc
    end interface writestat_nc
contains


  subroutine initstat_nc
    use modglobal, only : kmax,ifnamopt,fname_options,iexpnr
    use modmpi,    only : mpierr,mpi_logical,comm3d,myid
    implicit none

    integer             :: ierr

    namelist/NAMNETCDFSTATS/ &
    lnetcdf

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMNETCDFSTATS,iostat=ierr)
      write(6, NAMNETCDFSTATS)
      close(ifnamopt)
    end if

    call MPI_BCAST(lnetcdf    ,1,MPI_LOGICAL, 0,comm3d,mpierr)
    if(lnetcdf) stop'STOP: lnetcdf is true, but not supported by the DALES build.'

    
  end subroutine initstat_nc
!
! ----------------------------------------------------------------------
!> Subroutine Open_NC: Opens a NetCDF File and identifies starting record
!
  subroutine open_nc (fname, ncid,n1, n2, n3)
    use modglobal, only : author,version
    implicit none
    integer, intent (out) :: ncid
    integer, optional, intent (in) :: n1, n2, n3
    character (len=40), intent (in) :: fname

  end subroutine open_nc

  !
  ! ----------------------------------------------------------------------
  !> Subroutine Define_NC: Defines the structure of the nc file (if not
  !! already open)
  !
  subroutine define_nc(ncID, nVar, sx)
    implicit none
    integer, intent (in) :: nVar, ncID
    character (*), intent (in) :: sx(nVar,4)

  end subroutine define_nc

  subroutine redefine_nc(ncid)
  implicit none
    integer, intent(in) :: ncid
  
  end subroutine redefine_nc

  subroutine exitstat_nc(ncid)

   implicit none
   integer, intent(in) :: ncid
 
  end subroutine exitstat_nc
  subroutine writestat_dims_nc(ncid)
    use modglobal, only : dx,dy,zf,zh
    implicit none
    integer, intent(in) :: ncid
   
  end subroutine writestat_dims_nc
  subroutine writestat_time_nc(ncid,nvar,ncname,vars,nrec,lraise)
    implicit none
    integer, intent(in)                      :: ncid,nvar
    integer, intent(inout)                   :: nrec
    real,dimension(nvar),intent(in)          :: vars
    character(*), dimension(:,:),intent(in)  :: ncname
    logical, intent(in)                      :: lraise


  end subroutine writestat_time_nc

  subroutine writestat_1D_nc(ncid,nvar,ncname,vars,nrec,dim1)
    implicit none
    integer, intent(in)                      :: ncid,nvar,dim1
    integer, intent(in)                      :: nrec
    real,dimension(dim1,nvar),intent(in)     :: vars
    character(*), dimension(:,:),intent(in)  :: ncname


  end subroutine writestat_1D_nc

  subroutine writestat_2D_nc(ncid,nvar,ncname,vars,nrec,dim1,dim2)
    implicit none
    integer, intent(in)                      :: ncid,nvar,dim1,dim2
    integer, intent(in)                      :: nrec
    real,dimension(:,:,:),intent(in)         :: vars
    character(*), dimension(:,:),intent(in)  :: ncname

  end subroutine writestat_2D_nc
  subroutine writestat_3D_nc(ncid,nvar,ncname,vars,nrec,dim1,dim2,dim3)
    implicit none
    integer, intent(in)                      :: ncid,nvar,dim1,dim2,dim3
    integer, intent(in)                      :: nrec
    real,dimension(dim1,dim2,dim3,nvar),intent(in)       :: vars
    character(*), dimension(:,:),intent(in)  :: ncname

  end subroutine writestat_3D_nc
  subroutine writestat_3D_short_nc(ncid,nvar,ncname,vars,nrec,dim1,dim2,dim3)
    implicit none
    integer, intent(in)                      :: ncid,nvar,dim1,dim2,dim3
    integer, intent(in)                      :: nrec
    integer(KIND=selected_int_kind(4)),dimension(dim1,dim2,dim3,nvar),intent(in)       :: vars
    character(*), dimension(:,:),intent(in)  :: ncname

  end subroutine writestat_3D_short_nc


  subroutine ncinfo(out,in1,in2,in3,in4)

    implicit none
    character(*), dimension(4),intent(out) ::out
    character(*), intent(in) ::in1,in2,in3,in4
   
  end subroutine ncinfo

end module modstat_nc