!----------------------------------------------------------------------------
! This file is part of DALES.
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
! Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!----------------------------------------------------------------------------
!
!
module modstat_nc
    use typeSizes
    use netcdf
    integer :: ncid,tid,zid
    logical :: lnetcdf
contains
  subroutine initstat_nc
    use modglobal, only : kmax,ifnamopt,fname_options,iexpnr
    use modmpi,    only : mpierr,mpi_logical,comm3d,myid
    implicit none

    integer             :: status
    character(len = 20) :: ncfile
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

    if(.not.(lnetcdf)) return
    if (myid==0) then

      ncfile = 'stat.nc.xxx'
      write(ncfile(6:8),'(i3.3)') iexpnr
      write(6,*) "NETCDFSTATS: Creating: ", ncfile

      !create file
      status = nf90_create(ncfile, nf90_clobber, ncid)
      if (status /= nf90_noerr) call nchandle_error(status)

      !create dimensions
      status = nf90_def_dim(ncid, "z", kmax, zid)
      if (status /= nf90_noerr) call nchandle_error(status)
      status = nf90_def_dim(ncid, "t", nf90_unlimited, tid)
      if (status /= nf90_noerr) call nchandle_error(status)
      status = nf90_enddef(ncid)
      if (status /= nf90_noerr) call nchandle_error(status)


    end if


  return
  end subroutine initstat_nc

  subroutine exitstat_nc
    use modmpi, only : myid

    implicit none

    integer status

    if(lnetcdf .and. myid==0) then

      status = nf90_close(ncid)
      if (status /= nf90_noerr) call nchandle_error(status)
      status = nf90_close(ncid)
      if (status /= nf90_noerr) call nchandle_error(status)
    end if

  return
  end subroutine exitstat_nc


  subroutine initprof_nc(nvar,varid,varnames)
    use modglobal, only : kmax
    use modmpi,    only : myid
    implicit none
    integer :: status,n
    integer,intent(in) :: nvar
    character(len=*),dimension(nvar),intent(in) :: varnames
    integer,dimension(nvar),intent(out)       :: varid
    if (myid /=0) return
    status = nf90_redef(ncid)

    if (status /= nf90_noerr) call nchandle_error(status)
    do n=1,nvar
      status = nf90_def_var(ncid, varnames(n), nf90_float, (/zid, tid/), varid(n))
      if (status /= nf90_noerr) call nchandle_error(status)
    end do

    status = nf90_enddef(ncid)
    if (status /= nf90_noerr) call nchandle_error(status)

    return
  end subroutine initprof_nc

  subroutine writeprof_nc(nvar,varid,vars,nccall)
    use modglobal, only : kmax
    use modmpi,    only : myid
    implicit none
    integer, intent(in)             :: nvar,nccall
    real,dimension(kmax,nvar),intent(in)  :: vars
    integer,dimension(nvar),intent(in) :: varid
    integer :: status,n
    if (myid/=0) return

    do n=1,nvar
      status = nf90_put_var(ncid, varid(n), vars(1:kmax,n))
      if(status /= nf90_noerr) call nchandle_error(status)
    end do
  return
  end subroutine writeprof_nc

  subroutine inittstat_nc(nvar,varid,varnames)
    use modmpi,    only : myid
    implicit none
    integer,intent(in) :: nvar
    character(len=*),dimension(nvar),intent(in) :: varnames
    integer,dimension(nvar),intent(out)       :: varid
    integer :: status,n
    if (myid /=0) return

    status = nf90_redef(ncid)

    if (status /= nf90_noerr) call nchandle_error(status)
    do n=1,nvar
      status = nf90_def_var(ncid, varnames(n), nf90_float, tid, varid(n))
      if (status /= nf90_noerr) call nchandle_error(status)
    end do

    status = nf90_enddef(ncid)
    if (status /= nf90_noerr) call nchandle_error(status)


    return
  end subroutine inittstat_nc

  subroutine writetstat_nc(nvar,varid,vars,nccall)
    use modmpi,    only : myid
    implicit none
    integer, intent(in)             :: nvar,nccall
    real,dimension(nvar),intent(in)  :: vars
    integer,dimension(nvar),intent(in) :: varid
    integer :: status,n
    if (myid/=0) return

    do n=1,nvar
      status = nf90_put_var(ncid, varid(n), vars(n),(/nccall/))
      if(status /= nf90_noerr) call nchandle_error(status)
    end do
  return
  return
  end subroutine writetstat_nc

  subroutine nchandle_error(status)
    use typeSizes
    use netcdf
    implicit none

    integer, intent(in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if

  end subroutine nchandle_error

end module modstat_nc