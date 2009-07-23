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
module modnetcdfmovie

use modglobal

implicit none
private
public :: initnetcdfmovie, netcdfmovie, exitnetcdfmovie

save
  !namelist variables
  real                 :: dtmovie       = 1.e6
  real                 :: tnext
  logical              :: lnetcdfmovie  = .false.
  integer              :: ncklimit      = 60
  logical              :: lmoviez        = .true.
  integer              :: slicex        = 10
  integer              :: slicez        = 10

  integer              :: ndtmovie


  !VARIABLES FOR MOVIE
  !id of netcdf
  integer :: ncid

  !id of dimensions
  integer :: xid, yid, zid, tid

  !id of variables
  integer :: uid, vid, wid, thlid, qtid


  !COMMON VARIABLES
  !constants
  !integer ncklimit
  !parameter (ncklimit = 50)
  integer :: nccall = 1

  !integer :: ncfieldflag
  !parameter (ncfieldflag = 31)
  ! 1 = u
  ! 2 = v
  ! 4 = w
  ! 8 = thl
  ! 16 = qt

contains
  subroutine initnetcdfmovie

!     use typeSizes
    use netcdf, only : nf90_def_var,nf90_unlimited,nf90_enddef,nf90_float, nf90_clobber, nf90_create, nf90_def_dim, nf90_noerr
    use modmpi, only : myid,comm3d,mpierr, my_real,mpi_logical, mpi_integer
    use modglobal, only : dtmax,ladaptive, dt_lim

    implicit none

    integer             :: status, npes
    character(len = 20) :: ncfile
    integer             :: n, ierr
    character(20)       :: name

    namelist/NAMNETCDFMOVIE/ &
    dtmovie,lnetcdfmovie, ncklimit, lmoviez, slicex, slicez

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMNETCDFMOVIE,iostat=ierr)
      write(6, NAMNETCDFMOVIE)
      close(ifnamopt)
    end if

    call MPI_BCAST(dtmovie         ,1,MY_REAL,      0, comm3d, mpierr)
    call MPI_BCAST(lnetcdfmovie    ,1,MPI_LOGICAL,  0, comm3d, mpierr)
    call MPI_BCAST(ncklimit        ,1,MPI_INTEGER,  0, comm3d, mpierr)
    call MPI_BCAST(lmoviez         ,1,MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(slicex          ,1,MPI_INTEGER,  0, comm3d, mpierr)
    call MPI_BCAST(slicez          ,1,MPI_INTEGER,  0, comm3d, mpierr)
    tnext = dtmovie-1e-6



    if(.not.(lnetcdfmovie)) return
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtmovie/dtmax-nint(dtmovie/dtmax))>1e-4) then
      stop 'NETCDF: dtav should be a integer multiple of dtmax'
    end if

    ncfile = 'data/movie123.nc'
    write(ncfile(11:13),'(i3.3)') myid
    write(6,*) "NETCDFMOVIE: Creating: ", ncfile

    !create file
    status = nf90_create(ncfile, nf90_clobber, ncid)
    if (status /= nf90_noerr) call nchandle_error(status)

    if(lmoviez .eqv. .true.) then
      !create dimensions
      status = nf90_def_dim(ncid, "x", 1, xid)
      if (status /= nf90_noerr) call nchandle_error(status)
      status = nf90_def_dim(ncid, "y", jmax, yid)
      if (status /= nf90_noerr) call nchandle_error(status)
      status = nf90_def_dim(ncid, "z", ncklimit, zid)
      if (status /= nf90_noerr) call nchandle_error(status)
      status = nf90_def_dim(ncid, "t", nf90_unlimited, tid)
      if (status /= nf90_noerr) call nchandle_error(status)
    else
      !create dimensions
      status = nf90_def_dim(ncid, "x", imax, xid)
      if (status /= nf90_noerr) call nchandle_error(status)
      status = nf90_def_dim(ncid, "y", jmax, yid)
      if (status /= nf90_noerr) call nchandle_error(status)
      status = nf90_def_dim(ncid, "z", 1, zid)
      if (status /= nf90_noerr) call nchandle_error(status)
      status = nf90_def_dim(ncid, "t", nf90_unlimited, tid)
      if (status /= nf90_noerr) call nchandle_error(status)
    end if

    !create variables
    status = nf90_def_var(ncid, "u", nf90_float, (/xid, yid, zid, tid/), uid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "v", nf90_float, (/xid, yid, zid, tid/), vid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "w", nf90_float, (/xid, yid, zid, tid/), wid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "thl", nf90_float, (/xid, yid, zid, tid/), thlid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "qt", nf90_float, (/xid, yid, zid, tid/), qtid)
    if (status /= nf90_noerr) call nchandle_error(status)

    !turn off define mode
    status = nf90_enddef(ncid)
    if (status /= nf90_noerr) call nchandle_error(status)

  end subroutine initnetcdfmovie


  subroutine netcdfmovie

    use modglobal, only : rk3step,timee,dt_lim
    implicit none

    if (.not. lnetcdfmovie) return
    if (rk3step/=3) return
    if (timee==0) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+dtmovie
    dt_lim = min(dt_lim,tnext)

    call do_netcdfmovie
    nccall = nccall + 1

  end subroutine netcdfmovie


  subroutine do_netcdfmovie

    use typeSizes
    use netcdf
    use modfields
    use modmpi

    implicit none

    real, dimension(jmax, ncklimit)   :: uslicex, vslicex, wslicex, thlslicex, qtslicex
    real, dimension(imax, jmax)       :: uslicez, vslicez, wslicez, thlslicez, qtslicez
    integer                           :: status

    if(lmoviez .eqv. .true.) then
      ! Select a slice, in this case 10
      uslicex   = u0(slicex, 2:j1, 1:ncklimit)
      vslicex   = v0(slicex, 2:j1, 1:ncklimit)
      wslicex   = w0(slicex, 2:j1, 1:ncklimit)
      thlslicex = thl0(slicex, 2:j1, 1:ncklimit)
      qtslicex  = qt0(slicex, 2:j1, 1:ncklimit)

      status = nf90_put_var(ncid, uid, uslicex, (/1,1,1,nccall/), (/1, jmax, ncklimit , 1/))
      if(status /= nf90_noerr) call nchandle_error(status)
      status = nf90_put_var(ncid, vid, vslicex, (/1,1,1,nccall/), (/1, jmax, ncklimit , 1/))
      if(status /= nf90_noerr) call nchandle_error(status)
      status = nf90_put_var(ncid, wid, wslicex, (/1,1,1,nccall/), (/1, jmax, ncklimit , 1/))
      if(status /= nf90_noerr) call nchandle_error(status)
      status = nf90_put_var(ncid, thlid, thlslicex, (/1,1,1,nccall/), (/1, jmax, ncklimit , 1/))
      if(status /= nf90_noerr) call nchandle_error(status)
      status = nf90_put_var(ncid, qtid, qtslicex, (/1,1,1,nccall/), (/1, jmax, ncklimit , 1/))

    else
      ! Select a slice, in this case 10
      uslicez   = u0(2:i1, 2:j1, slicez)
      vslicez   = v0(2:i1, 2:j1, slicez)
      wslicez   = w0(2:i1, 2:j1, slicez)
      thlslicez = thl0(2:i1, 2:j1, slicez)
      qtslicez  = qt0(2:i1, 2:j1, slicez)

      status = nf90_put_var(ncid, uid, uslicez, (/1,1,1,nccall/), (/imax, jmax, 1 , 1/))
      if(status /= nf90_noerr) call nchandle_error(status)
      status = nf90_put_var(ncid, vid, vslicez, (/1,1,1,nccall/), (/imax, jmax, 1 , 1/))
      if(status /= nf90_noerr) call nchandle_error(status)
      status = nf90_put_var(ncid, wid, wslicez, (/1,1,1,nccall/), (/imax, jmax, 1 , 1/))
      if(status /= nf90_noerr) call nchandle_error(status)
      status = nf90_put_var(ncid, thlid, thlslicez, (/1,1,1,nccall/), (/imax, jmax, 1 , 1/))
      if(status /= nf90_noerr) call nchandle_error(status)
      status = nf90_put_var(ncid, qtid, qtslicez, (/1,1,1,nccall/), (/imax, jmax, 1 , 1/))

    end if

  end subroutine do_netcdfmovie

  subroutine exitnetcdfmovie

    use typeSizes
    use netcdf
    use modmpi, only : myid

    implicit none

    integer status

    if(.not.(lnetcdfmovie)) return

    write(6,*) "NETCDFMOVIE: Closing: ", myid

    status = nf90_close(ncid)
    if (status /= nf90_noerr) call nchandle_error(status)

  end subroutine exitnetcdfmovie

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

end module modnetcdfmovie
