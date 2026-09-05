!!> \file modnestingio.f90
!!!  reads the one-way nesting input file nesting.inp.<expnr>.nc
!
!>
!!  Input only: this module knows the file format (nesting spec section 5) and
!!  nothing about the nesting scheme. It deliberately does not use modnesting,
!!  so that it compiles and can be tested standalone.
!!
!!  The file is opened read-only on every rank with the serial netCDF library.
!!  Each rank then reads its own hyperslab of the decomposed (outermost spatial)
!!  index with nf90_get_var. Note that the CDL dimension order of the spec is
!!  reversed in Fortran: u_west(time, yf, zf, nzh) is seen here as
!!  u_west(nzh, zf, yf, time), so the decomposed index is dimension 3.
!!
!!  \author Maarten van Reeuwijk, Imperial College London
!
!  This file is part of uDALES.
!
! uDALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! uDALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1993-2009 Delft University of Technology, Wageningen University,
! Utrecht University, KNMI
!
module modnestingio
  use mpi,    only : MPI_Wtime
  use netcdf
  use modmpi, only : myid
  implicit none
  save
  private
  public :: nestio_open, nestio_validate, nestio_read, nestio_close
  public :: nestio_hdr, nestio_header_type, nestio_tread

  !> Schema version this reader understands (global attribute udales_nesting_schema).
  integer, parameter :: NESTIO_SCHEMA = 1

  !> Relative tolerance used when validating the header against the run.
  real, parameter :: nestio_tol = 1.e-10

  type nestio_header_type
    integer :: schema = 0, itot = 0, jtot = 0, ktot = 0, nzone = 0, ntime = 0
    real    :: xlen = 0., ylen = 0., rotation_deg = 0.
    logical :: divergence_corrected = .false.
    real, allocatable :: time(:), xf(:), xh(:), yf(:), yh(:), zf(:), zh(:)
    real, allocatable :: rhobf(:), rhobh(:), net_volume_flux(:)
  end type nestio_header_type

  type(nestio_header_type) :: nestio_hdr

  !> Cumulative wall time spent inside nf90_get_var in nestio_read [s].
  real :: nestio_tread = 0.

  integer            :: ncid    = -1
  logical            :: lopen   = .false.
  character(len=256) :: ncfname = ''

contains

  !> Open the nesting file read-only on every rank and populate nestio_hdr.
  !! Called from modnesting::nesting_init. ierr /= 0 on failure; the netCDF
  !! error string, the offending variable and the file name are printed first.
  subroutine nestio_open(fname, ierr)
    character(len=*), intent(in)  :: fname
    integer,          intent(out) :: ierr

    integer :: idum, nz, nzh

    if (lopen) call nestio_close()

    ncfname = fname
    nestio_tread = 0.

    ierr = nf90_open(trim(fname), NF90_NOWRITE, ncid)
    if (nestio_failed(ierr, 'nf90_open')) then
      ncid  = -1
      return
    end if
    lopen = .true.

    ! --- global attributes (all required, spec section 5) ---
    call nestio_get_att_int('udales_nesting_schema', nestio_hdr%schema, ierr)
    if (ierr /= nf90_noerr) return
    call nestio_get_att_int('itot',  nestio_hdr%itot,  ierr)
    if (ierr /= nf90_noerr) return
    call nestio_get_att_int('jtot',  nestio_hdr%jtot,  ierr)
    if (ierr /= nf90_noerr) return
    call nestio_get_att_int('ktot',  nestio_hdr%ktot,  ierr)
    if (ierr /= nf90_noerr) return
    call nestio_get_att_int('nzone', nestio_hdr%nzone, ierr)
    if (ierr /= nf90_noerr) return
    call nestio_get_att_int('divergence_corrected', idum, ierr)
    if (ierr /= nf90_noerr) return
    nestio_hdr%divergence_corrected = (idum /= 0)

    call nestio_get_att_real('xlen', nestio_hdr%xlen, ierr)
    if (ierr /= nf90_noerr) return
    call nestio_get_att_real('ylen', nestio_hdr%ylen, ierr)
    if (ierr /= nf90_noerr) return
    call nestio_get_att_real('rotation_deg', nestio_hdr%rotation_deg, ierr)
    if (ierr /= nf90_noerr) return

    ! --- coordinate and header variables ---
    call nestio_get_var1d('time', nestio_hdr%time, ierr)
    if (ierr /= nf90_noerr) return
    nestio_hdr%ntime = size(nestio_hdr%time)

    call nestio_get_var1d('xf', nestio_hdr%xf, ierr)
    if (ierr /= nf90_noerr) return
    call nestio_get_var1d('xh', nestio_hdr%xh, ierr)
    if (ierr /= nf90_noerr) return
    call nestio_get_var1d('yf', nestio_hdr%yf, ierr)
    if (ierr /= nf90_noerr) return
    call nestio_get_var1d('yh', nestio_hdr%yh, ierr)
    if (ierr /= nf90_noerr) return
    call nestio_get_var1d('zf', nestio_hdr%zf, ierr)
    if (ierr /= nf90_noerr) return
    call nestio_get_var1d('zh', nestio_hdr%zh, ierr)
    if (ierr /= nf90_noerr) return
    call nestio_get_var1d('rhobf', nestio_hdr%rhobf, ierr)
    if (ierr /= nf90_noerr) return
    call nestio_get_var1d('rhobh', nestio_hdr%rhobh, ierr)
    if (ierr /= nf90_noerr) return
    call nestio_get_var1d('net_volume_flux', nestio_hdr%net_volume_flux, ierr)
    if (ierr /= nf90_noerr) return

    ! --- internal consistency of the file itself ---
    call nestio_get_dim('nz',  nz,  ierr)
    if (ierr /= nf90_noerr) return
    call nestio_get_dim('nzh', nzh, ierr)
    if (ierr /= nf90_noerr) return

    if (nz /= nestio_hdr%nzone .or. nzh /= nestio_hdr%nzone + 1) then
      if (myid == 0) then
        write(*,'(a,a)') ' modnestingio: inconsistent zone dimensions in ', trim(ncfname)
        write(*,'(a,i0,a,i0,a,i0)') '   nzone = ', nestio_hdr%nzone, &
                                    ', nz = ', nz, ', nzh = ', nzh
      end if
      ierr = -1
      return
    end if

    if (size(nestio_hdr%net_volume_flux) /= nestio_hdr%ntime) then
      if (myid == 0) then
        write(*,'(a,a)') ' modnestingio: net_volume_flux has the wrong length in ', trim(ncfname)
        write(*,'(a,i0,a,i0)') '   size = ', size(nestio_hdr%net_volume_flux), &
                               ', ntime = ', nestio_hdr%ntime
      end if
      ierr = -1
      return
    end if

    ierr = nf90_noerr

  end subroutine nestio_open


  !> Compare the header against the run (modglobal) and the schema version.
  !! Called from modnesting::nesting_init, right after nestio_open. Any
  !! mismatch is reported by name with both values and then aborts with stop 1.
  !! Only rank 0 prints.
  subroutine nestio_validate()
    use modglobal, only : itot, jtot, ktot, xlen, ylen, xf, xh, yf, yh, zf, zh

    integer :: nerr

    nerr = 0

    if (.not. lopen) then
      if (myid == 0) write(*,'(a)') ' modnestingio: nestio_validate called before nestio_open'
      stop 1
    end if

    call chk_int('udales_nesting_schema', nestio_hdr%schema, NESTIO_SCHEMA)
    call chk_int('itot', nestio_hdr%itot, itot)
    call chk_int('jtot', nestio_hdr%jtot, jtot)
    call chk_int('ktot', nestio_hdr%ktot, ktot)

    call chk_real('xlen', nestio_hdr%xlen, xlen, xlen)
    call chk_real('ylen', nestio_hdr%ylen, ylen, ylen)

    call chk_arr('xf', nestio_hdr%xf, xf(1:itot),     xlen)
    call chk_arr('xh', nestio_hdr%xh, xh(1:itot + 1), xlen)
    call chk_arr('yf', nestio_hdr%yf, yf(1:jtot),     ylen)
    call chk_arr('yh', nestio_hdr%yh, yh(1:jtot + 1), ylen)
    call chk_arr('zf', nestio_hdr%zf, zf(1:ktot),     zh(ktot + 1))
    call chk_arr('zh', nestio_hdr%zh, zh(1:ktot + 1), zh(ktot + 1))

    ! Stagger tags. The grid can match while the data are laid out at the wrong
    ! staggered location, which would otherwise be read silently and wrongly.
    call chk_stagger('u', 'xh yf zf')
    call chk_stagger('v', 'xf yh zf')
    call chk_stagger('w', 'xf yf zh')

    if (nerr > 0) then
      if (myid == 0) then
        write(*,'(a,i0,a,a,a)') ' modnestingio: ', nerr, &
          ' mismatch(es) between ', trim(ncfname), ' and the current run - aborting'
      end if
      stop 1
    end if

    if (myid == 0) then
      write(*,'(a,a,a)') ' modnestingio: ', trim(ncfname), ' validated against the run grid'
    end if

  contains

    !> Check the stagger attribute of one velocity component on all four
    !! lateral slabs against the layout this reader assumes.
    subroutine chk_stagger(comp, expect)
      character(len=*), intent(in) :: comp, expect

      character(len=4), parameter :: faces(4) = (/ 'west', 'east', 'sout', 'nort' /)
      character(len=32) :: vname, got
      integer :: n, varid, status

      do n = 1, 4
        select case (n)
        case (1); vname = comp//'_west'
        case (2); vname = comp//'_east'
        case (3); vname = comp//'_south'
        case (4); vname = comp//'_north'
        end select

        status = nf90_inq_varid(ncid, trim(vname), varid)
        if (status /= nf90_noerr) cycle   ! absent slab: not this check's business

        got = ''
        status = nf90_get_att(ncid, varid, 'stagger', got)
        if (status /= nf90_noerr) then
          nerr = nerr + 1
          if (myid == 0) write(*,'(a,a,a)') ' modnestingio: MISMATCH ', trim(vname), &
            ' has no stagger attribute'
          cycle
        end if

        if (trim(got) /= expect) then
          nerr = nerr + 1
          if (myid == 0) write(*,'(a,a,a,a,a,a)') ' modnestingio: MISMATCH ', trim(vname), &
            ' stagger: file = "', trim(got), '", expected = "', expect//'"'
        end if
      end do

    end subroutine chk_stagger


    subroutine chk_int(name, ifile, irun)
      character(len=*), intent(in) :: name
      integer,          intent(in) :: ifile, irun

      if (ifile /= irun) then
        nerr = nerr + 1
        if (myid == 0) then
          write(*,'(a,a,a,i0,a,i0)') ' modnestingio: mismatch in ', name, &
            ': file = ', ifile, ', run = ', irun
        end if
      end if

    end subroutine chk_int

    subroutine chk_real(name, rfile, rrun, scale)
      character(len=*), intent(in) :: name
      real,             intent(in) :: rfile, rrun, scale

      if (rneq(rfile, rrun, scale)) then
        nerr = nerr + 1
        if (myid == 0) then
          write(*,'(a,a,a,es22.14,a,es22.14)') ' modnestingio: mismatch in ', name, &
            ': file = ', rfile, ', run = ', rrun
        end if
      end if

    end subroutine chk_real

    subroutine chk_arr(name, afile, arun, scale)
      character(len=*),  intent(in) :: name
      real, allocatable, intent(in) :: afile(:)
      real,              intent(in) :: arun(:), scale

      integer :: i

      if (.not. allocated(afile)) then
        nerr = nerr + 1
        if (myid == 0) then
          write(*,'(a,a,a)') ' modnestingio: mismatch in ', name, ': absent from the file'
        end if
        return
      end if

      if (size(afile) /= size(arun)) then
        nerr = nerr + 1
        if (myid == 0) then
          write(*,'(a,a,a,i0,a,i0)') ' modnestingio: mismatch in size of ', name, &
            ': file = ', size(afile), ', run = ', size(arun)
        end if
        return
      end if

      do i = 1, size(arun)
        if (rneq(afile(i), arun(i), scale)) then
          nerr = nerr + 1
          if (myid == 0) then
            write(*,'(a,a,a,i0,a,es22.14,a,es22.14)') ' modnestingio: mismatch in ', name, &
              '(', i, '): file = ', afile(i), ', run = ', arun(i)
          end if
          return
        end if
      end do

    end subroutine chk_arr

    !> Relative comparison; scale is a characteristic magnitude of the field so
    !! that coordinates at or near the origin are not compared against zero.
    logical function rneq(a, b, scale)
      real, intent(in) :: a, b, scale

      rneq = abs(a - b) > nestio_tol*max(abs(a), abs(b), abs(scale))

    end function rneq

  end subroutine nestio_validate


  !> Read one time level of one variable for this rank's range of the
  !! decomposed (outermost spatial, i.e. third Fortran) index. Called from
  !! modnesting when a new parent time level is needed. varname is e.g.
  !! 'u_west'; it, start2 and count2 are 1-based; buf is
  !! (n_zone_dim, n_z_dim, count2). ierr /= 0 on failure.
  subroutine nestio_read(varname, it, start2, count2, buf, ierr)
    character(len=*), intent(in)  :: varname
    integer,          intent(in)  :: it, start2, count2
    real,             intent(out) :: buf(:,:,:)
    integer,          intent(out) :: ierr

    integer :: varid, ndims, i
    integer :: dimids(NF90_MAX_VAR_DIMS)
    integer :: dlen(4), start(4), count(4)
    real    :: t0

    ierr = nf90_noerr

    if (.not. lopen) then
      write(*,'(a,i0,a,a)') ' modnestingio (rank ', myid, &
        '): nestio_read called before nestio_open, variable ', trim(varname)
      ierr = -1
      return
    end if

    if (count2 <= 0) return

    ierr = nf90_inq_varid(ncid, trim(varname), varid)
    if (nestio_failed(ierr, 'nf90_inq_varid('//trim(varname)//')')) return

    ierr = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
    if (nestio_failed(ierr, 'nf90_inquire_variable('//trim(varname)//')')) return

    if (ndims /= 4) then
      call nestio_abortmsg(varname, 'expected a 4-dimensional variable')
      ierr = -1
      return
    end if

    do i = 1, 4
      ierr = nf90_inquire_dimension(ncid, dimids(i), len=dlen(i))
      if (nestio_failed(ierr, 'nf90_inquire_dimension('//trim(varname)//')')) return
    end do

    ! dlen is in Fortran order: (zone, z, decomposed, time)
    if (size(buf,1) /= dlen(1) .or. size(buf,2) /= dlen(2) .or. size(buf,3) /= count2) then
      call nestio_abortmsg(varname, 'buffer shape does not match the file')
      if (myid == 0) then
        write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') '   buf = ', size(buf,1), ' x ', size(buf,2), &
          ' x ', size(buf,3), ', file = ', dlen(1), ' x ', dlen(2), ' x ', count2
      end if
      ierr = -1
      return
    end if

    if (it < 1 .or. it > dlen(4)) then
      call nestio_abortmsg(varname, 'time index out of range')
      ierr = -1
      return
    end if

    if (start2 < 1 .or. start2 + count2 - 1 > dlen(3)) then
      call nestio_abortmsg(varname, 'decomposed index range out of bounds')
      ierr = -1
      return
    end if

    start = (/ 1, 1, start2, it /)
    count = (/ dlen(1), dlen(2), count2, 1 /)

    t0 = MPI_Wtime()
    ierr = nf90_get_var(ncid, varid, buf, start=start, count=count)
    nestio_tread = nestio_tread + (MPI_Wtime() - t0)
    if (nestio_failed(ierr, 'nf90_get_var('//trim(varname)//')')) return

  end subroutine nestio_read


  !> Close the nesting file. Called from modnesting::nesting_finalize.
  !! The header contents are retained; a subsequent nestio_open replaces them.
  subroutine nestio_close()

    integer :: iret

    if (.not. lopen) return

    iret = nf90_close(ncid)
    if (nestio_failed(iret, 'nf90_close')) then
      ! reported above; carry on tearing down regardless
    end if

    ncid  = -1
    lopen = .false.

  end subroutine nestio_close


  ! ------------------------------------------------------------------ private

  !> .true. when status signals a netCDF error, which is then reported with the
  !! netCDF error string, the offending variable/call and the file name.
  logical function nestio_failed(status, what)
    integer,          intent(in) :: status
    character(len=*), intent(in) :: what

    nestio_failed = (status /= nf90_noerr)

    if (nestio_failed) then
      write(*,'(a,i0,a,a,a,a)') ' modnestingio (rank ', myid, '): ', trim(what), &
        ' failed on file ', trim(ncfname)
      write(*,'(a,a)') '   netCDF: ', trim(nf90_strerror(status))
    end if

  end function nestio_failed

  !> Report a non-netCDF failure against a named variable of the nesting file.
  subroutine nestio_abortmsg(varname, message)
    character(len=*), intent(in) :: varname, message

    write(*,'(a,i0,a,a,a,a,a,a)') ' modnestingio (rank ', myid, '): ', trim(message), &
      ' for variable ', trim(varname), ' in file ', trim(ncfname)

  end subroutine nestio_abortmsg

  !> Read a required integer global attribute.
  subroutine nestio_get_att_int(aname, ival, ierr)
    character(len=*), intent(in)  :: aname
    integer,          intent(out) :: ival
    integer,          intent(out) :: ierr

    ival = 0
    ierr = nf90_get_att(ncid, NF90_GLOBAL, aname, ival)
    if (nestio_failed(ierr, 'nf90_get_att('//aname//')')) return

  end subroutine nestio_get_att_int

  !> Read a required real global attribute.
  subroutine nestio_get_att_real(aname, rval, ierr)
    character(len=*), intent(in)  :: aname
    real,             intent(out) :: rval
    integer,          intent(out) :: ierr

    rval = 0.
    ierr = nf90_get_att(ncid, NF90_GLOBAL, aname, rval)
    if (nestio_failed(ierr, 'nf90_get_att('//aname//')')) return

  end subroutine nestio_get_att_real

  !> Length of a required dimension.
  subroutine nestio_get_dim(dname, n, ierr)
    character(len=*), intent(in)  :: dname
    integer,          intent(out) :: n
    integer,          intent(out) :: ierr

    integer :: dimid

    n = 0
    ierr = nf90_inq_dimid(ncid, dname, dimid)
    if (nestio_failed(ierr, 'nf90_inq_dimid('//dname//')')) return
    ierr = nf90_inquire_dimension(ncid, dimid, len=n)
    if (nestio_failed(ierr, 'nf90_inquire_dimension('//dname//')')) return

  end subroutine nestio_get_dim

  !> Read a required 1-D variable in full, allocating it to its own length.
  subroutine nestio_get_var1d(vname, arr, ierr)
    character(len=*),  intent(in)  :: vname
    real, allocatable, intent(out) :: arr(:)
    integer,           intent(out) :: ierr

    integer :: varid, ndims, n
    integer :: dimids(NF90_MAX_VAR_DIMS)

    ierr = nf90_inq_varid(ncid, vname, varid)
    if (nestio_failed(ierr, 'nf90_inq_varid('//vname//')')) return

    ierr = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
    if (nestio_failed(ierr, 'nf90_inquire_variable('//vname//')')) return

    if (ndims /= 1) then
      call nestio_abortmsg(vname, 'expected a 1-dimensional variable')
      ierr = -1
      return
    end if

    ierr = nf90_inquire_dimension(ncid, dimids(1), len=n)
    if (nestio_failed(ierr, 'nf90_inquire_dimension('//vname//')')) return

    allocate(arr(n))

    ierr = nf90_get_var(ncid, varid, arr)
    if (nestio_failed(ierr, 'nf90_get_var('//vname//')')) return

  end subroutine nestio_get_var1d

end module modnestingio
