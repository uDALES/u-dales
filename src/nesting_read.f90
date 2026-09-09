!!> \file nesting_read.f90
!!!  reads the one-way nesting input file nesting.inp.<expnr>.nc
!
!>
!!  Input only: this module knows the file format (nesting spec section 5) and
!!  nothing about the nesting scheme. It deliberately does not use nesting_scheme,
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
module nesting_read
  use mpi,    only : MPI_Wtime
  use netcdf
  use modmpi, only : myid
  implicit none
  save
  private
  public :: nestio_open, nestio_validate, nestio_read, nestio_read_block, nestio_close
  public :: nestio_hdr, nestio_header_type, nestio_tread
  public :: nestio_check_values, nestio_fill_value

  !> Schema versions this reader understands (global attribute
  !! udales_nesting_schema). Version 1 is the original file; version 2 adds the
  !! per-time-level post-correction flux_residual (with the fluid_lateral_area
  !! it was summed over) and the OPTIONAL full-domain initial-condition block
  !! u_init/v_init/w_init. A version 1 file must keep loading and running
  !! exactly as before, so everything version 2 adds is optional on read.
  integer, parameter :: NESTIO_SCHEMA_MIN = 1
  integer, parameter :: NESTIO_SCHEMA_MAX = 2

  !> Relative tolerance used when validating the header against the run.
  real, parameter :: nestio_tol = 1.e-10

  type nestio_header_type
    integer :: schema = 0, itot = 0, jtot = 0, ktot = 0, nzone = 0, ntime = 0
    real    :: xlen = 0., ylen = 0., rotation_deg = 0.
    logical :: divergence_corrected = .false.
    !> Schema 2: the file carries flux_residual(time), the residual of the data
    !! AS STORED, and fluid_lateral_area, the area it was summed over.
    logical :: has_flux_residual = .false.
    real    :: fluid_lateral_area = 0.
    !> Schema 2: the file carries u_init/v_init/w_init on the whole child grid.
    logical :: has_initial_condition = .false.
    real, allocatable :: time(:), xf(:), xh(:), yf(:), yh(:), zf(:), zh(:)
    real, allocatable :: rhobf(:), rhobh(:), net_volume_flux(:), flux_residual(:)
  end type nestio_header_type

  type(nestio_header_type) :: nestio_hdr

  !> Cumulative wall time spent in the read path of nestio_read /
  !! nestio_read_block [s]. This covers the WHOLE path -- the variable and
  !! dimension lookup as well as nf90_get_var -- because timing only
  !! nf90_get_var understated the cost: every call used to re-issue six
  !! metadata operations (inq_varid, inquire_variable and four
  !! inquire_dimension) outside the timer, on a shared-filesystem file, from
  !! every rank, for a varid and dimension lengths that cannot change while
  !! the file is open. V0c measured 46-49 % of child runtime in this path with
  !! only nf90_get_var counted (design section 10.5, V7).
  real :: nestio_tread = 0.

  integer            :: ncid    = -1
  logical            :: lopen   = .false.
  character(len=256) :: ncfname = ''

  ! --- variable metadata cache ------------------------------------------- #
  ! The file is opened read-only and never redefined, so a variable's id and
  ! its dimension lengths are fixed for the lifetime of the open file. There
  ! are 12 slab variables (4 faces x 3 components) plus at most 3 init
  ! variables, so a linear scan over a handful of names is far cheaper than
  ! one filesystem metadata round trip.
  integer, parameter :: NVCACHE = 32
  integer            :: nvcached = 0
  character(len=64)  :: vc_name(NVCACHE) = ''
  integer            :: vc_varid(NVCACHE) = -1
  integer            :: vc_dlen(4, NVCACHE) = 0

contains

  !> Open the nesting file read-only on every rank and populate nestio_hdr.
  !! Called from nesting::nesting_init. ierr /= 0 on failure; the netCDF
  !! error string, the offending variable and the file name are printed first.
  subroutine nestio_open(fname, ierr)
    character(len=*), intent(in)  :: fname
    integer,          intent(out) :: ierr

    integer :: idum, nz, nzh

    if (lopen) call nestio_close()

    ncfname = fname
    nestio_tread = 0.
    nvcached = 0
    vc_name = ''
    vc_varid = -1
    vc_dlen = 0

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

    ! --- schema 2 additions, all optional so that schema 1 still loads ---
    nestio_hdr%has_flux_residual     = .false.
    nestio_hdr%has_initial_condition = .false.
    nestio_hdr%fluid_lateral_area    = 0.
    if (allocated(nestio_hdr%flux_residual)) deallocate(nestio_hdr%flux_residual)

    if (nestio_hdr%schema >= 2) then
      call nestio_get_var1d('flux_residual', nestio_hdr%flux_residual, ierr)
      if (ierr /= nf90_noerr) return
      nestio_hdr%has_flux_residual = .true.

      call nestio_get_att_real('fluid_lateral_area', nestio_hdr%fluid_lateral_area, ierr)
      if (ierr /= nf90_noerr) return

      call nestio_get_att_int('has_initial_condition', idum, ierr)
      if (ierr /= nf90_noerr) return
      nestio_hdr%has_initial_condition = (idum /= 0)

      if (size(nestio_hdr%flux_residual) /= nestio_hdr%ntime) then
        if (myid == 0) then
          write(*,'(a,a)') ' nesting_read: flux_residual has the wrong length in ', trim(ncfname)
          write(*,'(a,i0,a,i0)') '   size = ', size(nestio_hdr%flux_residual), &
                                 ', ntime = ', nestio_hdr%ntime
        end if
        ierr = -1
        return
      end if
    end if

    ! --- internal consistency of the file itself ---
    call nestio_get_dim('nz',  nz,  ierr)
    if (ierr /= nf90_noerr) return
    call nestio_get_dim('nzh', nzh, ierr)
    if (ierr /= nf90_noerr) return

    if (nz /= nestio_hdr%nzone .or. nzh /= nestio_hdr%nzone + 1) then
      if (myid == 0) then
        write(*,'(a,a)') ' nesting_read: inconsistent zone dimensions in ', trim(ncfname)
        write(*,'(a,i0,a,i0,a,i0)') '   nzone = ', nestio_hdr%nzone, &
                                    ', nz = ', nz, ', nzh = ', nzh
      end if
      ierr = -1
      return
    end if

    if (size(nestio_hdr%net_volume_flux) /= nestio_hdr%ntime) then
      if (myid == 0) then
        write(*,'(a,a)') ' nesting_read: net_volume_flux has the wrong length in ', trim(ncfname)
        write(*,'(a,i0,a,i0)') '   size = ', size(nestio_hdr%net_volume_flux), &
                               ', ntime = ', nestio_hdr%ntime
      end if
      ierr = -1
      return
    end if

    ierr = nf90_noerr

  end subroutine nestio_open


  !> Compare the header against the run (modglobal) and the schema version.
  !! Called from nesting::nesting_init, right after nestio_open. Any
  !! mismatch is reported by name with both values and then aborts with stop 1.
  !! Only rank 0 prints.
  subroutine nestio_validate()
    use modglobal, only : itot, jtot, ktot, xlen, ylen, xf, xh, yf, yh, zf, zh

    integer :: nerr

    nerr = 0

    if (.not. lopen) then
      if (myid == 0) write(*,'(a)') ' nesting_read: nestio_validate called before nestio_open'
      stop 1
    end if

    if (nestio_hdr%schema < NESTIO_SCHEMA_MIN .or. nestio_hdr%schema > NESTIO_SCHEMA_MAX) then
      nerr = nerr + 1
      if (myid == 0) then
        write(*,'(a,a,a,i0,a,i0,a,i0)') ' nesting_read: mismatch in ', &
          'udales_nesting_schema', ': file = ', nestio_hdr%schema, &
          ', this reader supports ', NESTIO_SCHEMA_MIN, ' to ', NESTIO_SCHEMA_MAX
      end if
    end if
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

    ! Schema 2 initial-condition block: same stagger contract, plus the shape,
    ! since a full-domain array read at the wrong stagger would otherwise be one
    ! plane out in exactly the direction that matters.
    if (nestio_hdr%has_initial_condition) then
      call chk_init('u_init', 'xh yf zf', ktot,     jtot,     itot + 1)
      call chk_init('v_init', 'xf yh zf', ktot,     jtot + 1, itot)
      call chk_init('w_init', 'xf yf zh', ktot + 1, jtot,     itot)
    end if

    if (nerr > 0) then
      if (myid == 0) then
        write(*,'(a,i0,a,a,a)') ' nesting_read: ', nerr, &
          ' mismatch(es) between ', trim(ncfname), ' and the current run - aborting'
      end if
      stop 1
    end if

    if (myid == 0) then
      write(*,'(a,a,a)') ' nesting_read: ', trim(ncfname), ' validated against the run grid'
    end if

  contains

    !> Check one initial-condition variable: present, three dimensions in the
    !! Fortran order (z, y, x) with the lengths the stagger implies, and the
    !! matching stagger tag.
    subroutine chk_init(vname, expect, n1, n2, n3)
      character(len=*), intent(in) :: vname, expect
      integer,          intent(in) :: n1, n2, n3

      character(len=32) :: got
      integer :: varid, ndims, i, status, dlen(3)
      integer :: dimids(NF90_MAX_VAR_DIMS)
      integer :: want(3)

      want = (/ n1, n2, n3 /)

      status = nf90_inq_varid(ncid, vname, varid)
      if (status /= nf90_noerr) then
        nerr = nerr + 1
        if (myid == 0) write(*,'(a,a,a)') ' nesting_read: mismatch in ', vname, &
          ': has_initial_condition is set but the variable is absent'
        return
      end if

      status = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
      if (status /= nf90_noerr .or. ndims /= 3) then
        nerr = nerr + 1
        if (myid == 0) write(*,'(a,a,a,i0)') ' nesting_read: mismatch in ', vname, &
          ': expected a 3-dimensional variable, got ndims = ', ndims
        return
      end if

      do i = 1, 3
        status = nf90_inquire_dimension(ncid, dimids(i), len=dlen(i))
        if (status /= nf90_noerr) dlen(i) = -1
      end do

      do i = 1, 3
        if (dlen(i) /= want(i)) then
          nerr = nerr + 1
          if (myid == 0) write(*,'(a,a,a,i0,a,i0,a,i0)') ' nesting_read: mismatch in ', &
            vname, ': dimension ', i, ' is ', dlen(i), ', expected ', want(i)
          return
        end if
      end do

      got = ''
      status = nf90_get_att(ncid, varid, 'stagger', got)
      if (status /= nf90_noerr) then
        nerr = nerr + 1
        if (myid == 0) write(*,'(a,a,a)') ' nesting_read: MISMATCH ', vname, &
          ' has no stagger attribute'
        return
      end if
      if (trim(got) /= expect) then
        nerr = nerr + 1
        if (myid == 0) write(*,'(a,a,a,a,a,a)') ' nesting_read: MISMATCH ', vname, &
          ' stagger: file = "', trim(got), '", expected = "', expect//'"'
      end if

    end subroutine chk_init

    !> Check the stagger attribute of one velocity component on all four
    !! lateral slabs against the layout this reader assumes.
    subroutine chk_stagger(comp, expect)
      character(len=*), intent(in) :: comp, expect

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
          if (myid == 0) write(*,'(a,a,a)') ' nesting_read: MISMATCH ', trim(vname), &
            ' has no stagger attribute'
          cycle
        end if

        if (trim(got) /= expect) then
          nerr = nerr + 1
          if (myid == 0) write(*,'(a,a,a,a,a,a)') ' nesting_read: MISMATCH ', trim(vname), &
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
          write(*,'(a,a,a,i0,a,i0)') ' nesting_read: mismatch in ', name, &
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
          write(*,'(a,a,a,es22.14,a,es22.14)') ' nesting_read: mismatch in ', name, &
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
          write(*,'(a,a,a)') ' nesting_read: mismatch in ', name, ': absent from the file'
        end if
        return
      end if

      if (size(afile) /= size(arun)) then
        nerr = nerr + 1
        if (myid == 0) then
          write(*,'(a,a,a,i0,a,i0)') ' nesting_read: mismatch in size of ', name, &
            ': file = ', size(afile), ', run = ', size(arun)
        end if
        return
      end if

      do i = 1, size(arun)
        if (rneq(afile(i), arun(i), scale)) then
          nerr = nerr + 1
          if (myid == 0) then
            write(*,'(a,a,a,i0,a,es22.14,a,es22.14)') ' nesting_read: mismatch in ', name, &
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
  !! nesting when a new parent time level is needed. varname is e.g.
  !! 'u_west'; it, start2 and count2 are 1-based; buf is
  !! (n_zone_dim, n_z_dim, count2). ierr /= 0 on failure.
  !> Variable id and dimension lengths for ``varname``, from the cache.
  !! Populated on first use; a hit costs a string compare, a miss costs the
  !! six metadata calls this exists to stop repeating.
  subroutine nestio_varinfo(varname, varid, ndims, dlen, ierr)
    character(len=*), intent(in)  :: varname
    integer,          intent(out) :: varid, ndims, dlen(4)
    integer,          intent(out) :: ierr

    integer :: i, k
    integer :: dimids(NF90_MAX_VAR_DIMS)

    ierr = nf90_noerr
    dlen = 0

    do i = 1, nvcached
      if (trim(vc_name(i)) == trim(varname)) then
        varid = vc_varid(i)
        dlen  = vc_dlen(:, i)
        ndims = 4
        return
      end if
    end do

    ierr = nf90_inq_varid(ncid, trim(varname), varid)
    if (nestio_failed(ierr, 'nf90_inq_varid('//trim(varname)//')')) return

    ierr = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
    if (nestio_failed(ierr, 'nf90_inquire_variable('//trim(varname)//')')) return

    if (ndims /= 4) return          ! caller reports; nothing cached

    do k = 1, 4
      ierr = nf90_inquire_dimension(ncid, dimids(k), len=dlen(k))
      if (nestio_failed(ierr, 'nf90_inquire_dimension('//trim(varname)//')')) return
    end do

    if (nvcached < NVCACHE) then
      nvcached = nvcached + 1
      vc_name(nvcached)   = varname
      vc_varid(nvcached)  = varid
      vc_dlen(:, nvcached) = dlen
    end if
  end subroutine nestio_varinfo

  subroutine nestio_read(varname, it, start2, count2, buf, ierr)
    character(len=*), intent(in)  :: varname
    integer,          intent(in)  :: it, start2, count2
    real,             intent(out) :: buf(:,:,:)
    integer,          intent(out) :: ierr

    integer :: varid, ndims
    integer :: dlen(4), start(4), count(4)
    real    :: t0

    ierr = nf90_noerr

    if (.not. lopen) then
      write(*,'(a,i0,a,a)') ' nesting_read (rank ', myid, &
        '): nestio_read called before nestio_open, variable ', trim(varname)
      ierr = -1
      return
    end if

    if (count2 <= 0) return

    t0 = MPI_Wtime()
    call nestio_varinfo(varname, varid, ndims, dlen, ierr)
    if (ierr /= nf90_noerr) then
      nestio_tread = nestio_tread + (MPI_Wtime() - t0)
      return
    end if

    if (ndims /= 4) then
      nestio_tread = nestio_tread + (MPI_Wtime() - t0)
      call nestio_abortmsg(varname, 'expected a 4-dimensional variable')
      ierr = -1
      return
    end if

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

    ierr = nf90_get_var(ncid, varid, buf, start=start, count=count)
    nestio_tread = nestio_tread + (MPI_Wtime() - t0)
    if (nestio_failed(ierr, 'nf90_get_var('//trim(varname)//')')) return

    ierr = nestio_check_values(varname, it, buf)

  end subroutine nestio_read


  !> Read a rectangular block of one time-independent, full-domain variable
  !! (u_init/v_init/w_init, schema 2) for this rank. varname is e.g. 'u_init';
  !! start2/count2 index the y dimension and start3/count3 the x dimension,
  !! 1-based, and the whole vertical is read. buf is (nz, count2, count3),
  !! matching the file's Fortran dimension order (z, y, x). ierr /= 0 on failure.
  !! Called from nesting::nesting_init for nest_linitfromparent.
  subroutine nestio_read_block(varname, start2, count2, start3, count3, buf, ierr)
    character(len=*), intent(in)  :: varname
    integer,          intent(in)  :: start2, count2, start3, count3
    real,             intent(out) :: buf(:,:,:)
    integer,          intent(out) :: ierr

    integer :: varid, ndims, i
    integer :: dimids(NF90_MAX_VAR_DIMS)
    integer :: dlen(3), start(3), count(3)
    real    :: t0

    ierr = nf90_noerr

    if (.not. lopen) then
      write(*,'(a,i0,a,a)') ' nesting_read (rank ', myid, &
        '): nestio_read_block called before nestio_open, variable ', trim(varname)
      ierr = -1
      return
    end if

    if (count2 <= 0 .or. count3 <= 0) return

    ierr = nf90_inq_varid(ncid, trim(varname), varid)
    if (nestio_failed(ierr, 'nf90_inq_varid('//trim(varname)//')')) return

    ierr = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
    if (nestio_failed(ierr, 'nf90_inquire_variable('//trim(varname)//')')) return

    if (ndims /= 3) then
      call nestio_abortmsg(varname, 'expected a 3-dimensional variable')
      ierr = -1
      return
    end if

    do i = 1, 3
      ierr = nf90_inquire_dimension(ncid, dimids(i), len=dlen(i))
      if (nestio_failed(ierr, 'nf90_inquire_dimension('//trim(varname)//')')) return
    end do

    if (size(buf,1) /= dlen(1) .or. size(buf,2) /= count2 .or. size(buf,3) /= count3) then
      call nestio_abortmsg(varname, 'buffer shape does not match the file')
      ierr = -1
      return
    end if

    if (start2 < 1 .or. start2 + count2 - 1 > dlen(2) .or. &
        start3 < 1 .or. start3 + count3 - 1 > dlen(3)) then
      call nestio_abortmsg(varname, 'requested block is out of bounds')
      ierr = -1
      return
    end if

    start = (/ 1, start2, start3 /)
    count = (/ dlen(1), count2, count3 /)

    t0 = MPI_Wtime()
    ierr = nf90_get_var(ncid, varid, buf, start=start, count=count)
    nestio_tread = nestio_tread + (MPI_Wtime() - t0)
    if (nestio_failed(ierr, 'nf90_get_var('//trim(varname)//')')) return

    ierr = nestio_check_values(varname, 0, buf)

  end subroutine nestio_read_block


  !> Validate a buffer just read from variable varname at time level it (0
  !! for a time-independent block): every element must be finite and must not
  !! equal the variable's fill value. Missing or NaN data are an error, not a
  !! sentinel (spec section 5), and a NaN would otherwise pass straight
  !! through the flux assertion -- NaN makes Phi NaN and abs(NaN) > tol is
  !! false. Returns 0 when clean; otherwise reports (on the reading rank: the
  !! reads are per rank) the count and the first offending element and
  !! returns -1. Called by nestio_read and nestio_read_block on every read;
  !! public so the unit tests can poison a buffer and check it fires.
  integer function nestio_check_values(varname, it, buf) result(ierr)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

    character(len=*), intent(in) :: varname
    integer,          intent(in) :: it
    real,             intent(in) :: buf(:,:,:)

    integer :: i, j, k, nbad, ibad(3)
    real    :: fill, ftol
    logical :: lbad

    fill = nestio_fill_value(varname)
    ftol = 1.e-6*abs(fill)
    nbad = 0
    ibad = 0
    do k = 1, size(buf,3)
      do j = 1, size(buf,2)
        do i = 1, size(buf,1)
          lbad = .not. ieee_is_finite(buf(i,j,k))
          if (.not. lbad) lbad = abs(buf(i,j,k) - fill) <= ftol
          if (lbad) then
            if (nbad == 0) ibad = (/ i, j, k /)
            nbad = nbad + 1
          end if
        end do
      end do
    end do

    ierr = 0
    if (nbad == 0) return

    ierr = -1
    write(*,'(a,i0,a,i0,a,a,a,i0,a,a)') ' nesting_read (rank ', myid, '): ', nbad, &
      ' non-finite or fill-value element(s) in ', trim(varname), ' at time level ', it, &
      ' of file ', trim(ncfname)
    write(*,'(a,3i6,a,es12.4,a,es12.4)') '   first at buffer index ', ibad, &
      ', value ', buf(ibad(1),ibad(2),ibad(3)), ', fill value ', fill

  end function nestio_check_values


  !> The fill value the file declares for varname (_FillValue attribute), or
  !! netCDF's default fill for doubles when it declares none -- which is what
  !! an unwritten region of a variable reads back as.
  real function nestio_fill_value(varname)
    character(len=*), intent(in) :: varname

    integer :: varid, status
    real    :: fill

    nestio_fill_value = NF90_FILL_DOUBLE
    if (.not. lopen) return
    status = nf90_inq_varid(ncid, trim(varname), varid)
    if (status /= nf90_noerr) return
    status = nf90_get_att(ncid, varid, '_FillValue', fill)
    if (status == nf90_noerr) nestio_fill_value = fill

  end function nestio_fill_value


  !> Close the nesting file. Called from nesting::nesting_finalize.
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
      write(*,'(a,i0,a,a,a,a)') ' nesting_read (rank ', myid, '): ', trim(what), &
        ' failed on file ', trim(ncfname)
      write(*,'(a,a)') '   netCDF: ', trim(nf90_strerror(status))
    end if

  end function nestio_failed

  !> Report a non-netCDF failure against a named variable of the nesting file.
  subroutine nestio_abortmsg(varname, message)
    character(len=*), intent(in) :: varname, message

    write(*,'(a,i0,a,a,a,a,a,a)') ' nesting_read (rank ', myid, '): ', trim(message), &
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

end module nesting_read
