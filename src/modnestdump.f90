!> \file modnestdump.f90
!!  Parent-side zone dump for one-way nesting (design section 6.2, plan item D1).
!!
!!  A nested child needs the parent's velocity only in a band just inside the
!!  child's lateral boundary -- the guard strip plus the relaxation ramp -- and,
!!  once, over the whole child box for its initial condition.  At the cadence
!!  the scheme requires (C_dump <= 2, i.e. a parent level every ~2 dx/U_max)
!!  full-domain field dumps are ~15x too large for a production parent, so this
!!  module writes ONLY that band, at its own cadence, in single precision.
!!
!!  Geometry.  The child box [x0, x0+xsize] x [y0, y0+ysize] is given in parent
!!  coordinates and must coincide with parent cell faces; the run aborts
!!  otherwise.  The band is the nestdump_nzone parent cells inside each lateral
!!  face of the box -- the child's guard + ramp expressed on the parent grid,
!!  rounded up, plus one cell so that the tangential slopes of a refined child's
!!  prolongation (udprep.nesting.conservative_interpolate) have a neighbour on
!!  the inner side.  Each rank writes, for each of the four strips (west, east,
!!  south, north) it intersects, the rectangular intersection of that strip with
!!  its own subdomain.  The strips overlap at the box corners, so a corner rank
!!  writes those cells twice; that is a few nzone^2 columns and is accepted for
!!  the simplicity of four rectangles.  A rank whose subdomain misses the band
!!  writes nothing at all.
!!
!!  Staggering.  For a strip covering global parent cells i1..i2, j1..j2 (1-based,
!!  as in the namelist itot/jtot), u is written at the i1..i2+1 x-faces (xh),
!!  v at the j1..j2+1 y-faces (yh) and w at the ktot+1 z-faces (zh, the lid
!!  value included), so the strip carries the complete staggered set of the
!!  cells it covers.  The upper faces come from the halo, which is exchanged
!!  right before this routine runs (program.f90: halos, then fielddump, then
!!  nestdump).  The full-box initial block follows the same rule.
!!
!!  Files.  nestdump.<ipx>.<ipy>.<expnr>.nc, one per rank that writes, with
!!  time as the unlimited dimension; nestdump_init.<ipx>.<ipy>.<expnr>.nc, one
!!  per rank whose subdomain meets the box, written once at the first dump
!!  time.  All index ranges are attributes, so any tool can assemble the box
!!  from the rank files without knowing the decomposition.  The format is
!!  specified in docs/udales-nesting-spec.md (section "nestdump files").
!!
!!  Cadence.  The next dump is due at tnextnestdump = btime + tnestdump, then
!!  every tnestdump, tested at rk3step == 3 exactly like modfielddump, so that
!!  lfielddump and lnestdump at the same interval write the same instants.
!!  When tnestdump is at or below the timestep every step is written.
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
! Copyright (C) 2016- the uDALES Team, Imperial College London.
!
module modnestdump
  use mpi
  use netcdf
  implicit none
  private
  public :: initnestdump, nestdump, exitnestdump
  public :: lnestdump, tnestdump, nestdump_x0, nestdump_y0, nestdump_xsize, &
            nestdump_ysize, nestdump_nzone, nestdump_linit
  save

  ! --- namelist &NESTDUMP ------------------------------------------------- !
  logical :: lnestdump      = .false. !< switch for the parent-side zone dump
  real    :: tnestdump      = 1.      !< dump interval [s]; <= dt means every step
  real    :: nestdump_x0    = 0.      !< child box origin, parent coordinates [m]
  real    :: nestdump_y0    = 0.
  real    :: nestdump_xsize = -1.     !< child box size [m]; must span whole parent cells
  real    :: nestdump_ysize = -1.
  integer :: nestdump_nzone = 0       !< band thickness inside each lateral face, parent cells
  logical :: nestdump_linit = .true.  !< write the whole box once, at the first dump

  ! --- module state -------------------------------------------------------- !
  integer, parameter :: NESTDUMP_SCHEMA = 1
  character(len=5), parameter :: facename(4) = (/ 'west ', 'east ', 'south', 'north' /)
  real(kind=4), parameter :: nestdump_fill = -999.

  !> One strip of the band, restricted to this rank.  Index ranges are global
  !! 1-based parent cells; u covers faces i1..i2+1, v faces j1..j2+1.
  type piece_type
    logical :: active = .false.
    integer :: i1 = 0, i2 = 0, j1 = 0, j2 = 0
    integer :: vid_u = -1, vid_v = -1, vid_w = -1
  end type piece_type
  type(piece_type) :: piece(4)

  integer :: ilo = 0, ihi = 0, jlo = 0, jhi = 0  !< box cells, global
  integer :: bi1 = 0, bi2 = 0, bj1 = 0, bj2 = 0  !< this rank's part of the box
  logical :: lband = .false.  !< this rank writes a band file
  logical :: lbox  = .false.  !< this rank meets the box (initial block)
  logical :: linitdone = .false.
  integer :: ncid = -1, nrec = 0, vid_time = -1
  real    :: tnextnestdump = 0.
  integer :: ndump = 0
  real(kind=8) :: bytes_run = 0.d0, twrite_run = 0.d0
  character(len=80) :: fname = 'nestdump.xxx.xxx.xxx.nc'
  character(len=80) :: fname_init = 'nestdump_init.xxx.xxx.xxx.nc'

contains

  !> Locate the box on the parent grid, work out this rank's pieces and open
  !! the band file.  Called after readinitfiles (timee/btime are known).
  subroutine initnestdump
    use modglobal, only : cexpnr, itot, jtot, xh, yh, dx, dy, btime, kb, ke
    use modmpi,    only : myid, cmyidx, cmyidy
    use decomp_2d, only : zstart, zend
    implicit none
    integer :: ig1, ig2, jg1, jg2, ni, nj, n, ierr

    if (.not. lnestdump) return

    if (nestdump_xsize <= 0. .or. nestdump_ysize <= 0.) call nestdump_abort( &
      'nestdump_xsize and nestdump_ysize must be > 0')
    if (nestdump_nzone < 1) call nestdump_abort('nestdump_nzone must be >= 1')
    if (tnestdump <= 0.) call nestdump_abort('tnestdump must be > 0')

    ! Box faces must be parent faces.  xh(i) is the west face of cell i.
    ilo = face_index(xh, itot, dx, nestdump_x0, 'x0')
    ihi = face_index(xh, itot, dx, nestdump_x0 + nestdump_xsize, 'x0 + xsize') - 1
    jlo = face_index(yh, jtot, dy, nestdump_y0, 'y0')
    jhi = face_index(yh, jtot, dy, nestdump_y0 + nestdump_ysize, 'y0 + ysize') - 1
    ni = ihi - ilo + 1
    nj = jhi - jlo + 1
    if (ni < 1 .or. nj < 1) call nestdump_abort('the child box is empty')
    if (2*nestdump_nzone > min(ni, nj)) call nestdump_abort( &
      'the band would cover the whole box: 2 * nestdump_nzone exceeds the box size in cells')

    ! This rank's subdomain, global cells (2DECOMP z-pencil).
    ig1 = zstart(1); ig2 = zend(1)
    jg1 = zstart(2); jg2 = zend(2)
    bi1 = max(ig1, ilo); bi2 = min(ig2, ihi)
    bj1 = max(jg1, jlo); bj2 = min(jg2, jhi)
    lbox = (bi1 <= bi2) .and. (bj1 <= bj2)

    ! The four strips, each intersected with this rank.
    call set_piece(piece(1), ilo, ilo + nestdump_nzone - 1, jlo, jhi, ig1, ig2, jg1, jg2)
    call set_piece(piece(2), ihi - nestdump_nzone + 1, ihi, jlo, jhi, ig1, ig2, jg1, jg2)
    call set_piece(piece(3), ilo, ihi, jlo, jlo + nestdump_nzone - 1, ig1, ig2, jg1, jg2)
    call set_piece(piece(4), ilo, ihi, jhi - nestdump_nzone + 1, jhi, ig1, ig2, jg1, jg2)
    lband = any(piece(:)%active)

    tnextnestdump = btime + tnestdump
    linitdone = .not. nestdump_linit
    ndump = 0

    fname(10:12) = cmyidx
    fname(14:16) = cmyidy
    fname(18:20) = cexpnr
    fname_init(15:17) = cmyidx
    fname_init(19:21) = cmyidy
    fname_init(23:25) = cexpnr

    if (myid == 0) then
      write(*, '(a)') 'nestdump: parent-side zone dump enabled'
      write(*, '(a,i0,a,i0,a,i0,a,i0,a)') '   child box: parent cells i = ', ilo, '..', ihi, &
        ', j = ', jlo, '..', jhi, ' (1-based)'
      write(*, '(a,i0,a,f0.4,a)') '   band thickness: ', nestdump_nzone, &
        ' parent cells; cadence ', tnestdump, ' s'
      write(*, '(a,l1)') '   initial block over the whole box: ', nestdump_linit
    end if

    if (lband) call open_band_file(ni, nj, kb, ke)

    n = 0
    if (lband) n = 1
    call MPI_ALLREDUCE(MPI_IN_PLACE, n, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (n == 0) call nestdump_abort('no rank intersects the band; check the box')
    if (myid == 0) write(*, '(a,i0,a)') '   ', n, ' rank(s) write a band file'
  end subroutine initnestdump

  !> Write the band (and, once, the initial block) when a dump is due.
  !! Placed right after fielddump in the time loop: same instant, same
  !! criterion, so both dumps of one run are the same field.
  subroutine nestdump
    use modglobal, only : rk3step, timee, kb, ke
    use modfields, only : u0, v0, w0
    use modmpi,    only : myid
    use decomp_2d, only : zstart
    implicit none
    integer :: n, ierr, iret, l1, l2, m1, m2, nk
    real(kind=8) :: t0, nbytes, nbytes_all, twall, twall_max

    if (.not. lnestdump) return
    if (rk3step /= 3) return
    if (timee < tnextnestdump) return
    tnextnestdump = tnextnestdump + tnestdump

    t0 = MPI_Wtime()
    nbytes = 0.d0
    nk = ke - kb + 1

    if (.not. linitdone) then
      if (lbox) call write_init_file(nbytes)
      linitdone = .true.
    end if

    if (lband) then
      nrec = nrec + 1
      iret = nf90_put_var(ncid, vid_time, real(timee, kind=4), start=(/ nrec /))
      call check(iret, 'put time')
      do n = 1, 4
        if (.not. piece(n)%active) cycle
        l1 = piece(n)%i1 - zstart(1) + 1
        l2 = piece(n)%i2 - zstart(1) + 1
        m1 = piece(n)%j1 - zstart(2) + 1
        m2 = piece(n)%j2 - zstart(2) + 1
        iret = nf90_put_var(ncid, piece(n)%vid_u, u0(l1:l2+1, m1:m2, kb:ke), &
                            start=(/ 1, 1, 1, nrec /), count=(/ l2-l1+2, m2-m1+1, nk, 1 /))
        call check(iret, 'put u_'//trim(facename(n)))
        iret = nf90_put_var(ncid, piece(n)%vid_v, v0(l1:l2, m1:m2+1, kb:ke), &
                            start=(/ 1, 1, 1, nrec /), count=(/ l2-l1+1, m2-m1+2, nk, 1 /))
        call check(iret, 'put v_'//trim(facename(n)))
        iret = nf90_put_var(ncid, piece(n)%vid_w, w0(l1:l2, m1:m2, kb:ke+1), &
                            start=(/ 1, 1, 1, nrec /), count=(/ l2-l1+1, m2-m1+1, nk+1, 1 /))
        call check(iret, 'put w_'//trim(facename(n)))
        nbytes = nbytes + block_bytes(l2 - l1 + 1, m2 - m1 + 1, nk)
      end do
      iret = nf90_sync(ncid)
      call check(iret, 'sync')
    end if

    twall = MPI_Wtime() - t0
    ndump = ndump + 1
    bytes_run = bytes_run + nbytes
    twrite_run = twrite_run + twall

    ! The cost report: first dump and every 100th, reduced over the ranks.
    if (ndump == 1 .or. mod(ndump, 100) == 0) then
      call MPI_REDUCE(nbytes, nbytes_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(twall, twall_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
      if (myid == 0) write(*, '(a,i0,a,f0.3,a,f0.3,a,f0.4,a)') 'nestdump: dump ', ndump, &
        ' at t = ', timee, ' s: ', nbytes_all/1.d6, ' MB written (all ranks) in ', &
        twall_max, ' s (slowest rank)'
    end if
  end subroutine nestdump

  !> Close the band file and report the run total.
  subroutine exitnestdump
    use modmpi, only : myid
    implicit none
    integer :: ierr, iret
    real(kind=8) :: bytes_all, twrite_max

    if (.not. lnestdump) return
    if (lband .and. ncid >= 0) then
      iret = nf90_close(ncid)
      call check(iret, 'close')
      ncid = -1
    end if
    call MPI_REDUCE(bytes_run, bytes_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(twrite_run, twrite_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    if (myid == 0) write(*, '(a,i0,a,f0.3,a,f0.3,a)') 'nestdump: ', ndump, ' dumps, ', &
      bytes_all/1.d6, ' MB in total (all ranks), ', twrite_max, &
      ' s in the write calls (slowest rank)'
  end subroutine exitnestdump

  ! ======================================================================= !
  ! helpers
  ! ======================================================================= !

  !> Bytes of one single-precision (u, v, w) block over ni x nj cells, nk levels.
  function block_bytes(ni, nj, nk) result(b)
    integer, intent(in) :: ni, nj, nk
    real(kind=8) :: b
    b = 4.d0*real((ni + 1)*nj*nk + ni*(nj + 1)*nk + ni*nj*(nk + 1), kind=8)
  end function block_bytes

  !> Index i with h(i) == x within a tolerance of the spacing; aborts otherwise.
  function face_index(h, ntot, spacing, x, label) result(i)
    real,             intent(in) :: h(:)
    integer,          intent(in) :: ntot
    real,             intent(in) :: spacing, x
    character(len=*), intent(in) :: label
    integer :: i, k
    real :: tol

    tol = 1.e-6*spacing
    i = -1
    do k = 1, ntot + 1
      if (abs(h(k) - x) <= tol) then
        i = k
        exit
      end if
    end do
    if (i < 0) call nestdump_abort('nestdump_'//trim(label)//' does not coincide with a parent '// &
      'cell face; the child box must be aligned with the parent grid')
  end function face_index

  !> Intersect the strip [s1,s2]x[t1,t2] with the rank [ig1,ig2]x[jg1,jg2].
  subroutine set_piece(p, s1, s2, t1, t2, ig1, ig2, jg1, jg2)
    type(piece_type), intent(inout) :: p
    integer,          intent(in)    :: s1, s2, t1, t2, ig1, ig2, jg1, jg2
    p%i1 = max(s1, ig1); p%i2 = min(s2, ig2)
    p%j1 = max(t1, jg1); p%j2 = min(t2, jg2)
    p%active = (p%i1 <= p%i2) .and. (p%j1 <= p%j2)
  end subroutine set_piece

  !> Global attributes shared by the band and the init files.
  subroutine put_global_atts(id, ni, nj)
    use modglobal, only : author, version, itot, jtot, ktot, dx, dy
    use modmpi,    only : myidx, myidy, nprocx, nprocy
    use decomp_2d, only : zstart, zend
    integer, intent(in) :: id, ni, nj
    character(len=12) :: cdate = '', ctime = ''
    integer :: iret

    call date_and_time(cdate, ctime)
    iret = nf90_put_att(id, NF90_GLOBAL, 'history', 'Created on '//trim(cdate)//' at '//trim(ctime))
    iret = nf90_put_att(id, NF90_GLOBAL, 'Source', trim(version))
    iret = nf90_put_att(id, NF90_GLOBAL, 'Author', trim(author))
    iret = nf90_put_att(id, NF90_GLOBAL, 'udales_nestdump_schema', NESTDUMP_SCHEMA)
    iret = nf90_put_att(id, NF90_GLOBAL, 'itot', itot)
    iret = nf90_put_att(id, NF90_GLOBAL, 'jtot', jtot)
    iret = nf90_put_att(id, NF90_GLOBAL, 'ktot', ktot)
    iret = nf90_put_att(id, NF90_GLOBAL, 'dx', dx)
    iret = nf90_put_att(id, NF90_GLOBAL, 'dy', dy)
    iret = nf90_put_att(id, NF90_GLOBAL, 'box_x0', nestdump_x0)
    iret = nf90_put_att(id, NF90_GLOBAL, 'box_y0', nestdump_y0)
    iret = nf90_put_att(id, NF90_GLOBAL, 'box_xsize', nestdump_xsize)
    iret = nf90_put_att(id, NF90_GLOBAL, 'box_ysize', nestdump_ysize)
    iret = nf90_put_att(id, NF90_GLOBAL, 'box_i_start', ilo)
    iret = nf90_put_att(id, NF90_GLOBAL, 'box_i_end', ihi)
    iret = nf90_put_att(id, NF90_GLOBAL, 'box_j_start', jlo)
    iret = nf90_put_att(id, NF90_GLOBAL, 'box_j_end', jhi)
    iret = nf90_put_att(id, NF90_GLOBAL, 'box_ni', ni)
    iret = nf90_put_att(id, NF90_GLOBAL, 'box_nj', nj)
    iret = nf90_put_att(id, NF90_GLOBAL, 'nzone', nestdump_nzone)
    iret = nf90_put_att(id, NF90_GLOBAL, 'tnestdump', tnestdump)
    iret = nf90_put_att(id, NF90_GLOBAL, 'myidx', myidx)
    iret = nf90_put_att(id, NF90_GLOBAL, 'myidy', myidy)
    iret = nf90_put_att(id, NF90_GLOBAL, 'nprocx', nprocx)
    iret = nf90_put_att(id, NF90_GLOBAL, 'nprocy', nprocy)
    iret = nf90_put_att(id, NF90_GLOBAL, 'rank_i_start', zstart(1))
    iret = nf90_put_att(id, NF90_GLOBAL, 'rank_i_end', zend(1))
    iret = nf90_put_att(id, NF90_GLOBAL, 'rank_j_start', zstart(2))
    iret = nf90_put_att(id, NF90_GLOBAL, 'rank_j_end', zend(2))
    iret = nf90_put_att(id, NF90_GLOBAL, 'index_convention', &
      'global 1-based parent cell indices; u(i) is the x-face xh(i) west of cell i, '// &
      'v(j) the y-face yh(j) south of cell j, w(k) the z-face zh(k) below cell k; '// &
      'each block covers cells i_start..i_end, j_start..j_end and the faces '// &
      'i_start..i_end+1, j_start..j_end+1, 1..ktot+1')
    call check(iret, 'global attributes')
  end subroutine put_global_atts

  !> Define one (u, v, w) block over global cells i1..i2, j1..j2, plus its
  !! four horizontal coordinate variables.  ``suffix`` is '_west' etc. for
  !! the strips and '' for the initial block; ``timeid`` < 0 means no time
  !! dimension.  Coordinate values are written by write_block_coords after
  !! enddef.
  subroutine define_block(id, suffix, i1, i2, j1, j2, ztid, zmid, timeid, vu, vv, vw)
    integer,          intent(in)  :: id, i1, i2, j1, j2, ztid, zmid, timeid
    character(len=*), intent(in)  :: suffix
    integer,          intent(out) :: vu, vv, vw
    integer :: xtid, xmid, ytid, ymid, iret

    iret = nf90_def_dim(id, 'xt'//suffix, i2 - i1 + 1, xtid)
    iret = nf90_def_dim(id, 'xm'//suffix, i2 - i1 + 2, xmid)
    iret = nf90_def_dim(id, 'yt'//suffix, j2 - j1 + 1, ytid)
    iret = nf90_def_dim(id, 'ym'//suffix, j2 - j1 + 2, ymid)
    call check(iret, 'def dims'//suffix)
    call def_coord(id, 'xt'//suffix, xtid, 'West-East displacement of cell centers')
    call def_coord(id, 'xm'//suffix, xmid, 'West-East displacement of cell edges')
    call def_coord(id, 'yt'//suffix, ytid, 'South-North displacement of cell centers')
    call def_coord(id, 'ym'//suffix, ymid, 'South-North displacement of cell edges')
    if (timeid >= 0) then
      iret = nf90_def_var(id, 'u'//suffix, NF90_FLOAT, (/ xmid, ytid, ztid, timeid /), vu)
      call check(iret, 'def u'//suffix)
      iret = nf90_def_var(id, 'v'//suffix, NF90_FLOAT, (/ xtid, ymid, ztid, timeid /), vv)
      call check(iret, 'def v'//suffix)
      iret = nf90_def_var(id, 'w'//suffix, NF90_FLOAT, (/ xtid, ytid, zmid, timeid /), vw)
      call check(iret, 'def w'//suffix)
    else
      iret = nf90_def_var(id, 'u'//suffix, NF90_FLOAT, (/ xmid, ytid, ztid /), vu)
      call check(iret, 'def u'//suffix)
      iret = nf90_def_var(id, 'v'//suffix, NF90_FLOAT, (/ xtid, ymid, ztid /), vv)
      call check(iret, 'def v'//suffix)
      iret = nf90_def_var(id, 'w'//suffix, NF90_FLOAT, (/ xtid, ytid, zmid /), vw)
      call check(iret, 'def w'//suffix)
    end if
    call put_var_atts(id, vu, 'West-East velocity', 'xh yf zf', i1, i2, j1, j2)
    call put_var_atts(id, vv, 'South-North velocity', 'xf yh zf', i1, i2, j1, j2)
    call put_var_atts(id, vw, 'Vertical velocity', 'xf yf zh', i1, i2, j1, j2)
  end subroutine define_block

  subroutine put_var_atts(id, vid, longname, stagger, i1, i2, j1, j2)
    integer,          intent(in) :: id, vid, i1, i2, j1, j2
    character(len=*), intent(in) :: longname, stagger
    integer :: iret
    iret = nf90_put_att(id, vid, 'longname', longname)
    iret = nf90_put_att(id, vid, 'units', 'm/s')
    iret = nf90_put_att(id, vid, 'stagger', stagger)
    iret = nf90_put_att(id, vid, 'i_start', i1)
    iret = nf90_put_att(id, vid, 'i_end', i2)
    iret = nf90_put_att(id, vid, 'j_start', j1)
    iret = nf90_put_att(id, vid, 'j_end', j2)
    iret = nf90_put_att(id, vid, '_FillValue', nestdump_fill)
    call check(iret, 'variable attributes')
  end subroutine put_var_atts

  !> A coordinate variable (single precision, like modstat_nc's).
  subroutine def_coord(id, name, dimid, longname)
    integer,          intent(in) :: id, dimid
    character(len=*), intent(in) :: name, longname
    integer :: vid, iret
    iret = nf90_def_var(id, name, NF90_FLOAT, (/ dimid /), vid)
    call check(iret, 'def '//name)
    iret = nf90_put_att(id, vid, 'longname', longname)
    iret = nf90_put_att(id, vid, 'units', 'm')
  end subroutine def_coord

  subroutine put_coord(id, name, values)
    integer,          intent(in) :: id
    character(len=*), intent(in) :: name
    real,             intent(in) :: values(:)
    integer :: vid, iret
    iret = nf90_inq_varid(id, name, vid)
    call check(iret, 'inq '//name)
    iret = nf90_put_var(id, vid, values)
    call check(iret, 'put '//name)
  end subroutine put_coord

  !> After enddef: the horizontal coordinates of one block.
  subroutine write_block_coords(id, suffix, i1, i2, j1, j2)
    use modglobal, only : xf, xh, yf, yh
    integer,          intent(in) :: id, i1, i2, j1, j2
    character(len=*), intent(in) :: suffix
    call put_coord(id, 'xt'//suffix, xf(i1:i2))
    call put_coord(id, 'xm'//suffix, xh(i1:i2+1))
    call put_coord(id, 'yt'//suffix, yf(j1:j2))
    call put_coord(id, 'ym'//suffix, yh(j1:j2+1))
  end subroutine write_block_coords

  !> Create (or reopen on a warm start) the band file of this rank.
  subroutine open_band_file(ni, nj, kb, ke)
    use modglobal, only : zf, zh, timee
    integer, intent(in) :: ni, nj, kb, ke
    integer :: iret, n, ztid, zmid, timeid, recid, ntimes
    real(kind=4), allocatable :: xtimes(:)
    logical :: exans

    inquire(file=trim(fname), exist=exans)
    if (exans) then
      ! Warm start of a phase that already wrote this file: continue after the
      ! last record before the restart time, as modstat_nc::open_nc does.
      iret = nf90_open(trim(fname), NF90_WRITE, ncid)
      call check(iret, 'reopen '//trim(fname))
      iret = nf90_inq_varid(ncid, 'time', vid_time)
      call check(iret, 'inq time')
      iret = nf90_inquire(ncid, unlimitedDimId=recid)
      iret = nf90_inquire_dimension(ncid, recid, len=ntimes)
      call check(iret, 'inquire time')
      nrec = 0
      if (ntimes > 0) then
        allocate(xtimes(ntimes))
        iret = nf90_get_var(ncid, vid_time, xtimes)
        call check(iret, 'get time')
        do while (nrec < ntimes)
          if (xtimes(nrec + 1) >= real(timee, kind=4) - spacing(1._4)) exit
          nrec = nrec + 1
        end do
        deallocate(xtimes)
      end if
      do n = 1, 4
        if (.not. piece(n)%active) cycle
        iret = nf90_inq_varid(ncid, 'u_'//trim(facename(n)), piece(n)%vid_u)
        call check(iret, 'inq u_'//trim(facename(n)))
        iret = nf90_inq_varid(ncid, 'v_'//trim(facename(n)), piece(n)%vid_v)
        call check(iret, 'inq v_'//trim(facename(n)))
        iret = nf90_inq_varid(ncid, 'w_'//trim(facename(n)), piece(n)%vid_w)
        call check(iret, 'inq w_'//trim(facename(n)))
      end do
      return
    end if

    iret = nf90_create(trim(fname), NF90_NETCDF4, ncid)
    call check(iret, 'create '//trim(fname))
    iret = nf90_put_att(ncid, NF90_GLOBAL, 'title', trim(fname))
    call put_global_atts(ncid, ni, nj)
    iret = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, timeid)
    iret = nf90_def_dim(ncid, 'zt', ke - kb + 1, ztid)
    iret = nf90_def_dim(ncid, 'zm', ke - kb + 2, zmid)
    call check(iret, 'def z dims')
    iret = nf90_def_var(ncid, 'time', NF90_FLOAT, (/ timeid /), vid_time)
    call check(iret, 'def time')
    iret = nf90_put_att(ncid, vid_time, 'longname', 'Time')
    iret = nf90_put_att(ncid, vid_time, 'units', 's')
    call def_coord(ncid, 'zt', ztid, 'Vertical displacement of cell centers')
    call def_coord(ncid, 'zm', zmid, 'Vertical displacement of cell edges')
    do n = 1, 4
      if (.not. piece(n)%active) cycle
      call define_block(ncid, '_'//trim(facename(n)), piece(n)%i1, piece(n)%i2, &
                        piece(n)%j1, piece(n)%j2, ztid, zmid, timeid, &
                        piece(n)%vid_u, piece(n)%vid_v, piece(n)%vid_w)
    end do
    iret = nf90_enddef(ncid)
    call check(iret, 'enddef')
    call put_coord(ncid, 'zt', zf(kb:ke))
    call put_coord(ncid, 'zm', zh(kb:ke+1))
    do n = 1, 4
      if (.not. piece(n)%active) cycle
      call write_block_coords(ncid, '_'//trim(facename(n)), piece(n)%i1, piece(n)%i2, &
                              piece(n)%j1, piece(n)%j2)
    end do
    nrec = 0
  end subroutine open_band_file

  !> The whole box on this rank, once.  Overwrites an existing file.
  subroutine write_init_file(nbytes)
    use modglobal, only : zf, zh, timee, kb, ke
    use modfields, only : u0, v0, w0
    use decomp_2d, only : zstart
    real(kind=8), intent(inout) :: nbytes
    integer :: id, iret, ztid, zmid, vu, vv, vw, vt, l1, l2, m1, m2, nk

    nk = ke - kb + 1
    iret = nf90_create(trim(fname_init), IOR(NF90_NETCDF4, NF90_CLOBBER), id)
    call check(iret, 'create '//trim(fname_init))
    iret = nf90_put_att(id, NF90_GLOBAL, 'title', trim(fname_init))
    call put_global_atts(id, ihi - ilo + 1, jhi - jlo + 1)
    iret = nf90_def_dim(id, 'zt', nk, ztid)
    iret = nf90_def_dim(id, 'zm', nk + 1, zmid)
    call check(iret, 'def z dims (init)')
    iret = nf90_def_var(id, 'time', NF90_FLOAT, vt)
    call check(iret, 'def time (init)')
    iret = nf90_put_att(id, vt, 'longname', 'Time')
    iret = nf90_put_att(id, vt, 'units', 's')
    call def_coord(id, 'zt', ztid, 'Vertical displacement of cell centers')
    call def_coord(id, 'zm', zmid, 'Vertical displacement of cell edges')
    call define_block(id, '', bi1, bi2, bj1, bj2, ztid, zmid, -1, vu, vv, vw)
    iret = nf90_enddef(id)
    call check(iret, 'enddef (init)')
    call put_coord(id, 'zt', zf(kb:ke))
    call put_coord(id, 'zm', zh(kb:ke+1))
    call write_block_coords(id, '', bi1, bi2, bj1, bj2)
    iret = nf90_put_var(id, vt, real(timee, kind=4))
    call check(iret, 'put time (init)')
    l1 = bi1 - zstart(1) + 1
    l2 = bi2 - zstart(1) + 1
    m1 = bj1 - zstart(2) + 1
    m2 = bj2 - zstart(2) + 1
    iret = nf90_put_var(id, vu, u0(l1:l2+1, m1:m2, kb:ke))
    call check(iret, 'put u (init)')
    iret = nf90_put_var(id, vv, v0(l1:l2, m1:m2+1, kb:ke))
    call check(iret, 'put v (init)')
    iret = nf90_put_var(id, vw, w0(l1:l2, m1:m2, kb:ke+1))
    call check(iret, 'put w (init)')
    iret = nf90_close(id)
    call check(iret, 'close (init)')
    nbytes = nbytes + block_bytes(l2 - l1 + 1, m2 - m1 + 1, nk)
  end subroutine write_init_file

  subroutine check(iret, what)
    integer,          intent(in) :: iret
    character(len=*), intent(in) :: what
    if (iret /= nf90_noerr) call nestdump_abort('netCDF error ('//trim(what)//'): '// &
      trim(nf90_strerror(iret)))
  end subroutine check

  subroutine nestdump_abort(msg)
    use modmpi, only : myid
    character(len=*), intent(in) :: msg
    integer :: ierr
    write(0, '(a,i0,a)') 'ERROR (nestdump, rank ', myid, '): '//trim(msg)
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end subroutine nestdump_abort

end module modnestdump
