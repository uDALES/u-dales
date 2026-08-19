!> \file tests.f90
!> Module for testing functionality
!
! This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.

module tests
  !> In-solver test routines, executed via special runmode values defined in
  !! modglobal (TEST_*).  Each public entry point is dispatched from
  !! execute_runmode_actions in program.f90 and exercises solver
  !! infrastructure that is only reachable after full MPI/2DECOMP
  !! initialization.
  use decomp_2d
  use modmpi, only : myid, comm3d, mpierr, my_real, avexy_ibm, avey_ibm, sumx_ibm, sumy_ibm

  implicit none
  save
  public :: tests_read_sparse_ijk, tests_2decomp_init_exit, tests_mpi_operators, &
            tests_ibm_cell_lookup, tests_nudge, tests_ibm_wallfun

contains

  !> Report the currently initialized 2DECOMP layout and exit cleanly.
  !! This runmode is dispatched after the normal startup path has already
  !! called initmpi and init2decomp, so it must not initialize MPI or the
  !! decomposition a second time.
  subroutine tests_2decomp_init_exit
    write(*,*) xstart
    write(*,*) ystart
    write(*,*) zstart
    write(*,*) xend
    write(*,*) yend
    write(*,*) zend
    write(*,*) xsize
    write(*,*) ysize
    write(*,*) zsize

  end subroutine tests_2decomp_init_exit

  !> Test read_sparse_ijk by comparing with actual IBM initialization
  !> This test calls initibm which populates all global arrays,
  !> then compares with the generic read_sparse_ijk routine
  !> Returns .true. if all tests pass, .false. otherwise
  logical function tests_read_sparse_ijk()
    use modglobal,    only : runmode
    use readinput, only : read_sparse_ijk
    use modibm,       only : initibm
    use modibm,       only : solid_info_u, solid_info_v, solid_info_w, solid_info_c
    use modibm,       only : bound_info_u, bound_info_v, bound_info_w, bound_info_c
    use modibm,       only : nsolpts_u, nsolpts_v, nsolpts_w, nsolpts_c
    use modibm,       only : nbndpts_u, nbndpts_v, nbndpts_w, nbndpts_c
    use initfac,      only : readfacetfiles

    implicit none
    
    integer :: npts_loc_new
    integer, allocatable :: ids_loc_new(:)
    integer, allocatable :: pts_loc_new(:,:)
    logical :: all_passed
    
    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_read_sparse_ijk: SPARSE INPUT FILE TEST'
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'Testing sparse solid_*.txt and fluid_boundary_*.txt files'
    end if

    ! Read facet files and initialize IBM - populates all global arrays
    call readfacetfiles
    call initibm
    
    all_passed = .true.
    
    ! Test solid_u
    call read_sparse_ijk('solid_u.txt', nsolpts_u, npts_loc_new, ids_loc_new, pts_loc_new, 1)
    if (.not. compare_solid(solid_info_u, npts_loc_new, ids_loc_new, pts_loc_new, 'solid_u')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test solid_v
    call read_sparse_ijk('solid_v.txt', nsolpts_v, npts_loc_new, ids_loc_new, pts_loc_new, 1)
    if (.not. compare_solid(solid_info_v, npts_loc_new, ids_loc_new, pts_loc_new, 'solid_v')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test solid_w
    call read_sparse_ijk('solid_w.txt', nsolpts_w, npts_loc_new, ids_loc_new, pts_loc_new, 1)
    if (.not. compare_solid(solid_info_w, npts_loc_new, ids_loc_new, pts_loc_new, 'solid_w')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test solid_c
    call read_sparse_ijk('solid_c.txt', nsolpts_c, npts_loc_new, ids_loc_new, pts_loc_new, 1)
    if (.not. compare_solid(solid_info_c, npts_loc_new, ids_loc_new, pts_loc_new, 'solid_c')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test fluid_boundary_u
    call read_sparse_ijk('fluid_boundary_u.txt', nbndpts_u, npts_loc_new, ids_loc_new, pts_loc_new, 1)
    if (.not. compare_boundary(bound_info_u, npts_loc_new, ids_loc_new, pts_loc_new, 'fluid_boundary_u')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test fluid_boundary_v
    call read_sparse_ijk('fluid_boundary_v.txt', nbndpts_v, npts_loc_new, ids_loc_new, pts_loc_new, 1)
    if (.not. compare_boundary(bound_info_v, npts_loc_new, ids_loc_new, pts_loc_new, 'fluid_boundary_v')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test fluid_boundary_w
    call read_sparse_ijk('fluid_boundary_w.txt', nbndpts_w, npts_loc_new, ids_loc_new, pts_loc_new, 1)
    if (.not. compare_boundary(bound_info_w, npts_loc_new, ids_loc_new, pts_loc_new, 'fluid_boundary_w')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test fluid_boundary_c
    call read_sparse_ijk('fluid_boundary_c.txt', nbndpts_c, npts_loc_new, ids_loc_new, pts_loc_new, 1)
    if (.not. compare_boundary(bound_info_c, npts_loc_new, ids_loc_new, pts_loc_new, 'fluid_boundary_c')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    if (all_passed .and. myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'ALL TESTS PASSED: tests_read_sparse_ijk'
      write(*, '(A)') '  Tested 8 files successfully'
      write(*, '(A)') '  All results match IBM initialization code'
      write(*, '(A)') '================================================'
    else if (.not. all_passed .and. myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'TESTS FAILED: tests_read_sparse_ijk'
      write(*, '(A)') '  One or more tests did not pass'
      write(*, '(A)') '================================================'
    end if
    
    tests_read_sparse_ijk = all_passed
    
  end function tests_read_sparse_ijk

  !> Compare sparse solid points returned by read_sparse_ijk against
  !! the reference arrays populated by initibm.
  function compare_solid(solid_info, npts_loc_new, ids_loc_new, pts_loc_new, label) result(passed)
    use modmpi, only : myid
    use modibm, only : solid_info_type

    type(solid_info_type), intent(in) :: solid_info
    integer, intent(in) :: npts_loc_new
    integer, intent(in) :: ids_loc_new(:), pts_loc_new(:,:)
    character(len=*), intent(in) :: label
    logical :: passed
    integer :: m
    
    passed = .true.
    
    if (npts_loc_new /= solid_info%nsolptsrank) then
      write(*, '(A,I0,A,A,A)') 'FAIL on rank ', myid, ': ', trim(label), ' count mismatch'
      passed = .false.
      return
    end if
    
    do m = 1, solid_info%nsolptsrank
      if (ids_loc_new(m) /= solid_info%solptsrank(m)) then
        write(*, '(A,I0,A,A,A)') 'FAIL on rank ', myid, ': ', trim(label), ' index mismatch'
        passed = .false.
        return
      end if
      if (pts_loc_new(m,1) /= solid_info%solpts_loc(m,1) .or. &
          pts_loc_new(m,2) /= solid_info%solpts_loc(m,2) .or. &
          pts_loc_new(m,3) /= solid_info%solpts_loc(m,3)) then
        write(*, '(A,I0,A,A,A)') 'FAIL on rank ', myid, ': ', trim(label), ' coordinate mismatch'
        passed = .false.
        return
      end if
    end do
    
  end function compare_solid
  
  !> Compare sparse fluid-boundary points returned by read_sparse_ijk
  !! against the reference arrays populated by initibm.
  function compare_boundary(bound_info, npts_loc_new, ids_loc_new, pts_loc_new, label) result(passed)
    use modmpi, only : myid
    use modibm, only : bound_info_type

    type(bound_info_type), intent(in) :: bound_info
    integer, intent(in) :: npts_loc_new
    integer, intent(in) :: ids_loc_new(:), pts_loc_new(:,:)
    character(len=*), intent(in) :: label
    logical :: passed
    integer :: m
    
    passed = .true.
    
    if (npts_loc_new /= bound_info%nbndptsrank) then
      write(*, '(A,I0,A,A,A)') 'FAIL on rank ', myid, ': ', trim(label), ' count mismatch'
      passed = .false.
      return
    end if
    
    do m = 1, bound_info%nbndptsrank
      if (ids_loc_new(m) /= bound_info%bndptsrank(m)) then
        write(*, '(A,I0,A,A,A)') 'FAIL on rank ', myid, ': ', trim(label), ' index mismatch'
        passed = .false.
        return
      end if
      if (pts_loc_new(m,1) /= bound_info%bndpts_loc(m,1) .or. &
          pts_loc_new(m,2) /= bound_info%bndpts_loc(m,2) .or. &
          pts_loc_new(m,3) /= bound_info%bndpts_loc(m,3)) then
        write(*, '(A,I0,A,A,A)') 'FAIL on rank ', myid, ': ', trim(label), ' coordinate mismatch'
        passed = .false.
        return
      end if
    end do
    
  end function compare_boundary

  !> Validate the IBM-aware MPI reduction operators (avexy_ibm, avey_ibm,
  !! sumx_ibm, sumy_ibm) against brute-force local reference sums.
  !! Requires a case with IBM geometry (e.g. case 100) so that the
  !! mask arrays are non-trivial.
  logical function tests_mpi_operators()
    use mpi
    use modglobal, only : ib, ie, jb, je, kb, ke, khc, runmode
    use modfields, only : initfields, IIc, IIu, IIv, IIw, IIuw, IIvw, IIuv, &
                          IIct, IIut, IIvt, IIwt, IIuwt, &
                          IIcs, IIus, IIvs, IIws, IIuws, IIvws, IIuvs
    use modibm, only : initibm, createmasks
    use initfac, only : readfacetfiles

    implicit none

    logical :: all_passed

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_mpi_operators: MODMPI IBM OPERATOR TEST'
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'Using case 100 masks to validate avexy_ibm/avey_ibm/sumx_ibm/sumy_ibm'
    end if

    call initfields
    call readfacetfiles
    call initibm
    call createmasks

    all_passed = .true.

    if (.not. check_loc_xy('C', 1, IIc, IIcs)) all_passed = .false.
    if (.not. check_loc_xy('U', 2, IIu, IIus)) all_passed = .false.
    if (.not. check_loc_xy('V', 3, IIv, IIvs)) all_passed = .false.
    if (.not. check_loc_xy('W', 4, IIw, IIws)) all_passed = .false.
    if (.not. check_loc_xy('WU', 5, IIuw, IIuws)) all_passed = .false.
    if (.not. check_loc_xy('VW', 6, IIvw, IIvws)) all_passed = .false.
    if (.not. check_loc_xy('UV', 7, IIuv, IIuvs)) all_passed = .false.

    if (.not. check_loc_y('C', 1, IIc, IIct)) all_passed = .false.
    if (.not. check_loc_y('U', 2, IIu, IIut)) all_passed = .false.
    if (.not. check_loc_y('V', 3, IIv, IIvt)) all_passed = .false.
    if (.not. check_loc_y('W', 4, IIw, IIwt)) all_passed = .false.
    if (.not. check_loc_y('WU', 5, IIuw, IIuwt)) all_passed = .false.

    if (all_passed .and. myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'ALL TESTS PASSED: tests_mpi_operators'
      write(*, '(A)') '================================================'
    else if ((.not. all_passed) .and. myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'TESTS FAILED: tests_mpi_operators'
      write(*, '(A)') '================================================'
    end if

    tests_mpi_operators = all_passed

  contains

    logical function check_loc_xy(label, loc_id, mask_3d, mask_1d)
      implicit none
      character(len=*), intent(in) :: label
      integer, intent(in) :: loc_id
      integer, intent(in) :: mask_3d(:,:,:)
      integer, intent(in) :: mask_1d(:)
      real, allocatable :: var_clean(:,:,:)
      real, allocatable :: got(:), exp(:), sum_local(:), sum_global(:)
      integer :: i, j, k

      allocate(var_clean(ib:ie,jb:je,kb:ke+khc))
      allocate(got(kb:ke+khc), exp(kb:ke+khc), sum_local(kb:ke+khc), sum_global(kb:ke+khc))

      do k = kb, ke + khc
        do j = jb, je
          do i = ib, ie
            var_clean(i,j,k) = 0.25 * real(loc_id) + 0.13 * real(i) - 0.07 * real(j) + 0.011 * real(k)
          end do
        end do
      end do

      got = 0.
      call avexy_ibm(got, var_clean, ib, ie, jb, je, kb, ke, khc, mask_3d(ib:ie,jb:je,kb:ke+khc), mask_1d, .true.)

      do k = kb, ke + khc
        sum_local(k) = sum(var_clean(ib:ie,jb:je,k) * real(mask_3d(ib:ie,jb:je,k)))
      end do
      call MPI_ALLREDUCE(sum_local, sum_global, size(sum_local), MY_REAL, MPI_SUM, comm3d, mpierr)

      do k = kb, ke + khc
        if (mask_1d(k) == 0) then
          exp(k) = -999.
        else
          exp(k) = sum_global(k) / real(mask_1d(k))
        end if
      end do

      check_loc_xy = compare_real_1d('avexy_ibm '//trim(label), got, exp)

      deallocate(var_clean, got, exp, sum_local, sum_global)
    end function check_loc_xy

    logical function check_loc_y(label, loc_id, mask_3d, mask_2d)
      implicit none
      character(len=*), intent(in) :: label
      integer, intent(in) :: loc_id
      integer, intent(in) :: mask_3d(:,:,:)
      integer, intent(in) :: mask_2d(:,:)
      real, allocatable :: var_clean(:,:,:)
      real, allocatable :: got_avg(:,:), got_sum_y(:,:), got_sum_x(:,:)
      real, allocatable :: exp_avg(:,:), exp_sum_y(:,:), exp_sum_x(:,:)
      real, allocatable :: sum_local(:,:), sum_global(:,:)
      real, allocatable :: sumx_local(:,:), sumx_global(:,:)
      integer :: i, j, k

      allocate(var_clean(ib:ie,jb:je,kb:ke))
      allocate(got_avg(ib:ie,kb:ke), got_sum_y(ib:ie,kb:ke), got_sum_x(jb:je,kb:ke))
      allocate(exp_avg(ib:ie,kb:ke), exp_sum_y(ib:ie,kb:ke), exp_sum_x(jb:je,kb:ke))
      allocate(sum_local(ib:ie,kb:ke), sum_global(ib:ie,kb:ke))
      allocate(sumx_local(jb:je,kb:ke), sumx_global(jb:je,kb:ke))

      do k = kb, ke
        do j = jb, je
          do i = ib, ie
            var_clean(i,j,k) = 0.5 * real(loc_id) + 0.21 * real(i) - 0.03 * real(j) + 0.017 * real(k)
          end do
        end do
      end do

      got_avg = 0.
      got_sum_y = 0.
      got_sum_x = 0.
      call avey_ibm(got_avg, var_clean, ib, ie, jb, je, kb, ke, mask_3d(ib:ie,jb:je,kb:ke), mask_2d)
      call sumy_ibm(got_sum_y, var_clean, ib, ie, jb, je, kb, ke, mask_3d(ib:ie,jb:je,kb:ke))
      call sumx_ibm(got_sum_x, var_clean, ib, ie, jb, je, kb, ke, mask_3d(ib:ie,jb:je,kb:ke))

      do k = kb, ke
        do i = ib, ie
          sum_local(i,k) = sum(var_clean(i,jb:je,k) * real(mask_3d(i,jb:je,k)))
        end do
      end do
      call MPI_ALLREDUCE(sum_local, sum_global, size(sum_local), MY_REAL, MPI_SUM, comm3d, mpierr)

      do k = kb, ke
        do i = ib, ie
          exp_sum_y(i,k) = sum_global(i,k)
          if (mask_2d(i,k) == 0) then
            exp_avg(i,k) = -999.
          else
            exp_avg(i,k) = sum_global(i,k) / real(mask_2d(i,k))
          end if
        end do
      end do

      do k = kb, ke
        do j = jb, je
          sumx_local(j,k) = sum(var_clean(ib:ie,j,k) * real(mask_3d(ib:ie,j,k)))
        end do
      end do
      call MPI_ALLREDUCE(sumx_local, sumx_global, size(sumx_local), MY_REAL, MPI_SUM, comm3d, mpierr)
      exp_sum_x = sumx_global

      check_loc_y = compare_real_2d('avey_ibm '//trim(label), got_avg, exp_avg)
      if (.not. compare_real_2d('sumy_ibm '//trim(label), got_sum_y, exp_sum_y)) check_loc_y = .false.
      if (.not. compare_real_2d_jk('sumx_ibm '//trim(label), got_sum_x, exp_sum_x)) check_loc_y = .false.

      deallocate(var_clean, got_avg, got_sum_y, got_sum_x, exp_avg, exp_sum_y, exp_sum_x, &
                 sum_local, sum_global, sumx_local, sumx_global)
    end function check_loc_y

    logical function compare_real_1d(label, got, exp)
      implicit none
      character(len=*), intent(in) :: label
      real, intent(in) :: got(:), exp(:)
      real :: max_abs
      integer :: imax(1)

      max_abs = maxval(abs(got - exp))
      compare_real_1d = max_abs <= 1.e-9
      if ((.not. compare_real_1d) .and. myid == 0) then
        imax = maxloc(abs(got - exp))
        write(*,'(A,1X,A,1X,ES12.4,1X,A,I0,1X,A,ES12.4,1X,A,ES12.4)') &
             'FAIL', trim(label), max_abs, 'idx', imax(1), 'got', got(imax(1)), 'exp', exp(imax(1))
      end if
    end function compare_real_1d

    logical function compare_real_2d(label, got, exp)
      implicit none
      character(len=*), intent(in) :: label
      real, intent(in) :: got(:,:), exp(:,:)
      real :: max_abs
      integer :: imax(2)

      max_abs = maxval(abs(got - exp))
      compare_real_2d = max_abs <= 1.e-9
      if ((.not. compare_real_2d) .and. myid == 0) then
        imax = maxloc(abs(got - exp))
        write(*,'(A,1X,A,1X,ES12.4,1X,A,I0,A,I0,1X,A,ES12.4,1X,A,ES12.4)') &
             'FAIL', trim(label), max_abs, 'idx', imax(1), ',', imax(2), 'got', got(imax(1),imax(2)), 'exp', exp(imax(1),imax(2))
      end if
    end function compare_real_2d

    logical function compare_real_2d_jk(label, got, exp)
      implicit none
      character(len=*), intent(in) :: label
      real, intent(in) :: got(:,:), exp(:,:)
      real :: max_abs
      integer :: imax(2)

      max_abs = maxval(abs(got - exp))
      compare_real_2d_jk = max_abs <= 1.e-9
      if ((.not. compare_real_2d_jk) .and. myid == 0) then
        imax = maxloc(abs(got - exp))
        write(*,'(A,1X,A,1X,ES12.4,1X,A,I0,A,I0,1X,A,ES12.4,1X,A,ES12.4)') &
             'FAIL', trim(label), max_abs, 'idx', imax(1), ',', imax(2), 'got', got(imax(1),imax(2)), 'exp', exp(imax(1),imax(2))
      end if
    end function compare_real_2d_jk

  end function tests_mpi_operators

  !> Validate the grid-cell lookup used by the IBM wall-function reconstruction.
  !!
  !! modibm::cell_index locates the cell containing a reconstruction point with
  !!     count(p >= grid)
  !! which is only equivalent to the backward search it replaced
  !!     findloc(p >= grid, .true., 1, back=.true.)
  !! when `grid` is monotonically increasing and fully initialised over its whole
  !! declared extent. (findloc is not used because the NVHPC runtime aborts on
  !! logical arrays with "FINDLOC: unimplemented for data type".)
  !!
  !! This test therefore checks three things:
  !!   1. the precondition - every grid array really is strictly increasing across
  !!      its full extent, including the halo elements, so no uninitialised entry
  !!      can break the leading-run property that count() relies on;
  !!   2. the equivalence - cell_index agrees with an explicit backward search on a
  !!      probe sweep covering below-first, on-node, one ULP either side of a node,
  !!      midpoint, and above-last positions;
  !!   3. the real path - every reconstruction point produced by initibm for the
  !!      loaded case resolves to the same indices as the backward search.
  !!
  !! Step 3 only means anything on a geometry that actually enters the
  !! reconstruction branch, so the test fails if it finds no such sections rather
  !! than reporting a vacuous pass. Run it against a case with non-grid-aligned or
  !! near-wall facets (tests/cases/100).
  !! Returns .true. if all checks pass, .false. otherwise.
  logical function tests_ibm_cell_lookup()
    use modglobal, only : runmode, xh, xf, yh, yf, zh, zf
    use modibm,    only : initibm, bound_info_u, bound_info_v, bound_info_w, bound_info_c
    use initfac,   only : readfacetfiles

    implicit none

    logical :: all_passed
    integer :: nrec, nrec_total

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_ibm_cell_lookup: IBM CELL LOOKUP TEST'
      write(*, '(A)') '------------------------------------------------'
    end if

    all_passed = .true.

    ! 1. Precondition: grid arrays strictly increasing over their full extent.
    if (.not. check_monotonic('xh', xh)) all_passed = .false.
    if (.not. check_monotonic('xf', xf)) all_passed = .false.
    if (.not. check_monotonic('yh', yh)) all_passed = .false.
    if (.not. check_monotonic('yf', yf)) all_passed = .false.
    if (.not. check_monotonic('zh', zh)) all_passed = .false.
    if (.not. check_monotonic('zf', zf)) all_passed = .false.

    ! 2. Equivalence with an explicit backward search over a probe sweep.
    if (.not. check_probe_sweep('xh', xh)) all_passed = .false.
    if (.not. check_probe_sweep('xf', xf)) all_passed = .false.
    if (.not. check_probe_sweep('yh', yh)) all_passed = .false.
    if (.not. check_probe_sweep('yf', yf)) all_passed = .false.
    if (.not. check_probe_sweep('zh', zh)) all_passed = .false.
    if (.not. check_probe_sweep('zf', zf)) all_passed = .false.

    ! 3. The real reconstruction points computed by initibm.
    call readfacetfiles
    call initibm

    nrec_total = 0
    if (.not. check_recids('bound_info_u', bound_info_u, nrec)) all_passed = .false.
    nrec_total = nrec_total + nrec
    if (.not. check_recids('bound_info_v', bound_info_v, nrec)) all_passed = .false.
    nrec_total = nrec_total + nrec
    if (.not. check_recids('bound_info_w', bound_info_w, nrec)) all_passed = .false.
    nrec_total = nrec_total + nrec
    if (.not. check_recids('bound_info_c', bound_info_c, nrec)) all_passed = .false.
    nrec_total = nrec_total + nrec

    ! A pass on a geometry that never reconstructs would be vacuous: the branch
    ! containing the lookup is simply never executed. Treat that as a failure of
    ! the fixture rather than a success of the code.
    if (nrec_total == 0) then
      if (myid == 0) then
        write(*, '(A)') 'FAIL: no facet sections required reconstruction'
        write(*, '(A)') '  This case never executes the cell lookup, so the test proves nothing.'
        write(*, '(A)') '  Use a geometry with non-grid-aligned or near-wall facets.'
      end if
      all_passed = .false.
    end if

    if (myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      if (all_passed) then
        write(*, '(A)')      'ALL TESTS PASSED: tests_ibm_cell_lookup'
        write(*, '(A,I0,A)') '  Checked 6 grid arrays and ', nrec_total, ' reconstruction sections'
      else
        write(*, '(A)') 'TESTS FAILED: tests_ibm_cell_lookup'
        write(*, '(A)') '  One or more checks did not pass'
      end if
      write(*, '(A)') '================================================'
    end if

    tests_ibm_cell_lookup = all_passed

  contains

    !> Reference implementation of the lookup: an explicit backward search for the
    !! last element of grid that p is greater than or equal to, 0 if there is none.
    !! Written out longhand so that it depends on no intrinsic beyond comparison,
    !! and so it stays valid even if the monotonicity precondition is violated.
    integer function ref_cell_index(p, grid)
      real, intent(in) :: p
      real, intent(in) :: grid(:)
      integer :: i

      ref_cell_index = 0
      do i = size(grid), 1, -1
        if (p >= grid(i)) then
          ref_cell_index = i
          return
        end if
      end do

    end function ref_cell_index

    !> Check that a grid array is strictly increasing across its whole extent.
    !! A single uninitialised halo element shows up here, and would otherwise
    !! silently break the leading-run property that count() depends on.
    logical function check_monotonic(label, grid)
      character(len=*), intent(in) :: label
      real, intent(in) :: grid(:)
      integer :: i

      check_monotonic = .true.
      do i = 2, size(grid)
        if (.not. (grid(i) > grid(i-1))) then
          if (myid == 0) then
            write(*,'(A,A,A,I0,A,I0,A,ES23.15,A,ES23.15)') &
                 'FAIL: ', trim(label), ' not strictly increasing between i=', i-1, &
                 ' and i=', i, ': ', grid(i-1), ' -> ', grid(i)
          end if
          check_monotonic = .false.
          return
        end if
      end do

    end function check_monotonic

    !> Compare cell_index against the backward search, and optionally against an
    !! independently known expected index.
    logical function check_lookup(label, p, grid, expected)
      use modibm, only : cell_index
      character(len=*), intent(in) :: label
      real, intent(in) :: p
      real, intent(in) :: grid(:)
      integer, intent(in), optional :: expected
      integer :: got, ref

      check_lookup = .true.
      got = cell_index(p, grid)
      ref = ref_cell_index(p, grid)

      if (got /= ref) then
        write(*,'(A,I0,A,A,A,ES23.15,A,I0,A,I0)') 'FAIL on rank ', myid, ': ', trim(label), &
             ' p=', p, ' cell_index=', got, ' backward search=', ref
        check_lookup = .false.
      end if

      if (present(expected)) then
        if (got /= expected) then
          write(*,'(A,I0,A,A,A,ES23.15,A,I0,A,I0)') 'FAIL on rank ', myid, ': ', trim(label), &
               ' p=', p, ' cell_index=', got, ' expected=', expected
          check_lookup = .false.
        end if
      end if

    end function check_lookup

    !> Sweep probe points across a grid array. Covers the cases the reconstruction
    !! points can actually land on: outside either end of the grid, exactly on a
    !! node (the floating-point equality edge), one ULP either side of a node, and
    !! midway between nodes.
    logical function check_probe_sweep(label, grid)
      character(len=*), intent(in) :: label
      real, intent(in) :: grid(:)
      integer :: i, n
      real :: span

      check_probe_sweep = .true.
      n = size(grid)
      span = max(grid(n) - grid(1), 1.)

      ! Below the first node the lookup must return 0, above the last it must
      ! return n. These are the two results the caller's bounds check relies on.
      if (.not. check_lookup(label//' below-first', grid(1) - span, grid, 0)) check_probe_sweep = .false.
      if (.not. check_lookup(label//' above-last',  grid(n) + span, grid, n)) check_probe_sweep = .false.

      do i = 1, n
        if (.not. check_lookup(label//' on-node',    grid(i),                grid, i))   check_probe_sweep = .false.
        if (.not. check_lookup(label//' below-node', nearest(grid(i), -1.),  grid, i-1)) check_probe_sweep = .false.
        if (.not. check_lookup(label//' above-node', nearest(grid(i),  1.),  grid, i))   check_probe_sweep = .false.
        if (i < n) then
          if (.not. check_lookup(label//' midpoint', 0.5*(grid(i) + grid(i+1)), grid, i)) check_probe_sweep = .false.
        end if
      end do

    end function check_probe_sweep

    !> Re-derive the reconstruction indices stored by initibm and compare.
    !!
    !! initibm deallocates the global per-section arrays once it has distributed
    !! them, so this works from the surviving per-rank _loc copies.
    !!
    !! lskipsec is assigned for every section, but lcomprec is only assigned when
    !! lskipsec is .false.. Sections with lskipsec=.false. and lcomprec=.false. are
    !! exactly those that went through the reconstruction branch, so those are the
    !! ones whose recpts entry is defined and whose recids came from cell_index.
    !! Note initibm additionally overrides lcomprec_loc to .true. for sections whose
    !! reconstruction cell falls outside the rank's halo range; that only shrinks the
    !! set checked here, so the caller's requirement that some sections were found
    !! is what guards against the set collapsing to nothing.
    !! nchecked reports how many were found, so the caller can detect a fixture
    !! that never exercises the path.
    logical function check_recids(label, bound_info, nchecked)
      use modibm, only : bound_info_type
      character(len=*), intent(in) :: label
      type(bound_info_type), intent(in) :: bound_info
      integer, intent(out) :: nchecked
      integer :: n

      check_recids = .true.
      nchecked = 0

      do n = 1, bound_info%nfctsecsrank
        if (bound_info%lskipsec_loc(n)) cycle
        if (bound_info%lcomprec_loc(n)) cycle
        nchecked = nchecked + 1

        if (.not. check_lookup(label//' recids_u(1)', bound_info%recpts_loc(n,1), xh, bound_info%recids_u_loc(n,1))) check_recids = .false.
        if (.not. check_lookup(label//' recids_u(2)', bound_info%recpts_loc(n,2), yf, bound_info%recids_u_loc(n,2))) check_recids = .false.
        if (.not. check_lookup(label//' recids_u(3)', bound_info%recpts_loc(n,3), zf, bound_info%recids_u_loc(n,3))) check_recids = .false.

        if (.not. check_lookup(label//' recids_v(1)', bound_info%recpts_loc(n,1), xf, bound_info%recids_v_loc(n,1))) check_recids = .false.
        if (.not. check_lookup(label//' recids_v(2)', bound_info%recpts_loc(n,2), yh, bound_info%recids_v_loc(n,2))) check_recids = .false.
        if (.not. check_lookup(label//' recids_v(3)', bound_info%recpts_loc(n,3), zf, bound_info%recids_v_loc(n,3))) check_recids = .false.

        if (.not. check_lookup(label//' recids_w(1)', bound_info%recpts_loc(n,1), xf, bound_info%recids_w_loc(n,1))) check_recids = .false.
        if (.not. check_lookup(label//' recids_w(2)', bound_info%recpts_loc(n,2), yf, bound_info%recids_w_loc(n,2))) check_recids = .false.
        if (.not. check_lookup(label//' recids_w(3)', bound_info%recpts_loc(n,3), zh, bound_info%recids_w_loc(n,3))) check_recids = .false.

        if (.not. check_lookup(label//' recids_c(1)', bound_info%recpts_loc(n,1), xf, bound_info%recids_c_loc(n,1))) check_recids = .false.
        if (.not. check_lookup(label//' recids_c(2)', bound_info%recpts_loc(n,2), yf, bound_info%recids_c_loc(n,2))) check_recids = .false.
        if (.not. check_lookup(label//' recids_c(3)', bound_info%recpts_loc(n,3), zf, bound_info%recids_c_loc(n,3))) check_recids = .false.

        if (.not. check_recids) return
      end do

    end function check_recids

  end function tests_ibm_cell_lookup

  !> Validate the profile-nudging tendency applied by modforces::nudge.
  !!
  !! nudge relaxes the upper part of the domain towards prescribed profiles at
  !! rate 1/tnudge, starting at level kb+nnudge. Two properties matter and
  !! neither is obvious from the call site:
  !!   - the relaxation starts at kb+nnudge, so everything below must be
  !!     untouched;
  !!   - the CPU branch uses whole-array assignment, so the tendency is applied
  !!     across each array's entire declared extent, halo cells included. The
  !!     GPU branch reproduces that with explicit loops, and this test pins the
  !!     behaviour the two branches have to agree on.
  !!
  !! This exercises the host branch. The device kernels are covered by
  !! tests_cuda.f90::test_nudge_profiles, which runs under
  !! UDALES_RUN_CUDA_SELFTEST on a Debug GPU build.
  !! Returns .true. if all checks pass, .false. otherwise.
  logical function tests_nudge()
    use modglobal, only : runmode, ih, ihc, jh, jhc, kb, ke, nsv, &
                          ltempeq, lmoist, lnudge, lnudgevel, nnudge, tnudge
    use modfields, only : initfields, up, vp, thlp, qtp, svp, &
                          u0av, v0av, thl0av, qt0av, sv0av, &
                          uprof, vprof, thlprof, qtprof, svprof
    use modforces, only : nudge

    implicit none

    real, parameter :: initial_p = 2., av_value = 3., prof_value = 1., t_relax = 4.

    logical :: all_passed, saved_lnudge, saved_lnudgevel, saved_ltempeq, saved_lmoist
    integer :: saved_nnudge, ktest, n
    real    :: saved_tnudge, expected, tol

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_nudge: PROFILE NUDGING TENDENCY TEST'
      write(*, '(A)') '------------------------------------------------'
    end if

#if defined(_GPU)
    ! nudge writes the device-resident tendencies in a GPU build, so the host
    ! arrays this test inspects would never change. Say so explicitly rather
    ! than reporting an unexplained numeric mismatch, and do not pass vacuously.
    if (myid == 0) then
      write(*, '(A)') 'FAIL: runmode 1007 exercises the host branch of nudge.'
      write(*, '(A)') '  Run it against a CPU build (build/cpu/<type>/u-dales).'
      write(*, '(A)') '  The device kernels are covered by tests_cuda.f90::test_nudge_profiles,'
      write(*, '(A)') '  run with UDALES_RUN_CUDA_SELFTEST=1 on a Debug GPU build.'
    end if
    tests_nudge = .false.
    return
#endif

    call initfields

    all_passed = .true.

    if (ke - kb < 2) then
      if (myid == 0) write(*, '(A)') 'FAIL: need at least three vertical levels to test'
      tests_nudge = .false.
      return
    end if

    saved_lnudge    = lnudge
    saved_lnudgevel = lnudgevel
    saved_ltempeq   = ltempeq
    saved_lmoist    = lmoist
    saved_nnudge    = nnudge
    saved_tnudge    = tnudge

    lnudge    = .true.
    lnudgevel = .true.
    ltempeq   = .true.
    lmoist    = .true.
    nnudge    = 1
    tnudge    = t_relax
    ktest     = kb + nnudge

    up = initial_p ; vp = initial_p ; thlp = initial_p ; qtp = initial_p ; svp = initial_p
    u0av = av_value ; v0av = av_value ; thl0av = av_value ; qt0av = av_value ; sv0av = av_value
    uprof = prof_value ; vprof = prof_value ; thlprof = prof_value
    qtprof = prof_value ; svprof = prof_value

    call nudge

    expected = initial_p - (av_value - prof_value) / t_relax
    tol      = 64. * epsilon(1.) * max(1., abs(expected))

    ! Nudged levels, over the full declared extent including halos.
    if (.not. check_block('up',   up  (:,:,ktest:ke), expected, tol)) all_passed = .false.
    if (.not. check_block('vp',   vp  (:,:,ktest:ke), expected, tol)) all_passed = .false.
    if (.not. check_block('thlp', thlp(:,:,ktest:ke), expected, tol)) all_passed = .false.
    if (.not. check_block('qtp',  qtp (:,:,ktest:ke), expected, tol)) all_passed = .false.
    do n = 1, nsv
      if (.not. check_block('svp', svp(:,:,ktest:ke,n), expected, tol)) all_passed = .false.
    end do

    ! Levels below kb+nnudge must be untouched.
    if (.not. check_block('up below level',   up  (:,:,kb:ktest-1), initial_p, tol)) all_passed = .false.
    if (.not. check_block('vp below level',   vp  (:,:,kb:ktest-1), initial_p, tol)) all_passed = .false.
    if (.not. check_block('thlp below level', thlp(:,:,kb:ktest-1), initial_p, tol)) all_passed = .false.
    if (.not. check_block('qtp below level',  qtp (:,:,kb:ktest-1), initial_p, tol)) all_passed = .false.
    do n = 1, nsv
      if (.not. check_block('svp below level', svp(:,:,kb:ktest-1,n), initial_p, tol)) all_passed = .false.
    end do

    ! Halo cells must be nudged too, not just the interior. All four lateral
    ! faces are checked for every field, because an explicit-loop port that
    ! covers only part of the halo - ib-ih..ie, say - still reproduces every
    ! interior value and would pass an interior-only comparison.
    if (.not. check_halos('up',   up  (:,:,ktest:ke),   ih,  jh,  expected, tol)) all_passed = .false.
    if (.not. check_halos('vp',   vp  (:,:,ktest:ke),   ih,  jh,  expected, tol)) all_passed = .false.
    if (.not. check_halos('thlp', thlp(:,:,ktest:ke),   ih,  jh,  expected, tol)) all_passed = .false.
    if (.not. check_halos('qtp',  qtp (:,:,ktest:ke),   ih,  jh,  expected, tol)) all_passed = .false.
    do n = 1, nsv
      if (.not. check_halos('svp', svp(:,:,ktest:ke,n), ihc, jhc, expected, tol)) all_passed = .false.
    end do

    ! lnudge off must be a no-op.
    lnudge = .false.
    up = initial_p
    call nudge
    if (.not. check_block('disabled', up, initial_p, tol)) all_passed = .false.

    lnudge    = saved_lnudge
    lnudgevel = saved_lnudgevel
    ltempeq   = saved_ltempeq
    lmoist    = saved_lmoist
    nnudge    = saved_nnudge
    tnudge    = saved_tnudge

    if (myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      if (all_passed) then
        write(*, '(A)')      'ALL TESTS PASSED: tests_nudge'
        write(*, '(A,I0,A)') '  Checked momentum, temperature, moisture and ', nsv, ' scalar(s)'
      else
        write(*, '(A)') 'TESTS FAILED: tests_nudge'
        write(*, '(A)') '  One or more checks did not pass'
      end if
      write(*, '(A)') '================================================'
    end if

    tests_nudge = all_passed

  contains

    !> Assert every element of a block equals want, and report the worst offender.
    logical function check_block(label, got, want, tol)
      character(len=*), intent(in) :: label
      real, intent(in) :: got(:,:,:)
      real, intent(in) :: want, tol
      real :: worst

      check_block = .true.
      if (size(got) == 0) return

      worst = maxval(abs(got - want))
      if (worst > tol) then
        write(*,'(A,I0,A,A,A,ES12.4,A,ES12.4,A,ES12.4)') 'FAIL on rank ', myid, ': nudge ', &
             trim(label), ' worst deviation ', worst, ' (expected ', want, ', tol ', tol
        check_block = .false.
      end if

    end function check_block

    !> Assert all four lateral halo faces of a block equal want.
    !!
    !! got must be the full declared i and j extent of the field; the dummy is
    !! assumed-shape, so the faces are its first and last hx (hy) planes.
    logical function check_halos(label, got, hx, hy, want, tol)
      character(len=*), intent(in) :: label
      real,    intent(in) :: got(:,:,:)
      integer, intent(in) :: hx, hy
      real,    intent(in) :: want, tol
      integer :: ni, nj

      check_halos = .true.
      ni = size(got, 1)
      nj = size(got, 2)

      ! The extent must hold two disjoint faces plus an interior, otherwise the
      ! slices below would overlap and a missed face could be covered by the
      ! opposite one.
      if (hx > 0 .and. ni > 2*hx) then
        if (.not. check_block(label//' halo-x low',  got(1:hx,       :, :), want, tol)) check_halos = .false.
        if (.not. check_block(label//' halo-x high', got(ni-hx+1:ni, :, :), want, tol)) check_halos = .false.
      end if
      if (hy > 0 .and. nj > 2*hy) then
        if (.not. check_block(label//' halo-y low',  got(:, 1:hy,       :), want, tol)) check_halos = .false.
        if (.not. check_block(label//' halo-y high', got(:, nj-hy+1:nj, :), want, tol)) check_halos = .false.
      end if

    end function check_halos

  end function tests_nudge


  !> Check the host IBM wall functions against expectations derived
  !! independently of how they are written.
  !!
  !! tests_cuda.f90 already compares the device port against these routines, but
  !! that comparison is symmetric: a mistake present in both sides passes it, and
  !! it cannot run at all without a GPU. What follows instead asserts properties
  !! the routines must have whatever their implementation:
  !!
  !!   - a correction that must vanish (no solid neighbours, or a constant field)
  !!   - a single face whose value is worked out by hand
  !!   - accumulation rather than assignment
  !!   - locality: only listed boundary points move
  !!   - conservation: the per-facet totals and the field tendencies account for
  !!     the same flux
  !!
  !! Returns .true. if all checks pass.
  logical function tests_ibm_wallfun()
#if defined(_GPU)
    ! The host wall functions this checks are not compiled into a GPU release
    ! build, and on a GPU debug build the solver does not use them, so there is
    ! nothing here to exercise. The device port is covered by
    ! tests_cuda.f90::test_ibm_wallfunmom and test_ibm_wallfunheat.
    implicit none

    if (myid == 0) then
      write(*, '(A)') 'FAIL: runmode 1008 checks the host IBM wall functions.'
      write(*, '(A)') '  Run it against a CPU build (build/cpu/<type>/u-dales).'
      write(*, '(A)') '  The device port is covered by the CUDA self-tests.'
    end if
    tests_ibm_wallfun = .false.
#else
    use modglobal, only : runmode, ib, ie, ih, jb, je, jh, kb, ke, kh, &
                          libm, iwallmom, ltempeq, lmoist, lwritefac, iwallmoist, &
                          dx, dy, dzh, dx2i, nfcts, xhat, &
                          totheatflux, totqflux
    use modfields, only : initfields, u0, v0, w0, thl0, qt0, qtp, pres0, up, thlp
    use initfac,   only : faclGR, facqsat, fachurel, facf, facT
    use modsubgrid,     only : initsubgrid
    use modsubgriddata, only : ekh
    use initfac,   only : readfacetfiles
    use modibm,    only : initibm, createmasks, mask_c, &
                          diffc_corr, wallfunmom, wallfunheat, &
                          fac_tau_raw, fac_pres_raw, &
                          bound_info_u, bound_info_c

    implicit none

    real, parameter :: ek_uniform = 0.25
    logical :: all_passed
    integer :: i, j, k
    real    :: expected, got, tol, vol, flux_sum, delta_sum

    real, allocatable :: mask_c_save(:,:,:)
    real, allocatable :: thlp_before(:,:,:), up_before(:,:,:)
    real, allocatable :: qtp_before(:,:,:)
    logical, allocatable :: faclGR_save(:)
    real, allocatable :: facqsat_save(:), fachurel_save(:), facf_save(:,:)
    logical :: saved_lmoist
    integer :: saved_iwallmoist

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_ibm_wallfun: HOST IBM WALL FUNCTION TESTS'
      write(*, '(A)') '------------------------------------------------'
    end if

    all_passed = .true.

    if (.not. libm) then
      if (myid == 0) write(*, '(A)') 'FAIL: this runmode needs a case with libm = .true.'
      tests_ibm_wallfun = .false.
      return
    end if

    ! execute_runmode_actions runs before the solver builds any of this, so the
    ! geometry has to be set up here, as the other IBM runmodes do.
    call initfields
    call initsubgrid      ! allocates ekm and ekh, which the corrections read
    call readfacetfiles
    call initibm
    call createmasks

    if (bound_info_c%nbndptsrank < 1) then
      if (myid == 0) write(*, '(A)') 'FAIL: no cell-centred boundary points on this rank'
      tests_ibm_wallfun = .false.
      return
    end if

    allocate(mask_c_save, source=mask_c)
    allocate(thlp_before(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(up_before(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))

    tol = 1.e-10

    ! ---------------------------------------------------------------- 1
    ! No solid neighbours anywhere: nothing to cancel, so no correction.
    ! This is what a mask test written the wrong way round fails.
    mask_c = 1.
    ekh = ek_uniform
    call fill_ramp(thl0)
    thlp = 0.
    call diffc_corr(thl0, thlp, ih, jh, kh)
    if (maxval(abs(thlp)) > tol) then
      call report('correction applied with no solid neighbours', maxval(abs(thlp)))
      all_passed = .false.
    end if

    ! ---------------------------------------------------------------- 2
    ! Every neighbour solid but a constant field: every difference is zero, so
    ! the correction is zero whatever the coefficients are.
    mask_c = 0.
    thl0 = 17.
    thlp = 0.
    call diffc_corr(thl0, thlp, ih, jh, kh)
    if (maxval(abs(thlp)) > tol) then
      call report('correction applied to a constant field', maxval(abs(thlp)))
      all_passed = .false.
    end if

    ! ---------------------------------------------------------------- 3
    ! One isolated solid neighbour in +x, uniform diffusivity, unit jump.
    ! By hand: rhs = -0.5*(ekh + ekh)*(1 - 0)*dx2i = -ek_uniform*dx2i.
    i = bound_info_c%bndpts_loc(1,1)
    j = bound_info_c%bndpts_loc(1,2)
    k = bound_info_c%bndpts_loc(1,3)
    mask_c = 1.
    mask_c(i+1,j,k) = 0.
    ekh = ek_uniform
    thl0 = 0.
    thl0(i+1,j,k) = 1.
    thlp = 0.
    call diffc_corr(thl0, thlp, ih, jh, kh)
    expected = -ek_uniform * dx2i
    got = thlp(i,j,k)
    if (abs(got - expected) > tol * max(1., abs(expected))) then
      if (myid == 0) write(*,'(A,ES14.6,A,ES14.6)') &
        'FAIL: single-face correction got ', got, ' expected ', expected
      all_passed = .false.
    end if

    ! ---------------------------------------------------------------- 4
    ! Accumulates into rhs rather than assigning to it.
    thlp = 0.
    call diffc_corr(thl0, thlp, ih, jh, kh)
    call diffc_corr(thl0, thlp, ih, jh, kh)
    if (abs(thlp(i,j,k) - 2.*expected) > tol * max(1., abs(expected))) then
      call report('correction assigns instead of accumulating', thlp(i,j,k))
      all_passed = .false.
    end if

    ! ---------------------------------------------------------------- 5
    ! Only listed boundary points move. Marks every listed cell, then checks
    ! that nothing outside the marks changed.
    mask_c = 0.
    call fill_ramp(thl0)
    thlp = 0.
    call diffc_corr(thl0, thlp, ih, jh, kh)
    if (.not. only_boundary_points_touched()) all_passed = .false.

    mask_c = mask_c_save

    ! ---------------------------------------------------------------- 6
    ! Momentum: a fluid at rest exerts no stress, and reports none.
    if (iwallmom > 1 .and. bound_info_u%nfctsecsrank > 0) then
      u0 = 0.; v0 = 0.; w0 = 0.
      up = 0.
      call wallfunmom(xhat, up, bound_info_u)
      if (maxval(abs(up)) > tol) then
        call report('stress produced by a fluid at rest', maxval(abs(up)))
        all_passed = .false.
      end if
      if (allocated(fac_tau_raw)) then
        if (maxval(abs(fac_tau_raw)) > tol) then
          call report('facet stress reported for a fluid at rest', maxval(abs(fac_tau_raw)))
          all_passed = .false.
        end if
      end if

      ! ------------------------------------------------------------- 7
      ! The wall stress opposes the flow, and the momentum it removes from the
      ! field equals the momentum it reports on the facets. The second is a
      ! cross-check between two outputs written under different indexing, so an
      ! error in either scatter breaks it.
      call fill_ramp(u0); call fill_ramp(v0); call fill_ramp(w0)
      u0 = u0 + 1.
      call fill_ramp(thl0)
      thl0 = thl0 + 290.
      up = 0.
      up_before = up
      call wallfunmom(xhat, up, bound_info_u)

      if (.not. stress_opposes_flow()) all_passed = .false.

      if (allocated(fac_tau_raw)) then
        delta_sum = 0.
        do k = kb, ke+kh
          vol = dx*dy*dzh(k)
          do j = jb, je
            do i = ib, ie
              delta_sum = delta_sum + (up(i,j,k) - up_before(i,j,k)) * vol
            end do
          end do
        end do
        flux_sum = sum(fac_tau_raw)
        if (abs(delta_sum + flux_sum) > 1.e-6 * max(1., abs(flux_sum))) then
          if (myid == 0) write(*,'(A,ES14.6,A,ES14.6)') &
            'FAIL: momentum removed ', -delta_sum, ' but facets report ', flux_sum
          all_passed = .false.
        end if
      end if
    end if

    ! ---------------------------------------------------------------- 8
    ! Heat: the same conservation cross-check, plus the pressure accumulator
    ! recomputed here directly from the section list.
    if ((ltempeq .or. lmoist .or. lwritefac) .and. bound_info_c%nfctsecsrank > 0) then
      call fill_ramp(u0); call fill_ramp(v0); call fill_ramp(w0)
      u0 = u0 + 1.
      call fill_ramp(thl0); thl0 = thl0 + 290.
      call fill_ramp(pres0)
      if (lmoist) call fill_ramp(qt0)
      thlp = 0.
      thlp_before = thlp
      totheatflux = 0.
      call wallfunheat

      if (ltempeq) then
        delta_sum = 0.
        do k = kb, ke+kh
          vol = dx*dy*dzh(k)
          do j = jb, je
            do i = ib, ie
              delta_sum = delta_sum + (thlp(i,j,k) - thlp_before(i,j,k)) * vol
            end do
          end do
        end do
        if (abs(delta_sum + totheatflux) > 1.e-6 * max(1., abs(totheatflux))) then
          if (myid == 0) write(*,'(A,ES14.6,A,ES14.6)') &
            'FAIL: heat removed ', -delta_sum, ' but totheatflux is ', totheatflux
          all_passed = .false.
        end if
      end if

      if (allocated(fac_pres_raw)) then
        if (.not. pressure_accumulator_matches()) all_passed = .false.
      end if
    end if

    ! ---------------------------------------------------------------- 9
    ! The moisture wall function. Reached only for green-roof facets with
    ! iwallmoist = 2, a combination no case in the suite has, so it is set up
    ! here. What is asserted are properties of the flux itself rather than its
    ! formula: it only ever adds moisture, it is exactly zero at equilibrium,
    ! and what leaves the field is what totqflux reports.
    if (allocated(faclGR) .and. allocated(facqsat) .and. bound_info_c%nfctsecsrank > 0) then
      allocate(faclGR_save, source=faclGR)
      allocate(facqsat_save, source=facqsat)
      allocate(fachurel_save, source=fachurel)
      allocate(facf_save, source=facf)
      allocate(qtp_before(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
      saved_lmoist     = lmoist
      saved_iwallmoist = iwallmoist

      lmoist     = .true.
      iwallmoist = 2
      faclGR     = .true.
      fachurel   = 1.0     ! saturated soil, so the bare-soil term is qtair-qwall
      facf(:,4)  = 200.
      facf(:,5)  = 50.
      facqsat    = 0.02

      ! Air and surface far enough apart that every section gets a non-zero
      ! heat transfer coefficient, which is what opens the moisture branch.
      call fill_ramp(u0); call fill_ramp(v0); call fill_ramp(w0)
      u0 = u0 + 1.
      call fill_ramp(thl0); thl0 = thl0 + 290.
      if (allocated(facT)) facT(:,1) = 300.

      ! (a) Equilibrium: air already at the surface humidity, saturated soil.
      !     Both terms of the flux vanish, so nothing may move.
      qt0 = 0.02
      qtp = 0.
      totqflux = 0.
      call wallfunheat
      if (maxval(abs(qtp)) > tol .or. abs(totqflux) > tol) then
        call report('moisture flux at equilibrium', max(maxval(abs(qtp)), abs(totqflux)))
        all_passed = .false.
      end if

      ! (b) Air wetter than the surface. The flux is a min(0., .) so it must
      !     clamp: a green roof never condenses moisture out of the air here.
      qt0 = 0.05
      qtp = 0.
      totqflux = 0.
      call wallfunheat
      if (maxval(abs(qtp)) > tol .or. abs(totqflux) > tol) then
        call report('moisture removed from air wetter than the surface', maxval(abs(qtp)))
        all_passed = .false.
      end if

      ! (c) Air drier than the surface: evaporation, so the tendency may only
      !     add moisture, and it must balance what totqflux reports.
      qt0 = 0.001
      qtp = 0.
      qtp_before = qtp
      totqflux = 0.
      call wallfunheat
      if (maxval(abs(qtp)) <= tol) then
        call report('no moisture flux from a drying surface', maxval(abs(qtp)))
        all_passed = .false.
      end if
      if (minval(qtp) < -tol) then
        call report('evaporation removed moisture somewhere', minval(qtp))
        all_passed = .false.
      end if
      delta_sum = 0.
      do k = kb, ke+kh
        vol = dx*dy*dzh(k)
        do j = jb, je
          do i = ib, ie
            delta_sum = delta_sum + (qtp(i,j,k) - qtp_before(i,j,k)) * vol
          end do
        end do
      end do
      if (abs(delta_sum + totqflux) > 1.e-6 * max(1., abs(totqflux))) then
        if (myid == 0) write(*,'(A,ES14.6,A,ES14.6)') &
          'FAIL: moisture added ', -delta_sum, ' but totqflux is ', totqflux
        all_passed = .false.
      end if

      ! (d) With no green-roof facets the branch must not run at all.
      faclGR = .false.
      qtp = 0.
      totqflux = 0.
      call wallfunheat
      if (maxval(abs(qtp)) > tol .or. abs(totqflux) > tol) then
        call report('moisture flux without a green roof', maxval(abs(qtp)))
        all_passed = .false.
      end if

      lmoist     = saved_lmoist
      iwallmoist = saved_iwallmoist
      faclGR     = faclGR_save
      facqsat    = facqsat_save
      fachurel   = fachurel_save
      facf       = facf_save
      deallocate(faclGR_save, facqsat_save, fachurel_save, facf_save, qtp_before)
    end if

    mask_c = mask_c_save
    deallocate(mask_c_save, thlp_before, up_before)

    if (myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      if (all_passed) then
        write(*, '(A)') 'ALL TESTS PASSED: tests_ibm_wallfun'
      else
        write(*, '(A)') 'TESTS FAILED: tests_ibm_wallfun'
      end if
      write(*, '(A)') '================================================'
    end if

    tests_ibm_wallfun = all_passed

    return

  contains

    subroutine report(what, value)
      character(len=*), intent(in) :: what
      real, intent(in) :: value
      if (myid == 0) write(*,'(A,A,A,ES12.4)') 'FAIL: ', what, ', worst ', value
    end subroutine report

    !> A field that varies in all three directions, so no stencil neighbour
    !! coincides with another.
    subroutine fill_ramp(a)
      real, intent(out) :: a(:,:,:)
      integer :: p, q, r
      do r = 1, size(a,3)
        do q = 1, size(a,2)
          do p = 1, size(a,1)
            a(p,q,r) = 0.125*p - 0.0625*q + 0.03125*r
          end do
        end do
      end do
    end subroutine fill_ramp

    logical function only_boundary_points_touched()
      integer :: m, p, q, r
      integer, allocatable :: listed(:,:,:)

      only_boundary_points_touched = .true.
      allocate(listed(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
      listed = 0
      do m = 1, bound_info_c%nbndptsrank
        listed(bound_info_c%bndpts_loc(m,1), &
               bound_info_c%bndpts_loc(m,2), &
               bound_info_c%bndpts_loc(m,3)) = 1
      end do

      do r = kb, ke+kh
        do q = jb-jh, je+jh
          do p = ib-ih, ie+ih
            if (listed(p,q,r) == 0 .and. abs(thlp(p,q,r)) > tol) then
              if (myid == 0) write(*,'(A,3I6)') &
                'FAIL: correction written outside the boundary point list at ', p, q, r
              only_boundary_points_touched = .false.
              deallocate(listed)
              return
            end if
          end do
        end do
      end do
      deallocate(listed)
    end function only_boundary_points_touched

    logical function stress_opposes_flow()
      integer :: p, q, r
      real    :: delta

      stress_opposes_flow = .true.
      do r = kb, ke+kh
        do q = jb, je
          do p = ib, ie
            delta = up(p,q,r) - up_before(p,q,r)
            if (abs(delta) <= tol) cycle
            if (delta * u0(p,q,r) > 0.) then
              if (myid == 0) write(*,'(A,3I6,A,ES12.4,A,ES12.4)') &
                'FAIL: wall stress adds momentum at ', p, q, r, ' delta ', delta, ' u ', u0(p,q,r)
              stress_opposes_flow = .false.
              return
            end if
          end do
        end do
      end do
    end function stress_opposes_flow

    !> fac_pres is a plain area-weighted sum of pres0 over the sections, with no
    !! wall-function physics in it, so it can be recomputed here outright.
    logical function pressure_accumulator_matches()
      use decomp_2d, only : zstart
      integer :: sec, p, q, r
      real, allocatable :: ref(:)
      real :: worst

      pressure_accumulator_matches = .true.
      allocate(ref(1:nfcts))
      ref = 0.
      do sec = 1, bound_info_c%nfctsecsrank
        p = bound_info_c%secbndpts_loc(sec,1) - zstart(1) + 1
        q = bound_info_c%secbndpts_loc(sec,2) - zstart(2) + 1
        r = bound_info_c%secbndpts_loc(sec,3) - zstart(3) + 1
        ref(bound_info_c%secfacids_loc(sec)) = ref(bound_info_c%secfacids_loc(sec)) &
                                             + pres0(p,q,r) * bound_info_c%secareas_loc(sec)
      end do

      worst = maxval(abs(ref - fac_pres_raw))
      if (worst > 1.e-8 * max(1., maxval(abs(ref)))) then
        if (myid == 0) write(*,'(A,ES12.4)') 'FAIL: facet pressure accumulator, worst ', worst
        pressure_accumulator_matches = .false.
      end if
      deallocate(ref)
    end function pressure_accumulator_matches

#endif

  end function tests_ibm_wallfun

end module tests
