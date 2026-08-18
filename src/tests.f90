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
            tests_ibm_cell_lookup

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

end module tests
