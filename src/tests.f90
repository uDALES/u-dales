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
            tests_ibm_cell_lookup, tests_nudge, tests_ibm_wallfun, &
            tests_periodic_ebcorr, tests_masscorr, tests_ibmnorm

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
                          dx, dy, dzf, dzh, dxdydzfi, dxdydzhi, dx2i, nfcts, xhat, &
                          zhat, vec0, eps1, totheatflux, totqflux
    use modfields, only : initfields, u0, v0, w0, thl0, qt0, qtp, pres0, up, thlp
    use initfac,   only : faclGR, facqsat, fachurel, facf, facT
    use modsubgrid,     only : initsubgrid
    use modsubgriddata, only : ekh
    use initfac,   only : readfacetfiles
    use modibm,    only : initibm, createmasks, mask_c, &
                          diffc_corr, wallfunmom, wallfunheat, local_coords, &
                          check_wallfun_cache, &
                          fac_tau_raw, fac_pres_raw, &
                          bound_info_u, bound_info_c

    implicit none

    real, parameter :: ek_uniform = 0.25
    logical :: all_passed
    integer :: i, j, k
    real    :: expected, got, tol, vol, flux_sum, delta_sum
    real    :: span(3), strm(3), span2(3), strm2(3), uv(3)
    logical :: frame_valid, frame_valid2
    character(len=128) :: cache_problem

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
          ! dzf, matching wallfunmom. This read dzh until the cell-volume
          ! reciprocals went in; the two are equal on the uniform vertical grid
          ! every case in the suite uses, so the check could not tell them apart.
          vol = dx*dy*dzf(k)
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

    ! ---------------------------------------------------------------- 10
    ! The cell-volume reciprocals. Each is pinned to its own spacing, which no
    ! other check here can do: dzf and dzh are equal on the uniform vertical
    ! grid every case in the suite has, so a swap between dxdydzfi and dxdydzhi
    ! is invisible to the conservation tests above and to the CPU-GPU
    ! comparison. Asserted as a definition rather than by recomputing the same
    ! product, so the two sides cannot share a mistake.
    do k = kb, ke+kh
      if (abs(dxdydzfi(k)*(dx*dy*dzf(k)) - 1.) > 1.e-13) then
        call report('dxdydzfi is not 1/(dx*dy*dzf)', dxdydzfi(k)*(dx*dy*dzf(k)) - 1.)
        all_passed = .false.
        exit
      end if
    end do
    do k = kb, ke+kh
      if (abs(dxdydzhi(k)*(dx*dy*dzh(k)) - 1.) > 1.e-13) then
        call report('dxdydzhi is not 1/(dx*dy*dzh)', dxdydzhi(k)*(dx*dy*dzh(k)) - 1.)
        all_passed = .false.
        exit
      end if
    end do

    ! ---------------------------------------------------------------- 11
    ! local_coords normalises by a reciprocal rather than by three divisions.
    ! What has to hold is not the arithmetic but the frame: span and strm unit
    ! length and mutually orthogonal, both orthogonal to the normal, and both
    ! unchanged when the velocity is scaled. The last one is what a normalisation
    ! that silently stopped normalising would fail.
    uv = (/1.7, -0.9, 0.4/)
    call local_coords(uv, zhat, span, strm, frame_valid)
    if (.not. frame_valid) then
      call report('local_coords rejected a velocity across the wall', 0.)
      all_passed = .false.
    else
      if (abs(norm2(span) - 1.) > 1.e-13) then
        call report('local_coords span is not a unit vector', norm2(span) - 1.)
        all_passed = .false.
      end if
      if (abs(norm2(strm) - 1.) > 1.e-13) then
        call report('local_coords strm is not a unit vector', norm2(strm) - 1.)
        all_passed = .false.
      end if
      if (abs(dot_product(span, strm)) > 1.e-13 .or. &
          abs(dot_product(span, zhat)) > 1.e-13 .or. &
          abs(dot_product(strm, zhat)) > 1.e-13) then
        call report('local_coords frame is not orthogonal', &
                    max(abs(dot_product(span, strm)), abs(dot_product(span, zhat))))
        all_passed = .false.
      end if

      call local_coords(1000.*uv, zhat, span2, strm2, frame_valid2)
      if (.not. frame_valid2 .or. maxval(abs(span2 - span)) > 1.e-13 &
                             .or. maxval(abs(strm2 - strm)) > 1.e-13) then
        call report('local_coords frame changed when the velocity was scaled', &
                    max(maxval(abs(span2 - span)), maxval(abs(strm2 - strm))))
        all_passed = .false.
      end if

      ! A velocity straight into the wall has no tangential direction at all.
      call local_coords(zhat, zhat, span2, strm2, frame_valid2)
      if (frame_valid2) then
        call report('local_coords accepted a wall-normal velocity', 0.)
        all_passed = .false.
      end if
    end if

    ! ---------------------------------------------------------------- 12
    ! The per-section cache the wall functions read, against the expressions it
    ! replaced. The comparison lives in modibm, next to the arrays; this drives
    ! it. It is exact, so any difference is a plumbing error - a stencil on the
    ! wrong staggered grid, a cache column read under the wrong name, an index
    ! left global. On a GPU build the device mirrors are checked against the
    ! same cache by tests_cuda.f90; this half needs no GPU and so runs where the
    ! device tests cannot.
    call check_wallfun_cache(cache_problem)
    if (cache_problem /= '') then
      call report(trim(cache_problem), 0.)
      all_passed = .false.
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


  !> Check the periodic energy-balance correction against the physics it claims.
  !!
  !! periodicEBcorr removes, over the volume above the canopy and through the
  !! top level, the heat and moisture the facets put into the air over a
  !! periodic domain. The routine writes that out as three coupled scalars
  !! (R_theta_scaled, phi_theta_t and the M count behind them), so comparing
  !! against those expressions would only restate the code. What follows
  !! instead asserts the properties they exist to produce:
  !!
  !!   - the column-integrated tendency equals the domain-mean flux
  !!     tot_Tflux/(xlen*ylen), whatever sinkbase and fraction are. This is the
  !!     closure the scaling by ke/M exists to preserve, and it is what breaks
  !!     if M miscounts the levels;
  !!   - of that total, exactly (1-fraction) leaves through the top level and
  !!     the rest through the volume sink - the Grylls (2021) split;
  !!   - the flux entering is the sum over ranks, not the local contribution,
  !!     so the MPI_ALLREDUCE has to be there;
  !!   - support: nothing at or below sinkbase moves, and nothing outside
  !!     ib:ie, jb:je moves, halo cells included;
  !!   - the tendency is accumulated, is linear in the flux, and is switched off
  !!     by lperiodicEBcorr, ltempeq and lmoist.
  !!
  !! The column integral is only the domain-mean flux on a uniform vertical
  !! grid: the ke/M scaling weights levels by count, not by depth. That is a
  !! property of the scheme rather than of this port, so the test checks the
  !! grid is uniform and says so if it is not, instead of hiding the
  !! assumption.
  !!
  !! This exercises the host branch. The device kernels are covered by
  !! tests_cuda.f90::test_periodic_ebcorr, which runs under
  !! UDALES_RUN_CUDA_SELFTEST on a Debug GPU build.
  !! Returns .true. if all checks pass, .false. otherwise.
  logical function tests_periodic_ebcorr()
    use modglobal, only : runmode, ib, ie, jb, je, kb, ke, ih, jh, kh, &
                          ltempeq, lmoist, lperiodicEBcorr, sinkbase, fraction, &
                          totheatflux, totqflux, xlen, ylen, dzf, zh
    use modfields, only : initfields, thlp, qtp
    use modforces, only : periodicEBcorr
    use modmpi,    only : nprocs

    implicit none

    real, parameter :: initial_p = 2.

    logical :: all_passed, saved_lperiodicEBcorr, saved_ltempeq, saved_lmoist
    integer :: saved_sinkbase, k
    real    :: saved_fraction, saved_totheatflux, saved_totqflux
    real    :: unit_flux, mean_flux, tol

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_periodic_ebcorr: PERIODIC ENERGY BALANCE CORRECTION TEST'
      write(*, '(A)') '------------------------------------------------'
    end if

#if defined(_GPU)
    ! periodicEBcorr writes the device-resident tendencies in a GPU build, so
    ! the host arrays this test inspects would never change. Say so explicitly
    ! rather than reporting an unexplained numeric mismatch, and do not pass
    ! vacuously.
    if (myid == 0) then
      write(*, '(A)') 'FAIL: runmode 1009 exercises the host branch of periodicEBcorr.'
      write(*, '(A)') '  Run it against a CPU build (build/cpu/<type>/u-dales).'
      write(*, '(A)') '  The device kernels are covered by tests_cuda.f90::test_periodic_ebcorr,'
      write(*, '(A)') '  run with UDALES_RUN_CUDA_SELFTEST=1 on a Debug GPU build.'
    end if
    tests_periodic_ebcorr = .false.
    return
#endif

    call initfields

    all_passed = .true.

    ! Two levels have to sit above sinkbase for the top-level split to be
    ! separable from the volume sink, and one below it for the support check.
    if (ke - kb < 3) then
      if (myid == 0) write(*, '(A)') 'FAIL: need at least four vertical levels to test'
      tests_periodic_ebcorr = .false.
      return
    end if

    if (kb /= 1) then
      if (myid == 0) write(*, '(A)') 'FAIL: the ke/M scaling assumes kb = 1'
      tests_periodic_ebcorr = .false.
      return
    end if

    if (maxval(abs(dzf(kb:ke) - dzf(kb))) > 1.e-12 * dzf(kb) .or. &
        abs(zh(ke+1) - (ke - kb + 1) * dzf(kb)) > 1.e-12 * zh(ke+1)) then
      if (myid == 0) then
        write(*, '(A)') 'FAIL: this runmode needs a uniform vertical grid.'
        write(*, '(A)') '  The ke/M scaling weights levels by count, not by depth,'
        write(*, '(A)') '  so the column integral below is only the mean flux when dzf is constant.'
      end if
      tests_periodic_ebcorr = .false.
      return
    end if

    saved_lperiodicEBcorr = lperiodicEBcorr
    saved_ltempeq         = ltempeq
    saved_lmoist          = lmoist
    saved_sinkbase        = sinkbase
    saved_fraction        = fraction
    saved_totheatflux     = totheatflux
    saved_totqflux        = totqflux

    lperiodicEBcorr = .true.
    ltempeq         = .true.
    lmoist          = .true.

    ! Rank-dependent, so a routine that used the local flux instead of the
    ! reduced one gets a different answer on every rank but the first. The
    ! ranks contribute 1, 2, ... nprocs units, summing to nprocs(nprocs+1)/2.
    unit_flux = xlen * ylen
    mean_flux = 0.5 * real(nprocs) * real(nprocs + 1)

    ! 1. Column closure: whatever sinkbase and fraction are, the vertically
    !    integrated tendency is the domain-mean flux.
    if (.not. check_closure(2,        0.5, mean_flux)) all_passed = .false.
    if (.not. check_closure(ke/2,     0.5, mean_flux)) all_passed = .false.
    if (.not. check_closure(ke-2,     0.5, mean_flux)) all_passed = .false.
    if (.not. check_closure(2,        1.0, mean_flux)) all_passed = .false.
    if (.not. check_closure(2,        0.1, mean_flux)) all_passed = .false.

    ! 2. Linear in the flux: doubling what came in doubles what goes out. An
    !    additive constant anywhere in the scaling survives check 1 only if it
    !    is zero, but it would survive a single-amplitude check trivially.
    if (.not. check_closure_amplitude(2, 0.5, 2., 2. * mean_flux)) all_passed = .false.
    if (.not. check_closure_amplitude(2, 0.5, 0., 0.)) all_passed = .false.

    ! 3. The Grylls split: (1-fraction) of the total leaves through the top
    !    level, the rest through the volume sink. Read off as the excess the
    !    top level carries over the level below it, which the routine never
    !    forms.
    if (.not. check_split(2,    0.5)) all_passed = .false.
    if (.not. check_split(ke/2, 0.25)) all_passed = .false.
    if (.not. check_split(2,    1.0)) all_passed = .false.

    ! 4. Support and shape.
    sinkbase = ke/2
    fraction = 0.5
    call apply(unit_flux * real(myid + 1), unit_flux * real(myid + 1))

    tol = 1.e-12 * max(1., abs(mean_flux))

    ! Levels at and below sinkbase are untouched.
    do k = kb, sinkbase
      if (abs(thlp(ib,jb,k) - initial_p) > tol) then
        call report_scalar('thlp below sinkbase', thlp(ib,jb,k) - initial_p)
        all_passed = .false.
      end if
      if (abs(qtp(ib,jb,k) - initial_p) > tol) then
        call report_scalar('qtp below sinkbase', qtp(ib,jb,k) - initial_p)
        all_passed = .false.
      end if
    end do

    ! Every level above it is touched, so the loop bound is not off by one at
    ! the bottom either.
    do k = sinkbase+1, ke
      if (abs(thlp(ib,jb,k) - initial_p) <= tol) then
        call report_scalar('thlp level above sinkbase left untouched', real(k))
        all_passed = .false.
      end if
    end do

    ! Horizontally uniform over the interior, and nothing outside ib:ie, jb:je
    ! moves. The host branch loops over the interior only, so an explicit-loop
    ! port that spanned the halos - the opposite of the mistake nudge invites -
    ! would diverge here.
    do k = kb, ke
      if (maxval(abs(thlp(ib:ie,jb:je,k) - thlp(ib,jb,k))) > tol) then
        call report_scalar('thlp not horizontally uniform', &
                           maxval(abs(thlp(ib:ie,jb:je,k) - thlp(ib,jb,k))))
        all_passed = .false.
      end if
      if (maxval(abs(qtp(ib:ie,jb:je,k) - qtp(ib,jb,k))) > tol) then
        call report_scalar('qtp not horizontally uniform', &
                           maxval(abs(qtp(ib:ie,jb:je,k) - qtp(ib,jb,k))))
        all_passed = .false.
      end if
    end do

    if (.not. check_untouched('thlp halo', thlp)) all_passed = .false.
    if (.not. check_untouched('qtp halo',  qtp )) all_passed = .false.

    ! ke+kh sits above the loop bound and must stay put.
    if (kh > 0) then
      if (maxval(abs(thlp(ib:ie,jb:je,ke+1:ke+kh) - initial_p)) > tol) then
        call report_scalar('thlp above ke', &
                           maxval(abs(thlp(ib:ie,jb:je,ke+1:ke+kh) - initial_p)))
        all_passed = .false.
      end if
    end if

    ! 5. Heat and moisture follow the same algebra, so equal fluxes must give
    !    equal tendencies. Catches a copy-paste slip between the two branches.
    if (maxval(abs(thlp(ib:ie,jb:je,kb:ke) - qtp(ib:ie,jb:je,kb:ke))) > tol) then
      call report_scalar('heat and moisture tendencies differ', &
                         maxval(abs(thlp(ib:ie,jb:je,kb:ke) - qtp(ib:ie,jb:je,kb:ke))))
      all_passed = .false.
    end if

    ! 6. Switches. Each must silence its own field and leave the other alone.
    ltempeq = .false.
    lmoist  = .true.
    call apply(unit_flux, unit_flux)
    if (maxval(abs(thlp - initial_p)) > tol) then
      call report_scalar('ltempeq off still moved thlp', maxval(abs(thlp - initial_p)))
      all_passed = .false.
    end if
    if (maxval(abs(qtp - initial_p)) <= tol) then
      call report_scalar('ltempeq off also silenced qtp', 0.)
      all_passed = .false.
    end if

    ltempeq = .true.
    lmoist  = .false.
    call apply(unit_flux, unit_flux)
    if (maxval(abs(qtp - initial_p)) > tol) then
      call report_scalar('lmoist off still moved qtp', maxval(abs(qtp - initial_p)))
      all_passed = .false.
    end if
    if (maxval(abs(thlp - initial_p)) <= tol) then
      call report_scalar('lmoist off also silenced thlp', 0.)
      all_passed = .false.
    end if

    lmoist          = .true.
    lperiodicEBcorr = .false.
    call apply(unit_flux, unit_flux)
    if (maxval(abs(thlp - initial_p)) > tol .or. maxval(abs(qtp - initial_p)) > tol) then
      call report_scalar('disabled', max(maxval(abs(thlp - initial_p)), &
                                         maxval(abs(qtp - initial_p))))
      all_passed = .false.
    end if

    lperiodicEBcorr = saved_lperiodicEBcorr
    ltempeq         = saved_ltempeq
    lmoist          = saved_lmoist
    sinkbase        = saved_sinkbase
    fraction        = saved_fraction
    totheatflux     = saved_totheatflux
    totqflux        = saved_totqflux

    if (myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      if (all_passed) then
        write(*, '(A)') 'ALL TESTS PASSED: tests_periodic_ebcorr'
        write(*, '(A,I0,A)') '  Column closure, Grylls split, support and switches over ', &
                             nprocs, ' rank(s)'
      else
        write(*, '(A)') 'TESTS FAILED: tests_periodic_ebcorr'
        write(*, '(A)') '  One or more checks did not pass'
      end if
      write(*, '(A)') '================================================'
    end if

    tests_periodic_ebcorr = all_passed

  contains

    !> Reset the tendencies, set the per-rank fluxes and run the correction.
    subroutine apply(heat, moist)
      real, intent(in) :: heat, moist

      thlp = initial_p
      qtp  = initial_p
      totheatflux = heat
      totqflux    = moist
      call periodicEBcorr

    end subroutine apply

    !> Depth-weighted column integral of the tendency in an interior column.
    real function column_integral(field)
      real, intent(in) :: field(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
      integer :: kk

      column_integral = 0.
      do kk = kb, ke
        column_integral = column_integral + (field(ib,jb,kk) - initial_p) * dzf(kk)
      end do

    end function column_integral

    !> The column integral must be the domain-mean flux, at unit amplitude.
    logical function check_closure(base, frac, want)
      integer, intent(in) :: base
      real,    intent(in) :: frac, want

      check_closure = check_closure_amplitude(base, frac, 1., want)

    end function check_closure

    !> As check_closure, with the per-rank flux scaled by amplitude.
    logical function check_closure_amplitude(base, frac, amplitude, want)
      integer, intent(in) :: base
      real,    intent(in) :: frac, amplitude, want
      real :: got, ctol

      check_closure_amplitude = .true.
      sinkbase = base
      fraction = frac
      call apply(amplitude * unit_flux * real(myid + 1), &
                 amplitude * unit_flux * real(myid + 1))
      ctol = 1.e-12 * max(1., abs(want))

      got = column_integral(thlp)
      if (abs(got - want) > ctol) then
        write(*,'(A,I0,A,I0,A,F6.3,A,ES23.15,A,ES23.15)') 'FAIL on rank ', myid, &
             ': periodicEBcorr heat closure at sinkbase ', base, ', fraction ', frac, &
             ': integral ', got, ' expected ', want
        check_closure_amplitude = .false.
      end if

      got = column_integral(qtp)
      if (abs(got - want) > ctol) then
        write(*,'(A,I0,A,I0,A,F6.3,A,ES23.15,A,ES23.15)') 'FAIL on rank ', myid, &
             ': periodicEBcorr moisture closure at sinkbase ', base, ', fraction ', frac, &
             ': integral ', got, ' expected ', want
        check_closure_amplitude = .false.
      end if

    end function check_closure_amplitude

    !> The top level must carry (1-fraction) of the total on its own.
    !!
    !! Both ke and ke-1 receive the volume sink, so their difference isolates
    !! the top-level term without reference to how the routine forms it.
    logical function check_split(base, frac)
      integer, intent(in) :: base
      real,    intent(in) :: frac
      real :: got, want, stol

      check_split = .true.
      sinkbase = base
      fraction = frac
      call apply(unit_flux * real(myid + 1), unit_flux * real(myid + 1))

      want = (1. - frac) * mean_flux
      stol = 1.e-12 * max(1., abs(mean_flux))

      got = (thlp(ib,jb,ke) - thlp(ib,jb,ke-1)) * dzf(ke)
      if (abs(got - want) > stol) then
        write(*,'(A,I0,A,F6.3,A,ES23.15,A,ES23.15)') 'FAIL on rank ', myid, &
             ': periodicEBcorr heat top-level share at fraction ', frac, &
             ': got ', got, ' expected ', want
        check_split = .false.
      end if

      got = (qtp(ib,jb,ke) - qtp(ib,jb,ke-1)) * dzf(ke)
      if (abs(got - want) > stol) then
        write(*,'(A,I0,A,F6.3,A,ES23.15,A,ES23.15)') 'FAIL on rank ', myid, &
             ': periodicEBcorr moisture top-level share at fraction ', frac, &
             ': got ', got, ' expected ', want
        check_split = .false.
      end if

    end function check_split

    !> Assert the lateral halo faces still hold their initial value.
    logical function check_untouched(label, field)
      character(len=*), intent(in) :: label
      real, intent(in) :: field(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
      real :: worst

      check_untouched = .true.
      worst = 0.

      if (ih > 0) then
        worst = max(worst, maxval(abs(field(ib-ih:ib-1,   :, kb:ke) - initial_p)))
        worst = max(worst, maxval(abs(field(ie+1:ie+ih,   :, kb:ke) - initial_p)))
      end if
      if (jh > 0) then
        worst = max(worst, maxval(abs(field(:, jb-jh:jb-1, kb:ke) - initial_p)))
        worst = max(worst, maxval(abs(field(:, je+1:je+jh, kb:ke) - initial_p)))
      end if

      if (worst > tol) then
        call report_scalar(label, worst)
        check_untouched = .false.
      end if

    end function check_untouched

    !> Print one failing quantity with its rank.
    subroutine report_scalar(label, value)
      character(len=*), intent(in) :: label
      real,             intent(in) :: value

      write(*,'(A,I0,A,A,A,ES23.15)') 'FAIL on rank ', myid, ': periodicEBcorr ', &
           trim(label), ' = ', value

    end subroutine report_scalar

  end function tests_periodic_ebcorr


  !> Check the flow-rate controllers against the rate they are supposed to hit.
  !!
  !! masscorr diagnoses how far the predicted velocity is from the prescribed
  !! flow rate and adds a uniform defect to the tendency to close the gap. The
  !! diagnosis is spread over three quantities the routine forms internally
  !! (uoutflow, uflowrateold, udef), so comparing against those would restate
  !! the code. What follows instead measures the rate the routine claims to
  !! set - independently, with its own reduction - and asserts:
  !!
  !!   - after the call, the achieved rate equals uflowrate (or vflowrate).
  !!     This is the whole contract, and it holds for every branch;
  !!   - a second call is a no-op, because the first one already met the
  !!     target. This catches a defect that is right once but not consistent
  !!     with what the correction actually applies;
  !!   - the rate is the global one, so each rank is given a different field
  !!     and every rank is checked against the same answer;
  !!   - u and v are independent, the correction is uniform over ib:ie, jb:je,
  !!     kb:ke and touches nothing else, and linoutflow or all four switches
  !!     off is a genuine no-op.
  !!
  !! Two mask phases. The volume branches average over the fluid cells, so they
  !! run first against the real IBM masks. The outflow branches divide by an
  !! outlet area that calcfluidvolumes derives from IIc while the branch itself
  !! masks with IIu or IIv; the two agree only when the outlet plane is clear,
  !! so the second phase sets every mask to one and recomputes the areas. That
  !! is the configuration those branches are written for - "Assumes ie=itot",
  !! with no mention of blocks at the outlet.
  !!
  !! This exercises the host branch. The device kernels are covered by
  !! tests_cuda.f90::test_masscorr, which runs under UDALES_RUN_CUDA_SELFTEST
  !! on a Debug GPU build.
  !! Returns .true. if all checks pass, .false. otherwise.
  logical function tests_masscorr()
    use modglobal, only : runmode, ib, ie, ih, ihc, jb, je, jh, jhc, kb, ke, kh, khc, &
                          libm, dy, dxf, dzf, zh, rslabs, &
                          linoutflow, luoutflowr, lvoutflowr, luvolflowr, lvvolflowr, &
                          uflowrate, vflowrate, rk3coef, rk3coefi
    use modfields, only : initfields, um, up, vm, vp, udef, vdef, uouttot, &
                          uoutarea, voutarea, IIc, IIu, IIv, IIcs, IIus, IIvs
    use modforces, only : masscorr, calcfluidvolumes
    use initfac,   only : readfacetfiles
    use modibm,    only : initibm, createmasks
    use modmpi,    only : nprocy, MPI_SUM

    implicit none

    real, parameter :: rk3_test = 3., target_u = 1.25, target_v = -0.75

    logical :: all_passed
    integer :: k
    real    :: tol, got
    real, allocatable :: up_save(:,:,:), vp_save(:,:,:)
    integer, allocatable :: IIc_save(:,:,:), IIu_save(:,:,:), IIv_save(:,:,:)
    integer, allocatable :: IIcs_save(:), IIus_save(:), IIvs_save(:)

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_masscorr: FLOW RATE CORRECTION TEST'
      write(*, '(A)') '------------------------------------------------'
    end if

#if defined(_GPU)
    ! masscorr writes the device-resident tendencies in a GPU build, so the
    ! host arrays this test inspects would never change. Say so explicitly
    ! rather than reporting an unexplained numeric mismatch, and do not pass
    ! vacuously.
    if (myid == 0) then
      write(*, '(A)') 'FAIL: runmode 1010 exercises the host branch of masscorr.'
      write(*, '(A)') '  Run it against a CPU build (build/cpu/<type>/u-dales).'
      write(*, '(A)') '  The device kernels are covered by tests_cuda.f90::test_masscorr,'
      write(*, '(A)') '  run with UDALES_RUN_CUDA_SELFTEST=1 on a Debug GPU build.'
    end if
    tests_masscorr = .false.
    return
#endif

    call initfields
    if (libm) then
      call readfacetfiles
      call initibm
    end if
    call createmasks
    call calcfluidvolumes

    all_passed = .true.
    tol        = 1.e-10

    ! The volume branches divide each level by its fluid-cell count, and
    ! avexy_ibm returns a sentinel where that count is zero. A fully solid
    ! level would put -999. into the target, so say so rather than chase the
    ! resulting mismatch.
    do k = kb, ke
      if (IIus(k) == 0 .or. IIvs(k) == 0) then
        if (myid == 0) write(*, '(A,I0,A)') &
          'FAIL: level ', k, ' is entirely solid; this runmode needs a case with none'
        tests_masscorr = .false.
        return
      end if
    end do

    rk3coef  = rk3_test
    rk3coefi = 1. / rk3_test
    uflowrate = target_u
    vflowrate = target_v
    linoutflow = .false.

    allocate(up_save(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(vp_save(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))

    ! ---------------------------------------------------------------- 1
    ! Volume branches, against the real masks.
    call set_flags(.false., .true., .false., .false.)
    call fill_fields
    call masscorr
    got = volume_rate(um, up, IIu, IIus)
    if (abs(got - target_u) > tol) then
      call report('u volume rate not reached', got - target_u)
      all_passed = .false.
    end if
    ! The v tendency must not have moved: the u branch owns only u.
    if (maxval(abs(vp - vp_save)) > tol) then
      call report('u volume branch moved vp', maxval(abs(vp - vp_save)))
      all_passed = .false.
    end if
    ! Support: outside ib:ie, jb:je, kb:ke nothing is touched.
    if (.not. check_support('u volume', up, up_save)) all_passed = .false.
    ! Uniform: the same defect everywhere in the interior.
    if (.not. check_uniform('u volume', up, up_save)) all_passed = .false.

    ! ---------------------------------------------------------------- 2
    ! Idempotent: the target is already met, so the second call adds nothing.
    up_save = up
    call masscorr
    if (abs(udef) > tol) then
      call report('u volume defect not zero on the second call', udef)
      all_passed = .false.
    end if
    if (maxval(abs(up - up_save)) > tol) then
      call report('u volume second call moved up', maxval(abs(up - up_save)))
      all_passed = .false.
    end if

    ! ---------------------------------------------------------------- 3
    call set_flags(.false., .false., .false., .true.)
    call fill_fields
    call masscorr
    got = volume_rate(vm, vp, IIv, IIvs)
    if (abs(got - target_v) > tol) then
      call report('v volume rate not reached', got - target_v)
      all_passed = .false.
    end if
    if (maxval(abs(up - up_save)) > tol) then
      call report('v volume branch moved up', maxval(abs(up - up_save)))
      all_passed = .false.
    end if
    if (.not. check_support('v volume', vp, vp_save)) all_passed = .false.
    if (.not. check_uniform('v volume', vp, vp_save)) all_passed = .false.

    vp_save = vp
    call masscorr
    if (abs(vdef) > tol) then
      call report('v volume defect not zero on the second call', vdef)
      all_passed = .false.
    end if

    ! ---------------------------------------------------------------- 4
    ! Outflow branches, with the outlet planes clear. Save the masks first so
    ! the switch checks below still run against something realistic.
    allocate(IIc_save,  source=IIc)
    allocate(IIu_save,  source=IIu)
    allocate(IIv_save,  source=IIv)
    allocate(IIcs_save, source=IIcs)
    allocate(IIus_save, source=IIus)
    allocate(IIvs_save, source=IIvs)
    IIc = 1 ; IIu = 1 ; IIv = 1
    IIcs = nint(rslabs) ; IIus = nint(rslabs) ; IIvs = nint(rslabs)
    call calcfluidvolumes

    call set_flags(.true., .false., .false., .false.)
    call fill_fields
    call masscorr
    got = outlet_rate_u(um, up)
    if (abs(got - target_u) > tol) then
      call report('u outflow rate not reached', got - target_u)
      all_passed = .false.
    end if
    if (.not. check_support('u outflow', up, up_save)) all_passed = .false.
    if (.not. check_uniform('u outflow', up, up_save)) all_passed = .false.

    ! uouttot is read by the outflow boundary conditions and is the outlet mass
    ! flow of the tendency before the correction, so it is checked against the
    ! saved field rather than the corrected one.
    got = rk3_test * plane_sum_u(up_save) 
    if (abs(uouttot - got) > tol * max(1., abs(got))) then
      call report('uouttot', uouttot - got)
      all_passed = .false.
    end if

    up_save = up
    call masscorr
    if (abs(udef) > tol) then
      call report('u outflow defect not zero on the second call', udef)
      all_passed = .false.
    end if

    ! ---------------------------------------------------------------- 5
    ! The v outlet is a row of constant j, so it lives on one rank once the
    ! domain is split in y and the branch, which all-reduces over every rank,
    ! would add up one row per rank. That is a property of the scheme, so the
    ! check runs only where the row is whole.
    if (nprocy == 1) then
      call set_flags(.false., .false., .true., .false.)
      call fill_fields
      call masscorr
      got = outlet_rate_v(vm, vp)
      if (abs(got - target_v) > tol) then
        call report('v outflow rate not reached', got - target_v)
        all_passed = .false.
      end if
      if (.not. check_support('v outflow', vp, vp_save)) all_passed = .false.
      if (.not. check_uniform('v outflow', vp, vp_save)) all_passed = .false.

      vp_save = vp
      call masscorr
      if (abs(vdef) > tol) then
        call report('v outflow defect not zero on the second call', vdef)
        all_passed = .false.
      end if
    else if (myid == 0) then
      write(*, '(A)') '  skipped: v outflow branch needs nprocy = 1'
    end if

    IIc = IIc_save ; IIu = IIu_save ; IIv = IIv_save
    IIcs = IIcs_save ; IIus = IIus_save ; IIvs = IIvs_save
    call calcfluidvolumes

    ! ---------------------------------------------------------------- 6
    ! Switches. Every branch off, and linoutflow on with every branch set.
    call set_flags(.false., .false., .false., .false.)
    call fill_fields
    call masscorr
    if (maxval(abs(up - up_save)) > tol .or. maxval(abs(vp - vp_save)) > tol) then
      call report('all switches off', max(maxval(abs(up - up_save)), &
                                          maxval(abs(vp - vp_save))))
      all_passed = .false.
    end if

    call set_flags(.true., .true., .true., .true.)
    linoutflow = .true.
    call fill_fields
    call masscorr
    if (maxval(abs(up - up_save)) > tol .or. maxval(abs(vp - vp_save)) > tol) then
      call report('linoutflow', max(maxval(abs(up - up_save)), &
                                    maxval(abs(vp - vp_save))))
      all_passed = .false.
    end if
    linoutflow = .false.
    call set_flags(.false., .false., .false., .false.)

    deallocate(up_save, vp_save)
    deallocate(IIc_save, IIu_save, IIv_save, IIcs_save, IIus_save, IIvs_save)

    if (myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      if (all_passed) then
        write(*, '(A)') 'ALL TESTS PASSED: tests_masscorr'
        write(*, '(A)') '  Volume and outflow branches, both directions, over the fluid cells'
      else
        write(*, '(A)') 'TESTS FAILED: tests_masscorr'
        write(*, '(A)') '  One or more checks did not pass'
      end if
      write(*, '(A)') '================================================'
    end if

    tests_masscorr = all_passed

  contains

    subroutine set_flags(uout, uvol, vout, vvol)
      logical, intent(in) :: uout, uvol, vout, vvol

      luoutflowr = uout
      luvolflowr = uvol
      lvoutflowr = vout
      lvvolflowr = vvol

    end subroutine set_flags

    !> Give every rank a different field, so a missing all-reduce shows up.
    !!
    !! zstart(2) offsets the rank-local j onto the global grid, so the pattern
    !! is a single global field split across ranks rather than the same field
    !! repeated - the latter would let a rank-local reduction pass.
    subroutine fill_fields
      integer :: i, j, kk

      do kk = kb-kh, ke+kh
        do j = jb-jh, je+jh
          do i = ib-ih, ie+ih
            um(i,j,kk) = 0.5 + 0.01*real(i) - 0.003*real(j + zstart(2) - 1) + 0.007*real(kk)
            vm(i,j,kk) = -0.2 + 0.004*real(i) + 0.006*real(j + zstart(2) - 1) - 0.002*real(kk)
          end do
        end do
      end do
      do kk = kb, ke+kh
        do j = jb-jh, je+jh
          do i = ib-ih, ie+ih
            up(i,j,kk) = 0.02*real(i) + 0.05*real(j + zstart(2) - 1) - 0.01*real(kk)
            vp(i,j,kk) = -0.03*real(i) + 0.01*real(j + zstart(2) - 1) + 0.04*real(kk)
          end do
        end do
      end do

      up_save = up
      vp_save = vp

    end subroutine fill_fields

    !> Fluid-volume mean of m + rk3coef*p, the rate the volume branch targets.
    !!
    !! Formed here from the definition rather than by calling avexy_ibm, so a
    !! change of sign or of weighting inside masscorr cannot cancel out.
    real function volume_rate(m, p, II, IIs)
      real,    intent(in) :: m(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      real,    intent(in) :: p(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
      integer, intent(in) :: II(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc)
      integer, intent(in) :: IIs(kb:ke+khc)
      real    :: local(kb:ke), total(kb:ke)
      integer :: i, j, kk

      local = 0.
      do kk = kb, ke
        do j = jb, je
          do i = ib, ie
            local(kk) = local(kk) + (m(i,j,kk) + rk3coef*p(i,j,kk)) * II(i,j,kk)
          end do
        end do
      end do

      call MPI_ALLREDUCE(local, total, ke-kb+1, MY_REAL, MPI_SUM, comm3d, mpierr)

      volume_rate = 0.
      do kk = kb, ke
        volume_rate = volume_rate + (total(kk) / IIs(kk)) * dzf(kk)
      end do
      volume_rate = volume_rate / zh(ke+1)

    end function volume_rate

    !> Outlet-plane mean of m + rk3coef*p for the u branch.
    real function outlet_rate_u(m, p)
      real, intent(in) :: m(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      real, intent(in) :: p(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
      real    :: local(kb:ke), total(kb:ke)
      integer :: j, kk

      local = 0.
      do kk = kb, ke
        do j = jb, je
          local(kk) = local(kk) + (m(ie,j,kk) + rk3coef*p(ie,j,kk)) * IIu(ie,j,kk) * dy
        end do
      end do

      call MPI_ALLREDUCE(local, total, ke-kb+1, MY_REAL, MPI_SUM, comm3d, mpierr)

      outlet_rate_u = 0.
      do kk = kb, ke
        outlet_rate_u = outlet_rate_u + total(kk) * dzf(kk)
      end do
      outlet_rate_u = outlet_rate_u / uoutarea

    end function outlet_rate_u

    !> Outlet-plane sum of one field for the u branch, without the area divide.
    real function plane_sum_u(p)
      real, intent(in) :: p(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
      real    :: local(kb:ke), total(kb:ke)
      integer :: j, kk

      local = 0.
      do kk = kb, ke
        do j = jb, je
          local(kk) = local(kk) + p(ie,j,kk) * IIu(ie,j,kk) * dy
        end do
      end do

      call MPI_ALLREDUCE(local, total, ke-kb+1, MY_REAL, MPI_SUM, comm3d, mpierr)

      plane_sum_u = 0.
      do kk = kb, ke
        plane_sum_u = plane_sum_u + total(kk) * dzf(kk)
      end do

    end function plane_sum_u

    !> Outlet-row mean of m + rk3coef*p for the v branch.
    real function outlet_rate_v(m, p)
      real, intent(in) :: m(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      real, intent(in) :: p(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
      real    :: local(kb:ke), total(kb:ke)
      integer :: i, kk

      local = 0.
      do kk = kb, ke
        do i = ib, ie
          local(kk) = local(kk) + (m(i,je,kk) + rk3coef*p(i,je,kk)) * IIv(i,je,kk) * dxf(1)
        end do
      end do

      call MPI_ALLREDUCE(local, total, ke-kb+1, MY_REAL, MPI_SUM, comm3d, mpierr)

      outlet_rate_v = 0.
      do kk = kb, ke
        outlet_rate_v = outlet_rate_v + total(kk) * dzf(kk)
      end do
      outlet_rate_v = outlet_rate_v / voutarea

    end function outlet_rate_v

    !> Assert that nothing outside ib:ie, jb:je, kb:ke moved.
    logical function check_support(label, got, before)
      character(len=*), intent(in) :: label
      real, intent(in) :: got   (ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
      real, intent(in) :: before(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
      real    :: worst
      integer :: i, j, kk

      check_support = .true.
      worst = 0.
      do kk = kb, ke+kh
        do j = jb-jh, je+jh
          do i = ib-ih, ie+ih
            if (i < ib .or. i > ie .or. j < jb .or. j > je .or. kk > ke) then
              worst = max(worst, abs(got(i,j,kk) - before(i,j,kk)))
            end if
          end do
        end do
      end do

      if (worst > tol) then
        call report(label//' correction reached outside the interior', worst)
        check_support = .false.
      end if

    end function check_support

    !> Assert the interior moved by one and the same amount, and that it moved.
    logical function check_uniform(label, got, before)
      character(len=*), intent(in) :: label
      real, intent(in) :: got   (ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
      real, intent(in) :: before(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
      real    :: delta, worst
      integer :: i, j, kk

      check_uniform = .true.
      delta = got(ib,jb,kb) - before(ib,jb,kb)
      worst = 0.
      do kk = kb, ke
        do j = jb, je
          do i = ib, ie
            worst = max(worst, abs((got(i,j,kk) - before(i,j,kk)) - delta))
          end do
        end do
      end do

      if (worst > tol) then
        call report(label//' correction is not uniform', worst)
        check_uniform = .false.
      end if
      if (abs(delta) <= tol) then
        call report(label//' correction was zero, so nothing was tested', delta)
        check_uniform = .false.
      end if

    end function check_uniform

    !> Print one failing quantity with its rank.
    subroutine report(label, value)
      character(len=*), intent(in) :: label
      real,             intent(in) :: value

      write(*,'(A,I0,A,A,A,ES23.15)') 'FAIL on rank ', myid, ': masscorr ', &
           trim(label), ' = ', value

    end subroutine report

  end function tests_masscorr


  !> Check the host routines ibmnorm is built from against the conditions they
  !! are supposed to impose.
  !!
  !! ibmnorm does three things: it pins the velocity and its tendency to zero
  !! inside the solid, it pins the scalars to a fixed interior value except at
  !! solid cells with fluid neighbours, where it imposes zero flux instead, and
  !! it removes the advective flux the second-order scheme sent across solid
  !! faces. Each is stated as a property here rather than compared against the
  !! expression that produces it:
  !!
  !!   - a solid velocity point ends at exactly zero, and nothing else moves;
  !!   - a solid scalar point walled in on all six sides takes the interior
  !!     value, and one with fluid neighbours takes their mean instead - so a
  !!     constant field stays constant across the boundary, which is what zero
  !!     flux means, and the interior value cannot influence it;
  !!   - the liberal advection correction vanishes on a constant field, for any
  !!     velocity, because that is the flux it exists to cancel;
  !!   - both corrections vanish when the velocity does, are linear in the
  !!     scalar, and move only boundary points.
  !!
  !! The parallel port relies on the point lists holding no cell twice, since a
  !! repeat would be two threads writing one cell, so that is checked too.
  !!
  !! This exercises the host branch. The device kernels are covered by
  !! tests_cuda.f90::test_ibmnorm, which runs under UDALES_RUN_CUDA_SELFTEST on
  !! a Debug GPU build.
  !! Returns .true. if all checks pass, .false. otherwise.
  logical function tests_ibmnorm()
    use modglobal, only : runmode, ib, ie, ih, jb, je, jh, kb, ke, kh, &
                          libm, eps1, dxi, dyi, dxi5, dyi5, dzf, dzhi, dzfi, dzfi5
    use modfields, only : initfields, u0, v0, w0, thl0, thlp, thlm
    use initfac,   only : readfacetfiles
    use modibm,    only : initibm, createmasks, mask_c, mask_u, mask_v, mask_w, solid, &
                          advecc2nd_corr_conservative, advecc2nd_corr_liberal, &
                          solid_info_type, solid_info_u, solid_info_c, &
                          bound_info_c

    implicit none

    real, parameter :: interior = 300., constant_field = 7.5

    logical :: all_passed
    integer :: i, j, k, n, nsol, buried
    real    :: worst, tol
    type(solid_info_type) :: si
    real, allocatable :: mask_save(:,:,:), thlp_save(:,:,:)
    real, allocatable :: corr_a(:,:,:), corr_b(:,:,:)

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_ibmnorm: HOST IBM NORMAL-VELOCITY TESTS'
      write(*, '(A)') '------------------------------------------------'
    end if

#if defined(_GPU)
    ! The host routines this checks are not compiled into a GPU release build,
    ! and on a GPU debug build the solver does not use them, so there is
    ! nothing here to exercise. The device port is covered by
    ! tests_cuda.f90::test_ibmnorm.
    if (myid == 0) then
      write(*, '(A)') 'FAIL: runmode 1011 checks the host IBM normal-velocity routines.'
      write(*, '(A)') '  Run it against a CPU build (build/cpu/<type>/u-dales).'
      write(*, '(A)') '  The device port is covered by the CUDA self-tests.'
    end if
    tests_ibmnorm = .false.
    return
#endif

    if (.not. libm) then
      if (myid == 0) write(*, '(A)') 'FAIL: this runmode needs a case with libm = .true.'
      tests_ibmnorm = .false.
      return
    end if

    call initfields
    call readfacetfiles
    call initibm
    call createmasks

    all_passed = .true.
    tol = 1.e-12

    ! ---------------------------------------------------------------- 1
    ! No cell may appear twice in a solid list. The device port runs one
    ! thread per point, so a repeat would be two threads writing one cell -
    ! harmless in a sequential loop, a race in a parallel one.
    if (.not. check_distinct('solid_info_u', solid_info_u)) all_passed = .false.
    if (.not. check_distinct('solid_info_c', solid_info_c)) all_passed = .false.

    ! And the mask has to agree with the list: a cell is solid exactly when it
    ! is listed. The masked branch reads its neighbours in place, so the host's
    ! sequential loop and a one-thread-per-point port agree only while no
    ! listed cell is another's neighbour - which this invariant guarantees.
    worst = 0.
    do n = 1, solid_info_c%nsolptsrank
      i = solid_info_c%solpts_loc(n,1)
      j = solid_info_c%solpts_loc(n,2)
      k = solid_info_c%solpts_loc(n,3)
      worst = max(worst, abs(mask_c(i,j,k)))
    end do
    if (worst > 0.) then
      call report('a listed solid point is not masked solid', worst)
      all_passed = .false.
    end if

    ! ---------------------------------------------------------------- 2
    ! Velocity: solid points go to exactly zero, and nothing else moves.
    call fill(u0, 1.); call fill(thlp, 2.)
    allocate(thlp_save, source=thlp)
    call solid(solid_info_u, u0, thlp, 0., ih, jh, kh)

    worst = 0.
    do n = 1, solid_info_u%nsolptsrank
      i = solid_info_u%solpts_loc(n,1)
      j = solid_info_u%solpts_loc(n,2)
      k = solid_info_u%solpts_loc(n,3)
      worst = max(worst, abs(u0(i,j,k)), abs(thlp(i,j,k)))
    end do
    if (worst > 0.) then
      call report('solid velocity point not exactly zero', worst)
      all_passed = .false.
    end if
    if (.not. check_local('solid velocity', thlp, thlp_save, solid_info_u)) all_passed = .false.
    deallocate(thlp_save)

    ! ---------------------------------------------------------------- 3
    ! Scalars, on a mask the test owns so that a walled-in point exists. The
    ! interior of the block is solid, and (ib+1,jb+1,kb+1) sits inside it with
    ! all six neighbours solid as well.
    allocate(mask_save, source=mask_c)
    call build_block(si, buried)

    if (buried == 0) then
      if (myid == 0) write(*, '(A)') 'FAIL: could not place a walled-in point on this rank'
      tests_ibmnorm = .false.
      return
    end if

    ! A constant field must survive: every solid point with a fluid neighbour
    ! takes their mean, which for a constant field is the constant. This is
    ! what zero flux across the boundary means.
    thlm = constant_field
    thlp = 0.
    call solid(si, thlm, thlp, interior, ih, jh, kh, mask_c)

    worst = 0.
    do n = 1, si%nsolptsrank
      i = si%solpts_loc(n,1)
      j = si%solpts_loc(n,2)
      k = si%solpts_loc(n,3)
      if (n == buried) cycle
      if (fluid_neighbours(i,j,k) > 0) worst = max(worst, abs(thlm(i,j,k) - constant_field))
    end do
    if (worst > tol) then
      call report('constant scalar not preserved across the boundary', worst)
      all_passed = .false.
    end if

    ! The walled-in point has no fluid neighbour, so it takes the interior
    ! value instead - exactly, and with a zero tendency.
    i = si%solpts_loc(buried,1)
    j = si%solpts_loc(buried,2)
    k = si%solpts_loc(buried,3)
    if (abs(thlm(i,j,k) - interior) > 0. .or. abs(thlp(i,j,k)) > 0.) then
      call report('walled-in point did not take the interior value', &
                  max(abs(thlm(i,j,k) - interior), abs(thlp(i,j,k))))
      all_passed = .false.
    end if

    ! ---------------------------------------------------------------- 4
    ! A solid point with fluid neighbours takes their mean, so the interior
    ! value must not reach it. Run again with a different one and compare.
    thlm = constant_field
    thlp = 0.
    call solid(si, thlm, thlp, interior + 100., ih, jh, kh, mask_c)

    worst = 0.
    do n = 1, si%nsolptsrank
      i = si%solpts_loc(n,1)
      j = si%solpts_loc(n,2)
      k = si%solpts_loc(n,3)
      if (fluid_neighbours(i,j,k) > 0) worst = max(worst, abs(thlm(i,j,k) - constant_field))
    end do
    if (worst > tol) then
      call report('interior value leaked into a point with fluid neighbours', worst)
      all_passed = .false.
    end if

    ! And the mean itself, on a field that varies, worked out independently.
    call fill(thlm, 3.)
    thlp = 0.
    if (.not. check_neighbour_mean(si)) all_passed = .false.

    mask_c = mask_save
    deallocate(mask_save)
    if (allocated(si%solpts_loc)) deallocate(si%solpts_loc)

    ! ---------------------------------------------------------------- 5
    ! Advection corrections, against the case's own boundary points.
    if (bound_info_c%nbndptsrank < 1) then
      if (myid == 0) write(*, '(A)') 'FAIL: no cell-centred boundary points on this rank'
      tests_ibmnorm = .false.
      return
    end if

    allocate(corr_a(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(corr_b(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))

    ! No velocity, no flux to cancel: both variants must do nothing at all.
    u0 = 0.; v0 = 0.; w0 = 0.
    call fill(thl0, 4.)
    thlp = 0.
    call advecc2nd_corr_conservative(thl0, thlp)
    if (maxval(abs(thlp)) > 0.) then
      call report('conservative correction at zero velocity', maxval(abs(thlp)))
      all_passed = .false.
    end if
    thlp = 0.
    call advecc2nd_corr_liberal(thl0, thlp)
    if (maxval(abs(thlp)) > 0.) then
      call report('liberal correction at zero velocity', maxval(abs(thlp)))
      all_passed = .false.
    end if

    ! A constant scalar has no advective flux to cancel whatever the velocity
    ! is, because the liberal correction replaces the solid value with the
    ! fluid one and the two are equal.
    call fill(u0, 5.); call fill(v0, 6.); call fill(w0, 7.)
    thl0 = constant_field
    thlp = 0.
    call advecc2nd_corr_liberal(thl0, thlp)
    if (maxval(abs(thlp)) > tol) then
      call report('liberal correction on a constant field', maxval(abs(thlp)))
      all_passed = .false.
    end if

    ! Linear in the scalar: every term carries one factor of it.
    call fill(thl0, 8.)
    thlp = 0.
    call advecc2nd_corr_conservative(thl0, thlp)
    corr_a = thlp
    thl0 = 2. * thl0
    thlp = 0.
    call advecc2nd_corr_conservative(thl0, thlp)
    corr_b = thlp
    worst = maxval(abs(corr_b - 2.*corr_a))
    if (worst > 1.e-10 * max(1., maxval(abs(corr_a)))) then
      call report('conservative correction is not linear in the scalar', worst)
      all_passed = .false.
    end if
    if (maxval(abs(corr_a)) == 0.) then
      call report('conservative correction was zero, so nothing was tested', 0.)
      all_passed = .false.
    end if

    call fill(thl0, 8.)
    thlp = 0.
    call advecc2nd_corr_liberal(thl0, thlp)
    corr_a = thlp
    thl0 = 2. * thl0
    thlp = 0.
    call advecc2nd_corr_liberal(thl0, thlp)
    worst = maxval(abs(thlp - 2.*corr_a))
    if (worst > 1.e-10 * max(1., maxval(abs(corr_a)))) then
      call report('liberal correction is not linear in the scalar', worst)
      all_passed = .false.
    end if
    if (maxval(abs(corr_a)) == 0.) then
      call report('liberal correction was zero, so nothing was tested', 0.)
      all_passed = .false.
    end if

    ! Only boundary points move.
    if (.not. check_local_bnd('conservative', corr_b)) all_passed = .false.
    if (.not. check_local_bnd('liberal', corr_a)) all_passed = .false.

    ! ---------------------------------------------------------------- 6
    ! The value itself, from a closed form the routines do not use. The
    ! liberal variant replaces the solid neighbour with the cell's own value,
    ! so what it adds across a solid face is the velocity there times the jump
    ! the scheme saw - which collapses to one term per face instead of two.
    call fill(thl0, 8.)
    thlp = 0.
    call advecc2nd_corr_liberal(thl0, thlp)
    if (.not. check_liberal_value(thl0, thlp)) all_passed = .false.

    ! And the conservative variant removes the flux itself, so on a constant
    ! field it leaves that constant times the velocity divergence over the
    ! solid faces alone.
    thl0 = constant_field
    thlp = 0.
    call advecc2nd_corr_conservative(thl0, thlp)
    if (.not. check_conservative_value(constant_field, thlp)) all_passed = .false.

    deallocate(corr_a, corr_b)

    if (myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      if (all_passed) then
        write(*, '(A)') 'ALL TESTS PASSED: tests_ibmnorm'
        write(*, '(A)') '  Velocity pinning, the scalar boundary condition and both corrections'
      else
        write(*, '(A)') 'TESTS FAILED: tests_ibmnorm'
        write(*, '(A)') '  One or more checks did not pass'
      end if
      write(*, '(A)') '================================================'
    end if

    tests_ibmnorm = all_passed

  contains

    !> Fill a field with a value that varies in all three directions.
    subroutine fill(a, offset)
      real, intent(inout) :: a(:,:,:)
      real, intent(in)    :: offset
      integer :: ii, jj, kk

      do kk = 1, size(a,3)
        do jj = 1, size(a,2)
          do ii = 1, size(a,1)
            a(ii,jj,kk) = offset + 0.125*ii - 0.0625*jj + 0.03125*kk
          end do
        end do
      end do

    end subroutine fill

    !> Number of fluid neighbours of a cell, read off mask_c.
    integer function fluid_neighbours(ii, jj, kk)
      integer, intent(in) :: ii, jj, kk

      fluid_neighbours = 0
      if (abs(mask_c(ii,jj+1,kk) - 1.) < eps1) fluid_neighbours = fluid_neighbours + 1
      if (abs(mask_c(ii,jj-1,kk) - 1.) < eps1) fluid_neighbours = fluid_neighbours + 1
      if (abs(mask_c(ii,jj,kk+1) - 1.) < eps1) fluid_neighbours = fluid_neighbours + 1
      if (abs(mask_c(ii,jj,kk-1) - 1.) < eps1) fluid_neighbours = fluid_neighbours + 1
      if (abs(mask_c(ii+1,jj,kk) - 1.) < eps1) fluid_neighbours = fluid_neighbours + 1
      if (abs(mask_c(ii-1,jj,kk) - 1.) < eps1) fluid_neighbours = fluid_neighbours + 1

    end function fluid_neighbours

    !> Mark a 3x3x3 block solid and return the index of its walled-in centre.
    !!
    !! The case's own geometry need not contain a cell with six solid
    !! neighbours, and the interior value is only reachable through one, so the
    !! test builds a block rather than hoping to find one.
    subroutine build_block(info, centre)
      type(solid_info_type), intent(out) :: info
      integer, intent(out) :: centre
      integer :: i0, j0, k0, ii, jj, kk, m

      centre = 0
      i0 = (ib + ie)/2
      j0 = (jb + je)/2
      k0 = kb + 2
      if (i0-1 < ib .or. i0+1 > ie .or. j0-1 < jb .or. j0+1 > je .or. k0+1 > ke) return

      mask_c = 1.
      info%nsolptsrank = 27
      allocate(info%solpts_loc(27,3))
      m = 0
      do kk = k0-1, k0+1
        do jj = j0-1, j0+1
          do ii = i0-1, i0+1
            m = m + 1
            info%solpts_loc(m,1) = ii
            info%solpts_loc(m,2) = jj
            info%solpts_loc(m,3) = kk
            mask_c(ii,jj,kk) = 0.
            if (ii == i0 .and. jj == j0 .and. kk == k0) centre = m
          end do
        end do
      end do

    end subroutine build_block

    !> The mean over fluid neighbours, worked out here rather than read back.
    logical function check_neighbour_mean(info)
      type(solid_info_type), intent(in) :: info
      real, allocatable :: before(:,:,:)
      real    :: want, got, count, deviation
      integer :: m, ii, jj, kk

      check_neighbour_mean = .true.
      allocate(before, source=thlm)
      call solid(info, thlm, thlp, interior, ih, jh, kh, mask_c)

      deviation = 0.
      do m = 1, info%nsolptsrank
        ii = info%solpts_loc(m,1)
        jj = info%solpts_loc(m,2)
        kk = info%solpts_loc(m,3)
        count = fluid_neighbours(ii,jj,kk)
        if (count == 0) cycle
        want = 0.
        if (abs(mask_c(ii,jj+1,kk) - 1.) < eps1) want = want + before(ii,jj+1,kk)
        if (abs(mask_c(ii,jj-1,kk) - 1.) < eps1) want = want + before(ii,jj-1,kk)
        if (abs(mask_c(ii,jj,kk+1) - 1.) < eps1) want = want + before(ii,jj,kk+1)
        if (abs(mask_c(ii,jj,kk-1) - 1.) < eps1) want = want + before(ii,jj,kk-1)
        if (abs(mask_c(ii+1,jj,kk) - 1.) < eps1) want = want + before(ii+1,jj,kk)
        if (abs(mask_c(ii-1,jj,kk) - 1.) < eps1) want = want + before(ii-1,jj,kk)
        want = want / count
        got  = thlm(ii,jj,kk)
        deviation = max(deviation, abs(got - want))
      end do

      if (deviation > 1.e-10) then
        call report('solid point is not the mean of its fluid neighbours', deviation)
        check_neighbour_mean = .false.
      end if
      deallocate(before)

    end function check_neighbour_mean

    !> Assert no cell appears twice in a solid list.
    logical function check_distinct(label, info)
      character(len=*), intent(in) :: label
      type(solid_info_type), intent(in) :: info
      integer, allocatable :: seen(:,:,:)
      integer :: m

      check_distinct = .true.
      if (info%nsolptsrank <= 1) return

      allocate(seen(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
      seen = 0
      do m = 1, info%nsolptsrank
        if (seen(info%solpts_loc(m,1),info%solpts_loc(m,2),info%solpts_loc(m,3)) /= 0) then
          call report(label//' lists a cell twice', real(m))
          check_distinct = .false.
          exit
        end if
        seen(info%solpts_loc(m,1),info%solpts_loc(m,2),info%solpts_loc(m,3)) = m
      end do
      deallocate(seen)

    end function check_distinct

    !> Assert only the listed solid points differ from the saved field.
    logical function check_local(label, got, before, info)
      character(len=*), intent(in) :: label
      real, intent(in) :: got   (ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
      real, intent(in) :: before(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
      type(solid_info_type), intent(in) :: info
      real, allocatable :: diff(:,:,:)
      integer :: m

      check_local = .true.
      allocate(diff(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
      diff = got - before
      do m = 1, info%nsolptsrank
        if (info%solpts_loc(m,3) > ke+kh) cycle
        diff(info%solpts_loc(m,1),info%solpts_loc(m,2),info%solpts_loc(m,3)) = 0.
      end do
      if (maxval(abs(diff)) > 0.) then
        call report(label//' moved a cell that is not a solid point', maxval(abs(diff)))
        check_local = .false.
      end if
      deallocate(diff)

    end function check_local

    !> The liberal correction, one term per solid face rather than two.
    logical function check_liberal_value(var, got)
      real, intent(in) :: var(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      real, intent(in) :: got(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
      real    :: want, worst, scale
      integer :: m, ii, jj, kk

      check_liberal_value = .true.
      worst = 0.
      scale = 0.
      do m = 1, bound_info_c%nbndptsrank
        ii = bound_info_c%bndpts_loc(m,1)
        jj = bound_info_c%bndpts_loc(m,2)
        kk = bound_info_c%bndpts_loc(m,3)

        want = 0.
        if (abs(mask_c(ii+1,jj,kk)) < eps1) &
          want = want + u0(ii+1,jj,kk)*(var(ii+1,jj,kk) - var(ii,jj,kk))*dxi5
        if (abs(mask_c(ii-1,jj,kk)) < eps1) &
          want = want - u0(ii,jj,kk)*(var(ii-1,jj,kk) - var(ii,jj,kk))*dxi5
        if (abs(mask_c(ii,jj+1,kk)) < eps1) &
          want = want + v0(ii,jj+1,kk)*(var(ii,jj+1,kk) - var(ii,jj,kk))*dyi5
        if (abs(mask_c(ii,jj-1,kk)) < eps1) &
          want = want - v0(ii,jj,kk)*(var(ii,jj-1,kk) - var(ii,jj,kk))*dyi5
        if (abs(mask_c(ii,jj,kk+1)) < eps1) &
          want = want + w0(ii,jj,kk+1)*(var(ii,jj,kk+1) - var(ii,jj,kk))*dzf(kk)*dzhi(kk+1)*dzfi5(kk)
        if (abs(mask_c(ii,jj,kk-1)) < eps1) &
          want = want - w0(ii,jj,kk)*(var(ii,jj,kk-1) - var(ii,jj,kk))*dzf(kk)*dzhi(kk)*dzfi5(kk)

        worst = max(worst, abs(got(ii,jj,kk) - want))
        scale = max(scale, abs(want))
      end do

      if (worst > 1.e-10 * max(1., scale)) then
        call report('liberal correction does not match the flux it replaces', worst)
        check_liberal_value = .false.
      end if
      if (scale == 0.) then
        call report('liberal correction was zero, so nothing was tested', 0.)
        check_liberal_value = .false.
      end if

    end function check_liberal_value

    !> The conservative correction on a constant field: the constant times the
    !! velocity divergence taken over the solid faces alone.
    logical function check_conservative_value(c, got)
      real, intent(in) :: c
      real, intent(in) :: got(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
      real    :: want, worst, scale
      integer :: m, ii, jj, kk

      check_conservative_value = .true.
      worst = 0.
      scale = 0.
      do m = 1, bound_info_c%nbndptsrank
        ii = bound_info_c%bndpts_loc(m,1)
        jj = bound_info_c%bndpts_loc(m,2)
        kk = bound_info_c%bndpts_loc(m,3)

        want = 0.
        if ((abs(mask_u(ii+1,jj,kk)) < eps1) .or. (abs(mask_c(ii+1,jj,kk)) < eps1)) &
          want = want + c*u0(ii+1,jj,kk)*dxi
        if ((abs(mask_u(ii,jj,kk)) < eps1) .or. (abs(mask_c(ii-1,jj,kk)) < eps1)) &
          want = want - c*u0(ii,jj,kk)*dxi
        if ((abs(mask_v(ii,jj+1,kk)) < eps1) .or. (abs(mask_c(ii,jj+1,kk)) < eps1)) &
          want = want + c*v0(ii,jj+1,kk)*dyi
        if ((abs(mask_v(ii,jj,kk)) < eps1) .or. (abs(mask_c(ii,jj-1,kk)) < eps1)) &
          want = want - c*v0(ii,jj,kk)*dyi
        if ((abs(mask_w(ii,jj,kk+1)) < eps1) .or. (abs(mask_c(ii,jj,kk+1)) < eps1)) &
          want = want + c*w0(ii,jj,kk+1)*dzfi(kk)
        if ((abs(mask_w(ii,jj,kk)) < eps1) .or. (abs(mask_c(ii,jj,kk-1)) < eps1)) &
          want = want - c*w0(ii,jj,kk)*dzfi(kk)

        worst = max(worst, abs(got(ii,jj,kk) - want))
        scale = max(scale, abs(want))
      end do

      if (worst > 1.e-10 * max(1., scale)) then
        call report('conservative correction is not the solid-face flux', worst)
        check_conservative_value = .false.
      end if
      if (scale == 0.) then
        call report('conservative correction was zero, so nothing was tested', 0.)
        check_conservative_value = .false.
      end if

    end function check_conservative_value

    !> Assert only cell-centred boundary points carry a correction.
    logical function check_local_bnd(label, corr)
      character(len=*), intent(in) :: label
      real, intent(in) :: corr(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
      real, allocatable :: diff(:,:,:)
      integer :: m

      check_local_bnd = .true.
      allocate(diff(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
      diff = corr
      do m = 1, bound_info_c%nbndptsrank
        diff(bound_info_c%bndpts_loc(m,1),bound_info_c%bndpts_loc(m,2),bound_info_c%bndpts_loc(m,3)) = 0.
      end do
      if (maxval(abs(diff)) > 0.) then
        call report(label//' correction reached beyond the boundary points', maxval(abs(diff)))
        check_local_bnd = .false.
      end if
      deallocate(diff)

    end function check_local_bnd

    !> Print one failing quantity with its rank.
    subroutine report(label, value)
      character(len=*), intent(in) :: label
      real,             intent(in) :: value

      write(*,'(A,I0,A,A,A,ES23.15)') 'FAIL on rank ', myid, ': ibmnorm ', &
           trim(label), ' = ', value

    end subroutine report

  end function tests_ibmnorm

end module tests
