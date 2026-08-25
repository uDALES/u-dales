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
            tests_periodic_ebcorr, tests_masscorr, tests_ibmnorm, tests_eb, &
            tests_vegetation, tests_checksim, tests_driver_planes, &
            tests_thermodynamics, tests_tstep, tests_timedep

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
                          zhat, totheatflux, totqflux
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
    integer :: i, j, k, n, buried
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


  !> Properties of the facet energy balance, modEB::EB and its helpers.
  !!
  !! EB is the one routine in this part of the time loop that deliberately
  !! stays on the host, so there is no device branch to compare against and no
  !! symmetric test to write. What is checked here instead are properties the
  !! routine has to hold on its own, chosen so that each one fails for a
  !! different reason:
  !!
  !!   - intqH keeps the facet fluxes on the third Runge-Kutta stage only, and
  !!     clears them on every stage. This is the contract that lets modcuda's
  !!     updateFacFluxHost skip two crossings in three, so it is checked here
  !!     rather than assumed there. Over two ranks it also shows the flux is
  !!     the reduced one and the time integral lives on rank 0 alone.
  !!   - The energy balance fires only on the third stage and only once the
  !!     next EB time has arrived, and sets lfacetprops_dirty exactly when it
  !!     rewrites the facet properties. That flag is what the GPU mirror now
  !!     follows, so a run that never sets it would leave the wall functions on
  !!     initial-condition temperatures for the whole simulation.
  !!   - A facet whose surface is in equilibrium does not move. Uniform layer
  !!     temperatures with the net surface flux exactly balancing the emitted
  !!     longwave is a fixed point of the scheme, and almost any slip in the
  !!     matrix assembly, the inverse or the two matmuls destroys it.
  !!   - The surface energy balance closes: the conducted flux out of the
  !!     surface equals what went in, minus what was radiated. This uses only
  !!     faclam, facTdash and the forcing the test itself set, so it is the
  !!     physics rather than a second copy of the algebra.
  !!   - The innermost layer is held fixed, which is what the last row of the
  !!     system is for.
  !!   - calclw agrees with itself: the sparse and dense view-factor paths are
  !!     two separate loops that must produce the same incoming longwave.
  !!
  !! Every check has a negative control that deliberately breaks it, and the
  !! controls are themselves asserted to fire - a control that cannot fail
  !! proves nothing about the check beside it.
  logical function tests_eb()
    use modglobal, only : runmode, nfcts, nfaclyrs, lEB, lvfsparse, nnz, boltz, &
                          timee, tnextEB, tEB, rk3step, dt, cp, rhoa, &
                          lwriteEBfiles, lconstW
    use modEB,     only : EB, initEB, intqH, calclw
    use initfac,   only : readfacetfiles
    use initfac,   only : facT, facTdash, faclam, faccp, facd, facem, faca, &
                          facets, faclGR, netsw, facLWin, svf, vf, vfsparse, &
                          ivfsparse, jvfsparse, fachf, facef, fachfi, facefi, &
                          fachfsum, facefsum, facwsoil, lfacetprops_dirty
    use modmpi,    only : nprocs

    implicit none

    real, parameter :: T0     = 300.        !< uniform starting temperature [K]
    real, parameter :: probeH = 0.375       !< per-rank sensible flux probe
    real, parameter :: probeE = 0.125       !< per-rank latent flux probe
    real, parameter :: dt_lcl = 0.25        !< time step used for the integral
    real, parameter :: span   = 20.         !< seconds between energy balances
    real, parameter :: extraQ = 250.        !< net heating on top of balance
    real, parameter :: H0     = 40.         !< sensible flux fed through fachfi
    real, parameter :: E0     = 15.         !< latent flux fed through facefi

    logical :: all_passed, control_fired, layered
    integer :: n, j, stage, nlive
    real    :: rank_sum, want, got, tol, resid, scale, worst
    real    :: moved, lhs, rhs, dsum, dshare, ab_n, bb1

    ! Saved state. EB is destructive: it rewrites the temperatures, resets the
    ! integrated fluxes and advances the energy balance clock.
    real,    allocatable :: facT_s(:,:), facTdash_s(:,:), netsw_s(:), svf_s(:)
    real,    allocatable :: facLWin_s(:), fachfi_s(:), facefi_s(:), facwsoil_s(:)
    real,    allocatable :: fachf_s(:), facef_s(:), vfsparse_s(:)
    real,    allocatable :: faclam_s(:,:), faccp_s(:,:), facd_s(:,:)
    logical, allocatable :: faclGR_s(:)
    real,    allocatable :: facT_old(:,:), lw_sparse(:), lw_dense(:)
    real,    allocatable :: told(:), tdold(:)
    real                 :: timee_s, tnextEB_s, tEB_s, dt_s
    integer              :: rk3step_s
    logical              :: lvfsparse_s, lwriteEBfiles_s, lconstW_s

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_eb: FACET ENERGY BALANCE TEST'
      write(*, '(A)') '------------------------------------------------'
    end if

    if (.not. lEB) then
      if (myid == 0) then
        write(*, '(A)') 'FAIL: this runmode needs lEB = .true. in the namelist.'
        write(*, '(A)') '  Without it EB returns immediately and every check below is vacuous.'
      end if
      tests_eb = .false.
      return
    end if

    ! The runmode dispatch happens before the facet files are read and before
    ! initEB builds the layer system, so this entry point does both itself.
    ! The netCDF side is switched off first: writing facT.xxx.nc and
    ! facEB.xxx.nc from a test would leave meaningless records behind, and
    ! none of the properties below depend on the output.
    lwriteEBfiles = .false.
    call readfacetfiles
    call initEB

    if (nfcts < 1 .or. .not. allocated(facT)) then
      if (myid == 0) write(*, '(A)') 'FAIL: no facets to test'
      tests_eb = .false.
      return
    end if

    all_passed = .true.

    ! The radiative environment, the integrated fluxes and the temperature
    ! gradient are allocated on rank 0 alone - initfac keeps them there because
    ! that is where the energy balance is solved. Everything the wall functions
    ! read, facT among it, exists on every rank. So the setup below happens on
    ! rank 0, EB itself is called collectively for its broadcasts, and the
    ! checks split accordingly.
    allocate(facT_s,   source=facT)
    allocate(fachf_s,  source=fachf)
    allocate(facef_s,  source=facef)
    allocate(faclGR_s, source=faclGR)
    allocate(facT_old(1:nfcts,1:nfaclyrs+1))
    allocate(told(1:nfaclyrs+1), tdold(1:nfaclyrs+1))
    allocate(fachfi_s, source=fachfi)
    allocate(facefi_s, source=facefi)
    if (myid == 0) then
      allocate(facTdash_s, source=facTdash)
      allocate(netsw_s,    source=netsw)
      allocate(svf_s,      source=svf)
      allocate(facLWin_s,  source=facLWin)
      allocate(facwsoil_s, source=facwsoil)
      allocate(faclam_s,   source=faclam)
      allocate(faccp_s,    source=faccp)
      allocate(facd_s,     source=facd)
      if (allocated(vfsparse)) allocate(vfsparse_s, source=vfsparse)
    end if

    ! Give the wall a real profile before anything else runs.
    !
    ! Every wall type in the shipped case has identical layers - same
    ! thickness, same heat capacity, same conductivity all the way through - so
    ! a routine that indexes the wrong layer produces exactly the right answer
    ! and no check here can tell. That is not hypothetical: a mutation reading
    ! faclam(n,2) where the linearisation wants faclam(n,1) survived every
    ! check below until this was added.
    !
    ! The assertion afterwards is the part that matters. It is what makes the
    ! layer-indexed checks able to fail at all, so if the profile ever comes
    ! out flat the test says so rather than passing on a degenerate wall.
    if (myid == 0) then
      do n = 1, nfcts
        do j = 1, nfaclyrs
          facd(n,j)  = facd(n,j)  * (1. + 0.30 * real(j - 1))
          faccp(n,j) = faccp(n,j) * (1. + 0.20 * real(j - 1))
        end do
        do j = 1, nfaclyrs + 1
          faclam(n,j) = faclam(n,j) * (1. + 0.40 * real(j - 1))
        end do
      end do

      layered = .false.
      do n = 1, nfcts
        if (.not. live(n)) cycle
        if (faclam(n,1) /= faclam(n,2) .and. facd(n,1) /= facd(n,2) .and. &
            faccp(n,1) /= faccp(n,2)) layered = .true.
      end do
      if (.not. layered) then
        write(*, '(A)') 'FAIL: the wall layers are identical, so nothing here can'
        write(*, '(A)') '  detect a routine that reads the wrong one.'
        all_passed = .false.
      end if
    end if

    timee_s         = timee
    tnextEB_s       = tnextEB
    tEB_s           = tEB
    dt_s            = dt
    rk3step_s       = rk3step
    lvfsparse_s     = lvfsparse
    lwriteEBfiles_s = lwriteEBfiles
    lconstW_s       = lconstW

    ! ------------------------------------------------------------------
    ! 1. intqH: which stage is kept, what is cleared, and what is NOT summed
    ! ------------------------------------------------------------------
    ! intqH communicates nothing. Each rank integrates its own partial flux,
    ! and the sum over ranks happens once per energy balance instead of twice
    ! per step - see phase 2. So the probe is rank-dependent and every rank
    ! checks its OWN number: a routine that reduced here would leave the same
    ! total on all of them, which the two-rank pass catches on both.
    rank_sum = 0.5 * real(nprocs) * real(nprocs + 1)
    dt = dt_lcl
    fachfi = 0.
    facefi = 0.

    do stage = 1, 3
      fachf(1:nfcts) = probeH * real(myid + 1)
      facef(1:nfcts) = probeE * real(myid + 1)
      rk3step = stage
      call intqH

      ! Cleared on every stage, whichever stage it was. The device side relies
      ! on this to justify clearing fachf_d three times a step while
      ! integrating once - and letting it run on instead is what delivers the
      ! facet heat flux three times too large.
      if (maxval(abs(fachf(1:nfcts))) /= 0. .or. maxval(abs(facef(1:nfcts))) /= 0.) then
        call fail('intqH left a facet flux behind on stage', stage)
      end if

      ! Integrated on the third stage alone, and locally.
      if (stage < 3) then
        want = 0.
      else
        want = dt_lcl * probeH * real(myid + 1)
      end if
      got = maxval(abs(fachfi(1:nfcts)))
      if (abs(got - want) > 1.e-12 * max(1., want)) then
        write(*,'(A,I0,A,I0,A,ES23.15,A,ES23.15)') &
          'FAIL on rank ', myid, ': fachfi after stage ', stage, ' is ', got, &
          ', expected ', want
        all_passed = .false.
      end if
      want = 0.
      if (stage == 3) want = dt_lcl * probeE * real(myid + 1)
      got = maxval(abs(facefi(1:nfcts)))
      if (abs(got - want) > 1.e-12 * max(1., want)) then
        write(*,'(A,I0,A,I0,A,ES23.15,A,ES23.15)') &
          'FAIL on rank ', myid, ': facefi after stage ', stage, ' is ', got, &
          ', expected ', want
        all_passed = .false.
      end if
    end do

    ! Control: had intqH integrated on every stage it would now hold three
    ! times as much. Assert the two numbers are distinguishable, so a silently
    ! zero probe cannot make the check above pass for free.
    control_fired = abs(3. * dt_lcl * probeH - dt_lcl * probeH) > 1.e-9 * dt_lcl * probeH
    if (.not. control_fired) call fail('control: the stage-count probe cannot fire', 0)

    ! ------------------------------------------------------------------
    ! 2. EB sums the per-rank integrals, once, when it fires
    ! ------------------------------------------------------------------
    ! This is the collective that intqH no longer does. Give every rank a
    ! distinct local integral and check what comes out the other side:
    !
    !   - on rank 0, normalised to a flux, since that is where the balance is
    !     solved and where the division by tEB and faca happens;
    !   - elsewhere the raw reduced total, which also pins that the reduction
    !     is an all-reduce and that the normalisation is rank 0's alone;
    !   - and the local integrals are reset everywhere, or the next interval
    !     would double-count.
    call quiesce
    fachfi(1:nfcts) = probeH * real(myid + 1)
    facefi(1:nfcts) = probeE * real(myid + 1)
    call fire

    if (myid == 0) then
      worst = 0.
      do n = 1, nfcts
        if (.not. live(n)) cycle
        want = (probeH * rank_sum) / span / faca(n) * rhoa * cp
        worst = max(worst, abs(fachfsum(n) - want) / max(1., abs(want)))
      end do
      if (worst > 1.e-12) then
        write(*,'(A,ES23.15)') 'FAIL: the reduced facet heat flux is wrong, worst ', worst
        all_passed = .false.
      end if
    else
      want = probeH * rank_sum
      got  = maxval(abs(fachfsum(1:nfcts)))
      if (abs(got - want) > 1.e-12 * max(1., want)) then
        write(*,'(A,I0,A,ES23.15,A,ES23.15)') &
          'FAIL on rank ', myid, ': reduced facet heat flux is ', got, ', expected ', want
        all_passed = .false.
      end if
    end if

    if (maxval(abs(fachfi(1:nfcts))) /= 0. .or. maxval(abs(facefi(1:nfcts))) /= 0.) then
      call fail('EB did not reset the per-rank flux integrals', 0)
    end if

    ! Control: with one rank the reduced total and the local one coincide, so
    ! the reduction is untested. Say so rather than implying otherwise.
    if (myid == 0 .and. nprocs == 1) then
      write(*, '(A)') 'NOTE: on one rank the reduced and local facet fluxes coincide;'
      write(*, '(A)') '      the two-rank pass of run_test.sh is what separates them.'
    end if

    ! ------------------------------------------------------------------
    ! 3. When EB fires, and what it tells the GPU mirror
    ! ------------------------------------------------------------------
    call quiesce
    lfacetprops_dirty = .false.

    ! Not on the first two stages, however late it is.
    timee   = 1000.
    tnextEB = 100.
    tEB     = timee - span
    facT(1:nfcts,1) = T0
    do stage = 1, 2
      rk3step = stage
      call EB
      if (maxval(abs(facT(1:nfcts,1) - T0)) /= 0.) then
        call fail('EB moved the surface temperature on stage', stage)
      end if
      if (lfacetprops_dirty) then
        call fail('EB marked the facet properties dirty on stage', stage)
      end if
    end do

    ! Not on the third stage either, until the next EB time has arrived.
    rk3step = 3
    tnextEB = timee + 100.
    call EB
    if (maxval(abs(facT(1:nfcts,1) - T0)) /= 0.) then
      call fail('EB ran before its next time had arrived', 0)
    end if
    if (lfacetprops_dirty) then
      call fail('EB marked the facet properties dirty before it ran', 0)
    end if

    ! And then it does, and says so.
    tnextEB = timee - 1.
    tEB     = timee - span
    call EB
    if (.not. lfacetprops_dirty) then
      call fail('EB rewrote the facet properties without marking them dirty', 0)
    end if

    ! ------------------------------------------------------------------
    ! 4. Equilibrium is a fixed point
    ! ------------------------------------------------------------------
    ! With no sky and no view factors the incoming longwave is zero, so a
    ! facet radiating boltz*em*T^4 and absorbing exactly that in shortwave is
    ! in balance; the layers are uniform, so nothing conducts either.
    nlive = 0
    do n = 1, nfcts
      if (live(n)) nlive = nlive + 1
    end do
    if (nlive < 1) then
      if (myid == 0) write(*, '(A)') 'FAIL: no facet with a positive area and a real wall type'
      all_passed = .false.
    end if

    call quiesce
    worst = 0.
    if (myid == 0) then
      do n = 1, nfcts
        netsw(n) = boltz * facem(n) * T0**4
      end do
    end if
    call fire
    do n = 1, nfcts
      if (.not. live(n)) cycle
      worst = max(worst, maxval(abs(facT(n,1:nfaclyrs+1) - T0)))
    end do
    if (worst > 1.e-8) then
      if (myid == 0) write(*,'(A,ES23.15,A)') &
        'FAIL: a balanced facet drifted by ', worst, ' K, expected a fixed point'
      all_passed = .false.
    end if

    ! Control: the same setup with the surface flux 1% out must move, or the
    ! check above would pass on a routine that ignored its forcing entirely.
    call quiesce
    if (myid == 0) then
      do n = 1, nfcts
        netsw(n) = 1.01 * boltz * facem(n) * T0**4
      end do
    end if
    call fire
    control_fired = .false.
    do n = 1, nfcts
      if (.not. live(n)) cycle
      if (abs(facT(n,1) - T0) > 1.e-6) control_fired = .true.
    end do
    if (.not. control_fired) then
      call fail('control: a 1% flux imbalance left the facets unmoved', 0)
    end if

    ! ------------------------------------------------------------------
    ! 5. The surface energy balance closes, and the inside is held
    ! ------------------------------------------------------------------
    ! Heat the facet hard enough that the surface really moves, and drive all
    ! four terms of the net flux: shortwave, longwave, sensible and latent.
    call quiesce
    if (myid == 0) then
      do n = 1, nfcts
        netsw(n) = boltz * facem(n) * T0**4 + extraQ
      end do
    end if
    ! The whole contribution is put on rank 0 and the other ranks contribute
    ! nothing, so the reduced total is exactly H0 once normalised. Phase 2 is
    ! where the reduction itself is tested; here it only has to be predictable,
    ! and faca lives on rank 0 alone so no other rank could size its share.
    if (myid == 0) then
      do n = 1, nfcts
        if (faca(n) > 0.) fachfi(n) = H0 * span * faca(n) / (rhoa * cp)
        facefi(n) = E0
      end do
    end if
    facT_old = facT(1:nfcts,1:nfaclyrs+1)
    call fire

    worst = 0.
    do n = 1, nfcts
      if (.not. live(n)) cycle
      worst = max(worst, abs(facT(n,nfaclyrs+1) - T0))
    end do
    if (worst > 1.e-9) then
      if (myid == 0) write(*,'(A,ES23.15,A)') &
        'FAIL: the innermost layer moved by ', worst, ' K, it is meant to be held'
      all_passed = .false.
    end if

    control_fired = .false.
    do n = 1, nfcts
      if (.not. live(n)) cycle
      if (facT(n,1) - T0 > 1.e-3) control_fired = .true.
      if (facT(n,1) < T0 - 1.e-9) then
        call fail('a facet cooled under a net heat gain, facet', n)
      end if
    end do
    if (.not. control_fired) then
      call fail('control: net heating did not warm any surface', 0)
    end if

    ! facTdash is not broadcast, so the closure is rank 0's to check.
    if (myid == 0) then
      worst = 0.
      do n = 1, nfcts
        if (.not. live(n)) cycle
        resid = surface_residual(n, 1, 0.)
        scale = max(abs(netsw(n)), abs(faclam(n,1) * facTdash(n,1)), 1.)
        worst = max(worst, abs(resid) / scale)
      end do
      if (worst > 1.e-10) then
        write(*,'(A,ES23.15)') 'FAIL: the surface energy balance does not close, worst ', worst
        all_passed = .false.
      end if

      ! Control 1: the linearisation is around the OLD surface temperature.
      ! Using the new one has to break the closure, which also shows the
      ! surface moved far enough for the two to be told apart.
      worst = 0.
      do n = 1, nfcts
        if (.not. live(n)) cycle
        resid = surface_residual(n, 2, 0.)
        scale = max(abs(netsw(n)), abs(faclam(n,1) * facTdash(n,1)), 1.)
        worst = max(worst, abs(resid) / scale)
      end do
      if (worst < 1.e-8) then
        call fail('control: old and new surface temperatures are indistinguishable', 0)
      end if

      ! Control 2: dropping the sensible heat term must break it too, or the
      ! closure would be blind to how fachfi enters.
      worst = 0.
      do n = 1, nfcts
        if (.not. live(n)) cycle
        resid = surface_residual(n, 1, H0)
        scale = max(abs(netsw(n)), abs(faclam(n,1) * facTdash(n,1)), 1.)
        worst = max(worst, abs(resid) / scale)
      end do
      if (worst < 1.e-8) then
        call fail('control: the sensible heat term does not reach the closure', 0)
      end if

      ! --------------------------------------------------------------
      ! The conduction relation, layer by layer
      ! --------------------------------------------------------------
      ! Rows two and below of the system say the mean gradient across a layer
      ! is the temperature jump over its thickness:
      !
      !     (T'(j) + T'(j+1))/2 = (T(j+1) - T(j))/facd(n,j)
      !
      ! One equation per layer, naming facd(n,j) explicitly, needing neither
      ! the inverse nor any of C, D or E. Everything above this constrains the
      ! surface row and the held inner row and nothing in between, so a routine
      ! that assembled every layer from facd(n,1) passed the whole file. This
      ! is what catches that.
      worst = 0.
      moved = 0.
      do n = 1, nfcts
        if (.not. live(n)) cycle
        do j = 1, nfaclyrs
          if (facd(n,j) <= 0.) cycle
          resid = 0.5*(facTdash(n,j) + facTdash(n,j+1)) &
                  - (facT(n,j+1) - facT(n,j))/facd(n,j)
          scale = max(abs(facTdash(n,j)), abs(facTdash(n,j+1)), 1.)
          worst = max(worst, abs(resid)/scale)
        end do
        do j = 2, nfaclyrs
          moved = max(moved, abs(facT(n,j) - facT_old(n,j)))
        end do
      end do
      if (worst > 1.e-10) then
        write(*,'(A,ES23.15)') 'FAIL: the conduction relation does not hold, worst ', worst
        all_passed = .false.
      end if
      if (moved < 1.e-6) then
        call fail('control: no interior layer moved, so the relation above is trivial', 0)
      end if

      ! --------------------------------------------------------------
      ! Energy conservation over the step
      ! --------------------------------------------------------------
      ! Summing the layer equations telescopes the conductive terms and leaves
      !
      !   d(stored heat)/dt = flux in at the surface - flux out at the back
      !
      ! which in the scheme's own discrete form is
      !
      !   sum_j cp_j d_j (dT_j + dT_j+1)/2 + sum_j cp_j d_j^2 (dT'_j - dT'_j+1)/12
      !     = tEB ( -lam_1 T'_1 + lam_end T'_end )
      !
      ! The left side is the stored heat written as a sum over layers rather
      ! than as matrix rows, so it names faccp(n,j) and facd(n,j) for every j.
      ! This is the only check here that sees C and D at all.
      !
      ! T' at the old temperatures comes from the same two relations by forward
      ! substitution, with this step's forcing - not from a matrix inverse.
      worst  = 0.
      dshare = 0.
      do n = 1, nfcts
        if (.not. live(n)) cycle

        told(1:nfaclyrs+1) = facT_old(n,1:nfaclyrs+1)
        ab_n = boltz*facem(n)*told(1)**3/faclam(n,1)
        bb1  = -(netsw(n) + facLWin(n) + H0 + E0)/faclam(n,1)
        tdold(1) = bb1 + ab_n*told(1)
        do j = 1, nfaclyrs
          if (facd(n,j) <= 0.) then
            tdold(j+1) = -tdold(j)
          else
            tdold(j+1) = 2.*(told(j+1) - told(j))/facd(n,j) - tdold(j)
          end if
        end do

        lhs  = 0.
        dsum = 0.
        do j = 1, nfaclyrs
          lhs = lhs + faccp(n,j)*facd(n,j)/2. &
                      * ((facT(n,j) - told(j)) + (facT(n,j+1) - told(j+1)))
          dsum = dsum + faccp(n,j)*facd(n,j)**2/12. &
                      * ((facTdash(n,j) - tdold(j)) - (facTdash(n,j+1) - tdold(j+1)))
        end do

        rhs = 0.
        do j = 1, nfaclyrs
          rhs = rhs - faclam(n,j)*facTdash(n,j) + faclam(n,j+1)*facTdash(n,j+1)
        end do
        rhs = rhs * span

        scale = max(abs(lhs + dsum), abs(rhs), 1.)
        worst = max(worst, abs((lhs + dsum) - rhs)/scale)
        dshare = max(dshare, abs(dsum)/max(abs(lhs), 1.e-30))
      end do
      if (worst > 1.e-9) then
        write(*,'(A,ES23.15)') 'FAIL: the step does not conserve energy, worst ', worst
        all_passed = .false.
      end if
      ! The gradient correction has to carry enough weight for a mistake in it
      ! to show. If it is negligible the identity above is really only testing
      ! the heat capacity term.
      if (dshare < 1.e-3) then
        call fail('control: the gradient correction is negligible in the balance', 0)
      end if
    end if

    ! ------------------------------------------------------------------
    ! 6. calclw: the sparse and dense paths are the same calculation
    ! ------------------------------------------------------------------
    if (myid == 0 .and. lvfsparse .and. nnz > 0 .and. allocated(vfsparse)) then
      allocate(lw_sparse(1:nfcts), lw_dense(1:nfcts))
      ! quiesce blanked the radiative environment; put the case's own view
      ! factors and sky view back, or both paths would agree on zero.
      svf(1:nfcts)      = svf_s(1:nfcts)
      vfsparse(1:nnz)   = vfsparse_s(1:nnz)
      do n = 1, nfcts
        facT(n,1) = 280. + real(mod(n, 11))
      end do

      call calclw
      lw_sparse = facLWin(1:nfcts)

      allocate(vf(1:nfcts,1:nfcts)); vf = 0.
      do n = 1, nnz
        vf(ivfsparse(n), jvfsparse(n)) = vf(ivfsparse(n), jvfsparse(n)) + vfsparse(n)
      end do
      lvfsparse = .false.
      call calclw
      lw_dense = facLWin(1:nfcts)

      tol = 1.e-11 * max(1., maxval(abs(lw_sparse)))
      if (maxval(abs(lw_dense - lw_sparse)) > tol) then
        if (myid == 0) write(*,'(A,ES23.15)') &
          'FAIL: sparse and dense calclw disagree by ', maxval(abs(lw_dense - lw_sparse))
        all_passed = .false.
      end if

      ! Control: drop the largest view factor from the dense matrix. If that
      ! does not change the answer, the comparison above is between two arrays
      ! that never depended on the view factors at all.
      n = maxloc(abs(vfsparse(1:nnz)), 1)
      vf(ivfsparse(n), jvfsparse(n)) = 0.
      call calclw
      if (maxval(abs(facLWin(1:nfcts) - lw_sparse)) <= tol) then
        call fail('control: removing a view factor changed no incoming longwave', 0)
      end if

      lvfsparse = .true.
      deallocate(vf)
      deallocate(lw_sparse, lw_dense)
    else if (myid == 0) then
      write(*, '(A)') 'FAIL: this runmode needs sparse view factors.'
      write(*, '(A)') '  Set lvfsparse = .true. and nnz to the line count of'
      write(*, '(A)') '  vfsparse.inp.xxx, or the calclw cross-check between the'
      write(*, '(A)') '  sparse and dense paths cannot run at all.'
      all_passed = .false.
    end if

    ! ------------------------------------------------------------------
    ! Restore everything EB touched.
    ! ------------------------------------------------------------------
    facT     = facT_s
    fachf    = fachf_s
    facef    = facef_s
    faclGR   = faclGR_s
    fachfi   = fachfi_s
    facefi   = facefi_s
    fachfsum = 0.
    facefsum = 0.
    if (myid == 0) then
      facTdash = facTdash_s
      netsw    = netsw_s
      svf      = svf_s
      facLWin  = facLWin_s
      facwsoil = facwsoil_s
      faclam   = faclam_s
      faccp    = faccp_s
      facd     = facd_s
      if (allocated(vfsparse_s)) vfsparse = vfsparse_s
    end if

    timee             = timee_s
    tnextEB           = tnextEB_s
    tEB               = tEB_s
    dt                = dt_s
    rk3step           = rk3step_s
    lvfsparse         = lvfsparse_s
    lwriteEBfiles     = lwriteEBfiles_s
    lconstW           = lconstW_s
    ! Left dirty on purpose: this test moved the facet properties around, so a
    ! GPU build must refresh its mirror before the wall functions read them.
    lfacetprops_dirty = .true.

    deallocate(facT_s, fachf_s, facef_s, faclGR_s, facT_old, told, tdold)
    deallocate(fachfi_s, facefi_s)
    if (myid == 0) then
      deallocate(facTdash_s, netsw_s, svf_s, facLWin_s)
      deallocate(facwsoil_s, faclam_s, faccp_s, facd_s)
      if (allocated(vfsparse_s)) deallocate(vfsparse_s)
    end if

    if (myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      if (all_passed) then
        write(*, '(A)') 'ALL TESTS PASSED: tests_eb'
        write(*, '(A,I0,A,I0,A)') '  Stage contract, fixed point, surface closure and calclw over ', &
                                  nfcts, ' facet(s) on ', nprocs, ' rank(s)'
      else
        write(*, '(A)') 'TESTS FAILED: tests_eb'
        write(*, '(A)') '  One or more checks did not pass'
      end if
      write(*, '(A)') '================================================'
    end if

    tests_eb = all_passed

  contains

    !> A facet the energy balance actually solves for.
    !!
    !! facets and faca are on every rank; the wall properties are not, so the
    !! conductivity test can only be applied where faclam exists.
    logical function live(idx)
      integer, intent(in) :: idx
      live = (facets(idx) >= -100) .and. (faca(idx) > 0.)
      if (live .and. myid == 0) live = (faclam(idx,1) > 0.)
    end function live

    subroutine fail(msg, idx)
      character(len=*), intent(in) :: msg
      integer,          intent(in) :: idx
      if (idx == 0) then
        write(*,'(A,I0,A,A)') 'FAIL on rank ', myid, ': ', msg
      else
        write(*,'(A,I0,A,A,1X,I0)') 'FAIL on rank ', myid, ': ', msg, idx
      end if
      all_passed = .false.
    end subroutine fail

    !> Put the facets into a clean, quiet state: uniform temperature, no
    !! radiative environment, no green roofs, no accumulated flux.
    subroutine quiesce
      integer :: nn
      do nn = 1, nfcts
        facT(nn,1:nfaclyrs+1) = T0
      end do
      faclGR   = .false.
      lconstW  = .true.
      fachf    = 0.
      facef    = 0.
      fachfi   = 0.
      facefi   = 0.
      fachfsum = 0.
      facefsum = 0.
      if (myid == 0) then
        svf     = 0.
        netsw   = 0.
        facLWin = 0.
        if (allocated(vfsparse)) vfsparse = 0.
      end if
      dt      = dt_lcl
      timee   = 1000.
      rk3step = 3
    end subroutine quiesce

    !> Run one energy balance over `span` seconds of accumulated flux.
    subroutine fire
      tEB     = timee - span
      tnextEB = timee - 1.
      rk3step = 3
      call EB
    end subroutine fire

    !> Residual of the surface energy balance for facet n.
    !!
    !! Row one of the system is lam1*dT/dz = LWout - (SW + LWin + H + E), with
    !! the outgoing longwave linearised as boltz*em*Told^3*Tnew. Everything on
    !! the right is either set by this test or read back from the facet, so
    !! this is the physics and not a second copy of the matrix algebra.
    !!
    !! variant 2 uses the new surface temperature in the linearisation instead
    !! of the old one; drop_h subtracts a sensible heat flux that should be
    !! there. Both are negative controls.
    real function surface_residual(idx, variant, drop_h)
      integer, intent(in) :: idx, variant
      real,    intent(in) :: drop_h
      real :: lwout, influx, tcube

      if (variant == 2) then
        tcube = facT(idx,1)**3
      else
        tcube = facT_old(idx,1)**3
      end if
      lwout  = boltz * facem(idx) * tcube * facT(idx,1)
      influx = netsw(idx) + facLWin(idx) + (H0 - drop_h) + E0   ! fachfsum, facefsum
      surface_residual = faclam(idx,1) * facTdash(idx,1) - lwout + influx
    end function surface_residual

  end function tests_eb

  !> Validate the vegetation forcing against the equations it implements.
  !!
  !! Two things happened to vegetation_forcing when it was ported: every loop
  !! moved onto the device, and the time-independent part of the canopy energy
  !! balance was folded into a startup cache so that the per-step loop carries
  !! no exponentials and three divisions instead of six. Neither may move the
  !! physics, so what this test compares against is not the new code rearranged
  !! but the original expressions written out from the raw vegetation
  !! properties: the Beer-Lambert extinction with both its exponentials, the
  !! saturation-curve slope, the aerodynamic resistance in the form
  !! r_a = 130*sqrt(lsize/sqrt(wind2)), and the decoupling factor as
  !! omega = 1/(1 + 2*(gam/(s + 2*gam))*(rs/r_a)). If the cache or the folded
  !! divisions changed anything at all, it shows here.
  !!
  !! It also pins:
  !!   - the momentum drag on all three staggered face lists, and that nothing
  !!     outside those lists moves, which is what an index slip would do;
  !!   - qt = qtR + qtA, so the diagnostic split stays exhaustive;
  !!   - qe + qh = q_av_leaf, the absorbed radiation going somewhere;
  !!   - the mode dispatch, by running legacy and sveg over the same fields and
  !!     requiring each to match its own reference and to differ from the other;
  !!   - drag-only mode leaving heat and moisture untouched;
  !!   - scalar deposition, per scalar component.
  !!
  !! The shipped vegetation cases give every point the same properties, taken
  !! from one set of namelist values, so a per-point mix-up would reproduce the
  !! right answer everywhere. The setup below therefore gives lad, rs, lsize,
  !! ud, dec and sveg a bounded per-point profile and asserts it is not flat
  !! before trusting any of the comparisons.
  !!
  !! This exercises the host branch. The device kernels are covered by
  !! tests_cuda.f90::test_vegetation_forcing, which runs under
  !! UDALES_RUN_CUDA_SELFTEST on a Debug GPU build.
  !! Returns .true. if all checks pass, .false. otherwise.
  logical function tests_vegetation()
    use modglobal,  only : runmode, ib, ie, jb, je, kb, ke, ih, jh, kh, nsv, &
                           ltrees, itree_mode, TREE_MODE_DRAG_ONLY, &
                           TREE_MODE_SVEG, TREE_MODE_LEGACY_SEB, &
                           ltempeq, lmoist, Qstar, dzf, pref0, rlv, cp, rv, rd, rhoa
    use modfields,  only : initfields, um, vm, wm, up, vp, wp, &
                           thlm, qtm, thlp, qtp, svm, svp
    use vegetation, only : init_vegetation, init_vegetation_cache, vegetation_forcing, &
                           veg, vegp, sveg, vegetation_ready, &
                           npts_u, npts_v, npts_w, ijk_u, ijk_v, ijk_w, &
                           dcoef_u, dcoef_v, dcoef_w, veg_up, veg_vp, veg_wp
    use modmpi,     only : MPI_SUM

    implicit none

    real, parameter :: reltol = 1.e-11   !< reassociation only; the arithmetic is double

    logical :: all_passed
    integer :: m, n, i, j, k, mf, npts
    real    :: mine(1), everyone(1)
    real    :: want, got, spread
    real    :: ref_qt, ref_thl, ref_qtR, ref_qtA, ref_om, ref_q
    real    :: drag, uu, vv, ww
    real, allocatable :: legacy_qt(:), legacy_thl(:), legacy_om(:)
    logical, allocatable :: hit(:,:,:)

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_vegetation: CANOPY DRAG AND ENERGY BALANCE TEST'
      write(*, '(A)') '------------------------------------------------'
    end if

#if defined(_GPU)
    ! vegetation_forcing writes the device-resident tendencies in a GPU build,
    ! so the host arrays this test inspects would never change. Say so rather
    ! than reporting an unexplained mismatch, and do not pass vacuously.
    if (myid == 0) then
      write(*, '(A)') 'FAIL: runmode 1013 exercises the host branch of vegetation_forcing.'
      write(*, '(A)') '  Run it against a CPU build (build/cpu/<type>/u-dales).'
      write(*, '(A)') '  The device kernels are covered by tests_cuda.f90::test_vegetation_forcing,'
      write(*, '(A)') '  run with UDALES_RUN_CUDA_SELFTEST=1 on a Debug GPU build.'
    end if
    tests_vegetation = .false.
    return
#endif

    if (.not. ltrees) then
      if (myid == 0) then
        write(*, '(A)') 'FAIL: this runmode needs ltrees = .true. in the namelist.'
        write(*, '(A)') '  Without it init_vegetation returns immediately and every check is vacuous.'
      end if
      tests_vegetation = .false.
      return
    end if

    if (nsv < 1) then
      if (myid == 0) then
        write(*, '(A)') 'FAIL: this runmode needs nsv >= 1 in the namelist.'
        write(*, '(A)') '  The scalar deposition loop is one of the things under test.'
      end if
      tests_vegetation = .false.
      return
    end if

    all_passed = .true.

    ! The runmode dispatch happens before initfields and before
    ! init_vegetation, so this entry point does both itself. Heat and moisture
    ! are forced on: the canopy energy balance is skipped without them.
    ltempeq = .true.
    lmoist  = .true.
    call initfields

    itree_mode = TREE_MODE_LEGACY_SEB
    call init_vegetation

    ! init_vegetation broadcasts and barriers, and it is called once per mode
    ! below, so no rank may leave early. A rank whose subdomain holds no
    ! vegetation runs every collective and skips only the per-point checks; it
    ! is the global count that has to be non-zero.
    npts = veg%npts
    if (.not. vegetation_ready) npts = 0
    ! Counted through a real buffer: every other reduction in this file uses
    ! one, and MPI_ALLREDUCE has no explicit interface here.
    mine(1) = real(npts)
    call MPI_ALLREDUCE(mine, everyone, 1, MY_REAL, MPI_SUM, comm3d, mpierr)
    if (everyone(1) < 1.) then
      if (myid == 0) write(*, '(A)') 'FAIL: no vegetation points anywhere in the domain'
      tests_vegetation = .false.
      return
    end if

    allocate(legacy_qt(max(1,npts)), legacy_thl(max(1,npts)), legacy_om(max(1,npts)))
    allocate(hit(ib:ie, jb:je, kb:ke))

    call vary_properties
    call seed_fields
    call vegetation_forcing

    ! ----------------------------------------------------------------------
    ! 1. Momentum drag on the three staggered face lists.
    ! ----------------------------------------------------------------------
    mine(1) = real(npts_u + npts_v + npts_w)
    call MPI_ALLREDUCE(mine, everyone, 1, MY_REAL, MPI_SUM, comm3d, mpierr)
    if (everyone(1) < 1.) call vfail('control: no staggered vegetation faces in the domain', 0)

    do mf = 1, npts_u
      i = ijk_u(mf,1); j = ijk_u(mf,2); k = ijk_u(mf,3)
      uu = um(i,j,k)
      vv = 0.25*(vm(i,j,k) + vm(i,j+1,k) + vm(i-1,j,k) + vm(i-1,j+1,k))
      ww = 0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i-1,j,k) + wm(i-1,j,k+1))
      drag = -dcoef_u(mf)*uu*sqrt(uu*uu + vv*vv + ww*ww)
      if (.not. agrees(veg_up(mf), drag)) call vfail('u drag', mf)
      if (.not. agrees(up(i,j,k), drag)) call vfail('u drag scattered into up', mf)
    end do

    do mf = 1, npts_v
      i = ijk_v(mf,1); j = ijk_v(mf,2); k = ijk_v(mf,3)
      vv = vm(i,j,k)
      uu = 0.25*(um(i,j,k) + um(i+1,j,k) + um(i,j-1,k) + um(i+1,j-1,k))
      ww = 0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i,j-1,k) + wm(i,j-1,k+1))
      drag = -dcoef_v(mf)*vv*sqrt(vv*vv + uu*uu + ww*ww)
      if (.not. agrees(veg_vp(mf), drag)) call vfail('v drag', mf)
      if (.not. agrees(vp(i,j,k), drag)) call vfail('v drag scattered into vp', mf)
    end do

    do mf = 1, npts_w
      i = ijk_w(mf,1); j = ijk_w(mf,2); k = ijk_w(mf,3)
      ww = wm(i,j,k)
      uu = 0.25*(um(i,j,k) + um(i+1,j,k) + um(i,j,k-1) + um(i+1,j,k-1))
      vv = 0.25*(vm(i,j,k) + vm(i,j+1,k) + vm(i,j,k-1) + vm(i,j+1,k-1))
      drag = -dcoef_w(mf)*ww*sqrt(ww*ww + uu*uu + vv*vv)
      if (.not. agrees(veg_wp(mf), drag)) call vfail('w drag', mf)
      if (.not. agrees(wp(i,j,k), drag)) call vfail('w drag scattered into wp', mf)
    end do

    ! Nothing outside the face lists moved. Checked over the whole declared
    ! extent, halos included: a list built or indexed one cell off would still
    ! reproduce every value the loops above compare.
    call untouched_outside(up, npts_u, ijk_u, 'up')
    call untouched_outside(vp, npts_v, ijk_v, 'vp')
    call untouched_outside(wp, npts_w, ijk_w, 'wp')

    ! ----------------------------------------------------------------------
    ! 2. Canopy energy balance, legacy Beer-Lambert mode.
    ! ----------------------------------------------------------------------
    call check_canopy('legacy')

    do m = 1, npts
      legacy_qt(m)  = vegp%qt(m)
      legacy_thl(m) = vegp%thl(m)
      legacy_om(m)  = vegp%omega(m)
    end do

    ! ----------------------------------------------------------------------
    ! 3. Scalar deposition, one check per component.
    ! ----------------------------------------------------------------------
    do m = 1, npts
      i = veg%ijk(m,1); j = veg%ijk(m,2); k = veg%ijk(m,3)
      do n = 1, nsv
        want = -svm(i,j,k,n) * veg%lad(m) * veg%ud(m)
        if (.not. agrees(vegp%sv(m,n), want)) call vfail('scalar deposition', m)
        if (.not. agrees(svp(i,j,k,n), want)) call vfail('deposition scattered into svp', m)
      end do
    end do

    ! ----------------------------------------------------------------------
    ! 4. sveg mode: the radiation comes from sveg.inp, not from Qstar and the
    !    cumulative LAI. Both must match their own reference, and they must not
    !    agree with each other - otherwise the mode switch does nothing and
    !    check_canopy would pass either way.
    ! ----------------------------------------------------------------------
    itree_mode = TREE_MODE_SVEG
    call init_vegetation
    call vary_properties
    call seed_fields
    call vegetation_forcing
    call check_canopy('sveg')

    spread = 0.
    do m = 1, npts
      spread = max(spread, abs(vegp%qt(m) - legacy_qt(m)))
      spread = max(spread, abs(vegp%thl(m) - legacy_thl(m)))
    end do
    if (npts > 0 .and. spread <= 0.) then
      call vfail('control: legacy and sveg radiation give identical tendencies', 0)
    end if

    ! ----------------------------------------------------------------------
    ! 5. Drag-only mode: the canopy loop must not run at all, but the momentum
    !    drag must still be applied.
    ! ----------------------------------------------------------------------
    itree_mode = TREE_MODE_DRAG_ONLY
    call init_vegetation
    call vary_properties
    call seed_fields
    call vegetation_forcing

    do m = 1, npts
      if (vegp%qt(m) /= 0.)  call vfail('drag-only mode moved qt', m)
      if (vegp%thl(m) /= 0.) call vfail('drag-only mode moved thl', m)
    end do
    if (any(thlp /= 0.)) call vfail('drag-only mode wrote thlp', 0)
    if (any(qtp  /= 0.)) call vfail('drag-only mode wrote qtp', 0)

    spread = 0.
    do mf = 1, npts_u
      spread = max(spread, abs(veg_up(mf)))
    end do
    if (npts_u > 0 .and. spread <= 0.) call vfail('control: drag-only mode applied no drag either', 0)

    deallocate(legacy_qt, legacy_thl, legacy_om, hit)

    if (myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      if (all_passed) then
        write(*, '(A)') 'tests_vegetation: PASSED'
      else
        write(*, '(A)') 'tests_vegetation: FAILED'
      end if
      write(*, '(A)') '================================================'
    end if

    tests_vegetation = all_passed

  contains

    subroutine vfail(msg, idx)
      character(len=*), intent(in) :: msg
      integer,          intent(in) :: idx
      if (idx == 0) then
        write(*,'(A,I0,A,A)') 'FAIL on rank ', myid, ': ', msg
      else
        write(*,'(A,I0,A,A,1X,I0)') 'FAIL on rank ', myid, ': ', msg, idx
      end if
      all_passed = .false.
    end subroutine vfail

    !> Relative comparison, absolute near zero.
    logical function agrees(a, b)
      real, intent(in) :: a, b
      agrees = abs(a - b) <= reltol*max(1.e-30, abs(b))
    end function agrees

    !> Give the vegetation points properties that actually differ.
    !!
    !! Bounded and non-monotone, so no point ends up with a degenerate lad or a
    !! resistance large enough to drive the whole canopy term to zero. The face
    !! drag coefficients are deliberately left alone: they were built from the
    !! original lad in init_vegetation and are read straight out of dcoef_u/v/w
    !! by both the routine and the checks above, so the two paths stay
    !! independent.
    subroutine vary_properties
      integer :: mm
      real    :: lo, hi

      Qstar = 550.
      do mm = 1, veg%npts
        veg%lad(mm)   = 0.60 + 0.40*sin(0.90*real(mm))
        veg%rs(mm)    = 60.0 + 30.0*sin(0.53*real(mm) + 1.)
        veg%lsize(mm) = 0.15 + 0.08*sin(0.31*real(mm) + 2.)
        veg%ud(mm)    = 2.0e-4*(1.5 + sin(0.77*real(mm) + 3.))
        veg%dec(mm)   = 0.35 + 0.15*sin(0.19*real(mm) + 4.)
        sveg(mm)      = 250. + 120.*sin(0.61*real(mm) + 5.)
      end do
      call init_vegetation_cache

      ! The comparisons below are only meaningful if the points differ. A case
      ! whose vegetation is uniform would let a per-point index slip reproduce
      ! every value, so refuse to report a pass on one.
      if (veg%npts > 1) then
        lo = minval(veg%lad(1:veg%npts)); hi = maxval(veg%lad(1:veg%npts))
        if (hi - lo <= 0.) call vfail('control: leaf area density is flat across points', 0)
        lo = minval(veg%rs(1:veg%npts)); hi = maxval(veg%rs(1:veg%npts))
        if (hi - lo <= 0.) call vfail('control: stomatal resistance is flat across points', 0)
        lo = minval(sveg(1:veg%npts)); hi = maxval(sveg(1:veg%npts))
        if (hi - lo <= 0.) call vfail('control: absorbed shortwave is flat across points', 0)
      end if
      if (Qstar <= 0.) call vfail('control: Qstar is zero, so legacy radiation vanishes', 0)
    end subroutine vary_properties

    !> Spatially varying previous-step fields, and cleared tendencies.
    subroutine seed_fields
      integer :: nn
      call seed3(um,   0.80, 0.60)
      call seed3(vm,  -0.30, 0.50)
      call seed3(wm,   0.20, 0.40)
      call seed3(thlm, 292.0, 3.0)
      call seed3(qtm,  0.009, 0.002)
      do nn = 1, nsv
        call seed3(svm(:,:,:,nn), 1.0 + 0.5*real(nn), 0.4)
      end do
      up = 0.; vp = 0.; wp = 0.; thlp = 0.; qtp = 0.; svp = 0.
    end subroutine seed_fields

    !> Fill a field with a bounded, index-dependent pattern.
    !!
    !! Assumed-shape, so the index origin is lost - which does not matter: the
    !! checks read the same array back through the same indices, and all this
    !! has to do is make neighbouring cells differ.
    subroutine seed3(a, base, amp)
      real, intent(inout) :: a(:,:,:)
      real, intent(in)    :: base, amp
      integer :: i1, i2, i3
      do i3 = 1, size(a,3)
        do i2 = 1, size(a,2)
          do i1 = 1, size(a,1)
            a(i1,i2,i3) = base + amp*sin(0.70*real(i1) + 0.30*real(i2) + 0.11*real(i3))
          end do
        end do
      end do
    end subroutine seed3

    !> Every cell outside the point list must still hold exactly zero.
    subroutine untouched_outside(fld, nlist, list, label)
      real,             intent(in) :: fld(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh)
      integer,          intent(in) :: nlist
      integer,          intent(in) :: list(:,:)
      character(len=*), intent(in) :: label
      integer :: ii, jj, kk, mm

      hit = .false.
      do mm = 1, nlist
        hit(list(mm,1), list(mm,2), list(mm,3)) = .true.
      end do

      do kk = kb, ke+kh
        do jj = jb-jh, je+jh
          do ii = ib-ih, ie+ih
            if (ii >= ib .and. ii <= ie .and. jj >= jb .and. jj <= je .and. &
                kk >= kb .and. kk <= ke) then
              if (hit(ii,jj,kk)) cycle
            end if
            if (fld(ii,jj,kk) /= 0.) then
              call vfail(label//' moved outside the vegetation face list at k =', kk)
              return
            end if
          end do
        end do
      end do
    end subroutine untouched_outside

    !> Compare the canopy tendencies against the original expressions.
    subroutine check_canopy(label)
      character(len=*), intent(in) :: label
      integer :: mm

      do mm = 1, veg%npts
        call canopy_reference(mm, ref_qt, ref_thl, ref_qtR, ref_qtA, ref_om, ref_q)

        if (.not. agrees(vegp%omega(mm), ref_om)) call vfail(label//': decoupling factor', mm)
        if (.not. agrees(vegp%qt(mm),    ref_qt))  call vfail(label//': moisture tendency', mm)
        if (.not. agrees(vegp%thl(mm),   ref_thl)) call vfail(label//': heat tendency', mm)
        if (.not. agrees(vegp%qtR(mm),   ref_qtR)) call vfail(label//': radiative moisture share', mm)
        if (.not. agrees(vegp%qtA(mm),   ref_qtA)) call vfail(label//': aerodynamic moisture share', mm)

        ! The split is exhaustive.
        if (.not. agrees(vegp%qt(mm), vegp%qtR(mm) + vegp%qtA(mm))) then
          call vfail(label//': qt is not qtR + qtA', mm)
        end if

        ! And the absorbed radiation all goes somewhere: qe + qh = q_av_leaf,
        ! recovered from the two tendencies through their own scalings.
        want = ref_q
        got  = vegp%qt(mm)*(rhoa*rlv)/veg%lad(mm) + vegp%thl(mm)*(rhoa*cp)/veg%lad(mm)
        if (.not. agrees(got, want)) call vfail(label//': canopy energy does not close', mm)
      end do

      ! A canopy that absorbs nothing would satisfy every line above trivially.
      spread = 0.
      do mm = 1, veg%npts
        spread = max(spread, abs(vegp%thl(mm)))
      end do
      if (veg%npts > 0 .and. spread <= 0.) then
        call vfail('control: '//label//' produced no heat tendency at all', 0)
      end if
    end subroutine check_canopy

    !> The canopy energy balance as it was written before the cache.
    !!
    !! Deliberately verbose and division-heavy: this is the reference the
    !! folded form has to reproduce, so it keeps slope_sat, r_a and gam
    !! explicit rather than sharing anything with vegetation.f90.
    subroutine canopy_reference(mm, o_qt, o_thl, o_qtR, o_qtA, o_om, o_q)
      integer, intent(in)  :: mm
      real,    intent(out) :: o_qt, o_thl, o_qtR, o_qtA, o_om, o_q
      integer :: ii, jj, kk
      real :: ladv, decv, clai, rn_top, rn_bot, q_av_leaf
      real :: e_sat, e_vap, d_vap, slope_sat, r_a, omega, qe, qh, gam
      real :: lsizev, rsv, wind2

      ii = veg%ijk(mm,1); jj = veg%ijk(mm,2); kk = veg%ijk(mm,3)

      ladv = veg%lad(mm)
      gam  = (cp*pref0*rv)/(rlv*rd)

      if (itree_mode == TREE_MODE_SVEG) then
        q_av_leaf = sveg(mm) / max(ladv, 1.0e-12)
      else
        decv = veg%dec(mm)
        clai = veg%laiv(mm)
        rn_top = Qstar * exp(-decv * (clai - ladv * dzf(kk)))
        rn_bot = Qstar * exp(-decv * clai)
        q_av_leaf = (rn_top - rn_bot) / (dzf(kk) * max(ladv, 1.0e-12))
      end if

      lsizev = max(veg%lsize(mm), 1.0e-6)
      rsv    = max(veg%rs(mm), 1.0e-6)

      e_sat = 610.8*exp((17.27*(thlm(ii,jj,kk)-273.15))/(thlm(ii,jj,kk)-35.85))
      e_vap = (qtm(ii,jj,kk) * pref0) / (0.378 * qtm(ii,jj,kk) + 0.622)
      d_vap = max(e_sat - e_vap, 0.)
      slope_sat = (4098*e_sat)/((thlm(ii,jj,kk)-35.85)**2)

      wind2 = max((0.5*(um(ii,jj,kk)+um(ii+1,jj,kk)))**2 &
                +(0.5*(vm(ii,jj,kk)+vm(ii,jj+1,kk)))**2 &
                +(0.5*(wm(ii,jj,kk)+wm(ii,jj,kk+1)))**2, 1.0e-12)
      r_a = 130*sqrt(lsizev / sqrt(wind2))

      omega = 1/(1 + 2*(gam/(slope_sat+2*gam)) * (rsv/r_a))
      qe = omega*(slope_sat/(slope_sat+2*gam))*q_av_leaf + (1-omega)*(1/(gam*rsv))*rhoa*cp*d_vap
      qh = q_av_leaf - qe

      o_om  = omega
      o_qtR = ladv*(omega*(slope_sat/(slope_sat+2*gam))*q_av_leaf)/(rhoa*rlv)
      o_qtA = ladv*((1-omega)*(1/(gam*rsv))*rhoa*cp*d_vap)/(rhoa*rlv)
      o_qt  = ladv*qe/(rhoa*rlv)
      o_thl = ladv*qh/(rhoa*cp)
      o_q   = q_av_leaf
    end subroutine canopy_reference

  end function tests_vegetation

  !> Check the host branch of the checksim diagnostics.
  !!
  !! checksim itself only prints, so what is tested here are the three
  !! rank-local reductions it prints - courant_local, diffnr_local and
  !! div_local - plus the diffusion-number geometry cache they were rewritten
  !! around. Each is driven with a state whose answer follows from the grid
  !! rather than from re-running the loop, so a transcription error in the port
  !! cannot be reproduced by the expectation.
  !!
  !! What each group of checks is for:
  !!
  !!   - diffnrgeom is asserted equal, bit for bit, to the expression it
  !!     replaced. The point of caching it was that the printed diffusion
  !!     number must not move, and only exact equality says that.
  !!
  !!   - A single spike in one field, with everything else zero, isolates one
  !!     term of a sum: the maximum is then that term alone and can be written
  !!     down. Placed at each bound of the reduction box in turn, it also pins
  !!     the loop extent, and a matching spike placed just outside the box
  !!     asserts the loop stops where it should. A port that reduces over the
  !!     halo passes an interior-only comparison.
  !!
  !!   - div_local reads i+1, j+1 and k+1, so its spikes are planted at ie+1,
  !!     je+1 and ke+1: the cell that sees them is the last interior one. A
  !!     port written u0(i) - u0(i-1), or one that stops at ie-1, reads zero
  !!     there and returns nothing.
  !!
  !!   - The ramp gives every cell the same divergence, so the volume integral
  !!     is a cell count times a known value. That is the check that fails when
  !!     the sum covers one plane too few, which no single-spike test can see.
  !!
  !! This exercises the host branch. The device kernels are covered by
  !! tests_cuda.f90::test_checksim_reductions, which runs under
  !! UDALES_RUN_CUDA_SELFTEST on a Debug GPU build.
  !! Returns .true. if all checks pass, .false. otherwise.
  logical function tests_checksim()
    use modglobal,      only : runmode, ib, ie, jb, je, kb, ke, &
                               dx, dy, dxi, dyi, dy2i, dxhi, dzhi, dzf, dzfi, dzh, dxh2i, dvcell
    use modfields,      only : initfields, um, vm, wm, u0, v0, w0
    use modsubgrid,     only : initsubgrid
    use modsubgriddata, only : ekm, ekh
    use modchecksim,    only : initchecksim, courant_local, diffnr_local, div_local, diffnrgeom

    implicit none

    real, parameter :: spike = 2., alpha = 0.5, dtm = 0.25

    logical :: all_passed
    integer :: i, j, k, ni, nj
    real    :: got, want, gotmax, gottot, dzfsum
    real    :: dxi_save, dyi_save, dx_save, dy_save
    real, allocatable :: dxhi_save(:), dzhi_save(:), dzf_save(:), dzfi_save(:), dv_save(:), geom_save(:,:)

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_checksim: COURANT, DIFFUSION AND DIVERGENCE TEST'
      write(*, '(A)') '------------------------------------------------'
    end if

#if defined(_GPU)
    ! Under _GPU the three reductions read the device fields, so seeding the
    ! host arrays this test writes would leave them looking at whatever the
    ! device happens to hold. Say so rather than reporting an unexplained
    ! mismatch, and do not pass vacuously.
    if (myid == 0) then
      write(*, '(A)') 'FAIL: runmode 1014 exercises the host branch of the checksim reductions.'
      write(*, '(A)') '  Run it against a CPU build (build/cpu/<type>/u-dales).'
      write(*, '(A)') '  The device kernels are covered by tests_cuda.f90::test_checksim_reductions,'
      write(*, '(A)') '  run with UDALES_RUN_CUDA_SELFTEST=1 on a Debug GPU build.'
    end if
    tests_checksim = .false.
    return
#endif

    call initfields
    call initsubgrid
    call initchecksim

    all_passed = .true.
    ni = ie - ib + 1
    nj = je - jb + 1

    ! ---- the diffusion-number geometry cache --------------------------------
    ! Exact equality, not a tolerance: the whole justification for precomputing
    ! this was that it is the same arithmetic in the same order.
    do k = kb, ke
      do i = ib, ie
        want = 1/dzh(k)**2 + dxh2i(i) + dy2i
        if (diffnrgeom(i,k) /= want) then
          write(*,'(A,I0,A,I0,A,I0,A,ES22.14,A,ES22.14)') 'FAIL on rank ', myid, &
               ': diffnrgeom(', i, ',', k, ') = ', diffnrgeom(i,k), ' want ', want
          all_passed = .false.
        end if
      end do
    end do

    ! Same, for the cell volume modglobal builds. This one is not claimed to be
    ! bit-identical to the div*dx*dy*dzf(k) it replaces - see div_local - but it
    ! must be the product of exactly those three, so a cache built from dzh, or
    ! from dzfi, is caught.
    do k = kb, ke
      want = dx*dy*dzf(k)
      if (dvcell(k) /= want) then
        write(*,'(A,I0,A,I0,A,ES22.14,A,ES22.14)') 'FAIL on rank ', myid, &
             ': dvcell(', k, ') = ', dvcell(k), ' want ', want
        all_passed = .false.
      end if
    end do

    ! ---- make every grid factor distinct ------------------------------------
    ! Every case in the repository has dx = dy and a uniform grid in x, and all
    ! but one has a uniform grid in z. On such a grid dxhi(i) equals dyi,
    ! dxhi(i) equals dxhi(ib), and dzhi, dzfi and dzf are constant in k - so a
    ! loop that pairs the v term with the x spacing, or reads any of these at a
    ! fixed index, returns exactly the right answer and no assertion below can
    ! fail. Substituting values that differ in every position is what makes the
    ! pairing and the indexing load-bearing.
    !
    ! dzfi is deliberately not the reciprocal of dzf, and dxi not that of dx:
    ! div_local uses one of each pair for the gradient and the other for the
    ! cell volume, and on a consistent grid a kernel that confuses them is
    ! invisible.
    !
    ! This runs after the exactness check above, which has to see the grid
    ! initchecksim actually built, and everything is restored before returning.
    allocate(dxhi_save(lbound(dxhi,1):ubound(dxhi,1))) ; dxhi_save = dxhi
    allocate(dzhi_save(lbound(dzhi,1):ubound(dzhi,1))) ; dzhi_save = dzhi
    allocate(dzf_save (lbound(dzf ,1):ubound(dzf ,1))) ; dzf_save  = dzf
    allocate(dzfi_save(lbound(dzfi,1):ubound(dzfi,1))) ; dzfi_save = dzfi
    allocate(dv_save  (lbound(dvcell,1):ubound(dvcell,1))) ; dv_save   = dvcell
    allocate(geom_save(lbound(diffnrgeom,1):ubound(diffnrgeom,1), &
                       lbound(diffnrgeom,2):ubound(diffnrgeom,2)))
    geom_save = diffnrgeom
    dxi_save = dxi ; dyi_save = dyi ; dx_save = dx ; dy_save = dy

    ! Strictly increasing, and exact in binary: no two indices can collide, so
    ! reading any of these one index off is always visible. A cyclic pattern
    ! would not be - with a period of seven and a 64-wide domain, dxhi(ib) and
    ! dxhi(ie) come out equal and a kernel pinned to ib passes.
    do i = lbound(dxhi,1), ubound(dxhi,1)
      dxhi(i) = 0.375 + 0.0078125*real(i - lbound(dxhi,1))
    end do
    do k = lbound(dzhi,1), ubound(dzhi,1)
      dzhi(k) = 0.75 + 0.0078125*real(k - lbound(dzhi,1))
    end do
    do k = lbound(dzf,1), ubound(dzf,1)
      dzf(k) = 1.5 + 0.0078125*real(k - lbound(dzf,1))
    end do
    do k = lbound(dzfi,1), ubound(dzfi,1)
      dzfi(k) = 0.25 + 0.00390625*real(k - lbound(dzfi,1))
    end do
    do k = lbound(dvcell,1), ubound(dvcell,1)
      dvcell(k) = 2.5 + 0.0078125*real(k - lbound(dvcell,1))
    end do
    do k = lbound(diffnrgeom,2), ubound(diffnrgeom,2)
      do i = lbound(diffnrgeom,1), ubound(diffnrgeom,1)
        diffnrgeom(i,k) = 1.25 + 0.0078125*real( &
          (k - lbound(diffnrgeom,2))*(ubound(diffnrgeom,1) - lbound(diffnrgeom,1) + 1) &
          + (i - lbound(diffnrgeom,1)))
      end do
    end do
    dxi = 0.25 ; dyi = 0.5 ; dx = 3. ; dy = 5.

    ! ---- the Courant number -------------------------------------------------
    um = 0. ; vm = 0. ; wm = 0.
    call courant_local(dtm, got)
    if (.not. same('courant of a field at rest', got, 0.)) all_passed = .false.

    ! One term at a time, so each grid factor is pinned separately. A port that
    ! pairs vm with dzhi, say, reproduces neither of these.
    if (.not. courant_spike('courant um', 1, ie, je, ke, dxhi(ie))) all_passed = .false.
    if (.not. courant_spike('courant vm', 2, ie, je, ke, dyi     )) all_passed = .false.
    if (.not. courant_spike('courant wm', 3, ie, je, ke, dzhi(ke))) all_passed = .false.

    ! The same spike at every bound of the reduction box, so no face of it can
    ! be missing.
    if (.not. courant_spike('courant corner ib,jb,kb', 1, ib, jb, kb, dxhi(ib))) all_passed = .false.
    if (.not. courant_spike('courant corner ie,jb,kb', 1, ie, jb, kb, dxhi(ie))) all_passed = .false.
    if (.not. courant_spike('courant corner ib,je,kb', 1, ib, je, kb, dxhi(ib))) all_passed = .false.
    if (.not. courant_spike('courant corner ib,jb,ke', 1, ib, jb, ke, dxhi(ib))) all_passed = .false.

    ! ... and just outside every one of those bounds, where it must be ignored.
    if (.not. courant_poison('courant below ib', ib-1, jb,   kb  )) all_passed = .false.
    if (.not. courant_poison('courant above ie', ie+1, jb,   kb  )) all_passed = .false.
    if (.not. courant_poison('courant below jb', ib,   jb-1, kb  )) all_passed = .false.
    if (.not. courant_poison('courant above je', ib,   je+1, kb  )) all_passed = .false.
    if (.not. courant_poison('courant below kb', ib,   jb,   kb-1)) all_passed = .false.
    if (.not. courant_poison('courant above ke', ib,   jb,   ke+1)) all_passed = .false.

    ! The time step has to multiply the whole expression, not one term of it.
    um = 0. ; vm = 0. ; wm = 0.
    um(ie,je,ke) = spike
    call courant_local(2.*dtm, got)
    if (.not. same('courant scales with dt', got, spike*dxhi(ie)*(2.*dtm))) all_passed = .false.

    ! ---- the diffusion number -----------------------------------------------
    ekm = 0. ; ekh = 0.
    call diffnr_local(dtm, got)
    if (.not. same('diffusion number of zero viscosity', got, 0.)) all_passed = .false.

    ! ekm and ekh are checked one at a time with the other zeroed, so dropping
    ! either from the max is caught. A run whose Prandtl number is below one is
    ! limited by ekh, and that is the term a port is most likely to lose.
    if (.not. diffnr_spike('diffusion ekm', .true.,  ie, je, ke)) all_passed = .false.
    if (.not. diffnr_spike('diffusion ekh', .false., ie, je, ke)) all_passed = .false.
    if (.not. diffnr_spike('diffusion ekm at ib,jb,kb', .true.,  ib, jb, kb)) all_passed = .false.
    if (.not. diffnr_spike('diffusion ekh at ib,jb,kb', .false., ib, jb, kb)) all_passed = .false.

    if (.not. diffnr_poison('diffusion below ib', ib-1, jb,   kb  )) all_passed = .false.
    if (.not. diffnr_poison('diffusion above ie', ie+1, jb,   kb  )) all_passed = .false.
    if (.not. diffnr_poison('diffusion below jb', ib,   jb-1, kb  )) all_passed = .false.
    if (.not. diffnr_poison('diffusion above je', ib,   je+1, kb  )) all_passed = .false.
    if (.not. diffnr_poison('diffusion below kb', ib,   jb,   kb-1)) all_passed = .false.
    if (.not. diffnr_poison('diffusion above ke', ib,   jb,   ke+1)) all_passed = .false.

    ekm = 0. ; ekh = 0.
    ekm(ie,je,ke) = spike
    call diffnr_local(2.*dtm, got)
    if (.not. same('diffusion number scales with dt', got, spike*diffnrgeom(ie,ke)*(2.*dtm))) all_passed = .false.

    ! ---- the divergence -----------------------------------------------------
    u0 = 0. ; v0 = 0. ; w0 = 0.
    call div_local(gotmax, gottot)
    if (.not. same('divmax of a divergence-free field', gotmax, 0.)) all_passed = .false.
    if (.not. same('divtot of a divergence-free field', gottot, 0.)) all_passed = .false.

    ! The forward neighbour, planted one cell past the top of the box. The cell
    ! that sees it is the last interior one, so a backward difference or a loop
    ! that stops one short both return zero here.
    u0 = 0. ; v0 = 0. ; w0 = 0.
    u0(ie+1,je,ke) = spike
    call div_local(gotmax, gottot)
    if (.not. same('divmax from u0(ie+1)', gotmax, spike*dxi)) all_passed = .false.
    if (.not. same('divtot from u0(ie+1)', gottot, (spike*dxi)*dvcell(ke))) all_passed = .false.

    u0 = 0. ; v0 = 0. ; w0 = 0.
    v0(ie,je+1,ke) = spike
    call div_local(gotmax, gottot)
    if (.not. same('divmax from v0(je+1)', gotmax, spike*dyi)) all_passed = .false.
    if (.not. same('divtot from v0(je+1)', gottot, (spike*dyi)*dvcell(ke))) all_passed = .false.

    u0 = 0. ; v0 = 0. ; w0 = 0.
    w0(ie,je,ke+1) = spike
    call div_local(gotmax, gottot)
    if (.not. same('divmax from w0(ke+1)', gotmax, spike*dzfi(ke))) all_passed = .false.
    if (.not. same('divtot from w0(ke+1)', gottot, (spike*dzfi(ke))*dvcell(ke))) all_passed = .false.

    ! The sign is not free: a divergence built the other way round must come
    ! back negative in the integral, even though the maximum is an absolute one.
    ! The spike goes at ib, not ie: at ie the cell below it sees an equal and
    ! opposite divergence and the two cancel in the integral, whereas at ib
    ! that partner cell is outside the box. Which also pins the lower i bound.
    u0 = 0. ; v0 = 0. ; w0 = 0.
    u0(ib,je,ke) = spike
    call div_local(gotmax, gottot)
    if (.not. same('divmax from u0(ib)',  gotmax,  spike*dxi)) all_passed = .false.
    if (.not. same('divtot from u0(ib)',  gottot, -(spike*dxi)*dvcell(ke))) all_passed = .false.

    ! A ramp in x gives every cell the same divergence, so the volume integral
    ! is a cell count times a known value: the check that the sum covers the
    ! whole box and not one plane fewer.
    do k = lbound(u0,3), ubound(u0,3)
      do j = lbound(u0,2), ubound(u0,2)
        do i = lbound(u0,1), ubound(u0,1)
          u0(i,j,k) = alpha*(i - ib)
        end do
      end do
    end do
    v0 = 0. ; w0 = 0.
    dzfsum = 0.
    do k = kb, ke
      dzfsum = dzfsum + dvcell(k)
    end do
    call div_local(gotmax, gottot)
    ! Every cell holds the same divergence exactly - alpha and the index
    ! differences are both exact in binary - so the maximum is exact too. The
    ! integral is not: it is a running sum of ni*nj*nk terms, and the loop adds
    ! them in an order the closed form below does not reproduce. Hence a
    ! tolerance that scales with the term count. It is still four orders of
    ! magnitude tighter than the 1/nk error a missing plane would cause.
    if (.not. same('divmax of a uniform ramp', gotmax, alpha*dxi)) all_passed = .false.
    if (.not. nearly('divtot of a uniform ramp', gottot, &
                    (alpha*dxi)*dzfsum*real(ni)*real(nj), &
                    real(ni)*real(nj)*real(ke-kb+1))) all_passed = .false.

    dxhi = dxhi_save ; dzhi = dzhi_save ; dzf = dzf_save ; dzfi = dzfi_save
    diffnrgeom = geom_save ; dvcell = dv_save
    dxi = dxi_save ; dyi = dyi_save ; dx = dx_save ; dy = dy_save
    deallocate(geom_save, dv_save, dzfi_save, dzf_save, dzhi_save, dxhi_save)

    if (myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      if (all_passed) then
        write(*, '(A)')      'ALL TESTS PASSED: tests_checksim'
        write(*, '(A,I0,A,I0,A,I0,A)') '  Reduction box ', ni, ' x ', nj, ' x ', ke-kb+1, ' cells per rank'
      else
        write(*, '(A)') 'TESTS FAILED: tests_checksim'
        write(*, '(A)') '  One or more checks did not pass'
      end if
      write(*, '(A)') '================================================'
    end if

    tests_checksim = all_passed

  contains

    !> Assert a scalar matches to within a few ulps of its own magnitude.
    logical function same(label, got, want)
      character(len=*), intent(in) :: label
      real,             intent(in) :: got, want
      real :: tol

      tol  = 64. * epsilon(1.) * max(1., abs(want))
      same = abs(got - want) <= tol
      if (.not. same) then
        write(*,'(A,I0,3A,ES22.14,A,ES22.14)') 'FAIL on rank ', myid, ': checksim ', &
             trim(label), ' got ', got, ' want ', want
      end if

    end function same

    !> Assert a scalar matches within the rounding a sum of terms can accrue.
    !!
    !! Sequential summation of n terms carries a worst-case relative error of
    !! order n*epsilon, so the tolerance is scaled by the term count rather
    !! than fixed: a bound that holds on a 64^3 box and on a 256^3 one.
    logical function nearly(label, got, want, terms)
      character(len=*), intent(in) :: label
      real,             intent(in) :: got, want, terms
      real :: tol

      tol    = 8. * terms * epsilon(1.) * max(1., abs(want))
      nearly = abs(got - want) <= tol
      if (.not. nearly) then
        write(*,'(A,I0,3A,ES22.14,A,ES22.14)') 'FAIL on rank ', myid, ': checksim ', &
             trim(label), ' got ', got, ' want ', want
      end if

    end function nearly

    !> Put a single positive value in one momentum field - which selects um,
    !! vm or wm - and assert the maximum Courant number is that value times
    !! want_factor, the grid factor that term is supposed to carry.
    logical function courant_spike(label, which, is, js, ks, want_factor)
      character(len=*), intent(in) :: label
      integer,          intent(in) :: which, is, js, ks
      real,             intent(in) :: want_factor
      real :: got

      um = 0. ; vm = 0. ; wm = 0.
      select case (which)
      case (1) ; um(is,js,ks) = spike
      case (2) ; vm(is,js,ks) = spike
      case (3) ; wm(is,js,ks) = spike
      end select

      call courant_local(dtm, got)
      courant_spike = same(label, got, spike*want_factor*dtm)

    end function courant_spike

    !> Put a large value outside the reduction box and assert it is not seen.
    logical function courant_poison(label, is, js, ks)
      character(len=*), intent(in) :: label
      integer,          intent(in) :: is, js, ks
      real :: got

      um = 0. ; vm = 0. ; wm = 0.
      um(is,js,ks) = 1.e6
      vm(is,js,ks) = 1.e6
      wm(is,js,ks) = 1.e6

      call courant_local(dtm, got)
      courant_poison = same(label, got, 0.)

    end function courant_poison

    !> Put a single positive value in ekm (lekm) or ekh and assert the maximum
    !! diffusion number is that value times the cached geometry factor.
    logical function diffnr_spike(label, lekm, is, js, ks)
      character(len=*), intent(in) :: label
      logical,          intent(in) :: lekm
      integer,          intent(in) :: is, js, ks
      real :: got

      ekm = 0. ; ekh = 0.
      if (lekm) then
        ekm(is,js,ks) = spike
      else
        ekh(is,js,ks) = spike
      end if

      call diffnr_local(dtm, got)
      diffnr_spike = same(label, got, spike*diffnrgeom(is,ks)*dtm)

    end function diffnr_spike

    !> Put a large viscosity outside the reduction box and assert it is not seen.
    logical function diffnr_poison(label, is, js, ks)
      character(len=*), intent(in) :: label
      integer,          intent(in) :: is, js, ks
      real :: got

      ekm = 0. ; ekh = 0.
      ekm(is,js,ks) = 1.e6
      ekh(is,js,ks) = 1.e6

      call diffnr_local(dtm, got)
      diffnr_poison = same(label, got, 0.)

    end function diffnr_poison

  end function tests_checksim


  !> drivergen must leave every driver inlet plane current, in the one call.
  !!
  !! This is the host half of the invariant the GPU driver inlet rests on.
  !! initCUDA mirrors the twelve planes onto the device and fills them there
  !! and then - see modcuda.f90::allocDriverPlanesDevice - because by the time
  !! it runs, initdriver has read the driver file and the boundary call in the
  !! initialisation sequence has run drivergen over it. One upload at that
  !! point is a complete upload only if drivergen really does write all twelve
  !! planes, which is what the passes below pin down.
  !!
  !! The last pass is the reason the fill cannot be left to the time loop.
  !! drivergen writes the m planes only on rk3step 0 or 3, and boundary_device
  !! uploads on exactly those stages for that reason, so nothing uploads a
  !! complete set until the end of the first step. Allocating the device
  !! planes without filling them left the first two stages of the first step
  !! reading whatever the allocation landed on: case 452 came out of its first
  !! step with a divergence of 4.33 against 2.7E-13 on the host.
  !!
  !! Both readers are covered. With lchunkread off, drivergen indexes the
  !! whole record set directly. With it on, only a window of chunkread_size
  !! records is resident, slot 0 carries the last record of the previous
  !! window over, and drivergen maps the global record x to a slot with
  !!
  !!     xc = mod(x, chunkread_size)
  !!     if (xc == 0) xc = x - (chunkreadctr - 2)*chunkread_size
  !!
  !! which is its own way to be wrong. The chunk passes below put the window
  !! on the second chunk as well as the first, so the carry-over slot and the
  !! non-zero (chunkreadctr - 2) term both have to be right - on the first
  !! chunk that term is zero and the correction is indistinguishable from
  !! xc = x. Nothing about the transfer changes between the two: the record
  !! store never goes to the device, only the twelve planes drivergen builds
  !! from it, so the chunk boundary is invisible on the far side.
  !!
  !! The store is built here rather than read from a file, so the expected
  !! values follow from the interpolation arithmetic instead of from a fixture
  !! that would have to be regenerated with the code it checks. Every record
  !! is linear in the record index and every time below is a dyadic fraction,
  !! so the interpolated answer is exact and the tolerance only absorbs the
  !! rotation.
  !!
  !! Returns .true. if all checks pass, .false. otherwise.
  logical function tests_driver_planes()
    use modglobal,    only : runmode, jb, je, jh, kb, ke, kh, jhc, khc, nsv, &
                             ltempeq, lmoist, idriver, driverstore, lchunkread, &
                             lhdriver, lqdriver, lsdriver, timee, btime, runtime, &
                             rk3step, lwarmstart, driverid, chunkread_size
    use modinletdata, only : storetdriver, &
                             storeu0driver, storeumdriver, storev0driver, storevmdriver, &
                             storew0driver, storewmdriver, storethl0driver, storethlmdriver, &
                             storeqt0driver, storeqtmdriver, storesv0driver, storesvmdriver, &
                             u0driver, umdriver, v0driver, vmdriver, w0driver, wmdriver, &
                             thl0driver, thlmdriver, qt0driver, qtmdriver, &
                             sv0driver, svmdriver, u0driverrot, v0driverrot, &
                             iangle, nstepreaddriver, &
                             chunkreadctr, chunkread_s, chunkread_e
    use moddriver,    only : drivergen
    use decomp_2d,    only : zstart

    implicit none

    ! Eight records in two windows of four, so that the chunk passes can put
    ! the window on either one.
    integer, parameter :: nrec = 8, chunk = 4
    real,    parameter :: t_first = 0.25, dt_rec = 0.5
    real,    parameter :: c_u0 = 1., c_v0 = 3., c_w0 = 5., &
                          c_thl0 = 7., c_qt0 = 9., c_sv0 = 11.
    real,    parameter :: sentinel = -8888., poison = -7777.
    real,    parameter :: tol = 1e-12

    logical :: all_passed
    integer :: nreport

    ! A wrong plane is wrong at every point of it, and there are thousands.
    ! Name the first few and count the rest.
    nreport = 0

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_driver_planes: DRIVER INLET PLANE GENERATION TEST'
      write(*, '(A)') '------------------------------------------------'
    end if

    ! Do not pass vacuously on a namelist that switches the optional groups
    ! off - those are three of the twelve planes.
    if (.not. ltempeq .or. .not. lmoist .or. nsv < 1) then
      if (myid == 0) then
        write(*, '(A)') 'FAIL: this runmode needs ltempeq, lmoist and nsv > 0 to cover'
        write(*, '(A)') '  the temperature, moisture and scalar driver planes.'
      end if
      tests_driver_planes = .false.
      return
    end if

    ! Stand in for initdriver, which program.f90 runs long after this point
    ! and which would want a driver file on disk.
    idriver     = 2
    lchunkread  = .false.
    lhdriver    = .true.
    lqdriver    = .true.
    lsdriver    = .true.
    lwarmstart  = .false.
    driverid    = 0
    btime       = 0.
    runtime     = 4.
    nstepreaddriver = 0
    iangle      = 0.
    chunkread_size = chunk

    all_passed = .true.

    ! ---- the whole record set resident -----------------------------------
    lchunkread = .false.
    call alloc_store(1, nrec)
    call fill_store(0)

    ! rk3step 0 is the stage the initialisation-sequence call to boundary runs
    ! drivergen on, and the one whose result initCUDA uploads.
    rk3step = 0
    all_passed = one_pass('on a record',            1.25,  3.  ) .and. all_passed
    all_passed = one_pass('before the first',       0.125, 1.  ) .and. all_passed
    all_passed = one_pass('interpolated forward',   0.375, 1.25) .and. all_passed
    all_passed = one_pass('interpolated backward',  0.625, 1.75) .and. all_passed

    ! With an inflow angle, which is applied to u and v after the
    ! interpolation and before the m planes are taken. A copy made ahead of
    ! the rotation would leave um and u0 disagreeing.
    iangle = 0.375
    all_passed = one_pass('rotated inlet', 1.25, 3.) .and. all_passed
    iangle = 0.

    ! rk3step 3 is the stage the time loop uploads on.
    rk3step = 3
    all_passed = one_pass('last stage', 1.75, 4.) .and. all_passed

    all_passed = m_planes_held() .and. all_passed
    call drop_store

    ! ---- a window of the record set resident -----------------------------
    ! readdriverfile_chunk leaves chunkreadctr one past the window it just
    ! read, which is what drivergen's slot arithmetic is written against, so
    ! that is what is set here. The file reading itself is not exercised -
    ! what it produces is, which is the state below.
    lchunkread = .true.
    rk3step    = 0

    ! First window: records 1 to 4 in slots 1 to 4. Slot 0 would hold the
    ! record before the window, which does not exist here, so it is poisoned -
    ! no pass below has any business reading it.
    call alloc_store(0, chunk)
    call fill_store(0)
    call poison_slot(0)
    chunkreadctr = 2
    chunkread_s  = 1
    chunkread_e  = chunk
    all_passed = one_pass('chunk 1, on a record',           1.25,  3.  ) .and. all_passed
    all_passed = one_pass('chunk 1, before the first',      0.125, 1.  ) .and. all_passed
    all_passed = one_pass('chunk 1, interpolated forward',  0.375, 1.25) .and. all_passed
    all_passed = one_pass('chunk 1, interpolated backward', 0.625, 1.75) .and. all_passed
    ! x = 4 is a multiple of chunkread_size, so mod gives 0 and the correction
    ! runs. On this window its (chunkreadctr - 2) term is zero.
    all_passed = one_pass('chunk 1, last of the window',    1.75,  4.  ) .and. all_passed
    call drop_store

    ! Second window: records 5 to 8 in slots 1 to 4, with record 4 carried
    ! over into slot 0.
    call alloc_store(0, chunk)
    call fill_store(chunk)
    chunkreadctr = 3
    chunkread_s  = chunk + 1
    chunkread_e  = nrec
    all_passed = one_pass('chunk 2, on a record',           2.25,  5.  ) .and. all_passed
    ! Interpolating back across the window boundary: the older of the two
    ! records is the carry-over in slot 0, and nothing else reaches it.
    all_passed = one_pass('chunk 2, across the boundary',   2.125, 4.75) .and. all_passed
    all_passed = one_pass('chunk 2, interpolated forward',  2.375, 5.25) .and. all_passed
    ! x = 8, so mod gives 0 again - but here the correction's
    ! (chunkreadctr - 2)*chunkread_size term is 4, not 0, and dropping it
    ! would index slot 8 of a window that only has four.
    all_passed = one_pass('chunk 2, last of the window',    3.75,  8.  ) .and. all_passed
    call drop_store
    lchunkread = .false.

    if (myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      if (all_passed) then
        write(*, '(A)') 'tests_driver_planes: PASSED'
      else
        write(*, '(A,I0,A)') 'tests_driver_planes: FAILED (', nreport, ' points reported)'
      end if
      write(*, '(A)') '================================================'
    end if

    tests_driver_planes = all_passed

  contains

    !> The value stored for field c at record n, at (j,k).
    !!
    !! Linear in n, so linear interpolation between any two records is exact
    !! and the expected answer is the same expression at a fractional index.
    !! The global j is what stops a plane gathered for the other rank from
    !! passing on a two-rank run.
    real function stored(c, j, k, xn)
      real,    intent(in) :: c, xn
      integer, intent(in) :: j, k

      stored = 0.125*real(j + zstart(2) - 1) + 0.03125*real(k) + c + xn
    end function stored

    !> The record store, with slots lo..hi. For the whole-set reader that is
    !! 1..driverstore; for the chunk reader it is 0..chunkread_size, where
    !! slot 0 carries the record before the window.
    subroutine alloc_store(lo, hi)
      integer, intent(in) :: lo, hi
      integer :: m

      allocate(storetdriver(1:nrec))
      do m = 1, nrec
        storetdriver(m) = t_first + real(m - 1)*dt_rec
      end do
      driverstore = nrec

      allocate(storeu0driver (jb-jh:je+jh, kb-kh:ke+kh, lo:hi))
      allocate(storev0driver (jb-jh:je+jh, kb-kh:ke+kh, lo:hi))
      allocate(storew0driver (jb-jh:je+jh, kb-kh:ke+kh, lo:hi))
      allocate(storethl0driver(jb-jh:je+jh, kb-kh:ke+kh, lo:hi))
      allocate(storeqt0driver(jb-jh:je+jh, kb-kh:ke+kh, lo:hi))
      allocate(storesv0driver(jb-jhc:je+jhc, kb-khc:ke+khc, 1:nsv, lo:hi))

      ! drivergen takes the m planes from the 0 planes it has just built, not
      ! from the store. Poison the m store so that sourcing them from the file
      ! instead would show up as an m plane that does not match its 0 plane.
      allocate(storeumdriver (jb-jh:je+jh, kb-kh:ke+kh, lo:hi))
      allocate(storevmdriver (jb-jh:je+jh, kb-kh:ke+kh, lo:hi))
      allocate(storewmdriver (jb-jh:je+jh, kb-kh:ke+kh, lo:hi))
      allocate(storethlmdriver(jb-jh:je+jh, kb-kh:ke+kh, lo:hi))
      allocate(storeqtmdriver(jb-jh:je+jh, kb-kh:ke+kh, lo:hi))
      allocate(storesvmdriver(jb-jhc:je+jhc, kb-khc:ke+khc, 1:nsv, lo:hi))
      storeumdriver  = poison
      storevmdriver  = poison
      storewmdriver  = poison
      storethlmdriver = poison
      storeqtmdriver = poison
      storesvmdriver = poison

      call alloc_planes
    end subroutine alloc_store

    !> Put global record base + m into slot m, over whatever slots exist.
    !!
    !! base is 0 for the whole-set reader and for the first window, and
    !! chunkread_s - 1 for a later window - so slot 0 of the second window
    !! holds the last record of the first, which is what
    !! readdriverfile_chunk carries over.
    subroutine fill_store(base)
      integer, intent(in) :: base
      integer :: j, k, m, n

      do m = lbound(storeu0driver, 3), ubound(storeu0driver, 3)
        do k = kb-kh, ke+kh
          do j = jb-jh, je+jh
            storeu0driver  (j,k,m) = stored(c_u0,   j, k, real(base + m))
            storev0driver  (j,k,m) = stored(c_v0,   j, k, real(base + m))
            storew0driver  (j,k,m) = stored(c_w0,   j, k, real(base + m))
            storethl0driver(j,k,m) = stored(c_thl0, j, k, real(base + m))
            storeqt0driver (j,k,m) = stored(c_qt0,  j, k, real(base + m))
          end do
        end do
        do n = 1, nsv
          do k = kb-khc, ke+khc
            do j = jb-jhc, je+jhc
              storesv0driver(j,k,n,m) = stored(c_sv0 + 0.5*real(n), j, k, real(base + m))
            end do
          end do
        end do
      end do
    end subroutine fill_store

    !> Make one slot unreadable, for a slot no pass should reach.
    subroutine poison_slot(m)
      integer, intent(in) :: m

      storeu0driver  (:,:,m)   = poison
      storev0driver  (:,:,m)   = poison
      storew0driver  (:,:,m)   = poison
      storethl0driver(:,:,m)   = poison
      storeqt0driver (:,:,m)   = poison
      storesv0driver (:,:,:,m) = poison
    end subroutine poison_slot

    subroutine alloc_planes
      allocate(u0driver  (jb-jh:je+jh, kb-kh:ke+kh))
      allocate(umdriver  (jb-jh:je+jh, kb-kh:ke+kh))
      allocate(v0driver  (jb-jh:je+jh, kb-kh:ke+kh))
      allocate(vmdriver  (jb-jh:je+jh, kb-kh:ke+kh))
      allocate(w0driver  (jb-jh:je+jh, kb-kh:ke+kh))
      allocate(wmdriver  (jb-jh:je+jh, kb-kh:ke+kh))
      allocate(thl0driver(jb-jh:je+jh, kb-kh:ke+kh))
      allocate(thlmdriver(jb-jh:je+jh, kb-kh:ke+kh))
      allocate(qt0driver (jb-jh:je+jh, kb-kh:ke+kh))
      allocate(qtmdriver (jb-jh:je+jh, kb-kh:ke+kh))
      allocate(sv0driver (jb-jhc:je+jhc, kb-khc:ke+khc, 1:nsv))
      allocate(svmdriver (jb-jhc:je+jhc, kb-khc:ke+khc, 1:nsv))
      allocate(u0driverrot(jb-jh:je+jh, kb-kh:ke+kh))
      allocate(v0driverrot(jb-jh:je+jh, kb-kh:ke+kh))
    end subroutine alloc_planes

    subroutine drop_store
      deallocate(storetdriver)
      deallocate(storeu0driver, storev0driver, storew0driver, &
                 storethl0driver, storeqt0driver, storesv0driver)
      deallocate(storeumdriver, storevmdriver, storewmdriver, &
                 storethlmdriver, storeqtmdriver, storesvmdriver)
      deallocate(u0driver, umdriver, v0driver, vmdriver, w0driver, wmdriver, &
                 thl0driver, thlmdriver, qt0driver, qtmdriver, &
                 sv0driver, svmdriver, u0driverrot, v0driverrot)
    end subroutine drop_store

    subroutine sentinel_planes
      u0driver = sentinel ; umdriver = sentinel
      v0driver = sentinel ; vmdriver = sentinel
      w0driver = sentinel ; wmdriver = sentinel
      thl0driver = sentinel ; thlmdriver = sentinel
      qt0driver = sentinel ; qtmdriver = sentinel
      sv0driver = sentinel ; svmdriver = sentinel
    end subroutine sentinel_planes

    !> One drivergen call at time t, whose answer is record index xn.
    !!
    !! Three separate claims: that the 0 planes hold the interpolated record,
    !! that the m planes hold the same thing, and that nothing was left at the
    !! sentinel. The last is what makes a plane nobody thought to check here
    !! still able to fail.
    logical function one_pass(label, t, xn)
      character(len=*), intent(in) :: label
      real,             intent(in) :: t, xn

      real    :: ca, sa, eu, ev
      integer :: j, k, n

      call sentinel_planes
      timee = t
      call drivergen

      ca = cos(iangle)
      sa = sin(iangle)
      one_pass = .true.

      do k = kb-kh, ke+kh
        do j = jb-jh, je+jh
          eu = stored(c_u0, j, k, xn)
          ev = stored(c_v0, j, k, xn)
          one_pass = near(label, 'u0driver',   u0driver(j,k),   eu*ca - ev*sa) .and. one_pass
          one_pass = near(label, 'v0driver',   v0driver(j,k),   ev*ca + eu*sa) .and. one_pass
          one_pass = near(label, 'w0driver',   w0driver(j,k),   stored(c_w0,   j, k, xn)) .and. one_pass
          one_pass = near(label, 'thl0driver', thl0driver(j,k), stored(c_thl0, j, k, xn)) .and. one_pass
          one_pass = near(label, 'qt0driver',  qt0driver(j,k),  stored(c_qt0,  j, k, xn)) .and. one_pass

          ! The m planes are a straight copy, so exact equality is the claim.
          one_pass = same(label, 'umdriver',  umdriver(j,k),  u0driver(j,k))   .and. one_pass
          one_pass = same(label, 'vmdriver',  vmdriver(j,k),  v0driver(j,k))   .and. one_pass
          one_pass = same(label, 'wmdriver',  wmdriver(j,k),  w0driver(j,k))   .and. one_pass
          one_pass = same(label, 'thlmdriver', thlmdriver(j,k), thl0driver(j,k)) .and. one_pass
          one_pass = same(label, 'qtmdriver', qtmdriver(j,k), qt0driver(j,k))  .and. one_pass
        end do
      end do

      do n = 1, nsv
        do k = kb-khc, ke+khc
          do j = jb-jhc, je+jhc
            one_pass = near(label, 'sv0driver', sv0driver(j,k,n), &
                            stored(c_sv0 + 0.5*real(n), j, k, xn)) .and. one_pass
            one_pass = same(label, 'svmdriver', svmdriver(j,k,n), sv0driver(j,k,n)) .and. one_pass
          end do
        end do
      end do

      one_pass = no_sentinel(label) .and. one_pass
    end function one_pass

    !> Every one of the twelve planes must have been written.
    logical function no_sentinel(label)
      character(len=*), intent(in) :: label

      no_sentinel = .true.
      no_sentinel = written(label, 'u0driver',   any(u0driver   == sentinel)) .and. no_sentinel
      no_sentinel = written(label, 'umdriver',   any(umdriver   == sentinel)) .and. no_sentinel
      no_sentinel = written(label, 'v0driver',   any(v0driver   == sentinel)) .and. no_sentinel
      no_sentinel = written(label, 'vmdriver',   any(vmdriver   == sentinel)) .and. no_sentinel
      no_sentinel = written(label, 'w0driver',   any(w0driver   == sentinel)) .and. no_sentinel
      no_sentinel = written(label, 'wmdriver',   any(wmdriver   == sentinel)) .and. no_sentinel
      no_sentinel = written(label, 'thl0driver', any(thl0driver == sentinel)) .and. no_sentinel
      no_sentinel = written(label, 'thlmdriver', any(thlmdriver == sentinel)) .and. no_sentinel
      no_sentinel = written(label, 'qt0driver',  any(qt0driver  == sentinel)) .and. no_sentinel
      no_sentinel = written(label, 'qtmdriver',  any(qtmdriver  == sentinel)) .and. no_sentinel
      no_sentinel = written(label, 'sv0driver',  any(sv0driver  == sentinel)) .and. no_sentinel
      no_sentinel = written(label, 'svmdriver',  any(svmdriver  == sentinel)) .and. no_sentinel
    end function no_sentinel

    logical function written(label, name, left)
      character(len=*), intent(in) :: label, name
      logical,          intent(in) :: left

      written = .not. left
      if (left) then
        nreport = nreport + 1
        if (nreport <= 20 .and. myid == 0) then
          write(*,'(A,I0,A,A,A,A)') 'FAIL on rank ', myid, ': ', name, &
               ' left unwritten at ', label
        end if
      end if
    end function written

    !> The m planes move only on rk3step 0 and 3.
    !!
    !! Which is why boundary_device uploads on exactly those stages, and why
    !! the first fill has to come from initCUDA rather than from the loop:
    !! the first two stages of the first step never upload at all.
    logical function m_planes_held()
      real, allocatable :: keep_um(:,:), keep_thlm(:,:)
      real, allocatable :: keep_svm(:,:,:)

      m_planes_held = .true.

      rk3step = 3
      timee   = 1.25
      call drivergen
      allocate(keep_um,   source=umdriver)
      allocate(keep_thlm, source=thlmdriver)
      allocate(keep_svm,  source=svmdriver)

      rk3step = 1
      timee   = 1.75
      call drivergen

      if (any(umdriver /= keep_um) .or. any(thlmdriver /= keep_thlm) .or. &
          any(svmdriver /= keep_svm)) then
        nreport = nreport + 1
        if (myid == 0) write(*,'(A,I0,A)') 'FAIL on rank ', myid, &
             ': the m driver planes moved on rk3step 1'
        m_planes_held = .false.
      end if

      ! And the 0 planes did move, so the pass above is not passing because
      ! drivergen did nothing at all.
      if (all(u0driver == keep_um)) then
        nreport = nreport + 1
        if (myid == 0) write(*,'(A,I0,A)') 'FAIL on rank ', myid, &
             ': the 0 driver planes did not move on rk3step 1'
        m_planes_held = .false.
      end if

      deallocate(keep_um, keep_thlm, keep_svm)
      rk3step = 0
    end function m_planes_held

    logical function near(label, name, got, want)
      character(len=*), intent(in) :: label, name
      real,             intent(in) :: got, want

      near = abs(got - want) <= tol*max(1., abs(want))
      if (.not. near) then
        nreport = nreport + 1
        if (nreport <= 20 .and. myid == 0) then
          write(*,'(A,I0,A,A,A,A,A,ES22.14,A,ES22.14)') 'FAIL on rank ', myid, &
               ': ', name, ' at ', label, ' = ', got, ' want ', want
        end if
      end if
    end function near

    logical function same(label, name, got, want)
      character(len=*), intent(in) :: label, name
      real,             intent(in) :: got, want

      same = got == want
      if (.not. same) then
        nreport = nreport + 1
        if (nreport <= 20 .and. myid == 0) then
          write(*,'(A,I0,A,A,A,A,A,ES22.14,A,ES22.14)') 'FAIL on rank ', myid, &
               ': ', name, ' at ', label, ' = ', got, ' want ', want
        end if
      end if
    end function same

  end function tests_driver_planes

  !> Pin the facts about thermodynamics that the GPU port had to reproduce.
  !!
  !! Deliberately not a re-derivation of the physics. tests_cuda.f90's
  !! test_thermodynamics_device already drives the host routine and the device
  !! routine over one seed and requires everything they produce to agree, so
  !! it catches any disagreement between the two - and nothing that both get
  !! wrong together. This is exactly that blind spot: five facts about the
  !! host routine that are odd enough that someone would reasonably change
  !! one, where the device kernels are written to match rather than to be
  !! right.
  !!
  !! The first is the one that matters. thermo writes ql0 one k level low.
  !! Its ql dummy is declared kb:ke+kh while its qt and thl dummies are
  !! declared kb-kh:ke+kh - the commented-out declarations above the routine
  !! show ql being left behind when the other two were widened - and the first
  !! call passes the whole of ql0, whose own lower bound is kb-kh. Explicit
  !! shape means sequence association, so dummy level kb lands on actual level
  !! kb-kh: the saturation computed from level k is stored at level k-kh, and
  !! ql0's top level is never written at all. ql0h is unaffected, because its
  !! lower bound is kb and the second call is aligned.
  !!
  !! That shift reaches ql0av, calthv and the ql field in a fielddump, so it
  !! is not cosmetic - but correcting it is a physics change to the CPU solver
  !! and belongs in its own commit. Until then thermo_device declares its
  !! dummies character for character the same way, so the two shift together
  !! and stay comparable. If this test starts failing because the host was
  !! corrected, the device declaration has to be corrected in the same commit.
  logical function tests_thermodynamics()
    use modglobal,  only : runmode, ib, ie, ih, jb, je, jh, kb, ke, kh, &
                           ltempeq, lmoist, timee, &
                           dzf, dzh, eps1, es0, at, bt, rd, rv, rlv, cp, tmelt
    use modfields,  only : initfields, &
                           u0, v0, thl0, qt0, ql0, ql0h, thl0h, qt0h, thv0h, dthvdz, &
                           presf, exnf
    use modsurfdata,only : thls, qts
    use modthermodynamics, only : thermodynamics, initthermodynamics
    use initfac,    only : readfacetfiles
    use modibm,     only : initibm, createmasks

    implicit none

    real, parameter :: sentinel = -6543.
    ! Every comparison below is between two evaluations of the same expression
    ! on the same host, so the slack is for the compiler's freedom to keep an
    ! intermediate in a wider register, not for a numerical difference.
    real, parameter :: tol = 1e-12

    real, allocatable :: presf_seed(:), exnf_seed(:)
    logical :: all_passed
    integer :: i, j, k, nreport
    real    :: want, got

    nreport = 0

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_thermodynamics: HOST THERMODYNAMICS CONTRACT TEST'
      write(*, '(A)') '------------------------------------------------'
    end if

    ! Four of the five checks are about the moist branches. Fail rather than
    ! pass vacuously on a namelist that switches them off.
    if (.not. ltempeq .or. .not. lmoist) then
      if (myid == 0) then
        write(*, '(A)') 'FAIL: this runmode needs ltempeq and lmoist - four of'
        write(*, '(A)') '  the five checks are about the saturated branches.'
      end if
      tests_thermodynamics = .false.
      return
    end if
    if (qts <= 0.) then
      if (myid == 0) then
        write(*, '(A)') 'FAIL: qts must be set in &BC. fromztop divides by a'
        write(*, '(A)') '  surface virtual temperature built from it.'
      end if
      tests_thermodynamics = .false.
      return
    end if

    call initfields
    call initthermodynamics
    call readfacetfiles
    call initibm
    call createmasks

    allocate(presf_seed(kb:ke+kh), exnf_seed(kb:ke+kh))

    all_passed = .true.

    ! ---- the state the first saturation pass sees -------------------------
    !
    ! timee is set away from zero so that the leading diagfld does not run.
    ! That is what makes the reference computable: the first thermo then reads
    ! exactly the column seeded here, where after a diagfld it would read one
    ! this test cannot reproduce without reimplementing fromztop.
    timee = 1.
    do k = kb, ke+kh
      presf_seed(k) = 101500. - 100.*real(k-kb)
      exnf_seed(k)  = 1. - 0.0009765625*real(k-kb)
    end do
    presf = presf_seed
    exnf  = exnf_seed

    call seed_fields

    ql0  = sentinel
    ql0h = sentinel
    ! The halo columns of thv0h, which the saturated branch of calthv does not
    ! write - see check 5.
    thv0h = sentinel

    call thermodynamics

    if (.not. check_saturation_shift()) all_passed = .false.
    if (.not. check_dim_is_max())       all_passed = .false.
    if (.not. check_half_levels())      all_passed = .false.
    if (.not. check_dthvdz())           all_passed = .false.
    if (.not. check_thv0h_halo_moist()) all_passed = .false.

    ! ---- and the branch a moist namelist cannot otherwise reach -----------
    lmoist = .false.
    thv0h  = sentinel
    call thermodynamics
    if (.not. check_thv0h_halo_dry()) all_passed = .false.
    lmoist = .true.

    if (myid == 0) then
      if (all_passed) then
        write(*, '(A)') 'tests_thermodynamics: PASS'
      else
        write(*, '(A)') 'tests_thermodynamics: FAIL'
      end if
      write(*, '(A)') '================================================'
    end if

    deallocate(presf_seed, exnf_seed)
    tests_thermodynamics = all_passed

  contains

    !> A profile with a real vertical gradient and a horizontal pattern.
    !!
    !! qt0 alternates by level between one value no saturation specific
    !! humidity reaches and one every plausible value sits below, so both arms
    !! of the saturated test are taken and neither is near its threshold.
    subroutine seed_fields
      implicit none
      integer :: ii, jj, kk

      do kk = kb-kh, ke+kh
        do jj = jb-jh, je+jh
          do ii = ib-ih, ie+ih
            u0(ii,jj,kk)   =  1.0 + wiggle(ii,jj,kk)
            v0(ii,jj,kk)   = -0.5 + wiggle(ii,jj,kk)
            thl0(ii,jj,kk) = 290. + 0.25*real(kk-kb) + wiggle(ii,jj,kk)
            if (modulo(kk,2) == 0) then
              qt0(ii,jj,kk) = 0.03125    + 0.00390625*wiggle(ii,jj,kk)
            else
              qt0(ii,jj,kk) = 0.00390625 + 0.00048828125*wiggle(ii,jj,kk)
            end if
          end do
        end do
      end do

    end subroutine seed_fields

    real function wiggle(ii, jj, kk)
      implicit none
      integer, intent(in) :: ii, jj, kk

      wiggle = 0.0625*real(modulo(ii*7 + jj*13 + kk*29, 9)) - 0.25

    end function wiggle

    !> thermo's saturation deficit, spelled the way thermo spells it.
    real function sat_deficit(thl, qt, pressure, exner)
      implicit none
      real, intent(in) :: thl, qt, pressure, exner

      real :: tl, es, qs, qsl, b1

      tl = thl*exner
      if (tl<100.0) tl = 100.0
      es  = es0*exp(at*(tl-tmelt)/(tl-bt))
      qsl = rd/rv*es/(pressure-(1-rd/rv)*es)
      b1  = rlv**2/(tl**2*cp*rv)
      qs  = qsl*(1.+b1*qt)/(1.+b1*qsl)
      sat_deficit = dim(qt-qs, 0.)

    end function sat_deficit

    !> The saturation computed from level k is stored at level k-kh.
    logical function check_saturation_shift()
      implicit none

      check_saturation_shift = .true.
      nreport = 0

      do k = kb, ke+kh
        do j = jb, je
          do i = ib, ie
            want = sat_deficit(thl0(i,j,k), qt0(i,j,k), presf_seed(k), exnf_seed(k))
            got  = ql0(i,j,k-kh)
            if (abs(got - want) > tol*max(1e-6, abs(want))) then
              call flag('ql0 at the shifted level', i, j, k-kh, got, want)
              check_saturation_shift = .false.
            end if
          end do
        end do
      end do

      ! The top level is the other half of the same fact: the loop runs to
      ! ke+kh on a dummy that starts one level low, so nothing reaches it.
      do j = jb, je
        do i = ib, ie
          if (ql0(i,j,ke+kh) /= sentinel) then
            call flag('ql0 top level written', i, j, ke+kh, ql0(i,j,ke+kh), sentinel)
            check_saturation_shift = .false.
          end if
        end do
      end do

      ! And ql0h is not shifted, so every level of it is written.
      do k = kb, ke+kh
        do j = jb, je
          do i = ib, ie
            if (ql0h(i,j,k) == sentinel) then
              call flag('ql0h level not written', i, j, k, ql0h(i,j,k), sentinel)
              check_saturation_shift = .false.
            end if
          end do
        end do
      end do

      call verdict('saturation is stored one level low', check_saturation_shift)

    end function check_saturation_shift

    !> dim(x,0.) and max(x,0.) agree on the values this field produces.
    !!
    !! thermo_device spells the clamp as max because DIM has no guaranteed
    !! device runtime, and the standard defines DIM(x,y) as x-y when x is
    !! greater and zero otherwise - which is what MAX(x-y,0) evaluates to.
    !! Written down here so the substitution is not a claim about the compiler.
    logical function check_dim_is_max()
      implicit none

      real :: d, x

      check_dim_is_max = .true.
      nreport = 0

      do k = kb, ke+kh
        do j = jb, je
          do i = ib, ie
            x = qt0(i,j,k) - (qt0(i,j,k) - ql0(i,j,k-kh))
            d = dim(x, 0.)
            if (d /= max(x, 0.)) then
              call flag('dim and max disagree', i, j, k, d, max(x, 0.))
              check_dim_is_max = .false.
            end if
          end do
        end do
      end do

      call verdict('dim(x,0.) is max(x,0.)', check_dim_is_max)

    end function check_dim_is_max

    !> calc_halflev interpolates by cell height and overrides the surface.
    logical function check_half_levels()
      implicit none

      check_half_levels = .true.
      nreport = 0

      do j = jb, je
        do i = ib, ie
          if (thl0h(i,j,kb) /= thls) then
            call flag('thl0h surface override', i, j, kb, thl0h(i,j,kb), thls)
            check_half_levels = .false.
          end if
          if (qt0h(i,j,kb) /= qts) then
            call flag('qt0h surface override', i, j, kb, qt0h(i,j,kb), qts)
            check_half_levels = .false.
          end if
        end do
      end do

      do k = kb+1, ke+kh
        do j = jb, je
          do i = ib, ie
            want = (thl0(i,j,k)*dzf(k-1)+thl0(i,j,k-1)*dzf(k))/(2*dzh(k))
            if (abs(thl0h(i,j,k) - want) > tol*max(1e-6, abs(want))) then
              call flag('thl0h interpolation', i, j, k, thl0h(i,j,k), want)
              check_half_levels = .false.
            end if
            want = (qt0(i,j,k)*dzf(k-1)+qt0(i,j,k-1)*dzf(k))/(2*dzh(k))
            if (abs(qt0h(i,j,k) - want) > tol*max(1e-6, abs(want))) then
              call flag('qt0h interpolation', i, j, k, qt0h(i,j,k), want)
              check_half_levels = .false.
            end if
          end do
        end do
      end do

      call verdict('half levels interpolate and take the surface at kb', check_half_levels)

    end function check_half_levels

    !> The eps1 clamp reaches every level it is applied to, and no other.
    !!
    !! The lowest level is where this bites: calthv writes a hard zero there
    !! and the clamp then turns it into +eps1, because sign(eps1, 0.) is
    !! +eps1. subgrid divides by dthvdz, so a level left at zero would be a
    !! division by zero rather than a small number - which is the whole point
    !! of the clamp and the reason the device kernel applies it over the same
    !! range rather than over the range it writes.
    logical function check_dthvdz()
      implicit none

      check_dthvdz = .true.
      nreport = 0

      do j = jb, je
        do i = ib, ie
          if (dthvdz(i,j,kb) /= eps1) then
            call flag('dthvdz at kb is not the clamped zero', i, j, kb, dthvdz(i,j,kb), eps1)
            check_dthvdz = .false.
          end if
        end do
      end do

      do k = kb, ke
        do j = jb, je
          do i = ib, ie
            if (abs(dthvdz(i,j,k)) < eps1) then
              call flag('dthvdz below the clamp', i, j, k, dthvdz(i,j,k), eps1)
              check_dthvdz = .false.
            end if
          end do
        end do
      end do

      ! ke+kh is outside the clamp's range and outside the range calthv
      ! writes, so it keeps the whole-array zero.
      do j = jb, je
        do i = ib, ie
          if (dthvdz(i,j,ke+kh) /= 0.) then
            call flag('dthvdz above the clamp', i, j, ke+kh, dthvdz(i,j,ke+kh), 0.)
            check_dthvdz = .false.
          end if
        end do
      end do

      call verdict('dthvdz is clamped over kb to ke and nowhere else', check_dthvdz)

    end function check_dthvdz

    !> The saturated branch writes the interior of thv0h and not its halo.
    logical function check_thv0h_halo_moist()
      implicit none

      check_thv0h_halo_moist = .true.
      nreport = 0

      do k = kb, ke+kh
        do j = jb-jh, je+jh
          do i = ib-ih, ie+ih
            if (i >= ib .and. i <= ie .and. j >= jb .and. j <= je) cycle
            if (thv0h(i,j,k) /= sentinel) then
              call flag('moist calthv wrote a thv0h halo cell', i, j, k, thv0h(i,j,k), sentinel)
              check_thv0h_halo_moist = .false.
            end if
          end do
        end do
      end do

      call verdict('the saturated branch leaves the thv0h halo alone', check_thv0h_halo_moist)

    end function check_thv0h_halo_moist

    !> The unsaturated branch assigns thv0h whole, halo included.
    logical function check_thv0h_halo_dry()
      implicit none

      check_thv0h_halo_dry = .true.
      nreport = 0

      do k = kb, ke+kh
        do j = jb-jh, je+jh
          do i = ib-ih, ie+ih
            if (thv0h(i,j,k) /= thl0h(i,j,k)) then
              call flag('dry calthv did not copy thl0h', i, j, k, thv0h(i,j,k), thl0h(i,j,k))
              check_thv0h_halo_dry = .false.
            end if
          end do
        end do
      end do

      call verdict('the unsaturated branch copies thl0h including the halo', check_thv0h_halo_dry)

    end function check_thv0h_halo_dry

    !> Name the first few failures and count the rest: a wrong level is wrong
    !! at every point of it, and there are tens of thousands.
    subroutine flag(what, ii, jj, kk, g, w)
      implicit none
      character(len=*), intent(in) :: what
      integer,          intent(in) :: ii, jj, kk
      real,             intent(in) :: g, w

      nreport = nreport + 1
      if (nreport <= 20 .and. myid == 0) then
        write(*, '(A,A,A,I0,A,I0,A,I0,A,ES22.14,A,ES22.14)') '  FAIL: ', trim(what), &
             ' at (', ii, ',', jj, ',', kk, ') got ', g, ' want ', w
      end if

    end subroutine flag

    subroutine verdict(what, passed)
      implicit none
      character(len=*), intent(in) :: what
      logical,          intent(in) :: passed

      if (myid /= 0) return
      if (passed) then
        write(*, '(A,A)') '  ok: ', trim(what)
      else
        write(*, '(A,A,A,I0,A)') '  FAILED: ', trim(what), ' (', nreport, ' points)'
      end if

    end subroutine verdict

  end function tests_thermodynamics


  !> The two rank-local reductions that choose the adaptive time step.
  !!
  !! tstep_limits is what tstep_update evaluates over the interior of the
  !! domain before it picks dt; everything else in that routine is scalar
  !! bookkeeping over clocks and counters. Driving the reduction directly is
  !! what lets a test see the answer at all - through tstep_update it is
  !! visible only as a shift in dt, after an all-reduce, a division by a
  !! namelist Courant number and a min against two other limits.
  !!
  !! What it pins, over and above the equivalent checks on modchecksim's pair,
  !! which look like the same two numbers and are not:
  !!
  !!   - The absolute values. tstep_update takes |u|, |v| and |w|; checksim's
  !!     Courant number does not. A port that drops them reports a negative
  !!     maximum as zero, which on a flow with a strong downdraught and a weak
  !!     updraught gives a Courant number several times too small and a dt
  !!     several times too large - a run that goes unstable rather than one
  !!     that is quietly a little wrong.
  !!
  !!   - The grid factors this pair carries, which are the cell-centred dxi
  !!     and dyi and the half-level dzh, not checksim's dxhi and dzhi.
  !!
  !!   - That the geometry cache is exactly the expression it replaced, and
  !!     that the diffusion term reads it at k.
  !!
  !!   - That the two numbers are independent: a velocity must not move the
  !!     diffusion number and a diffusivity must not move the Courant number.
  !!     They come out of one fused loop, so a reduction variable used for both
  !!     is a real way to get this wrong and would otherwise show up only as an
  !!     over-small dt on a run with a large viscosity.
  logical function tests_tstep()
    use modglobal,      only : runmode, ib, ie, jb, je, kb, ke, &
                               dxi, dyi, dx2i, dy2i, dzh, dzh2i
    use modfields,      only : initfields, um, vm, wm
    use modsubgrid,     only : initsubgrid
    use modsubgriddata, only : ekm, ekh
    use modtstep,       only : inittstep, tstep_limits, diffgeom

    implicit none

    real, parameter :: spike = 2., dtm = 0.25, poison = 1.e6

    logical :: all_passed
    integer :: k
    real    :: cour, diff, want
    real    :: dxi_save, dyi_save
    real, allocatable :: dzh_save(:), geom_save(:)

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_tstep: ADAPTIVE TIME STEP LIMITS TEST'
      write(*, '(A)') '------------------------------------------------'
    end if

#if defined(_GPU)
    ! Under _GPU tstep_limits reduces over the device fields, so seeding the
    ! host arrays this test writes would leave it looking at whatever the
    ! device happens to hold. Say so rather than reporting an unexplained
    ! mismatch, and do not pass vacuously.
    if (myid == 0) then
      write(*, '(A)') 'FAIL: runmode 1017 exercises the host branch of tstep_limits.'
      write(*, '(A)') '  Run it against a CPU build (build/cpu/<type>/u-dales).'
      write(*, '(A)') '  The device kernel is covered by tests_cuda.f90::test_tstep_limits,'
      write(*, '(A)') '  run with UDALES_RUN_CUDA_SELFTEST=1 on a Debug GPU build.'
    end if
    tests_tstep = .false.
    return
#endif

    call initfields
    call initsubgrid
    call inittstep

    all_passed = .true.

    ! ---- the geometry cache -------------------------------------------------
    ! Exact equality, not a tolerance: the whole justification for precomputing
    ! this was that it is the same arithmetic in the same order, so that the
    ! diffusion number, and with it dt, does not move by a single bit.
    do k = kb, ke
      want = dzh2i(k) + dx2i + dy2i
      if (diffgeom(k) /= want) then
        write(*,'(A,I0,A,I0,A,ES22.14,A,ES22.14)') 'FAIL on rank ', myid, &
             ': diffgeom(', k, ') = ', diffgeom(k), ' want ', want
        all_passed = .false.
      end if
    end do

    ! ---- make every grid factor distinct ------------------------------------
    ! Every case in the repository has dx = dy, and all but one a uniform grid
    ! in z. There dxi equals dyi and dzh is constant in k, so a loop that pairs
    ! the v term with the x spacing, or reads dzh at a fixed index, returns
    ! exactly the right answer and no assertion below can fail. Substituting
    ! values that differ in every position is what makes the pairing and the
    ! indexing load-bearing.
    !
    ! The cache is substituted rather than rebuilt from the new dzh2i, so that
    ! it varies in k independently of dzh: a diffusion term that reached for
    ! the Courant number's spacing instead would otherwise still agree.
    !
    ! This runs after the exactness check above, which has to see the grid
    ! inittstep actually built, and everything is restored before returning.
    allocate(dzh_save (lbound(dzh,1):ubound(dzh,1)))           ; dzh_save  = dzh
    allocate(geom_save(lbound(diffgeom,1):ubound(diffgeom,1))) ; geom_save = diffgeom
    dxi_save = dxi ; dyi_save = dyi

    ! Strictly increasing, and exact in binary: no two indices can collide, so
    ! reading either of these one index off is always visible. A cyclic pattern
    ! would not be - with a period of seven and a 64-deep domain, dzh(kb) and
    ! dzh(ke) come out equal and a kernel pinned to kb passes.
    do k = lbound(dzh,1), ubound(dzh,1)
      dzh(k) = 0.75 + 0.0078125*real(k - lbound(dzh,1))
    end do
    do k = lbound(diffgeom,1), ubound(diffgeom,1)
      diffgeom(k) = 1.25 + 0.0078125*real(k - lbound(diffgeom,1))
    end do
    dxi = 0.25 ; dyi = 0.5

    ! ---- the Courant number -------------------------------------------------
    call zero_all
    call tstep_limits(dtm, cour, diff)
    if (.not. same('courant of a field at rest', cour, 0.)) all_passed = .false.
    if (.not. same('diffusion number of zero viscosity', diff, 0.)) all_passed = .false.

    ! One term at a time, so each grid factor is pinned separately. A port that
    ! pairs vm with 1/dzh, say, reproduces none of these.
    if (.not. cour_spike('courant um', 1, ie, je, ke,  spike*dxi       )) all_passed = .false.
    if (.not. cour_spike('courant vm', 2, ie, je, ke,  spike*dyi       )) all_passed = .false.
    if (.not. cour_spike('courant wm', 3, ie, je, ke,  spike/dzh(ke)   )) all_passed = .false.

    ! The absolute values, one per term. Without them each of these comes back
    ! zero, because the reduction starts at zero and a negative contribution
    ! never displaces it.
    if (.not. cour_spike('courant |um|', 1, ie, je, ke, spike*dxi,     -spike)) all_passed = .false.
    if (.not. cour_spike('courant |vm|', 2, ie, je, ke, spike*dyi,     -spike)) all_passed = .false.
    if (.not. cour_spike('courant |wm|', 3, ie, je, ke, spike/dzh(ke), -spike)) all_passed = .false.

    ! The same spike at every bound of the reduction box, so no face of it can
    ! be missing. wm, because its factor is the one that varies with k.
    if (.not. cour_spike('courant corner ib,jb,kb', 3, ib, jb, kb, spike/dzh(kb))) all_passed = .false.
    if (.not. cour_spike('courant corner ie,jb,kb', 3, ie, jb, kb, spike/dzh(kb))) all_passed = .false.
    if (.not. cour_spike('courant corner ib,je,kb', 3, ib, je, kb, spike/dzh(kb))) all_passed = .false.
    if (.not. cour_spike('courant corner ib,jb,ke', 3, ib, jb, ke, spike/dzh(ke))) all_passed = .false.

    ! ... and just outside every one of those bounds, where it must be ignored.
    if (.not. poison_at('below ib', ib-1, jb,   kb  )) all_passed = .false.
    if (.not. poison_at('above ie', ie+1, jb,   kb  )) all_passed = .false.
    if (.not. poison_at('below jb', ib,   jb-1, kb  )) all_passed = .false.
    if (.not. poison_at('above je', ib,   je+1, kb  )) all_passed = .false.
    if (.not. poison_at('below kb', ib,   jb,   kb-1)) all_passed = .false.
    if (.not. poison_at('above ke', ib,   jb,   ke+1)) all_passed = .false.

    ! ---- the diffusion number -----------------------------------------------
    ! ekm and ekh one at a time with the other zeroed, so dropping either from
    ! the max is caught. A run whose Prandtl number is below one is limited by
    ! ekh, and that is the term a port is most likely to lose.
    if (.not. diff_spike('diffusion ekm', .true.,  ie, je, ke)) all_passed = .false.
    if (.not. diff_spike('diffusion ekh', .false., ie, je, ke)) all_passed = .false.
    if (.not. diff_spike('diffusion ekm at ib,jb,kb', .true.,  ib, jb, kb)) all_passed = .false.
    if (.not. diff_spike('diffusion ekh at ib,jb,kb', .false., ib, jb, kb)) all_passed = .false.

    ! ---- the two numbers are independent ------------------------------------
    ! One loop produces both, so a single reduction variable serving both terms
    ! passes every check above and fails these two.
    call zero_all
    um(ie,je,ke) = spike
    call tstep_limits(dtm, cour, diff)
    if (.not. same('a velocity does not move the diffusion number', diff, 0.)) all_passed = .false.
    call zero_all
    ekm(ie,je,ke) = spike
    call tstep_limits(dtm, cour, diff)
    if (.not. same('a diffusivity does not move the courant number', cour, 0.)) all_passed = .false.

    ! ---- the time step ------------------------------------------------------
    ! It multiplies the whole of each expression, and both of them. The
    ! solver's loop used to multiply every cell by dt inside the maximum and
    ! now multiplies the maximum once; these say the answer is the same, and
    ! at exact equality, because that claim is a bit-for-bit one.
    call zero_all
    um(ie,je,ke) = spike
    ekh(ib,jb,kb) = spike
    call tstep_limits(2.*dtm, cour, diff)
    if (cour /= (spike*dxi)*(2.*dtm)) then
      write(*,'(A,I0,A,ES22.14,A,ES22.14)') 'FAIL on rank ', myid, &
           ': tstep courant scales with dt, got ', cour, ' want ', (spike*dxi)*(2.*dtm)
      all_passed = .false.
    end if
    if (diff /= (spike*diffgeom(kb))*(2.*dtm)) then
      write(*,'(A,I0,A,ES22.14,A,ES22.14)') 'FAIL on rank ', myid, &
           ': tstep diffusion scales with dt, got ', diff, ' want ', (spike*diffgeom(kb))*(2.*dtm)
      all_passed = .false.
    end if

    dzh = dzh_save ; diffgeom = geom_save
    dxi = dxi_save ; dyi = dyi_save
    deallocate(geom_save, dzh_save)

    if (myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      if (all_passed) then
        write(*, '(A)') 'ALL TESTS PASSED: tests_tstep'
        write(*, '(A,I0,A,I0,A,I0,A)') '  Reduction box ', ie-ib+1, ' x ', je-jb+1, &
             ' x ', ke-kb+1, ' cells per rank'
      else
        write(*, '(A)') 'TESTS FAILED: tests_tstep'
        write(*, '(A)') '  One or more checks did not pass'
      end if
      write(*, '(A)') '================================================'
    end if

    tests_tstep = all_passed

  contains

    !> Clear everything the reduction reads, halos included.
    subroutine zero_all
      um = 0. ; vm = 0. ; wm = 0. ; ekm = 0. ; ekh = 0.
    end subroutine zero_all

    !> Assert a scalar matches to within a few ulps of its own magnitude.
    logical function same(label, got, want)
      character(len=*), intent(in) :: label
      real,             intent(in) :: got, want
      real :: tol

      tol  = 64. * epsilon(1.) * max(1., abs(want))
      same = abs(got - want) <= tol
      if (.not. same) then
        write(*,'(A,I0,3A,ES22.14,A,ES22.14)') 'FAIL on rank ', myid, ': tstep ', &
             trim(label), ' got ', got, ' want ', want
      end if

    end function same

    !> Put one value in um, vm or wm - which selects it - and assert the
    !! Courant number is want, the grid factor that term is supposed to carry
    !! times the value planted. val defaults to spike; pass a negative one to
    !! drive the absolute value.
    logical function cour_spike(label, which, is, js, ks, want, val)
      character(len=*), intent(in) :: label
      integer,          intent(in) :: which, is, js, ks
      real,             intent(in) :: want
      real, optional,   intent(in) :: val
      real :: put

      put = spike
      if (present(val)) put = val

      call zero_all
      select case (which)
      case (1) ; um(is,js,ks) = put
      case (2) ; vm(is,js,ks) = put
      case (3) ; wm(is,js,ks) = put
      end select

      call tstep_limits(dtm, cour, diff)
      cour_spike = same(label, cour, want*dtm)

    end function cour_spike

    !> Put a single positive value in ekm (lekm) or ekh and assert the
    !! diffusion number is that value times the cached geometry factor at ks.
    logical function diff_spike(label, lekm, is, js, ks)
      character(len=*), intent(in) :: label
      logical,          intent(in) :: lekm
      integer,          intent(in) :: is, js, ks

      call zero_all
      if (lekm) then
        ekm(is,js,ks) = spike
      else
        ekh(is,js,ks) = spike
      end if

      call tstep_limits(dtm, cour, diff)
      diff_spike = same(label, diff, spike*diffgeom(ks)*dtm)

    end function diff_spike

    !> Put a large value in all five fields outside the reduction box and
    !! assert neither number sees it.
    logical function poison_at(label, is, js, ks)
      character(len=*), intent(in) :: label
      integer,          intent(in) :: is, js, ks

      call zero_all
      um(is,js,ks) = poison ; vm(is,js,ks) = poison ; wm(is,js,ks) = poison
      ekm(is,js,ks) = poison ; ekh(is,js,ks) = poison

      call tstep_limits(dtm, cour, diff)
      poison_at = same('courant ' // label, cour, 0.)
      poison_at = same('diffusion ' // label, diff, 0.) .and. poison_at

    end function poison_at

  end function tests_tstep

  !> Time interpolation of the prescribed surface fluxes and nudging profiles.
  !!
  !! timedep itself is four interpolations and a bracket search, and every one
  !! of them is invisible from the outside: the surface fluxes reach the
  !! solution through a wall function, the nudging profiles through a tendency
  !! divided by tnudge, and both are then integrated. A parity case can tell
  !! that something moved; only this can tell that it moved to the right value
  !! at the right time.
  !!
  !! Neither routine is called through timedep here. timedep is a dispatcher
  !! over four switches and calling it would mean staging all four tables to
  !! test any one of them, with the two energy-balance branches needing an
  !! initialised facet set on top of that.
  !!
  !! There is no input file in a runmode test, so the tables are installed
  !! directly. That is also what makes the checks below sharp: the tables read
  !! from a file are smooth in k and in t, and a smooth table hides a
  !! transposed index. Every entry here differs from every other one.
  !!
  !! What it pins:
  !!
  !!   - The bracket. The search walks down from the last entry and takes the
  !!     first whose time has passed, so an interval boundary belongs to the
  !!     interval above it. Off by one and the whole forcing runs one interval
  !!     early or late, which on a diurnal file is hours.
  !!
  !!   - The endpoints, exactly. At a table time fac is zero and the result is
  !!     the table entry itself, bit for bit - a + 0*(b - a) is a in IEEE
  !!     arithmetic for any finite a and b. So these are checked with /=, and
  !!     an interpolation rewritten as (1-fac)*a + fac*b, which is the obvious
  !!     tidy-up and is not the same expression, fails them.
  !!
  !!   - The freeze, and where it starts. The bracket search returns the last
  !!     index once timee reaches the last table time, and the guard then
  !!     writes nothing - so the last column is never installed, and what
  !!     stands from there on is the interpolated value from the step before.
  !!     In a run that is the last column to within one step, which is why
  !!     this has never mattered; it is still not the same thing, and a port
  !!     that clamped to the last column instead would pass every other check
  !!     here. A run longer than its forcing file is the normal case, not an
  !!     edge case.
  !!
  !!   - That each of the five fluxes and each of the four profiles carries its
  !!     own column of the table, at its own level. Five scalars of the same
  !!     magnitude interpolated in one routine is exactly the shape a
  !!     copy-paste error survives in.
  logical function tests_timedep()
    use modglobal,  only : runmode, kb, ke, kh, timee
    use modfields,  only : initfields, thlprof, qtprof, uprof, vprof
    use modibmdata, only : bctfxm, bctfxp, bctfym, bctfyp, bctfz
    use modtimedep, only : timedepsurf, timedepnudge, &
                           ltimedepsurf, ntimedepsurf, ltimedepnudge, ntimedepnudge, &
                           timeflux, bctfxmt, bctfxpt, bctfymt, bctfypt, bctfzt, &
                           timenudge, thlproft, qtproft, uproft, vproft

    implicit none

    ! Three entries: enough for two intervals, so that the bracket search has
    ! to choose rather than having one answer. The times are not evenly spaced
    ! and not integers, so a fac formed from the wrong pair of entries lands
    ! somewhere visibly wrong rather than nearly right.
    integer, parameter :: nt = 3
    real,    parameter :: t1 = 0., t2 = 10., t3 = 25.

    logical :: all_passed
    real    :: fac

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_timedep: TIME-DEPENDENT FORCING TEST'
      write(*, '(A)') '------------------------------------------------'
    end if

    call initfields

    all_passed = .true.

    ! ---- surface fluxes -----------------------------------------------------
    ltimedepsurf = .true.
    ntimedepsurf = nt
    allocate(timeflux(nt), bctfxmt(nt), bctfxpt(nt), bctfymt(nt), bctfypt(nt), bctfzt(nt))

    timeflux = [ t1, t2, t3 ]
    ! One decade apart between the five, so a flux that ends up in the wrong
    ! variable is off by a factor of ten rather than by a few percent, and
    ! every entry differs from every other one in the table.
    bctfxmt = [  1.,  2.,  4. ]
    bctfxpt = [ 10., 20., 40. ]
    bctfymt = [ 100., 200., 400. ]
    bctfypt = [ 1000., 2000., 4000. ]
    bctfzt  = [ 10000., 20000., 40000. ]

    ! At the first table time. fac is zero, so this is the exact first entry
    ! and also the first bracket - the search must not run off the bottom.
    timee = t1
    call timedepsurf
    all_passed = surf_exact('t = t1', 1) .and. all_passed

    ! Inside the first interval, at a fraction that is not a half: a fac
    ! formed as 0.5, or from the second interval's width, misses.
    timee = t1 + 0.25*(t2 - t1)
    call timedepsurf
    all_passed = surf_interp('inside first interval', 1, 0.25) .and. all_passed

    ! Exactly on the interior table time. The search takes the interval above,
    ! so this is entry 2 with fac = 0, not entry 1 with fac = 1. Both give the
    ! same number here, which is the point: it is the one time the two
    ! bracketings agree, so the check that separates them is the one below.
    timee = t2
    call timedepsurf
    all_passed = surf_exact('t = t2', 2) .and. all_passed

    ! Inside the second interval. The two intervals have different widths, so
    ! a fac divided by the wrong one is wrong here even at the same offset.
    timee = t2 + 0.25*(t3 - t2)
    call timedepsurf
    all_passed = surf_interp('inside second interval', 2, 0.25) .and. all_passed

    ! At the last table time, and past it, nothing is written at all - the
    ! search returns ntimedepsurf and the guard stops there. Poison first, so
    ! a routine that ran anyway and happened to land on the last column is not
    ! mistaken for one that correctly did nothing.
    bctfxm = -1. ; bctfxp = -1. ; bctfym = -1. ; bctfyp = -1. ; bctfz = -1.
    timee = t3
    call timedepsurf
    all_passed = surf_frozen('at last table time') .and. all_passed

    timee = t3 + 1000.
    call timedepsurf
    all_passed = surf_frozen('past last table time') .and. all_passed

    ! ---- nudging profiles ---------------------------------------------------
    ltimedepnudge = .true.
    ntimedepnudge = nt
    allocate(timenudge(nt))
    allocate(thlproft(kb:ke+kh,nt), qtproft(kb:ke+kh,nt), &
             uproft  (kb:ke+kh,nt), vproft (kb:ke+kh,nt))

    timenudge = [ t1, t2, t3 ]
    call fill_tables

    timee = t1
    call timedepnudge
    all_passed = prof_exact('t = t1', 1) .and. all_passed

    timee = t1 + 0.25*(t2 - t1)
    call timedepnudge
    all_passed = prof_interp('inside first interval', 1, 0.25) .and. all_passed

    timee = t2
    call timedepnudge
    all_passed = prof_exact('t = t2', 2) .and. all_passed

    ! Not a round fraction and not the same one as above, so the second
    ! interval is genuinely a second case.
    fac = 0.75
    timee = t2 + fac*(t3 - t2)
    call timedepnudge
    all_passed = prof_interp('inside second interval', 2, fac) .and. all_passed

    ! Frozen from the last table time onwards, as for the fluxes, and poisoned
    ! first for the same reason.
    thlprof = -1. ; qtprof = -2. ; uprof = -3. ; vprof = -4.
    timee = t3
    call timedepnudge
    all_passed = prof_frozen('at last table time') .and. all_passed

    timee = t3 + 1000.
    call timedepnudge
    all_passed = prof_frozen('past last table time') .and. all_passed

    ! ---- the switches are respected ----------------------------------------
    ! Both routines return before touching anything when their own switch is
    ! off. A run that sets one of the four is the normal way to use this
    ! module, so a port that interpolates unconditionally would overwrite the
    ! profiles readinitfiles produced with whatever the table held.
    ltimedepsurf  = .false.
    ltimedepnudge = .false.
    bctfxm = 7. ; thlprof = 7. ; uprof = 7.
    timee = t1 + 0.25*(t2 - t1)
    call timedepsurf
    call timedepnudge
    all_passed = same('bctfxm untouched with ltimedepsurf off', bctfxm, 7.) .and. all_passed
    all_passed = same('thlprof untouched with ltimedepnudge off', thlprof(kb), 7.) .and. all_passed
    all_passed = same('uprof untouched with ltimedepnudge off', uprof(kb), 7.) .and. all_passed

    deallocate(timeflux, bctfxmt, bctfxpt, bctfymt, bctfypt, bctfzt)
    deallocate(timenudge, thlproft, qtproft, uproft, vproft)

    if (myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      if (all_passed) then
        write(*, '(A)') 'ALL TESTS PASSED: tests_timedep'
        write(*, '(A,I0,A,I0,A)') '  ', nt, ' table entries over ', ke+kh-kb+1, ' levels'
      else
        write(*, '(A)') 'TESTS FAILED: tests_timedep'
        write(*, '(A)') '  One or more checks did not pass'
      end if
      write(*, '(A)') '================================================'
    end if

    tests_timedep = all_passed

  contains

    !> A table entry that is unique in k, in t and across the four fields.
    !!
    !! Linear in k so that a level read at a fixed index, or off by one, is
    !! wrong everywhere rather than only in the interior; the four fields are
    !! separated by a decade so a swapped pair is unmistakable; and the t
    !! dependence is not the same shape as the k dependence, so a transposed
    !! subscript cannot land on the right number.
    real function entry_at(field, k, t)
      integer, intent(in) :: field, k, t

      entry_at = (10.**field) * (1. + 0.125*(k - kb)) * (1. + 0.5*(t - 1))

    end function entry_at

    subroutine fill_tables
      integer :: kk, tt
      do tt = 1, nt
        do kk = kb, ke+kh
          thlproft(kk,tt) = entry_at(1, kk, tt)
          qtproft (kk,tt) = entry_at(2, kk, tt)
          uproft  (kk,tt) = entry_at(3, kk, tt)
          vproft  (kk,tt) = entry_at(4, kk, tt)
        end do
      end do
    end subroutine fill_tables

    !> Assert a scalar matches to within a few ulps of its own magnitude.
    logical function same(label, got, want)
      character(len=*), intent(in) :: label
      real,             intent(in) :: got, want
      real :: tol

      tol  = 64. * epsilon(1.) * max(1., abs(want))
      same = abs(got - want) <= tol
      if (.not. same) then
        write(*,'(A,I0,3A,ES22.14,A,ES22.14)') 'FAIL on rank ', myid, ': timedep ', &
             trim(label), ' got ', got, ' want ', want
      end if

    end function same

    !> Assert a scalar is bit for bit the value given.
    logical function identical(label, got, want)
      character(len=*), intent(in) :: label
      real,             intent(in) :: got, want

      identical = (got == want)
      if (.not. identical) then
        write(*,'(A,I0,3A,ES22.14,A,ES22.14)') 'FAIL on rank ', myid, ': timedep ', &
             trim(label), ' got ', got, ' want exactly ', want
      end if

    end function identical

    !> Every flux is exactly table entry t.
    logical function surf_exact(label, t)
      character(len=*), intent(in) :: label
      integer,          intent(in) :: t

      surf_exact = identical('bctfxm ' // label, bctfxm, bctfxmt(t))
      surf_exact = identical('bctfxp ' // label, bctfxp, bctfxpt(t)) .and. surf_exact
      surf_exact = identical('bctfym ' // label, bctfym, bctfymt(t)) .and. surf_exact
      surf_exact = identical('bctfyp ' // label, bctfyp, bctfypt(t)) .and. surf_exact
      surf_exact = identical('bctfz '  // label, bctfz,  bctfzt(t))  .and. surf_exact

    end function surf_exact

    !> Nothing was written: every flux still holds the poison put there.
    logical function surf_frozen(label)
      character(len=*), intent(in) :: label

      surf_frozen = same('bctfxm frozen ' // label, bctfxm, -1.)
      surf_frozen = same('bctfxp frozen ' // label, bctfxp, -1.) .and. surf_frozen
      surf_frozen = same('bctfym frozen ' // label, bctfym, -1.) .and. surf_frozen
      surf_frozen = same('bctfyp frozen ' // label, bctfyp, -1.) .and. surf_frozen
      surf_frozen = same('bctfz frozen '  // label, bctfz,  -1.) .and. surf_frozen

    end function surf_frozen

    !> Every flux is entry t blended a fraction f of the way towards t+1.
    logical function surf_interp(label, t, f)
      character(len=*), intent(in) :: label
      integer,          intent(in) :: t
      real,             intent(in) :: f

      surf_interp = same('bctfxm ' // label, bctfxm, bctfxmt(t) + f*(bctfxmt(t+1) - bctfxmt(t)))
      surf_interp = same('bctfxp ' // label, bctfxp, bctfxpt(t) + f*(bctfxpt(t+1) - bctfxpt(t))) .and. surf_interp
      surf_interp = same('bctfym ' // label, bctfym, bctfymt(t) + f*(bctfymt(t+1) - bctfymt(t))) .and. surf_interp
      surf_interp = same('bctfyp ' // label, bctfyp, bctfypt(t) + f*(bctfypt(t+1) - bctfypt(t))) .and. surf_interp
      surf_interp = same('bctfz '  // label, bctfz,  bctfzt(t)  + f*(bctfzt(t+1)  - bctfzt(t)))  .and. surf_interp

    end function surf_interp

    !> Every profile is exactly column t of its own table, level by level.
    logical function prof_exact(label, t)
      character(len=*), intent(in) :: label
      integer,          intent(in) :: t
      integer :: kk

      prof_exact = .true.
      do kk = kb, ke+kh
        prof_exact = identical('thlprof ' // label, thlprof(kk), thlproft(kk,t)) .and. prof_exact
        prof_exact = identical('qtprof '  // label, qtprof(kk),  qtproft(kk,t))  .and. prof_exact
        prof_exact = identical('uprof '   // label, uprof(kk),   uproft(kk,t))   .and. prof_exact
        prof_exact = identical('vprof '   // label, vprof(kk),   vproft(kk,t))   .and. prof_exact
      end do

    end function prof_exact

    !> Nothing was written: every profile still holds the poison put there.
    logical function prof_frozen(label)
      character(len=*), intent(in) :: label
      integer :: kk

      prof_frozen = .true.
      do kk = kb, ke+kh
        prof_frozen = same('thlprof frozen ' // label, thlprof(kk), -1.) .and. prof_frozen
        prof_frozen = same('qtprof frozen '  // label, qtprof(kk),  -2.) .and. prof_frozen
        prof_frozen = same('uprof frozen '   // label, uprof(kk),   -3.) .and. prof_frozen
        prof_frozen = same('vprof frozen '   // label, vprof(kk),   -4.) .and. prof_frozen
      end do

    end function prof_frozen

    !> Every profile is column t blended a fraction f towards column t+1.
    logical function prof_interp(label, t, f)
      character(len=*), intent(in) :: label
      integer,          intent(in) :: t
      real,             intent(in) :: f
      integer :: kk

      prof_interp = .true.
      do kk = kb, ke+kh
        prof_interp = same('thlprof ' // label, thlprof(kk), &
                           thlproft(kk,t) + f*(thlproft(kk,t+1) - thlproft(kk,t))) .and. prof_interp
        prof_interp = same('qtprof ' // label, qtprof(kk), &
                           qtproft(kk,t) + f*(qtproft(kk,t+1) - qtproft(kk,t))) .and. prof_interp
        prof_interp = same('uprof ' // label, uprof(kk), &
                           uproft(kk,t) + f*(uproft(kk,t+1) - uproft(kk,t))) .and. prof_interp
        prof_interp = same('vprof ' // label, vprof(kk), &
                           vproft(kk,t) + f*(vproft(kk,t+1) - vproft(kk,t))) .and. prof_interp
      end do

    end function prof_interp

  end function tests_timedep


end module tests
