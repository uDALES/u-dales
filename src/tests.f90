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
  !!
  !! Exact real comparisons (`==` / `/=` on reals) are deliberate throughout
  !! this file: the assertions here are of the form "bitwise unchanged",
  !! "exactly zero", "exactly one", or "this integer count came back through an
  !! MPI reduction as a real".  Replacing them with tolerance comparisons would
  !! remove the property being tested, so gfortran's -Wcompare-reals is
  !! disabled for this file only, in the top-level CMakeLists.txt.  Nothing
  !! else in src/ is exempt -- a real-equality comparison in solver code is
  !! still flagged, which is the point.
  use decomp_2d
  use modmpi, only : myid, comm3d, mpierr, my_real, avexy_ibm, avey_ibm, sumx_ibm, sumy_ibm

  implicit none
  save
  public :: tests_read_sparse_ijk, tests_2decomp_init_exit, tests_mpi_operators
  public :: tests_nesting_weights, tests_nesting_geometry, tests_nesting_io, &
            tests_nesting_flux, tests_nesting_update, tests_nesting_init

  !> Synthetic solid box used by the nesting IBM subtests (U11-U13, U34).
  !! Defined on GLOBAL indices so that it marks the same physical cells on
  !! every decomposition.
  logical            :: nest_solid_on = .false.
  integer, parameter :: NSOL_I1 = 5, NSOL_I2 = 7
  integer, parameter :: NSOL_J1 = 5, NSOL_J2 = 7
  integer, parameter :: NSOL_K1 = 3, NSOL_K2 = 5

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

  ! ===================================================================
  !  Nesting unit tests, docs/udales-nesting-design.md section 10.1.
  !
  !  Every subtest calls the PRODUCTION routines. Where the quantity of
  !  interest lives in a private variable of modnesting (the zone weights
  !  and the time-interpolated target), it is recovered through the public
  !  interface rather than recomputed: see nest_probe_weight and
  !  nest_target_now below.
  ! ===================================================================

  !> Shape function W(s), the union of face weights, the C1/C2 joints, the
  !! union bounds and the integral identity of design section 1.4(b).
  !! Pure algebra: no grid, no MPI, no I/O. Covers U1-U7.
  logical function tests_nesting_weights()
    use modnesting, only : nest_shape_fn, nest_union

    implicit none

    real, parameter :: limp = 3., lrel = 8.
    logical :: all_passed
    integer :: ish

    call nest_banner('tests_nesting_weights', 'SHAPE FUNCTION AND UNION (U1-U7)')

    all_passed = .true.

    do ish = 1, 2
      if (.not. u1_values(ish))       all_passed = .false.
      if (.not. u2_monotone(ish))     all_passed = .false.
      if (.not. u3_joints(ish))       all_passed = .false.
      if (.not. u7_integral(ish))     all_passed = .false.
    end do
    if (.not. u4_curvature())         all_passed = .false.
    if (.not. u5_union_bounds())      all_passed = .false.
    if (.not. u6_union_degenerate())  all_passed = .false.

    call nest_verdict('tests_nesting_weights', all_passed)
    tests_nesting_weights = nest_all_ranks(all_passed)

  contains

    !> U1: W = 1 for s <= L_imp (and at s = 0), W = 0 for s >= L_imp + L_rel.
    logical function u1_values(ishape)
      integer, intent(in) :: ishape
      integer :: n
      real    :: s
      character(len=32) :: lbl

      write(lbl,'(a,i0,a)') 'U1 values (shape ', ishape, ')'
      u1_values = .true.

      if (nest_shape_fn(0., limp, lrel, ishape) /= 1.) u1_values = .false.
      do n = 0, 100
        s = limp*real(n)/100.
        if (nest_shape_fn(s, limp, lrel, ishape) /= 1.) u1_values = .false.
      end do
      do n = 0, 100
        s = limp + lrel + 5.*real(n)/100.
        if (nest_shape_fn(s, limp, lrel, ishape) /= 0.) u1_values = .false.
      end do
      ! degenerate ramp: a pure guard strip is a step
      if (nest_shape_fn(limp, limp, 0., ishape) /= 1.) u1_values = .false.
      if (nest_shape_fn(limp + 1.e-8, limp, 0., ishape) /= 0.) u1_values = .false.

      call nest_report(trim(lbl), u1_values)
    end function u1_values

    !> U2: 0 <= W <= 1 and W non-increasing over 10^4 samples.
    logical function u2_monotone(ishape)
      integer, intent(in) :: ishape
      integer, parameter  :: nsamp = 10000
      integer :: n
      real    :: s, w, wprev, ds
      character(len=32) :: lbl

      write(lbl,'(a,i0,a)') 'U2 monotone (shape ', ishape, ')'
      u2_monotone = .true.
      ds = (limp + lrel + 2.)/real(nsamp)
      wprev = 1.
      do n = 0, nsamp
        s = -1. + real(n)*ds
        w = nest_shape_fn(s, limp, lrel, ishape)
        if (w < 0. .or. w > 1.) u2_monotone = .false.
        if (w > wprev + 1.e-15) u2_monotone = .false.
        wprev = w
      end do

      call nest_report(trim(lbl), u2_monotone)
    end function u2_monotone

    !> U3: the numerical slope is continuous and ZERO at both joints, i.e.
    !! it vanishes relative to the slope in the middle of the ramp.
    logical function u3_joints(ishape)
      integer, intent(in) :: ishape
      real    :: h, dmid, dlo_in, dhi_in, dlo_out, dhi_out
      character(len=32) :: lbl

      write(lbl,'(a,i0,a)') 'U3 zero slope (shape ', ishape, ')'
      h = 1.e-4*lrel

      dmid    = (nest_shape_fn(limp + 0.5*lrel + h, limp, lrel, ishape) -   &
                 nest_shape_fn(limp + 0.5*lrel - h, limp, lrel, ishape))/(2.*h)
      dlo_in  = (nest_shape_fn(limp + h, limp, lrel, ishape) -              &
                 nest_shape_fn(limp,     limp, lrel, ishape))/h
      dlo_out = (nest_shape_fn(limp,     limp, lrel, ishape) -              &
                 nest_shape_fn(limp - h, limp, lrel, ishape))/h
      dhi_in  = (nest_shape_fn(limp + lrel,     limp, lrel, ishape) -       &
                 nest_shape_fn(limp + lrel - h, limp, lrel, ishape))/h
      dhi_out = (nest_shape_fn(limp + lrel + h, limp, lrel, ishape) -       &
                 nest_shape_fn(limp + lrel,     limp, lrel, ishape))/h

      u3_joints = (abs(dlo_in)  <= 1.e-3*abs(dmid)) .and.                   &
                  (abs(dlo_out) <= 1.e-3*abs(dmid)) .and.                   &
                  (abs(dhi_in)  <= 1.e-3*abs(dmid)) .and.                   &
                  (abs(dhi_out) <= 1.e-3*abs(dmid)) .and.                   &
                  (abs(dlo_in - dlo_out) <= 1.e-3*abs(dmid)) .and.          &
                  (abs(dhi_in - dhi_out) <= 1.e-3*abs(dmid))

      if (.not. u3_joints .and. myid == 0) then
        write(*,'(a,4es12.4,a,es12.4)') '   slopes in/out at both joints ',  &
          dlo_out, dlo_in, dhi_in, dhi_out, ' vs mid-ramp ', dmid
      end if
      call nest_report(trim(lbl), u3_joints)
    end function u3_joints

    !> Second difference of W at s, step h.
    real function d2w(s, ishape, h)
      real,    intent(in) :: s, h
      integer, intent(in) :: ishape
      d2w = (nest_shape_fn(s + h, limp, lrel, ishape)                       &
             - 2.*nest_shape_fn(s, limp, lrel, ishape)                      &
             + nest_shape_fn(s - h, limp, lrel, ishape))/(h*h)
    end function d2w

    !> U4: the quintic is C2 at the joints; the raised cosine deliberately
    !! is not, and the test documents that difference rather than hiding it.
    logical function u4_curvature()
      real :: h, c1lo, c2lo, c2hi, cmid

      h = 1.e-4*lrel
      c1lo = d2w(limp, 1, h)
      c2lo = d2w(limp, 2, h)
      c2hi = d2w(limp + lrel, 2, h)
      cmid = d2w(limp + 0.25*lrel, 2, h)

      u4_curvature = (abs(c2lo) <= 1.e-3*abs(cmid)) .and.                   &
                     (abs(c2hi) <= 1.e-3*abs(cmid)) .and.                   &
                     (abs(c1lo) >  1.e-2*abs(cmid))

      if (myid == 0) then
        write(*,'(a,3es12.4)') '   d2W/ds2 at joints: quintic lo/hi, cosine lo = ', &
          c2lo, c2hi, c1lo
      end if
      call nest_report('U4 quintic C2 (cosine only C1)', u4_curvature)
    end function u4_curvature

    !> U5: the bounded union is <= 1, >= max face weight, equals 1 exactly
    !! when a face demands 1, and is invariant under face permutation.
    logical function u5_union_bounds()
      real, parameter :: vals(5) = (/ 0., 0.2, 0.5, 0.7, 1. /)
      integer :: i1, i2, i3, i4, p
      real    :: w(4), wp(4), u, up_, wmax
      logical :: anyone

      u5_union_bounds = .true.

      do i1 = 1, 5
        do i2 = 1, 5
          do i3 = 1, 5
            do i4 = 1, 5
              w = (/ vals(i1), vals(i2), vals(i3), vals(i4) /)
              u = nest_union(w, 4)
              wmax = maxval(w)
              anyone = any(w >= 1.)
              if (u > 1. .or. u < 0.) u5_union_bounds = .false.
              if (u < wmax - 1.e-15) u5_union_bounds = .false.
              if (anyone .and. u /= 1.) u5_union_bounds = .false.
              if ((.not. anyone) .and. u >= 1.) u5_union_bounds = .false.
              do p = 1, 3
                wp = cshift(w, p)
                up_ = nest_union(wp, 4)
                if (abs(up_ - u) > 1.e-15) u5_union_bounds = .false.
              end do
            end do
          end do
        end do
      end do

      call nest_report('U5 union bounds and symmetry', u5_union_bounds)
    end function u5_union_bounds

    !> U6: with a single active face the union reduces to that face weight.
    !! 1 - (1 - a) is not bitwise a in binary floating point, so the bound is
    !! one unit in the last place rather than exact equality.
    logical function u6_union_degenerate()
      integer :: n, m
      real    :: a, w(4), u, dmax

      u6_union_degenerate = .true.
      dmax = 0.
      do m = 1, 4
        do n = 0, 1000
          a = real(n)/1000.
          w = 0.
          w(m) = a
          u = nest_union(w, 4)
          dmax = max(dmax, abs(u - a))
          if (abs(u - a) > epsilon(1.)) u6_union_degenerate = .false.
        end do
      end do

      if (myid == 0) write(*,'(a,es12.4,a,es12.4)')                        &
        '   max |W_union - W_f| = ', dmax, ', eps = ', epsilon(1.)
      call nest_report('U6 union degeneracy', u6_union_degenerate)
    end function u6_union_degenerate

    !> U7: int W ds = L_imp + L_rel/2, the optical depth of design 1.4(b).
    !! The two panels meet exactly at the joint, so the only error is the
    !! Simpson truncation of a smooth integrand.
    logical function u7_integral(ishape)
      integer, intent(in) :: ishape
      integer, parameter  :: np = 2000
      integer :: n
      real    :: total, h, s, wgt
      character(len=32) :: lbl

      write(lbl,'(a,i0,a)') 'U7 integral (shape ', ishape, ')'

      total = limp
      h = lrel/real(np)
      do n = 0, np
        s = limp + real(n)*h
        if (n == 0 .or. n == np) then
          wgt = 1.
        else if (mod(n,2) == 1) then
          wgt = 4.
        else
          wgt = 2.
        end if
        total = total + wgt*nest_shape_fn(s, limp, lrel, ishape)*h/3.
      end do

      u7_integral = abs(total - (limp + 0.5*lrel)) <= 1.e-10*(limp + lrel)
      if (myid == 0) write(*,'(a,es22.14,a,es22.14)') '   int W ds = ', total, &
        ', expected ', limp + 0.5*lrel
      call nest_report(trim(lbl), u7_integral)
    end function u7_integral

  end function tests_nesting_weights


  ! ------------------------------------------------------------------
  !  Shared helpers for the nesting tests
  ! ------------------------------------------------------------------

  !> Test banner, rank 0 only.
  subroutine nest_banner(name, what)
    use modglobal, only : runmode

    character(len=*), intent(in) :: name, what

    if (myid == 0) then
      write(*,'(a)') '================================================'
      write(*,'(a,i8)') 'runmode = ', runmode
      write(*,'(a,a,a)') trim(name), ': ', trim(what)
      write(*,'(a)') '------------------------------------------------'
    end if

  end subroutine nest_banner

  !> One subtest result line, rank 0 only.
  subroutine nest_report(label, ok)
    character(len=*), intent(in) :: label
    logical,          intent(in) :: ok

    if (myid == 0) then
      if (ok) then
        write(*,'(a,a)') ' PASS  ', trim(label)
      else
        write(*,'(a,a)') ' FAIL  ', trim(label)
      end if
    end if

  end subroutine nest_report

  !> Closing verdict line, rank 0 only.
  subroutine nest_verdict(name, ok)
    character(len=*), intent(in) :: name
    logical,          intent(in) :: ok

    if (myid == 0) then
      write(*,'(a)') '------------------------------------------------'
      if (ok) then
        write(*,'(a,a)') 'ALL TESTS PASSED: ', trim(name)
      else
        write(*,'(a,a)') 'TESTS FAILED: ', trim(name)
      end if
      write(*,'(a)') '================================================'
    end if

  end subroutine nest_verdict

  !> .true. only when every rank passed, so all ranks return the same code.
  logical function nest_all_ranks(ok)
    use mpi

    logical, intent(in) :: ok

    integer :: il, ig

    il = 0
    if (.not. ok) il = 1
    call MPI_ALLREDUCE(il, ig, 1, MPI_INTEGER, MPI_SUM, comm3d, mpierr)
    nest_all_ranks = (ig == 0)

  end function nest_all_ranks

  !> The reference analytic velocity field of the Python fixture writer,
  !! udprep.nesting.analytic_field:
  !!    f = sin(a x + b y) cos(c z) (1 + d t) + e x y z
  !! with the coefficients of udprep.nesting.ANALYTIC_COEFFS. ic: 1 = u,
  !! 2 = v, 3 = w. The build promotes literals to double (-r8), so these
  !! are the same doubles the writer used.
  real function nest_analytic(ic, x, y, z, t)
    integer, intent(in) :: ic
    real,    intent(in) :: x, y, z, t

    real :: a, b, c, d, e

    select case (ic)
    case (1)
      a = 0.017; b = 0.011; c = 0.023; d = 0.05;  e =  1.0e-5
    case (2)
      a = 0.013; b = 0.019; c = 0.029; d = 0.03;  e = -7.0e-6
    case default
      a = 0.023; b = 0.007; c = 0.013; d = 0.07;  e =  4.0e-6
    end select

    nest_analytic = sin(a*x + b*y)*cos(c*z)*(1. + d*t) + e*(x*y*z)

  end function nest_analytic

  !> Synthetic solid box, on GLOBAL indices so that it is the same set on
  !! every decomposition. Used by U11, U12, U13 and U34.
  logical function nest_ref_solid(ig, jg, k)
    integer, intent(in) :: ig, jg, k

    nest_ref_solid = (ig >= NSOL_I1 .and. ig <= NSOL_I2) .and.             &
                     (jg >= NSOL_J1 .and. jg <= NSOL_J2) .and.             &
                     (k  >= NSOL_K1 .and. k  <= NSOL_K2)

  end function nest_ref_solid

  !> Solid, or within nwall cells (Chebyshev, three dimensions) of a solid.
  !! This is the independent dilation the erosion of design section 5 item 2
  !! has to reproduce.
  logical function nest_ref_blocked(ig, jg, k, nwall)
    integer, intent(in) :: ig, jg, k, nwall

    integer :: di, dj, dk

    nest_ref_blocked = .false.
    if (.not. nest_solid_on) return

    do dk = -nwall, nwall
      do dj = -nwall, nwall
        do di = -nwall, nwall
          if (nest_ref_solid(ig + di, jg + dj, k + dk)) nest_ref_blocked = .true.
        end do
      end do
    end do

  end function nest_ref_blocked

  !> Reference weight of component ivar at GLOBAL indices, built from an
  !! independent coordinate formula (x = (ig-1) dx for a face, (ig-1/2) dx
  !! for a centre) but from the PRODUCTION nest_shape_fn and nest_union.
  !! Includes the IBM mask, the wall erosion and the ground-plane rule that
  !! modnesting applies when it builds the zone lists.
  real function nest_ref_weight(ivar, ig, jg, k)
    use modglobal,  only : dx, dy, xlen, ylen, kb, libm
    use modnesting, only : nest_shape_fn, nest_union, nest_guardwidth,     &
                           nest_zonewidth, nest_shape, nest_lateral, nest_nwall

    integer, intent(in) :: ivar, ig, jg, k

    real :: x, y, wf(4)

    nest_ref_weight = 0.

    ! w on the ground plane is set by the bottom BC, not by nesting
    if (ivar == 3 .and. k == kb) return

    select case (ivar)
    case (1)
      x = real(ig - 1)*dx;    y = (real(jg) - 0.5)*dy
    case (2)
      x = (real(ig) - 0.5)*dx; y = real(jg - 1)*dy
    case default
      x = (real(ig) - 0.5)*dx; y = (real(jg) - 0.5)*dy
    end select

    wf = 0.
    if (nest_lateral(1)) wf(1) = nest_shape_fn(x, nest_guardwidth, nest_zonewidth, nest_shape)
    if (nest_lateral(2)) wf(2) = nest_shape_fn(xlen - x, nest_guardwidth, nest_zonewidth, nest_shape)
    if (nest_lateral(3)) wf(3) = nest_shape_fn(y, nest_guardwidth, nest_zonewidth, nest_shape)
    if (nest_lateral(4)) wf(4) = nest_shape_fn(ylen - y, nest_guardwidth, nest_zonewidth, nest_shape)

    nest_ref_weight = nest_union(wf, 4)

    if (nest_ref_weight <= 0.) return
    if (.not. nest_solid_on) return
    if (nest_ref_solid(ig, jg, k)) then
      nest_ref_weight = 0.
    else if (libm .and. nest_nwall >= 1) then
      if (nest_ref_blocked(ig, jg, k, nest_nwall)) nest_ref_weight = 0.
    end if

  end function nest_ref_weight

  !> Install (or remove) the synthetic solid box in IIu/IIv/IIw. The box is
  !! defined on global indices, so every rank marks the same physical cells.
  subroutine nest_set_solids(on)
    use modglobal, only : ib, ie, jb, je, kb, ke, kh, libm
    use modfields, only : IIu, IIv, IIw

    logical, intent(in) :: on

    integer :: i, j, k, ig, jg

    nest_solid_on = on
    libm = on

    IIu = 1; IIv = 1; IIw = 1
    if (.not. on) return

    do k = kb, ke + kh
      do j = jb, je
        do i = ib, ie
          ig = i + zstart(1) - 1
          jg = j + zstart(2) - 1
          if (nest_ref_solid(ig, jg, k)) then
            IIu(i,j,k) = 0
            IIv(i,j,k) = 0
            IIw(i,j,k) = 0
          end if
        end do
      end do
    end do

  end subroutine nest_set_solids

  !> Close any open parent file and re-open the given one with the given
  !! settings. Everything here is a public namelist variable of modnesting,
  !! so the tests reconfigure the production initialisation rather than
  !! duplicating it.
  subroutine nest_reinit(fname, t0, tau, timeinterp, nwall, fluxtol, lassert)
    use modglobal,  only : timee
    use modnesting, only : lnesting, nestfile, nest_tau, nest_timeinterp,  &
                           nest_nwall, nest_fluxtol, nest_lfluxassert,     &
                           nest_lfluxcheckall, nest_linitfromparent,       &
                           nesting_init, nesting_finalize

    character(len=*), intent(in) :: fname
    real,             intent(in) :: t0, tau, fluxtol
    integer,          intent(in) :: timeinterp, nwall
    logical,          intent(in) :: lassert

    call nesting_finalize

    lnesting         = .true.
    nestfile         = fname
    nest_tau         = tau
    nest_timeinterp  = timeinterp
    nest_nwall       = nwall
    nest_fluxtol     = fluxtol
    nest_lfluxassert = lassert
    ! pinned, so a subtest exercises the path it names rather than whatever the
    ! namelist happened to select
    nest_lfluxcheckall   = .false.
    nest_linitfromparent = .false.
    timee            = t0

    call nesting_init

  end subroutine nest_reinit

  !> Recover the production zone weight W at every local point.
  !!
  !! nesting_apply computes  up = (tgt - um)(1 - exp(-W dt_s/tau))/dt_s, so
  !! running it twice from um = 0 and um = 1 removes the (unknown) target:
  !!    up(um=0) - up(um=1) = 1 - exp(-W dt_s/tau).
  !! With rk3step = 0 (dt_s = 1) and tau = 1 that is 1 - exp(-W), from which
  !! W follows exactly. Off-zone points are left untouched by nesting_apply,
  !! so they come back as W = 0.
  subroutine nest_probe_weight(wu, wv, ww)
    use modglobal,  only : ib, ie, jb, je, kb, ke, rk3step
    use modfields,  only : um, up, vm, vp, wm, wp
    use modnesting, only : nest_tau, nesting_apply

    real, intent(out) :: wu(ib:ie,jb:je,kb:ke)
    real, intent(out) :: wv(ib:ie,jb:je,kb:ke)
    real, intent(out) :: ww(ib:ie,jb:je,kb:ke)

    real :: tausave

    tausave  = nest_tau
    nest_tau = 1.
    rk3step  = 0

    um = 0.; vm = 0.; wm = 0.
    up = 0.; vp = 0.; wp = 0.
    call nesting_apply
    wu = up(ib:ie,jb:je,kb:ke)
    wv = vp(ib:ie,jb:je,kb:ke)
    ww = wp(ib:ie,jb:je,kb:ke)

    um = 1.; vm = 1.; wm = 1.
    up = 0.; vp = 0.; wp = 0.
    call nesting_apply
    wu = wu - up(ib:ie,jb:je,kb:ke)
    wv = wv - vp(ib:ie,jb:je,kb:ke)
    ww = ww - wp(ib:ie,jb:je,kb:ke)

    ! wu is now 1 - exp(-W); invert it
    where (wu > 0.) wu = -log(1. - wu)
    where (wv > 0.) wv = -log(1. - wv)
    where (ww > 0.) ww = -log(1. - ww)
    where (wu < 0.) wu = 0.
    where (wv < 0.) wv = 0.
    where (ww < 0.) ww = 0.

    nest_tau = tausave

  end subroutine nest_probe_weight

  !> Recover the CURRENT time-interpolated parent target at every zone point,
  !! without advancing the buffer. tau <= 0 makes the update Dirichlet, so
  !! nesting_apply leaves up = (tgt - um)/dt_s; with um = 0 and dt_s = 1 that
  !! is the target itself. Off-zone points come back as 0.
  subroutine nest_target_now(tu, tv, tw)
    use modglobal,  only : ib, ie, jb, je, kb, ke, rk3step
    use modfields,  only : um, up, vm, vp, wm, wp
    use modnesting, only : nest_tau, nesting_apply

    real, intent(out) :: tu(ib:ie,jb:je,kb:ke)
    real, intent(out) :: tv(ib:ie,jb:je,kb:ke)
    real, intent(out) :: tw(ib:ie,jb:je,kb:ke)

    real :: tausave

    tausave  = nest_tau
    nest_tau = 0.
    rk3step  = 0

    um = 0.; vm = 0.; wm = 0.
    up = 0.; vp = 0.; wp = 0.
    call nesting_apply

    tu = up(ib:ie,jb:je,kb:ke)
    tv = vp(ib:ie,jb:je,kb:ke)
    tw = wp(ib:ie,jb:je,kb:ke)

    nest_tau = tausave

  end subroutine nest_target_now

  !> Advance the production buffer to time t and recover the target there.
  subroutine nest_target_at(t, tu, tv, tw)
    use modglobal,  only : ib, ie, jb, je, kb, ke, timee
    use modnesting, only : nesting_update_target

    real, intent(in)  :: t
    real, intent(out) :: tu(ib:ie,jb:je,kb:ke)
    real, intent(out) :: tv(ib:ie,jb:je,kb:ke)
    real, intent(out) :: tw(ib:ie,jb:je,kb:ke)

    timee = t
    call nesting_update_target
    call nest_target_now(tu, tv, tw)

  end subroutine nest_target_at

  !> Physical coordinates of component ivar at LOCAL indices, from the same
  !! independent formula as nest_ref_weight (used to predict the target).
  subroutine nest_ref_coord(ivar, i, j, k, x, y, z)
    use modglobal, only : dx, dy, zf, zh

    integer, intent(in)  :: ivar, i, j, k
    real,    intent(out) :: x, y, z

    integer :: ig, jg

    ig = i + zstart(1) - 1
    jg = j + zstart(2) - 1

    select case (ivar)
    case (1)
      x = real(ig - 1)*dx;     y = (real(jg) - 0.5)*dy; z = zf(k)
    case (2)
      x = (real(ig) - 0.5)*dx; y = real(jg - 1)*dy;     z = zf(k)
    case default
      x = (real(ig) - 0.5)*dx; y = (real(jg) - 0.5)*dy; z = zh(k)
    end select

  end subroutine nest_ref_coord


  !> Stagger coordinates, zone membership and weights across decompositions,
  !! IBM masking, wall erosion, the building-free rule and the width report.
  !! Covers U8-U14. Run on 1x1, 2x1, 1x2 and 2x2.
  logical function tests_nesting_geometry()
    use mpi
    use modglobal,  only : ib, ie, jb, je, kb, ke, itot, jtot, ktot,       &
                           dx, dy, dzf, zh, xlen, ylen, cexpnr
    use modfields,  only : initfields
    use modibm,     only : createmasks
    use modnesting, only : nest_stagger_coord, nest_guardwidth,            &
                           nest_zonewidth, nest_lparentgeom, nest_nwall,   &
                           nesting_init, nestfile, lnesting

    implicit none

    logical :: all_passed
    real    :: gsave, zsave
    real, allocatable :: wu(:,:,:), wv(:,:,:), ww(:,:,:)

    call nest_banner('tests_nesting_geometry',                             &
                     'STAGGER, ZONE, IBM MASK, EROSION (U8-U14)')

    call initfields
    call createmasks

    allocate(wu(ib:ie,jb:je,kb:ke), wv(ib:ie,jb:je,kb:ke), ww(ib:ie,jb:je,kb:ke))

    ! The namelist selects which half of U13 this invocation runs: the
    ! error path aborts inside nesting_init and can only be observed from
    ! outside the process, so the driver runs it as a separate case.
    if (.not. nest_lparentgeom) then
      tests_nesting_geometry = u13_must_abort()
      deallocate(wu, wv, ww)
      return
    end if

    all_passed = .true.

    if (.not. u8_stagger()) all_passed = .false.

    call nest_set_solids(.false.)
    call nest_reinit('nesting_analytic.'//cexpnr//'.nc', 0., 4., 2, 1, 1.e30, .false.)
    call nest_probe_weight(wu, wv, ww)
    if (.not. u9_membership(wu, wv, ww)) all_passed = .false.
    if (.not. u10_weightsum(wu, wv, ww)) all_passed = .false.

    ! ---- IBM: a solid box inside the zone, no erosion (U11) ----
    call nest_set_solids(.true.)
    call nest_reinit('nesting_analytic.'//cexpnr//'.nc', 0., 4., 2, 0, 1.e30, .false.)
    call nest_probe_weight(wu, wv, ww)
    if (.not. u11_mask(wu, wv, ww)) all_passed = .false.

    ! ---- wall erosion (U12) ----
    call nest_reinit('nesting_analytic.'//cexpnr//'.nc', 0., 4., 2, 1, 1.e30, .false.)
    call nest_probe_weight(wu, wv, ww)
    if (.not. u12_erosion(wu, wv, ww, 1)) all_passed = .false.

    call nest_reinit('nesting_analytic.'//cexpnr//'.nc', 0., 4., 2, 2, 1.e30, .false.)
    call nest_probe_weight(wu, wv, ww)
    if (.not. u12_erosion(wu, wv, ww, 2)) all_passed = .false.

    ! U13, warning half: solids inside the zone are tolerated when the
    ! parent is declared to resolve the child geometry. Reaching this line
    ! at all means the three re-inits above warned instead of aborting.
    call nest_report('U13 building-free rule warns with nest_lparentgeom', .true.)

    call nest_set_solids(.false.)

    ! ---- width reporting and the two width warnings (U14) ----
    gsave = nest_guardwidth
    zsave = nest_zonewidth
    if (.not. u14_width()) all_passed = .false.
    nest_guardwidth = gsave
    nest_zonewidth  = zsave

    deallocate(wu, wv, ww)

    call nest_verdict('tests_nesting_geometry', all_passed)
    tests_nesting_geometry = nest_all_ranks(all_passed)

  contains

    !> U8: nest_stagger_coord against an independent global formula. This is
    !! the half-cell-error test: u sits on xh, v on yh, w on zh, and the
    !! local-to-global shift is zstart - 1.
    logical function u8_stagger()
      integer :: ivar, i, j, k, ig, jg
      real    :: x, y, z, xr, yr, zr, dzu, dmax, tol

      dzu = zh(ktot + 1)/real(ktot)
      dmax = 0.
      tol  = 1.e-12*max(xlen, ylen, zh(ktot + 1))

      ! the independent z formula assumes the fixture's uniform grid
      do k = kb, ke
        dmax = max(dmax, abs(dzf(k) - dzu))
      end do
      if (dmax > tol) then
        call nest_report('U8 stagger (fixture grid is not uniform in z)', .false.)
        u8_stagger = .false.
        return
      end if

      dmax = 0.
      do ivar = 1, 3
        do k = kb, ke
          do j = jb, je
            do i = ib, ie
              call nest_stagger_coord(ivar, i, j, k, x, y, z)
              ig = i + zstart(1) - 1
              jg = j + zstart(2) - 1
              select case (ivar)
              case (1)
                xr = real(ig - 1)*dx;     yr = (real(jg) - 0.5)*dy; zr = (real(k) - 0.5)*dzu
              case (2)
                xr = (real(ig) - 0.5)*dx; yr = real(jg - 1)*dy;     zr = (real(k) - 0.5)*dzu
              case default
                xr = (real(ig) - 0.5)*dx; yr = (real(jg) - 0.5)*dy; zr = real(k - 1)*dzu
              end select
              dmax = max(dmax, abs(x - xr), abs(y - yr), abs(z - zr))
            end do
          end do
        end do
      end do

      ! the top w face, which lives on zh(ke+1)
      do j = jb, je
        do i = ib, ie
          call nest_stagger_coord(3, i, j, ke + 1, x, y, z)
          dmax = max(dmax, abs(z - real(ke)*dzu))
        end do
      end do

      dmax = nest_maxall(dmax)
      u8_stagger = dmax <= tol
      if (myid == 0) write(*,'(a,es12.4)') '   max stagger coordinate error = ', dmax
      call nest_report('U8 stagger coordinates', u8_stagger)
    end function u8_stagger

    !> U9: the set of points carrying W > 0 is exactly the reference set,
    !! and its MPI-reduced size equals the serial global count. Because the
    !! reference is a function of GLOBAL indices only, agreeing with it on
    !! any layout means the layouts agree with each other.
    logical function u9_membership(wu, wv, ww)
      real, intent(in) :: wu(ib:ie,jb:je,kb:ke), wv(ib:ie,jb:je,kb:ke)
      real, intent(in) :: ww(ib:ie,jb:je,kb:ke)

      integer :: ivar, i, j, k, ig, jg, nbad
      real    :: wref, cl(3), cg(3), cref(3), gbad

      nbad = 0
      cl = 0.
      do ivar = 1, 3
        do k = kb, ke
          do j = jb, je
            do i = ib, ie
              ig = i + zstart(1) - 1
              jg = j + zstart(2) - 1
              wref = nest_ref_weight(ivar, ig, jg, k)
              if (wref > 0.) cl(ivar) = cl(ivar) + 1.
              if ((wref > 0.) .neqv. (nest_pick(ivar, wu, wv, ww, i, j, k) > 0.)) then
                nbad = nbad + 1
                if (nbad <= 5 .and. myid == 0) write(*,'(a,4i5,2es12.4)')  &
                  '   membership mismatch ivar,ig,jg,k,got,ref ', ivar, ig, jg, k, &
                  nest_pick(ivar, wu, wv, ww, i, j, k), wref
              end if
            end do
          end do
        end do
      end do

      call MPI_ALLREDUCE(cl, cg, 3, MY_REAL, MPI_SUM, comm3d, mpierr)

      cref = 0.
      do ivar = 1, 3
        do k = kb, ke
          do jg = 1, jtot
            do ig = 1, itot
              if (nest_ref_weight(ivar, ig, jg, k) > 0.) cref(ivar) = cref(ivar) + 1.
            end do
          end do
        end do
      end do

      gbad = nest_sumall(real(nbad))
      u9_membership = (gbad == 0.) .and. all(abs(cg - cref) < 0.5)

      if (myid == 0) write(*,'(a,3i8,a,3i8)') '   zone points u/v/w ',      &
        nint(cg(1)), nint(cg(2)), nint(cg(3)), ' reference ',               &
        nint(cref(1)), nint(cref(2)), nint(cref(3))
      call nest_report('U9 zone membership', u9_membership)
    end function u9_membership

    !> U10: the MPI-reduced sum of W equals the serial global sum, so the
    !! total forcing the zone applies does not depend on the decomposition.
    logical function u10_weightsum(wu, wv, ww)
      real, intent(in) :: wu(ib:ie,jb:je,kb:ke), wv(ib:ie,jb:je,kb:ke)
      real, intent(in) :: ww(ib:ie,jb:je,kb:ke)

      integer :: ivar, k, ig, jg
      real    :: sl(3), sg(3), sref(3), dmax

      sl(1) = sum(wu); sl(2) = sum(wv); sl(3) = sum(ww)
      call MPI_ALLREDUCE(sl, sg, 3, MY_REAL, MPI_SUM, comm3d, mpierr)

      sref = 0.
      do ivar = 1, 3
        do k = kb, ke
          do jg = 1, jtot
            do ig = 1, itot
              sref(ivar) = sref(ivar) + nest_ref_weight(ivar, ig, jg, k)
            end do
          end do
        end do
      end do

      dmax = maxval(abs(sg - sref)/max(sref, 1.))
      u10_weightsum = dmax <= 1.e-12

      if (myid == 0) write(*,'(a,3es22.14)') '   sum W  u/v/w = ', sg
      if (myid == 0) write(*,'(a,3es22.14)') '   reference    = ', sref
      call nest_report('U10 weight sum invariance', u10_weightsum)
    end function u10_weightsum

    !> U11: no solid point carries W > 0, and no fluid point lost its weight.
    logical function u11_mask(wu, wv, ww)
      real, intent(in) :: wu(ib:ie,jb:je,kb:ke), wv(ib:ie,jb:je,kb:ke)
      real, intent(in) :: ww(ib:ie,jb:je,kb:ke)

      integer :: ivar, i, j, k, ig, jg, nsol, nbad
      real    :: wref, wgot, gbad, gsol

      nsol = 0
      nbad = 0
      do ivar = 1, 3
        do k = kb, ke
          do j = jb, je
            do i = ib, ie
              ig = i + zstart(1) - 1
              jg = j + zstart(2) - 1
              wref = nest_ref_weight(ivar, ig, jg, k)
              wgot = nest_pick(ivar, wu, wv, ww, i, j, k)
              if (nest_ref_solid(ig, jg, k)) then
                nsol = nsol + 1
                if (wgot /= 0.) nbad = nbad + 1
              else if (abs(wgot - wref) > 1.e-12) then
                nbad = nbad + 1
              end if
            end do
          end do
        end do
      end do

      gbad = nest_sumall(real(nbad))
      gsol = nest_sumall(real(nsol))
      u11_mask = (gbad == 0.) .and. (gsol > 0.)
      if (myid == 0) write(*,'(a,i0,a,i0)') '   solid points visited ',     &
        nint(gsol), ', mismatches ', nint(gbad)
      call nest_report('U11 IBM masking', u11_mask)
    end function u11_mask

    !> U12: no point within nest_nwall cells of a solid carries W > 0, and the
    !! eroded set matches an independently computed Chebyshev dilation.
    logical function u12_erosion(wu, wv, ww, nwall)
      real,    intent(in) :: wu(ib:ie,jb:je,kb:ke), wv(ib:ie,jb:je,kb:ke)
      real,    intent(in) :: ww(ib:ie,jb:je,kb:ke)
      integer, intent(in) :: nwall

      integer :: ivar, i, j, k, ig, jg, nbad, nero
      real    :: wref, wgot, gbad, gero
      character(len=40) :: lbl

      write(lbl,'(a,i0,a)') 'U12 wall erosion (nest_nwall = ', nwall, ')'

      nbad = 0
      nero = 0
      do ivar = 1, 3
        do k = kb, ke
          do j = jb, je
            do i = ib, ie
              ig = i + zstart(1) - 1
              jg = j + zstart(2) - 1
              wref = nest_ref_weight(ivar, ig, jg, k)
              wgot = nest_pick(ivar, wu, wv, ww, i, j, k)
              if (nest_ref_blocked(ig, jg, k, nwall)) then
                nero = nero + 1
                if (wgot /= 0.) then
                  nbad = nbad + 1
                  if (nbad <= 5 .and. myid == 0) write(*,'(a,4i5,es12.4)') &
                    '   eroded point still forced ivar,ig,jg,k,W ', ivar, ig, jg, k, wgot
                end if
              else if (abs(wgot - wref) > 1.e-12) then
                nbad = nbad + 1
                if (nbad <= 5 .and. myid == 0) write(*,'(a,4i5,2es12.4)')  &
                  '   weight mismatch ivar,ig,jg,k,got,ref ', ivar, ig, jg, k, wgot, wref
              end if
            end do
          end do
        end do
      end do

      gbad = nest_sumall(real(nbad))
      gero = nest_sumall(real(nero))
      u12_erosion = (gbad == 0.) .and. (gero > 0.)
      if (myid == 0) write(*,'(a,i0,a,i0)') '   dilated points ',           &
        nint(gero), ', mismatches ', nint(gbad)
      call nest_report(trim(lbl), u12_erosion)
    end function u12_erosion

    !> U13, error half: with solid points inside the zone and
    !! nest_lparentgeom = .false., nesting_init must abort. Returning from it
    !! is the failure.
    logical function u13_must_abort()
      call nest_set_solids(.true.)
      lnesting = .true.
      nestfile = 'nesting_analytic.'//cexpnr//'.nc'
      nest_nwall = 1
      if (myid == 0) write(*,'(a)')                                        &
        ' U13: nesting_init must abort (solid points in the zone, nest_lparentgeom = F)'
      call nesting_init
      call nest_report('U13 building-free enforcement (init did NOT abort)', .false.)
      call nest_verdict('tests_nesting_geometry', .false.)
      u13_must_abort = .false.
    end function u13_must_abort

    !> U14: the cell-equivalent width reported from a width in metres, and
    !! the two width warnings. The warning text itself is checked by the
    !! Python driver, which greps this run's output.
    logical function u14_width()
      real :: ltot, ncell

      nest_guardwidth = 2.
      nest_zonewidth  = 3.5
      ltot  = nest_guardwidth + nest_zonewidth
      ncell = ltot/dx

      u14_width = (abs(ncell - 5.5) <= 1.e-12) .and. (ncell < 6.) .and.     &
                  (ltot > 0.15*xlen) .and. (abs(dx - xlen/real(itot)) <= 1.e-12)

      if (myid == 0) write(*,'(a,f8.3,a,f8.3,a)') '   L_imp + L_rel = ',    &
        ltot, ' m = ', ncell, ' cells; expecting both width warnings below'
      call nest_reinit('nesting_analytic.'//cexpnr//'.nc', 0., 4., 2, 1, 1.e30, .false.)

      call nest_report('U14 width reporting', u14_width)
    end function u14_width

  end function tests_nesting_geometry

  !> Pick component ivar out of the three probe arrays.
  real function nest_pick(ivar, wu, wv, ww, i, j, k)
    use modglobal, only : ib, ie, jb, je, kb, ke

    integer, intent(in) :: ivar, i, j, k
    real,    intent(in) :: wu(ib:ie,jb:je,kb:ke), wv(ib:ie,jb:je,kb:ke)
    real,    intent(in) :: ww(ib:ie,jb:je,kb:ke)

    select case (ivar)
    case (1)
      nest_pick = wu(i,j,k)
    case (2)
      nest_pick = wv(i,j,k)
    case default
      nest_pick = ww(i,j,k)
    end select

  end function nest_pick

  !> MPI sum of a scalar over comm3d.
  real function nest_sumall(v)
    use mpi

    real, intent(in) :: v

    real :: g

    call MPI_ALLREDUCE(v, g, 1, MY_REAL, MPI_SUM, comm3d, mpierr)
    nest_sumall = g

  end function nest_sumall

  !> MPI maximum of a scalar over comm3d.
  real function nest_maxall(v)
    use mpi

    real, intent(in) :: v

    real :: g

    call MPI_ALLREDUCE(v, g, 1, MY_REAL, MPI_MAX, comm3d, mpierr)
    nest_maxall = g

  end function nest_maxall


  !> Maximum, over all local points and all three components, of the
  !! difference between a probed target and the analytic field of the
  !! fixture (zero off the zone). MPI-reduced.
  real function nest_target_error(t, tu, tv, tw)
    use modglobal, only : ib, ie, jb, je, kb, ke

    real, intent(in) :: t
    real, intent(in) :: tu(ib:ie,jb:je,kb:ke), tv(ib:ie,jb:je,kb:ke)
    real, intent(in) :: tw(ib:ie,jb:je,kb:ke)

    integer :: ivar, i, j, k, ig, jg
    real    :: x, y, z, ex, dmax

    dmax = 0.
    do ivar = 1, 3
      do k = kb, ke
        do j = jb, je
          do i = ib, ie
            ig = i + zstart(1) - 1
            jg = j + zstart(2) - 1
            if (nest_ref_weight(ivar, ig, jg, k) > 0.) then
              call nest_ref_coord(ivar, i, j, k, x, y, z)
              ex = nest_analytic(ivar, x, y, z, t)
            else
              ex = 0.
            end if
            dmax = max(dmax, abs(nest_pick(ivar, tu, tv, tw, i, j, k) - ex))
          end do
        end do
      end do
    end do

    nest_target_error = nest_maxall(dmax)

  end function nest_target_error

  !> Number of elements of two probed target sets that are not bitwise equal.
  real function nest_target_ndiff(au, av, aw, bu, bv, bw)
    use modglobal, only : ib, ie, jb, je, kb, ke

    real, intent(in) :: au(ib:ie,jb:je,kb:ke), av(ib:ie,jb:je,kb:ke), aw(ib:ie,jb:je,kb:ke)
    real, intent(in) :: bu(ib:ie,jb:je,kb:ke), bv(ib:ie,jb:je,kb:ke), bw(ib:ie,jb:je,kb:ke)

    nest_target_ndiff = nest_sumall(real(count(au /= bu) + count(av /= bv) + count(aw /= bw)))

  end function nest_target_ndiff


  !> Reader, writer contract and time buffer. Covers U15-U22.
  !! U15 and U16 are the only tests that exercise the Python writer and the
  !! Fortran reader end to end, so they read the file through modnestingio
  !! directly and compare against udprep.nesting.analytic_field evaluated at
  !! each variable's own stagger.
  logical function tests_nesting_io()
    use mpi
    use modglobal,    only : ib, ie, jb, je, kb, ke, itot, jtot, ktot,      &
                             xf, xh, yf, yh, zf, zh, cexpnr
    use modfields,    only : initfields
    use modibm,       only : createmasks
    use modnestingio, only : nestio_open, nestio_validate, nestio_read,     &
                             nestio_close, nestio_hdr
    use modnesting,   only : nestfile, nesting_finalize

    implicit none

    character(len=5), parameter :: fnam(4) = (/ 'west ', 'east ', 'south', 'north' /)
    character(len=1), parameter :: cnam(3) = (/ 'u', 'v', 'w' /)

    logical :: all_passed
    integer :: ierr, c, f, nzone, ntime
    real, allocatable :: au(:,:,:), av(:,:,:), aw(:,:,:)
    real, allocatable :: bu(:,:,:), bv(:,:,:), bw(:,:,:)

    call nest_banner('tests_nesting_io', 'READER, WRITER CONTRACT, TIME BUFFER (U15-U22)')

    call initfields
    call createmasks
    call nest_set_solids(.false.)

    ! The header-validation cases (U22) abort inside nestio_validate, so the
    ! driver runs each corrupted file as its own invocation and selects it
    ! through the nestfile namelist entry.
    if (index(nestfile, 'bad_') > 0) then
      tests_nesting_io = u22_header()
      return
    end if

    all_passed = .true.

    ! ---- U15/U16: straight through modnestingio, no scheme involved ----
    call nestio_open('nesting_analytic.'//cexpnr//'.nc', ierr)
    if (ierr /= 0) then
      call nest_report('U15 open nesting_analytic', .false.)
      call nest_verdict('tests_nesting_io', .false.)
      tests_nesting_io = .false.
      return
    end if
    call nestio_validate()
    nzone = nestio_hdr%nzone
    ntime = nestio_hdr%ntime

    do f = 1, 4
      do c = 1, 3
        if (.not. u15_u16_slab(c, f)) all_passed = .false.
      end do
    end do
    call nestio_close()

    ! ---- the scheme-level tests need the zone machinery ----
    allocate(au(ib:ie,jb:je,kb:ke), av(ib:ie,jb:je,kb:ke), aw(ib:ie,jb:je,kb:ke))
    allocate(bu(ib:ie,jb:je,kb:ke), bv(ib:ie,jb:je,kb:ke), bw(ib:ie,jb:je,kb:ke))

    if (.not. u17_time_exact(1)) all_passed = .false.
    if (.not. u17_time_exact(2)) all_passed = .false.
    if (.not. u18_hermite_c1())  all_passed = .false.
    if (.not. u44_interp_linear(1)) all_passed = .false.
    if (.not. u44_interp_linear(2)) all_passed = .false.
    if (.not. u19_u21_buffer())  all_passed = .false.
    if (.not. u20_restart())     all_passed = .false.

    deallocate(au, av, aw, bu, bv, bw)
    call nesting_finalize

    call nest_verdict('tests_nesting_io', all_passed)
    tests_nesting_io = nest_all_ranks(all_passed)

  contains

    !> U15: every stored value equals analytic_field at that variable's own
    !! stagger, for every time level, every face and every component.
    !! U16: this rank's hyperslab is bitwise identical to the same range of a
    !! full read, so gathering the buffers from any decomposition gives the
    !! same bits.
    logical function u15_u16_slab(c, f)
      integer, intent(in) :: c, f

      integer :: n1, n2, n3, base, ax, m, kk, d, it, ds, de, nh, igx, jgy
      real    :: x, y, z, ex, dmax, t
      real, allocatable :: full(:,:,:), part(:,:,:)
      logical :: lface_dim, ok16
      character(len=24) :: vname

      vname = cnam(c)//'_'//trim(fnam(f))

      if (f <= 2) then
        lface_dim = (c == 1)
        ax = 2
        n3 = jtot
        if (c == 2) n3 = jtot + 1
        base = 0
        if (f == 2) base = itot - nzone
      else
        lface_dim = (c == 2)
        ax = 1
        n3 = itot
        if (c == 1) n3 = itot + 1
        base = 0
        if (f == 4) base = jtot - nzone
      end if
      n1 = nzone
      if (lface_dim) n1 = nzone + 1
      n2 = ktot
      if (c == 3) n2 = ktot + 1

      allocate(full(n1,n2,n3))

      ds = min(max(zstart(ax), 1), n3)
      de = min(zend(ax) + 1, n3)
      if (de < ds) de = ds
      nh = de - ds + 1
      allocate(part(n1,n2,nh))

      dmax = 0.
      ok16 = .true.

      do it = 1, ntime
        t = nestio_hdr%time(it)

        call nestio_read(trim(vname), it, 1, n3, full, ierr)
        if (ierr /= 0) then
          call nest_report('U15 read '//trim(vname), .false.)
          deallocate(full, part)
          u15_u16_slab = .false.
          return
        end if

        do d = 1, n3
          do kk = 1, n2
            do m = 1, n1
              if (f <= 2) then
                igx = base + m
                jgy = d
              else
                jgy = base + m
                igx = d
              end if
              select case (c)
              case (1)
                x = xh(igx); y = yf(jgy); z = zf(kk)
              case (2)
                x = xf(igx); y = yh(jgy); z = zf(kk)
              case default
                x = xf(igx); y = yf(jgy); z = zh(kk)
              end select
              ex = nest_analytic(c, x, y, z, t)
              dmax = max(dmax, abs(full(m,kk,d) - ex))
            end do
          end do
        end do

        call nestio_read(trim(vname), it, ds, nh, part, ierr)
        if (ierr /= 0) then
          call nest_report('U16 hyperslab read '//trim(vname), .false.)
          deallocate(full, part)
          u15_u16_slab = .false.
          return
        end if
        if (any(part /= full(:,:,ds:de))) ok16 = .false.
      end do

      dmax = nest_maxall(dmax)
      u15_u16_slab = (dmax <= 1.e-13) .and. nest_all_ranks(ok16)

      if (myid == 0) write(*,'(a,a8,a,es12.4,a,l1)') '   ', trim(vname),    &
        ': max |file - analytic| = ', dmax, ', hyperslab bitwise = ', ok16
      call nest_report('U15/U16 '//trim(vname), u15_u16_slab)

      deallocate(full, part)
    end function u15_u16_slab

    !> U17: the fixture is linear in t, so both interpolants must reproduce
    !! it exactly at any time.
    logical function u17_time_exact(interp)
      integer, intent(in) :: interp

      integer :: n
      real    :: t, dmax, err
      character(len=40) :: lbl

      write(lbl,'(a,i0,a)') 'U17 time interpolation (mode ', interp, ')'
      call nest_reinit('nesting_analytic.'//cexpnr//'.nc', 0., 4., interp, 1, 1.e30, .false.)

      dmax = 0.
      do n = 1, 20
        t = 50.*real(n - 1)/19.
        call nest_target_at(t, au, av, aw)
        err = nest_target_error(t, au, av, aw)
        dmax = max(dmax, err)
      end do

      u17_time_exact = dmax <= 1.e-13
      if (myid == 0) write(*,'(a,es12.4)') '   max |target - analytic| over 20 times = ', dmax
      call nest_report(trim(lbl), u17_time_exact)
    end function u17_time_exact

    !> One-sided d(target)/dt jump at a parent node, for interpolation mode
    !! interp, maximised over all zone points and components.
    real function u18_jump(interp)
      integer, intent(in) :: interp

      real, parameter :: tc = 20., h = 1.e-3
      real    :: dmax
      real, allocatable :: du(:,:,:), dv(:,:,:), dw(:,:,:)

      allocate(du(ib:ie,jb:je,kb:ke), dv(ib:ie,jb:je,kb:ke), dw(ib:ie,jb:je,kb:ke))

      call nest_reinit('nesting_nonlinear.'//cexpnr//'.nc', 0., 4., interp, 1, 1.e30, .false.)

      call nest_target_at(tc - h, au, av, aw)
      call nest_target_at(tc,     bu, bv, bw)
      du = (bu - au)/h
      dv = (bv - av)/h
      dw = (bw - aw)/h

      call nest_target_at(tc + h, au, av, aw)
      du = (au - bu)/h - du
      dv = (av - bv)/h - dv
      dw = (aw - bw)/h - dw

      dmax = max(maxval(abs(du)), maxval(abs(dv)), maxval(abs(dw)))
      u18_jump = nest_maxall(dmax)

      deallocate(du, dv, dw)
    end function u18_jump

    !> U18: across a parent-interval crossing the Hermite target has a
    !! continuous time derivative and the linear one does not. The fixture
    !! for this test is deliberately NON-linear in time; a field linear in t
    !! is reproduced exactly by both schemes and so cannot separate them.
    logical function u18_hermite_c1()
      real :: jlin, jher

      jlin = u18_jump(1)
      jher = u18_jump(2)

      u18_hermite_c1 = (jlin > 1.e-4) .and. (jher <= 1.e-2*jlin)
      if (myid == 0) write(*,'(a,es12.4,a,es12.4)')                         &
        '   d(target)/dt jump at t = 20: linear ', jlin, ', Hermite ', jher
      call nest_report('U18 Hermite C1 across an interval crossing', u18_hermite_c1)
    end function u18_hermite_c1

    !> U44: the time interpolant must be LINEAR IN THE DATA.
    !!
    !! This is the property design section 3.1(2) rests on. The net boundary
    !! flux Phi is a linear functional of the boundary values, and the offline
    !! correction makes it zero at every stored level; so if -- and only if --
    !! the interpolated value is a fixed linear combination of the stored
    !! levels, with coefficients that depend on the times but NOT on the data,
    !! does Phi stay zero between levels. Superposition is exactly that
    !! statement, so test it directly:
    !!
    !!     H(alpha*a + beta*b) == alpha*H(a) + beta*H(b)
    !!
    !! An unlimited Hermite passes. A monotone (Fritsch-Carlson) limiter fails
    !! by O(1), which is how the original implementation of this branch broke:
    !! Phi reached 5.2e-5 at run time and divtot 1.7, and no other test in the
    !! matrix could see it -- U28 checks the linear mode's result rather than
    !! its linearity, and I3/I4/I5 all use time-constant boundary data.
    !!
    !! Deliberately uses NON-uniform level spacing and sign-changing data:
    !! a limiter is only active where the one-sided slopes disagree, so data
    !! that happens to be monotone would let it pass.
    logical function u44_interp_linear(mode)
      use modnesting, only : nest_time_interp
      integer, intent(in) :: mode

      integer, parameter :: NTRY = 64
      real,    parameter :: h1 = 0.7, h2 = 1.9, h3 = 0.4
      real,    parameter :: alpha = 1.7, beta = -0.9

      integer :: n, it
      real    :: a(4), b(4), c(4), th, lhs, rhs, sc, worst

      worst = 0.
      do n = 1, NTRY
         ! Deterministic, spread over sign changes and magnitudes so the
         ! limiter's s1*s2 <= 0 branch is exercised.
         do it = 1, 4
            a(it) = sin(1.7*n + 2.3*it) * (1. + 0.5*cos(0.9*n))
            b(it) = cos(0.6*n - 1.1*it) * (1. + 0.5*sin(1.3*n))
         end do
         c = alpha*a + beta*b

         do it = 0, 10
            th  = 0.1*real(it)
            lhs = nest_time_interp(c(1), c(2), c(3), c(4), h1, h2, h3, th, mode)
            rhs = alpha*nest_time_interp(a(1), a(2), a(3), a(4), h1, h2, h3, th, mode) &
                + beta *nest_time_interp(b(1), b(2), b(3), b(4), h1, h2, h3, th, mode)
            sc  = max(1., abs(lhs), abs(rhs))
            worst = max(worst, abs(lhs - rhs)/sc)
         end do
      end do

      u44_interp_linear = (worst <= 1.e-12)
      if (myid == 0) write(*,'(a,i0,a,es12.4)')                             &
        '   nest_timeinterp = ', mode, ': worst relative superposition '//  &
        'error ', worst
      call nest_report('U44 time interpolant is linear in the data', &
                       u44_interp_linear)
    end function u44_interp_linear

    !> U19: stepping across several parent intervals gives the same target as
    !! a fresh initialisation at the same time (which reads every level
    !! eagerly). U21: the incremental roll path and the full-reload path,
    !! which is what a prefetch would bypass, agree bitwise. modnesting has
    !! no separate prefetch switch; the two code paths in set_interval are
    !! the thing to compare.
    logical function u19_u21_buffer()
      real, parameter :: tq = 47.3
      real    :: t, nd1, nd2
      integer :: n

      ! eager: fresh init at tq reloads all NSLOT slots
      call nest_reinit('nesting_nonlinear.'//cexpnr//'.nc', tq, 4., 2, 1, 1.e30, .false.)
      call nest_target_now(au, av, aw)

      ! incremental: init at 0, walk forward in small steps (roll path)
      call nest_reinit('nesting_nonlinear.'//cexpnr//'.nc', 0., 4., 2, 1, 1.e30, .false.)
      do n = 1, 20
        t = min(tq, 2.5*real(n))
        call nest_target_at(t, bu, bv, bw)
      end do
      call nest_target_at(tq, bu, bv, bw)
      nd1 = nest_target_ndiff(au, av, aw, bu, bv, bw)

      ! one jump of NSLOT intervals, which takes the full-reload branch
      call nest_reinit('nesting_nonlinear.'//cexpnr//'.nc', 0., 4., 2, 1, 1.e30, .false.)
      call nest_target_at(tq, bu, bv, bw)
      nd2 = nest_target_ndiff(au, av, aw, bu, bv, bw)

      u19_u21_buffer = (nd1 == 0.) .and. (nd2 == 0.)
      if (myid == 0) write(*,'(a,i0,a,i0)') '   differing targets: rolled ', &
        nint(nd1), ', reloaded ', nint(nd2)
      call nest_report('U19/U21 buffer roll and full reload', u19_u21_buffer)
    end function u19_u21_buffer

    !> U20: a restart mid-interval reproduces the target of a continuous run
    !! at the same time, bitwise.
    logical function u20_restart()
      use modnesting, only : nesting_restart_write, nesting_restart_read

      real, parameter :: tstar = 33.7
      integer :: iu, n
      real    :: t, nd

      call nest_reinit('nesting_nonlinear.'//cexpnr//'.nc', 0., 4., 2, 1, 1.e30, .false.)
      do n = 1, 20
        t = min(tstar, 2.5*real(n))
        call nest_target_at(t, au, av, aw)
      end do
      call nest_target_at(tstar, au, av, aw)

      open(newunit=iu, form='unformatted', status='scratch')
      call nesting_restart_write(iu)
      rewind(iu)

      call nest_reinit('nesting_nonlinear.'//cexpnr//'.nc', 0., 4., 2, 1, 1.e30, .false.)
      call nesting_restart_read(iu)
      close(iu)

      call nest_target_now(bu, bv, bw)
      nd = nest_target_ndiff(au, av, aw, bu, bv, bw)

      u20_restart = (nd == 0.)
      if (myid == 0) write(*,'(a,i0)') '   differing targets after restart = ', nint(nd)
      call nest_report('U20 restart repositioning', u20_restart)
    end function u20_restart

    !> U22: a header field that disagrees with the run must abort with a
    !! message naming the field. One corrupted file per invocation.
    logical function u22_header()
      integer :: ie2

      if (myid == 0) write(*,'(a,a)') ' U22: nestio_validate must abort on ', trim(nestfile)

      call nestio_open(trim(nestfile), ie2)
      if (ie2 /= 0) then
        call nest_report('U22 open '//trim(nestfile), .false.)
        call nest_verdict('tests_nesting_io', .false.)
        u22_header = .false.
        return
      end if

      ! bad_stagger is now handled by the common path below: nestio_validate
      ! checks the per-variable stagger attribute (spec section 6), so a wrong
      ! tag must abort exactly like a wrong grid.

      call nestio_validate()

      call nest_report('U22 header validation ('//trim(nestfile)//                &
                       '): validate did NOT abort', .false.)
      call nest_verdict('tests_nesting_io', .false.)
      u22_header = .false.
    end function u22_header

  end function tests_nesting_io


  !> Fill an array with a deterministic pattern (assumed shape, so halos are
  !! filled too and an unintended write anywhere shows up as a difference).
  subroutine nest_fill(a, s)
    real, intent(out) :: a(:,:,:)
    real, intent(in)  :: s

    integer :: i, j, k

    do k = 1, size(a,3)
      do j = 1, size(a,2)
        do i = 1, size(a,1)
          a(i,j,k) = s*(0.017*real(i) - 0.011*real(j) + 0.0031*real(k) + 0.25)
        end do
      end do
    end do

  end subroutine nest_fill


  !> The mass-compatibility machinery: Phi against an analytically known
  !! value, its decomposition invariance, the fluid-face mask, a corrected
  !! input file, the assertion firing on an uncorrected one, and the
  !! linearity of Phi in time. Covers U23-U28.
  logical function tests_nesting_flux()
    use mpi
    use modglobal,  only : ib, ie, ih, jb, je, jh, kb, ke, kh, itot, jtot,  &
                           dx, dy, dzf, ibrank, ierank, jbrank, jerank,     &
                           cexpnr, timee
    use modfields,  only : initfields, IIu, IIv, IIw, rhobf, rhobh
    use modibm,     only : createmasks
    use modnesting,   only : nest_flux_residual, nest_flux_split, nesting_bcpup, &
                             nesting_update_target, nestfile, lnesting,          &
                             nest_fluxtol, nest_lfluxassert, nesting_init,       &
                             nesting_finalize
    use modnestingio, only : nestio_hdr, nestio_read

    implicit none

    real, parameter :: rk3c = 0.7

    logical :: all_passed, lmask
    real, allocatable :: pup(:,:,:), pvp(:,:,:), pwp(:,:,:)

    call nest_banner('tests_nesting_flux', 'BOUNDARY FLUX RESIDUAL (U23-U28)')

    call initfields
    call createmasks
    call nest_set_solids(.false.)

    allocate(pup(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
    allocate(pvp(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
    allocate(pwp(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))

    ! U27 is an abort case and gets its own invocation, selected by nestfile.
    if (index(nestfile, 'assertfire') > 0) then
      tests_nesting_flux = u27_assert_fires()
      deallocate(pup, pvp, pwp)
      return
    end if

    all_passed = .true.

    lmask = .false.
    if (.not. u23_phi(lmask)) all_passed = .false.

    lmask = .true.
    call set_boundary_mask(.true.)
    if (.not. u23_phi(lmask)) all_passed = .false.
    call set_boundary_mask(.false.)

    if (.not. u26_corrected()) all_passed = .false.
    if (.not. u28_linear())    all_passed = .false.

    ! design section 10.6 item 5 -- the lid contribution to Phi
    if (.not. u35_lid_split())  all_passed = .false.
    if (.not. u36_closed_lid()) all_passed = .false.

    ! design section 10.6 item 3 -- the stored post-correction residual
    if (.not. u37_stored_residual()) all_passed = .false.
    if (.not. u38_schema1_file())    all_passed = .false.
    if (.not. u39_lying_file())      all_passed = .false.

    deallocate(pup, pvp, pwp)
    call nesting_finalize

    call nest_verdict('tests_nesting_flux', all_passed)
    tests_nesting_flux = nest_all_ranks(all_passed)

  contains

    !> Imposed normal velocity on the six domain-boundary faces, as a
    !! function of GLOBAL indices. iface: 1 west, 2 east, 3 south, 4 north,
    !! 5 bottom, 6 top.
    real function face_val(iface, ga, gb)
      integer, intent(in) :: iface, ga, gb

      select case (iface)
      case (1)
        face_val =  0.31*sin(0.21*real(ga)) + 0.07*real(gb)
      case (2)
        face_val = -0.19*cos(0.13*real(ga)) + 0.05*real(gb)
      case (3)
        face_val =  0.23*sin(0.17*real(ga)) - 0.03*real(gb)
      case (4)
        face_val =  0.11*cos(0.29*real(ga)) + 0.09*real(gb)
      case (5)
        face_val =  0.013*real(ga) - 0.007*real(gb)
      case default
        face_val = -0.005*real(ga) + 0.011*real(gb)
      end select

    end function face_val

    !> Synthetic fluid/solid pattern on the boundary faces (U25), on global
    !! indices. .true. = fluid.
    logical function face_fluid(iface, ga, gb, lmasked)
      integer, intent(in) :: iface, ga, gb
      logical, intent(in) :: lmasked

      face_fluid = .true.
      if (.not. lmasked) return

      select case (iface)
      case (1)
        face_fluid = .not. (mod(ga,4) == 0 .and. gb <= 4)
      case (2)
        face_fluid = .not. (mod(ga,5) == 0)
      case (3)
        face_fluid = .not. (mod(ga,3) == 0 .and. gb >= 12)
      case (4)
        face_fluid = .not. (mod(ga,7) == 0)
      case (5)
        face_fluid = .not. (mod(ga + gb, 6) == 0)
      case default
        face_fluid = .true.
      end select

    end function face_fluid

    !> Apply (or clear) the synthetic boundary-face mask in IIu/IIv/IIw.
    subroutine set_boundary_mask(on)
      logical, intent(in) :: on

      integer :: i, j, k, ig, jg

      IIu = 1; IIv = 1; IIw = 1
      if (.not. on) return

      if (ibrank) then
        do k = kb, ke
          do j = jb, je
            jg = j + zstart(2) - 1
            if (.not. face_fluid(1, jg, k, .true.)) IIu(ib,j,k) = 0
          end do
        end do
      end if
      if (ierank) then
        do k = kb, ke
          do j = jb, je
            jg = j + zstart(2) - 1
            if (.not. face_fluid(2, jg, k, .true.)) IIu(ie+1,j,k) = 0
          end do
        end do
      end if
      if (jbrank) then
        do k = kb, ke
          do i = ib, ie
            ig = i + zstart(1) - 1
            if (.not. face_fluid(3, ig, k, .true.)) IIv(i,jb,k) = 0
          end do
        end do
      end if
      if (jerank) then
        do k = kb, ke
          do i = ib, ie
            ig = i + zstart(1) - 1
            if (.not. face_fluid(4, ig, k, .true.)) IIv(i,je+1,k) = 0
          end do
        end do
      end if
      do j = jb, je
        do i = ib, ie
          ig = i + zstart(1) - 1
          jg = j + zstart(2) - 1
          if (.not. face_fluid(5, ig, jg, .true.)) IIw(i,j,kb) = 0
        end do
      end do

    end subroutine set_boundary_mask

    !> Impose the six-face field of face_val on the predicted velocity.
    !! ltop selects whether the lid carries data or is closed (w* = 0), which is
    !! what a rigid lid does in bcpup.
    subroutine fill_boundary(ltop)
      logical, intent(in) :: ltop

      integer :: i, j, k, ig, jg

      ! interior values must not enter Phi
      call nest_fill(pup, 3.1)
      call nest_fill(pvp, -2.4)
      call nest_fill(pwp, 1.7)

      if (ibrank) then
        do k = kb, ke
          do j = jb - jh, je + jh
            jg = j + zstart(2) - 1
            pup(ib,j,k) = face_val(1, jg, k)/rk3c
          end do
        end do
      end if
      if (ierank) then
        do k = kb, ke
          do j = jb - jh, je + jh
            jg = j + zstart(2) - 1
            pup(ie+1,j,k) = face_val(2, jg, k)/rk3c
          end do
        end do
      end if
      if (jbrank) then
        do k = kb, ke
          do i = ib - ih, ie + ih
            ig = i + zstart(1) - 1
            pvp(i,jb,k) = face_val(3, ig, k)/rk3c
          end do
        end do
      end if
      if (jerank) then
        do k = kb, ke
          do i = ib - ih, ie + ih
            ig = i + zstart(1) - 1
            pvp(i,je+1,k) = face_val(4, ig, k)/rk3c
          end do
        end do
      end if
      do j = jb - jh, je + jh
        do i = ib - ih, ie + ih
          ig = i + zstart(1) - 1
          jg = j + zstart(2) - 1
          pwp(i,j,kb) = face_val(5, ig, jg)/rk3c
          if (ltop) then
            pwp(i,j,ke+1) = face_val(6, ig, jg)/rk3c
          else
            pwp(i,j,ke+1) = 0.
          end if
        end do
      end do

    end subroutine fill_boundary

    !> Reference Phi over global indices. iset: 1 all six faces, 2 the lid
    !! alone, 3 the five closed faces. All three are normalised by the SAME
    !! total fluid boundary area, as nest_flux_split is.
    real function ref_phi(iset, lmasked)
      integer, intent(in) :: iset
      logical, intent(in) :: lmasked

      integer :: ig, jg, k
      real    :: num, area, af

      num  = 0.
      area = 0.
      do k = kb, ke
        af = dy*dzf(k)
        do jg = 1, jtot
          if (face_fluid(1, jg, k, lmasked)) then
            if (iset /= 2) num = num - rhobf(k)*face_val(1, jg, k)*af
            area = area + af
          end if
          if (face_fluid(2, jg, k, lmasked)) then
            if (iset /= 2) num = num + rhobf(k)*face_val(2, jg, k)*af
            area = area + af
          end if
        end do
        af = dx*dzf(k)
        do ig = 1, itot
          if (face_fluid(3, ig, k, lmasked)) then
            if (iset /= 2) num = num - rhobf(k)*face_val(3, ig, k)*af
            area = area + af
          end if
          if (face_fluid(4, ig, k, lmasked)) then
            if (iset /= 2) num = num + rhobf(k)*face_val(4, ig, k)*af
            area = area + af
          end if
        end do
      end do
      af = dx*dy
      do jg = 1, jtot
        do ig = 1, itot
          if (face_fluid(5, ig, jg, lmasked)) then
            if (iset /= 2) num = num - rhobh(kb)*face_val(5, ig, jg)*af
            area = area + af
          end if
          if (face_fluid(6, ig, jg, lmasked)) then
            if (iset /= 3) num = num + rhobh(ke+1)*face_val(6, ig, jg)*af
            area = area + af
          end if
        end do
      end do

      ref_phi = num/area

    end function ref_phi

    !> U23/U24/U25: Phi for an imposed face field with an analytically known
    !! net flux. The reference is summed over GLOBAL indices, so agreeing
    !! with it on any layout is decomposition invariance.
    logical function u23_phi(lmasked)
      logical, intent(in) :: lmasked

      real    :: phi, ref, tol
      character(len=48) :: lbl

      if (lmasked) then
        lbl = 'U25 Phi over fluid faces only (masked)'
      else
        lbl = 'U23/U24 Phi against an analytic net flux'
      end if

      call fill_boundary(.true.)

      phi = nest_flux_residual(pup, pvp, pwp, rk3c)

      ! independent global reference
      ref = ref_phi(1, lmasked)

      tol = 1.e-13*max(abs(ref), 1.e-3)
      u23_phi = abs(phi - ref) <= tol

      if (myid == 0) then
        write(*,'(a,es24.16)') ' NESTFLUX phi      = ', phi
        write(*,'(a,es24.16)') ' NESTFLUX phi_ref  = ', ref
      end if
      call nest_report(trim(lbl), u23_phi)
    end function u23_phi

    !> U26: every stored level of a writer-corrected file is flux balanced
    !! (nesting_init checks all of them and aborts otherwise), and the faces
    !! nesting_bcpup imposes give Phi = 0 at run time as well.
    logical function u26_corrected()
      integer :: n
      real    :: t, phi, pmax

      call nest_reinit('nesting_corrected.'//cexpnr//'.nc', 0., 4., 2, 1, 1.e-10, .true.)

      pmax = 0.
      do n = 1, 7
        t = 50.*real(n - 1)/6.
        timee = t
        call nesting_update_target
        pup = 0.; pvp = 0.; pwp = 0.
        call nesting_bcpup(pup, pvp, pwp, rk3c)
        phi = nest_flux_residual(pup, pvp, pwp, rk3c)
        pmax = max(pmax, abs(phi))
      end do

      u26_corrected = pmax <= 1.e-10
      if (myid == 0) write(*,'(a,es12.4)') '   max |Phi| over 7 times = ', pmax
      call nest_report('U26 corrected input gives Phi = 0', u26_corrected)
    end function u26_corrected

    !> U27: an uncorrected file must trip the flux assertion. Reaching the
    !! line after nesting_init is the failure.
    logical function u27_assert_fires()
      lnesting         = .true.
      nest_fluxtol     = 1.e-10
      nest_lfluxassert = .true.
      timee            = 0.
      if (myid == 0) write(*,'(a,a)')                                       &
        ' U27: nesting_init must abort on the uncorrected file ', trim(nestfile)

      call nesting_init

      call nest_report('U27 flux assertion (init did NOT abort)', .false.)
      call nest_verdict('tests_nesting_flux', .false.)
      u27_assert_fires = .false.
    end function u27_assert_fires


    !> U35: the lid split of design section 10.6 item 5. phi_all must still be
    !! the six-face Phi, phi_lid the top face alone, and the difference the five
    !! closed faces -- each against an independently summed global reference.
    !! This is what lets the flux assertion stay on under BCtopm_pressure, where
    !! the lid flux is a free response and not an imposed datum.
    logical function u35_lid_split()
      real :: phi_all, phi_lid, r_all, r_lid, r_closed, tol

      call fill_boundary(.true.)
      call nest_flux_split(pup, pvp, pwp, rk3c, phi_all, phi_lid)

      r_all    = ref_phi(1, .false.)
      r_lid    = ref_phi(2, .false.)
      r_closed = ref_phi(3, .false.)

      tol = 1.e-13*max(abs(r_all), abs(r_lid), 1.e-3)

      u35_lid_split = abs(phi_all - r_all) <= tol .and.                     &
                      abs(phi_lid - r_lid) <= tol .and.                     &
                      abs((phi_all - phi_lid) - r_closed) <= tol

      if (myid == 0) then
        write(*,'(a,2es22.14)') ' NESTFLUX phi_all, ref = ', phi_all, r_all
        write(*,'(a,2es22.14)') ' NESTFLUX phi_lid, ref = ', phi_lid, r_lid
        write(*,'(a,2es22.14)') ' NESTFLUX closed,  ref = ', phi_all - phi_lid, r_closed
      end if
      call nest_report('U35 Phi splits into lid and closed faces', u35_lid_split)
    end function u35_lid_split

    !> U36: with a rigid lid bcpup sets w* = 0 there, so the lid carries no flux
    !! at all and the asserted quantity is exactly the old six-face Phi. The
    !! split must therefore change nothing for design case A.
    logical function u36_closed_lid()
      real :: phi_all, phi_lid

      call fill_boundary(.false.)
      call nest_flux_split(pup, pvp, pwp, rk3c, phi_all, phi_lid)

      u36_closed_lid = (phi_lid == 0.) .and. (abs(phi_all - ref_phi(3, .false.)) <= 1.e-13)

      if (myid == 0) write(*,'(a,es12.4,a,es12.4)') '   phi_lid = ', phi_lid, &
        ', phi_all = ', phi_all
      call nest_report('U36 a closed lid contributes nothing to Phi', u36_closed_lid)
    end function u36_closed_lid

    !> U37: the schema 2 flux_residual is the residual of the data AS STORED.
    !! Recomputed here from the four boundary-normal slabs, read straight through
    !! modnestingio, which is exactly the quantity the cheap init-time check
    !! trusts (design section 10.6 item 3). Rank 0 reads the whole range, so the
    !! comparison does not depend on the decomposition.
    logical function u37_stored_residual()
      use modglobal, only : cexpnr, ktot

      integer :: it, kk, d, ierr, nz
      real    :: num, dmax, aref
      real, allocatable :: bw(:,:,:), be(:,:,:), bs(:,:,:), bn(:,:,:)
      logical :: ok

      call nest_reinit('nesting_corrected.'//cexpnr//'.nc', 0., 4., 2, 1, 1.e-10, .true.)

      ok = nestio_hdr%has_flux_residual
      dmax = 0.
      aref = 0.

      if (ok .and. myid == 0) then
        nz = nestio_hdr%nzone
        allocate(bw(nz + 1, ktot, jtot), be(nz + 1, ktot, jtot))
        allocate(bs(nz + 1, ktot, itot), bn(nz + 1, ktot, itot))
        do it = 1, nestio_hdr%ntime
          call nestio_read('u_west',  it, 1, jtot, bw, ierr)
          call nestio_read('u_east',  it, 1, jtot, be, ierr)
          call nestio_read('v_south', it, 1, itot, bs, ierr)
          call nestio_read('v_north', it, 1, itot, bn, ierr)
          num = 0.
          do d = 1, jtot
            do kk = 1, ktot
              num = num + (be(nz + 1, kk, d) - bw(1, kk, d))*rhobf(kk)*dy*dzf(kk)
            end do
          end do
          do d = 1, itot
            do kk = 1, ktot
              num = num + (bn(nz + 1, kk, d) - bs(1, kk, d))*rhobf(kk)*dx*dzf(kk)
            end do
          end do
          dmax = max(dmax, abs(num - nestio_hdr%flux_residual(it)))
        end do
        deallocate(bw, be, bs, bn)

        ! the same slab data give the fluid lateral area the file advertises
        aref = 2.*real(jtot)*dy*sum(dzf(kb:ke)) + 2.*real(itot)*dx*sum(dzf(kb:ke))
        ok = ok .and. abs(aref - nestio_hdr%fluid_lateral_area) <= 1.e-12*aref
      end if

      if (myid == 0) then
        ok = ok .and. (dmax <= 1.e-12*max(1., aref))
        write(*,'(a,es12.4)') '   max |stored residual - recomputed| = ', dmax
        write(*,'(a,2es22.14)') '   fluid lateral area file/ref = ', &
          nestio_hdr%fluid_lateral_area, aref
      end if

      u37_stored_residual = nest_all_ranks(ok .or. myid /= 0)
      call nest_report('U37 stored residual matches the stored slabs', u37_stored_residual)
    end function u37_stored_residual

    !> U38: a schema 1 file has no stored residual and must still load and pass,
    !! by the full recompute, with a warning. Backwards compatibility of the file
    !! format is a hard requirement, so it gets its own subtest.
    logical function u38_schema1_file()
      use modglobal, only : cexpnr

      real :: phi

      call nest_reinit('nesting_v1.'//cexpnr//'.nc', 0., 4., 2, 1, 1.e-10, .true.)

      pup = 0.; pvp = 0.; pwp = 0.
      call nesting_bcpup(pup, pvp, pwp, rk3c)
      phi = nest_flux_residual(pup, pvp, pwp, rk3c)

      u38_schema1_file = (nestio_hdr%schema == 1) .and.                     &
                         (.not. nestio_hdr%has_flux_residual) .and.         &
                         (abs(phi) <= 1.e-10)

      if (myid == 0) write(*,'(a,i0,a,l1,a,es12.4)') '   schema = ', nestio_hdr%schema, &
        ', has_flux_residual = ', nestio_hdr%has_flux_residual, ', Phi = ', phi
      call nest_report('U38 a schema 1 file still loads and is checked', u38_schema1_file)
    end function u38_schema1_file

    !> U39: the cheap check believes the writer. A file whose stored residual
    !! claims zero while its data do not is accepted at init -- and caught by the
    !! per-substep assertion in nesting_bcpup instead. The subtest pins that
    !! trade-off rather than leaving it implied; nest_lfluxcheckall is what
    !! closes it, and its abort case is in the driver.
    logical function u39_lying_file()
      use modglobal, only : cexpnr

      real :: phi

      call nest_reinit('assertfire_lying.'//cexpnr//'.nc', 0., 4., 2, 1, 1.e-10, .false.)

      pup = 0.; pvp = 0.; pwp = 0.
      call nesting_bcpup(pup, pvp, pwp, rk3c)
      phi = nest_flux_residual(pup, pvp, pwp, rk3c)

      u39_lying_file = abs(phi) > 1.e-10

      if (myid == 0) write(*,'(a,es12.4)') &
        '   run-time Phi of the file the cheap check accepted = ', phi
      call nest_report('U39 the run-time assertion backs up the cheap check', &
                       u39_lying_file)
    end function u39_lying_file

    !> U28: Phi is a linear functional of the boundary data, so with a target
    !! that is linear in time Phi at the midpoint equals the mean of the
    !! endpoints. This is design section 3.1 item 2, checked in the code.
    logical function u28_linear()
      real :: p1, p2, pm, t1, t2

      call nest_reinit('nesting_analytic.'//cexpnr//'.nc', 0., 4., 2, 1, 1.e30, .false.)

      t1 = 5.
      t2 = 45.
      p1 = phi_at(t1)
      p2 = phi_at(t2)
      pm = phi_at(0.5*(t1 + t2))

      u28_linear = abs(pm - 0.5*(p1 + p2)) <= 1.e-12*max(abs(p1), abs(p2), 1.e-6)
      if (myid == 0) write(*,'(a,3es22.14)') '   Phi(t1), Phi(t2), Phi(tm) = ', p1, p2, pm
      call nest_report('U28 linearity of Phi in time', u28_linear)
    end function u28_linear

    real function phi_at(t)
      real, intent(in) :: t

      timee = t
      call nesting_update_target
      pup = 0.; pvp = 0.; pwp = 0.
      call nesting_bcpup(pup, pvp, pwp, rk3c)
      phi_at = nest_flux_residual(pup, pvp, pwp, rk3c)

    end function phi_at

  end function tests_nesting_flux


  !> The substep-implicit relaxation update of design section 1.2, on the
  !! production zone. Covers U29-U34.
  logical function tests_nesting_update()
    use mpi
    use modglobal,  only : ib, ie, ih, jb, je, jh, kb, ke, kh, dt, rk3step, &
                           cexpnr, timee
    use modfields,  only : initfields, um, up, vm, vp, wm, wp
    use modibm,     only : createmasks
    use modnesting, only : nest_tau, nest_lparentgeom, nesting_apply,       &
                           nesting_finalize

    implicit none

    logical :: all_passed
    real, allocatable :: wu(:,:,:), wv(:,:,:), ww(:,:,:)
    real, allocatable :: su(:,:,:), sv(:,:,:), sw(:,:,:)
    real, allocatable :: mu(:,:,:), mv(:,:,:), mw(:,:,:)

    call nest_banner('tests_nesting_update', 'RELAXATION UPDATE (U29-U34)')

    call initfields
    call createmasks
    call nest_set_solids(.false.)

    allocate(wu(ib:ie,jb:je,kb:ke), wv(ib:ie,jb:je,kb:ke), ww(ib:ie,jb:je,kb:ke))
    allocate(su(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
    allocate(sv(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
    allocate(sw(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
    allocate(mu(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
    allocate(mv(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
    allocate(mw(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))

    call nest_reinit('nesting_analytic.'//cexpnr//'.nc', 0., 4., 2, 1, 1.e30, .false.)
    call nest_probe_weight(wu, wv, ww)

    all_passed = .true.
    if (.not. u29_noop())        all_passed = .false.
    if (.not. u30_dirichlet(0.))     all_passed = .false.
    if (.not. u30_dirichlet(1.e-30)) all_passed = .false.
    if (.not. u31_linear())      all_passed = .false.
    if (.not. u32_fullstep())    all_passed = .false.
    if (.not. u33_stability())   all_passed = .false.
    if (.not. u34_solid())       all_passed = .false.

    deallocate(wu, wv, ww, su, sv, sw, mu, mv, mw)
    call nesting_finalize

    call nest_verdict('tests_nesting_update', all_passed)
    tests_nesting_update = nest_all_ranks(all_passed)

  contains

    !> Set up a deterministic (um, up) state and remember it.
    subroutine set_state(s)
      real, intent(in) :: s

      call nest_fill(um, s)
      call nest_fill(vm, 1.3*s)
      call nest_fill(wm, -0.7*s)
      call nest_fill(up, 0.5*s)
      call nest_fill(vp, -0.9*s)
      call nest_fill(wp, 1.1*s)
      mu = um(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh)
      mv = vm(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh)
      mw = wm(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh)
      su = up
      sv = vp
      sw = wp

    end subroutine set_state

    !> Count local points where up changed but the reference weight is zero,
    !! i.e. where nesting_apply was not entitled to touch anything.
    real function offzone_changes()
      integer :: i, j, k, ig, jg, n

      n = 0
      do k = kb, ke + kh
        do j = jb - jh, je + jh
          do i = ib - ih, ie + ih
            ig = i + zstart(1) - 1
            jg = j + zstart(2) - 1
            if (in_zone(1, i, j, k)) cycle
            if (up(i,j,k) /= su(i,j,k)) n = n + 1
          end do
        end do
      end do
      do k = kb, ke + kh
        do j = jb - jh, je + jh
          do i = ib - ih, ie + ih
            if (in_zone(2, i, j, k)) cycle
            if (vp(i,j,k) /= sv(i,j,k)) n = n + 1
          end do
        end do
      end do
      do k = kb, ke + kh
        do j = jb - jh, je + jh
          do i = ib - ih, ie + ih
            if (in_zone(3, i, j, k)) cycle
            if (wp(i,j,k) /= sw(i,j,k)) n = n + 1
          end do
        end do
      end do

      offzone_changes = nest_sumall(real(n))

    end function offzone_changes

    !> .true. when (ivar, i, j, k) is a point the probe found in the zone.
    logical function in_zone(ivar, i, j, k)
      integer, intent(in) :: ivar, i, j, k

      in_zone = .false.
      if (i < ib .or. i > ie) return
      if (j < jb .or. j > je) return
      if (k < kb .or. k > ke) return
      in_zone = (nest_pick(ivar, wu, wv, ww, i, j, k) > 0.)

    end function in_zone

    !> U29: W = 0 leaves up bitwise unchanged, halo planes included.
    logical function u29_noop()
      real :: n

      nest_tau = 4.
      dt       = 0.5
      rk3step  = 2
      timee    = 17.
      call nest_target_at(17., wu, wv, ww)   ! advance the buffer
      call nest_probe_weight(wu, wv, ww)
      nest_tau = 4.
      dt       = 0.5
      rk3step  = 2

      call set_state(1.)
      call nesting_apply
      n = offzone_changes()

      u29_noop = (n == 0.)
      if (myid == 0) write(*,'(a,i0)') '   off-zone points changed = ', nint(n)
      call nest_report('U29 no-op outside the zone', u29_noop)
    end function u29_noop

    !> U30: tau -> 0 (or W dt/tau -> infinity) gives q_new = target exactly.
    logical function u30_dirichlet(tau)
      real, intent(in) :: tau

      integer :: ivar, i, j, k
      real    :: rk3coef, x, y, z, tgt, qnew, dmax
      character(len=40) :: lbl

      write(lbl,'(a,es9.2,a)') 'U30 Dirichlet limit (tau = ', tau, ')'

      nest_tau = tau
      dt       = 0.5
      rk3step  = 2
      rk3coef  = dt/2.

      call set_state(1.)
      call nesting_apply

      dmax = 0.
      do ivar = 1, 3
        do k = kb, ke
          do j = jb, je
            do i = ib, ie
              if (.not. in_zone(ivar, i, j, k)) cycle
              call nest_ref_coord(ivar, i, j, k, x, y, z)
              tgt = nest_analytic(ivar, x, y, z, 17.)
              select case (ivar)
              case (1)
                qnew = mu(i,j,k) + rk3coef*up(i,j,k)
              case (2)
                qnew = mv(i,j,k) + rk3coef*vp(i,j,k)
              case default
                qnew = mw(i,j,k) + rk3coef*wp(i,j,k)
              end select
              dmax = max(dmax, abs(qnew - tgt))
            end do
          end do
        end do
      end do
      dmax = nest_maxall(dmax)

      u30_dirichlet = dmax <= 1.e-13
      if (myid == 0) write(*,'(a,es12.4)') '   max |q_new - target| = ', dmax
      call nest_report(trim(lbl), u30_dirichlet)
    end function u30_dirichlet

    !> U31: for W dt_s / tau <= 1e-4 the update agrees with the explicit
    !! linear forcing qp += (W/tau)(target - q*) to second order in that
    !! small parameter.
    logical function u31_linear()
      integer :: ivar, i, j, k
      real    :: rk3coef, x, y, z, tgt, qstar, w, eps, lin, got, err, bound, worst

      nest_tau = 1.e6
      dt       = 1.
      rk3step  = 2
      rk3coef  = dt/2.

      call set_state(1.)
      call nesting_apply

      worst = 0.
      do ivar = 1, 3
        do k = kb, ke
          do j = jb, je
            do i = ib, ie
              if (.not. in_zone(ivar, i, j, k)) cycle
              w = nest_pick(ivar, wu, wv, ww, i, j, k)
              call nest_ref_coord(ivar, i, j, k, x, y, z)
              tgt = nest_analytic(ivar, x, y, z, 17.)
              select case (ivar)
              case (1)
                qstar = mu(i,j,k) + rk3coef*su(i,j,k)
                lin   = su(i,j,k) + (w/nest_tau)*(tgt - qstar)
                got   = up(i,j,k)
              case (2)
                qstar = mv(i,j,k) + rk3coef*sv(i,j,k)
                lin   = sv(i,j,k) + (w/nest_tau)*(tgt - qstar)
                got   = vp(i,j,k)
              case default
                qstar = mw(i,j,k) + rk3coef*sw(i,j,k)
                lin   = sw(i,j,k) + (w/nest_tau)*(tgt - qstar)
                got   = wp(i,j,k)
              end select
              eps   = w*rk3coef/nest_tau
              err   = abs(got - lin)
              bound = abs(tgt - qstar)*eps*w/nest_tau + 1.e-14
              worst = max(worst, err/bound)
            end do
          end do
        end do
      end do
      worst = nest_maxall(worst)

      u31_linear = worst <= 1.
      if (myid == 0) write(*,'(a,es12.4)')                                  &
        '   max |update - linear form| / O(eps^2) bound = ', worst
      call nest_report('U31 linear limit', u31_linear)
    end function u31_linear

    !> U32: driving a point through all three RK3 substeps with N = 0 gives
    !! exactly target + (q^n - target) exp(-W dt / tau) -- the C10 property.
    logical function u32_fullstep()
      integer :: ivar, i, j, k, s
      real    :: rk3coef, x, y, z, tgt, w, dmax, qend, expect

      nest_tau = 4.
      dt       = 1.2

      call set_state(1.)
      up = 0.; vp = 0.; wp = 0.
      su = 0.; sv = 0.; sw = 0.

      dmax = 0.
      do s = 1, 3
        rk3step = s
        rk3coef = dt/(4. - real(s))
        up = 0.; vp = 0.; wp = 0.        ! tstep_integrate resets the tendencies
        call nesting_apply
        if (s < 3) cycle
        do ivar = 1, 3
          do k = kb, ke
            do j = jb, je
              do i = ib, ie
                if (.not. in_zone(ivar, i, j, k)) cycle
                w = nest_pick(ivar, wu, wv, ww, i, j, k)
                call nest_ref_coord(ivar, i, j, k, x, y, z)
                tgt = nest_analytic(ivar, x, y, z, 17.)
                select case (ivar)
                case (1)
                  qend   = mu(i,j,k) + rk3coef*up(i,j,k)
                  expect = tgt + (mu(i,j,k) - tgt)*exp(-w*dt/nest_tau)
                case (2)
                  qend   = mv(i,j,k) + rk3coef*vp(i,j,k)
                  expect = tgt + (mv(i,j,k) - tgt)*exp(-w*dt/nest_tau)
                case default
                  qend   = mw(i,j,k) + rk3coef*wp(i,j,k)
                  expect = tgt + (mw(i,j,k) - tgt)*exp(-w*dt/nest_tau)
                end select
                dmax = max(dmax, abs(qend - expect))
              end do
            end do
          end do
        end do
      end do
      dmax = nest_maxall(dmax)

      u32_fullstep = dmax <= 1.e-13
      if (myid == 0) write(*,'(a,es12.4)')                                  &
        '   max |q after 3 substeps - target - (q^n - target) exp(-W dt/tau)| = ', dmax
      call nest_report('U32 full-step composition', u32_fullstep)
    end function u32_fullstep

    !> U33: no overshoot, no sign change and no non-finite value for
    !! dt/tau between 1e-3 and 1e6.
    logical function u33_stability()
      integer :: n, ivar, i, j, k
      real    :: ratio, x, y, z, tgt, qnew, q0, nbad, ov, ovmax, sgmax, tolabs

      q0    = 0.3
      nbad  = 0.
      ovmax = 0.
      sgmax = 0.
      dt      = 1.
      rk3step = 3

      do n = 0, 18
        ratio    = 10.**(-3. + 0.5*real(n))
        nest_tau = dt/ratio

        um = q0; vm = q0; wm = q0
        up = 0.; vp = 0.; wp = 0.
        call nesting_apply

        do ivar = 1, 3
          do k = kb, ke
            do j = jb, je
              do i = ib, ie
                if (.not. in_zone(ivar, i, j, k)) cycle
                call nest_ref_coord(ivar, i, j, k, x, y, z)
                tgt = nest_analytic(ivar, x, y, z, 17.)
                select case (ivar)
                case (1)
                  qnew = q0 + dt*up(i,j,k)
                case (2)
                  qnew = q0 + dt*vp(i,j,k)
                case default
                  qnew = q0 + dt*wp(i,j,k)
                end select
                tolabs = 8.*epsilon(1.)*max(abs(q0), abs(tgt), abs(qnew))
                ov = abs(qnew - tgt) - abs(q0 - tgt)
                ovmax = max(ovmax, ov)
                if (.not. (abs(qnew) < huge(1.))) nbad = nbad + 1.
                if ((qnew - tgt)*(q0 - tgt) < 0.) then
                  sgmax = max(sgmax, abs(qnew - tgt))
                  if (abs(qnew - tgt) > tolabs) nbad = nbad + 1.
                end if
                if (ov > tolabs) nbad = nbad + 1.
              end do
            end do
          end do
        end do
      end do

      nbad  = nest_sumall(nbad)
      ovmax = nest_maxall(ovmax)
      sgmax = nest_maxall(sgmax)

      u33_stability = (nbad == 0.)
      if (myid == 0) write(*,'(a,es12.4,a,es12.4)')                         &
        '   worst overshoot = ', ovmax, ', worst wrong-sign residual = ', sgmax
      if (myid == 0) write(*,'(a,i0)')                                      &
        '   overshoot / sign-change / non-finite violations = ', nint(nbad)
      call nest_report('U33 stability over dt/tau in [1e-3, 1e6]', u33_stability)
    end function u33_stability

    !> U34: solid points are not touched by nesting_apply.
    logical function u34_solid()
      integer :: i, j, k, ig, jg, n
      real    :: gn

      call nest_set_solids(.true.)
      nest_lparentgeom = .true.
      call nest_reinit('nesting_analytic.'//cexpnr//'.nc', 17., 4., 2, 1, 1.e30, .false.)

      nest_tau = 4.
      dt       = 0.5
      rk3step  = 2
      call set_state(1.)
      call nesting_apply

      n = 0
      do k = kb, ke
        do j = jb, je
          do i = ib, ie
            ig = i + zstart(1) - 1
            jg = j + zstart(2) - 1
            if (.not. nest_ref_solid(ig, jg, k)) cycle
            if (up(i,j,k) /= su(i,j,k)) n = n + 1
            if (vp(i,j,k) /= sv(i,j,k)) n = n + 1
            if (wp(i,j,k) /= sw(i,j,k)) n = n + 1
          end do
        end do
      end do

      gn = nest_sumall(real(n))
      u34_solid = (gn == 0.)
      if (myid == 0) write(*,'(a,i0)') '   solid points changed = ', nint(gn)
      call nest_report('U34 solid points untouched', u34_solid)

      call nest_set_solids(.false.)
    end function u34_solid

  end function tests_nesting_update

  !> Cold-start initialisation of the interior from the parent, design section
  !! 10.6 item 4. Covers U40-U43.
  !!
  !! The fixture nesting_initial.<expnr>.nc carries the analytic field of
  !! udprep.nesting.analytic_field as its full-domain block, written WITHOUT the
  !! divergence correction, so every value the reader is supposed to place can be
  !! predicted here exactly and a transposed, off-by-one or wrongly staggered
  !! read cannot pass. The abort cases (no block, wrong shape, wrong stagger) are
  !! separate invocations selected through the nestfile namelist entry, since
  !! they stop the process.
  logical function tests_nesting_init()
    use mpi
    use modglobal,  only : ib, ie, ih, jb, je, jh, kb, ke, kh, cexpnr,      &
                           timee, lwarmstart, ierank, jerank
    use modfields,  only : initfields, u0, um, v0, vm, w0, wm
    use modibm,     only : createmasks
    use modnesting, only : lnesting, nestfile, nest_tau, nest_fluxtol,       &
                           nest_lfluxassert, nest_linitfromparent,           &
                           nesting_init, nesting_finalize

    implicit none

    real, parameter :: SENTINEL_U = 3.75, SENTINEL_V = -1.25, SENTINEL_W = 0.5

    logical :: all_passed, lwarm_save
    real, allocatable :: su(:,:,:), sv(:,:,:), sw(:,:,:)

    call nest_banner('tests_nesting_init', 'COLD-START INIT FROM THE PARENT (U40-U43)')

    call initfields
    call createmasks
    call nest_set_solids(.false.)

    ! The failure modes -- no block, wrong shape, wrong stagger -- stop the
    ! process, so the driver runs each of them as its own invocation and selects
    ! it by pointing nestfile somewhere other than the good fixture.
    if (index(nestfile, 'nesting_initial') == 0) then
      tests_nesting_init = u43_init_aborts()
      return
    end if

    allocate(su(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
    allocate(sv(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
    allocate(sw(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))

    lwarm_save = lwarmstart
    all_passed = .true.

    if (.not. u40_fills_from_the_block()) all_passed = .false.
    if (.not. u41_off_leaves_the_fields()) all_passed = .false.
    if (.not. u42_warmstart_is_not_touched()) all_passed = .false.

    lwarmstart = lwarm_save
    deallocate(su, sv, sw)
    call nesting_finalize

    call nest_verdict('tests_nesting_init', all_passed)
    tests_nesting_init = nest_all_ranks(all_passed)

  contains

    !> Fill the velocity fields with a value nothing in the file can produce and
    !! remember them, so any point the reader does not write is recognisable.
    subroutine plant_sentinel()
      u0 = SENTINEL_U; um = SENTINEL_U
      v0 = SENTINEL_V; vm = SENTINEL_V
      w0 = SENTINEL_W; wm = SENTINEL_W
      su = u0; sv = v0; sw = w0
    end subroutine plant_sentinel

    !> Re-initialise nesting with the initial-condition switch in a given state.
    subroutine reinit(fname, lfrom, lwarm)
      character(len=*), intent(in) :: fname
      logical,          intent(in) :: lfrom, lwarm

      call nesting_finalize
      lnesting             = .true.
      nestfile             = fname
      nest_tau             = 4.
      nest_fluxtol         = 1.e30
      nest_lfluxassert     = .false.
      nest_linitfromparent = lfrom
      lwarmstart           = lwarm
      timee                = 0.
      call nesting_init
    end subroutine reinit

    !> Largest deviation of one component from the analytic field over the
    !! index range the reader is contractually required to fill, MPI-reduced.
    !! iu/ju/ku extend the local range by one where this rank owns the far face.
    real function block_error(ivar, q)
      integer, intent(in) :: ivar
      real,    intent(in) :: q(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)

      integer :: i, j, k, iend, jend, kend
      real    :: x, y, z, d, dmax

      iend = ie
      jend = je
      kend = ke
      if (ivar == 1 .and. ierank) iend = ie + 1
      if (ivar == 2 .and. jerank) jend = je + 1
      if (ivar == 3) kend = ke + 1

      dmax = 0.
      do k = kb, kend
        do j = jb, jend
          do i = ib, iend
            call nest_ref_coord(ivar, i, j, k, x, y, z)
            d = abs(q(i,j,k) - nest_analytic(ivar, x, y, z, 0.))
            dmax = max(dmax, d)
          end do
        end do
      end do

      block_error = nest_maxall(dmax)

    end function block_error

    !> U43: nesting_init must stop on a file that cannot serve the switch --
    !! it has no block, or one at the wrong shape or the wrong stagger. Reaching
    !! the line after nesting_init is the failure.
    logical function u43_init_aborts()
      lnesting             = .true.
      nest_tau             = 4.
      nest_fluxtol         = 1.e30
      nest_lfluxassert     = .false.
      nest_linitfromparent = .true.
      lwarmstart           = .false.
      timee                = 0.
      if (myid == 0) write(*,'(a,a)')                                       &
        ' U43: nesting_init must abort on ', trim(nestfile)

      call nesting_init

      call nest_report('U43 cold-start init (init did NOT abort)', .false.)
      call nest_verdict('tests_nesting_init', .false.)
      u43_init_aborts = .false.
    end function u43_init_aborts

    !> U40: with nest_linitfromparent the three components come back as the
    !! stored block, and the m-level copies are bitwise identical to the 0-level.
    logical function u40_fills_from_the_block()
      real :: du, dv, dw
      logical :: lm

      call plant_sentinel()
      call reinit('nesting_initial.'//cexpnr//'.nc', .true., .false.)

      du = block_error(1, u0)
      dv = block_error(2, v0)
      dw = block_error(3, w0)
      lm = nest_sumall(real(count(um /= u0) + count(vm /= v0) + count(wm /= w0))) == 0.

      u40_fills_from_the_block = (max(du, dv, dw) <= 1.e-14) .and. lm

      if (myid == 0) then
        write(*,'(a,3es12.4)') '   max |q - analytic| u/v/w = ', du, dv, dw
        if (.not. lm) write(*,'(a)') '   the m-level fields differ from the 0-level ones'
      end if
      call nest_report('U40 cold start filled from the parent block', &
                       u40_fills_from_the_block)
    end function u40_fills_from_the_block

    !> U41: with the switch off nothing is touched, bitwise. The same file is
    !! used, so this isolates the switch and not the file.
    logical function u41_off_leaves_the_fields()
      real :: n

      call plant_sentinel()
      call reinit('nesting_initial.'//cexpnr//'.nc', .false., .false.)

      n = nest_sumall(real(count(u0 /= su) + count(v0 /= sv) + count(w0 /= sw) + &
                           count(um /= su) + count(vm /= sv) + count(wm /= sw)))
      u41_off_leaves_the_fields = (n == 0.)
      if (myid == 0) write(*,'(a,i0)') '   points changed with the switch off = ', nint(n)
      call nest_report('U41 switch off leaves the fields bitwise unchanged', &
                       u41_off_leaves_the_fields)
    end function u41_off_leaves_the_fields

    !> U42: a warm start already holds a state consistent with the parent, and
    !! overwriting it would break restart parity (I6). The switch must be
    !! ignored, not honoured.
    logical function u42_warmstart_is_not_touched()
      real :: n

      call plant_sentinel()
      call reinit('nesting_initial.'//cexpnr//'.nc', .true., .true.)

      n = nest_sumall(real(count(u0 /= su) + count(v0 /= sv) + count(w0 /= sw) + &
                           count(um /= su) + count(vm /= sv) + count(wm /= sw)))
      u42_warmstart_is_not_touched = (n == 0.)
      if (myid == 0) write(*,'(a,i0)') '   points changed on a warm start = ', nint(n)
      call nest_report('U42 warm start is not overwritten', &
                       u42_warmstart_is_not_touched)
    end function u42_warmstart_is_not_touched

  end function tests_nesting_init

end module tests
