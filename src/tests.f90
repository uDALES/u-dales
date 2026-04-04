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
  use decomp_2d
  use modmpi, only : myid
  use architecture, only : comm3d, mpierr, my_real, nprocs, bcast, allreduce, global_sum, global_max
  
  implicit none
  save
  public :: tests_read_sparse_ijk, tests_2decomp_init_exit, tests_mpi_fieldops, tests_mpi_architecture

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
  !> then compares with the new generic read_sparse_ijk function
  !> Returns .true. if all tests pass, .false. otherwise
  logical function tests_read_sparse_ijk()
    use modglobal,    only : libm, cexpnr, runmode
    use readinput, only : read_sparse_ijk
    use ibm,       only : initibm
    use ibm,       only : solid_info_u, solid_info_v, solid_info_w, solid_info_c
    use ibm,       only : bound_info_u, bound_info_v, bound_info_w, bound_info_c
    use ibm,       only : nsolpts_u, nsolpts_v, nsolpts_w, nsolpts_c
    use ibm,       only : nbndpts_u, nbndpts_v, nbndpts_w, nbndpts_c
    use initfac,      only : readfacetfiles
    use decomp_2d,    only : zstart, zend, xstart, xend, ystart, yend
    
    implicit none
    
    integer :: npts_loc_new
    integer, allocatable :: ids_loc_new(:)
    integer, allocatable :: pts_loc_new(:,:)
    integer :: m
    logical :: all_passed
    character(200) :: filename
    
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
    call read_test_output('solid_u', npts_loc_new, ids_loc_new, pts_loc_new)
    if (.not. compare_solid(solid_info_u, npts_loc_new, ids_loc_new, pts_loc_new, 'solid_u')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test solid_v
    call read_test_output('solid_v', npts_loc_new, ids_loc_new, pts_loc_new)
    if (.not. compare_solid(solid_info_v, npts_loc_new, ids_loc_new, pts_loc_new, 'solid_v')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test solid_w
    call read_test_output('solid_w', npts_loc_new, ids_loc_new, pts_loc_new)
    if (.not. compare_solid(solid_info_w, npts_loc_new, ids_loc_new, pts_loc_new, 'solid_w')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test solid_c
    call read_test_output('solid_c', npts_loc_new, ids_loc_new, pts_loc_new)
    if (.not. compare_solid(solid_info_c, npts_loc_new, ids_loc_new, pts_loc_new, 'solid_c')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test fluid_boundary_u
    call read_test_output('fluid_boundary_u', npts_loc_new, ids_loc_new, pts_loc_new)
    if (.not. compare_boundary(bound_info_u, npts_loc_new, ids_loc_new, pts_loc_new, 'fluid_boundary_u')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test fluid_boundary_v
    call read_test_output('fluid_boundary_v', npts_loc_new, ids_loc_new, pts_loc_new)
    if (.not. compare_boundary(bound_info_v, npts_loc_new, ids_loc_new, pts_loc_new, 'fluid_boundary_v')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test fluid_boundary_w
    call read_test_output('fluid_boundary_w', npts_loc_new, ids_loc_new, pts_loc_new)
    if (.not. compare_boundary(bound_info_w, npts_loc_new, ids_loc_new, pts_loc_new, 'fluid_boundary_w')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test fluid_boundary_c
    call read_test_output('fluid_boundary_c', npts_loc_new, ids_loc_new, pts_loc_new)
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

  !> Read test output from file with rank coordinates in filename
  subroutine read_test_output(label, npts_loc, ids_loc, pts_loc)
    use modmpi, only : myidx, myidy
    
    implicit none
    character(len=*), intent(in) :: label
    integer, intent(out) :: npts_loc
    integer, allocatable, intent(out) :: ids_loc(:), pts_loc(:,:)
    
    character(200) :: filename
    character(500) :: headerline
    integer :: funit, m, ierr, colon_pos
    
    ! Create filename with cartesian rank coordinates: label_X_Y.txt
    write(filename, '(A,A,I0,A,I0,A)') trim(label), '_', myidx, '_', myidy, '.txt'
    
    ! Open file and read data
    open(newunit=funit, file=trim(filename), status='old', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*, '(A,A)') 'ERROR: Cannot open test output file: ', trim(filename)
      stop 1
    end if
    
    ! Read header lines
    read(funit, '(A500)', iostat=ierr) headerline  ! Line 1: Test output for...
    read(funit, '(A500)', iostat=ierr) headerline  ! Line 2: Rank cartesian coordinates
    read(funit, '(A500)', iostat=ierr) headerline  ! Line 3: Number of local points
    if (ierr /= 0) then
      write(*, '(A)') 'ERROR: Cannot read header line 3'
      stop 1
    end if
    
    ! Extract npts_loc from line 3
    ! Format is: "# Number of local points: N"
    colon_pos = index(headerline, ':')
    if (colon_pos > 0) then
      read(headerline(colon_pos+1:), *, iostat=ierr) npts_loc
      if (ierr /= 0) then
        write(*, '(A,A)') 'ERROR: Cannot parse number of points from: ', trim(headerline)
        stop 1
      end if
    else
      write(*, '(A,A,A)') 'ERROR: Cannot find colon in header: "', trim(headerline), '"'
      stop 1
    end if
    
    read(funit, '(A500)', iostat=ierr) headerline  ! Line 4: Format description
    
    ! Allocate arrays
    allocate(ids_loc(npts_loc))
    allocate(pts_loc(npts_loc, 3))
    
    ! Read data
    if (npts_loc > 0) then
      do m = 1, npts_loc
        read(funit, *, iostat=ierr) ids_loc(m), pts_loc(m,1), pts_loc(m,2), pts_loc(m,3)
        if (ierr /= 0) then
          write(*, '(A,I0,A,A)') 'ERROR: Cannot read data line ', m, ' from ', trim(filename)
          stop 1
        end if
      end do
    end if
    
    close(funit)
    
  end subroutine read_test_output
  
  function compare_solid(solid_info, npts_loc_new, ids_loc_new, pts_loc_new, label) result(passed)
    use modmpi, only : myid
    use ibm, only : solid_info_type
    use decomp_2d, only : zstart, zend
    
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
      stop 1
    end if
    
    do m = 1, solid_info%nsolptsrank
      if (ids_loc_new(m) /= solid_info%solptsrank(m)) then
        write(*, '(A,I0,A,A,A)') 'FAIL on rank ', myid, ': ', trim(label), ' index mismatch'
        passed = .false.
        stop 1
      end if
      if (pts_loc_new(m,1) /= solid_info%solpts_loc(m,1) .or. &
          pts_loc_new(m,2) /= solid_info%solpts_loc(m,2) .or. &
          pts_loc_new(m,3) /= solid_info%solpts_loc(m,3)) then
        write(*, '(A,I0,A,A,A)') 'FAIL on rank ', myid, ': ', trim(label), ' coordinate mismatch'
        passed = .false.
        stop 1
      end if
    end do
    
  end function compare_solid
  
  function compare_boundary(bound_info, npts_loc_new, ids_loc_new, pts_loc_new, label) result(passed)
    use modmpi, only : myid
    use ibm, only : bound_info_type
    use decomp_2d, only : zstart, zend
    
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
      stop 1
    end if
    
    do m = 1, bound_info%nbndptsrank
      if (ids_loc_new(m) /= bound_info%bndptsrank(m)) then
        write(*, '(A,I0,A,A,A)') 'FAIL on rank ', myid, ': ', trim(label), ' index mismatch'
        passed = .false.
        stop 1
      end if
      if (pts_loc_new(m,1) /= bound_info%bndpts_loc(m,1) .or. &
          pts_loc_new(m,2) /= bound_info%bndpts_loc(m,2) .or. &
          pts_loc_new(m,3) /= bound_info%bndpts_loc(m,3)) then
        write(*, '(A,I0,A,A,A)') 'FAIL on rank ', myid, ': ', trim(label), ' coordinate mismatch'
        passed = .false.
        stop 1
      end if
    end do
    
  end function compare_boundary

  logical function tests_mpi_architecture()
    use mpi, only : MPI_SUM, MPI_MAX
    use modglobal, only : runmode
    implicit none

    logical :: all_passed

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_mpi_architecture: DIRECT MPI ABSTRACTION TEST'
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'Testing generic bcast/allreduce/global_sum/global_max interfaces'
    end if

    all_passed = .true.

    if (.not. check_bcast_scalars()) all_passed = .false.
    if (.not. check_bcast_arrays()) all_passed = .false.
    if (.not. check_allreduce_scalars()) all_passed = .false.
    if (.not. check_allreduce_arrays()) all_passed = .false.
    if (.not. check_global_sum_scalars()) all_passed = .false.
    if (.not. check_global_sum_arrays()) all_passed = .false.
    if (.not. check_global_max_scalars()) all_passed = .false.

    if (all_passed .and. myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'ALL TESTS PASSED: tests_mpi_architecture'
      write(*, '(A)') '================================================'
    else if ((.not. all_passed) .and. myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'TESTS FAILED: tests_mpi_architecture'
      write(*, '(A)') '================================================'
    end if

    tests_mpi_architecture = all_passed

  contains

    logical function check_bcast_scalars()
      implicit none
      real :: real_val
      integer :: int_val
      logical :: logical_val
      character(16) :: char_val

      real_val = -1.
      int_val = -1
      logical_val = .false.
      char_val = 'unset'
      if (myid == 0) then
        real_val = 3.25
        int_val = 42
        logical_val = .true.
        char_val = 'arch-bcast-ok'
      end if

      call bcast(real_val, 0)
      call bcast(int_val, 0)
      call bcast(logical_val, 0)
      call bcast(char_val, 0)

      check_bcast_scalars = .true.
      if (.not. nearly_equal_scalar_arch(real_val, 3.25)) then
        call report_scalar_failure_arch('bcast real scalar', real_val, 3.25)
        check_bcast_scalars = .false.
      end if
      if (int_val /= 42) then
        call report_int_scalar_failure('bcast int scalar', int_val, 42)
        check_bcast_scalars = .false.
      end if
      if (.not. logical_val) then
        call report_logical_failure('bcast logical scalar', logical_val, .true.)
        check_bcast_scalars = .false.
      end if
      if (trim(char_val) /= 'arch-bcast-ok') then
        call report_char_failure('bcast char scalar', char_val, 'arch-bcast-ok')
        check_bcast_scalars = .false.
      end if
    end function check_bcast_scalars

    logical function check_bcast_arrays()
      implicit none
      real :: real_vec(3)
      integer :: int_vec(4)
      logical :: logical_vec(3)

      real_vec = -1.
      int_vec = -1
      logical_vec = .false.
      if (myid == 0) then
        real_vec = (/ 1.5, -2.5, 4.0 /)
        int_vec = (/ 2, 4, 6, 8 /)
        logical_vec = (/ .true., .false., .true. /)
      end if

      call bcast(real_vec, 0)
      call bcast(int_vec, 0)
      call bcast(logical_vec, 0)

      check_bcast_arrays = .true.
      if (.not. compare_real_1d_arch('bcast real 1d', real_vec, (/ 1.5, -2.5, 4.0 /))) check_bcast_arrays = .false.
      if (.not. compare_int_1d('bcast int 1d', int_vec, (/ 2, 4, 6, 8 /))) check_bcast_arrays = .false.
      if (.not. compare_logical_1d('bcast logical 1d', logical_vec, (/ .true., .false., .true. /))) check_bcast_arrays = .false.
    end function check_bcast_arrays

    logical function check_allreduce_scalars()
      implicit none
      real :: real_local, real_sum, real_max
      integer :: int_local, int_sum, int_max

      real_local = 0.5 * real(myid + 1)
      int_local = myid + 1
      real_sum = allreduce(real_local, MPI_SUM)
      real_max = allreduce(real_local, MPI_MAX)
      int_sum = allreduce(int_local, MPI_SUM)
      int_max = allreduce(int_local, MPI_MAX)

      check_allreduce_scalars = .true.
      if (.not. nearly_equal_scalar_arch(real_sum, real(nprocs * (nprocs + 1)) / 4.0)) then
        call report_scalar_failure_arch('allreduce real scalar sum', real_sum, real(nprocs * (nprocs + 1)) / 4.0)
        check_allreduce_scalars = .false.
      end if
      if (.not. nearly_equal_scalar_arch(real_max, 0.5 * real(nprocs))) then
        call report_scalar_failure_arch('allreduce real scalar max', real_max, 0.5 * real(nprocs))
        check_allreduce_scalars = .false.
      end if
      if (int_sum /= nprocs * (nprocs + 1) / 2) then
        call report_int_scalar_failure('allreduce int scalar sum', int_sum, nprocs * (nprocs + 1) / 2)
        check_allreduce_scalars = .false.
      end if
      if (int_max /= nprocs) then
        call report_int_scalar_failure('allreduce int scalar max', int_max, nprocs)
        check_allreduce_scalars = .false.
      end if
    end function check_allreduce_scalars

    logical function check_allreduce_arrays()
      implicit none
      real :: real_local(3), real_sum(3), real_max(3), real_exp_sum(3), real_exp_max(3)
      integer :: int_local(3), int_sum(3), int_max(3), int_exp_sum(3), int_exp_max(3)
      integer :: p

      real_local = (/ real(myid + 1), real(2 * (myid + 1)), -real(myid + 1) /)
      int_local = (/ myid + 1, 2 * (myid + 1), -(myid + 1) /)

      real_sum = allreduce(real_local, MPI_SUM)
      real_max = allreduce(real_local, MPI_MAX)
      int_sum = allreduce(int_local, MPI_SUM)
      int_max = allreduce(int_local, MPI_MAX)

      real_exp_sum = 0.
      int_exp_sum = 0
      do p = 1, nprocs
        real_exp_sum = real_exp_sum + (/ real(p), real(2 * p), -real(p) /)
        int_exp_sum = int_exp_sum + (/ p, 2 * p, -p /)
      end do
      real_exp_max = (/ real(nprocs), real(2 * nprocs), -1.0 /)
      int_exp_max = (/ nprocs, 2 * nprocs, -1 /)

      check_allreduce_arrays = .true.
      if (.not. compare_real_1d_arch('allreduce real 1d sum', real_sum, real_exp_sum)) check_allreduce_arrays = .false.
      if (.not. compare_real_1d_arch('allreduce real 1d max', real_max, real_exp_max)) check_allreduce_arrays = .false.
      if (.not. compare_int_1d('allreduce int 1d sum', int_sum, int_exp_sum)) check_allreduce_arrays = .false.
      if (.not. compare_int_1d('allreduce int 1d max', int_max, int_exp_max)) check_allreduce_arrays = .false.
    end function check_allreduce_arrays

    logical function check_global_sum_scalars()
      implicit none
      real :: real_got
      integer :: int_got

      real_got = global_sum(0.25 * real(myid + 1))
      int_got = global_sum(myid + 1)

      check_global_sum_scalars = .true.
      if (.not. nearly_equal_scalar_arch(real_got, 0.125 * real(nprocs * (nprocs + 1)))) then
        call report_scalar_failure_arch('global_sum real scalar', real_got, 0.125 * real(nprocs * (nprocs + 1)))
        check_global_sum_scalars = .false.
      end if
      if (int_got /= nprocs * (nprocs + 1) / 2) then
        call report_int_scalar_failure('global_sum int scalar', int_got, nprocs * (nprocs + 1) / 2)
        check_global_sum_scalars = .false.
      end if
    end function check_global_sum_scalars

    logical function check_global_sum_arrays()
      implicit none
      real :: real_local(3), real_got(3), real_exp(3)
      integer :: int_local(3), int_got(3), int_exp(3)
      integer :: p

      real_local = (/ real(myid + 1), real(myid + 2), real(2 * myid + 1) /)
      int_local = (/ myid + 1, myid + 2, 2 * myid + 1 /)
      real_got = global_sum(real_local)
      int_got = global_sum(int_local)

      real_exp = 0.
      int_exp = 0
      do p = 0, nprocs - 1
        real_exp = real_exp + (/ real(p + 1), real(p + 2), real(2 * p + 1) /)
        int_exp = int_exp + (/ p + 1, p + 2, 2 * p + 1 /)
      end do

      check_global_sum_arrays = .true.
      if (.not. compare_real_1d_arch('global_sum real 1d', real_got, real_exp)) check_global_sum_arrays = .false.
      if (.not. compare_int_1d('global_sum int 1d', int_got, int_exp)) check_global_sum_arrays = .false.
    end function check_global_sum_arrays

    logical function check_global_max_scalars()
      implicit none
      real :: real_got
      integer :: int_got

      real_got = global_max(1.5 * real(myid + 1))
      int_got = global_max(myid + 1)

      check_global_max_scalars = .true.
      if (.not. nearly_equal_scalar_arch(real_got, 1.5 * real(nprocs))) then
        call report_scalar_failure_arch('global_max real scalar', real_got, 1.5 * real(nprocs))
        check_global_max_scalars = .false.
      end if
      if (int_got /= nprocs) then
        call report_int_scalar_failure('global_max int scalar', int_got, nprocs)
        check_global_max_scalars = .false.
      end if
    end function check_global_max_scalars

    logical function compare_real_1d_arch(label, got, exp)
      implicit none
      character(len=*), intent(in) :: label
      real, intent(in)             :: got(:), exp(:)
      real                         :: max_abs
      integer                      :: imax(1)

      max_abs = maxval(abs(got - exp))
      compare_real_1d_arch = max_abs <= 1.e-9
      if ((.not. compare_real_1d_arch) .and. myid == 0) then
        imax = maxloc(abs(got - exp))
        write(*,'(A,1X,A,1X,ES12.4,1X,A,I0,1X,A,ES12.4,1X,A,ES12.4)') &
             'FAIL', trim(label), max_abs, 'idx', imax(1), 'got', got(imax(1)), 'exp', exp(imax(1))
      end if
    end function compare_real_1d_arch

    logical function compare_int_1d(label, got, exp)
      implicit none
      character(len=*), intent(in) :: label
      integer, intent(in)          :: got(:), exp(:)
      integer                      :: max_abs
      integer                      :: imax(1)

      max_abs = maxval(abs(got - exp))
      compare_int_1d = max_abs == 0
      if ((.not. compare_int_1d) .and. myid == 0) then
        imax = maxloc(abs(got - exp))
        write(*,'(A,1X,A,1X,I0,1X,A,I0,1X,A,I0,1X,A,I0)') &
             'FAIL', trim(label), max_abs, 'idx', imax(1), 'got', got(imax(1)), 'exp', exp(imax(1))
      end if
    end function compare_int_1d

    logical function compare_logical_1d(label, got, exp)
      implicit none
      character(len=*), intent(in) :: label
      logical, intent(in)          :: got(:), exp(:)
      integer :: i

      compare_logical_1d = .true.
      do i = 1, size(got)
        if (got(i) .neqv. exp(i)) then
          compare_logical_1d = .false.
          if (myid == 0) then
            write(*,'(A,1X,A,1X,A,I0)') 'FAIL', trim(label), 'idx', i
          end if
          exit
        end if
      end do
    end function compare_logical_1d

    logical function nearly_equal_scalar_arch(got, exp)
      implicit none
      real, intent(in) :: got, exp
      nearly_equal_scalar_arch = abs(got - exp) <= 1.e-9
    end function nearly_equal_scalar_arch

    subroutine report_scalar_failure_arch(label, got, exp)
      implicit none
      character(len=*), intent(in) :: label
      real, intent(in)             :: got, exp

      if (myid == 0) then
        write(*,'(A,1X,A,1X,ES12.4,1X,ES12.4)') 'FAIL', trim(label), got, exp
      end if
    end subroutine report_scalar_failure_arch

    subroutine report_int_scalar_failure(label, got, exp)
      implicit none
      character(len=*), intent(in) :: label
      integer, intent(in)          :: got, exp

      if (myid == 0) then
        write(*,'(A,1X,A,1X,I0,1X,I0)') 'FAIL', trim(label), got, exp
      end if
    end subroutine report_int_scalar_failure

    subroutine report_logical_failure(label, got, exp)
      implicit none
      character(len=*), intent(in) :: label
      logical, intent(in)          :: got, exp

      if (myid == 0) then
        write(*,'(A,1X,A,1X,L1,1X,L1)') 'FAIL', trim(label), got, exp
      end if
    end subroutine report_logical_failure

    subroutine report_char_failure(label, got, exp)
      implicit none
      character(len=*), intent(in) :: label, got, exp

      if (myid == 0) then
        write(*,'(A,1X,A,1X,A,1X,A)') 'FAIL', trim(label), trim(got), trim(exp)
      end if
    end subroutine report_char_failure

  end function tests_mpi_architecture

  logical function tests_mpi_fieldops()
    use mpi
    use modglobal,   only : ib, ie, jb, je, kb, ke, kh, runmode
    use readinput,   only : read_sparse_ijk
    use ibm,         only : initibm_support, createmasks, &
                            solid_info_type, solid_info_u, solid_info_v, solid_info_w, solid_info_c, &
                            nsolpts_u, nsolpts_v, nsolpts_w, nsolpts_c
    use ibmmasks,    only : IIc, IIu, IIv, IIw, IIuv, IIuw, IIvw, &
                            IIct, IIut, IIvt, IIwt, IIuwt, &
                            IIcs, IIus, IIvs, IIws, IIuvs, IIuws, IIvws
    use definitions, only : LOC_C, LOC_U, LOC_V, LOC_W, LOC_UV, LOC_WU, LOC_VW
    use operators,    only : reduce_xy_sum, reduce_yz_sum, avg_xy_fluid, avg_y_fluid, sum_x_fluid, sum_y_fluid

    implicit none

    logical :: all_passed

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_mpi_fieldops: DIRECT MPI FIELD OPERATOR TEST'
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'Using Xie/Castro geometry from case 100 with explicit masked references'
      write(*, '(A,6(1X,I0))') 'domain bounds ib ie jb je kb ke =', ib, ie, jb, je, kb, ke
      write(*, '(A,1X,I0)') 'kh =', kh
    end if

    call initibm_support
    call load_solid_info_from_sparse('solid_u.txt', nsolpts_u, solid_info_u)
    call load_solid_info_from_sparse('solid_v.txt', nsolpts_v, solid_info_v)
    call load_solid_info_from_sparse('solid_w.txt', nsolpts_w, solid_info_w)
    call load_solid_info_from_sparse('solid_c.txt', nsolpts_c, solid_info_c)
    call createmasks

    all_passed = .true.

    if (.not. check_global_reductions()) all_passed = .false.
    if (.not. check_reduce_xy_sum()) all_passed = .false.
    if (.not. check_reduce_yz_sum()) all_passed = .false.

    if (.not. check_fluid_loc('LOC_C', LOC_C, IIc, IIcs)) all_passed = .false.
    if (.not. check_fluid_loc('LOC_U', LOC_U, IIu, IIus)) all_passed = .false.
    if (.not. check_fluid_loc('LOC_V', LOC_V, IIv, IIvs)) all_passed = .false.
    if (.not. check_fluid_loc('LOC_W', LOC_W, IIw, IIws)) all_passed = .false.
    if (.not. check_fluid_loc('LOC_UV', LOC_UV, IIuv, IIuvs)) all_passed = .false.
    if (.not. check_fluid_loc('LOC_WU', LOC_WU, IIuw, IIuws)) all_passed = .false.
    if (.not. check_fluid_loc('LOC_VW', LOC_VW, IIvw, IIvws)) all_passed = .false.

    if (.not. check_fluid_y_loc('LOC_C', LOC_C, IIc, IIct)) all_passed = .false.
    if (.not. check_fluid_y_loc('LOC_U', LOC_U, IIu, IIut)) all_passed = .false.
    if (.not. check_fluid_y_loc('LOC_V', LOC_V, IIv, IIvt)) all_passed = .false.
    if (.not. check_fluid_y_loc('LOC_W', LOC_W, IIw, IIwt)) all_passed = .false.
    if (.not. check_fluid_y_loc('LOC_WU', LOC_WU, IIuw, IIuwt)) all_passed = .false.

    if (all_passed .and. myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'ALL TESTS PASSED: tests_mpi_fieldops'
      write(*, '(A)') '================================================'
    else if ((.not. all_passed) .and. myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'TESTS FAILED: tests_mpi_fieldops'
      write(*, '(A)') '================================================'
    end if

    tests_mpi_fieldops = all_passed

  contains

    logical function check_global_reductions()
      implicit none
      real :: local_val, got_sum, got_max, exp_sum, exp_max

      local_val = real(myid + 1)
      got_sum = global_sum(local_val)
      got_max = global_max(local_val)
      exp_sum = real(nprocs * (nprocs + 1)) / 2.
      exp_max = real(nprocs)

      check_global_reductions = .true.
      if (.not. nearly_equal_scalar(got_sum, exp_sum)) then
        call report_scalar_failure('global_sum', got_sum, exp_sum)
        check_global_reductions = .false.
      end if
      if (.not. nearly_equal_scalar(got_max, exp_max)) then
        call report_scalar_failure('global_max', got_max, exp_max)
        check_global_reductions = .false.
      end if
    end function check_global_reductions

    logical function check_reduce_xy_sum()
      use modglobal, only : ib, ie, jb, je, kb, ke
      implicit none
      real, allocatable :: var(:,:,:), got(:), exp_local(:), exp(:)
      integer :: i, j, k

      allocate(var(ib:ie, jb:je, kb:ke))
      allocate(got(kb:ke), exp_local(kb:ke), exp(kb:ke))

      do k = kb, ke
        do j = jb, je
          do i = ib, ie
            var(i,j,k) = real(i) + 10.*real(j) + 100.*real(k)
          end do
        end do
      end do

      got = 0.
      call reduce_xy_sum(got, var)

      do k = kb, ke
        exp_local(k) = sum(var(:,:,k))
      end do
      call MPI_ALLREDUCE(exp_local, exp, size(exp), MY_REAL, MPI_SUM, comm3d, mpierr)

      check_reduce_xy_sum = compare_real_1d('reduce_xy_sum', got, exp)

      deallocate(var, got, exp_local, exp)
    end function check_reduce_xy_sum

    logical function check_reduce_yz_sum()
      use modglobal, only : ib, ie, jb, je, kb, ke
      implicit none
      real, allocatable :: var(:,:,:), got(:), exp_local(:), exp(:)
      integer :: i, j, k

      allocate(var(ib:ie, jb:je, kb:ke))
      allocate(got(ib:ie), exp_local(ib:ie), exp(ib:ie))

      do k = kb, ke
        do j = jb, je
          do i = ib, ie
            var(i,j,k) = 2.*real(i) - 3.*real(j) + 5.*real(k)
          end do
        end do
      end do

      got = 0.
      call reduce_yz_sum(got, var)

      do i = ib, ie
        exp_local(i) = sum(var(i,:,:))
      end do
      call MPI_ALLREDUCE(exp_local, exp, size(exp), MY_REAL, MPI_SUM, comm3d, mpierr)

      check_reduce_yz_sum = compare_real_1d('reduce_yz_sum', got, exp)

      deallocate(var, got, exp_local, exp)
    end function check_reduce_yz_sum

    logical function check_fluid_loc(label, var_loc, mask_3d, mask_1d)
      use modglobal, only : kb, ke, kh
      implicit none
      character(len=*), intent(in) :: label
      integer, intent(in)          :: var_loc
      integer, intent(in)          :: mask_3d(:,:,:), mask_1d(:)
      real, allocatable            :: var_clean(:,:,:), var_halo(:,:,:)
      real, allocatable            :: got_clean(:), got_halo(:), exp(:), sum_local(:), sum_global(:)
      integer                      :: i, j, k
      integer                      :: i1, i2, j1, j2, k1, k2
      integer                      :: ii1, ii2, ij1, ij2, ik1, ik2
      integer                      :: is1, is2

      call select_fluid_xy_reference_windows(mask_real_compatible(mask_3d), mask_3d, mask_1d, zsize(3) + kh, kh, &
           i1, i2, j1, j2, k1, k2, ii1, ii2, ij1, ij2, ik1, ik2, is1, is2)

      allocate(var_clean(lbound(mask_3d,1):ubound(mask_3d,1), &
                         lbound(mask_3d,2):ubound(mask_3d,2), &
                         lbound(mask_3d,3):ubound(mask_3d,3)))
      allocate(var_halo(lbound(mask_3d,1):ubound(mask_3d,1), &
                        lbound(mask_3d,2):ubound(mask_3d,2), &
                        lbound(mask_3d,3):ubound(mask_3d,3)))
      allocate(got_clean(kb:ke+kh), got_halo(kb:ke+kh), exp(kb:ke+kh), sum_local(kb:ke+kh), sum_global(kb:ke+kh))

      var_clean = -9999.
      var_halo = 7777.
      do k = k1, k2
        do j = j1, j2
          do i = i1, i2
            var_clean(i,j,k) = real(var_loc) + 0.13*real(i) - 0.07*real(j) + 0.011*real(k)
            var_halo(i,j,k) = var_clean(i,j,k)
          end do
        end do
      end do

      got_clean = 0.
      got_halo = 0.
      call avg_xy_fluid(got_clean, var_clean, var_loc, kh, .true.)
      call avg_xy_fluid(got_halo, var_halo, var_loc, kh, .true.)

      do k = kb, ke + kh
        sum_local(k) = sum(var_clean(i1:i2,j1:j2,k1+k-1) * real(mask_3d(ii1:ii2,ij1:ij2,ik1+k-1)))
      end do
      call MPI_ALLREDUCE(sum_local, sum_global, size(sum_local), MY_REAL, MPI_SUM, comm3d, mpierr)

      do k = kb, ke + kh
        if (mask_1d(is1+k-1) == 0) then
          exp(k) = -999.
        else
          exp(k) = sum_global(k) / real(mask_1d(is1+k-1))
        end if
      end do

      check_fluid_loc = compare_real_1d('avg_xy_fluid '//trim(label), got_clean, exp)
      if (.not. compare_real_1d('avg_xy_fluid '//trim(label)//' halo invariance', got_halo, exp)) check_fluid_loc = .false.
      if (.not. compare_real_1d('avg_xy_fluid '//trim(label)//' clean vs halo', got_halo, got_clean)) check_fluid_loc = .false.

      if (.not. check_sum_y_loc(label, var_loc, mask_3d, var_clean, var_halo)) check_fluid_loc = .false.
      if (.not. check_sum_x_loc(label, var_loc, mask_3d, var_clean, var_halo)) check_fluid_loc = .false.

      deallocate(var_clean, var_halo, got_clean, got_halo, exp, sum_local, sum_global)
    end function check_fluid_loc

    logical function check_fluid_y_loc(label, var_loc, mask_3d, mask_2d)
      use modglobal, only : kb, ke
      implicit none
      character(len=*), intent(in) :: label
      integer, intent(in)          :: var_loc
      integer, intent(in)          :: mask_3d(:,:,:), mask_2d(:,:)
      real, allocatable            :: var_clean(:,:,:), var_halo(:,:,:)
      real, allocatable            :: got_clean(:,:), got_halo(:,:), exp(:,:), sum_local(:,:), sum_global(:,:)
      integer                      :: i, j, k
      integer                      :: i1, i2, j1, j2, k1, k2
      integer                      :: ii1, ii2, ij1, ij2, ik1, ik2
      integer                      :: it1, it2, kt1, kt2

      call select_fluid_y_reference_windows(mask_real_compatible(mask_3d), mask_3d, mask_2d, &
           i1, i2, j1, j2, k1, k2, ii1, ii2, ij1, ij2, ik1, ik2, it1, it2, kt1, kt2)

      allocate(var_clean(lbound(mask_3d,1):ubound(mask_3d,1), &
                         lbound(mask_3d,2):ubound(mask_3d,2), &
                         lbound(mask_3d,3):ubound(mask_3d,3)))
      allocate(var_halo(lbound(mask_3d,1):ubound(mask_3d,1), &
                        lbound(mask_3d,2):ubound(mask_3d,2), &
                        lbound(mask_3d,3):ubound(mask_3d,3)))
      allocate(got_clean(1:zsize(1),1:zsize(3)), got_halo(1:zsize(1),1:zsize(3)), exp(1:zsize(1),1:zsize(3)), &
               sum_local(1:zsize(1),1:zsize(3)), sum_global(1:zsize(1),1:zsize(3)))

      var_clean = -4321.
      var_halo = 6543.
      do k = k1, k2
        do j = j1, j2
          do i = i1, i2
            var_clean(i,j,k) = 0.25*real(var_loc) + 0.21*real(i) - 0.03*real(j) + 0.017*real(k)
            var_halo(i,j,k) = var_clean(i,j,k)
          end do
        end do
      end do

      call avg_y_fluid(got_clean, var_clean, var_loc)
      call avg_y_fluid(got_halo, var_halo, var_loc)

      do k = 1, zsize(3)
        do i = 1, zsize(1)
          sum_local(i,k) = sum(var_clean(i1+i-1,j1:j2,k1+k-1) * real(mask_3d(ii1+i-1,ij1:ij2,ik1+k-1)))
        end do
      end do
      call MPI_ALLREDUCE(sum_local, sum_global, size(sum_local), MY_REAL, MPI_SUM, comm3d, mpierr)

      do k = 1, zsize(3)
        do i = 1, zsize(1)
          if (mask_2d(it1+i-1,kt1+k-1) == 0) then
            exp(i,k) = -999.
          else
            exp(i,k) = sum_global(i,k) / real(mask_2d(it1+i-1,kt1+k-1))
          end if
        end do
      end do

      check_fluid_y_loc = compare_real_2d('avg_y_fluid '//trim(label), got_clean, exp)
      if (.not. compare_real_2d('avg_y_fluid '//trim(label)//' halo invariance', got_halo, exp)) check_fluid_y_loc = .false.
      if (.not. compare_real_2d('avg_y_fluid '//trim(label)//' clean vs halo', got_halo, got_clean)) check_fluid_y_loc = .false.

      deallocate(var_clean, var_halo, got_clean, got_halo, exp, sum_local, sum_global)
    end function check_fluid_y_loc

    logical function check_sum_y_loc(label, var_loc, mask_3d, var_clean, var_halo)
      implicit none
      character(len=*), intent(in) :: label
      integer, intent(in)          :: var_loc
      integer, intent(in)          :: mask_3d(:,:,:)
      real, intent(in)             :: var_clean(:,:,:), var_halo(:,:,:)
      real, allocatable            :: got_clean(:,:), got_halo(:,:), exp_local(:,:), exp(:,:)
      integer                      :: i, k
      integer                      :: i1, i2, j1, j2, k1, k2
      integer                      :: ii1, ii2, ij1, ij2, ik1, ik2

      allocate(got_clean(1:zsize(1),1:zsize(3)), got_halo(1:zsize(1),1:zsize(3)), &
               exp_local(1:zsize(1),1:zsize(3)), exp(1:zsize(1),1:zsize(3)))

      call sum_y_fluid(got_clean, var_clean, var_loc)
      call sum_y_fluid(got_halo, var_halo, var_loc)

      call select_fluid_sum_reference_windows(var_clean, mask_3d, &
           i1, i2, j1, j2, k1, k2, ii1, ii2, ij1, ij2, ik1, ik2)

      do k = 1, zsize(3)
        do i = 1, zsize(1)
          exp_local(i,k) = sum(var_clean(i1+i-1,j1:j2,k1+k-1) * real(mask_3d(ii1+i-1,ij1:ij2,ik1+k-1)))
        end do
      end do
      call MPI_ALLREDUCE(exp_local, exp, size(exp_local), MY_REAL, MPI_SUM, comm3d, mpierr)

      check_sum_y_loc = compare_real_2d('sum_y_fluid '//trim(label), got_clean, exp)
      if (.not. compare_real_2d('sum_y_fluid '//trim(label)//' halo invariance', got_halo, exp)) check_sum_y_loc = .false.
      if (.not. compare_real_2d('sum_y_fluid '//trim(label)//' clean vs halo', got_halo, got_clean)) check_sum_y_loc = .false.

      deallocate(got_clean, got_halo, exp_local, exp)
    end function check_sum_y_loc

    logical function check_sum_x_loc(label, var_loc, mask_3d, var_clean, var_halo)
      implicit none
      character(len=*), intent(in) :: label
      integer, intent(in)          :: var_loc
      integer, intent(in)          :: mask_3d(:,:,:)
      real, intent(in)             :: var_clean(:,:,:), var_halo(:,:,:)
      real, allocatable            :: got_clean(:,:), got_halo(:,:), exp_local(:,:), exp(:,:)
      integer                      :: j, k
      integer                      :: i1, i2, j1, j2, k1, k2
      integer                      :: ii1, ii2, ij1, ij2, ik1, ik2

      allocate(got_clean(1:zsize(2),1:zsize(3)), got_halo(1:zsize(2),1:zsize(3)), &
               exp_local(1:zsize(2),1:zsize(3)), exp(1:zsize(2),1:zsize(3)))

      call sum_x_fluid(got_clean, var_clean, var_loc)
      call sum_x_fluid(got_halo, var_halo, var_loc)

      call select_fluid_sum_reference_windows(var_clean, mask_3d, &
           i1, i2, j1, j2, k1, k2, ii1, ii2, ij1, ij2, ik1, ik2)

      do k = 1, zsize(3)
        do j = 1, zsize(2)
          exp_local(j,k) = sum(var_clean(i1:i2,j1+j-1,k1+k-1) * real(mask_3d(ii1:ii2,ij1+j-1,ik1+k-1)))
        end do
      end do
      call MPI_ALLREDUCE(exp_local, exp, size(exp_local), MY_REAL, MPI_SUM, comm3d, mpierr)

      check_sum_x_loc = compare_real_2d('sum_x_fluid '//trim(label), got_clean, exp)
      if (.not. compare_real_2d('sum_x_fluid '//trim(label)//' halo invariance', got_halo, exp)) check_sum_x_loc = .false.
      if (.not. compare_real_2d('sum_x_fluid '//trim(label)//' clean vs halo', got_halo, got_clean)) check_sum_x_loc = .false.

      deallocate(got_clean, got_halo, exp_local, exp)
    end function check_sum_x_loc

    logical function compare_real_1d(label, got, exp)
      implicit none
      character(len=*), intent(in) :: label
      real, intent(in)             :: got(:), exp(:)
      real                         :: max_abs
      integer                      :: imax(1)

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
      real, intent(in)             :: got(:,:), exp(:,:)
      real                         :: max_abs
      integer                      :: imax(2)

      max_abs = maxval(abs(got - exp))
      compare_real_2d = max_abs <= 1.e-9
      if ((.not. compare_real_2d) .and. myid == 0) then
        imax = maxloc(abs(got - exp))
        write(*,'(A,1X,A,1X,ES12.4,1X,A,I0,A,I0,1X,A,ES12.4,1X,A,ES12.4)') &
             'FAIL', trim(label), max_abs, 'idx', imax(1), ',', imax(2), 'got', got(imax(1),imax(2)), 'exp', exp(imax(1),imax(2))
      end if
    end function compare_real_2d

    logical function nearly_equal_scalar(got, exp)
      implicit none
      real, intent(in) :: got, exp
      nearly_equal_scalar = abs(got - exp) <= 1.e-9
    end function nearly_equal_scalar

    subroutine report_scalar_failure(label, got, exp)
      implicit none
      character(len=*), intent(in) :: label
      real, intent(in)             :: got, exp

      if (myid == 0) then
        write(*,'(A,1X,A,1X,ES12.4,1X,ES12.4)') 'FAIL', trim(label), got, exp
      end if
    end subroutine report_scalar_failure

    function mask_real_compatible(mask_3d) result(var)
      implicit none
      integer, intent(in) :: mask_3d(:,:,:)
      real, allocatable   :: var(:,:,:)

      allocate(var(lbound(mask_3d,1):ubound(mask_3d,1), &
                   lbound(mask_3d,2):ubound(mask_3d,2), &
                   lbound(mask_3d,3):ubound(mask_3d,3)))
      var = 0.
    end function mask_real_compatible

    subroutine select_fluid_xy_reference_windows(var, mask_3d, mask_1d, interior_z, lower_halo_z, &
         i1, i2, j1, j2, k1, k2, ii1, ii2, ij1, ij2, ik1, ik2, is1, is2)
      implicit none
      real, intent(in)    :: var(:,:,:)
      integer, intent(in) :: mask_3d(:,:,:), mask_1d(:)
      integer, intent(in) :: interior_z, lower_halo_z
      integer, intent(out):: i1, i2, j1, j2, k1, k2
      integer, intent(out):: ii1, ii2, ij1, ij2, ik1, ik2
      integer, intent(out):: is1, is2

      call reference_bounds_real3d(var, zsize(1), zsize(2), interior_z, lower_halo_z, 0, i1, i2, j1, j2, k1, k2)
      call reference_bounds_int3d(mask_3d, zsize(1), zsize(2), interior_z, lower_halo_z, 0, ii1, ii2, ij1, ij2, ik1, ik2)
      call reference_bounds_int1d(mask_1d, interior_z, lower_halo_z, 0, is1, is2)
    end subroutine select_fluid_xy_reference_windows

    subroutine select_fluid_y_reference_windows(var, mask_3d, mask_2d, &
         i1, i2, j1, j2, k1, k2, ii1, ii2, ij1, ij2, ik1, ik2, it1, it2, kt1, kt2)
      implicit none
      real, intent(in)    :: var(:,:,:)
      integer, intent(in) :: mask_3d(:,:,:), mask_2d(:,:)
      integer, intent(out):: i1, i2, j1, j2, k1, k2
      integer, intent(out):: ii1, ii2, ij1, ij2, ik1, ik2
      integer, intent(out):: it1, it2, kt1, kt2

      call reference_bounds_real3d(var, zsize(1), zsize(2), zsize(3), 0, 1, i1, i2, j1, j2, k1, k2)
      call reference_bounds_int3d(mask_3d, zsize(1), zsize(2), zsize(3), 0, 1, ii1, ii2, ij1, ij2, ik1, ik2)
      call reference_bounds_int2d_xz(mask_2d, zsize(1), zsize(3), it1, it2, kt1, kt2)
    end subroutine select_fluid_y_reference_windows

    subroutine select_fluid_sum_reference_windows(var, mask_3d, &
         i1, i2, j1, j2, k1, k2, ii1, ii2, ij1, ij2, ik1, ik2)
      implicit none
      real, intent(in)    :: var(:,:,:)
      integer, intent(in) :: mask_3d(:,:,:)
      integer, intent(out):: i1, i2, j1, j2, k1, k2
      integer, intent(out):: ii1, ii2, ij1, ij2, ik1, ik2

      call reference_bounds_real3d(var, zsize(1), zsize(2), zsize(3), 0, 1, i1, i2, j1, j2, k1, k2)
      call reference_bounds_int3d(mask_3d, zsize(1), zsize(2), zsize(3), 0, 1, ii1, ii2, ij1, ij2, ik1, ik2)
    end subroutine select_fluid_sum_reference_windows

    subroutine reference_bounds_real3d(arr, interior_x, interior_y, interior_z, lower_halo_z, upper_halo_z, &
         i1, i2, j1, j2, k1, k2)
      implicit none
      real, intent(in)    :: arr(:,:,:)
      integer, intent(in) :: interior_x, interior_y, interior_z, lower_halo_z, upper_halo_z
      integer, intent(out):: i1, i2, j1, j2, k1, k2

      call reference_bounds_xy(size(arr,1), interior_x, i1, i2)
      call reference_bounds_xy(size(arr,2), interior_y, j1, j2)
      call reference_bounds_z(size(arr,3), interior_z, lower_halo_z, upper_halo_z, k1, k2)
    end subroutine reference_bounds_real3d

    subroutine reference_bounds_int3d(arr, interior_x, interior_y, interior_z, lower_halo_z, upper_halo_z, &
         i1, i2, j1, j2, k1, k2)
      implicit none
      integer, intent(in) :: arr(:,:,:)
      integer, intent(in) :: interior_x, interior_y, interior_z, lower_halo_z, upper_halo_z
      integer, intent(out):: i1, i2, j1, j2, k1, k2

      call reference_bounds_xy(size(arr,1), interior_x, i1, i2)
      call reference_bounds_xy(size(arr,2), interior_y, j1, j2)
      call reference_bounds_z(size(arr,3), interior_z, lower_halo_z, upper_halo_z, k1, k2)
    end subroutine reference_bounds_int3d

    subroutine reference_bounds_int2d_xz(arr, interior_x, interior_z, i1, i2, k1, k2)
      implicit none
      integer, intent(in) :: arr(:,:)
      integer, intent(in) :: interior_x, interior_z
      integer, intent(out):: i1, i2, k1, k2

      call reference_bounds_xy(size(arr,1), interior_x, i1, i2)
      call reference_bounds_z(size(arr,2), interior_z, 0, 0, k1, k2)
    end subroutine reference_bounds_int2d_xz

    subroutine reference_bounds_int1d(arr, interior_z, lower_halo_z, upper_halo_z, k1, k2)
      implicit none
      integer, intent(in) :: arr(:)
      integer, intent(in) :: interior_z, lower_halo_z, upper_halo_z
      integer, intent(out):: k1, k2

      call reference_bounds_z(size(arr), interior_z, lower_halo_z, upper_halo_z, k1, k2)
    end subroutine reference_bounds_int1d

    subroutine reference_bounds_xy(actual_size, interior_size, start_idx, end_idx)
      implicit none
      integer, intent(in)  :: actual_size, interior_size
      integer, intent(out) :: start_idx, end_idx
      integer              :: halo_total

      halo_total = actual_size - interior_size
      if (halo_total < 0 .or. mod(halo_total, 2) /= 0) stop 1
      start_idx = halo_total / 2 + 1
      end_idx = start_idx + interior_size - 1
    end subroutine reference_bounds_xy

    subroutine reference_bounds_z(actual_size, interior_size, lower_halo, upper_halo, start_idx, end_idx)
      implicit none
      integer, intent(in)  :: actual_size, interior_size, lower_halo, upper_halo
      integer, intent(out) :: start_idx, end_idx

      if (actual_size == interior_size) then
        start_idx = 1
      elseif (actual_size == interior_size + lower_halo + upper_halo) then
        start_idx = lower_halo + 1
      elseif (actual_size == interior_size + lower_halo) then
        start_idx = lower_halo + 1
      elseif (actual_size == interior_size + upper_halo) then
        start_idx = 1
      else
        stop 1
      end if

      end_idx = start_idx + interior_size - 1
    end subroutine reference_bounds_z

    subroutine load_solid_info_from_sparse(filename, npts, solid_info)
      implicit none
      character(len=*), intent(in)   :: filename
      integer, intent(in)            :: npts
      type(solid_info_type), intent(inout) :: solid_info

      integer, allocatable :: ids_loc(:), pts_loc(:,:)

      call read_sparse_ijk(filename, npts, solid_info%nsolptsrank, ids_loc, pts_loc)
      solid_info%nsolpts = npts
      call move_alloc(ids_loc, solid_info%solptsrank)
      call move_alloc(pts_loc, solid_info%solpts_loc)
    end subroutine load_solid_info_from_sparse

  end function tests_mpi_fieldops

end module tests
