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
  use modmpi, only : myid, comm3d, mpierr, my_real, avexy_ibm, avey_ibm, sumx_ibm, sumy_ibm
  
  implicit none
  save
  public :: tests_read_sparse_ijk, tests_2decomp_init_exit, tests_mpi_operators

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
    use modibm,       only : initibm
    use modibm,       only : solid_info_u, solid_info_v, solid_info_w, solid_info_c
    use modibm,       only : bound_info_u, bound_info_v, bound_info_w, bound_info_c
    use modibm,       only : nsolpts_u, nsolpts_v, nsolpts_w, nsolpts_c
    use modibm,       only : nbndpts_u, nbndpts_v, nbndpts_w, nbndpts_c
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
    use modibm, only : solid_info_type
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
    use modibm, only : bound_info_type
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
      call avexy_ibm(got, var_clean, ib, ie, jb, je, kb, ke, 0, 0, khc, mask_3d(ib:ie,jb:je,kb:ke+khc), mask_1d, .true.)

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

end module tests
