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
  use MPI
  use decomp_2d
  use modmpi, only : myid
  
  implicit none
  save
  public :: tests_read_sparse_ijk, tests_2decomp_init_exit, tests_trees_sparse_compare

contains

  !> Initialize MPI and 2DECOMP for testing
  subroutine tests_2decomp_init_exit
    integer :: nx=64, ny=64, nz=64
    integer :: p_row=0, p_col=0
    integer :: ierror

    call MPI_INIT(ierror)
    call decomp_2d_init(nx,ny,nz,p_row,p_col)

    write(*,*) xstart
    write(*,*) ystart
    write(*,*) zstart
    write(*,*) xend
    write(*,*) yend
    write(*,*) zend
    write(*,*) xsize
    write(*,*) ysize
    write(*,*) zsize

    call decomp_2d_finalize
    call MPI_FINALIZE(ierror)

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
      if (pts_loc_new(m,1) /= solid_info%solpts_loc(m,1) - zstart(1) + 1 .or. &
          pts_loc_new(m,2) /= solid_info%solpts_loc(m,2) - zstart(2) + 1 .or. &
          pts_loc_new(m,3) /= solid_info%solpts_loc(m,3) - zstart(3) + 1) then
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
      if (pts_loc_new(m,1) /= bound_info%bndpts_loc(m,1) - zstart(1) + 1 .or. &
          pts_loc_new(m,2) /= bound_info%bndpts_loc(m,2) - zstart(2) + 1 .or. &
          pts_loc_new(m,3) /= bound_info%bndpts_loc(m,3) - zstart(3) + 1) then
        write(*, '(A,I0,A,A,A)') 'FAIL on rank ', myid, ': ', trim(label), ' coordinate mismatch'
        passed = .false.
        stop 1
      end if
    end do
    
  end function compare_boundary

  !> Compare block-based tree forcing with sparse vegetation masks
  logical function tests_trees_sparse_compare()
    use modglobal, only : ib, ie, jb, je, kb, ke, itot, jtot, ntree_max, cd, cexpnr, ltrees
    use modfields, only : um, vm, wm, tr_u, tr_v, tr_w, up, vp, wp, ladzh, ladzf
    use modtrees,  only : trees, createtrees, trees_sparse
    use modmpi,    only : myid, myidx, myidy, comm3d, mpierr, my_real
    use decomp_2d, only : xstart, ystart, zstart
    implicit none

    integer :: i, j, k
    real, allocatable :: block_u(:,:,:), block_v(:,:,:), block_w(:,:,:)
    real :: diff_local, diff_global
    real :: diff_u_local, diff_v_local, diff_w_local
    real :: diff_u_global, diff_v_global, diff_w_global
    real(kind=8) :: pair_u(2), pair_v(2), pair_w(2)
    real(kind=8) :: pair_u_g(2), pair_v_g(2), pair_w_g(2)
    integer :: rank_u, rank_v, rank_w
    integer :: idx_u(3), idx_v(3), idx_w(3)
    integer :: idx_u_root(3), idx_v_root(3), idx_w_root(3)

    tests_trees_sparse_compare = .false.

    if (.not. ltrees) then
      if (myid == 0) write(*,*) 'ltrees is false; skipping tests_trees_sparse_compare'
      return
    end if

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A)') 'tests_trees_sparse_compare: SPARSE TREE INPUT TEST'
      write(*, '(A)') '------------------------------------------------'
    end if

    ! Ensure tree properties are initialized (also reads veg_lad)
    call createtrees

    ! Build a simple deterministic velocity field so drag is non-zero
    um = 0.0
    vm = 0.0
    wm = 0.0
    do k = kb, ke
      do j = jb, je
        do i = ib, ie
          um(i,j,k) = real(0.01, kind=kind(um)) * real(i, kind=kind(um)) + real(0.002, kind=kind(um)) * real(j, kind=kind(um)) + real(0.0001, kind=kind(um)) * real(k, kind=kind(um))
          vm(i,j,k) = real(0.015, kind=kind(vm)) * real(i, kind=kind(vm)) + real(0.003, kind=kind(vm)) * real(j, kind=kind(vm)) + real(0.0002, kind=kind(vm)) * real(k, kind=kind(vm))
          wm(i,j,k) = real(0.02, kind=kind(wm)) * real(i, kind=kind(wm)) + real(0.004, kind=kind(wm)) * real(j, kind=kind(wm)) + real(0.0003, kind=kind(wm)) * real(k, kind=kind(wm))
        end do
      end do
    end do

    tr_u = 0.0
    tr_v = 0.0
    tr_w = 0.0
    up = 0.0
    vp = 0.0
    wp = 0.0

    call trees

    allocate(block_u(ib:ie,jb:je,kb:ke))
    allocate(block_v(ib:ie,jb:je,kb:ke))
    allocate(block_w(ib:ie,jb:je,kb:ke))
    block_u = tr_u
    block_v = tr_v
    block_w = tr_w
    tr_u = 0.0
    tr_v = 0.0
    tr_w = 0.0

    ! Apply sparse vegetation forcing via modtrees routine
    call trees_sparse
    diff_local = 0.0
    diff_u_local = 0.0
    diff_v_local = 0.0
    diff_w_local = 0.0
    idx_u = (/ ib, jb, kb /)
    idx_v = (/ ib, jb, kb /)
    idx_w = (/ ib, jb, kb /)
    do k = kb, ke
      do j = jb, je
        do i = ib, ie
          diff_u_local = max(diff_u_local, abs(tr_u(i,j,k) - block_u(i,j,k)))
          diff_v_local = max(diff_v_local, abs(tr_v(i,j,k) - block_v(i,j,k)))
          diff_w_local = max(diff_w_local, abs(tr_w(i,j,k) - block_w(i,j,k)))
          if (abs(tr_u(i,j,k) - block_u(i,j,k)) >= diff_u_local) then
            idx_u = (/ i, j, k /)
          end if
          if (abs(tr_v(i,j,k) - block_v(i,j,k)) >= diff_v_local) then
            idx_v = (/ i, j, k /)
          end if
          if (abs(tr_w(i,j,k) - block_w(i,j,k)) >= diff_w_local) then
            idx_w = (/ i, j, k /)
          end if
          diff_local = max(diff_local, diff_u_local)
          diff_local = max(diff_local, diff_v_local)
          diff_local = max(diff_local, diff_w_local)
        end do
      end do
    end do

    call MPI_ALLREDUCE(diff_local, diff_global, 1, MY_REAL, MPI_MAX, comm3d, mpierr)

    pair_u(1) = real(diff_u_local, kind=8)
    pair_u(2) = real(myid, kind=8)
    pair_v(1) = real(diff_v_local, kind=8)
    pair_v(2) = real(myid, kind=8)
    pair_w(1) = real(diff_w_local, kind=8)
    pair_w(2) = real(myid, kind=8)

    call MPI_REDUCE(pair_u, pair_u_g, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, 0, comm3d, mpierr)
    call MPI_REDUCE(pair_v, pair_v_g, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, 0, comm3d, mpierr)
    call MPI_REDUCE(pair_w, pair_w_g, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, 0, comm3d, mpierr)

    if (myid == 0) then
      rank_u = nint(pair_u_g(2))
      rank_v = nint(pair_v_g(2))
      rank_w = nint(pair_w_g(2))
    end if

    call MPI_BCAST(rank_u, 1, MPI_INTEGER, 0, comm3d, mpierr)
    call MPI_BCAST(rank_v, 1, MPI_INTEGER, 0, comm3d, mpierr)
    call MPI_BCAST(rank_w, 1, MPI_INTEGER, 0, comm3d, mpierr)

    if (myid == rank_u .and. myid /= 0) call MPI_SEND(idx_u, 3, MPI_INTEGER, 0, 101, comm3d, mpierr)
    if (myid == rank_v .and. myid /= 0) call MPI_SEND(idx_v, 3, MPI_INTEGER, 0, 102, comm3d, mpierr)
    if (myid == rank_w .and. myid /= 0) call MPI_SEND(idx_w, 3, MPI_INTEGER, 0, 103, comm3d, mpierr)

    if (myid == 0) then
      if (rank_u == 0) then
        idx_u_root = idx_u
      else
        call MPI_RECV(idx_u_root, 3, MPI_INTEGER, rank_u, 101, comm3d, MPI_STATUS_IGNORE, mpierr)
      end if
      if (rank_v == 0) then
        idx_v_root = idx_v
      else
        call MPI_RECV(idx_v_root, 3, MPI_INTEGER, rank_v, 102, comm3d, MPI_STATUS_IGNORE, mpierr)
      end if
      if (rank_w == 0) then
        idx_w_root = idx_w
      else
        call MPI_RECV(idx_w_root, 3, MPI_INTEGER, rank_w, 103, comm3d, MPI_STATUS_IGNORE, mpierr)
      end if
      diff_u_global = pair_u_g(1)
      diff_v_global = pair_v_g(1)
      diff_w_global = pair_w_g(1)
    end if

    call MPI_BCAST(diff_u_global, 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(diff_v_global, 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(diff_w_global, 1, MY_REAL, 0, comm3d, mpierr)

    if (myid == 0) then
      write(*,*) 'Max diff (sparse vs block) tr_u = ', diff_u_global, &
             ' at local (', idx_u_root(1), ',', idx_u_root(2), ',', idx_u_root(3), &
             ') rank ', rank_u, ' (myidx,myidy)=(', myidx, ',', myidy, ')'
      write(*,*) 'Max diff (sparse vs block) tr_v = ', diff_v_global, &
             ' at local (', idx_v_root(1), ',', idx_v_root(2), ',', idx_v_root(3), &
             ') rank ', rank_v, ' (myidx,myidy)=(', myidx, ',', myidy, ')'
      write(*,*) 'Max diff (sparse vs block) tr_w = ', diff_w_global, &
             ' at local (', idx_w_root(1), ',', idx_w_root(2), ',', idx_w_root(3), &
             ') rank ', rank_w, ' (myidx,myidy)=(', myidx, ',', myidy, ')'
      if (diff_global < 1.0e-12) then
        write(*, '(A)') '------------------------------------------------'
        write(*, '(A)') 'ALL TESTS PASSED: tests_trees_sparse_compare'
        write(*, '(A)') '  Sparse and block tree forcing match'
        write(*, '(A)') '================================================'
      else
        write(*, '(A)') '------------------------------------------------'
        write(*, '(A)') 'TESTS FAILED: tests_trees_sparse_compare'
        write(*, '(A)') '  Sparse and block tree forcing differ'
        write(*, '(A)') '================================================'
      end if
    end if

    tests_trees_sparse_compare = diff_global < 1.0e-12

    deallocate(block_u, block_v, block_w)

  end function tests_trees_sparse_compare

end module tests
