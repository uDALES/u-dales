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
  
  implicit none
  save
  public :: test_read_sparse_ijk, inittest, exittest

contains

  !> Initialize MPI and 2DECOMP for testing
  subroutine tests_2decomp_init_exit
    integer :: nx=64, ny=64, nz=64
    !integer :: p_row=2, p_col=2
    integer :: p_row=0, p_col=0
    integer :: ierror

    call MPI_INIT(ierror)
    !call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierr)
    !call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    !call decomp_2d_init(nx,ny,nz,p_row,p_col)
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

contains

  !> Test the read_sparse_ijk routine by comparing with existing manual reading
  subroutine test_read_sparse_ijk(filename, npts, nskip)
    use modglobal,    only : ifinput
    use modmpi,       only : myid, comm3d, mpierr
    use modfileinput, only : read_sparse_ijk
    use decomp_2d,    only : zstart, zend
    use mpi
    
    implicit none
    
    character(len=*), intent(in) :: filename
    integer, intent(in)           :: npts
    integer, intent(in), optional :: nskip
    
    ! Variables for old method
    integer, allocatable :: pts_glob_old(:,:)
    logical, allocatable :: lpts_old(:)
    integer, allocatable :: ids_loc_old(:)
    integer, allocatable :: pts_loc_old(:,:)
    integer :: npts_loc_old
    
    ! Variables for new method
    integer :: npts_loc_new
    integer, allocatable :: ids_loc_new(:)
    integer, allocatable :: pts_loc_new(:,:)
    
    ! Temporary variables
    integer :: n, m, nh, ierr, i, j, k
    character(80) :: chmess
    logical :: test_passed
    
    ! Set number of header lines
    nh = 1
    if (present(nskip)) nh = nskip
    
    !---------------------------------------------------------------------------
    ! OLD METHOD: Manual reading (existing code pattern)
    !---------------------------------------------------------------------------
    allocate(pts_glob_old(npts, 3))
    allocate(lpts_old(npts))
    
    ! Read file on rank 0
    if (myid == 0) then
      open(ifinput, file=filename, status='old', iostat=ierr)
      if (ierr /= 0) then
        write(*, '(A,A)') 'ERROR: Cannot open file: ', trim(filename)
        stop 1
      end if
      
      ! Skip header lines
      do n = 1, nh
        read(ifinput, '(a80)') chmess
      end do
      
      ! Read coordinates
      do n = 1, npts
        read(ifinput, *) pts_glob_old(n, 1), pts_glob_old(n, 2), pts_glob_old(n, 3)
      end do
      
      close(ifinput)
    end if
    
    ! Broadcast to all ranks
    call MPI_BCAST(pts_glob_old, npts*3, MPI_INTEGER, 0, comm3d, mpierr)
    
    ! Determine local points and count them
    npts_loc_old = 0
    do n = 1, npts
      if ((pts_glob_old(n,1) >= zstart(1) .and. pts_glob_old(n,1) <= zend(1)) .and. &
          (pts_glob_old(n,2) >= zstart(2) .and. pts_glob_old(n,2) <= zend(2))) then
        lpts_old(n) = .true.
        npts_loc_old = npts_loc_old + 1
      else
        lpts_old(n) = .false.
      end if
    end do
    
    ! Pack into compact arrays with index remapping
    allocate(ids_loc_old(npts_loc_old))
    allocate(pts_loc_old(npts_loc_old, 3))
    
    m = 0
    do n = 1, npts
      if (lpts_old(n)) then
        m = m + 1
        ids_loc_old(m) = n
        pts_loc_old(m,1) = pts_glob_old(n,1) - zstart(1) + 1
        pts_loc_old(m,2) = pts_glob_old(n,2) - zstart(2) + 1
        pts_loc_old(m,3) = pts_glob_old(n,3) - zstart(3) + 1
      end if
    end do
    
    deallocate(pts_glob_old, lpts_old)
    
    !---------------------------------------------------------------------------
    ! NEW METHOD: Using read_sparse_ijk
    !---------------------------------------------------------------------------
    if (present(nskip)) then
      call read_sparse_ijk(filename, npts, npts_loc_new, ids_loc_new, pts_loc_new, nskip=nh)
    else
      call read_sparse_ijk(filename, npts, npts_loc_new, ids_loc_new, pts_loc_new)
    end if
    
    !---------------------------------------------------------------------------
    ! COMPARISON
    !---------------------------------------------------------------------------
    test_passed = .true.
    
    ! Check 1: Number of local points
    if (npts_loc_new /= npts_loc_old) then
      write(*, '(A,I0,A)') 'FAIL on rank ', myid, ': Local point count mismatch'
      write(*, '(A,I0)') '  Old method: ', npts_loc_old
      write(*, '(A,I0)') '  New method: ', npts_loc_new
      test_passed = .false.
    end if
    
    ! Check 2: Global indices
    if (test_passed) then
      do m = 1, npts_loc_old
        if (ids_loc_new(m) /= ids_loc_old(m)) then
          write(*, '(A,I0,A,I0)') 'FAIL on rank ', myid, ': Index mismatch at local position ', m
          write(*, '(A,I0)') '  Old: ', ids_loc_old(m)
          write(*, '(A,I0)') '  New: ', ids_loc_new(m)
          test_passed = .false.
          exit
        end if
      end do
    end if
    
    ! Check 3: Local coordinates (remapped)
    if (test_passed) then
      do m = 1, npts_loc_old
        if (pts_loc_new(m,1) /= pts_loc_old(m,1) .or. &
            pts_loc_new(m,2) /= pts_loc_old(m,2) .or. &
            pts_loc_new(m,3) /= pts_loc_old(m,3)) then
          write(*, '(A,I0,A,I0)') 'FAIL on rank ', myid, ': Coordinate mismatch at local position ', m
          write(*, '(A,3I5)') '  Old: ', pts_loc_old(m,:)
          write(*, '(A,3I5)') '  New: ', pts_loc_new(m,:)
          test_passed = .false.
          exit
        end if
      end do
    end if
    
    !---------------------------------------------------------------------------
    ! REPORT RESULTS
    !---------------------------------------------------------------------------
    if (test_passed) then
      if (myid == 0) then
        write(*, '(A)') '================================================'
        write(*, '(A,A)') 'TEST PASSED: read_sparse_ijk for ', trim(filename)
        write(*, '(A,I0,A)') '  Total points in file: ', npts
        write(*, '(A)') '  All ranks verified successfully'
        write(*, '(A)') '================================================'
      end if
    else
      write(*, '(A,A)') 'TEST FAILED for ', trim(filename)
      stop 1
    end if
    
    ! Cleanup
    deallocate(ids_loc_old, pts_loc_old)
    deallocate(ids_loc_new, pts_loc_new)
    
  end subroutine test_read_sparse_ijk

end module tests
