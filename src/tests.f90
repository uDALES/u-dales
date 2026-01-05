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
  public :: tests_read_sparse_ijk, tests_2decomp_init_exit

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
  subroutine tests_read_sparse_ijk()
    use modglobal,    only : libm, cexpnr
    use modmpi,       only : myid
    use modfileinput, only : read_sparse_ijk
    use modibm,       only : initibm
    use modibm,       only : solid_info_u, solid_info_v, solid_info_w, solid_info_c
    use modibm,       only : bound_info_u, bound_info_v, bound_info_w, bound_info_c
    use modibm,       only : nsolpts_u, nsolpts_v, nsolpts_w, nsolpts_c
    use modibm,       only : nbndpts_u, nbndpts_v, nbndpts_w, nbndpts_c
    use initfac,      only : readfacetfiles
    use decomp_2d,    only : zstart, zend
    
    implicit none
    
    integer :: npts_loc_new
    integer, allocatable :: ids_loc_new(:)
    integer, allocatable :: pts_loc_new(:,:)
    integer :: m
    logical :: all_passed
    
    if (.not. libm) then
      if (myid == 0) write(*, '(A)') 'Skipping read_sparse_ijk test: libm=false'
      return
    end if
    
    ! Read facet files and initialize IBM - populates all global arrays
    call readfacetfiles
    call initibm
    
    all_passed = .true.
    
    ! Test solid_u
    call read_sparse_ijk('solid_u.' // cexpnr, nsolpts_u, npts_loc_new, ids_loc_new, pts_loc_new)
    if (.not. compare_solid(solid_info_u, npts_loc_new, ids_loc_new, pts_loc_new, 'solid_u')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test solid_v
    call read_sparse_ijk('solid_v.' // cexpnr, nsolpts_v, npts_loc_new, ids_loc_new, pts_loc_new)
    if (.not. compare_solid(solid_info_v, npts_loc_new, ids_loc_new, pts_loc_new, 'solid_v')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test solid_w
    call read_sparse_ijk('solid_w.' // cexpnr, nsolpts_w, npts_loc_new, ids_loc_new, pts_loc_new)
    if (.not. compare_solid(solid_info_w, npts_loc_new, ids_loc_new, pts_loc_new, 'solid_w')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test solid_c
    call read_sparse_ijk('solid_c.' // cexpnr, nsolpts_c, npts_loc_new, ids_loc_new, pts_loc_new)
    if (.not. compare_solid(solid_info_c, npts_loc_new, ids_loc_new, pts_loc_new, 'solid_c')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test fluid_boundary_u
    call read_sparse_ijk('fluid_boundary_u.' // cexpnr, nbndpts_u, npts_loc_new, ids_loc_new, pts_loc_new)
    if (.not. compare_boundary(bound_info_u, npts_loc_new, ids_loc_new, pts_loc_new, 'fluid_boundary_u')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test fluid_boundary_v
    call read_sparse_ijk('fluid_boundary_v.' // cexpnr, nbndpts_v, npts_loc_new, ids_loc_new, pts_loc_new)
    if (.not. compare_boundary(bound_info_v, npts_loc_new, ids_loc_new, pts_loc_new, 'fluid_boundary_v')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test fluid_boundary_w
    call read_sparse_ijk('fluid_boundary_w.' // cexpnr, nbndpts_w, npts_loc_new, ids_loc_new, pts_loc_new)
    if (.not. compare_boundary(bound_info_w, npts_loc_new, ids_loc_new, pts_loc_new, 'fluid_boundary_w')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    ! Test fluid_boundary_c
    call read_sparse_ijk('fluid_boundary_c.' // cexpnr, nbndpts_c, npts_loc_new, ids_loc_new, pts_loc_new)
    if (.not. compare_boundary(bound_info_c, npts_loc_new, ids_loc_new, pts_loc_new, 'fluid_boundary_c')) all_passed = .false.
    deallocate(ids_loc_new, pts_loc_new)
    
    if (all_passed .and. myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A)') 'ALL TESTS PASSED: read_sparse_ijk'
      write(*, '(A)') '  Tested 8 files successfully'
      write(*, '(A)') '  All results match IBM initialization code'
      write(*, '(A)') '================================================'
    end if
    
  end subroutine tests_read_sparse_ijk
  
  function compare_solid(solid_info, npts_loc_new, ids_loc_new, pts_loc_new, label) result(passed)
    use modmpi, only : myid
    use modibmdata, only : solid_info_type
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
    use modibmdata, only : bound_info_type
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

end module tests
