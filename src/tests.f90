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
  public :: tests_read_sparse_ijk, tests_2decomp_init_exit, tests_mpi_operators, tests_basestate

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

  !> Check the derived hydrostatic base state (modbasestate, issue #302) against
  !! closed-form results.
  !!
  !! The profiles are built here rather than read from a case, so this checks the
  !! arithmetic itself and not merely that a run reproduces its own previous
  !! output. It runs at the runmode dispatch point, after initglobal has set up
  !! the vertical grid but before readinitfiles would call initbasestate.
  !!
  !! For a column of uniform thv the discrete integration in initbasestate
  !! telescopes, because dzh(k) = zf(k)-zf(k-1) and dzf(k) = zh(k+1)-zh(k), and
  !! every level then has a closed form:
  !!
  !!   p(z)**(Rd/cp) = ps**(Rd/cp) - g*pref0**(Rd/cp)*z/(cp*thv)
  !!
  !! evaluated at zf for full levels and zh for half levels. That is an
  !! independent expectation, so a wrong hydrostatic scheme fails here even
  !! though it would reproduce itself perfectly in a regression comparison.
  logical function tests_basestate()
    use modglobal,    only : kb, ke, kh, zf, zh, dzf, dzh, grav, cp, rd, rv, pref0
    use modbasestate, only : initbasestate, exitbasestate, ps, thv_b, pf_b, ph_b, &
                             exnf_b, exnh_b

    implicit none

    real, allocatable :: thlprof(:), qtprof(:), thv_expect(:)
    real    :: rdocp, thv0, gamma, qt0, dqtdz, tol_exact, tol_trunc
    integer :: k
    logical :: ok

    ok        = .true.
    tol_exact = 1.e-10  ! vs the telescoped sum: same arithmetic, different order
    tol_trunc = 1.e-5   ! vs the continuous solution: the scheme is 2nd order in dz
    rdocp     = rd/cp

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A)') 'tests_basestate: hydrostatic base state (#302)'
      write(*, '(A)') '================================================'
    end if

    allocate(thlprof(kb:ke), qtprof(kb:ke), thv_expect(kb:ke+kh))

    ! initbasestate integrates upward from ps at z = 0, so every expectation
    ! below is only the right answer if the domain bottom really is at z = 0.
    if (abs(zh(kb)) > 1.e-12) then
      if (myid == 0) write(*, '(A,E12.4)') '  FAIL: zh(kb) is not 0:', zh(kb)
      ok = .false.
    end if

    ! ---- A. dry, uniform thv ------------------------------------------------
    ! The simplest regime, and the one where a closed form exists outright:
    ! p(z)**rdocp = ps**rdocp - g*pref0**rdocp*z/(cp*thv).
    thv0    = 290.
    thlprof = thv0
    qtprof  = 0.
    call initbasestate(thlprof, qtprof)

    ! thv_b(kb) is what the retired thls sentinel used to supply. #302's bug gave
    ! thvs = thls*(1+(rv/rd-1)*qts) with qts unset (-1), i.e. ~0.39*thls; anything
    ! but thl(kb) here means that class of error is back.
    if (.not. close_rel(thv_b(kb), thv0, tol_exact)) then
      if (myid == 0) write(*, '(A,F12.5,A,F12.5)') '  FAIL: thv_b(kb) =', thv_b(kb), ' expected', thv0
      ok = .false.
    end if
    if (ph_b(kb) /= ps) then
      if (myid == 0) write(*, '(A,F12.3,A,F12.3)') '  FAIL: ph_b(kb) =', ph_b(kb), ' expected ps =', ps
      ok = .false.
    end if

    do k = kb, ke + kh
      if (.not. close_rel(pf_b(k), (ps**rdocp - grav*(pref0**rdocp)*zf(k)/(cp*thv0))**(1./rdocp), tol_exact)) then
        if (myid == 0) write(*, '(A,I4,A,F14.4)') '  FAIL: uniform pf_b(', k, ') =', pf_b(k)
        ok = .false.
      end if
      if (.not. close_rel(ph_b(k), (ps**rdocp - grav*(pref0**rdocp)*zh(k)/(cp*thv0))**(1./rdocp), tol_exact)) then
        if (myid == 0) write(*, '(A,I4,A,F14.4)') '  FAIL: uniform ph_b(', k, ') =', ph_b(k)
        ok = .false.
      end if
      if (.not. close_rel(exnf_b(k), (pf_b(k)/pref0)**rdocp, tol_exact)) then
        if (myid == 0) write(*, '(A,I4)') '  FAIL: exnf_b inconsistent with pf_b at k =', k
        ok = .false.
      end if
      if (.not. close_rel(exnh_b(k), (ph_b(k)/pref0)**rdocp, tol_exact)) then
        if (myid == 0) write(*, '(A,I4)') '  FAIL: exnh_b inconsistent with ph_b at k =', k
        ok = .false.
      end if
    end do
    do k = kb + 1, ke + kh
      if (pf_b(k) >= pf_b(k-1)) then
        if (myid == 0) write(*, '(A,I4)') '  FAIL: pf_b not decreasing at k =', k
        ok = .false.
      end if
    end do
    call exitbasestate

    ! ---- B. dry, linearly stratified ---------------------------------------
    ! thv(z) = thv0 + gamma*z. A uniform column cannot tell a correct thvh
    ! weighting from a wrong one -- every average of equal values is the same
    ! value -- so the interface weighting is only really exercised once thv
    ! varies with height. Two independent expectations are checked:
    !
    !   (i)  the telescoped sum, exact for any profile;
    !   (ii) the continuous solution, which for linear thv integrates to a log:
    !        p(z)**rdocp = ps**rdocp - (g*pref0**rdocp/(cp*gamma))*ln(1+gamma*z/thv0)
    !
    ! (i) pins the discretisation; (ii) says the discretisation is actually
    ! solving hydrostatic balance, which (i) alone would not catch if the scheme
    ! were self-consistently wrong.
    gamma = 0.01                                  ! K/m, a stable lapse
    do k = kb, ke
      thlprof(k) = thv0 + gamma*zf(k)
    end do
    qtprof = 0.
    call initbasestate(thlprof, qtprof)

    do k = kb, ke
      thv_expect(k) = thv0 + gamma*zf(k)
    end do
    thv_expect(ke+kh) = thv_expect(ke)
    do k = kb, ke + kh
      if (.not. close_rel(thv_b(k), thv_expect(k), tol_exact)) then
        if (myid == 0) write(*, '(A,I4,A,F12.5,A,F12.5)') '  FAIL: stratified thv_b(', k, ') =', thv_b(k), ' expected', thv_expect(k)
        ok = .false.
      end if
    end do

    if (.not. check_telescoped_sum(thv_expect, rdocp, tol_exact)) ok = .false.

    ! (ii) continuous solution. Only meaningful where thv really is linear, i.e.
    ! below ke; thv_b is flattened at the top ghost level by construction.
    do k = kb, ke
      if (.not. close_rel(pf_b(k), &
            (ps**rdocp - (grav*(pref0**rdocp)/(cp*gamma))*log(1. + gamma*zf(k)/thv0))**(1./rdocp), &
            tol_trunc)) then
        if (myid == 0) write(*, '(A,I4,A,F14.4)') '  FAIL: stratified pf_b vs continuous solution at k =', k, ' got', pf_b(k)
        ok = .false.
      end if
    end do
    call exitbasestate

    ! ---- C. moist, both thl and qt varying with height ----------------------
    ! thv = thl*(1 + (rv/rd - 1)*qt) with ql = 0. Varying qt as well as thl means
    ! a qt that is dropped, or mis-scaled, or read at the wrong level shows up
    ! here, where a constant qt would hide it in a single offset.
    qt0   = 8.e-3
    dqtdz = -5.e-5                                ! kg/kg per m, drying upward
    do k = kb, ke
      thlprof(k) = thv0 + gamma*zf(k)
      qtprof(k)  = qt0 + dqtdz*zf(k)
    end do
    call initbasestate(thlprof, qtprof)

    do k = kb, ke
      thv_expect(k) = thlprof(k)*(1. + (rv/rd - 1.)*qtprof(k))
    end do
    thv_expect(ke+kh) = thlprof(ke)*(1. + (rv/rd - 1.)*qtprof(ke))
    do k = kb, ke + kh
      if (.not. close_rel(thv_b(k), thv_expect(k), tol_exact)) then
        if (myid == 0) write(*, '(A,I4,A,F12.5,A,F12.5)') '  FAIL: moist thv_b(', k, ') =', thv_b(k), ' expected', thv_expect(k)
        ok = .false.
      end if
    end do

    if (.not. check_telescoped_sum(thv_expect, rdocp, tol_exact)) ok = .false.
    call exitbasestate

    deallocate(thlprof, qtprof, thv_expect)

    if (myid == 0) then
      if (ok) then
        write(*, '(A)') 'tests_basestate: PASS'
      else
        write(*, '(A)') 'tests_basestate: FAIL'
      end if
    end if

    tests_basestate = ok

  contains

    logical function close_rel(got, want, rtol)
      real, intent(in) :: got, want, rtol
      close_rel = abs(got - want) <= rtol*max(abs(want), tiny(want))
    end function close_rel

    !> pf_b/ph_b written as a running sum instead of the module's step-by-step
    !! recurrence. Exact for any profile, because each step of the recurrence
    !! subtracts an independent increment in p**rdocp space:
    !!
    !!   pf_b(k)**rdocp = ps**rdocp
    !!                  - (g*pref0**rdocp/cp)*[ zf(kb)/thv(kb)
    !!                                        + sum_{j=kb+1..k} dzh(j)/thvh(j) ]
    !!
    !! with thvh(j) the module's own interface weighting. This catches a wrong
    !! increment, a wrong dz, a sign slip or a mis-weighted interface, none of
    !! which a uniform-thv column can distinguish.
    logical function check_telescoped_sum(thv, rdocp_in, rtol)
      real, intent(in) :: thv(kb:ke+kh), rdocp_in, rtol
      real    :: accf, acch, thvh, fac, expect
      integer :: j
      check_telescoped_sum = .true.
      fac  = grav*(pref0**rdocp_in)/cp
      accf = zf(kb)/thv(kb)
      acch = 0.
      do j = kb + 1, ke + kh
        thvh = (thv(j)*dzf(j-1) + thv(j-1)*dzf(j))/(2.*dzh(j))
        accf = accf + dzh(j)/thvh
        acch = acch + dzf(j-1)/thv(j-1)
        expect = (ps**rdocp_in - fac*accf)**(1./rdocp_in)
        if (.not. close_rel(pf_b(j), expect, rtol)) then
          if (myid == 0) write(*, '(A,I4,A,F14.4,A,F14.4)') '  FAIL: telescoped pf_b(', j, ') =', pf_b(j), ' expected', expect
          check_telescoped_sum = .false.
        end if
        expect = (ps**rdocp_in - fac*acch)**(1./rdocp_in)
        if (.not. close_rel(ph_b(j), expect, rtol)) then
          if (myid == 0) write(*, '(A,I4,A,F14.4,A,F14.4)') '  FAIL: telescoped ph_b(', j, ') =', ph_b(j), ' expected', expect
          check_telescoped_sum = .false.
        end if
      end do
    end function check_telescoped_sum

  end function tests_basestate

end module tests
