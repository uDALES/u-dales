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
  public :: tests_read_sparse_ijk, tests_2decomp_init_exit, tests_mpi_operators, tests_zgrid_metrics, tests_sgs_statistics, &
            tests_sgs_oneeqn, tests_lstend, tests_thermo_levels

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

      check_loc_xy = compare_real_1d('avexy_ibm '//trim(label), got, exp, 1.e-9)

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

      check_loc_y = compare_real_2d('avey_ibm '//trim(label), got_avg, exp_avg, 1.e-9)
      if (.not. compare_real_2d('sumy_ibm '//trim(label), got_sum_y, exp_sum_y, 1.e-9)) check_loc_y = .false.
      if (.not. compare_real_2d_jk('sumx_ibm '//trim(label), got_sum_x, exp_sum_x, 1.e-9)) check_loc_y = .false.

      deallocate(var_clean, got_avg, got_sum_y, got_sum_x, exp_avg, exp_sum_y, exp_sum_x, &
                 sum_local, sum_global, sumx_local, sumx_global)
    end function check_loc_y

  end function tests_mpi_operators

  !> Verify the vertical grid metrics after initglobal (TEST_ZGRID_METRICS).
  !!
  !! The solver reads only the cell centres zf from prof.inp and derives every
  !! other vertical metric (zh, dzf, dzh, their inverses and the kappa-scheme
  !! copies). The rest of the code relies on identities of that construction,
  !! most importantly that each centre is the midpoint of its own faces (which
  !! makes the plain 0.5 half-to-full-level averages exact on stretched grids)
  !! and that dzh is the mean of the adjacent dzf. This test asserts those
  !! identities for the grid actually loaded, so it validates uniform and
  !! stretched prof.inp files alike. On a stretched grid it also requires a
  !! Poisson configuration that solves the vertical direction by Gaussian
  !! elimination (ipoiss=0 with BCzp=1, or ipoiss=3); the other branches
  !! assume equidistant spacing in their vertical eigenvalues and are already
  !! rejected by initglobal before this test runs.
  logical function tests_zgrid_metrics()
    use modglobal, only : kb, ke, kh, zf, zh, dzf, dzh, dzfi, dzhi, dzfiq, &
                          dzfi5, dzhiq, dzf2, dzh2i, dzfc, dzfci, dzhci, &
                          delta, dxf, dy, ipoiss, POISS_FFT2D, POISS_FFT2D_2DECOMP, BCzp
    implicit none

    real :: tol, stretch
    integer :: k, i, nfail
    logical :: stretched

    nfail = 0
    tol = 1000.*epsilon(zh(kb))*max(1., zh(ke+1))

    call expect(abs(zh(kb)) <= tol, 'zh(kb) is 0')

    do k = kb, ke
      call expect(dzf(k) > 0., 'dzf positive')
      call expect(abs(dzf(k) - (zh(k+1) - zh(k))) <= tol, 'dzf = zh(k+1)-zh(k)')
      ! the midpoint identity: 0.5*(w(k)+w(k+1)) style averages depend on it
      call expect(abs(zf(k) - 0.5*(zh(k) + zh(k+1))) <= tol, 'zf is the midpoint of its faces')
      if (k > kb) call expect(zf(k) > zf(k-1), 'zf strictly increasing')
    end do

    call expect(abs(dzh(kb) - dzf(kb)) <= tol, 'dzh(kb) ghost convention')
    do k = kb + 1, ke
      call expect(abs(dzh(k) - (zf(k) - zf(k-1))) <= tol, 'dzh = zf(k)-zf(k-1)')
      call expect(abs(dzh(k) - 0.5*(dzf(k-1) + dzf(k))) <= tol, 'dzh is the mean of adjacent dzf')
    end do
    call expect(abs(dzf(kb-1) - dzf(kb)) <= tol, 'dzf bottom ghost convention')
    call expect(abs(dzf(ke+1) - dzf(ke)) <= tol, 'dzf top ghost convention')

    do k = kb, ke + kh
      call expect(abs(dzfi(k)*dzf(k) - 1.) <= tol, 'dzfi = 1/dzf')
      call expect(abs(dzhi(k)*dzh(k) - 1.) <= tol, 'dzhi = 1/dzh')
      call expect(abs(dzfiq(k) - 0.25*dzfi(k)) <= tol, 'dzfiq = dzfi/4')
      call expect(abs(dzhiq(k) - 0.25*dzhi(k)) <= tol, 'dzhiq = dzhi/4')
      call expect(abs(dzfi5(k) - 0.5*dzfi(k)) <= tol, 'dzfi5 = dzfi/2')
      call expect(abs(dzf2(k) - dzf(k)*dzf(k)) <= tol*dzf(k), 'dzf2 = dzf**2')
      call expect(abs(dzh2i(k) - dzhi(k)*dzhi(k)) <= tol*dzhi(k)*dzhi(k), 'dzh2i = dzhi**2')
    end do

    ! kappa-scheme copies exist only when a kappa advection scheme is selected
    if (allocated(dzfc)) then
      do k = kb, ke
        call expect(abs(dzfc(k) - dzf(k)) <= tol, 'dzfc matches dzf on the interior')
        call expect(abs(dzfci(k)*dzfc(k) - 1.) <= tol, 'dzfci = 1/dzfc')
        call expect(abs(dzhci(k) - dzhi(k)) <= tol*dzhi(k), 'dzhci matches dzhi on the interior')
      end do
    end if

    do k = kb, ke
      i = lbound(delta, 1)
      call expect(abs(delta(i,k) - (dxf(i)*dy*dzf(k))**(1./3.)) <= tol, &
                  'delta = (dxf*dy*dzf)**(1/3)')
    end do

    ! a stretched grid must only be combined with a Poisson solve that treats the
    ! vertical direction by Gaussian elimination; the FFT branches hard-code
    ! dzfi(1) in their vertical eigenvalues (modpois: "Assumes equidistant")
    stretch = maxval(dzf(kb:ke))/minval(dzf(kb:ke))
    stretched = stretch > 1. + 100.*epsilon(stretch)
    if (stretched) then
      call expect((ipoiss == POISS_FFT2D .and. BCzp == 1) .or. ipoiss == POISS_FFT2D_2DECOMP, &
                  'stretched grid requires ipoiss=0 with BCzp=1 or ipoiss=3 (vertical Gaussian elimination)')
    end if

    if (myid == 0) then
      write(*,'(A,ES10.2,A,L1)') 'zgrid metrics: max/min dzf = ', stretch, ', stretched = ', stretched
      if (nfail == 0) then
        write(*,*) 'PASS: all vertical grid metric checks passed'
      else
        write(*,'(A,I0,A)') 'FAIL: ', nfail, ' vertical grid metric checks failed'
      end if
    end if

    tests_zgrid_metrics = (nfail == 0)

  contains

    subroutine expect(condition, label)
      logical, intent(in) :: condition
      character(*), intent(in) :: label
      if (.not. condition) then
        nfail = nfail + 1
        if (myid == 0) write(*,'(A,A,A,I0)') 'FAIL zgrid metric: ', label, ' at k=', k
      end if
    end subroutine expect

  end function tests_zgrid_metrics

  !> Validate the SGS fluxes written by the statistics output (modstatsdump).
  !! Uses the stretched-grid case 300 (libm = .false.) so that the vertical
  !! metric factors are genuinely non-uniform and dzf/dzh confusion is visible.
  !!
  !! All prescribed fields are pure functions of the GLOBAL coordinates, so the
  !! halos can be filled analytically and the test needs no halo exchange and is
  !! decomposition-invariant.
  !!
  !! Three independent assertions:
  !!   (a) manufactured solution: u = a*z^2, w = b*x, thl = c*z^2 with constant
  !!       ekm/ekh must reproduce -nu*(2*a*z + b) and -kappa*2*c*z, both exactly
  !!       in the discrete sense and to within the expected O(dz) truncation
  !!       error of the analytic profile;
  !!   (b) consistency with the solver: the vertical divergence of the SGS
  !!       fluxes must reproduce the tendency diffu/diffw/diffc add to
  !!       up/wp/thlp, to machine precision. This is what pins wsgs to a stress
  !!       (m^2/s^2) rather than a tendency (m/s^2);
  !!   (c) sign regression: for a monotone shear the SGS momentum flux must be
  !!       negative (down-gradient), guarding the sign convention.
  logical function tests_sgs_statistics()
    use modglobal,    only : ib, ie, ih, jb, je, jh, kb, ke, kh, runmode, &
                             dx, dzf, dzfi, dzh, dzhi, zf, zh
    use modfields,    only : initfields, um, vm, wm, thlm, u0, v0, w0, thl0, up, wp, thlp
    use modsubgrid,   only : initsubgrid, diffc, diffu, diffw, ekh, ekm
    use modboundary,  only : initboundary
    use modibm,       only : initibm, createmasks
    use modstatsdump, only : compute_sgs_fluxes
    use initfac,      only : readfacetfiles
    use decomp_2d,    only : zstart

    implicit none

    real, parameter :: acoef = 0.37   ! u   = acoef*z^2
    real, parameter :: bcoef = 0.21   ! w   = bcoef*x + dcoef*z
    real, parameter :: ccoef = 0.53   ! thl = ccoef*z^2
    real, parameter :: dcoef = 0.43
    real, parameter :: nuconst = 0.17 ! constant ekm for the manufactured solution
    real, parameter :: kaconst = 0.29 ! constant ekh for the manufactured solution

    real, allocatable :: usgs(:,:,:), vsgs(:,:,:), wsgs(:,:,:)
    real, allocatable :: thlsgs(:,:,:), qtsgs(:,:,:)
    real, allocatable :: sv1sgs(:,:,:), sv2sgs(:,:,:), sv3sgs(:,:,:), sv4sgs(:,:,:)
    real, allocatable :: zc(:), zw(:), got(:,:,:), exp(:,:,:)
    logical :: all_passed
    real    :: trunc_tol, scale
    integer :: i, j, k, gi

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_sgs_statistics: SGS STATISTICS TEST'
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'Using stretched-grid case 300 to validate modstatsdump SGS fluxes'
    end if

    call initfields
    call initboundary
    call initsubgrid
    call readfacetfiles
    call initibm
    call createmasks

    allocate(usgs  (ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(vsgs  (ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(wsgs  (ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(thlsgs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(qtsgs (ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(sv1sgs(ib:ie,jb:je,kb:ke+kh))
    allocate(sv2sgs(ib:ie,jb:je,kb:ke+kh))
    allocate(sv3sgs(ib:ie,jb:je,kb:ke+kh))
    allocate(sv4sgs(ib:ie,jb:je,kb:ke+kh))
    allocate(got(ib:ie,jb:je,kb:ke+kh), exp(ib:ie,jb:je,kb:ke+kh))

    ! Cell-centre heights, extended one level below the surface by mirroring
    ! (zh(kb) = 0), which is the level the kb stencils actually reach into.
    allocate(zc(kb-kh:ke+kh), zw(kb-kh:ke+kh))
    do k = kb, ke+kh
      zc(k) = zf(k)
      zw(k) = zh(k)
    end do
    zc(kb-kh) = zf(kb) - dzh(kb)
    zw(kb-kh) = zh(kb) - dzf(kb)

    all_passed = .true.

    !----------------------------------------------------------------------
    ! (a) Manufactured solution on the stretched grid
    !----------------------------------------------------------------------
    ekm = nuconst
    ekh = kaconst

    do k = kb-kh, ke+kh
      do j = jb-jh, je+jh
        do i = ib-ih, ie+ih
          gi = zstart(1) + i - 1
          um(i,j,k)   = acoef * zc(k)**2
          vm(i,j,k)   = 0.
          wm(i,j,k)   = bcoef * (real(gi) - 0.5) * dx + dcoef * zw(k)
          thlm(i,j,k) = ccoef * zc(k)**2
        end do
      end do
    end do

    call compute_sgs_fluxes(usgs,vsgs,wsgs,thlsgs,qtsgs,sv1sgs,sv2sgs,sv3sgs,sv4sgs)

    ! (a1) exact discrete identity: the second-order stencils are exact for a
    !      quadratic in z and a linear in x, so no truncation error enters here.
    do k = kb, ke+kh
      got(:,:,k) = usgs(ib:ie,jb:je,k)
      exp(:,:,k) = -nuconst * (acoef*(zc(k) + zc(k-1)) + bcoef)
    end do
    if (.not. compare_real_3d('usgs discrete', got, exp, 1.e-10)) all_passed = .false.

    do k = kb, ke+kh
      got(:,:,k) = thlsgs(ib:ie,jb:je,k)
      exp(:,:,k) = -kaconst * ccoef * (zc(k) + zc(k-1))
    end do
    if (.not. compare_real_3d('thlsgs discrete', got, exp, 1.e-10)) all_passed = .false.

    ! wsgs is the normal stress tau_33 = -2*K_m*dw/dz at the CELL CENTRE, so with
    ! dw/dz = dcoef it is exactly -2*nu*dcoef everywhere. If wsgs were the
    ! vertical part of the diffw tendency instead (units m/s^2, not m^2/s^2) this
    ! would be identically zero for a linear w and the check would fail.
    do k = kb, ke
      got(:,:,k) = wsgs(ib:ie,jb:je,k)
      exp(:,:,k) = -2. * nuconst * dcoef
    end do
    got(:,:,ke+kh) = 0.; exp(:,:,ke+kh) = 0.
    if (.not. compare_real_3d('wsgs discrete', got, exp, 1.e-10)) all_passed = .false.

    ! (a2) second-order consistency with the analytic flux at the w-levels.
    !      The stretched grid gives zf(k)+zf(k-1) - 2*zh(k) = (dzf(k)-dzf(k-1))/2,
    !      so anything worse than that bound (a dzf/dzh swap, a missing factor,
    !      a sign error) is caught here.
    trunc_tol = 0.
    do k = kb+1, ke+kh
      trunc_tol = max(trunc_tol, abs(dzf(k) - dzf(k-1)))
    end do
    trunc_tol = 0.5 * trunc_tol

    do k = kb, ke+kh
      got(:,:,k) = usgs(ib:ie,jb:je,k)
      exp(:,:,k) = -nuconst * (2.*acoef*zh(k) + bcoef)
    end do
    if (.not. compare_real_3d('usgs vs analytic', got, exp, &
                              1.05*nuconst*acoef*trunc_tol + 1.e-10)) all_passed = .false.

    do k = kb, ke+kh
      got(:,:,k) = thlsgs(ib:ie,jb:je,k)
      exp(:,:,k) = -kaconst * 2.*ccoef * zh(k)
    end do
    if (.not. compare_real_3d('thlsgs vs analytic', got, exp, &
                              1.05*kaconst*ccoef*trunc_tol + 1.e-10)) all_passed = .false.

    !----------------------------------------------------------------------
    ! (c) Sign regression: monotone positive shear must give a negative
    !     (down-gradient) SGS momentum flux. One-line guard on the convention.
    !----------------------------------------------------------------------
    if (maxval(usgs(ib:ie,jb:je,kb+1:ke)) >= 0.) then
      if (myid == 0) write(*,'(A,1X,A,1X,ES12.4)') 'FAIL', 'usgs sign', &
           maxval(usgs(ib:ie,jb:je,kb+1:ke))
      all_passed = .false.
    end if
    if (maxval(thlsgs(ib:ie,jb:je,kb+1:ke)) >= 0.) then
      if (myid == 0) write(*,'(A,1X,A,1X,ES12.4)') 'FAIL', 'thlsgs sign', &
           maxval(thlsgs(ib:ie,jb:je,kb+1:ke))
      all_passed = .false.
    end if

    !----------------------------------------------------------------------
    ! (b) Consistency with the solver's diffusion terms.
    !     u, v and thl vary in z only and w varies in x only, so the x- and
    !     y-terms of diffu/diffc vanish identically and what is left is exactly
    !     minus the vertical divergence of the SGS fluxes.
    !----------------------------------------------------------------------
    do k = kb-kh, ke+kh
      do j = jb-jh, je+jh
        do i = ib-ih, ie+ih
          gi = zstart(1) + i - 1
          ekm(i,j,k)  = 0.5 + 0.4*sin(1.7*real(gi) + 0.9*real(j) + 0.31*real(k))
          ekh(i,j,k)  = 0.6 + 0.3*cos(0.7*real(gi) - 1.1*real(j) + 0.53*real(k))
          um(i,j,k)   = sin(0.41*real(k)) + 0.3*cos(0.17*real(k))
          vm(i,j,k)   = cos(0.23*real(k)) - 0.2*sin(0.61*real(k))
          wm(i,j,k)   = 0.7*sin(0.29*(real(gi) - 0.5)*dx)
          thlm(i,j,k) = 288. + sin(0.37*real(k)) + 0.4*cos(0.13*real(k))
        end do
      end do
    end do
    u0   = um
    v0   = vm
    w0   = wm
    thl0 = thlm

    call compute_sgs_fluxes(usgs,vsgs,wsgs,thlsgs,qtsgs,sv1sgs,sv2sgs,sv3sgs,sv4sgs)

    up   = 0.
    thlp = 0.
    call diffu(up)
    call diffc(ih,jh,kh,thl0,thlp)

    do k = kb, ke
      got(:,:,k) = up(ib:ie,jb:je,k)
      exp(:,:,k) = -(usgs(ib:ie,jb:je,k+1) - usgs(ib:ie,jb:je,k)) * dzfi(k)
    end do
    got(:,:,ke+kh) = 0.; exp(:,:,ke+kh) = 0.
    scale = max(maxval(abs(exp)), 1.)
    if (.not. compare_real_3d('diffu vs div(usgs)', got, exp, 1.e-10*scale)) all_passed = .false.

    do k = kb, ke
      got(:,:,k) = thlp(ib:ie,jb:je,k)
      exp(:,:,k) = -(thlsgs(ib:ie,jb:je,k+1) - thlsgs(ib:ie,jb:je,k)) * dzfi(k)
    end do
    got(:,:,ke+kh) = 0.; exp(:,:,ke+kh) = 0.
    scale = max(maxval(abs(exp)), 1.)
    if (.not. compare_real_3d('diffc vs div(thlsgs)', got, exp, 1.e-10*scale)) all_passed = .false.

    ! Same idea for the normal stress: with u = v = 0 and w a function of z only,
    ! the horizontal terms of diffw vanish and what is left is exactly minus the
    ! vertical divergence of wsgs. This is the check that pins wsgs to a flux
    ! (m^2/s^2) rather than a tendency (m/s^2).
    do k = kb-kh, ke+kh
      do j = jb-jh, je+jh
        do i = ib-ih, ie+ih
          um(i,j,k) = 0.
          vm(i,j,k) = 0.
          wm(i,j,k) = sin(0.33*real(k)) + 0.5*cos(0.19*real(k))
        end do
      end do
    end do
    u0 = um
    v0 = vm
    w0 = wm

    call compute_sgs_fluxes(usgs,vsgs,wsgs,thlsgs,qtsgs,sv1sgs,sv2sgs,sv3sgs,sv4sgs)

    wp = 0.
    call diffw(wp)

    got = 0.; exp = 0.
    do k = kb+1, ke
      got(:,:,k) = wp(ib:ie,jb:je,k)
      exp(:,:,k) = -(wsgs(ib:ie,jb:je,k) - wsgs(ib:ie,jb:je,k-1)) * dzhi(k)
    end do
    scale = max(maxval(abs(exp)), 1.)
    if (.not. compare_real_3d('diffw vs div(wsgs)', got, exp, 1.e-10*scale)) all_passed = .false.

    deallocate(usgs, vsgs, wsgs, thlsgs, qtsgs, sv1sgs, sv2sgs, sv3sgs, sv4sgs, got, exp, zc, zw)

    if (all_passed .and. myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'ALL TESTS PASSED: tests_sgs_statistics'
      write(*, '(A)') '================================================'
    else if ((.not. all_passed) .and. myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'TESTS FAILED: tests_sgs_statistics'
      write(*, '(A)') '================================================'
    end if

    tests_sgs_statistics = all_passed

  end function tests_sgs_statistics

  !> One-equation SGS TKE model (Task 10, aligned with DALES dev 075ac4e modsubgrid.f90).
  !! Serial, on a synthetic stratified state; the closure and sources are called directly.
  !!   (a) closure: Kh/Km (molecular parts removed) is 1/Prandtl = 3 where dthvdz <= 0 and
  !!       ch1 + ch2*zlt/delta where dthvdz > 0, with zlt from grav/thvf(k) (not thvs);
  !!   (b) sources without the MO surface: no TKE source at kb, set explicitly (the
  !!       work arrays are pre-filled with garbage), and interior dissipation without
  !!       the old factor 2, buoyancy with grav/thvf(k);
  !!   (c) sources with the MO surface: the kb sources use dudz_mo/dvdz_mo and dthvdz(kb);
  !!   (d) a synthetic fully solid level k0 > kb (IIcs(k0) = 0, thvf(k0) = -999, the
  !!       avexy_ibm fill) with dthvdz > 0: neutral length scale, no buoyancy source,
  !!       everything finite (before the guard: sqrt of a negative number).
  !! The namelist must set loneeqn = .true. and ldelta = .false.; Prandtl must be 1/3 to
  !! full precision (the default).
  logical function tests_sgs_oneeqn()
    use modglobal,      only : ib, ie, jb, je, kb, ke, runmode, delta, grav, eps1, numol, &
                               prandtlmoli, pi, BCbotm, BCbotm_mo
    use modfields,      only : initfields, u0, v0, w0, e120, e12p, dthvdz, thvf, IIcs
    use modsurfdata,    only : thvs
    use modsubgrid,     only : initsubgrid, closure, sources
    use modsubgriddata, only : ekm, ekh, zlt, sbshr, sbbuo, sbdiss, loneeqn, ldelta, Prandtl, &
                               cn, cf, Rigc, alpha_kolm
    use modsurface,     only : initsurface_mo, lmo_surface, dudz_mo, dvdz_mo
    use modboundary,    only : initboundary
    use modibm,         only : initibm, createmasks
    use initfac,        only : readfacetfiles
    use decomp_2d,      only : zstart
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan

    implicit none

    real, parameter :: tol = 1.e-10
    real, parameter :: dudz0 = 0.05, dvdz0 = -0.03, dthvdz0 = -0.01
    real :: kmt, kht, ratio, zexp, cm_l, ceps_l, ce1_l, ce2_l, expd, expb, exps, worst
    integer :: i, j, k, gi, nstable, nunstable, nlimited, BCbotm_save, k0, IIcs_save
    real :: thvf_save
    logical :: all_passed

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_sgs_oneeqn: ONE-EQUATION SGS MODEL TEST'
      write(*, '(A)') '------------------------------------------------'
    end if

    call initfields
    call initboundary
    call initsubgrid
    call initsurface_mo
    call readfacetfiles
    call initibm
    call createmasks

    all_passed = .true.
    if ((.not. loneeqn) .or. ldelta .or. abs(1./Prandtl - 3.) > 1.e-12) then
      if (myid == 0) write(*,'(A)') 'FAIL namelist must set loneeqn, .not. ldelta and Prandtl = 1/3'
      tests_sgs_oneeqn = .false.
      return
    end if

    ! Synthetic stratified state. thvs is set far from thvf so a regression to grav/thvs fails.
    thvs = 250.
    do k = kb, ke
      thvf(k) = 300. + 0.01*real(k)
    end do
    do k = lbound(e120,3), ubound(e120,3)
      do j = lbound(e120,2), ubound(e120,2)
        do i = lbound(e120,1), ubound(e120,1)
          gi = zstart(1) + i - 1
          e120(i,j,k) = 0.3 + 0.1*sin(0.7*real(gi) + 1.1*real(j) + 0.3*real(k))
          u0(i,j,k)   = 2. + 0.5*sin(0.37*real(gi) + 0.21*real(k))
          v0(i,j,k)   = -1. + 0.4*cos(0.29*real(j) + 0.17*real(k))
          w0(i,j,k)   = 0.1*sin(0.53*real(gi) - 0.41*real(j) + 0.11*real(k))
        end do
      end do
    end do
    do k = kb, ke
      do j = jb, je
        do i = ib, ie
          gi = zstart(1) + i - 1
          dthvdz(i,j,k) = 0.02*sin(1.3*real(gi) + 0.7*real(j) + 0.4*real(k))
        end do
      end do
    end do

    !----------------------------------------------------------------------
    ! (a) closure: Kh/Km and the stable length scale
    !----------------------------------------------------------------------
    call closure
    nstable = 0; nunstable = 0; nlimited = 0; worst = 0.
    do k = kb, ke
      do j = jb, je
        do i = ib, ie
          kmt = ekm(i,j,k) - numol
          kht = ekh(i,j,k) - numol*prandtlmoli
          ratio = kht / kmt
          if (dthvdz(i,j,k) <= 0.) then
            nunstable = nunstable + 1
            worst = max(worst, abs(zlt(i,j,k) - delta(i,k))/delta(i,k), abs(ratio - 3.))
          else
            nstable = nstable + 1
            zexp = min(delta(i,k), cn*e120(i,j,k)/sqrt(grav/thvf(k)*max(abs(dthvdz(i,j,k)),eps1)))
            if (zexp < delta(i,k)) nlimited = nlimited + 1
            worst = max(worst, abs(zlt(i,j,k) - zexp)/zexp, &
                        abs(ratio - (1. + (3. - 1.)*zlt(i,j,k)/delta(i,k))))
          end if
        end do
      end do
    end do
    if (myid == 0) write(*,'(A,3(1X,I0),1X,A,ES12.4)') 'closure cells unstable/stable/limited', &
         nunstable, nstable, nlimited, 'worst', worst
    if (worst > tol) then
      if (myid == 0) write(*,'(A,ES12.4)') 'FAIL closure Kh/Km or zlt, worst deviation', worst
      all_passed = .false.
    end if
    if (nunstable == 0 .or. nstable == 0 .or. nlimited == 0) then
      if (myid == 0) write(*,'(A)') 'FAIL closure test is vacuous (need unstable, stable and limited cells)'
      all_passed = .false.
    end if

    !----------------------------------------------------------------------
    ! (b) sources without the MO surface
    !----------------------------------------------------------------------
    if (lmo_surface()) then
      if (myid == 0) write(*,'(A)') 'FAIL namelist must not select the MO surface (BCbotm = 4/5)'
      all_passed = .false.
    end if
    sbshr(:,:,kb) = 1.e30; sbbuo(:,:,kb) = -1.e30; sbdiss(:,:,kb) = 1.e30
    e12p = 0.
    call sources
    worst = maxval(abs(e12p(ib:ie,jb:je,kb)))
    if (worst /= 0.) then
      if (myid == 0) write(*,'(A,ES12.4)') 'FAIL e12p(kb) without MO surface is not zero', worst
      all_passed = .false.
    end if
    if (maxval(abs(e12p(ib:ie,jb:je,kb+1))) == 0.) then
      if (myid == 0) write(*,'(A)') 'FAIL e12p(kb+1) is zero: sources test is vacuous'
      all_passed = .false.
    end if

    ! interior dissipation and buoyancy, constants derived independently of initsubgrid
    cm_l   = cf / (2. * pi) * (1.5*alpha_kolm)**(-1.5)
    ceps_l = 2. * pi / cf * (1.5*alpha_kolm)**(-1.5)
    ce1_l  = cn**2 * (cm_l/Rigc - cm_l)
    ce2_l  = ceps_l - ce1_l
    worst = 0.
    do k = kb+1, ke
      do j = jb, je
        do i = ib, ie
          expd = -(ce1_l + ce2_l*zlt(i,j,k)/delta(i,k)) * e120(i,j,k)**2 / (2.*zlt(i,j,k))
          expb = -(ekh(i,j,k) - numol*prandtlmoli)*grav/thvf(k)*dthvdz(i,j,k)/(2.*e120(i,j,k))
          worst = max(worst, abs(sbdiss(i,j,k) - expd)/abs(expd), &
                      abs(sbbuo(i,j,k) - expb)/max(abs(expb), 1.e-30))
        end do
      end do
    end do
    if (myid == 0) write(*,'(A,ES12.4)') 'interior sbdiss/sbbuo worst relative deviation', worst
    if (worst > tol) then
      if (myid == 0) write(*,'(A,ES12.4)') 'FAIL interior sbdiss or sbbuo', worst
      all_passed = .false.
    end if

    !----------------------------------------------------------------------
    ! (c) sources with the MO surface: uniform kb plane, so the only shear is the MO one
    !----------------------------------------------------------------------
    BCbotm_save = BCbotm
    BCbotm = BCbotm_mo
    u0(:,:,kb) = 3.; v0(:,:,kb) = -2.; w0(:,:,kb) = 0.; w0(:,:,kb+1) = 0.
    dudz_mo = dudz0; dvdz_mo = dvdz0
    dthvdz(:,:,kb) = dthvdz0
    call closure
    e12p = 0.
    call sources
    worst = 0.
    do j = jb, je
      do i = ib, ie
        exps = (ekm(i,j,kb) - numol)*(dudz0**2 + dvdz0**2)/(2.*e120(i,j,kb)) &
             - (ekh(i,j,kb) - numol*prandtlmoli)*grav/thvf(kb)*dthvdz0/(2.*e120(i,j,kb)) &
             - (ce1_l + ce2_l*zlt(i,j,kb)/delta(i,kb)) * e120(i,j,kb)**2 / (2.*zlt(i,j,kb))
        worst = max(worst, abs(e12p(i,j,kb) - exps)/abs(exps))
      end do
    end do
    BCbotm = BCbotm_save
    if (myid == 0) write(*,'(A,ES12.4)') 'MO kb sources worst relative deviation', worst
    if (worst > tol) then
      if (myid == 0) write(*,'(A,ES12.4)') 'FAIL MO kb sources', worst
      all_passed = .false.
    end if

    !----------------------------------------------------------------------
    ! (d) fully solid level: the grav/thvf guard
    !----------------------------------------------------------------------
    k0 = kb + 2
    IIcs_save = IIcs(k0); thvf_save = thvf(k0)
    IIcs(k0) = 0
    thvf(k0) = -999.
    dthvdz(:,:,k0) = 0.01
    call closure
    e12p = 0.
    call sources
    worst = 0.
    do j = jb, je
      do i = ib, ie
        worst = max(worst, abs(zlt(i,j,k0) - delta(i,k0))/delta(i,k0), abs(sbbuo(i,j,k0)))
        if (ieee_is_nan(zlt(i,j,k0)) .or. ieee_is_nan(ekm(i,j,k0)) .or. ieee_is_nan(ekh(i,j,k0)) &
            .or. ieee_is_nan(sbshr(i,j,k0)) .or. ieee_is_nan(sbdiss(i,j,k0)) &
            .or. ieee_is_nan(e12p(i,j,k0))) worst = huge(1.)
      end do
    end do
    IIcs(k0) = IIcs_save; thvf(k0) = thvf_save
    if (myid == 0) write(*,'(A,ES12.4)') 'fully solid level guard worst deviation', worst
    if (worst > tol) then
      if (myid == 0) write(*,'(A,ES12.4)') 'FAIL fully solid level guard', worst
      all_passed = .false.
    end if

    if (all_passed .and. myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'ALL TESTS PASSED: tests_sgs_oneeqn'
      write(*, '(A)') '================================================'
    else if ((.not. all_passed) .and. myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'TESTS FAILED: tests_sgs_oneeqn'
      write(*, '(A)') '================================================'
    end if

    tests_sgs_oneeqn = all_passed

  end function tests_sgs_oneeqn

  !> Large-scale subsidence (modforces:lstend) against DALES dev 075ac4e lstend, cell by
  !! cell. thl0, qt0, sv0(:,:,:,1), u0 and v0 vary in x and z, whls changes sign with
  !! height, lmomsubs is on and the horizontal large-scale gradients are zero, so the
  !! tendency is the upwind subsidence term only:
  !!   whls(k+1) < 0:     -whls(k+1)*(phi(k+1) - phi(k))/dzh(k+1)
  !!   else, k > kb:      -whls(k)  *(phi(k)   - phi(k-1))/dzh(k)
  !!   else, k = kb:      0 (mean ascent neglected at the first level)
  !! Two passes with opposite whls cover both kb rules. A slab-mean stencil fails the
  !! test, because the fields vary in x. The namelist must set nsv >= 1.
  logical function tests_lstend()
    use modglobal, only : ib, ie, jb, je, kb, ke, dzh, nsv, lmomsubs, itot, runmode, pi
    use modfields, only : initfields, thl0, qt0, sv0, u0, v0, thlp, qtp, svp, up, vp, whls, &
                          u0av, v0av
    use modforces, only : lstend
    use modmpi,    only : mpi_max
    use decomp_2d, only : zstart

    implicit none

    real, parameter :: tol = 1.e-12
    integer :: i, j, k, gi, ipass, n, ndown, nup, nkb
    real :: worst, scale, sgn, e, s, xvar, loc(2), glob(2)
    logical :: all_passed, lmomsubs_save
    real, allocatable :: expt(:,:,:,:), got(:,:,:,:), phi(:,:,:,:)

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_lstend: LOCAL-FIELD SUBSIDENCE TEST'
      write(*, '(A)') '------------------------------------------------'
    end if

    call initfields
    all_passed = .true.
    if (nsv < 1) then
      if (myid == 0) write(*,'(A)') 'FAIL namelist must set nsv >= 1'
      tests_lstend = .false.
      return
    end if
    lmomsubs_save = lmomsubs
    lmomsubs = .true.
    u0av = 0.; v0av = 0.

    do k = lbound(thl0,3), ubound(thl0,3)
      do j = lbound(thl0,2), ubound(thl0,2)
        do i = lbound(thl0,1), ubound(thl0,1)
          gi = zstart(1) + i - 1
          s = sin(2.*pi*real(gi)/real(itot))
          thl0(i,j,k) = 290. + 0.004*real(k) + 0.3*s*cos(0.2*real(k)) + 0.01*real(j)
          qt0(i,j,k)  = 1.e-2 - 2.e-5*real(k) + 1.e-4*s*sin(0.3*real(k))
          u0(i,j,k)   = 5. + 0.05*real(k) + 0.8*s*cos(0.15*real(k))
          v0(i,j,k)   = -2. + 0.03*real(k) - 0.6*s*sin(0.25*real(k) + 0.1*real(j))
        end do
      end do
    end do
    do k = lbound(sv0,3), ubound(sv0,3)
      do j = lbound(sv0,2), ubound(sv0,2)
        do i = lbound(sv0,1), ubound(sv0,1)
          gi = zstart(1) + i - 1
          s = sin(2.*pi*real(gi)/real(itot))
          sv0(i,j,k,1) = 1. + 0.02*real(k) + 0.5*s*cos(0.35*real(k))
        end do
      end do
    end do

    allocate(expt(ib:ie,jb:je,kb:ke,5), got(ib:ie,jb:je,kb:ke,5), phi(ib:ie,jb:je,kb-1:ke+1,5))
    phi(:,:,:,1) = thl0(ib:ie,jb:je,kb-1:ke+1)
    phi(:,:,:,2) = qt0 (ib:ie,jb:je,kb-1:ke+1)
    phi(:,:,:,3) = sv0 (ib:ie,jb:je,kb-1:ke+1,1)
    phi(:,:,:,4) = u0  (ib:ie,jb:je,kb-1:ke+1)
    phi(:,:,:,5) = v0  (ib:ie,jb:je,kb-1:ke+1)

    ndown = 0; nup = 0; nkb = 0; worst = 0.; xvar = 0.
    do ipass = 1, 2
      sgn = merge(1., -1., ipass == 1)
      ! descent aloft and ascent below (pass 1), the opposite in pass 2; zero at one level
      do k = lbound(whls,1), ubound(whls,1)
        whls(k) = sgn * 0.01 * real((kb + ke)/2 - k) / real(ke - kb + 1)
      end do
      thlp = 0.; qtp = 0.; svp = 0.; up = 0.; vp = 0.
      call lstend
      got(:,:,:,1) = thlp(ib:ie,jb:je,kb:ke)
      got(:,:,:,2) = qtp (ib:ie,jb:je,kb:ke)
      got(:,:,:,3) = svp (ib:ie,jb:je,kb:ke,1)
      got(:,:,:,4) = up  (ib:ie,jb:je,kb:ke)
      got(:,:,:,5) = vp  (ib:ie,jb:je,kb:ke)
      do k = kb, ke
        if (whls(k+1) < 0.) then
          ndown = ndown + 1
          expt(:,:,k,:) = -whls(k+1) * (phi(:,:,k+1,:) - phi(:,:,k,:)) / dzh(k+1)
        else if (k > kb) then
          nup = nup + 1
          expt(:,:,k,:) = -whls(k) * (phi(:,:,k,:) - phi(:,:,k-1,:)) / dzh(k)
        else
          nkb = nkb + 1
          expt(:,:,k,:) = 0.
        end if
      end do
      do n = 1, 5
        scale = max(maxval(abs(expt(:,:,:,n))), tiny(1.))
        e = maxval(abs(got(:,:,:,n) - expt(:,:,:,n))) / scale
        worst = max(worst, e)
        xvar = max(xvar, maxval(abs(got(ie,jb,:,n) - got(ib,jb,:,n))) / scale)
      end do
    end do
    loc = (/worst, xvar/)
    call MPI_ALLREDUCE(loc, glob, 2, MY_REAL, MPI_MAX, comm3d, mpierr)
    worst = glob(1); xvar = glob(2)
    lmomsubs = lmomsubs_save
    deallocate(expt, got, phi)

    if (myid == 0) write(*,'(A,3(1X,I0),1X,A,ES12.4,1X,A,ES12.4)') &
         'lstend levels descent/ascent/kb-ascent', ndown, nup, nkb, &
         'worst relative deviation', worst, 'x variation', xvar
    if (worst > tol) then
      if (myid == 0) write(*,'(A,ES12.4)') 'FAIL lstend differs from the DALES upwind formula', worst
      all_passed = .false.
    end if
    if (ndown == 0 .or. nup == 0 .or. nkb == 0 .or. xvar < 1.e-3) then
      if (myid == 0) write(*,'(A)') 'FAIL lstend test is vacuous (need descent, ascent, kb ascent, x variation)'
      all_passed = .false.
    end if

    if (all_passed .and. myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'ALL TESTS PASSED: tests_lstend'
      write(*, '(A)') '================================================'
    else if ((.not. all_passed) .and. myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'TESTS FAILED: tests_lstend'
      write(*, '(A)') '================================================'
    end if

    tests_lstend = all_passed

  end function tests_lstend

  !> Saturation adjustment at full and half levels (#403). thermodynamics is called on a
  !! saturated state that varies in x and z; then ql0(k) must be the one-step adjustment
  !! of thl0(k), qt0(k) at the pressure that thermo saw (presf/exnf before the call), and
  !! ql0h(k) that of thl0h(k), qt0h(k) at presh/exnh (thermo runs after diagfld and
  !! calc_halflev for the half levels). Before #403 ql0(k) held the ql of level k+1.
  !! The namelist must set lmoist = .true. and lqlnr = .false. (the default).
  logical function tests_thermo_levels()
    use modglobal,  only : ib, ie, jb, je, kb, ke, kh, runmode, lmoist, timee, itot, pi, &
                           es0, at, bt, rd, rv, rlv, cp, tmelt
    use modfields,  only : initfields, thl0, qt0, ql0, thl0h, qt0h, ql0h, presf, presh, exnf, exnh
    use modthermodynamics, only : initthermodynamics, thermodynamics, lqlnr
    use modmpi,     only : mpi_max, mpi_min
    use modboundary, only : initboundary
    use modibm,     only : initibm, createmasks
    use initfac,    only : readfacetfiles
    use decomp_2d,  only : zstart

    implicit none

    real, parameter :: tol = 1.e-12
    integer :: i, j, k, gi, nsat_f, nsat_h
    real :: s, worst_f, worst_h, qlmax, loc(4), glob(4)
    real, allocatable :: pf(:), ef(:)
    logical :: all_passed

    if (myid == 0) then
      write(*, '(A)') '================================================'
      write(*, '(A, I8)') 'runmode = ', runmode
      write(*, '(A)') 'tests_thermo_levels: SATURATION ADJUSTMENT LEVELS'
      write(*, '(A)') '------------------------------------------------'
    end if

    call initfields
    call initboundary
    call initthermodynamics
    call readfacetfiles
    call initibm
    call createmasks

    all_passed = .true.
    if ((.not. lmoist) .or. lqlnr) then
      if (myid == 0) write(*,'(A)') 'FAIL namelist must set lmoist = .true. and lqlnr = .false.'
      tests_thermo_levels = .false.
      return
    end if

    ! saturated mixed layer with a lapse in thl and a dip in qt; varies in x
    do k = lbound(thl0,3), ubound(thl0,3)
      do j = lbound(thl0,2), ubound(thl0,2)
        do i = lbound(thl0,1), ubound(thl0,1)
          gi = zstart(1) + i - 1
          s = sin(2.*pi*real(gi)/real(itot))
          thl0(i,j,k) = 287. + 0.004*real(k) + 0.05*s
          qt0(i,j,k)  = 1.2e-2 - 1.e-5*real(k) + 2.e-4*s
        end do
      end do
    end do

    timee = 0.
    call thermodynamics            ! sets presf/exnf (diagfld) for the next call
    timee = 1.                     ! no extra diagfld before thermo
    allocate(pf(kb:ke+kh), ef(kb:ke+kh))
    pf = presf; ef = exnf
    call thermodynamics

    worst_f = 0.; worst_h = 0.; nsat_f = 0; nsat_h = 0; qlmax = 0.
    do k = kb, ke
      do j = jb, je
        do i = ib, ie
          s = ql_ref(thl0(i,j,k), qt0(i,j,k), pf(k), ef(k))
          if (s > 0.) nsat_f = nsat_f + 1
          qlmax = max(qlmax, s)
          worst_f = max(worst_f, abs(ql0(i,j,k) - s))
          if (k > kb) then
            s = ql_ref(thl0h(i,j,k), qt0h(i,j,k), presh(k), exnh(k))
            if (s > 0.) nsat_h = nsat_h + 1
            worst_h = max(worst_h, abs(ql0h(i,j,k) - s))
          end if
        end do
      end do
    end do
    loc = (/worst_f, worst_h, qlmax, real(min(nsat_f, nsat_h))/)
    call MPI_ALLREDUCE(loc(1:3), glob(1:3), 3, MY_REAL, MPI_MAX, comm3d, mpierr)
    call MPI_ALLREDUCE(loc(4), glob(4), 1, MY_REAL, MPI_MIN, comm3d, mpierr)
    worst_f = glob(1) / max(glob(3), tiny(1.)); worst_h = glob(2) / max(glob(3), tiny(1.))
    deallocate(pf, ef)

    if (myid == 0) write(*,'(A,ES12.4,1X,A,ES12.4,1X,A,ES12.4,1X,A,1X,I0)') &
         'ql0 full-level worst', worst_f, 'ql0h half-level worst', worst_h, 'max ql', glob(3), &
         'min saturated cells', nint(glob(4))
    if (worst_f > tol) then
      if (myid == 0) write(*,'(A,ES12.4)') 'FAIL ql0 does not match thl0/qt0 at the same level', worst_f
      all_passed = .false.
    end if
    if (worst_h > tol) then
      if (myid == 0) write(*,'(A,ES12.4)') 'FAIL ql0h does not match thl0h/qt0h at the same level', worst_h
      all_passed = .false.
    end if
    if (glob(4) < 1. .or. glob(3) < 1.e-4) then
      if (myid == 0) write(*,'(A)') 'FAIL thermo levels test is vacuous (state not saturated)'
      all_passed = .false.
    end if

    if (all_passed .and. myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'ALL TESTS PASSED: tests_thermo_levels'
      write(*, '(A)') '================================================'
    else if ((.not. all_passed) .and. myid == 0) then
      write(*, '(A)') '------------------------------------------------'
      write(*, '(A)') 'TESTS FAILED: tests_thermo_levels'
      write(*, '(A)') '================================================'
    end if
    tests_thermo_levels = all_passed

  contains

    !> uDALES modthermodynamics:thermo, one linearised step (lqlnr = .false.)
    real function ql_ref(thl, qt, p, exn)
      real, intent(in) :: thl, qt, p, exn
      real :: tl, es, qsl, b1, qs
      tl  = max(thl*exn, 100.)
      es  = es0*exp(at*(tl-tmelt)/(tl-bt))
      qsl = rd/rv*es/(p-(1-rd/rv)*es)
      b1  = rlv**2/(tl**2*cp*rv)
      qs  = qsl*(1.+b1*qt)/(1.+b1*qsl)
      ql_ref = max(qt-qs, 0.)
    end function ql_ref

  end function tests_thermo_levels

  !> Shared comparators. Report the largest absolute deviation and where it is.
  logical function compare_real_1d(label, got, exp, tol)
    implicit none
    character(len=*), intent(in) :: label
    real, intent(in) :: got(:), exp(:)
    real, intent(in) :: tol
    real :: max_abs
    integer :: imax(1)

    max_abs = maxval(abs(got - exp))
    compare_real_1d = max_abs <= tol
    if ((.not. compare_real_1d) .and. myid == 0) then
      imax = maxloc(abs(got - exp))
      write(*,'(A,1X,A,1X,ES12.4,1X,A,I0,1X,A,ES12.4,1X,A,ES12.4)') &
           'FAIL', trim(label), max_abs, 'idx', imax(1), 'got', got(imax(1)), 'exp', exp(imax(1))
    end if
  end function compare_real_1d

  logical function compare_real_2d(label, got, exp, tol)
    implicit none
    character(len=*), intent(in) :: label
    real, intent(in) :: got(:,:), exp(:,:)
    real, intent(in) :: tol
    real :: max_abs
    integer :: imax(2)

    max_abs = maxval(abs(got - exp))
    compare_real_2d = max_abs <= tol
    if ((.not. compare_real_2d) .and. myid == 0) then
      imax = maxloc(abs(got - exp))
      write(*,'(A,1X,A,1X,ES12.4,1X,A,I0,A,I0,1X,A,ES12.4,1X,A,ES12.4)') &
           'FAIL', trim(label), max_abs, 'idx', imax(1), ',', imax(2), 'got', got(imax(1),imax(2)), 'exp', exp(imax(1),imax(2))
    end if
  end function compare_real_2d

  logical function compare_real_2d_jk(label, got, exp, tol)
    implicit none
    character(len=*), intent(in) :: label
    real, intent(in) :: got(:,:), exp(:,:)
    real, intent(in) :: tol
    real :: max_abs
    integer :: imax(2)

    max_abs = maxval(abs(got - exp))
    compare_real_2d_jk = max_abs <= tol
    if ((.not. compare_real_2d_jk) .and. myid == 0) then
      imax = maxloc(abs(got - exp))
      write(*,'(A,1X,A,1X,ES12.4,1X,A,I0,A,I0,1X,A,ES12.4,1X,A,ES12.4)') &
           'FAIL', trim(label), max_abs, 'idx', imax(1), ',', imax(2), 'got', got(imax(1),imax(2)), 'exp', exp(imax(1),imax(2))
    end if
  end function compare_real_2d_jk

  logical function compare_real_3d(label, got, exp, tol)
    implicit none
    character(len=*), intent(in) :: label
    real, intent(in) :: got(:,:,:), exp(:,:,:)
    real, intent(in) :: tol
    real :: max_abs
    integer :: imax(3)

    max_abs = maxval(abs(got - exp))
    compare_real_3d = max_abs <= tol
    if ((.not. compare_real_3d) .and. myid == 0) then
      imax = maxloc(abs(got - exp))
      write(*,'(A,1X,A,1X,ES12.4,1X,A,I0,A,I0,A,I0,1X,A,ES12.4,1X,A,ES12.4,1X,A,ES12.4)') &
           'FAIL', trim(label), max_abs, 'idx', imax(1), ',', imax(2), ',', imax(3), &
           'got', got(imax(1),imax(2),imax(3)), 'exp', exp(imax(1),imax(2),imax(3)), 'tol', tol
    end if
  end function compare_real_3d

end module tests
