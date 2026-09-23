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
  public :: tests_read_sparse_ijk, tests_2decomp_init_exit, tests_mpi_operators, tests_sgs_statistics

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
