!> \file tests_cuda.f90
!! Opt-in CUDA tests driven by the external GPU integration-test harness.

module tests_cuda
#if defined(_GPU) && defined(UDALES_DEBUG)
   use cudafor
   use modadvection, only : advecc_upw_cuda, rlim_cuda
   use modboundary,  only : bcpup_pup_BCxm_driver_cuda
   use modcuda,  only : blockdim, griddim, checkCUDA, initfield, &
                        thlptothlpc_cuda, thlpctothlp_cuda, &
                        dxfci_d, dzfci_d, u0_d, v0_d, w0_d, &
                        thlp_d, thlpc_d, pup_d, up_d, u0driver_d
   use modfields, only : u0, v0, w0
   use modglobal, only : ib, ie, jb, je, kb, ke, ih, jh, kh, &
                         ihc, jhc, khc, dxfci, dyi, dzfci, eps1
   use modinletdata, only : u0driver
   use modmpi,   only : myid

   implicit none
   private

   public :: run_cuda_selftests_if_requested

contains

   !> Run CUDA self-tests only when requested by the external test harness.
   subroutine run_cuda_selftests_if_requested
      implicit none

      character(len=32) :: request
      integer :: status

      request = ''
      call get_environment_variable('UDALES_RUN_CUDA_SELFTEST', value=request, status=status)
      if (status /= 0) return

      select case (trim(adjustl(request)))
      case ('1', 'true', 'TRUE', 'yes', 'YES', 'on', 'ON')
         call test_initfield_extended_halos
         call test_kappa_limiter
         call test_upwind_scalar_advection
         call test_temperature_tendency_copy
         call test_driver_inlet_boundary
         write(*,'(A,I0)') 'CUDA device self-tests passed. rank=', myid
      case ('', '0', 'false', 'FALSE', 'no', 'NO', 'off', 'OFF')
         return
      case default
         if (myid == 0) then
            write(*,'(A,A)') 'Invalid UDALES_RUN_CUDA_SELFTEST value: ', trim(request)
         end if
         error stop 1
      end select
   end subroutine run_cuda_selftests_if_requested

   !> Verify that initfield writes every cell in the extended scalar halos.
   subroutine test_initfield_extended_halos
      implicit none

      real, device, allocatable :: test_d(:, :, :)
      real, allocatable :: test_h(:, :, :)

      allocate(test_d(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb:ke+khc))
      allocate(test_h(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb:ke+khc))

      test_d = 1.
      call initfield<<<griddim,blockdim>>>(test_d, 0., ihc, jhc, khc)
      call checkCUDA(cudaGetLastError(), 'extended-halo initfield self-test launch')
      call checkCUDA(cudaDeviceSynchronize(), 'extended-halo initfield self-test synchronization')
      test_h = test_d

      if (any(test_h /= 0.)) then
         call fail_cuda_selftest('extended-halo initfield')
      end if

      deallocate(test_h, test_d)
   end subroutine test_initfield_extended_halos

   !> Check the device limiter against the same algebra evaluated on the host.
   subroutine test_kappa_limiter
      implicit none

      integer, parameter :: nvalues = 6
      real, parameter :: d1_h(nvalues) = [1., 1., -1., 2., 0.25, -0.5]
      real, parameter :: d2_h(nvalues) = [1., -1., -1., 0.5, 0.5, 0.125]
      real, device, allocatable :: d1_d(:), d2_d(:), result_d(:)
      real :: expected(nvalues), result_h(nvalues)
      real :: phir, ri, tolerance
      integer :: n

      allocate(d1_d(nvalues), d2_d(nvalues), result_d(nvalues))
      d1_d = d1_h
      d2_d = d2_h
      result_d = 0.

      call evaluate_rlim_cuda<<<1,32>>>(nvalues, d1_d, d2_d, result_d)
      call checkCUDA(cudaGetLastError(), 'Kappa limiter self-test launch')
      call checkCUDA(cudaDeviceSynchronize(), 'Kappa limiter self-test synchronization')
      result_h = result_d

      do n = 1, nvalues
         ri = (d2_h(n) + eps1)/(d1_h(n) + eps1)
         phir = max(0., min(2.*ri, min(1./3. + 2./3.*ri, 2.)))
         expected(n) = 0.5*phir*d1_h(n)
      end do
      tolerance = 128.*epsilon(1.)*max(1., maxval(abs(expected)))
      if (maxval(abs(result_h - expected)) > tolerance) then
         call fail_cuda_selftest('Kappa limiter')
      end if

      deallocate(result_d, d2_d, d1_d)
   end subroutine test_kappa_limiter

   !> Exercise the otherwise-unselected first-order scalar CUDA kernel.
   subroutine test_upwind_scalar_advection
      implicit none

      real, device, allocatable :: input_d(:, :, :), output_d(:, :, :)
      real, allocatable :: input_h(:, :, :), output_h(:, :, :), expected(:, :, :)
      real :: tolerance
      integer :: i, j, k

      ! The extended-grid metrics are allocated for Kappa or upwind scalar cases.
      if (.not. allocated(dxfci_d)) return

      allocate(input_d(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb-khc:ke+khc))
      allocate(output_d(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb:ke+khc))
      allocate(input_h(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb-khc:ke+khc))
      allocate(output_h(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb:ke+khc))
      allocate(expected(ib:ie, jb:je, kb:ke))

      do k = kb-khc, ke+khc
         do j = jb-jhc, je+jhc
            do i = ib-ihc, ie+ihc
               input_h(i,j,k) = real(i + 2*j + 3*k)
            end do
         end do
      end do
      input_d = input_h
      output_d = 0.
      u0_d = 1.
      v0_d = -0.5
      w0_d = 0.25

      call advecc_upw_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, input_d, output_d)
      call checkCUDA(cudaGetLastError(), 'upwind scalar self-test launch')
      call checkCUDA(cudaDeviceSynchronize(), 'upwind scalar self-test synchronization')
      output_h = output_d

      do k = kb, ke
         do j = jb, je
            do i = ib, ie
               expected(i,j,k) = &
                  -(input_h(i,j,k) - input_h(i-1,j,k))*dxfci(i) &
                  +0.5*(input_h(i,j+1,k) - input_h(i,j,k))*dyi &
                  -0.25*(input_h(i,j,k) - input_h(i,j,k-1))*dzfci(k)
            end do
         end do
      end do
      tolerance = 256.*epsilon(1.)*max(1., maxval(abs(expected)))
      if (maxval(abs(output_h(ib:ie,jb:je,kb:ke) - expected)) > tolerance) then
         call fail_cuda_selftest('upwind scalar advection')
      end if

      ! The normal first timestep refreshes these arrays too; restoring them here
      ! keeps the self-test independent of that ordering contract.
      u0_d = u0
      v0_d = v0
      w0_d = w0
      deallocate(expected, output_h, input_h, output_d, input_d)
   end subroutine test_upwind_scalar_advection

   !> Verify both conversion kernels used by Kappa temperature advection.
   subroutine test_temperature_tendency_copy
      implicit none

      real, allocatable :: thlp_h(:, :, :), thlpc_h(:, :, :)
      real :: expected, tolerance
      integer :: i, j, k

      if (.not. allocated(thlpc_d)) return

      allocate(thlp_h(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
      allocate(thlpc_h(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb:ke+khc))
      thlp_h = 0.
      do k = kb, ke
         do j = jb, je
            do i = ib, ie
               thlp_h(i,j,k) = real(i + 10*j + 100*k)/1000.
            end do
         end do
      end do
      thlp_d = thlp_h
      thlpc_d = -999.

      call thlptothlpc_cuda<<<griddim,blockdim>>>
      call checkCUDA(cudaGetLastError(), 'temperature copy-to-extended self-test launch')
      call checkCUDA(cudaDeviceSynchronize(), 'temperature copy-to-extended self-test synchronization')
      thlpc_h = thlpc_d
      tolerance = 64.*epsilon(1.)
      do k = kb, ke
         do j = jb, je
            do i = ib, ie
               expected = real(i + 10*j + 100*k)/1000.
               if (abs(thlpc_h(i,j,k) - expected) > tolerance) then
                  call fail_cuda_selftest('temperature copy to extended grid')
               end if
            end do
         end do
      end do

      thlp_d = -777.
      call thlpctothlp_cuda<<<griddim,blockdim>>>
      call checkCUDA(cudaGetLastError(), 'temperature copy-from-extended self-test launch')
      call checkCUDA(cudaDeviceSynchronize(), 'temperature copy-from-extended self-test synchronization')
      thlp_h = thlp_d
      do k = kb, ke
         do j = jb, je
            do i = ib, ie
               expected = real(i + 10*j + 100*k)/1000.
               if (abs(thlp_h(i,j,k) - expected) > tolerance) then
                  call fail_cuda_selftest('temperature copy from extended grid')
               end if
            end do
         end do
      end do

      deallocate(thlpc_h, thlp_h)
   end subroutine test_temperature_tendency_copy

   !> Test the driver-inlet boundary kernel without requiring a driver file.
   subroutine test_driver_inlet_boundary
      implicit none

      real, allocatable :: pup_h(:, :, :), up_h(:, :, :)
      real, parameter :: driver_value = 3.25, rk3coefi_test = 2.
      real :: tolerance
      integer :: j, k

      allocate(pup_h(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
      allocate(up_h(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
      u0driver_d = driver_value
      pup_d = -999.
      up_d = -999.

      call bcpup_pup_BCxm_driver_cuda<<<griddim,blockdim>>>( &
         .true., .false., rk3coefi_test, 1., 0.)
      call checkCUDA(cudaGetLastError(), 'driver inlet boundary self-test launch')
      call checkCUDA(cudaDeviceSynchronize(), 'driver inlet boundary self-test synchronization')
      pup_h = pup_d
      up_h = up_d
      tolerance = 64.*epsilon(1.)
      do k = kb, ke
         do j = jb-1, je+1
            if (abs(pup_h(ib,j,k) - driver_value*rk3coefi_test) > tolerance) then
               call fail_cuda_selftest('driver inlet boundary tendency')
            end if
            if (abs(up_h(ib,j,k)) > tolerance) then
               call fail_cuda_selftest('driver inlet boundary increment')
            end if
         end do
      end do

      u0driver_d = u0driver
      deallocate(up_h, pup_h)
   end subroutine test_driver_inlet_boundary

   !> Evaluate the device-only limiter for a compact vector of test inputs.
   attributes(global) subroutine evaluate_rlim_cuda(nvalues, d1, d2, result)
      implicit none

      integer, value, intent(in) :: nvalues
      real, intent(in) :: d1(nvalues), d2(nvalues)
      real, intent(out) :: result(nvalues)
      integer :: n

      n = (blockIdx%x - 1)*blockDim%x + threadIdx%x
      if (n <= nvalues) result(n) = rlim_cuda(d1(n), d2(n))
   end subroutine evaluate_rlim_cuda

   subroutine fail_cuda_selftest(name)
      implicit none

      character(len=*), intent(in) :: name

      write(*,'(A,A,A,I0)') 'CUDA device self-tests failed: ', trim(name), '. rank=', myid
      error stop 1
   end subroutine fail_cuda_selftest
#else
   implicit none
#endif

end module tests_cuda
