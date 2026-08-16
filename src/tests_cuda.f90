!> \file tests_cuda.f90
!! Opt-in CUDA tests driven by the external GPU integration-test harness.

module tests_cuda
#if defined(_GPU) && defined(UDALES_DEBUG)
   use cudafor
   use modadvection, only : advecc_kappa_reset_cuda, advecc_kappa_ducdx_cuda, &
                            advecc_kappa_dvcdy_cuda, advecc_kappa_dwcdz_cuda, &
                            advecc_kappa_add_cuda, advecc_upw_cuda, rlim_cuda
   use modboundary,  only : bcpup_pup_BCxm_driver_cuda, &
                            xm_periodic_device, xT_periodic_device, &
                            xq_periodic_device, xs_periodic_device, &
                            ym_periodic_device, yT_periodic_device, &
                            yq_periodic_device, ys_periodic_device
   use modcuda,  only : blockdim, griddim, checkCUDA, initfield, &
                        thlptothlpc_cuda, thlpctothlp_cuda, &
                        dxhci_d, dxfc_d, dxfci_d, dzhci_d, dzfc_d, dzfci_d, &
                        u0_d, v0_d, w0_d, um_d, vm_d, wm_d, e120_d, e12m_d, &
                        thl0_d, thlm_d, thl0c_d, qt0_d, qtm_d, sv0_d, svm_d, &
                        thlp_d, thlpc_d, pup_d, up_d, u0driver_d, wp_d
   use modfields, only : u0, v0, w0, um, vm, wm, e120, e12m, &
                         thl0, thlm, thl0c, qt0, qtm, sv0, svm, thlp, wp
   use modglobal, only : ib, ie, jb, je, kb, ke, ih, jh, kh, nsv, &
                         ihc, jhc, khc, dxhci, dxfc, dxfci, dxi, dyi, &
                         dzhci, dzfc, dzfci, dzfi, eps1, &
                         lheatpump, lfan_hp, nhppoints, ltempeq
   use modheatpump, only : heatpump, nhppoints_local, idhppts_local_d, &
                           thl_dot_hp, w_hp_exhaust
   use modinletdata, only : u0driver
   use modmpi,   only : myid

   implicit none
   private

   public :: run_cuda_selftests_if_requested

   real, parameter :: periodic_halo_sentinel = -12345678.

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
         call test_kappa_scalar_advection
         call test_upwind_scalar_advection
         call test_temperature_tendency_copy
         call test_driver_inlet_boundary
         call test_periodic_device_halos
         call test_heatpump_scatter
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

   !> Check the directional Kappa kernels against the host flux-divergence algebra.
   subroutine test_kappa_scalar_advection
      implicit none

      real, device, allocatable :: input_d(:, :, :), output_d(:, :, :)
      real, allocatable :: input_h(:, :, :), output_initial_h(:, :, :), output_h(:, :, :)
      real, allocatable :: expected(:, :, :), u_test_h(:, :, :), v_test_h(:, :, :), w_test_h(:, :, :)
      real :: cf, d1, d2, flux_upper, flux_lower, tendency, tolerance
      integer :: i, j, k

      if (.not. allocated(dxhci_d) .or. .not. allocated(dxfc_d) .or. &
          .not. allocated(dzhci_d) .or. .not. allocated(dzfc_d)) return

      allocate(input_d(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb-khc:ke+khc))
      allocate(output_d(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb:ke+khc))
      allocate(input_h(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb-khc:ke+khc))
      allocate(output_initial_h(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb:ke+khc))
      allocate(output_h(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb:ke+khc))
      allocate(expected(ib:ie, jb:je, kb:ke))
      allocate(u_test_h(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
      allocate(v_test_h(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
      allocate(w_test_h(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))

      do k = kb-khc, ke+khc
         do j = jb-jhc, je+jhc
            do i = ib-ihc, ie+ihc
               input_h(i,j,k) = real(13*i*i + 7*j*j - 5*k*k + 3*i*j - 2*j*k)/1000.
            end do
         end do
      end do
      do k = kb, ke+khc
         do j = jb-jhc, je+jhc
            do i = ib-ihc, ie+ihc
               output_initial_h(i,j,k) = real(i - 2*j + 3*k)/10000.
            end do
         end do
      end do
      do k = kb-kh, ke+kh
         do j = jb-jh, je+jh
            do i = ib-ih, ie+ih
               if (mod(i + j + k, 2) == 0) then
                  u_test_h(i,j,k) = 0.75
                  v_test_h(i,j,k) = -0.5
                  w_test_h(i,j,k) = 0.25
               else
                  u_test_h(i,j,k) = -0.625
                  v_test_h(i,j,k) = 0.375
                  w_test_h(i,j,k) = -0.875
               end if
            end do
         end do
      end do

      input_d = input_h
      output_d = output_initial_h
      u0_d = u_test_h
      v0_d = v_test_h
      w0_d = w_test_h

      call advecc_kappa_reset_cuda<<<griddim,blockdim>>>(ihc, jhc, khc)
      call checkCUDA(cudaGetLastError(), 'Kappa X reset self-test launch')
      call advecc_kappa_ducdx_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, input_d)
      call checkCUDA(cudaGetLastError(), 'Kappa X flux self-test launch')
      call advecc_kappa_add_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, output_d)
      call checkCUDA(cudaGetLastError(), 'Kappa X add self-test launch')

      call advecc_kappa_reset_cuda<<<griddim,blockdim>>>(ihc, jhc, khc)
      call checkCUDA(cudaGetLastError(), 'Kappa Y reset self-test launch')
      call advecc_kappa_dvcdy_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, input_d)
      call checkCUDA(cudaGetLastError(), 'Kappa Y flux self-test launch')
      call advecc_kappa_add_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, output_d)
      call checkCUDA(cudaGetLastError(), 'Kappa Y add self-test launch')

      call advecc_kappa_reset_cuda<<<griddim,blockdim>>>(ihc, jhc, khc)
      call checkCUDA(cudaGetLastError(), 'Kappa Z reset self-test launch')
      call advecc_kappa_dwcdz_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, input_d)
      call checkCUDA(cudaGetLastError(), 'Kappa Z flux self-test launch')
      call advecc_kappa_add_cuda<<<griddim,blockdim>>>(ihc, jhc, khc, output_d)
      call checkCUDA(cudaGetLastError(), 'Kappa Z add self-test launch')
      call checkCUDA(cudaDeviceSynchronize(), 'Kappa advection self-test synchronization')
      output_h = output_d

      do k = kb, ke
         do j = jb, je
            do i = ib, ie
               tendency = output_initial_h(i,j,k)

               if (u_test_h(i+1,j,k) > 0.) then
                  d1 = (input_h(i,j,k) - input_h(i-1,j,k))*dxhci(i)
                  d2 = (input_h(i+1,j,k) - input_h(i,j,k))*dxhci(i+1)
                  cf = input_h(i,j,k)
               else
                  d1 = (input_h(i+1,j,k) - input_h(i+2,j,k))*dxhci(i+2)
                  d2 = (input_h(i,j,k) - input_h(i+1,j,k))*dxhci(i+1)
                  cf = input_h(i+1,j,k)
               end if
               cf = cf + dxfc(i+1)*kappa_rlim_host(d1, d2)
               flux_upper = -cf*u_test_h(i+1,j,k)*dxfci(i)

               if (u_test_h(i,j,k) > 0.) then
                  d1 = (input_h(i-1,j,k) - input_h(i-2,j,k))*dxhci(i-1)
                  d2 = (input_h(i,j,k) - input_h(i-1,j,k))*dxhci(i)
                  cf = input_h(i-1,j,k)
               else
                  d1 = (input_h(i,j,k) - input_h(i+1,j,k))*dxhci(i+1)
                  d2 = (input_h(i-1,j,k) - input_h(i,j,k))*dxhci(i)
                  cf = input_h(i,j,k)
               end if
               cf = cf + dxfc(i)*kappa_rlim_host(d1, d2)
               flux_lower = cf*u_test_h(i,j,k)*dxfci(i)
               tendency = tendency + flux_upper + flux_lower

               if (v_test_h(i,j+1,k) > 0.) then
                  d1 = input_h(i,j,k) - input_h(i,j-1,k)
                  d2 = input_h(i,j+1,k) - input_h(i,j,k)
                  cf = input_h(i,j,k)
               else
                  d1 = input_h(i,j+1,k) - input_h(i,j+2,k)
                  d2 = input_h(i,j,k) - input_h(i,j+1,k)
                  cf = input_h(i,j+1,k)
               end if
               cf = cf + kappa_rlim_host(d1, d2)
               flux_upper = -cf*v_test_h(i,j+1,k)*dyi

               if (v_test_h(i,j,k) > 0.) then
                  d1 = input_h(i,j-1,k) - input_h(i,j-2,k)
                  d2 = input_h(i,j,k) - input_h(i,j-1,k)
                  cf = input_h(i,j-1,k)
               else
                  d1 = input_h(i,j,k) - input_h(i,j+1,k)
                  d2 = input_h(i,j-1,k) - input_h(i,j,k)
                  cf = input_h(i,j,k)
               end if
               cf = cf + kappa_rlim_host(d1, d2)
               flux_lower = cf*v_test_h(i,j,k)*dyi
               tendency = tendency + flux_upper + flux_lower

               if (w_test_h(i,j,k+1) > 0.) then
                  d1 = (input_h(i,j,k) - input_h(i,j,k-1))*dzhci(k)
                  d2 = (input_h(i,j,k+1) - input_h(i,j,k))*dzhci(k+1)
                  cf = input_h(i,j,k)
               else
                  d1 = (input_h(i,j,k+1) - input_h(i,j,k+2))*dzhci(k+2)
                  d2 = (input_h(i,j,k) - input_h(i,j,k+1))*dzhci(k+1)
                  cf = input_h(i,j,k+1)
               end if
               cf = cf + dzfc(k+1)*kappa_rlim_host(d1, d2)
               flux_upper = -cf*w_test_h(i,j,k+1)*dzfci(k)

               flux_lower = 0.
               if (k > kb) then
                  if (w_test_h(i,j,k) > 0.) then
                     d1 = (input_h(i,j,k-1) - input_h(i,j,k-2))*dzhci(k-1)
                     d2 = (input_h(i,j,k) - input_h(i,j,k-1))*dzhci(k)
                     cf = input_h(i,j,k-1)
                  else
                     d1 = (input_h(i,j,k) - input_h(i,j,k+1))*dzhci(k+1)
                     d2 = (input_h(i,j,k-1) - input_h(i,j,k))*dzhci(k)
                     cf = input_h(i,j,k)
                  end if
                  cf = cf + dzfc(k)*kappa_rlim_host(d1, d2)
                  flux_lower = cf*w_test_h(i,j,k)*dzfci(k)
               end if
               expected(i,j,k) = tendency + flux_upper + flux_lower
            end do
         end do
      end do

      tolerance = 4096.*epsilon(1.)*max(1., maxval(abs(expected)))
      if (maxval(abs(output_h(ib:ie,jb:je,kb:ke) - expected)) > tolerance) then
         call fail_cuda_selftest('Kappa scalar advection')
      end if

      u0_d = u0
      v0_d = v0
      w0_d = w0
      deallocate(w_test_h, v_test_h, u_test_h, expected, output_h, output_initial_h, input_h)
      deallocate(output_d, input_d)
   end subroutine test_kappa_scalar_advection

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

   !> Check every periodic device routine against independently generated halo values.
   subroutine test_periodic_device_halos
      implicit none

      call test_x_momentum_periodic
      call test_y_momentum_periodic
      call test_x_temperature_periodic
      call test_y_temperature_periodic
      call test_x_moisture_periodic
      call test_y_moisture_periodic
      call test_x_scalar_periodic
      call test_y_scalar_periodic

      call restore_periodic_device_fields
      call checkCUDA(cudaDeviceSynchronize(), &
                     'periodic halo self-test field restoration')
   end subroutine test_periodic_device_halos

   subroutine test_x_momentum_periodic
      implicit none

      call prepare_standard_device(u0_d,   10., 'x')
      call prepare_standard_device(um_d,   20., 'x')
      call prepare_standard_device(v0_d,   30., 'x')
      call prepare_standard_device(vm_d,   40., 'x')
      call prepare_standard_device(w0_d,   50., 'x')
      call prepare_standard_device(wm_d,   60., 'x')
      if (allocated(e120_d)) then
         call prepare_standard_device(e120_d, 70., 'x')
         call prepare_standard_device(e12m_d, 80., 'x')
      end if

      call xm_periodic_device
      call checkCUDA(cudaDeviceSynchronize(), 'x momentum periodic halo self-test')

      call assert_standard_periodic(u0_d,   10., 'x', 'x-periodic u0')
      call assert_standard_periodic(um_d,   20., 'x', 'x-periodic um')
      call assert_standard_periodic(v0_d,   30., 'x', 'x-periodic v0')
      call assert_standard_periodic(vm_d,   40., 'x', 'x-periodic vm')
      call assert_standard_periodic(w0_d,   50., 'x', 'x-periodic w0')
      call assert_standard_periodic(wm_d,   60., 'x', 'x-periodic wm')
      if (allocated(e120_d)) then
         call assert_standard_periodic(e120_d, 70., 'x', 'x-periodic e120')
         call assert_standard_periodic(e12m_d, 80., 'x', 'x-periodic e12m')
      end if
   end subroutine test_x_momentum_periodic

   subroutine test_y_momentum_periodic
      implicit none

      call prepare_standard_device(u0_d,   10., 'y')
      call prepare_standard_device(um_d,   20., 'y')
      call prepare_standard_device(v0_d,   30., 'y')
      call prepare_standard_device(vm_d,   40., 'y')
      call prepare_standard_device(w0_d,   50., 'y')
      call prepare_standard_device(wm_d,   60., 'y')
      if (allocated(e120_d)) then
         call prepare_standard_device(e120_d, 70., 'y')
         call prepare_standard_device(e12m_d, 80., 'y')
      end if

      call ym_periodic_device
      call checkCUDA(cudaDeviceSynchronize(), 'y momentum periodic halo self-test')

      call assert_standard_periodic(u0_d,   10., 'y', 'y-periodic u0')
      call assert_standard_periodic(um_d,   20., 'y', 'y-periodic um')
      call assert_standard_periodic(v0_d,   30., 'y', 'y-periodic v0')
      call assert_standard_periodic(vm_d,   40., 'y', 'y-periodic vm')
      call assert_standard_periodic(w0_d,   50., 'y', 'y-periodic w0')
      call assert_standard_periodic(wm_d,   60., 'y', 'y-periodic wm')
      if (allocated(e120_d)) then
         call assert_standard_periodic(e120_d, 70., 'y', 'y-periodic e120')
         call assert_standard_periodic(e12m_d, 80., 'y', 'y-periodic e12m')
      end if
   end subroutine test_y_momentum_periodic

   subroutine test_x_temperature_periodic
      implicit none

      if (.not. allocated(thl0_d)) return

      call prepare_standard_device(thl0_d, 110., 'x')
      call prepare_standard_device(thlm_d, 120., 'x')
      if (allocated(thl0c_d)) call prepare_extended_device(thl0c_d, 130., 'x')

      call xT_periodic_device
      call checkCUDA(cudaDeviceSynchronize(), 'x temperature periodic halo self-test')

      call assert_standard_periodic(thl0_d, 110., 'x', 'x-periodic thl0')
      call assert_standard_periodic(thlm_d, 120., 'x', 'x-periodic thlm')
      if (allocated(thl0c_d)) then
         call assert_extended_periodic(thl0c_d, 130., 'x', 'x-periodic thl0c')
      end if
   end subroutine test_x_temperature_periodic

   subroutine test_y_temperature_periodic
      implicit none

      if (.not. allocated(thl0_d)) return

      call prepare_standard_device(thl0_d, 110., 'y')
      call prepare_standard_device(thlm_d, 120., 'y')
      if (allocated(thl0c_d)) call prepare_extended_device(thl0c_d, 130., 'y')

      call yT_periodic_device
      call checkCUDA(cudaDeviceSynchronize(), 'y temperature periodic halo self-test')

      call assert_standard_periodic(thl0_d, 110., 'y', 'y-periodic thl0')
      call assert_standard_periodic(thlm_d, 120., 'y', 'y-periodic thlm')
      if (allocated(thl0c_d)) then
         call assert_extended_periodic(thl0c_d, 130., 'y', 'y-periodic thl0c')
      end if
   end subroutine test_y_temperature_periodic

   subroutine test_x_moisture_periodic
      implicit none

      if (.not. allocated(qt0_d)) return

      call prepare_standard_device(qt0_d, 210., 'x')
      call prepare_standard_device(qtm_d, 220., 'x')
      call xq_periodic_device
      call checkCUDA(cudaDeviceSynchronize(), 'x moisture periodic halo self-test')
      call assert_standard_periodic(qt0_d, 210., 'x', 'x-periodic qt0')
      call assert_standard_periodic(qtm_d, 220., 'x', 'x-periodic qtm')
   end subroutine test_x_moisture_periodic

   subroutine test_y_moisture_periodic
      implicit none

      if (.not. allocated(qt0_d)) return

      call prepare_standard_device(qt0_d, 210., 'y')
      call prepare_standard_device(qtm_d, 220., 'y')
      call yq_periodic_device
      call checkCUDA(cudaDeviceSynchronize(), 'y moisture periodic halo self-test')
      call assert_standard_periodic(qt0_d, 210., 'y', 'y-periodic qt0')
      call assert_standard_periodic(qtm_d, 220., 'y', 'y-periodic qtm')
   end subroutine test_y_moisture_periodic

   subroutine test_x_scalar_periodic
      implicit none

      if (.not. allocated(sv0_d)) return

      call prepare_scalar_device(sv0_d, 310., 'x')
      call prepare_scalar_device(svm_d, 320., 'x')
      call xs_periodic_device
      call checkCUDA(cudaDeviceSynchronize(), 'x scalar periodic halo self-test')
      call assert_scalar_periodic(sv0_d, 310., 'x', 'x-periodic sv0')
      call assert_scalar_periodic(svm_d, 320., 'x', 'x-periodic svm')
   end subroutine test_x_scalar_periodic

   subroutine test_y_scalar_periodic
      implicit none

      if (.not. allocated(sv0_d)) return

      call prepare_scalar_device(sv0_d, 310., 'y')
      call prepare_scalar_device(svm_d, 320., 'y')
      call ys_periodic_device
      call checkCUDA(cudaDeviceSynchronize(), 'y scalar periodic halo self-test')
      call assert_scalar_periodic(sv0_d, 310., 'y', 'y-periodic sv0')
      call assert_scalar_periodic(svm_d, 320., 'y', 'y-periodic svm')
   end subroutine test_y_scalar_periodic

   subroutine prepare_standard_device(field_d, offset, direction)
      implicit none

      real, device, intent(out) :: field_d(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
      real, intent(in) :: offset
      character(len=1), intent(in) :: direction
      real, allocatable :: field_h(:, :, :)
      integer i, j, k

      allocate(field_h(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
      do k = kb-kh, ke+kh
         do j = jb-jh, je+jh
            do i = ib-ih, ie+ih
               field_h(i,j,k) = periodic_value(i, j, k, 0, offset)
            end do
         end do
      end do
      select case (direction)
      case ('x')
         field_h(ib-ih:ib-1,:,:) = periodic_halo_sentinel
         field_h(ie+1:ie+ih,:,:) = periodic_halo_sentinel
      case ('y')
         field_h(:,jb-jh:jb-1,:) = periodic_halo_sentinel
         field_h(:,je+1:je+jh,:) = periodic_halo_sentinel
      case default
         call fail_cuda_selftest('invalid standard periodic direction')
      end select
      field_d = field_h
      deallocate(field_h)
   end subroutine prepare_standard_device

   subroutine prepare_extended_device(field_d, offset, direction)
      implicit none

      real, device, intent(out) :: field_d(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb-khc:ke+khc)
      real, intent(in) :: offset
      character(len=1), intent(in) :: direction
      real, allocatable :: field_h(:, :, :)
      integer i, j, k

      allocate(field_h(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb-khc:ke+khc))
      do k = kb-khc, ke+khc
         do j = jb-jhc, je+jhc
            do i = ib-ihc, ie+ihc
               field_h(i,j,k) = periodic_value(i, j, k, 0, offset)
            end do
         end do
      end do
      select case (direction)
      case ('x')
         field_h(ib-ihc:ib-1,:,:) = periodic_halo_sentinel
         field_h(ie+1:ie+ihc,:,:) = periodic_halo_sentinel
      case ('y')
         field_h(:,jb-jhc:jb-1,:) = periodic_halo_sentinel
         field_h(:,je+1:je+jhc,:) = periodic_halo_sentinel
      case default
         call fail_cuda_selftest('invalid extended periodic direction')
      end select
      field_d = field_h
      deallocate(field_h)
   end subroutine prepare_extended_device

   subroutine prepare_scalar_device(field_d, offset, direction)
      implicit none

      real, device, intent(out) :: field_d(ib-ihc:ie+ihc, jb-jhc:je+jhc, &
                                             kb-khc:ke+khc, 1:nsv)
      real, intent(in) :: offset
      character(len=1), intent(in) :: direction
      real, allocatable :: field_h(:, :, :, :)
      integer i, j, k, n

      allocate(field_h(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb-khc:ke+khc, nsv))
      do n = 1, nsv
         do k = kb-khc, ke+khc
            do j = jb-jhc, je+jhc
               do i = ib-ihc, ie+ihc
                  field_h(i,j,k,n) = periodic_value(i, j, k, n, offset)
               end do
            end do
         end do
      end do
      select case (direction)
      case ('x')
         field_h(ib-ihc:ib-1,:,:,:) = periodic_halo_sentinel
         field_h(ie+1:ie+ihc,:,:,:) = periodic_halo_sentinel
      case ('y')
         field_h(:,jb-jhc:jb-1,:,:) = periodic_halo_sentinel
         field_h(:,je+1:je+jhc,:,:) = periodic_halo_sentinel
      case default
         call fail_cuda_selftest('invalid scalar periodic direction')
      end select
      field_d = field_h
      deallocate(field_h)
   end subroutine prepare_scalar_device

   subroutine assert_standard_periodic(field_d, offset, direction, name)
      implicit none

      real, device, intent(in) :: field_d(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
      real, intent(in) :: offset
      character(len=1), intent(in) :: direction
      character(len=*), intent(in) :: name
      real, allocatable :: field_h(:, :, :)
      real expected
      integer i, j, k, source_i, source_j

      allocate(field_h(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
      field_h = field_d
      do k = kb-kh, ke+kh
         do j = jb-jh, je+jh
            do i = ib-ih, ie+ih
               source_i = i
               source_j = j
               call periodic_source_index(i, j, ih, jh, direction, source_i, source_j)
               expected = periodic_value(source_i, source_j, k, 0, offset)
               if (field_h(i,j,k) /= expected) call fail_cuda_selftest(name)
            end do
         end do
      end do
      deallocate(field_h)
   end subroutine assert_standard_periodic

   subroutine assert_extended_periodic(field_d, offset, direction, name)
      implicit none

      real, device, intent(in) :: field_d(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb-khc:ke+khc)
      real, intent(in) :: offset
      character(len=1), intent(in) :: direction
      character(len=*), intent(in) :: name
      real, allocatable :: field_h(:, :, :)
      real expected
      integer i, j, k, source_i, source_j

      allocate(field_h(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb-khc:ke+khc))
      field_h = field_d
      do k = kb-khc, ke+khc
         do j = jb-jhc, je+jhc
            do i = ib-ihc, ie+ihc
               source_i = i
               source_j = j
               call periodic_source_index(i, j, ihc, jhc, direction, source_i, source_j)
               expected = periodic_value(source_i, source_j, k, 0, offset)
               if (field_h(i,j,k) /= expected) call fail_cuda_selftest(name)
            end do
         end do
      end do
      deallocate(field_h)
   end subroutine assert_extended_periodic

   subroutine assert_scalar_periodic(field_d, offset, direction, name)
      implicit none

      real, device, intent(in) :: field_d(ib-ihc:ie+ihc, jb-jhc:je+jhc, &
                                            kb-khc:ke+khc, 1:nsv)
      real, intent(in) :: offset
      character(len=1), intent(in) :: direction
      character(len=*), intent(in) :: name
      real, allocatable :: field_h(:, :, :, :)
      real expected
      integer i, j, k, n, source_i, source_j

      allocate(field_h(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb-khc:ke+khc, nsv))
      field_h = field_d
      do n = 1, nsv
         do k = kb-khc, ke+khc
            do j = jb-jhc, je+jhc
               do i = ib-ihc, ie+ihc
                  source_i = i
                  source_j = j
                  call periodic_source_index(i, j, ihc, jhc, direction, source_i, source_j)
                  expected = periodic_value(source_i, source_j, k, n, offset)
                  if (field_h(i,j,k,n) /= expected) call fail_cuda_selftest(name)
               end do
            end do
         end do
      end do
      deallocate(field_h)
   end subroutine assert_scalar_periodic

   subroutine periodic_source_index(i, j, halo_i, halo_j, direction, source_i, source_j)
      implicit none

      integer, intent(in) :: i, j, halo_i, halo_j
      character(len=1), intent(in) :: direction
      integer, intent(inout) :: source_i, source_j

      select case (direction)
      case ('x')
         if (i < ib) source_i = ie + 1 - (ib - i)
         if (i > ie) source_i = ib - 1 + (i - ie)
      case ('y')
         if (j < jb) source_j = je + 1 - (jb - j)
         if (j > je) source_j = jb - 1 + (j - je)
      case default
         call fail_cuda_selftest('invalid periodic assertion direction')
      end select

      ! Ensure callers passed the halo depth matching the tested array shape.
      if ((direction == 'x') .and. ((i < ib-halo_i) .or. (i > ie+halo_i))) then
         call fail_cuda_selftest('invalid x periodic halo extent')
      end if
      if ((direction == 'y') .and. ((j < jb-halo_j) .or. (j > je+halo_j))) then
         call fail_cuda_selftest('invalid y periodic halo extent')
      end if
   end subroutine periodic_source_index

   real function periodic_value(i, j, k, n, offset)
      implicit none

      integer, intent(in) :: i, j, k, n
      real, intent(in) :: offset

      periodic_value = offset + real(i + 101*j + 10007*k + 1000003*n)
   end function periodic_value

   subroutine restore_periodic_device_fields
      implicit none

      u0_d = u0
      v0_d = v0
      w0_d = w0
      um_d = um
      vm_d = vm
      wm_d = wm
      if (allocated(e120_d)) then
         e120_d = e120
         e12m_d = e12m
      end if
      if (allocated(thl0_d)) then
         thl0_d = thl0
         thlm_d = thlm
      end if
      if (allocated(thl0c_d)) thl0c_d = thl0c
      if (allocated(qt0_d)) then
         qt0_d = qt0
         qtm_d = qtm
      end if
      if (allocated(sv0_d)) then
         sv0_d = sv0
         svm_d = svm
      end if
   end subroutine restore_periodic_device_fields

   !> Verify the heatpump OpenACC kernel writes the correct cells and no others.
   subroutine test_heatpump_scatter
      implicit none

      integer, parameter :: test_i = 3, test_j = 3, test_k = 3
      real,    parameter :: initial_thlp = 5., thl_dot = 0.25, w_exhaust = 0.75

      logical :: saved_lheatpump, saved_ltempeq, saved_lfan_hp
      integer :: saved_nhppoints, saved_nhppoints_local
      real    :: saved_thl_dot, saved_w_exhaust
      real, allocatable :: thlp_h(:,:,:), wm_h(:,:,:), w0_h(:,:,:), wp_h(:,:,:)
      real    :: expected_thlp, tolerance
      integer :: i, j, k

      if (.not. allocated(thlp_d)) return

      saved_lheatpump       = lheatpump
      saved_ltempeq         = ltempeq
      saved_lfan_hp         = lfan_hp
      saved_nhppoints       = nhppoints
      saved_nhppoints_local = nhppoints_local
      saved_thl_dot         = thl_dot_hp
      saved_w_exhaust       = w_hp_exhaust

      lheatpump       = .true.
      ltempeq         = .true.
      lfan_hp         = .true.
      nhppoints       = 1
      nhppoints_local = 1
      thl_dot_hp      = thl_dot
      w_hp_exhaust    = w_exhaust
      if (allocated(idhppts_local_d)) deallocate(idhppts_local_d)
      allocate(idhppts_local_d(1, 3))
      idhppts_local_d(1, 1) = test_i
      idhppts_local_d(1, 2) = test_j
      idhppts_local_d(1, 3) = test_k

      thlp_d = initial_thlp
      wm_d   = 0.
      w0_d   = 0.
      wp_d   = 0.

      call heatpump
      call checkCUDA(cudaDeviceSynchronize(), 'heatpump scatter self-test synchronization')

      allocate(thlp_h(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
      allocate(wm_h  (ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
      allocate(w0_h  (ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
      allocate(wp_h  (ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
      thlp_h = thlp_d
      wm_h   = wm_d
      w0_h   = w0_d
      wp_h   = wp_d

      expected_thlp = initial_thlp - thl_dot * dxi * dyi * dzfi(test_k)
      tolerance = 64. * epsilon(1.) * max(1., abs(expected_thlp))
      if (abs(thlp_h(test_i, test_j, test_k) - expected_thlp) > tolerance) then
         call fail_cuda_selftest('heatpump temperature tendency')
      end if
      do k = kb, ke
         do j = jb, je
            do i = ib, ie
               if (i == test_i .and. j == test_j .and. k == test_k) cycle
               if (abs(thlp_h(i,j,k) - initial_thlp) > tolerance) then
                  call fail_cuda_selftest('heatpump temperature tendency side-effect')
               end if
            end do
         end do
      end do

      tolerance = 64. * epsilon(1.) * max(1., abs(w_exhaust))
      if (abs(wm_h(test_i, test_j, test_k+1) - w_exhaust) > tolerance) then
         call fail_cuda_selftest('heatpump fan wm')
      end if
      if (abs(w0_h(test_i, test_j, test_k+1) - w_exhaust) > tolerance) then
         call fail_cuda_selftest('heatpump fan w0')
      end if
      if (abs(wp_h(test_i, test_j, test_k+1)) > tolerance) then
         call fail_cuda_selftest('heatpump fan wp zeroed')
      end if

      deallocate(wp_h, w0_h, wm_h, thlp_h)
      deallocate(idhppts_local_d)

      lheatpump       = saved_lheatpump
      ltempeq         = saved_ltempeq
      lfan_hp         = saved_lfan_hp
      nhppoints       = saved_nhppoints
      nhppoints_local = saved_nhppoints_local
      thl_dot_hp      = saved_thl_dot
      w_hp_exhaust    = saved_w_exhaust

      thlp_d = thlp
      wm_d   = wm
      w0_d   = w0
      wp_d   = wp
   end subroutine test_heatpump_scatter

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

   real function kappa_rlim_host(d1, d2)
      implicit none

      real, intent(in) :: d1, d2
      real :: phir, ri

      ri = (d2 + eps1)/(d1 + eps1)
      phir = max(0., min(2.*ri, min(1./3. + 2./3.*ri, 2.)))
      kappa_rlim_host = 0.5*phir*d1
   end function kappa_rlim_host

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
