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
                        thlp_d, thlpc_d, pup_d, up_d, u0driver_d, wp_d, &
                        vp_d, qtp_d, svp_d, u0av_d, v0av_d, thl0av_d, qt0av_d, sv0av_d, &
                        uprof_d, vprof_d, thlprof_d, qtprof_d, svprof_d, &
                        ekm_d, ekh_d, pres0_d, fachf_d, facef_d, updateFacFluxHost, &
                        facT1_d, facqsat_d, fachurel_d, facf_d, updateFacetPropsDevice, &
                        dxdydzfi_d, dxdydzhi_d
   use modfields, only : u0, v0, w0, um, vm, wm, e120, e12m, pres0, &
                         thl0, thlm, thl0c, qt0, qtm, sv0, svm, thlp, wp, &
                         up, vp, qtp, svp, &
                         u0av, v0av, thl0av, qt0av, sv0av, &
                         uprof, vprof, thlprof, qtprof, svprof
   use modglobal, only : ib, ie, jb, je, kb, ke, ih, jh, kh, nsv, &
                         ihc, jhc, khc, dxhci, dxfc, dxfci, dxi, dyi, &
                         dzhci, dzfc, dzfci, dzfi, eps1, &
                         dx, dy, dzf, dzh, dxdydzfi, dxdydzhi, &
                         lheatpump, lfan_hp, nhppoints, ltempeq, &
                         lmoist, lnudge, lnudgevel, tnudge, nnudge, &
                         iwallmom, nfcts, xhat, yhat, zhat, &
                         totheatflux, totqflux, lEB, rk3step, iwallmoist, &
                         lperiodicEBcorr, sinkbase, fraction, xlen, ylen, zh
   use modforces,   only : nudge, periodicEBcorr
   use modheatpump, only : heatpump, nhppoints_local, idhppts_local_d, &
                           thl_dot_hp, w_hp_exhaust
   use initfac,     only : fachf, facef, facT, facqsat, fachurel, facf, faclGR
   use modibm,      only : mask_u, mask_v, mask_w, mask_c, &
                           mask_u_d, mask_v_d, mask_w_d, mask_c_d, &
                           bndpts_u_d, bndpts_v_d, bndpts_w_d, bndpts_c_d, &
                           bound_info_u, bound_info_v, bound_info_w, bound_info_c, &
                           diffu_corr, diffv_corr, diffw_corr, diffc_corr, &
                           diffu_corr_device, diffv_corr_device, diffw_corr_device, &
                           diffc_corr_device, &
                           wallfunmom, wallfunmom_dir_device, fac_tau_d, fac_tau_raw, &
                           bound_info_type, &
                           wallfunheat, wallfunheat_dir_device, bound_info_c, faclGR_d, &
                           fac_htc_raw, fac_cth_raw, fac_pres_raw, fac_pres2_raw, &
                           fac_htc_d, fac_cth_d, fac_pres_d, fac_pres2_d, &
                           check_ibm_section_cache
   use modsubgriddata, only : ekm, ekh
   use modinletdata, only : u0driver
   use modmpi,   only : myid, nprocs

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
         call test_nudge_profiles
         call test_periodic_ebcorr
         call test_ibm_device_geometry
         call test_cell_volume_reciprocals
         call test_ibm_section_cache
         call test_ibm_diff_corr
         call test_ibm_wallfunmom
         call test_ibm_wallfunheat
         call test_facflux_handover
         call test_facet_props_refresh
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
      integer, allocatable :: saved_ids_h(:,:)
      logical :: had_ids
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
      ! The run's real point list is kept in place and only its first row is
      ! borrowed. Neither destroying it nor swapping it for another allocation
      ! is safe: init_heatpump is the only thing that builds it and does not run
      ! again, and re-pointing the module array between kernel launches leaves
      ! the accelerator using the allocation it saw first. Writing through the
      ! existing allocation keeps that identity fixed.
      had_ids = allocated(idhppts_local_d)
      if (had_ids) then
         allocate(saved_ids_h(size(idhppts_local_d,1), size(idhppts_local_d,2)))
         saved_ids_h = idhppts_local_d
      else
         allocate(idhppts_local_d(1, 3))
      end if
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
      if (had_ids) then
         idhppts_local_d = saved_ids_h
         deallocate(saved_ids_h)
      else
         ! Nothing to put back: with no points on this rank nhppoints_local
         ! restores to zero and heatpump returns before reading the list.
         deallocate(idhppts_local_d)
      end if

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

   !> Check the OpenACC nudge kernels against the same algebra on the host.
   !!
   !! Covers the two things most likely to go wrong in this port: that the
   !! relaxation is applied only from kb+nnudge upwards, and that it spans each
   !! array's full declared extent including the halos, matching the whole-array
   !! assignment used by the CPU branch. Cells below the nudging level must be
   !! left untouched.
   subroutine test_nudge_profiles
      implicit none

      real, parameter :: initial_p = 2., av_value = 3., prof_value = 1., t_relax = 4.

      logical :: saved_lnudge, saved_lnudgevel, saved_ltempeq, saved_lmoist
      integer :: saved_nnudge
      real    :: saved_tnudge
      real, allocatable :: up_h(:,:,:), vp_h(:,:,:), thlp_h(:,:,:), qtp_h(:,:,:)
      real, allocatable :: svp_h(:,:,:,:)
      real    :: expected, tolerance
      integer :: i, j, k, n, ktest

      if (.not. allocated(up_d)) return

      ! Need at least one nudged and one un-nudged level to test both.
      if (ke - kb < 2) return

      saved_lnudge    = lnudge
      saved_lnudgevel = lnudgevel
      saved_ltempeq   = ltempeq
      saved_lmoist    = lmoist
      saved_nnudge    = nnudge
      saved_tnudge    = tnudge

      lnudge    = .true.
      lnudgevel = .true.
      nnudge    = 1
      tnudge    = t_relax
      ktest     = kb + nnudge

      up_d      = initial_p
      vp_d      = initial_p
      u0av_d    = av_value
      v0av_d    = av_value
      uprof_d   = prof_value
      vprof_d   = prof_value
      if (ltempeq) then
         thlp_d    = initial_p
         thl0av_d  = av_value
         thlprof_d = prof_value
      end if
      if (lmoist) then
         qtp_d    = initial_p
         qt0av_d  = av_value
         qtprof_d = prof_value
      end if
      if (nsv > 0) then
         svp_d    = initial_p
         sv0av_d  = av_value
         svprof_d = prof_value
      end if

      call nudge
      call checkCUDA(cudaDeviceSynchronize(), 'nudge self-test synchronization')

      expected  = initial_p - (av_value - prof_value) / t_relax
      tolerance = 64. * epsilon(1.) * max(1., abs(expected))

      allocate(up_h (ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
      allocate(vp_h (ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
      up_h = up_d
      vp_h = vp_d

      ! Nudged levels, checked across the full halo extent.
      do k = ktest, ke
         do j = jb-jh, je+jh
            do i = ib-ih, ie+ih
               if (abs(up_h(i,j,k) - expected) > tolerance) call fail_cuda_selftest('nudge up')
               if (abs(vp_h(i,j,k) - expected) > tolerance) call fail_cuda_selftest('nudge vp')
            end do
         end do
      end do

      ! Levels below kb+nnudge must be untouched.
      do k = kb, ktest-1
         do j = jb-jh, je+jh
            do i = ib-ih, ie+ih
               if (abs(up_h(i,j,k) - initial_p) > tolerance) call fail_cuda_selftest('nudge up below level')
               if (abs(vp_h(i,j,k) - initial_p) > tolerance) call fail_cuda_selftest('nudge vp below level')
            end do
         end do
      end do
      deallocate(vp_h, up_h)

      if (ltempeq) then
         allocate(thlp_h(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
         thlp_h = thlp_d
         do k = ktest, ke
            do j = jb-jh, je+jh
               do i = ib-ih, ie+ih
                  if (abs(thlp_h(i,j,k) - expected) > tolerance) call fail_cuda_selftest('nudge thlp')
               end do
            end do
         end do
         do k = kb, ktest-1
            do j = jb-jh, je+jh
               do i = ib-ih, ie+ih
                  if (abs(thlp_h(i,j,k) - initial_p) > tolerance) call fail_cuda_selftest('nudge thlp below level')
               end do
            end do
         end do
         deallocate(thlp_h)
      end if

      if (lmoist) then
         allocate(qtp_h(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
         qtp_h = qtp_d
         do k = ktest, ke
            do j = jb-jh, je+jh
               do i = ib-ih, ie+ih
                  if (abs(qtp_h(i,j,k) - expected) > tolerance) call fail_cuda_selftest('nudge qtp')
               end do
            end do
         end do
         do k = kb, ktest-1
            do j = jb-jh, je+jh
               do i = ib-ih, ie+ih
                  if (abs(qtp_h(i,j,k) - initial_p) > tolerance) call fail_cuda_selftest('nudge qtp below level')
               end do
            end do
         end do
         deallocate(qtp_h)
      end if

      if (nsv > 0) then
         allocate(svp_h(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb:ke+khc, nsv))
         svp_h = svp_d
         do n = 1, nsv
            do k = ktest, ke
               do j = jb-jhc, je+jhc
                  do i = ib-ihc, ie+ihc
                     if (abs(svp_h(i,j,k,n) - expected) > tolerance) call fail_cuda_selftest('nudge svp')
                  end do
               end do
            end do
            do k = kb, ktest-1
               do j = jb-jhc, je+jhc
                  do i = ib-ihc, ie+ihc
                     if (abs(svp_h(i,j,k,n) - initial_p) > tolerance) call fail_cuda_selftest('nudge svp below level')
                  end do
               end do
            end do
         end do
         deallocate(svp_h)
      end if

      ! lnudge off must be a no-op.
      lnudge = .false.
      up_d = initial_p
      call nudge
      call checkCUDA(cudaDeviceSynchronize(), 'nudge disabled self-test synchronization')
      allocate(up_h(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
      up_h = up_d
      if (any(abs(up_h - initial_p) > tolerance)) call fail_cuda_selftest('nudge disabled')
      deallocate(up_h)

      lnudge    = saved_lnudge
      lnudgevel = saved_lnudgevel
      ltempeq   = saved_ltempeq
      lmoist    = saved_lmoist
      nnudge    = saved_nnudge
      tnudge    = saved_tnudge

      ! Restore the device state the rest of the run expects.
      up_d = up
      vp_d = vp
      if (ltempeq) thlp_d = thlp
      if (lmoist)  qtp_d  = qtp
      if (nsv > 0) svp_d  = svp
      u0av_d = u0av
      v0av_d = v0av
      if (ltempeq) thl0av_d = thl0av
      if (lmoist)  qt0av_d  = qt0av
      if (nsv > 0) sv0av_d  = sv0av
      uprof_d = uprof
      vprof_d = vprof
      if (ltempeq) thlprof_d = thlprof
      if (lmoist)  qtprof_d  = qtprof
      if (nsv > 0) svprof_d  = svprof
   end subroutine test_nudge_profiles

   !> Verify the IBM geometry mirrored onto the device by init_ibm_device.
   !!
   !! The mirrors are written once and then only read, so the failure this
   !! guards against is a silent one: a mask or point list that was never
   !! transferred reads as zeros on the device and simply switches off the
   !! wall-function corrections at those points instead of producing an error.
   !! Bounds are compared as well as values: the kernels index these from ib-ih
   !! and kb-kh, and a device array allocated over the wrong extent still copies
   !! and still runs, it just reads the wrong cells.
   subroutine test_ibm_device_geometry
      implicit none

      ! Bounds are inquired on the module arrays themselves. Passing them to a
      ! helper with declared bounds would make the comparison vacuous, because
      ! the dummy would impose the bounds being checked for.
      if (allocated(mask_u_d)) then
         call check_bounds('mask_u', lbound(mask_u_d), ubound(mask_u_d), lbound(mask_u), ubound(mask_u))
         call check_mask('mask_u', mask_u, mask_u_d)
      end if
      if (allocated(mask_v_d)) then
         call check_bounds('mask_v', lbound(mask_v_d), ubound(mask_v_d), lbound(mask_v), ubound(mask_v))
         call check_mask('mask_v', mask_v, mask_v_d)
      end if
      if (allocated(mask_w_d)) then
         call check_bounds('mask_w', lbound(mask_w_d), ubound(mask_w_d), lbound(mask_w), ubound(mask_w))
         call check_mask('mask_w', mask_w, mask_w_d)
      end if
      if (allocated(mask_c_d)) then
         call check_bounds('mask_c', lbound(mask_c_d), ubound(mask_c_d), lbound(mask_c), ubound(mask_c))
         call check_mask('mask_c', mask_c, mask_c_d)
      end if

      if (allocated(bndpts_u_d)) &
         call check_bndpts('bndpts_u', bound_info_u%bndpts_loc, bndpts_u_d, bound_info_u%nbndptsrank)
      if (allocated(bndpts_v_d)) &
         call check_bndpts('bndpts_v', bound_info_v%bndpts_loc, bndpts_v_d, bound_info_v%nbndptsrank)
      if (allocated(bndpts_w_d)) &
         call check_bndpts('bndpts_w', bound_info_w%bndpts_loc, bndpts_w_d, bound_info_w%nbndptsrank)
      if (allocated(bndpts_c_d)) &
         call check_bndpts('bndpts_c', bound_info_c%bndpts_loc, bndpts_c_d, bound_info_c%nbndptsrank)

   contains

      subroutine check_bounds(name, lo_d, hi_d, lo_h, hi_h)
         implicit none
         character(len=*), intent(in) :: name
         integer, intent(in) :: lo_d(3), hi_d(3), lo_h(3), hi_h(3)

         if (any(lo_d /= lo_h) .or. any(hi_d /= hi_h)) then
            call fail_cuda_selftest('ibm device geometry bounds '//name)
         end if
      end subroutine check_bounds

      subroutine check_mask(name, host, dev)
         implicit none
         character(len=*), intent(in) :: name
         real, intent(in) :: host(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
         real, device, intent(in) :: dev(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
         real, allocatable :: back(:,:,:)

         allocate(back(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         back = dev
         if (any(back /= host)) call fail_cuda_selftest('ibm device geometry '//name)
         deallocate(back)
      end subroutine check_mask

      subroutine check_bndpts(name, host, dev, npts)
         implicit none
         character(len=*), intent(in) :: name
         integer, intent(in) :: npts
         integer, intent(in) :: host(npts,3)
         integer, device, intent(in) :: dev(npts,3)
         integer, allocatable :: back(:,:)

         if (npts == 0) return
         allocate(back(npts,3))
         back = dev
         if (any(back /= host)) call fail_cuda_selftest('ibm device geometry '//name)
         deallocate(back)
      end subroutine check_bndpts

   end subroutine test_ibm_device_geometry

   !> Compare the device diff*_corr routines against the host originals.
   !!
   !! The masks are synthesised rather than taken from the case geometry. Real
   !! geometries are not neutral: on the IBM test case only the i+, i- and k-
   !! neighbour branches ever fire, so a port that mishandled the j or k+ terms
   !! would compare equal and pass. A checkerboard of solid cells makes all six
   !! branches live for a large share of the boundary points, independently of
   !! which case this runs under.
   !!
   !! Also asserts what lets the device loops omit atomics: that the entries of
   !! each bndpts list are distinct cells. If a cell were listed twice the host
   !! loop would apply both corrections and the device loop would lose one.
   !> The mirrored cell-volume reciprocals.
   !!
   !! Two things, because a mirror comparison alone would not be enough. First
   !! that each device array carries its host array's values. Second that each
   !! host array really is the reciprocal of its own cell volume - dzf and dzh
   !! are equal on the uniform vertical grid every case in the suite uses, so a
   !! swap between the two would survive both the mirror check and the whole
   !! CPU-GPU parity matrix. The definition is asserted rather than the product
   !! recomputed, so the two sides cannot share a mistake.
   subroutine test_cell_volume_reciprocals
      implicit none

      real, allocatable :: back(:)
      integer :: k

      allocate(back(kb-kh:ke+kh))
      back = dxdydzfi_d
      if (any(back /= dxdydzfi)) call fail_cuda_selftest('dxdydzfi mirror')
      do k = kb-kh, ke+kh
         if (abs(back(k)*(dx*dy*dzf(k)) - 1.) > 1.e-13) &
            call fail_cuda_selftest('dxdydzfi is not 1/(dx*dy*dzf)')
      end do
      deallocate(back)

      allocate(back(kb:ke+kh))
      back = dxdydzhi_d
      if (any(back /= dxdydzhi)) call fail_cuda_selftest('dxdydzhi mirror')
      do k = kb, ke+kh
         if (abs(back(k)*(dx*dy*dzh(k)) - 1.) > 1.e-13) &
            call fail_cuda_selftest('dxdydzhi is not 1/(dx*dy*dzh)')
      end do
      deallocate(back)
   end subroutine test_cell_volume_reciprocals

   !> The wall-function section caches against the expressions they replaced.
   !!
   !! The comparison lives in modibm, next to the arrays, because they are not
   !! exported; this drives it and reports in the usual way. It is exact: both
   !! sides evaluate the same expression on the same host, so any difference is
   !! a plumbing error - a stencil on the wrong staggered grid, a cache column
   !! read under the wrong name, an index left global.
   subroutine test_ibm_section_cache
      implicit none

      character(len=128) :: problem

      call check_ibm_section_cache(problem)
      if (problem /= '') call fail_cuda_selftest(trim(problem))
   end subroutine test_ibm_section_cache

   subroutine test_ibm_diff_corr
      implicit none

      real, allocatable :: u0_s(:,:,:), v0_s(:,:,:), w0_s(:,:,:)
      real, allocatable :: up_s(:,:,:), vp_s(:,:,:), wp_s(:,:,:)
      real, allocatable :: thl0_s(:,:,:), thlp_s(:,:,:), ekm_s(:,:,:), ekh_s(:,:,:)
      real, allocatable :: mask_u_s(:,:,:), mask_v_s(:,:,:), mask_w_s(:,:,:), mask_c_s(:,:,:)

      if (.not. allocated(mask_u_d)) return
      if (.not. allocated(bndpts_u_d)) return

      call check_distinct('bndpts_u', bound_info_u%bndpts_loc, bound_info_u%nbndptsrank)
      call check_distinct('bndpts_v', bound_info_v%bndpts_loc, bound_info_v%nbndptsrank)
      call check_distinct('bndpts_w', bound_info_w%bndpts_loc, bound_info_w%nbndptsrank)
      if (allocated(bndpts_c_d)) &
         call check_distinct('bndpts_c', bound_info_c%bndpts_loc, bound_info_c%nbndptsrank)

      ! The self-tests run before the time loop, so everything touched here has
      ! to be handed back exactly as it was found, on both host and device.
      allocate(u0_s, source=u0); allocate(v0_s, source=v0); allocate(w0_s, source=w0)
      allocate(up_s, source=up); allocate(vp_s, source=vp); allocate(wp_s, source=wp)
      allocate(ekm_s, source=ekm); allocate(ekh_s, source=ekh)
      allocate(mask_u_s, source=mask_u)
      allocate(mask_v_s, source=mask_v)
      allocate(mask_w_s, source=mask_w)

      call checkerboard(mask_u); call checkerboard(mask_v); call checkerboard(mask_w)
      mask_u_d = mask_u; mask_v_d = mask_v; mask_w_d = mask_w

      call seed(u0, 1.); call seed(v0, 2.); call seed(w0, 3.)
      call seed(up, 4.); call seed(vp, 5.); call seed(wp, 6.)
      call seed(ekm, 7.); call seed(ekh, 8.)
      ekm = 0.5 + 0.01*ekm   ! diffusivities are positive and O(1)
      ekh = 0.5 + 0.01*ekh

      u0_d = u0; v0_d = v0; w0_d = w0
      up_d = up; vp_d = vp; wp_d = wp
      ekm_d = ekm; ekh_d = ekh

      call diffu_corr
      call diffu_corr_device
      call compare('diffu_corr', up, up_d, up_s, ih, kh)

      call diffv_corr
      call diffv_corr_device
      call compare('diffv_corr', vp, vp_d, vp_s, ih, kh)

      call diffw_corr
      call diffw_corr_device
      call compare('diffw_corr', wp, wp_d, wp_s, ih, kh)

      if (allocated(bndpts_c_d) .and. ltempeq) then
         allocate(thl0_s, source=thl0); allocate(thlp_s, source=thlp)
         allocate(mask_c_s, source=mask_c)
         call checkerboard(mask_c)
         mask_c_d = mask_c
         call seed(thl0, 9.); call seed(thlp, 10.)
         thl0_d = thl0; thlp_d = thlp
         call diffc_corr(thl0, thlp, ih, jh, kh)
         call diffc_corr_device(thl0_d, thlp_d, ih, jh, kh)
         call compare('diffc_corr thl', thlp, thlp_d, thlp_s, ih, kh)
         thl0 = thl0_s; thlp = thlp_s; mask_c = mask_c_s
         thl0_d = thl0; thlp_d = thlp; mask_c_d = mask_c
         deallocate(thl0_s, thlp_s, mask_c_s)
      end if

      u0 = u0_s; v0 = v0_s; w0 = w0_s
      up = up_s; vp = vp_s; wp = wp_s
      ekm = ekm_s; ekh = ekh_s
      mask_u = mask_u_s; mask_v = mask_v_s; mask_w = mask_w_s
      u0_d = u0; v0_d = v0; w0_d = w0
      up_d = up; vp_d = vp; wp_d = wp
      ekm_d = ekm; ekh_d = ekh
      mask_u_d = mask_u; mask_v_d = mask_v; mask_w_d = mask_w
      deallocate(u0_s, v0_s, w0_s, up_s, vp_s, wp_s, ekm_s, ekh_s)
      deallocate(mask_u_s, mask_v_s, mask_w_s)

   contains

      !> Mark a repeating subset of cells solid, so that every one of the six
      !! neighbour tests in the diff*_corr routines is true somewhere.
      subroutine checkerboard(msk)
         implicit none
         real, intent(out) :: msk(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
         integer :: i, j, k

         do k = kb-kh, ke+kh
            do j = jb-jh, je+jh
               do i = ib-ih, ie+ih
                  if (modulo(i + j + k, 3) == 0) then
                     msk(i,j,k) = 0.
                  else
                     msk(i,j,k) = 1.
                  end if
               end do
            end do
         end do
      end subroutine checkerboard

      !> Fill a with a value that varies in all three directions, so that every
      !! stencil neighbour a correction reads is distinct.
      subroutine seed(a, offset)
         implicit none
         real, intent(inout) :: a(:,:,:)
         real, intent(in)    :: offset
         integer :: i, j, k

         do k = 1, size(a,3)
            do j = 1, size(a,2)
               do i = 1, size(a,1)
                  a(i,j,k) = offset + 0.125*i - 0.0625*j + 0.03125*k
               end do
            end do
         end do
      end subroutine seed

      !> Compare device against host, and refuse to pass if neither moved.
      subroutine compare(name, host, dev, before, hi, hk)
         implicit none
         character(len=*), intent(in) :: name
         integer, intent(in) :: hi, hk
         real, intent(in) :: host(ib-hi:ie+hi,jb-jh:je+jh,kb:ke+hk)
         real, device, intent(in) :: dev(ib-hi:ie+hi,jb-jh:je+jh,kb:ke+hk)
         real, intent(in) :: before(ib-hi:ie+hi,jb-jh:je+jh,kb:ke+hk)
         real, allocatable :: back(:,:,:)
         real :: worst, scale

         allocate(back(ib-hi:ie+hi,jb-jh:je+jh,kb:ke+hk))
         back = dev
         if (any(back /= back)) call fail_cuda_selftest('ibm diff_corr '//name//' device produced NaN')

         ! A routine that corrected nothing would agree with any port of itself.
         if (maxval(abs(host - before)) == 0.) then
            call fail_cuda_selftest('ibm diff_corr applied no correction: '//name)
         end if

         worst = maxval(abs(back - host))
         scale = max(1., maxval(abs(host)))
         ! Not required to be exact: the device compiler may contract these
         ! expressions into FMAs where the host compiler does not.
         if (worst > 1.e-10 * scale) then
            write(*,'(A,A,A,ES12.4)') 'ibm diff_corr mismatch ', name, ' worst ', worst
            call fail_cuda_selftest('ibm diff_corr '//name)
         end if
         deallocate(back)
      end subroutine compare

      !> Assert no cell appears twice in a boundary-point list.
      subroutine check_distinct(name, pts, npts)
         implicit none
         character(len=*), intent(in) :: name
         integer, intent(in) :: npts
         integer, intent(in) :: pts(npts,3)
         integer, allocatable :: seen(:,:,:)
         integer :: n

         if (npts <= 1) return
         allocate(seen(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         seen = 0
         do n = 1, npts
            if (seen(pts(n,1),pts(n,2),pts(n,3)) /= 0) then
               call fail_cuda_selftest('ibm boundary points not distinct: '//name)
            end if
            seen(pts(n,1),pts(n,2),pts(n,3)) = n
         end do
         deallocate(seen)
      end subroutine check_distinct

   end subroutine test_ibm_diff_corr

   !> Compare the device wallfunmom against the host original.
   !!
   !! Run for all three staggered directions on a seeded velocity and
   !! temperature field. The per-facet stress accumulator is compared as well as
   !! the tendency, because it is the one output written under a different
   !! indexing (by facet, not by cell) and so has its own way of going wrong.
   !!
   !! The comparison tolerance is looser than for diff_corr: the device sums
   !! into rhs and fac_tau with atomics, so contributions from sections sharing
   !! a cell or a facet arrive in an order that varies between runs.
   subroutine test_ibm_wallfunmom
      implicit none

      real, allocatable :: u0_s(:,:,:), v0_s(:,:,:), w0_s(:,:,:), thl0_s(:,:,:)
      real, allocatable :: up_s(:,:,:), vp_s(:,:,:), wp_s(:,:,:)

      if (.not. allocated(fac_tau_d)) return
      if (iwallmom <= 1) return

      allocate(u0_s, source=u0); allocate(v0_s, source=v0); allocate(w0_s, source=w0)
      allocate(thl0_s, source=thl0)
      allocate(up_s, source=up); allocate(vp_s, source=vp); allocate(wp_s, source=wp)

      call seed_wf(u0, 1.); call seed_wf(v0, 2.); call seed_wf(w0, 3.)
      thl0 = 290. + 0.5*thl0_s*0.  ! uniform enough to keep the stability branch finite
      call seed_wf(thl0, 290.)
      up = 0.; vp = 0.; wp = 0.

      u0_d = u0; v0_d = v0; w0_d = w0; thl0_d = thl0
      up_d = up; vp_d = vp; wp_d = wp

      call one_direction(1, xhat, bound_info_u, up, up_d)
      call one_direction(2, yhat, bound_info_v, vp, vp_d)
      call one_direction(3, zhat, bound_info_w, wp, wp_d)

      u0 = u0_s; v0 = v0_s; w0 = w0_s; thl0 = thl0_s
      up = up_s; vp = vp_s; wp = wp_s
      u0_d = u0; v0_d = v0; w0_d = w0; thl0_d = thl0
      up_d = up; vp_d = vp; wp_d = wp
      deallocate(u0_s, v0_s, w0_s, thl0_s, up_s, vp_s, wp_s)

   contains

      subroutine one_direction(dirsel, dirvec, bi, rhs, rhs_d)
         implicit none
         integer, intent(in) :: dirsel
         real, intent(in) :: dirvec(3)
         type(bound_info_type), intent(inout) :: bi
         real, intent(inout) :: rhs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
         real, device, intent(inout) :: rhs_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)

         real, allocatable :: back(:,:,:), host_tau(:), dev_tau(:)
         real :: worst, scale
         character(len=1) :: tag

         if (bi%nfctsecsrank < 1) return
         tag = 'u'
         if (dirsel == 2) tag = 'v'
         if (dirsel == 3) tag = 'w'

         rhs = 0.
         rhs_d = rhs
         call wallfunmom(dirvec, rhs, bi)
         call wallfunmom_dir_device(dirsel, rhs_d)
         call checkCUDA(cudaDeviceSynchronize(), 'wallfunmom self-test synchronization')

         ! A direction that produced no stress anywhere would agree with any
         ! port of itself, so refuse to pass on that.
         if (maxval(abs(rhs)) == 0.) then
            call fail_cuda_selftest('ibm wallfunmom produced no tendency: '//tag)
         end if

         allocate(back(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
         back = rhs_d
         if (any(back /= back)) call fail_cuda_selftest('ibm wallfunmom '//tag//' device produced NaN')
         if (any(rhs /= rhs)) call fail_cuda_selftest('ibm wallfunmom '//tag//' host produced NaN')
         worst = maxval(abs(back - rhs))
         scale = max(1., maxval(abs(rhs)))
         if (worst > 1.e-8 * scale) then
            write(*,'(A,A,A,ES12.4,A,ES12.4)') 'ibm wallfunmom mismatch ', tag, ' worst ', worst, &
                 ' scale ', scale
            call fail_cuda_selftest('ibm wallfunmom '//tag)
         end if
         deallocate(back)

         allocate(host_tau(nfcts), dev_tau(nfcts))
         host_tau = fac_tau_raw
         dev_tau  = fac_tau_d
         worst = maxval(abs(dev_tau - host_tau))
         scale = max(1., maxval(abs(host_tau)))
         if (worst > 1.e-8 * scale) then
            write(*,'(A,A,A,ES12.4)') 'ibm wallfunmom fac_tau mismatch ', tag, ' worst ', worst
            call fail_cuda_selftest('ibm wallfunmom fac_tau '//tag)
         end if
         deallocate(host_tau, dev_tau)
      end subroutine one_direction

      subroutine seed_wf(a, offset)
         implicit none
         real, intent(inout) :: a(:,:,:)
         real, intent(in)    :: offset
         integer :: i, j, k

         do k = 1, size(a,3)
            do j = 1, size(a,2)
               do i = 1, size(a,1)
                  a(i,j,k) = offset + 0.03125*i - 0.015625*j + 0.0078125*k
               end do
            end do
         end do
      end subroutine seed_wf

   end subroutine test_ibm_wallfunmom

   !> Compare the device wallfunheat against the host original.
   !!
   !! Checks both tendencies, the two whole-domain flux sums, and all four
   !! per-facet accumulators. The pressure pair is worth its own comparison
   !! because it is accumulated before the lskipsec test, so it covers sections
   !! the rest of the routine never reaches.
   subroutine test_ibm_wallfunheat
      implicit none

      real, allocatable :: u0_s(:,:,:), v0_s(:,:,:), w0_s(:,:,:)
      real, allocatable :: thl0_s(:,:,:), qt0_s(:,:,:), pres0_s(:,:,:)
      real, allocatable :: thlp_s(:,:,:), qtp_s(:,:,:)
      real, allocatable :: back(:,:,:)
      real, allocatable :: h_htc(:), h_cth(:), h_pres(:), h_pres2(:)
      real :: host_totheat, host_totq, dev_totheat, dev_totq
      real :: save_totheat, save_totq
      real, allocatable :: fachf_s(:), facef_s(:)
      logical, allocatable :: faclGR_s(:)
      real, allocatable :: facqsat_s(:), fachurel_s(:)
      integer :: saved_iwallmoist

      if (.not. allocated(fac_htc_d)) return
      if (bound_info_c%nfctsecsrank < 1) return

      allocate(u0_s, source=u0); allocate(v0_s, source=v0); allocate(w0_s, source=w0)
      allocate(thl0_s, source=thl0); allocate(pres0_s, source=pres0)
      allocate(thlp_s, source=thlp)
      if (lmoist) then
         allocate(qt0_s, source=qt0); allocate(qtp_s, source=qtp)
      end if
      save_totheat = totheatflux
      save_totq    = totqflux
      ! Both wall functions add into the energy balance accumulators when lEB is
      ! set, on the host and on the device. Neither is reset by anything before
      ! the run proper, so the test has to put both back itself: the device side
      ! now carries its total across the Runge-Kutta stages and would otherwise
      ! deliver this test's contribution into the first real time step.
      if (allocated(fachf)) allocate(fachf_s, source=fachf)
      if (allocated(facef)) allocate(facef_s, source=facef)

      ! No case in the suite has a green-roof facet, so the moisture wall
      ! function would otherwise compare zero against zero. Turning every facet
      ! into a green roof, on both sides, makes the branch actually run.
      saved_iwallmoist = iwallmoist
      if (lmoist .and. allocated(faclGR) .and. allocated(facqsat)) then
         allocate(faclGR_s, source=faclGR)
         allocate(facqsat_s, source=facqsat)
         allocate(fachurel_s, source=fachurel)
         iwallmoist = 2
         faclGR   = .true.
         facqsat  = 0.02
         fachurel = 1.0
         if (allocated(faclGR_d))   faclGR_d   = faclGR
         if (allocated(facqsat_d))  facqsat_d  = facqsat
         if (allocated(fachurel_d)) fachurel_d = fachurel
         if (allocated(facf_d))     facf_d     = facf
      end if

      call seed_wh(u0, 1.); call seed_wh(v0, 2.); call seed_wh(w0, 3.)
      call seed_wh(thl0, 290.); call seed_wh(pres0, 5.)
      thlp = 0.
      if (lmoist) then
         call seed_wh(qt0, 0.01)
         qtp = 0.
      end if

      u0_d = u0; v0_d = v0; w0_d = w0
      thl0_d = thl0; pres0_d = pres0; thlp_d = thlp
      if (lmoist) then
         qt0_d = qt0; qtp_d = qtp
      end if

      totheatflux = 0.
      totqflux    = 0.
      call wallfunheat
      host_totheat = totheatflux
      host_totq    = totqflux
      allocate(h_htc, source=fac_htc_raw)
      allocate(h_cth, source=fac_cth_raw)
      allocate(h_pres, source=fac_pres_raw)
      allocate(h_pres2, source=fac_pres2_raw)

      totheatflux = 0.
      totqflux    = 0.
      call wallfunheat_dir_device
      call checkCUDA(cudaDeviceSynchronize(), 'wallfunheat self-test synchronization')
      dev_totheat = totheatflux
      dev_totq    = totqflux

      ! The pressure accumulator is non-zero for any geometry at all, so it is
      ! the one that proves the loop ran.
      if (maxval(abs(h_pres)) == 0.) call fail_cuda_selftest('ibm wallfunheat saw no sections')

      allocate(back(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
      back = thlp_d
      call cmp3('wallfunheat thlp', back, thlp)
      if (lmoist) then
         ! A moisture branch that produced nothing would agree with any port of
         ! itself, so refuse to pass on that.
         if (allocated(faclGR_s) .and. maxval(abs(qtp)) == 0.) then
            call fail_cuda_selftest('ibm wallfunheat produced no moisture tendency')
         end if
         back = qtp_d
         call cmp3('wallfunheat qtp', back, qtp)
      end if
      deallocate(back)

      call cmp1('wallfunheat fac_pres',  fac_pres_d,  h_pres)
      call cmp1('wallfunheat fac_pres2', fac_pres2_d, h_pres2)
      call cmp1('wallfunheat fac_htc',   fac_htc_d,   h_htc)
      call cmp1('wallfunheat fac_cth',   fac_cth_d,   h_cth)
      call cmp0('wallfunheat totheatflux', dev_totheat, host_totheat)
      call cmp0('wallfunheat totqflux',    dev_totq,    host_totq)

      u0 = u0_s; v0 = v0_s; w0 = w0_s; thl0 = thl0_s; pres0 = pres0_s; thlp = thlp_s
      u0_d = u0; v0_d = v0; w0_d = w0
      thl0_d = thl0; pres0_d = pres0; thlp_d = thlp
      if (lmoist) then
         qt0 = qt0_s; qtp = qtp_s
         qt0_d = qt0; qtp_d = qtp
         deallocate(qt0_s, qtp_s)
      end if
      totheatflux = save_totheat
      totqflux    = save_totq
      if (allocated(fachf_s)) then
         fachf = fachf_s
         deallocate(fachf_s)
      end if
      if (allocated(facef_s)) then
         facef = facef_s
         deallocate(facef_s)
      end if
      if (allocated(fachf_d)) fachf_d = 0.
      if (allocated(facef_d)) facef_d = 0.
      iwallmoist = saved_iwallmoist
      if (allocated(faclGR_s)) then
         faclGR = faclGR_s; facqsat = facqsat_s; fachurel = fachurel_s
         if (allocated(faclGR_d))   faclGR_d   = faclGR
         if (allocated(facqsat_d))  facqsat_d  = facqsat
         if (allocated(fachurel_d)) fachurel_d = fachurel
         deallocate(faclGR_s, facqsat_s, fachurel_s)
      end if
      deallocate(u0_s, v0_s, w0_s, thl0_s, pres0_s, thlp_s)
      deallocate(h_htc, h_cth, h_pres, h_pres2)

   contains

      subroutine seed_wh(a, offset)
         implicit none
         real, intent(inout) :: a(:,:,:)
         real, intent(in)    :: offset
         integer :: i, j, k

         do k = 1, size(a,3)
            do j = 1, size(a,2)
               do i = 1, size(a,1)
                  a(i,j,k) = offset + 0.03125*i - 0.015625*j + 0.0078125*k
               end do
            end do
         end do
      end subroutine seed_wh

      subroutine cmp3(name, dev, host)
         implicit none
         character(len=*), intent(in) :: name
         real, intent(in) :: dev(:,:,:), host(:,:,:)
         real :: worst, scale

         if (any(dev /= dev)) call fail_cuda_selftest('ibm '//name//' device produced NaN')
         if (any(host /= host)) call fail_cuda_selftest('ibm '//name//' host produced NaN')
         worst = maxval(abs(dev - host))
         scale = max(1., maxval(abs(host)))
         if (worst > 1.e-8 * scale) then
            write(*,'(A,A,A,ES12.4)') 'ibm ', name, ' mismatch worst ', worst
            call fail_cuda_selftest('ibm '//name)
         end if
      end subroutine cmp3

      subroutine cmp1(name, dev, host)
         implicit none
         character(len=*), intent(in) :: name
         real, device, intent(in) :: dev(:)
         real, intent(in) :: host(:)
         real, allocatable :: b(:)
         real :: worst, scale

         allocate(b(size(host)))
         b = dev
         if (any(b /= b)) call fail_cuda_selftest('ibm '//name//' device produced NaN')
         if (any(host /= host)) call fail_cuda_selftest('ibm '//name//' host produced NaN')
         worst = maxval(abs(b - host))
         scale = max(1., maxval(abs(host)))
         if (worst > 1.e-8 * scale) then
            write(*,'(A,A,A,ES12.4)') 'ibm ', name, ' mismatch worst ', worst
            call fail_cuda_selftest('ibm '//name)
         end if
         deallocate(b)
      end subroutine cmp1

      subroutine cmp0(name, dev, host)
         implicit none
         character(len=*), intent(in) :: name
         real, intent(in) :: dev, host

         if (dev /= dev .or. host /= host) call fail_cuda_selftest('ibm '//name//' is NaN')
         if (abs(dev - host) > 1.e-8 * max(1., abs(host))) then
            write(*,'(A,A,A,ES14.6,A,ES14.6)') 'ibm ', name, ' mismatch dev ', dev, ' host ', host
            call fail_cuda_selftest('ibm '//name)
         end if
      end subroutine cmp0

   end subroutine test_ibm_wallfunheat

   !> Verify how the energy balance accumulators reach the host.
   !!
   !! The handover must happen on EVERY Runge-Kutta stage, and must clear the
   !! device copy as it goes. It is tempting to accumulate on the device and
   !! deliver once per step, but intqH in modEB resets fachf and facef at the
   !! end of every stage, outside its own rk3step == 3 test, so the host array
   !! holds one stage when the reduction reads it. Delivering three stages at
   !! once makes the facet heat flux three times too large.
   !!
   !! That is not a hypothetical: it shipped, and a manual CPU/GPU comparison
   !! found it. This test and the facEB comparison in the surface-energy-balance
   !! parity case are the two things that now stand between it and a rerun.
   subroutine test_facflux_handover
      implicit none

      real, parameter :: probe = 3.25
      real, allocatable :: fachf_s(:), facef_s(:), back(:)
      integer :: saved_rk3step, stage

      if (.not. lEB) return
      if (.not. allocated(fachf_d)) return
      if (.not. allocated(fachf)) return

      allocate(fachf_s, source=fachf)
      if (allocated(facef)) allocate(facef_s, source=facef)
      saved_rk3step = rk3step
      allocate(back(0:nfcts))

      do stage = 1, 3
         rk3step = stage

         ! One stage's contribution must arrive on the host, whichever stage.
         fachf = 0.
         fachf_d = probe
         call updateFacFluxHost
         if (maxval(abs(fachf(0:nfcts) - probe)) > 64.*epsilon(1.)*probe) then
            call fail_cuda_selftest('facflux not handed over on every stage')
         end if

         ! And the device copy must be cleared, so the next stage starts from
         ! zero instead of carrying this stage forward.
         back = fachf_d
         if (any(back /= 0.)) call fail_cuda_selftest('facflux device not reset on handover')

         ! With nothing new on the device, a further handover adds nothing.
         call updateFacFluxHost
         if (maxval(abs(fachf(0:nfcts) - probe)) > 64.*epsilon(1.)*probe) then
            call fail_cuda_selftest('facflux counted twice')
         end if
      end do

      deallocate(back)
      rk3step = saved_rk3step
      fachf = fachf_s
      deallocate(fachf_s)
      if (allocated(facef_s)) then
         facef = facef_s
         deallocate(facef_s)
      end if
      fachf_d = 0.
      if (allocated(facef_d)) facef_d = 0.

   end subroutine test_facflux_handover

   !> Verify that the facet properties on the device follow the host.
   !!
   !! facT, facqsat, fachurel and facf are all rewritten by the energy balance
   !! during the run. Mirroring them once at initialisation leaves the wall
   !! functions running the whole simulation on initial-condition values, which
   !! no single-call comparison can detect: the device and host agree perfectly
   !! on any one call and disagree only because a refresh never happened.
   !!
   !! So this perturbs the host arrays and checks the refresh picks the change
   !! up, which is the property that actually matters.
   subroutine test_facet_props_refresh
      implicit none

      real, allocatable :: facT_s(:,:), facqsat_s(:), fachurel_s(:), facf_s(:,:)
      real, allocatable :: back(:)
      integer :: saved_rk3step

      if (.not. lEB) return
      if (.not. allocated(facT1_d)) return

      allocate(facT_s, source=facT)
      allocate(back(0:nfcts))
      saved_rk3step = rk3step
      rk3step = 1

      ! Move the host value, then ask for the refresh the time loop would do.
      facT(0:nfcts,1) = facT(0:nfcts,1) + 7.5
      call updateFacetPropsDevice
      back = facT1_d
      if (maxval(abs(back - facT(0:nfcts,1))) > 64.*epsilon(1.)*maxval(abs(facT(0:nfcts,1)))) then
         call fail_cuda_selftest('facet property facT is stale on the device')
      end if

      if (allocated(facqsat_d) .and. allocated(facqsat)) then
         allocate(facqsat_s, source=facqsat)
         allocate(fachurel_s, source=fachurel)
         allocate(facf_s, source=facf)

         facqsat(0:nfcts)  = facqsat(0:nfcts) + 0.125
         fachurel(0:nfcts) = fachurel(0:nfcts) + 0.0625
         facf(0:nfcts,4)   = facf(0:nfcts,4) + 11.
         facf(0:nfcts,5)   = facf(0:nfcts,5) + 13.
         call updateFacetPropsDevice

         back = facqsat_d
         if (maxval(abs(back - facqsat(0:nfcts))) > 1.e-10) &
            call fail_cuda_selftest('facet property facqsat is stale on the device')
         back = fachurel_d
         if (maxval(abs(back - fachurel(0:nfcts))) > 1.e-10) &
            call fail_cuda_selftest('facet property fachurel is stale on the device')
         call check_facf_column(4)
         call check_facf_column(5)

         facqsat = facqsat_s; fachurel = fachurel_s; facf = facf_s
         deallocate(facqsat_s, fachurel_s, facf_s)
      end if

      facT = facT_s
      deallocate(facT_s)
      rk3step = 1
      call updateFacetPropsDevice
      rk3step = saved_rk3step
      deallocate(back)

   contains

      subroutine check_facf_column(col)
         implicit none
         integer, intent(in) :: col
         real, allocatable :: b2(:,:)
         allocate(b2(0:nfcts,5))
         b2 = facf_d
         if (maxval(abs(b2(0:nfcts,col) - facf(0:nfcts,col))) > 1.e-10) then
            call fail_cuda_selftest('facet property facf is stale on the device')
         end if
         deallocate(b2)
      end subroutine check_facf_column

   end subroutine test_facet_props_refresh

   !> Check the OpenACC periodicEBcorr kernels against the physics they encode.
   !!
   !! The host branch is covered by tests.f90::tests_periodic_ebcorr, which runs
   !! the same closure argument on a CPU build. Repeating it here rather than
   !! diffing the two branches keeps the check independent of the host: a slip
   !! shared by both sides would pass a symmetric comparison.
   !!
   !! The three mistakes this port invites are covered explicitly. The host
   !! branch loops over the interior only, so a device loop that spanned the
   !! halos - the opposite of what nudge needs - is caught by the halo check.
   !! The tendency is applied by two separate kernels, so dropping the top-level
   !! one, or aiming it at the wrong level, is caught by the closure and the
   !! split. And both kernels accumulate, so an assignment is caught by starting
   !! from a non-zero baseline.
   subroutine test_periodic_ebcorr
      implicit none

      real, parameter :: initial_p = 2.

      logical :: saved_lperiodicEBcorr, uniform
      integer :: saved_sinkbase, i, j, k
      real    :: saved_fraction, saved_totheatflux, saved_totqflux
      real    :: unit_flux, mean_flux, tol, colsum, want
      real, allocatable :: thlp_h(:,:,:), qtp_h(:,:,:)

      if (.not. ltempeq) return
      if (.not. allocated(thlp_d)) return
      if (kb /= 1) return

      ! Two levels above sinkbase to separate the top-level term from the
      ! volume sink, one below it for the support check.
      if (ke - kb < 3) return

      ! The ke/M scaling weights levels by count rather than by depth, so the
      ! column integral is the domain-mean flux only on a uniform grid. The
      ! support and halo checks below do not depend on that, so they run either
      ! way.
      uniform = maxval(abs(dzf(kb:ke) - dzf(kb))) <= 1.e-12 * dzf(kb) .and. &
                abs(zh(ke+1) - real(ke - kb + 1) * dzf(kb)) <= 1.e-12 * zh(ke+1)

      saved_lperiodicEBcorr = lperiodicEBcorr
      saved_sinkbase        = sinkbase
      saved_fraction        = fraction
      saved_totheatflux     = totheatflux
      saved_totqflux        = totqflux

      allocate(thlp_h(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
      if (lmoist) allocate(qtp_h(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))

      lperiodicEBcorr = .true.
      sinkbase        = ke/2
      fraction        = 0.5

      ! Rank-dependent, so a correction built from the local flux instead of
      ! the reduced one lands on a different answer on every rank but the
      ! first. The ranks contribute 1, 2, ... nprocs units.
      unit_flux = xlen * ylen
      mean_flux = 0.5 * real(nprocs) * real(nprocs + 1)
      tol       = 1.e-12 * max(1., abs(mean_flux))

      call run_correction(1.)

      if (uniform) then
         ! The column-integrated tendency is the domain-mean flux: what the
         ! facets put in is what the correction takes out.
         colsum = 0.
         do k = kb, ke
            colsum = colsum + (thlp_h(ib,jb,k) - initial_p) * dzf(k)
         end do
         if (abs(colsum - mean_flux) > tol) call fail_cuda_selftest('periodicEBcorr column closure')

         ! Of that, (1-fraction) leaves through the top level. Both ke and ke-1
         ! carry the volume sink, so their difference isolates the top-level
         ! term without restating how the kernel forms it.
         want = (1. - fraction) * mean_flux
         if (abs((thlp_h(ib,jb,ke) - thlp_h(ib,jb,ke-1)) * dzf(ke) - want) > tol) &
            call fail_cuda_selftest('periodicEBcorr top-level share')
      end if

      ! Nothing at or below sinkbase moves, anywhere in the plane.
      do k = kb, sinkbase
         do j = jb-jh, je+jh
            do i = ib-ih, ie+ih
               if (abs(thlp_h(i,j,k) - initial_p) > tol) &
                  call fail_cuda_selftest('periodicEBcorr below sinkbase')
            end do
         end do
      end do

      ! Nothing above ke moves either: the kernels stop at the top level.
      do k = ke+1, ke+kh
         do j = jb-jh, je+jh
            do i = ib-ih, ie+ih
               if (abs(thlp_h(i,j,k) - initial_p) > tol) &
                  call fail_cuda_selftest('periodicEBcorr above ke')
            end do
         end do
      end do

      ! Above sinkbase: the interior is updated and horizontally uniform, the
      ! halo is untouched.
      do k = sinkbase+1, ke
         do j = jb-jh, je+jh
            do i = ib-ih, ie+ih
               if (i < ib .or. i > ie .or. j < jb .or. j > je) then
                  if (abs(thlp_h(i,j,k) - initial_p) > tol) &
                     call fail_cuda_selftest('periodicEBcorr halo')
               else
                  if (abs(thlp_h(i,j,k) - initial_p) <= tol) &
                     call fail_cuda_selftest('periodicEBcorr level not updated')
                  if (abs(thlp_h(i,j,k) - thlp_h(ib,jb,k)) > tol) &
                     call fail_cuda_selftest('periodicEBcorr horizontal uniformity')
               end if
            end do
         end do
      end do

      ! Heat and moisture run the same algebra on equal fluxes, so a copy-paste
      ! slip between the two kernel pairs shows up as a difference.
      if (lmoist) then
         if (maxval(abs(thlp_h - qtp_h)) > tol) &
            call fail_cuda_selftest('periodicEBcorr heat and moisture differ')
      end if

      ! Linear in the flux: doubling what came in doubles what goes out.
      if (uniform) then
         call run_correction(2.)
         colsum = 0.
         do k = kb, ke
            colsum = colsum + (thlp_h(ib,jb,k) - initial_p) * dzf(k)
         end do
         if (abs(colsum - 2. * mean_flux) > tol) &
            call fail_cuda_selftest('periodicEBcorr flux linearity')
      end if

      ! No flux, no tendency: nothing is added on top of the scaling.
      call run_correction(0.)
      if (maxval(abs(thlp_h - initial_p)) > tol) &
         call fail_cuda_selftest('periodicEBcorr zero flux')

      ! Switched off, the routine must not touch the device fields at all.
      lperiodicEBcorr = .false.
      call run_correction(1.)
      if (maxval(abs(thlp_h - initial_p)) > tol) &
         call fail_cuda_selftest('periodicEBcorr disabled')
      if (lmoist) then
         if (maxval(abs(qtp_h - initial_p)) > tol) &
            call fail_cuda_selftest('periodicEBcorr disabled moisture')
      end if

      deallocate(thlp_h)
      if (allocated(qtp_h)) deallocate(qtp_h)

      lperiodicEBcorr = saved_lperiodicEBcorr
      sinkbase        = saved_sinkbase
      fraction        = saved_fraction
      totheatflux     = saved_totheatflux
      totqflux        = saved_totqflux

      ! Restore the device state the rest of the run expects.
      thlp_d = thlp
      if (lmoist) qtp_d = qtp

   contains

      !> Reset the device tendencies, set the per-rank flux and run the
      !! correction, leaving the result in thlp_h and qtp_h.
      subroutine run_correction(amplitude)
         real, intent(in) :: amplitude

         thlp_d = initial_p
         if (lmoist) qtp_d = initial_p
         totheatflux = amplitude * unit_flux * real(myid + 1)
         totqflux    = amplitude * unit_flux * real(myid + 1)
         call periodicEBcorr
         call checkCUDA(cudaDeviceSynchronize(), 'periodicEBcorr self-test synchronization')
         thlp_h = thlp_d
         if (lmoist) qtp_h = qtp_d

      end subroutine run_correction

   end subroutine test_periodic_ebcorr

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
