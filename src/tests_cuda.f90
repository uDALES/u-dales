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
                        dxdydzfi_d, dxdydzhi_d, &
                        col_d, col_stage, IIu_d, IIv_d
   use modfields, only : u0, v0, w0, um, vm, wm, e120, e12m, pres0, &
                         thl0, thlm, thl0c, qt0, qtm, sv0, svm, thlp, wp, &
                         up, vp, qtp, svp, &
                         u0av, v0av, thl0av, qt0av, sv0av, &
                         uprof, vprof, thlprof, qtprof, svprof, &
                         IIc, IIu, IIv, IIcs, IIus, IIvs, &
                         uoutarea, voutarea, udef, vdef, uouttot
   use modglobal, only : ib, ie, jb, je, kb, ke, ih, jh, kh, nsv, &
                         ihc, jhc, khc, dxhci, dxfc, dxfci, dxi, dyi, &
                         dzhci, dzfc, dzfci, dzfi, eps1, &
                         dx, dy, dxf, dzf, dzh, zh, rslabs, dxdydzfi, dxdydzhi, &
                         lheatpump, lfan_hp, nhppoints, ltempeq, &
                         lmoist, lnudge, lnudgevel, tnudge, nnudge, &
                         iwallmom, nfcts, xhat, yhat, zhat, &
                         totheatflux, totqflux, lEB, rk3step, iwallmoist, &
                         lperiodicEBcorr, sinkbase, fraction, xlen, ylen, &
                         linoutflow, luoutflowr, lvoutflowr, luvolflowr, lvvolflowr, &
                         uflowrate, vflowrate, rk3coef, rk3coefi
   use modforces,   only : nudge, periodicEBcorr, masscorr, calcfluidvolumes
   use modheatpump, only : heatpump, nhppoints_local, idhppts_local_d, &
                           thl_dot_hp, w_hp_exhaust
   use initfac,     only : fachf, facef, facT, facqsat, fachurel, facf, faclGR, lfacetprops_dirty
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
                           check_ibm_section_cache, &
                           solid, advecc2nd_corr_conservative, advecc2nd_corr_liberal, &
                           solid_info_type, solid_info_u, solid_info_v, solid_info_w, solid_info_c, &
                           solpts_u_d, solpts_v_d, solpts_w_d, solpts_c_d, &
                           solid_device, solid_masked_device, &
                           advecc2nd_corr_conservative_device, advecc2nd_corr_liberal_device
   use modsubgriddata, only : ekm, ekh
   use modinletdata, only : u0driver
   use modmpi,   only : myid, nprocs, nprocy, comm3d, mpierr, MY_REAL, MPI_SUM, MPI_INTEGER
   use decomp_2d, only : zstart

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
         call test_masscorr
         call test_ibm_device_geometry
         call test_cell_volume_reciprocals
         call test_ibm_section_cache
         call test_ibm_diff_corr
         call test_ibmnorm
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
   !! There are two properties here and keeping them apart is the whole point,
   !! because conflating them is what made this handover go wrong once already.
   !!
   !! The device accumulators must be cleared on EVERY Runge-Kutta stage. The
   !! wall functions add into them each stage; let them run on and the third
   !! stage delivers the sum of all three, which is the facet heat flux three
   !! times too large. That bug shipped, and a manual CPU/GPU comparison found
   !! it. This test and the facEB comparison in the surface-energy-balance
   !! parity case are what stand between it and a rerun.
   !!
   !! The copy belongs on the third stage alone. intqH in modEB reduces and
   !! time-integrates there, and clears its host arrays at the end of every
   !! stage regardless of which one it was - so whatever the first two stages
   !! hand over is discarded untouched, and copying it is two crossings in
   !! three spent on a value nothing reads. That claim about intqH is not
   !! taken on trust here: tests.f90::tests_eb checks it against the routine
   !! itself on a CPU build.
   !!
   !! The last block is the one that ties the two together. It plays three
   !! stages through the way a time step does - the wall functions add their
   !! share to whatever the device already holds, the handover follows, intqH
   !! clears the host - and asserts that what the third stage leaves for the
   !! reduction is one stage's worth. Drop the clear from the discarded stages
   !! and it arrives at three times that.
   subroutine test_facflux_handover
      implicit none

      real, parameter :: probe = 3.25
      real, allocatable :: fachf_s(:), facef_s(:), back(:)
      integer :: saved_rk3step, stage
      real    :: consumed

      if (.not. lEB) return
      if (.not. allocated(fachf_d)) return
      if (.not. allocated(fachf)) return

      allocate(fachf_s, source=fachf)
      if (allocated(facef)) allocate(facef_s, source=facef)
      saved_rk3step = rk3step
      allocate(back(0:nfcts))

      do stage = 1, 3
         rk3step = stage

         fachf = 0.
         fachf_d = probe
         call updateFacFluxHost

         if (stage == 3) then
            ! The consumed stage has to arrive.
            if (maxval(abs(fachf(0:nfcts) - probe)) > 64.*epsilon(1.)*probe) then
               call fail_cuda_selftest('facflux not handed over on the consumed stage')
            end if
         else
            ! The discarded stages must not cross the bus at all.
            if (maxval(abs(fachf(0:nfcts))) /= 0.) then
               call fail_cuda_selftest('facflux crossed on a stage intqH discards')
            end if
         end if

         ! And the device copy must be cleared on every stage, so the next one
         ! starts from zero instead of carrying this one forward.
         back = fachf_d
         if (any(back /= 0.)) call fail_cuda_selftest('facflux device not reset on handover')

         ! With nothing new on the device, a further handover adds nothing.
         call updateFacFluxHost
         if (stage == 3) then
            if (maxval(abs(fachf(0:nfcts) - probe)) > 64.*epsilon(1.)*probe) then
               call fail_cuda_selftest('facflux counted twice')
            end if
         end if
      end do

      ! One time step, played through.
      fachf   = 0.
      fachf_d = 0.
      consumed = -1.
      do stage = 1, 3
         rk3step = stage
         back = fachf_d               ! whatever the handover left there
         back = back + probe          ! what the wall functions add this stage
         fachf_d = back
         call updateFacFluxHost
         if (stage == 3) consumed = maxval(abs(fachf(0:nfcts)))
         fachf = 0.                   ! intqH clears the host at every stage
      end do
      if (abs(consumed - probe) > 64.*epsilon(1.)*probe) then
         call fail_cuda_selftest('the reduction would read the wrong number of stages')
      end if

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
   !!
   !! The refresh is now gated on lfacetprops_dirty, which modEB sets on the
   !! steps where the energy balance actually rewrites these four arrays - one
   !! step in dtEB/dt, rather than every step as before. That makes two
   !! properties to check, not one: the refresh must follow the host when the
   !! flag is set, and it must genuinely skip when it is not. A gate that never
   !! skips saves nothing; a gate that never follows leaves the wall functions
   !! on initial-condition temperatures for the whole run. The half of the
   !! contract that lives in modEB - that the flag is set exactly when those
   !! arrays move - is checked by tests.f90::tests_eb on a CPU build.
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

      ! Move the host value, then ask for the refresh the time loop would do.
      facT(0:nfcts,1) = facT(0:nfcts,1) + 7.5
      lfacetprops_dirty = .true.
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
         lfacetprops_dirty = .true.
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

      ! The gate is real: a change made without the flag must not cross.
      facT(0:nfcts,1) = facT(0:nfcts,1) + 3.25
      lfacetprops_dirty = .false.
      call updateFacetPropsDevice
      back = facT1_d
      if (maxval(abs(back - facT(0:nfcts,1))) < 1.) then
         call fail_cuda_selftest('the facet property refresh ignores its dirty flag')
      end if

      facT = facT_s
      deallocate(facT_s)
      lfacetprops_dirty = .true.
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

   !> Check the OpenACC flow-rate correction against the rate it must produce.
   !!
   !! The host branch is covered by tests.f90::tests_masscorr on a CPU build.
   !! Repeating the argument here rather than diffing the two branches keeps
   !! the check independent: a mistake shared by both sides passes a symmetric
   !! comparison.
   !!
   !! Two phases, as on the host. The first uses whatever masks the case has,
   !! so a device reduction that dropped the fluid-cell mask is caught wherever
   !! the case has solid cells. The second clears the masks, which is the only
   !! configuration in which the outlet-plane branches can hit their target -
   !! the area they divide by comes from IIc while the branch masks with IIu or
   !! IIv, and those agree only when the outlet is clear.
   !!
   !! Every case gets all four branches, whatever its namelist asks for: the
   !! test allocates the mask mirrors and the staging column that initCUDA
   !! leaves out when a controller is off, and gives back exactly what it took.
   subroutine test_masscorr
      implicit none

      real, parameter :: rk3_test = 3., target_u = 1.25, target_v = -0.75

      logical :: saved_linoutflow, saved_uout, saved_vout, saved_uvol, saved_vvol
      logical :: mine_col, mine_IIu, mine_IIv
      real    :: saved_rk3coef, saved_rk3coefi, saved_uflowrate, saved_vflowrate
      real    :: saved_udef, saved_vdef, saved_uouttot
      real    :: got
      integer :: k
      real, allocatable :: um_h(:,:,:), vm_h(:,:,:)
      real, allocatable :: up_h(:,:,:), vp_h(:,:,:)
      real, allocatable :: up_ref(:,:,:), vp_ref(:,:,:)
      integer, allocatable :: IIc_save(:,:,:), IIu_save(:,:,:), IIv_save(:,:,:)
      integer, allocatable :: IIcs_save(:), IIus_save(:), IIvs_save(:)

      if (.not. allocated(up_d)) return

      ! initCUDA allocates these only for the controllers the namelist turns
      ! on. Borrow or create whatever is missing so every branch is reachable,
      ! and remember which ones to hand back.
      mine_col = .not. allocated(col_d)
      if (mine_col) then
         allocate(col_d(kb:ke+kh))
         allocate(col_stage(kb:ke+kh))
      end if
      mine_IIu = .not. allocated(IIu_d)
      if (mine_IIu) then
         allocate(IIu_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
         IIu_d = IIu
      end if
      mine_IIv = .not. allocated(IIv_d)
      if (mine_IIv) then
         allocate(IIv_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
         IIv_d = IIv
      end if

      saved_linoutflow = linoutflow
      saved_uout       = luoutflowr
      saved_vout       = lvoutflowr
      saved_uvol       = luvolflowr
      saved_vvol       = lvvolflowr
      saved_rk3coef    = rk3coef
      saved_rk3coefi   = rk3coefi
      saved_uflowrate  = uflowrate
      saved_vflowrate  = vflowrate
      ! masscorr leaves these behind, and uouttot is read by the outflow
      ! boundary conditions during the first step - boundary only recomputes
      ! it at the end of that step, after the damage would be done.
      saved_udef       = udef
      saved_vdef       = vdef
      saved_uouttot    = uouttot

      linoutflow = .false.
      rk3coef    = rk3_test
      rk3coefi   = 1. / rk3_test
      uflowrate  = target_u
      vflowrate  = target_v

      allocate(um_h(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
      allocate(vm_h(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
      allocate(up_h(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
      allocate(vp_h(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
      allocate(up_ref(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
      allocate(vp_ref(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))

      allocate(IIc_save,  source=IIc)
      allocate(IIu_save,  source=IIu)
      allocate(IIv_save,  source=IIv)
      allocate(IIcs_save, source=IIcs)
      allocate(IIus_save, source=IIus)
      allocate(IIvs_save, source=IIvs)

      ! Phase 1: a mask the test punches itself rather than whatever the case
      ! happens to carry. Most GPU cases run with libm off, where every mask is
      ! one and a device reduction that ignored the mask entirely would agree
      ! with a correct one - the check has to make its own solid cells to have
      ! anything to say.
      call punch_mask
      IIu_d = IIu
      IIv_d = IIv

      call set_flags(.false., .true., .false., .false.)
      call run_correction
      got = volume_rate(um_h, up_h, IIu, IIus)
      if (abs(got - target_u) > tolerance_for(target_u)) &
         call fail_cuda_selftest('masscorr u volume rate, masked')
      if (maxval(abs(vp_h - vp_ref)) > 1.e-12) &
         call fail_cuda_selftest('masscorr u volume branch moved vp')
      call check_shape('masscorr u volume, masked')

      call set_flags(.false., .false., .false., .true.)
      call run_correction
      got = volume_rate(vm_h, vp_h, IIv, IIvs)
      if (abs(got - target_v) > tolerance_for(target_v)) &
         call fail_cuda_selftest('masscorr v volume rate, masked')
      if (maxval(abs(up_h - up_ref)) > 1.e-12) &
         call fail_cuda_selftest('masscorr v volume branch moved up')
      call check_shape_v('masscorr v volume, masked')

      ! Phase 2: outlet planes clear, so every branch can reach its target.
      IIc = 1 ; IIu = 1 ; IIv = 1
      IIcs = nint(rslabs) ; IIus = nint(rslabs) ; IIvs = nint(rslabs)
      IIu_d = IIu
      IIv_d = IIv
      call calcfluidvolumes

      call set_flags(.false., .true., .false., .false.)
      call run_correction
      got = volume_rate(um_h, up_h, IIu, IIus)
      if (abs(got - target_u) > tolerance_for(target_u)) &
         call fail_cuda_selftest('masscorr u volume rate, clear masks')
      call check_shape('masscorr u volume, clear masks')

      ! Already on target, so nothing more may be added. A defect that is right
      ! once but inconsistent with the correction applied fails here.
      up_ref = up_h
      call rerun_correction
      if (abs(udef) > 1.e-12) call fail_cuda_selftest('masscorr u volume second call defect')
      if (maxval(abs(up_h - up_ref)) > 1.e-12) &
         call fail_cuda_selftest('masscorr u volume second call moved up')

      call set_flags(.false., .false., .false., .true.)
      call run_correction
      got = volume_rate(vm_h, vp_h, IIv, IIvs)
      if (abs(got - target_v) > tolerance_for(target_v)) &
         call fail_cuda_selftest('masscorr v volume rate, clear masks')
      call check_shape_v('masscorr v volume, clear masks')

      call set_flags(.true., .false., .false., .false.)
      call run_correction
      got = outlet_rate_u(um_h, up_h)
      if (abs(got - target_u) > tolerance_for(target_u)) &
         call fail_cuda_selftest('masscorr u outflow rate')
      call check_shape('masscorr u outflow')
      ! uouttot feeds the outflow boundary conditions and is the outlet mass
      ! flow of the tendency before the correction, not after it.
      got = rk3_test * plane_sum_u(up_ref)
      if (abs(uouttot - got) > tolerance_for(got)) &
         call fail_cuda_selftest('masscorr uouttot')

      ! The v outlet is a row of constant j, which lives on one rank once the
      ! domain is split in y while the branch all-reduces over every rank.
      if (nprocy == 1) then
         call set_flags(.false., .false., .true., .false.)
         call run_correction
         got = outlet_rate_v(vm_h, vp_h)
         if (abs(got - target_v) > tolerance_for(target_v)) &
            call fail_cuda_selftest('masscorr v outflow rate')
         call check_shape_v('masscorr v outflow')
      end if

      ! Switches: nothing on, and linoutflow overriding everything.
      call set_flags(.false., .false., .false., .false.)
      call run_correction
      if (maxval(abs(up_h - up_ref)) > 1.e-12 .or. maxval(abs(vp_h - vp_ref)) > 1.e-12) &
         call fail_cuda_selftest('masscorr all switches off')

      call set_flags(.true., .true., .true., .true.)
      linoutflow = .true.
      call run_correction
      if (maxval(abs(up_h - up_ref)) > 1.e-12 .or. maxval(abs(vp_h - vp_ref)) > 1.e-12) &
         call fail_cuda_selftest('masscorr linoutflow')

      IIc = IIc_save ; IIu = IIu_save ; IIv = IIv_save
      IIcs = IIcs_save ; IIus = IIus_save ; IIvs = IIvs_save
      call calcfluidvolumes
      deallocate(IIc_save, IIu_save, IIv_save, IIcs_save, IIus_save, IIvs_save)

      linoutflow = saved_linoutflow
      luoutflowr = saved_uout
      lvoutflowr = saved_vout
      luvolflowr = saved_uvol
      lvvolflowr = saved_vvol
      rk3coef    = saved_rk3coef
      rk3coefi   = saved_rk3coefi
      uflowrate  = saved_uflowrate
      vflowrate  = saved_vflowrate
      udef       = saved_udef
      vdef       = saved_vdef
      uouttot    = saved_uouttot

      deallocate(um_h, vm_h, up_h, vp_h, up_ref, vp_ref)

      ! Give back exactly what was borrowed, and restore the device state the
      ! rest of the run expects.
      if (mine_IIv) deallocate(IIv_d)
      if (mine_IIu) deallocate(IIu_d)
      if (.not. mine_IIu) IIu_d = IIu
      if (.not. mine_IIv) IIv_d = IIv
      if (mine_col) deallocate(col_d, col_stage)
      up_d = up
      vp_d = vp
      um_d = um
      vm_d = vm

   contains

      !> Zero a deterministic third of the interior cells and recount.
      !!
      !! IIus and IIvs are the global fluid-cell count per level, so they come
      !! from an all-reduce of the rank-local count: getting them from anywhere
      !! else would make the target agree with a wrong reduction.
      subroutine punch_mask
         integer :: i, j, kk, jg
         integer :: countu(kb:ke), countv(kb:ke)
         integer :: totalu(kb:ke), totalv(kb:ke)

         IIu = 1
         IIv = 1
         countu = 0
         countv = 0

         do kk = kb, ke
            do j = jb, je
               jg = j + zstart(2) - 1
               do i = ib, ie
                  if (mod(i + 2*jg + 3*kk, 3) == 0) IIu(i,j,kk) = 0
                  if (mod(i + 3*jg + 2*kk, 3) == 0) IIv(i,j,kk) = 0
                  countu(kk) = countu(kk) + IIu(i,j,kk)
                  countv(kk) = countv(kk) + IIv(i,j,kk)
               end do
            end do
         end do

         call MPI_ALLREDUCE(countu, totalu, ke-kb+1, MPI_INTEGER, MPI_SUM, comm3d, mpierr)
         call MPI_ALLREDUCE(countv, totalv, ke-kb+1, MPI_INTEGER, MPI_SUM, comm3d, mpierr)

         IIus = nint(rslabs)
         IIvs = nint(rslabs)
         do kk = kb, ke
            if (totalu(kk) == 0 .or. totalv(kk) == 0) &
               call fail_cuda_selftest('masscorr test mask left a level empty')
            IIus(kk) = totalu(kk)
            IIvs(kk) = totalv(kk)
         end do

      end subroutine punch_mask

      subroutine set_flags(uout, uvol, vout, vvol)
         logical, intent(in) :: uout, uvol, vout, vvol

         luoutflowr = uout
         luvolflowr = uvol
         lvoutflowr = vout
         lvvolflowr = vvol

      end subroutine set_flags

      real function tolerance_for(want)
         real, intent(in) :: want

         tolerance_for = 1.e-10 * max(1., abs(want))

      end function tolerance_for

      !> Fill the device fields, run masscorr and read the tendencies back.
      !!
      !! zstart(2) puts the rank-local j onto the global grid, so the ranks
      !! hold different slabs of one field rather than the same field each: a
      !! reduction that stayed rank-local would otherwise pass.
      subroutine run_correction
         integer :: i, j, kk

         do kk = kb-kh, ke+kh
            do j = jb-jh, je+jh
               do i = ib-ih, ie+ih
                  um_h(i,j,kk) = 0.5 + 0.01*real(i) - 0.003*real(j + zstart(2) - 1) + 0.007*real(kk)
                  vm_h(i,j,kk) = -0.2 + 0.004*real(i) + 0.006*real(j + zstart(2) - 1) - 0.002*real(kk)
               end do
            end do
         end do
         do kk = kb, ke+kh
            do j = jb-jh, je+jh
               do i = ib-ih, ie+ih
                  up_ref(i,j,kk) = 0.02*real(i) + 0.05*real(j + zstart(2) - 1) - 0.01*real(kk)
                  vp_ref(i,j,kk) = -0.03*real(i) + 0.01*real(j + zstart(2) - 1) + 0.04*real(kk)
               end do
            end do
         end do

         um_d = um_h
         vm_d = vm_h
         up_d = up_ref
         vp_d = vp_ref

         call rerun_correction

      end subroutine run_correction

      !> Run masscorr again on whatever the device already holds.
      subroutine rerun_correction

         call masscorr
         call checkCUDA(cudaDeviceSynchronize(), 'masscorr self-test synchronization')
         up_h = up_d
         vp_h = vp_d

      end subroutine rerun_correction

      !> Fluid-volume mean of m + rk3coef*p, formed from the definition.
      real function volume_rate(m, p, II, IIs)
         real,    intent(in) :: m(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
         real,    intent(in) :: p(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
         integer, intent(in) :: II(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc)
         integer, intent(in) :: IIs(kb:ke+khc)
         real    :: local(kb:ke), total(kb:ke)
         integer :: i, j, kk

         local = 0.
         do kk = kb, ke
            do j = jb, je
               do i = ib, ie
                  local(kk) = local(kk) + (m(i,j,kk) + rk3coef*p(i,j,kk)) * II(i,j,kk)
               end do
            end do
         end do

         call MPI_ALLREDUCE(local, total, ke-kb+1, MY_REAL, MPI_SUM, comm3d, mpierr)

         volume_rate = 0.
         do kk = kb, ke
            volume_rate = volume_rate + (total(kk) / IIs(kk)) * dzf(kk)
         end do
         volume_rate = volume_rate / zh(ke+1)

      end function volume_rate

      !> Outlet-plane mean of m + rk3coef*p for the u branch.
      real function outlet_rate_u(m, p)
         real, intent(in) :: m(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
         real, intent(in) :: p(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
         real    :: local(kb:ke), total(kb:ke)
         integer :: j, kk

         local = 0.
         do kk = kb, ke
            do j = jb, je
               local(kk) = local(kk) + (m(ie,j,kk) + rk3coef*p(ie,j,kk)) * IIu(ie,j,kk) * dy
            end do
         end do

         call MPI_ALLREDUCE(local, total, ke-kb+1, MY_REAL, MPI_SUM, comm3d, mpierr)

         outlet_rate_u = 0.
         do kk = kb, ke
            outlet_rate_u = outlet_rate_u + total(kk) * dzf(kk)
         end do
         outlet_rate_u = outlet_rate_u / uoutarea

      end function outlet_rate_u

      !> Outlet-plane sum of one field for the u branch, without the area.
      real function plane_sum_u(p)
         real, intent(in) :: p(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
         real    :: local(kb:ke), total(kb:ke)
         integer :: j, kk

         local = 0.
         do kk = kb, ke
            do j = jb, je
               local(kk) = local(kk) + p(ie,j,kk) * IIu(ie,j,kk) * dy
            end do
         end do

         call MPI_ALLREDUCE(local, total, ke-kb+1, MY_REAL, MPI_SUM, comm3d, mpierr)

         plane_sum_u = 0.
         do kk = kb, ke
            plane_sum_u = plane_sum_u + total(kk) * dzf(kk)
         end do

      end function plane_sum_u

      !> Outlet-row mean of m + rk3coef*p for the v branch.
      real function outlet_rate_v(m, p)
         real, intent(in) :: m(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
         real, intent(in) :: p(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
         real    :: local(kb:ke), total(kb:ke)
         integer :: i, kk

         local = 0.
         do kk = kb, ke
            do i = ib, ie
               local(kk) = local(kk) + (m(i,je,kk) + rk3coef*p(i,je,kk)) * IIv(i,je,kk) * dxf(1)
            end do
         end do

         call MPI_ALLREDUCE(local, total, ke-kb+1, MY_REAL, MPI_SUM, comm3d, mpierr)

         outlet_rate_v = 0.
         do kk = kb, ke
            outlet_rate_v = outlet_rate_v + total(kk) * dzf(kk)
         end do
         outlet_rate_v = outlet_rate_v / voutarea

      end function outlet_rate_v

      !> The u correction must be one constant over ib:ie, jb:je, kb:ke, and
      !! must reach nothing outside it.
      subroutine check_shape(label)
         character(len=*), intent(in) :: label
         real    :: delta
         integer :: i, j, kk

         delta = up_h(ib,jb,kb) - up_ref(ib,jb,kb)
         if (abs(delta) <= 1.e-12) call fail_cuda_selftest(label//': correction was zero')
         do kk = kb, ke+kh
            do j = jb-jh, je+jh
               do i = ib-ih, ie+ih
                  if (i < ib .or. i > ie .or. j < jb .or. j > je .or. kk > ke) then
                     if (abs(up_h(i,j,kk) - up_ref(i,j,kk)) > 1.e-12) &
                        call fail_cuda_selftest(label//': reached outside the interior')
                  else
                     if (abs((up_h(i,j,kk) - up_ref(i,j,kk)) - delta) > 1.e-12) &
                        call fail_cuda_selftest(label//': not uniform')
                  end if
               end do
            end do
         end do

      end subroutine check_shape

      !> As check_shape, for the v tendency.
      subroutine check_shape_v(label)
         character(len=*), intent(in) :: label
         real    :: delta
         integer :: i, j, kk

         delta = vp_h(ib,jb,kb) - vp_ref(ib,jb,kb)
         if (abs(delta) <= 1.e-12) call fail_cuda_selftest(label//': correction was zero')
         do kk = kb, ke+kh
            do j = jb-jh, je+jh
               do i = ib-ih, ie+ih
                  if (i < ib .or. i > ie .or. j < jb .or. j > je .or. kk > ke) then
                     if (abs(vp_h(i,j,kk) - vp_ref(i,j,kk)) > 1.e-12) &
                        call fail_cuda_selftest(label//': reached outside the interior')
                  else
                     if (abs((vp_h(i,j,kk) - vp_ref(i,j,kk)) - delta) > 1.e-12) &
                        call fail_cuda_selftest(label//': not uniform')
                  end if
               end do
            end do
         end do

      end subroutine check_shape_v

   end subroutine test_masscorr

   !> Compare the device ibmnorm kernels against the host originals.
   !!
   !! The host routines are property-tested separately by tests.f90's runmode
   !! 1011, so what is left for the port is that it reproduces them, which is
   !! what this checks - the same shape as test_ibm_diff_corr, on the same
   !! checkerboard geometry so that every branch in the six-neighbour tests is
   !! taken somewhere.
   !!
   !! The solid lists are checked for repeats first. A sequential loop tolerates
   !! a cell appearing twice; one thread per point does not.
   subroutine test_ibmnorm
      implicit none

      real, allocatable :: u0_s(:,:,:), v0_s(:,:,:), w0_s(:,:,:)
      real, allocatable :: um_s(:,:,:), vm_s(:,:,:), wm_s(:,:,:)
      real, allocatable :: up_s(:,:,:), vp_s(:,:,:), wp_s(:,:,:)
      real, allocatable :: thl0_s(:,:,:), thlm_s(:,:,:), thlp_s(:,:,:)
      real, allocatable :: mask_c_s(:,:,:)

      if (.not. allocated(solpts_u_d)) return

      call check_distinct_solid('solpts_u', solid_info_u)
      call check_distinct_solid('solpts_v', solid_info_v)
      call check_distinct_solid('solpts_w', solid_info_w)
      if (allocated(solpts_c_d)) call check_distinct_solid('solpts_c', solid_info_c)

      ! The self-tests run before the time loop, so everything touched here has
      ! to be handed back exactly as it was found, on both host and device.
      allocate(u0_s, source=u0); allocate(v0_s, source=v0); allocate(w0_s, source=w0)
      allocate(um_s, source=um); allocate(vm_s, source=vm); allocate(wm_s, source=wm)
      allocate(up_s, source=up); allocate(vp_s, source=vp); allocate(wp_s, source=wp)

      call seed_field(um, 1.); call seed_field(vm, 2.); call seed_field(wm, 3.)
      call seed_field(up, 4.); call seed_field(vp, 5.); call seed_field(wp, 6.)
      um_d = um; vm_d = vm; wm_d = wm
      up_d = up; vp_d = vp; wp_d = wp

      call solid(solid_info_u, um, up, 0., ih, jh, kh)
      call solid_device(solpts_u_d, solid_info_u%nsolptsrank, um_d, up_d, 0., ih, jh, kh)
      call cmp_m('solid u var', um, um_d, um_s)
      call cmp_p('solid u rhs', up, up_d, up_s)

      call solid(solid_info_v, vm, vp, 0., ih, jh, kh)
      call solid_device(solpts_v_d, solid_info_v%nsolptsrank, vm_d, vp_d, 0., ih, jh, kh)
      call cmp_m('solid v var', vm, vm_d, vm_s)
      call cmp_p('solid v rhs', vp, vp_d, vp_s)

      call solid(solid_info_w, wm, wp, 0., ih, jh, kh)
      call solid_device(solpts_w_d, solid_info_w%nsolptsrank, wm_d, wp_d, 0., ih, jh, kh)
      call cmp_m('solid w var', wm, wm_d, wm_s)
      call cmp_p('solid w rhs', wp, wp_d, wp_s)

      if (allocated(solpts_c_d) .and. ltempeq) then
         allocate(thl0_s, source=thl0)
         allocate(thlm_s, source=thlm)
         allocate(thlp_s, source=thlp)
         allocate(mask_c_s, source=mask_c)

         ! The mask has to agree with the solid list: a cell is solid exactly
         ! when it is listed. The solver builds it that way - createmasks runs
         ! solid over the same list - and the port depends on it, because a
         ! solid point whose neighbour were also listed would be read by the
         ! host after that neighbour had been updated in place and by the
         ! device before. An arbitrary mask would make the two differ for a
         ! reason that cannot arise in a real run.
         call mask_from_solid(mask_c, solid_info_c)
         mask_c_d = mask_c

         call seed_field(thlm, 7.); call seed_field(thlp, 8.)
         thlm_d = thlm; thlp_d = thlp
         call solid(solid_info_c, thlm, thlp, 300., ih, jh, kh, mask_c)
         call solid_masked_device(solpts_c_d, solid_info_c%nsolptsrank, thlm_d, thlp_d, 300., ih, jh, kh)
         call cmp_m('solid masked var', thlm, thlm_d, thlm_s)
         call cmp_p('solid masked rhs', thlp, thlp_d, thlp_s)

         if (allocated(bndpts_c_d)) then
            ! The advection corrections only read the mask and write the
            ! tendency at boundary points, so they carry none of the in-place
            ! neighbour dependency solid has and a checkerboard is safe here.
            ! It is also necessary: case geometries are often uniform in one
            ! direction, and then the faces normal to it are never solid and
            ! two of the six branches are never taken.
            call checkerboard_mask(mask_c)
            mask_c_d = mask_c

            call seed_field(u0, 9.); call seed_field(v0, 10.); call seed_field(w0, 11.)
            u0_d = u0; v0_d = v0; w0_d = w0

            call seed_field(thl0, 12.)
            thl0_d = thl0
            call seed_field(thlp, 13.)
            thlp_d = thlp
            call advecc2nd_corr_conservative(thl0, thlp)
            call advecc2nd_corr_conservative_device(thl0_d, thlp_d)
            call cmp_p('advecc2nd_corr_conservative', thlp, thlp_d, thlp_s)

            call seed_field(thlp, 13.)
            thlp_d = thlp
            call advecc2nd_corr_liberal(thl0, thlp)
            call advecc2nd_corr_liberal_device(thl0_d, thlp_d)
            call cmp_p('advecc2nd_corr_liberal', thlp, thlp_d, thlp_s)
         end if

         mask_c = mask_c_s
         thl0 = thl0_s; thlm = thlm_s; thlp = thlp_s
         mask_c_d = mask_c
         thl0_d = thl0; thlm_d = thlm; thlp_d = thlp
         deallocate(thl0_s, thlm_s, thlp_s, mask_c_s)
      end if

      u0 = u0_s; v0 = v0_s; w0 = w0_s
      um = um_s; vm = vm_s; wm = wm_s
      up = up_s; vp = vp_s; wp = wp_s
      u0_d = u0; v0_d = v0; w0_d = w0
      um_d = um; vm_d = vm; wm_d = wm
      up_d = up; vp_d = vp; wp_d = wp
      deallocate(u0_s, v0_s, w0_s, um_s, vm_s, wm_s, up_s, vp_s, wp_s)

   contains

      !> Fill a field with a value that varies in all three directions.
      subroutine seed_field(a, offset)
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
      end subroutine seed_field

      !> Mark a repeating subset of cells solid, so that every one of the six
      !! neighbour tests in the advection corrections is true somewhere.
      subroutine checkerboard_mask(msk)
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
      end subroutine checkerboard_mask

      !> Build the mask the solid list implies: zero at a listed cell, one
      !! elsewhere. This is what createmasks produces in a real run.
      subroutine mask_from_solid(msk, info)
         implicit none
         real, intent(out) :: msk(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
         type(solid_info_type), intent(in) :: info
         integer :: n

         msk = 1.
         msk(:,:,kb-kh) = 0.
         do n = 1, info%nsolptsrank
            msk(info%solpts_loc(n,1),info%solpts_loc(n,2),info%solpts_loc(n,3)) = 0.
         end do
      end subroutine mask_from_solid

      !> Compare a previous-step field, which spans kb-kh upwards.
      subroutine cmp_m(name, host, dev, before)
         implicit none
         character(len=*), intent(in) :: name
         real, intent(in) :: host  (ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
         real, device, intent(in) :: dev(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
         real, intent(in) :: before(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
         real, allocatable :: back(:,:,:)

         allocate(back(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         back = dev
         call verify(name, host, back, before)
         deallocate(back)
      end subroutine cmp_m

      !> Compare a tendency, which spans kb upwards.
      subroutine cmp_p(name, host, dev, before)
         implicit none
         character(len=*), intent(in) :: name
         real, intent(in) :: host  (ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
         real, device, intent(in) :: dev(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
         real, intent(in) :: before(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
         real, allocatable :: back(:,:,:)

         allocate(back(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
         back = dev
         call verify(name, host, back, before)
         deallocate(back)
      end subroutine cmp_p

      !> Device against host, refusing to pass if the routine did nothing.
      subroutine verify(name, host, back, before)
         implicit none
         character(len=*), intent(in) :: name
         real, intent(in) :: host(:,:,:), back(:,:,:), before(:,:,:)
         real :: worst, scale

         if (any(back /= back)) call fail_cuda_selftest('ibmnorm '//name//' device produced NaN')

         ! A routine that changed nothing would agree with any port of itself.
         if (maxval(abs(host - before)) == 0.) then
            call fail_cuda_selftest('ibmnorm changed nothing: '//name)
         end if

         worst = maxval(abs(back - host))
         scale = max(1., maxval(abs(host)))
         ! Not required to be exact: the device compiler may contract these
         ! expressions into FMAs where the host compiler does not.
         if (worst > 1.e-10 * scale) then
            write(*,'(A,A,A,ES12.4)') 'ibmnorm mismatch ', name, ' worst ', worst
            call fail_cuda_selftest('ibmnorm '//name)
         end if
      end subroutine verify

      !> Assert no cell appears twice in a solid list.
      subroutine check_distinct_solid(name, info)
         implicit none
         character(len=*), intent(in) :: name
         type(solid_info_type), intent(in) :: info
         integer, allocatable :: seen(:,:,:)
         integer :: n

         if (info%nsolptsrank <= 1) return
         allocate(seen(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         seen = 0
         do n = 1, info%nsolptsrank
            if (seen(info%solpts_loc(n,1),info%solpts_loc(n,2),info%solpts_loc(n,3)) /= 0) then
               call fail_cuda_selftest('ibmnorm solid points not distinct: '//name)
            end if
            seen(info%solpts_loc(n,1),info%solpts_loc(n,2),info%solpts_loc(n,3)) = n
         end do
         deallocate(seen)
      end subroutine check_distinct_solid

   end subroutine test_ibmnorm

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
