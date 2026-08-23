!> \file tests_cuda.f90
!! Opt-in CUDA tests driven by the external GPU integration-test harness.

module tests_cuda
#if defined(_GPU) && defined(UDALES_DEBUG)
   use cudafor
   use modadvection, only : advecc_kappa_reset_cuda, advecc_kappa_ducdx_cuda, &
                            advecc_kappa_dvcdy_cuda, advecc_kappa_dwcdz_cuda, &
                            advecc_kappa_add_cuda, advecc_upw_cuda, rlim_cuda
   use modboundary,  only : boundary, boundary_device, &
                            bcpup_pup_BCxm_driver_cuda, &
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
                        uprof_d, vprof_d, thlprof_d, qtprof_d, svprof_d, e12prof_d, &
                        allocDriverPlanesDevice, deallocDriverPlanesDevice, &
                        updateHostForDriverDump, plane_d, plane_stage, &
                        planec_d, planec_stage, &
                        ekm_d, ekh_d, pres0_d, fachf_d, facef_d, fachfi_d, facefi_d, &
                        integrateFacFluxDevice, updateFacIntegralsHost, &
                        facT1_d, facqsat_d, fachurel_d, facf_d, updateFacetPropsDevice, &
                        dxdydzfi_d, dxdydzhi_d, dvcell_d, dxhi_d, dzhi_d, dzf_d, dzfi_d, &
                        col_d, col_stage, IIu_d, IIv_d, updateDevice, &
                        updateHostForFielddump, updateHostForStatsdump, &
                        updateHostForUnportedRoutines, tau_x_d, tau_y_d, tau_z_d, &
                        thl_flux_d, momfluxb_d, tfluxb_d
   use modfields, only : u0, v0, w0, um, vm, wm, e120, e12m, pres0, &
                         thl0, thlm, thl0c, qt0, qtm, sv0, svm, thlp, wp, &
                         up, vp, qtp, svp, &
                         u0av, v0av, thl0av, qt0av, sv0av, &
                         uprof, vprof, thlprof, qtprof, svprof, e12prof, vouttot, &
                         IIc, IIu, IIv, IIcs, IIus, IIvs, &
                         tau_x, tau_y, tau_z, thl_flux, momfluxb, tfluxb, &
                         uoutarea, voutarea, udef, vdef, uouttot
   use modglobal, only : ib, ie, jb, je, kb, ke, ih, jh, kh, nsv, &
                         ihc, jhc, khc, dxhci, dxfc, dxfci, dxi, dyi, &
                         dzhci, dzfc, dzfci, dzfi, eps1, &
                         dx, dy, dxf, dxhi, dzf, dzh, dzhi, zh, rslabs, dxdydzfi, dxdydzhi, dvcell, &
                         lheatpump, lfan_hp, nhppoints, ltempeq, fielddump_wants, &
                         lmoist, lnudge, lnudgevel, tnudge, nnudge, &
                         iwallmom, nfcts, xhat, yhat, zhat, &
                         totheatflux, totqflux, lEB, rk3step, iwallmoist, &
                         lperiodicEBcorr, sinkbase, fraction, xlen, ylen, &
                         linoutflow, luoutflowr, lvoutflowr, luvolflowr, lvvolflowr, &
                         uflowrate, vflowrate, rk3coef, rk3coefi, dt, &
                         ltrees, ltreedump, &
                         BCxm, BCym, BCxm_periodic, BCxm_profile, BCxm_driver, &
                         BCym_periodic, BCym_profile, &
                         BCxT, BCxT_periodic, BCxT_profile, BCxT_driver, &
                         BCyT, BCyT_periodic, BCyT_profile, &
                         BCxq, BCxq_periodic, BCxq_profile, BCxq_driver, &
                         BCyq, BCyq_periodic, BCyq_profile, &
                         BCxs, BCxs_periodic, BCxs_profile, BCxs_driver, BCxs_custom, &
                         BCys, BCys_periodic, &
                         BCtopm, BCtopm_freeslip, BCtopm_noslip, BCtopm_pressure, &
                         BCtopT, BCtopT_flux, BCtopT_value, &
                         BCtopq, BCtopq_flux, BCtopq_value, &
                         BCtops, BCtops_flux, BCtops_value, &
                         Uinf, Vinf, idriver, lchunkread, iadv_thl, iadv_kappa
   use modforces,   only : nudge, periodicEBcorr, masscorr, calcfluidvolumes
   use modchecksim, only : courant_local, diffnr_local, div_local, &
                           diffnrgeom, diffnrgeom_d
   use modheatpump, only : heatpump, nhppoints_local, idhppts_local_d, &
                           thl_dot_hp, w_hp_exhaust
   use initfac,     only : fachf, facef, fachfi, facefi, facT, facqsat, fachurel, facf, faclGR, &
                           lfacetprops_dirty
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
   use modsubgriddata, only : ekm, ekh, loneeqn
   use modsurfdata,  only : wttop, wqtop, wsvtop, thl_top, qt_top, sv_top
   use vegetation,   only : veg, vegp, vegetation_ready, &
                            npts_u, npts_v, npts_w, veg_up, veg_vp, veg_wp, &
                            veg_up_d, veg_vp_d, veg_wp_d, &
                            vegp_qt_d, vegp_qtR_d, vegp_qtA_d, vegp_thl_d, &
                            vegp_omega_d, vegp_sv_d, &
                            vegetation_forcing_device, vegetation_forcing_host, &
                            updateVegDiagHost
   use modinletdata, only : u0driver, umdriver, v0driver, vmdriver, w0driver, wmdriver, &
                            thl0driver, thlmdriver, qt0driver, qtmdriver, &
                            sv0driver, svmdriver, irecydriver
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
         call test_vegetation_forcing
         call test_bc_profile_upload
         call test_boundary_conditions
         call test_driver_plane_handover
         call test_post_poisson_handover
         call test_checksim_reductions
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

      ! Guarded: initdriver allocates the host array under idriver alone, so
      ! on every case in the matrix there is nothing here to copy from.
      if (allocated(u0driver)) then
         u0driver_d = u0driver
      else
         u0driver_d = 0.
      end if
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

      ! The volume itself, which chkdiv integrates over. Exact equality, not a
      ! tolerance: it is meant to be the same product in the same association
      ! as the expression it replaced, not merely close to it.
      allocate(back(kb-kh:ke+kh))
      back = dvcell_d
      if (any(back /= dvcell)) call fail_cuda_selftest('dvcell mirror')
      do k = kb-kh, ke+kh
         if (back(k) /= dx*dy*dzf(k)) call fail_cuda_selftest('dvcell is not dx*dy*dzf')
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

   !> Verify how the energy balance accumulators become a time integral.
   !!
   !! There are three properties and keeping them apart is the whole point,
   !! because conflating them is what made this handover go wrong once already.
   !!
   !! The accumulators must be cleared on EVERY Runge-Kutta stage. The wall
   !! functions add into them each stage; let them run on and the third stage
   !! integrates the sum of all three, which is the facet heat flux three times
   !! too large. That bug shipped, and a manual CPU/GPU comparison found it.
   !!
   !! Only the third stage is integrated, because that is the stage intqH kept
   !! back when it did this on the host. tests.f90::tests_eb checks that claim
   !! against intqH itself on a CPU build.
   !!
   !! And the integration costs no traffic. The time integral is per-rank -
   !! summing over steps and summing over ranks commute, and dt is the same on
   !! every rank - so it accumulates on the device across the hundreds of steps
   !! between energy balances, and only the total comes down. What this test
   !! asserts about that is the arithmetic: after three stages the integral is
   !! one stage's worth times dt, not three.
   subroutine test_facflux_handover
      implicit none

      real, parameter :: probe = 3.25
      real, allocatable :: fachf_s(:), facef_s(:), fachfi_s(:), facefi_s(:)
      real, allocatable :: back(:)
      integer :: saved_rk3step, stage
      real    :: saved_dt, want

      if (.not. lEB) return
      if (.not. allocated(fachf_d)) return
      if (.not. allocated(fachfi_d)) return
      if (.not. allocated(fachfi)) return

      allocate(fachf_s, source=fachf)
      if (allocated(facef)) allocate(facef_s, source=facef)
      allocate(fachfi_s, source=fachfi)
      if (allocated(facefi)) allocate(facefi_s, source=facefi)
      saved_rk3step = rk3step
      saved_dt      = dt
      dt = 0.25
      allocate(back(0:nfcts))

      ! One time step, played through the way the loop runs it: the wall
      ! functions add their share to whatever is on the device, the integration
      ! follows, and nothing crosses the bus.
      fachfi_d = 0.
      facefi_d = 0.
      do stage = 1, 3
         rk3step = stage
         back = fachf_d              ! whatever the last stage left there
         back = back + probe         ! what the wall functions add this stage
         fachf_d = back
         call integrateFacFluxDevice

         ! Cleared on every stage, whichever it was.
         back = fachf_d
         if (any(back /= 0.)) call fail_cuda_selftest('facflux device not reset after integration')
      end do

      ! Three stages in, one stage out. A missing clear lands here at three
      ! times the value; an ungated integration lands at six.
      back = fachfi_d
      want = dt * probe
      if (maxval(abs(back(1:nfcts) - want)) > 64.*epsilon(1.)*want) then
         call fail_cuda_selftest('the facet flux integral counted the wrong number of stages')
      end if

      ! And the integral reaches the host only when it is asked for, exactly
      ! once, leaving the device side clean for the next interval.
      fachfi = 0.
      call updateFacIntegralsHost
      if (maxval(abs(fachfi(1:nfcts) - want)) > 64.*epsilon(1.)*want) then
         call fail_cuda_selftest('the facet flux integral did not reach the host')
      end if
      back = fachfi_d
      if (any(back /= 0.)) call fail_cuda_selftest('facflux integral not reset on handover')

      ! A second call with nothing new must not add to what the host holds:
      ! the download assigns rather than accumulates, so the interval cannot
      ! be counted twice.
      call updateFacIntegralsHost
      if (maxval(abs(fachfi(1:nfcts))) /= 0.) then
         call fail_cuda_selftest('facflux integral counted twice')
      end if

      deallocate(back)
      rk3step = saved_rk3step
      dt      = saved_dt
      fachf   = fachf_s
      fachfi  = fachfi_s
      deallocate(fachf_s, fachfi_s)
      if (allocated(facef_s)) then
         facef = facef_s
         deallocate(facef_s)
      end if
      if (allocated(facefi_s)) then
         facefi = facefi_s
         deallocate(facefi_s)
      end if
      fachf_d  = 0.
      fachfi_d = 0.
      if (allocated(facef_d))  facef_d  = 0.
      if (allocated(facefi_d)) facefi_d = 0.

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


   !> Compare the device vegetation forcing against the host branch it replaced.
   !!
   !! Both branches are compiled into a Debug GPU build, so this feeds them the
   !! same previous-step fields and requires the same answer - the drag on all
   !! three staggered face lists, the canopy tendencies, the scalar deposition,
   !! and the scatter into up/vp/wp/thlp/qtp/svp, compared over the whole
   !! declared extent so a misplaced write cannot hide in a halo.
   !!
   !! It also pins two contracts the port depends on:
   !!   - the device tendency diagnostics are assigned, not accumulated. The
   !!     host clears them through reset_vegetation_sources every step; the
   !!     device has no such pass, and relies on every loop writing every
   !!     element. Running the kernels twice without clearing must therefore
   !!     change nothing.
   !!   - updateVegDiagHost lands exactly what the device holds, which is what
   !!     statsdump reads under ltreedump.
   !!
   !! Skipped on a case without vegetation, like the other configuration-gated
   !! self-tests here.
   subroutine test_vegetation_forcing
      implicit none

      real, parameter :: tol = 1.e-10

      real, allocatable :: um_s(:,:,:), vm_s(:,:,:), wm_s(:,:,:)
      real, allocatable :: thlm_s(:,:,:), qtm_s(:,:,:), svm_s(:,:,:,:)
      real, allocatable :: up_s(:,:,:), vp_s(:,:,:), wp_s(:,:,:)
      real, allocatable :: thlp_s(:,:,:), qtp_s(:,:,:), svp_s(:,:,:,:)
      real, allocatable :: dev3(:,:,:), dev4(:,:,:,:)
      real, allocatable :: host_up(:), host_vp(:), host_wp(:)
      real, allocatable :: host_qt(:), host_qtR(:), host_qtA(:), host_thl(:), host_om(:)
      real, allocatable :: host_sv(:,:)
      real, allocatable :: dev1(:), dev1b(:), dev2(:,:)
      integer :: m, n, npts

      if (.not. ltrees) return
      if (.not. vegetation_ready) return
      npts = veg%npts
      if (npts < 1) return

      allocate(um_s,   source=um)
      allocate(vm_s,   source=vm)
      allocate(wm_s,   source=wm)
      allocate(up_s,   source=up)
      allocate(vp_s,   source=vp)
      allocate(wp_s,   source=wp)
      if (ltempeq) then
         allocate(thlm_s, source=thlm)
         allocate(thlp_s, source=thlp)
      end if
      if (lmoist) then
         allocate(qtm_s, source=qtm)
         allocate(qtp_s, source=qtp)
      end if
      if (nsv > 0) then
         allocate(svm_s, source=svm)
         allocate(svp_s, source=svp)
      end if

      call seed_vegetation_inputs

      ! ---- device branch -------------------------------------------------
      um_d = um
      vm_d = vm
      wm_d = wm
      up_d = up
      vp_d = vp
      wp_d = wp
      if (ltempeq) then
         thlm_d = thlm
         thlp_d = thlp
      end if
      if (lmoist) then
         qtm_d = qtm
         qtp_d = qtp
      end if
      if (nsv > 0) then
         svm_d = svm
         svp_d = svp
      end if

      call vegetation_forcing_device

      allocate(dev1(max(1,npts)), dev1b(max(1,npts)))
      allocate(host_qt(max(1,npts)), host_qtR(max(1,npts)), host_qtA(max(1,npts)))
      allocate(host_thl(max(1,npts)), host_om(max(1,npts)))
      allocate(host_sv(max(1,npts), max(1,nsv)))
      if (npts_u > 0) allocate(host_up(npts_u))
      if (npts_v > 0) allocate(host_vp(npts_v))
      if (npts_w > 0) allocate(host_wp(npts_w))

      ! The idempotence contract: a second pass over unchanged inputs, with the
      ! fields put back but the diagnostics deliberately left alone, must leave
      ! the diagnostics where they were. An accumulating kernel doubles here.
      dev1 = vegp_qt_d
      up_d = up
      vp_d = vp
      wp_d = wp
      if (ltempeq) thlp_d = thlp
      if (lmoist)  qtp_d  = qtp
      if (nsv > 0) svp_d  = svp
      call vegetation_forcing_device
      ! Brought down before comparing: a device array may only appear in a
      ! whole-array assignment here, never inside an expression.
      dev1b = vegp_qt_d
      if (maxval(abs(dev1b - dev1)) > tol*max(1.e-30, maxval(abs(dev1)))) then
         call fail_cuda_selftest('the device vegetation tendencies accumulate across calls')
      end if

      ! ---- host branch, same inputs --------------------------------------
      call vegetation_forcing_host

      ! ---- per-point comparisons -----------------------------------------
      if (npts_u > 0) then
         host_up = veg_up
         deallocate(dev1b); allocate(dev1b(npts_u))
         dev1b = veg_up_d
         call same1('u drag', dev1b, host_up)
      end if
      if (npts_v > 0) then
         host_vp = veg_vp
         deallocate(dev1b); allocate(dev1b(npts_v))
         dev1b = veg_vp_d
         call same1('v drag', dev1b, host_vp)
      end if
      if (npts_w > 0) then
         host_wp = veg_wp
         deallocate(dev1b); allocate(dev1b(npts_w))
         dev1b = veg_wp_d
         call same1('w drag', dev1b, host_wp)
      end if
      deallocate(dev1b); allocate(dev1b(max(1,npts)))

      host_qt  = vegp%qt
      host_qtR = vegp%qtR
      host_qtA = vegp%qtA
      host_thl = vegp%thl
      host_om  = vegp%omega
      dev1 = vegp_qt_d;    call same1('canopy qt',    dev1, host_qt)
      dev1 = vegp_qtR_d;   call same1('canopy qtR',   dev1, host_qtR)
      dev1 = vegp_qtA_d;   call same1('canopy qtA',   dev1, host_qtA)
      dev1 = vegp_thl_d;   call same1('canopy thl',   dev1, host_thl)
      dev1 = vegp_omega_d; call same1('canopy omega', dev1, host_om)

      if (nsv > 0) then
         host_sv = vegp%sv
         allocate(dev2(npts, max(1,nsv)))
         dev2 = vegp_sv_d
         do n = 1, nsv
            call same1('scalar deposition', dev2(:,n), host_sv(:,n))
         end do
         deallocate(dev2)
      end if

      ! ---- scattered fields, whole declared extent ------------------------
      ! The tendencies start at kb, not kb-kh, on both sides.
      allocate(dev3(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh))
      dev3 = up_d; call same3('up', dev3, up)
      dev3 = vp_d; call same3('vp', dev3, vp)
      dev3 = wp_d; call same3('wp', dev3, wp)
      if (ltempeq) then
         dev3 = thlp_d; call same3('thlp', dev3, thlp)
      end if
      if (lmoist) then
         dev3 = qtp_d; call same3('qtp', dev3, qtp)
      end if
      deallocate(dev3)

      if (nsv > 0) then
         allocate(dev4(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb:ke+khc, nsv))
         dev4 = svp_d
         do n = 1, nsv
            if (maxval(abs(dev4(:,:,:,n) - svp(:,:,:,n))) > &
                tol*max(1.e-30, maxval(abs(svp(:,:,:,n))))) then
               call fail_cuda_selftest('device and host svp disagree after vegetation forcing')
            end if
         end do
         deallocate(dev4)
      end if

      ! ---- the diagnostic drain ------------------------------------------
      if (ltreedump) then
         vegp%qt = 0.; vegp%qtR = 0.; vegp%qtA = 0.; vegp%thl = 0.; vegp%omega = 0.
         if (nsv > 0) vegp%sv = 0.
         if (npts_u > 0) veg_up = 0.
         if (npts_v > 0) veg_vp = 0.
         if (npts_w > 0) veg_wp = 0.

         call updateVegDiagHost

         if (lmoist) then
            call same1('drained qt', vegp%qt, host_qt)
            call same1('drained qtR', vegp%qtR, host_qtR)
            call same1('drained qtA', vegp%qtA, host_qtA)
            call same1('drained omega', vegp%omega, host_om)
         end if
         if (ltempeq) call same1('drained thl', vegp%thl, host_thl)
         do n = 1, nsv
            call same1('drained scalar deposition', vegp%sv(:,n), host_sv(:,n))
         end do
         if (npts_u > 0) call same1('drained u drag', veg_up, host_up)
         if (npts_v > 0) call same1('drained v drag', veg_vp, host_vp)
         if (npts_w > 0) call same1('drained w drag', veg_wp, host_wp)
      end if

      ! ---- restore -------------------------------------------------------
      um = um_s; vm = vm_s; wm = wm_s
      up = up_s; vp = vp_s; wp = wp_s
      if (ltempeq) then
         thlm = thlm_s; thlp = thlp_s
      end if
      if (lmoist) then
         qtm = qtm_s; qtp = qtp_s
      end if
      if (nsv > 0) then
         svm = svm_s; svp = svp_s
      end if

   contains

      !> Bounded, index-dependent previous-step fields and cleared tendencies.
      subroutine seed_vegetation_inputs
         integer :: i, j, k, nn
         do k = kb-kh, ke+kh
            do j = jb-jh, je+jh
               do i = ib-ih, ie+ih
                  um(i,j,k) = 0.80 + 0.60*sin(0.70*real(i) + 0.30*real(j) + 0.11*real(k))
                  vm(i,j,k) = -0.30 + 0.50*sin(0.23*real(i) + 0.61*real(j) + 0.17*real(k))
                  wm(i,j,k) = 0.20 + 0.40*sin(0.41*real(i) + 0.13*real(j) + 0.53*real(k))
                  if (ltempeq) thlm(i,j,k) = 292. + 3.0*sin(0.19*real(i) + 0.29*real(j) + 0.37*real(k))
                  if (lmoist)  qtm(i,j,k)  = 0.009 + 0.002*sin(0.31*real(i) + 0.47*real(j) + 0.07*real(k))
               end do
            end do
         end do
         do nn = 1, nsv
            do k = kb-khc, ke+khc
               do j = jb-jhc, je+jhc
                  do i = ib-ihc, ie+ihc
                     svm(i,j,k,nn) = 1.0 + 0.5*real(nn) &
                                   + 0.4*sin(0.43*real(i) + 0.11*real(j) + 0.27*real(k))
                  end do
               end do
            end do
         end do
         up = 0.; vp = 0.; wp = 0.
         if (ltempeq) thlp = 0.
         if (lmoist)  qtp  = 0.
         if (nsv > 0) svp  = 0.
      end subroutine seed_vegetation_inputs

      subroutine same1(label, a, b)
         character(len=*), intent(in) :: label
         real,             intent(in) :: a(:), b(:)
         if (size(a) < 1) return
         if (maxval(abs(a - b)) > tol*max(1.e-30, maxval(abs(b)))) then
            call fail_cuda_selftest('device and host vegetation disagree: '//label)
         end if
      end subroutine same1

      subroutine same3(label, a, b)
         character(len=*), intent(in) :: label
         real,             intent(in) :: a(:,:,:), b(:,:,:)
         if (maxval(abs(a - b)) > tol*max(1.e-30, maxval(abs(b)))) then
            call fail_cuda_selftest('device and host vegetation disagree: '//label)
         end if
      end subroutine same3

   end subroutine test_vegetation_forcing


   !> Pin the inlet uploads that moved out of updateDevicePriorPoiss.
   !!
   !! That routine was the only place copying uprof, vprof and u0driver for the
   !! profile and driver inlets: the nudging guard in updateDevice reaches them
   !! only when velocity nudging happens to be on as well. Deleting the routine
   !! therefore had to move those three, and this is what says so - it puts a
   !! sentinel in each host array, runs the real updateDevice, and reads the
   !! mirror back.
   !!
   !! The switches are set here rather than taken from the case, and lnudge is
   !! forced off so the nudging arm cannot mask the BC one. Each guard is
   !! checked in the configuration that needs it and in one that does not, so
   !! a guard that is simply always true fails rather than passing.
   !!
   !! uprof and vprof travel together, for either inlet. That is not symmetry
   !! for its own sake: xmi_profile sets v at the inlet from vprof and
   !! ymi_profile sets u from uprof, so a guard that uploaded only the profile
   !! named after its own axis would leave the other one at whatever the last
   !! nudging step put there.
   !!
   !! The driver arm matters most. No case in the GPU matrix sets idriver, so
   !! u0driver_d has no parity coverage at all and this is the only thing that
   !! exercises its upload. initdriver allocates the host array under idriver
   !! alone, so where it is missing the test allocates it and gives it back.
   subroutine test_bc_profile_upload
      implicit none

      real, parameter :: sentinel_u = 7.125, sentinel_v = -3.5, sentinel_d = 11.75

      real, parameter :: sentinel_T = 297.25, sentinel_q = 0.015625

      real, allocatable :: uprof_s(:), vprof_s(:), u0driver_s(:,:)
      real, allocatable :: thlprof_s(:), qtprof_s(:)
      real, allocatable :: back1(:), back2(:,:)
      integer :: saved_BCxm, saved_BCym, saved_BCxT, saved_BCyT, saved_BCxq, saved_BCyq
      logical :: saved_lnudge, faked_driver

      if (.not. allocated(uprof_d)) return
      if (.not. allocated(vprof_d)) return

      allocate(uprof_s, source=uprof)
      allocate(vprof_s, source=vprof)
      allocate(back1(lbound(uprof,1):ubound(uprof,1)))
      allocate(back2(jb-jh:je+jh, kb-kh:ke+kh))
      if (allocated(thlprof_d)) allocate(thlprof_s, source=thlprof)
      if (allocated(qtprof_d))  allocate(qtprof_s,  source=qtprof)

      saved_BCxm   = BCxm
      saved_BCym   = BCym
      saved_BCxT   = BCxT
      saved_BCyT   = BCyT
      saved_BCxq   = BCxq
      saved_BCyq   = BCyq
      saved_lnudge = lnudge
      lnudge       = .false.
      BCxT         = BCxT_periodic
      BCyT         = BCyT_periodic
      BCxq         = BCxq_periodic
      BCyq         = BCyq_periodic

      ! --- neither inlet: nothing crosses ----------------------------------
      uprof   = sentinel_u
      vprof   = sentinel_v
      uprof_d = 0.
      vprof_d = 0.
      BCxm    = BCxm_periodic
      BCym    = BCym_periodic
      call updateDevice

      back1 = uprof_d
      if (any(back1 /= 0.)) then
         call fail_cuda_selftest('uprof was uploaded with no profile inlet and no nudging')
      end if
      back1 = vprof_d
      if (any(back1 /= 0.)) then
         call fail_cuda_selftest('vprof was uploaded with no profile inlet and no nudging')
      end if

      ! --- x profile inlet: both profiles cross ----------------------------
      uprof_d = 0.
      vprof_d = 0.
      BCxm    = BCxm_profile
      BCym    = BCym_periodic
      call updateDevice

      back1 = uprof_d
      if (any(back1 /= sentinel_u)) then
         call fail_cuda_selftest('updateDevice did not upload uprof for a profile inlet')
      end if
      back1 = vprof_d
      if (any(back1 /= sentinel_v)) then
         call fail_cuda_selftest('the x inlet guard left vprof behind, which xmi_profile reads')
      end if

      ! --- y profile inlet: the mirror image -------------------------------
      uprof_d = 0.
      vprof_d = 0.
      BCxm    = BCxm_periodic
      BCym    = BCym_profile
      call updateDevice

      back1 = vprof_d
      if (any(back1 /= sentinel_v)) then
         call fail_cuda_selftest('updateDevice did not upload vprof for a profile inlet')
      end if
      back1 = uprof_d
      if (any(back1 /= sentinel_u)) then
         call fail_cuda_selftest('the y inlet guard left uprof behind, which ymi_profile reads')
      end if
      BCym = BCym_periodic

      ! --- the temperature and moisture inlet profiles ---------------------
      if (allocated(thlprof_d)) then
         thlprof   = sentinel_T
         thlprof_d = 0.
         call updateDevice
         back1 = thlprof_d
         if (any(back1 /= 0.)) then
            call fail_cuda_selftest('thlprof was uploaded with no profile inlet and no nudging')
         end if

         BCxT = BCxT_profile
         call updateDevice
         back1 = thlprof_d
         if (any(back1 /= sentinel_T)) then
            call fail_cuda_selftest('updateDevice did not upload thlprof for an x profile inlet')
         end if

         BCxT = BCxT_periodic
         BCyT = BCyT_profile
         thlprof_d = 0.
         call updateDevice
         back1 = thlprof_d
         if (any(back1 /= sentinel_T)) then
            call fail_cuda_selftest('updateDevice did not upload thlprof for a y profile inlet')
         end if
         BCyT = BCyT_periodic
      end if

      if (allocated(qtprof_d)) then
         qtprof   = sentinel_q
         qtprof_d = 0.
         call updateDevice
         back1 = qtprof_d
         if (any(back1 /= 0.)) then
            call fail_cuda_selftest('qtprof was uploaded with no profile inlet and no nudging')
         end if

         BCxq = BCxq_profile
         call updateDevice
         back1 = qtprof_d
         if (any(back1 /= sentinel_q)) then
            call fail_cuda_selftest('updateDevice did not upload qtprof for an x profile inlet')
         end if

         BCxq = BCxq_periodic
         BCyq = BCyq_profile
         qtprof_d = 0.
         call updateDevice
         back1 = qtprof_d
         if (any(back1 /= sentinel_q)) then
            call fail_cuda_selftest('updateDevice did not upload qtprof for a y profile inlet')
         end if
         BCxq = BCxq_periodic
         BCyq = BCyq_periodic
      end if

      ! --- driver inlet -----------------------------------------------------
      faked_driver = .not. allocated(u0driver)
      if (faked_driver) then
         allocate(u0driver(jb-jh:je+jh, kb-kh:ke+kh))
      else
         allocate(u0driver_s, source=u0driver)
      end if

      u0driver   = sentinel_d
      u0driver_d = 0.
      BCxm       = BCxm_driver
      BCym       = BCym_periodic
      call updateDevice

      back2 = u0driver_d
      if (any(back2 /= sentinel_d)) then
         call fail_cuda_selftest('updateDevice did not upload u0driver for a driver inlet')
      end if

      ! --- put everything back ---------------------------------------------
      if (faked_driver) then
         deallocate(u0driver)
         u0driver_d = 0.
      else
         u0driver   = u0driver_s
         u0driver_d = u0driver
         deallocate(u0driver_s)
      end if

      uprof   = uprof_s
      vprof   = vprof_s
      if (allocated(thlprof_s)) then
         thlprof   = thlprof_s
         thlprof_d = thlprof
         deallocate(thlprof_s)
      end if
      if (allocated(qtprof_s)) then
         qtprof   = qtprof_s
         qtprof_d = qtprof
         deallocate(qtprof_s)
      end if
      uprof_d = uprof
      vprof_d = vprof
      BCxm    = saved_BCxm
      BCym    = saved_BCym
      BCxT    = saved_BCxT
      BCyT    = saved_BCyT
      BCxq    = saved_BCxq
      BCyq    = saved_BCyq
      lnudge  = saved_lnudge

      ! updateDevice cleared the facet dirty flag on the way through. Set it
      ! again so the first real call copies the facet properties exactly as it
      ! would have without this test.
      lfacetprops_dirty = .true.

      deallocate(uprof_s, vprof_s, back1, back2)

   end subroutine test_bc_profile_upload

   !> Compare boundary against boundary_device, branch by branch.
   !!
   !! The port of boundary is a large one - thirty-odd small kernels, most of
   !! them writing a single plane - and the ways it can be wrong are all quiet:
   !! a ghost level left at the previous step's value, an index range that
   !! covers the momentum halo where the host covers the scalar one, a
   !! sequential dependence parallelised away. None of that shows up as a
   !! crash, and a parity case only sees it if the branch is one that case
   !! happens to select.
   !!
   !! So this drives both implementations from the same seeded fields and
   !! requires them to agree, walking the branches rather than the
   !! configurations: six passes cover every top condition, every lateral
   !! inlet and every outflow, in a build whose namelist selects one of them.
   !!
   !! Only the host side is saved and put back. The device copies do not need
   !! restoring, because every field written here is one updateDevice uploads
   !! from the host on the first stage of the loop, which is the next thing
   !! that runs.
   subroutine test_boundary_conditions
      implicit none

      integer :: s_BCxm, s_BCym, s_BCxT, s_BCyT, s_BCxq, s_BCyq, s_BCxs, s_BCys
      integer :: s_BCtopm, s_BCtopT, s_BCtopq, s_BCtops, s_rk3step, s_idriver
      real    :: s_rk3coef, s_uouttot, s_vouttot
      real    :: s_wttop, s_wqtop, s_thl_top, s_qt_top, s_Uinf, s_Vinf
      logical :: s_lchunkread, driver_here

      real, allocatable :: h_u0(:,:,:), h_um(:,:,:), h_v0(:,:,:), h_vm(:,:,:), &
                           h_w0(:,:,:), h_wm(:,:,:), h_thl0(:,:,:), h_thlm(:,:,:), &
                           h_thl0c(:,:,:), h_qt0(:,:,:), h_qtm(:,:,:), &
                           h_e120(:,:,:), h_e12m(:,:,:)
      real, allocatable :: h_sv0(:,:,:,:), h_svm(:,:,:,:)
      real, allocatable :: h_ekm(:,:,:), h_ekh(:,:,:)
      real, allocatable :: s_wsvtop(:), s_sv_top(:)
      real, allocatable :: s_uprof(:), s_vprof(:), s_thlprof(:), s_qtprof(:), &
                           s_e12prof(:), s_svprof(:,:)
      integer :: k, n

      call snapshot_host

      s_BCxm = BCxm ; s_BCym = BCym ; s_BCxT = BCxT ; s_BCyT = BCyT
      s_BCxq = BCxq ; s_BCyq = BCyq ; s_BCxs = BCxs ; s_BCys = BCys
      s_BCtopm = BCtopm ; s_BCtopT = BCtopT ; s_BCtopq = BCtopq ; s_BCtops = BCtops
      s_rk3step = rk3step ; s_rk3coef = rk3coef ; s_idriver = idriver
      s_uouttot = uouttot ; s_vouttot = vouttot ; s_lchunkread = lchunkread
      s_wttop = wttop ; s_wqtop = wqtop ; s_thl_top = thl_top ; s_qt_top = qt_top
      s_Uinf = Uinf ; s_Vinf = Vinf
      if (allocated(wsvtop)) then
         allocate(s_wsvtop(size(wsvtop))) ; s_wsvtop = wsvtop
      end if
      if (allocated(sv_top)) then
         allocate(s_sv_top(size(sv_top))) ; s_sv_top = sv_top
      end if

      ! Non-zero fluxes and values throughout, so that a branch that silently
      ! reduced to "copy the level below" would be caught rather than agree.
      rk3step = 3
      rk3coef = 0.375
      idriver = 0
      lchunkread = .false.
      wttop = -2.5e-3 ; wqtop = 1.25e-4
      thl_top = 291.5 ; qt_top = 6.25e-3
      Uinf = 1.5 ; Vinf = -0.75
      if (allocated(wsvtop)) wsvtop = 3.125e-3
      if (allocated(sv_top)) sv_top = 0.625

      ! The inlet profiles the lateral conditions read. They are given values
      ! here rather than taken from the case, because a case that never uses a
      ! profile inlet leaves them at zero and a routine that read the wrong one
      ! would then agree with one that read the right one.
      allocate(s_uprof, source=uprof) ; allocate(s_vprof, source=vprof)
      do k = lbound(uprof,1), ubound(uprof,1)
         uprof(k) =  1.25 + 0.03125*real(k)
         vprof(k) = -0.50 + 0.015625*real(k)
      end do
      uprof_d = uprof ; vprof_d = vprof

      if (allocated(thlprof_d)) then
         allocate(s_thlprof, source=thlprof)
         do k = lbound(thlprof,1), ubound(thlprof,1)
            thlprof(k) = 289.0 + 0.0625*real(k)
         end do
         thlprof_d = thlprof
      end if
      if (allocated(qtprof_d)) then
         allocate(s_qtprof, source=qtprof)
         do k = lbound(qtprof,1), ubound(qtprof,1)
            qtprof(k) = 0.25 + 0.0078125*real(k)
         end do
         qtprof_d = qtprof
      end if
      if (allocated(e12prof_d)) then
         allocate(s_e12prof, source=e12prof)
         do k = lbound(e12prof,1), ubound(e12prof,1)
            e12prof(k) = 0.125 + 0.00390625*real(k)
         end do
         e12prof_d = e12prof
      end if
      if (allocated(svprof_d)) then
         allocate(s_svprof, source=svprof)
         do n = 1, nsv
            do k = lbound(svprof,1), ubound(svprof,1)
               svprof(k,n) = 0.5*real(n) + 0.03125*real(k)
            end do
         end do
         svprof_d = svprof
      end if

      ! Both eddy diffusivities are read by the flux-type top conditions and
      ! written by nothing here. They are given values rather than inherited:
      ! initCUDA allocates the device copies without seeding them, and the host
      ! copies hold whatever initsubgrid left, which on some cases is small
      ! enough to make the flux term swamp the field it is added to.
      do k = kb-kh, ke+kh
         ekm(:,:,k) = 0.25  + 0.0078125*real(k - kb)
         ekh(:,:,k) = 0.125 + 0.00390625*real(k - kb)
      end do
      ekm_d = ekm
      ekh_d = ekh

      ! 1. Zero-flux top, flux-type scalar tops, everything periodic sideways.
      BCtopm = BCtopm_freeslip
      BCtopT = BCtopT_flux ; BCtopq = BCtopq_flux ; BCtops = BCtops_flux
      BCxm = BCxm_periodic ; BCxT = BCxT_periodic ; BCxq = BCxq_periodic ; BCxs = BCxs_periodic
      BCym = BCym_periodic ; BCyT = BCyT_periodic ; BCyq = BCyq_periodic ; BCys = BCys_periodic
      call run_case('freeslip-flux-tops')

      ! 2. Fixed-velocity top and fixed-value scalar tops.
      BCtopm = BCtopm_noslip
      BCtopT = BCtopT_value ; BCtopq = BCtopq_value ; BCtops = BCtops_value
      call run_case('noslip-value-tops')

      ! 3. Pressure top, which leaves w to modpois, with the x profile inlet
      !    and the convective outflow it implies.
      BCtopm = BCtopm_pressure
      BCtopT = BCtopT_flux ; BCtopq = BCtopq_flux ; BCtops = BCtops_flux
      BCxm = BCxm_profile ; BCxT = BCxT_profile ; BCxq = BCxq_profile ; BCxs = BCxs_profile
      call run_case('x-profile-inlet')

      ! 4. The same in y.
      BCxm = BCxm_periodic ; BCxT = BCxT_periodic ; BCxq = BCxq_periodic ; BCxs = BCxs_periodic
      BCym = BCym_profile ; BCyT = BCyT_profile ; BCyq = BCyq_profile ; BCys = 2
      call run_case('y-profile-inlet')

      ! 5. The single-column scalar inlet, which is the one branch that has to
      !    find its own j rather than being handed a range.
      BCym = BCym_periodic ; BCyT = BCyT_periodic ; BCyq = BCyq_periodic ; BCys = BCys_periodic
      BCxs = BCxs_custom
      call run_case('x-custom-scalar-inlet')

      ! 6. The driver inlet, if this run is not itself reading a driver file.
      !    When it is, the planes hold data the first step needs and the run
      !    exercises the branch anyway.
      BCxs = BCxs_periodic
      driver_here = .not. allocated(u0driver)
      if (driver_here) then
         call open_driver_planes
         BCxm = BCxm_driver ; BCxT = BCxT_driver ; BCxq = BCxq_driver ; BCxs = BCxs_driver
         call run_case('x-driver-inlet')
         call close_driver_planes
      end if

      BCxm = s_BCxm ; BCym = s_BCym ; BCxT = s_BCxT ; BCyT = s_BCyT
      BCxq = s_BCxq ; BCyq = s_BCyq ; BCxs = s_BCxs ; BCys = s_BCys
      BCtopm = s_BCtopm ; BCtopT = s_BCtopT ; BCtopq = s_BCtopq ; BCtops = s_BCtops
      rk3step = s_rk3step ; rk3coef = s_rk3coef ; idriver = s_idriver
      uouttot = s_uouttot ; vouttot = s_vouttot ; lchunkread = s_lchunkread
      wttop = s_wttop ; wqtop = s_wqtop ; thl_top = s_thl_top ; qt_top = s_qt_top
      Uinf = s_Uinf ; Vinf = s_Vinf
      if (allocated(s_wsvtop)) wsvtop = s_wsvtop
      if (allocated(s_sv_top)) sv_top = s_sv_top

      uprof = s_uprof ; vprof = s_vprof
      uprof_d = uprof ; vprof_d = vprof
      if (allocated(s_thlprof)) then
         thlprof = s_thlprof ; thlprof_d = thlprof
      end if
      if (allocated(s_qtprof)) then
         qtprof = s_qtprof ; qtprof_d = qtprof
      end if
      if (allocated(s_e12prof)) then
         e12prof = s_e12prof ; e12prof_d = e12prof
      end if
      if (allocated(s_svprof)) then
         svprof = s_svprof ; svprof_d = svprof
      end if

      call restore_host

   contains

      subroutine snapshot_host
         implicit none

         allocate(h_u0, source=u0) ; allocate(h_um, source=um)
         allocate(h_v0, source=v0) ; allocate(h_vm, source=vm)
         allocate(h_w0, source=w0) ; allocate(h_wm, source=wm)
         allocate(h_thl0, source=thl0) ; allocate(h_thlm, source=thlm)
         allocate(h_thl0c, source=thl0c)
         allocate(h_qt0, source=qt0) ; allocate(h_qtm, source=qtm)
         allocate(h_ekm, source=ekm) ; allocate(h_ekh, source=ekh)
         if (allocated(e120)) then
            allocate(h_e120, source=e120) ; allocate(h_e12m, source=e12m)
         end if
         if (nsv > 0) then
            allocate(h_sv0, source=sv0) ; allocate(h_svm, source=svm)
         end if
      end subroutine snapshot_host

      subroutine restore_host
         implicit none

         u0 = h_u0 ; um = h_um ; v0 = h_v0 ; vm = h_vm ; w0 = h_w0 ; wm = h_wm
         thl0 = h_thl0 ; thlm = h_thlm ; thl0c = h_thl0c
         qt0 = h_qt0 ; qtm = h_qtm

         ! The device copies of the diffusivities go back too. Everything else
         ! written here is uploaded again by the first updateDevice of the
         ! loop; these two are not - nothing uploads them, so a seeded value
         ! left behind would survive into the run wherever subgrid does not
         ! write, which is every halo cell.
         ekm = h_ekm ; ekh = h_ekh
         ekm_d = ekm ; ekh_d = ekh
         if (allocated(h_e120)) then
            e120 = h_e120 ; e12m = h_e12m
         end if
         if (allocated(h_sv0)) then
            sv0 = h_sv0 ; svm = h_svm
         end if
      end subroutine restore_host

      !> Deterministic, dyadic and distinct per field.
      !!
      !! Increments are negative powers of two so that the halving and doubling
      !! the value-type conditions do is exact, and every field gets its own
      !! offset so that a routine writing into the wrong one cannot pass.
      subroutine seed_fields
         implicit none
         integer :: i, j, k, n

         do k = kb-kh, ke+kh
            do j = jb-jh, je+jh
               do i = ib-ih, ie+ih
                  u0(i,j,k)   =  1.0  + at(i,j,k)
                  um(i,j,k)   =  2.0  + at(i,j,k)
                  v0(i,j,k)   =  3.0  + at(i,j,k)
                  vm(i,j,k)   =  4.0  + at(i,j,k)
                  w0(i,j,k)   =  5.0  + at(i,j,k)
                  wm(i,j,k)   =  6.0  + at(i,j,k)
                  thl0(i,j,k) = 290.0 + at(i,j,k)
                  thlm(i,j,k) = 291.0 + at(i,j,k)
                  qt0(i,j,k)  =  0.5  + at(i,j,k)
                  qtm(i,j,k)  =  0.75 + at(i,j,k)
               end do
            end do
         end do

         if (allocated(e120)) then
            do k = kb-kh, ke+kh
               do j = jb-jh, je+jh
                  do i = ib-ih, ie+ih
                     e120(i,j,k) = 7.0 + at(i,j,k)
                     e12m(i,j,k) = 8.0 + at(i,j,k)
                  end do
               end do
            end do
         end if

         do k = kb-khc, ke+khc
            do j = jb-jhc, je+jhc
               do i = ib-ihc, ie+ihc
                  thl0c(i,j,k) = 292.0 + at(i,j,k)
                  do n = 1, nsv
                     sv0(i,j,k,n) = 0.25*real(n) + 9.0  + at(i,j,k)
                     svm(i,j,k,n) = 0.25*real(n) + 10.0 + at(i,j,k)
                  end do
               end do
            end do
         end do

         u0_d = u0 ; um_d = um ; v0_d = v0 ; vm_d = vm ; w0_d = w0 ; wm_d = wm
         if (allocated(thl0_d)) then
            thl0_d = thl0 ; thlm_d = thlm
         end if
         if (allocated(thl0c_d)) thl0c_d = thl0c
         if (allocated(qt0_d)) then
            qt0_d = qt0 ; qtm_d = qtm
         end if
         if (allocated(e120_d)) then
            e120_d = e120 ; e12m_d = e12m
         end if
         if (allocated(sv0_d)) then
            sv0_d = sv0 ; svm_d = svm
         end if
      end subroutine seed_fields

      real function at(i, j, k)
         implicit none
         integer, intent(in) :: i, j, k

         at = 0.125*real(i) + 0.0625*real(j + zstart(2) - 1) + 0.03125*real(k)
      end function at

      subroutine run_case(label)
         implicit none
         character(len=*), intent(in) :: label
         integer :: n

         call seed_fields

         call boundary
         call boundary_device
         call checkCUDA(cudaDeviceSynchronize(), 'boundary device self-test: '//label)

         call check3(label, 'u0', u0, u0_d, ib-ih, jb-jh, kb-kh)
         call check3(label, 'um', um, um_d, ib-ih, jb-jh, kb-kh)
         call check3(label, 'v0', v0, v0_d, ib-ih, jb-jh, kb-kh)
         call check3(label, 'vm', vm, vm_d, ib-ih, jb-jh, kb-kh)
         call check3(label, 'w0', w0, w0_d, ib-ih, jb-jh, kb-kh)
         call check3(label, 'wm', wm, wm_d, ib-ih, jb-jh, kb-kh)
         ! On the physics switches, not on what is allocated. The host applies
         ! the top conditions to temperature and moisture whether or not the
         ! run solves for them; the device routines return instead, which is
         ! the same choice every other _device twin in modboundary makes and
         ! costs nothing because nothing reads a field the run does not solve.
         ! Allocation is not the same question: an earlier self-test assigns
         ! into thl0_d unguarded, and Fortran allocates on assignment.
         if (ltempeq) then
            call check3(label, 'thl0', thl0, thl0_d, ib-ih, jb-jh, kb-kh)
            call check3(label, 'thlm', thlm, thlm_d, ib-ih, jb-jh, kb-kh)
            if (iadv_thl == iadv_kappa) &
               call check3(label, 'thl0c', thl0c, thl0c_d, ib-ihc, jb-jhc, kb-khc)
         end if
         if (lmoist) then
            call check3(label, 'qt0', qt0, qt0_d, ib-ih, jb-jh, kb-kh)
            call check3(label, 'qtm', qtm, qtm_d, ib-ih, jb-jh, kb-kh)
         end if
         if (loneeqn) then
            call check3(label, 'e120', e120, e120_d, ib-ih, jb-jh, kb-kh)
            call check3(label, 'e12m', e12m, e12m_d, ib-ih, jb-jh, kb-kh)
         end if
         if (nsv > 0) then
            do n = 1, nsv
               call check3(label, 'sv0', sv0(:,:,:,n), sv0_d(:,:,:,n), ib-ihc, jb-jhc, kb-khc)
               call check3(label, 'svm', svm(:,:,:,n), svm_d(:,:,:,n), ib-ihc, jb-jhc, kb-khc)
            end do
         end if
      end subroutine run_case

      !> Require host and device to agree, and say where they do not.
      !!
      !! The lower bounds travel as arguments because an assumed-shape dummy
      !! renumbers them from one, and a report that named the wrong cell would
      !! send the next reader to the wrong branch.
      subroutine check3(label, name, host, dev, lbi, lbj, lbk)
         implicit none
         character(len=*), intent(in) :: label, name
         real,             intent(in) :: host(:,:,:)
         real, device,     intent(in) :: dev(:,:,:)
         integer,          intent(in) :: lbi, lbj, lbk

         real, allocatable :: back(:,:,:)
         real :: worst, scale
         integer :: at3(3)

         allocate(back(size(dev,1), size(dev,2), size(dev,3)))
         back = dev

         worst = maxval(abs(back - host))
         scale = max(1., maxval(abs(host)))
         if (worst > 1.e-12*scale) then
            at3 = maxloc(abs(back - host))
            write(*,'(4A)') 'boundary parity: ', trim(name), ' differs at ', trim(label)
            write(*,'(A,I0,A,I0,A,I0,A,ES22.14,A,ES22.14,A,I0)') &
                 '  worst cell (i,j,k) = (', at3(1) + lbi - 1, ', ', at3(2) + lbj - 1, &
                 ', ', at3(3) + lbk - 1, ')  host ', host(at3(1), at3(2), at3(3)), &
                 '  device ', back(at3(1), at3(2), at3(3)), '  differing cells ', &
                 count(abs(back - host) > 1.e-12*scale)
            deallocate(back)
            call fail_cuda_selftest('boundary '//trim(label)//' '//trim(name))
         end if

         deallocate(back)
      end subroutine check3

      !> Give the driver branch planes to read, on both sides.
      subroutine open_driver_planes
         implicit none
         integer :: j, k, n

         allocate(u0driver(jb-jh:je+jh,kb-kh:ke+kh))
         allocate(umdriver(jb-jh:je+jh,kb-kh:ke+kh))
         allocate(v0driver(jb-jh:je+jh,kb-kh:ke+kh))
         allocate(vmdriver(jb-jh:je+jh,kb-kh:ke+kh))
         allocate(w0driver(jb-jh:je+jh,kb-kh:ke+kh))
         allocate(wmdriver(jb-jh:je+jh,kb-kh:ke+kh))
         if (ltempeq) then
            allocate(thl0driver(jb-jh:je+jh,kb-kh:ke+kh))
            allocate(thlmdriver(jb-jh:je+jh,kb-kh:ke+kh))
         end if
         if (lmoist) then
            allocate(qt0driver(jb-jh:je+jh,kb-kh:ke+kh))
            allocate(qtmdriver(jb-jh:je+jh,kb-kh:ke+kh))
         end if
         if (nsv > 0) then
            allocate(sv0driver(jb-jhc:je+jhc,kb-khc:ke+khc,1:nsv))
            allocate(svmdriver(jb-jhc:je+jhc,kb-khc:ke+khc,1:nsv))
         end if

         do k = kb-kh, ke+kh
            do j = jb-jh, je+jh
               u0driver(j,k) = 1.5  + at(0,j,k)
               umdriver(j,k) = 2.5  + at(0,j,k)
               v0driver(j,k) = 3.5  + at(0,j,k)
               vmdriver(j,k) = 4.5  + at(0,j,k)
               w0driver(j,k) = 5.5  + at(0,j,k)
               wmdriver(j,k) = 6.5  + at(0,j,k)
               if (ltempeq) then
                  thl0driver(j,k) = 293.5 + at(0,j,k)
                  thlmdriver(j,k) = 294.5 + at(0,j,k)
               end if
               if (lmoist) then
                  qt0driver(j,k) = 0.875 + at(0,j,k)
                  qtmdriver(j,k) = 0.9375 + at(0,j,k)
               end if
            end do
         end do
         do n = 1, nsv
            do k = kb-khc, ke+khc
               do j = jb-jhc, je+jhc
                  sv0driver(j,k,n) = 11.0 + 0.25*real(n) + at(0,j,k)
                  svmdriver(j,k,n) = 12.0 + 0.25*real(n) + at(0,j,k)
               end do
            end do
         end do

         call allocDriverPlanesDevice
      end subroutine open_driver_planes

      subroutine close_driver_planes
         implicit none

         call deallocDriverPlanesDevice
         deallocate(u0driver, umdriver, v0driver, vmdriver, w0driver, wmdriver)
         if (allocated(thl0driver)) deallocate(thl0driver, thlmdriver)
         if (allocated(qt0driver))  deallocate(qt0driver, qtmdriver)
         if (allocated(sv0driver))  deallocate(sv0driver, svmdriver)
      end subroutine close_driver_planes

   end subroutine test_boundary_conditions

   !> Check that the driver recording plane comes down, and the right one.
   !!
   !! A driver-generation run records one i-plane of each prognostic field
   !! every dtdriver, and the plane is now gathered on the device rather than
   !! the whole field being fetched. Nothing else covers it: the parity
   !! harness compares netCDF outputs and a driver file is raw binary, and no
   !! case in the matrix sets idriver at all.
   !!
   !! What it would catch is the pair of mistakes this rewrite could make -
   !! gathering the plane below the one writedriverfile writes, or handing
   !! back a buffer that was never filled. The device planes are given values
   !! the host does not have, so a routine that quietly did nothing leaves the
   !! sentinel behind.
   subroutine test_driver_plane_handover
      implicit none

      integer :: i, j, k, n, saved_irecy
      logical :: made_plane, made_planec
      real, parameter :: sentinel = -98765.
      real, allocatable :: want(:,:), got(:,:)
      real, allocatable :: s_u0(:,:), s_v0(:,:), s_w0(:,:), s_thl0(:,:), s_qt0(:,:)
      real, allocatable :: s_sv0(:,:,:)

      ! The recycle plane and the one below it both have to be inside the
      ! domain, which they are for any iplane a real run would use.
      saved_irecy = irecydriver
      irecydriver = ib + 1

      made_plane = .not. allocated(plane_d)
      if (made_plane) then
         allocate(plane_d(jb-jh:je+jh, kb-kh:ke+kh))
         allocate(plane_stage(jb-jh:je+jh, kb-kh:ke+kh))
      end if
      made_planec = (nsv > 0) .and. (.not. allocated(planec_d))
      if (made_planec) then
         allocate(planec_d(jb-jhc:je+jhc, kb-khc:ke+khc))
         allocate(planec_stage(jb-jhc:je+jhc, kb-khc:ke+khc))
      end if

      allocate(want(jb-jh:je+jh, kb-kh:ke+kh))
      allocate(got (jb-jh:je+jh, kb-kh:ke+kh))

      ! Keep what is about to be overwritten on the host.
      allocate(s_u0,  source=u0(irecydriver, :, :))
      allocate(s_v0,  source=v0(irecydriver-1, :, :))
      allocate(s_w0,  source=w0(irecydriver-1, :, :))
      if (ltempeq) allocate(s_thl0, source=thl0(irecydriver-1, :, :))
      if (lmoist)  allocate(s_qt0,  source=qt0(irecydriver-1, :, :))
      if (nsv > 0) allocate(s_sv0,  source=sv0(irecydriver-1, :, :, :))

      call check_plane(irecydriver,   10., 'u0')
      call check_plane(irecydriver-1, 20., 'v0')
      call check_plane(irecydriver-1, 30., 'w0')
      if (ltempeq) call check_plane(irecydriver-1, 40., 'thl0')
      if (lmoist)  call check_plane(irecydriver-1, 50., 'qt0')
      do n = 1, nsv
         call check_plane_sv(irecydriver-1, 60. + real(n), n)
      end do

      ! Put the host back. The device copies are restored by the first
      ! updateDevice of the loop, which uploads every field touched here.
      u0(irecydriver, :, :)   = s_u0
      v0(irecydriver-1, :, :) = s_v0
      w0(irecydriver-1, :, :) = s_w0
      if (allocated(s_thl0)) thl0(irecydriver-1, :, :) = s_thl0
      if (allocated(s_qt0))  qt0(irecydriver-1, :, :)  = s_qt0
      if (allocated(s_sv0))  sv0(irecydriver-1, :, :, :) = s_sv0

      if (made_planec) deallocate(planec_d, planec_stage)
      if (made_plane)  deallocate(plane_d, plane_stage)
      deallocate(want, got)
      irecydriver = saved_irecy

   contains

      !> Put a pattern on one device plane, sentinel the host, and fetch.
      !!
      !! Each field is written and read one at a time, so a gather that named
      !! the wrong array would return the previous field's pattern rather than
      !! one that happens to look plausible.
      subroutine check_plane(iplane, base, name)
         implicit none
         integer,          intent(in) :: iplane
         real,             intent(in) :: base
         character(len=*), intent(in) :: name

         do k = kb - kh, ke + kh
            do j = jb - jh, je + jh
               want(j, k) = base + 0.125*real(j + zstart(2) - 1) + 0.03125*real(k)
            end do
         end do

         ! The pattern goes onto the device through the host copy and a
         ! whole-array upload, then the host plane is overwritten. Writing the
         ! device plane directly would be a strided host-to-device section
         ! assignment, which is exactly the transfer this rewrite exists to
         ! avoid and not something the test should depend on.
         select case (name)
         case ('u0')
            u0(iplane, :, :) = want ; u0_d = u0 ; u0(iplane, :, :) = sentinel
         case ('v0')
            v0(iplane, :, :) = want ; v0_d = v0 ; v0(iplane, :, :) = sentinel
         case ('w0')
            w0(iplane, :, :) = want ; w0_d = w0 ; w0(iplane, :, :) = sentinel
         case ('thl0')
            thl0(iplane, :, :) = want ; thl0_d = thl0 ; thl0(iplane, :, :) = sentinel
         case ('qt0')
            qt0(iplane, :, :) = want ; qt0_d = qt0 ; qt0(iplane, :, :) = sentinel
         case default
            call fail_cuda_selftest('driver plane self-test: unknown field '//trim(name))
         end select

         call updateHostForDriverDump

         select case (name)
         case ('u0')
            got = u0(iplane, :, :)
         case ('v0')
            got = v0(iplane, :, :)
         case ('w0')
            got = w0(iplane, :, :)
         case ('thl0')
            got = thl0(iplane, :, :)
         case ('qt0')
            got = qt0(iplane, :, :)
         end select

         if (any(got /= want)) then
            call fail_cuda_selftest('driver plane '//trim(name))
         end if
      end subroutine check_plane

      subroutine check_plane_sv(iplane, base, species)
         implicit none
         integer, intent(in) :: iplane, species
         real,    intent(in) :: base

         real, allocatable :: wantc(:,:), gotc(:,:)

         allocate(wantc(jb-jhc:je+jhc, kb-khc:ke+khc))
         allocate(gotc (jb-jhc:je+jhc, kb-khc:ke+khc))

         do k = kb - khc, ke + khc
            do j = jb - jhc, je + jhc
               wantc(j, k) = base + 0.125*real(j + zstart(2) - 1) + 0.03125*real(k)
            end do
         end do

         sv0(iplane, :, :, species) = wantc
         sv0_d = sv0
         sv0(iplane, :, :, species) = sentinel

         call updateHostForDriverDump

         gotc = sv0(iplane, :, :, species)
         if (any(gotc /= wantc)) then
            call fail_cuda_selftest('driver plane sv0')
         end if

         deallocate(wantc, gotc)
      end subroutine check_plane_sv

   end subroutine test_driver_plane_handover


   !> Pin the two contracts left behind by folding updateHost away.
   !!
   !! updateHost used to run before the pressure step and bring down eighteen
   !! fields. Ten of them had no host reader at all in a GPU build and were
   !! dropped; the six that do have one now come down in the post-Poisson
   !! handover,
   !! which is what the second half of this checks. Nothing writes them on the
   !! device in between, so the values have to be the ones the device holds.
   !!
   !! The first half is the trap. momfluxb and tfluxb are accumulators - the
   !! wall functions add into them and nothing ever clears them, on either side
   !! - so the download/upload round trip was what carried their running total
   !! across stages. Dropping the download alone would have frozen them at
   !! whatever the host last held. Both halves went, and initCUDA seeds the
   !! device copies once, which leaves them accumulating on the device exactly
   !! as they do on the host. That only holds while updateDevice keeps its
   !! hands off them, which is what this asserts.
   subroutine test_post_poisson_handover
      implicit none

      real, parameter :: sent_mom = 41.5, sent_tf  = -17.25
      real, parameter :: sent_ekm = 3.125, sent_ekh = 6.375
      real, parameter :: sent_tx = 11.5, sent_ty = -22.25, sent_tz = 33.75
      real, parameter :: sent_hf = -44.125

      real, allocatable :: ekm_s(:,:,:), ekh_s(:,:,:)
      real, allocatable :: tau_x_s(:,:,:), tau_y_s(:,:,:), tau_z_s(:,:,:)
      real, allocatable :: thl_flux_s(:,:,:)
      real, allocatable :: momfluxb_s(:,:,:), tfluxb_s(:,:,:)
      real, allocatable :: mom_dev(:,:,:), tf_dev(:,:,:), back(:,:,:)

      if (.not. allocated(momfluxb_d)) return
      if (.not. allocated(ekm_d)) return

      allocate(back(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))

      allocate(ekm_s,      source=ekm)
      allocate(ekh_s,      source=ekh)
      allocate(tau_x_s,    source=tau_x)
      allocate(tau_y_s,    source=tau_y)
      allocate(tau_z_s,    source=tau_z)
      allocate(momfluxb_s, source=momfluxb)
      if (ltempeq) then
         allocate(thl_flux_s, source=thl_flux)
         allocate(tfluxb_s,   source=tfluxb)
      end if

      ! What the device currently holds, so the accumulators can be handed back
      ! untouched: the run continues from here.
      allocate(mom_dev(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
      mom_dev = momfluxb_d
      if (ltempeq) then
         allocate(tf_dev(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
         tf_dev = tfluxb_d
      end if

      ! Make every mirror updateDevice uploads agree with the host, so that the
      ! handover calls below are an identity for everything except
      ! the six fields under test.
      call updateDevice

      ! --- the accumulators survive updateDevice ---------------------------
      momfluxb_d = sent_mom
      momfluxb   = sent_mom + 1.
      if (ltempeq) then
         tfluxb_d = sent_tf
         tfluxb   = sent_tf + 1.
      end if

      call updateDevice

      back = momfluxb_d
      if (any(back /= sent_mom)) then
         call fail_cuda_selftest('updateDevice overwrote the momentum flux accumulator')
      end if
      if (ltempeq) then
         back = tfluxb_d
         if (any(back /= sent_tf)) then
            call fail_cuda_selftest('updateDevice overwrote the heat flux accumulator')
         end if
      end if

      ! --- the handover brings the six survivors down -----------------------
      ekm_d   = sent_ekm ; ekm   = 0.
      ekh_d   = sent_ekh ; ekh   = 0.
      tau_x_d = sent_tx  ; tau_x = 0.
      tau_y_d = sent_ty  ; tau_y = 0.
      tau_z_d = sent_tz  ; tau_z = 0.
      if (ltempeq) then
         thl_flux_d = sent_hf ; thl_flux = 0.
      end if

      ! The post-Poisson handover is four routines now, one per reader, and
      ! between them they still have to bring the whole set down. Called in the
      ! order the time loop calls them, so the pull bookkeeping is exercised
      ! the way it runs: whichever routine reaches a field first fetches it and
      ! the others leave it alone.
      call updateHostForFielddump
      call updateHostForStatsdump
      call updateHostForUnportedRoutines

      if (any(ekm   /= sent_ekm)) call fail_cuda_selftest('the handover did not bring down ekm')
      if (any(ekh   /= sent_ekh)) call fail_cuda_selftest('the handover did not bring down ekh')

      ! The stresses and the heat flux are pulled only when fielddump names
      ! them, so both directions are asserted: brought down when there is a
      ! reader, and left alone when there is not. The host sentinel is already
      ! zero, so the negative half costs nothing, and it is the only check
      ! anywhere that the guard really skips rather than merely being written
      ! down. The parity matrix exercises both branches - fielddump-stress-
      ! fields asks for all four, every other case asks for none.
      call check_gated('tx', tau_x, sent_tx, 'tau_x')
      call check_gated('ty', tau_y, sent_ty, 'tau_y')
      call check_gated('tz', tau_z, sent_tz, 'tau_z')
      if (ltempeq) call check_gated('hf', thl_flux, sent_hf, 'thl_flux')

      ! --- put everything back ---------------------------------------------
      ekm   = ekm_s   ; ekm_d   = ekm
      ekh   = ekh_s   ; ekh_d   = ekh
      tau_x = tau_x_s ; tau_x_d = tau_x
      tau_y = tau_y_s ; tau_y_d = tau_y
      tau_z = tau_z_s ; tau_z_d = tau_z
      momfluxb   = momfluxb_s
      momfluxb_d = mom_dev
      if (ltempeq) then
         thl_flux   = thl_flux_s ; thl_flux_d = thl_flux
         tfluxb     = tfluxb_s
         tfluxb_d   = tf_dev
         deallocate(thl_flux_s, tfluxb_s, tf_dev)
      end if

      ! updateDevice cleared the facet dirty flag on the way through. Set it
      ! again so the first real call copies the facet properties exactly as it
      ! would have without this test.
      lfacetprops_dirty = .true.

      deallocate(ekm_s, ekh_s, tau_x_s, tau_y_s, tau_z_s, momfluxb_s, mom_dev, back)

   contains

      !> Assert a fielddump-gated field was pulled exactly when it has a reader.
      subroutine check_gated(code, host, sentinel, name)
         implicit none
         character(len=2), intent(in) :: code
         real,             intent(in) :: host(:,:,:)
         real,             intent(in) :: sentinel
         character(len=*), intent(in) :: name

         if (fielddump_wants(code)) then
            if (any(host /= sentinel)) &
               call fail_cuda_selftest('the handover did not bring down '//name)
         else
            if (any(host /= 0.)) &
               call fail_cuda_selftest('the handover brought down '//name//' with no reader')
         end if

      end subroutine check_gated

   end subroutine test_post_poisson_handover


   !> Check the device branch of the checksim reductions against a host
   !! reference built from the same seed.
   !!
   !! The three reductions are the whole of what modchecksim computes, and
   !! under _GPU they read the device fields, so runmode 1014 - which drives
   !! the host branch - cannot reach them. This is their only coverage outside
   !! a full parity run.
   !!
   !! The reference is an explicit host loop over the same seed function, not a
   !! call to the host branch, which is not compiled in this build at all. What
   !! it pins:
   !!
   !!   - Every grid factor. dxhi, dzhi and dzfi vary with the index only on a
   !!     stretched grid, which no GPU case uses, so the seed varies in all
   !!     three directions instead and the reference multiplies each term by
   !!     the factor it is supposed to carry.
   !!
   !!   - The reduction box. The Courant and diffusion fields are poisoned
   !!     everywhere outside ib..ie, jb..je, kb..ke, so a kernel that reduces
   !!     over the halo returns the poison rather than a slightly wrong number.
   !!
   !!   - The forward reach in div_local, which is the one place a device loop
   !!     reads past the box on purpose. It is checked twice: once with the
   !!     poison confined to the cells div never reads, and once with a single
   !!     spike at ie+1, je+1 and ke+1 in turn.
   !!
   !!   - That diffnrgeom_d holds the cache initchecksim built. The reference
   !!     reads the host diffnrgeom, so a mirror that was never uploaded, or
   !!     uploaded transposed, gives the wrong answer rather than none.
   !!
   !! The device fields are saved and restored, because the self-tests run
   !! between initCUDA and the time loop and the run continues afterwards.
   subroutine test_checksim_reductions
      implicit none

      real, parameter :: dtm = 0.25, spike = 2., poison = 1.e6

      real, allocatable :: sa(:,:,:), sb(:,:,:), sc(:,:,:)   ! saved device state
      real, allocatable :: ta(:,:,:), tb(:,:,:), tc(:,:,:)   ! staging for the seed
      real, allocatable :: dxhi_save(:), dzhi_save(:), dzf_save(:), dzfi_save(:), dv_save(:)
      real, allocatable :: geom_save(:,:)
      real    :: dxi_save, dyi_save, dx_save, dy_save
      real    :: got, gotmax, gottot, ref, refmax, reftot
      integer :: i, j, k

      if (.not. allocated(um_d)) return
      if (.not. allocated(diffnrgeom)) call fail_cuda_selftest('checksim geometry cache not built')

      allocate(sa(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
      allocate(sb(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
      allocate(sc(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
      allocate(ta(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
      allocate(tb(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
      allocate(tc(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))

      ! ---- make every grid factor distinct ------------------------------------
      ! Every GPU parity case runs on a grid that is uniform in all three
      ! directions and has dx = dy. There dxhi(i) equals dyi and every one of
      ! dxhi, dzhi, dzfi, dzf and the geometry cache is constant, so a kernel
      ! that pairs the v term with the x spacing, or reads any of them at a
      ! fixed index, gives exactly the right answer. Substituting values that
      ! differ in every position is what makes the pairing and the indexing
      ! load-bearing here, in the absence of an anisotropic or stretched case.
      !
      ! dzfi is deliberately not the reciprocal of dzf, and dxi not that of dx,
      ! because div_local uses one of each pair for the gradient and the other
      ! for the cell volume.
      !
      ! The scalars reach the kernels as implicit firstprivate copies taken
      ! from the host at launch, so setting the host value is enough; the
      ! arrays need their device mirrors set too. Both are restored below.
      allocate(dxhi_save(lbound(dxhi,1):ubound(dxhi,1))) ; dxhi_save = dxhi
      allocate(dzhi_save(lbound(dzhi,1):ubound(dzhi,1))) ; dzhi_save = dzhi
      allocate(dzf_save (lbound(dzf ,1):ubound(dzf ,1))) ; dzf_save  = dzf
      allocate(dzfi_save(lbound(dzfi,1):ubound(dzfi,1))) ; dzfi_save = dzfi
      allocate(dv_save  (lbound(dvcell,1):ubound(dvcell,1))) ; dv_save   = dvcell
      allocate(geom_save(lbound(diffnrgeom,1):ubound(diffnrgeom,1), &
                         lbound(diffnrgeom,2):ubound(diffnrgeom,2)))
      geom_save = diffnrgeom
      dxi_save = dxi ; dyi_save = dyi ; dx_save = dx ; dy_save = dy

      ! Strictly increasing, and exact in binary: no two indices can collide,
      ! so reading any of these one index off is always visible. A cyclic
      ! pattern would not be - with a period of seven and a 64-wide domain,
      ! dxhi(ib) and dxhi(ie) come out equal and a kernel pinned to ib passes.
      do i = lbound(dxhi,1), ubound(dxhi,1)
         dxhi(i) = 0.375 + 0.0078125*real(i - lbound(dxhi,1))
      end do
      do k = lbound(dzhi,1), ubound(dzhi,1)
         dzhi(k) = 0.75 + 0.0078125*real(k - lbound(dzhi,1))
      end do
      do k = lbound(dzf,1), ubound(dzf,1)
         dzf(k) = 1.5 + 0.0078125*real(k - lbound(dzf,1))
      end do
      do k = lbound(dzfi,1), ubound(dzfi,1)
         dzfi(k) = 0.25 + 0.00390625*real(k - lbound(dzfi,1))
      end do
      do k = lbound(dvcell,1), ubound(dvcell,1)
         dvcell(k) = 2.5 + 0.0078125*real(k - lbound(dvcell,1))
      end do
      do k = lbound(diffnrgeom,2), ubound(diffnrgeom,2)
         do i = lbound(diffnrgeom,1), ubound(diffnrgeom,1)
            diffnrgeom(i,k) = 1.25 + 0.0078125*real( &
              (k - lbound(diffnrgeom,2))*(ubound(diffnrgeom,1) - lbound(diffnrgeom,1) + 1) &
              + (i - lbound(diffnrgeom,1)))
         end do
      end do
      dxi = 0.25 ; dyi = 0.5 ; dx = 3. ; dy = 5.

      dxhi_d = dxhi ; dzhi_d = dzhi ; dzf_d = dzf ; dzfi_d = dzfi
      diffnrgeom_d = diffnrgeom ; dvcell_d = dvcell

      ! ---- the Courant number -------------------------------------------------
      sa = um_d ; sb = vm_d ; sc = wm_d

      call seed_box(ta, 1, poison)
      call seed_box(tb, 2, poison)
      call seed_box(tc, 3, poison)
      um_d = ta ; vm_d = tb ; wm_d = tc

      ref = 0.
      do k = kb, ke
         do j = jb, je
            do i = ib, ie
               ref = max(ref, (ta(i,j,k)*dxhi(i) + tb(i,j,k)*dyi + tc(i,j,k)*dzhi(k))*dtm)
            end do
         end do
      end do

      call courant_local(dtm, got)
      if (.not. matches(got, ref, 1.)) call fail_cuda_selftest('checksim courant_local')

      um_d = sa ; vm_d = sb ; wm_d = sc

      ! ---- the diffusion number -----------------------------------------------
      sa = ekm_d ; sb = ekh_d

      ! Positive, as a diffusivity is, and different in ekm and ekh so that
      ! whichever term the kernel drops changes the answer.
      call seed_box(ta, 4, poison)
      call seed_box(tb, 5, poison)
      do k = kb, ke
         do j = jb, je
            do i = ib, ie
               ta(i,j,k) = abs(ta(i,j,k)) + 0.5
               tb(i,j,k) = abs(tb(i,j,k)) + 1.5
            end do
         end do
      end do
      ekm_d = ta ; ekh_d = tb

      ref = 0.
      do k = kb, ke
         do j = jb, je
            do i = ib, ie
               ref = max(ref, ta(i,j,k)*diffnrgeom(i,k)*dtm, &
                              tb(i,j,k)*diffnrgeom(i,k)*dtm)
            end do
         end do
      end do

      call diffnr_local(dtm, got)
      if (.not. matches(got, ref, 1.)) call fail_cuda_selftest('checksim diffnr_local')

      ! ekh alone, with ekm zeroed, so a kernel that keeps only ekm returns
      ! zero here rather than a number that happens to be close.
      ekm_d = 0.
      ref = 0.
      do k = kb, ke
         do j = jb, je
            do i = ib, ie
               ref = max(ref, tb(i,j,k)*diffnrgeom(i,k)*dtm)
            end do
         end do
      end do
      call diffnr_local(dtm, got)
      if (.not. matches(got, ref, 1.)) call fail_cuda_selftest('checksim diffnr_local ekh term')

      ekm_d = sa ; ekh_d = sb

      ! ---- the divergence -----------------------------------------------------
      sa = u0_d ; sb = v0_d ; sc = w0_d

      ! No poison this time: div_local reads i+1, j+1 and k+1 on purpose, so
      ! the seed covers the whole declared extent and the reference reads the
      ! same cells the kernel is entitled to.
      call seed_all(ta, 6)
      call seed_all(tb, 7)
      call seed_all(tc, 8)
      u0_d = ta ; v0_d = tb ; w0_d = tc

      call div_reference(ta, tb, tc, refmax, reftot)
      call div_local(gotmax, gottot)
      if (.not. matches(gotmax, refmax, 1.)) call fail_cuda_selftest('checksim div_local divmax')
      if (.not. sums_to(gottot, reftot)) call fail_cuda_selftest('checksim div_local divtot')

      ! Now poison only what div_local must never read - one cell short of the
      ! forward neighbours it must. Both results have to be unchanged.
      call poison_unread(ta, poison)
      call poison_unread(tb, poison)
      call poison_unread(tc, poison)
      u0_d = ta ; v0_d = tb ; w0_d = tc
      call div_local(gotmax, gottot)
      if (.not. matches(gotmax, refmax, 1.)) call fail_cuda_selftest('checksim div_local reads outside its box')
      if (.not. sums_to(gottot, reftot)) call fail_cuda_selftest('checksim div_local sums outside its box')

      ! The forward reach itself: one spike per direction, everything else
      ! zero, so the cell that sees it is the last interior one. A backward
      ! difference, or a loop stopping at ie-1, returns nothing here.
      if (.not. div_spike(1, spike*dxi     )) call fail_cuda_selftest('checksim div_local u0 at ie+1')
      if (.not. div_spike(2, spike*dyi     )) call fail_cuda_selftest('checksim div_local v0 at je+1')
      if (.not. div_spike(3, spike*dzfi(ke))) call fail_cuda_selftest('checksim div_local w0 at ke+1')

      u0_d = sa ; v0_d = sb ; w0_d = sc

      dxhi = dxhi_save ; dzhi = dzhi_save ; dzf = dzf_save ; dzfi = dzfi_save
      diffnrgeom = geom_save ; dvcell = dv_save
      dxi = dxi_save ; dyi = dyi_save ; dx = dx_save ; dy = dy_save
      dxhi_d = dxhi ; dzhi_d = dzhi ; dzf_d = dzf ; dzfi_d = dzfi
      diffnrgeom_d = diffnrgeom ; dvcell_d = dvcell

      deallocate(geom_save, dv_save, dzfi_save, dzf_save, dzhi_save, dxhi_save)
      deallocate(tc, tb, ta, sc, sb, sa)

   contains

      !> Deterministic, exactly representable seed. Varies in all three
      !! directions so no grid factor can be swapped for another unnoticed.
      real function seed_at(n, i, j, k)
         implicit none
         integer, intent(in) :: n, i, j, k

         seed_at = 0.125*real(modulo(i*7 + j*13 + k*29 + n*5, 17)) - 1.

      end function seed_at

      !> Seed the reduction box and poison everything outside it.
      subroutine seed_box(a, n, bad)
         implicit none
         real,    intent(out) :: a(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
         integer, intent(in)  :: n
         real,    intent(in)  :: bad
         integer :: i, j, k

         a = bad
         do k = kb, ke
            do j = jb, je
               do i = ib, ie
                  a(i,j,k) = seed_at(n, i, j, k)
               end do
            end do
         end do

      end subroutine seed_box

      !> Seed the whole declared extent, halos included.
      subroutine seed_all(a, n)
         implicit none
         real,    intent(out) :: a(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
         integer, intent(in)  :: n
         integer :: i, j, k

         do k = kb-kh, ke+kh
            do j = jb-jh, je+jh
               do i = ib-ih, ie+ih
                  a(i,j,k) = seed_at(n, i, j, k)
               end do
            end do
         end do

      end subroutine seed_all

      !> Overwrite every cell div_local has no business reading, leaving the
      !! forward neighbours at ie+1, je+1 and ke+1 alone.
      subroutine poison_unread(a, bad)
         implicit none
         real, intent(inout) :: a(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
         real, intent(in)    :: bad
         integer :: i, j, k

         do k = kb-kh, ke+kh
            do j = jb-jh, je+jh
               do i = ib-ih, ie+ih
                  if (i < ib .or. i > ie+1 .or. &
                      j < jb .or. j > je+1 .or. &
                      k < kb .or. k > ke+1) a(i,j,k) = bad
               end do
            end do
         end do

      end subroutine poison_unread

      !> The divergence reduction, written out here rather than called, so a
      !! transcription error in the kernel cannot be reproduced by its own
      !! reference.
      subroutine div_reference(au, av, aw, dmax, dtot)
         implicit none
         real, intent(in)  :: au(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
         real, intent(in)  :: av(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
         real, intent(in)  :: aw(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
         real, intent(out) :: dmax, dtot
         real    :: d
         integer :: i, j, k

         dmax = 0.
         dtot = 0.
         do k = kb, ke
            do j = jb, je
               do i = ib, ie
                  d = (au(i+1,j,k) - au(i,j,k))*dxi + &
                      (av(i,j+1,k) - av(i,j,k))*dyi + &
                      (aw(i,j,k+1) - aw(i,j,k))*dzfi(k)
                  dmax = max(dmax, abs(d))
                  dtot = dtot + d*dvcell(k)
               end do
            end do
         end do

      end subroutine div_reference

      !> One spike just past the top of the box in the direction given, with
      !! every field otherwise zero. want is the divergence the single cell
      !! that sees it must report.
      logical function div_spike(which, want)
         implicit none
         integer, intent(in) :: which
         real,    intent(in) :: want

         u0_d = 0. ; v0_d = 0. ; w0_d = 0.
         ta = 0.
         select case (which)
         case (1) ; ta(ie+1,je,ke) = spike ; u0_d = ta
         case (2) ; ta(ie,je+1,ke) = spike ; v0_d = ta
         case (3) ; ta(ie,je,ke+1) = spike ; w0_d = ta
         end select

         call div_local(gotmax, gottot)
         div_spike = matches(gotmax, abs(want), 1.) .and. &
                     matches(gottot, want*dvcell(ke), 1.)

      end function div_spike

      !> Equal to within a few ulps of the larger of want and floor.
      logical function matches(got, want, floor)
         implicit none
         real, intent(in) :: got, want, floor

         matches = abs(got - want) <= 64.*epsilon(1.)*max(floor, abs(want))
         if (.not. matches) then
            write(*,'(A,I0,A,ES22.14,A,ES22.14)') 'checksim self-test rank ', myid, &
                 ': got ', got, ' want ', want
         end if

      end function matches

      !> Equal to within the rounding a sum over the reduction box can accrue.
      !! The device adds the terms in scheduler order and the host reference
      !! sequentially, so the two agree to order n*epsilon, not to the last bit.
      logical function sums_to(got, want)
         implicit none
         real, intent(in) :: got, want
         real :: terms

         terms   = real(ie-ib+1)*real(je-jb+1)*real(ke-kb+1)
         sums_to = abs(got - want) <= 8.*terms*epsilon(1.)*max(1., abs(want))
         if (.not. sums_to) then
            write(*,'(A,I0,A,ES22.14,A,ES22.14)') 'checksim self-test rank ', myid, &
                 ': sum got ', got, ' want ', want
         end if

      end function sums_to

   end subroutine test_checksim_reductions

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
