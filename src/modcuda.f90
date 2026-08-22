module modcuda
#if defined(_GPU)        
   use cudafor
#if defined(UDALES_DEBUG)
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#endif
   use modglobal,      only: itot, ib, ie, jb, je, kb, ke, ih, jh, kh, ihc, jhc, khc, &
                             dx2, dxi, dx2i, dxi5, dxiq, dy2, dyi, dy2i, dyi5, dyiq, dxf, dxhi, &
                             dzf, dzf2, dzfi, dzfi5, dzfiq, dzh, dzhi, dzh2i, dzhiq, &
                             dxdydzfi, dxdydzhi, dvcell, &
                             dzfc, dzfci, dzhci, dxfc, dxfci, dxhci, delta, &
                             ltempeq, lmoist, nsv, lchem, lles, lbuoyancy, ltrees, lscasrc, lscasrcl, &
                             BCxm, BCxm_profile, BCxm_driver, BCym, BCym_profile, &
                             lnudge, lnudgevel, libm, nfcts, &
                             linoutflow, luoutflowr, lvoutflowr, luvolflowr, lvvolflowr, &
                             iadv_sv, iadv_thl, iadv_kappa, iadv_upw, &
                             xh, &
                             eps1, numol, prandtlmoli, prandtlturb, grav, fkar2
   use modfields,      only: u0, v0, w0, pres0, e120, thl0, thl0c, qt0, sv0, &
                             up, vp, wp, e12p, thlp, thlpc, qtp, svp, &
                             um, vm, wm, e12m, thlm, qtm, svm, &
                             tau_x, tau_y, tau_z, thl_flux, momfluxb, tfluxb, &
                             u0av, v0av, thl0av, qt0av, sv0av, dthvdz, ug, vg, whls, tsc, &
                             dpdxl, dpdyl, thv0h, thvh, thlpcar, &
                             dudxls, dudyls, dvdxls, dvdyls, dthldxls, dthldyls, dqtdxls, dqtdyls, dqtdtls, &
                             uprof, vprof, thlprof, qtprof, svprof, &
                             IIc, IIu, IIv, scalar_source_tendency
   use modsubgriddata, only: lsmagorinsky, lvreman, loneeqn, ldelta, lbuoycorr, &
                             ekm, ekh, &
                             sbshr, sbbuo, sbdiss, zlt, damp, csz, &
                             cn, cm, ch1, ch2, ce1, ce2, dampmin, prandtli, c_vreman
   use modsurfdata,    only: thvs
   use initfac,        only: facT, fachf, facef, fachfi, facefi, facqsat, fachurel, facf
   use modinletdata,   only: u0driver
   use decomp_2d,      only: zstart
   implicit none
   save

   type(dim3) :: griddim, blockdim

   integer, device :: ib_d, ie_d, jb_d, je_d, kb_d, ke_d, ih_d, jh_d, kh_d
   
   logical, device :: ltempeq_d, lles_d, lsmagorinsky_d, lvreman_d, loneeqn_d, ldelta_d, lbuoyancy_d, lbuoycorr_d

   real,    device :: dx2_d, dxi_d, dx2i_d, dxi5_d, dxiq_d, dy2_d, dyi_d, dy2i_d, dyi5_d, dyiq_d, &
                      eps1_d, numol_d, prandtlmoli_d, prandtlturb_d, prandtli_d, grav_d, dampmin_d, c_vreman_d, fkar2_d, &
                      cn_d, cm_d, ch1_d, ch2_d, ce1_d, ce2_d, &
                      thvs_d, bcTfluxA_d

   integer, device, dimension(3) :: zstart_d

   real, device, allocatable :: dzf_d(:), dzf2_d(:), dzfi_d(:), dzfi5_d(:), dzfiq_d(:), dzh_d(:), dzhi_d(:), dzh2i_d(:), dzhiq_d(:), &
                                dzfc_d(:), dzfci_d(:), dzhci_d(:), dxfc_d(:), dxfci_d(:), dxhci_d(:), &
                                dxf_d(:), dxhi_d(:), dxdydzfi_d(:), dxdydzhi_d(:), dvcell_d(:), &
                                xh_d(:), &
                                u0av_d(:), v0av_d(:), ug_d(:), vg_d(:), whls_d(:), thl0av_d(:), qt0av_d(:), tsc_d(:), &
                                dpdxl_d(:), dpdyl_d(:), thvh_d(:), thlpcar_d(:), &
                                dudxls_d(:), dudyls_d(:), dvdxls_d(:), dvdyls_d(:), dthldxls_d(:), dthldyls_d(:), dqtdxls_d(:), dqtdyls_d(:), dqtdtls_d(:), &
                                uprof_d(:), vprof_d(:), thlprof_d(:), qtprof_d(:)

   real, device, allocatable :: delta_d(:, :), csz_d(:,:), sv0av_d(:,:), svprof_d(:,:), u0driver_d(:,:)

   real, device, allocatable :: u0_d(:,:,:), v0_d(:,:,:), w0_d(:,:,:), pres0_d(:,:,:), e120_d(:,:,:), thl0_d(:,:,:), thl0c_d(:,:,:), qt0_d(:,:,:), sv0_d(:,:,:,:)
   real, device, allocatable :: up_d(:,:,:), vp_d(:,:,:), wp_d(:,:,:), e12p_d(:,:,:), thlp_d(:,:,:), thlpc_d(:,:,:), qtp_d(:,:,:), svp_d(:,:,:,:)
   real, device, allocatable :: um_d(:,:,:), vm_d(:,:,:), wm_d(:,:,:), e12m_d(:,:,:), thlm_d(:,:,:), qtm_d(:,:,:), svm_d(:,:,:,:)
   real, device, allocatable :: p_d(:,:,:), pup_d(:,:,:), pvp_d(:,:,:), pwp_d(:,:,:)
   real, device, allocatable :: tau_x_d(:,:,:), tau_y_d(:,:,:), tau_z_d(:,:,:), thl_flux_d(:,:,:)
   real, device, allocatable :: momfluxb_d(:,:,:), tfluxb_d(:,:,:)
   real, device, allocatable :: dthvdz_d(:,:,:)
   real, device, allocatable :: ekm_d(:,:,:), ekh_d(:,:,:), sbshr_d(:,:,:), sbbuo_d(:,:,:), sbdiss_d(:,:,:), zlt_d(:,:,:), damp_d(:,:,:)
   real, device, allocatable :: thv0h_d(:,:,:)
   real, device, allocatable :: dumu_d(:,:,:), duml_d(:,:,:)
   real, device, allocatable :: dummyNO_d(:,:,:), dummyNO2_d(:,:,:), dummyO3_d(:,:,:)

   integer, device, allocatable :: IIc_d(:,:,:), IIu_d(:,:,:), IIv_d(:,:,:)

   real, device, allocatable :: scalar_source_tendency_d(:,:,:,:)

   ! IBM facet quantities that cross the bus inside the time loop. They are
   ! declared here, rather than in modibm, so that the transfers can live in the
   ! four update routines below: modibm already uses modcuda, so the reverse
   ! dependency would be circular. The IBM geometry that never changes stays in
   ! modibm and is transferred once by init_ibm_device.
   real, device, allocatable :: facT1_d(:)                  ! outdoor surface temperature
   ! Green-roof properties. The energy balance rewrites all of these whenever it
   ! solves, so they are refreshed alongside facT rather than mirrored once.
   real, device, allocatable :: facqsat_d(:), fachurel_d(:), facf_d(:,:)
   real, device, allocatable :: fachf_d(:), facef_d(:)      ! energy balance accumulators
   ! ... and their per-rank time integrals, which stay on the device between
   ! energy balances so that nothing has to cross the bus per step.
   real, device, allocatable :: fachfi_d(:), facefi_d(:)
   real, device, allocatable :: fac_tau_d(:)                ! momentum, one direction at a time
   real, device, allocatable :: fac_htc_d(:), fac_cth_d(:), fac_pres_d(:), fac_pres2_d(:)

   ! Pinned staging for the facet transfers. Automatic arrays cannot be pinned,
   ! so these are held at module level to keep every loop-time copy pinned.
   real, allocatable, pinned :: fac_stage(:), fachf_stage(:), facef_stage(:)

   ! The mass correction reduces one column of per-level sums per call. MPI is
   ! not CUDA-aware, so that column has to reach the host before the
   ! all-reduce: col_d holds it on the device and col_stage is the pinned
   ! landing site. Automatic arrays cannot be pinned, hence module scope.
   real, device, allocatable :: col_d(:)
   real, allocatable, pinned :: col_stage(:)
   logical :: lfacets_on_device = .false.

   contains
      subroutine initCUDA
         implicit none
         integer :: threadnumx, threadnumy, threadnumz
         integer :: blocknumx, blocknumy, blocknumz
         type(cudaDeviceProp) :: props
         integer :: deviceID, SMcount
         
         call checkCUDA( cudaGetDevice( deviceID ), 'cudaGetDevice_in_modcuda' )
         call checkCUDA( cudaGetDeviceProperties(props, deviceId), 'cudaGetDeviceProperties_in_modcuda' )
         
         !! SMcount = props%multiProcessorCount
         !! write(*,*) props%warpsize, props%multiProcessorCount, props%maxThreadsPerBlock, props%maxThreadsPerMultiProcessor

         threadnumx = props%warpsize/2
         threadnumy = props%warpsize/8
         threadnumz = 2
         if (threadnumx*threadnumy*threadnumz > props%maxThreadsPerBlock) then
            write(*,*) "Incorrect block dimension configuration."
            stop 1
         end if

         blocknumx = min( max( floor(real( (ie - ib + 1)/threadnumx ) ), 1 ), props%maxGridSize(1) )
         blocknumy = min( max( floor(real( (je - jb + 1)/threadnumy ) ), 1 ), props%maxGridSize(2) )
         blocknumz = 32 ! min( max( floor(real( (ke - kb + 1)/threadnumz ) ), 1 ), props%maxGridSize(3) )

         write(*,*) "CUDA block dimension: (", threadnumx, ",", threadnumy, ",", threadnumz, ")"
         write(*,*) "CUDA grid dimension:  (", blocknumx, ",", blocknumy, ",", blocknumz, ")"

         blockdim  = dim3(threadnumx,threadnumy,threadnumz)
         griddim   = dim3(blocknumx,blocknumy,blocknumz)

         ib_d = ib
         ie_d = ie
         ih_d = ih
         jb_d = jb
         je_d = je
         jh_d = jh
         kb_d = kb
         ke_d = ke
         kh_d = kh

         dx2_d  = dx2
         dxi_d  = dxi
         dxi5_d = dxi5
         dxiq_d = dxiq
         dx2i_d = dx2i

         dy2_d  = dy2
         dyi_d  = dyi
         dyi5_d = dyi5
         dyiq_d = dyiq
         dy2i_d = dy2i

         allocate (dxf_d(ib-ih:itot+ih))
         allocate (dxhi_d(ib:itot+ih))
         dxf_d  = dxf
         dxhi_d = dxhi

         allocate (dzf_d(kb - kh:ke + kh))
         allocate (dzf2_d(kb - kh:ke + kh))
         allocate (dzfi_d(kb - kh:ke + kh))
         allocate (dzfi5_d(kb - kh:ke + kh))
         allocate (dzfiq_d(kb - kh:ke + kh))
         allocate (dzh_d(kb:ke + kh))
         allocate (dzhi_d(kb:ke + kh))
         allocate (dzh2i_d(kb:ke + kh))
         allocate (dxdydzfi_d(kb - kh:ke + kh))
         allocate (dxdydzhi_d(kb:ke + kh))
         allocate (dvcell_d(kb - kh:ke + kh))
         allocate (dzhiq_d(kb:ke + kh))
         dzf_d   = dzf
         dzf2_d  = dzf2
         dzfi_d  = dzfi
         dzfi5_d = dzfi5
         dzfiq_d = dzfiq
         dzh_d   = dzh
         dzhi_d  = dzhi
         dxdydzfi_d = dxdydzfi
         dxdydzhi_d = dxdydzhi
         dvcell_d   = dvcell
         dzh2i_d = dzh2i
         dzhiq_d = dzhiq

         allocate(delta_d(ib-ih:itot+ih, kb:ke + kh))
         delta_d = delta

         allocate(u0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         allocate(v0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         allocate(w0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         allocate(pres0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         
         allocate(um_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         allocate(vm_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         allocate(wm_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))

         allocate(up_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
         allocate(vp_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
         allocate(wp_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))

         ! Seeded once here rather than every stage. bottom rewrites the three
         ! stresses whole before anything reads them, so their seed only has to
         ! be defined - but momfluxb is different: the wall functions add into
         ! it and nothing ever clears it, on either the host or the device, so
         ! it carries a running total for the whole run. Uploading the host
         ! copy every stage used to be what carried that total across; now it
         ! simply stays on the device, which is what the host does too.
         allocate(tau_x_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh));    tau_x_d    = tau_x
         allocate(tau_y_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh));    tau_y_d    = tau_y
         allocate(tau_z_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh));    tau_z_d    = tau_z
         allocate(momfluxb_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)); momfluxb_d = momfluxb

         allocate(u0av_d(kb:ke+kh))
         allocate(v0av_d(kb:ke+kh))

         allocate(ug_d(kb:ke+kh))
         allocate(vg_d(kb:ke+kh))
         allocate(whls_d(kb:ke+kh))
         ug_d = ug
         vg_d = vg
         whls_d = whls

         allocate(tsc_d(kb:ke+kh))
         tsc_d = tsc

         allocate(dpdxl_d(kb:ke+kh))
         allocate(dpdyl_d(kb:ke+kh))
         dpdxl_d = dpdxl
         dpdyl_d = dpdyl

         allocate(dudxls_d(kb:ke+kh))
         allocate(dudyls_d(kb:ke+kh))
         allocate(dvdxls_d(kb:ke+kh))
         allocate(dvdyls_d(kb:ke+kh))
         dudxls_d = dudxls
         dudyls_d = dudyls
         dvdxls_d = dvdxls
         dvdyls_d = dvdyls

         allocate(uprof_d(kb:ke+kh))
         allocate(vprof_d(kb:ke+kh))
         allocate(u0driver_d(jb-jh:je+jh,kb-kh:ke+kh))

         if (loneeqn) then
            allocate(e120_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
            allocate(e12m_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
            allocate(e12p_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
            allocate(sbshr_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
            allocate(sbbuo_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
            allocate(sbdiss_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
            allocate(zlt_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
         end if

         if (ltempeq) then
            allocate(thl0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
            allocate(thlm_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
            allocate(thlp_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
            if (iadv_thl == iadv_kappa) then
               allocate(thl0c_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb-khc:ke+khc))
               allocate(thlpc_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
            end if
            allocate(thl_flux_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)); thl_flux_d = thl_flux
            allocate(thv0h_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
            allocate(thvh_d(kb:ke+kh))
            allocate(thl0av_d(kb:ke+kh))
            allocate(thlprof_d(kb:ke+kh))

            allocate(thlpcar_d(kb:ke+kh))
            allocate(dthldxls_d(kb:ke+kh))
            allocate(dthldyls_d(kb:ke+kh))
            thlpcar_d = thlpcar
            dthldxls_d = dthldxls
            dthldyls_d = dthldyls

            ! An accumulator like momfluxb - see the comment above.
            allocate(tfluxb_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)); tfluxb_d = tfluxb
         end if

         if (lmoist) then
            allocate(qt0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
            allocate(qtm_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
            allocate(qtp_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
            allocate(qt0av_d(kb:ke+kh))
            allocate(qtprof_d(kb:ke+kh))

            allocate(dqtdxls_d(kb:ke+kh))
            allocate(dqtdyls_d(kb:ke+kh))
            allocate(dqtdtls_d(kb:ke+kh))
            dqtdxls_d = dqtdxls
            dqtdyls_d = dqtdyls
            dqtdtls_d = dqtdtls
         end if

         if (nsv>0) then
            allocate(sv0_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb-khc:ke+khc,nsv))
            allocate(svm_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb-khc:ke+khc,nsv))
            allocate(svp_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc,nsv))
            allocate(sv0av_d(kb:ke+khc,nsv))
            allocate(svprof_d(kb:ke+kh,nsv))
            if (nsv==3 .and. lchem) then
               allocate(dummyNO_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
               allocate(dummyNO2_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
               allocate(dummyO3_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
               allocate(IIc_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
               IIc_d = IIc
            end if
            if (lscasrc .or. lscasrcl) then
               allocate(scalar_source_tendency_d(ib:ie,jb:je,kb:ke,nsv))
               scalar_source_tendency_d = scalar_source_tendency
            end if
         end if

         allocate(ekm_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         allocate(ekh_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         allocate(dthvdz_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))

         if (any(iadv_sv(1:nsv) == iadv_kappa) .or. (iadv_thl == iadv_kappa)) then
            allocate(dumu_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
            allocate(duml_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
            allocate (dxhci_d(ib - 1:itot + ihc))
            allocate (dxfc_d(ib - ihc:itot + ihc))
            allocate (dzhci_d(kb - 1:ke + khc))
            allocate (dzfc_d(kb - khc:ke + khc))
            dxhci_d = dxhci
            dxfc_d  = dxfc
            dzhci_d = dzhci
            dzfc_d  = dzfc
         end if
         if (any(iadv_sv(1:nsv) == iadv_kappa) .or. any(iadv_sv(1:nsv) == iadv_upw) .or. (iadv_thl == iadv_kappa)) then
            allocate (dxfci_d(ib - ihc:itot + ihc))
            allocate (dzfci_d(kb - khc:ke + khc))
            dxfci_d = dxfci
            dzfci_d = dzfci
         end if
         
         eps1_d   = eps1
         zstart_d = zstart
         ! Only xh is mirrored: the wall-function reconstruction used to need
         ! the full set, but its trilinear offsets are now precomputed per
         ! section in modibm, so the grid never enters those kernels.
         allocate (xh_d(ib:itot+ih))
         xh_d = xh

         ltempeq_d      = ltempeq
         lles_d         = lles
         lsmagorinsky_d = lsmagorinsky
         lvreman_d      = lvreman
         loneeqn_d      = loneeqn
         ldelta_d       = ldelta
         lbuoyancy_d    = lbuoyancy
         lbuoycorr_d    = lbuoycorr

         numol_d       = numol
         prandtlmoli_d = prandtlmoli
         prandtlturb_d = prandtlturb
         grav_d        = grav
         thvs_d        = thvs
         cm_d          = cm
         cn_d          = cn
         ch1_d         = ch1
         ch2_d         = ch2
         ce1_d         = ce1
         ce2_d         = ce2
         prandtli_d    = prandtli
         c_vreman_d    = c_vreman
         fkar2_d       = fkar2

         if (lsmagorinsky .or. loneeqn) then
            allocate(damp_d(ib:ie,jb:je,kb:ke))
            dampmin_d = dampmin
         end if

         if (lsmagorinsky) then
            allocate(csz_d(ib-ih:ie+ih,kb:ke+kh))
            csz_d = csz
         end if

         if (libm .and. nfcts > 0) then
            allocate(fac_stage(1:nfcts))
            allocate(fachf_stage(0:nfcts))
            allocate(facef_stage(1:nfcts))
         end if

         ! Only the flow-rate controllers masscorr can actually run need the
         ! masks and the staging column, so an unforced run pays for neither.
         if ((.not. linoutflow) .and. &
             (luoutflowr .or. luvolflowr .or. lvoutflowr .or. lvvolflowr)) then
            allocate(col_d(kb:ke+kh))
            allocate(col_stage(kb:ke+kh))
            if (luoutflowr .or. luvolflowr) then
               allocate(IIu_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
               IIu_d = IIu
            end if
            if (lvoutflowr .or. lvvolflowr) then
               allocate(IIv_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
               IIv_d = IIv
            end if
         end if

      end subroutine initCUDA

      subroutine exitCUDA
         implicit none
         deallocate(dxf_d, dxhi_d, dzf_d, dzf2_d, dzfi_d, dzfi5_d, dzfiq_d, dzh_d, dzhi_d, dzh2i_d, dzhiq_d, delta_d)
         deallocate(dxdydzfi_d, dxdydzhi_d, dvcell_d)
         deallocate(u0_d, v0_d, w0_d, pres0_d, um_d, vm_d, wm_d, up_d, vp_d, wp_d)
         deallocate(tau_x_d, tau_y_d, tau_z_d, momfluxb_d)
         deallocate(u0av_d, v0av_d, ug_d, vg_d, whls_d, tsc_d)
         deallocate(dpdxl_d, dpdyl_d, dudxls_d, dudyls_d, dvdxls_d, dvdyls_d)
         deallocate(uprof_d, vprof_d, u0driver_d)
         if (loneeqn) deallocate(e120_d, e12m_d, e12p_d, sbshr_d, sbbuo_d, sbdiss_d, zlt_d)
         if (ltempeq) then
            deallocate(thl0_d, thlm_d, thlp_d, thl_flux_d)
            if (iadv_thl == iadv_kappa) deallocate(thl0c_d, thlpc_d)
            deallocate(thv0h_d, thvh_d, thl0av_d, thlprof_d, thlpcar_d)
            deallocate(dthldxls_d, dthldyls_d)
            deallocate(tfluxb_d)
         end if
         if (lmoist) deallocate(qt0_d, qtm_d, qtp_d, qt0av_d, qtprof_d, dqtdxls_d, dqtdyls_d, dqtdtls_d)
         if (nsv>0) then
            deallocate(sv0_d, svm_d, svp_d, sv0av_d, svprof_d)
            if (lscasrc .or. lscasrcl) deallocate(scalar_source_tendency_d)
            if (nsv==3 .and. lchem) deallocate(dummyNO_d, dummyNO2_d, dummyO3_d, IIc_d)
         end if
         if (any(iadv_sv(1:nsv) == iadv_kappa) .or. (iadv_thl == iadv_kappa)) then
            deallocate(dumu_d, duml_d, dxhci_d, dxfc_d, dzhci_d, dzfc_d)
         end if
         if (any(iadv_sv(1:nsv) == iadv_kappa) .or. any(iadv_sv(1:nsv) == iadv_upw) .or. (iadv_thl == iadv_kappa)) then
            deallocate(dxfci_d, dzfci_d)
         end if
         deallocate(xh_d)
         deallocate(ekm_d, ekh_d)
         if (lsmagorinsky .or. loneeqn) then
            deallocate(damp_d)
         end if
         if (lsmagorinsky) deallocate(csz_d)
         deallocate(dthvdz_d)
         if (allocated(col_d)) deallocate(col_d, col_stage)
         if (allocated(IIu_d)) deallocate(IIu_d)
         if (allocated(IIv_d)) deallocate(IIv_d)
      end subroutine exitCUDA

      !> Refresh the facet properties the IBM wall functions read.
      !!
      !! facT, facqsat, fachurel and facf are rewritten by modEB, and only by
      !! modEB, on the steps where the energy balance fires - one step in
      !! dtEB/dt, which is a few hundred at a typical dtEB of 10 s. This used
      !! to upload all four on every step that reached rk3step == 1, so on the
      !! order of 99% of the traffic re-sent values the device already held.
      !!
      !! lfacetprops_dirty is how modEB says which steps are the real ones. It
      !! starts set, so the initial-condition values still make the crossing
      !! whether or not lEB is on, and the flag is cleared here rather than at
      !! the call site so that no caller can consume it without doing the copy.
      subroutine updateFacetPropsDevice
         use modglobal, only : nfcts
         use initfac,   only : lfacetprops_dirty
         implicit none

         if (lfacets_on_device .and. .not. lfacetprops_dirty) return

         if (allocated(facT1_d))    facT1_d    = facT(0:nfcts,1)
         if (allocated(facqsat_d))  facqsat_d  = facqsat(0:nfcts)
         if (allocated(fachurel_d)) fachurel_d = fachurel(0:nfcts)
         if (allocated(facf_d))     facf_d     = facf(0:nfcts,1:5)

         lfacets_on_device = .true.
         lfacetprops_dirty = .false.
      end subroutine updateFacetPropsDevice

      !> Integrate the energy balance accumulators on the device. Nothing copies.
      !!
      !! The wall functions add into fachf_d and facef_d on every Runge-Kutta
      !! stage; only the third is kept, because that is the one intqH used to
      !! reduce and time-integrate before clearing its host arrays regardless
      !! of the stage. So the third stage is integrated into fachfi_d and every
      !! stage is cleared.
      !!
      !! Clearing on every stage is not optional. Let fachf_d run on across the
      !! stages and the third one carries the sum of all three: the facet heat
      !! flux three times too large. That bug shipped once and a manual CPU/GPU
      !! comparison found it.
      !!
      !! The point of integrating here rather than on the host is that this
      !! whole routine is now free of traffic. The time integral is per-rank -
      !! summing over steps and summing over ranks commute, and dt is the same
      !! on every rank - so it can accumulate on the device across the hundreds
      !! of steps between energy balances, and only the total needs to come
      !! down. See updateFacIntegralsHost.
      subroutine integrateFacFluxDevice
         use modglobal, only : nfcts, lEB, rk3step, dt
         implicit none
         integer :: n

         if (.not. lEB) return
         if (.not. allocated(fachf_d)) return

         if (rk3step == 3) then
            !$acc parallel loop default(present)
            do n = 1, nfcts
               fachfi_d(n) = fachfi_d(n) + dt*fachf_d(n)
               facefi_d(n) = facefi_d(n) + dt*facef_d(n)
            end do
            !$acc end parallel loop
         end if

         fachf_d = 0.
         facef_d = 0.
      end subroutine integrateFacFluxDevice

      !> Bring the per-rank facet flux integrals down, for the energy balance.
      !!
      !! The one crossing the energy balance needs, and it happens only on the
      !! steps where the balance actually fires - one in dtEB/dt, a few hundred
      !! at a typical dtEB. Call it under modEB::eb_will_run so the decision
      !! here and the decision inside EB cannot drift apart.
      !!
      !! The device copies are cleared as they are read: nothing looks at them
      !! again until the next accumulation, and the host side is assigned
      !! rather than added to, so the interval starts clean on both sides.
      subroutine updateFacIntegralsHost
         use modglobal, only : nfcts, lEB
         implicit none

         if (.not. lEB) return

         if (allocated(fachfi_d) .and. allocated(fachfi)) then
            fachf_stage = fachfi_d
            fachfi(0:nfcts) = fachf_stage
            fachfi_d = 0.
         end if
         if (allocated(facefi_d) .and. allocated(facefi)) then
            facef_stage = facefi_d
            facefi(1:nfcts) = facef_stage
            facefi_d = 0.
         end if
      end subroutine updateFacIntegralsHost

      subroutine updateDevice
         implicit none
         integer :: n

         call updateFacetPropsDevice

         u0_d = u0
         v0_d = v0
         w0_d = w0
         pres0_d = pres0

         ! The previous-step fields. masscorr reduces um and vm, and ibmnorm
         ! pins all six inside the solid, both before updateHost hands the
         ! stage over - so the mirrors have to be current by this point. The
         ! host is the last writer: boundary applies the lateral and top
         ! conditions at the end of the previous step, for a profile or driver
         ! inlet at the interior plane i = ib and not only in the halo.
         um_d = um
         vm_d = vm
         wm_d = wm
         if (ltempeq) thlm_d = thlm
         if (lmoist)  qtm_d  = qtm
         if (nsv > 0) svm_d  = svm

         call initfield<<<griddim,blockdim>>>(up_d, 0., ih, jh, kh)
         call checkCUDA( cudaGetLastError(), 'initfield up_d' )

         call initfield<<<griddim,blockdim>>>(vp_d, 0., ih, jh, kh)
         call checkCUDA( cudaGetLastError(), 'initfield vp_d' )

         call initfield<<<griddim,blockdim>>>(wp_d, 0., ih, jh, kh)
         call checkCUDA( cudaGetLastError(), 'initfield wp_d' )

         u0av_d = u0av
         v0av_d = v0av

         ! The inlet profiles and the driver plane. updateDevicePriorPoiss
         ! used to be the only place that copied these for the BC paths: the
         ! nudging guard on its own leaves them behind whenever a profile or
         ! driver inlet is used without velocity nudging.
         !
         ! Both host writers run earlier in this same iteration: timedep fills
         ! uprof and vprof at the top of the loop, and drivergen fills
         ! u0driver from inside boundary at the end of the previous one.
         ! Nothing rewrites them between here and bcpup, which is the only
         ! reader, so the value that lands is the one that used to.
         if ((lnudge .and. lnudgevel) .or. BCxm == BCxm_profile) uprof_d = uprof
         if ((lnudge .and. lnudgevel) .or. BCym == BCym_profile) vprof_d = vprof
         if (BCxm == BCxm_driver) u0driver_d = u0driver

         if (loneeqn) then
            e12m_d = e12m
            e120_d = e120
            call initfield<<<griddim,blockdim>>>(e12p_d, 0., ih, jh, kh)
            call checkCUDA( cudaGetLastError(), 'initfield e12p_d' )
         end if

         if (ltempeq) then
            thl0_d = thl0
            call initfield<<<griddim,blockdim>>>(thlp_d, 0., ih, jh, kh)
            call checkCUDA( cudaGetLastError(), 'initfield thlp_d' )

            if (iadv_thl == iadv_kappa) then
               thl0c_d = thl0c
               call initfield<<<griddim,blockdim>>>(thlpc_d, 0., ihc, jhc, khc)
               call checkCUDA( cudaGetLastError(), 'initfield thlpc_d' )
            end if

            thv0h_d    = thv0h
            thvh_d     = thvh
            thl0av_d   = thl0av
            if (lnudge) thlprof_d = thlprof
            if (ltrees .and. lmoist) then
               thlpcar_d = thlpcar
            end if
         end if

         if (lmoist) then
            qt0_d = qt0
            call initfield<<<griddim,blockdim>>>(qtp_d, 0., ih, jh, kh)
            call checkCUDA( cudaGetLastError(), 'initfield qtp_d' )
            qt0av_d = qt0av
            if (lnudge) qtprof_d = qtprof
         end if

         if (nsv>0) then
            sv0_d = sv0
            do n = 1, nsv
               call initfield<<<griddim,blockdim>>>(svp_d(:, :, :, n), 0., ihc, jhc, khc)
               call checkCUDA( cudaGetLastError(), 'initfield svp_d' )
            end do
            sv0av_d = sv0av
            if (lnudge) svprof_d = svprof
         end if

         dthvdz_d = dthvdz
      end subroutine updateDevice

      subroutine updateHostAfterPoiss
         implicit none
         u0 = u0_d
         v0 = v0_d
         w0 = w0_d
         um = um_d
         vm = vm_d
         wm = wm_d
         pres0 = pres0_d
         if (ltempeq) then
            thl0 = thl0_d
            thlm = thlm_d
            if (iadv_thl == iadv_kappa) thl0c = thl0c_d
         end if
         if (lmoist) then
            qt0 = qt0_d
            qtm = qtm_d
         end if
         if (nsv>0) then
            sv0 = sv0_d
            svm = svm_d
#if defined(UDALES_DEBUG)
            if (.not. all(ieee_is_finite(sv0(ib:ie, jb:je, kb:ke, 1:nsv)))) then
               write(*,*) 'Non-finite scalar value detected after the GPU pressure step.'
               error stop 1
            end if
#endif
         end if
         if (loneeqn) then
            e120 = e120_d
            e12m = e12m_d
         end if

         ! What updateHost used to bring down half a stage earlier. It ran
         ! before poisson because vegetation_forcing needed the fields there;
         ! once that moved onto the device, everything left in that routine
         ! either had no host reader at all - the tendencies, momfluxb, tfluxb -
         ! or had one that runs past this point:
         !
         !   ekm, ekh           statsdump, boundary's fluxtop,
         !                      writerestartfiles
         !   tau_x/y/z          fielddump, when named in fieldvars
         !   thl_flux           the same
         !
         ! Nothing writes any of them on the device between there and here -
         ! subgrid and bottom produce them before, poisson and tstep_integrate
         ! do not touch them, and halos_device exchanges only the prognostic
         ! fields - so the values are the ones updateHost used to hand over.
         ekm = ekm_d
         ekh = ekh_d
         tau_x = tau_x_d
         tau_y = tau_y_d
         tau_z = tau_z_d
         if (ltempeq) thl_flux = thl_flux_d
      end subroutine updateHostAfterPoiss

      !> Rank-local fluid-cell column sums of a device field, for avexy_ibm.
      !!
      !! Forms the same averl the host branch of avexy_ibm forms, on the
      !! device, and then hands it to avexy_ibm_finish. Only the column crosses
      !! the bus: MPI is not CUDA-aware, so the all-reduce has to happen on the
      !! host, and the normalization stays with avexy_ibm_finish so the two
      !! paths cannot drift apart.
      !!
      !! kzb is the lower z bound of var, because the tendencies start at kb
      !! and the previous-step fields at kb-kh. An explicit-shape device dummy
      !! takes the address of the first element, so getting this wrong would
      !! silently read the array shifted by kh levels; call sites pass
      !! lbound(var,3).
      subroutine avexy_ibm_device(aver, var, kzb, II, IIs, lnan)
         use modmpi, only : avexy_ibm_finish
         implicit none

         integer,         intent(in)  :: kzb
         real   , device, intent(in)  :: var(ib-ih:ie+ih,jb-jh:je+jh,kzb:ke+kh)
         integer, device, intent(in)  :: II (ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc)
         real   ,         intent(out) :: aver(kb:ke+kh)
         integer,         intent(in)  :: IIs(kb:ke+kh)
         logical,         intent(in)  :: lnan

         real    :: s, averl_kb_nomask
         integer :: i, j, k

         !$acc parallel loop gang default(present) private(s)
         do k = kb, ke+kh
           s = 0.
           !$acc loop vector collapse(2) reduction(+:s)
           do j = jb, je
             do i = ib, ie
               s = s + var(i,j,k)*II(i,j,k)
             end do
           end do
           col_d(k) = s
         end do
         !$acc end parallel loop

         ! Only needed when the lowest level is entirely solid, exactly as on
         ! the host, so the extra pass is not paid for otherwise.
         averl_kb_nomask = 0.
         if ((.not. lnan) .and. (IIs(kb) == 0)) then
           s = 0.
           !$acc parallel loop collapse(2) default(present) reduction(+:s)
           do j = jb, je
             do i = ib, ie
               s = s + var(i,j,kb)
             end do
           end do
           !$acc end parallel loop
           averl_kb_nomask = s
         end if

         col_stage = col_d
         call avexy_ibm_finish(aver, col_stage, averl_kb_nomask, kb, ke, kh, IIs, lnan)

      end subroutine avexy_ibm_device

      !> Rank-local fluid-cell sum over j on the i = iplane column, weighted.
      !!
      !! The device counterpart of the sumy_ibm call in the u outflow branch of
      !! masscorr, whose i range there is the single outlet plane. weight is
      !! applied to var before the mask, in the association the host expression
      !! uses. kzb is the lower z bound of var; see avexy_ibm_device.
      subroutine sumy_ibm_column_device(total, var, kzb, II, iplane, weight)
         use modmpi, only : sum_ibm_reduce
         implicit none

         integer,         intent(in)  :: kzb
         real   , device, intent(in)  :: var(ib-ih:ie+ih,jb-jh:je+jh,kzb:ke+kh)
         integer, device, intent(in)  :: II (ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc)
         integer,         intent(in)  :: iplane
         real   ,         intent(in)  :: weight
         real   ,         intent(out) :: total(kb:ke)

         real    :: s
         integer :: j, k

         !$acc parallel loop gang default(present) private(s)
         do k = kb, ke
           s = 0.
           !$acc loop vector reduction(+:s)
           do j = jb, je
             s = s + (var(iplane,j,k)*weight)*II(iplane,j,k)
           end do
           col_d(k) = s
         end do
         !$acc end parallel loop

         col_stage(kb:ke) = col_d(kb:ke)
         call sum_ibm_reduce(total, col_stage(kb:ke), ke-kb+1)

      end subroutine sumy_ibm_column_device

      !> Rank-local fluid-cell sum over i on the j = jplane column, weighted.
      !!
      !! The device counterpart of the sumx_ibm call in the v outflow branch of
      !! masscorr. kzb is the lower z bound of var; see avexy_ibm_device.
      subroutine sumx_ibm_column_device(total, var, kzb, II, jplane, weight)
         use modmpi, only : sum_ibm_reduce
         implicit none

         integer,         intent(in)  :: kzb
         real   , device, intent(in)  :: var(ib-ih:ie+ih,jb-jh:je+jh,kzb:ke+kh)
         integer, device, intent(in)  :: II (ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc)
         integer,         intent(in)  :: jplane
         real   ,         intent(in)  :: weight
         real   ,         intent(out) :: total(kb:ke)

         real    :: s
         integer :: i, k

         !$acc parallel loop gang default(present) private(s)
         do k = kb, ke
           s = 0.
           !$acc loop vector reduction(+:s)
           do i = ib, ie
             s = s + (var(i,jplane,k)*weight)*II(i,jplane,k)
           end do
           col_d(k) = s
         end do
         !$acc end parallel loop

         col_stage(kb:ke) = col_d(kb:ke)
         call sum_ibm_reduce(total, col_stage(kb:ke), ke-kb+1)

      end subroutine sumx_ibm_column_device

      subroutine checkCUDA(istat, kernelname)
         implicit none
         integer, intent(in)          :: istat
         character(len=*), intent(in) :: kernelname
         if(istat /= cudaSuccess) then
            write(*,*) "Error in ", trim(kernelname), ": ", cudaGetErrorString(istat)
            error stop 1
         end if
      end subroutine checkCUDA

      attributes(device) subroutine tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
         implicit none
         integer :: tidx, tidy, tidz, stridex, stridey, stridez
         tidx = (blockIDx%x - 1) * blockDim%x + threadIdx%x
         tidy = (blockIDx%y - 1) * blockDim%y + threadIdx%y
         tidz = (blockIDx%z - 1) * blockDim%z + threadIdx%z
         stridex = gridDim%x * blockDim%x
         stridey = gridDim%y * blockDim%y
         stridez = gridDim%z * blockDim%z
      end subroutine tidandstride

      attributes(global) subroutine initfield(var, varvalue, halo_i, halo_j, halo_k)
         implicit none
         integer, value, intent(in) :: halo_i, halo_j, halo_k
         real, dimension(ib_d-halo_i:ie_d+halo_i, jb_d-halo_j:je_d+halo_j, kb_d:ke_d+halo_k), intent(inout) :: var
         real, value, intent(in) :: varvalue

         integer :: tidx, tidy, tidz, stridex, stridey, stridez
         integer :: i, j, k

         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

         do i = tidx-halo_i, ie_d+halo_i, stridex
            do j = tidy-halo_j, je_d+halo_j, stridey
               do k = tidz, ke_d+halo_k, stridez
                  var(i,j,k) = varvalue
               end do
            end do
         end do
      end subroutine initfield

      ! copy routines called inside advection, for kappa scheme of thlp
      attributes(global) subroutine thlptothlpc_cuda
         implicit none
         integer :: i, j, k, tidx, tidy, tidz, stridex, stridey, stridez
         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
         do k = tidz, ke_d, stridez
            do j = tidy, je_d, stridey
               do i = tidx, ie_d, stridex
                  thlpc_d(i,j,k) = thlp_d(i,j,k)
               end do
            end do
         end do
      end subroutine thlptothlpc_cuda

      attributes(global) subroutine thlpctothlp_cuda
         implicit none
         integer :: i, j, k, tidx, tidy, tidz, stridex, stridey, stridez
         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
         do k = tidz, ke_d, stridez
            do j = tidy, je_d, stridey
               do i = tidx, ie_d, stridex
                  thlp_d(i,j,k) = thlpc_d(i,j,k)
               end do
            end do
         end do
      end subroutine thlpctothlp_cuda

#endif
end module modcuda
