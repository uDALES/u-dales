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
                             eps1, numol, prandtlmoli, prandtlturb, grav, fkar2, &
                             fielddump_wants
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

   ! Bookkeeping for the post-Poisson handover.
   !
   ! Four routines below pull fields down for four different readers, and the
   ! sets overlap: fielddump and the restart both want u0, statsdump and the
   ! unported host routines both want ekm. Each field is pulled at most once
   ! per step, and which routine got there first does not matter - so rather
   ! than each routine carrying conditions about what the others already did,
   ! every transfer goes through a pull_ routine and is skipped if the field is
   ! already current. The dedup is then structural: listing a field in two
   ! routines is harmless, and adding a reader cannot silently double a copy.
   !
   ! updateDevice clears the flags, which is the one point in the step where
   ! every host field has just been handed to the device.
   ! The rule that governs which of them may be conditional: updateDevice
   ! uploads sixteen of these fields back from the host on every stage, so if
   ! one of those is not pulled this step, the host's stale copy is written
   ! over the device's current one and the solution itself is wrong - not just
   ! a dump. Those sixteen therefore belong in updateHostForUnportedRoutines,
   ! which runs unconditionally, and only ekm, ekh, tau_x, tau_y, tau_z and
   ! thl_flux - the six updateDevice does not upload - can be pulled on demand.
   ! assertRoundTripPulled below holds the two lists to that.
   integer, parameter :: F_U0 = 1, F_V0 = 2, F_W0 = 3, &
                         F_UM = 4, F_VM = 5, F_WM = 6, &
                         F_PRES0 = 7, F_THL0 = 8, F_THLM = 9, F_THL0C = 10, &
                         F_QT0 = 11, F_QTM = 12, F_SV0 = 13, F_SVM = 14, &
                         F_E120 = 15, F_E12M = 16, F_EKM = 17, F_EKH = 18, &
                         F_TAUX = 19, F_TAUY = 20, F_TAUZ = 21, F_THLFLUX = 22, &
                         F_COUNT = 22
   logical :: pulled(F_COUNT) = .false.
#if defined(UDALES_DEBUG)
   ! The round-trip invariant only holds inside the time loop. Initialisation
   ! and the device self-tests call updateDevice with nothing pulled, quite
   ! legitimately, so program.f90 arms the check once the loop is about to run.
   logical :: lcheck_round_trip = .false.
#endif

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

#if defined(UDALES_DEBUG)
      !> Every field updateDevice is about to upload must have been pulled.
      !!
      !! updateDevice writes the host copy of sixteen fields onto the device.
      !! If one of them was not brought down during the step just finishing,
      !! that upload replaces a current device value with a stale host one, and
      !! the error is in the solution rather than in an output file - it shows
      !! up as a slow divergence with nothing pointing at the cause.
      !!
      !! So the list here is the list of uploads in updateDevice, and it has to
      !! be kept beside it: adding an upload without adding the field to
      !! updateHostForUnportedRoutines is exactly the mistake this catches.
      !! Which is not hypothetical - pres0 was left out of that routine on the
      !! first attempt at this split, on the grounds that its only readers are
      !! the three that have their own routines.
      !!
      !! Inert until enableRoundTripCheck arms it, because initialisation and
      !! the device self-tests both call updateDevice with nothing pulled.
      subroutine assertRoundTripPulled
         use modmpi, only : myid
         implicit none

         if (.not. lcheck_round_trip) return

         call need(F_U0, 'u0') ; call need(F_V0, 'v0') ; call need(F_W0, 'w0')
         call need(F_UM, 'um') ; call need(F_VM, 'vm') ; call need(F_WM, 'wm')
         call need(F_PRES0, 'pres0')
         if (ltempeq) then
            call need(F_THL0, 'thl0') ; call need(F_THLM, 'thlm')
            if (iadv_thl == iadv_kappa) call need(F_THL0C, 'thl0c')
         end if
         if (lmoist) then
            call need(F_QT0, 'qt0') ; call need(F_QTM, 'qtm')
         end if
         if (nsv > 0) then
            call need(F_SV0, 'sv0') ; call need(F_SVM, 'svm')
         end if
         if (loneeqn) then
            call need(F_E120, 'e120') ; call need(F_E12M, 'e12m')
         end if

      contains

         subroutine need(id, name)
            implicit none
            integer,          intent(in) :: id
            character(len=*), intent(in) :: name

            if (.not. pulled(id)) then
               write(*,'(A,I0,3A)') 'Round-trip field not pulled on rank ', myid, ': ', &
                    trim(name), ' is about to be uploaded from a host copy this step never fetched.'
               error stop 1
            end if

         end subroutine need

      end subroutine assertRoundTripPulled

      !> Arm the round-trip check. Called once, as the time loop starts.
      subroutine enableRoundTripCheck
         implicit none
         lcheck_round_trip = .true.
      end subroutine enableRoundTripCheck
#endif

      subroutine updateDevice
         implicit none
         integer :: n

#if defined(UDALES_DEBUG)
         call assertRoundTripPulled
#endif

         ! Everything the host holds is about to be on the device, so nothing
         ! is pulled yet for the step this starts. See the declaration.
         pulled = .false.

         call updateFacetPropsDevice

         um_d = um
         vm_d = vm
         wm_d = wm

         u0_d = u0
         v0_d = v0
         w0_d = w0
         pres0_d = pres0

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
            thlm_d = thlm
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
            qtm_d  = qtm
            qt0_d = qt0
            call initfield<<<griddim,blockdim>>>(qtp_d, 0., ih, jh, kh)
            call checkCUDA( cudaGetLastError(), 'initfield qtp_d' )
            qt0av_d = qt0av
            if (lnudge) qtprof_d = qtprof
         end if

         if (nsv>0) then
            svm_d  = svm
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

      ! One routine per field rather than one generic helper taking the pair as
      ! arguments. The generic version was written first and cost 130 seconds on
      ! a 256^3 run: an assumed-shape dummy turns a whole-array device-to-host
      ! assignment into a descriptor copy, and a 137 MB field cannot afford
      ! that. With no dummies the compiler sees exactly the assignment it saw
      ! before this split existed.

      subroutine pull_u0
         implicit none
         if (pulled(F_U0)) return
         u0 = u0_d
         pulled(F_U0) = .true.
      end subroutine pull_u0

      subroutine pull_v0
         implicit none
         if (pulled(F_V0)) return
         v0 = v0_d
         pulled(F_V0) = .true.
      end subroutine pull_v0

      subroutine pull_w0
         implicit none
         if (pulled(F_W0)) return
         w0 = w0_d
         pulled(F_W0) = .true.
      end subroutine pull_w0

      subroutine pull_um
         implicit none
         if (pulled(F_UM)) return
         um = um_d
         pulled(F_UM) = .true.
      end subroutine pull_um

      subroutine pull_vm
         implicit none
         if (pulled(F_VM)) return
         vm = vm_d
         pulled(F_VM) = .true.
      end subroutine pull_vm

      subroutine pull_wm
         implicit none
         if (pulled(F_WM)) return
         wm = wm_d
         pulled(F_WM) = .true.
      end subroutine pull_wm

      subroutine pull_pres0
         implicit none
         if (pulled(F_PRES0)) return
         pres0 = pres0_d
         pulled(F_PRES0) = .true.
      end subroutine pull_pres0

      subroutine pull_thl0
         implicit none
         if (pulled(F_THL0)) return
         thl0 = thl0_d
         pulled(F_THL0) = .true.
      end subroutine pull_thl0

      subroutine pull_thlm
         implicit none
         if (pulled(F_THLM)) return
         thlm = thlm_d
         pulled(F_THLM) = .true.
      end subroutine pull_thlm

      subroutine pull_thl0c
         implicit none
         if (pulled(F_THL0C)) return
         thl0c = thl0c_d
         pulled(F_THL0C) = .true.
      end subroutine pull_thl0c

      subroutine pull_qt0
         implicit none
         if (pulled(F_QT0)) return
         qt0 = qt0_d
         pulled(F_QT0) = .true.
      end subroutine pull_qt0

      subroutine pull_qtm
         implicit none
         if (pulled(F_QTM)) return
         qtm = qtm_d
         pulled(F_QTM) = .true.
      end subroutine pull_qtm

      subroutine pull_sv0
         implicit none
         if (pulled(F_SV0)) return
         sv0 = sv0_d
         pulled(F_SV0) = .true.
      end subroutine pull_sv0

      subroutine pull_svm
         implicit none
         if (pulled(F_SVM)) return
         svm = svm_d
         pulled(F_SVM) = .true.
      end subroutine pull_svm

      subroutine pull_e120
         implicit none
         if (pulled(F_E120)) return
         e120 = e120_d
         pulled(F_E120) = .true.
      end subroutine pull_e120

      subroutine pull_e12m
         implicit none
         if (pulled(F_E12M)) return
         e12m = e12m_d
         pulled(F_E12M) = .true.
      end subroutine pull_e12m

      subroutine pull_ekm
         implicit none
         if (pulled(F_EKM)) return
         ekm = ekm_d
         pulled(F_EKM) = .true.
      end subroutine pull_ekm

      subroutine pull_ekh
         implicit none
         if (pulled(F_EKH)) return
         ekh = ekh_d
         pulled(F_EKH) = .true.
      end subroutine pull_ekh

      subroutine pull_tau_x
         implicit none
         if (pulled(F_TAUX)) return
         tau_x = tau_x_d
         pulled(F_TAUX) = .true.
      end subroutine pull_tau_x

      subroutine pull_tau_y
         implicit none
         if (pulled(F_TAUY)) return
         tau_y = tau_y_d
         pulled(F_TAUY) = .true.
      end subroutine pull_tau_y

      subroutine pull_tau_z
         implicit none
         if (pulled(F_TAUZ)) return
         tau_z = tau_z_d
         pulled(F_TAUZ) = .true.
      end subroutine pull_tau_z

      subroutine pull_thl_flux
         implicit none
         if (pulled(F_THLFLUX)) return
         thl_flux = thl_flux_d
         pulled(F_THLFLUX) = .true.
      end subroutine pull_thl_flux

      !> Bring down what fielddump is about to read.
      !!
      !! Called only on the steps fielddump writes, so the fields that exist
      !! for it alone - the wall stresses and the heat flux - cross the bus a
      !! handful of times a run rather than on every stage.
      !!
      !! u0, v0 and w0 are unconditional because fielddump forms div from them
      !! before it looks at fieldvars at all. Everything else is asked for by
      !! name, through fielddump's own predicate, so a field cannot be pulled
      !! under one condition and read under another.
      subroutine updateHostForFielddump
         implicit none

         call pull_u0
         call pull_v0
         call pull_w0

         if (fielddump_wants('p0')) call pull_pres0
         if (ltempeq .and. fielddump_wants('th')) call pull_thl0
         if (lmoist  .and. fielddump_wants('qt')) call pull_qt0

         ! sv0 comes down whole, so any one of the five scalar codes is enough.
         if (nsv > 0) then
            if (fielddump_wants('s1') .or. fielddump_wants('s2') .or. &
                fielddump_wants('s3') .or. fielddump_wants('s4') .or. &
                fielddump_wants('s5')) call pull_sv0
         end if

         if (fielddump_wants('tx')) call pull_tau_x
         if (fielddump_wants('ty')) call pull_tau_y
         if (fielddump_wants('tz')) call pull_tau_z
         if (ltempeq .and. fielddump_wants('hf')) call pull_thl_flux

      end subroutine updateHostForFielddump

      !> Bring down what statsdump is about to read.
      !!
      !! Called only on the steps statsdump samples. statsdump reads the
      !! previous-stage fields, not the current ones - um rather than u0 - plus
      !! both eddy diffusivities and the pressure. It does not read u0, v0 or
      !! w0: the only statistics that would are in tkestats, which is reached
      !! from genstats, whose one call site is commented out.
      !!
      !! The tree diagnostics statsdump also reads are handed over by
      !! vegetation's own updateVegDiagHost, called beside this one under the
      !! same predicate. They cannot be folded in here, because vegetation
      !! reaches modcuda and the reverse would be a cycle.
      subroutine updateHostForStatsdump
         implicit none

         call pull_um
         call pull_vm
         call pull_wm
         call pull_pres0
         call pull_ekm
         call pull_ekh

         if (ltempeq) call pull_thlm
         if (lmoist)  call pull_qtm
         if (nsv > 0) call pull_svm

      end subroutine updateHostForStatsdump

      !> Bring down what writerestartfiles is about to write.
      !!
      !! Called from inside that routine's own guard rather than from the time
      !! loop, so the restart schedule is evaluated once. Reproducing it here
      !! would mean repeating a filesystem inquire and an MPI broadcast every
      !! step to answer a question the writer already asks.
      !!
      !! By that point boundary has written the top and lateral planes of these
      !! fields on the host, which is what belongs in a restart - and the pull
      !! routines leave them alone, because updateHostForUnportedRoutines
      !! already marked them pulled earlier in the same step. Everything here
      !! is therefore a no-op today; the routine exists so that when that
      !! unconditional pull goes away, the restart keeps working. ql0 and ql0h
      !! are also written to the restart and are absent on purpose:
      !! thermodynamics produces them on the host and there is no device copy
      !! to fetch.
      subroutine updateHostForRestart
         implicit none

         call pull_u0
         call pull_v0
         call pull_w0
         call pull_pres0
         call pull_ekm

         if (ltempeq) call pull_thl0
         if (lmoist)  call pull_qt0
         if (nsv > 0) call pull_sv0
         if (loneeqn) call pull_e120

      end subroutine updateHostForRestart

      !> Bring down what the host routines still in the loop are about to use.
      !!
      !! This is not a bin of leftovers, it is the correctness path. boundary
      !! writes the top, bottom and lateral planes of every prognostic field on
      !! the host, and thermodynamics reads thl0 and qt0 and writes ql0, thv0h
      !! and the slab averages - and the next updateDevice uploads all of it
      !! back. Skip a field here and the device gets the previous step's value
      !! for it, so this runs on every stage and cannot be made conditional.
      !!
      !! It is also the one routine in this group that shrinks. Every entry
      !! below leaves when the routine that needs it is ported, and once
      !! boundary and thermodynamics are on the device the routine goes with
      !! them - unlike the three above, which survive a fully ported loop
      !! because their readers write files from host memory.
      !!
      !! ekm and ekh are here for boundary's fluxtop and fluxtopscal, which
      !! take them as the diffusivity for a flux-type top boundary.
      !!
      !! pres0 is here for a different reason. Its only host readers are
      !! statsdump, fielddump and the restart, all of which have their own
      !! routine - but updateDevice uploads it, so leaving it to them would put
      !! a stale pressure back on the device on every step between samples.
      !! The same argument keeps thlm, qtm and svm here even though boundary's
      !! writes to them are the only thing downstream that touches them.
      subroutine updateHostForUnportedRoutines
         implicit none

         call pull_u0
         call pull_v0
         call pull_w0
         call pull_um
         call pull_vm
         call pull_wm
         call pull_pres0
         call pull_ekm
         call pull_ekh

         if (ltempeq) then
            call pull_thl0
            call pull_thlm
            if (iadv_thl == iadv_kappa) call pull_thl0c
         end if

         if (lmoist) then
            call pull_qt0
            call pull_qtm
         end if

         if (nsv > 0) then
            call pull_sv0
            call pull_svm
#if defined(UDALES_DEBUG)
            if (.not. all(ieee_is_finite(sv0(ib:ie, jb:je, kb:ke, 1:nsv)))) then
               write(*,*) 'Non-finite scalar value detected after the GPU pressure step.'
               error stop 1
            end if
#endif
         end if

         if (loneeqn) then
            call pull_e120
            call pull_e12m
         end if

      end subroutine updateHostForUnportedRoutines

#if defined(UDALES_DEBUG)
      !> Abort if a host field a reader is about to use has drifted from the
      !! device copy it mirrors.
      !!
      !! Called from the time loop immediately before fielddump and before
      !! statsdump, and asks the question that matters at that moment: is every
      !! field this reader will look at the one the device holds?
      !!
      !! The lists below are written from what the reader reads. They are a
      !! second derivation of the lists in updateHostForFielddump and
      !! updateHostForStatsdump, on purpose - a check that shares its
      !! derivation with the thing it checks cannot catch a wrong list. The two
      !! have to agree, and when they do not this aborts naming the field.
      !!
      !! What it is defending against: those transfers are conditional now, and
      !! every way a condition can be wrong has the same symptom - a dump
      !! quietly written from the previous step, on a configuration nobody
      !! runs, with no error anywhere.
      !!
      !! Placement is load-bearing. It has to run before boundary, which writes
      !! the top and lateral planes on the host and so makes host and device
      !! legitimately differ until the next updateDevice.
      !!
      !! Debug builds only, and it copies everything down a second time, so it
      !! is not cheap - but the GPU parity cases are 64^3 and run in a second.
      subroutine assertHostMatchesDevice(label)
         use modmpi, only : myid
         implicit none

         character(len=*), intent(in) :: label
         integer :: n

         select case (label)

         case ('fielddump')
            ! div is formed from these before fieldvars is consulted at all.
            call check3(label, 'u0', u0, u0_d)
            call check3(label, 'v0', v0, v0_d)
            call check3(label, 'w0', w0, w0_d)

            if (fielddump_wants('p0')) call check3(label, 'pres0', pres0, pres0_d)
            if (ltempeq .and. fielddump_wants('th')) call check3(label, 'thl0', thl0, thl0_d)
            if (lmoist  .and. fielddump_wants('qt')) call check3(label, 'qt0',  qt0,  qt0_d)

            if (nsv > 0) then
               if (fielddump_wants('s1') .or. fielddump_wants('s2') .or. &
                   fielddump_wants('s3') .or. fielddump_wants('s4') .or. &
                   fielddump_wants('s5')) then
                  do n = 1, nsv
                     call check3(label, 'sv0', sv0(:,:,:,n), sv0_d(:,:,:,n))
                  end do
               end if
            end if

            if (fielddump_wants('tx')) call check3(label, 'tau_x', tau_x, tau_x_d)
            if (fielddump_wants('ty')) call check3(label, 'tau_y', tau_y, tau_y_d)
            if (fielddump_wants('tz')) call check3(label, 'tau_z', tau_z, tau_z_d)
            if (ltempeq .and. fielddump_wants('hf')) &
               call check3(label, 'thl_flux', thl_flux, thl_flux_d)

         case ('statsdump')
            ! The previous-stage fields, both diffusivities and the pressure.
            call check3(label, 'um', um, um_d)
            call check3(label, 'vm', vm, vm_d)
            call check3(label, 'wm', wm, wm_d)
            call check3(label, 'pres0', pres0, pres0_d)
            call check3(label, 'ekm', ekm, ekm_d)
            call check3(label, 'ekh', ekh, ekh_d)

            if (ltempeq) call check3(label, 'thlm', thlm, thlm_d)
            if (lmoist)  call check3(label, 'qtm',  qtm,  qtm_d)
            if (nsv > 0) then
               do n = 1, nsv
                  call check3(label, 'svm', svm(:,:,:,n), svm_d(:,:,:,n))
               end do
            end if

         case default
            write(*,'(3A)') 'assertHostMatchesDevice: no field list for reader ', trim(label), '.'
            error stop 1

         end select

      contains

         !> Bring one device field down into scratch and require exact equality.
         !!
         !! Exact, not a tolerance: the host copy was produced by copying this
         !! same device array, so anything other than equality means the copy
         !! did not happen this step. A tolerance would let a field that is
         !! merely one step stale, and therefore close, pass unnoticed.
         subroutine check3(where, name, host, dev)
            implicit none
            character(len=*), intent(in) :: where, name
            real,             intent(in) :: host(:,:,:)
            real, device,     intent(in) :: dev(:,:,:)

            real, allocatable :: back(:,:,:)
            integer :: bad
            real    :: worst

            allocate(back(size(dev,1), size(dev,2), size(dev,3)))
            back = dev

            bad = count(back /= host)
            if (bad > 0) then
               worst = maxval(abs(back - host))
               write(*,'(A,I0,7A,I0,A,I0,A,ES12.4)') 'Stale host field on rank ', myid, &
                    ': ', trim(name), ' differs from ', trim(name), '_d at ', trim(where), &
                    ' - ', bad, ' of ', size(host), ' elements, worst ', worst
               deallocate(back)
               error stop 1
            end if

            deallocate(back)

         end subroutine check3

      end subroutine assertHostMatchesDevice
#endif

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
