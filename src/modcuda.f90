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
                             ltempeq, lmoist, nsv, lchem, lles, lbuoyancy, lscasrc, lscasrcl, &
                             BCxm, BCxm_driver, libm, nfcts, idriver, &
                             iadv_sv, iadv_thl, iadv_kappa, iadv_upw, &
                             xh, &
                             eps1, numol, prandtlmoli, prandtlturb, grav, fkar2, &
                             fielddump_wants
   use modfields,      only: u0, v0, w0, pres0, e120, thl0, thl0c, qt0, sv0, &
                             ql0, ql0h, thl0h, qt0h, presf, presh, exnf, exnh, IIw, &
                             up, vp, wp, e12p, thlp, thlpc, qtp, svp, &
                             um, vm, wm, e12m, thlm, qtm, svm, &
                             tau_x, tau_y, tau_z, thl_flux, momfluxb, tfluxb, &
                             u0av, v0av, thl0av, qt0av, sv0av, dthvdz, ug, vg, whls, tsc, &
                             dpdxl, dpdyl, thv0h, thvh, thlpcar, &
                             dudxls, dudyls, dvdxls, dvdyls, dthldxls, dthldyls, dqtdxls, dqtdyls, dqtdtls, &
                             uprof, vprof, thlprof, qtprof, svprof, e12prof, &
                             IIc, IIu, IIv, scalar_source_tendency
   use modsubgriddata, only: lsmagorinsky, lvreman, loneeqn, ldelta, lbuoycorr, &
                             ekm, ekh, &
                             sbshr, sbbuo, sbdiss, zlt, damp, csz, &
                             cn, cm, ch1, ch2, ce1, ce2, dampmin, prandtli, c_vreman
   use modsurfdata,    only: thvs
   use initfac,        only: facT, fachf, facef, fachfi, facefi, facqsat, fachurel, facf
   use modinletdata,   only: u0driver, umdriver, v0driver, vmdriver, w0driver, wmdriver, &
                             thl0driver, thlmdriver, qt0driver, qtmdriver, &
                             sv0driver, svmdriver
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
                                uprof_d(:), vprof_d(:), thlprof_d(:), qtprof_d(:), e12prof_d(:)

   real, device, allocatable :: delta_d(:, :), csz_d(:,:), sv0av_d(:,:), svprof_d(:,:), u0driver_d(:,:)

   ! The rest of the driver inlet planes. Small - one j-k plane each - and
   ! only allocated for a run that reads a driver file, but the boundary
   ! conditions that use them are on the device now, so they have to be too.
   real, device, allocatable :: umdriver_d(:,:), v0driver_d(:,:), vmdriver_d(:,:), &
                                w0driver_d(:,:), wmdriver_d(:,:), &
                                thl0driver_d(:,:), thlmdriver_d(:,:), &
                                qt0driver_d(:,:), qtmdriver_d(:,:)
   real, device, allocatable :: sv0driver_d(:,:,:), svmdriver_d(:,:,:)

   ! Staging for the recycle plane a driver-generation run records. The
   ! plane is gathered on the device into one of these and crosses the bus
   ! contiguously; the alternative, a strided section copy, would be a
   ! two-dimensional transfer eight bytes wide. Two shapes because the
   ! scalars carry a wider halo than the momentum fields.
   real, device, allocatable :: plane_d(:,:), planec_d(:,:)
   real, allocatable, pinned :: plane_stage(:,:), planec_stage(:,:)

   real, device, allocatable :: u0_d(:,:,:), v0_d(:,:,:), w0_d(:,:,:), pres0_d(:,:,:), e120_d(:,:,:), thl0_d(:,:,:), thl0c_d(:,:,:), qt0_d(:,:,:), sv0_d(:,:,:,:)
   real, device, allocatable :: up_d(:,:,:), vp_d(:,:,:), wp_d(:,:,:), e12p_d(:,:,:), thlp_d(:,:,:), thlpc_d(:,:,:), qtp_d(:,:,:), svp_d(:,:,:,:)
   real, device, allocatable :: um_d(:,:,:), vm_d(:,:,:), wm_d(:,:,:), e12m_d(:,:,:), thlm_d(:,:,:), qtm_d(:,:,:), svm_d(:,:,:,:)
   real, device, allocatable :: p_d(:,:,:), pup_d(:,:,:), pvp_d(:,:,:), pwp_d(:,:,:)
   real, device, allocatable :: tau_x_d(:,:,:), tau_y_d(:,:,:), tau_z_d(:,:,:), thl_flux_d(:,:,:)
   real, device, allocatable :: momfluxb_d(:,:,:), tfluxb_d(:,:,:)
   real, device, allocatable :: dthvdz_d(:,:,:)
   real, device, allocatable :: ekm_d(:,:,:), ekh_d(:,:,:), sbshr_d(:,:,:), sbbuo_d(:,:,:), sbdiss_d(:,:,:), zlt_d(:,:,:), damp_d(:,:,:)
   real, device, allocatable :: thv0h_d(:,:,:)

   ! What thermodynamics derives on the device, and the hydrostatic columns
   ! its kernels read. ql0_d and qt0_d exist for a dry-but-thermal run too,
   ! because the virtual temperature the thvf reduction forms is the full
   ! moist expression whether or not lmoist is on: nothing writes ql0 or qt0
   ! then, but thl0 moves and the product does with it.
   real, device, allocatable :: ql0_d(:,:,:), ql0h_d(:,:,:), thl0h_d(:,:,:), qt0h_d(:,:,:)
   real, device, allocatable :: presf_d(:), presh_d(:), exnf_d(:), exnh_d(:)
   real, device, allocatable :: dumu_d(:,:,:), duml_d(:,:,:)
   real, device, allocatable :: dummyNO_d(:,:,:), dummyNO2_d(:,:,:), dummyO3_d(:,:,:)

   integer, device, allocatable :: IIc_d(:,:,:), IIu_d(:,:,:), IIv_d(:,:,:), IIw_d(:,:,:)

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
   !
   ! Every transfer here is on demand now. It was not: updateDevice used to
   ! upload sixteen of these fields back from the host on every stage, so a
   ! field not pulled during the step had its stale host copy written over the
   ! device's current one and the error landed in the solution rather than in
   ! a file. That was the price of thermodynamics running on the host, and it
   ! is what porting thermodynamics bought back - updateDevice no longer
   ! uploads a single prognostic field, so nothing here is load-bearing for
   ! correctness beyond the reader that asks for it.
   integer, parameter :: F_U0 = 1, F_V0 = 2, F_W0 = 3, &
                         F_UM = 4, F_VM = 5, F_WM = 6, &
                         F_PRES0 = 7, F_THL0 = 8, F_THLM = 9, F_THL0C = 10, &
                         F_QT0 = 11, F_QTM = 12, F_SV0 = 13, F_SVM = 14, &
                         F_E120 = 15, F_E12M = 16, F_EKM = 17, F_EKH = 18, &
                         F_TAUX = 19, F_TAUY = 20, F_TAUZ = 21, F_THLFLUX = 22, &
                         F_QL0 = 23, F_QL0H = 24, &
                         F_COUNT = 24
   logical :: pulled(F_COUNT) = .false.
#if defined(UDALES_DEBUG)
   ! Neither assertion below means anything before the time loop starts.
   ! Initialisation and the device self-tests reach the pull routines with the
   ! bitmap uncleared and with host and device deliberately out of step - the
   ! self-tests poison one side and compare - so both checks stay inert until
   ! program.f90 arms them, which it does once the loop is about to run.
   !
   ! The asserts test this themselves rather than their call sites testing it,
   ! so a routine reached from a self-test cannot fire one by accident.
   logical :: lchecks_armed = .false.
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

         ! Seeded here and then owned by the device, like e12prof_d below.
         ! readinitfiles fills the host pair and inittimedep interpolates it
         ! to t = 0, both before this point; from the first stage onwards the
         ! only writer is timedepnudge_device, and it writes these.
         !
         ! Unconditional rather than under the readers' guards. updateDevice
         ! used to upload the pair per stage under
         ! (lnudge .and. lnudgevel) .or. BCxm_profile .or. BCym_profile, and
         ! that list had already been wrong once - it has to name every
         ! consumer of a profile inlet, and there is no reason to keep paying
         ! attention to it for one column.
         allocate(uprof_d(kb:ke+kh)); uprof_d = uprof
         allocate(vprof_d(kb:ke+kh)); vprof_d = vprof
         allocate(u0driver_d(jb-jh:je+jh,kb-kh:ke+kh))

         ! The remaining driver planes, for a run that reads one, filled as
         ! they are allocated - see allocDriverPlanesDevice for why the two
         ! are one step. u0driver_d above stays unconditional because the
         ! device self-tests use it as scratch for the driver pressure kernel.
         if (BCxm == BCxm_driver) call allocDriverPlanesDevice

         ! And the staging for a run that records one.
         if (idriver == 1) then
            allocate(plane_d(jb-jh:je+jh,kb-kh:ke+kh))
            allocate(plane_stage(jb-jh:je+jh,kb-kh:ke+kh))
            if (nsv > 0) then
               allocate(planec_d(jb-jhc:je+jhc,kb-khc:ke+khc))
               allocate(planec_stage(jb-jhc:je+jhc,kb-khc:ke+khc))
            end if
         end if

         if (loneeqn) then
            ! Set once: the inlet turbulence profile is read at startup and
            ! nothing in the loop writes it, unlike uprof and vprof which
            ! timedep can move.
            allocate(e12prof_d(kb:ke+kh))
            e12prof_d = e12prof
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
            allocate(thlprof_d(kb:ke+kh)); thlprof_d = thlprof

            ! readinitfiles is the only writer of thlpcar, and it has run by
            ! now, so this seed is the whole of it. updateDevice used to
            ! repeat the copy per stage under ltrees .and. lmoist, which was
            ! dead the moment it was written.
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
            allocate(qtm_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
            allocate(qtp_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
            allocate(qt0av_d(kb:ke+kh))
            allocate(qtprof_d(kb:ke+kh)); qtprof_d = qtprof

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
            ! Set once, for the same reason as e12prof_d above: timedep has
            ! no scalar profile to interpolate, so this was being uploaded
            ! unchanged on every stage of every step.
            allocate(svprof_d(kb:ke+kh,nsv))
            svprof_d = svprof
            if (nsv==3 .and. lchem) then
               allocate(dummyNO_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
               allocate(dummyNO2_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
               allocate(dummyO3_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
            end if
            if (lscasrc .or. lscasrcl) then
               allocate(scalar_source_tendency_d(ib:ie,jb:je,kb:ke,nsv))
               scalar_source_tendency_d = scalar_source_tendency
            end if
         end if

         ! Seeded, not left as whatever the allocation returned. subgrid
         ! writes the interior before anything reads it, but nothing writes the
         ! outer halo columns until closurebc and the halo exchange do, and the
         ! host copies are a defined starting point for both.
         allocate(ekm_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)); ekm_d = ekm
         allocate(ekh_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)); ekh_d = ekh
         allocate(dthvdz_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))

         ! What thermodynamics owns on the device, seeded from the host copies
         ! rather than left as whatever the allocation returned. readinitfiles
         ! runs the host thermodynamics once, before this point, so every one
         ! of these already holds the value the first stage of the first step
         ! expects to read - which for presf and exnf is not a convenience but
         ! the definition: thermodynamics' first saturation pass runs before
         ! diagfld rebuilds them, so it reads the previous call's profiles, and
         ! at the first stage the previous call is the one in readinitfiles.
         !
         ! dthvdz_d and thv0h_d are seeded for a narrower reason: subgrid and
         ! the buoyancy term read them on the first stage of the first step,
         ! which is before the first calthv writes them.
         dthvdz_d = dthvdz
         if (ltempeq) then
            allocate(thl0h_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)); thl0h_d = thl0h
            allocate(qt0h_d (ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)); qt0h_d  = qt0h
            allocate(ql0_d  (ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)); ql0_d   = ql0
            allocate(ql0h_d (ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh));    ql0h_d  = ql0h
            allocate(presf_d(kb:ke+kh)); presf_d = presf
            allocate(presh_d(kb:ke+kh)); presh_d = presh
            allocate(exnf_d (kb:ke+kh)); exnf_d  = exnf
            allocate(exnh_d (kb:ke+kh)); exnh_d  = exnh
            thv0h_d = thv0h
         end if

         ! qt0 is read by the thvf reduction whether or not there is moisture,
         ! and nothing writes it when there is none - so it is allocated for
         ! the wider condition and seeded here, once, for the narrower case.
         if (ltempeq .or. lmoist) then
            allocate(qt0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
            if (.not. lmoist) qt0_d = qt0
         end if

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

         ! The fluid-cell masks and the staging column for every rank-local
         ! reduction on the device. These used to be conditional on masscorr's
         ! flow-rate controllers, which were the only readers; diagfld's slab
         ! averages are on the device now and they run on every stage of every
         ! configuration, so u, v and c are unconditional. Four integer masks
         ! is real memory - about 70 MB each at 256^3 - but the host already
         ! carries seven of them, and the alternative is bringing three
         ! prognostic fields down per stage to average them there.
         allocate(col_d(kb:ke+kh))
         allocate(col_stage(kb:ke+kh))
         allocate(IIu_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc)); IIu_d = IIu
         allocate(IIv_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc)); IIv_d = IIv
         allocate(IIc_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc)); IIc_d = IIc
         ! The w mask has one reader, the thvh average, and that only exists
         ! when there is a temperature to average.
         if (ltempeq) then
            allocate(IIw_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc)); IIw_d = IIw
         end if

         call uploadPrognosticFieldsDevice

         ! The tendencies start empty, exactly as initfields starts the host
         ! ones - allocate(up(...)) ; up=0. and so on. This used to be done by
         ! the updateDevice at the top of the first stage; now that the clear
         ! belongs to the end of a stage, the first stage has no stage before
         ! it to have done it, and advection would accumulate into whatever
         ! the allocation returned.
         call resetTendenciesDevice

      end subroutine initCUDA

      !> The one and only upload of the prognostic fields.
      !!
      !! updateDevice used to do this at the top of every stage. That is what
      !! made the host copies authoritative and forced the unconditional pull
      !! that answered it, and it was necessary for exactly as long as
      !! thermodynamics ran on the host. Every writer in the loop is on the
      !! device now, so the handover happens once, from initCUDA, and the
      !! device holds these from there until exitCUDA - nothing uploads them
      !! again during the run.
      !!
      !! Which also means a field left out of this list is not merely stale, it
      !! is never initialised at all, on a run of any length and in every
      !! configuration that allocates it. That is a loud failure rather than a
      !! quiet one, and it is why one upload is enough.
      !!
      !! It is a routine rather than a block inside initCUDA because the device
      !! self-tests seed these arrays and have to put them back, and there is
      !! no longer a per-stage upload to do it for them.
      subroutine uploadPrognosticFieldsDevice
         implicit none

         um_d = um
         vm_d = vm
         wm_d = wm
         u0_d = u0
         v0_d = v0
         w0_d = w0
         pres0_d = pres0
         if (loneeqn) then
            e120_d = e120
            e12m_d = e12m
         end if
         if (ltempeq) then
            thl0_d = thl0
            thlm_d = thlm
            if (iadv_thl == iadv_kappa) thl0c_d = thl0c
         end if
         ! qt0_d is allocated for a run with temperature and no moisture too;
         ! qtm_d is not.
         if (ltempeq .or. lmoist) qt0_d = qt0
         if (lmoist) qtm_d = qtm
         if (nsv > 0) then
            sv0_d = sv0
            svm_d = svm
         end if

      end subroutine uploadPrognosticFieldsDevice

      subroutine exitCUDA
         implicit none
         deallocate(dxf_d, dxhi_d, dzf_d, dzf2_d, dzfi_d, dzfi5_d, dzfiq_d, dzh_d, dzhi_d, dzh2i_d, dzhiq_d, delta_d)
         deallocate(dxdydzfi_d, dxdydzhi_d, dvcell_d)
         deallocate(u0_d, v0_d, w0_d, pres0_d, um_d, vm_d, wm_d, up_d, vp_d, wp_d)
         deallocate(tau_x_d, tau_y_d, tau_z_d, momfluxb_d)
         deallocate(u0av_d, v0av_d, ug_d, vg_d, whls_d, tsc_d)
         deallocate(dpdxl_d, dpdyl_d, dudxls_d, dudyls_d, dvdxls_d, dvdyls_d)
         deallocate(uprof_d, vprof_d, u0driver_d)
         if (BCxm == BCxm_driver) call deallocDriverPlanesDevice
         if (allocated(plane_d)) deallocate(plane_d, plane_stage)
         if (allocated(planec_d)) deallocate(planec_d, planec_stage)
         if (loneeqn) deallocate(e120_d, e12m_d, e12p_d, sbshr_d, sbbuo_d, sbdiss_d, zlt_d, e12prof_d)
         if (ltempeq) then
            deallocate(thl0_d, thlm_d, thlp_d, thl_flux_d)
            if (iadv_thl == iadv_kappa) deallocate(thl0c_d, thlpc_d)
            deallocate(thv0h_d, thvh_d, thl0av_d, thlprof_d, thlpcar_d)
            deallocate(dthldxls_d, dthldyls_d)
            deallocate(tfluxb_d)
         end if
         if (lmoist) deallocate(qtm_d, qtp_d, qt0av_d, qtprof_d, dqtdxls_d, dqtdyls_d, dqtdtls_d)
         if (ltempeq .or. lmoist) deallocate(qt0_d)
         if (ltempeq) deallocate(thl0h_d, qt0h_d, ql0_d, ql0h_d, presf_d, presh_d, exnf_d, exnh_d)
         if (nsv>0) then
            deallocate(sv0_d, svm_d, svp_d, sv0av_d, svprof_d)
            if (lscasrc .or. lscasrcl) deallocate(scalar_source_tendency_d)
            if (nsv==3 .and. lchem) deallocate(dummyNO_d, dummyNO2_d, dummyO3_d)
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
         deallocate(col_d, col_stage)
         deallocate(IIu_d, IIv_d, IIc_d)
         if (allocated(IIw_d)) deallocate(IIw_d)
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
      !> Arm the host/device assertions. Called once, as the loop starts.
      subroutine enableDeviceHostChecks
         implicit none
         lchecks_armed = .true.
      end subroutine enableDeviceHostChecks
#endif

      !> Send the six slab averages the host reduces back up to the device.
      !!
      !! No prognostic field crosses here any more. These six are the only
      !! host-to-device transfers left in the time loop, and not because
      !! something is unported: avexy_ibm_device reduces over the domain on the
      !! device, but the rank sum in the middle of it is an MPI_ALLREDUCE, so
      !! the finished column only ever exists on the host.
      subroutine updateDevice
         implicit none

         ! The host copies are about to stop being the ones that matter for
         ! this step: whatever a reader fetched last step is stale from here
         ! on. See the declaration.
         pulled = .false.

         call updateFacetPropsDevice

         u0av_d = u0av
         v0av_d = v0av

         if (ltempeq) then
            thvh_d     = thvh
            thl0av_d   = thl0av
         end if

         if (lmoist) then
            qt0av_d = qt0av
         end if

         if (nsv>0) then
            sv0av_d = sv0av
         end if
      end subroutine updateDevice

      !> Clear the tendency accumulators, once per Runge-Kutta stage.
      !!
      !! Called from the end of tstep_integrate, which is where the host does
      !! the same thing, and for the same reason: advection and everything
      !! after it add into these, so the stage that is finishing has to leave
      !! them empty for the stage that follows. Having the two implementations
      !! zero them at the same point in the step is the whole reason this is a
      !! routine - it used to happen at the top of the next stage, inside
      !! updateDevice, which is the same instant in the sequence and a
      !! different place in the source.
      !!
      !! Safe to move because nothing between the end of tstep_integrate and
      !! the next advection touches a tendency: exchange_halos, checksim, the
      !! dumps, boundary_conditions, thermodynamics and the restart writer all
      !! read prognostic fields. Every routine that does add into one -
      !! coriolis, forces, lstend, nudge, the wall functions, masscorr,
      !! grwdamp and bcpup - runs between updateDevice and poisson.
      !!
      !! thlpc_d is here and is not in the host's list, which is not an
      !! oversight on either side. advection assigns thlpc from thlp before it
      !! accumulates, on both, so the interior never needs clearing; what the
      !! host relies on instead is the thlpc = 0. initfields does at
      !! allocation, because advecc_kappa_add writes the halo as well and only
      !! the interior is copied back. Zeroing it per stage is what the device
      !! has always done and costs one kernel on a kappa run.
      subroutine resetTendenciesDevice
         implicit none
         integer :: n

         call initfield<<<griddim,blockdim>>>(up_d, 0., ih, jh, kh)
         call checkCUDA( cudaGetLastError(), 'initfield up_d' )

         call initfield<<<griddim,blockdim>>>(vp_d, 0., ih, jh, kh)
         call checkCUDA( cudaGetLastError(), 'initfield vp_d' )

         call initfield<<<griddim,blockdim>>>(wp_d, 0., ih, jh, kh)
         call checkCUDA( cudaGetLastError(), 'initfield wp_d' )

         if (loneeqn) then
            call initfield<<<griddim,blockdim>>>(e12p_d, 0., ih, jh, kh)
            call checkCUDA( cudaGetLastError(), 'initfield e12p_d' )
         end if

         if (ltempeq) then
            call initfield<<<griddim,blockdim>>>(thlp_d, 0., ih, jh, kh)
            call checkCUDA( cudaGetLastError(), 'initfield thlp_d' )

            if (iadv_thl == iadv_kappa) then
               call initfield<<<griddim,blockdim>>>(thlpc_d, 0., ihc, jhc, khc)
               call checkCUDA( cudaGetLastError(), 'initfield thlpc_d' )
            end if
         end if

         if (lmoist) then
            call initfield<<<griddim,blockdim>>>(qtp_d, 0., ih, jh, kh)
            call checkCUDA( cudaGetLastError(), 'initfield qtp_d' )
         end if

         if (nsv>0) then
            do n = 1, nsv
               call initfield<<<griddim,blockdim>>>(svp_d(:, :, :, n), 0., ihc, jhc, khc)
               call checkCUDA( cudaGetLastError(), 'initfield svp_d' )
            end do
         end if

      end subroutine resetTendenciesDevice

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

      subroutine pull_ql0
         implicit none
         if (pulled(F_QL0)) return
         ql0 = ql0_d
         pulled(F_QL0) = .true.
      end subroutine pull_ql0

      subroutine pull_ql0h
         implicit none
         if (pulled(F_QL0H)) return
         ql0h = ql0h_d
         pulled(F_QL0H) = .true.
      end subroutine pull_ql0h

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
      subroutine updateHostForFielddump
         implicit none

         call pull_u0
         call pull_v0
         call pull_w0

         if (fielddump_wants('p0')) call pull_pres0
         if (ltempeq .and. fielddump_wants('th')) call pull_thl0
         if (lmoist  .and. fielddump_wants('qt')) call pull_qt0
         if (ltempeq .and. fielddump_wants('ql')) call pull_ql0

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

#if defined(UDALES_DEBUG)
         call assertHostMatchesDevice('fielddump')
#endif

      end subroutine updateHostForFielddump

      !> Bring down what statsdump is about to read.
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

#if defined(UDALES_DEBUG)
         call assertHostMatchesDevice('statsdump')
#endif

      end subroutine updateHostForStatsdump

      !> Bring down what writerestartfiles is about to write.
      subroutine updateHostForRestart
         implicit none

         call invalidateHostFieldPulls

         call pull_u0
         call pull_v0
         call pull_w0
         call pull_pres0
         call pull_ekm

         if (ltempeq) then
            call pull_thl0
            call pull_ql0
            call pull_ql0h
         end if
         if (lmoist)  call pull_qt0
         if (nsv > 0) call pull_sv0
         if (loneeqn) call pull_e120

#if defined(UDALES_DEBUG)
         call assertHostMatchesDevice('restart')
#endif

      end subroutine updateHostForRestart

      !> Send the hydrostatic column diagfld just rebuilt up to the device.
      !!
      !! Four columns of ke-kb+2 doubles, which is nothing, and they are what
      !! the saturation and virtual-temperature kernels read. Uploaded from
      !! inside diagfld rather than from updateDevice because thermodynamics'
      !! first saturation pass runs before diagfld in the same call and has to
      !! see the previous call's profiles, exactly as the host does.
      subroutine updateThermoProfilesDevice
         implicit none

         if (.not. allocated(presf_d)) return

         presf_d = presf
         presh_d = presh
         exnf_d  = exnf
         exnh_d  = exnh

      end subroutine updateThermoProfilesDevice

      !> Mirror the driver inlet planes that boundary_device reads, and fill them.
      !!
      !! Allocated from initCUDA, which runs after initdriver, so the host
      !! shapes are already known. The optional groups follow the host: a run
      !! can drive momentum from a file without driving temperature, moisture
      !! or the scalars, and then those planes do not exist on either side.
      !!
      !! The fill is part of allocating rather than a second call at the call
      !! site because there is no point in the run at which these may be read
      !! unfilled, and one caller forgetting is not a state worth being able
      !! to represent. By the time initCUDA runs the host planes are already
      !! complete: initdriver has read the driver file and the boundary call
      !! in the initialisation sequence has run drivergen over it at rk3step
      !! 0, which is the step that writes the m planes as well as the 0 ones.
      !!
      !! Leaving the fill to the loop does not work. boundary_device calls
      !! updateDriverPlanesDevice only on the stage drivergen itself runs on,
      !! the last of the three, so the first two stages of the first step
      !! would read whatever the allocation happened to land on. updateDevice
      !! uploads u0driver_d for bcpup and nothing else, so u0 was the one
      !! inlet field that came out right: case 452 left its first step with a
      !! divergence of 4.33 against 2.7E-13 for the same case on the host, a
      !! 40 K temperature anomaly at the inlet, and never recovered.
      subroutine allocDriverPlanesDevice
         implicit none

         ! initdriver allocates the host planes under idriver .and. ibrank,
         ! and xmi_driver_device only runs under ibrank too. On a rank that
         ! does not own the inlet there is nothing to mirror and nothing that
         ! would read the mirror.
         if (.not. allocated(u0driver)) return

         allocate(umdriver_d(jb-jh:je+jh,kb-kh:ke+kh))
         allocate(v0driver_d(jb-jh:je+jh,kb-kh:ke+kh))
         allocate(vmdriver_d(jb-jh:je+jh,kb-kh:ke+kh))
         allocate(w0driver_d(jb-jh:je+jh,kb-kh:ke+kh))
         allocate(wmdriver_d(jb-jh:je+jh,kb-kh:ke+kh))

         if (allocated(thl0driver)) then
            allocate(thl0driver_d(jb-jh:je+jh,kb-kh:ke+kh))
            allocate(thlmdriver_d(jb-jh:je+jh,kb-kh:ke+kh))
         end if
         if (allocated(qt0driver)) then
            allocate(qt0driver_d(jb-jh:je+jh,kb-kh:ke+kh))
            allocate(qtmdriver_d(jb-jh:je+jh,kb-kh:ke+kh))
         end if
         if (allocated(sv0driver)) then
            allocate(sv0driver_d(jb-jhc:je+jhc,kb-khc:ke+khc,1:nsv))
            allocate(svmdriver_d(jb-jhc:je+jhc,kb-khc:ke+khc,1:nsv))
         end if

         call updateDriverPlanesDevice

      end subroutine allocDriverPlanesDevice

      subroutine deallocDriverPlanesDevice
         implicit none

         ! Nothing was allocated on a rank that does not own the inlet.
         if (.not. allocated(umdriver_d)) return

         deallocate(umdriver_d, v0driver_d, vmdriver_d, w0driver_d, wmdriver_d)
         if (allocated(thl0driver_d)) deallocate(thl0driver_d, thlmdriver_d)
         if (allocated(qt0driver_d))  deallocate(qt0driver_d, qtmdriver_d)
         if (allocated(sv0driver_d))  deallocate(sv0driver_d, svmdriver_d)

      end subroutine deallocDriverPlanesDevice

      !> Send the driver inlet planes up, after drivergen has rewritten them.
      !!
      !! Called from boundary_device at the point drivergen returns, which is
      !! the only thing that writes them, and from allocDriverPlanesDevice for
      !! the first fill. Twelve j-k planes is a few hundred kilobytes against
      !! the hundred-odd megabytes a prognostic field costs, and it happens on
      !! the first and last stage of a step rather than on all three - so the
      !! driver inlet stays on the device for what it is worth to keep it there.
      subroutine updateDriverPlanesDevice
         implicit none

         ! The mirror only exists on the ranks that own the inlet.
         if (.not. allocated(umdriver_d)) return

         u0driver_d = u0driver
         umdriver_d = umdriver
         v0driver_d = v0driver
         vmdriver_d = vmdriver
         w0driver_d = w0driver
         wmdriver_d = wmdriver

         if (allocated(thl0driver_d)) then
            thl0driver_d = thl0driver
            thlmdriver_d = thlmdriver
         end if
         if (allocated(qt0driver_d)) then
            qt0driver_d = qt0driver
            qtmdriver_d = qtmdriver
         end if
         if (allocated(sv0driver_d)) then
            sv0driver_d = sv0driver
            svmdriver_d = svmdriver
         end if

      end subroutine updateDriverPlanesDevice

      !> Bring down the recycle plane writedriverfile is about to write.
      !!
      !! The plane, not the fields. Fetching u0, v0, w0, thl0, qt0 and sv0
      !! whole is what this did first, and on a 256^3 generation run with
      !! three scalars that is 1.1 GB a record - 3.6 seconds over a run that
      !! records thirty-six of them, which is more than the whole of boundary
      !! costs. Each plane is gathered on the device into a contiguous buffer
      !! instead and crosses as half a megabyte.
      !!
      !! Called from inside writedriverfile's own guards, the way
      !! updateHostForRestart is - but forced, rather than routed through the
      !! pull routines. It runs from the middle of boundary_device, after the
      !! top conditions have been applied on the device and before the lateral
      !! ones, because that is where the host path read these fields. The
      !! bitmap has no way to say "current as of here", so the flags are left
      !! alone and boundary_device clears all of them when it finishes.
      !!
      !! The i indices are writedriverfile's: it takes u from the recycle
      !! plane itself and everything else from the plane below it.
      subroutine updateHostForDriverDump
         use modinletdata, only : irecydriver
         implicit none
         integer :: i, n

         i = irecydriver
         call gatherPlane_u0(i)
         u0(i, :, :) = plane_stage

         i = irecydriver - 1
         call gatherPlane_v0(i)
         v0(i, :, :) = plane_stage
         call gatherPlane_w0(i)
         w0(i, :, :) = plane_stage

         if (ltempeq) then
            call gatherPlane_thl0(i)
            thl0(i, :, :) = plane_stage
         end if
         if (lmoist) then
            call gatherPlane_qt0(i)
            qt0(i, :, :) = plane_stage
         end if
         do n = 1, nsv
            call gatherPlane_sv0(i, n)
            sv0(i, :, :, n) = planec_stage
         end do

      end subroutine updateHostForDriverDump

      ! One gather per field, for the reason the pull routines are one per
      ! field: a device dummy argument would put a descriptor between the
      ! kernel and the array.

      subroutine gatherPlane_u0(i)
         implicit none
         integer, intent(in) :: i
         integer :: j, k

         !$acc parallel loop collapse(2) default(present)
         do k = kb - kh, ke + kh
            do j = jb - jh, je + jh
               plane_d(j, k) = u0_d(i, j, k)
            end do
         end do
         !$acc end parallel loop
         plane_stage = plane_d
      end subroutine gatherPlane_u0

      subroutine gatherPlane_v0(i)
         implicit none
         integer, intent(in) :: i
         integer :: j, k

         !$acc parallel loop collapse(2) default(present)
         do k = kb - kh, ke + kh
            do j = jb - jh, je + jh
               plane_d(j, k) = v0_d(i, j, k)
            end do
         end do
         !$acc end parallel loop
         plane_stage = plane_d
      end subroutine gatherPlane_v0

      subroutine gatherPlane_w0(i)
         implicit none
         integer, intent(in) :: i
         integer :: j, k

         !$acc parallel loop collapse(2) default(present)
         do k = kb - kh, ke + kh
            do j = jb - jh, je + jh
               plane_d(j, k) = w0_d(i, j, k)
            end do
         end do
         !$acc end parallel loop
         plane_stage = plane_d
      end subroutine gatherPlane_w0

      subroutine gatherPlane_thl0(i)
         implicit none
         integer, intent(in) :: i
         integer :: j, k

         !$acc parallel loop collapse(2) default(present)
         do k = kb - kh, ke + kh
            do j = jb - jh, je + jh
               plane_d(j, k) = thl0_d(i, j, k)
            end do
         end do
         !$acc end parallel loop
         plane_stage = plane_d
      end subroutine gatherPlane_thl0

      subroutine gatherPlane_qt0(i)
         implicit none
         integer, intent(in) :: i
         integer :: j, k

         !$acc parallel loop collapse(2) default(present)
         do k = kb - kh, ke + kh
            do j = jb - jh, je + jh
               plane_d(j, k) = qt0_d(i, j, k)
            end do
         end do
         !$acc end parallel loop
         plane_stage = plane_d
      end subroutine gatherPlane_qt0

      subroutine gatherPlane_sv0(i, n)
         implicit none
         integer, intent(in) :: i, n
         integer :: j, k

         !$acc parallel loop collapse(2) default(present)
         do k = kb - khc, ke + khc
            do j = jb - jhc, je + jhc
               planec_d(j, k) = sv0_d(i, j, k, n)
            end do
         end do
         !$acc end parallel loop
         planec_stage = planec_d
      end subroutine gatherPlane_sv0

      !> The device has written past what the host holds.
      !!
      !! boundary_device applies the top, bottom and lateral conditions on the
      !! device, and it runs after fielddump and statsdump have already pulled
      !! their fields down. Those copies were current when they were made and
      !! are not any more - they are missing the boundary planes. Clearing the
      !! bitmap is what makes the restart fetch them again rather than skip
      !! them as already pulled, and skipping them is exactly the failure the
      !! bitmap would otherwise introduce: a restart written on that step would
      !! record velocities with no boundary condition in them.
      !!
      !! All of them, not only the fifteen boundary_device writes. A list of
      !! exactly those fifteen would spare pres0, ekm and ekh a second fetch on
      !! the steps statsdump samples - about a third of a second on a 256^3 run
      !! of five thousand steps - at the price of a list that has to be revised
      !! every time a branch is added to boundary, with a stale device field
      !! and no error message as the failure mode. The blanket clear is the
      !! statement that the device has moved on, which is what actually
      !! happened.
      subroutine invalidateHostFieldPulls
         implicit none
         pulled = .false.
      end subroutine invalidateHostFieldPulls

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
      !! Placement is load-bearing. It has to run before boundary, which
      !! applies the top and lateral planes on the device and so makes host and
      !! device legitimately differ until the pull that follows it.
      !!
      !! Debug builds only, and it copies everything down a second time, so it
      !! is not cheap - but the GPU parity cases are 64^3 and run in a second.
      subroutine assertHostMatchesDevice(label)
         use modmpi, only : myid
         implicit none

         character(len=*), intent(in) :: label
         integer :: n

         if (.not. lchecks_armed) return

         select case (label)

         case ('fielddump')
            ! div is formed from these before fieldvars is consulted at all.
            call check3(label, 'u0', u0, u0_d)
            call check3(label, 'v0', v0, v0_d)
            call check3(label, 'w0', w0, w0_d)

            if (fielddump_wants('p0')) call check3(label, 'pres0', pres0, pres0_d)
            if (ltempeq .and. fielddump_wants('th')) call check3(label, 'thl0', thl0, thl0_d)
            if (lmoist  .and. fielddump_wants('qt')) call check3(label, 'qt0',  qt0,  qt0_d)
            if (ltempeq .and. fielddump_wants('ql')) call check3(label, 'ql0',  ql0,  ql0_d)

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

         case ('restart')
            ! Exactly what updateHostForRestart fetches, in the same order and
            ! under the same guards. The two lists are the same list written
            ! twice, which is the one thing wrong with this check - a field
            ! added to the pull without being added here is not caught.
            call check3(label, 'u0', u0, u0_d)
            call check3(label, 'v0', v0, v0_d)
            call check3(label, 'w0', w0, w0_d)
            call check3(label, 'pres0', pres0, pres0_d)
            call check3(label, 'ekm', ekm, ekm_d)

            if (ltempeq) then
               call check3(label, 'thl0', thl0, thl0_d)
               call check3(label, 'ql0',  ql0,  ql0_d)
               call check3(label, 'ql0h', ql0h, ql0h_d)
            end if
            if (lmoist) call check3(label, 'qt0', qt0, qt0_d)
            if (nsv > 0) then
               do n = 1, nsv
                  call check3(label, 'sv0', sv0(:,:,:,n), sv0_d(:,:,:,n))
               end do
            end if
            if (loneeqn) call check3(label, 'e120', e120, e120_d)

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

      !> avexy_ibm_device for a field carrying the scalar halo.
      !!
      !! sv0 is the only one, and it needs its own routine rather than a wider
      !! dummy on the one above because an explicit-shape device dummy takes
      !! the address of the first element: pass sv0_d to a dummy declared with
      !! the momentum halo and, whenever ihc differs from ih, every index is
      !! wrong by a margin that grows with k.
      !!
      !! The k loop still stops at ke+kh, not ke+khc, which is what the host
      !! does. There the truncation is an accident of sequence association -
      !! avexy_ibm's aver dummy is declared kb:ke+kh and the caller hands it a
      !! kb:ke+khc slice - and the tail level of sv0av keeps the zero diagfld
      !! wrote before the reduction. Here it is written down.
      subroutine avexy_ibm_c_device(aver, var, kzb, II, IIs, lnan)
         use modmpi, only : avexy_ibm_finish
         implicit none

         integer,         intent(in)  :: kzb
         real   , device, intent(in)  :: var(ib-ihc:ie+ihc,jb-jhc:je+jhc,kzb:ke+khc)
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

      end subroutine avexy_ibm_c_device

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
