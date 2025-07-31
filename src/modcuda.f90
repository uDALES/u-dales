module modcuda
#if defined(_GPU)        
   use cudafor
   use modglobal,      only: itot, ib, ie, jb, je, kb, ke, ih, jh, kh, ihc, jhc, khc, &
                             dx2, dxi, dx2i, dxi5, dxiq, dy2, dyi, dy2i, dyi5, dyiq, dxf, dxhi, &
                             dzf, dzf2, dzfi, dzfi5, dzfiq, dzh, dzhi, dzh2i, dzhiq, &
                             dzfc, dzfci, dzhci, dxfc, dxfci, dxhci, delta, &
                             ltempeq, lmoist, nsv, lles, lbuoyancy, &
                             BCxm, BCxm_periodic, BCym, BCym_periodic, &
                             BCtopm, BCtopm_freeslip, BCtopm_pressure, BCtopm_noslip, &
                             iadv_sv, iadv_thl, iadv_kappa, iadv_upw, &
                             xlen, ds, xh, &
                             pi, eps1, numol, prandtlmoli, prandtlturb, grav, fkar2
   use modfields,      only: u0, v0, w0, pres0, e120, thl0, thl0c, qt0, sv0, &
                             up, vp, wp, e12p, thlp, thlpc, qtp, svp, &
                             e12m, &
                             tau_x, tau_y, tau_z, thl_flux, &
                             u0av, dthvdz, ug
   use modsubgriddata, only: lsmagorinsky, lvreman, loneeqn, ldelta, lbuoycorr, &
                             ekm, ekh, &
                             sbshr, sbbuo, sbdiss, zlt, damp, csz, &
                             cn, cm, ch1, ch2, ce1, ce2, dampmin, prandtli, c_vreman
   use modsurfdata,    only: thvs
   use decomp_2d,      only: zstart
   implicit none
   save

   type(dim3) :: griddim, blockdim

   integer, device :: itot_d, ib_d, ie_d, jb_d, je_d, kb_d, ke_d, ih_d, jh_d, kh_d, &
                      BCxm_d, BCxm_periodic_d, BCym_d, BCym_periodic_d, &
                      BCtopm_d, BCtopm_freeslip_d, BCtopm_pressure_d, BCtopm_noslip_d
   
   logical, device :: ltempeq_d, lles_d, lsmagorinsky_d, lvreman_d, loneeqn_d, ldelta_d, lbuoyancy_d, lbuoycorr_d

   real,    device :: dx2_d, dxi_d, dx2i_d, dxi5_d, dxiq_d, dy2_d, dyi_d, dy2i_d, dyi5_d, dyiq_d, &
                      xlen_d, ds_d, &
                      eps1_d, pi_d, numol_d, prandtlmoli_d, prandtlturb_d, prandtli_d, grav_d, dampmin_d, c_vreman_d, fkar2_d, &
                      cn_d, cm_d, ch1_d, ch2_d, ce1_d, ce2_d, &
                      thvs_d

   integer, device, dimension(3) :: zstart_d

   real, device, allocatable :: dzf_d(:), dzf2_d(:), dzfi_d(:), dzfi5_d(:), dzfiq_d(:), dzh_d(:), dzhi_d(:), dzh2i_d(:), dzhiq_d(:), &
                                dzfc_d(:), dzfci_d(:), dzhci_d(:), dxfc_d(:), dxfci_d(:), dxhci_d(:), &
                                dxf_d(:), dxhi_d(:), &
                                xh_d(:), u0av_d(:), ug_d(:)

   real, device, allocatable :: delta_d(:, :), csz_d(:,:)

   real, device, allocatable :: u0_d(:,:,:), v0_d(:,:,:), w0_d(:,:,:), pres0_d(:,:,:), e120_d(:,:,:), thl0_d(:,:,:), thl0c_d(:,:,:), qt0_d(:,:,:), sv0_d(:,:,:,:)
   real, device, allocatable :: up_d(:,:,:), vp_d(:,:,:), wp_d(:,:,:), e12p_d(:,:,:), thlp_d(:,:,:), thlpc_d(:,:,:), qtp_d(:,:,:), svp_d(:,:,:,:)
   real, device, allocatable :: e12m_d(:,:,:)
   real, device, allocatable :: tau_x_d(:,:,:), tau_y_d(:,:,:), tau_z_d(:,:,:), thl_flux_d(:,:,:)
   real, device, allocatable :: dthvdz_d(:,:,:)
   real, device, allocatable :: ekm_d(:,:,:), ekh_d(:,:,:), sbshr_d(:,:,:), sbbuo_d(:,:,:), sbdiss_d(:,:,:), zlt_d(:,:,:), damp_d(:,:,:)

   real, device, allocatable :: dumu_d(:,:,:), duml_d(:,:,:)

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

         itot_d = itot
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
         allocate (dzhiq_d(kb:ke + kh))
         dzf_d   = dzf
         dzf2_d  = dzf2
         dzfi_d  = dzfi
         dzfi5_d = dzfi5
         dzfiq_d = dzfiq
         dzh_d   = dzh
         dzhi_d  = dzhi
         dzh2i_d = dzh2i
         dzhiq_d = dzhiq

         allocate(delta_d(ib-ih:itot+ih, kb:ke + kh))
         delta_d = delta

         BCxm_d            = BCxm
         BCxm_periodic_d   = BCxm_periodic
         BCym_d            = BCym
         BCym_periodic_d   = BCym_periodic
         BCtopm_d          = BCtopm
         BCtopm_freeslip_d = BCtopm_freeslip
         BCtopm_pressure_d = BCtopm_pressure
         BCtopm_noslip_d   = BCtopm_noslip

         allocate(u0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         allocate(v0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         allocate(w0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         allocate(pres0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         
         allocate(up_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
         allocate(vp_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
         allocate(wp_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))

         allocate(tau_x_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         allocate(tau_y_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         allocate(tau_z_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))

         if (loneeqn) then
            allocate(e120_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
            allocate(e12p_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
            allocate(e12m_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
            allocate(sbshr_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
            allocate(sbbuo_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
            allocate(sbdiss_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
            allocate(zlt_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
         end if

         if (ltempeq) then
            allocate(thl0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
            allocate(thlp_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
            if (iadv_thl == iadv_kappa) then
               allocate(thl0c_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb-khc:ke+khc))
               allocate(thlpc_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
            end if
            allocate(thl_flux_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         end if

         if (lmoist) then
            allocate(qt0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
            allocate(qtp_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
         end if

         if (nsv>0) then
            allocate(sv0_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb-khc:ke+khc,nsv))
            allocate(svp_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc,nsv))
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
         
         eps1_d = eps1
         pi_d     = pi
         xlen_d   = xlen
         ds_d     = ds
         zstart_d = zstart
         allocate (xh_d(ib:itot+ih))
         xh_d = xh

         allocate(u0av_d(kb:ke+kh))

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

         allocate(ug_d(kb:ke+kh))
         ug_d = ug

      end subroutine initCUDA

      subroutine exitCUDA
         implicit none
         deallocate(dxf_d, dxhi_d, dzf_d, dzf2_d, dzfi_d, dzfi5_d, dzfiq_d, dzh_d, dzhi_d, dzh2i_d, dzhiq_d, delta_d)
         deallocate(u0_d, v0_d, w0_d, pres0_d)
         deallocate(up_d, vp_d, wp_d)
         deallocate(tau_x_d, tau_y_d, tau_z_d)
         if (loneeqn) deallocate(e120_d, e12p_d, e12m_d, sbshr_d, sbbuo_d, sbdiss_d, zlt_d)
         if (ltempeq) then
            deallocate(thl0_d, thlp_d, thl_flux_d)
            if (iadv_thl == iadv_kappa) deallocate(thl0c_d, thlpc_d)
         end if
         if (lmoist) deallocate(qt0_d, qtp_d)
         if (nsv>0) deallocate(sv0_d, svp_d)
         if (any(iadv_sv(1:nsv) == iadv_kappa) .or. (iadv_thl == iadv_kappa)) then
            deallocate(dumu_d, duml_d, dxhci_d, dxfc_d, dzhci_d, dzfc_d)
         end if
         if (any(iadv_sv(1:nsv) == iadv_kappa) .or. any(iadv_sv(1:nsv) == iadv_upw) .or. (iadv_thl == iadv_kappa)) then
            deallocate(dxfci_d, dzfci_d)
         end if
         deallocate(xh_d, u0av_d)
         deallocate(ekm_d, ekh_d)
         if (lsmagorinsky .or. loneeqn) then
            deallocate(damp_d)
         end if
         if (lsmagorinsky) deallocate(csz_d)
         deallocate(dthvdz_d)
         deallocate(ug_d)
      end subroutine exitCUDA

      subroutine updateDevice
         implicit none
         integer :: n

         u0_d = u0
         v0_d = v0
         w0_d = w0
         pres0_d = pres0

         call initfield<<<griddim,blockdim>>>(up_d, 0.)
         call checkCUDA( cudaGetLastError(), 'initfield up_d' )

         call initfield<<<griddim,blockdim>>>(vp_d, 0.)
         call checkCUDA( cudaGetLastError(), 'initfield vp_d' )

         call initfield<<<griddim,blockdim>>>(wp_d, 0.)
         call checkCUDA( cudaGetLastError(), 'initfield wp_d' )

         tau_x_d = tau_x
         tau_y_d = tau_y
         tau_z_d = tau_z

         if (loneeqn) then
            e12m_d = e12m
            e120_d = e120
            call initfield<<<griddim,blockdim>>>(e12p_d, 0.)
            call checkCUDA( cudaGetLastError(), 'initfield e12p_d' )
         end if

         if (ltempeq) then
            thl0_d = thl0
            call initfield<<<griddim,blockdim>>>(thlp_d, 0.)
            call checkCUDA( cudaGetLastError(), 'initfield thlp_d' )

            if (iadv_thl == iadv_kappa) then
               thl0c_d = thl0c
               call initfield<<<griddim,blockdim>>>(thlpc_d, 0.)
               call checkCUDA( cudaGetLastError(), 'initfield thlpc_d' )
            end if

            thl_flux_d = thl_flux
         end if

         if (lmoist) then
            qt0_d = qt0
            call initfield<<<griddim,blockdim>>>(qtp_d, 0.)
            call checkCUDA( cudaGetLastError(), 'initfield qtp_d' )
         end if

         if (nsv>0) then
            sv0_d = sv0
            do n = 1, nsv
               call initfield<<<griddim,blockdim>>>(svp_d(:, :, :, n), 0.)
               call checkCUDA( cudaGetLastError(), 'initfield svp_d' )
            end do
         end if

         u0av_d = u0av
         dthvdz_d = dthvdz
      end subroutine updateDevice

      subroutine updateHost
         implicit none
         up = up_d
         vp = vp_d
         wp = wp_d
         tau_x = tau_x_d
         tau_y = tau_y_d
         tau_z = tau_z_d
         if (loneeqn) then
            e12m = e12m_d
            e120 = e120_d
            e12p = e12p_d
         end if
         if (ltempeq) then
            thlp = thlp_d
            if (iadv_thl == iadv_kappa) thlpc = thlpc_d
            thl_flux = thl_flux_d
         end if
         if (lmoist) then
            qtp = qtp_d
         end if
         if (nsv>0) then
            svp = svp_d
         end if
         ekm = ekm_d
         ekh = ekh_d
      end subroutine updateHost

      subroutine checkCUDA(istat, kernelname)
         implicit none
         integer, intent(in)          :: istat
         character(len=*), intent(in) :: kernelname
         if(istat /= cudaSuccess) then
            write(*,*) "Error in ", trim(kernelname), ": ", cudaGetErrorString(istat)
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

      attributes(global) subroutine initfield(var, varvalue)
         implicit none
         real, dimension(ib_d-ih_d:ie_d+ih_d, jb_d-jh_d:je_d+jh_d, kb_d:ke_d+kh_d), intent(inout) :: var
         real, value, intent(in) :: varvalue

         integer :: tidx, tidy, tidz, stridex, stridey, stridez
         integer :: i, j, k

         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

         do i = tidx-ih_d, ie_d+ih_d, stridex
            do j = tidy-jh_d, je_d+jh_d, stridey
               do k = tidz, ke_d+kh_d, stridez
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
