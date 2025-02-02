module modcuda
#if defined(_GPU)        
   use cudafor
   use modglobal,      only: itot, ib, ie, jb, je, kb, ke, ih, jh, kh, ihc, jhc, khc, &
                             dxi, dx2i, dxi5, dxiq, dyi, dy2i, dyi5, dyiq, dzf, dzfi, dzfi5, dzhi, dzh2i, dzhiq, &
                             dzfc, dzfci, dzhci, dxfc, dxfci, dxhci, delta, &
                             ltempeq, lmoist, nsv, lles, &
                             iadv_sv, iadv_thl, iadv_kappa, iadv_upw, &
                             xlen, ds, xh, &
                             rk3step, dt, &
                             pi, eps1, numol, prandtlmoli, grav
   use modfields,      only: u0, v0, w0, pres0, e120, thl0, thl0c, qt0, sv0, &
                             up, vp, wp, e12p, thlp, thlpc, qtp, svp, &
                             u0av, dthvdz
   use modsubgriddata, only: lsmagorinsky, lvreman, loneeqn, ldelta, &
                             ekm, ekh, &
                             sbshr, sbbuo, sbdiss, zlt, damp, csz, &
                             cn, cm, ch1, ch2, ce1, ce2, dampmin, prandtli
   use modsurfdata,    only: thvs
   use decomp_2d,      only: zstart
   implicit none
   save

   type(dim3) :: griddim, blockdim

   integer, device :: itot_d, ib_d, ie_d, jb_d, je_d, kb_d, ke_d, ih_d, jh_d, kh_d, rk3step_d
   logical, device :: lles_d, lsmagorinsky_d, lvreman_d, loneeqn_d, ldelta_d
   real,    device :: dxi_d, dx2i_d, dxi5_d, dxiq_d, dyi_d, dy2i_d, dyi5_d, dyiq_d, &
                      xlen_d, ds_d, &
                      eps1_d, pi_d, numol_d, prandtlmoli_d, prandtli_d, grav_d, dampmin_d, &
                      thvs_d, cn_d, cm_d, ch1_d, ch2_d, ce1_d, ce2_d, &
                      dt_d

   integer, device, dimension(3) :: zstart_d

   real, device, allocatable :: dzf_d(:), dzfi_d(:), dzhi_d(:), dzh2i_d(:), dzfi5_d(:), dzhiq_d(:), &
                                dzfc_d(:), dzfci_d(:), dzhci_d(:), dxfc_d(:), dxfci_d(:), dxhci_d(:), &
                                xh_d(:), u0av_d(:)

   real, device, allocatable :: delta_d(:, :), csz_d(:,:)

   real, device, allocatable :: u0_d(:,:,:), v0_d(:,:,:), w0_d(:,:,:), pres0_d(:,:,:), e120_d(:,:,:), thl0_d(:,:,:), thl0c_d(:,:,:), qt0_d(:,:,:), sv0_d(:,:,:,:)
   real, device, allocatable :: up_d(:,:,:), vp_d(:,:,:), wp_d(:,:,:), e12p_d(:,:,:), thlp_d(:,:,:), thlpc_d(:,:,:), qtp_d(:,:,:), svp_d(:,:,:,:)
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

         threadnumx = props%warpsize/4
         threadnumy = props%warpsize/4
         threadnumz = 1
         if (threadnumx*threadnumy*threadnumz > props%maxThreadsPerBlock) then
            write(*,*) "Incorrect block dimension configuration."
            stop 1
         end if

         blocknumx = min( max( floor(real( (ie - ib + 1)/threadnumx ) ), 1 ), props%maxGridSize(1) )
         blocknumy = min( max( floor(real( (je - jb + 1)/threadnumy ) ), 1 ), props%maxGridSize(2) )
         blocknumz = min( max( floor(real( (ke - kb + 1)/threadnumz ) ), 1 ), props%maxGridSize(3) )

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
         dxi5_d = dxi5
         dyi5_d = dyi5
         dxiq_d = dxiq
         dyiq_d = dyiq
         dxi_d = dxi
         dyi_d = dyi
         allocate (dzf_d(kb - kh:ke + kh))
         allocate (dzhi_d(kb:ke + kh))
         allocate (dzfi5_d(kb - kh:ke + kh))
         allocate (dzhiq_d(kb:ke + kh))
         dzf_d = dzf
         dzhi_d = dzhi
         dzfi5_d = dzfi5
         dzhiq_d = dzhiq

         allocate(u0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         allocate(v0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         allocate(w0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         allocate(pres0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         
         allocate(up_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
         allocate(vp_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
         allocate(wp_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))

         if (loneeqn) then
            allocate(e120_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
            allocate(e12p_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
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
         end if

         if (lmoist) then
            allocate(qt0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
            allocate(qtp_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
         end if

         if (nsv>0) then
            allocate(sv0_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb-khc:ke+khc,nsv))
            allocate(svp_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc,nsv))
         end if

         if (any(iadv_sv(1:nsv) == iadv_kappa) .or. any(iadv_sv(1:nsv) == iadv_upw) .or. (iadv_thl == iadv_kappa)) then
            allocate(dumu_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
            allocate(duml_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
            allocate (dxhci_d(ib - 1:itot + ihc))
            allocate (dxfc_d(ib - ihc:itot + ihc))
            allocate (dxfci_d(ib - ihc:itot + ihc))
            allocate (dzhci_d(kb - 1:ke + khc))
            allocate (dzfc_d(kb - khc:ke + khc))
            allocate (dzfci_d(kb - khc:ke + khc))
            dxhci_d = dxhci
            dxfc_d  = dxfc
            dxfci_d = dxfci
            dzhci_d = dzhci
            dzfc_d  = dzfc
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

         lles_d = lles
         lsmagorinsky_d = lsmagorinsky
         lvreman_d = lvreman
         loneeqn_d = loneeqn
         ldelta_d  = ldelta
         allocate(ekm_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         allocate(ekh_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
         dx2i_d = dx2i
         allocate (dzfi_d(kb - kh:ke + kh))
         dzfi_d = dzfi
         numol_d = numol
         prandtlmoli_d = prandtlmoli
         
         dy2i_d = dy2i
         allocate (dzh2i_d(kb:ke + kh))
         dzh2i_d = dzh2i
         grav_d = grav
         thvs_d = thvs
         allocate(dthvdz_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
         cm_d = cm
         cn_d = cn
         ch1_d = ch1
         ch2_d = ch2
         ce1_d = ce1
         ce2_d = ce2
         prandtli_d = prandtli

         if (lsmagorinsky .or. loneeqn) then
            allocate(damp_d(ib:ie,jb:je,kb:ke))
            dampmin_d = dampmin
         end if

         allocate(delta_d(ib-ih:itot+ih, kb:ke + kh))
         delta_d = delta

         if (lsmagorinsky) then
            allocate(csz_d(ib-ih:ie+ih,kb:ke+kh))
            csz_d = csz
         end if

      end subroutine initCUDA

      subroutine exitCUDA
         implicit none
         deallocate(dzf_d, dzfi_d, dzh2i_d, dzhi_d, dzfi5_d, dzhiq_d, delta_d)
         deallocate(u0_d, v0_d, w0_d, pres0_d)
         deallocate(up_d, vp_d, wp_d)
         if (loneeqn) deallocate(e120_d, e12p_d, sbshr_d, sbbuo_d, sbdiss_d, zlt_d)
         if (ltempeq) then
            deallocate(thl0_d, thlp_d)
            if (iadv_thl == iadv_kappa) deallocate(thl0c_d, thlpc_d)
         end if
         if (lmoist) deallocate(qt0_d, qtp_d)
         if (nsv>0) deallocate(sv0_d, svp_d)
         if (any(iadv_sv(1:nsv) == iadv_kappa) .or. any(iadv_sv(1:nsv) == iadv_upw) .or. (iadv_thl == iadv_kappa)) then
            deallocate(dumu_d, duml_d, dxhci_d, dxfc_d, dxfci_d, dzhci, dzfc, dzfci)
         end if
         deallocate(xh_d, u0av_d)
         deallocate(ekm_d, ekh_d)
         if (lsmagorinsky .or. loneeqn) then
            deallocate(damp_d)
         end if
         if (lsmagorinsky) deallocate(csz_d)
         deallocate(dthvdz_d)
      end subroutine exitCUDA

      subroutine updateDevice
         implicit none
         dt_d = dt
         rk3step_d = rk3step
         u0_d = u0
         v0_d = v0
         w0_d = w0
         pres0_d = pres0
         up_d = up
         vp_d = vp
         wp_d = wp
         if (loneeqn) then
            e120_d = e120
            e12p_d = e12p
         end if
         if (ltempeq) then
            thl0_d = thl0
            thlp_d = thlp
            if (iadv_thl == iadv_kappa) then
               thl0c_d = thl0c
               thlpc_d = thlpc
            end if
         end if
         if (lmoist) then
            qt0_d = qt0
            qtp_d = qtp
         end if
         if (nsv>0) then
            sv0_d = sv0
            svp_d = svp
         end if
         u0av_d = u0av
         dthvdz_d = dthvdz
      end subroutine updateDevice

      subroutine updateHost
         implicit none
     !    u0 = u0_d
     !    v0 = v0_d
     !    w0 = w0_d
     !    pres0 = pres0_d
         up = up_d
         vp = vp_d
         wp = wp_d
         if (loneeqn) then
     !       e120 = e120_d
            e12p = e12p_d
         end if
         if (ltempeq) then
     !       thl0 = thl0_d
            thlp = thlp_d
            if (iadv_thl == iadv_kappa) thlpc = thlpc_d
         end if
         if (lmoist) then
     !       qt0 = qt0_d
            qtp = qtp_d
         end if
         if (nsv>0) then
     !       sv0 = sv0_d
            svp = svp_d
         end if
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
      

      !> Advection at cell center central difference
      attributes(global) subroutine advecc_2nd_cuda(hi, hj, hk, putin, putout)
         implicit none
         
         integer, value, intent(in) :: hi, hj, hk

         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d - hk:ke_d + hk), intent(in)    :: putin
         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d:ke_d + hk)     , intent(inout) :: putout
         
         integer :: tidx, tidy, tidz, stridex, stridey, stridez
         integer :: i, j, k, im, ip, jm, jp, km, kp
         
         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

         do i = tidx, ie_d, stridex
            im = i - 1
            ip = i + 1
            do j = tidy, je_d, stridey
               jm = j - 1
               jp = j + 1
               do k = tidz, ke_d, stridez
                  putout(i, j, k) = putout(i, j, k) - ( &
                                    ( &
                                    u0_d(ip, j, k)*(putin(ip, j, k) + putin(i, j, k)) &
                                  - u0_d(i, j, k)*(putin(im, j, k) + putin(i, j, k)) & ! d(uc)/dx
                                    )*dxi5_d &
                                  + ( & !
                                    v0_d(i, jp, k)*(putin(i, jp, k) + putin(i, j, k)) &
                                  - v0_d(i, j, k)*(putin(i, jm, k) + putin(i, j, k)) & ! d(vc)/dy
                                    )*dyi5_d)
               end do
            end do
         end do

         do i = tidx, ie_d, stridex
            im = i - 1
            ip = i + 1
            do j = tidy, je_d, stridey
               jm = j - 1
               jp = j + 1
               do k = tidz, ke_d, stridez
                  km = k - 1
                  kp = k + 1
                  putout(i, j, k) = putout(i, j, k) - ( &
                                    w0_d(i, j, kp)*(putin(i, j, kp)*dzf_d(k) + putin(i, j, k)*dzf_d(kp))*dzhi_d(kp) &
                                  - w0_d(i, j, k)*(putin(i, j, km)*dzf_d(k) + putin(i, j, k)*dzf_d(km))*dzhi_d(k) &
                                    )*dzfi5_d(k)
               end do
            end do
         end do

      end subroutine advecc_2nd_cuda
     

      !> Advection at the u point. central difference
      attributes(global) subroutine advecu_2nd_cuda(putin, putout)
         implicit none
             
         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d - kh_d:ke_d + kh_d), intent(in)    :: putin
         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d:ke_d + kh_d)       , intent(inout) :: putout

         integer :: tidx, tidy, tidz, stridex, stridey, stridez         
         integer :: i, j, k, im, ip, jm, jp, km, kp

         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

         do i = tidx, ie_d, stridex
            im = i - 1
            ip = i + 1
            do j = tidy, je_d, stridey
               jm = j - 1
               jp = j + 1
               do k = tidz, ke_d, stridez
                  km = k - 1
                  kp = k + 1
                  putout(i, j, k) = putout(i, j, k) - ( &
                                    ( &
                                    (putin(i, j, k) + putin(ip, j, k))*(u0_d(i, j, k) + u0_d(ip, j, k)) &
                                  - (putin(i, j, k) + putin(im, j, k))*(u0_d(i, j, k) + u0_d(im, j, k)) & ! d(uu)/dx
                                    )*dxiq_d &
                                  + ( &
                                    (putin(i, j, k) + putin(i, jp, k))*(v0_d(i, jp, k) + v0_d(im, jp, k)) &
                                  - (putin(i, j, k) + putin(i, jm, k))*(v0_d(i, j, k) + v0_d(im, j, k)) & ! d(vu)/dy
                                    )*dyiq_d) &
                                  - ((pres0_d(i, j, k) - pres0_d(im, j, k))*dxi_d) ! - dp/dx
               end do
            end do
         end do

         do k = tidz, ke_d, stridez
            km = k - 1
            kp = k + 1
            do j = tidy, je_d, stridey
               jm = j - 1
               jp = j + 1
               do i = tidx, ie_d, stridex
                  im = i - 1
                  ip = i + 1
                  putout(i, j, k) = putout(i, j, k) - ( &
                                    (putin(i, j, kp)*dzf_d(k) + putin(i, j, k)*dzf_d(kp))*dzhi_d(kp) &
                                  * (w0_d(i, j, kp) + w0_d(im, j, kp)) &
                                  - (putin(i, j, k)*dzf_d(km) + putin(i, j, km)*dzf_d(k))*dzhi_d(k) &
                                  * (w0_d(i, j, k) + w0_d(im, j, k)) &
                                    )*0.5*dzfi5_d(k)
               end do
            end do
         end do

      end subroutine advecu_2nd_cuda


      !> Advection at the v point. central difference
      attributes(global) subroutine advecv_2nd_cuda(putin, putout)
         implicit none
         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d - kh_d:ke_d + kh_d), intent(in)    :: putin
         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d:ke_d + kh_d)       , intent(inout) :: putout

         integer :: tidx, tidy, tidz, stridex, stridey, stridez
         integer :: i, j, k, im, ip, jm, jp, km, kp

         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

         do i = tidx, ie_d, stridex
            im = i - 1
            ip = i + 1
            do j = tidy, je_d, stridey
               jm = j - 1
               jp = j + 1
               do k = tidz, ke_d, stridez
                  km = k - 1
                  kp = k + 1
                  putout(i, j, k) = putout(i, j, k) - ( &
                                    ( &
                                    (u0_d(ip, j, k) + u0_d(ip, jm, k))*(putin(i, j, k) + putin(ip, j, k)) &
                                  - (u0_d(i, j, k)  + u0_d(i, jm, k)) *(putin(i, j, k) + putin(im, j, k)) & ! d(uv)/dx
                                    )*dxiq_d &
                                  + ( &
                                    (v0_d(i, jp, k) + v0_d(i, j, k))*(putin(i, j, k) + putin(i, jp, k)) &
                                  - (v0_d(i, jm, k) + v0_d(i, j, k))*(putin(i, j, k) + putin(i, jm, k)) & ! d(vv)/dy
                                    )*dyiq_d &
                                    ) &
                                  - ((pres0_d(i, j, k) - pres0_d(i, jm, k))*dyi_d) ! - dp/dy
               end do
            end do
         end do

         do k = tidz, ke_d, stridez
            km = k - 1
            kp = k + 1
            do j = tidy, je_d, stridey
               jm = j - 1
               jp = j + 1
               do i = tidx, ie_d, stridex
                  im = i - 1
                  ip = i + 1
                  putout(i, j, k) = putout(i, j, k) - ( &
                                    (w0_d(i, j, kp) + w0_d(i, jm, kp)) &
                                  * (putin(i, j, kp)*dzf_d(k) + putin(i, j, k)*dzf_d(kp))*dzhi_d(kp) &
                                  - (w0_d(i, j, k) + w0_d(i, jm, k)) &
                                  * (putin(i, j, km)*dzf_d(k) + putin(i, j, k)*dzf_d(km))*dzhi_d(k) &
                                    )*0.5*dzfi5_d(k)
               end do
            end do
         end do

      end subroutine advecv_2nd_cuda


      !> Advection at the w point. central difference
      attributes(global) subroutine advecw_2nd_cuda(putin, putout)
         implicit none
         
         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d - kh_d:ke_d + kh_d), intent(in)    :: putin
         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d:ke_d + kh_d)       , intent(inout) :: putout

         integer :: tidx, tidy, tidz, stridex, stridey, stridez
         integer :: i, j, k, im, ip, jm, jp, km, kp

         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

         do i = tidx, ie_d, stridex
            im = i - 1
            ip = i + 1
            do j = tidy, je_d, stridey
               jm = j - 1
               jp = j + 1
               do k = tidz, ke_d, stridez
                  km = k - 1
                  kp = k + 1
                  putout(i, j, k) = putout(i, j, k) - ( &
                                    ( &
                                    (putin(ip, j, k) + putin(i, j, k))*(dzf_d(km)*u0_d(ip, j, k) + dzf_d(k)*u0_d(ip, j, km)) &
                                  - (putin(i, j, k)  + putin(im, j, k))*(dzf_d(km)*u0_d(i, j, k) + dzf_d(k)*u0_d(i, j, km)) &
                                    )*dxiq_d*dzhi_d(k) & ! d(uw)/dx
                                  + ( &
                                    (putin(i, jp, k) + putin(i, j, k))*(dzf_d(km)*v0_d(i, jp, k) + dzf_d(k)*v0_d(i, jp, km)) &
                                  - (putin(i, j, k) + putin(i, jm, k))*(dzf_d(km)*v0_d(i, j, k) + dzf_d(k)*v0_d(i, j, km)) &
                                    )*dyiq_d*dzhi_d(k) & ! d(vw)/dy
                                  + ( &
                                    (putin(i, j, k) + putin(i, j, kp))*(w0_d(i, j, k) + w0_d(i, j, kp)) &
                                  - (putin(i, j, k) + putin(i, j, km))*(w0_d(i, j, k) + w0_d(i, j, km)) &
                                    )*dzhiq_d(k) & ! d(ww)/dz
                                    ) &
                                  - ((pres0_d(i, j, k) - pres0_d(i, j, km))*dzhi_d(k)) ! - dp/dz
               end do
            end do
         end do

      end subroutine advecw_2nd_cuda


      !> Advection at cell center through kappa scheme
      attributes(global) subroutine advecc_kappa_reset_cuda(hi, hj, hk)
         implicit none
         integer, value, intent(in) :: hi, hj, hk
         integer :: i, j, k, tidx, tidy, tidz, stridex, stridey, stridez
         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
         do k = tidz, ke_d + hk, stridez
            do j = tidy - hj, je_d + hj, stridey
               do i = tidx - hi, ie_d + hi, stridex
                  dumu_d(i,j,k) = 0.
                  duml_d(i,j,k) = 0.
               end do
            end do
         end do
      end subroutine advecc_kappa_reset_cuda

      ! -d(uc)/dx (stretched grid)
      attributes(global) subroutine advecc_kappa_ducdx_cuda(hi, hj, hk, var)
         implicit none
         integer, value, intent(in) :: hi, hj, hk
         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d - hk:ke_d + hk), intent(in)    :: var
         real    :: cf, d1, d2
         integer :: i, j, k, tidx, tidy, tidz, stridex, stridey, stridez

         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

         do k = tidz, ke_d, stridez
            do j = tidy, je_d, stridey
               do i = tidx, ie_d + 1, stridex
                  if (u0_d(i, j, k) > 0) then
                     d1 = (var(i - 1, j, k) - var(i - 2, j, k))*dxhci_d(i - 1)
                     d2 = (var(i, j, k) - var(i - 1, j, k))*dxhci_d(i)
                     cf = var(i - 1, j, k)
                  else
                     d1 = (var(i, j, k) - var(i + 1, j, k))*dxhci_d(i + 1)
                     d2 = (var(i - 1, j, k) - var(i, j, k))*dxhci_d(i)
                     cf = var(i, j, k)
                  end if
                  cf = cf + dxfc_d(i)*rlim_cuda(d1, d2)
                  dumu_d(i - 1, j, k) = -cf*u0_d(i, j, k)*dxfci_d(i - 1)
                  duml_d(i, j, k) = cf*u0_d(i, j, k)*dxfci_d(i)
               end do
            end do
         end do
      end subroutine advecc_kappa_ducdx_cuda

      ! -d(vc)/dy (no stretched grid)
      attributes(global) subroutine advecc_kappa_dvcdy_cuda(hi, hj, hk, var)
         implicit none
         integer, value, intent(in) :: hi, hj, hk
         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d - hk:ke_d + hk), intent(in)    :: var
         real    :: cf, d1, d2
         integer :: i, j, k, tidx, tidy, tidz, stridex, stridey, stridez

         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

         do k = tidz, ke_d, stridez
            do j = tidy, je_d + 1, stridey
               do i = tidx, ie_d, stridex
                  if (v0_d(i, j, k) > 0) then
                     d1 = var(i, j - 1, k) - var(i, j - 2, k)
                     d2 = var(i, j, k) - var(i, j - 1, k)
                     cf = var(i, j - 1, k)
                  else
                     d1 = var(i, j, k) - var(i, j + 1, k)
                     d2 = var(i, j - 1, k) - var(i, j, k)
                     cf = var(i, j, k)
                  end if
                  cf = cf + rlim_cuda(d1, d2)
                  duml_d(i, j, k) = cf*v0_d(i, j, k)*dyi_d
                  dumu_d(i, j - 1, k) = -cf*v0_d(i, j, k)*dyi_d
               end do
            end do
         end do
      end subroutine advecc_kappa_dvcdy_cuda

      ! -d(wc)/dz (stretched grid)
      attributes(global) subroutine advecc_kappa_dwcdz_cuda(hi, hj, hk, var)
         implicit none
         integer, value, intent(in) :: hi, hj, hk
         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d - hk:ke_d + hk), intent(in)    :: var
         real    :: cf, d1, d2
         integer :: i, j, k, tidx, tidy, tidz, stridex, stridey, stridez

         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

         do k = tidz + 1, ke_d + 1, stridez
            do j = tidy, je_d, stridey
               do i = tidx, ie_d, stridex
                  if (w0_d(i, j, k) > 0) then
                     d1 = (var(i, j, k - 1) - var(i, j, k - 2))*dzhci_d(k - 1)
                     d2 = (var(i, j, k) - var(i, j, k - 1))*dzhci_d(k)
                     cf = var(i, j, k - 1)
                  else
                     d1 = (var(i, j, k) - var(i, j, k + 1))*dzhci_d(k + 1)
                     d2 = (var(i, j, k - 1) - var(i, j, k))*dzhci_d(k)
                     cf = var(i, j, k)
                  end if
                  cf = cf + dzfc_d(k)*rlim_cuda(d1, d2)
                  duml_d(i, j, k) = cf*w0_d(i, j, k)*dzfci_d(k)
                  dumu_d(i, j, k - 1) = -cf*w0_d(i, j, k)*dzfci_d(k - 1)
               end do
            end do
         end do
      end subroutine advecc_kappa_dwcdz_cuda

      attributes(global) subroutine advecc_kappa_add_cuda(hi, hj, hk, varp)
         implicit none
         integer, value, intent(in) :: hi, hj, hk
         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d:ke_d + hk), intent(inout) :: varp
         integer :: i, j, k, tidx, tidy, tidz, stridex, stridey, stridez
         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
         do k = tidz, ke_d + hk, stridez
            do j = tidy - hj, je_d + hj, stridey
               do i = tidx - hi, ie_d + hi, stridex
                  varp(i,j,k) = varp(i,j,k) + dumu_d(i,j,k) + duml_d(i,j,k)
               end do
            end do
         end do
      end subroutine advecc_kappa_add_cuda

      !> Determination of the limiter function
      attributes(device) real function rlim_cuda(d1, d2)
         implicit none
         real, intent(in) :: d1 !< Scalar flux at 1.5 cells upwind
         real, intent(in) :: d2 !< Scalar flux at 0.5 cells upwind
         real :: ri, phir
         ri = (d2 + eps1_d)/(d1 + eps1_d)
         phir = max(0., min(2.*ri, min(1./3.+2./3.*ri, 2.)))
         rlim_cuda = 0.5*phir*d1
      end function rlim_cuda

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


      !> Advection at cell center through upwind scheme
      attributes(global) subroutine advecc_upw_cuda(hi, hj, hk, putin, putout)
         implicit none
         integer, value, intent(in) :: hi, hj, hk
         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d - hk:ke_d + hk), intent(in)    :: putin
         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d:ke_d + hk)     , intent(inout) :: putout
         integer :: i, j, k, tidx, tidy, tidz, stridex, stridey, stridez
         real    :: fluxr, fluxl, fluxb, fluxf, fluxu, fluxd
         
         call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)
         
         do k = tidz, ke_d, stridez
            do j = tidy, je_d, stridey
               do i = tidx, ie_d, stridex
                   if (u0_d(i+1, j, k) > 0) then   
                      fluxr = putin(i, j, k)
                   else
                      fluxr = putin(i + 1, j, k)
                   end if
 
                   if (u0_d(i, j, k) > 0) then
                      fluxl = putin(i - 1, j, k)
                   else
                      fluxl = putin(i, j, k)
                   end if
 
                   if (v0_d(i, j+1, k) > 0) then
                      fluxb = putin(i, j, k)
                   else
                      fluxb = putin(i, j + 1, k)
                   end if
 
                   if (v0_d(i, j, k) > 0) then
                      fluxf = putin(i, j - 1, k)
                   else
                      fluxf = putin(i, j, k)
                   end if
 
                   if (w0_d(i, j, k+1) > 0) then
                      fluxu = putin(i, j, k)
                   else
                      fluxu = putin(i, j, k + 1)
                   end if
 
                   if (w0_d(i, j, k) > 0) then
                      fluxd = putin(i, j, k - 1)
                   else
                      fluxd = putin(i, j, k)
                   end if
 
                   putout(i, j, k) = putout(i, j, k) &
                                     - (u0_d(i + 1, j, k)*fluxr - u0_d(i, j, k)*fluxl)*dxfci_d(i) & ! -d(uc)/dx (stretched grid)
                                     - (v0_d(i, j + 1, k)*fluxb - v0_d(i, j, k)*fluxf)*dyi_d &      ! -d(vc)/dy (no stretched grid)
                                     - (w0_d(i, j, k + 1)*fluxu - w0_d(i, j, k)*fluxd)*dzfci_d(k)   ! -d(wc)/dz (stretched grid)
               end do
            end do
         end do
      end subroutine advecc_upw_cuda


      attributes(global) subroutine shiftedPBCs_cuda
         implicit none
         integer :: tidx, tidy, tidz, stridex, stridey, stridez
         integer :: i, j, k, ig
         real :: vs, rk3coef

         if (ds_d > 0) then

            call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

            rk3coef = dt_d / (4. - dble(rk3step_d))

            do i = tidx, ie_d, stridex
               ig = i + zstart_d(1) - 1 ! global i position
               if (ig > int(itot_d/2)) then
                  do j = tidy, je_d, stridey
                     do k = tidz, ke_d, stridez
                        vs = 0.5 * pi_d * ds_d / (0.5*xlen_d) * u0av_d(k) * sin(pi_d*(xh_d(ig)-xh_d(int(itot_d/2))) / (0.5*xlen_d))
                        up_d(i,j,k) = up_d(i,j,k) - vs * (u0_d(i,j,k) - u0_d(i,j-1,k)) * dyi_d
                        vp_d(i,j,k) = vp_d(i,j,k) - vs * (v0_d(i,j,k) - v0_d(i,j-1,k)) * dyi_d
                        wp_d(i,j,k) = wp_d(i,j,k) - vs * (w0_d(i,j,k) - w0_d(i,j-1,k)) * dyi_d
                     end do
                  end do
               end if
            end do

         end if
      end subroutine shiftedPBCs_cuda

#endif
end module modcuda
