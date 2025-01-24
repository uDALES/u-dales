module cudamodule
#if defined(_GPU)        
   use cudafor
   use modglobal, only : ib, ie, jb, je, kb, ke, ih, jh, kh, ihc, jhc, khc, dxi5, dyi5, dzf, dzhi, dzfi5, dzhiq, dxiq, dyiq, dxi, dyi, &
                         loneeqn, ltempeq, lmoist, nsv
   use modfields, only : u0, v0, w0, pres0, e120, thl0, qt0, sv0, up, vp, wp, e12p, thlp, qtp, svp
   implicit none
 !  public :: initCUDA, updateDeviceVariables, checkCUDA, saxpy, advecc_2nd_cuda, advecu_2nd_cuda, advecv_2nd_cuda
   save

   type(dim3) :: griddim, blockdim

   integer, device :: ib_d, ie_d, jb_d, je_d, kb_d, ke_d, ih_d, jh_d, kh_d
   real, device    :: dxi5_d, dyi5_d, dxiq_d, dyiq_d, dxi_d, dyi_d
   real, device, allocatable :: dzf_d(:), dzhi_d(:), dzfi5_d(:), dzhiq_d(:)
   real, device, allocatable :: u0_d(:,:,:), v0_d(:,:,:), w0_d(:,:,:), pres0_d(:,:,:), e120_d(:,:,:), thl0_d(:,:,:), qt0_d(:,:,:), sv0_d(:,:,:,:)
   real, device, allocatable :: up_d(:,:,:), vp_d(:,:,:), wp_d(:,:,:), e12p_d(:,:,:), thlp_d(:,:,:), qtp_d(:,:,:), svp_d(:,:,:,:)

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
         end if

         if (ltempeq) then
            allocate(thl0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
            allocate(thlp_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
         end if

         if (lmoist) then
            allocate(qt0_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
            allocate(qtp_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
         end if

         if (nsv>0) then
            allocate(sv0_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb-khc:ke+khc,nsv))
            allocate(svp_d(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc,nsv))
         end if

      end subroutine initCUDA

      subroutine exitCUDA
         implicit none
         deallocate(dzf_d, dzhi_d, dzfi5_d, dzhiq_d)
         deallocate(u0_d, v0_d, w0_d, pres0_d)
         deallocate(up_d, vp_d, wp_d)
         if (loneeqn) deallocate(e120_d, e12p_d)
         if (ltempeq) deallocate(thl0_d, thlp_d)
         if (lmoist) deallocate(qt0_d, qtp_d)
         if (nsv>0) deallocate(sv0_d, svp_d)
      end subroutine exitCUDA

      subroutine updateDevice
         implicit none
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
         end if
         if (lmoist) then
            qt0_d = qt0
            qtp_d = qtp
         end if
     !    if (nsv>0) then
     !       sv0_d = sv0
     !       svp_d = svp
     !    end if
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
         end if
         if (lmoist) then
     !       qt0 = qt0_d
            qtp = qtp_d
         end if
     !    if (nsv>0) then
     !       sv0 = sv0_d
     !       svp = svp_d
     !    end if
      end subroutine updateHost

      subroutine checkCUDA(istat, kernelname)
         implicit none
         integer, intent(in)          :: istat
         character(len=*), intent(in) :: kernelname
         if(istat /= cudaSuccess) then
            write(*,*) "Error in ", trim(kernelname), ": ", cudaGetErrorString(istat)
         !!   STOP 1
         end if
      end subroutine checkCUDA

      attributes(global) subroutine saxpy(c)
         implicit none
         real :: c(ib_d:ie_d,jb_d:je_d,kb_d:ke_d)
         integer :: i, j, k, tidx, tidy, stridex, stridey, kk
         tidx = (blockIDx%x - 1) * blockDim%x + threadIdx%x
         tidy = (blockIDx%y - 1) * blockDim%y + threadIdx%y
         stridex = gridDim%x * blockDim%x
         stridey = gridDim%y * blockDim%y
         
         do k = kb_d, ke_d
            do j = tidy, je_d, stridey
               do i = tidx, ie_d, stridex
                  c(i,j,k) = ih_d + jh_d + kh_d
               end do
            end do
         end do
      end subroutine saxpy
      
      !> Advection at cell center
      attributes(global) subroutine advecc_2nd_cuda(hi, hj, hk, putin, putout)
         implicit none
         
         integer, value, intent(in) :: hi, hj, hk

         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d - hk:ke_d + hk), intent(in)    :: putin
         real, dimension(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d:ke_d + hk)     , intent(inout) :: putout
         
         integer :: tidx, tidy, tidz, stridex, stridey, stridez
         integer :: i, j, k, im, ip, jm, jp, km, kp
               
         tidx = (blockIDx%x - 1) * blockDim%x + threadIdx%x
         tidy = (blockIDx%y - 1) * blockDim%y + threadIdx%y
         tidz = (blockIDx%z - 1) * blockDim%z + threadIdx%z
         stridex = gridDim%x * blockDim%x
         stridey = gridDim%y * blockDim%y
         stridez = gridDim%z * blockDim%z

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
     
      !> Advection at the u point.
      attributes(global) subroutine advecu_2nd_cuda(putin, putout)
         implicit none
             
         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d - kh_d:ke_d + kh_d), intent(in)    :: putin
         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d:ke_d + kh_d)       , intent(inout) :: putout

         integer :: tidx, tidy, tidz, stridex, stridey, stridez         
         integer :: i, j, k, im, ip, jm, jp, km, kp
                                                               
         tidx = (blockIDx%x - 1) * blockDim%x + threadIdx%x
         tidy = (blockIDx%y - 1) * blockDim%y + threadIdx%y                                               
         tidz = (blockIDx%z - 1) * blockDim%z + threadIdx%z                                             
         stridex = gridDim%x * blockDim%x                              
         stridey = gridDim%y * blockDim%y                              
         stridez = gridDim%z * blockDim%z

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

      !> Advection at the v point.
      attributes(global) subroutine advecv_2nd_cuda(putin, putout)
         implicit none
         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d - kh_d:ke_d + kh_d), intent(in)    :: putin
         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d:ke_d + kh_d)       , intent(inout) :: putout

         integer :: tidx, tidy, tidz, stridex, stridey, stridez
         integer :: i, j, k, im, ip, jm, jp, km, kp

         tidx = (blockIDx%x - 1) * blockDim%x + threadIdx%x
         tidy = (blockIDx%y - 1) * blockDim%y + threadIdx%y
         tidz = (blockIDx%z - 1) * blockDim%z + threadIdx%z
         stridex = gridDim%x * blockDim%x
         stridey = gridDim%y * blockDim%y
         stridez = gridDim%z * blockDim%z

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

      !> Advection at the w point.
      attributes(global) subroutine advecw_2nd_cuda(putin, putout)
         implicit none
         
         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d - kh_d:ke_d + kh_d), intent(in)    :: putin
         real, dimension(ib_d - ih_d:ie_d + ih_d, jb_d - jh_d:je_d + jh_d, kb_d:ke_d + kh_d)       , intent(inout) :: putout

         integer :: tidx, tidy, tidz, stridex, stridey, stridez
         integer :: i, j, k, im, ip, jm, jp, km, kp

         tidx = (blockIDx%x - 1) * blockDim%x + threadIdx%x
         tidy = (blockIDx%y - 1) * blockDim%y + threadIdx%y
         tidz = (blockIDx%z - 1) * blockDim%z + threadIdx%z
         stridex = gridDim%x * blockDim%x
         stridey = gridDim%y * blockDim%y
         stridez = gridDim%z * blockDim%z

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

#endif
end module cudamodule
