module modgpu
        
#if defined(_GPU)        
   use cudafor
#endif
   implicit none
   save
   public :: exchange_halo_z_gpu, exchange_halo_z_gpu_alloc_z, exchange_halo_z_gpu_ihc, &
             exchange_halo_z_gpu_scalar, exchange_halo_z_gpu_kb_opt_ih, exchange_halo_z_gpu_kbke_opt_ih
   
   contains
   
   subroutine exchange_halo_z_gpu(dummyvar)
     use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh
     use decomp_2d, only : exchange_halo_z
     implicit none

     real, dimension(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh), intent(inout) :: dummyvar

#if defined(_GPU)
     attributes(device) :: dummyvar
#endif

     call exchange_halo_z(dummyvar)
     
   end subroutine exchange_halo_z_gpu
   
   
   subroutine exchange_halo_z_gpu_alloc_z(dummyvar)
     use decomp_2d, only : exchange_halo_z
     implicit none
     
     real, intent(inout) :: dummyvar(:,:,:)

#if defined(_GPU)
     attributes(device) :: dummyvar
#endif

     call exchange_halo_z(dummyvar)
     
   end subroutine exchange_halo_z_gpu_alloc_z
   

   subroutine exchange_halo_z_gpu_ihc(dummyvar)
     use modglobal, only : ib, ie, ihc, jb, je, jhc, kb, ke, khc
     use decomp_2d, only : exchange_halo_z
     implicit none

     real, dimension(ib-ihc:ie+ihc, jb-jhc:je+jhc, kb-khc:ke+khc), intent(inout) :: dummyvar

#if defined(_GPU)
     attributes(device) :: dummyvar
#endif

     call exchange_halo_z(dummyvar, opt_zlevel=(/ihc,jhc,khc/))

   end subroutine exchange_halo_z_gpu_ihc


   subroutine exchange_halo_z_gpu_scalar(sv0, svm)

     use modglobal, only : ib, ie, ihc, jb, je, jhc, kb, ke, khc, nsv
     use decomp_2d, only : exchange_halo_z
     implicit none

     integer n

     real, dimension(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb-khc:ke+khc,nsv), intent(inout) :: sv0    !<  temporary copy of sv0
     real, dimension(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb-khc:ke+khc,nsv), intent(inout) :: svm    !<  temporary copy of svm

#if defined(_GPU)
     attributes(device) :: sv0, svm
#endif

     do n = 1, nsv
        call exchange_halo_z(sv0(:, :, :, n), opt_zlevel=(/ihc,jhc,khc/))
        call exchange_halo_z(svm(:, :, :, n), opt_zlevel=(/ihc,jhc,khc/))
     enddo

   end subroutine exchange_halo_z_gpu_scalar

   
   subroutine exchange_halo_z_gpu_kb_opt_ih(dummyvar)
     use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh
     use decomp_2d, only : exchange_halo_z
     implicit none

     real, dimension(ib-ih:ie+ih, jb-jh:je+jh, kb:ke+kh), intent(inout) :: dummyvar

#if defined(_GPU)
     attributes(device) :: dummyvar
#endif

     call exchange_halo_z(dummyvar, opt_zlevel=(/ih,jh,0/))
     
   end subroutine exchange_halo_z_gpu_kb_opt_ih
   
   
   subroutine exchange_halo_z_gpu_kbke_opt_ih(dummyvar)
     use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke
     use decomp_2d, only : exchange_halo_z
     implicit none

     real, dimension(ib-1:ie+1,jb-1:je+1,kb:ke), intent(inout) :: dummyvar

#if defined(_GPU)
     attributes(device) :: dummyvar
#endif

     call exchange_halo_z(dummyvar, opt_zlevel=(/ih,jh,0/))
      
   end subroutine exchange_halo_z_gpu_kbke_opt_ih
   
   
end module modgpu
