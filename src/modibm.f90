!!> \file modibm.f90
!!!  adds forcing terms for immersed boundaries
!
!>
!!  \author Jasper Thomas TU Delft / Ivo Suter Imperial College London
!
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!
module modibm
   use mpi
#if defined(_GPU)
   use modcuda, only : facT1_d, fachf_d, facef_d, fac_tau_d, &
                       fac_htc_d, fac_cth_d, fac_pres_d, fac_pres2_d, fac_stage, &
                       facqsat_d, fachurel_d, facf_d, updateFacetPropsDevice
#endif
   use modibmdata
   !use wf_uno
   implicit none
   save
   public :: initibm, ibmnorm, ibmwallfun, bottom, createmasks, cell_index, &
             nsolpts_u, nsolpts_v, nsolpts_w, nsolpts_c, &
             nbndpts_u, nbndpts_v, nbndpts_w, nbndpts_c, &
             nfctsecs_u, nfctsecs_v, nfctsecs_w, nfctsecs_c, &
             mask_u, mask_v, mask_w, mask_c
   ! The host wall functions and the geometry they read, exposed so that both
   ! tests.f90 (which checks them against independent expectations on the CPU)
   ! and tests_cuda.f90 (which checks the device port against them) can drive
   ! them directly.
#if !defined(_GPU) || defined(UDALES_DEBUG)
   public :: wallfunmom, wallfunheat, local_coords, check_wallfun_cache, &
             diffu_corr, diffv_corr, diffw_corr, diffc_corr, &
             fac_tau_raw, fac_htc_raw, fac_cth_raw, fac_pres_raw, fac_pres2_raw
#endif
   public :: bound_info_type, bound_info_u, bound_info_v, bound_info_w, bound_info_c
#if defined(_GPU)
   ! Re-exported from modcuda so that the self-tests can reach them by their
   ! usual name.
   public :: fac_tau_d, fac_htc_d, fac_cth_d, fac_pres_d, fac_pres2_d, fachf_d, facef_d
   public :: mask_u_d, mask_v_d, mask_w_d, mask_c_d, &
             bndpts_u_d, bndpts_v_d, bndpts_w_d, bndpts_c_d, &
             init_ibm_device, &
             diffu_corr_device, diffv_corr_device, diffw_corr_device, diffc_corr_device, &
             wallfunmom_dir_device, &
             wallfunheat_dir_device, faclGR_d, check_ibm_section_cache
#endif

    abstract interface
      function interp_velocity(i,j,k)
        real :: interp_velocity(3)
        integer, intent(in) :: i, j, k
      end function interp_velocity
    end interface

    abstract interface
      real function interp_temperature(i,j,k)
        integer, intent(in) :: i, j, k
      end function interp_temperature
    end interface

   logical :: lnorec = .false.
   
   ! read from namoptions
   integer :: nsolpts_u, nsolpts_v, nsolpts_w, nsolpts_c, &
              nbndpts_u, nbndpts_v, nbndpts_w, nbndpts_c, &
              nfctsecs_u, nfctsecs_v, nfctsecs_w, nfctsecs_c

   real, allocatable, target, dimension(:,:,:) :: mask_u, mask_v, mask_w, mask_c

   ! Column layout of the two per-section caches build_wallfun_cache fills.
   ! Both the host wall functions and, through mirror_sections, the device
   ! kernels index them with these.
   integer, parameter :: RECW_U = 1, RECW_V = 4, RECW_W = 7, RECW_C = 10, RECW_N = 12
   integer, parameter :: WFC_LOGDZ = 1, WFC_LOGZH = 2, WFC_COEF = 3, &
                         WFC_CM = 4, WFC_CH = 5, WFC_CTMNEU = 6, WFC_N = 6

#if defined(_GPU)
   ! Device mirrors of the IBM geometry. All of it is built once in initibm and
   ! never changes afterwards, so a single transfer at the end of initibm is
   ! enough - there is no per-step counterpart to this.
   real,    device, allocatable :: mask_u_d(:,:,:), mask_v_d(:,:,:), mask_w_d(:,:,:), mask_c_d(:,:,:)
   ! The boundary-point lists live in bound_info_type, whose allocatable
   ! components a kernel cannot reach, so they are mirrored as flat arrays in
   ! the same way heatpump mirrors idhppts_local.
   integer, device, allocatable :: bndpts_u_d(:,:), bndpts_v_d(:,:), bndpts_w_d(:,:), bndpts_c_d(:,:)

   ! Facet-section geometry, one set per staggered direction. The prefix is the
   ! section set; recids_<x> is the cell on the <x> grid that the reconstruction
   ! point falls in, so u_recids_v_d is "for the u sections, the v-grid cell".
   !
   ! None of this changes once the facet files are read, so the arrays hold the
   ! precomputed form the kernels use rather than the form the host keeps. See
   ! mirror_sections, which builds them, for what each one holds. recpts is not
   ! mirrored at all: it is only needed to build _secdist_d and _recw_d.
   !
   real,    device, allocatable :: u_secareas_d(:), u_secdist_d(:)
   real,    device, allocatable :: u_recw_d(:,:), u_wfc_d(:,:)
   integer, device, allocatable :: u_secfacids_d(:), u_secbndloc_d(:,:)
   integer, device, allocatable :: u_recids_u_d(:,:), u_recids_v_d(:,:), u_recids_w_d(:,:), u_recids_c_d(:,:)
   logical, device, allocatable :: u_lskipsec_d(:), u_lcomprec_d(:)

   real,    device, allocatable :: v_secareas_d(:), v_secdist_d(:)
   real,    device, allocatable :: v_recw_d(:,:), v_wfc_d(:,:)
   integer, device, allocatable :: v_secfacids_d(:), v_secbndloc_d(:,:)
   integer, device, allocatable :: v_recids_u_d(:,:), v_recids_v_d(:,:), v_recids_w_d(:,:), v_recids_c_d(:,:)
   logical, device, allocatable :: v_lskipsec_d(:), v_lcomprec_d(:)

   real,    device, allocatable :: w_secareas_d(:), w_secdist_d(:)
   real,    device, allocatable :: w_recw_d(:,:), w_wfc_d(:,:)
   integer, device, allocatable :: w_secfacids_d(:), w_secbndloc_d(:,:)
   integer, device, allocatable :: w_recids_u_d(:,:), w_recids_v_d(:,:), w_recids_w_d(:,:), w_recids_c_d(:,:)
   logical, device, allocatable :: w_lskipsec_d(:), w_lcomprec_d(:)

   real,    device, allocatable :: c_secareas_d(:), c_secdist_d(:)
   real,    device, allocatable :: c_recw_d(:,:), c_wfc_d(:,:)
   integer, device, allocatable :: c_secfacids_d(:), c_secbndloc_d(:,:)
   integer, device, allocatable :: c_recids_u_d(:,:), c_recids_v_d(:,:), c_recids_w_d(:,:), c_recids_c_d(:,:)
   logical, device, allocatable :: c_lskipsec_d(:), c_lcomprec_d(:)

   ! Facet properties. facz0 and facz0h are indexed from 0 on the host. facT1_d
   ! lives in modcuda instead, because it is refreshed inside the time loop.
   real, device, allocatable :: facnorm_d(:,:), facz0_d(:), facz0h_d(:)
   ! facqsat_d, fachurel_d and facf_d live in modcuda: the energy balance
   ! rewrites them during the run, so they are refreshed there with facT.
   ! faclGR is fixed geometry and stays here.
   logical, device, allocatable :: faclGR_d(:)

   ! The per-facet accumulators live in modcuda: fachf_d and facef_d cross the
   ! bus inside the time loop, and the rest are kept beside them so that all the
   ! facet transfers are declared in one place.
   ! Per-section moisture tendency increment, scattered into qtp in a second
   ! pass. Keeping it out of the main kernel means that kernel never takes qtp
   ! as an argument, which matters because qtp_d does not exist when lmoist is
   ! off - and passing an unallocated device array faults the kernel even where
   ! no branch reads it.
   real, device, allocatable :: qflux_sec_d(:), hflux_sec_d(:)

   ! Scratch copy of a tendency taken before a wall function runs, so that its
   ! contribution alone can be added to tau_x/tau_y/tau_z or thl_flux.
   real, device, allocatable :: rhs_ibm_d(:,:,:)
#endif

   TYPE solid_info_type
     integer :: nsolpts
     integer, allocatable :: solpts(:,:)
     logical, allocatable :: lsolptsrank(:) !
     integer, allocatable :: solptsrank(:) ! indices of points on current rank
     integer :: nsolptsrank
     integer, allocatable :: solpts_loc(:,:)
   end TYPE solid_info_type

   type(solid_info_type) :: solid_info_u, solid_info_v, solid_info_w, solid_info_c

   TYPE bound_info_type
     integer :: nbndpts
     integer, allocatable :: bndpts(:,:) ! ijk location of fluid boundary point
     !real, allocatable    :: intpts(:,:) ! xyz location of boundary intercept point
     !real, allocatable    :: bndvec(:,:) ! vector from boundary to fluid point (normalised)
     real, allocatable    :: recpts(:,:) ! xyz location of reconstruction point
     integer, allocatable :: recids_u(:,:) ! ijk location of u grid cell that rec point is in
     integer, allocatable :: recids_v(:,:) ! ijk location of u grid cell that rec point is in
     integer, allocatable :: recids_w(:,:) ! ijk location of u grid cell that rec point is in
     integer, allocatable :: recids_c(:,:) ! ijk location of u grid cell that rec point is in
     real, allocatable    :: bnddst(:) ! distance between surface & bound point
     integer, allocatable :: bndptsrank(:) ! indices of points on current rank
     !integer, allocatable :: bndpts_loc(:,:) ! indices of points on current rank
     logical, allocatable :: lcomprec(:) ! Switch whether reconstruction point is a computational point
     logical, allocatable :: lskipsec(:) ! Switch whether to skip finding the shear stress at this point
     integer :: nbndptsrank
     integer, allocatable :: bndpts_loc(:,:) ! ijk location of fluid boundary point on rank

     integer :: nfctsecs
     integer, allocatable :: secbndptids(:)
     integer, allocatable :: secfacids(:)
     real,    allocatable :: secareas(:)
     integer, allocatable :: fctsecsrank(:)
     integer :: nfctsecsrank
     integer, allocatable :: secfacids_loc(:)
     real   , allocatable :: secareas_loc(:)
     integer, allocatable :: secbndpts_loc(:,:)
     real   , allocatable :: bnddst_loc(:)
     real   , allocatable :: recpts_loc(:,:)
     integer, allocatable :: recids_u_loc(:,:)
     integer, allocatable :: recids_v_loc(:,:)
     integer, allocatable :: recids_w_loc(:,:)
     integer, allocatable :: recids_c_loc(:,:)
     logical, allocatable :: lcomprec_loc(:)
     logical, allocatable :: lskipsec_loc(:)

     ! Per-section quantities that follow from the geometry alone, built once by
     ! build_wallfun_cache. The host wall functions read these directly and the
     ! device mirrors are a straight copy of them, so there is one derivation
     ! rather than two that have to be kept in step.
     real   , allocatable :: secdist_loc(:)   ! wall distance, reconstruction offset included
     integer, allocatable :: bndloc_loc(:,:)  ! boundary point, rank-local ijk
     integer, allocatable :: recloc_u(:,:)    ! reconstruction cells, rank-local ijk
     integer, allocatable :: recloc_v(:,:)
     integer, allocatable :: recloc_w(:,:)
     integer, allocatable :: recloc_c(:,:)
     real   , allocatable :: recw_loc(:,:)    ! trilinear offsets, RECW_N columns
     real   , allocatable :: wfc_loc(:,:)     ! roughness terms, WFC_N columns
   end TYPE bound_info_type

   type(bound_info_type) :: bound_info_u, bound_info_v, bound_info_w, bound_info_c

   integer :: nstatfac=7, ncidfac, nrecfac=0
   character(80), allocatable :: ncstatfac(:,:)
   character(80) :: facname = 'fac.xxx.nc'
   character(80),dimension(1,4) :: tncstatfac
   real, allocatable :: varsfac(:,:)
#if !defined(_GPU) || defined(UDALES_DEBUG)
   ! Per-facet totals as the host wall functions accumulate them, before the
   ! area normalisation and MPI reduction that only happen under lwritefac.
   ! They exist so the tests have something to compare against, so they follow
   ! the host routines out of a GPU release build.
   real, allocatable :: fac_tau_raw(:)
   real, allocatable :: fac_htc_raw(:), fac_cth_raw(:), fac_pres_raw(:), fac_pres2_raw(:)
#endif
   real, allocatable :: fac_tau_x(:)
   real, allocatable :: fac_tau_y(:)
   real, allocatable :: fac_tau_z(:)
   real, allocatable :: fac_pres(:)
   real, allocatable :: fac_pres2(:)
   real, allocatable :: fac_htc(:)
   real, allocatable :: fac_cth(:)
   real, allocatable :: fac_tau_x_av(:)
   real, allocatable :: fac_tau_y_av(:)
   real, allocatable :: fac_tau_z_av(:)
   real, allocatable :: fac_pres_av(:)
   real, allocatable :: fac_pres2_av(:)
   real, allocatable :: fac_htc_av(:)
   real, allocatable :: fac_cth_av(:)

   contains

   !> Index of the cell containing coordinate p on a monotonically increasing grid.
   !!
   !! Returns the position (1-based, as the grid arrays all start at index 1) of the
   !! last element of `grid` that p is greater than or equal to, and 0 when p lies
   !! below grid(1). This is the same result as
   !!     findloc(p >= grid, .true., 1, back=.true.)
   !! which is how it used to be written. Because `grid` is monotonically increasing,
   !! the mask (p >= grid) is .true. for a leading run of elements only, so the last
   !! .true. position is simply the number of .true. values and count() is exact.
   !!
   !! count() is used rather than findloc() because the NVHPC runtime does not
   !! implement FINDLOC for logical arrays - a GPU build aborts at run time with
   !! "FINDLOC: unimplemented for data type" the first time this line is reached.
   !! That only happens for geometries with facet sections close enough to a wall to
   !! need full reconstruction, which is why it went unnoticed on grid-aligned cases.
   !!
   !! The monotonicity requirement is asserted by tests_ibm_cell_lookup.
   pure integer function cell_index(p, grid)
     real, intent(in) :: p
     real, intent(in) :: grid(:)

     cell_index = count(p >= grid)

   end function cell_index

   subroutine initibm
     use modglobal, only : libm, xhat, yhat, zhat, vec0, &
                           ib, ie, ih, jb, je, jh, kb, ke, kh, nsv, &
                           iwallmom, lmoist, ltempeq, cexpnr, nfcts, lwritefac
     use m_halo, only : halo_exchange
     use modmpi,    only : myid
     use modstat_nc,only: open_nc, define_nc, ncinfo, writestat_dims_nc

     real, allocatable :: rhs(:,:,:)

     if (.not. libm) return

     solid_info_u%nsolpts = nsolpts_u
     solid_info_v%nsolpts = nsolpts_v
     solid_info_w%nsolpts = nsolpts_w
     call initibmnorm('solid_u.txt', solid_info_u)
     call initibmnorm('solid_v.txt', solid_info_v)
     call initibmnorm('solid_w.txt', solid_info_w)

     ! Define (real) masks
     ! Hopefully this can be removed eventually if (integer) IIx halos can be communicated
     ! These are only used in modibm, to cancel subgrid term across solid boundaries
     allocate(mask_u(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)); mask_u = 1.
     allocate(mask_v(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)); mask_v = 1.
     allocate(mask_w(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)); mask_w = 1.
     mask_w(:,:,kb) = 0.     ! In future this shouldn't be needed?
     mask_u(:,:,kb-kh) = 0.
     mask_v(:,:,kb-kh) = 0.
     mask_w(:,:,kb-kh) = 0.

     allocate(rhs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
     call solid(solid_info_u, mask_u, rhs, 0., ih, jh, kh)
     call solid(solid_info_v, mask_v, rhs, 0., ih, jh, kh)
     call solid(solid_info_w, mask_w, rhs, 0., ih, jh, kh)
!$acc data create(mask_u, mask_v, mask_w)
!$acc update device(mask_u, mask_v, mask_w)
     call halo_exchange(mask_u, 3)!, opt_zlevel=(/ih,jh,0/))
     call halo_exchange(mask_v, 3)!, opt_zlevel=(/ih,jh,0/))
     call halo_exchange(mask_w, 3)!, opt_zlevel=(/ih,jh,0/))
!$acc update host(mask_u, mask_v, mask_w)
!$acc end data

#if !defined(_GPU) || defined(UDALES_DEBUG)
     if (iwallmom > 1) allocate(fac_tau_raw(1:nfcts))
     if (ltempeq .or. lmoist .or. lwritefac) then
       allocate(fac_htc_raw(1:nfcts), fac_cth_raw(1:nfcts))
       allocate(fac_pres_raw(1:nfcts), fac_pres2_raw(1:nfcts))
     end if
#endif

     if (iwallmom > 1) then
       bound_info_u%nbndpts = nbndpts_u
       bound_info_v%nbndpts = nbndpts_v
       bound_info_w%nbndpts = nbndpts_w
       bound_info_u%nfctsecs = nfctsecs_u
       bound_info_v%nfctsecs = nfctsecs_v
       bound_info_w%nfctsecs = nfctsecs_w
       call initibmwallfun('fluid_boundary_u.txt', 'facet_sections_u.txt', xhat, bound_info_u)
       call initibmwallfun('fluid_boundary_v.txt', 'facet_sections_v.txt', yhat, bound_info_v)
       call initibmwallfun('fluid_boundary_w.txt', 'facet_sections_w.txt', zhat, bound_info_w)
     end if

     if (ltempeq .or. lmoist .or. nsv>0 .or. lwritefac) then
       solid_info_c%nsolpts = nsolpts_c
       call initibmnorm('solid_c.txt', solid_info_c)

       bound_info_c%nbndpts = nbndpts_c
       bound_info_c%nfctsecs = nfctsecs_c
       call initibmwallfun('fluid_boundary_c.txt', 'facet_sections_c.txt', vec0, bound_info_c)

       allocate(mask_c(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)); mask_c = 1.
       mask_c(:,:,kb-kh) = 0.
       call solid(solid_info_c, mask_c, rhs, 0., ih, jh, kh)
!$acc data create(mask_c)
!$acc update device(mask_c)
       call halo_exchange(mask_c, 3)!, opt_zlevel=(/ih,jh,0/))
!$acc update host(mask_c)
!$acc end data
     end if

     deallocate(rhs)

     ! write facet stresses and pressure to fac.xxx.nc
     if (lwritefac) then
       allocate(fac_tau_x(1:nfcts))
       allocate(fac_tau_y(1:nfcts))
       allocate(fac_tau_z(1:nfcts))
       allocate(fac_pres(1:nfcts))
       allocate(fac_pres2(1:nfcts))
       allocate(fac_htc(1:nfcts))
       allocate(fac_cth(1:nfcts))
       fac_tau_x = 0.
       fac_tau_y = 0.
       fac_tau_z = 0.
       fac_pres = 0.
       fac_pres2=0.
       fac_htc = 0.
       fac_cth = 0.
       allocate(fac_tau_x_av(1:nfcts))
       allocate(fac_tau_y_av(1:nfcts))
       allocate(fac_tau_z_av(1:nfcts))
       allocate(fac_pres_av(1:nfcts))
       allocate(fac_pres2_av(1:nfcts))
       allocate(fac_htc_av(1:nfcts))
       allocate(fac_cth_av(1:nfcts))
       fac_tau_x_av = 0.
       fac_tau_y_av = 0.
       fac_tau_z_av = 0.
       fac_pres_av = 0.
       fac_pres2_av=0.
       fac_htc_av = 0.
       fac_cth_av = 0.

       facname(5:7) = cexpnr
       allocate(ncstatfac(nstatfac,4))
       call ncinfo(tncstatfac(1,:),'t', 'Time', 's', 'time')
       call ncinfo(ncstatfac( 1,:),'tau_x', 'tau_x', 'm^2/s^2','ft')
       call ncinfo(ncstatfac( 2,:),'tau_y', 'tau_y', 'm^2/s^2','ft')
       call ncinfo(ncstatfac( 3,:),'tau_z', 'tau_z', 'm^2/s^2','ft')
       call ncinfo(ncstatfac( 4,:),'pres', 'pressure', 'm^2/s^2','ft')
       call ncinfo(ncstatfac( 5,:),'htc', 'heat transfer coefficient', '','ft')
       call ncinfo(ncstatfac( 6,:),'cth', 'heat transfer coefficient (Ivo)', '','ft')
       call ncinfo(ncstatfac( 7,:),'pres_flc', 'pressure fluctuation', '','ft')

       if (myid==0) then
         call open_nc(facname, ncidfac, nrecfac, nfcts=nfcts)
         if (nrecfac==0) then
           call define_nc(ncidfac, 1, tncstatfac)
           call writestat_dims_nc(ncidfac)
         end if
         call define_nc(ncidfac, nstatfac, ncstatfac)
       end if
     end if

#if defined(_GPU)
     ! Last, so that the mirrored masks include the halos filled by the
     ! exchanges above.
     call init_ibm_device
#endif

   end subroutine initibm

#if defined(_GPU)
   !> Mirror the IBM geometry onto the device.
   !!
   !! Which of these arrays exist depends on iwallmom, ltempeq, lmoist, nsv and
   !! lwritefac. Rather than restating those conditions - and drifting from them
   !! later - this mirrors whatever initibm actually allocated.
   subroutine init_ibm_device
     use modglobal, only : libm, ib, ie, ih, jb, je, jh, kb, ke, kh, iwallmom, nfcts
     use initfac,   only : facnorm, facz0, facz0h, facT, faclGR, facqsat, fachurel, facf
     implicit none

     if (.not. libm) return

     ! Bounds are stated explicitly, as everywhere else device memory is
     ! allocated, so that the halo-relative extents the kernels index from
     ! (ib-ih, kb-kh) are visible here rather than inherited from the host
     ! array. test_ibm_device_geometry asserts they match.
     if (allocated(mask_u)) then
       allocate(mask_u_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
       mask_u_d = mask_u
     end if
     if (allocated(mask_v)) then
       allocate(mask_v_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
       mask_v_d = mask_v
     end if
     if (allocated(mask_w)) then
       allocate(mask_w_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
       mask_w_d = mask_w
     end if
     if (allocated(mask_c)) then
       allocate(mask_c_d(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
       mask_c_d = mask_c
     end if

     ! The point lists are 1-based on the host, so these bounds match directly.
     if (allocated(bound_info_u%bndpts_loc) .and. bound_info_u%nbndptsrank > 0) then
       allocate(bndpts_u_d(bound_info_u%nbndptsrank,3))
       bndpts_u_d = bound_info_u%bndpts_loc
     end if
     if (allocated(bound_info_v%bndpts_loc) .and. bound_info_v%nbndptsrank > 0) then
       allocate(bndpts_v_d(bound_info_v%nbndptsrank,3))
       bndpts_v_d = bound_info_v%bndpts_loc
     end if
     if (allocated(bound_info_w%bndpts_loc) .and. bound_info_w%nbndptsrank > 0) then
       allocate(bndpts_w_d(bound_info_w%nbndptsrank,3))
       bndpts_w_d = bound_info_w%bndpts_loc
     end if
     if (allocated(bound_info_c%bndpts_loc) .and. bound_info_c%nbndptsrank > 0) then
       allocate(bndpts_c_d(bound_info_c%nbndptsrank,3))
       bndpts_c_d = bound_info_c%bndpts_loc
     end if

     if (iwallmom > 1) then
     call mirror_sections(bound_info_u, u_secareas_d, u_secdist_d, u_recw_d, u_wfc_d, &
                          u_secfacids_d, u_secbndloc_d, &
                          u_recids_u_d, u_recids_v_d, u_recids_w_d, u_recids_c_d, &
                          u_lskipsec_d, u_lcomprec_d)
     call mirror_sections(bound_info_v, v_secareas_d, v_secdist_d, v_recw_d, v_wfc_d, &
                          v_secfacids_d, v_secbndloc_d, &
                          v_recids_u_d, v_recids_v_d, v_recids_w_d, v_recids_c_d, &
                          v_lskipsec_d, v_lcomprec_d)
     call mirror_sections(bound_info_w, w_secareas_d, w_secdist_d, w_recw_d, w_wfc_d, &
                          w_secfacids_d, w_secbndloc_d, &
                          w_recids_u_d, w_recids_v_d, w_recids_w_d, w_recids_c_d, &
                          w_lskipsec_d, w_lcomprec_d)

       allocate(facnorm_d(nfcts,3));  facnorm_d = facnorm
       allocate(facz0_d(0:nfcts));    facz0_d   = facz0
       allocate(facz0h_d(0:nfcts));   facz0h_d  = facz0h
       if (allocated(facT)) allocate(facT1_d(0:nfcts))
       allocate(fac_tau_d(nfcts))

       allocate(rhs_ibm_d(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))

     end if

     ! The heat wall function works on the cell-centred sections, which exist
     ! under a different condition from the momentum ones.
     if (allocated(bound_info_c%secareas_loc) .and. bound_info_c%nfctsecsrank > 0) then
       call mirror_sections(bound_info_c, c_secareas_d, c_secdist_d, c_recw_d, c_wfc_d, &
                            c_secfacids_d, c_secbndloc_d, &
                            c_recids_u_d, c_recids_v_d, c_recids_w_d, c_recids_c_d, &
                            c_lskipsec_d, c_lcomprec_d)

       if (.not. allocated(facnorm_d)) then
         allocate(facnorm_d(nfcts,3)); facnorm_d = facnorm
         allocate(facz0_d(0:nfcts));   facz0_d   = facz0
         allocate(facz0h_d(0:nfcts));  facz0h_d  = facz0h
       end if
       if (.not. allocated(facT1_d) .and. allocated(facT)) allocate(facT1_d(0:nfcts))

       allocate(faclGR_d(0:nfcts)); faclGR_d = faclGR
       if (allocated(facqsat)) then
         allocate(facqsat_d(0:nfcts))
         allocate(fachurel_d(0:nfcts))
         allocate(facf_d(0:nfcts,5))
       end if

       ! Filled straight away, not left to the first updateDevice. The solver
       ! would be fine either way, but the CUDA self-tests run before the time
       ! loop and would otherwise read an uninitialised surface temperature -
       ! zero, which divides through in the stability functions and hands the
       ! comparison NaN on the device and finite values on the host.
       call updateFacetPropsDevice

       allocate(fac_htc_d(nfcts), fac_cth_d(nfcts))
       allocate(fac_pres_d(nfcts), fac_pres2_d(nfcts))
       allocate(fachf_d(0:nfcts), facef_d(nfcts))
       fachf_d = 0.
       facef_d = 0.
       if (bound_info_c%nfctsecsrank > 0) then
         allocate(qflux_sec_d(bound_info_c%nfctsecsrank))
         allocate(hflux_sec_d(bound_info_c%nfctsecsrank))
       end if
     end if

   end subroutine init_ibm_device

   !> Copy one bound_info's per-rank section arrays onto the device.
   !!
   !! Everything the kernels need in a derived form - the wall distance, the
   !! rank-local indices, the trilinear offsets and the roughness terms - is
   !! built once by build_wallfun_cache, which the host wall functions read as
   !! well. Nothing is recomputed here, so the two paths cannot drift apart.
   subroutine mirror_sections(bi, areas, dst, recw, wfc, facids, bndloc, &
                              rids_u, rids_v, rids_w, rids_c, skip, comp)
     implicit none
     type(bound_info_type), intent(in) :: bi
     real,    device, allocatable, intent(out) :: areas(:), dst(:), recw(:,:), wfc(:,:)
     integer, device, allocatable, intent(out) :: facids(:), bndloc(:,:)
     integer, device, allocatable, intent(out) :: rids_u(:,:), rids_v(:,:), rids_w(:,:), rids_c(:,:)
     logical, device, allocatable, intent(out) :: skip(:), comp(:)
     integer :: n

     n = bi%nfctsecsrank

     ! A rank can own no sections at all. Allocating zero bytes on the device is
     ! a fatal error under nvfortran, so leave these unallocated and let the
     ! callers skip the direction entirely.
     if (n < 1) return

     allocate(areas(n));       areas  = bi%secareas_loc
     allocate(dst(n));         dst    = bi%secdist_loc
     allocate(recw(n,RECW_N)); recw   = bi%recw_loc
     allocate(wfc(n,WFC_N));   wfc    = bi%wfc_loc
     allocate(facids(n));      facids = bi%secfacids_loc
     allocate(bndloc(n,3));    bndloc = bi%bndloc_loc
     allocate(rids_u(n,3));    rids_u = bi%recloc_u
     allocate(rids_v(n,3));    rids_v = bi%recloc_v
     allocate(rids_w(n,3));    rids_w = bi%recloc_w
     allocate(rids_c(n,3));    rids_c = bi%recloc_c
     allocate(skip(n));        skip   = bi%lskipsec_loc
     allocate(comp(n));        comp   = bi%lcomprec_loc

   end subroutine mirror_sections


   !> Check that the device mirrors carry the host cache unchanged.
   !!
   !! The cache itself is checked against the expressions it replaced by
   !! check_wallfun_cache, which is host code and runs on a CPU build too. All
   !! that is left here is that mirror_sections copied it faithfully.
   subroutine check_ibm_section_cache(problem)
     implicit none
     character(len=*), intent(out) :: problem

     call check_wallfun_cache(problem)
     if (problem /= '') return

     if (allocated(u_secdist_d)) &
       call check_one('u', bound_info_u, bound_info_u%nfctsecsrank, &
                      u_secdist_d, u_recw_d, u_wfc_d, u_secbndloc_d, &
                      u_recids_u_d, u_recids_v_d, u_recids_w_d, u_recids_c_d)
     if (allocated(v_secdist_d)) &
       call check_one('v', bound_info_v, bound_info_v%nfctsecsrank, &
                      v_secdist_d, v_recw_d, v_wfc_d, v_secbndloc_d, &
                      v_recids_u_d, v_recids_v_d, v_recids_w_d, v_recids_c_d)
     if (allocated(w_secdist_d)) &
       call check_one('w', bound_info_w, bound_info_w%nfctsecsrank, &
                      w_secdist_d, w_recw_d, w_wfc_d, w_secbndloc_d, &
                      w_recids_u_d, w_recids_v_d, w_recids_w_d, w_recids_c_d)
     if (allocated(c_secdist_d)) &
       call check_one('c', bound_info_c, bound_info_c%nfctsecsrank, &
                      c_secdist_d, c_recw_d, c_wfc_d, c_secbndloc_d, &
                      c_recids_u_d, c_recids_v_d, c_recids_w_d, c_recids_c_d)

   contains

     subroutine check_one(tag, bi, n, dst, recw, wfc, bndloc, ru, rv, rw, rc)
       implicit none
       character(len=*), intent(in) :: tag
       type(bound_info_type), intent(in) :: bi
       integer, intent(in) :: n
       real,    device, intent(in) :: dst(n), recw(n,RECW_N), wfc(n,WFC_N)
       integer, device, intent(in) :: bndloc(n,3)
       integer, device, intent(in) :: ru(n,3), rv(n,3), rw(n,3), rc(n,3)

       real,    allocatable :: dst_b(:), recw_b(:,:), wfc_b(:,:)
       integer, allocatable :: bnd_b(:,:), ru_b(:,:), rv_b(:,:), rw_b(:,:), rc_b(:,:)

       if (problem /= '') return
       if (n < 1) return

       allocate(dst_b(n), recw_b(n,RECW_N), wfc_b(n,WFC_N))
       allocate(bnd_b(n,3), ru_b(n,3), rv_b(n,3), rw_b(n,3), rc_b(n,3))
       dst_b = dst; recw_b = recw; wfc_b = wfc
       bnd_b = bndloc; ru_b = ru; rv_b = rv; rw_b = rw; rc_b = rc

       if (any(dst_b  /= bi%secdist_loc)) call flag(tag, 'wall distance mirror')
       if (any(recw_b /= bi%recw_loc))    call flag(tag, 'trilinear offset mirror')
       if (any(wfc_b  /= bi%wfc_loc))     call flag(tag, 'roughness term mirror')
       if (any(bnd_b  /= bi%bndloc_loc))  call flag(tag, 'boundary point mirror')
       if (any(ru_b   /= bi%recloc_u))    call flag(tag, 'u reconstruction cell mirror')
       if (any(rv_b   /= bi%recloc_v))    call flag(tag, 'v reconstruction cell mirror')
       if (any(rw_b   /= bi%recloc_w))    call flag(tag, 'w reconstruction cell mirror')
       if (any(rc_b   /= bi%recloc_c))    call flag(tag, 'c reconstruction cell mirror')

       deallocate(dst_b, recw_b, wfc_b, bnd_b, ru_b, rv_b, rw_b, rc_b)

     end subroutine check_one

     subroutine flag(tag, what)
       implicit none
       character(len=*), intent(in) :: tag, what

       if (problem == '') problem = 'ibm section cache ('//tag//'): '//what
     end subroutine flag

   end subroutine check_ibm_section_cache



   !> Device counterpart of wallfunmom, for one staggered direction.
   !!
   !! align is 1, 2 or 3 for the u, v and w sections, and replaces the procedure
   !! pointers the host version uses to pick an interpolation stencil: those
   !! cannot be called from a kernel. Because the host passes dir = xhat, yhat
   !! or zhat, every dot_product with dir reduces to a single component, which
   !! is what align indexes here. Those reductions are exact, not approximate:
   !! multiplying by an exact 0 or 1 and summing changes no bits.
   !!
   !! Both scatters need atomics. Several facet sections can share one boundary
   !! point, and many share a facet, so neither rhs nor fac_tau is written by a
   !! single iteration only.
   !!
   !! Everything the host recomputes per section from fixed geometry - the local
   !! indices, the wall distance, the trilinear offsets and the roughness terms
   !! - is read from the caches mirror_sections builds. What is left depends on
   !! the solution: the interpolation, the local frame and the stress.
   subroutine wallfunmom_device(align, rhs, nsec, nfct, &
                                secareas, secfacids, secbndloc, secdist, recw, wfc, &
                                recids_u, recids_v, recids_w, recids_c, &
                                lskipsec, lcomprec, fac_tau)
     use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh, &
                           iwallmom, eps1, grav, prandtlturb
     use modcuda,   only : u0_d, v0_d, w0_d, thl0_d, dxdydzfi_d
     implicit none

     integer, intent(in) :: align, nsec, nfct
     real,    device, intent(inout) :: rhs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
     real,    device, intent(in)    :: secareas(nsec), secdist(nsec)
     real,    device, intent(in)    :: recw(nsec,RECW_N), wfc(nsec,WFC_N)
     integer, device, intent(in)    :: secfacids(nsec), secbndloc(nsec,3)
     integer, device, intent(in)    :: recids_u(nsec,3), recids_v(nsec,3)
     integer, device, intent(in)    :: recids_w(nsec,3), recids_c(nsec,3)
     logical, device, intent(in)    :: lskipsec(nsec), lcomprec(nsec)
     real,    device, intent(inout) :: fac_tau(nfct)

     integer :: sec, fac, i, j, k, iwm
     logical :: norec
     real    :: area, dist, stress, stress_dir, momvol, Tair, utan, ctm
     real    :: ux, uy, uz, nx, ny, nz, sx, sy, sz, tx, ty, tz, smagi
     real    :: a_is, sig, logdz

     if (nsec < 1) return

     ! Hoisted: reading these from inside the region would pull whole modules,
     ! and bound_info in particular, into the kernel parameter space.
     norec = lnorec
     iwm = iwallmom

     !$acc parallel loop default(present) &
     !$acc   private(fac, i, j, k, area, dist, stress, stress_dir, momvol, &
     !$acc           Tair, utan, ctm, ux, uy, uz, nx, ny, nz, sx, sy, sz, &
     !$acc           tx, ty, tz, smagi, a_is, sig, logdz)
     do sec = 1, nsec
       if (lskipsec(sec)) cycle

       ! log(dist/z0) is fixed geometry, so it is read rather than computed, and
       ! the test it gates is hoisted above the interpolation: the host tests it
       ! afterwards, but the interpolation has no side effect, so a section that
       ! fails here skips the whole gather instead of only the stress.
       logdz = wfc(sec,WFC_LOGDZ)
       if (logdz <= 1.) cycle

       area = secareas(sec)
       fac  = secfacids(sec)
       nx   = facnorm_d(fac,1)
       ny   = facnorm_d(fac,2)
       nz   = facnorm_d(fac,3)

       i = secbndloc(sec,1)
       j = secbndloc(sec,2)
       k = secbndloc(sec,3)

       dist = secdist(sec)
       Tair = 0.
       if (lcomprec(sec) .or. norec) then
         ! Velocity at the boundary point itself.
         select case (align)
         case (1)
           ux = u0_d(i,j,k)
           uy = 0.25 * (v0_d(i,j,k) + v0_d(i,j+1,k) + v0_d(i-1,j,k) + v0_d(i-1,j+1,k))
           uz = 0.25 * (w0_d(i,j,k) + w0_d(i,j,k+1) + w0_d(i-1,j,k) + w0_d(i-1,j,k+1))
           if (iwm == 2) Tair = 0.5 * (thl0_d(i  ,j,k)*mask_c_d(i  ,j,k)*(2.-mask_c_d(i-1,j,k)) &
                                    +  thl0_d(i-1,j,k)*mask_c_d(i-1,j,k)*(2.-mask_c_d(i  ,j,k)))
         case (2)
           ux = 0.25 * (u0_d(i,j,k) + u0_d(i+1,j,k) + u0_d(i,j-1,k) + u0_d(i+1,j-1,k))
           uy = v0_d(i,j,k)
           uz = 0.25 * (w0_d(i,j,k) + w0_d(i,j,k+1) + w0_d(i,j-1,k) + w0_d(i,j-1,k+1))
           if (iwm == 2) Tair = 0.5 * (thl0_d(i,j  ,k)*mask_c_d(i,j  ,k)*(2.-mask_c_d(i,j-1,k)) &
                                     + thl0_d(i,j-1,k)*mask_c_d(i,j-1,k)*(2.-mask_c_d(i,j  ,k)))
         case default
           ! Deliberately the same stencil as case (2): interp_velocity_w on the
           ! host is a copy of interp_velocity_v. Reproduced rather than
           ! corrected, so that the two paths agree; see the note in the commit.
           ux = 0.25 * (u0_d(i,j,k) + u0_d(i+1,j,k) + u0_d(i,j-1,k) + u0_d(i+1,j-1,k))
           uy = v0_d(i,j,k)
           uz = 0.25 * (w0_d(i,j,k) + w0_d(i,j,k+1) + w0_d(i,j-1,k) + w0_d(i,j-1,k+1))
           if (iwm == 2) Tair = 0.5 * (thl0_d(i,j,k  )*mask_c_d(i,j,k  )*(2.-mask_c_d(i,j,k-1)) &
                                    +  thl0_d(i,j,k-1)*mask_c_d(i,j,k-1)*(2.-mask_c_d(i,j,k  )))
         end select
       else
         ux = trilinear_device(u0_d, recids_u(sec,1), recids_u(sec,2), recids_u(sec,3), &
                               recw(sec,RECW_U), recw(sec,RECW_U+1), recw(sec,RECW_U+2))
         uy = trilinear_device(v0_d, recids_v(sec,1), recids_v(sec,2), recids_v(sec,3), &
                               recw(sec,RECW_V), recw(sec,RECW_V+1), recw(sec,RECW_V+2))
         uz = trilinear_device(w0_d, recids_w(sec,1), recids_w(sec,2), recids_w(sec,3), &
                               recw(sec,RECW_W), recw(sec,RECW_W+1), recw(sec,RECW_W+2))
         if (iwm == 2) &
           Tair = trilinear_device(thl0_d, recids_c(sec,1), recids_c(sec,2), recids_c(sec,3), &
                                   recw(sec,RECW_C), recw(sec,RECW_C+1), recw(sec,RECW_C+2))
       end if

       if (abs(ux) < eps1 .and. abs(uy) < eps1 .and. abs(uz) < eps1) cycle

       ! local_coords: span = norm x uvec, normalised; strm = span x norm.
       ! One reciprocal and three multiplies rather than three divisions, which
       ! is what local_coords does on the host too - in double precision on a
       ! 1:64 FP64 device a division is not a cheap operation.
       sx = ny*uz - nz*uy
       sy = nz*ux - nx*uz
       sz = nx*uy - ny*ux
       if (abs(sx) < eps1 .and. abs(sy) < eps1 .and. abs(sz) < eps1) cycle
       smagi = 1./sqrt(sx*sx + sy*sy + sz*sz)
       sx = sx * smagi
       sy = sy * smagi
       sz = sz * smagi
       tx = sy*nz - sz*ny
       ty = sz*nx - sx*nz
       tz = sx*ny - sy*nx

       utan = ux*tx + uy*ty + uz*tz

       ctm = 0.
       if (iwm == 2) then
         ctm = mom_transfer_coef_stability_device(utan, dist, Tair, facT1_d(fac), &
                                                  logdz, wfc(sec,WFC_LOGZH), wfc(sec,WFC_COEF), &
                                                  wfc(sec,WFC_CM), wfc(sec,WFC_CH), &
                                                  grav, prandtlturb)
       else if (iwm == 3) then
         ctm = wfc(sec,WFC_CTMNEU)
       end if

       stress = ctm * utan**2

       ! dir is a unit axis vector, so dot_product(dir, strm) is strm(align) and
       ! dot_product(xhat, norm) is norm(1), and so on.
       if (align == 1) then
         a_is = tx
         sig  = ux
       else if (align == 2) then
         a_is = ty
         sig  = uy
       else
         a_is = tz
         sig  = uz
       end if

       if (lcomprec(sec)) then
         stress_dir = a_is * stress
       else
         stress_dir = sqrt((a_is*nx*stress)**2 + (a_is*ny*stress)**2 + (a_is*nz*stress)**2)
       end if

       stress_dir = sign(stress_dir, sig)

       momvol = stress_dir * area * dxdydzfi_d(k)

       !$acc atomic update
       rhs(i,j,k) = rhs(i,j,k) - momvol

       !$acc atomic update
       fac_tau(fac) = fac_tau(fac) + stress_dir * area
     end do
     !$acc end parallel loop

   end subroutine wallfunmom_device
   !> Trilinear interpolation of a device field at a reconstruction point.
   !!
   !! i, j, k are the rank-local indices of the containing cell and xd, yd, zd
   !! the normalised offsets within it. Both are precomputed per section, so
   !! this is an eight-point gather and nothing else: the grid arrays, the
   !! global-to-local shift and the three divisions the offsets cost are all
   !! resolved once, at initialisation.
   real function trilinear_device(var, i, j, k, xd, yd, zd)
     !$acc routine seq
     use modcuda, only : ib_d, ie_d, jb_d, je_d, kb_d, ke_d, ih_d, jh_d, kh_d
     implicit none

     real, device, intent(in) :: var(ib_d-ih_d:ie_d+ih_d, jb_d-jh_d:je_d+jh_d, kb_d-kh_d:ke_d+kh_d)
     integer, value, intent(in) :: i, j, k
     real,    value, intent(in) :: xd, yd, zd

     trilinear_device = var(i  ,j  ,k  ) * (1-xd)*(1-yd)*(1-zd) + &
                        var(i+1,j  ,k  ) * (  xd)*(1-yd)*(1-zd) + &
                        var(i  ,j+1,k  ) * (1-xd)*(  yd)*(1-zd) + &
                        var(i+1,j+1,k  ) * (  xd)*(  yd)*(1-zd) + &
                        var(i  ,j  ,k+1) * (1-xd)*(1-yd)*(  zd) + &
                        var(i+1,j  ,k+1) * (  xd)*(1-yd)*(  zd) + &
                        var(i  ,j+1,k+1) * (1-xd)*(  yd)*(  zd) + &
                        var(i+1,j+1,k+1) * (  xd)*(  yd)*(  zd)

   end function trilinear_device

   !> Device copy of mom_transfer_coef_stability.
   !!
   !! The roughness terms are passed in already evaluated, from the per-section
   !! cache built at initialisation: logdz = log(dist/z0), logzh = log(z0/z0h),
   !! coef = fkar2/logdz**2, and cmfac, chfac the two (d*fkar2)/logdz**2*b1*sqdz
   !! products. None of them depends on the solution, and together they are the
   !! only logarithms and the only square root the original routine evaluated.
   !! The remaining module constants are passed rather than used, so the routine
   !! carries no module state into the kernel.
   real function mom_transfer_coef_stability_device(utan, dist, Tair, Tsurf, &
                                                    logdz, logzh, coef, cmfac, chfac, &
                                                    grav_in, prandtlturb_in)
     !$acc routine seq
     implicit none
     real, value, intent(in) :: utan, dist, Tair, Tsurf
     real, value, intent(in) :: logdz, logzh, coef, cmfac, chfac
     real, value, intent(in) :: grav_in, prandtlturb_in

     real, parameter :: b1 = 9.4
     real, parameter :: b2 = 4.7
     real :: dT, Ribl0, Ribl1, Fm, Fh, M

     dT = Tair - Tsurf
     Ribl0 = grav_in * dist * dT / (Tsurf * utan**2)

     if (Ribl0 > 0.) then
        Fm = 1./(1. + b2*Ribl0)**2
        Fh = Fm
     else
        Fm = 1. - (b1*Ribl0)/(1. + cmfac*sqrt(abs(Ribl0)))
        Fh = 1. - (b1*Ribl0)/(1. + chfac*sqrt(abs(Ribl0)))
     end if

     M = prandtlturb_in*logdz*sqrt(Fm)/Fh

     Ribl1 = Ribl0 - Ribl0*prandtlturb_in*logzh/(prandtlturb_in*logzh + M)

     if (Ribl1 > 0.) then
        Fm = 1./(1. + b2*Ribl1)**2
     else
        Fm = 1. - (b1*Ribl1)/(1. + cmfac*sqrt(abs(Ribl1)))
     end if

     mom_transfer_coef_stability_device = coef*Fm

   end function mom_transfer_coef_stability_device

   !> Run wallfunmom_device for one direction: 1 = u, 2 = v, 3 = w.
   !!
   !! fac_tau_d is zeroed here rather than by the caller, matching the host
   !! version where fac_tau_loc is a local reset on entry.
   subroutine wallfunmom_dir_device(dirsel, rhs)
     use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh, nfcts
     implicit none

     integer, intent(in) :: dirsel
     real, device, intent(inout) :: rhs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)

     select case (dirsel)
     case (1)
       if (.not. allocated(u_secareas_d)) return
     case (2)
       if (.not. allocated(v_secareas_d)) return
     case (3)
       if (.not. allocated(w_secareas_d)) return
     end select

     fac_tau_d = 0.

     select case (dirsel)
     case (1)
       call wallfunmom_device(1, rhs, bound_info_u%nfctsecsrank, nfcts, &
                              u_secareas_d, u_secfacids_d, u_secbndloc_d, u_secdist_d, &
                              u_recw_d, u_wfc_d, &
                              u_recids_u_d, u_recids_v_d, u_recids_w_d, u_recids_c_d, &
                              u_lskipsec_d, u_lcomprec_d, fac_tau_d)
     case (2)
       call wallfunmom_device(2, rhs, bound_info_v%nfctsecsrank, nfcts, &
                              v_secareas_d, v_secfacids_d, v_secbndloc_d, v_secdist_d, &
                              v_recw_d, v_wfc_d, &
                              v_recids_u_d, v_recids_v_d, v_recids_w_d, v_recids_c_d, &
                              v_lskipsec_d, v_lcomprec_d, fac_tau_d)
     case (3)
       call wallfunmom_device(3, rhs, bound_info_w%nfctsecsrank, nfcts, &
                              w_secareas_d, w_secfacids_d, w_secbndloc_d, w_secdist_d, &
                              w_recw_d, w_wfc_d, &
                              w_recids_u_d, w_recids_v_d, w_recids_w_d, w_recids_c_d, &
                              w_lskipsec_d, w_lcomprec_d, fac_tau_d)
     end select

   end subroutine wallfunmom_dir_device


   !> Device counterpart of wallfunheat.
   !!
   !! Structured like wallfunmom_device: the cell-centred sections, atomics on
   !! every per-facet and per-cell scatter, and the two whole-domain sums taken
   !! as loop reductions rather than running totals.
   !!
   !! flux, htc and cth are reset at the top of each section. On the host they
   !! are procedure-local and therefore survive between iterations, which the
   !! code relies on in configurations it does not really support: with
   !! iwalltemp = 1 and a facet normal matching none of the five axis cases,
   !! flux keeps the previous section's value, and the moisture branch reads htc
   !! from the sensible-heat block that only runs under ltempeq. Neither is
   !! reproducible in a parallel loop. Where the host is well defined - ltempeq
   !! with iwalltemp = 2, which sets both every iteration - the two agree; where
   !! it is not, this contributes zero instead of an arbitrary carried value.
   subroutine wallfunheat_device(nsec, nfct, hflux_sec, qflux_sec, &
                                 secareas, secfacids, secbndloc, secdist, recw, wfc, &
                                 recids_u, recids_v, recids_w, recids_c, &
                                 lskipsec, lcomprec, &
                                 f_htc, f_cth, f_pres, f_pres2, f_hf, f_ef, &
                                 totheat_out, totq_out)
     use modglobal,  only : ib, ie, ih, jb, je, jh, kb, ke, kh, &
                            eps1, grav, prandtlturb, &
                            ltempeq, lmoist, iwalltemp, iwallmoist, lEB
     use modcuda,    only : u0_d, v0_d, w0_d, thl0_d, qt0_d, pres0_d, dxdydzhi_d
     use modibmdata, only : bctfxm, bctfxp, bctfyp, bctfz, &
                            bcqfxm, bcqfxp, bcqfyp, bcqfym, bcqfz
     implicit none

     integer, intent(in) :: nsec, nfct
     real,    device, intent(out)   :: hflux_sec(nsec), qflux_sec(nsec)
     real,    device, intent(in)    :: secareas(nsec), secdist(nsec)
     real,    device, intent(in)    :: recw(nsec,RECW_N), wfc(nsec,WFC_N)
     integer, device, intent(in)    :: secfacids(nsec), secbndloc(nsec,3)
     integer, device, intent(in)    :: recids_u(nsec,3), recids_v(nsec,3)
     integer, device, intent(in)    :: recids_w(nsec,3), recids_c(nsec,3)
     logical, device, intent(in)    :: lskipsec(nsec), lcomprec(nsec)
     real,    device, intent(inout) :: f_htc(nfct), f_cth(nfct), f_pres(nfct), f_pres2(nfct)
     real,    device, intent(inout) :: f_hf(0:nfct), f_ef(nfct)
     real,    intent(out) :: totheat_out, totq_out

     integer :: sec, fac, i, j, k, iwt, iwq
     logical :: norec, ltemp, lmoi, leb_l
     real    :: area, dist, Tair, qtair, utan, cth, htc, flux
     real    :: ux, uy, uz, nx, ny, nz, sx, sy, sz, tx, ty, tz, smagi
     real    :: p0
     real    :: cveg, hurel, qwall, resa, resc, ress
     real    :: dT, Ribl0, Ribl1, logdz, logzh, coef, cmfac, chfac
     real    :: Fm, Fh, Mst, dTrough
     real    :: totheat, totq

     real, parameter :: b1 = 9.4
     real, parameter :: b2 = 4.7

     totheat_out = 0.
     totq_out    = 0.
     if (nsec < 1) return

     norec = lnorec
     iwt   = iwalltemp
     iwq   = iwallmoist
     ltemp = ltempeq
     lmoi  = lmoist
     leb_l = lEB

     totheat = 0.
     totq    = 0.

     !$acc parallel loop default(present)
     do sec = 1, nsec
       hflux_sec(sec) = 0.
       qflux_sec(sec) = 0.
     end do
     !$acc end parallel loop

     !$acc parallel loop default(present) reduction(+:totheat,totq) &
     !$acc   private(fac, i, j, k, area, dist, Tair, qtair, utan, cth, htc, flux, &
     !$acc           ux, uy, uz, nx, ny, nz, sx, sy, sz, tx, ty, tz, smagi, &
     !$acc           p0, cveg, hurel, qwall, resa, resc, ress, &
     !$acc           dT, Ribl0, Ribl1, logdz, logzh, coef, cmfac, chfac, &
     !$acc           Fm, Fh, Mst, dTrough)
     do sec = 1, nsec
       fac  = secfacids(sec)
       area = secareas(sec)

       i = secbndloc(sec,1)
       j = secbndloc(sec,2)
       k = secbndloc(sec,3)

       ! Accumulated for every section, including the skipped ones, as on the host.
       p0 = pres0_d(i,j,k)
       !$acc atomic update
       f_pres(fac) = f_pres(fac) + p0 * area
       !$acc atomic update
       f_pres2(fac) = f_pres2(fac) + p0 * p0 * area

       if (lskipsec(sec)) cycle

       ! log(dist/z0) is fixed geometry, so it is read rather than computed, and
       ! the test it gates is hoisted above the interpolation: the host tests it
       ! afterwards, but the interpolation has no side effect, so a section that
       ! fails here skips the whole gather instead of only the flux.
       logdz = wfc(sec,WFC_LOGDZ)
       if (logdz <= 1.) cycle

       nx = facnorm_d(fac,1)
       ny = facnorm_d(fac,2)
       nz = facnorm_d(fac,3)

       dist  = secdist(sec)
       qtair = 0.
       if (lcomprec(sec) .or. norec) then
         ux = 0.5 * (u0_d(i,j,k) + u0_d(i+1,j,k))
         uy = 0.5 * (v0_d(i,j,k) + v0_d(i,j+1,k))
         uz = 0.5 * (w0_d(i,j,k) + w0_d(i,j,k+1))
         Tair  = thl0_d(i,j,k)
         ! Guarded because qt0_d exists only when lmoist is set, while the host
         ! qt0 is always allocated and so is read unconditionally there. qtair
         ! is used only inside the moisture branch, so the two still agree.
         if (lmoi) qtair = qt0_d(i,j,k)
       else
         ux = trilinear_device(u0_d, recids_u(sec,1), recids_u(sec,2), recids_u(sec,3), &
                               recw(sec,RECW_U), recw(sec,RECW_U+1), recw(sec,RECW_U+2))
         uy = trilinear_device(v0_d, recids_v(sec,1), recids_v(sec,2), recids_v(sec,3), &
                               recw(sec,RECW_V), recw(sec,RECW_V+1), recw(sec,RECW_V+2))
         uz = trilinear_device(w0_d, recids_w(sec,1), recids_w(sec,2), recids_w(sec,3), &
                               recw(sec,RECW_W), recw(sec,RECW_W+1), recw(sec,RECW_W+2))
         Tair = trilinear_device(thl0_d, recids_c(sec,1), recids_c(sec,2), recids_c(sec,3), &
                                 recw(sec,RECW_C), recw(sec,RECW_C+1), recw(sec,RECW_C+2))
         if (lmoi) &
           qtair = trilinear_device(qt0_d, recids_c(sec,1), recids_c(sec,2), recids_c(sec,3), &
                                    recw(sec,RECW_C), recw(sec,RECW_C+1), recw(sec,RECW_C+2))
       end if

       if (abs(ux) < eps1 .and. abs(uy) < eps1 .and. abs(uz) < eps1) cycle

       sx = ny*uz - nz*uy
       sy = nz*ux - nx*uz
       sz = nx*uy - ny*ux
       if (abs(sx) < eps1 .and. abs(sy) < eps1 .and. abs(sz) < eps1) cycle
       smagi = 1./sqrt(sx*sx + sy*sy + sz*sz)
       sx = sx * smagi
       sy = sy * smagi
       sz = sz * smagi
       tx = sy*nz - sz*ny
       ty = sz*nx - sx*nz
       tz = sx*ny - sy*nx

       utan = ux*tx + uy*ty + uz*tz

       flux = 0.
       htc  = 0.
       cth  = 0.

       if (ltemp) then
         if (iwt == 1) then
           if      (abs(nx-1.) < eps1 .and. abs(ny) < eps1 .and. abs(nz) < eps1) then
             flux = bctfxp
           else if (abs(nx+1.) < eps1 .and. abs(ny) < eps1 .and. abs(nz) < eps1) then
             flux = bctfxm
           else if (abs(ny-1.) < eps1 .and. abs(nx) < eps1 .and. abs(nz) < eps1) then
             flux = bctfyp
           else if (abs(ny+1.) < eps1 .and. abs(nx) < eps1 .and. abs(nz) < eps1) then
             ! bctfxm, not bctfym: reproduced from the host, where it looks like
             ! a typo but changing it here would break agreement between them.
             flux = bctfxm
           else if (abs(nz-1.) < eps1 .and. abs(nx) < eps1 .and. abs(ny) < eps1) then
             flux = bctfz
           end if
         else if (iwt == 2) then
           ! heat_transfer_coef_flux, inlined. It returns three values, which is
           ! awkward to carry across a routine seq boundary, and it has a single
           ! call site. The four roughness terms below are read from wfc rather
           ! than rebuilt: they cost two logarithms, a square root and three
           ! divisions, and none of them depends on the solution.
           logzh = wfc(sec,WFC_LOGZH)
           coef  = wfc(sec,WFC_COEF)
           cmfac = wfc(sec,WFC_CM)
           chfac = wfc(sec,WFC_CH)

           dT = Tair - facT1_d(fac)
           Ribl0 = grav * dist * dT / (facT1_d(fac) * utan**2)
           if (Ribl0 > 0.) then
             Fm = 1./(1. + b2*Ribl0)**2
             Fh = Fm
           else
             Fm = 1. - (b1*Ribl0)/(1. + cmfac*sqrt(abs(Ribl0)))
             Fh = 1. - (b1*Ribl0)/(1. + chfac*sqrt(abs(Ribl0)))
           end if
           Mst = prandtlturb*logdz*sqrt(Fm)/Fh
           Ribl1 = Ribl0 - Ribl0*prandtlturb*logzh/(prandtlturb*logzh + Mst)
           if (Ribl1 > 0.) then
             Fm = 1./(1. + b2*Ribl1)**2
             Fh = Fm
           else
             Fm = 1. - (b1*Ribl1)/(1. + cmfac*sqrt(abs(Ribl1)))
             Fh = 1. - (b1*Ribl1)/(1. + chfac*sqrt(abs(Ribl1)))
           end if
           Mst = prandtlturb*logdz*sqrt(Fm)/Fh
           dTrough = dT*1./(prandtlturb*logzh/Mst + 1.)
           cth = coef*Fh/prandtlturb
           flux = abs(utan)*cth*dTrough
           if (abs(abs(utan)*dT) > 0.) then
             htc = flux / (abs(utan)*dT)
           else
             htc = 0.
           end if

           !$acc atomic update
           f_cth(fac) = f_cth(fac) + cth * area
           !$acc atomic update
           f_htc(fac) = f_htc(fac) + htc * area
         end if

         hflux_sec(sec) = flux * area * dxdydzhi_d(k)

         totheat = totheat + flux*area

         if (leb_l) then
           !$acc atomic update
           f_hf(fac) = f_hf(fac) + flux * area
         end if
       end if

       if (lmoi) then
         if (faclGR_d(fac)) then
           if (iwq == 1) then
             if      (abs(nx-1.) < eps1 .and. abs(ny) < eps1 .and. abs(nz) < eps1) then
               flux = bcqfxp
             else if (abs(nx+1.) < eps1 .and. abs(ny) < eps1 .and. abs(nz) < eps1) then
               flux = bcqfxm
             else if (abs(ny-1.) < eps1 .and. abs(nx) < eps1 .and. abs(nz) < eps1) then
               flux = bcqfyp
             else if (abs(ny+1.) < eps1 .and. abs(nx) < eps1 .and. abs(nz) < eps1) then
               flux = bcqfym
             else if (abs(nz-1.) < eps1 .and. abs(nx) < eps1 .and. abs(ny) < eps1) then
               flux = bcqfz
             end if
           else if (iwq == 2) then
             if (abs(htc*abs(utan)) > 0.) then
               qwall = facqsat_d(fac)
               hurel = fachurel_d(fac)
               resa  = 1./(htc*abs(utan))
               resc  = facf_d(fac,4)
               ress  = facf_d(fac,5)
               cveg  = 0.8
               flux  = min(0., cveg * (qtair - qwall)         / (resa + resc) + &
                           (1 - cveg)* (qtair - qwall * hurel) / (resa + ress))
             end if
           end if

           totq = totq + flux*area

           qflux_sec(sec) = flux * area * dxdydzhi_d(k)

           if (leb_l) then
             !$acc atomic update
             f_ef(fac) = f_ef(fac) + flux * area
           end if
         end if
       end if
     end do
     !$acc end parallel loop

     totheat_out = totheat
     totq_out    = totq

   end subroutine wallfunheat_device

   !> Scatter per-section tendency increments into a field.
   !!
   !! A separate pass so that wallfunheat_device is handed neither thlp nor qtp.
   !! Neither exists on the device unless its own switch is set, and passing an
   !! unallocated device array faults the kernel even where no branch reads it.
   subroutine scatter_flux_device(nsec, qtp_dev, qflux_sec, secbndloc)
     use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh
     implicit none

     integer, intent(in) :: nsec
     real,    device, intent(inout) :: qtp_dev(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
     real,    device, intent(in)    :: qflux_sec(nsec)
     integer, device, intent(in)    :: secbndloc(nsec,3)

     integer :: sec, i, j, k

     if (nsec < 1) return

     !$acc parallel loop default(present) private(i, j, k)
     do sec = 1, nsec
       if (qflux_sec(sec) == 0.) cycle
       i = secbndloc(sec,1)
       j = secbndloc(sec,2)
       k = secbndloc(sec,3)
       !$acc atomic update
       qtp_dev(i,j,k) = qtp_dev(i,j,k) - qflux_sec(sec)
     end do
     !$acc end parallel loop

   end subroutine scatter_flux_device


   !> Run wallfunheat_device against the cell-centred sections.
   !!
   !! Zeroes the per-call accumulators, runs the kernel and folds the two
   !! whole-domain sums into their host counterparts, so the caller sees the
   !! same state the host wallfunheat would leave behind.
   subroutine wallfunheat_dir_device
     use modglobal, only : nfcts, totheatflux, totqflux, lmoist, ltempeq
     use modcuda,   only : thlp_d, qtp_d
     implicit none

     real :: totheat, totq

     if (.not. allocated(c_secareas_d)) return

     fac_htc_d = 0.
     fac_cth_d = 0.
     fac_pres_d = 0.
     fac_pres2_d = 0.
     ! fachf_d and facef_d are deliberately not reset here: modEB consumes them
     ! once per time step, so they accumulate across the Runge-Kutta stages and
     ! updateFacFluxHost zeroes them when it hands them over.

     call wallfunheat_device(bound_info_c%nfctsecsrank, nfcts, hflux_sec_d, qflux_sec_d, &
                             c_secareas_d, c_secfacids_d, c_secbndloc_d, c_secdist_d, &
                             c_recw_d, c_wfc_d, &
                             c_recids_u_d, c_recids_v_d, c_recids_w_d, c_recids_c_d, &
                             c_lskipsec_d, c_lcomprec_d, &
                             fac_htc_d, fac_cth_d, fac_pres_d, fac_pres2_d, fachf_d, facef_d, &
                             totheat, totq)

     if (ltempeq) call scatter_flux_device(bound_info_c%nfctsecsrank, thlp_d, &
                                          hflux_sec_d, c_secbndloc_d)
     if (lmoist)  call scatter_flux_device(bound_info_c%nfctsecsrank, qtp_d, &
                                          qflux_sec_d, c_secbndloc_d)

     totheatflux = totheatflux + totheat
     totqflux    = totqflux + totq

   end subroutine wallfunheat_dir_device


   !> dst = src over the tendency extent.
   subroutine copy_tendency_device(dst, src)
     use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh
     implicit none
     real, device, intent(out) :: dst(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
     real, device, intent(in)  :: src(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
     integer :: i, j, k

     !$acc parallel loop collapse(3) default(present)
     do k = kb, ke+kh
       do j = jb-jh, je+jh
         do i = ib-ih, ie+ih
           dst(i,j,k) = src(i,j,k)
         end do
       end do
     end do
     !$acc end parallel loop

   end subroutine copy_tendency_device

   !> acc = acc + (cur - prev), over the same k range the host slices.
   subroutine accumulate_delta_device(acc, cur, prev)
     use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh
     implicit none
     real, device, intent(inout) :: acc(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
     real, device, intent(in)    :: cur(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
     real, device, intent(in)    :: prev(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
     integer :: i, j, k

     !$acc parallel loop collapse(3) default(present)
     do k = kb, ke+kh
       do j = jb-jh, je+jh
         do i = ib-ih, ie+ih
           acc(i,j,k) = acc(i,j,k) + (cur(i,j,k) - prev(i,j,k))
         end do
       end do
     end do
     !$acc end parallel loop

   end subroutine accumulate_delta_device

   !> Facet-stress output for one direction, mirroring the tail of wallfunmom.
   subroutine reduce_fac_tau_device(dirsel)
     use modglobal, only : nfcts, lwritefac, rk3step
     use initfac,   only : faca
     use modmpi,    only : comm3d, mpi_sum, mpierr, my_real
     implicit none

     integer, intent(in) :: dirsel
     real :: tot(1:nfcts)

     if (.not. (lwritefac .and. rk3step == 3)) return

     fac_stage = fac_tau_d
     fac_stage = fac_stage / faca(1:nfcts)
     call MPI_ALLREDUCE(fac_stage, tot, nfcts, MY_REAL, MPI_SUM, comm3d, mpierr)

     select case (dirsel)
     case (1)
       fac_tau_x = tot
     case (2)
       fac_tau_y = tot
     case (3)
       fac_tau_z = tot
     end select

   end subroutine reduce_fac_tau_device

   !> Facet heat output, mirroring the tail of wallfunheat, plus the energy
   !! balance accumulators, which are added on every call rather than only at
   !! the last Runge-Kutta stage.
   subroutine reduce_fac_heat_device
     use modglobal, only : nfcts, lwritefac, rk3step
     use initfac,   only : faca
     use modmpi,    only : comm3d, mpi_sum, mpierr, my_real
     implicit none

     real :: tot(1:nfcts)

     ! The energy balance accumulators are handed over by updateFacFluxHost, in
     ! modcuda, so that every loop-time transfer is declared in one place.
     if (.not. (lwritefac .and. rk3step == 3)) return

     fac_stage = fac_cth_d
     fac_stage = fac_stage / faca(1:nfcts)
     call MPI_ALLREDUCE(fac_stage, tot, nfcts, MY_REAL, MPI_SUM, comm3d, mpierr)
     fac_cth = tot

     fac_stage = fac_htc_d
     fac_stage = fac_stage / faca(1:nfcts)
     call MPI_ALLREDUCE(fac_stage, tot, nfcts, MY_REAL, MPI_SUM, comm3d, mpierr)
     fac_htc = tot

     fac_stage = fac_pres_d
     fac_stage = fac_stage / faca(1:nfcts)
     call MPI_ALLREDUCE(fac_stage, tot, nfcts, MY_REAL, MPI_SUM, comm3d, mpierr)
     fac_pres = tot

     fac_stage = fac_pres2_d
     fac_stage = fac_stage / faca(1:nfcts)
     call MPI_ALLREDUCE(fac_stage, tot, nfcts, MY_REAL, MPI_SUM, comm3d, mpierr)
     fac_pres2 = tot

   end subroutine reduce_fac_heat_device


   !> Device counterparts of diffu_corr, diffv_corr, diffw_corr and diffc_corr.
   !!
   !! Each cancels the subgrid tendency that modsubgrid added across a solid
   !! boundary, at this rank's fluid boundary points for one staggered location.
   !!
   !! The scatter into the tendency uses no atomics. That is only sound because
   !! the entries of a bndpts list are distinct cells, so no two iterations of
   !! the loop touch the same (i,j,k). Were a cell listed twice, the host loop
   !! would apply both corrections and the device loop would lose one.
   !! test_ibm_diff_corr asserts the distinctness rather than trusting it.
   subroutine diffu_corr_device
     use modglobal, only : eps1, dy2i
     use modcuda,   only : up_d, u0_d, ekm_d, dzf_d, dzfi_d, dzhi_d, dzhiq_d
     implicit none

     real    :: empo, emmo, emop, emom
     integer :: i, j, k, n, npts

     ! Held in a local: referencing bound_info_u%nbndptsrank inside the region
     ! pulls the whole derived type, allocatable descriptors and all, into the
     ! kernel parameter space and overflows its 4096-byte limit.
     if (.not. allocated(bndpts_u_d)) return
     npts = bound_info_u%nbndptsrank
     if (npts < 1) return

     !$acc parallel loop default(present) private(i, j, k, empo, emmo, emop, emom)
     do n = 1, npts
        i = bndpts_u_d(n,1)
        j = bndpts_u_d(n,2)
        k = bndpts_u_d(n,3)

        if (abs(mask_u_d(i,j+1,k)) < eps1) then
          empo = 0.25 * ((ekm_d(i,j,k) + ekm_d(i,j+1,k)) + (ekm_d(i-1,j,k) + ekm_d(i-1,j+1,k)))
          up_d(i,j,k) = up_d(i,j,k) - empo * (u0_d(i,j+1,k) - u0_d(i,j,k))*dy2i
        end if

        if (abs(mask_u_d(i,j-1,k)) < eps1) then
          emmo = 0.25 * ((ekm_d(i,j,k) + ekm_d(i,j-1,k)) + (ekm_d(i-1,j-1,k) + ekm_d(i-1,j,k)))
          up_d(i,j,k) = up_d(i,j,k) + emmo * (u0_d(i,j,k) - u0_d(i,j-1,k))*dy2i
        end if

        if (abs(mask_u_d(i,j,k+1)) < eps1) then
          emop = (dzf_d(k+1) * ( ekm_d(i,j,k)   + ekm_d(i-1,j,k  ))  + &
                  dzf_d(k)   * ( ekm_d(i,j,k+1) + ekm_d(i-1,j,k+1))) * dzhiq_d(k+1)
          up_d(i,j,k) = up_d(i,j,k) - emop * (u0_d(i,j,k+1) - u0_d(i,j,k))*dzhi_d(k+1)*dzfi_d(k)
        end if

        if (abs(mask_u_d(i,j,k-1)) < eps1) then
          emom = (dzf_d(k-1) * (ekm_d(i,j,k  ) + ekm_d(i-1,j,k  ))  + &
                  dzf_d(k)   * (ekm_d(i,j,k-1) + ekm_d(i-1,j,k-1))) * dzhiq_d(k)
          up_d(i,j,k) = up_d(i,j,k) + emom * (u0_d(i,j,k) - u0_d(i,j,k-1))*dzhi_d(k)*dzfi_d(k)
        end if
     end do
     !$acc end parallel loop

   end subroutine diffu_corr_device

   subroutine diffv_corr_device
     use modglobal, only : eps1, dx2i
     use modcuda,   only : vp_d, v0_d, ekm_d, dzf_d, dzfi_d, dzhi_d, dzhiq_d
     implicit none

     real    :: epmo, emmo, eomp, eomm
     integer :: i, j, k, n, npts

     ! Held in a local: referencing bound_info_v%nbndptsrank inside the region
     ! pulls the whole derived type, allocatable descriptors and all, into the
     ! kernel parameter space and overflows its 4096-byte limit.
     if (.not. allocated(bndpts_v_d)) return
     npts = bound_info_v%nbndptsrank
     if (npts < 1) return

     !$acc parallel loop default(present) private(i, j, k, epmo, emmo, eomp, eomm)
     do n = 1, npts
        i = bndpts_v_d(n,1)
        j = bndpts_v_d(n,2)
        k = bndpts_v_d(n,3)

        if (abs(mask_v_d(i+1,j,k)) < eps1) then
          epmo = 0.25 * (ekm_d(i,j,k) + ekm_d(i,j-1,k) + ekm_d(i+1,j-1,k) + ekm_d(i+1,j,k))
          vp_d(i,j,k) = vp_d(i,j,k) - epmo * (v0_d(i+1,j,k) - v0_d(i,j,k))*dx2i
        end if

        if (abs(mask_v_d(i-1,j,k)) < eps1) then
          emmo = 0.25 * (ekm_d(i,j,k) + ekm_d(i,j-1,k) + ekm_d(i-1,j-1,k) + ekm_d(i-1,j,k))
          vp_d(i,j,k) = vp_d(i,j,k) + emmo * (v0_d(i,j,k) - v0_d(i-1,j,k))*dx2i
        end if

        if (abs(mask_v_d(i,j,k+1)) < eps1) then
          eomp = ( dzf_d(k+1) * ( ekm_d(i,j,k)   + ekm_d(i,j-1,k)  )  + &
                   dzf_d(k  ) * ( ekm_d(i,j,k+1) + ekm_d(i,j-1,k+1))) * dzhiq_d(k+1)
          vp_d(i,j,k) = vp_d(i,j,k) - eomp * (v0_d(i,j,k+1) - v0_d(i,j,k))*dzhi_d(k+1)*dzfi_d(k)
        end if

        if (abs(mask_v_d(i,j,k-1)) < eps1) then
          eomm = ( dzf_d(k-1) * ( ekm_d(i,j,k  )  + ekm_d(i,j-1,k)   ) + &
                   dzf_d(k)   * ( ekm_d(i,j,k-1)  + ekm_d(i,j-1,k-1))) * dzhiq_d(k)
          vp_d(i,j,k) = vp_d(i,j,k) + eomm * (v0_d(i,j,k) - v0_d(i,j,k-1))*dzhi_d(k)*dzfi_d(k)
        end if
     end do
     !$acc end parallel loop

   end subroutine diffv_corr_device

   subroutine diffw_corr_device
     use modglobal, only : eps1, dx2i, dy2i
     use modcuda,   only : wp_d, w0_d, ekm_d, dzf_d, dzhiq_d
     implicit none

     real    :: epom, emom, eopm, eomm
     integer :: i, j, k, n, npts

     ! Held in a local: referencing bound_info_w%nbndptsrank inside the region
     ! pulls the whole derived type, allocatable descriptors and all, into the
     ! kernel parameter space and overflows its 4096-byte limit.
     if (.not. allocated(bndpts_w_d)) return
     npts = bound_info_w%nbndptsrank
     if (npts < 1) return

     !$acc parallel loop default(present) private(i, j, k, epom, emom, eopm, eomm)
     do n = 1, npts
        i = bndpts_w_d(n,1)
        j = bndpts_w_d(n,2)
        k = bndpts_w_d(n,3)

        if (abs(mask_w_d(i+1,j,k)) < eps1) then
          epom = ( dzf_d(k-1) * ( ekm_d(i,j,k  ) + ekm_d(i+1,j,k  ))    + &
                   dzf_d(k  ) * ( ekm_d(i,j,k-1) + ekm_d(i+1,j,k-1))) * dzhiq_d(k)
          wp_d(i,j,k) = wp_d(i,j,k) - epom * (w0_d(i+1,j,k) - w0_d(i,j,k))*dx2i
        end if

        if (abs(mask_w_d(i-1,j,k)) < eps1) then
          emom = ( dzf_d(k-1) * ( ekm_d(i,j,k  ) + ekm_d(i-1,j,k  ))  + &
                   dzf_d(k  ) * ( ekm_d(i,j,k-1) + ekm_d(i-1,j,k-1))) * dzhiq_d(k)
          wp_d(i,j,k) = wp_d(i,j,k) + emom * (w0_d(i,j,k) - w0_d(i-1,j,k))*dx2i
        end if

        if (abs(mask_w_d(i,j+1,k)) < eps1) then
          eopm = ( dzf_d(k-1) * ( ekm_d(i,j,k  ) + ekm_d(i,j+1,k  ))  + &
                   dzf_d(k  ) * ( ekm_d(i,j,k-1) + ekm_d(i,j+1,k-1))) * dzhiq_d(k)
          wp_d(i,j,k) = wp_d(i,j,k) - eopm * (w0_d(i,j+1,k) - w0_d(i,j,k))*dy2i
        end if

        if (abs(mask_w_d(i,j-1,k)) < eps1) then
          eomm = ( dzf_d(k-1) * ( ekm_d(i,j,k  ) + ekm_d(i,j-1,k  ))  + &
                   dzf_d(k  ) * ( ekm_d(i,j,k-1) + ekm_d(i,j-1,k-1))) * dzhiq_d(k)
          wp_d(i,j,k) = wp_d(i,j,k) + eomm * (w0_d(i,j,k) - w0_d(i,j-1,k))*dy2i
        end if
     end do
     !$acc end parallel loop

   end subroutine diffw_corr_device

   !> Scalar counterpart. var and rhs are device arrays; hi/hj/hk are their halo
   !! widths, which differ between thl/qt (ih,jh,kh) and sv (ihc,jhc,khc).
   subroutine diffc_corr_device(var, rhs, hi, hj, hk)
     use modglobal, only : eps1, ib, ie, jb, je, kb, ke, dx2i, dy2i
     use modcuda,   only : ekh_d, dzf_d, dzfi_d, dzh2i_d
     implicit none

     integer, intent(in) :: hi, hj, hk
     real, device, intent(in)    :: var(ib-hi:ie+hi,jb-hj:je+hj,kb-hk:ke+hk)
     real, device, intent(inout) :: rhs(ib-hi:ie+hi,jb-hj:je+hj,kb   :ke+hk)

     integer :: i, j, k, n, npts

     if (.not. allocated(bndpts_c_d)) return
     npts = bound_info_c%nbndptsrank
     if (npts < 1) return

     !$acc parallel loop default(present) private(i, j, k)
     do n = 1, npts
        i = bndpts_c_d(n,1)
        j = bndpts_c_d(n,2)
        k = bndpts_c_d(n,3)

        if (abs(mask_c_d(i+1,j,k)) < eps1) then
          rhs(i,j,k) = rhs(i,j,k) - 0.5 * (ekh_d(i+1,j,k) + ekh_d(i,j,k)) * (var(i+1,j,k) - var(i,j,k))*dx2i
        end if

        if (abs(mask_c_d(i-1,j,k)) < eps1) then
          rhs(i,j,k) = rhs(i,j,k) + 0.5 * (ekh_d(i,j,k) + ekh_d(i-1,j,k)) * (var(i,j,k) - var(i-1,j,k))*dx2i
        end if

        if (abs(mask_c_d(i,j+1,k)) < eps1) then
          rhs(i,j,k) = rhs(i,j,k) - 0.5 * (ekh_d(i,j+1,k) + ekh_d(i,j,k)) * (var(i,j+1,k) - var(i,j,k))*dy2i
        end if

        if (abs(mask_c_d(i,j-1,k)) < eps1) then
          rhs(i,j,k) = rhs(i,j,k) + 0.5 * (ekh_d(i,j,k) + ekh_d(i,j-1,k)) * (var(i,j,k) - var(i,j-1,k))*dy2i
        end if

        if (abs(mask_c_d(i,j,k+1)) < eps1) then
          rhs(i,j,k) = rhs(i,j,k) - 0.5 * (dzf_d(k+1)*ekh_d(i,j,k) + dzf_d(k)*ekh_d(i,j,k+1)) &
                                        * (var(i,j,k+1) - var(i,j,k))*dzh2i_d(k+1)*dzfi_d(k)
        end if

        if (abs(mask_c_d(i,j,k-1)) < eps1) then
          rhs(i,j,k) = rhs(i,j,k) + 0.5 * (dzf_d(k-1)*ekh_d(i,j,k) + dzf_d(k)*ekh_d(i,j,k-1)) &
                                        * (var(i,j,k) - var(i,j,k-1))*dzh2i_d(k)*dzfi_d(k)
        end if
     end do
     !$acc end parallel loop

   end subroutine diffc_corr_device
#endif


   subroutine initibmnorm(fname, solid_info)
     use readinput, only : read_sparse_ijk

     character(11), intent(in) :: fname

     type(solid_info_type), intent(inout) :: solid_info

     integer, allocatable :: ids_loc(:), pts_loc(:,:)

     ! Use generic read_sparse_ijk to read and distribute solid points
     call read_sparse_ijk(fname, solid_info%nsolpts, solid_info%nsolptsrank, ids_loc, pts_loc)

     ! Transfer ownership of arrays (no copying, no conversion needed)
     call move_alloc(ids_loc, solid_info%solptsrank)
     call move_alloc(pts_loc, solid_info%solpts_loc)

     !write(*,*) "rank ", myid, " has ", solid_info%nsolptsrank, " solid points from ", fname

   end subroutine initibmnorm


   subroutine initibmwallfun(fname_bnd, fname_sec, dir, bound_info)
     use modglobal, only : ifinput, ib, itot, ih, jb, jtot, jh, kb, ke, ktot, kh, &
                           xf, yf, zf, xh, yh, zh, dx, dy, dzf, xhat, yhat, zhat, eps1
     use modmpi,    only : myid, comm3d, MY_REAL, mpierr
     use initfac,   only : facnorm, facz0
     use decomp_2d, only : zstart, zend
     use readinput, only : read_sparse_ijk

     character(20), intent(in) :: fname_bnd, fname_sec
     type(bound_info_type) :: bound_info
     real, intent(in), dimension(3) :: dir
     real, dimension(ib:itot+ih) :: xgrid
     real, dimension(jb:jtot+jh) :: ygrid
     real, dimension(kb:ktot+kh) :: zgrid
     logical, dimension(bound_info%nfctsecs) :: lfctsecsrank
     logical, dimension(:), allocatable :: lbndptsrank
     real, dimension(3) :: norm, p0, p1, pxl, pxu, pyl, pyu, pzl, pzu
     integer, dimension(6) :: check
     integer, dimension(1) :: pos_min_dist
     real, dimension(6,3) :: inter
     real, dimension(6) :: inter_dists
     real :: xc, yc, zc, xl, yl, zl, xu, yu, zu
     integer n, m, norm_align, dir_align, pos, p
     character(80) :: chmess

     integer, dimension(:), allocatable :: ids_loc
     integer, dimension(:,:), allocatable :: pts_loc

     ! Read boundary points using generic read_sparse_ijk routine (skips 1 header line)
     ! Request both local and global arrays since sections need global indices
     call read_sparse_ijk(fname_bnd, bound_info%nbndpts, bound_info%nbndptsrank, ids_loc, pts_loc, nskip=1, pts_glob_out=bound_info%bndpts)

     ! Transfer ownership of local arrays (no copying, no conversion needed)
     call move_alloc(ids_loc, bound_info%bndptsrank)
     call move_alloc(pts_loc, bound_info%bndpts_loc)

     ! Build lbndptsrank lookup array for determining which sections are on this rank
     allocate(lbndptsrank(bound_info%nbndpts))
     lbndptsrank = .false.
     do m = 1, bound_info%nbndptsrank
       lbndptsrank(bound_info%bndptsrank(m)) = .true.
     end do

     allocate(bound_info%secfacids(bound_info%nfctsecs))
     allocate(bound_info%secareas(bound_info%nfctsecs))
     allocate(bound_info%secbndptids(bound_info%nfctsecs))
     !allocate(bound_info%intpts(bound_info%nfctsecs,3))
     allocate(bound_info%bnddst(bound_info%nfctsecs))
     !allocate(bound_info%bndvec(bound_info%nfctsecs,3))
     allocate(bound_info%recpts(bound_info%nfctsecs,3))
     allocate(bound_info%recids_u(bound_info%nfctsecs,3))
     allocate(bound_info%recids_v(bound_info%nfctsecs,3))
     allocate(bound_info%recids_w(bound_info%nfctsecs,3))
     allocate(bound_info%recids_c(bound_info%nfctsecs,3))
     allocate(bound_info%lcomprec(bound_info%nfctsecs))
     allocate(bound_info%lskipsec(bound_info%nfctsecs))

     dir_align = alignment(dir)
     select case(dir_align)
     case(1)
       xgrid = xh
       ygrid = yf
       zgrid = zf
     case(2)
       xgrid = xf
       ygrid = yh
       zgrid = zf
     case(3)
       xgrid = xf
       ygrid = yf
       zgrid = zh
     case(0)
       xgrid = xf
       ygrid = yf
       zgrid = zf
     end select

     if (myid == 0) then
       open (ifinput, file=fname_sec)
       read (ifinput, '(a80)') chmess
       do n = 1, bound_info%nfctsecs
         read (ifinput, *) bound_info%secfacids(n), bound_info%secareas(n), bound_info%secbndptids(n), bound_info%bnddst(n)
                           !bound_info%intpts(n,1),  bound_info%intpts(n,2), bound_info%intpts(n,3)
       end do
       close (ifinput)

       do n = 1,bound_info%nfctsecs
         m = bound_info%secbndptids(n)
         !bound_info%bndvec(n,1) = xgrid(bound_info%bndpts(m,1)) - bound_info%intpts(n,1)
         !bound_info%bndvec(n,2) = ygrid(bound_info%bndpts(m,2)) - bound_info%intpts(n,2)
         !bound_info%bndvec(n,3) = zgrid(bound_info%bndpts(m,3)) - bound_info%intpts(n,3)
         !bound_info%bnddst(n) = norm2(bound_info%bndvec(n,:))
         !write(*,*) bound_info%bnddst(n)
         !bound_info%bndvec(n,:) = bound_info%bndvec(n,:) / bound_info%bnddst(n)

         norm = facnorm(bound_info%secfacids(n),:)
         norm_align = alignment(norm)

         if ((dir_align /= 0 .and. dir_align == norm_align) .or. (facz0(bound_info%secfacids(n)) < eps1)) then
           ! (for velocities) if the facet is aligned with the grid AND in the same direction as the current velocity grid direction
           ! therefore no tangential component, don't need to calculate shear stress
           bound_info%lskipsec(n) = .true.
           cycle
         else
            bound_info%lskipsec(n) = .false.
         end if

         if (log(bound_info%bnddst(n)/facz0(bound_info%secfacids(n))) > 1. .or. lnorec) then ! the wall function is well-defined
            bound_info%lcomprec(n) = .true. ! do simple reconstruction
         else ! need to reconstruct
           bound_info%lcomprec(n) = .false.
           ! Find reconstruction point
           ! cell centre (of current grid)
           xc = xgrid(bound_info%bndpts(m,1))
           yc = ygrid(bound_info%bndpts(m,2))
           zc = zgrid(bound_info%bndpts(m,3))

           ! cell edges
           xl = xc - dx/2.
           xu = xc + dx/2.
           yl = yc - dy/2.
           yu = yc + dy/2.
           zl = zc - dzf(1)/2. ! assumes equidistant
           zu = zc + dzf(1)/2. ! assumes equidistant

           ! points on planes
           pxl = (/xl, yc, zc/)
           pxu = (/xu, yc, zc/)
           pyl = (/xc, yl, zc/)
           pyu = (/xc, yu, zc/)
           pzl = (/xc, yc, zl/)
           pzu = (/xc, yc, zu/)

           p0 = (/xc, yc, zc/)
           p1 = p0 + norm * sqrt(3.)*(dx*dy*dzf(1))**(1./3.)

           call plane_line_intersection(xhat, pxl, p0, p1, inter(1,:), check(1), inter_dists(1))
           call plane_line_intersection(xhat, pxu, p0, p1, inter(2,:), check(2), inter_dists(2))
           call plane_line_intersection(yhat, pyl, p0, p1, inter(3,:), check(3), inter_dists(3))
           call plane_line_intersection(yhat, pyu, p0, p1, inter(4,:), check(4), inter_dists(4))
           call plane_line_intersection(zhat, pzl, p0, p1, inter(5,:), check(5), inter_dists(5))
           call plane_line_intersection(zhat, pzu, p0, p1, inter(6,:), check(6), inter_dists(6))

           pos_min_dist = minloc(inter_dists, mask=check==1)
           pos = pos_min_dist(1)

           if (pos == 0) then
             write(*,*) "ERROR: no intersection found"
             stop 1
           else
             bound_info%recpts(n,:) = inter(pos,:) ! x y z
           end if

           ! find which cell the point lies in
           bound_info%recids_u(n,1) = cell_index(bound_info%recpts(n,1), xh)
           bound_info%recids_u(n,2) = cell_index(bound_info%recpts(n,2), yf)
           bound_info%recids_u(n,3) = cell_index(bound_info%recpts(n,3), zf)

           bound_info%recids_v(n,1) = cell_index(bound_info%recpts(n,1), xf)
           bound_info%recids_v(n,2) = cell_index(bound_info%recpts(n,2), yh)
           bound_info%recids_v(n,3) = cell_index(bound_info%recpts(n,3), zf)

           bound_info%recids_w(n,1) = cell_index(bound_info%recpts(n,1), xf)
           bound_info%recids_w(n,2) = cell_index(bound_info%recpts(n,2), yf)
           bound_info%recids_w(n,3) = cell_index(bound_info%recpts(n,3), zh)

           bound_info%recids_c(n,1) = cell_index(bound_info%recpts(n,1), xf)
           bound_info%recids_c(n,2) = cell_index(bound_info%recpts(n,2), yf)
           bound_info%recids_c(n,3) = cell_index(bound_info%recpts(n,3), zf)

           ! check to see if recids is inside the domain
           if (bound_info%recids_u(n,1) < ib .or. bound_info%recids_u(n,1)+1 > itot+ih .or. &
               bound_info%recids_u(n,2) < jb .or. bound_info%recids_u(n,2)+1 > jtot+jh .or. &
               bound_info%recids_u(n,3) < kb .or. bound_info%recids_u(n,3)+1 > ke+kh) then
               ! if (myid == 0) then
               !   write(*,*) "DEBUG: skipping section n=", n, "- u-velocity reconstruction index out of bounds:", bound_info%recids_u(n,1), bound_info%recids_u(n,1)+1, bound_info%recids_u(n,2), bound_info%recids_u(n,2)+1, bound_info%recids_u(n,3), bound_info%recids_u(n,3)+1
               ! end if
               bound_info%lskipsec(n) = .true.
               cycle
           end if
           if (bound_info%recids_v(n,1) < ib .or. bound_info%recids_v(n,1)+1 > itot+ih .or. &
               bound_info%recids_v(n,2) < jb .or. bound_info%recids_v(n,2)+1 > jtot+jh .or. &
               bound_info%recids_v(n,3) < kb .or. bound_info%recids_v(n,3)+1 > ke+kh) then
               ! if (myid == 0) then
               !   write(*,*) "DEBUG: skipping section n=", n, "- v-velocity reconstruction index out of bounds:", bound_info%recids_v(n,1), bound_info%recids_v(n,1)+1, bound_info%recids_v(n,2), bound_info%recids_v(n,2)+1, bound_info%recids_v(n,3), bound_info%recids_v(n,3)+1
               ! end if
               bound_info%lskipsec(n) = .true.
             cycle
           end if
           if (bound_info%recids_w(n,1) < ib .or. bound_info%recids_w(n,1)+1 > itot+ih .or. &
               bound_info%recids_w(n,2) < jb .or. bound_info%recids_w(n,2)+1 > jtot+jh .or. &
               bound_info%recids_w(n,3) < kb .or. bound_info%recids_w(n,3)+1 > ke+kh) then
               ! if (myid == 0) then
               !   write(*,*) "DEBUG: skipping section n=", n, "- w-velocity reconstruction index out of bounds:", bound_info%recids_w(n,1), bound_info%recids_w(n,1)+1, bound_info%recids_w(n,2), bound_info%recids_w(n,2)+1, bound_info%recids_w(n,3), bound_info%recids_w(n,3)+1
               ! end if
               bound_info%lskipsec(n) = .true.
               cycle
           end if
           if (bound_info%recids_c(n,1) < ib .or. bound_info%recids_c(n,1)+1 > itot+ih .or. &
               bound_info%recids_c(n,2) < jb .or. bound_info%recids_c(n,2)+1 > jtot+jh .or. &
               bound_info%recids_c(n,3) < kb .or. bound_info%recids_c(n,3)+1 > ke+kh) then
               ! if (myid == 0) then
               !   write(*,*) "DEBUG: skipping section n=", n, "- c-velocity index out of bounds:", bound_info%recids_c(n,1), bound_info%recids_c(n,1)+1, bound_info%recids_c(n,2), bound_info%recids_c(n,2)+1, bound_info%recids_c(n,3), bound_info%recids_c(n,3)+1
               ! end if
               bound_info%lskipsec(n) = .true.
             cycle
           end if

           !check recpts is inside the box defined by the corners
           ! u
           if ((bound_info%recpts(n,1) < xh(bound_info%recids_u(n,1))) .or. &
               (bound_info%recpts(n,1) > xh(bound_info%recids_u(n,1)+1))) then
             write(*,*) "ERROR: x out of bounds"
             stop 1
           end if
           if ((bound_info%recpts(n,2) < yf(bound_info%recids_u(n,2))) .or. &
               (bound_info%recpts(n,2) > yf(bound_info%recids_u(n,2)+1))) then
             write(*,*) "ERROR: y out of bounds"
             stop 1
           end if
           if ((bound_info%recpts(n,3) < zf(bound_info%recids_u(n,3))) .or. &
               (bound_info%recpts(n,3) > zf(bound_info%recids_u(n,3)+1))) then
             write(*,*) "ERROR: z out of bounds"
             stop 1
           end if

           ! v
           if ((bound_info%recpts(n,1) < xf(bound_info%recids_v(n,1))) .or. &
               (bound_info%recpts(n,1) > xf(bound_info%recids_v(n,1)+1))) then
             write(*,*) "ERROR: x out of bounds"
             stop 1
           end if
           if ((bound_info%recpts(n,2) < yh(bound_info%recids_v(n,2))) .or. &
               (bound_info%recpts(n,2) > yh(bound_info%recids_v(n,2)+1))) then
             write(*,*) "ERROR: y out of bounds"
             stop 1
           end if
           if ((bound_info%recpts(n,3) < zf(bound_info%recids_v(n,3))) .or. &
               (bound_info%recpts(n,3) > zf(bound_info%recids_v(n,3)+1))) then
             write(*,*) "ERROR: z out of bounds"
             stop 1
           end if

           ! w
           if ((bound_info%recpts(n,1) < xf(bound_info%recids_w(n,1))) .or. &
               (bound_info%recpts(n,1) > xf(bound_info%recids_w(n,1)+1))) then
             write(*,*) "ERROR: x out of bounds"
             stop 1
           end if
           if ((bound_info%recpts(n,2) < yf(bound_info%recids_w(n,2))) .or. &
               (bound_info%recpts(n,2) > yf(bound_info%recids_w(n,2)+1))) then
             write(*,*) "ERROR: y out of bounds"
             stop 1
           end if
           if ((bound_info%recpts(n,3) < zh(bound_info%recids_w(n,3))) .or. &
               (bound_info%recpts(n,3) > zh(bound_info%recids_w(n,3)+1))) then
             write(*,*) "ERROR: z out of bounds"
             stop 1
           end if
         end if
       end do
     end if ! myid==0

     call MPI_BCAST(bound_info%secfacids,   bound_info%nfctsecs,   MPI_INTEGER, 0, comm3d, mpierr)
     call MPI_BCAST(bound_info%secareas,    bound_info%nfctsecs,   MY_REAL,     0, comm3d, mpierr)
     call MPI_BCAST(bound_info%secbndptids, bound_info%nfctsecs,   MPI_INTEGER, 0, comm3d, mpierr)
     !call MPI_BCAST(bound_info%intpts,      bound_info%nfctsecs*3, MY_REAL,     0, comm3d, mpierr)
     !call MPI_BCAST(bound_info%bndvec,      bound_info%nfctsecs*3, MY_REAL,     0, comm3d, mpierr)
     call MPI_BCAST(bound_info%bnddst,      bound_info%nfctsecs,   MY_REAL,     0, comm3d, mpierr)
     call MPI_BCAST(bound_info%recpts,      bound_info%nfctsecs*3, MY_REAL,     0, comm3d, mpierr)
     call MPI_BCAST(bound_info%recids_u,    bound_info%nfctsecs*3, MPI_INTEGER, 0, comm3d, mpierr)
     call MPI_BCAST(bound_info%recids_v,    bound_info%nfctsecs*3, MPI_INTEGER, 0, comm3d, mpierr)
     call MPI_BCAST(bound_info%recids_w,    bound_info%nfctsecs*3, MPI_INTEGER, 0, comm3d, mpierr)
     call MPI_BCAST(bound_info%recids_c,    bound_info%nfctsecs*3, MPI_INTEGER, 0, comm3d, mpierr)
     call MPI_BCAST(bound_info%lskipsec,    bound_info%nfctsecs,   MPI_LOGICAL, 0, comm3d, mpierr)
     call MPI_BCAST(bound_info%lcomprec,    bound_info%nfctsecs,   MPI_LOGICAL, 0, comm3d, mpierr)

     ! Determine whether section needs to be updated by this rank
     bound_info%nfctsecsrank = 0
     do n = 1, bound_info%nfctsecs
       if (lbndptsrank(bound_info%secbndptids(n))) then
          lfctsecsrank(n) = .true.
          bound_info%nfctsecsrank = bound_info%nfctsecsrank + 1
        else
          lfctsecsrank(n) = .false.
       end if
     end do

     ! Store indices of sections on current rank - only loop through these sections
     allocate(bound_info%fctsecsrank(bound_info%nfctsecsrank))
     ! allocate local arrays
     allocate(bound_info%secfacids_loc(bound_info%nfctsecsrank))
     allocate(bound_info%secareas_loc(bound_info%nfctsecsrank))
     allocate(bound_info%secbndpts_loc(bound_info%nfctsecsrank,3))
     allocate(bound_info%bnddst_loc(bound_info%nfctsecsrank))
     allocate(bound_info%recpts_loc(bound_info%nfctsecsrank,3))
     allocate(bound_info%recids_u_loc(bound_info%nfctsecsrank,3))
     allocate(bound_info%recids_v_loc(bound_info%nfctsecsrank,3))
     allocate(bound_info%recids_w_loc(bound_info%nfctsecsrank,3))
     allocate(bound_info%recids_c_loc(bound_info%nfctsecsrank,3))
     allocate(bound_info%lcomprec_loc(bound_info%nfctsecsrank))
     allocate(bound_info%lskipsec_loc(bound_info%nfctsecsrank))

     m = 0
     p = 0  ! counter for sections overridden to simple reconstruction
     do n = 1, bound_info%nfctsecs
       if (lfctsecsrank(n)) then
          m = m + 1
          bound_info%fctsecsrank(m) = n
          bound_info%secfacids_loc(m) = bound_info%secfacids(n) ! facet id
          bound_info%secareas_loc(m) = bound_info%secareas(n)
          bound_info%secbndpts_loc(m,:) = bound_info%bndpts(bound_info%secbndptids(n),:) ! boundary point location (in global coordinates)

          if (bound_info%bndpts(bound_info%secbndptids(n),1) < zstart(1) .or. bound_info%bndpts(bound_info%secbndptids(n),1) > zend(1)) then
            write(*,*) "problem in x boundary points on : ", myid, n, bound_info%secbndptids(n), bound_info%bndpts(bound_info%secbndptids(n),1), zstart(1), zend(1)
          end if
          if (bound_info%bndpts(bound_info%secbndptids(n),2) < zstart(2) .or. bound_info%bndpts(bound_info%secbndptids(n),2) > zend(2)) then
             write(*,*) "problem in y boundary points on rank: ", myid, n, bound_info%secbndptids(n), bound_info%bndpts(bound_info%secbndptids(n),2), zstart(2), zend(2)
          end if

          bound_info%bnddst_loc(m) = bound_info%bnddst(n)
          bound_info%recpts_loc(m,:) = bound_info%recpts(n,:)
          bound_info%recids_u_loc(m,:) = bound_info%recids_u(n,:)
          bound_info%recids_v_loc(m,:) = bound_info%recids_v(n,:)
          bound_info%recids_w_loc(m,:) = bound_info%recids_w(n,:)
          bound_info%recids_c_loc(m,:) = bound_info%recids_c(n,:)
          bound_info%lcomprec_loc(m) = bound_info%lcomprec(n)

          ! If any reconstruction cell index falls outside this rank's halo-accessible range,
          ! override to simple reconstruction to avoid out-of-bounds access in trilinear_interp_var.
          !
          ! Skipped sections (lskipsec(n)) never reach trilinear_interp_var
          ! so they are harmless even if their reconstruction cell is out of range.
          !
          ! These sections then self-skip at runtime via the log(dist/z0) <= 1 check, contributing zero stress — a small, localized inaccuracy
          if (.not. bound_info%lcomprec_loc(m) .and. .not. bound_info%lskipsec(n)) then
            if (bound_info%recids_u_loc(m,1) < zstart(1)-1 .or. bound_info%recids_u_loc(m,1) > zend(1)+1 .or. &
                bound_info%recids_u_loc(m,2) < zstart(2)-1 .or. bound_info%recids_u_loc(m,2) > zend(2)+1 .or. &
                bound_info%recids_v_loc(m,1) < zstart(1)-1 .or. bound_info%recids_v_loc(m,1) > zend(1)+1 .or. &
                bound_info%recids_v_loc(m,2) < zstart(2)-1 .or. bound_info%recids_v_loc(m,2) > zend(2)+1 .or. &
                bound_info%recids_w_loc(m,1) < zstart(1)-1 .or. bound_info%recids_w_loc(m,1) > zend(1)+1 .or. &
                bound_info%recids_w_loc(m,2) < zstart(2)-1 .or. bound_info%recids_w_loc(m,2) > zend(2)+1 .or. &
                bound_info%recids_c_loc(m,1) < zstart(1)-1 .or. bound_info%recids_c_loc(m,1) > zend(1)+1 .or. &
                bound_info%recids_c_loc(m,2) < zstart(2)-1 .or. bound_info%recids_c_loc(m,2) > zend(2)+1) then
              bound_info%lcomprec_loc(m) = .true.
              p = p + 1
            end if
          end if

          bound_info%lskipsec_loc(m) = bound_info%lskipsec(n)
       end if
     end do

     if (p > 0) then
       write(*,*) "WARNING initibmwallfun: MPI rank", myid, "overrode", p, &
                  "facet section(s) to simple reconstruction because reconstruction cell falls outside halo range."
     end if

     call build_wallfun_cache(bound_info)

     deallocate(bound_info%bndpts)
     deallocate(bound_info%secfacids)
     deallocate(bound_info%secbndptids)
     deallocate(bound_info%bnddst)
     deallocate(bound_info%recpts)
     deallocate(bound_info%recids_u)
     deallocate(bound_info%recids_v)
     deallocate(bound_info%recids_w)
     deallocate(bound_info%recids_c)
     deallocate(bound_info%lcomprec)
     deallocate(bound_info%lskipsec)
     deallocate(lbndptsrank)

   end subroutine initibmwallfun


   !> Derive, once, everything the wall functions would otherwise recompute per
   !! section per Runge-Kutta stage.
   !!
   !! The facet section geometry is fixed as soon as the facet files are read,
   !! so none of the following depends on the solution:
   !!
   !!   secdist_loc  the wall distance. For a reconstructed section the offset
   !!                from the boundary point to the reconstruction point costs a
   !!                square root and three grid lookups.
   !!   bndloc_loc   the boundary point as rank-local (i,j,k), and recloc_* the
   !!                reconstruction cells likewise, so zstart is applied here
   !!                rather than per section.
   !!   recw_loc     the trilinear offsets xd, yd, zd for the u, v, w and c
   !!                stencils: three divisions and six grid lookups each.
   !!   wfc_loc      the roughness terms - log(dist/z0), log(z0/z0h),
   !!                fkar2/logdz**2, the two b1*sqdz products the unstable
   !!                branch needs, and the whole neutral coefficient
   !!                (fkar/logdz)**2. Between them these are every logarithm
   !!                and every square root the transfer coefficients evaluated.
   !!
   !! Both the host wall functions and the device kernels read this, the latter
   !! through mirror_sections, which only copies. Deriving it in one place is
   !! what keeps them bit-identical: the alternative, deriving it twice, is what
   !! the earlier device-only version did and it made every change here a choice
   !! between speed on one side and agreement between the two.
   !!
   !! The expressions are written as the wall functions wrote them, in the same
   !! association, so building them here changes no bits.
   subroutine build_wallfun_cache(bi)
     use modglobal, only : ib, ie, jb, je, xf, yf, zf, xh, yh, zh, eps1, fkar, fkar2
     use initfac,   only : facz0, facz0h
     use decomp_2d, only : zstart
     use modmpi,    only : myid
     implicit none
     type(bound_info_type), intent(inout) :: bi

     real, parameter :: b1 = 9.4
     real, parameter :: dm = 7.4
     real, parameter :: dh = 5.3

     integer :: n, sec, fac, i, j
     real    :: d, z0, z0h, logdz, sqdz

     n = bi%nfctsecsrank
     if (n < 1) return

     allocate(bi%secdist_loc(n))
     allocate(bi%bndloc_loc(n,3))
     allocate(bi%recloc_u(n,3), bi%recloc_v(n,3), bi%recloc_w(n,3), bi%recloc_c(n,3))
     allocate(bi%recw_loc(n,RECW_N))
     allocate(bi%wfc_loc(n,WFC_N))

     bi%recw_loc = 0.
     bi%wfc_loc  = 0.

     do sec = 1, n
       fac = bi%secfacids_loc(sec)

       bi%bndloc_loc(sec,1) = bi%secbndpts_loc(sec,1) - zstart(1) + 1
       bi%bndloc_loc(sec,2) = bi%secbndpts_loc(sec,2) - zstart(2) + 1
       bi%bndloc_loc(sec,3) = bi%secbndpts_loc(sec,3) - zstart(3) + 1

       bi%recloc_u(sec,:) = bi%recids_u_loc(sec,:) - zstart + 1
       bi%recloc_v(sec,:) = bi%recids_v_loc(sec,:) - zstart + 1
       bi%recloc_w(sec,:) = bi%recids_w_loc(sec,:) - zstart + 1
       bi%recloc_c(sec,:) = bi%recids_c_loc(sec,:) - zstart + 1

       ! The wall functions used to test this on every section on every stage,
       ! and a kernel cannot stop. Both conditions are on fixed indices, so they
       ! are settled here instead.
       i = bi%bndloc_loc(sec,1)
       j = bi%bndloc_loc(sec,2)
       if ((i < ib) .or. (i > ie) .or. (j < jb) .or. (j > je)) then
         write(*,*) 'problem in wallfunmom ', myid, bi%secbndpts_loc(sec,1), bi%secbndpts_loc(sec,2)
         stop 1
       end if

       if (bi%lcomprec_loc(sec) .or. lnorec) then
         ! The reconstruction point is not used, and recids_* may never have
         ! been set, so the offsets are left at zero.
         d = bi%bnddst_loc(sec)
       else
         call validate_reconstruction_cell('sec', bi%recids_u_loc(sec,:))
         call validate_reconstruction_cell('sec', bi%recids_v_loc(sec,:))
         call validate_reconstruction_cell('sec', bi%recids_w_loc(sec,:))
         call validate_reconstruction_cell('sec', bi%recids_c_loc(sec,:))

         d = bi%bnddst_loc(sec) + norm2((/bi%recpts_loc(sec,1) - xf(bi%secbndpts_loc(sec,1)), &
                                          bi%recpts_loc(sec,2) - yf(bi%secbndpts_loc(sec,2)), &
                                          bi%recpts_loc(sec,3) - zf(bi%secbndpts_loc(sec,3))/))

         call cell_offsets(bi%recids_u_loc(sec,:), xh, yf, zf, bi%recpts_loc(sec,:), &
                           bi%recw_loc(sec,RECW_U:RECW_U+2))
         call cell_offsets(bi%recids_v_loc(sec,:), xf, yh, zf, bi%recpts_loc(sec,:), &
                           bi%recw_loc(sec,RECW_V:RECW_V+2))
         call cell_offsets(bi%recids_w_loc(sec,:), xf, yf, zh, bi%recpts_loc(sec,:), &
                           bi%recw_loc(sec,RECW_W:RECW_W+2))
         call cell_offsets(bi%recids_c_loc(sec,:), xf, yf, zf, bi%recpts_loc(sec,:), &
                           bi%recw_loc(sec,RECW_C:RECW_C+2))
       end if
       bi%secdist_loc(sec) = d

       ! Skipped sections never reach the wall function, and lskipsec already
       ! covers facz0 < eps1, so the logarithm below is always well defined.
       if (bi%lskipsec_loc(sec)) cycle
       z0 = facz0(fac)
       if (z0 < eps1) cycle

       logdz = log(d/z0)
       bi%wfc_loc(sec,WFC_LOGDZ) = logdz
       ! The wall functions test logdz > 1. before touching anything else, so
       ! for a section that fails it the rest of the row stays zero.
       if (logdz <= 1.) cycle

       ! Only the stability paths read logzh, and a zero heat roughness makes it
       ! non-finite there in any case; left at zero rather than propagating an
       ! infinity out of initialisation.
       z0h = facz0h(fac)
       if (z0h > 0.) bi%wfc_loc(sec,WFC_LOGZH) = log(z0/z0h)

       sqdz = sqrt(d/z0)
       bi%wfc_loc(sec,WFC_COEF)   = fkar2/(logdz**2)
       bi%wfc_loc(sec,WFC_CM)     = (dm*fkar2)/(logdz**2)*b1*sqdz
       bi%wfc_loc(sec,WFC_CH)     = (dh*fkar2)/(logdz**2)*b1*sqdz
       ! mom_transfer_coef_neutral in full. Its own column rather than a reuse
       ! of WFC_COEF: the two are the same number mathematically but not
       ! bitwise, one squaring a rounded quotient and the other dividing by a
       ! rounded square.
       bi%wfc_loc(sec,WFC_CTMNEU) = (fkar/logdz)**2
     end do

   end subroutine build_wallfun_cache


   !> Normalised offsets of a point within its containing grid cell.
   !!
   !! cell holds global indices and the grid arrays are indexed by them; the
   !! expression is the one trilinear_interp used to evaluate per call.
   subroutine cell_offsets(cell, xgrid, ygrid, zgrid, p, d)
     use modglobal, only : ib, jb, kb, itot, jtot, ktot, ih, jh, kh
     implicit none
     integer, intent(in)  :: cell(3)
     real,    intent(in)  :: xgrid(ib:itot+ih), ygrid(jb:jtot+jh), zgrid(kb:ktot+kh)
     real,    intent(in)  :: p(3)
     real,    intent(out) :: d(3)

     d(1) = (p(1) - xgrid(cell(1))) / (xgrid(cell(1)+1) - xgrid(cell(1)))
     d(2) = (p(2) - ygrid(cell(2))) / (ygrid(cell(2)+1) - ygrid(cell(2)))
     d(3) = (p(3) - zgrid(cell(3))) / (zgrid(cell(3)+1) - zgrid(cell(3)))

   end subroutine cell_offsets


   subroutine validate_reconstruction_cell(tag, cell)
     use modglobal, only : ib, ie, jb, je
     use decomp_2d, only : zstart
     implicit none
     character(len=*), intent(in) :: tag
     integer, intent(in) :: cell(3)
     integer :: i, j, k

     i = cell(1) - zstart(1) + 1
     j = cell(2) - zstart(2) + 1
     k = cell(3) - zstart(3) + 1
     if ((i < ib-1) .or. (i > ie+1) .or. (j < jb-1) .or. (j > je+1)) then
       write(*,*) 'problem in trilinear_interp_var ', tag, i, j, k
       stop 1
     end if

   end subroutine validate_reconstruction_cell


   !> Compare the per-section cache against the expressions it replaced.
   !!
   !! Everything build_wallfun_cache precomputes was, until it was hoisted out
   !! of the time loop, evaluated inside the wall functions on every Runge-Kutta
   !! stage. The expressions below are those, written out again here rather than
   !! shared with the producer, so that the two can disagree: a swapped stencil
   !! grid, a mislabelled cache column or an index left global then shows up as
   !! a mismatch instead of as a quietly wrong wall stress. They are compared
   !! exactly, because both sides evaluate the same expression on the same host.
   !!
   !! problem comes back empty when everything agrees, and otherwise names the
   !! section set and the quantity that did not.
   subroutine check_wallfun_cache(problem)
     implicit none
     character(len=*), intent(out) :: problem

     problem = ''
     call check_one('u', bound_info_u)
     call check_one('v', bound_info_v)
     call check_one('w', bound_info_w)
     call check_one('c', bound_info_c)

   contains

     subroutine check_one(tag, bi)
       use modglobal, only : xf, yf, zf, xh, yh, zh, eps1, fkar, fkar2
       use initfac,   only : facz0, facz0h
       use decomp_2d, only : zstart
       implicit none
       character(len=*), intent(in) :: tag
       type(bound_info_type), intent(in) :: bi

       real, parameter :: b1 = 9.4
       real, parameter :: dm = 7.4
       real, parameter :: dh = 5.3

       integer :: n, sec, fac
       real    :: d, z0, z0h, logdz, sqdz, want(WFC_N), wantw(RECW_N)

       if (problem /= '') return
       if (.not. allocated(bi%wfc_loc)) return
       n = bi%nfctsecsrank
       if (n < 1) return

       do sec = 1, n
         fac = bi%secfacids_loc(sec)

         if (any(bi%bndloc_loc(sec,:) /= bi%secbndpts_loc(sec,:) - zstart + 1)) &
           call flag(tag, 'boundary point')
         if (any(bi%recloc_u(sec,:) /= bi%recids_u_loc(sec,:) - zstart + 1)) call flag(tag, 'u reconstruction cell')
         if (any(bi%recloc_v(sec,:) /= bi%recids_v_loc(sec,:) - zstart + 1)) call flag(tag, 'v reconstruction cell')
         if (any(bi%recloc_w(sec,:) /= bi%recids_w_loc(sec,:) - zstart + 1)) call flag(tag, 'w reconstruction cell')
         if (any(bi%recloc_c(sec,:) /= bi%recids_c_loc(sec,:) - zstart + 1)) call flag(tag, 'c reconstruction cell')

         wantw = 0.
         if (bi%lcomprec_loc(sec) .or. lnorec) then
           d = bi%bnddst_loc(sec)
         else
           d = bi%bnddst_loc(sec) + norm2((/bi%recpts_loc(sec,1) - xf(bi%secbndpts_loc(sec,1)), &
                                            bi%recpts_loc(sec,2) - yf(bi%secbndpts_loc(sec,2)), &
                                            bi%recpts_loc(sec,3) - zf(bi%secbndpts_loc(sec,3))/))
           ! The grid each stencil is staggered on, restated.
           wantw(RECW_U:RECW_U+2) = offsets(bi%recids_u_loc(sec,:), xh, yf, zf, bi%recpts_loc(sec,:))
           wantw(RECW_V:RECW_V+2) = offsets(bi%recids_v_loc(sec,:), xf, yh, zf, bi%recpts_loc(sec,:))
           wantw(RECW_W:RECW_W+2) = offsets(bi%recids_w_loc(sec,:), xf, yf, zh, bi%recpts_loc(sec,:))
           wantw(RECW_C:RECW_C+2) = offsets(bi%recids_c_loc(sec,:), xf, yf, zf, bi%recpts_loc(sec,:))
         end if
         if (bi%secdist_loc(sec) /= d) call flag(tag, 'wall distance')
         if (any(bi%recw_loc(sec,:) /= wantw)) call flag(tag, 'trilinear offsets')

         want = 0.
         z0 = facz0(fac)
         if (.not. bi%lskipsec_loc(sec) .and. z0 >= eps1) then
           logdz = log(d/z0)
           want(WFC_LOGDZ) = logdz
           if (logdz > 1.) then
             z0h = facz0h(fac)
             if (z0h > 0.) want(WFC_LOGZH) = log(z0/z0h)
             sqdz = sqrt(d/z0)
             want(WFC_COEF)   = fkar2/(logdz**2)
             want(WFC_CM)     = (dm*fkar2)/(logdz**2)*b1*sqdz
             want(WFC_CH)     = (dh*fkar2)/(logdz**2)*b1*sqdz
             want(WFC_CTMNEU) = (fkar/logdz)**2
           end if
         end if
         if (any(bi%wfc_loc(sec,:) /= want)) call flag(tag, 'roughness terms')

         if (problem /= '') exit
       end do

     end subroutine check_one

     !> Normalised offsets within the containing cell, as trilinear_interp had them.
     function offsets(cell, xgrid, ygrid, zgrid, p)
       use modglobal, only : ib, jb, kb, itot, jtot, ktot, ih, jh, kh
       implicit none
       integer, intent(in) :: cell(3)
       real,    intent(in) :: xgrid(ib:itot+ih), ygrid(jb:jtot+jh), zgrid(kb:ktot+kh)
       real,    intent(in) :: p(3)
       real :: offsets(3)

       offsets(1) = (p(1) - xgrid(cell(1))) / (xgrid(cell(1)+1) - xgrid(cell(1)))
       offsets(2) = (p(2) - ygrid(cell(2))) / (ygrid(cell(2)+1) - ygrid(cell(2)))
       offsets(3) = (p(3) - zgrid(cell(3))) / (zgrid(cell(3)+1) - zgrid(cell(3)))
     end function offsets

     subroutine flag(tag, what)
       implicit none
       character(len=*), intent(in) :: tag, what

       if (problem == '') problem = 'wall function cache ('//tag//'): '//what
     end subroutine flag

   end subroutine check_wallfun_cache



   subroutine plane_line_intersection(norm, V0, P0, P1, I, check, dist)
     use modglobal, only : vec0, eps1
     implicit none
     ! determines the intersection of a plane and a line segment
     ! norm: plane normal
     ! V0: point on the plane
     ! P0: start of line segment
     ! P1: end of line segment
     ! I: intersection point
     ! dist: distance from P0 to intersection point
     ! check: 0 if no intersection
     !        1 if unique intersection
     !        2 if line segment is in the plane
     !        3 if intersection is outside line segment
     real, intent(in),  dimension(3) :: norm, V0, P0, P1
     real, intent(out), dimension(3) :: I
     integer, intent(out) :: check
     real, intent(out) :: dist
     real, dimension(3) :: u, w
     real :: D, N, sI

     I = vec0
     w = P0 - V0
     u = P1 - P0
     D = dot_product(norm, u)
     N =-dot_product(norm, w)

     if (abs(D) < eps1) then ! line orthogonal to plane normal -> segment parallel to plane
       if (abs(N) < eps1) then ! start point is on the plane -> segment lies in the plane
         check = 2
         return
       else
         check = 0
         return
       end if
     end if

     sI = N / D
     I = P0 + sI * u
     dist = norm2(I - P0)

     if ((sI < 0.) .or. (sI > 1.)) then
       check = 3
     else
       check = 1
     end if

   end subroutine plane_line_intersection


   subroutine ibmnorm
     use modglobal,   only : ih, jh, kh, ihc, jhc, khc, nsv, dzf, zh, kb, ke, kh, nsv, libm, ltempeq, lmoist, iadv_sv, iadv_cd2, iadv_thl, lconservativeibm
     use modfields,   only : um, vm, wm, thlm, qtm, svm, up, vp, wp, thlp, qtp, svp, thl0, qt0, sv0, thl0av
     use modboundary, only : halos

     integer n

     if (.not. libm) return

     ! Set internal velocities to zero
     call solid(solid_info_u, um, up, 0., ih, jh, kh)
     call solid(solid_info_v, vm, vp, 0., ih, jh, kh)
     call solid(solid_info_w, wm, wp, 0., ih, jh, kh)

     ! Scalars
     ! Solid value does not matter when using second order scheme
     ! Set interior to a constant and boundary to average of fluid neighbours
     if (ltempeq) then
        call solid(solid_info_c, thlm, thlp, sum(thl0av(kb:ke)*dzf(kb:ke))/zh(ke+1), ih, jh, kh, mask_c)
        if (iadv_thl == iadv_cd2) then
          if (lconservativeibm) then
            call advecc2nd_corr_conservative(thl0, thlp)
          else
            call advecc2nd_corr_liberal(thl0, thlp)
          end if
        end if
     end if

     if (lmoist) then
       call solid(solid_info_c, qtm, qtp, 0., ih, jh, kh, mask_c)
       if (lconservativeibm) then
         call advecc2nd_corr_conservative(qt0, qtp)
       else
         call advecc2nd_corr_liberal(qt0, qtp)
       end if
     end if

     do n=1,nsv
        call solid(solid_info_c, svm(:,:,:,n), svp(:,:,:,n), 0., ihc, jhc, khc, mask_c)
        if (iadv_sv(n) == iadv_cd2) then
          if (lconservativeibm) then
            call advecc2nd_corr_conservative(sv0(:,:,:,n), svp(:,:,:,n))
          else
            call advecc2nd_corr_liberal(sv0(:,:,:,n), svp(:,:,:,n))
          end if
        end if
     end do

   end subroutine ibmnorm


   subroutine solid(solid_info, var, rhs, val, hi, hj, hk, mask)
     use modglobal, only : ib, ie, jb, je, kb, ke, ih, jh, kh, eps1

     type(solid_info_type), intent(in) :: solid_info
     integer, intent(in) :: hi, hj, hk
     real, intent(inout) :: var(ib-hi:ie+hi,jb-hj:je+hj,kb-hk:ke+hk)
     real, intent(inout) :: rhs(ib-hi:ie+hi,jb-hj:je+hj,kb   :ke+hk)
     real, intent(in) :: val
     real, intent(in), optional :: mask(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
     real :: count
     integer :: i, j, k, n

     if (present(mask) .eqv. .false.) then
        do n=1,solid_info%nsolptsrank
           !n = solid_info%solptsrank(m)
           i = solid_info%solpts_loc(n,1)
           j = solid_info%solpts_loc(n,2)
           k = solid_info%solpts_loc(n,3)
           var(i,j,k) = val
           rhs(i,j,k) = 0.
        end do

     else
        do n=1,solid_info%nsolptsrank
           !n = solid_info%solptsrank(m)
           i = solid_info%solpts_loc(n,1)
           j = solid_info%solpts_loc(n,2)
           k = solid_info%solpts_loc(n,3)
           var(i,j,k) = val
           rhs(i,j,k) = 0.
           count = 0

           ! Attempt to set zero flux BC
           if (abs(mask(i,j+1,k) - 1.) < eps1) then ! fluid neighbour
             count = count + 1
             var(i,j,k) = var(i,j,k) + var(i,j+1,k)
             rhs(i,j,k) = rhs(i,j,k) + rhs(i,j+1,k)
          end if

          if (abs(mask(i,j-1,k) - 1.) < eps1) then
             count = count + 1
             var(i,j,k) = var(i,j,k) + var(i,j-1,k)
             rhs(i,j,k) = rhs(i,j,k) + rhs(i,j-1,k)
          end if

          if (abs(mask(i,j,k+1) - 1.) < eps1) then
             count = count + 1
             var(i,j,k) = var(i,j,k) + var(i,j,k+1)
             rhs(i,j,k) = rhs(i,j,k) + rhs(i,j,k+1)
          end if

          if (abs(mask(i,j,k-1) - 1.) < eps1) then
             count = count + 1
             var(i,j,k) = var(i,j,k) + var(i,j,k-1)
             rhs(i,j,k) = rhs(i,j,k) + rhs(i,j,k-1)
          end if

          if (abs(mask(i+1,j,k) - 1.) < eps1) then
             count = count + 1
             var(i,j,k) = var(i,j,k) + var(i+1,j,k)
             rhs(i,j,k) = rhs(i,j,k) + rhs(i+1,j,k)
          end if

          if (abs(mask(i-1,j,k) - 1.) < eps1) then
             count = count + 1
             var(i,j,k) = var(i,j,k) + var(i-1,j,k)
             rhs(i,j,k) = rhs(i,j,k) + rhs(i-1,j,k)
          end if

          if (count > 0) then
             var(i,j,k) = (var(i,j,k) - val) / count
             rhs(i,j,k) = rhs(i,j,k) / count
          end if

       end do

   end if

   end subroutine solid


   ! subroutine solid_boundary(bound_info, mask, var, rhs, hi, hj, hk)
   !   use modglobal, only : eps1, ib, ie, ih, jb, je, jh, kb, ke, kh
   !   use decomp_2d, only : zstart
   !   ! uDALES 1 approach
   !   ! Not truly conservative
   !   type(bound_info_type), intent(in) :: bound_info
   !   integer, intent(in) :: hi, hj, hk
   !   real, intent(in)    :: mask(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
   !   real, intent(inout) :: var(ib-hi:ie+hi,jb-hj:je+hj,kb-hk:ke+hk)
   !   real, intent(inout) :: rhs(ib-hi:ie+hi,jb-hj:je+hj,kb   :ke+hk)
   !
   !   integer i, j, k, n, m
   !
   !   do m = 1,bound_info%nbndptsrank
   !     n = bound_info%bndptsrank(m)
   !     i = bound_info%bndpts(n,1) - zstart(1) + 1
   !     j = bound_info%bndpts(n,2) - zstart(2) + 1
   !     k = bound_info%bndpts(n,3) - zstart(3) + 1
   !
   !     if (abs(mask(i,j+1,k)) < eps1) then
   !       ! rhs(i,j+1,k) = 0.
   !       rhs(i,j+1,k) = rhs(i,j,k)
   !       var(i,j+1,k) = var(i,j,k)
   !     end if
   !
   !     if (abs(mask(i,j-1,k)) < eps1) then
   !       ! rhs(i,j-1,k) = 0.
   !       rhs(i,j-1,k) = rhs(i,j,k)
   !       var(i,j-1,k) = var(i,j,k)
   !     end if
   !
   !     if (abs(mask(i,j,k+1)) < eps1) then
   !       ! rhs(i,j,k+1) = 0.
   !       rhs(i,j,k+1) = rhs(i,j,k)
   !       var(i,j,k+1) = var(i,j,k)
   !     end if
   !
   !     if (abs(mask(i,j,k-1)) < eps1) then
   !       ! rhs(i,j,k-1) = 0.
   !       rhs(i,j,k-1) = rhs(i,j,k)
   !       var(i,j,k-1) = var(i,j,k)
   !     end if
   !
   !     if (abs(mask(i+1,j,k)) < eps1) then
   !       ! rhs(i+1,j,k) = 0.
   !       rhs(i+1,j,k) = rhs(i,j,k)
   !       var(i+1,j,k) = var(i,j,k)
   !     end if
   !
   !     if (abs(mask(i-1,j,k)) < eps1) then
   !       ! rhs(i-1,j,k) = 0.
   !       rhs(i-1,j,k) = rhs(i,j,k)
   !       var(i-1,j,k) = var(i,j,k)
   !     end if
   !
   !   end do
   !
   ! end subroutine


   subroutine advecc2nd_corr_conservative(var, rhs)
     ! Removes the advection contribution from solid velocities, which should be
     ! close to zero but are not necessarily due to pressure correction.
     ! Has a fairly drastic effect on the initial flow, but the scalar is
     ! conserved throughout the simulation.
     use modglobal,      only : eps1, ib, ie, ih, jb, je, jh, kb, ke, kh, &
                                dxi5, dyi5, dzf, dzhi, dzfi5
     use modfields,      only : u0, v0, w0
     real, intent(in)    :: var(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
     real, intent(inout) :: rhs(ib-ih:ie+ih,jb-jh:je+jh,kb   :ke+kh)
     integer :: i, j, k, n

     do n = 1,bound_info_c%nbndptsrank
      !n = bound_info_c%bndptsrank(m)
         i = bound_info_c%bndpts_loc(n,1)
         j = bound_info_c%bndpts_loc(n,2)
         k = bound_info_c%bndpts_loc(n,3)

         if ((abs(mask_u(i+1,j,k)) < eps1) .or. (abs(mask_c(i+1,j,k)) < eps1)) then
           rhs(i,j,k) = rhs(i,j,k) + u0(i+1,j,k)*(var(i+1,j,k) + var(i,j,k))*dxi5
         end if

         if ((abs(mask_u(i,j,k)) < eps1) .or. (abs(mask_c(i-1,j,k)) < eps1)) then
           rhs(i,j,k) = rhs(i,j,k) - u0(i,j,k)*(var(i-1,j,k) + var(i,j,k))*dxi5
         end if

         if ((abs(mask_v(i,j+1,k)) < eps1) .or. (abs(mask_c(i,j+1,k)) < eps1)) then
           rhs(i,j,k) = rhs(i,j,k) + v0(i,j+1,k)*(var(i,j+1,k) + var(i,j,k))*dyi5
         end if

         if ((abs(mask_v(i,j,k)) < eps1) .or. (abs(mask_c(i,j-1,k)) < eps1)) then
           rhs(i,j,k) = rhs(i,j,k) - v0(i,j,k)*(var(i,j-1,k) + var(i,j,k))*dyi5
         end if

         if ((abs(mask_w(i,j,k+1)) < eps1) .or. (abs(mask_c(i,j,k+1)) < eps1)) then
           rhs(i,j,k) = rhs(i,j,k) + w0(i,j,k+1)*(var(i,j,k+1)*dzf(k) + var(i,j,k)*dzf(k+1))*dzhi(k+1)*dzfi5(k)
         end if

         if ((abs(mask_w(i,j,k)) < eps1) .or. (abs(mask_c(i,j,k-1)) < eps1)) then
           rhs(i,j,k) = rhs(i,j,k) - w0(i,j,k)*(var(i,j,k-1)*dzf(k) + var(i,j,k)*dzf(k-1))*dzhi(k)*dzfi5(k)
         end if

     end do

   end subroutine advecc2nd_corr_conservative


   subroutine advecc2nd_corr_liberal(var, rhs)
     ! Removes the advection contribution from solid scalar points as calculated
     ! by the 2nd order scheme, and replaces it with a contribution in which the
     ! value inside the solid is equal to the value outside, thereby modelling
     ! a zero (advective) flux condition.
     ! Due to potentially nonzero solid velocities due to the pressure correction,
     ! the IBM will not be conservative.
     use modglobal,      only : eps1, ib, ie, ih, jb, je, jh, kb, ke, kh, &
                                dxi5, dyi5, dzf, dzhi, dzfi5
     use modfields,      only : u0, v0, w0
     real, intent(in)    :: var(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
     real, intent(inout) :: rhs(ib-ih:ie+ih,jb-jh:je+jh,kb   :ke+kh)
     integer :: i, j, k, n

     do n = 1,bound_info_c%nbndptsrank
      !n = bound_info_c%bndptsrank(m)
         i = bound_info_c%bndpts_loc(n,1)
         j = bound_info_c%bndpts_loc(n,2)
         k = bound_info_c%bndpts_loc(n,3)

         if (abs(mask_c(i+1,j,k)) < eps1) then ! var(i+1) is solid
           rhs(i,j,k) = rhs(i,j,k) + u0(i+1,j,k)*(var(i+1,j,k) + var(i,j,k))*dxi5 & ! negate contribution added in advection using var(i+1)
                                   - u0(i+1,j,k)*(var(i  ,j,k) + var(i,j,k))*dxi5   ! add corresponding contribution with var(i+1) = var(i)
         end if

         if (abs(mask_c(i-1,j,k)) < eps1) then ! var(i-1) is solid
           rhs(i,j,k) = rhs(i,j,k) - u0(i,j,k)*(var(i-1,j,k) + var(i,j,k))*dxi5 & ! negate contribution added in advection using var(i-1)
                                   + u0(i,j,k)*(var(i  ,j,k) + var(i,j,k))*dxi5   ! add corresponding contribution with var(i-1) = var(i)
         end if

         if (abs(mask_c(i,j+1,k)) < eps1) then ! var(j+1) is solid
           rhs(i,j,k) = rhs(i,j,k) + v0(i,j+1,k)*(var(i,j+1,k) + var(i,j,k))*dyi5 & ! negate contribution added in advection using var(j+1)
                                   - v0(i,j+1,k)*(var(i,j  ,k) + var(i,j,k))*dyi5   ! add corresponding contribution with var(j+1) = var(j)
         end if

         if (abs(mask_c(i,j-1,k)) < eps1) then ! var(j-1) is solid
           rhs(i,j,k) = rhs(i,j,k) - v0(i,j,k)*(var(i,j-1,k) + var(i,j,k))*dyi5 & ! negate contribution added in advection using var(j-1)
                                   + v0(i,j,k)*(var(i,j  ,k) + var(i,j,k))*dyi5   ! add corresponding contribution with var(j-1) = var(j)
         end if

         if (abs(mask_c(i,j,k+1)) < eps1) then ! var(k+1) is solid
           rhs(i,j,k) = rhs(i,j,k) + w0(i,j,k+1)*(var(i,j,k+1)*dzf(k) + var(i,j,k)*dzf(k+1))*dzhi(k+1)*dzfi5(k) & ! negate contribution added in advection using var(k+1)
                                   - w0(i,j,k+1)*(var(i,j,k  )*dzf(k) + var(i,j,k)*dzf(k+1))*dzhi(k+1)*dzfi5(k)   ! add corresponding contribution with var(k+1) = var(k)
         end if

         if (abs(mask_c(i,j,k-1)) < eps1) then ! var(k-1) is solid
           rhs(i,j,k) = rhs(i,j,k) - w0(i,j,k)*(var(i,j,k-1)*dzf(k) + var(i,j,k)*dzf(k-1))*dzhi(k)*dzfi5(k) & ! negate contribution added in advection using var(k-1)
                                   + w0(i,j,k)*(var(i,j,k  )*dzf(k) + var(i,j,k)*dzf(k-1))*dzhi(k)*dzfi5(k)   ! add corresponding contribution with var(k-1) = var(k)
         end if

     end do
   end subroutine advecc2nd_corr_liberal


#if !defined(_GPU) || defined(UDALES_DEBUG)
! The host wall functions below, and the helpers only they use, are the CPU
! implementation that the device port mirrors. A GPU release build has no use
! for them, so they are compiled out of it. A GPU debug build keeps them,
! because tests_cuda.f90 verifies the port by running both on the same state.
   subroutine diffu_corr
     ! Negate subgrid rhs contributions from solid points (added by diffu in modsubgrid)
     use modglobal,      only : eps1, &
                                dy2i, dzf, dzfi, dzhi, dzhiq
     use modfields,      only : u0, up
     use modsubgriddata, only : ekm
     real :: empo, emmo, emop, emom
     integer :: i, j, k, n

     do n = 1,bound_info_u%nbndptsrank
      !n = bound_info_u%bndptsrank(m)
         i = bound_info_u%bndpts_loc(n,1)
         j = bound_info_u%bndpts_loc(n,2)
         k = bound_info_u%bndpts_loc(n,3)

         if (abs(mask_u(i,j+1,k)) < eps1) then
           empo = 0.25 * ((ekm(i,j,k) + ekm(i,j+1,k)) + (ekm(i-1,j,k) + ekm(i-1,j+1,k)))
           up(i,j,k) = up(i,j,k) - empo * (u0(i,j+1,k) - u0(i,j,k))*dy2i
         end if

         if (abs(mask_u(i,j-1,k)) < eps1) then
           emmo = 0.25 * ((ekm(i,j,k) + ekm(i,j-1,k)) + (ekm(i-1,j-1,k) + ekm(i-1,j,k)))
           up(i,j,k) = up(i,j,k) + emmo * (u0(i,j,k) - u0(i,j-1,k))*dy2i
         end if

         if (abs(mask_u(i,j,k+1)) < eps1) then
           emop = (dzf(k+1) * ( ekm(i,j,k)   + ekm(i-1,j,k  ))  + &
                   dzf(k)   * ( ekm(i,j,k+1) + ekm(i-1,j,k+1))) * dzhiq(k+1)
           up(i,j,k) = up(i,j,k) - emop * (u0(i,j,k+1) - u0(i,j,k))*dzhi(k+1)*dzfi(k)
         end if

         if (abs(mask_u(i,j,k-1)) < eps1) then
           emom = (dzf(k-1) * (ekm(i,j,k  ) + ekm(i-1,j,k  ))  + &
                   dzf(k)   * (ekm(i,j,k-1) + ekm(i-1,j,k-1))) * dzhiq(k)
           up(i,j,k) = up(i,j,k) + emom * (u0(i,j,k) - u0(i,j,k-1))*dzhi(k)*dzfi(k)
         end if

     end do


   end subroutine diffu_corr


   subroutine diffv_corr
     ! Negate subgrid rhs contributions from solid points (added by diffv in modsubgrid)
     use modglobal,      only : eps1, &
                                dx2i, dzf, dzfi, dzhi, dzhiq
     use modfields,      only : v0, vp
     use modsubgriddata, only : ekm
     real :: epmo, emmo, eomp, eomm
     integer :: i, j, k, n

     do n = 1,bound_info_v%nbndptsrank
      !n = bound_info_v%bndptsrank(m)
         i = bound_info_v%bndpts_loc(n,1)
         j = bound_info_v%bndpts_loc(n,2)
         k = bound_info_v%bndpts_loc(n,3)

         if (abs(mask_v(i+1,j,k)) < eps1) then
           epmo = 0.25 * (ekm(i,j,k) + ekm(i,j-1,k) + ekm(i+1,j-1,k) + ekm(i+1,j,k))
           vp(i,j,k) = vp(i,j,k) - epmo * (v0(i+1,j,k) - v0(i,j,k))*dx2i
         end if

         if (abs(mask_v(i-1,j,k)) < eps1) then
           emmo = 0.25 * (ekm(i,j,k) + ekm(i,j-1,k) + ekm(i-1,j-1,k) + ekm(i-1,j,k))
           vp(i,j,k) = vp(i,j,k) + emmo * (v0(i,j,k) - v0(i-1,j,k))*dx2i
         end if

         if (abs(mask_v(i,j,k+1)) < eps1) then
           eomp = ( dzf(k+1) * ( ekm(i,j,k)   + ekm(i,j-1,k)  )  + &
                    dzf(k  ) * ( ekm(i,j,k+1) + ekm(i,j-1,k+1))) * dzhiq(k+1)
           vp(i,j,k) = vp(i,j,k) - eomp * (v0(i,j,k+1) - v0(i,j,k))*dzhi(k+1)*dzfi(k)
         end if

         if (abs(mask_v(i,j,k-1)) < eps1) then
           eomm = ( dzf(k-1) * ( ekm(i,j,k  )  + ekm(i,j-1,k)   ) + &
                    dzf(k)   * ( ekm(i,j,k-1)  + ekm(i,j-1,k-1))) * dzhiq(k)
           vp(i,j,k) = vp(i,j,k) + eomm * (v0(i,j,k) - v0(i,j,k-1))*dzhi(k)*dzfi(k)
         end if

     end do

   end subroutine diffv_corr


   subroutine diffw_corr
     ! Negate subgrid rhs contributions from solid points (added by diffw in modsubgrid)
     use modglobal,      only : eps1, &
                                dx2i, dy2i, dzf, dzhiq
     use modfields,      only : w0, wp
     use modsubgriddata, only : ekm
     real :: epom, emom, eopm, eomm
     integer :: i, j, k, n

     do n = 1,bound_info_w%nbndptsrank
      !n = bound_info_w%bndptsrank(m)
         i = bound_info_w%bndpts_loc(n,1)
         j = bound_info_w%bndpts_loc(n,2)
         k = bound_info_w%bndpts_loc(n,3)

         ! Account for solid w points
         if (abs(mask_w(i+1,j,k)) < eps1) then
           epom = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i+1,j,k  ))    + &
                    dzf(k  ) * ( ekm(i,j,k-1) + ekm(i+1,j,k-1))) * dzhiq(k)
           wp(i,j,k) = wp(i,j,k) - epom * (w0(i+1,j,k) - w0(i,j,k))*dx2i
         end if

         if (abs(mask_w(i-1,j,k)) < eps1) then
           emom = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i-1,j,k  ))  + &
                    dzf(k  ) * ( ekm(i,j,k-1) + ekm(i-1,j,k-1))) * dzhiq(k)
           wp(i,j,k) = wp(i,j,k) + emom * (w0(i,j,k) - w0(i-1,j,k))*dx2i
         end if

         if (abs(mask_w(i,j+1,k)) < eps1) then
           eopm = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i,j+1,k  ))  + &
                    dzf(k  ) * ( ekm(i,j,k-1) + ekm(i,j+1,k-1))) * dzhiq(k)
           wp(i,j,k) = wp(i,j,k) - eopm * (w0(i,j+1,k) - w0(i,j,k))*dy2i
         end if

         if (abs(mask_w(i,j-1,k)) < eps1) then
           eomm = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i,j-1,k  ))  + &
                    dzf(k  ) * ( ekm(i,j,k-1) + ekm(i,j-1,k-1))) * dzhiq(k)
           wp(i,j,k) = wp(i,j,k) + eomm * (w0(i,j,k) - w0(i,j-1,k))*dy2i
         end if

     end do

   end subroutine diffw_corr


   subroutine diffc_corr(var, rhs, hi, hj, hk)
     ! Negate subgrid rhs contributions from solid points (added by diffc in modsubgrid)
     use modglobal,      only : eps1, ib, ie, jb, je, kb, ke, &
                                dx2i, dy2i, dzf, dzh2i, dzfi
     use modsubgriddata, only : ekh
     integer, intent(in) :: hi, hj, hk
     real, intent(in)    :: var(ib-hi:ie+hi,jb-hj:je+hj,kb-hk:ke+hk)
     real, intent(inout) :: rhs(ib-hi:ie+hi,jb-hj:je+hj,kb   :ke+hk)
     integer :: i, j, k, n

     do n = 1,bound_info_c%nbndptsrank
      !n = bound_info_c%bndptsrank(m)
         i = bound_info_c%bndpts_loc(n,1)
         j = bound_info_c%bndpts_loc(n,2)
         k = bound_info_c%bndpts_loc(n,3)

         if (abs(mask_c(i+1,j,k)) < eps1) then
           rhs(i,j,k) = rhs(i,j,k) - 0.5 * (ekh(i+1,j,k) + ekh(i,j,k)) * (var(i+1,j,k) - var(i,j,k))*dx2i
         end if

         if (abs(mask_c(i-1,j,k)) < eps1) then
           rhs(i,j,k) = rhs(i,j,k) + 0.5 * (ekh(i,j,k) + ekh(i-1,j,k)) * (var(i,j,k) - var(i-1,j,k))*dx2i
         end if

         if (abs(mask_c(i,j+1,k)) < eps1) then
           rhs(i,j,k) = rhs(i,j,k) - 0.5 * (ekh(i,j+1,k) + ekh(i,j,k)) * (var(i,j+1,k) - var(i,j,k))*dy2i
         end if

         if (abs(mask_c(i,j-1,k)) < eps1) then
           rhs(i,j,k) = rhs(i,j,k) + 0.5 * (ekh(i,j,k) + ekh(i,j-1,k)) * (var(i,j,k) - var(i,j-1,k))*dy2i
         end if

         if (abs(mask_c(i,j,k+1)) < eps1) then
           rhs(i,j,k) = rhs(i,j,k) - 0.5 * (dzf(k+1)*ekh(i,j,k) + dzf(k)*ekh(i,j,k+1)) &
                                         * (var(i,j,k+1) - var(i,j,k))*dzh2i(k+1)*dzfi(k)
         end if

         if (abs(mask_c(i,j,k-1)) < eps1) then
           rhs(i,j,k) = rhs(i,j,k) + 0.5 * (dzf(k-1)*ekh(i,j,k) + dzf(k)*ekh(i,j,k-1)) &
                                         * (var(i,j,k) - var(i,j,k-1))*dzh2i(k)*dzfi(k)
         end if

     end do

   end subroutine diffc_corr
#endif


   subroutine ibmwallfun
     use modglobal, only : libm, iwallmom, xhat, yhat, zhat, ltempeq, lmoist, &
                           ib, ie, ih, ihc, jb, je, jh, jhc, kb, ke, kh, khc, nsv, totheatflux, totqflux, nfcts, rk3step, timee, nfcts, lwritefac, dt, dtfac, tfac, tnextfac
     use modfields, only : thl0, qt0, sv0, up, vp, wp, thlp, qtp, svp, &
                           tau_x, tau_y, tau_z, thl_flux
     use modmpi, only : myid
     use modstat_nc, only : writestat_nc, writestat_1D_nc, writestat_2D_nc
#if defined(_GPU)
     use modcuda, only : up_d, vp_d, wp_d, thlp_d, qtp_d, svp_d, &
                         thl0_d, qt0_d, sv0_d, tau_x_d, tau_y_d, tau_z_d, thl_flux_d
#endif

     real, allocatable :: rhs(:,:,:)
     integer n

      if (.not. libm) return

#if defined(_GPU)
      if (iwallmom > 1) then
        call copy_tendency_device(rhs_ibm_d, up_d)
        call wallfunmom_dir_device(1, up_d)
        call accumulate_delta_device(tau_x_d, up_d, rhs_ibm_d)
        call reduce_fac_tau_device(1)

        call copy_tendency_device(rhs_ibm_d, vp_d)
        call wallfunmom_dir_device(2, vp_d)
        call accumulate_delta_device(tau_y_d, vp_d, rhs_ibm_d)
        call reduce_fac_tau_device(2)

        call copy_tendency_device(rhs_ibm_d, wp_d)
        call wallfunmom_dir_device(3, wp_d)
        call accumulate_delta_device(tau_z_d, wp_d, rhs_ibm_d)
        call reduce_fac_tau_device(3)
      end if

      call diffu_corr_device
      call diffv_corr_device
      call diffw_corr_device

      if (ltempeq .or. lmoist .or. lwritefac) then
        if (ltempeq) call copy_tendency_device(rhs_ibm_d, thlp_d)
        totheatflux = 0.
        totqflux = 0.
        call wallfunheat_dir_device
        if (ltempeq) call accumulate_delta_device(thl_flux_d, thlp_d, rhs_ibm_d)
        if (ltempeq) call diffc_corr_device(thl0_d, thlp_d, ih, jh, kh)
        if (lmoist)  call diffc_corr_device(qt0_d, qtp_d, ih, jh, kh)
        call reduce_fac_heat_device
      end if

      do n = 1,nsv
        call diffc_corr_device(sv0_d(:,:,:,n), svp_d(:,:,:,n), ihc, jhc, khc)
      end do
#else

      allocate(rhs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))

      if (iwallmom > 1) then
        rhs = up
        call wallfunmom(xhat, up, bound_info_u)
        tau_x(:,:,kb:ke+kh) = tau_x(:,:,kb:ke+kh) + (up - rhs)

        rhs = vp
        call wallfunmom(yhat, vp, bound_info_v)
        tau_y(:,:,kb:ke+kh) = tau_y(:,:,kb:ke+kh) + (vp - rhs)

        rhs = wp
        call wallfunmom(zhat, wp, bound_info_w)
        tau_z(:,:,kb:ke+kh) = tau_z(:,:,kb:ke+kh) + (wp - rhs)

        ! mom_flux_sum = sum(tau_x(ib:ie,jb:je,kb+1:ke) + tau_y(ib:ie,jb:je,kb+1:ke) + tau_z(ib:ie,jb:je,kb+1:ke))
        ! call MPI_ALLREDUCE(mom_flux_sum, mom_flux_tot, 1, MY_REAL, MPI_SUM, comm3d, mpierr)
        ! if (myid == 0) then
        !    if (rk3step == 3) then
        !         inquire(file="mom_flux.txt", exist=mom_flux_file_exists)
        !         if (mom_flux_file_exists) then
        !           open(12, file="mom_flux.txt", status="old", position="append", action="write")
        !         else
        !           open(12, file="mom_flux.txt", status="new", action="write")
        !         end if
        !         write(12, *) timee, -mom_flux_tot
        !         close(12)
        !    end if
        ! end if
      end if

      call diffu_corr
      call diffv_corr
      call diffw_corr

      if (ltempeq .or. lmoist .or. lwritefac) then
        rhs = thlp
        totheatflux = 0 ! Reset total heat flux to zero so we only account for that in this step.
        totqflux = 0
        call wallfunheat
        thl_flux(:,:,kb:ke+kh) = thl_flux(:,:,kb:ke+kh) + (thlp - rhs)
        if (ltempeq) call diffc_corr(thl0, thlp, ih, jh, kh)
        if (lmoist)  call diffc_corr(qt0, qtp, ih, jh, kh)

        ! thl_flux_sum = sum(thl_flux(ib:ie,jb:je,kb+1:ke))
        ! call MPI_ALLREDUCE(thl_flux_sum, thl_flux_tot, 1, MY_REAL, MPI_SUM, comm3d, mpierr)
        ! if (myid == 0) then
        !    if (rk3step == 3) then
        !         inquire(file="thl_flux.txt", exist=thl_flux_file_exists)
        !         if (thl_flux_file_exists) then
        !           open(12, file="thl_flux.txt", status="old", position="append", action="write")
        !         else
        !           open(12, file="thl_flux.txt", status="new", action="write")
        !         end if
        !         write(12, *) timee, thl_flux_tot
        !         close(12)
        !    end if
        ! end if
      end if

      do n = 1,nsv
        call diffc_corr(sv0(:,:,:,n), svp(:,:,:,n), ihc, jhc, khc)
      end do

      deallocate(rhs)
#endif

      if (lwritefac .and. rk3step==3) then
        if (myid == 0) then
            fac_tau_x_av = fac_tau_x_av + dt*fac_tau_x
            fac_tau_y_av = fac_tau_y_av + dt*fac_tau_y
            fac_tau_z_av = fac_tau_z_av + dt*fac_tau_z
            fac_pres_av = fac_pres_av + dt*fac_pres
            fac_pres2_av = fac_pres2_av + dt*fac_pres2
            fac_htc_av = fac_htc_av + dt*fac_htc
            fac_cth_av = fac_cth_av + dt*fac_cth

            if (timee >= tnextfac) then
               tfac = timee - tfac
               allocate(varsfac(nfcts,nstatfac))
               varsfac(:,1) = fac_tau_x_av(1:nfcts)/tfac
               varsfac(:,2) = fac_tau_y_av(1:nfcts)/tfac
               varsfac(:,3) = fac_tau_z_av(1:nfcts)/tfac
               varsfac(:,4) = fac_pres_av(1:nfcts)/tfac
               varsfac(:,5) = fac_htc_av(1:nfcts)/tfac
               varsfac(:,6) = fac_cth_av(1:nfcts)/tfac
               varsfac(:,7) = fac_pres2_av(1:nfcts)/tfac - (fac_pres_av(1:nfcts)/dtfac * fac_pres_av(1:nfcts)/tfac)
               call writestat_nc(ncidfac,1,tncstatfac,(/timee/),nrecfac,.true.)
               call writestat_1D_nc(ncidfac,nstatfac,ncstatfac,varsfac,nrecfac,nfcts)
               deallocate(varsfac)

               tfac = timee
               tnextfac = NINT((timee + dtfac))*1.0

               fac_tau_x_av = 0.
               fac_tau_y_av = 0.
               fac_tau_z_av = 0.
               fac_pres_av = 0.
               fac_pres2_av= 0.
               fac_htc_av = 0.
               fac_cth_av = 0.
            end if
        end if !myid
      end if
   end subroutine ibmwallfun


#if !defined(_GPU) || defined(UDALES_DEBUG)
   subroutine wallfunmom(dir, rhs, bound_info)
     use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh, &
                           dxdydzfi, iwallmom, xhat, yhat, zhat, vec0, nfcts, lwritefac, rk3step
     use modfields, only : u0, v0, w0, thl0
     use initfac,   only : facT, facnorm, faca
     use decomp_2d, only : zstart
     use modmpi,    only : comm3d, mpi_sum, mpierr, my_real

     real, intent(in)    :: dir(3)
     real, intent(inout) :: rhs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
     type(bound_info_type) :: bound_info

     integer i, j, k, sec, fac
     real dist, stress, stress_dir, area, momvol, Tair, logdz, &
          utan, ctm, a, a_is, a_xn, a_yn, a_zn, stress_ix, stress_iy, stress_iz
     real, dimension(3) :: uvec, norm, strm, span, stressvec
     logical :: valid
     real, dimension(1:nfcts) :: fac_tau_loc, fac_tau
     !real, dimension(:), allocatable :: fac_tau, fac_pres

     procedure(interp_velocity), pointer :: interp_velocity_ptr => null()
     procedure(interp_temperature), pointer :: interp_temperature_ptr => null()

     select case(alignment(dir))
     case(1)
       interp_velocity_ptr => interp_velocity_u
       interp_temperature_ptr => interp_temperature_u
     case(2)
       interp_velocity_ptr => interp_velocity_v
       interp_temperature_ptr => interp_temperature_v
     case(3)
       interp_velocity_ptr => interp_velocity_w
       interp_temperature_ptr => interp_temperature_w
     end select

     fac_tau_loc = 0.

     do sec = 1,bound_info%nfctsecsrank
       !sec = bound_info%fctsecsrank(m) ! index of section
       !n = bound_info%secbndptids(sec) ! index of boundary point
       area = bound_info%secareas_loc(sec) ! area of section
       fac = bound_info%secfacids_loc(sec) ! index of facet
       norm = facnorm(fac,:) ! facet normal

       if (bound_info%lskipsec_loc(sec)) cycle

       ! Fixed geometry, so read rather than computed; see build_wallfun_cache.
       ! The test is made before the interpolation rather than after it, as it
       ! once was: the interpolation has no side effect, so a section that fails
       ! here can skip the whole stencil instead of only the stress.
       logdz = bound_info%wfc_loc(sec,WFC_LOGDZ)
       if (logdz <= 1.) cycle

       i = bound_info%bndloc_loc(sec,1)
       j = bound_info%bndloc_loc(sec,2)
       k = bound_info%bndloc_loc(sec,3)

       dist = bound_info%secdist_loc(sec)
       if (bound_info%lcomprec_loc(sec) .or. lnorec) then
         uvec = interp_velocity_ptr(i, j, k)
         if (iwallmom == 2) then
           Tair = interp_temperature_ptr(i, j, k)
         end if
       else
         uvec(1) = trilinear_interp_cached(u0, bound_info%recloc_u(sec,:), &
                                           bound_info%recw_loc(sec,RECW_U:RECW_U+2))
         uvec(2) = trilinear_interp_cached(v0, bound_info%recloc_v(sec,:), &
                                           bound_info%recw_loc(sec,RECW_V:RECW_V+2))
         uvec(3) = trilinear_interp_cached(w0, bound_info%recloc_w(sec,:), &
                                           bound_info%recw_loc(sec,RECW_W:RECW_W+2))
         if (iwallmom == 2) &
           Tair = trilinear_interp_cached(thl0, bound_info%recloc_c(sec,:), &
                                          bound_info%recw_loc(sec,RECW_C:RECW_C+2))
       end if

       if (is_equal(uvec, vec0)) cycle

       call local_coords(uvec, norm, span, strm, valid)
       if (.not. valid) cycle

       utan = dot_product(uvec, strm)
       !utan = max(0.01, utan) ! uDALES 1

       ! calcualate momentum transfer coefficient
       ! make into interface somehow? because iwallmom doesn't change in the loop
       if (iwallmom == 2) then ! stability included
         ctm = mom_transfer_coef_stability(utan, dist, Tair, facT(fac,1), &
                                           logdz, bound_info%wfc_loc(sec,WFC_LOGZH), &
                                           bound_info%wfc_loc(sec,WFC_COEF), &
                                           bound_info%wfc_loc(sec,WFC_CM), &
                                           bound_info%wfc_loc(sec,WFC_CH))
       else if (iwallmom == 3) then ! neutral
         ctm = bound_info%wfc_loc(sec,WFC_CTMNEU)
       end if

       stress = ctm * utan**2

       if (bound_info%lcomprec_loc(sec)) then
         a = dot_product(dir, strm)
         stress_dir = a * stress
       else
         ! Rotation from local (strm,span,norm) to global (xhat,yhat,zhat) basis
         ! \tau'_ij = a_ip a_jq \tau_pq
         ! \tau_pq in local coordinates is something like \tau \delta_13, because we only have \tau_{strm,norm})
         a_is = dot_product(dir, strm)
         a_xn = dot_product(xhat, norm)
         a_yn = dot_product(yhat, norm)
         a_zn = dot_product(zhat, norm)

         stress_ix = a_is * a_xn * stress
         stress_iy = a_is * a_yn * stress
         stress_iz = a_is * a_zn * stress

         stressvec(1) = stress_ix
         stressvec(2) = stress_iy
         stressvec(3) = stress_iz
         stress_dir = norm2(stressvec)
       end if

       stress_dir = sign(stress_dir, dot_product(uvec, dir))

       ! dxdydzfi is 1/(dx*dy*dzf(k)), formed once in initglobal. The device
       ! kernel multiplies by the same array, so the two stay bit-identical.
       momvol = stress_dir * area * dxdydzfi(k)
       rhs(i,j,k) = rhs(i,j,k) - momvol
       fac_tau_loc(fac) = fac_tau_loc(fac) + stress_dir * area ! output stresses on facets
     end do

     if (allocated(fac_tau_raw)) fac_tau_raw = fac_tau_loc(1:nfcts)

     if (lwritefac .and. rk3step==3) then
        fac_tau_loc(1:nfcts) = fac_tau_loc(1:nfcts) / faca(1:nfcts)
        call MPI_ALLREDUCE(fac_tau_loc(1:nfcts), fac_tau(1:nfcts), nfcts, MY_REAL, MPI_SUM, comm3d, mpierr)

        select case(alignment(dir))
        case(1)
           fac_tau_x = fac_tau
        case(2)
           fac_tau_y = fac_tau
        case(3)
           fac_tau_z = fac_tau
        end select
     end if

     ! Do time-averaging like in modEB

   end subroutine wallfunmom


   subroutine wallfunheat
     use modglobal, only : ib, ie, jb, je, dxdydzhi, &
                           xhat, yhat, zhat, vec0, ltempeq, lmoist, iwalltemp, iwallmoist, lEB, lwritefac, nfcts, rk3step, totheatflux, totqflux
     use modfields, only : u0, v0, w0, thl0, thlp, qt0, qtp, pres0
     use initfac,   only : facT, facnorm, fachf, facef, facqsat, fachurel, facf, faclGR, faca
     use modmpi,    only : comm3d, mpi_sum, mpierr, my_real
     use modibmdata, only : bctfxm, bctfxp, bctfyp, bctfz

     integer i, j, k, sec, fac
     real :: dist, flux, area, Tair, utan, cth, htc, logdz, cveg, hurel, qtair, qwall, resa, resc, ress
     real, dimension(3) :: uvec, norm, span, strm
     real, dimension(1:nfcts) :: fac_htc_loc, fac_cth_loc, fac_pres_loc, fac_pres2_loc
     logical :: valid

     fac_htc_loc = 0.
     fac_cth_loc = 0.
     fac_pres_loc = 0.
     fac_pres2_loc = 0.

     do sec = 1,bound_info_c%nfctsecsrank
       ! sec = bound_info_c%fctsecsrank(m) ! index of section
       !n =   bound_info_c%secbndptids(sec) ! index of boundary point
       fac = bound_info_c%secfacids_loc(sec) ! index of facet
       area = bound_info_c%secareas_loc(sec) ! area
       norm = facnorm(fac,:)

       i = bound_info_c%bndloc_loc(sec,1)
       j = bound_info_c%bndloc_loc(sec,2)
       k = bound_info_c%bndloc_loc(sec,3)

       fac_pres_loc(fac) = fac_pres_loc(fac) + pres0(i,j,k) * area ! output pressure on facets
       fac_pres2_loc(fac) = fac_pres2_loc(fac) +  pres0(i,j,k)* pres0(i,j,k) * area

       if (bound_info_c%lskipsec_loc(sec)) cycle

       ! Fixed geometry; see build_wallfun_cache. Tested before the
       ! interpolation rather than after it, which is free of side effects.
       logdz = bound_info_c%wfc_loc(sec,WFC_LOGDZ)
       if (logdz <= 1.) cycle

       dist = bound_info_c%secdist_loc(sec)
       if (bound_info_c%lcomprec_loc(sec) .or. lnorec) then ! section aligned with grid - use this cell's velocity
         uvec = interp_velocity_c(i, j, k)
         Tair = thl0(i,j,k)
         qtair = qt0(i,j,k)
       else ! use velocity at reconstruction point
         uvec(1) = trilinear_interp_cached(u0, bound_info_c%recloc_u(sec,:), &
                                           bound_info_c%recw_loc(sec,RECW_U:RECW_U+2))
         uvec(2) = trilinear_interp_cached(v0, bound_info_c%recloc_v(sec,:), &
                                           bound_info_c%recw_loc(sec,RECW_V:RECW_V+2))
         uvec(3) = trilinear_interp_cached(w0, bound_info_c%recloc_w(sec,:), &
                                           bound_info_c%recw_loc(sec,RECW_W:RECW_W+2))
         Tair  = trilinear_interp_cached(thl0, bound_info_c%recloc_c(sec,:), &
                                         bound_info_c%recw_loc(sec,RECW_C:RECW_C+2))
         qtair = trilinear_interp_cached( qt0, bound_info_c%recloc_c(sec,:), &
                                         bound_info_c%recw_loc(sec,RECW_C:RECW_C+2))
       end if

       if (is_equal(uvec, vec0)) cycle

       call local_coords(uvec, norm, span, strm, valid)
       if (.not. valid) cycle
       utan = dot_product(uvec, strm)
       !utan = max(0.01, utan) ! uDALES 1

       ! Sensible heat
       if (ltempeq) then
         if (iwalltemp == 1) then ! probably remove this eventually, only relevant to grid-aligned facets
           !if     (all(abs(norm - xhat) < eps1)) then
           if     (is_equal(norm, xhat)) then
             flux = bctfxp
           !elseif (all(abs(norm + xhat) < eps1)) then
           elseif (is_equal(norm, -xhat)) then
             flux = bctfxm
           !elseif (all(abs(norm - yhat) < eps1)) then
           elseif (is_equal(norm, yhat)) then
             flux = bctfyp
           !elseif (all(abs(norm + yhat) < eps1)) then
           elseif (is_equal(norm, -yhat)) then
             flux = bctfxm
           !elseif (all(abs(norm - zhat) < eps1)) then
           elseif (is_equal(norm, zhat)) then
             flux = bctfz
           end if

         elseif (iwalltemp == 2) then
           call heat_transfer_coef_flux(utan, dist, Tair, facT(fac, 1), &
                                        logdz, bound_info_c%wfc_loc(sec,WFC_LOGZH), &
                                        bound_info_c%wfc_loc(sec,WFC_COEF), &
                                        bound_info_c%wfc_loc(sec,WFC_CM), &
                                        bound_info_c%wfc_loc(sec,WFC_CH), &
                                        cth, flux, htc)
           fac_cth_loc(fac) = fac_cth_loc(fac) + cth * area ! output heat transfer coefficients on facets
           fac_htc_loc(fac) = fac_htc_loc(fac) + htc * area ! output heat transfer coefficients on facets
         end if

         ! flux [Km/s]
         ! fluid volumetric sensible heat source/sink = flux * area / volume [K/s]
         ! facet sensible heat flux = volumetric heat capacity of air * flux * sectionarea / facetarea [W/m^2]
         thlp(i,j,k) = thlp(i,j,k) - flux * area * dxdydzhi(k)

         totheatflux = totheatflux + flux*area ! [Km^3s^-1] This sums the flux over all facets (unconditional, mirrors totqflux; decouples periodicEBcorr from lEB)

         if (lEB) then
           fachf(fac) = fachf(fac) + flux * area ! [Km^2/s] (will be divided by facetarea(fac) in modEB)
         end if !fachf=[Km/s]
       end if

       ! Latent heat
       if (lmoist .and. faclGR(fac)) then
         if (iwallmoist == 1) then ! probably remove this eventually, only relevant to grid-aligned facets
           if     (is_equal(norm, xhat)) then
             flux = bcqfxp
           elseif (is_equal(norm, -xhat)) then
             flux = bcqfxm
           elseif (is_equal(norm, yhat)) then
             flux = bcqfyp
           elseif (is_equal(norm, -yhat)) then
             flux = bcqfym
           elseif (is_equal(norm, zhat)) then
             flux = bcqfz
           end if

         elseif (iwallmoist == 2) then
            if (abs(htc*abs(utan)) > 0.) then
               qwall = facqsat(fac) ! saturation humidity
               hurel = fachurel(fac) ! relative humidity
               resa = 1./(htc*abs(utan)) ! aerodynamic resistance
               resc = facf(fac,4) ! canopy resistance
               ress = facf(fac,5) ! soil resistance
               cveg = 0.8 ! vegetation fraction
               flux = moist_flux(cveg, resa, qtair, qwall, hurel, resc, ress)
            end if
         end if

         ! flux [kg/kg m/s]
         ! fluid volumetric latent heat source/sink = flux * area / volume [kg/kg / s]
         ! facet latent heat flux = volumetric heat capacity of air * flux * sectionarea / facetarea [W/m^2]
         totqflux = totqflux + flux*area ! [Km^3s^-1] This sums the flux over all facets
         qtp(i,j,k) = qtp(i,j,k) - flux * area * dxdydzhi(k)

         if (lEB) then
           facef(fac) = facef(fac) + flux * area ! [Km^2/s] (will be divided by facetarea(fac) in modEB)
         end if
       end if

     end do

     ! Published before the area normalisation, so that the device port has the
     ! raw accumulators to be compared against. See fac_tau_raw.
     if (allocated(fac_htc_raw)) then
       fac_htc_raw   = fac_htc_loc(1:nfcts)
       fac_cth_raw   = fac_cth_loc(1:nfcts)
       fac_pres_raw  = fac_pres_loc(1:nfcts)
       fac_pres2_raw = fac_pres2_loc(1:nfcts)
     end if

     if (lwritefac .and. rk3step==3) then
        fac_cth_loc(1:nfcts) = fac_cth_loc(1:nfcts) / faca(1:nfcts)
        fac_htc_loc(1:nfcts) = fac_htc_loc(1:nfcts) / faca(1:nfcts)
        fac_pres_loc(1:nfcts) = fac_pres_loc(1:nfcts) / faca(1:nfcts)
        fac_pres2_loc(1:nfcts) = fac_pres2_loc(1:nfcts) / faca(1:nfcts)
        
        call MPI_ALLREDUCE(fac_cth_loc(1:nfcts), fac_cth(1:nfcts), nfcts, MY_REAL, MPI_SUM, comm3d, mpierr)
        call MPI_ALLREDUCE(fac_htc_loc(1:nfcts), fac_htc(1:nfcts), nfcts, MY_REAL, MPI_SUM, comm3d, mpierr)
        call MPI_ALLREDUCE(fac_pres_loc(1:nfcts), fac_pres(1:nfcts), nfcts, MY_REAL, MPI_SUM, comm3d, mpierr)
        call MPI_ALLREDUCE(fac_pres2_loc(1:nfcts), fac_pres2(1:nfcts), nfcts, MY_REAL, MPI_SUM, comm3d, mpierr)
     end if

   end subroutine wallfunheat


   !> Trilinear interpolation of a field at a reconstruction point.
   !!
   !! cell holds rank-local indices and d the normalised offsets within it, both
   !! precomputed per section by build_wallfun_cache. What is left is an
   !! eight-point gather: the grid arrays, the global-to-local shift and the
   !! three divisions the offsets cost are all settled at initialisation, as is
   !! the off-rank check this used to make on every call.
   real function trilinear_interp_cached(var, cell, d)
     use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh
     implicit none
     real,    intent(in) :: var(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
     integer, intent(in) :: cell(3)
     real,    intent(in) :: d(3)
     real, dimension(8)  :: corners
     real :: xd, yd, zd

     corners = eval_corners(var, cell(1), cell(2), cell(3))
     xd = d(1); yd = d(2); zd = d(3)

     trilinear_interp_cached = corners(1) * (1-xd)*(1-yd)*(1-zd) + & ! c000
                               corners(2) * (  xd)*(1-yd)*(1-zd) + & ! c100
                               corners(3) * (1-xd)*(  yd)*(1-zd) + & ! c010
                               corners(4) * (  xd)*(  yd)*(1-zd) + & ! c110
                               corners(5) * (1-xd)*(1-yd)*(  zd) + & ! c001
                               corners(6) * (  xd)*(1-yd)*(  zd) + & ! c101
                               corners(7) * (1-xd)*(  yd)*(  zd) + & ! c011
                               corners(8) * (  xd)*(  yd)*(  zd)     ! c111

   end function trilinear_interp_cached

   function eval_corners(var, i, j, k)
     use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh
     integer, intent(in) :: i, j, k ! LOCAL indices
     real, intent(in)    :: var(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
     real, dimension(8)  :: eval_corners(8)

     eval_corners(1) = var(i  ,j  ,k  ) !c000
     eval_corners(2) = var(i+1,j  ,k  ) !c100
     eval_corners(3) = var(i  ,j+1,k  ) !c010
     eval_corners(4) = var(i+1,j+1,k  ) !c110
     eval_corners(5) = var(i  ,j  ,k+1) !c001
     eval_corners(6) = var(i+1,j  ,k+1) !c101
     eval_corners(7) = var(i  ,j+1,k+1) !c011
     eval_corners(8) = var(i+1,j+1,k+1) !c111

   end function eval_corners
#endif


   integer function alignment(n)
     ! returns an integer determining whether a unit vector n is aligned with the
     ! coordinates axes.
     use modglobal, only : xhat, yhat, zhat
     implicit none
     real, dimension(3), intent(in) :: n ! must be unit vector

     if     (is_equal(n, xhat)) then
       alignment = 1
     elseif (is_equal(n, yhat)) then
       alignment = 2
     elseif (is_equal(n, zhat)) then
       alignment = 3
     elseif (is_equal(n, -xhat)) then
       alignment = -1
     elseif (is_equal(n, -yhat)) then
       alignment = -2
     elseif (is_equal(n, -zhat)) then
       alignment = -3
     else
       alignment = 0
     end if

   end function alignment

   ! overload this with real dimension(1)
   ! could put somewhere else because not specific to ibm
   logical function is_equal(a,b)
     ! determines whether two vectors are equal to each other within a tolerance of eps1
     use modglobal, only : eps1
     implicit none
     real, dimension(3), intent(in) :: a, b

     if (all(abs(a - b) < eps1)) then
       is_equal = .true.
     else
       is_equal = .false.
     end if

   end function is_equal


#if !defined(_GPU) || defined(UDALES_DEBUG)
   function cross_product(a,b)
     ! Calculate the cross product (a x b)
     implicit none
     real, dimension(3) :: cross_product
     real, dimension(3), intent(in) :: a, b

     cross_product(1) = a(2)*b(3) - a(3)*b(2)
     cross_product(2) = a(3)*b(1) - a(1)*b(3)
     cross_product(3) = a(1)*b(2) - a(2)*b(1)

   end function cross_product


   function interp_velocity_u(i, j, k)
     ! interpolates the velocity at u-grid location i,j,k
     use modfields, only :  u0, v0, w0
     real ::  interp_velocity_u(3)
     integer, intent(in) :: i, j, k

     interp_velocity_u(1) = u0(i,j,k)
     interp_velocity_u(2) = 0.25 * (v0(i,j,k) + v0(i,j+1,k) + v0(i-1,j,k) + v0(i-1,j+1,k))
     interp_velocity_u(3) = 0.25 * (w0(i,j,k) + w0(i,j,k+1) + w0(i-1,j,k) + w0(i-1,j,k+1)) !only for equidistant grid!

     return
   end function interp_velocity_u


   function interp_velocity_v(i, j, k)
     ! interpolates the velocity at v-grid location i,j,k
     use modfields, only :  u0, v0, w0
     real ::  interp_velocity_v(3)
     integer, intent(in) :: i, j, k

     interp_velocity_v(1) = 0.25 * (u0(i,j,k) + u0(i+1,j,k) + u0(i,j-1,k) + u0(i+1,j-1,k))
     interp_velocity_v(2) = v0(i,j,k)
     interp_velocity_v(3) = 0.25 * (w0(i,j,k) + w0(i,j,k+1) + w0(i,j-1,k) + w0(i,j-1,k+1)) !only for equidistant grid!

     return
   end function interp_velocity_v


   function interp_velocity_w(i, j, k)
     ! interpolates the velocity at w-grid location i,j,k
     use modfields, only :  u0, v0, w0
     real ::  interp_velocity_w(3)
     integer, intent(in) :: i, j, k

     interp_velocity_w(1) = 0.25 * (u0(i,j,k) + u0(i+1,j,k) + u0(i,j-1,k) + u0(i+1,j-1,k))
     interp_velocity_w(2) = v0(i,j,k)
     interp_velocity_w(3) = 0.25 * (w0(i,j,k) + w0(i,j,k+1) + w0(i,j-1,k) + w0(i,j-1,k+1)) !only for equidistant grid!

     return
   end function interp_velocity_w


   function interp_velocity_c(i, j, k)
     ! interpolates the velocity at c-grid location i,j,k
     use modfields, only :  u0, v0, w0
     real ::  interp_velocity_c(3)
     integer, intent(in) :: i, j, k

     interp_velocity_c(1) = 0.5 * (u0(i,j,k) + u0(i+1,j,k))
     interp_velocity_c(2) = 0.5 * (v0(i,j,k) + v0(i,j+1,k))
     interp_velocity_c(3) = 0.5 * (w0(i,j,k) + w0(i,j,k+1))

     return
   end function interp_velocity_c


   real function interp_temperature_u(i, j, k)
     ! interpolates the temperature at u-grid location i,j,k
     use modfields, only :  thl0
     integer, intent(in) :: i, j, k

     !interp_temperature_u = 0.5 * (thl0(i,j,k) + thl0(i-1,j,k))
     interp_temperature_u = 0.5 * (thl0(i  ,j,k)*mask_c(i  ,j,k)*(2.-mask_c(i-1,j,k)) &
                                +  thl0(i-1,j,k)*mask_c(i-1,j,k)*(2.-mask_c(i  ,j,k)))

     return
   end function interp_temperature_u


   real function interp_temperature_v(i, j, k)
     ! interpolates the temperature at v-grid location i,j,k
     use modfields, only :  thl0
     integer, intent(in) :: i, j, k

     !interp_temperature_v = 0.5 * (thl0(i,j,k) + thl0(i,j-1,k))
     interp_temperature_v = 0.5 * (thl0(i,j  ,k)*mask_c(i,j  ,k)*(2.-mask_c(i,j-1,k)) &
                                 + thl0(i,j-1,k)*mask_c(i,j-1,k)*(2.-mask_c(i,j  ,k)))

     return
   end function interp_temperature_v


   real function interp_temperature_w(i, j, k)
     ! interpolates the temperature at w-grid location i,j,k
     use modfields, only :  thl0
     integer, intent(in) :: i, j, k

     !interp_temperature_w = 0.5 * (thl0(i,j,k) + thl0(i,j,k-1))
     interp_temperature_w = 0.5 * (thl0(i,j,k  )*mask_c(i,j,k  )*(2.-mask_c(i,j,k-1)) &
                                +  thl0(i,j,k-1)*mask_c(i,j,k-1)*(2.-mask_c(i,j,k  )))

     return
   end function interp_temperature_w


   subroutine local_coords(uvec, norm, span, strm, valid)
     ! returns the local streamwise (strm) and spanwise vectors (span) in the
     ! plane normal to (norm) containing the velocity vector (uvec)
     use modglobal, only : vec0
     real, intent(in),  dimension(3) :: uvec, norm
     real, intent(out), dimension(3) :: span, strm
     logical, intent(out) :: valid

     span = cross_product(norm, uvec)
     !if (is_equal(span, (/0.,0.,0./))) then
     ! velocity is pointing into or outof the surface
     if (is_equal(span, vec0)) then
       strm = 0.
       valid = .false.
     else
       ! One reciprocal and three multiplies rather than three divisions; the
       ! device kernel normalises the same way, so the two agree bit for bit.
       span = span * (1./norm2(span))
       valid = .true.
     end if
     strm = cross_product(span, norm)

   end subroutine local_coords


   real function mom_transfer_coef_stability(utan, dist, Tair, Tsurf, &
                                             logdz, logzh, coef, cmfac, chfac)
     ! By Ivo Suter. calculates the momentum transfer coefficient based on the
     ! surface tangential velocity 'utan' at a distance 'dist' from the surface.
     ! Stability are included using the air temperature Tair and surface temperature Tsurf.
     !
     ! The roughness terms arrive already evaluated, from the per-section cache
     ! build_wallfun_cache fills: logdz = log(dist/z0), logzh = log(z0/z0h),
     ! coef = fkar2/logdz**2 and cmfac, chfac the two (d*fkar2)/logdz**2*b1*sqdz
     ! products. None of them depends on the solution, and together they were
     ! every logarithm and the square root this routine evaluated per call. The
     ! neutral coefficient (fkar/logdz)**2 is cached whole, so the routine
     ! mom_transfer_coef_neutral used to provide it is gone.
     use modglobal, only : grav, prandtlturb

      implicit none
      real, intent(in) :: dist, Tsurf, Tair, utan
      real, intent(in) :: logdz, logzh, coef, cmfac, chfac
      real, parameter :: b1 = 9.4 !parameters from uno1995
      real, parameter :: b2 = 4.7
      real :: dT, Ribl0, Ribl1, Fm, Fh, M

      dT = Tair - Tsurf
      Ribl0 = grav * dist * dT / (Tsurf * utan**2) !Eq. 6, guess initial Ri

      IF (Ribl0 > 0.) THEN !0.25 approx critical for bulk Richardson number  => stable
         Fm = 1./(1. + b2*Ribl0)**2 !Eq. 4
         Fh = Fm !Eq. 4
      ELSE ! => unstable
         Fm = 1. - (b1*Ribl0)/(1. + cmfac*SQRT(ABS(Ribl0))) !Eq. 3
         Fh = 1. - (b1*Ribl0)/(1. + chfac*SQRT(ABS(Ribl0))) !Eq. 3
      END IF

      M = prandtlturb*logdz*SQRT(Fm)/Fh !Eq. 14

      Ribl1 = Ribl0 - Ribl0*prandtlturb*logzh/(prandtlturb*logzh + M) !Eq. 17

      !interate to get new Richardson number
      IF (Ribl1 > 0.) THEN !0.25 approx critical for bulk Richardson number  => stable
         Fm = 1./(1. + b2*Ribl1)**2 !Eq. 4
      ELSE ! => unstable
         Fm = 1. - (b1*Ribl1)/(1. + cmfac*SQRT(ABS(Ribl1))) !Eq. 3
      END IF

      mom_transfer_coef_stability = coef*Fm !Eq. 7

   end function mom_transfer_coef_stability


   subroutine heat_transfer_coef_flux(utan, dist, Tair, Tsurf, &
                                      logdz, logzh, coef, cmfac, chfac, cth, flux, htc)
     ! The roughness terms arrive already evaluated; see
     ! mom_transfer_coef_stability for what they are and why.
     use modglobal, only : grav, prandtlturb

      implicit none
      real, intent(in)  :: dist, Tsurf, Tair, utan
      real, intent(in)  :: logdz, logzh, coef, cmfac, chfac
      real, intent(out) :: cth, flux, htc
      real, parameter :: b1 = 9.4 !parameters from Uno1995
      real, parameter :: b2 = 4.7
      !real :: Pr
      real :: dT, Ribl0, Ribl1, Fm, Fh, M, dTrough

      !Pr = 1.
      !Pr = prandtlmol
      dT = Tair - Tsurf
      Ribl0 = grav * dist * dT / (Tsurf * utan**2) !Eq. 6, guess initial Ri

      cth = 0.
      flux = 0.
      if (Ribl0 > 0.) then
         Fm = 1./(1. + b2*Ribl0)**2 !Eq. 4
         Fh = Fm !Eq. 4
      else ! => unstable
         Fm = 1. - (b1*Ribl0)/(1. + cmfac*sqrt(abs(Ribl0))) !Eq. 3
         Fh = 1. - (b1*Ribl0)/(1. + chfac*sqrt(abs(Ribl0))) !Eq. 3
      end if

      M = prandtlturb*logdz*sqrt(Fm)/Fh !Eq. 14
      Ribl1 = Ribl0 - Ribl0*prandtlturb*logzh/(prandtlturb*logzh + M) !Eq. 17

      !interate to get new Richardson number
      if (Ribl1 > 0.) then
         Fm = 1./(1. + b2*Ribl1)**2 !Eq. 4
         Fh = Fm !Eq. 4
      else ! => unstable
         Fm = 1. - (b1*Ribl1)/(1. + cmfac*sqrt(abs(Ribl1))) !Eq. 3
         Fh = 1. - (b1*Ribl1)/(1. + chfac*sqrt(abs(Ribl1))) !Eq. 3
      end if

      ! Uno (2)
      M = prandtlturb*logdz*sqrt(Fm)/Fh !Eq. 14
      dTrough = dT*1./(prandtlturb*logzh/M + 1.) !Eq. 13a
      cth = coef*Fh/prandtlturb ! Ivo's heat transfer coefficient
      flux = abs(utan)*cth*dTrough

      if (abs(abs(utan)*dT) > 0.) then
         htc = flux / (abs(utan)*dT)
      else
         htc = 0.
      end if

      ! ! Uno (8)
      ! cth = abs(utan)*fkar2/(logdz*logdzh)*Fh/prandtlturb !Eq. 8
      ! flux = cth*dT !Eq. 2, Eq. 8

   end subroutine heat_transfer_coef_flux


   real function moist_flux(cveg, resa, qtair, qwall, hurel, resc, ress)
     real, intent(in) :: cveg, resa, qtair, qwall, hurel, resc, ress

     moist_flux = min(0., cveg * (qtair - qwall)         / (resa + resc) + &
                     (1 - cveg)* (qtair - qwall * hurel) / (resa + ress))

   end function moist_flux
#endif


#if defined(_GPU)
   attributes(global) subroutine bottom_set_tau_cuda
      use modcuda, only : ih_d, jh_d, kh_d, ie_d, je_d, kb_d, ke_d, &
                          loneeqn_d, ltempeq_d, &
                          e120_d, e12m_d, up_d, vp_d, wp_d, thlp_d, &
                          tau_x_d, tau_y_d, tau_z_d, thl_flux_d, &
                          tidandstride
      implicit none

      integer :: i, j, k, tidx, tidy, tidz, stridex, stridey, stridez

      call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

      if (loneeqn_d) then
         if (tidz == kb_d) then
            do j = tidy - jh_d, je_d + jh_d, stridey
               do i = tidx - ih_d, ie_d + ih_d, stridex
                  e120_d(i, j, kb_d - 1) = e120_d(i, j, kb_d)
                  e12m_d(i, j, kb_d - 1) = e12m_d(i, j, kb_d)
               end do
            end do
         end if
      end if

      do k = tidz, ke_d + kh_d, stridez
         do j = tidy - jh_d, je_d + jh_d, stridey
            do i = tidx - ih_d, ie_d + ih_d, stridex
               tau_x_d(i,j,k) = up_d(i,j,k)
               tau_y_d(i,j,k) = vp_d(i,j,k)
               tau_z_d(i,j,k) = wp_d(i,j,k)
            end do
         end do
      end do

      if (ltempeq_d) then
         do k = tidz, ke_d + kh_d, stridez
            do j = tidy - jh_d, je_d + jh_d, stridey
               do i = tidx - ih_d, ie_d + ih_d, stridex
                  thl_flux_d(i,j,k) = thlp_d(i,j,k)
               end do
            end do
         end do
      end if
   end subroutine bottom_set_tau_cuda

   attributes(global) subroutine bottom_update_thlp_cuda(wtsurf)
      use modcuda, only : ie_d, je_d, kb_d, dzf_d, dzfi_d, dzh2i_d, &
                          ekh_d, thl0_d, thlp_d, &
                          tidandstride
      implicit none
      real, value, intent(in) :: wtsurf
      integer :: i, j, tidx, tidy, tidz, stridex, stridey, stridez

      call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

      if (tidz == kb_d) then
         do j = tidy, je_d, stridey
            do i = tidx, ie_d, stridex
               thlp_d(i, j, kb_d) = thlp_d(i, j, kb_d) + ( &
                                    0.5*(dzf_d(kb_d - 1)*ekh_d(i, j, kb_d) + dzf_d(kb_d)*ekh_d(i, j, kb_d - 1)) &
                                    *(thl0_d(i, j, kb_d) - thl0_d(i, j, kb_d - 1)) &
                                    *dzh2i_d(kb_d) &
                                    - wtsurf &
                                    )*dzfi_d(kb_d)
            end do
         end do
     end if
   end subroutine bottom_update_thlp_cuda

   attributes(global) subroutine bottom_update_qtp_cuda(wqsurf)
      use modcuda, only : ie_d, je_d, kb_d, dzf_d, dzfi_d, dzh2i_d, &
                          ekh_d, qt0_d, qtp_d, &
                          tidandstride
      implicit none
      real, value, intent(in) :: wqsurf
      integer :: i, j, tidx, tidy, tidz, stridex, stridey, stridez

      call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

      if (tidz == kb_d) then  
         do j = tidy, je_d, stridey
            do i = tidx, ie_d, stridex
               qtp_d(i, j, kb_d) = qtp_d(i, j, kb_d) + ( &
                                   0.5*(dzf_d(kb_d - 1)*ekh_d(i, j, kb_d) + dzf_d(kb_d)*ekh_d(i, j, kb_d - 1)) &
                                   *(qt0_d(i, j, kb_d) - qt0_d(i, j, kb_d - 1)) &
                                   *dzh2i_d(kb_d) &
                                   + wqsurf &
                                   )*dzfi_d(kb_d)
            end do
         end do
      end if
   end subroutine bottom_update_qtp_cuda

   attributes(global) subroutine bottom_update_svp_cuda(m)
      use modcuda, only : ie_d, je_d, kb_d, dzf_d, dzfi_d, dzh2i_d, &
                          ekh_d, sv0_d, svp_d, &
                          tidandstride
      implicit none
      integer, value, intent(in) :: m
      integer :: i, j, tidx, tidy, tidz, stridex, stridey, stridez

      call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

      if (tidz == kb_d) then
         do j = tidy, je_d, stridey
            do i = tidx, ie_d, stridex
               svp_d(i, j, kb_d, m) = svp_d(i, j, kb_d, m) + ( &
                                      0.5*(dzf_d(kb_d - 1)*ekh_d(i, j, kb_d) + dzf_d(kb_d)*ekh_d(i, j, kb_d - 1)) &
                                      *(sv0_d(i, j, kb_d, m) - sv0_d(i, j, kb_d - 1, m)) &
                                      *dzh2i_d(kb_d) &
                                      + 0. &
                                      )*dzfi_d(kb_d)
            end do
         end do
      end if
   end subroutine bottom_update_svp_cuda

   attributes(global) subroutine bottom_update_tau_cuda
      use modcuda, only : ih_d, jh_d, kh_d, ie_d, je_d, ke_d, &
                          ltempeq_d, &
                          up_d, vp_d, wp_d, thlp_d, &
                          tau_x_d, tau_y_d, tau_z_d, thl_flux_d, &
                          tidandstride
      implicit none

      integer :: i, j, k, tidx, tidy, tidz, stridex, stridey, stridez

      call tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

      do k = tidz, ke_d + kh_d, stridez
         do j = tidy - jh_d, je_d + jh_d, stridey
            do i = tidx - ih_d, ie_d + ih_d, stridex
               tau_x_d(i,j,k) = up_d(i,j,k) - tau_x_d(i,j,k)
               tau_y_d(i,j,k) = vp_d(i,j,k) - tau_y_d(i,j,k)
               tau_z_d(i,j,k) = wp_d(i,j,k) - tau_z_d(i,j,k)
            end do
         end do
      end do

      if (ltempeq_d) then
         do k = tidz, ke_d + kh_d, stridez
            do j = tidy - jh_d, je_d + jh_d, stridey
               do i = tidx - ih_d, ie_d + ih_d, stridex
                  thl_flux_d(i,j,k) = thlp_d(i,j,k) - thl_flux_d(i,j,k)
               end do
            end do
         end do
      end if
   end subroutine bottom_update_tau_cuda
#endif

   subroutine bottom
      ! By Ivo Suter.
      !kind of obsolete when road facets are being used
      !vegetated floor not added (could simply be copied from vegetated horizontal facets)
      use modglobal,   only: ih, jh, kh, lbottom, ltempeq, lmoist, nsv, BCbotm, BCbotT, BCbotq, BCbots
      use modsurfdata, only: thls, z0, z0h, wtsurf, wqsurf
#if defined(_GPU)
      use cudafor
      use modcuda,              only: griddim, blockdim, checkCUDA, u0_d, v0_d, thl0_d, up_d, vp_d, thlp_d, &
                                      momfluxb_d, tfluxb_d, bcTfluxA_d
      use modwallfunctions,     only: wfmneutral_cuda, wfuno_cuda
#else
      use modglobal,            only: ib, ie, jb, je, kb, ke, dzf, dzfi, dzh2i
      use modfields,            only: u0, v0, e120, thl0, qt0, sv0, e12m, up, vp, wp, thlp, qtp, svp, &
                                      momfluxb, tfluxb, tau_x, tau_y, tau_z, thl_flux
      use modsubgriddata,       only: ekh
      use modwallfunctions,     only: wfmneutral, wfuno
#endif
      implicit none
      integer :: i, j, m

#if defined(_GPU)
      call bottom_set_tau_cuda<<<griddim,blockdim>>>
      call checkCUDA( cudaGetLastError(), 'bottom_set_tau_cuda' )
#else
      e120(:, :, kb - 1) = e120(:, :, kb)
      e12m(:, :, kb - 1) = e12m(:, :, kb)
      ! wm(:, :, kb) = 0. ! SO moved to modboundary
      ! w0(:, :, kb) = 0.
      tau_x(:,:,kb:ke+kh) = up
      tau_y(:,:,kb:ke+kh) = vp
      tau_z(:,:,kb:ke+kh) = wp
      thl_flux(:,:,kb:ke+kh) = thlp
#endif

      if (lbottom) then
         !momentum
         if (BCbotm.eq.2) then
#if defined(_GPU)
            bcTfluxA_d = 0.
            call wfuno_cuda<<<griddim,blockdim>>>(ih, jh, kh, up_d, vp_d, thlp_d, &
                                                  momfluxb_d, tfluxb_d, bcTfluxA_d, &
                                                  u0_d, v0_d, thl0_d, thls, z0, z0h, 91)
            call checkCUDA( cudaGetLastError(), 'wfuno_cuda under BCbotm' )
            bcTfluxA = bcTfluxA_d
#else
            call wfuno(ih, jh, kh, up, vp, thlp, momfluxb, tfluxb, bcTfluxA, &
                       u0, v0, thl0, thls, z0, z0h, 91)
#endif
         elseif (BCbotm.eq.3) then
#if defined(_GPU)
            call wfmneutral_cuda<<<griddim,blockdim>>>(ih, jh, kh, up_d, vp_d, momfluxb_d, u0_d, v0_d, z0, 91)
            call checkCUDA( cudaGetLastError(), 'wfmneutral_cuda' )
#else
            call wfmneutral(ih, jh, kh, up, vp, momfluxb, u0, v0, z0, 91)
#endif
         else
            write(0, *) "ERROR: bottom boundary type for momentum undefined"
            stop 1
         end if


         if (ltempeq) then
            if (BCbotT.eq.1) then !neumann/fixed flux bc for temperature
#if defined(_GPU)
               call bottom_update_thlp_cuda<<<griddim,blockdim>>>(wtsurf)
               call checkCUDA( cudaGetLastError(), 'bottom_update_thlp_cuda' )
#else
               do j = jb, je
                  do i = ib, ie
                     thlp(i, j, kb) = thlp(i, j, kb) &
                                      + ( &
                                      0.5*(dzf(kb - 1)*ekh(i, j, kb) + dzf(kb)*ekh(i, j, kb - 1)) &
                                      *(thl0(i, j, kb) - thl0(i, j, kb - 1)) &
                                      *dzh2i(kb) &
                                      - wtsurf &
                                      )*dzfi(kb)
                  end do
               end do
#endif
            else if (BCbotT.eq.2) then !wall function bc for temperature (fixed temperature)
#if defined(_GPU)
               bcTfluxA_d = 0.
               call wfuno_cuda<<<griddim,blockdim>>>(ih, jh, kh, up_d, vp_d, thlp_d, &
                                                     momfluxb_d, tfluxb_d, bcTfluxA_d, &
                                                     u0_d, v0_d, thl0_d, thls, z0, z0h, 92)
               call checkCUDA( cudaGetLastError(), 'wfuno_cuda under BCbotT' )
               bcTfluxA = bcTfluxA_d
#else
               call wfuno(ih, jh, kh, up, vp, thlp, momfluxb, tfluxb, bcTfluxA, &
                          u0, v0, thl0, thls, z0, z0h, 92)
#endif
            else
               write(0, *) "ERROR: bottom boundary type for temperature undefined"
               stop 1
            end if
         end if ! ltempeq

         if (lmoist) then
            if (BCbotq.eq.1) then !neumann/fixed flux bc for moisture
#if defined(_GPU)
               call bottom_update_qtp_cuda<<<griddim,blockdim>>>(wqsurf)
               call checkCUDA( cudaGetLastError(), 'bottom_update_qtp_cuda' )
#else
               do j = jb, je
                  do i = ib, ie
                     qtp(i, j, kb) = qtp(i, j, kb) + ( &
                                     0.5*(dzf(kb - 1)*ekh(i, j, kb) + dzf(kb)*ekh(i, j, kb - 1)) &
                                     *(qt0(i, j, kb) - qt0(i, j, kb - 1)) &
                                     *dzh2i(kb) &
                                     + wqsurf &
                                     )*dzfi(kb)
                  end do
               end do
#endif
            else
               write(0, *) "ERROR: bottom boundary type for moisture undefined"
               stop 1
            end if
         end if !lmoist

         if (nsv>0) then
            if (BCbots.eq.1) then !neumann/fixed flux bc for moisture
               do m = 1, nsv
#if defined(_GPU)
                  call bottom_update_svp_cuda<<<griddim,blockdim>>>(m)
                  call checkCUDA( cudaGetLastError(), 'bottom_update_svp_cuda' )
#else
                  do j = jb, je
                     do i = ib, ie
                        svp(i, j, kb, m) = svp(i, j, kb, m) + ( &
                                           0.5*(dzf(kb - 1)*ekh(i, j, kb) + dzf(kb)*ekh(i, j, kb - 1)) &
                                           *(sv0(i, j, kb, m) - sv0(i, j, kb - 1, m)) &
                                           *dzh2i(kb) &
                                           + 0. &
                                           )*dzfi(kb)
                     end do
                  end do
#endif
               end do
            else
               write(0, *) "ERROR: bottom boundary type for scalars undefined"
               stop 1
            end if
         end if

      end if

#if defined(_GPU)
      call bottom_update_tau_cuda<<<griddim,blockdim>>>
      call checkCUDA( cudaGetLastError(), 'bottom_update_tau_cuda' )
#else
      tau_x(:,:,kb:ke+kh) = up - tau_x(:,:,kb:ke+kh)
      tau_y(:,:,kb:ke+kh) = vp - tau_y(:,:,kb:ke+kh)
      tau_z(:,:,kb:ke+kh) = wp - tau_z(:,:,kb:ke+kh)
      thl_flux(:,:,kb:ke+kh) = thlp - thl_flux(:,:,kb:ke+kh)
#endif

      return
   end subroutine bottom


   subroutine createmasks
      use modglobal, only : libm, ib, ie, jb, je, kb, ke, khc, jtot, rslabs
      use modfields, only : IIc,  IIu,  IIv,  IIw,  IIuw,  IIvw,  IIuv,  &
                            IIcs, IIus, IIvs, IIws, IIuws, IIvws, IIuvs, &
                            IIct, IIut, IIvt, IIwt, IIuwt
      use modmpi,    only : comm3d, mpierr
      use m_halo,    only : halo_exchange

      integer :: IIcl(kb:ke + khc), IIul(kb:ke + khc), IIvl(kb:ke + khc), IIwl(kb:ke + khc), IIuwl(kb:ke + khc), IIvwl(kb:ke + khc), IIuvl(kb:ke + khc)
      integer :: IIcd(ib:ie, kb:ke)
      integer :: IIwd(ib:ie, kb:ke)
      integer :: IIuwd(ib:ie, kb:ke)
      integer :: IIud(ib:ie, kb:ke)
      integer :: IIvd(ib:ie, kb:ke)
      integer :: i, j, k, n

      ! II*l needn't be defined up to ke_khc, but for now would require large scale changes in modstatsdump so if works leave as is ! tg3315 04/07/18

      if (.not. libm) then
         IIc(:, :, :) = 1
         IIu(:, :, :) = 1
         IIv(:, :, :) = 1
         IIw(:, :, :) = 1
         IIuw(:, :, :) = 1
         IIvw(:, :, :) = 1
         IIuv(:, :, :) = 1
         IIcs(:) = nint(rslabs)
         IIus(:) = nint(rslabs)
         IIvs(:) = nint(rslabs)
         IIws(:) = nint(rslabs)
         IIuws(:) = nint(rslabs)
         IIvws(:) = nint(rslabs)
         IIuvs(:) = nint(rslabs)
         IIct(:, :) = jtot
         IIut(:, :) = jtot
         IIvt(:, :) = jtot
         IIwt(:, :) = jtot
         IIuwt(:, :) = jtot
         return
      end if
      ! Create masking matrices
      IIc = 1; IIu = 1; IIv = 1; IIct = 1; IIw = 1; IIuw = 1; IIvw = 1; IIuv = 1; IIwt = 1; IIut = 1; IIvt = 1; IIuwt = 1; IIcs = 1; IIus = 1; IIvs = 1; IIws = 1; IIuws = 1; IIvws = 1; IIuvs = 1

      do n = 1,solid_info_u%nsolptsrank
       !n = solid_info_u%solptsrank(m)
          i = solid_info_u%solpts_loc(n,1)
          j = solid_info_u%solpts_loc(n,2)
          k = solid_info_u%solpts_loc(n,3)
          IIu(i,j,k) = 0
      end do

      do n = 1,solid_info_v%nsolptsrank
       !n = solid_info_v%solptsrank(m)
          i = solid_info_v%solpts_loc(n,1)
          j = solid_info_v%solpts_loc(n,2)
          k = solid_info_v%solpts_loc(n,3)
          IIv(i,j,k) = 0
      end do

      do n = 1,solid_info_w%nsolptsrank
       !n = solid_info_w%solptsrank(m)
          i = solid_info_w%solpts_loc(n,1)
          j = solid_info_w%solpts_loc(n,2)
          k = solid_info_w%solpts_loc(n,3)
          IIw(i,j,k) = 0
      end do

      do n = 1,solid_info_c%nsolptsrank
       !n = solid_info_c%solptsrank(m)
          i = solid_info_c%solpts_loc(n,1)
          j = solid_info_c%solpts_loc(n,2)
          k = solid_info_c%solpts_loc(n,3)
          IIc(i,j,k) = 0
      end do

      IIw(:, :, kb) = 0; IIuw(:, :, kb) = 0; IIvw(:, :, kb) = 0

      do i=ib,ie
        do j=jb,je
          IIuv(i,j,kb) = IIu(i,j,kb) * IIu(i,j-1,kb) * IIv(i,j,kb) * IIv(i-1,j,kb)
          do k=kb+1,ke
            ! Classed as solid (set to zero) unless ALL points in the stencil are fluid
            IIuv(i,j,k) = IIu(i,j,k) * IIu(i,j-1,k) * IIv(i,j,k) * IIv(i-1,j,k)
            IIuw(i,j,k) = IIu(i,j,k) * IIu(i,j,k-1) * IIw(i,j,k) * IIw(i-1,j,k)
            IIvw(i,j,k) = IIv(i,j,k) * IIv(i,j,k-1) * IIw(i,j,k) * IIw(i,j-1,k)
          end do
        end do
      end do

      ! Can't do this because no interface for integers
      ! call halo_exchange(IIuv, 3, opt_zlevel=(/ihc,jhc,0/))
      ! call halo_exchange(IIuv, 3, opt_zlevel=(/ihc,jhc,0/))
      ! call halo_exchange(IIvw, 3, opt_zlevel=(/ihc,jhc,0/))

      do k = kb, ke + khc
         IIcl(k) = sum(IIc(ib:ie, jb:je, k))
         IIul(k) = sum(IIu(ib:ie, jb:je, k))
         IIvl(k) = sum(IIv(ib:ie, jb:je, k))
         IIwl(k) = sum(IIw(ib:ie, jb:je, k))
         IIuwl(k) = sum(IIuw(ib:ie, jb:je, k))
         IIvwl(k) = sum(IIvw(ib:ie, jb:je, k))
         IIuvl(k) = sum(IIuv(ib:ie, jb:je, k))
      enddo

      call MPI_ALLREDUCE(IIcl, IIcs, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIul, IIus, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIvl, IIvs, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIwl, IIws, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIuwl, IIuws, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIvwl, IIvws, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIuvl, IIuvs, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)

      IIcd(ib:ie, kb:ke) = sum(IIc(ib:ie, jb:je, kb:ke), DIM=2)
      IIwd(ib:ie, kb:ke) = sum(IIw(ib:ie, jb:je, kb:ke), DIM=2)
      IIuwd(ib:ie, kb:ke) = sum(IIuw(ib:ie, jb:je, kb:ke), DIM=2)
      IIud(ib:ie, kb:ke) = sum(IIu(ib:ie, jb:je, kb:ke), DIM=2)
      IIvd(ib:ie, kb:ke) = sum(IIv(ib:ie, jb:je, kb:ke), DIM=2)

      call MPI_ALLREDUCE(IIwd(ib:ie, kb:ke), IIwt(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIcd(ib:ie, kb:ke), IIct(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIuwd(ib:ie, kb:ke), IIuwt(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIud(ib:ie, kb:ke), IIut(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIvd(ib:ie, kb:ke), IIvt(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, comm3d, mpierr)

   end subroutine createmasks

end module modibm
