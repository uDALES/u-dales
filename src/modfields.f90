!> \file modfields.f90
!!  Declares, allocates and initializes the 3D fields

!>
!!  Declares, allocates and initializes the 3D fields
!>


module modfields

  implicit none
  save

  ! Prognostic variables

  real, allocatable :: worksave(:)      !<   Used in POISR!
  real, allocatable :: um(:,:,:)        !<   x-component of velocity at time step t-1
  real, allocatable :: vm(:,:,:)        !<   y-component of velocity at time step t-1
  real, allocatable :: wm(:,:,:)        !<   z-component of velocity at time step t-1
  real, allocatable :: thlm(:,:,:)      !<   liq. water pot. temperature at time step t-1
  real, allocatable :: e12m(:,:,:)      !<   turb. kin. energy at time step t-1
  real, allocatable :: qtm(:,:,:)       !<   total specific humidity at time step t
  real, allocatable, target :: u0(:,:,:)        !<   x-component of velocity at time step t
  real, allocatable, target :: v0(:,:,:)        !<   y-component of velocity at time step t
  real, allocatable, target :: w0(:,:,:)        !<   z-component of velocity at time step t
  real, allocatable, target :: pres0(:,:,:)     !<   pressure at time step t
  real, allocatable, target :: thl0(:,:,:)      !<   liq. water pot. temperature at time step t
  real, allocatable :: thl0h(:,:,:)     !<   3d-field of theta_l at half levels for kappa scheme

  real, allocatable :: qt0h(:,:,:)      !<  3d-field of q_tot   at half levels for kappa scheme
  real, allocatable :: e120(:,:,:)      !<   turb. kin. energy at time step t
  real, allocatable, target :: qt0(:,:,:)       !<   total specific humidity at time step t

  real, allocatable :: up(:,:,:)        !<   tendency of um
  real, allocatable :: vp(:,:,:)        !<   tendency of vm
  real, allocatable :: wp(:,:,:)        !<   tendency of wm
  real, allocatable :: thlp(:,:,:)      !<   tendency of thlm
  real, allocatable :: e12p(:,:,:)      !<   tendency of e12m
  real, allocatable :: qtp(:,:,:)       !<   tendency of qtm

  real, allocatable :: svm(:,:,:,:)     !<  scalar sv(n) at time step t-1
  real, allocatable, target :: sv0(:,:,:,:)     !<  scalar sv(n) at time step t
  real, allocatable :: svp(:,:,:,:)     !<  tendency of sv(n)
  real, allocatable :: svpp(:,:,:,:)
 
  ! Diagnostic variables
  real, allocatable :: mindist(:,:,:)   !< minimal distance of cell center to a wall

  real, allocatable :: shear(:,:,:,:)   !<   wall shear (last rank indicates the type of shear componenten (uym, uyp, etc.)

   real, allocatable :: momfluxb(:,:,:) !< fields for the wallfluxes of total momentum  
   real, allocatable :: tfluxb(:,:,:)  !< heat
   real, allocatable :: qfluxb(:,:,:)  !< and moisture
   real, allocatable :: cth(:,:,:)     !< heat transfer coefficient

  !tg3315 added variables (statistics, masking and others) 
  integer, allocatable :: IIc(:,:,:)        !< Masking matrix for blocks at cell centres
  integer, allocatable :: IIu(:,:,:)        !< Masking matrix for blocks at x-direction half cells
  integer, allocatable :: IIv(:,:,:)        !< Masking matrix for blocks at y-direction half cells
  integer, allocatable :: IIw(:,:,:)        !< Masking matrix for blocks at z-direction half cells
  integer, allocatable :: IIuw(:,:,:)       !< Masking matrix for blocks at x-and z-direction half cells
  integer, allocatable :: IIvw(:,:,:)       !< Masking matrix for blocks at y- and z-direction half cells
  integer, allocatable :: IIuv(:,:,:)       !< Masking matrix for blocks at x- and y-direction half cells  
  integer, allocatable :: IIct(:,:)         !< 2-D Masking matrix for blocks at cell centre that span 1:jtot
  integer, allocatable :: IIwt(:,:)         !< 2-D Masking matrix for blocks at z-direction half cells that span 1:jtot
  integer, allocatable :: IIuwt(:,:)        !< 2-D Masking matrix for blocks at x- and z-direction half cells that span 1:jtot
  integer, allocatable :: IIut(:,:)         !< 2-D Masking matrix for blocks at x-direction half cells that span 1:jtot
  integer, allocatable :: IIvt(:,:)         !< 2-D Masking matrix for blocks at y-direction half cells that span 1:jtot
  integer, allocatable :: IIcs(:)           !< 1-D Masking matrix for blocks at cell centres that span ib:ie and 1:jtot
  integer, allocatable :: IIus(:)           !< 1-D Masking matrix for blocks at x-direction half cells that span ib:ie and 1:jtot
  integer, allocatable :: IIvs(:)           !< 1-D Masking matrix for blocks at y-direction half cells that span ib:ie and 1:jtot
  integer, allocatable :: IIws(:)           !< 1-D Masking matrix for blocks at z-direction half cells that span ib:ie and 1:jtot
  integer, allocatable :: IIuws(:)          !< 1-D Masking matrix for blocks at x- and z-direction half cells that span ib:ie and 1:jtot
  integer, allocatable :: IIvws(:)          !< 1-D Masking matrix for blocks at y- and z-direction half cells that span ib:ie and 1:jtot
  integer, allocatable :: IIuvs(:)          !< 1-D Masking matrix for blocks at x- and y-direction half cells that span ib:ie and 1:jtot

!  integer              :: IIbl = 1          !< Switch for if layer at kb is all blocks

  ! statistical fields following notation "[statistical name][averaging directions - x,y,z,t][position in grid - i,j,k]"
  real, allocatable :: uyt(:,:) 
  real, allocatable :: uytik(:,:)
  real, allocatable :: vyt(:,:)        
  real, allocatable :: wyt(:,:)        
  real, allocatable :: wytik(:,:)
  real, allocatable :: thlyt(:,:)      
  real, allocatable :: qtyt(:,:)      
  real, allocatable :: thlytk(:,:)      
  real, allocatable :: sca1yt(:,:)     
  real, allocatable :: sca2yt(:,:)     
  real, allocatable :: sca3yt(:,:)     
  real, allocatable :: thlsgsyt(:,:)     
  real, allocatable :: qtsgsyt(:,:)
  real, allocatable :: usgsyt(:,:) 
  real, allocatable :: wsgsyt(:,:)    
  real, allocatable :: sv1sgsyt(:,:)
  real, allocatable :: sv2sgsyt(:,:)
  real, allocatable :: sv3sgsyt(:,:)

  real, allocatable :: uxyt(:)        
  real, allocatable :: vxyt(:)        
  real, allocatable :: wxyt(:)        
  real, allocatable :: thlxyt(:)        
  real, allocatable :: qtxyt(:)
  real, allocatable :: pxyt(:)  ! bss116      
  real, allocatable :: usgsxyt(:)
  real, allocatable :: thlsgsxyt(:)
  real, allocatable :: vsgsxyt(:)

  real, allocatable :: uwtik(:,:,:)
  real, allocatable :: wthltk(:,:,:)
  real, allocatable :: wqttk(:,:,:)
  real, allocatable :: thlthlt(:,:,:)
  real, allocatable :: qtqtt(:,:,:)
  real, allocatable :: sv1sv1t(:,:,:)
  real, allocatable :: sv2sv2t(:,:,:)
  real, allocatable :: sv3sv3t(:,:,:)
  real, allocatable :: sv4sv4t(:,:,:)
  ! real, allocatable :: sv1max(:,:,:)
  ! real, allocatable :: sv2max(:,:,:)
  ! real, allocatable :: sv3max(:,:,:)
  ! real, allocatable :: sv4max(:,:,:)
  real, allocatable :: PSSt(:,:,:)
  real, allocatable :: uutc(:,:,:)
  real, allocatable :: vvtc(:,:,:)
  real, allocatable :: wwtc(:,:,:)
  real, allocatable :: vwtjk(:,:,:)
  real, allocatable :: uvtij(:,:,:)
  real, allocatable :: utik(:,:,:)
  real, allocatable :: wtik(:,:,:)
  real, allocatable :: vtjk(:,:,:)
  real, allocatable :: wtjk(:,:,:)
  real, allocatable :: utij(:,:,:)
  real, allocatable :: vtij(:,:,:)
  real, allocatable :: wmt(:,:,:)
  real, allocatable :: thltk(:,:,:)
  real, allocatable :: qttk(:,:,:)
  real, allocatable :: thlt(:,:,:)
  real, allocatable :: utc(:,:,:)
  real, allocatable :: vtc(:,:,:)
  real, allocatable :: wtc(:,:,:)

  real, allocatable :: vmt(:,:,:)
  real, allocatable :: umt(:,:,:)
  real, allocatable :: sv1t(:,:,:)
  real, allocatable :: sv2t(:,:,:)
  real, allocatable :: sv3t(:,:,:)
  real, allocatable :: sv4t(:,:,:)
  real, allocatable :: sv1tk(:,:,:)
  real, allocatable :: sv2tk(:,:,:)
  real, allocatable :: sv3tk(:,:,:)
  real, allocatable :: sv4tk(:,:,:)
  real, allocatable :: wsv1tk(:,:,:)
  real, allocatable :: wsv2tk(:,:,:)
  real, allocatable :: wsv3tk(:,:,:) 
  real, allocatable :: wsv4tk(:,:,:)
  real, allocatable :: sv1sgst(:,:,:)
  real, allocatable :: sv2sgst(:,:,:)
  real, allocatable :: sv3sgst(:,:,:)
  real, allocatable :: sv4sgst(:,:,:)
  real, allocatable :: qtt(:,:,:) ! bss116
  real, allocatable :: pt(:,:,:)  ! bss116

  real, allocatable :: slice(:,:) 
  real, allocatable :: slice2(:,:) 
  real, allocatable :: slice3(:,:) 
  real, allocatable :: slice4(:,:) 
  real, allocatable :: slice5(:,:) 
  real, allocatable :: slice6(:,:) 
  real, allocatable :: slice7(:,:) 
  real, allocatable :: slice8(:,:)

  ! fields for scalar sources
  real, allocatable :: scar(:,:)
  real, allocatable :: scarl(:,:)

  real, allocatable :: uav(:,:,:)       !<   time-averaged u-velocity
  real, allocatable :: vav(:,:,:)       !<   time-averaged u-velocity
  real, allocatable :: wav(:,:,:)       !<   time-averaged u-velocity 
  real, allocatable :: thlav(:,:,:)     !<   time-averaged liquid temperature
  real, allocatable :: qtav(:,:,:)      !<   time-averaged specific humidity
  real, allocatable :: qlav(:,:,:)      !<   time-averaged liquid water
  real, allocatable :: presav(:,:,:)    !<   time-averaged pressure
  real, allocatable :: svav(:,:,:,:)    !<   time-averaged scalar concentration 
  real, allocatable :: viscratioav(:,:,:)    !<   time-averaged viscosity ratio; turb viscosity / molecular viscosity
  real, allocatable :: umint(:,:,:)     !<   um interpolated to cell-center
  real, allocatable :: vmint(:,:,:)     !<   vm interpolated to cell-center
  real, allocatable :: wmint(:,:,:)     !<   wm interpolated to cell-center
  real, allocatable :: thl2av(:,:,:)    !<   time-average: liquid temperature squared
  real, allocatable :: ql2av(:,:,:)    !<   time-average: liquid temperature squared
  real, allocatable :: qt2av(:,:,:)    !<   time-average: liquid temperature squared
  real, allocatable :: sv2av(:,:,:,:)   !<   time-average: scalar concentration squared
  real, allocatable :: uuav(:,:,:)      !<   time-average: u-velocity squared
  real, allocatable :: vvav(:,:,:)      !<   time-average: v-velocity squared
  real, allocatable :: wwav(:,:,:)      !<   time-average: w-velocity squared
  real, allocatable :: uvav(:,:,:)      !<   time-average: u-velocity times v-velocity
  real, allocatable :: uwav(:,:,:)      !<   time-average: u-velocity times fluctuation
  real, allocatable :: vwav(:,:,:)      !<   time-average: v-velocity times w-velocity
  real, allocatable :: thluav(:,:,:)    !<   time-average: thl times u-velocity
  real, allocatable :: thlvav(:,:,:)    !<   time-average: thl times v-velocity
  real, allocatable :: thlwav(:,:,:)    !<   time-average: thl times w-velocity
  real, allocatable :: thlthlav(:,:,:)  !<   time-average: thl times thl
  real, allocatable :: qluav(:,:,:)    !<   time-average: ql times u-velocity
  real, allocatable :: qlvav(:,:,:)    !<   time-average: ql times v-velocity
  real, allocatable :: qlwav(:,:,:)    !<   time-average: ql times w-velocity
  real, allocatable :: qtuav(:,:,:)    !<   time-average: qt times u-velocity
  real, allocatable :: qtvav(:,:,:)    !<   time-average: qt times v-velocity
  real, allocatable :: qtwav(:,:,:)    !<   time-average: qt times w-velocity
  real, allocatable :: svuav(:,:,:,:)   !<   time-average: sv times u-velocity
  real, allocatable :: svvav(:,:,:,:)   !<   time-average: sv times v-velocity
  real, allocatable :: svwav(:,:,:,:)   !<   time-average: sv times w-velocity
!  real, allocatable :: tekm(:,:,:)     !tekm = ekm - numol !tg3315

  real, allocatable :: upupav(:,:,:)    !<   time-average: u'u'
  real, allocatable :: vpvpav(:,:,:)    !<   time-average: v'v'
  real, allocatable :: wpwpav(:,:,:)    !<   time-average: w'w'
  real, allocatable :: thlpthlpav(:,:,:)!<   time-average: thl'thl'
  real, allocatable :: qlpqlpav(:,:,:)  !<   time-average: ql'ql'
  real, allocatable :: qtpqtpav(:,:,:)!<   time-average: thl'thl'
  real, allocatable :: svpsvpav(:,:,:,:)!<   time-average: sv'sv'
  real, allocatable :: upvpav(:,:,:)    !<   time-average: u'v'
  real, allocatable :: upwpav(:,:,:)    !<   time-average: u'w'
  real, allocatable :: vpwpav(:,:,:)    !<   time-average: v'w'
  real, allocatable :: thlpupav(:,:,:)  !<   time-average: thl'u'
  real, allocatable :: thlpvpav(:,:,:)  !<   time-average: thl'v'
  real, allocatable :: thlpwpav(:,:,:)  !<   time-average: thl'w'
  real, allocatable :: qlpupav(:,:,:)  !<   time-average: ql'u'
  real, allocatable :: qlpvpav(:,:,:)  !<   time-average: ql'v'
  real, allocatable :: qlpwpav(:,:,:)  !<   time-average: ql'w'
  real, allocatable :: qtpupav(:,:,:)  !<   time-average: qt'u'
  real, allocatable :: qtpvpav(:,:,:)  !<   time-average: qt'v'
  real, allocatable :: qtpwpav(:,:,:)  !<   time-average: qt'w'
  real, allocatable :: svpupav(:,:,:,:) !<   time-average: sv'u'
  real, allocatable :: svpvpav(:,:,:,:) !<   time-average: sv'v'
  real, allocatable :: svpwpav(:,:,:,:) !<   time-average: sv'w'

! SGS fields
  real, allocatable :: uusgsav(:,:,:)    !<   time-average subgrid contribution (estimate)
  real, allocatable :: vvsgsav(:,:,:)    !<   time-average subgrid contribution (estimate)
  real, allocatable :: wwsgsav(:,:,:)    !<   time-average subgrid contribution (estimate)
  real, allocatable :: uwsgsav(:,:,:)    !<   time-average subgrid contribution (estimate)
  real, allocatable :: thlusgsav(:,:,:)  !<   time-average subgrid contribution (estimate)
  real, allocatable :: thlwsgsav(:,:,:)  !<   time-average subgrid contribution (estimate)
  real, allocatable :: qlusgsav(:,:,:)  !<   time-average subgrid contribution (estimate)
  real, allocatable :: qlwsgsav(:,:,:)  !<   time-average subgrid contribution (estimate)
  real, allocatable :: qtusgsav(:,:,:)  !<   time-average subgrid contribution (estimate)
  real, allocatable :: qtwsgsav(:,:,:)  !<   time-average subgrid contribution (estimate)
  real, allocatable :: svusgsav(:,:,:,:) !<   time-average subgrid contribution (estimate)
  real, allocatable :: svwsgsav(:,:,:,:) !<   time-average subgrid contribution (estimate)
  real, allocatable :: tkesgsav(:,:,:)   !<   time-average subgrid turbulence kinetic energy
  real, allocatable :: nusgsav(:,:,:)    !<   time-average subgrid viscosity

! Resolved dissipation 'terms'
  real, allocatable :: strain2av(:,:,:)  !<   <Sij*Sij> used to compute <Sij'*Sij'> = <Sij*Sij> - <S>ij*<S>ij
  real, allocatable :: disssgsav(:,:,:)  !<   mean subgrid dissipation: <nu_sgs*2.*Sij*Sij>
                                         !<   which is used for resolved dissipation = nu*2*<Sij'*Sij'> 
! TKE budget terms:
  real, allocatable :: tvmx(:,:,:)        !<   needed for viscous transport: <u*d/dxj(2*nu*S1j)>
  real, allocatable :: tvmy(:,:,:)        !<   needed for viscous transport: <v*d/dxj(2*nu*S2j)>
  real, allocatable :: tvmz(:,:,:)        !<   needed for viscous transport: <w*d/dxj(2*nu*S3j)>
  real, allocatable :: tpm (:,:,:)        !<   needed for transport by pressure fluctuations
  real, allocatable :: ttmx(:,:,:)        !<   needed for transport by turb. vel. fluctuations
  real, allocatable :: ttmy(:,:,:)        !<   needed for transport by turb. vel. fluctuations
  real, allocatable :: ttmz(:,:,:)        !<   needed for transport by turb. vel. fluctuations
  real, allocatable :: tsgsmx1(:,:,:)     !<   needed for transport by subgrid x = <u*d/dxj(2*nu_t*S1j)>
  real, allocatable :: tsgsmy1(:,:,:)     !<   needed for transport by subgrid y = <v*d/dxj(2*nu_t*S2j)>
  real, allocatable :: tsgsmz1(:,:,:)     !<   needed for transport by subgrid z = <w*d/dxj(2*nu_t*S3j)>
  real, allocatable :: tsgsmx2(:,:,:)     !<   needed for transport by subgrid x = <d/dxj(2*nu_t*S1j)>
  real, allocatable :: tsgsmy2(:,:,:)     !<   needed for transport by subgrid y = <d/dxj(2*nu_t*S2j)>
  real, allocatable :: tsgsmz2(:,:,:)     !<   needed for transport by subgrid z = <d/dxj(2*nu_t*S3j)>
! TKE budget results (written to files):
  real, allocatable :: t_vav(:,:,:)        !<   viscous transport
  real, allocatable :: t_sgsav(:,:,:)      !<   transport by subgrid
  real, allocatable :: t_pav(:,:,:)        !<   transport by pressure fluctuations
  real, allocatable :: t_tav(:,:,:)        !<   transport by by turb. vel. fluctuations
  real, allocatable :: p_tav(:,:,:)        !<   production by shear
  real, allocatable :: p_bav(:,:,:)        !<   production/destruction by buoyancy
  real, allocatable :: d_sgsav(:,:,:)      !<   dissipation by subgrid
  real, allocatable :: tkeadv(:,:,:)       !<   advection of tke 
  
! TKE budget results (written to tkedump): !tg3315
  real, allocatable :: t_v(:)            !<   viscous transport
  real, allocatable :: t_sgs(:)          !<   transport by subgrid
  real, allocatable :: t_p(:)            !<   transport by pressurefluctuations
  real, allocatable :: t_t(:)            !<   transport by by turb. vel. fluctuations
  real, allocatable :: p_t(:)            !<   production by shear
  real, allocatable :: p_b(:)            !<   production/destruction by buoyancy
  real, allocatable :: d_sgs(:)          !<   dissipation by subgrid
  real, allocatable :: adv(:)            !<   advection of tke

  real, allocatable, target :: ql0(:,:,:)       !<   liquid water content

  real, allocatable :: thv0h(:,:,:)     !<   theta_v at half level

  real, allocatable :: whls(:)          !<   large scale vert velocity at half levels

  real, allocatable :: presf(:)         !<   hydrostatic pressure at full level
  real, allocatable :: presh(:)         !<   hydrostatic pressure at half level
  real, allocatable :: exnf(:)          !<   hydrostatic exner function at full level
  real, allocatable :: exnh(:)          !<   hydrostatic exner function at half level
  real, allocatable :: thvf(:)          !<   hydrostatic exner function at full level
  real, allocatable :: thvh(:)          !<   hydrostatic exner function at half level
  real, allocatable :: rhof(:)          !<   slab averaged density at full level
  real, allocatable :: qt0av(:)         !<   slab averaged q_tot
  real, allocatable :: ql0av(:)         !<   slab averaged q_liq

  real, allocatable :: thl0av(:)        !<   slab averaged th_liq
  real, allocatable :: u0av(:)          !<   slab averaged u
  real, allocatable :: v0av(:)          !<   slab averaged v
  real, allocatable :: ug(:)            !<   geostrophic u-wind
  real, allocatable :: vg(:)            !<   geostrophic v-wind

  real, allocatable :: pgx(:)            !<   driving pressure gradient in x, this is dp/dx [(\Delta p) / (\Delta x)] across one cell, already divided by \rho -> in units of [m/s^2]
  real, allocatable :: pgy(:)            !<   driving pressure gradient in y [m/s^2]

  real, allocatable :: dpdxl(:)                      !<   large scale pressure x-gradient [m/s^2]
  real, allocatable :: dpdyl(:)                      !<   large scale pressure y-gradient [m/s^2]

  real, allocatable :: dthldxls(:)                   !<   large scale x-gradient of th_liq
  real, allocatable :: dthldyls(:)                   !<   large scale y-gradient of th_liq
  real, allocatable :: dqtdxls(:)                    !<   large scale x-gradient of q_tot
  real, allocatable :: dqtdyls(:)                    !<   large scale y-gradient of q_tot
  real, allocatable :: dqtdtls(:)                    !<   large scale y-gradient of q_tot
  real, allocatable :: dudxls(:)                     !<   large scale x-gradient of u

  real, allocatable :: dudyls(:)                     !<   large scale y-gradient of u
  real, allocatable :: dvdxls(:)                     !<   large scale x-gradient of v
  real, allocatable :: dvdyls(:)                     !<   large scale y-gradient of v
  real, allocatable :: wfls  (:)                     !<   large scale y-gradient of v
  real, allocatable :: ql0h(:,:,:)
  real, allocatable :: dthvdz(:,:,:)                 !<   theta_v at half level

  real, allocatable :: thlprof(:)                    !<   initial thl-profile
  real, allocatable :: qtprof(:)                     !<   initial qt-profile
  real, allocatable :: uprof(:)                      !<   initial u-profile
  real, allocatable :: vprof(:)                      !<   initial v-profile
  real, allocatable :: e12prof(:)                    !<   initial subgrid TKE profile
  real, allocatable :: sv0av(:,:)                    !<   slab average of sv(n)
  real, allocatable :: svprof(:,:)                   !<   initial sv(n)-profile
  real, allocatable :: qlprof(:)

  real, allocatable :: thlpcar(:)                    !< prescribed radiatively forced thl tendency
  real, allocatable :: SW_up_TOA(:,:), SW_dn_TOA(:,:), LW_up_TOA(:,:), LW_dn_TOA(:,:)
  real, allocatable :: uout(:)                      !< height average outlet velocity (used in convective outflow BC)
  real, allocatable :: wout(:)                      !< j-averaged top velocity
  real, allocatable :: friction(:)                  !< skin-friction coeff: from y-line-averaged shear 
  real, allocatable :: momthick(:)                  !< momentum thickness: y-line average 
  real, allocatable :: displthick(:)                !< displacement thickness: y-line average 
  real              :: uouttot                      !< area-averaged outflow velocity (used in convective outflow BC) 
  real              :: wouttot                      !< area-averaveraged top velocity
  real              :: udef
  real              :: vdef
  real, allocatable :: vout(:)
  real              :: vouttot

  real              :: thlsrcdt                     ! thlsrc -> thlsrcdt is used to solve 1-order ODE for thlsrc 
  real              :: dgdt                         ! g = dp/dx -> dgdt is used to solve 1-order ODE for dpdx 
  real              :: dpdx = 0.                   ! dpdx given in namoptions

  real              :: uoutarea                     !< area of domain u-outlet
  real              :: voutarea                     !< area of domain v-outlet
  real              :: fluidvol                     !< fluid volume (excluding blocks)

  character(80), allocatable :: ncname(:,:)
  character(80), allocatable :: ncstaty(:,:)
  character(80), allocatable :: ncstatyt(:,:)
  character(80), allocatable :: ncstattke(:,:)
  character(80), allocatable :: ncstatxy(:,:)
  character(80), allocatable :: ncstatxyt(:,:)
  character(80), allocatable :: ncstatslice(:,:)
  character(80), allocatable :: ncstatt(:,:)

  integer, allocatable :: wall(:,:,:,:)             !< wall(ic,jc,kc,1-5) gives the global indices of the wall closest to cell center ic,jc,kc. The 4th and 5th integer gives the corresponding shear components

contains
  !> Allocate and initialize the prognostic variables
  subroutine initfields

    use modglobal, only : ib,ie,jb,je,ih,jh,kb,ke,kh,nsv,jtot,imax,jmax,kmax,&
         ihc,jhc,khc!, iadv_kappa,iadv_sv
    ! Allocation of prognostic variables
    implicit none

    allocate(worksave(2*imax*jmax*kmax))
    allocate(um(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(vm(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(wm(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(thlm(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(e12m(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(qtm(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(u0(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(v0(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(w0(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(pres0(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(thl0(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(thl0h(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(qt0h(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(e120(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(qt0(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(ql0(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(up(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(vp(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(wp(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(thlp(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(e12p(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(qtp(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(svm(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb-khc:ke+khc,nsv))
    allocate(sv0(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb-khc:ke+khc,nsv))
    allocate(svp(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc,nsv))
    allocate(svpp(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc,1)) 

    ! Allocation of diagnostic variables
    allocate(mindist(ib:ie,jb:je,kb:ke))
    allocate(thv0h(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(whls(kb:ke+kh))
    allocate(presf(kb:ke+kh))
    allocate(presh(kb:ke+kh))
    allocate(exnf(kb:ke+kh))
    allocate(exnh(kb:ke+kh))
    allocate(thvf(kb:ke+kh))
    allocate(thvh(kb:ke+kh))
    allocate(rhof(kb:ke+kh))
    allocate(qt0av(kb:ke+kh))
    allocate(ql0av(kb:ke+kh))
    allocate(thl0av(kb:ke+kh))
    allocate(u0av(kb:ke+kh))
    allocate(v0av(kb:ke+kh))
    allocate(ug(kb:ke+kh))
    allocate(vg(kb:ke+kh))
    allocate(pgx(kb:ke+kh))
    allocate(pgy(kb:ke+kh))
    allocate(dpdxl(kb:ke+kh))
    allocate(dpdyl(kb:ke+kh))
    allocate(dthldxls(kb:ke+kh))
    allocate(dthldyls(kb:ke+kh))
    allocate(dqtdxls(kb:ke+kh))
    allocate(dqtdyls(kb:ke+kh))
    allocate(dqtdtls(kb:ke+kh))
    allocate(dudxls(kb:ke+kh))
    allocate(dudyls(kb:ke+kh))
    allocate(dvdxls(kb:ke+kh))
    allocate(dvdyls(kb:ke+kh))
    allocate(wfls  (kb:ke+kh))
    allocate(ql0h(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(dthvdz(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(thlprof(kb:ke+kh))
    allocate(qtprof(kb:ke+kh))
    allocate(qlprof(kb:ke+kh))
    allocate(uprof(kb:ke+kh))
    allocate(vprof(kb:ke+kh))
    allocate(e12prof(kb:ke+kh))
    allocate(sv0av(kb:ke+khc,nsv))
    allocate(svprof(kb:ke+kh,nsv))
    allocate(thlpcar(kb:ke+kh))
    allocate(uout(kb:ke))         ! height average outlet velocity (used in convective outflow BC)
    allocate(vout(kb:ke))
    allocate(wout(ib:ie))         ! j -averaged top velocity
    allocate(friction(ib:ie))     ! line-averaged (along j) skin friction 
    allocate(momthick(ib:ie))     ! line-averaged (along j) momentum thickness 
    allocate(displthick(ib:ie))   ! line-averaged (along j) displacement thickness 
    allocate(SW_up_TOA(ib-ih:ie+ih,jb-jh:je+jh))
    allocate(SW_dn_TOA(ib-ih:ie+ih,jb-jh:je+jh))
    allocate(LW_up_TOA(ib-ih:ie+ih,jb-jh:je+jh))
    allocate(LW_dn_TOA(ib-ih:ie+ih,jb-jh:je+jh))

    ! allocate averaged variables
    allocate(uav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(vav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(wav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(thlav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(qtav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(qlav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(presav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(svav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh,1:nsv))
    allocate(viscratioav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))

    allocate(IIc(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
    allocate(IIu(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
    allocate(IIv(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
    allocate(IIw(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
    allocate(IIuw(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
    allocate(IIvw(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
    allocate(IIuv(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc))
    allocate(IIct(ib:ie,kb:ke))
    allocate(IIwt(ib:ie,kb:ke))
    allocate(IIuwt(ib:ie,kb:ke))
    allocate(IIut(ib:ie,kb:ke))
    allocate(IIvt(ib:ie,kb:ke))
    allocate(IIcs(kb:ke+khc))
    allocate(IIus(kb:ke+khc))
    allocate(IIvs(kb:ke+khc))
    allocate(IIws(kb:ke+khc))
    allocate(IIuws(kb:ke+khc))
    allocate(IIvws(kb:ke+khc))
    allocate(IIuvs(kb:ke+khc))

    allocate(uyt(ib:ie,kb:ke))
    allocate(uytik(ib:ie,kb:ke))
    allocate(vyt(ib:ie,kb:ke))
    allocate(wyt(ib:ie,kb:ke))
    allocate(wytik(ib:ie,kb:ke))
    allocate(thlyt(ib:ie,kb:ke))
    allocate(qtyt(ib:ie,kb:ke))
    allocate(thlytk(ib:ie,kb:ke))
    allocate(sca1yt(ib:ie,kb:ke))
    allocate(sca2yt(ib:ie,kb:ke))
    allocate(sca3yt(ib:ie,kb:ke))
    allocate(usgsyt(ib:ie,kb:ke))
    allocate(thlsgsyt(ib:ie,kb:ke))
    allocate(qtsgsyt(ib:ie,kb:ke))
    allocate(wsgsyt(ib:ie,kb:ke))
    allocate(sv1sgsyt(ib:ie,kb:ke))
    allocate(sv2sgsyt(ib:ie,kb:ke))
    allocate(sv3sgsyt(ib:ie,kb:ke))

    allocate(uxyt(kb:ke+kh))
    allocate(vxyt(kb:ke+kh))
    allocate(wxyt(kb:ke+kh))
    allocate(thlxyt(kb:ke+kh))
    allocate(qtxyt(kb:ke+kh))
    allocate(pxyt(kb:ke+kh))
    allocate(usgsxyt(kb:ke+kh))
    allocate(thlsgsxyt(kb:ke+kh))
    allocate(vsgsxyt(kb:ke+kh))

    allocate(uwtik(ib:ie,jb:je,kb:ke+kh))
    allocate(wthltk(ib:ie,jb:je,kb:ke+kh))
    allocate(wqttk(ib:ie,jb:je,kb:ke+kh))
    allocate(thlthlt(ib:ie,jb:je,kb:ke+kh))
    allocate(qtqtt(ib:ie,jb:je,kb:ke+kh))
    allocate(sv1sv1t(ib:ie,jb:je,kb:ke+kh))
    allocate(sv2sv2t(ib:ie,jb:je,kb:ke+kh))
    allocate(sv3sv3t(ib:ie,jb:je,kb:ke+kh))
    allocate(sv4sv4t(ib:ie,jb:je,kb:ke+kh))
    ! allocate(sv1max(ib:ie,jb:je,kb:ke+kh))
    ! allocate(sv2max(ib:ie,jb:je,kb:ke+kh))
    ! allocate(sv3max(ib:ie,jb:je,kb:ke+kh))
    ! allocate(sv4max(ib:ie,jb:je,kb:ke+kh))
    allocate(PSSt(ib:ie,jb:je,kb:ke+kh))

    allocate(uutc(ib:ie,jb:je,kb:ke+kh))
    allocate(vvtc(ib:ie,jb:je,kb:ke+kh))
    allocate(wwtc(ib:ie,jb:je,kb:ke+kh))
    allocate(vwtjk(ib:ie,jb:je,kb:ke+kh))
    allocate(uvtij(ib:ie,jb:je,kb:ke+kh))
    allocate(utik(ib:ie,jb:je,kb:ke+kh))
    allocate(wtik(ib:ie,jb:je,kb:ke+kh))
    allocate(vtjk(ib:ie,jb:je,kb:ke+kh))
    allocate(wtjk(ib:ie,jb:je,kb:ke+kh))
    allocate(utij(ib:ie,jb:je,kb:ke+kh))
    allocate(vtij(ib:ie,jb:je,kb:ke+kh))
    allocate(wmt(ib:ie,jb:je,kb:ke+kh))
    allocate(thltk(ib:ie,jb:je,kb:ke+kh))
    allocate(qttk(ib:ie,jb:je,kb:ke+kh))
    allocate(thlt(ib:ie,jb:je,kb:ke+kh))
    allocate(utc(ib:ie,jb:je,kb:ke+kh))
    allocate(vtc(ib:ie,jb:je,kb:ke+kh))
    allocate(wtc(ib:ie,jb:je,kb:ke+kh))

    allocate(umt(ib:ie,jb:je,kb:ke+kh))
    allocate(vmt(ib:ie,jb:je,kb:ke+kh))
    allocate(sv1t(ib:ie,jb:je,kb:ke+kh))
    allocate(sv2t(ib:ie,jb:je,kb:ke+kh))
    allocate(sv3t(ib:ie,jb:je,kb:ke+kh))
    allocate(sv4t(ib:ie,jb:je,kb:ke+kh))
    allocate(sv1tk(ib:ie,jb:je,kb:ke+kh))
    allocate(sv2tk(ib:ie,jb:je,kb:ke+kh))
    allocate(sv3tk(ib:ie,jb:je,kb:ke+kh))
    allocate(sv4tk(ib:ie,jb:je,kb:ke+kh))
    allocate(wsv1tk(ib:ie,jb:je,kb:ke+kh))
    allocate(wsv2tk(ib:ie,jb:je,kb:ke+kh))
    allocate(wsv3tk(ib:ie,jb:je,kb:ke+kh))
    allocate(wsv4tk(ib:ie,jb:je,kb:ke+kh))
    allocate(sv1sgst(ib:ie,jb:je,kb:ke+kh))
    allocate(sv2sgst(ib:ie,jb:je,kb:ke+kh))
    allocate(sv3sgst(ib:ie,jb:je,kb:ke+kh))
    allocate(sv4sgst(ib:ie,jb:je,kb:ke+kh))
    allocate(qtt(ib:ie,jb:je,kb:ke+kh))
    allocate(pt(ib:ie,jb:je,kb:ke+kh))

    allocate(slice(ib:ie,jb:je))
    allocate(slice2(ib:ie,jb:je))
    allocate(slice3(ib:ie,jb:je))  
    allocate(slice4(ib:ie,jb:je))
    allocate(slice5(ib:ie,jb:je)) 
    allocate(slice6(ib:ie,jb:je))  
    allocate(slice7(ib:ie,jb:je))
    allocate(slice8(ib:ie,jb:je))

    allocate(scar(ib:ie,jb:jtot))
    allocate(scarl(ib:ie,jb:je))

    allocate(thl2av(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(ql2av(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(qt2av(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(sv2av(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh,1:nsv))
    allocate(umint(ib:ie,jb:je,kb:ke))
    allocate(vmint(ib:ie,jb:je,kb:ke))
    allocate(wmint(ib:ie,jb:je,kb:ke))
    allocate(uuav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(vvav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(wwav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(uvav(ib:ie+ih,jb:je+jh,kb:ke   ))
    allocate(uwav(ib:ie+ih,jb:je   ,kb:ke+kh))
    allocate(vwav(ib:ie   ,jb:je+jh,kb:ke+kh))

    allocate(thluav(ib:ie+ih,jb:je   ,kb:ke   ))
    allocate(thlvav(ib:ie   ,jb:je+jh,kb:ke   ))
    allocate(thlwav(ib:ie   ,jb:je,   kb:ke+kh))
    allocate(thlthlav(ib:ie   ,jb:je,   kb:ke+kh))
    allocate(qluav(ib:ie+ih,jb:je   ,kb:ke   ))
    allocate(qlvav(ib:ie   ,jb:je+jh,kb:ke   ))
    allocate(qlwav(ib:ie   ,jb:je,   kb:ke+kh))
    allocate(qtuav(ib:ie+ih,jb:je   ,kb:ke   ))
    allocate(qtvav(ib:ie   ,jb:je+jh,kb:ke   ))
    allocate(qtwav(ib:ie   ,jb:je,   kb:ke+kh))
    allocate(svuav (ib:ie+ih,jb:je   ,kb:ke   ,1:nsv))
    allocate(svvav (ib:ie   ,jb:je+jh,kb:ke   ,1:nsv))
    allocate(svwav (ib:ie   ,jb:je   ,kb:ke+kh,1:nsv))

    ! <x'x> ( = <xx> -<x><x> )
    allocate(upupav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(vpvpav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(wpwpav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(upvpav(ib:ie+ih,jb:je+jh,kb:ke   ))
    allocate(upwpav(ib:ie+ih,jb:je   ,kb:ke+kh))
    allocate(vpwpav(ib:ie   ,jb:je+jh,kb:ke+kh))

    allocate(thlpthlpav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(thlpupav(ib:ie+ih,jb:je   ,kb:ke   ))
    allocate(thlpvpav(ib:ie   ,jb:je+jh,kb:ke   ))
    allocate(thlpwpav(ib:ie   ,jb:je   ,kb:ke+kh))
    allocate(qlpqlpav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(qlpupav(ib:ie+ih,jb:je   ,kb:ke   ))
    allocate(qlpvpav(ib:ie   ,jb:je+jh,kb:ke   ))
    allocate(qlpwpav(ib:ie   ,jb:je   ,kb:ke+kh))
    allocate(qtpqtpav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(qtpupav(ib:ie+ih,jb:je   ,kb:ke   ))
    allocate(qtpvpav(ib:ie   ,jb:je+jh,kb:ke   ))
    allocate(qtpwpav(ib:ie   ,jb:je   ,kb:ke+kh))
    allocate(svpsvpav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh,1:nsv))
    allocate(svpupav(ib:ie+ih,jb:je   ,kb:ke   ,1:nsv))
    allocate(svpvpav(ib:ie   ,jb:je+jh,kb:ke   ,1:nsv))
    allocate(svpwpav(ib:ie   ,jb:je   ,kb:ke+kh,1:nsv))
 
! Subgrid

    allocate(uusgsav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(vvsgsav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(wwsgsav(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(uwsgsav(ib:ie+ih,jb:je   ,kb:ke+kh))
    allocate(thlusgsav(ib:ie+ih,jb:je   ,kb:ke   ))
    allocate(thlwsgsav(ib:ie   ,jb:je,   kb:ke+kh))
    allocate(qlusgsav(ib:ie+ih,jb:je   ,kb:ke   ))
    allocate(qlwsgsav(ib:ie   ,jb:je,   kb:ke+kh))
    allocate(qtusgsav(ib:ie+ih,jb:je   ,kb:ke   ))
    allocate(qtwsgsav(ib:ie   ,jb:je,   kb:ke+kh))
    allocate(tkesgsav (ib:ie   ,jb:je   ,kb:ke   ))
    allocate(svusgsav (ib:ie+ih,jb:je   ,kb:ke   ,1:nsv))
    allocate(svwsgsav (ib:ie   ,jb:je   ,kb:ke+kh,1:nsv))
    allocate(nusgsav  (ib:ie   ,jb:je   ,kb:ke   ))

! resolved dissipation
    allocate(strain2av(ib:ie,jb:je,kb:ke))
    allocate(disssgsav(ib:ie,jb:je,kb:ke))

! TKE budget terms
    allocate(tvmx  (ib:ie+1,jb:je,kb:ke))
    allocate(tvmy  (ib:ie,jb-1:je+1,kb:ke))
    allocate(tvmz  (ib:ie,jb:je,kb:ke+1))
    allocate(tpm   (ib:ie,jb:je,kb:ke))
    allocate(ttmx  (ib:ie+1,jb:je,    kb:ke))
    allocate(ttmy  (ib:ie,  jb-1:je+1,kb:ke))
    allocate(ttmz  (ib:ie,  jb:je,    kb:ke+1))
    allocate(tsgsmx1(ib:ie+1,jb:je,kb:ke))
    allocate(tsgsmy1(ib:ie,jb-1:je+1,kb:ke))
    allocate(tsgsmz1(ib:ie,jb:je,kb:ke+1))
    allocate(tsgsmx2(ib:ie+1,jb:je,kb:ke))
    allocate(tsgsmy2(ib:ie,jb-1:je+1,kb:ke))
    allocate(tsgsmz2(ib:ie,jb:je,kb:ke+1))
    
    allocate(t_pav  (ib:ie,jb:je,kb:ke))
    allocate(t_vav  (ib:ie,jb:je,kb:ke))
    allocate(t_tav  (ib:ie,jb:je,kb:ke))
    allocate(t_sgsav(ib:ie,jb:je,kb:ke))
    allocate(p_tav  (ib:ie,jb:je,kb:ke))
    allocate(p_bav  (ib:ie,jb:je,kb:ke))
    allocate(d_sgsav(ib:ie,jb:je,kb:ke))
    allocate(tkeadv(ib:ie,jb:je,kb:ke))

    allocate(t_p    (kb:ke))
    allocate(t_v    (kb:ke))
    allocate(t_t    (kb:ke))
    allocate(t_sgs  (kb:ke))
    allocate(p_t    (kb:ke))
    allocate(p_b    (kb:ke))
    allocate(d_sgs  (kb:ke))
    allocate(adv    (kb:ke))

    ! allocate wall shear-stress terms (immersed boundaries)
    allocate(shear(ib-1:ie+1,jb-1:je+1,kb-1:ke+1,0:12))    ! halo is set to 1
!     allocate(shear(ib:ie,jb-1:je+1,kb:ke,12)
    allocate(momfluxb(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
  !  write(*,*) "allocate momfluxb, indeces:",ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh 
   allocate(tfluxb(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(qfluxb(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
   allocate(cth(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))


    momfluxb=0.;tfluxb=0.;qfluxb=0.;cth=0.


    allocate(wall(ib:ie,jb:je,kb:ke,5))              ! 

    um=0.;u0=0.;up=0.
    vm=0.;v0=0.;vp=0.
    wm=0.;w0=0.;wp=0.
    pres0=0.;
    thlm=0.;thl0=0.;thlp=0.
    qtm=0.;qt0=0.;qtp=0.
    e12m=0.;e120=0.;e12p=0.
    svm=0.;sv0=0.;svp=0.;svpp=0.

    ql0=0.;qt0h=0.;
    thv0h=0.;thl0h=0.;
    mindist=1.0e10;
    presf=0.;presh=0.;exnf=1.;exnh=0.;thvf=0.;thvh=0.;rhof=0.    ! OG   
    !Exner function should be called in startup and just be initialised here
    qt0av=0.;ql0av=0.;thl0av=0.;u0av=0.;v0av=0.;sv0av=0.
    thlprof=0.;qtprof=0.;qlprof=0.;uprof=0.;vprof=0.;e12prof=0.;svprof=0.
    ug=0.;vg=0.;pgx=0.;pgy=0.;dpdxl=0.;dpdyl=0.;wfls=0.;whls=0.;thlpcar = 0.;uout=0.;vout=0.;wout=0.;udef=0.;vdef=0.;uouttot=0.;wouttot=0.;vouttot=0.
    dthldxls=0.;dthldyls=0.;dqtdxls=0.;dqtdyls=0.;dudxls=0.;dudyls=0.;dvdxls=0.;dvdyls=0.
    dthvdz=0.
    SW_up_TOA=0.;SW_dn_TOA=0.;LW_up_TOA=0.;LW_dn_TOA=0.

    uyt=0.;uytik=0.;vyt=0.;wyt=0.;wytik=0.;thlyt=0.;qtyt=0.;thlytk=0.;sca1yt=0.;sca2yt=0.;sca3yt=0.;thlsgsyt=0.;wsgsyt=0.;qtsgsyt=0.;sv1sgsyt=0.;sv2sgsyt=0.;sv3sgsyt=0.
    usgsyt=0.
    uxyt=0.;vxyt=0.;wxyt=0.;thlxyt=0.;qtxyt=0.;pxyt=0.;usgsxyt=0.;vsgsxyt=0.;thlsgsxyt=0.;
    uwtik=0.;wthltk=0.;wqttk=0.;thlthlt=0.;qtqtt=0.;sv1sv1t=0.;sv2sv2t=0.;sv3sv3t=0.;sv4sv4t=0.;uutc=0.;vvtc=0.;wwtc=0.;vwtjk=0.;uvtij=0.;utik=0.;wtik=0.;wtjk=0.;vtjk=0.;utij=0.;vtij=0.;
    wmt=0.;thltk=0.;qttk=0.;thlt=0.;slice=0.;slice2=0.;slice3=0.;slice4=0.;slice5=0.;utc=0.;vtc=0.;wtc=0.
    slice6=0.;slice7=0.;slice8=0.;umt=0.;vmt=0.;sv1t=0.;sv2t=0.;sv3t=0.;sv4t=0.;sv1tk=0.;sv2tk=0.;sv3tk=0.;sv4tk=0.
    wsv1tk=0.;wsv2tk=0.;wsv3tk=0.;wsv4tk=0.;sv1sgst=0.;sv2sgst=0.;sv3sgst=0.;sv4sgst=0.;qtt=0.;pt=0.
    PSSt = 0. !sv1max = 0.; sv2max = 0.; sv3max = 0.; sv4max = 0.

    scar=0.;scarl=0.

    IIc=1;IIu=1;IIv=1;IIct=1;IIw=1;IIuw=1;IIvw=1;IIuwt=1;IIut=1;IIvt=1;IIwt=1;IIcs=1;IIus=1;IIvs=1;IIws=1;IIuws=1;IIvws=1;IIuw=1;IIuvs=1

    uav=0.;vav=0.;wav=0.;thlav=0.;qtav=0.;svav=0.;viscratioav=0.;uuav=0.;vvav=0.
    wwav=0.;uvav=0.;uwav=0.;vwav=0.;sv2av=0.;thl2av=0.;ql2av=0.;qt2av=0.;presav=0.
    thluav=0.;thlvav=0.;thlwav=0.;thlthlav=0.;svuav=0.;svvav=0.;svwav=0.
    shear=0.
    upupav=0.;vpvpav=0.;wpwpav=0.;thlpthlpav=0.;qlpqlpav=0.;qtpqtpav=0.;svpsvpav=0.;upvpav=0.;upwpav=0.;vpwpav=0.
    thlpupav=0.;thlpvpav=0.;thlpwpav=0.;qlpupav=0.;qlpvpav=0.;qlpwpav=0.;qtpwpav=0.;qtpvpav=0.;qtpupav=0.;svpupav=0.;svpvpav=0.;svpwpav=0.
    umint=0.;vmint=0.;wmint=0.
! SGS
    uusgsav=0.;vvsgsav=0.;wwsgsav=0.;uwsgsav=0.;thlusgsav=0.;thlwsgsav=0.;qlusgsav=0.;qlwsgsav=0.;qtwsgsav=0.;qtusgsav=0.;
    svusgsav=0.;svwsgsav=0.;tkesgsav=0.;nusgsav=0.
! Resolved dissipation 
    strain2av=0.
! Subgrid dissipation 
    disssgsav=0.
! TKE budget
    t_vav=0.;tvmx=0.;tvmy=0.;tvmz=0.;tpm=0.;ttmx=0.;ttmy=0.;ttmz=0.;t_sgsav=0.;p_tav=0.
    tsgsmx1=0.;tsgsmy1=0.;tsgsmz1=0.;tsgsmx2=0.;tsgsmy2=0.;tsgsmz2=0.
    t_pav=0.;t_tav=0.;p_bav=0.;d_sgsav=0.;tkeadv=0.;t_p=0.;t_v=0.;t_t=0.;t_sgs=0.;p_t=0.;p_b=0.;d_sgs=0.;adv=0.
    ! domain fluid volume and area calculations
    uoutarea=0.;voutarea=0.;fluidvol=0.
  end subroutine initfields

  !> Deallocate the fields
  subroutine exitfields
    implicit none

    deallocate(um,vm,wm,thlm,e12m,qtm,u0,v0,w0,pres0,thl0,thl0h,qt0h,e120,qt0)
    deallocate(up,vp,wp,thlp,e12p,qtp)
    deallocate(svm,sv0,svp,svpp)
    deallocate(ql0,ql0h,thv0h,dthvdz,whls,presf,presh,exnf,exnh,thvf,thvh,rhof,qt0av,ql0av,thl0av,u0av,v0av)
    deallocate(ug,vg,pgx,pgy,dpdxl,dpdyl,dthldxls,dthldyls,dqtdxls,dqtdyls,dqtdtls,dudxls,dudyls,dvdxls,dvdyls,wfls)
    deallocate(thlprof,qtprof,uprof,vprof,e12prof,sv0av,svprof)
    deallocate(thlpcar)
    deallocate(momfluxb,tfluxb,qfluxb,cth)
    deallocate(SW_up_TOA,SW_dn_TOA,LW_up_TOA,LW_dn_TOA)

  end subroutine exitfields

end module modfields
