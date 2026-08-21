!> Sparse vegetation module: reads per-cell vegetation properties from
!! input files (veg.inp, veg_params.inp, sveg.inp), builds face-averaged
!! drag coefficient lists for the staggered u/v/w grids, and applies
!! momentum drag, canopy energy balance, and scalar deposition forcing
!! each timestep.
module vegetation
  use mpi
#if defined(_GPU)
  use cudafor
#endif
  implicit none
  save

  ! Vegetation properties.
  type veg_type
    integer :: npts = 0               ! number of sparse points on this rank
    integer, allocatable :: id(:)     ! vegetation ID / class label per point
    integer, allocatable :: gidx(:)   ! global index (position in veg_params.inp / veg.inp)
    integer, allocatable :: ijk(:,:)  ! local grid indices (i,j,k)
    real,    allocatable :: lad(:)    ! leaf area density (m2/m3)
    real,    allocatable :: cd(:)     ! volumetric drag coefficient (-)
    real,    allocatable :: ud(:)     ! deposition velocity for scalars (m/s)
    real,    allocatable :: dec(:)    ! extinction coefficient
    real,    allocatable :: lsize(:)  ! characteristic leaf size (m)
    real,    allocatable :: rs(:)     ! stomatal resistance (s/m)
    real,    allocatable :: laiv(:)  ! cumulative LAI from domain top to top of this element
  end type veg_type

  type(veg_type) :: veg

  type vegp_type
    real, allocatable :: qt(:)        ! qt tendency contribution on vegetation cells
    real, allocatable :: qtR(:)       ! radiation-driven qt tendency contribution
    real, allocatable :: qtA(:)       ! aerodynamic qt tendency contribution
    real, allocatable :: thl(:)       ! thl tendency contribution on vegetation cells
    real, allocatable :: omega(:)     ! decoupling factor on vegetation cells
    real, allocatable :: sv(:,:)      ! scalar tendency contribution on vegetation cells
  end type vegp_type

  type(vegp_type) :: vegp

  logical :: has_sveg = .false.
  real, allocatable :: sveg(:)        ! absorbed shortwave on vegetation cells (W/m3)

  integer :: npts_u = 0, npts_v = 0, npts_w = 0
  integer, allocatable :: ijk_u(:,:), ijk_v(:,:), ijk_w(:,:)
  real, allocatable :: dcoef_u(:), dcoef_v(:), dcoef_w(:)
  real, allocatable :: veg_up(:), veg_vp(:), veg_wp(:)

  real, allocatable :: lad_3d(:,:,:)  ! cell-centered LAD with halos for face averaging
  real, allocatable :: dcoef_3d(:,:,:) ! cell-centered (lad*cd) with halos for face averaging

  logical :: vegetation_ready = .false.

  ! ==========================================================================
  ! Time-independent factors of the canopy energy balance.
  !
  ! Everything in the heat and moisture loop that depends only on the
  ! vegetation properties is folded in here once, at startup: the two
  ! Beer-Lambert exponentials of the legacy radiation extinction, the leaf-size
  ! and stomatal-resistance clamps, and the divisions by rhoa*rlv, rhoa*cp,
  ! gam*rs and 130*sqrt(lsize). What is left in the loop is the part that
  ! genuinely moves with the flow. These arrays are veg%npts long - a handful
  ! of words per vegetation cell - so they cost nothing against the 3D fields.
  !
  ! They also collapse the two canopy modes into one routine. The only thing
  ! that separated legacy from sveg was how q_av_leaf is formed, and that is
  ! now a table lookup times a mode-dependent scalar.
  ! ==========================================================================
  real, allocatable :: veg_qleaf(:)   !< absorbed radiation per unit leaf area.
                                      !!  sveg mode:   W/m2, complete.
                                      !!  legacy mode: per unit Qstar, so the
                                      !!  loop multiplies by Qstar. Qstar is a
                                      !!  namelist constant today, but leaving
                                      !!  it outside the cache keeps this valid
                                      !!  if dQdt is ever wired up.
  real, allocatable :: veg_rs_ra(:)   !< rs/(130*sqrt(lsize)); rs/r_a = veg_rs_ra*wind2**0.25
  real, allocatable :: veg_aero(:)    !< rhoa*cp/(gam*rs)
  real, allocatable :: veg_lad_rlv(:) !< lad/(rhoa*rlv)
  real, allocatable :: veg_lad_cp(:)  !< lad/(rhoa*cp)
  real, allocatable :: veg_lad_ud(:)  !< lad*ud, the scalar deposition rate
  real :: veg_twogam = 0.             !< 2*(cp*pref0*rv)/(rlv*rd)

#if defined(_GPU)
  ! Device mirrors of the point lists, the cache and the tendency diagnostics.
  ! They live here rather than in modcuda because this module reads modcuda's
  ! field mirrors, so the dependency has to run one way only; the facet arrays
  ! sit in modcuda for the mirror-image reason.
  integer, device, allocatable :: ijk_u_d(:,:), ijk_v_d(:,:), ijk_w_d(:,:), veg_ijk_d(:,:)
  real,    device, allocatable :: dcoef_u_d(:), dcoef_v_d(:), dcoef_w_d(:)
  real,    device, allocatable :: veg_up_d(:), veg_vp_d(:), veg_wp_d(:)
  real,    device, allocatable :: veg_qleaf_d(:), veg_rs_ra_d(:), veg_aero_d(:), &
                                  veg_lad_rlv_d(:), veg_lad_cp_d(:), veg_lad_ud_d(:)
  real,    device, allocatable :: vegp_qt_d(:), vegp_qtR_d(:), vegp_qtA_d(:), &
                                  vegp_thl_d(:), vegp_omega_d(:), vegp_sv_d(:,:)

  ! Pinned landing sites for the tree diagnostics, which are the only thing
  ! vegetation puts on the bus inside the time loop - and only under ltreedump,
  ! on the steps statsdump actually samples. Automatic arrays cannot be pinned,
  ! hence module scope.
  real, allocatable, pinned :: veg_stage_u(:), veg_stage_v(:), veg_stage_w(:), veg_stage_c(:)
#endif
contains

  !> Read sparse vegetation data, distribute across MPI ranks, build 3D
  !! LAD/drag fields with halo exchanges, and precompute staggered face
  !! lists used by vegetation_forcing.
  subroutine init_vegetation
    use modglobal,  only : ltrees,ltreedump,itree_mode,TREE_MODE_SVEG,TREE_MODE_LEGACY_SEB,ib,ie,jb,je,kb,ke,ih,jh,kh,cexpnr,nsv
    use modmpi,     only : myid,comm3d,mpierr,MY_REAL
    use readinput,  only : read_sparse_ijk, read_sparse_real
    use m_halo,     only : halo_exchange
    implicit none
    integer :: i,j,k,m
    integer :: npts
    integer, allocatable :: ids_loc(:)
    integer, allocatable :: pts_in(:,:)
    character(200) :: filename
    integer :: ierr, ifinput
    character(256) :: line
    integer, allocatable :: id_all(:)
    real, allocatable :: lad_all(:), cd_all(:), ud_all(:), dec_all(:), lsize_all(:), rs_all(:), sveg_loc(:)
    integer :: nread
    integer :: idx
    real, allocatable :: dcoef_f(:,:,:)
    logical, allocatable :: mask_f(:,:,:)
    logical :: sveg_exists
#if defined(_GPU)
    real, allocatable, device :: lad_3d_d(:,:,:), dcoef_3d_d(:,:,:)
#endif

    if (.not. ltrees) then
      ltreedump = .false.
      return
    end if
    vegetation_ready = .false.
    has_sveg = .false.

    ! count points in veg.inp.<expnr> to set npts
    ! TEMPORARY WHILE WE RUN TREES AND VEG TOGETHER
    ! npts is contained in ntrees namelist variable
    write(filename, '(A,A)') 'veg.inp.', trim(cexpnr)
    ifinput = 99
    npts = 0
    if (myid == 0) then
      open(ifinput, file=filename, status='old', iostat=ierr)
      if (ierr /= 0) then
        write(*, '(A,A)') 'ERROR: Cannot open file: ', trim(filename)
        stop 1
      end if
      read(ifinput, '(a256)', iostat=ierr) line  ! skip header
      do
        read(ifinput, *, iostat=ierr) i, j, k
        if (ierr /= 0) exit
        npts = npts + 1
      end do
      close(ifinput)
    end if
    call MPI_BCAST(npts, 1, MPI_INTEGER, 0, comm3d, mpierr)

    ! Read veg_params.<expnr> (aligned with veg.inp ordering)
    if (allocated(id_all)) deallocate(id_all, lad_all, cd_all, ud_all, dec_all, lsize_all, rs_all)
    allocate(id_all(npts))
    allocate(lad_all(npts))
    allocate(cd_all(npts))
    allocate(ud_all(npts))
    allocate(dec_all(npts))
    allocate(lsize_all(npts))
    allocate(rs_all(npts))
    if (myid == 0) then
      write(filename, '(A,A)') 'veg_params.inp.', trim(cexpnr)
      open(ifinput, file=filename, status='old', iostat=ierr)
      if (ierr /= 0) then
        write(*, '(A,A)') 'ERROR: Cannot open file: ', trim(filename)
        stop 1
      end if
      read(ifinput, '(a256)', iostat=ierr) line  ! skip header
      do nread = 1, npts
        read(ifinput, *, iostat=ierr) id_all(nread), lad_all(nread), cd_all(nread), ud_all(nread), dec_all(nread), lsize_all(nread), rs_all(nread)
        if (ierr /= 0) then
          write(*, '(A,I0)') 'ERROR reading veg_params line ', nread
          stop 1
        end if
      end do
      close(ifinput)

      write(filename, '(A,A)') 'sveg.inp.', trim(cexpnr)
      inquire(file=filename, exist=sveg_exists)
      if (itree_mode /= TREE_MODE_SVEG) then
        if (sveg_exists) then
          write(*,'(A,A)') 'NOTE: Found optional vegetation shortwave file: ', trim(filename)
        end if
      else
        if (.not. sveg_exists) then
          write(*,'(A,A)') 'ERROR: Missing vegetation shortwave file: ', trim(filename)
          stop 1
        end if
      end if

    end if

    call MPI_BCAST(id_all, npts, MPI_INTEGER, 0, comm3d, mpierr)
    call MPI_BCAST(lad_all, npts, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(cd_all, npts, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(ud_all, npts, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(dec_all, npts, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(lsize_all, npts, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(rs_all, npts, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(sveg_exists, 1, MPI_LOGICAL, 0, comm3d, mpierr)

    call read_sparse_ijk('veg.inp.'//trim(cexpnr), npts, veg%npts, ids_loc, pts_in, nskip=1)
    if (sveg_exists) then
      call read_sparse_real('sveg.inp.'//trim(cexpnr), npts, ids_loc, sveg_loc, nskip=1)
    end if

    if (allocated(veg%id)) then
      deallocate(veg%id, veg%gidx, veg%ijk, veg%lad, veg%cd, veg%ud, veg%dec, veg%lsize, veg%rs, veg%laiv)
    end if
    if (allocated(sveg)) deallocate(sveg)
    if (allocated(vegp%qt)) deallocate(vegp%qt, vegp%qtR, vegp%qtA, vegp%thl, vegp%omega, vegp%sv)

    allocate(veg%id(veg%npts))
    allocate(veg%gidx(veg%npts))
    allocate(veg%ijk(veg%npts,3))
    allocate(veg%lad(veg%npts))
    allocate(veg%cd(veg%npts))
    allocate(veg%ud(veg%npts))
    allocate(veg%dec(veg%npts))
    allocate(veg%lsize(veg%npts))
    allocate(veg%rs(veg%npts))
    allocate(veg%laiv(veg%npts))
    allocate(sveg(veg%npts))
    allocate(vegp%qt(veg%npts))
    allocate(vegp%qtR(veg%npts))
    allocate(vegp%qtA(veg%npts))
    allocate(vegp%thl(veg%npts))
    allocate(vegp%omega(veg%npts))
    allocate(vegp%sv(veg%npts,max(1,nsv)))

    veg%ijk = 0
    veg%laiv = 0.
    sveg = 0.
    vegp%qt = 0.
    vegp%qtR = 0.
    vegp%qtA = 0.
    vegp%thl = 0.
    vegp%omega = 0.
    vegp%sv = 0.
    has_sveg = sveg_exists

    do m = 1, veg%npts
      veg%id(m)    = id_all(ids_loc(m))
      veg%gidx(m)  = ids_loc(m)
      veg%lad(m)   = lad_all(ids_loc(m))
      veg%cd(m)    = cd_all(ids_loc(m))
      veg%ud(m)    = ud_all(ids_loc(m))
      veg%dec(m)   = dec_all(ids_loc(m))
      veg%lsize(m) = lsize_all(ids_loc(m))
      veg%rs(m)    = rs_all(ids_loc(m))
      veg%ijk(m,1:3) = pts_in(m,1:3)
      if (has_sveg) sveg(m) = sveg_loc(m)
    end do

    call MPI_BARRIER(comm3d, mpierr)

    deallocate(ids_loc)
    deallocate(pts_in)
    if (allocated(sveg_loc)) deallocate(sveg_loc)
    deallocate(id_all, lad_all, cd_all, ud_all, dec_all, lsize_all, rs_all)

    if (allocated(lad_3d)) deallocate(lad_3d)
    if (allocated(dcoef_3d)) deallocate(dcoef_3d)
    allocate(lad_3d(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
    allocate(dcoef_3d(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
    lad_3d = 0.
    dcoef_3d = 0.

    ! Build 3D LAD and drag parameter fields from sparse vegetation points.
    do m = 1, veg%npts
      i = veg%ijk(m,1)
      j = veg%ijk(m,2)
      k = veg%ijk(m,3)
      lad_3d(i, j, k) = veg%lad(m)
      dcoef_3d(i, j, k) = veg%lad(m) * veg%cd(m)
    end do

    ! Exchange halos so face averaging has neighbor values (2D decomposition)
#if defined(_GPU)
    allocate(lad_3d_d(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
    allocate(dcoef_3d_d(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
    lad_3d_d = lad_3d
    dcoef_3d_d = dcoef_3d
    call halo_exchange(lad_3d_d, 3)
    call halo_exchange(dcoef_3d_d, 3)
    lad_3d = lad_3d_d
    dcoef_3d = dcoef_3d_d
    deallocate(lad_3d_d, dcoef_3d_d)
#else
    call halo_exchange(lad_3d, 3)
    call halo_exchange(dcoef_3d, 3)
#endif

    if (itree_mode == TREE_MODE_LEGACY_SEB) then
      call init_vegetation_legacy
    end if

    ! ========================================================================
    ! Precompute face lists for staggered u/v/w points
    ! ========================================================================
    if (allocated(ijk_u)) deallocate(ijk_u, dcoef_u, veg_up)
    if (allocated(ijk_v)) deallocate(ijk_v, dcoef_v, veg_vp)
    if (allocated(ijk_w)) deallocate(ijk_w, dcoef_w, veg_wp)
    npts_u = 0
    npts_v = 0
    npts_w = 0

    allocate(dcoef_f(ib:ie, jb:je, kb:ke))
    allocate(mask_f(ib:ie, jb:je, kb:ke))

    ! u-faces (i-1/i)
    dcoef_f = 0.5*(dcoef_3d(ib-1:ie-1, jb:je, kb:ke) + dcoef_3d(ib:ie, jb:je, kb:ke))
    mask_f = (dcoef_f > 0.0)

    npts_u = count(mask_f)
    if (npts_u > 0) then
      allocate(ijk_u(npts_u,3), dcoef_u(npts_u), veg_up(npts_u))
      idx = 0
      do k = kb, ke
        do j = jb, je
          do i = ib, ie
            if (mask_f(i,j,k)) then
              idx = idx + 1
              ijk_u(idx,1:3) = (/ i, j, k /)
              dcoef_u(idx) = dcoef_f(i,j,k)
            end if
          end do
        end do
      end do
      veg_up = 0.
    end if

    ! v-faces (j-1/j)
    dcoef_f = 0.5*(dcoef_3d(ib:ie, jb-1:je-1, kb:ke) + dcoef_3d(ib:ie, jb:je, kb:ke))
    mask_f = (dcoef_f > 0.0)
    npts_v = count(mask_f)
    if (npts_v > 0) then
      allocate(ijk_v(npts_v,3), dcoef_v(npts_v), veg_vp(npts_v))
      idx = 0
      do k = kb, ke
        do j = jb, je
          do i = ib, ie
            if (mask_f(i,j,k)) then
              idx = idx + 1
              ijk_v(idx,1:3) = (/ i, j, k /)
              dcoef_v(idx) = dcoef_f(i,j,k)
            end if
          end do
        end do
      end do
      veg_vp = 0.
    end if

    ! w-faces (k-1/k); no z-halos, so start at kb+1
    dcoef_f = 0.0
    dcoef_f(:,:,kb+1:ke) = 0.5*(dcoef_3d(ib:ie, jb:je, kb:ke-1) + dcoef_3d(ib:ie, jb:je, kb+1:ke))
    mask_f = (dcoef_f > 0.0)
    npts_w = count(mask_f)
    if (npts_w > 0) then
      allocate(ijk_w(npts_w,3), dcoef_w(npts_w), veg_wp(npts_w))
      idx = 0
      do k = kb+1, ke
        do j = jb, je
          do i = ib, ie
            if (mask_f(i,j,k)) then
              idx = idx + 1
              ijk_w(idx,1:3) = (/ i, j, k /)
              dcoef_w(idx) = dcoef_f(i,j,k)
            end if
          end do
        end do
      end do
      veg_wp = 0.
    end if

    deallocate(dcoef_f, mask_f)

    ! Fold away everything in the forcing that does not move with the flow.
    ! This has to follow init_vegetation_legacy, which fills veg%laiv.
    call init_vegetation_cache

#if defined(_GPU)
    call init_vegetation_device
#endif

    vegetation_ready = .true.

  end subroutine init_vegetation

  !> Precompute the time-independent factors of the canopy energy balance and
  !! the scalar deposition rate. See the declarations above for what each
  !! array holds; vegetation_canopy_host documents the algebra that puts them
  !! in that form.
  subroutine init_vegetation_cache
    use modglobal, only : pref0, rlv, cp, rv, rd, rhoa, dzf, itree_mode, &
                          TREE_MODE_SVEG, TREE_MODE_LEGACY_SEB
    implicit none
    integer :: m, k
    real    :: gam, ladv, ladd, rsv, decv, clai, dzfk

    if (allocated(veg_qleaf)) then
      deallocate(veg_qleaf, veg_rs_ra, veg_aero, veg_lad_rlv, veg_lad_cp, veg_lad_ud)
    end if
    allocate(veg_qleaf(veg%npts))
    allocate(veg_rs_ra(veg%npts))
    allocate(veg_aero(veg%npts))
    allocate(veg_lad_rlv(veg%npts))
    allocate(veg_lad_cp(veg%npts))
    allocate(veg_lad_ud(veg%npts))

    gam = (cp*pref0*rv)/(rlv*rd)
    veg_twogam = 2.*gam

    ! Left at zero for TREE_MODE_DRAG_ONLY, where the canopy loop never runs.
    veg_qleaf = 0.

    do m = 1, veg%npts
      k    = veg%ijk(m,3)
      ladv = veg%lad(m)
      ladd = max(ladv, 1.0e-12)
      rsv  = max(veg%rs(m), 1.0e-6)

      veg_rs_ra(m)   = rsv / (130.*sqrt(max(veg%lsize(m), 1.0e-6)))
      veg_aero(m)    = rhoa*cp / (gam*rsv)
      veg_lad_rlv(m) = ladv / (rhoa*rlv)
      veg_lad_cp(m)  = ladv / (rhoa*cp)
      veg_lad_ud(m)  = ladv * veg%ud(m)

      if (itree_mode == TREE_MODE_SVEG) then
        veg_qleaf(m) = sveg(m) / ladd
      else if (itree_mode == TREE_MODE_LEGACY_SEB) then
        decv = veg%dec(m)
        clai = veg%laiv(m)
        dzfk = dzf(k)
        veg_qleaf(m) = (exp(-decv*(clai - ladv*dzfk)) - exp(-decv*clai)) / (dzfk*ladd)
      end if
    end do
  end subroutine init_vegetation_cache

  !> Compute cumulative LAI from domain top for legacy SEB radiation extinction.
  subroutine init_vegetation_legacy
    use modglobal, only : ib, ie, jb, je, kb, ke, dzf
    implicit none
    integer :: i, j, k, m
    real, allocatable :: lai_3d(:,:,:)

    allocate(lai_3d(ib:ie, jb:je, kb:ke+1))
    lai_3d = 0.

    do k = ke, kb, -1
      lai_3d(:, :, k) = lai_3d(:, :, k+1) + lad_3d(ib:ie, jb:je, k) * dzf(k)
    end do

    do m = 1, veg%npts
      i = veg%ijk(m,1)
      j = veg%ijk(m,2)
      k = veg%ijk(m,3)
      veg%laiv(m) = lai_3d(i, j, k)
    end do

    deallocate(lai_3d)
  end subroutine init_vegetation_legacy

  !> Apply vegetation forcing: momentum drag on staggered faces, canopy
  !! energy balance (heat/moisture), and scalar deposition.  Called once
  !! per timestep from the main time loop.
  !!
  !! On a GPU build every one of those loops runs on the device, against the
  !! field mirrors modcuda already holds. Nothing crosses the bus here: the
  !! tendencies are added straight into up_d/vp_d/wp_d/thlp_d/qtp_d/svp_d, and
  !! the per-cell diagnostics stay on the device until statsdump asks for them
  !! - see updateVegDiagHost.
  subroutine vegetation_forcing
    implicit none

    if (.not. vegetation_ready) return

    ! Kept exactly as it was: a rank with no vegetation points of its own does
    ! nothing at all, even where its halo carries a neighbour's vegetation and
    ! npts_u/v/w is therefore non-zero at the shared face.
    if (veg%npts <= 0) return

#if defined(_GPU)
    call vegetation_forcing_device
#else
    call vegetation_forcing_host
#endif

  end subroutine vegetation_forcing

#if !defined(_GPU) || defined(UDALES_DEBUG)
  !> Host branch of vegetation_forcing.
  !!
  !! Compiled into a Debug GPU build as well, so that
  !! tests_cuda.f90::test_vegetation_forcing can run it against the device
  !! branch on identical inputs rather than against hand-written expectations.
  subroutine vegetation_forcing_host
    use modglobal,  only : lmoist,ltempeq,nsv,itree_mode,TREE_MODE_SVEG,TREE_MODE_LEGACY_SEB,Qstar
    use modfields,  only : um,vm,wm,up,vp,wp,svp,svm
    implicit none
    integer :: i,j,k,m,n
    integer :: mf
    real :: dcoefv

    call reset_vegetation_sources()

    ! ========================================================================
    ! Loop 1: Momentum drag forces (precomputed staggered faces, no branching)
    ! ========================================================================
    do mf = 1, npts_u
      i = ijk_u(mf,1)
      j = ijk_u(mf,2)
      k = ijk_u(mf,3)
      dcoefv = dcoef_u(mf)
      veg_up(mf) = - dcoefv * um(i,j,k) * &
                         sqrt( um(i,j,k)**2 &
                         + (0.25*(vm(i,j,k) + vm(i,j+1,k) + vm(i-1,j,k) + vm(i-1,j+1,k)))**2 &
                         + (0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i-1,j,k) + wm(i-1,j,k+1)))**2 )
      up(i,j,k) = up(i,j,k) + veg_up(mf)
    end do

    do mf = 1, npts_v
      i = ijk_v(mf,1)
      j = ijk_v(mf,2)
      k = ijk_v(mf,3)
      dcoefv = dcoef_v(mf)
      veg_vp(mf) = - dcoefv * vm(i,j,k) * &
                         sqrt( vm(i,j,k)**2 &
                         + (0.25*(um(i,j,k) + um(i+1,j,k) + um(i,j-1,k) + um(i+1,j-1,k)))**2 &
                         + (0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i,j-1,k) + wm(i,j-1,k+1)))**2 )
      vp(i,j,k) = vp(i,j,k) + veg_vp(mf)
    end do

    do mf = 1, npts_w
      i = ijk_w(mf,1)
      j = ijk_w(mf,2)
      k = ijk_w(mf,3)
      dcoefv = dcoef_w(mf)
      veg_wp(mf) = - dcoefv * wm(i,j,k) * &
                         sqrt( wm(i,j,k)**2 &
                         + (0.25*(um(i,j,k) + um(i+1,j,k) + um(i,j,k-1) + um(i+1,j,k-1)))**2 &
                         + (0.25*(vm(i,j,k) + vm(i,j+1,k) + vm(i,j,k-1) + vm(i,j+1,k-1)))**2 )
      wp(i,j,k) = wp(i,j,k) + veg_wp(mf)
    end do

    ! ========================================================================
    ! Loop 2: Canopy energy balance (heat and moisture fluxes)
    ! ========================================================================
    if (lmoist .and. ltempeq) then
      if (itree_mode == TREE_MODE_LEGACY_SEB) then
        call vegetation_canopy_host(Qstar)
      else if (itree_mode == TREE_MODE_SVEG) then
        call vegetation_canopy_host(1.)
      end if
    end if

    ! ========================================================================
    ! Loop 3: Scalar deposition
    ! ========================================================================
    if (nsv > 0) then
      do m = 1, veg%npts
        i = veg%ijk(m,1)
        j = veg%ijk(m,2)
        k = veg%ijk(m,3)

        do n = 1, nsv
          vegp%sv(m,n) = vegp%sv(m,n) - svm(i,j,k,n) * veg_lad_ud(m)
          svp(i,j,k,n) = svp(i,j,k,n) + vegp%sv(m,n)
        end do
      end do
    end if

  end subroutine vegetation_forcing_host

  !> Canopy energy balance for the vegetation cells on this rank.
  !!
  !! qrad scales the cached radiation term: Qstar in the legacy Beer-Lambert
  !! mode, where veg_qleaf holds the extinction difference per unit Qstar, and
  !! 1 in sveg mode, where veg_qleaf already holds W/m2 of leaf. That single
  !! argument is all that ever distinguished the two modes.
  !!
  !! The algebra is the original, with the divisions folded away. Writing s for
  !! the saturation-curve slope and g for the psychrometric constant:
  !!
  !!   s = 4098*e_sat/(T - 35.85)**2, and the loop only ever needs the ratio
  !!       f = s/(s + 2g) = 4098*e_sat / (4098*e_sat + 2g*(T - 35.85)**2),
  !!   which removes the slope itself along with its division;
  !!
  !!   g/(s + 2g) = (1 - f)/2, so
  !!       omega = 1/(1 + 2*(g/(s+2g))*(rs/r_a)) = 1/(1 + (1 - f)*rs/r_a);
  !!
  !!   r_a = 130*sqrt(lsize/sqrt(wind2)) = 130*sqrt(lsize)*wind2**(-1/4), so
  !!       rs/r_a = veg_rs_ra * wind2**0.25.
  !!
  !! Three divisions are left, all of them genuinely per-step.
  subroutine vegetation_canopy_host(qrad)
    use modglobal, only : pref0
    use modfields, only : thlm, qtm, qtp, thlp, um, vm, wm
    implicit none
    real, intent(in) :: qrad
    integer :: i, j, k, m
    real :: q_av_leaf, e_sat, e_vap, d_vap, esat4098, f, omega, qe, qh, wind2

    do m = 1, veg%npts
      i = veg%ijk(m,1)
      j = veg%ijk(m,2)
      k = veg%ijk(m,3)

      q_av_leaf = qrad * veg_qleaf(m)

      e_sat = 610.8*exp((17.27*(thlm(i,j,k)-273.15))/(thlm(i,j,k)-35.85))
      e_vap = (qtm(i,j,k) * pref0) / (0.378 * qtm(i,j,k) + 0.622)
      d_vap = max(e_sat - e_vap, 0.)

      esat4098 = 4098.*e_sat
      f = esat4098 / (esat4098 + veg_twogam*(thlm(i,j,k)-35.85)**2)

      wind2 = max((0.5*(um(i,j,k)+um(i+1,j,k)))**2 &
                +(0.5*(vm(i,j,k)+vm(i,j+1,k)))**2 &
                +(0.5*(wm(i,j,k)+wm(i,j,k+1)))**2, 1.0e-12)

      omega = 1./(1. + (1.-f)*veg_rs_ra(m)*sqrt(sqrt(wind2)))

      qe = omega*f*q_av_leaf + (1.-omega)*veg_aero(m)*d_vap
      qh = q_av_leaf - qe

      vegp%omega(m) = omega
      vegp%qtR(m) = veg_lad_rlv(m)*omega*f*q_av_leaf
      vegp%qtA(m) = veg_lad_rlv(m)*(1.-omega)*veg_aero(m)*d_vap
      vegp%qt(m) = vegp%qt(m) + veg_lad_rlv(m)*qe
      vegp%thl(m) = vegp%thl(m) + veg_lad_cp(m)*qh
      qtp(i,j,k) = qtp(i,j,k) + vegp%qt(m)
      thlp(i,j,k) = thlp(i,j,k) + vegp%thl(m)
    end do
  end subroutine vegetation_canopy_host
#endif

#if defined(_GPU)
  !> Mirror the point lists and the startup cache onto the device.
  !!
  !! Called from init_vegetation, alongside modibm's init_ibm_device: none of
  !! this changes after startup, so it crosses the bus exactly once.
  subroutine init_vegetation_device
    use modglobal, only : nsv, ltreedump
    implicit none

    if (veg%npts > 0) then
      allocate(veg_ijk_d(veg%npts,3));   veg_ijk_d     = veg%ijk
      allocate(veg_qleaf_d(veg%npts));   veg_qleaf_d   = veg_qleaf
      allocate(veg_rs_ra_d(veg%npts));   veg_rs_ra_d   = veg_rs_ra
      allocate(veg_aero_d(veg%npts));    veg_aero_d    = veg_aero
      allocate(veg_lad_rlv_d(veg%npts)); veg_lad_rlv_d = veg_lad_rlv
      allocate(veg_lad_cp_d(veg%npts));  veg_lad_cp_d  = veg_lad_cp
      allocate(veg_lad_ud_d(veg%npts));  veg_lad_ud_d  = veg_lad_ud

      ! The diagnostics are written whole by whichever loops are enabled, so
      ! they need no per-step reset - but a disabled loop never writes its
      ! array at all, and statsdump still reads it. Zero here, as the host's
      ! reset_vegetation_sources does on the first step.
      allocate(vegp_qt_d(veg%npts));     vegp_qt_d     = 0.
      allocate(vegp_qtR_d(veg%npts));    vegp_qtR_d    = 0.
      allocate(vegp_qtA_d(veg%npts));    vegp_qtA_d    = 0.
      allocate(vegp_thl_d(veg%npts));    vegp_thl_d    = 0.
      allocate(vegp_omega_d(veg%npts));  vegp_omega_d  = 0.
      allocate(vegp_sv_d(veg%npts,max(1,nsv)));       vegp_sv_d = 0.
    end if

    if (npts_u > 0) then
      allocate(ijk_u_d(npts_u,3)); ijk_u_d   = ijk_u
      allocate(dcoef_u_d(npts_u)); dcoef_u_d = dcoef_u
      allocate(veg_up_d(npts_u));  veg_up_d  = 0.
    end if
    if (npts_v > 0) then
      allocate(ijk_v_d(npts_v,3)); ijk_v_d   = ijk_v
      allocate(dcoef_v_d(npts_v)); dcoef_v_d = dcoef_v
      allocate(veg_vp_d(npts_v));  veg_vp_d  = 0.
    end if
    if (npts_w > 0) then
      allocate(ijk_w_d(npts_w,3)); ijk_w_d   = ijk_w
      allocate(dcoef_w_d(npts_w)); dcoef_w_d = dcoef_w
      allocate(veg_wp_d(npts_w));  veg_wp_d  = 0.
    end if

    ! Only ltreedump ever moves any of this back, so a run without it never
    ! allocates the staging at all.
    if (ltreedump) then
      if (npts_u    > 0) allocate(veg_stage_u(npts_u))
      if (npts_v    > 0) allocate(veg_stage_v(npts_v))
      if (npts_w    > 0) allocate(veg_stage_w(npts_w))
      if (veg%npts  > 0) allocate(veg_stage_c(veg%npts))
    end if

  end subroutine init_vegetation_device

  !> Device branch of vegetation_forcing.
  !!
  !! One kernel per loop, over the same sparse lists the host walks. The face
  !! lists are built from a mask over the 3D grid, so each (i,j,k) appears at
  !! most once in ijk_u/v/w and the momentum scatter needs no atomics. The
  !! cell list comes from veg.inp and carries no such guarantee, so the writes
  !! into thlp/qtp/svp are atomic - which is also what reproduces the host, a
  !! repeated point being added twice there.
  subroutine vegetation_forcing_device
    use modglobal, only : lmoist,ltempeq,nsv,itree_mode,TREE_MODE_SVEG,TREE_MODE_LEGACY_SEB,Qstar
    use modcuda,   only : um_d, vm_d, wm_d, up_d, vp_d, wp_d
    implicit none
    integer :: i, j, k, mf
    real :: dcoefv

    !$acc parallel loop default(present) private(i, j, k, dcoefv)
    do mf = 1, npts_u
      i = ijk_u_d(mf,1)
      j = ijk_u_d(mf,2)
      k = ijk_u_d(mf,3)
      dcoefv = dcoef_u_d(mf)
      veg_up_d(mf) = - dcoefv * um_d(i,j,k) * &
                         sqrt( um_d(i,j,k)**2 &
                         + (0.25*(vm_d(i,j,k) + vm_d(i,j+1,k) + vm_d(i-1,j,k) + vm_d(i-1,j+1,k)))**2 &
                         + (0.25*(wm_d(i,j,k) + wm_d(i,j,k+1) + wm_d(i-1,j,k) + wm_d(i-1,j,k+1)))**2 )
      up_d(i,j,k) = up_d(i,j,k) + veg_up_d(mf)
    end do
    !$acc end parallel loop

    !$acc parallel loop default(present) private(i, j, k, dcoefv)
    do mf = 1, npts_v
      i = ijk_v_d(mf,1)
      j = ijk_v_d(mf,2)
      k = ijk_v_d(mf,3)
      dcoefv = dcoef_v_d(mf)
      veg_vp_d(mf) = - dcoefv * vm_d(i,j,k) * &
                         sqrt( vm_d(i,j,k)**2 &
                         + (0.25*(um_d(i,j,k) + um_d(i+1,j,k) + um_d(i,j-1,k) + um_d(i+1,j-1,k)))**2 &
                         + (0.25*(wm_d(i,j,k) + wm_d(i,j,k+1) + wm_d(i,j-1,k) + wm_d(i,j-1,k+1)))**2 )
      vp_d(i,j,k) = vp_d(i,j,k) + veg_vp_d(mf)
    end do
    !$acc end parallel loop

    !$acc parallel loop default(present) private(i, j, k, dcoefv)
    do mf = 1, npts_w
      i = ijk_w_d(mf,1)
      j = ijk_w_d(mf,2)
      k = ijk_w_d(mf,3)
      dcoefv = dcoef_w_d(mf)
      veg_wp_d(mf) = - dcoefv * wm_d(i,j,k) * &
                         sqrt( wm_d(i,j,k)**2 &
                         + (0.25*(um_d(i,j,k) + um_d(i+1,j,k) + um_d(i,j,k-1) + um_d(i+1,j,k-1)))**2 &
                         + (0.25*(vm_d(i,j,k) + vm_d(i,j+1,k) + vm_d(i,j,k-1) + vm_d(i,j+1,k-1)))**2 )
      wp_d(i,j,k) = wp_d(i,j,k) + veg_wp_d(mf)
    end do
    !$acc end parallel loop

    if (lmoist .and. ltempeq) then
      if (itree_mode == TREE_MODE_LEGACY_SEB) then
        call vegetation_canopy_device(Qstar)
      else if (itree_mode == TREE_MODE_SVEG) then
        call vegetation_canopy_device(1.)
      end if
    end if

    if (nsv > 0) call vegetation_deposition_device(nsv)

  end subroutine vegetation_forcing_device

  !> Device counterpart of vegetation_canopy_host. Same algebra, same cache.
  subroutine vegetation_canopy_device(qrad)
    use modglobal, only : pref0
    use modcuda,   only : thlm_d, qtm_d, thlp_d, qtp_d, um_d, vm_d, wm_d
    implicit none
    real, intent(in) :: qrad
    integer :: i, j, k, m, npts
    real :: q_av_leaf, e_sat, e_vap, d_vap, esat4098, f, omega, qe, qh, wind2
    real :: twogam

    ! Copied out of the derived type first: naming veg%npts in the loop bound
    ! makes the compiler ask for veg itself on the device, and its allocatable
    ! components have no business being there.
    npts   = veg%npts
    twogam = veg_twogam

    !$acc parallel loop default(present) &
    !$acc   private(i, j, k, q_av_leaf, e_sat, e_vap, d_vap, esat4098, f, omega, qe, qh, wind2)
    do m = 1, npts
      i = veg_ijk_d(m,1)
      j = veg_ijk_d(m,2)
      k = veg_ijk_d(m,3)

      q_av_leaf = qrad * veg_qleaf_d(m)

      e_sat = 610.8*exp((17.27*(thlm_d(i,j,k)-273.15))/(thlm_d(i,j,k)-35.85))
      e_vap = (qtm_d(i,j,k) * pref0) / (0.378 * qtm_d(i,j,k) + 0.622)
      d_vap = max(e_sat - e_vap, 0.)

      esat4098 = 4098.*e_sat
      f = esat4098 / (esat4098 + twogam*(thlm_d(i,j,k)-35.85)**2)

      wind2 = max((0.5*(um_d(i,j,k)+um_d(i+1,j,k)))**2 &
                +(0.5*(vm_d(i,j,k)+vm_d(i,j+1,k)))**2 &
                +(0.5*(wm_d(i,j,k)+wm_d(i,j,k+1)))**2, 1.0e-12)

      omega = 1./(1. + (1.-f)*veg_rs_ra_d(m)*sqrt(sqrt(wind2)))

      qe = omega*f*q_av_leaf + (1.-omega)*veg_aero_d(m)*d_vap
      qh = q_av_leaf - qe

      vegp_omega_d(m) = omega
      vegp_qtR_d(m) = veg_lad_rlv_d(m)*omega*f*q_av_leaf
      vegp_qtA_d(m) = veg_lad_rlv_d(m)*(1.-omega)*veg_aero_d(m)*d_vap
      vegp_qt_d(m)  = veg_lad_rlv_d(m)*qe
      vegp_thl_d(m) = veg_lad_cp_d(m)*qh

      !$acc atomic update
      qtp_d(i,j,k) = qtp_d(i,j,k) + vegp_qt_d(m)
      !$acc atomic update
      thlp_d(i,j,k) = thlp_d(i,j,k) + vegp_thl_d(m)
    end do
    !$acc end parallel loop

  end subroutine vegetation_canopy_device

  !> Device counterpart of the scalar deposition loop.
  subroutine vegetation_deposition_device(nsvl)
    use modcuda, only : svm_d, svp_d
    implicit none
    integer, intent(in) :: nsvl
    integer :: i, j, k, m, n, npts

    npts = veg%npts

    !$acc parallel loop default(present) private(i, j, k, n)
    do m = 1, npts
      i = veg_ijk_d(m,1)
      j = veg_ijk_d(m,2)
      k = veg_ijk_d(m,3)
      do n = 1, nsvl
        vegp_sv_d(m,n) = - svm_d(i,j,k,n) * veg_lad_ud_d(m)
        !$acc atomic update
        svp_d(i,j,k,n) = svp_d(i,j,k,n) + vegp_sv_d(m,n)
      end do
    end do
    !$acc end parallel loop

  end subroutine vegetation_deposition_device

  !> Bring the tree diagnostics down for statsdump.
  !!
  !! This is the only vegetation transfer inside the time loop, and it does not
  !! happen at all unless ltreedump is set. The caller gates it on the step
  !! statsdump will actually sample - modstatsdump::statsdump_will_sample -
  !! rather than on every step, so the cadence follows tsample and not dt.
  !!
  !! Nothing on the device reads these arrays, so there is no upload to match.
  subroutine updateVegDiagHost
    use modglobal, only : ltreedump, ltempeq, lmoist, nsv
    implicit none
    integer :: n

    if (.not. ltreedump) return
    if (.not. vegetation_ready) return
    if (veg%npts <= 0) return

    if (npts_u > 0) then
      veg_stage_u = veg_up_d
      veg_up      = veg_stage_u
    end if
    if (npts_v > 0) then
      veg_stage_v = veg_vp_d
      veg_vp      = veg_stage_v
    end if
    if (npts_w > 0) then
      veg_stage_w = veg_wp_d
      veg_wp      = veg_stage_w
    end if

    if (ltempeq) then
      veg_stage_c = vegp_thl_d
      vegp%thl    = veg_stage_c
    end if

    if (lmoist) then
      veg_stage_c = vegp_qt_d
      vegp%qt     = veg_stage_c
      veg_stage_c = vegp_qtR_d
      vegp%qtR    = veg_stage_c
      veg_stage_c = vegp_qtA_d
      vegp%qtA    = veg_stage_c
      veg_stage_c = vegp_omega_d
      vegp%omega  = veg_stage_c
    end if

    do n = 1, nsv
      veg_stage_c   = vegp_sv_d(:,n)
      vegp%sv(:,n)  = veg_stage_c
    end do

  end subroutine updateVegDiagHost
#endif

  !> Zero all vegetation source/tendency arrays before accumulation.
  subroutine reset_vegetation_sources()
    if (allocated(veg_up)) veg_up = 0.
    if (allocated(veg_vp)) veg_vp = 0.
    if (allocated(veg_wp)) veg_wp = 0.
    if (allocated(vegp%qt)) vegp%qt = 0.
    if (allocated(vegp%qtR)) vegp%qtR = 0.
    if (allocated(vegp%qtA)) vegp%qtA = 0.
    if (allocated(vegp%thl)) vegp%thl = 0.
    if (allocated(vegp%omega)) vegp%omega = 0.
    if (allocated(vegp%sv)) vegp%sv = 0.
  end subroutine reset_vegetation_sources

end module vegetation
