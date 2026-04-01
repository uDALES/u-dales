!> Sparse vegetation support (veg_params-based masking and drag application)
module vegetation
  use mpi
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
    real, allocatable :: qtp(:)       ! qt tendency contribution on vegetation cells
    real, allocatable :: qtpR(:)      ! radiation-driven qt tendency contribution
    real, allocatable :: qtpA(:)      ! aerodynamic qt tendency contribution
    real, allocatable :: thlp(:)      ! thl tendency contribution on vegetation cells
    real, allocatable :: omega(:)     ! decoupling factor on vegetation cells
    real, allocatable :: svp(:,:)     ! scalar tendency contribution on vegetation cells
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
contains

  subroutine init_vegetation
    use modglobal,  only : ltrees,itree_mode,TREE_MODE_DRAG_ONLY,TREE_MODE_SVEG,TREE_MODE_LEGACY_SEB,ib,ie,jb,je,kb,ke,ih,jh,kh,cexpnr,dzf,nsv
    use modmpi,     only : myid,comm3d,mpierr,MY_REAL
    use readinput,  only : read_sparse_ijk, read_sparse_real
    use decomp_2d,  only : exchange_halo_x, exchange_halo_y, exchange_halo_z
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

    if (.not. ltrees) return
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
    if (allocated(vegp%qtp)) deallocate(vegp%qtp, vegp%qtpR, vegp%qtpA, vegp%thlp, vegp%omega, vegp%svp)

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
    allocate(vegp%qtp(veg%npts))
    allocate(vegp%qtpR(veg%npts))
    allocate(vegp%qtpA(veg%npts))
    allocate(vegp%thlp(veg%npts))
    allocate(vegp%omega(veg%npts))
    allocate(vegp%svp(veg%npts,max(1,nsv)))

    veg%ijk = 0
    veg%laiv = 0.
    sveg = 0.
    vegp%qtp = 0.
    vegp%qtpR = 0.
    vegp%qtpA = 0.
    vegp%thlp = 0.
    vegp%omega = 0.
    vegp%svp = 0.
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
    call exchange_halo_z(lad_3d)
    call exchange_halo_z(dcoef_3d)

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

    vegetation_ready = .true.

  end subroutine init_vegetation

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

  subroutine apply_vegetation
    use modglobal,  only : ib,ie,jb,je,kb,ke,lmoist,ltempeq,nsv,itree_mode,TREE_MODE_DRAG_ONLY,TREE_MODE_SVEG,TREE_MODE_LEGACY_SEB
    use modfields,  only : um,vm,wm,up,vp,wp,svp,svm
    implicit none
    integer :: i,j,k,m,n,npts
    integer :: mf
    real :: dcoefv, ladv, udv

    if (.not. vegetation_ready) return

    npts = veg%npts

    if (npts <= 0) return

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
        call apply_vegetation_legacy
      else if (itree_mode == TREE_MODE_SVEG) then
        call apply_vegetation_sveg
      end if
    end if

    ! ========================================================================
    ! Loop 3: Scalar deposition
    ! ========================================================================
    if (nsv > 0) then
      do m = 1, npts
        i = veg%ijk(m,1)
        j = veg%ijk(m,2)
        k = veg%ijk(m,3)

        ladv = veg%lad(m)
        udv  = veg%ud(m)

        do n = 1, nsv
          vegp%svp(m,n) = vegp%svp(m,n) - svm(i,j,k,n) * ladv * udv
          svp(i,j,k,n) = svp(i,j,k,n) + vegp%svp(m,n)
        end do
      end do
    end if

  end subroutine apply_vegetation

  subroutine apply_vegetation_legacy
    use modglobal, only : Qstar, dzf, pref0, rlv, cp, rv, rd, rhoa
    use modfields, only : thlm, qtm, qtp, thlp, um, vm, wm
    implicit none
    integer :: i, j, k, m
    real :: ladv, decv, clai, rn_top, rn_bot, q_av_leaf
    real :: e_sat, e_vap, d_vap, slope_sat, r_a, omega, qe, qh, gam
    real :: lsizev, rsv, wind2

    gam = (cp*pref0*rv)/(rlv*rd)
    do m = 1, veg%npts
      i = veg%ijk(m,1)
      j = veg%ijk(m,2)
      k = veg%ijk(m,3)

      ladv = veg%lad(m)
      decv = veg%dec(m)
      clai = veg%laiv(m)

      rn_top = Qstar * exp(-decv * (clai - ladv * dzf(k)))
      rn_bot = Qstar * exp(-decv * clai)
      q_av_leaf = (rn_top - rn_bot) / (dzf(k) * max(ladv, 1.0e-12))

      lsizev = max(veg%lsize(m), 1.0e-6)
      rsv = max(veg%rs(m), 1.0e-6)

      e_sat = 610.8*exp((17.27*(thlm(i,j,k)-273.15))/(thlm(i,j,k)-35.85))
      e_vap = (qtm(i,j,k) * pref0) / (0.378 * qtm(i,j,k) + 0.622)
      d_vap = max(e_sat - e_vap, 0.)
      slope_sat = (4098*e_sat)/((thlm(i,j,k)-35.85)**2)

      wind2 = max((0.5*(um(i,j,k)+um(i+1,j,k)))**2 &
                +(0.5*(vm(i,j,k)+vm(i,j+1,k)))**2 &
                +(0.5*(wm(i,j,k)+wm(i,j,k+1)))**2, 1.0e-12)
      r_a = 130*sqrt(lsizev / sqrt(wind2))

      omega = 1/(1 + 2*(gam/(slope_sat+2*gam)) * (rsv/r_a))
      qe = omega*(slope_sat/(slope_sat+2*gam))*q_av_leaf + (1-omega)*(1/(gam*rsv))*rhoa*cp*d_vap
      qh = q_av_leaf - qe

      vegp%omega(m) = omega
      vegp%qtpR(m) = ladv*(omega*(slope_sat/(slope_sat+2*gam))*q_av_leaf)/(rhoa*rlv)
      vegp%qtpA(m) = ladv*((1-omega)*(1/(gam*rsv))*rhoa*cp*d_vap)/(rhoa*rlv)
      vegp%qtp(m) = vegp%qtp(m) + ladv*qe/(rhoa*rlv)
      vegp%thlp(m) = vegp%thlp(m) + ladv*qh/(rhoa*cp)
      qtp(i,j,k) = qtp(i,j,k) + vegp%qtp(m)
      thlp(i,j,k) = thlp(i,j,k) + vegp%thlp(m)
    end do
  end subroutine apply_vegetation_legacy

  subroutine apply_vegetation_sveg
    use modglobal, only : pref0, rlv, cp, rv, rd, rhoa
    use modfields, only : thlm, qtm, qtp, thlp, um, vm, wm
    implicit none
    integer :: i, j, k, m
    real :: ladv, q_av_leaf
    real :: e_sat, e_vap, d_vap, slope_sat, r_a, omega, qe, qh, gam
    real :: lsizev, rsv, wind2

    gam = (cp*pref0*rv)/(rlv*rd)
    do m = 1, veg%npts
      i = veg%ijk(m,1)
      j = veg%ijk(m,2)
      k = veg%ijk(m,3)

      ladv = veg%lad(m)
      q_av_leaf = sveg(m) / max(ladv, 1.0e-12)

      lsizev = max(veg%lsize(m), 1.0e-6)
      rsv = max(veg%rs(m), 1.0e-6)

      e_sat = 610.8*exp((17.27*(thlm(i,j,k)-273.15))/(thlm(i,j,k)-35.85))
      e_vap = (qtm(i,j,k) * pref0) / (0.378 * qtm(i,j,k) + 0.622)
      d_vap = max(e_sat - e_vap, 0.)
      slope_sat = (4098*e_sat)/((thlm(i,j,k)-35.85)**2)

      wind2 = max((0.5*(um(i,j,k)+um(i+1,j,k)))**2 &
                +(0.5*(vm(i,j,k)+vm(i,j+1,k)))**2 &
                +(0.5*(wm(i,j,k)+wm(i,j,k+1)))**2, 1.0e-12)
      r_a = 130*sqrt(lsizev / sqrt(wind2))

      omega = 1/(1 + 2*(gam/(slope_sat+2*gam)) * (rsv/r_a))
      qe = omega*(slope_sat/(slope_sat+2*gam))*q_av_leaf + (1-omega)*(1/(gam*rsv))*rhoa*cp*d_vap
      qh = q_av_leaf - qe

      vegp%omega(m) = omega
      vegp%qtpR(m) = ladv*(omega*(slope_sat/(slope_sat+2*gam))*q_av_leaf)/(rhoa*rlv)
      vegp%qtpA(m) = ladv*((1-omega)*(1/(gam*rsv))*rhoa*cp*d_vap)/(rhoa*rlv)
      vegp%qtp(m) = vegp%qtp(m) + ladv*qe/(rhoa*rlv)
      vegp%thlp(m) = vegp%thlp(m) + ladv*qh/(rhoa*cp)
      qtp(i,j,k) = qtp(i,j,k) + vegp%qtp(m)
      thlp(i,j,k) = thlp(i,j,k) + vegp%thlp(m)
    end do
  end subroutine apply_vegetation_sveg

  subroutine reset_vegetation_sources()
    if (allocated(veg_up)) veg_up = 0.
    if (allocated(veg_vp)) veg_vp = 0.
    if (allocated(veg_wp)) veg_wp = 0.
    if (allocated(vegp%qtp)) vegp%qtp = 0.
    if (allocated(vegp%qtpR)) vegp%qtpR = 0.
    if (allocated(vegp%qtpA)) vegp%qtpA = 0.
    if (allocated(vegp%thlp)) vegp%thlp = 0.
    if (allocated(vegp%omega)) vegp%omega = 0.
    if (allocated(vegp%svp)) vegp%svp = 0.
  end subroutine reset_vegetation_sources

  subroutine scatter_vegetation_sources(tr_up, tr_vp, tr_wp, tr_qtp, tr_qtpR, tr_qtpA, tr_thlp, tr_sv, tr_omega)
    use modglobal, only : ib, ie, jb, je, kb, ke, nsv
    implicit none
    real, intent(out) :: tr_up(ib:ie,jb:je,kb:ke)
    real, intent(out) :: tr_vp(ib:ie,jb:je,kb:ke)
    real, intent(out) :: tr_wp(ib:ie,jb:je,kb:ke)
    real, intent(out) :: tr_qtp(ib:ie,jb:je,kb:ke)
    real, intent(out) :: tr_qtpR(ib:ie,jb:je,kb:ke)
    real, intent(out) :: tr_qtpA(ib:ie,jb:je,kb:ke)
    real, intent(out) :: tr_thlp(ib:ie,jb:je,kb:ke)
    real, intent(out) :: tr_sv(ib:ie,jb:je,kb:ke,max(1,nsv))
    real, intent(out) :: tr_omega(ib:ie,jb:je,kb:ke)
    integer :: m, mf, i, j, k

    tr_up = 0.
    tr_vp = 0.
    tr_wp = 0.
    tr_qtp = 0.
    tr_qtpR = 0.
    tr_qtpA = 0.
    tr_thlp = 0.
    tr_sv = 0.
    tr_omega = 0.

    do mf = 1, npts_u
      i = ijk_u(mf,1); j = ijk_u(mf,2); k = ijk_u(mf,3)
      tr_up(i,j,k) = veg_up(mf)
    end do
    do mf = 1, npts_v
      i = ijk_v(mf,1); j = ijk_v(mf,2); k = ijk_v(mf,3)
      tr_vp(i,j,k) = veg_vp(mf)
    end do
    do mf = 1, npts_w
      i = ijk_w(mf,1); j = ijk_w(mf,2); k = ijk_w(mf,3)
      tr_wp(i,j,k) = veg_wp(mf)
    end do

    do m = 1, veg%npts
      i = veg%ijk(m,1); j = veg%ijk(m,2); k = veg%ijk(m,3)
      tr_qtp(i,j,k) = vegp%qtp(m)
      tr_qtpR(i,j,k) = vegp%qtpR(m)
      tr_qtpA(i,j,k) = vegp%qtpA(m)
      tr_thlp(i,j,k) = vegp%thlp(m)
      tr_omega(i,j,k) = vegp%omega(m)
      if (nsv > 0) tr_sv(i,j,k,1:nsv) = vegp%svp(m,1:nsv)
    end do
  end subroutine scatter_vegetation_sources

end module vegetation
