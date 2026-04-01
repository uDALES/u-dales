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
    real,    allocatable :: sveg(:)   ! absorbed shortwave on vegetation cells (W/m3)
    real,    allocatable :: tr_qtp(:)     ! qt tendency contribution on vegetation cells
    real,    allocatable :: tr_qtpR(:)    ! radiation-driven qt tendency contribution
    real,    allocatable :: tr_qtpA(:)    ! aerodynamic qt tendency contribution
    real,    allocatable :: tr_thlp(:)    ! thl tendency contribution on vegetation cells
    real,    allocatable :: tr_omega(:)   ! decoupling factor on vegetation cells
    real,    allocatable :: tr_svp(:,:)   ! scalar tendency contribution on vegetation cells
    logical :: has_sveg = .false.
  end type veg_type

  type(veg_type) :: veg

  type veg_face_type
    integer :: npts = 0               ! number of face points on this rank
    integer, allocatable :: ijk(:,:)  ! local grid indices (i,j,k) at face
    real,    allocatable :: dcoef(:)  ! face-centered (lad*cd)
    real,    allocatable :: tr_p(:)   ! momentum tendency contribution on face points
  end type veg_face_type

  type(veg_face_type) :: veg_u, veg_v, veg_w

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
    veg%has_sveg = .false.

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
      deallocate(veg%id, veg%gidx, veg%ijk, veg%lad, veg%cd, veg%ud, veg%dec, veg%lsize, veg%rs, veg%laiv, veg%sveg, &
                 veg%tr_qtp, veg%tr_qtpR, veg%tr_qtpA, veg%tr_thlp, veg%tr_omega, veg%tr_svp)
    end if

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
    allocate(veg%sveg(veg%npts))
    allocate(veg%tr_qtp(veg%npts))
    allocate(veg%tr_qtpR(veg%npts))
    allocate(veg%tr_qtpA(veg%npts))
    allocate(veg%tr_thlp(veg%npts))
    allocate(veg%tr_omega(veg%npts))
    allocate(veg%tr_svp(veg%npts,max(1,nsv)))

    veg%ijk = 0
    veg%laiv = 0.
    veg%sveg = 0.
    veg%tr_qtp = 0.
    veg%tr_qtpR = 0.
    veg%tr_qtpA = 0.
    veg%tr_thlp = 0.
    veg%tr_omega = 0.
    veg%tr_svp = 0.
    veg%has_sveg = sveg_exists

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
      if (veg%has_sveg) veg%sveg(m) = sveg_loc(m)
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
    if (allocated(veg_u%ijk)) then
      deallocate(veg_u%ijk, veg_u%dcoef, veg_u%tr_p)
      veg_u%npts = 0
    end if
    if (allocated(veg_v%ijk)) then
      deallocate(veg_v%ijk, veg_v%dcoef, veg_v%tr_p)
      veg_v%npts = 0
    end if
    if (allocated(veg_w%ijk)) then
      deallocate(veg_w%ijk, veg_w%dcoef, veg_w%tr_p)
      veg_w%npts = 0
    end if

    allocate(dcoef_f(ib:ie, jb:je, kb:ke))
    allocate(mask_f(ib:ie, jb:je, kb:ke))

    ! u-faces (i-1/i)
    dcoef_f = 0.5*(dcoef_3d(ib-1:ie-1, jb:je, kb:ke) + dcoef_3d(ib:ie, jb:je, kb:ke))
    mask_f = (dcoef_f > 0.0)

    veg_u%npts = count(mask_f)
    if (veg_u%npts > 0) then
      allocate(veg_u%ijk(veg_u%npts,3), veg_u%dcoef(veg_u%npts), veg_u%tr_p(veg_u%npts))
      idx = 0
      do k = kb, ke
        do j = jb, je
          do i = ib, ie
            if (mask_f(i,j,k)) then
              idx = idx + 1
              veg_u%ijk(idx,1:3) = (/ i, j, k /)
              veg_u%dcoef(idx) = dcoef_f(i,j,k)
            end if
          end do
        end do
      end do
      veg_u%tr_p = 0.
    end if

    ! v-faces (j-1/j)
    dcoef_f = 0.5*(dcoef_3d(ib:ie, jb-1:je-1, kb:ke) + dcoef_3d(ib:ie, jb:je, kb:ke))
    mask_f = (dcoef_f > 0.0)
    veg_v%npts = count(mask_f)
    if (veg_v%npts > 0) then
      allocate(veg_v%ijk(veg_v%npts,3), veg_v%dcoef(veg_v%npts), veg_v%tr_p(veg_v%npts))
      idx = 0
      do k = kb, ke
        do j = jb, je
          do i = ib, ie
            if (mask_f(i,j,k)) then
              idx = idx + 1
              veg_v%ijk(idx,1:3) = (/ i, j, k /)
              veg_v%dcoef(idx) = dcoef_f(i,j,k)
            end if
          end do
        end do
      end do
      veg_v%tr_p = 0.
    end if

    ! w-faces (k-1/k); no z-halos, so start at kb+1
    dcoef_f = 0.0
    dcoef_f(:,:,kb+1:ke) = 0.5*(dcoef_3d(ib:ie, jb:je, kb:ke-1) + dcoef_3d(ib:ie, jb:je, kb+1:ke))
    mask_f = (dcoef_f > 0.0)
    veg_w%npts = count(mask_f)
    if (veg_w%npts > 0) then
      allocate(veg_w%ijk(veg_w%npts,3), veg_w%dcoef(veg_w%npts), veg_w%tr_p(veg_w%npts))
      idx = 0
      do k = kb+1, ke
        do j = jb, je
          do i = ib, ie
            if (mask_f(i,j,k)) then
              idx = idx + 1
              veg_w%ijk(idx,1:3) = (/ i, j, k /)
              veg_w%dcoef(idx) = dcoef_f(i,j,k)
            end if
          end do
        end do
      end do
      veg_w%tr_p = 0.
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
    do mf = 1, veg_u%npts
      i = veg_u%ijk(mf,1)
      j = veg_u%ijk(mf,2)
      k = veg_u%ijk(mf,3)
      dcoefv = veg_u%dcoef(mf)
      veg_u%tr_p(mf) = - dcoefv * um(i,j,k) * &
                         sqrt( um(i,j,k)**2 &
                         + (0.25*(vm(i,j,k) + vm(i,j+1,k) + vm(i-1,j,k) + vm(i-1,j+1,k)))**2 &
                         + (0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i-1,j,k) + wm(i-1,j,k+1)))**2 )
      up(i,j,k) = up(i,j,k) + veg_u%tr_p(mf)
    end do

    do mf = 1, veg_v%npts
      i = veg_v%ijk(mf,1)
      j = veg_v%ijk(mf,2)
      k = veg_v%ijk(mf,3)
      dcoefv = veg_v%dcoef(mf)
      veg_v%tr_p(mf) = - dcoefv * vm(i,j,k) * &
                         sqrt( vm(i,j,k)**2 &
                         + (0.25*(um(i,j,k) + um(i+1,j,k) + um(i,j-1,k) + um(i+1,j-1,k)))**2 &
                         + (0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i,j-1,k) + wm(i,j-1,k+1)))**2 )
      vp(i,j,k) = vp(i,j,k) + veg_v%tr_p(mf)
    end do

    do mf = 1, veg_w%npts
      i = veg_w%ijk(mf,1)
      j = veg_w%ijk(mf,2)
      k = veg_w%ijk(mf,3)
      dcoefv = veg_w%dcoef(mf)
      veg_w%tr_p(mf) = - dcoefv * wm(i,j,k) * &
                         sqrt( wm(i,j,k)**2 &
                         + (0.25*(um(i,j,k) + um(i+1,j,k) + um(i,j,k-1) + um(i+1,j,k-1)))**2 &
                         + (0.25*(vm(i,j,k) + vm(i,j+1,k) + vm(i,j,k-1) + vm(i,j+1,k-1)))**2 )
      wp(i,j,k) = wp(i,j,k) + veg_w%tr_p(mf)
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
          veg%tr_svp(m,n) = veg%tr_svp(m,n) - svm(i,j,k,n) * ladv * udv
          svp(i,j,k,n) = svp(i,j,k,n) + veg%tr_svp(m,n)
        end do
      end do
    end if

  end subroutine apply_vegetation

  subroutine apply_vegetation_legacy
    use modglobal, only : Qstar, dzf
    implicit none
    integer :: i, j, k, m
    real :: ladv, decv, clai, rn_top, rn_bot, q_av_leaf

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

      call apply_vegetation_energy_point(m, i, j, k, ladv, veg%lsize(m), veg%rs(m), q_av_leaf)
    end do
  end subroutine apply_vegetation_legacy

  subroutine apply_vegetation_sveg
    implicit none
    integer :: i, j, k, m
    real :: ladv, q_av_leaf

    do m = 1, veg%npts
      i = veg%ijk(m,1)
      j = veg%ijk(m,2)
      k = veg%ijk(m,3)

      ladv = veg%lad(m)
      q_av_leaf = veg%sveg(m) / max(ladv, 1.0e-12)

      call apply_vegetation_energy_point(m, i, j, k, ladv, veg%lsize(m), veg%rs(m), q_av_leaf)
    end do
  end subroutine apply_vegetation_sveg

  subroutine apply_vegetation_energy_point(m, i, j, k, ladv, lsize_in, rs_in, q_av_leaf)
    use modglobal, only : pref0, rlv, cp, rv, rd, rhoa
    use modfields, only : thlm, qtm, qtp, thlp, um, vm, wm
    implicit none
    integer, intent(in) :: m, i, j, k
    real, intent(in) :: ladv, lsize_in, rs_in, q_av_leaf
    real :: e_sat, e_vap, d_vap, slope_sat, r_a, omega, qe, qh, gam
    real :: lsizev, rsv, wind2

    lsizev = max(lsize_in, 1.0e-6)
    rsv = max(rs_in, 1.0e-6)
    gam = (cp*pref0*rv)/(rlv*rd)

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

    veg%tr_omega(m) = omega
    veg%tr_qtpR(m) = ladv*(omega*(slope_sat/(slope_sat+2*gam))*q_av_leaf)/(rhoa*rlv)
    veg%tr_qtpA(m) = ladv*((1-omega)*(1/(gam*rsv))*rhoa*cp*d_vap)/(rhoa*rlv)
    veg%tr_qtp(m) = veg%tr_qtp(m) + ladv*qe/(rhoa*rlv)
    veg%tr_thlp(m) = veg%tr_thlp(m) + ladv*qh/(rhoa*cp)
    qtp(i,j,k) = qtp(i,j,k) + veg%tr_qtp(m)
    thlp(i,j,k) = thlp(i,j,k) + veg%tr_thlp(m)
  end subroutine apply_vegetation_energy_point

  subroutine reset_vegetation_sources()
    if (allocated(veg_u%tr_p)) veg_u%tr_p = 0.
    if (allocated(veg_v%tr_p)) veg_v%tr_p = 0.
    if (allocated(veg_w%tr_p)) veg_w%tr_p = 0.
    if (allocated(veg%tr_qtp)) veg%tr_qtp = 0.
    if (allocated(veg%tr_qtpR)) veg%tr_qtpR = 0.
    if (allocated(veg%tr_qtpA)) veg%tr_qtpA = 0.
    if (allocated(veg%tr_thlp)) veg%tr_thlp = 0.
    if (allocated(veg%tr_omega)) veg%tr_omega = 0.
    if (allocated(veg%tr_svp)) veg%tr_svp = 0.
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

    do mf = 1, veg_u%npts
      i = veg_u%ijk(mf,1); j = veg_u%ijk(mf,2); k = veg_u%ijk(mf,3)
      tr_up(i,j,k) = veg_u%tr_p(mf)
    end do
    do mf = 1, veg_v%npts
      i = veg_v%ijk(mf,1); j = veg_v%ijk(mf,2); k = veg_v%ijk(mf,3)
      tr_vp(i,j,k) = veg_v%tr_p(mf)
    end do
    do mf = 1, veg_w%npts
      i = veg_w%ijk(mf,1); j = veg_w%ijk(mf,2); k = veg_w%ijk(mf,3)
      tr_wp(i,j,k) = veg_w%tr_p(mf)
    end do

    do m = 1, veg%npts
      i = veg%ijk(m,1); j = veg%ijk(m,2); k = veg%ijk(m,3)
      tr_qtp(i,j,k) = veg%tr_qtp(m)
      tr_qtpR(i,j,k) = veg%tr_qtpR(m)
      tr_qtpA(i,j,k) = veg%tr_qtpA(m)
      tr_thlp(i,j,k) = veg%tr_thlp(m)
      tr_omega(i,j,k) = veg%tr_omega(m)
      if (nsv > 0) tr_sv(i,j,k,1:nsv) = veg%tr_svp(m,1:nsv)
    end do
  end subroutine scatter_vegetation_sources

end module vegetation
