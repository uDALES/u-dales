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

  type veg_face_type
    integer :: npts = 0               ! number of face points on this rank
    integer, allocatable :: ijk(:,:)  ! local grid indices (i,j,k) at face
    real,    allocatable :: dcoef(:)  ! face-centered (lad*cd)
  end type veg_face_type

  type(veg_face_type) :: veg_u, veg_v, veg_w

  real, allocatable :: lad_3d(:,:,:)  ! cell-centered LAD with halos for face averaging
  real, allocatable :: dcoef_3d(:,:,:) ! cell-centered (lad*cd) with halos for face averaging

  logical :: vegetation_ready = .false.
contains

  subroutine init_vegetation
    use modglobal,  only : ltrees,ib,ie,jb,je,kb,ke,ih,jh,kh,cexpnr,dzf
    use modmpi,     only : myid,comm3d,mpierr,MY_REAL
    use readinput,  only : read_sparse_ijk
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
    real, allocatable :: lad_all(:), cd_all(:), ud_all(:), dec_all(:), lsize_all(:), rs_all(:)
    integer :: nread
    real, allocatable :: lai_3d(:,:,:)      ! LAI (ib:ie, jb:je, kb:ke+1)
    integer :: idx
    real, allocatable :: dcoef_f(:,:,:)
    logical, allocatable :: mask_f(:,:,:)

    if (.not. ltrees) return
    vegetation_ready = .false.

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
    end if

    call MPI_BCAST(id_all, npts, MPI_INTEGER, 0, comm3d, mpierr)
    call MPI_BCAST(lad_all, npts, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(cd_all, npts, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(ud_all, npts, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(dec_all, npts, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(lsize_all, npts, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(rs_all, npts, MY_REAL, 0, comm3d, mpierr)

    call read_sparse_ijk('veg.inp.'//trim(cexpnr), npts, veg%npts, ids_loc, pts_in, nskip=1)

    if (allocated(veg%id)) then
      deallocate(veg%id, veg%gidx, veg%ijk, veg%lad, veg%cd, veg%ud, veg%dec, veg%lsize, veg%rs, veg%laiv)
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

    veg%ijk = 0
    veg%laiv = 0.

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
    end do

    call MPI_BARRIER(comm3d, mpierr)

    deallocate(ids_loc)
    deallocate(pts_in)
    deallocate(id_all, lad_all, cd_all, ud_all, dec_all, lsize_all, rs_all)

    ! ========================================================================
    ! Compute LAI from top of domain downward
    ! Note: Domain is decomposed in x and y (2D pencil decomposition)
    ! Each rank has the full vertical extent; we need x/y halos for face averaging.
    
    ! Allocate LAD with halos and LAI without halos
    if (allocated(lad_3d)) deallocate(lad_3d)
    if (allocated(dcoef_3d)) deallocate(dcoef_3d)
    allocate(lad_3d(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
    allocate(dcoef_3d(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh))
    lad_3d = 0.
    dcoef_3d = 0.
    allocate(lai_3d(ib:ie, jb:je, kb:ke+1))
    
    ! Step 1: Build 3D LAD and dcoef (drag parameter) field from 
    ! sparse vegetation points
    do m = 1, veg%npts
      i = veg%ijk(m,1)
      j = veg%ijk(m,2)
      k = veg%ijk(m,3)
      lad_3d(i, j, k) = veg%lad(m)
      dcoef_3d(i, j, k) = veg%lad(m) * veg%cd(m)
    end do

    ! Exchange halos so face averaging has neighbor values (2D decomposition)
    call exchange_halo_x(lad_3d)
    call exchange_halo_y(lad_3d)
    call exchange_halo_x(dcoef_3d)
    call exchange_halo_y(dcoef_3d)
    
    ! Step 2: Compute cumulative LAI by integrating downward from top
    ! Start at ke+1 (top of domain) with zero
    lai_3d = 0.
    do k = ke, kb, -1
      lai_3d(:, :, k) = lai_3d(:, :, k+1) + lad_3d(ib:ie, jb:je, k) * dzf(k)
    end do
    
    ! Step 3: Sample cumulative LAI at each vegetation point location
    do m = 1, veg%npts
      i = veg%ijk(m,1)
      j = veg%ijk(m,2)
      k = veg%ijk(m,3)
      
      veg%laiv(m) = lai_3d(i, j, k)
    end do
    
    ! Clean up temporary arrays
    deallocate(lai_3d)

    ! ========================================================================
    ! Precompute face lists for staggered u/v/w points 
    ! ========================================================================
    if (allocated(veg_u%ijk)) then
      deallocate(veg_u%ijk, veg_u%dcoef)
      veg_u%npts = 0
    end if
    if (allocated(veg_v%ijk)) then
      deallocate(veg_v%ijk, veg_v%dcoef)
      veg_v%npts = 0
    end if
    if (allocated(veg_w%ijk)) then
      deallocate(veg_w%ijk, veg_w%dcoef)
      veg_w%npts = 0
    end if

    allocate(dcoef_f(ib:ie, jb:je, kb:ke))
    allocate(mask_f(ib:ie, jb:je, kb:ke))

    ! u-faces (i-1/i)
    dcoef_f = 0.5*(dcoef_3d(ib-1:ie-1, jb:je, kb:ke) + dcoef_3d(ib:ie, jb:je, kb:ke))
    mask_f = (dcoef_f > 0.0)

    veg_u%npts = count(mask_f)
    if (veg_u%npts > 0) then
      allocate(veg_u%ijk(veg_u%npts,3), veg_u%dcoef(veg_u%npts))
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
    end if

    ! v-faces (j-1/j)
    dcoef_f = 0.5*(dcoef_3d(ib:ie, jb-1:je-1, kb:ke) + dcoef_3d(ib:ie, jb:je, kb:ke))
    mask_f = (dcoef_f > 0.0)
    veg_v%npts = count(mask_f)
    if (veg_v%npts > 0) then
      allocate(veg_v%ijk(veg_v%npts,3), veg_v%dcoef(veg_v%npts))
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
    end if

    ! w-faces (k-1/k); no z-halos, so start at kb+1
    dcoef_f = 0.0
    dcoef_f(:,:,kb+1:ke) = 0.5*(dcoef_3d(ib:ie, jb:je, kb:ke-1) + dcoef_3d(ib:ie, jb:je, kb+1:ke))
    mask_f = (dcoef_f > 0.0)
    veg_w%npts = count(mask_f)
    if (veg_w%npts > 0) then
      allocate(veg_w%ijk(veg_w%npts,3), veg_w%dcoef(veg_w%npts))
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
    end if

    deallocate(dcoef_f, mask_f)

    vegetation_ready = .true.

  end subroutine init_vegetation

  subroutine apply_vegetation
    use modglobal,  only : ib,ie,jb,je,kb,ke,lmoist,ltempeq,nsv,pref0,rlv,cp,rv,rd,rhoa,dzf,&
                 Qstar,dQdt,dec,dy,ih,jh,lad
    use modfields,  only : um,vm,wm,tr_u,tr_v,tr_w,tr_qt,tr_thl,tr_sv,thlm,qtm,qtp,thlp,svm,svp,&
                           tr_qtR,tr_qtA,tr_omega,Rn,qa
    use modibmdata, only : bctfz
    use modmpi,     only : myid
    implicit none
    integer :: i,j,k,m,n,npts
    integer :: mf
    real :: dcoefv, ladv, udv, lsizev, rsv, decv
    real :: e_sat, e_vap, D, s, r_a, omega, qe, qh, gam
    real :: clai, Rn_top, Rn_bot, q_av, Rq, shade

    if (.not. vegetation_ready) return

    ! Compute gamma parameter for psychrometric calculations
    gam = (cp*pref0*rv)/(rlv*rd)

    npts = veg%npts

    if (npts <= 0) return

    ! ========================================================================
    ! Loop 1: Momentum drag forces (precomputed staggered faces, no branching)
    ! ========================================================================
    do mf = 1, veg_u%npts
      i = veg_u%ijk(mf,1)
      j = veg_u%ijk(mf,2)
      k = veg_u%ijk(mf,3)
      dcoefv = veg_u%dcoef(mf)
      tr_u(i,j,k) = - dcoefv * um(i,j,k) * &
                      sqrt( um(i,j,k)**2 &
                      + (0.25*(vm(i,j,k) + vm(i,j+1,k) + vm(i-1,j,k) + vm(i-1,j+1,k)))**2 &
                      + (0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i-1,j,k) + wm(i-1,j,k+1)))**2 )
    end do

    do mf = 1, veg_v%npts
      i = veg_v%ijk(mf,1)
      j = veg_v%ijk(mf,2)
      k = veg_v%ijk(mf,3)
      dcoefv = veg_v%dcoef(mf)
      tr_v(i,j,k) = - dcoefv * vm(i,j,k) * &
                      sqrt( vm(i,j,k)**2 &
                      + (0.25*(um(i,j,k) + um(i+1,j,k) + um(i,j-1,k) + um(i+1,j-1,k)))**2 &
                      + (0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i,j-1,k) + wm(i,j-1,k+1)))**2 )
    end do

    do mf = 1, veg_w%npts
      i = veg_w%ijk(mf,1)
      j = veg_w%ijk(mf,2)
      k = veg_w%ijk(mf,3)
      dcoefv = veg_w%dcoef(mf)
      tr_w(i,j,k) = - dcoefv * wm(i,j,k) * &
                      sqrt( wm(i,j,k)**2 &
                      + (0.25*(um(i,j,k) + um(i+1,j,k) + um(i,j,k-1) + um(i+1,j,k-1)))**2 &
                      + (0.25*(vm(i,j,k) + vm(i,j+1,k) + vm(i,j,k-1) + vm(i,j+1,k-1)))**2 )
    end do

    ! ========================================================================
    ! Loop 2: Canopy energy balance (heat and moisture fluxes)
    ! ========================================================================
    if (lmoist .and. ltempeq) then
      do m = 1, npts
        i = veg%ijk(m,1)
        j = veg%ijk(m,2)
        k = veg%ijk(m,3)

        ladv   = veg%lad(m)
        lsizev = max(veg%lsize(m), 1.0e-6)
        rsv    = max(veg%rs(m), 1.0e-6)
        decv   = veg%dec(m)

        ! Use pre-computed cumulative LAI (from top down to bottom of this element)
        clai = veg%laiv(m)
        
        ! Net radiation at bottom and top of vegetation element
        ! clai includes this element, so use it for bottom radiation
        Rn_top = Qstar * exp(-decv * (clai - ladv * dzf(k)))
        Rn_bot = Qstar * exp(-decv * clai)
        
        ! Available energy (W/m²) - radiation absorbed by this element
        q_av = Rn_top - Rn_bot
        
        ! Saturation vapour pressure
        e_sat = 610.8*exp((17.27*(thlm(i,j,k)-273.15))/(thlm(i,j,k)-35.85))

        ! Water vapour partial pressure
        e_vap = (qtm(i,j,k) * pref0) / (0.378 * qtm(i,j,k) + 0.622)

        ! Vapour pressure deficit
        D = max(e_sat - e_vap, 0.)

        ! Slope of saturation vapour pressure curve
        s = (4098*e_sat)/((thlm(i,j,k)-35.85)**2)

        ! Aerodynamic resistance
        r_a = 130*sqrt(max(lsizev, 1.0e-6)/(sqrt(max((0.5*(um(i,j,k)+um(i+1,j,k)))**2 &
                 +(0.5*(vm(i,j,k)+vm(i,j+1,k)))**2 &
                 +(0.5*(wm(i,j,k)+wm(i,j,k+1)))**2, 1.0e-12))))

        ! Decoupling factor
        omega = 1/(1 + 2*(gam/(s+2*gam)) * (rsv/r_a))

        ! Latent heat flux (Penman-Monteith formulation)
        ! Available energy per unit volume: q_av/(dzf(k)*ladv)
        qe = omega*(s/(s+2*gam))*(q_av/(dzf(k)*max(ladv, 1.0e-12))) + (1-omega)*(1/(gam*rsv))*rhoa*cp*D

        ! Sensible heat flux (energy balance residual)
        qh = q_av/(dzf(k)*max(ladv, 1.0e-12)) - qe

        ! Store components for diagnostics
        tr_omega(i,j,k) = omega
        tr_qtR(i,j,k) = ladv*(omega*(s/(s+2*gam))*(q_av/(dzf(k)*ladv)))/(rhoa*rlv)
        tr_qtA(i,j,k) = ladv*((1-omega)*(1/(gam*rsv))*rhoa*cp*D)/(rhoa*rlv)

        ! Volumetric sources/sinks
        tr_qt(i,j,k) = tr_qt(i,j,k) + ladv*qe/(rhoa*rlv)
        tr_thl(i,j,k) = tr_thl(i,j,k) + ladv*qh/(rhoa*cp)

        ! Add to RHS
        qtp(i,j,k) = qtp(i,j,k) + ladv*qe/(rhoa*rlv)
        thlp(i,j,k) = thlp(i,j,k) + ladv*qh/(rhoa*cp)
      end do
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
          tr_sv(i,j,k,n) = tr_sv(i,j,k,n) - svm(i,j,k,n) * ladv * udv
          svp(i,j,k,n) = svp(i,j,k,n) - svm(i,j,k,n) * ladv * udv
        end do
      end do
    end if

  end subroutine apply_vegetation

end module vegetation
