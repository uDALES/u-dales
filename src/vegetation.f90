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

  logical :: trees_sparse_ready = .false.
contains

  subroutine init_vegetation
    use modglobal,  only : ltrees,ib,ie,jb,je,kb,ke,cexpnr,dzf
    use modmpi,     only : myid,comm3d,mpierr,MY_REAL
    use readinput,  only : read_sparse_ijk
    use decomp_2d,  only : DECOMP_2D_COMM_CART_Z
    implicit none
    integer :: i,j,k,m,kk
    integer :: npts
    integer, allocatable :: ids_loc(:)
    integer, allocatable :: pts_in(:,:)
    character(200) :: filename
    integer :: ierr, ifinput
    character(256) :: line
    integer, allocatable :: id_all(:)
    real, allocatable :: lad_all(:), cd_all(:), ud_all(:), dec_all(:), lsize_all(:), rs_all(:)
    integer :: nread
    real, allocatable :: lad_3d(:,:,:)      ! 3D LAD field (ib:ie, jb:je, kb:ke)
    real, allocatable :: lai_3d(:,:,:)  ! Cumulative LAI (ib:ie, jb:je, kb:ke+1)
    real, allocatable :: lai_from_above(:,:) ! LAI from ranks above (ib:ie, jb:je)
    real, allocatable :: lai_local_total(:,:) ! Total LAI through this rank (ib:ie, jb:je)
    real, allocatable :: lai_sum_to_here(:,:) ! Cumulative sum from rank 0 to this rank (ib:ie, jb:je)
    real, allocatable :: lai_global_total(:,:) ! Total LAI across all ranks (ib:ie, jb:je)

    if (.not. ltrees) return

    trees_sparse_ready = .false.

    ! count points in veg.inp.<expnr> to set npts
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
    ! Compute cumulative LAI from top of domain downward
    ! ========================================================================
    
    ! Allocate temporary arrays for LAI calculation
    allocate(lad_3d(ib:ie, jb:je, kb:ke))
    allocate(lai_3d(ib:ie, jb:je, kb:ke+1))
    allocate(lai_from_above(ib:ie, jb:je))
    allocate(lai_local_total(ib:ie, jb:je))
    allocate(lai_sum_to_here(ib:ie, jb:je))
    allocate(lai_global_total(ib:ie, jb:je))
    
    lad_3d = 0.
    lai_3d = 0.
    lai_from_above = 0.
    lai_local_total = 0.
    lai_sum_to_here = 0.
    lai_global_total = 0.
    
    ! Step 1: Build 3D LAD field from sparse vegetation points
    do m = 1, veg%npts
      i = veg%ijk(m,1)
      j = veg%ijk(m,2)
      k = veg%ijk(m,3)
      lad_3d(i, j, k) = veg%lad(m)
    end do
    
    ! Step 2: Compute local cumulative LAI by integrating downward from top
    ! Start at ke+1 (top boundary) with zero
    lai_3d(:, :, ke+1) = 0.
    do k = ke, kb, -1
      lai_3d(:, :, k) = lai_3d(:, :, k+1) + lad_3d(:, :, k) * dzf(k)
    end do
    
    ! Step 3: Get cumulative LAI from all ranks above using collective operations
    ! Store total LAI through this rank (at bottom, kb)
    lai_local_total = lai_3d(:, :, kb)
    
    ! Get cumulative sum from rank 0 through this rank (inclusive)
    call MPI_Scan(lai_local_total, lai_sum_to_here, (ie-ib+1)*(je-jb+1), MY_REAL, &
                  MPI_SUM, DECOMP_2D_COMM_CART_Z, mpierr)
    
    ! Get total LAI across all ranks in z-direction
    call MPI_Allreduce(lai_local_total, lai_global_total, (ie-ib+1)*(je-jb+1), MY_REAL, &
                       MPI_SUM, DECOMP_2D_COMM_CART_Z, mpierr)
    
    ! LAI from ranks above = total - sum through here
    ! (Assumes rank 0 is at bottom, higher ranks at top)
    lai_from_above = lai_global_total - lai_sum_to_here
    
    ! Step 4: Add contribution from above to all levels in this rank
    do k = kb, ke+1
      lai_3d(:, :, k) = lai_3d(:, :, k) + lai_from_above(:, :)
    end do
    
    ! Step 5: Sample cumulative LAI at each vegetation point location
    do m = 1, veg%npts
      i = veg%ijk(m,1)
      j = veg%ijk(m,2)
      k = veg%ijk(m,3)
      
      ! LAI above this element (at level k, before adding this element's contribution)
      veg%laiv(m) = lai_3d(i, j, k)
    end do
    
    ! Clean up temporary arrays
    deallocate(lad_3d, lai_3d, lai_from_above, lai_local_total, lai_sum_to_here, lai_global_total)

    trees_sparse_ready = .true.

  end subroutine init_vegetation

  subroutine apply_vegetation
    use modglobal,  only : ib,ie,jb,je,kb,ke,lmoist,ltempeq,nsv,pref0,rlv,cp,rv,rd,rhoa,dzf,&
                           Qstar,dQdt,dec,dy,ih,jh
    use modfields,  only : um,vm,wm,tr_u,tr_v,tr_w,tr_qt,tr_thl,tr_sv,thlm,qtm,qtp,thlp,svm,svp,&
                           tr_qtR,tr_qtA,tr_omega,Rn,qa
    use modibmdata, only : bctfz
    use modmpi,     only : myid
    implicit none
    integer :: i,j,k,m,n,npts
    real :: cdv, ladv, udv, lsizev, rsv, decv
    real :: e_sat, e_vap, D, s, r_a, omega, qe, qh, gam
    real :: clai, Rn_top, Rn_bot, qc, qa_local, Rq, shade

    if (.not. trees_sparse_ready) return

    ! Compute gamma parameter for psychrometric calculations
    gam = (cp*pref0*rv)/(rlv*rd)

    npts = veg%npts

    ! ========================================================================
    ! Loop 1: Momentum drag forces
    ! ========================================================================
    ! The below is not strictly correct since it does not fully account for staggering. 
    ! However, it is a simple method to avoid complications requiring halos in both the 
    ! 3d tr_* arrays and the sparse vegetation data.
    ! Moreover, it is more or less consistent with the original implementation in modtrees.
    do m = 1, npts
      i = veg%ijk(m,1)
      j = veg%ijk(m,2)
      k = veg%ijk(m,3)

      cdv  = veg%cd(m)
      ladv = veg%lad(m)

      ! u face (x-direction) at i
      tr_u(i,j,k) = - cdv * ladv * um(i,j,k) * &
                      sqrt( um(i,j,k)**2 &
                      + (0.25*(vm(i,j,k) + vm(i,j+1,k) + vm(i-1,j,k) + vm(i-1,j+1,k)))**2 &
                      + (0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i-1,j,k) + wm(i-1,j,k+1)))**2 )

      ! v face (y-direction) at j
      tr_v(i,j,k) = - cdv * ladv * vm(i,j,k) * &
                      sqrt( vm(i,j,k)**2 &
                      + (0.25*(um(i,j,k) + um(i+1,j,k) + um(i,j-1,k) + um(i+1,j-1,k)))**2 &
                      + (0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i,j-1,k) + wm(i,j-1,k+1)))**2 )

      ! w face (vertical) at k
      tr_w(i,j,k) = - cdv * ladv * wm(i,j,k) * &
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
        lsizev = veg%lsize(m)
        rsv    = veg%rs(m)
        decv   = veg%dec(m)

        ! Use pre-computed cumulative LAI from top of domain to this point
        clai = veg%laiv(m)
        
        ! Net radiation at top and bottom of vegetation element
        Rn_top = Qstar * exp(-decv * clai)
        Rn_bot = Qstar * exp(-decv * (clai + ladv * dzf(k)))
        
        ! Change in radiation (absorbed by this element)
        qc = Rn_top - Rn_bot
        
        ! Available energy (W/m²) - simplified, no storage term
        qa_local = qc
        
        ! Saturation vapour pressure
        e_sat = 610.8*exp((17.27*(thlm(i,j,k)-273.15))/(thlm(i,j,k)-35.85))

        ! Water vapour partial pressure
        e_vap = (qtm(i,j,k) * pref0) / (0.378 * qtm(i,j,k) + 0.622)

        ! Vapour pressure deficit
        D = max(e_sat - e_vap, 0.)

        ! Slope of saturation vapour pressure curve
        s = (4098*e_sat)/((thlm(i,j,k)-35.85)**2)

        ! Aerodynamic resistance
        r_a = 130*sqrt(lsizev/(sqrt((0.5*(um(i,j,k)+um(i+1,j,k)))**2 &
                                   +(0.5*(vm(i,j,k)+vm(i,j+1,k)))**2 &
                                   +(0.5*(wm(i,j,k)+wm(i,j,k+1)))**2)))

        ! Decoupling factor
        omega = 1/(1 + 2*(gam/(s+2*gam)) * (rsv/r_a))

        ! Latent heat flux (Penman-Monteith formulation)
        ! Available energy per unit volume: qa_local/(dzf(k)*ladv)
        qe = omega*(s/(s+2*gam))*(qa_local/(dzf(k)*ladv)) + (1-omega)*(1/(gam*rsv))*rhoa*cp*D

        ! Sensible heat flux (energy balance residual)
        qh = qa_local/(dzf(k)*ladv) - qe

        ! Store components for diagnostics
        tr_omega(i,j,k) = omega
        tr_qtR(i,j,k) = ladv*(omega*(s/(s+2*gam))*(qa_local/(dzf(k)*ladv)))/(rhoa*rlv)
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
