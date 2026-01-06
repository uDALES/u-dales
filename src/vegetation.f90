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
  end type veg_type

  type(veg_type) :: veg

  logical :: trees_sparse_ready = .false.

contains

  subroutine createtrees
    use modglobal,  only : ltrees,ib,ie,jb,je,kb,ke,cexpnr
    use modmpi,     only : myid,comm3d,mpierr,MY_REAL
    use readinput,  only : read_sparse_ijk
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
      deallocate(veg%id, veg%gidx, veg%ijk, veg%lad, veg%cd, veg%ud, veg%dec, veg%lsize, veg%rs)
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
    veg%ijk = 0

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

    trees_sparse_ready = .true.

  end subroutine createtrees

  subroutine trees
    use modglobal,  only : ib,ie,jb,je,kb,ke
    use modfields,  only : um,vm,wm,tr_u,tr_v,tr_w
    implicit none
    integer :: i,j,k,m,npts
    real :: cdv, ladv

    if (.not. trees_sparse_ready) return

    npts = veg%npts

    do m = 1, npts
      i = veg%ijk(m,1)
      j = veg%ijk(m,2)
      k = veg%ijk(m,3)
      if (i < ib .or. i > ie .or. j < jb .or. j > je .or. k < kb .or. k > ke) cycle

      cdv   = veg%cd(m)
      ladv  = veg%lad(m)

      ! w face (vertical) at k
      tr_w(i,j,k) = - cdv * ladv * wm(i,j,k) * &
                    sqrt( wm(i,j,k)**2 &
                    + (0.25*(um(i,j,k) + um(i+1,j,k) + um(i,j,k-1) + um(i+1,j,k-1)))**2 &
                    + (0.25*(vm(i,j,k) + vm(i,j+1,k) + vm(i,j,k-1) + vm(i,j+1,k-1)))**2 )

      ! v face (y-direction) at j
      tr_v(i,j,k) = - cdv * ladv * vm(i,j,k) * &
                    sqrt( vm(i,j,k)**2 &
                    + (0.25*(um(i,j,k) + um(i+1,j,k) + um(i,j-1,k) + um(i+1,j-1,k)))**2 &
                    + (0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i,j-1,k) + wm(i,j-1,k+1)))**2 )

      ! u face (x-direction) at i
      tr_u(i,j,k) = - cdv * ladv * um(i,j,k) * &
                    sqrt( um(i,j,k)**2 &
                    + (0.25*(vm(i,j,k) + vm(i,j+1,k) + vm(i-1,j,k) + vm(i-1,j+1,k)))**2 &
                    + (0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i-1,j,k) + wm(i-1,j,k+1)))**2 )
    end do

  end subroutine trees

end module vegetation
