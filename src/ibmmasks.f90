module ibmmasks
  use mpi
  use modglobal, only : libm, ib, ie, ihc, jb, je, jhc, kb, ke, khc, jtot, rslabs
  use ibmdata, only : solid_info_u, solid_info_v, solid_info_w, solid_info_c
  implicit none
  save

  integer, allocatable :: IIc(:,:,:)        !< Masking matrix for blocks at cell centres
  integer, allocatable :: IIu(:,:,:)        !< Masking matrix for blocks at x-direction half cells
  integer, allocatable :: IIv(:,:,:)        !< Masking matrix for blocks at y-direction half cells
  integer, allocatable :: IIw(:,:,:)        !< Masking matrix for blocks at z-direction half cells
  integer, allocatable :: IIuw(:,:,:)       !< Masking matrix for blocks at x-and z-direction half cells
  integer, allocatable :: IIvw(:,:,:)       !< Masking matrix for blocks at y- and z-direction half cells
  integer, allocatable :: IIuv(:,:,:)       !< Masking matrix for blocks at x- and y-direction half cells
  integer, allocatable :: IIct(:,:)         !< 2-D masking matrix for cell centres spanning 1:jtot
  integer, allocatable :: IIwt(:,:)         !< 2-D masking matrix for w-locations spanning 1:jtot
  integer, allocatable :: IIuwt(:,:)        !< 2-D masking matrix for uw-locations spanning 1:jtot
  integer, allocatable :: IIut(:,:)         !< 2-D masking matrix for u-locations spanning 1:jtot
  integer, allocatable :: IIvt(:,:)         !< 2-D masking matrix for v-locations spanning 1:jtot
  integer, allocatable :: IIcs(:)           !< 1-D masking matrix for cell centres spanning ib:ie and 1:jtot
  integer, allocatable :: IIus(:)           !< 1-D masking matrix for u-locations spanning ib:ie and 1:jtot
  integer, allocatable :: IIvs(:)           !< 1-D masking matrix for v-locations spanning ib:ie and 1:jtot
  integer, allocatable :: IIws(:)           !< 1-D masking matrix for w-locations spanning ib:ie and 1:jtot
  integer, allocatable :: IIuws(:)          !< 1-D masking matrix for uw-locations spanning ib:ie and 1:jtot
  integer, allocatable :: IIvws(:)          !< 1-D masking matrix for vw-locations spanning ib:ie and 1:jtot
  integer, allocatable :: IIuvs(:)          !< 1-D masking matrix for uv-locations spanning ib:ie and 1:jtot

contains

  subroutine initibm_support
    implicit none

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

    IIc = 1; IIu = 1; IIv = 1; IIw = 1
    IIuw = 1; IIvw = 1; IIuv = 1
    IIct = 1; IIwt = 1; IIuwt = 1; IIut = 1; IIvt = 1
    IIcs = 1; IIus = 1; IIvs = 1; IIws = 1
    IIuws = 1; IIvws = 1; IIuvs = 1
  end subroutine initibm_support

  subroutine createmasks
    implicit none

    integer :: IIcl(kb:ke + khc), IIul(kb:ke + khc), IIvl(kb:ke + khc), IIwl(kb:ke + khc)
    integer :: IIuwl(kb:ke + khc), IIvwl(kb:ke + khc), IIuvl(kb:ke + khc)
    integer :: IIcd(ib:ie, kb:ke), IIwd(ib:ie, kb:ke), IIuwd(ib:ie, kb:ke)
    integer :: IIud(ib:ie, kb:ke), IIvd(ib:ie, kb:ke)
    integer :: i, j, k, n, ierr

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

    IIc = 1; IIu = 1; IIv = 1; IIw = 1
    IIuw = 1; IIvw = 1; IIuv = 1
    IIct = 1; IIwt = 1; IIut = 1; IIvt = 1; IIuwt = 1
    IIcs = 1; IIus = 1; IIvs = 1; IIws = 1
    IIuws = 1; IIvws = 1; IIuvs = 1

    do n = 1, solid_info_u%nsolptsrank
      i = solid_info_u%solpts_loc(n,1)
      j = solid_info_u%solpts_loc(n,2)
      k = solid_info_u%solpts_loc(n,3)
      IIu(i,j,k) = 0
    end do

    do n = 1, solid_info_v%nsolptsrank
      i = solid_info_v%solpts_loc(n,1)
      j = solid_info_v%solpts_loc(n,2)
      k = solid_info_v%solpts_loc(n,3)
      IIv(i,j,k) = 0
    end do

    do n = 1, solid_info_w%nsolptsrank
      i = solid_info_w%solpts_loc(n,1)
      j = solid_info_w%solpts_loc(n,2)
      k = solid_info_w%solpts_loc(n,3)
      IIw(i,j,k) = 0
    end do

    do n = 1, solid_info_c%nsolptsrank
      i = solid_info_c%solpts_loc(n,1)
      j = solid_info_c%solpts_loc(n,2)
      k = solid_info_c%solpts_loc(n,3)
      IIc(i,j,k) = 0
    end do

    IIw(:, :, kb) = 0
    IIuw(:, :, kb) = 0
    IIvw(:, :, kb) = 0

    do i = ib, ie
      do j = jb, je
        IIuv(i,j,kb) = IIu(i,j,kb) * IIu(i,j-1,kb) * IIv(i,j,kb) * IIv(i-1,j,kb)
        do k = kb+1, ke
          IIuv(i,j,k) = IIu(i,j,k) * IIu(i,j-1,k) * IIv(i,j,k) * IIv(i-1,j,k)
          IIuw(i,j,k) = IIu(i,j,k) * IIu(i,j,k-1) * IIw(i,j,k) * IIw(i-1,j,k)
          IIvw(i,j,k) = IIv(i,j,k) * IIv(i,j,k-1) * IIw(i,j,k) * IIw(i,j-1,k)
        end do
      end do
    end do

    do k = kb, ke + khc
      IIcl(k) = sum(IIc(ib:ie, jb:je, k))
      IIul(k) = sum(IIu(ib:ie, jb:je, k))
      IIvl(k) = sum(IIv(ib:ie, jb:je, k))
      IIwl(k) = sum(IIw(ib:ie, jb:je, k))
      IIuwl(k) = sum(IIuw(ib:ie, jb:je, k))
      IIvwl(k) = sum(IIvw(ib:ie, jb:je, k))
      IIuvl(k) = sum(IIuv(ib:ie, jb:je, k))
    end do

    call MPI_ALLREDUCE(IIcl, IIcs, ke + khc - kb + 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(IIul, IIus, ke + khc - kb + 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(IIvl, IIvs, ke + khc - kb + 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(IIwl, IIws, ke + khc - kb + 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(IIuwl, IIuws, ke + khc - kb + 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(IIvwl, IIvws, ke + khc - kb + 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(IIuvl, IIuvs, ke + khc - kb + 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    IIcd(ib:ie, kb:ke) = sum(IIc(ib:ie, jb:je, kb:ke), DIM=2)
    IIwd(ib:ie, kb:ke) = sum(IIw(ib:ie, jb:je, kb:ke), DIM=2)
    IIuwd(ib:ie, kb:ke) = sum(IIuw(ib:ie, jb:je, kb:ke), DIM=2)
    IIud(ib:ie, kb:ke) = sum(IIu(ib:ie, jb:je, kb:ke), DIM=2)
    IIvd(ib:ie, kb:ke) = sum(IIv(ib:ie, jb:je, kb:ke), DIM=2)

    call MPI_ALLREDUCE(IIwd(ib:ie, kb:ke), IIwt(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(IIcd(ib:ie, kb:ke), IIct(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(IIuwd(ib:ie, kb:ke), IIuwt(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(IIud(ib:ie, kb:ke), IIut(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(IIvd(ib:ie, kb:ke), IIvt(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  end subroutine createmasks

end module ibmmasks
