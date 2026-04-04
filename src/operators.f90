module operators
  use mpi
  use architecture, only : comm3d, mpierr, my_real, abort_run
  use definitions, only : LOC_C, LOC_U, LOC_V, LOC_W, LOC_UV, LOC_WU, LOC_VW
  use decomp_2d, only : zsize
  use modglobal, only : ib, ie, jb, je, kb, ke, dxfi, dyi, dzhi
  implicit none
  save

  ! Keep instantaneous field operators here. This currently covers reductions
  ! and intrinsic averages, and can grow toward sparkle-style divergence,
  ! gradient, and interpolation operators.

contains

  subroutine reduce_xy_sum(aver,var)
    implicit none

    real, intent(inout) :: aver(:)
    real, intent(in)    :: var(:,:,:)
    real                :: averl(size(aver))
    real                :: avers(size(aver))
    integer             :: k

    averl = 0.
    avers = 0.

    if (size(var,3) /= size(aver)) then
      write(*,*) 'reduce_xy_sum: output length ', size(aver), ' does not match z extent ', size(var,3)
      call abort_run(1)
    end if

    do k=1,size(aver)
      averl(k) = sum(var(:,:,k))
    end do

    call MPI_ALLREDUCE(averl, avers, size(aver), MY_REAL, MPI_SUM, comm3d, mpierr)
    aver = aver + avers
  end subroutine reduce_xy_sum

  subroutine reduce_yz_sum(aver,var)
    implicit none

    real, intent(inout) :: aver(:)
    real, intent(in)    :: var(:,:,:)
    real                :: averl(size(aver))
    real                :: avers(size(aver))
    integer             :: i

    averl = 0.
    avers = 0.

    if (size(var,1) /= size(aver)) then
      write(*,*) 'reduce_yz_sum: output length ', size(aver), ' does not match x extent ', size(var,1)
      call abort_run(1)
    end if

    do i=1,size(aver)
      averl(i) = sum(var(i,:,:))
    end do

    call MPI_ALLREDUCE(averl, avers, size(aver), MY_REAL, MPI_SUM, comm3d, mpierr)
    aver = aver + avers
  end subroutine reduce_yz_sum

  subroutine av_intr(aver,var,var_loc,kh,lnan)
    use ibmmasks, only : IIc, IIu, IIv, IIw, IIuv, IIuw, IIvw, IIcs, IIus, IIvs, IIws, IIuvs, IIuws, IIvws
    implicit none

    integer, intent(in) :: kh
    integer, intent(in) :: var_loc
    real, intent(out)   :: aver(:)
    real, intent(in)    :: var(:,:,:)
    logical, intent(in) :: lnan

    select case (var_loc)
    case (LOC_C)
      call compute_avg_xy_intrinsic(aver,var,IIc,IIcs,kh,lnan,'av_intr')
    case (LOC_U)
      call compute_avg_xy_intrinsic(aver,var,IIu,IIus,kh,lnan,'av_intr')
    case (LOC_V)
      call compute_avg_xy_intrinsic(aver,var,IIv,IIvs,kh,lnan,'av_intr')
    case (LOC_W)
      call compute_avg_xy_intrinsic(aver,var,IIw,IIws,kh,lnan,'av_intr')
    case (LOC_UV)
      call compute_avg_xy_intrinsic(aver,var,IIuv,IIuvs,kh,lnan,'av_intr')
    case (LOC_WU)
      call compute_avg_xy_intrinsic(aver,var,IIuw,IIuws,kh,lnan,'av_intr')
    case (LOC_VW)
      call compute_avg_xy_intrinsic(aver,var,IIvw,IIvws,kh,lnan,'av_intr')
    case default
      write(*,*) 'av_intr: unknown location selector ', var_loc
      call abort_run(1)
    end select
  end subroutine av_intr

  subroutine av_y_intr(aver,var,var_loc)
    use ibmmasks, only : IIc, IIu, IIv, IIw, IIct, IIut, IIvt, IIwt, IIuw, IIuwt
    implicit none

    real, intent(out)   :: aver(:,:)
    real, intent(in)    :: var(:,:,:)
    integer, intent(in) :: var_loc

    select case (var_loc)
    case (LOC_C)
      call compute_avg_y_intrinsic(aver,var,IIc,IIct,'av_y_intr')
    case (LOC_U)
      call compute_avg_y_intrinsic(aver,var,IIu,IIut,'av_y_intr')
    case (LOC_V)
      call compute_avg_y_intrinsic(aver,var,IIv,IIvt,'av_y_intr')
    case (LOC_W)
      call compute_avg_y_intrinsic(aver,var,IIw,IIwt,'av_y_intr')
    case (LOC_WU)
      call compute_avg_y_intrinsic(aver,var,IIuw,IIuwt,'av_y_intr')
    case default
      write(*,*) 'av_y_intr: unknown location selector ', var_loc
      call abort_run(1)
    end select
  end subroutine av_y_intr

  subroutine sum_y_intr(sumy,var,var_loc)
    use ibmmasks, only : IIc, IIu, IIv, IIw, IIuv, IIuw, IIvw
    implicit none

    real, intent(out)   :: sumy(:,:)
    real, intent(in)    :: var(:,:,:)
    integer, intent(in) :: var_loc

    select case (var_loc)
    case (LOC_C)
      call compute_sum_y_intrinsic(sumy,var,IIc,'sum_y_intr')
    case (LOC_U)
      call compute_sum_y_intrinsic(sumy,var,IIu,'sum_y_intr')
    case (LOC_V)
      call compute_sum_y_intrinsic(sumy,var,IIv,'sum_y_intr')
    case (LOC_W)
      call compute_sum_y_intrinsic(sumy,var,IIw,'sum_y_intr')
    case (LOC_UV)
      call compute_sum_y_intrinsic(sumy,var,IIuv,'sum_y_intr')
    case (LOC_WU)
      call compute_sum_y_intrinsic(sumy,var,IIuw,'sum_y_intr')
    case (LOC_VW)
      call compute_sum_y_intrinsic(sumy,var,IIvw,'sum_y_intr')
    case default
      write(*,*) 'sum_y_intr: unknown location selector ', var_loc
      call abort_run(1)
    end select
  end subroutine sum_y_intr

  subroutine sum_x_intr(sumx,var,var_loc)
    use ibmmasks, only : IIc, IIu, IIv, IIw, IIuv, IIuw, IIvw
    implicit none

    real, intent(out)   :: sumx(:,:)
    real, intent(in)    :: var(:,:,:)
    integer, intent(in) :: var_loc

    select case (var_loc)
    case (LOC_C)
      call compute_sum_x_intrinsic(sumx,var,IIc,'sum_x_intr')
    case (LOC_U)
      call compute_sum_x_intrinsic(sumx,var,IIu,'sum_x_intr')
    case (LOC_V)
      call compute_sum_x_intrinsic(sumx,var,IIv,'sum_x_intr')
    case (LOC_W)
      call compute_sum_x_intrinsic(sumx,var,IIw,'sum_x_intr')
    case (LOC_UV)
      call compute_sum_x_intrinsic(sumx,var,IIuv,'sum_x_intr')
    case (LOC_WU)
      call compute_sum_x_intrinsic(sumx,var,IIuw,'sum_x_intr')
    case (LOC_VW)
      call compute_sum_x_intrinsic(sumx,var,IIvw,'sum_x_intr')
    case default
      write(*,*) 'sum_x_intr: unknown location selector ', var_loc
      call abort_run(1)
    end select
  end subroutine sum_x_intr

  subroutine compute_avg_xy_intrinsic(aver,var,mask_3d,mask_1d,kh,lnan,label_prefix)
    implicit none

    real, intent(out)            :: aver(:)
    real, intent(in)             :: var(:,:,:)
    integer, intent(in)          :: mask_3d(:,:,:)
    integer, intent(in)          :: mask_1d(:)
    integer, intent(in)          :: kh
    logical, intent(in)          :: lnan
    character(len=*), intent(in) :: label_prefix

    integer :: i1, i2, j1, j2, k1, k2
    integer :: ii1, ii2, ij1, ij2, ik1, ik2
    integer :: is1, is2, k, nz_expected
    integer :: fluid_count(size(aver))
    real    :: averl(size(aver))
    real    :: avers(size(aver))

    nz_expected = zsize(3) + kh
    if (size(aver) /= nz_expected) then
      write(*,*) trim(label_prefix), ': unexpected output size ', size(aver), ' expected ', nz_expected
      call abort_run(1)
    end if

    call select_intrinsic_xy_windows(var, mask_3d, mask_1d, nz_expected, kh, trim(label_prefix), &
         i1, i2, j1, j2, k1, k2, ii1, ii2, ij1, ij2, ik1, ik2, is1, is2)

    fluid_count = mask_1d(is1:is2)
    averl = 0.
    avers = 0.

    do k=1,nz_expected
      averl(k) = sum(var(i1:i2,j1:j2,k1+k-1) * mask_3d(ii1:ii2,ij1:ij2,ik1+k-1))
    end do

    if ((.not. lnan) .and. (fluid_count(1)==0)) then
      averl(1) = sum(var(i1:i2,j1:j2,k1))
      fluid_count(1) = fluid_count(size(fluid_count)-kh)
    end if

    call MPI_ALLREDUCE(averl, avers, size(averl), MY_REAL, MPI_SUM, comm3d, mpierr)

    where (fluid_count==0)
      aver = -999.
    elsewhere
      aver = avers / fluid_count
    endwhere
  end subroutine compute_avg_xy_intrinsic

  subroutine compute_avg_y_intrinsic(aver,var,mask_3d,mask_2d,label_prefix)
    implicit none

    real, intent(out)            :: aver(:,:)
    real, intent(in)             :: var(:,:,:)
    integer, intent(in)          :: mask_3d(:,:,:)
    integer, intent(in)          :: mask_2d(:,:)
    character(len=*), intent(in) :: label_prefix

    real    :: avero(size(aver,1),size(aver,2))
    integer :: i1, i2, j1, j2, k1, k2
    integer :: ii1, ii2, ij1, ij2, ik1, ik2
    integer :: it1, it2, kt1, kt2

    if (size(aver,1) /= zsize(1) .or. size(aver,2) /= zsize(3)) then
      write(*,*) trim(label_prefix), ': unexpected output shape ', shape(aver)
      call abort_run(1)
    end if

    call select_intrinsic_y_windows(var, mask_3d, mask_2d, trim(label_prefix), &
         i1, i2, j1, j2, k1, k2, ii1, ii2, ij1, ij2, ik1, ik2, it1, it2, kt1, kt2)

    avero = sum(var(i1:i2,j1:j2,k1:k2) * mask_3d(ii1:ii2,ij1:ij2,ik1:ik2), DIM=2)
    aver = 0.
    call MPI_ALLREDUCE(avero, aver, size(avero), MY_REAL, MPI_SUM, comm3d, mpierr)

    where (mask_2d(it1:it2,kt1:kt2)==0)
      aver = -999.
    elsewhere
      aver = aver / mask_2d(it1:it2,kt1:kt2)
    endwhere
  end subroutine compute_avg_y_intrinsic

  subroutine compute_sum_y_intrinsic(sumy,var,mask_3d,label_prefix)
    implicit none

    real, intent(out)            :: sumy(:,:)
    real, intent(in)             :: var(:,:,:)
    integer, intent(in)          :: mask_3d(:,:,:)
    character(len=*), intent(in) :: label_prefix

    real    :: sumproc(size(sumy,1),size(sumy,2))
    integer :: i1, i2, j1, j2, k1, k2
    integer :: ii1, ii2, ij1, ij2, ik1, ik2

    if (size(sumy,1) /= zsize(1) .or. size(sumy,2) /= zsize(3)) then
      write(*,*) trim(label_prefix), ': unexpected output shape ', shape(sumy), ' expected ', zsize(1), zsize(3)
      call abort_run(1)
    end if

    call select_intrinsic_sum_windows(var, mask_3d, trim(label_prefix), &
         i1, i2, j1, j2, k1, k2, ii1, ii2, ij1, ij2, ik1, ik2)

    sumproc = sum(var(i1:i2,j1:j2,k1:k2) * mask_3d(ii1:ii2,ij1:ij2,ik1:ik2), DIM=2)
    sumy = 0.
    call MPI_ALLREDUCE(sumproc, sumy, size(sumproc), MY_REAL, MPI_SUM, comm3d, mpierr)
  end subroutine compute_sum_y_intrinsic

  subroutine compute_sum_x_intrinsic(sumx,var,mask_3d,label_prefix)
    implicit none

    real, intent(out)            :: sumx(:,:)
    real, intent(in)             :: var(:,:,:)
    integer, intent(in)          :: mask_3d(:,:,:)
    character(len=*), intent(in) :: label_prefix

    real    :: sumproc(size(sumx,1),size(sumx,2))
    integer :: i1, i2, j1, j2, k1, k2
    integer :: ii1, ii2, ij1, ij2, ik1, ik2

    if (size(sumx,1) /= zsize(2) .or. size(sumx,2) /= zsize(3)) then
      write(*,*) trim(label_prefix), ': unexpected output shape ', shape(sumx), ' expected ', zsize(2), zsize(3)
      call abort_run(1)
    end if

    call select_intrinsic_sum_windows(var, mask_3d, trim(label_prefix), &
         i1, i2, j1, j2, k1, k2, ii1, ii2, ij1, ij2, ik1, ik2)

    sumproc = sum(var(i1:i2,j1:j2,k1:k2) * mask_3d(ii1:ii2,ij1:ij2,ik1:ik2), DIM=1)
    sumx = 0.
    call MPI_ALLREDUCE(sumproc, sumx, size(sumproc), MY_REAL, MPI_SUM, comm3d, mpierr)
  end subroutine compute_sum_x_intrinsic

  subroutine select_intrinsic_xy_windows(var, mask_3d, mask_1d, interior_z, lower_halo_z, label_prefix, &
       i1, i2, j1, j2, k1, k2, ii1, ii2, ij1, ij2, ik1, ik2, is1, is2)
    implicit none

    real, intent(in)             :: var(:,:,:)
    integer, intent(in)          :: mask_3d(:,:,:)
    integer, intent(in)          :: mask_1d(:)
    integer, intent(in)          :: interior_z, lower_halo_z
    character(len=*), intent(in) :: label_prefix
    integer, intent(out)         :: i1, i2, j1, j2, k1, k2
    integer, intent(out)         :: ii1, ii2, ij1, ij2, ik1, ik2
    integer, intent(out)         :: is1, is2

    call interior_bounds_real3d(var, zsize(1), zsize(2), interior_z, lower_halo_z, 0, &
         i1, i2, j1, j2, k1, k2, trim(label_prefix)//' var')
    call interior_bounds_int3d(mask_3d, zsize(1), zsize(2), interior_z, lower_halo_z, 0, &
         ii1, ii2, ij1, ij2, ik1, ik2, trim(label_prefix)//' mask')
    call interior_bounds_int1d(mask_1d, interior_z, lower_halo_z, 0, is1, is2, trim(label_prefix)//' count')
  end subroutine select_intrinsic_xy_windows

  subroutine select_intrinsic_y_windows(var, mask_3d, mask_2d, label_prefix, &
       i1, i2, j1, j2, k1, k2, ii1, ii2, ij1, ij2, ik1, ik2, it1, it2, kt1, kt2)
    implicit none

    real, intent(in)             :: var(:,:,:)
    integer, intent(in)          :: mask_3d(:,:,:)
    integer, intent(in)          :: mask_2d(:,:)
    character(len=*), intent(in) :: label_prefix
    integer, intent(out)         :: i1, i2, j1, j2, k1, k2
    integer, intent(out)         :: ii1, ii2, ij1, ij2, ik1, ik2
    integer, intent(out)         :: it1, it2, kt1, kt2

    call interior_bounds_real3d(var, zsize(1), zsize(2), zsize(3), 0, 1, &
         i1, i2, j1, j2, k1, k2, trim(label_prefix)//' var')
    call interior_bounds_int3d(mask_3d, zsize(1), zsize(2), zsize(3), 0, 1, &
         ii1, ii2, ij1, ij2, ik1, ik2, trim(label_prefix)//' mask')
    call interior_bounds_int2d_xz(mask_2d, zsize(1), zsize(3), it1, it2, kt1, kt2, trim(label_prefix)//' plane')
  end subroutine select_intrinsic_y_windows

  subroutine select_intrinsic_sum_windows(var, mask_3d, label_prefix, &
       i1, i2, j1, j2, k1, k2, ii1, ii2, ij1, ij2, ik1, ik2)
    implicit none

    real, intent(in)             :: var(:,:,:)
    integer, intent(in)          :: mask_3d(:,:,:)
    character(len=*), intent(in) :: label_prefix
    integer, intent(out)         :: i1, i2, j1, j2, k1, k2
    integer, intent(out)         :: ii1, ii2, ij1, ij2, ik1, ik2

    call interior_bounds_real3d(var, zsize(1), zsize(2), zsize(3), 0, 1, &
         i1, i2, j1, j2, k1, k2, trim(label_prefix)//' var')
    call interior_bounds_int3d(mask_3d, zsize(1), zsize(2), zsize(3), 0, 1, &
         ii1, ii2, ij1, ij2, ik1, ik2, trim(label_prefix)//' mask')
  end subroutine select_intrinsic_sum_windows

  subroutine interior_bounds_real3d(arr, interior_x, interior_y, interior_z, lower_halo_z, upper_halo_z, &
       i1, i2, j1, j2, k1, k2, label_prefix)
    implicit none

    real, intent(in)             :: arr(:,:,:)
    integer, intent(in)          :: interior_x, interior_y, interior_z
    integer, intent(in)          :: lower_halo_z, upper_halo_z
    integer, intent(out)         :: i1, i2, j1, j2, k1, k2
    character(len=*), intent(in) :: label_prefix

    call interior_bounds_xy(size(arr,1), interior_x, i1, i2, trim(label_prefix)//' x')
    call interior_bounds_xy(size(arr,2), interior_y, j1, j2, trim(label_prefix)//' y')
    call interior_bounds_z(size(arr,3), interior_z, lower_halo_z, upper_halo_z, k1, k2, trim(label_prefix)//' z')
  end subroutine interior_bounds_real3d

  subroutine interior_bounds_int3d(arr, interior_x, interior_y, interior_z, lower_halo_z, upper_halo_z, &
       i1, i2, j1, j2, k1, k2, label_prefix)
    implicit none

    integer, intent(in)          :: arr(:,:,:)
    integer, intent(in)          :: interior_x, interior_y, interior_z
    integer, intent(in)          :: lower_halo_z, upper_halo_z
    integer, intent(out)         :: i1, i2, j1, j2, k1, k2
    character(len=*), intent(in) :: label_prefix

    call interior_bounds_xy(size(arr,1), interior_x, i1, i2, trim(label_prefix)//' x')
    call interior_bounds_xy(size(arr,2), interior_y, j1, j2, trim(label_prefix)//' y')
    call interior_bounds_z(size(arr,3), interior_z, lower_halo_z, upper_halo_z, k1, k2, trim(label_prefix)//' z')
  end subroutine interior_bounds_int3d

  subroutine interior_bounds_int2d_xz(arr, interior_x, interior_z, i1, i2, k1, k2, label_prefix)
    implicit none

    integer, intent(in)          :: arr(:,:)
    integer, intent(in)          :: interior_x, interior_z
    integer, intent(out)         :: i1, i2, k1, k2
    character(len=*), intent(in) :: label_prefix

    call interior_bounds_xy(size(arr,1), interior_x, i1, i2, trim(label_prefix)//' x')
    call interior_bounds_z(size(arr,2), interior_z, 0, 0, k1, k2, trim(label_prefix)//' z')
  end subroutine interior_bounds_int2d_xz

  subroutine interior_bounds_int1d(arr, interior_z, lower_halo_z, upper_halo_z, k1, k2, label_prefix)
    implicit none

    integer, intent(in)          :: arr(:)
    integer, intent(in)          :: interior_z, lower_halo_z, upper_halo_z
    integer, intent(out)         :: k1, k2
    character(len=*), intent(in) :: label_prefix

    call interior_bounds_z(size(arr), interior_z, lower_halo_z, upper_halo_z, k1, k2, trim(label_prefix)//' z')
  end subroutine interior_bounds_int1d

  subroutine interior_bounds_xy(actual_size, interior_size, start_idx, end_idx, label)
    implicit none

    integer, intent(in)          :: actual_size, interior_size
    integer, intent(out)         :: start_idx, end_idx
    character(len=*), intent(in) :: label
    integer                      :: halo_total

    halo_total = actual_size - interior_size
    if (halo_total < 0 .or. mod(halo_total, 2) /= 0) then
      write(*,*) trim(label), ': unsupported extent ', actual_size, ' interior ', interior_size
      call abort_run(1)
    end if

    start_idx = halo_total / 2 + 1
    end_idx = start_idx + interior_size - 1
  end subroutine interior_bounds_xy

  subroutine interior_bounds_z(actual_size, interior_size, lower_halo, upper_halo, start_idx, end_idx, label)
    implicit none

    integer, intent(in)          :: actual_size, interior_size, lower_halo, upper_halo
    integer, intent(out)         :: start_idx, end_idx
    character(len=*), intent(in) :: label

    if (actual_size == interior_size) then
      start_idx = 1
    elseif (actual_size == interior_size + lower_halo + upper_halo) then
      start_idx = lower_halo + 1
    elseif (actual_size == interior_size + lower_halo) then
      start_idx = lower_halo + 1
    elseif (actual_size == interior_size + upper_halo) then
      start_idx = 1
    else
      write(*,*) trim(label), ': unsupported extent ', actual_size, ' interior ', interior_size, ' halos ', lower_halo, upper_halo
      call abort_run(1)
    end if

    end_idx = start_idx + interior_size - 1
  end subroutine interior_bounds_z

end module operators
