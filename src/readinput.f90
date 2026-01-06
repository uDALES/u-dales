module readinput
  !> Generic file input routines for uDALES
  !> Provides standardized reading of commonly-used input file formats
  use mpi
  use modglobal, only : ifinput
  use modmpi,    only : myid, comm3d, mpierr
  use decomp_2d, only : zstart, zend
  implicit none

  private

  public :: read_sparse_ijk

contains

  !> Generic routine to read files containing sparse i,j,k grid point coordinates
  !> and return only the points that belong to the current MPI rank
  !> 
  !> Used by:
  !>  - heatpump.inp (2 header lines)
  !>  - solid_*.txt files (1 header line)
  !>  - fluid_boundary_*.txt files (1 header line)
  !>
  !> @param filename      Input filename to read
  !> @param npts          Number of data rows in file
  !> @param npts_loc      Output: number of points on this rank
  !> @param ids_loc       Output: indices of local points in global array (npts_loc)
  !> @param pts_loc       Output: local point coordinates remapped to local indices (npts_loc, 3)
  !> @param nskip         Number of header/comment lines to skip (default=1)
  
  subroutine read_sparse_ijk(filename, npts, npts_loc, ids_loc, pts_loc, nskip, pts_glob_out)
    implicit none

    character(len=*), intent(in)              :: filename      
    integer, intent(in)                       :: npts          
    integer, intent(out)                      :: npts_loc      
    integer, allocatable, intent(out)         :: ids_loc(:)    
    integer, allocatable, intent(out)         :: pts_loc(:,:)  
    integer, intent(in), optional             :: nskip   
    integer, allocatable, intent(out), optional :: pts_glob_out(:,:)   

    integer, allocatable :: pts_glob(:,:)
    logical, allocatable :: lpts(:)
    integer :: n, m, nh, ierr
    character(80) :: chmess

    ! Set default number of header lines to skip
    nh = 1
    if (present(nskip)) nh = nskip

    ! Allocate temporary arrays
    allocate(pts_glob(npts, 3))
    allocate(lpts(npts))

    ! Read file on rank 0 only
    if (myid == 0) then
      open(ifinput, file=filename, status='old', iostat=ierr)
      if (ierr /= 0) then
        write(*, '(A,A)') 'ERROR: Cannot open file: ', trim(filename)
        write(*, '(A,I0)') '  iostat error code: ', ierr
        stop 1
      end if

      ! Skip header lines
      do n = 1, nh
        read(ifinput, '(a80)', iostat=ierr) chmess
        if (ierr /= 0) then
          write(*, '(A,I0,A,A)') 'ERROR: Cannot read header line ', n, ' from ', trim(filename)
          stop 1
        end if
      end do

      ! Read integer triplets (i, j, k coordinates)
      do n = 1, npts
        read(ifinput, *, iostat=ierr) pts_glob(n, 1), pts_glob(n, 2), pts_glob(n, 3)
        if (ierr /= 0) then
          write(*, '(A,I0,A,A)') 'ERROR: Cannot read data line ', n, ' from ', trim(filename)
          write(*, '(A,I0)') '  iostat error code: ', ierr
          stop 1
        end if
      end do

      close(ifinput)
    end if

    ! Broadcast to all MPI ranks
    call MPI_BCAST(pts_glob, npts*3, MPI_INTEGER, 0, comm3d, mpierr)

    ! Determine whether points are on this rank and count them
    m = 0
    do n = 1, npts
      if ((pts_glob(n,1) >= zstart(1) .and. pts_glob(n,1) <= zend(1)) .and. &
          (pts_glob(n,2) >= zstart(2) .and. pts_glob(n,2) <= zend(2))) then
        lpts(n) = .true.
        m = m + 1
      else
        lpts(n) = .false.
      end if
    end do

    ! Pack local points into compact arrays
    npts_loc = m
    allocate(ids_loc(npts_loc))
    allocate(pts_loc(npts_loc, 3))
    
    m = 0
    do n = 1, npts
      if (lpts(n)) then
        m = m + 1
        ids_loc(m) = n
        ! Store global coordinates (no conversion)
        pts_loc(m,1) = pts_glob(n,1) - zstart(1) + 1
        pts_loc(m,2) = pts_glob(n,2) - zstart(2) + 1
        pts_loc(m,3) = pts_glob(n,3) - zstart(3) + 1
      end if
    end do

    ! Clean up or return temporary arrays
    if (present(pts_glob_out)) then
      ! Transfer ownership of global array to caller
      call move_alloc(pts_glob, pts_glob_out)
    else
      ! Clean up global array
      deallocate(pts_glob)
    end if
    deallocate(lpts)

  end subroutine read_sparse_ijk

end module readinput