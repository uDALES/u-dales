module tests
  use MPI
  use decomp_2d
  use modglobal
  use modmpi, only : comm3d, my_real, mpierr, myid
  use readparameters, only : readnamelists, writenamelists
  implicit none

  contains
    subroutine inittest
      integer :: nx=64, ny=64, nz=64
      integer :: p_row=0, p_col=0
      integer :: ierror

      call MPI_INIT(ierror)
      call decomp_2d_init(nx,ny,nz,p_row,p_col)

      write(*,*) xstart
      write(*,*) ystart
      write(*,*) zstart
      write(*,*) xend
      write(*,*) yend
      write(*,*) zend
      write(*,*) xsize
      write(*,*) ysize
      write(*,*) zsize

    end subroutine inittest

    subroutine exittest
      integer :: ierror

      call decomp_2d_finalize
      call MPI_FINALIZE(ierror)

    end subroutine exittest

    subroutine tests_roundtrip
      implicit none

      integer :: ierr, last_proc, nprocs, myid

      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
      last_proc = nprocs - 1

      if (myid == 0) then
        write(*,*) '========================================='
        write(*,*) 'INPUT ROUND-TRIP TEST'
        write(*,*) '========================================='
        write(*,*) 'Total processes:', nprocs
        write(*,*) 'Last process ID:', last_proc
        write(*,*) '========================================='
      end if

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if (myid == last_proc) then
        write(*,*) '  Last process (', myid, ') received ALL INPUT VALUES:'
        call writenamelists
      end if

    end subroutine tests_roundtrip

end module tests
