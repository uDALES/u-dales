module tests
  use MPI
  use decomp_2d
  use modglobal
  use modmpi, only : comm3d, my_real, mpierr, myid
  use modstartup, only : init2decomp
  use readparameters, only : writenamelists
  implicit none

  contains
    subroutine init_tests
      implicit none

      integer :: ierr, myid

      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

      if (nprocx * nprocy > 32) then
        if (myid == 0) then
          write(*,*) 'ERROR: input round-trip test refuses to run with more than 32 MPI ranks.'
          write(*,*) 'Requested nprocx*nprocy =', nprocx * nprocy
        end if
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      end if

      call init2decomp

    end subroutine init_tests

    subroutine exit_tests
      implicit none

      integer :: ierr

      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
    end subroutine exit_tests

    subroutine tests_roundtrip
      implicit none

      integer :: ierr, nprocs, myid, last_proc

      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
      last_proc = nprocs - 1

      if (myid == 0) then
        write(*,*) '========================================='
        write(*,*) 'INPUT ROUND-TRIP TEST'
        write(*,*) '========================================='
        write(*,*) 'Total processes:', nprocs
        write(*,*) 'Requested nprocx*nprocy:', nprocx * nprocy
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
