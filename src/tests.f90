module tests
  use MPI
  use decomp_2d
  use modglobal  ! Import all global variables for comprehensive output
  use modmpi, only : comm3d, my_real, mpierr, myid
  use json_module
  use modstartup, only : readnamelists, readjsonconfig, RUN, DOMAIN, PHYSICS, DYNAMICS, BC, INLET, DRIVER, WALLS, ENERGYBALANCE, SCALARS, CHEMISTRY, OUTPUT, TREES, PURIFS, HEATPUMP

  implicit none

  contains
    subroutine inittest
      integer :: nx=64, ny=64, nz=64
      !integer :: p_row=2, p_col=2
      integer :: p_row=0, p_col=0
      integer :: ierror

      call MPI_INIT(ierror)
      !call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierr)
      !call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      !call decomp_2d_init(nx,ny,nz,p_row,p_col)
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

    subroutine tests_json
      !-----------------------------------------------------------------|
      !                                                                 |
      !     JSON input/output tests                                     |
      !     1) Load namoptions file and output namelists (last proc)   |
      !     2) Load config.json file and output namelists (last proc)  |
      !                                                                 |
      !-----------------------------------------------------------------|
      
      implicit none
      
      integer :: ierr, last_proc, nprocs, myid
      
      ! Get MPI information
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
      last_proc = nprocs - 1
      
      if (myid == 0) then
        write(*,*) '========================================='
        write(*,*) 'JSON INPUT/OUTPUT TESTS'
        write(*,*) '========================================='
        write(*,*) 'Total processes:', nprocs
        write(*,*) 'Last process ID:', last_proc
        write(*,*) '========================================='
      end if
      
      ! Barrier to ensure all processes have completed the broadcasts
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      
      ! Last process outputs all namelist values to verify broadcast worked
      if (myid == last_proc) then
        write(*,*) '  Last process (', myid, ') received ALL JSON VALUES:'
        call output_all_namelists
      end if
      
      if (myid == 0) then
        write(*,*) '========================================='
        write(*,*) 'JSON TESTS COMPLETED'
        write(*,*) '========================================='
      end if
      
    end subroutine tests_json

    subroutine output_all_namelists
      !-----------------------------------------------------------------|
      !     Write all namelists to file to verify broadcast worked     |
      !-----------------------------------------------------------------|
      
      implicit none
      integer :: iunit = 99
      character(len=100) :: filename
      
      ! Create filename with process ID for debugging
      write(filename,'(a,i0,a)') 'namoptions_json_test.', iexpnr
      
      write(*,*) '  Writing all namelists to file: ', trim(filename)
      
      open(unit=iunit, file=filename, status='replace', action='write')
      
      ! Write each namelist section using Fortran's built-in namelist write
      write(iunit, nml=RUN)
      write(iunit, nml=DOMAIN) 
      write(iunit, nml=PHYSICS)
      write(iunit, nml=DYNAMICS)
      write(iunit, nml=BC)
      write(iunit, nml=INLET)
      write(iunit, nml=DRIVER)
      write(iunit, nml=WALLS)
      write(iunit, nml=ENERGYBALANCE)
      write(iunit, nml=SCALARS)
      write(iunit, nml=CHEMISTRY)
      write(iunit, nml=OUTPUT)
      write(iunit, nml=TREES)
      write(iunit, nml=PURIFS)
      write(iunit, nml=HEATPUMP)
      
      close(iunit)
      
      write(*,*) '  All namelists written to file successfully!'
      
    end subroutine output_all_namelists

end module tests
