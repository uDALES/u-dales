module tests
  use MPI
  use decomp_2d
  use modglobal  ! Import all global variables for comprehensive output
  use modmpi, only : comm3d, my_real, mpierr, myid, myidx, myidy, myid1dy, nprocs, nprocx, nprocy, cmyidx
!  use stats, only : stats_createnc_tavg1d, ncidt1d, nrect
  use modstat_nc, only : open_nc, define_nc, writestat_dims_nc, ncinfo, writeoffset
  !use json_module
  !use readparameters, only : readnamelists, readjsonconfig, writenamelists
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


    subroutine tests_io
      implicit none
      integer :: i, xdim, ydim, zdim, nrect, ncid, time_step
      real :: coord_value
      real, allocatable :: test_var(:,:,:)
      
      character(80) :: filename
      character(80) :: timeVar(1,4)
      character(80) :: varinfo(1,4)

      nrect = 0
      xdim = ie-ib+1
      ydim = je-jb+1
      zdim = ke-kb+1
      
      call MPI_BARRIER(comm3d, mpierr)

      do i = 0, nprocs - 1
         if (myid == i) then
          write(*,*) '========================================='
          write(*,*) ' Process myid = ', myid, ' of ', nprocs - 1
          write(*,*) '-----------------------------------------'
          write(*,*) ' Grid position: (', myidx, ',', myidy, ')  [nprocx=', nprocx, ', nprocy=', nprocy, ']'
      end if
      call MPI_BARRIER(comm3d, mpierr)
      end do
! ==================================== !

      !call stats_createnc_tavg1d
        if (myid1dy == 0) then
        
        filename = 'tests_t_out.xxx.xxx.nc'
        filename(13:15) = cmyidx  ! Only x-processor ID
        filename(17:19) = cexpnr  ! Experiment number

        call open_nc(filename, ncid, nrect, n1=xdim, n2=jtot, n3=zdim)

        ! Define time variable and dimensions if new file
        if (nrect == 0) then
          call ncinfo(timeVar(1,:), 'time', 'Time', 's', 'time')
          call define_nc(ncid, 1, timeVar)
          call writestat_dims_nc(ncid)
        end if

        ! Define test variable
        call ncinfo(varinfo(1,:), 'test_value', 'Test coordinate value', '1', 'tttt')
        call define_nc(ncid, 1, varinfo)

        if (myid == 0) then
          write(*,*) 'Created file:', trim(filename)
        end if
      end if

      if (myid == 0) then
        write(*,*) 'NetCDF file created!'
        !write(*,*) 'Array dimensions:', xdim, 'x', ydim, 'x', zdim
      end if     
! ============================================= !

      allocate(test_var(ib:ie, jb:je, kb:ke))
       coord_value = 100.0 * myidx + myidy
       test_var = coord_value


      call MPI_BARRIER(comm3d, mpierr)

      if (myid == 0) then
        write(*,*) 'Writing test data with coord_value = 10*myid1dx + myid1dy + 0.1*time_step...'
      end if

      do time_step  = 1, 5
         nrect = time_step
           call writeoffset(ncid, 'test_value', test_var, nrect, xdim, ydim, zdim)
           test_var = coord_value + 0.1 * time_step
           
           call MPI_BARRIER(comm3d, mpierr)
      end do

    end subroutine tests_io


end module tests
