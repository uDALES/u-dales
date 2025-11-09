module tests
  use MPI
  use decomp_2d
  use modglobal  ! Import all global variables for comprehensive output
  use modmpi, only : comm3d, my_real, mpierr, myid, myidx, myidy, myid1dy, nprocs, nprocx, nprocy, cmyidx, cmyidy
!  use stats, only : stats_createnc_tavg1d, ncidt1d, nrect
  use modstat_nc, only : open_nc, define_nc, writestat_dims_nc, ncinfo, writeoffset, writeoffset_1dx
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
      use netcdf
      implicit none
      integer :: i, k, xdim, ydim, zdim, nrect, ncid, time_step
      real :: coord_value
      real, allocatable :: test_var(:,:,:)
      
      character(80) :: filename
      character(80) :: timeVar(1,4)
      character(80) :: varinfo(1,4)
      
      ! islice variables (use global islice from modglobal)
      integer :: local_nislice, ncid_islice, nrec_islice, local_idx, ii
      real, allocatable :: islice_data(:,:,:)
      character(80) :: islice_filename
      integer :: nislice_count  ! Actual number of islices
      
      ! kslice variables (use global kslice from modglobal)
      integer :: ncid_kslice, nrec_kslice
      real, allocatable :: kslice_data(:,:,:)
      character(80) :: kslice_filename
      integer :: nkslice_count  ! Actual number of kslices
      
      ! jslice variables (use global jslice from modglobal)
      integer :: local_njslice, ncid_jslice, nrec_jslice, local_idy, jj
      real, allocatable :: jslice_data(:,:,:)
      character(80) :: jslice_filename
      integer :: njslice_count  ! Actual number of jslices

      nrect = 0
      xdim = ie-ib+1
      ydim = je-jb+1
      zdim = ke-kb+1
      
      ! Get number of islices, jslices and kslices from modglobal
      nislice_count = nislice
      njslice_count = njslice
      nkslice_count = nkslice
      
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
        
        filename = ' test_t_out.xxx.xxx.nc'
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
! Create islice files
! ============================================= !
      if (nislice_count > 0) then
        if (myid == 0) write(*,*) 'Creating islice files for ', nislice_count, ' slices at x=', islice(1:nislice_count)
        
        ! Count local islices
        local_nislice = 0
        do i = 1, nislice_count
          if ( (islice(i)-1)/nprocx == myidx) then
            local_nislice = local_nislice + 1
          end if
        end do
        
        if (local_nislice > 0 .and. myidy == 0) then
          islice_filename = 'test_islice.xxx.xxx.nc'
          islice_filename(13:15) = cmyidx
          islice_filename(17:19) = cexpnr
          
          nrec_islice = 0
          call open_nc(islice_filename, ncid_islice, nrec_islice, n1=local_nislice, n2=jtot, n3=zdim)
          
          if (nrec_islice == 0) then
            call ncinfo(timeVar(1,:), 'time', 'Time', 's', 'time')
            call define_nc(ncid_islice, 1, timeVar)
            call writestat_dims_nc(ncid_islice)
            call write_islice_coords(ncid_islice, local_nislice, nislice_count, islice)
          end if
          
          call ncinfo(varinfo(1,:), 'test_value', 'Test coordinate value', '1', 'tttt')
          call define_nc(ncid_islice, 1, varinfo)
          
          write(*,'(A,I3,A,I2,A)') 'Processor myidx=', myidx, ' created islice file with ', local_nislice, ' slices'
        end if
      end if
      
! ============================================= !
! Create jslice files
! ============================================= !
      if (njslice_count > 0) then
        if (myid == 0) write(*,*) 'Creating jslice files for ', njslice_count, ' slices at y=', jslice(1:njslice_count)
        
        ! Count local jslices (round-robin by myidy)
        local_njslice = 0
        do i = 1, njslice_count
          if ( (jslice(i)-1)/nprocy == myidy) then
            local_njslice = local_njslice + 1
          end if
        end do
        
        if (local_njslice > 0 .and. myidx == 0) then
          jslice_filename = 'test_jslice.xxx.xxx.nc'
          jslice_filename(13:15) = cmyidy
          jslice_filename(17:19) = cexpnr
          
          nrec_jslice = 0
          call open_nc(jslice_filename, ncid_jslice, nrec_jslice, n1=itot, n2=local_njslice, n3=zdim)
          
          if (nrec_jslice == 0) then
            call ncinfo(timeVar(1,:), 'time', 'Time', 's', 'time')
            call define_nc(ncid_jslice, 1, timeVar)
            call writestat_dims_nc(ncid_jslice)
            call write_jslice_coords(ncid_jslice, local_njslice, njslice_count, jslice)
          end if
          
          call ncinfo(varinfo(1,:), 'test_value', 'Test coordinate value', '1', 'tttt')
          call define_nc(ncid_jslice, 1, varinfo)
          
          write(*,'(A,I3,A,I2,A)') 'Processor myidy=', myidy, ' created jslice file with ', local_njslice, ' slices'
        end if
      end if
      
! ============================================= !
! Create kslice files
! ============================================= !
      if (nkslice_count > 0) then
        if (myid == 0) write(*,*) 'Creating kslice files for ', nkslice_count, ' slices at k=', kslice(1:nkslice_count)
        
        if (myidy == 0) then
          kslice_filename = 'test_kslice.xxx.xxx.nc'
          kslice_filename(13:15) = cmyidx
          kslice_filename(17:19) = cexpnr
          
          nrec_kslice = 0
          call open_nc(kslice_filename, ncid_kslice, nrec_kslice, n1=xdim, n2=jtot, n3=nkslice_count)
          
          if (nrec_kslice == 0) then
            call ncinfo(timeVar(1,:), 'time', 'Time', 's', 'time')
            call define_nc(ncid_kslice, 1, timeVar)
            call writestat_dims_nc(ncid_kslice)
            call write_kslice_coords(ncid_kslice, nkslice_count, kslice)
          end if
          
          call ncinfo(varinfo(1,:), 'test_value', 'Test coordinate value', '1', 'tttt')
          call define_nc(ncid_kslice, 1, varinfo)
          
          write(*,'(A,I3,A,I2,A)') 'Processor myidx=', myidx, ' created kslice file with ', nkslice_count, ' slices'
        end if
      end if
      
      call MPI_BARRIER(comm3d, mpierr)

! ============================================= !
! Allocate and fill test data
! ============================================= !
      allocate(test_var(ib:ie, jb:je, kb:ke))
      coord_value = 100.0 * myidx + myidy
      test_var = coord_value

      if (nislice_count > 0 .and. local_nislice > 0) allocate(islice_data(local_nislice, jb:je, kb:ke))
      if (njslice_count > 0 .and. local_njslice > 0) allocate(jslice_data(ib:ie, local_njslice, kb:ke))
      if (nkslice_count > 0) allocate(kslice_data(ib:ie, jb:je, nkslice_count))

      call MPI_BARRIER(comm3d, mpierr)

      if (myid == 0) then
        write(*,*) 'Writing test data for 5 time steps...'
      end if

! ============================================= !
! Write multiple time steps
! ============================================= !
      do time_step = 1, 5
         nrect = time_step
         
         ! Update test data
         test_var = coord_value + 0.1 * time_step
         
         ! Write full 3D field
         call writeoffset(ncid, 'test_value', test_var, nrect, xdim, ydim, zdim)
         
         ! Write islice
         if (nislice_count > 0 .and. local_nislice > 0) then
           islice_data = 0.0
           local_idx = 0
           do i = 1, nislice_count
             if ( (islice(i)-1)/nprocx == myidx) then
               local_idx = local_idx + 1
               ii = islice(i)
               if (ii >= xstart(1) .and. ii <= xend(1)) then
                 islice_data(local_idx, :, :) = test_var(ii - xstart(1) + ib, :, :)
               end if
             end if
           end do
           
           nrec_islice = time_step
           call writeoffset(ncid_islice, 'test_value', islice_data, nrec_islice, local_nislice, ydim, zdim)
         end if
         
         ! Write jslice
         if (njslice_count > 0 .and. local_njslice > 0) then
           jslice_data = 0.0
           local_idy = 0
           do i = 1, njslice_count
             if ( (jslice(i)-1)/nprocy == myidy) then
               local_idy = local_idy + 1
               jj = jslice(i)
               if (jj >= ystart(2) .and. jj <= yend(2)) then
                 jslice_data(:, local_idy, :) = test_var(:, jj - ystart(2) + jb, :)
               end if
             end if
           end do
           
           nrec_jslice = time_step
           call writeoffset_1dx(ncid_jslice, 'test_value', jslice_data, nrec_jslice, xdim, local_njslice, zdim)
         end if
         
         ! Write kslice
         if (nkslice_count > 0) then
           kslice_data = 0.0
           do k = 1, nkslice_count
             kslice_data(:, :, k) = test_var(:, :, kslice(k))
           end do
           
           nrec_kslice = time_step
           call writeoffset(ncid_kslice, 'test_value', kslice_data, nrec_kslice, xdim, ydim, nkslice_count)
         end if
         
         call MPI_BARRIER(comm3d, mpierr)
      end do
      
      ! Cleanup
      deallocate(test_var)
      if (allocated(islice_data)) deallocate(islice_data)
      if (allocated(jslice_data)) deallocate(jslice_data)
      if (allocated(kslice_data)) deallocate(kslice_data)
      
      if (myid == 0) then
        write(*,*) 'Test completed! Files written:'
        write(*,*) '  -  test_t_out.xxx.xxx.nc (full 3D)'
        if (nislice_count > 0) write(*,*) '  - test_islice.xxx.xxx.nc (', nislice_count, ' islices)'
        if (njslice_count > 0) write(*,*) '  - test_jslice.xxx.xxx.nc (', njslice_count, ' jslices)'
        if (nkslice_count > 0) write(*,*) '  - test_kslice.xxx.xxx.nc (', nkslice_count, ' kslices)'
      end if

    end subroutine tests_io


    subroutine write_islice_coords(ncid, local_nislice, nislice_total, islice_positions)
      use modglobal, only : xf, xh
      use netcdf
      implicit none
      integer, intent(in) :: ncid, local_nislice, nislice_total
      integer, dimension(nislice_total), intent(in) :: islice_positions
      integer :: varid, ierr, i, local_idx
      real, allocatable :: x_f(:), x_h(:)
      integer, allocatable :: indices(:)
      
      allocate(x_f(local_nislice))
      allocate(x_h(local_nislice))
      allocate(indices(local_nislice))
      
      local_idx = 0
      do i = 1, nislice_total
        if ( (islice(i)-1)/nprocx == myidx) then
          local_idx = local_idx + 1
          x_f(local_idx) = xf(islice_positions(i))
          x_h(local_idx) = xh(islice_positions(i))
          indices(local_idx) = islice_positions(i)
        end if
      end do
      
      ierr = nf90_inq_varid(ncid, 'xt', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'islice_indices', indices)
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, x_f)
      end if
      
      ierr = nf90_inq_varid(ncid, 'xm', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'islice_indices', indices)
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, x_h)
      end if
      
      deallocate(x_f, x_h, indices)
    end subroutine write_islice_coords


    subroutine write_kslice_coords(ncid, nkslice, kslice_positions)
      use modglobal, only : zf, zh
      use netcdf
      implicit none
      integer, intent(in) :: ncid, nkslice
      integer, dimension(nkslice), intent(in) :: kslice_positions
      integer :: varid, ierr, k
      real, allocatable :: z_f(:), z_h(:)
      
      allocate(z_f(nkslice))
      allocate(z_h(nkslice))
      
      do k = 1, nkslice
        z_f(k) = zf(kslice_positions(k))
        z_h(k) = zh(kslice_positions(k))
      end do
      
      ierr = nf90_inq_varid(ncid, 'zt', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'kslice_indices', kslice_positions)
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, z_f)
      end if
      
      ierr = nf90_inq_varid(ncid, 'zm', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'kslice_indices', kslice_positions)
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, z_h)
      end if
      
      deallocate(z_f, z_h)
    end subroutine write_kslice_coords


    subroutine write_jslice_coords(ncid, local_njslice, njslice_total, jslice_positions)
      use modglobal, only : yf, yh
      use netcdf
      implicit none
      integer, intent(in) :: ncid, local_njslice, njslice_total
      integer, dimension(njslice_total), intent(in) :: jslice_positions
      integer :: varid, ierr, j, local_idy
      real, allocatable :: y_f(:), y_h(:)
      integer, allocatable :: indices(:)
      
      allocate(y_f(local_njslice))
      allocate(y_h(local_njslice))
      allocate(indices(local_njslice))
      
      local_idy = 0
      do j = 1, njslice_total
        if ( (jslice(j)-1)/nprocy == myidy) then
          local_idy = local_idy + 1
          y_f(local_idy) = yf(jslice_positions(j))
          y_h(local_idy) = yh(jslice_positions(j))
          indices(local_idy) = jslice_positions(j)
        end if
      end do
      
      ierr = nf90_inq_varid(ncid, 'yt', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'jslice_indices', indices)
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, y_f)
      end if
      
      ierr = nf90_inq_varid(ncid, 'ym', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'jslice_indices', indices)
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, y_h)
      end if
      
      deallocate(y_f, y_h, indices)
    end subroutine write_jslice_coords


end module tests
