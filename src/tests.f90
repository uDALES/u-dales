module tests
  use MPI
  use decomp_2d
  use modglobal  ! Import all global variables for comprehensive output
  use modmpi, only : comm3d, my_real, mpierr, myid, myidx, myidy, myid1dy, nprocs, nprocx, nprocy, cmyidx, cmyidy
  use modstat_nc, only : open_nc, define_nc, writestat_dims_nc, ncinfo, writeoffset, writeoffset_1dx, writestat_nc
  use instant, only : write_islice_xcoord_local, write_jslice_ycoord_local, write_kslice_zcoord
  implicit none

  contains
    subroutine inittest
      integer :: nx=64, ny=64, nz=64
      !integer :: p_row=2, p_col=2
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


    subroutine tests_io
      use netcdf
      implicit none
      integer :: i, j, k, n, xdim, ydim, zdim, nrect, ncid, time_step
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

      ! Probe variables (New Implementation)
      integer :: nprobe_count
      integer :: ncid_probe, nrec_probe
      integer :: point_dimid, time_dimid, varid_x, varid_y, varid_z, varid_val, varid_time
      character(80) :: probe_filename
      real, allocatable :: probe_send_buf(:), probe_recv_buf(:)
      real, allocatable :: probe_x(:), probe_y(:), probe_z(:)
      integer :: pn
      integer :: gi, gj, gk, li, lj, lk
      integer :: target_idx, target_idy
      logical :: is_local

      integer :: ii_local, jj_local

      nrect = 0
      xdim = ie-ib+1
      ydim = je-jb+1
      zdim = ke-kb+1
      
      ! Get number of islices, jslices and kslices from modglobal
      nislice_count = nislice
      njslice_count = njslice
      nkslice_count = nkslice
      
      call MPI_BARRIER(comm3d, mpierr)

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
      end if     

! ============================================= !
! Create islice files
! ============================================= !
      if (nislice_count > 0) then
        if (myid == 0) write(*,*) 'Creating islice files for ', nislice_count, ' slices at x=', islice(1:nislice_count)
        
        ! Count local islices
        local_nislice = 0
        do i = 1, nislice_count
          !if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
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
            call write_islice_xcoord_local(ncid_islice, local_nislice, nislice_count, islice, myidx, nprocx, itot)
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
          !if ( (jslice(i)-1)/nprocy == myidy) then
          if (jslice(i) >= zstart(2) .and. jslice(i) <= zend(2)) then
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
            call write_jslice_ycoord_local(ncid_jslice, local_njslice, njslice_count, jslice, myidy, nprocy, jtot)
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
            call write_kslice_zcoord(ncid_kslice, nkslice_count, kslice)
          end if
          
          call ncinfo(varinfo(1,:), 'test_value', 'Test coordinate value', '1', 'tttt')
          call define_nc(ncid_kslice, 1, varinfo)
          
          write(*,'(A,I3,A,I2,A)') 'Processor myidx=', myidx, ' created kslice file with ', nkslice_count, ' slices'
        end if
      end if

! ============================================= !
! Create probe files (Master only)
! ============================================= !
      nprobe_count = nprobe
      if (nprobe_count > 0) then

        ! Only master creates the file
        if (myid == 0) then
           write(*,*) 'Creating single probe file for ', nprobe_count, ' probes'
           
           probe_filename = 'test_probe.xxx.nc'
           probe_filename(12:14) = cexpnr
           
           ! Create NetCDF file
           call check( nf90_create(trim(probe_filename), NF90_CLOBBER, ncid_probe) )
           
           ! Define Dimensions
           call check( nf90_def_dim(ncid_probe, 'point', nprobe_count, point_dimid) )
           call check( nf90_def_dim(ncid_probe, 'time', NF90_UNLIMITED, time_dimid) )
           
           ! Define Variables
           call check( nf90_def_var(ncid_probe, 'x', NF90_FLOAT, (/point_dimid/), varid_x) )
           call check( nf90_def_var(ncid_probe, 'y', NF90_FLOAT, (/point_dimid/), varid_y) )
           call check( nf90_def_var(ncid_probe, 'z', NF90_FLOAT, (/point_dimid/), varid_z) )
           call check( nf90_def_var(ncid_probe, 'time', NF90_FLOAT, (/time_dimid/), varid_time) )
           call check( nf90_def_var(ncid_probe, 'test_var', NF90_FLOAT, (/point_dimid, time_dimid/), varid_val) )
           
           call check( nf90_put_att(ncid_probe, varid_time, 'units', 'seconds since start') )
           call check( nf90_put_att(ncid_probe, varid_val, 'units', 'test_units') )

           call check( nf90_enddef(ncid_probe) )
           
           ! Write Coordinate Data (x, y, z)
           allocate(probe_x(nprobe_count))
           allocate(probe_y(nprobe_count))
           allocate(probe_z(nprobe_count))
           
           do n = 1, nprobe_count
              probe_x(n) = xf(iprobe(n))
              probe_y(n) = yf(jprobe(n))
              probe_z(n) = zf(kprobe(n))
           end do
           
           call check( nf90_put_var(ncid_probe, varid_x, probe_x) )
           call check( nf90_put_var(ncid_probe, varid_y, probe_y) )
           call check( nf90_put_var(ncid_probe, varid_z, probe_z) )
           
           deallocate(probe_x, probe_y, probe_z)
           
        endif
      endif
      
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
             !if ( (islice(i)-1)/(itot/nprocx) == myidx) then
             if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
               local_idx = local_idx + 1
               ii = islice(i)
               !ii_local = ii - (myidx * nprocx)
               ii_local = ii - zstart(1) + 1
                 islice_data(local_idx,:, :) = test_var(ii_local,jb:je, kb:ke)
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
             !if ( (jslice(i)-1)/(jtot/nprocy) == myidy) then
             if (jslice(i) >= zstart(2) .and. jslice(i) <= zend(2)) then
               local_idy = local_idy + 1
               jj = jslice(i)
               !jj_local = jj - (myidy * (jtot/nprocy))
               jj_local = jj - zstart(2) + 1
                 jslice_data(:,local_idy,:) = test_var(ib:ie, jj_local,kb:ke)
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
         
         ! Write probe data (New Implementation)
         if (nprobe_count > 0) then
         
           allocate(probe_send_buf(nprobe_count))
           probe_send_buf = 0.0
           
           do n = 1, nprobe_count
             gi = iprobe(n)
             gj = jprobe(n)
             gk = kprobe(n)
             
             ! Check if point is on this processor using zstart/zend
             is_local = .false.
             
             !target_idx = (gi - 1) / (itot / nprocx)
             !target_idy = (gj - 1) / (jtot / nprocy)
             !if (target_idx == myidx .and. target_idy == myidy) then
             !   if (gk >= 1 .and. gk <= ktot) then
             !       is_local = .true.
             !   endif
             !endif
             
             if (gi >= zstart(1) .and. gi <= zend(1)) then
                if (gj >= zstart(2) .and. gj <= zend(2)) then
                   if (gk >= zstart(3) .and. gk <= zend(3)) then
                       is_local = .true.
                   endif
                endif
             endif
             
             if (is_local) then
               !li = gi - (myidx * (itot / nprocx))
               !lj = gj - (myidy * (jtot / nprocy))
               !lk = gk
               li = gi - zstart(1) + 1
               lj = gj - zstart(2) + 1
               lk = gk - zstart(3) + 1
               
               probe_send_buf(n) = test_var(li, lj, lk)
               
             endif
           end do
           
           allocate(probe_recv_buf(nprobe_count))
           probe_recv_buf = 0.0 ! Initialize to zero
           
           ! Debug print before reduce
           if (myid == 0) then
              write(*,*) 'Calling MPI_REDUCE with nprobe_count=', nprobe_count
           endif
           
           call MPI_REDUCE(probe_send_buf, probe_recv_buf, nprobe_count, MY_REAL, MPI_SUM, 0, comm3d, mpierr)
           
           if (myid == 0) then
             ! Debug: Print received values
             write(*,*) 'Step', time_step, ' Probe values:', probe_recv_buf(1:nprobe_count)
           
             ! Write time
             call check( nf90_put_var(ncid_probe, varid_time, real(time_step), start=(/time_step/)) )
             
             ! Write probe values
             call check( nf90_put_var(ncid_probe, varid_val, probe_recv_buf, start=(/1, time_step/), count=(/nprobe_count, 1/)) )
           endif
           
           deallocate(probe_send_buf)
           deallocate(probe_recv_buf)
           
         endif
         
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
        if (nprobe_count > 0) write(*,*) '  - test_probe.nc (', nprobe_count, ' probes, single file)'
      end if
      
      ! Close Probe NC file
      if (myid == 0 .and. nprobe_count > 0) then
         call check( nf90_close(ncid_probe) )
      endif

    end subroutine tests_io

    subroutine check(status)
      use netcdf
      integer, intent ( in) :: status
      if(status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        stop 2
      end if
    end subroutine check

end module tests
