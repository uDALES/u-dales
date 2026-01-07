module tests
  use MPI
  use decomp_2d
  use modglobal  ! Import all global variables for comprehensive output
  use modmpi, only : comm3d, my_real, mpierr, myid, myidx, myidy, myid1dy, nprocs, nprocx, nprocy, cmyidx, cmyidy
  use modstat_nc, only : open_nc, define_nc, writestat_dims_nc, ncinfo, writeoffset, writeoffset_1dx, writestat_nc
  use instant_slice, only : write_islice_xcoord_local, write_jslice_ycoord_local, write_kslice_zcoord
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
      
      ! Probe variables
      integer :: nprobe_count
      integer :: ncid_probe, nrec_probe
      logical, allocatable :: probe_on_rank(:)
      integer, allocatable :: probe_local_i(:), probe_local_j(:)
      integer :: local_nprobe
      integer, allocatable :: local_probe_map(:)
      real, allocatable :: probe_data(:)
      character(80) :: probe_filename

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
          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
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
      
      call MPI_BARRIER(comm3d, mpierr)

! ============================================= !
! Create probe files
! ============================================= !
      nprobe_count = nprobe
      if (nprobe_count > 0) then
        if (myid == 0) write(*,*) 'Creating probe files for ', nprobe_count, ' probes'
        
        ! Allocate probe tracking arrays
        allocate(probe_on_rank(nprobe_count))
        allocate(probe_local_i(nprobe_count))
        allocate(probe_local_j(nprobe_count))
        
        ! Determine which probes are on this processor
        probe_on_rank = .false.
        local_nprobe = 0
        
        do n = 1, nprobe_count
          call get_probe_processor(iprobe(n), jprobe(n), probe_on_rank(n), &
                                   probe_local_i(n), probe_local_j(n))
          if (probe_on_rank(n)) then
            local_nprobe = local_nprobe + 1
          end if
        end do
        
        ! Create mapping from local to global probe index
        if (local_nprobe > 0) then
          allocate(local_probe_map(local_nprobe))
          allocate(probe_data(local_nprobe))
          local_idx = 0
          do n = 1, nprobe_count
            if (probe_on_rank(n)) then
              local_idx = local_idx + 1
              local_probe_map(local_idx) = n
            end if
          end do
          
          ! Create NetCDF file
          probe_filename = 'test_probe.xxx.xxx.nc'
          probe_filename(12:14) = cmyidx
          probe_filename(16:18) = cmyidy
          
          nrec_probe = 0
          call open_nc(probe_filename, ncid_probe, nrec_probe, n3=local_nprobe)
          
          if (nrec_probe == 0) then
            call ncinfo(timeVar(1,:), 'time', 'Time', 's', 'time')
            call define_nc(ncid_probe, 1, timeVar)
            call writestat_dims_nc(ncid_probe)
            call write_probe_coords(ncid_probe, local_nprobe, local_probe_map, nprobe_count)
          end if
          
          call ncinfo(varinfo(1,:), 'test_value', 'Test coordinate value', '1', 'tttt')
          call define_nc(ncid_probe, 1, varinfo)
          
          write(*,'(A,I3,A,I3,A,I2,A)') 'Processor ', myid, ' created probe file with ', &
                                         local_nprobe, ' of ', nprobe_count, ' probes'
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
             if ( (islice(i)-1)/(itot/nprocx) == myidx) then
               local_idx = local_idx + 1
               ii = islice(i)
               ii_local = ii - (myidx * nprocx)
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
             if ( (jslice(i)-1)/(jtot/nprocy) == myidy) then
               local_idy = local_idy + 1
               jj = jslice(i)
               jj_local = jj - (myidy * (jtot/nprocy))
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
         
         ! Write probe data
         if (nprobe_count > 0 .and. local_nprobe > 0) then
           do local_idx = 1, local_nprobe
             n = local_probe_map(local_idx)
             i = probe_local_i(n)
             j = probe_local_j(n)
             k = kprobe(n)
             probe_data(local_idx) = test_var(i, j, k)
           end do
           
           nrec_probe = time_step
           call writestat_nc(ncid_probe, 1, nrec_probe, varinfo(1,:), probe_data, local_nprobe)
         end if
         
         call MPI_BARRIER(comm3d, mpierr)
      end do
      
      ! Cleanup
      deallocate(test_var)
      if (allocated(islice_data)) deallocate(islice_data)
      if (allocated(jslice_data)) deallocate(jslice_data)
      if (allocated(kslice_data)) deallocate(kslice_data)
      if (allocated(probe_on_rank)) deallocate(probe_on_rank)
      if (allocated(probe_local_i)) deallocate(probe_local_i)
      if (allocated(probe_local_j)) deallocate(probe_local_j)
      if (allocated(local_probe_map)) deallocate(local_probe_map)
      if (allocated(probe_data)) deallocate(probe_data)
      
      if (myid == 0) then
        write(*,*) 'Test completed! Files written:'
        write(*,*) '  -  test_t_out.xxx.xxx.nc (full 3D)'
        if (nislice_count > 0) write(*,*) '  - test_islice.xxx.xxx.nc (', nislice_count, ' islices)'
        if (njslice_count > 0) write(*,*) '  - test_jslice.xxx.xxx.nc (', njslice_count, ' jslices)'
        if (nkslice_count > 0) write(*,*) '  - test_kslice.xxx.xxx.nc (', nkslice_count, ' kslices)'
        if (nprobe_count > 0) write(*,*) '  - test_probe.xxx.xxx.nc (', nprobe_count, ' probes)'
      end if

    end subroutine tests_io
    
    subroutine get_probe_processor(ig, jg, on_rank, local_i, local_j)
      ! Determine if a probe point is on this processor
      ! and calculate local indices
      implicit none
      integer, intent(in)  :: ig, jg  ! Global indices
      logical, intent(out) :: on_rank
      integer, intent(out) :: local_i, local_j
      
      integer :: istart, iend, jstart, jend
      
      ! Calculate global index range for this processor
      istart = myidx * (itot / nprocx) + 1
      iend = (myidx + 1) * (itot / nprocx)
      jstart = myidy * (jtot / nprocy) + 1
      jend = (myidy + 1) * (jtot / nprocy)
      
      ! Check if probe is within this processor's domain
      on_rank = (ig >= istart .and. ig <= iend .and. jg >= jstart .and. jg <= jend)
      
      if (on_rank) then
        local_i = ig - istart + ib
        local_j = jg - jstart + jb
      else
        local_i = -1
        local_j = -1
      end if
      
    end subroutine get_probe_processor
    
    subroutine write_probe_coords(ncid, local_nprobe, local_probe_map, nprobe_total)
      ! Write probe coordinates as variable attributes (only for local probes)
      use netcdf
      implicit none
      integer, intent(in) :: ncid, local_nprobe, nprobe_total
      integer, intent(in) :: local_probe_map(local_nprobe)
      
      integer :: varid, ierr, n, local_idx
      integer :: local_iprobe(local_nprobe), local_jprobe(local_nprobe), local_kprobe(local_nprobe)
      
      ! Get local probe coordinates from global arrays
      do local_idx = 1, local_nprobe
        n = local_probe_map(local_idx)
        local_iprobe(local_idx) = iprobe(n)
        local_jprobe(local_idx) = jprobe(n)
        local_kprobe(local_idx) = kprobe(n)
      end do
      
      ! Get variable ID for test_value
      ierr = nf90_inq_varid(ncid, 'test_value', varid)
      
      ! Write coordinates as attributes
      ierr = nf90_put_att(ncid, varid, 'probe_i_indices', local_iprobe)
      ierr = nf90_put_att(ncid, varid, 'probe_j_indices', local_jprobe)
      ierr = nf90_put_att(ncid, varid, 'probe_k_indices', local_kprobe)
      
    end subroutine write_probe_coords


end module tests
