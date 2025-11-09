module instant_slice
  use modglobal,  only : cexpnr, rk3step, ltempeq, lmoist, nsv, tsample, dt, timee, runtime, &
                         slicevars, lislicedump, islice, ljslicedump, jslice, lkslicedump, kslice, &
                         nkslice, nislice, njslice, &
                         ib, ie, jb, je, kb, ke, dzfi, dzh, &
                         timee, tstatstart, jtot, itot, kmax
  use modfields,  only : um, vm, wm, thlm, qtm, svm
  use modmpi,     only : myid, cmyidx, cmyidy
  use decomp_2d,  only : zstart, zend, xstart, xend, ystart, yend
  use modstat_nc, only : ncinfo, open_nc, define_nc, writestat_dims_nc, writestat_nc, writeoffset, writeoffset_1dx
  implicit none
  private
  public :: instant_init, instant_main
  save

  integer :: xdim, ydim, zdim, kdim, idim, jdim  ! Added kdim, idim, jdim for multiple slices
  real    :: tsampleslice

  logical, allocatable :: islicerank(:)   ! array of flags for each islice on this core
  integer, allocatable :: isliceloc(:)    ! local islice positions on this core
  logical, allocatable :: jslicerank(:)   ! array of flags for each jslice on this core
  integer, allocatable :: jsliceloc(:)    ! local jslice positions on this core
  
  integer :: local_nislice, local_njslice  ! Number of islices on this X-processor (saved from create phase)
  integer, allocatable :: local_islice_map(:)  ! Mapping from local index to global islice index
  
  integer                    :: nslicevars
  character(80)              :: islicename
  character(80)              :: jslicename
  character(80)              :: kslicename
  character(80)              :: kslicename1d    ! Separate filename for 1D output
  
  character(80)              :: slicetimeVar(1,4)
  character(80), allocatable :: isliceVars(:,:)
  character(80), allocatable :: jsliceVars(:,:)
  character(80), allocatable :: ksliceVars(:,:)
  integer                    :: ncidislice, nrecislice, ncidjslice, nrecjslice, ncidkslice, nreckslice
  integer                    :: ncidkslice1d, nreckslice1d    ! Separate file ID for 1D output

  contains

    subroutine instant_init
      implicit none

      if(.not.(lislicedump .or. ljslicedump .or. lkslicedump)) return

      nslicevars = (LEN(trim(slicevars))+1)/3

      if (nslicevars == 0) then
        lislicedump = .false.
        ljslicedump = .false.
        lkslicedump = .false.
        print *, "NOTE: invalid 'slicevars' therefore all of lislicedump, ljslicedump &
                  &and lkslicedump are set to false, and no instantaneous slice will be outputted !!"
        return
      end if

      if(runtime <= tstatstart) then
        if(myid==0) then
          write(*,*) "ERROR: no instantaneous slice will be written as runtime <= tstatstart. Note that runtime &
                      &must be greater than tstatstart for wiriting slice files."
          write(*,*) "You have used runtime = ", runtime, ", tstatstart = ", tstatstart
          write(*,*) "Either correct the time settings or change all the slice writing flags to false."
          stop 1
        end if
      end if
      
      xdim = ie-ib+1
      ydim = je-jb+1
      zdim = ke-kb+1
      kdim = nkslice  ! Set kdim to number of kslices
      idim = nislice  ! Set idim to number of islices
      jdim = njslice  ! Set jdim to number of jslices
      tsampleslice = 0.
      call ncinfo(slicetimeVar( 1,:), 'time', 'Time', 's', 'time')

      if (lislicedump) then
        allocate(isliceVars(nslicevars,4))
        allocate(islicerank(nislice))
        allocate(isliceloc(nislice))
        allocate(local_islice_map(nislice))  ! Allocate mapping array
        call instant_ncdescription_islice
        call instant_create_ncislice    !> Generate sliced NetCDF: islice.xxx.xxx.nc
        deallocate(isliceVars)
      end if
      
      if (ljslicedump) then
        allocate(jsliceVars(nslicevars,4))
        allocate(jslicerank(njslice))
        allocate(jsliceloc(njslice))
        call instant_ncdescription_jslice
        call instant_create_ncjslice    ! Unified jslice creation (per myidy file, myidx==0)
        deallocate(jsliceVars)
      end if

      if (lkslicedump) then
        allocate(ksliceVars(nslicevars,4))
        call instant_ncdescription_kslice
        call instant_create_nckslice1d    !> Generate sliced NetCDF: kslice.xxx.xxx.xxx.nc
        deallocate(ksliceVars)
      end if
    end subroutine instant_init


    subroutine instant_main
      implicit none
      if (timee < tstatstart) return
      if (.not. rk3step==3)  return
      if(.not.(lislicedump .or. ljslicedump .or. lkslicedump)) return

      if (tsampleslice > tsample) then
        if (lislicedump) call instant_write_islice
        if (ljslicedump) call instant_write_jslice
        if (lkslicedump) call instant_write_kslice
        tsampleslice = dt
      else
        tsampleslice = tsampleslice + dt
      endif
    end subroutine instant_main


    subroutine instant_ncdescription_islice
      implicit none
      integer :: n
      do n=1,nslicevars
        select case(slicevars(3*n-2:3*n-1))
          case('u0')
            call ncinfo( isliceVars(n,:), 'u' , 'Streamwise velocity' , 'm/s' , 'tmtt' )
          case('v0')
            call ncinfo( isliceVars(n,:), 'v' , 'Spanwise velocity'   , 'm/s' , 'tttt' )
          case('w0')
            call ncinfo( isliceVars(n,:), 'w' , 'Vertical velocity'   , 'm/s' , 'ttmt' )
          case('th')
            if (ltempeq) call ncinfo( isliceVars(n,:), 'thl' , 'Potential temperature' , 'K'     , 'tttt' )
          case('qt')
            if (lmoist)  call ncinfo( isliceVars(n,:), 'qt'  , 'Specific humidity'     , 'kg/kg' , 'tttt' )
          case('s1')
            if (nsv>0)   call ncinfo( isliceVars(n,:), 's1'  , 'Concentration field 1' , 'g/m^3' , 'tttt' )
          case('s2')
            if (nsv>1)   call ncinfo( isliceVars(n,:), 's2'  , 'Concentration field 2' , 'g/m^3' , 'tttt' )
          case('s3')
            if (nsv>2)   call ncinfo( isliceVars(n,:), 's3'  , 'Concentration field 3' , 'g/m^3' , 'tttt' )
          case('s4')
            if (nsv>3)   call ncinfo( isliceVars(n,:), 's4'  , 'Concentration field 4' , 'g/m^3' , 'tttt' )
          case default
            print *, "Invalid slice variables name. Check namoptions variable 'slicevars'. &
                      &There should not be any space in the string."
            STOP 1
        end select
      end do
    end subroutine instant_ncdescription_islice

    subroutine instant_ncdescription_jslice
      implicit none
      integer :: n
      do n=1,nslicevars
        select case(slicevars(3*n-2:3*n-1))
          case('u0')
            call ncinfo( jsliceVars(n,:), 'u' , 'Streamwise velocity' , 'm/s' , 'mttt' )
          case('v0')
            call ncinfo( jsliceVars(n,:), 'v' , 'Spanwise velocity '  , 'm/s' , 'tmtt' )
          case('w0')
            call ncinfo( jsliceVars(n,:), 'w' , 'Vertical velocity'   , 'm/s' , 'ttmt' )
          case('th')
            if (ltempeq) call ncinfo( jsliceVars(n,:), 'thl' , 'Potential temperature' , 'K'     , 'tttt' )
          case('qt')
            if (lmoist)  call ncinfo( jsliceVars(n,:), 'qt'  , 'Specific humidity'     , 'kg/kg' , 'tttt' )
          case('s1')
            if (nsv>0)   call ncinfo( jsliceVars(n,:), 's1'  , 'Concentration field 1' , 'g/m^3' , 'tttt' )
          case('s2')
            if (nsv>1)   call ncinfo( jsliceVars(n,:), 's2'  , 'Concentration field 2' , 'g/m^3' , 'tttt' )
          case('s3')
            if (nsv>2)   call ncinfo( jsliceVars(n,:), 's3'  , 'Concentration field 3' , 'g/m^3' , 'tttt' )
          case('s4')
            if (nsv>3)   call ncinfo( jsliceVars(n,:), 's4'  , 'Concentration field 4' , 'g/m^3' , 'tttt' )
          case default
            print *, "Invalid slice variables name. Check namoptions variable 'slicevars'. &
                      &There should not be any space in the string."
            STOP 1
        end select
      end do
    end subroutine instant_ncdescription_jslice

    subroutine instant_ncdescription_kslice
      implicit none
      integer :: n
      do n=1,nslicevars
        select case(slicevars(3*n-2:3*n-1))
          case('u0')
            call ncinfo( ksliceVars(n,:), 'u' , 'Streamwise velocity' , 'm/s' , 'mttt' )
          case('v0')
            call ncinfo( ksliceVars(n,:), 'v' , 'Spanwise velocity'   , 'm/s' , 'tmtt' )
          case('w0')
            call ncinfo( ksliceVars(n,:), 'w' , 'Vertical velocity'   , 'm/s' , 'ttmt' )
          case('th')
            if (ltempeq) call ncinfo( ksliceVars( n,:), 'thl' , 'Potential temperature' , 'K'     , 'tttt' )
          case('qt')
            if (lmoist)  call ncinfo( ksliceVars( n,:), 'qt'  , 'Specific humidity'     , 'kg/kg' , 'tttt' )
          case('s1')
            if (nsv>0)   call ncinfo( ksliceVars( n,:), 's1'  , 'Concentration field 1' , 'g/m^3' , 'tttt' )
          case('s2')
            if (nsv>1)   call ncinfo( ksliceVars( n,:), 's2'  , 'Concentration field 2' , 'g/m^3' , 'tttt' )
          case('s3')
            if (nsv>2)   call ncinfo( ksliceVars( n,:), 's3'  , 'Concentration field 3' , 'g/m^3' , 'tttt' )
          case('s4')
            if (nsv>3)   call ncinfo( ksliceVars( n,:), 's4'  , 'Concentration field 4' , 'g/m^3' , 'tttt' )
          case default
            print *, "Invalid slice variables name. Check namoptions variable 'slicevars'. &
                      &There should not be any space in the string."
            STOP 1
        end select
      end do
    end subroutine instant_ncdescription_kslice


    !! ## %% 1D parallel output creation for islice (y-direction processes only)
    subroutine instant_create_ncislice
      use modglobal, only : jtot
      use modmpi, only : myidy, myidx, nprocx
      implicit none
      integer :: i
      logical :: has_islice
      
      if (nislice == 0) then
        if (myid == 0) write(*,*) "WARNING: nislice=0, no i-slices will be output"
        return
      end if
      
      ! Count how many islices are in this X-processor's range (save to module variable)
      local_nislice = 0
      has_islice = .false.
      do i = 1, nislice
        if ( (islice(i)-1)/nprocx == myidx) then
          local_nislice = local_nislice + 1
          has_islice = .true.
        end if
      end do
      
      ! Only processors with islices AND myidy==0 create files
      if (has_islice .and. myidy == 0) then
        islicename = 'islice.xxx.xxx.nc'
        islicename(8:10) = cmyidx   ! X-processor ID
        islicename(12:14) = cexpnr  ! Experiment number

        nrecislice = 0
        ! Use local_nislice (module variable) as first dimension, jtot as global y-dimension
        call open_nc(islicename, ncidislice, nrecislice, n1=local_nislice, n2=jtot, n3=zdim)
        if (nrecislice == 0) then
          call define_nc(ncidislice, 1, slicetimeVar)
          call writestat_dims_nc(ncidislice)
          ! Add islice-specific x coordinate information
          call write_islice_xcoord_local
        end if
        call define_nc(ncidislice, nslicevars, isliceVars)
        
        write(*,'(A,A,A,I2,A)') '  Processor (myidx=', cmyidx, ') created islice file with ', local_nislice, ' slices'
      end if
    end subroutine instant_create_ncislice

    subroutine instant_create_ncjslice
      use modglobal, only : itot
      use modmpi, only : myidx, myidy, nprocy
      implicit none
      integer :: j
      logical :: has_jslice
      
      if (njslice == 0) then
        if (myid == 0) write(*,*) "WARNING: njslice=0, no j-slices will be output"
        return
      end if

      ! Count how many jslices are assigned to this Y-processor (round-robin by myidy)
      local_njslice = 0
      has_jslice = .false.
      do j = 1, njslice
        if ( (jslice(j)-1)/nprocy == myidy) then
          local_njslice = local_njslice + 1
          has_jslice = .true.
        end if
      end do
      
      ! Only processors with myidx==0 create files for their myidy column
      if (has_jslice .and. myidx == 0) then
        jslicename = 'jslice.xxx.xxx.nc'
        jslicename(8:10) = cmyidy   ! Y-processor ID
        jslicename(12:14) = cexpnr  ! Experiment number
        
        nrecjslice = 0
        ! Use itot as x-dimension (global), local_n as jslice-dimension, zdim as z
        call open_nc(jslicename, ncidjslice, nrecjslice, n1=itot, n2=local_njslice, n3=zdim)
        if (nrecjslice==0) then
          call define_nc(ncidjslice, 1, slicetimeVar)
          call writestat_dims_nc(ncidjslice)
          call write_jslice_ycoord_local
        end if
        call define_nc(ncidjslice, nslicevars, jsliceVars)
        write(*,'(A,A,A,I2,A)') '  Processor (myidy=', cmyidy, ') created islice file with ', local_njslice, ' slices'
      end if
    end subroutine instant_create_ncjslice

    subroutine instant_create_nckslice
      implicit none
      
      if (nkslice == 0) then
        if (myid == 0) write(*,*) "WARNING: nkslice=0, no k-slices will be output"
        return
      end if
      
      kslicename = 'kslice.xxx.xxx.xxx.nc'
      kslicename(8:10)  = cmyidx
      kslicename(12:14) = cmyidy
      kslicename(16:18) = cexpnr

      nreckslice = 0
      ! Use kdim (nkslice) as the z-dimension instead of default
      call open_nc(kslicename, ncidkslice, nreckslice, n1=xdim, n2=ydim, n3=kdim)
      if (nreckslice==0) then
        call define_nc(ncidkslice, 1, slicetimeVar)
        ! Use standard writestat_dims_nc - it will write x, y, z dimensions
        call writestat_dims_nc(ncidkslice)
        ! Now add kslice-specific z coordinate information
        call write_kslice_zcoord(ncidkslice)
      end if
      call define_nc(ncidkslice, nslicevars, ksliceVars)
    end subroutine instant_create_nckslice

    !! ## %% 1D parallel output creation (y-direction processes only)
    subroutine instant_create_nckslice1d
      use modglobal, only : jtot
      use modmpi, only : myidy
      implicit none
      
      if (nkslice == 0) then
        if (myid == 0) write(*,*) "WARNING: nkslice=0, no k-slices will be output"
        return
      end if

      if (myidy == 0) then
        kslicename1d = 'kslice.xxx.xxx.nc'
        kslicename1d(8:10) = cmyidx   ! Only x-processor ID
        kslicename1d(12:14) = cexpnr  ! Experiment number

        nreckslice1d = 0
        ! Use kdim (nkslice) as the z-dimension, jtot as global y-dimension
        call open_nc(kslicename1d, ncidkslice1d, nreckslice1d, n1=xdim, n2=jtot, n3=kdim)
        if (nreckslice1d == 0) then
          call define_nc(ncidkslice1d, 1, slicetimeVar)
          call writestat_dims_nc(ncidkslice1d)
          ! Add kslice-specific z coordinate information
          call write_kslice_zcoord(ncidkslice1d)
        end if
        call define_nc(ncidkslice1d, nslicevars, ksliceVars)
      end if
    end subroutine instant_create_nckslice1d

    subroutine instant_write_islice
      use modmpi, only : myidy, myidx, nprocx
      use modglobal, only : jtot
      implicit none
      real, allocatable :: tmp_slice(:,:,:)
      integer :: i, ii, ii_local, local_idx
      
      if (nislice == 0) return
      if (local_nislice == 0) return  ! No islices on this processor
      
      ! Allocate temporary array with local_nislice (module variable)
      allocate(tmp_slice(local_nislice, jb:je, kb:ke))
      tmp_slice = 0.0
      
      ! Write time variable (only myidy==0 writes)
      if (myidy == 0) then
        call writestat_nc(ncidislice, 'time', timee, nrecislice, .true.)
      end if
      
      ! u velocity (interpolated to cell centers in x-direction)
      if (present('u0')) then
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/nprocx == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)  ! Global X-position
            ii_local = ii - xstart(1) + ib
            if (ii_local >= ib .and. ii_local < ie) then
              tmp_slice(local_idx,:,:) = 0.5*(um(ii_local,:,:) + um(ii_local+1,:,:))
            end if
          end if
        end do
        call writeoffset(ncidislice, 'u', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      ! v velocity
      if (present('v0')) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/nprocx == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - xstart(1) + ib
            if (ii_local >= ib .and. ii_local <= ie) then
              tmp_slice(local_idx,:,:) = vm(ii_local,:,:)
            end if
          end if
        end do
        call writeoffset(ncidislice, 'v', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      ! w velocity
      if (present('w0')) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/nprocx == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - xstart(1) + ib
            if (ii_local >= ib .and. ii_local <= ie) then
              tmp_slice(local_idx,:,:) = wm(ii_local,:,:)
            end if
          end if
        end do
        call writeoffset(ncidislice, 'w', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      ! temperature
      if (present('th') .and. ltempeq) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/nprocx == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - xstart(1) + ib
            if (ii_local >= ib .and. ii_local <= ie) then
              tmp_slice(local_idx,:,:) = thlm(ii_local,:,:)
            end if
          end if
        end do
        call writeoffset(ncidislice, 'thl', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      ! moisture
      if (present('qt') .and. lmoist) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/nprocx == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - xstart(1) + ib
            if (ii_local >= ib .and. ii_local <= ie) then
              tmp_slice(local_idx,:,:) = qtm(ii_local,:,:)
            end if
          end if
        end do
        call writeoffset(ncidislice, 'qt', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      ! scalars s1-s4
      if (present('s1') .and. nsv>0) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/nprocx == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - xstart(1) + ib
            if (ii_local >= ib .and. ii_local <= ie) then
              tmp_slice(local_idx,:,:) = svm(ii_local,:,:,1)
            end if
          end if
        end do
        call writeoffset(ncidislice, 's1', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      if (present('s2') .and. nsv>1) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/nprocx == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - xstart(1) + ib
            if (ii_local >= ib .and. ii_local <= ie) then
              tmp_slice(local_idx,:,:) = svm(ii_local,:,:,2)
            end if
          end if
        end do
        call writeoffset(ncidislice, 's2', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      if (present('s3') .and. nsv>2) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/nprocx == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - xstart(1) + ib
            if (ii_local >= ib .and. ii_local <= ie) then
              tmp_slice(local_idx,:,:) = svm(ii_local,:,:,3)
            end if
          end if
        end do
        call writeoffset(ncidislice, 's3', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      if (present('s4') .and. nsv>3) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/nprocx == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - xstart(1) + ib
            if (ii_local >= ib .and. ii_local <= ie) then
              tmp_slice(local_idx,:,:) = svm(ii_local,:,:,4)
            end if
          end if
        end do
        call writeoffset(ncidislice, 's4', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      deallocate(tmp_slice)
    end subroutine instant_write_islice

    subroutine instant_write_jslice
      use modmpi, only : myidx, myidy, nprocy
      use modglobal, only : itot
      implicit none
      real, allocatable :: tmp_slice(:,:,:)
      integer :: j, jj, jj_local, local_idy
      
      if (njslice == 0) return
      if (local_njslice == 0) return

      allocate(tmp_slice(ib:ie, local_njslice, kb:ke))
      tmp_slice = 0.0

      ! Write time variable (only myidx==0 writes)
      if (myidx == 0) then
        call writestat_nc(ncidjslice, 'time', timee, nrecjslice, .true.)
      end if

      ! u velocity
      if (present('u0')) then
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/nprocy == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - ystart(2) + jb
            if (jj_local >= jb .and. jj_local <= je) then
              tmp_slice(:, local_idy, :) = um(ib:ie, jj_local, kb:ke)
            end if
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'u', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      ! v velocity (interpolated to cell centers in y-direction)
      if (present('v0')) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/nprocy == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - ystart(2) + jb
            if (jj_local >= jb .and. jj_local < je) then
              tmp_slice(:, local_idy, :) = 0.5*(vm(ib:ie, jj_local, kb:ke) + vm(ib:ie, jj_local+1, kb:ke))
            end if
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'v', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      ! w velocity
      if (present('w0')) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/nprocy == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - ystart(2) + jb
            if (jj_local >= jb .and. jj_local <= je) then
              tmp_slice(:, local_idy, :) = wm(ib:ie, jj_local, kb:ke)
            end if
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'w', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      ! temperature
      if (present('th') .and. ltempeq) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/nprocy == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - ystart(2) + jb
            if (jj_local >= jb .and. jj_local <= je) then
              tmp_slice(:, local_idy, :) = thlm(ib:ie, jj_local, kb:ke)
            end if
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'thl', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      ! moisture
      if (present('qt') .and. lmoist) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/nprocy == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - ystart(2) + jb
            if (jj_local >= jb .and. jj_local <= je) then
              tmp_slice(:, local_idy, :) = qtm(ib:ie, jj_local, kb:ke)
            end if
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'qt', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      ! scalars s1-s4
      if (present('s1') .and. nsv>0) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/nprocy == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - ystart(2) + jb
            if (jj_local >= jb .and. jj_local <= je) then
              tmp_slice(:, local_idy, :) = svm(ib:ie, jj_local, kb:ke, 1)
            end if
          end if
        end do
        call writeoffset_1dx(ncidjslice, 's1', tmp_slice, nrecjslice, xdim,  local_njslice, zdim)
      end if

      if (present('s2') .and. nsv>1) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/nprocy == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - ystart(2) + jb
            if (jj_local >= jb .and. jj_local <= je) then
              tmp_slice(:, local_idy, :) = svm(ib:ie, jj_local, kb:ke, 2)
            end if
          end if
        end do
        call writeoffset_1dx(ncidjslice, 's2', tmp_slice, nrecjslice, xdim,  local_njslice, zdim)
      end if

      if (present('s3') .and. nsv>2) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/nprocy == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - ystart(2) + jb
            if (jj_local >= jb .and. jj_local <= je) then
              tmp_slice(:, local_idy, :) = svm(ib:ie, jj_local, kb:ke, 3)
            end if
          end if
        end do
        call writeoffset_1dx(ncidjslice, 's3', tmp_slice, nrecjslice, xdim,  local_njslice, zdim)
      end if

      if (present('s4') .and. nsv>3) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/nprocy == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - ystart(2) + jb
            if (jj_local >= jb .and. jj_local <= je) then
              tmp_slice(:, local_idy, :) = svm(ib:ie, jj_local, kb:ke, 4)
            end if
          end if
        end do
        call writeoffset_1dx(ncidjslice, 's4', tmp_slice, nrecjslice, xdim,  local_njslice, zdim)
      end if

      deallocate(tmp_slice)
    end subroutine instant_write_jslice

    subroutine write_islice_xcoord_local
      ! Write x-coordinates for LOCAL islice positions only
      use modglobal, only : xf, xh
      use modmpi, only : myidx, nprocx
      use netcdf
      implicit none
      integer :: varid, ierr, i, local_idx
      real, allocatable :: x_islice_f(:), x_islice_h(:)
      integer, allocatable :: islice_indices(:)
      
      ! Allocate and fill x coordinates for LOCAL islices only (use module variable)
      allocate(x_islice_f(local_nislice))
      allocate(x_islice_h(local_nislice))
      allocate(islice_indices(local_nislice))
      
      local_idx = 0
      do i = 1, nislice
        if ( (islice(i)-1)/nprocx == myidx) then
          local_idx = local_idx + 1
          x_islice_f(local_idx) = xf(islice(i))  ! full level (cell center)
          x_islice_h(local_idx) = xh(islice(i))  ! half level (cell edge)
          islice_indices(local_idx) = islice(i)
        end if
      end do
      
      ! Write xt (full level)
      ierr = nf90_inq_varid(ncidislice, 'xt', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncidislice)
        ierr = nf90_put_att(ncidislice, varid, 'islice_indices', islice_indices)
        ierr = nf90_put_att(ncidislice, varid, 'long_name', 'x-coordinate of i-slices (cell center)')
        ierr = nf90_enddef(ncidislice)
        ierr = nf90_put_var(ncidislice, varid, x_islice_f)
      end if
      
      ! Write xm (half level)
      ierr = nf90_inq_varid(ncidislice, 'xm', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncidislice)
        ierr = nf90_put_att(ncidislice, varid, 'islice_indices', islice_indices)
        ierr = nf90_put_att(ncidislice, varid, 'long_name', 'x-coordinate of i-slices (cell edge)')
        ierr = nf90_enddef(ncidislice)
        ierr = nf90_put_var(ncidislice, varid, x_islice_h)
      end if
      
      deallocate(x_islice_f)
      deallocate(x_islice_h)
      deallocate(islice_indices)
    end subroutine write_islice_xcoord_local

    subroutine write_jslice_ycoord_local
      ! Write y-coordinates for LOCAL jslice positions (round-robin assignment over myidy)
      use modglobal, only : yf, yh
      use modmpi, only : myidy, nprocy
      use netcdf
      implicit none
      integer :: varid, ierr, j, local_idy
      real, allocatable :: y_jslice_f(:), y_jslice_h(:)
      integer, allocatable :: jslice_indices(:)

      allocate(y_jslice_f(local_njslice))
      allocate(y_jslice_h(local_njslice))
      allocate(jslice_indices(local_njslice))
      
      local_idy = 0
      do j = 1, njslice
        if ( (jslice(j)-1)/nprocy == myidy) then
          local_idy = local_idy + 1
          y_jslice_f(local_idy) = yf(jslice(j))
          y_jslice_h(local_idy) = yh(jslice(j))
          jslice_indices(local_idy) = jslice(j)
        end if
      end do

      ierr = nf90_inq_varid(ncidjslice, 'yt', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncidjslice)
        ierr = nf90_put_att(ncidjslice, varid, 'jslice_indices', jslice_indices)
        ierr = nf90_put_att(ncidjslice, varid, 'long_name', 'y-coordinate of j-slices (cell center)')
        ierr = nf90_enddef(ncidjslice)
        ierr = nf90_put_var(ncidjslice, varid, y_jslice_f)
      end if

      ierr = nf90_inq_varid(ncidjslice, 'ym', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncidjslice)
        ierr = nf90_put_att(ncidjslice, varid, 'jslice_indices', jslice_indices)
        ierr = nf90_put_att(ncidjslice, varid, 'long_name', 'y-coordinate of j-slices (cell edge)')
        ierr = nf90_enddef(ncidjslice)
        ierr = nf90_put_var(ncidjslice, varid, y_jslice_h)
      end if

      deallocate(y_jslice_f)
      deallocate(y_jslice_h)
      deallocate(jslice_indices)
    end subroutine write_jslice_ycoord_local

    subroutine write_kslice_zcoord(ncid)
      ! Update z-coordinate to reflect kslice levels instead of full z levels
      use modglobal, only : zf, zh
      use netcdf
      implicit none
      integer, intent(in) :: ncid
      integer :: varid, ierr, k
      real, allocatable :: z_kslice_f(:), z_kslice_h(:)
      
      ! Allocate and fill z coordinates for kslices
      allocate(z_kslice_f(nkslice))
      allocate(z_kslice_h(nkslice))
      do k = 1, nkslice
        z_kslice_f(k) = zf(kslice(k))  ! full level (cell center)
        z_kslice_h(k) = zh(kslice(k))  ! half level (cell edge)
      end do
      
      ! Write zt (full level)
      ierr = nf90_inq_varid(ncid, 'zt', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'kslice_indices', kslice(1:nkslice))
        ierr = nf90_put_att(ncid, varid, 'long_name', 'z-coordinate of k-slices (cell center)')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, z_kslice_f)
        
        if (ierr == nf90_noerr .and. myid == 0) then
          write(*,*) 'Updated zt for kslices:', z_kslice_f
        end if
      end if
      
      ! Write zm (half level)
      ierr = nf90_inq_varid(ncid, 'zm', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'kslice_indices', kslice(1:nkslice))
        ierr = nf90_put_att(ncid, varid, 'long_name', 'z-coordinate of k-slices (cell edge)')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, z_kslice_h)
        
        if (ierr == nf90_noerr .and. myid == 0) then
          write(*,*) 'Updated zm for kslices:', z_kslice_h
        end if
      end if
      
      deallocate(z_kslice_f)
      deallocate(z_kslice_h)
    end subroutine write_kslice_zcoord

    subroutine instant_write_kslice
      use modmpi, only : myidy
      use modglobal, only : jtot
      implicit none
      real, allocatable :: tmp_slice(:,:,:)
      integer :: k, kk
      
      if (nkslice == 0) return
      
      ! Allocate temporary array to hold all kslice levels: (x, y, nkslice)
      allocate(tmp_slice(ib:ie, jb:je, nkslice))
      
      ! Write time variable (only myidy==0 writes)
      if (myidy == 0) then
        call writestat_nc(ncidkslice1d, 'time', timee, nreckslice1d, .true.)
      end if
      
      ! u velocity
      if (present('u0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = um(ib:ie, jb:je, kk)
        end do
        call writeoffset(ncidkslice1d, 'u', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if
      
      ! v velocity
      if (present('v0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = vm(ib:ie, jb:je, kk)
        end do
        call writeoffset(ncidkslice1d, 'v', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if
      
      ! w velocity (interpolated to cell centers)
      if (present('w0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = 0.5 * (wm(ib:ie, jb:je, kk) + wm(ib:ie, jb:je, kk+1))
        end do
        call writeoffset(ncidkslice1d, 'w', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if
      
      ! temperature
      if (present('th') .and. ltempeq) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = thlm(ib:ie, jb:je, kk)
        end do
        call writeoffset(ncidkslice1d, 'thl', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if
      
      ! moisture
      if (present('qt') .and. lmoist) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = qtm(ib:ie, jb:je, kk)
        end do
        call writeoffset(ncidkslice1d, 'qt', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if
      
      ! scalars s1-s4
      if (present('s1') .and. nsv>0) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie, jb:je, kk, 1)
        end do
        call writeoffset(ncidkslice1d, 's1', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if
      
      if (present('s2') .and. nsv>1) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie, jb:je, kk, 2)
        end do
        call writeoffset(ncidkslice1d, 's2', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if
      
      if (present('s3') .and. nsv>2) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie, jb:je, kk, 3)
        end do
        call writeoffset(ncidkslice1d, 's3', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if
      
      if (present('s4') .and. nsv>3) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie, jb:je, kk, 4)
        end do
        call writeoffset(ncidkslice1d, 's4', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if
      
      deallocate(tmp_slice)
    end subroutine instant_write_kslice

    logical function present(str)
      character(len=*), intent(in) :: str
      integer :: pos

      pos = index(slicevars, str)   ! returns the position of the substring str in slicevars. If not found, it returns 0.
      present = (pos > 0)
    end function present

end module instant_slice