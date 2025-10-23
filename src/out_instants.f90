module instant_slice
  use modglobal,  only : cexpnr, rk3step, ltempeq, lmoist, nsv, tsample, dt, timee, runtime, &
                         slicevars, lislicedump, islice, ljslicedump, jslice, lkslicedump, kslice, &
                         nkslice, nislice, njslice, &
                         ib, ie, jb, je, kb, ke, dzfi, dzh, &
                         timee, tstatstart
  use modfields,  only : um, vm, wm, thlm, qtm, svm
  use modmpi,     only : myid, cmyidx, cmyidy
  use decomp_2d,  only : zstart, zend, xstart, xend, ystart, yend
  use modstat_nc, only : ncinfo, open_nc, define_nc, writestat_dims_nc, writestat_nc
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
  
  integer                    :: nslicevars
  character(80)              :: islicename
  character(80)              :: jslicename
  character(80)              :: kslicename
  character(80)              :: slicetimeVar(1,4)
  character(80), allocatable :: isliceVars(:,:)
  character(80), allocatable :: jsliceVars(:,:)
  character(80), allocatable :: ksliceVars(:,:)
  integer                    :: ncidislice, nrecislice, ncidjslice, nrecjslice, ncidkslice, nreckslice

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
        call instant_ncdescription_islice
        call instant_create_ncislice    !> Generate sliced NetCDF: islice.xxx.xxx.xxx.nc
        deallocate(isliceVars)
      end if
      
      if (ljslicedump) then
        allocate(jsliceVars(nslicevars,4))
        allocate(jslicerank(njslice))
        allocate(jsliceloc(njslice))
        call instant_ncdescription_jslice
        call instant_create_ncjslice    !> Generate sliced NetCDF: jslice.xxx.xxx.xxx.nc
        deallocate(jsliceVars)
      end if

      if (lkslicedump) then
        allocate(ksliceVars(nslicevars,4))
        call instant_ncdescription_kslice
        call instant_create_nckslice    !> Generate sliced NetCDF: kslice.xxx.xxx.xxx.nc
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


    subroutine instant_create_ncislice
      implicit none
      integer :: i
      
      if (nislice == 0) then
        if (myid == 0) write(*,*) "WARNING: nislice=0, no i-slices will be output"
        return
      end if
      
      ! Check which islices are on this processor and store their local positions
      do i = 1, nislice
        if ((islice(i) >= xstart(1)) .and. (islice(i) <= xend(1))) then
          islicerank(i) = .true.
          isliceloc(i) = islice(i) - xstart(1) + 1
        else
          islicerank(i) = .false.
          isliceloc(i) = -1
        end if
      end do
      
      ! All processors create one file: islice.xxx.xxx.xxx.nc
      ! Use idim as the first dimension (number of islices)
      islicename = 'islice.' // cmyidx // '.' // cmyidy // '.' // cexpnr // '.nc'
      call open_nc(islicename, ncidislice, nrecislice, n1=idim, n2=ydim, n3=zdim)
      if (nrecislice==0) then
        call define_nc(ncidislice, 1, slicetimeVar)
        call writestat_dims_nc(ncidislice)
      end if
      call define_nc(ncidislice, nslicevars, isliceVars)
      
      ! Write x-coordinates for the islice positions
      call write_islice_xcoord
    end subroutine instant_create_ncislice

    subroutine instant_create_ncjslice
      implicit none
      integer :: j
      
      if (njslice == 0) then
        if (myid == 0) write(*,*) "WARNING: njslice=0, no j-slices will be output"
        return
      end if
      
      ! Check which jslices are on this processor and store their local positions
      do j = 1, njslice
        if ((jslice(j) >= ystart(2)) .and. (jslice(j) <= yend(2))) then
          jslicerank(j) = .true.
          jsliceloc(j) = jslice(j) - ystart(2) + 1
        else
          jslicerank(j) = .false.
          jsliceloc(j) = -1
        end if
      end do
      
      ! All processors create one file: jslice.xxx.xxx.xxx.nc
      ! Use jdim as the second dimension (number of jslices)
      jslicename = 'jslice.' // cmyidx // '.' // cmyidy // '.' // cexpnr // '.nc'
      call open_nc(jslicename, ncidjslice, nrecjslice, n1=xdim, n2=jdim, n3=zdim)
      if (nrecjslice==0) then
        call define_nc(ncidjslice, 1, slicetimeVar)
        call writestat_dims_nc(ncidjslice)
      end if
      call define_nc(ncidjslice, nslicevars, jsliceVars)
      
      ! Write y-coordinates for the jslice positions
      call write_jslice_ycoord
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


    subroutine instant_write_islice
      implicit none
      real, allocatable :: tmp_slice(:,:,:)
      integer :: i, ii
      
      if (nislice == 0) return
      
      ! Allocate temporary array to hold all islice levels: (nislice, y, z)
      allocate(tmp_slice(nislice, jb:je, kb:ke))
      
      ! Write time variable
      call writestat_nc(ncidislice, 'time', timee, nrecislice, .true.)
      
      ! u velocity (interpolated to cell centers in x-direction)
      if (present('u0')) then
        tmp_slice = 0.0  ! Initialize with zeros for processors without the slice
        do i = 1, nislice
          if (islicerank(i)) then
            ii = isliceloc(i)
            tmp_slice(i,:,:) = 0.5*(um(ii,jb:je,kb:ke) + um(ii+1,jb:je,kb:ke))
          end if
        end do
        call writestat_nc(ncidislice, 'u', tmp_slice, nrecislice, idim, ydim, zdim)
      end if
      
      ! v velocity
      if (present('v0')) then
        tmp_slice = 0.0
        do i = 1, nislice
          if (islicerank(i)) then
            ii = isliceloc(i)
            tmp_slice(i,:,:) = vm(ii,jb:je,kb:ke)
          end if
        end do
        call writestat_nc(ncidislice, 'v', tmp_slice, nrecislice, idim, ydim, zdim)
      end if
      
      ! w velocity
      if (present('w0')) then
        tmp_slice = 0.0
        do i = 1, nislice
          if (islicerank(i)) then
            ii = isliceloc(i)
            tmp_slice(i,:,:) = wm(ii,jb:je,kb:ke)
          end if
        end do
        call writestat_nc(ncidislice, 'w', tmp_slice, nrecislice, idim, ydim, zdim)
      end if
      
      ! temperature
      if (present('th') .and. ltempeq) then
        tmp_slice = 0.0
        do i = 1, nislice
          if (islicerank(i)) then
            ii = isliceloc(i)
            tmp_slice(i,:,:) = thlm(ii,jb:je,kb:ke)
          end if
        end do
        call writestat_nc(ncidislice, 'thl', tmp_slice, nrecislice, idim, ydim, zdim)
      end if
      
      ! moisture
      if (present('qt') .and. lmoist) then
        tmp_slice = 0.0
        do i = 1, nislice
          if (islicerank(i)) then
            ii = isliceloc(i)
            tmp_slice(i,:,:) = qtm(ii,jb:je,kb:ke)
          end if
        end do
        call writestat_nc(ncidislice, 'qt', tmp_slice, nrecislice, idim, ydim, zdim)
      end if
      
      ! scalars s1-s4
      if (present('s1') .and. nsv>0) then
        tmp_slice = 0.0
        do i = 1, nislice
          if (islicerank(i)) then
            ii = isliceloc(i)
            tmp_slice(i,:,:) = svm(ii,jb:je,kb:ke,1)
          end if
        end do
        call writestat_nc(ncidislice, 's1', tmp_slice, nrecislice, idim, ydim, zdim)
      end if
      
      if (present('s2') .and. nsv>1) then
        tmp_slice = 0.0
        do i = 1, nislice
          if (islicerank(i)) then
            ii = isliceloc(i)
            tmp_slice(i,:,:) = svm(ii,jb:je,kb:ke,2)
          end if
        end do
        call writestat_nc(ncidislice, 's2', tmp_slice, nrecislice, idim, ydim, zdim)
      end if
      
      if (present('s3') .and. nsv>2) then
        tmp_slice = 0.0
        do i = 1, nislice
          if (islicerank(i)) then
            ii = isliceloc(i)
            tmp_slice(i,:,:) = svm(ii,jb:je,kb:ke,3)
          end if
        end do
        call writestat_nc(ncidislice, 's3', tmp_slice, nrecislice, idim, ydim, zdim)
      end if
      
      if (present('s4') .and. nsv>3) then
        tmp_slice = 0.0
        do i = 1, nislice
          if (islicerank(i)) then
            ii = isliceloc(i)
            tmp_slice(i,:,:) = svm(ii,jb:je,kb:ke,4)
          end if
        end do
        call writestat_nc(ncidislice, 's4', tmp_slice, nrecislice, idim, ydim, zdim)
      end if
      
      deallocate(tmp_slice)
    end subroutine instant_write_islice

    subroutine instant_write_jslice
      implicit none
      real, allocatable :: tmp_slice(:,:,:)
      integer :: j, jj
      
      if (njslice == 0) return
      
      ! Allocate temporary array to hold all jslice levels: (x, njslice, z)
      allocate(tmp_slice(ib:ie, njslice, kb:ke))
      
      ! Write time variable
      call writestat_nc(ncidjslice, 'time', timee, nrecjslice, .true.)
      
      ! u velocity
      if (present('u0')) then
        tmp_slice = 0.0  ! Initialize with zeros for processors without the slice
        do j = 1, njslice
          if (jslicerank(j)) then
            jj = jsliceloc(j)
            tmp_slice(:,j,:) = um(ib:ie,jj,kb:ke)
          end if
        end do
        call writestat_nc(ncidjslice, 'u', tmp_slice, nrecjslice, xdim, jdim, zdim)
      end if
      
      ! v velocity (interpolated to cell centers in y-direction)
      if (present('v0')) then
        tmp_slice = 0.0
        do j = 1, njslice
          if (jslicerank(j)) then
            jj = jsliceloc(j)
            tmp_slice(:,j,:) = 0.5*(vm(ib:ie,jj,kb:ke) + vm(ib:ie,jj+1,kb:ke))
          end if
        end do
        call writestat_nc(ncidjslice, 'v', tmp_slice, nrecjslice, xdim, jdim, zdim)
      end if
      
      ! w velocity
      if (present('w0')) then
        tmp_slice = 0.0
        do j = 1, njslice
          if (jslicerank(j)) then
            jj = jsliceloc(j)
            tmp_slice(:,j,:) = wm(ib:ie,jj,kb:ke)
          end if
        end do
        call writestat_nc(ncidjslice, 'w', tmp_slice, nrecjslice, xdim, jdim, zdim)
      end if
      
      ! temperature
      if (present('th') .and. ltempeq) then
        tmp_slice = 0.0
        do j = 1, njslice
          if (jslicerank(j)) then
            jj = jsliceloc(j)
            tmp_slice(:,j,:) = thlm(ib:ie,jj,kb:ke)
          end if
        end do
        call writestat_nc(ncidjslice, 'thl', tmp_slice, nrecjslice, xdim, jdim, zdim)
      end if
      
      ! moisture
      if (present('qt') .and. lmoist) then
        tmp_slice = 0.0
        do j = 1, njslice
          if (jslicerank(j)) then
            jj = jsliceloc(j)
            tmp_slice(:,j,:) = qtm(ib:ie,jj,kb:ke)
          end if
        end do
        call writestat_nc(ncidjslice, 'qt', tmp_slice, nrecjslice, xdim, jdim, zdim)
      end if
      
      ! scalars s1-s4
      if (present('s1') .and. nsv>0) then
        tmp_slice = 0.0
        do j = 1, njslice
          if (jslicerank(j)) then
            jj = jsliceloc(j)
            tmp_slice(:,j,:) = svm(ib:ie,jj,kb:ke,1)
          end if
        end do
        call writestat_nc(ncidjslice, 's1', tmp_slice, nrecjslice, xdim, jdim, zdim)
      end if
      
      if (present('s2') .and. nsv>1) then
        tmp_slice = 0.0
        do j = 1, njslice
          if (jslicerank(j)) then
            jj = jsliceloc(j)
            tmp_slice(:,j,:) = svm(ib:ie,jj,kb:ke,2)
          end if
        end do
        call writestat_nc(ncidjslice, 's2', tmp_slice, nrecjslice, xdim, jdim, zdim)
      end if
      
      if (present('s3') .and. nsv>2) then
        tmp_slice = 0.0
        do j = 1, njslice
          if (jslicerank(j)) then
            jj = jsliceloc(j)
            tmp_slice(:,j,:) = svm(ib:ie,jj,kb:ke,3)
          end if
        end do
        call writestat_nc(ncidjslice, 's3', tmp_slice, nrecjslice, xdim, jdim, zdim)
      end if
      
      if (present('s4') .and. nsv>3) then
        tmp_slice = 0.0
        do j = 1, njslice
          if (jslicerank(j)) then
            jj = jsliceloc(j)
            tmp_slice(:,j,:) = svm(ib:ie,jj,kb:ke,4)
          end if
        end do
        call writestat_nc(ncidjslice, 's4', tmp_slice, nrecjslice, xdim, jdim, zdim)
      end if
      
      deallocate(tmp_slice)
    end subroutine instant_write_jslice

    subroutine write_islice_xcoord
      ! Update x-coordinate to reflect islice positions
      use modglobal, only : dx, xf
      use netcdf
      implicit none
      integer :: varid, ierr, i
      real, allocatable :: x_islice(:)
      
      ! Allocate and fill x coordinates for islices
      allocate(x_islice(nislice))
      do i = 1, nislice
        x_islice(i) = xf(islice(i))
      end do
      
      ! Find the x variable and overwrite it
      ierr = nf90_inq_varid(ncidislice, 'xt', varid)
      if (ierr /= nf90_noerr) then
        if (myid == 0) write(*,*) 'WARNING: Cannot find xt variable, trying xm'
        ierr = nf90_inq_varid(ncidislice, 'xm', varid)
      end if
      
      if (ierr == nf90_noerr) then
        ! Add attributes
        ierr = nf90_redef(ncidislice)
        ierr = nf90_put_att(ncidislice, varid, 'islice_indices', islice(1:nislice))
        ierr = nf90_put_att(ncidislice, varid, 'long_name', 'x-coordinate of i-slices')
        ierr = nf90_enddef(ncidislice)
        
        ! Write the islice x-coordinates
        ierr = nf90_put_var(ncidislice, varid, x_islice)
      end if
      
      deallocate(x_islice)
    end subroutine write_islice_xcoord

    subroutine write_jslice_ycoord
      ! Update y-coordinate to reflect jslice positions
      use modglobal, only : dy, yf
      use netcdf
      implicit none
      integer :: varid, ierr, j
      real, allocatable :: y_jslice(:)
      
      ! Allocate and fill y coordinates for jslices
      allocate(y_jslice(njslice))
      do j = 1, njslice
        y_jslice(j) = yf(jslice(j))
      end do
      
      ! Find the y variable and overwrite it
      ierr = nf90_inq_varid(ncidjslice, 'yt', varid)
      if (ierr /= nf90_noerr) then
        if (myid == 0) write(*,*) 'WARNING: Cannot find yt variable, trying ym'
        ierr = nf90_inq_varid(ncidjslice, 'ym', varid)
      end if
      
      if (ierr == nf90_noerr) then
        ! Add attributes
        ierr = nf90_redef(ncidjslice)
        ierr = nf90_put_att(ncidjslice, varid, 'jslice_indices', jslice(1:njslice))
        ierr = nf90_put_att(ncidjslice, varid, 'long_name', 'y-coordinate of j-slices')
        ierr = nf90_enddef(ncidjslice)
        
        ! Write the jslice y-coordinates
        ierr = nf90_put_var(ncidjslice, varid, y_jslice)
      end if
      
      deallocate(y_jslice)
    end subroutine write_jslice_ycoord

    subroutine write_kslice_zcoord(ncid)
      ! Update z-coordinate to reflect kslice levels instead of full z levels
      use modglobal, only : zh
      use netcdf
      implicit none
      integer, intent(in) :: ncid
      integer :: varid, ierr, k
      real, allocatable :: z_kslice(:)
      
      ! Allocate and fill z coordinates for kslices
      allocate(z_kslice(nkslice))
      do k = 1, nkslice
        z_kslice(k) = zh(kslice(k))
      end do
      
      ! Find the z variable and overwrite it with kslice heights
      ierr = nf90_inq_varid(ncid, 'zt', varid)
      if (ierr /= nf90_noerr) then
        if (myid == 0) write(*,*) 'WARNING: Cannot find zt variable, trying zm'
        ierr = nf90_inq_varid(ncid, 'zm', varid)
      end if
      
      if (ierr == nf90_noerr) then
        ! Add attribute to indicate these are kslice levels
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'kslice_indices', kslice(1:nkslice))
        ierr = nf90_put_att(ncid, varid, 'long_name', 'z-coordinate of k-slices')
        ierr = nf90_enddef(ncid)
        
        ! Write the kslice z-coordinates
        ierr = nf90_put_var(ncid, varid, z_kslice)
        
        if (ierr == nf90_noerr .and. myid == 0) then
          write(*,*) 'Updated z-coordinates for kslices:'
          write(*,*) '  z values (m):', z_kslice
          write(*,*) '  k indices:', kslice(1:nkslice)
        end if
      else
        if (myid == 0) then
          write(*,*) 'ERROR: Cannot find z coordinate variable'
        end if
      end if
      
      deallocate(z_kslice)
    end subroutine write_kslice_zcoord

    subroutine instant_write_kslice
      implicit none
      real, allocatable :: tmp_slice(:,:,:)
      integer :: k, kk
      
      if (nkslice == 0) return
      
      ! Allocate temporary array to hold all kslice levels: (x, y, nkslice)
      allocate(tmp_slice(ib:ie, jb:je, nkslice))
      
      ! Write time variable
      call writestat_nc(ncidkslice, 'time', timee, nreckslice, .true.)
      
      ! u velocity
      if (present('u0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = um(ib:ie, jb:je, kk)
        end do
        call writestat_nc(ncidkslice, 'u', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      ! v velocity
      if (present('v0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = vm(ib:ie, jb:je, kk)
        end do
        call writestat_nc(ncidkslice, 'v', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      ! w velocity (interpolated to cell centers)
      if (present('w0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = 0.5 * (wm(ib:ie, jb:je, kk) + wm(ib:ie, jb:je, kk+1))
        end do
        call writestat_nc(ncidkslice, 'w', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      ! temperature
      if (present('th') .and. ltempeq) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = thlm(ib:ie, jb:je, kk)
        end do
        call writestat_nc(ncidkslice, 'thl', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      ! moisture
      if (present('qt') .and. lmoist) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = qtm(ib:ie, jb:je, kk)
        end do
        call writestat_nc(ncidkslice, 'qt', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      ! scalars s1-s4
      if (present('s1') .and. nsv>0) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie, jb:je, kk, 1)
        end do
        call writestat_nc(ncidkslice, 's1', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      if (present('s2') .and. nsv>1) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie, jb:je, kk, 2)
        end do
        call writestat_nc(ncidkslice, 's2', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      if (present('s3') .and. nsv>2) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie, jb:je, kk, 3)
        end do
        call writestat_nc(ncidkslice, 's3', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      if (present('s4') .and. nsv>3) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie, jb:je, kk, 4)
        end do
        call writestat_nc(ncidkslice, 's4', tmp_slice, nreckslice, xdim, ydim, kdim)
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