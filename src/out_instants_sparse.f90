module instant_slice_sparse
  use modglobal,  only : cexpnr, rk3step, ltempeq, lmoist, nsv, tsample, dt, timee, runtime, &
                         slicevars, lislicedump, islice, ljslicedump, jslice, lkslicedump, kslice, sparse_factor,&
                         nkslice, nislice, njslice, itot,jtot, &
                         ib, ie, jb, je, kb, ke, dzfi, dzh, &
                         timee, tstatstart, jtot, itot, kmax
  use modfields,  only : um, vm, wm, thlm, qtm, svm, pres0
  use modmpi,     only : myid, cmyidx, cmyidy
  use decomp_2d,  only : zstart, zend, xstart, xend, ystart, yend
  use modstat_nc, only : ncinfo, open_nc, define_nc, writestat_dims_nc, writestat_nc, writeoffset, writeoffset_1dx
  implicit none
  private
  public :: instant_init, instant_main
  save

  !> Sparse sampling factor: write every sparse_factor grid points
  !> Set to 1 to match original output exactly
  !> Now 3D array: (/x_factor, y_factor, z_factor/)
  integer :: xdim, ydim, zdim, kdim, idim, jdim
  real    :: tsampleslice

  logical, allocatable :: islicerank(:)   
  integer, allocatable :: isliceloc(:)    
  logical, allocatable :: jslicerank(:)   
  integer, allocatable :: jsliceloc(:)    
  
  integer :: local_nislice, local_njslice
  integer, allocatable :: local_islice_map(:)
  
  integer                    :: nslicevars
  character(80)              :: islicename
  character(80)              :: jslicename
  character(80)              :: kslicename
  character(80)              :: kslicename1d
  
  character(80)              :: slicetimeVar(1,4)
  character(80), allocatable :: isliceVars(:,:)
  character(80), allocatable :: jsliceVars(:,:)
  character(80), allocatable :: ksliceVars(:,:)
  integer                    :: ncidislice, nrecislice, ncidjslice, nrecjslice, ncidkslice, nreckslice
  integer                    :: ncidkslice1d, nreckslice1d

  contains

    subroutine instant_init
      implicit none

      if(.not.(lislicedump .or. ljslicedump .or. lkslicedump)) return

      nslicevars = (LEN(trim(slicevars))+1)/3

      if (nslicevars == 0) then
        lislicedump = .false.
        ljslicedump = .false.
        lkslicedump = .false.
        print *, "NOTE: invalid 'slicevars'..."
        return
      end if

      if(runtime <= tstatstart) then
        if(myid==0) then
          write(*,*) "ERROR: runtime <= tstatstart..."
          stop 1
        end if
      end if
      
      ! Correct dimension calculation for strided array: (N-1)/S + 1
      ! This ensures correct size even if not perfectly divisible
      xdim = (ie-ib)/sparse_factor(1) + 1
      ydim = (je-jb)/sparse_factor(2) + 1
      zdim = (ke-kb)/sparse_factor(3) + 1
      kdim = nkslice
      idim = nislice
      jdim = njslice
      tsampleslice = 0.
      
      if(myid==0) then
        write(*,*) "=== Sparse Output Init ==="
        write(*,*) "Factors (x,y,z):", sparse_factor
        write(*,*) "X-dim:", ie-ib+1, "->", xdim
        write(*,*) "Y-dim:", je-jb+1, "->", ydim
        write(*,*) "Z-dim:", ke-kb+1, "->", zdim
        write(*,*) "=========================="
      end if
      
      call ncinfo(slicetimeVar( 1,:), 'time', 'Time', 's', 'time')

      if (lislicedump) then
        allocate(isliceVars(nslicevars,4))
        allocate(islicerank(nislice))
        allocate(isliceloc(nislice))
        allocate(local_islice_map(nislice))
        call instant_ncdescription_islice
        call instant_create_ncislice
        deallocate(isliceVars)
      end if
      
      if (ljslicedump) then
        allocate(jsliceVars(nslicevars,4))
        allocate(jslicerank(njslice))
        allocate(jsliceloc(njslice))
        call instant_ncdescription_jslice
        call instant_create_ncjslice
        deallocate(jsliceVars)
      end if

      if (lkslicedump) then
        allocate(ksliceVars(nslicevars,4))
        call instant_ncdescription_kslice
        call instant_create_nckslice
        deallocate(ksliceVars)
      end if
    end subroutine instant_init


    subroutine instant_main
      implicit none
      if (timee < tstatstart) return
      if (.not. rk3step==3) return
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
            call ncinfo(isliceVars(n,:), 'u', 'Streamwise velocity', 'm/s', 'tmtt')
          case('v0')
            call ncinfo(isliceVars(n,:), 'v', 'Spanwise velocity', 'm/s', 'tttt')
          case('w0')
            call ncinfo(isliceVars(n,:), 'w', 'Vertical velocity', 'm/s', 'ttmt')
          case('p0')
            call ncinfo(isliceVars(n,:), 'pres', 'pressure', 'M', 'tttt')
          case('th')
            if (ltempeq) call ncinfo(isliceVars(n,:), 'thl', 'Potential temperature', 'K', 'tttt')
          case('qt')
            if (lmoist) call ncinfo(isliceVars(n,:), 'qt', 'Specific humidity', 'kg/kg', 'tttt')
          case('s1')
            if (nsv>0) call ncinfo(isliceVars(n,:), 's1', 'Concentration field 1', 'g/m^3', 'tttt')
          case('s2')
            if (nsv>1) call ncinfo(isliceVars(n,:), 's2', 'Concentration field 2', 'g/m^3', 'tttt')
          case('s3')
            if (nsv>2) call ncinfo(isliceVars(n,:), 's3', 'Concentration field 3', 'g/m^3', 'tttt')
          case('s4')
            if (nsv>3) call ncinfo(isliceVars(n,:), 's4', 'Concentration field 4', 'g/m^3', 'tttt')
          case default
            if(myid==0) print *, "Invalid slice variable: ", slicevars(3*n-2:3*n-1)
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
            call ncinfo(jsliceVars(n,:), 'u', 'Streamwise velocity', 'm/s', 'mttt')
          case('v0')
            call ncinfo(jsliceVars(n,:), 'v', 'Spanwise velocity', 'm/s', 'tmtt')
          case('w0')
            call ncinfo(jsliceVars(n,:), 'w', 'Vertical velocity', 'm/s', 'ttmt')
          case('p0')
            call ncinfo(jsliceVars(n,:), 'pres', 'pressure', 'M', 'tttt')
          case('th')
            if (ltempeq) call ncinfo(jsliceVars(n,:), 'thl', 'Potential temperature', 'K', 'tttt')
          case('qt')
            if (lmoist) call ncinfo(jsliceVars(n,:), 'qt', 'Specific humidity', 'kg/kg', 'tttt')
          case('s1')
            if (nsv>0) call ncinfo(jsliceVars(n,:), 's1', 'Concentration field 1', 'g/m^3', 'tttt')
          case('s2')
            if (nsv>1) call ncinfo(jsliceVars(n,:), 's2', 'Concentration field 2', 'g/m^3', 'tttt')
          case('s3')
            if (nsv>2) call ncinfo(jsliceVars(n,:), 's3', 'Concentration field 3', 'g/m^3', 'tttt')
          case('s4')
            if (nsv>3) call ncinfo(jsliceVars(n,:), 's4', 'Concentration field 4', 'g/m^3', 'tttt')
          case default
            if(myid==0) print *, "Invalid slice variable: ", slicevars(3*n-2:3*n-1)
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
            call ncinfo(ksliceVars(n,:), 'u', 'Streamwise velocity', 'm/s', 'mttt')
          case('v0')
            call ncinfo(ksliceVars(n,:), 'v', 'Spanwise velocity', 'm/s', 'tmtt')
          case('w0')
            call ncinfo(ksliceVars(n,:), 'w', 'Vertical velocity', 'm/s', 'ttmt')
          case('p0')
            call ncinfo(ksliceVars(n,:), 'pres', 'pressure', 'M', 'tttt')
          case('th')
            if (ltempeq) call ncinfo(ksliceVars(n,:), 'thl', 'Potential temperature', 'K', 'tttt')
          case('qt')
            if (lmoist) call ncinfo(ksliceVars(n,:), 'qt', 'Specific humidity', 'kg/kg', 'tttt')
          case('s1')
            if (nsv>0) call ncinfo(ksliceVars(n,:), 's1', 'Concentration field 1', 'g/m^3', 'tttt')
          case('s2')
            if (nsv>1) call ncinfo(ksliceVars(n,:), 's2', 'Concentration field 2', 'g/m^3', 'tttt')
          case('s3')
            if (nsv>2) call ncinfo(ksliceVars(n,:), 's3', 'Concentration field 3', 'g/m^3', 'tttt')
          case('s4')
            if (nsv>3) call ncinfo(ksliceVars(n,:), 's4', 'Concentration field 4', 'g/m^3', 'tttt')
          case default
            if(myid==0) print *, "Invalid slice variable: ", slicevars(3*n-2:3*n-1)
            STOP 1
        end select
      end do
    end subroutine instant_ncdescription_kslice


    subroutine instant_create_ncislice
      use modglobal, only : jtot
      use modmpi, only : myidy, myidx, nprocx
      implicit none
      integer :: i
      logical :: has_islice
      integer :: jtot_sparse
      
      if (nislice == 0) return
      
      local_nislice = 0
      has_islice = .false.
      do i = 1, nislice
        if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          local_nislice = local_nislice + 1
          has_islice = .true.
        end if
      end do
      
      if (has_islice .and. myidy == 0) then
        islicename = 'islice.xxx.xxx.nc'
        islicename(8:10) = cmyidx
        islicename(12:14) = cexpnr

        nrecislice = 0
        ! FIX: Use correct integer arithmetic for strided size: (N-1)/S + 1
        jtot_sparse = (jtot - 1) / sparse_factor(2) + 1
        
        call open_nc(islicename, ncidislice, nrecislice, n1=local_nislice, n2=jtot_sparse, n3=zdim)
        if (nrecislice == 0) then
          call define_nc(ncidislice, 1, slicetimeVar)
          call writestat_dims_nc(ncidislice)
          call write_sparse_ycoord(ncidislice)
          call write_sparse_zcoord(ncidislice)
          call write_islice_xcoord_local
        end if
        call define_nc(ncidislice, nslicevars, isliceVars)
        
        write(*,'(A,A,A,I2,A)') '  Processor (myidx=', cmyidx, ') created islice file'
      end if
    end subroutine instant_create_ncislice

    subroutine instant_create_ncjslice
      use modglobal, only : itot
      use modmpi, only : myidx, myidy, nprocy
      implicit none
      integer :: j
      logical :: has_jslice
      integer :: itot_sparse
      
      if (njslice == 0) return

      local_njslice = 0
      has_jslice = .false.
      do j = 1, njslice
        if ( (jslice(j)-1)/(jtot/nprocy) == myidy) then
          local_njslice = local_njslice + 1
          has_jslice = .true.
        end if
      end do
      
      if (has_jslice .and. myidx == 0) then
        jslicename = 'jslice.xxx.xxx.nc'
        jslicename(8:10) = cmyidy
        jslicename(12:14) = cexpnr
        
        nrecjslice = 0
        ! FIX: Use correct integer arithmetic
        itot_sparse = (itot - 1) / sparse_factor(1) + 1
        
        call open_nc(jslicename, ncidjslice, nrecjslice, n1=itot_sparse, n2=local_njslice, n3=zdim)
        if (nrecjslice==0) then
          call define_nc(ncidjslice, 1, slicetimeVar)
          call writestat_dims_nc(ncidjslice)
          call write_sparse_xcoord(ncidjslice)
          call write_sparse_zcoord(ncidjslice)
          call write_jslice_ycoord_local
        end if
        call define_nc(ncidjslice, nslicevars, jsliceVars)
        write(*,'(A,A,A,I2,A)') '  Processor (myidy=', cmyidy, ') created jslice file'
      end if
    end subroutine instant_create_ncjslice

    subroutine instant_create_nckslice
      use modglobal, only : jtot
      use modmpi, only : myidy
      implicit none
      integer :: jtot_sparse
      
      if (nkslice == 0) return

      if (myidy == 0) then
        kslicename1d = 'kslice.xxx.xxx.nc'
        kslicename1d(8:10) = cmyidx
        kslicename1d(12:14) = cexpnr

        nreckslice1d = 0
        ! FIX: Use correct integer arithmetic
        jtot_sparse = (jtot - 1) / sparse_factor(2) + 1
        
        call open_nc(kslicename1d, ncidkslice1d, nreckslice1d, n1=xdim, n2=jtot_sparse, n3=kdim)
        if (nreckslice1d == 0) then
          call define_nc(ncidkslice1d, 1, slicetimeVar)
          call writestat_dims_nc(ncidkslice1d)
          call write_sparse_ycoord(ncidkslice1d)
          call write_kslice_zcoord(ncidkslice1d)
        end if
        call define_nc(ncidkslice1d, nslicevars, ksliceVars)
      end if
    end subroutine instant_create_nckslice

    subroutine instant_write_islice
      use modmpi, only : myidy, myidx, nprocx
      use modglobal, only : jtot
      implicit none
      real, allocatable :: tmp_slice(:,:,:)
      integer :: i, ii, ii_local, local_idx
      
      if (nislice == 0) return
      if (local_nislice == 0) return
      
      allocate(tmp_slice(local_nislice, ydim, zdim))
      tmp_slice = 0.0
      
      if (myidy == 0) then
        call writestat_nc(ncidislice, 'time', timee, nrecislice, .true.)
      end if
      
      if (present('u0')) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - (myidx * (itot/nprocx))
              tmp_slice(local_idx,:,:) = 0.5*(um(ii_local,jb:je:sparse_factor(2),kb:ke:sparse_factor(3)) + um(ii_local+1,jb:je:sparse_factor(2),kb:ke:sparse_factor(3)))
          end if
        end do
        call writeoffset(ncidislice, 'u', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if

      if (present('v0')) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - (myidx * (itot/nprocx))
              tmp_slice(local_idx,:,:) = vm(ii_local,jb:je:sparse_factor(2),kb:ke:sparse_factor(3))
          end if
        end do
        call writeoffset(ncidislice, 'v', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      if (present('w0')) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - (myidx * (itot/nprocx))           
              tmp_slice(local_idx,:,:) = wm(ii_local,jb:je:sparse_factor(2),kb:ke:sparse_factor(3))
          end if
        end do
        call writeoffset(ncidislice, 'w', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      if (present('p0')) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - (myidx * (itot/nprocx))           
              tmp_slice(local_idx,:,:) = pres0(ii_local,jb:je:sparse_factor(2),kb:ke:sparse_factor(3))
          end if
        end do
        call writeoffset(ncidislice, 'pres', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if

      if (present('th') .and. ltempeq) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - (myidx * (itot/nprocx))
              tmp_slice(local_idx,:,:) = thlm(ii_local,jb:je:sparse_factor(2),kb:ke:sparse_factor(3))
          end if
        end do
        call writeoffset(ncidislice, 'thl', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      if (present('qt') .and. lmoist) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - (myidx * (itot/nprocx))
              tmp_slice(local_idx,:,:) = qtm(ii_local,jb:je:sparse_factor(2),kb:ke:sparse_factor(3))
          end if
        end do
        call writeoffset(ncidislice, 'qt', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      if (present('s1') .and. nsv>0) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - (myidx * (itot/nprocx))
              tmp_slice(local_idx,:,:) = svm(ii_local,jb:je:sparse_factor(2),kb:ke:sparse_factor(3),1)
          end if
        end do
        call writeoffset(ncidislice, 's1', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      if (present('s2') .and. nsv>1) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - (myidx * (itot/nprocx))
              tmp_slice(local_idx,:,:) = svm(ii_local,jb:je:sparse_factor(2),kb:ke:sparse_factor(3),2)
          end if
        end do
        call writeoffset(ncidislice, 's2', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      if (present('s3') .and. nsv>2) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - (myidx * (itot/nprocx))
              tmp_slice(local_idx,:,:) = svm(ii_local,jb:je:sparse_factor(2),kb:ke:sparse_factor(3),3)
          end if
        end do
        call writeoffset(ncidislice, 's3', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      if (present('s4') .and. nsv>3) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
            local_idx = local_idx + 1
            ii = islice(i)
            ii_local = ii - (myidx * (itot/nprocx))
              tmp_slice(local_idx,:,:) = svm(ii_local,jb:je:sparse_factor(2),kb:ke:sparse_factor(3),4)
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

      allocate(tmp_slice(xdim, local_njslice, zdim))
      tmp_slice = 0.0

      if (myidx == 0) then
        call writestat_nc(ncidjslice, 'time', timee, nrecjslice, .true.)
      end if

      if (present('u0')) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/(jtot/nprocy) == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - (myidy * (jtot/nprocy))
              tmp_slice(:, local_idy, :) = um(ib:ie:sparse_factor(1), jj_local, kb:ke:sparse_factor(3))
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'u', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      if (present('v0')) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/(jtot/nprocy) == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - (myidy * (jtot/nprocy))
              tmp_slice(:, local_idy, :) = 0.5*(vm(ib:ie:sparse_factor(1), jj_local, kb:ke:sparse_factor(3)) + vm(ib:ie:sparse_factor(1), jj_local+1, kb:ke:sparse_factor(3)))
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'v', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      if (present('w0')) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/(jtot/nprocy) == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - (myidy * (jtot/nprocy))
              tmp_slice(:, local_idy, :) = wm(ib:ie:sparse_factor(1), jj_local, kb:ke:sparse_factor(3))
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'w', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      if (present('p0')) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/(jtot/nprocy) == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - (myidy * (jtot/nprocy))
              tmp_slice(:, local_idy, :) = pres0(ib:ie:sparse_factor(1), jj_local, kb:ke:sparse_factor(3))
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'pres', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      if (present('th') .and. ltempeq) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/(jtot/nprocy) == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - (myidy * (jtot/nprocy))
              tmp_slice(:, local_idy, :) = thlm(ib:ie:sparse_factor(1), jj_local, kb:ke:sparse_factor(3))
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'thl', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      if (present('qt') .and. lmoist) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/(jtot/nprocy) == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - (myidy * (jtot/nprocy))
              tmp_slice(:, local_idy, :) = qtm(ib:ie:sparse_factor(1), jj_local, kb:ke:sparse_factor(3))
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'qt', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      if (present('s1') .and. nsv>0) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/(jtot/nprocy) == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - (myidy * (jtot/nprocy))
              tmp_slice(:, local_idy, :) = svm(ib:ie:sparse_factor(1), jj_local, kb:ke:sparse_factor(3), 1)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 's1', tmp_slice, nrecjslice, xdim,  local_njslice, zdim)
      end if

      if (present('s2') .and. nsv>1) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/(jtot/nprocy) == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - (myidy * (jtot/nprocy))
              tmp_slice(:, local_idy, :) = svm(ib:ie:sparse_factor(1), jj_local, kb:ke:sparse_factor(3), 2)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 's2', tmp_slice, nrecjslice, xdim,  local_njslice, zdim)
      end if

      if (present('s3') .and. nsv>2) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/(jtot/nprocy) == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - (myidy * (jtot/nprocy))
              tmp_slice(:, local_idy, :) = svm(ib:ie:sparse_factor(1), jj_local, kb:ke:sparse_factor(3), 3)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 's3', tmp_slice, nrecjslice, xdim,  local_njslice, zdim)
      end if

      if (present('s4') .and. nsv>3) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
          if ( (jslice(j)-1)/(jtot/nprocy) == myidy) then
            local_idy = local_idy + 1
            jj = jslice(j)
            jj_local = jj - (myidy * (jtot/nprocy))
              tmp_slice(:, local_idy, :) = svm(ib:ie:sparse_factor(1), jj_local, kb:ke:sparse_factor(3), 4)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 's4', tmp_slice, nrecjslice, xdim,  local_njslice, zdim)
      end if

      deallocate(tmp_slice)
    end subroutine instant_write_jslice

    subroutine write_islice_xcoord_local
      use modglobal, only : xf, xh
      use modmpi, only : myidx, nprocx
      use netcdf
      implicit none
      integer :: varid, ierr, i, local_idx
      real, allocatable :: x_islice_f(:), x_islice_h(:)
      integer, allocatable :: islice_indices(:)
      
      allocate(x_islice_f(local_nislice))
      allocate(x_islice_h(local_nislice))
      allocate(islice_indices(local_nislice))
      
      local_idx = 0
      do i = 1, nislice
        if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          local_idx = local_idx + 1
          x_islice_f(local_idx) = xf(islice(i))
          x_islice_h(local_idx) = xh(islice(i))
          islice_indices(local_idx) = islice(i)
        end if
      end do
      write(*,*) "local_idx=", local_idx, "myidx=", myidx
      ierr = nf90_inq_varid(ncidislice, 'xt', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncidislice)
        ierr = nf90_put_att(ncidislice, varid, 'islice_indices', islice_indices)
        ierr = nf90_put_att(ncidislice, varid, 'long_name', 'x-coordinate of i-slices (cell center)')
        ierr = nf90_enddef(ncidislice)
        ierr = nf90_put_var(ncidislice, varid, x_islice_f)
      end if
      
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
        if ( (jslice(j)-1)/(jtot/nprocy) == myidy) then
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
      use modglobal, only : zf, zh
      use netcdf
      implicit none
      integer, intent(in) :: ncid
      integer :: varid, ierr, k
      real, allocatable :: z_kslice_f(:), z_kslice_h(:)
      
      allocate(z_kslice_f(nkslice))
      allocate(z_kslice_h(nkslice))
      do k = 1, nkslice
        z_kslice_f(k) = zf(kslice(k))
        z_kslice_h(k) = zh(kslice(k))
      end do
      
      ierr = nf90_inq_varid(ncid, 'zt', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'kslice_indices', kslice(1:nkslice))
        ierr = nf90_put_att(ncid, varid, 'long_name', 'z-coordinate of k-slices (cell center)')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, z_kslice_f)
      end if
      
      ierr = nf90_inq_varid(ncid, 'zm', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'kslice_indices', kslice(1:nkslice))
        ierr = nf90_put_att(ncid, varid, 'long_name', 'z-coordinate of k-slices (cell edge)')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, z_kslice_h)
      end if
      
      deallocate(z_kslice_f)
      deallocate(z_kslice_h)
    end subroutine write_kslice_zcoord

    subroutine write_sparse_xcoord(ncid)
      use modglobal, only : xf, xh, itot
      use netcdf
      implicit none
      integer, intent(in) :: ncid
      integer :: varid, ierr
      real, allocatable :: x_sparse_f(:), x_sparse_h(:)
      integer :: nx_sparse
      
      nx_sparse = size(xf(1:itot:sparse_factor(1)))
      
      allocate(x_sparse_f(nx_sparse))
      allocate(x_sparse_h(nx_sparse))
      
      x_sparse_f = xf(1:itot:sparse_factor(1))
      x_sparse_h = xh(1:itot:sparse_factor(1))
      
      ierr = nf90_inq_varid(ncid, 'xt', varid)
      if (ierr == nf90_noerr) then
         ierr = nf90_put_var(ncid, varid, x_sparse_f)
      end if
      
      ierr = nf90_inq_varid(ncid, 'xm', varid)
      if (ierr == nf90_noerr) then
         ierr = nf90_put_var(ncid, varid, x_sparse_h)
      end if
      
      deallocate(x_sparse_f)
      deallocate(x_sparse_h)
    end subroutine write_sparse_xcoord

    subroutine write_sparse_ycoord(ncid)
      use modglobal, only : yf, yh, jtot
      use netcdf
      implicit none
      integer, intent(in) :: ncid
      integer :: varid, ierr
      real, allocatable :: y_sparse_f(:), y_sparse_h(:)
      integer :: ny_sparse
      
      ny_sparse = size(yf(1:jtot:sparse_factor(2)))
      
      allocate(y_sparse_f(ny_sparse))
      allocate(y_sparse_h(ny_sparse))
      
      y_sparse_f = yf(1:jtot:sparse_factor(2))
      y_sparse_h = yh(1:jtot:sparse_factor(2))
      
      ierr = nf90_inq_varid(ncid, 'yt', varid)
      if (ierr == nf90_noerr) then
         ierr = nf90_put_var(ncid, varid, y_sparse_f)
      end if
      
      ierr = nf90_inq_varid(ncid, 'ym', varid)
      if (ierr == nf90_noerr) then
         ierr = nf90_put_var(ncid, varid, y_sparse_h)
      end if
      
      deallocate(y_sparse_f)
      deallocate(y_sparse_h)
    end subroutine write_sparse_ycoord

    subroutine write_sparse_zcoord(ncid)
      use modglobal, only : zf, zh, kmax
      use netcdf
      implicit none
      integer, intent(in) :: ncid
      integer :: varid, ierr
      real, allocatable :: z_sparse_f(:), z_sparse_h(:)
      integer :: nz_sparse
      
      nz_sparse = size(zf(1:kmax:sparse_factor(3)))
      
      allocate(z_sparse_f(nz_sparse))
      allocate(z_sparse_h(nz_sparse))
      
      z_sparse_f = zf(1:kmax:sparse_factor(3))
      z_sparse_h = zh(1:kmax:sparse_factor(3))
      
      ierr = nf90_inq_varid(ncid, 'zt', varid)
      if (ierr == nf90_noerr) then
         ierr = nf90_put_var(ncid, varid, z_sparse_f)
      end if
      
      ierr = nf90_inq_varid(ncid, 'zm', varid)
      if (ierr == nf90_noerr) then
         ierr = nf90_put_var(ncid, varid, z_sparse_h)
      end if
      
      deallocate(z_sparse_f)
      deallocate(z_sparse_h)
    end subroutine write_sparse_zcoord

    subroutine instant_write_kslice
      use modmpi, only : myidy
      use modglobal, only : jtot
      implicit none
      real, allocatable :: tmp_slice(:,:,:)
      integer :: k, kk
      
      if (nkslice == 0) return
      
      allocate(tmp_slice(xdim, ydim, nkslice))
      
      if (myidy == 0) then
        call writestat_nc(ncidkslice1d, 'time', timee, nreckslice1d, .true.)
      end if
      
      if (present('u0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = um(ib:ie:sparse_factor(1), jb:je:sparse_factor(2), kk)
        end do
        call writeoffset(ncidkslice1d, 'u', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if
      
      if (present('v0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = vm(ib:ie:sparse_factor(1), jb:je:sparse_factor(2), kk)
        end do
        call writeoffset(ncidkslice1d, 'v', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if
      
      if (present('w0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = 0.5 * (wm(ib:ie:sparse_factor(1), jb:je:sparse_factor(2), kk) + wm(ib:ie:sparse_factor(1), jb:je:sparse_factor(2), kk+1))
        end do
        call writeoffset(ncidkslice1d, 'w', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if

      if (present('p0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = pres0(ib:ie:sparse_factor(1), jb:je:sparse_factor(2), kk)
        end do
        call writeoffset(ncidkslice1d, 'pres', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if

      if (present('th') .and. ltempeq) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = thlm(ib:ie:sparse_factor(1), jb:je:sparse_factor(2), kk)
        end do
        call writeoffset(ncidkslice1d, 'thl', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if
      
      if (present('qt') .and. lmoist) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = qtm(ib:ie:sparse_factor(1), jb:je:sparse_factor(2), kk)
        end do
        call writeoffset(ncidkslice1d, 'qt', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if
      
      if (present('s1') .and. nsv>0) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie:sparse_factor(1), jb:je:sparse_factor(2), kk, 1)
        end do
        call writeoffset(ncidkslice1d, 's1', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if
      
      if (present('s2') .and. nsv>1) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie:sparse_factor(1), jb:je:sparse_factor(2), kk, 2)
        end do
        call writeoffset(ncidkslice1d, 's2', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if
      
      if (present('s3') .and. nsv>2) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie:sparse_factor(1), jb:je:sparse_factor(2), kk, 3)
        end do
        call writeoffset(ncidkslice1d, 's3', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if
      
      if (present('s4') .and. nsv>3) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie:sparse_factor(1), jb:je:sparse_factor(2), kk, 4)
        end do
        call writeoffset(ncidkslice1d, 's4', tmp_slice, nreckslice1d, xdim, ydim, kdim)
      end if
      
      deallocate(tmp_slice)
    end subroutine instant_write_kslice

    logical function present(str)
      character(len=*), intent(in) :: str
      integer :: pos

      pos = index(slicevars, str)
      present = (pos > 0)
    end function present

end module instant_slice_sparse
