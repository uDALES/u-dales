!> \file out_instants.f90
!! Writes instantaneous planes of various fields to different NetCDF files.
!>
!
!! Inspired from the slice routines in uDALES v2.2.0 modstatsdump.f90 written by Sam O. Owens, ICL (2024).
!! \author Jingzi Huang, ICL (2024-2026)
!! \author Dipanjan Majumdar, ICL (2023-2026)
!! \todo documentation
!
! This file is part of uDALES (https://github.com/uDALES/u-dales).
!
! uDALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! uDALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) 2016- the uDALES Team, Imperial College London.
!
module instant
  use modglobal,  only : cexpnr, rk3step, ltempeq, lmoist, nsv, tsample, dt, timee, runtime, &
                         slicevars, lislicedump, islice, ljslicedump, jslice, lkslicedump, kslice, &
                         nkslice, nislice, njslice, itot,jtot, &
                         ib, ie, jb, je, kb, ke, &
                         tstatstart, &
                         probevars, lprobedump, iprobe, jprobe, kprobe, nprobe, &
                         xf, yf, zf, xh, yh, zh
  use modfields,  only : um, vm, wm, thlm, qtm, svm, pres0
  use modmpi,     only : myid, myidy, cmyidx, cmyidy, comm3d, mpierr, my_real, NPROCY, NPROCX,myidx
  use decomp_2d,  only : zstart, zend
  use modstat_nc, only : ncinfo, open_nc, define_nc, writestat_dims_nc, writestat_nc, writeoffset, writeoffset_1dx
  use mpi
  use netcdf
  implicit none
  private
  public :: slice_init, slice_main, probe_init, probe_main
  save

  integer :: xdim, ydim, zdim, kdim  ! Added kdim for multiple slices
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
  
  character(80)              :: slicetimeVar(1,4)
  character(80), allocatable :: isliceVars(:,:)
  character(80), allocatable :: jsliceVars(:,:)
  character(80), allocatable :: ksliceVars(:,:)
  integer                    :: ncidislice, nrecislice, ncidjslice, nrecjslice, ncidkslice, nreckslice

  ! Probe variables
  real    :: tsampleprobe
  integer :: nprobevars
  integer :: ncid_probe
  integer :: point_dimid, time_dimid
  integer :: varid_xt, varid_xm, varid_yt, varid_ym, varid_zt, varid_zm, varid_time
  integer, allocatable :: varid_vals(:)
  integer :: nrecprobe = 0

  contains

    subroutine slice_init
      implicit none

      if(.not.(lislicedump .or. ljslicedump .or. lkslicedump)) return

      ! Validate slice inputs before proceeding
      call instant_validate_out_time('slices')
      call instant_validate_output_vars(nslicevars, slicevars, 'slicevars', 'one or more of lislicedump, ljslicedump and lkslicedump')  ! Validate slicevars only if at least one slice dump is true
      call instant_validate_slice_inputs(lkslicedump, nkslice, kslice, 'kslice')
      call instant_validate_slice_inputs(lislicedump, nislice, islice, 'islice') 
      call instant_validate_slice_inputs(ljslicedump, njslice, jslice, 'jslice')
      
      xdim = ie-ib+1
      ydim = je-jb+1
      zdim = ke-kb+1
      kdim = nkslice  ! Set kdim to number of kslices
      tsampleslice = 0.
      call ncinfo(slicetimeVar( 1,:), 'time', 'Time', 's', 'time')

      if (lislicedump) then
        allocate(isliceVars(nslicevars,4))
        allocate(islicerank(nislice))
        allocate(isliceloc(nislice))
        allocate(local_islice_map(nislice))  ! Allocate mapping array
        call instant_ncdescription_islice
        call instant_create_ncislice    !> Generate sliced NetCDF: ins_islice.xxx.xxx.nc
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
        call instant_create_nckslice    !> Generate sliced NetCDF: ins_kslice.xxx.xxx.nc
        deallocate(ksliceVars)
      end if
    end subroutine slice_init


    subroutine slice_main
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
    end subroutine slice_main


    subroutine instant_ncdescription_islice
      implicit none
      integer :: n
      do n=1,nslicevars
        select case(slicevars(3*n-2:3*n-1))
          case('u0')
            call ncinfo( isliceVars(n,:), 'u' , 'Streamwise velocity' , 'm/s' , 'tttt' )
          case('v0')
            call ncinfo( isliceVars(n,:), 'v' , 'Spanwise velocity'   , 'm/s' , 'tmtt' )
          case('w0')
            call ncinfo( isliceVars(n,:), 'w' , 'Vertical velocity'   , 'm/s' , 'ttmt' )
          case('p0')
            call ncinfo( isliceVars(n,:), 'p' , 'Kinematic Pressure'  , 'm^2/s^2' , 'tttt' )
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
            print *, "Invalid slice variables name. Check namoptions setting for 'slicevars'. &
                      &The variables can be only from the list: u0,v0,w0,p0,th,qt,s1,s2,s3,s4. &
                      &There should not be any space in the string. "
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
            call ncinfo( jsliceVars(n,:), 'v' , 'Spanwise velocity '  , 'm/s' , 'tttt' )
          case('w0')
            call ncinfo( jsliceVars(n,:), 'w' , 'Vertical velocity'   , 'm/s' , 'ttmt' )
          case('p0')
            call ncinfo( jsliceVars(n,:), 'p' , 'Kinematic Pressure'  , 'm^2/s^2' , 'tttt' )
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
            print *, "Invalid slice variables name. Check namoptions setting for 'slicevars'. &
                      &The variables can be only from the list: u0,v0,w0,p0,th,qt,s1,s2,s3,s4. &
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
            call ncinfo( ksliceVars(n,:), 'w' , 'Vertical velocity'   , 'm/s' , 'tttt' )
          case('p0')
            call ncinfo( ksliceVars(n,:), 'p' , 'Kinematic Pressure'  , 'm^2/s^2' , 'tttt' )
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
            print *, "Invalid slice variables name. Check namoptions setting for 'slicevars'. &
                      &The variables can be only from the list: u0,v0,w0,p0,th,qt,s1,s2,s3,s4. &
                      &There should not be any space in the string."
            STOP 1
        end select
      end do
    end subroutine instant_ncdescription_kslice


    !! ## %% 1D parallel output creation for islice (y-direction processes only)
    subroutine instant_create_ncislice
      implicit none
      integer :: i
      logical :: has_islice
      
      ! Count how many islices are in this X-processor's range (save to module variable)
      local_nislice = 0
      has_islice = .false.
      do i = 1, nislice
!        if ( (islice(i)-1)/(itot/nprocx) == myidx) then
        if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
          local_nislice = local_nislice + 1
          has_islice = .true.
        end if
      end do
      
      ! Only processors with islices AND myidy==0 create files
      if (has_islice .and. myidy == 0) then
        islicename = 'ins_islice.xxx.xxx.nc'
        islicename(12:14) = cmyidx   ! X-processor ID
        islicename(16:18) = cexpnr  ! Experiment number

        nrecislice = 0
        ! Use local_nislice (module variable) as first dimension, jtot as global y-dimension
        call open_nc(islicename, ncidislice, nrecislice, n1=local_nislice, n2=jtot, n3=zdim)
        if (nrecislice == 0) then
          call define_nc(ncidislice, 1, slicetimeVar)
          call writestat_dims_nc(ncidislice)
          ! Add islice-specific x coordinate information
          call write_islice_xcoord_local(ncidislice, local_nislice, nislice, islice, myidx, nprocx, itot)
        end if
        call define_nc(ncidislice, nslicevars, isliceVars)
        
        ! write(*,'(A,A,A,I2,A)') '  Processor (myidx=', cmyidx, ') created islice file with ', local_nislice, ' slices'
      end if
    end subroutine instant_create_ncislice

    subroutine instant_create_ncjslice
      implicit none
      integer :: j
      logical :: has_jslice

      ! Count how many jslices are assigned to this Y-processor
      local_njslice = 0
      has_jslice = .false.
      do j = 1, njslice
!        if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
        if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
          local_njslice = local_njslice + 1
          has_jslice = .true.
        end if
      end do
      
      ! Only processors with myidx==0 create files for their myidy column
      if (has_jslice .and. myidx == 0) then
        jslicename = 'ins_jslice.xxx.xxx.nc'
        jslicename(12:14) = cmyidy   ! Y-processor ID
        jslicename(16:18) = cexpnr  ! Experiment number
        
        nrecjslice = 0
        ! Use itot as x-dimension (global), local_n as jslice-dimension, zdim as z
        call open_nc(jslicename, ncidjslice, nrecjslice, n1=itot, n2=local_njslice, n3=zdim)
        if (nrecjslice==0) then
          call define_nc(ncidjslice, 1, slicetimeVar)
          call writestat_dims_nc(ncidjslice)
          call write_jslice_ycoord_local(ncidjslice, local_njslice, njslice, jslice, myidy, nprocy, jtot)
        end if
        call define_nc(ncidjslice, nslicevars, jsliceVars)
        ! write(*,'(A,A,A,I2,A)') '  Processor (myidy=', cmyidy, ') created islice file with ', local_njslice, ' slices'
      end if
    end subroutine instant_create_ncjslice

    !! ## %% 1D parallel output creation (y-direction processes only)
    subroutine instant_create_nckslice
      implicit none

      if (myidy == 0) then
        kslicename = 'ins_kslice.xxx.xxx.nc'
        kslicename(12:14) = cmyidx   ! Only x-processor ID
        kslicename(16:18) = cexpnr  ! Experiment number

        nreckslice = 0
        ! Use kdim (nkslice) as the z-dimension, jtot as global y-dimension
        call open_nc(kslicename, ncidkslice, nreckslice, n1=xdim, n2=jtot, n3=kdim)
        if (nreckslice == 0) then
          call define_nc(ncidkslice, 1, slicetimeVar)
          call writestat_dims_nc(ncidkslice)
          ! Add kslice-specific z coordinate information
          call write_kslice_zcoord(ncidkslice, nkslice, kslice)
        end if
        call define_nc(ncidkslice, nslicevars, ksliceVars)
      end if
    end subroutine instant_create_nckslice

    subroutine instant_write_islice
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
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
!          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)  ! Global X-position
!            ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = 0.5*(um(ii_local,jb:je,kb:ke) + um(ii_local+1,jb:je,kb:ke))
          end if
        end do
        call writeoffset(ncidislice, 'u', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if

      ! v velocity
      if (present('v0')) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
!          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
!            ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = vm(ii_local,jb:je,kb:ke)
          end if
        end do
        call writeoffset(ncidislice, 'v', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      ! w velocity
      if (present('w0')) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
!          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
!            ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = wm(ii_local,jb:je,kb:ke)
          end if
        end do
        call writeoffset(ncidislice, 'w', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      ! p pressure
      if (present('p0')) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
!          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
!            ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = pres0(ii_local,jb:je,kb:ke)
          end if
        end do
        call writeoffset(ncidislice, 'p', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if

      ! temperature
      if (present('th') .and. ltempeq) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
!          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
!            ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = thlm(ii_local,jb:je,kb:ke)
          end if
        end do
        call writeoffset(ncidislice, 'thl', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      ! moisture
      if (present('qt') .and. lmoist) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
!          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
!            ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = qtm(ii_local,jb:je,kb:ke)
          end if
        end do
        call writeoffset(ncidislice, 'qt', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      ! scalars s1-s4
      if (present('s1') .and. nsv>0) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
!          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
!            ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = svm(ii_local,jb:je,kb:ke,1)
          end if
        end do
        call writeoffset(ncidislice, 's1', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      if (present('s2') .and. nsv>1) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
!          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
!            ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = svm(ii_local,jb:je,kb:ke,2)
          end if
        end do
        call writeoffset(ncidislice, 's2', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      if (present('s3') .and. nsv>2) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
!          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
!            ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = svm(ii_local,jb:je,kb:ke,3)
          end if
        end do
        call writeoffset(ncidislice, 's3', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      if (present('s4') .and. nsv>3) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
!          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
!            ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = svm(ii_local,jb:je,kb:ke,4)
          end if
        end do
        call writeoffset(ncidislice, 's4', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      deallocate(tmp_slice)
    end subroutine instant_write_islice

    subroutine instant_write_jslice
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
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = um(ib:ie, jj_local, kb:ke)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'u', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      ! v velocity (interpolated to cell centers in y-direction)
      if (present('v0')) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = 0.5*(vm(ib:ie, jj_local, kb:ke) + vm(ib:ie, jj_local+1, kb:ke))
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'v', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      ! w velocity
      if (present('w0')) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = wm(ib:ie, jj_local, kb:ke)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'w', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      ! p pressure
      if (present('p0')) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = pres0(ib:ie, jj_local, kb:ke)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'p', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      ! temperature
      if (present('th') .and. ltempeq) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = thlm(ib:ie, jj_local, kb:ke)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'thl', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      ! moisture
      if (present('qt') .and. lmoist) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = qtm(ib:ie, jj_local, kb:ke)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'qt', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      ! scalars s1-s4
      if (present('s1') .and. nsv>0) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = svm(ib:ie, jj_local, kb:ke, 1)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 's1', tmp_slice, nrecjslice, xdim,  local_njslice, zdim)
      end if

      if (present('s2') .and. nsv>1) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = svm(ib:ie, jj_local, kb:ke, 2)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 's2', tmp_slice, nrecjslice, xdim,  local_njslice, zdim)
      end if

      if (present('s3') .and. nsv>2) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = svm(ib:ie, jj_local, kb:ke, 3)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 's3', tmp_slice, nrecjslice, xdim,  local_njslice, zdim)
      end if

      if (present('s4') .and. nsv>3) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = svm(ib:ie, jj_local, kb:ke, 4)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 's4', tmp_slice, nrecjslice, xdim,  local_njslice, zdim)
      end if

      deallocate(tmp_slice)
    end subroutine instant_write_jslice

    subroutine write_islice_xcoord_local(ncid, local_n, nislice_total, islice_positions, myidx_in, nprocx_in, itot_in)
      ! Write x-coordinates for LOCAL islice positions only
      use modglobal, only : xf, xh
      implicit none
      integer, intent(in) :: ncid, local_n, nislice_total
      integer, dimension(nislice_total), intent(in) :: islice_positions
      integer, intent(in) :: myidx_in, nprocx_in, itot_in
      integer :: varid, ierr, i, local_idx
      real, allocatable :: x_islice_f(:), x_islice_h(:)
      integer, allocatable :: islice_indices(:)
      
      ! Allocate and fill x coordinates for LOCAL islices only
      allocate(x_islice_f(local_n))
      allocate(x_islice_h(local_n))
      allocate(islice_indices(local_n))
      
      local_idx = 0
      do i = 1, nislice_total
!        if ( (islice_positions(i)-1)/(itot_in/nprocx_in) == myidx_in) then
        if (islice_positions(i) >= zstart(1) .and. islice_positions(i) <= zend(1)) then
          local_idx = local_idx + 1
          x_islice_f(local_idx) = xf(islice_positions(i))  ! full level (cell center)
          x_islice_h(local_idx) = xh(islice_positions(i))  ! half level (cell edge)
          islice_indices(local_idx) = islice_positions(i)
        end if
      end do
      
      ! Write xt (full level)
      ierr = nf90_inq_varid(ncid, 'xt', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'islice_indices', islice_indices)
        ierr = nf90_put_att(ncid, varid, 'long_name', 'x-coordinate of i-slices (cell center)')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, x_islice_f)
      end if
      
      ! Write xm (half level)
      ierr = nf90_inq_varid(ncid, 'xm', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'islice_indices', islice_indices)
        ierr = nf90_put_att(ncid, varid, 'long_name', 'x-coordinate of i-slices (cell edge)')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, x_islice_h)
      end if
      
      deallocate(x_islice_f)
      deallocate(x_islice_h)
      deallocate(islice_indices)
    end subroutine write_islice_xcoord_local

    subroutine write_jslice_ycoord_local(ncid, local_n, njslice_total, jslice_positions, myidy_in, nprocy_in, jtot_in)
      ! Write y-coordinates for LOCAL jslice positions (round-robin assignment over myidy)
      use modglobal, only : yf, yh
      implicit none
      integer, intent(in) :: ncid, local_n, njslice_total
      integer, dimension(njslice_total), intent(in) :: jslice_positions
      integer, intent(in) :: myidy_in, nprocy_in, jtot_in
      integer :: varid, ierr, j, local_idy
      real, allocatable :: y_jslice_f(:), y_jslice_h(:)
      integer, allocatable :: jslice_indices(:)

      allocate(y_jslice_f(local_n))
      allocate(y_jslice_h(local_n))
      allocate(jslice_indices(local_n))
      
      local_idy = 0
      do j = 1, njslice_total
!        if (jslice_positions(j) >= ystart(2) .and. jslice_positions(j) <= yend(2)) then
        if (jslice_positions(j) >= zstart(2) .and. jslice_positions(j) <= zend(2)) then
          local_idy = local_idy + 1
          y_jslice_f(local_idy) = yf(jslice_positions(j))
          y_jslice_h(local_idy) = yh(jslice_positions(j))
          jslice_indices(local_idy) = jslice_positions(j)
        end if
      end do

      ierr = nf90_inq_varid(ncid, 'yt', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'jslice_indices', jslice_indices)
        ierr = nf90_put_att(ncid, varid, 'long_name', 'y-coordinate of j-slices (cell center)')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, y_jslice_f)
      end if

      ierr = nf90_inq_varid(ncid, 'ym', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'jslice_indices', jslice_indices)
        ierr = nf90_put_att(ncid, varid, 'long_name', 'y-coordinate of j-slices (cell edge)')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, y_jslice_h)
      end if

      deallocate(y_jslice_f)
      deallocate(y_jslice_h)
      deallocate(jslice_indices)
    end subroutine write_jslice_ycoord_local

    subroutine write_kslice_zcoord(ncid, nkslice_count, kslice_positions)
      ! Update z-coordinate to reflect kslice levels instead of full z levels
      use modglobal, only : zf, zh
      implicit none
      integer, intent(in) :: ncid, nkslice_count
      integer, dimension(nkslice_count), intent(in) :: kslice_positions
      integer :: varid, ierr, k
      real, allocatable :: z_kslice_f(:), z_kslice_h(:)
      
      ! Allocate and fill z coordinates for kslices
      allocate(z_kslice_f(nkslice_count))
      allocate(z_kslice_h(nkslice_count))
      do k = 1, nkslice_count
        z_kslice_f(k) = zf(kslice_positions(k))  ! full level (cell center)
        z_kslice_h(k) = zh(kslice_positions(k))  ! half level (cell edge)
      end do
      
      ! Write zt (full level)
      ierr = nf90_inq_varid(ncid, 'zt', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'kslice_indices', kslice_positions)
        ierr = nf90_put_att(ncid, varid, 'long_name', 'z-coordinate of k-slices (cell center)')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, z_kslice_f)
      end if
      
      ! Write zm (half level)
      ierr = nf90_inq_varid(ncid, 'zm', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'kslice_indices', kslice_positions)
        ierr = nf90_put_att(ncid, varid, 'long_name', 'z-coordinate of k-slices (cell edge)')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, z_kslice_h)
      end if
      
      deallocate(z_kslice_f)
      deallocate(z_kslice_h)
    end subroutine write_kslice_zcoord

    subroutine instant_write_kslice
      implicit none
      real, allocatable :: tmp_slice(:,:,:)
      integer :: k, kk
      
      if (nkslice == 0) return
      
      ! Allocate temporary array to hold all kslice levels: (x, y, nkslice)
      allocate(tmp_slice(ib:ie, jb:je, nkslice))
      
      ! Write time variable (only myidy==0 writes)
      if (myidy == 0) then
        call writestat_nc(ncidkslice, 'time', timee, nreckslice, .true.)
      end if
      
      ! u velocity
      if (present('u0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = um(ib:ie, jb:je, kk)
        end do
        call writeoffset(ncidkslice, 'u', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      ! v velocity
      if (present('v0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = vm(ib:ie, jb:je, kk)
        end do
        call writeoffset(ncidkslice, 'v', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      ! w velocity (interpolated to cell centers)
      if (present('w0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = 0.5 * (wm(ib:ie, jb:je, kk) + wm(ib:ie, jb:je, kk+1))
        end do
        call writeoffset(ncidkslice, 'w', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if

      ! pressure
      if (present('p0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = pres0(ib:ie, jb:je, kk)
        end do
        call writeoffset(ncidkslice, 'p', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if


      ! temperature
      if (present('th') .and. ltempeq) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = thlm(ib:ie, jb:je, kk)
        end do
        call writeoffset(ncidkslice, 'thl', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      ! moisture
      if (present('qt') .and. lmoist) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = qtm(ib:ie, jb:je, kk)
        end do
        call writeoffset(ncidkslice, 'qt', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      ! scalars s1-s4
      if (present('s1') .and. nsv>0) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie, jb:je, kk, 1)
        end do
        call writeoffset(ncidkslice, 's1', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      if (present('s2') .and. nsv>1) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie, jb:je, kk, 2)
        end do
        call writeoffset(ncidkslice, 's2', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      if (present('s3') .and. nsv>2) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie, jb:je, kk, 3)
        end do
        call writeoffset(ncidkslice, 's3', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      if (present('s4') .and. nsv>3) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie, jb:je, kk, 4)
        end do
        call writeoffset(ncidkslice, 's4', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      deallocate(tmp_slice)
    end subroutine instant_write_kslice

    logical function present(str)
      character(len=*), intent(in) :: str
      integer :: pos

      pos = index(slicevars, str)   ! returns the position of the substring str in slicevars. If not found, it returns 0.
      present = (pos > 0)
    end function present

    subroutine probe_init
      implicit none
      integer :: n, ierr, vn
      character(80) :: probe_filename
      character(2) :: varname
      real, allocatable :: xt(:), xm(:), yt(:), ym(:), zt(:), zm(:)
      
      if (.not. lprobedump) return

      ! Validate probe inputs before proceeding
      call instant_validate_out_time('probes')
      call instant_validate_output_vars(nprobevars, probevars, 'probevars', 'lprobedump')

      if (nprobe <= 0) then
        write(0, *) 'ERROR: lprobedump=.true. but nprobe=', nprobe, ' (must be > 0)'
        stop 1
      end if
      if (count(iprobe(1:nprobe) > 0) /= nprobe .or. &
          count(jprobe(1:nprobe) > 0) /= nprobe .or. &
          count(kprobe(1:nprobe) > 0) /= nprobe) then
        write(0, *) 'ERROR: nprobe=', nprobe, ' but iprobe/jprobe/kprobe do not all have', &
                    ' nprobe valid (>0) entries'
        write(0, *) 'Check that iprobe, jprobe, kprobe each have exactly nprobe positive values'
        stop 1
      end if
      ! write(*, *) 'probe output enabled for', nprobe, 'points'
      ! write(*, *) 'probevars:', trim(probevars)
      
      tsampleprobe = 0.
      
      if (myid == 0) then
          allocate(varid_vals(nprobevars))
          
          ! write(*,*) "=== Probe Output Init (Gather to Rank 0) ==="
          ! write(*,*) "Number of probes:", nprobe
          ! write(*,*) "Variables:", trim(probevars)
          
          probe_filename = 'ins_probe.xxx.nc'
          probe_filename(11:13) = cexpnr
          
          ! Create NetCDF file
          call check( nf90_create(trim(probe_filename), NF90_CLOBBER, ncid_probe) )
          
          ! Define Dimensions
          call check( nf90_def_dim(ncid_probe, 'point', nprobe, point_dimid) )
          call check( nf90_def_dim(ncid_probe, 'time', NF90_UNLIMITED, time_dimid) )
          
          ! Define Variables Dimensions
          ! Coordinates (cell center and cell edge)
          call check( nf90_def_var(ncid_probe, 'xt', NF90_FLOAT, (/point_dimid/), varid_xt) )
          call check( nf90_put_att(ncid_probe, varid_xt, 'long_name', 'x-coordinate (cell center)') )
          call check( nf90_put_att(ncid_probe, varid_xt, 'units', 'm') )
          
          call check( nf90_def_var(ncid_probe, 'xm', NF90_FLOAT, (/point_dimid/), varid_xm) )
          call check( nf90_put_att(ncid_probe, varid_xm, 'long_name', 'x-coordinate (cell edge)') )
          call check( nf90_put_att(ncid_probe, varid_xm, 'units', 'm') )
          
          call check( nf90_def_var(ncid_probe, 'yt', NF90_FLOAT, (/point_dimid/), varid_yt) )
          call check( nf90_put_att(ncid_probe, varid_yt, 'long_name', 'y-coordinate (cell center)') )
          call check( nf90_put_att(ncid_probe, varid_yt, 'units', 'm') )
          
          call check( nf90_def_var(ncid_probe, 'ym', NF90_FLOAT, (/point_dimid/), varid_ym) )
          call check( nf90_put_att(ncid_probe, varid_ym, 'long_name', 'y-coordinate (cell edge)') )
          call check( nf90_put_att(ncid_probe, varid_ym, 'units', 'm') )

          call check( nf90_def_var(ncid_probe, 'zt', NF90_FLOAT, (/point_dimid/), varid_zt) )
          call check( nf90_put_att(ncid_probe, varid_zt, 'long_name', 'z-coordinate (cell center)') )
          call check( nf90_put_att(ncid_probe, varid_zt, 'units', 'm') )
          
          call check( nf90_def_var(ncid_probe, 'zm', NF90_FLOAT, (/point_dimid/), varid_zm) )
          call check( nf90_put_att(ncid_probe, varid_zm, 'long_name', 'z-coordinate (cell edge)') )
          call check( nf90_put_att(ncid_probe, varid_zm, 'units', 'm') )

          ! Time
          call check( nf90_def_var(ncid_probe, 'time', NF90_FLOAT, (/time_dimid/), varid_time) )
          call check( nf90_put_att(ncid_probe, varid_time, 'long_name', 'Time') )
          call check( nf90_put_att(ncid_probe, varid_time, 'units', 's') )

          ! Data Variables
          do vn = 1, nprobevars
             varname = probevars(3*vn-2 : 3*vn-1)
             call check( nf90_def_var(ncid_probe, trim(get_nc_varname(varname)), NF90_FLOAT, &
                                      (/point_dimid, time_dimid/), varid_vals(vn)) )
             call add_var_atts(ncid_probe, varid_vals(vn), varname)
          end do
          
          call check( nf90_enddef(ncid_probe) )
          
          ! Write Static Coordinates immediately
          allocate(xt(nprobe), xm(nprobe))
          allocate(yt(nprobe), ym(nprobe))
          allocate(zt(nprobe), zm(nprobe))
          do n=1, nprobe
             xt(n) = xf(iprobe(n))  ! cell center
             xm(n) = xh(iprobe(n))  ! cell edge
             yt(n) = yf(jprobe(n))  ! cell center
             ym(n) = yh(jprobe(n))  ! cell edge
             zt(n) = zf(kprobe(n))  ! cell center
             zm(n) = zh(kprobe(n))  ! cell edge
          end do
          
          call check( nf90_put_var(ncid_probe, varid_xt, xt, start=(/1/), count=(/nprobe/)) )
          call check( nf90_put_var(ncid_probe, varid_xm, xm, start=(/1/), count=(/nprobe/)) )
          call check( nf90_put_var(ncid_probe, varid_yt, yt, start=(/1/), count=(/nprobe/)) )
          call check( nf90_put_var(ncid_probe, varid_ym, ym, start=(/1/), count=(/nprobe/)) )
          call check( nf90_put_var(ncid_probe, varid_zt, zt, start=(/1/), count=(/nprobe/)) )
          call check( nf90_put_var(ncid_probe, varid_zm, zm, start=(/1/), count=(/nprobe/)) )
          call check( nf90_sync(ncid_probe) )
          
          deallocate(xt, xm, yt, ym, zt, zm)
      end if
      
    end subroutine probe_init


    subroutine probe_main
      implicit none
      real, allocatable :: send_buf(:), recv_buf(:)
      integer :: vn
      character(80) :: varname
      
      if (timee < tstatstart) return
      if (.not. rk3step==3) return
      if (.not. lprobedump) return

      if (tsampleprobe > tsample) then
         
         ! Only rank 0 tracks record number
         if (myid == 0) nrecprobe = nrecprobe + 1
         
         allocate(send_buf(nprobe))
         allocate(recv_buf(nprobe))
         
         ! Always write time first (Rank 0)
         if (myid == 0) then
             call check( nf90_put_var(ncid_probe, varid_time, timee, start=(/nrecprobe/)) )
         end if

         ! Loop over all requested variables
         do vn = 1, nprobevars
            varname = probevars(3*vn-2 : 3*vn-1)
            send_buf = 0.0
            recv_buf = 0.0
            
            call gather_probe_var(varname, send_buf)
            
            ! Reduce to Rank 0
            call MPI_REDUCE(send_buf, recv_buf, nprobe, my_real, MPI_SUM, 0, comm3d, mpierr)
            
            if (myid == 0) then
                ! Use start/count to append to time dimension
               call check( nf90_put_var(ncid_probe, varid_vals(vn), recv_buf, &
                                        start=(/1, nrecprobe/), count=(/nprobe, 1/)) )
            end if
         end do
         
         if (myid == 0) then
            call check( nf90_sync(ncid_probe) )
         end if
         
         deallocate(send_buf, recv_buf)
         
         tsampleprobe = dt ! Reset timer
      else
        tsampleprobe = tsampleprobe + dt
      endif
    end subroutine probe_main


    subroutine gather_probe_var(vname, buf)
      character(len=*), intent(in) :: vname
      real, intent(out) :: buf(:)
      integer :: n, gi, gj, gk, li, lj, lk
      
      do n = 1, nprobe
         gi = iprobe(n)
         gj = jprobe(n)
         gk = kprobe(n)
         
         ! Check if local using decomp_2d zstart/zend
         if (gi >= zstart(1) .and. gi <= zend(1) .and. &
             gj >= zstart(2) .and. gj <= zend(2) .and. &
             gk >= zstart(3) .and. gk <= zend(3)) then
             
             li = gi - zstart(1) + 1
             lj = gj - zstart(2) + 1
             lk = gk - zstart(3) + 1
             
             select case (vname)
               case ('u0')
                  buf(n) = 0.5 * (um(li,lj,lk) + um(li+1,lj,lk))
               case ('v0')
                  buf(n) = 0.5 * (vm(li,lj,lk) + vm(li,lj+1,lk))
               case ('w0')
                  buf(n) = 0.5 * (wm(li,lj,lk) + wm(li,lj,lk+1))
               case ('p0')
                  buf(n) = pres0(li,lj,lk)
               case ('th')
                  if (ltempeq) buf(n) = thlm(li,lj,lk)
               case ('qt')
                  if (lmoist) buf(n) = qtm(li,lj,lk)
               case ('s1')
                  if (nsv>0) buf(n) = svm(li,lj,lk,1)
               case ('s2')
                  if (nsv>1) buf(n) = svm(li,lj,lk,2)
               case ('s3')
                  if (nsv>2) buf(n) = svm(li,lj,lk,3)
               case ('s4')
                  if (nsv>3) buf(n) = svm(li,lj,lk,4)
             end select
         endif
      end do
    end subroutine gather_probe_var

    function get_nc_varname(vname) result(ncname)
       character(len=*), intent(in) :: vname
       character(len=10) :: ncname
       select case(vname)
         case('u0')
            ncname = 'u'
         case('v0')
            ncname = 'v'
         case('w0')
            ncname = 'w'
         case('p0')
            ncname = 'p'
         case('th')
            ncname = 'thl'
         case('qt')
            ncname = 'qt'
         case('s1','s2','s3','s4')
            ncname = vname
         case default
            ncname = 'unknown'
       end select
    end function
    
    subroutine add_var_atts(ncid, vid, vname)
       integer, intent(in) :: ncid, vid
       character(len=*), intent(in) :: vname
       integer :: ierr
       select case(vname)
         case('u0')
            ierr = nf90_put_att(ncid, vid, 'long_name', 'Streamwise velocity')
            ierr = nf90_put_att(ncid, vid, 'units', 'm/s')
         case('v0')
            ierr = nf90_put_att(ncid, vid, 'long_name', 'Spanwise velocity')
            ierr = nf90_put_att(ncid, vid, 'units', 'm/s')
         case('w0')
             ierr = nf90_put_att(ncid, vid, 'long_name', 'Vertical velocity')
             ierr = nf90_put_att(ncid, vid, 'units', 'm/s')
         case('p0')
             ierr = nf90_put_att(ncid, vid, 'long_name', 'Pressure')
             ierr = nf90_put_att(ncid, vid, 'units', 'Pa')
         case('th')
             ierr = nf90_put_att(ncid, vid, 'long_name', 'Potential temperature')
             ierr = nf90_put_att(ncid, vid, 'units', 'K')
         case('qt')
             ierr = nf90_put_att(ncid, vid, 'long_name', 'Specific humidity')
             ierr = nf90_put_att(ncid, vid, 'units', 'kg/kg')
         case('s1','s2','s3','s4')
             ierr = nf90_put_att(ncid, vid, 'long_name', 'Scalar concentration')
             ierr = nf90_put_att(ncid, vid, 'units', 'g/m^3')
       end select
    end subroutine add_var_atts

    subroutine check(status)
      integer, intent ( in) :: status
      if(status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        stop 2
      end if
    end subroutine check


    subroutine instant_validate_out_time(var_name)
      implicit none
      character(len=*), intent(in)  :: var_name
      if(runtime <= tstatstart) then
        if(myid==0) then
          write(*,*) "ERROR: no instantaneous ", trim(var_name), " will be written as runtime <= tstatstart. Note that runtime &
                      &must be greater than tstatstart for wiriting ", trim(var_name), " files."
          write(*,*) "You have used runtime = ", runtime, ", tstatstart = ", tstatstart
          write(*,*) "Either correct the time settings or change all the ", trim(var_name), " writing flags to false."
          stop 1
        end if
      end if
    end subroutine instant_validate_out_time

    subroutine instant_validate_output_vars(nvars, vars, var_name, var_flag)
      ! Count variables from vars string; stop with error if none found
      implicit none
      integer,          intent(out) :: nvars
      character(len=*), intent(in)  :: vars
      character(len=*), intent(in)  :: var_name
      character(len=*), intent(in)  :: var_flag

      nvars = (LEN(trim(vars))+1)/3
      if (nvars == 0) then
        if(myid==0) then
          print *, "ERROR: no '", trim(var_name), "' entries found, although ", trim(var_flag), " are set to true !!"
          stop 1
        end if
      end if
    end subroutine instant_validate_output_vars

    subroutine instant_validate_slice_inputs(slice_enabled, n_slice, slice_array, slice_name)
      ! Validate slice arrays: check count and non-zero entries
      implicit none
      logical, intent(in) :: slice_enabled
      integer, intent(in) :: n_slice
      integer, intent(in) :: slice_array(:)
      character(len=*), intent(in) :: slice_name
      
      if (slice_enabled) then
         if (n_slice <= 0) then
            write(0, *) 'ERROR: l', trim(slice_name), 'dump=.true. but n', trim(slice_name), '=', n_slice, ' (must be > 0)'
            stop 1
         end if
         if (count(slice_array(1:n_slice) > 0) /= n_slice) then
            write(0, *) 'ERROR: n', trim(slice_name), '=', n_slice, ' but only', count(slice_array(1:n_slice) > 0), &
                         ' valid (>0) entries found in ', trim(slice_name), ' array'
            write(0, *) 'Check that ', trim(slice_name), ' has exactly n', trim(slice_name), ' positive values'
            stop 1
         end if
      end if
    end subroutine instant_validate_slice_inputs

end module instant