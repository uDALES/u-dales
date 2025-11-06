module instant_slice
  use modglobal,  only : cexpnr, rk3step, ltempeq, lmoist, nsv, tsample, dt, timee, runtime, &
                         slicevars, lislicedump, islice, ljslicedump, jslice, lkslicedump, kslice, &
                         ib, ie, jb, je, kb, ke, dzfi, dzh, &
                         timee, tstatstart
  use modfields,  only : um, vm, wm, thlm, qtm, svm
  use modmpi,     only : myid, cmyidx, cmyidy
  use decomp_2d,  only : zstart, zend
  use modstat_nc, only : ncinfo, open_nc, define_nc, writestat_dims_nc, writestat_nc
  implicit none
  private
  public :: instant_init, instant_main
  save

  integer :: xdim, ydim, zdim
  real    :: tsampleslice

  integer :: isliceloc    ! local islice on core
  logical :: islicerank   ! cpu that islice is on
  integer :: jsliceloc    ! local jslice on core
  logical :: jslicerank   ! cpu that jslice is on
  
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
      tsampleslice = 0.
      call ncinfo(slicetimeVar( 1,:), 'time', 'Time', 's', 'time')

      if (lislicedump) then
        allocate(isliceVars(nslicevars,4))
        call instant_ncdescription_islice
        call instant_create_ncislice    !> Generate sliced NetCDF: islice.xxx.xxx.nc
        deallocate(isliceVars)
      end if
      
      if (ljslicedump) then
        allocate(jsliceVars(nslicevars,4))
        call instant_ncdescription_jslice
        call instant_create_ncjslice    !> Generate sliced NetCDF: jslice.xxx.xxx.nc
        deallocate(jsliceVars)
      end if

      if (lkslicedump) then
        allocate(ksliceVars(nslicevars,4))
        call instant_ncdescription_kslice
        call instant_create_nckslice    !> Generate sliced NetCDF: kslice.xxx.xxx.nc
        deallocate(ksliceVars)
      end if
    end subroutine instant_init


    subroutine instant_main
      implicit none
      if (timee < tstatstart) return
      if (.not. rk3step==3)  return
      if(.not.(lislicedump .or. ljslicedump .or. lkslicedump)) return

      if (tsampleslice > tsample) then
        if (lislicedump .and. islicerank) call instant_write_islice
        if (ljslicedump .and. jslicerank) call instant_write_jslice
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
            call ncinfo( isliceVars(n,:), 'u' , 'Streamwise velocity' , 'm/s' , '0ttt' )
          case('v0')
            call ncinfo( isliceVars(n,:), 'v' , 'Spanwise velocity'   , 'm/s' , '0mtt' )
          case('w0')
            call ncinfo( isliceVars(n,:), 'w' , 'Vertical velocity'   , 'm/s' , '0tmt' )
          case('th')
            if (ltempeq) call ncinfo( isliceVars(n,:), 'thl' , 'Potential temperature' , 'K'     , '0ttt' )
          case('qt')
            if (lmoist)  call ncinfo( isliceVars(n,:), 'qt'  , 'Specific humidity'     , 'kg/kg' , '0ttt' )
          case('s1')
            if (nsv>0)   call ncinfo( isliceVars(n,:), 's1'  , 'Concentration field 1' , 'g/m^3' , '0ttt' )
          case('s2')
            if (nsv>1)   call ncinfo( isliceVars(n,:), 's2'  , 'Concentration field 2' , 'g/m^3' , '0ttt' )
          case('s3')
            if (nsv>2)   call ncinfo( isliceVars(n,:), 's3'  , 'Concentration field 3' , 'g/m^3' , '0ttt' )
          case('s4')
            if (nsv>3)   call ncinfo( isliceVars(n,:), 's4'  , 'Concentration field 4' , 'g/m^3' , '0ttt' )
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
            call ncinfo( jsliceVars(n,:), 'u' , 'Streamwise velocity' , 'm/s' , 'm0tt' )
          case('v0')
            call ncinfo( jsliceVars(n,:), 'v' , 'Spanwise velocity '  , 'm/s' , 't0tt' )
          case('w0')
            call ncinfo( jsliceVars(n,:), 'w' , 'Vertical velocity'   , 'm/s' , 't0mt' )
          case('th')
            if (ltempeq) call ncinfo( jsliceVars(n,:), 'thl' , 'Potential temperature' , 'K'     , 't0tt' )
          case('qt')
            if (lmoist)  call ncinfo( jsliceVars(n,:), 'qt'  , 'Specific humidity'     , 'kg/kg' , 't0tt' )
          case('s1')
            if (nsv>0)   call ncinfo( jsliceVars(n,:), 's1'  , 'Concentration field 1' , 'g/m^3' , 't0tt' )
          case('s2')
            if (nsv>1)   call ncinfo( jsliceVars(n,:), 's2'  , 'Concentration field 2' , 'g/m^3' , 't0tt' )
          case('s3')
            if (nsv>2)   call ncinfo( jsliceVars(n,:), 's3'  , 'Concentration field 3' , 'g/m^3' , 't0tt' )
          case('s4')
            if (nsv>3)   call ncinfo( jsliceVars(n,:), 's4'  , 'Concentration field 4' , 'g/m^3' , 't0tt' )
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
            call ncinfo( ksliceVars(n,:), 'u' , 'Streamwise velocity' , 'm/s' , 'mt0t' )
          case('v0')
            call ncinfo( ksliceVars(n,:), 'v' , 'Spanwise velocity'   , 'm/s' , 'tm0t' )
          case('w0')
            call ncinfo( ksliceVars(n,:), 'w' , 'Vertical velocity'   , 'm/s' , 'tt0t' )
          case('th')
            if (ltempeq) call ncinfo( ksliceVars( n,:), 'thl' , 'Potential temperature' , 'K'     , 'tt0t' )
          case('qt')
            if (lmoist)  call ncinfo( ksliceVars( n,:), 'qt'  , 'Specific humidity'     , 'kg/kg' , 'tt0t' )
          case('s1')
            if (nsv>0)   call ncinfo( ksliceVars( n,:), 's1'  , 'Concentration field 1' , 'g/m^3' , 'tt0t' )
          case('s2')
            if (nsv>1)   call ncinfo( ksliceVars( n,:), 's2'  , 'Concentration field 2' , 'g/m^3' , 'tt0t' )
          case('s3')
            if (nsv>2)   call ncinfo( ksliceVars( n,:), 's3'  , 'Concentration field 3' , 'g/m^3' , 'tt0t' )
          case('s4')
            if (nsv>3)   call ncinfo( ksliceVars( n,:), 's4'  , 'Concentration field 4' , 'g/m^3' , 'tt0t' )
          case default
            print *, "Invalid slice variables name. Check namoptions variable 'slicevars'. &
                      &There should not be any space in the string."
            STOP 1
        end select
      end do
    end subroutine instant_ncdescription_kslice


    subroutine instant_create_ncislice
      implicit none
      islicename = 'islice.xxx.xxx.xxx.nc'
      islicename(8:10)  = cmyidx
      islicename(12:14) = cmyidy
      islicename(16:18) = cexpnr

      if ((islice >= zstart(1)) .and. (islice <= zend(1))) then
        islicerank = .true.
        isliceloc = islice - zstart(1) + 1
      else
        islicerank = .false.
      end if

      nrecislice = 0
      if (islicerank) then
        call open_nc(islicename, ncidislice, nrecislice, n2=ydim, n3=zdim)
        if (nrecislice==0) then
          call define_nc(ncidislice, 1, slicetimeVar)
          call writestat_dims_nc(ncidislice)
        end if
        call define_nc(ncidislice, nslicevars, isliceVars)
      end if
    end subroutine instant_create_ncislice

    subroutine instant_create_ncjslice
      implicit none
      jslicename = 'jslice.xxx.xxx.xxx.nc'
      jslicename(8:10)  = cmyidx
      jslicename(12:14) = cmyidy
      jslicename(16:18) = cexpnr

      if ((jslice >= zstart(2)) .and. (jslice <= zend(2))) then
        jslicerank = .true.
        jsliceloc = jslice - zstart(2) + 1
      else
        jslicerank = .false.
      end if

      nrecjslice = 0
      if (jslicerank) then
         call open_nc(jslicename, ncidjslice, nrecjslice, n1=xdim, n3=zdim)
         if (nrecjslice==0) then
            call define_nc(ncidjslice, 1, slicetimeVar)
            call writestat_dims_nc(ncidjslice)
         end if
         call define_nc(ncidjslice, nslicevars, jsliceVars)
      end if
    end subroutine instant_create_ncjslice

    subroutine instant_create_nckslice
      implicit none
      kslicename = 'kslice.xxx.xxx.xxx.nc'
      kslicename(8:10)  = cmyidx
      kslicename(12:14) = cmyidy
      kslicename(16:18) = cexpnr

      nreckslice = 0
      call open_nc(kslicename, ncidkslice, nreckslice, n1=xdim, n2=ydim)
      if (nreckslice==0) then
        call define_nc(ncidkslice, 1, slicetimeVar)
        call writestat_dims_nc(ncidkslice)
      end if
      call define_nc(ncidkslice, nslicevars, ksliceVars)
    end subroutine instant_create_nckslice


    subroutine instant_write_islice
      implicit none
      call writestat_nc(ncidislice, 'time', timee, nrecislice, .true.)
      if (present('u0')) call writestat_nc(ncidislice, 'u' , 0.5*(um(isliceloc,jb:je,kb:ke)+um(isliceloc+1,jb:je,kb:ke)) , nrecislice, ydim, zdim)
      if (present('v0')) call writestat_nc(ncidislice, 'v' , vm(isliceloc,jb:je,kb:ke)                                   , nrecislice, ydim, zdim)
      if (present('w0')) call writestat_nc(ncidislice, 'w' , wm(isliceloc,jb:je,kb:ke)                                   , nrecislice, ydim, zdim)
      if (present('th') .and. ltempeq) call writestat_nc(ncidislice, 'thl' , thlm(isliceloc,jb:je,kb:ke)  , nrecislice, ydim, zdim)
      if (present('qt') .and. lmoist)  call writestat_nc(ncidislice, 'qt'  , qtm(isliceloc,jb:je,kb:ke)   , nrecislice, ydim, zdim)
      if (present('s1') .and. nsv>0)   call writestat_nc(ncidislice, 's1'  , svm(isliceloc,jb:je,kb:ke,1) , nrecislice, ydim, zdim)
      if (present('s2') .and. nsv>1)   call writestat_nc(ncidislice, 's2'  , svm(isliceloc,jb:je,kb:ke,2) , nrecislice, ydim, zdim)
      if (present('s3') .and. nsv>2)   call writestat_nc(ncidislice, 's3'  , svm(isliceloc,jb:je,kb:ke,3) , nrecislice, ydim, zdim)
      if (present('s4') .and. nsv>3)   call writestat_nc(ncidislice, 's4'  , svm(isliceloc,jb:je,kb:ke,4) , nrecislice, ydim, zdim)
    end subroutine instant_write_islice

    subroutine instant_write_jslice
      implicit none
      call writestat_nc(ncidjslice, 'time', timee, nrecjslice, .true.)
      if (present('u0')) call writestat_nc(ncidjslice, 'u' , um(ib:ie,jsliceloc,kb:ke)                                   , nrecjslice, xdim, zdim)
      if (present('v0')) call writestat_nc(ncidjslice, 'v' , 0.5*(vm(ib:ie,jsliceloc,kb:ke)+vm(ib:ie,jsliceloc+1,kb:ke)) , nrecjslice, xdim, zdim)
      if (present('w0')) call writestat_nc(ncidjslice, 'w' , wm(ib:ie,jsliceloc,kb:ke)                                   , nrecjslice, xdim, zdim)
      if (present('th') .and. ltempeq) call writestat_nc(ncidjslice, 'thl' , thlm(ib:ie,jsliceloc,kb:ke)  , nrecjslice, xdim, zdim)
      if (present('qt') .and. lmoist)  call writestat_nc(ncidjslice, 'qt'  , qtm(ib:ie,jsliceloc,kb:ke)   , nrecjslice, xdim, zdim)
      if (present('s1') .and. nsv>0)   call writestat_nc(ncidjslice, 's1'  , svm(ib:ie,jsliceloc,kb:ke,1) , nrecjslice, xdim, zdim)
      if (present('s2') .and. nsv>1)   call writestat_nc(ncidjslice, 's2'  , svm(ib:ie,jsliceloc,kb:ke,2) , nrecjslice, xdim, zdim)
      if (present('s3') .and. nsv>2)   call writestat_nc(ncidjslice, 's3'  , svm(ib:ie,jsliceloc,kb:ke,3) , nrecjslice, xdim, zdim)
      if (present('s4') .and. nsv>3)   call writestat_nc(ncidjslice, 's4'  , svm(ib:ie,jsliceloc,kb:ke,4) , nrecjslice, xdim, zdim)
    end subroutine instant_write_jslice

    subroutine instant_write_kslice
      implicit none
      call writestat_nc(ncidkslice, 'time', timee, nreckslice, .true.)
      if (present('u0')) call writestat_nc(ncidkslice, 'u' , um(ib:ie,jb:je,kslice)                                , nreckslice, xdim, ydim)
      if (present('v0')) call writestat_nc(ncidkslice, 'v' , vm(ib:ie,jb:je,kslice)                                , nreckslice, xdim, ydim)
      if (present('w0')) call writestat_nc(ncidkslice, 'w' , 0.5*(wm(ib:ie,jb:je,kslice)+wm(ib:ie,jb:je,kslice+1)) , nreckslice, xdim, ydim)
      if (present('th') .and. ltempeq) call writestat_nc(ncidkslice, 'thl' , thlm(ib:ie,jb:je,kslice)  , nreckslice, xdim, ydim)
      if (present('qt') .and. lmoist)  call writestat_nc(ncidkslice, 'qt'  , qtm(ib:ie,jb:je,kslice)   , nreckslice, xdim, ydim)
      if (present('s1') .and. nsv>0)   call writestat_nc(ncidkslice, 's1'  , svm(ib:ie,jb:je,kslice,1) , nreckslice, xdim, ydim)
      if (present('s2') .and. nsv>1)   call writestat_nc(ncidkslice, 's2'  , svm(ib:ie,jb:je,kslice,2) , nreckslice, xdim, ydim)
      if (present('s3') .and. nsv>2)   call writestat_nc(ncidkslice, 's3'  , svm(ib:ie,jb:je,kslice,3) , nreckslice, xdim, ydim)
      if (present('s4') .and. nsv>3)   call writestat_nc(ncidkslice, 's4'  , svm(ib:ie,jb:je,kslice,4) , nreckslice, xdim, ydim)
    end subroutine instant_write_kslice

    logical function present(str)
      character(len=*), intent(in) :: str
      integer :: pos

      pos = index(slicevars, str)   ! returns the position of the substring str in slicevars. If not found, it returns 0.
      present = (pos > 0)
    end function present
end module instant_slice