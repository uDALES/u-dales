!> \file modslicedump.f90
!!  Dumps slices of several variables
module modslicedump

  use modglobal, only : longint
  use modfields, only : ncname_slice
  implicit none
  private
  PUBLIC :: initslicedump, slicedump,exitslicedump
  save
  !NetCDF variables
  integer :: ncid,ncid1,ncid2,nrec = 0
!  real, pointer :: point
  type domainptr
    real, pointer :: point(:,:)
  end type domainptr
  type(domainptr), dimension(30) :: pfields

  character(80) :: fname = 'slicedump.xxx.xxx.xxx.nc'
  !dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname

  integer :: ilow,ihigh,jlow,jhigh,klow,khigh,nvar
  logical :: ldiracc   = .false. !< switch for doing direct access writing (on/off)
  logical :: lbinary   = .false. !
  logical :: lhalos

contains
 !> Initializing slicedump. Read out the namelist, initializing the variables
  subroutine initslicedump
    use modmpi,   only   :myid,my_real,mpierr,comm3d,mpi_logical,mpi_integer,cmyidx,cmyidy,mpi_character
    use modglobal,only   :imax,jmax,kmax,imax1,jmax1,kmax1,imax2,jmax2,kmax2,cexpnr,ifnamopt,fname_options,dtmax,kb,ke, ladaptive,dt_lim,btime,nsv,slicevars,ib,ie,jb,je,kb,ke, ih,jh,lfielddump, lslice,jslice,ktot,kh
    use modstat_nc,only  : open_nc, define_nc,ncinfo,writestat_dims_nc
    use modfields, only  : u0,v0,w0,thl0,sv0,ql0,qt0,pres0,div,dudx,dvdy,dwdz,ru,rv,rw,tau_x, tau_y, tau_z, thl_flux, up,wp,vp,thlp
    use modpois, only : p, pup,pvp,pwp, rhs, dpupdx, dpvpdy, dpwpdz, xyzrt, Fxy, Fxyz
    use modibm, only : mask_u, mask_v, mask_w, mask_c
    use decomp_2d, only: zstart
    implicit none
    integer :: ierr, n, jeg,jbg,jslicelocal
    real, allocatable, target :: thl0_slice(:,:)
    allocate(thl0_slice(ib:ie,kb:ke)) ; thl0_slice=0.
    !type(domainptr), dimension(nvar) :: pfields
    nvar = (LEN(trim(slicevars))+1)/3
    if (nvar == 0) then
      lslice = .false.
      print *, 'empty slicevars therefore lfielddump = .false. and no instantaneous fields outputted'
      return
    else
      allocate(ncname_slice(nvar,4))
    end if

      ilow = ib
      ihigh = ie
      jlow = jb
      jhigh = je
      klow=kb
      khigh=ke



!ils13 13.08.18: why is this broadcast, doesn't every processor do it anyway?
    call MPI_BCAST(klow        ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(khigh       ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(lslice  ,1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(ldiracc     ,1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(lbinary     ,1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(ncname_slice     ,80,MPI_CHARACTER,0,comm3d,mpierr)
    call MPI_BCAST(nvar       ,1,MPI_INTEGER,0,comm3d,mpierr)

    !    dt_lim = min(dt_lim,tnext)

    if(.not.(lslice)) return

    fname(11:13) = cmyidx
    fname(15:17) = cmyidy
    fname(19:21) = cexpnr
    call ncinfo(tncname(1,:),'time','Time','s','time')

    ! tg3315 reads in fields specified by fieldvars
      if (lslice) then
        jbg = jb + zstart(2)-1
        jeg = je + zstart(2)-1
        jslicelocal = jslice - zstart(2)+1
      end if
      if (jbg <= jslice .and. jeg >= jslice) then
        thl0_slice = thl0(ib:ie,jslice,kb:ke)
      do n=1,nvar
        select case(slicevars(3*n-2:3*n-1))
        case('ts')
          call ncinfo(ncname_slice( n,:),'thl_slice','Liquid water potential temperature slice','K','t0tt')
            pfields(n)%point => thl0_slice(ib:ie,kb:ke)
        end select
      end do
    end if


    !call ncinfo(ncname( n,:),'u','West-East velocity','m/s','mttt')
    !call open_nc( fname, ncid, nrec, n1=imax+2, n2=jmax+2, n3=khigh-klow+1)
    call open_nc( fname, ncid, nrec, n1=ihigh-ilow+1, n3=khigh-klow+1)
    if (nrec==0) then
      call define_nc( ncid, 1, tncname)
      call writestat_dims_nc(ncid)
    end if

    call define_nc( ncid, nvar, ncname_slice)
    !call open_nc( fname, ncid, nrec, n1=imax+2, n2=jmax+2, n3=khigh-klow+1)
    call open_nc( fname, ncid, nrec, n1=ihigh-ilow+1, n3=khigh-klow+1)

    if (nrec==0) then
      call define_nc( ncid, 1, tncname)
      call writestat_dims_nc(ncid)
    end if
    call define_nc( ncid, nvar, ncname_slice)

  end subroutine initslicedump

  !> Do fielddump. Collect data to truncated (2 byte) integers, and write them to file
  subroutine slicedump
    use modfields, only : u0,v0,w0,thl0,qt0,ql0,sv0,pres0,u01,u02,u0h,um,div,dudx,dvdy,dwdz,tau_x  !ILS13 21.04.2015 changed to u0 from um  etc
    use modsurfdata,only : thls,qts,thvs
    use modglobal, only : ib,ie,ih,jb,je,jh,ke,kb,kh,rk3step,timee,dt_lim,cexpnr,ifoutput,imax,jmax,&
                          tfielddump, tnextfielddump,nsv, lfielddump,lslice,jslice, ktot,fieldvars, imax1,jmax1,kmax1,imax2,jmax2,kmax2,rk3step,dyi,dxfi,dzhi
    !use modmpi,    only : myid,cmyid
    !use modsubgriddata, only : ekm,sbshr
    use modstat_nc, only : writestat_nc
    use modmpi, only : myid,cmyid
    use decomp_2d, only: zstart
    implicit none
    real, allocatable :: vars(:,:,:), vars1(:,:,:,:), vars2(:,:,:,:)
    integer i,j,k,n, jbg,jeg,jslicelocal
    integer :: writecounter = 1
    integer, parameter :: sz = 64
    real, dimension(:,:), allocatable :: ones_array
    allocate(ones_array(sz, sz))
    ones_array = 1.0
    if (.not. ((timee>=tnextfielddump) .or. (rk3step==0))) return
    if (.not. lslice) return
    if (rk3step/=3 .and. rk3step/=0) return

    if (rk3step == 3) tnextfielddump=tnextfielddump+tfielddump
    jbg = jb + zstart(2)-1
    jeg = je + zstart(2)-1
    jslicelocal = jslice - zstart(2)+1
      if (jbg <= jslice .and. jeg >= jslice) then
        allocate(vars(ib:ie,kb:ke,nvar)); vars = 0;
        do n=1,nvar
          !vars(ib:ie,kb:ke,n) = pfields(n)%point
          vars(ib:ie,kb:ke,n) = ones_array
        end do
      end if

   call writestat_nc(ncid,1,tncname,(/timee/),nrec,.true.)
   !call writestat_nc(ncid,nvar,ncname,vars,nrec,imax+2,jmax+2,khigh-klow+1)

   call writestat_nc(ncid,nvar,ncname_slice,vars,nrec,ihigh-ilow+1,khigh-klow+1)

   deallocate(vars)

 end subroutine slicedump

  !> Clean up when leaving the run
  subroutine exitslicedump
      use modglobal, only : lslice
      use modstat_nc, only : exitstat_nc
    implicit none

       if (lslice) call exitstat_nc(ncid)
  end subroutine exitslicedump

! Remove all data from cores not needed for slice.
 subroutine tdumpslice
    use netcdf
    use modmpi,   only : cmyidx,cmyidy
    use modglobal,only : cexpnr, jb, je, jslice, lslice
    use decomp_2d, only: zstart
    implicit none
    integer :: ncid, status, jbg,jeg, nf90_inq_nvars
    ! Specify the name of your NetCDF file
    character(80) :: filename = 'tdump.xxx.xxx.xxx.nc'
    jbg = jb + zstart(2)-1
    jeg = je + zstart(2)-1
    if (.not. lslice) return
    if ((jbg <= jslice .and. jeg >= jslice)) return !   We don't edit the data on the course we want.
    filename(11:13) = cmyidx
    filename(15:17) = cmyidy
    filename(19:21) = cexpnr
    ! Open the NetCDF file in write mode
    status = nf90_open(filename, nf90_write, ncid)

    if (status /= nf90_noerr) then
        print *, "Error opening NetCDF file"
        stop
    end if

    ! Call a subroutine to clear the data
    call clear_data(ncid)

    ! Close the NetCDF file
    status = nf90_close(ncid)

    if (status /= nf90_noerr) then
        print *, "Error closing NetCDF file"
        stop
    end if
 end subroutine tdumpslice

  subroutine clear_data(ncid)
     ! Remove data from all variables in the NetCDF file
     USE netcdf
     integer, intent(in) :: ncid
     integer :: varid, dimid, status, nvarstdump,i, ndims, nf90_inq_nvars, nf90_inq_varndims, nf90_put_vara
     INTEGER, ALLOCATABLE :: start(:), count(:)
     ! Inquire about the total number of variables
    status = nf90_inq_nvars(ncid, nvarstdump)

    !F (status /= nf90_noerr) THEN
      !  PRINT *, "Error inquiring about the number of variables"
    !    STOP
    !END IF

    ! Loop through variables
    DO i = 1, nvarstdump
        ! Inquire about variable ID
        status = nf90_inq_varid(ncid, 'thlt', varid)

        !IF (status /= nf90_noerr) THEN
        !    PRINT *, "Error inquiring about variable ID"
        !    STOP
        !END IF

        ! Inquire about variable dimensions
        status = nf90_inq_varndims(ncid, varid, ndims)

        !IF (status /= nf90_noerr) THEN
        !    PRINT *, "Error inquiring about variable dimensions"
        !    STOP
        !END IF

        ! Define a hyperslab to write an empty array

        ALLOCATE(start(ndims), count(ndims))

        start = 1
        count = 0

        ! Write the empty array to the variable
        status = nf90_put_vara_real(ncid, varid, start, count, [1])

        !IF (status /= nf90_noerr) THEN
        !    PRINT *, "Error writing empty array to variable"
        !    STOP
        !END IF
    END DO

    ! Close the NetCDF file
    status = nf90_close(ncid)

    !IF (status /= nf90_noerr) THEN
    !    PRINT *, "Error closing NetCDF file"
    !    STOP
    !END IF

   end subroutine clear_data

end module modslicedump
