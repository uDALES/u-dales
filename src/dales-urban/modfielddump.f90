!> \file modfielddump.f90
!!  Dumps 3D fields of several variables
!>
!!  Dumps 3D fields of several variables Written to wb*.myid.expnr
!! If netcdf is true, this module leads the fielddump.myid.expnr.nc output
!!  \author Jasper Tomas, TU Delft Match 31 2014
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!! Possibility to write tecplot (formatted!!) file 
!
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!
module modfielddump

  use modglobal, only : longint
  use modfields, only : ncname
  implicit none
  private
  PUBLIC :: initfielddump, fielddump,exitfielddump
  save
  !NetCDF variables
  integer :: ncid,nrec = 0
!  real, pointer :: point
  type domainptr
    real, pointer :: point(:,:,:)
  end type domainptr
  type(domainptr), dimension(30) :: pfields

  character(80) :: fname = 'fielddump.xxx.xxx.nc'
  !dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname
  
  real    :: dtav,tnext
  integer :: klow,khigh,n,nvar
  logical :: lfielddump= .true.  !< switch to enable the fielddump (on/off)
  logical :: ldiracc   = .false. !< switch for doing direct access writing (on/off)
  logical :: lbinary   = .false. ! 

contains
 !> Initializing fielddump. Read out the namelist, initializing the variables
  subroutine initfielddump
    use modmpi,   only   :myid,my_real,mpierr,comm3d,mpi_logical,mpi_integer,cmyid, mpi_character
    use modglobal,only   :imax,jmax,kmax,cexpnr,ifnamopt,fname_options,dtmax,dtav_glob,kb,ke, ladaptive,dt_lim,btime,nsv,fieldvars,ib,ie,jb,je,kb,ke
    use modstat_nc,only  : open_nc, define_nc,ncinfo,writestat_dims_nc
    use modfields, only  : u0,v0,w0,thl0,sv0,ql0,qt0,pres0
    implicit none
    integer :: ierr

  !  type(domainptr), dimension(nvar) :: pfields

    namelist/NAMFIELDDUMP/ &
         dtav,lfielddump,ldiracc,lbinary,klow,khigh

    write(*,*), '(LEN(trim(fieldvars))+1)/3', (LEN(trim(fieldvars))+1)/3

    nvar = (LEN(trim(fieldvars))+1)/3
    
    if (nvar == 0) then
      lfielddump = .false.
      print *, 'empty fieldvars therefore lfielddump = .false. and no instantaneous fields outputted'
      return
    else
      allocate(ncname(nvar,4)) 
    end if

    dtav=dtav_glob
    klow=kb
    khigh=ke
    if(myid==0)then
       open(ifnamopt,file=fname_options,status='old',iostat=ierr)
       read (ifnamopt,NAMFIELDDUMP,iostat=ierr)
       if (ierr > 0) then
          print *, 'Problem in namoptions NAMFIELDDUMP'
          print *, 'iostat error: ', ierr
          stop 'ERROR: Problem in namoptions NAMFIELDDUMP'
       endif
       write(6 ,NAMFIELDDUMP)
       close(ifnamopt)
    end if

    call MPI_BCAST(klow        ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(khigh       ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(dtav        ,1,MY_REAL    ,0,comm3d,ierr)
    call MPI_BCAST(lfielddump  ,1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(ldiracc     ,1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(lbinary     ,1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(ncname     ,80,MPI_CHARACTER,0,comm3d,mpierr) 
    call MPI_BCAST(nvar       ,1,MPI_INTEGER,0,comm3d,mpierr)

    tnext      = dtav   +btime
    !    dt_lim = min(dt_lim,tnext)

    if(.not.(lfielddump)) return
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then !tg3315 this is not relevant to these fields, dtav is stats sampling interval...
       stop 'dtav should be a integer multiple of dtmax'
    end if

    fname(11:13) = cmyid
    fname(15:17) = cexpnr
    call ncinfo(tncname(1,:),'time','Time','s','time')

    ! tg3315 reads in fields specified by fieldvars
    do n=1,nvar
      select case(fieldvars(3*n-2:3*n-1))
        case('u0')
        call ncinfo(ncname( n,:),'u','West-East velocity','m/s','mttt')
        pfields(n)%point => u0(ib:ie,jb:je,kb:ke)
        case('v0')
        call ncinfo(ncname( n,:),'v','South-North velocity','m/s','tmtt')
        pfields(n)%point => v0(ib:ie,jb:je,kb:ke)  
        case('w0')
        call ncinfo(ncname( n,:),'w','Vertical velocity','m/s','ttmt')
        pfields(n)%point => w0(ib:ie,jb:je,kb:ke)  
        case('th')
        call ncinfo(ncname( n,:),'thl','Liquid water potential temperature','K','tttt')
        pfields(n)%point => thl0(ib:ie,jb:je,kb:ke)  
        case('ql')
        call ncinfo(ncname( n,:),'ql','Liquid water mixing ratio','1e-5kg/kg','tttt')
        pfields(n)%point => ql0(ib:ie,jb:je,kb:ke) 
        case('qt')
        call ncinfo(ncname( n,:),'qt','Total water mixing ratio','1e-5kg/kg','tttt')
        pfields(n)%point => qt0(ib:ie,jb:je,kb:ke)
        case('s1')
        call ncinfo(ncname( n,:),'sca1','scalar 1','M','tttt')
        pfields(n)%point => sv0(ib:ie,jb:je,kb:ke,1)
        case('s2')
        call ncinfo(ncname( n,:),'sca2','scalar 2','M','tttt')
        pfields(n)%point => sv0(ib:ie,jb:je,kb:ke,2)
        case('s3')
        call ncinfo(ncname( n,:),'sca3','scalar 3','M','tttt')
        pfields(n)%point => sv0(ib:ie,jb:je,kb:ke,3)
        case('s4')
        call ncinfo(ncname( n,:),'sca4','scalar 4','M','tttt')
        pfields(n)%point => sv0(ib:ie,jb:je,kb:ke,4)
        case('s5')
        call ncinfo(ncname( n,:),'sca5','scalar 5','M','tttt')
        pfields(n)%point => sv0(ib:ie,jb:je,kb:ke,5)
        case('p0')
        call ncinfo(ncname( n,:),'pres','pressure field','M','tttt')
        pfields(n)%point => pres0(ib:ie,jb:je,kb:ke)
        case default
        call ncinfo(ncname( n,:),'u','West-East velocity','m/s','mttt')
        pfields(n)%point => u0(ib:ie,jb:je,kb:ke)
      end select
    end do
 
    call open_nc( fname, ncid, nrec, n1=imax, n2=jmax, n3=khigh-klow+1)
    if (nrec==0) then
      call define_nc( ncid, 1, tncname)
      call writestat_dims_nc(ncid)  
    end if
    call define_nc( ncid, nvar, ncname)
      
           call open_nc( fname, ncid, nrec, n1=imax, n2=jmax, n3=khigh-klow+1)
!        call open_nc( fname, ncid, nrec, n1=imax+2, n2=jmax+2, n3=khigh-klow+1)  !if want to print ghostcells

        if (nrec==0) then
          write(*,*) "calling define_nc"
          call define_nc( ncid, 1, tncname)
          write(*,*) "calling writestat_dims_nc"
          call writestat_dims_nc(ncid)  
        end if
       call define_nc( ncid, nvar, ncname)

  end subroutine initfielddump

  !> Do fielddump. Collect data to truncated (2 byte) integers, and write them to file
  subroutine fielddump
    use modfields, only : u0,v0,w0,thl0,qt0,ql0,sv0  !ILS13 21.04.2015 changed to u0 from um  etc
    use modsurfdata,only : thls,qts,thvs
    use modglobal, only : ib,ie,ih,jb,je,jh,ke,kb,kh,rk3step,timee,dt_lim,cexpnr,ifoutput,imax,jmax,&
                          tfielddump, tnextfielddump,nsv
    !use modmpi,    only : myid,cmyid
    !use modsubgriddata, only : ekm,sbshr
    use modstat_nc, only : writestat_nc
    use modmpi, only : myid,cmyid
    implicit none

    real, allocatable :: vars(:,:,:,:)
    integer i,j,k
    integer :: writecounter = 1
 
    if (.not. (timee>=tnextfielddump)) return 

    if (.not. lfielddump) return

    if (rk3step/=3) return

    tnextfielddump=tnextfielddump+tfielddump

   allocate(vars(ib:ie,jb:je,kb:ke,nvar))

    do n=1,nvar
      vars(ib:ie,jb:je,kb:ke,n) = pfields(n)%point
    end do

   call writestat_nc(ncid,1,tncname,(/timee/),nrec,.true.)
         
   call writestat_nc(ncid,nvar,ncname,vars,nrec,imax,jmax,khigh-klow+1)
   deallocate(vars)

  end subroutine fielddump

  !> Clean up when leaving the run
  subroutine exitfielddump
      use modstat_nc, only : exitstat_nc
    implicit none

       if (lfielddump) call exitstat_nc(ncid)
  end subroutine exitfielddump

end module modfielddump
