!> \file modfielddump.f90
!!  Dumps 3D fields of several variables

!>
!!  Dumps 3D fields of several variables
!>
!!  Dumps 3D fields of several variables Written to wb*.myid.expnr
!! If netcdf is true, this module leads the fielddump.myid.expnr.nc output
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
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


implicit none
private
PUBLIC :: initfielddump, fielddump,exitfielddump
save
!NetCDF variables
  integer,parameter :: nvar = 6
  integer :: ncid,nrec = 0
  character(80) :: fname = 'fielddump.xxx.xxx.nc'
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav,tnext
  integer :: klow,khigh
  logical :: lfielddump= .false. !< switch to enable the fielddump (on/off)
  logical :: ldiracc   = .false. !< switch for doing direct access writing (on/off)

contains
!> Initializing fielddump. Read out the namelist, initializing the variables
  subroutine initfielddump
    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_logical,mpi_integer,cmyid
    use modglobal,only :imax,jmax,kmax,cexpnr,ifnamopt,fname_options,dtmax,dtav_glob,kmax, ladaptive,dt_lim,btime
    use modstat_nc,only : lnetcdf,open_nc, define_nc, redefine_nc,ncinfo,writestat_dims_nc
    implicit none
    integer :: ierr


    namelist/NAMFIELDDUMP/ &
    dtav,lfielddump,ldiracc,klow,khigh

    dtav=dtav_glob
    klow=1
    khigh=kmax
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMFIELDDUMP,iostat=ierr)
      write(6 ,NAMFIELDDUMP)
      close(ifnamopt)
    end if
    call MPI_BCAST(klow        ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(khigh       ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(dtav        ,1,MY_REAL   ,0,comm3d,ierr)
    call MPI_BCAST(lfielddump  ,1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(ldiracc     ,1,MPI_LOGICAL,0,comm3d,ierr)

    tnext   = dtav-1e-3+btime
    if(.not.(lfielddump)) return
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if
    if (lnetcdf) then
      dtav = dtav_glob
      fname(11:13) = cmyid
      fname(15:17) = cexpnr
      call ncinfo(tncname(1,:),'time','Time','s','time')
      call ncinfo(ncname( 1,:),'u','West-East velocity','m/s','mttt')
      call ncinfo(ncname( 2,:),'v','South-North velocity','m/s','tmtt')
      call ncinfo(ncname( 3,:),'w','Vertical velocity','m/s','ttmt')
      call ncinfo(ncname( 4,:),'qt','Total water mixing ratio','1e-5kg/kg','tttt')
      call ncinfo(ncname( 5,:),'ql','Liquid water mixing ratio','1e-5kg/kg','tttt')
      call ncinfo(ncname( 6,:),'thl','Liquid water potential temperature above 300K','K','tttt')

      call open_nc(fname,  ncid,n1=imax,n2=jmax,n3=khigh-klow+1)
      call define_nc( ncid, 1, tncname)
      call writestat_dims_nc(ncid)
      call redefine_nc(ncid)
      call define_nc( ncid, NVar, ncname)
    end if

  end subroutine initfielddump

!> Do fielddump. Collect data to truncated (2 byte) integers, and write them to file
  subroutine fielddump
    use modfields, only : um,vm,wm,thlm,qtm,ql0
    use modsurfdata,only : thls,qts,thvs
    use modglobal, only : imax,i1,ih,jmax,j1,jh,kmax,k1,rk3step,&
                          timee,dt_lim,cexpnr,ifoutput
    use modmpi,    only : myid,cmyid
    use modstat_nc, only : lnetcdf, writestat_nc
    implicit none

    integer(KIND=selected_int_kind(4)), allocatable :: field(:,:,:),vars(:,:,:,:)
    integer i,j,k
    integer :: writecounter = 1
    integer :: reclength


    if (.not. lfielddump) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+dtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    allocate(field(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(vars(imax,jmax,khigh-klow+1,nvar))

    reclength = imax*jmax*(khigh-klow+1)*2

    field = NINT(1.0E3*um,2)
    if (lnetcdf) vars(:,:,:,1) = field(2:i1,2:j1,klow:khigh)
    if (ldiracc) then
      open (ifoutput,file='wbuu.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
      write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
    else
      open  (ifoutput,file='wbuu.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
      write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
    end if
    close (ifoutput)

    field = NINT(1.0E3*vm,2)
    if (lnetcdf) vars(:,:,:,2) = field(2:i1,2:j1,klow:khigh)
    if (ldiracc) then
      open (ifoutput,file='wbvv.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
      write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
    else
      open  (ifoutput,file='wbvv.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
      write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
    end if
    close (ifoutput)

    field = NINT(1.0E3*wm,2)
    if (lnetcdf) vars(:,:,:,3) = field(2:i1,2:j1,klow:khigh)
    if (ldiracc) then
      open (ifoutput,file='wbww.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
      write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
    else
      open  (ifoutput,file='wbww.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
      write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
    end if
    close (ifoutput)

    field = NINT(1.0E5*qtm,2)
    if (lnetcdf) vars(:,:,:,4) = field(2:i1,2:j1,klow:khigh)
    if (ldiracc) then
      open (ifoutput,file='wbqt.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
      write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
    else
      open  (ifoutput,file='wbqt.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
      write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
    end if
    close (ifoutput)

    field = NINT(1.0E5*ql0,2)
    if (lnetcdf) vars(:,:,:,5) = field(2:i1,2:j1,klow:khigh)
    if (ldiracc) then
      open (ifoutput,file='wbql.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
      write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
    else
      open  (ifoutput,file='wbql.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
      write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
    end if
    close (ifoutput)

    field = NINT(1.0E3*(thlm-300),2)
    if (lnetcdf) vars(:,:,:,6) = field(2:i1,2:j1,klow:khigh)
    if (ldiracc) then
      open (ifoutput,file='wbtl.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
      write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
    else
      open  (ifoutput,file='wbtl.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
      write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
    end if
    close (ifoutput)
    if(lnetcdf) then
      call writestat_nc(ncid,1,tncname,(/timee/),nrec,.true.)
      call writestat_nc(ncid,nvar,ncname,vars,nrec,imax,jmax,khigh-klow+1)
    end if

!     if (myid==0) then
!       open(ifoutput, file='wbthls.'//cexpnr,form='formatted',position='append')
!       write(ifoutput,'(F12.1 3F12.5)') timee,thls, qts,thvs
!       close(ifoutput)
!     end if
    writecounter=writecounter+1

    deallocate(field,vars)

  end subroutine fielddump
!> Clean up when leaving the run
  subroutine exitfielddump
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(lfielddump .and. lnetcdf) call exitstat_nc(ncid)
  end subroutine exitfielddump

end module modfielddump
