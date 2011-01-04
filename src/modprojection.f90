!----------------------------------------------------------------------------
! This file is part of DALES.
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
! Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!----------------------------------------------------------------------------
!
!
module modprojection

    !-----------------------------------------------------------------|
    !                                                                 |
    !*** *projection*  dumps an instantenous integrated projection    |
    !                                                   of the field  |
    !                                                                 |
    !     purpose.                                                    |
    !     --------                                                    |
    !                                                                 |
    !                                                                 |
    !    projections in the yz-plane and in the xy-plane            |
    !        of thl,thv,qt,ql                                   |
    !____________________SETTINGS_AND_SWITCHES________________________|
    !                     IN &NAMprojection                         |
    !                                                                 |
    !    dtav        SAMPLING INTERVAL                                |
    !                                                                 |
    !    lprojection SWITCH TO ENABLE projection                    |
    !                                                                 |
    !    projectheight HEIGHT OF THE XY-projection                    |
    !                                                                 |
    !    projectplane  LOCATION OF THE YZ-PLANE ON EVERY PROCESSOR      |
    !-----------------------------------------------------------------|
  use modglobal, only : longint

implicit none
private
PUBLIC :: initprojection, projection
save
!NetCDF variables
  integer,parameter :: nvar = 8
  integer :: ncid = 0
  integer :: nrec = 0
  character(80) :: fname = 'proj.xxx.xxx.nc'
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav
  integer(kind=longint)  :: idtav,tnext
  logical:: lproject = .false. ! switch for conditional sampling cloud (on/off)
  real    :: projectheight = 0.
  integer :: ksplit = 2  !lowest integration boundary

contains

  subroutine initprojection
    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_logical,mpi_integer,cmyid
    use modglobal,only :imax,jmax,ifnamopt,fname_options,dtmax,rk3step, dtav_glob,ladaptive,j1,kmax,dt_lim,tres,btime,cexpnr,zf
    use modstat_nc, only : lnetcdf, open_nc,define_nc,ncinfo, writestat_dims_nc
    implicit none

    integer :: ierr

    namelist/NAMprojection/ &
    lproject, dtav, projectheight

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMprojection,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMprojection'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMprojection'
      endif
      write(6 ,NAMprojection)
      close(ifnamopt)
    end if
    do ksplit=1,kmax
      if (zf(ksplit)>=projectheight) exit
    end do
    ksplit = ksplit - 1
    call MPI_BCAST(dtav       ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(lproject     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(ksplit,1,MPI_INTEGER,0,comm3d,mpierr)
    idtav = dtav/tres

    tnext      = idtav   +btime
    if(.not.(lproject)) return
    dt_lim = min(dt_lim,tnext)

    if(ksplit>kmax) then
      stop 'projection: projection out of range'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'projection: dtav should be a integer multiple of dtmax'
    end if
    fname(6:8) = cmyid
    fname(10:12) = cexpnr
    call ncinfo(tncname(1,:),'time','Time','s','time')
    call ncinfo(ncname( 1,:),'thlxylow','Subcloud Integrated liquid water potential temperature path','K/m','tt0t')
    call ncinfo(ncname( 2,:),'thlxyhigh','Cloudlayer Integrated liquid water potential temperature path','K/m','tt0t')
    call ncinfo(ncname( 3,:),'thvxylow','Subcloud Integrated virtual potential temperature path','K/m','tt0t')
    call ncinfo(ncname( 4,:),'thvxyhigh','Cloudlayer Integrated virtual potential temperature path','K/m','tt0t')
    call ncinfo(ncname( 5,:),'qtxylow','Subcloud Integrated total water path','kg/m','tt0t')
    call ncinfo(ncname( 6,:),'qtxyhigh','Cloudlayer Integrated total water path','kg/m','tt0t')
    call ncinfo(ncname( 7,:),'qlxylow','Subcloud Integrate liquid waterpath ','kg/m','tt0t')
    call ncinfo(ncname( 8,:),'qlxyhigh','Cloudlayer Integrate liquid waterpath ','kg/m','tt0t')
    call open_nc(fname,  ncid,nrec,n1=imax,n2=jmax)
    if (nrec==0) then
      call define_nc( ncid, 1, tncname)
      call writestat_dims_nc(ncid)
    end if
    call define_nc( ncid, NVar, ncname)


  end subroutine initprojection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine projection
    use modglobal, only : rk3step,timee,dt_lim
    implicit none


    if (.not. lproject) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))
     call wrthorz


  end subroutine projection
  subroutine wrthorz
    use modglobal, only : i1,j1,imax,jmax,kmax,nsv,cexpnr,ifoutput,rtimee ,dzf
    use modfields, only : thlm,qtm,svm,thv0h,ql0,thl0av,qt0av,rhof
    use modmpi,    only : cmyid
    use modstat_nc,  only : lnetcdf, writestat_nc
   implicit none


    ! LOCAL
    integer i,j
!     character(20) :: name
    real,allocatable,dimension(:,:,:) :: vars
!     
!     open(ifoutput,file='projh_thl.'//cmyid//'.'//cexpnr,position='append',action='write')
!     write(ifoutput,'(es12.5)') ((sum(thlm(i,j,projectheight:kmax)-thl0av(projectheight:kmax)),i=2,i1),j=2,j1)
!     close(ifoutput)
! 
!     open(ifoutput,file='projh_thv.'//cmyid//'.'//cexpnr,position='append',action='write')
!     write(ifoutput,'(es12.5)') ((sum(thv0h(i,j,projectheight:kmax)-thl0av(projectheight:kmax)),i=2,i1),j=2,j1)
!     close(ifoutput)
! 
!     open(ifoutput,file='projh_qt.'//cmyid//'.'//cexpnr,position='append',action='write')
!     write(ifoutput,'(es12.5)') ((1.e3*sum(qtm(i,j,projectheight:kmax)-qt0av(projectheight:kmax)),i=2,i1),j=2,j1)
!     close(ifoutput)
! 
!     open(ifoutput,file='projh_ql.'//cmyid//'.'//cexpnr,position='append',action='write')
!     write(ifoutput,'(es12.5)') ((1.e3*sum(ql0(i,j,projectheight:kmax)),i=2,i1),j=2,j1)
!     close(ifoutput)
! 
!     do n = 1,nsv
!       name = 'projh_snn.'//cmyid//'.'//cexpnr
!       write(name(7:8),'(i2.2)') n
!       open(ifoutput,file=name,position='append',action='write')
!       write(ifoutput,'(es12.5)') ((sum(svm(i,j,projectheight:kmax,n)),i=2,i1),j=2,j1)
!       close(ifoutput)
!     end do

    if (lnetcdf) then
      allocate(vars(1:imax,1:jmax,nvar))
      do j = 2,j1
        do i = 2,i1
          vars(i-1,j-1,1) = sum(thlm (i,j,1:ksplit)*rhof(1:ksplit)*dzf(1:ksplit))
          vars(i-1,j-1,2) = sum(thlm (i,j,ksplit+1:kmax)*rhof(ksplit+1:kmax)*dzf(ksplit+1:kmax))
          vars(i-1,j-1,3) = sum(thv0h(i,j,1:ksplit)*rhof(1:ksplit)*dzf(1:ksplit))
          vars(i-1,j-1,4) = sum(thv0h(i,j,ksplit+1:kmax)*rhof(ksplit+1:kmax)*dzf(ksplit+1:kmax))
          vars(i-1,j-1,5) = sum(qtm  (i,j,1:ksplit)*rhof(1:ksplit)*dzf(1:ksplit))
          vars(i-1,j-1,6) = sum(qtm  (i,j,ksplit+1:kmax)*rhof(ksplit+1:kmax)*dzf(ksplit+1:kmax))
          vars(i-1,j-1,7) = sum(ql0  (i,j,1:ksplit)*rhof(1:ksplit)*dzf(1:ksplit))
          vars(i-1,j-1,8) = sum(ql0  (i,j,ksplit+1:kmax)*rhof(ksplit+1:kmax)*dzf(ksplit+1:kmax))
        end do
      end do
      
      
      call writestat_nc(ncid,1,tncname,(/rtimee/),nrec,.true.)
      call writestat_nc(ncid,nvar,ncname(1:nvar,:),vars,nrec,imax,jmax)
      deallocate(vars)
    end if
    

  end subroutine wrthorz
end module modprojection
