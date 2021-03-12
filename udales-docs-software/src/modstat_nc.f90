!> \file modstat_nc.f90
!!  Background routines to write NetCDF output

!>
!!  Background routines to write NetCDF output.
!>
!! All calls to the netcdf library should be directed through here.
!! Inspired on the UCLA-LES routine by Bjorn Stevens.
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  \todo documentation
!!   \todo restartfiles in NetCDF?
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
module modstat_nc
    use netcdf
    use modmpi, only : myid
    implicit none
    integer, save :: timeID=0, ztID=0, zmID=0, xtID=0, xmID=0, ytID=0, ymID=0, ztsID=0, fctID=0, lyrID=0
    real(kind=4) :: nc_fillvalue = -999.
!> The only interface necessary to write data to netcdf, regardless of the dimensions.
    interface writestat_nc
      module procedure writestat_time_nc
      module procedure writestat_1D_nc
      module procedure writestat_2D_nc
      module procedure writestat_3D_nc
      module procedure writestat_3D_short_nc
    end interface writestat_nc
contains


  subroutine initstat_nc
    use modglobal, only : kmax,ifnamopt,fname_options,iexpnr
    use modmpi,    only : mpierr,mpi_logical,comm3d,myid
    implicit none

    integer             :: ierr

  end subroutine initstat_nc
!
! ----------------------------------------------------------------------
!> Subroutine Open_NC: Opens a NetCDF File and identifies starting record
!
  subroutine open_nc (fname, ncid,nrec,n1, n2, n3, ns, nfcts, nlyrs)
    use modglobal, only : author,version,timee
    implicit none
    integer, intent (out) :: ncid,nrec
    integer, optional, intent (in) :: n1, n2, n3, ns, nfcts, nlyrs
    character (len=40), intent (in) :: fname

    character (len=12):: date='',time=''
    integer :: iret,varid,ncall,RecordDimID
    real, allocatable :: xtimes(:)
    logical :: exans

    inquire(file=trim(fname),exist=exans)
    !write(*,*) 'opennc'
    ncall = 0
    if (.not.exans) then

      call date_and_time(date,time)
      iret = nf90_create(fname,NF90_SHARE,ncid)
      iret = nf90_put_att(ncid,NF90_GLOBAL,'title',fname)
      iret = nf90_put_att(ncid,NF90_GLOBAL,'history','Created on '//trim(date)//' at '//trim(time))
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'Source',trim(version))
      iret = nf90_put_att(ncid, NF90_GLOBAL, 'Author',trim(author))
      iret = nf90_def_dim(ncID, 'time', NF90_UNLIMITED, timeID)
      if (present(n1)) then
        iret = nf90_def_dim(ncID, 'xt', n1, xtID)
        iret = nf90_def_dim(ncID, 'xm', n1, xmID)
        iret = nf90_def_var(ncID,'xt',NF90_FLOAT,(/xtID/) ,VarID)
        iret=nf90_put_att(ncID,VarID,'longname','West-East displacement of cell centers')
        iret=nf90_put_att(ncID,VarID,'units','m')
        iret = nf90_def_var(ncID,'xm',NF90_FLOAT,(/xmID/),VarID)
        iret=nf90_put_att(ncID,VarID,'longname','West-East displacement of cell edges')
        iret=nf90_put_att(ncID,VarID,'units','m')
      end if
      if (present(n2)) then
        iret = nf90_def_dim(ncID, 'yt', n2, ytID)
        iret = nf90_def_dim(ncID, 'ym', n2, ymID)
        iret = nf90_def_var(ncID,'yt',NF90_FLOAT,ytID ,VarID)
        iret=nf90_put_att(ncID,VarID,'longname','South-North displacement of cell centers')
        iret=nf90_put_att(ncID,VarID,'units','m')
        iret = nf90_def_var(ncID,'ym',NF90_FLOAT,ymID,VarID)
        iret=nf90_put_att(ncID,VarID,'longname','South-North displacement of cell edges')
        iret=nf90_put_att(ncID,VarID,'units','m')
      end if
      if (present(n3)) then
        iret = nf90_def_dim(ncID, 'zt', n3, ztID)
        iret = nf90_def_dim(ncID, 'zm', n3, zmID)
        iret = nf90_def_var(ncID,'zt',NF90_FLOAT,(/ztID/) ,VarID)
        iret=nf90_put_att(ncID,VarID,'longname','Vertical displacement of cell centers')
        iret=nf90_put_att(ncID,VarID,'units','m')
        iret = nf90_def_var(ncID,'zm',NF90_FLOAT,(/zmID/),VarID)
        iret=nf90_put_att(ncID,VarID,'longname','Vertical displacement of cell edges')
        iret=nf90_put_att(ncID,VarID,'units','m')
      end if
      if (present(ns)) then
        iret = nf90_def_dim(ncID, 'zts', ns, ztsID)
        iret = nf90_def_var(ncID,'zts',NF90_FLOAT,(/ztsID/) ,VarID)
        iret=nf90_put_att(ncID,VarID,'longname','Soil level depth of cell centers')
        iret=nf90_put_att(ncID,VarID,'units','m')
      end if
      if (present(nfcts)) then
        iret = nf90_def_dim(ncID, 'fct', nfcts, fctID)
        iret = nf90_def_var(ncID, 'fct',NF90_INT,(/fctID/) ,VarID)
        iret=nf90_put_att(ncID,VarID,'longname','Facet number')
      end if
      if (present(nlyrs)) then
        iret = nf90_def_dim(ncID, 'lyr', nlyrs, lyrID)
        iret = nf90_def_var(ncID, 'lyr',NF90_INT,(/lyrID/) ,VarID)
        iret=nf90_put_att(ncID,VarID,'longname','Number of wall layers')
      end if

    else
       nrec = 0
       ncall= 0
       iret = nf90_open (trim(fname), NF90_WRITE, ncid)
       iret = nf90_inquire(ncid, unlimitedDimId = RecordDimID)
       iret = nf90_inquire_dimension(ncid, RecordDimID, len=nrec)
       if (nrec>0) then
        iret = nf90_inq_varid(ncid,'time',timeID)
        allocate (xtimes(nrec))
        iret = nf90_get_var(ncid, timeId, xtimes(1:nrec))

        do while(xtimes(ncall+1) < timee - spacing(1.))
            ncall=ncall+1
            if (ncall >= nrec) exit 
        end do
        deallocate(xtimes)
       end if
       if (present(n1)) then
         iret = nf90_inq_dimid(ncid,'xt',xtId)
         iret = nf90_inq_dimid(ncid,'xm',xmId)
       end if
       if (present(n2)) then
         iret = nf90_inq_dimid(ncid,'yt',ytId)
         iret = nf90_inq_dimid(ncid,'ym',ymId)
       end if
       if (present(n3)) then
         iret = nf90_inq_dimid(ncid,'zt',ztId)
         iret = nf90_inq_dimid(ncid,'zm',zmId)
       end if
       if (present(ns)) then
         iret = nf90_inq_dimid(ncid,'zts',ztsId)
       end if
       if (present(nfcts)) then
         iret = nf90_inq_dimid(ncid,'fct',fctId)
       end if
    end if
    nrec = ncall
    iret = nf90_sync(ncid)

    iret= nf90_enddef(ncID)

  end subroutine open_nc

  !
  ! ----------------------------------------------------------------------
  !> Subroutine Define_NC: Defines the structure of the nc file (if not
  !! already open)
  !

  subroutine define_nc(ncID, nVar, sx)
    implicit none
    integer, intent (in) :: nVar, ncID
    character (*), intent (in) :: sx(nVar,4)

    integer, save ::  dim_mttt(4) = 0, dim_tmtt(4) = 0, dim_ttmt(4) = 0, dim_tttt(4) = 0, &
                      dim_tt(2)= 0, dim_mt(2)= 0,dim_t0tt(3)=0,dim_m0tt(3)=0,dim_t0mt(3)=0,&
                      dim_m0mt(3)=0, dim_tt0t(3)=0, &
                      dim_mt0t(3)=0,dim_tm0t(3)=0,dim_0ttt(3)=0,dim_0mtt(3)=0,dim_0tmt(3)=0,&
                      dim_tts(2)=0,dim_t0tts(3)=0,dim_0ttts(3)=0,dim_tttts(4)=0,dim_ttt0(3)=0,& !tg3315 added last one
                      dim_mtmt(4),dim_tmmt(4),dim_mmtt(4),& !bss116
                      dim_ft(2), dim_flt(3)!SO

    integer :: iret, n, VarID
    !write(*,*) 'definenc'
    iret = nf90_inq_dimid(ncid,'time',timeId)
    iret = nf90_inq_dimid(ncid,'xt',xtId)
    iret = nf90_inq_dimid(ncid,'xm',xmId)
    iret = nf90_inq_dimid(ncid,'yt',ytId)
    iret = nf90_inq_dimid(ncid,'ym',ymId)
    iret = nf90_inq_dimid(ncid,'zt',ztId)
    iret = nf90_inq_dimid(ncid,'zm',zmId)
    iret = nf90_inq_dimid(ncid,'zts',ztsId)
    iret = nf90_inq_dimid(ncid,'fct',fctId) ! so4718 for energy balance output
    iret = nf90_inq_dimid(ncid,'lyr',lyrId) ! so4718 for energy balance output
    
    iret = nf90_redef(ncid) 
    dim_tt = (/ztId,timeId/)
    dim_mt = (/zmId,timeId/)

    dim_t0tt= (/xtID,ztID,timeId/)! thermo point
    dim_t0mt= (/xtID,zmID,timeId/)! zpoint
    dim_m0tt= (/xmID,ztID,timeId/)! upoint
    dim_m0mt= (/xmID,ztID,timeId/)! uw stats point
    dim_tt0t= (/xtID,ytID,timeId/)! thermo point
    dim_tm0t= (/xtID,ymID,timeId/)! vpoint
    dim_mt0t= (/xmID,ytID,timeId/)! upoint
    dim_0ttt= (/ytID,ztID,timeId/)! thermo point
    dim_0tmt= (/ytID,zmID,timeId/)! wpoint
    dim_0mtt= (/ymID,ztID,timeId/)! vpoint

    dim_tttt= (/xtID,ytID,ztID,timeId/)! thermo point
    dim_ttmt= (/xtID,ytID,zmID,timeId/)! zpoint
    dim_mttt= (/xmID,ytID,ztID,timeId/)! upoint
    dim_tmtt= (/xtID,ymID,ztId,timeId/)! ypoint
    dim_mtmt= (/xmID,ytID,zmId,timeId/)! uw stats point bss116
    dim_tmmt= (/xtID,ymID,zmId,timeId/)! vw stats point bss116
    dim_mmtt= (/xmID,ymID,ztId,timeId/)! uv stats point bss116

    dim_ttt0= (/xtID,ytID,ztID/)! stats point tg3315

    dim_tts = (/ztsId,timeId/)
    dim_t0tts= (/xtID,ztsID,timeId/)! thermo soil point
    dim_0ttts= (/ytID,ztsID,timeId/)! thermo point
    dim_tttts= (/xtID,ytID,ztsID,timeId/)! thermo point
    
    dim_ft = (/fctID,timeId/)
    dim_flt = (/fctID,lyrID,timeId/)
    
    do n=1,nVar
!      write(*,*) 'n', n
!      write(*,*) "dummyline1"
!      write(*,*) 'sx1', sx(1,:)
!      write(*,*) 'sx2', sx(2,:)
!      write(*,*) "dummyline2"
!      write(*,*) 'trim(sx(n,1))', trim(sx(n,1))
!      write(*,*) 'trim(sx(n,2))', trim(sx(n,2))
!      write(*,*) 'trim(sx(n,3))', trim(sx(n,3))
!      write(*,*) 'trim(sx(n,4))', trim(sx(n,4))
!      write (*,*) 'ncID', ncID
      iret = nf90_inq_varid(ncid, trim(sx(n,1)), VarID)
      if (iret == 0) cycle
      select case(trim(sx(n,4)))
        case ('time')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,(/timeID/) ,VarID)
        case ('tt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_tt ,VarID)
        case ('mt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_mt,VarID)
  !2D Fields
        case ('t0tt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_t0tt,VarID)
        case ('t0mt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_t0mt,VarID)
        case ('m0tt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_m0tt,VarID)
        case ('m0mt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_m0mt,VarID)
        case ('tt0t')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_tt0t,VarID)
        case ('tm0t')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_tm0t,VarID)
        case ('mt0t')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_mt0t,VarID)
        case ('0ttt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_0ttt,VarID)
        case ('0tmt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_0tmt,VarID)
        case ('0mtt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_0mtt,VarID)
        case ('ttt0')                  !tg3315 for uav,vav,wav etc.
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_ttt0,VarID)
  !3D Fields
        case ('tttt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_tttt,VarID)
        case ('mttt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_mttt,VarID)
        case ('tmtt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_tmtt,VarID)
        case ('ttmt')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_ttmt,VarID)
        case ('mtmt')                                                   ! bss116
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_mtmt,VarID)
        case ('tmmt')                                                   ! bss116
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_tmmt,VarID)
        case ('mmtt')                                                   ! bss116
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_mmtt,VarID)
!Soil fields
        case ('tts')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_tts ,VarID)
        case ('t0tts')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_t0tts,VarID)
        case ('0ttts')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_0ttts,VarID)
        case ('tttts')
          iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_tttts,VarID)
       
!Facet information
		case ('ft')
		  iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_ft,VarID)
		case ('flt')
		  iret=nf90_def_var(ncID,sx(n,1),NF90_FLOAT,dim_flt,VarID)

        case default
        write(0, *) 'nvar', nvar, sx(n,:)
        write(0, *) 'ERROR: Bad dimensional information ',sx(n,:)
        stop 1
        ! call appl_abort(0)
      end select
      if (iret/=0) then
!        write (*,*) 'nvar', nvar, sx(n,:)
!        write (*,*) 'ncID', ncID
        call nchandle_error(iret)
      end if
      iret = nf90_put_att(ncID,VarID,'longname',sx(n,2))
      iret = nf90_put_att(ncID,VarID,'units',sx(n,3))
      iret = nf90_put_att(ncid, VarID, '_FillValue',nc_fillvalue)

    end do
    iret= nf90_enddef(ncID)
  end subroutine define_nc

  subroutine redefine_nc(ncid)
  implicit none
    integer, intent(in) :: ncid
    integer :: iret
    iret = nf90_redef(ncid)
  end subroutine redefine_nc

  subroutine exitstat_nc(ncid)

   implicit none
   integer, intent(in) :: ncid
   integer status

   status = nf90_close(ncid)
   if (status /= nf90_noerr) call nchandle_error(status)
 end subroutine exitstat_nc
  subroutine writestat_dims_nc(ncid)
    use modglobal, only : xf,xh,dy,zf,zh,jmax
    use modmpi, only : myid
    implicit none
    integer, intent(in) :: ncid
    integer             :: i=0,iret,length,varid
    write(*,*) 'writestat_dims_nc'
    iret = nf90_inq_varid(ncid, 'xt', VarID)
    if (iret==0) iret=nf90_inquire_dimension(ncid, xtID, len=length)
    !if (iret==0) iret = nf90_put_var(ncid, varID, zf(1:length),(/1/))
    if (iret==0) iret = nf90_put_var(ncid, varID, xf(1:length),(/1/))
    
    iret = nf90_inq_varid(ncid, 'xm', VarID)
    if (iret==0) iret=nf90_inquire_dimension(ncid, xmID, len=length)
    if (iret==0) iret = nf90_put_var(ncid, varID, xh(1:length),(/1/))

    iret = nf90_inq_varid(ncid, 'yt', VarID)
    if (iret==0) iret=nf90_inquire_dimension(ncid, ytID, len=length)
    if (iret==0) iret = nf90_put_var(ncid, varID, (/(dy*(0.5+i)+myid*jmax*dy,i=0,length-1)/),(/1/))
    iret = nf90_inq_varid(ncid, 'ym', VarID)
    if (iret==0) iret=nf90_inquire_dimension(ncid, ymID, len=length)
    if (iret==0) iret = nf90_put_var(ncid, varID, (/(dy*i+myid*jmax*dy,i=0,length-1)/),(/1/))

    iret = nf90_inq_varid(ncid, 'zt', VarID)
    if (iret==0) iret=nf90_inquire_dimension(ncid,ztID, len=length)
    if (iret==0) iret = nf90_put_var(ncid, varID, zf(0:length-1),(/1/))  !ils13, 29.06.2017 zf starts at 0, not at 1
    iret = nf90_inq_varid(ncid, 'zm', VarID)
    if (iret==0) iret=nf90_inquire_dimension(ncid, zmID, len=length)
    if (iret==0) iret = nf90_put_var(ncid, varID, zh(0:length-1),(/1/))  !same for zh
    !if (isurf==1) then
      !iret = nf90_inq_varid(ncid, 'zts', VarID)
      !if (iret==0) iret = nf90_inquire_dimension(ncid, ztsID, len=length)
      !if (iret==0) iret = nf90_put_var(ncid, varID, zsoilc(1:length),(/1/))
    !end if

  end subroutine writestat_dims_nc

  subroutine writestat_time_nc(ncid,nvar,ncname,vars,nrec,lraise)
    implicit none
    integer, intent(in)                      :: ncid,nvar
    integer, intent(inout)                   :: nrec
    real,dimension(nvar),intent(in)          :: vars
    character(*), dimension(:,:),intent(in)  :: ncname
    logical, intent(in)                      :: lraise

    integer :: iret,n,varid
    !write(*,*) 'time-nc'
    if(lraise) then
      nrec = nrec+1
    end if
    do n=1,nvar
       iret = nf90_inq_varid(ncid, ncname(n,1), VarID)
       iret = nf90_put_var(ncid, VarID, vars(n), start=(/nrec/))
    end do
    iret = nf90_sync(ncid)

  end subroutine writestat_time_nc

  subroutine writestat_1D_nc(ncid,nvar,ncname,vars,nrec,dim1)
    implicit none
    integer, intent(in)                      :: ncid,nvar,dim1
    integer, intent(in)                      :: nrec
    real,dimension(dim1,nvar),intent(in)     :: vars
    character(*), dimension(:,:),intent(in)  :: ncname

    integer :: iret,n,varid
    do n=1,nvar
      iret = nf90_inq_varid(ncid, ncname(n,1), VarID)
      iret = nf90_put_var(ncid, VarID, vars(1:dim1,n),(/1,nrec/),(/dim1,1/))
    end do
    iret = nf90_sync(ncid)

  end subroutine writestat_1D_nc

  subroutine writestat_2D_nc(ncid,nvar,ncname,vars,nrec,dim1,dim2)
    implicit none
    integer, intent(in)                      :: ncid,nvar,dim1,dim2
    integer, intent(in)                      :: nrec
    real,dimension(:,:,:),intent(in)         :: vars
    character(*), dimension(:,:),intent(in)  :: ncname

    integer :: iret,n,varid
    do n=1,nvar
      iret = nf90_inq_varid(ncid, ncname(n,1), VarID)
      iret = nf90_put_var(ncid, VarID, vars(1:dim1,1:dim2,n),(/1,1,nrec/),(/dim1,dim2,1/))
    end do
    iret = nf90_sync(ncid)

  end subroutine writestat_2D_nc
  subroutine writestat_3D_nc(ncid,nvar,ncname,vars,nrec,dim1,dim2,dim3)
    implicit none
    integer, intent(in)                      :: ncid,nvar,dim1,dim2,dim3
    integer, intent(in)                      :: nrec
    real,dimension(dim1,dim2,dim3,nvar),intent(in)       :: vars
    character(*), dimension(:,:),intent(in)  :: ncname

    integer :: iret,n,varid
    !write(*,*) 'write 3Dnc'
    do n=1,nvar
      iret = nf90_inq_varid(ncid, ncname(n,1), VarID)
           !write(*,*) 'MYID,.', myid
           !write(*,*) "nth netcdf variable",n
           !write(*,*) 'ncid',ncid
           !write(*,*) "n",n
           !write(*,*) 'ncname(n,1)',ncname(n,1)
           !write(*,*) 'VarID',VarID
           !write(*,*) 'nrec',nrec
           !write(*,*) 'imax',dim1
           !write(*,*) 'jmax',dim2
           !write(*,*) 'kmax',dim3
           !write(*,*) 'shape(vars)',shape(vars)
      iret = nf90_put_var(ncid, VarID, vars(1:dim1,1:dim2,1:dim3,n),(/1,1,1,nrec/),(/dim1,dim2,dim3,1/))
    end do
    iret = nf90_sync(ncid)

  end subroutine writestat_3D_nc
  subroutine writestat_3D_short_nc(ncid,nvar,ncname,vars,nrec,dim1,dim2,dim3)
    implicit none
    integer, intent(in)                      :: ncid,nvar,dim1,dim2,dim3
    integer, intent(in)                      :: nrec
    integer(KIND=selected_int_kind(4)),dimension(dim1,dim2,dim3,nvar),intent(in)       :: vars
    character(*), dimension(:,:),intent(in)  :: ncname

    integer :: iret,n,varid
    !write(*,*) '3Dnc short'
    do n=1,nvar
      iret = nf90_inq_varid(ncid, ncname(n,1), VarID)
      iret = nf90_put_var(ncid, VarID, vars(1:dim1,1:dim2,1:dim3,n),(/1,1,1,nrec/),(/dim1,dim2,dim3,1/))
    end do
    iret = nf90_sync(ncid)

  end subroutine writestat_3D_short_nc


  subroutine ncinfo(out,in1,in2,in3,in4)

    implicit none
    character(*), dimension(4),intent(out) ::out
    character(*), intent(in) ::in1,in2,in3,in4
    out(1) = in1
    out(2) = in2
    out(3) = in3
    out(4) = in4
  end subroutine ncinfo

  subroutine nchandle_error(status)
    use netcdf
    implicit none

    integer, intent(in) :: status

    if(status /= nf90_noerr) then
      write(0, *) trim(nf90_strerror(status))
      write(0, *) 'ERROR: status',  status
      stop 1
    end if

  end subroutine nchandle_error

end module modstat_nc
