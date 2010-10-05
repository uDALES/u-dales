!> \file modlsmcrosssection.f90
!!   Dumps an instantenous lsmcrosssection of the field

!>
!! Dumps an instantenous lsmcrosssection of the field.
!>
!! lsmcrosssections in the yz-plane and in the xy-plane            |
    !        of u,v,w,thl,thv,qt,ql. Written to movv_*.expnr and movh_*.expnr
!! If netcdf is true, this module leads the cross.myid.expnr.nc output

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
module modlsmcrosssection


  use modglobal, only : longint
  use modsurfdata,only : ksoilmax

implicit none
private
PUBLIC :: initlsmcrosssection, lsmcrosssection,exitlsmcrosssection
save
!NetCDF variables
  integer,parameter :: nvar = 2,nvar3=12
  integer :: ncid1 = 0
  integer :: ncid2 = 0
  integer :: ncid3 = 0
  integer :: nrec1 = 0
  integer :: nrec2 = 0
  integer :: nrec3 = 0
  integer :: crossheight
!   integer :: nxy = 0
!   integer :: cross
!   integer :: nrc
  character(4) :: cheight
  character(80) :: fname1 = 'lsmcrossxz.xxx.xxx.nc'
  character(80) :: fname2 = 'lsmcrossxy.xxxx.xxx.xxx.nc'
  character(80) :: fname3 = 'surfcross.xxx.xxx.nc'
  character(80),dimension(nvar,4) :: ncname1
  character(80),dimension(1,4) :: tncname1
  character(80),dimension(nvar,4) :: ncname2
  character(80),dimension(1,4) :: tncname2
  character(80),dimension(nvar3,4) :: ncname3
  character(80),dimension(1,4) :: tncname3


  real    :: dtav
  integer(kind=longint) :: idtav,tnext
  logical :: lcross = .false. !< switch for doing the lsmcrosssection (on/off)
  integer :: crossplane = 2 !< Location of the xz lsmcrosssection

contains
!> Initializing lsmcrosssection. Read out the namelist, initializing the variables
  subroutine initlsmcrosssection
    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_logical,mpi_integer,cmyid
    use modglobal,only :imax,jmax,ifnamopt,fname_options,dtmax,rk3step, dtav_glob,ladaptive,j1,i1,dt_lim,cexpnr,tres,btime
    use modstat_nc,only : lnetcdf,open_nc, define_nc, redefine_nc,ncinfo,writestat_dims_nc
   implicit none

    integer :: ierr

    namelist/NAMLSMCROSSSECTION/ &
    lcross, dtav, crossheight, crossplane

    crossheight=2
    ncid2=2
    nrec2=0

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMLSMCROSSSECTION,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMLSMCROSSSECTION'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMLSMCROSSSECTION'
      endif
      write(6 ,NAMLSMCROSSSECTION)
      close(ifnamopt)
    end if

    call MPI_BCAST(dtav       ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(lcross     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(crossheight,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(crossplane ,1,MPI_INTEGER,0,comm3d,mpierr)


    idtav = dtav/tres
    tnext   = idtav+btime
    if(.not.(lcross)) return
    dt_lim = min(dt_lim,tnext)

    if((crossheight>ksoilmax) .or. crossplane>j1) then
      stop 'lsmcrosssection: lsmcrosssection out of range'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'lsmcrosssection: dtav should be a integer multiple of dtmax'
    end if
    if (lnetcdf) then
      if (myid==0) then
        fname1(12:14) = cmyid
        fname1(16:18) = cexpnr
        call ncinfo(tncname1(1,:),'time','Time','s','time')
        call ncinfo(ncname1( 1,:),'tsoil', 'xz crosssection of the Soil temperature','K','t0tts')
        call ncinfo(ncname1( 2,:),'phiw', 'xz crosssection of the Soil moisture','m3/m3','t0tts')
        call open_nc(fname1,  ncid1,nrec1,n1=imax,ns=ksoilmax)
        if (nrec1==0) then
          call define_nc( ncid1, 1, tncname1)
          call writestat_dims_nc(ncid1)
          call redefine_nc(ncid1)
          call define_nc( ncid1, NVar, ncname1)
        end if
      end if
        write(cheight,'(i4.4)') crossheight
        fname2(12:15) = cheight
        fname2(17:19) = cmyid
        fname2(21:23) = cexpnr
        call ncinfo(tncname2(1,:),'time','Time','s','time')
        call ncinfo(ncname2( 1,:),'tsoil', 'xy crosssection of the Soil temperature','K','tt0t')
        call ncinfo(ncname2( 2,:),'phiw', 'xy crosssection of the Soil moisture','m3/m3','tt0t')
        call open_nc(fname2,  ncid2,nrec2,n1=imax,n2=jmax)
        if (nrec2==0) then
          call define_nc( ncid2, 1, tncname2)
          call writestat_dims_nc(ncid2)
          call redefine_nc(ncid2)
          call define_nc( ncid2, NVar, ncname2)
        end if
! 
! !Surface values
        fname3(11:13) = cmyid
        fname3(15:17) = cexpnr
        call ncinfo(tncname3(1,:),'time','Time','s','time')
        call ncinfo(ncname3( 1,:),'Qnet','Net radiation','W/m^2','tt0t')
        call ncinfo(ncname3( 2,:),'H','Sensible heat flux','W/m^2','tt0t')
        call ncinfo(ncname3( 3,:),'LE','Latent heat flux','W/m^2','tt0t')
        call ncinfo(ncname3( 4,:),'G0','Ground heat flux','W/m^2','tt0t')
        call ncinfo(ncname3( 5,:),'tskin','Skin temperature','K','tt0t')
        call ncinfo(ncname3( 6,:),'tendskin','Skin tendency','W/m^2','tt0t')
        call ncinfo(ncname3( 7,:),'rs','Surface resistance','s/m','tt0t')
        call ncinfo(ncname3( 8,:),'ra','Aerodynamic resistance','s/m','tt0t')
        call ncinfo(ncname3( 9,:),'cliq','Fraction of vegetated surface covered with liquid water','-','tt0t')
        call ncinfo(ncname3(10,:),'Wl','Liquid water reservoir','m','tt0t')
        call ncinfo(ncname3(11,:),'rssoil','Soil evaporation resistance','s/m','tt0t')
        call ncinfo(ncname3(12,:),'rsveg','Vegitation resistance','s/m','tt0t')
        call open_nc(fname3,  ncid3,nrec3,n1=imax,n2=jmax)
        if (nrec3==0) then
          call define_nc( ncid3, 1, tncname3)
          call writestat_dims_nc(ncid3)
        end if
        call redefine_nc(ncid3)
        call define_nc( ncid3, NVar3, ncname3)
    end if


  end subroutine initlsmcrosssection
!>Run lsmcrosssection. Mainly timekeeping
  subroutine lsmcrosssection
    use modglobal, only : rk3step,timee,rtimee,dt_lim
    use modstat_nc, only : lnetcdf, writestat_nc
    implicit none


    if (.not. lcross) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))
   call wrtvert
   call wrthorz
   call wrtsurf

  end subroutine lsmcrosssection


!> Do the xz lsmcrosssections and dump them to file
  subroutine wrtvert
  use modglobal, only : imax,i1,j1,cexpnr,ifoutput,rtimee
  use modsurfdata, only : tsoil, phiw
  use modmpi,    only : myid
  use modstat_nc, only : lnetcdf, writestat_nc
  implicit none

  integer i,k

  real, allocatable :: vars(:,:,:)

  if( myid /= 0 ) return

    open(ifoutput,file='movv_tsoil.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((tsoil(i,crossplane,k),i=2,i1),k=1,ksoilmax)
    close(ifoutput)

    open(ifoutput,file='movv_phiw.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((phiw(i,crossplane,k),i=2,i1),k=1,ksoilmax)
    close(ifoutput)
    if (lnetcdf) then
      allocate(vars(1:imax,1:ksoilmax,2))
      vars(:,:,1) = tsoil(2:i1,crossplane,1:ksoilmax)
      vars(:,:,2) = phiw(2:i1,crossplane,1:ksoilmax)
      call writestat_nc(ncid1,1,tncname1,(/rtimee/),nrec1,.true.)
      call writestat_nc(ncid1,2,ncname1(1:2,:),vars,nrec1,imax,ksoilmax)
      deallocate(vars)
    end if

  end subroutine wrtvert

!> Do the xy lsmcrosssections and dump them to file
  subroutine wrthorz
    use modglobal, only : imax,jmax,i1,j1,cexpnr,ifoutput,rtimee
    use modsurfdata, only : tsoil,phiw
    use modmpi,    only : cmyid
    use modstat_nc, only : lnetcdf, writestat_nc
    implicit none


    ! LOCAL
    integer i,j
    real, allocatable :: vars(:,:,:)




    write(cheight,'(i4.4)') crossheight
    open(ifoutput,file='movh_tsoil.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((tsoil(i,j,crossheight),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_phiw.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((phiw(i,j,crossheight),i=2,i1),j=2,j1)
    close(ifoutput)


    if (lnetcdf) then
        allocate(vars(1:imax,1:jmax,2))
        vars(:,:,1) = tsoil(2:i1,2:j1,crossheight)
        vars(:,:,2) = phiw(2:i1,2:j1,crossheight)
        call writestat_nc(ncid2,1,tncname2,(/rtimee/),nrec2,.true.)
        call writestat_nc(ncid2,2,ncname2(1:2,:),vars,nrec2,imax,jmax)
        deallocate(vars)
    end if


  end subroutine wrthorz
  !> Do the xy lsmcrosssections and dump them to file
  subroutine wrtsurf
    use modglobal, only : imax,jmax,i1,j1,cexpnr,ifoutput,rtimee
    use modsurfdata, only : Qnet, H, LE, G0, rs, ra, tskin, tendskin, &
                           cliq,rsveg,rssoil,Wl
    use modmpi,    only : cmyid
    use modstat_nc, only : lnetcdf, writestat_nc
    implicit none


    ! LOCAL
    integer i,j
    real, allocatable :: vars(:,:,:)




    open(ifoutput,file='movh_qnet.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((qnet(i,j),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_h.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((h(i,j),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_le.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((le(i,j),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_g0.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((g0(i,j),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_rs.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((rs(i,j),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_ra.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((ra(i,j),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_tskin.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((tskin(i,j),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_tendskin.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((tendskin(i,j),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_cliq.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((cliq(i,j),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_rsveg.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((rsveg(i,j),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_rssoil.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((rssoil(i,j),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_wl.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((wl(i,j),i=2,i1),j=2,j1)
    close(ifoutput)


    if (lnetcdf) then
        allocate(vars(1:imax,1:jmax,nvar3))
        vars(:,:,1) = qnet(2:i1,2:j1)
        vars(:,:,2) = h(2:i1,2:j1)
        vars(:,:,3) = le(2:i1,2:j1)
        vars(:,:,4) = g0(2:i1,2:j1)
        vars(:,:,5) = tskin(2:i1,2:j1)
        vars(:,:,6) = tendskin(2:i1,2:j1)
        vars(:,:,7) = rs(2:i1,2:j1)
        vars(:,:,8) = ra(2:i1,2:j1)
        vars(:,:,9) = cliq(2:i1,2:j1)
        vars(:,:,10) = Wl(2:i1,2:j1)
        vars(:,:,11) = rssoil(2:i1,2:j1)
        vars(:,:,12) = rsveg(2:i1,2:j1)
        call writestat_nc(ncid3,1,tncname3,(/rtimee/),nrec3,.true.)
        call writestat_nc(ncid3,nvar3,ncname3(1:nvar3,:),vars,nrec3,imax,jmax)
        deallocate(vars)
    end if


  end subroutine wrtsurf
!> Clean up when leaving the run
  subroutine exitlsmcrosssection
    use modstat_nc, only : exitstat_nc,lnetcdf
    use modmpi, only : myid
    implicit none

    if(lcross .and. lnetcdf) then
    if (myid==0) then
      call exitstat_nc(ncid1)
    end if
        call exitstat_nc(ncid2)
    end if

  end subroutine exitlsmcrosssection

end module modlsmcrosssection
