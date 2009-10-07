!> \file modcrosssection.f90
!!   Dumps an instantenous crosssection of the field

!>
!! Dumps an instantenous crosssection of the field.
!>
!! Crosssections in the yz-plane and in the xy-plane            |
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
module modcrosssection



implicit none
private
PUBLIC :: initcrosssection, crosssection,exitcrosssection
save
!NetCDF variables
  integer,parameter :: nvar = 14
  integer :: ncid,nrec = 0
  character(80) :: fname = 'cross.xxx.xxx.nc'
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav,tnext
  logical :: lcross = .false. !< switch for doing the crosssection (on/off)
  integer :: crossheight = 2 !< Height of the xy crosssection
  integer :: crossplane = 2 !< Location of the xz crosssection

contains
!> Initializing Crosssection. Read out the namelist, initializing the variables
  subroutine initcrosssection
    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_logical,mpi_integer,cmyid
    use modglobal,only :imax,jmax,ifnamopt,fname_options,dtmax,rk3step, dtav_glob,ladaptive,j1,kmax,dt_lim,cexpnr
    use modstat_nc,only : lnetcdf,open_nc, define_nc, redefine_nc,ncinfo,writestat_dims_nc
   implicit none

    integer :: ierr

    namelist/NAMCROSSSECTION/ &
    lcross, dtav, crossheight, crossplane

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMCROSSSECTION,iostat=ierr)
      write(6 ,NAMCROSSSECTION)
      close(ifnamopt)
    end if

    call MPI_BCAST(dtav       ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(lcross     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(crossheight,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(crossplane ,1,MPI_INTEGER,0,comm3d,mpierr)
    tnext   = dtav-1e-3
    if(.not.(lcross)) return
    dt_lim = min(dt_lim,tnext)

    if(crossheight>kmax .or. crossplane>j1) then
      stop 'CROSSSECTION: crosssection out of range'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'CROSSSECTION: dtav should be a integer multiple of dtmax'
    end if
    if (lnetcdf) then
      dtav = dtav_glob
      fname(7:9) = cmyid
      fname(11:13) = cexpnr
      call ncinfo(tncname(1,:),'time','Time','s','time')
      call ncinfo(ncname( 1,:),'uxz', 'xz crosssection of the West-East velocity','m/s','m0tt')
      call ncinfo(ncname( 2,:),'vxz', 'xz crosssection of the South-North velocity','m/s','t0tt')
      call ncinfo(ncname( 3,:),'wxz', 'xz crosssection of the Vertical velocity','m/s','t0mt')
      call ncinfo(ncname( 4,:),'thlxz','xz crosssection of the Liquid water potential temperature','K','t0tt')
      call ncinfo(ncname( 5,:),'thvxz','xz crosssection of the Virtual potential temperature','K','t0tt')
      call ncinfo(ncname( 6,:),'qtxz','xz crosssection of the Total water mixing ratio','kg/kg','t0tt')
      call ncinfo(ncname( 7,:),'qlxz','xz crosssection of the Liquid water mixing ratio','kg/kg','t0tt')
      call ncinfo(ncname( 8,:),'uxy','xy crosssection of the West-East velocity','m/s','mt0t')
      call ncinfo(ncname( 9,:),'vxy','xy crosssection of the South-North velocity','m/s','tm0t')
      call ncinfo(ncname(10,:),'wxy','xy crosssection of the Vertical velocity','m/s','tt0t')
      call ncinfo(ncname(11,:),'thlxy','xy crosssection of the Liquid water potential temperature','K','tt0t')
      call ncinfo(ncname(12,:),'thvxy','xy crosssection of the Virtual potential temperature','K','tt0t')
      call ncinfo(ncname(13,:),'qtxy','xy crosssection of the Total water mixing ratio','kg/kg','tt0t')
      call ncinfo(ncname(14,:),'qlxy','xy crosssection of the Liquid water mixing ratio','kg/kg','tt0t')

      call open_nc(fname,  ncid,n1=imax,n2=jmax,n3=kmax)
      call define_nc( ncid, 1, tncname)
      call writestat_dims_nc(ncid)
      call redefine_nc(ncid)
      call define_nc( ncid, NVar, ncname)
    end if


  end subroutine initcrosssection
!>Run crosssection. Mainly timekeeping
  subroutine crosssection
    use modglobal, only : rk3step,timee,dt_lim
    use modstat_nc, only : lnetcdf, writestat_nc
    implicit none


    if (.not. lcross) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+dtav
    dt_lim = minval((/dt_lim,tnext-timee/))
     if(lnetcdf) call writestat_nc(ncid,1,tncname,(/timee/),nrec,.true.)
   call wrtvert
    call wrthorz


  end subroutine crosssection


!> Do the xz crosssections and dump them to file
  subroutine wrtvert
  use modglobal, only : imax,i1,j1,kmax,nsv,rlv,cp,rv,rd,cu,cv,cexpnr,ifoutput
  use modfields, only : um,vm,wm,thlm,qtm,svm,thl0,qt0,ql0,exnf
  use modmpi,    only : myid
  use modstat_nc, only : lnetcdf, writestat_nc
  implicit none

  integer i,j,k,n
  character(20) :: name

  real, allocatable :: thv0(:,:,:),vars(:,:,:)

  if( myid /= 0 ) return

  allocate(thv0(2:i1,2:j1,1:kmax))


    do  j=2,j1
    do  i=2,i1
    do  k=1,kmax
      thv0(i,j,k) = (thl0(i,j,k)+rlv*ql0(i,j,k)/(cp*exnf(k))) &
                    *(1+(rv/rd-1)*qt0(i,j,k)-rv/rd*ql0(i,j,k))
    enddo
    enddo
    enddo

    open(ifoutput,file='movv_u.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es11.4)') ((um(i,crossplane,k)+cu,i=2,i1),k=1,kmax)
    close(ifoutput)

    open(ifoutput,file='movv_v.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es11.4)') ((vm(i,crossplane,k)+cv,i=2,i1),k=1,kmax)
    close(ifoutput)

    open(ifoutput,file='movv_w.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es11.4)') ((wm(i,crossplane,k),i=2,i1),k=1,kmax)
    close(ifoutput)

    open(ifoutput,file='movv_thl.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es11.4)') ((thlm(i,crossplane,k),i=2,i1),k=1,kmax)
    close(ifoutput)

    open(ifoutput,file='movv_thv.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es11.4)') ((thv0(i,crossplane,k),i=2,i1),k=1,kmax)
    close(ifoutput)

    open(ifoutput,file='movv_qt.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es11.4)') ((1.e3*qtm(i,crossplane,k),i=2,i1),k=1,kmax)
    close(ifoutput)

    open(ifoutput,file='movv_ql.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es11.4)') ((1.e3*ql0(i,crossplane,k),i=2,i1),k=1,kmax)
    close(ifoutput)

    do n = 1,nsv
      name = 'movh_tnn.'//cexpnr
      write(name(7:8),'(i2.2)') n
      open(ifoutput,file=name,position='append',action='write')
      write(ifoutput,'(es11.4)') ((svm(i,crossplane,k,n),i=2,i1),k=1,kmax)
      close(ifoutput)
    end do
    if (lnetcdf) then
      allocate(vars(1:imax,1:kmax,7))
      vars(:,:,1) = um(2:i1,crossplane,1:kmax)+cu
      vars(:,:,2) = vm(2:i1,crossplane,1:kmax)+cv
      vars(:,:,3) = wm(2:i1,crossplane,1:kmax)
      vars(:,:,4) = thlm(2:i1,crossplane,1:kmax)
      vars(:,:,5) = thv0(2:i1,crossplane,1:kmax)
      vars(:,:,6) = qtm(2:i1,crossplane,1:kmax)
      vars(:,:,7) = ql0(2:i1,crossplane,1:kmax)
      call writestat_nc(ncid,7,ncname(1:7,:),vars,nrec,imax,kmax)
      deallocate(vars)
    end if
    deallocate(thv0)

  end subroutine wrtvert

!> Do the xy crosssections and dump them to file
  subroutine wrthorz
    use modglobal, only : imax,jmax,i1,j1,nsv,rlv,cp,rv,rd,cu,cv,cexpnr,ifoutput
    use modfields, only : um,vm,wm,thlm,qtm,svm,thl0,qt0,ql0,exnf
    use modmpi,    only : cmyid
    use modstat_nc, only : lnetcdf, writestat_nc
    implicit none


    ! LOCAL
    integer i,j,n
    character(20) :: name

    real, allocatable :: thv0(:,:),vars(:,:,:)

    allocate(thv0(2:i1,2:j1))

    do  j=2,j1
    do  i=2,i1
      thv0(i,j) =&
       (thl0(i,j,crossheight)+&
       rlv*ql0(i,j,crossheight)/&
       (cp*exnf(crossheight))) &
                    *(1+(rv/rd-1)*qt0(i,j,crossheight)&
                    -rv/rd*ql0(i,j,crossheight))
    enddo
    enddo
    open(ifoutput,file='movh_u.'//cmyid//'.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((um(i,j,crossheight)+cu,i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_v.'//cmyid//'.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((vm(i,j,crossheight)+cv,i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_w.'//cmyid//'.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((wm(i,j,crossheight),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_thl.'//cmyid//'.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((thlm(i,j,crossheight),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_thv.'//cmyid//'.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((thv0(i,j),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_qt.'//cmyid//'.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((1.e3*qtm(i,j,crossheight),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_ql.'//cmyid//'.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((1.e3*ql0(i,j,crossheight),i=2,i1),j=2,j1)
    close(ifoutput)

    do n = 1,nsv
      name = 'movh_tnn.'//cmyid//'.'//cexpnr
      write(name(7:8),'(i2.2)') n
      open(ifoutput,file=name,position='append',action='write')
      write(ifoutput,'(es12.5)') ((svm(i,j,crossheight,n),i=2,i1),j=2,j1)
      close(ifoutput)
    end do
    if (lnetcdf) then
      allocate(vars(1:imax,1:jmax,7))
      vars(:,:,1) = um(2:i1,2:j1,crossheight)+cu
      vars(:,:,2) = vm(2:i1,2:j1,crossheight)+cv
      vars(:,:,3) = wm(2:i1,2:j1,crossheight)
      vars(:,:,4) = thlm(2:i1,2:j1,crossheight)
      vars(:,:,5) = thv0(2:i1,2:j1)
      vars(:,:,6) = qtm(2:i1,2:j1,crossheight)
      vars(:,:,7) = ql0(2:i1,2:j1,crossheight)
      call writestat_nc(ncid,7,ncname(8:14,:),vars,nrec,imax,jmax)
      deallocate(vars)
    end if

    deallocate(thv0)

  end subroutine wrthorz
!> Clean up when leaving the run
  subroutine exitcrosssection
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(lcross .and. lnetcdf) call exitstat_nc(ncid)
  end subroutine exitcrosssection

end module modcrosssection