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
module modcrosssection

    !-----------------------------------------------------------------|
    !                                                                 |
    !*** *cross*  dumps an instantenous crosssection of the field     |
    !                                                                 |
    !     purpose.                                                    |
    !     --------                                                    |
    !                                                                 |
    !                                                                 |
    !    Crosssections in the yz-plane and in the xy-plane            |
    !        of u,v,w,thl,thv,qt,ql                                   |
    !____________________SETTINGS_AND_SWITCHES________________________|
    !                     IN &NAMCROSSSECTION                         |
    !                                                                 |
    !    dtav        SAMPLING INTERVAL                                |
    !                                                                 |
    !    lcross      SWITCH TO ENABLE CROSSSECTION                    |
    !                                                                 |
    !    crossheight HEIGHT OF THE XY-CROSSSECTION                    |
    !                                                                 |
    !    crossplane  LOCATION OF THE YZ-PLANE ON EVERY PROCESSOR      |
    !-----------------------------------------------------------------|

implicit none
private
PUBLIC :: initcrosssection, crosssection
save

  real    :: dtav,tnext
  logical :: lcross = .false. ! switch for conditional sampling cloud (on/off)
  integer :: crossheight = 2
  integer :: crossplane = 2

contains

  subroutine initcrosssection
    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_logical,mpi_integer
    use modglobal,only :ifnamopt,fname_options,dtmax,rk3step, dtav_glob,ladaptive,j1,kmax,dt_lim
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


  end subroutine initcrosssection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine crosssection
    use modglobal, only : rk3step,timee,dt_lim
    implicit none


    if (.not. lcross) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+dtav
    dt_lim = minval((/dt_lim,tnext-timee/))
    call wrtvert
    call wrthorz


  end subroutine crosssection
  !*****************************************************************************
  subroutine wrtvert
  use modglobal, only : i1,j1,kmax,nsv,rlv,cp,rv,rd,cu,cv,cexpnr,ifoutput
  use modfields, only : um,vm,wm,thlm,qtm,svm,thl0,qt0,ql0,exnf
  use modmpi,    only : myid
  implicit none

  integer i,j,k,n
  character(20) :: name

  real, allocatable :: thv0(:,:,:)

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

    deallocate(thv0)

  end subroutine wrtvert
!*****************************************************************************
  subroutine wrthorz
    use modglobal, only : i1,j1,nsv,rlv,cp,rv,rd,cu,cv,cexpnr,ifoutput
    use modfields, only : um,vm,wm,thlm,qtm,svm,thl0,qt0,ql0,exnf
    use modmpi,    only : cmyid
    implicit none


    ! LOCAL
    integer i,j,n
    character(20) :: name

    real, allocatable :: thv0(:,:)

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

    deallocate(thv0)

  end subroutine wrthorz
end module modcrosssection