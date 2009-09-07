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

implicit none
private
PUBLIC :: initprojection, projection
save

  real    :: dtav,tnext
  logical :: lproject = .false. ! switch for conditional sampling cloud (on/off)
  integer :: projectheight = 2  !lowest integration boundary
  integer :: projectplane = 2

contains

  subroutine initprojection
    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_logical,mpi_integer
    use modglobal,only :ifnamopt,fname_options,dtmax,rk3step, dtav_glob,ladaptive,j1,kmax,dt_lim
    implicit none

    integer :: ierr

    namelist/NAMprojection/ &
    lproject, dtav, projectheight, projectplane

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMprojection,iostat=ierr)
      write(6 ,NAMprojection)
      close(ifnamopt)
    end if

    call MPI_BCAST(dtav       ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(lproject     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(projectheight,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(projectplane ,1,MPI_INTEGER,0,comm3d,mpierr)
    tnext   = dtav-1e-3
    if(.not.(lproject)) return
    dt_lim = min(dt_lim,tnext)

    if(projectheight>kmax .or. projectplane>j1) then
      stop 'projection: projection out of range' 
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'projection: dtav should be a integer multiple of dtmax'
    end if


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
    tnext = tnext+dtav
    dt_lim = minval((/dt_lim,tnext-timee/))
    call wrtvert
    call wrthorz


  end subroutine projection
  !*****************************************************************************
  subroutine wrtvert
  use modglobal, only : i1,j1,kmax,nsv,cu,cv,cexpnr,ifoutput
  use modfields, only : thv0h,thlm,qtm,svm,ql0
  use modmpi,    only : myid
  implicit none

  integer i,j,k,n
  character(20) :: name

 
  if( myid /= 0 ) return

    open(ifoutput,file='projv_thl.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es11.4)') ((sum(thlm(i,projectplane:j1,k)),i=2,i1),k=1,kmax)
    close(ifoutput)

    open(ifoutput,file='projv_thv.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es11.4)') ((sum(thv0h(i,projectplane:j1,k)),i=2,i1),k=1,kmax)
    close(ifoutput)

    open(ifoutput,file='projv_qt.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es11.4)') ((1.e3*sum(qtm(i,projectplane:j1,k)),i=2,i1),k=1,kmax)
    close(ifoutput)

    open(ifoutput,file='projv_ql.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es11.4)') ((1.e3*sum(ql0(i,projectplane:j1,k)),i=2,i1),k=1,kmax)
    close(ifoutput)

    do n = 1,nsv
      name = 'projh_tnn.'//cexpnr
      write(name(7:8),'(i2.2)') n
      open(ifoutput,file=name,position='append',action='write')
      write(ifoutput,'(es11.4)') ((sum(svm(i,projectplane:j1,k,n)),i=2,i1),k=1,kmax)
      close(ifoutput)
    end do


  end subroutine wrtvert
!*****************************************************************************
  subroutine wrthorz
    use modglobal, only : i1,j1,kmax,nsv,cexpnr,ifoutput
    use modfields, only : thlm,qtm,svm,thv0h,ql0,thl0av,qt0av
    use modmpi,    only : cmyid
    implicit none


    ! LOCAL
    integer i,j,n
    character(20) :: name

    open(ifoutput,file='projh_thl.'//cmyid//'.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((sum(thlm(i,j,projectheight:kmax)-thl0av(projectheight:kmax)),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='projh_thv.'//cmyid//'.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((sum(thv0h(i,j,projectheight:kmax)-thl0av(projectheight:kmax)),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='projh_qt.'//cmyid//'.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((1.e3*sum(qtm(i,j,projectheight:kmax)-qt0av(projectheight:kmax)),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='projh_ql.'//cmyid//'.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((1.e3*sum(ql0(i,j,projectheight:kmax)),i=2,i1),j=2,j1)
    close(ifoutput)

    do n = 1,nsv
      name = 'projh_snn.'//cmyid//'.'//cexpnr
      write(name(7:8),'(i2.2)') n
      open(ifoutput,file=name,position='append',action='write')
      write(ifoutput,'(es12.5)') ((sum(svm(i,j,projectheight:kmax,n)),i=2,i1),j=2,j1)
      close(ifoutput)
    end do

  end subroutine wrthorz
end module modprojection