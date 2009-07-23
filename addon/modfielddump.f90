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
module modfielddump

    !-----------------------------------------------------------------|
    !                                                                 |
    !*** *fielddump*  dumps complete 3d fields in 2-byte integers     |
    !                                                                 |
    !      Thijs Heus                   19/06/2007                    |
    !                                                                 |
    !     purpose.                                                    |
    !     --------                                                    |
    !                                                                 |
    !                                                                 |
    !    Dumps fields of:                                             |
    !                u,v,w, thl,thv,qt and ql                         |
    !____________________SETTINGS_AND_SWITCHES________________________|
    !                     IN &NAMFIELDDUMP                             |
    !                                                                 |
    !    dtav            SAMPLING INTERVAL                             |
    !    lfielddump      SWITCH TO ENABLE FIELDDUMP                    |
    !    ldiracc         SWITCH TO DUMP IN DIRECT ACCESS FILES         |
    !-----------------------------------------------------------------|

implicit none
private
PUBLIC :: initfielddump, fielddump
save

  real    :: dtav,tnext
  integer :: klow,khigh
  logical :: lfielddump= .false. ! switch for conditional sampling cloud (on/off)
  logical :: ldiracc   = .false. ! switch for conditional sampling cloud (on/off)

contains

  subroutine initfielddump
    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_logical,mpi_integer
    use modglobal,only :ifnamopt,fname_options,dtmax,dtav_glob,kmax, ladaptive,dt_lim,btime
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


  end subroutine initfielddump
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fielddump
    use modfields, only : um,vm,wm,thlm,qtm,ql0
    use modsurfdata,only : thls,qts,thvs
    use modglobal, only : imax,i1,ih,jmax,j1,jh,kmax,k1,rk3step,&
                          timee,dt_lim,cexpnr,ifoutput
    use modmpi,    only : myid,cmyid
    implicit none

    integer(KIND=selected_int_kind(4)), allocatable :: field(:,:,:)
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


    reclength = imax*jmax*(khigh-klow+1)*2

    field = NINT(1.0E3*um,2)
    if (ldiracc) then
      open (ifoutput,file='wbuu.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
      write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
    else
      open  (ifoutput,file='wbuu.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
      write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
    end if
    close (ifoutput)

    field = NINT(1.0E3*vm,2)
    if (ldiracc) then
      open (ifoutput,file='wbvv.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
      write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
    else
      open  (ifoutput,file='wbvv.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
      write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
    end if
    close (ifoutput)

    field = NINT(1.0E3*wm,2)
    if (ldiracc) then
      open (ifoutput,file='wbww.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
      write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
    else
      open  (ifoutput,file='wbww.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
      write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
    end if
    close (ifoutput)

    field = NINT(1.0E5*qtm,2)
    if (ldiracc) then
      open (ifoutput,file='wbqt.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
      write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
    else
      open  (ifoutput,file='wbqt.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
      write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
    end if
    close (ifoutput)

    field = NINT(1.0E5*ql0,2)
    if (ldiracc) then
      open (ifoutput,file='wbql.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
      write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
    else
      open  (ifoutput,file='wbql.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
      write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
    end if
    close (ifoutput)

    field = NINT(1.0E3*(thlm-thls),2)
    if (ldiracc) then
      open (ifoutput,file='wbtl.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
      write (ifoutput, rec=writecounter) field(2:i1,2:j1,klow:khigh)
    else
      open  (ifoutput,file='wbtl.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
      write (ifoutput) (((field(i,j,k),i=2,i1),j=2,j1),k=klow,khigh)
    end if
    close (ifoutput)

    if (myid==0) then
      open(ifoutput, file='wbthls.'//cexpnr,form='formatted',position='append')
      write(ifoutput,'(F12.1 3F12.5)') timee,thls, qts,thvs
      close(ifoutput)
    end if
    writecounter=writecounter+1

    deallocate(field)

  end subroutine fielddump

end module modfielddump