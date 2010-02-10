!> \file modtilt.f90
!!  Calculates flow over a slope

!>
!!  Calculates flow over a slope
!>
!!  \author Simon Axelsen, UU
!! \todo documentation
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

    !-----------------------------------------------------------------|
    !                                                                 |
    !* module tilt calculates prognostics of momentum and temperature |
    !     due to a sloping surface                                    |
    !                                                                 |
    !      Simon Axelsen                19/11/2007                    |
    !                                                                 |
    !____________________SETTINGS_AND_SWITCHES________________________|
    !                      IN &MODTILT                                |
    !                                                                 |
    !    ltilted        SWITCHES THE ROUTINE ON/OFF                   |
    !                                                                 |
    !    alfa           SLOPE ANGLE                                   |
    !                                                                 |
    !    dtav           TIME INTERVAL FOR SAMPLING OF STATISTICS      |
    !                                                                 |
    !    timeav         TIME INTERVAL FOR WRITING OF STATISTICS       |
    !                                                                 |
    !                                                                 |
    !  CONTENTS:                                                      |
    !          1.1 inittilt: initialize variables                     !
    !          2.1 initthla: construct or read thla                   !
    !          3.1 gravityforce: change in momentum due to slope      !
    !          4.1 caldefh: calculate deficit at half-level           !
    !          5.1 thetatodelta: pot temp to pot temp deficit         !
    !          5.2 thetatodelta: pot temp deficit to pot temp         !
    !          6.1 writethla: write the background pot.temp           !
    !          6.2 readthla: read the background pot.temp             !
    !          7.1 tiltstat: make call to dostat if time=timeav       !
    !          7.2 dostat: calcualte statistics on deficit            !
    !          7.3 writestat: write the statistics                    !
    !          8.1 exittilt: deacllocte variables                     !
    !-----------------------------------------------------------------|
module modtilt
implicit none
PRIVATE
PUBLIC :: thldefm,thldef0, ltilted,&
         inittilt,initthla,tiltedgravity,tiltedboundary,adjustbudget,&
         tiltstat,deltatotheta,thetatodelta,exittilt,thla
SAVE

  real    :: alfa     = 0.
  logical :: ltilted  = .false.
  logical :: lstat    = .true. !default = true, but subroutine only entered if ltilted=true
  real    :: dtav, timeav
  integer :: idtav,itimeav,tnext,tnextwrite
  integer :: nsamples



  !Prognostic / diagnostic variables
  real, allocatable :: thldefm(:,:,:)   !  liq. water pot. temp. deficit at time step t-1
  real, allocatable :: thldef0(:,:,:)   !  liq. water pot. temp. deficit at time step t
  real, allocatable :: thldef0h(:,:,:)  !  thl deficit at half levels for gravity x-direction
  real, allocatable :: thla(:,:,:)      !  liq.water pot.temp of ambient air (only if tilted)

  !Statistical variables
  real, allocatable :: thldefmn(:)      ! time averaged temp. deficit
  real, allocatable :: thldefav(:)      ! slab averaged temp. deficit
  real, allocatable :: thlaav(:)        ! slab average of ambient pot.temp


contains

!*****************************************************************************
!  1.1 INITTILT
!*****************************************************************************
  subroutine inittilt
    use modmpi,   only : myid,my_real,mpierr,comm3d,mpi_integer,mpi_logical,slabsum
    use modglobal,only : ifnamopt,ifoutput,fname_options,ifinput,cexpnr,&
                         i1,j1,ih,jh,k1,iexpnr,cu,lmoist,timeav_glob, dtav_glob, dt_lim,btime,tres
    use modstartup,only : irandom,randthl,krand,randomnize

    implicit none


  ! LOCAL
    integer :: k, ierr

    namelist/NAMTILT/  &
         ltilted,alfa,lstat,dtav,timeav

    dtav=dtav_glob;timeav=timeav_glob
    if(myid==0)then
        open(ifnamopt,file=fname_options,status='old',iostat=ierr)
        read (ifnamopt,NAMTILT,iostat=ierr)
        if (ierr > 0) then
          print *, 'Problem in namoptions NAMTILT'
          print *, 'iostat error: ', ierr
          stop 'ERROR: Problem in namoptions NAMTILT'
        endif
        write(6 ,NAMTILT)
      close(ifnamopt)
    end if



    call MPI_BCAST(ltilted ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lstat   ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(alfa    ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(dtav    ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(timeav  ,1,MY_REAL    ,0,comm3d,mpierr)
    idtav = dtav/tres
    dtav  = idtav*tres
    itimeav = timeav/tres
    timeav  = itimeav*tres

    tnext      = dtav   +btime
    tnextwrite = timeav +btime
    nsamples = itimeav/idtav
    if(.not.(ltilted)) return

    if(cu/=0.) stop 'cu/=0 not allowed in a tilted environment'
    if(lmoist) stop 'no humidity allowed in a tilted environment'

    allocate(thldefm(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thldef0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thldef0h(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thla(2-ih:i1+ih,2-jh:j1+jh,k1))

    allocate(thldefmn(k1))
    allocate(thldefav(k1))
    allocate(thlaav(k1))

    thldefm=0.;thldef0=0.;thldef0h=0.;thla=0.
    thldefmn=0.;thldefav=0.
    call initthla
    do k=1,krand
      call randomnize(thldefm,k,randthl,irandom,ih,jh)
      call randomnize(thldef0,k,randthl,irandom,ih,jh)
    enddo
    call deltatotheta
    call initthla
    if(lstat) then
      if(myid==0)then
        open (ifoutput,file='fielddef.'//cexpnr,status='replace')
        close (ifoutput)
      endif

      dt_lim = min(dt_lim,tnext)

    endif


  end subroutine inittilt

!*****************************************************************************
!  2.1 INITTHLA
!*****************************************************************************
  subroutine initthla

    use modfields, only : thlprof
    use modglobal, only : i1,j1,ih,jh,k1,imax,kmax,dzf,dx,lwarmstart
    use modmpi,    only : myid

    implicit none
    integer i,j,k
    real dthl0dz !individual for all k

    if(.not.(ltilted)) return

    if(.not.(lwarmstart)) then
       do j=2-jh,j1+jh
       do i=2-ih,i1+ih
          !k=1
          dthl0dz = (thlprof(2)-thlprof(1))/dzf(1)
          thla(i,j,1) = thlprof(1)-(i-1)*dx*dthl0dz*sin(alfa)
          do k=2,kmax
             dthl0dz = (thlprof(k)-thlprof(k-1))/dzf(k)
             thla(i,j,k) = thlprof(1)-(i-1)*dx*dthl0dz*sin(alfa)
             thla(i,j,k) =  thla(i,j,k-1)+dzf(k)*dthl0dz*cos(alfa)
          enddo
          thla(i,j,k1) = thla(i,j,kmax)
       enddo
       enddo


       else !.not.lwarmstart
          call readthla
          call thetatodelta

       endif

       !Save thla field to unformatted file, used for restart
       call writethla

       if(lstat)then
          !thla is independent of time. Calculate slab-average thla
          do k=1,k1
             thlaav(k) = thla(imax/2,2,k)
          enddo
       endif

  end subroutine initthla


  subroutine tiltedboundary
  use modmpi,    only : excjs
  use modglobal, only : i1,j1,k1,ih,jh,i2,kmax

  implicit none
  integer m
  if(ltilted)then
    call thetatodelta
    do m=1,ih
       thldef0(2-m,2:j1,1:k1)   = thldef0(i2-m,2:j1,1:k1)
       thldef0(i1+m,2:j1,1:k1)  = thldef0(1+m,2:j1,1:k1)
       thldefm(2-m,2:j1,1:k1)   = thldefm(i2-m,2:j1,1:k1)
       thldefm(i1+m,2:j1,1:k1)  = thldefm(1+m,2:j1,1:k1)
    enddo
    call excjs( thldef0           , 2,i1,2,j1,1,k1,ih,jh)
    call excjs( thldefm           , 2,i1,2,j1,1,k1,ih,jh)
    thldef0(2:i1,2:j1,k1) = thldef0(2:i1,2:j1,kmax)
    thldefm(2:i1,2:j1,k1) = thldefm(2:i1,2:j1,kmax)
    call deltatotheta !Revert back to pot temp.
  endif
  !CvH


  end subroutine tiltedboundary

!*****************************************************************************
!  3.1 GRAVITYFORCE
!*****************************************************************************
  subroutine tiltedgravity

    use modglobal,  only : grav,kmax,i1,j1
    use modsurfdata,only : thvs
    use modfields,  only : thv0h,up,wp

    implicit none
    integer i,j,k

    if(.not.(ltilted)) return

    call caldefh
    do k=1,kmax
    do j=2,j1
    do i=2,i1
       up(i,j,k) = up(i,j,k) - (grav/thvs)*thldef0h(i,j,k)*sin(alfa)
       wp(i,j,k) = wp(i,j,k) - (grav/thvs)*thv0h(i,j,k) !remove contribution from forces
       wp(i,j,k) = wp(i,j,k) + (grav/thvs)*thldef0h(i,j,k)*cos(alfa)
       !cos(alfa)-1: remove contribution from modforces and add tilted contribution
    enddo
    enddo
    enddo
    !Overwrite tilting correction to wp at the surface
    wp(2:i1,2:j1,1) = 0.0

  end subroutine tiltedgravity
!*****************************************************************************
!  4.1 CALDEFH
!*****************************************************************************

  subroutine caldefh

    use modglobal, only : dzf,dzh,&
                          i1,j1,k1,&
                          iadv_kappa,iadv_thl

    implicit none
    integer i,j,k

    if(.not.(ltilted)) return

    if (iadv_thl==iadv_kappa) then
       !Special ifattention temperature is advected using kappa scheme
       call halflev_kappa(thldef0,thldef0h)
    else
       do  j=2,j1
       do  i=2,i1
          do  k=2,k1
             thldef0h(i,j,k) = (thldef0(i,j,k)*dzf(k-1)+thldef0(i,j,k-1)*dzf(k))/(2*dzh(k))
          end do
          thldef0h(i,j,1) = thldef0(i,j,1)
       end do
       end do

    endif



  end subroutine caldefh


!*****************************************************************************
!  5. THETATODELTA AND DELTATOTHETA
!*****************************************************************************
!------------------------
! 5.1 THETATODELTA
!------------------------
  subroutine thetatodelta
    use modfields, only : thl0,thlm
    use modglobal, only : i1,j1,ih,jh,k1

    implicit none
    integer i,j,k

    if(.not.(ltilted)) return

       do k=1,k1
       do i=2-ih,i1+ih
       do j=2-jh,j1+jh
          thldef0(i,j,k) = thl0(i,j,k)-thla(i,j,k)
          thldefm(i,j,k) = thlm(i,j,k)-thla(i,j,k)
       enddo
       enddo
       enddo

     end subroutine thetatodelta

!------------------------
! 5.2 DELTTOTHETA
!------------------------
  subroutine deltatotheta

    use modfields, only : thl0,thlm
    use modglobal, only : i1,j1,ih,jh,k1

    implicit none
    integer i,j,k

    if(.not.(ltilted)) return

       do k=1,k1
       do i=2-jh,i1+ih
       do j=2-jh,j1+jh
          thl0(i,j,k) = thldef0(i,j,k)+thla(i,j,k)
          thlm(i,j,k) = thldefm(i,j,k)+thla(i,j,k)
       enddo
       enddo
       enddo

     end subroutine deltatotheta


!*****************************************************************************
!  6. READTHLA AND WRITETHLA
!*****************************************************************************
!------------------------
! 6.1 READTHLA
!------------------------

!-------------- Read thla--------------------
  subroutine readthla

    use modglobal, only : i1,i2,ih,j1,j2,jh,k1,iexpnr,ifinput
    use modmpi,    only : myid
    implicit none
   !thla homogeneous in y-direction, e.g. same input for all nprocs
    !variable myid only used for writing to log

    character(20),parameter :: fname = 'thlafield'
    integer i,j,k,ios
    logical :: lexists

    !Checking the existance of the file thlafield
    inquire(file=fname, EXIST=lexists)
    if(lexists) then
       write(6,*) 'Opening, myid ',fname,myid
       open(unit=ifinput,file=fname,form='unformatted', status='old',IOSTAT=ios)
       if(ios==0)then
          read(ifinput)  (((thla    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
          read(ifinput)  alfa
       else !Error opening thlafield, cannot continue
          stop 'Thlafield does exist, but could not be opened. Stopping'
       endif
    endif

    close(ifinput)
  end subroutine readthla

!------------------------
! 6.2 WRITETHLA
!------------------------

  subroutine writethla

    use modglobal, only : i1,i2,ih,j1,j2,jh,k1,iexpnr,ifoutput
    use modmpi,    only : myid
    implicit none
    !thla homogeneous in y-direction, e.g. same input for all nprocs
    !Write field only if myid==0, and if file does not already exist

    character(20),parameter :: fname = 'thlafield'
    integer i,j,k,ios
    logical :: lexists

    if(myid==0)then
       inquire(file=fname, EXIST=lexists)
       if(lexists) then
          return !file already exists
       else
          write(*,*) 'Thlafield does not yet exist, open it'
          open(unit=ifoutput,file=fname,form='unformatted', status='new',IOSTAT=ios)
          if(ios==0) then
             write(ifoutput)  (((thla    (i,j,k),i=2-ih,i1+ih),j=2-jh,j1+jh),k=1,k1)
             write(ifoutput)  alfa
          else !problems opening output file, stopping
             stop 'An error occured when trying to open thlafield for writing, stopping'
          endif
       end if
       close(ifoutput)
    end if

  end subroutine writethla

!*****************************************************************************
!  7. TILTSTAT & DO_STAT & WRITE_STAT
!*****************************************************************************
!------------------------
! 7.1 Tiltstat
!     Tiltstat is called from program
!------------------------

  subroutine tiltstat
    use modglobal, only : rk3step,timee,dt_lim
    implicit none
    if (.not.(ltilted)) return
    if (.not.(lstat))   return
    if (rk3step/=3)     return
    if(timee<tnext .and. timee<tnextwrite) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if
    if (timee>=tnext) then
      tnext = tnext+idtav
      call do_tiltstat
    end if
    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+itimeav
      call writetiltstat
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
  end subroutine tiltstat

!------------------------
! 7.2 Do_stat
!------------------------

  subroutine do_tiltstat

    use modglobal, only : i1,ih,j1,jh,k1,kmax,rslabs
    use modmpi,    only : nprocs,comm3d,nprocs,my_real, mpi_sum,mpierr, slabsum


    implicit none

    !Calculate slab average
    call slabsum(thldefav,1,k1,thldefm,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    thldefav = thldefav/rslabs

    !Add slab average to time mean
    thldefmn  = thldefmn + thldefav




  end subroutine do_tiltstat

!------------------------
! 7.3 Write_tiltstat
!------------------------
  subroutine writetiltstat

    use modglobal, only : kmax,k1,zf,rtimee,cexpnr,ifoutput
    use modmpi,    only : myid

    implicit none
    integer nsecs, nhrs, nminut,k
    real convt
    nsecs   = nint(rtimee)
    nhrs    = int(nsecs/3600)
    nminut  = int(nsecs/60)-nhrs*60
    nsecs   = mod(nsecs,60)
    convt   = 86400.

    thldefmn  = thldefmn  *nsamples


    if(myid==0)then
       open (ifoutput,file='fielddef.'//cexpnr,position='append')
       write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
            '#--------------------------------------------------------'      &
            ,'#',(timeav),'--- AVERAGING TIMESTEP --- '      &
            ,nhrs,':',nminut,':',nsecs      &
            ,'   HRS:MIN:SEC AFTER INITIALIZATION '
       write (ifoutput,'(A/A/A)') &
            '#--------------------------------------------------------' &
            ,'#LEV  HGHT     THLDEF        THLA     ' &
            ,'#      (M)    (---- (KELVIN) --_------)    '
       do k=1,kmax
          write(ifoutput,'(I3,F8.2,2F10.4)') &
               k,          &
               zf      (k), &
               thldefmn(k),&
               thlaav  (k)
       end do
       close (ifoutput)
    endif

    thldefmn=0.

  end subroutine writetiltstat




!*****************************************************************************
!  8.1 ADJUSTBUDGET
!*****************************************************************************
! The resolved tke budget is altered by the slope, e.g. the buoyancy term

  subroutine adjustbudget(buoyav)
    use modglobal,  only : j1,i1,rslabs,k1,cu,grav
    use modfields,  only : u0av,u0
    use modsurface, only : thvs

    implicit none
    real buox
    real buoyav(k1)
    integer i,j,k

    do k=2,k1
       buox=0.
       do j=2,j1
       do i=2,i1
          buox = buox + u0(i,j,k)*thldef0(i,j,k)
       enddo
       enddo
       buoyav(k) = buoyav(k)*cos(alfa) & !Tilt term in modbudget
                 + (grav/thvs)*(-buox/rslabs+(u0av(k)-cu)*thldefav(k))*sin(alfa)
    enddo

  end subroutine adjustbudget

!*****************************************************************************
!  9.1 EXITTILT
!*****************************************************************************
  subroutine exittilt
    implicit none

    if (.not.(ltilted)) return

    deallocate(thldefm,thldef0,thldef0h,thla)
    if(lstat) deallocate(thlaav,thldefav,thldefmn)

  end subroutine exittilt

!*****************************************************************************


!*****************************************************************************
end module modtilt

