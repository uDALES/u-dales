!> \file modtimedep.f90
!!  Prescribes surface values, fluxes and LS forcings at certain times

!>
!!  Prescribes surface values, fluxes and LS forcings at certain times
!>
!!  \author Roel Neggers, KNMI
!!  \author Thijs Heus,MPI-M
!!  \author Stephan de Roode, TU Delft
!!  \author Simon Axelsen, UU
!!  \par Revision list
!! \todo documentation
!  This file is part of DALES.
!
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

module modtimedep

use mpi

implicit none
private
public :: inittimedep, timedep, ltimedep, ltimedepsurf, ltimedepnudge, ltimedeplw, ltimedepsw, &
          ntimedepsurf, ntimedepnudge, ntimedeplw, ntimedepsw, exittimedep

save
! switches for timedependent surface fluxes and large scale forcings
  logical       :: ltimedep      = .false.  !< Overall switch
  logical       :: ltimedepsurf  = .false.  !< Switch for fluid BC fluxes
  logical       :: ltimedepnudge = .false.  !< Switch for nudging profiles
  logical       :: ltimedeplw   = .false.  !< Switch for longwave radiative fluxes
  logical       :: ltimedepsw   = .false.  !< Switch for shortwave radiative fluxes

  integer    :: ntimedepsurf
  integer    :: ntimedepnudge
  integer    :: ntimedeplw
  integer    :: ntimedepsw

  real, allocatable     :: timeflux (:)
  real, allocatable     :: bctfxmt (:)
  real, allocatable     :: bctfxpt (:)
  real, allocatable     :: bctfymt (:)
  real, allocatable     :: bctfypt (:)
  real, allocatable     :: bctfzt (:)
  !real, allocatable     :: bctfzft (:)

  real, allocatable     :: timenudge (:)
  real, allocatable     :: thlproft (:,:)
  real, allocatable     :: qtproft (:,:)
  real, allocatable     :: uproft (:,:)
  real, allocatable     :: vproft (:,:)

  real, allocatable     :: timelw (:)
  real, allocatable     :: skyLWt (:)
  real, allocatable     :: timesw (:)
  real, allocatable     :: netswt (:,:)


contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine inittimedep
    use modmpi,    only :myid,my_real,mpi_logical,mpierr,comm3d
    use modglobal, only :cexpnr,kb,ke,kh,kmax,ifinput,runtime,zf,skyLW,nfcts
    use modibmdata, only : bctfxm, bctfxp, bctfym, bctfyp, bctfz!, bctfzf
    use modfields, only: thlprof
    !use initfac, only : netsw !Should probably be moved to somewhere else

    implicit none

    character (80):: chmess
    character (1) :: chmess1
    integer :: k,t,n,ierr
    real :: dummyr
    real, allocatable, dimension (:) :: height

    ltimedep = (ltimedepsurf .or. ltimedepnudge) .or. (ltimedeplw .or. ltimedepsw)

    if (.not. ltimedep) return

    if (ltimedepsurf) then
      allocate(timeflux (1:ntimedepsurf))
      allocate(bctfxmt (1:ntimedepsurf))
      allocate(bctfxpt (1:ntimedepsurf))
      allocate(bctfymt (1:ntimedepsurf))
      allocate(bctfypt (1:ntimedepsurf))
      allocate(bctfzt (1:ntimedepsurf))
      !allocate(bctfzft (1:ntimedepsurf))

      timeflux = 0.
      bctfxmt = bctfxm
      bctfxpt = bctfxp
      bctfymt = bctfym
      bctfypt = bctfyp
      bctfzt = bctfz
      !bctfzft = bctfzf

      if (myid==0) then
        open(ifinput,file='timedepsurf.inp.'//cexpnr)
        read(ifinput,'(a80)') chmess
        !write(6,*) chmess
        read(ifinput,'(a80)') chmess
        !write(6,*) chmess

        !--- load fluxes---
        !t    = 1
        ierr = 0
        do t = 1,ntimedepsurf
          read(ifinput,*, iostat = ierr) timeflux(t), bctfxmt(t), bctfxpt(t), bctfymt(t), bctfypt(t), bctfzt(t)!, bctfzft(t)
          !write(*,*) t, timeflux(t), bctfxmt(t), bctfxpt(t), bctfymt(t), bctfypt(t), bctfzt(t)!, bctfzft(t)
          !if (ierr < 0) then
            !stop 'STOP: No time dependend data for end of run (surface fluxes)'
          !end if
        end do
        !if(timeflux(1)>runtime) then
         !write(6,*) 'Time dependent surface variables do not change before end of simulation'
         !ltimedepsurf=.false.
        !endif
        ! flush to the end of fluxlist
        !do while (ierr ==0)
          !read (ifinput,*,iostat=ierr) dummyr
        !end do
        !backspace (ifinput)
        !close(ifinput)

      end if !myid==0

      call MPI_BCAST(timeflux, ntimedepsurf, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctfxmt , ntimedepsurf, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctfxpt , ntimedepsurf, MY_REAL, 0, comm3d,mpierr)
      call MPI_BCAST(bctfymt , ntimedepsurf, MY_REAL, 0, comm3d,mpierr)
      call MPI_BCAST(bctfypt , ntimedepsurf, MY_REAL, 0, comm3d,mpierr)
      call MPI_BCAST(bctfzt  , ntimedepsurf, MY_REAL, 0, comm3d,mpierr)
      !call MPI_BCAST(bctfzft , ntimedepsurf, MY_REAL, 0, comm3d,mpierr)

    end if

    if (ltimedepnudge) then

      allocate(timenudge (1:ntimedepnudge))
      allocate(height    (kb:ke+kh))
      allocate(uproft    (kb:ke+kh,ntimedepnudge))
      allocate(vproft    (kb:ke+kh,ntimedepnudge))
      allocate(thlproft  (kb:ke+kh,ntimedepnudge))
      allocate(qtproft   (kb:ke+kh,ntimedepnudge))

      timenudge = 0
      thlproft = 0
      qtproft = 0
      uproft = 0
      vproft = 0

      if (myid == 0) then
        !---load nudging profiles----
        open(ifinput,file='timedepnudge.inp.'//cexpnr)
        read(ifinput,'(a80)') chmess
        !write(6,*) chmess

        !t = 0
        do t = 1,ntimedepnudge
          !t = t + 1
          chmess1 = "#"
          ierr = 1 ! not zero
          do while (.not.(chmess1 == "#" .and. ierr == 0)) !search for the next line consisting of "# time", from there onwards the profiles will be read
            read(ifinput,*,iostat=ierr) chmess1, timenudge(t)
            !if (ierr < 0) then
              !stop 'STOP: No time dependend data for end of run'
            !end if
          end do

          !write (*,*) 'timenudge = ',timenudge(t)
          !write(*,*) 'Nudging profiles'
          do k=kb,ke
            read (ifinput,*) &
            height   (k)  , &
            thlproft (k,t), &
            qtproft  (k,t), &
            uproft   (k,t), &
            vproft   (k,t)

            !write(*,*) height(k), thlproft (k,t), qtproft(k,t), uproft(k,t), vproft(k,t)
          end do
        end do
      end if !myid == 0

      call MPI_BCAST(timenudge, ntimedepnudge, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(thlproft, (kmax+1)*ntimedepnudge,MY_REAL,0,comm3d,mpierr)
      call MPI_BCAST(qtproft,  (kmax+1)*ntimedepnudge,MY_REAL,0,comm3d,mpierr)
      call MPI_BCAST(uproft,   (kmax+1)*ntimedepnudge,MY_REAL,0,comm3d,mpierr)
      call MPI_BCAST(vproft,   (kmax+1)*ntimedepnudge,MY_REAL,0,comm3d,mpierr)

      deallocate(height)

    end if

    if (ltimedeplw) then
      allocate(timelw (1:ntimedeplw))
      allocate(skyLWt (1:ntimedeplw))

      ! Read longwave
      timelw = 0.
      skyLWt = skyLW

      if (myid==0) then
       open(ifinput,file='timedeplw.inp.'//cexpnr)
       read(ifinput,'(a80)') chmess
       !write(6,*) chmess
       read(ifinput,'(a80)') chmess
       !write(6,*) chmess

       !--- load fluxes---
       !t    = 1
       ierr = 0
       do t = 1,ntimedeplw
          read(ifinput,*, iostat = ierr) timelw(t), skyLWt(t)
          !write(*,*) t, timelw(t), skyLWt(t)
          !if (ierr < 0) then
            !stop 'STOP: No time dependend data for end of run (surface fluxes)'
          !end if
       end do
       !if(timelw(1)>runtime) then
         !write(6,*) 'Time dependent surface variables do not change before end of simulation'
         !ltimedeplw=.false.
       !endif
       ! flush to the end of fluxlist
       !do while (ierr ==0)
          !read (ifinput,*,iostat=ierr) dummyr
       !end do
       !backspace (ifinput)
       !close(ifinput)

      end if ! myid = 0

      call MPI_BCAST(timelw, ntimedeplw, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(skyLWt, ntimedeplw, MY_REAL, 0, comm3d, mpierr)

    end if !ltimedeplw

    if (ltimedepsw) then
      allocate(timesw(1:ntimedepsw))
      allocate(netswt(1:nfcts, 1:ntimedepsw))

      timesw = 0.
      netswt = 0.

      if (myid == 0) then
        ! Read shortwave
        open (ifinput,file='timedepsw.inp.'//cexpnr)
        read (ifinput,'(a80)') chmess ! first line is a description of the file
        read (ifinput, *) (timesw(t), t=1,ntimedepsw) ! second line is the times

        do n = 1,nfcts
         read (ifinput, *) (netswt(n,t), t=1,ntimedepsw)
        end do

        !write(*,*) "read timedepsw"

      end if !myid==0

      call MPI_BCAST(timesw, ntimedepsw, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(netswt, ntimedepsw*nfcts, MY_REAL, 0, comm3d, mpierr)

    end if ! timedepsw

    call MPI_BCAST(ltimedepsurf ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(ltimedepnudge,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(ltimedeplw,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(ltimedepsw,1,MPI_LOGICAL,0,comm3d,mpierr)

    call timedep

  end subroutine inittimedep

  subroutine timedep

!-----------------------------------------------------------------|
!                                                                 |
!*** *timedep*  calculates ls forcings and surface forcings       |
!               case as a funtion of timee                        |
!                                                                 |
!      Roel Neggers    K.N.M.I.     01/05/2001                    |
!                                                                 |
!                                                                 |
!    calls                                                        |
!    * timedepz                                                   |
!      calculation of large scale advection, radiation and        |
!      surface fluxes by interpolation between prescribed         |
!      values at certain times                                    |
!                                                                 |
!    * timedepsurf                                                |
!      calculation  surface fluxes by interpolation               |
!      between prescribed values at certain times                 |
!                                                                 |
!                                                                 |
!-----------------------------------------------------------------|
    implicit none

    if (.not. ltimedep) return

    call timedepsurf
    call timedepnudge
    call timedeplw
    call timedepsw

  end subroutine timedep

  subroutine timedepsurf
    use modmpi,     only : myid
    use modglobal,  only : timee
    use modibmdata, only : bctfxm, bctfxp, bctfym, bctfyp, bctfz!, bctfzf
    implicit none
    integer t
    real fac

    if(.not.(ltimedepsurf)) return

    !     --- interpolate! ----
    do t=ntimedepsurf,1,-1
      if (timee .ge. timeflux(t)) then
        exit
      endif
    end do

    ! if ((myid == 0) .or. (myid == 1)) then
    !   write(*, *) "myid", myid, "t", t, "timee", timee, "timeflux(t)", timeflux(t)
    ! end if

    if (t .ne. ntimedepsurf) then
      fac = (timee-timeflux(t)) / (timeflux(t+1)-timeflux(t))
      bctfxm = bctfxmt(t) + fac * (bctfxmt(t+1) - bctfxmt(t))
      bctfxp = bctfxpt(t) + fac * (bctfxpt(t+1) - bctfxpt(t))
      bctfym = bctfymt(t) + fac * (bctfymt(t+1) - bctfymt(t))
      bctfyp = bctfypt(t) + fac * (bctfypt(t+1) - bctfypt(t))
      bctfz = bctfzt(t) + fac * (bctfzt(t+1) - bctfzt(t))
      !bctfzf = bctfzft(t) + fac * (bctfzft(t+1) - bctfzft(t))
   end if

    ! if ((myid == 0) .or. (myid == 1)) then
    !   write(*, *) "myid", myid, "bctfz", bctfz, "bctfzf", bctfzf
    ! end if

    return
  end subroutine timedepsurf

  subroutine timedepnudge
    use modfields,   only : thlprof, qtprof, uprof, vprof
    use modglobal,   only : timee,dzf,dzh,kb,ke,kh,kmax
    use modmpi,      only : myid

    implicit none
    integer t,k
    real fac

    if(.not.(ltimedepnudge)) return

    !---- interpolate ----
    do t=ntimedepnudge,1,-1
      if (timee .ge. timenudge(t)) then
        exit
      endif
    end do

    if (t .ne. ntimedepnudge) then
      fac = (timee - timenudge(t)) / (timenudge(t+1) - timenudge(t))
      thlprof = thlproft(:,t) + fac * (thlproft(:,t+1) - thlproft(:,t))
      qtprof  = qtproft (:,t) + fac * (qtproft (:,t+1) - qtproft (:,t))
      uprof   = uproft  (:,t) + fac * (uproft  (:,t+1) - uproft  (:,t))
      vprof   = vproft  (:,t) + fac * (vproft  (:,t+1) - vproft  (:,t))
    end if

    ! if ((myid == 0) .or. (myid == 1)) then
    !   write(*, *) "myid, t, timee, timenudge(t), thlproft(ke,t), thlproft(ke,t+1), thlprof(ke)"
    !   write(*, *) myid, t, timee, timenudge(t), thlproft(ke,t), thlproft(ke,t+1), thlprof(ke)
    ! end if
    !write(*, *) "myid, thlproft(k,t), thlprof(k)"
    !do k = kb,ke
        !write(*,*) myid, thlproft(k, t), thlprof(k)
    !end do

  return

  end subroutine timedepnudge

  subroutine timedeplw
    use modglobal,    only : timee, skyLW, rk3step, tnextEB
    use modmpi,       only : myid

    implicit none
    integer t,k
    real fac

    if(.not.(ltimedeplw)) return

    if ((rk3step .eq. 3) .and. (timee .ge. tnextEB)) then

    ! if (myid == 0) then
    !  write(*,*) "EB coming up so changing longwave forcing"
    ! end if

    !---- interpolate ----
    do t=ntimedeplw,1,-1
      if (timee .ge. timelw(t)) then
        exit
      endif
    end do

    if (t .ne. ntimedeplw) then
      fac = (timee - timelw(t)) / (timelw(t+1) - timelw(t))
      skyLW = skyLWt(t) + fac * (skyLWt(t+1) - skyLWt(t))
    end if

    end if

  end subroutine timedeplw

  subroutine timedepsw
    use modglobal, only : timee, nfcts, rk3step, tnextEB
    use initfac, only : netsw
    use modmpi, only : myid

    implicit none
    integer t,n
    real fac

    if(.not.(ltimedepsw .and. myid==0)) return

    if ((rk3step .eq. 3) .and. (timee .ge. tnextEB)) then

    ! if (myid == 0) then
    !  write(*,*) "EB coming up so changing solar position"
    ! end if

    !---- interpolate ----
    do t=ntimedepsw,1,-1
      if (timee .ge. timesw(t)) then
        exit
      endif
    end do

    if (t .ne. ntimedepsw) then
      fac = (timee - timesw(t)) / (timesw(t+1) - timesw(t))
      do n=1,nfcts
         netsw(n) = netswt(n,t) + fac * (netswt(n,t+1) - netswt(n,t))
      end do
    end if

    end if

  end subroutine timedepsw

  subroutine exittimedep

    implicit none

    if (.not. ltimedep) return

    if (ltimedepsurf) then
      deallocate(timeflux, bctfxmt, bctfxpt, bctfymt, bctfypt, bctfzt)!, bctfzft)
    end if

    if (ltimedepnudge) then
      deallocate(timenudge, thlproft, qtproft, uproft, vproft)
    end if

  end subroutine

end module modtimedep
