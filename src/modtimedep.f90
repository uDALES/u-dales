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
#if defined(_GPU)
use cudafor
#endif

implicit none
private
public :: inittimedep, timedep, timedep_step, ltimedep, ltimedepsurf, ltimedepnudge, ltimedeplw, ltimedepsw, &
          ntimedepsurf, ntimedepnudge, ntimedeplw, ntimedepsw, exittimedep
! The interpolation tables and the two routines that read them, for
! tests.f90::tests_timedep. A runmode test has no input file to read, so the
! tables are what it has to stage; the interpolation is what it is testing.
! Public unconditionally rather than under UDALES_DEBUG because tests.f90 is
! compiled into the Release build too.
public :: timedepsurf, timedepnudge
public :: timeflux, bctfxmt, bctfxpt, bctfymt, bctfypt, bctfzt
public :: timenudge, thlproft, qtproft, uproft, vproft
#if defined(_GPU) && defined(UDALES_DEBUG)
public :: timedepnudge_device
public :: thlproft_d, qtproft_d, uproft_d, vproft_d
#endif

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
#if defined(_GPU)
  ! Pinned because they are handed to the device, once, from inittimedep. The
  ! bracket search stays on the host and reads timenudge there, so that one is
  ! not mirrored and not pinned.
  real, allocatable, pinned :: thlproft (:,:)
  real, allocatable, pinned :: qtproft (:,:)
  real, allocatable, pinned :: uproft (:,:)
  real, allocatable, pinned :: vproft (:,:)
  real, device, allocatable :: thlproft_d (:,:)
  real, device, allocatable :: qtproft_d (:,:)
  real, device, allocatable :: uproft_d (:,:)
  real, device, allocatable :: vproft_d (:,:)
#else
  real, allocatable     :: thlproft (:,:)
  real, allocatable     :: qtproft (:,:)
  real, allocatable     :: uproft (:,:)
  real, allocatable     :: vproft (:,:)
#endif

  real, allocatable     :: timelw (:)
  real, allocatable     :: skyLWt (:)
  real, allocatable     :: timesw (:)
  real, allocatable     :: netswt (:,:)


contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine inittimedep
    use modmpi,    only :myid,my_real,mpi_logical,mpierr,comm3d
    use modglobal, only :cexpnr,kb,ke,kh,kmax,ifinput,skyLW,nfcts
    use modibmdata, only : bctfxm, bctfxp, bctfym, bctfyp, bctfz!, bctfzf
    !use initfac, only : netsw !Should probably be moved to somewhere else

    implicit none

    character (80):: chmess
    character (1) :: chmess1
    integer :: k,t,n,ierr
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

#if defined(_GPU)
      ! The tables go up once and stay. Nothing writes them after this point -
      ! they are the file, broadcast - so the per-stage interpolation reads
      ! them where it runs and the loop moves nothing across the bus for the
      ! nudging profiles at all. Four columns of kmax+1 per table entry, which
      ! at the sizes these files are written for is tens of kilobytes.
      !
      ! Allocated here rather than from initCUDA because they are this
      ! module's state and initCUDA has not run yet - the device context is
      ! created on first use, exactly as inittstep's geometry cache relies on.
      allocate(thlproft_d(kb:ke+kh,ntimedepnudge))
      allocate(qtproft_d (kb:ke+kh,ntimedepnudge))
      allocate(uproft_d  (kb:ke+kh,ntimedepnudge))
      allocate(vproft_d  (kb:ke+kh,ntimedepnudge))
      thlproft_d = thlproft
      qtproft_d  = qtproft
      uproft_d   = uproft
      vproft_d   = vproft
#endif

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

    ! The host one, deliberately, in a GPU build too: initCUDA has not run,
    ! so uprof_d and the rest do not exist yet, and what this call is for is
    ! to leave the host profiles holding the t = 0 interpolation for initCUDA
    ! to hand over. From the first stage onwards timedep_step takes over and
    ! writes the device copies instead.
    call timedep

  end subroutine inittimedep

  !> The one the time loop calls.
  subroutine timedep_step
    implicit none
#if defined(_GPU)
    call timedep_device
#else
    call timedep
#endif
  end subroutine timedep_step

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
    use modglobal,   only : timee

    implicit none
    integer t
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

#if defined(_GPU)
  !> timedep for a run whose profiles live on the device.
  !!
  !! Three of the four parts are the host ones, unchanged, and that is not an
  !! omission.
  !!
  !! timedepsurf writes five scalars in modibmdata. The wall-function kernel
  !! that reads them is an OpenACC region in a host routine, so the scalars
  !! reach it as implicit firstprivate, captured at launch - the value a stage
  !! sees is the one this routine wrote at the top of that stage, with nothing
  !! to copy and no second copy to keep in step. The compiler says so: the
  !! wall-function kernel's -Minfo line reads
  !! "Generating implicit firstprivate(...,bctfxm,...,bctfxp,bctfyp,...,bctfz,...)".
  !! Deleting the call here would not break a copy; it would stop the flux
  !! changing with time, and runmode 1018 is what covers the interpolation -
  !! no case in the repository combines iwalltemp = 1 with libm, so the wall
  !! function never reads a moving flux in any test on either device.
  !!
  !! timedeplw and timedepsw feed the energy balance. That runs on the host,
  !! reads skyLW and netsw there, and both routines already restrict
  !! themselves to the stage the balance fires on.
  !!
  !! Only the nudging profiles have a device copy, so only they differ.
  subroutine timedep_device
    implicit none

    if (.not. ltimedep) return

    call timedepsurf
    call timedepnudge_device
    call timedeplw
    call timedepsw

  end subroutine timedep_device

  !> timedepnudge writing the device profiles.
  !!
  !! The split follows the shape of the host routine rather than fighting it.
  !! The bracket search is a walk over ntimedepnudge scalars against timee and
  !! stays here, on the host, where timee lives; t and fac then reach the
  !! kernels as implicit firstprivate, so no clause names them - a private
  !! clause on either would hand every thread an uninitialised copy.
  !!
  !! The arithmetic is the host expression, operand for operand, so the only
  !! room left between the two is fused-multiply-add contraction. Nothing is
  !! precomputed out of it: the obvious candidates are the table differences,
  !! which cost as much memory again to save four subtractions per level, and
  !! the reciprocal of the bracket width, which would turn a division into a
  !! multiplication and move the profile - and the profile is what the nudging
  !! tendency and the inlet are built from.
  !!
  !! Three kernels rather than one because thlprof_d and qtprof_d only exist
  !! under ltempeq and lmoist. Their launches are what this routine costs -
  !! the work itself is one column - and the trade is against the four
  !! unconditional column uploads updateDevice used to do on every stage of
  !! every run, ltimedep or not, which are now gone.
  subroutine timedepnudge_device
    use modglobal, only : timee, kb, ke, kh, ltempeq, lmoist
    use modcuda,   only : uprof_d, vprof_d, thlprof_d, qtprof_d

    implicit none
    integer t, k
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

      !$acc parallel loop default(present)
      do k = kb, ke+kh
        uprof_d(k) = uproft_d(k,t) + fac * (uproft_d(k,t+1) - uproft_d(k,t))
        vprof_d(k) = vproft_d(k,t) + fac * (vproft_d(k,t+1) - vproft_d(k,t))
      end do
      !$acc end parallel loop

      if (ltempeq) then
        !$acc parallel loop default(present)
        do k = kb, ke+kh
          thlprof_d(k) = thlproft_d(k,t) + fac * (thlproft_d(k,t+1) - thlproft_d(k,t))
        end do
        !$acc end parallel loop
      end if

      if (lmoist) then
        !$acc parallel loop default(present)
        do k = kb, ke+kh
          qtprof_d(k) = qtproft_d(k,t) + fac * (qtproft_d(k,t+1) - qtproft_d(k,t))
        end do
        !$acc end parallel loop
      end if
    end if

  return

  end subroutine timedepnudge_device
#endif

  subroutine timedeplw
    use modglobal,    only : timee, skyLW, rk3step, tnextEB

    implicit none
    integer t
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
#if defined(_GPU)
      deallocate(thlproft_d, qtproft_d, uproft_d, vproft_d)
#endif
    end if

  end subroutine

end module modtimedep
