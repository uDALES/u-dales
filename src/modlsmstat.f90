!> \file modlsmstat.f90
!!  Calculates the land surface model statistics


!>
!!  Calculates the land surface model statistics
!>
!! Profiles of the LSM statistics. Written to lsmstat.expnr
!! If netcdf is true, this module also writes in the profiles.expnr.nc output
!!  \author Stephan de Roode, TU Delft
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
module modlsmstat

  use modglobal, only : longint

implicit none
!private
PUBLIC :: initlsmstat, lsmstat, exitlsmstat
save
!NetCDF variables
  integer,parameter :: nvar = 5
  character(80),dimension(nvar,4) :: ncname

  real    :: dtav, timeav
  integer(kind=longint) :: idtav,itimeav,tnext,tnextwrite
  integer :: nsamples
  logical :: lstat= .false. !< switch to enable the lsmiative statistics (on/off)

!     ------

!   --------------
  real, allocatable :: gammasav(:)
  real, allocatable :: phiwav(:)
  real, allocatable :: tsoilav(:)
  real, allocatable :: lambdaav(:)
  real, allocatable :: lambdasav(:)

!
  real, allocatable :: gammasmn(:)
  real, allocatable :: phiwmn(:)
  real, allocatable :: tsoilmn(:)
  real, allocatable :: lambdamn(:)
  real, allocatable :: lambdasmn(:)

contains
!> Initialization routine, reads namelists and inits variables
  subroutine initlsmstat
    use modmpi,    only : myid,mpierr, comm3d,my_real, mpi_logical
    use modglobal, only : dtmax, ifnamopt,fname_options, ifoutput, cexpnr,dtav_glob,timeav_glob,ladaptive,dt_lim,btime,tres
    use modstat_nc, only : lnetcdf, redefine_nc,define_nc,ncinfo
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav,ncid_prof=>ncid
    use modsurfdata,only : ksoilmax,isurf
    implicit none

    integer ierr
    namelist/NAMLSMSTAT/ &
    dtav,timeav,lstat

    dtav=dtav_glob;timeav=timeav_glob
    lstat = .false.
    if (isurf /=1) return
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMLSMSTAT,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMLSMSTAT'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMLSMSTAT'
      endif
      write(6 ,NAMLSMSTAT)
      close(ifnamopt)
    end if

    call MPI_BCAST(timeav     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dtav       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(lstat   ,1,MPI_LOGICAL,0,comm3d,mpierr)
    idtav = dtav/tres
    itimeav = timeav/tres

    tnext      = idtav   +btime
    tnextwrite = itimeav +btime
    nsamples = itimeav/idtav


    if(.not.(lstat)) return
    dt_lim = min(dt_lim,tnext)

    if (abs(timeav/dtav-nsamples)>1e-4) then
      stop 'timeav must be a integer multiple of dtav'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if

    allocate(phiwav(ksoilmax))
    allocate(tsoilav(ksoilmax))
    allocate(lambdaav(ksoilmax))
    allocate(lambdasav(ksoilmax))
    allocate(gammasav(ksoilmax))

    allocate(phiwmn(ksoilmax))
    allocate(tsoilmn(ksoilmax))
    allocate(lambdamn(ksoilmax))
    allocate(lambdasmn(ksoilmax))
    allocate(gammasmn(ksoilmax))

    phiwmn = 0.0
    tsoilmn = 0.0
    lambdamn = 0.0
    lambdasmn = 0.0
    gammasmn = 0.0

    if(myid==0)then
      open (ifoutput,file='lsmstat.'//cexpnr,status='replace')
      close (ifoutput)
    end if
    if (lnetcdf) then
      idtav = idtav_prof
      itimeav = itimeav_prof
      tnext      = idtav+btime
      tnextwrite = itimeav+btime
      nsamples = itimeav/idtav

      if (myid==0) then
        call ncinfo(ncname( 1,:),'tsoil','Soil temperature','W/m^2','tts')
        call ncinfo(ncname( 2,:),'phiw','Soil moisture content','W/m^2','tts')
        call ncinfo(ncname( 3,:),'lambda','Heat conductivity soil layer','W/m/K','tts')
        call ncinfo(ncname( 4,:),'lambdas','Soil moisture diffusivity soil layer','m^2/s','tts')
        call ncinfo(ncname( 5,:),'gammas','Soil moisture conductivity soil layer','M/s','tts')

        call redefine_nc(ncid_prof)
        call define_nc( ncid_prof, NVar, ncname)
      end if

   end if

  end subroutine initlsmstat
!> General routine, does the timekeeping
  subroutine lsmstat
    use modglobal, only : rk3step,timee,dt_lim
    implicit none
    if (.not. lstat) return
    if (rk3step/=3) return
    if(timee<tnext .and. timee<tnextwrite) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if
    if (timee>=tnext) then
      tnext = tnext+idtav
      call do_lsmstat
    end if
    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+itimeav
      call writelsmstat
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))

  end subroutine lsmstat

!> Calculates the statistics
  subroutine do_lsmstat

    use modmpi,    only :  slabsum
    use modglobal, only : rslabs,dzf,i1,j1,i2,j2
    use modsurfdata, only : ksoilmax,tsoil,phiw,lambda,lambdas,gammas

    implicit none
    integer :: k

    tsoilav  = 0.
    phiwav  = 0.
    lambdaav  = 0.
    lambdasav  = 0.
    gammasav = 0.

    call slabsum(tsoilav ,1,ksoilmax,tsoil ,1,i2,1,j2,1,ksoilmax,2,i1,2,j1,1,ksoilmax)
    call slabsum(phiwav ,1,ksoilmax,phiw ,1,i2,1,j2,1,ksoilmax,2,i1,2,j1,1,ksoilmax)
    call slabsum(lambdaav ,1,ksoilmax,lambda ,1,i2,1,j2,1,ksoilmax,2,i1,2,j1,1,ksoilmax)
    call slabsum(lambdasav ,1,ksoilmax,lambdas ,1,i2,1,j2,1,ksoilmax,2,i1,2,j1,1,ksoilmax)
    call slabsum(gammasav ,1,ksoilmax,gammas ,1,i2,1,j2,1,ksoilmax,2,i1,2,j1,1,ksoilmax)
 !    ADD SLAB AVERAGES TO TIME MEAN

    phiwmn = phiwmn + phiwav/rslabs
    tsoilmn = tsoilmn + tsoilav/rslabs
    lambdamn = lambdamn + lambdaav/rslabs
    lambdasmn = lambdasmn + lambdasav/rslabs
    gammasmn = gammasmn + gammasav/rslabs

  end subroutine do_lsmstat


!> Write the statistics to file
  subroutine writelsmstat
      use modmpi,    only : myid
      use modglobal, only : cexpnr,ifoutput,zf,zh,rtimee
      use modstat_nc, only: lnetcdf, writestat_nc
      use modgenstat, only: ncid_prof=>ncid,nrec_prof=>nrec
      use modsurfdata,only : ksoilmax,zsoil
      implicit none
      real,dimension(ksoilmax,nvar) :: vars
      integer nsecs, nhrs, nminut,k


      nsecs   = nint(rtimee)
      nhrs    = int(nsecs/3600)
      nminut  = int(nsecs/60)-nhrs*60
      nsecs   = mod(nsecs,60)

      phiwmn   = phiwmn    /nsamples
      tsoilmn   = tsoilmn    /nsamples
      lambdamn   = lambdamn    /nsamples
      lambdasmn   = lambdasmn    /nsamples
      gammasmn   = gammasmn   /nsamples
  !     ----------------------
  !     2.0  write the fields
  !           ----------------

    if(myid==0)then
      open (ifoutput,file='lsmstat.'//cexpnr,position='append')
      write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
      '#--------------------------------------------------------'      &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '      &
      ,nhrs,':',nminut,':',nsecs      &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
      write (ifoutput,'(A/A/A)') &
          '#--------------------------------------------------------------------------' &
          ,'#LEV HEIGHT  T_SOIL      SOIL MOIST   HEAT COND.   MOIST DIFF.  MOIST COND.' &
          ,'#    (M)    (K)          (M^3/M^3)    (W/M/K)      (M^2/S)       (M/S)      '
      do k=1,ksoilmax
        write(ifoutput,'(I3,F8.2,F10.4,4E13.4)') &
            k, zsoil(k),&
            tsoilmn(k),&
            phiwmn(k),&
            lambdamn(k),&
            lambdasmn(k),&
            gammasmn(k)
      end do
      close (ifoutput)
      if (lnetcdf) then
        vars(:, 1) = tsoilmn
        vars(:, 2) = phiwmn
        vars(:, 3) = lambdamn
        vars(:, 4) = lambdasmn
        vars(:, 5) = gammasmn
       call writestat_nc(ncid_prof,nvar,ncname,vars(1:ksoilmax,:),nrec_prof,ksoilmax)
      end if
    end if ! end if(myid==0)

    phiwmn = 0.0
    tsoilmn = 0.0
    lambdamn = 0.0
    lambdasmn = 0.0
    gammasmn  = 0.0


  end subroutine writelsmstat

!> Cleans up after the run
  subroutine exitlsmstat
    implicit none

    !deallocate variables that are needed in modlsmiation

    if(.not.(lstat)) return

    deallocate(phiwav,tsoilav,lambdaav,lambdasav,gammasav)
    deallocate(phiwmn,tsoilmn,lambdamn,lambdasmn,gammasmn)



  end subroutine exitlsmstat


end module modlsmstat
