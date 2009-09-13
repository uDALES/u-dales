!> \file modradstat.f90
!!  Calculates the radiative statistics


!>
!!  Calculates the radiative statistics
!>
!! Profiles of the radiative statistics. Written to radstat.expnr
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
module modradstat


implicit none
!private
PUBLIC :: initradstat, radstat, exitradstat
save
!NetCDF variables
  integer,parameter :: nvar = 7
  character(80),dimension(nvar,4) :: ncname

  real    :: dtav, timeav,tnext,tnextwrite
  integer :: nsamples
  logical :: lstat= .false. !< switch to enable the radiative statistics (on/off)

!     ------

  real, allocatable :: tllwtendavl(:)
  real, allocatable :: tlswtendavl(:)
  real, allocatable :: lwuavl(:)
  real, allocatable :: lwdavl(:)
  real, allocatable :: swdavl(:)
  real, allocatable :: swuavl(:)

!   --------------
  real, allocatable :: tllwtendav(:)
  real, allocatable :: tlswtendav(:)
  real, allocatable :: lwuav(:)
  real, allocatable :: lwdav(:)
  real, allocatable :: swdav(:)
  real, allocatable :: swuav(:)

!
  real, allocatable :: tllwtendmn(:)
  real, allocatable :: tlswtendmn(:)
  real, allocatable :: lwumn(:)
  real, allocatable :: lwdmn(:)
  real, allocatable :: swdmn(:)
  real, allocatable :: swumn(:)
  real, allocatable :: tlradlsmn(:)

contains
!> Initialization routine, reads namelists and inits variables
  subroutine initradstat
    use modmpi,    only : myid,mpierr, comm3d,my_real, mpi_logical
    use modglobal, only : dtmax, k1,kmax, ifnamopt,fname_options, ifoutput, cexpnr,dtav_glob,timeav_glob,ladaptive,dt_lim,btime
    use modstat_nc, only : lnetcdf, redefine_nc,define_nc,ncinfo
    use modgenstat, only : dtav_prof=>dtav, timeav_prof=>timeav,ncid_prof=>ncid

    implicit none

    integer ierr
    namelist/NAMRADSTAT/ &
    dtav,timeav,lstat

    dtav=dtav_glob;timeav=timeav_glob
    lstat = .false.

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMRADSTAT,iostat=ierr)
      write(6 ,NAMRADSTAT)
      close(ifnamopt)
    end if

    call MPI_BCAST(timeav     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dtav       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(lstat   ,1,MPI_LOGICAL,0,comm3d,mpierr)

    tnext      = dtav-1e-3+btime
    tnextwrite = timeav-1e-3+btime
    nsamples = nint(timeav/dtav)

   !allocate variables that are needed in modradiation
    allocate(lwdavl(k1))
    allocate(lwuavl(k1))
    allocate(swdavl(k1))
    allocate(swuavl(k1))
    allocate (tllwtendavl(k1))
    allocate (tlswtendavl(k1))

    if(.not.(lstat)) return
    dt_lim = min(dt_lim,tnext)

    if (abs(timeav/dtav-nsamples)>1e-4) then
      stop 'timeav must be a integer multiple of dtav'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if

    allocate(lwuav(k1))
    allocate(lwdav(k1))
    allocate(swdav(k1))
    allocate(swuav(k1))
    allocate(tllwtendav(k1))
    allocate(tlswtendav(k1))

    allocate(lwumn(k1))
    allocate(lwdmn(k1))
    allocate(swdmn(k1))
    allocate(swumn(k1))
    allocate(tllwtendmn(k1))
    allocate(tlswtendmn(k1))
    allocate(tlradlsmn(k1))

    lwumn = 0.0
    lwdmn = 0.0
    swdmn = 0.0
    swumn = 0.0
    tllwtendmn = 0.0
    tlswtendmn = 0.0
    tlradlsmn  = 0.0

    if(myid==0)then
      open (ifoutput,file='radstat.'//cexpnr,status='replace')
      close (ifoutput)
    end if
    if (lnetcdf) then
      dtav = dtav_prof
      timeav = timeav_prof

      tnext      = dtav-1e-3+btime
      tnextwrite = timeav-1e-3+btime
      nsamples = nint(timeav/dtav)

      if (myid==0) then
        call ncinfo(ncname( 1,:),'tllwtend','Long wave radiative tendency','K/s','tt')
        call ncinfo(ncname( 2,:),'tlswtend','Short wave radiative tendency','K/s','tt')
        call ncinfo(ncname( 3,:),'tlradls','Large scale radiative tendency','K/s','tt')
        call ncinfo(ncname( 4,:),'lwu','Long wave upward radiative flux','W/m^2','mt')
        call ncinfo(ncname( 5,:),'lwd','Long wave downward radiative flux','W/m^2','mt')
        call ncinfo(ncname( 6,:),'swd','Short wave downward radiative flux','W/m^2','mt')

        call redefine_nc(ncid_prof)
        call define_nc( ncid_prof, NVar, ncname)
      end if

   end if

  end subroutine initradstat
!> General routine, does the timekeeping
  subroutine radstat
    use modglobal, only : rk3step,timee,dt_lim
    implicit none
    if (.not. lstat) return
    if (rk3step/=3) return
    if(timee<tnext .and. timee<tnextwrite) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if
    if (timee>=tnext) then
      tnext = tnext+dtav
      call do_radstat
    end if
    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+timeav
      call writeradstat
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))

  end subroutine radstat

!> Calculates the statistics
  subroutine do_radstat

    use modmpi,    only : nprocs,comm3d,nprocs,my_real, mpi_sum,mpierr, slabsum
    use modglobal, only : kmax,rslabs,cp,dzf,i1,j1,k1
    use modfields, only : thlpcar
    use modraddata, only : lwd,lwu,swd,swu, rho_air_mn

    implicit none
    integer :: k

    lwdav  = 0.
    lwuav  = 0.
    swdav  = 0.
    tllwtendav = 0.
    tlswtendav = 0.
    lwdavl  = 0.
    lwuavl  = 0.
    swdavl  = 0.
    swuavl  = 0.
    tllwtendavl = 0.
    tlswtendavl = 0.

    do k=1,k1
      lwdavl(k) = sum(lwd(2:i1,2:j1,k))
      lwuavl(k) = sum(lwu(2:i1,2:j1,k))
      swdavl(k) = sum(swd(2:i1,2:j1,k))
      swuavl(k) = sum(swu(2:i1,2:j1,k))
    end do

    do k=1,kmax
      tllwtendavl(k) = -(lwdavl(k+1)-lwdavl(k)+lwuavl(k+1)-lwuavl(k))/(rho_air_mn*cp*dzf(k))
      tlswtendavl(k) =  (swdavl(k+1)-swdavl(k))/(rho_air_mn*cp*dzf(k))
    end do
    !swdavl = swdav

    call MPI_ALLREDUCE(lwdavl, lwdav, kmax,    MY_REAL, &
                         MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(lwuavl, lwuav, kmax,    MY_REAL, &
                         MPI_SUM, comm3d,mpierr)
   call MPI_ALLREDUCE(swdavl, swdav, kmax,    MY_REAL, &
                         MPI_SUM, comm3d,mpierr)
   call MPI_ALLREDUCE(swuavl, swuav, kmax,    MY_REAL, &
                         MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(tllwtendavl, tllwtendav, kmax,    MY_REAL, &
                         MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(tlswtendavl, tlswtendav, kmax,    MY_REAL, &
                         MPI_SUM, comm3d,mpierr)

 !    ADD SLAB AVERAGES TO TIME MEAN

    lwumn = lwumn + lwuav/rslabs
    lwdmn = lwdmn + lwdav/rslabs
    swdmn = swdmn + swdav/rslabs
    swumn = swumn + swuav/rslabs
    tllwtendmn = tllwtendmn + tllwtendav/rslabs
    tlswtendmn = tlswtendmn + tlswtendav/rslabs
    tlradlsmn  = tlradlsmn  + thlpcar

  end subroutine do_radstat

!> Write the statistics to file
  subroutine writeradstat
      use modmpi,    only : myid
      use modglobal, only : cexpnr,ifoutput,kmax,k1,zf,timee
      use modstat_nc, only: lnetcdf, writestat_nc
      use modgenstat, only: ncid_prof=>ncid,nrec_prof=>nrec

      implicit none
      real,dimension(k1,nvar) :: vars
      integer nsecs, nhrs, nminut,k


      nsecs   = nint(timee)
      nhrs    = int(nsecs/3600)
      nminut  = int(nsecs/60)-nhrs*60
      nsecs   = mod(nsecs,60)

      lwumn   = lwumn    /nsamples
      lwdmn   = lwdmn    /nsamples
      swdmn   = swdmn    /nsamples
      swumn   = swumn    /nsamples
      tllwtendmn = tllwtendmn /nsamples
      tlswtendmn = tlswtendmn /nsamples
      tlradlsmn  = tlradlsmn  /nsamples

  !     ----------------------
  !     2.0  write the fields
  !           ----------------

    if(myid==0)then
      open (ifoutput,file='radstat.'//cexpnr,position='append')
      write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
      '#--------------------------------------------------------'      &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '      &
      ,nhrs,':',nminut,':',nsecs      &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
      write (ifoutput,'(A/2A/2A)') &
          '#--------------------------------------------------------------------------' &
          ,'#LEV HGHT     LW_UP        LW_DN        SW_UP       SW_DN       ' &
          ,'TL_LW_TEND   TL_SW_TEND   TL_LS_TEND' &
          ,'#    (M)      (W/M^2)      (W/M^2)      (W/M^2)      (W/M^2)      ' &
          ,'(K/H)         (K/H)        (K/H)'
      do k=1,kmax
        write(ifoutput,'(I3,F8.2,7E13.4)') &
            k,zf(k),&
            lwumn(k),&
            lwdmn(k),&
            swdmn(k),&
            swumn(k),&
            tllwtendmn(k)*3600,&
            tlswtendmn(k)*3600,&
            tlradlsmn(k) *3600
      end do
      close (ifoutput)
      if (lnetcdf) then
        vars(:, 1) = tllwtendmn
        vars(:, 2) = tlswtendmn
        vars(:, 3) = tlradlsmn
        vars(:, 4) = lwumn
        vars(:, 5) = lwdmn
        vars(:, 6) = swumn
        vars(:, 7) = swdmn
        call writestat_nc(ncid_prof,nvar,ncname,vars(1:kmax,:),nrec_prof+1,kmax)
      end if
    end if ! end if(myid==0)

    lwumn = 0.0
    lwdmn = 0.0
    swdmn = 0.0
    swumn = 0.0
    tllwtendmn = 0.0
    tlswtendmn = 0.0
    tlradlsmn  = 0.0


  end subroutine writeradstat

!> Cleans up after the run
  subroutine exitradstat
    implicit none

    !deallocate variables that are needed in modradiation
    deallocate(lwdavl,lwuavl,swdavl,swuavl,tllwtendavl,tlswtendavl)

    if(.not.(lstat)) return

    deallocate(lwuav,lwdav,swdav,swuav)
    deallocate(tllwtendav,tlswtendav)
    deallocate(lwumn,lwdmn,swdmn,swumn)
    deallocate(tllwtendmn,tlswtendmn,tlradlsmn)



  end subroutine exitradstat


end module modradstat
