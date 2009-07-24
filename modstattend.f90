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
module modstattend

    !-----------------------------------------------------------------|
    !                                                                 |
    !*** *stattend*  calculates tendencies of prognostic variables    |
    !                                                                 |
    !      Thijs Heus      TU Delft     15/11/2007                    |
    !                                                                 |
    !     purpose.                                                    |
    !     --------                                                    |
    !                                                                 |
    !                                                                 |
    !    NOTE: A STATEMENT call stattend NEEDS TO BE SET BEFORE       |
    !    ADVECTION AND AFTER ADVECTION, SUBGRID, THE ADDONS, DAMPING  |
    !    AND PRESSURE SECTIONS IN program.f90                          |
    !                                                                 |
    !                                                                 |
    !_________________________ON OUTPUT_______________________________|
    !                                                                 |
    !     Various slab averages                                       |
    !                                                                 |
    !____________________SETTINGS_AND_SWITCHES________________________|
    !                     IN &NAMSTATTEND                             |
    !                                                                 |
    !    dtav           SAMPLING INTERVAL                             |
    !    timeav         WRITING INTERVAL                              |
    !    lstattend      SWITCH TO ENABLE TENDENCY STATISTICS          |
    !-----------------------------------------------------------------|
  implicit none
!   private
!   public :: initstattend, stattend, exitstattend
  save

  real    :: dtav, timeav,tnext,tnextwrite
  integer,parameter :: tend_tot=1,tend_start=1,tend_adv=2,tend_subg=3,tend_force=4,&
                       tend_rad=5,tend_ls=6,tend_micro=7, tend_topbound=8,tend_pois=9,tend_addon=10, tend_coriolis=11
  integer,parameter :: nrfields = 11
  integer :: nsamples
  logical :: ltend = .false.

  real, allocatable :: upmn(:,:),vpmn(:,:),wpmn(:,:),thlpmn(:,:),qtpmn(:,:)
  real, allocatable :: upav(:,:),vpav(:,:),wpav(:,:),thlpav(:,:),qtpav(:,:)
contains
  subroutine initstattend
    use modmpi,   only : mpierr,my_real,mpi_logical,comm3d,myid
    use modglobal,only : dtmax,imax,jmax,ifnamopt,fname_options,k1,dtav_glob,timeav_glob,ladaptive, dt_lim,btime
    implicit none
    integer :: ierr


    namelist/NAMSTATTEND/ &
    timeav,dtav,ltend

    dtav=dtav_glob;timeav=timeav_glob

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMSTATTEND,iostat=ierr)
      write(6 ,NAMSTATTEND)
      close(ifnamopt)
    end if

    call MPI_BCAST(dtav       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(timeav     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ltend      ,1,MPI_LOGICAL,0,comm3d,mpierr)

    tnext   = dtav - 1e-3 + btime
    tnextwrite = timeav - 1e-3 + btime
    nsamples = nint(timeav/dtav)

    if(.not.(ltend)) return
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if
    if (abs(timeav/dtav-nsamples)>1e-4) then
      stop 'timeav should be a integer multiple of dtav'
    end if

    allocate (upmn(k1,nrfields),vpmn(k1,nrfields),wpmn(k1,nrfields),thlpmn(k1,nrfields),qtpmn(k1,nrfields))
    allocate (upav(k1,nrfields),vpav(k1,nrfields),wpav(k1,nrfields),thlpav(k1,nrfields),qtpav(k1,nrfields))
    upmn = 0
    vpmn = 0
    wpmn = 0
    thlpmn = 0
    qtpmn = 0
    upav = 0
    vpav = 0
    wpav = 0
    thlpav = 0
    qtpav = 0

  end subroutine initstattend

  subroutine stattend(tendterm,lastterm)
    use modmpi,    only : myid,slabsum
    use modglobal, only : ih,jh,i1,j1,k1,rk3step,timee,dt_lim,rslabs
    use modfields, only : up,vp,wp,thlp,qtp
    implicit none
    integer, intent(in)           :: tendterm
    logical, intent(in), optional :: lastterm
    real, dimension(:),allocatable :: avfield

    if (.not.(ltend)) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if
    allocate(avfield(k1))

    avfield = 0.0
    call slabsum(avfield  ,1,k1,up  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    upav(:,tendterm) = avfield  - upav(:,tend_tot)
    upav(:,tend_tot)  = avfield

    avfield = 0.0
    call slabsum(avfield  ,1,k1,vp  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    vpav(:,tendterm) = avfield  - vpav(:,tend_tot)
    vpav(:,tend_tot)  = avfield

    avfield = 0.0
    call slabsum(avfield  ,1,k1,wp  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    wpav(:,tendterm) = avfield  - wpav(:,tend_tot)
    wpav(:,tend_tot)  = avfield

    avfield = 0.0
    call slabsum(avfield  ,1,k1,thlp  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    thlpav(:,tendterm) = avfield  - thlpav(:,tend_tot)
    thlpav(:,tend_tot)  = avfield

    avfield = 0.0
    call slabsum(avfield  ,1,k1,qtp  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    qtpav(:,tendterm) = avfield  - qtpav(:,tend_tot)
    qtpav(:,tend_tot)  = avfield

    if (present(lastterm)) then
    if (lastterm) then
      tnext = tnext+dtav
      upmn  = upmn  + upav /nsamples/rslabs
      vpmn  = vpmn  + vpav /nsamples/rslabs
      wpmn  = wpmn  + wpav /nsamples/rslabs
      qtpmn  = qtpmn  + qtpav /nsamples/rslabs
      thlpmn  = thlpmn  + thlpav /nsamples/rslabs
      upav  = 0.0
      vpav  = 0.0
      wpav  = 0.0
      qtpav  = 0.0
      thlpav  = 0.0
      if (timee>=tnextwrite) then
        tnextwrite = tnextwrite+timeav
        call writestattend
      end if
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
    end if
    end if

    deallocate(avfield)
  end subroutine stattend


  subroutine writestattend
    use modglobal, only : timee,ifoutput,kmax, zf, cexpnr
    use modfields, only : presf
    implicit none

    integer nsecs, nhrs, nminut,k

    nsecs   = nint(timee)
    nhrs    = int(nsecs/3600)
    nminut  = int(nsecs/60)-nhrs*60
    nsecs   = mod(nsecs,60)

    open(ifoutput,file='utend.'//cexpnr,position='append')
     write(ifoutput,'(//A,/A,I4,A,I2,A,I2,A)') &
         '#--------------------------------------------------------'      &
         ,'#--- AVERAGING TIMESTEP --- '      &
         ,nhrs,':',nminut,':',nsecs      &
         ,'   HRS:MIN:SEC AFTER INITIALIZATION '

      write (ifoutput,'(2A/A/A/2A/2A/2A/2A)') &
           '#-----------------------------------------------------------' &
           , '------------' &
           ,'#               --------   T E N D E N C I E S  U --------    ' &
           ,'#                                                           ' &
           ,'#                  |                               VELOCITY ' &
           ,' TENDENCIES                                    |           ' &
           ,'# LEV HEIGHT   PRES  |  ADVECT       DIFF         FORCES       CORIOLIS     LARGE SCALE    DAMP' &
           ,'         PRESSURE         ADD ON        TOTAL' &
           ,'#      (M)   (MB)  |  ---------   (M/S^2)     ----------' &
           ,'      ' &
           ,'#----------------------------------------------------------' &
           ,'----------------------------------------------------------'
      write(ifoutput,'(I4,F8.2,F7.1,9E13.5)') &
           (k, &
            zf       (k), &
            presf    (k)/100., &
            upmn   (k,tend_adv), &
            upmn   (k,tend_subg), &
            upmn   (k,tend_force), &
            upmn   (k,tend_coriolis), &
            upmn   (k,tend_ls), &
            upmn   (k,tend_topbound), &
            upmn   (k,tend_pois), &
            upmn   (k,tend_addon), &
            upmn   (k,tend_tot), &
                         k=1,kmax)
      close(ifoutput)

    open(ifoutput,file='vtend.'//cexpnr,position='append')
     write(ifoutput,'(//A,/A,I4,A,I2,A,I2,A)') &
         '#--------------------------------------------------------'      &
         ,'#--- AVERAGING TIMESTEP --- '      &
         ,nhrs,':',nminut,':',nsecs      &
         ,'   HRS:MIN:SEC AFTER INITIALIZATION '

      write (ifoutput,'(2A/A/A/2A/2A/2A/2A)') &
           '#-----------------------------------------------------------' &
           , '------------' &
           ,'#               --------   T E N D E N C I E S  V --------    ' &
           ,'#                                                           ' &
           ,'#                  |                               VELOCITY ' &
           ,' TENDENCIES                                    |           ' &
           ,'# LEV HEIGHT   PRES  |  ADVECT       DIFF         FORCES       CORIOLIS       LARGE SCALE    DAMP' &
           ,'         PRESSURE         ADD ON      TOTAL' &
           ,'#      (M)   (MB)  |  ---------   (M/S^2)     ----------' &
           ,'      ' &
           ,'#----------------------------------------------------------' &
           ,'----------------------------------------------------------'
      write(ifoutput,'(I4,F8.2,F7.1,9E13.5)') &
           (k, &
            zf       (k), &
            presf    (k)/100., &
            vpmn   (k,tend_adv), &
            vpmn   (k,tend_subg), &
            vpmn   (k,tend_force), &
            vpmn   (k,tend_coriolis), &
            vpmn   (k,tend_ls), &
            vpmn   (k,tend_topbound), &
            vpmn   (k,tend_pois), &
            vpmn   (k,tend_addon), &
            vpmn   (k,tend_tot), &
                         k=1,kmax)
      close(ifoutput)

    open(ifoutput,file='wtend.'//cexpnr,position='append')
     write(ifoutput,'(//A,/A,I4,A,I2,A,I2,A)') &
         '#--------------------------------------------------------'      &
         ,'#--- AVERAGING TIMESTEP --- '      &
         ,nhrs,':',nminut,':',nsecs      &
         ,'   HRS:MIN:SEC AFTER INITIALIZATION '

      write (ifoutput,'(2A/A/A/2A/2A/2A/2A)') &
           '#-----------------------------------------------------------' &
           , '------------' &
           ,'#               --------   T E N D E N C I E S  W --------    ' &
           ,'#                                                           ' &
           ,'#                  |                               VELOCITY ' &
           ,' TENDENCIES                                    |           ' &
           ,'# LEV HEIGHT   PRES  |  ADVECT       DIFF         FORCES       CORIOLIS       LARGE SCALE    DAMP' &
           ,'         PRESSURE         ADD ON      TOTAL' &
           ,'#      (M)   (MB)  |  ---------   (M/S^2)     ----------' &
           ,'      ' &
           ,'#----------------------------------------------------------' &
           ,'----------------------------------------------------------'
      write(ifoutput,'(I4,F8.2,F7.1,9E13.5)') &
           (k, &
            zf       (k), &
            presf    (k)/100., &
            wpmn   (k,tend_adv), &
            wpmn   (k,tend_subg), &
            wpmn   (k,tend_force), &
            wpmn   (k,tend_coriolis), &
            wpmn   (k,tend_ls), &
            wpmn   (k,tend_topbound), &
            wpmn   (k,tend_pois), &
            wpmn   (k,tend_addon), &
            wpmn   (k,tend_tot), &
                         k=1,kmax)
      close(ifoutput)

    open(ifoutput,file='tltend.'//cexpnr,position='append')
     write(ifoutput,'(//A,/A,I4,A,I2,A,I2,A)') &
         '#--------------------------------------------------------'      &
         ,'#--- AVERAGING TIMESTEP --- '      &
         ,nhrs,':',nminut,':',nsecs      &
         ,'   HRS:MIN:SEC AFTER INITIALIZATION '

      write (ifoutput,'(2A/A/A/2A/2A/2A/2A)') &
           '#-----------------------------------------------------------' &
           , '------------' &
           ,'#               --------   T E N D E N C I E S  TL --------    ' &
           ,'#                                                           ' &
           ,'#                  |                               THETA_L ' &
           ,' TENDENCIES                                    |           ' &
           ,'# LEV HEIGHT   PRES  |  ADVECT       DIFF         RADIATION     MICRO       BOUND' &
           ,'      LARGE SCALE         ADD ON      TOTAL' &
           ,'#      (M)   (MB)  |  ---------   (K/S)        ----------' &
           ,'      ' &
           ,'#----------------------------------------------------------' &
           ,'----------------------------------------------------------'
      write(ifoutput,'(I4,F8.2,F7.1,8E13.5)') &
           (k, &
            zf       (k), &
            presf    (k)/100., &
            thlpmn   (k,tend_adv), &
            thlpmn   (k,tend_subg), &
            thlpmn   (k,tend_rad), &
            thlpmn   (k,tend_micro), &
            thlpmn   (k,tend_topbound), &
            thlpmn   (k,tend_ls), &
            thlpmn   (k,tend_addon), &
            thlpmn   (k,tend_tot), &
                          k=1,kmax)
     close(ifoutput)

    open(ifoutput,file='qttend.'//cexpnr,position='append')
     write(ifoutput,'(//A,/A,I4,A,I2,A,I2,A)') &
         '#--------------------------------------------------------'      &
         ,'#--- AVERAGING TIMESTEP --- '      &
         ,nhrs,':',nminut,':',nsecs      &
         ,'   HRS:MIN:SEC AFTER INITIALIZATION '

      write (ifoutput,'(2A/A/A/2A/2A/2A/2A)') &
           '#-----------------------------------------------------------' &
           , '------------' &
           ,'#               --------   T E N D E N C I E S  QT --------    ' &
           ,'#                                                           ' &
           ,'#                  |                               VELOCITY ' &
           ,' TENDENCIES                                    |           ' &
           ,'# LEV HEIGHT   PRES  |  ADVECT       DIFF     MICRO       BOUND' &
           ,'      LARGE SCALE         ADD ON      TOTAL' &
           ,'#      (M)   (MB)  |  ---------   (KG/KG/S)     ----------' &
           ,'      ' &
           ,'#----------------------------------------------------------' &
           ,'----------------------------------------------------------'
      write(ifoutput,'(I4,F8.2,F7.1,7E13.5)') &
           (k, &
            zf       (k), &
            presf    (k)/100., &
            qtpmn   (k,tend_adv), &
            qtpmn   (k,tend_subg), &
            qtpmn   (k,tend_micro), &
            qtpmn   (k,tend_topbound), &
            qtpmn   (k,tend_ls), &
            qtpmn   (k,tend_addon), &
            qtpmn   (k,tend_tot), &
                          k=1,kmax)
     close(ifoutput)

  end subroutine writestattend

  subroutine exitstattend
  implicit none
   if(.not.(ltend)) return

    deallocate (upmn,vpmn,wpmn,thlpmn,qtpmn)
    deallocate (upav,vpav,wpav,thlpav,qtpav)

  end subroutine exitstattend
end module