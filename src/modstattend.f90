!> \file modstattend.f90
!!  Calculates the tendencies of the main fields


!>
!!  Calculates the tendencies of the main fields
!>
!! Profiles of the individual terms of the prognostic equations.  Written to *tend.expnr
!! If netcdf is true, this module also writes in the profiles.expnr.nc output
!!  \author Thijs Heus, MPI
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
module modstattend
  use modglobal, only : longint

  implicit none
!   private
!   public :: initstattend, stattend, exitstattend
  save
!NetCDF variables
  integer,parameter :: nvar = 43
  integer :: ncid,nrec = 0
  character(80) :: fname = 'stattend.xxx.nc'
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav, timeav
  integer(kind=longint) :: idtav,itimeav,tnext,tnextwrite
  integer,parameter :: tend_tot=1,tend_start=1,tend_adv=2,tend_subg=3,tend_force=4,&
                       tend_rad=5,tend_ls=6,tend_micro=7, tend_topbound=8,tend_pois=9,tend_addon=10, tend_coriolis=11
  integer,parameter :: nrfields = 11
  integer :: nsamples
  logical :: ltend = .false.

  real, allocatable :: upmn(:,:),vpmn(:,:),wpmn(:,:),thlpmn(:,:),qtpmn(:,:)
  real, allocatable :: upav(:,:),vpav(:,:),wpav(:,:),thlpav(:,:),qtpav(:,:)
contains
!> Initialization routine, reads namelists and inits variables
subroutine initstattend
    use modmpi,   only : mpierr,my_real,mpi_logical,comm3d,myid
    use modglobal,only : cexpnr,dtmax,imax,jmax,kmax,ifnamopt,fname_options,k1,dtav_glob,timeav_glob,ladaptive, dt_lim,btime,kmax,tres,ifoutput
    use modstat_nc, only : lnetcdf, open_nc,define_nc,ncinfo,writestat_dims_nc
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav,ncid_prof=>ncid

    implicit none
    integer :: ierr


    namelist/NAMSTATTEND/ &
    timeav,dtav,ltend

    dtav=dtav_glob;timeav=timeav_glob

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMSTATTEND,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMSTATTEND'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMSTATTEND'
      endif
      write(6 ,NAMSTATTEND)
      close(ifnamopt)
    end if

    call MPI_BCAST(dtav       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(timeav     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ltend      ,1,MPI_LOGICAL,0,comm3d,mpierr)

    idtav = dtav/tres
    itimeav = timeav/tres

    tnext      = idtav   +btime
    tnextwrite = itimeav +btime
    nsamples = itimeav/idtav
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

    if(myid==0)then
      open (ifoutput,file='utend.'//cexpnr,status='replace')  
      close (ifoutput)
      open (ifoutput,file='vtend.'//cexpnr,status='replace')  
      close (ifoutput)
      open (ifoutput,file='wtend.'//cexpnr,status='replace')  
      close (ifoutput)
      open (ifoutput,file='tltend.'//cexpnr,status='replace')  
      close (ifoutput)
      open (ifoutput,file='qttend.'//cexpnr,status='replace')  
      close (ifoutput)
    endif

    if (lnetcdf) then
    idtav = dtav/tres
    itimeav = timeav/tres

    tnext      = idtav   +btime
    tnextwrite = itimeav +btime
    nsamples = itimeav/idtav
     if (myid==0) then
        fname(10:12) = cexpnr
        call ncinfo(tncname(1,:),'time','Time','s','time')
        call ncinfo(ncname( 1,:),'utendadv','U advective tendency','m/s^2','tt')
        call ncinfo(ncname( 2,:),'utenddif','U diffusive tendency','m/s^2','tt')
        call ncinfo(ncname( 3,:),'utendfor','U tendency due to other forces','m/s^2','tt')
        call ncinfo(ncname( 4,:),'utendcor','U coriolis tendency','m/s^2','tt')
        call ncinfo(ncname( 5,:),'utendls','U large scale tendency','m/s^2','tt')
        call ncinfo(ncname( 6,:),'utendtop','U top boundary tendency','m/s^2','tt')
        call ncinfo(ncname( 7,:),'utendpois','U pressure gradient tendency','m/s^2','tt')
        call ncinfo(ncname( 8,:),'utendaddon','U in addons tendency','m/s^2','tt')
        call ncinfo(ncname( 9,:),'utendtot','U total tendency','m/s^2','tt')
        call ncinfo(ncname(10,:),'vtendadv','V advective tendency','m/s^2','tt')
        call ncinfo(ncname(11,:),'vtenddif','V diffusive tendency','m/s^2','tt')
        call ncinfo(ncname(12,:),'vtendfor','V tendency due to other forces','m/s^2','tt')
        call ncinfo(ncname(13,:),'vtendcor','V coriolis tendency','m/s^2','tt')
        call ncinfo(ncname(14,:),'vtendls','V large scale tendency','m/s^2','tt')
        call ncinfo(ncname(15,:),'vtendtop','V top boundary tendency','m/s^2','tt')
        call ncinfo(ncname(16,:),'vtendpois','V pressure gradient tendency','m/s^2','tt')
        call ncinfo(ncname(17,:),'vtendaddon','V in addons tendency','m/s^2','tt')
        call ncinfo(ncname(18,:),'vtendtot','V total tendency','m/s^2','tt')
        call ncinfo(ncname(19,:),'wtendadv','W advective tendency','m/s^2','mt')
        call ncinfo(ncname(20,:),'wtenddif','W diffusive tendency','m/s^2','mt')
        call ncinfo(ncname(21,:),'wtendfor','W tendency due to other forces','m/s^2','mt')
        call ncinfo(ncname(22,:),'wtendcor','W coriolis tendency','m/s^2','mt')
        call ncinfo(ncname(23,:),'wtendls','W large scale tendency','m/s^2','mt')
        call ncinfo(ncname(24,:),'wtendtop','W top boundary tendency','m/s^2','mt')
        call ncinfo(ncname(25,:),'wtendpois','W pressure gradient tendency','m/s^2','mt')
        call ncinfo(ncname(26,:),'wtendaddon','W in addons tendency','m/s^2','mt')
        call ncinfo(ncname(27,:),'wtendtot','W total tendency','m/s^2','mt')
        call ncinfo(ncname(28,:),'tltendadv','theta_l advective tendency','K/s','tt')
        call ncinfo(ncname(29,:),'tltenddif','theta_l diffusive tendency','K/s','tt')
        call ncinfo(ncname(30,:),'tltendrad','theta_l radiative tendency','K/s','tt')
        call ncinfo(ncname(31,:),'tltendmicro','theta_l microphysical tendency','K/s','tt')
        call ncinfo(ncname(32,:),'tltendls','theta_l large scale tendency','K/s','tt')
        call ncinfo(ncname(33,:),'tltendtop','theta_l  top boundary tendency','K/s','tt')
        call ncinfo(ncname(34,:),'tltendaddon','theta_l in addons tendency','K/s','tt')
        call ncinfo(ncname(35,:),'tltendtot','theta_l total tendency','K/s','tt')
        call ncinfo(ncname(36,:),'qttendadv','total water content advective tendency','kg/kg/s','tt')
        call ncinfo(ncname(37,:),'qttenddif','total water content diffusive tendency','kg/kg/s','tt')
        call ncinfo(ncname(38,:),'qttendrad','total water content radiative tendency','kg/kg/s','tt')
        call ncinfo(ncname(39,:),'qttendmicro','total water content microphysical tendency','kg/kg/s','tt')
        call ncinfo(ncname(40,:),'qttendls','total water content large scale tendency','kg/kg/s','tt')
        call ncinfo(ncname(41,:),'qttendtop','total water content  top boundary tendency','kg/kg/s','tt')
        call ncinfo(ncname(42,:),'qttendaddon','total water content in addons tendency','kg/kg/s','tt')
        call ncinfo(ncname(43,:),'qttendtot','total water content total tendency','kg/kg/s','tt')


        call define_nc( ncid_prof, NVar, ncname)
     end if

   end if

  end subroutine initstattend

!> Performs the statistics, keeps track of what the tendencies were last time, and what they are this time.
  subroutine stattend(tendterm,lastterm)
    use modmpi,    only : myid,slabsum
    use modglobal, only : ih,jh,i1,j1,kmax,k1,rk3step,timee,dt_lim,rslabs,btime
    use modfields, only : up,vp,wp,thlp,qtp
    implicit none
    integer, intent(in)           :: tendterm !< name of the term to write down
    logical, intent(in), optional :: lastterm !< true if this is the last term of the equations; the write routine is entered.
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
      tnext = tnext+idtav
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
        tnextwrite = tnextwrite+itimeav
        call writestattend
      end if
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
    end if
    end if

    deallocate(avfield)
  end subroutine stattend

!> Write the statistics to file
  subroutine writestattend
    use modglobal, only : rtimee,ifoutput,kmax,k1, zf, cexpnr
    use modmpi,    only : myid
    use modfields, only : presf
      use modstat_nc, only: lnetcdf, writestat_nc
      use modgenstat, only: ncid_prof=>ncid,nrec_prof=>nrec
    implicit none
    real,dimension(k1,nvar) :: vars
    integer nsecs, nhrs, nminut,k

    nsecs   = nint(rtimee)
    nhrs    = int(nsecs/3600)
    nminut  = int(nsecs/60)-nhrs*60
    nsecs   = mod(nsecs,60)

    if(myid == 0) then
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
      if (lnetcdf) then
        vars(:, 1) =upmn(:,tend_adv)
        vars(:, 2) =upmn(:,tend_subg)
        vars(:, 3) =upmn(:,tend_force)
        vars(:, 4) =upmn(:,tend_coriolis)
        vars(:, 5) =upmn(:,tend_ls)
        vars(:, 6) =upmn(:,tend_topbound)
        vars(:, 7) =upmn(:,tend_pois)
        vars(:, 8) =upmn(:,tend_addon)
        vars(:, 9) =upmn(:,tend_tot)
        vars(:,10) =vpmn(:,tend_adv)
        vars(:,11) =vpmn(:,tend_subg)
        vars(:,12) =vpmn(:,tend_force)
        vars(:,13) =vpmn(:,tend_coriolis)
        vars(:,14) =vpmn(:,tend_ls)
        vars(:,15) =vpmn(:,tend_topbound)
        vars(:,16) =vpmn(:,tend_pois)
        vars(:,17) =vpmn(:,tend_addon)
        vars(:,18) =vpmn(:,tend_tot)
        vars(:,19) =wpmn(:,tend_adv)
        vars(:,20) =wpmn(:,tend_subg)
        vars(:,21) =wpmn(:,tend_force)
        vars(:,22) =wpmn(:,tend_coriolis)
        vars(:,23) =wpmn(:,tend_ls)
        vars(:,24) =wpmn(:,tend_topbound)
        vars(:,25) =wpmn(:,tend_pois)
        vars(:,26) =wpmn(:,tend_addon)
        vars(:,27) =wpmn(:,tend_tot)
        vars(:,28) =thlpmn(:,tend_adv)
        vars(:,29) =thlpmn(:,tend_subg)
        vars(:,30) =thlpmn(:,tend_rad)
        vars(:,31) =thlpmn(:,tend_micro)
        vars(:,32) =thlpmn(:,tend_ls)
        vars(:,33) =thlpmn(:,tend_topbound)
        vars(:,34) =thlpmn(:,tend_addon)
        vars(:,35) =thlpmn(:,tend_tot)
        vars(:,36) =qtpmn(:,tend_adv)
        vars(:,37) =qtpmn(:,tend_subg)
        vars(:,38) =qtpmn(:,tend_rad)
        vars(:,39) =qtpmn(:,tend_micro)
        vars(:,40) =qtpmn(:,tend_ls)
        vars(:,41) =qtpmn(:,tend_topbound)
        vars(:,42) =qtpmn(:,tend_addon)
        vars(:,43) =qtpmn(:,tend_tot)
        call writestat_nc(ncid_prof,nvar,ncname,vars(1:kmax,:),nrec_prof+1,kmax)
      end if
    end if

  end subroutine writestattend
!> Cleans up after the run
  subroutine exitstattend
  implicit none
   if(.not.(ltend)) return

    deallocate (upmn,vpmn,wpmn,thlpmn,qtpmn)
    deallocate (upav,vpav,wpav,thlpav,qtpav)

  end subroutine exitstattend
end module
