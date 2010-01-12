!> \file modbulkmicrostat.f90
!!  Calculates profiles coming from the bulkmicrophysics


!>
!!  Calculates profiles coming from the bulkmicrophysics
!>
!! Profiles coming from the bulkmicrophysics. Written to precep.expnr for the
!! rain rates etc., and to qlptend.expnr, nptend.expnr and qtptend.expnr for the
!! tendencies is rain water content, droplet number, and total water content,
!! respectively.
!! If netcdf is true, this module also writes in the profiles.expnr.nc output
!!  \author Olivier Geoffroy, KNMI
!!  \author Johan van de Dussen, TU Delft
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
module modbulkmicrostat

implicit none
private
PUBLIC  :: initbulkmicrostat, bulkmicrostat, exitbulkmicrostat, bulkmicrotend
save
!NetCDF variables
  integer,parameter :: nvar = 23
  integer :: ncid,nrec = 0
  character(80) :: fname = 'bulkmicrostat.xxx.nc'
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname

  real          :: dtav, timeav, tnext, tnextwrite
  integer          :: nsamples
  logical          :: lmicrostat = .false.
  integer, parameter      :: nrfields = 5  , &
                 iauto    = 2 , &
                  iaccr    = 3 , &
               ievap    = 4 , &
               ised      = 5
  real, allocatable, dimension(:,:)  :: Npav    , &
               Npmn    , &
               qlpav  , &
               qlpmn  , &
               qtpav  , &
               qtpmn
  real, allocatable, dimension(:)    :: precavl  , &
               precav  , &
               precmn  , &
               preccountavl  , &
               preccountav  , &
               preccountmn  , &
               prec_prcavl  , &
               prec_prcav  , &
               prec_prcmn  , &
               cloudcountavl, &
               cloudcountav  , &
               cloudcountmn  , &
               raincountavl  , &
               raincountav  , &
               raincountmn  , &
               Nrrainavl  , &
               Nrrainav  , &
               Nrrainmn  , &
               qravl  , &
               qrav    , &
               qrmn    , &
               Dvravl  , &
               Dvrav  , &
               Dvrmn

contains
!> Initialization routine, reads namelists and inits variables
subroutine initbulkmicrostat
    use modmpi,    only  : myid, mpi_logical, my_real, comm3d, mpierr
    use modglobal, only  : ifnamopt, fname_options, cexpnr, ifoutput, &
              dtav_glob, timeav_glob, ladaptive, k1,kmax, dtmax,btime
    use modstat_nc, only : lnetcdf, open_nc,define_nc,redefine_nc,ncinfo,writestat_dims_nc
    use modgenstat, only : dtav_prof=>dtav, timeav_prof=>timeav,ncid_prof=>ncid
    use modmicrodata,only: imicro, imicro_bulk
    implicit none
    integer      :: ierr

    namelist/NAMBULKMICROSTAT/ &
    lmicrostat, dtav, timeav

    if (imicro /=imicro_bulk) return

    dtav  = dtav_glob
    timeav  = timeav_glob
    if(myid==0)then
      open (ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMBULKMICROSTAT,iostat=ierr)
      write(6,NAMBULKMICROSTAT)
      close(ifnamopt)
    end if

    call MPI_BCAST(lmicrostat  ,1,MPI_LOGICAL  ,0,comm3d,mpierr)
    call MPI_BCAST(dtav    ,1,MY_REAL  ,0,comm3d,mpierr)
    call MPI_BCAST(timeav    ,1,MY_REAL  ,0,comm3d,mpierr)

    tnext    = dtav - 1e-3+btime
    tnextwrite  = timeav - 1e-3+btime
    nsamples  = nint(timeav/dtav)

    if (.not. lmicrostat) return
    if (abs(timeav/dtav - nsamples) > 1e-4) then
      stop 'timeav must be an integer multiple of dtav (NAMBULKMICROSTAT)'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax - nint(dtav/dtmax)) > 1e-4) then
      stop 'dtav must be an integer multiple of dtmax (NAMBULKMICROSTAT)'
    end if

    allocate(Npav    (k1, nrfields)  , &
       Npmn    (k1, nrfields)  , &
       qlpav    (k1, nrfields)  , &
       qlpmn    (k1, nrfields)  , &
       qtpav    (k1, nrfields)  , &
       qtpmn    (k1, nrfields)  )
    allocate(precavl  (k1)    , &
       precav    (k1)    , &
       precmn    (k1)    , &
       preccountavl  (k1)    , &
       preccountav  (k1)    , &
       preccountmn  (k1)    , &
       prec_prcavl  (k1)    , &
       prec_prcav  (k1)    , &
       prec_prcmn  (k1)    , &
       cloudcountavl  (k1)    , &
       cloudcountav  (k1)    , &
       cloudcountmn  (k1)    , &
       raincountavl  (k1)    , &
       raincountav  (k1)    , &
       raincountmn  (k1)    , &
       Nrrainavl  (k1)    , &
       Nrrainav  (k1)    , &
       Nrrainmn  (k1)    , &
       qravl    (k1)    , &
       qrav    (k1)    , &
       qrmn    (k1)    , &
       Dvravl    (k1)    , &
       Dvrav    (k1)    , &
       Dvrmn    (k1)    )
    Npmn    = 0.0
    qlpmn    = 0.0
    qtpmn    = 0.0
    precmn    = 0.0
    preccountmn  = 0.0
    prec_prcmn  = 0.0
    cloudcountmn  = 0.0
    raincountmn  = 0.0
    Nrrainmn  = 0.0
    qrmn    = 0.0
    Dvrmn    = 0.0

    if (myid == 0) then
      open (ifoutput,file = 'precep.'//cexpnr ,status = 'replace')
      close(ifoutput)
      open (ifoutput,file = 'nptend.'//cexpnr ,status = 'replace')
      close(ifoutput)
      open (ifoutput,file = 'qlptend.'//cexpnr,status = 'replace')
      close(ifoutput)
      open (ifoutput,file = 'qtptend.'//cexpnr,status = 'replace')
      close(ifoutput)
    end if
    if (lnetcdf) then
      dtav = dtav_prof
      timeav = timeav_prof
      if (myid==0) then
        fname(15:17) = cexpnr
        call ncinfo(tncname(1,:),'time','Time','s','time')
        call ncinfo(ncname( 1,:),'cfrac','Cloud fraction','-','tt')
        call ncinfo(ncname( 2,:),'rainrate','Echo Rain Rate','W/m^2','tt')
        call ncinfo(ncname( 3,:),'preccount','Preccount','W/m^2','mt')
        call ncinfo(ncname( 4,:),'nrrain','nrrain','W/m^2','mt')
        call ncinfo(ncname( 5,:),'raincount','raincount','W/m^2','tt')
        call ncinfo(ncname( 6,:),'precmn','precmn','W/m^2','mt')
        call ncinfo(ncname( 7,:),'dvrmn','dvrmn','W/m^2','tt')
        call ncinfo(ncname( 8,:),'qrmn','qrmn','W/m^2','tt')
        call ncinfo(ncname( 9,:),'npauto','Autoconversion rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(10,:),'npaccr','Accretion rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(11,:),'npsed','Sedimentation rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(12,:),'npevap','Evaporation rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(13,:),'nptot','Total rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(14,:),'qrpauto','Autoconversion rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(15,:),'qrpaccr','Accretion rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(16,:),'qrpsed','Sedimentation rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(17,:),'qrpevap','Evaporation rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(18,:),'qrptot','Total rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(19,:),'qtpauto','Autoconversion total water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(20,:),'qtpaccr','Accretion total water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(21,:),'qtpsed','Sedimentation total water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(22,:),'qtpevap','Evaporation total water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(23,:),'qtptot','Total total water content tendency','kg/kg/s','tt')
        call open_nc(fname,ncid,n3=kmax)
        call define_nc( ncid, 1, tncname)
        call writestat_dims_nc(ncid)
        call redefine_nc(ncid_prof)
        call define_nc( ncid_prof, NVar, ncname)
      end if

   end if

  end subroutine initbulkmicrostat

!------------------------------------------------------------------------------!
!> General routine, does the timekeeping
  subroutine bulkmicrostat
    use modmpi,    only  : myid
    use modglobal,    only  : rk3step, timee, dt_lim
    implicit none
    if (.not. lmicrostat)  return
    if (rk3step /= 3)  return
    if (timee == 0)    return
    if (timee < tnext .and. timee < tnextwrite) then
      dt_lim  = minval((/dt_lim, tnext - timee, tnextwrite - timee/))
      return
    end if
    if (timee >= tnext) then
      tnext = tnext + dtav
      call dobulkmicrostat
    end if
    if (timee >= tnextwrite) then
      tnextwrite = tnextwrite + timeav
      call writebulkmicrostat
    end if

  end subroutine bulkmicrostat

!------------------------------------------------------------------------------!
!> Performs the calculations for rainrate etc.
  subroutine dobulkmicrostat
    use modmpi,    only  : my_real, mpi_sum, comm3d, mpierr
    use modglobal,    only  : i1, j1, k1, rslabs
    use modmicrodata,  only  : qc, qr, precep, Dvr, Nr, epscloud, epsqr, epsprec
    implicit none

    integer      :: k

    precav      = 0.0
    preccountav    = 0.0
    prec_prcav    = 0.0
    cloudcountav    = 0.0
    raincountav    = 0.0
    Nrrainav    = 0.0
    qrav      = 0.0
    Dvrav      = 0.0

    do k = 1,k1
      cloudcountavl(k)  = count(qc      (2:i1,2:j1,k) > epscloud)
      raincountavl (k)  = count(qr      (2:i1,2:j1,k) > epsqr)
      preccountavl (k)  = count(precep  (2:i1,2:j1,k) > epsprec)
      prec_prcavl  (k)  = sum  (precep  (2:i1,2:j1,k)  , precep(2:i1,2:j1,k) > epsprec)
      Dvravl       (k)  = sum  (Dvr  (2:i1,2:j1,k)  , qr  (2:i1,2:j1,k) > epsqr)
      Nrrainavl    (k)  = sum  (Nr  (2:i1,2:j1,k))
      precavl      (k)  = sum  (precep  (2:i1,2:j1,k))
      qravl        (k)  = sum  (qr  (2:i1,2:j1,k))
    end do

    call MPI_ALLREDUCE(cloudcountavl,cloudcountav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(raincountavl ,raincountav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(preccountavl  ,preccountav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(prec_prcavl  ,prec_prcav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(Dvravl  ,Dvrav    ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(Nrrainavl  ,Nrrainav  ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(precavl  ,precav    ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(qravl  ,qrav    ,k1,MY_REAL,MPI_SUM,comm3d,mpierr)

    cloudcountmn  = cloudcountmn  +  cloudcountav  /rslabs
    raincountmn  = raincountmn  +  raincountav  /rslabs
    preccountmn  = preccountmn  +  preccountav  /rslabs
    prec_prcmn  = prec_prcmn  +  prec_prcav  /rslabs
    Dvrmn    = Dvrmn    +  Dvrav    /rslabs
    Nrrainmn  = Nrrainmn  +  Nrrainav  /rslabs
    precmn    = precmn  +  precav    /rslabs
    qrmn    = qrmn    +  qrav    /rslabs

  end subroutine dobulkmicrostat

!------------------------------------------------------------------------------!
!> Performs the calculations for the tendencies etc.
  subroutine bulkmicrotend
    use modmpi,    only  : slabsum
    use modglobal,    only  : rk3step, timee, dt_lim, k1, ih, i1, jh, j1, rslabs
    use modfields,    only  : qtp
    use modmicrodata,  only  : qrp, Nrp
    implicit none

    real, dimension(:), allocatable  :: avfield
    integer        :: ifield = 0

    if (.not. lmicrostat)  return
    if (rk3step /= 3)  return
    if (timee == 0)    return
    if (timee < tnext .and. timee < tnextwrite) then
      dt_lim  = minval((/dt_lim, tnext - timee, tnextwrite - timee/))
      return
    end if
!    tnext = tnext+dtav

    allocate(avfield(k1))

    ifield    = mod(ifield, nrfields) + 1

    avfield    = 0.0
    call slabsum(avfield  ,1,k1,Nrp  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    Npav(:,ifield)  = avfield - sum(Npav  (:,1:ifield-1),2)

    avfield    = 0.0
    call slabsum(avfield  ,1,k1,qrp  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    qlpav(:,ifield) = avfield - sum(qlpav  (:,1:ifield-1),2)

    avfield    = 0.0
    call slabsum(avfield  ,1,k1,qtp  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    qtpav(:,ifield) = avfield - sum(qtpav  (:,1:ifield-1),2)

    if (ifield == nrfields) then
      Npmn    = Npmn  + Npav  /nsamples/rslabs
      qlpmn    = qlpmn  + qlpav /nsamples/rslabs
      qtpmn    = qtpmn  + qtpav /nsamples/rslabs
      Npav    = 0.0
      qlpav    = 0.0
      qtpav    = 0.0
    end if

    deallocate(avfield)
  end subroutine bulkmicrotend

!------------------------------------------------------------------------------!
!> Write the stats to file
  subroutine writebulkmicrostat
    use modmpi,    only  : myid
    use modglobal,    only  : timee, ifoutput, cexpnr, k1,kmax, &
              rlv, zf
    use modfields,    only  : presf
    use modmicrodata,  only  : rhoz
      use modstat_nc, only: lnetcdf, writestat_nc
      use modgenstat, only: ncid_prof=>ncid,nrec_prof=>nrec

      implicit none
      real,dimension(k1,nvar) :: vars

    integer    :: nsecs, nhrs, nminut
    integer    :: k

    nsecs    = nint(timee)
    nhrs    = int (nsecs/3600)
    nminut    = int (nsecs/60)-nhrs*60
    nsecs    = mod (nsecs,60)

    cloudcountmn    = cloudcountmn  /nsamples
                raincountmn     = raincountmn   /nsamples
                preccountmn     = preccountmn   /nsamples
                prec_prcmn      = prec_prcmn    /nsamples
                Dvrmn           = Dvrmn         /nsamples
                Nrrainmn        = Nrrainmn      /nsamples
                precmn          = precmn        /nsamples
                qrmn            = qrmn          /nsamples

    where (raincountmn > 0.)
      Dvrmn        = Dvrmn / raincountmn
    elsewhere
      Dvrmn = 0.0
    end where
    where (preccountmn > 0.)
      prec_prcmn = prec_prcmn/preccountmn
    elsewhere
      prec_prcmn = 0.0
    end where

    if (myid == 0) then
    open (ifoutput,file='precep.'//cexpnr,position='append')
    write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
      '#-------------------------------------------------------------'   &
      ,'---------------------------------)'           &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
      ,nhrs,':',nminut,':',nsecs             &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
    write (ifoutput,'(2A/A/A/2A/2A/2A)')             &
      '#------------------------------------------------------------'     &
      ,'------------'                 &
      ,'#               --------   PRECIPITATION ------    '       &
      ,'#                                                           '     &
      ,'# LEV HEIGHT   RHO(k)  PRES  |CLOUDCOVER  ECHORAINRATE  PRECCOUNT '   &
      ,'    NRRAIN      RAINCOUNT     PREC(k)     <Dvr(k)>     <qr(k)>'   &
      ,'#      (M)             (MB)  |----------  ---W/M2----   --------- '   &
      ,'    ------      ---------     -------     --------    ---------'   &
      ,'#-----------------------------------------------------------------'   &
      ,'---------------------------------------------------------------'
    write(ifoutput,'(I4,F8.2,F8.3,F7.1,8E13.5)') &
      (k          , &
      zf    (k)      , &
      rhoz    (2,2,k)      , &
      presf    (k)/100.    , &
      cloudcountmn  (k)      , &
      prec_prcmn  (k)*rhoz(2,2,k)*rlv  , &
      preccountmn  (k)      , &
      Nrrainmn  (k)      , &
      raincountmn  (k)      , &
      precmn    (k)*rhoz(2,2,k)*rlv  , &
      Dvrmn    (k)      , &
      qrmn    (k)      , &
      k=1,kmax)
    close(ifoutput)

    open (ifoutput,file='nptend.'//cexpnr,position='append')
    write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
      '#-------------------------------------------------------------'   &
      ,'---------------------------------)'           &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
      ,nhrs,':',nminut,':',nsecs             &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
    write (ifoutput,'(2A/A/A/2A/A/A)')             &
      '#------------------------------------------------------------'     &
      , '------------'               &
      ,'#               --------   T E N D E N C I E S NRAIN ------    '     &
      ,'#                                                           '     &
      ,'# LEV HEIGHT   PRES  |  AUTO         ACCR          SEDIM    '     &
      ,'     EVAP         TOT '             &
      ,'#      (M)   (MB)  |  ---------   (#/M3/S)      ----------'     &
      ,'#-----------------------------------------------------------'
    write(ifoutput,'(I4,F8.2,F7.1,5E13.5)') &
      (k          , &
      zf    (k)      , &
      presf    (k)/100.    , &
      Npmn    (k,iauto)    , &
      Npmn    (k,iaccr)    , &
      Npmn    (k,ised)    , &
      Npmn    (k,ievap)    , &
      sum(Npmn  (k,2:nrfields))    , &
      k=1,kmax)
    close(ifoutput)

    open (ifoutput,file='qlptend.'//cexpnr,position='append')
    write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
      '#-------------------------------------------------------------'   &
      ,'---------------------------------)'           &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
      ,nhrs,':',nminut,':',nsecs             &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
    write (ifoutput,'(2A/A/A/2A/A/A)')             &
      '#------------------------------------------------------------'     &
      , '------------'               &
      ,'#               --------   T E N D E N C I E S QRAIN ------    '   &
      ,'#                                                           '     &
      ,'# LEV HEIGHT   PRES  |  AUTO         ACCR          SEDIM    '     &
      ,'     EVAP         TOT '             &
      ,'#      (M)   (MB)  |  ---------   (KG/KG/S)      ----------'     &
      ,'#-----------------------------------------------------------'
    write(ifoutput,'(I4,F8.2,F7.1,5E13.5)') &
      (k          , &
      zf    (k)      , &
      presf    (k)/100.    , &
      qlpmn    (k,iauto)    , &
      qlpmn    (k,iaccr)    , &
      qlpmn    (k,ised)    , &
      qlpmn    (k,ievap)    , &
      sum(qlpmn  (k,2:nrfields))    , &
                        k=1,kmax)
    close(ifoutput)

    open (ifoutput,file='qtptend.'//cexpnr,position='append')
    write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
      '#-------------------------------------------------------------'   &
      ,'---------------------------------)'           &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
      ,nhrs,':',nminut,':',nsecs             &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
    write (ifoutput,'(2A/A/A/2A/A/A)')             &
      '#------------------------------------------------------------'     &
      , '------------'               &
      ,'#               --------   T E N D E N C I E S QTP ------    '   &
      ,'#                                                           '     &
      ,'# LEV HEIGHT   PRES  |  AUTO         ACCR          SEDIM    '     &
      ,'     EVAP         TOT '             &
      ,'#      (M)   (MB)  |  ---------   (KG/KG/S)      ----------'     &
      ,'#-----------------------------------------------------------'
    write(ifoutput,'(I4,F8.2,F7.1,5E13.5)') &
      (k          , &
      zf    (k)      , &
      presf    (k)/100.    , &
      qtpmn    (k,iauto)    , &
      qtpmn    (k,iaccr)    , &
      qtpmn    (k,ised)    , &
      qtpmn    (k,ievap)    , &
      sum    (qtpmn(k,2:nrfields))  , &
      k=1,kmax)
      close(ifoutput)
      if (lnetcdf) then
        vars(:, 1) = cloudcountmn
        vars(:, 2) = prec_prcmn  (:)*rhoz(2,2,:)*rlv
        vars(:, 3) = preccountmn  (:)
        vars(:, 4) = Nrrainmn  (:)
        vars(:, 5) = raincountmn  (:)
        vars(:, 6) = precmn    (:)*rhoz(2,2,:)*rlv
        vars(:, 7) = Dvrmn    (:)
        vars(:, 8) = qrmn    (:)
        vars(:, 9) =Npmn    (k,iauto)
        vars(:,10) =Npmn    (k,iaccr)
        vars(:,11) =Npmn    (k,ised)
        vars(:,12) =Npmn    (k,ievap)
        vars(:,13) =sum(Npmn  (k,2:nrfields))
        vars(:,14) =qlpmn    (k,iauto)
        vars(:,15) =qlpmn    (k,iaccr)
        vars(:,16) =qlpmn    (k,ised)
        vars(:,17) =qlpmn    (k,ievap)
        vars(:,18) =sum(qlpmn  (k,2:nrfields))
        vars(:,19) =qtpmn    (k,iauto)
        vars(:,20) =qtpmn    (k,iaccr)
        vars(:,21) =qtpmn    (k,ised)
        vars(:,22) =qtpmn    (k,ievap)
        vars(:,23) = sum    (qtpmn(k,2:nrfields))
        call writestat_nc(ncid_prof,nvar,ncname,vars(1:kmax,:),nrec_prof,kmax)
      end if

    end if

    cloudcountmn    = 0.0
    raincountmn    = 0.0
    preccountmn    = 0.0
    prec_prcmn    = 0.0
    Dvrmn      = 0.0
    Nrrainmn    = 0.0
    precmn      = 0.0
    qrmn      = 0.0
    Npmn      = 0.0
    qlpmn      = 0.0
    qtpmn      = 0.0

  end subroutine writebulkmicrostat

!------------------------------------------------------------------------------!

  subroutine exitbulkmicrostat
    implicit none
    if (.not. lmicrostat)  return

    deallocate(Npav      , &
         Npmn      , &
         qlpav    , &
         qlpmn    , &
         qtpav    , &
         qtpmn    )
    deallocate(precavl    , &
         precav    , &
         precmn    , &
         preccountavl    , &
         preccountav    , &
         preccountmn    , &
         prec_prcavl    , &
         prec_prcav    , &
         prec_prcmn    , &
         cloudcountavl  , &
         cloudcountav    , &
         cloudcountmn    , &
         raincountavl    , &
         raincountav    , &
         raincountmn    , &
         Nrrainavl    , &
         Nrrainav    , &
         Nrrainmn    , &
         qravl    , &
         qrav      , &
         qrmn      , &
         Dvravl    , &
         Dvrav    , &
         Dvrmn    )
  end subroutine exitbulkmicrostat

!------------------------------------------------------------------------------!

end module
