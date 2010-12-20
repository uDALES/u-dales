!> \file modgenstat.f90
!!  Genstat calculates slab averages of several variables

!>
!!  Genstat calculates slab averages of several variables
!>
!!  Written to fields.expnr, moments.expnr and flux1.expnr and flux2.expnr
!! If netcdf is true, this module leads the profiles.expnr.nc output
!!  \author Hans Cuijpers, K.N.M.I.
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  \todo Documentation
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



module modgenstat

    !-----------------------------------------------------------------|
    !                                                                 |
    !*** *stattend*  calculates generic slabaveraged statistics       |
    !                                                                 |
    !      Pier Siebesma   K.N.M.I.     12/05/1995                    |
    !      Hans Cuijpers   I.M.A.U.     21/06/1995                    |
    !                                                                 |
    !     purpose.                                                    |
    !     --------                                                    |
    !                                                                 |
    !     stattend.f calculates:                                      |
    !                                                                 |
    !     * time averaged fieldsof theta_l, theta, theta_v, qt, qv   |
    !        ql, u and v                                              |
    !     * time averaged tendencies  of theta_l, theta, qt, qv, ql,  |
    !        u and v                                                  |
    !     * time averaged turbulent fluxes of theta_l, theta_v,       |
    !       theta, qt, qv, ql, u and v                                |
    !     * variances of qt, w, u, theta_l and theta_v                |
    !     * skewness of qt and w                                      |
    !*** *genstat*  calculates timeseries of several variables       |
    !                                                                 |
    !____________________SETTINGS_AND_SWITCHES________________________|
    !                     IN &NAMTIMESTAT                             |
    !                                                                 |
    !    dtav           SAMPLING INTERVAL                             |
    !                                                                 |
    !    timeav         INTERVAL OF WRITING                           |
    !                                                                 |
    !    lstat      SWITCH TO ENABLE TIMESERIES                       |
    !-----------------------------------------------------------------|
  use modglobal, only : longint

implicit none
! private
PUBLIC :: initgenstat, genstat, exitgenstat
save

!NetCDF variables
  integer :: nvar = 37
  integer :: ncid,nrec = 0
  character(80) :: fname = 'profiles.xxx.nc'
  character(80),allocatable, dimension(:,:) :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav, timeav
  integer(kind=longint) :: idtav,itimeav,tnext,tnextwrite
  logical :: lstat= .false. ! switch for conditional sampling cloud (on/off)
  integer :: nsamples
!     ----  total fields  ---

  real, allocatable  :: umn   (:)       ,vmn   (:)
  real, allocatable  :: thlmn (:)       ,thvmn (:)
  real, allocatable  :: qtmn  (:)       ,qlmn  (:),  qlhmn(:),cfracmn(:)

! real, allocatable  ::     --- fluxes (resolved, subgrid and total) ---
  real, allocatable  :: wtlsmn (:),wtlrmn (:),wtltmn(:)
  real, allocatable  :: wtvsmn (:),wtvrmn (:),wtvtmn(:)
  real, allocatable  :: wqlsmn (:),wqlrmn (:),wqltmn(:)
  real, allocatable  :: wqtsmn (:),wqtrmn (:),wqttmn(:)
  real, allocatable  :: wsvsmn (:,:),wsvrmn(:,:),wsvtmn(:,:)
  real, allocatable  :: uwtmn  (:),vwtmn  (:) !total    uw, vw
  real, allocatable  :: uwsmn  (:),vwsmn  (:) !resolved uw, vw
  real, allocatable  :: uwrmn  (:),vwrmn  (:) !subgrid  uw, vw
! real, allocatable  ::     --- various moments ---

!   real, allocatable  :: rmn        (:), r2mn   (:), r3mn (:), rhmn (:)
  real, allocatable  :: w2mn       (:), skewmn (:)
  real, allocatable  :: w2submn    (:)
  real, allocatable  :: u2mn       (:), v2mn  (:),     qt2mn(:)
  real, allocatable  :: thl2mn     (:), thv2mn(:),     th2mn(:),     ql2mn(:)
!   real, allocatable  :: qs2mn      (:), qsmn  (:)
  real, allocatable  :: svmmn(:,:),svptmn(:,:),svplsmn(:,:),svpmn(:,:)
  real, allocatable  :: sv2mn(:,:)

 real, allocatable :: umav (:)     ! slab averaged ql_0    at full level
 real, allocatable :: vmav (:)     ! slab averaged ql_0    at full level
 real, allocatable :: thlmav (:)     ! slab averaged ql_0    at full level
 real, allocatable :: thmav (:)     ! slab averaged ql_0    at full level
 real, allocatable :: qtmav (:)     ! slab averaged ql_0    at full level
 real, allocatable :: qlmav (:)     ! slab averaged ql_0    at full level
 real, allocatable :: cfracav (:)     ! slab averaged ql_0    at full level
 real, allocatable :: svmav (:,:)     ! slab averaged ql_0    at full level
  real, allocatable :: svpav(:,:)                  !  slab average total tendency of sv(n)
  real, allocatable :: svptav(:,:)                 !  slab average tendency of sv(n) due to turb.

  real, allocatable :: uptav(:)                      !  slab averaged tendency for u
  real, allocatable :: vptav(:)                      !  slab averaged tendency for v

 real, allocatable :: thptav(:)     ! slab averaged turbulence tendency of theta
 real, allocatable :: qlptav(:)     ! slab averaged turbulence tendency of q_liq
 real, allocatable :: uwtot (:)     ! slab averaged tot w-u flux at half levels
 real, allocatable :: vwtot (:)     ! slab averaged tot w-v flux at half levels
 real, allocatable :: uwsub (:)     ! slab averaged sgs w-u flux at half levels
 real, allocatable :: vwsub (:)     ! slab averaged sgs w-v flux at half levels
 real, allocatable :: uwres (:)     ! slab averaged res w-u flux at half levels
 real, allocatable :: vwres (:)     ! slab averaged res w-v flux at half levels


  real, allocatable :: thvhav(:)
  real, allocatable :: th0av(:)
 real, allocatable :: wtlsub(:)     ! slab averaged sub w-theta_l flux at half levels
 real, allocatable :: wtlres(:)     ! slab averaged res w-theta_l flux at half levels
 real, allocatable :: wtltot(:)     ! slab averaged tot w-theta_l flux at half levels

 real, allocatable :: wqtsub(:)    ! slab averaged sub w-qtot    flux at half levels
 real, allocatable :: wqtres(:)    ! slab averaged res w-qtot    flux at half levels
 real, allocatable :: wqttot(:)    ! slab averaged tot w-qtot    flux at half levels

 real, allocatable :: wqlsub(:)    ! slab averaged sub w-ql      flux at half levels
 real, allocatable :: wqlres(:)    ! slab averaged tot w-ql      flux at half levels
 real, allocatable :: wqltot(:)    ! slab averaged tot w-ql      flux at half levels

 real, allocatable :: wtvtot(:)    !  slab averaged total wthv-flux
 real, allocatable :: wtvres(:)    !  slab averaged res   wthv-flux
 real, allocatable :: wtvsub(:)    !  slab averaged sub   wthv-flux

 real, allocatable :: wsvsub(:,:)! slab averaged sub w-sv(n)  flux
 real, allocatable :: wsvres(:,:)! slab averaged res w-sv(n)  flux
 real, allocatable :: wsvtot(:,:)! slab averaged tot w-sv(n)  flux

 real, allocatable :: cszmn(:)    ! Smagorinsky constant
 real, allocatable :: cszav(:)    ! Smagorinsky constant
contains

  subroutine initgenstat
    use modmpi,    only : myid,mpierr, comm3d,my_real, mpi_logical
    use modglobal, only : dtmax, kmax,k1, nsv,ifnamopt,fname_options, ifoutput, cexpnr,dtav_glob,timeav_glob,ladaptive,dt_lim,btime,tres
    use modstat_nc, only : lnetcdf, open_nc,define_nc,redefine_nc,ncinfo,writestat_dims_nc
    use modsurfdata, only : isurf,ksoilmax

    implicit none

    integer n, ierr
    character(40) :: name
    character(3) :: csvname
    namelist/NAMGENSTAT/ &
    dtav,timeav,lstat

    dtav=dtav_glob;timeav=timeav_glob

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMGENSTAT,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMGENSTAT'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMGENSTAT'
      endif
      write(6 ,NAMGENSTAT)
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

    allocate(umn(k1)       ,vmn   (k1))
    allocate(thlmn (k1)       ,thvmn (k1))
    allocate(qtmn  (k1)       ,qlmn  (k1),  qlhmn(k1),cfracmn(k1))
    allocate(wtlsmn (k1),wtlrmn (k1),wtltmn(k1))
    allocate(wtvsmn (k1),wtvrmn (k1),wtvtmn(k1))
    allocate(wqlsmn (k1),wqlrmn (k1),wqltmn(k1))
    allocate(wqtsmn (k1),wqtrmn (k1),wqttmn(k1))
    allocate(wsvsmn (k1,nsv),wsvrmn(k1,nsv),wsvtmn(k1,nsv))
    allocate(uwtmn  (k1),vwtmn  (k1),uwrmn  (k1),  vwrmn(k1),uwsmn(k1),vwsmn(k1))
!     allocate(rmn(k1), r2mn   (k1), r3mn (k1), rhmn (k1))
    allocate(w2mn       (k1), skewmn (k1))
    allocate(w2submn    (k1))
    allocate(u2mn       (k1), v2mn  (k1),     qt2mn(k1))
    allocate(thl2mn     (k1), thv2mn(k1),     th2mn(k1),     ql2mn(k1))
!     allocate(qs2mn      (k1), qsmn  (k1))
    allocate(svmmn(k1,nsv),svptmn(k1,nsv),svplsmn(k1,nsv),svpmn(k1,nsv))
    allocate(sv2mn(k1,nsv))

    allocate(umav (k1))
    allocate(vmav (k1))
    allocate(thlmav (k1))
    allocate(thmav (k1))
    allocate(qtmav (k1))
    allocate(qlmav (k1))
    allocate(cfracav(k1))
    allocate(svmav (k1,nsv))
    allocate(uptav(k1))
    allocate(vptav(k1))

    allocate(thptav(k1))
    allocate(qlptav(k1))
    allocate(uwtot (k1))
    allocate(vwtot (k1))
    allocate(uwsub (k1))
    allocate(vwsub (k1))
    allocate(uwres (k1))
    allocate(vwres (k1))
    allocate(wtlsub(k1))
    allocate(wtlres(k1))
    allocate(wtltot(k1))
    allocate(wqtsub(k1))
    allocate(wqtres(k1))
    allocate(wqttot(k1))
    allocate(wqlsub(k1))
    allocate(wqlres(k1))
    allocate(wqltot(k1))
    allocate(wtvtot(k1))
    allocate(wtvres(k1))
    allocate(wtvsub(k1))
    allocate(wsvsub(k1,nsv))
    allocate(wsvres(k1,nsv))
    allocate(wsvtot(k1,nsv))
    allocate(thvhav(k1))
    allocate(th0av(k1))
    allocate(svptav(k1,nsv))
    allocate(svpav(k1,nsv))

    allocate(cszmn(k1), cszav(k1))

      umn      = 0.
      vmn      = 0.
      thlmn    = 0.
      thvmn    = 0.
      qtmn     = 0.
      qlmn     = 0.
      qlhmn    = 0.
      cfracmn  = 0.

      wtlsmn =  0.
      wtlrmn =  0.
      wtltmn =  0.

      wtvsmn =  0.
      wtvrmn =  0.
      wtvtmn =  0.

      wqtsmn =  0.
      wqtrmn =  0.
      wqttmn =  0.

      wqlsmn =  0.
      wqlrmn =  0.
      wqltmn =  0.

      uwtmn  = 0.
      vwtmn  = 0.
      uwrmn  = 0.
      vwrmn  = 0.
      uwsmn  = 0.
      vwsmn  = 0.


      u2mn     = 0.
      v2mn     = 0.
      w2mn     = 0.
      w2submn  = 0.
      skewmn   = 0.
      qt2mn    = 0.
      thl2mn   = 0.
      thv2mn   = 0.
      th2mn    = 0.
      ql2mn    = 0.
!       qs2mn    = 0.
!       qsmn     = 0.
!       rhmn     = 0.
!       rmn      = 0.
!       r2mn     = 0.
!       r3mn     = 0.

      svmmn   = 0.
      svpmn   = 0.
      svpav   = 0.
!       svplsmn = 0.
      svptmn  = 0.
      svptav  = 0.

      sv2mn = 0.

      wsvsmn = 0.
      wsvrmn = 0.
      wsvtmn = 0.

      cszav = 0.
      cszmn = 0.

      if(myid==0)then
        open (ifoutput,file='field.'//cexpnr,status='replace')
        close (ifoutput)
        open (ifoutput,file='flux1.'//cexpnr,status='replace')
        close (ifoutput)
        open (ifoutput,file='flux2.'//cexpnr,status='replace')
        close (ifoutput)
        open (ifoutput,file='moments.'//cexpnr,status='replace')
        close (ifoutput)
        do n=1,nsv
            name = 'svnnnfld.'//cexpnr
            write (name(3:5),'(i3.3)') n
            open (ifoutput,file=name,status='replace')
            close (ifoutput)
        end do
        do n=1,nsv
            name = 'svnnnflx.'//cexpnr
            write (name(3:5),'(i3.3)') n
            open (ifoutput,file=name,status='replace')
            close (ifoutput)
        end do
      if (lnetcdf) then
        fname(10:12) = cexpnr
        nvar = nvar + 7*nsv
        allocate(ncname(nvar,4))
        call ncinfo(tncname(1,:),'time','Time','s','time')
        call ncinfo(ncname( 1,:),'dn0','Base-state density','kg/m^3','tt')
        call ncinfo(ncname( 2,:),'presh','Pressure at cell center','Pa','tt')
        call ncinfo(ncname( 3,:),'u','West-East velocity','m/s','tt')
        call ncinfo(ncname( 4,:),'v','South-North velocity','m/s','tt')
        call ncinfo(ncname( 5,:),'thl','Liquid water potential temperature','K','tt')
        call ncinfo(ncname( 6,:),'thv','Virtual potential temperature','K','tt')
        call ncinfo(ncname( 7,:),'qt','Total water mixing ratio','kg/kg','tt')
        call ncinfo(ncname( 8,:),'ql','Liquid water mixing ratio','kg/kg','tt')
        call ncinfo(ncname( 9,:),'wtls','SFS-Theta_l flux','Km/s','mt')
        call ncinfo(ncname(10,:),'wtlr','Resolved Theta_l flux','Km/s','mt')
        call ncinfo(ncname(11,:),'wtlt','Total Theta_l flux','Km/s','mt')
        call ncinfo(ncname(12,:),'wtvs','SFS-buoyancy flux','Km/s','mt')
        call ncinfo(ncname(13,:),'wtvr','Resolved buoyancy flux','Km/s','mt')
        call ncinfo(ncname(14,:),'wtvt','Total buoyancy flux','Km/s','mt')
        call ncinfo(ncname(15,:),'wqts','SFS-moisture flux','kg/kg m/s','mt')
        call ncinfo(ncname(16,:),'wqtr','Resolved moisture flux','kg/kg m/s','mt')
        call ncinfo(ncname(17,:),'wqtt','Total moisture flux','kg/kg m/s','mt')
        call ncinfo(ncname(18,:),'wqls','SFS-liquid water flux','kg/kg m/s','mt')
        call ncinfo(ncname(19,:),'wqlr','Resolved liquid water flux','kg/kg m/s','mt')
        call ncinfo(ncname(20,:),'wqlt','Total liquid water flux','kg/kg m/s','mt')
        call ncinfo(ncname(21,:),'uws','SFS-momentum flux (uw)','m^2/s^2','mt')
        call ncinfo(ncname(22,:),'uwr','Resolved momentum flux (uw)','m^2/s^2','mt')
        call ncinfo(ncname(23,:),'uwt','Total momentum flux (vw)','m^2/s^2','mt')
        call ncinfo(ncname(24,:),'vws','SFS-momentum flux (vw)','m^2/s^2','mt')
        call ncinfo(ncname(25,:),'vwr','Resolved momentum flux (vw)','m^2/s^2','mt')
        call ncinfo(ncname(26,:),'vwt','Total momentum flux (vw)','m^2/s^2','mt')
        call ncinfo(ncname(27,:),'w2s','SFS-TKE','m^2/s^2','mt')
        call ncinfo(ncname(28,:),'w2r','Resolved vertical velocity variance','m^2/s^2','mt')
        !call ncinfo(ncname(29,:),'w2t','Total vertical velocity variance','m^2/s^2','mt')
        call ncinfo(ncname(29,:),'skew','vertical velocity skewness','-','mt')
        call ncinfo(ncname(30,:),'u2r','Resolved horizontal velocity variance (u)','m^2/s^2','tt')
        call ncinfo(ncname(31,:),'v2r','Resolved horizontal velocity variance (v)','m^2/s^2','tt')
        call ncinfo(ncname(32,:),'thl2r','Resolved theta_l variance','K^2','tt')
        call ncinfo(ncname(33,:),'thv2r','Resolved buoyancy variance','K^2','tt')
        call ncinfo(ncname(34,:),'th2r','Resolved theta variance','K^2','tt')
        call ncinfo(ncname(35,:),'qt2r','Resolved total water variance','(kg/kg)^2','tt')
        call ncinfo(ncname(36,:),'ql2r','Resolved liquid water variance','(kg/kg)^2','tt')
        call ncinfo(ncname(37,:),'cs','Smagorinsky constant','-','tt')
        do n=1,nsv
          write (csvname(1:3),'(i3.3)') n
          call ncinfo(ncname(37+7*(n-1),:),'sv'//csvname,'Scalar '//csvname//' mixing ratio','(kg/kg)','tt')
          call ncinfo(ncname(37+7*(n-1)+1,:),'svp'//csvname,'Scalar '//csvname//' tendency','(kg/kg/s)','tt')
          call ncinfo(ncname(37+7*(n-1)+2,:),'svpt'//csvname,'Scalar '//csvname//' turbulence tendency','(kg/kg/s)','tt')
          call ncinfo(ncname(37+7*(n-1)+3,:),'sv'//csvname//'2r','Resolved scalar '//csvname//' variance','(kg/kg)^2','tt')
          call ncinfo(ncname(37+7*(n-1)+4,:),'wsv'//csvname//'s','SFS scalar '//csvname//' flux','kg/kg m/s','mt')
          call ncinfo(ncname(37+7*(n-1)+5,:),'wsv'//csvname//'r','Resolved scalar '//csvname//' flux','kg/kg m/s','mt')
          call ncinfo(ncname(37+7*(n-1)+6,:),'wsv'//csvname//'t','Total scalar '//csvname//' flux','kg/kg m/s','mt')
        end do

        if (isurf==1) then
          call open_nc(fname,  ncid,n3=kmax,ns=ksoilmax)
        else
          call open_nc(fname,  ncid,n3=kmax,ns=ksoilmax)
        endif
        call define_nc( ncid, 1, tncname)
        call writestat_dims_nc(ncid)
        call redefine_nc(ncid)
        call define_nc( ncid, NVar, ncname)
      end if


      end if

  end subroutine initgenstat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine genstat

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
      call do_genstat
    end if
    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+itimeav
      call writestat
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
  end subroutine genstat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine do_genstat

    use modfields, only : u0,v0,w0,um,vm,wm,qtm,thlm,thl0,qt0,qt0h, &
                          ql0,ql0h,thl0h,thv0h,sv0, svm, e12m,exnf,exnh
    use modsurfdata,only: thls,qts,svs,ustar,thlflux,qtflux,svflux
    use modsubgriddata,only : ekm, ekh, csz
    use modglobal, only : i1,ih,j1,jh,k1,kmax,nsv,dzf,dzh,rlv,rv,rd,cp, &
                          rslabs,cu,cv,iadv_thl,iadv_kappa,eps1,dxi,dyi
    use modmpi,    only : nprocs,comm3d,nprocs,my_real,mpi_sum,mpierr,slabsum
    implicit none



    real cthl,cqt,den

    real,allocatable, dimension(:) :: &
        qlhavl , & ! slab averaged ql_0 at half level &
        u2avl    , &
        v2avl    , &
        w2avl    , &
        w3avl    , &
        w2subavl , &
        qt2avl   , &
        thl2avl  , &
        thv2avl  , &
        ql2avl   , &
        th2avl
    real,allocatable, dimension(:,:) :: &
        wsvsubl,&   ! slab averaged sub w-sv(n)  flux &
        wsvresl,&   ! slab averaged res w-sv(n)  flux &
        sv2avl,&
        sv2av

    real,allocatable, dimension(:) :: wqlsubl
    real,allocatable, dimension(:) :: wqlresl

    real,allocatable, dimension(:):: wtlsubl
    real,allocatable, dimension(:):: wtlresl

    real,allocatable, dimension(:):: wqtsubl
    real,allocatable, dimension(:):: wqtresl

    real,allocatable, dimension(:):: wtvsubl
    real,allocatable, dimension(:) ::wtvresl

    real,allocatable, dimension(:):: cfracavl ! cloudfraction    at full level


    real,allocatable, dimension(:):: qlptavl   ! slab averaged turbulence tendency of q_liq
    real,allocatable, dimension(:):: uwsubl
    real,allocatable, dimension(:):: vwsubl
    real,allocatable, dimension(:):: uwresl
    real,allocatable, dimension(:):: vwresl
    real,allocatable, dimension(:):: qlhav
    real,allocatable, dimension(:):: u2av    , &
              v2av    , &
              w2av    , &
              w3av    , &
              w2subav , &
              qt2av   , &
              thl2av  , &
              thv2av  , &
              th2av   , &
              ql2av
    real,allocatable, dimension(:,:,:)::  thv0
    real,allocatable, dimension(:)::   thvmav
    real ,allocatable, dimension(:,:,:):: sv0h

    integer i, j, k, n, km
    real    tsurf, qsat, c1, c2
    real    qs0h, t0h, ekhalf, euhalf, evhalf
    real    wtls, wtlr, wqts, wqtr, wqls, wqlr, wtvs, wtvr
    real    uws,vws,uwr,vwr
    real    upcu, vpcv
    real    qls
    allocate( &
        qlhavl (k1), & ! slab averaged ql_0 at half level &
        wsvsubl(k1,nsv),&   ! slab averaged sub w-sv(n)  flux &
        wsvresl(k1,nsv),&   ! slab averaged res w-sv(n)  flux &
        u2avl    (k1), &
        v2avl    (k1), &
        w2avl    (k1), &
        w3avl    (k1), &
        w2subavl (k1), &
        qt2avl   (k1), &
        thl2avl  (k1), &
        thv2avl  (k1), &
        th2avl   (k1))
    allocate( &
        ql2avl   (k1), &
        sv2avl   (k1,nsv))


    allocate( wqlsubl    (k1))
    allocate( wqlresl    (k1))

    allocate( wtlsubl    (k1))
    allocate( wtlresl    (k1))

    allocate( wqtsubl    (k1))
    allocate( wqtresl    (k1))

    allocate( wtvsubl    (k1))
    allocate( wtvresl    (k1))

    allocate( cfracavl(k1))  ! slab averaged cloud fraction


    allocate( qlptavl(k1))   ! slab averaged turbulence tendency of q_liq
    allocate( uwsubl(k1))
    allocate( vwsubl(k1))
    allocate( uwresl(k1))
    allocate( vwresl(k1))
    allocate( qlhav(k1))
    allocate( u2av    (k1), &
              v2av    (k1), &
              w2av    (k1), &
              w3av    (k1), &
              w2subav (k1), &
              qt2av   (k1), &
              thl2av  (k1), &
              thv2av  (k1), &
              th2av   (k1), &
              ql2av   (k1), &
              sv2av   (k1,nsv))
    allocate(thv0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thvmav(k1))
    allocate(sv0h(2-ih:i1+ih,2-jh:j1+jh,k1))



  !-----------------------------------------------------------------------
  !     1.    INITIALISE LOCAL CONSTANTS
  !     --    --------------------------

  !     --------------------------------------------------------
  !     3.0    RESET ARRAYS FOR SLAB AVERAGES
  !     ---    ------------------------------
  !     --------------------------------------------------------
    qlhavl      = 0.0
    cfracavl    = 0.0


    qlptavl     = 0.0

    wqlsubl     = 0.0
    wqlresl     = 0.0
    wqltot      = 0.0

    wtlsubl     = 0.0
    wtlresl     = 0.0
    wtltot      = 0.0

    wqtsubl     = 0.0
    wqtresl     = 0.0
    wqttot      = 0.0

    wtvsubl     = 0.0
    wtvresl     = 0.0
    wtvtot      = 0.0


    wsvsubl = 0.
    wsvresl = 0.
    sv2avl  = 0.

    uwresl  = 0.
    vwresl  = 0.
    uwtot   = 0.
    uwsubl  = 0.
    vwsubl  = 0.
    vwtot   = 0.

    u2avl     = 0.0
    v2avl     = 0.0
    w2avl     = 0.0
    w3avl     = 0.0
    w2subavl  = 0.0
    qt2avl    = 0.0
    thl2avl   = 0.0
    thv2avl   = 0.0
    th2avl    = 0.0
    ql2avl    = 0.0
    thvmav    = 0.0

    sv2av   = 0.0

    umav = 0.0
    vmav = 0.0
    thlmav = 0.0
    thmav  = 0.0
    qtmav  = 0.0
    qlmav  = 0.0
    cfracav= 0.0
    svmav = 0.

    cszav = 0.

    do  k=1,k1
      do  j=2,j1
        do  i=2,i1
          thv0(i,j,k) = (thl0(i,j,k)+rlv*ql0(i,j,k)/(cp*exnf(k))) &
                        *(1+(rv/rd-1)*qt0(i,j,k)-rv/rd*ql0(i,j,k))
        enddo
      enddo
    enddo

    do k=1,k1
      cfracavl(k)    = cfracavl(k)+count(ql0(2:i1,2:j1,k)>0)
    end do

    call MPI_ALLREDUCE(cfracavl,cfracav,k1,MY_REAL,MPI_SUM,comm3d,mpierr)

    call slabsum(umav  ,1,k1,um  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(vmav  ,1,k1,vm  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(thlmav,1,k1,thlm,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(qtmav ,1,k1,qtm ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(qlmav ,1,k1,ql0 ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(thvmav,1,k1,thv0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)

    umav  = umav  /rslabs + cu
    vmav  = vmav  /rslabs + cv
    thlmav = thlmav/rslabs
    qtmav = qtmav /rslabs
    qlmav = qlmav /rslabs
    cfracav = cfracav / rslabs
    thmav  = thlmav + (rlv/cp)*qlmav/exnf
    thvmav = thvmav/rslabs

    cszav  = csz
  !

    do n=1,nsv
      call slabsum(svmav(1,n),1,k1,svm(1,1,1,n),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    enddo
    svmav = svmav/rslabs
  !------------------------------------------------------------------
  !     4     CALCULATE SLAB AVERAGED OF FLUXES AND SEVERAL MOMENTS
  !     -------------------------------------------------------------
!         4.1 special treatment for lowest level
  !     -------------------------------------------------

    qls   = 0.0 ! hj: no liquid water at the surface
    tsurf = thls*exnh(1)+(rlv/cp)*qls
    qsat  = qts - qls
    if (qls==0) then
      c1  = 1.+(rv/rd-1)*qts
      c2  = (rv/rd-1)
    else
      c1    = (1.-qts+rv/rd*qsat*(1.+rd/(rv*tsurf))) &
                /(1.+rlv/(rv*tsurf)*rlv/(cp*tsurf)*qsat)
      c2    = c1*rlv/(tsurf*cp)-1.
    end if
    den   = 1. + (rlv**2)*qsat/(rv*cp*(tsurf**2))
    cthl  = (exnh(1)*cp/rlv)*((1-den)/den)
    cqt   = 1./den
    do j=2,j1
    do i=2,i1
      qlhavl(1) = qlhavl(1) + ql0h(i,j,1)
  !     thv(1) = thlm(i,j,1) * (1.+(rv/rd-1)*qtm(i,j,1))

      wtlsubl(1) = wtlsubl(1) + thlflux(i,j)
      wqtsubl(1) = wqtsubl(1) + qtflux (i,j)
      wqlsubl(1) = 0
      wtvsubl(1) = wtvsubl(1) + ( c1*thlflux(i,j)+c2*thls*qtflux(i,j) ) !hj: thv0 replaced by thls

      !Momentum flux
      if (abs(um(i,j,1)+cu)<eps1) then
        upcu = sign(eps1,um(i,j,1)+cu)
      else
        upcu = um(i,j,1)+cu
      end if
      uwsubl(1) = uwsubl(1) - ( 0.5*( ustar(i,j)+ustar(i-1,j) ) )**2  * &
                upcu/sqrt(upcu**2  + &
          ((vm(i,j,1)+vm(i-1,j,1)+vm(i,j+1,1)+vm(i-1,j+1,1))/4.+cv)**2)

      if (abs(vm(i,j,1)+cv)<eps1) then
        vpcv = sign(eps1,vm(i,j,1)+cv)
      else
        vpcv = vm(i,j,1)+cv
      end if
      vwsubl(1) = vwsubl(1) - ( 0.5*( ustar(i,j)+ustar(i,j-1) ) )**2  * &
                vpcv/sqrt(vpcv**2  + &
          ((um(i,j,1)+um(i+1,j,1)+um(i,j-1,1)+um(i+1,j-1,1))/4.+cu)**2)


      !Higher order moments
      u2avl    (1) = u2avl    (1) + (um (i,j,1)+cu - umav(1))**2
      v2avl    (1) = v2avl    (1) + (vm (i,j,1)+cv - vmav(1))**2
      w2avl    (1) = w2avl    (1) + (wm  (i,j,1)**2)
      w3avl    (1) = w3avl    (1) + (wm  (i,j,1)**3)
      w2subavl (1) = w2subavl (1) + (e12m(i,j,1)**2)
      qt2avl   (1) = qt2avl   (1) + (qtm (i,j,1) - qtmav (1))**2
      thl2avl  (1) = thl2avl  (1) + (thlm(i,j,1) - thlmav(1))**2
      thv2avl  (1) = thv2avl  (1) + (thv0(i,j,1) - thvmav(1))**2
      th2avl   (1) = th2avl   (1) + (thlm(i,j,1) - thmav (1))**2
      ql2avl   (1) = ql2avl   (1) + (ql0(i,j,1)  - qlmav (1))**2
!       qs2avl   (1) = qs2avl   (1) + qs0**2
!       qsavl    (1) = qsavl    (1) + qs0
!       rhavl    (1) = rhavl    (1) + qtm (i,j,1)/qs0
!       ravl     (1) = ravl     (1) + ((qt0(i,j,1) - qs0))
!       r2avl    (1) = r2avl    (1) + ((qt0(i,j,1) - qs0)**2)
!       r3avl    (1) = r3avl    (1) + ((qt0(i,j,1) - qs0)**3)

      do n=1,nsv
        wsvsubl(1,n) = wsvsubl(1,n) + svflux(i,j,n)
        sv2avl(1,n)  = sv2avl(1,n) + (svm(i,j,1,n)-svmav(1,n))**2
      end do
    end do
    end do

  !      --------------------------
  !      4.2 higher levels
  !      --------------------------

    do j=2,j1
    do i=2,i1

  !     --------------------------------------------------------
  !      Calculate half level fields for thl and qt consistent
  !      with the used advection scheme ( kappa or cent. diff.)
  !     ---------------

  !     --------------------------------------------------------

      do k=2,kmax
        km = k-1


    !     ------------------------------------------------------
    !     calculate ql and thv at time t0 at full and half level
    !      ----------------------------------------------------
        qlhavl(k) = qlhavl(k)  + ql0h(i,j,k)

    !     -----------------------------------------------------------
    !     calculate prefactors for subgrid wthv and wql fluxes
    !      at half levels
    !     -----------------------------------------------------------
        qs0h  =  (qt0h(i,j,k) - ql0h(i,j,k))
        t0h   =  exnh(k)*thl0h(i,j,k) + (rlv/cp)*ql0h(i,j,k)

        den   = 1. + (rlv**2)*qs0h/(rv*cp*(t0h**2))
        cthl  = (exnh(k)*cp/rlv)*((1-den)/den)
        cqt   =  (1./den)
        if (ql0h(i,j,k)>0) then
          c1    = (1.-qt0h(i,j,k)+rv/rd*qs0h &
                  * (1.+rd/rv*rlv/(rd*t0h)))/den
          c2    =  c1*rlv/(t0h*cp)-1.
        else
          c1 = 1. + (rv/rd-1)*qt0h(i,j,k)
          c2 = (rv/rd-1)
        end if

    !     -----------------------------------------------------------
    !     calculate resolved and subgrid fluxes at half levels
    !     -----------------------------------------------------------

        ekhalf  = (ekh(i,j,k)*dzf(km)+ekh(i,j,km)*dzf(k))/(2*dzh(k))
        euhalf = ( dzf(km) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i-1,j,km) ) ) / &
                    ( 4.   * dzh(k) )
        evhalf = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,j-1,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i,j-1,km) ) ) / &
                    ( 4.   * dzh(k) )

        wtls    = -ekhalf*(thl0(i,j,k)-thl0(i,j,km))/dzh(k)
        wtlr    = w0(i,j,k)*thl0h(i,j,k)

        wqts    = -ekhalf*(qt0(i,j,k)-qt0(i,j,km))/dzh(k)
        wqtr    = w0(i,j,k)*qt0h(i,j,k)

        wqls    = cthl*wtls+ cqt*wqts
        wqlr    = w0(i,j,k)*ql0h(i,j,k)

        wtvs    = c1*wtls + c2*thl0h(i,j,k)*wqts
        wtvr    = w0(i,j,k)*thv0h(i,j,k)

        uwr     = (w0(i,j,k)+w0(i-1,j,k)) &
                  *((u0(i,j,k-1)+cu)*dzf(k)+(u0(i,j,k)+cu)*dzf(k-1))/(4*dzh(k))
        vwr     = (w0(i,j,k)+w0(i,j-1,k)) &
                  *((v0(i,j,k-1)+cv)*dzf(k)+(v0(i,j,k)+cv)*dzf(k-1))/(4*dzh(k))
        uws     = -euhalf &
                  *((u0(i,j,k)-u0(i,j,k-1))/dzh(k)+(w0(i,j,k)-w0(i-1,j,k))*dxi)
        vws     = -evhalf &
                  *((v0(i,j,k)-v0(i,j,k-1))/dzh(k)+(w0(i,j,k)-w0(i,j-1,k))*dyi)


        if (ql0h(i,j,k)>0) then
          wqlsubl(k) = wqlsubl(k) + wqls
        end if

        wqlresl(k) = wqlresl(k) + wqlr

        wtlsubl(k) = wtlsubl(k) + wtls
        wtlresl(k) = wtlresl(k) + wtlr

        wtvsubl(k) = wtvsubl(k) + wtvs
        wtvresl(k) = wtvresl(k) + wtvr

        wqtsubl(k) = wqtsubl(k) + wqts
        wqtresl(k) = wqtresl(k) + wqtr

        uwresl(k) = uwresl(k) + uwr
        vwresl(k) = vwresl(k) + vwr
        uwsubl(k) = uwsubl(k) + uws
        vwsubl(k) = vwsubl(k) + vws
    !     -----------------------------------------------------------
    !     calculate various moments
    !     -----------------------------------------------------------

        u2avl    (k) = u2avl    (k) + (um (i,j,k)+cu - umav(k))**2
        v2avl    (k) = v2avl    (k) + (vm (i,j,k)+cv - vmav(k))**2
        w2avl    (k) = w2avl    (k) + (wm  (i,j,k)**2)
        w3avl    (k) = w3avl    (k) + (wm  (i,j,k)**3)
        w2subavl (k) = w2subavl (k) + (e12m(i,j,k)**2)
        qt2avl   (k) = qt2avl   (k) + (qtm (i,j,k) - qtmav (k))**2
        thl2avl  (k) = thl2avl  (k) + (thlm(i,j,k) - thlmav(k))**2
        thv2avl  (k) = thv2avl  (k) + (thv0(i,j,k) - thvmav(k))**2
        th2avl   (k) = th2avl   (k) + (thlm(i,j,k) - thmav (k))**2 !thlm, no thm !?!
        ql2avl   (k) = ql2avl   (k) + (ql0(i,j,k)  - qlmav (k))**2
!         qs2avl   (k) = qs2avl   (k) + qs0**2
!         qsavl    (k) = qsavl    (k) + qs0
!         rhavl    (k) = rhavl    (k) + qtm (i,j,k)/qs0
!         ravl     (k) = ravl     (k) + (qt0(i,j,k) - qs0)
!         r2avl    (k) = r2avl    (k) + (qt0(i,j,k) - qs0)**2
!         r3avl    (k) = r3avl    (k) + (qt0(i,j,k) - qs0)**3

      end do
    end do
    end do
  !     -------------------

    do n=1,nsv
      do k=2,kmax
      do j=2,j1
      do i=2,i1
        sv2avl(k,n)  = sv2avl(k,n) + (svm(i,j,k,n)-svmav(k,n))**2
      end do
      end do
      end do
    end do

    do n=1,nsv
      if (iadv_thl==iadv_kappa) then
         call halflev_kappa(sv0(2-ih:i1+ih,2-jh:j1+jh,1:k1,n),sv0h)
      else
        do  j=2,j1
        do  i=2,i1
        do  k=2,k1
          sv0h(i,j,k) = (sv0(i,j,k,n)*dzf(k-1)+sv0(i,j,k-1,n)*dzf(k))/(2*dzh(k))
        enddo
        enddo
        enddo
        sv0h(2:i1,2:j1,1) = svs(n)

      end if
      do  j=2,j1
      do  i=2,i1
      do  k=2,kmax
        wsvresl(k,n) = wsvresl(k,n) + w0(i,j,k)*sv0h(i,j,k)
      end do
      end do
      end do

      do k=2,kmax
        km = k-1
      do j=2,j1
      do i=2,i1
        ekhalf      = (ekh(i,j,k)*dzf(km)+ekh(i,j,km)*dzf(k))/(2*dzh(k))
        wsvsubl(k,n)= wsvsubl(k,n)-ekhalf*(sv0(i,j,k,n)-sv0(i,j,km,n)) &
                                                        /dzh(k)
      end do
      end do
      end do

    end do


  !     -------------------------------
  !     5   CALCULATE MOMENTUM FLUXES
  !     -------------------------------

  !     5.1 special treatment for lowest level
  !     -------------------------------------------------
  !      DEPRECATED


  !     5.2 higher levels by vert. integr. of the mom. tendencies
  !     ---------------------------------------------------------
  !         DEPRECATED

  ! MPI communication
    call MPI_ALLREDUCE(qlhavl, qlhav, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(wqlsubl, wqlsub, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(wqlresl, wqlres, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(wtlsubl, wtlsub, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(wtlresl, wtlres, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(wqtsubl, wqtsub, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(wqtresl, wqtres, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(wtvsubl, wtvsub, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(wtvresl, wtvres, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(uwsubl, uwsub, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(vwsubl, vwsub, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(uwresl, uwres, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(vwresl, vwres, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(u2avl, u2av, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(v2avl, v2av, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(w2avl, w2av, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(w3avl, w3av, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(w2subavl, w2subav, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(qt2avl, qt2av, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(thl2avl, thl2av, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(thv2avl, thv2av, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(th2avl, th2av, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(ql2avl, ql2av, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)
!     call MPI_ALLREDUCE(qs2avl, qs2av, k1,    MY_REAL, &
!                       MPI_SUM, comm3d,mpierr)
!     call MPI_ALLREDUCE(qsavl, qsav, k1,    MY_REAL, &
!                       MPI_SUM, comm3d,mpierr)
!     call MPI_ALLREDUCE(rhavl, rhav, k1,    MY_REAL, &
!                       MPI_SUM, comm3d,mpierr)
!     call MPI_ALLREDUCE(ravl, rav, k1,    MY_REAL, &
!                       MPI_SUM, comm3d,mpierr)
!     call MPI_ALLREDUCE(r2avl, r2av, k1,    MY_REAL, &
!                       MPI_SUM, comm3d,mpierr)
!     call MPI_ALLREDUCE(r3avl, r3av, k1,    MY_REAL, &
!                       MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(qlptavl, qlptav, k1,    MY_REAL, &
                      MPI_SUM, comm3d,mpierr)

    do n=1,nsv
  call MPI_ALLREDUCE(sv2avl(:,n),sv2av(:,n),k1,MY_REAL, &
                        MPI_SUM, comm3d,mpierr)
  call MPI_ALLREDUCE(wsvsubl(:,n),wsvsub(:,n), k1,    MY_REAL, &
      MPI_SUM, comm3d,mpierr)
  call MPI_ALLREDUCE(wsvresl(:,n),wsvres(:,n), k1,    MY_REAL, &
      MPI_SUM, comm3d,mpierr)
    end do

  !     -----------------------------------------------
  !     6   NORMALIZATION OF THE FIELDS AND FLUXES
  !     -----------------------------------------------

      qlhav   = qlhav  /rslabs

      wqlsub  = wqlsub /rslabs
      wqlres  = wqlres /rslabs

      wtlsub  = wtlsub /rslabs
      wtlres  = wtlres /rslabs

      wqtsub  = wqtsub /rslabs
      wqtres  = wqtres /rslabs

      wtvsub  = wtvsub /rslabs
      wtvres  = wtvres /rslabs

      wqttot  = wqtres + wqtsub
      wqltot  = wqlres + wqlsub
      wtltot  = wtlres + wtlsub
      wtvtot  = wtvres + wtvsub


        wsvsub = wsvsub /rslabs
        wsvres = wsvres /rslabs
        wsvtot = wsvsub + wsvres

      uwres    = uwres    /rslabs
      vwres    = vwres    /rslabs
      uwsub    = uwsub    /rslabs
      vwsub    = vwsub    /rslabs
      uwtot    = uwres + uwsub
      vwtot    = vwres + vwsub

      u2av     = u2av     /rslabs
      v2av     = v2av     /rslabs
      w2av     = w2av     /rslabs
      w3av     = w3av     /rslabs
      w2subav  = w2subav  /rslabs
      qt2av    = qt2av    /rslabs
      thl2av   = thl2av   /rslabs
      thv2av   = thv2av   /rslabs
      th2av    = th2av    /rslabs
      ql2av    = ql2av    /rslabs
!       qs2av    = qs2av    /rslabs
!       qsav     = qsav     /rslabs
!       qs2av    = qs2av    - qsav**2
!       rhav     = rhav     /rslabs
!       rav      = rav      /rslabs
!       r2av     = r2av     /rslabs
!       r3av     = r3av     /rslabs


        sv2av = sv2av/rslabs

  !********************************************************************

  !     4.0   ADD SLAB AVERAGES TO TIME MEAN
  !           ------------------------------
      umn    = umn   + umav
      vmn    = vmn   + vmav
      thvmn  = thvmn + thvmav
      thlmn  = thlmn + thlmav
      qtmn   = qtmn  + qtmav
      qlmn   = qlmn  + qlmav
      cfracmn= cfracmn+cfracav
      qlhmn  = qlhmn + qlhav

      wtlsmn = wtlsmn + wtlsub
      wtlrmn = wtlrmn + wtlres
      wtltmn = wtltmn + wtltot
      wtvsmn = wtvsmn + wtvsub
      wtvrmn = wtvrmn + wtvres
      wtvtmn = wtvtmn + wtvtot
      wqtsmn = wqtsmn + wqtsub
      wqtrmn = wqtrmn + wqtres
      wqttmn = wqttmn + wqttot
      wqlsmn = wqlsmn + wqlsub
      wqlrmn = wqlrmn + wqlres
      wqltmn = wqltmn + wqltot
      uwtmn  = uwtmn + uwtot
      vwtmn  = vwtmn + vwtot
      uwrmn  = uwrmn + uwres
      vwrmn  = vwrmn + vwres
      uwsmn  = uwsmn + uwsub
      vwsmn  = vwsmn + vwsub
      u2mn     = u2mn     + u2av
      v2mn     = v2mn     + v2av
      w2mn     = w2mn     + w2av
      w2submn  = w2submn  + w2subav
      qt2mn    = qt2mn    + qt2av
      thl2mn   = thl2mn   + thl2av
      thv2mn   = thv2mn   + thv2av
      th2mn    = th2mn    + th2av
!       ql2mn    = ql2mn    + ql2av
!       qs2mn    = qs2mn    + qs2av
!       qsmn     = qsmn     + qsav
!       rhmn     = rhmn     + rhav
!       rmn      = rmn      + rav
!       r2mn     = r2mn     + r2av
!       r3mn     = r3mn     + r3av

        svmmn   = svmmn  + svmav
        svpmn   = svpmn  + svpav
!         svplsmn = svplsmn+ svplsav
        svptmn  = svptmn + svptav

        sv2mn  = sv2mn + sv2av

        wsvsmn = wsvsmn + wsvsub
        wsvrmn = wsvrmn + wsvres
        wsvtmn = wsvtmn + wsvtot
      skewmn   = skewmn   + w3av/max(w2av**1.5,epsilon(w2av(1)))

      cszmn = cszmn + cszav
      
    deallocate( &
        qlhavl , & ! slab averaged ql_0 at half level &
        wsvsubl,&   ! slab averaged sub w-sv(n)  flux &
        wsvresl,&   ! slab averaged res w-sv(n)  flux &
        u2avl    , &
        v2avl    , &
        w2avl    , &
        w3avl    , &
        w2subavl , &
        qt2avl   , &
        thl2avl  , &
        thv2avl  , &
        th2avl   )
    deallocate( &
        ql2avl   , &
        sv2avl   )


    deallocate( wqlsubl    )
    deallocate( wqlresl    )

    deallocate( wtlsubl    )
    deallocate( wtlresl    )

    deallocate( wqtsubl    )
    deallocate( wqtresl    )

    deallocate( wtvsubl    )
    deallocate( wtvresl    )

    deallocate( cfracavl )


    deallocate( qlptavl)   ! slab averaged turbulence tendency of q_liq
    deallocate( uwsubl)
    deallocate( vwsubl)
    deallocate( uwresl)
    deallocate( vwresl)

    deallocate( qlhav)
    deallocate( u2av    , &
              v2av    , &
              w2av    , &
              w3av    , &
              w2subav , &
              qt2av   , &
              thl2av  , &
              thv2av  , &
              th2av   , &
              ql2av   , &
              sv2av   )
    deallocate(thv0)
    deallocate(thvmav)
    deallocate(sv0h)
  end subroutine do_genstat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writestat
      use modglobal, only : kmax,k1,nsv, zh,zf,rtimee,rlv,cp,cexpnr,ifoutput
      use modfields, only : presf,presh,exnf,exnh,rhof
      use modsubgriddata, only : csz
      use modmpi,    only : myid
      use modstat_nc, only: lnetcdf, writestat_nc
      implicit none


      real,dimension(k1,nvar) :: vars

      real,allocatable, dimension(:) :: tmn, thmn
      integer nsecs, nhrs, nminut,k,n
      real convt, convq
      character(20) :: name
      nsecs   = nint(rtimee)
      nhrs    = int(nsecs/3600)
      nminut  = int(nsecs/60)-nhrs*60
      nsecs   = mod(nsecs,60)
      convt   = 86400.
      convq   = 86400*1000.
      allocate(tmn   (k1), thmn  (k1))

      umn    = umn    /nsamples
      vmn    = vmn    /nsamples
      thvmn  = thvmn  /nsamples
      thlmn  = thlmn  /nsamples
      qtmn   = qtmn   /nsamples
      qlmn   = qlmn   /nsamples
      cfracmn= cfracmn/nsamples
      qlhmn  = qlhmn  /nsamples


      wtlsmn = wtlsmn/nsamples
      wtlrmn = wtlrmn/nsamples
      wtltmn = wtltmn/nsamples

      wqtsmn = wqtsmn/nsamples
      wqtrmn = wqtrmn/nsamples
      wqttmn = wqttmn/nsamples

      wqlsmn = wqlsmn/nsamples
      wqlrmn = wqlrmn/nsamples
      wqltmn = wqltmn/nsamples

      wtvsmn = wtvsmn/nsamples
      wtvrmn = wtvrmn/nsamples
      wtvtmn = wtvtmn/nsamples

      uwtmn  = uwtmn /nsamples
      vwtmn  = vwtmn /nsamples
      uwrmn  = uwrmn /nsamples
      vwrmn  = vwrmn /nsamples
      uwsmn  = uwsmn /nsamples
      vwsmn  = vwsmn /nsamples

!       rmn      = rmn    /nsamples
!       r2mn     = r2mn   /nsamples
!       rhmn     = rhmn   /nsamples
      w2mn     = w2mn   /nsamples
      skewmn   = skewmn /nsamples
      w2submn  = w2submn/nsamples
      qt2mn    = qt2mn  /nsamples
      v2mn     = v2mn   /nsamples
      u2mn     = u2mn   /nsamples
      thl2mn   = thl2mn /nsamples
      thv2mn   = thv2mn /nsamples
      th2mn    = th2mn  /nsamples
      ql2mn    = ql2mn  /nsamples
!       qs2mn    = qs2mn  /nsamples
!       qsmn     = qsmn   /nsamples


        svmmn   = svmmn  /nsamples
        svpmn   = svpmn  /nsamples
!         svplsmn = svplsmn/nsamples
        svptmn  = svptmn /nsamples

        sv2mn = sv2mn/nsamples

        wsvsmn = wsvsmn/nsamples
        wsvrmn = wsvrmn/nsamples
        wsvtmn = wsvtmn/nsamples

        cszmn = cszmn / nsamples


  !     ------------------------------------------
  !     2.0  Construct other time averaged fields
  !     ------------------------------------------


      thmn = thlmn + (rlv/cp)*qlmn/exnf
      tmn  = thmn*exnf

  !     ----------------------
  !     2.0  write the fields
  !           ----------------

    if(myid==0)then
      open (ifoutput,file='field.'//cexpnr,position='append')
      write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
      '#--------------------------------------------------------'      &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '      &
      ,nhrs,':',nminut,':',nsecs      &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
      write (ifoutput,'(A/2A/2A)') &
          '#--------------------------------------------------------' &
          ,'#LEV  HGHT    PRES    TEMP       TH_L     THETA      TH_V     ' &
          ,'  QT_AV      QL_AV      U       V   CLOUD FRACTION  CS' &
          ,'#      (M)    (MB)   (----------- (KELVIN) ---------------)    ' &
          ,'(----(G/KG)------)  (--- (M/S ---)   (-----------)  (---)'
      do k=1,kmax
        write(ifoutput,'(I3,F8.2,F7.1,5F10.4,F12.5,3F11.4,F11.5)') &
            k, &
            zf    (k), &
            presf (k)/100., &
            tmn   (k), &
            thlmn (k), &
            thmn  (k), &
            thvmn (k), &
            qtmn  (k)*1000., &
            qlmn  (k)*1000., &
            umn   (k), &
            vmn   (k), &
            cfracmn(k), &
            cszmn(k)
      end do
      close (ifoutput)

  !     -------------------------
  !     6.5   write the fluxes
  !     -----------------------------------------------------------------------

      open (ifoutput,file='flux1.'//cexpnr,position='append')

      write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
            '#-------------------------------------------------------------' &
          ,'---------------------------------)' &
            ,'#',(timeav),'--- AVERAGING TIMESTEP --- ' &
            ,nhrs,':',nminut,':',nsecs &
            ,'   HRS:MIN:SEC AFTER INITIALIZATION '
      write (ifoutput,'(2A/A/2A/2A/2A/2A)') &
            '#---------------------------------------------------------------' &
          ,'---------------------------------)' &
          ,'#                                                               ' &
          ,'#                  |                                  TURBULENT ' &
          ,'FLUXES                           |                             ' &
          ,'#LEV HEIGHT  PRES  |   WTL_SUB     WTL_RES     WTL_TOT         WQ' &
          ,'T_SUB      WQT_RES     WQT_TOT   ' &
          ,'#     (M)    (MB) | (----------   (K M/S)   --------------)     ' &
          ,'(---------- (M/S)   -------------)                             ' &
          ,'#---------------------------------------------------------------' &
          ,'---------------------------------)'

      write(ifoutput,'(I3,F8.2,F7.1,6E13.5)') &
            (k, &
            zh       (k), &
            presh     (k)/100., &
            wtlsmn   (k)                              , &
            wtlrmn   (k)                              , &
            wtltmn   (k)                              , &
            wqtsmn   (k)                              , &
            wqtrmn   (k)                              , &
            wqttmn   (k),&
             k=1,kmax)



      close(ifoutput)

      open (ifoutput,file='flux2.'//cexpnr,position='append')

      write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
            '#--------------------------------------------------------' &
            ,'#',(timeav),'--- AVERAGING TIMESTEP --- ' &
            ,nhrs,':',nminut,':',nsecs &
            ,'   HRS:MIN:SEC AFTER INITIALIZATION '
      write (ifoutput,'(A/A/A/3A/4A/2A)') &
            '#(------------------------------------------------------------)' &
          ,'#                                                             ' &
          ,'                                      TURBULENT FLUXES        ' &
          ,'#LEV HEIGHT  PRES  |    UW_TOT      VW_TOT       UW_SGS   ' &
          ,'   VW_SGS       UW_RES       VW_RES ' &
          ,'      WTH_TOT       WQ_L       WTHV_SUB     WTHV_RES     WTHV_TOT' &
          ,'#     (M)    (MB) |   ' &
          ,'(---------------------------- (M/S)^2  ---------------------------------)' &
          ,'   (-(K M/S)-)   (-(' &
          ,'M/S)-)     (----------  (K M/S)  -----------)' &
          ,'#(------------------------------------------------------------' &
          ,'-----------------------------------------)'

      write(ifoutput,'(I3,F8.2,F7.1,11E13.5)') &
          (k, &
            zh       (k), &
            presh     (k)/100., &
            uwtmn    (k)                              , &
            vwtmn    (k)                              , &
            uwsmn    (k)                              , &
            vwsmn    (k)                              , &
            uwrmn    (k)                              , &
            vwrmn    (k)                              , &
            wtltmn   (k) + wqltmn(k)*(rlv/cp)/exnh(k) , &
            wqltmn   (k)                              , &
            wtvsmn   (k)                              , &
            wtvrmn   (k)                              , &
            wtvtmn   (k)                              , &
            k=1,kmax)


      close(ifoutput)

      open (ifoutput,file='moments.'//cexpnr,position='append')

      write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
        '#--------------------------------------------------------'      &
        ,'#',(timeav),'--- AVERAGING TIMESTEP --- '      &
        ,nhrs,':',nminut,':',nsecs      &
        ,'   HRS:MIN:SEC AFTER INITIALIZATION '

      write (ifoutput,'(3A/3A/3A)') &
            '#----------------------------------------------------' &
            ,'---------------------------------------------------' &
            ,'------------------------------' &
            ,'#  LEV   HGHT   PRES     THL**2       THV**2         ' &
            ,'TH**2       QT**2       U*U    ' &
            ,'  V*V     HGHT     W*W     SKEWW     SFS-TKE' &
            ,'#        (M)   (MB)     (--------------(K*K)--------' &
            ,'------)     (G/KG)^2     (-(M/S)' &
            ,'^2)-)      (M)     (M/S)^2  ()        (M/S)^2'

      write(ifoutput,'(I5,F8.2,F7.1,4E13.5,2F9.4,F8.2,3F9.4)') &
            (k, &
            zf    (k), &
            presf (k)/100., &
            thl2mn(k), &
            thv2mn(k), &
            th2mn (k), &
            qt2mn (k)*1e6, &
!             qs2mn (k)*1e6, &
            u2mn  (k), &
            v2mn  (k), &
            zh    (k), &
            w2mn  (k), &
            skewmn(k), &
            w2submn(k), &
            k=1,kmax)

      close(ifoutput)


  !----   Write information about scalar field and its tendencies ------

      do n=1,nsv
        name = 'svnnnfld.'//cexpnr
        write (name(3:5),'(i3.3)') n
        open (ifoutput,file=name,position='append')

        write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
            '#-------------------------------------------------------------' &
            ,'---------------------------------)' &
            ,'#',(timeav),'--- AVERAGING TIMESTEP --- ' &
            ,nhrs,':',nminut,':',nsecs &
            ,'   HRS:MIN:SEC AFTER INITIALIZATION '

        write (ifoutput,'(2A/A/A/A,I1,A,/A/A/2A)') &
            '#--------------------------------------------------------- ' &
            ,'--------------' &
            ,'#               --FIELD &  T E N D E N C I E S  --------    ' &
            ,'#                                                           ' &
            ,'#                  |        SCALAR(' &
            ,n &
            ,')                  |' &
            ,'# LEV HEIGHT  PRES      SV(1)      TURB    TOTAL  |' &
            ,'#      (M)   (MB)  |       -----  (KG/KG/DAY) ----- |' &
            ,'#----------------------------------------------------------' &
            ,'-------------'
        write(ifoutput,'(I4,2F8.2,E14.5e3,2F8.2,E14.5E3)') &
            (k, &
              zf       (k), &
              presf    (k)/100., &
              svmmn    (k,n), &
              svptmn   (k,n)  *convt, &
              svpmn    (k,n)  *convt, &
              sv2mn    (k,n), &
              k=1,kmax)

        close(ifoutput)



  !        -----------------------
  !        Write the scalar fluxes
  !        -----------------------
        name = 'svnnnflx.'//cexpnr
        write (name(3:5),'(i3.3)') n
        open (ifoutput,file=name,position='append')

        write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
            '#-------------------------------------------------------------' &
            ,'---------------------------------)' &
            ,'#',(timeav),'--- AVERAGING TIMESTEP --- ' &
            ,nhrs,':',nminut,':',nsecs &
            ,'   HRS:MIN:SEC AFTER INITIALIZATION '
        write (ifoutput,'(2A/A/A,I1,A,/A/A/2A)') &
            '#---------------------------------------------------------------' &
            ,'---------------------------------)' &
            ,'#                                                               ' &
            ,'#                 |         TURBULENT FLUXES (SV=',n,')       |' &
            ,'#LEV HEIGHT  PRES       WSV_SUB     WSV_RES     WSV_TOT' &
            ,'#     (M)    (MB) | (----------   (KG/KG M/S)   --------------) ' &
            ,'#---------------------------------------------------------------' &
            ,'---------------------------------)'

        write(ifoutput,'(I3,2F8.2,3E14.5E3)') &
            (k, &
              zh       (k), &
              presh    (k)/100., &
              wsvsmn   (k,n)                            , &
              wsvrmn   (k,n)                            , &
              wsvtmn   (k,n)                            , &
              k=1,kmax)

         close(ifoutput)

      end do
      if (lnetcdf) then
        vars(:, 1)=rhof
        vars(:, 2)=presh
        vars(:, 3)=umn
        vars(:, 4)=vmn
        vars(:, 5)=thlmn
        vars(:, 6)=thvmn
        vars(:, 7)=qtmn
        vars(:, 8)=qlmn
        vars(:, 9)=wtlsmn
        vars(:,10)=wtlrmn
        vars(:,11)=wtltmn
        vars(:,12)=wtvsmn
        vars(:,13)=wtvrmn
        vars(:,14)=wtvtmn
        vars(:,15)=wqtsmn
        vars(:,16)=wqtrmn
        vars(:,17)=wqttmn
        vars(:,18)=wqlsmn
        vars(:,19)=wqlrmn
        vars(:,20)=wqltmn
        vars(:,21)=uwsmn
        vars(:,22)=uwrmn
        vars(:,23)=uwtmn
        vars(:,24)=vwsmn
        vars(:,25)=vwrmn
        vars(:,26)=vwtmn
        vars(:,27)=w2submn
        vars(:,28)=w2mn
        !vars(:,29)=w2submn+w2mn
        vars(:,29)=skewmn
        vars(:,30)=u2mn
        vars(:,31)=v2mn
        vars(:,32)=thl2mn
        vars(:,33)=thv2mn
        vars(:,34)=th2mn
        vars(:,35)=qt2mn
        vars(:,36)=ql2mn
        vars(:,37)=csz
        do n=1,nsv
          vars(:,37+7*(n-1)  )=svmmn(:,n)
          vars(:,37+7*(n-1)+1)=svpmn(:,n)
          vars(:,37+7*(n-1)+2)=svptmn(:,n)
          vars(:,37+7*(n-1)+3)=sv2mn(:,n)
          vars(:,37+7*(n-1)+4)=wsvsmn(:,n)
          vars(:,37+7*(n-1)+5)=wsvrmn(:,n)
          vars(:,37+7*(n-1)+6)=wsvtmn(:,n)
        end do
        call writestat_nc(ncid,1,tncname,(/rtimee/),nrec,.true.)
        call writestat_nc(ncid,nvar,ncname,vars(1:kmax,:),nrec,kmax)
      end if


    end if ! end if(myid==0)
      umn      = 0.
      vmn      = 0.
      thlmn    = 0.
      thvmn    = 0.
      qtmn     = 0.
      qlmn     = 0.
      qlhmn    = 0.
      cfracmn  = 0.

      wtlsmn =  0.
      wtlrmn =  0.
      wtltmn =  0.

      wtvsmn =  0.
      wtvrmn =  0.
      wtvtmn =  0.

      wqtsmn =  0.
      wqtrmn =  0.
      wqttmn =  0.

      wqlsmn =  0.
      wqlrmn =  0.
      wqltmn =  0.

      uwtmn  =  0.
      vwtmn  =  0.
      uwrmn  =  0.
      vwrmn  =  0.
      uwsmn  =  0.
      vwsmn  =  0.


      u2mn     = 0.
      v2mn     = 0.
      w2mn     = 0.
      w2submn  = 0.
      skewmn   = 0.
      qt2mn    = 0.
      thl2mn   = 0.
      thv2mn   = 0.
      th2mn    = 0.
      ql2mn    = 0.
!       qs2mn    = 0.
!       qsmn     = 0.
!       rhmn     = 0.
!       rmn      = 0.
!       r2mn     = 0.
!       r3mn     = 0.

      svmmn   = 0.
      svpmn   = 0.
!       svplsmn = 0.
      svptmn  = 0.

      sv2mn = 0.

      wsvsmn = 0.
      wsvrmn = 0.
      wsvtmn = 0.

      cszmn  = 0.

      deallocate(tmn, thmn)


  end subroutine writestat
  subroutine exitgenstat
    use modmpi, only : myid
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(.not.(lstat)) return

    if(lnetcdf .and. myid==0) call exitstat_nc(ncid)

    deallocate(umn       ,vmn   )
    deallocate(thlmn        ,thvmn )
    deallocate(qtmn         ,qlmn  ,  qlhmn, cfracmn)
    deallocate(wtlsmn ,wtlrmn ,wtltmn)
    deallocate(wtvsmn ,wtvrmn ,wtvtmn)
    deallocate(wqlsmn ,wqlrmn ,wqltmn)
    deallocate(wqtsmn ,wqtrmn ,wqttmn)
    deallocate(wsvsmn ,wsvrmn,wsvtmn)
    deallocate(uwtmn,vwtmn,uwrmn,vwrmn,uwsmn,vwsmn )
!     deallocate(rmn, r2mn   , r3mn , rhmn )
    deallocate(w2mn       , skewmn )
    deallocate(w2submn    )
    deallocate(u2mn       , v2mn  ,     qt2mn)
    deallocate(thl2mn     , thv2mn,     th2mn,     ql2mn)
!     deallocate(qs2mn      , qsmn  )
    deallocate(svmmn,svptmn,svplsmn,svpmn)
    deallocate(sv2mn)

    deallocate(umav )
    deallocate(vmav )
    deallocate(thlmav )
    deallocate(thmav )
    deallocate(qtmav )
    deallocate(qlmav )
    deallocate(cfracav )
    deallocate(svmav )
    deallocate(uptav)
    deallocate(vptav)

    deallocate(thptav)
    deallocate(qlptav)
    deallocate(uwtot )
    deallocate(vwtot )
    deallocate(uwres )
    deallocate(vwres )
    deallocate(uwsub )
    deallocate(vwsub )
    deallocate(wtlsub)
    deallocate(wtlres)
    deallocate(wtltot)
    deallocate(wqtsub)
    deallocate(wqtres)
    deallocate(wqttot)
    deallocate(wqlsub)
    deallocate(wqlres)
    deallocate(wqltot)
    deallocate(wtvtot)
    deallocate(wtvres)
    deallocate(wtvsub)
    deallocate(wsvsub)
    deallocate(wsvres)
    deallocate(wsvtot)
    deallocate(thvhav)
    deallocate(th0av)
    deallocate(svpav)
    deallocate(svptav)
    
    deallocate(cszmn)
    deallocate(cszav)

  end subroutine exitgenstat


end module modgenstat
