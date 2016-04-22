!> \file modstartup.f90
!!  Initializes the run

!>
!! Initializes the run.
!>
!! Modstartup reads the namelists and initial data, sets the fields and calls
!! the inits of the other routines where necessary. Reading of the
!! restart files also live in this module.
!!  \author Jasper Tomas, TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \todo documentation
!!  \par Revision list
!

module modstartup


implicit none
! private
! public :: startup,trestart
save

  integer (KIND=selected_int_kind(6)) :: irandom= 0     !    * number to seed the randomnizer with
  integer :: krand = huge(0)
  real :: randu = 2.0,randthl= 0.1,randqt=1e-5         !    * uvw,thl and qt amplitude of randomnization

contains
  subroutine startup

      !-----------------------------------------------------------------|
      !                                                                 |
      !     Reads all general options from namoptions                   |
      !                                                                 |
      !      Jasper Tomas                 31/03/2014                    |
      !      Chiel van Heerwaarden        15/06/2007                    |
      !      Thijs Heus                   15/06/2007                    |
      !-----------------------------------------------------------------|

    use modglobal,         only : initglobal,iexpnr,runtime, dtmax,dtav_glob,timeav_glob,&
                                  lwarmstart,lreadscal,startfile,tfielddump,trestart,trestart2,&
                                  nsv,nvar,imax,jtot,kmax,xsize,ysize,xlat,xlon,xday,xtime,lwalldist,&
                                  lmoist,lcoriol,igrw_damp,geodamptime,lmomsubs,cu, cv,ifnamopt,fname_options,llsadv,&
                                  xS,yS,zS,SS,sigS,iadv_mom,iadv_tke,iadv_thl,iadv_qt,iadv_sv,courant,peclet,ladaptive,author,&
                                  linoutflow,lper2inout,libm,lles,lwallfunc,lmassflowr,massflowrate,lstoreplane,iplane,&
                                  lreadmean,linletgen,inletav,lreadminl,Uinf,linletRA,nxwalls,nywalls,nzwalls,nblocks,&
                                  ltec3d,lscalinout,lscalrec,lscasrc,lzerogradtop,lzerogradtopscal,lbuoyancy,lneutral,ltempeq,numol,prandtlmol,numoli,prandtlmoli,&
                                  lfixinlet,lfixutauin,lstore3d,lstorexz,startmean,pi,lMOST,ltfluxtop,ltfluxbot,ltec2d,&
                                  thlsrc,nkplane,kplane,nsvl,nsvp,ifixuinf,tscale,ltempinout,lmoistinout,lnetcdf
    use modsurfdata,       only : z0,ustin,wtsurf,wttop,wqsurf,wsvsurf,ps,thvs,thls,thl_top,qt_top,qts,isurf
    use modsurface,        only : initsurface
    use modfields,         only : initfields,dpdx,ncname
    use modpois,           only : initpois
    use modboundary,       only : initboundary,ksp
    use modthermodynamics, only : initthermodynamics,lqlnr, chi_half
    use modsubgrid,        only : initsubgrid
    use modsubgriddata,    only : ldelta, cf,cn,Rigc,Prandtl,lmason
    use modmpi,            only : comm3d,myid, mpi_integer,mpi_logical,my_real,mpierr, mpi_character
    use modforces,          only : initforces
    use modinlet,          only : initinlet
    use modinletdata,      only : di,dr,di_test,dti,iangledeg,iangle

    implicit none
    integer :: ierr

  !declare namelists

    namelist/RUN/ &
        iexpnr,lwarmstart,lreadscal,startfile, runtime, dtmax,dtav_glob,timeav_glob,&
        trestart,trestart2,tfielddump,irandom,randthl,randqt,krand,nsv,courant,peclet,ladaptive,&
        author,linoutflow,lper2inout,libm,lles,lwallfunc,lMOST,lmassflowr,lreadmean,& 
        startmean,ltec3d,lscalinout,lscalrec,lscasrc,ltempinout,lmoistinout,lzerogradtop,lzerogradtopscal,lwalldist,lstore3d,lstorexz,&
        ltfluxtop,ltfluxbot,randu,ltec2d,nkplane,kplane,nsvl,nsvp,ifixuinf,tscale,dpdx,lnetcdf
    namelist/DOMAIN/ &
        imax,jtot,kmax,&
        xsize,ysize,&
        xlat,xlon,xday,xtime,ksp,&
        nxwalls,nywalls,nzwalls,nblocks
    namelist/INLET/ &
        linletgen,iangledeg,Uinf,di,dti,iplane,lstoreplane,inletav,lreadminl,linletRA,lfixinlet,lfixutauin,xS,yS,zS,SS,sigS
    namelist/PHYSICS/ &
         z0,ustin,ps,wtsurf,wttop,wqsurf,wsvsurf,thvs,thls,thl_top,qt_top,qts,lmoist,isurf,&
        lcoriol,igrw_damp,geodamptime,massflowrate,numol,prandtlmol,&
        lbuoyancy,ltempeq,thlsrc,lneutral
    namelist/DYNAMICS/ &
        llsadv, lqlnr, cu, cv, iadv_mom, iadv_tke, iadv_thl, iadv_qt, iadv_sv
    

    !add nsv to nvar

    nvar=nvar+nsv
    allocate(ncname(nvar,4))
    write(*,*) 'nvar', nvar
    !read namelists
    

     

    if(myid==0)then
      if (command_argument_count() >=1) then
        call get_command_argument(1,fname_options)
      end if
      write (*,*) fname_options

      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      if (ierr /= 0) then
        stop 'ERROR:Namoptions does not exist'
      end if
      read (ifnamopt,RUN,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions RUN'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions RUN'
      endif
      write(6 ,RUN)
      rewind(ifnamopt)
      read (ifnamopt,DOMAIN,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions DOMAIN'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions DOMAIN'
      endif
      write(6 ,DOMAIN)
      rewind(ifnamopt)
      read (ifnamopt,INLET,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions INLET'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions INLET'
      endif
      write(6 ,INLET)
      rewind(ifnamopt)
      read (ifnamopt,PHYSICS,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions PHYSICS'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions PHYSICS'
      endif
      write(6 ,PHYSICS)
      rewind(ifnamopt)
      read (ifnamopt,DYNAMICS,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions DYNAMICS'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions DYNAMICS'
      endif
      write(6 ,DYNAMICS)
      close(ifnamopt)
    end if
  
  
  !broadcast namelists
    call MPI_BCAST(iexpnr     ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(lwarmstart ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lreadscal  ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch to read scalar pollutant fields (warm start)
    call MPI_BCAST(linoutflow ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for inflow/outflow in i-direction
    call MPI_BCAST(lscalinout ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for inflow/outflow in i-direction for scalar
    call MPI_BCAST(lscalrec   ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lscasrc    ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(ltempinout ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for inflow/outflow in i-direction for temperature
    call MPI_BCAST(lmoistinout,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for inflow/outflow in i-direction for moisture
    call MPI_BCAST(lzerogradtop ,1,MPI_LOGICAL,0,comm3d,mpierr)      ! J.Tomas: added switch for zero-gradient BC at top boundary
    call MPI_BCAST(lzerogradtopscal ,1,MPI_LOGICAL,0,comm3d,mpierr)      ! ILS13
    call MPI_BCAST(ltfluxtop ,1,MPI_LOGICAL,0,comm3d,mpierr)      ! ILS13
    call MPI_BCAST(ltfluxbot ,1,MPI_LOGICAL,0,comm3d,mpierr)      ! ILS13
    call MPI_BCAST(lbuoyancy  ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for buoyancy force in modforces
    call MPI_BCAST(lneutral  ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! ILS13: added switch for neutral atmosphere (no stability correction)   
    call MPI_BCAST(ltempeq    ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for solving adv/diff equation for temperature
    call MPI_BCAST(lper2inout ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for restart periodic flow to inoutflow
    call MPI_BCAST(libm       ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for turning on/off IBM method
    call MPI_BCAST(lwalldist  ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for computing wall distances
    call MPI_BCAST(lles       ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for turning on/off LES functionality (subgrid model)
    call MPI_BCAST(ltec3d     ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for turning on/off writing formatted tec3d file
    call MPI_BCAST(ltec2d     ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for turning on/off writing formatted tec3d file
    call MPI_BCAST(linletgen  ,1,MPI_INTEGER,0,comm3d,mpierr)        ! J.Tomas: 0: no inletgen, 1: Lund (1998), 2: read inlet from file
    call MPI_BCAST(linletRA   ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for turning on/off Running Average in inletgenerator
    call MPI_BCAST(lfixinlet  ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for keeping average inlet velocit and temp fixed at inlet (linletgen=1,2)
    call MPI_BCAST(lfixutauin ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for keeping utau fixed at inlet (linletgen=1,2)
    call MPI_BCAST(xS ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(yS ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(zS ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(SS ,1,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(sigS ,1,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(lreadminl,  1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for reading mean inlet/recycle plane profiles (Uinl,Urec,Wrec)
    call MPI_BCAST(lwallfunc  ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for turning on/off wall function
    call MPI_BCAST(lMOST  ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for turning on/off wall function
    call MPI_BCAST(lmassflowr ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for turning on/off u-velocity correction for fixed mass flow rate
    call MPI_BCAST(lstoreplane,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for turning on/off for storing i-plane data to serve as inlet for future sim. 
    call MPI_BCAST(lstore3d   ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for turning on/off storing 3d fields
    call MPI_BCAST(lstorexz   ,1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for turning on/off storing xz fields
    call MPI_BCAST(lreadmean,  1,MPI_LOGICAL,0,comm3d,mpierr)        ! J.Tomas: added switch for reading mean variables from means#MYID#.#EXPNR#
    call MPI_BCAST(iplane     ,1,MPI_INTEGER,0,comm3d,mpierr)        ! J.Tomas: ib+iplane is the i-plane that is stored if lstoreplane is .true.
    call MPI_BCAST(startfile  ,50,MPI_CHARACTER,0,comm3d,mpierr)
    call MPI_BCAST(author     ,80,MPI_CHARACTER,0,comm3d,mpierr)
    call MPI_BCAST(runtime    ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(trestart   ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(trestart2  ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(tfielddump ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dtmax      ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dtav_glob  ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(timeav_glob,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(nsv        ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(nvar       ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(ncname     ,80,MPI_CHARACTER,0,comm3d,mpierr) 
    call MPI_BCAST(lnetcdf    ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(ifixuinf   ,1,MPI_INTEGER,0,comm3d,mpierr) 
    call MPI_BCAST(dpdx       ,1,MY_REAL,0,comm3d,mpierr) 
    call MPI_BCAST(nsvp        ,1,MPI_INTEGER,0,comm3d,mpierr) 
    call MPI_BCAST(nsvl        ,1,MPI_INTEGER,0,comm3d,mpierr) 
    call MPI_BCAST(tscale     ,1,MY_REAL,0,comm3d,mpierr)     
    !nsv=nsvl+nsvp
    call MPI_BCAST(imax       ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(jtot       ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(kmax       ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(xsize      ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ysize      ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(xlat       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(xlon       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(xday       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(xtime      ,1,MY_REAL   ,0,comm3d,mpierr)

    call MPI_BCAST(z0         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ustin      ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wtsurf     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wqsurf     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wsvsurf(1:nsv),nsv,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ps         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(thvs       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(thls       ,1,MY_REAL   ,0,comm3d,mpierr)      
    call MPI_BCAST(thl_top    ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(qt_top     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(qts        ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(lmoist     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lcoriol    ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(igrw_damp  ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(geodamptime,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wttop      ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(thlsrc     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(startmean  ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(massflowrate,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(Uinf       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(numol      ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(prandtlmol ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(di         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dti        ,1,MY_REAL   ,0,comm3d,mpierr)
    dr      = di                  ! initial value is needed
    di_test = di                  ! initial value is needed
    numoli = 1./numol
    prandtlmoli = 1./prandtlmol

    call MPI_BCAST(numoli     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(prandtlmoli,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(iangledeg,1,MY_REAL   ,0,comm3d,mpierr)
    iangle = iangledeg * pi / 180.  
     call MPI_BCAST(inletav     ,1,MY_REAL   ,0,comm3d,mpierr)

    call MPI_BCAST(llsadv     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lqlnr      ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(cu         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(cv         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ksp        ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(nxwalls    ,1,MPI_INTEGER,0,comm3d,mpierr)    ! no. of x-walls used in IBM
    call MPI_BCAST(nywalls    ,1,MPI_INTEGER,0,comm3d,mpierr)    ! no. of y-walls used in IBM
    call MPI_BCAST(nzwalls    ,1,MPI_INTEGER,0,comm3d,mpierr)    ! no. of z-walls used in IBM
    call MPI_BCAST(nblocks    ,1,MPI_INTEGER,0,comm3d,mpierr)    ! no. of blocks used in IBM
    call MPI_BCAST(irandom    ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(krand      ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(nkplane    ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(kplane    ,nkplane,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(randthl    ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(randu    ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(randqt     ,1,MY_REAL   ,0,comm3d,mpierr)

    call MPI_BCAST(ladaptive  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(courant,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(peclet,1,MY_REAL   ,0,comm3d,mpierr)

    call MPI_BCAST(isurf   ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(iadv_mom,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(iadv_tke,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(iadv_thl,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(iadv_qt ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(iadv_sv(1:nsv) ,nsv,MPI_INTEGER,0,comm3d,mpierr)

    ! Allocate and initialize core modules
    call initglobal
    call initfields

    call initboundary
    call initthermodynamics
    call initsurface
    call initsubgrid
    call initpois
    call initinlet   ! added by J. Tomas: initialize inlet generator

    call checkinitvalues
    call readinitfiles
    call initforces

  end subroutine startup


  subroutine checkinitvalues
  !-----------------------------------------------------------------|
  !                                                                 |
  !      Thijs Heus   TU Delft  9/2/2006                            |
  !                                                                 |
  !     purpose.                                                    |
  !     --------                                                    |
  !                                                                 |
  !      checks whether crucial parameters are set correctly        |
  !                                                                 |
  !     interface.                                                  |
  !     ----------                                                  |
  !                                                                 |
  !     *checkinitvalues* is called from *program*.                 |
  !                                                                 |
  !-----------------------------------------------------------------|

    use modsurfdata,only : wtsurf,wqsurf,ustin,thvs,qts,isurf,ps
    use modglobal, only : imax,jtot, ysize,xsize,dtmax,runtime, startfile,lwarmstart
    use modmpi,    only : myid, nprocs,mpierr


      if(mod(jtot,nprocs) /= 0) then
        if(myid==0)then
          write(6,*)'STOP ERROR IN NUMBER OF PROCESSORS'
          write(6,*)'nprocs must divide jtot!!! '
          write(6,*)'nprocs and jtot are: ',nprocs, jtot
        end if
        call MPI_FINALIZE(mpierr)
        stop
      end if


      !ILS13, 18.04.16  is this necessary? no parallelisation in x....

      !if(mod(imax,nprocs)/=0)then
      !  if(myid==0)then
      !    write(6,*)'STOP ERROR IN NUMBER OF PROCESSORS'
      !    write(6,*)'nprocs must divide imax!!! '
      !    write(6,*)'nprocs and imax are: ',nprocs,imax
      !  end if
      !  call MPI_FINALIZE(mpierr)
      !  stop
      !end if

  !Check Namroptions


    if (runtime < 0)stop 'runtime out of range/not set'
    if (dtmax < 0)  stop 'dtmax out of range/not set '
    if (ps < 0)     stop 'psout of range/not set'
    !if (thvs < 0)   stop 'thvs out of range/not set'
    if (xsize < 0)  stop 'xsize out of range/not set'
    if (ysize < 0)  stop 'ysize out of range/not set '

    if (lwarmstart) then
      if (startfile == '') stop 'no restartfile set'
    end if
  !isurf
    if (myid == 0) then
      select case (isurf)
      case(1)
      case(2,10)
      case(3:4)
        if (wtsurf == -1)  stop 'wtsurf not set'
        if (wqsurf == -1)  stop 'wqsurf not set'
      case default
        stop 'isurf out of range/not set'
      end select
      if (isurf ==3) then
        if (ustin < 0)  stop 'ustin out of range/not set'
      end if
    end if

  end subroutine checkinitvalues

  subroutine readinitfiles
    use modfields,         only : u0,v0,w0,um,vm,wm,thlm,thl0,thl0h,qtm,qt0,qt0h,&
                                  ql0,ql0h,thv0h,sv0,svm,e12m,e120,&
                                  dudxls,dudyls,dvdxls,dvdyls,dthldxls,dthldyls,&
                                  dqtdxls,dqtdyls,dqtdtls,dpdx,dpdxl,dpdyl,&
                                  wfls,whls,ug,vg,pgx,pgy,uprof,vprof,thlprof, qtprof,e12prof, svprof,&
                                  v0av,u0av,qt0av,ql0av,thl0av,qt0av,sv0av,exnf,exnh,presf,presh,rhof,&
                                  thlpcar,uav,thvh,thvf
    use modglobal,         only : ib,ie,ih,jb,je,jh,kb,ke,kh,kmax,dtmax,dt,runtime,timeleft,timee,ntimee,ntrun,btime,dt_lim,nsv,&
                                  zf,zh,dzf,dzh,rv,rd,grav,cp,rlv,pref0,om23_gs,jgb,jge,Uinf,dy,&
                                  rslabs,cu,cv,e12min,dzh,dtheta,dqt,dsv,cexpnr,ifinput,lwarmstart,trestart,numol,&
                                  ladaptive,llsadv,tnextrestart,jmax,linoutflow,lper2inout,linletgen,lreadminl,&
                                  massflowrate,trestart2,tnextrestart2,ltempeq,prandtlmoli,freestreamav,tnextfielddump,tfielddump
    use modsubgriddata,    only : ekm,ekh
    use modsurfdata,       only : wtsurf,wqsurf,wsvsurf, &
                                  thls,thvs,ustin,ps,qts,isurf,svs,obl,oblav
    use modsurface,        only : surface,dthldz
    use modboundary,       only : boundary,bottom,tqaver
    use modmpi,            only : slabsum,myid,comm3d,mpierr,my_real
    use modthermodynamics, only : thermodynamics,calc_halflev
    use modinletdata,      only : Uinl,Urec,Wrec,u0inletbc,v0inletbc,w0inletbc,ubulk,irecy,Utav,Ttav,&
                                  uminletbc,vminletbc,wminletbc,u0inletbcold,v0inletbcold,w0inletbcold,&
                                  storeu0inletbc,storev0inletbc,storew0inletbc,nstepread,nfile,Tinl,&
                                  Trec,tminletbc,t0inletbcold,t0inletbc,storet0inletbc,utaui,ttaui,iangle
    use modinlet,          only : readinletfile

    integer i,j,k,n

    real, allocatable :: height(:), th0av(:)
    real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh) :: thv0
    real, dimension(kb:ke) :: uaverage        ! volume averaged u-velocity
    real, dimension(kb:ke) :: uaverager       ! recycle plane
    real, dimension(kb:ke) :: uaveragei       ! inlet plane
    real, dimension(kb:ke) :: taverager       ! recycle plane
    real, dimension(kb:ke) :: taveragei       ! inlet plane
    real, dimension(kb:ke+1) :: waverage
    real, dimension(kb:ke+1) :: uprofrot
    real, dimension(kb:ke+1) :: vprofrot
    real tv,ran,ran1


    character(80) chmess

    allocate (height(kb:ke+kh))
    allocate (th0av(kb:ke+kh))



    if (.not. lwarmstart) then

    !********************************************************************

    !    1.0 prepare initial fields from files 'prof.inp' and 'scalar.inp'
    !    ----------------------------------------------------------------

    !--------------------------------------------------------------------
    !    1.1 read fields
    !-----------------------------------------------------------------

      dt = dtmax / 100.
      timee = 0
      if (myid==0) then
        open (ifinput,file='prof.inp.'//cexpnr)
        read (ifinput,'(a80)') chmess
        write(*,     '(a80)') chmess
        read (ifinput,'(a80)') chmess

        do k=kb,ke
          read (ifinput,*) &
                height (k), &
                thlprof(k), &
                qtprof (k), &
                uprof  (k), &
                vprof  (k), &
                e12prof(k)
        end do

        ! Apply rotation in horizontal
        write(6,*) 'iangle = ',iangle
        
        uprofrot = uprof * cos(iangle) - vprof * sin(iangle)
        vprofrot = vprof * cos(iangle) + uprof * sin(iangle)
        uprof    = uprofrot
        vprof    = vprofrot

        close(ifinput)
        write(*,*) 'height    thl     qt      u      v     e12'
        do k=ke,kb,-1
          write (*,'(f7.1,2f8.1,3f7.1)') &
                height (k), &
                thlprof(k), &
                qtprof (k), &
                uprof  (k), &
                vprof  (k), &
                e12prof(k)

        end do

        if (minval(e12prof(kb:ke)) < e12min) then
          write(*,*)  'e12 value is zero (or less) in prof.inp'
          do k=kb,ke
            e12prof(k) = max(e12prof(k),e12min)
          end do
        end if

      end if ! end if myid==0
    ! MPI broadcast numbers reading
      call MPI_BCAST(thlprof,kmax,MY_REAL   ,0,comm3d,mpierr)
      call MPI_BCAST(qtprof ,kmax,MY_REAL    ,0,comm3d,mpierr)
      call MPI_BCAST(uprof  ,kmax,MY_REAL   ,0,comm3d,mpierr)
      call MPI_BCAST(vprof  ,kmax,MY_REAL   ,0,comm3d,mpierr)
      call MPI_BCAST(e12prof,kmax,MY_REAL   ,0,comm3d,mpierr)
      do k=kb,ke
      do j=jb-1,je+1
      do i=ib-1,ie+1
        thl0(i,j,k) = thlprof(k)
        thlm(i,j,k) = thlprof(k)
        qt0 (i,j,k) = qtprof (k)
        qtm (i,j,k) = qtprof (k)
        u0  (i,j,k) = uprof  (k) - cu
        um  (i,j,k) = uprof  (k) - cu
        v0  (i,j,k) = vprof  (k) - cv
        vm  (i,j,k) = vprof  (k) - cv
        w0  (i,j,k) = 0.0
        wm  (i,j,k) = 0.0
        e120(i,j,k) = e12prof(k)
        e12m(i,j,k) = e12prof(k)
!        ekm (i,j,k) = 0.0
!        ekh (i,j,k) = 0.0
        ekm (i,j,k) = numol
        ekh (i,j,k) = numol
      end do
      end do
      end do

!! add random fluctuations
       krand  = min(krand,ke)
        do k = kb,krand
          call randomnize(um,k,randu,irandom,ih,jh)
        end do
        do k = kb,krand
          call randomnize(vm,k,randu,irandom,ih,jh)
        end do
        do k = kb,krand
          call randomnize(wm,k,randu,irandom,ih,jh)
        end do




!       do k=kb+1,ke-1
!       do j=jb,je
!       do i=ib+1,ie-1
!         call random_number(ran)
!         ran1 = -1. +2.*ran
!         wm(i,j,k)=wm(i,j,k)+ 0.1*Uinf*ran1
!       end do
!       end do
!       end do
!
!       do k=kb+1,ke-1
!       do j=jb,je
!       do i=ib+2,ie-1
!         call random_number(ran)
!         ran1 = -1. +2.*ran
!         um(i,j,k)=um(i,j,k)+ 0.1*Uinf*ran1
!       end do
!       end do
!       end do
!
!       do k=kb+1,ke-1
!       do j=jb,je
!       do i=ib+1,ie-1
!         call random_number(ran)
!         ran1 = -1. +2.*ran
!         vm(i,j,k)=vm(i,j,k)+ 0.1*Uinf*ran1
!       end do
!       end do
!       end do

       u0=um
       v0=vm
       w0=wm


        uaverage=0.
        call slabsum(uaverage ,kb,ke,um,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,ib,ie,jb,je,kb,ke)       
        do k=kb,ke
          uaverage(k) = uprof(k)*dzf(k)
        end do
        ubulk = sum(uaverage(kb:ke))/(zh(ke+1)-zh(kb))      ! volume-averaged u-velocity
        massflowrate = ubulk*(jge-jgb+1)*dy*(zh(ke+1)-zh(kb))*1000.
        write(6,*) 'Modstartup: massflowrate=',massflowrate
        

! Set average inlet profile to initial inlet profile in case of inletgenerator mode
      if (linletgen==1) then
        uaverage=0.
        call slabsum(uaverage ,kb,ke,um,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,ib,ie,jb,je,kb,ke)       
        do k=kb,ke
          uaverage(k) = uprof(k)*dzf(k)
        end do
        ubulk = sum(uaverage(kb:ke))/(zh(ke+1)-zh(kb))      ! volume-averaged u-velocity
        write(6,*) 'Modstartup: ubulk=',ubulk
        Utav(ib:ie,kb:ke) = um(ib:ie,jb,kb:ke) 
        Uinl = um(ib,jb,kb:ke)               ! set the initial time-averaged inlet profile equal to um
        Urec = um(ib,jb,kb:ke)               ! set the initial time-averaged inlet profile equal to um
        Wrec(kb:ke+1) = wm(ib,jb,kb:ke+1)    ! set the initial time-averaged inlet profile equal to mean w-profile
        u0inletbcold(jb:je,kb:ke) = um(ib,jb:je,kb:ke)
        v0inletbcold(jb:je,kb:ke) = vm(ib-1,jb:je,kb:ke)
        w0inletbcold(jb:je,kb:ke+1) = wm(ib-1,jb:je,kb:ke+1)
        uminletbc(jb:je,kb:ke) = um(ib,jb:je,kb:ke)
        vminletbc(jb:je,kb:ke) = vm(ib-1,jb:je,kb:ke)
        wminletbc(jb:je,kb:ke) = wm(ib-1,jb:je,kb:ke)
        u0inletbc(jb:je,kb:ke) = um(ib,jb:je,kb:ke)
        v0inletbc(jb:je,kb:ke) = vm(ib-1,jb:je,kb:ke)
        w0inletbc(jb:je,kb:ke+1) = wm(ib-1,jb:je,kb:ke+1)
        utaui = sqrt(abs(2*numol* Uinl(kb)/dzf(kb)))     ! average streamwise friction at inlet (need for first time step)

        if (ltempeq == .true.) then  
          Ttav(ib:ie,kb:ke) = thlm(ib:ie,jb,kb:ke)             ! set the initial time-averaged inlet profile equal to thlm
          Tinl = thlm(ib,jb,kb:ke)             ! set the initial time-averaged inlet profile equal to thlm
          Trec = thlm(ib,jb,kb:ke)             ! set the initial time-averaged inlet profile equal to thlm
          t0inletbcold(jb:je,kb:ke) = thlm(ib-1,jb:je,kb:ke)
          t0inletbc(jb:je,kb:ke) = thl0(ib-1,jb:je,kb:ke)
          tminletbc(jb:je,kb:ke) = thlm(ib-1,jb:je,kb:ke)
          ttaui = numol*prandtlmoli*2.*(Tinl(kb)-thls)/(dzf(kb)*utaui) ! average friction temp. at inlet (need for first time step)
        end if

       ! add random perturbations
       if (myid==0) then
       call random_number(ran)
       ran1 = -1. +2.*ran
       write(6,*) 'random=', ran,ran1        
       call random_number(ran)
       ran1 = -1. +2.*ran
       write(6,*) 'random=', ran,ran1        
       call random_number(ran)
       ran1 = -1. +2.*ran
       write(6,*) 'random=', ran,ran1        
       call random_number(ran)
       ran1 = -1. +2.*ran
       write(6,*) 'random=', ran,ran1        
       call random_number(ran)
       ran1 = -1. +2.*ran
       write(6,*) 'random=', ran,ran1        
       call random_number(ran)
       ran1 = -1. +2.*ran
       write(6,*) 'random=', ran,ran1        
       end if

       do k=kb+1,kb+48
       do j=jb,je
       do i=ib+1,ie-1
         call random_number(ran)
         ran1 = -1. +2.*ran
         wm(i,j,k)=wm(i,j,k)+ 0.1*Uinf*ran1
       end do
       end do
       end do

!       do k=kb+1,ke-1
       do k=kb+1,kb+48
       do j=jb,je
       do i=ib+2,ie-1
         call random_number(ran)
         ran1 = -1. +2.*ran
         um(i,j,k)=um(i,j,k)+ 0.1*Uinf*ran1
       end do
       end do
       end do

!       do k=kb+1,ke-1
       do k=kb+1,kb+48
       do j=jb,je
       do i=ib+1,ie-1
         call random_number(ran)
         ran1 = -1. +2.*ran
         vm(i,j,k)=vm(i,j,k)+ 0.1*Uinf*ran1
       end do
       end do
       end do

       u0=um
       v0=vm
       w0=wm

      else if (linletgen == 2) then 

        nfile = nfile + 1
        call readinletfile
        u0inletbc(:,:) = storeu0inletbc(:,:,nstepread)
        v0inletbc(:,:) = storev0inletbc(:,:,nstepread)
        w0inletbc(:,:) = storew0inletbc(:,:,nstepread)
        uminletbc(:,:) = storeu0inletbc(:,:,nstepread)
        vminletbc(:,:) = storev0inletbc(:,:,nstepread)
        wminletbc(:,:) = storew0inletbc(:,:,nstepread)
        ! determine bulk velocity
        call slabsum(uaverage ,kb,ke,u0,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,ib,ie,jb,je,kb,ke)
        uaverage = uaverage / ((ie-ib+1)*(jge-jgb+1))  ! this gives the i-j-averaged velocity (only correct for equidistant grid?)        
        do k=kb,ke
          uaverage(k) = uaverage(k)*dzf(k)
        end do
        ubulk = sum(uaverage(kb:ke))/(zh(ke+1)-zh(kb))      ! volume-averaged u-velocity
        write(6,*) 'Modstartup: ubulk=',ubulk
      end if ! linletgen
 

    !---------------------------------------------------------------
    !  1.2 randomnize fields
    !---------------------------------------------------------------
!     if (linletgen /= 2 .and. linletgen /= 1) then
!       write(6,*) 'randomnizing temperature!' 
!       krand  = min(krand,ke)
!        do k = kb,krand
!          call randomnize(thlm,k,randthl,irandom,ih,jh)
!          call randomnize(thl0,k,randthl,irandom,ih,jh)
!        end do
!      end if

      svprof = 0.
      if(myid==0)then
        if (nsv>0) then
          open (ifinput,file='scalar.inp.'//cexpnr)
          read (ifinput,'(a80)') chmess
          read (ifinput,'(a80)') chmess
          do k=kb,ke
            read (ifinput,*) &
                  height (k), &
                  (svprof (k,n),n=1,nsv)
          end do
          open (ifinput,file='scalar.inp.'//cexpnr)
          write (6,*) 'height   sv(1) --------- sv(nsv) '
          do k=ke,kb,-1
            write (6,*) &
                  height (k), &
                (svprof (k,n),n=1,nsv)
          end do

        end if
      end if ! end if myid==0

      call MPI_BCAST(wsvsurf,nsv   ,MY_REAL   ,0,comm3d,mpierr)
     
      call MPI_BCAST(svprof ,(ke+kh-(kb-kh))*nsv,MY_REAL   ,0,comm3d,mpierr)
      do k=kb,ke
        do j=jb-1,je+1
          do i=ib-1,ie+1
            do n=1,nsv
              sv0(i,j,k,n) = svprof(k,n)
              svm(i,j,k,n) = svprof(k,n)
            end do
          end do
        end do
      end do

!-----------------------------------------------------------------
!    2.2 Initialize surface layer
!-----------------------------------------------------------------
      obl=-0.1
      oblav=-0.1
      thvs = thls * (1. + (rv/rd -1.) * qts) 
      call boundary
      call thermodynamics         ! turned of when pot. temp = temp.
      call surface


      call boundary
      call thermodynamics        ! turned of when pot. temp = temp.

    else !if lwarmstart

      call readrestartfiles
      um   = u0
      vm   = v0
      wm   = w0
      thlm = thl0
      qtm  = qt0
      svm  = sv0
      e12m = e120

!ILS13 reintroduced thv
      call calc_halflev
     ! exnf = (presf/pref0)**(rd/cp)  !exner functions not in restart files
     ! anymore.. or at least not read
     ! exnh = (presh/pref0)**(rd/cp)

     !   write(*,*) "exnf",enf
     !   write(*,*) "exnh",exnh

      do  j=jb,je
      do  i=ib,ie
      do  k=kb,ke+kh
        !write(*,*) "thl0h",thl0h(i,j,k)
        thv0h(i,j,k) = (thl0h(i,j,k)+rlv*ql0h(i,j,k)/(cp)) &
                      *(1+(rv/rd-1)*qt0h(i,j,k)-rv/rd*ql0h(i,j,k))
      end do
      end do
      end do

      do  j=j,je
      do  i=ib,ie
      do  k=kb,ke+kh
        thv0(i,j,k) = (thl0(i,j,k)+rlv*ql0(i,j,k)/(cp)) &
                      *(1+(rv/rd-1)*qt0(i,j,k)-rv/rd*ql0(i,j,k))
      end do
      end do
      end do

      thvh=0.
      call slabsum(thvh,kb,ke,thv0h,ib-ih,ie+ih,jb-jh,je+jh,kb-kh,ke+kh,ib,ie,jb,je,kb,ke) !redefine halflevel thv using calculated thv
      thvh = thvh/rslabs

      thvf = 0.0
      call slabsum(thvf,kb,ke,thv0,ib-ih,ie+ih,jb-jh,je+jh,kb-kh,ke+kh,ib,ie,jb,je,kb,ke)
      thvf = thvf/rslabs





! Set average inlet profile to initial inlet profile in case of inletgenerator mode
      uaverage=0.
      uaveragei=0.
      uaverager=0.
      waverage=0.
      taveragei=0.
      taverager=0.
      if (linletgen==1) then
        call slabsum(uaveragei ,kb,ke,u0,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,ib,ib,jb,je,kb,ke)
        call slabsum(uaverager ,kb,ke,u0,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,irecy,irecy,jb,je,kb,ke)
        call slabsum(waverage ,kb,ke+1,w0,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,ib,ie,jb,je,kb,ke+1)
        call slabsum(uaverage ,kb,ke,u0,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,ib,ie,jb,je,kb,ke)
        uaverage = uaverage / ((ie-ib+1)*(jge-jgb+1))  ! this gives the i-j-averaged velocity (only correct for equidistant grid?)            
        uaveragei = uaveragei / (jge-jgb+1)            ! this gives the j-averaged u-velocity at the inlet
        uaverager = uaverager / (jge-jgb+1)            ! this gives the j-averaged u-velocity at the recycle plane
        waverage = waverage / ((ie-ib+1)*(jge-jgb+1))  ! this gives the i-j-averaged w-velocity (only correct for equidistant grid?)
        if (ltempeq == .true.) then
          call slabsum(taveragei ,kb,ke,thl0,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,ib,ie,jb,je,kb,ke)
          call slabsum(taverager ,kb,ke,thl0,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,irecy-1,irecy-1,jb,je,kb,ke)
          taveragei = taveragei / ((ie-ib+1)*(jge-jgb+1))            ! this gives the j-averaged temperature at the inlet
          taverager = taverager / (jge-jgb+1)            ! this gives the j-averaged temperature at the recycle plane
        end if
        if (lreadminl==.false.) then
          if (myid==0) then
            write(6,*) 'uaverage(kb)=',uaverage(kb)
            write(6,*) 'uaverage(ke)=',uaverage(ke)
            write(6,*) 'waverage(ke)=',waverage(ke)            
            write(6,*) 'waverage(ke-20)=',waverage(ke-20)            
            write(6,*) 'taveragei(kb)=',taveragei(kb)
            write(6,*) 'taveragei(ke)=',taveragei(ke)
          end if
          
          Utav=0.
          do i=ib,ie
            Utav(i,:) = uaverage
          end do
          
          Uinl = uaverage         ! set the initial time-averaged inlet profile equal to mean u-profile read from means
          write(6,*) 'Uinl(kb+10)=',Uinl(kb+10)
          utaui = sqrt(abs(2*numol* Uinl(kb)/dzf(kb)))     ! average streamwise friction at inlet (need for first time step)
          Urec = uaverage         ! set the initial time-averaged inlet profile equal to mean u-profile

          Wrec(kb:ke+1) = waverage(kb:ke+1)    ! set the initial time-averaged inlet profile equal to mean w-profile
          Wrec(kb) = 0.            ! set the initial time-averaged inlet profile equal to zero
          if (ltempeq == .true.) then
            Ttav=0.
            do i=ib,ie
              Ttav(i,:) = taveragei(:)
            end do
            Tinl = taveragei
            Trec = taveragei
            ttaui = numol*prandtlmoli*2.*(Tinl(kb)-thls)/(dzf(kb)*utaui)  ! friction temp. at inlet (need at first time step)
          end if
        else  ! -> lreadminl==.true. -> Uinl, Urec, Wrec already read
          call slabsum(uaverage ,kb,ke,u0,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,ib,ie,jb,je,kb,ke)
          uaverage = uaverage / ((ie-ib+1)*(jge-jgb+1))  ! this gives the i-j-averaged velocity (only correct for equidistant grid?)            
        end if

   ! determine bulk velocity
        do k=kb,ke
          uaverage(k) = uaverage(k)*dzf(k)
        end do
        ubulk = sum(uaverage(kb:ke))/(zh(ke+1)-zh(kb))      ! volume-averaged u-velocity
        write(6,*) 'Modstartup: ubulk=',ubulk
        
        do k=kb,ke
        do j=jb,je
        uminletbc(j,k) = um(ib,j,k)
        vminletbc(j,k) = vm(ib-1,j,k)
        u0inletbcold(j,k) = um(ib,j,k)
        v0inletbcold(j,k) = vm(ib-1,j,k)
        u0inletbc(j,k) = um(ib,j,k)
        v0inletbc(j,k) = vm(ib-1,j,k)
        end do
        end do 

        do k=kb,ke+1
        do j=jb,je
        wminletbc(j,k) = wm(ib-1,j,k)
        w0inletbcold(j,k) = wm(ib-1,j,k)
        w0inletbc(j,k) = wm(ib-1,j,k)
        end do
        end do
         
        if (ltempeq == .true.) then
          do k=kb,ke
          do j=jb,je
            tminletbc(j,k) = thlm(ib-1,j,k)
            t0inletbcold(j,k) = thlm(ib-1,j,k)
            t0inletbc(j,k) = thlm(ib-1,j,k)
          end do
          end do
        end if 
 
        write(6,*) 'uminletbc(jb,kb),um(ib,jb,kb)=',uminletbc(jb,kb),um(ib,jb,kb)
        write(6,*) 'uminletbc(jb+1,kb+10),um(ib,jb+1,kb+10)=',uminletbc(jb+1,kb+10),um(ib,jb+1,kb+10)
        write(6,*) 'uminletbc(je,kb+10),um(ib,je,kb+10)=',uminletbc(je,kb+10),um(ib,je,kb+10)

      else if (linletgen == 2) then

        nfile = nfile + 1
        write(6,*) 'Loading inletfile' 
        call readinletfile
        u0inletbc(:,:) = storeu0inletbc(:,:,nstepread)
        v0inletbc(:,:) = storev0inletbc(:,:,nstepread)
        w0inletbc(:,:) = storew0inletbc(:,:,nstepread)
        uminletbc(:,:) = storeu0inletbc(:,:,nstepread)
        vminletbc(:,:) = storev0inletbc(:,:,nstepread)
        wminletbc(:,:) = storew0inletbc(:,:,nstepread)
        if (ltempeq == .true.) then
          t0inletbc(:,:) = storet0inletbc(:,:,nstepread)
          tminletbc(:,:) = storet0inletbc(:,:,nstepread)
        end if
        ! determine bulk velocity
        call slabsum(uaverage ,kb,ke,u0,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,ib,ie,jb,je,kb,ke)
        uaverage = uaverage / ((ie-ib+1)*(jge-jgb+1))  ! this gives the i-j-averaged velocity (only correct for equidistant grid?)        
        do k=kb,ke
          uaverage(k) = uaverage(k)*dzf(k)
        end do
        ubulk = sum(uaverage(kb:ke))/(zh(ke+1)-zh(kb))      ! volume-averaged u-velocity
        write(6,*) 'Modstartup: ubulk=',ubulk
      end if  ! linletgen!

      if (lper2inout == .true.) then               ! if the restart starts from a periodic simulation to in/outflow, lper2inout should be set to .true.
        if (myid==0) then
          write(6,*) 'per2inout=.true. -> reading inlet profile from prof.inp.XXX and scalar.inp.XXX'
          open (ifinput,file='prof.inp.'//cexpnr) !  read the inlet profile from prof.inp
          read (ifinput,'(a80)') chmess
          write(*,     '(a80)') chmess
          read (ifinput,'(a80)') chmess

          do k=kb,ke
            read (ifinput,*) &
                  height (k), &
                  thlprof(k), &
                  qtprof (k), &
                  uprof  (k), &
                  vprof  (k), &
                  e12prof(k)
          end do
          svprof = 0.
          if (nsv>0) then
            open (ifinput,file='scalar.inp.'//cexpnr)
            read (ifinput,'(a80)') chmess
            read (ifinput,'(a80)') chmess
            do k=kb,ke
              read (ifinput,*) &
                    height (k), &
                    (svprof (k,n),n=1,nsv)
            end do
            open (ifinput,file='scalar.inp.'//cexpnr)
            write (6,*) 'height   sv(1) --------- sv(nsv) '
            do k=ke,kb,-1
              write (6,*) &
                    height (k), &
                  (svprof (k,n),n=1,nsv)
            end do

          end if
        end if ! end if myid==0

      ! MPI broadcast numbers reading
        call MPI_BCAST(thlprof,kmax,MY_REAL   ,0,comm3d,mpierr)
        call MPI_BCAST(uprof  ,kmax,MY_REAL   ,0,comm3d,mpierr)
        call MPI_BCAST(vprof  ,kmax,MY_REAL   ,0,comm3d,mpierr)
        call MPI_BCAST(e12prof,kmax,MY_REAL   ,0,comm3d,mpierr)
        call MPI_BCAST(qtprof ,kmax,MY_REAL   ,0,comm3d,mpierr)
        call MPI_BCAST(svprof ,(ke+kh-(kb-kh))*nsv,MY_REAL   ,0,comm3d,mpierr)
 
      else if (linoutflow == .true.) then   ! restart of inoutflow simulation: reproduce inlet boundary condition from restartfile
        do j=jb-1,je+1
          do k=kb,ke+1          
            uprof(k)   = u0(ib,j,k)
            vprof(k)   = (v0(ib-1,j,k)+v0(ib,j,k))/2
            thlprof(k) = (thl0(ib-1,j,k)+thl0(ib,j,k))/2
            qtprof(k)  = (qt0(ib-1,j,k)+qt0(ib,j,k))/2
            e12prof(k) = (e120(ib-1,j,k)+e120(ib,j,k))/2
            do n=1,nsv
              svprof(k,n)=(sv0(ib-1,j,k,n)+sv0(ib,j,k,n))/2
            enddo
          enddo
        enddo 
        ! outlet bulk velocity
        call slabsum(uaverage ,kb,ke,u0,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,ib,ie,jb,je,kb,ke)
        uaverage = uaverage / ((ie-ib+1)*(jge-jgb+1))  ! this gives the i-j-averaged velocity (only correct for equidistant grid?)
        ! determine bulk velocity
        do k=kb,ke
          uaverage(k) = uaverage(k)*dzf(k)
        end do
        ubulk = sum(uaverage(kb:ke))/(zh(ke+1)-zh(kb))      ! volume-averaged u-velocity
        write(6,*) 'Modstartup: ubulk=',ubulk
      endif  ! end if lper2inout


      u0av = 0.0
      v0av = 0.0
      thl0av = 0.0
      qt0av   = 0.0
      th0av  = 0.0
      sv0av = 0.

      call slabsum(u0av  ,kb,ke+kh,u0(:,:,kb:ke+kh)  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
      call slabsum(v0av  ,kb,ke+kh,v0(:,:,kb:ke+kh)  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
      call slabsum(thl0av,kb,ke+kh,thl0(:,:,kb:ke+kh),ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
      call slabsum(qt0av,kb,ke+kh,qt0(:,:,kb:ke+kh),ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
      do n=1,nsv
        call slabsum(sv0av(kb,n),kb,ke+kh,sv0(ib-ih,jb-jh,kb,n),ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
      end do

      u0av  = u0av  /rslabs + cu
      v0av  = v0av  /rslabs + cv
      thl0av = thl0av/rslabs
      qt0av = qt0av/rslabs
      sv0av = sv0av /rslabs

      ! CvH - only do this for fixed timestepping. In adaptive dt comes from restartfile
      if(ladaptive .eqv. .false.) dt=dtmax
    !  call boundary
    end if

!-----------------------------------------------------------------
!    2.1 read and initialise fields
!-----------------------------------------------------------------


    if(myid==0)then
      open (ifinput,file='lscale.inp.'//cexpnr)
      read (ifinput,'(a80)') chmess
      read (ifinput,'(a80)') chmess
      write(6,*) ' height  u_geo  v_geo  pgx  pgy  subs     ' &
                    ,'   dqtdx      dqtdy        dqtdtls     thl_rad '
      do  k=kb,ke
        read (ifinput,*) &
              height (k), &
              ug     (k), &
              vg     (k), &
              pgx    (k), &
              pgy    (k), &
              wfls   (k), &
              dqtdxls(k), &
              dqtdyls(k), &
              dqtdtls(k), &
              thlpcar(k) 
      end do
      close(ifinput)

      do k=ke,kb,-1
        write (6,'(3f7.1,5e12.4)') &
              height (k), &
              ug     (k), &
              vg     (k), &
              pgx    (k), &
              pgy    (k), &
              wfls   (k), &
              dqtdxls(k), &
              dqtdyls(k), &
              dqtdtls(k), &
              thlpcar(k)
      end do


    end if ! end myid==0

! MPI broadcast variables read in

    call MPI_BCAST(ug       ,kmax,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(vg       ,kmax,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(pgx      ,kmax,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(pgy      ,kmax,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wfls     ,kmax,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dqtdxls  ,kmax,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dqtdyls  ,kmax,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dqtdtls  ,kmax,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(thlpcar  ,kmax,MY_REAL   ,0,comm3d,mpierr)
    

!-----------------------------------------------------------------
!    2.3 make large-scale horizontal pressure gradient
!-----------------------------------------------------------------

!******include rho if rho = rho(z) /= 1.0 ***********

    do k=kb,ke
      !dpdxl(k) =  om23_gs*vg(k)
      !dpdyl(k) = -om23_gs*ug(k)
      !dpdxl(k) =  -ug(k) 
      !dpdyl(k) =  -vg(k)
      dpdxl(k) =  om23_gs*vg(k)-pgx(k)-dpdx   !corriolis forcing and pressure gradient
      dpdyl(k) =  -om23_gs*ug(k)-pgy(k)
 end do


!-----------------------------------------------------------------
!    2.4 large-scale subsidence, reintroduced ILS13 05.06.2014
!-----------------------------------------------------------------

whls(kb) = 0.0
do k=kb+1,ke
  whls(k) =  (wfls(k)*dzf(k-1) +  wfls(k-1)*dzf(k) )/(2*dzh(k))
end do
whls(ke+1) = (wfls(ke)+0.5*dzf(ke)*(wfls(ke)-wfls(ke-1)) /dzh(ke))



!    idtmax = floor(dtmax/tres)
    btime   = timee
!    timeleft=ceiling(runtime/tres)
    timeleft=runtime
    dt_lim = timeleft
!    write(6,*) 'real(dt)*tres= ',rdt, ' dtmax/100= ',dtmax/100 
    ntrun   = 0
    ntimee  = nint(timee/dtmax)
    tnextrestart = btime + trestart
    tnextrestart2 = btime + trestart2
    tnextfielddump = btime + tfielddump
    deallocate (height,th0av)

!    call boundary
    call bottom   ! this is done to have set the values of 'shear' at the bottom: needed for near-wall damping function at first time step
 
  end subroutine readinitfiles

  subroutine readrestartfiles

    use modsurfdata, only : ustar,thlflux,qtflux,svflux,dudz,dvdz,dthldz,dqtdz,ps,thls,qts,thvs,oblav,&
                           isurf,wtsurf
    use modfields,  only : u0,v0,w0,thl0,qt0,ql0,ql0h,qtav,qlav,e120,dthvdz,presf,presh,sv0,mindist,wall,&
                           uav,vav,wav,uuav,vvav,wwav,uvav,uwav,vwav,svav,thlav,thl2av,sv2av,pres0,svm,&
                           svprof,viscratioav,thluav,thlvav,thlwav,svuav,svvav,svwav,presav,&
                           uusgsav,vvsgsav,wwsgsav,uwsgsav,thlusgsav,thlwsgsav,svusgsav,svwsgsav,tkesgsav,&
                           strain2av,nusgsav
    use modglobal,  only : ib,ie,ih,jb,je,jh,kb,ke,kh,dtheta,dqt,dsv,startfile,timee,totavtime,runavtime,&
                           iexpnr,ntimee,rk3step,ifinput,nsv,nvar,runtime,dt,cu,cv,cexpnr,lreadmean,lreadminl, &
                           totinletav,lreadscal,ltempeq,dzf,numol,prandtlmoli
    use modmpi,     only : cmyid, myid
    use modsubgriddata, only : ekm
    use modinlet,       only : zinterpolate1d,zinterpolatet1d,zinterpolatew1d,zinterpolate2d
    use modinletdata,   only : Uinl,Urec,Wrec,Utav,Tinl,Trec,linuf,linuh,&
                               kbin,kein,lzinzsim,utaui,Ttav,ttaui

    real, dimension(ib:ie,jb:je,kb:ke)  ::  dummy3d
    real, dimension(ib:ie,kbin:kein)    ::  Utavin
    real, dimension(ib:ie,kbin:kein)    ::  Ttavin
    real, dimension(kbin:kein)          ::  Uinlin
    real, dimension(kbin:kein)          ::  Urecin
    real, dimension(kbin:kein)          ::  Tinlin
    real, dimension(kbin:kein)          ::  Trecin
    real, dimension(kbin:kein+1)        ::  Wrecin
    character(50) :: name,name2,name4
    real dummy
    integer i,j,k,n
    !********************************************************************

  !    1.0 Read initfiles
  !-----------------------------------------------------------------
    name = startfile
    name(5:5) = 'd'
    name(15:17)=cmyid
    write(6,*) 'loading ',name
    open(unit=ifinput,file=name,form='unformatted', status='old')

      read(ifinput)  (((mindist(i,j,k),i=ib,ie     ),j=jb,je      ),k=kb,ke   )
      read(ifinput)  ((((wall(i,j,k,n),i=ib,ie     ),j=jb,je      ),k=kb,ke   ),n=1,5)
      read(ifinput)  (((u0    (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      read(ifinput)  (((v0    (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      read(ifinput)  (((w0    (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      read(ifinput)  (((pres0 (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      read(ifinput)  (((thl0  (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      read(ifinput)  (((e120  (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      read(ifinput)  (((ekm   (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      read(ifinput)  (((qt0   (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      read(ifinput)  (((ql0   (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      read(ifinput)  (((ql0h   (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      read(ifinput)  timee,dt
    close(ifinput)
    write(6,*) 'finished loading ',name

    if (nsv>0 .and. lreadscal==.true.) then
      name(5:5) = 's'
      write(6,*) 'loading ',name
      open(unit=ifinput,file=name,form='unformatted')
      read(ifinput) ((((sv0(i,j,k,n),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh),n=1,nsv)
      read(ifinput)  timee
      close(ifinput)
      write(6,*) 'finished loading ',name
    elseif (nsv>0 .and. lreadscal==.false.) then
      sv0 = 0.
      svprof = 0.
    end if

! read mean variables if asked for by lreadmean
    name2 = 'means   .'
    name2( 6:8 ) = cmyid
    name2(10:12) = cexpnr
    if (lreadmean == .true.) then
      write(6,*) 'Reading meansXXX.XXX, proc = ',myid
      open(unit=ifinput,file=name2,form='unformatted')
      read(ifinput)  totavtime, nsv
      read(ifinput)  (((uav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      read(ifinput)  (((vav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      read(ifinput)  (((wav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      read(ifinput)  (((thlav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      read(ifinput)  (((qtav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )    
      read(ifinput)  (((qlav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )    
      read(ifinput)  (((presav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      read(ifinput) ((((svav(i,j,k,n),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   ),n=1,nsv)
      read(ifinput)  (((viscratioav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      read(ifinput)  (((uuav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      read(ifinput)  (((vvav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      read(ifinput)  (((wwav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      read(ifinput)  (((thl2av(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      read(ifinput) ((((sv2av(i,j,k,n),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   ),n=1,nsv)
      read(ifinput)  (((uvav(i,j,k),i=ib,ie+ih  ),j=jb,je+jh      ),k=kb,ke   )
      read(ifinput)  (((uwav(i,j,k),i=ib,ie+ih  ),j=jb,je         ),k=kb,ke+kh   )
      read(ifinput)  (((vwav(i,j,k),i=ib,ie     ),j=jb,je+jh      ),k=kb,ke+kh   )
      read(ifinput)  (((thluav(i,j,k),i=ib,ie     ),j=jb,je         ),k=kb,ke      )
      read(ifinput)  (((thlvav(i,j,k),i=ib,ie     ),j=jb,je+jh      ),k=kb,ke      )
      read(ifinput)  (((thlwav(i,j,k),i=ib,ie     ),j=jb,je         ),k=kb,ke+kh   )
      read(ifinput) ((((svuav(i,j,k,n),i=ib,ie     ),j=jb,je         ),k=kb,ke      ),n=1,nsv)
      read(ifinput) ((((svvav(i,j,k,n),i=ib,ie     ),j=jb,je+jh      ),k=kb,ke      ),n=1,nsv)
      read(ifinput) ((((svwav(i,j,k,n),i=ib,ie     ),j=jb,je         ),k=kb,ke+kh   ),n=1,nsv)
      close(ifinput)
      write(6,*) 'Total averaging time so far: ', totavtime
     
!     read <x'y'>_SGS to file.
      name2 = 'SGS__   .'
      name2( 6:8 ) = cmyid
      name2(10:12) = cexpnr
      open(unit=ifinput,file=name2,form='unformatted')
      read(ifinput)  dummy, dummy
      read(ifinput)  (((uusgsav(i,j,k)   ,i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
      read(ifinput)  (((vvsgsav(i,j,k)   ,i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
      read(ifinput)  (((wwsgsav(i,j,k)   ,i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
      read(ifinput)  (((uwsgsav(i,j,k)   ,i=ib   ,ie+ih),j=jb   ,je   ),k=kb   ,ke+kh)
      read(ifinput)  (((dummy3d(i,j,k)   ,i=ib   ,ie   ),j=jb   ,je   ),k=kb   ,ke   )   ! this is dissresav, which will be computed using other mean quantities
      read(ifinput)  (((tkesgsav(i,j,k)  ,i=ib   ,ie   ),j=jb   ,je   ),k=kb   ,ke   )
      read(ifinput)  (((dummy3d(i,j,k)   ,i=ib   ,ie   ),j=jb   ,je   ),k=kb   ,ke   )   ! this is disssgsav, which will be computed using other mean quantities
      read(ifinput)  (((strain2av(i,j,k) ,i=ib   ,ie   ),j=jb   ,je   ),k=kb   ,ke   )   ! <SijSij> (NOT <Sij><Sij> !!) (average over time)
      read(ifinput)  (((nusgsav(i,j,k)   ,i=ib   ,ie   ),j=jb   ,je   ),k=kb   ,ke   )   ! <nu_sgs> (average over time)
      read(ifinput)  (((thlusgsav(i,j,k) ,i=ib,ie+ih  ),j=jb,je         ),k=kb,ke      )
      read(ifinput)  (((thlwsgsav(i,j,k) ,i=ib,ie     ),j=jb,je         ),k=kb,ke+kh   )
      read(ifinput) ((((svusgsav(i,j,k,n),i=ib,ie+ih  ),j=jb,je         ),k=kb,ke      ),n=1,nsv)
      read(ifinput) ((((svwsgsav(i,j,k,n),i=ib,ie     ),j=jb,je         ),k=kb,ke+kh   ),n=1,nsv)
      close (ifinput)

    end if
 
    
! read mean profiles for inlet generator       
    if (lreadminl==.true.) then
      if (lzinzsim==.false.) then 
        name4 = 'meaninlet.   '
        name4(11:13) = cexpnr
        open (unit=ifinput,file=name4,form='unformatted')
        read(ifinput) totinletav      ! interval of time-average
        read(ifinput) (Uinlin(k),k=kbin,kein   )
        read(ifinput) (Urecin(k),k=kbin,kein   )
        read(ifinput) (Wrecin(k),k=kbin,kein+1   )
        read(ifinput) ((Utavin(i,k),i=ib,ie ),k=kbin,kein   )
        close(ifinput)

        call zinterpolate1d (Uinlin,Uinl)   ! interpolate inlet profile to zgrid
        call zinterpolate1d (Urecin,Urec)  
        call zinterpolatew1d(Wrecin,Wrec)   
        call zinterpolate2d (Utavin,Utav)   

        if (ltempeq == .true.) then
          name4 = 'tempinlet.   '
          name4(11:13) = cexpnr
          open (unit=ifinput,file=name4,form='unformatted')
          read(ifinput) totinletav      ! interval of time-average
          read(ifinput) (Tinlin(k),k=kbin,kein   )
          read(ifinput) (Trecin(k),k=kbin,kein   )
          read(ifinput) ((Ttavin(i,k),i=ib,ie ),k=kbin,kein   )
          close(ifinput)

          call zinterpolatet1d(Tinlin,Tinl)  
          call zinterpolatet1d(Trecin,Trec)  
          call zinterpolate2d (Ttavin,Ttav)   
        end if ! ltempeq

      else  !lzinzsim=.true. -> inlet grid equals sim grid
        name4 = 'meaninlet.   '
        name4(11:13) = cexpnr
        open (unit=ifinput,file=name4,form='unformatted')
        read(ifinput) totinletav      ! interval of time-average
        read(ifinput) (Uinl(k),k=kb,ke   )
        read(ifinput) (Urec(k),k=kb,ke   )
        read(ifinput) (Wrec(k),k=kb,ke+1   ) 
        read(ifinput) ((Utav(i,k),i=ib,ie ),k=kb,ke   )
           
        close(ifinput)
  
        if (ltempeq == .true.) then
          name4 = 'tempinlet.   '
          name4(11:13) = cexpnr
          open (unit=ifinput,file=name4,form='unformatted')
          read(ifinput) totinletav      ! interval of time-average
          read(ifinput) (Tinl(k),k=kb,ke   )
          read(ifinput) (Trec(k),k=kb,ke   )
          read(ifinput) ((Ttav(i,k),i=ib,ie ),k=kb,ke   )

          close(ifinput)
        end if ! ltempeq
      end if ! lzinzsim
      
      utaui = sqrt(abs(2*numol* Uinl(kb)/dzf(kb)))     ! average streamwise friction at inlet (need for first time step)
      if (ltempeq == .true.) then
        ttaui = numol*prandtlmoli*2.*(Tinl(kb)-thls)/(dzf(kb)*utaui)
      end if
    end if !(lreadminl==.true.)


  end subroutine readrestartfiles

  subroutine exitmodules
    use modfields,         only : exitfields
    use modglobal,         only : exitglobal
    use modmpi,            only : exitmpi
    use modboundary,       only : exitboundary
    use modpois,           only : exitpois
    use modsubgrid,        only : exitsubgrid
    use modsurface,        only : exitsurface
    use modthermodynamics, only : exitthermodynamics
    use modinlet,          only : exitinlet

    call exitthermodynamics
    call exitsurface
    call exitsubgrid
    call exitpois
    call exitboundary
    call exitfields
    call exitglobal
    call exitinlet
    call exitmpi

 end subroutine exitmodules
!----------------------------------------------------------------
  subroutine randomnize(field,klev,ampl,ir,ihl,jhl)

    use modmpi,    only :  myid,nprocs
    use modglobal, only : ib,ie,imax,jmax,jb,je,kb,ke,kh
    integer (KIND=selected_int_kind(6)):: imm, ia, ic,ir
    integer ihl, jhl
    integer i,j,klev
    integer m,mfac
    real ran,ampl
    real field(ib-ihl:ie+ihl,jb-jhl:je+jhl,kb-kh:ke+kh)
    parameter (imm = 134456, ia = 8121, ic = 28411)

    if (myid>0) then
      mfac = myid * jmax * imax
      do m =1,mfac
        ir=mod((ir)*ia+ic,imm)

      end do
    end if
    do j=jb,je
    do i=ib,ie
      ir=mod((ir)*ia+ic,imm)
      ran=real(ir)/real(imm)
      field(i,j,klev) = field(i,j,klev) + (ran-0.5)*2.0*ampl
    end do
    end do

    if (nprocs-1-myid > 0) then
      mfac = (nprocs-1-myid) * imax * jmax
      do m=1,mfac
        ir=mod((ir)*ia+ic,imm)
      end do
    end if

    return
  end subroutine randomnize

  subroutine readscalpointsource  
! 12-12-2014: Jasper Tomas
!
! This routine eventually only produces scalpointsrc(i,j,k,1:nsvp) that gives the 
! scalar source distribution on this proc.
! The routine reads scalptloc.inp.XXX which contains the locations and standard
! deviations of the point sources. It sets the source to zero when an obstacle is
! present and rescales the total source to 1.

    use modglobal, only :  ifinput,cexpnr,scalptloc,scalptsrc,nsv,ib,ie,jb,&
                           je,kb,ke,xf,zf,xh,zh,dxf,dzf,nblocks,libm,&
                           nsvl,nsvp,scalptsrc,dy,jmax
    use modibm,    only :  block
    use modmpi,    only :  myid,MY_REAL,comm3d,mpierr,MPI_SUM
    real scalsum
    integer i,j,k,n,il,iu,jl,ju,kl,ku
    character(80) chmess
    real, dimension(1:nsvp) :: scalptsuml ! sum of scalptsrc on this proc
    real, dimension(1:nsvp) :: scalptsum  ! total sum of scalptsrc 

    allocate(scalptloc(1:nsvp,4))
    allocate(scalptsrc(ib:ie,jb:je,kb:ke,1:nsvp)) ! local point source distr.

    scalptsrc(:,:,:,:) = 0.

  if(nsvp > 0) then
    if(myid==0) then
      open (ifinput,file='scalptloc.inp.'//cexpnr)
      read (ifinput,'(a80)') chmess
      read (ifinput,'(a80)') chmess
      do n=1,nsvp
        read (ifinput,*) &
              scalptloc(n,1), &        ! x-location of scalar source
              scalptloc(n,2), &        ! y-location of scalar source
              scalptloc(n,3), &        ! z-location of scalar source
              scalptloc(n,4)           ! standard deviation (in all directions)
      end do
      write (6,*) 'Point source number,  x-loc, y-loc, z-loc, sigma'
      do n=1,nsvp
        write (6,*) &
              n , &
              scalptloc(n,1), &        ! x-location of scalar source
              scalptloc(n,2), &        ! y-location of scalar source
              scalptloc(n,3), &        ! z-location of scalar source
              scalptloc(n,4)           ! standard deviation (in all directions)
      end do
    end if ! end if myid==0
    call MPI_BCAST(scalptloc ,4*nsvp,MY_REAL ,0,comm3d,mpierr)   ! share with other procs

    do n=1,nsvp
    do k=kb,ke
    do j=jb,je
    do i=ib,ie
      scalptsrc(i,j,k,n) = exp( - ( ((xf(i)-scalptloc(n,1))**2)/(2.*(scalptloc(n,4)**2)) &
                                  +((((j+myid*jmax)*dy-0.5*dy)-scalptloc(n,2))**2)/(2.*(scalptloc(n,4)**2)) &
                                  +((zf(k)-scalptloc(n,3))**2)/(2.*(scalptloc(n,4)**2))))
    end do
    end do
    end do
    end do

    ! set scalptsrc to zero when an obstacle is inside
    if (libm==.true.) then
      do n=1,nblocks
        il = block(n,1)
        iu = block(n,2)
        kl = block(n,5)
        ku = block(n,6)
        jl = block(n,3)-myid*jmax
        ju = block(n,4)-myid*jmax
        if (ju < jb .or. jl > je) then
          cycle
        else
          if (ju > je) ju=je
          if (jl < jb) jl=jb
          scalptsrc(il:iu,jl:ju,kl:ku,:) = 0.
        end if
      end do
    end if
    
    ! compute the sum of scalptsrc on this proc
    scalptsuml(:) = 0.0
    do n=1,nsvp
      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        scalptsuml(n) = scalptsuml(n) + dxf(i)*dy*dzf(k)*scalptsrc(i,j,k,n) 
      end do
      end do
      end do
    end do
  call MPI_ALLREDUCE(scalptsuml(1:nsvp),scalptsum(1:nsvp),nsvp,MY_REAL,MPI_SUM,comm3d,mpierr)
!    write(6,*) 'scalptsuml,scalptsum = ', scalptsuml(1),scalptsum(1) 
    do n=1,nsvp
      scalptsrc(:,:,:,n) = scalptsrc(:,:,:,n)/scalptsum(n) ! normalize the source term
!      write(6,*) 'sum(scalptsrc) = ', sum(sum(sum(scalptsrc(:,:,:,n),1),1),1)
    end do
  end if ! nsvp>0   
!    ! test 
!    scalptsuml(:) = 0.0
!    do n=1,nsvp
!      do k=kb,ke
!      do j=jb,je
!      do i=ib,ie
!        scalptsuml(n) = scalptsuml(n) + dxf(i)*dy*dzf(k)*scalptsrc(i,j,k,n) 
!      end do
!      end do
!      end do
!      write(6,*) 'total point source [kg/sec] = ',scalptsuml(n)
!    end do

  end subroutine readscalpointsource  

  subroutine readscalsource
! 31-03-2014: Jasper Tomas
!
! This routine eventually only produces scalsrc(i,k,nsv) that gives the scalar
! source distribution in the x-z plane.
! The routine reads scalarloc.inp.XXX which contains the locations and standard
! deviations of the line sources. It sets the source to zero when an obstacle is
! present (anywhere along j) and rescales the total source to 1.
    use modglobal, only :  ifinput,cexpnr,scalloc,scalsrc,nsv,ib,ie,jb,& 
                           je,kb,ke,xf,zf,xh,zh,dxf,dzf,nblocks,libm,nsvl,&
                           nsvp,jtot,ysize
    use modibm,    only :  block
    use modmpi,    only :  myid,MY_REAL,comm3d,mpierr
    integer i,j,k,n,nb,m
    character(80) chmess
    real scalsum
    allocate(scalloc(1:nsv,4))
    allocate(scalsrc(ib:ie,kb:ke,1:nsv))

    scalsrc(:,:,:) = 0.

    if(myid==0 .and. nsv > 0) then
      open (ifinput,file='scalarloc.inp.'//cexpnr)
      read (ifinput,'(a80)') chmess
      read (ifinput,'(a80)') chmess
      do n=1,nsvl
        read (ifinput,*) &
              scalloc(n,1), &        ! x-location of scalar source
              scalloc(n,2), &        ! z-location of scalar source
              scalloc(n,3), &        ! standard deviation in x
              scalloc(n,4)           ! standard deviation in z
      end do
      write (6,*) 'Line source number,  x-loc, z-loc, sigma_x, sigma_z'
      do n=1,nsvl
        write (6,*) &
              n , &
              scalloc(n,1), &        ! x-location of scalar source
              scalloc(n,2), &        ! z-location of scalar source
              scalloc(n,3), &        ! standard deviation in x
              scalloc(n,4)           ! standard deviation in z
      end do
    end if ! end if myid==0
    call MPI_BCAST(scalloc ,4*nsvl,MY_REAL ,0,comm3d,mpierr)   ! share with other procs
     
    do n=1,nsvl
    do k=kb,ke
    do i=ib,ie

    ! commented by tg3315 !undone
    scalsrc(i,k,n) = exp( - ( ((xf(i)-scalloc(n,1))**2)/(2.*(scalloc(n,3)**2)) & 
                             +((zf(k)-scalloc(n,2))**2)/(2.*(scalloc(n,4)**2)) )) 
    end do
    end do
    end do

    ! set scalsrc to zero when an obstacle is inside (aywhere along the
    ! j-direction)
    if (libm==.true.) then
      do nb=1,nblocks ! loop over blocks
      do k=block(nb,5),block(nb,6)
      do i=block(nb,1),block(nb,2)
        scalsrc(i,k,:) = 0.   ! set to zero for all scalars (no nsv loop necessary)
      end do
      end do
      end do
      end if

! Rescale the total source to 1.
    do n=1,nsvl
      scalsum = 0.0
      do k=kb,ke
      do i=ib,ie
!        scalsum = scalsum + dxf(i)*dzf(k)*scalsrc(i,k,n)  ! this is probably thecorrect way...
        scalsum = scalsum + dxf(i)*ysize*dzf(k)*scalsrc(i,k,n)  ! this is probably thecorrect way...
      end do
      end do
      scalsrc(:,:,n) = scalsrc(:,:,n)/scalsum     ! correct the source to have a total of 1.
    end do

!Added by tg3315 !undone
   ! scalsrc = 0

  end subroutine readscalsource


end module modstartup
