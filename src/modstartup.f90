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

   integer(KIND=selected_int_kind(6)) :: irandom = 0 !    * number to seed the randomnizer with
   integer :: krand = huge(0)  ! returns the largest integer that is not an infinity
   real :: randu = 0.0, randthl = 0.0, randqt = 0.0 !    * uvw,thl and qt amplitude of randomnization

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

      use modglobal, only:initglobal, iexpnr, runtime, dtmax,  &
         lwarmstart, lstratstart, lfielddump, lreadscal, startfile, tfielddump, fieldvars, tsample, tstatsdump, trestart, &
         nsv, imax, jtot, kmax, xsize, ysize, xlat, xlon, xday, xtime, lwalldist, &
         lmoist, lcoriol, igrw_damp, geodamptime, ifnamopt, fname_options, &
         xS,yS,zS,SS,sigS,iwallmom,iwalltemp,iwallmoist,iwallscal,ipoiss,iadv_mom,iadv_tke,iadv_thl,iadv_qt,iadv_sv,courant,diffnr,ladaptive,author,&
         linoutflow, lper2inout, libm, lnudge, tnudge, nnudge, lles, luoutflowr, lvoutflowr, luvolflowr, lvvolflowr, &
         uflowrate, vflowrate, lstoreplane, iplane, &
         lreadmean, iinletgen, inletav, lreadminl, Uinf, Vinf, linletRA, nblocks, &
         lscalrec,lSIRANEinout,lscasrc,lscasrcl,lscasrcr,lydump,lytdump,lxydump,lxytdump,lslicedump,ltdump,ltkedump,lzerogradtop,&
         lzerogradtopscal, lbuoyancy, ltempeq, &
         lfixinlet, lfixutauin, pi, &
         thlsrc, ifixuinf, lvinf, tscale, ltempinout, lmoistinout,  &
         lwallfunc,lprofforc,lchem,k1,JNO2,rv,rd,tnextEB,tEB,dtEB,bldT,wsoil,wgrmax,wwilt,wfc,skyLW,GRLAI,rsmin,nfcts,lEB,lwriteEBfiles,nwalllayers,lconstW, &
         BCxm,BCxT,BCxq,BCxs,BCym,BCyT,BCyq,BCys, &
         BCtopm,BCtopT,BCtopq,BCtops,BCbotm,BCbotT,BCbotq,BCbots, &
         idriver,tdriverstart,driverjobnr,dtdriver,driverstore
      use modsurfdata, only:z0, z0h,  wtsurf, wttop, wqtop, wqsurf, wsvsurf, wsvtop, wsvsurfdum, wsvtopdum, ps, thvs, thls, thl_top, qt_top, qts 
      ! use modsurface,        only : initsurface
      use modfields, only:initfields, dpdx, ncname
      use modpois, only:initpois
      use modboundary, only:initboundary, ksp
      use modthermodynamics, only:initthermodynamics, lqlnr, chi_half
      use modsubgrid, only:initsubgrid
      use modmpi, only:comm3d, myid, mpi_integer, mpi_logical, my_real, mpierr, mpi_character
      use modinlet, only:initinlet
      use modinletdata, only:di, dr, di_test, dti, iangledeg, iangle
      use modibmdata, only:bctfxm, bctfxp, bctfym, bctfyp, bctfz
      use modforces, only: calcfluidvolumes
      use moddriver, only: initdriver

      implicit none
      integer :: ierr

      !declare namelists

      namelist/RUN/ &
         iexpnr, lwarmstart, lstratstart, startfile, &
         runtime, dtmax, trestart, ladaptive, &
         irandom, randu, randthl, randqt, krand, &
         courant, diffnr, author, &
         libm, lles, &
         lper2inout, lwalldist, &
         lreadmean
      namelist/DOMAIN/ &
         imax, jtot, kmax, xsize, ysize, &
         xlat, xlon, xday, xtime, ksp 
      namelist/PHYSICS/ &
         ps, igrw_damp, lmoist, lcoriol, lbuoyancy, ltempeq, &
         lprofforc, ifixuinf, lvinf, tscale, dpdx, &
         luoutflowr, lvoutflowr, luvolflowr, lvvolflowr, &
         uflowrate, vflowrate, &
         lnudge, tnudge, nnudge
      namelist/DYNAMICS/ &
         lqlnr, ipoiss, &
         iadv_mom, iadv_tke, iadv_thl, iadv_qt, iadv_sv
      namelist/BC/ &
         BCxm, BCxT, BCxq, BCxs, &
         BCym, BCyT, BCyq, BCys, &
         BCtopm, BCtopT, BCtopq, BCtops, &
         BCbotm, BCbotT, BCbotq, BCbots, &
         bctfxm, bctfxp, bctfym, bctfyp, bctfz, &
         wttop, thl_top, qt_top, qts, wsvsurfdum, wsvtopdum, &
         wtsurf, wqsurf, thls, z0, z0h
      namelist/INLET/ &
         Uinf, Vinf, di, dti, inletav, linletRA, &
         lstoreplane, lreadminl, lfixinlet, lfixutauin, &
         lwallfunc
      namelist/DRIVER/ &
         idriver, tdriverstart, driverjobnr, dtdriver, &
         driverstore, iplane
      namelist/WALLS/ &
         nblocks, nfcts, iwallmom, iwalltemp, iwallmoist, iwallscal
      namelist/ENERGYBALANCE/ &
         lEB, lwriteEBfiles, lconstW, dtEB, bldT, wsoil, wgrmax, wwilt, wfc, &
         skyLW, GRLAI, rsmin, nwalllayers
      namelist/SCALARS/ &
         lreadscal, lscasrc, lscasrcl, lscasrcr, &
         nsv, xS, yS, zS, SS, sigS
      namelist/CHEMISTRY/ &
         lchem, k1, JNO2
      namelist/OUTPUT/ &
         lfielddump, tfielddump, fieldvars, &
         ltdump, lydump, lytdump, lxydump, lxytdump, &
         lslicedump, ltkedump, tstatsdump, tsample

      if (myid == 0) then
         if (command_argument_count() >= 1) then
            call get_command_argument(1, fname_options)
         end if
         write (*, *) fname_options

         open (ifnamopt, file=fname_options, status='old', iostat=ierr)
         if (ierr /= 0) then
            write(0, *) 'ERROR: Namoptions does not exist'
            write(0, *) 'iostat error: ', ierr
            stop 1
         end if

         read (ifnamopt, RUN, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions RUN'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         write (6, RUN)
         rewind (ifnamopt)

         read (ifnamopt, DOMAIN, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions DOMAIN'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         write (6, DOMAIN)
         rewind (ifnamopt)

         read (ifnamopt, PHYSICS, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions PHYSICS'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         write (6, PHYSICS)
         rewind (ifnamopt)

         read (ifnamopt, DYNAMICS, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions DYNAMICS'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         write (6, DYNAMICS)
         rewind (ifnamopt)

         read (ifnamopt, BC, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions BC'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         write (6, BC)
         rewind (ifnamopt)

         read (ifnamopt, INLET, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions INLET'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         write (6, INLET)
         rewind (ifnamopt)

         read (ifnamopt, DRIVER, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'Problem in namoptions DRIVER'
            write(0, *) 'iostat error: ', ierr
            stop 'ERROR: Problem in namoptions DRIVER'
         endif
         write (6, DRIVER)
         rewind (ifnamopt)

         read (ifnamopt, WALLS, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions WALLS'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         write (6, WALLS)
         rewind (ifnamopt)

         read (ifnamopt, ENERGYBALANCE, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions EB'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         write (6, ENERGYBALANCE)
         rewind (ifnamopt)

         read (ifnamopt, SCALARS, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions SCALARS'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         write (6, SCALARS)
         rewind (ifnamopt)

         read (ifnamopt, CHEMISTRY, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions CHEMISTRY'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         write (6, CHEMISTRY)
         rewind (ifnamopt)

         read (ifnamopt, OUTPUT, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions OUTPUT'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         write (6, OUTPUT)
         close (ifnamopt)
      end if

      thvs = thls*(1.+(rv/rd - 1.)*qts)

      write (*, *) "starting broadcast"
      !broadcast namelists
      call MPI_BCAST(iexpnr, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(lwarmstart, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lstratstart, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lfielddump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lreadscal, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch to read scalar pollutant fields (warm start)
      call MPI_BCAST(lscasrc, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! tg3315
      call MPI_BCAST(lscasrcl, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! tg3315
      call MPI_BCAST(lscasrcr, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! tg3315

      call MPI_BCAST(lbuoyancy, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for buoyancy force in modforces
      call MPI_BCAST(ltempeq, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for solving adv/diff equation for temperature
      call MPI_BCAST(lper2inout, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for restart periodic flow to inoutflow
      call MPI_BCAST(libm, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for turning on/off IBM method
      call MPI_BCAST(lnudge, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(nnudge, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(tnudge, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(lwalldist, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for computing wall distances
      call MPI_BCAST(lles, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for turning on/off LES functionality (subgrid model)
      call MPI_BCAST(linletRA, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for turning on/off Running Average in inletgenerator
      call MPI_BCAST(lfixinlet, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for keeping average inlet velocit and temp fixed at inlet (iinletgen=1,2)
      call MPI_BCAST(lfixutauin, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for keeping utau fixed at inlet (iinletgen=1,2)
      call MPI_BCAST(xS, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(yS, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(zS, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(SS, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(sigS, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(idriver    ,1,MPI_INTEGER,0,comm3d,mpierr)        ! ae1212: Added switch for driver inlet simulation                                                                         
      call MPI_BCAST(tdriverstart,1,MY_REAL   ,0,comm3d,mpierr)        ! ae1212
      call MPI_BCAST(driverjobnr,1,MPI_INTEGER,0,comm3d,mpierr)        ! ae1212
      call MPI_BCAST(dtdriver   ,1,MY_REAL    ,0,comm3d,mpierr)        ! ae1212
      call MPI_BCAST(driverstore,1,MPI_INTEGER ,0,comm3d,mpierr)
      write (*, *) "sec BC"
         call MPI_BCAST(BCxm, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCxT, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCxq, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCxs, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCym, 1, MPI_INTEGER, 0, comm3d, mpierr)                                                                                                                                                                             
         call MPI_BCAST(BCyT, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCyq, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCys, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCtopm, 1, MPI_INTEGER, 0, comm3d, mpierr)                                                         
         call MPI_BCAST(BCtopT, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCtopq, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCtops, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCbotm, 1, MPI_INTEGER, 0, comm3d, mpierr) 
         call MPI_BCAST(BCbotT, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCbotq, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCbots, 1, MPI_INTEGER, 0, comm3d, mpierr)

      write (*, *) "sec c"

      call MPI_BCAST(lwallfunc, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for reading mean inlet/recycle plane profiles (Uinl,Urec,Wrec)
      call MPI_BCAST(lreadminl, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for reading mean inlet/recycle plane profiles (Uinl,Urec,Wrec)
      call MPI_BCAST(iwalltemp, 1, MPI_INTEGER, 0, comm3d, mpierr) ! case (integer) for wall treatment for temperature (1=no wall function/fixed flux, 2=no wall function/fixed value, 3=uno)
      call MPI_BCAST(iwallmoist, 1, MPI_INTEGER, 0, comm3d, mpierr) ! case (integer) for wall treatment for moisture (1=no wall function/fixed flux, 2=no wall function/fixed value, 3=uno)
      call MPI_BCAST(iwallscal, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iwallmom, 1, MPI_INTEGER, 0, comm3d, mpierr) ! case (integer) for wall treatment for momentum (1=no wall function, 2=werner-wengle, 3=uno)
      write (*, *) "sec d"
      call MPI_BCAST(luoutflowr, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for turning on/off u-velocity correction for fixed mass outflow rate
      call MPI_BCAST(lvoutflowr, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! tg3315: added switch for turning on/off v-velocity correction for fixed mass outflow rate
      call MPI_BCAST(luvolflowr, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! bss166: added switch for turning on/off u-velocity correction for fixed volume flow rate
      call MPI_BCAST(lvvolflowr, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! bss116: added switch for turning on/off v-velocity correction for fixed volume flow rate
      call MPI_BCAST(lstoreplane, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for turning on/off for storing i-plane data to serve as inlet for future sim.
      call MPI_BCAST(lreadmean, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for reading mean variables from means#MYID#.#EXPNR#
      call MPI_BCAST(lydump, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! tg3315 added switch for writing statistics files
      call MPI_BCAST(lytdump, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! tg3315 added switch for writing statistics files
      call MPI_BCAST(lxydump, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! tg3315 added switch for writing statistics files
      call MPI_BCAST(lxytdump, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! tg3315 added switch for writing statistics files
      call MPI_BCAST(lslicedump, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! tg3315 added switch for writing statistics files
      call MPI_BCAST(ltdump, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! tg3315 added switch for writing statistics files
      call MPI_BCAST(ltkedump, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! tg3315 added switch for writing tke budget files
      call MPI_BCAST(iplane, 1, MPI_INTEGER, 0, comm3d, mpierr) ! J.Tomas: ib+iplane is the i-plane that is stored if lstoreplane is .true.
      call MPI_BCAST(startfile, 50, MPI_CHARACTER, 0, comm3d, mpierr)
      call MPI_BCAST(author, 80, MPI_CHARACTER, 0, comm3d, mpierr)
      call MPI_BCAST(runtime, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(trestart, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(tfielddump, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(tsample, 1, MY_REAL, 0, comm3d, mpierr) !tg3315
      call MPI_BCAST(tstatsdump, 1, MY_REAL, 0, comm3d, mpierr) !tg3315

      call MPI_BCAST(tEB, 1, MY_REAL, 0, comm3d, mpierr)
      tnextEB = dtEB
      call MPI_BCAST(tnextEB, 1, MY_REAL, 0, comm3d, mpierr)

      call MPI_BCAST(dtmax, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(nsv, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(fieldvars, 50, MPI_CHARACTER, 0, comm3d, mpierr)
      !call MPI_BCAST(nstat      ,1,MPI_INTEGER,0,comm3d,mpierr) !tg3315
      !call MPI_BCAST(ncstat     ,80,MPI_CHARACTER,0,comm3d,mpierr) !tg3315
      call MPI_BCAST(ifixuinf, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(lvinf, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(dpdx, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(tscale, 1, MY_REAL, 0, comm3d, mpierr)

      call MPI_BCAST(imax, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(jtot, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(kmax, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(xsize, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(ysize, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(xlat, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(xlon, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(xday, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(xtime, 1, MY_REAL, 0, comm3d, mpierr)
      write (*, *) "sec f"

      call MPI_BCAST(z0, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(z0h, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctfxm, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctfxp, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctfym, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctfyp, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctfz, 1, MY_REAL, 0, comm3d, mpierr)

      call MPI_BCAST(wtsurf, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wqsurf, 1, MY_REAL, 0, comm3d, mpierr)
      allocate (wsvsurf(1:nsv))
      wsvsurf = wsvsurfdum(1:nsv)
      call MPI_BCAST(wsvsurf(1:nsv), nsv, MY_REAL, 0, comm3d, mpierr)
      allocate (wsvtop(1:nsv))
      wsvtop = wsvtopdum(1:nsv)                                                                   
      call MPI_BCAST(wsvtop(1:nsv), nsv, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(ps, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(thvs, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(thls, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(thl_top, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(qt_top, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(qts, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(lmoist, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lcoriol, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lprofforc, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lchem, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(k1, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(JNO2, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(igrw_damp, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(geodamptime, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wttop, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wqtop, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(thlsrc, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(uflowrate, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(vflowrate, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(Uinf, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(Vinf, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(di, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(dti, 1, MY_REAL, 0, comm3d, mpierr)
      dr = di ! initial value is needed
      di_test = di ! initial value is needed
      write (*, *) "sec g"
      call MPI_BCAST(iangledeg, 1, MY_REAL, 0, comm3d, mpierr)
      iangle = iangledeg*pi/180.
      call MPI_BCAST(inletav, 1, MY_REAL, 0, comm3d, mpierr)

      call MPI_BCAST(lqlnr, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ksp, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nblocks, 1, MPI_INTEGER, 0, comm3d, mpierr) ! no. of blocks used in IBM

      call MPI_BCAST(nfcts, 1, MPI_INTEGER, 0, comm3d, mpierr)

      call MPI_BCAST(lconstW, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lEB, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lwriteEBfiles, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(wsoil, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wgrmax, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wwilt, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wfc, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(dtEB, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bldT, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(skyLW, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(GRLAI, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(rsmin, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(nwalllayers, 1, MPI_INTEGER, 0, comm3d, mpierr)

      call MPI_BCAST(irandom, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(krand, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(randthl, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(randu, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(randqt, 1, MY_REAL, 0, comm3d, mpierr)

      call MPI_BCAST(ladaptive, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(courant, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(diffnr, 1, MY_REAL, 0, comm3d, mpierr)

      call MPI_BCAST(ipoiss, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iadv_mom, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iadv_tke, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iadv_thl, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iadv_qt, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iadv_sv(1:nsv), nsv, MPI_INTEGER, 0, comm3d, mpierr)
      !write(*,*) "sec h"

      ! Allocate and initialize core modules
      call initglobal
      write (*, *) "done initglobal"
      call initfields
      write (*, *) "done initfields"
      call initboundary
      write (*, *) "done initboundayi"
      call initthermodynamics
      write (*, *) "done initthermodynamics"
      !!depreated!!
      ! call initsurface
      write (*, *) "done initsurface"
      call initsubgrid
      write (*, *) "done initsubgrid"
      call initpois
      write (*, *) "done initpois"
      call initinlet ! added by J. Tomas: initialize inlet generator
      write (*, *) "done initinlet"
      call initdriver  ! added by ae1212: initialise driver inlet                             
      write(*,*) "done initdriver"
      call checkinitvalues
      write (*, *) "done checkinitvalues"
      write (6, *) 'Determine masking matrices'
      call createmasks ! determine walls/blocks
      write (6, *) 'Finished determining masking matrices'
      ! calculate fluid volume and outlet areas, needs masking matrices
      call calcfluidvolumes

      call readinitfiles
      write (*, *) "done readinitfiles"
      write (*, *) "done startup"

      call createscals
      write (*, *) "done create scals"

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

      use modsurfdata, only : wtsurf, wqsurf, qts, ps
      use modglobal, only   : imax,kmax,jtot,ysize,xsize,dxf,ib,ie,&
                              dtmax,runtime,startfile,lwarmstart,lstratstart,&
                              BCxm,BCxT,BCxq,BCxs,BCtopm,BCbotm,&
                              iinletgen,linoutflow,ltempeq,iwalltemp,iwallmom,&
                              ipoiss,POISS_FFT,POISS_CYC
      use modmpi, only      : myid, nprocs, mpierr, comm3d, MPI_INTEGER, MPI_LOGICAL
      use modglobal, only   : idriver
      implicit none
      real :: d(1:imax-1)
      logical :: inequi

      if (mod(jtot, nprocs) /= 0) then
         if (myid == 0) then
            write (0, *) 'STOP ERROR IN NUMBER OF PROCESSORS'
            write (0, *) 'nprocs must divide jtot!!! '
            write (0, *) 'nprocs and jtot are: ', nprocs, jtot
         end if
         call MPI_FINALIZE(mpierr)
         stop 1
      end if

      if (ipoiss==POISS_FFT) then
        if(mod(imax,nprocs)/=0)then
          if(myid==0)then
            write(0,*)'STOP ERROR IN NUMBER OF PROCESSORS'
            write(0,*)'nprocs must divide imax!!! '
            write(0,*)'nprocs and imax are: ',nprocs,imax
          end if
          call MPI_FINALIZE(mpierr)
          stop 1
        end if
      end if

      if (mod(kmax, nprocs) /= 0) then
         if (myid == 0) then
            write (0, *) 'STOP ERROR IN NUMBER OF PROCESSORS'
            write (0, *) 'nprocs must divide kmax!!! '
            write (0, *) 'nprocs and kmax are: ', nprocs, kmax
         end if
         call MPI_FINALIZE(mpierr)
         stop 1
      end if

      !Check Namoptions
      if (runtime < 0) then
         write(0, *) 'ERROR: runtime out of range/not set'
         stop 1
      end if
      if (dtmax < 0) then
         write(0, *) 'ERROR: dtmax out of range/not set'
         stop 1
      end if
      if (ps < 0) then
         write(0, *) 'ERROR: psout of range/not set'
         stop 1
      end if
      if (xsize < 0) then
         write(0, *) 'ERROR: xsize out of range/not set'
         stop 1
      end if
      if (ysize < 0) then
         write(0, *) 'ERROR: ysize out of range/not set'
         stop 1
      end if

      if ((lwarmstart) .or. (lstratstart)) then
         if (startfile == '') then 
            write(0, *) 'ERROR: no restartfile set'
            stop 1
         end if
      end if

      ! Switch to ensure that neutral wall function is called when ltempeq=false and if iwalltemp==1 (constant flux and therefore wall temp is not resolved.
      if ((ltempeq .eqv. .false.) .or. (iwalltemp==1)) then
         iwallmom = 3
         BCbotm = 3
      end if

      ! choosing inoutflow in x requires switches to be set
      ! tg3315 - these could be moved to init boundary
      if (BCxm .eq. 2) then
         write (*, *) "inoutflow conditions, setting appropriate switches (1)"
         iinletgen = 1
         BCxT = 3 !temperature is considered in inletgen & iolet
         BCxq = 3 !humidity is considered in iolet
         BCxs = 3 !scalars are considered in iolet
         BCtopm = 3 !velocity at top determined by topm
         linoutflow = .true.
         call MPI_BCAST(iinletgen, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCxT, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCxq, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCxs, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCtopm, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(linoutflow, 1, MPI_LOGICAL, 0, comm3d, mpierr)

      else if (BCxm .eq. 3) then
         write (*, *) "inoutflow conditions, setting appropriate switches (2)"

         iinletgen = 2
         ! see modstartup for conditions that apply with inletgenerators
         ! move to modstartup
         BCxT = 3 !temperature is considered in inletgen & iolet
         BCxq = 3 !humidity is considered in iolet
         BCxs = 3 !scalars are considered in iolet
         BCtopm = 3 !velocity at top determined by topm
         linoutflow = .true.
         call MPI_BCAST(iinletgen, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCxT, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCxq, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCxs, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCtopm, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(linoutflow, 1, MPI_LOGICAL, 0, comm3d, mpierr)

      else if (BCxm .eq. 4) then
         write (*, *) "inoutflow conditions, setting appropriate switches (0)"

         iinletgen = 0
         ! see modstartup for conditions that apply with inletgenerators
         ! move to modstartup
         BCxT = 3 !temperature is considered in inletgen & iolet
         BCxq = 3 !humidity is considered in iolet
         BCxs = 3 !scalars are considered in iolet
         BCtopm = 3 !velocity at top determined by topm
         linoutflow = .true.
         call MPI_BCAST(iinletgen, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCxT, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCxq, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCxs, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCtopm, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(linoutflow, 1, MPI_LOGICAL, 0, comm3d, mpierr)

      else if (BCxm .eq. 5) then

         write (*, *) "inoutflow conditions and idriver, setting appropriate switches (0)"

         iinletgen = 0
         idriver = 2
         BCxT = 3 !temperature is considered in inletgen & iolet
         BCxq = 3 !humidity is considered in iolet
         BCxs = 3 !scalars are considered in iolet
         BCtopm = 3 !velocity at top determined by topm
         linoutflow = .true.
         call MPI_BCAST(iinletgen, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(idriver, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCxT, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCxq, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCxs, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(BCtopm, 1, MPI_INTEGER, 0, comm3d, mpierr)
         call MPI_BCAST(linoutflow, 1, MPI_LOGICAL, 0, comm3d, mpierr)

      end if

      ! check the Poisson solver setting w.r.t. x-grid
      d(1:imax-1) = dxf(ib+1:ie) - dxf(ib:ie-1)
      inequi = any(abs(d)>dxf(ib)*1e-5)

      if ((.not. inequi) .and. (ipoiss == POISS_CYC) .and. (.not. linoutflow)) then
         write(*, *) "WARNING: consider using FFT poisson solver for better performance!"
      end if

      if ((ipoiss == POISS_FFT) .and. (inequi)) then
         write(*, *) "ERROR: POISS_FFT requires equidistant grid. Aborting..."
         call MPI_FINALIZE(mpierr)
         stop 1
      end if

   end subroutine checkinitvalues

   subroutine readinitfiles
      use modfields, only:u0, v0, w0, um, vm, wm, thlm, thl0, thl0h, qtm, qt0, qt0h, &
         ql0, ql0h, thv0h, sv0, svm, e12m, e120, &
         dudxls, dudyls, dvdxls, dvdyls, dthldxls, dthldyls, &
         dqtdxls, dqtdyls, dqtdtls, dpdx, dpdxl, dpdyl, &
         wfls, whls, ug, vg, pgx, pgy, uprof, vprof, thlprof, qtprof, e12prof, svprof, &
         v0av, u0av, qt0av, ql0av, thl0av, qt0av, sv0av, exnf, exnh, presf, presh, rhof, &
         thlpcar, uav, thvh, thvf, IIc, IIcs, IIu, IIus, IIv, IIvs, IIw, IIws
            use modglobal,         only : ib,ie,ih,jb,je,jh,kb,ke,kh,khc,kmax,dtmax,dt,runtime,timeleft,timee,ntimee,ntrun,btime,dt_lim,nsv,&
         zf, zh, dzf, dzh, rv, rd, grav, cp, rlv, pref0, om23_gs, jgb, jge, Uinf, Vinf, dy, &
         rslabs, e12min, dzh, dtheta, dqt, dsv, cexpnr, ifinput, lwarmstart, lstratstart, trestart, numol, &
         ladaptive, tnextrestart, jmax, linoutflow, lper2inout, iinletgen, lreadminl, &
         uflowrate, vflowrate,ltempeq, prandtlmoli, freestreamav, &
         tnextfielddump, tfielddump, tsample, tstatsdump, startfile, lprofforc, lchem, k1, JNO2,&
         idriver,dtdriver,driverstore,tdriverstart,tdriverdump        
      use modsubgriddata, only:ekm, ekh
      use modsurfdata, only:wtsurf, wqsurf, wsvsurf, &
         thls, thvs, ps, qts, svs, sv_top
      ! use modsurface,        only : surface,dthldz
      use modboundary, only:boundary, tqaver
      use modmpi, only:slabsum, myid, comm3d, mpierr, my_real, avexy_ibm
      use modthermodynamics, only:thermodynamics, calc_halflev
      use modinletdata, only:Uinl, Urec, Wrec, u0inletbc, v0inletbc, w0inletbc, ubulk, irecy, Utav, Ttav, &
         uminletbc, vminletbc, wminletbc, u0inletbcold, v0inletbcold, w0inletbcold, &
         storeu0inletbc, storev0inletbc, storew0inletbc, nstepread, nfile, Tinl, &
         Trec, tminletbc, t0inletbcold, t0inletbc, storet0inletbc, utaui, ttaui, iangle,&
         u0driver,v0driver,w0driver,e120driver,tdriver,thl0driver,qt0driver,storetdriver,&
         storeu0driver,storev0driver,storew0driver,storee120driver,storethl0driver,storeqt0driver,&
         nstepreaddriver
      use modinlet, only:readinletfile
      use moddriver, only: readdriverfile,initdriver,drivergen

      integer i, j, k, n

      real, allocatable :: height(:), th0av(:)
      real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh) :: thv0
      real, dimension(kb:ke) :: uaverage ! volume averaged u-velocity
      real, dimension(kb:ke) :: vaverage ! volume averaged v-velocity
      real, dimension(kb:ke) :: uaverager ! recycle plane
      real, dimension(kb:ke) :: uaveragei ! inlet plane
      real, dimension(kb:ke) :: taverager ! recycle plane
      real, dimension(kb:ke) :: taveragei ! inlet plane
      real, dimension(kb:ke + 1) :: waverage
      real, dimension(kb:ke + 1) :: uprofrot
      real, dimension(kb:ke + 1) :: vprofrot
      real tv, ran, ran1, vbulk

      character(80) chmess

      allocate (height(kb:ke + kh))
      allocate (th0av(kb:ke + kh))

      if (lstratstart) then ! Switch

         ! Read restart files as in lwarmstart
         call readrestartfiles
         um = u0
         vm = v0
         wm = w0
         thlm = thl0 !do this before or just not?
         qtm = qt0
         svm = sv0
         e12m = e120

         ! Overwrite thlm, thl0, qtm, qt0 from prof.inp.xxx
         if (myid == 0) then
            open (ifinput, file='prof.inp.'//cexpnr)
            read (ifinput, '(a80)') chmess
            write (*, '(a80)') chmess
            read (ifinput, '(a80)') chmess

            do k = kb, ke
               read (ifinput, *) &
                  height(k), &
                  thlprof(k), &
                  qtprof(k), &
                  uprof(k), &
                  vprof(k), &
                  e12prof(k)
            end do
            close (ifinput)

            write (*, *) 'height    thl     qt      u      v     e12'
            do k = ke, kb, -1
               write (*, '(f7.1,2f8.1,3f7.1)') &
                  height(k), &
                  thlprof(k), &
                  qtprof(k), &
                  uprof(k), &
                  vprof(k), &
                  e12prof(k)

            end do
         end if !myid=0

         ! MPI broadcast thl and qt
         call MPI_BCAST(thlprof, kmax, MY_REAL, 0, comm3d, mpierr)
         call MPI_BCAST(qtprof, kmax, MY_REAL, 0, comm3d, mpierr)

         do k = kb, ke
            do j = jb - 1, je + 1
               do i = ib - 1, ie + 1
                  thl0(i, j, k) = thlprof(k)
                  thlm(i, j, k) = thlprof(k)
                  qt0(i, j, k) = qtprof(k)
                  qtm(i, j, k) = qtprof(k)
               end do
            end do
         end do

         !ILS13 reintroduced thv !tg3315 this part may wrong, could need to use
         call calc_halflev
         ! exnf = (presf/pref0)**(rd/cp)  !exner functions not in restart files
         ! anymore.. or at least not read
         ! exnh = (presh/pref0)**(rd/cp)

         !   write(*,*) "exnf",enf
         !   write(*,*) "exnh",exnh
         do k = kb, ke + kh
            do j = jb, je
               do i = ib, ie
                  !write(*,*) "thl0h",thl0h(i,j,k)
                  thv0h(i, j, k) = (thl0h(i, j, k) + rlv*ql0h(i, j, k)/(cp)) &
                                   *(1 + (rv/rd - 1)*qt0h(i, j, k) - rv/rd*ql0h(i, j, k))
               end do
            end do
         end do

         do j = j, je
            do i = ib, ie
               do k = kb, ke + kh
                  thv0(i, j, k) = (thl0(i, j, k) + rlv*ql0(i, j, k)/(cp)) &
                                  *(1 + (rv/rd - 1)*qt0(i, j, k) - rv/rd*ql0(i, j, k))
               end do
            end do
         end do

         thvh = 0.
         ! call slabsum(thvh,kb,ke,thv0h,ib-ih,ie+ih,jb-jh,je+jh,kb-kh,ke+kh,ib,ie,jb,je,kb,ke) ! redefine halflevel thv using calculated thv
         call avexy_ibm(thvh(kb:ke+kh),thv0h(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
         ! thvh = thvh/rslabs

         thvf = 0.0
         call avexy_ibm(thvf(kb:ke+kh),thv0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
         ! call slabsum(thvf,kb,ke,thv0,ib-ih,ie+ih,jb-jh,je+jh,kb-kh,ke+kh,ib,ie,jb,je,kb,ke)
         ! thvf = thvf/rslabs

      else !if not lstratstart

         if (.not. lwarmstart) then

            !********************************************************************

            !    1.0 prepare initial fields from files 'prof.inp' and 'scalar.inp'
            !    ----------------------------------------------------------------

            !--------------------------------------------------------------------
            !    1.1 read fields
            !-----------------------------------------------------------------

            dt = dtmax/100.
            timee = 0.
            if (myid == 0) then
               open (ifinput, file='prof.inp.'//cexpnr)
               read (ifinput, '(a80)') chmess
               write (*, '(a80)') chmess
               read (ifinput, '(a80)') chmess

               do k = kb, ke
                  read (ifinput, *) &
                     height(k), &
                     thlprof(k), &
                     qtprof(k), &
                     uprof(k), &
                     vprof(k), &
                     e12prof(k)
               end do

               ! Apply rotation in horizontal
               !write (6, *) 'iangle = ', iangle

               !uprofrot = uprof*cos(iangle) - vprof*sin(iangle)
               !vprofrot = vprof*cos(iangle) + uprof*sin(iangle)
               !uprof = uprofrot
               !vprof = vprofrot

               close (ifinput)
               write (*, *) 'height    thl     qt      u      v     e12'
               do k = ke, kb, -1
                  write (*, '(f7.1,2f8.1,3f7.1)') &
                     height(k), &
                     thlprof(k), &
                     qtprof(k), &
                     uprof(k), &
                     vprof(k), &
                     e12prof(k)

               end do

               if (minval(e12prof(kb:ke)) < e12min) then
                  write (*, *) 'e12 value is zero (or less) in prof.inp'
                  do k = kb, ke
                     e12prof(k) = max(e12prof(k), e12min)
                  end do
               end if

            end if ! end if myid==0
            ! MPI broadcast numbers reading
            call MPI_BCAST(thlprof, kmax, MY_REAL, 0, comm3d, mpierr)
            call MPI_BCAST(qtprof, kmax, MY_REAL, 0, comm3d, mpierr)
            call MPI_BCAST(uprof, kmax, MY_REAL, 0, comm3d, mpierr)
            call MPI_BCAST(vprof, kmax, MY_REAL, 0, comm3d, mpierr)
            call MPI_BCAST(e12prof, kmax, MY_REAL, 0, comm3d, mpierr)
            do k = kb, ke
            do j = jb - 1, je + 1
            do i = ib - 1, ie + 1
               thl0(i, j, k) = thlprof(k)
               thlm(i, j, k) = thlprof(k)
               qt0(i, j, k) = qtprof(k)
               qtm(i, j, k) = qtprof(k)
               u0(i, j, k) = uprof(k)
               um(i, j, k) = uprof(k)
               v0(i, j, k) = vprof(k)
               vm(i, j, k) = vprof(k)
               w0(i, j, k) = 0.0
               wm(i, j, k) = 0.0
               e120(i, j, k) = e12prof(k)
               e12m(i, j, k) = e12prof(k)
               !        ekm (i,j,k) = 0.0
               !        ekh (i,j,k) = 0.0
               ekm(i, j, k) = numol
               ekh(i, j, k) = numol
            end do
            end do
            end do

            ekh(:, :, ke + 1) = ekh(:, :, ke) ! also for start up

            ! ILS13 30.11.17, added, not sure if necessary
            ! ILS13 30.11.1, commented
            do j = jb - jh, je + jh
               do i = ib - ih, ie + ih
                  thl0(i, j, ke + 1) = thl0(i, j, ke)
                  thl0(i, j, kb - 1) = thl0(i, j, kb)
               end do
            end do

            !! add random fluctuations
            krand = min(krand, ke)
            do k = kb, krand
               call randomnize(um, k, randu, irandom, ih, jh)
            end do
            do k = kb, krand
               call randomnize(vm, k, randu, irandom, ih, jh)
            end do
            do k = kb, krand
               call randomnize(wm, k, randu, irandom, ih, jh)
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

            u0 = um
            v0 = vm
            w0 = wm

            uaverage = 0.
            ! call slabsum(uaverage, kb, ke, um, ib - 1, ie + 1, jb - 1, je + 1, kb - 1, ke + 1, ib, ie, jb, je, kb, ke)
            do k = kb, ke
               uaverage(k) = uprof(k)*dzf(k)
            end do
            ubulk = sum(uaverage(kb:ke))/(zh(ke + 1) - zh(kb)) ! averaged u-velocity inflow profile

            vaverage = 0.
            ! call slabsum(vaverage, kb, ke, vm, ib - 1, ie + 1, jb - 1, je + 1, kb - 1, ke + 1, ib, ie, jb, je, kb, ke)
            do k = kb, ke
               vaverage(k) = vprof(k)*dzf(k)
            end do
            vbulk = sum(vaverage(kb:ke))/(zh(ke + 1) - zh(kb)) ! averaged u-velocity inflow profile

            ! Set average inlet profile to initial inlet profile in case of inletgenerator mode
            if (iinletgen == 1) then

               uaverage = 0.
               call slabsum(uaverage, kb, ke, um, ib - 1, ie + 1, jb - 1, je + 1, kb - 1, ke + 1, ib, ie, jb, je, kb, ke)
               do k = kb, ke
                  uaverage(k) = uprof(k)*dzf(k)
               end do
               ubulk = sum(uaverage(kb:ke))/(zh(ke + 1) - zh(kb)) ! volume-averaged u-velocity
               write (6, *) 'Modstartup: ubulk=', ubulk
               Utav(ib:ie, kb:ke) = um(ib:ie, jb, kb:ke)
               Uinl = um(ib, jb, kb:ke) ! set the initial time-averaged inlet profile equal to um
               Urec = um(ib, jb, kb:ke) ! set the initial time-averaged inlet profile equal to um
               Wrec(kb:ke + 1) = wm(ib, jb, kb:ke + 1) ! set the initial time-averaged inlet profile equal to mean w-profile
               u0inletbcold(jb:je, kb:ke) = um(ib, jb:je, kb:ke)
               v0inletbcold(jb:je, kb:ke) = vm(ib - 1, jb:je, kb:ke)
               w0inletbcold(jb:je, kb:ke + 1) = wm(ib - 1, jb:je, kb:ke + 1)
               uminletbc(jb:je, kb:ke) = um(ib, jb:je, kb:ke)
               vminletbc(jb:je, kb:ke) = vm(ib - 1, jb:je, kb:ke)
               wminletbc(jb:je, kb:ke) = wm(ib - 1, jb:je, kb:ke)
               u0inletbc(jb:je, kb:ke) = um(ib, jb:je, kb:ke)
               v0inletbc(jb:je, kb:ke) = vm(ib - 1, jb:je, kb:ke)
               w0inletbc(jb:je, kb:ke + 1) = wm(ib - 1, jb:je, kb:ke + 1)
               utaui = sqrt(abs(2*numol*Uinl(kb)/dzf(kb))) ! average streamwise friction at inlet (need for first time step)

               if (ltempeq) then
                  Ttav(ib:ie, kb:ke) = thlm(ib:ie, jb, kb:ke) ! set the initial time-averaged inlet profile equal to thlm
                  Tinl = thlm(ib, jb, kb:ke) ! set the initial time-averaged inlet profile equal to thlm
                  Trec = thlm(ib, jb, kb:ke) ! set the initial time-averaged inlet profile equal to thlm
                  t0inletbcold(jb:je, kb:ke) = thlm(ib - 1, jb:je, kb:ke)
                  t0inletbc(jb:je, kb:ke) = thl0(ib - 1, jb:je, kb:ke)
                  tminletbc(jb:je, kb:ke) = thlm(ib - 1, jb:je, kb:ke)
                  ttaui = numol*prandtlmoli*2.*(Tinl(kb) - thls)/(dzf(kb)*utaui) ! average friction temp. at inlet (need for first time step)
               end if

               ! add random perturbations
               if (myid == 0) then
                  call random_number(ran)
                  ran1 = -1.+2.*ran
                  write (6, *) 'random=', ran, ran1
                  call random_number(ran)
                  ran1 = -1.+2.*ran
                  write (6, *) 'random=', ran, ran1
                  call random_number(ran)
                  ran1 = -1.+2.*ran
                  write (6, *) 'random=', ran, ran1
                  call random_number(ran)
                  ran1 = -1.+2.*ran
                  write (6, *) 'random=', ran, ran1
                  call random_number(ran)
                  ran1 = -1.+2.*ran
                  write (6, *) 'random=', ran, ran1
                  call random_number(ran)
                  ran1 = -1.+2.*ran
                  write (6, *) 'random=', ran, ran1
               end if

               do k = kb + 1, kb + 48
               do j = jb, je
               do i = ib + 1, ie - 1
                  call random_number(ran)
                  ran1 = -1.+2.*ran
                  wm(i, j, k) = wm(i, j, k) + 0.1*Uinf*ran1
               end do
               end do
               end do

               !       do k=kb+1,ke-1
               do k = kb + 1, kb + 48
               do j = jb, je
               do i = ib + 2, ie - 1
                  call random_number(ran)
                  ran1 = -1.+2.*ran
                  um(i, j, k) = um(i, j, k) + 0.1*Uinf*ran1
               end do
               end do
               end do

               !       do k=kb+1,ke-1
               do k = kb + 1, kb + 48
               do j = jb, je
               do i = ib + 1, ie - 1
                  call random_number(ran)
                  ran1 = -1.+2.*ran
                  vm(i, j, k) = vm(i, j, k) + 0.1*Uinf*ran1
               end do
               end do
               end do

               u0 = um
               v0 = vm
               w0 = wm

            else if (iinletgen == 2) then

               nfile = nfile + 1
               call readinletfile
               u0inletbc(:, :) = storeu0inletbc(:, :, nstepread)
               v0inletbc(:, :) = storev0inletbc(:, :, nstepread)
               w0inletbc(:, :) = storew0inletbc(:, :, nstepread)
               uminletbc(:, :) = storeu0inletbc(:, :, nstepread)
               vminletbc(:, :) = storev0inletbc(:, :, nstepread)
               wminletbc(:, :) = storew0inletbc(:, :, nstepread)
               ! determine bulk velocity
               call slabsum(uaverage, kb, ke, u0, ib - 1, ie + 1, jb - 1, je + 1, kb - 1, ke + 1, ib, ie, jb, je, kb, ke)
               uaverage = uaverage/((ie - ib + 1)*(jge - jgb + 1)) ! this gives the i-j-averaged velocity (only correct for equidistant grid?)
               do k = kb, ke
                  uaverage(k) = uaverage(k)*dzf(k)
               end do
               ubulk = sum(uaverage(kb:ke))/(zh(ke + 1) - zh(kb)) ! volume-averaged u-velocity
               write (6, *) 'Modstartup: ubulk=', ubulk
            elseif (idriver==2) then ! idriver

               call readdriverfile

               ! if(myid==0) then
                 ! write(*,*) 'Driver inlet velocity'
                 ! do n=1,driverstore
                   ! write (*,'(f9.2,e20.12)') storetdriver(n),     storeu0driver(1,32,n)
                 ! end do
               ! endif

              ! call slabsum(uaverage,kb,ke,u0,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,ib,ie,jb,je,kb,ke)
              ! uaverage = uaverage / ((ie-ib+1)*(jge-jgb+1))  ! this gives the i-j-averaged velocity (only correct for equidistant grid?)

              call avexy_ibm(uaverage(kb:ke),u0(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,ih,jh,kh,IIu(ib:ie,jb:je,kb:ke),IIus(kb:ke),.false.)
              do k=kb,ke
                uaverage(k) = uaverage(k)*dzf(k)
              end do
              ubulk = sum(uaverage(kb:ke))/(zh(ke+1)-zh(kb)) !volume-averaged u-velocity

              if (myid==0) then
                 write(6,*) 'Modstartup: ubulk=',ubulk
              end if

            elseif (idriver==1) then

              call drivergen

            end if

            !---------------------------------------------------------------
            !  1.2 randomnize fields
            !---------------------------------------------------------------
            !     if (iinletgen /= 2 .and. iinletgen /= 1) then
            !       write(6,*) 'randomnizing temperature!'
            !       krand  = min(krand,ke)
            !        do k = kb,ke !edited tg3315 krand --> ke
            !          call randomnize(thlm,k,randthl,irandom,ih,jh)
            !          call randomnize(thl0,k,randthl,irandom,ih,jh)
            !        end do
            !       end if

            svprof = 0.
            if (myid == 0) then
               if (nsv > 0) then
                  open (ifinput, file='scalar.inp.'//cexpnr)
                  read (ifinput, '(a80)') chmess
                  read (ifinput, '(a80)') chmess
                  do k = kb, ke
                     read (ifinput, *) &
                        height(k), &
                        (svprof(k, n), n=1, nsv)
                  end do
                  open (ifinput, file='scalar.inp.'//cexpnr)
                  write (6, *) 'height   sv(1) --------- sv(nsv) '
                  do k = ke, kb, -1
                     write (6, *) &
                        height(k), &
                        (svprof(k, n), n=1, nsv)
                  end do

               end if
            end if ! end if myid==0

            call MPI_BCAST(svprof, (ke + kh - (kb - kh))*nsv, MY_REAL, 0, comm3d, mpierr)
            do k = kb, ke
               do j = jb - 1, je + 1
                  do i = ib - 1, ie + 1
                     do n = 1, nsv
                        sv0(i, j, k, n) = svprof(k, n)
                        svm(i, j, k, n) = svprof(k, n)
                     end do
                  end do
               end do
            end do
          
            if (nsv>0) then !tg3315 set these variables here for now and repeat for warmstart

              allocate(sv_top(1:nsv))
              sv_top(:) = svprof(ke,1:nsv)

              call MPI_BCAST(sv_top, nsv, MY_REAL, 0, comm3d, mpierr)

              write(*,*) 'svprof', svprof
              write(*,*) 'sv_top', sv_top

            end if
 
            !do n = 1,nsv
            !  do j = jb - jhc, je + jhc
            !    do i = ib - ihc, ie + ihc
            !      svm(i, j, ke + 1, n) = svm(i, j, ke)
            !      sv0(i, j, kb - 1, n) = sv0(i, j, kb)
            !    end do
            !  end do
            !end do

            !-----------------------------------------------------------------
            !    2.2 Initialize surface layer
            !-----------------------------------------------------------------

         !ILS13 reintroduced thv !tg3315 this part may wrong, could need to use
         call calc_halflev
         ! exnf = (presf/pref0)**(rd/cp)  !exner functions not in restart files
         ! anymore.. or at least not read
         ! exnh = (presh/pref0)**(rd/cp)

            call boundary ! tg3315 17.10.17 having this in startup was causing issues for running with lmoist ! turned of when pot. temp = temp.
            call thermodynamics ! turned off when pot. temp = temp.

            call boundary
            call thermodynamics ! turned off when pot. temp = temp.

         else !if lwarmstart
            write (*, *) "doing warmstart"
            call readrestartfiles

            um = u0
            vm = v0
            wm = w0
            thlm = thl0
            qtm = qt0
            svm = sv0
            e12m = e120
            ekm(:, :, :) = numol
            ekh(:, :, :) = numol*prandtlmoli !tg3315 added because wttop using ekh in modboundary which is called in startup

            ekh(:, :, ke + 1) = ekh(:, :, ke) ! also for start up

            if (idriver==1) then                                                                   
              !driverstore = (timeleft - tdriverstart)/dtdriver + 1
              !if(myid==0) then
              !  write(*,*) 'driverstore: ', driverstore
              !end if
              call drivergen
              tdriverdump = tdriverstart
            endif

            !ILS13 reintroduced thv
            call calc_halflev
            ! exnf = (presf/pref0)**(rd/cp)  !exner functions not in restart files
            ! anymore.. or at least not read
            ! exnh = (presh/pref0)**(rd/cp)

            !   write(*,*) "exnf",enf
            !   write(*,*) "exnh",exnh

            do j = jb, je
            do i = ib, ie
            do k = kb, ke + kh
               !write(*,*) "thl0h",thl0h(i,j,k)
               thv0h(i, j, k) = (thl0h(i, j, k) + rlv*ql0h(i, j, k)/(cp)) &
                                *(1 + (rv/rd - 1)*qt0h(i, j, k) - rv/rd*ql0h(i, j, k))
            end do
            end do
            end do

            do j = j, je
            do i = ib, ie
            do k = kb, ke + kh
               thv0(i, j, k) = (thl0(i, j, k) + rlv*ql0(i, j, k)/(cp)) &
                               *(1 + (rv/rd - 1)*qt0(i, j, k) - rv/rd*ql0(i, j, k))
            end do
            end do
            end do

            thvh = 0.
            ! call slabsum(thvh,kb,ke,thv0h,ib-ih,ie+ih,jb-jh,je+jh,kb-kh,ke+kh,ib,ie,jb,je,kb,ke) ! redefine halflevel thv using calculated thv
            call avexy_ibm(thvh(kb:ke+kh),thv0h(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
            ! thvh = thvh/rslabs

            thvf = 0.0
            call avexy_ibm(thvf(kb:ke+kh),thv0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
            ! call slabsum(thvf,kb,ke,thv0,ib-ih,ie+ih,jb-jh,je+jh,kb-kh,ke+kh,ib,ie,jb,je,kb,ke)
            ! thvf = thvf/rslabs

            ! Set average inlet profile to initial inlet profile in case of inletgenerator mode
            uaverage = 0.
            uaveragei = 0.
            uaverager = 0.
            waverage = 0.
            taveragei = 0.
            taverager = 0.
            if (iinletgen == 1) then
               call slabsum(uaveragei, kb, ke, u0, ib - 1, ie + 1, jb - 1, je + 1, kb - 1, ke + 1, ib, ib, jb, je, kb, ke)
               call slabsum(uaverager, kb, ke, u0, ib - 1, ie + 1, jb - 1, je + 1, kb - 1, ke + 1, irecy, irecy, jb, je, kb, ke)
               call slabsum(waverage, kb, ke + 1, w0, ib - 1, ie + 1, jb - 1, je + 1, kb - 1, ke + 1, ib, ie, jb, je, kb, ke + 1)
               call slabsum(uaverage, kb, ke, u0, ib - 1, ie + 1, jb - 1, je + 1, kb - 1, ke + 1, ib, ie, jb, je, kb, ke)
               uaverage = uaverage/((ie - ib + 1)*(jge - jgb + 1)) ! this gives the i-j-averaged velocity (only correct for equidistant grid?)
               uaveragei = uaveragei/(jge - jgb + 1) ! this gives the j-averaged u-velocity at the inlet
               uaverager = uaverager/(jge - jgb + 1) ! this gives the j-averaged u-velocity at the recycle plane
               waverage = waverage/((ie - ib + 1)*(jge - jgb + 1)) ! this gives the i-j-averaged w-velocity (only correct for equidistant grid?)
               if (ltempeq) then
                  call slabsum(taveragei, kb, ke, thl0, ib - 1, ie + 1, jb - 1, je + 1, kb - 1, ke + 1, ib, ie, jb, je, kb, ke)
                  call slabsum(taverager, kb, ke, thl0, ib - 1, ie + 1, jb - 1, je + 1, kb - 1, ke + 1, irecy - 1, irecy - 1, jb, je, kb, ke)
                  taveragei = taveragei/((ie - ib + 1)*(jge - jgb + 1)) ! this gives the j-averaged temperature at the inlet
                  taverager = taverager/(jge - jgb + 1) ! this gives the j-averaged temperature at the recycle plane
               end if
               if (.not. lreadminl) then
                  if (myid == 0) then
                     write (6, *) 'uaverage(kb)=', uaverage(kb)
                     write (6, *) 'uaverage(ke)=', uaverage(ke)
                     write (6, *) 'waverage(ke)=', waverage(ke)
                     write (6, *) 'waverage(ke-20)=', waverage(ke - 20)
                     write (6, *) 'taveragei(kb)=', taveragei(kb)
                     write (6, *) 'taveragei(ke)=', taveragei(ke)
                  end if

                  Utav = 0.
                  do i = ib, ie
                     Utav(i, :) = uaverage
                  end do

                  Uinl = uaverage ! set the initial time-averaged inlet profile equal to mean u-profile read from means
                  write (6, *) 'Uinl(kb+10)=', Uinl(kb + 10)
                  utaui = sqrt(abs(2*numol*Uinl(kb)/dzf(kb))) ! average streamwise friction at inlet (need for first time step)
                  Urec = uaverage ! set the initial time-averaged inlet profile equal to mean u-profile

                  Wrec(kb:ke + 1) = waverage(kb:ke + 1) ! set the initial time-averaged inlet profile equal to mean w-profile
                  Wrec(kb) = 0. ! set the initial time-averaged inlet profile equal to zero
                  if (ltempeq) then
                     Ttav = 0.
                     do i = ib, ie
                        Ttav(i, :) = taveragei(:)
                     end do
                     Tinl = taveragei
                     Trec = taveragei
                     ttaui = numol*prandtlmoli*2.*(Tinl(kb) - thls)/(dzf(kb)*utaui) ! friction temp. at inlet (need at first time step)
                  end if
               else ! -> lreadminl -> Uinl, Urec, Wrec already read
                  call slabsum(uaverage, kb, ke, u0, ib - 1, ie + 1, jb - 1, je + 1, kb - 1, ke + 1, ib, ie, jb, je, kb, ke)
                  uaverage = uaverage/((ie - ib + 1)*(jge - jgb + 1)) ! this gives the i-j-averaged velocity (only correct for equidistant grid?)
               end if

               ! determine bulk velocity
               do k = kb, ke
                  uaverage(k) = uaverage(k)*dzf(k)
               end do
               ubulk = sum(uaverage(kb:ke))/(zh(ke + 1) - zh(kb)) ! volume-averaged u-velocity
               write (6, *) 'Modstartup: ubulk=', ubulk

               do k = kb, ke
               do j = jb, je
                  uminletbc(j, k) = um(ib, j, k)
                  vminletbc(j, k) = vm(ib - 1, j, k)
                  u0inletbcold(j, k) = um(ib, j, k)
                  v0inletbcold(j, k) = vm(ib - 1, j, k)
                  u0inletbc(j, k) = um(ib, j, k)
                  v0inletbc(j, k) = vm(ib - 1, j, k)
               end do
               end do

               do k = kb, ke + 1
               do j = jb, je
                  wminletbc(j, k) = wm(ib - 1, j, k)
                  w0inletbcold(j, k) = wm(ib - 1, j, k)
                  w0inletbc(j, k) = wm(ib - 1, j, k)
               end do
               end do

               if (ltempeq) then
                  do k = kb, ke
                  do j = jb, je
                     tminletbc(j, k) = thlm(ib - 1, j, k)
                     t0inletbcold(j, k) = thlm(ib - 1, j, k)
                     t0inletbc(j, k) = thlm(ib - 1, j, k)
                  end do
                  end do
               end if

               write (6, *) 'uminletbc(jb,kb),um(ib,jb,kb)=', uminletbc(jb, kb), um(ib, jb, kb)
               write (6, *) 'uminletbc(jb+1,kb+10),um(ib,jb+1,kb+10)=', uminletbc(jb + 1, kb + 10), um(ib, jb + 1, kb + 10)
               write (6, *) 'uminletbc(je,kb+10),um(ib,je,kb+10)=', uminletbc(je, kb + 10), um(ib, je, kb + 10)

            else if (iinletgen == 2) then

               nfile = nfile + 1
               write (6, *) 'Loading inletfile'
               call readinletfile
               u0inletbc(:, :) = storeu0inletbc(:, :, nstepread)
               v0inletbc(:, :) = storev0inletbc(:, :, nstepread)
               w0inletbc(:, :) = storew0inletbc(:, :, nstepread)
               uminletbc(:, :) = storeu0inletbc(:, :, nstepread)
               vminletbc(:, :) = storev0inletbc(:, :, nstepread)
               wminletbc(:, :) = storew0inletbc(:, :, nstepread)
               if (ltempeq) then
                  t0inletbc(:, :) = storet0inletbc(:, :, nstepread)
                  tminletbc(:, :) = storet0inletbc(:, :, nstepread)
               end if
               ! determine bulk velocity
               call slabsum(uaverage, kb, ke, u0, ib - 1, ie + 1, jb - 1, je + 1, kb - 1, ke + 1, ib, ie, jb, je, kb, ke)
               uaverage = uaverage/((ie - ib + 1)*(jge - jgb + 1)) ! this gives the i-j-averaged velocity (only correct for equidistant grid?)
               do k = kb, ke
                  uaverage(k) = uaverage(k)*dzf(k)
               end do
               ubulk = sum(uaverage(kb:ke))/(zh(ke + 1) - zh(kb)) ! volume-averaged u-velocity
               write (6, *) 'Modstartup: ubulk=', ubulk

            elseif (idriver==2) then ! idriver

               call readdriverfile
               call drivergen

              !call slabsum(uaverage,kb,ke,u0,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,ib,ie,jb,je,kb,ke)
              !uaverage = uaverage / ((ie-ib+1)*(jge-jgb+1))  ! this gives the i-j-averaged velocity (only correct for equidistant grid?)
              call avexy_ibm(uaverage(kb:ke),u0(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,ih,jh,kh,IIu(ib:ie,jb:je,kb:ke),IIus(kb:ke),.false.)
              do k=kb,ke
                uaverage(k) = uaverage(k)*dzf(k)
              end do
              ubulk = sum(uaverage(kb:ke))/(zh(ke+1)-zh(kb)) !volume-averaged u-velocity
              if (myid==0) then
                 write(6,*) 'Modstartup: ubulk=',ubulk
              end if

            end if ! iinletgen/idriver

            if (lper2inout) then ! if the restart starts from a periodic simulation to in/outflow, lper2inout should be set to .true.
               if (myid == 0) then
                  write (6, *) 'per2inout=.true. -> reading inlet profile from prof.inp.XXX and scalar.inp.XXX'
                  open (ifinput, file='prof.inp.'//cexpnr) !  read the inlet profile from prof.inp
                  read (ifinput, '(a80)') chmess
                  write (*, '(a80)') chmess
                  read (ifinput, '(a80)') chmess

                  do k = kb, ke
                     read (ifinput, *) &
                        height(k), &
                        thlprof(k), &
                        qtprof(k), &
                        uprof(k), &
                        vprof(k), &
                        e12prof(k)
                  end do
                  svprof = 0.
                  if (nsv > 0) then
                     open (ifinput, file='scalar.inp.'//cexpnr)
                     read (ifinput, '(a80)') chmess
                     read (ifinput, '(a80)') chmess
                     do k = kb, ke
                        read (ifinput, *) &
                           height(k), &
                           (svprof(k, n), n=1, nsv)
                     end do
                     open (ifinput, file='scalar.inp.'//cexpnr)
                     write (6, *) 'height   sv(1) --------- sv(nsv) '
                     do k = ke, kb, -1
                        write (6, *) &
                           height(k), &
                           (svprof(k, n), n=1, nsv)
                     end do

                  end if
               end if ! end if myid==0

               ! MPI broadcast numbers reading
               call MPI_BCAST(thlprof, kmax, MY_REAL, 0, comm3d, mpierr)
               call MPI_BCAST(uprof, kmax, MY_REAL, 0, comm3d, mpierr)
               call MPI_BCAST(vprof, kmax, MY_REAL, 0, comm3d, mpierr)
               call MPI_BCAST(e12prof, kmax, MY_REAL, 0, comm3d, mpierr)
               call MPI_BCAST(qtprof, kmax, MY_REAL, 0, comm3d, mpierr)
               call MPI_BCAST(svprof, (ke + kh - (kb - kh))*nsv, MY_REAL, 0, comm3d, mpierr)

            else if (linoutflow) then ! restart of inoutflow simulation: reproduce inlet boundary condition from restartfile
               do j = jb - 1, je + 1
                  do k = kb, ke + 1
                     uprof(k) = u0(ib, j, k)
                     vprof(k) = (v0(ib - 1, j, k) + v0(ib, j, k))/2
                     thlprof(k) = (thl0(ib - 1, j, k) + thl0(ib, j, k))/2
                     qtprof(k) = (qt0(ib - 1, j, k) + qt0(ib, j, k))/2
                     e12prof(k) = (e120(ib - 1, j, k) + e120(ib, j, k))/2
                     do n = 1, nsv
                        svprof(k, n) = (sv0(ib - 1, j, k, n) + sv0(ib, j, k, n))/2
                     enddo
                  enddo
               enddo
               ! outlet bulk velocity
               call slabsum(uaverage, kb, ke, u0, ib - 1, ie + 1, jb - 1, je + 1, kb - 1, ke + 1, ib, ie, jb, je, kb, ke)
               uaverage = uaverage/((ie - ib + 1)*(jge - jgb + 1)) ! this gives the i-j-averaged velocity (only correct for equidistant grid?)
               ! determine bulk velocity
               do k = kb, ke
                  uaverage(k) = uaverage(k)*dzf(k)
               end do
               ubulk = sum(uaverage(kb:ke))/(zh(ke + 1) - zh(kb)) ! volume-averaged u-velocity
               write (6, *) 'Modstartup: ubulk=', ubulk
            else ! else per2per... read svprof regardless...

            ! tg3315 read svprof (but do not use regardless of above...)
              svprof = 0.
              if (myid == 0) then
                 if (nsv > 0) then
                    open (ifinput, file='scalar.inp.'//cexpnr)
                    read (ifinput, '(a80)') chmess
                    read (ifinput, '(a80)') chmess
                    do k = kb, ke
                       read (ifinput, *) &
                          height(k), &
                          (svprof(k, n), n=1, nsv)
                    end do
                    open (ifinput, file='scalar.inp.'//cexpnr)
                    write (6, *) 'height   sv(1) --------- sv(nsv) '
                    do k = ke, kb, -1
                       write (6, *) &
                          height(k), &
                          (svprof(k, n), n=1, nsv)
                    end do

                 end if
              end if ! end if myid==0

              call MPI_BCAST(svprof, (ke + kh - (kb - kh))*nsv, MY_REAL, 0, comm3d, mpierr)
          
              if (nsv>0) then !tg3315 set these variables here for now and repeat for warmstart

                allocate(sv_top(1:nsv))
                sv_top(:) = svprof(ke,1:nsv)

                call MPI_BCAST(sv_top, nsv, MY_REAL, 0, comm3d, mpierr)

              end if

            end if ! end if lper2inout

            u0av = 0.0
            v0av = 0.0
            thl0av = 0.0
            qt0av = 0.0
            th0av = 0.0
            sv0av = 0.

            ! call slabsum(u0av  ,kb,ke+kh,u0(:,:,kb:ke+kh)  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
            call avexy_ibm(u0av(kb:ke+kh),u0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIu(ib:ie,jb:je,kb:ke+kh),IIus(kb:ke+kh),.false.)
            ! call slabsum(v0av  ,kb,ke+kh,v0(:,:,kb:ke+kh)  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
            call avexy_ibm(v0av(kb:ke+kh),v0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIv(ib:ie,jb:je,kb:ke+kh),IIvs(kb:ke+kh),.false.)
            ! call slabsum(thl0av,kb,ke+kh,thl0(:,:,kb:ke+kh),ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
            call avexy_ibm(thl0av(kb:ke+kh),thl0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
            ! call slabsum(qt0av,kb,ke+kh,qt0(:,:,kb:ke+kh),ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
            call avexy_ibm(qt0av(kb:ke+kh),qt0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
            do n = 1, nsv
               ! call slabsum(sv0av(kb,n),kb,ke+kh,sv0(ib-ih,jb-jh,kb,n),ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ib,ie,jb,je,kb,ke+kh)
               call avexy_ibm(sv0av(kb:ke+khc,n),sv0(ib:ie,jb:je,kb:ke+khc,n),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+khc),IIcs(kb:ke+khc),.false.)
            end do

            ! CvH - only do this for fixed timestepping. In adaptive dt comes from restartfile
            if (.not. ladaptive) dt = dtmax
            !  call boundary
         end if ! lwarmstart
      end if ! not lstratstart
      !-----------------------------------------------------------------
      !    2.1 read and initialise fields
      !-----------------------------------------------------------------

      if (myid == 0) then
         open (ifinput, file='lscale.inp.'//cexpnr)
         read (ifinput, '(a80)') chmess
         read (ifinput, '(a80)') chmess
         write (6, *) ' height  u_geo  v_geo  pgx  pgy  subs     ' &
            , '   dqtdx      dqtdy        dqtdtls     thl_rad '
         do k = kb, ke
            read (ifinput, *) &
               height(k), &
               ug(k), &
               vg(k), &
               pgx(k), &
               pgy(k), &
               wfls(k), &
               dqtdxls(k), &
               dqtdyls(k), &
               dqtdtls(k), &
               thlpcar(k)
         end do
         close (ifinput)

         do k = ke, kb, -1
            write (6, '(3f7.1,5e12.4)') &
               height(k), &
               ug(k), &
               vg(k), &
               pgx(k), &
               pgy(k), &
               wfls(k), &
               dqtdxls(k), &
               dqtdyls(k), &
               dqtdtls(k), &
               thlpcar(k)
         end do

      end if ! end myid==0

      ! MPI broadcast variables read in

      call MPI_BCAST(ug, kmax, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(vg, kmax, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(pgx, kmax, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(pgy, kmax, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wfls, kmax, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(dqtdxls, kmax, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(dqtdyls, kmax, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(dqtdtls, kmax, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(thlpcar, kmax, MY_REAL, 0, comm3d, mpierr)

      !-----------------------------------------------------------------
      !    2.3 make large-scale horizontal pressure gradient
      !-----------------------------------------------------------------

      !******include rho if rho = rho(z) /= 1.0 ***********

      if (lprofforc) then !tg3315
         do k = kb, ke
            dpdxl(k) = -pgx(k) - dpdx !-om23_gs*ug(k)-pgx(k)-dpdx
            dpdyl(k) = -pgy(k)
         end do
      else
         do k = kb, ke
            !dpdxl(k) =  om23_gs*vg(k)
            !dpdyl(k) = -om23_gs*ug(k)
            !dpdxl(k) =  -ug(k)
            !dpdyl(k) =  -vg(k)
            dpdxl(k) = om23_gs*vg(k) - pgx(k) - dpdx !corriolis forcing and pressure gradient
            dpdyl(k) = -om23_gs*ug(k) - pgy(k)
         end do
      endif

      !-----------------------------------------------------------------
      !    2.4 large-scale subsidence, reintroduced ILS13 05.06.2014
      !-----------------------------------------------------------------

      whls(kb) = 0.0
      do k = kb + 1, ke
         whls(k) = (wfls(k)*dzf(k - 1) + wfls(k - 1)*dzf(k))/(2*dzh(k))
      end do
      whls(ke + 1) = (wfls(ke) + dzf(ke)*(wfls(ke) - wfls(ke - 1))/dzh(ke)) ! tg3315 31/07/18 removed a 0.5

      !    idtmax = floor(dtmax/tres)
      btime = timee
      !    timeleft=ceiling(runtime/tres)
      timeleft = runtime
      dt_lim = timeleft
      !    write(6,*) 'real(dt)*tres= ',rdt, ' dtmax/100= ',dtmax/100
      if ((lwarmstart) .or. (lstratstart)) then ! tg3315 to have cumulative number on restart files
         read (startfile(6:13), '(i4)') ntrun
         ! ntrun = ichar(startfile(6:13))
      else
         ntrun = 0
      end if
      ntimee = nint(timee/dtmax)
      tnextrestart = btime + trestart
      tnextfielddump = btime + tfielddump
      deallocate (height, th0av)

      !    call boundary

   end subroutine readinitfiles

   subroutine readrestartfiles

      use modsurfdata, only:ustar, thlflux, qtflux, svflux, dudz, dvdz, dthldz, dqtdz, ps, thls, qts, thvs, oblav, &
         wtsurf
      use modfields, only:u0, v0, w0, thl0, qt0, ql0, ql0h, qtav, qlav, e120, dthvdz, presf, presh, sv0, mindist, wall, &
         uav, vav, wav, uuav, vvav, wwav, uvav, uwav, vwav, svav, thlav, thl2av, sv2av, pres0, svm, &
         svprof, viscratioav, thluav, thlvav, thlwav, svuav, svvav, svwav, presav, &
         uusgsav, vvsgsav, wwsgsav, uwsgsav, thlusgsav, thlwsgsav, svusgsav, svwsgsav, tkesgsav, &
         strain2av, nusgsav
      use modglobal, only:ib, ie, ih, jb, je, jh, kb, ke, kh, dtheta, dqt, dsv, startfile, timee, totavtime, runavtime, &
         iexpnr, ntimee, rk3step, ifinput, nsv, runtime, dt, cexpnr, lreadmean, lreadminl, &
         totinletav, lreadscal, ltempeq, dzf, numol, prandtlmoli
      use modmpi, only:cmyid, myid
      use modsubgriddata, only:ekm
      use modinlet, only:zinterpolate1d, zinterpolatet1d, zinterpolatew1d, zinterpolate2d
      use modinletdata, only:Uinl, Urec, Wrec, Utav, Tinl, Trec, linuf, linuh, &
         kbin, kein, lzinzsim, utaui, Ttav, ttaui

      real, dimension(ib:ie, jb:je, kb:ke)  ::  dummy3d
      real, dimension(ib:ie, kbin:kein)    ::  Utavin
      real, dimension(ib:ie, kbin:kein)    ::  Ttavin
      real, dimension(kbin:kein)          ::  Uinlin
      real, dimension(kbin:kein)          ::  Urecin
      real, dimension(kbin:kein)          ::  Tinlin
      real, dimension(kbin:kein)          ::  Trecin
      real, dimension(kbin:kein + 1)        ::  Wrecin
      character(50) :: name, name2, name4
      real dummy
      integer i, j, k, n
      !********************************************************************

      !    1.0 Read initfiles
      !-----------------------------------------------------------------
      name = startfile
      name(5:5) = 'd'
      name(15:17) = cmyid
      write (6, *) 'loading ', name
      open (unit=ifinput, file=name, form='unformatted', status='old')

      read (ifinput) (((mindist(i, j, k), i=ib, ie), j=jb, je), k=kb, ke)
      read (ifinput) ((((wall(i, j, k, n), i=ib, ie), j=jb, je), k=kb, ke), n=1, 5)
      read (ifinput) (((u0(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb, ke + kh)
      read (ifinput) (((v0(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb, ke + kh)
      read (ifinput) (((w0(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb, ke + kh)
      read (ifinput) (((pres0(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb, ke + kh)
      read (ifinput) (((thl0(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb, ke + kh)
      read (ifinput) (((e120(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb, ke + kh)
      read (ifinput) (((ekm(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb, ke + kh)
      read (ifinput) (((qt0(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb, ke + kh)
      read (ifinput) (((ql0(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb, ke + kh)
      read (ifinput) (((ql0h(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb, ke + kh)
      read (ifinput) timee, dt
      close (ifinput)
      write (6, *) 'finished loading ', name

      if ((nsv > 0) .and. (lreadscal)) then
         name(5:5) = 's'
         write (6, *) 'loading ', name
         open (unit=ifinput, file=name, form='unformatted')
         read (ifinput) ((((sv0(i, j, k, n), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb, ke + kh), n=1, nsv)
         read (ifinput) timee
         close (ifinput)
         write (6, *) 'finished loading ', name
      elseif ((nsv > 0) .and. (.not. lreadscal)) then
         sv0 = 0.
         svprof = 0.
      end if

      ! read mean variables if asked for by lreadmean
      name2 = 'means   .'
      name2(6:8) = cmyid
      name2(10:12) = cexpnr
      if (lreadmean) then
         write (6, *) 'Reading meansXXX.XXX, proc = ', myid
         open (unit=ifinput, file=name2, form='unformatted')
         read (ifinput) totavtime, nsv
         read (ifinput) (((uav(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb - kh, ke + kh)
         read (ifinput) (((vav(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb - kh, ke + kh)
         read (ifinput) (((wav(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb - kh, ke + kh)
         read (ifinput) (((thlav(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb - kh, ke + kh)
         read (ifinput) (((qtav(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb - kh, ke + kh)
         read (ifinput) (((qlav(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb - kh, ke + kh)
         read (ifinput) (((presav(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb - kh, ke + kh)
         read (ifinput) ((((svav(i, j, k, n), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb - kh, ke + kh), n=1, nsv)
         read (ifinput) (((viscratioav(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb - kh, ke + kh)
         read (ifinput) (((uuav(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb - kh, ke + kh)
         read (ifinput) (((vvav(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb - kh, ke + kh)
         read (ifinput) (((wwav(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb - kh, ke + kh)
         read (ifinput) (((thl2av(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb - kh, ke + kh)
         read (ifinput) ((((sv2av(i, j, k, n), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb - kh, ke + kh), n=1, nsv)
         read (ifinput) (((uvav(i, j, k), i=ib, ie + ih), j=jb, je + jh), k=kb, ke)
         read (ifinput) (((uwav(i, j, k), i=ib, ie + ih), j=jb, je), k=kb, ke + kh)
         read (ifinput) (((vwav(i, j, k), i=ib, ie), j=jb, je + jh), k=kb, ke + kh)
         read (ifinput) (((thluav(i, j, k), i=ib, ie), j=jb, je), k=kb, ke)
         read (ifinput) (((thlvav(i, j, k), i=ib, ie), j=jb, je + jh), k=kb, ke)
         read (ifinput) (((thlwav(i, j, k), i=ib, ie), j=jb, je), k=kb, ke + kh)
         read (ifinput) ((((svuav(i, j, k, n), i=ib, ie), j=jb, je), k=kb, ke), n=1, nsv)
         read (ifinput) ((((svvav(i, j, k, n), i=ib, ie), j=jb, je + jh), k=kb, ke), n=1, nsv)
         read (ifinput) ((((svwav(i, j, k, n), i=ib, ie), j=jb, je), k=kb, ke + kh), n=1, nsv)
         close (ifinput)
         write (6, *) 'Total averaging time so far: ', totavtime

         ! read <x'y'>_SGS to file.
         name2 = 'SGS__   .'
         name2(6:8) = cmyid
         name2(10:12) = cexpnr
         open (unit=ifinput, file=name2, form='unformatted')
         read (ifinput) dummy, dummy
         read (ifinput) (((uusgsav(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb - kh, ke + kh)
         read (ifinput) (((vvsgsav(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb - kh, ke + kh)
         read (ifinput) (((wwsgsav(i, j, k), i=ib - ih, ie + ih), j=jb - jh, je + jh), k=kb - kh, ke + kh)
         read (ifinput) (((uwsgsav(i, j, k), i=ib, ie + ih), j=jb, je), k=kb, ke + kh)
         read (ifinput) (((dummy3d(i, j, k), i=ib, ie), j=jb, je), k=kb, ke) ! this is dissresav, which will be computed using other mean quantities
         read (ifinput) (((tkesgsav(i, j, k), i=ib, ie), j=jb, je), k=kb, ke)
         read (ifinput) (((dummy3d(i, j, k), i=ib, ie), j=jb, je), k=kb, ke) ! this is disssgsav, which will be computed using other mean quantities
         read (ifinput) (((strain2av(i, j, k), i=ib, ie), j=jb, je), k=kb, ke) ! <SijSij> (NOT <Sij><Sij> !!) (average over time)
         read (ifinput) (((nusgsav(i, j, k), i=ib, ie), j=jb, je), k=kb, ke) ! <nu_sgs> (average over time)
         read (ifinput) (((thlusgsav(i, j, k), i=ib, ie + ih), j=jb, je), k=kb, ke)
         read (ifinput) (((thlwsgsav(i, j, k), i=ib, ie), j=jb, je), k=kb, ke + kh)
         read (ifinput) ((((svusgsav(i, j, k, n), i=ib, ie + ih), j=jb, je), k=kb, ke), n=1, nsv)
         read (ifinput) ((((svwsgsav(i, j, k, n), i=ib, ie), j=jb, je), k=kb, ke + kh), n=1, nsv)
         close (ifinput)

      end if

      ! read mean profiles for inlet generator
      if (lreadminl) then
         if (.not. lzinzsim) then
            name4 = 'meaninlet.   '
            name4(11:13) = cexpnr
            open (unit=ifinput, file=name4, form='unformatted')
            read (ifinput) totinletav ! interval of time-average
            read (ifinput) (Uinlin(k), k=kbin, kein)
            read (ifinput) (Urecin(k), k=kbin, kein)
            read (ifinput) (Wrecin(k), k=kbin, kein + 1)
            read (ifinput) ((Utavin(i, k), i=ib, ie), k=kbin, kein)
            close (ifinput)

            call zinterpolate1d(Uinlin, Uinl) ! interpolate inlet profile to zgrid
            call zinterpolate1d(Urecin, Urec)
            call zinterpolatew1d(Wrecin, Wrec)
            call zinterpolate2d(Utavin, Utav)

            if (ltempeq) then
               name4 = 'tempinlet.   '
               name4(11:13) = cexpnr
               open (unit=ifinput, file=name4, form='unformatted')
               read (ifinput) totinletav ! interval of time-average
               read (ifinput) (Tinlin(k), k=kbin, kein)
               read (ifinput) (Trecin(k), k=kbin, kein)
               read (ifinput) ((Ttavin(i, k), i=ib, ie), k=kbin, kein)
               close (ifinput)

               call zinterpolatet1d(Tinlin, Tinl)
               call zinterpolatet1d(Trecin, Trec)
               call zinterpolate2d(Ttavin, Ttav)
            end if ! ltempeq

         else !lzinzsim=.true. -> inlet grid equals sim grid
            name4 = 'meaninlet.   '
            name4(11:13) = cexpnr
            open (unit=ifinput, file=name4, form='unformatted')
            read (ifinput) totinletav ! interval of time-average
            read (ifinput) (Uinl(k), k=kb, ke)
            read (ifinput) (Urec(k), k=kb, ke)
            read (ifinput) (Wrec(k), k=kb, ke + 1)
            read (ifinput) ((Utav(i, k), i=ib, ie), k=kb, ke)

            close (ifinput)

            if (ltempeq) then
               name4 = 'tempinlet.   '
               name4(11:13) = cexpnr
               open (unit=ifinput, file=name4, form='unformatted')
               read (ifinput) totinletav ! interval of time-average
               read (ifinput) (Tinl(k), k=kb, ke)
               read (ifinput) (Trec(k), k=kb, ke)
               read (ifinput) ((Ttav(i, k), i=ib, ie), k=kb, ke)

               close (ifinput)
            end if ! ltempeq
         end if ! lzinzsim

         utaui = sqrt(abs(2*numol*Uinl(kb)/dzf(kb))) ! average streamwise friction at inlet (need for first time step)
         if (ltempeq) then
            ttaui = numol*prandtlmoli*2.*(Tinl(kb) - thls)/(dzf(kb)*utaui)
         end if
      end if !(lreadminl)

   end subroutine readrestartfiles

   subroutine exitmodules
      use modfields, only:exitfields
      use modglobal, only:exitglobal
      use modmpi, only:exitmpi
      use modpois, only:exitpois
      use modsubgrid, only:exitsubgrid
      use modthermodynamics, only:exitthermodynamics
      use modinlet, only:exitinlet

      call exitthermodynamics
      call exitsubgrid
      call exitpois
      call exitfields
      call exitglobal
      call exitinlet
      call exitmpi

   end subroutine exitmodules

   subroutine randomnize(field, klev, ampl, ir, ihl, jhl)

      use modmpi, only:myid, nprocs
      use modglobal, only:ib, ie, imax, jmax, jb, je, kb, ke, kh
      integer(KIND=selected_int_kind(6)):: imm, ia, ic, ir
      integer ihl, jhl
      integer i, j, klev
      integer m, mfac
      real ran, ampl
      real field(ib - ihl:ie + ihl, jb - jhl:je + jhl, kb - kh:ke + kh)
      parameter(imm=134456, ia=8121, ic=28411)

      if (myid > 0) then
         mfac = myid*jmax*imax
         do m = 1, mfac
            ir = mod((ir)*ia + ic, imm)

         end do
      end if
      do j = jb, je
      do i = ib, ie
         ir = mod((ir)*ia + ic, imm)
         ran = real(ir)/real(imm)
         field(i, j, klev) = field(i, j, klev) + (ran - 0.5)*2.0*ampl
      end do
      end do

      if (nprocs - 1 - myid > 0) then
         mfac = (nprocs - 1 - myid)*imax*jmax
         do m = 1, mfac
            ir = mod((ir)*ia + ic, imm)
         end do
      end if

      return
   end subroutine randomnize

   subroutine createmasks
      use modglobal, only:ib, ie, ih, ihc, jb, je, jh, jhc, kb, ke, kh, khc, rslabs, jmax, nblocks,&
         ifinput, cexpnr, libm, jtot, block
      use modfields, only:IIc, IIu, IIv, IIw, IIuw, IIvw, IIuv, IIct, IIwt, IIut, IIuwt, IIvt,&
         IIcs, IIus, IIuws, IIvws, IIuvs, IIvs, IIws, &
         um, u0, vm, v0, wm, w0
      use modmpi, only:myid, comm3d, mpierr, MPI_INTEGER, MPI_DOUBLE_PRECISION, MY_REAL, nprocs, &
         cmyid, MPI_REAL8, MPI_REAL4, MPI_SUM, excjs
      ! use initfac, only:block
      integer k, n, il, iu, jl, ju, kl, ku
      integer :: IIcl(kb:ke + khc), IIul(kb:ke + khc), IIvl(kb:ke + khc), IIwl(kb:ke + khc), IIuwl(kb:ke + khc), IIvwl(kb:ke + khc), IIuvl(kb:ke + khc)
      integer :: IIcd(ib:ie, kb:ke)
      integer :: IIwd(ib:ie, kb:ke)
      integer :: IIuwd(ib:ie, kb:ke)
      integer :: IIud(ib:ie, kb:ke)
      integer :: IIvd(ib:ie, kb:ke)
      character(80) chmess, name2

      ! II*l needn't be defined up to ke_khc, but for now would require large scale changes in modstatsdump so if works leave as is ! tg3315 04/07/18

      if (.not. libm) then
         IIc(:, :, :) = 1
         IIu(:, :, :) = 1
         IIv(:, :, :) = 1
         IIw(:, :, :) = 1
         IIuw(:, :, :) = 1
         IIvw(:, :, :) = 1
         IIuv(:, :, :) = 1
         IIcs(:) = nint(rslabs)
         IIus(:) = nint(rslabs)
         IIvs(:) = nint(rslabs)
         IIws(:) = nint(rslabs)
         IIuws(:) = nint(rslabs)
         IIvws(:) = nint(rslabs)
         IIuvs(:) = nint(rslabs)
         IIct(:, :) = jtot 
         IIut(:, :) = jtot
         IIvt(:, :) = jtot
         IIwt(:, :) = jtot
         IIuwt(:, :) = jtot
         return
      end if

      allocate (block(1:nblocks, 1:11))

      if (myid == 0) then
         if (nblocks > 0) then
            open (ifinput, file='blocks.inp.'//cexpnr)
            read (ifinput, '(a80)') chmess
            read (ifinput, '(a80)') chmess
            do n = 1, nblocks
               read (ifinput, *) &
                  block(n, 1), &
                  block(n, 2), &
                  block(n, 3), &
                  block(n, 4), &
                  block(n, 5), &
                  block(n, 6), &
                  block(n, 7), &
                  block(n, 8), &
                  block(n, 9), &
                  block(n, 10), &
                  block(n, 11)
            end do
            close (ifinput)

            do n = 1, nblocks
               write (6, *) &
                  n, &
                  block(n, 1), &
                  block(n, 2), &
                  block(n, 3), &
                  block(n, 4), &
                  block(n, 5), &
                  block(n, 6)
            end do
         end if !nblocks>0
      end if !myid

      call MPI_BCAST(block, 11*nblocks, MPI_INTEGER, 0, comm3d, mpierr)

      ! Create masking matrices
      IIc = 1; IIu = 1; IIv = 1; IIct = 1; IIw = 1; IIuw = 1; IIvw = 1; IIuv = 1; IIwt = 1; IIut = 1; IIvt = 1; IIuwt = 1; IIcs = 1; IIus = 1; IIvs = 1; IIws = 1; IIuws = 1; IIvws = 1; IIuvs = 1

      do n = 1, nblocks
         il = block(n, 1)
         iu = block(n, 2)
         !kl = block(n, 5)
         kl = kb ! tg3315 changed as buildings for lEB must start at kb+1 not kb with no block below
         ku = block(n, 6)
         jl = block(n, 3) - myid*jmax
         ju = block(n, 4) - myid*jmax
         if (ju < jb - 1 .or. jl > je) then
            cycle
         else
            if (ju >= je) then !tg3315 04/07/18 to avoid ju+1 when is last cell...
               if (jl < jb) jl = jb
               ju = je

               ! Masking matrices !tg3315
               IIc(il:iu, jl:ju, kl:ku) = 0
               IIu(il:iu + 1, jl:ju, kl:ku) = 0
               IIv(il:iu, jl:ju, kl:ku) = 0
               IIw(il:iu, jl:ju, kl:ku + 1) = 0
               IIuw(il:iu + 1, jl:ju, kl:ku + 1) = 0
               IIvw(il:iu, jl:ju, kl:ku + 1) = 0
               IIuv(il:iu + 1, jl:ju, kl:ku) = 0

            else if (ju == jb - 1) then ! if end of block is in cell before proc

               IIv(il:iu, jb, kl:ku) = 0
               IIvw(il:iu, jb, kl:ku + 1) = 0
               IIuv(il:iu + 1, jb, kl:ku) = 0

            else ! ju is in this proc...
               if (jl < jb) jl = jb

               ! Masking matrices !tg3315
               IIc(il:iu, jl:ju, kl:ku) = 0
               IIu(il:iu + 1, jl:ju, kl:ku) = 0
               IIv(il:iu, jl:ju + 1, kl:ku) = 0
               IIw(il:iu, jl:ju, kl:ku + 1) = 0
               IIuw(il:iu + 1, jl:ju, kl:ku + 1) = 0
               IIvw(il:iu, jl:ju + 1, kl:ku + 1) = 0
               IIuv(il:iu + 1, jl:ju + 1, kl:ku) = 0

            end if

            ! ensure that ghost cells know where blocks are !tg3315 this is not necessary
            ! if (jl<jb+jh)  IIc(il:iu,je+jh,kl:ku) = 0
            ! if (jl<jb+jhc) IIc(il:iu,je+jhc,kl:ku) = 0
            ! if (ju>je-jh)  IIc(il:iu,jb-jh,kl:ku) = 0
            ! if (ju>je-jhc) IIc(il:iu,jb-jhc,kl:ku) = 0

            ! if (il<ib+ih)  IIc(ie+ih,jl:ju,kl:ku) = 0
            ! if (il<ib+ihc) IIc(ie+ihc,jl:ju,kl:ku) = 0
            ! if (iu>ie-ih)  IIc(ib-ih,jl:ju,kl:ku) = 0
            ! if (iu>ie-ihc) IIc(ib-ihc,jl:ju,kl:ku) = 0

         end if
      end do

      IIw(:, :, kb) = 0; IIuw(:, :, kb) = 0; IIvw(:, :, kb) = 0

      ! for correct ghost cells from adjacent processors !tg3315 ?unsure if this is correct
      ! tg3315 22/11/17 does not work because II is an integer and needs real numbers... !tg3315 not necessary
      !call excjs( IIc  , ib,ie,jb,je,kb,ke+khc,ihc,jhc)
      !call excjs( IIu  , ib,ie,jb,je,kb,ke+khc,ihc,jhc)
      !call excjs( IIv  , ib,ie,jb,je,kb,ke+khc,ihc,jhc)
      !call excjs( IIw  , ib,ie,jb,je,kb,ke+khc,ihc,jhc)

      do k = kb, ke + khc
         IIcl(k) = sum(IIc(ib:ie, jb:je, k))
         IIul(k) = sum(IIu(ib:ie, jb:je, k))
         IIvl(k) = sum(IIv(ib:ie, jb:je, k))
         IIwl(k) = sum(IIw(ib:ie, jb:je, k))
         IIuwl(k) = sum(IIuw(ib:ie, jb:je, k))
         IIvwl(k) = sum(IIvw(ib:ie, jb:je, k))
         IIuvl(k) = sum(IIuv(ib:ie, jb:je, k))
      enddo

      call MPI_ALLREDUCE(IIcl, IIcs, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIul, IIus, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIvl, IIvs, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIwl, IIws, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIuwl, IIuws, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIvwl, IIvws, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIuvl, IIuvs, ke + khc - kb + 1, MPI_INTEGER, &
                         MPI_SUM, comm3d, mpierr)

      IIcd(ib:ie, kb:ke) = sum(IIc(ib:ie, jb:je, kb:ke), DIM=2)
      IIwd(ib:ie, kb:ke) = sum(IIw(ib:ie, jb:je, kb:ke), DIM=2)
      IIuwd(ib:ie, kb:ke) = sum(IIuw(ib:ie, jb:je, kb:ke), DIM=2)
      IIud(ib:ie, kb:ke) = sum(IIu(ib:ie, jb:je, kb:ke), DIM=2)
      IIvd(ib:ie, kb:ke) = sum(IIv(ib:ie, jb:je, kb:ke), DIM=2)

      call MPI_ALLREDUCE(IIwd(ib:ie, kb:ke), IIwt(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIcd(ib:ie, kb:ke), IIct(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIuwd(ib:ie, kb:ke), IIuwt(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIud(ib:ie, kb:ke), IIut(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(IIvd(ib:ie, kb:ke), IIvt(ib:ie, kb:ke), (ke - kb + 1)*(ie - ib + 1), MPI_INTEGER, MPI_SUM, comm3d, mpierr)

      ! masking matrix for switch if entire slab is blocks
      !if (IIcs(kb) == 0) then
      !  IIbl = 0
      !else
      !  IIbl = 1
      !end if

      !where (IIcs == 0)
      !IIcs = nint(rslabs)
      !endwhere
      !where (IIus == 0)
      !IIus = nint(rslabs)
      !endwhere
      !where (IIvs == 0)
      !IIvs = nint(rslabs)
      !endwhere
      !where (IIws == 0)
      !IIws = nint(rslabs)
      !endwhere
      !where (IIuws == 0)
      !IIuws = nint(rslabs)
      !endwhere
      !where (IIvws == 0)
      !IIvws = nint(rslabs)
      !endwhere

      ! use masking matrices to set 0 in blocks from start? tg3315 13/12/17
      ! um(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh) = IIu(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)*um(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
      ! vm(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh) = IIv(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)*vm(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
      ! wm(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh) = IIw(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)*wm(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)

      ! u0 = um
      ! v0 = vm
      ! w0 = wm

   end subroutine createmasks

end module modstartup
