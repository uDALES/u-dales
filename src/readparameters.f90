!> \file readparameters.f90
!!  Reads and broadcasts configuration parameters for uDALES

module readparameters

   use mpi
   implicit none
   save

   public :: readnamelists
   public :: irandom, krand, randu, randthl, randqt

   integer(KIND=selected_int_kind(6)) :: irandom = 43
   integer :: krand = huge(0)
   real :: randu = 0.01, randthl = 0.0, randqt = 0.0

contains
   subroutine readnamelists

      !-----------------------------------------------------------------|
      !                                                                 |
      !     Reads all general options from namoptions                   |
      !                                                                 |
      !      Jasper Tomas                 31/03/2014                    |
      !      Chiel van Heerwaarden        15/06/2007                    |
      !      Thijs Heus                   15/06/2007                    |
      !-----------------------------------------------------------------|

      use modglobal,         only : initglobal, iexpnr, runtime, dtmax,  &
                                    lwarmstart, lstratstart, lfielddump, lreadscal, startfile, tfielddump, fieldvars, slicevars, tsample, tstatsdump, tstatstart, trestart, &
                                    nsv, itot, jtot, ktot, xlen, ylen, xlat, xlon, xday, xtime, lwalldist, &
                                    lmoist, lcoriol, igrw_damp, geodamptime, ifnamopt, fname_options, &
                                    nscasrc,nscasrcl,iwallmom,iwalltemp,iwallmoist,iwallscal,ipoiss,iadv_mom,iadv_tke,iadv_thl,iadv_qt,iadv_sv,courant,diffnr,ladaptive,author,&
                                    linoutflow, lper2inout, libm, lnudge, lnudgevel, tnudge, nnudge, lles, luoutflowr, lvoutflowr, luvolflowr, lvvolflowr, &
                                    uflowrate, vflowrate, lstoreplane, iplane, &
                                    lreadmean, iinletgen, inletav, lreadminl, Uinf, Vinf, linletRA, nblocks, &
                                    lscalrec,lSIRANEinout,lscasrc,lscasrcl,lscasrcr,lydump,lytdump,lxydump,lxytdump,ltdump,lmintdump,ltkedump,lzerogradtop,&
                                    lkslicedump,lislicedump,ljslicedump,kslice,islice,jslice,&
                                    lzerogradtopscal, lbuoyancy, ltempeq, &
                                    lfixinlet, lfixutauin, pi, &
                                    thlsrc, ifixuinf, lvinf, tscale, ltempinout, lmoistinout,  &
                                    lwallfunc,lprofforc,lchem,k1,JNO2,rv,rd,tnextEB,tEB,dtEB,bldT,flrT, lperiodicEBcorr, fraction,sinkbase,wsoil,wgrmax,wwilt,wfc,skyLW,GRLAI,rsmin,nfcts,lEB,lwriteEBfiles,nfaclyrs,lconstW,lvfsparse,nnz,lfacTlyrs, &
                                    BCxm,BCxT,BCxq,BCxs,BCym,BCyT,BCyq,BCys,BCzp,ds, &
                                    BCtopm,BCtopT,BCtopq,BCtops,BCbotm,BCbotT,BCbotq,BCbots, &
                                    BCxm_periodic, BCym_periodic, &
                                    idriver,tdriverstart,driverjobnr,dtdriver,driverstore,lchunkread,chunkread_size, &
                                    lrandomize, prandtlturb, fkar, lwritefac, dtfac, tfac, tnextfac, &
                                    ltrees,ntrees,Qstar,dQdt,lad,lsize,r_s,cd,dec,ud,ltreedump, &
                                    lpurif,npurif,Qpu,epu, &
                                    lheatpump,lfan_hp,nhppoints,Q_dot_hp,QH_dot_hp
      use modsurfdata,       only : z0, z0h,  wtsurf, wttop, wqtop, wqsurf, wsvsurf, wsvtop, wsvsurfdum, wsvtopdum, ps, thvs, thls, thl_top, qt_top, qts
      use modfields,         only : initfields, dpdx, ncname
      use modpois,           only : initpois
      use modboundary,       only : initboundary, ksp
      use modthermodynamics, only : initthermodynamics, lqlnr, chi_half
      use modsubgrid,        only : initsubgrid
      use modmpi,            only : comm3d, myid, myidx, myidy, cmyid, cmyidx, cmyidy, mpi_integer, mpi_logical, my_real, mpierr, mpi_character, nprocx, nprocy, nbreast, nbrwest, nbrnorth, nbrsouth
      use modinlet,          only : initinlet
      use modinletdata,      only : di, dr, di_test, dti, iangledeg, iangle
      use modibmdata,        only : bctfxm, bctfxp, bctfym, bctfyp, bctfz, bcqfxm, bcqfxp, bcqfym, bcqfyp, bcqfz
      use modforces,         only : calcfluidvolumes
      use moddriver,         only : initdriver
      use modtimedep,        only : ltimedepsurf, ntimedepsurf, ltimedepnudge, ntimedepnudge, &
                                    ltimedeplw, ntimedeplw, ltimedepsw, ntimedepsw
      use modibm,            only : nsolpts_u, nsolpts_v, nsolpts_w, nsolpts_c, &
                                    nbndpts_u, nbndpts_v, nbndpts_w, nbndpts_c, &
                                    nfctsecs_u, nfctsecs_v, nfctsecs_w, nfctsecs_c, &
                                    createmasks, lbottom, lnorec
      use decomp_2d

      implicit none
      integer :: ierr
      logical, dimension(3) :: periodic_bc
      integer, dimension(2) :: myids

      !declare namelists

      namelist/RUN/ &
         iexpnr, lwarmstart, lstratstart, startfile, &
         runtime, dtmax, trestart, ladaptive, &
         irandom, randu, randthl, randqt, krand, &
         courant, diffnr, author, &
         libm, lles, &
         lper2inout, lwalldist, &
         lreadmean, &
         nprocx, nprocy, &
         lrandomize
      namelist/DOMAIN/ &
         itot, jtot, ktot, xlen, ylen, &
         xlat, xlon, xday, xtime, ksp
      namelist/PHYSICS/ &
         ps, igrw_damp, lmoist, lcoriol, lbuoyancy, ltempeq, &
         lprofforc, ifixuinf, lvinf, tscale, dpdx, &
         luoutflowr, lvoutflowr, luvolflowr, lvvolflowr, &
         uflowrate, vflowrate, &
         lnudge, lnudgevel, tnudge, nnudge, &
         ltimedepsurf, ntimedepsurf, ltimedepnudge, ntimedepnudge, &
         ltimedeplw, ntimedeplw, ltimedepsw, ntimedepsw
      namelist/DYNAMICS/ &
         lqlnr, ipoiss, &
         iadv_mom, iadv_tke, iadv_thl, iadv_qt, iadv_sv
      namelist/BC/ &
         BCxm, BCxT, BCxq, BCxs, &
         BCym, BCyT, BCyq, BCys, &
         BCtopm, BCtopT, BCtopq, BCtops, &
         BCbotm, BCbotT, BCbotq, BCbots, &
         bctfxm, bctfxp, bctfym, bctfyp, bctfz, &
         bcqfxm, bcqfxp, bcqfym, bcqfyp, bcqfz, &
         wttop, thl_top, qt_top, qts, wsvsurfdum, wsvtopdum, &
         wtsurf, wqsurf, thls, z0, z0h, BCzp, ds
      namelist/INLET/ &
         Uinf, Vinf, di, dti, inletav, linletRA, &
         lstoreplane, lreadminl, lfixinlet, lfixutauin, &
         lwallfunc
      namelist/DRIVER/ &
         idriver, tdriverstart, driverjobnr, dtdriver, &
         driverstore, iplane, iangledeg, &
         lchunkread, chunkread_size
      namelist/WALLS/ &
         nblocks, nfcts, iwallmom, iwalltemp, iwallmoist, iwallscal, &
         nsolpts_u, nsolpts_v, nsolpts_w, nsolpts_c, &
         nbndpts_u, nbndpts_v, nbndpts_w, nbndpts_c, &
         nfctsecs_u, nfctsecs_v, nfctsecs_w, nfctsecs_c, lbottom, lnorec, &
         prandtlturb, fkar, lwritefac, dtfac
      namelist/ENERGYBALANCE/ &
         lEB, lwriteEBfiles, lperiodicEBcorr, sinkbase, lconstW, dtEB, bldT, flrT, wsoil, wgrmax, wwilt, wfc, &
         skyLW, GRLAI, rsmin, nfaclyrs, lfacTlyrs, lvfsparse, nnz, fraction
      namelist/SCALARS/ &
         lreadscal, lscasrc, lscasrcl, lscasrcr, &
         nsv, nscasrc, nscasrcl !!xS, yS, zS, SS, sigS
      namelist/CHEMISTRY/ &
         lchem, k1, JNO2
      namelist/OUTPUT/ &
         lfielddump, tfielddump, fieldvars, &
         ltdump, lydump, lytdump, lxydump, lxytdump, lmintdump, ltkedump, &
         slicevars, lkslicedump, kslice, lislicedump, islice, ljslicedump, jslice, &
         tstatsdump, tsample, tstatstart
      namelist/TREES/ &
         ltrees, ntrees, cd, dec, ud, lad, Qstar, dQdt, lsize, r_s, ltreedump
      namelist/PURIFS/ &
         lpurif, npurif, Qpu, epu
      namelist/HEATPUMP/ &
         lheatpump, lfan_hp, nhppoints, Q_dot_hp, QH_dot_hp

      if (myid == 0) then
         if (command_argument_count() >= 1) then
            call get_command_argument(1, fname_options)
         end if
         !write (*, *) fname_options

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
         !write (6, RUN)
         rewind (ifnamopt)

         read (ifnamopt, DOMAIN, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions DOMAIN'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         !write (6, DOMAIN)
         rewind (ifnamopt)

         read (ifnamopt, PHYSICS, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions PHYSICS'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         !write (6, PHYSICS)
         rewind (ifnamopt)

         read (ifnamopt, DYNAMICS, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions DYNAMICS'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         !write (6, DYNAMICS)
         rewind (ifnamopt)

         read (ifnamopt, BC, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions BC'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         !write (6, BC)
         rewind (ifnamopt)

         read (ifnamopt, INLET, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions INLET'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         !write (6, INLET)
         rewind (ifnamopt)

         read (ifnamopt, DRIVER, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'Problem in namoptions DRIVER'
            write(0, *) 'iostat error: ', ierr
            stop 'ERROR: Problem in namoptions DRIVER'
         endif
         !write (6, DRIVER)
         rewind (ifnamopt)

         read (ifnamopt, WALLS, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions WALLS'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         !write (6, WALLS)
         rewind (ifnamopt)

         read (ifnamopt, ENERGYBALANCE, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions EB'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         !write (6, ENERGYBALANCE)
         rewind (ifnamopt)

         read (ifnamopt, SCALARS, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions SCALARS'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         !write (6, SCALARS)
         rewind (ifnamopt)

         read (ifnamopt, CHEMISTRY, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions CHEMISTRY'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         !write (6, CHEMISTRY)
         rewind (ifnamopt)

         read (ifnamopt, TREES, iostat=ierr)
         if (ierr > 0) then
            print *, 'ERROR: Problem in namoptions TREES'
            print *, 'iostat error: ', ierr
            stop 1
         endif
         !write (6, TREES)
         rewind (ifnamopt)

         read (ifnamopt, PURIFS, iostat=ierr)
         if (ierr > 0) then
            print *, 'ERROR: Problem in namoptions PURIFS'
            print *, 'iostat error: ', ierr
            stop 1
         endif
         !write (6, PURIFS)
         rewind (ifnamopt)

         read (ifnamopt, HEATPUMP, iostat=ierr)
         if (ierr > 0) then
            print *, 'ERROR: Problem in namoptions HEATPUMP'
            print *, 'iostat error: ', ierr
            stop 1
         endif
         !write (6, HEATPUMP)
         rewind (ifnamopt)

         read (ifnamopt, OUTPUT, iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions OUTPUT'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         !write (6, OUTPUT)
         close (ifnamopt)
      end if

      call MPI_BCAST(itot,1,MPI_INTEGER,0,comm3d,mpierr)
      call MPI_BCAST(jtot,1,MPI_INTEGER,0,comm3d,mpierr)
      call MPI_BCAST(ktot,1,MPI_INTEGER,0,comm3d,mpierr)
      call MPI_BCAST(nprocx,1,MPI_INTEGER,0,comm3d,mpierr)
      call MPI_BCAST(nprocy,1,MPI_INTEGER,0,comm3d,mpierr)
      call MPI_BCAST(BCxm, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(BCym, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(BCzp, 1, MPI_INTEGER, 0, comm3d, mpierr)

      ! if (BCxm .eq. BCxm_periodic .and. nprocx > 1) then
      !   periodic_bc(1) = .true.
      ! else
      !   periodic_bc(1) = .false.
      ! end if
      !
      ! if (BCym .eq. BCym_periodic .and. nprocy > 1) then
      !   periodic_bc(2) = .true.
      ! else
      !   periodic_bc(2) = .false.
      ! end if
      !
      ! periodic_bc(3) = .false.
      ! call decomp_2d_init(itot,jtot,ktot,nprocx,nprocy,periodic_bc)
      ! !myid = nrank
      ! !write(cmyid,'(i3.3)') myid
      !
      ! comm3d = DECOMP_2D_COMM_CART_Z
      ! !write(*,*) "myid", myid
      ! call MPI_CART_COORDS(comm3d,myid,2,myids,mpierr)
      ! !write(*,*) "myids", myids
      ! myidx = myids(1)
      ! myidy = myids(2)
      ! ! write(*,*) "myid", " myids", myid, myids'
      !
      ! write(cmyidx,'(i3.3)') myidx
      ! write(cmyidy,'(i3.3)') myidy
      !
      ! call MPI_CART_SHIFT(comm3d, 0,  1, nbrwest,  nbreast ,   mpierr)
      ! call MPI_CART_SHIFT(comm3d, 1,  1, nbrsouth, nbrnorth,   mpierr)

      !call init2decomp

      !write (*, *) "starting broadcast"
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
      call MPI_BCAST(lnudgevel, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(nnudge, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(tnudge, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(ltimedepsurf, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ltimedepnudge, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ltimedeplw, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ltimedepsw, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ntimedepsurf, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(ntimedepnudge, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(ntimedeplw, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(ntimedepsw, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(lwalldist, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for computing wall distances
      call MPI_BCAST(lles, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for turning on/off LES functionality (subgrid model)
      call MPI_BCAST(linletRA, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for turning on/off Running Average in inletgenerator
      call MPI_BCAST(lfixinlet, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for keeping average inlet velocit and temp fixed at inlet (iinletgen=1,2)
      call MPI_BCAST(lfixutauin, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for keeping utau fixed at inlet (iinletgen=1,2)
      !call MPI_BCAST(xS, 1, MY_REAL, 0, comm3d, mpierr)
      !call MPI_BCAST(yS, 1, MY_REAL, 0, comm3d, mpierr)
      !call MPI_BCAST(zS, 1, MY_REAL, 0, comm3d, mpierr)
      !call MPI_BCAST(SS, 1, MY_REAL, 0, comm3d, mpierr)
      !call MPI_BCAST(sigS, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(idriver    ,1,MPI_INTEGER,0,comm3d,mpierr)        ! ae1212: Added switch for driver inlet simulation
      call MPI_BCAST(tdriverstart,1,MY_REAL   ,0,comm3d,mpierr)        ! ae1212
      call MPI_BCAST(driverjobnr,1,MPI_INTEGER,0,comm3d,mpierr)        ! ae1212
      call MPI_BCAST(dtdriver   ,1,MY_REAL    ,0,comm3d,mpierr)        ! ae1212
      call MPI_BCAST(driverstore,1,MPI_INTEGER ,0,comm3d,mpierr)
      call MPI_BCAST(lchunkread ,1,MPI_LOGICAL,0,comm3d,mpierr)
      call MPI_BCAST(chunkread_size,1,MPI_INTEGER,0,comm3d,mpierr)
      !call MPI_BCAST(BCxm, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(BCxT, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(BCxq, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(BCxs, 1, MPI_INTEGER, 0, comm3d, mpierr)
      !call MPI_BCAST(BCym, 1, MPI_INTEGER, 0, comm3d, mpierr)
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
      call MPI_BCAST(ds, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(lwallfunc, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for reading mean inlet/recycle plane profiles (Uinl,Urec,Wrec)
      call MPI_BCAST(lreadminl, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! J.Tomas: added switch for reading mean inlet/recycle plane profiles (Uinl,Urec,Wrec)
      call MPI_BCAST(iwalltemp, 1, MPI_INTEGER, 0, comm3d, mpierr) ! case (integer) for wall treatment for temperature (1=no wall function/fixed flux, 2=no wall function/fixed value, 3=uno)
      call MPI_BCAST(iwallmoist, 1, MPI_INTEGER, 0, comm3d, mpierr) ! case (integer) for wall treatment for moisture (1=no wall function/fixed flux, 2=no wall function/fixed value, 3=uno)
      call MPI_BCAST(iwallscal, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iwallmom, 1, MPI_INTEGER, 0, comm3d, mpierr) ! case (integer) for wall treatment for momentum (1=no wall function, 2=werner-wengle, 3=uno)
      call MPI_BCAST(nsolpts_u, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nsolpts_v, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nsolpts_w, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nsolpts_c, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nbndpts_u, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nbndpts_v, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nbndpts_w, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nbndpts_c, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nfctsecs_u, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nfctsecs_v, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nfctsecs_w, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nfctsecs_c, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(lbottom, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lnorec, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lwritefac, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(tfac, 1, MY_REAL, 0, comm3d, mpierr)
      tnextfac = dtfac
      call MPI_BCAST(tnextfac, 1, MY_REAL, 0, comm3d, mpierr)
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
      call MPI_BCAST(lkslicedump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lislicedump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ljslicedump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(kslice, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(islice, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(jslice, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(ltdump, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! tg3315 added switch for writing statistics files
      call MPI_BCAST(lmintdump, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! tg3315 added switch for writing statistics files
      call MPI_BCAST(ltkedump, 1, MPI_LOGICAL, 0, comm3d, mpierr) ! tg3315 added switch for writing tke budget files
      call MPI_BCAST(iplane, 1, MPI_INTEGER, 0, comm3d, mpierr) ! J.Tomas: ib+iplane is the i-plane that is stored if lstoreplane is .true.
      call MPI_BCAST(startfile, 50, MPI_CHARACTER, 0, comm3d, mpierr)
      call MPI_BCAST(author, 80, MPI_CHARACTER, 0, comm3d, mpierr)
      call MPI_BCAST(runtime, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(trestart, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(tfielddump, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(tsample, 1, MY_REAL, 0, comm3d, mpierr) !tg3315
      call MPI_BCAST(tstatsdump, 1, MY_REAL, 0, comm3d, mpierr) !tg3315
      call MPI_BCAST(tstatstart, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(tEB, 1, MY_REAL, 0, comm3d, mpierr)
      tnextEB = dtEB
      call MPI_BCAST(tnextEB, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(dtmax, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(nsv, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nscasrc,1,MPI_INTEGER,0,comm3d,mpierr)
      call MPI_BCAST(nscasrcl,1,MPI_INTEGER,0,comm3d,mpierr)
      call MPI_BCAST(fieldvars, 50, MPI_CHARACTER, 0, comm3d, mpierr)
      call MPI_BCAST(slicevars, 50, MPI_CHARACTER, 0, comm3d, mpierr)
      !call MPI_BCAST(nstat      ,1,MPI_INTEGER,0,comm3d,mpierr) !tg3315
      !call MPI_BCAST(ncstat     ,80,MPI_CHARACTER,0,comm3d,mpierr) !tg3315
      call MPI_BCAST(ifixuinf, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(lvinf, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(dpdx, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(tscale, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(itot, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(jtot, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(ktot, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(xlen, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(ylen, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(xlat, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(xlon, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(xday, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(xtime, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(z0, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(z0h, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctfxm, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctfxp, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctfym, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctfyp, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctfz, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bcqfxm, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bcqfxp, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bcqfym, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bcqfyp, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bcqfz, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wtsurf, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wqsurf, 1, MY_REAL, 0, comm3d, mpierr)
      allocate (wsvsurf(1:nsv))
      wsvsurf = wsvsurfdum(1:nsv)
      call MPI_BCAST(wsvsurf(1:nsv), nsv, MY_REAL, 0, comm3d, mpierr)
      allocate (wsvtop(1:nsv))
      wsvtop = wsvtopdum(1:nsv)
      call MPI_BCAST(wsvtop(1:nsv), nsv, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(ps, 1, MY_REAL, 0, comm3d, mpierr)
      thvs = thls*(1.+(rv/rd - 1.)*qts)
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
      call MPI_BCAST(flrT, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(lperiodicEBcorr, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(sinkbase, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(fraction, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(skyLW, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(GRLAI, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(rsmin, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(nfaclyrs, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(lfacTlyrs, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lvfsparse, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(nnz, 1, MPI_INTEGER, 0, comm3d, mpierr)
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
      call MPI_BCAST(lrandomize, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(prandtlturb, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(fkar, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(ltrees, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ntrees, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(ltreedump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(Qstar, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(dQdt, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(lsize, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(lad, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(r_s, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(cd, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(dec, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(ud, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(lpurif, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(npurif, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(Qpu, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(epu, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(lheatpump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lfan_hp, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(nhppoints, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(Q_dot_hp, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(QH_dot_hp, 1, MY_REAL, 0, comm3d, mpierr)

      ! ! Allocate and initialize core modules
      ! call initglobal
      ! !write (*, *) "done initglobal"
      ! call initfields
      ! !write (*, *) "done initfields"
      ! call initboundary
      ! !write (*, *) "done initboundary"
      ! call initthermodynamics
      ! ! write (*, *) "done initthermodynamics"
      ! ! !!depreated!!
      ! ! ! call initsurface
      ! ! write (*, *) "done initsurface"
      ! call initsubgrid
      ! !write (*, *) "done initsubgrid"
      ! ! call initpois
      ! ! write (*, *) "done initpois"
      ! ! call initinlet ! added by J. Tomas: initialize inlet generator
      ! ! write (*, *) "done initinlet"
      ! call initdriver  ! added by ae1212: initialise driver inlet
      ! ! write(*,*) "done initdriver"
      ! call checkinitvalues
      ! !write (*, *) "done checkinitvalues"
      ! call initpois
      ! !write (*, *) "done initpois"
      ! ! write (6, *) 'Determine masking matrices'
      ! call createmasks ! determine walls/blocks
      ! ! write (6, *) 'Finished determining masking matrices'
      ! ! ! calculate fluid volume and outlet areas, needs masking matrices
      ! call calcfluidvolumes
      ! !
      ! call readinitfiles
      ! !write (*, *) "done readinitfiles"
      ! ! write (*, *) "done startup"
      ! !
      ! ! call createscals
      ! ! write (*, *) "done create scals"

   end subroutine readnamelists

end module readparameters
