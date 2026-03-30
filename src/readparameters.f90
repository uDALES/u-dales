!> \file readparameters.f90
!!  Reads and broadcasts configuration parameters for uDALES

!>
!! uDALES: Urbanized Dutch Atmospheric Large-Eddy Simulation
!>
!! This module handles reading of configuration parameters from Fortran
!! namelists and broadcasting them across all MPI processes. It supports
!! modular reading and broadcasting for each configuration section, and
!! provides routines for writing current configuration to file for
!! verification.
!!  \author Maarten van Reeuwijk, Imperial College London
!!  \author uDALES contributors
!!  \par Revision list
!
!  This file is part of uDALES.
!
! uDALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! uDALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 2017-2024 uDALES contributors, Imperial College London
!

module readparameters

use mpi
   use modglobal, mg_rv => rv
   use modsurfdata
   use modfields
   use modpois
   use modboundary
   use modthermodynamics
   use modsubgrid
   use modmpi
   use modinlet
   use modinletdata
   use modibmdata
   use modforces
   use moddriver
   use modtimedep
   use modibm
   use decomp_2d
   use modsubgriddata, sg_cs => cs

   implicit none
   ! private
   ! public :: startup,trestart
!   public :: RUN, DOMAIN, PHYSICS, DYNAMICS, BC, INLET, DRIVER, WALLS, ENERGYBALANCE, SCALARS, CHEMISTRY, OUTPUT, TREES, PURIFS, HEATPUMP
   public :: writenamelists
   save

   integer(KIND=selected_int_kind(6)) :: irandom = 43 !    * number to seed the randomnizer with
   integer :: krand = huge(0)  ! returns the largest integer that is not an infinity
   real :: randu = 0.01, randthl = 0.0, randqt = 0.0 !    * uvw,thl and qt amplitude of randomnization

   ! Module-level namelist definitions
   namelist/RUN/ &
      iexpnr, lwarmstart, lstratstart, startfile, &
      runtime, dtmax, trestart, ladaptive, &
      irandom, randu, randthl, randqt, krand, &
      courant, diffnr, author, &
      libm, lles, &
      lper2inout, lwalldist, &
      lreadmean, &
      nprocx, nprocy, &
      lrandomize, runmode
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
      nsv, nscasrc, nscasrcl
   namelist/CHEMISTRY/ &
      lchem, k1, JNO2
   namelist/OUTPUT/ &
      lfielddump, tfielddump, fieldvars, &
      ltdump, lydump, lytdump, lxydump, lxytdump, lmintdump, ltkedump, &
      slicevars, lkslicedump, kslice, lislicedump, islice, ljslicedump, jslice, &
      tstatsdump, tsample, tstatstart, tcheck
   namelist/TREES/ &
      ltrees, ntrees, cd, dec, ud, lad, Qstar, dQdt, lsize, r_s, ltreedump
   namelist/PURIFS/ &
      lpurif, npurif, Qpu, epu
   namelist/HEATPUMP/ &
      lheatpump, lfan_hp, nhppoints, Q_dot_hp, QH_dot_hp
   namelist/NAMSUBGRID/ &
         ldelta,lmason, cf,cn,Rigc,Prandtl,lsmagorinsky,lvreman,loneeqn,c_vreman,sg_cs,nmason,lbuoycorr
   namelist/INFO/ &
         nprocsinl,jgtotinl,kmaxin,dtin,wtop,totalreadu

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

      implicit none
      integer :: ierr
      logical, dimension(3) :: periodic_bc
      integer, dimension(2) :: myids

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
         rewind (ifnamopt)

         read (ifnamopt,NAMSUBGRID,iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions NAMSUBGRID'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         rewind (ifnamopt)
         !write(6 ,NAMSUBGRID)
 
         if (has_namelist_group('INFO')) then
            read (ifnamopt,INFO,iostat=ierr)
            if (ierr > 0) then
               write(0, *) 'Problem in namoptions INFO'
               write(0, *) 'iostat error: ', ierr
               stop 1
            endif
            !write(6,INFO)
         endif

         if (lrandomize) then
            if (.not. lles) then
               write(0, *) 'Cannot use lrandomize=.true. with lles=.false.'
               stop 1
            end if
            if (irandom == 0) then
               write(0, *) 'Please define RUN.irandom when lrandomize=.true.'
               stop 1
            end if
         end if

         if (any(iadv_sv(1:nsv) == -1)) then
            write(0, *) 'ERROR: iadv_sv contains -1, which is invalid. Please check your configuration.'
            stop 1
         end if

         close(ifnamopt)
      end if
   end subroutine readnamelists

   subroutine writenamelists
      !-----------------------------------------------------------------|
      !     Write all namelists to file to verify broadcast worked     |
      !-----------------------------------------------------------------|
      
      implicit none
      integer :: iunit = 99
      character(len=100) :: filename
      
      ! Create filename with process ID for debugging
      write(filename,'(a)') 'namoptions_roundtrip'
      write(*,*) '  Writing all namelists to file: ', trim(filename)
      
      open(unit=iunit, file=filename, status='replace', action='write')
      
      ! Write each namelist section using Fortran's built-in namelist write
      write(iunit, nml=RUN)
      write(iunit, nml=DOMAIN) 
      write(iunit, nml=PHYSICS)
      write(iunit, nml=DYNAMICS)
      write(iunit, nml=BC)
      write(iunit, nml=INLET)
      write(iunit, nml=DRIVER)
      write(iunit, nml=WALLS)
      write(iunit, nml=ENERGYBALANCE)
      write(iunit, nml=SCALARS)
      write(iunit, nml=CHEMISTRY)
      write(iunit, nml=OUTPUT)
      write(iunit, nml=TREES)
      write(iunit, nml=PURIFS)
      write(iunit, nml=HEATPUMP)
      write(iunit, nml=INFO)
      write(iunit, nml=NAMSUBGRID)
      close(iunit)   
     end subroutine writenamelists
   subroutine broadcast_config_parameters
      !-----------------------------------------------------------------|
      !                                                                 |
      !     Broadcast all configuration parameters using modular funcs  |
      !     Called after the root process reads the namelists           |
      !                                                                 |
      !-----------------------------------------------------------------|
      
      implicit none
      
      ! Use modular broadcast functions for better organization
      call broadcast_run_parameters()
      call broadcast_domain_parameters()
      call broadcast_physics_parameters()
      call broadcast_scalars_parameters()
      call broadcast_dynamics_parameters()
      call broadcast_bc_parameters()
      call broadcast_driver_parameters()
      call broadcast_energybalance_parameters()
      call broadcast_heatpump_parameters()
      call broadcast_inlet_parameters()
      call broadcast_output_parameters()
      call broadcast_purifs_parameters()
      call broadcast_chemistry_parameters()
      call broadcast_trees_parameters()
      call broadcast_walls_parameters()
      if (has_namelist_group('INFO')) call broadcast_info_parameters()
      call broadcast_namsubgrid_parameters()
      
   end subroutine broadcast_config_parameters

   !-----------------------------------------------------------------|
   !              MODULAR PARAMETER BROADCAST FUNCTIONS             |
   !-----------------------------------------------------------------|
   subroutine broadcast_run_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast RUN namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(iexpnr, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(lwarmstart, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lstratstart, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(startfile, 50, MPI_CHARACTER, 0, comm3d, mpierr)
      call MPI_BCAST(runtime, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(dtmax, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(trestart, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(ladaptive, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(irandom, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(randu, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(randthl, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(randqt, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(krand, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(courant, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(diffnr, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(author, 50, MPI_CHARACTER, 0, comm3d, mpierr)
      call MPI_BCAST(libm, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lles, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lper2inout, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lwalldist, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lreadmean, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(nprocx, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nprocy, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(lrandomize, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(runmode, 1, MPI_INTEGER, 0, comm3d, mpierr)
   end subroutine broadcast_run_parameters
   subroutine broadcast_domain_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast DOMAIN namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(itot, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(jtot, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(ktot, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(xlen, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(ylen, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(xlat, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(xlon, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(xday, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(xtime, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(ksp, 1, MPI_INTEGER, 0, comm3d, mpierr)
   end subroutine broadcast_domain_parameters
   subroutine broadcast_physics_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast PHYSICS namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(ps, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(igrw_damp, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(lmoist, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lcoriol, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lbuoyancy, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ltempeq, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lprofforc, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ifixuinf, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(lvinf, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(tscale, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(dpdx, 1, MY_REAL, 0, comm3d, mpierr)
      ! Broadcast PHYSICS variables that are part of the input interface
      call MPI_BCAST(lnudge, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lnudgevel, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ltimedeplw, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ltimedepnudge, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ltimedepsurf, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ltimedepsw, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(luoutflowr, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(luvolflowr, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lvoutflowr, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lvvolflowr, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(nnudge, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(ntimedeplw, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(ntimedepnudge, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(ntimedepsurf, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(ntimedepsw, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(tnudge, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(uflowrate, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(vflowrate, 1, MY_REAL, 0, comm3d, mpierr)
   end subroutine broadcast_physics_parameters
   subroutine broadcast_dynamics_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast DYNAMICS namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(lqlnr, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(iadv_mom, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iadv_tke, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iadv_thl, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iadv_qt, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iadv_sv, size(iadv_sv), MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(ipoiss, 1, MPI_INTEGER, 0, comm3d, mpierr)
   end subroutine broadcast_dynamics_parameters
   subroutine broadcast_bc_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast BC namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(bcbotm, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(bcbotq, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(bcbots, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(bcbott, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(bcqfxm, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bcqfxp, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bcqfym, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bcqfyp, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bcqfz, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctfxm, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctfxp, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctfym, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctfyp, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctfz, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(bctopm, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(bctopq, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(bctops, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(bctopt, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(bcxm, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(bcxq, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(bcxs, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(bcxt, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(bcym, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(bcyq, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(bcys, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(bcyt, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(bczp, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(ds, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(qt_top, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(qts, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(thl_top, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(thls, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wqsurf, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wsvsurfdum, size(wsvsurfdum), MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wsvtopdum, size(wsvtopdum), MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wtsurf, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wttop, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(z0, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(z0h, 1, MY_REAL, 0, comm3d, mpierr)
   end subroutine broadcast_bc_parameters
   subroutine broadcast_driver_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast DRIVER namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(idriver, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(chunkread_size, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(driverjobnr, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(driverstore, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(dtdriver, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(iangledeg, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(iplane, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(lchunkread, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(tdriverstart, 1, MY_REAL, 0, comm3d, mpierr)
   end subroutine broadcast_driver_parameters
   subroutine broadcast_energybalance_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast ENERGYBALANCE namelist parameters to all processes
      !-----------------------------------------------------------------|
      ! Broadcast all ENERGYBALANCE input variables
      call MPI_BCAST(bldt, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(dteb, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(flrt, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(fraction, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(grlai, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(lconstw, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(leb, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lfactlyrs, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lperiodicebcorr, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lvfsparse, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lwriteebfiles, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(nfaclyrs, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nnz, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(rsmin, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(sinkbase, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(skylw, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wfc, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wgrmax, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wsoil, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wwilt, 1, MY_REAL, 0, comm3d, mpierr)
   end subroutine broadcast_energybalance_parameters
   subroutine broadcast_heatpump_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast HEATPUMP namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(lheatpump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lfan_hp, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(nhppoints, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(Q_dot_hp, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(QH_dot_hp, 1, MY_REAL, 0, comm3d, mpierr)
   end subroutine broadcast_heatpump_parameters
   subroutine broadcast_inlet_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast INLET namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(di, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(dti, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(inletav, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(lfixinlet, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lfixutauin, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(linletra, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lreadminl, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lstoreplane, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lwallfunc, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(uinf, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(vinf, 1, MY_REAL, 0, comm3d, mpierr)
   end subroutine broadcast_inlet_parameters
   subroutine broadcast_output_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast OUTPUT namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(lfielddump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(tfielddump, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(fieldvars, 50, MPI_CHARACTER, 0, comm3d, mpierr)
      call MPI_BCAST(tsample, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(tstatsdump, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(tstatstart, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(tcheck, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(islice, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(jslice, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(kslice, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(slicevars, 50, MPI_CHARACTER, 0, comm3d, mpierr)
      call MPI_BCAST(lislicedump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ljslicedump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lkslicedump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lmintdump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ltdump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ltkedump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lxydump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lxytdump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lydump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lytdump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
   end subroutine broadcast_output_parameters
   subroutine broadcast_purifs_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast PURIFS namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(lpurif, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(npurif, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(epu, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(qpu, 1, MY_REAL, 0, comm3d, mpierr)
   end subroutine broadcast_purifs_parameters
   subroutine broadcast_scalars_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast SCALARS namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(lreadscal, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lscasrc, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lscasrcl, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lscasrcr, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(nscasrc, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nscasrcl, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nsv, 1, MPI_INTEGER, 0, comm3d, mpierr)
   end subroutine broadcast_scalars_parameters
   subroutine broadcast_chemistry_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast CHEMISTRY namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(lchem, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(k1, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(JNO2, 1, MY_REAL, 0, comm3d, mpierr)
   end subroutine broadcast_chemistry_parameters
   subroutine broadcast_trees_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast TREES namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(ltrees, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ltreedump, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(ntrees, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(cd, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(dec, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(dqdt, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(lad, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(lsize, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(qstar, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(r_s, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(ud, 1, MY_REAL, 0, comm3d, mpierr)
   end subroutine broadcast_trees_parameters
   subroutine broadcast_walls_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast WALLS namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(dtfac, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(fkar, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(iwallmoist, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iwallmom, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iwallscal, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iwalltemp, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(lbottom, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lnorec, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lwritefac, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(nblocks, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nbndpts_c, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nbndpts_u, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nbndpts_v, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nbndpts_w, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nfcts, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nfctsecs_c, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nfctsecs_u, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nfctsecs_v, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nfctsecs_w, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nsolpts_c, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nsolpts_u, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nsolpts_v, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(nsolpts_w, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(prandtlturb, 1, MY_REAL, 0, comm3d, mpierr)
   end subroutine broadcast_walls_parameters
   subroutine broadcast_info_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast INFO namelist parameters to all processes
      !-----------------------------------------------------------------|
      use modinletdata, only: nprocsinl, jgtotinl, kmaxin, dtin, wtop, totalreadu
      
      call MPI_BCAST(nprocsinl, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(jgtotinl, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(kmaxin, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(dtin, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wtop, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(totalreadu, 1, MY_REAL, 0, comm3d, mpierr)
   end subroutine broadcast_info_parameters
   subroutine broadcast_namsubgrid_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast NAMSUBGRID namelist parameters to all processes
      !-----------------------------------------------------------------|
      use modsubgrid, only: ldelta, lmason, cf, cn, Rigc, Prandtl, lsmagorinsky, &
                                lvreman, loneeqn, c_vreman, cs, nmason, lbuoycorr
      
      call MPI_BCAST(ldelta, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lmason, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(cf, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(cn, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(Rigc, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(Prandtl, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(lsmagorinsky, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(lvreman, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(loneeqn, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(c_vreman, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(cs, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(nmason, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(lbuoycorr, 1, MPI_LOGICAL, 0, comm3d, mpierr)
   end subroutine broadcast_namsubgrid_parameters

   logical function has_namelist_group(group_name)
      implicit none
      character(len=*), intent(in) :: group_name
      integer :: scan_unit, ioerr
      character(len=256) :: line

      has_namelist_group = .false.
      scan_unit = 98

      open(unit=scan_unit, file=fname_options, status='old', action='read', iostat=ioerr)
      if (ioerr /= 0) return

      do
         read(scan_unit, '(A)', iostat=ioerr) line
         if (ioerr /= 0) exit
         if (index(adjustl(line), '&' // trim(group_name)) == 1) then
            has_namelist_group = .true.
            exit
         end if
      end do

      close(scan_unit)
   end function has_namelist_group
end module readparameters
