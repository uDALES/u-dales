!> \file modstartup.f90
!!  Initializes the run

!>
!! Initializes the run.
!>
!! Modstartup reads the namelists and initial data, sets the fields and calls
!! the inits of the other routines where necessary. Reading of the
!! restart files also live in this module. Also supports JSON input format
!! when compiled with USE_JSON_INPUT.
!!  \author Maarten van Reeuwijk, Imperial College London
!!  \author Jasper Tomas, TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \todo documentation
!!  \par Revision list
!
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

module readparameters

use mpi
   use json_module, wp => json_RK, CK => json_CK
   use iso_fortran_env, only: real64
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
   public :: readnamelists, readjsonconfig, broadcast_config_parameters, writenamelists
   save

   !integer(KIND=selected_int_kind(6)) :: irandom = 43 !    * number to seed the randomnizer with
   !integer :: krand = huge(0)  ! returns the largest integer that is not an infinity
   !real :: randu = 0.01, randthl = 0.0, randqt = 0.0 !    * uvw,thl and qt amplitude of randomnization

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
         close (ifnamopt)

         open(ifnamopt,file=fname_options,status='old',iostat=ierr)
         read (ifnamopt,NAMSUBGRID,iostat=ierr)
         if (ierr > 0) then
            write(0, *) 'ERROR: Problem in namoptions NAMSUBGRID'
            write(0, *) 'iostat error: ', ierr
            stop 1
         endif
         !write(6 ,NAMSUBGRID)
         close(ifnamopt)
 
         read (ifinput,INFO,iostat=ierr)
         if (ierr > 0) then
           write(0, *) 'Problem in zgrid.inf INFO'
           write(0, *) 'iostat error: ', ierr
           stop 1
         endif
         !write(6,INFO)
         close(ifinput)
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
      write(filename,'(a,i0,a)') 'namoptions_json_test.', iexpnr
      
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
      
      write(*,*) '  All namelists written to file successfully!'
      
    end subroutine writenamelists

   subroutine readjsonconfig

      !-----------------------------------------------------------------|
      !                                                                 |
      !     Reads configuration from JSON file using modular functions |
      !                                                                 |
      !     Uses existing module defaults when parameters not found    |
      !-----------------------------------------------------------------|

      use modmpi,            only : myid, comm3d, mpierr, my_real, mpi_logical, mpi_integer

      type(json_file) :: json
      logical :: found_any
      character(len=80) :: json_fname

      if (myid == 0) then
         if (command_argument_count() >= 1) then
            call get_command_argument(1, json_fname)
         end if

         ! Initialize JSON parser
         call json%initialize()
         
         ! Load JSON file
         call json%load_file(json_fname)
         if (json%failed()) then
            write(0,*) 'ERROR: Could not load JSON file: ', trim(json_fname)
            call json%print_error_message()
            stop 1
         end if

         ! Read configuration using modular functions
         call read_run_json(json)
         call read_domain_json(json)
         call read_physics_json(json)
         call read_dynamics_json(json)
         call read_bc_json(json)
         call read_driver_json(json)
         call read_energybalance_json(json)
         call read_heatpump_json(json)
         call read_inlet_json(json)
         call read_output_json(json)
         call read_purifs_json(json)
         call read_scalars_json(json)
         call read_chemistry_json(json)
         call read_trees_json(json)
         call read_walls_json(json)
         call read_info_json(json)
         !call read_namchecksim_json(json)
         !call read_namstatsdump_json(json)
         call read_namsubgrid_json(json)

         ! Clean up
         call json%destroy()
      end if
   end subroutine readjsonconfig

   subroutine broadcast_config_parameters
      !-----------------------------------------------------------------|
      !                                                                 |
      !     Broadcast all configuration parameters using modular funcs  |
      !     Called by both readnamelists and readjsonconfig             |
      !                                                                 |
      !-----------------------------------------------------------------|
      
      implicit none
      
      ! Use modular broadcast functions for better organization
      call broadcast_run_parameters()
      call broadcast_domain_parameters()
      call broadcast_physics_parameters()
      call broadcast_dynamics_parameters()
      call broadcast_bc_parameters()
      call broadcast_driver_parameters()
      call broadcast_energybalance_parameters()
      call broadcast_heatpump_parameters()
      call broadcast_inlet_parameters()
      call broadcast_output_parameters()
      call broadcast_purifs_parameters()
      call broadcast_scalars_parameters()
      call broadcast_chemistry_parameters()
      call broadcast_trees_parameters()
      call broadcast_walls_parameters()
      call broadcast_info_parameters()
      !call broadcast_namchecksim_parameters()
      !call broadcast_namstatsdump_parameters()
      call broadcast_namsubgrid_parameters()
      
      ! Legacy broadcasts for remaining parameters (to be modularized)
      ! TODO: Modularize remaining namelist broadcasts
      
   end subroutine broadcast_config_parameters

   !-----------------------------------------------------------------|
   !                 MODULAR JSON READING FUNCTIONS                 |
   !-----------------------------------------------------------------|
   
   subroutine read_run_json(json)
      !-----------------------------------------------------------------|
      ! Read RUN namelist parameters from JSON
      !-----------------------------------------------------------------|
      type(json_file), intent(inout) :: json
      logical :: found
      integer :: temp_int
      real :: temp_real
      logical :: temp_logical
      character(kind=CK,len=:), allocatable :: temp_string

      call json%get('RUN.iexpnr', temp_int, found)
      if (found) iexpnr = temp_int
      
      call json%get('RUN.lwarmstart', temp_logical, found)
      if (found) lwarmstart = temp_logical

      call json%get('RUN.lstratstart', temp_logical, found)
      if (found) lstratstart = temp_logical

      call json%get('RUN.startfile', temp_string, found)
      if (found) startfile = temp_string
      
      call json%get('RUN.runtime', temp_real, found)
      if (found) runtime = temp_real

      call json%get('RUN.trestart', temp_real, found)
      if (found) trestart = temp_real

      call json%get('RUN.dtmax', temp_real, found)
      if (found) dtmax = temp_real
      
      call json%get('RUN.ladaptive', temp_logical, found)
      if (found) ladaptive = temp_logical

      call json%get('RUN.courant', temp_real, found)
      if (found) courant = temp_real

      call json%get('RUN.diffnr', temp_real, found)
      if (found) diffnr = temp_real
      
      call json%get('RUN.author', temp_string, found)
      if (found) author = temp_string

      call json%get('RUN.lrandomize', temp_logical, found)
      if (found) lrandomize = temp_logical

      call json%get('RUN.irandom', temp_int, found)
      if (found) irandom = temp_int

      call json%get('RUN.randu', temp_real, found)
      if (found) randu = temp_real

      call json%get('RUN.randthl', temp_real, found)
      if (found) randthl = temp_real
      
      call json%get('RUN.randqt', temp_real, found)
      if (found) randqt = temp_real

      call json%get('RUN.krand', temp_int, found)
      if (found) krand = temp_int

      call json%get('RUN.libm', temp_logical, found)
      if (found) libm = temp_logical
      
      call json%get('RUN.lles', temp_logical, found)
      if (found) lles = temp_logical

      call json%get('RUN.lper2inout', temp_logical, found)
      if (found) lper2inout = temp_logical

      call json%get('RUN.lreadmean', temp_logical, found)
      if (found) lreadmean = temp_logical
      
      call json%get('RUN.lwalldist', temp_logical, found)
      if (found) lwalldist = temp_logical

      ! Required parameters with error handling
      call json%get('RUN.nprocx', temp_int, found)
      if (found) then
         nprocx = temp_int
      else
         write(0,*) 'ERROR: nprocx is required but not found in JSON'
         stop 1
      end if
      
      call json%get('RUN.nprocy', temp_int, found)
      if (found) then
         nprocy = temp_int
      else
         write(0,*) 'ERROR: nprocy is required but not found in JSON'
         stop 1
      end if
      
      call json%get('RUN.runmode', temp_int, found)
      if (found) runmode = temp_int

   end subroutine read_run_json

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

   subroutine read_domain_json(json)
      !-----------------------------------------------------------------|
      ! Read DOMAIN namelist parameters from JSON
      !-----------------------------------------------------------------|
      type(json_file), intent(inout) :: json
      logical :: found
      integer :: temp_int
      real :: temp_real

      call json%get('DOMAIN.itot', temp_int, found)
      if (found) itot = temp_int
      
      call json%get('DOMAIN.jtot', temp_int, found)
      if (found) jtot = temp_int

      call json%get('DOMAIN.ktot', temp_int, found)
      if (found) ktot = temp_int

      call json%get('DOMAIN.xlen', temp_real, found)
      if (found) xlen = temp_real

      call json%get('DOMAIN.ylen', temp_real, found)
      if (found) ylen = temp_real

      call json%get('DOMAIN.xlat', temp_real, found)
      if (found) xlat = temp_real
      
      call json%get('DOMAIN.xlon', temp_real, found)
      if (found) xlon = temp_real

      call json%get('DOMAIN.xday', temp_real, found)
      if (found) xday = temp_real

      call json%get('DOMAIN.xtime', temp_real, found)
      if (found) xtime = temp_real
      
      call json%get('DOMAIN.ksp', temp_int, found)
      if (found) ksp = temp_int

   end subroutine read_domain_json

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

   subroutine read_physics_json(json)
      !-----------------------------------------------------------------|
      ! Read PHYSICS namelist parameters from JSON
      !-----------------------------------------------------------------|
      type(json_file), intent(inout) :: json
      logical :: found
      integer :: temp_int
      real :: temp_real
      logical :: temp_logical

      call json%get('PHYSICS.ps', temp_real, found)
      if (found) ps = temp_real
      
      call json%get('PHYSICS.ltempeq', temp_logical, found)
      if (found) ltempeq = temp_logical

      call json%get('PHYSICS.lbuoyancy', temp_logical, found)
      if (found) lbuoyancy = temp_logical

      call json%get('PHYSICS.lmoist', temp_logical, found)
      if (found) lmoist = temp_logical
      
      call json%get('PHYSICS.lcoriol', temp_logical, found)
      if (found) lcoriol = temp_logical

      call json%get('PHYSICS.igrw_damp', temp_int, found)
      if (found) igrw_damp = temp_int

      call json%get('PHYSICS.lprofforc', temp_logical, found)
      if (found) lprofforc = temp_logical
      
      call json%get('PHYSICS.ifixuinf', temp_int, found)
      if (found) ifixuinf = temp_int

      call json%get('PHYSICS.lvinf', temp_logical, found)
      if (found) lvinf = temp_logical

      call json%get('PHYSICS.tscale', temp_real, found)
      if (found) tscale = temp_real
      
      call json%get('PHYSICS.dpdx', temp_real, found)
      if (found) dpdx = temp_real

      call json%get('PHYSICS.lnudge', temp_logical, found)
      if (found) lnudge = temp_logical

      call json%get('PHYSICS.lnudgevel', temp_logical, found)
      if (found) lnudgevel = temp_logical
      
      call json%get('PHYSICS.ltimedeplw', temp_logical, found)
      if (found) ltimedeplw = temp_logical

      call json%get('PHYSICS.ltimedepnudge', temp_logical, found)
      if (found) ltimedepnudge = temp_logical

      call json%get('PHYSICS.ltimedepsurf', temp_logical, found)
      if (found) ltimedepsurf = temp_logical
      
      call json%get('PHYSICS.ltimedepsw', temp_logical, found)
      if (found) ltimedepsw = temp_logical

      call json%get('PHYSICS.luoutflowr', temp_logical, found)
      if (found) luoutflowr = temp_logical

      call json%get('PHYSICS.luvolflowr', temp_logical, found)
      if (found) luvolflowr = temp_logical
      
      call json%get('PHYSICS.lvoutflowr', temp_logical, found)
      if (found) lvoutflowr = temp_logical

      call json%get('PHYSICS.lvvolflowr', temp_logical, found)
      if (found) lvvolflowr = temp_logical

      call json%get('PHYSICS.nnudge', temp_int, found)
      if (found) nnudge = temp_int
      
      call json%get('PHYSICS.ntimedeplw', temp_int, found)
      if (found) ntimedeplw = temp_int

      call json%get('PHYSICS.ntimedepnudge', temp_int, found)
      if (found) ntimedepnudge = temp_int

      call json%get('PHYSICS.ntimedepsurf', temp_int, found)
      if (found) ntimedepsurf = temp_int
      
      call json%get('PHYSICS.ntimedepsw', temp_int, found)
      if (found) ntimedepsw = temp_int

      call json%get('PHYSICS.tnudge', temp_real, found)
      if (found) tnudge = temp_real

      call json%get('PHYSICS.uflowrate', temp_real, found)
      if (found) uflowrate = temp_real
      
      call json%get('PHYSICS.vflowrate', temp_real, found)
      if (found) vflowrate = temp_real

   end subroutine read_physics_json

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
      ! Add missing broadcasts for PHYSICS variables with JSON support
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

   subroutine read_dynamics_json(json)
      !-----------------------------------------------------------------|
      ! Read DYNAMICS namelist parameters from JSON
      !-----------------------------------------------------------------|
      type(json_file), intent(inout) :: json
      logical :: found
      integer :: temp_int
      logical :: temp_logical

      call json%get('DYNAMICS.iadv_mom', temp_int, found)
      if (found) iadv_mom = temp_int
      
      call json%get('DYNAMICS.iadv_tke', temp_int, found)
      if (found) iadv_tke = temp_int

      call json%get('DYNAMICS.iadv_thl', temp_int, found)
      if (found) iadv_thl = temp_int

      call json%get('DYNAMICS.iadv_qt', temp_int, found)
      if (found) iadv_qt = temp_int
      
      call json%get('DYNAMICS.iadv_sv', temp_int, found)
      if (found) iadv_sv = temp_int

      call json%get('DYNAMICS.ipoiss', temp_int, found)
      if (found) ipoiss = temp_int

      call json%get('DYNAMICS.lqlnr', temp_logical, found)
      if (found) lqlnr = temp_logical

   end subroutine read_dynamics_json

   subroutine broadcast_dynamics_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast DYNAMICS namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(lqlnr, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(iadv_mom, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iadv_tke, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iadv_thl, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iadv_qt, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iadv_sv, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(ipoiss, 1, MPI_INTEGER, 0, comm3d, mpierr)
   end subroutine broadcast_dynamics_parameters

   subroutine read_bc_json(json)
      !-----------------------------------------------------------------|
      ! Read BC namelist parameters from JSON  
      !-----------------------------------------------------------------|
      type(json_file), intent(inout) :: json
      logical :: found
      integer :: temp_int
      real :: temp_real
      logical :: temp_logical

      call json%get('BC.bcbotm', temp_int, found)
      if (found) bcbotm = temp_int

      call json%get('BC.bcbotq', temp_int, found)
      if (found) bcbotq = temp_int
      
      call json%get('BC.bcbots', temp_int, found)
      if (found) bcbots = temp_int

      call json%get('BC.bcbott', temp_int, found)
      if (found) bcbott = temp_int
      
      call json%get('BC.bcqfxm', temp_real, found)
      if (found) bcqfxm = temp_real

      call json%get('BC.bcqfxp', temp_real, found)
      if (found) bcqfxp = temp_real

      call json%get('BC.bcqfym', temp_real, found)
      if (found) bcqfym = temp_real
      
      call json%get('BC.bcqfyp', temp_real, found)
      if (found) bcqfyp = temp_real

      call json%get('BC.bcqfz', temp_real, found)
      if (found) bcqfz = temp_real

      call json%get('BC.bctfxm', temp_real, found)
      if (found) bctfxm = temp_real
      
      call json%get('BC.bctfxp', temp_real, found)
      if (found) bctfxp = temp_real

      call json%get('BC.bctfym', temp_real, found)
      if (found) bctfym = temp_real

      call json%get('BC.bctfyp', temp_real, found)
      if (found) bctfyp = temp_real
      
      call json%get('BC.bctfz', temp_real, found)
      if (found) bctfz = temp_real

      call json%get('BC.bctopm', temp_int, found)
      if (found) bctopm = temp_int

      call json%get('BC.bctopq', temp_int, found)
      if (found) bctopq = temp_int
      
      call json%get('BC.bctops', temp_int, found)
      if (found) bctops = temp_int

      call json%get('BC.bctopt', temp_int, found)
      if (found) bctopt = temp_int

      call json%get('BC.bcxm', temp_int, found)
      if (found) bcxm = temp_int
      
      call json%get('BC.bcxq', temp_int, found)
      if (found) bcxq = temp_int

      call json%get('BC.bcxs', temp_int, found)
      if (found) bcxs = temp_int

      call json%get('BC.bcxt', temp_int, found)
      if (found) bcxt = temp_int
      
      call json%get('BC.bcym', temp_int, found)
      if (found) bcym = temp_int

      call json%get('BC.bcyq', temp_int, found)
      if (found) bcyq = temp_int

      call json%get('BC.bcys', temp_int, found)
      if (found) bcys = temp_int
      
      call json%get('BC.bcyt', temp_int, found)
      if (found) bcyt = temp_int

      call json%get('BC.bczp', temp_int, found)
      if (found) bczp = temp_int

      call json%get('BC.ds', temp_real, found)
      if (found) ds = temp_real
      
      call json%get('BC.qt_top', temp_real, found)
      if (found) qt_top = temp_real

      call json%get('BC.qts', temp_real, found)
      if (found) qts = temp_real

      call json%get('BC.thl_top', temp_real, found)
      if (found) thl_top = temp_real
      
      call json%get('BC.thls', temp_real, found)
      if (found) thls = temp_real

      call json%get('BC.wqsurf', temp_real, found)
      if (found) wqsurf = temp_real

      call json%get('BC.wtsurf', temp_real, found)
      if (found) wtsurf = temp_real
      
      call json%get('BC.wttop', temp_real, found)
      if (found) wttop = temp_real

      call json%get('BC.z0', temp_real, found)
      if (found) z0 = temp_real

      call json%get('BC.z0h', temp_real, found)
      if (found) z0h = temp_real
      
      call json%get('BC.wsvsurfdum', temp_real, found)
      if (found) wsvsurfdum = temp_real

      call json%get('BC.wsvtopdum', temp_real, found)
      if (found) wsvtopdum = temp_real

   end subroutine read_bc_json

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
      call MPI_BCAST(wsvsurfdum, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wsvtopdum, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wtsurf, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(wttop, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(z0, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(z0h, 1, MY_REAL, 0, comm3d, mpierr)
   end subroutine broadcast_bc_parameters

   subroutine read_driver_json(json)
      !-----------------------------------------------------------------|
      ! Read DRIVER namelist parameters from JSON
      !-----------------------------------------------------------------|
      type(json_file), intent(inout) :: json
      logical :: found
      integer :: temp_int
      real :: temp_real
      logical :: temp_logical

      call json%get('DRIVER.chunkread_size', temp_int, found)
      if (found) chunkread_size = temp_int
      
      call json%get('DRIVER.driverjobnr', temp_int, found)
      if (found) driverjobnr = temp_int

      call json%get('DRIVER.driverstore', temp_logical, found)
      if (found) driverstore = temp_logical

      call json%get('DRIVER.dtdriver', temp_real, found)
      if (found) dtdriver = temp_real
      
      call json%get('DRIVER.iangledeg', temp_int, found)
      if (found) iangledeg = temp_int

      call json%get('DRIVER.idriver', temp_int, found)
      if (found) idriver = temp_int

      call json%get('DRIVER.iplane', temp_int, found)
      if (found) iplane = temp_int
      
      call json%get('DRIVER.lchunkread', temp_logical, found)
      if (found) lchunkread = temp_logical

      call json%get('DRIVER.tdriverstart', temp_real, found)
      if (found) tdriverstart = temp_real

   end subroutine read_driver_json

   subroutine broadcast_driver_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast DRIVER namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(idriver, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(chunkread_size, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(driverjobnr, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(driverstore, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(dtdriver, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(iangledeg, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(iplane, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(lchunkread, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(tdriverstart, 1, MY_REAL, 0, comm3d, mpierr)
   end subroutine broadcast_driver_parameters

   subroutine read_energybalance_json(json)
      !-----------------------------------------------------------------|
      ! Read ENERGYBALANCE namelist parameters from JSON
      !-----------------------------------------------------------------|
      type(json_file), intent(inout) :: json
      logical :: found
      integer :: temp_int
      real :: temp_real
      logical :: temp_logical

      call json%get('ENERGYBALANCE.bldt', temp_real, found)
      if (found) bldt = temp_real

      call json%get('ENERGYBALANCE.dteb', temp_real, found)
      if (found) dteb = temp_real
      
      call json%get('ENERGYBALANCE.flrt', temp_real, found)
      if (found) flrt = temp_real

      call json%get('ENERGYBALANCE.fraction', temp_real, found)
      if (found) fraction = temp_real

      call json%get('ENERGYBALANCE.grlai', temp_real, found)
      if (found) grlai = temp_real
      
      call json%get('ENERGYBALANCE.lconstw', temp_logical, found)
      if (found) lconstw = temp_logical

      call json%get('ENERGYBALANCE.leb', temp_logical, found)
      if (found) leb = temp_logical

      call json%get('ENERGYBALANCE.lfactlyrs', temp_logical, found)
      if (found) lfactlyrs = temp_logical
      
      call json%get('ENERGYBALANCE.lperiodicebcorr', temp_logical, found)
      if (found) lperiodicebcorr = temp_logical

      call json%get('ENERGYBALANCE.lvfsparse', temp_logical, found)
      if (found) lvfsparse = temp_logical

      call json%get('ENERGYBALANCE.lwriteebfiles', temp_logical, found)
      if (found) lwriteebfiles = temp_logical
      
      call json%get('ENERGYBALANCE.nfaclyrs', temp_int, found)
      if (found) nfaclyrs = temp_int

      call json%get('ENERGYBALANCE.nnz', temp_int, found)
      if (found) nnz = temp_int

      call json%get('ENERGYBALANCE.rsmin', temp_real, found)
      if (found) rsmin = temp_real
      
      call json%get('ENERGYBALANCE.sinkbase', temp_real, found)
      if (found) sinkbase = temp_real

      call json%get('ENERGYBALANCE.skylw', temp_real, found)
      if (found) skylw = temp_real

      call json%get('ENERGYBALANCE.wfc', temp_real, found)
      if (found) wfc = temp_real
      
      call json%get('ENERGYBALANCE.wgrmax', temp_real, found)
      if (found) wgrmax = temp_real

      call json%get('ENERGYBALANCE.wsoil', temp_real, found)
      if (found) wsoil = temp_real

      call json%get('ENERGYBALANCE.wwilt', temp_real, found)
      if (found) wwilt = temp_real

   end subroutine read_energybalance_json

   subroutine broadcast_energybalance_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast ENERGYBALANCE namelist parameters to all processes
      !-----------------------------------------------------------------|
      ! Add all ENERGYBALANCE variables with JSON support
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

   subroutine read_heatpump_json(json)
      !-----------------------------------------------------------------|
      ! Read HEATPUMP namelist parameters from JSON
      !-----------------------------------------------------------------|
      type(json_file), intent(inout) :: json
      logical :: found
      integer :: temp_int
      real :: temp_real
      logical :: temp_logical

      call json%get('HEATPUMP.lheatpump', temp_logical, found)
      if (found) lheatpump = temp_logical

      call json%get('HEATPUMP.lfan_hp', temp_logical, found)
      if (found) lfan_hp = temp_logical

      call json%get('HEATPUMP.nhppoints', temp_int, found)
      if (found) nhppoints = temp_int

      call json%get('HEATPUMP.Q_dot_hp', temp_real, found)
      if (found) Q_dot_hp = temp_real
      
      call json%get('HEATPUMP.QH_dot_hp', temp_real, found)
      if (found) QH_dot_hp = temp_real

   end subroutine read_heatpump_json

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

   subroutine read_inlet_json(json)
      !-----------------------------------------------------------------|
      ! Read INLET namelist parameters from JSON
      !-----------------------------------------------------------------|
      type(json_file), intent(inout) :: json
      logical :: found
      integer :: temp_int
      real :: temp_real
      logical :: temp_logical

      call json%get('INLET.di', temp_real, found)
      if (found) di = temp_real

      call json%get('INLET.dti', temp_real, found)
      if (found) dti = temp_real
      
      call json%get('INLET.inletav', temp_int, found)
      if (found) inletav = temp_int

      call json%get('INLET.lfixinlet', temp_logical, found)
      if (found) lfixinlet = temp_logical

      call json%get('INLET.lfixutauin', temp_logical, found)
      if (found) lfixutauin = temp_logical
      
      call json%get('INLET.linletra', temp_logical, found)
      if (found) linletra = temp_logical

      call json%get('INLET.lreadminl', temp_logical, found)
      if (found) lreadminl = temp_logical

      call json%get('INLET.lstoreplane', temp_logical, found)
      if (found) lstoreplane = temp_logical
      
      call json%get('INLET.lwallfunc', temp_logical, found)
      if (found) lwallfunc = temp_logical

      call json%get('INLET.uinf', temp_real, found)
      if (found) uinf = temp_real

      call json%get('INLET.vinf', temp_real, found)
      if (found) vinf = temp_real

   end subroutine read_inlet_json

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

   subroutine read_output_json(json)
      !-----------------------------------------------------------------|
      ! Read OUTPUT namelist parameters from JSON
      !-----------------------------------------------------------------|
      type(json_file), intent(inout) :: json
      logical :: found
      integer :: temp_int
      real :: temp_real
      logical :: temp_logical
      character(kind=CK,len=:), allocatable :: temp_string

      call json%get('OUTPUT.lfielddump', temp_logical, found)
      if (found) lfielddump = temp_logical
      
      call json%get('OUTPUT.tfielddump', temp_real, found)
      if (found) tfielddump = temp_real

      call json%get('OUTPUT.fieldvars', temp_string, found)
      if (found) fieldvars = temp_string
      
      call json%get('OUTPUT.tsample', temp_real, found)
      if (found) tsample = temp_real

      call json%get('OUTPUT.tstatsdump', temp_real, found)
      if (found) tstatsdump = temp_real

      call json%get('OUTPUT.tstatstart', temp_real, found)
      if (found) tstatstart = temp_real
      
      call json%get('OUTPUT.tcheck', temp_real, found)
      if (found) tcheck = temp_real

      call json%get('OUTPUT.islice', temp_int, found)
      if (found) islice = temp_int

      call json%get('OUTPUT.jslice', temp_int, found)
      if (found) jslice = temp_int

      call json%get('OUTPUT.kslice', temp_int, found)
      if (found) kslice = temp_int
      
      call json%get('OUTPUT.lislicedump', temp_logical, found)
      if (found) lislicedump = temp_logical

      call json%get('OUTPUT.ljslicedump', temp_logical, found)
      if (found) ljslicedump = temp_logical

      call json%get('OUTPUT.lkslicedump', temp_logical, found)
      if (found) lkslicedump = temp_logical

      call json%get('OUTPUT.slicevars', temp_string, found)
      if (found) slicevars = temp_string
      
      call json%get('OUTPUT.lmintdump', temp_logical, found)
      if (found) lmintdump = temp_logical

      call json%get('OUTPUT.ltdump', temp_logical, found)
      if (found) ltdump = temp_logical

      call json%get('OUTPUT.ltkedump', temp_logical, found)
      if (found) ltkedump = temp_logical
      
      call json%get('OUTPUT.lxydump', temp_logical, found)
      if (found) lxydump = temp_logical

      call json%get('OUTPUT.lxytdump', temp_logical, found)
      if (found) lxytdump = temp_logical

      call json%get('OUTPUT.lydump', temp_logical, found)
      if (found) lydump = temp_logical
      
      call json%get('OUTPUT.lytdump', temp_logical, found)
      if (found) lytdump = temp_logical

   end subroutine read_output_json

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

   subroutine read_purifs_json(json)
      !-----------------------------------------------------------------|
      ! Read PURIFS namelist parameters from JSON
      !-----------------------------------------------------------------|
      type(json_file), intent(inout) :: json
      logical :: found
      integer :: temp_int
      real :: temp_real
      logical :: temp_logical

      call json%get('PURIFS.lpurif', temp_logical, found)
      if (found) lpurif = temp_logical

      call json%get('PURIFS.npurif', temp_int, found)
      if (found) npurif = temp_int

      call json%get('PURIFS.epu', temp_real, found)
      if (found) epu = temp_real
      
      call json%get('PURIFS.qpu', temp_real, found)
      if (found) qpu = temp_real

   end subroutine read_purifs_json

   subroutine broadcast_purifs_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast PURIFS namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(lpurif, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(npurif, 1, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(epu, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(qpu, 1, MY_REAL, 0, comm3d, mpierr)
   end subroutine broadcast_purifs_parameters

   subroutine read_scalars_json(json)
      !-----------------------------------------------------------------|
      ! Read SCALARS namelist parameters from JSON
      !-----------------------------------------------------------------|
      type(json_file), intent(inout) :: json
      logical :: found
      integer :: temp_int
      logical :: temp_logical

      call json%get('SCALARS.lreadscal', temp_logical, found)
      if (found) lreadscal = temp_logical

      call json%get('SCALARS.lscasrc', temp_logical, found)
      if (found) lscasrc = temp_logical

      call json%get('SCALARS.lscasrcl', temp_logical, found)
      if (found) lscasrcl = temp_logical

      call json%get('SCALARS.lscasrcr', temp_logical, found)
      if (found) lscasrcr = temp_logical
      
      call json%get('SCALARS.nscasrc', temp_int, found)
      if (found) nscasrc = temp_int

      call json%get('SCALARS.nscasrcl', temp_int, found)
      if (found) nscasrcl = temp_int

      call json%get('SCALARS.nsv', temp_int, found)
      if (found) nsv = temp_int

   end subroutine read_scalars_json

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

   subroutine read_chemistry_json(json)
      !-----------------------------------------------------------------|
      ! Read CHEMISTRY namelist parameters from JSON
      !-----------------------------------------------------------------|
      use json_module
      type(json_file), intent(inout) :: json
      logical :: found
      real :: temp_real
      logical :: temp_logical

      call json%get('CHEMISTRY.lchem', temp_logical, found)
      if (found) lchem = temp_logical

      call json%get('CHEMISTRY.k1', temp_real, found)
      if (found) k1 = temp_real

      call json%get('CHEMISTRY.jno2', temp_real, found)
      if (found) JNO2 = temp_real

   end subroutine read_chemistry_json

   subroutine broadcast_chemistry_parameters()
      !-----------------------------------------------------------------|
      ! Broadcast CHEMISTRY namelist parameters to all processes
      !-----------------------------------------------------------------|
      call MPI_BCAST(lchem, 1, MPI_LOGICAL, 0, comm3d, mpierr)
      call MPI_BCAST(k1, 1, MY_REAL, 0, comm3d, mpierr)
      call MPI_BCAST(JNO2, 1, MY_REAL, 0, comm3d, mpierr)
   end subroutine broadcast_chemistry_parameters

   subroutine read_trees_json(json)
      !-----------------------------------------------------------------|
      ! Read TREES namelist parameters from JSON
      !-----------------------------------------------------------------|
      type(json_file), intent(inout) :: json
      logical :: found
      integer :: temp_int
      real :: temp_real
      logical :: temp_logical

      call json%get('TREES.ltrees', temp_logical, found)
      if (found) ltrees = temp_logical
      
      call json%get('TREES.ltreedump', temp_logical, found)
      if (found) ltreedump = temp_logical

      call json%get('TREES.ntrees', temp_int, found)
      if (found) ntrees = temp_int

      call json%get('TREES.cd', temp_real, found)
      if (found) cd = temp_real
      
      call json%get('TREES.dec', temp_real, found)
      if (found) dec = temp_real

      call json%get('TREES.dqdt', temp_real, found)
      if (found) dqdt = temp_real

      call json%get('TREES.lad', temp_real, found)
      if (found) lad = temp_real
      
      call json%get('TREES.lsize', temp_real, found)
      if (found) lsize = temp_real

      call json%get('TREES.qstar', temp_real, found)
      if (found) qstar = temp_real

      call json%get('TREES.r_s', temp_real, found)
      if (found) r_s = temp_real
      
      call json%get('TREES.ud', temp_real, found)
      if (found) ud = temp_real

   end subroutine read_trees_json

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

   subroutine read_walls_json(json)
      !-----------------------------------------------------------------|
      ! Read WALLS namelist parameters from JSON
      !-----------------------------------------------------------------|
      type(json_file), intent(inout) :: json
      logical :: found
      integer :: temp_int
      real :: temp_real
      logical :: temp_logical

      call json%get('WALLS.dtfac', temp_real, found)
      if (found) dtfac = temp_real

      call json%get('WALLS.fkar', temp_real, found)
      if (found) fkar = temp_real
      
      call json%get('WALLS.iwallmoist', temp_int, found)
      if (found) iwallmoist = temp_int

      call json%get('WALLS.iwallmom', temp_int, found)
      if (found) iwallmom = temp_int

      call json%get('WALLS.iwallscal', temp_int, found)
      if (found) iwallscal = temp_int

      call json%get('WALLS.iwalltemp', temp_int, found)
      if (found) iwalltemp = temp_int
      
      call json%get('WALLS.lbottom', temp_logical, found)
      if (found) lbottom = temp_logical

      call json%get('WALLS.lnorec', temp_logical, found)
      if (found) lnorec = temp_logical

      call json%get('WALLS.lwritefac', temp_logical, found)
      if (found) lwritefac = temp_logical
      
      call json%get('WALLS.nblocks', temp_int, found)
      if (found) nblocks = temp_int

      call json%get('WALLS.nbndpts_c', temp_int, found)
      if (found) nbndpts_c = temp_int

      call json%get('WALLS.nbndpts_u', temp_int, found)
      if (found) nbndpts_u = temp_int
      
      call json%get('WALLS.nbndpts_v', temp_int, found)
      if (found) nbndpts_v = temp_int

      call json%get('WALLS.nbndpts_w', temp_int, found)
      if (found) nbndpts_w = temp_int

      call json%get('WALLS.nfcts', temp_int, found)
      if (found) nfcts = temp_int
      
      call json%get('WALLS.nfctsecs_c', temp_int, found)
      if (found) nfctsecs_c = temp_int

      call json%get('WALLS.nfctsecs_u', temp_int, found)
      if (found) nfctsecs_u = temp_int

      call json%get('WALLS.nfctsecs_v', temp_int, found)
      if (found) nfctsecs_v = temp_int
      
      call json%get('WALLS.nfctsecs_w', temp_int, found)
      if (found) nfctsecs_w = temp_int

      call json%get('WALLS.nsolpts_c', temp_int, found)
      if (found) nsolpts_c = temp_int

      call json%get('WALLS.nsolpts_u', temp_int, found)
      if (found) nsolpts_u = temp_int
      
      call json%get('WALLS.nsolpts_v', temp_int, found)
      if (found) nsolpts_v = temp_int

      call json%get('WALLS.nsolpts_w', temp_int, found)
      if (found) nsolpts_w = temp_int

      call json%get('WALLS.prandtlturb', temp_real, found)
      if (found) prandtlturb = temp_real

   end subroutine read_walls_json

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

   subroutine read_info_json(json)
      !-----------------------------------------------------------------|
      ! Read INFO namelist parameters from JSON
      !-----------------------------------------------------------------|
      use modinletdata, only: nprocsinl, jgtotinl, kmaxin, dtin, wtop, totalreadu
      
      type(json_file), intent(inout) :: json
      logical :: found
      integer :: temp_int
      real :: temp_real
      character(len=100) :: temp_string

      call json%get('INFO.nprocsinl', temp_int, found)
      if (found) nprocsinl = temp_int

      call json%get('INFO.jgtotinl', temp_int, found)
      if (found) jgtotinl = temp_int

      call json%get('INFO.kmaxin', temp_int, found)
      if (found) kmaxin = temp_int

      call json%get('INFO.dtin', temp_real, found)
      if (found) dtin = temp_real
      
      call json%get('INFO.wtop', temp_real, found)
      if (found) wtop = temp_real

      call json%get('INFO.totalreadu', temp_real, found)
      if (found) totalreadu = temp_real

      ! Note: namezinlet and zgrid are file-specific parameters
      ! that don't translate to JSON configuration

   end subroutine read_info_json

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

   subroutine read_namsubgrid_json(json)
      !-----------------------------------------------------------------|
      ! Read NAMSUBGRID namelist parameters from JSON
      !-----------------------------------------------------------------|
      use modsubgrid, only: ldelta, lmason, cf, cn, Rigc, Prandtl, &
                                lsmagorinsky, lvreman, loneeqn, c_vreman, &
                                cs, nmason, lbuoycorr 

      type(json_file), intent(inout) :: json
      logical :: found
      logical :: temp_logical
      real :: temp_real

      call json%get('NAMSUBGRID.ldelta', temp_logical, found)
      if (found) ldelta = temp_logical

      call json%get('NAMSUBGRID.lmason', temp_logical, found)
      if (found) lmason = temp_logical

      call json%get('NAMSUBGRID.lsmagorinsky', temp_logical, found)
      if (found) lsmagorinsky = temp_logical

      call json%get('NAMSUBGRID.lvreman', temp_logical, found)
      if (found) lvreman = temp_logical

      call json%get('NAMSUBGRID.loneeqn', temp_logical, found)
      if (found) loneeqn = temp_logical

      call json%get('NAMSUBGRID.lbuoycorr', temp_logical, found)
      if (found) lbuoycorr = temp_logical

      call json%get('NAMSUBGRID.cf', temp_real, found)
      if (found) cf = temp_real

      call json%get('NAMSUBGRID.cn', temp_real, found)
      if (found) cn = temp_real

      call json%get('NAMSUBGRID.Rigc', temp_real, found)
      if (found) Rigc = temp_real

      call json%get('NAMSUBGRID.Prandtl', temp_real, found)
      if (found) Prandtl = temp_real

      call json%get('NAMSUBGRID.c_vreman', temp_real, found)
      if (found) c_vreman = temp_real

      call json%get('NAMSUBGRID.cs', temp_real, found)
      if (found) cs = temp_real

      call json%get('NAMSUBGRID.nmason', temp_real, found)
      if (found) nmason = temp_real

   end subroutine read_namsubgrid_json

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
end module readparameters