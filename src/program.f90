!> \file program.f90
!! Main program

!! \section License License
!!  This file is part of DALES.
!!
!!  DALES is free software; you can redistribute it and/or modify it under the
!! terms of the GNU General Public License as published by the Free Software
!! Foundation; either version 3 of the License, or (at your option) any later
!! version.
!!
!!  DALES is distributed in the hope that it will be useful, but WITHOUT ANY
!! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License along with
!! this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!!  Copyright 1993-2009 Delft University of Technology, Wageningen University,
!! Utrecht University, KNMI
!!
program uDALES 

!!----------------------------------------------------------------
!!     0.0    USE STATEMENTS FOR CORE MODULES
!!----------------------------------------------------------------
  use modmpi,            only : initmpi,exitmpi,myid,starttimer
  use modglobal,         only : initglobal,rk3step,timeleft
  use modglobal,         only : runmode,RUN_COLDSTART,RUN_WARMSTART,RUN_DRIVER,RUN_STRATSTART,TEST_ROUNDTRIP,TEST_IO,TEST_SPARSE_IJK,TEST_2DCOMP_INIT_EXIT,TEST_MPI_OPERATORS
  use modstartup,        only : readconfig,init2decomp,checkinitvalues,readinitfiles,exitmodules
  use modfields,         only : initfields
  use modsave,           only : writerestartfiles
  use modboundary,       only : initboundary,boundary,grwdamp,halos
  use modthermodynamics, only : initthermodynamics,thermodynamics
  use modsubgrid,        only : initsubgrid,subgrid
  use modforces,         only : calcfluidvolumes,forces,coriolis,lstend,fixuinf1,fixuinf2,fixthetainf,nudge,masscorr,shiftedPBCs,periodicEBcorr
  use modpois,           only : initpois,poisson
  use modibm,            only : initibm,createmasks,ibmwallfun,ibmnorm,bottom
  use vegetation,        only : init_vegetation, vegetation_forcing
  use modpurifiers,      only : createpurifiers,purifiers
  use modheatpump,       only : init_heatpump,heatpump,exit_heatpump
  use initfac,           only : readfacetfiles
  use modEB,             only : initEB,EB
  use moddriver,         only : initdriver
  use modchecksim,       only : initchecksim,checksim
  use modtimedep,        only : inittimedep,timedep

!------------------------------------------------------------------------------
!     0.1     USE STATEMENTS FOR STATISTICAL AND INSTANTANEOUS OUTPUT ROUTINES
!------------------------------------------------------------------------------
  use modstatsdump,      only : initstatsdump,statsdump,exitstatsdump    !tg3315
  use stats,             only : stats_init,stats_main,stats_exit
  use instant,           only : instant_init,instant_main,instant_exit
  
!------------------------------------------------------------------------------
!     0.2     USE STATEMENTS FOR TESTS
!------------------------------------------------------------------------------
  use tests,             only : tests_read_sparse_ijk,tests_2decomp_init_exit,tests_mpi_operators
  use tests,             only : init_tests,tests_roundtrip,exit_tests

  implicit none

!----------------------------------------------------------------
!     0      READ NAMELISTS,INITIALISE GRID, CONSTANTS AND FIELDS
!----------------------------------------------------------------
  call initmpi

  !call startup
  call readconfig

  ! TEST_ROUNDTRIP is dispatched before init2decomp because it manages its
  ! own decomposition setup and teardown in init_tests/exit_tests.
  call execute_early_runmode_actions

  call init2decomp

  call checkinitvalues

  call initglobal

  ! Execute tests if needed
  call execute_runmode_actions

  call initfields

  call initboundary

  call initthermodynamics

  call initsubgrid

  ! call initinlet

  call initdriver

  call initpois

  call readfacetfiles
  ! These should be combined once file format is sorted
  call initibm

  call createmasks

  call calcfluidvolumes

  call readinitfiles

  call createscals

!---------------------------------------------------------
!      2     INITIALIZE STATISTICAL ROUTINES AND ADD-ONS
!---------------------------------------------------------
  call initchecksim ! Could be deprecated

  call initstatsdump
  call stats_init
  call instant_init

  call initEB

  call inittimedep

  call boundary

  call init_vegetation

  call createpurifiers

  call init_heatpump

!------------------------------------------------------
!   3.0   MAIN TIME LOOP
!------------------------------------------------------
  !write(*,*) 'Starting rank ', myid
  call starttimer
  do while ((timeleft>0) .or. (rk3step < 3))

    call tstep_update

    call timedep

!-----------------------------------------------------
!   3.2   ADVECTION AND DIFFUSION
!-----------------------------------------------------

    call advection ! includes predicted pressure gradient term

    call shiftedPBCs

    call subgrid

!-----------------------------------------------------
!   3.3   THE SURFACE LAYER
!-----------------------------------------------------

    call bottom

!-----------------------------------------------------
!   3.4   REMAINING TERMS
!-----------------------------------------------------

    call coriolis       !remaining terms of ns equation

    call forces         !remaining terms of ns equation

    call lstend         !large scale forcings

    call nudge          ! nudge top cells of fields to enforce steady-state

    call ibmwallfun     ! immersed boundary forcing: only shear forces.
    call periodicEBcorr

    call masscorr       ! correct pred. velocity pup to get correct mass flow

    call ibmnorm        ! immersed boundary forcing: set normal velocities to zero

    call EB

    call vegetation_forcing

    call heatpump

    call scalsource     ! adds continuous forces in specified region of domain

!------------------------------------------------------
!   3.4   EXECUTE ADD ONS
!------------------------------------------------------
    call fixuinf2

    call fixuinf1

!-----------------------------------------------------------------------
!   3.5  PRESSURE FLUCTUATIONS, TIME INTEGRATION AND BOUNDARY CONDITIONS
!-----------------------------------------------------------------------
    call grwdamp        !damping at top of the model

    call poisson

    call purifiers      !placing of purifiers here may need to be checked

    call tstep_integrate

    call halos

    call checksim

    call statsdump     ! will depricate soon; contains tke budget only(not working)
    call stats_main
    call instant_main

    call boundary

    !call fixthetainf ! deprecated

!-----------------------------------------------------
!   3.6   LIQUID WATER CONTENT AND DIAGNOSTIC FIELDS
!-----------------------------------------------------
    call thermodynamics

!-----------------------------------------------------
!   3.7  WRITE RESTARTFILES AND DO STATISTICS
!------------------------------------------------------

    call writerestartfiles

  end do
!-------------------------------------------------------
!             END OF TIME LOOP
!-------------------------------------------------------

!--------------------------------------------------------
!    4    FINALIZE ADD ONS AND THE MAIN PROGRAM
!-------------------------------------------------------
  call exitstatsdump     !tg3315
  call exit_heatpump
  call stats_exit
  call instant_exit
  !call exitmodules
  call exitmpi

contains
  subroutine execute_early_runmode_actions
    select case (runmode)
      case (TEST_ROUNDTRIP)
        call init_tests
        call tests_roundtrip
        call exit_tests
        stop
      case default
        return
    end select
  end subroutine execute_early_runmode_actions

  subroutine execute_runmode_actions
    logical :: test_failed
    logical :: invalid_runmode

    test_failed = .false.
    invalid_runmode = .false.
    select case (runmode)
      case (RUN_COLDSTART, RUN_WARMSTART, RUN_DRIVER, RUN_STRATSTART)
        return
        ! Normal execution mode, do nothing special here
      case (TEST_SPARSE_IJK)
        ! Execute tests for reading sparse arrays
        test_failed = .not. tests_read_sparse_ijk()
      case (TEST_MPI_OPERATORS)
        test_failed = .not. tests_mpi_operators()
      case (TEST_2DCOMP_INIT_EXIT)
        call tests_2decomp_init_exit
      case (TEST_IO)
        write(*,*) 'TEST_IO mode not yet implemented'
        invalid_runmode = .true.
      case default
        write(*,*) 'Unknown runmode:', runmode
        invalid_runmode = .true.
    end select

    call exitmpi

    if (invalid_runmode) then
      stop 1
    end if

    ! Return appropriate exit code for unit tests:
    ! 0 = success, 1 = failure
    if (test_failed) then
      stop 1
    else
      stop 0
    end if

  end subroutine execute_runmode_actions

end program uDALES
