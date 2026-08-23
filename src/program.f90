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
  use mpi
  use modmpi,            only : initmpi,exitmpi,starttimer,myid
#if defined(_GPU)
  use cudafor
  use modcuda,           only : initCUDA, updateDevice, &
                                updateHostForFielddump, updateHostForStatsdump, &
                                updateHostForUnportedRoutines, integrateFacFluxDevice, &
                                updateFacIntegralsHost, checkCUDA, exitCUDA
#endif
#if defined(_GPU) && defined(UDALES_DEBUG)
  use modcuda,           only : assertHostMatchesDevice, enableRoundTripCheck
  use tests_cuda,        only : run_cuda_selftests_if_requested
#endif
  use modglobal,         only : initglobal,rk3step,timeleft
#if defined(_GPU) && defined(UDALES_DEBUG)
  use modglobal,         only : ladaptive
#endif
  use modglobal,         only : runmode,RUN_COLDSTART,RUN_WARMSTART,RUN_DRIVER,RUN_STRATSTART,TEST_SPARSE_IJK,TEST_2DCOMP_INIT_EXIT, &
                                TEST_MPI_OPERATORS,TEST_IBM_CELL_LOOKUP,TEST_NUDGE,TEST_IBM_WALLFUN, &
                                TEST_PERIODIC_EBCORR,TEST_MASSCORR,TEST_IBMNORM,TEST_EB, &
                                TEST_VEGETATION,TEST_CHECKSIM,TEST_DRIVER_PLANES
  use modstartup,        only : readnamelists,init2decomp,checkinitvalues,readinitfiles,exitmodules
  use modfields,         only : initfields
  use modsave,           only : writerestartfiles
  use modboundary,       only : initboundary,boundary,boundary_conditions,grwdamp,exchange_halos
  use modthermodynamics, only : initthermodynamics,thermodynamics
  use modsubgrid,        only : initsubgrid,subgrid
  use modforces,         only : calcfluidvolumes,forces,coriolis,lstend,fixuinf1,fixuinf2,nudge,masscorr,shiftedPBCs,periodicEBcorr
  use modpois,           only : initpois,poisson
  use modibm,            only : initibm,createmasks,ibmwallfun,ibmnorm,bottom
  use vegetation,        only : init_vegetation, vegetation_forcing
#if defined(_GPU)
  use vegetation,        only : updateVegDiagHost
#endif
  ! Purifier forcing is not yet supported by the GPU execution path.
  ! use modpurifiers,      only : createpurifiers,purifiers
  use modheatpump,       only : init_heatpump,heatpump,exit_heatpump
  use initfac,           only : readfacetfiles
  use modEB,             only : initEB,EB,eb_will_run
  use moddriver,         only : initdriver
  use modadvection,      only : advection
  use modtstep,          only : tstep_update,tstep_integrate
  use modscalsource,     only : createscals,scalsource,exitscals

!----------------------------------------------------------------
!     0.1     USE STATEMENTS FOR ADDONS STATISTICAL ROUTINES
!----------------------------------------------------------------
  use modchecksim,     only : initchecksim,checksim
  use modstat_nc,      only : initstat_nc
  use modfielddump,    only : initfielddump,fielddump,exitfielddump
  use modstatsdump,    only : initstatsdump,statsdump,exitstatsdump
#if defined(_GPU)
  use modfielddump,    only : fielddump_will_sample
  use modstatsdump,    only : statsdump_will_sample
#endif
  use modtimedep,      only : inittimedep,timedep
  use tests,           only : tests_read_sparse_ijk,tests_2decomp_init_exit,tests_mpi_operators,tests_ibm_cell_lookup,tests_nudge,tests_ibm_wallfun, &
                            tests_periodic_ebcorr,tests_masscorr,tests_ibmnorm,tests_eb, &
                            tests_vegetation,tests_checksim,tests_driver_planes
  implicit none

  real    :: stime

!----------------------------------------------------------------
!     0      READ NAMELISTS,INITIALISE GRID, CONSTANTS AND FIELDS
!----------------------------------------------------------------
  call initmpi

  !call startup
  call readnamelists

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

  call initstat_nc ! Could be deprecated

  call initstatsdump

  call initEB

  call inittimedep

  call initfielddump

  call boundary

  call init_vegetation

  ! call createpurifiers

  call init_heatpump

  !call fielddump

#if defined(_GPU)
  call initCUDA
#endif
#if defined(_GPU) && defined(UDALES_DEBUG)
  call run_cuda_selftests_if_requested
#endif

!------------------------------------------------------
!   3.0   MAIN TIME LOOP
!------------------------------------------------------
#if defined(_GPU) && defined(UDALES_DEBUG)
  ! From here on, every field updateDevice uploads has to have been brought
  ! down during the step. Armed now rather than at initialisation, because the
  ! device self-tests call updateDevice with nothing pulled.
  call enableRoundTripCheck
#endif

  call starttimer
  do while ((timeleft>0) .or. (rk3step < 3))

    ! Any routine added to this loop must have its GPU data transfers addressed:
    ! nothing between updateDevice and the post-Poisson handover may read host
    ! fields.

    stime = MPI_Wtime()

#if defined(_GPU) && defined(UDALES_DEBUG)
    ! tstep_update is the one host routine left above updateDevice that reads
    ! fields rather than clocks, and with an adaptive timestep what it reads
    ! decides dt - so a stale copy here does not corrupt an output, it moves
    ! the whole solution. ladaptive is the condition because that is when it
    ! reads them at all; which fields it reads is stated in modcuda, so the
    ! two lists can disagree and be caught.
    if (ladaptive) call assertHostMatchesDevice('tstep_update')
#endif

    call tstep_update
    call print_time('tstep_update')

    call timedep
    call print_time('timedep')

!-----------------------------------------------------
!   3.2   ADVECTION AND DIFFUSION
!-----------------------------------------------------
#if defined(_GPU)
    call updateDevice
#endif
    call print_time('updateDevice')

    call advection ! includes predicted pressure gradient term
    call print_time('advection')

    call shiftedPBCs
    call print_time('shiftedPBCs')

    call subgrid
    call print_time('subgrid')

!-----------------------------------------------------
!   3.3   THE SURFACE LAYER
!-----------------------------------------------------

    call bottom
    call print_time('bottom')

!-----------------------------------------------------
!   3.4   REMAINING TERMS
!-----------------------------------------------------

    call coriolis       !remaining terms of ns equation
    call print_time('coriolis')

    call forces         !remaining terms of ns equation
    call print_time('forces')

    call lstend         !large scale forcings
    call print_time('lstend')

    call nudge          ! nudge top cells of fields to enforce steady-state
    call print_time('nudge')

    call ibmwallfun     ! immersed boundary forcing: only shear forces.
    call print_time('ibmwallfun')

    call periodicEBcorr
    call print_time('periodicEBcorr')

    call masscorr       ! correct pred. velocity pup to get correct mass flow
    call print_time('masscorr')

    call ibmnorm        ! immersed boundary forcing: set normal velocities to zero
    call print_time('ibmnorm')

#if defined(_GPU)
    ! The facet flux accumulators the wall functions filled on the device.
    ! Integrating them costs no traffic at all: the time integral is per-rank,
    ! so it accumulates on the device across the hundreds of steps between
    ! energy balances and only the total ever comes down.
    call integrateFacFluxDevice

    ! And that total comes down only on the steps where the balance fires.
    ! Both this and EB's own guard read eb_will_run, so they cannot drift: if
    ! the copy were skipped on a firing step the balance would silently run on
    ! zero facet heat flux, on the GPU and nowhere else, which is what the
    ! facEB comparison in the surface-energy-balance parity case catches.
    if (eb_will_run()) call updateFacIntegralsHost
#endif
    call print_time('Additional EB traffic')

    call EB
    call print_time('EB')

    call vegetation_forcing
    call print_time('vegetation_forcing')

    call heatpump
    call print_time('heatpump')

    call scalsource     ! adds continuous forces in specified region of domain
    call print_time('scalsource')

!------------------------------------------------------
!   3.4   EXECUTE ADD ONS
!------------------------------------------------------
    call fixuinf2
    call print_time('fixuinf2')
    call fixuinf1
    call print_time('fixuinf1')

!-----------------------------------------------------------------------
!   3.5  PRESSURE FLUCTUATIONS, TIME INTEGRATION AND BOUNDARY CONDITIONS
!-----------------------------------------------------------------------
    call grwdamp        !damping at top of the model
    call print_time('grwdamp')

    call poisson
    call print_time('poisson')

    ! call purifiers      !placing need to be checked; Not GPU compatible yet

    call tstep_integrate
    call print_time('tstep_integrate')

    call exchange_halos
    call print_time('exchange_halos')

    call checksim
    call print_time('checksim')

#if defined(_GPU)
    if (fielddump_will_sample()) call updateHostForFielddump
#endif
    call print_time('updateHostForFielddump')

#if defined(_GPU) && defined(UDALES_DEBUG)
    ! Every host field fielddump can read must still be the one the device
    ! holds. The transfers above are conditional, and the failure mode of a
    ! wrong condition is a dump quietly written from the previous step, so the
    ! check is here rather than left to a parity comparison that only some
    ! configurations would notice.
    if (fielddump_will_sample()) call assertHostMatchesDevice('fielddump')
#endif

    call fielddump
    call print_time('fielddump')

#if defined(_GPU)
    if (statsdump_will_sample()) then
       call updateHostForStatsdump
       call updateVegDiagHost
    end if
#endif
    call print_time('updateHostForStatsdump')

#if defined(_GPU) && defined(UDALES_DEBUG)
    if (statsdump_will_sample()) call assertHostMatchesDevice('statsdump')
#endif

    call statsdump
    call print_time('statsdump')

    ! Above the handover, not below it, now that it runs on the device: the
    ! boundary planes go on before the host copies are taken, rather than being
    ! written into host memory afterwards and carried back up by the next
    ! updateDevice. It ends by clearing the pull bitmap, because fielddump and
    ! statsdump have already fetched some of these fields further up.
    call boundary_conditions
    call print_time('boundary_conditions')

#if defined(_GPU)
    ! And what the host routines still left in the loop need, which is every
    ! prognostic field, on every stage. thermodynamics reads and derives from
    ! them, and the next updateDevice uploads the result - so this one is
    ! correctness, not output.
    call updateHostForUnportedRoutines
#endif
    call print_time('updateHostForUnportedRoutines')

!-----------------------------------------------------
!   3.6   LIQUID WATER CONTENT AND DIAGNOSTIC FIELDS
!-----------------------------------------------------
    call thermodynamics
    call print_time('thermodynamics')

!-----------------------------------------------------
!   3.7  WRITE RESTARTFILES AND DO STATISTICS
!------------------------------------------------------

    call writerestartfiles
    call print_time('writerestartfiles')
  end do
!-------------------------------------------------------
!             END OF TIME LOOP
!-------------------------------------------------------

!--------------------------------------------------------
!    4    FINALIZE ADD ONS AND THE MAIN PROGRAM
!-------------------------------------------------------
  call exitfielddump
  call exitstatsdump     !tg3315
  call exit_heatpump
  !call exitmodules
  !call exittest

  call exitscals

#if defined(_GPU)
  call exitCUDA
#endif

  call exitmpi

contains
  subroutine print_time(routine_name)
    implicit none
    character(len=*), intent(in) :: routine_name
#if defined(_GPU)
    call checkCUDA( cudaDeviceSynchronize(), 'cudaDeviceSynchronize in program print_time for ' // trim(routine_name) )
#endif
    write(6,'(A,I0,3A,F10.6,A)')'Rank:',myid,': Time taken by ', trim(routine_name), ' : ', MPI_Wtime() - stime, ' seconds'

    stime = MPI_Wtime()
  end subroutine print_time

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
      case (TEST_IBM_CELL_LOOKUP)
        test_failed = .not. tests_ibm_cell_lookup()
      case (TEST_NUDGE)
        test_failed = .not. tests_nudge()
      case (TEST_IBM_WALLFUN)
        test_failed = .not. tests_ibm_wallfun()
      case (TEST_PERIODIC_EBCORR)
        test_failed = .not. tests_periodic_ebcorr()
      case (TEST_MASSCORR)
        test_failed = .not. tests_masscorr()
      case (TEST_IBMNORM)
        test_failed = .not. tests_ibmnorm()
      case (TEST_EB)
        test_failed = .not. tests_eb()
      case (TEST_VEGETATION)
        test_failed = .not. tests_vegetation()
      case (TEST_CHECKSIM)
        test_failed = .not. tests_checksim()
      case (TEST_DRIVER_PLANES)
        test_failed = .not. tests_driver_planes()
      case (TEST_2DCOMP_INIT_EXIT)
        call tests_2decomp_init_exit
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
