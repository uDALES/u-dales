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
  use modmpi,            only : initmpi,exitmpi,starttimer
#if defined(_GPU)
  use cudafor
  use modcuda,           only : initCUDA, updateDevice, updateHost, &
                                updateHostAfterPoiss, integrateFacFluxDevice, &
                                updateFacIntegralsHost, checkCUDA, exitCUDA
#endif
#if defined(_GPU) && defined(UDALES_DEBUG)
  use tests_cuda,        only : run_cuda_selftests_if_requested
#endif
  use modglobal,         only : initglobal,rk3step,timeleft
  use modglobal,         only : runmode,RUN_COLDSTART,RUN_WARMSTART,RUN_DRIVER,RUN_STRATSTART,TEST_SPARSE_IJK,TEST_2DCOMP_INIT_EXIT, &
                                TEST_MPI_OPERATORS,TEST_IBM_CELL_LOOKUP,TEST_NUDGE,TEST_IBM_WALLFUN, &
                                TEST_PERIODIC_EBCORR,TEST_MASSCORR,TEST_IBMNORM,TEST_EB, &
                                TEST_VEGETATION
  use modstartup,        only : readnamelists,init2decomp,checkinitvalues,readinitfiles,exitmodules
  use modfields,         only : initfields
  use modsave,           only : writerestartfiles
  use modboundary,       only : initboundary,boundary,grwdamp
#if defined(_GPU)
  use modboundary,       only : halos_device
#else
  use modboundary,       only : halos
#endif
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
  use modstatsdump,    only : initstatsdump,statsdump,exitstatsdump    !tg3315
#if defined(_GPU)
  use modstatsdump,    only : statsdump_will_sample
#endif
  use modtimedep,      only : inittimedep,timedep
  use tests,           only : tests_read_sparse_ijk,tests_2decomp_init_exit,tests_mpi_operators,tests_ibm_cell_lookup,tests_nudge,tests_ibm_wallfun, &
                            tests_periodic_ebcorr,tests_masscorr,tests_ibmnorm,tests_eb, &
                            tests_vegetation
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
  call starttimer
  do while ((timeleft>0) .or. (rk3step < 3))

    call tstep_update

    call timedep

!-----------------------------------------------------
!   3.2   ADVECTION AND DIFFUSION
!-----------------------------------------------------
    stime = MPI_Wtime()
#if defined(_GPU)
    call updateDevice
#endif
    write(6,'(A,F10.6)')'updateDevice time = ', MPI_Wtime() - stime

    stime = MPI_Wtime()

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

#if defined(_GPU)
    call checkCUDA( cudaDeviceSynchronize(), 'cudaDeviceSynchronize in program' )
#endif
    write(6,'(A,F10.6)')'(advection to ibmnorm) time = ', MPI_Wtime() - stime

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
    !
    ! Neither of these depends on updateHost, so EB survives unchanged when
    ! vegetation_forcing is ported and updateHost goes away.
    if (eb_will_run()) call updateFacIntegralsHost
#endif

    call EB

    ! Ahead of updateHost now, not behind it. On a GPU build every loop in here
    ! runs on the device against the mirrors ibmnorm has just pinned, so the
    ! previous-step fields no longer have to come down for it - which is what
    ! the ltrees block updateHost used to carry. On a CPU build nothing moved:
    ! updateHost does not exist there, and the order against EB is unchanged.
    call vegetation_forcing

    stime = MPI_Wtime()
#if defined(_GPU)
    call updateHost
#endif
    write(6,'(A,F10.6)')'updateHost time = ', MPI_Wtime() - stime

    ! updateDevicePriorPoiss used to sit here, re-uploading what updateHost had
    ! just brought down plus the fields updateDevice had already sent at the
    ! top of the stage. Once vegetation_forcing moved onto the device there was
    ! no host work left between the two to justify it: nothing on the host
    ! writes any of those fields in this window, and nothing on the device
    ! writes u0/v0/w0/pres0/thl0/thl0c/qt0/sv0 either, so both halves of the
    ! round trip carried the value they already held. The inlet profiles and
    ! the driver plane were the exception - that routine was their only upload
    ! - and they are copied by updateDevice now.

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

    stime = MPI_Wtime()
    call poisson
#if defined(_GPU)
    call checkCUDA( cudaDeviceSynchronize(), 'cudaDeviceSynchronize in program' )
#endif
    write(6,'(A,F10.6)')'poisson time = ', MPI_Wtime() - stime

    ! call purifiers      !placing need to be checked; Not GPU compatible yet

    call tstep_integrate

#if defined(_GPU)
    call halos_device
#else
    call halos
#endif

#if defined(_GPU)
    call updateHostAfterPoiss
#endif

    call checksim

    call fielddump

#if defined(_GPU)
    ! The tree tendencies stay on the device between samples. statsdump reads
    ! them under ltreedump only, and only on the steps it samples, so ask it
    ! rather than draining every step - the same shape as the facet integrals
    ! above. updateVegDiagHost is itself a no-op without ltreedump, so the
    ! predicate never has to know about trees.
    if (statsdump_will_sample()) call updateVegDiagHost
#endif

    call statsdump

    call boundary

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
