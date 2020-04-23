!> \file program.f90
!! Main program

!>
!! \mainpage
!! Dutch Atmospheric Large Eddy Simulation -URBAN
!! \section DALES Dutch Atmospheric Large Eddy Simulation -URBAN
!!
!! @version 48
!!
!
!! \todo
!!
!! \section License License
!!  This file is part of DALESURBAN.
!!  Copyright 1993-2014 Delft University of Technology
!!
program DALESURBAN      !Version 48

!!----------------------------------------------------------------
!!     0.0    USE STATEMENTS FOR CORE MODULES
!!----------------------------------------------------------------
  use modmpi,            only : myid, initmpi
  use modglobal,         only : rk3step,timeleft,ib,jb,kb,ke
  use modstartup,        only : startup,exitmodules
  use modsave,           only : writerestartfiles
  use modboundary,       only : boundary, grwdamp,tqaver
  use modthermodynamics, only : thermodynamics
!  use modsurface,        only : surface
  use modsubgrid,        only : subgrid
  use modforces,         only : forces,coriolis,lstend,fixuinf1,fixuinf2,fixthetainf,nudge, masscorr
  use modpois,           only : poisson
  use modibm,            only : createwalls,ibmwallfun,ibmnorm,nearwall,bottom
  use initfac,           only : readfacetfiles
  use modEB,             only : initEB,EB

!----------------------------------------------------------------
!     0.1     USE STATEMENTS FOR ADDONS STATISTICAL ROUTINES
!----------------------------------------------------------------
  use modchecksim,     only : initchecksim, checksim
  use modstat_nc,      only : initstat_nc
  use modfielddump,    only : initfielddump, fielddump,exitfielddump
  use modstatsdump,    only : initstatsdump,statsdump,exitstatsdump    !tg3315
  !use modbudget,       only : initbudget, budgetstat, exitbudget
  implicit none


!----------------------------------------------------------------
!     1      READ NAMELISTS,INITIALISE GRID, CONSTANTS AND FIELDS
!----------------------------------------------------------------
  call initmpi
  write(*,*) "done initmpi"
  call startup
  write(*,*) "done startup"

!---------------------------------------------------------
!      2     INITIALIZE STATISTICAL ROUTINES AND ADD-ONS
!---------------------------------------------------------
  call initchecksim
  call initstat_nc

  call initfielddump
  call initstatsdump !tg3315

  call readfacetfiles
  call initEB
  write(*,*) "done init stuff"

  write(6,*) 'Determine immersed walls'
  call createwalls    ! determine walls/blocks
 ! call nearwall       ! determine minimum distance and corresponding shear components, ils13 10.07.17, commented, not functional at the moment, not needed for vreman but for smag., fix in modibm
  write(6,*) 'Finished determining immersed walls'

  call boundary  !ils13 22.06.2017 inserted boundary here to get values at ghost cells before iteration starts

!  not necessary but abates the fact that temp field is randomised by randomisation of just velocity fields
!  (because advection at start of time loop without being divergence free)
!  call poisson

!------------------------------------------------------
!   3.0   MAIN TIME LOOP
!------------------------------------------------------
  write(*,*)'START myid ', myid
  do while ((timeleft>0) .or. (rk3step < 3))
    call tstep_update

!-----------------------------------------------------
!   3.2   ADVECTION AND DIFFUSION
!-----------------------------------------------------
  
    call advection                ! now also includes predicted pressure gradient term  

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

    call masscorr       ! correct pred. velocity pup to get correct mass flow
                                                                                         
    call ibmnorm        ! immersed boundary forcing: set normal velocities to zero  

    call EB

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

    call tstep_integrate

    call boundary

    call fixthetainf

!-----------------------------------------------------
!   3.6   LIQUID WATER CONTENT AND DIAGNOSTIC FIELDS
!-----------------------------------------------------
    call thermodynamics

!-----------------------------------------------------
!   3.7  WRITE RESTARTFILES AND DO STATISTICS
!------------------------------------------------------

    call checksim
   ! call writedatafiles   ! write data files for later analysis
    call writerestartfiles
    call fielddump
    call statsdump        ! tg3315
 
  end do

!-------------------------------------------------------
!             END OF TIME LOOP
!-------------------------------------------------------

!--------------------------------------------------------
!    4    FINALIZE ADD ONS AND THE MAIN PROGRAM
!-------------------------------------------------------
  call exitfielddump  
  call exitstatsdump     !tg3315
  call exitmodules

end program DALESURBAN
