!> \file program.f90
!! Main program

!>
!! \mainpage
!! Dutch Atmospheric Large Eddy Simulation -URBAN
!! \section DALES Dutch Atmospheric Large Eddy Simulation -URBAN
!!
!! @version 48
!!
!!
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
  use modglobal,         only : rk3step,timee,btime,runtime,timeleft,ib,jb,kb,libm,ntrun
  use modstartup,        only : startup, readscalsource,readscalpointsource,exitmodules
  use modsave,           only : writerestartfiles,writedatafiles
  use modboundary,       only : boundary, grwdamp,tqaver,masscorr
  use modthermodynamics, only : thermodynamics
  use modsurface,        only : surface
  use modsubgrid,        only : subgrid
  use modforces,         only : forces, coriolis,lstend,fixuinf1,fixuinf2,fixthetainf
  use modpois,           only : poisson
  use modibm,            only : createwalls,ibmshear,ibmnorm,nearwall,ibmforce
 
!----------------------------------------------------------------
!     0.1     USE STATEMENTS FOR ADDONS STATISTICAL ROUTINES
!----------------------------------------------------------------
  use modchecksim,     only : initchecksim, checksim
  use modstat_nc,      only : initstat_nc
  use modgenstat,      only : initgenstat, exitgenstat,average1homog,average2homog,interpolate
  use modfielddump,    only : initfielddump, fielddump,exitfielddump,tec3d
  !use modbudget,       only : initbudget, budgetstat, exitbudget
  implicit none


!----------------------------------------------------------------
!     1      READ NAMELISTS,INITIALISE GRID, CONSTANTS AND FIELDS
!----------------------------------------------------------------
  call initmpi

  call startup

!---------------------------------------------------------
!      2     INITIALIZE STATISTICAL ROUTINES AND ADD-ONS
!---------------------------------------------------------
  call initchecksim
  call initstat_nc
  call initgenstat   ! Genstat must preceed all other statistics that could write in the same netCDF file (unless stated otherwise
  call initfielddump
  !call initbudget

  write(6,*) 'Determine immersed walls'
  call createwalls    ! determine walls/blocks
  call nearwall       ! determine minimum distance and corresponding shear components
  write(6,*) 'Finished determining immersed walls'
  
!  call readscalsource      ! reads scalar line sources taking into account obstacles
!  call readscalpointsource ! reads scalar point sources taking into account obstacles
  
   ! tg3315 makes instantaneous
   ! call scalsource

!------------------------------------------------------
!   3.0   MAIN TIME LOOP
!------------------------------------------------------
  write(*,*)'START myid ', myid

  do while (timeleft>0 .or. rk3step < 3)
    call tstep_update

!-----------------------------------------------------
!   3.2   THE SURFACE LAYER
!-----------------------------------------------------
    call surface
!-----------------------------------------------------
!   3.3   ADVECTION AND DIFFUSION
!-----------------------------------------------------
    call advection                ! now also includes predicted pressure gradient term 
    call subgrid

!-----------------------------------------------------
!   3.4   REMAINING TERMS
!-----------------------------------------------------

    call coriolis       !remaining terms of ns equation
    call forces         !remaining terms of ns equation
    call lstend         !large scale forcings
    call ibmforce
    call ibmshear       ! immersed boundary forcing: only shear forces.
    call masscorr       ! correct pred. velocity pup to get correct mass flow
    call ibmnorm        ! immersed boundary forcing: set normal velocities to zero

    call scalsource     ! adds continuous forces in specified region of domain

!------------------------------------------------------
!   3.4   EXECUTE ADD ONS
!------------------------------------------------------
    call fixuinf2
!-----------------------------------------------------------------------
!   3.5  PRESSURE FLUCTUATIONS, TIME INTEGRATION AND BOUNDARY CONDITIONS
!-----------------------------------------------------------------------
    call grwdamp        !damping at top of the model
    call poisson
    
    call tstep_integrate
    call boundary
   ! call fixuinf1
   ! call fixthetainf
    
!-----------------------------------------------------
!   3.6   LIQUID WATER CONTENT AND DIAGNOSTIC FIELDS
!-----------------------------------------------------
    call thermodynamics

!-----------------------------------------------------
!   3.7  WRITE RESTARTFILES AND DO STATISTICS
!------------------------------------------------------

    call checksim
    call average1homog    ! average in j-direction and time
    call fielddump
    call writedatafiles   ! write data files for later analysis
    call writerestartfiles
   ! call budgetstat
  end do



!-------------------------------------------------------
!             END OF TIME LOOP
!-------------------------------------------------------


!--------------------------------------------------------
!    4    FINALIZE ADD ONS AND THE MAIN PROGRAM
!-------------------------------------------------------
  call exitgenstat
  !call exitbudget
  call exitfielddump
  call exitmodules

end program DALESURBAN
