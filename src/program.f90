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
program DALESURBAN      !Version 48

!!----------------------------------------------------------------
!!     0.0    USE STATEMENTS FOR CORE MODULES
!!----------------------------------------------------------------
  use modmpi,            only : initmpi, exitmpi, myid, barrou
  !use modtest,           only : inittest, exittest
  use modglobal,         only : rk3step,timeleft,ib,jb,kb,ke,driverid,ibrank,timee
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
  use modfields, only : u0, um
  !use modbudget,       only : initbudget, budgetstat, exitbudget
  implicit none

!----------------------------------------------------------------
!     1      READ NAMELISTS,INITIALISE GRID, CONSTANTS AND FIELDS
!----------------------------------------------------------------
  call initmpi
  call startup
  write(*,*) "done startup"
  !call inittest

!---------------------------------------------------------
!      2     INITIALIZE STATISTICAL ROUTINES AND ADD-ONS
!---------------------------------------------------------
  call initchecksim ! Could be deprecated
  call initstat_nc ! Could be deprecated
  call initfielddump
  !call fielddump
  !call initstatsdump !tg3315

  !call readfacetfiles
  !call initEB
  !write(*,*) "done init stuff"

  !write(6,*) 'Determine immersed walls'
  !call createwalls    ! determine walls/blocks
 ! call nearwall       ! determine minimum distance and corresponding shear components, ils13 10.07.17, commented, not functional at the moment, not needed for vreman but for smag., fix in modibm
  !write(6,*) 'Finished determining immersed walls'

  call boundary  !ils13 22.06.2017 inserted boundary here to get values at ghost cells before iteration starts
  ! call fielddump
  ! write(*,*) "done fielddump after boundary"
!  not necessary but abates the fact that temp field is randomised by randomisation of just velocity fields
!  (because advection at start of time loop without being divergence free)
  !call poisson
  write(*,*) "Done initial boundary"
  call fielddump
  write(*,*) "Done initial fielddump"
!------------------------------------------------------
!   3.0   MAIN TIME LOOP
!------------------------------------------------------
  !write(*,*)'START myid ', myid
  !do while (.false.)
  do while ((timeleft>0) .or. (rk3step < 3))
    !write(*,*) timeleft
    call tstep_update
    !if (driverid==0 .and. ibrank) write(*,*) "--------------- start of timestep. rk3step, time: ", rk3step, timee
!-----------------------------------------------------
!   3.2   ADVECTION AND DIFFUSION
!-----------------------------------------------------
    !call boundary
    call advection                ! now also includes predicted pressure gradient term
    call subgrid
!-----------------------------------------------------
!   3.3   THE SURFACE LAYER
!-----------------------------------------------------

    call bottom
!-----------------------------------------------------
!   3.4   REMAINING TERMS
!-----------------------------------------------------

    !call coriolis       !remaining terms of ns equation

    !call forces         !remaining terms of ns equation

    !call lstend         !large scale forcings

    !call nudge          ! nudge top cells of fields to enforce steady-state

    !call ibmwallfun     ! immersed boundary forcing: only shear forces.

    !call masscorr       ! correct pred. velocity pup to get correct mass flow

    !call ibmnorm        ! immersed boundary forcing: set normal velocities to zero

    !call EB

    !call scalsource     ! adds continuous forces in specified region of domain

!------------------------------------------------------
!   3.4   EXECUTE ADD ONS
!------------------------------------------------------
    !call fixuinf2

    !call fixuinf1

!-----------------------------------------------------------------------
!   3.5  PRESSURE FLUCTUATIONS, TIME INTEGRATION AND BOUNDARY CONDITIONS
!-----------------------------------------------------------------------
    !call grwdamp        !damping at top of the model

    call poisson
    call tstep_integrate
    call checksim
    call fielddump
    call boundary
    !call fixthetainf

!-----------------------------------------------------
!   3.6   LIQUID WATER CONTENT AND DIAGNOSTIC FIELDS
!-----------------------------------------------------
    !call thermodynamics

!-----------------------------------------------------
!   3.7  WRITE RESTARTFILES AND DO STATISTICS
!------------------------------------------------------

    !call checksim
   ! call writedatafiles   ! write data files for later analysis
    !call writerestartfiles
    !call fielddump
    !call statsdump        ! tg3315

  end do
!-------------------------------------------------------
!             END OF TIME LOOP
!-------------------------------------------------------

!--------------------------------------------------------
!    4    FINALIZE ADD ONS AND THE MAIN PROGRAM
!-------------------------------------------------------
  call exitfielddump
  !call exitstatsdump     !tg3315
  !call exitmodules
  !call exittest
  call exitmpi

end program DALESURBAN
