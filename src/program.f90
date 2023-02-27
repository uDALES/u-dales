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
  use modmpi,            only : initmpi, exitmpi, myid, starttimer
  use modglobal,         only : rk3step,timeleft
  use modstartup,        only : startup,exitmodules
  use modsave,           only : writerestartfiles
  use modboundary,       only : boundary, grwdamp,tqaver,halos
  use modthermodynamics, only : thermodynamics
!  use modsurface,        only : surface
  use modsubgrid,        only : subgrid
  use modforces,         only : forces,coriolis,lstend,fixuinf1,fixuinf2,fixthetainf,nudge, masscorr
  use modpois,           only : poisson
  use modibm,            only : ibmwallfun,ibmnorm,bottom,initibm
  use initfac,           only : readfacetfiles
  use modEB,             only : initEB,EB

!----------------------------------------------------------------
!     0.1     USE STATEMENTS FOR ADDONS STATISTICAL ROUTINES
!----------------------------------------------------------------
  use modchecksim,     only : initchecksim, checksim
  use modstat_nc,      only : initstat_nc
  use modfielddump,    only : initfielddump, fielddump,exitfielddump
  use modstatsdump,    only : initstatsdump,statsdump,exitstatsdump    !tg3315
  use modtimedep,      only : inittimedep, timedep
  implicit none

!----------------------------------------------------------------
!     1      READ NAMELISTS,INITIALISE GRID, CONSTANTS AND FIELDS
!----------------------------------------------------------------
  call initmpi
  call startup
  !write(*,*) "done startup"
  !call inittest
  !write(*,*) myid, "done initibm"
!---------------------------------------------------------
!      2     INITIALIZE STATISTICAL ROUTINES AND ADD-ONS
!---------------------------------------------------------
  call initchecksim ! Could be deprecated
  call initstat_nc ! Could be deprecated

  !write(*,*) myid, "done initfielddump"

  call initstatsdump !tg3315
  !write(*,*) myid, "done initstatsdump"

  call readfacetfiles
  call initEB
  call inittimedep
  call initibm
  call initfielddump
  call boundary

  !call fielddump
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

    call subgrid

!-----------------------------------------------------
!   3.3   THE SURFACE LAYER
!-----------------------------------------------------

    call bottom
!-----------------------------------------------------
!   3.4   REMAINING TERMS
!-----------------------------------------------------

    !call coriolis       !remaining terms of ns equation

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
    !call fixuinf2

    !call fixuinf1

!-----------------------------------------------------------------------
!   3.5  PRESSURE FLUCTUATIONS, TIME INTEGRATION AND BOUNDARY CONDITIONS
!-----------------------------------------------------------------------
    !call grwdamp        !damping at top of the model

    call poisson

    call tstep_integrate

    call halos

    call checksim

    call fielddump

    call statsdump

    call boundary

    !call fixthetainf

!-----------------------------------------------------
!   3.6   LIQUID WATER CONTENT AND DIAGNOSTIC FIELDS
!-----------------------------------------------------
    call thermodynamics

!-----------------------------------------------------
!   3.7  WRITE RESTARTFILES AND DO STATISTICS
!------------------------------------------------------

    !call writerestartfiles

  end do
!-------------------------------------------------------
!             END OF TIME LOOP
!-------------------------------------------------------

!--------------------------------------------------------
!    4    FINALIZE ADD ONS AND THE MAIN PROGRAM
!-------------------------------------------------------
  call exitfielddump
  call exitstatsdump     !tg3315
  !call exitmodules
  !call exittest
  call exitmpi

end program DALESURBAN
