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
  use modmpi,            only : initmpi,exitmpi,myid,starttimer
  use modglobal,         only : initglobal,rk3step,timeleft
  use modstartup,        only : readnamelists,init2decomp,checkinitvalues,readinitfiles,exitmodules
  use modfields,         only : initfields
  use modsave,           only : writerestartfiles
  use modboundary,       only : initboundary,boundary,grwdamp,halos
  use modthermodynamics, only : initthermodynamics,thermodynamics
  use modsubgrid,        only : initsubgrid,subgrid
  use modforces,         only : calcfluidvolumes,forces,coriolis,lstend,fixuinf1,fixuinf2,fixthetainf,nudge,masscorr,shiftedPBCs,periodicEBcorr
  use modpois,           only : initpois,poisson
  use modibm,            only : initibm,createmasks,ibmwallfun,ibmnorm,bottom
  use modtrees,          only : createtrees,trees
  use modpurifiers,      only : createpurifiers,purifiers
  use modheatpump,       only : init_heatpump,heatpump,exit_heatpump
  use initfac,           only : readfacetfiles
  use modEB,             only : initEB,EB
  use moddriver,         only : initdriver

!----------------------------------------------------------------
!     0.1     USE STATEMENTS FOR ADDONS STATISTICAL ROUTINES
!----------------------------------------------------------------
  use modchecksim,     only : initchecksim,checksim
  use modstat_nc,      only : initstat_nc
  use modfielddump,    only : initfielddump,fielddump,exitfielddump
  use modstatsdump,    only : initstatsdump,statsdump,exitstatsdump    !tg3315
  use modtimedep,      only : inittimedep,timedep
  implicit none

!----------------------------------------------------------------
!     1      READ NAMELISTS,INITIALISE GRID, CONSTANTS AND FIELDS
!----------------------------------------------------------------
  call initmpi

  !call startup
  call readnamelists

  call init2decomp

  call checkinitvalues

  call initglobal

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

  call createtrees

  call createpurifiers

  call init_heatpump

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

    call trees

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

    call fielddump

    call statsdump

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
  call exitfielddump
  call exitstatsdump     !tg3315
  call exit_heatpump
  !call exitmodules
  !call exittest
  call exitmpi

end program DALESURBAN
