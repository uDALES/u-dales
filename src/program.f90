!> \file program.f90
!! Main program

!>
!! \mainpage
!! Dutch Atmospheric Large Eddy Simulation
!! \section DALES Dutch Atmospheric Large Eddy Simulation
!!
!! @version 3.2git
!!
!! @author
!! Stephan de Roode
!! (TU Delft)
!! \author
!! Chiel van Heerwaarden
!! (Wageningen University)
!! \author
!! Thijs Heus
!! (Max Planck Institute Hamburg)
!>
!! \section Log Change log
!! \par New Features
!!  - Land Surface model
!!  - Chemistry
!!  - NetCDF output
!!  - Gitourius repository, including case database
!!  - Optional Smagorinsky SFS-scheme
!! \par Changes
!! - Newton-Raphson Ql calculations now by default true\n
!! - Subsidence by default off for momentum. Changable with lmomsubs switch in
!! physics
!!
!! \todo
!! - Radiation
!! - Review documentation
!! - Cleaning up for release : Removal of debug statements, cleaner code, less
!! warnings in chemistry
!!
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
program DALES      !Version 3.2 Beta 1

!!----------------------------------------------------------------
!!     0.0    USE STATEMENTS FOR CORE MODULES
!!----------------------------------------------------------------
  use modmpi,            only : myid, initmpi
  use modglobal,         only : rk3step,timee,btime,runtime
  use modfields,         only : thl0
  use modstartup,        only : startup, writerestartfiles,exitmodules
  use modtimedep,        only : timedep
  use modboundary,       only : boundary, grwdamp,tqaver
  use modthermodynamics, only : thermodynamics
  use modmicrophysics,   only : microsources
  use modsurface,        only : surface
  use modsubgrid,        only : subgrid
  use modforces,         only : forces, coriolis, lstend
  use modradiation,      only : radiation
  use modpois,           only : poisson


!----------------------------------------------------------------
!     0.1     USE STATEMENTS FOR ADDONS STATISTICAL ROUTINES
!----------------------------------------------------------------
  use modchecksim,     only : initchecksim, checksim
  use modstat_nc,      only : initstat_nc
  use modtimestat,     only : inittimestat, timestat
  use modgenstat,      only : initgenstat, genstat, exitgenstat
  use modradstat,      only : initradstat ,radstat, exitradstat
  use modsampling,     only : initsampling, sampling,exitsampling
  use modcrosssection, only : initcrosssection, crosssection,exitcrosssection
  use modprojection,   only : initprojection, projection
  use modcloudfield,   only : initcloudfield, cloudfield
  use modfielddump,    only : initfielddump, fielddump,exitfielddump
  use modstattend,     only : initstattend, stattend,exitstattend, tend_start,tend_adv,tend_subg,tend_force,&
                              tend_rad,tend_ls,tend_micro, tend_topbound,tend_pois,tend_addon, tend_coriolis

  use modbulkmicrostat,only : initbulkmicrostat, bulkmicrostat,exitbulkmicrostat
  use modbudget,       only : initbudget, budgetstat, exitbudget

  use modtilt,         only : inittilt, tiltedgravity, tiltedboundary, exittilt
  use modparticles,    only : initparticles, particles, exitparticles
  use modnudge,        only : initnudge, nudge, exitnudge
!   use modnetcdfstats,  only : initnetcdfstats, netcdfstats, exitnetcdfstats
!   use modnetcdfmovie,  only : initnetcdfmovie, netcdfmovie, exitnetcdfmovie
  use modchem,         only : initchem,inputchem, twostep
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
  call initstat_nc  ! should be called before stat-routines that might do netCDF
  call inittimestat  !Timestat must preceed all other timeseries that could write in the same netCDF file (unless stated otherwise
  call initgenstat  !Genstat must preceed all other statistics that could write in the same netCDF file (unless stated otherwise
  call inittilt
  call initsampling
  call initcrosssection
  call initprojection
  call initcloudfield
  call initfielddump
  call initstattend
  call initradstat
  call initparticles
  call initnudge
  ! call initnetcdfstats
  ! call initnetcdfmovie
  call initbulkmicrostat
  call initbudget

  call initchem
  call inputchem

!------------------------------------------------------
!   3.0   MAIN TIME LOOP
!------------------------------------------------------
  write(*,*)'START myid ', myid
  do while ((timee-btime)<runtime .or. rk3step < 3)
    call tstep_update
    call timedep
    call stattend(tend_start)

!-----------------------------------------------------
!   3.1   RADIATION         
!-----------------------------------------------------
    call radiation !radiation scheme
    call stattend(tend_rad)

!-----------------------------------------------------
!   3.2   THE SURFACE LAYER
!-----------------------------------------------------
    call surface

!-----------------------------------------------------
!   3.3   ADVECTION AND DIFFUSION
!-----------------------------------------------------
    call advection
    call stattend(tend_adv)
    call subgrid
    call stattend(tend_subg)

!-----------------------------------------------------
!   3.4   REMAINING TERMS
!-----------------------------------------------------
    call coriolis !remaining terms of ns equation
    call stattend(tend_coriolis)
    call forces !remaining terms of ns equation
    call stattend(tend_force)

    call lstend !large scale forcings
    call stattend(tend_ls)
    call microsources !Drizzle etc.
    call stattend(tend_micro)

!------------------------------------------------------
!   3.4   EXECUTE ADD ONS
!------------------------------------------------------
    call nudge
    call tiltedgravity
    call stattend(tend_addon)

!-----------------------------------------------------------------------
!   3.5  PRESSURE FLUCTUATIONS, TIME INTEGRATION AND BOUNDARY CONDITIONS
!-----------------------------------------------------------------------
    call grwdamp !damping at top of the model
    call tqaver !set thl, qt and sv(n) equal to slab average at level kmax
    call stattend(tend_topbound)
    call poisson
    call stattend(tend_pois,lastterm=.true.)

    call tstep_integrate
    call boundary
    call tiltedboundary
!-----------------------------------------------------
!   3.6   LIQUID WATER CONTENT AND DIAGNOSTIC FIELDS
!-----------------------------------------------------
    call thermodynamics


!-----------------------------------------------------
!   3.7  WRITE RESTARTFILES AND DO STATISTICS
!------------------------------------------------------
    ! Chemistry
    call twostep

    call checksim
    call timestat  !Timestat must preceed all other timeseries that could write in the same netCDF file (unless stated otherwise
    call genstat  !Genstat must preceed all other statistics that could write in the same netCDF file (unless stated otherwise
    call radstat
    call sampling
    call crosssection
    call projection
    call cloudfield
    call fielddump
    call particles

    call bulkmicrostat
    call budgetstat

    ! call netcdfmovie

    call writerestartfiles
  end do

!-------------------------------------------------------
!             END OF TIME LOOP
!-------------------------------------------------------


!--------------------------------------------------------
!    4    FINALIZE ADD ONS AND THE MAIN PROGRAM
!-------------------------------------------------------
  call exitgenstat
  call exitradstat
  call exitparticles
  call exitnudge
  call exitsampling
  call exitstattend
  call exitbulkmicrostat
  call exitbudget
  call exitcrosssection
  call exitfielddump
  ! call exitnetcdfmovie
  call exitmodules

end program DALES
