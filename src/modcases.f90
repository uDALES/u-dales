!> \file modcases.f90
!! Registry of runtime case hooks: case-specific physics that lives in
!! src/cases/ instead of in the core modules.
!>
!! The namelist &NAMCASES lists the hooks to switch on, `case_hooks = 'a,b'`.
!! Each name maps (select case in register_case) to a register_<name> routine
!! that reads and checks the case's namelist and stores the case's procedures
!! in a hook set. The core modules call the call_hooks_* routines at fixed
!! points; without hooks they do nothing. Hooks run in case_hooks order.
!!
!! Hook points:
!! - init      once after startup (after readinitfiles: cold and warm start);
!!             checks that need the initial fields
!! - radiation in modradiation:radiation at each radiation update (rk3step 1,
!!             every timerad), after the generic radpar tendency and before
!!             thlprad is added to thlp; may add to thlprad (held until the next update)
!! - forcing   in the time loop after lstend, every RK substep; may add to the tendencies
!! - stats     in the time loop after thermostat, every RK substep; a stats hook
!!             gates itself on rk3step == 3 and its own output time
!!
!! See docs/udales-case-hooks.md for how to add a case.
!
! This file is part of uDALES (https://github.com/uDALES/u-dales).
!
! uDALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! uDALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) 2016- the uDALES Team, Imperial College London.
!
module modcases
  implicit none
  private
  public :: initcases, call_hooks_init, call_hooks_radiation, call_hooks_forcing, call_hooks_stats
  save

  abstract interface
    subroutine hook_proc()
    end subroutine hook_proc
  end interface

  !> One case's procedures; a component left null is a hook point the case does not use.
  type hook_set
    character(32) :: name = ''
    procedure(hook_proc), pointer, nopass :: init      => null()
    procedure(hook_proc), pointer, nopass :: radiation => null()
    procedure(hook_proc), pointer, nopass :: forcing   => null()
    procedure(hook_proc), pointer, nopass :: stats     => null()
  end type hook_set

  integer, parameter :: max_hooks = 8
  character(*), parameter :: valid_names = 'dycoms_rf01'
  type(hook_set) :: hooks(max_hooks)
  integer :: nhooks = 0

  character(256) :: case_hooks = ''   !< comma-separated hook names; '' = none

contains

  !> Read &NAMCASES and register the listed hooks (all ranks).
  subroutine initcases
    use modglobal, only : ifnamopt, fname_options
    use modmpi,    only : myid, comm3d, mpierr, mpi_character
    implicit none
    integer :: ierr
    namelist/NAMCASES/ case_hooks

    if (myid == 0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read (ifnamopt, NAMCASES, iostat=ierr)
      if (ierr > 0) then
        write(0, *) 'ERROR: Problem in namoptions NAMCASES'
        write(0, *) 'iostat error: ', ierr
        stop 1
      end if
      close(ifnamopt)
    end if
    call MPI_BCAST(case_hooks, len(case_hooks), MPI_CHARACTER, 0, comm3d, mpierr)
    call register_list(case_hooks)
  end subroutine initcases

  !> Abort on all ranks with one message (two lines) from rank 0.
  subroutine fail(msg1, msg2)
    use modmpi, only : myid, abort_all_ranks
    implicit none
    character(*), intent(in) :: msg1, msg2
    if (myid == 0) then
      write(0, '(2a)') 'ERROR: ', msg1
      if (len(msg2) > 0) write(0, '(a)') msg2
    end if
    call abort_all_ranks(1)
  end subroutine fail

  !> Split a comma-separated list and register each non-empty name.
  subroutine register_list(list)
    implicit none
    character(*), intent(in) :: list
    integer :: i0, i1

    i0 = 1
    do while (i0 <= len_trim(list))
      i1 = index(list(i0:), ',')
      if (i1 == 0) then
        i1 = len_trim(list) + 1
      else
        i1 = i0 + i1 - 1
      end if
      if (len_trim(list(i0:i1-1)) > 0) call register_case(trim(adjustl(list(i0:i1-1))))
      i0 = i1 + 1
    end do
  end subroutine register_list

  !> Map a hook name to its register routine; abort on unknown or repeated names.
  subroutine register_case(name)
    implicit none
    character(*), intent(in) :: name
    integer :: n

    do n = 1, nhooks
      if (trim(hooks(n)%name) == name) call fail('case hook listed twice in case_hooks: ' // name, '')
    end do
    if (nhooks == max_hooks) call fail('too many case hooks in case_hooks', '')
    nhooks = nhooks + 1
    hooks(nhooks)%name = name
    select case (name)
    case ('dycoms_rf01')
      call register_dycoms_rf01(hooks(nhooks))
    case default
      call fail('unknown case hook in case_hooks: ' // name, 'Valid case hooks: ' // valid_names)
    end select
  end subroutine register_case

  !> DYCOMS-II RF01 free-troposphere radiative term (src/cases/case_dycoms_rf01.f90).
  subroutine register_dycoms_rf01(h)
    use case_dycoms_rf01, only : dycoms_rf01_config, dycoms_rf01_init, dycoms_rf01_radiation
    implicit none
    type(hook_set), intent(inout) :: h
    call dycoms_rf01_config          ! read and check &NAMCASE_DYCOMS now, before the rest of startup
    h%init      => dycoms_rf01_init
    h%radiation => dycoms_rf01_radiation
  end subroutine register_dycoms_rf01

  subroutine call_hooks_init
    implicit none
    integer :: n
    do n = 1, nhooks
      if (associated(hooks(n)%init)) call hooks(n)%init()
    end do
  end subroutine call_hooks_init

  subroutine call_hooks_radiation
    implicit none
    integer :: n
    do n = 1, nhooks
      if (associated(hooks(n)%radiation)) call hooks(n)%radiation()
    end do
  end subroutine call_hooks_radiation

  subroutine call_hooks_forcing
    implicit none
    integer :: n
    do n = 1, nhooks
      if (associated(hooks(n)%forcing)) call hooks(n)%forcing()
    end do
  end subroutine call_hooks_forcing

  subroutine call_hooks_stats
    implicit none
    integer :: n
    do n = 1, nhooks
      if (associated(hooks(n)%stats)) call hooks(n)%stats()
    end do
  end subroutine call_hooks_stats

end module modcases
