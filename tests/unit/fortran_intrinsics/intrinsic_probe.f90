!> Probe runtime support for a single Fortran intrinsic.
!!
!! Some compilers accept an intrinsic at compile time but have no runtime
!! implementation for a particular argument type, and abort only when the call is
!! actually reached. NVHPC 24.11 does this with FINDLOC on logical arrays:
!!
!!     0: FINDLOC: unimplemented for data type
!!
!! A build gives no warning, so the gap only surfaces when a code path that uses
!! it happens to execute. Each probe therefore runs in its own process, selected
!! by a command-line argument, so that an abort in one does not hide the others.
!!
!! Exit 0 means the intrinsic produced the correct answer, exit 1 means it
!! produced a wrong answer, and any other termination (signal, runtime abort)
!! means the runtime does not implement it.
program intrinsic_probe
   implicit none

   character(len=64) :: which
   integer :: argl, status

   call get_command_argument(1, which, argl, status)
   if (status /= 0) then
      write(*,'(A)') 'usage: intrinsic_probe <check-name>'
      stop 2
   end if

   select case (trim(which))
   case ('findloc-logical')
      call probe_findloc_logical
   case ('findloc-real')
      call probe_findloc_real
   case ('findloc-integer')
      call probe_findloc_integer
   case ('count-logical')
      call probe_count_logical
   case default
      write(*,'(A,A)') 'unknown check: ', trim(which)
      stop 2
   end select

contains

   subroutine expect(got, want)
      integer, intent(in) :: got, want
      if (got /= want) then
         write(*,'(A,I0,A,I0)') 'WRONG ANSWER: got ', got, ' want ', want
         stop 1
      end if
      write(*,'(A)') 'ok'
   end subroutine expect

   !> findloc over a logical mask - the form modibm used to locate the grid cell
   !! containing a reconstruction point. Unimplemented in the NVHPC runtime.
   subroutine probe_findloc_logical
      real :: grid(5) = [0., 0.5, 1.0, 1.5, 2.0]
      call expect(findloc(1.2 >= grid, .true., 1, back=.true.), 3)
   end subroutine probe_findloc_logical

   subroutine probe_findloc_real
      real :: a(5) = [0., 0.5, 1.0, 1.5, 2.0]
      call expect(findloc(a, 1.0, 1), 3)
   end subroutine probe_findloc_real

   subroutine probe_findloc_integer
      integer :: a(5) = [1, 2, 3, 4, 5]
      call expect(findloc(a, 3, 1), 3)
   end subroutine probe_findloc_integer

   !> The replacement modibm::cell_index uses.
   subroutine probe_count_logical
      real :: grid(5) = [0., 0.5, 1.0, 1.5, 2.0]
      call expect(count(1.2 >= grid), 3)
   end subroutine probe_count_logical

end program intrinsic_probe
