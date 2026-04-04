module architecture
  use mpi
  implicit none
  save

  integer :: comm3d
  integer :: nbrnorth
  integer :: nbrsouth
  integer :: nbreast
  integer :: nbrwest
  integer :: myid
  integer :: myidx
  integer :: myidy
  integer :: nprocs
  integer :: nprocx
  integer :: nprocy
  integer :: mpierr
  integer :: my_real
  real    :: CPU_program
  real    :: CPU_program0
  character(3) :: cmyid
  character(3) :: cmyidx
  character(3) :: cmyidy

  interface bcast
    module procedure bcast_real0d
    module procedure bcast_real1d
    module procedure bcast_int0d
    module procedure bcast_int1d
    module procedure bcast_logical0d
    module procedure bcast_logical1d
    module procedure bcast_char0d
  end interface

  interface allreduce
    module procedure allreduce_real0d
    module procedure allreduce_real1d
    module procedure allreduce_int0d
    module procedure allreduce_int1d
  end interface

  interface global_sum
    module procedure global_sum_real0d
    module procedure global_sum_real1d
    module procedure global_sum_int0d
    module procedure global_sum_int1d
  end interface

  interface global_max
    module procedure global_max_real0d
    module procedure global_max_int0d
  end interface

contains

  subroutine barrier_domain()
    implicit none

    call MPI_BARRIER(comm3d, mpierr)
  end subroutine barrier_domain

  subroutine abort_run(err_code)
    implicit none
    integer, intent(in) :: err_code

    call MPI_ABORT(comm3d, err_code, mpierr)
  end subroutine abort_run

  subroutine bcast_real0d(var, root)
    implicit none
    real, intent(inout) :: var
    integer, intent(in) :: root

    call MPI_BCAST(var, 1, MY_REAL, root, comm3d, mpierr)
  end subroutine bcast_real0d

  subroutine bcast_real1d(var, root)
    implicit none
    real, intent(inout) :: var(:)
    integer, intent(in) :: root

    if (size(var) == 0) return
    call MPI_BCAST(var, size(var), MY_REAL, root, comm3d, mpierr)
  end subroutine bcast_real1d

  subroutine bcast_int0d(var, root)
    implicit none
    integer, intent(inout) :: var
    integer, intent(in)    :: root

    call MPI_BCAST(var, 1, MPI_INTEGER, root, comm3d, mpierr)
  end subroutine bcast_int0d

  subroutine bcast_int1d(var, root)
    implicit none
    integer, intent(inout) :: var(:)
    integer, intent(in)    :: root

    if (size(var) == 0) return
    call MPI_BCAST(var, size(var), MPI_INTEGER, root, comm3d, mpierr)
  end subroutine bcast_int1d

  subroutine bcast_logical0d(var, root)
    implicit none
    logical, intent(inout) :: var
    integer, intent(in)    :: root

    call MPI_BCAST(var, 1, MPI_LOGICAL, root, comm3d, mpierr)
  end subroutine bcast_logical0d

  subroutine bcast_logical1d(var, root)
    implicit none
    logical, intent(inout) :: var(:)
    integer, intent(in)    :: root

    if (size(var) == 0) return
    call MPI_BCAST(var, size(var), MPI_LOGICAL, root, comm3d, mpierr)
  end subroutine bcast_logical1d

  subroutine bcast_char0d(var, root)
    implicit none
    character(len=*), intent(inout) :: var
    integer, intent(in)             :: root

    call MPI_BCAST(var, len(var), MPI_CHARACTER, root, comm3d, mpierr)
  end subroutine bcast_char0d

  function allreduce_real0d(lvar, op) result(gvar)
    implicit none
    real, intent(in)    :: lvar
    integer, intent(in) :: op
    real                :: gvar

    call MPI_ALLREDUCE(lvar, gvar, 1, MY_REAL, op, comm3d, mpierr)
  end function allreduce_real0d

  function allreduce_real1d(lvar, op) result(gvar)
    implicit none
    real, intent(in)    :: lvar(:)
    integer, intent(in) :: op
    real                :: gvar(size(lvar))

    if (size(lvar) == 0) return
    call MPI_ALLREDUCE(lvar, gvar, size(lvar), MY_REAL, op, comm3d, mpierr)
  end function allreduce_real1d

  function allreduce_int0d(lvar, op) result(gvar)
    implicit none
    integer, intent(in) :: lvar
    integer, intent(in) :: op
    integer             :: gvar

    call MPI_ALLREDUCE(lvar, gvar, 1, MPI_INTEGER, op, comm3d, mpierr)
  end function allreduce_int0d

  function allreduce_int1d(lvar, op) result(gvar)
    implicit none
    integer, intent(in) :: lvar(:)
    integer, intent(in) :: op
    integer             :: gvar(size(lvar))

    if (size(lvar) == 0) return
    call MPI_ALLREDUCE(lvar, gvar, size(lvar), MPI_INTEGER, op, comm3d, mpierr)
  end function allreduce_int1d

  function global_sum_real0d(lvar) result(gvar)
    implicit none
    real, intent(in) :: lvar
    real             :: gvar

    call MPI_ALLREDUCE(lvar, gvar, 1, MY_REAL, MPI_SUM, comm3d, mpierr)
  end function global_sum_real0d

  function global_sum_real1d(lvar) result(gvar)
    implicit none
    real, intent(in) :: lvar(:)
    real             :: gvar(size(lvar))

    if (size(lvar) == 0) return
    call MPI_ALLREDUCE(lvar, gvar, size(lvar), MY_REAL, MPI_SUM, comm3d, mpierr)
  end function global_sum_real1d

  function global_sum_int0d(lvar) result(gvar)
    implicit none
    integer, intent(in) :: lvar
    integer             :: gvar

    call MPI_ALLREDUCE(lvar, gvar, 1, MPI_INTEGER, MPI_SUM, comm3d, mpierr)
  end function global_sum_int0d

  function global_sum_int1d(lvar) result(gvar)
    implicit none
    integer, intent(in) :: lvar(:)
    integer             :: gvar(size(lvar))

    if (size(lvar) == 0) return
    call MPI_ALLREDUCE(lvar, gvar, size(lvar), MPI_INTEGER, MPI_SUM, comm3d, mpierr)
  end function global_sum_int1d

  function global_max_real0d(lvar) result(gvar)
    implicit none
    real, intent(in) :: lvar
    real             :: gvar

    call MPI_ALLREDUCE(lvar, gvar, 1, MY_REAL, MPI_MAX, comm3d, mpierr)
  end function global_max_real0d

  function global_max_int0d(lvar) result(gvar)
    implicit none
    integer, intent(in) :: lvar
    integer             :: gvar

    call MPI_ALLREDUCE(lvar, gvar, 1, MPI_INTEGER, MPI_MAX, comm3d, mpierr)
  end function global_max_int0d

end module architecture
