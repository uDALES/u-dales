!> \file modmpi.f90
!!  Layer to deal with the parallelization.

!>
!!  Layer to deal with the parallelization.
!>
!!  \author Matthieu Pourquie, TU Delft
!!  \par Revision list
!!  \todo Documentation
!!  \todo 2D/3D parallelization
!!  \todo Include interfaces for MPI_ALLREDUCE, MPI_ALLTOALL, MPI_BCAST,
!! MPI_SENDRECV to get rid of pure mpi calls in the code
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!



module modmpi
use mpi
implicit none
save
  integer comm3d
  integer nbrtop
  integer nbrbottom
  integer myid
  integer nprocs
  integer mpierr
  integer my_real
  real    CPU_program    !end time
  real    CPU_program0   !start time

   character(3) :: cmyid

contains
  subroutine initmpi
    implicit none
    integer dims(1)
    logical periods(1)
    integer periods2(1)

    call MPI_INIT(mpierr)
    MY_REAL = MPI_DOUBLE_PRECISION
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, mpierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, mpierr )
! Specify the # procs in each direction.
! specifying a 0 means that MPI will try to find a useful # procs in
! the corresponding  direction,

! specifying 1 means only 1 processor in this direction, meaning that
! we have in fact a grid of (at most) 2 dimensions left. This is used
! when we want the array index range in 1 particular direction to be
! present on all processors in the grid

    dims(1) = 0


! directions 1 and 2 are chosen periodic


    periods(1) = .true.
! Soares 20080115
    periods2(1) = 1

! find suitable # procs in each direction

    call MPI_DIMS_CREATE( nprocs, 1, dims, mpierr )

! create the Cartesian communicator, denoted by the integer comm3d

    ! BUG - Thijs, Harm
    !call MPI_CART_CREATE(MPI_COMM_WORLD, 1, dims, periods,.false., &
    !                    comm3d, ierr )

    call MPI_CART_CREATE(MPI_COMM_WORLD, 1, dims, periods,.true., &
                        comm3d, mpierr )

! Soares 20080115
!     call MPI_CART_CREATE(MPI_COMM_WORLD, 1, dims, periods2,1, &
!                         comm3d, mpierr )

! Get my processor number in this communicator

    call MPI_COMM_RANK( comm3d, myid, mpierr )


! when applying boundary conditions, we need to know which processors
! are neighbours in all 3 directions


! these are determined with the aid of the MPI routine MPI_CART_SHIFT,

    call MPI_CART_SHIFT( comm3d, 0,  1, nbrbottom, nbrtop,   mpierr )

! determine some useful MPI datatypes for sending/receiving data

     write(cmyid,'(i3.3)') myid

    if(myid==0)then
      CPU_program0 = MPI_Wtime()
    end if

    write(*,*)'nprocs = ', nprocs
  end subroutine initmpi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exitmpi
    implicit none


    if(myid==0)then
      CPU_program = MPI_Wtime() - CPU_program0
      write(6,*)'TOTAL CPU time = ', CPU_program
    end if

    call MPI_Comm_free( comm3d, mpierr )
    call MPI_FINALIZE(mpierr)
  end subroutine exitmpi

  subroutine barrou()
    implicit none
    call MPI_BARRIER(comm3d,mpierr)

  return
  end subroutine barrou



  subroutine excj( a, sx, ex, sy, ey, sz,ez)
    implicit none

  integer sx, ex, sy, ey, sz, ez
  real a(sx:ex, sy:ey, sz:ez)
  integer iiget, status(MPI_STATUS_SIZE)
  integer ii, i, k
  real,allocatable, dimension(:) :: buffj1,buffj2,buffj3,buffj4
  iiget = (ex - sx + 1)*(ez - sz + 1)

  allocate( buffj1(iiget),&
            buffj2(iiget),&
            buffj3(iiget),&
            buffj4(iiget))



  if(nbrtop/=MPI_PROC_NULL)then
    do k=sz,ez
    do i=sx,ex
      ii = i - sx + 1 + (k - sz )*(ex - sx + 1)
      buffj1(ii) = a(i,ey-1,k)
    enddo
    enddo
  endif
  call MPI_SENDRECV(  buffj1,  ii    , MY_REAL, nbrtop, 4, &
                      buffj2,  iiget , MY_REAL, nbrbottom,  4,      &
                      comm3d, status, mpierr )
  if(nbrbottom/=MPI_PROC_NULL)then
    do k=sz,ez
    do i=sx,ex
      ii = i - sx + 1 + (k - sz )*(ex - sx + 1)
      a(i,sy,k) = buffj2(ii)
    enddo
    enddo
  endif

!   call barrou()

  if(nbrbottom/=MPI_PROC_NULL)then
    do k=sz,ez
    do i=sx,ex
      ii = i - sx + 1 + (k - sz )*(ex - sx + 1)
      buffj3(ii) = a(i,sy+1,k)
    enddo
    enddo
  endif
  call MPI_SENDRECV(  buffj3,  ii    , MY_REAL, nbrbottom,  5, &
                          buffj4,  iiget , MY_REAL, nbrtop, 5, &
                          comm3d, status, mpierr )
  if(nbrtop/=MPI_PROC_NULL)then
    do k=sz,ez
    do i=sx,ex
      ii = i - sx + 1 + (k - sz )*(ex - sx + 1)
      a(i,ey,k) = buffj4(ii)
    enddo
    enddo
  endif

!   call barrou()
  deallocate (buffj1,buffj2,buffj3,buffj4)

  return
  end subroutine excj

  subroutine excjs( a, sx, ex, sy, ey, sz,ez,ih,jh)
    implicit none
  integer sx, ex, sy, ey, sz, ez, ih, jh
  real a(sx-ih:ex+ih, sy-jh:ey+jh, sz:ez)
  integer status(MPI_STATUS_SIZE), iiget
  integer ii, i, j, k
  real,allocatable, dimension(:) :: buffj1,buffj2,buffj3,buffj4
  iiget = jh*(ex - sx + 1 + 2*ih)*(ez - sz + 1)

  allocate( buffj1(iiget),&
            buffj2(iiget),&
            buffj3(iiget),&
            buffj4(iiget))


  if(nbrtop/=MPI_PROC_NULL)then
    ii = 0
    do j=1,jh
    do k=sz,ez
    do i=sx-ih,ex+ih
      ii = ii + 1
      buffj1(ii) = a(i,ey-j+1,k)
    enddo
    enddo
    enddo
  endif

  call MPI_SENDRECV(  buffj1,  ii    , MY_REAL, nbrtop, 4, &
                           buffj2,  iiget , MY_REAL, nbrbottom,  4, &
                           comm3d, status, mpierr )
  if(nbrbottom/=MPI_PROC_NULL)then
    ii = 0
    do j=1,jh
    do k=sz,ez
    do i=sx-ih,ex+ih
      ii = ii + 1
      a(i,sy-j,k) = buffj2(ii)
    enddo
    enddo
    enddo
  endif

!   call barrou()

  if(nbrbottom/=MPI_PROC_NULL)then
    ii = 0
    do j=1,jh
    do k=sz,ez
    do i=sx-ih,ex+ih
      ii = ii + 1
      buffj3(ii) = a(i,sy+j-1,k)
    enddo
    enddo
    enddo
  endif
  call MPI_SENDRECV(  buffj3,  ii    , MY_REAL, nbrbottom,  5, &
                          buffj4,  iiget , MY_REAL, nbrtop, 5, &
                          comm3d, status, mpierr )
  if(nbrtop/=MPI_PROC_NULL)then
    ii = 0
    do j=1,jh
    do k=sz,ez
    do i=sx-ih,ex+ih
      ii = ii + 1
      a(i,ey+j,k) = buffj4(ii)
    enddo
    enddo
    enddo
  endif

!   call barrou()
  deallocate (buffj1,buffj2,buffj3,buffj4)

  return
  end subroutine excjs

  subroutine slabsum(aver,ks,kf,var,ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes)
    implicit none

    integer :: ks,kf
    integer :: ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes
    real    :: aver(ks:kf)
    real    :: var (ib:ie,jb:je,kb:ke)
    real    :: averl(ks:kf)
    real    :: avers(ks:kf)
    integer :: k

    averl       = 0.
    avers       = 0.

    do k=kbs,kes
      averl(k) = sum(var(ibs:ies,jbs:jes,k))
    enddo

    call MPI_ALLREDUCE(averl, avers, kf-ks+1,  MY_REAL, &
                          MPI_SUM, comm3d,mpierr)

    aver = aver + avers

    return
  end subroutine slabsum


end module
