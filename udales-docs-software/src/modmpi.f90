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
  !integer nbrtop
  integer nbrnorth
  !integer nbrbottom
  integer nbrsouth
  integer nbreast
  integer nbrwest
  integer myid
  integer myidx
  integer myidy
  integer nprocs ! Remove
  integer nprocx
  integer nprocy
  integer mpierr
  integer my_real
  real    CPU_program    !end time
  real    CPU_program0   !start time
  character(3) :: cmyid
  character(3) :: cmyidx
  character(3) :: cmyidy

contains
  subroutine initmpi
    use decomp_2d, only : nrank, nproc
    implicit none
    ! integer dims(1)
    ! logical periods(1)
    ! integer periods2(1)

     call MPI_INIT(mpierr)
     MY_REAL = MPI_DOUBLE_PRECISION  !MPI_REAL8 should be the same..
     comm3d = MPI_COMM_WORLD
     call MPI_COMM_RANK( MPI_COMM_WORLD, nrank, mpierr )
     call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, mpierr )
     myid = nrank
     nprocs = nproc
     write(cmyid,'(i3.3)') myid
! ! Specify the # procs in each direction.
! ! specifying a 0 means that MPI will try to find a useful # procs in
! ! the corresponding  direction,
!
! ! specifying 1 means only 1 processor in this direction, meaning that
! ! we have in fact a grid of (at most) 2 dimensions left. This is used
! ! when we want the array index range in 1 particular direction to be
! ! present on all processors in the grid
!
!     dims(1) = 0
!
!
! ! directions 1 and 2 are chosen periodic
!
!
!     periods(1) = .true.
! ! Soares 20080115
!     periods2(1) = 1
!
! ! find suitable # procs in each direction
!
!     call MPI_DIMS_CREATE( nprocs, 1, dims, mpierr )
!
! ! create the Cartesian communicator, denoted by the integer comm3d
!
!     ! BUG - Thijs, Harm
!     !call MPI_CART_CREATE(MPI_COMM_WORLD, 1, dims, periods,.false., &
!     !                    comm3d, ierr )
!
!     call MPI_CART_CREATE(MPI_COMM_WORLD, 1, dims, periods,.true., &
!                         comm3d, mpierr )
!
! ! Soares 20080115
! !     call MPI_CART_CREATE(MPI_COMM_WORLD, 1, dims, periods2,1, &
! !                         comm3d, mpierr )
!
! ! Get my processor number in this communicator
!
!     call MPI_COMM_RANK( comm3d, myid, mpierr )
!
!
! ! when applying boundary conditions, we need to know which processors
! ! are neighbours in all 3 directions
!
!
! ! these are determined with the aid of the MPI routine MPI_CART_SHIFT,
!
!     call MPI_CART_SHIFT( comm3d, 0,  1, nbrbottom, nbrtop,   mpierr )
!
! ! determine some useful MPI datatypes for sending/receiving data
!
!      write(cmyid,'(i3.3)') myid
!
     ! if(myid==0)then
     !   CPU_program0 = MPI_Wtime()
     ! end if
!
!     write(*,*)'nprocs = ', nprocs
  end subroutine initmpi

  subroutine starttimer

    if(myid==0)then
      CPU_program0 = MPI_Wtime()
    end if

  end subroutine starttimer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exitmpi
    use decomp_2d
    implicit none

    if(myid==0)then
      CPU_program = MPI_Wtime() - CPU_program0
      write(6,*)'TOTAL CPU time = ', CPU_program
    end if

    !call MPI_Comm_free( comm3d, mpierr )
    !call MPI_FINALIZE(mpierr)
    call decomp_2d_finalize
    call MPI_FINALIZE(mpierr)
  end subroutine exitmpi

  subroutine barrou()
    implicit none
    call MPI_BARRIER(comm3d,mpierr)

  return
  end subroutine barrou


!   subroutine excj( a, sx, ex, sy, ey, sz,ez)
!     implicit none
!
!   integer sx, ex, sy, ey, sz, ez
!   real a(sx:ex, sy:ey, sz:ez)
!   integer iiget, status(MPI_STATUS_SIZE)
!   integer ii, i, k
!   real,allocatable, dimension(:) :: buffj1,buffj2,buffj3,buffj4
!   iiget = (ex - sx + 1)*(ez - sz + 1)
!
!   allocate( buffj1(iiget),&
!             buffj2(iiget),&
!             buffj3(iiget),&
!             buffj4(iiget))
!
!
!
!   if(nbrnorth/=MPI_PROC_NULL)then
!     do k=sz,ez
!     do i=sx,ex
!       ii = i - sx + 1 + (k - sz )*(ex - sx + 1)
!       buffj1(ii) = a(i,ey-1,k)
!     enddo
!     enddo
!   endif
!   call MPI_SENDRECV(  buffj1,  ii    , MY_REAL, nbrnorth, 4, &
!                       buffj2,  iiget , MY_REAL, nbrsouth,  4,      &
!                       comm3d, status, mpierr )
!   if(nbrsouth/=MPI_PROC_NULL)then
!     do k=sz,ez
!     do i=sx,ex
!       ii = i - sx + 1 + (k - sz )*(ex - sx + 1)
!       a(i,sy,k) = buffj2(ii)
!     enddo
!     enddo
!   endif
!
! !   call barrou()
!
!   if(nbrsouth/=MPI_PROC_NULL)then
!     do k=sz,ez
!     do i=sx,ex
!       ii = i - sx + 1 + (k - sz )*(ex - sx + 1)
!       buffj3(ii) = a(i,sy+1,k)
!     enddo
!     enddo
!   endif
!   call MPI_SENDRECV(  buffj3,  ii    , MY_REAL, nbrnorth,  5, &
!                           buffj4,  iiget , MY_REAL, nbrsouth, 5, &
!                           comm3d, status, mpierr )
!   if(nbrnorth/=MPI_PROC_NULL)then
!     do k=sz,ez
!     do i=sx,ex
!       ii = i - sx + 1 + (k - sz )*(ex - sx + 1)
!       a(i,ey,k) = buffj4(ii)
!     enddo
!     enddo
!   endif
!
! !   call barrou()
!   deallocate (buffj1,buffj2,buffj3,buffj4)
!
!   return
!   end subroutine excj

!   ! tg3315 - does some transfer of ghost cell information in the y-direction.
!   ! nbrbottom is cpu below and nbrtop is cpu above
!   ! MPI_PROC_NULL = -2???
!   subroutine excjs( a, sx, ex, sy, ey, sz,ez,ih,jh)
!     implicit none
!   integer sx, ex, sy, ey, sz, ez, ih, jh
!   real a(sx-ih:ex+ih, sy-jh:ey+jh, sz:ez)
!   integer status(MPI_STATUS_SIZE), iiget
!   integer ii, i, j, k
!   real,allocatable, dimension(:) :: buffj1,buffj2,buffj3,buffj4
!   iiget = jh*(ex - sx + 1 + 2*ih)*(ez - sz + 1)
!
!   allocate( buffj1(iiget),&
!             buffj2(iiget),&
!             buffj3(iiget),&
!             buffj4(iiget))
!
!   if(nbrtop/=MPI_PROC_NULL)then
!     ii = 0
!     do j=1,jh
!     do k=sz,ez
!     do i=sx-ih,ex+ih
!       ii = ii + 1
!       buffj1(ii) = a(i,ey-j+1,k) ! tg3315 buffj1 is je-jhc ghost cells on myid
!     enddo
!     enddo
!     enddo
!   endif
!
!   call MPI_SENDRECV(  buffj1,  ii    , MY_REAL, nbrtop, 4, &
!                            buffj2,  iiget , MY_REAL, nbrbottom,  4, &
!                            comm3d, status, mpierr )
!
!   ! tg3315 sends this to nbrtop and pulls buffj2 from nbrbottom! send and receive process that is good for executing chain shifts.
!
!   if(nbrbottom/=MPI_PROC_NULL)then
!     ii = 0
!     do j=1,jh
!     do k=sz,ez
!     do i=sx-ih,ex+ih
!       ii = ii + 1
!       a(i,sy-j,k) = buffj2(ii) ! set the previous ghost cells to buffj2 (last cells of nbrbottom I think)
!     enddo
!     enddo
!     enddo
!   endif
!
! !   call barrou()
!
!   ! repeats this process for other way round
!
!   if(nbrbottom/=MPI_PROC_NULL)then
!     ii = 0
!     do j=1,jh
!     do k=sz,ez
!     do i=sx-ih,ex+ih
!       ii = ii + 1
!       buffj3(ii) = a(i,sy+j-1,k)
!     enddo
!     enddo
!     enddo
!   endif
!   call MPI_SENDRECV(  buffj3,  ii    , MY_REAL, nbrbottom,  5, &
!                           buffj4,  iiget , MY_REAL, nbrtop, 5, &
!                           comm3d, status, mpierr )
!   if(nbrtop/=MPI_PROC_NULL)then
!     ii = 0
!     do j=1,jh
!     do k=sz,ez
!     do i=sx-ih,ex+ih
!       ii = ii + 1
!       a(i,ey+j,k) = buffj4(ii)
!     enddo
!     enddo
!     enddo
!   endif
!
! !   call barrou()
!   deallocate (buffj1,buffj2,buffj3,buffj4)
!
!   return
!   end subroutine excjs

! subroutine exci(a,sx,ex,sy,ey,sz,ez)
!   implicit none
!   integer sx, ex, sy, ey, sz, ez
!   real a(sx:ex, sy:ey, sz:ez)
!   integer status(MPI_STATUS_SIZE)
!   integer ii, i, j, k
!   integer reqe, reqw
!   integer ewsize
!   real,allocatable, dimension(:) :: sende,recve
!   real,allocatable, dimension(:) :: sendw,recvw

! !   Calculate buffer size
!   ewsize = (ey - sy + 1)*(ez - sz + 1)

! !   Allocate send / receive buffers
!   allocate(sende(ewsize),sendw(ewsize))
!   allocate(recve(ewsize),recvw(ewsize))

!   if(nprocx .gt. 1)then
!     !   Send east/west
!     ii = 0
!     do k=sz,ez
!     do j=sy,ey
!       ii = ii + 1
!       sende(ii) = a(ex-i+1,j,k)
!       sendw(ii) = a(sx+i-1,j,k)
!     enddo
!     enddo

!     call MPI_ISEND(sende, ewsize, MY_REAL, nbreast, 6, comm3d, reqe, mpierr)
!     call MPI_ISEND(sendw, ewsize, MY_REAL, nbrwest, 7, comm3d, reqw, mpierr)

!     !   Receive west/east
!     call MPI_RECV(recvw, ewsize, MY_REAL, nbrwest, 6, comm3d, status, mpierr)
!     call MPI_RECV(recve, ewsize, MY_REAL, nbreast, 7, comm3d, status, mpierr)

!     ii = 0
!     do k=sz,ez
!     do j=sy,ey
!       ii = ii + 1
!       a(sx-i,j,k) = recvw(ii)
!       a(ex+i,j,k) = recve(ii)
!     enddo
!     enddo
!   else
!     ! Single processor, make sure the field is periodic
!     do k=sz,ez
!     do j=sy,ey
!       a(sx-i,j,k) = a(ex-i+1,j,k)
!       a(ex+i,j,k) = a(sx+i-1,j,k)
!     enddo
!     enddo
!   endif

!   if(nprocx.gt.1)then
!     call MPI_WAIT(reqe, status, mpierr)
!     call MPI_WAIT(reqw, status, mpierr)
!   endif

!   deallocate (sende, sendw)
!   deallocate (recve, recvw)

!   return
! end subroutine exci

subroutine excis(a,sx,ex,sy,ey,sz,ez,ih,jh)
  implicit none
  integer sx, ex, sy, ey, sz, ez, ih, jh
  real a(sx-ih:ex+ih, sy-jh:ey+jh, sz:ez)
  integer status(MPI_STATUS_SIZE)
  integer ii, i, j, k
  integer reqe, reqw
  integer ewsize
  real,allocatable, dimension(:) :: sende,recve
  real,allocatable, dimension(:) :: sendw,recvw

!   Calculate buffer size
  ewsize = ih*(ey - sy + 1 + 2*jh)*(ez - sz + 1)

!   Allocate send / receive buffers
  allocate(sende(ewsize),sendw(ewsize))
  allocate(recve(ewsize),recvw(ewsize))

  if(nprocx .gt. 1)then
    !   Send east/west
    ii = 0
    do i=1,ih
    do k=sz,ez
    do j=sy-jh,ey+jh
      ii = ii + 1
      sende(ii) = a(ex-i+1,j,k)
      sendw(ii) = a(sx+i-1,j,k)
    enddo
    enddo
    enddo

    call MPI_ISEND(sende, ewsize, MY_REAL, nbreast, 6, comm3d, reqe, mpierr)
    call MPI_ISEND(sendw, ewsize, MY_REAL, nbrwest, 7, comm3d, reqw, mpierr)

    !   Receive west/east
    call MPI_RECV(recvw, ewsize, MY_REAL, nbrwest, 6, comm3d, status, mpierr)
    call MPI_RECV(recve, ewsize, MY_REAL, nbreast, 7, comm3d, status, mpierr)

    ii = 0
    do i=1,ih
    do k=sz,ez
    do j=sy-jh,ey+jh
      ii = ii + 1
      a(sx-i,j,k) = recvw(ii)
      a(ex+i,j,k) = recve(ii)
    enddo
    enddo
    enddo
  else
    ! Single processor, make sure the field is periodic
    do i=1,ih
    do k=sz,ez
    do j=sy-jh,ey+jh
      a(sx-i,j,k) = a(ex-i+1,j,k)
      a(ex+i,j,k) = a(sx+i-1,j,k)
    enddo
    enddo
    enddo
  endif

  if(nprocx.gt.1)then
    call MPI_WAIT(reqe, status, mpierr)
    call MPI_WAIT(reqw, status, mpierr)
  endif

  deallocate (sende, sendw)
  deallocate (recve, recvw)

  return
end subroutine excis

! subroutine excj(a,sx,ex,sy,ey,sz,ez)
!   implicit none
!   integer sx, ex, sy, ey, sz, ez
!   real a(sx:ex, sy:ey, sz:ez)
!   integer status(MPI_STATUS_SIZE)
!   integer ii, i, j, k
!   integer reqn, reqs
!   integer nssize
!   real,allocatable, dimension(:) :: sendn,recvn
!   real,allocatable, dimension(:) :: sends,recvs

! !   Calculate buffer size
!   nssize = (ex - sx + 1)*(ez - sz + 1)

! !   Allocate send / receive buffers
!   allocate(sendn(nssize),sends(nssize))
!   allocate(recvn(nssize),recvs(nssize))

!   if(nprocy .gt. 1)then
!     !   Send north/south
!     ii = 0
!     do k=sz,ez
!     do i=sx,ex
!       ii = ii + 1
!       sendn(ii) = a(i,ey-j+1,k)
!       sends(ii) = a(i,sy+j-1,k)
!     enddo
!     enddo

!     call MPI_ISEND(sendn, nssize, MY_REAL, nbrnorth, 4, comm3d, reqn, mpierr)
!     call MPI_ISEND(sends, nssize, MY_REAL, nbrsouth, 5, comm3d, reqs, mpierr)

!     !   Receive south/north
!     call MPI_RECV(recvs, nssize, MY_REAL, nbrsouth, 4, comm3d, status, mpierr)
!     call MPI_RECV(recvn, nssize, MY_REAL, nbrnorth, 5, comm3d, status, mpierr)

!     ii = 0
!     do k=sz,ez
!     do i=sx,ex
!       ii = ii + 1
!       a(i,sy-j,k) = recvs(ii)
!       a(i,ey+j,k) = recvn(ii)
!     enddo
!     enddo
!   else
!     ! Single processor, make sure the field is periodic
!     do k=sz,ez
!     do i=sx,ex
!       a(i,sy-j,k) = a(i,ey-j+1,k)
!       a(i,ey+j,k) = a(i,sy+j-1,k)
!     enddo
!     enddo
!   endif

!   if(nprocy.gt.1)then
!     call MPI_WAIT(reqn, status, mpierr)
!     call MPI_WAIT(reqs, status, mpierr)
!   endif

!   deallocate (sendn, sends)
!   deallocate (recvn, recvs)

!   return
!   end subroutine excj


subroutine excjs(a,sx,ex,sy,ey,sz,ez,ih,jh)
  implicit none
  integer sx, ex, sy, ey, sz, ez, ih, jh
  real a(sx-ih:ex+ih, sy-jh:ey+jh, sz:ez)
  integer status(MPI_STATUS_SIZE)
  integer ii, i, j, k
  integer reqn, reqs
  integer nssize
  real,allocatable, dimension(:) :: sendn,recvn
  real,allocatable, dimension(:) :: sends,recvs

!   Calculate buffer size
  nssize = jh*(ex - sx + 1 + 2*ih)*(ez - sz + 1)

!   Allocate send / receive buffers
  allocate(sendn(nssize),sends(nssize))
  allocate(recvn(nssize),recvs(nssize))

  if(nprocy .gt. 1)then
    !   Send north/south
    ii = 0
    do j=1,jh
    do k=sz,ez
    do i=sx-ih,ex+ih
      ii = ii + 1
      sendn(ii) = a(i,ey-j+1,k)
      sends(ii) = a(i,sy+j-1,k)
    enddo
    enddo
    enddo

    call MPI_ISEND(sendn, nssize, MY_REAL, nbrnorth, 4, comm3d, reqn, mpierr)
    call MPI_ISEND(sends, nssize, MY_REAL, nbrsouth, 5, comm3d, reqs, mpierr)

    !   Receive south/north
    call MPI_RECV(recvs, nssize, MY_REAL, nbrsouth, 4, comm3d, status, mpierr)
    call MPI_RECV(recvn, nssize, MY_REAL, nbrnorth, 5, comm3d, status, mpierr)

    ii = 0
    do j=1,jh
    do k=sz,ez
    do i=sx-ih,ex+ih
      ii = ii + 1
      a(i,sy-j,k) = recvs(ii)
      a(i,ey+j,k) = recvn(ii)
    enddo
    enddo
    enddo
  else
    ! Single processor, make sure the field is periodic
    do j=1,jh
    do k=sz,ez
    do i=sx-ih,ex+ih
      a(i,sy-j,k) = a(i,ey-j+1,k)
      a(i,ey+j,k) = a(i,sy+j-1,k)
    enddo
    enddo
    enddo
  endif

  if(nprocy.gt.1)then
    call MPI_WAIT(reqn, status, mpierr)
    call MPI_WAIT(reqs, status, mpierr)
  endif

  deallocate (sendn, sends)
  deallocate (recvn, recvs)

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

  subroutine avexy_ibm(aver,var,ib,ie,jb,je,kb,ke,ih,jh,kh,II,IIs,lnan)
    implicit none

    integer :: ib,ie,jb,je,kb,ke,ih,jh,kh
    real    :: aver(kb:ke+kh)
    real    :: var(ib:ie,jb:je,kb:ke+kh)
    integer :: II(ib:ie,jb:je,kb:ke+kh)
    integer :: IIs(kb:ke+kh)
    integer :: IId(kb:ke+kh)
    real    :: averl(kb:ke+kh)
    real    :: avers(kb:ke+kh)
    integer :: k
    logical :: lnan

    averl       = 0.
    avers       = 0.

    do k=kb,ke+kh
      averl(k) = sum(var(ib:ie,jb:je,k)*II(ib:ie,jb:je,k))
    enddo

    IId = IIs

    ! tg3315 22.03.19 - if not calculating stats and all blocks on lowest layer...
    ! should not be necessary but value at kb is used in modthermo so reasonable value must
    ! be assigned. Potentially should leave as before and only account for in modthermo...
    if ((.not. lnan) .and. (IId(kb)==0)) then
      averl(kb) = sum(var(ib:ie,jb:je,kb))
      IId(kb) = IId(ke)
    end if

    call MPI_ALLREDUCE(averl, avers, ke+kh-kb+1,  MY_REAL, &
                          MPI_SUM, comm3d,mpierr)

    where (IId==0)
      aver = -999.
    elsewhere
      aver = avers/IId
    endwhere

    return
  end subroutine avexy_ibm

  subroutine slabsumi(aver,iis,iif,var,ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes)
    implicit none

    integer :: iis,iif
    integer :: ib,ie,jb,je,kb,ke,ibs,ies,jbs,jes,kbs,kes
    real    :: aver(iis:iif)
    real    :: var (ib:ie,jb:je,kb:ke)
    real    :: averl(iis:iif)
    real    :: avers(iis:iif)
    integer :: i

    averl       = 0.
    avers       = 0.

    do i=ibs,ies
      averl(i) = sum(var(i,jbs:jes,kbs:kes))
    enddo

    call MPI_ALLREDUCE(averl, avers, iif-iis+1,  MY_REAL, &
                          MPI_SUM, comm3d,mpierr)

    aver = aver + avers

    return
  end subroutine slabsumi

  !Could make this so it has cases like if the input is 1,2 or 3D...
  subroutine avey_ibm(aver,var,ib,ie,jb,je,kb,ke,II,IIt)
  implicit none
  integer                 :: ib,ie,jb,je,kb,ke
  real                    :: aver(ib:ie,kb:ke)
  real                    :: avero(ib:ie,kb:ke)
  real                    :: var(ib:ie,jb:je,kb:ke)
  integer                 :: II(ib:ie,jb:je,kb:ke)
  integer                 :: IIt(ib:ie,kb:ke)
  logical                 :: lytdump,lnan

  avero = 0.
  aver  = 0.

  avero = sum(var(ib:ie,jb:je,kb:ke)*II(ib:ie,jb:je,kb:ke), DIM=2)

  call MPI_ALLREDUCE(avero(ib:ie,kb:ke), aver(ib:ie,kb:ke), (ke-kb+1)*(ie-ib+1), MY_REAL,MPI_SUM, comm3d,mpierr)

  where (IIt==0)
    aver = -999.
  elsewhere
    aver = aver/IIt
  endwhere

  end subroutine avey_ibm

  subroutine sumy_ibm(sumy,var,ib,ie,jb,je,kb,ke,II)
    ! This routine sums up a variable over the y direction,
    ! only including the fluid cells.
    implicit none
    integer                 :: ib,ie,jb,je,kb,ke
    real                    :: sumy(ib:ie,kb:ke)
    real                    :: sumproc(ib:ie,kb:ke)
    real                    :: var(ib:ie,jb:je,kb:ke)
    integer                 :: II(ib:ie,jb:je,kb:ke)

    sumproc = 0.
    sumy  = 0.

    sumproc = sum(var(ib:ie,jb:je,kb:ke)*II(ib:ie,jb:je,kb:ke), DIM=2)

    call MPI_ALLREDUCE(sumproc(ib:ie,kb:ke), sumy(ib:ie,kb:ke), (ke-kb+1)*(ie-ib+1), MY_REAL,MPI_SUM, comm3d,mpierr)

    end subroutine sumy_ibm


    subroutine sumx_ibm(sumx,var,ib,ie,jb,je,kb,ke,II)
     ! This routine sums up a variable over the x direction,
     ! only including the fluid cells.
     implicit none
     integer                 :: ib,ie,jb,je,kb,ke
     real                    :: sumx(jb:je,kb:ke)
     real                    :: sumproc(jb:je,kb:ke)
     real                    :: var(ib:ie,jb:je,kb:ke)
     integer                 :: II(ib:ie,jb:je,kb:ke)

     sumproc = 0.
     sumx  = 0.

     sumproc = sum(var(ib:ie,jb:je,kb:ke)*II(ib:ie,jb:je,kb:ke), DIM=1)

     call MPI_ALLREDUCE(sumproc(jb:je,kb:ke), sumx(jb:je,kb:ke), (ke-kb+1)*(je-jb+1), MY_REAL,MPI_SUM, comm3d,mpierr)

     end subroutine sumx_ibm

end module
