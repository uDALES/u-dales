!> \file modchecksim.f90
!!  Monitors Courant and Peclet numbers, and divergence.

!>
!!  Monitors Courant and Peclet numbers, and divergence.
!>
!!  These numbers are put out to screen either every tcheck seconds, or every time step (if tcheck=0).
!!  \autor Jasper Tomas, TU Delft, June 4th 2015
!!  \author Thijs Heus,MPI-M
!!  \author Hans Cuijpers, KNMI
!!  \par Revision list
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
!
module modchecksim
  use modglobal, only : longint
#if defined(_GPU)
  use cudafor
#endif

  implicit none
  private
  public initchecksim,checksim
  ! The three rank-local reductions are the only part of this module that has
  ! both a host and a device implementation, so the tests drive them directly:
  ! going through checksim would mean reproducing its timekeeping state and
  ! reading numbers back off a screen line. diffnrgeom is exposed for the same
  ! reason - tests.f90 asserts the cache against the expression it replaced.
  ! Public on every build rather than behind the usual test guard, because
  ! tests.f90 compiles its host-branch checks in a GPU build too, where they
  ! only report that the runmode belongs on a CPU binary.
  public :: courant_local, diffnr_local, div_local, diffnrgeom
#if defined(_GPU) && defined(UDALES_DEBUG)
  ! The device mirror as well, so tests_cuda.f90 can substitute a cache whose
  ! value differs at every (i,k). The grids the GPU cases run on are uniform,
  ! which makes the real cache constant and a kernel that reads it at a fixed
  ! index indistinguishable from one that indexes it properly. Test-only, so
  ! this one stays behind the guard.
  public :: diffnrgeom_d
#endif

  real    :: tcheck = 0.
  !integer(kind=longint) :: tnext = 3600.,itcheck
  real    :: tnext = 0.
  real    :: dtmn =0.,ndt =0.

  !> Geometry factor of the diffusion number, 1/dzh(k)**2 + dxh2i(i) + dy2i.
  !!
  !! It depends on the grid alone, so it is built once in initchecksim rather
  !! than per cell per call: on a 256^3 grid the expression it replaces is 16.7
  !! million double-precision divisions, three adds each, every time checksim
  !! fires. The terms are combined in the order the original expression used,
  !! so the cached value is bit-identical to it and the printed diffusion
  !! number does not move.
  !!
  !! It is derived once for both paths - the host loop reads it and the device
  !! mirror is a straight copy - so the two cannot drift apart, and the device
  !! kernel needs no mirror of dxh2i.
  !!
  !! The i extent is the loop's own, ib:ie, so it is indexed exactly as dxh2i
  !! was. It costs (ie-ib+1)*(ke-kb+1) reals: half a megabyte on that 256^3
  !! grid, against the 269 megabytes of ekm and ekh the same loop streams.
#if defined(_GPU)
  real, allocatable, pinned :: diffnrgeom(:,:)
  real, device, allocatable :: diffnrgeom_d(:,:)
#else
  real, allocatable :: diffnrgeom(:,:)
#endif

  save
contains
!> Initializing Checksim. Read out the namelist, initializing the variables
  subroutine initchecksim
    use modglobal, only : ifnamopt, fname_options,dtmax,ladaptive,btime, &
                          ib,ie,kb,ke,dxh2i,dy2i,dzh
    use modmpi,    only : myid,my_real,comm3d,mpierr
    implicit none
    integer :: ierr
    integer :: i, k
    namelist/NAMCHECKSIM/ &
    tcheck

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMCHECKSIM,iostat=ierr)
      if (ierr > 0) then
        write(0, *) 'ERROR: Problem in namoptions NAMCHECKSIM'
        write(0, *) 'iostat error: ', ierr
        stop 1
      endif
      !write(6 ,NAMCHECKSIM)
      close(ifnamopt)

      if ((.not. ladaptive) .and. (tcheck < dtmax)) then
        tcheck = dtmax
      end if
    end if

    call MPI_BCAST(tcheck     ,1,MY_REAL   ,0,comm3d,mpierr)
!    itcheck = floor(tcheck/tres)
    tnext = tcheck+btime

    ! The diffusion-number geometry, once. See the declaration for why. The
    ! cell volume the divergence integral needs is the same kind of quantity
    ! but a generic one, so it sits with the other grid-only arrays in
    ! modglobal as dvcell, next to the reciprocal dxdydzfi.
    if (.not. allocated(diffnrgeom)) allocate(diffnrgeom(ib:ie,kb:ke))
    do k=kb,ke
      do i=ib,ie
        diffnrgeom(i,k) = 1/dzh(k)**2 + dxh2i(i) + dy2i
      end do
    end do
#if defined(_GPU)
    if (.not. allocated(diffnrgeom_d)) allocate(diffnrgeom_d(ib:ie,kb:ke))
    diffnrgeom_d = diffnrgeom
#endif

  end subroutine initchecksim
!>Run checksim. Timekeeping, and output
  subroutine checksim
    use modglobal, only : timee, rk3step, dt
    use modmpi,    only : myid
    implicit none
    character(20) :: timeday
    if (timee ==0.0) return
    if (rk3step/=3) return
    dtmn = dtmn +dt; ndt =ndt+1.
    if(timee<tnext) return
    tnext = tnext+tcheck
    dtmn  = dtmn / ndt
    if (myid==0) then
      call date_and_time(time=timeday)
      write (*,*) '================================================================='
      write (*,'(3A,F15.5,A,F12.9)') 'Time of Day: ', timeday(1:10),'    Time of Simulation: ', timee, '    dt: ',dtmn
    end if
    call calccourant
    call calcdiffnr
    call chkdiv
    dtmn  = 0.
    ndt   = 0.

  end subroutine checksim
!>      Calculates the courant number as in max(w)*deltat/deltaz
  subroutine calccourant
    use modmpi,    only : myid,comm3d,mpierr,mpi_max,my_real
    implicit none

    real          :: courtotl, courtot

    courtot  = 0.0
    call courant_local(dtmn, courtotl)

    call MPI_ALLREDUCE(courtotl,courtot,1,MY_REAL,MPI_MAX,comm3d,mpierr)
    if (myid==0) then
      write(*,'(A,ES10.2)') 'Courant numbers (x,y,z,tot):',courtot
    end if

    return
  end subroutine calccourant

!> Rank-local maximum Courant number over the cells this rank owns.
!!
!! Split out of calccourant so that the loop, which is the only part with a
!! host and a device form, can be driven by the tests without an MPI
!! communicator, and so that the reduction it feeds stays in one place.
!! dtm is passed rather than read from dtmn so the result depends on nothing
!! but its arguments.
!!
!! dtm multiplies the maximum rather than every cell, which is one multiply
!! instead of one per cell and gives the same answer to the last bit: rounding
!! is monotonic and dtm is a mean time step, so it is never negative, and
!! max(fl(a*dtm), fl(b*dtm)) is fl(max(a,b)*dtm). The zero the reduction starts
!! from survives the move too, since 0*dtm is 0.
  subroutine courant_local(dtm, courtotl)
    use modglobal, only : ib,ie,jb,je,kb,ke,dyi
#if defined(_GPU)
    use modcuda,   only : um_d,vm_d,wm_d,dxhi_d,dzhi_d
#else
    use modglobal, only : dxhi,dzhi
    use modfields, only : um,vm,wm
#endif
    implicit none

    real, intent(in)  :: dtm
    real, intent(out) :: courtotl
    integer           :: i, j, k

    courtotl = 0.0
#if defined(_GPU)
    !$acc parallel loop collapse(3) default(present) reduction(max:courtotl)
    do k=kb,ke
      do j=jb,je
        do i=ib,ie
          courtotl = max(courtotl, um_d(i,j,k)*dxhi_d(i) + vm_d(i,j,k)*dyi + wm_d(i,j,k)*dzhi_d(k))
        end do
      end do
    end do
    !$acc end parallel loop
#else
    do k=kb,ke
      do j=jb,je
        do i=ib,ie
          courtotl = max(courtotl, um(i,j,k)*dxhi(i) + vm(i,j,k)*dyi + wm(i,j,k)*dzhi(k))
        end do
      end do
    end do
#endif
    courtotl = courtotl*dtm

  end subroutine courant_local

!> Calculates the diffusion number as max(ekm) *deltat/deltax**2
  subroutine calcdiffnr
    use modmpi,         only : myid,comm3d,mpierr,mpi_max,my_real
    implicit none

    real diffnrtotl,diffnrtot

    diffnrtot  = 0.
    call diffnr_local(dtmn, diffnrtotl)

    call MPI_ALLREDUCE(diffnrtotl,diffnrtot,1,MY_REAL,MPI_MAX,comm3d,mpierr)
    if (myid==0) then
      write(6,'(A,ES10.2)') 'Diffusion number:',diffnrtot
    end if

    return
  end subroutine calcdiffnr

!> Rank-local maximum diffusion number over the cells this rank owns.
!!
!! diffnrgeom(i,k) is 1/dzh(k)**2 + dxh2i(i) + dy2i, built once in
!! initchecksim; see its declaration. ekh is taken alongside ekm because the
!! stability limit is set by whichever diffusivity is larger, and on a run with
!! a Prandtl number below one that is ekh.
!!
!! dtm multiplies the maximum rather than every cell, for the reason given in
!! courant_local.
  subroutine diffnr_local(dtm, diffnrtotl)
    use modglobal,      only : ib,ie,jb,je,kb,ke
#if defined(_GPU)
    use modcuda,        only : ekm_d,ekh_d
#else
    use modsubgriddata, only : ekm,ekh
#endif
    implicit none

    real, intent(in)  :: dtm
    real, intent(out) :: diffnrtotl
    integer           :: i, j, k

    diffnrtotl = 0.
#if defined(_GPU)
    !$acc parallel loop collapse(3) default(present) reduction(max:diffnrtotl)
    do k=kb,ke
    do j=jb,je
    do i=ib,ie
      diffnrtotl = max(diffnrtotl,  ekm_d(i,j,k)*diffnrgeom_d(i,k), &
                                    ekh_d(i,j,k)*diffnrgeom_d(i,k) )  ! or should I interpolate ekm to the correct position?
    end do
    end do
    end do
    !$acc end parallel loop
#else
    do k=kb,ke
    do j=jb,je
    do i=ib,ie
      diffnrtotl = max(diffnrtotl,  ekm(i,j,k)*diffnrgeom(i,k), &
                                    ekh(i,j,k)*diffnrgeom(i,k) )  ! or should I interpolate ekm to the correct position?
    end do
    end do
    end do
#endif
    diffnrtotl = diffnrtotl*dtm

  end subroutine diffnr_local

!> Checks local and total divergence
  subroutine chkdiv
    use modmpi,    only : myid,comm3d,mpi_sum,mpi_max,my_real,mpierr
    implicit none

    real divmax, divtot
    real divmaxl, divtotl

    divmax = 0.
    divtot = 0.
    call div_local(divmaxl, divtotl)

    call MPI_ALLREDUCE(divtotl, divtot, 1,    MY_REAL, &
                          MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(divmaxl, divmax, 1,    MY_REAL, &
                          MPI_MAX, comm3d,mpierr)

    if(myid==0)then
      write(6,'(A,2ES11.2)')'divmax, divtot = ', divmax, divtot
    end if

    return
  end subroutine chkdiv

!> Rank-local peak and volume-integrated divergence over the cells this rank owns.
!!
!! The i+1, j+1 and k+1 neighbours reach one cell past ie, je and ke, so the
!! halo has to be current when this runs - which after the halo exchange at the
!! end of the step it is.
!!
!! dvcell(k) is the cell volume dx*dy*dzf(k), built once in initglobal beside
!! the reciprocal the tendency routines use; see its declaration in modglobal.
!! Unlike diffnrgeom it is not bit-identical to the expression it replaces -
!! div*dx*dy*dzf(k) associates left to right, so folding the three grid factors
!! into one operand rounds differently. That costs nothing here: divtot is a
!! printed diagnostic, compared nowhere, and a sum that very nearly cancels. On
!! a 256^3 run two invocations of the same binary already disagree about every
!! one of its printed values, because the wall-function atomics leave the
!! fields themselves differing in their last bits.
  subroutine div_local(divmaxl, divtotl)
    use modglobal, only : ib,ie,jb,je,ke,kb,dxi,dyi
#if defined(_GPU)
    use modcuda,   only : u0_d,v0_d,w0_d,dzfi_d,dvcell_d
#else
    use modglobal, only : dzfi,dvcell
    use modfields, only : u0,v0,w0!,divergentie
#endif
    implicit none

    real, intent(out) :: divmaxl, divtotl
    real    :: div
    integer :: i, j, k

    divmaxl= 0.
    divtotl= 0.
#if defined(_GPU)
    !$acc parallel loop collapse(3) default(present) private(div) &
    !$acc          reduction(max:divmaxl) reduction(+:divtotl)
    do k=kb,ke
    do j=jb,je
    do i=ib,ie
      div = &
                (u0_d(i+1,j,k) - u0_d(i,j,k) )*dxi + &
                (v0_d(i,j+1,k) - v0_d(i,j,k) )*dyi + &
                (w0_d(i,j,k+1) - w0_d(i,j,k) )*dzfi_d(k)
      divmaxl = max(divmaxl,abs(div))
      divtotl = divtotl + div*dvcell_d(k)
    end do
    end do
    end do
    !$acc end parallel loop
#else
    do k=kb,ke
    do j=jb,je
    do i=ib,ie
      div = &
                (u0(i+1,j,k) - u0(i,j,k) )*dxi + &
                (v0(i,j+1,k) - v0(i,j,k) )*dyi + &
                (w0(i,j,k+1) - w0(i,j,k) )*dzfi(k)
!      divergentie(i,j,k)=div
      divmaxl = max(divmaxl,abs(div))
      divtotl = divtotl + div*dvcell(k)
    end do
    end do
    end do
#endif

  end subroutine div_local

end module modchecksim
