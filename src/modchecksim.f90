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
  use modglobal, only : longint, tcheck

  implicit none
  private
  public initchecksim,checksim

  !integer(kind=longint) :: tnext = 3600.,itcheck
  real    :: tnext = 0.
  real    :: dtmn =0.,ndt =0.

  save
contains
!> Initializing Checksim. Read out the namelist, initializing the variables
  subroutine initchecksim
    use modglobal, only : dtmax,ladaptive,btime
    implicit none

    if ((.not. ladaptive) .and. (tcheck < dtmax)) then
      tcheck = dtmax
    end if

    tnext = tcheck+btime


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
    use mpi
    use modglobal, only : ib,ie,jb,je,kb,ke,dxhi,dyi,dzhi
    use modfields, only : um,vm,wm,thl0
    use modmpi,    only : myid,comm3d,mpierr,mpi_max,my_real
    use decomp_2d, only : zstart
    implicit none


    real          :: courtotl, courtot, cour
    real          :: maxl(2), maxg(2)
    integer       :: i, j, k, im, jm, km

    courtotl = 0.0
    courtot  = 0.0
    im = ib; jm = jb; km = kb
    do k=kb,ke
      do j=jb,je
        do i=ib,ie
          cour = (um(i,j,k)*dxhi(i) + vm(i,j,k)*dyi + wm(i,j,k)*dzhi(k))*dtmn
          if (cour > courtotl) then
            courtotl = cour
            im = i; jm = j; km = k
          end if
        end do
      end do
    end do

    call MPI_ALLREDUCE(courtotl,courtot,1,MY_REAL,MPI_MAX,comm3d,mpierr)
    if (myid==0) then
      write(*,'(A,ES10.2)') 'Courant numbers (x,y,z,tot):',courtot
    end if

    ! Diagnostic: report where the domain-wide Courant maximum lives.
    maxl(1) = courtotl
    maxl(2) = real(myid)
    call MPI_ALLREDUCE(maxl,maxg,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,comm3d,mpierr)
    if (myid == nint(maxg(2))) then
      write(*,'(A,I6,A,3I5,A,4ES11.3)') 'maxcour: rank', myid, &
        '  glob(i,j,k)=', im+zstart(1)-1, jm+zstart(2)-1, km+zstart(3)-1, &
        '  u,v,w,thl=', um(im,jm,km), vm(im,jm,km), wm(im,jm,km), thl0(im,jm,km)
    end if

    return
  end subroutine calccourant

!> Calculates the diffusion number as max(ekm) *deltat/deltax**2
  subroutine calcdiffnr

    use modglobal,      only : ib,ie,jb,je,kb,ke,dxh2i,dy2i,dzh
    use modsubgriddata, only : ekm,ekh
    use modmpi,         only : myid,comm3d,mpierr,mpi_max,my_real
    implicit none


    real diffnrtotl,diffnrtot
    integer       :: i,j,k

    diffnrtotl = 0.
    diffnrtot  = 0.
    do k=kb,ke
    do j=jb,je
    do i=ib,ie
!      diffnrtotl = max(diffnrtotl,  ekm(i,j,k)*(1/dzh(k)**2 + dxh2i(i) + dy2i)*dtmn )  ! or should I interpolate ekm to the correct position?
      diffnrtotl = max(diffnrtotl,  ekm(i,j,k)*(1/dzh(k)**2 + dxh2i(i) + dy2i)*dtmn, &
                                    ekh(i,j,k)*(1/dzh(k)**2 + dxh2i(i) + dy2i)*dtmn )  ! or should I interpolate ekm to the correct position?
    end do
    end do
    end do

    call MPI_ALLREDUCE(diffnrtotl,diffnrtot,1,MY_REAL,MPI_MAX,comm3d,mpierr)
    if (myid==0) then
      write(6,'(A,ES10.2)') 'Diffusion number:',diffnrtot
    end if

    return
  end subroutine calcdiffnr

!> Checks local and total divergence
  subroutine chkdiv

    use modglobal, only : ib,ie,jb,je,ke,kb,dx,dxi,dy,dyi,dzf,dzfi
    use modfields, only : u0,v0,w0!,divergentie
    use modmpi,    only : myid,comm3d,mpi_sum,mpi_max,my_real,mpierr
    implicit none



    real div, divmax, divtot
    real divmaxl, divtotl
    integer i, j, k

    divmax = 0.
    divtot = 0.
    divmaxl= 0.
    divtotl= 0.

    do k=kb,ke
    do j=jb,je
    do i=ib,ie
      div = &
                (u0(i+1,j,k) - u0(i,j,k) )*dxi + &
                (v0(i,j+1,k) - v0(i,j,k) )*dyi + &
                (w0(i,j,k+1) - w0(i,j,k) )*dzfi(k)
!      divergentie(i,j,k)=div
      divmaxl = max(divmaxl,abs(div))
      divtotl = divtotl + div*dx*dy*dzf(k)
    end do
    end do
    end do

    call MPI_ALLREDUCE(divtotl, divtot, 1,    MY_REAL, &
                          MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(divmaxl, divmax, 1,    MY_REAL, &
                          MPI_MAX, comm3d,mpierr)

    if(myid==0)then
      write(6,'(A,2ES11.2)')'divmax, divtot = ', divmax, divtot
    end if

    return
  end subroutine chkdiv

end module modchecksim
