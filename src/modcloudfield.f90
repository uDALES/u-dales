!> \file modcloudfield.f90
!!  Dumps all the wet points in the field

!>
!! Dumps all the wet points in the field to clouds.myid.expnr
!>
!!  \author Harm Jonker, TU Delft
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

module modcloudfield


implicit none
private
PUBLIC :: initcloudfield, cloudfield
save

  real    :: dtav,tnext
  logical :: lcloudfield= .false. !< switch for writing cloud field (on/off)
  logical :: laddinfo   = .false. !< switch to write ql and w values (on/off)

contains
!> Initializing Cloudfield. Read out the namelist, initializing the variables
  subroutine initcloudfield

    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_logical
    use modglobal,only :ifnamopt,fname_options,dtmax,dtav_glob,dt_lim,btime
    implicit none
    integer :: ierr
    namelist/NAMCLOUDFIELD/ &
    dtav,lcloudfield, laddinfo

    dtav = dtav_glob

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMCLOUDFIELD,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMCLOUDFIELD'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMCLOUDFIELD'
      endif
      write(6 ,NAMCLOUDFIELD)
      close(ifnamopt)
    end if

    call MPI_BCAST(dtav         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(lcloudfield  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(laddinfo     ,1,MPI_LOGICAL,0,comm3d,mpierr)

    if(.not.(lcloudfield)) return
    tnext   = dtav-1e-3+btime
    dt_lim = min(dt_lim,tnext)

    if (abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if


  end subroutine initcloudfield
!>Run cloudfield. Dump the coordinates to file
   subroutine cloudfield
    use modglobal, only : imax,i1,jmax,j1,kmax, rk3step,dt_lim,timee, cexpnr,ifoutput
    use modfields, only : w0,ql0
    use modmpi,    only : cmyid
    implicit none

    integer  :: ncl
    integer,allocatable ::   ipos(:), jpos(:), kpos(:)
    integer,allocatable ::   wcl(:), qlcl(:)
    integer i,j,k,n,cldpoints

    if (.not. lcloudfield) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+dtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    ncl  = imax*jmax*kmax
    allocate (ipos(ncl),jpos(ncl),kpos(ncl),wcl(ncl),qlcl(ncl))

    if(laddinfo) then
      n = 0
      do  i=2,i1
      do  j=2,j1
      do  k=1,kmax
        if (ql0(i,j,k) > 0 ) then
          n = n + 1
          ipos(n) = i-1
          jpos(n) = j-1
          kpos(n) = k
          qlcl(n) = nint(1.0e6 * ql0(i,j,k))
          wcl(n)  = nint(100.0  * (w0(i,j,k)+w0(i,j,k+1))/2)
        end if
      end do
      end do
      end do
    else
      n = 0
      do  i=2,i1
      do  j=2,j1
      do  k=1,kmax
        if (ql0(i,j,k) > 0 ) then
          n = n + 1
          ipos(n) = i-1
          jpos(n) = j-1
          kpos(n) = k
        end if
      end do
      end do
      end do
    end if

    cldpoints = n

    open(ifoutput,file='clouds'//cmyid//'.'//cexpnr,position='append')

    if(laddinfo) then
      write(ifoutput,'(2i10)') nint(timee),cldpoints
      if( cldpoints > 0 ) then
      write(ifoutput,'(i3,2i4,2i5)') &
          (ipos(n),jpos(n),kpos(n),qlcl(n),wcl(n),n=1,cldpoints)
      end if
    else
      write(ifoutput,'(2i10)') nint(timee),cldpoints
      if( cldpoints > 0 ) then
      write(ifoutput,'(i3,2i4)') &
          (ipos(n),jpos(n),kpos(n),n=1,cldpoints)
      end if
    end if

    close(ifoutput)

  end subroutine cloudfield

end module modcloudfield
