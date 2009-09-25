!> \file modcloudfld.f90
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

module modcloudfield1


implicit none
! private
! PUBLIC :: initcloudfld
save

  real    :: dtav,tnext
  logical :: lcloudfld= .false. !< switch for writing cloud field (on/off)
  logical :: laddinfo   = .false. !< switch to write ql and w values (on/off)

contains
!> Initializing Cloudfield. Read out the namelist, initializing the variables
  subroutine initcloudfield

    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_logical
    use modglobal,only :ifnamopt,fname_options,dtmax,dtav_glob,dt_lim,btime
    implicit none
    integer :: ierr
    namelist/NAMCLOUDFIELD/ &
    dtav,lcloudfld, laddinfo

    dtav = dtav_glob

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMCLOUDFIELD,iostat=ierr)
      write(6 ,NAMCLOUDFIELD)
      close(ifnamopt)
    end if

    call MPI_BCAST(dtav         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(lcloudfld  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(laddinfo     ,1,MPI_LOGICAL,0,comm3d,mpierr)

    if(.not.(lcloudfld)) return
    tnext   = dtav-1e-3+btime
    dt_lim = min(dt_lim,tnext)

    if (abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if


  end subroutine initcloudfield
!>Run cloudfld. Dump the coordinates to file
   subroutine cloudfield
    use modglobal, only : imax,ih,i1,jmax,jtot,jh,j1,kmax,k1, rk3step,dt_lim,timee, cexpnr,ifoutput,zf
    use modfields, only : w0,ql0
    use modmpi,    only : MPI_INTEGER,MPI_Wtime,myid,comm3d,mpierr
    implicit none

    integer,allocatable,dimension (:,:,:)  :: cloudfldl,cloudfld
    integer,allocatable,dimension (:,:)  :: projcloudnr
    integer,allocatable,dimension (:,:)  :: projcloudhght,projcloudminhght
    integer,allocatable,dimension (:)   :: cloudminhght
    real :: profminhght(kmax),cfrac(kmax)
    real :: cpu_time,cpu_time0
    integer i,j,k,n,nmax,status,cc
      integer nsecs, nhrs, nminut


      nsecs   = nint(timee)
      nhrs    = int(nsecs/3600)
      nminut  = int(nsecs/60)-nhrs*60
      nsecs   = mod(nsecs,60)
    if (.not. lcloudfld) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+dtav
    dt_lim = minval((/dt_lim,tnext-timee/))

if(myid==0) CPU_time0 = MPI_Wtime()

    allocate(cloudfldl(2-ih:i1+ih,2-jh:j1+jh,1:k1))
    allocate(cloudfld(1:imax,1:jtot,1:kmax))
    allocate(projcloudhght(1:imax,1:jtot),projcloudminhght(1:imax,1:jtot),projcloudnr(1:imax,1:jtot))

    cloudfldl = 0
    where(ql0>0.) cloudfldl = -1

    call MPI_Gather( cloudfldl(2:i1,2:j1,1:kmax), imax*jmax*kmax, MPI_INTEGER, cloudfld,imax*jmax*kmax, MPI_INTEGER,0,comm3d, status, mpierr)
    if (myid==0) then
      n=0
      do i=1,imax
      do j=1,jmax
      do k=1,kmax
        if (cloudfld(i,j,k)<0) then
          n=n+1
          call findneighbour(cloudfld,i,j,k,n)
        end if
      end do
      end do
      end do
      nmax=maxval(cloudfld)
write (*,*) 'CLOUDFIELD lost clouds ', count(cloudfld<0), ' number of clouds', nmax
      allocate(cloudminhght(nmax))
      cloudminhght=kmax
      do i=1,imax
      do j=1,jmax
      do k=1,kmax
        if (cloudfld(i,j,k)>0) then
          projcloudnr(i,j) = cloudfld(i,j,k)
          projcloudhght(i,j) = k
          if (k<cloudminhght(n)) then
            cloudminhght(n) = k
          end if
          exit
        end if
      end do
      end do
      end do
      cc = count(projcloudnr>0)

      do i=1,imax
      do j=1,jmax
        if (projcloudnr(i,j)>0) then
          projcloudminhght(i,j) = cloudminhght(projcloudnr(i,j))
        end if
      end do
      end do
      do k=1,kmax
        cfrac(k) = 1.0*count(cloudfld<0)
        profminhght(k) = 1.0*count(projcloudminhght==k)
      end do
      cfrac       = cfrac/cc
      profminhght = profminhght/cc

      open (ifoutput,file='cloudoverlap.'//cexpnr,position='append')
      write(ifoutput,'(//A,/A,I4,A,I2,A,I2,A)') &
      '#--------------------------------------------------------'      &
      ,'#--- AVERAGING TIMESTEP --- '      &
      ,nhrs,':',nminut,':',nsecs      &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
      write (ifoutput,'(A/A)') &
          '#--------------------------------------------------------------------------' &
          ,'#LEV HGHT    #CFrac     #CBases'
      do k=1,kmax
        write(ifoutput,'(I3,F8.2,2E13.4)') &
            k,zf(k),&
            cfrac(k), &
            profminhght(k)
      end do
      close (ifoutput)


 CPU_time = MPI_Wtime() - CPU_time0
write(6,*)'CLOUDFIELD CPU time = ', CPU_time

    end if

  deallocate(cloudfld, cloudfldl,projcloudhght,projcloudnr)
  end subroutine cloudfield
  recursive subroutine findneighbour(cloudfld,i,j,k,n)
    use modglobal, only : imax,jtot,kmax
    integer, intent(inout) ::cloudfld(:,:,:)
    integer, intent(in) :: i,j,k,n
    integer :: ii,iii,jj,jjj,kk,kkk
    cloudfld(i,j,k) = n
    do ii=-1,1
    do jj=-1,1
    do kk=-1,1
      iii = modulo(i+ii-1,imax)+1
      jjj = modulo(j+jj-1,jtot)+1
      kkk = k+kk
      if (cloudfld(iii,jjj,kkk)<0) then
        call findneighbour(cloudfld,iii,jjj,kkk,n)
      end if
    end do
    end do
    end do
  end subroutine findneighbour
end module modcloudfield1