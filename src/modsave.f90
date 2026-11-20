!> \file modsave.f90
!! Writes restart and data files.
!>
!
!! \author Jasper Tomas, June 4th 2015
!! \par Revision list
!!   Dipanjan Majumdar, ICL (2025-2026)
!! \todo documentation
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
module modsave

implicit none
save

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writerestartfiles

    use mpi
    use modfields, only : u0,v0,w0,thl0,qt0,ql0,ql0h,e120,pres0,sv0,mindist,wall
    use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh,trestart,tnextrestart,timee,&
                          cexpnr,rk3step,ifoutput,nsv,dt,ntrun,&
                          iinletgen,timee,nstore,&
                          nfcts,nfaclyrs,lEB,iwallmom,iwalltemp,iwallmoist,ifacrestartversion,&
                          tEB,tnextEB
    use modmpi,    only : cmyidx,cmyidy,myid,slabsum,excjs,comm3d
    use modsubgriddata, only : ekm
    use modinletdata, only   : nstepread
    use initfac,   only : facT,facwsoil,facqsat,fachurel,facf,fachfi,facefi

    implicit none
    logical :: lexitnow = .false.
    integer :: i,j,k,n
    integer :: ierr, err_code
    character(len=25) :: name

    if (timee == 0) return
!    if (rk3step /=3) return
    if ((iinletgen==2) .and. (nstepread==nstore)) then                ! This overrules the need for rk3step to be 3 in case of reading inletfiles
      write(6,*) 'Writing restartfiles after reading in new inletfiles'
    else
      if (rk3step /=3) return   ! Normal check
    end if

    if (myid == 0) then
      name = 'exit_now.'//cexpnr
      inquire(file=trim(name), EXIST=lexitnow)
    end if
    call MPI_Bcast(lexitnow, 1, MPI_LOGICAL, 0, comm3d, ierr)
    if (ierr /= 0) then
      if (myid == 0) then
        print *, "Error in MPI Broadcast!"
      end if
      err_code = ierr
      call MPI_Abort(MPI_COMM_WORLD, err_code, ierr)
    end if

    if (((timee>=tnextrestart)) .or. ((lexitnow) .or. (nstepread == nstore+1))) then
      tnextrestart = tnextrestart+trestart

      name = 'initd        _   _   .'
      write (name(6:13)  ,'(i8.8)') ntrun
      name(15:17)= cmyidx
      name(19:21)= cmyidy
      name(23:25)= cexpnr
      open  (ifoutput,file=name,form='unformatted',status='replace')

      write(ifoutput)  (((mindist(i,j,k),i=ib,ie  ),j=jb,je      ),k=kb,ke   )
      write(ifoutput)  ((((wall(i,j,k,n),i=ib,ie  ),j=jb,je      ),k=kb,ke   ),n=1,5)
      write(ifoutput)  (((u0 (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((v0 (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((w0    (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((pres0 (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((thl0  (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((e120  (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((ekm   (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((qt0   (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((ql0   (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((ql0h  (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  timee,  dt

      if (myid==0) then
        write(*,*) '-------------------------'
        write(*,*) 'Saving initd restart file'
        write(*,*) 'ntrun ::: ', ntrun
        write(*,*) 'timee ::: ', timee
        write(*,*) '-------------------------'
      endif

      close (ifoutput)

      if (nsv>0) then
        name  = 'inits        _   _   .'
        write (name(6:13) ,'(i8.8)') ntrun
        name(15:17) = cmyidx
        name(19:21) = cmyidy
        name(23:25) = cexpnr
        open  (ifoutput,file=name,form='unformatted')
        write(ifoutput) ((((sv0(i,j,k,n),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh),n=1,nsv)
        write(ifoutput)  timee

        close (ifoutput)
      end if

      ! Facet state (surface energy balance). The facet arrays are global rather than
      ! domain-decomposed, so this file is written by rank 0 only and carries no rank
      ! fields in its name. Same guard as the facet allocation in initfac.
      if ((nfcts > 0) .and. (myid == 0) .and. &
          ((lEB) .or. (iwalltemp == 2) .or. (iwallmom == 2) .or. (iwallmoist == 2))) then
        name  = 'initf        .'
        write (name(6:13) ,'(i8.8)') ntrun
        name(15:17) = cexpnr
        open  (ifoutput,file=name,form='unformatted',status='replace')
        write(ifoutput)  ifacrestartversion, nfcts, nfaclyrs
        write(ifoutput)  ((facT(n,j),n=1,nfcts),j=1,nfaclyrs+1)
        if (lEB) then
          write(ifoutput)  (facwsoil(n),n=1,nfcts)   ! rank-0 array, only allocated for lEB
        end if
        ! facqsat, fachurel and facf are derived quantities (from facT, facwsoil and
        ! netsw) that the energy balance recomputes at every EB step. They are stored
        ! anyway so that the wall functions, which read them directly, see the
        ! checkpoint's exact values in the steps between the restart and the first EB
        ! update; facf in particular depends on netsw at the time of the last update,
        ! which cannot be recovered from the inputs without replaying that update.
        write(ifoutput)  (facqsat(n),n=1,nfcts)
        write(ifoutput)  (fachurel(n),n=1,nfcts)
        write(ifoutput)  ((facf(n,j),n=1,nfcts),j=1,5)
        write(ifoutput)  timee
        if (lEB) then
          ! Energy-balance interval state. The energy balance is
          ! evaluated every dtEB; in between, intqH accumulates the facet heat and
          ! moisture fluxes into fachfi/facefi (rank-0 arrays), which the next EB step
          ! divides by the elapsed interval timee - tEB. Storing the interval clock and
          ! the accumulators lets a restart that falls between two EB steps continue
          ! the interval in progress instead of restarting it at btime and dropping
          ! the flux accumulated since the last update.
          write(ifoutput)  tEB, tnextEB
          write(ifoutput)  (fachfi(n),n=1,nfcts)
          write(ifoutput)  (facefi(n),n=1,nfcts)
        end if

        write(*,*) '-------------------------'
        write(*,*) 'Saving initf restart file'
        write(*,*) 'ntrun ::: ', ntrun
        write(*,*) 'timee ::: ', timee
        write(*,*) '-------------------------'

        close (ifoutput)
      end if

      if (myid==0) then
        write(*,'(A,F15.7,A,I4)') 'dump at time = ',timee,' unit = ',ifoutput
      end if

    end if

  end subroutine writerestartfiles

end module modsave
