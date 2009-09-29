!> \file modsampling.f90
!!  Calculates statistics under conditional criteria


!>
!!  Calculates statistics under conditional criteria
!>
!!  Calculates statistics under conditional criteria, for instance over updrafts or cloudy parts of the domain. Written to $sampname_fld.expnr and to $sampname_flx.expnr.
!! If netcdf is true, this module also writes in the profiles.expnr.nc output
!! Currently implemented criteria for sampling are:
!! - Updraft (w>0)
!! - Cloud (ql>0)
!! - Buoyant updraft (w>0 and thv>0)
!! - Cloud core (ql>0 and thv>0)
!!
!!  \author Thijs Heus, MPI-M
!!  \author Pier Siebesma, KNMI
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
module modsampling



implicit none
private
PUBLIC :: initsampling, sampling, exitsampling
save
!NetCDF variables
  integer,parameter :: nvar = 13
  character(80),allocatable,dimension(:,:,:) :: ncname
  real    :: dtav, timeav,tnext,tnextwrite
  integer :: nsamples,isamp,isamptot
  character(20),dimension(10) :: samplname,longsamplname
  logical :: lsampcl  = .false. !< switch for conditional sampling cloud (on/off)
  logical :: lsampco  = .false. !< switch for conditional sampling core (on/off)
  logical :: lsampup  = .false. !< switch for conditional sampling updraft (on/off)
  logical :: lsampbuup  = .false. !< switch for conditional sampling buoyant updraft (on/off)
  real, allocatable, dimension(:,:) ::  wavl,tlavl,tvavl,qtavl,qlavl,nrsampl,massflxavl, &
                                        wtlavl,wtvavl,wqtavl,wqlavl,uwavl,vwavl

contains
!> Initialization routine, reads namelists and inits variables
  subroutine initsampling

    use modmpi,    only : comm3d, my_real,mpierr,myid,mpi_logical
    use modglobal, only : ladaptive, dtmax,rk3step,k1,ifnamopt,fname_options,   &
                           dtav_glob,timeav_glob,dt_lim,btime
    use modstat_nc, only : lnetcdf, redefine_nc,define_nc,ncinfo
    use modgenstat, only : dtav_prof=>dtav, timeav_prof=>timeav,ncid_prof=>ncid
    implicit none

    integer :: ierr

    namelist/NAMSAMPLING/ &
    dtav,timeav,lsampcl,lsampco,lsampup,lsampbuup

    dtav=dtav_glob;timeav=timeav_glob

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMSAMPLING,iostat=ierr)
      write(6 ,NAMSAMPLING)
      close(ifnamopt)
    end if

    call MPI_BCAST(timeav     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dtav       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampco,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampup,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampbuup,1,MPI_LOGICAL,0,comm3d,mpierr)

    isamptot = 0
    if (lsampup) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'upd'
      longsamplname(isamptot) = 'Updraft '
    end if
    if (lsampbuup) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'buupd'
      longsamplname(isamptot) = 'Buoyant Updraft '
    end if
    if (lsampcl) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cld'
      longsamplname(isamptot) = 'Cloud '
    end if
    if (lsampco) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cldcr'
      longsamplname(isamptot) = 'Cloud Core '
    end if

    if(isamptot == 0) return
    tnext   = dtav-1e-3+btime
    tnextwrite = timeav-1e-3+btime
    nsamples = nint(timeav/dtav)
    dt_lim = min(dt_lim,tnext)

    if (abs(timeav/dtav-nint(timeav/dtav))>1e-4) then
      stop 'timeav must be a integer multiple of dtav'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if

    allocate( wavl(k1,isamptot),tlavl(k1,isamptot),tvavl(k1,isamptot), &
          qtavl(k1,isamptot),qlavl(k1,isamptot),nrsampl(k1,isamptot),massflxavl(k1,isamptot), &
          wtlavl(k1,isamptot),wtvavl(k1,isamptot),wqtavl(k1,isamptot), &
          wqlavl(k1,isamptot),uwavl(k1,isamptot),vwavl(k1,isamptot))

!initialize variables
    nrsampl    = 0.0
    wavl       = 0.0
    tlavl      = 0.0
    tvavl      = 0.0
    qtavl      = 0.0
    qlavl      = 0.0
    massflxavl = 0.0
    wtlavl     = 0.0
    wtvavl     = 0.0
    wqtavl     = 0.0
    wqlavl     = 0.0
    uwavl      = 0.0
    vwavl      = 0.0
    if (lnetcdf) then
      dtav = dtav_prof
      timeav = timeav_prof
       tnext      = dtav-1e-3+btime
      tnextwrite = timeav-1e-3+btime
      nsamples = nint(timeav/dtav)
     if (myid==0) then
        allocate(ncname(nvar,4,isamptot))
        do isamp=1,isamptot
          call ncinfo(ncname( 1,:,isamp),'nrsamp'//samplname(isamp),trim(longsamplname(isamp))//' '//'number of points','-','tt')
          call ncinfo(ncname( 2,:,isamp),'w'//samplname(isamp),trim(longsamplname(isamp))//' '//'mean vertical velocity','m/s','mt')
          call ncinfo(ncname( 3,:,isamp),'tl'//samplname(isamp),trim(longsamplname(isamp))//' '//'mean liquid water potential temperature','K','tt')
          call ncinfo(ncname( 4,:,isamp),'qt'//samplname(isamp),trim(longsamplname(isamp))//' '//'mean total water content','kg/kg','tt')
          call ncinfo(ncname( 5,:,isamp),'ql'//samplname(isamp),trim(longsamplname(isamp))//' '//'mean liquid water content','kg/kg','tt')
          call ncinfo(ncname( 6,:,isamp),'tv'//samplname(isamp),trim(longsamplname(isamp))//' '//'mean virtual potential temperature','K','tt')
          call ncinfo(ncname( 7,:,isamp),'massflx'//samplname(isamp),trim(longsamplname(isamp))//' '//'mass flux','m^3/s','mt')
          call ncinfo(ncname( 8,:,isamp),'wtl'//samplname(isamp),trim(longsamplname(isamp))//' '//'theta_l flux','K m/s','mt')
          call ncinfo(ncname( 9,:,isamp),'wqt'//samplname(isamp),trim(longsamplname(isamp))//' '//'total water flux','kg/kg m/s','mt')
          call ncinfo(ncname(10,:,isamp),'wql'//samplname(isamp),trim(longsamplname(isamp))//' '//'liquid water flux','kg/kg m/s','mt')
          call ncinfo(ncname(11,:,isamp),'wtv'//samplname(isamp),trim(longsamplname(isamp))//' '//'theta_v flux','K m/s','mt')
          call ncinfo(ncname(12,:,isamp),'uw'//samplname(isamp),trim(longsamplname(isamp))//' '//'uw flux','m^2/s^2','mt')
          call ncinfo(ncname(13,:,isamp),'vw'//samplname(isamp),trim(longsamplname(isamp))//' '//'vw flux','m^2/s^2','mt')

          call redefine_nc(ncid_prof)
          call define_nc( ncid_prof, NVar, ncname(:,:,isamp))
        end do
     end if

   end if


  end subroutine initsampling
!> Cleans up after the run
  subroutine exitsampling
    use modstat_nc, only : lnetcdf
    use modmpi,     only : myid
    implicit none

    if (isamptot>0) then
       deallocate( wavl,tlavl,tvavl,qtavl,qlavl,nrsampl,massflxavl, &
                wtlavl,wtvavl,wqtavl,wqlavl,uwavl,vwavl)
       if (lnetcdf .and. myid==0) deallocate(ncname)
    end if
  end subroutine exitsampling
!> General routine, does the timekeeping
  subroutine sampling
    use modglobal, only : rk3step,timee,dt_lim
    implicit none
   if (isamptot==0) return
    if (rk3step/=3) return
    if(timee<tnext .and. timee<tnextwrite) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if
    if (timee>=tnext) then
      tnext = tnext+dtav
      do isamp = 1,isamptot
        call dosampling
      end do
    end if
    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+timeav
      call writesampling
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))

  return
  end subroutine sampling
!> Performs the actual sampling
  subroutine dosampling
    use modglobal, only : i1,i2,j1,j2,kmax,k1,ih,jh,&
                          dx,dy,dzh,dzf,cp,rv,rlv,rd,rslabs
    use modfields, only : u0,v0,w0,thl0,thl0h,qt0,qt0h,ql0,ql0h,thv0h,exnf,exnh
    use modsubgrid,only : ekh,ekm
    use modmpi    ,only : slabsum
    implicit none

    logical, allocatable, dimension(:,:,:) :: mask
    real, allocatable, dimension(:,:,:) :: wtlt,wqtt,wqlt,wtvt,uwt,vwt
    real, allocatable, dimension(:,:) :: w0f
    real, allocatable, dimension(:,:,:) :: thv0
    real, allocatable, dimension(:) :: thvav

    integer :: i,j,k,km,kp
    real :: cqt,cthl,den,ekhalf,c2,c1,t0h,qs0h,ekav
    real :: uws,uwr,vws,vwr,wtvs,wtvr,wtls,wtlr,wqtr,wqts,wqlr,wqls
    allocate(thvav(k1))
    allocate(mask(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate( wtlt(2-ih:i1+ih,2-jh:j1+jh,k1),&
              wtvt(2-ih:i1+ih,2-jh:j1+jh,k1),&
              wqtt(2-ih:i1+ih,2-jh:j1+jh,k1),&
              wqlt(2-ih:i1+ih,2-jh:j1+jh,k1),&
              uwt(2-ih:i1+ih,2-jh:j1+jh,k1),&
              vwt(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thv0(2-ih:i1+ih,2-jh:j1+jh,k1),&
              w0f(2-ih:i1+ih,2-jh:j1+jh))

    do k=1,k1
      thv0(2:i1,2:j1,k) = (thl0(2:i1,2:j1,k)+rlv*ql0(2:i1,2:j1,k)/(cp*exnf(k))) &
                  *(1+(rv/rd-1)*qt0(2:i1,2:j1,k)-rv/rd*ql0(2:i1,2:j1,k))
    end do

    mask = .false.
    select case (samplname(isamp))
    case ('upd')
      where(w0>0.0) mask = .true.
    case ('buup')
      thvav = 0.0
      call slabsum(thvav,1,k1,thv0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      thvav = thvav/rslabs
      do k=1,kmax
        where(w0(2:i1,2:j1,k)>0.0 .and. thv0(2:i1,2:j1,k) > thvav(k)) mask(2:i1,2:j1,k) = .true.
      end do
    case ('cld')
      where(ql0>epsilon(1.0)) mask = .true.
    case ('cldcr')
      thvav = 0.0
      call slabsum(thvav,1,k1,thv0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      thvav = thvav/rslabs
      do k=1,kmax
        where(ql0(2:i1,2:j1,k)>0.0 .and. thv0(2:i1,2:j1,k) > thvav(k)) mask(2:i1,2:j1,k) = .true.
      end do
    end select


!calculate fluxes at half levels
    wtlt = 0.0
    wtvt = 0.0
    wqtt = 0.0
    wqlt = 0.0
    uwt  = 0.0
    vwt  = 0.0
    do j=2,j1
    do i=2,i1
    do k=2,kmax
      km = k-1
      kp = k+1

!     -----------------------------------------------------------
!     calculate prefactors c1 and c2 for subgrid wthv and cthl
!     cqt for wql fluxes at half levels
!     -----------------------------------------------------------
      qs0h  =  (qt0h(i,j,k) - ql0h(i,j,k))
      t0h   =  exnh(k)*thl0h(i,j,k) + (rlv/cp)*ql0h(i,j,k)
      den   = 1. + (rlv**2)*qs0h/(rv*cp*(t0h**2))
      cthl  = (exnh(k)*cp/rlv)*((1-den)/den)
      cqt   =  (1./den)
      if (ql0h(i,j,k)>0) then
        c1    = (1.-qt0h(i,j,k)+rv/rd*qs0h*(1.+rlv/(rv*t0h)))/den
        c2    = c1*rlv/(t0h*cp)-1.
      else
        c1    = 1. + (rv/rd-1)*qt0h(i,j,k)
        c2    = rv/rd-1
      end if

!     -----------------------------------------------------------
!     calculate resolved and subgrid fluxes at half levels
!     -----------------------------------------------------------

      ekhalf = (ekh(i,j,k) *dzf(km)+ekh(i,j,km)*dzf(k) )/(2*dzh(k) )

      wtls   = -ekhalf*(thl0(i,j,k)-thl0(i,j,km))/dzh(k)
      wtlr   = w0(i,j,k)*thl0h(i,j,k)
      wtlt(i,j,k)   = wtls + wtlr

      wqts   = -ekhalf*(qt0(i,j,k)-qt0(i,j,km))/dzh(k)
      wqtr   = w0(i,j,k)*qt0h(i,j,k)
      wqtt(i,j,k)   = wqts + wqtr

      if (ql0h(i,j,k)>0) then
        wqls   = cthl*wtls+ cqt*wqts
      else
        wqls   = 0.
      end if

      wqlr   = w0(i,j,k)*ql0h(i,j,k)
      wqlt(i,j,k)   = wqlr + wqls

      wtvs   = c1*wtls + c2*thl0h(i,j,k)*wqts
      wtvr   = w0(i,j,k)*thv0h(i,j,k)
      wtvt(i,j,k)   = wtvr + wtvs

      ekav = 0.25 * &
              ( ekm(i,j,k-1)+ekm(i,j,k)+ekm(i-1,j,k-1)+ekm(i-1,j,k) )
      uws  =  - ekav * &
            ((u0(i,j,k)-u0(i,j,k-1))/dzf(k)+(w0(i,j,k)-w0(i-1,j,k))/dx)
      uwr = (w0(i,j,k)+w0(i-1,j,k))*(u0(i,j,k-1)+u0(i,j,k))/4.
      uwt(i,j,k) = uws + uwr

      ekav = 0.25 * &
          ( ekm(i,j,k-1)+ekm(i,j,k)+ekm(i,j-1,k-1)+ekm(i,j-1,k) )
      vws  = - ekav * &
          ((v0(i,j,k)-v0(i,j,k-1))/dzf(k)+(w0(i,j,k)-w0(i,j-1,k))/dy)
      vwr = (w0(i,j,k)+w0(i,j-1,k))*(v0(i,j,k-1)+v0(i,j,k))/4.
      vwt(i,j,k) = vws + vwr
    end do
    end do
    end do

!add fields and fluxes to mean
!     1)       fields on full levels
    do k=1,kmax
      w0f(2:i1,2:j1) = 0.5*(w0  (2:i1,2:j1,k) + w0  (2:i1,2:j1,k+1))
      nrsampl(k,isamp)= nrsampl(k,isamp)+count(mask(2:i1,2:j1,k))
      wavl (k,isamp) = wavl (k,isamp)+sum(w0f (2:i1,2:j1  ),mask(2:i1,2:j1,k))
      tlavl(k,isamp) = tlavl(k,isamp)+sum(thl0(2:i1,2:j1,k),mask(2:i1,2:j1,k))
      tvavl(k,isamp) = tvavl(k,isamp)+sum(thv0(2:i1,2:j1,k),mask(2:i1,2:j1,k))
      qtavl(k,isamp) = qtavl(k,isamp)+sum(qt0(2:i1,2:j1,k),mask(2:i1,2:j1,k))
      qlavl(k,isamp) = qlavl(k,isamp)+sum(ql0(2:i1,2:j1,k),mask(2:i1,2:j1,k))
    end do

!     2)       fluxes on half levels
    massflxavl(1,isamp) =0.
    wtlavl(1,isamp) = 0.
    wtvavl(1,isamp) = 0.
    wqtavl(1,isamp) = 0.
    wqlavl(1,isamp) = 0.
    uwavl (1,isamp) = 0.
    vwavl (1,isamp) = 0.
    do k=2,kmax
      massflxavl(k,isamp) = massflxavl(k,isamp)+0.5*(sum(w0 (2:i1,2:j1,k),mask(2:i1,2:j1,k)) &
                                                    +sum(w0 (2:i1,2:j1,k),mask(2:i1,2:j1,k-1)))
      wtlavl(k,isamp) = wtlavl(k,isamp)+0.5*(sum(wtlt(2:i1,2:j1,k),mask(2:i1,2:j1,k)) &
                                            +sum(wtlt(2:i1,2:j1,k),mask(2:i1,2:j1,k-1)))
      wtvavl(k,isamp) = wtvavl(k,isamp)+0.5*(sum(wtvt(2:i1,2:j1,k),mask(2:i1,2:j1,k)) &
                                            +sum(wtvt(2:i1,2:j1,k),mask(2:i1,2:j1,k-1)))
      wqtavl(k,isamp) = wqtavl(k,isamp)+0.5*(sum(wqtt(2:i1,2:j1,k),mask(2:i1,2:j1,k)) &
                                            +sum(wqtt(2:i1,2:j1,k),mask(2:i1,2:j1,k-1)))
      wqlavl(k,isamp) = wqlavl(k,isamp)+0.5*(sum(wqlt(2:i1,2:j1,k),mask(2:i1,2:j1,k)) &
                                            +sum(wqlt(2:i1,2:j1,k),mask(2:i1,2:j1,k-1)))
      uwavl (k,isamp) = uwavl (k,isamp)+0.5*(sum(uwt (2:i1,2:j1,k),mask(2:i1,2:j1,k)) &
                                            +sum(uwt (2:i1,2:j1,k),mask(2:i1,2:j1,k-1)))
      vwavl (k,isamp) = vwavl (k,isamp)+0.5*(sum(vwt (2:i1,2:j1,k),mask(2:i1,2:j1,k)) &
                                            +sum(vwt (2:i1,2:j1,k),mask(2:i1,2:j1,k-1)))
    end do

    deallocate(mask,wtlt,wqtt,wqlt,wtvt,uwt,vwt,thvav,w0f,thv0)


  end subroutine dosampling
!> Write the statistics to file
  subroutine writesampling

    use modglobal, only : timee,k1,kmax,zf,zh,cexpnr,ifoutput,rslabs
    use modfields, only : presf,presh
    use modmpi,    only : myid,my_real,comm3d,mpierr,mpi_sum
    use modstat_nc, only: lnetcdf, writestat_nc
    use modgenstat, only: ncid_prof=>ncid,nrec_prof=>nrec
    implicit none
    real,dimension(k1,nvar) :: vars

    real, allocatable, dimension(:)  :: wmn,tlmn,tvmn,qtmn,qlmn,nrsampmn,massflxmn, &
                                        wtlmn,wtvmn,wqtmn,wqlmn,uwmn,vwmn
    real, allocatable, dimension(:,:):: wav,tlav,tvav,qtav,qlav,nrsamp,massflxav, &
                                        wtlav,wtvav,wqtav,wqlav,uwav,vwav

    integer :: nsecs, nhrs, nminut, k
    integer :: inorm
    allocate( wmn(k1),tlmn(k1),tvmn(k1),qtmn(k1),qlmn(k1),nrsampmn(k1),massflxmn(k1), &
              wtlmn(k1),wtvmn(k1),wqtmn(k1),wqlmn(k1),uwmn(k1),vwmn(k1))
    allocate( wav(k1,isamptot),tlav(k1,isamptot),tvav(k1,isamptot), &
              qtav(k1,isamptot),qlav(k1,isamptot),nrsamp(k1,isamptot),massflxav(k1,isamptot), &
              wtlav(k1,isamptot),wtvav(k1,isamptot),wqtav(k1,isamptot), &
              wqlav(k1,isamptot),uwav(k1,isamptot),vwav(k1,isamptot))
    nsecs   = nint(timee)
    nhrs    = int(nsecs/3600)
    nminut  = int(nsecs/60)-nhrs*60
    nsecs   = mod(nsecs,60)
    inorm   = nint(rslabs*timeav/dtav)

    call MPI_ALLREDUCE(nrsampl   ,nrsamp   ,isamptot*kmax,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wavl      ,wav      ,isamptot*kmax,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(tlavl     ,tlav     ,isamptot*kmax,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(tvavl     ,tvav     ,isamptot*kmax,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(qtavl     ,qtav     ,isamptot*kmax,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(qlavl     ,qlav     ,isamptot*kmax,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(massflxavl,massflxav,isamptot*kmax,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wtlavl    ,wtlav    ,isamptot*kmax,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wtvavl    ,wtvav    ,isamptot*kmax,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wqtavl    ,wqtav    ,isamptot*kmax,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wqlavl    ,wqlav    ,isamptot*kmax,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(uwavl     ,uwav     ,isamptot*kmax,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(vwavl     ,vwav     ,isamptot*kmax,MY_REAL,MPI_SUM,comm3d,mpierr)
!reset variables
    nrsampl    = 0.0
    wavl       = 0.0
    tlavl      = 0.0
    tvavl      = 0.0
    qtavl      = 0.0
    qlavl      = 0.0
    massflxavl = 0.0
    wtlavl     = 0.0
    wtvavl     = 0.0
    wqtavl     = 0.0
    wqlavl     = 0.0
    uwavl      = 0.0
    vwavl      = 0.0

    if (myid==0) then
      do isamp = 1,isamptot
!normalize variables
        where (nint(nrsamp(1:kmax,isamp))>0)
          wmn  = wav (1:kmax,isamp)/nrsamp(1:kmax,isamp)
          tlmn = tlav(1:kmax,isamp)/nrsamp(1:kmax,isamp)
          tvmn = tvav(1:kmax,isamp)/nrsamp(1:kmax,isamp)
          qtmn = qtav(1:kmax,isamp)/nrsamp(1:kmax,isamp)
          qlmn = qlav(1:kmax,isamp)/nrsamp(1:kmax,isamp)
        elsewhere
          wmn  = 0.
          tlmn = 0.
          tvmn = 0.
          qtmn = 0.
          qlmn = 0.
        end where
        nrsampmn  = nrsamp   (1:kmax,isamp)/inorm
        massflxmn = massflxav(1:kmax,isamp)/inorm
        wtlmn = wtlav(1:kmax,isamp)/inorm
        wtvmn = wtvav(1:kmax,isamp)/inorm
        wqtmn = wqtav(1:kmax,isamp)/inorm
        wqlmn = wqlav(1:kmax,isamp)/inorm
        uwmn  = uwav (1:kmax,isamp)/inorm
        vwmn  = vwav (1:kmax,isamp)/inorm


!write files
        open (ifoutput,file=trim(samplname(isamp))//'fld.'//cexpnr,position='append')
        write(ifoutput,'(//3A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
          '#------------------------------ ',trim(longsamplname(isamp)),' -------------------------------'      &
          ,'#',timeav,'--- AVERAGING TIMESTEP --- '      &
          ,nhrs,':',nminut,':',nsecs      &
          ,'   HRS:MIN:SEC AFTER INITIALIZATION '

        write (ifoutput,'(2A/2A)') &
           '#------------------------------------------------------' &
           ,'------------------------------' &
           ,'  LEV HGHT   PRES       COV       W       THL      QT      ' &
           ,'QL       THV   '
        do k=1,kmax
          write(ifoutput,'(i5,F6.0,F7.1,F9.4,5F9.3)') &
              k, &
              zf      (k), &
              presf   (k)/100., &
              nrsampmn(k), &
              wmn     (k), &
              tlmn    (k), &
              qtmn    (k)*1000., &
              qlmn    (k)*1000., &
              tvmn    (k)
        end do
        close(ifoutput)

        open (ifoutput,file=trim(samplname(isamp))//'flx.'//cexpnr,position='append')
        write(ifoutput,'(//3A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
          '#------------------------------- ',trim(longsamplname(isamp)),' ----------------------------------------'      &
          ,'#',timeav,'--- AVERAGING TIMESTEP --- '      &
          ,nhrs,':',nminut,':',nsecs      &
          ,'   HRS:MIN:SEC AFTER INITIALIZATION '

        write (ifoutput,'(2A/2A)') &
            '#------------------------------------------------------' &
            ,'------------------------------' &
            ,'   LEV  HGHT  PRES       AW                WTHL                WQT                ' &
            ,'WQL                WTHV                UW                VW'
        do k=1,kmax
          write(ifoutput,'(i5,F6.0,F7.1,7E16.8)') &
                k, &
                zh       (k), &
                presh    (k)/100., &
                massflxmn(k), &
                wtlmn    (k), &
                wqtmn    (k), &
                wqlmn    (k), &
                wtvmn    (k), &
                uwmn     (k), &
                vwmn     (k)
        end do
        close(ifoutput)
      if (lnetcdf) then
        vars(:, 1) =nrsampmn
        vars(:, 2) =wmn
        vars(:, 3) =tlmn
        vars(:, 4) =qtmn
        vars(:, 5) =qlmn
        vars(:, 6) =tvmn
        vars(:, 7) =massflxmn
        vars(:, 8) =wtlmn
        vars(:, 9) =wqtmn
        vars(:,10) =wqlmn
        vars(:,11) =wtvmn
        vars(:,12) =uwmn
        vars(:,13) =vwmn

        call writestat_nc(ncid_prof,nvar,ncname(:,:,isamp),vars(1:kmax,:),nrec_prof,kmax)
      end if

      end do
    end if
    deallocate( wmn,tlmn,tvmn,qtmn,qlmn,nrsampmn,massflxmn, &
                wtlmn,wtvmn,wqtmn,wqlmn,uwmn,vwmn)
    deallocate( wav,tlav,tvav,qtav,qlav,nrsamp,massflxav, &
                wtlav,wtvav,wqtav,wqlav,uwav,vwav)

  end subroutine writesampling

end module modsampling
