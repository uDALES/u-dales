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


  use modglobal, only : longint

implicit none
private
PUBLIC :: initsampling, sampling, exitsampling
save
!NetCDF variables
  integer,parameter :: nvar = 13
  character(80),allocatable,dimension(:,:,:) :: ncname
  real    :: dtav, timeav
  integer(kind=longint) :: idtav,itimeav,tnext,tnextwrite
  integer :: nsamples,isamp,isamptot
  character(20),dimension(10) :: samplname,longsamplname
  logical :: lsampcl  = .false. !< switch for conditional sampling cloud (on/off)
  logical :: lsampco  = .false. !< switch for conditional sampling core (on/off)
  logical :: lsampup  = .false. !< switch for conditional sampling updraft (on/off)
  logical :: lsampbuup  = .false. !< switch for conditional sampling buoyant updraft (on/off)
  real, allocatable, dimension(:,:) ::  wavl,tlavl,tvavl,qtavl,qlavl,nrsampl,massflxavl, &
                                        wtlavl,wtvavl,wqtavl,wqlavl,uwavl,vwavl
  real, allocatable, dimension(:,:) :: wwavl,pavl,dwdthavl,dwwdzhavl,dpdzhavl, &
                                       duwdxhavl,dtaudxhavl,dtaudzhavl,thvhavl, &
                                       fcoravl,nrsamphl
  real,allocatable, dimension(:,:) :: w_el,sig_el


contains
!> Initialization routine, reads namelists and inits variables
  subroutine initsampling

    use modmpi,    only : comm3d, my_real,mpierr,myid,mpi_logical
    use modglobal, only : ladaptive, dtmax,rk3step,k1,ifnamopt,fname_options,   &
                           dtav_glob,timeav_glob,dt_lim,btime,tres
    use modstat_nc, only : lnetcdf, redefine_nc,define_nc,ncinfo
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav,ncid_prof=>ncid
    implicit none

    integer :: ierr

    namelist/NAMSAMPLING/ &
    dtav,timeav,lsampcl,lsampco,lsampup,lsampbuup

    dtav=dtav_glob;timeav=timeav_glob

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMSAMPLING,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMSAMPLING'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMSAMPLING'
      endif
      write(6 ,NAMSAMPLING)
      close(ifnamopt)
    end if

    call MPI_BCAST(timeav     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dtav       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampco,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampup,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampbuup,1,MPI_LOGICAL,0,comm3d,mpierr)

    isamptot = 1
    samplname (isamptot) = 'all'
    longsamplname(isamptot) = 'All '
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
     idtav = dtav/tres
    itimeav = timeav/tres

    tnext      = idtav   +btime
    tnextwrite = itimeav +btime
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

    allocate (  nrsamphl(k1,isamptot),wwavl(k1,isamptot), &
                pavl(k1,isamptot),dwdthavl(k1,isamptot),dwwdzhavl(k1,isamptot),&
                dpdzhavl(k1,isamptot),duwdxhavl(k1,isamptot),dtaudxhavl(k1,isamptot),&
                dtaudzhavl(k1,isamptot),thvhavl(k1,isamptot),fcoravl(k1,isamptot))

    allocate (w_el(k1,isamptot),sig_el(k1,isamptot))


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

    nrsamphl   = 0.0
    wwavl      = 0.0
    pavl       = 0.0
    dwdthavl   = 0.0
    dwwdzhavl  = 0.0
    dpdzhavl   = 0.0
    duwdxhavl  = 0.0
    dtaudxhavl = 0.0
    dtaudzhavl = 0.0
    thvhavl    = 0.0
    fcoravl    = 0.0

    w_el       = 0.0
    sig_el     = 0.0

    if (lnetcdf) then
      idtav = idtav_prof
      itimeav = itimeav_prof
      tnext      = idtav+btime
      tnextwrite = itimeav+btime
      nsamples = itimeav/idtav
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
       deallocate( nrsamphl,wwavl, &
                pavl,dwdthavl,dwwdzhavl,dpdzhavl,duwdxhavl,dtaudxhavl,dtaudzhavl,thvhavl,fcoravl, &
                w_el,sig_el)
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
      tnext = tnext+idtav
      do isamp = 1,isamptot
        call dosampling
      end do
    end if
    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+itimeav
      call writesampling
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))

  return
  end subroutine sampling
!> Performs the actual sampling
  subroutine dosampling
    use modglobal, only : i1,i2,j1,j2,kmax,k1,ih,jh,&
                          dx,dy,dzh,dzf,cp,rv,rlv,rd,rslabs, &
                          grav,om22,cu,timee
    use modfields, only : u0,v0,w0,thl0,thl0h,qt0,qt0h,ql0,ql0h,thv0h,exnf,exnh, &
                          wp
    use modsubgriddata,only : ekh,ekm
    use modmpi    ,only : slabsum
    use modpois,   only : p
    use modsurfdata,only: thvs
    implicit none

    logical, allocatable, dimension(:,:,:) :: mask
    real, allocatable, dimension(:,:,:) :: wtlt,wqtt,wqlt,wtvt,uwt,vwt
    real, allocatable, dimension(:,:,:) :: w0f
    real, allocatable, dimension(:,:,:) :: thv0
    real, allocatable, dimension(:) :: thvav

    logical, allocatable, dimension(:,:,:) :: maskh
    real, allocatable, dimension(:,:,:) :: uws,vws,uwr,vwr,wwr,wwsf
    real, allocatable, dimension(:) :: thvhav


    integer :: i,j,k,km,kp
    real :: cqt,cthl,den,ekhalf,c2,c1,t0h,qs0h,ekav
    real :: wtvs,wtvr,wtls,wtlr,wqtr,wqts,wqlr,wqls
    real :: beta

    allocate(thvav(k1))
    allocate(mask(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate( wtlt(2-ih:i1+ih,2-jh:j1+jh,k1),&
              wtvt(2-ih:i1+ih,2-jh:j1+jh,k1),&
              wqtt(2-ih:i1+ih,2-jh:j1+jh,k1),&
              wqlt(2-ih:i1+ih,2-jh:j1+jh,k1),&
              uwt(2-ih:i1+ih,2-jh:j1+jh,k1),&
              vwt(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thv0(2-ih:i1+ih,2-jh:j1+jh,k1),&
              w0f(2-ih:i1+ih,2-jh:j1+jh,k1))      !Change to 3D array
    allocate(maskh(2-ih:i1+ih,2-jh:j1+jh,k1))     !New
    allocate(uws  (2-ih:i1+ih,2-jh:j1+jh,k1),&
             vws  (2-ih:i1+ih,2-jh:j1+jh,k1),&
             uwr  (2-ih:i1+ih,2-jh:j1+jh,k1),&
             vwr  (2-ih:i1+ih,2-jh:j1+jh,k1),&
             wwr  (2-ih:i1+ih,2-jh:j1+jh,k1),&
             wwsf (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thvhav(k1))
    
    beta = grav/thvs

    do k=1,k1
      thv0(2:i1,2:j1,k) = (thl0(2:i1,2:j1,k)+rlv*ql0(2:i1,2:j1,k)/(cp*exnf(k))) &
                  *(1+(rv/rd-1)*qt0(2:i1,2:j1,k)-rv/rd*ql0(2:i1,2:j1,k))
    enddo
    do k=1,kmax
      w0f (2:i1,2:j1,k) = 0.5*(w0 (2:i1,2:j1,k) + w0  (2:i1,2:j1,k+1))
    end do

    mask = .false.
    maskh= .false.

    thvav = 0.0
    call slabsum(thvav,1,k1,thv0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    thvav = thvav/rslabs

    thvhav = 0.0
    call slabsum(thvhav,1,k1,thv0h,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    thvhav = thvhav/rslabs


    select case (samplname(isamp))
    case ('upd')
      do i=2,i1
      do j=2,j1
      do k=1,kmax
         if (w0f(i,j,k).gt.0.) then
             mask(i,j,k) = .true.
         endif
         if (w0(i,j,k).gt.0.) then
             maskh(i,j,k) = .true.
         endif
      enddo
      enddo
      enddo
    case ('buup')
      do i=2,i1
      do j=2,j1
      do k=1,kmax
         if (w0f(i,j,k)>0.0.and.thv0 (i,j,k) > thvav(k)) then
             mask(i,j,k) = .true.
         endif
         if (w0(i,j,k).gt.0.0.and.thv0h(i,j,k) > thvhav(k)) then
             maskh(i,j,k) = .true.
         endif
      enddo
      enddo
      enddo

    case ('cld')
      do i=2,i1
      do j=2,j1
      do k=1,kmax
         if (ql0(i,j,k)>epsilon(1.0)) then
             mask(i,j,k) = .true.
         endif
         if (ql0h(i,j,k)>epsilon(1.0)) then
             maskh(i,j,k) = .true.
         endif
      enddo
      enddo
      enddo
    case ('cldcr')
      do i=2,i1
      do j=2,j1
      do k=1,kmax
         if (ql0(i,j,k)>epsilon(1.0).and.thv0(i,j,k) > thvav(k)) then
             mask(i,j,k) = .true.
         endif
         if (ql0h(i,j,k)>epsilon(1.0).and.thv0h(i,j,k) > thvhav(k)) then
             maskh(i,j,k) = .true.
         endif
      enddo
      enddo
      enddo

    case ('all')
        mask  = .true.
        maskh = .true.
    end select


!calculate fluxes at half levels
    wtlt = 0.0
    wtvt = 0.0
    wqtt = 0.0
    wqlt = 0.0
    uwt  = 0.0
    vwt  = 0.0

    uws  = 0.0
    vws  = 0.0
    uwr  = 0.0
    vwr  = 0.0
    wwr  = 0.0
    wwsf = 0.0

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
      uws(i,j,k)  =  - ekav * &
            ((u0(i,j,k)-u0(i,j,k-1))/dzf(k)+(w0(i,j,k)-w0(i-1,j,k))/dx)
      uwr(i,j,k) = (w0(i,j,k)+w0(i-1,j,k))*(u0(i,j,k-1)+u0(i,j,k))/4.
      uwt(i,j,k) = uws(i,j,k) + uwr(i,j,k)

      ekav = 0.25 * &
          ( ekm(i,j,k-1)+ekm(i,j,k)+ekm(i,j-1,k-1)+ekm(i,j-1,k) )
      vws(i,j,k)  = - ekav * &
          ((v0(i,j,k)-v0(i,j,k-1))/dzf(k)+(w0(i,j,k)-w0(i,j-1,k))/dy)
      vwr(i,j,k) = (w0(i,j,k)+w0(i,j-1,k))*(v0(i,j,k-1)+v0(i,j,k))/4.
      vwt(i,j,k) = vws(i,j,k) + vwr(i,j,k)

      wwr  (i,j,k) = w0(i,j,k)**2
      wwsf (i,j,k-1) = -2 * ekm(i,j,k-1) * &
               (w0(i,j,k)-w0(i,j,k-1)) / dzf(k-1)
    end do
    end do
    end do

    do k=2,kmax
       i=i2
       do j=2,j1
          ekav = 0.25 * &
               ( ekm(i,j,k-1)+ekm(i,j,k)+ekm(i-1,j,k-1)+ekm(i-1,j,k) )
          uws(i,j,k)  =  - ekav * &
            ((u0(i,j,k)-u0(i,j,k-1))/dzf(k)+(w0(i,j,k)-w0(i-1,j,k))/dx)
          uwr (i,j,k) = (w0(i,j,k)+w0(i-1,j,k))*(u0(i,j,k-1)+u0(i,j,k))/4.
       end do
       j=j2
       do i=2,i1
          ekav = 0.25 * &
             ( ekm(i,j,k-1)+ekm(i,j,k)+ekm(i,j-1,k-1)+ekm(i,j-1,k) )
          vws (i,j,k) = - ekav * &
             ((v0(i,j,k)-v0(i,j,k-1))/dzf(k)+(w0(i,j,k)-w0(i,j-1,k))/dy)
          vwr (i,j,k) = (w0(i,j,k)+w0(i,j-1,k))*(v0(i,j,k-1)+v0(i,j,k))/4.
       end do
    end do


!add fields and fluxes to mean
!     1)       fields on full levels
    do k=1,kmax
      nrsampl(k,isamp)= nrsampl(k,isamp)+count(mask(2:i1,2:j1,k))
      nrsamphl(k,isamp) = nrsamphl(k,isamp)+count(maskh(2:i1,2:j1,k))
      wavl (k,isamp) = wavl (k,isamp)+sum(w0f (2:i1,2:j1,k),mask(2:i1,2:j1,k))
      tlavl(k,isamp) = tlavl(k,isamp)+sum(thl0(2:i1,2:j1,k),mask(2:i1,2:j1,k))
      tvavl(k,isamp) = tvavl(k,isamp)+sum(thv0(2:i1,2:j1,k),mask(2:i1,2:j1,k))
      qtavl(k,isamp) = qtavl(k,isamp)+sum(qt0(2:i1,2:j1,k),mask(2:i1,2:j1,k))
      qlavl(k,isamp) = qlavl(k,isamp)+sum(ql0(2:i1,2:j1,k),mask(2:i1,2:j1,k))
      pavl (k,isamp) = pavl (k,isamp)+sum(p  (2:i1,2:j1,k),mask(2:i1,2:j1,k))
    end do

!     2)       fluxes on half levels
!    massflxavl (1,isamp) =0.  !resetting not necessary
!    wtlavl     (1,isamp) = 0.
!    wtvavl     (1,isamp) = 0.
!    wqtavl     (1,isamp) = 0.
!    wqlavl     (1,isamp) = 0.
!    uwavl      (1,isamp) = 0.
!    vwavl      (1,isamp) = 0.
!    wwavl      (1,isamp) = 0.
!    thvhavl    (1,isamp) = 0.
!    dpdzhavl   (1,isamp) = 0.
!    dwdthavl   (1,isamp) = 0.
!    dwwdzhavl  (1,isamp) = 0.
!    duwdxhavl  (1,isamp) = 0.
!    dtaudxhavl (1,isamp) = 0.
!    dtaudzhavl (1,isamp) = 0.
!    fcoravl    (1,isamp) = 0.

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
      wwavl (k,isamp) = wwavl (k,isamp)+     sum(wwr (2:i1,2:j1,k),maskh(2:i1,2:j1,k))
      w_el  (k,isamp) = sum(w0 (2:i1,2:j1,k),maskh(2:i1,2:j1,k))
      sig_el(k,isamp) = count(maskh(2:i1,2:j1,k))

      dwdthavl(k,isamp) = dwdthavl(k,isamp)+sum(wp(2:i1,2:j1,k),maskh(2:i1,2:j1,k))

      do j=2,j1
      do i=2,i1
       if (maskh(i,j,k)) then
         thvhavl   (k,isamp) = thvhavl   (k,isamp) + thv0h(i,j,k) - thvhav(k)
         dpdzhavl  (k,isamp) = dpdzhavl  (k,isamp) - (p(i,j,k)-p(i,j,k-1))/dzh(k) &
                                                   + beta * thvhav(k)
         dwwdzhavl (k,isamp) = dwwdzhavl (k,isamp) - (w0f(i,j,k)**2- w0f(i,j,k-1)**2)/dzh(k)
         duwdxhavl (k,isamp) = duwdxhavl (k,isamp) - (uwr(i+1,j,k) - uwr (i,j,k))/dx - &
                                                     (vwr(i,j+1,k) - vwr (i,j,k))/dy
         dtaudxhavl(k,isamp) = dtaudxhavl(k,isamp) - (uws(i+1,j,k) - uws (i,j,k))/dx - &
                                                     (vws(i,j+1,k) - vws (i,j,k))/dy   
         dtaudzhavl(k,isamp) = dtaudzhavl(k,isamp) - (wwsf(i,j,k)  - wwsf(i,j,k-1))/dzh(k)
         fcoravl(k,isamp) = fcoravl(k,isamp) + om22 * cu  &
                                    +( (dzf(k-1) * (u0(i,j,k)   + u0(i+1,j,k) )    &
                                    +    dzf(k)  * (u0(i,j,k-1) + u0(i+1,j,k-1))  ) / dzh(k) ) &
                                    * om22*0.25
       endif
      end do
      end do
    enddo



    deallocate(mask,wtlt,wqtt,wqlt,wtvt,uwt,vwt,thvav,w0f,thv0)
    deallocate(maskh,uws,vws,uwr,vwr)
    deallocate(wwr,wwsf,thvhav)



  end subroutine dosampling
!> Write the statistics to file
  subroutine writesampling

    use modglobal, only : rtimee,k1,kmax,zf,zh,cexpnr,ifoutput,rslabs,grav
    use modfields, only : presf,presh
    use modmpi,    only : myid,my_real,comm3d,mpierr,mpi_sum
    use modstat_nc, only: lnetcdf, writestat_nc,nc_fillvalue
    use modgenstat, only: ncid_prof=>ncid,nrec_prof=>nrec
    use modsurfdata, only: thvs

    implicit none
    real,dimension(k1,nvar) :: vars

    real, allocatable, dimension(:)  :: wmn,tlmn,tvmn,qtmn,qlmn,nrsampmn,massflxmn, &
                                        wtlmn,wtvmn,wqtmn,wqlmn,uwmn,vwmn
    real, allocatable, dimension(:,:):: wav,tlav,tvav,qtav,qlav,nrsamp,massflxav, &
                                        wtlav,wtvav,wqtav,wqlav,uwav,vwav
    real, allocatable, dimension(:)  :: nrsamphmn, wwmn,&
                                        pmn,dwdthmn,dwwdzhmn,dpdzhmn,duwdxhmn,&
                                        dtaudxhmn,dtaudzhmn,thvhmn, &
                                        fcormn
    real, allocatable, dimension(:,:):: nrsamph,wwav, & 
                                        pav,dwdthav,dwwdzhav,dpdzhav,duwdxhav,&
                                        dtaudxhav,dtaudzhav,thvhav,&
                                        fcorav,w_e,sig_e

    integer :: nsecs, nhrs, nminut, k
    integer :: inorm
    allocate( wmn(k1),tlmn(k1),tvmn(k1),qtmn(k1),qlmn(k1),nrsampmn(k1),massflxmn(k1), &
              wtlmn(k1),wtvmn(k1),wqtmn(k1),wqlmn(k1),uwmn(k1),vwmn(k1))
    allocate( wav(k1,isamptot),tlav(k1,isamptot),tvav(k1,isamptot), &
              qtav(k1,isamptot),qlav(k1,isamptot),nrsamp(k1,isamptot),massflxav(k1,isamptot), &
              wtlav(k1,isamptot),wtvav(k1,isamptot),wqtav(k1,isamptot), &
              wqlav(k1,isamptot),uwav(k1,isamptot),vwav(k1,isamptot))
    allocate (nrsamphmn(k1),wwmn(k1))
    allocate( pmn(k1),dwdthmn(k1),dwwdzhmn(k1),dpdzhmn(k1),duwdxhmn(k1),&
              dtaudxhmn(k1),dtaudzhmn(k1),thvhmn(k1),fcormn(k1))
    allocate( nrsamph(k1,isamptot),wwav(k1,isamptot))
    allocate( pav(k1,isamptot),dwdthav(k1,isamptot),dwwdzhav(k1,isamptot),&
              dpdzhav(k1,isamptot),duwdxhav(k1,isamptot),&
              dtaudxhav(k1,isamptot),dtaudzhav(k1,isamptot),thvhav(k1,isamptot),fcorav(k1,isamptot))
    allocate (w_e(k1,isamptot),sig_e(k1,isamptot))

    nsecs   = nint(rtimee)
    nhrs    = int(nsecs/3600)
    nminut  = int(nsecs/60)-nhrs*60
    nsecs   = mod(nsecs,60)
    inorm   = nint(rslabs*timeav/dtav)

    call MPI_ALLREDUCE(nrsampl   ,nrsamp   ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wavl      ,wav      ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(tlavl     ,tlav     ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(tvavl     ,tvav     ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(qtavl     ,qtav     ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(qlavl     ,qlav     ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(massflxavl,massflxav,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wtlavl    ,wtlav    ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wtvavl    ,wtvav    ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wqtavl    ,wqtav    ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wqlavl    ,wqlav    ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(uwavl     ,uwav     ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(vwavl     ,vwav     ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)

    call MPI_ALLREDUCE(nrsamphl  ,nrsamph  ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wwavl     ,wwav     ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(pavl      ,pav      ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(dwdthavl  ,dwdthav  ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(dwwdzhavl ,dwwdzhav ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(dpdzhavl  ,dpdzhav  ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(duwdxhavl ,duwdxhav ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(dtaudxhavl,dtaudxhav,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(dtaudzhavl,dtaudzhav,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(thvhavl   ,thvhav   ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(fcoravl   ,fcorav   ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)

    call MPI_ALLREDUCE(w_el      ,w_e      ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(sig_el    ,sig_e    ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)



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

    nrsamphl   = 0.0
    wwavl      = 0.0
    pavl       = 0.0
    dwdthavl   = 0.0
    dwwdzhavl  = 0.0
    dpdzhavl   = 0.0
    duwdxhavl  = 0.0
    dtaudxhavl = 0.0
    dtaudzhavl = 0.0
    thvhavl    = 0.0
    fcoravl    = 0.0

    w_el = 0.
    sig_el = 0.


    if (myid==0) then
      do isamp = 1,isamptot

         wmn  = 0.
         tlmn = 0.
         tvmn = 0.
         qtmn = 0.
         qlmn = 0.

         pmn  = 0.
         dwdthmn  = 0.
         dwwdzhmn = 0.
         dpdzhmn  = 0.
         duwdxhmn = 0.
         dtaudxhmn= 0.
         dtaudzhmn= 0.
         thvhmn   = 0.
         fcormn   = 0.
         wwmn     = 0.

!normalize variables

        do k=1,kmax
          if (nrsamp(k,isamp)>0) then
            wmn  (k) = wav (k,isamp)/nrsamp(k,isamp)
            tlmn (k) = tlav(k,isamp)/nrsamp(k,isamp)
            tvmn (k) = tvav(k,isamp)/nrsamp(k,isamp)
            qtmn (k) = qtav(k,isamp)/nrsamp(k,isamp)
            qlmn (k) = qlav(k,isamp)/nrsamp(k,isamp)
          endif

          if (sig_e(k,isamp).gt.0) then
            w_e  (k,isamp) = w_e  (k,isamp) / sig_e(k,isamp)
            sig_e(k,isamp) = sig_e(k,isamp) / rslabs
          endif

          if (nrsamph(k,isamp).gt.0) then
            pmn      (k) = pav      (k,isamp)/nrsamph(k,isamp)
            dwdthmn  (k) = dwdthav  (k,isamp)/nrsamph(k,isamp)
            dwwdzhmn (k) = dwwdzhav (k,isamp)/nrsamph(k,isamp)
            dpdzhmn  (k) = dpdzhav  (k,isamp)/nrsamph(k,isamp)
            duwdxhmn (k) = duwdxhav (k,isamp)/nrsamph(k,isamp)
            dtaudxhmn(k) = dtaudxhav(k,isamp)/nrsamph(k,isamp)
            dtaudzhmn(k) = dtaudzhav(k,isamp)/nrsamph(k,isamp)
            thvhmn   (k) = thvhav   (k,isamp)/nrsamph(k,isamp)
            fcormn   (k) = fcorav   (k,isamp)/nrsamph(k,isamp)
            wwmn     (k) = wwav     (k,isamp)/nrsamph(k,isamp)
          endif

          nrsampmn  (k) = nrsamp   (k,isamp)/inorm
          massflxmn (k) = massflxav(k,isamp)/inorm
          wtlmn     (k) = wtlav    (k,isamp)/inorm
          wtvmn     (k) = wtvav    (k,isamp)/inorm
          wqtmn     (k) = wqtav    (k,isamp)/inorm
          wqlmn     (k) = wqlav    (k,isamp)/inorm
          uwmn      (k) = uwav     (k,isamp)/inorm
          vwmn      (k) = vwav     (k,isamp)/inorm

          nrsamphmn (k) = nrsamph  (k,isamp)/inorm


        enddo

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
           ,'  LEV  HGHT_F HGHT_H   PRES   COV_F  COV_H       W       THL      QT      ' &
           ,'QL       THV     P   WW_H '
        do k=1,kmax
          write(ifoutput,'(i5,2F6.0,F7.1,2F10.5,5F11.5,E14.5,F14.5)') &
              k, &
              zf      (k), &
              zh      (k), &
              presf   (k)/100., &
              nrsampmn(k), &
              nrsamphmn(k),&
              wmn     (k), &
              tlmn    (k), &
              qtmn    (k)*1000., &
              qlmn    (k)*1000., &
              tvmn    (k),&
              pmn     (k), &
              wwmn    (k)

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

        open (ifoutput,file=trim(samplname(isamp))//'wbudg.'//cexpnr,position='append')
        write(ifoutput,'(//3A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
          '#------------------------------ ',trim(samplname(isamp)),' -------------------------------'      &
          ,'#',timeav,'--- AVERAGING TIMESTEP --- '      &
          ,nhrs,':',nminut,':',nsecs      &
          ,'   HRS:MIN:SEC AFTER INITIALIZATION '


        write (ifoutput,'(2A/3A)') &
           '#------------------------------------------------------' &
           ,'------------------------------' &
           ,'  LEV HGHT   PRES     COVER     DWDTMN        BUO          DPDZMN       DWWDZHMN    DUWDXHMN ' &
           ,'    DTAUDZHMN     DTAUDXHMN    CORIOLIS     RESIDUAL    WS_END   SIG_END ' &
           ,'  '
        do k=1,kmax
          write(ifoutput,'(i5,F6.0,F7.1,F9.4,11F13.7)') &
              k, &
              zh      (k), &
              presh   (k)/100., &
              nrsamphmn(k),&
              dwdthmn  (k),&
              grav/thvs*thvhmn(k),&
              dpdzhmn  (k),&
              dwwdzhmn (k),&
              duwdxhmn (k),&
              dtaudzhmn(k),&
              dtaudxhmn(k),&
              fcormn   (k),&
              dwwdzhmn(k) +  grav/thvs*thvhmn(k) + dpdzhmn  (k) +  duwdxhmn(k) + dtaudxhmn(k) + dtaudzhmn(k) + &
                fcormn(k) -  dwdthmn  (k),&
              w_e(k,isamp),&
              sig_e(k,isamp)
        enddo
        close(ifoutput)

      if (lnetcdf) then
        vars(:, 1) =nrsampmn
        if (any(nrsampmn>0)) then
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
        else
          vars(:,2:13)=nc_fillvalue
        end if

        call writestat_nc(ncid_prof,nvar,ncname(:,:,isamp),vars(1:kmax,:),nrec_prof,kmax)
      end if

      end do
    end if
    deallocate( wmn,tlmn,tvmn,qtmn,qlmn,nrsampmn,massflxmn, &
                wtlmn,wtvmn,wqtmn,wqlmn,uwmn,vwmn)
    deallocate( wav,tlav,tvav,qtav,qlav,nrsamp,massflxav, &
                wtlav,wtvav,wqtav,wqlav,uwav,vwav)
    deallocate (nrsamphmn,wwmn,pmn,dwdthmn,dwwdzhmn,dpdzhmn,duwdxhmn,&
                 dtaudzhmn,dtaudxhmn,thvhmn,fcormn)
    deallocate (nrsamph,wwav)
    deallocate( pav,dwdthav,dwwdzhav,dpdzhav,duwdxhav,dtaudxhav,dtaudzhav,thvhav,fcorav)
    deallocate (w_e,sig_e)


  end subroutine writesampling

end module modsampling
