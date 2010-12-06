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
  integer,parameter :: nvar = 28
  character(80),allocatable,dimension(:,:,:) :: ncname
  real :: dtav, timeav
  integer(kind=longint) :: idtav,itimeav,tnext,tnextwrite
  integer :: nsamples,isamp,isamptot
  character(20),dimension(10) :: samplname,longsamplname
  logical :: lsampall = .false. !< switch for sampling (on/off)
  logical :: lsampcl  = .false. !< switch for conditional sampling cloud (on/off)
  logical :: lsampco  = .false. !< switch for conditional sampling core (on/off)
  logical :: lsampup  = .false. !< switch for conditional sampling updraft (on/off)
  logical :: lsampbuup  = .false. !< switch for conditional sampling buoyant updraft (on/off)
  logical :: lsampcldup = .false. !<switch for condtional sampling cloudy updraft (on/off)
  real, allocatable, dimension(:,:) ::  wfavl,thlfavl,thvfavl,qtfavl,qlfavl,nrsampfl,massflxhavl, &
                                        wtlthavl,wtvthavl,wqtthavl,wqlthavl,uwthavl,vwthavl
  real, allocatable, dimension(:,:) :: wwrhavl,wwshavl,pfavl,dwdthavl,dwwdzhavl,dpdzhavl, &
                                       duwdxhavl,dtaudxhavl,dtaudzhavl,thvhavl, &
                                       fcorhavl,nrsamphl
  real,allocatable, dimension(:,:) :: wh_el,sigh_el


contains
!> Initialization routine, reads namelists and inits variables
  subroutine initsampling

    use modmpi,    only : comm3d, my_real,mpierr,myid,mpi_logical
    use modglobal, only : ladaptive, dtmax,rk3step,k1,ifnamopt,fname_options,   &
                           dtav_glob,timeav_glob,dt_lim,btime,tres,cexpnr,ifoutput
    use modstat_nc, only : lnetcdf, redefine_nc,define_nc,ncinfo
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav,ncid_prof=>ncid
    implicit none

    integer :: ierr

    namelist/NAMSAMPLING/ &
    dtav,timeav,lsampcl,lsampco,lsampup,lsampbuup,lsampcldup

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


    call MPI_BCAST(timeav    ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dtav      ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(lsampall  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcl   ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampco   ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampup   ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampbuup ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lsampcldup,1,MPI_LOGICAL,0,comm3d,mpierr)

    isamptot = 0
    if (lsampall) then
      isamptot = isamptot + 1
      samplname (isamptot) = 'all'
      longsamplname(isamptot) = 'All '
    endif
    if (lsampup) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'upd'
      longsamplname(isamptot) = 'Updraft '
    end if
    if (lsampbuup) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'buup'
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
    if (lsampcldup) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cldup'
      longsamplname(isamptot) = 'Cloud Updraft '
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

    allocate( wfavl     (k1,isamptot),thlfavl  (k1,isamptot),thvfavl   (k1,isamptot), &
              qtfavl    (k1,isamptot),qlfavl   (k1,isamptot),nrsampfl  (k1,isamptot), &
              wtlthavl  (k1,isamptot),wtvthavl (k1,isamptot),wqtthavl  (k1,isamptot), &
              wqlthavl  (k1,isamptot),uwthavl  (k1,isamptot),vwthavl   (k1,isamptot))

    allocate( nrsamphl  (k1,isamptot),wwrhavl  (k1,isamptot),wwshavl   (k1,isamptot), &
              pfavl     (k1,isamptot),dwdthavl (k1,isamptot),dwwdzhavl (k1,isamptot), &
              dpdzhavl  (k1,isamptot),duwdxhavl(k1,isamptot),dtaudxhavl(k1,isamptot), &
              dtaudzhavl(k1,isamptot),thvhavl  (k1,isamptot),fcorhavl  (k1,isamptot))

    allocate (wh_el     (k1,isamptot),sigh_el  (k1,isamptot),massflxhavl(k1,isamptot))


!initialize variables
    nrsampfl    = 0.0
    wfavl       = 0.0
    thlfavl     = 0.0
    thvfavl     = 0.0
    qtfavl      = 0.0
    qlfavl      = 0.0
    massflxhavl = 0.0
    wtlthavl    = 0.0
    wtvthavl    = 0.0
    wqtthavl    = 0.0
    wqlthavl    = 0.0
    uwthavl     = 0.0
    vwthavl     = 0.0

    nrsamphl    = 0.0
    wwrhavl     = 0.0
    wwshavl     = 0.0
    pfavl       = 0.0
    dwdthavl    = 0.0
    dwwdzhavl   = 0.0
    dpdzhavl    = 0.0
    duwdxhavl   = 0.0
    dtaudxhavl  = 0.0
    dtaudzhavl  = 0.0
    thvhavl     = 0.0
    fcorhavl    = 0.0

    wh_el       = 0.0
    sigh_el     = 0.0

    if(myid==0)then
      do isamp = 1,isamptot
        open (ifoutput,file=trim(samplname(isamp))//'wbudg.'//cexpnr,status='replace')
        close (ifoutput)
        open (ifoutput,file=trim(samplname(isamp))//'fld.'//cexpnr,status='replace')
        close (ifoutput)
        open (ifoutput,file=trim(samplname(isamp))//'flx.'//cexpnr,status='replace')
        close (ifoutput)
      enddo
    endif

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
          call ncinfo(ncname(14,:,isamp),'nrsamph'//samplname(isamp),trim(longsamplname(isamp))//' '//'number of points at halflevel','-','mt')
          call ncinfo(ncname(15,:,isamp),'pf'//samplname(isamp),trim(longsamplname(isamp))//' '// 'pressure','kg/m/s^2','tt')
          call ncinfo(ncname(16,:,isamp),'wwrh'//samplname(isamp),trim(longsamplname(isamp))//' '//'ww res flux','m^2/s^2','mt')
          call ncinfo(ncname(17,:,isamp),'wwsh'//samplname(isamp),trim(longsamplname(isamp))//' '//'ww sub flux','m^2/s^2','mt')
          call ncinfo(ncname(18,:,isamp),'dwdth'//samplname(isamp),trim(longsamplname(isamp))//' '//'dwdt at sampled point','m/s^2','mt')
          call ncinfo(ncname(19,:,isamp),'buoyh'//samplname(isamp),trim(longsamplname(isamp))//' '//'buoyancy force','m/s^2','mt')
          call ncinfo(ncname(20,:,isamp),'dpdzh'//samplname(isamp),trim(longsamplname(isamp))//' '//'vertical pressure force','m/s^2','mt')
          call ncinfo(ncname(21,:,isamp),'dwwdzh'//samplname(isamp),trim(longsamplname(isamp))//' '//'resolved vertical w advection','m/s^2','mt')
          call ncinfo(ncname(22,:,isamp),'duwdxh'//samplname(isamp),trim(longsamplname(isamp))//' '//'resolved horizontal w advection','m/s^2','mt')
          call ncinfo(ncname(23,:,isamp),'dtaudzh'//samplname(isamp),trim(longsamplname(isamp))//' '//'sg vertical w advection','m/s^2','mt')
          call ncinfo(ncname(24,:,isamp),'dtaudxh'//samplname(isamp),trim(longsamplname(isamp))//' '//'sg horizontal w advection','m/s^2','mt')
          call ncinfo(ncname(25,:,isamp),'fcorh'//samplname(isamp),trim(longsamplname(isamp))//' '//'coriolis force on w','m/s^2','mt')
          call ncinfo(ncname(26,:,isamp),'resid'//samplname(isamp),trim(longsamplname(isamp))//' '//'residual term in sampled budget eqn','m/s^2','mt')
          call ncinfo(ncname(27,:,isamp),'whend'//samplname(isamp),trim(longsamplname(isamp))//' '//'ws at end of sampling period','m/s','mt')
          call ncinfo(ncname(28,:,isamp),'sighend'//samplname(isamp),trim(longsamplname(isamp))//' '//'sigma at end of period','-','mt')
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
       deallocate( wfavl   ,thlfavl ,thvfavl ,qtfavl   ,qlfavl   ,nrsampfl  ,massflxhavl, &
                   wtlthavl,wtvthavl,wqtthavl,wqlthavl ,uwthavl  ,vwthavl)
       deallocate( nrsamphl,wwrhavl ,wwshavl , &
                   pfavl   ,dwdthavl,dwwdzhavl,dpdzhavl,duwdxhavl,dtaudxhavl,dtaudzhavl,  &
                   thvhavl ,fcorhavl,wh_el,sigh_el) 
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
                          wp_store
    use modsubgrid,only : ekh,ekm
    use modmpi    ,only : slabsum
    use modpois,   only : p
    use modsurfdata,only: thvs
    implicit none

    logical, allocatable, dimension(:,:,:) :: maskf
    real, allocatable, dimension(:,:,:) :: wtlth,wqtth,wqlth,wtvth,uwth,vwth
    real, allocatable, dimension(:,:,:) :: w0f
    real, allocatable, dimension(:,:,:) :: thv0
    real, allocatable, dimension(:) :: thvav

    logical, allocatable, dimension(:,:,:) :: maskh
    real, allocatable, dimension(:,:,:) :: uwsh,vwsh,uwrh,vwrh,wwrh,wwsf
    real, allocatable, dimension(:) :: thvhav


    integer :: i,j,k,km,kp
    real :: cqt,cthl,den,ekhalf,c2,c1,t0h,qs0h,ekav
    real :: wtvsh,wtvrh,wtlsh,wtlrh,wqtrh,wqtsh,wqlrh,wqls
    real :: beta

    allocate(thvav(k1))
    allocate( maskf(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate( wtlth(2-ih:i1+ih,2-jh:j1+jh,k1),&
              wtvth(2-ih:i1+ih,2-jh:j1+jh,k1),&
              wqtth(2-ih:i1+ih,2-jh:j1+jh,k1),&
              wqlth(2-ih:i1+ih,2-jh:j1+jh,k1),&
              uwth (2-ih:i1+ih,2-jh:j1+jh,k1),&
              vwth (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate( thv0 (2-ih:i1+ih,2-jh:j1+jh,k1),&
              w0f  (2-ih:i1+ih,2-jh:j1+jh,k1))      !Change to 3D array
    allocate(maskh (2-ih:i1+ih,2-jh:j1+jh,k1))     !New
    allocate(uwsh  (2-ih:i1+ih,2-jh:j1+jh,k1),&
             vwsh  (2-ih:i1+ih,2-jh:j1+jh,k1),&
             uwrh  (2-ih:i1+ih,2-jh:j1+jh,k1),&
             vwrh  (2-ih:i1+ih,2-jh:j1+jh,k1),&
             wwrh  (2-ih:i1+ih,2-jh:j1+jh,k1),&
             wwsf  (2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thvhav(k1))
    
    beta = grav/thvs

    do k=1,k1
      thv0(2:i1,2:j1,k) = (thl0(2:i1,2:j1,k)+rlv*ql0(2:i1,2:j1,k)/(cp*exnf(k))) &
                        *(1+(rv/rd-1)*qt0(2:i1,2:j1,k)-rv/rd*ql0(2:i1,2:j1,k))
    enddo
    do k=1,kmax
      w0f (2:i1,2:j1,k) = 0.5*(w0 (2:i1,2:j1,k) + w0  (2:i1,2:j1,k+1))
    end do

    maskf = .false.
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
             maskf(i,j,k) = .true.
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
             maskf(i,j,k) = .true.
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
             maskf(i,j,k) = .true.
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
             maskf(i,j,k) = .true.
         endif
         if (ql0h(i,j,k)>epsilon(1.0).and.thv0h(i,j,k) > thvhav(k)) then
             maskh(i,j,k) = .true.
         endif
      enddo
      enddo
      enddo

    case ('cldup')
      do i=2,i1
      do j=2,j1
      do k=1,kmax
         if (ql0(i,j,k)>epsilon(1.0).and.w0f(i,j,k).gt.0.) then
             maskf(i,j,k) = .true.
         endif
         if (ql0h(i,j,k)>epsilon(1.0).and.w0(i,j,k).gt.0.) then
             maskh(i,j,k) = .true.
         endif
      enddo
      enddo
      enddo


    case ('all')
        maskf  = .true.
        maskh = .true.
    end select


!calculate fluxes at half levels
    wtlth = 0.0
    wtvth = 0.0
    wqtth = 0.0
    wqlth = 0.0
    uwth  = 0.0
    vwth  = 0.0

    uwsh  = 0.0
    vwsh  = 0.0
    uwrh  = 0.0
    vwrh  = 0.0
    wwrh  = 0.0
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

      wtlsh   = -ekhalf*(thl0(i,j,k)-thl0(i,j,km))/dzh(k)
      wtlrh   = w0(i,j,k)*thl0h(i,j,k)
      wtlth (i,j,k)   = wtlsh + wtlrh

      wqtsh   = -ekhalf*(qt0(i,j,k)-qt0(i,j,km))/dzh(k)
      wqtrh   = w0(i,j,k)*qt0h(i,j,k)
      wqtth(i,j,k)   = wqtsh + wqtrh

      if (ql0h(i,j,k)>0) then
        wqls   = cthl*wtlsh+ cqt*wqtsh
      else
        wqls   = 0.
      end if

      wqlrh   = w0(i,j,k)*ql0h(i,j,k)
      wqlth(i,j,k)   = wqlrh + wqls

      wtvsh   = c1*wtlsh + c2*thl0h(i,j,k)*wqtsh
      wtvrh   = w0(i,j,k)*thv0h(i,j,k)
      wtvth(i,j,k)   = wtvrh + wtvsh

      ekav = 0.25 * &
              ( ekm(i,j,k-1)+ekm(i,j,k)+ekm(i-1,j,k-1)+ekm(i-1,j,k) )
      uwsh(i,j,k)  =  - ekav * &
            ((u0(i,j,k)-u0(i,j,k-1))/dzf(k)+(w0(i,j,k)-w0(i-1,j,k))/dx)
      uwrh(i,j,k) = (w0(i,j,k)+w0(i-1,j,k))*(u0(i,j,k-1)+u0(i,j,k))/4.
      uwth(i,j,k) = uwsh(i,j,k) + uwrh(i,j,k)

      ekav = 0.25 * &
          ( ekm(i,j,k-1)+ekm(i,j,k)+ekm(i,j-1,k-1)+ekm(i,j-1,k) )
      vwsh(i,j,k)  = - ekav * &
          ((v0(i,j,k)-v0(i,j,k-1))/dzf(k)+(w0(i,j,k)-w0(i,j-1,k))/dy)
      vwrh(i,j,k) = (w0(i,j,k)+w0(i,j-1,k))*(v0(i,j,k-1)+v0(i,j,k))/4.
      vwth(i,j,k) = vwsh(i,j,k) + vwrh(i,j,k)

      wwrh  (i,j,k) = w0(i,j,k)**2
    end do
    end do
    end do

    do k=2,kmax
       i=i2
       do j=2,j1
          ekav = 0.25 * &
               ( ekm(i,j,k-1)+ekm(i,j,k)+ekm(i-1,j,k-1)+ekm(i-1,j,k) )
          uwsh(i,j,k)  =  - ekav * &
            ((u0(i,j,k)-u0(i,j,k-1))/dzf(k)+(w0(i,j,k)-w0(i-1,j,k))/dx)
          uwrh (i,j,k) = (w0(i,j,k)+w0(i-1,j,k))*(u0(i,j,k-1)+u0(i,j,k))/4.
       end do
       j=j2
       do i=2,i1
          ekav = 0.25 * &
             ( ekm(i,j,k-1)+ekm(i,j,k)+ekm(i,j-1,k-1)+ekm(i,j-1,k) )
          vwsh (i,j,k) = - ekav * &
             ((v0(i,j,k)-v0(i,j,k-1))/dzf(k)+(w0(i,j,k)-w0(i,j-1,k))/dy)
          vwrh (i,j,k) = (w0(i,j,k)+w0(i,j-1,k))*(v0(i,j,k-1)+v0(i,j,k))/4.
       end do
    end do

!
!  full levels
!
    do j=2,j1
    do i=2,i1
    do k=2,kmax
      km = k-1
      wwsf (i,j,km) = -2 * ekm(i,j,km) * &
               (w0(i,j,k)-w0(i,j,km)) / dzf(km)

    end do
    end do
    end do


!add fields and fluxes to mean
!     1)       fields on full levels
    do k=1,kmax
      nrsampfl(k,isamp) = nrsampfl(k,isamp)+count(maskf(2:i1,2:j1,k))
      nrsamphl(k,isamp) = nrsamphl(k,isamp)+count(maskh(2:i1,2:j1,k))
      wfavl   (k,isamp) = wfavl   (k,isamp)+sum  (w0f  (2:i1,2:j1,k),maskf(2:i1,2:j1,k))
      thlfavl (k,isamp) = thlfavl (k,isamp)+sum  (thl0 (2:i1,2:j1,k),maskf(2:i1,2:j1,k))
      thvfavl (k,isamp) = thvfavl (k,isamp)+sum  (thv0 (2:i1,2:j1,k),maskf(2:i1,2:j1,k))
      qtfavl  (k,isamp) = qtfavl  (k,isamp)+sum  (qt0  (2:i1,2:j1,k),maskf(2:i1,2:j1,k))
      qlfavl  (k,isamp) = qlfavl  (k,isamp)+sum  (ql0  (2:i1,2:j1,k),maskf(2:i1,2:j1,k))
      pfavl   (k,isamp) = pfavl   (k,isamp)+sum  (p    (2:i1,2:j1,k),maskf(2:i1,2:j1,k))
    end do

!     2)       fluxes on half levels

    do k=2,kmax
      massflxhavl(k,isamp) = massflxhavl(k,isamp) + sum(w0   (2:i1,2:j1,k),maskh(2:i1,2:j1,k)) 
      wtlthavl   (k,isamp) = wtlthavl   (k,isamp) + sum(wtlth(2:i1,2:j1,k),maskh(2:i1,2:j1,k)) 
      wtvthavl   (k,isamp) = wtvthavl   (k,isamp) + sum(wtvth(2:i1,2:j1,k),maskh(2:i1,2:j1,k)) 
      wqtthavl   (k,isamp) = wqtthavl   (k,isamp) + sum(wqtth(2:i1,2:j1,k),maskh(2:i1,2:j1,k)) 
      wqlthavl   (k,isamp) = wqlthavl   (k,isamp) + sum(wqlth(2:i1,2:j1,k),maskh(2:i1,2:j1,k)) 
      uwthavl    (k,isamp) = uwthavl    (k,isamp) + sum(uwth (2:i1,2:j1,k),maskh(2:i1,2:j1,k)) 
      vwthavl    (k,isamp) = vwthavl    (k,isamp) + sum(vwth (2:i1,2:j1,k),maskh(2:i1,2:j1,k)) 
      wwrhavl    (k,isamp) = wwrhavl    (k,isamp) + sum(wwrh (2:i1,2:j1,k),maskh(2:i1,2:j1,k))
      dwdthavl   (k,isamp) = dwdthavl   (k,isamp) + sum(wp_store(2:i1,2:j1,k),maskh(2:i1,2:j1,k))
      wwshavl    (k,isamp) = wwshavl    (k,isamp) +(sum(wwsf (2:i1,2:j1,k),maskf(2:i1,2:j1,k))  & 
                                                  + sum(wwsf (2:i1,2:j1,k),maskf(2:i1,2:j1,k-1))) * 0.5
                                            
      wh_el      (k,isamp) =                        sum(w0 (2:i1,2:j1,k),maskh(2:i1,2:j1,k))
      sigh_el    (k,isamp) =                        count(maskh(2:i1,2:j1,k))


      do j=2,j1
      do i=2,i1
       if (maskh(i,j,k)) then
         thvhavl   (k,isamp) = thvhavl   (k,isamp) + thv0h(i,j,k) - thvhav(k)
         dpdzhavl  (k,isamp) = dpdzhavl  (k,isamp) - (p(i,j,k)-p(i,j,k-1))/dzh(k) &
                                                   + beta * thvhav(k)
         dwwdzhavl (k,isamp) = dwwdzhavl (k,isamp) - (w0f(i,j,k)**2- w0f(i,j,k-1)**2)/dzh(k)
         duwdxhavl (k,isamp) = duwdxhavl (k,isamp) - (uwrh(i+1,j,k) - uwrh (i,j,k))/dx - &
                                                     (vwrh(i,j+1,k) - vwrh (i,j,k))/dy
         dtaudxhavl(k,isamp) = dtaudxhavl(k,isamp) - (uwsh(i+1,j,k) - uwsh (i,j,k))/dx - &
                                                     (vwsh(i,j+1,k) - vwsh (i,j,k))/dy   
         dtaudzhavl(k,isamp) = dtaudzhavl(k,isamp) - (wwsf(i,j,k)   - wwsf(i,j,k-1))/dzh(k)
         fcorhavl  (k,isamp) = fcorhavl(k,isamp) + om22 * cu  &
                                    +( (dzf(k-1) * (u0(i,j,k)   + u0(i+1,j,k) )    &
                                    +    dzf(k)  * (u0(i,j,k-1) + u0(i+1,j,k-1))  ) / dzh(k) ) &
                                    * om22*0.25
       endif
      end do
      end do
    enddo



    deallocate(maskf,wtlth,wqtth,wqlth,wtvth,uwth,vwth,thvav,w0f,thv0)
    deallocate(maskh,uwsh,vwsh,uwrh,vwrh)
    deallocate(wwrh,wwsf,thvhav)



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

    real, allocatable, dimension(:)  :: wfmn,thlfmn,thvfmn,qtfmn,qlfmn,nrsampfmn,massflxhmn, &
                                        wtlthmn,wtvthmn,wqtthmn,wqlthmn,uwthmn,vwthmn
    real, allocatable, dimension(:,:):: wfav,thlfav,thvfav,qtfav,qlfav,nrsampf,massflxhav, &
                                        wtlthav,wtvthav,wqtthav,wqlthav,uwthav,vwthav
    real, allocatable, dimension(:)  :: nrsamphmn, wwrhmn,wwshmn,&
                                        pfmn,dwdthmn,dwwdzhmn,dpdzhmn,duwdxhmn,&
                                        dtaudxhmn,dtaudzhmn,thvhmn, &
                                        fcorhmn
    real, allocatable, dimension(:,:):: nrsamph,wwrhav,wwshav, & 
                                        pfav,dwdthav,dwwdzhav,dpdzhav,duwdxhav,&
                                        dtaudxhav,dtaudzhav,thvhav,&
                                        fcorhav,wh_e,sigh_e

    integer :: nsecs, nhrs, nminut, k
    integer :: inorm
    allocate( wfmn(k1),thlfmn(k1),thvfmn(k1),qtfmn(k1),qlfmn(k1),nrsampfmn(k1),massflxhmn(k1), &
              wtlthmn(k1),wtvthmn(k1),wqtthmn(k1),wqlthmn(k1),uwthmn(k1),vwthmn(k1))
    allocate( wfav(k1,isamptot),thlfav(k1,isamptot),thvfav(k1,isamptot), &
              qtfav(k1,isamptot),qlfav(k1,isamptot),nrsampf(k1,isamptot),massflxhav(k1,isamptot), &
              wtlthav(k1,isamptot),wtvthav(k1,isamptot),wqtthav(k1,isamptot), &
              wqlthav(k1,isamptot),uwthav(k1,isamptot),vwthav(k1,isamptot))
    allocate (nrsamphmn(k1),wwrhmn(k1),wwshmn(k1))
    allocate( pfmn(k1),dwdthmn(k1),dwwdzhmn(k1),dpdzhmn(k1),duwdxhmn(k1),&
              dtaudxhmn(k1),dtaudzhmn(k1),thvhmn(k1),fcorhmn(k1))
    allocate( nrsamph(k1,isamptot),wwrhav(k1,isamptot),wwshav(k1,isamptot))
    allocate( pfav(k1,isamptot),dwdthav(k1,isamptot),dwwdzhav(k1,isamptot),&
              dpdzhav(k1,isamptot),duwdxhav(k1,isamptot),&
              dtaudxhav(k1,isamptot),dtaudzhav(k1,isamptot),thvhav(k1,isamptot),fcorhav(k1,isamptot))
    allocate (wh_e(k1,isamptot),sigh_e(k1,isamptot))

    nsecs   = nint(rtimee)
    nhrs    = int(nsecs/3600)
    nminut  = int(nsecs/60)-nhrs*60
    nsecs   = mod(nsecs,60)
    inorm   = nint(rslabs*timeav/dtav)

    call MPI_ALLREDUCE(nrsampfl   ,nrsampf    ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wfavl      ,wfav       ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(thlfavl    ,thlfav     ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(thvfavl    ,thvfav     ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(qtfavl     ,qtfav      ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(qlfavl     ,qlfav      ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(massflxhavl,massflxhav ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wtlthavl   ,wtlthav    ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wtvthavl   ,wtvthav    ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wqtthavl   ,wqtthav    ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wqlthavl   ,wqlthav    ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(uwthavl    ,uwthav     ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(vwthavl    ,vwthav     ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)

    call MPI_ALLREDUCE(nrsamphl   ,nrsamph    ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wwrhavl    ,wwrhav     ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(wwshavl    ,wwshav     ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(pfavl      ,pfav       ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(dwdthavl   ,dwdthav    ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(dwwdzhavl  ,dwwdzhav   ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(dpdzhavl   ,dpdzhav    ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(duwdxhavl  ,duwdxhav   ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(dtaudxhavl ,dtaudxhav  ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(dtaudzhavl ,dtaudzhav  ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(thvhavl    ,thvhav     ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(fcorhavl   ,fcorhav    ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)

    call MPI_ALLREDUCE(wh_el      ,wh_e       ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(sigh_el    ,sigh_e     ,isamptot*k1,MY_REAL,MPI_SUM,comm3d,mpierr)



!reset variables
    nrsampfl    = 0.0
    wfavl       = 0.0
    thlfavl     = 0.0
    thvfavl     = 0.0
    qtfavl      = 0.0
    qlfavl      = 0.0
    massflxhavl = 0.0
    wtlthavl    = 0.0
    wtvthavl    = 0.0
    wqtthavl    = 0.0
    wqlthavl    = 0.0
    uwthavl     = 0.0
    vwthavl     = 0.0

    nrsamphl    = 0.0
    wwrhavl     = 0.0
    wwshavl     = 0.0
    pfavl       = 0.0
    dwdthavl    = 0.0
    dwwdzhavl   = 0.0
    dpdzhavl    = 0.0
    duwdxhavl   = 0.0
    dtaudxhavl  = 0.0
    dtaudzhavl  = 0.0
    thvhavl     = 0.0
    fcorhavl    = 0.0

    wh_el       = 0.
    sigh_el     = 0.


    if (myid==0) then
      do isamp = 1,isamptot

         wfmn     = 0.
         thlfmn   = 0.
         thvfmn   = 0.
         qtfmn    = 0.
         qlfmn    = 0.

         pfmn     = 0.
         dwdthmn  = 0.
         dwwdzhmn = 0.
         dpdzhmn  = 0.
         duwdxhmn = 0.
         dtaudxhmn= 0.
         dtaudzhmn= 0.
         thvhmn   = 0.
         fcorhmn  = 0.
         wwrhmn   = 0.
         wwshmn   = 0.

!normalize variables

        do k=1,kmax
          if (nrsampf(k,isamp)>0) then
            wfmn   (k) = wfav  (k,isamp)/nrsampf(k,isamp)
            thlfmn (k) = thlfav(k,isamp)/nrsampf(k,isamp)
            thvfmn (k) = thvfav(k,isamp)/nrsampf(k,isamp)
            qtfmn  (k) = qtfav (k,isamp)/nrsampf(k,isamp)
            qlfmn  (k) = qlfav (k,isamp)/nrsampf(k,isamp)
          endif

          if (sigh_e(k,isamp).gt.0) then
            wh_e  (k,isamp) = wh_e  (k,isamp) / sigh_e(k,isamp)
            sigh_e(k,isamp) = sigh_e(k,isamp) / rslabs
          endif

          if (nrsamph(k,isamp).gt.0) then
            pfmn     (k) = pfav      (k,isamp)/nrsamph(k,isamp)
            dwdthmn  (k) = dwdthav   (k,isamp)/nrsamph(k,isamp)
            dwwdzhmn (k) = dwwdzhav  (k,isamp)/nrsamph(k,isamp)
            dpdzhmn  (k) = dpdzhav   (k,isamp)/nrsamph(k,isamp)
            duwdxhmn (k) = duwdxhav  (k,isamp)/nrsamph(k,isamp)
            dtaudxhmn(k) = dtaudxhav (k,isamp)/nrsamph(k,isamp)
            dtaudzhmn(k) = dtaudzhav (k,isamp)/nrsamph(k,isamp)
            thvhmn   (k) = thvhav    (k,isamp)/nrsamph(k,isamp)
            fcorhmn  (k) = fcorhav   (k,isamp)/nrsamph(k,isamp)
            wwrhmn   (k) = wwrhav    (k,isamp)/nrsamph(k,isamp)
            wwshmn   (k) = wwshav    (k,isamp)/nrsamph(k,isamp)
          endif

          nrsampfmn  (k) = nrsampf   (k,isamp)/inorm
          massflxhmn (k) = massflxhav(k,isamp)/inorm
          wtlthmn    (k) = wtlthav   (k,isamp)/inorm
          wtvthmn    (k) = wtvthav   (k,isamp)/inorm
          wqtthmn    (k) = wqtthav   (k,isamp)/inorm
          wqlthmn    (k) = wqlthav   (k,isamp)/inorm
          uwthmn     (k) = uwthav    (k,isamp)/inorm
          vwthmn     (k) = vwthav    (k,isamp)/inorm

          nrsamphmn  (k) = nrsamph   (k,isamp)/inorm


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
           ,'QL       THV     P   WW_RES_H   WW_SUB_H'
        do k=1,kmax
          write(ifoutput,'(i5,2F6.0,F7.1,2F10.5,5F11.5,E14.5,2F14.5)') &
              k, &
              zf       (k), &
              zh       (k), &
              presf    (k)/100., &
              nrsampfmn(k), &
              nrsamphmn(k),&
              wfmn     (k), &
              thlfmn   (k), &
              qtfmn    (k)*1000., &
              qlfmn    (k)*1000., &
              thvfmn   (k),&
              pfmn     (k), &
              wwrhmn   (k),&
              wwshmn   (k)

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
                zh         (k), &
                presh      (k)/100., &
                massflxhmn (k), &
                wtlthmn    (k), &
                wqtthmn    (k), &
                wqlthmn    (k), &
                wtvthmn    (k), &
                uwthmn     (k), &
                vwthmn     (k)
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
              zh        (k), &
              presh     (k)/100., &
              nrsamphmn (k),&
              dwdthmn   (k),&
              grav/thvs*thvhmn(k),&
              dpdzhmn   (k),&
              dwwdzhmn  (k),&
              duwdxhmn  (k),&
              dtaudzhmn (k),&
              dtaudxhmn (k),&
              fcorhmn   (k),&
              dwwdzhmn(k) +  grav/thvs*thvhmn(k) + dpdzhmn  (k) +  duwdxhmn(k) + dtaudxhmn(k) + dtaudzhmn(k) + &
                fcorhmn(k) -  dwdthmn  (k),&
              wh_e      (k,isamp),&
              sigh_e    (k,isamp)
        enddo
        close(ifoutput)

      if (lnetcdf) then
        vars(:, 1)=nrsampfmn
        if (any(nrsampfmn>0)) then
        vars(:, 2)=wfmn
        vars(:, 3)=thlfmn
        vars(:, 4)=qtfmn
        vars(:, 5)=qlfmn
        vars(:, 6)=thvfmn
        vars(:, 7)=massflxhmn
        vars(:, 8)=wtlthmn
        vars(:, 9)=wqtthmn
        vars(:,10)=wqlthmn
        vars(:,11)=wtvthmn
        vars(:,12)=uwthmn
        vars(:,13)=vwthmn
        vars(:,14)=nrsamphmn
        vars(:,15)=pfmn
        vars(:,16)=wwrhmn
        vars(:,17)=wwshmn
        vars(:,18)=dwdthmn
        vars(:,19)=grav/thvs*thvhmn
        vars(:,20)=dpdzhmn
        vars(:,21)=dwwdzhmn
        vars(:,22)=duwdxhmn
        vars(:,23)=dtaudzhmn
        vars(:,24)=dtaudxhmn
        vars(:,25)=fcorhmn
        vars(:,26)=dwwdzhmn(:)+grav/thvs*thvhmn(:)+dpdzhmn(:)+duwdxhmn(:)+dtaudxhmn(:)+dtaudzhmn(:)+fcorhmn(:)-dwdthmn(:)
        vars(:,27)=wh_e(:,isamp)
        vars(:,28)=sigh_e(:,isamp)
        else
          vars(:,2:28)=nc_fillvalue
        end if
        call writestat_nc(ncid_prof,nvar,ncname(:,:,isamp),vars(1:kmax,:),nrec_prof,kmax)
      end if

      end do
    end if
    deallocate( wfmn,thlfmn,thvfmn,qtfmn,qlfmn,nrsampfmn,massflxhmn, &
                wtlthmn,wtvthmn,wqtthmn,wqlthmn,uwthmn,vwthmn)
    deallocate( wfav,thlfav,thvfav,qtfav,qlfav,nrsampf,massflxhav, &
                wtlthav,wtvthav,wqtthav,wqlthav,uwthav,vwthav)
    deallocate (nrsamphmn,wwrhmn,wwshmn,pfmn,dwdthmn,dwwdzhmn,dpdzhmn,duwdxhmn,&
                 dtaudzhmn,dtaudxhmn,thvhmn,fcorhmn)
    deallocate (nrsamph,wwrhav,wwshav)
    deallocate( pfav,dwdthav,dwwdzhav,dpdzhav,duwdxhav,dtaudxhav,dtaudzhav,thvhav,fcorhav)
    deallocate (wh_e,sigh_e)


  end subroutine writesampling

end module modsampling
