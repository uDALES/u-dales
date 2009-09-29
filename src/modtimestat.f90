!> \file modtimestat.f90
!!  Timestat calculates timeseries of several variables

!>
!! Timestat calculates timeseries of several variables
!>
!! Timeseries of the most relevant parameters. Written to tmser1.expnr and tmsurf.expnr
!! If netcdf is true, this module leads the tmser.expnr.nc output
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
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


module modtimestat



implicit none
! private
! PUBLIC :: inittimestat, timestat
save
!NetCDF variables
  integer,parameter :: nvar = 21
  integer :: ncid,nrec = 0
  character(80) :: fname = 'tmser.xxx.nc'
  character(80),dimension(nvar,4) :: ncname

  real    :: dtav,tnext
  logical :: ltimestat= .false. !<switch for timestatistics (on/off)
  real    :: zi,ziold=-1, we
  integer, parameter :: iblh_flux = 1, iblh_grad = 2, iblh_thres = 3
  integer, parameter :: iblh_thv = -1,iblh_thl = -2, iblh_qt = -3
  integer :: iblh_meth = iblh_grad, iblh_var = iblh_thv
  integer :: blh_nsamp = 4
  real    :: blh_thres=-1 ,blh_sign=1.0
  real   :: zbaseav, ztopav, ztopmax,zbasemin
  real   :: qlintav, qlintmax, tke_tot
  real   :: cc, wmax, qlmax
  real   :: qlint

contains
!> Initializing Timestat. Read out the namelist, initializing the variables
  subroutine inittimestat
    use modmpi,    only : my_real,myid,comm3d,mpi_logical,mpierr,mpi_integer
    use modglobal, only : ifnamopt, fname_options,cexpnr,dtmax,ifoutput,dtav_glob,ladaptive,k1,kmax,rd,rv,dt_lim,btime
    use modfields, only : thlprof,qtprof,svprof
    use modsurfdata, only : isurf
    use modstat_nc, only : lnetcdf, open_nc,define_nc,ncinfo
    implicit none
    integer :: ierr,k,location = 1
    real :: gradient = 0.0
    real, allocatable,dimension(:) :: profile


    namelist/NAMTIMESTAT/ & !< namelist
    dtav,ltimestat,blh_thres,iblh_meth,iblh_var,blh_nsamp,blh_thres !! namelist contents
!!bla
!!dibla

    dtav=dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMTIMESTAT,iostat=ierr)
      write(6 ,NAMTIMESTAT)
      close(ifnamopt)
    end if

    call MPI_BCAST(dtav       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ltimestat  ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(blh_thres  ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(iblh_meth  ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(iblh_var   ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(blh_nsamp  ,1,MPI_INTEGER,0,comm3d,mpierr)

    tnext = dtav-1e-3+btime

    if(.not.(ltimestat)) return
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'TIMESTAT: dtav should be a integer multiple of dtmax'
    end if

    allocate(profile(1:k1))
    select case (iblh_var)
    case(iblh_qt)
      profile = qtprof
    case(iblh_thl)
      profile = thlprof
    case(iblh_thv)
      do k=1,k1
        profile(k) = thlprof(k)*(1+(rv/rd-1)*qtprof(k))
      end do
    case(1:)
      profile = svprof(:,iblh_var)
    end select
    blh_sign = sign(1.0,profile(kmax)-profile(1))

    select case(iblh_meth)
    case (iblh_flux)
    case (iblh_grad)
    case (iblh_thres)
      if (blh_thres<0) then
        do k=kmax,2,-1
          if (blh_sign*(profile(k+1) - profile(k-1)) > gradient) then
            location = k
            gradient = blh_sign*(profile(k+1) - profile(k-1))
          endif
        enddo
        blh_thres=profile(location)
        if (myid==0) write (*,*) 'TIMESTAT: blh_tres =',blh_thres
      end if
    case default
      stop 'TIMESTAT: Incorrect iblh_meth'
    end select
    deallocate(profile)
    if(myid==0) then
      !tmser1
      open (ifoutput,file='tmser1.'//cexpnr,status='replace',position='append')
      write(ifoutput,'(2a)') &
             '#  time      cc     z_cbase    z_ctop_avg  z_ctop_max      zi         we', &
             '   <<ql>>  <<ql>>_max   w_max   tke     ql_max'
      close(ifoutput)
      !tmsurf
      open (ifoutput,file='tmsurf.'//cexpnr,status='replace',position='append')
      write(ifoutput,'(2a)') &
             '#  time        ust        tst        qst         obukh', &
             '      thls        z0        wthls      wthvs      wqls '
      close(ifoutput)
      if(isurf == 1) then
        open (ifoutput,file='tmlsm.'//cexpnr,status='replace',position='append')
        write(ifoutput,'(2a)') &
               '#     time      Qnet        H          LE         G0  ', &
               '   tendskin       rs         ra        tskin'
        write(ifoutput,'(2a)') &
               '#      [s]     [W/m2]     [W/m2]     [W/m2]     [W/m2]  ', &
               '   [W/m2]      [s/m]      [s/m]        [K]'
        close(ifoutput)
      end if
      if (lnetcdf) then
        fname(7:9) = cexpnr
        call ncinfo(ncname( 1,:),'time','Time','s','time')
        call ncinfo(ncname( 2,:),'cfrac','Cloud fraction','-','time')
        call ncinfo(ncname( 3,:),'zb','Cloud-base height','m','time')
        call ncinfo(ncname( 4,:),'zc_av','Average Cloud-top height','m','time')
        call ncinfo(ncname( 5,:),'zc_max','Maximum Cloud-top height','m','time')
        call ncinfo(ncname( 6,:),'zi','Boundary layer height','m','time')
        call ncinfo(ncname( 7,:),'we','Entrainment velocity','m/s','time')
        call ncinfo(ncname( 8,:),'lwp_bar','Liquid-water path','kg/m^2','time')
        call ncinfo(ncname( 9,:),'lwp_max','Maximum Liquid-water path','kg/m^2','time')
        call ncinfo(ncname(10,:),'wmax','Maximum vertical velocity','m/s','time')
        call ncinfo(ncname(11,:),'vtke','Vertical integral of total TKE','kg/s','time')
        call ncinfo(ncname(12,:),'lmax','Maximum liquid water mixing ratio','kg/kg','time')
        call ncinfo(ncname(13,:),'ustar','Surface friction velocity','m/s','time')
        call ncinfo(ncname(14,:),'tstr','Turbulent temperature scale','K','time')
        call ncinfo(ncname(15,:),'qtstr','Turbulent humidity scale','K','time')
        call ncinfo(ncname(16,:),'obukh','Obukhov Length','m','time')
        call ncinfo(ncname(17,:),'tsrf','Surface liquid water potential temperature','K','time')
        call ncinfo(ncname(18,:),'z0','Roughness height','m','time')
        call ncinfo(ncname(19,:),'shf_bar','Sensible heat flux','K m/s','time')
        call ncinfo(ncname(20,:),'sfcbflx','Surface Buoyancy Flux','K m/s','time')
        call ncinfo(ncname(21,:),'lhf_bar','Latent heat flux','kg/kg m/s','time')
        call open_nc(fname,  ncid)
        call define_nc( ncid, NVar, ncname)
      end if
    end if

  end subroutine inittimestat

!>Run timestat. Calculate and write the statistics
  subroutine timestat

    use modglobal,  only : i1,j1,kmax,zf,dzf,cu,cv,rv,rd,&
                          rslabs,timee,dt_lim,rk3step,cexpnr,ifoutput
!
    use modfields,  only : um,vm,wm,e12m,ql0,u0av,v0av,rhof
    use modsurfdata,only : wtsurf, wqsurf, isurf,ustar,tstar,qstar,z0,oblav,qts,thls,&
                           Qnet, H, LE, G0, rs, ra, tskin, tendskin
    use modmpi,     only : my_real,mpi_sum,mpi_max,mpi_min,comm3d,mpierr,myid
    use modstat_nc,  only : lnetcdf, writestat_nc
    implicit none

    real   :: zbaseavl, ztopavl, ztopmaxl, ztop,zbaseminl
    real   :: qlintavl, qlintmaxl, tke_totl
    real   :: ccl, wmaxl, qlmaxl
    real   :: ust,tst,qst,ustl,tstl,qstl
    real   :: usttst, ustqst, usttstl, ustqstl
    real   :: wts, wqls,wtvs
    real   :: c1,c2 !Used to calculate wthvs
    real,dimension(nvar) :: vars

    ! lsm variables
    real   :: Qnetavl, Havl, LEavl, G0avl, tendskinavl, rsavl, raavl, tskinavl
    real   :: Qnetav, Hav, LEav, G0av, tendskinav, rsav, raav, tskinav
    integer:: i, j, k

    if (.not.(ltimestat)) return
    if (rk3step/=3) return
    if (timee==0) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+dtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    !      -----------------------------------------------------------
  !     1     EVALUATION OF CLOUD COVER, CLOUD BASE, ETC.
  !    -----------------------------------------------------------

  !     -----------------------------------------------------
  !     1.   Set A:  entrainment and time evolution
  !     -----------------------------------------------------

    zbaseavl = 0.0
    ztopavl = 0.0
    zbaseminl = zf(kmax)

    call calcblheight

  !     --------------------------------------------------------------
  !     9.2  liq. waterpath, cloudcover, cloudbase and cloudtop
  !     --------------------------------------------------------------

    ccl      = 0.0
    qlintavl = 0.0
    qlintmaxl= 0.0
    tke_totl = 0.0

    do j=2,j1
    do i=2,i1
      qlint     = 0.0
      do k=1,kmax
        qlint = qlint + ql0(i,j,k)*rhof(k)*dzf(k)
      end do
      if (qlint>0.) then
        ccl      = ccl      + 1.0
        qlintavl = qlintavl + qlint
        qlintmaxl = max(qlint,qlintmaxl)
      end if

      do k=1,kmax
        if (ql0(i,j,k) > 0.) then
        zbaseavl = zbaseavl + zf(k)
        zbaseminl = min(zf(k),zbaseminl)
        exit
        end if
      end do

    end do
    end do

    call MPI_ALLREDUCE(ccl   , cc   , 1,    MY_REAL, &
                          MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(qlintavl, qlintav, 1,    MY_REAL, &
                          MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(qlintmaxl, qlintmax, 1,    MY_REAL, &
                          MPI_MAX, comm3d,mpierr)
    call MPI_ALLREDUCE(zbaseavl, zbaseav, 1,    MY_REAL, &
                          MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(zbaseminl, zbasemin, 1,    MY_REAL, &
                          MPI_MIN, comm3d,mpierr)
  !     ---------------------------------------
  !     9.3  determine maximum ql_max and w_max
  !     ---------------------------------------

    wmaxl  = 0.0
    qlmaxl = 0.0
    ztopavl = 0.0
    ztopmaxl = 0.0

    do  i=2,i1
    do  j=2,j1
      ztop  = 0.0

      do  k=1,kmax
        if (ql0(i,j,k) > 0) ztop = zf(k)
        wmaxl = max(wm(i,j,k),wmaxl)
        qlmaxl = max(ql0(i,j,k),qlmaxl)
      end do

      ztopavl = ztopavl + ztop
      if (ztop > ztopmaxl) ztopmaxl = ztop
    end do
    end do

    call MPI_ALLREDUCE(wmaxl   , wmax   , 1,    MY_REAL, &
                          MPI_MAX, comm3d,mpierr)
    call MPI_ALLREDUCE(qlmaxl, qlmax, 1,    MY_REAL, &
                          MPI_MAX, comm3d,mpierr)
    call MPI_ALLREDUCE(ztopavl, ztopav, 1,    MY_REAL, &
                          MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(ztopmaxl, ztopmax, 1,    MY_REAL, &
                          MPI_MAX, comm3d,mpierr)
  !     -------------------------
  !     9.4  normalise the fields
  !     -------------------------

    if (cc > 0.0) then
      zbaseav = zbaseav/cc
      ztopav  = ztopav/cc
    else
      zbaseav = 0.0
      ztopav = 0.0
    end if

    cc      = cc/rslabs
    qlintav = qlintav / rslabs !domain averaged liquid water path

  !     -------------------------
  !     9.5  Domain Averaged TKE
  !     -------------------------

    do  i=2,i1
    do  j=2,j1
    do  k=1,kmax

      tke_totl = tke_totl + 0.5*( &
                          (0.5*(um(i,j,k)+um(i+1,j,k))+cu-u0av(k))**2 &
                          +(0.5*(vm(i,j,k)+vm(i,j+1,k))+cv-v0av(k))**2 &
                          +(0.5*(wm(i,j,k)+wm(i,j,k+1))           )**2 &
                                ) + e12m(i,j,k)**2
    end do
    end do
    end do

    call MPI_ALLREDUCE(tke_totl, tke_tot, 1,    MY_REAL, &
                          MPI_SUM, comm3d,mpierr)

    tke_tot = tke_tot/rslabs
!     -------------------------
!     9.6  Horizontally  Averaged ustar, tstar and obl
!     -------------------------
    ustl=sum(ustar(2:i1,2:j1))
    tstl=sum(tstar(2:i1,2:j1))
    qstl=sum(qstar(2:i1,2:j1))

    usttstl=sum(ustar(2:i1,2:j1)*tstar(2:i1,2:j1))
    ustqstl=sum(ustar(2:i1,2:j1)*qstar(2:i1,2:j1))

!       ustl=min(1e-20, ustl)
!       tstl=min(1e-20, tstl)
!       qstl=min(1e-20, qstl)

    call MPI_ALLREDUCE(ustl, ust, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(tstl, tst, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(qstl, qst, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)

    ust = ust / rslabs
    tst = tst / rslabs
    qst = qst / rslabs

    if(isurf < 3) then
      call MPI_ALLREDUCE(usttstl, usttst, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(ustqstl, ustqst, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)

      usttst = usttst / rslabs
      ustqst = ustqst / rslabs
    end if

    !Constants c1 and c2
    c1   = 1.+(rv/rd-1)*qts
    c2   = (rv/rd-1)

    if(isurf >= 3) then
      wts  = wtsurf
      wqls = wqsurf
      wtvs = c1*wts + c2*thls*wqls
    else
      wts  = -usttst
      wqls = -ustqst
      wtvs = c1*wts + c2*thls*wqls
    end if

  !  9.7  Create statistics for the land surface scheme
    if(isurf == 1) then
      Qnetavl      = sum(Qnet(2:i1,2:j1))
      Havl         = sum(H(2:i1,2:j1))
      LEavl        = sum(LE(2:i1,2:j1))
      G0avl        = sum(G0(2:i1,2:j1))
      tendskinavl  = sum(tendskin(2:i1,2:j1))
      rsavl        = sum(rs(2:i1,2:j1))
      raavl        = sum(ra(2:i1,2:j1))
      tskinavl     = sum(tskin(2:i1,2:j1))

      call MPI_ALLREDUCE(Qnetavl,     Qnetav,     1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(Havl,        Hav,        1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(LEavl,       LEav,       1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(G0avl,       G0av,       1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(tendskinavl, tendskinav, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(rsavl,       rsav,       1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(raavl,       raav,       1,  MY_REAL,MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(tskinavl,    tskinav,    1,  MY_REAL,MPI_SUM, comm3d,mpierr)

      Qnetav        = Qnetav      / rslabs
      Hav           = Hav         / rslabs
      LEav          = LEav        / rslabs
      G0av          = G0av        / rslabs
      tendskinav    = tendskinav  / rslabs
      rsav          = rsav        / rslabs
      raav          = raav        / rslabs
      tskinav       = tskinav     / rslabs
    end if

  !  9.8  write the results to output file
  !     ---------------------------------------

    if(myid==0)then
       !tmser1
      open (ifoutput,file='tmser1.'//cexpnr,position='append')
      write( ifoutput,'(f10.2,f6.3,4f12.3,f10.4,5f9.3)') &
          timee, &
          cc, &
          zbaseav, &
          ztopav, &
          ztopmax, &
          zi, &
          we, &
          qlintav*1000., &
          qlintmax*1000., &
          wmax, &
          tke_tot*dzf(1), &
          qlmax*1000.
      close(ifoutput)

      !tmsurf
      open (ifoutput,file='tmsurf.'//cexpnr,position='append')
      write( ifoutput,'(f10.2,3e11.3,2f11.3,4e11.3)') &
          timee   ,&
          ust     ,&
          tst     ,&
          qst     ,&
          oblav   ,&
          thls    ,&
          z0      ,&
          wts     ,&
          wtvs    ,&
          wqls
      close(ifoutput)

      if (isurf == 1) then
        !tmlsm
        open (ifoutput,file='tmlsm.'//cexpnr,position='append')
        write(ifoutput,'(f10.2,8f11.3)') &
            timee       ,&
            Qnetav      ,&
            Hav         ,&
            LEav        ,&
            G0av        ,&
            tendskinav  ,&
            rsav        ,&
            raav        ,&
            tskinav
        close(ifoutput)
      end if
      if (lnetcdf) then
        vars( 1) = timee
        vars( 2) = cc
        vars( 3) = zbaseav
        vars( 4) = ztopav
        vars( 5) = ztopmax
        vars( 6) = zi
        vars( 7) = we
        vars( 8) = qlintav
        vars( 9) = qlintmax
        vars(10) = wmax
        vars(11) = tke_tot*dzf(1)
        vars(12) = qlmax
        vars(13) = ust
        vars(14) = tst
        vars(15) = qst
        vars(16) = oblav
        vars(17) = thls
        vars(18) = z0
        vars(19) = wts
        vars(20) = wtvs
        vars(21) = wqls
        call writestat_nc(ncid,nvar,ncname,vars,nrec,.true.)
      end if
    end if

  end subroutine timestat

!>Calculate the boundary layer height
!!
!! There are 3 available ways to calculate the boundary layer height:
!! - By determining the minimum flux in some scalar, e.g. buoyancy
!! - By determining the minimum local gradient of some scalar, averaged over a definable number of columns
!! - By monitoring a threshold value of some scalar, averaged over a definable number of columns
  subroutine calcblheight
    use modglobal,  only : ih,i1,jh,j1,kmax,k1,cp,rlv,imax,rd,zh,dzh,zf,dzf,rv,rslabs,iadv_sv,iadv_kappa
    use modfields,  only : w0,qt0,qt0h,ql0,thl0,thl0h,thv0h,sv0,exnf,whls
    use modsurfdata,only :svs
    use modmpi,     only : mpierr, comm3d,mpi_sum,my_real
    implicit none
    real   :: zil, dhdt,locval,oldlocval
    integer :: location,i,j,k,nsamp,stride
    real, allocatable,dimension(:,:,:) :: blh_fld,  sv0h
    real, allocatable, dimension(:) :: profile, gradient, dgrad
    allocate(blh_fld(2-ih:i1+ih,2-jh:j1+jh,k1),sv0h(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(profile(k1),gradient(k1),dgrad(k1))
    zil = 0.0
    gradient = 0.0
    dgrad = 0.0
    select case (iblh_meth)
    case (iblh_flux)
      select case (iblh_var)
      case(iblh_qt)
        blh_fld = w0*qt0h
      case(iblh_thl)
        blh_fld = w0*thl0h
      case(iblh_thv)
        blh_fld = w0*thv0h
      case(1:)
        if (iadv_sv(iblh_var)==iadv_kappa) then
          call halflev_kappa(sv0(:,:,:,iblh_var),sv0h)
        else
          do  j=2,j1
          do  i=2,i1
          do  k=2,k1
            sv0h(i,j,k) = (sv0(i,j,k,iblh_var)*dzf(k-1)+sv0(i,j,k-1,iblh_var)*dzf(k))/(2*dzh(k))
          enddo
          enddo
          enddo
          sv0h(2:i1,2:j1,1) = svs(iblh_var)
          blh_fld = w0*sv0h
        end if
      end select
    case (iblh_grad,iblh_thres)
      select case (iblh_var)
      case(iblh_qt)
        blh_fld = qt0
      case(iblh_thl)
        blh_fld = thl0
      case(iblh_thv)
        do k=1,k1
          blh_fld(2:i1,2:j1,k) = (thl0(2:i1,2:j1,k)+rlv*ql0(2:i1,2:j1,k)/(cp*exnf(k))) &
                      *(1+(rv/rd-1)*qt0(2:i1,2:j1,k)-rv/rd*ql0(2:i1,2:j1,k))
        end do
      case(1:)
        blh_fld = sv0(2:i1,2:j1,1:k1,iblh_var)
      end select
    end select

    select case (iblh_meth)

    case (iblh_flux)
      stride = ceiling(real(imax)/real(blh_nsamp))
      do i=2,stride+1
        nsamp =  ceiling(real(i1-i+1)/real(stride))
        do j=2,j1
          zil = zil + nsamp*zh(minloc(sum(blh_fld(i:i1:stride,j,:),1),1))
        end do
      end do

    case (iblh_grad)
      stride = ceiling(real(imax)/real(blh_nsamp))
      do i=2,stride+1
        nsamp =  ceiling(real(i1-i+1)/real(stride))
        do j=2,j1
          profile  = sum(blh_fld(i:i1:stride,j,:),1)
          gradient(2:k1) = abs(profile(2:k1) - profile(1:kmax))/dzh(2:k1)
          dgrad(2:kmax)    = (gradient(3:k1) - gradient(2:kmax))/dzf(2:kmax)
          location = maxloc(gradient,1)
          zil  = zil + nsamp*(zh(location-1) - dzh(location)*dgrad(location-1)/(dgrad(location)-dgrad(location-1)))
        enddo
      enddo
    case (iblh_thres)
      stride = ceiling(real(imax)/real(blh_nsamp))
      do i=2,stride+1
        nsamp =  ceiling(real(i1-i+1)/real(stride))
        do j=2,j1
          locval = 0.0
          do k=kmax,1,-1
            oldlocval = locval
            locval = blh_sign*sum(blh_fld(i:i1:stride,j,k))/nsamp
            if (locval < blh_sign*blh_thres) then
              zil = zil + nsamp *(zf(k) +  (blh_sign*blh_thres-locval) &
                        *dzh(k+1)/(oldlocval-locval))
              exit
            endif
          enddo
        enddo
      enddo

    end select

    call MPI_ALLREDUCE(zil, zi, 1, MY_REAL, MPI_SUM, comm3d,mpierr)
    zi = zi / rslabs

    if (ziold< 0) ziold = zi
    dhdt = (zi-ziold)/dtav
    ziold = zi

    k=2
    do while (zh(k)<zi .and. k < kmax)
      k=k+1
    end do
    we = dhdt - whls (k)   !include for large-scale vertical velocity

    deallocate(blh_fld,sv0h)
    deallocate(profile,gradient,dgrad)

  end subroutine calcblheight

!> Clean up when leaving the run
  subroutine exittimestat
    use modmpi, only : myid
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(ltimestat .and. lnetcdf .and. myid==0) call exitstat_nc(ncid)
  end subroutine exittimestat

end module modtimestat
