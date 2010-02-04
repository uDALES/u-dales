!> \file modbudget.f90
!!  Calculates the turbulent budgets


!>
!! Calculates the turbulent budgets
!>
!! Profiles of the resolved and SFS TKE budgets. Written to budget.expnr and sbbudget.expnr
!! If netcdf is true, this module also writes in the profiles.expnr.nc output
!!  \author Hans Cuijpers, IMAU
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
module modbudget


  implicit none
  PRIVATE
  PUBLIC :: initbudget,budgetstat,exitbudget
  save
!NetCDF variables
  integer,parameter :: nvar = 18
  integer :: ncid,nrec = 0
  character(80) :: fname = 'budget.xxx.nc'
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav, timeav
  integer :: idtav, itimeav,tnext,tnextwrite
  integer :: nsamples
  logical :: lbudget= .false. ! switch for turbulent TKE budget

  !time averaged fields, resolved TKE
  real, allocatable :: tkemn(:)   !< Resolved TKE
  real, allocatable :: shrmn(:)   !< Shear
  real, allocatable :: buomn(:)   !< Buoyancy
  real, allocatable :: trspmn(:)  !< Transport
  real, allocatable :: ptrspmn(:) !< Pressure transport (redistribution)
  real, allocatable :: dissmn(:)  !< Dissipation
  real, allocatable :: stormn(:)  !< Storage = dE/dt
  real, allocatable :: budgmn(:)  !< Budget = sum of contributions excl storage
  real, allocatable :: residmn(:) !< Residual = budget - storage
  real, allocatable :: tkeb(:)    !< TKE at beginning of averaging periode. Used to calculate stormn
  real, allocatable :: tkeav(:)   !< Must be module global to calculate dE/dt in combination with tkeb


  !time averaged fields, subgrid TKE
  real, allocatable :: sbtkemn(:)   !< subgrid TKE
  real, allocatable :: sbshrmn(:)   !< subgrid Shear
  real, allocatable :: sbbuomn(:)   !< subgrid Buoyancy
  real, allocatable :: sbdissmn(:)  !< subgrid Dissipation
  real, allocatable :: sbstormn(:)  !< subgrid Storage = dE/dt
  real, allocatable :: sbbudgmn(:)  !< subgrid Budget = sum of contributions excl storage
  real, allocatable :: sbresidmn(:) !< subgrid Residual = budget - storage
  real, allocatable :: sbtkeb(:)    ! TKE at beginning of averaging periode. Used to calculate stormn
  real, allocatable :: sbtkeav(:)   ! Must be module global to calculate dE/dt in combination with tkeb
  real, allocatable :: ekmmn(:)     !< Turbulent exchange coefficient momentum
  real, allocatable :: khkmmn(:)    !< Kh / Km, in post-processing used to determine filter-grid ratio


  logical :: ltkeb     !Switch to tell if the tke   at beg of av periode has been stored
  logical :: lsbtkeb   !Switch to tell if the sbtke at beg of av periode has been stored

contains
!> Initialization routine, reads namelists and inits variables
  subroutine initbudget
    use modmpi,    only : myid,mpierr, comm3d,my_real, mpi_logical
    use modglobal, only : dtmax,idtmax, k1,ifnamopt,fname_options, ifoutput, cexpnr,dtav_glob,timeav_glob,ladaptive,dt_lim,btime,kmax,tres
    use modstat_nc, only : lnetcdf, open_nc,define_nc,redefine_nc,ncinfo,writestat_dims_nc
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav,ncid_prof=>ncid


    implicit none

    integer ierr
    namelist/NAMBUDGET/ &
         dtav,timeav,lbudget

    dtav=dtav_glob;timeav=timeav_glob

    if(myid==0)then
       open(ifnamopt,file=fname_options,status='old',iostat=ierr)
       read (ifnamopt,NAMBUDGET,iostat=ierr)
       if (ierr > 0) then
         print *, 'Problem in namoptions NAMBUDGET'
         print *, 'iostat error: ', ierr
         stop 'ERROR: Problem in namoptions NAMBUDGET'
       endif
       write(6 ,NAMBUDGET)
       close(ifnamopt)
    end if

    call MPI_BCAST(timeav     ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(dtav       ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(lbudget    ,1,MPI_LOGICAL,0,comm3d,mpierr)
    idtav = dtav/tres
    dtav  = idtav*tres
    itimeav = timeav/tres
    timeav  = itimeav*tres

    tnext      = dtav   +btime
    tnextwrite = timeav +btime
    nsamples = itimeav/idtav
    if(.not.(lbudget)) return
    dt_lim = min(dt_lim,tnext)

    if (abs(timeav/dtav-nsamples)>1e-4) then
      stop 'timeav must be a integer multiple of dtav'
    end if
    if (.not. ladaptive .and.abs( dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if

    !time averaged fields, resolved TKE
    allocate(tkemn(k1),tkeb(k1),tkeav(k1),shrmn(k1),buomn(k1),trspmn(k1))
    allocate(ptrspmn(k1),stormn(k1),budgmn(k1),residmn(k1),dissmn(k1))
    !horizontal mean fields, resolved TKE
    !time averaged fields, subgrid TKE
    allocate(sbtkemn(k1),sbtkeb(k1),sbtkeav(k1),sbshrmn(k1),sbbuomn(k1))
    allocate(sbstormn(k1),sbbudgmn(k1),sbresidmn(k1),&
         sbdissmn(k1),ekmmn(k1),khkmmn(k1))


    !Setting time mean variables to zero
    tkemn=0.;shrmn=0.;buomn=0.;trspmn=0.;dissmn=0.
    ptrspmn=0.;stormn=0.;budgmn=0.;residmn=0.
    tkeb=0. !tke at start of averaging period
    sbtkemn=0.;sbshrmn=0.;sbbuomn=0.;sbstormn=0.;
    sbbudgmn=0.;sbresidmn=0.;sbdissmn=0.;
       ekmmn=0.;khkmmn=0.
    sbtkeb=0.;
    ltkeb=.false. ; lsbtkeb=.false.

   !Preparing output files
    if(myid==0)then
       open (ifoutput,file='budget.'//cexpnr,status='replace')
       close (ifoutput)
       open (ifoutput,file='sbbudget.'//cexpnr,status='replace')
       close (ifoutput)
    endif
    if (lnetcdf) then
      idtav = idtav_prof
      itimeav = itimeav_prof
      tnext      = idtav+btime
      tnextwrite = itimeav+btime
      nsamples = itimeav/idtav
     if (myid==0) then
        fname(8:10) = cexpnr
        call ncinfo(tncname(1,:),'time','Time','s','time')
        call ncinfo(ncname( 1,:),'tker','Resolved TKE','m/s^2','tt')
        call ncinfo(ncname( 2,:),'shr','Resolved Shear','m/s^2','tt')
        call ncinfo(ncname( 3,:),'buo','Resolved Buoyancy','m/s^2','tt')
        call ncinfo(ncname( 4,:),'trsp','Resolved Transport','m/s^2','tt')
        call ncinfo(ncname( 5,:),'ptrsp','Resolved Pressure transport (redistribution)','m/s^2','tt')
        call ncinfo(ncname( 6,:),'diss','Resolved Dissipation','m/s^2','tt')
        call ncinfo(ncname( 7,:),'budg','Resolved Storage = dE/dt','m/s^2','tt')
        call ncinfo(ncname( 8,:),'stor','Resolved Budget = sum of contributions excl storage','m/s^2','tt')
        call ncinfo(ncname( 9,:),'resid','Resolved Residual = budget - storage','m/s^2','tt')
        call ncinfo(ncname(10,:),'sbtke','Subgrid TKE','m/s^2','tt')
        call ncinfo(ncname(11,:),'sbshr','Subgrid Shear','m/s^2','tt')
        call ncinfo(ncname(12,:),'sbbuo','Subgrid Buoyancy','m/s^2','tt')
        call ncinfo(ncname(13,:),'sbdiss','Subgrid Dissipation','m/s^2','tt')
        call ncinfo(ncname(14,:),'sbstor','Subgrid Storage = dE/dt','m/s^2','tt')
        call ncinfo(ncname(15,:),'sbbudg','Subgrid Budget = sum of contributions excl storage','m/s^2','tt')
        call ncinfo(ncname(16,:),'sbresid','Subgrid Residual = budget - storage','m/s^2','tt')
        call ncinfo(ncname(17,:),'ekm','Turbulent exchange coefficient momentum','m/s^2','tt')
        call ncinfo(ncname(18,:),'khkm   ','Kh / Km, in post-processing used to determine filter-grid ratio','m/s^2','tt')
        call open_nc(fname,ncid,n3=kmax)
        call define_nc( ncid, 1, tncname)
        call writestat_dims_nc(ncid)
        call redefine_nc(ncid_prof)
        call define_nc( ncid_prof, NVar, ncname)
     end if

   end if

  end subroutine initbudget


!> General routine, does the timekeeping
  subroutine budgetstat

    use modglobal, only : rk3step,timee, dt_lim
    implicit none

    if (.not. lbudget) return
    if (rk3step/=3) return
    if(timee<tnext .and. timee<tnextwrite) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if
    if (timee>=tnext) then
      tnext = tnext+idtav
      call do_genbudget
      call do_gensbbudget
    end if
    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+itimeav
      call writebudget
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
  end subroutine budgetstat

!> Performs the resolved budget calculations
  subroutine do_genbudget
    use modglobal,  only : i1,j1,k1,kmax,dzf,dzh, &
                          rslabs,cu,cv,iadv_thl,grav, &
                          dxi,dyi,dx2i,dy2i
    use modsurfdata,only : thvs, ustar
    use modsubgrid, only : ekm
    use modpois,    only : p
    use modfields,  only : u0,v0,w0,thl0h,thv0h,u0av,v0av
!cstep    use modtilt,    only : adjustbudget,ltilted
    use modmpi,     only : nprocs,comm3d,nprocs,my_real, mpi_sum,mpierr

    implicit none
    integer :: i,j,k,km,kp,jm,jp
    real    :: uwrs,vwrs,egp,ph
    real    :: emom,emmo,emop,epom,eopm,epmo,eomp,eomm,empo
    real    :: fu,fv

    !horizontal mean fields, resolved TKE
    real, allocatable,dimension(:) :: &
          shrav, buoav,trspav, ptrspav,dissav
    real,allocatable,dimension(:) :: &
         tkeavl,shravl,buoavl,trspavl,ptrspavl
    real, allocatable, dimension(:) :: &
         buoz,weres,weresl
    real,allocatable, dimension(:) :: &
         tau1m,tau2m,tau3m,tau1ml,tau2ml,tau3ml
    real,allocatable, dimension(:,:,:) :: &
         tau1j,tau2j,tau3j
    real,allocatable, dimension(:) :: &
         subx,suby,subz,subxl,subyl,subzl



    !help variables

    !**************************************************
    ! 1.1  ALLOCATE VARIABLES
    !**************************************************
    allocate(shrav(k1),buoav(k1),trspav(k1))
    allocate(ptrspav(k1),dissav(k1))
    allocate(tkeavl(k1),shravl(k1),buoavl(k1),trspavl(k1),ptrspavl(k1))

    !help variables
    allocate(buoz(k1),weres(k1),weresl(k1))
    allocate(tau1m(k1),tau2m(k1),tau3m(k1),tau1ml(k1),tau2ml(k1),tau3ml(k1))
    allocate(tau1j(i1+1,j1+1,k1),tau2j(i1+1,j1+1,k1),tau3j(i1+1,j1+1,k1))
    allocate(subx(k1),suby(k1),subz(k1),subxl(k1),subyl(k1),subzl(k1))

    !**************************************************
    ! 2.1 RESET ARRAYS FOR SLAB AVERAGES
    !**************************************************

    tkeav=0.;shrav=0.;buoav=0.;trspav=0.;ptrspav=0.;dissav=0.;
    tkeavl=0.;shravl=0.;buoavl=0.;trspavl=0.; ptrspavl=0.



    !**************************************
    ! 3: Calculation of individual terms
    !**************************************

    !-------------------------------------------------
    ! 3.1  Buoyancy
    !-------------------------------------------------

    do k=2,k1
       km = k-1
       buoz(k) = 0.0
       do j=2,j1
          do i=2,i1
             buoz(k) = buoz(k) + w0(i,j,k)*thv0h(i,j,k)
             !buoyancy in x-dir calculated in modtilt.
          end do
       end do
       buoz(k) = grav/thvs * buoz(k)/rslabs
       buoavl(k) = buoz(k)
    end do
    !No buoyancy at the surface
    buoavl(1) = 0.
    !If tilted surface, add buoyancy in x-dir.
!    if(ltilted) call adjustbudget(buoavl)
    call MPI_ALLREDUCE(buoavl, buoav, k1,    MY_REAL, &
         MPI_SUM, comm3d,mpierr)


    !-------------------------------------------------
    ! 3.2  Shear
    !-------------------------------------------------
    do k=2,kmax
       uwrs = 0.0
       vwrs = 0.0
       do j=2,j1
       do i=2,i1
          uwrs = uwrs+(w0(i,j,k)+w0(i-1,j,k))*(u0(i,j,k-1)+u0(i,j,k))/4.
          vwrs = vwrs+(w0(i,j,k)+w0(i,j-1,k))*(v0(i,j,k-1)+v0(i,j,k))/4.
       end do
       end do
       uwrs = uwrs/rslabs
       vwrs = vwrs/rslabs
       shravl(k) = -uwrs*(u0av(k)-u0av(k-1))/dzh(k) &
                   -vwrs*(v0av(k)-v0av(k-1))/dzh(k)
    end do
    !No shear at the surface
    shravl(1)=0.
    shravl(k1)=shravl(kmax)
    call MPI_ALLREDUCE(shravl, shrav, k1,    MY_REAL, &
         MPI_SUM, comm3d,mpierr)

    !-------------------------------------------------
    ! 3.3  Transport, resolved TKE and storage
    !-------------------------------------------------

  do k=1,kmax
    tkeavl(k)   = 0.0
    weresl(k) = 0.0
    do j=2,j1
    do i=2,i1
      egp = 0.5*( (0.5*(u0(i,j,k)+u0(i+1,j,k))-(u0av(k)-cu))**2 &
                 +(0.5*(v0(i,j,k)+v0(i,j+1,k))-(v0av(k)-cv))**2 &
                 +(0.5*(w0(i,j,k)+w0(i,j,k+1))             )**2 )

      tkeavl(k) = tkeavl(k) + egp
      weresl(k) = weresl(k) + egp * 0.5*(w0(i,j,k)+w0(i,j,k+1))
    end do
    end do
    tkeavl(k)   = tkeavl(k)/rslabs
    weresl(k) = weresl(k)/rslabs
  end do

  call MPI_ALLREDUCE(tkeavl, tkeav, k1,    MY_REAL, &
         MPI_SUM, comm3d,mpierr)
  call MPI_ALLREDUCE(weresl, weres, k1,    MY_REAL, &
         MPI_SUM, comm3d,mpierr)

  do k=2,kmax
    trspav(k) = -(weres(k)-weres(k-1))/dzh(k)
  end do

  trspav(k1) = 0.
  trspav(1)  = -weres(1)/(0.5*dzh(1))

  !storage: tke at beginning of averaging period
    if (.not.(ltkeb)) then
       tkeb = tkeav
       ltkeb = .true.
    endif


  !-------------------------------------------------
  ! 3.4  Pressure redistribution (pressure transport)
  !-------------------------------------------------
  do k=2,kmax
     do j=2,j1
     do i=2,i1
        ph = (p(i,j,k)*dzf(k-1)+p(i,j,k-1)*dzf(k))/(2*dzh(k))
        ptrspavl(k) = ptrspavl(k) + w0(i,j,k)*ph
     end do
     end do
     ptrspavl(k) = ptrspavl(k)/rslabs
  enddo

  ptrspavl(1)  = 0.
  call MPI_ALLREDUCE(ptrspavl, ptrspav, k1,    MY_REAL, &
         MPI_SUM, comm3d,mpierr)

  ptrspav(k1) = ptrspav(kmax)
  ptrspav(1)  = 0.
  do k=1,kmax
    ptrspav(k) = -(ptrspav(k+1)-ptrspav(k))/dzf(k)
  end do

  !-------------------------------------------------
  ! 3.5  Dissipation
  !-------------------------------------------------

  do k=2,kmax
    kp = k+1
    km = k-1
    tau1ml(k) = 0.0
    do j=2,j1
      jp = j+1
      jm = j-1
    do i=2,i1
      emom = ( dzf(km) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                   dzf(k)  * ( ekm(i,j,km) + ekm(i-1,j,km) ) ) / &
                 ( 4.   * dzh(k) )

      emop = ( dzf(kp) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                   dzf(k)  * ( ekm(i,j,kp) + ekm(i-1,j,kp) ) ) / &
                 ( 4.   * dzh(kp) )

      empo = 0.25 * ( &
                 ekm(i,j,k)+ekm(i,jp,k)+ekm(i-1,jp,k)+ekm(i-1,j,k) )

      emmo = 0.25 * ( &
                 ekm(i,j,k)+ekm(i,jm,k)+ekm(i-1,jm,k)+ekm(i-1,j,k) )

      tau1j(i,j,k) = ( &
               ( ekm(i,j,k)  * (u0(i+1,j,k)-u0(i,j,k)) &
                -ekm(i-1,j,k)* (u0(i,j,k)-u0(i-1,j,k)) ) * 2. * dx2i &
               + &
               ( empo * ( (u0(i,jp,k)-u0(i,j,k))   *dyi &
                         +(v0(i,jp,k)-v0(i-1,jp,k))*dxi) &
                -emmo * ( (u0(i,j,k)-u0(i,jm,k))   *dyi &
                         +(v0(i,j,k)-v0(i-1,j,k))  *dxi) ) * dyi &
               + &
               ( emop * ( (u0(i,j,kp)-u0(i,j,k))   /dzh(kp) &
                         +(w0(i,j,kp)-w0(i-1,j,kp))*dxi) &
                -emom * ( (u0(i,j,k)-u0(i,j,km))   /dzh(k) &
                         +(w0(i,j,k)-w0(i-1,j,k))  *dxi)   ) / dzf(k) &
               )

      tau1ml(k) = tau1ml(k) + tau1j(i,j,k)

    enddo
    enddo
    tau1ml(k) = tau1ml(k)/rslabs
  enddo

  do k=2,kmax
    kp = k+1
    km = k-1
    tau2ml(k) = 0.0
    do j=2,j1
      jp = j+1
      jm = j-1
    do i=2,i1


      eomm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
          dzf(k)  * ( ekm(i,j,km) + ekm(i,jm,km) ) ) / &
          ( 4.   * dzh(k) )

      eomp = ( dzf(kp) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
        dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) / &
        ( 4.   * dzh(kp) )

      emmo = 0.25  * ( &
        ekm(i,j,k)+ekm(i,jm,k)+ekm(i-1,jm,k)+ekm(i-1,j,k)  )

      epmo = 0.25  * ( &
        ekm(i,j,k)+ekm(i,jm,k)+ekm(i+1,jm,k)+ekm(i+1,j,k)  )


      tau2j(i,j,k) = ( &
               ( epmo * ( (v0(i+1,j,k)-v0(i,j,k))   *dxi &
                         +(u0(i+1,j,k)-u0(i+1,jm,k))*dyi) &
                -emmo * ( (v0(i,j,k)-v0(i-1,j,k))   *dxi &
                         +(u0(i,j,k)-u0(i,jm,k))    *dyi)  ) * dxi &
               + &
               ( ekm(i,j,k) * (v0(i,jp,k)-v0(i,j,k)) &
                -ekm(i,jm,k)* (v0(i,j,k)-v0(i,jm,k))  ) * 2. * dy2i &
               + &
               ( eomp * ( (v0(i,j,kp)-v0(i,j,k))    /dzh(kp) &
                         +(w0(i,j,kp)-w0(i,jm,kp))  *dyi) &
                -eomm * ( (v0(i,j,k)-v0(i,j,km))    /dzh(k) &
                         +(w0(i,j,k)-w0(i,jm,k))    *dyi)  ) / dzf(k) &
               )

      tau2ml(k) = tau2ml(k) + tau2j(i,j,k)
    enddo
    enddo
    tau2ml(k) = tau2ml(k)/rslabs
  enddo

  do k=2,kmax
    kp = k+1
    km = k-1
    tau3ml(k) = 0.0
    do j=2,j1
      jp = j+1
      jm = j-1
    do i=2,i1

      emom = ( dzf(km) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
          dzf(k)  * ( ekm(i,j,km) + ekm(i-1,j,km) ) ) / &
          ( 4.   * dzh(k) )

      eomm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
        dzf(k)  * ( ekm(i,j,km) + ekm(i,jm,km) ) ) / &
        ( 4.   * dzh(k) )

      eopm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jp,k)  )  + &
        dzf(k)  * ( ekm(i,j,km) + ekm(i,jp,km) ) ) / &
        ( 4.   * dzh(k) )

      epom = ( dzf(km) * ( ekm(i,j,k)  + ekm(i+1,j,k)  )  + &
        dzf(k)  * ( ekm(i,j,km) + ekm(i+1,j,km) ) ) / &
        ( 4.   * dzh(k) )


      tau3j(i,j,k) = ( &
                 ( epom * ( (w0(i+1,j,k)-w0(i,j,k))    *dxi &
                           +(u0(i+1,j,k)-u0(i+1,j,km)) /dzh(k) ) &
                  -emom * ( (w0(i,j,k)-w0(i-1,j,k))    *dxi &
                           +(u0(i,j,k)-u0(i,j,km))     /dzh(k) ))*dxi &
                 + &
                 ( eopm * ( (w0(i,jp,k)-w0(i,j,k))     *dyi &
                           +(v0(i,jp,k)-v0(i,jp,km))   /dzh(k) ) &
                  -eomm * ( (w0(i,j,k)-w0(i,jm,k))     *dyi &
                           +(v0(i,j,k)-v0(i,j,km))     /dzh(k) ))*dyi &
                 + &
                  (ekm(i,j,k) * (w0(i,j,kp)-w0(i,j,k)) /dzf(k) &
                  -ekm(i,j,km)* (w0(i,j,k)-w0(i,j,km)) /dzf(km) ) * 2. &
                                                       / dzh(k) &
                 )

      tau3ml(k) = tau3ml(k) + tau3j(i,j,k)

    enddo
    enddo
    tau3ml(k) = tau3ml(k)/rslabs
  enddo


!     --------------------------------------------
!     special treatment for lowest full level: k=1
!     --------------------------------------------

   tau1ml(1) = 0.0
   tau2ml(1) = 0.0
   tau3ml(1) = 0.0

!     ----------------------------------------------start j-loop
  do j=2,j1
    jp = j+1
    jm = j-1
!     ----------------------------------------------start i-loop
  do i=2,i1

    empo = 0.25 * ( &
              ekm(i,j,1)+ekm(i,jp,1)+ekm(i-1,jp,1)+ekm(i-1,j,1)  )

    emmo = 0.25 * ( &
              ekm(i,j,1)+ekm(i,jm,1)+ekm(i-1,jm,1)+ekm(i-1,j,1)  )

    emop = ( dzf(2) * ( ekm(i,j,1) + ekm(i-1,j,1) )  + &
                 dzf(1) * ( ekm(i,j,2) + ekm(i-1,j,2) ) ) / &
               ( 4.   * dzh(2) )

    fu = ( 0.5*( ustar(i,j)+ustar(i-1,j) ) )**2  * &
                   (u0(i,j,1)+cu)/sqrt((u0(i,j,1)+cu)**2  + &
                   ((v0(i,j,1)+v0(i-1,j,1)+v0(i,jp,1)+v0(i-1,jp,1))/4. &
                   +cv)**2)


    tau1j(i,j,1) = &
             ( ekm(i,j,1)  * (u0(i+1,j,1)-u0(i,j,1)) &
              -ekm(i-1,j,1)* (u0(i,j,1)-u0(i-1,j,1)) ) * 2. * dx2i &
             + &
             ( empo * ( (u0(i,jp,1)-u0(i,j,1))   *dyi &
                       +(v0(i,jp,1)-v0(i-1,jp,1))*dxi) &
             - emmo * ( (u0(i,j,1)-u0(i,jm,1))   *dyi &
                       +(v0(i,j,1)-v0(i-1,j,1))  *dxi)   ) * dyi &
             + &
             ( emop * ( (u0(i,j,2)-u0(i,j,1))    /dzh(2) &
                       +(w0(i,j,2)-w0(i-1,j,2))  *dxi) &
             - fu   ) / dzf(1)


    tau1ml(1) = tau1ml(1) + tau1j(i,j,1)

    epmo = 0.25  * ( &
              ekm(i,j,1)+ekm(i,jm,1)+ekm(i+1,jm,1)+ekm(i+1,j,1)  )

    eomp = ( dzf(2) * ( ekm(i,j,1) + ekm(i,jm,1)  )  + &
                 dzf(1) * ( ekm(i,j,2) + ekm(i,jm,2) ) ) / &
               ( 4.   * dzh(2) )


    fv = ( 0.5*( ustar(i,j)+ustar(i,j-1) ) )**2  * &
              (v0(i,j,1)+cv)/sqrt((v0(i,j,1)+cv)**2  + &
             ((u0(i,j,1)+u0(i+1,j,1)+u0(i,jm,1)+u0(i+1,jm,1))/4.+cu)**2)


    tau2j(i,j,1) = &
            ( epmo * ( (v0(i+1,j,1)-v0(i,j,1))   *dxi &
                      +(u0(i+1,j,1)-u0(i+1,jm,1))*dyi) &
             -emmo * ( (v0(i,j,1)-v0(i-1,j,1))   *dxi &
                      +(u0(i,j,1)-u0(i,jm,1))    *dyi)   ) * dxi &
            + &
            ( ekm(i,j,1) * (v0(i,jp,1)-v0(i,j,1)) &
             -ekm(i,jm,1)* (v0(i,j,1)-v0(i,jm,1))  ) * 2. * dy2i &
            + &
            ( eomp * ( (v0(i,j,2)-v0(i,j,1))     /dzh(2) &
                      +(w0(i,j,2)-w0(i,jm,2))    *dyi) &
           - fv   ) / dzf(1)


    tau2ml(1) = tau2ml(1) + tau2j(i,j,1)

    tau3j(i,j,1) = 0.0

  end do
  end do

  tau1ml(1) = tau1ml(1)/rslabs
  tau2ml(1) = tau2ml(1)/rslabs

  call MPI_ALLREDUCE(tau1ml, tau1m, k1,    MY_REAL, &
         MPI_SUM, comm3d,mpierr)
  call MPI_ALLREDUCE(tau2ml, tau2m, k1,    MY_REAL, &
         MPI_SUM, comm3d,mpierr)
  call MPI_ALLREDUCE(tau3ml, tau3m, k1,    MY_REAL, &
         MPI_SUM, comm3d,mpierr)

  do k=1,kmax
    subxl(k) = 0.0
    subyl(k) = 0.0
    subzl(k) = 0.0
    do j=2,j1
    do i=2,i1
      subxl(k) = subxl(k) + &
                    (u0(i,j,k)-u0av(k)) * (tau1j(i,j,k)-tau1m(k))
      subyl(k) = subyl(k) + &
                    (v0(i,j,k)-v0av(k)) * (tau2j(i,j,k)-tau2m(k))
      subzl(k) = subzl(k) + &
                    (w0(i,j,k)        ) * (tau3j(i,j,k)-tau3m(k))
    end do
    end do
    subxl(k) = subxl(k)/rslabs
    subyl(k) = subyl(k)/rslabs
    subzl(k) = subzl(k)/rslabs
  end do
  call MPI_ALLREDUCE(subxl, subx, k1,    MY_REAL, &
         MPI_SUM, comm3d,mpierr)
  call MPI_ALLREDUCE(subyl, suby, k1,    MY_REAL, &
         MPI_SUM, comm3d,mpierr)
  call MPI_ALLREDUCE(subzl, subz, k1,    MY_REAL, &
         MPI_SUM, comm3d,mpierr)

  subz(k1) = subz(kmax)

  do  k=1,kmax
    dissav(k) = subx(k)+suby(k)+0.5*(subz(k+1)+subz(k))
  end do


  !**************************************************
  ! 4.1 Add slab averages to time mean
  !**************************************************
  do k=1,kmax
     tkemn(k)   = tkemn(k)   + tkeav(k)
     shrmn(k)   = shrmn(k)   + 0.5*(shrav(k+1)+shrav(k))
     buomn(k)   = buomn(k)   + 0.5*(buoav(k+1)+buoav(k))
     trspmn(k)  = trspmn(k)  + 0.5*(trspav(k+1)+trspav(k))
     ptrspmn(k) = ptrspmn(k) + ptrspav(k)
     dissmn(k)  = dissmn(k)  + dissav(k)
  enddo

  !**************************************************
  ! 5.1 Deallocate temporary arrays
  !**************************************************

  deallocate(shrav,buoav,trspav,ptrspav,dissav)
  deallocate(tkeavl,shravl,buoavl,trspavl,ptrspavl)
  deallocate(buoz,weres,weresl)
  deallocate(tau1m,tau2m,tau3m,tau1ml,tau2ml,tau3ml)
  deallocate(tau1j,tau2j,tau3j)
  deallocate(subx,suby,subz,subxl,subyl,subzl)
end subroutine do_genbudget

!> Performs the SFS - budget calculations
  subroutine do_gensbbudget
    use modglobal,  only : i1,j1,ih,jh,k1,kmax,rslabs
    use modsubgrid, only : ekm,ekh,sbdiss,sbshr,sbbuo
    use modfields,  only : e120
    use modmpi,     only : slabsum,nprocs,comm3d,nprocs,my_real, mpi_sum,mpierr
    !----------------------------
    ! 1.1 Declare allocatable
    !----------------------------
    integer i,j,k
    real, allocatable :: sbshrav(:),sbbuoav(:) ,sbdissav(:),&
         ekmav(:),khkmav(:)
    real, allocatable :: sbtkeavl(:),khkmavl(:)

    !----------------------------
    ! 1.2 Allocate variables
    !----------------------------
    allocate(sbshrav(k1),sbbuoav(k1),sbdissav(k1),ekmav(k1),khkmav(k1))
    allocate(sbtkeavl(k1),khkmavl(k1))

    !----------------------------
    ! 2. Reset variables
    !----------------------------
    sbtkeav=0.;sbshrav=0.;sbbuoav=0.;sbdissav=0;ekmav=0.;khkmav=0.
    sbtkeavl=0.;khkmavl=0.

    !----------------------------
    ! 3. Calculate sbtke and kh/km
    !----------------------------
    do k=1,k1
       do j=2,j1
       do i=2,i1
          sbtkeavl(k) = sbtkeavl(k) + e120(i,j,k)*e120(i,j,k)
          khkmavl(k)  = khkmavl(k)  + ekh(i,j,k)/ekm(i,j,k)
       enddo
       enddo
    enddo

    !----------------------------
    ! 4. Calculate slab-averaged values
    !----------------------------
    call slabsum(sbshrav ,1,k1,sbshr ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(sbbuoav ,1,k1,sbbuo ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(sbdissav,1,k1,sbdiss,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(ekmav   ,1,k1,ekm   ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call MPI_ALLREDUCE(khkmavl, khkmav, k1, MY_REAL, MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(sbtkeavl,sbtkeav,k1, MY_REAL, MPI_SUM, comm3d,mpierr)
    sbshrav  = sbshrav  / rslabs
    sbbuoav  = sbbuoav / rslabs
    sbdissav = sbdissav/ rslabs
    sbtkeav  = sbtkeav  / rslabs
    ekmav    = ekmav    / rslabs
    khkmav   = khkmav   / rslabs
    !storage: tke at beginning of averaging period
    if (.not.(lsbtkeb)) then
       sbtkeb = sbtkeav
       lsbtkeb=.true.
    endif

    !----------------------------
    ! 5. Add to time mean
    !----------------------------

  do k=1,k1
     sbshrmn(k)   = sbshrmn(k)   + sbshrav(k)
     sbbuomn(k)   = sbbuomn(k)   + sbbuoav(k)
     sbdissmn(k)  = sbdissmn(k)  + sbdissav(k)
     sbtkemn(k)   = sbtkemn(k)   + sbtkeav(k)
     ekmmn(k)     = ekmmn(k)     + ekmav(k)
     khkmmn(k)    = khkmmn(k)    + khkmav(k)
  enddo

    !Deallocate
    deallocate(sbshrav,sbbuoav,sbdissav,ekmav,khkmav)
    deallocate(sbtkeavl,khkmavl)
  end subroutine do_gensbbudget

!> Write the budgets to file
  subroutine writebudget
    use modglobal, only : kmax,k1,zf,rtimee,cexpnr,ifoutput
    use modmpi,    only : myid
    use modstat_nc,only : writestat_nc,lnetcdf
      use modgenstat, only: ncid_prof=>ncid,nrec_prof=>nrec
    implicit none
    real,dimension(k1,nvar) :: vars
    integer nsecs, nhrs, nminut,k
    nsecs   = nint(rtimee)
    nhrs    = int(nsecs/3600)
    nminut  = int(nsecs/60)-nhrs*60
    nsecs   = mod(nsecs,60)

    !Resolved TKE
    tkemn  = tkemn    /nsamples
    shrmn   = shrmn   /nsamples
    buomn   = buomn   /nsamples
    trspmn  = trspmn  /nsamples
    ptrspmn = ptrspmn /nsamples
    dissmn  = dissmn  /nsamples
    stormn  = (tkeav-tkeb)/timeav
    budgmn  = shrmn+buomn+trspmn+ptrspmn+dissmn
    residmn = budgmn - stormn

    !Subgrid TKE
    sbtkemn   = sbtkemn   /nsamples
    sbshrmn   = sbshrmn   /nsamples
    sbbuomn   = sbbuomn   /nsamples
    sbdissmn  = sbdissmn  /nsamples
    ekmmn     = ekmmn     /nsamples
    khkmmn    = khkmmn    /nsamples
    sbstormn  = (sbtkeav-sbtkeb)/timeav
    sbbudgmn  = sbshrmn+sbbuomn+sbdissmn
    sbresidmn = sbbudgmn - sbstormn


    if(myid==0) then
       open (ifoutput,file='budget.'//cexpnr,position='append')

       write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
           '#-------------------------------------------------------------------' &
            ,'#',(timeav),'--- AVERAGING TIMESTEP --- ' &
            ,nhrs,':',nminut,':',nsecs &
            ,'   HRS:MIN:SEC AFTER INITIALIZATION '
       write (ifoutput,'(A/2A/3A)') &
            '#-------------------------------------------------------------------' &
          ,'#LEV HEIGHT   |   TKE        SHEAR      BUOYANCY   ' &
          ,'  TRANSP     PRES_TRSP     DISS      BUDGET      STORAGE      RESID'&
          ,'#     (M)    |   ' &
          ,'(---------------------------------- (M/S)^2  ------------------' &
          ,'-------)'


       write(ifoutput,'(I3,F9.3,9E12.4)') &
            (k, &
            zf      (k), &
            tkemn  (k), & !!!
            shrmn   (k), &
            buomn   (k), &
            trspmn  (k), &
            ptrspmn (k), &
            dissmn  (k), &!!!
            budgmn  (k), &!!!
            stormn  (k), &
            residmn (k), &!!!
            k=1,kmax)
       close(ifoutput)

       open (ifoutput,file='sbbudget.'//cexpnr,position='append')

       write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
            '#---------------------------------------------------------------' &
            ,'#',(timeav),'--- AVERAGING TIMESTEP --- ' &
            ,nhrs,':',nminut,':',nsecs &
            ,'   HRS:MIN:SEC AFTER INITIALIZATION '
       write (ifoutput,'(A/2A/3A)') &
            '#---------------------------------------------------------------' &
          ,'#LEV HEIGHT   |   SBTKE     SBSHEAR     BUOYANCY     SBDISS' &
          ,'     SBSTORAGE    SBBUDGET   SBRESID    EKM          KH/KM '&
          ,'#       (M)   | ' &
          ,'(-------------------------------------------------- (M/S)^2 -------------' &
          ,'-----------------------------------)'


       write(ifoutput,'(I3,F9.3,9E12.4)') &
            (k, &
            zf      (k), &
            sbtkemn  (k), &
            sbshrmn   (k), &
            sbbuomn   (k), &
            sbdissmn  (k), &
            sbstormn  (k), &
            sbbudgmn  (k), &
            sbresidmn (k), &
            ekmmn     (k), & !!!
            khkmmn    (k), &
            k=1,kmax)
       close(ifoutput)
    endif !endif myid==0
      if (lnetcdf) then
        vars(:, 1) =tkemn
        vars(:, 2) =shrmn
        vars(:, 3) =buomn
        vars(:, 4) =trspmn
        vars(:, 5) =ptrspmn
        vars(:, 6) =dissmn
        vars(:, 7) =budgmn
        vars(:, 8) =stormn
        vars(:, 9) =residmn
        vars(:,10) =sbtkemn
        vars(:,11) =sbshrmn
        vars(:,12) =sbbuomn
        vars(:,13) =sbdissmn
        vars(:,14) =sbstormn
        vars(:,15) =sbbudgmn
        vars(:,16) =sbresidmn
        vars(:,17) =ekmmn
        vars(:,18) =khkmmn
        call writestat_nc(ncid_prof,nvar,ncname,vars(1:kmax,:),nrec_prof,kmax)
      end if
      !Reset time mean variables; resolved TKE
      tkemn=0.;tkeb=0.;shrmn=0.;buomn=0.;trspmn=0.;ptrspmn=0.;
      dissmn=0.;stormn=0.;budgmn=0.;residmn=0.
      ltkeb=.false.
      !Reset time mean variables; subgrid TKE
      sbtkemn=0.;sbtkeb=0.;sbshrmn=0.;sbbuomn=0.;sbdissmn=0.;sbtkeb=0.
      sbstormn=0.;sbbudgmn=0.;sbresidmn=0.;ekmmn=0.;khkmmn=0.;
      lsbtkeb=.false.
  end subroutine writebudget


!> Cleans up after the run
  subroutine exitbudget
  implicit none

    if(.not.(lbudget)) return

    deallocate(tkemn,tkeb,tkeav,shrmn,buomn,trspmn,ptrspmn,stormn,budgmn,residmn,dissmn)
    deallocate(sbtkemn,sbshrmn,sbbuomn,sbstormn,sbbudgmn,sbresidmn,sbdissmn,sbtkeb,sbtkeav)
    deallocate(ekmmn,khkmmn)
  end subroutine exitbudget

end module modbudget

