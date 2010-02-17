!> \file modstress.f90
!!  Calculates the momentum stress budget


!>
!! Calculates the mometum stress budgets
!>
!! Profiles of the resolved and SFS momentum budgets. Written to stressbudget.expnr
!! If netcdf is true, this module also writes in the profiles.expnr.nc output
!!  \author Arnold Moene, WU and David Pino UPC
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

module modstress
  use modglobal, only : longint

  implicit none
  PRIVATE
  PUBLIC :: initstressbudget,stressbudgetstat,exitstressbudget
  save

!NetCDF variables
  integer,parameter :: nvar = 11
  character(80),dimension(nvar,4) :: ncname

  real    :: dtav, timeav
  integer(kind=longint) :: idtav, itimeav,tnext,tnextwrite
  integer :: nsamples
  logical :: lstress= .false. ! switch for turbulent stress budget
  logical :: lstressb !Switch tot tell  if the stress  at beg of av periode has been stored
 
  !time averaged fields

  real, allocatable, dimension (:,:,:) :: str_res  !Residual stress
  real, allocatable, dimension (:,:,:) :: str_stor !Storage stress
  real, allocatable, dimension (:,:,:) :: str_budg !Sum of budget

  real, allocatable, dimension (:,:,:) :: tstr   !stress
  real, allocatable, dimension (:,:,:) :: tshr   !shear
  real, allocatable, dimension (:,:,:) :: tttr   !turbulent transport
  real, allocatable, dimension (:,:,:) :: tbuo   !buoyancy production
  real, allocatable, dimension (:,:,:) :: tcor   !Coriolis production
  real, allocatable, dimension (:,:,:) :: tptr   !pressure-velocity interaction
  real, allocatable, dimension (:,:,:) :: tdis   !disspation
  real, allocatable, dimension (:,:,:) :: tadv   !advection

  real, allocatable, dimension (:,:,:) :: rstrb, rstre ! local storage for the stress at the beginning and at the end of averaging period


contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initstressbudget
    use modmpi,    only : myid, mpierr, comm3d, my_real, mpi_logical
    use modglobal, only : dtmax,idtmax, i2, j2, k1, ifnamopt,fname_options, ifoutput, cexpnr,dtav_glob,timeav_glob,ladaptive,dt_lim,btime,tres
    use modstat_nc, only : lnetcdf, redefine_nc,define_nc,ncinfo
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav,ncid_prof=>ncid

    implicit none

    integer ierr
    namelist/NAMSTRESS/ &
         dtav,timeav,lstress

    dtav=dtav_glob;timeav=timeav_glob

    if(myid==0)then
       open(ifnamopt,file=fname_options,status='old',iostat=ierr)
       read (ifnamopt,NAMSTRESS,iostat=ierr)
       write(6 ,NAMSTRESS)
       close(ifnamopt)
    end if

    call MPI_BCAST(timeav     ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(dtav       ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(lstress    ,1,MPI_LOGICAL,0,comm3d,mpierr)
    idtav = dtav/tres
    itimeav = timeav/tres

    tnext      = idtav   +btime
    tnextwrite = itimeav +btime
    nsamples = itimeav/idtav

    if(.not.(lstress)) return
    dt_lim = min(dt_lim,tnext)

    if (abs(timeav/dtav-nsamples)>1e-4) then
      stop 'timeav must be a integer multiple of dtav'
    end if
    if (.not. ladaptive .and.abs( dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if

! Time averaged budget terms

    allocate(tstr(k1,3,3),tshr(k1,3,3),tttr(k1,3,3),tbuo(k1,3,3),tcor(k1,3,3),tptr(k1,3,3),&
    tdis(k1,3,3),tadv(k1,3,3))

    allocate(rstrb(k1,3,3),rstre(k1,3,3))

    allocate(str_stor(k1,3,3),str_budg(k1,3,3),str_res(k1,3,3))

!Setting time mean variables to zero

    tstr=0.;tshr=0.;tttr=0.;tbuo=0.;tcor=0.;tptr=0.;tdis=0.;tadv=0.
    rstrb=0.;rstre=0.;str_stor=0.;str_budg=0.;str_res=0.
    lstressb=.false.

   !Preparing output files
    if(myid==0)then
       open (ifoutput,file='stressbudget.'//cexpnr,status='replace')
       close (ifoutput)
    endif

    if (lnetcdf) then
      idtav = idtav_prof
      itimeav = itimeav_prof
      tnext      = idtav+btime
      tnextwrite = itimeav+btime
      nsamples = itimeav/idtav
      if (myid==0) then

        call ncinfo(ncname( 1,:),'str_res','Residual stress','m/s^2','tttt')
        call ncinfo(ncname( 2,:),'str_sto','Storage stress','m/s^2','tttt')
        call ncinfo(ncname( 3,:),'str_budg','Budget stress','m/s^2','tttt')
        call ncinfo(ncname( 4,:),'str_shr','Shear stress','m/s^2','tttt')
        call ncinfo(ncname( 5,:),'str_tr','Transport stress','m/s^2','tttt')
        call ncinfo(ncname( 6,:),'str_buo','Buoyancy stress','m/s^2','tttt')
        call ncinfo(ncname( 7,:),'str_ptr','Pressure stress','m/s^2','tttt')
        call ncinfo(ncname( 8,:),'str_dis','Dissipation stress','m/s^2','tttt')
        call ncinfo(ncname( 9,:),'str','Stress','m/s^2','tttt')
        call ncinfo(ncname( 10,:),'str_cor','Coriolis stress','m/s^2','tttt')
        call ncinfo(ncname( 11,:),'str_adv','Advection stress','m/s^2','tttt')

        call redefine_nc(ncid_prof)
        call define_nc(ncid_prof, NVar, ncname)
     end if

   end if

  end subroutine initstressbudget


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> General routine, does the timekeeping

  subroutine stressbudgetstat

    use modglobal, only : rk3step,timee, dt_lim
    implicit none

    if (.not. lstress) return
    if (rk3step/=3) return
    if(timee<tnext .and. timee<tnextwrite) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if
    if (timee>=tnext) then
      tnext = tnext+idtav
      call do_stressbudget
    end if
    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+itimeav
      call writestressbudget
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
  end subroutine stressbudgetstat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Performs the resolved budget calculations

  subroutine do_stressbudget
    use modglobal,  only : i1,i2,j1,j2,k1,ih,jh,imax,jmax,kmax,dzf,dzh, &
                          rslabs,cu,cv,iadv_thl,grav, &
                          dxi,dyi,dx2i,dy2i,om22,om23,iadv_cd2,iadv_5th,iadv_cd6,iadv_mom
    use modsurfdata, only : thvs, ustar
    use modsubgrid, only : ekm, diffu, diffv, diffw
    use modpois,    only : p
    use modfields,  only : u0,v0,w0,um,vm,wm,up,vp,wp,thl0h,thv0h,u0av,v0av,thl0
    use modmpi,     only : nprocs,comm3d,nprocs,my_real,mpi_sum,mpierr,slabsum, myid

    implicit none
    integer m,n,k

    real, allocatable, dimension (:)     :: w0av          ! mean w0

    real, allocatable, dimension (:,:,:) :: u0_dev         ! container for a u-field (deviation)
    real, allocatable, dimension (:,:,:) :: v0_dev         ! container for a v-field (deviation)
    real, allocatable, dimension (:,:,:) :: w0_dev         ! container for a w-field (deviation)

    real, allocatable, dimension (:,:,:) :: u_mean         ! container for a u-field (mean)
    real, allocatable, dimension (:,:,:) :: v_mean         ! container for a v-field (mean)
    real, allocatable, dimension (:,:,:) :: w_mean         ! container for a w-field (mean)
    real, allocatable, dimension (:,:,:) :: dumfield       ! temp storage

    real, allocatable, dimension (:,:,:) :: u0_stor    ! temporary storage for u0-field
    real, allocatable, dimension (:,:,:) :: v0_stor    ! temporary storage for v0-field
    real, allocatable, dimension (:,:,:) :: w0_stor    ! temporary storage for w0-field
    real, allocatable, dimension (:,:,:) :: u_term     ! container for a u-term in budget
    real, allocatable, dimension (:,:,:) :: v_term     ! container for a v-term in budget
    real, allocatable, dimension (:,:,:) :: w_term     ! container for a w-term in budget

    real, allocatable, dimension (:) :: uterm_avl, uterm_av
    real, allocatable, dimension (:) :: vterm_avl, vterm_av
    real, allocatable, dimension (:) :: wterm_avl, wterm_av
    real, allocatable, dimension (:) :: thv0h_avl, thv0h_av

    real, allocatable, dimension (:,:,:) :: fstrl   !stress
    real, allocatable, dimension (:,:,:) :: fshrl   !shear
    real, allocatable, dimension (:,:,:) :: fttrl   !turbulent transport
    real, allocatable, dimension (:,:,:) :: fbuol   !buoyancy production
    real, allocatable, dimension (:,:,:) :: fcorl   !Coriolis production
    real, allocatable, dimension (:,:,:) :: fptrl   !pressure-velocity interaction
    real, allocatable, dimension (:,:,:) :: fdisl   !disspation
    real, allocatable, dimension (:,:,:) :: fadvl   !advection

    real, allocatable, dimension (:,:,:) :: fstr   !stress
    real, allocatable, dimension (:,:,:) :: fshr   !shear
    real, allocatable, dimension (:,:,:) :: fttr   !turbulent transport
    real, allocatable, dimension (:,:,:) :: fbuo   !buoyancy production
    real, allocatable, dimension (:,:,:) :: fcor   !Coriolis production
    real, allocatable, dimension (:,:,:) :: fptr   !pressure-velocity interaction
    real, allocatable, dimension (:,:,:) :: fdis   !disspation
    real, allocatable, dimension (:,:,:) :: fadv   !advection

    real, allocatable, dimension (:,:,:) :: rstrel  ! local storage for the stress at end of averaging period



    !********************************************************
    ! 0.    CALCULATE SLAB AVERAGED HORIZONTAL WIND SPEED
    !********************************************************

    !***************************************************************
    ! 0.0 Allocate variables
    !***************************************************************

    allocate(dumfield(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(uterm_avl(k1),vterm_avl(k1),wterm_avl(k1),thv0h_avl(k1))
    allocate(uterm_av(k1),vterm_av(k1),wterm_av(k1),thv0h_av(k1))
    allocate(w0av(k1))
    allocate(u0_dev(2-ih:i1+ih,2-jh:j1+jh,k1),v0_dev(2-ih:i1+ih,2-jh:j1+jh,k1),w0_dev(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(u_mean(2-ih:i1+ih,2-jh:j1+jh,k1),v_mean(2-ih:i1+ih,2-jh:j1+jh,k1),w_mean(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(u0_stor(2-ih:i1+ih,2-jh:j1+jh,k1),v0_stor(2-ih:i1+ih,2-jh:j1+jh,k1),w0_stor(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(u_term(2-ih:i1+ih,2-jh:j1+jh,k1),v_term(2-ih:i1+ih,2-jh:j1+jh,k1),w_term(2-ih:i1+ih,2-jh:j1+jh,k1))

    allocate(fstr(k1,3,3),fshr(k1,3,3),fttr(k1,3,3),fbuo(k1,3,3),fcor(k1,3,3),fptr(k1,3,3),&
    fdis(k1,3,3),fadv(k1,3,3))
    allocate(fstrl(k1,3,3),fshrl(k1,3,3),fttrl(k1,3,3),fbuol(k1,3,3),fcorl(k1,3,3),fptrl(k1,3,3),&
    fdisl(k1,3,3),fadvl(k1,3,3))

    allocate(rstrel(k1,3,3))


    !***************************************************************
    ! 0.1 Reset arrays for slab averages and calculate slab averages
    !***************************************************************

    w0av=0.
    fstrl=0.;fbuol=0.;fadvl=0.;fcorl=0.;fdisl=0.;fshrl=0.;fptrl=0.;fttrl=0.
    fstr=0.;fbuo=0.;fadv=0.;fcor=0.;fdis=0.;fshr=0.;fptr=0.;fttr=0.
    dumfield=0.
    rstrel=0.

    call slabsum(w0av  ,1,k1,w0  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)

    w0av = w0av / rslabs

    !**************************************************
    ! 0.2 Store the wind fields in u0, v0 and w0
    !**************************************************

    u0_stor = u0
    v0_stor = v0
    w0_stor = w0

    !*************************************************
    ! 0.3 Deviation wind field
    !*************************************************

    do k=1,k1
       u0_dev(:,:,k) = u0(:,:,k) - (u0av(k) - cu)
       v0_dev(:,:,k) = v0(:,:,k) - (v0av(k) - cv)
       w0_dev(:,:,k) = w0(:,:,k) - w0av(k)
    enddo

    ! Mean wind field

    do k=1,k1
       u_mean(:,:,k) = u0av(k)
       v_mean(:,:,k) = v0av(k)
       w_mean(:,:,k) = w0av(k)
    enddo

    !******************************************************************
    ! 1.    RESET ARRAYS CONTAINING TIME AVERAGES EACH TIMEAV SECONDS
    !*****************************************************************

    !***********************************************************
    ! 1.2  Calculate stress
    !***********************************************************


    do k=2,k1
    !   uu-stress
    fstrl(k,1,1) = sum((u0(2:i1,2:j1,k)-(u0av(k)-cu))**2)
    !   vv-stress
    fstrl(k,2,2) = sum((v0(2:i1,2:j1,k)-(v0av(k)-cv))**2)
    !   ww-stress
    fstrl(k,3,3) = sum((w0(2:i1,2:j1,k)-w0av(k))**2)
    !   uv-stress
    fstrl(k,1,2) = sum(0.5*(u0(2:i1,2:j1,k)-(u0av(k)-cu) +     &
                            u0(2:i1,1:jmax,k)-(u0av(k)-cu))*     &
                       0.5*(v0(2:i1,2:j1,k)-(v0av(k)-cv) +     &
                            v0(1:imax,2:j1,k)-(v0av(k)-cv)))
    !   vu-stress
    fstrl(k,2,1) = fstrl(k,1,2)
    !   uw-stress
    fstrl(k,1,3) = sum(0.5*(u0(2:i1,2:j1,k)-(u0av(k)-cu) +     &
                            u0(2:i1,2:j1,k-1)-(u0av(k-1)-cu))* &
                       0.5*(w0(2:i1,2:j1,k)-w0av(k) +     &
                            w0(1:imax,2:j1,k)-w0av(k)))
    !   wu-stress
    fstrl(k,3,1) = fstrl(k,1,3)
    !   vw-stress
    fstrl(k,2,3) = sum(0.5*(v0(2:i1,2:j1,k)-(v0av(k)-cv) +     &
                            v0(2:i1,2:j1,k-1)-(v0av(k-1)-cv))* &
                       0.5*(w0(2:i1,2:j1,k)-w0av(k) +     &
                            w0(2:i1,1:jmax,k)-w0av(k)))
    !   wv-stress
    fstrl(k,3,2) = fstrl(k,2,3)

    enddo

    fstrl(1,1,1) = fstrl(2,1,1)
    fstrl(1,2,2) = fstrl(2,2,2)
    fstrl(1,3,3) = 0.
    fstrl(1,1,2) = fstrl(2,1,2)
    fstrl(1,2,1) = fstrl(1,1,2)
    fstrl(1,1,3) = 0.
    fstrl(1,3,1) = fstrl(1,1,3)
    fstrl(1,2,3) = 0.
    fstrl(1,3,2) = fstrl(1,2,3)

    fstrl = fstrl/rslabs

    do m=1,3
      do n=1,3
        call MPI_ALLREDUCE(fstrl(:,m,n), fstr(:,m,n), k1,  MY_REAL, &
                           MPI_SUM, comm3d,mpierr)
      enddo
    enddo


    !***************************************************************
    ! 1.3  Calculate stress at the beginning of the averaging period
    !***************************************************************

    if (.not.(lstressb)) then
       rstrb = fstr
       lstressb = .true.
    endif


    !************************************************************
    ! 2.0    CALCULATE SLAB AVERAGES OF ALL TERMS
    !************************************************************

    !************************************************************
    ! 2.1    Buoyancy production terms
    !(see routine 'forces' for how term in w-momentum equation is treated)
    !**************************************************************

    call cyclicx(thv0h)

    do k=1,kmax
       thv0h_avl(k) =  sum(thv0h(2:i1,2:j1,k))/rslabs
    enddo

    call MPI_ALLREDUCE(thv0h_avl, thv0h_av, k1,  MY_REAL, &
                           MPI_SUM, comm3d,mpierr)

    do k=1,k1
      w_term(:,:,k) = (grav/thvs) * (thv0h(2-ih:i1+ih,2-jh:j1+jh,k) - thv0h_av(k))
    enddo

    ! Lateral BC; we don't need lower BC

    call cyclicx(w_term)

    dumfield(2:i1,2:j1,2:k1) =  0.25*(u0_dev(2:i1,2:j1,1:kmax)+u0_dev(2:i1,2:j1,2:k1))*  &
                                     (w_term(2:i1,2:j1,2:k1)+w_term(1:imax,2:j1,2:k1))

    dumfield(:,:,1) = 0.

    do k=1,k1
       fbuol(k,1,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    fbuol(:,3,1) = fbuol(:,1,3)

    dumfield(2:i1,2:j1,2:k1) = 0.25*(v0_dev(2:i1,2:j1,1:kmax)+v0_dev(2:i1,2:j1,2:k1))*   &
                                    (w_term(2:i1,2:j1,2:k1)+w_term(2:i1,1:jmax,2:k1))
    dumfield(:,:,1) = 0.

    do k=1,k1
       fbuol(k,2,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    fbuol(:,3,2) = fbuol(:,2,3)

    dumfield = 2.*(w0_dev * w_term)

    do k=1,k1
       fbuol(k,3,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    do m=1,3
      do n=1,3
        call MPI_ALLREDUCE(fbuol(:,m,n), fbuo(:,m,n), k1,  MY_REAL, &
                           MPI_SUM, comm3d,mpierr)
      enddo
    enddo


    !*********************************************************************
    !     2.2    Shear production terms
    !*********************************************************************

    ! Note that the minus sign is already taken care of in the advecx_ routines

    u0 = u0_dev
    v0 = v0_dev
    w0 = w0_dev

    u_term=0.
    v_term=0.
    w_term=0.

    select case(iadv_mom)
      case(iadv_cd2)
        call advecu_2nd(u_mean,u_term)
        call advecv_2nd(v_mean,v_term)
        call advecw_2nd(w_mean,w_term)
      case(iadv_5th)
        call advecu_5th(u_mean,u_term)
        call advecv_5th(v_mean,v_term)
        call advecw_5th(w_mean,w_term)
      case(iadv_cd6)
        call advecu_6th(u_mean,u_term)
        call advecv_6th(v_mean,v_term)
        call advecw_6th(w_mean,w_term)
      case default
          stop "Unknown advection scheme "
    end select

    u0 = u0_stor
    v0 = v0_stor
    w0 = w0_stor

    call cyclicx(u_term)
    call cyclicx(v_term)
    call cyclicx(w_term)

    dumfield = 2.*u0_dev*u_term

    do k=1,k1
       fshrl(k,1,1) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    dumfield = 2.*v0_dev*v_term

    do k=1,k1
       fshrl(k,2,2) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    dumfield = 2.*w0_dev*w_term

    do k=1,k1
       fshrl(k,3,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    dumfield(2:i1,2:j1,:)    = 0.25*(u0_dev(2:i1,2:j1,:)+u0_dev(2:i1,1:jmax,:))*          &
                                    (v_term(2:i1,2:j1,:)+v_term(1:imax,2:j1,:))   +     &
                               0.25*(u_term(2:i1,2:j1,:)+u_term(2:i1,1:jmax,:))*        &
                                    (v0_dev(2:i1,2:j1,:)+v0_dev(1:imax,2:j1,:))

    do k=1,k1
       fshrl(k,1,2) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    fshrl(:,2,1) = fshrl(:,1,2)

    dumfield(2:i1,2:j1,2:k1) = 0.25*(u0_dev(2:i1,2:j1,2:k1)+u0_dev(2:i1,2:j1,1:kmax))*       &
                                    (w_term(2:i1,2:j1,2:k1)+w_term(1:imax,2:j1,2:k1))  +   &
                               0.25*(u_term(2:i1,2:j1,2:k1)+u_term(2:i1,2:j1,1:kmax))*     &
                                    (w0_dev(2:i1,2:j1,2:k1)+w0_dev(1:imax,2:j1,2:k1))

    dumfield(:,:,1) = 0.

    do k=1,k1
       fshrl(k,1,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    fshrl(:,3,1) = fshrl(:,1,3)

    dumfield(2:i1,2:j1,2:k1) = 0.25*(v0_dev(2:i1,2:j1,2:k1)+v0_dev(2:i1,2:j1,1:kmax))*          &
                                    (w_term(2:i1,2:j1,2:k1)+w_term(2:i1,1:jmax,2:k1))  +      &
                               0.25*(v_term(2:i1,2:j1,2:k1)+v_term(2:i1,2:j1,1:kmax))*        &
                                    (w0_dev(2:i1,2:j1,2:k1)+w0_dev(2:i1,1:jmax,2:k1))

    dumfield(:,:,1) = 0.

    do k=1,k1
       fshrl(k,2,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    fshrl(:,3,2) = fshrl(:,2,3)

    do m=1,3
        do n=1,3
          call MPI_ALLREDUCE(fshrl(:,m,n), fshr(:,m,n), k1,  MY_REAL, &
                             MPI_SUM, comm3d,mpierr)
        enddo
    enddo

   !*********************************************************
   !     2.3    Turbulent transport terms
   !*************************************************

   ! Note that the minus sign is already taken care of in the advecx routines

   u0 = u0_dev
   v0 = v0_dev
   w0 = w0_dev

   u_term=0.
   v_term=0.
   w_term=0.

   select case(iadv_mom)
      case(iadv_cd2)
        call advecu_2nd(u0_dev,u_term)
        call advecv_2nd(v0_dev,v_term)
        call advecw_2nd(w0_dev,w_term)
      case(iadv_5th)
        call advecu_5th(u0_dev,u_term)
        call advecv_5th(v0_dev,v_term)
        call advecw_5th(w0_dev,w_term)
      case(iadv_cd6)
        call advecu_6th(u0_dev,u_term)
        call advecv_6th(v0_dev,v_term)
        call advecw_6th(w0_dev,w_term)
      case default
          stop "Unknown advection scheme "
   end select

   u0 = u0_stor
   v0 = v0_stor
   w0 = w0_stor

   call cyclicx(u_term)
   call cyclicx(v_term)
   call cyclicx(w_term)

   dumfield = 2.*u0_dev*u_term

   do k=1,k1
     fttrl(k,1,1) = sum(dumfield(2:i1,2:j1,k))/rslabs
   enddo

   dumfield = 2.*v0_dev*v_term

   do k=1,k1
     fttrl(k,2,2) = sum(dumfield(2:i1,2:j1,k))/rslabs
   enddo

   dumfield = 2.*w0_dev*w_term

   do k=1,k1
     fttrl(k,3,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
   enddo

   dumfield(2:i1,2:j1,:) = 0.25*(u0_dev(2:i1,2:j1,:)+u0_dev(2:i1,1:jmax,:))*  &
                                (v_term(2:i1,2:j1,:)+v_term(1:imax,2:j1,:))   +     &
                           0.25*(u_term(2:i1,2:j1,:)+u_term(2:i1,1:jmax,:))*        &
                                (v0_dev(2:i1,2:j1,:)+v0_dev(1:imax,2:j1,:))

   do k=1,k1
     fttrl(k,1,2) = sum(dumfield(2:i1,2:j1,k))/rslabs
   enddo

   fttrl(:,2,1) = fttrl(:,1,2)

   dumfield(2:i1,2:j1,2:k1) = 0.25*(u0_dev(2:i1,2:j1,2:k1)+u0_dev(2:i1,2:j1,1:kmax))* &
                                   (w_term(2:i1,2:j1,2:k1)+w_term(1:imax,2:j1,2:k1))  +   &
                              0.25*(u_term(2:i1,2:j1,2:k1)+u_term(2:i1,2:j1,1:kmax))*     &
                                   (w0_dev(2:i1,2:j1,2:k1)+w0_dev(1:imax,2:j1,2:k1))

   dumfield(:,:,1) = 0.

   do k=1,k1
     fttrl(k,1,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
   enddo

   fttrl(:,3,1) = fttrl(:,1,3)

   dumfield(2:i1,2:j1,2:k1) = 0.25*(v0_dev(2:i1,2:j1,2:k1)+v0_dev(2:i1,2:j1,1:kmax))*   &
                                   (w_term(2:i1,2:j1,2:k1)+w_term(2:i1,1:jmax,2:k1))   +     &
                              0.25*(v_term(2:i1,2:j1,2:k1)+v_term(2:i1,2:j1,1:kmax))*        &
                                   (w0_dev(2:i1,2:j1,2:k1)+w0_dev(2:i1,1:jmax,2:k1))

   dumfield(:,:,1) = 0.

   do k=1,k1
     fttrl(k,2,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
   enddo

   fttrl(:,3,2) = fttrl(:,2,3)

   do m=1,3
      do n=1,3
        call MPI_ALLREDUCE(fttrl(:,m,n), fttr(:,m,n), k1,  MY_REAL, &
                          MPI_SUM, comm3d, mpierr)
      enddo
   enddo

    !********************************************************
    !     2.3b    Advection (should be zero)
    !********************************************************

! Note that the minus sign is already taken care of in the advecx routines

    u0 = u_mean
    v0 = v_mean
    w0 = w_mean

    u_term=0.
    v_term=0.
    w_term=0.

    select case(iadv_mom)
      case(iadv_cd2)
        call advecu_2nd(u0_dev,u_term)
        call advecv_2nd(v0_dev,v_term)
        call advecw_2nd(w0_dev,w_term)
      case(iadv_5th)
        call advecu_5th(u0_dev,u_term)
        call advecv_5th(v0_dev,v_term)
        call advecw_5th(w0_dev,w_term)
      case(iadv_cd6)
        call advecu_6th(u0_dev,u_term)
        call advecv_6th(v0_dev,v_term)
        call advecw_6th(w0_dev,w_term)
      case default
          stop "Unknown advection scheme "
    end select

    u0 = u0_stor
    v0 = v0_stor
    w0 = w0_stor

    call cyclicx(u_term)
    call cyclicx(v_term)
    call cyclicx(w_term)

    dumfield = 2.*u0_dev*u_term

    do k=1,k1
       fadvl(k,1,1) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    dumfield = 2.*v0_dev*v_term

    do k=1,k1
       fadvl(k,2,2) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    dumfield = 2.*w0_dev*w_term

    do k=1,k1
       fadvl(k,3,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    dumfield(2:i1,2:j1,:) =  0.25*(u0_dev(2:i1,2:j1,:)+u0_dev(2:i1,1:jmax,:))*          &
                                  (v_term(2:i1,2:j1,:)+v_term(1:imax,2:j1,:))   +     &
                             0.25*(u_term(2:i1,2:j1,:)+u_term(2:i1,1:jmax,:))*        &
                                  (v0_dev(2:i1,2:j1,:)+v0_dev(1:imax,2:j1,:))

    do k=1,k1
       fadvl(k,1,2) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    fadvl(:,2,1) = fadvl(:,1,2)

    dumfield(2:i1,2:j1,2:k1) = 0.25*(u0_dev(2:i1,2:j1,2:k1)+u0_dev(2:i1,2:j1,1:kmax))*       &
                                    (w_term(2:i1,2:j1,2:k1)+w_term(1:imax,2:j1,2:k1))  +   &
                               0.25*(u_term(2:i1,2:j1,2:k1)+u_term(2:i1,2:j1,1:kmax))*     &
                                    (w0_dev(2:i1,2:j1,2:k1)+w0_dev(1:imax,2:j1,2:k1))

    dumfield(:,:,1) = 0.

    do k=1,k1
       fadvl(k,1,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    fadvl(:,3,1) = fadvl(:,1,3)

    dumfield(2:i1,2:j1,2:k1) = 0.25*(v0_dev(2:i1,2:j1,2:k1)+v0_dev(2:i1,2:j1,1:kmax))*          &
                                    (w_term(2:i1,2:j1,2:k1)+w_term(2:i1,1:jmax,2:k1))   +     &
                               0.25*(v_term(2:i1,2:j1,2:k1)+v_term(2:i1,2:j1,1:kmax))*        &
                                    (w0_dev(2:i1,2:j1,2:k1)+w0_dev(2:i1,1:jmax,2:k1))

    dumfield(:,:,1) = 0.

    do k=1,k1
       fadvl(k,2,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    fadvl(:,3,2) = fadvl(:,2,3)

    do m=1,3
        do n=1,3
          call MPI_ALLREDUCE(fadvl(:,m,n), fadv(:,m,n), k1,  MY_REAL, &
                             MPI_SUM, comm3d,mpierr)
        enddo
    enddo


    !*******************************************************************
    ! 2.4    Coriolis term
    !******************************************************************

    u_term(2:i1,2:j1,1:kmax) =   &
            +(v0_dev(2:i1,2:j1,1:kmax)+v0_dev(2:i1,3:j2,1:kmax)+v0_dev(1:imax,2:j1,1:kmax)+v0_dev(1:imax,3:j2,1:kmax))*om23*0.25 &
            -(w0_dev(2:i1,2:j1,1:kmax)+w0_dev(2:i1,2:j1,2:k1)+w0_dev(1:imax,2:j1,2:k1)+w0_dev(1:imax,2:j1,1:kmax))*om22*0.25

    v_term(2:i1,2:j1,:) =  &
            -(u0_dev(2:i1,2:j1,:)+u0_dev(2:i1,1:j1-1,:)+u0_dev(3:i2,1:j1-1,:)+u0_dev(3:i2,2:j1,:))*om23*0.25

    do k=2,k1
       w_term(2:i1,2:j1,k) =   &
                    om22 * 0.25*( (dzf(k-1) * (u0_dev(2:i1,2:j1,k-1)  + u0_dev(3:i2,2:j1,k-1) )     &
                                 + dzf(k)  *  (u0_dev(2:i1,2:j1,k) + u0_dev(3:i2,2:j1,k))  ) / dzh(k) )
    enddo

    w_term(:,:,1) = 0.

    call cyclicx(u_term)
    call cyclicx(v_term)
    call cyclicx(w_term)

    dumfield = 2.* u0_dev * u_term

    do k=1,k1
       fcorl(k,1,1) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    dumfield = 2.* v0_dev * v_term

    do k=1,k1
       fcorl(k,2,2) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    dumfield = 2.* w0_dev * w_term

    do k=1,k1
       fcorl(k,3,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    dumfield(2:i1,2:j1,:) =     0.25*(u0_dev(2:i1,2:j1,:)+u0_dev(2:i1,1:jmax,:))*          &
                                     (v_term(2:i1,2:j1,:)+v_term(1:imax,2:j1,:))   +     &
                                0.25*(u_term(2:i1,2:j1,:)+u_term(2:i1,1:jmax,:))*        &
                                     (v0_dev(2:i1,2:j1,:)+v0_dev(1:imax,2:j1,:))

    do k=1,k1
       fcorl(k,1,2) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    fcorl(:,2,1) = fcorl(:,1,2)

    dumfield(2:i1,2:j1,2:k1) =    0.25*(u0_dev(2:i1,2:j1,2:k1)+u0_dev(2:i1,2:j1,1:kmax))*       &
                                       (w_term(2:i1,2:j1,2:k1)+w_term(1:imax,2:j1,2:k1))  +   &
                                  0.25*(u_term(2:i1,2:j1,2:k1)+u_term(2:i1,2:j1,1:kmax))*     &
                                       (w0_dev(2:i1,2:j1,2:k1)+w0_dev(1:imax,2:j1,2:k1))

    do k=1,k1
       fcorl(k,1,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    fcorl(:,3,1) = fcorl(:,1,3)

    dumfield(2:i1,2:j1,2:k1) =    0.25*(v0_dev(2:i1,2:j1,2:k1)+v0_dev(2:i1,2:j1,1:kmax))*          &
                                       (w_term(2:i1,2:j1,2:k1)+w_term(2:i1,1:jmax,2:k1))  +      &
                                  0.25*(v_term(2:i1,2:j1,2:k1)+v_term(2:i1,2:j1,1:kmax))*        &
                                       (w0_dev(2:i1,2:j1,2:k1)+w0_dev(2:i1,1:jmax,2:k1))

    do k=1,k1
       fcorl(k,2,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    fcorl(:,3,2) = fcorl(:,2,3)

    do m=1,3
        do n=1,3
          call MPI_ALLREDUCE(fcorl(:,m,n), fcor(:,m,n), k1,  MY_REAL, &
                             MPI_SUM, comm3d,mpierr)
        enddo
    enddo

    !***************************************************************************
    !   2.4    Pressure fluctuation terms
    !**************************************************************************


    u_term(2:i1,2:j1,1:k1) = -(p(2:i1,2:j1,1:k1)-p(1:imax,2:j1,1:k1))*dxi
    v_term(2:i1,2:j1,1:k1) = -(p(2:i1,2:j1,1:k1)-p(2:i1,1:jmax,1:k1))*dyi

    do k=2,k1
       w_term(2:i1,2:j1,k) = -(p(2:i1,2:j1,k)-p(2:i1,2:j1,k-1))/dzh(k)
    enddo

    call cyclicx(u_term)
    call cyclicx(v_term)
    call cyclicx(w_term)

    w_term(:,:,1) = 0.

    dumfield = 2.*u_term*u0_dev
    do k=1,k1
       fptrl(k,1,1) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    dumfield = 2.*v_term*v0_dev
    do k=1,k1
       fptrl(k,2,2) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    dumfield = 2.*w_term*w0_dev
    do k=1,k1
       fptrl(k,3,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    dumfield(2:i1,2:j1,2:k1) = 0.25*(u0_dev(2:i1,2:j1,2:k1) + u0_dev(2:i1,1:jmax,2:k1))*      &
                                    (v_term(2:i1,2:j1,2:k1) + v_term(1:imax,2:j1,2:k1)) +   &
                               0.25*(u_term(2:i1,2:j1,2:k1) + u_term(2:i1,1:jmax,2:k1))*      &
                                    (v0_dev(2:i1,2:j1,2:k1) + v0_dev(1:imax,2:j1,2:k1))
    do k=1,k1
       fptrl(k,1,2) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo
    fptrl(:,2,1) = fptrl(:,1,2)

    dumfield(2:i1,2:j1,2:k1) = 0.25*(u0_dev(2:i1,2:j1,2:k1) + u0_dev(2:i1,2:j1,1:kmax))*      &
                                    (w_term(2:i1,2:j1,2:k1) + w_term(1:imax,2:j1,2:k1)) +   &
                               0.25*(u_term(2:i1,2:j1,2:k1) + u_term(2:i1,2:j1,1:kmax))*      &
                                    (w0_dev(2:i1,2:j1,2:k1) + w0_dev(1:imax,2:j1,2:k1))
    do k=1,k1
       fptrl(k,1,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo
    fptrl(:,3,1) = fptrl(:,1,3)

    dumfield(2:i1,2:j1,2:k1) = 0.25*(v0_dev(2:i1,2:j1,2:k1) + v0_dev(2:i1,2:j1,1:kmax))*      &
                                    (w_term(2:i1,2:j1,2:k1) + w_term(2:i1,1:jmax,2:k1)) +   &
                               0.25*(v_term(2:i1,2:j1,2:k1) + v_term(2:i1,2:j1,1:kmax))*      &
                                    (w0_dev(2:i1,2:j1,2:k1) + w0_dev(2:i1,1:jmax,2:k1))
    do k=1,k1
       fptrl(k,2,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo
    fptrl(:,3,2) = fptrl(:,2,3)

    do m=1,3
        do n=1,3
           call MPI_ALLREDUCE(fptrl(:,m,n), fptr(:,m,n), k1,  MY_REAL, &
                           MPI_SUM, comm3d, mpierr)
       enddo
    enddo

    !********************************************************************************
    !   2.5    Dissipation terms
    !******************************************************************************


    u_term = 0.
    v_term = 0.
    w_term = 0.

    call diffu(u_term)
    call diffv(v_term)
    call diffw(w_term)

    call cyclicx(u_term)
    call cyclicx(v_term)
    call cyclicx(w_term)

    do k=1,k1
       uterm_avl(k) = sum(u_term(2:i1,2:j1,k))/rslabs
       vterm_avl(k) = sum(v_term(2:i1,2:j1,k))/rslabs
       wterm_avl(k) = sum(w_term(2:i1,2:j1,k))/rslabs
    enddo

    call MPI_ALLREDUCE(uterm_avl, uterm_av, k1,  MY_REAL, &
                           MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(vterm_avl, vterm_av, k1,  MY_REAL, &
                           MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(wterm_avl, wterm_av, k1,  MY_REAL, &
                           MPI_SUM, comm3d, mpierr)

    do k=1,k1
       u_term(:,:,k) = u_term(:,:,k) - uterm_av(k)
       v_term(:,:,k) = v_term(:,:,k) - vterm_av(k)
       w_term(:,:,k) = w_term(:,:,k) - wterm_av(k)
    enddo


    dumfield = 2.*u_term*u0_dev
    do k=1,k1
       fdisl(k,1,1) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    dumfield = 2.*v_term*v0_dev
    do k=1,k1
       fdisl(k,2,2) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    dumfield = 2.*w_term*w0_dev
    do k=1,k1
       fdisl(k,3,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo

    dumfield(2:i1,2:j1,2:k1) = 0.25*(u0_dev(2:i1,2:j1,2:k1) + u0_dev(2:i1,1:jmax,2:k1))*      &
                                    (v_term(2:i1,2:j1,2:k1) + v_term(1:imax,2:j1,2:k1)) +   &
                               0.25*(u_term(2:i1,2:j1,2:k1) + u_term(2:i1,1:jmax,2:k1))*      &
                                    (v0_dev(2:i1,2:j1,2:k1) + v0_dev(1:imax,2:j1,2:k1))
    do k=1,k1
       fdisl(k,1,2) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo
    fdisl(:,2,1) = fdisl(:,1,2)

    dumfield(2:i1,2:j1,2:k1) = 0.25*(u0_dev(2:i1,2:j1,2:k1) + u0_dev(2:i1,2:j1,1:kmax))*      &
                                    (w_term(2:i1,2:j1,2:k1) + w_term(1:imax,2:j1,2:k1)) +   &
                               0.25*(u_term(2:i1,2:j1,2:k1) + u_term(2:i1,2:j1,1:kmax))*      &
                                    (w0_dev(2:i1,2:j1,2:k1) + w0_dev(1:imax,2:j1,2:k1))
    do k=1,k1
       fdisl(k,1,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo
    fdisl(:,3,1) = fdisl(:,1,3)

    dumfield(2:i1,2:j1,2:k1) = 0.25*(v0_dev(2:i1,2:j1,2:k1) + v0_dev(2:i1,2:j1,1:kmax))*      &
                                    (w_term(2:i1,2:j1,2:k1) + w_term(2:i1,1:jmax,2:k1)) +   &
                               0.25*(v_term(2:i1,2:j1,2:k1) + v_term(2:i1,2:j1,1:kmax))*      &
                                    (w0_dev(2:i1,2:j1,2:k1) + w0_dev(2:i1,1:jmax,2:k1))
    do k=1,k1
       fdisl(k,2,3) = sum(dumfield(2:i1,2:j1,k))/rslabs
    enddo
    fdisl(:,3,2) = fdisl(:,2,3)

    do m=1,3
        do n=1,3
           call MPI_ALLREDUCE(fdisl(:,m,n), fdis(:,m,n), k1,  MY_REAL, &
                           MPI_SUM, comm3d, mpierr)
       enddo
    enddo


    !*******************************************************************
    ! 3  CALCULATION OF TIME MEAN
    !*******************************************************************

    do k=1,kmax
      tstr(k,:,:) = tstr(k,:,:) + fstr(k,:,:) * dtav
      tbuo(k,:,:) = tbuo(k,:,:) + fbuo(k,:,:) * dtav
      tshr(k,:,:) = tshr(k,:,:) + fshr(k,:,:) * dtav
      tttr(k,:,:) = tttr(k,:,:) + fttr(k,:,:) * dtav
      tcor(k,:,:) = tcor(k,:,:) + fcor(k,:,:) * dtav
      tptr(k,:,:) = tptr(k,:,:) + fptr(k,:,:) * dtav
      tdis(k,:,:) = tdis(k,:,:) + fdis(k,:,:) * dtav
      tadv(k,:,:) = tadv(k,:,:) + fadv(k,:,:) * dtav
    end do

    !*********************************************************************
    ! 3.1 Calculate stress at the beginning of the averaging period
    !*********************************************************************

    do k=2,k1
  !   uu-stress
  rstrel(k,1,1) = sum((u0(2:i1,2:j1,k)-(u0av(k)-cu))**2)
  !   vv-stress
  rstrel(k,2,2) = sum((v0(2:i1,2:j1,k)-(v0av(k)-cv))**2)
  !   ww-stress
  rstrel(k,3,3) = sum((w0(2:i1,2:j1,k)-w0av(k))**2)
  !   uv-stress
  rstrel(k,1,2) = sum(0.5*(u0(2:i1,2:j1,k)-(u0av(k)-cu) +     &
        u0(2:i1,1:jmax,k)-(u0av(k)-cu))*     &
        0.5*(v0(2:i1,2:j1,k)-(v0av(k)-cv) +     &
        v0(1:imax,2:j1,k)-(v0av(k)-cv)))
  !   vu-stress
  rstrel(k,2,1) = rstrel(k,1,2)
  !   uw-stress
  rstrel(k,1,3) = sum(0.5*(u0(2:i1,2:j1,k)-(u0av(k)-cu) +     &
        u0(2:i1,2:j1,k-1)-(u0av(k-1)-cu))* &
        0.5*(w0(2:i1,2:j1,k)-w0av(k) +     &
        w0(1:imax,2:j1,k)-w0av(k)))
  !   wu-stress
  rstrel(k,3,1) = rstrel(k,1,3)
  !   vw-stress
  rstrel(k,2,3) = sum(0.5*(v0(2:i1,2:j1,k)-(v0av(k)-cv) +     &
        v0(2:i1,2:j1,k-1)-(v0av(k-1)-cv))* &
        0.5*(w0(2:i1,2:j1,k)-w0av(k) +     &
        w0(2:i1,1:jmax,k)-w0av(k)))
  !   wv-stress
  rstrel(k,3,2) = rstrel(k,2,3)
    enddo
    rstrel(1,1,1) = rstrel(2,1,1)
    rstrel(1,2,2) = rstrel(2,2,2)
    rstrel(1,3,3) = 0.
    rstrel(1,1,2) = rstrel(2,1,2)
    rstrel(1,2,1) = rstrel(1,1,2)
    rstrel(1,1,3) = 0.
    rstrel(1,3,1) = rstrel(1,1,3)
    rstrel(1,2,3) = 0.
    rstrel(1,3,2) = rstrel(1,2,3)

    rstrel = rstrel/rslabs

    do m=1,3
      do n=1,3
        call MPI_ALLREDUCE(rstrel(:,m,n), rstre(:,m,n), k1,  MY_REAL, &
                           MPI_SUM, comm3d,mpierr)
      enddo
    enddo

    !**************************************************
    ! 4 DEALLOCATE TEMPORARY VARIABLES
    !**************************************************

    deallocate(u0_dev,v0_dev,w0_dev)
    deallocate(w0av)
    deallocate(u_mean,v_mean,w_mean)
    deallocate(dumfield)
    deallocate(u0_stor,v0_stor,w0_stor)
    deallocate(u_term,v_term,w_term)
    deallocate(uterm_av,vterm_av,wterm_av,thv0h_av)
    deallocate(uterm_avl,vterm_avl,wterm_avl,thv0h_avl)
    deallocate(fshr,fttr,fbuo,fptr,fdis,fstr,fcor,fadv)
    deallocate(fshrl,fttrl,fbuol,fptrl,fdisl,fstrl,fcorl,fadvl)
    deallocate(rstrel)

  end subroutine do_stressbudget

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine writestressbudget

    use modglobal, only : kmax,k1,zf,zh,rtimee,cexpnr,ifoutput,rslabs,i1,j1
    use modmpi,    only : myid
    use modstat_nc,only : writestat_nc,lnetcdf
    use modgenstat, only: ncid_prof=>ncid,nrec_prof=>nrec

    implicit none  
    
    real,dimension(k1,3,3,nvar) :: vars
    integer nsecs, nhrs, nminut, k, m, n

    nsecs   = nint(rtimee)
    nhrs    = int(nsecs/3600)
    nminut  = int(nsecs/60)-nhrs*60
    nsecs   = mod(nsecs,60)

    do k=kmax,1,-1
       tstr(k,:,:) = tstr(k,:,:)/timeav
       tshr(k,:,:) = tshr(k,:,:)/timeav
       tbuo(k,:,:) = tbuo(k,:,:)/timeav
       tttr(k,:,:) = tttr(k,:,:)/timeav
       tcor(k,:,:) = tcor(k,:,:)/timeav
       tptr(k,:,:) = tptr(k,:,:)/timeav
       tdis(k,:,:) = tdis(k,:,:)/timeav
       tadv(k,:,:) = tadv(k,:,:)/timeav
       str_stor (k,:,:) = (rstre(k,:,:) - rstrb(k,:,:))/timeav
       str_budg (k,:,:) = tshr(k,:,:)+tbuo(k,:,:)+tttr(k,:,:)+tcor(k,:,:)+tptr(k,:,:)+tdis(k,:,:)+tadv(k,:,:)
       str_res  (k,:,:) = str_budg(k,:,:) - str_stor(k,:,:)
     end do

    if(myid == 0)then
       
       open (ifoutput,file='stressbudget.'//cexpnr,position='append')       
       
       do m=1,3
          do n=1,3
            if (n>=m) then
              write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
              '#--------------------------------------------------------' &
              ,'#',(timeav),'--- AVERAGING TIMESTEP --- ' &
              ,nhrs,':',nminut,':',nsecs &
              ,'   HRS:MIN:SEC AFTER INITIALIZATION '
              write (ifoutput,'(I2,A,I2,A)') m,',',n,'-budget'
              write (ifoutput,'(A/3A)') &
              '#--------------------------------------------------------' &
              ,'#  zfull[m] zhalf[m]  resid       store         budg        shear' &
              ,'         trsp         buoy         prsc         diss' &
              ,'         stress       cor          adv'        
              do k=kmax,1,-1
               write(ifoutput,'(f7.1,f7.1,11(e13.5))') &
                  zf(k),zh(k), str_res(k,m,n),str_stor(k,m,n),str_budg(k,m,n),tshr(k,m,n),tttr(k,m,n), &
                  tbuo(k,m,n),tptr(k,m,n),tdis(k,m,n),tstr(k,m,n), tcor(k,m,n), tadv(k,m,n)
              enddo
            endif ! end if(n>=m)
          enddo ! end m-loop
       enddo ! end n-loop
              
       write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
       '#--------------------------------------------------------' &
       ,'#',(timeav),'--- AVERAGING TIMESTEP --- ' &
       ,nhrs,':',nminut,':',nsecs &
       ,'   HRS:MIN:SEC AFTER INITIALIZATION '
       write (ifoutput,'(A)') 'stress-budget'
       write (ifoutput,'(A/3A)') &
       '#--------------------------------------------------------' &
       ,'#  zfull[m] zhalf[m]  resid       store         budg        shear' &
       ,'         trsp         buoy         prsc         diss' &
       ,'         stress       cor          adv'
       do k=kmax,1,-1
            write(ifoutput,'(f7.1,f7.1,11(e13.5))') &
               zf(k),zh(k),   0.5*(str_res(k,1,1) + str_res(k,2,2) + str_res(k,3,3)), &                
                              0.5*(str_stor(k,1,1) + str_stor(k,2,2) + str_stor(k,3,3)), &
                              0.5*(str_budg(k,1,1) + str_budg(k,2,2) + str_budg(k,3,3)), &
                              0.5*(tshr(k,1,1) + tshr(k,2,2) + tshr(k,3,3)), &
                              0.5*(tttr(k,1,1) + tttr(k,2,2) + tttr(k,3,3)), &
                              0.5*(tbuo(k,1,1) + tbuo(k,2,2) + tbuo(k,3,3)), & 
                              0.5*(tptr(k,1,1) + tptr(k,2,2) + tptr(k,3,3)), &
                              0.5*(tdis(k,1,1) + tdis(k,2,2) + tdis(k,3,3)), &
                              0.5*(tstr(k,1,1) + tstr(k,2,2) + tstr(k,3,3)), &
                              0.5*(tcor(k,1,1) + tcor(k,2,2) + tcor(k,3,3)), &
                              0.5*(tadv(k,1,1) + tadv(k,2,2) + tadv(k,3,3))
       enddo
       close(ifoutput)
     endif ! end if(myid==0)

  if (lnetcdf) then
        vars(:,:,:, 1) =str_res
        vars(:,:,:, 2) =str_stor
        vars(:,:,:, 3) =str_budg
        vars(:,:,:, 4) =tshr
        vars(:,:,:, 5) =tttr
        vars(:,:,:, 6) =tbuo
        vars(:,:,:, 7) =tptr
        vars(:,:,:, 8) =tdis
        vars(:,:,:, 9) =tstr
        vars(:,:,:,10) =tcor
        vars(:,:,:,11) =tadv
        call writestat_nc(ncid_prof,nvar,ncname,vars(1:kmax,:,:,:),nrec_prof,kmax,3,3)
  end if

!Setting time mean variables to zero

  tstr=0.;tshr=0.;tttr=0.;tbuo=0.;tcor=0.;tptr=0.;tdis=0.;tadv=0.
  rstrb=0.;rstre=0.;str_stor=0.;str_budg=0.;str_res=0.
  lstressb=.false.

  end subroutine writestressbudget


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine exitstressbudget

  use modmpi,    only : myid
  
  implicit none

    if(.not.(lstress)) return

    deallocate(str_res,str_stor,str_budg)
    deallocate(tshr,tttr,tbuo,tptr,tdis,tstr,tcor,tadv)
    deallocate(rstre,rstrb)

  end subroutine exitstressbudget


  subroutine cyclicx(field)

!-----------------------------------------------------------------|
!                                                                 |
!*** *cyclicx*  set lateral boundary conditions                  |
!      Based on cyclich                                           |
!                                                                 |
!      Hans Cuijpers   I.M.A.U.                                   |
!      Pier Siebesma   K.N.M.I.     22/02/1995                    |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!     Set cyclic boundary condition                               |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!             *cyclicx* is called from *modstress*               |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal, only : i1,i2,j1,j2,k1,ih,jh
  use modmpi,    only : excjs 
 
  implicit none

  real field(2-ih:i1+ih,2-jh:j1+jh,k1)
  integer m

  do m = 1, ih

  field(2-m,:,:)  = field(i2-m,:,:)
  field(i1+m,:,:) = field(1+m,:,:)

  end do

  call excjs( field          , 2,i1,2,j1,1,k1,ih,jh)

  return

  end subroutine cyclicx

  end module modstress
