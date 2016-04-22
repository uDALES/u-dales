!> \file modgenstat.f90
!!  Genstat calculates slab averages of several variables

!>
!!  Genstat calculates slab averages of several variables
!>
!!!!!!!  MODIFIED BY J.M. TOMAS       !!!!!!!
!! average1homog should be used when only the spanwise direction (j) is statist. homogeneous
!! average2homog should be used when the horizontal directions are statist. homogeneous
!!
!!  Written to fields.expnr, moments.expnr and flux1.expnr and flux2.expnr
!! If netcdf is true, this module leads the profiles.expnr.nc output
!!  \author Jasper Tomas, TU Delft, March 31 2014
!!  \author Hans Cuijpers, K.N.M.I.
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  \todo Documentation
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



module modgenstat

  !-----------------------------------------------------------------|
  !                                                                 |
  !*** *stattend*  calculates generic slabaveraged statistics       |
  !                                                                 |
  !      Pier Siebesma   K.N.M.I.     12/05/1995                    |
  !      Hans Cuijpers   I.M.A.U.     21/06/1995                    |
  !                                                                 |
  !     purpose.                                                    |
  !     --------                                                    |
  !                                                                 |
  !     stattend.f calculates:                                      |
  !                                                                 |
  !     * time averaged fieldsof theta_l, theta, theta_v, qt, qv   |
  !        ql, u and v                                              |
  !     * time averaged tendencies  of theta_l, theta, qt, qv, ql,  |
  !        u and v                                                  |
  !     * time averaged turbulent fluxes of theta_l, theta_v,       |
  !       theta, qt, qv, ql, u and v                                |
  !     * variances of qt, w, u, theta_l and theta_v                |
  !     * skewness of qt and w                                      |
  !*** *genstat*  calculates timeseries of several variables       |
  !                                                                 |
  !____________________SETTINGS_AND_SWITCHES________________________|
  !                     IN &NAMTIMESTAT                             |
  !                                                                 |
  !    dtav           SAMPLING INTERVAL                             |
  !                                                                 |
  !    timeav         INTERVAL OF WRITING                           |
  !                                                                 |
  !    lstat      SWITCH TO ENABLE TIMESERIES                       |
  !-----------------------------------------------------------------|
  use modglobal, only : longint

  implicit none
  ! private
  PUBLIC :: initgenstat, exitgenstat, average2homog, average1homog, interpolate
  save

  !NetCDF variables
  !integer :: nvar = 37
  integer :: ncid,nrec = 0
  character(80) :: fname = 'profiles.xxx.nc'
  character(80),allocatable, dimension(:,:) :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav, timeav,tnext,tnextwrite
  !  integer(kind=longint) :: idtav,itimeav,tnext,tnextwrite
  logical :: lstat= .false. ! switch for conditional sampling cloud (on/off)
  integer :: nsamples
  logical :: lcanopy= .false.    ! switch for computing canopy properties
  logical :: ltkebudget= .false. ! switch for computing tke budget equation
  integer :: istart         ! i-index where canopy starts
  integer :: ni             ! length of repeating unit i x-direction (until end of domain)
!     ----  total fields  ---


! J. Tomas:
    real, allocatable :: ucflux (:,:,:,:) ! turb u'c' flux
    real, allocatable :: wcflux (:,:,:,:) ! turb w'c' flux
    real, allocatable :: uwstress (:,:,:) ! Reynolds stress u'v'
    real, allocatable :: dissip   (:,:,:) ! Resolved dissipation u'v'*du/dz
    real, allocatable :: dissipsgs(:,:,:) ! Subgrid dissipation: nu_t * (du/dz)**2
    real, allocatable :: wthlflux (:,:,:) ! turbulent temp flux w'thl'
    real, allocatable :: tkeu (:,:,:)     ! tke u
    real, allocatable :: tkev (:,:,:)     ! tke v
    real, allocatable :: tkew (:,:,:)     ! tke w
    real, allocatable :: shearu (:,:,:)   ! shear stress in x-direction: nu*du/dz (on lower cell face)
    real, allocatable :: shearv (:,:,:)   ! shear stress in y-direction: nu*dv/dz (on lower cell face)

    real, allocatable :: slabu (:)     ! slab averaged u
    real, allocatable :: slabv (:)     ! slab averaged v
    real, allocatable :: slabw (:)     ! slab averaged w
    real, allocatable :: slabp (:)     ! slab averaged p
    real, allocatable :: slabthl (:)   ! slab averaged thl
    real, allocatable :: slabviscratio(:) ! slab averaged viscratio
    real, allocatable :: slabsv (:,:)     ! slab averaged sv
    real, allocatable :: slabucflux (:,:) ! slab averaged u'c'
    real, allocatable :: slabwcflux (:,:) ! slab averaged w'c'
    real, allocatable :: slabuwstress (:) ! slab averaged u'v'
    real, allocatable :: slabdissip   (:) ! slab averaged dissip u'v'*du/dz
    real, allocatable :: slabdissipsgs(:) ! slab averaged dissipsgs nu_t*(du/dz)**2
    real, allocatable :: slabwthlflux (:) ! slab averaged w'thl' 
    real, allocatable :: slabtkeu (:)     ! slab averaged tke u
    real, allocatable :: slabtkev (:)     ! slab averaged tke v
    real, allocatable :: slabtkew (:)     ! slab averaged tke w
    real, allocatable :: slabshearu (:)     ! slab averaged shear in x-direction (on lower cell face)
    real, allocatable :: slabshearv (:)     ! slab averaged shear in y direction (on lower cell face)
    real, allocatable :: ucfluxtot(:)       ! volume-averaged u'c'
    real, allocatable :: wcfluxtot(:)       ! volume-averaged w'c'
    real, allocatable :: ucfluxtotav(:)     ! volume-averaged and time-averaged u'c'
    real, allocatable :: wcfluxtotav(:)     ! volume-averaged and time-averaged w'c'

    real :: tketotu                       ! volume-averaged tke u
    real :: tketotv                       ! volume-averaged tke v
    real :: tketotw                       ! volume-averaged tke w
    real :: tketotuav                     ! volume-averaged and time-averaged tke u
    real :: tketotvav                     ! volume-averaged and time-averaged tke v
    real :: tketotwav                     ! volume-averaged and time-averaged tke w
    real :: timecompl                     ! completed time timeav (time interval for statistics)
    real :: timeaver                      ! completed part of statistics resolution


contains

  subroutine initgenstat
    use modmpi,    only : myid,mpierr, comm3d,my_real, mpi_logical, mpi_integer
    use modglobal, only : dtmax,kb,ke,kh,nsv,ifnamopt,fname_options, ifoutput, cexpnr,dtav_glob,timeav_glob,ladaptive,dt_lim,btime,ib,ie,ih,jb,je,jh,kb,ke,kh,totavtime
    use modsurfdata, only : isurf

    implicit none

    integer n, ierr
    character(40) :: name
    character(3) :: csvname
    namelist/NAMGENSTAT/ &
    dtav,timeav,lstat,lcanopy,ni,istart,ltkebudget

    !dtav=dtav_glob;timeav=timeav_glob

    if(myid==0)then
       open(ifnamopt,file=fname_options,status='old',iostat=ierr)
       read (ifnamopt,NAMGENSTAT,iostat=ierr)
       if (ierr > 0) then
          print *, 'Problem in namoptions NAMGENSTAT'
          print *, 'iostat error: ', ierr
          stop 'ERROR: Problem in namoptions NAMGENSTAT'
       endif
       write(6 ,NAMGENSTAT)
       close(ifnamopt)
    end if

    call MPI_BCAST(timeav     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dtav       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(lstat      ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lcanopy    ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(istart     ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(ni         ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(lstat   ,1,MPI_LOGICAL,0,comm3d,mpierr)
    !    idtav = dtav/tres
    !    itimeav = timeav/tres
    !    idtav = NINT(dtav/tres)           ! J. Tomas
    !    itimeav = NINT(timeav/tres)       ! J. Tomas

    tnext      = dtav   +btime
    tnextwrite = timeav +btime
    nsamples = NINT(timeav/dtav)
    if(.not.(lstat)) return
    !    dt_lim = min(dt_lim,tnext)
    if (abs(timeav/dtav-nsamples)>1e-4) then
       stop 'timeav must be a integer multiple of dtav'
    end if

    ! J. Tomas:
    allocate(slabu (kb:ke))
    allocate(slabv (kb:ke))
    allocate(slabw (kb:ke))
    allocate(slabp (kb:ke))
    allocate(slabthl(kb:ke))
    allocate(slabviscratio(kb:ke))
    allocate(slabsv (kb:ke,1:nsv))
    allocate(slabucflux (kb:ke,1:nsv))
    allocate(slabwcflux (kb:ke,1:nsv))
    allocate(slabuwstress (kb:ke))
    allocate(slabdissip (kb:ke))
    allocate(slabdissipsgs(kb:ke))
    allocate(slabwthlflux (kb:ke))
    allocate(slabtkeu (kb:ke))
    allocate(slabtkev (kb:ke))
    allocate(slabtkew (kb:ke))
    allocate(slabshearu (kb:ke))
    allocate(slabshearv (kb:ke))
    allocate(ucflux(ib-ih:ie+ih,jb-jh:je+jh,kb:ke,1:nsv))
    allocate(wcflux(ib-ih:ie+ih,jb-jh:je+jh,kb:ke,1:nsv))
    allocate(uwstress(ib-ih:ie+ih,jb-jh:je+jh,kb:ke))
    allocate(dissip(ib-ih:ie+ih,jb-jh:je+jh,kb:ke))
    allocate(dissipsgs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke))
    allocate(wthlflux(ib-ih:ie+ih,jb-jh:je+jh,kb:ke))
    allocate(tkeu (ib-ih:ie+ih,jb-jh:je+jh,kb:ke))
    allocate(tkev (ib-ih:ie+ih,jb-jh:je+jh,kb:ke))
    allocate(tkew (ib-ih:ie+ih,jb-jh:je+jh,kb:ke))
    allocate(shearu (ib-ih:ie+ih,jb-jh:je+jh,kb:ke))
    allocate(shearv (ib-ih:ie+ih,jb-jh:je+jh,kb:ke))
    allocate(ucfluxtot (1:nsv))
    allocate(wcfluxtot (1:nsv))
    allocate(ucfluxtotav (1:nsv))
    allocate(wcfluxtotav (1:nsv))


    slabu       = 0.
    slabv       = 0.
    slabw       = 0.
    slabp       = 0.
    slabthl     = 0.
    slabsv      = 0.
    slabucflux  = 0.
    slabwcflux  = 0.
    slabwthlflux = 0.
    slabuwstress= 0.
    slabdissip  = 0.
    slabdissipsgs= 0.
    slabviscratio = 0.
    slabtkeu    = 0.
    slabtkev    = 0.
    slabtkew    = 0.
    slabshearu  = 0.
    slabshearv  = 0.
    ucflux      = 0.
    wcflux      = 0.
    uwstress    = 0.
    dissip      = 0.
    dissipsgs   = 0.
    wthlflux    = 0.
    tkeu        = 0.
    tkev        = 0.
    tkew        = 0.
    shearu      = 0.
    shearv      = 0.
    ucfluxtot   = 0.
    wcfluxtot   = 0.
    ucfluxtotav = 0.
    wcfluxtotav = 0.
    timecompl   = totavtime
    timeaver    = 0.
    ! end J. Tomas

  end subroutine initgenstat

  subroutine exitgenstat
    use modmpi, only : myid
    implicit none

    if(.not.(lstat)) return
    deallocate(ucflux,wcflux,uwstress,wthlflux,tkeu,tkev,tkew,shearu,shearv,slabu, &
         slabv,slabw,slabp,slabthl,slabsv,slabucflux,slabwcflux,slabuwstress,& 
         slabwthlflux,slabtkeu,slabtkev,slabtkew,slabshearu,slabshearv,ucfluxtot,&
         wcfluxtot,ucfluxtotav,wcfluxtotav,slabviscratio,dissip,&
         dissipsgs,slabdissip,slabdissipsgs)


  end subroutine exitgenstat

  ! create time-averaged 3D fields for channel flow
  subroutine average2homog   
    use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh,jgb,jge,timee,btime,dt,rk3step,nsv,numol,numoli,&
         dzh,dzf,dzhi,dzhiq,dxf,dxhi,linoutflow,uinf,totavtime,linletgen,&
         dzfi,dxfi,dyi
    use modfields, only : um,vm,wm,umint,vmint,wmint,thlm,svm,uav,vav,wav,qtm,ql0,svm,uav,vav,wav,thlav,qtav,qlav,uuav,vvav,wwav,&
         uvav,uwav,vwav,svav,sv2av,thl2av,shear,friction, momthick,displthick, &
         upupav,vpvpav,wpwpav,upvpav,upwpav,vpwpav,pres0,viscratioav
    use modsubgriddata, only : ekm
    use modmpi, only : slabsum,myid,MY_REAL,comm3d,mpierr,MPI_SUM
    use modinletdata, only : upupavinl,vpvpavinl,wpwpavinl
    implicit none

    integer i,j,k,n,nhcells,snhcells,nk,kp,km,jp,jm,ip,im
    real timecompli

    !    real, allocatable :: slabu (:)     ! slab averaged u
    !    real, allocatable :: slabv (:)     ! slab averaged v
    !    real, allocatable :: slabw (:)     ! slab averaged w
    !    real, allocatable :: slabtkeu (:)     ! slab averaged tke u
    !    real, allocatable :: slabtkev (:)     ! slab averaged tke v
    !    real, allocatable :: slabtkew (:)     ! slab averaged tke w
    !    real, allocatable :: tkeu (:,:,:)     ! tke u
    !    real, allocatable :: tkev (:,:,:)     ! tke v
    !    real, allocatable :: tkew (:,:,:)     ! tke w
    real, allocatable :: slabwh (:)        ! slab averaged w at cell faces
    real, allocatable :: wmin (:,:,: )     ! w interpolated to cell center
    real, allocatable :: up (:,:,: )       ! u' at u-location
    real, allocatable :: vp (:,:,: )       ! v' at v-location
    real, allocatable :: wp (:,:,: )       ! w' at w-location
    real, allocatable :: uavj (:)          ! u averaged over j
    real, allocatable :: slabu_old (:)     ! old slabu
    real, allocatable :: slabv_old (:)     ! old slabv
    real, allocatable :: slabw_old (:)     ! old slabw
    real, allocatable :: slabp_old (:)     ! old slabp
    real, allocatable :: slabthl_old (:)   ! old slabthl
    real, allocatable :: slabviscratio_old (:)   ! old slabviscratio
    real, allocatable :: slabsv_old (:,:)  ! old slabsv
    real, allocatable :: slabucflux_old (:,:) ! old value
    real, allocatable :: slabwcflux_old (:,:) ! old value
    real, allocatable :: slabuwstress_old (:) ! old value
    real, allocatable :: slabdissip_old (:) ! old value
    real, allocatable :: slabdissipsgs_old (:) ! old value
    real, allocatable :: slabwthlflux_old (:)     ! old value
    real, allocatable :: slabtkeu_old (:)     ! old value
    real, allocatable :: slabtkev_old (:)     ! old value
    real, allocatable :: slabtkew_old (:)     ! old value
    real, allocatable :: slabshearu_old (:)   ! old value
    real, allocatable :: slabshearv_old (:)   ! old value

    real   :: myfriction(ib:ie)
    real   :: totfriction(ib:ie)
    real   :: mthick(kb:ke)
    real   :: dthick(kb:ke)
    real   :: mom(kb:ke)
    real   :: utopavj                        ! j-average top velocity at i-location  
    real   :: myutop                        ! sending utop (summed over jb:je 
    real   :: totutop                        ! receiving utop (summed over jgb:jge) 
    real   :: tketotuav_old                    ! old value of total averaged tke (u-component)
    real   :: tketotvav_old                    ! old value of total averaged tke (v-component)
    real   :: tketotwav_old                    ! old value of total averaged tke (w-component)
    real   :: ucfluxtotav_old(1:nsv)           ! old value of total averaged u'c'
    real   :: wcfluxtotav_old(1:nsv)           ! old value of total averaged u'c'
    real   :: wprime                           ! dummy variable
    real   :: eomp, eomm, emom, emop, dudz  

    if (lstat==.true. .AND. rk3step==3) then

       !    real :: tketotu
       !    real :: tketotv
       !    real :: tketotw

       !    allocate(slabu (kb-kh:ke+kh))
       !    allocate(slabv (kb-kh:ke+kh))
       !    allocate(slabw (kb-kh:ke+kh))
       !    allocate(slabtkeu (kb-kh:ke+kh))
       !    allocate(slabtkev (kb-kh:ke+kh))
       !    allocate(slabtkew (kb-kh:ke+kh))
       !    allocate(tkeu (ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
       !    allocate(tkev (ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
       !    allocate(tkew (ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
       allocate(wmin (ib-ih:ie+ih,jb:je,kb:ke))
       allocate(up (ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
       allocate(vp (ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
       allocate(wp (ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
       allocate(slabwh (kb:ke))
       allocate(uavj (kb:ke))
       allocate(slabu_old (kb:ke))
       allocate(slabv_old (kb:ke))
       allocate(slabw_old (kb:ke))
       allocate(slabp_old (kb:ke))
       allocate(slabthl_old (kb:ke))
       allocate(slabviscratio_old (kb:ke))
       allocate(slabsv_old (kb:ke,1:nsv))
       allocate(slabucflux_old (kb:ke,1:nsv))
       allocate(slabwcflux_old (kb:ke,1:nsv))
       allocate(slabuwstress_old (kb:ke))
       allocate(slabdissip_old (kb:ke))
       allocate(slabdissipsgs_old (kb:ke))
       allocate(slabwthlflux_old (kb:ke))
       allocate(slabtkeu_old (kb:ke))
       allocate(slabtkev_old (kb:ke))
       allocate(slabtkew_old (kb:ke))
       allocate(slabshearu_old (kb:ke))
       allocate(slabshearv_old (kb:ke))


       if (timecompl>timeav) then
          timecompl = dt
          call writestatistics2homog
       else
          timecompl = timecompl + dt
       end if

       if (timeaver> dtav) then

          timecompli = 1./timecompl
          slabwh = 0.

          ! remember old values
          slabu_old = slabu
          slabv_old = slabv
          slabw_old = slabw
          slabp_old = slabp
          slabthl_old = slabthl
          slabviscratio_old = slabviscratio
          slabsv_old  = slabsv
          slabucflux_old   = slabucflux
          slabwcflux_old   = slabwcflux
          slabuwstress_old = slabuwstress
          slabdissip_old   = slabdissip
          slabdissipsgs_old= slabdissipsgs
          slabwthlflux_old = slabwthlflux
          slabtkeu_old     = slabtkeu
          slabtkev_old     = slabtkev
          slabtkew_old     = slabtkew     
          slabshearu_old   = slabshearu
          slabshearv_old   = slabshearv
          tketotuav_old    = tketotuav
          tketotvav_old    = tketotvav
          tketotwav_old    = tketotwav
          ucfluxtotav_old  = ucfluxtotav
          wcfluxtotav_old  = wcfluxtotav
          ! slab-averaged velocity
          !    nhcells = (ie-ib+1+2*ih)*(jge-jgb+1+2*jh) ! no. of cells in horizontal plane: use global j-indices!!!
          nhcells = (ie-ib+1+2*ih)*(jge-jgb+1) ! no. of cells in horizontal plane: use global j-indices!!!
          snhcells = (ie-ib+1)*(jge-jgb+1)      ! no. of cells in horizontal small plane: use global j-indices!!!
          nk      = (ke-kb+1)                     ! no of cells in vertical direction  

          ! Interpolate wm to half cell height
          do k=kb,ke
             do j=jb,je
                do i=ib-ih,ie+ih
                   wmin(i,j,k)=0.5*(wm(i,j,k+1)+wm(i,j,k))
                end do
             end do
          end do

          call slabsum(slabu  ,kb,ke,um(:,:,kb:ke)   ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke)
          call slabsum(slabv  ,kb,ke,vm(:,:,kb:ke)   ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke)
          call slabsum(slabw  ,kb,ke,wmin            ,ib-ih,ie+ih,jb   ,je   ,kb,ke,ib-ih,ie+ih,jb,je,kb,ke)
          call slabsum(slabwh  ,kb,ke,wm(:,:,kb:ke)  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke)
          call slabsum(slabp  ,kb,ke,pres0(:,:,kb:ke),ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke)
          call slabsum(slabthl,kb,ke,thlm(:,:,kb:ke) ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke)
          call slabsum(slabviscratio,kb,ke,(ekm(:,:,kb:ke)-numol)*numoli,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke)
          do n=1,nsv
             call slabsum(slabsv(:,n),kb,ke,svm(ib-ih:ie+ih,jb-jh:je+jh,kb:ke,n) ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke)
          end do
          slabu   = slabu / nhcells
          slabv   = slabv / nhcells
          slabw   = slabw / nhcells
          slabwh  = slabwh / nhcells
          slabp   = slabp / nhcells
          slabthl = slabthl / nhcells
          slabviscratio = slabviscratio / nhcells
          slabsv  = slabsv  / nhcells

          ! Define fluctuating velocities
          do k = kb,ke
             do j = jb-jh,je+jh
                do i = ib-ih,ie+ih
                   up(i,j,k) = um(i,j,k) - slabu(k)
                   vp(i,j,k) = vm(i,j,k) - slabv(k)
                   wp(i,j,k) = wm(i,j,k) - slabwh(k)
                end do
             end do
          end do


          ! BC's
          up(:,:,kb-1)= -up(:,:,kb)   ! no slip
          up(:,:,ke+1)= -up(:,:,ke)   ! no slip
          vp(:,:,kb-1)= -vp(:,:,kb)   ! no slip
          vp(:,:,ke+1)= -vp(:,:,ke)   ! no slip
          wp(:,:,kb-1)= -wp(:,:,kb+1)
          wp(:,:,ke+1)= 0.

          ! Compute rate of strain tensor based on fluctuating components
          ! This will give SijSij, while dissipation = 2*nu*SijSij
          do k=kb,ke
             kp=k+1
             km=k-1
             do j=jb,je
                jp=j+1
                jm=j-1
                do i=ib,ie
                   ip=i+1
                   im=i-1
                   dissip(i,j,k) =  ( &
                        ((up(ip,j,k)-up(i,j,k))    *dxfi(i)        )**2    + &
                        ((vp(i,jp,k)-vp(i,j,k))    *dyi        )**2    + &
                        ((wp(i,j,kp)-wp(i,j,k))    *dzfi(k)     )**2    )

                   dissip(i,j,k) = dissip(i,j,k) + 0.125 * ( &
                        ((wp(i,j,kp)-wp(im,j,kp))   *dxhi(i)     + &
                        (up(i,j,kp)-up(i,j,k))      *dzhi(kp)  )**2    + &
                        ((wp(i,j,k)-wp(im,j,k))     *dxhi(i)     + &
                        (up(i,j,k)-up(i,j,km))      *dzhi(k)   )**2    + &
                        ((wp(ip,j,k)-wp(i,j,k))     *dxhi(ip)     + &
                        (up(ip,j,k)-up(ip,j,km))    *dzhi(k)   )**2    + &
                        ((wp(ip,j,kp)-wp(i,j,kp))   *dxhi(ip)     + &
                        (up(ip,j,kp)-up(ip,j,k))    *dzhi(kp)  )**2    )

                   dissip(i,j,k) = dissip(i,j,k) + 0.125 * ( &
                        ((up(i,jp,k)-up(i,j,k))     *dyi     + &
                        (vp(i,jp,k)-vp(im,jp,k))    *dxhi(i)        )**2    + &
                        ((up(i,j,k)-up(i,jm,k))     *dyi     + &
                        (vp(i,j,k)-vp(im,j,k))      *dxhi(i)        )**2    + &
                        ((up(ip,j,k)-up(ip,jm,k))   *dyi     + &
                        (vp(ip,j,k)-vp(i,j,k))      *dxhi(ip)       )**2    + &
                        ((up(ip,jp,k)-up(ip,j,k))   *dyi     + &
                        (vp(ip,jp,k)-vp(i,jp,k))    *dxhi(ip)       )**2    )

                   dissip(i,j,k)= dissip(i,j,k) + 0.125 * ( &
                        ((vp(i,j,kp)-vp(i,j,k))    *dzhi(kp) + &
                        (wp(i,j,kp)-wp(i,jm,kp))   *dyi        )**2    + &
                        ((vp(i,j,k)-vp(i,j,km))    *dzhi(k)+ &
                        (wp(i,j,k)-wp(i,jm,k))     *dyi        )**2    + &
                        ((vp(i,jp,k)-vp(i,jp,km))  *dzhi(k)+ &
                        (wp(i,jp,k)-wp(i,j,k))     *dyi        )**2    + &
                        ((vp(i,jp,kp)-vp(i,jp,k))  *dzhi(kp) + &
                        (wp(i,jp,kp)-wp(i,j,kp))   *dyi        )**2    )
                end do
             end do
          end do
          dissip(:,:,:) =2.*numol* dissip(:,:,:)  ! D = 2*nu*s_ij*s_ij
          ! add BC's because average is over ib-ih:ie+ih (just extrapolate)
          dissip(ib-ih,:,:) = dissip(ib,:,:)   
          dissip(ie+ih,:,:) = dissip(ie,:,:)   

          ! slab-averaged TKE        
          do k=kb,ke
             do j=jb,je   
                do i=ib-ih,ie
                   wprime          = (wmin(i,j,k)-slabw(k))
                   uwstress(i,j,k) = (   0.5*(um(i,j,k)+um(i+1,j,k)) -slabu(k)) * wprime     ! -u'w' (both u' and w' are interpolated to the c.c.
                   ! dissip approximately equals u'w' * du/dz

                   dudz            =   0.25*( (um(i,j,k+1)-um(i,j,k))*dzhi(k+1)  + (um(i,j,k)-um(i,j,k-1))*dzhi(k) + &            ! interpolate du/dz
                        (um(i+1,j,k+1)-um(i+1,j,k))*dzhi(k+1)  + (um(i+1,j,k)-um(i+1,j,k-1))*dzhi(k) )
                   !          dissip(i,j,k)   = uwstress(i,j,k) * dudz 
                   dissipsgs(i,j,k)= (ekm(i,j,k) - numol)* dudz**2
                   wthlflux(i,j,k) = wprime * (thlm(i,j,k) -slabthl(k))                      ! w'thl' 
                   do n=1,nsv
                      ucflux(i,j,k,n) = (   0.5*(um(i,j,k)+um(i+1,j,k)) -slabu(k)) * &
                           (svm(i,j,k,n)-slabsv(k,n))                      ! u'c' (u' is interpolated to the c.c.)
                      wcflux(i,j,k,n) = wprime * (svm(i,j,k,n)-slabsv(k,n))                   ! w'c' (w' is interpolated to the c.c.)
                   end do
                end do
                i=ie+ih        ! for that last cell
                wprime          = (wmin(i,j,k)-slabw(k))
                uwstress(i,j,k) = (   0.5*(um(i,j,k)+um(ib-ih,j,k))-slabu(k)) * wprime    ! -u'w' (both u' and w' are interpolated to the c.c.
                wthlflux(i,j,k) = wprime * (thlm(i,j,k) -slabthl(k))                      ! w'thl'  
                do n=1,nsv
                   ucflux(i,j,k,n) = (   0.5*(um(i,j,k)+um(ib-ih,j,k)) -slabu(k)) * &
                        (svm(i,j,k,n)-slabsv(k,n))                      ! u'c' (u' is interpolated to the c.c.)
                   wcflux(i,j,k,n) = wprime * (svm(i,j,k,n)-slabsv(k,n))                   ! w'c' (w' is interpolated to the c.c.)
                end do
             end do
          end do

          do k=kb,ke  
             do j=jb,je
                do i=ib-ih,ie+ih
                   tkeu(i,j,k)     = (um(i,j,k)-slabu(k))**2
                   tkev(i,j,k)     = (vm(i,j,k)-slabv(k))**2 
                   tkew(i,j,k)     = (wmin(i,j,k)-slabw(k))**2 
                end do
             end do
          end do

          !    call slabsum(slabtkeu  ,kb,ke,tkeu  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib,ie,jb,je,kb,ke) 
          !    call slabsum(slabtkev  ,kb,ke,tkev  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib,ie,jb,je,kb,ke) 
          !    call slabsum(slabtkew  ,kb,ke,tkew  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib,ie,jb,je,kb,ke) 
          call slabsum(slabuwstress,kb,ke,uwstress  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke) 
          call slabsum(slabdissip,kb,ke,dissip  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke) 
          call slabsum(slabdissipsgs,kb,ke,dissipsgs,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke) 
          call slabsum(slabwthlflux,kb,ke,wthlflux  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke) 
          call slabsum(slabtkeu  ,kb,ke,tkeu  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke) 
          call slabsum(slabtkev  ,kb,ke,tkev  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke) 
          call slabsum(slabtkew  ,kb,ke,tkew  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke) 
          do n=1,nsv
             call slabsum(slabucflux(:,n),kb,ke,ucflux(:,:,:,n),ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke) 
             call slabsum(slabwcflux(:,n),kb,ke,wcflux(:,:,:,n),ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke) 
          end do
          slabuwstress = slabuwstress / nhcells
          slabdissip   = slabdissip / nhcells
          slabdissipsgs= slabdissipsgs / nhcells
          slabwthlflux = slabwthlflux / nhcells
          slabtkeu     = slabtkeu / nhcells
          slabtkev     = slabtkev / nhcells
          slabtkew     = slabtkew / nhcells
          slabucflux  = slabucflux / nhcells
          slabwcflux  = slabwcflux / nhcells

          ! volume-averaged TKE
          tketotu=sum(slabtkeu(:))/nk
          tketotv=sum(slabtkev(:))/nk
          tketotw=sum(slabtkew(:))/nk
          do n=1,nsv
             ucfluxtot(n) = sum(slabucflux(:,n))/nk   
             wcfluxtot(n) = sum(slabwcflux(:,n))/nk   
          end do

          ! Computing ustar 
          do k=kb,ke
             kp = k+1
             km = k-1   
             do j=jb,je
                jp = j+1
                jm = j-1
                !        do i=ib-ih,ie+ih
                do i=ib,ie
                   emop = ( dzf(kp) * ( ekm(i,j,k)*dxf(i-1)  + ekm(i-1,j,k)*dxf(i) )  + &              ! dx is non-equidistant
                        dzf(k)  * ( ekm(i,j,kp)*dxf(i-1) + ekm(i-1,j,kp)*dxf(i) ) )*dxhi(i) * dzhiq(kp)
                   emom = ( dzf(km) * ( ekm(i,j,k)*dxf(i-1)  + ekm(i-1,j,k)*dxf(i) )  + &             ! dx is non-equidistant
                        dzf(k)  * ( ekm(i,j,km)*dxf(i-1) + ekm(i-1,j,km)*dxf(i) ) )*dxhi(i) * dzhiq(k)
                   eomp = ( dzf(kp) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                        dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
                   eomm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                        dzf(k)  * ( ekm(i,j,km) + ekm(i,jm,km) ) ) * dzhiq(k)

                   shearu(i,j,k)=0.5*( emop*(um(i,j,k+1)-um(i,j,k))*dzhi(k+1)  + emom*(um(i,j,k)-um(i,j,k-1))*dzhi(k))  ! average of upper and lower u-shear stress
                   shearv(i,j,k)=0.5*( eomp*(vm(i,j,k+1)-vm(i,j,k))*dzhi(k+1)  + eomm*(vm(i,j,k)-vm(i,j,k-1))*dzhi(k))  ! average of upper and lower v-shear stress
                end do
             end do
          end do
          !    call slabsum(slabshearu  ,kb,ke,shearu  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib,ie,jb,je,kb,ke)
          !    call slabsum(slabshearv  ,kb,ke,shearv  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib,ie,jb,je,kb,ke)
          call slabsum(slabshearu  ,kb,ke,shearu  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib,ie,jb,je,kb,ke)
          call slabsum(slabshearv  ,kb,ke,shearv  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib,ie,jb,je,kb,ke)
          slabshearu=slabshearu/snhcells
          slabshearv=slabshearv/snhcells

          ! compute the new time averaged values  
          !    timecompl  = runtime - timeleft
          !    write(6,*) 'totavtime=', totavtime
          !    timecompl  = totavtime + timee - btime   ! completed time in this simulation + totavtime (from previous sim)

          tketotuav=(tketotuav_old*(timecompl-timeaver) + tketotu*timeaver)*timecompli
          tketotvav=(tketotvav_old*(timecompl-timeaver) + tketotv*timeaver)*timecompli
          tketotwav=(tketotwav_old*(timecompl-timeaver) + tketotw*timeaver)*timecompli
          ucfluxtotav = (ucfluxtotav_old*(timecompl-timeaver) + ucfluxtot*timeaver)*timecompli
          wcfluxtotav = (wcfluxtotav_old*(timecompl-timeaver) + wcfluxtot*timeaver)*timecompli

          do k=kb,ke
             slabu(k)=(slabu_old(k)*(timecompl-timeaver) + slabu(k)*timeaver)*timecompli
             slabv(k)=(slabv_old(k)*(timecompl-timeaver) + slabv(k)*timeaver)*timecompli
             slabw(k)=(slabw_old(k)*(timecompl-timeaver) + slabw(k)*timeaver)*timecompli
             slabp(k)=(slabp_old(k)*(timecompl-timeaver) + slabp(k)*timeaver)*timecompli
             slabthl(k)=(slabthl_old(k)*(timecompl-timeaver) + slabthl(k)*timeaver)*timecompli
             slabviscratio(k)=(slabviscratio_old(k)*(timecompl-timeaver) + slabviscratio(k)*timeaver)*timecompli
             slabuwstress(k)=(slabuwstress_old(k)*(timecompl-timeaver) + slabuwstress(k)*timeaver)*timecompli
             slabdissip(k)=(slabdissip_old(k)*(timecompl-timeaver) + slabdissip(k)*timeaver)*timecompli
             slabdissipsgs(k)=(slabdissipsgs_old(k)*(timecompl-timeaver) + slabdissipsgs(k)*timeaver)*timecompli
             slabwthlflux(k)=(slabwthlflux_old(k)*(timecompl-timeaver) + slabwthlflux(k)*timeaver)*timecompli
             slabtkeu(k)=(slabtkeu_old(k)*(timecompl-timeaver) + slabtkeu(k)*timeaver)*timecompli
             slabtkev(k)=(slabtkev_old(k)*(timecompl-timeaver) + slabtkev(k)*timeaver)*timecompli
             slabtkew(k)=(slabtkew_old(k)*(timecompl-timeaver) + slabtkew(k)*timeaver)*timecompli
             slabshearu(k)=(slabshearu_old(k)*(timecompl-timeaver) + slabshearu(k)*timeaver)*timecompli
             slabshearv(k)=(slabshearv_old(k)*(timecompl-timeaver) + slabshearv(k)*timeaver)*timecompli
             ! scalar pollutants
             slabsv(k,:)=(slabsv_old(k,:)*(timecompl-timeaver) + slabsv(k,:)*timeaver)*timecompli
             slabucflux(k,:) = (slabucflux_old(k,:)*(timecompl-timeaver) + slabucflux(k,:)*timeaver)*timecompli
             slabwcflux(k,:) = (slabwcflux_old(k,:)*(timecompl-timeaver) + slabwcflux(k,:)*timeaver)*timecompli
          enddo

          do k=kb,ke  
             do j=jb,je
                do i=ib,ie
                   uav(i,j,k)=(uav(i,j,k)*(timecompl-timeaver)  + umint(i,j,k)*timeaver)*timecompli
                   vav(i,j,k)=(vav(i,j,k)*(timecompl-timeaver)  + vmint(i,j,k)*timeaver)*timecompli
                   wav(i,j,k)=(wav(i,j,k)*(timecompl-timeaver)  + wmint(i,j,k)*timeaver)*timecompli
                   thlav(i,j,k)=(thlav(i,j,k)*(timecompl-timeaver)  + thlm(i,j,k)*timeaver)*timecompli
                   qtav(i,j,k)=(qtav(i,j,k)*(timecompl-timeaver)  + qtm(i,j,k)*timeaver)*timecompli
                   qlav(i,j,k)=(qlav(i,j,k)*(timecompl-timeaver)  + ql0(i,j,k)*timeaver)*timecompli  !no qlm since ql is calculated from difference of qt and qsat
                   ! overline((u')^2) = overline(u^2) - uav^2
                   uuav(i,j,k)=(uuav(i,j,k)*(timecompl-timeaver) + (umint(i,j,k)**2)*timeaver)*timecompli    
                   vvav(i,j,k)=(vvav(i,j,k)*(timecompl-timeaver) + (vmint(i,j,k)**2)*timeaver)*timecompli
                   wwav(i,j,k)=(wwav(i,j,k)*(timecompl-timeaver) + (wmint(i,j,k)**2)*timeaver)*timecompli
                   uvav(i,j,k)=(uvav(i,j,k)*(timecompl-timeaver) + umint(i,j,k)*vmint(i,j,k)*timeaver)*timecompli
                   vwav(i,j,k)=(vwav(i,j,k)*(timecompl-timeaver) + vmint(i,j,k)*wmint(i,j,k)*timeaver)*timecompli
                   uwav(i,j,k)=(uwav(i,j,k)*(timecompl-timeaver) + umint(i,j,k)*wmint(i,j,k)*timeaver)*timecompli
                   thl2av(i,j,k)=(thl2av(i,j,k)*(timecompl-timeaver)  + (thlm(i,j,k)**2)*timeaver)*timecompli 
                   viscratioav(i,j,k)=(viscratioav(i,j,k)*(timecompl-timeaver)  + (ekm(i,j,k)-numol)*numoli*timeaver)*timecompli
                   do n=1,nsv   
                   svav(i,j,k,n)=(svav(i,j,k,n)*(timecompl-timeaver)  + svm(i,j,k,n)*timeaver)*timecompli
                   sv2av(i,j,k,n)=(sv2av(i,j,k,n)*(timecompl-timeaver)  + (svm(i,j,k,n)**2)*timeaver)*timecompli 
                      !            svav(i,j,k,n)=(svav(i,j,k,n)*(ntrun-1)  + svm(i,j,k,n))*ntruni
                      !            sv2av(i,j,k,n)=(sv2av(i,j,k,n)*(ntrun-1)  + svm(i,j,k,n)**2)*ntruni - svav(i,j,k,n)**2
                   enddo
                enddo
             enddo
          enddo

          upupav = uuav - uav**2             ! overline(u'u') = overline(uu) - U^2
          vpvpav = vvav - vav**2             ! overline(v'v') = overline(vv) - V^2
          wpwpav = wwav - wav**2             ! overline(w'w') = overline(ww) - W^2
          upvpav = uvav - uav*vav     ! overline(u'v') = overline(uv) - U*V
          upwpav = uwav - uav*wav     ! overline(u'w') = overline(uw) - U*W
          vpwpav = vwav - vav*wav     ! overline(v'w') = overline(vw) - V*W

          !     do k=kb,ke
          !     do j=jb,je
          !     do i=ib,ie
          !       upupav(i,j,k) = uuav(i,j,k) - uav(i,j,k)*uav(i,j,k)
          !       vpvpav(i,j,k) = vvav(i,j,k) - vav(i,j,k)*vav(i,j,k)
          !       wpwpav(i,j,k) = wwav(i,j,k) - wav(i,j,k)*wav(i,j,k)
          !       upvpav(i,j,k) = uvav(i,j,k) - uav(i,j,k)*vav(i,j,k)
          !       upwpav(i,j,k) = uwav(i,j,k) - uav(i,j,k)*wav(i,j,k)
          !       vpwpav(i,j,k) = vwav(i,j,k) - vav(i,j,k)*wav(i,j,k)
          !     end do
          !     end do
          !     end do


          if (linletgen == 1) then
             upupavinl=0. 
             call slabsum(upupavinl  ,kb,ke,upupav(ib:ib,jb:je,kb:ke),ib,ib,jb,je,kb,ke,ib,ib,jb,je,kb,ke)
             vpvpavinl=0. 
             call slabsum(vpvpavinl  ,kb,ke,vpvpav(ib:ib,jb:je,kb:ke),ib,ib,jb,je,kb,ke,ib,ib,jb,je,kb,ke)
             wpwpavinl=0. 
             call slabsum(wpwpavinl  ,kb,ke,wpwpav(ib:ib,jb:je,kb:ke),ib,ib,jb,je,kb,ke,ib,ib,jb,je,kb,ke)
             upupavinl = upupavinl /(jge-jgb+1)
             vpvpavinl = vpvpavinl /(jge-jgb+1)
             wpwpavinl = wpwpavinl /(jge-jgb+1)
          end if


          !! compute boundary layer properties
          !!   if (linoutflow==.true.) then
          !! Determine average skin friction
          !      do i = ib,ie
          !        myfriction(i) = sum(shear(i,jb:je,kb,3))
          !      end do
          !      call MPI_ALLREDUCE(myfriction, totfriction, ie-ib+1,  MY_REAL, &
          !                          MPI_SUM, comm3d,mpierr)
          !      friction = totfriction/(jge-jgb+1)                  ! Now we have the average shear at each i-location
          !      friction = 2.*friction/((uinf)**2)    ! compute skin friction Cf = tau/(0.5*rho*u_inf^2)  (shear is already tau/rho!)
          !! End of determining skin-friction
          !
          !! determine momentum thickness 
          !     do i=ib,ie
          !       call slabsum(uavj  ,kb,ke,um(i:i,jb:je,kb:ke) ,i,i,jb,je,kb,ke,i,i,jb,je,kb,ke) ! determine horizontal (j) average velocity (variable in k)
          !!       call slabsum(uavj  ,kb,ke,uav(i:i,jb:je,kb:ke) ,i,i,jb,je,kb,ke,i,i,jb,je,kb,ke) ! determine horizontal (j) average velocity (variable in k) use average flow field
          !       uavj = uavj / (jge-jgb+1)
          !!       myutop = sum(um(i,jb:je,ke))
          !!       myutop = sum(uav(i,:,ke))
          !!       call MPI_ALLREDUCE(myutop, totutop, 1,  MY_REAL, &
          !!                          MPI_SUM, comm3d,mpierr)
          !!       utopavj = totutop / (jge-jgb+1)
          !       utopavj = um(ib,jb,ke)  ! inflow is plug-flow 
          !!       momthick(i) = sum(   ((uavj(i,:)/uinf) - (uavj(i,:)**2)/(uinf**2))*dzf)  ! momentum thickness
          !!       momthick(i) = sum( ((uavj(:)/uinf) - (uavj(:)**2 / uinf**2) )*dzf)  ! momentum thickness 
          !       do k=kb,ke
          !         mthick(k) = ((uavj(k)/utopavj) - (uavj(k)**2 / utopavj**2) )*dzf(k)
          !         dthick(k) = (1. - (uavj(k) / utopavj) )*dzf(k)
          !!         mthick(k) = ((uavj(k)/0.226) - (uavj(k)**2 / 0.226**2) )*dzf(k)
          !!         dthick(k) = (1. - (uavj(k) / 0.226) )*dzf(k)
          !       end do
          !       momthick(i)   = sum(mthick)  ! momentum thickness 
          !       displthick(i)  = sum(dthick)  ! displacement thickness 
          !!       momthick(i)   = sum( ((uavj(:)/utopavj) - (uavj(:)**2 / utopavj**2) )*dzf)  ! momentum thickness 
          !!       displthick(i) = sum( (1. - (uavj(:) / utopavj) )*dzf)  ! displacement thickness 
          !     end do 
          !!    end if ! linoutflow

          timeaver = dt 
       else ! timeaver < dtav
          timeaver = timeaver + dt
       end if

       deallocate(wmin,uavj,slabu_old,slabv_old,slabw_old,slabp_old,slabthl_old,slabsv_old,&
            slabucflux_old,slabwcflux_old,slabuwstress_old,slabwthlflux_old,slabtkeu_old,&
            slabtkev_old,slabtkew_old,slabshearu_old,slabshearv_old,slabviscratio_old,&
            slabdissip_old,slabdissipsgs_old,up,vp,wp,slabwh)

    end if ! lstat==true .and. rk3step ==3
  end subroutine average2homog

  ! create time-averaged 3D fields for spatially-developing flow (1 homogeneous direction)
  subroutine average1homog   
    use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh,jgb,jge,timee,btime,dt,rk3step,nsv,numol,numoli,&
         dzh,dzf,dzhi,dzhiq,dxf,dxhi,linoutflow,uinf,totavtime,linletgen,&
         dzfi,dxfi,dyi,startmean,dy2i,dxfiq,dxhiq,dyiq,dzfi5
    use modfields, only : u0,v0,w0,thl0,qt0,ql0,sv0,uav,vav,wav,thlav,qtav,qlav,uuav,vvav,wwav,&
         uvav,uwav,vwav,svav,sv2av,thl2av,thluav,thlvav,thlwav,svuav,svvav,svwav, &
         pres0,presav,viscratioav,uusgsav,vvsgsav,wwsgsav,uwsgsav,thlusgsav,&
         thlwsgsav,svusgsav,svwsgsav,tkesgsav,strain2av,disssgsav,nusgsav,&
         tvmx,tvmy,tvmz,tpm,ttmx,ttmy,ttmz,tsgsmx1,tsgsmx2,tsgsmy1,tsgsmy2,&
         tsgsmz1,tsgsmz2
    use modsubgriddata, only : ekm,prandtli
    use modmpi, only : slabsum,myid,MY_REAL,comm3d,mpierr,MPI_SUM
    implicit none
 
    real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)  :: tekm  ! turbulent viscosity 
    !real, dimension(ib:ie+ih,jb:je+jh,kb:ke+kh)           :: tkeres  ! instantaneous TKE (resolved)
 
    integer i,j,k,n,nhcells,snhcells,nk,kpp,km,jpp,jm,ipp,im,ip,jp,kp
    real timecompli,uw,uv,vw,thlu,thlv,thlw,svu,svv,svw,&
         tkesgs,nusgs,diffsgs,uusgs,vvsgs,wwsgs,uwsgs,&
         thlusgs,thlwsgs,strain2,svusgs,svwsgs,dummy,&
         emom,eomm,eopm,epom,emmo,eomp,epmo,emop,empo
 
    if (lstat==.true. .AND. rk3step==3) then

    if (timecompl>timeav) then
      timecompl = dt
!      call writestatistics1homog
    elseif (timee-btime < startmean) then  ! only start taking averages after startmean has elapsed
!      timecompl = dt
      timecompl = timeaver
!      write(6,*) 'elapsed time is smaller than startmean ' 
    else
      timecompl = timecompl + dt
    end if
  
    if (timeaver> dtav) then
!    write(6,*) 'timecompl = ', timecompl
    timecompli = 1./timecompl
   
    ! remember old values
    nhcells = (ie-ib+1+2*ih)*(jge-jgb+1) ! no. of cells in horizontal plane: use global j-indices!!!
    nk      = (ke-kb+1)                     ! no of cells in vertical direction  

! turbulent viscosity 
    tekm(:,:,:) = ekm(:,:,:) - numol

! Compute stresses and fluxes at c.c.
    do k = kb,ke
      kp=k+1
      km=k-1
      do j = jb,je
        jp=j+1
        jm=j-1
        do i = ib,ie
          ip=i+1
          im=i-1
          
          strain2 =  ( &
            ((u0(ip,j,k)-u0(i,j,k))    *dxfi(i)     )**2    + &
            ((v0(i,jp,k)-v0(i,j,k))    *dyi         )**2    + &
            ((w0(i,j,kp)-w0(i,j,k))    *dzfi(k)     )**2    )

          strain2 = strain2 + 0.125 * ( &
            ((w0(i,j,kp)-w0(im,j,kp))   *dxhi(i)     + &
            (u0(i,j,kp)-u0(i,j,k))      *dzhi(kp)  )**2    + &
            ((w0(i,j,k)-w0(im,j,k))     *dxhi(i)     + &
            (u0(i,j,k)-u0(i,j,km))      *dzhi(k)   )**2    + &
            ((w0(ip,j,k)-w0(i,j,k))     *dxhi(ip)     + &
            (u0(ip,j,k)-u0(ip,j,km))    *dzhi(k)   )**2    + &
            ((w0(ip,j,kp)-w0(i,j,kp))   *dxhi(ip)     + &
            (u0(ip,j,kp)-u0(ip,j,k))    *dzhi(kp)  )**2    )

          strain2 = strain2 + 0.125 * ( &
            ((u0(i,jp,k)-u0(i,j,k))     *dyi     + &
            (v0(i,jp,k)-v0(im,jp,k))    *dxhi(i)        )**2    + &
            ((u0(i,j,k)-u0(i,jm,k))     *dyi     + &
            (v0(i,j,k)-v0(im,j,k))      *dxhi(i)        )**2    + &
            ((u0(ip,j,k)-u0(ip,jm,k))   *dyi     + &
            (v0(ip,j,k)-v0(i,j,k))      *dxhi(ip)       )**2    + &
            ((u0(ip,jp,k)-u0(ip,j,k))   *dyi     + &
            (v0(ip,jp,k)-v0(i,jp,k))    *dxhi(ip)       )**2    )

          strain2 = strain2 + 0.125 * ( &
            ((v0(i,j,kp)-v0(i,j,k))    *dzhi(kp) + &
            (w0(i,j,kp)-w0(i,jm,kp))   *dyi        )**2    + &
            ((v0(i,j,k)-v0(i,j,km))    *dzhi(k)+ &
            (w0(i,j,k)-w0(i,jm,k))     *dyi        )**2    + &
            ((v0(i,jp,k)-v0(i,jp,km))  *dzhi(k)+ &
            (w0(i,jp,k)-w0(i,j,k))     *dyi        )**2    + &
            ((v0(i,jp,kp)-v0(i,jp,k))  *dzhi(kp) + &
            (w0(i,jp,kp)-w0(i,j,kp))   *dyi        )**2    )

!          snorm(i,j,k) = sqrt(2*strain2)
          tkesgs         = 2.*tekm(i,j,k)*sqrt(2.*strain2) 
          strain2av(i,j,k)=(strain2av(i,j,k)*(timecompl-timeaver) + strain2*timeaver)*timecompli   ! update average strain2av
          tkesgsav(i,j,k) =(tkesgsav(i,j,k) *(timecompl-timeaver) + tkesgs*timeaver) *timecompli   ! update average tkesgsav
          disssgsav(i,j,k)=(disssgsav(i,j,k)*(timecompl-timeaver) + 2.*tekm(i,j,k)*strain2*timeaver)*timecompli              ! = 2*<nu_t*SijSij> 
          nusgsav(i,j,k) = (nusgsav(i,j,k)*(timecompl-timeaver) + tekm(i,j,k)*timeaver)*timecompli   ! update average subgrid viscosity
        end do
      end do
    end do

!    ! BC's for snorm  
!    snorm(ib-1,:,:) = snorm(ib,:,:)
!    snorm(ie+1,:,:) = snorm(ie,:,:)
!    snorm(:,:,kb-1) = snorm(:,:,kb)   ! not correct...
!    snorm(:,:,ke+1) = snorm(:,:,ke)
!    snorm(:,jb-1,:) = snorm(:,je,:)
!    snorm(:,je+1,:) = snorm(:,jb,:)


! Averages are at cell faces(u,v,w) and centers (temperature, viscosity)
! One ghost cell is saved
    do k=kb-1,ke+1  
      do j=jb-1,je+1
        do i=ib-1,ie+1
          uav(i,j,k)=(uav(i,j,k)*(timecompl-timeaver)      + u0(i,j,k)*timeaver)*timecompli
          vav(i,j,k)=(vav(i,j,k)*(timecompl-timeaver)      + v0(i,j,k)*timeaver)*timecompli
          wav(i,j,k)=(wav(i,j,k)*(timecompl-timeaver)      + w0(i,j,k)*timeaver)*timecompli
          thlav(i,j,k)=(thlav(i,j,k)*(timecompl-timeaver)  + thl0(i,j,k)*timeaver)*timecompli
          qtav(i,j,k)=(qtav(i,j,k)*(timecompl-timeaver)    + qt0(i,j,k)*timeaver)*timecompli
          qlav(i,j,k)=(qlav(i,j,k)*(timecompl-timeaver)    + ql0(i,j,k)*timeaver)*timecompli
          presav(i,j,k)=(presav(i,j,k)*(timecompl-timeaver)+ pres0(i,j,k)*timeaver)*timecompli
          uuav(i,j,k)=(uuav(i,j,k)*(timecompl-timeaver)    + (u0(i,j,k)**2)*timeaver)*timecompli    
          vvav(i,j,k)=(vvav(i,j,k)*(timecompl-timeaver)    + (v0(i,j,k)**2)*timeaver)*timecompli
          wwav(i,j,k)=(wwav(i,j,k)*(timecompl-timeaver)    + (w0(i,j,k)**2)*timeaver)*timecompli
          thl2av(i,j,k)=(thl2av(i,j,k)*(timecompl-timeaver)+ (thl0(i,j,k)**2)*timeaver)*timecompli 
          viscratioav(i,j,k)=(viscratioav(i,j,k)*(timecompl-timeaver)  + (ekm(i,j,k)-numol)*numoli*timeaver)*timecompli       
          do n=1,nsv
            svav(i,j,k,n)=(svav(i,j,k,n)*(timecompl-timeaver)  + sv0(i,j,k,n)*timeaver)*timecompli
            sv2av(i,j,k,n)=(sv2av(i,j,k,n)*(timecompl-timeaver)  + (sv0(i,j,k,n)**2)*timeaver)*timecompli 
          end do
        enddo
      enddo
    enddo

! Compute stresses and fluxes

! SGS normal stresses:
    do k=kb-1,ke+1
      km = k-1
      do j=jb-1,je-1
!        do i=ib,ie+1
        do i=ib,ie
          im = i-1
          nusgs = 0.5*(tekm(i,j,k)*dxf(i-1)+tekm(i-1,j,k)*dxf(i))*dxhi(i)                    ! turb. visc at uu locations
          uusgs = - 2.*nusgs*0.5*((u0(i+1,j,k)-u0(i,j,k))*dxfi(i)*dxf(i-1) + &
                                  (u0(i,j,k)-u0(i-1,j,k))*dxfi(i-1)*dxf(i))*dxhi(i) 
          uusgsav(i,j,k)=(uusgsav(i,j,k)*(timecompl-timeaver) + uusgs*timeaver)*timecompli   ! update average uusgsav
        end do
      end do
    end do

    do k=kb-1,ke+1
      km = k-1
      do j=jb,je
        do i=ib,ie+1
          im = i-1
          nusgs = 0.5*(tekm(i,j,k)+tekm(i,j-1,k))                                 ! turb. visc at vv locations
          vvsgs = - 2.*nusgs *0.5*((v0(i,j+1,k)-v0(i,j,k))*dyi + &
                                   (v0(i,j,k)-v0(i,j-1,k))*dyi)
          vvsgsav(i,j,k)=(vvsgsav(i,j,k)*(timecompl-timeaver) + vvsgs*timeaver)*timecompli   ! update average vvsgsav

        end do
      end do
    end do

!    do k=kb,ke+1
    do k=kb,ke
      km = k-1
      do j=jb-1,je-1
        do i=ib-1,ie+1
          im = i-1
          nusgs = 0.5*(tekm(i,j,k)*dzf(k-1)+tekm(i,j,k-1)*dzf(k))*dzhi(k)         ! turb. visc at ww locations
          wwsgs = - 2.*nusgs*0.5*((w0(i,j,k+1)-w0(i,j,k))*dzfi(k)*dzf(k-1) +    &
                                  (w0(i,j,k)-w0(i,j,k-1))*dzfi(k-1)*dzf(k))*dzhi(k)
          wwsgsav(i,j,k)=(wwsgsav(i,j,k)*(timecompl-timeaver) + wwsgs*timeaver)*timecompli   ! update average ww
        end do
      end do
    end do

! Resolved and normal stresses:
! uw: ib:ie+1 jb:je kb:ke+1
    do k=kb,ke+1 
      km = k-1 
      do j=jb,je
        do i=ib,ie+1
          im = i-1
          uw         = 0.25*(u0(i,j,km)*dzf(k) + u0(i,j,k)*dzf(km))*dzhi(k) * &     ! interpolate u to edges
                            (w0(im,j,k)*dxf(i) + w0(i,j,k)*dxf(im))*dxhi(i)         ! interpolate w to edges
          uwav(i,j,k)=(uwav(i,j,k)*(timecompl-timeaver) + uw*timeaver)*timecompli   ! update average vw
! subgrid:
          nusgs = 0.5*(0.5*(tekm(i,j,k)*dxf(i-1) + tekm(i-1,j,k)*dxf(i))*dxhi(i)*dzf(k-1) + &
                       0.5*(tekm(i,j,k-1)*dxf(i-1) + tekm(i-1,j,k-1)*dxf(i))*dxhi(i)*dzf(k))*dzhi(k)

          uwsgs = - nusgs*( (u0(i,j,k)-u0(i,j,k-1))*dzhi(k) + &
                            (w0(i,j,k)-w0(i-1,j,k))*dxhi(i))

          uwsgsav(i,j,k)=(uwsgsav(i,j,k)*(timecompl-timeaver) + uwsgs*timeaver)*timecompli   ! update average uwsgsav
        end do
      end do
    end do

! uv: ib:ie+1 jb:je kb:ke+1
    do k=kb,ke  
      do j=jb,je+1
        jm = j-1
        do i=ib,ie+1
          im = i-1
          uv         = 0.25*(u0(i,jm,k)  + u0(i,j,k)) * &                           ! interpolate u to edges
                            (v0(im,j,k)*dxf(i) + v0(i,j,k)*dxf(im))*dxhi(i)         ! interpolate v to edges
          uvav(i,j,k)=(uvav(i,j,k)*(timecompl-timeaver) + uv*timeaver)*timecompli   ! update average uv
        end do
      end do
    end do
 

! vw: ib:ie jb:je+1 kb:ke+1
    do k=kb,ke+1
      km = k-1  
      do j=jb,je+1
        jm = j-1
        do i=ib,ie
          vw         = 0.25*(w0(i,jm,k)  + w0(i,j,k)) * &                           ! interpolate w to edges
                            (v0(i,j,km)*dzf(k) + v0(i,j,k)*dzf(km))*dzhi(k)         ! interpolate v to edges
          vwav(i,j,k)=(vwav(i,j,k)*(timecompl-timeaver) + vw*timeaver)*timecompli   ! update average vw
        end do
      end do
    end do



! thlu and svu: ib:ie+1 jb:je kb:ke  (located on u-faces)
    do k=kb,ke  
      do j=jb,je
        do i=ib,ie+1
          im = i-1
          thlu         = 0.5* u0(i,j,k)  * &                                              ! no interpolation
                            (thl0(im,j,k)*dxf(i) + thl0(i,j,k)*dxf(im))*dxhi(i)           ! interpolate thl to u-faces
          thluav(i,j,k)=(thluav(i,j,k)*(timecompl-timeaver) + thlu*timeaver)*timecompli   ! update average thlu
! SGS
          diffsgs  = prandtli*0.5*(tekm(i,j,k)*dxf(i-1)+tekm(i-1,j,k)*dxf(i))*dxhi(i)
          thlusgs = - diffsgs * (thl0(i,j,k) - thl0(i-1,j,k))*dxhi(i)
          thlusgsav(i,j,k)=(thlusgsav(i,j,k)*(timecompl-timeaver) + thlusgs*timeaver)*timecompli   ! update average thlusgs

          do n=1,nsv
            svu         = 0.5* u0(i,j,k)  * &                                              ! no interpolation
                              (sv0(im,j,k,n)*dxf(i) + sv0(i,j,k,n)*dxf(im))*dxhi(i)        ! interpolate sv to u-faces
            svuav(i,j,k,n)=(svuav(i,j,k,n)*(timecompl-timeaver) + svu*timeaver)*timecompli   ! update average svu
   ! SGS
            svusgs = - diffsgs * (sv0(i,j,k,n) - sv0(i-1,j,k,n))*dxhi(i)
            svusgsav(i,j,k,n)=(svusgsav(i,j,k,n)*(timecompl-timeaver) + svusgs*timeaver)*timecompli   ! update average svusgs
          end do
        end do
      end do
    end do

! thlv and svv: ib:ie jb:je+1 kb:ke  (located on v-faces)
    do k=kb,ke 
      do j=jb,je+1
      jm = j-1
        do i=ib,ie
          thlv         = 0.5* v0(i,j,k) * &                                               ! no interpolation
                           (thl0(i,jm,k) + thl0(i,j,k))                                   ! interpolate thl to v-faces
          thlvav(i,j,k)=(thlvav(i,j,k)*(timecompl-timeaver) + thlv*timeaver)*timecompli   ! update average thlv
          do n=1,nsv
            svv         = 0.5* v0(i,j,k) * &                                              ! no interpolation
                             (sv0(i,jm,k,n) + sv0(i,j,k,n))                               ! interpolate sv to v-faces
            svvav(i,j,k,n)=(svvav(i,j,k,n)*(timecompl-timeaver) + svv*timeaver)*timecompli    ! update average svv
          end do
        end do
      end do
    end do

! thlw and svw: ib:ie jb:je kb:ke+1  (located on w-faces)
    do k=kb,ke+1  
      km = k-1
      do j=jb,je
        do i=ib,ie
          thlw         = 0.5 * w0(i,j,k) * &                                              ! no interpolation
                            (thl0(i,j,km)*dzf(k) + thl0(i,j,k)*dzf(km))*dzhi(k)           ! interpolate thl to w-faces
          thlwav(i,j,k)=(thlwav(i,j,k)*(timecompl-timeaver) + thlw*timeaver)*timecompli   ! update average thlw
! SGS flux
          diffsgs  = prandtli*0.5*(tekm(i,j,k)*dzf(k-1)+tekm(i,j,k-1)*dzf(k))*dzhi(k)
          thlwsgs = - diffsgs * (thl0(i,j,k) - thl0(i,j,k-1))*dzhi(k)
          thlwsgsav(i,j,k)=(thlwsgsav(i,j,k)*(timecompl-timeaver) + thlwsgs*timeaver)*timecompli   ! update average thlw
          do n=1,nsv
            svw         = 0.5 * w0(i,j,k) * &                                              ! no interpolation
                              (sv0(i,j,km,n)*dzf(k) + sv0(i,j,k,n)*dzf(km))*dzhi(k)        ! interpolate sv to w-faces
            svwav(i,j,k,n)=(svwav(i,j,k,n)*(timecompl-timeaver) + svw*timeaver)*timecompli ! update average svw
  ! SGS flux
            svwsgs = - diffsgs * (sv0(i,j,k,n) - sv0(i,j,k-1,n))*dzhi(k)
            svwsgsav(i,j,k,n)=(svwsgsav(i,j,k,n)*(timecompl-timeaver) + svwsgs*timeaver)*timecompli ! update average svw
          end do
        end do
      end do
    end do

    if (ltkebudget==.true.) then
!      do k=kb-kh,ke
!      kp =k+1
!      do j=jb-jh,je
!      jp =j+1
!      do i=ib-ih,ie
!        ip =i+1
!        tkeres(i,j,k) = 0.5*(0.5*(uuav(i,j,k)+ uuav(ip,j,k)) + &
!                             0.5*(vvav(i,j,k)+ vvav(i,jp,k)) + &
!                             0.5*(wwav(i,j,k)+ wwav(i,j,kp)))
!      end do
!      end do
!      end do

!    write(6,*) 'START TKE BUDGET COMPUTATION'

      do k=kb,ke
      kp =k+1
      km =k-1
      do j=jb,je
      jp =j+1
      jm =j-1
      do i=ib,ie
      ip =i+1
      im =i-1
      ! Transport by pressure fluctuations (mean: -<u_j*d/dxj(p)>, needed for -<u_j'*d/dxj(p')> =  -<u_j*d/dxj(p)> + <u_j>*d/dxj(<p>)
        dummy      = -  0.5*(u0(ip,j,k)*(pres0(ip,j,k) - pres0(i,j,k)) *dxhi(ip)   + &
                             u0(i,j,k) *(pres0(i,j,k)  - pres0(im,j,k))*dxhi(i))     &
                     -  0.5*(v0(i,jp,k)*(pres0(i,jp,k) - pres0(i,j,k)) *dyi        + &
                             v0(i,j,k) *(pres0(i,j,k)  - pres0(i,jm,k))*dyi)         &
                     -  0.5*(w0(i,j,kp)*(pres0(i,j,kp) - pres0(i,j,k)) *dzhi(kp)   + &
                             w0(i,j,k) *(pres0(i,j,k)  - pres0(i,j,km))*dzhi(k))
        tpm(i,j,k)=(tpm(i,j,k)*(timecompl-timeaver) + dummy*timeaver)*timecompli   ! update average uwsgsav
      end do
      end do
      end do


      do k=kb,ke
      kp =k+1
!      kpp=k+2
      km =k-1
      do j=jb,je
      jp =j+1
!      jpp=j+2
      jm =j-1
      do i=ib,ie
      ip =i+1
!      ipp=i+2
      im =i-1
 

!       <ui*d/dxj(ui*uj)> equals the advection routine times ui
          dummy  =      &
                u0(i,j,k)*(( &
                (u0(i,j,k) + u0(ip,j,k)) *( u0(i,j,k)+u0(ip,j,k) ) &
               -(u0(i,j,k) + u0(im,j,k)) *( u0(i,j,k)+u0(im,j,k) ) &     ! d(uu)/dx
                )*dxhiq(i) &
                  +(  &
                (u0(i,j,k)+u0(i,jp,k))*(v0(i,jp,k)+v0(im,jp ,k)) &
               -(u0(i,j,k)+u0(i,jm,k))*(v0(i,j  ,k)+v0(im,j  ,k)) &       ! d(vu)/dy
                )*dyiq  & 
                 +  &
                ( &
                ( u0(i,j,kp)*dzf(k) + u0(i,j,k)*dzf(kp) ) * dzhi(kp) &
                  *( w0(i,j,kp)+ w0(im,j,kp) ) &
               -( u0(i,j,k)*dzf(km) + u0(i,j,km)*dzf(k) ) * dzhi(k) &
                  *( w0(i,j,k)  + w0(im,j,k)   ) &
                )*0.5*dzfi5(k))        

        ttmx(i,j,k)=(ttmx(i,j,k)*(timecompl-timeaver) + dummy*timeaver)*timecompli   ! update 

        dummy =      & 
              v0(i,j,k)*(( &
              ( u0(ip,j,k)+u0(ip,jm,k)) &
                *(v0(i,j,k)*dxf(ip) + v0(ip,j,k)*dxf(i) ) * dxhi(ip) &
              -(u0(i ,j,k)+u0(i ,jm,k)) &
                *(v0(i,j,k)*dxf(im) + v0(im,j,k)*dxf(i) ) * dxhi(i)  &             ! d(uv)/dx
              )*dxfiq(i) &
              +( &
              ( v0(i,jp,k)+v0(i,j,k))*(v0(i,j,k)+v0(i,jp,k)) &
              -(v0(i,jm,k)+v0(i,j,k))*(v0(i,j,k)+v0(i,jm,k)) &               ! d(vv)/dy
              )*dyiq  +                                                     &
             ( &
            (w0(i,j,kp)+w0(i,jm,kp)) &
            *(v0(i,j,kp)*dzf(k)+v0(i,j,k)*dzf(kp) ) * dzhi(kp) &
            -(w0(i,j,k)+w0(i,jm,k)) &
            *(v0(i,j,km)*dzf(k)+v0(i,j,k)*dzf(km)) * dzhi(k) &
            ) * 0.5*dzfi5(k))                                              
 
        ttmy(i,j,k)=(ttmy(i,j,k)*(timecompl-timeaver) + dummy*timeaver)*timecompli   ! update 

        dummy    =   & 
               w0(i,j,k)*( ( &
                ( w0(ip,j,k)*dxf(i) + w0(i,j,k)*dxf(ip) )*dxhi(ip) &         ! d(uw)/dx
              *( dzf(km)*u0(ip,j,k) + dzf(k)*u0(ip,j,km) ) &
              -( w0(i,j,k)*dxf(im) + w0(im,j,k)*dxf(i) ) *dxhi(i)  &
              *( dzf(km)*u0(i,j,k)+dzf(k)*u0(i ,j,km) ) &
                )*dxfiq(i) * dzhi(k) &
              + &
                ( &
                ( w0(i,jp,k) + w0(i,j,k) ) &                                ! d(vw)/dy
              *( dzf(km)*v0(i,jp,k) + dzf(k)*v0(i,jp,km) ) &
              -( w0(i,j,k) + w0(i,j-1,k) ) &
              *( dzf(km)*v0(i,j,k) + dzf(k)*v0(i,j,km) ) &
                ) *dyiq * dzhi(k) &
              + &
                ( &
                ( w0(i,j,k)+w0(i,j,kp) ) * (w0(i,j,k) + w0(i,j,kp) ) &      ! d(ww)/dz 
              -( w0(i,j,k)+w0(i,j,km) ) * (w0(i,j,k) + w0(i,j,km) ) &
                ) * dzhiq(k) )                                              

        ttmz(i,j,k)=(ttmz(i,j,k)*(timecompl-timeaver) + dummy*timeaver)*timecompli   ! update 

 
! Tvmx at u-locations (ib:ih+ih:jb:je,kb:ke)
! This is similar to routine diffu time u_i
        dummy =  u0(i,j,k)*(                           &
                  ( numol  * (u0(i+1,j,k)-u0(i,j,k))*dxfi(i) &
                   -numol * (u0(i,j,k)-u0(i-1,j,k))*dxfi(i-1) ) * 2. * dxhi(i) &
                  + &
                  ( numol * ( (u0(i,jp,k)-u0(i,j,k))   *dyi &
                            +(v0(i,jp,k)-v0(i-1,jp,k))*dxhi(i)) &
                   -numol * ( (u0(i,j,k)-u0(i,jm,k))   *dyi &
                            +(v0(i,j,k)-v0(i-1,j,k))  *dxhi(i)) &
                                       ) * dyi &
                  + &
                  ( numol * ( (u0(i,j,kp)-u0(i,j,k))   *dzhi(kp) &
                            +(w0(i,j,kp)-w0(i-1,j,kp))*dxhi(i)) &
                   -numol * ( (u0(i,j,k)-u0(i,j,km))   *dzhi(k) &
                            +(w0(i,j,k)-w0(i-1,j,k))  *dxhi(i)) &
                                       ) *dzfi(k)  )
        tvmx(i,j,k)=(tvmx(i,j,k)*(timecompl-timeaver) + dummy*timeaver)*timecompli   ! update average tvmx

        emom = ( dzf(km) * ( tekm(i,j,k)*dxf(i-1)  + tekm(i-1,j,k)*dxf(i) )  + &             ! dx is non-equidistant
                 dzf(k)  * ( tekm(i,j,km)*dxf(i-1) + tekm(i-1,j,km)*dxf(i) ) )*dxhi(i) * dzhiq(k)
        emop = ( dzf(kp) * ( tekm(i,j,k)*dxf(i-1)  + tekm(i-1,j,k)*dxf(i) )  + &              ! dx is non-equidistant
                 dzf(k)  * ( tekm(i,j,kp)*dxf(i-1) + tekm(i-1,j,kp)*dxf(i) ) )*dxhi(i) * dzhiq(kp)
        empo = 0.25 * ( ( tekm(i,j,k)+tekm(i,jp,k))*dxf(i-1) + (tekm(i-1,j,k)+tekm(i-1,jp,k))*dxf(i) ) * dxhi(i)  ! dx is non-equidistant
        emmo = 0.25 * ( ( tekm(i,j,k)+tekm(i,jm,k))*dxf(i-1)  +(tekm(i-1,jm,k)+tekm(i-1,j,k))*dxf(i) ) * dxhi(i)  ! dx is non-equidistant

!        dummy =  u0(i,j,k)*(                           &
        dummy =  (                           &
                   ( tekm(i,j,k)  * (u0(i+1,j,k)-u0(i,j,k))*dxfi(i) &
                    -tekm(i-1,j,k)* (u0(i,j,k)-u0(i-1,j,k))*dxfi(i-1) ) * 2. * dxhi(i) &
                  + &
                  ( empo * ( (u0(i,jp,k)-u0(i,j,k))   *dyi &
                            +(v0(i,jp,k)-v0(i-1,jp,k))*dxhi(i)) &
                    -emmo * ( (u0(i,j,k)-u0(i,jm,k))   *dyi &
                            +(v0(i,j,k)-v0(i-1,j,k))  *dxhi(i)) &
                                       ) * dyi &
                  + &
                  ( emop * ( (u0(i,j,kp)-u0(i,j,k))   *dzhi(kp) &
                            +(w0(i,j,kp)-w0(i-1,j,kp))*dxhi(i)) &
                    -emom * ( (u0(i,j,k)-u0(i,j,km))   *dzhi(k) &
                            +(w0(i,j,k)-w0(i-1,j,k))  *dxhi(i)) &
                                       ) *dzfi(k) )
        tsgsmx1(i,j,k)=(tsgsmx1(i,j,k)*(timecompl-timeaver) + dummy*u0(i,j,k)*timeaver)*timecompli   ! update average tsgsmx1
        tsgsmx2(i,j,k)=(tsgsmx2(i,j,k)*(timecompl-timeaver) + dummy*timeaver)*timecompli             ! update average tsgsmx2


! Tvmv at v-locations (ib:ih:jb:je+1,kb:ke)
! This is similar to routine diffv time v
       dummy = v0(i,j,k) * (                            &
              ( numol * ( (v0(i+1,j,k)-v0(i,j,k))   *dxhi(i+1) &
                        +(u0(i+1,j,k)-u0(i+1,jm,k))*dyi) &
                -numol * ( (v0(i,j,k)-v0(i-1,j,k))   *dxhi(i) &
                        +(u0(i,j,k)-u0(i,jm,k))    *dyi) &
                           ) * dxfi(i) &        ! = d/dx( Km*(dv/dx + du/dy) )
                + &
              (numol * (v0(i,jp,k)-v0(i,j,k)) &
              -numol * (v0(i,j,k)-v0(i,jm,k))  ) * 2. * dy2i &        ! = d/dy( 2*Km*(dv/dy) )
                + &
              ( numol * ( (v0(i,j,kp)-v0(i,j,k))    *dzhi(kp) &
                        +(w0(i,j,kp)-w0(i,jm,kp))  *dyi) &
                -numol * ( (v0(i,j,k)-v0(i,j,km))    *dzhi(k) &
                        +(w0(i,j,k)-w0(i,jm,k))    *dyi)   &
                           ) * dzfi(k) )      ! = d/dz( Km*(dv/dz + dw/dy) )
        tvmy(i,j,k)=(tvmy(i,j,k)*(timecompl-timeaver) + dummy*timeaver)*timecompli   ! update average uwsgsav

          eomm = ( dzf(km) * ( tekm(i,j,k)  + tekm(i,jm,k)  )  + &
                      dzf(k)  * ( tekm(i,j,km) + tekm(i,jm,km) ) ) * dzhiq(k)
          eomp = ( dzf(kp) * ( tekm(i,j,k)  + tekm(i,jm,k)  )  + &
                      dzf(k)  * ( tekm(i,j,kp) + tekm(i,jm,kp) ) ) * dzhiq(kp)
          emmo = 0.25 * ( ( tekm(i,j,k)+tekm(i,jm,k))*dxf(i-1)  +(tekm(i-1,jm,k)+tekm(i-1,j,k))*dxf(i) ) * dxhi(i)  ! dx is non-equidistant
          epmo = 0.25 * ( ( tekm(i,j,k)+tekm(i,jm,k))*dxf(i+1) + (tekm(i+1,jm,k)+tekm(i+1,j,k))*dxf(i) ) * dxhi(i+1)  ! dx is non-equidistant

!       dummy = v0(i,j,k) * (                            &
       dummy = (                            &
               ( epmo * ( (v0(i+1,j,k)-v0(i,j,k))   *dxhi(i+1) &
                        +(u0(i+1,j,k)-u0(i+1,jm,k))*dyi) &
                -emmo * ( (v0(i,j,k)-v0(i-1,j,k))   *dxhi(i) &
                        +(u0(i,j,k)-u0(i,jm,k))    *dyi) &
                           ) * dxfi(i) &        ! = d/dx( Km*(dv/dx + du/dy) )
                + &
              (tekm(i,j,k) * (v0(i,jp,k)-v0(i,j,k)) &
              -tekm(i,jm,k)* (v0(i,j,k)-v0(i,jm,k))  ) * 2. * dy2i &        ! = d/dy( 2*Km*(dv/dy) )
                + &
              ( eomp * ( (v0(i,j,kp)-v0(i,j,k))    *dzhi(kp) &
                        +(w0(i,j,kp)-w0(i,jm,kp))  *dyi) &
                -eomm * ( (v0(i,j,k)-v0(i,j,km))    *dzhi(k) &
                        +(w0(i,j,k)-w0(i,jm,k))    *dyi)   &
                           ) * dzfi(k)  )     ! = d/dz( Km*(dv/dz + dw/dy) )
        tsgsmy1(i,j,k)=(tsgsmy1(i,j,k)*(timecompl-timeaver) + dummy*v0(i,j,k)*timeaver)*timecompli   ! update average tsgsmy1  = <v*d/dxj(2*nu*S2j)>
        tsgsmy2(i,j,k)=(tsgsmy2(i,j,k)*(timecompl-timeaver) + dummy*timeaver)*timecompli             ! update average tsgsmy2  = <d/dxj(2*nu*S2j)>

! Tvmz at w-locations (ib:ih:jb:je,kb:ke+kh)
 ! This is similar to routine diffw time w
        dummy = w0(i,j,k) * (                                         &
                          ( numol * ( (w0(i+1,j,k)-w0(i,j,k))    *dxhi(i+1) &
                            +(u0(i+1,j,k)-u0(i+1,j,km)) *dzhi(k) ) &
                    -numol * ( (w0(i,j,k)-w0(i-1,j,k))    *dxhi(i) &
                            +(u0(i,j,k)-u0(i,j,km))     *dzhi(k) ) &
                             )*dxfi(i) &
                + &
                  ( numol * ( (w0(i,jp,k)-w0(i,j,k))     *dyi &
                            +(v0(i,jp,k)-v0(i,jp,km))   *dzhi(k) ) &
                    -numol * ( (w0(i,j,k)-w0(i,jm,k))     *dyi &
                            +(v0(i,j,k)-v0(i,j,km))     *dzhi(k) ) &
                             )*dyi &
                + &
                  ( numol * (w0(i,j,kp)-w0(i,j,k)) *dzfi(k) &
                   -numol * (w0(i,j,k)-w0(i,j,km)) *dzfi(km) ) * 2. &
                                                              * dzhi(k))
        tvmz(i,j,k)=(tvmz(i,j,k)*(timecompl-timeaver) + dummy*timeaver)*timecompli   ! update average uwsgsav

          emom = ( dzf(km) * ( tekm(i,j,k)*dxf(i-1)  + tekm(i-1,j,k)*dxf(i) )*dxhi(i)  + &
                      dzf(k)  * ( tekm(i,j,km)*dxf(i-1) + tekm(i-1,j,km)*dxf(i) )*dxhi(i) ) * dzhiq(k)
          eomm = ( dzf(km) * ( tekm(i,j,k)  + tekm(i,jm,k)  )  + &
                      dzf(k)  * ( tekm(i,j,km) + tekm(i,jm,km) ) ) * dzhiq(k)
          eopm = ( dzf(km) * ( tekm(i,j,k)  + tekm(i,jp,k)  )  + &
                      dzf(k)  * ( tekm(i,j,km) + tekm(i,jp,km) ) ) * dzhiq(k)
          epom = ( dzf(km) * ( tekm(i,j,k)*dxf(i+1)  + tekm(i+1,j,k)*dxf(i) )*dxhi(i+1)  + &
                      dzf(k)  * ( tekm(i,j,km)*dxf(i+1) + tekm(i+1,j,km)*dxf(i) )*dxhi(i+1) ) * dzhiq(k)

!        dummy = w0(i,j,k) * (                                         &
        dummy =   (                                         &
                  ( epom * ( (w0(i+1,j,k)-w0(i,j,k))    *dxhi(i+1) &
                            +(u0(i+1,j,k)-u0(i+1,j,km)) *dzhi(k) ) &
                    -emom * ( (w0(i,j,k)-w0(i-1,j,k))    *dxhi(i) &
                            +(u0(i,j,k)-u0(i,j,km))     *dzhi(k) ) &
                             )*dxfi(i) &
                + &
                  ( eopm * ( (w0(i,jp,k)-w0(i,j,k))     *dyi &
                            +(v0(i,jp,k)-v0(i,jp,km))   *dzhi(k) ) &
                    -eomm * ( (w0(i,j,k)-w0(i,jm,k))     *dyi &
                            +(v0(i,j,k)-v0(i,j,km))     *dzhi(k) ) &
                             )*dyi &
                + &
                  ( tekm(i,j,k) * (w0(i,j,kp)-w0(i,j,k)) *dzfi(k) &
                  -tekm(i,j,km)* (w0(i,j,k)-w0(i,j,km)) *dzfi(km) ) * 2. &
                                                              * dzhi(k))
        tsgsmz1(i,j,k)=(tsgsmz1(i,j,k)*(timecompl-timeaver) + dummy*w0(i,j,k)*timeaver)*timecompli   ! update average tsgsmz1 = <w*d/dxj(2*nu*S3j)>
        tsgsmz2(i,j,k)=(tsgsmz2(i,j,k)*(timecompl-timeaver) + dummy*timeaver)*timecompli            ! update average tsgsmz2 = <d/dxj(2*nu*S3j)>
      end do
      end do
      end do


    end if



    call slabsum(slabu  ,kb,ke,u0(:,:,kb:ke)   ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke)
    call slabsum(slabv  ,kb,ke,v0(:,:,kb:ke)   ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke)
    call slabsum(slabw  ,kb,ke,w0(:,:,kb:ke)  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke)
    slabu   = slabu / nhcells
    slabv   = slabv / nhcells
    slabw   = slabw / nhcells

    
    do k=kb,ke  
      do j=jb,je
        do i=ib-ih,ie+ih
          tkeu(i,j,k)     = (u0(i,j,k)-slabu(k))**2
          tkev(i,j,k)     = (v0(i,j,k)-slabv(k))**2 
          tkew(i,j,k)     = (w0(i,j,k)-slabw(k))**2 
        end do
      end do
    end do

    call slabsum(slabtkeu  ,kb,ke,tkeu  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke) 
    call slabsum(slabtkev  ,kb,ke,tkev  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke) 
    call slabsum(slabtkew  ,kb,ke,tkew  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke,ib-ih,ie+ih,jb,je,kb,ke) 
    slabtkeu     = slabtkeu / nhcells
    slabtkev     = slabtkev / nhcells
    slabtkew     = slabtkew / nhcells

   ! volume-averaged TKE
   tketotu=sum(slabtkeu(:))/nk
   tketotv=sum(slabtkev(:))/nk
   tketotw=sum(slabtkew(:))/nk



          timeaver = dt 
       else ! timeaver < dtav
          timeaver = timeaver + dt
       end if

    end if ! lstat==true .and. rk3step ==3
  end subroutine average1homog

  subroutine writestatistics2homog
    use modglobal, only : ib,ie,jb,je,kb,ke,ifoutput,cexpnr,zf,numol,grav,nsv
    use modsurfdata, only: thls,thl_top
    use modmpi,    only : myid
    integer n,k    
    character(25) filename, scalno
   
    if (myid==0) then
      write(6,*) 'Writing vertical profiles to statistics.//cexpnr'
      write(6,*) '!!! thls,thl_top=',thls,thl_top
      open (ifoutput,file='statistics.'//cexpnr,position='append')
      write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
      '#--------------------------------------------------------'      &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '    
      write (ifoutput,'(A/2A)') &
          '#--------------------------------------------------------' &
          , '#LEV     HGHT          U          V          W          T-Ts          P          nu*dU/dz         ' &
          , ' nu*dV/dz        -<uw>        nu*dU/dz - <uw>         <wT>         Ri_flux       TKEu          ' &
          , ' TKEv          TKEw        TKETOT_U        TKETOT_V          TKETOT_W ,     viscratio        ' &
          , ' <2*nu*s_ij*s_ij>    <nu_t*(du/dz)^2>' 

      do k=kb,ke
        write(ifoutput,'(I3,21e14.6)') &
            k, &
            zf           (k), &
            slabu        (k), &
            slabv        (k), &
            slabw        (k), &
            slabthl(k)-thls , &
            slabp        (k), &
            slabshearu   (k), &
            slabshearv   (k), &
           -slabuwstress (k), &
            slabshearu(k)-slabuwstress(k), &
            slabwthlflux (k), &
            (grav*slabwthlflux(k)*numol/thls)/(slabuwstress(k)*slabshearu(k)) , &
            slabtkeu     (k), &
            slabtkev     (k), &
            slabtkew     (k), &
            tketotuav       , &
            tketotvav       , &
            tketotwav       , &
            slabviscratio(k), &
            slabdissip(k),    &
            slabdissipsgs(k)
 
      end do
      close (ifoutput)

      do n=1,nsv
      write(scalno,'(i2.2)') n
      filename = 'statistics_scal'//trim(scalno)//'.'//cexpnr
      write(6,*) 'Writing vertical profiles to statistics_scalar.//cexpnr'
      open (ifoutput,file=filename,position='append')
      write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
      '#--------------------------------------------------------'      &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '
      write (ifoutput,'(A/2A)') &
          '#--------------------------------------------------------' &
          , '#LEV     HGHT          C          u''c''          w''c''       (u''c'')tot      (w''c'')tot '   

      do k=kb,ke
!        write(ifoutput,'(I3,F8.2,10F12.6)') &
        write(ifoutput,'(I3,6e14.6)') &
            k, &
            zf           (k)  , &
            slabsv       (k,n), &
            slabucflux   (k,n), &
            slabwcflux   (k,n), &
            ucfluxtotav  (n)  , &
            wcfluxtotav  (n) 
      end do ! loop over k
      close (ifoutput)
      end do ! loop over n scalars

    end if 
  end subroutine writestatistics2homog
  
! create 3D velocity fields in cell-centers
  subroutine interpolate   
    use modglobal, only : ib,ie,jb,je,kb,ke,rk3step
    use modfields, only : um,vm,wm,umint,vmint,wmint
    implicit none

    integer i,j,k
    if (rk3step==3) then
    do k=kb,ke
    do j=jb,je
    do i=ib,ie
       umint(i,j,k) = 0.5*(um(i,j,k)+um(i+1,j,k) )
       vmint(i,j,k) = 0.5*(vm(i,j,k)+vm(i,j+1,k) )
       wmint(i,j,k) = 0.5*(wm(i,j,k)+wm(i,j,k+1) )
    enddo
    enddo
    enddo
    endif
  end subroutine interpolate
end module modgenstat
