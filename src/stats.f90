module stats
  use modglobal,  only : cexpnr, ltdump, ltempeq, lmoist, lchem, nsv, rk3step, &
                         ib, ie, jb, je, kb, ke, kh, &
                         dxf, dzf, dzfi, dxhi, dzhi, dzh2i, dyi, dzhiq, &
                         timee, tstatsdump, tsample, dt, &
                         k1, JNO2
  use modfields,  only : um, vm, wm, pres0, thlm, qtm, svm, IIc
  use modsubgrid, only : ekh, ekm
  use modmpi,     only : cmyidx, cmyidy, myid
  use modstat_nc, only : ncinfo, open_nc, define_nc, writestat_dims_nc, writestat_nc
  implicit none
  public :: stats_init, stats_main, stats_exit
            ! stats_allocate_vel, stats_allocate_temp, stats_allocate_moist, stats_allocate_scalar, stats_allocate_PSS
            ! stats_compute_vel,  stats_compute_temp,  stats_compute_moist,  stats_compute_scalar , stats_compute_PSS
            ! stats_write_vel,    stats_write_temp,    stats_write_moist,    stats_write_scalar   , stats_write_PSS
            ! stats_compute_tavg, stats_interpolate_k, stats_compute_sgs
  save

  integer :: xdim, ydim, zdim
  real    :: tsamplep, tstatsdumpp, tstatsdumppi
  
  character(80)              :: filenamet
  character(80), allocatable :: tVars(:,:)
  integer                    :: ctrt, dumpcount, ncidt, nrect

  !! Variables to store time averaged quantities
  real, allocatable :: ut(:,:,:)
  real, allocatable :: vt(:,:,:)
  real, allocatable :: wt(:,:,:)
  real, allocatable :: pt(:,:,:)

  real, allocatable :: uik(:,:,:)
  real, allocatable :: wik(:,:,:)
  real, allocatable :: vjk(:,:,:)
  real, allocatable :: wjk(:,:,:)
  real, allocatable :: uij(:,:,:)
  real, allocatable :: vij(:,:,:)
  real, allocatable :: utik(:,:,:)
  real, allocatable :: wtik(:,:,:)
  real, allocatable :: uwtik(:,:,:)
  real, allocatable :: vtjk(:,:,:)
  real, allocatable :: wtjk(:,:,:)
  real, allocatable :: vwtjk(:,:,:)
  real, allocatable :: utij(:,:,:)
  real, allocatable :: vtij(:,:,:)
  real, allocatable :: uvtij(:,:,:)

  real, allocatable :: uc(:,:,:)
  real, allocatable :: vc(:,:,:)
  real, allocatable :: wc(:,:,:)
  real, allocatable :: utc(:,:,:)
  real, allocatable :: vtc(:,:,:)
  real, allocatable :: wtc(:,:,:)
  real, allocatable :: uutc(:,:,:)
  real, allocatable :: vvtc(:,:,:)
  real, allocatable :: wwtc(:,:,:)

  real, allocatable :: usgs(:,:,:)
  real, allocatable :: vsgs(:,:,:)
  real, allocatable :: wsgs(:,:,:)
  real, allocatable :: usgst(:,:,:)
  real, allocatable :: vsgst(:,:,:)
  real, allocatable :: wsgst(:,:,:)

  real, allocatable :: thlt(:,:,:)
  real, allocatable :: thlk(:,:,:)
  real, allocatable :: thltk(:,:,:)
  real, allocatable :: wthltk(:,:,:)
  real, allocatable :: thlthlt(:,:,:)
  real, allocatable :: thlsgs(:,:,:)
  real, allocatable :: thlsgst(:,:,:)

  real, allocatable :: qtt(:,:,:)
  real, allocatable :: qtk(:,:,:)
  real, allocatable :: qttk(:,:,:)
  real, allocatable :: wqttk(:,:,:)
  real, allocatable :: qtqtt(:,:,:)
  real, allocatable :: qtsgs(:,:,:)
  real, allocatable :: qtsgst(:,:,:)

  character(10), allocatable :: svtname(:)
  character(20), allocatable :: wpsvptname(:)
  character(20), allocatable :: svpsvptname(:)
  character(10), allocatable :: svsgsname(:)
  real, allocatable :: svt(:,:,:,:)
  real, allocatable :: svk(:,:,:,:)
  real, allocatable :: svtk(:,:,:,:)
  real, allocatable :: wsvtk(:,:,:,:)
  real, allocatable :: svsvt(:,:,:,:)
  real, allocatable :: svsgs(:,:,:,:)
  real, allocatable :: svsgst(:,:,:,:)

  real, allocatable :: PSS(:,:,:)
  real, allocatable :: PSSt(:,:,:)

  contains
    subroutine stats_init
      implicit none
      integer :: tVarsCount
      character(80), dimension(1,4) :: tVar0
      
      xdim = ie-ib+1
      ydim = je-jb+1
      zdim = ke-kb+1
      tsamplep = 0.
      tstatsdumpp = 0.
      dumpcount = 0
      nrect = 0

      !> Generate time averaged NetCDF: stats.xxx.xxx.xxx.nc
      if (ltdump) then
        filenamet = 'sdump.xxx.xxx.xxx.nc'
        filenamet(7:9) = cmyidx
        filenamet(11:13) = cmyidy
        filenamet(15:17) = cexpnr
        
        !> Total numbers of variables to be written
        tVarsCount = 14
        if (ltempeq) tVarsCount = tVarsCount + 4
        if (lmoist) tVarsCount = tVarsCount + 4
        if (nsv>0) tVarsCount = tVarsCount + 4*nsv
        if ((lchem) .and. (nsv>2)) tVarsCount = tVarsCount + 1

        !> Array to store the variable description of the quantities to be written
        allocate(tVars(tVarsCount,4))
        
        ctrt = 0
        call stats_allocate_vel
        if (ltempeq) call stats_allocate_temp
        if (lmoist)  call stats_allocate_moist
        if (nsv>0)   call stats_allocate_scalar
        if ((lchem) .and. (nsv>2)) call stats_allocate_PSS
        
        call ncinfo(tVar0( 1,:), 'time', 'Time', 's', 'time')
        call open_nc(filenamet, ncidt, nrect, n1=xdim, n2=ydim, n3=zdim)
        if (nrect==0) then
          call define_nc(ncidt, 1, tVar0)
          call writestat_dims_nc(ncidt)
        end if
        call define_nc(ncidt, tVarsCount, tVars)
        
        deallocate(tVars)
      end if
    end subroutine stats_init

    subroutine stats_main
      implicit none

      if (.not. rk3step==3)  return
      
      if(ltdump) then
        if (tsamplep > tsample) then
          tstatsdumppi = 1./tstatsdumpp
          call stats_compute_vel
          if (ltempeq) call stats_compute_temp
          if (lmoist)  call stats_compute_moist
          if (nsv>0)   call stats_compute_scalar
          if ((lchem) .and. (nsv>2)) call stats_compute_PSS
          tsamplep = dt
        else
          tsamplep = tsamplep + dt
        endif

        if (tstatsdumpp > tstatsdump) then
          dumpcount = dumpcount + 1
          if (myid==0) then
            write(*,*) "---------------------------------"
            write(*,*) "Writing time average statistics"
            write(*,*) "dump count ::: ", dumpcount
            write(*,*) "timee ::: ", timee
            write(*,*) "---------------------------------"
          end if
          call stats_write_vel
          if (ltempeq) call stats_write_temp
          if (lmoist)  call stats_write_moist
          if (nsv>0)   call stats_write_scalar
          if ((lchem) .and. (nsv>2)) call stats_write_PSS
          tstatsdumpp = dt
        else
          tstatsdumpp = tstatsdumpp + dt
        endif
      end if
    end subroutine stats_main


    !! ## %% Time averaging initialization routines
    subroutine stats_allocate_vel
      implicit none
      !> allocate variables to compute time-averaged quantities
      allocate(ut(ib:ie,jb:je,kb:ke+kh))   ; ut = 0;
      allocate(vt(ib:ie,jb:je,kb:ke+kh))   ; vt = 0;
      allocate(wt(ib:ie,jb:je,kb:ke+kh))   ; wt = 0;
      allocate(pt(ib:ie,jb:je,kb:ke+kh))   ; pt = 0;

      allocate(uik(ib:ie,jb:je,kb:ke+kh))
      allocate(wik(ib:ie,jb:je,kb:ke+kh))
      allocate(vjk(ib:ie,jb:je,kb:ke+kh))
      allocate(wjk(ib:ie,jb:je,kb:ke+kh))
      allocate(uij(ib:ie,jb:je,kb:ke+kh))
      allocate(vij(ib:ie,jb:je,kb:ke+kh))
      allocate(utik(ib:ie,jb:je,kb:ke+kh)) ; utik  = 0;
      allocate(wtik(ib:ie,jb:je,kb:ke+kh)) ; wtik  = 0;
      allocate(uwtik(ib:ie,jb:je,kb:ke+kh)); uwtik = 0;
      allocate(vtjk(ib:ie,jb:je,kb:ke+kh)) ; vtjk  = 0;
      allocate(wtjk(ib:ie,jb:je,kb:ke+kh)) ; wtjk  = 0;
      allocate(vwtjk(ib:ie,jb:je,kb:ke+kh)); vwtjk = 0;
      allocate(utij(ib:ie,jb:je,kb:ke+kh)) ; utij  = 0;
      allocate(vtij(ib:ie,jb:je,kb:ke+kh)) ; vtij  = 0;
      allocate(uvtij(ib:ie,jb:je,kb:ke+kh)); uvtij = 0;

      allocate(uc(ib:ie,jb:je,kb:ke+kh))
      allocate(vc(ib:ie,jb:je,kb:ke+kh))
      allocate(wc(ib:ie,jb:je,kb:ke+kh))
      allocate(utc(ib:ie,jb:je,kb:ke+kh))  ; utc  = 0;
      allocate(vtc(ib:ie,jb:je,kb:ke+kh))  ; vtc  = 0;
      allocate(wtc(ib:ie,jb:je,kb:ke+kh))  ; wtc  = 0;
      allocate(uutc(ib:ie,jb:je,kb:ke+kh)) ; uutc = 0;
      allocate(vvtc(ib:ie,jb:je,kb:ke+kh)) ; vvtc = 0;
      allocate(wwtc(ib:ie,jb:je,kb:ke+kh)) ; wwtc = 0;

      allocate(usgs(ib:ie,jb:je,kb:ke+kh))
      allocate(vsgs(ib:ie,jb:je,kb:ke+kh))
      allocate(wsgs(ib:ie,jb:je,kb:ke+kh))
      allocate(usgst(ib:ie,jb:je,kb:ke+kh)); usgst = 0;
      allocate(vsgst(ib:ie,jb:je,kb:ke+kh)); vsgst = 0;
      allocate(wsgst(ib:ie,jb:je,kb:ke+kh)); wsgst = 0;

      !> Generate variable description for the quantities to be written in the time averaged NetCDF: stats.xxx.xxx.xxx.nc
      call ncinfo( tVars(ctrt+ 1,:), 'ut'       , 'Streamwise velocity'       , 'm/s'       , 'mttt' )
      call ncinfo( tVars(ctrt+ 2,:), 'vt'       , 'Spanwise velocity'         , 'm/s'       , 'tmtt' )
      call ncinfo( tVars(ctrt+ 3,:), 'wt'       , 'Vertical velocity'         , 'm/s'       , 'ttmt' )
      call ncinfo( tVars(ctrt+ 4,:), 'pt'       , 'Kinematic Pressure'        , 'm^2/s^2'   , 'tttt' )

      call ncinfo( tVars(ctrt+ 5,:), 'upwpt'    , 'Turbulent momentum flux'   , 'm^2/s^2'   , 'mtmt' )
      call ncinfo( tVars(ctrt+ 6,:), 'vpwpt'    , 'Turbulent momentum flux'   , 'm^2/s^2'   , 'tmmt' )
      call ncinfo( tVars(ctrt+ 7,:), 'upvpt'    , 'Turbulent momentum flux'   , 'm^2/s^2'   , 'mmtt' )

      call ncinfo( tVars(ctrt+ 8,:), 'upuptc'   , 'u variance'                , 'm^2/s^2'   , 'tttt' )
      call ncinfo( tVars(ctrt+ 9,:), 'vpvptc'   , 'v variance'                , 'm^2/s^2'   , 'tttt' )
      call ncinfo( tVars(ctrt+10,:), 'wpwptc'   , 'w variance'                , 'm^2/s^2'   , 'tttt' )
      call ncinfo( tVars(ctrt+11,:), 'tketc'    , 'TKE'                       , 'm^2/s^2'   , 'tttt' )

      call ncinfo( tVars(ctrt+12,:), 'usgst'    , 'SGS u flux'                , 'm^2/s^2'   , 'mtmt' )
      call ncinfo( tVars(ctrt+13,:), 'vsgst'    , 'SGS v flux'                , 'm^2/s^2'   , 'tmmt' )
      call ncinfo( tVars(ctrt+14,:), 'wsgst'    , 'SGS w flux'                , 'm^2/s^2'   , 'ttmt' )
      ctrt = ctrt+14
    end subroutine stats_allocate_vel

    subroutine stats_allocate_temp
      implicit none
      allocate(thlk(ib:ie,jb:je,kb:ke+kh))
      allocate(thlt(ib:ie,jb:je,kb:ke+kh))   ; thlt    = 0;
      allocate(thltk(ib:ie,jb:je,kb:ke+kh))  ; thltk   = 0;
      allocate(wthltk(ib:ie,jb:je,kb:ke+kh)) ; wthltk  = 0;
      allocate(thlthlt(ib:ie,jb:je,kb:ke+kh)); thlthlt = 0;
      allocate(thlsgs(ib:ie,jb:je,kb:ke+kh))
      allocate(thlsgst(ib:ie,jb:je,kb:ke+kh)); thlsgst = 0;
      call ncinfo( tVars(ctrt+1,:) , 'thlt'     , 'Temperature'               , 'K'         , 'tttt' )
      call ncinfo( tVars(ctrt+2,:) , 'wpthlpt'  , 'Turbulent heat flux'       , 'K m/s'     , 'ttmt' )
      call ncinfo( tVars(ctrt+3,:) , 'thlpthlpt', 'Temperature variance'      , 'K^2'       , 'tttt' )
      call ncinfo( tVars(ctrt+4,:) , 'thlsgst'  , 'SGS temperature flux'      , 'K m/s'     , 'ttmt' )
      ctrt = ctrt+4
    end subroutine stats_allocate_temp

    subroutine stats_allocate_moist
      implicit none
      allocate(qtk(ib:ie,jb:je,kb:ke+kh))
      allocate(qtt(ib:ie,jb:je,kb:ke+kh))   ; qtt    = 0;
      allocate(qttk(ib:ie,jb:je,kb:ke+kh))  ; qttk   = 0;
      allocate(wqttk(ib:ie,jb:je,kb:ke+kh)) ; wqttk  = 0;
      allocate(qtqtt(ib:ie,jb:je,kb:ke+kh)) ; qtqtt  = 0;
      allocate(qtsgs(ib:ie,jb:je,kb:ke+kh))
      allocate(qtsgst(ib:ie,jb:je,kb:ke+kh)); qtsgst = 0;
      call ncinfo( tVars(ctrt+1,:) , 'qtt'      , 'Moisture'                  , 'kg/kg'     , 'tttt' )
      call ncinfo( tVars(ctrt+2,:) , 'wpqtpt'   , 'Turbulent moisture flux'   , 'kg m/kg s' , 'ttmt' )
      call ncinfo( tVars(ctrt+3,:) , 'qtpqtpt'  , 'Moisture variance'         , 'kg^2/kg^2' , 'tttt' )
      call ncinfo( tVars(ctrt+4,:) , 'qtsgst'   , 'SGS moisture flux'         , 'kg m/kg s' , 'ttmt' )
      ctrt = ctrt+4
    end subroutine stats_allocate_moist

    subroutine stats_allocate_scalar
      integer :: n
      character(2) :: sid
      allocate(svtname(nsv))
      allocate(wpsvptname(nsv))
      allocate(svpsvptname(nsv))
      allocate(svsgsname(nsv))
      allocate(svk(ib:ie,jb:je,kb:ke+kh,nsv))
      allocate(svt(ib:ie,jb:je,kb:ke+kh,nsv))   ; svt    = 0;
      allocate(svtk(ib:ie,jb:je,kb:ke+kh,nsv))  ; svtk   = 0;
      allocate(wsvtk(ib:ie,jb:je,kb:ke+kh,nsv)) ; wsvtk  = 0;
      allocate(svsvt(ib:ie,jb:je,kb:ke+kh,nsv)) ; svsvt  = 0;
      allocate(svsgs(ib:ie,jb:je,kb:ke+kh,nsv))
      allocate(svsgst(ib:ie,jb:je,kb:ke+kh,nsv)); svsgst = 0;
      do n = 1, nsv
        write (sid, '(I0)') n
        svtname(n)     = 'sca'//trim(sid)//'t'                      ! sca1t       at n = 1
        wpsvptname(n)  = 'wpsca'//trim(sid)//'pt'                   ! wpsca1pt    at n = 1
        svpsvptname(n) = 'sca'//trim(sid)//'psca'//trim(sid)//'pt'  ! sca1psca1pt at n = 1
        svsgsname(n)   = 'sv'//trim(sid)//'sgs'                     ! sv1sgs      at n = 1
        call ncinfo(tVars(ctrt+n,:)      , trim(svtname(n))    , 'Concentration field '//trim(sid)   , 'g/m^3'  , 'tttt' )
        call ncinfo(tVars(ctrt+nsv+n,:)  , trim(wpsvptname(n)) , 'Turbulent scalar flux '//trim(sid) , 'g/m^2s' , 'ttmt' )
        call ncinfo(tVars(ctrt+2*nsv+n,:), trim(svpsvptname(n)), 'Concentration variance '//trim(sid), 'g^2/m^6', 'tttt' )
        call ncinfo(tVars(ctrt+3*nsv+n,:), trim(svsgsname(n))  , 'SGS scalar flux '//trim(sid)       , 'g/m^2s' , 'ttmt' )
      end do
      ctrt = ctrt+4*nsv
    end subroutine stats_allocate_scalar

    subroutine stats_allocate_PSS
      implicit none
      allocate(PSS(ib:ie,jb:je,kb:ke+kh))
      allocate(PSSt(ib:ie,jb:je,kb:ke+kh)); PSSt = 0;
      call ncinfo( tVars(ctrt+1,:) , 'PSSt'      , 'PSS defect'                , 'gm/s'      , 'tttt' )
      ctrt = ctrt+1
    end subroutine stats_allocate_PSS


    !! ## %% Time averaging computations routines
    subroutine stats_compute_vel
      implicit none
      integer :: i, j, k
      real    :: emom
      
      !> Perform required interpolations to cell centers
      do k=kb,ke+kh
        do j=jb,je
          do i=ib,ie
            uik(i,j,k) = 0.5*dzhi(k)*(um(i,j,k)*dzf(k-1) + um(i,j,k-1)*dzf(k))
            wik(i,j,k) = 0.5*dxhi(i)*(wm(i,j,k)*dxf(i-1) + wm(i-1,j,k)*dxf(i))
            vjk(i,j,k) = 0.5*dzhi(k)*(vm(i,j,k)*dzf(k-1) + vm(i,j,k-1)*dzf(k))
            wjk(i,j,k) = 0.5*        (wm(i,j,k)          + wm(i,j-1,k))
            uij(i,j,k) = 0.5*        (um(i,j,k)          + um(i,j-1,k))
            vij(i,j,k) = 0.5*dxhi(i)*(vm(i,j,k)*dxf(i-1) + vm(i-1,j,k)*dxf(i))
            uc(i,j,k)  = 0.5*dxhi(i)*(um(i,j,k)*dxf(i-1) + um(i-1,j,k)*dxf(i))
            vc(i,j,k)  = 0.5*        (vm(i,j,k)          + vm(i,j-1,k))
            wc(i,j,k)  = 0.5*dzhi(k)*(wm(i,j,k)*dzf(k-1) + wm(i,j,k-1)*dzf(k))

            ! SGS fluxes
            ! interps ekm to cell corner (uw)
            emom = ( dzf(k-1) * ( ekm(i,j,k)*dxf(i-1)  + ekm(i-1,j,k)*dxf(i) )  + &
                     dzf(k)   * ( ekm(i,j,k-1)*dxf(i-1) + ekm(i-1,j,k-1)*dxf(i) ) )*dxhi(i) * dzhiq(k)
            usgs(i,j,k)  = emom * ( (um(i,j,k)-um(i,j,k-1)) *dzhi(k) &
                            +(wm(i,j,k)-wm(i-1,j,k))  *dxhi(i))

            ! interps ekm to cell corner (vw)
            emom = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i,j-1,k) )  + &
                     dzf(k)   * ( ekm(i,j,k-1) + ekm(i,j-1,k-1) ) ) * dzhiq(k)

            vsgs(i,j,k)  = emom * ( (vm(i,j,k)-vm(i,j,k-1)) *dzhi(k) &
                            +(wm(i,j,k)-wm(i,j-1,k))  *dyi)

          end do
        end do
      end do

      do k=kb,ke
        do j=jb,je
          do i=ib,ie
            wsgs(i,j,k) = ( ekm(i,j,k) * (wm(i,j,k+1)-wm(i,j,k)) *dzfi(k) &
                            -ekm(i,j,k-1)* (wm(i,j,k)-wm(i,j,k-1)) *dzfi(k-1) ) * 2. &
                            * dzhi(k) ! tg3315 check this
          end do
        end do
      end do

      call stats_compute_tavg(ut, um(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(vt, vm(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(wt, wm(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(pt, pres0(ib:ie,jb:je,kb:ke+kh))

      call stats_compute_tavg(utik,  uik)
      call stats_compute_tavg(wtik,  wik)
      call stats_compute_tavg(uwtik, wik*uik)

      call stats_compute_tavg(vtjk,  vjk)
      call stats_compute_tavg(wtjk,  wjk)
      call stats_compute_tavg(vwtjk, vjk*wjk)

      call stats_compute_tavg(utij,  uij)
      call stats_compute_tavg(vtij,  vij)
      call stats_compute_tavg(uvtij, uij*vij)

      call stats_compute_tavg(utc,  uc(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(vtc,  vc(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(wtc,  wc(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(uutc, uc(ib:ie,jb:je,kb:ke+kh)*uc(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(vvtc, vc(ib:ie,jb:je,kb:ke+kh)*vc(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(wwtc, wc(ib:ie,jb:je,kb:ke+kh)*wc(ib:ie,jb:je,kb:ke+kh))

      call stats_compute_tavg(usgst,  usgs)
      call stats_compute_tavg(vsgst,  vsgs)
      call stats_compute_tavg(wsgst,  wsgs)
    end subroutine stats_compute_vel

    subroutine stats_compute_temp
      implicit none
      call stats_interpolate_k(thlk, thlm(ib:ie,jb:je,kb-kh:ke+kh))
      call stats_compute_sgs(thlsgs, thlm(ib:ie,jb:je,kb-kh:ke+kh), ekh(ib:ie,jb:je,kb-kh:ke+kh))

      call stats_compute_tavg(thlt,    thlm(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(thltk,   thlk)
      call stats_compute_tavg(wthltk,  wm(ib:ie,jb:je,kb:ke+kh)*thlk)
      call stats_compute_tavg(thlthlt, thlm(ib:ie,jb:je,kb:ke+kh)*thlm(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(thlsgst, thlsgs)
    end subroutine stats_compute_temp

    subroutine stats_compute_moist
      implicit none
      call stats_interpolate_k(qtk, qtm(ib:ie,jb:je,kb-kh:ke+kh))
      call stats_compute_sgs(qtsgs, qtm(ib:ie,jb:je,kb-kh:ke+kh), ekh(ib:ie,jb:je,kb-kh:ke+kh))

      call stats_compute_tavg(qtt,    qtm(ib:ie,jb:je,kb:ke+kh) )
      call stats_compute_tavg(qttk,   qtk)
      call stats_compute_tavg(wqttk,  wm(ib:ie,jb:je,kb:ke+kh)*qtk)
      call stats_compute_tavg(qtqtt,  qtm(ib:ie,jb:je,kb:ke+kh)*qtm(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(qtsgst, qtsgs)
    end subroutine stats_compute_moist

    subroutine stats_compute_scalar
      implicit none
      integer :: n
      do n = 1, nsv
        call stats_interpolate_k(svk(:,:,:,n), svm(ib:ie,jb:je,kb-kh:ke+kh,n))
        call stats_compute_sgs(svsgs(:,:,:,n), svm(ib:ie,jb:je,kb-kh:ke+kh,n), ekh(ib:ie,jb:je,kb-kh:ke+kh))

        call stats_compute_tavg(svt(:,:,:,n)   , svm(ib:ie,jb:je,kb:ke+kh,n) )
        call stats_compute_tavg(svtk(:,:,:,n)  , svk(:,:,:,n) )
        call stats_compute_tavg(wsvtk(:,:,:,n) , wm(ib:ie,jb:je,kb:ke+kh)*svk(:,:,:,n) )
        call stats_compute_tavg(svsvt(:,:,:,n) , svm(ib:ie,jb:je,kb:ke+kh,n)*svm(ib:ie,jb:je,kb:ke+kh,n) )
        call stats_compute_tavg(svsgst(:,:,:,n), svsgs(:,:,:,n) )
      end do
    end subroutine stats_compute_scalar

    subroutine stats_compute_PSS
      implicit none
      integer :: i, j, k
      do k=kb,ke
        do j=jb,je
          do i=ib,ie
            if ((ABS(svm(i,j,k,2)) .gt. 1.e-40) .and. (IIc(i,j,k)==1)) then
              PSS(i,j,k) = ( ( (k1*(svm(i,j,k,1)/30.)*(svm(i,j,k,3)/48.))/(JNO2*(svm(i,j,k,2)/46.)) ) - 1 ) * 100
            end if
          end do
        end do
      end do
      call stats_compute_tavg(PSSt, PSS)
    end subroutine stats_compute_PSS


    !! ## %% Low level routines
    subroutine stats_compute_tavg(vart,var)
      implicit none
      real, intent(inout) :: vart(:,:,:)
      real, intent(in)    :: var(:,:,:)
      vart = ( vart*(tstatsdumpp-tsamplep) + var*tsamplep )*tstatsdumppi
    end subroutine stats_compute_tavg

    subroutine stats_interpolate_k(vark,varm)
      implicit none
      real, intent(inout) :: vark(ib:ie,jb:je,kb:ke+kh)
      real, intent(in)    :: varm(ib:ie,jb:je,kb-kh:ke+kh)
      integer :: k
      do k=kb,ke+kh
        vark(:,:,k) = 0.5*dzhi(k)*(varm(:,:,k)*dzf(k-1) + varm(:,:,k-1)*dzf(k))
      end do
    end subroutine stats_interpolate_k

    subroutine stats_compute_sgs(varsgs,varm,ekvar)
      implicit none
      real, intent(inout) :: varsgs(ib:ie,jb:je,kb:ke+kh)
      real, intent(in)    :: varm(ib:ie,jb:je,kb-kh:ke+kh)
      real, intent(in)    :: ekvar(ib:ie,jb:je,kb-kh:ke+kh)
      integer :: k
      do k=kb,ke
        varsgs(:,:,k) = 0.5 * (dzf(k-1)*ekvar(:,:,k) + dzf(k)*ekvar(:,:,k-1)) &
                            * (varm(:,:,k)-varm(:,:,k-1)) * dzh2i(k)
      end do
    end subroutine stats_compute_sgs


    !! ## %% Time averaged statistics writing routines 
    subroutine stats_write_vel
      implicit none
      real, allocatable :: upuptc(:,:,:)
      real, allocatable :: vpvptc(:,:,:)
      real, allocatable :: wpwptc(:,:,:)
      real, allocatable :: tketc(:,:,:)

      call writestat_nc(ncidt, 'time', timee, nrect, .true.)
      call writestat_nc(ncidt, 'ut', ut(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'vt', vt(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'wt', wt(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'pt', pt(:,:,kb:ke), nrect, xdim, ydim, zdim)

      call writestat_nc(ncidt, 'upwpt', uwtik(:,:,kb:ke) - utik(:,:,kb:ke)*wtik(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'vpwpt', vwtjk(:,:,kb:ke) - vtjk(:,:,kb:ke)*wtjk(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'upvpt', uvtij(:,:,kb:ke) - utij(:,:,kb:ke)*vtij(:,:,kb:ke), nrect, xdim, ydim, zdim)

      allocate(upuptc(ib:ie,jb:je,kb:ke+kh))
      allocate(vpvptc(ib:ie,jb:je,kb:ke+kh))
      allocate(wpwptc(ib:ie,jb:je,kb:ke+kh))
      allocate(tketc(ib:ie,jb:je,kb:ke+kh))
      upuptc = uutc - utc*utc
      vpvptc = vvtc - vtc*vtc
      wpwptc = wwtc - wtc*wtc
      tketc = 0.5*(upuptc + vpvptc + wpwptc)
      call writestat_nc(ncidt, 'upuptc', upuptc(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'vpvptc', vpvptc(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'wpwptc', wpwptc(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'tketc' , tketc(:,:,kb:ke) , nrect, xdim, ydim, zdim)
      deallocate(upuptc,vpvptc,wpwptc,tketc)

      call writestat_nc(ncidt, 'usgst', usgst(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'vsgst', vsgst(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'wsgst', wsgst(:,:,kb:ke), nrect, xdim, ydim, zdim)
    end subroutine stats_write_vel

    subroutine stats_write_temp
      implicit none
      call writestat_nc(ncidt, 'thlt'     , thlt(:,:,kb:ke)                                     , nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'wpthlpt'  , wthltk(:,:,kb:ke) - wt(:,:,kb:ke)*thltk(:,:,kb:ke)  , nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'thlpthlpt', thlthlt(:,:,kb:ke) - thlt(:,:,kb:ke)*thlt(:,:,kb:ke), nrect, xdim, ydim, zdim)
    end subroutine stats_write_temp

    subroutine stats_write_moist
      implicit none
      call writestat_nc(ncidt, 'qtt'    , qtt(:,:,kb:ke)                                  , nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'wpqtpt' , wqttk(:,:,kb:ke) - wt(:,:,kb:ke)*qttk(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'qtpqtpt', qtqtt(:,:,kb:ke) - qtt(:,:,kb:ke)*qtt(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'qtsgst' , qtsgst(:,:,kb:ke)                               , nrect, xdim, ydim, zdim)
    end subroutine stats_write_moist

    subroutine stats_write_scalar
      implicit none
      integer :: n
      do n = 1, nsv
        call writestat_nc(ncidt, trim(svtname(n))    , svt(:,:,kb:ke,n)                                      , nrect, xdim, ydim, zdim)
        call writestat_nc(ncidt, trim(wpsvptname(n)) , wsvtk(:,:,kb:ke,n) - wt(:,:,kb:ke)*svtk(:,:,kb:ke,n)  , nrect, xdim, ydim, zdim)
        call writestat_nc(ncidt, trim(svpsvptname(n)), svsvt(:,:,kb:ke,n) - svt(:,:,kb:ke,n)*svt(:,:,kb:ke,n), nrect, xdim, ydim, zdim)
        call writestat_nc(ncidt, trim(svsgsname(n))  , svsgst(:,:,kb:ke,n)                                   , nrect, xdim, ydim, zdim)
      end do
    end subroutine stats_write_scalar

    subroutine stats_write_PSS
      implicit none
      call writestat_nc(ncidt, 'PSSt', PSSt(:,:,kb:ke), nrect, xdim, ydim, zdim)
    end subroutine stats_write_PSS


    subroutine stats_exit
      implicit none
      if (ltdump) then
        deallocate(ut,vt,wt,pt)
        deallocate(uik,wik,vjk,wjk,uij,vij,uc,vc,wc)
        deallocate(utik,wtik,uwtik,vtjk,wtjk,vwtjk,utij,vtij,uvtij)
        deallocate(utc,vtc,wtc,uutc,vvtc,wwtc)
        deallocate(usgs,vsgs,wsgs,usgst,vsgst,wsgst)
        if (ltempeq) deallocate(thlt,thlk,thltk,wthltk,thlthlt,thlsgs,thlsgst)
        if (lmoist)  deallocate(qtt,qtk,qttk,wqttk,qtqtt,qtsgs,qtsgst)
        if (nsv>0)   deallocate(svtname,wpsvptname,svpsvptname,svsgsname,svt,svk,svtk,wsvtk,svsvt,svsgs,svsgst)
        if ((lchem) .and. (nsv>2)) deallocate(PSS,PSSt)
      end if
    end subroutine stats_exit
end module stats