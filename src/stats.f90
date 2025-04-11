module stats
  use modglobal,  only : cexpnr, ltdump, lxytdump, lxydump, ltempeq, lmoist, lchem, nsv, rk3step, &
                         ib, ie, ih, jb, je, jh, kb, ke, kh, &
                         dxf, dzf, dzfi, dxhi, dzhi, dzh2i, dyi, dzhiq, &
                         timee, tstatsdump, tsample, dt, &
                         k1, JNO2
  use modfields,  only : um, vm, wm, pres0, thlm, qtm, svm, &
                         IIc, IIu, IIus, IIv, IIvs, IIw, IIws, IIc, IIcs, IIuw, IIuws, IIvw, IIvws, IIuv, IIuvs
  use modsubgrid, only : ekh, ekm
  use modmpi,     only : cmyidx, cmyidy, myid, avexy_ibm
  use modstat_nc, only : ncinfo, open_nc, define_nc, writestat_dims_nc, writestat_nc
  implicit none
  public :: stats_init, stats_main, stats_exit
  save

  integer :: xdim, ydim, zdim
  real    :: tsamplep, tstatsdumpp, tstatsdumppi
  
  character(80)              :: filenamet
  character(80)              :: filenamexy
  character(80)              :: filenamexyt
  character(80), allocatable :: tVars(:,:), xytVars(:,:), xyVars(:,:)
  integer                    :: ctrt, ncidt, nrect, &
                                ctrxyt, ncidxyt, nrecxyt, &
                                ctrxy, ncidxy, nrecxy, &
                                dumpcount

  !!> Variables to store time averaged quantities
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

  !!> Variables to store time, y and x averaged quantities
  real, allocatable :: uxyt(:)
  real, allocatable :: vxyt(:)
  real, allocatable :: wxyt(:)
  real, allocatable :: pxyt(:)
  real, allocatable :: upwpxytik(:)
  real, allocatable :: vpwpxytjk(:)
  real, allocatable :: upvpxytij(:)
  real, allocatable :: uwxytik(:)
  real, allocatable :: vwxytjk(:)
  real, allocatable :: uvxytij(:)
  real, allocatable :: uuxyti(:)
  real, allocatable :: vvxytj(:)
  real, allocatable :: wwxytk(:)
  real, allocatable :: upupxytc(:)
  real, allocatable :: vpvpxytc(:)
  real, allocatable :: wpwpxytc(:)
  real, allocatable :: tkexytc(:)
  real, allocatable :: usgsxyt(:)
  real, allocatable :: vsgsxyt(:)
  real, allocatable :: wsgsxyt(:)

  real, allocatable :: thlxyt(:)
  real, allocatable :: wpthlpxytk(:)
  real, allocatable :: wthlxytk(:)
  real, allocatable :: thlpthlpxyt(:)
  real, allocatable :: thlsgsxyt(:)

  real, allocatable :: qtxyt(:)
  real, allocatable :: wpqtpxytk(:)
  real, allocatable :: wqtxytk(:)
  real, allocatable :: qtpqtpxyt(:)
  real, allocatable :: qtsgsxyt(:)

  !!> Variables to store y and x averaged quantities
  real, allocatable :: uxy(:)
  real, allocatable :: vxy(:)
  real, allocatable :: wxy(:)
  real, allocatable :: pxy(:)
  real, allocatable :: upwpxyik(:)
  real, allocatable :: vpwpxyjk(:)
  real, allocatable :: upvpxyij(:)
  real, allocatable :: uwxyik(:)
  real, allocatable :: uxyik(:)
  real, allocatable :: wxyik(:)
  real, allocatable :: vwxyjk(:)
  real, allocatable :: vxyjk(:)
  real, allocatable :: wxyjk(:)
  real, allocatable :: uvxyij(:)
  real, allocatable :: uxyij(:)
  real, allocatable :: vxyij(:)
  real, allocatable :: uuxyi(:)
  real, allocatable :: vvxyj(:)
  real, allocatable :: wwxyk(:)
  real, allocatable :: usgsxy(:)
  real, allocatable :: vsgsxy(:)
  real, allocatable :: wsgsxy(:)

  real, allocatable :: thlxy(:)
  real, allocatable :: wpthlpxyk(:)
  real, allocatable :: wthlxyk(:)
  real, allocatable :: thlxyk(:)
  real, allocatable :: thlsgsxy(:)

  real, allocatable :: qtxy(:)
  real, allocatable :: wpqtpxyk(:)
  real, allocatable :: wqtxyk(:)
  real, allocatable :: qtxyk(:)
  real, allocatable :: qtsgsxy(:)

  interface stats_compute_tavg
    module procedure stats_compute_tavg_1D
    module procedure stats_compute_tavg_3D
  end interface stats_compute_tavg

  contains
    subroutine stats_init
      implicit none
      integer :: tVarsCount, xytVarsCount, xyVarsCount
      character(80), dimension(1,4) :: tVar0
      
      xdim = ie-ib+1
      ydim = je-jb+1
      zdim = ke-kb+1
      tsamplep = 0.
      tstatsdumpp = 0.
      dumpcount = 0
      nrect = 0
      nrecxyt = 0
      nrecxy = 0

      call ncinfo(tVar0( 1,:), 'time', 'Time', 's', 'time')

      if(ltdump .or. lxytdump) then
        !> allocate variables to compute time-averaged quantities
        call stats_allocate_tavg_vel
        if (ltempeq) call stats_allocate_tavg_temp
        if (lmoist)  call stats_allocate_tavg_moist
      end if

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

        allocate(tVars(tVarsCount,4))   !!> Array to store the variable description of the quantities to be written
        ctrt = 0
        call stats_ncdescription_tavg_vel
        if (ltempeq) call stats_ncdescription_tavg_temp
        if (lmoist)  call stats_ncdescription_tavg_moist
        if (nsv>0)   call stats_allocate_tavg_scalar
        if ((lchem) .and. (nsv>2)) call stats_allocate_tavg_PSS
        
        call open_nc(filenamet, ncidt, nrect, n1=xdim, n2=ydim, n3=zdim)
        if (nrect==0) then
          call define_nc(ncidt, 1, tVar0)
          call writestat_dims_nc(ncidt)
        end if
        call define_nc(ncidt, tVarsCount, tVars)
        deallocate(tVars)
      end if

      !> Generate time, y and x averaged NetCDF: xytdump.xxx.nc
      if (lxytdump) then
        filenamexyt = 'xysdump.xxx.nc'
        filenamexyt(9:11) = cexpnr

        xytVarsCount = 20
        if (ltempeq) xytVarsCount = xytVarsCount + 5
        if (lmoist)  xytVarsCount = xytVarsCount + 5

        allocate(xytVars(xytVarsCount,4))   !!> Array to store the variable description of the quantities to be written
        ctrxyt = 0
        call stats_allocate_xytavg_vel
        if (ltempeq) call stats_allocate_xytavg_temp
        if (lmoist)  call stats_allocate_xytavg_moist

        if (myid==0) then
          call open_nc(filenamexyt, ncidxyt, nrecxyt, n3=zdim)
          if (nrecxyt==0) then
            call define_nc( ncidxyt, 1,  tVar0)
            call writestat_dims_nc(ncidxyt)
          end if
          call define_nc( ncidxyt, xytVarsCount, xytVars)
        end if
      end if

      !> Generate y and x averaged NetCDF: xydump.xxx.nc
      if (lxydump) then
        filenamexy = 'xydumps.xxx.nc'
        filenamexy(9:11) = cexpnr

        xyVarsCount = 16
        if (ltempeq) xyVarsCount = xyVarsCount + 4
        if (lmoist)  xyVarsCount = xyVarsCount + 4

        allocate(xyVars(xyVarsCount,4))   !!> Array to store the variable description of the quantities to be written
        ctrxy = 0
        call stats_allocate_xyavg_vel
        if (ltempeq) call stats_allocate_xyavg_temp
        if (lmoist)  call stats_allocate_xyavg_moist

        if (myid==0) then
          call open_nc(filenamexy, ncidxy, nrecxy, n3=zdim)
          if (nrecxy==0) then
            call define_nc( ncidxy, 1,  tVar0)
            call writestat_dims_nc(ncidxy)
          end if
          call define_nc( ncidxy, xyVarsCount, xyVars)
        end if
      end if
    end subroutine stats_init

    subroutine stats_main
      implicit none

      if (.not. rk3step==3)  return

      if (tsamplep > tsample) then        ! at every stats sampling instance
        tstatsdumppi = 1./tstatsdumpp

        if(ltdump .or. lxytdump .or. lxydump) then
          call stats_interpolate_and_sgs_vel
          if (ltempeq) call stats_interpolate_and_sgs_temp
          if (lmoist)  call stats_interpolate_and_sgs_moist
        end if

        if(ltdump .or. lxytdump) then
          call stats_compute_tavg_vel
          if (ltempeq) call stats_compute_tavg_temp
          if (lmoist)  call stats_compute_tavg_moist
        end if
        if(ltdump) then
          if (nsv>0) then
            call stats_interpolate_and_sgs_scalar
            call stats_compute_tavg_scalar
          end if
          if ((lchem) .and. (nsv>2)) call stats_compute_tavg_PSS
        end if

        if(lxydump) then
          call stats_compute_xyavg_vel
          if (ltempeq) call stats_compute_xyavg_temp
          if (lmoist)  call stats_compute_xyavg_moist
          if (myid==0) then
            call writestat_nc(ncidxy, 'time', timee, nrecxy, .true.)
            call stats_write_xyavg_vel
            if (ltempeq) call stats_write_xyavg_temp
            if (lmoist)  call stats_write_xyavg_moist
          end if
        end if
        
        tsamplep = dt
      else
        tsamplep = tsamplep + dt
      endif

      if (tstatsdumpp > tstatsdump) then    ! at every stats dump time instance
        dumpcount = dumpcount + 1
        if (myid==0) then
          write(*,*) "---------------------------------"
          write(*,*) "Writing time average statistics"
          write(*,*) "stats dump count ::: ", dumpcount
          write(*,*) "timee ::: ", timee
          write(*,*) "---------------------------------"
        end if
        
        if(ltdump) then
          call writestat_nc(ncidt, 'time', timee, nrect, .true.)
          call stats_write_tavg_vel
          if (ltempeq) call stats_write_tavg_temp
          if (lmoist)  call stats_write_tavg_moist
          if (nsv>0)   call stats_write_tavg_scalar
          if ((lchem) .and. (nsv>2)) call stats_write_tavg_PSS
        end if

        if(lxytdump) then
          call stats_compute_xytavg_vel
          if (ltempeq) call stats_compute_xytavg_temp
          if (lmoist)  call stats_compute_xytavg_moist
          if (myid==0) then
            call writestat_nc(ncidxyt, 'time', timee, nrecxyt, .true.)
            call stats_write_xytavg_vel
            if (ltempeq) call stats_write_xytavg_temp
            if (lmoist)  call stats_write_xytavg_moist
          end if
        end if

        tstatsdumpp = dt
      else
        tstatsdumpp = tstatsdumpp + dt
      endif

    end subroutine stats_main


    !! ## %% Time averaging initialization routines
    subroutine stats_allocate_tavg_vel
      implicit none
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
    end subroutine stats_allocate_tavg_vel
    subroutine stats_ncdescription_tavg_vel
      implicit none
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
    end subroutine stats_ncdescription_tavg_vel

    subroutine stats_allocate_tavg_temp
      implicit none
      allocate(thlk(ib:ie,jb:je,kb:ke+kh))
      allocate(thlt(ib:ie,jb:je,kb:ke+kh))   ; thlt    = 0;
      allocate(thltk(ib:ie,jb:je,kb:ke+kh))  ; thltk   = 0;
      allocate(wthltk(ib:ie,jb:je,kb:ke+kh)) ; wthltk  = 0;
      allocate(thlthlt(ib:ie,jb:je,kb:ke+kh)); thlthlt = 0;
      allocate(thlsgs(ib:ie,jb:je,kb:ke+kh))
      allocate(thlsgst(ib:ie,jb:je,kb:ke+kh)); thlsgst = 0;
    end subroutine stats_allocate_tavg_temp
    subroutine stats_ncdescription_tavg_temp
      implicit none
      call ncinfo( tVars(ctrt+1,:) , 'thlt'     , 'Temperature'               , 'K'         , 'tttt' )
      call ncinfo( tVars(ctrt+2,:) , 'wpthlpt'  , 'Turbulent heat flux'       , 'K m/s'     , 'ttmt' )
      call ncinfo( tVars(ctrt+3,:) , 'thlpthlpt', 'Temperature variance'      , 'K^2'       , 'tttt' )
      call ncinfo( tVars(ctrt+4,:) , 'thlsgst'  , 'SGS temperature flux'      , 'K m/s'     , 'ttmt' )
      ctrt = ctrt+4
    end subroutine stats_ncdescription_tavg_temp

    subroutine stats_allocate_tavg_moist
      implicit none
      allocate(qtk(ib:ie,jb:je,kb:ke+kh))
      allocate(qtt(ib:ie,jb:je,kb:ke+kh))   ; qtt    = 0;
      allocate(qttk(ib:ie,jb:je,kb:ke+kh))  ; qttk   = 0;
      allocate(wqttk(ib:ie,jb:je,kb:ke+kh)) ; wqttk  = 0;
      allocate(qtqtt(ib:ie,jb:je,kb:ke+kh)) ; qtqtt  = 0;
      allocate(qtsgs(ib:ie,jb:je,kb:ke+kh))
      allocate(qtsgst(ib:ie,jb:je,kb:ke+kh)); qtsgst = 0;
    end subroutine stats_allocate_tavg_moist
    subroutine stats_ncdescription_tavg_moist
      implicit none
      call ncinfo( tVars(ctrt+1,:) , 'qtt'      , 'Moisture'                  , 'kg/kg'     , 'tttt' )
      call ncinfo( tVars(ctrt+2,:) , 'wpqtpt'   , 'Turbulent moisture flux'   , 'kg m/kg s' , 'ttmt' )
      call ncinfo( tVars(ctrt+3,:) , 'qtpqtpt'  , 'Moisture variance'         , 'kg^2/kg^2' , 'tttt' )
      call ncinfo( tVars(ctrt+4,:) , 'qtsgst'   , 'SGS moisture flux'         , 'kg m/kg s' , 'ttmt' )
      ctrt = ctrt+4
    end subroutine stats_ncdescription_tavg_moist

    subroutine stats_allocate_tavg_scalar
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
    end subroutine stats_allocate_tavg_scalar

    subroutine stats_allocate_tavg_PSS
      implicit none
      allocate(PSS(ib:ie,jb:je,kb:ke+kh))
      allocate(PSSt(ib:ie,jb:je,kb:ke+kh)); PSSt = 0;
      call ncinfo( tVars(ctrt+1,:) , 'PSSt'      , 'PSS defect'                , 'gm/s'      , 'tttt' )
      ctrt = ctrt+1
    end subroutine stats_allocate_tavg_PSS


    !! ## %% Time, y and x averaging initialization routines
    subroutine stats_allocate_xytavg_vel
      implicit none
      allocate(uxyt(kb:ke+kh))
      allocate(vxyt(kb:ke+kh))
      allocate(wxyt(kb:ke+kh))
      allocate(pxyt(kb:ke+kh))
      allocate(upwpxytik(kb:ke+kh))
      allocate(vpwpxytjk(kb:ke+kh))
      allocate(upvpxytij(kb:ke+kh))
      allocate(uwxytik(kb:ke+kh))
      allocate(vwxytjk(kb:ke+kh))
      allocate(uvxytij(kb:ke+kh))
      allocate(uuxyti(kb:ke+kh))
      allocate(vvxytj(kb:ke+kh))
      allocate(wwxytk(kb:ke+kh))
      allocate(upupxytc(kb:ke+kh))
      allocate(vpvpxytc(kb:ke+kh))
      allocate(wpwpxytc(kb:ke+kh))
      allocate(tkexytc(kb:ke+kh))
      allocate(usgsxyt(kb:ke+kh))
      allocate(vsgsxyt(kb:ke+kh))
      allocate(wsgsxyt(kb:ke+kh))

      call ncinfo( xytVars(ctrxyt+ 1,:), 'uxyt'        , 'Streamwise velocity'      , 'm/s'       , 'tt' )
      call ncinfo( xytVars(ctrxyt+ 2,:), 'vxyt'        , 'Spanwise velocity'        , 'm/s'       , 'tt' )
      call ncinfo( xytVars(ctrxyt+ 3,:), 'wxyt'        , 'Vertical velocity'        , 'm/s'       , 'mt' )
      call ncinfo( xytVars(ctrxyt+ 4,:), 'pxyt'        , 'Kinematic Pressure'       , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+ 5,:), 'upwpxyt'     , 'Turbulent mom. flux'      , 'm^2/s^2'   , 'mt' )
      call ncinfo( xytVars(ctrxyt+ 6,:), 'vpwpxyt'     , 'Turbulent mom. flux'      , 'm^2/s^2'   , 'mt' )
      call ncinfo( xytVars(ctrxyt+ 7,:), 'upvpxyt'     , 'Turbulent mom. flux'      , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+ 8,:), 'uwxyt'       , 'Kinematic mom. flux'      , 'm^2/s^2'   , 'mt' )
      call ncinfo( xytVars(ctrxyt+ 9,:), 'vwxyt'       , 'Kinematic mom. flux'      , 'm^2/s^2'   , 'mt' )
      call ncinfo( xytVars(ctrxyt+10,:), 'uvxyt'       , 'Kinematic mom. flux'      , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+11,:), 'uuxyt'       , 'Kinematic mom. flux'      , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+12,:), 'vvxyt'       , 'Kinematic mom. flux'      , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+13,:), 'wwxyt'       , 'Kinematic mom. flux'      , 'm^2/s^2'   , 'mt' )
      call ncinfo( xytVars(ctrxyt+14,:), 'upuptxyc'    , 'u variance'               , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+15,:), 'vpvptxyc'    , 'v variance'               , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+16,:), 'wpwptxyc'    , 'w variance'               , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+17,:), 'tketxyc'     , 'TKE'                      , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+18,:), 'usgsxyt'     , 'SGS mom. flux'            , 'm^2/s^2'   , 'mt' )
      call ncinfo( xytVars(ctrxyt+19,:), 'vsgsxyt'     , 'SGS mom. flux'            , 'm^2/s^2'   , 'mt' )
      call ncinfo( xytVars(ctrxyt+20,:), 'wsgsxyt'     , 'SGS mom. flux'            , 'm^2/s^2'   , 'mt' )
      ctrxyt = ctrxyt+20
    end subroutine stats_allocate_xytavg_vel

    subroutine stats_allocate_xytavg_temp
      implicit none
      allocate(thlxyt(kb:ke+kh))
      allocate(wpthlpxytk(kb:ke+kh))
      allocate(wthlxytk(kb:ke+kh))
      allocate(thlpthlpxyt(kb:ke+kh))
      allocate(thlsgsxyt(kb:ke+kh))
      call ncinfo( xytVars(ctrxyt+ 1,:), 'thlxyt'      , 'Temperature'              , 'K'         , 'tt' )
      call ncinfo( xytVars(ctrxyt+ 2,:), 'wpthlpxyt'   , 'Turbulent heat flux'      , 'K m/s'     , 'mt' )
      call ncinfo( xytVars(ctrxyt+ 3,:), 'wthlxyt'     , 'Kinematic heat flux'      , 'K m/s'     , 'mt' )
      call ncinfo( xytVars(ctrxyt+ 4,:), 'thlpthlptxy' , 'Temp. variance'           , 'K^2'       , 'tt' )
      call ncinfo( xytVars(ctrxyt+ 5,:), 'thlsgsxyt'   , 'SGS heat flux'            , 'K m/s'     , 'mt' )
      ctrxyt = ctrxyt+5
    end subroutine stats_allocate_xytavg_temp

    subroutine stats_allocate_xytavg_moist
      implicit none
      allocate(qtxyt(kb:ke+kh))
      allocate(wpqtpxytk(kb:ke+kh))
      allocate(wqtxytk(kb:ke+kh))
      allocate(qtpqtpxyt(kb:ke+kh))
      allocate(qtsgsxyt(kb:ke+kh))
      call ncinfo( xytVars(ctrxyt+ 1,:), 'qtxyt'       , 'Moisture'                 , 'kg/kg'     , 'tt' )
      call ncinfo( xytVars(ctrxyt+ 2,:), 'wpqtpxyt'    , 'Turbulent moisture flux'  , 'kg m/kg s' , 'mt' )
      call ncinfo( xytVars(ctrxyt+ 3,:), 'wqtxyt'      , 'Kinematic moisture flux'  , 'kg m/kg s' , 'mt' )
      call ncinfo( xytVars(ctrxyt+ 4,:), 'qtpqtptxy'   , 'Moisture variance'        , 'kg^2/kg^2' , 'tt' )
      call ncinfo( xytVars(ctrxyt+ 5,:), 'qtsgsxyt'    , 'SGS moisture flux'        , 'kg m/kg s' , 'mt' )
      ctrxyt = ctrxyt+5
    end subroutine stats_allocate_xytavg_moist


    !! ## %% y and x averaging initialization routines
    subroutine stats_allocate_xyavg_vel
      implicit none
      allocate(uxy(kb:ke+kh))
      allocate(vxy(kb:ke+kh))
      allocate(wxy(kb:ke+kh))
      allocate(pxy(kb:ke+kh))
      allocate(upwpxyik(kb:ke+kh))
      allocate(vpwpxyjk(kb:ke+kh))
      allocate(upvpxyij(kb:ke+kh))
      allocate(uwxyik(kb:ke+kh))
      allocate(uxyik(kb:ke+kh))
      allocate(wxyik(kb:ke+kh))
      allocate(vwxyjk(kb:ke+kh))
      allocate(vxyjk(kb:ke+kh))
      allocate(wxyjk(kb:ke+kh))
      allocate(uvxyij(kb:ke+kh))
      allocate(uxyij(kb:ke+kh))
      allocate(vxyij(kb:ke+kh))
      allocate(uuxyi(kb:ke+kh))
      allocate(vvxyj(kb:ke+kh))
      allocate(wwxyk(kb:ke+kh))
      allocate(usgsxy(kb:ke+kh))
      allocate(vsgsxy(kb:ke+kh))
      allocate(wsgsxy(kb:ke+kh))

      call ncinfo( xyVars(ctrxy+ 1,:), 'uxy'        , 'Streamwise velocity'      , 'm/s'       , 'tt' )
      call ncinfo( xyVars(ctrxy+ 2,:), 'vxy'        , 'Spanwise velocity'        , 'm/s'       , 'tt' )
      call ncinfo( xyVars(ctrxy+ 3,:), 'wxy'        , 'Vertical velocity'        , 'm/s'       , 'mt' )
      call ncinfo( xyVars(ctrxy+ 4,:), 'pxy'        , 'Kinematic Pressure'       , 'm^2/s^2'   , 'tt' )
      call ncinfo( xyVars(ctrxy+ 5,:), 'upwpxy'     , 'Turbulent mom. flux'      , 'm^2/s^2'   , 'mt' )
      call ncinfo( xyVars(ctrxy+ 6,:), 'vpwpxy'     , 'Turbulent mom. flux'      , 'm^2/s^2'   , 'mt' )
      call ncinfo( xyVars(ctrxy+ 7,:), 'upvpxy'     , 'Turbulent mom. flux'      , 'm^2/s^2'   , 'tt' )
      call ncinfo( xyVars(ctrxy+ 8,:), 'uwxy'       , 'Kinematic mom. flux'      , 'm^2/s^2'   , 'mt' )
      call ncinfo( xyVars(ctrxy+ 9,:), 'vwxy'       , 'Kinematic mom. flux'      , 'm^2/s^2'   , 'mt' )
      call ncinfo( xyVars(ctrxy+10,:), 'uvxy'       , 'Kinematic mom. flux'      , 'm^2/s^2'   , 'tt' )
      call ncinfo( xyVars(ctrxy+11,:), 'uuxy'       , 'Kinematic mom. flux'      , 'm^2/s^2'   , 'tt' )
      call ncinfo( xyVars(ctrxy+12,:), 'vvxy'       , 'Kinematic mom. flux'      , 'm^2/s^2'   , 'tt' )
      call ncinfo( xyVars(ctrxy+13,:), 'wwxy'       , 'Kinematic mom. flux'      , 'm^2/s^2'   , 'mt' )
      call ncinfo( xyVars(ctrxy+14,:), 'usgsxy'     , 'SGS mom. flux'            , 'm^2/s^2'   , 'mt' )
      call ncinfo( xyVars(ctrxy+15,:), 'vsgsxy'     , 'SGS mom. flux'            , 'm^2/s^2'   , 'mt' )
      call ncinfo( xyVars(ctrxy+16,:), 'wsgsxy'     , 'SGS mom. flux'            , 'm^2/s^2'   , 'mt' )
      ctrxy = ctrxy+16
    end subroutine stats_allocate_xyavg_vel

    subroutine stats_allocate_xyavg_temp
      implicit none
      allocate(thlxy(kb:ke+kh))
      allocate(wpthlpxyk(kb:ke+kh))
      allocate(wthlxyk(kb:ke+kh))
      allocate(thlxyk(kb:ke+kh))
      allocate(thlsgsxy(kb:ke+kh))
      call ncinfo( xyVars(ctrxy+ 1,:), 'thlxy'      , 'Temperature'              , 'K'         , 'tt' )
      call ncinfo( xyVars(ctrxy+ 2,:), 'wpthlpxy'   , 'Turbulent heat flux'      , 'K m/s'     , 'mt' )
      call ncinfo( xyVars(ctrxy+ 3,:), 'wthlxy'     , 'Advective heat flux'      , 'K m/s'     , 'mt' )
      call ncinfo( xyVars(ctrxy+ 4,:), 'thlsgsxy'   , 'SGS heat flux'            , 'K m/s'     , 'mt' )
      ctrxy = ctrxy+4
    end subroutine stats_allocate_xyavg_temp

    subroutine stats_allocate_xyavg_moist
      implicit none
      allocate(qtxy(kb:ke+kh))
      allocate(wpqtpxyk(kb:ke+kh))
      allocate(wqtxyk(kb:ke+kh))
      allocate(qtxyk(kb:ke+kh))
      allocate(qtsgsxy(kb:ke+kh))
      call ncinfo( xyVars(ctrxy+ 1,:), 'qtxy'       , 'Moisture'                 , 'kg/kg'     , 'tt' )
      call ncinfo( xyVars(ctrxy+ 2,:), 'wpqtpxy'    , 'Turbulent moisture flux'  , 'kg m/kg s' , 'mt' )
      call ncinfo( xyVars(ctrxy+ 3,:), 'wqtxy'      , 'Advective moisture flux'  , 'kg m/kg s' , 'mt' )
      call ncinfo( xyVars(ctrxy+ 4,:), 'qtsgsxy'    , 'SGS moisture flux'        , 'kg m/kg s' , 'mt' )
      ctrxy = ctrxy+4
    end subroutine stats_allocate_xyavg_moist


    !! ## %% Interpolate variables at cell faces and compute sgs fluxes
    subroutine stats_interpolate_and_sgs_vel
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
          end do
        end do
      end do

      !> Perform required interpolations to cell centers
      do k=kb,ke+kh
        do j=jb,je
          do i=ib,ie
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
    end subroutine stats_interpolate_and_sgs_vel

    subroutine stats_interpolate_and_sgs_temp
      implicit none
      call stats_interpolate_k(thlk, thlm(ib:ie,jb:je,kb-kh:ke+kh))
      call stats_compute_sgs(thlsgs, thlm(ib:ie,jb:je,kb-kh:ke+kh), ekh(ib:ie,jb:je,kb-kh:ke+kh))
    end subroutine stats_interpolate_and_sgs_temp

    subroutine stats_interpolate_and_sgs_moist
      implicit none
      call stats_interpolate_k(qtk, qtm(ib:ie,jb:je,kb-kh:ke+kh))
      call stats_compute_sgs(qtsgs, qtm(ib:ie,jb:je,kb-kh:ke+kh), ekh(ib:ie,jb:je,kb-kh:ke+kh))
    end subroutine stats_interpolate_and_sgs_moist

    subroutine stats_interpolate_and_sgs_scalar
      implicit none
      integer :: n
      do n = 1, nsv
        call stats_interpolate_k(svk(:,:,:,n), svm(ib:ie,jb:je,kb-kh:ke+kh,n))
        call stats_compute_sgs(svsgs(:,:,:,n), svm(ib:ie,jb:je,kb-kh:ke+kh,n), ekh(ib:ie,jb:je,kb-kh:ke+kh))
      end do
    end subroutine stats_interpolate_and_sgs_scalar

    !! Low level routines
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


    !! ## %% Time averaging computations routines
    subroutine stats_compute_tavg_vel
      implicit none 
      call stats_compute_tavg(ut, um(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(vt, vm(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(wt, wm(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(pt, pres0(ib:ie,jb:je,kb:ke+kh))

      call stats_compute_tavg(utik , uik)
      call stats_compute_tavg(wtik , wik)
      call stats_compute_tavg(uwtik, wik*uik)

      call stats_compute_tavg(vtjk , vjk)
      call stats_compute_tavg(wtjk , wjk)
      call stats_compute_tavg(vwtjk, vjk*wjk)

      call stats_compute_tavg(utij , uij)
      call stats_compute_tavg(vtij , vij)
      call stats_compute_tavg(uvtij, uij*vij)

      call stats_compute_tavg(utc , uc(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(vtc , vc(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(wtc , wc(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(uutc, uc(ib:ie,jb:je,kb:ke+kh)*uc(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(vvtc, vc(ib:ie,jb:je,kb:ke+kh)*vc(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(wwtc, wc(ib:ie,jb:je,kb:ke+kh)*wc(ib:ie,jb:je,kb:ke+kh))

      call stats_compute_tavg(usgst,  usgs)
      call stats_compute_tavg(vsgst,  vsgs)
      call stats_compute_tavg(wsgst,  wsgs)
    end subroutine stats_compute_tavg_vel

    subroutine stats_compute_tavg_temp
      implicit none
      call stats_compute_tavg(thlt   , thlm(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(thltk  , thlk)
      call stats_compute_tavg(wthltk , wm(ib:ie,jb:je,kb:ke+kh)*thlk)
      call stats_compute_tavg(thlthlt, thlm(ib:ie,jb:je,kb:ke+kh)*thlm(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(thlsgst, thlsgs)
    end subroutine stats_compute_tavg_temp

    subroutine stats_compute_tavg_moist
      implicit none
      call stats_compute_tavg(qtt   , qtm(ib:ie,jb:je,kb:ke+kh) )
      call stats_compute_tavg(qttk  , qtk)
      call stats_compute_tavg(wqttk , wm(ib:ie,jb:je,kb:ke+kh)*qtk)
      call stats_compute_tavg(qtqtt , qtm(ib:ie,jb:je,kb:ke+kh)*qtm(ib:ie,jb:je,kb:ke+kh))
      call stats_compute_tavg(qtsgst, qtsgs)
    end subroutine stats_compute_tavg_moist

    subroutine stats_compute_tavg_scalar
      implicit none
      integer :: n
      do n = 1, nsv
        call stats_compute_tavg(svt(:,:,:,n)   , svm(ib:ie,jb:je,kb:ke+kh,n) )
        call stats_compute_tavg(svtk(:,:,:,n)  , svk(:,:,:,n) )
        call stats_compute_tavg(wsvtk(:,:,:,n) , wm(ib:ie,jb:je,kb:ke+kh)*svk(:,:,:,n) )
        call stats_compute_tavg(svsvt(:,:,:,n) , svm(ib:ie,jb:je,kb:ke+kh,n)*svm(ib:ie,jb:je,kb:ke+kh,n) )
        call stats_compute_tavg(svsgst(:,:,:,n), svsgs(:,:,:,n) )
      end do
    end subroutine stats_compute_tavg_scalar

    subroutine stats_compute_tavg_PSS
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
    end subroutine stats_compute_tavg_PSS

    !! Low level tavg routines
    subroutine stats_compute_tavg_1D(vart,var)
      implicit none
      real, intent(inout) :: vart(:)
      real, intent(in)    :: var(:)
      vart = ( vart*(tstatsdumpp-tsamplep) + var*tsamplep )*tstatsdumppi
    end subroutine stats_compute_tavg_1D
    subroutine stats_compute_tavg_3D(vart,var)
      implicit none
      real, intent(inout) :: vart(:,:,:)
      real, intent(in)    :: var(:,:,:)
      vart = ( vart*(tstatsdumpp-tsamplep) + var*tsamplep )*tstatsdumppi
    end subroutine stats_compute_tavg_3D


    !! ## %% Time, y and x averaging computations routines
    subroutine stats_compute_xytavg_vel
      implicit none
      !> Mean
      uxyt = 0.; vxyt = 0.; wxyt = 0; pxyt = 0;
      call avexy_ibm(uxyt,ut(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIu(ib:ie,jb:je,kb:ke+kh),IIus(kb:ke+kh),.false.)
      call avexy_ibm(vxyt,vt(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIv(ib:ie,jb:je,kb:ke+kh),IIvs(kb:ke+kh),.false.)
      call avexy_ibm(wxyt,wt(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
      call avexy_ibm(pxyt,pt(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)

      !> Turbulent fluxes
      upwpxytik = 0; vpwpxytjk = 0; upvpxytij = 0;
      call avexy_ibm(upwpxytik,uwtik(ib:ie,jb:je,kb:ke+kh)-utik(ib:ie,jb:je,kb:ke+kh)*wtik(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.false.)
      call avexy_ibm(vpwpxytjk,vwtjk(ib:ie,jb:je,kb:ke+kh)-vtjk(ib:ie,jb:je,kb:ke+kh)*wtjk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.false.)
      call avexy_ibm(upvpxytij,uvtij(ib:ie,jb:je,kb:ke+kh)-utij(ib:ie,jb:je,kb:ke+kh)*vtij(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuv(ib:ie,jb:je,kb:ke+kh),IIuvs(kb:ke+kh),.false.)

      !> Advective fluxes
      uwxytik = 0; vwxytjk = 0; uvxytij = 0; uuxyti = 0; vvxytj = 0; wwxytk = 0;
      call avexy_ibm(uwxytik,utik(ib:ie,jb:je,kb:ke+kh)*wtik(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.false.)
      call avexy_ibm(vwxytjk,vtjk(ib:ie,jb:je,kb:ke+kh)*wtjk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.false.)
      call avexy_ibm(uvxytij,utij(ib:ie,jb:je,kb:ke+kh)*vtij(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuv(ib:ie,jb:je,kb:ke+kh),IIuvs(kb:ke+kh),.false.)
      call avexy_ibm(uuxyti,ut(ib:ie,jb:je,kb:ke+kh)*ut(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIu(ib:ie,jb:je,kb:ke+kh),IIus(kb:ke+kh),.false.)
      call avexy_ibm(vvxytj,vt(ib:ie,jb:je,kb:ke+kh)*vt(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIv(ib:ie,jb:je,kb:ke+kh),IIvs(kb:ke+kh),.false.)
      call avexy_ibm(wwxytk,wt(ib:ie,jb:je,kb:ke+kh)*wt(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)

      !> Variances and TKE
      upupxytc = 0; vpvpxytc = 0; wpwpxytc = 0; tkexytc = 0;
      call avexy_ibm(upupxytc,uutc(ib:ie,jb:je,kb:ke+kh)-utc(ib:ie,jb:je,kb:ke+kh)*utc(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
      call avexy_ibm(vpvpxytc,vvtc(ib:ie,jb:je,kb:ke+kh)-vtc(ib:ie,jb:je,kb:ke+kh)*vtc(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
      call avexy_ibm(wpwpxytc,wwtc(ib:ie,jb:je,kb:ke+kh)-wtc(ib:ie,jb:je,kb:ke+kh)*wtc(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
      call avexy_ibm(tkexytc,0.5*((wwtc(ib:ie,jb:je,kb:ke+kh)-wtc(ib:ie,jb:je,kb:ke+kh)*wtc(ib:ie,jb:je,kb:ke+kh))+(vvtc(ib:ie,jb:je,kb:ke+kh)-vtc(ib:ie,jb:je,kb:ke+kh)*vtc(ib:ie,jb:je,kb:ke+kh))+(uutc(ib:ie,jb:je,kb:ke+kh)-utc(ib:ie,jb:je,kb:ke+kh)*utc(ib:ie,jb:je,kb:ke+kh))),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)

      !> SGS fluxes
      usgsxyt = 0; vsgsxyt = 0; wsgsxyt = 0;
      call avexy_ibm(usgsxyt,usgst(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.false.)
      call avexy_ibm(vsgsxyt,vsgst(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.false.)
      call avexy_ibm(wsgsxyt,wsgst(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
    end subroutine stats_compute_xytavg_vel

    subroutine stats_compute_xytavg_temp
      implicit none
      thlxyt = 0.; wpthlpxytk = 0; wthlxytk = 0; thlpthlpxyt = 0; thlsgsxyt = 0;
      call avexy_ibm(thlxyt,thlt(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
      call avexy_ibm(wpthlpxytk,wthltk(ib:ie,jb:je,kb:ke+kh)-wt(ib:ie,jb:je,kb:ke+kh)*thltk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
      call avexy_ibm(wthlxytk,wt(ib:ie,jb:je,kb:ke+kh)*thltk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
      call avexy_ibm(thlpthlpxyt,thlthlt(ib:ie,jb:je,kb:ke+kh)-thlt(ib:ie,jb:je,kb:ke+kh)*thlt(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
      call avexy_ibm(thlsgsxyt,thlsgst(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
    end subroutine stats_compute_xytavg_temp

    subroutine stats_compute_xytavg_moist
      implicit none
      qtxyt = 0.; wpqtpxytk = 0; wqtxytk = 0; qtpqtpxyt = 0; qtsgsxyt = 0;
      call avexy_ibm(qtxyt,qtt(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
      call avexy_ibm(wpqtpxytk,wqttk(ib:ie,jb:je,kb:ke+kh)-wt(ib:ie,jb:je,kb:ke+kh)*qttk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
      call avexy_ibm(wqtxytk,wt(ib:ie,jb:je,kb:ke+kh)*qttk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
      call avexy_ibm(qtpqtpxyt,qtqtt(ib:ie,jb:je,kb:ke+kh)-qtt(ib:ie,jb:je,kb:ke+kh)*qtt(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
      call avexy_ibm(qtsgsxyt,qtsgst(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
    end subroutine stats_compute_xytavg_moist


    !! ## %% y and x averaging computations routines
    subroutine stats_compute_xyavg_vel
      implicit none
      !> Mean
      uxy = 0.; vxy = 0.; wxy = 0.; pxy = 0.;
      call avexy_ibm(uxy,um(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIu(ib:ie,jb:je,kb:ke+kh),IIus(kb:ke+kh),.false.)
      call avexy_ibm(vxy,vm(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIv(ib:ie,jb:je,kb:ke+kh),IIvs(kb:ke+kh),.false.)
      call avexy_ibm(wxy,wm(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
      call avexy_ibm(pxy,pres0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)

      !> Advective fluxes and some necesseary mean
      uwxyik = 0.; uxyik = 0.; wxyik = 0.;
      vwxyjk = 0.; vxyjk = 0.; wxyjk = 0.;
      uvxyij = 0.; uxyij = 0.; vxyij = 0.;
      uuxyi  = 0.; vvxyj = 0.; wwxyk = 0.;
      call avexy_ibm(uwxyik,uik(ib:ie,jb:je,kb:ke+kh)*wik(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.true.)
      call avexy_ibm(uxyik,uik(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.true.)
      call avexy_ibm(wxyik,wik(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.true.)
      call avexy_ibm(vwxyjk,vjk(ib:ie,jb:je,kb:ke+kh)*wjk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.true.)
      call avexy_ibm(vxyjk,vjk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.true.)
      call avexy_ibm(wxyjk,wjk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.true.)
      call avexy_ibm(uvxyij,uij(ib:ie,jb:je,kb:ke+kh)*vij(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuv(ib:ie,jb:je,kb:ke+kh),IIuvs(kb:ke+kh),.true.)
      call avexy_ibm(uxyij,uij(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuv(ib:ie,jb:je,kb:ke+kh),IIuvs(kb:ke+kh),.true.)
      call avexy_ibm(vxyij,vij(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuv(ib:ie,jb:je,kb:ke+kh),IIuvs(kb:ke+kh),.true.)

      call avexy_ibm(uuxyi,um(ib:ie,jb:je,kb:ke+kh)*um(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIu(ib:ie,jb:je,kb:ke+kh),IIus(kb:ke+kh),.true.)
      call avexy_ibm(vvxyj,vm(ib:ie,jb:je,kb:ke+kh)*vm(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIv(ib:ie,jb:je,kb:ke+kh),IIvs(kb:ke+kh),.true.)
      call avexy_ibm(wwxyk,wm(ib:ie,jb:je,kb:ke+kh)*wm(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.true.)

      !> Turbulent fluxes
      upwpxyik = uwxyik - uxyik*wxyik
      vpwpxyjk = vwxyjk - vxyjk*wxyjk
      upvpxyij = uvxyij - uxyij*vxyij

      !> SGS fluxes
      usgsxy = 0.; vsgsxy = 0.; wsgsxy = 0.;
      call avexy_ibm(usgsxy,usgs(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.false.)
      call avexy_ibm(vsgsxy,vsgs(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.false.)
      call avexy_ibm(wsgsxy,wsgs(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
    end subroutine stats_compute_xyavg_vel

    subroutine stats_compute_xyavg_temp
      implicit none
      thlxy = 0.; wthlxyk = 0.; thlxyk = 0.; wpthlpxyk = 0.; thlsgsxy = 0.;

      call avexy_ibm(thlxy,thlm(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
      
      call avexy_ibm(wthlxyk,wm(ib:ie,jb:je,kb:ke+kh)*thlk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.true.)
      call avexy_ibm(thlxyk,thlk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.true.)
      
      wpthlpxyk = wthlxyk - wxy*thlxyk

      call avexy_ibm(thlsgsxy,thlsgs(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
    end subroutine stats_compute_xyavg_temp

    subroutine stats_compute_xyavg_moist
      implicit none
      qtxy = 0.; wqtxyk = 0.; qtxyk = 0.; wpqtpxyk = 0.; qtsgsxy = 0.;

      call avexy_ibm(qtxy,qtm(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)

      call avexy_ibm(wqtxyk,wm(ib:ie,jb:je,kb:ke+kh)*qtk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.true.)
      call avexy_ibm(qtxyk,qtk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.true.)
      
      wpqtpxyk = wqtxyk - wxy*qtxyk

      call avexy_ibm(qtsgsxy,qtsgs(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
    end subroutine stats_compute_xyavg_moist


    !! ## %% Time averaged statistics writing routines 
    subroutine stats_write_tavg_vel
      implicit none
      call writestat_nc(ncidt, 'ut', ut(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'vt', vt(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'wt', wt(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'pt', pt(:,:,kb:ke), nrect, xdim, ydim, zdim)

      call writestat_nc(ncidt, 'upwpt', uwtik(:,:,kb:ke) - utik(:,:,kb:ke)*wtik(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'vpwpt', vwtjk(:,:,kb:ke) - vtjk(:,:,kb:ke)*wtjk(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'upvpt', uvtij(:,:,kb:ke) - utij(:,:,kb:ke)*vtij(:,:,kb:ke), nrect, xdim, ydim, zdim)
      
      call writestat_nc(ncidt, 'upuptc', uutc(:,:,kb:ke)-utc(:,:,kb:ke)*utc(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'vpvptc', vvtc(:,:,kb:ke)-vtc(:,:,kb:ke)*vtc(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'wpwptc', wwtc(:,:,kb:ke)-wtc(:,:,kb:ke)*wtc(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'tketc' , 0.5*( (uutc(:,:,kb:ke)-utc(:,:,kb:ke)*utc(:,:,kb:ke)) + (vvtc(:,:,kb:ke)-vtc(:,:,kb:ke)*vtc(:,:,kb:ke)) + (wwtc(:,:,kb:ke)-wtc(:,:,kb:ke)*wtc(:,:,kb:ke)) ) , nrect, xdim, ydim, zdim)

      call writestat_nc(ncidt, 'usgst', usgst(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'vsgst', vsgst(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'wsgst', wsgst(:,:,kb:ke), nrect, xdim, ydim, zdim)
    end subroutine stats_write_tavg_vel

    subroutine stats_write_tavg_temp
      implicit none
      call writestat_nc(ncidt, 'thlt'     , thlt(:,:,kb:ke)                                     , nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'wpthlpt'  , wthltk(:,:,kb:ke) - wt(:,:,kb:ke)*thltk(:,:,kb:ke)  , nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'thlpthlpt', thlthlt(:,:,kb:ke) - thlt(:,:,kb:ke)*thlt(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'thlsgst'  , thlsgst(:,:,kb:ke)                                  , nrect, xdim, ydim, zdim)
    end subroutine stats_write_tavg_temp

    subroutine stats_write_tavg_moist
      implicit none
      call writestat_nc(ncidt, 'qtt'    , qtt(:,:,kb:ke)                                  , nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'wpqtpt' , wqttk(:,:,kb:ke) - wt(:,:,kb:ke)*qttk(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'qtpqtpt', qtqtt(:,:,kb:ke) - qtt(:,:,kb:ke)*qtt(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'qtsgst' , qtsgst(:,:,kb:ke)                               , nrect, xdim, ydim, zdim)
    end subroutine stats_write_tavg_moist

    subroutine stats_write_tavg_scalar
      implicit none
      integer :: n
      do n = 1, nsv
        call writestat_nc(ncidt, trim(svtname(n))    , svt(:,:,kb:ke,n)                                      , nrect, xdim, ydim, zdim)
        call writestat_nc(ncidt, trim(wpsvptname(n)) , wsvtk(:,:,kb:ke,n) - wt(:,:,kb:ke)*svtk(:,:,kb:ke,n)  , nrect, xdim, ydim, zdim)
        call writestat_nc(ncidt, trim(svpsvptname(n)), svsvt(:,:,kb:ke,n) - svt(:,:,kb:ke,n)*svt(:,:,kb:ke,n), nrect, xdim, ydim, zdim)
        call writestat_nc(ncidt, trim(svsgsname(n))  , svsgst(:,:,kb:ke,n)                                   , nrect, xdim, ydim, zdim)
      end do
    end subroutine stats_write_tavg_scalar

    subroutine stats_write_tavg_PSS
      implicit none
      call writestat_nc(ncidt, 'PSSt', PSSt(:,:,kb:ke), nrect, xdim, ydim, zdim)
    end subroutine stats_write_tavg_PSS


    !! ## %% Time, y and x averaged statistics writing routines 
    subroutine stats_write_xytavg_vel
      implicit none
      call writestat_nc(ncidxyt, 'uxyt'       , uxyt(kb:ke)       , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'vxyt'       , vxyt(kb:ke)       , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'wxyt'       , wxyt(kb:ke)       , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'pxyt'       , pxyt(kb:ke)       , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'upwpxyt'    , upwpxytik(kb:ke)  , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'vpwpxyt'    , vpwpxytjk(kb:ke)  , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'upvpxyt'    , upvpxytij(kb:ke)  , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'uwxyt'      , uwxytik(kb:ke)    , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'vwxyt'      , vwxytjk(kb:ke)    , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'uvxyt'      , uvxytij(kb:ke)    , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'uuxyt'      , uuxyti(kb:ke)     , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'vvxyt'      , vvxytj(kb:ke)     , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'wwxyt'      , wwxytk(kb:ke)     , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'upuptxyc'   , upupxytc(kb:ke)   , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'vpvptxyc'   , vpvpxytc(kb:ke)   , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'wpwptxyc'   , wpwpxytc(kb:ke)   , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'tketxyc'    , tkexytc(kb:ke)    , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'usgsxyt'    , usgsxyt(kb:ke)    , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'vsgsxyt'    , vsgsxyt(kb:ke)    , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'wsgsxyt'    , wsgsxyt(kb:ke)    , nrecxyt, zdim)
    end subroutine stats_write_xytavg_vel
    
    subroutine stats_write_xytavg_temp
      implicit none
      call writestat_nc(ncidxyt, 'thlxyt'     , thlxyt(kb:ke)     , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'wpthlpxyt'  , wpthlpxytk(kb:ke) , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'wthlxyt'    , wthlxytk(kb:ke)   , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'thlpthlptxy', thlpthlpxyt(kb:ke), nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'thlsgsxyt'  , thlsgsxyt(kb:ke)  , nrecxyt, zdim)
    end subroutine stats_write_xytavg_temp
    
    subroutine stats_write_xytavg_moist
      implicit none
      call writestat_nc(ncidxyt, 'qtxyt'      , qtxyt(kb:ke)      , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'wpqtpxyt'   , wpqtpxytk(kb:ke)  , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'wqtxyt'     , wqtxytk(kb:ke)    , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'qtpqtptxy'  , qtpqtpxyt(kb:ke)  , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'qtsgsxyt'   , qtsgsxyt(kb:ke)   , nrecxyt, zdim)
    end subroutine stats_write_xytavg_moist


    !! ## %% y and x averaged statistics writing routines 
    subroutine stats_write_xyavg_vel
      implicit none
      call writestat_nc(ncidxy, 'uxy'       , uxy(kb:ke)       , nrecxy, zdim)
      call writestat_nc(ncidxy, 'vxy'       , vxy(kb:ke)       , nrecxy, zdim)
      call writestat_nc(ncidxy, 'wxy'       , wxy(kb:ke)       , nrecxy, zdim)
      call writestat_nc(ncidxy, 'pxy'       , pxy(kb:ke)       , nrecxy, zdim)
      call writestat_nc(ncidxy, 'upwpxy'    , upwpxyik(kb:ke)  , nrecxy, zdim)
      call writestat_nc(ncidxy, 'vpwpxy'    , vpwpxyjk(kb:ke)  , nrecxy, zdim)
      call writestat_nc(ncidxy, 'upvpxy'    , upvpxyij(kb:ke)  , nrecxy, zdim)
      call writestat_nc(ncidxy, 'uwxy'      , uwxyik(kb:ke)    , nrecxy, zdim)
      call writestat_nc(ncidxy, 'vwxy'      , vwxyjk(kb:ke)    , nrecxy, zdim)
      call writestat_nc(ncidxy, 'uvxy'      , uvxyij(kb:ke)    , nrecxy, zdim)
      call writestat_nc(ncidxy, 'uuxy'      , uuxyi(kb:ke)     , nrecxy, zdim)
      call writestat_nc(ncidxy, 'vvxy'      , vvxyj(kb:ke)     , nrecxy, zdim)
      call writestat_nc(ncidxy, 'wwxy'      , wwxyk(kb:ke)     , nrecxy, zdim)
      call writestat_nc(ncidxy, 'usgsxy'    , usgsxy(kb:ke)    , nrecxy, zdim)
      call writestat_nc(ncidxy, 'vsgsxy'    , vsgsxy(kb:ke)    , nrecxy, zdim)
      call writestat_nc(ncidxy, 'wsgsxy'    , wsgsxy(kb:ke)    , nrecxy, zdim)
    end subroutine stats_write_xyavg_vel

    subroutine stats_write_xyavg_temp
      implicit none
      call writestat_nc(ncidxy, 'thlxy'     , thlxy(kb:ke)     , nrecxy, zdim)
      call writestat_nc(ncidxy, 'wpthlpxy'  , wpthlpxyk(kb:ke) , nrecxy, zdim)
      call writestat_nc(ncidxy, 'wthlxy'    , wthlxyk(kb:ke)   , nrecxy, zdim)
      call writestat_nc(ncidxy, 'thlsgsxy'  , thlsgsxy(kb:ke)  , nrecxy, zdim)
    end subroutine stats_write_xyavg_temp

    subroutine stats_write_xyavg_moist
      implicit none
      call writestat_nc(ncidxy, 'qtxy'      , qtxy(kb:ke)      , nrecxy, zdim)
      call writestat_nc(ncidxy, 'wpqtpxy'   , wpqtpxyk(kb:ke)  , nrecxy, zdim)
      call writestat_nc(ncidxy, 'wqtxy'     , wqtxyk(kb:ke)    , nrecxy, zdim)
      call writestat_nc(ncidxy, 'qtsgsxy'   , qtsgsxy(kb:ke)   , nrecxy, zdim)
    end subroutine stats_write_xyavg_moist


    subroutine stats_exit
      implicit none
      if (ltdump .or. lxytdump) then
        deallocate(ut,vt,wt,pt)
        deallocate(uik,wik,vjk,wjk,uij,vij,uc,vc,wc)
        deallocate(utc,vtc,wtc,uutc,vvtc,wwtc)
        deallocate(utik,wtik,uwtik,vtjk,wtjk,vwtjk,utij,vtij,uvtij)
        deallocate(usgs,vsgs,wsgs,usgst,vsgst,wsgst)
        if (ltempeq) deallocate(thlt,thlk,thltk,wthltk,thlthlt,thlsgs,thlsgst)
        if (lmoist)  deallocate(qtt,qtk,qttk,wqttk,qtqtt,qtsgs,qtsgst)
      end if
      if (ltdump) then  
        if (nsv>0)   deallocate(svtname,wpsvptname,svpsvptname,svsgsname,svt,svk,svtk,wsvtk,svsvt,svsgs,svsgst)
        if ((lchem) .and. (nsv>2)) deallocate(PSS,PSSt)
      end if
      if (lxytdump) then
        deallocate(uxyt,vxyt,wxyt,pxyt,usgsxyt,vsgsxyt,wsgsxyt)
        deallocate(upwpxytik,vpwpxytjk,upvpxytij,upupxytc,vpvpxytc,wpwpxytc,tkexytc)
        deallocate(uwxytik,vwxytjk,uvxytij,uuxyti,vvxytj,wwxytk)
        if (ltempeq) deallocate(thlxyt,wpthlpxytk,wthlxytk,thlpthlpxyt,thlsgsxyt)
        if (lmoist)  deallocate(qtxyt,wpqtpxytk,wqtxytk,qtpqtpxyt,qtsgsxyt)
      end if
      if (lxydump) then
        deallocate(uxy,vxy,wxy,pxy,usgsxy,vsgsxy,wsgsxy)
        deallocate(upwpxyik,vpwpxyjk,upvpxyij)
        deallocate(uwxyik,uxyik,wxyik,vwxyjk,vxyjk,wxyjk,uvxyij,uxyij,vxyij,uuxyi,vvxyj,wwxyk)
        if (ltempeq) deallocate(thlxy,wpthlpxyk,wthlxyk,thlxyk,thlsgsxy)
        if (lmoist)  deallocate(qtxy,wpqtpxyk,wqtxyk,qtxyk,qtsgsxy)
      end if
    end subroutine stats_exit
end module stats