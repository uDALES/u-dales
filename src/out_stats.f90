module stats
  use modglobal,  only : cexpnr, ltempeq, lmoist, lchem, nsv, rk3step, &
                         ltdump, lxytdump, lxydump, lytdump, lydump, ltreedump, &
                         ib, ie, ih, jb, je, jh, kb, ke, kh, &
                         dxf, dzf, dzfi, dxhi, dzhi, dzh2i, dyi, dzhiq, &
                         timee, tstatsdump, tsample, dt, &
                         k1, JNO2
  use modfields,  only : um, vm, wm, pres0, thlm, qtm, svm, &
                         IIu, IIus, IIut, IIv, IIvs, IIvt, IIw, IIws, IIwt, IIc, IIcs, IIct, &
                         IIuw, IIuws, IIuwt, IIvw, IIvws, IIuv, IIuvs, &
                         tr_u, tr_v, tr_w, tr_thl, tr_qt, tr_qtR, tr_qtA, tr_omega, tr_sv
  use modsubgrid, only : ekh, ekm
  use modmpi,     only : cmyidx, cmyidy, myid, myidy, spatial_avg
  use modstat_nc, only : ncinfo, open_nc, define_nc, writestat_dims_nc, writestat_nc
  implicit none
  private
  public :: stats_init, stats_main, stats_exit
  save

  integer :: xdim, ydim, zdim
  real    :: tsamplep, tstatsdumpp, tstatsdumppi
  
  character(80)              :: filenamet
  character(80)              :: filenamexyt
  character(80)              :: filenamexy
  character(80)              :: filenameyt
  character(80)              :: filenamey
  character(80)              :: filenametree
  character(80)              :: timeVar(1,4)
  character(80), allocatable :: tVars(:,:), xytVars(:,:), xyVars(:,:), ytVars(:,:), yVars(:,:), treeVars(:,:)
  integer                    :: tVarsCount, xytVarsCount, xyVarsCount, ytVarsCount, yVarsCount, treeVarsCount
  integer                    :: ctrt,    ncidt,    nrect, &
                                ctrxyt,  ncidxyt,  nrecxyt, &
                                ctrxy,   ncidxy,   nrecxy, &
                                ctryt,   ncidyt,   nrecyt, &
                                ctry,    ncidy,    nrecy, &
                                ctrtree, ncidtree, nrectree, &
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

  !!> Variables to store time and y averaged quantities
  real, allocatable :: uyt(:,:)
  real, allocatable :: vyt(:,:)
  real, allocatable :: wyt(:,:)
  real, allocatable :: pyt(:,:)
  real, allocatable :: upwpytik(:,:)
  real, allocatable :: uwytik(:,:)
  real, allocatable :: upupytc(:,:)
  real, allocatable :: vpvpytc(:,:)
  real, allocatable :: wpwpytc(:,:)
  real, allocatable :: usgsyt(:,:)
  real, allocatable :: wsgsyt(:,:)

  real, allocatable :: thlyt(:,:)
  real, allocatable :: wpthlpytk(:,:)
  real, allocatable :: wthlytk(:,:)
  real, allocatable :: thlpthlpyt(:,:)
  real, allocatable :: thlsgsyt(:,:)

  real, allocatable :: qtyt(:,:)
  real, allocatable :: wpqtpytk(:,:)
  real, allocatable :: wqtytk(:,:)
  real, allocatable :: qtpqtpyt(:,:)
  real, allocatable :: qtsgsyt(:,:)

  character(10), allocatable :: svytname(:)
  character(20), allocatable :: wpsvpytname(:)
  character(20), allocatable :: wsvytname(:)
  character(20), allocatable :: svpsvpytname(:)
  character(10), allocatable :: svsgsytname(:)
  real, allocatable :: svyt(:,:,:)
  real, allocatable :: wpsvpytk(:,:,:)
  real, allocatable :: wsvytk(:,:,:)
  real, allocatable :: svpsvpyt(:,:,:)
  real, allocatable :: svsgsyt(:,:,:)

  !!> Variables to store y averaged quantities
  real, allocatable :: uy(:,:)
  real, allocatable :: vy(:,:)
  real, allocatable :: wy(:,:)
  real, allocatable :: py(:,:)
  real, allocatable :: upwpyik(:,:)
  real, allocatable :: uwyik(:,:)
  real, allocatable :: uyik(:,:)
  real, allocatable :: wyik(:,:)
  real, allocatable :: usgsy(:,:)
  real, allocatable :: wsgsy(:,:)

  real, allocatable :: thly(:,:)
  real, allocatable :: wpthlpyk(:,:)
  real, allocatable :: wthlyk(:,:)
  real, allocatable :: thlyk(:,:)
  real, allocatable :: thlsgsy(:,:)

  real, allocatable :: qty(:,:)
  real, allocatable :: wpqtpyk(:,:)
  real, allocatable :: wqtyk(:,:)
  real, allocatable :: qtyk(:,:)
  real, allocatable :: qtsgsy(:,:)

  character(10), allocatable :: svyname(:)
  character(20), allocatable :: wpsvpyname(:)
  character(20), allocatable :: wsvyname(:)
  character(10), allocatable :: svsgsyname(:)
  real, allocatable :: svy(:,:,:)
  real, allocatable :: wpsvpyk(:,:,:)
  real, allocatable :: wsvyk(:,:,:)
  real, allocatable :: svyk(:,:,:)
  real, allocatable :: svsgsy(:,:,:)

  !!> Variables to store time averaged tree quantities
  real, allocatable :: tr_ut(:,:,:)
  real, allocatable :: tr_vt(:,:,:)
  real, allocatable :: tr_wt(:,:,:)
  real, allocatable :: tr_thlt(:,:,:)
  real, allocatable :: tr_qtt(:,:,:)
  real, allocatable :: tr_qtRt(:,:,:)
  real, allocatable :: tr_qtAt(:,:,:)
  real, allocatable :: tr_omegat(:,:,:)
  character(10), allocatable :: svtreename(:)
  real, allocatable :: tr_svt(:,:,:,:)


  interface stats_compute_tavg
    module procedure stats_compute_tavg_1D
    module procedure stats_compute_tavg_3D
  end interface stats_compute_tavg


  contains
    subroutine stats_init
      implicit none
      
      if(.not.(ltdump .or. lxytdump .or. lxydump .or. lytdump .or. lydump .or. ltreedump)) return
      
      xdim = ie-ib+1
      ydim = je-jb+1
      zdim = ke-kb+1
      tsamplep = 0.
      tstatsdumpp = 0.
      dumpcount = 0
      call ncinfo(timeVar( 1,:), 'time', 'Time', 's', 'time')

      if(ltdump .or. lxytdump .or. lxydump .or. lytdump .or. lydump) then
        call stats_allocate_interp_and_sgs_vel
        if (ltempeq) call stats_allocate_interp_and_sgs_temp
        if (lmoist)  call stats_allocate_interp_and_sgs_moist
      end if
      if(ltdump .or. lytdump .or. lydump) then
        if (nsv>0)   call stats_allocate_interp_and_sgs_scalar
      end if

      if(ltdump .or. lxytdump .or. lytdump) then
        !> allocate variables to compute time-averaged quantities
        call stats_allocate_tavg_vel
        if (ltempeq) call stats_allocate_tavg_temp
        if (lmoist)  call stats_allocate_tavg_moist
      end if
      if(ltdump .or. lytdump) then
        if (nsv>0)   call stats_allocate_tavg_scalar
      end if

      !> Generate time averaged NetCDF: stats_t_out.xxx.xxx.xxx.nc
      if (ltdump) then
        !> Total numbers of variables to be written
        tVarsCount = 14
        if (ltempeq) tVarsCount = tVarsCount + 4
        if (lmoist)  tVarsCount = tVarsCount + 4
        if (nsv>0)   tVarsCount = tVarsCount + 4*nsv
        if ((lchem) .and. (nsv>2)) tVarsCount = tVarsCount + 1

        allocate(tVars(tVarsCount,4))   !!> Array to store the variable description of the quantities to be written
        ctrt = 0
        call stats_ncdescription_tavg_vel
        if (ltempeq) call stats_ncdescription_tavg_temp
        if (lmoist)  call stats_ncdescription_tavg_moist
        if (nsv>0)   call stats_ncdescription_tavg_scalar
        if ((lchem) .and. (nsv>2)) call stats_init_tavg_PSS
        
        call stats_createnc_tavg
        
        deallocate(tVars)
      end if

      !> Generate time, y and x averaged NetCDF: stats_xyt_out.xxx.nc
      if (lxytdump) then
        xytVarsCount = 20
        if (ltempeq) xytVarsCount = xytVarsCount + 5
        if (lmoist)  xytVarsCount = xytVarsCount + 5

        allocate(xytVars(xytVarsCount,4))   !!> Array to store the variable description of the quantities to be written
        ctrxyt = 0
        call stats_init_xytavg_vel
        if (ltempeq) call stats_init_xytavg_temp
        if (lmoist)  call stats_init_xytavg_moist

        call stats_createnc_xytavg
        
        deallocate(xytVars)
      end if

      !> Generate y and x averaged NetCDF: stats_xy_out.xxx.nc
      if (lxydump) then
        xyVarsCount = 16
        if (ltempeq) xyVarsCount = xyVarsCount + 4
        if (lmoist)  xyVarsCount = xyVarsCount + 4

        allocate(xyVars(xyVarsCount,4))   !!> Array to store the variable description of the quantities to be written
        ctrxy = 0
        call stats_init_xyavg_vel
        if (ltempeq) call stats_init_xyavg_temp
        if (lmoist)  call stats_init_xyavg_moist

        call stats_createnc_xyavg
        
        deallocate(xyVars)
      end if

      !> Generate time and y averaged NetCDF: stats_yt_out.xxx.xxx.nc
      if (lytdump) then
        ytVarsCount = 11
        if (ltempeq) ytVarsCount = ytVarsCount + 5
        if (lmoist)  ytVarsCount = ytVarsCount + 5
        if (nsv>0)   ytVarsCount = ytVarsCount + 5*nsv

        allocate(ytVars(ytVarsCount,4))   !!> Array to store the variable description of the quantities to be written
        ctryt = 0
        call stats_init_ytavg_vel
        if (ltempeq) call stats_init_ytavg_temp
        if (lmoist)  call stats_init_ytavg_moist
        if (nsv>0)   call stats_init_ytavg_scalar

        call stats_createnc_ytavg
        
        deallocate(ytVars)
      end if

      !> Generate y averaged NetCDF: stats_y_out.xxx.xxx.nc
      if (lydump) then
        yVarsCount = 8
        if (ltempeq) yVarsCount = yVarsCount + 4
        if (lmoist)  yVarsCount = yVarsCount + 4
        if (nsv>0)   yVarsCount = yVarsCount + 4*nsv

        allocate(yVars(yVarsCount,4))   !!> Array to store the variable description of the quantities to be written
        ctry = 0
        call stats_init_yavg_vel
        if (ltempeq) call stats_init_yavg_temp
        if (lmoist)  call stats_init_yavg_moist
        if (nsv>0)   call stats_init_yavg_scalar

        call stats_createnc_yavg

        deallocate(yVars)
      end if

      !> Generate time averaged tree data NetCDF: stats_tree_out.xxx.xxx.xxx.nc
      if (ltreedump) then
        treeVarsCount = 3
        if (ltempeq) treeVarsCount = treeVarsCount + 1
        if (lmoist)  treeVarsCount = treeVarsCount + 4
        if (nsv>0)   treeVarsCount = treeVarsCount + nsv
        allocate(treeVars(treeVarsCount,4))   !!> Array to store the variable description of the quantities to be written
        ctrtree = 0
        call stats_init_tree_vel
        if (ltempeq) call stats_init_tree_temp
        if (lmoist)  call stats_init_tree_moist
        if (nsv>0)   call stats_init_tree_scalar

        call stats_createnc_tree
        
        deallocate(treeVars)
      end if
    end subroutine stats_init


    subroutine stats_main
      implicit none

      if (.not. rk3step==3)  return
      if(.not.(ltdump .or. lxytdump .or. lxydump .or. lytdump .or. lydump .or. ltreedump)) return

      if (tsamplep > tsample) then        ! at every stats sampling instance
        tstatsdumppi = 1./tstatsdumpp

        if(ltdump .or. lxytdump .or. lxydump .or. lytdump .or. lydump) then
          call stats_interpolate_and_sgs_vel
          if (ltempeq) call stats_interpolate_and_sgs_temp
          if (lmoist)  call stats_interpolate_and_sgs_moist
        end if
        if(ltdump .or. lytdump .or. lydump) then
          if (nsv>0)   call stats_interpolate_and_sgs_scalar
        end if

        if(ltdump .or. lxytdump .or. lytdump) then
          call stats_compute_tavg_vel
          if (ltempeq) call stats_compute_tavg_temp
          if (lmoist)  call stats_compute_tavg_moist
        end if
        if(ltdump .or. lytdump) then
          if (nsv>0) call stats_compute_tavg_scalar
        end if
        if(ltdump) then
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

        if(lydump) then
          call stats_compute_yavg_vel
          if (ltempeq) call stats_compute_yavg_temp
          if (lmoist)  call stats_compute_yavg_moist
          if (nsv>0)   call stats_compute_yavg_scalar
          if (myidy==0) then
            call writestat_nc(ncidy, 'time', timee, nrecy, .true.)
            call stats_write_yavg_vel
            if (ltempeq) call stats_write_yavg_temp
            if (lmoist)  call stats_write_yavg_moist
            if (nsv>0)   call stats_write_yavg_scalar
          end if
        end if

        if(ltreedump) then
          call stats_compute_tree_vel
          if (ltempeq) call stats_compute_tree_temp
          if (lmoist)  call stats_compute_tree_moist
          if (nsv>0)   call stats_compute_tree_scalar
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

        if(lytdump) then
          call stats_compute_ytavg_vel
          if (ltempeq) call stats_compute_ytavg_temp
          if (lmoist)  call stats_compute_ytavg_moist
          if (nsv>0)   call stats_compute_ytavg_scalar
          if (myidy==0) then
            call writestat_nc(ncidyt, 'time', timee, nrecyt, .true.)
            call stats_write_ytavg_vel
            if (ltempeq) call stats_write_ytavg_temp
            if (lmoist)  call stats_write_ytavg_moist
            if (nsv>0)   call stats_write_ytavg_scalar
          end if
        end if

        if(ltreedump) then
          call writestat_nc(ncidtree, 'time', timee, nrectree, .true.)
          call stats_write_tree_vel
          if (ltempeq) call stats_write_tree_temp
          if (lmoist)  call stats_write_tree_moist
          if (nsv>0)   call stats_write_tree_scalar
        end if

        tstatsdumpp = dt
      else
        tstatsdumpp = tstatsdumpp + dt
      endif

    end subroutine stats_main


    !! ## %% Interpolated and sgs fields initialization routines
    subroutine stats_allocate_interp_and_sgs_vel
      implicit none
      allocate(uik(ib:ie,jb:je,kb:ke+kh))
      allocate(wik(ib:ie,jb:je,kb:ke+kh))
      allocate(vjk(ib:ie,jb:je,kb:ke+kh))
      allocate(wjk(ib:ie,jb:je,kb:ke+kh))
      allocate(uij(ib:ie,jb:je,kb:ke+kh))
      allocate(vij(ib:ie,jb:je,kb:ke+kh))
      allocate(uc(ib:ie,jb:je,kb:ke+kh))
      allocate(vc(ib:ie,jb:je,kb:ke+kh))
      allocate(wc(ib:ie,jb:je,kb:ke+kh))
      allocate(usgs(ib:ie,jb:je,kb:ke+kh))
      allocate(vsgs(ib:ie,jb:je,kb:ke+kh))
      allocate(wsgs(ib:ie,jb:je,kb:ke+kh))
    end subroutine stats_allocate_interp_and_sgs_vel

    subroutine stats_allocate_interp_and_sgs_temp
      implicit none
      allocate(thlk(ib:ie,jb:je,kb:ke+kh))
      allocate(thlsgs(ib:ie,jb:je,kb:ke+kh))
    end subroutine stats_allocate_interp_and_sgs_temp

    subroutine stats_allocate_interp_and_sgs_moist
      implicit none
      allocate(qtk(ib:ie,jb:je,kb:ke+kh))
      allocate(qtsgs(ib:ie,jb:je,kb:ke+kh))
    end subroutine stats_allocate_interp_and_sgs_moist

    subroutine stats_allocate_interp_and_sgs_scalar
      implicit none
      allocate(svk(ib:ie,jb:je,kb:ke+kh,nsv))
      allocate(svsgs(ib:ie,jb:je,kb:ke+kh,nsv))
    end subroutine stats_allocate_interp_and_sgs_scalar
      

    !! ## %% Time averaging initialization routines
    subroutine stats_allocate_tavg_vel
      implicit none
      allocate(ut(ib:ie,jb:je,kb:ke+kh))   ; ut    = 0;
      allocate(vt(ib:ie,jb:je,kb:ke+kh))   ; vt    = 0;
      allocate(wt(ib:ie,jb:je,kb:ke+kh))   ; wt    = 0;
      allocate(pt(ib:ie,jb:je,kb:ke+kh))   ; pt    = 0;
      allocate(utik(ib:ie,jb:je,kb:ke+kh)) ; utik  = 0;
      allocate(wtik(ib:ie,jb:je,kb:ke+kh)) ; wtik  = 0;
      allocate(uwtik(ib:ie,jb:je,kb:ke+kh)); uwtik = 0;
      allocate(vtjk(ib:ie,jb:je,kb:ke+kh)) ; vtjk  = 0;
      allocate(wtjk(ib:ie,jb:je,kb:ke+kh)) ; wtjk  = 0;
      allocate(vwtjk(ib:ie,jb:je,kb:ke+kh)); vwtjk = 0;
      allocate(utij(ib:ie,jb:je,kb:ke+kh)) ; utij  = 0;
      allocate(vtij(ib:ie,jb:je,kb:ke+kh)) ; vtij  = 0;
      allocate(uvtij(ib:ie,jb:je,kb:ke+kh)); uvtij = 0;
      allocate(utc(ib:ie,jb:je,kb:ke+kh))  ; utc   = 0;
      allocate(vtc(ib:ie,jb:je,kb:ke+kh))  ; vtc   = 0;
      allocate(wtc(ib:ie,jb:je,kb:ke+kh))  ; wtc   = 0;
      allocate(uutc(ib:ie,jb:je,kb:ke+kh)) ; uutc  = 0;
      allocate(vvtc(ib:ie,jb:je,kb:ke+kh)) ; vvtc  = 0;
      allocate(wwtc(ib:ie,jb:je,kb:ke+kh)) ; wwtc  = 0;      
      allocate(usgst(ib:ie,jb:je,kb:ke+kh)); usgst = 0;
      allocate(vsgst(ib:ie,jb:je,kb:ke+kh)); vsgst = 0;
      allocate(wsgst(ib:ie,jb:je,kb:ke+kh)); wsgst = 0;
    end subroutine stats_allocate_tavg_vel
    subroutine stats_ncdescription_tavg_vel
      implicit none
      !> Generate variable description for the quantities to be written in the time averaged NetCDF: stats.xxx.xxx.xxx.nc
      call ncinfo( tVars(ctrt+ 1,:), 'u'       , 'Streamwise velocity'        , 'm/s'       , 'mttt' )
      call ncinfo( tVars(ctrt+ 2,:), 'v'       , 'Spanwise velocity'          , 'm/s'       , 'tmtt' )
      call ncinfo( tVars(ctrt+ 3,:), 'w'       , 'Vertical velocity'          , 'm/s'       , 'ttmt' )
      call ncinfo( tVars(ctrt+ 4,:), 'p'       , 'Kinematic Pressure'         , 'm^2/s^2'   , 'tttt' )

      call ncinfo( tVars(ctrt+ 5,:), 'upwp'    , 'Turbulent momentum flux'    , 'm^2/s^2'   , 'mtmt' )
      call ncinfo( tVars(ctrt+ 6,:), 'vpwp'    , 'Turbulent momentum flux'    , 'm^2/s^2'   , 'tmmt' )
      call ncinfo( tVars(ctrt+ 7,:), 'upvp'    , 'Turbulent momentum flux'    , 'm^2/s^2'   , 'mmtt' )

      call ncinfo( tVars(ctrt+ 8,:), 'upup'    , 'u variance - cell centered' , 'm^2/s^2'   , 'tttt' )
      call ncinfo( tVars(ctrt+ 9,:), 'vpvp'    , 'v variance - cell centered' , 'm^2/s^2'   , 'tttt' )
      call ncinfo( tVars(ctrt+10,:), 'wpwp'    , 'w variance - cell centered' , 'm^2/s^2'   , 'tttt' )
      call ncinfo( tVars(ctrt+11,:), 'tke'     , 'TKE - cell centered'        , 'm^2/s^2'   , 'tttt' )

      call ncinfo( tVars(ctrt+12,:), 'usgs'    , 'SGS u flux'                 , 'm^2/s^2'   , 'mtmt' )
      call ncinfo( tVars(ctrt+13,:), 'vsgs'    , 'SGS v flux'                 , 'm^2/s^2'   , 'tmmt' )
      call ncinfo( tVars(ctrt+14,:), 'wsgs'    , 'SGS w flux'                 , 'm^2/s^2'   , 'ttmt' )
      ctrt = ctrt+14
    end subroutine stats_ncdescription_tavg_vel

    subroutine stats_allocate_tavg_temp
      implicit none
      allocate(thlt(ib:ie,jb:je,kb:ke+kh))   ; thlt    = 0;
      allocate(thltk(ib:ie,jb:je,kb:ke+kh))  ; thltk   = 0;
      allocate(wthltk(ib:ie,jb:je,kb:ke+kh)) ; wthltk  = 0;
      allocate(thlthlt(ib:ie,jb:je,kb:ke+kh)); thlthlt = 0;
      allocate(thlsgst(ib:ie,jb:je,kb:ke+kh)); thlsgst = 0;
    end subroutine stats_allocate_tavg_temp
    subroutine stats_ncdescription_tavg_temp
      implicit none
      call ncinfo( tVars(ctrt+1,:) , 'thl'     , 'Temperature'               , 'K'         , 'tttt' )
      call ncinfo( tVars(ctrt+2,:) , 'wpthlp'  , 'Turbulent heat flux'       , 'K m/s'     , 'ttmt' )
      call ncinfo( tVars(ctrt+3,:) , 'thlpthlp', 'Temperature variance'      , 'K^2'       , 'tttt' )
      call ncinfo( tVars(ctrt+4,:) , 'thlsgs'  , 'SGS temperature flux'      , 'K m/s'     , 'ttmt' )
      ctrt = ctrt+4
    end subroutine stats_ncdescription_tavg_temp

    subroutine stats_allocate_tavg_moist
      implicit none
      allocate(qtt(ib:ie,jb:je,kb:ke+kh))   ; qtt    = 0;
      allocate(qttk(ib:ie,jb:je,kb:ke+kh))  ; qttk   = 0;
      allocate(wqttk(ib:ie,jb:je,kb:ke+kh)) ; wqttk  = 0;
      allocate(qtqtt(ib:ie,jb:je,kb:ke+kh)) ; qtqtt  = 0;
      allocate(qtsgst(ib:ie,jb:je,kb:ke+kh)); qtsgst = 0;
    end subroutine stats_allocate_tavg_moist
    subroutine stats_ncdescription_tavg_moist
      implicit none
      call ncinfo( tVars(ctrt+1,:) , 'qt'      , 'Moisture'                  , 'kg/kg'     , 'tttt' )
      call ncinfo( tVars(ctrt+2,:) , 'wpqtp'   , 'Turbulent moisture flux'   , 'kg m/kg s' , 'ttmt' )
      call ncinfo( tVars(ctrt+3,:) , 'qtpqtp'  , 'Moisture variance'         , 'kg^2/kg^2' , 'tttt' )
      call ncinfo( tVars(ctrt+4,:) , 'qtsgs'   , 'SGS moisture flux'         , 'kg m/kg s' , 'ttmt' )
      ctrt = ctrt+4
    end subroutine stats_ncdescription_tavg_moist

    subroutine stats_allocate_tavg_scalar
      implicit none
      allocate(svt(ib:ie,jb:je,kb:ke+kh,nsv))   ; svt    = 0;
      allocate(svtk(ib:ie,jb:je,kb:ke+kh,nsv))  ; svtk   = 0;
      allocate(wsvtk(ib:ie,jb:je,kb:ke+kh,nsv)) ; wsvtk  = 0;
      allocate(svsvt(ib:ie,jb:je,kb:ke+kh,nsv)) ; svsvt  = 0;
      allocate(svsgst(ib:ie,jb:je,kb:ke+kh,nsv)); svsgst = 0;
    end subroutine stats_allocate_tavg_scalar
    subroutine stats_ncdescription_tavg_scalar
      implicit none
      integer :: n
      character(2) :: sid
      allocate(svtname(nsv))
      allocate(wpsvptname(nsv))
      allocate(svpsvptname(nsv))
      allocate(svsgsname(nsv))
      do n = 1, nsv
        write (sid, '(I0)') n
        svtname(n)     = 's'//trim(sid)                        ! s1       at n = 1
        wpsvptname(n)  = 'wps'//trim(sid)//'p'                 ! wps1p    at n = 1
        svpsvptname(n) = 's'//trim(sid)//'ps'//trim(sid)//'p'  ! s1ps1p   at n = 1
        svsgsname(n)   = 's'//trim(sid)//'sgs'                 ! s1sgs    at n = 1
        call ncinfo(tVars(ctrt+n,:)      , trim(svtname(n))    , 'Concentration field '//trim(sid)   , 'g/m^3'  , 'tttt' )
        call ncinfo(tVars(ctrt+nsv+n,:)  , trim(wpsvptname(n)) , 'Turbulent scalar flux '//trim(sid) , 'g/m^2s' , 'ttmt' )
        call ncinfo(tVars(ctrt+2*nsv+n,:), trim(svpsvptname(n)), 'Concentration variance '//trim(sid), 'g^2/m^6', 'tttt' )
        call ncinfo(tVars(ctrt+3*nsv+n,:), trim(svsgsname(n))  , 'SGS scalar flux '//trim(sid)       , 'g/m^2s' , 'ttmt' )
      end do
      ctrt = ctrt+4*nsv
    end subroutine stats_ncdescription_tavg_scalar

    subroutine stats_init_tavg_PSS
      implicit none
      allocate(PSS(ib:ie,jb:je,kb:ke+kh))
      allocate(PSSt(ib:ie,jb:je,kb:ke+kh)); PSSt = 0;
      call ncinfo( tVars(ctrt+1,:) , 'PSS'      , 'PSS defect'                , 'gm/s'      , 'tttt' )
      ctrt = ctrt+1
    end subroutine stats_init_tavg_PSS

    subroutine stats_createnc_tavg
      implicit none
      filenamet = 'stats_t_out.xxx.xxx.xxx.nc'
      filenamet(13:15) = cmyidx
      filenamet(17:19) = cmyidy
      filenamet(21:23) = cexpnr

      nrect = 0
      call open_nc(filenamet, ncidt, nrect, n1=xdim, n2=ydim, n3=zdim)
      if (nrect==0) then
        call define_nc(ncidt, 1, timeVar)
        call writestat_dims_nc(ncidt)
      end if
      call define_nc(ncidt, tVarsCount, tVars)
    end subroutine stats_createnc_tavg


    !! ## %% Time, y and x averaging initialization routines
    subroutine stats_init_xytavg_vel
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

      call ncinfo( xytVars(ctrxyt+ 1,:), 'u'        , 'Streamwise velocity'        , 'm/s'       , 'tt' )
      call ncinfo( xytVars(ctrxyt+ 2,:), 'v'        , 'Spanwise velocity'          , 'm/s'       , 'tt' )
      call ncinfo( xytVars(ctrxyt+ 3,:), 'w'        , 'Vertical velocity'          , 'm/s'       , 'mt' )
      call ncinfo( xytVars(ctrxyt+ 4,:), 'p'        , 'Kinematic Pressure'         , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+ 5,:), 'upwp'     , 'Turbulent mom. flux'        , 'm^2/s^2'   , 'mt' )
      call ncinfo( xytVars(ctrxyt+ 6,:), 'vpwp'     , 'Turbulent mom. flux'        , 'm^2/s^2'   , 'mt' )
      call ncinfo( xytVars(ctrxyt+ 7,:), 'upvp'     , 'Turbulent mom. flux'        , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+ 8,:), 'uw'       , 'Dispersive mom. flux'       , 'm^2/s^2'   , 'mt' )
      call ncinfo( xytVars(ctrxyt+ 9,:), 'vw'       , 'Dispersive mom. flux'       , 'm^2/s^2'   , 'mt' )
      call ncinfo( xytVars(ctrxyt+10,:), 'uv'       , 'Dispersive mom. flux'       , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+11,:), 'uu'       , 'Dispersive mom. flux'       , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+12,:), 'vv'       , 'Dispersive mom. flux'       , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+13,:), 'ww'       , 'Dispersive mom. flux'       , 'm^2/s^2'   , 'mt' )
      call ncinfo( xytVars(ctrxyt+14,:), 'upup'     , 'u variance - cell centered' , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+15,:), 'vpvp'     , 'v variance - cell centered' , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+16,:), 'wpwp'     , 'w variance - cell centered' , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+17,:), 'tke'      , 'TKE - cell centered'        , 'm^2/s^2'   , 'tt' )
      call ncinfo( xytVars(ctrxyt+18,:), 'usgs'     , 'SGS mom. flux'              , 'm^2/s^2'   , 'mt' )
      call ncinfo( xytVars(ctrxyt+19,:), 'vsgs'     , 'SGS mom. flux'              , 'm^2/s^2'   , 'mt' )
      call ncinfo( xytVars(ctrxyt+20,:), 'wsgs'     , 'SGS mom. flux'              , 'm^2/s^2'   , 'mt' )
      ctrxyt = ctrxyt+20
    end subroutine stats_init_xytavg_vel

    subroutine stats_init_xytavg_temp
      implicit none
      allocate(thlxyt(kb:ke+kh))
      allocate(wpthlpxytk(kb:ke+kh))
      allocate(wthlxytk(kb:ke+kh))
      allocate(thlpthlpxyt(kb:ke+kh))
      allocate(thlsgsxyt(kb:ke+kh))
      call ncinfo( xytVars(ctrxyt+ 1,:), 'thl'      , 'Temperature'              , 'K'         , 'tt' )
      call ncinfo( xytVars(ctrxyt+ 2,:), 'wpthlp'   , 'Turbulent heat flux'      , 'K m/s'     , 'mt' )
      call ncinfo( xytVars(ctrxyt+ 3,:), 'wthl'     , 'Dispersive heat flux'     , 'K m/s'     , 'mt' )
      call ncinfo( xytVars(ctrxyt+ 4,:), 'thlpthlp' , 'Temp. variance'           , 'K^2'       , 'tt' )
      call ncinfo( xytVars(ctrxyt+ 5,:), 'thlsgs'   , 'SGS heat flux'            , 'K m/s'     , 'mt' )
      ctrxyt = ctrxyt+5
    end subroutine stats_init_xytavg_temp

    subroutine stats_init_xytavg_moist
      implicit none
      allocate(qtxyt(kb:ke+kh))
      allocate(wpqtpxytk(kb:ke+kh))
      allocate(wqtxytk(kb:ke+kh))
      allocate(qtpqtpxyt(kb:ke+kh))
      allocate(qtsgsxyt(kb:ke+kh))
      call ncinfo( xytVars(ctrxyt+ 1,:), 'qt'       , 'Moisture'                 , 'kg/kg'     , 'tt' )
      call ncinfo( xytVars(ctrxyt+ 2,:), 'wpqtp'    , 'Turbulent moisture flux'  , 'kg m/kg s' , 'mt' )
      call ncinfo( xytVars(ctrxyt+ 3,:), 'wqt'      , 'Dispersive moisture flux' , 'kg m/kg s' , 'mt' )
      call ncinfo( xytVars(ctrxyt+ 4,:), 'qtpqtp'   , 'Moisture variance'        , 'kg^2/kg^2' , 'tt' )
      call ncinfo( xytVars(ctrxyt+ 5,:), 'qtsgs'    , 'SGS moisture flux'        , 'kg m/kg s' , 'mt' )
      ctrxyt = ctrxyt+5
    end subroutine stats_init_xytavg_moist

    subroutine stats_createnc_xytavg
      implicit none
      filenamexyt = 'stats_xyt_out.xxx.nc'
      filenamexyt(15:17) = cexpnr

      nrecxyt = 0
      if (myid==0) then
        call open_nc(filenamexyt, ncidxyt, nrecxyt, n3=zdim)
        if (nrecxyt==0) then
          call define_nc(ncidxyt, 1,  timeVar)
          call writestat_dims_nc(ncidxyt)
        end if
        call define_nc(ncidxyt, xytVarsCount, xytVars)
      end if
    end subroutine stats_createnc_xytavg


    !! ## %% y and x averaging initialization routines
    subroutine stats_init_xyavg_vel
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

      call ncinfo( xyVars(ctrxy+ 1,:), 'u'        , 'Streamwise velocity'      , 'm/s'       , 'tt' )
      call ncinfo( xyVars(ctrxy+ 2,:), 'v'        , 'Spanwise velocity'        , 'm/s'       , 'tt' )
      call ncinfo( xyVars(ctrxy+ 3,:), 'w'        , 'Vertical velocity'        , 'm/s'       , 'mt' )
      call ncinfo( xyVars(ctrxy+ 4,:), 'p'        , 'Kinematic Pressure'       , 'm^2/s^2'   , 'tt' )
      call ncinfo( xyVars(ctrxy+ 5,:), 'upwp'     , 'Turbulent mom. flux'      , 'm^2/s^2'   , 'mt' )
      call ncinfo( xyVars(ctrxy+ 6,:), 'vpwp'     , 'Turbulent mom. flux'      , 'm^2/s^2'   , 'mt' )
      call ncinfo( xyVars(ctrxy+ 7,:), 'upvp'     , 'Turbulent mom. flux'      , 'm^2/s^2'   , 'tt' )
      call ncinfo( xyVars(ctrxy+ 8,:), 'uw'       , 'Advective mom. flux'      , 'm^2/s^2'   , 'mt' )
      call ncinfo( xyVars(ctrxy+ 9,:), 'vw'       , 'Advective mom. flux'      , 'm^2/s^2'   , 'mt' )
      call ncinfo( xyVars(ctrxy+10,:), 'uv'       , 'Advective mom. flux'      , 'm^2/s^2'   , 'tt' )
      call ncinfo( xyVars(ctrxy+11,:), 'uu'       , 'Advective mom. flux'      , 'm^2/s^2'   , 'tt' )
      call ncinfo( xyVars(ctrxy+12,:), 'vv'       , 'Advective mom. flux'      , 'm^2/s^2'   , 'tt' )
      call ncinfo( xyVars(ctrxy+13,:), 'ww'       , 'Advective mom. flux'      , 'm^2/s^2'   , 'mt' )
      call ncinfo( xyVars(ctrxy+14,:), 'usgs'     , 'SGS mom. flux'            , 'm^2/s^2'   , 'mt' )
      call ncinfo( xyVars(ctrxy+15,:), 'vsgs'     , 'SGS mom. flux'            , 'm^2/s^2'   , 'mt' )
      call ncinfo( xyVars(ctrxy+16,:), 'wsgs'     , 'SGS mom. flux'            , 'm^2/s^2'   , 'mt' )
      ctrxy = ctrxy+16
    end subroutine stats_init_xyavg_vel

    subroutine stats_init_xyavg_temp
      implicit none
      allocate(thlxy(kb:ke+kh))
      allocate(wpthlpxyk(kb:ke+kh))
      allocate(wthlxyk(kb:ke+kh))
      allocate(thlxyk(kb:ke+kh))
      allocate(thlsgsxy(kb:ke+kh))
      call ncinfo( xyVars(ctrxy+ 1,:), 'thl'      , 'Temperature'              , 'K'         , 'tt' )
      call ncinfo( xyVars(ctrxy+ 2,:), 'wpthlp'   , 'Turbulent heat flux'      , 'K m/s'     , 'mt' )
      call ncinfo( xyVars(ctrxy+ 3,:), 'wthl'     , 'Advective heat flux'      , 'K m/s'     , 'mt' )
      call ncinfo( xyVars(ctrxy+ 4,:), 'thlsgs'   , 'SGS heat flux'            , 'K m/s'     , 'mt' )
      ctrxy = ctrxy+4
    end subroutine stats_init_xyavg_temp

    subroutine stats_init_xyavg_moist
      implicit none
      allocate(qtxy(kb:ke+kh))
      allocate(wpqtpxyk(kb:ke+kh))
      allocate(wqtxyk(kb:ke+kh))
      allocate(qtxyk(kb:ke+kh))
      allocate(qtsgsxy(kb:ke+kh))
      call ncinfo( xyVars(ctrxy+ 1,:), 'qt'       , 'Moisture'                 , 'kg/kg'     , 'tt' )
      call ncinfo( xyVars(ctrxy+ 2,:), 'wpqtp'    , 'Turbulent moisture flux'  , 'kg m/kg s' , 'mt' )
      call ncinfo( xyVars(ctrxy+ 3,:), 'wqt'      , 'Advective moisture flux'  , 'kg m/kg s' , 'mt' )
      call ncinfo( xyVars(ctrxy+ 4,:), 'qtsgs'    , 'SGS moisture flux'        , 'kg m/kg s' , 'mt' )
      ctrxy = ctrxy+4
    end subroutine stats_init_xyavg_moist

    subroutine stats_createnc_xyavg
      implicit none
      filenamexy = 'stats_xy_out.xxx.nc'
      filenamexy(14:16) = cexpnr

      nrecxy = 0
      if (myid==0) then
        call open_nc(filenamexy, ncidxy, nrecxy, n3=zdim)
        if (nrecxy==0) then
          call define_nc(ncidxy, 1,  timeVar)
          call writestat_dims_nc(ncidxy)
        end if
        call define_nc(ncidxy, xyVarsCount, xyVars)
      end if
    end subroutine stats_createnc_xyavg


    !! ## %% Time and y averaging initialization routines
    subroutine stats_init_ytavg_vel
      implicit none
      allocate(uyt(ib:ie,kb:ke))
      allocate(vyt(ib:ie,kb:ke))
      allocate(wyt(ib:ie,kb:ke))
      allocate(pyt(ib:ie,kb:ke))
      allocate(upwpytik(ib:ie,kb:ke))
      allocate(uwytik(ib:ie,kb:ke))
      allocate(upupytc(ib:ie,kb:ke))
      allocate(vpvpytc(ib:ie,kb:ke))
      allocate(wpwpytc(ib:ie,kb:ke))
      allocate(usgsyt(ib:ie,kb:ke))
      allocate(wsgsyt(ib:ie,kb:ke))

      call ncinfo( ytVars(ctryt+ 1,:), 'u'        , 'Streamwise velocity'        , 'm/s'       , 'm0tt' )
      call ncinfo( ytVars(ctryt+ 2,:), 'v'        , 'Spanwise velocity'          , 'm/s'       , 't0tt' )
      call ncinfo( ytVars(ctryt+ 3,:), 'w'        , 'Vertical velocity'          , 'm/s'       , 't0mt' )
      call ncinfo( ytVars(ctryt+ 4,:), 'p'        , 'Kinematic Pressure'         , 'm^2/s^2'   , 't0tt' )
      call ncinfo( ytVars(ctryt+ 5,:), 'upwp'     , 'Turbulent mom. flux'        , 'm^2/s^2'   , 'm0mt' )
      call ncinfo( ytVars(ctryt+ 6,:), 'uw'       , 'Dispersive mom. flux'       , 'm^2/s^2'   , 'm0mt' )
      call ncinfo( ytVars(ctryt+ 7,:), 'upup'     , 'u variance - cell centered' , 'm^2/s^2'   , 't0tt' )
      call ncinfo( ytVars(ctryt+ 8,:), 'vpvp'     , 'v variance - cell centered' , 'm^2/s^2'   , 't0tt' )
      call ncinfo( ytVars(ctryt+ 9,:), 'wpwp'     , 'w variance - cell centered' , 'm^2/s^2'   , 't0tt' )
      call ncinfo( ytVars(ctryt+10,:), 'usgs'     , 'SGS mom. flux'              , 'm^2/s^2'   , 'm0mt' )
      call ncinfo( ytVars(ctryt+11,:), 'wsgs'     , 'SGS mom. flux'              , 'm^2/s^2'   , 't0mt' )
      ctryt = ctryt+11
    end subroutine stats_init_ytavg_vel

    subroutine stats_init_ytavg_temp
      implicit none
      allocate(thlyt(ib:ie,kb:ke))
      allocate(wpthlpytk(ib:ie,kb:ke))
      allocate(wthlytk(ib:ie,kb:ke))
      allocate(thlpthlpyt(ib:ie,kb:ke))
      allocate(thlsgsyt(ib:ie,kb:ke))
      call ncinfo( ytVars(ctryt+ 1,:), 'thl'      , 'Temperature'              , 'K'         , 't0tt' )
      call ncinfo( ytVars(ctryt+ 2,:), 'wpthlp'   , 'Turbulent heat flux'      , 'K m/s'     , 't0mt' )
      call ncinfo( ytVars(ctryt+ 3,:), 'wthl'     , 'Dispersive heat flux'     , 'K m/s'     , 't0mt' )
      call ncinfo( ytVars(ctryt+ 4,:), 'thlpthlp' , 'Temp. variance'           , 'K^2'       , 't0tt' )
      call ncinfo( ytVars(ctryt+ 5,:), 'thlsgs'   , 'SGS heat flux'            , 'K m/s'     , 't0mt' )
      ctryt = ctryt+5
    end subroutine stats_init_ytavg_temp

    subroutine stats_init_ytavg_moist
      implicit none
      allocate(qtyt(ib:ie,kb:ke))
      allocate(wpqtpytk(ib:ie,kb:ke))
      allocate(wqtytk(ib:ie,kb:ke))
      allocate(qtpqtpyt(ib:ie,kb:ke))
      allocate(qtsgsyt(ib:ie,kb:ke))
      call ncinfo( ytVars(ctryt+ 1,:), 'qt'       , 'Moisture'                 , 'kg/kg'     , 't0tt' )
      call ncinfo( ytVars(ctryt+ 2,:), 'wpqtp'    , 'Turbulent moisture flux'  , 'kg m/kg s' , 't0mt' )
      call ncinfo( ytVars(ctryt+ 3,:), 'wqt'      , 'Dispersive moisture flux' , 'kg m/kg s' , 't0mt' )
      call ncinfo( ytVars(ctryt+ 4,:), 'qtpqtp'   , 'Moisture variance'        , 'kg^2/kg^2' , 't0tt' )
      call ncinfo( ytVars(ctryt+ 5,:), 'qtsgs'    , 'SGS moisture flux'        , 'kg m/kg s' , 't0mt' )
      ctryt = ctryt+5
    end subroutine stats_init_ytavg_moist

    subroutine stats_init_ytavg_scalar
      integer :: n
      character(2) :: sid
      allocate(svytname(nsv))
      allocate(wpsvpytname(nsv))
      allocate(wsvytname(nsv))
      allocate(svpsvpytname(nsv))
      allocate(svsgsytname(nsv))
      allocate(svyt(ib:ie,kb:ke,nsv))
      allocate(wpsvpytk(ib:ie,kb:ke,nsv))
      allocate(wsvytk(ib:ie,kb:ke,nsv))
      allocate(svpsvpyt(ib:ie,kb:ke,nsv))
      allocate(svsgsyt(ib:ie,kb:ke,nsv))
      do n = 1, nsv
        write (sid, '(I0)') n
        svytname(n)     = 's'//trim(sid)                        ! s1       at n = 1
        wpsvpytname(n)  = 'wps'//trim(sid)//'p'                 ! wps1p    at n = 1
        wsvytname(n)    = 'ws'//trim(sid)                       ! ws1      at n = 1
        svpsvpytname(n) = 's'//trim(sid)//'ps'//trim(sid)//'p'  ! s1ps1p   at n = 1
        svsgsytname(n)  = 's'//trim(sid)//'sgs'                 ! s1sgs    at n = 1
        call ncinfo(ytVars(ctryt+n,:)      , trim(svytname(n))    , 'Concentration field '//trim(sid)   , 'g/m^3'  , 't0tt' )
        call ncinfo(ytVars(ctryt+nsv+n,:)  , trim(wpsvpytname(n)) , 'Turbulent scalar flux '//trim(sid) , 'g/m^2s' , 't0mt' )
        call ncinfo(ytVars(ctryt+2*nsv+n,:), trim(wsvytname(n))   , 'Dispersive scalar flux '//trim(sid), 'g/m^2s' , 't0mt' )
        call ncinfo(ytVars(ctryt+3*nsv+n,:), trim(svpsvpytname(n)), 'Concentration variance '//trim(sid), 'g^2/m^6', 't0tt' )
        call ncinfo(ytVars(ctryt+4*nsv+n,:), trim(svsgsytname(n)) , 'SGS scalar flux '//trim(sid)       , 'g/m^2s' , 't0mt' )
      end do
      ctryt = ctryt+5*nsv
    end subroutine stats_init_ytavg_scalar

    subroutine stats_createnc_ytavg
      implicit none
      filenameyt = 'stats_yt_out.xxx.xxx.nc'
      filenameyt(14:16)  = cmyidx
      filenameyt(18:20) = cexpnr

      nrecyt = 0
      if (myidy==0) then
        call open_nc(filenameyt, ncidyt, nrecyt, n1=xdim, n3=zdim)
        if (nrecyt==0) then
          call define_nc(ncidyt, 1,  timeVar)
          call writestat_dims_nc(ncidyt)
        end if
        call define_nc(ncidyt, ytVarsCount, ytVars)
      end if
    end subroutine stats_createnc_ytavg


    !! ## %% y averaging initialization routines
    subroutine stats_init_yavg_vel
      implicit none
      allocate(uy(ib:ie,kb:ke))
      allocate(vy(ib:ie,kb:ke))
      allocate(wy(ib:ie,kb:ke))
      allocate(py(ib:ie,kb:ke))
      allocate(upwpyik(ib:ie,kb:ke))
      allocate(uwyik(ib:ie,kb:ke))
      allocate(uyik(ib:ie,kb:ke))
      allocate(wyik(ib:ie,kb:ke))
      allocate(usgsy(ib:ie,kb:ke))
      allocate(wsgsy(ib:ie,kb:ke))
      call ncinfo( yVars(ctry+ 1,:), 'u'        , 'Streamwise velocity'      , 'm/s'       , 'm0tt' )
      call ncinfo( yVars(ctry+ 2,:), 'v'        , 'Spanwise velocity'        , 'm/s'       , 't0tt' )
      call ncinfo( yVars(ctry+ 3,:), 'w'        , 'Vertical velocity'        , 'm/s'       , 't0mt' )
      call ncinfo( yVars(ctry+ 4,:), 'p'        , 'Kinematic Pressure'       , 'm^2/s^2'   , 't0tt' )
      call ncinfo( yVars(ctry+ 5,:), 'upwp'     , 'Turbulent mom. flux'      , 'm^2/s^2'   , 'm0mt' )
      call ncinfo( yVars(ctry+ 6,:), 'uw'       , 'Advective mom. flux'      , 'm^2/s^2'   , 'm0mt' )
      call ncinfo( yVars(ctry+ 7,:), 'usgs'     , 'SGS mom. flux'            , 'm^2/s^2'   , 'm0mt' )
      call ncinfo( yVars(ctry+ 8,:), 'wsgs'     , 'SGS mom. flux'            , 'm^2/s^2'   , 't0mt' )
      ctry = ctry+8
    end subroutine stats_init_yavg_vel

    subroutine stats_init_yavg_temp
      implicit none
      allocate(thly(ib:ie,kb:ke))
      allocate(wpthlpyk(ib:ie,kb:ke))
      allocate(wthlyk(ib:ie,kb:ke))
      allocate(thlyk(ib:ie,kb:ke))
      allocate(thlsgsy(ib:ie,kb:ke))
      call ncinfo( yVars(ctry+ 1,:), 'thl'      , 'Temperature'              , 'K'         , 't0tt' )
      call ncinfo( yVars(ctry+ 2,:), 'wpthlp'   , 'Turbulent heat flux'      , 'K m/s'     , 't0mt' )
      call ncinfo( yVars(ctry+ 3,:), 'wthl'     , 'Advective heat flux'      , 'K m/s'     , 't0mt' )
      call ncinfo( yVars(ctry+ 4,:), 'thlsgs'   , 'SGS heat flux'            , 'K m/s'     , 't0mt' )
      ctry = ctry+4
    end subroutine stats_init_yavg_temp

    subroutine stats_init_yavg_moist
      implicit none
      allocate(qty(ib:ie,kb:ke))
      allocate(wpqtpyk(ib:ie,kb:ke))
      allocate(wqtyk(ib:ie,kb:ke))
      allocate(qtyk(ib:ie,kb:ke))
      allocate(qtsgsy(ib:ie,kb:ke))
      call ncinfo( yVars(ctry+ 1,:), 'qt'       , 'Moisture'                 , 'kg/kg'     , 't0tt' )
      call ncinfo( yVars(ctry+ 2,:), 'wpqtp'    , 'Turbulent moisture flux'  , 'kg m/kg s' , 't0mt' )
      call ncinfo( yVars(ctry+ 3,:), 'wqt'      , 'Advective moisture flux'  , 'kg m/kg s' , 't0mt' )
      call ncinfo( yVars(ctry+ 4,:), 'qtsgs'    , 'SGS moisture flux'        , 'kg m/kg s' , 't0mt' )
      ctry = ctry+4
    end subroutine stats_init_yavg_moist

    subroutine stats_init_yavg_scalar
      integer :: n
      character(2) :: sid
      allocate(svyname(nsv))
      allocate(wpsvpyname(nsv))
      allocate(wsvyname(nsv))
      allocate(svsgsyname(nsv))
      allocate(svy(ib:ie,kb:ke,nsv))
      allocate(wpsvpyk(ib:ie,kb:ke,nsv))
      allocate(wsvyk(ib:ie,kb:ke,nsv))
      allocate(svyk(ib:ie,kb:ke,nsv))
      allocate(svsgsy(ib:ie,kb:ke,nsv))
      do n = 1, nsv
        write (sid, '(I0)') n
        svyname(n)     = 's'//trim(sid)                      ! s1       at n = 1
        wpsvpyname(n)  = 'wps'//trim(sid)//'p'               ! wps1p    at n = 1
        wsvyname(n)    = 'ws'//trim(sid)                     ! ws1      at n = 1
        svsgsyname(n)  = 's'//trim(sid)//'sgs'               ! s1sgs    at n = 1
        call ncinfo(yVars(ctry+n,:)      , trim(svyname(n))    , 'Concentration field '//trim(sid)   , 'g/m^3'  , 't0tt' )
        call ncinfo(yVars(ctry+nsv+n,:)  , trim(wpsvpyname(n)) , 'Turbulent scalar flux '//trim(sid) , 'g/m^2s' , 't0mt' )
        call ncinfo(yVars(ctry+2*nsv+n,:), trim(wsvyname(n))   , 'Advective scalar flux '//trim(sid) , 'g/m^2s' , 't0mt' )
        call ncinfo(yVars(ctry+3*nsv+n,:), trim(svsgsyname(n)) , 'SGS scalar flux '//trim(sid)       , 'g/m^2s' , 't0mt' )
      end do
      ctry = ctry+4*nsv
    end subroutine stats_init_yavg_scalar

    subroutine stats_createnc_yavg
      implicit none
      filenamey = 'stats_y_out.xxx.xxx.nc'
      filenamey(13:15)  = cmyidx
      filenamey(17:19) = cexpnr
      
      nrecy = 0
      if (myidy==0) then
        call open_nc(filenamey, ncidy, nrecy, n1=xdim, n3=zdim)
        if (nrecy==0) then
          call define_nc(ncidy, 1,  timeVar)
          call writestat_dims_nc(ncidy)
        end if
        call define_nc(ncidy, yVarsCount, yVars)
      end if
    end subroutine stats_createnc_yavg


    !! ## %% time averaging tree data initialization routines
    subroutine stats_init_tree_vel
      implicit none
      allocate(tr_ut(ib:ie,jb:je,kb:ke))     ; tr_ut     = 0;
      allocate(tr_vt(ib:ie,jb:je,kb:ke))     ; tr_vt     = 0;
      allocate(tr_wt(ib:ie,jb:je,kb:ke))     ; tr_wt     = 0;
      call ncinfo( treeVars(ctrtree+ 1,:), 'tr_u'      , 'Drag in x'            , 'm/s^2'   , 'tttt' )
      call ncinfo( treeVars(ctrtree+ 2,:), 'tr_v'      , 'Drag in y'            , 'm/s^2'   , 'tttt' )
      call ncinfo( treeVars(ctrtree+ 3,:), 'tr_w'      , 'Drag in z'            , 'm/s^2'   , 'ttmt' )
      ctrtree = ctrtree+3
    end subroutine stats_init_tree_vel

    subroutine stats_init_tree_temp
      implicit none
      allocate(tr_thlt(ib:ie,jb:je,kb:ke))   ; tr_thlt   = 0;
      call ncinfo( treeVars(ctrtree+ 1,:), 'tr_thl'    , 'Temp source/ sink'    , 'K/s'     , 'tttt' )
      ctrtree = ctrtree+1
    end subroutine stats_init_tree_temp

    subroutine stats_init_tree_moist
      implicit none
      allocate(tr_qtt(ib:ie,jb:je,kb:ke))    ; tr_qtt    = 0;
      allocate(tr_qtRt(ib:ie,jb:je,kb:ke))   ; tr_qtRt   = 0;
      allocate(tr_qtAt(ib:ie,jb:je,kb:ke))   ; tr_qtAt   = 0;
      allocate(tr_omegat(ib:ie,jb:je,kb:ke)) ; tr_omegat = 0;
      call ncinfo( treeVars(ctrtree+ 1,:), 'tr_qt'     , 'Moisture source sink' , '1/s'     , 'tttt' )
      call ncinfo( treeVars(ctrtree+ 2,:), 'tr_qtR'    , 'Moisture source sink' , '1/s'     , 'tttt' )
      call ncinfo( treeVars(ctrtree+ 3,:), 'tr_qtA'    , 'Moisture source sink' , '1/s'     , 'tttt' )
      call ncinfo( treeVars(ctrtree+ 4,:), 'tr_omega'  , 'Decoupling factor'    , '-'       , 'tttt' )
      ctrtree = ctrtree+4
    end subroutine stats_init_tree_moist

    subroutine stats_init_tree_scalar
      implicit none
      integer :: n
      character(2) :: sid
      allocate(svtreename(nsv))
      allocate(tr_svt(ib:ie,jb:je,kb:ke,nsv)); tr_svt = 0;
      do n = 1, nsv
        write (sid, '(I0)') n
        svtreename(n) = 'tr_sv'//trim(sid)         ! tr_sv1       at n = 1
        call ncinfo( treeVars(ctrtree+ n,:), trim(svtreename(n)) , 'Scalar source sink '//trim(sid) , 'kg/m^3s' , 'tttt' )
      end do
      ctrtree = ctrtree+nsv
    end subroutine stats_init_tree_scalar

    subroutine stats_createnc_tree
      implicit none
      filenametree = 'stats_tree_out.xxx.xxx.xxx.nc'
      filenametree(16:18) = cmyidx
      filenametree(20:22) = cmyidy
      filenametree(24:26) = cexpnr

      nrectree = 0
      call open_nc(filenametree, ncidtree, nrectree, n1=xdim, n2=ydim, n3=zdim)
      if (nrectree==0) then
        call define_nc(ncidtree, 1, timeVar)
        call writestat_dims_nc(ncidtree)
      end if
      call define_nc(ncidtree, treeVarsCount, treeVars)
    end subroutine stats_createnc_tree


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
            wc(i,j,k)  = 0.5*dzhi(k)*(wm(i,j,k)*dzf(k-1) + wm(i,j,k-1)*dzf(k))  !! Needs careful checking DMajumdar
            ! wc(i,j,k)  = 0.5*( wm(i,j,k+1) + wm(i,j,k) ) 
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
      !> Means
      call spatial_avg(uxyt,ut(ib:ie,jb:je,kb:ke+kh),kb,ke,IIu(ib:ie,jb:je,kb:ke+kh),IIus(kb:ke+kh),.false.)
      call spatial_avg(vxyt,vt(ib:ie,jb:je,kb:ke+kh),kb,ke,IIv(ib:ie,jb:je,kb:ke+kh),IIvs(kb:ke+kh),.false.)
      call spatial_avg(wxyt,wt(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
      call spatial_avg(pxyt,pt(ib:ie,jb:je,kb:ke+kh),kb,ke,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)

      !> Turbulent fluxes
      call spatial_avg(upwpxytik,uwtik(ib:ie,jb:je,kb:ke+kh)-utik(ib:ie,jb:je,kb:ke+kh)*wtik(ib:ie,jb:je,kb:ke+kh),kb,ke,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.false.)
      call spatial_avg(vpwpxytjk,vwtjk(ib:ie,jb:je,kb:ke+kh)-vtjk(ib:ie,jb:je,kb:ke+kh)*wtjk(ib:ie,jb:je,kb:ke+kh),kb,ke,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.false.)
      call spatial_avg(upvpxytij,uvtij(ib:ie,jb:je,kb:ke+kh)-utij(ib:ie,jb:je,kb:ke+kh)*vtij(ib:ie,jb:je,kb:ke+kh),kb,ke,IIuv(ib:ie,jb:je,kb:ke+kh),IIuvs(kb:ke+kh),.false.)

      !> Dispersive fluxes
      call spatial_avg(uwxytik,utik(ib:ie,jb:je,kb:ke+kh)*wtik(ib:ie,jb:je,kb:ke+kh),kb,ke,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.false.)
      call spatial_avg(vwxytjk,vtjk(ib:ie,jb:je,kb:ke+kh)*wtjk(ib:ie,jb:je,kb:ke+kh),kb,ke,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.false.)
      call spatial_avg(uvxytij,utij(ib:ie,jb:je,kb:ke+kh)*vtij(ib:ie,jb:je,kb:ke+kh),kb,ke,IIuv(ib:ie,jb:je,kb:ke+kh),IIuvs(kb:ke+kh),.false.)
      call spatial_avg(uuxyti,ut(ib:ie,jb:je,kb:ke+kh)*ut(ib:ie,jb:je,kb:ke+kh),kb,ke,IIu(ib:ie,jb:je,kb:ke+kh),IIus(kb:ke+kh),.false.)
      call spatial_avg(vvxytj,vt(ib:ie,jb:je,kb:ke+kh)*vt(ib:ie,jb:je,kb:ke+kh),kb,ke,IIv(ib:ie,jb:je,kb:ke+kh),IIvs(kb:ke+kh),.false.)
      call spatial_avg(wwxytk,wt(ib:ie,jb:je,kb:ke+kh)*wt(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)

      !> Variances and TKE
      call spatial_avg(upupxytc,uutc(ib:ie,jb:je,kb:ke+kh)-utc(ib:ie,jb:je,kb:ke+kh)*utc(ib:ie,jb:je,kb:ke+kh),kb,ke,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
      call spatial_avg(vpvpxytc,vvtc(ib:ie,jb:je,kb:ke+kh)-vtc(ib:ie,jb:je,kb:ke+kh)*vtc(ib:ie,jb:je,kb:ke+kh),kb,ke,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
      call spatial_avg(wpwpxytc,wwtc(ib:ie,jb:je,kb:ke+kh)-wtc(ib:ie,jb:je,kb:ke+kh)*wtc(ib:ie,jb:je,kb:ke+kh),kb,ke,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
      call spatial_avg(tkexytc,0.5*((wwtc(ib:ie,jb:je,kb:ke+kh)-wtc(ib:ie,jb:je,kb:ke+kh)*wtc(ib:ie,jb:je,kb:ke+kh))+(vvtc(ib:ie,jb:je,kb:ke+kh)-vtc(ib:ie,jb:je,kb:ke+kh)*vtc(ib:ie,jb:je,kb:ke+kh))+(uutc(ib:ie,jb:je,kb:ke+kh)-utc(ib:ie,jb:je,kb:ke+kh)*utc(ib:ie,jb:je,kb:ke+kh))),kb,ke,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)

      !> SGS fluxes
      call spatial_avg(usgsxyt,usgst(ib:ie,jb:je,kb:ke+kh),kb,ke,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.false.)
      call spatial_avg(vsgsxyt,vsgst(ib:ie,jb:je,kb:ke+kh),kb,ke,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.false.)
      call spatial_avg(wsgsxyt,wsgst(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
    end subroutine stats_compute_xytavg_vel

    subroutine stats_compute_xytavg_temp
      implicit none
      call spatial_avg(thlxyt,thlt(ib:ie,jb:je,kb:ke+kh),kb,ke,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
      call spatial_avg(wpthlpxytk,wthltk(ib:ie,jb:je,kb:ke+kh)-wt(ib:ie,jb:je,kb:ke+kh)*thltk(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
      call spatial_avg(wthlxytk,wt(ib:ie,jb:je,kb:ke+kh)*thltk(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
      call spatial_avg(thlpthlpxyt,thlthlt(ib:ie,jb:je,kb:ke+kh)-thlt(ib:ie,jb:je,kb:ke+kh)*thlt(ib:ie,jb:je,kb:ke+kh),kb,ke,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
      call spatial_avg(thlsgsxyt,thlsgst(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
    end subroutine stats_compute_xytavg_temp

    subroutine stats_compute_xytavg_moist
      implicit none
      call spatial_avg(qtxyt,qtt(ib:ie,jb:je,kb:ke+kh),kb,ke,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
      call spatial_avg(wpqtpxytk,wqttk(ib:ie,jb:je,kb:ke+kh)-wt(ib:ie,jb:je,kb:ke+kh)*qttk(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
      call spatial_avg(wqtxytk,wt(ib:ie,jb:je,kb:ke+kh)*qttk(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
      call spatial_avg(qtpqtpxyt,qtqtt(ib:ie,jb:je,kb:ke+kh)-qtt(ib:ie,jb:je,kb:ke+kh)*qtt(ib:ie,jb:je,kb:ke+kh),kb,ke,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
      call spatial_avg(qtsgsxyt,qtsgst(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
    end subroutine stats_compute_xytavg_moist


    !! ## %% y and x averaging computations routines
    subroutine stats_compute_xyavg_vel
      implicit none
      !> Mean
      call spatial_avg(uxy,um(ib:ie,jb:je,kb:ke+kh),kb,ke,IIu(ib:ie,jb:je,kb:ke+kh),IIus(kb:ke+kh),.false.)
      call spatial_avg(vxy,vm(ib:ie,jb:je,kb:ke+kh),kb,ke,IIv(ib:ie,jb:je,kb:ke+kh),IIvs(kb:ke+kh),.false.)
      call spatial_avg(wxy,wm(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
      call spatial_avg(pxy,pres0(ib:ie,jb:je,kb:ke+kh),kb,ke,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)

      !> Advective fluxes and some necesseary mean
      call spatial_avg(uwxyik,uik(ib:ie,jb:je,kb:ke+kh)*wik(ib:ie,jb:je,kb:ke+kh),kb,ke,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.true.)
      call spatial_avg(uxyik,uik(ib:ie,jb:je,kb:ke+kh),kb,ke,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.true.)
      call spatial_avg(wxyik,wik(ib:ie,jb:je,kb:ke+kh),kb,ke,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.true.)
      call spatial_avg(vwxyjk,vjk(ib:ie,jb:je,kb:ke+kh)*wjk(ib:ie,jb:je,kb:ke+kh),kb,ke,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.true.)
      call spatial_avg(vxyjk,vjk(ib:ie,jb:je,kb:ke+kh),kb,ke,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.true.)
      call spatial_avg(wxyjk,wjk(ib:ie,jb:je,kb:ke+kh),kb,ke,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.true.)
      call spatial_avg(uvxyij,uij(ib:ie,jb:je,kb:ke+kh)*vij(ib:ie,jb:je,kb:ke+kh),kb,ke,IIuv(ib:ie,jb:je,kb:ke+kh),IIuvs(kb:ke+kh),.true.)
      call spatial_avg(uxyij,uij(ib:ie,jb:je,kb:ke+kh),kb,ke,IIuv(ib:ie,jb:je,kb:ke+kh),IIuvs(kb:ke+kh),.true.)
      call spatial_avg(vxyij,vij(ib:ie,jb:je,kb:ke+kh),kb,ke,IIuv(ib:ie,jb:je,kb:ke+kh),IIuvs(kb:ke+kh),.true.)

      call spatial_avg(uuxyi,um(ib:ie,jb:je,kb:ke+kh)*um(ib:ie,jb:je,kb:ke+kh),kb,ke,IIu(ib:ie,jb:je,kb:ke+kh),IIus(kb:ke+kh),.true.)
      call spatial_avg(vvxyj,vm(ib:ie,jb:je,kb:ke+kh)*vm(ib:ie,jb:je,kb:ke+kh),kb,ke,IIv(ib:ie,jb:je,kb:ke+kh),IIvs(kb:ke+kh),.true.)
      call spatial_avg(wwxyk,wm(ib:ie,jb:je,kb:ke+kh)*wm(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.true.)

      !> Turbulent fluxes
      upwpxyik = uwxyik - uxyik*wxyik
      vpwpxyjk = vwxyjk - vxyjk*wxyjk
      upvpxyij = uvxyij - uxyij*vxyij

      !> SGS fluxes
      call spatial_avg(usgsxy,usgs(ib:ie,jb:je,kb:ke+kh),kb,ke,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.false.)
      call spatial_avg(vsgsxy,vsgs(ib:ie,jb:je,kb:ke+kh),kb,ke,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.false.)
      call spatial_avg(wsgsxy,wsgs(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
    end subroutine stats_compute_xyavg_vel

    subroutine stats_compute_xyavg_temp
      implicit none
      call spatial_avg(thlxy,thlm(ib:ie,jb:je,kb:ke+kh),kb,ke,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
      call spatial_avg(wthlxyk,wm(ib:ie,jb:je,kb:ke+kh)*thlk(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.true.)
      call spatial_avg(thlxyk,thlk(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.true.)
      wpthlpxyk = wthlxyk - wxy*thlxyk
      call spatial_avg(thlsgsxy,thlsgs(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
    end subroutine stats_compute_xyavg_temp

    subroutine stats_compute_xyavg_moist
      implicit none
      call spatial_avg(qtxy,qtm(ib:ie,jb:je,kb:ke+kh),kb,ke,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)
      call spatial_avg(wqtxyk,wm(ib:ie,jb:je,kb:ke+kh)*qtk(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.true.)
      call spatial_avg(qtxyk,qtk(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.true.)
      wpqtpxyk = wqtxyk - wxy*qtxyk
      call spatial_avg(qtsgsxy,qtsgs(ib:ie,jb:je,kb:ke+kh),kb,ke,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.false.)
    end subroutine stats_compute_xyavg_moist


    !! ## %% Time and y averaging computations routines
    subroutine stats_compute_ytavg_vel
      implicit none
      !> Mean
      call spatial_avg(uyt,ut(ib:ie,jb:je,kb:ke),IIu(ib:ie,jb:je,kb:ke),IIut)
      call spatial_avg(vyt,vt(ib:ie,jb:je,kb:ke),IIv(ib:ie,jb:je,kb:ke),IIvt)
      call spatial_avg(wyt,wt(ib:ie,jb:je,kb:ke),IIw(ib:ie,jb:je,kb:ke),IIwt)
      call spatial_avg(pyt,pt(ib:ie,jb:je,kb:ke),IIc(ib:ie,jb:je,kb:ke),IIct)

      !> Turbulent fluxes
      call spatial_avg(upwpytik,uwtik(ib:ie,jb:je,kb:ke)-utik(ib:ie,jb:je,kb:ke)*wtik(ib:ie,jb:je,kb:ke),IIuw(ib:ie,jb:je,kb:ke),IIuwt)
      
      !> Dispersive fluxes
      call spatial_avg(uwytik,utik(ib:ie,jb:je,kb:ke)*wtik(ib:ie,jb:je,kb:ke),IIuw(ib:ie,jb:je,kb:ke),IIuwt)

      !> Variances
      call spatial_avg(upupytc,uutc(ib:ie,jb:je,kb:ke)-utc(ib:ie,jb:je,kb:ke)*utc(ib:ie,jb:je,kb:ke),IIc(ib:ie,jb:je,kb:ke),IIct)
      call spatial_avg(vpvpytc,vvtc(ib:ie,jb:je,kb:ke)-vtc(ib:ie,jb:je,kb:ke)*vtc(ib:ie,jb:je,kb:ke),IIc(ib:ie,jb:je,kb:ke),IIct)
      call spatial_avg(wpwpytc,wwtc(ib:ie,jb:je,kb:ke)-wtc(ib:ie,jb:je,kb:ke)*wtc(ib:ie,jb:je,kb:ke),IIc(ib:ie,jb:je,kb:ke),IIct)

      !> SGS fluxes
      call spatial_avg(usgsyt,usgst(ib:ie,jb:je,kb:ke),IIuw(ib:ie,jb:je,kb:ke),IIuwt)
      call spatial_avg(wsgsyt,wsgst(ib:ie,jb:je,kb:ke),IIw(ib:ie,jb:je,kb:ke),IIwt)
    end subroutine stats_compute_ytavg_vel

    subroutine stats_compute_ytavg_temp
      implicit none
      call spatial_avg(thlyt,thlt(ib:ie,jb:je,kb:ke),IIc(ib:ie,jb:je,kb:ke),IIct)
      call spatial_avg(wpthlpytk,wthltk(ib:ie,jb:je,kb:ke)-wt(ib:ie,jb:je,kb:ke)*thltk(ib:ie,jb:je,kb:ke),IIw(ib:ie,jb:je,kb:ke),IIwt)
      call spatial_avg(wthlytk,wt(ib:ie,jb:je,kb:ke)*thltk(ib:ie,jb:je,kb:ke),IIw(ib:ie,jb:je,kb:ke),IIwt)
      call spatial_avg(thlpthlpyt,thlthlt(ib:ie,jb:je,kb:ke)-thlt(ib:ie,jb:je,kb:ke)*thlt(ib:ie,jb:je,kb:ke),IIc(ib:ie,jb:je,kb:ke),IIct)
      call spatial_avg(thlsgsyt,thlsgst(ib:ie,jb:je,kb:ke),IIw(ib:ie,jb:je,kb:ke),IIwt)
    end subroutine stats_compute_ytavg_temp

    subroutine stats_compute_ytavg_moist
      implicit none
      call spatial_avg(qtyt,qtt(ib:ie,jb:je,kb:ke),IIc(ib:ie,jb:je,kb:ke),IIct)
      call spatial_avg(wpqtpytk,wqttk(ib:ie,jb:je,kb:ke)-wt(ib:ie,jb:je,kb:ke)*qttk(ib:ie,jb:je,kb:ke),IIw(ib:ie,jb:je,kb:ke),IIwt)
      call spatial_avg(wqtytk,wt(ib:ie,jb:je,kb:ke)*qttk(ib:ie,jb:je,kb:ke),IIw(ib:ie,jb:je,kb:ke),IIwt)
      call spatial_avg(qtpqtpyt,qtqtt(ib:ie,jb:je,kb:ke)-qtt(ib:ie,jb:je,kb:ke)*qtt(ib:ie,jb:je,kb:ke),IIc(ib:ie,jb:je,kb:ke),IIct)
      call spatial_avg(qtsgsyt,qtsgst(ib:ie,jb:je,kb:ke),IIw(ib:ie,jb:je,kb:ke),IIwt)
    end subroutine stats_compute_ytavg_moist

    subroutine stats_compute_ytavg_scalar
      implicit none
      integer :: n
      do n = 1, nsv
        call spatial_avg(svyt(:,:,n),svt(ib:ie,jb:je,kb:ke,n),IIc(ib:ie,jb:je,kb:ke),IIct)
        call spatial_avg(wpsvpytk(:,:,n),wsvtk(ib:ie,jb:je,kb:ke,n)-wt(ib:ie,jb:je,kb:ke)*svtk(ib:ie,jb:je,kb:ke,n),IIw(ib:ie,jb:je,kb:ke),IIwt)
        call spatial_avg(wsvytk(:,:,n),wt(ib:ie,jb:je,kb:ke)*svtk(ib:ie,jb:je,kb:ke,n),IIw(ib:ie,jb:je,kb:ke),IIwt)
        call spatial_avg(svpsvpyt(:,:,n),svsvt(ib:ie,jb:je,kb:ke,n)-svt(ib:ie,jb:je,kb:ke,n)*svt(ib:ie,jb:je,kb:ke,n),IIc(ib:ie,jb:je,kb:ke),IIct)
        call spatial_avg(svsgsyt(:,:,n),svsgst(ib:ie,jb:je,kb:ke,n),IIw(ib:ie,jb:je,kb:ke),IIwt)
      end do
    end subroutine stats_compute_ytavg_scalar


    !! ## %% y averaging computations routines
    subroutine stats_compute_yavg_vel
      implicit none
      !> Mean
      call spatial_avg(uy,um(ib:ie,jb:je,kb:ke),IIu(ib:ie,jb:je,kb:ke),IIut)
      call spatial_avg(vy,vm(ib:ie,jb:je,kb:ke),IIv(ib:ie,jb:je,kb:ke),IIvt)
      call spatial_avg(wy,wm(ib:ie,jb:je,kb:ke),IIw(ib:ie,jb:je,kb:ke),IIwt)
      call spatial_avg(py,pres0(ib:ie,jb:je,kb:ke),IIc(ib:ie,jb:je,kb:ke),IIct)

      !> Advective fluxes and some necesseary mean
      call spatial_avg(uwyik,uik(ib:ie,jb:je,kb:ke)*wik(ib:ie,jb:je,kb:ke),IIuw(ib:ie,jb:je,kb:ke),IIuwt)
      call spatial_avg(uyik,uik(ib:ie,jb:je,kb:ke),IIuw(ib:ie,jb:je,kb:ke),IIuwt)
      call spatial_avg(wyik,wik(ib:ie,jb:je,kb:ke),IIuw(ib:ie,jb:je,kb:ke),IIuwt)

      !> Turbulent fluxes
      upwpyik = uwyik - uyik*wyik
      where (IIuwt==0)
        upwpyik    = -999.0
      endwhere

      !> SGS fluxes
      call spatial_avg(usgsy,usgs(ib:ie,jb:je,kb:ke),IIuw(ib:ie,jb:je,kb:ke),IIuwt)
      call spatial_avg(wsgsy,wsgs(ib:ie,jb:je,kb:ke),IIw(ib:ie,jb:je,kb:ke),IIwt)
    end subroutine stats_compute_yavg_vel

    subroutine stats_compute_yavg_temp
      implicit none
      call spatial_avg(thly,thlm(ib:ie,jb:je,kb:ke),IIc(ib:ie,jb:je,kb:ke),IIct)
      call spatial_avg(wthlyk,wm(ib:ie,jb:je,kb:ke)*thlk(ib:ie,jb:je,kb:ke),IIw(ib:ie,jb:je,kb:ke),IIwt)
      call spatial_avg(thlyk,thlk(ib:ie,jb:je,kb:ke),IIw(ib:ie,jb:je,kb:ke),IIwt)
      
      wpthlpyk = wthlyk - wy*thlyk
      where (IIwt==0)
        wpthlpyk  = -999.0
      endwhere
      
      call spatial_avg(thlsgsy,thlsgs(ib:ie,jb:je,kb:ke),IIw(ib:ie,jb:je,kb:ke),IIwt)
    end subroutine stats_compute_yavg_temp

    subroutine stats_compute_yavg_moist
      implicit none
      call spatial_avg(qty,qtm(ib:ie,jb:je,kb:ke),IIc(ib:ie,jb:je,kb:ke),IIct)
      call spatial_avg(wqtyk,wm(ib:ie,jb:je,kb:ke)*qtk(ib:ie,jb:je,kb:ke),IIw(ib:ie,jb:je,kb:ke),IIwt)
      call spatial_avg(qtyk,thlk(ib:ie,jb:je,kb:ke),IIw(ib:ie,jb:je,kb:ke),IIwt)
      
      wpqtpyk = wqtyk - wy*qtyk
      where (IIwt==0)
        wpqtpyk  = -999.0
      endwhere

      call spatial_avg(qtsgsy,qtsgs(ib:ie,jb:je,kb:ke),IIw(ib:ie,jb:je,kb:ke),IIwt)
    end subroutine stats_compute_yavg_moist

    subroutine stats_compute_yavg_scalar
      implicit none
      integer :: n
      do n = 1, nsv
        call spatial_avg(svy(:,:,n),svm(ib:ie,jb:je,kb:ke,n),IIc(ib:ie,jb:je,kb:ke),IIct)

        call spatial_avg(wsvyk(:,:,n),wm(ib:ie,jb:je,kb:ke)*svk(ib:ie,jb:je,kb:ke,n),IIw(ib:ie,jb:je,kb:ke),IIwt)
        call spatial_avg(svyk(:,:,n),svk(ib:ie,jb:je,kb:ke,n),IIw(ib:ie,jb:je,kb:ke),IIwt)

        wpsvpyk(:,:,n) = wsvyk(:,:,n) - wy*svyk(:,:,n)
        where (IIwt==0)
          wpsvpyk(:,:,n)  = -999.0
        endwhere

        call spatial_avg(svsgsy(:,:,n),svsgs(ib:ie,jb:je,kb:ke,n),IIw(ib:ie,jb:je,kb:ke),IIwt)
      end do
    end subroutine stats_compute_yavg_scalar


    !! ## %% Time averaging tree data computations routines
    subroutine stats_compute_tree_vel
      implicit none 
      call stats_compute_tavg(tr_ut,     tr_u(ib:ie,jb:je,kb:ke))
      call stats_compute_tavg(tr_vt,     tr_v(ib:ie,jb:je,kb:ke))
      call stats_compute_tavg(tr_wt,     tr_w(ib:ie,jb:je,kb:ke))
    end subroutine stats_compute_tree_vel

    subroutine stats_compute_tree_temp
      implicit none 
      call stats_compute_tavg(tr_thlt,   tr_thl(ib:ie,jb:je,kb:ke))
    end subroutine stats_compute_tree_temp

    subroutine stats_compute_tree_moist
      implicit none 
      call stats_compute_tavg(tr_qtt,    tr_qt(ib:ie,jb:je,kb:ke))
      call stats_compute_tavg(tr_qtRt,   tr_qtR(ib:ie,jb:je,kb:ke))
      call stats_compute_tavg(tr_qtAt,   tr_qtA(ib:ie,jb:je,kb:ke))
      call stats_compute_tavg(tr_omegat, tr_omega(ib:ie,jb:je,kb:ke))
    end subroutine stats_compute_tree_moist

    subroutine stats_compute_tree_scalar
      implicit none
      integer :: n
      do n = 1, nsv
        call stats_compute_tavg(tr_svt(:,:,:,n), tr_sv(ib:ie,jb:je,kb:ke,n))
      end do
    end subroutine stats_compute_tree_scalar


    !! ## %% Time averaged statistics writing routines 
    subroutine stats_write_tavg_vel
      implicit none
      call writestat_nc(ncidt, 'u', ut(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'v', vt(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'w', wt(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'p', pt(:,:,kb:ke), nrect, xdim, ydim, zdim)

      call writestat_nc(ncidt, 'upwp', uwtik(:,:,kb:ke) - utik(:,:,kb:ke)*wtik(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'vpwp', vwtjk(:,:,kb:ke) - vtjk(:,:,kb:ke)*wtjk(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'upvp', uvtij(:,:,kb:ke) - utij(:,:,kb:ke)*vtij(:,:,kb:ke), nrect, xdim, ydim, zdim)
      
      call writestat_nc(ncidt, 'upup', uutc(:,:,kb:ke)-utc(:,:,kb:ke)*utc(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'vpvp', vvtc(:,:,kb:ke)-vtc(:,:,kb:ke)*vtc(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'wpwp', wwtc(:,:,kb:ke)-wtc(:,:,kb:ke)*wtc(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'tke' , 0.5*( (uutc(:,:,kb:ke)-utc(:,:,kb:ke)*utc(:,:,kb:ke)) + (vvtc(:,:,kb:ke)-vtc(:,:,kb:ke)*vtc(:,:,kb:ke)) + (wwtc(:,:,kb:ke)-wtc(:,:,kb:ke)*wtc(:,:,kb:ke)) ) , nrect, xdim, ydim, zdim)

      call writestat_nc(ncidt, 'usgs', usgst(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'vsgs', vsgst(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'wsgs', wsgst(:,:,kb:ke), nrect, xdim, ydim, zdim)
    end subroutine stats_write_tavg_vel

    subroutine stats_write_tavg_temp
      implicit none
      call writestat_nc(ncidt, 'thl'     , thlt(:,:,kb:ke)                                     , nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'wpthlp'  , wthltk(:,:,kb:ke) - wt(:,:,kb:ke)*thltk(:,:,kb:ke)  , nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'thlpthlp', thlthlt(:,:,kb:ke) - thlt(:,:,kb:ke)*thlt(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'thlsgs'  , thlsgst(:,:,kb:ke)                                  , nrect, xdim, ydim, zdim)
    end subroutine stats_write_tavg_temp

    subroutine stats_write_tavg_moist
      implicit none
      call writestat_nc(ncidt, 'qt'    , qtt(:,:,kb:ke)                                  , nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'wpqtp' , wqttk(:,:,kb:ke) - wt(:,:,kb:ke)*qttk(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'qtpqtp', qtqtt(:,:,kb:ke) - qtt(:,:,kb:ke)*qtt(:,:,kb:ke), nrect, xdim, ydim, zdim)
      call writestat_nc(ncidt, 'qtsgs' , qtsgst(:,:,kb:ke)                               , nrect, xdim, ydim, zdim)
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
      call writestat_nc(ncidt, 'PSS', PSSt(:,:,kb:ke), nrect, xdim, ydim, zdim)
    end subroutine stats_write_tavg_PSS


    !! ## %% Time, y and x averaged statistics writing routines 
    subroutine stats_write_xytavg_vel
      implicit none
      call writestat_nc(ncidxyt, 'u'       , uxyt(kb:ke)       , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'v'       , vxyt(kb:ke)       , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'w'       , wxyt(kb:ke)       , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'p'       , pxyt(kb:ke)       , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'upwp'    , upwpxytik(kb:ke)  , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'vpwp'    , vpwpxytjk(kb:ke)  , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'upvp'    , upvpxytij(kb:ke)  , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'uw'      , uwxytik(kb:ke)    , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'vw'      , vwxytjk(kb:ke)    , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'uv'      , uvxytij(kb:ke)    , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'uu'      , uuxyti(kb:ke)     , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'vv'      , vvxytj(kb:ke)     , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'ww'      , wwxytk(kb:ke)     , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'upup'    , upupxytc(kb:ke)   , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'vpvp'    , vpvpxytc(kb:ke)   , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'wpwp'    , wpwpxytc(kb:ke)   , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'tke'     , tkexytc(kb:ke)    , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'usgs'    , usgsxyt(kb:ke)    , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'vsgs'    , vsgsxyt(kb:ke)    , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'wsgs'    , wsgsxyt(kb:ke)    , nrecxyt, zdim)
    end subroutine stats_write_xytavg_vel
    
    subroutine stats_write_xytavg_temp
      implicit none
      call writestat_nc(ncidxyt, 'thl'     , thlxyt(kb:ke)     , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'wpthlp'  , wpthlpxytk(kb:ke) , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'wthl'    , wthlxytk(kb:ke)   , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'thlpthlp', thlpthlpxyt(kb:ke), nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'thlsgs'  , thlsgsxyt(kb:ke)  , nrecxyt, zdim)
    end subroutine stats_write_xytavg_temp
    
    subroutine stats_write_xytavg_moist
      implicit none
      call writestat_nc(ncidxyt, 'qt'      , qtxyt(kb:ke)      , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'wpqtp'   , wpqtpxytk(kb:ke)  , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'wqt'     , wqtxytk(kb:ke)    , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'qtpqtp'  , qtpqtpxyt(kb:ke)  , nrecxyt, zdim)
      call writestat_nc(ncidxyt, 'qtsgs'   , qtsgsxyt(kb:ke)   , nrecxyt, zdim)
    end subroutine stats_write_xytavg_moist


    !! ## %% y and x averaged statistics writing routines 
    subroutine stats_write_xyavg_vel
      implicit none
      call writestat_nc(ncidxy, 'u'       , uxy(kb:ke)       , nrecxy, zdim)
      call writestat_nc(ncidxy, 'v'       , vxy(kb:ke)       , nrecxy, zdim)
      call writestat_nc(ncidxy, 'w'       , wxy(kb:ke)       , nrecxy, zdim)
      call writestat_nc(ncidxy, 'p'       , pxy(kb:ke)       , nrecxy, zdim)
      call writestat_nc(ncidxy, 'upwp'    , upwpxyik(kb:ke)  , nrecxy, zdim)
      call writestat_nc(ncidxy, 'vpwp'    , vpwpxyjk(kb:ke)  , nrecxy, zdim)
      call writestat_nc(ncidxy, 'upvp'    , upvpxyij(kb:ke)  , nrecxy, zdim)
      call writestat_nc(ncidxy, 'uw'      , uwxyik(kb:ke)    , nrecxy, zdim)
      call writestat_nc(ncidxy, 'vw'      , vwxyjk(kb:ke)    , nrecxy, zdim)
      call writestat_nc(ncidxy, 'uv'      , uvxyij(kb:ke)    , nrecxy, zdim)
      call writestat_nc(ncidxy, 'uu'      , uuxyi(kb:ke)     , nrecxy, zdim)
      call writestat_nc(ncidxy, 'vv'      , vvxyj(kb:ke)     , nrecxy, zdim)
      call writestat_nc(ncidxy, 'ww'      , wwxyk(kb:ke)     , nrecxy, zdim)
      call writestat_nc(ncidxy, 'usgs'    , usgsxy(kb:ke)    , nrecxy, zdim)
      call writestat_nc(ncidxy, 'vsgs'    , vsgsxy(kb:ke)    , nrecxy, zdim)
      call writestat_nc(ncidxy, 'wsgs'    , wsgsxy(kb:ke)    , nrecxy, zdim)
    end subroutine stats_write_xyavg_vel

    subroutine stats_write_xyavg_temp
      implicit none
      call writestat_nc(ncidxy, 'thl'     , thlxy(kb:ke)     , nrecxy, zdim)
      call writestat_nc(ncidxy, 'wpthlp'  , wpthlpxyk(kb:ke) , nrecxy, zdim)
      call writestat_nc(ncidxy, 'wthl'    , wthlxyk(kb:ke)   , nrecxy, zdim)
      call writestat_nc(ncidxy, 'thlsgs'  , thlsgsxy(kb:ke)  , nrecxy, zdim)
    end subroutine stats_write_xyavg_temp

    subroutine stats_write_xyavg_moist
      implicit none
      call writestat_nc(ncidxy, 'qt'      , qtxy(kb:ke)      , nrecxy, zdim)
      call writestat_nc(ncidxy, 'wpqtp'   , wpqtpxyk(kb:ke)  , nrecxy, zdim)
      call writestat_nc(ncidxy, 'wqt'     , wqtxyk(kb:ke)    , nrecxy, zdim)
      call writestat_nc(ncidxy, 'qtsgs'   , qtsgsxy(kb:ke)   , nrecxy, zdim)
    end subroutine stats_write_xyavg_moist


    !! ## %% Time and y averaged statistics writing routines 
    subroutine stats_write_ytavg_vel
      implicit none
      call writestat_nc(ncidyt, 'u'    , uyt     , nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'v'    , vyt     , nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'w'    , wyt     , nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'p'    , pyt     , nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'upwp' , upwpytik, nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'uw'   , uwytik  , nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'upup' , upupytc , nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'vpvp' , vpvpytc , nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'wpwp' , wpwpytc , nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'usgs' , usgsyt  , nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'wsgs' , wsgsyt  , nrecyt, xdim, zdim)
    end subroutine stats_write_ytavg_vel

    subroutine stats_write_ytavg_temp
      implicit none
      call writestat_nc(ncidyt, 'thl'     , thlyt     , nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'wpthlp'  , wpthlpytk , nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'wthl'    , wthlytk   , nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'thlpthlp', thlpthlpyt, nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'thlsgs'  , thlsgsyt  , nrecyt, xdim, zdim)
    end subroutine stats_write_ytavg_temp

    subroutine stats_write_ytavg_moist
      implicit none
      call writestat_nc(ncidyt, 'qt'     , qtyt     , nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'wpqtp'  , wpqtpytk , nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'wqt'    , wqtytk   , nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'qtpqtp' , qtpqtpyt , nrecyt, xdim, zdim)
      call writestat_nc(ncidyt, 'qtsgs'  , qtsgsyt  , nrecyt, xdim, zdim)
    end subroutine stats_write_ytavg_moist

    subroutine stats_write_ytavg_scalar
      implicit none
      integer :: n
      do n = 1, nsv
        call writestat_nc(ncidyt, trim(svytname(n))    , svyt(:,:,n)     , nrecyt, xdim, zdim)
        call writestat_nc(ncidyt, trim(wpsvpytname(n)) , wpsvpytk(:,:,n) , nrecyt, xdim, zdim)
        call writestat_nc(ncidyt, trim(wsvytname(n))   , wsvytk(:,:,n)   , nrecyt, xdim, zdim)
        call writestat_nc(ncidyt, trim(svpsvpytname(n)), svpsvpyt(:,:,n) , nrecyt, xdim, zdim)
        call writestat_nc(ncidyt, trim(svsgsytname(n)) , svsgsyt(:,:,n)  , nrecyt, xdim, zdim)
      end do
    end subroutine stats_write_ytavg_scalar


    !! ## %% y averaged statistics writing routines 
    subroutine stats_write_yavg_vel
      implicit none
      call writestat_nc(ncidy, 'u'    , uy     , nrecy, xdim, zdim)
      call writestat_nc(ncidy, 'v'    , vy     , nrecy, xdim, zdim)
      call writestat_nc(ncidy, 'w'    , wy     , nrecy, xdim, zdim)
      call writestat_nc(ncidy, 'p'    , py     , nrecy, xdim, zdim)
      call writestat_nc(ncidy, 'upwp' , upwpyik, nrecy, xdim, zdim)
      call writestat_nc(ncidy, 'uw'   , uwyik  , nrecy, xdim, zdim)
      call writestat_nc(ncidy, 'usgs' , usgsy  , nrecy, xdim, zdim)
      call writestat_nc(ncidy, 'wsgs' , wsgsy  , nrecy, xdim, zdim)
    end subroutine stats_write_yavg_vel

    subroutine stats_write_yavg_temp
      implicit none
      call writestat_nc(ncidy, 'thl'     , thly     , nrecy, xdim, zdim)
      call writestat_nc(ncidy, 'wpthlp'  , wpthlpyk , nrecy, xdim, zdim)
      call writestat_nc(ncidy, 'wthl'    , wthlyk   , nrecy, xdim, zdim)
      call writestat_nc(ncidy, 'thlsgs'  , thlsgsy  , nrecy, xdim, zdim)
    end subroutine stats_write_yavg_temp

    subroutine stats_write_yavg_moist
      implicit none
      call writestat_nc(ncidy, 'qt'     , qty     , nrecy, xdim, zdim)
      call writestat_nc(ncidy, 'wpqtp'  , wpqtpyk , nrecy, xdim, zdim)
      call writestat_nc(ncidy, 'wqt'    , wqtyk   , nrecy, xdim, zdim)
      call writestat_nc(ncidy, 'qtsgs'  , qtsgsy  , nrecy, xdim, zdim)
    end subroutine stats_write_yavg_moist

    subroutine stats_write_yavg_scalar
      implicit none
      integer :: n
      do n = 1, nsv
        call writestat_nc(ncidy, trim(svyname(n))    , svy(:,:,n)     , nrecy, xdim, zdim)
        call writestat_nc(ncidy, trim(wpsvpyname(n)) , wpsvpyk(:,:,n) , nrecy, xdim, zdim)
        call writestat_nc(ncidy, trim(wsvyname(n))   , wsvyk(:,:,n)   , nrecy, xdim, zdim)
        call writestat_nc(ncidy, trim(svsgsyname(n)) , svsgsy(:,:,n)  , nrecy, xdim, zdim)
      end do
    end subroutine stats_write_yavg_scalar


    !! ## %% Time averaged tree data statistics writing routines 
    subroutine stats_write_tree_vel
      implicit none
      call writestat_nc(ncidtree, 'tr_u'    , tr_ut    , nrectree, xdim, ydim, zdim)
      call writestat_nc(ncidtree, 'tr_v'    , tr_vt    , nrectree, xdim, ydim, zdim)
      call writestat_nc(ncidtree, 'tr_w'    , tr_wt    , nrectree, xdim, ydim, zdim)
    end subroutine stats_write_tree_vel

    subroutine stats_write_tree_temp
      implicit none
      call writestat_nc(ncidtree, 'tr_thl'  , tr_thlt  , nrectree, xdim, ydim, zdim)
    end subroutine stats_write_tree_temp

    subroutine stats_write_tree_moist
      implicit none
      call writestat_nc(ncidtree, 'tr_qt'   , tr_qtt   , nrectree, xdim, ydim, zdim)
      call writestat_nc(ncidtree, 'tr_qtR'  , tr_qtRt  , nrectree, xdim, ydim, zdim)
      call writestat_nc(ncidtree, 'tr_qtA'  , tr_qtAt  , nrectree, xdim, ydim, zdim)
      call writestat_nc(ncidtree, 'tr_omega', tr_omegat, nrectree, xdim, ydim, zdim)
    end subroutine stats_write_tree_moist

    subroutine stats_write_tree_scalar
      implicit none
      integer :: n
      do n = 1, nsv
        call writestat_nc(ncidtree, trim(svtreename(n)), tr_svt(:,:,:,n), nrectree, xdim, ydim, zdim)
      end do
    end subroutine stats_write_tree_scalar

    subroutine stats_exit
      implicit none
      if(.not.(ltdump .or. lxytdump .or. lxydump .or. lytdump .or. lydump .or. ltreedump)) return

      if (ltdump .or. lxytdump .or. lxydump .or. lytdump .or. lydump) then
        deallocate(uik,wik,vjk,wjk,uij,vij,uc,vc,wc,usgs,vsgs,wsgs)
        if (ltempeq) deallocate(thlk,thlsgs)
        if (lmoist)  deallocate(qtk,qtsgs)
      end if
      if (ltdump .or. lytdump) then
        if (nsv>0)   deallocate(svk,svsgs)
      end if

      if (ltdump .or. lxytdump .or. lytdump) then
        deallocate(ut,vt,wt,pt)
        deallocate(utc,vtc,wtc,uutc,vvtc,wwtc)
        deallocate(utik,wtik,uwtik,vtjk,wtjk,vwtjk,utij,vtij,uvtij)
        deallocate(usgst,vsgst,wsgst)
        if (ltempeq) deallocate(thlt,thltk,wthltk,thlthlt,thlsgst)
        if (lmoist)  deallocate(qtt,qttk,wqttk,qtqtt,qtsgst)
      end if
      if (ltdump .or. lytdump) then
        if (nsv>0)   deallocate(svt,svtk,wsvtk,svsvt,svsgst)
      end if
      if (ltdump) then  
        if (nsv>0)   deallocate(svtname,wpsvptname,svpsvptname,svsgsname)
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

      if (lytdump) then
        deallocate(uyt,vyt,wyt,pyt,usgsyt,wsgsyt)
        deallocate(upwpytik,uwytik,upupytc,vpvpytc,wpwpytc)
        if (ltempeq) deallocate(thlyt,wpthlpytk,wthlytk,thlpthlpyt,thlsgsyt)
        if (lmoist)  deallocate(qtyt,wpqtpytk,wqtytk,qtpqtpyt,qtsgsyt)
        if (nsv>0)   deallocate(svytname,wpsvpytname,wsvytname,svpsvpytname,svsgsytname,svyt,wpsvpytk,wsvytk,svpsvpyt,svsgsyt)
      end if

      if (lytdump) then
        deallocate(uy,vy,wy,py,usgsy,wsgsy)
        deallocate(upwpyik,uwyik,uyik,wyik)
        if (ltempeq) deallocate(thly,wpthlpyk,wthlyk,thlyk,thlsgsy)
        if (lmoist)  deallocate(qty,wpqtpyk,wqtyk,qtyk,qtsgsy)
        if (nsv>0)   deallocate(svyname,wpsvpyname,wsvyname,svsgsyname,svy,wpsvpyk,wsvyk,svyk,svsgsy)
      end if

      if (ltreedump) then
        deallocate(tr_ut,tr_vt,tr_wt)
        if (ltempeq) deallocate(tr_thlt)
        if (lmoist)  deallocate(tr_qtt,tr_qtRt,tr_qtAt,tr_omegat)
        if (nsv>0)   deallocate(svtreename,tr_svt)
      end if
    end subroutine stats_exit
end module stats