!> Cloud, boundary-layer and radiation diagnostics (bldump): domain-mean time
!> series and slab profiles, the uDALES counterpart of DALES tmser/profiles.
!> Called from modstatsdump (initstatsdump, exitstatsdump) and sampled from the
!> time loop after thermodynamics, every tstatsdump from tstatstart on. Written
!> to bldump_t.<expnr>.nc and bldump_z.<expnr>.nc when lbldump is set.
module modbldiag
  use mpi   ! explicit MPI interfaces: gfortran rejects MPI_BCAST calls of different types in one file otherwise
  implicit none
  save
  private
  public :: initbldiag, bldiag, exitbldiag

  real    :: tnext = 0.
  integer :: ncidt, ncidz, nrect = 0, nrecz = 0
  integer, parameter :: nvtmax = 13, nvz = 11, ivt_zi_qt = 4, ivt_zi_sv = 5
  integer :: nvt = 0               !< number of time series written
  integer :: ivt(nvtmax)           !< indices into vt/nct of the time series written
  character(80) :: fnamet = 'bldump_t.xxx.nc', fnamez = 'bldump_z.xxx.nc'
  character(80) :: tnct(1,4), tncz(1,4), nct(nvtmax,4), ncz(nvz,4)

contains

  subroutine initbldiag
    use modglobal, only : cexpnr, kb, ke, timee, tstatsdump, lbldump, zi_qt_thres, zi_sv_index
    use modmpi,    only : myid
    use modstat_nc, only : ncinfo, open_nc, define_nc, writestat_dims_nc
    implicit none
    integer :: n

    call check_bldiag_options
    if (.not. lbldump) return

    tnext = timee + tstatsdump       ! timee is set by readinitfiles (restart time) before this
    fnamet(10:12) = cexpnr
    fnamez(10:12) = cexpnr
    call ncinfo(tnct(1,:), 'time', 'Time', 's', 'time')
    tncz = tnct
    call ncinfo(nct( 1,:), 'lwp_bar',  'Liquid water path',             'kg/m^2', 'time')
    call ncinfo(nct( 2,:), 'cfrac',    'Cloud cover',                   '-',      'time')
    call ncinfo(nct( 3,:), 'zi',       'BL height (max dthl/dz)',       'm',      'time')
    call ncinfo(nct( 4,:), 'zi_qt',    'BL height (qt = zi_qt_thres)',  'm',      'time')
    call ncinfo(nct( 5,:), 'zi_sv',    'BL height (scalar zi_sv_index = zi_sv_thres)', 'm', 'time')
    call ncinfo(nct( 6,:), 'wmax',     'Maximum vertical velocity',     'm/s',    'time')
    call ncinfo(nct( 7,:), 'qlmax',    'Maximum liquid water',          'kg/kg',  'time')
    call ncinfo(nct( 8,:), 'ustar',    'Friction velocity',             'm/s',    'time')
    call ncinfo(nct( 9,:), 'obl',      'Obukhov length',                'm',      'time')
    call ncinfo(nct(10,:), 'thls',     'Surface thl',                   'K',      'time')
    call ncinfo(nct(11,:), 'qts',      'Surface qt',                    'kg/kg',  'time')
    call ncinfo(nct(12,:), 'u1',       'Slab-mean u at first level',    'm/s',    'time')
    call ncinfo(nct(13,:), 'v1',       'Slab-mean v at first level',    'm/s',    'time')
    nvt = 0
    do n = 1, nvtmax
      if (n == ivt_zi_qt .and. zi_qt_thres == 0.) cycle   ! zi_qt only when configured
      if (n == ivt_zi_sv .and. zi_sv_index == 0) cycle    ! zi_sv only when configured
      nvt = nvt + 1
      ivt(nvt) = n
    end do
    call ncinfo(ncz( 1,:), 'thl',     'Liquid water potential temperature', 'K',      'tt')
    call ncinfo(ncz( 2,:), 'qt',      'Total water',                        'kg/kg',  'tt')
    call ncinfo(ncz( 3,:), 'ql',      'Liquid water',                       'kg/kg',  'tt')
    call ncinfo(ncz( 4,:), 'cfracz',  'Cloud fraction',                     '-',      'tt')
    call ncinfo(ncz( 5,:), 'rhof',    'Density',                            'kg/m^3', 'tt')
    call ncinfo(ncz( 6,:), 'exnf',    'Exner function',                     '-',      'tt')
    call ncinfo(ncz( 7,:), 'thlprad', 'Radiative thl tendency',             'K/s',    'tt')
    call ncinfo(ncz( 8,:), 'wthl',    'Resolved w thl flux',                'K m/s',  'mt')
    call ncinfo(ncz( 9,:), 'wqt',     'Resolved w qt flux',                 'm/s',    'mt')
    call ncinfo(ncz(10,:), 'lwu',     'Longwave flux, F1 term',             'W/m^2',  'mt')
    call ncinfo(ncz(11,:), 'lwd',     'Longwave flux, F0 term',             'W/m^2',  'mt')
    if (myid == 0) then
      call open_nc(fnamet, ncidt, nrect)
      if (nrect == 0) call define_nc(ncidt, 1, tnct)
      call define_nc(ncidt, nvt, nct(ivt(1:nvt),:))
      call open_nc(fnamez, ncidz, nrecz, n3=ke-kb+1)
      if (nrecz == 0) then
        call define_nc(ncidz, 1, tncz)
        call writestat_dims_nc(ncidz)      ! fill zt/zm, as modstatsdump does
      end if
      call define_nc(ncidz, nvz, ncz)
    end if
  end subroutine initbldiag

  subroutine bldiag
    use modglobal,  only : ib, ie, jb, je, kb, ke, timee, rk3step, zf, zh, dzf, dzh, &
                           lbldump, tstatsdump, tstatstart, zi_qt_thres, zi_sv_index, zi_sv_thres
    use modfields,  only : thl0, qt0, ql0, w0, sv0, rhof, exnf, IIc, IIcs, IIw, IIws, u0av, v0av
    use modsurfdata, only : thls, qts
    use modraddata, only : thlprad, lwu, lwd, iradiation, irad_none
    use modsurface, only : ustar_mo, obl_mo
    use modmpi,     only : myid, comm3d, mpierr, my_real, mpi_sum, mpi_max, avexy_ibm
    use modstat_nc, only : writestat_nc
    implicit none
    real :: vt(nvtmax), vz(kb:ke, nvz)
    real, allocatable :: tmp3(:,:,:)   ! allocatable: an automatic array is ~22 MB/rank on case 303
    real :: loc(8), glob(8), colwp, gmax, dthl
    integer :: i, j, k, kmx
    real :: ncol

    if (.not. lbldump) return
    if (rk3step /= 3 .or. timee < tstatstart .or. timee < tnext) return
    tnext = tnext + tstatsdump
    allocate(tmp3(ib:ie,jb:je,kb:ke+1))

    ! --- column diagnostics: lwp, cloud cover, three BL heights
    loc = 0.
    do j = jb, je
      do i = ib, ie
        if (sum(IIc(i,j,kb:ke)) == 0) cycle      ! column entirely inside a building
        loc(8) = loc(8) + 1.
        colwp = 0.
        do k = kb, ke
          colwp = colwp + rhof(k) * ql0(i,j,k) * dzf(k) * IIc(i,j,k)
        end do
        loc(1) = loc(1) + colwp
        if (any(ql0(i,j,kb:ke) * IIc(i,j,kb:ke) > 0.)) loc(2) = loc(2) + 1.
        gmax = -huge(1.); kmx = kb + 1
        do k = kb + 1, ke
          if (IIc(i,j,k) == 0 .or. IIc(i,j,k-1) == 0) cycle
          dthl = (thl0(i,j,k) - thl0(i,j,k-1)) / dzh(k)
          if (dthl > gmax) then
            gmax = dthl; kmx = k
          end if
        end do
        loc(3) = loc(3) + zh(kmx)
        if (zi_qt_thres > 0.) loc(4) = loc(4) + crossing(qt0(i,j,kb:ke), IIc(i,j,kb:ke), zi_qt_thres)
        if (zi_sv_index > 0) loc(5) = loc(5) + crossing(sv0(i,j,kb:ke,zi_sv_index), IIc(i,j,kb:ke), zi_sv_thres)
      end do
    end do
    loc(6) = maxval(w0(ib:ie,jb:je,kb:ke))
    loc(7) = maxval(ql0(ib:ie,jb:je,kb:ke))
    call MPI_ALLREDUCE(loc(1:5), glob(1:5), 5, MY_REAL, MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(loc(6:7), glob(6:7), 2, MY_REAL, MPI_MAX, comm3d, mpierr)
    call MPI_ALLREDUCE(loc(8), glob(8), 1, MY_REAL, MPI_SUM, comm3d, mpierr)
    ncol = max(glob(8), 1.)                  ! columns with at least one fluid cell
    vt(1:5) = glob(1:5) / ncol
    vt(6:7) = glob(6:7)
    vt(8)  = ustar_mo
    vt(9)  = obl_mo
    vt(10) = thls
    vt(11) = qts
    vt(12) = u0av(kb)
    vt(13) = v0av(kb)

    ! --- slab profiles (fluid cells only)
    call avexy_ibm(vz(kb:ke,1), thl0(ib:ie,jb:je,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIc(ib:ie,jb:je,kb:ke), IIcs(kb:ke), .false.)
    call avexy_ibm(vz(kb:ke,2), qt0 (ib:ie,jb:je,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIc(ib:ie,jb:je,kb:ke), IIcs(kb:ke), .false.)
    call avexy_ibm(vz(kb:ke,3), ql0 (ib:ie,jb:je,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIc(ib:ie,jb:je,kb:ke), IIcs(kb:ke), .false.)
    tmp3(:,:,kb:ke) = merge(1., 0., ql0(ib:ie,jb:je,kb:ke) > 0.)
    call avexy_ibm(vz(kb:ke,4), tmp3(:,:,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIc(ib:ie,jb:je,kb:ke), IIcs(kb:ke), .false.)
    vz(:,5) = rhof(kb:ke)
    vz(:,6) = exnf(kb:ke)
    if (iradiation /= irad_none) then
      call avexy_ibm(vz(kb:ke,7), thlprad(ib:ie,jb:je,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIc(ib:ie,jb:je,kb:ke), IIcs(kb:ke), .false.)
    else
      vz(:,7) = 0.
    end if
    ! resolved fluxes at half levels k (face between k-1 and k); kb face is 0
    do k = kb + 1, ke
      tmp3(:,:,k) = w0(ib:ie,jb:je,k) * 0.5 * (thl0(ib:ie,jb:je,k) + thl0(ib:ie,jb:je,k-1))
    end do
    tmp3(:,:,kb) = 0.
    call avexy_ibm(vz(kb:ke,8), tmp3(:,:,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIw(ib:ie,jb:je,kb:ke), IIws(kb:ke), .false.)
    do k = kb + 1, ke
      tmp3(:,:,k) = w0(ib:ie,jb:je,k) * 0.5 * (qt0(ib:ie,jb:je,k) + qt0(ib:ie,jb:je,k-1))
    end do
    tmp3(:,:,kb) = 0.
    call avexy_ibm(vz(kb:ke,9), tmp3(:,:,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIw(ib:ie,jb:je,kb:ke), IIws(kb:ke), .false.)
    if (iradiation /= irad_none) then
      ! F1 / F0 terms at faces kb:ke (the top face ke+1 is not written); mean over fluid-face columns
      call avexy_ibm(vz(kb:ke,10), lwu(ib:ie,jb:je,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIw(ib:ie,jb:je,kb:ke), IIws(kb:ke), .false.)
      call avexy_ibm(vz(kb:ke,11), lwd(ib:ie,jb:je,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIw(ib:ie,jb:je,kb:ke), IIws(kb:ke), .false.)
    else
      vz(:,10) = 0.
      vz(:,11) = 0.
    end if

    if (myid == 0) then
      call writestat_nc(ncidt, 1, tnct, (/timee/), nrect, .true.)
      call writestat_nc(ncidt, nvt, nct(ivt(1:nvt),:), vt(ivt(1:nvt)), nrect, .false.)
      call writestat_nc(ncidz, 1, tncz, (/timee/), nrecz, .true.)
      call writestat_nc(ncidz, nvz, ncz, vz, nrecz, ke-kb+1)
    end if
    deallocate(tmp3)

  contains

    !> Height of the highest downward crossing of `thres` in a column; 0 if none.
    real function crossing(col, mask, thres)
      real, intent(in) :: col(kb:ke), thres
      integer, intent(in) :: mask(kb:ke)
      integer :: kk
      crossing = 0.
      do kk = ke - 1, kb, -1
        if (mask(kk) == 1 .and. mask(kk+1) == 1 .and. col(kk) > thres .and. col(kk+1) <= thres) then
          crossing = zf(kk) + (thres - col(kk)) * (zf(kk+1) - zf(kk)) / (col(kk+1) - col(kk))
          return
        end if
      end do
    end function crossing

  end subroutine bldiag

  !> Abort on settings that cannot be meaningful; warn (rank 0) on settings that are inactive.
  subroutine check_bldiag_options
    use modglobal, only : nsv, lmoist, lbldump, tstatsdump, zi_qt_thres, zi_sv_index, zi_sv_thres
    use modmpi,    only : myid, abort_all_ranks
    implicit none
    if (zi_qt_thres < 0. .or. (zi_qt_thres > 0. .and. .not. lmoist)) then
      if (myid == 0) write(0, '(a,es12.4)') 'ERROR: &OUTPUT zi_qt_thres must be >= 0, and > 0 requires lmoist; got ', &
                                 zi_qt_thres
      call abort_all_ranks(1)
    end if
    if (zi_sv_index < 0 .or. zi_sv_index > nsv .or. (zi_sv_index > 0 .and. zi_sv_thres <= 0.)) then
      if (myid == 0) write(0, '(2a,i0,es12.4)') 'ERROR: &OUTPUT zi_sv_index > 0 requires zi_sv_index <= nsv ', &
                                 'and zi_sv_thres > 0; got ', zi_sv_index, zi_sv_thres
      call abort_all_ranks(1)
    end if
    if (lbldump .and. tstatsdump <= 0.) then
      if (myid == 0) write(0, '(a,es12.4)') 'ERROR: &OUTPUT lbldump = .true. requires tstatsdump > 0; got ', tstatsdump
      call abort_all_ranks(1)
    end if
    if (myid /= 0) return
    if (zi_sv_index == 0 .and. zi_sv_thres /= 0.) &
      write(0, '(a)') 'WARNING: &OUTPUT zi_sv_thres ignored: zi_sv_index = 0'
    if (.not. lbldump .and. (zi_qt_thres /= 0. .or. zi_sv_index /= 0)) &
      write(0, '(a)') 'WARNING: &OUTPUT zi_qt_thres/zi_sv_index ignored: lbldump = .false.'
  end subroutine check_bldiag_options

  subroutine exitbldiag
    use modglobal, only : lbldump
    use modmpi, only : myid
    use modstat_nc, only : exitstat_nc
    implicit none
    if (.not. lbldump) return
    if (myid == 0) then
      call exitstat_nc(ncidt)
      call exitstat_nc(ncidz)
    end if
  end subroutine exitbldiag

end module modbldiag
