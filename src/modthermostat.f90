!> Cloud, surface and radiation diagnostics: domain-mean time series and slab
!> profiles, the uDALES counterpart of DALES tmser/profiles. Kept separate from
!> modstatsdump so it can move into the statsout framework as one unit.
module modthermostat
  implicit none
  save
  private
  public :: initthermostat, thermostat, exitthermostat, lthermostat, tthermostat

  logical :: lthermostat = .false.
  real    :: tthermostat = 60.    !< output interval [s]
  real    :: tnext = 0.
  integer :: ncidt, ncidz, nrect = 0, nrecz = 0
  integer, parameter :: nvt = 11, nvz = 11
  character(80) :: fnamet = 'thermostat_t.xxx.nc', fnamez = 'thermostat_z.xxx.nc'
  character(80) :: tnct(1,4), tncz(1,4), nct(nvt,4), ncz(nvz,4)

contains

  subroutine initthermostat
    use modglobal, only : ifnamopt, fname_options, cexpnr, kb, ke
    use modmpi,    only : myid, comm3d, mpierr, my_real, mpi_logical
    use modstat_nc, only : ncinfo, open_nc, define_nc, writestat_dims_nc
    implicit none
    integer :: ierr
    namelist/NAMTHERMOSTAT/ lthermostat, tthermostat

    if (myid == 0) then
      open(ifnamopt, file=fname_options, status='old', iostat=ierr)
      read (ifnamopt, NAMTHERMOSTAT, iostat=ierr)
      if (ierr > 0) then
        write(0, *) 'ERROR: Problem in namoptions NAMTHERMOSTAT'
        write(0, *) 'iostat error: ', ierr
        stop 1
      end if
      close(ifnamopt)
    end if
    call MPI_BCAST(lthermostat, 1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(tthermostat, 1, MY_REAL,     0, comm3d, mpierr)
    if (.not. lthermostat) return

    tnext = tthermostat
    fnamet(14:16) = cexpnr
    fnamez(14:16) = cexpnr
    call ncinfo(tnct(1,:), 'time', 'Time', 's', 'time')
    tncz = tnct
    call ncinfo(nct( 1,:), 'lwp_bar',  'Liquid water path',             'kg/m^2', 'time')
    call ncinfo(nct( 2,:), 'cfrac',    'Cloud cover',                   '-',      'time')
    call ncinfo(nct( 3,:), 'zi',       'BL height (max dthl/dz)',       'm',      'time')
    call ncinfo(nct( 4,:), 'zi_qt8',   'BL height (qt = 8 g/kg)',       'm',      'time')
    call ncinfo(nct( 5,:), 'zi_smoke', 'BL height (smoke = 0.5)',       'm',      'time')
    call ncinfo(nct( 6,:), 'wmax',     'Maximum vertical velocity',     'm/s',    'time')
    call ncinfo(nct( 7,:), 'qlmax',    'Maximum liquid water',          'kg/kg',  'time')
    call ncinfo(nct( 8,:), 'ustar',    'Friction velocity',             'm/s',    'time')
    call ncinfo(nct( 9,:), 'obl',      'Obukhov length',                'm',      'time')
    call ncinfo(nct(10,:), 'thls',     'Surface thl',                   'K',      'time')
    call ncinfo(nct(11,:), 'qts',      'Surface qt',                    'kg/kg',  'time')
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
      call define_nc(ncidt, nvt, nct)
      call open_nc(fnamez, ncidz, nrecz, n3=ke-kb+1)
      if (nrecz == 0) then
        call define_nc(ncidz, 1, tncz)
        call writestat_dims_nc(ncidz)      ! fill zt/zm, as modstatsdump.f90:219-222
      end if
      call define_nc(ncidz, nvz, ncz)
    end if
  end subroutine initthermostat

  subroutine thermostat
    use modglobal,  only : ib, ie, jb, je, kb, ke, itot, jtot, timee, rk3step, zf, zh, dzf, dzh, &
                           nsv, lmoist
    use modfields,  only : thl0, qt0, ql0, w0, sv0, rhof, exnf, IIc, IIcs, IIw, IIws
    use modsurfdata, only : thls, qts
    use modmpi,     only : myid, comm3d, mpierr, my_real, mpi_sum, mpi_max, avexy_ibm
    use modstat_nc, only : writestat_nc
    implicit none
    real :: vt(nvt), vz(kb:ke, nvz)
    real, allocatable :: tmp3(:,:,:)   ! allocatable: an automatic array is ~22 MB/rank on case 303
    real :: loc(7), glob(7), colwp, gmax, dthl
    integer :: i, j, k, kmx
    real :: ncol

    if (.not. lthermostat) return
    if (rk3step /= 3 .or. timee < tnext) return
    tnext = tnext + tthermostat
    ncol = real(itot) * real(jtot)
    allocate(tmp3(ib:ie,jb:je,kb:ke+1))

    ! --- column diagnostics: lwp, cloud cover, three BL heights
    loc = 0.
    do j = jb, je
      do i = ib, ie
        colwp = 0.
        do k = kb, ke
          colwp = colwp + rhof(k) * ql0(i,j,k) * dzf(k) * IIc(i,j,k)
        end do
        loc(1) = loc(1) + colwp
        if (any(ql0(i,j,kb:ke) * IIc(i,j,kb:ke) > 0.)) loc(2) = loc(2) + 1.
        gmax = -huge(1.); kmx = kb + 1
        do k = kb + 1, ke
          dthl = (thl0(i,j,k) - thl0(i,j,k-1)) / dzh(k)
          if (dthl > gmax) then
            gmax = dthl; kmx = k
          end if
        end do
        loc(3) = loc(3) + zh(kmx)
        if (lmoist) loc(4) = loc(4) + crossing(qt0(i,j,kb:ke), 8.e-3)
        if (nsv > 0) loc(5) = loc(5) + crossing(sv0(i,j,kb:ke,1), 0.5)   ! Task 5: index -> isvsmoke
      end do
    end do
    loc(6) = maxval(w0(ib:ie,jb:je,kb:ke))
    loc(7) = maxval(ql0(ib:ie,jb:je,kb:ke))
    call MPI_ALLREDUCE(loc(1:5), glob(1:5), 5, MY_REAL, MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(loc(6:7), glob(6:7), 2, MY_REAL, MPI_MAX, comm3d, mpierr)
    vt(1:5) = glob(1:5) / ncol
    vt(6:7) = glob(6:7)
    vt(8)  = 0.          ! ustar: Task 3
    vt(9)  = 0.          ! obl:   Task 3
    vt(10) = thls
    vt(11) = qts

    ! --- slab profiles (fluid cells only)
    call avexy_ibm(vz(kb:ke,1), thl0(ib:ie,jb:je,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIc(ib:ie,jb:je,kb:ke), IIcs(kb:ke), .false.)
    call avexy_ibm(vz(kb:ke,2), qt0 (ib:ie,jb:je,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIc(ib:ie,jb:je,kb:ke), IIcs(kb:ke), .false.)
    call avexy_ibm(vz(kb:ke,3), ql0 (ib:ie,jb:je,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIc(ib:ie,jb:je,kb:ke), IIcs(kb:ke), .false.)
    tmp3(:,:,kb:ke) = merge(1., 0., ql0(ib:ie,jb:je,kb:ke) > 0.)
    call avexy_ibm(vz(kb:ke,4), tmp3(:,:,kb:ke), ib,ie,jb,je,kb,ke-1,1, IIc(ib:ie,jb:je,kb:ke), IIcs(kb:ke), .false.)
    vz(:,5) = rhof(kb:ke)
    vz(:,6) = exnf(kb:ke)
    vz(:,7) = 0.         ! thlprad: Task 5
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
    vz(:,10) = 0.        ! lwu: Task 5
    vz(:,11) = 0.        ! lwd: Task 5

    if (myid == 0) then
      call writestat_nc(ncidt, 1, tnct, (/timee/), nrect, .true.)
      call writestat_nc(ncidt, nvt, nct, vt, nrect, .false.)
      call writestat_nc(ncidz, 1, tncz, (/timee/), nrecz, .true.)
      call writestat_nc(ncidz, nvz, ncz, vz, nrecz, ke-kb+1)
    end if
    deallocate(tmp3)

  contains

    !> Height of the highest downward crossing of `thres` in a column; 0 if none.
    real function crossing(col, thres)
      real, intent(in) :: col(kb:ke), thres
      integer :: kk
      crossing = 0.
      do kk = ke - 1, kb, -1
        if (col(kk) > thres .and. col(kk+1) <= thres) then
          crossing = zf(kk) + (thres - col(kk)) * (zf(kk+1) - zf(kk)) / (col(kk+1) - col(kk))
          return
        end if
      end do
    end function crossing

  end subroutine thermostat

  subroutine exitthermostat
    use modmpi, only : myid
    use modstat_nc, only : exitstat_nc
    implicit none
    if (.not. lthermostat) return
    if (myid == 0) then
      call exitstat_nc(ncidt)
      call exitstat_nc(ncidz)
    end if
  end subroutine exitthermostat

end module modthermostat
