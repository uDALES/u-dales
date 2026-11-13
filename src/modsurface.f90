!> Monin-Obukhov surface layer for the flat bottom (lbottom), ported from
!> DALES dev 075ac4e src/surface/modsurface.f90 (isurf = 3/4, lmostlocal = .false.):
!> domain-mean Obukhov length and friction velocity, locally applied stress.
!>   BCbotm = 4 (BCbotm_mo):       u* from MO similarity            (DALES isurf = 4)
!>   BCbotm = 5 (BCbotm_mo_ustar): u* = ustin                       (DALES isurf = 3)
!> Heat and moisture fluxes stay prescribed by wtsurf/wqsurf (BCbotT/BCbotq = 1).
module modsurface
  implicit none
  save
  private
  public :: surface_mo, lmo_surface, ustar_mo, obl_mo, dthldz_mo, dqtdz_mo

  real :: ustar_mo  = 0.      !< friction velocity [m/s]
  real :: obl_mo    = -0.1    !< Obukhov length [m] (also the iteration's first guess)
  real :: dthldz_mo = 0.      !< MO surface-layer thl gradient at zf(kb) [K/m]
  real :: dqtdz_mo  = 0.      !< MO surface-layer qt gradient at zf(kb) [1/m]
  real, parameter :: min_horv = 0.1   !< minimum wind speed [m/s] (DALES)

contains

  logical function lmo_surface()
    use modglobal, only : BCbotm, BCbotm_mo, BCbotm_mo_ustar
    lmo_surface = (BCbotm == BCbotm_mo) .or. (BCbotm == BCbotm_mo_ustar)
  end function lmo_surface

  subroutine surface_mo
    use modglobal,  only : ib, ie, jb, je, kb, itot, jtot, zf, dzf, dzfi, dzhi, dzhiq, &
                           dxf, dxhi, fkar, rv, rd, lmoist, BCbotm, BCbotm_mo_ustar
    use modfields,  only : u0, v0, up, vp, thl0, qt0, u0av, v0av, thl0av, qt0av, momfluxb
    use modsurfdata, only : wtsurf, wqsurf, thls, qts, thvs, z0, z0h, ustin
    use modsubgriddata, only : ekm
    use modmpi,     only : comm3d, mpierr, my_real, mpi_sum
    implicit none
    integer :: i, j
    real :: z1, horvav, cm, cs, scaling, upcu, vpcv, horv, emom, utan2, bcmom
    real :: loc(2), glob(2), wthl_up, wqt_up

    z1 = zf(kb)
    ! DALES upward-positive kinematic fluxes; uDALES prescribes them positive in -z (D1)
    wthl_up = -wtsurf
    wqt_up  = -wqsurf

    ! drag coefficients (DALES calc_drag_coefficients, uniform here). DALES order for
    ! isurf 3/4: drag coefficients BEFORE get_obl, i.e. with the previous obl.
    cm = fkar**2 / (log(z1/z0) - psim(z1/obl_mo) + psim(z0/obl_mo))**2
    cs = fkar**2 / (log(z1/z0) - psim(z1/obl_mo) + psim(z0/obl_mo)) &
                 / (log(z1/z0h) - psih(z1/obl_mo) + psih(z0h/obl_mo))

    ! Obukhov length from the slab-mean first-level state and current surface values
    obl_mo = calc_obl_iter(thl0av(kb), qt0av(kb), thls, qts, z1, z0, z0h, u0av(kb), v0av(kb), obl_mo)

    ! friction velocity (DALES calc_friction_velocity / presc_friction_velocity)
    horvav = max(sqrt(u0av(kb)**2 + v0av(kb)**2), min_horv)
    if (BCbotm == BCbotm_mo_ustar) then
      ustar_mo = ustin
    else
      ustar_mo = sqrt(cm) * horvav
    end if
    ustar_mo = max(ustar_mo, 1.e-2)

    ! surface-layer gradients (DALES calc_surface_gradients; fluxes upward-positive)
    scaling   = 1. / (fkar * z1)
    dthldz_mo = - wthl_up / ustar_mo * phih(z1/obl_mo) * scaling
    if (lmoist) then
      dqtdz_mo = - wqt_up / ustar_mo * phih(z1/obl_mo) * scaling
    else
      dqtdz_mo = 0.
    end if

    ! surface scalars (DALES calc_surface_scalars): diagnose tskin/qskin per column, average
    loc = 0.
    do j = jb, je
      do i = ib, ie
        upcu = 0.5 * (u0(i,j,kb) + u0(i+1,j,kb))
        vpcv = 0.5 * (v0(i,j,kb) + v0(i,j+1,kb))
        horv = max(sqrt(upcu**2 + vpcv**2), min_horv)
        loc(1) = loc(1) + min(max(wthl_up / (cs * horv), -10.), 10.) + thl0(i,j,kb)
        loc(2) = loc(2) + min(max(wqt_up / (cs * horv), -5.e-2), 5.e-2) + qt0(i,j,kb)
      end do
    end do
    call MPI_ALLREDUCE(loc, glob, 2, MY_REAL, MPI_SUM, comm3d, mpierr)
    thls = glob(1) / (real(itot) * real(jtot))
    if (lmoist) qts = glob(2) / (real(itot) * real(jtot))
    thvs = thls * (1. + (rv/rd - 1.) * qts)

    ! momentum: replace the resolved diffusive flux through the bottom face by
    ! u'w' = -u*^2 u/|U| (pattern of modwallfunctions:wfmneutral, case 91)
    do j = jb, je
      do i = ib, ie
        utan2 = 0.25 * (v0(i,j,kb) + v0(i-1,j,kb) + v0(i,j+1,kb) + v0(i-1,j+1,kb))
        horv  = sqrt(max(min_horv**2, u0(i,j,kb)**2 + utan2**2))
        bcmom = ustar_mo**2 * u0(i,j,kb) / horv
        momfluxb(i,j,kb) = momfluxb(i,j,kb) + bcmom * dzfi(kb)
        emom = (dzf(kb-1)*(ekm(i,j,kb)*dxf(i-1) + ekm(i-1,j,kb)*dxf(i)) + &
                dzf(kb)*(ekm(i,j,kb-1)*dxf(i-1) + ekm(i-1,j,kb-1)*dxf(i)))*dxhi(i)*dzhiq(kb)
        up(i,j,kb) = up(i,j,kb) + (u0(i,j,kb) - u0(i,j,kb-1))*emom*dzhi(kb)*dzfi(kb) - bcmom*dzfi(kb)
      end do
    end do
    do j = jb, je
      do i = ib, ie
        utan2 = 0.25 * (u0(i,j,kb) + u0(i,j-1,kb) + u0(i+1,j-1,kb) + u0(i+1,j,kb))
        horv  = sqrt(max(min_horv**2, utan2**2 + v0(i,j,kb)**2))
        bcmom = ustar_mo**2 * v0(i,j,kb) / horv
        momfluxb(i,j,kb) = momfluxb(i,j,kb) + bcmom * dzfi(kb)
        emom = (dzf(kb-1)*(ekm(i,j,kb) + ekm(i,j-1,kb)) + dzf(kb)*(ekm(i,j,kb-1) + ekm(i,j-1,kb-1)))*dzhiq(kb)
        vp(i,j,kb) = vp(i,j,kb) + (v0(i,j,kb) - v0(i,j,kb-1))*emom*dzhi(kb)*dzfi(kb) - bcmom*dzfi(kb)
      end do
    end do
  end subroutine surface_mo

  !> Port of DALES dev 075ac4e modsurface.f90:calc_obl_iter (Newton iteration on Rib).
  real function calc_obl_iter(thl, qt, tskin, qskin, z, z0m, z0h, u, v, Lin) result(L)
    use modglobal, only : grav, rv, rd
    real, intent(in) :: thl, qt, tskin, qskin, z, z0m, z0h, u, v, Lin
    real :: horv2, rib, thv, thvsl, fx, fxdif, lold, lstart, lend
    integer :: iter
    L = Lin
    thv   = thl   * (1. + (rv/rd - 1.) * qt)
    thvsl = tskin * (1. + (rv/rd - 1.) * qskin)
    horv2 = max(u**2 + v**2, min_horv**2)
    rib   = grav / thvsl * z * (thv - thvsl) / horv2
    if (rib == 0.) then
      L = 1.e6
      return
    end if
    iter = 0
    lold = 0.
    if (rib * L < 0. .or. abs(L) == 1.e5) then
      if (rib > 0.) L =  0.01
      if (rib < 0.) L = -0.01
    end if
    do while (abs((L - lold) / L) > 1.e-4 .and. iter < 1000)
      iter   = iter + 1
      lold   = L
      lstart = L - 0.001 * L
      lend   = L + 0.001 * L
      fx     = rib - rib_from_obl(z, L, z0h, z0m)
      fxdif  = (- rib_from_obl(z, lstart, z0h, z0m) + rib_from_obl(z, lend, z0h, z0m)) / (lstart - lend)
      L      = L - fx / fxdif
      if (rib * L < 0. .or. abs(L) == 1.e5) then
        if (rib > 0.) L =  0.01
        if (rib < 0.) L = -0.01
      end if
      if (abs(L) > 1.e6) L = sign(1.e6, L)
    end do
  end function calc_obl_iter

  elemental real function rib_from_obl(z, L, z0h, z0m)
    real, intent(in) :: z, L, z0h, z0m
    rib_from_obl = z / L * (log(z/z0h) - psih(z/L) + psih(z0h/L)) &
                         / (log(z/z0m) - psim(z/L) + psim(z0m/L))**2
  end function rib_from_obl

  ! DALES dev 075ac4e modsurface.f90: psim, psih (Businger-Dyer / Beljaars-Holtslag), phim, phih
  elemental real function psim(zeta)
    real, intent(in) :: zeta
    real :: x
    if (zeta <= 0.) then
      x    = (1. - 16. * zeta) ** 0.25
      psim = 3.14159265 / 2. - 2. * atan(x) + log((1. + x)**2 * (1. + x**2) / 8.)
    else
      psim = -2./3. * (zeta - 5./0.35) * exp(-0.35 * zeta) - zeta - (10./3.) / 0.35
    end if
  end function psim

  elemental real function psih(zeta)
    real, intent(in) :: zeta
    real :: x
    if (zeta <= 0.) then
      x    = (1. - 16. * zeta) ** 0.25
      psih = 2. * log((1. + x**2) / 2.)
    else
      psih = -2./3. * (zeta - 5./0.35) * exp(-0.35 * zeta) - (1. + (2./3.) * zeta) ** 1.5 &
             - (10./3.) / 0.35 + 1.
    end if
  end function psih

  elemental real function phih(zeta)
    real, intent(in) :: zeta
    if (zeta < 0.) then
      phih = (1. - 16. * zeta) ** (-0.50)
    else if (zeta < 1.) then
      phih = 1. + 5. * zeta
    else
      phih = 6.
    end if
  end function phih

end module modsurface
