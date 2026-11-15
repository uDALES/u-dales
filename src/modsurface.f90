!> Monin-Obukhov surface-layer kernels, ported from DALES dev 075ac4e src/surface/modsurface.f90:
!> the stability functions psim, psih, phim, phih, the iterated Obukhov length (calc_obl_iter)
!> and the per-element surface update of isurf 3/4 with lmostlocal = .true. (mo_local_update).
!> Used by the Monin-Obukhov facets (iwallmom = 4, modibm:wallfunmom).
module modsurface
  implicit none
  save
  private
  public :: psim, psih, phim, phih, calc_obl_iter, rib_from_obl, mo_local_update   ! pure kernels, tested by tests.f90:tests_mo_kernels
  public :: min_horv, ustar_min

  real, parameter :: min_horv = 0.1   !< minimum wind speed [m/s] (DALES)
  real, parameter :: ustar_min = 1.e-2   !< floor of the friction velocity [m/s] (DALES)
  real, parameter :: dtskin_max = 10.    !< bound of |tskin - thl(kb)| [K] (DALES calc_surface_scalars)
  real, parameter :: dqskin_max = 5.e-2  !< bound of |qskin - qt(kb)| [kg/kg] (DALES calc_surface_scalars)
  real, parameter :: obl_inf = 1.e6      !< |L| of a neutral surface layer and its upper bound [m] (DALES)
  ! Stability functions of DALES dev 075ac4e modsurface.f90. Unstable: Businger-Dyer with
  ! (1 - gamma_bd zeta)^(1/4) (Businger et al. 1971; Dyer 1974; Paulson 1970 integrals). Stable:
  ! Beljaars and Holtslag (1991, J. Appl. Meteorol. 30, 327) with a = 1, b = 2/3, c = 5, d = 0.35
  ! for psim/psih, and the linear phi = 1 + beta_bd zeta capped at 6 (zeta >= 1) for phim/phih.
  real, parameter :: gamma_bd = 16.      !< Businger-Dyer unstable coefficient [-]
  real, parameter :: beta_bd  = 5.       !< Businger-Dyer stable slope of phi [-]
  real, parameter :: bh_c = 5., bh_d = 0.35   !< Beljaars-Holtslag c, d [-]
  real, parameter :: phi_max = 6.        !< cap of phim/phih for zeta >= 1 [-]
  !> pi as written in DALES psim (9 digits, not modglobal pi); kept so that psim is unchanged to the bit
  real, parameter :: pi_dales = 3.14159265

contains

  !> Local Monin-Obukhov update of one surface element (a facet section), lmostlocal = .true.
  !> Port of the per-column sequence of DALES dev 075ac4e modsurface.f90:surface, isurf = 4
  !> (isurf = 3 when ustin > 0): calc_drag_coefficients with the previous L, get_obl,
  !> calc_friction_velocity (presc_friction_velocity), calc_surface_scalars.
  !> thl, qt and horv are the air state at the reference point, a height z1 above the surface;
  !> wthl_up, wqt_up the prescribed kinematic fluxes, positive upward (away from the surface).
  !> obl, tskin and qskin carry the state from the previous call.
  subroutine mo_local_update(thl, qt, horv, z1, z0m, z0h, wthl_up, wqt_up, ustin, obl, tskin, qskin, ustar)
    use modglobal, only : fkar
    real, intent(in)    :: thl, qt, horv, z1, z0m, z0h, wthl_up, wqt_up, ustin
    real, intent(inout) :: obl, tskin, qskin
    real, intent(out)   :: ustar
    real :: horvf, cm, cs

    horvf = max(horv, min_horv)
    cm = fkar**2 / (log(z1/z0m) - psim(z1/obl) + psim(z0m/obl))**2
    cs = fkar**2 / (log(z1/z0m) - psim(z1/obl) + psim(z0m/obl)) &
                 / (log(z1/z0h) - psih(z1/obl) + psih(z0h/obl))
    obl = calc_obl_iter(thl, qt, tskin, qskin, z1, z0m, z0h, horv, 0., obl)
    if (ustin > 0.) then
      ustar = ustin
    else
      ustar = sqrt(cm) * horvf
    end if
    ustar = max(ustar, ustar_min)
    tskin = min(max(wthl_up / (cs * horvf), -dtskin_max), dtskin_max) + thl
    qskin = min(max(wqt_up  / (cs * horvf), -dqskin_max), dqskin_max) + qt
  end subroutine mo_local_update

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
      L = obl_inf
      return
    end if
    iter = 0
    lold = 0.
    ! The abs(L) == 1.e5 reset is inherited from DALES, whose first guess is 1e5;
    ! it is inactive with the -0.1 first guess used here.
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
      if (abs(L) > obl_inf) L = sign(obl_inf, L)
    end do
  end function calc_obl_iter

  elemental real function rib_from_obl(z, L, z0h, z0m)
    real, intent(in) :: z, L, z0h, z0m
    rib_from_obl = z / L * (log(z/z0h) - psih(z/L) + psih(z0h/L)) &
                         / (log(z/z0m) - psim(z/L) + psim(z0m/L))**2
  end function rib_from_obl

  !> Integrated stability functions psim, psih and gradient functions phim, phih
  !> (DALES dev 075ac4e modsurface.f90), constants as documented at their declaration.
  elemental real function psim(zeta)
    real, intent(in) :: zeta
    real :: x
    if (zeta <= 0.) then
      x    = (1. - gamma_bd * zeta) ** 0.25
      psim = pi_dales / 2. - 2. * atan(x) + log((1. + x)**2 * (1. + x**2) / 8.)
    else
      psim = -2./3. * (zeta - bh_c/bh_d) * exp(-bh_d * zeta) - zeta - (10./3.) / bh_d
    end if
  end function psim

  elemental real function psih(zeta)
    real, intent(in) :: zeta
    real :: x
    if (zeta <= 0.) then
      x    = (1. - gamma_bd * zeta) ** 0.25
      psih = 2. * log((1. + x**2) / 2.)
    else
      psih = -2./3. * (zeta - bh_c/bh_d) * exp(-bh_d * zeta) - (1. + (2./3.) * zeta) ** 1.5 &
             - (10./3.) / bh_d + 1.
    end if
  end function psih

  elemental real function phim(zeta)
    real, intent(in) :: zeta
    if (zeta < 0.) then
      phim = (1. - gamma_bd * zeta) ** (-0.25)
    else if (zeta < 1.) then
      phim = 1. + beta_bd * zeta
    else
      phim = phi_max
    end if
  end function phim

  elemental real function phih(zeta)
    real, intent(in) :: zeta
    if (zeta < 0.) then
      phih = (1. - gamma_bd * zeta) ** (-0.50)
    else if (zeta < 1.) then
      phih = 1. + beta_bd * zeta
    else
      phih = phi_max
    end if
  end function phih

end module modsurface
