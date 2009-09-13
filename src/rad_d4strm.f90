!> \file rad_d4strm.f90
!!  Delta 4-stream radiation

!>
!!  Delta 4-stream radiation
!>
!!  \author Robert Pincus
!!  \author Bjorn Stevens
!!  \author Thijs Heus
!!  \todo Documentation
!!  \par Revision list
!----------------------------------------------------------------------------
! This file is part of DALES.
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
! Copyright 1999-2008, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module fuliou

  use cldwtr, only : init_cldwtr, cloud_water
  use rad_solver, only : qft, nv,nv1,mb,totalpower
  use modglobal, only : pi,grav,rd,ep2
  use RandomNumbers
  use ckd

  implicit none

  logical, save :: Initialized = .False.
  type(randomNumberSequence), save :: randoms
  real, parameter :: minSolarZenithCosForVis = 1.e-4

contains
  !> Subroutine rad_init initialize data arrays for gases, ice model and water
  !> model on first call
  subroutine rad_init

    if (.not.Initialized) then
       call init_ckd
       randoms = new_RandomNumberSequence(1)
       call init_cldwtr
       Initialized = .True.
    end if

    Initialized = .True.
  end subroutine rad_init
  !> Subroutine rad: Computes radiative fluxes using a band structure
  !> defined by input ckd file
  subroutine rad (as, u0, ss, pts, ee, pp, pt, ph, po, fds, fus, fdir, fuir, &
       plwc, pre, useMcICA )

    real, intent (in)  :: pp (nv1) ! pressure at interfaces

    real, dimension(nv), intent (in)  :: &
         pt,   & ! temperature [K] at mid points
         ph,   & ! humidity mixing ratio in kg/kg
         po      ! ozone mixing ratio

    real, optional, dimension(nv), intent (in)  :: &
         plwc, & ! cloud liquid water content [g/m^3]
         pre!,  & ! effective radius of cloud droplets [microns]
!          piwc, & ! cloud ice water content [g/m^3]
!          pde,  & ! effective diameter of ice particles [microns]
!          prwc, & ! rain water content [g/m^3]
!          pgwc    ! graupel water content

    real, intent (in) :: &
         as, & ! broadband albedo (all visible bands given this value)
         ee, & ! broadband surface emissivity (all IR bands given this value)
         u0, & ! cosine of solar zenith angle
         ss, & ! Solar constant
         pts   ! Surface skin temperature

    logical, optional, intent (in ) :: useMcICA

    real, dimension(nv1), intent (out)::  &
         fds, fus,  & ! downward and upward solar flux
         fdir, fuir   ! downward and upward ir flux

    logical            :: McICA = .False.

    if(present(useMcICA)) McICA = useMcICA

    call rad_ir(pts, ee, pp, pt, ph, po, fdir, fuir, &
                 plwc, pre, McICA  )
    call rad_vis(as, u0, ss, pp, pt, ph, po, fds, fus,  &
                 plwc, pre, McICA  )

  end subroutine rad

  !> Subroutine rad_ir
  !> Computes IR radiative fluxes using a band structure
  !> defined by input ckd file
  !>
  subroutine rad_ir (pts, ee, pp, pt, ph, po, fdir, fuir, &
       plwc, pre, useMcICA  )

    real, intent (in)  :: pp (nv1) ! pressure at interfaces

    real, dimension(nv), intent (in)  :: &
         pt,   & ! temperature [K] at mid points
         ph,   & ! humidity mixing ratio in kg/kg
         po      ! ozone mixing ratio

    real, optional, dimension(nv), intent (in)  :: &
         plwc, & ! cloud liquid water content [g/m^3]
         pre!,  & ! effective radius of cloud droplets [microns]
!          piwc, & ! cloud ice water content [g/m^3]
!          pde,  & ! effective diameter of ice particles [microns]
!          prwc, & ! rain water content [g/m^3]
!          pgwc    ! graupel water content

    real, intent (in) :: &
         ee, & ! broadband surface emissivity (all IR bands given this value)
         pts   ! Surface skin temperature

    logical, optional, intent (in ) :: useMcICA

    real, dimension(nv1), intent (out)::  &
         fdir, fuir   ! downward and upward ir flux

    ! ----------------------------------------
    logical            :: McICA = .False.
    logical, parameter :: irWeighted = .False.

    real, dimension (nv)   :: tw,ww,tg,dz,tauNoGas, wNoGas, Tau, w
    real, dimension (nv1)  :: fu1, fd1, bf
    real, dimension (nv,4) :: www, pfNoGas, pf

    integer :: ib, ig, k, ig1, ig2, ibandloop, iblimit
    real :: fuq2, xir_norm
    real, dimension(:), allocatable, save :: bandWeights
    real :: randomNumber
    ! ----------------------------------------

    ib=0;ig2=0
    if (.not.Initialized) then
       call init_ckd
       randoms = new_RandomNumberSequence(1)
       call init_cldwtr
       Initialized = .True.
    end if

    if(.not. allocated(bandweights)) then
      allocate(bandweights(size(ir_bands)))
      call computeIRBandWeights(ir_bands, irWeighted, bandWeights)
    end if
    if(present(useMcICA)) McICA = useMcICA

    fdir(:) = 0.0; fuir(:) = 0.0

    call thicks(pp, pt, ph, dz)

    if (McICA) then
       !
       ! Select a single band and g-point (ib, ig1) and use these as the limits
       !   in the loop through the spectrum below.
       !
       randomNumber = getRandomReal(randoms)
       call select_bandg(ir_bands, bandweights, randomNumber, ib, ig1)
       ig2 = ig1
       iblimit = 1
    else
       iblimit = size(ir_bands)
    end if

    bandLoop: do ibandloop = 1, iblimit
      if (.not. McICA) then
         ib  = ibandloop
         ig1 = 1
         ig2 = kg(ir_bands(ib))
      end if
      !
      ! Water vapor continuum optical depth
      !
      call gascon ( center(ir_bands(ib)), pp, pt, ph, TauNoGas )
      wNoGas = 0.; pfNoGas  = 0.
      if (present(plwc)) then
        call cloud_water(ib + size(solar_bands), pre, plwc, dz, tw, ww, www)
        call combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, tw, ww, www)
      end if

      call planck(pt, pts, llimit(ir_bands(ib)), rlimit(ir_bands(ib)), bf)

      gPointLoop: do ig = ig1, ig2
         tau = TauNoGas; w = wNoGas; pf = pfNoGas
         call gases (ir_bands(ib), ig, pp, pt, ph, po, tg )
         call combineOpticalProperties(tau, w, pf, tg)

         !
         ! Solver expects cumulative optical depth
         !
         do k = 2, nv
           tau(k) = tau(k) + tau(k - 1)
         end do
         call qft (.False., ee, 0., 0., bf, tau, w, pf(:, 1), pf(:, 2),      &
              pf(:, 3), pf(:, 4), fu1, fd1)

         if (McICA) then
            xir_norm = 1./bandweights(ib)
         else
            xir_norm = gPointWeight(ir_bands(ib), ig)
         end if

         fdir(:) = fdir(:) + fd1(:) * xir_norm
         fuir(:) = fuir(:) + fu1(:) * xir_norm
      end do gPointLoop
    end do bandLoop
    !
    ! fuq2 is the surface emitted flux in the band 0 - 280 cm**-1 with a
    ! hk of 0.03.
    !
    fuq2 = bf(nv1) * 0.03 * pi * ee
    fuir(:) = fuir(:) + fuq2
  end subroutine rad_ir
  !> Subroutine rad_vis: Computes radiative fluxes using a band structure
  !> defined by input ckd file
  subroutine rad_vis (as, u0, ss, pp, pt, ph, po, fds, fus,  &
       plwc, pre, useMcICA  )

    real, intent (in)  :: pp (nv1) ! pressure at interfaces

    real, dimension(nv), intent (in)  :: &
         pt,   & ! temperature [K] at mid points
         ph,   & ! humidity mixing ratio in kg/kg
         po      ! ozone mixing ratio

    real, optional, dimension(nv), intent (in)  :: &
         plwc, & ! cloud liquid water content [g/m^3]
         pre!,  & ! effective radius of cloud droplets [microns]
!          piwc, & ! cloud ice water content [g/m^3]
!          pde,  & ! effective diameter of ice particles [microns]
!          prwc, & ! rain water content [g/m^3]
!          pgwc    ! graupel water content

    real, intent (in) :: &
         as, & ! broadband albedo (all visible bands given this value)
         u0, & ! cosine of solar zenith angle
         ss    ! Solar constant

    logical, optional, intent (in ) :: useMcICA

    real, dimension(nv1), intent (out)::  &
         fds, fus    ! downward and upward solar flux

    ! ----------------------------------------
    logical            :: McICA = .False.
    logical, parameter :: solarWeighted = .false.

    real, dimension (nv)   :: tw,ww,tg,tgm,dz, tauNoGas, wNoGas, tau, w
    real, dimension (nv1)  :: fu1, fd1, bf
    real, dimension (nv,4) :: www, pfNoGas, pf

    real, dimension(:), allocatable, save :: bandWeights

    integer :: ib, ig, k, ig1, ig2, ibandloop, iblimit
    real    :: fuq1, xs_norm
    real    :: randomNumber
    ! ----------------------------------------
    ig2=0
    if (.not.Initialized) call rad_init

    if (.not. allocated(bandweights)) then
      allocate(bandweights(size(solar_bands)))
      call computeSolarBandWeights(solar_bands, solarWeighted, bandWeights)
    end if
    if(present(useMcICA)) McICA = useMcICA

    fds(:)  = 0.0
    fus(:)  = 0.0
    bf(:)   = 0.0

    if(u0 > minSolarZenithCosForVis) then
      call thicks(pp, pt, ph, dz)

      if (McICA) then
         randomNumber = getRandomReal(randoms)
         !
         ! Select a single band and g-point (ib, ig1) and use these as the
         ! limits in the loop through the spectrum below.
         !
         call select_bandg(solar_bands, bandweights, randomNumber, ib, ig1)
         ig2 = ig1
         iblimit = 1
      else
         iblimit = size(solar_bands)
      end if

      bandLoop: do ibandloop =  1, iblimit
         !
         ! select g points either all, or one depending on McICA
         !
         if (.not. McICA) then
           ib  = ibandloop
           ig1 = 1
           ig2 = kg(solar_bands(ib))
         end if

         !
         ! Rayleigh scattering
         !
         call rayle ( ib, u0, power(solar_bands(ib)), pp, pt, dz, tauNoGas, &
              wNoGas, pfNoGas)
         !
         ! Water vapor continuum
         !
         call gascon ( center(solar_bands(ib)), pp, pt, ph, tgm )
         if(any(tgm > 0.)) &
           call combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, tgm)
         !
         ! Cloud water
         !
         if (present(plwc)) then
           call cloud_water(ib, pre, plwc, dz, tw, ww, www)
           call combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, tw,ww,www)
         end if

         gPointLoop: do ig =  ig1, ig2
           tau = tauNoGas; w = wNoGas; pf = pfNoGas
           call gases (solar_bands(ib), ig, pp, pt, ph, po, tg )
           call combineOpticalProperties(tau, w, pf, tg)

           !
           ! Solver expects cumulative optical depth
           !
           do k = 2, nv
             tau(k) = tau(k) + tau(k - 1)
           end do
           call qft (.true., 0., as, u0, bf, tau, w, pf(:, 1), pf(:, 2),    &
                pf(:, 3), pf(:, 4), fu1, fd1)
           if (McICA) then
              xs_norm =power(solar_bands(ib))/ bandweights(ib)
           else
              xs_norm =gPointWeight(solar_bands(ib), ig)*power(solar_bands(ib))
           end if
           fds(:) = fds(:) + fd1(:) * xs_norm
           fus(:) = fus(:) + fu1(:) * xs_norm
         end do gPointLoop
      end do bandLoop
      !
      ! In this model, we used the solar spectral irradiance determined by
      ! Thekaekara (1973), and 1340.0 W/m**2 is the solar energy contained
      ! in the spectral region 0.2 - 4.0 um., thus scale solar fluxes by
      ! fuq1
      !
      fuq1 = ss / totalpower
      fds(:)  = fds(:)*fuq1
      fus(:)  = fus(:)*fuq1
    end if
  end subroutine rad_vis
  !> Subroutine select_bandg
  !>
  !> selects the band (i) and the g point (j) based on the probability a
  !> photon would be found in the wavelengths covered by the band and g
  !> point, which is given by g_prob.  Note g_prob sums to unity for both
  !> the solar bands (power > 0.) and the infrared bands respectively.
  !>
  subroutine select_bandg(bands, bandweights, randomNumber, i, j)
    type(band_properties), &
          dimension(:),    &
             intent ( in) :: bands
    real, dimension(:), &
             intent ( in) :: bandweights
    real,    intent ( in) :: randomNumber
    integer, intent (out) :: i, j

    real :: cumulative

    i=1; j=1
    ! The probability contained in the first g point of the first band
    cumulative = gPointWeight(bands(i), j) * bandweights(i)

    do while (randomNumber > cumulative .and. cumulative < 1.0)
       j = j+1
       if (j > kg(bands(i)) ) then
          i=i+1
          j=1
       end if
       cumulative = cumulative + gPointWeight(bands(i), j) * bandweights(i)
    end do
  end subroutine select_bandg

  !> subroutine thicks: Integrates the hydrostatic equation to provide
  !> layer thicknesses
  !>
  subroutine thicks(pp, pt, ph, dz)

    real, intent (in) :: pp(nv1), pt(nv), ph(nv)
    real, intent (out):: dz(nv)

    integer :: i
    real    :: tv

    do  i = 1, nv
       tv = pt(i)*(1+0. + ep2*ph(i) )
       dz(i) = (rd/grav) * tv * alog( pp(i+1) / pp(i) )
    end do

  end subroutine thicks

    !> Adds optical properties to running sum
    !>   If ssa and/or w[1-4] are not present we assume the new medium is
    !> strictly absorbring
    !>
  subroutine combineOpticalProperties(tau,      ssa,      pF, &
                                      tauToAdd, ssaToAdd, pFtoAdd)
    real, dimension(:),    intent(inout) :: tau, ssa
    real, dimension(:, :), intent(inout) :: pF   ! Phs function (level, moment)
    real, dimension(:),    intent(in)    :: tautoAdd
    real, dimension(:),    optional, intent(in) :: ssaToAdd
    real, dimension(:, :), optional, intent(in) :: pFToAdd ! Phs function

    integer :: j
    if(present(ssaToAdd) .and. present(pfToAdd)) then
       do j = 1, size(pF, 2)
          where (ssa(:) * tau(:) + ssaToAdd(:) * tauToAdd(:) > 0.)
             pf(:, j) = (ssa(:)*tau(:)*pf(:, j) + ssaToAdd(:)*tauToAdd(:)     &
                  * pfToAdd(:, j))/(ssa(:)*tau(:) + ssaToAdd(:) * tauToAdd(:))
          elsewhere
             pf(:, j) = 0.
          end where
       end do
       where (tau(:) + tauToAdd(:) > 0.)
          ssa(:) = (ssa(:) * tau(:) + ssaToAdd(:) * tauToAdd(:)) /            &
               (tau(:) + tauToAdd(:))
       elsewhere
          ssa(:) = 0.
       end where
       tau(:) = tau(:) + tauToAdd(:)
    else
      !
      ! New medium is absorbing - phase function doesn't change
      !
       ssa(:) = (ssa(:) * tau(:)) / (tau(:) + tauToAdd(:))
       tau(:) = tau(:) + tauToAdd(:)
    end if

  end subroutine combineOpticalProperties
  !> Subroutine rayle:  computes optical properties associated with rayleigh
  !> scattering
  !>
  !> ri is the coefficient in Eq.(4.8) of Fu (1991) to compute the optical
  !> depth due to Rayleigh scattering in the solar bands.
  !>
  !> tr, wr, and wwr are the optical depth, single scattering albedo,
  !> and expansion coefficients of the phase function ( 1, 2, 3, and
  !> 4 ) due to the Rayleigh scattering for a given layer.
  !>
  subroutine rayle ( ib, u0, power, pp, pt, dz, tr, wr, wwr)
    integer, intent (in) :: ib
    real, intent (in)    :: u0, power, pp(nv1), pt(nv), dz(nv)
    real, intent (out)   :: tr(nv), wr(nv), wwr(nv,4)

    real, parameter :: ri(6)=(/ 0.9022e-5, 0.5282e-6, 0.5722e-7, &
         0.1433e-7, 0.4526e-8, 0.1529e-8 /)

    integer :: i
    real    :: x

    if ( ib == 1 ) then
       x = -3.902860e-6*u0*u0 + 6.120070e-6*u0 + 4.177440e-6
    else
       x = ri(ib)
    endif

    if(power > 0.) then
      do  i = 1, nv
        tr(i) = x * ( pp(i) + pp(i+1) ) * dz(i) * 0.5 / pt(i)
      end do
      wr(:) = 1.0
      wwr(:, :) = 0.
      wwr(:, 2) = 0.5
    else
      tr(:)     = 0.
      wr(:)     = 0.
      wwr(:, :) = 0.
    end if
  end subroutine rayle

  !> tgm(nv) are the optical depthes due to water vapor continuum absorp-
  !> tion in nv layers for a given band ib. We include continuum absorp-
  !> tion in the 280 to 1250 cm**-1 region. vv(11)-vv(17) are the central
  !> wavenumbers of each band in this region.
  subroutine gascon ( center, pp, pt, ph, tgm)
    real,  intent (in) :: center, pp(nv1), pt(nv), ph(nv)
    real, intent (out) :: tgm(nv)

    integer :: k
    real    :: ff, pe, s, pmid

    if ( center >= 280 .and. center <= 1250.) then
       s = ( 4.18 +  5577.8 * exp ( - 0.00787 * center ) ) / 1013.25
       do k = 1, nv
          pmid   = (pp(k) + pp(k+1))*0.5
          pe     = pmid * ph(k) / ( 0.622 + 0.378 * ph(k) )
          ff     = s*(pe +  0.002*pmid ) *exp (1800.0/pt(k) -6.08108)
          tgm(k) = ff * ph(k) * ( pp(k+1) - pp(k) ) * 1.019767
       end do
    else
       tgm(:) = 0.
    end if
  end subroutine gascon
  !> Subroutine planck:  Integrates planck function over band in xk sub
  !> intervals.  The temperatures at the interfaces are taken as the
  !> average of the mid-point temperatures, the temperature at the lowest
  !> pressure interface is set to the temperature at the first mid point,
  !> and surface temperatures are taken as the skin temperature.
  !>
  subroutine planck ( pt, tskin, llimit, rlimit, bf)
    real, intent (in)    :: pt(nv), tskin, llimit, rlimit
    real, intent (out)   :: bf(nv1) ! intensity [W/m^2/Sr]

    real, parameter :: xk = 10.

    integer :: k
    real    :: v1, v2, vmid, fq1, fq2, tk

    do k = 1, nv1
       bf(k) = 0.0
    end do

    v1 = llimit
    do while (v1 > rlimit+epsilon(rlimit))
       v2 = max(v1 - xk, rlimit)
       vmid = ( v1 + v2 ) * 0.5
       fq1 = 1.19107e-8 * vmid * vmid * vmid
       fq2 = 1.43884 * vmid
       do k = 2, nv
          tk = (pt(k)+pt(k-1))*0.5
          bf(k) = bf(k) + (fq1/(exp(fq2/tk) - 1.0))*(v1-v2)
       end do
       bf(1) = bf(1) + (fq1/(exp(fq2/pt(1)) - 1.0))*(v1-v2)
       bf(nv1) = bf(nv1) + (fq1/(exp(fq2/tskin) - 1.0))*(v1-v2)
       v1 = v2
    end do

  end subroutine planck
  !
    !>
    !> find the weighting for band points so that the probability of a photon
    !> existing in the g-point range of a band can be calculated, and used for
    !> McICA calculations.  This is the relative band width for the IR bands
    !>
  subroutine computeIRBandWeights(bands, weighted, bandweights)
    type(band_properties), &
          dimension(:), intent( in) :: bands
    logical,            intent( in) :: weighted
    real, dimension(:), intent(out) :: bandweights

    integer :: ib

    if(size(bands) /= size(bandweights)) &
      stop "Didn't provide the right amount of storage for band weights"
    if(any(isSolar(bands))) stop "Can't compute IR band weights for solar bands."

    if (weighted) then
       do ib = 1, size(bands)
         bandweights(ib) = (llimit(bands(ib)) - rlimit(bands(ib)))/(bllmx-brlmn)
       end do
    else
       bandweights(:) = 1./(real(size(bands)))
    end if
  end subroutine computeIRBandWeights
    !>
    !> find the weighting for band points so that the probability of a photon
    !> existing in the g-point range of a band can be calculated, and used for
    !> McICA calculations.  This is the solar power for the solar bands
    !>
  subroutine computeSolarBandWeights(bands, weighted, bandweights)
    type(band_properties), &
          dimension(:), intent( in) :: bands
    logical,            intent( in) :: weighted
    real, dimension(:), intent(out) :: bandweights

    integer :: i
    i=1
        if(size(bands) /= size(bandweights)) &
         stop "Didn't provide the right amount of storage for band weights"
    if(any(.not. isSolar(bands))) stop "Can't compute solar band weights in IR"

    if(weighted) then
       bandweights(:) = (/ (power(bands(i))/totalpower, i = 1, size(bands)) /)
    else
       bandweights(:) = 1./(real(size(bands)))
    end if
  end subroutine computeSolarBandWeights

end module fuliou
