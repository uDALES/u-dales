!> \file rad_corkds.f90
!! Correlated k-distribution

!>
!! Correlated k-distribution
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
module ckd

  use modglobal, only : mair,cexpnr
  use rad_solver, only :nv, nv1, mb, totalpower
  implicit none
  private

  character (len=20) :: gasfile = 'ckd.dat.'
  logical, save      :: Initialized = .False.
  real, save         :: bllmx, brlmn
  integer, save      :: ngases

  TYPE ckd_properties
     character (len=5) :: name
     integer           :: ng, nt, np, noverlap, iband
     real              :: mweight, default_conc, tbase
     real, allocatable :: hk(:),sp(:),xk(:,:,:,:)
  END TYPE ckd_properties

  TYPE band_properties
     private
     integer              :: kg, ngases
     real                 :: llimit, rlimit, center, power = 0.
     integer, allocatable :: gas_id(:)
     real, allocatable    :: hk(:)
  END TYPE band_properties

  TYPE (ckd_properties),  allocatable :: gas(:)
  TYPE (band_properties), allocatable :: band(:), solar_bands(:), ir_bands(:)

  public :: band_properties, bllmx, brlmn, &
            band, solar_bands, ir_bands, &
            init_ckd, gases, &
            kg, llimit, rlimit, center, power, gPointWeight, isSolar
contains
  !
  ! ---------------------------------------------------------------------------
  !> Subroutine ckd_init:  Reads the correlated K distribution data an assures
  !> that it conforms to expected properties
  !>
  subroutine init_ckd

    implicit none

    integer :: i, j, k, l, n, ib, ii, mbs, mbir
    logical :: check
    real    :: bllmx, brlmn
    real, dimension(2000) :: realVars

    real, allocatable :: gasesinband(:)

    OPEN ( unit = 66, file = trim(gasfile//cexpnr), status = 'old' )
    read (66, '(2I5)') mb, ngases
    allocate (band(mb), gas(ngases))
    read (66, '(300(6E12.4,/))') realVars(1:mb+1)
    do ib = 1, mb
       band(ib)%llimit = realVars(ib)
    end do
    band(mb)%rlimit =  realVars(mb+1)
    read (66, '(300(6E12.4,/))') realVars(1:mb)
    do ib = 1, mb
       band(ib)%power = realVars(ib)
    end do

    mbs = 0.
    do ib=1,mb-1
       band(ib)%rlimit = band(ib+1)%llimit
       band(ib)%center =(band(ib)%rlimit+band(ib)%llimit)*0.5
       if (band(ib)%power > 0.) mbs = mbs + 1
    end do
    band(mb)%center =(band(mb)%rlimit+band(mb)%llimit)*0.5

    if (band(mb)%power > 0.) mbs = mbs + 1
    mbir = mb - mbs
    print 600, trim(gasfile), ngases, mb, mbs, mbir, sum(band%power)

    do n=1,ngases
       read (66,'(A5,I4)') gas(n)%name,gas(n)%iband
       read (66,'(4I4)') gas(n)%noverlap,gas(n)%ng,gas(n)%np,gas(n)%nt
       read (66,'(3E13.5)') gas(n)%mweight, gas(n)%default_conc, gas(n)%tbase
       allocate (gas(n)%hk(gas(n)%ng))
       allocate (gas(n)%sp(gas(n)%np))
       allocate (gas(n)%xk(gas(n)%nt,gas(n)%np,gas(n)%ng,gas(n)%noverlap))
       read (66, '(300(6E12.4,/))') gas(n)%hk(1:gas(n)%ng)
       read (66,'(300(6E12.4,/))') gas(n)%sp(1:gas(n)%np)
       read (66,'(300(6E12.4,/))')                                            &
            realVars(:(gas(n)%ng * gas(n)%np * gas(n)%nt * gas(n)%noverlap))
       ib = 1
       do l = 1, gas(n)%noverlap
         do k = 1, gas(n)%nt
           do j = 1, gas(n)%np
             do i = 1, gas(n)%ng
               gas(n)%xk(k,j,i,l) = realVars(ib)
               ib = ib + 1
             end do
           end do
         end do
       end do

       if (abs(sum(gas(n)%hk) - 1.) <= 1.1 * spacing(1.) ) then
          print 601, gas(n)%name, gas(n)%iband, gas(n)%noverlap,              &
               gas(n)%ng, gas(n)%np, gas(n)%nt
       else
          print *, gas(n)%hk, sum(gas(n)%hk(:))
          stop 'TERMINATING: gas did not occur with probability one in band'
       end if
    end do

    bllmx = tiny(1.)
    brlmn = huge(1.)
    allocate (gasesinband(ngases))
    totalpower = 0.
    do ib=1,mb
       check = .False.
       ii = 0
       do n=1,ngases
          if (gas(n)%iband == ib) then
             check = .True.
             ii = ii+1
             gasesinband(ii) = n
             if (gas(n)%ng > 1 .and. .not. allocated(band(ib)%hk)) then
                allocate(band(ib)%hk(gas(n)%ng))
                band(ib)%hk = gas(n)%hk
                band(ib)%kg = gas(n)%ng
             end if
          end if
       end do
       if (.not.check) stop 'TERMINATING: Gases do not span bands'
       band(ib)%ngases = ii
       allocate(band(ib)%gas_id(ii))
       band(ib)%gas_id(:) = gasesinband(1:ii)
       if (band(ib)%power < epsilon(1.)) then
          bllmx = max(band(ib)%llimit, bllmx)
          brlmn = min(band(ib)%rlimit, brlmn)
       end if
       totalpower = totalpower+band(ib)%power
    end do
    deallocate (gasesinband)

    !
    ! Make separate solar and IR spectra
    !
    allocate(solar_bands(mbs), ir_bands(mbir))
    i = 1; j = 1
    do ib = 1, mb
       if(isSolar(band(ib))) then
          if(i > size(solar_bands)) stop 'TERMINATING: mismatch in solar bands'
          solar_bands(i) = copy_band_properties(band(ib))
          i = i + 1
       else
          if(j > size(ir_bands))    stop 'TERMINATING: mismatch in solar bands'
          ir_bands(j) = copy_band_properties(band(ib))
          j = j + 1
       end if
    end do

    do ib=1,mb
       print 602, ib, band(ib)%power, band(ib)%llimit, band(ib)%rlimit,    &
            band(ib)%ngases, band(ib)%kg
    end do
    print 604

    Initialized = .True.

600 format ('-----------------------------------------------------------', &
         /3x,'Computing with file: ',A20,' containing ',I3,' gases',   &
         /3x,'Number of Bands (Total/Solar/IR): ',I3,'/',I3,'/',I3     &
         /3x,'Total Power in Solar Bands: ',F8.2)
601 format ('                                                -----------', &
         /3x,'Reading Gas: ',A5,', Band ',I3,', Overlap Type ',I1,     &
         /3x,'# of g-points/ P Levels/ T Coeffs :',3I3)
602 format ('-----------                                                ', &
         /3x,'Band: ',I3,': ',F8.2,' Wm^-2',', between ',F6.0,' and ',F6.0,&
         ' cm^-1',/3x,I3,' gase(s): and ',I3,' g-points')
604 format ('---------------------------------------- Finished band init ')

  end subroutine init_ckd
  !
  ! ---------------------------------------------------------------------------
  !> Subroutine gases:  given an atmospheric state, a band and g point number
  !> this routroutine calculates the optical depth at that g-point within that
  !> band for all the gases (nongray gaseous absorbers) in the band.
  !>
  subroutine gases ( this_band, ig, pp, pt, ph, po, tg )

    implicit none

    type(band_properties), &
             intent (in) :: this_band
    integer, intent (in) :: ig
    real, intent(in)     :: pp(nv1), pt(nv), ph(nv), po(nv)
    real, intent (out)   :: tg(nv)

    real, dimension (nv) :: fkg, fkga, fkgb, pq

    integer :: k, n, igg, nn
    real    :: xfct

    if (.not.Initialized) stop 'TERMINATING:  ckd_gases not initialized'
    do k = 1, nv
       tg(k) = 0.
    end do
    !
    ! loop through gases finding the cumulative probability for the band and g
    ! point, weight this by the power if in the solar bands, for which power>0.
    !
    do nn=1,this_band%ngases

       n = this_band%gas_id(nn)
       igg=min(ig,gas(n)%ng)
       select case(gas(n)%noverlap)
       case (1)
          igg=min(ig,gas(n)%ng)
          call qk (gas(n)%nt, gas(n)%np, gas(n)%sp, gas(n)%tbase,            &
               gas(n)%xk(1,1,igg,1), pp, pt, fkg )
          call select_gas(gas(n)%name, gas(n)%default_conc, gas(n)%mweight,  &
               pp, ph, po, pq)
          xfct = (2.24e4/gas(n)%mweight) * 10./9.81
          do k = 1, nv
             tg(k) = tg(k) + fkg(k)*pq(k)*(pp(k+1)-pp(k))*xfct
          end do

       case (2)
          call qk (gas(n)%nt, gas(n)%np, gas(n)%sp, gas(n)%tbase,            &
               gas(n)%xk(1,1,igg,1), pp, pt, fkga )
          call qk (gas(n)%nt, gas(n)%np, gas(n)%sp, gas(n)%tbase,            &
               gas(n)%xk(1,1,igg,2), pp, pt, fkgb )
          call select_gas(gas(n)%name, gas(n)%default_conc, gas(n)%mweight,  &
               pp, ph, po, pq)
          do k = 1, nv
             fkg(k) = fkga(k) + pq(k) * fkgb(k)
          end do
          call select_gas('  CO2', gas(n)%default_conc, gas(n)%mweight,       &
               pp, ph, po, pq)
          xfct = (2.24e4/gas(n)%mweight) * 10./9.81
          do k = 1, nv
             tg(k) = tg(k) + fkg(k)*pq(k)*(pp(k+1)-pp(k))*xfct
          end do

       case default
          stop 'TERMINATING: overlap type not supported'
       end select
    end do

  end subroutine gases
  !> Subroutine select_gas determines the mixing ratio to use in the optical
  !> depth calculation.  For simple overlaps, this amounts to selecting the
  !> correct input array, ph for water vapor, po for ozone. For fixed gases
  !> this converts a concentration to a mixing ratio given the molecular weight
  !> of the gas and some specified background concentration.  For C02 the
  !> mixing ratio is chosen so that if conc = 330.e-6, the pq*xfct in subroutine
  !> po = 0.5
  !>
  subroutine select_gas (name, conc_x, mx, pp, ph, po, pq)

    character (len=5), intent (in) :: name
    real, intent (in) :: conc_x,mx
    real, intent (in) :: pp(nv1), ph(nv), po(nv)
    real, intent (out):: pq(nv)

    integer :: k
    real    :: xx

    select case(name)
    case ('  H2O')
       pq(:) = ph(:)
    case ('   O3')
       pq(:) = po(:)
    case ('  CO2')
       xx = conc_x/330.e-6 / ((2.24e4/mx) * 10./9.81)
       pq(:) = xx
    case ('OVRLP')
       do k=1,nv
          if ( pp(k) >= 63.1 ) then
             pq(k) = ph(k)
          else
             pq(k) = 0.0
          endif
       end do
    case default
       xx = conc_x*(mx/mair)
       pq(:) = xx
    end select

    return
  end subroutine select_gas
  !
  ! ---------------------------------------------------------------------------
  !> Subroutine qk: interpolates the gasesous absorption coefficients in units
  !> of (cm-atm)**-1 to the given temperature and pressure in each layer
  !> following: ln k = a + b * ( t - tbase ) + c * ( t - tbase ) ** 2 in
  !> temperature and  linear interpolation in pressure.
  !>
  subroutine qk (nt, np, stanp, tbase, coefki, pp, pt, fkg )

    implicit none

    integer, intent (in) :: np, nt
    real, intent(in)  :: pp(nv1), pt(nv), coefki(nt,np), stanp(np), tbase
    real, intent(out) :: fkg(nv)

    integer :: i1, k
    real    :: x1, x2, y1, y2, xp, pmid

    if (nt < 3) then
       do k = 1, nv
          fkg(k) = coefki(1,1)
       end do
    else
       i1 = 1
       x1 = 0.
       xp = 0.
       do k = 1, nv
          pmid = 0.5*(pp(k)+pp(k+1))
          do while ( pmid >= stanp(i1) .and. i1 < np)
             i1 = i1 + 1
          end do
          y1 = (pt(k)-tbase)
          y2 = y1 * y1
          x2 = exp (coefki(1,i1) + coefki(2,i1)*y1 + coefki(3,i1)*y2)
          if (i1 > 1) then
             x1 = exp (coefki(1,i1-1) + coefki(2,i1-1)*y1 + coefki(3,i1-1)*y2)
             xp = stanp(i1-1)
          endif
          fkg(k) = x1 + (x2-x1)*(pmid-xp)/(stanp(i1)-xp)
       end do
    end if
    return
  end subroutine qk
  !
  ! ----------------------------------------------------------------------
  !
  elemental integer function kg(thisBand)
    type(band_properties), intent(in) :: thisBand

    kg = thisBand%kg
  end function kg
  !
  ! ----------------------------------------------------------------------
  !
  elemental real function power(thisBand)
    type(band_properties), intent(in) :: thisBand

    power = thisBand%power
  end function power
  !
  ! ----------------------------------------------------------------------
  !
  elemental real function llimit(thisBand)
    type(band_properties), intent(in) :: thisBand

    llimit = thisBand%llimit
  end function llimit
  !
  ! ----------------------------------------------------------------------
  !
  elemental real function rlimit(thisBand)
    type(band_properties), intent(in) :: thisBand

    rlimit = thisBand%rlimit
  end function rlimit
  !
  ! ----------------------------------------------------------------------
  !
  elemental real function center(thisBand)
    type(band_properties), intent(in) :: thisBand

    center = (thisBand%llimit + thisBand%rlimit) / 2.
  end function center
  !
  ! ----------------------------------------------------------------------
  !
  elemental real function gPointWeight(thisBand, ig)
    type(band_properties), intent(in) :: thisBand
    integer,               intent(in) :: ig

    gPointWeight = thisBand%hk(ig)
  end function gPointWeight
  !> function isSolar:  Returns True if the band is in the Solar region
  elemental logical function isSolar(thisBand)
    type(band_properties), intent(in) :: thisBand

    isSolar = thisBand%power > 0.
  end function isSolar
  function copy_band_properties(original) result(copy)

    type(band_properties), intent(in) :: original
    type(band_properties)             :: copy

    copy%kg     = original%kg
    copy%ngases = original%ngases
    copy%llimit = original%llimit
    copy%rlimit = original%rlimit
    copy%center = original%center
    copy%power  = original%power

    if(allocated(original%gas_id)) then
       allocate(copy%gas_id(size(original%gas_id)))
       copy%gas_id(:) = original%gas_id(:)
    end if

    if(allocated(original%hk)) then
       allocate(copy%hk(size(original%hk)))
       copy%hk(:) = original%hk(:)
    end if
  end function copy_band_properties

end module ckd
