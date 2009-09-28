!> \file rad_cldwtr.f90
!>  The cloudwater radiation

!>
!>  The cloudwater radiation
!>
!>  \author Robert Pincus
!>  \author Bjorn Stevens
!>  \author Thijs Heus
!>  \todo Documentation
!>  \par Revision list
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
module cldwtr

  use rad_solver, only : nv, mb
  integer, save :: nsizes = 8
  logical, save :: Initialized = .False.

  real, allocatable    :: re(:), fl(:), bz(:,:), wz(:,:), gz(:,:)

contains
  !> Surbourine cloud_init initialize data arrays for the cloud model,
  !> checking for consistency between band structure of cloud model and CKD
  !>
  subroutine init_cldwtr
    use modglobal, only : cexpnr
    use ckd, only : band, center
    integer, parameter  :: nrec = 21600

    real, dimension(mb) :: cntrs

    integer             :: ib, i, nbands=18, ierr
    character (len=12)  :: frmt
    character (len=20)  :: filenm

    filenm = 'cldwtr.inp.'//cexpnr
    open ( unit = 71, file = filenm, status = 'old', recl=nrec,iostat=ierr)
    if (ierr==0) read (71,'(2I3)') nsizes, nbands
    if (nbands /= mb .or. nsizes*nbands*15 > nrec) &
         stop 'TERMINATING: incompatible cldwtr.dat file'

    allocate (re(nsizes),fl(nsizes),bz(nsizes,mb),wz(nsizes,mb),gz(nsizes,mb))
    write(frmt,'(A1,I2.2,A8)') '(',mb,'E15.7)    '
    if (ierr/=0) call initvar(cntrs,re,fl,bz,wz,gz)
    if (ierr==0) read (71,frmt) (cntrs(i), i=1,mb)
    do i=1,mb
       if (spacing(1.) < abs(cntrs(i)- center(band(i))) ) &
            stop 'TERMINATING: cloud properties not matched to band structure'
    end do

    write(frmt,'(A1,I2.2,A9)') '(',nsizes,'E15.7)   '
    if (ierr==0) read (71,frmt) (re(i), i=1,nsizes)
    if (ierr==0) read (71,frmt) (fl(i), i=1,nsizes)

    write(frmt,'(A1,I4.4,A7)') '(',nsizes*mb,'E15.7) '
    if (ierr==0) read (71,frmt) ((bz(i,ib), i=1,nsizes), ib=1,mb)
    if (ierr==0) read (71,frmt) ((wz(i,ib), i=1,nsizes), ib=1,mb)
    if (ierr==0) read (71,frmt) ((gz(i,ib), i=1,nsizes), ib=1,mb)
    if (ierr==0) close (71)

    if (minval((/bz,wz,gz/)) < 0.) &
         stop 'TERMINATING: cloud properties out of bounds'

    Initialized = .True.

  end subroutine init_cldwtr

  !> Subroutine cloud_water:  calculates the optical depth (tw), single
  !> scattering albedo (ww), and phase function (www(4)) given the cloud
  !> water [g/m^3] and effective radius [microns] by interpolating based on
  !> known optical properties at predefined sizes
  !>
  subroutine cloud_water ( ib, pre, pcw, dz, tw, ww, www )

    implicit none

    integer, intent (in) :: ib
    real, dimension (nv), intent (in) :: pre, pcw, dz
    real, intent (out) :: tw(nv), ww(nv), www(nv,4)

    integer :: k, j, j0, j1
    real    :: gg, wght, cwmks

    if (.not.Initialized) stop 'TERMINATING: Cloud not Initialized'

    do k = 1, nv
       cwmks = pcw(k)*1.e-3
       if ( cwmks .ge. 1.e-8) then
          j = 0
          do while (j<nsizes .and. pre(k) > re(j+1))
             j = j + 1
          end do
          if (j >= 1 .and. j < nsizes) then
             j1 = j+1
             wght = (pre(k)-re(j))/(re(j1)-re(j))
             tw(k) = dz(k) * cwmks * ( bz(j,ib) / fl(j) +   &
                  ( bz(j1,ib) / fl(j1) - bz(j,ib) / fl(j) ) /    &
                  ( 1.0 / re(j1) - 1.0 / re(j) ) * ( 1.0 / pre(k) &
                  - 1.0 / re(j) ) )
             ww(k) = wz(j,ib) + (wz(j1,ib) - wz(j,ib) ) * wght
             gg    = gz(j,ib) + (gz(j1,ib) - gz(j,ib) ) * wght
          else
             j0 = max(j,1)
             tw(k) = dz(k) * cwmks * (bz(j0,ib)/fl(j0))
             ww(k) = wz(j0,ib)
             gg    = gz(j0,ib)
          end if
          www(k,1) = 3.0 * gg
          do j=2,4
             wght = real(2*j+1)/real(2*j-1)
             www(k,j) = www(k,j-1) * gg * wght
          end do
       else
          www(k,:) = 0.0
          tw(k) = 0.0
          ww(k) = 0.0
          gg    = 0.
       end if
    end do

    return
  end subroutine cloud_water

  !> linear interpolation between two points, returns indicies of the
  !> interpolation points and weights
  !>
  subroutine interpolate(x,ny,y,i1,i2,alpha)

    integer, intent (in) :: ny
    real, intent (in)    :: x, y(ny)

    integer, intent (out) :: i1, i2
    real, intent (out)    :: alpha

    if (y(1) < y(2)) stop 'TERMINATING: band centers increasing'

    i2 = 1
    do while (x < y(i2) .and. i2 < ny)
       i2 = i2+1
    end do
    i1 = max(1,i2-1)
    alpha = 1.

    if(i2.ne.i1) alpha = (x-y(i1))/(y(i2)-y(i1))
    if (alpha <0 .or. alpha >1) print 600, x, y(1), y(ny), alpha

    return

600 format(/'CLOUD_INIT WARNING:  Extrapolating because data out of range', &
         /1x,'x = ',F8.1,', ymax = ',F7.0,', ymin =',F7.0,', alpha = ',F6.3)
  end subroutine interpolate
  subroutine initvar(cntrs,re,fl,bz,wz,gz)
  real, dimension(:),intent(out) :: cntrs,re,fl
  real, dimension(:,:),intent(out) :: bz,wz,gz
    cntrs = (/0.3225000E+05, 0.1110000E+05, 0.6475000E+04, 0.4625000E+04, 0.3425000E+04, 0.2675000E+04, 0.2200000E+04, 0.1800000E+04, 0.1550000E+04, 0.1325000E+04, 0.1175000E+04, 0.1040000E+04, 0.8900000E+03, 0.7350000E+03, 0.6050000E+03, 0.4700000E+03, 0.3400000E+03, 0.1400000E+03/)
    re = (/0.5000000E-01, 0.1400000E+00, 0.2200000E+00, 0.2800000E+00, 0.5000000E+00, 0.4700000E+00, 0.1000000E+01, 0.2500000E+01/)
    fl = (/0.5000000E-01, 0.1400000E+00, 0.2200000E+00, 0.2800000E+00, 0.5000000E+00, 0.4700000E+00, 0.1000000E+01, 0.2500000E+01/)
    bz = reshape((/0.1511000E+02, 0.4025000E+02, 0.5981000E+02, 0.7243000E+02, 0.8369000E+02, 0.7399000E+02, 0.1281700E+03, 0.1209100E+03, 0.1574000E+02, 0.4170000E+02, 0.6152000E+02, 0.7447000E+02, 0.8578000E+02, 0.7559000E+02, 0.1304600E+03, 0.1218400E+03, 0.1638000E+02, 0.4352000E+02, 0.6484000E+02, 0.7797000E+02, 0.8731000E+02, 0.7736000E+02, 0.1343000E+03, 0.1240600E+03, 0.1757000E+02, 0.4578000E+02, 0.6644000E+02, 0.8015000E+02, 0.9049000E+02, 0.7990000E+02, 0.1375600E+03, 0.1259200E+03, 0.1819000E+02, 0.4663000E+02, 0.6939000E+02, 0.8220000E+02, 0.9146000E+02, 0.7999000E+02, 0.1382100E+03, 0.1260800E+03, 0.2130000E+02, 0.5188000E+02, 0.7777000E+02, 0.8702000E+02, 0.9491000E+02, 0.8355000E+02, 0.1434600E+03, 0.1284500E+03, 0.2244000E+02, 0.5735000E+02, 0.8441000E+02, 0.1035000E+03, 0.1034900E+03, 0.8417000E+02, 0.1527700E+03, 0.1320700E+03, 0.1832000E+02, 0.5269000E+02, 0.7667000E+02, 0.1003100E+03, 0.1054600E+03, 0.9286000E+02, 0.1578200E+03, 0.1330300E+03, 0.1727000E+02, 0.5044000E+02, 0.7418000E+02, 0.9676000E+02, 0.1053200E+03, 0.9525000E+02, 0.1580700E+03, 0.1344800E+03, 0.1373000E+02, 0.4490000E+02, 0.6770000E+02, 0.9085000E+02, 0.1091600E+03, 0.1054800E+03, 0.1631100E+03, 0.1362100E+03, 0.1030000E+02, 0.3628000E+02, 0.5723000E+02, 0.7643000E+02, 0.1064500E+03, 0.1049000E+03, 0.1617300E+03, 0.1366200E+03, 0.7160000E+01, 0.2640000E+02, 0.4351000E+02, 0.5724000E+02, 0.9255000E+02, 0.9055000E+02, 0.1491000E+03, 0.1351300E+03, 0.6390000E+01, 0.2100000E+02, 0.3381000E+02, 0.4336000E+02, 0.6690000E+02, 0.6358000E+02, 0.1138300E+03, 0.1256500E+03, 0.1033000E+02, 0.3087000E+02, 0.4763000E+02, 0.6033000E+02, 0.7954000E+02, 0.7392000E+02, 0.1274600E+03, 0.1282100E+03, 0.1186000E+02, 0.3564000E+02, 0.5481000E+02, 0.6985000E+02, 0.9039000E+02, 0.8416000E+02, 0.1424900E+03, 0.1352500E+03, 0.1027000E+02, 0.3308000E+02, 0.5181000E+02, 0.6726000E+02, 0.9324000E+02, 0.8860000E+02, 0.1487100E+03, 0.1404200E+03, 0.6720000E+01, 0.2409000E+02, 0.3942000E+02, 0.5168000E+02, 0.8334000E+02, 0.8072000E+02, 0.1401400E+03, 0.1435700E+03, 0.3920000E+01, 0.1476000E+02, 0.2532000E+02, 0.3263000E+02, 0.6085000E+02, 0.5881000E+02, 0.1123000E+03, 0.1456200E+03/),(/nsizes,mb/))
    wz = reshape((/0.9999990E+00, 0.9999990E+00, 0.9999990E+00, 0.9999990E+00, 0.9999980E+00, 0.9999990E+00, 0.9999980E+00, 0.9999970E+00, 0.9997530E+00, 0.9997000E+00, 0.9996670E+00, 0.9996460E+00, 0.9994920E+00, 0.9994700E+00, 0.9993440E+00, 0.9986670E+00, 0.9959140E+00, 0.9949670E+00, 0.9943790E+00, 0.9938420E+00, 0.9913850E+00, 0.9907530E+00, 0.9889080E+00, 0.9748310E+00, 0.9837610E+00, 0.9789810E+00, 0.9765680E+00, 0.9747000E+00, 0.9634660E+00, 0.9599340E+00, 0.9538650E+00, 0.8976900E+00, 0.7029490E+00, 0.6832410E+00, 0.6797230E+00, 0.6690450E+00, 0.6426160E+00, 0.6329960E+00, 0.6297760E+00, 0.5888200E+00, 0.9473430E+00, 0.9296190E+00, 0.9248060E+00, 0.9145570E+00, 0.8771690E+00, 0.8670470E+00, 0.8536610E+00, 0.7374260E+00, 0.9193560E+00, 0.8962740E+00, 0.8859240E+00, 0.8810970E+00, 0.8127720E+00, 0.7816370E+00, 0.7754180E+00, 0.6373410E+00, 0.8747170E+00, 0.8611220E+00, 0.8478500E+00, 0.8516770E+00, 0.7871710E+00, 0.7729520E+00, 0.7531430E+00, 0.6186560E+00, 0.7647500E+00, 0.7524100E+00, 0.7365290E+00, 0.7434350E+00, 0.6712720E+00, 0.6593920E+00, 0.6394920E+00, 0.5499410E+00, 0.8075360E+00, 0.8087000E+00, 0.7959940E+00, 0.8054890E+00, 0.7505770E+00, 0.7555240E+00, 0.7094720E+00, 0.5719890E+00, 0.7533460E+00, 0.7720260E+00, 0.7672730E+00, 0.7770790E+00, 0.7512640E+00, 0.7609730E+00, 0.7125360E+00, 0.5682860E+00, 0.6327220E+00, 0.6763320E+00, 0.6846310E+00, 0.6935520E+00, 0.7079860E+00, 0.7177240E+00, 0.6824300E+00, 0.5528670E+00, 0.2888850E+00, 0.3484890E+00, 0.3716530E+00, 0.3803670E+00, 0.4545400E+00, 0.4657690E+00, 0.4754090E+00, 0.4938810E+00, 0.2618270E+00, 0.3062830E+00, 0.3213400E+00, 0.3330510E+00, 0.3929170E+00, 0.4068760E+00, 0.4174500E+00, 0.4845930E+00, 0.2958040E+00, 0.3399290E+00, 0.3524940E+00, 0.3655020E+00, 0.4162290E+00, 0.4303690E+00, 0.4352670E+00, 0.4913560E+00, 0.3012140E+00, 0.3547460E+00, 0.3693460E+00, 0.3819060E+00, 0.4336020E+00, 0.4473970E+00, 0.4474060E+00, 0.4869680E+00, 0.2437140E+00, 0.3187610E+00, 0.3446420E+00, 0.3527700E+00, 0.4279060E+00, 0.4389790E+00, 0.4459720E+00, 0.4772640E+00, 0.1090120E+00, 0.1872300E+00, 0.2268490E+00, 0.2249760E+00, 0.3313820E+00, 0.3359170E+00, 0.3748820E+00, 0.4570670E+00/),(/nsizes,mb/))
    gz = reshape((/0.8380000E+00, 0.8390000E+00, 0.8440000E+00, 0.8470000E+00, 0.8490000E+00, 0.8600000E+00, 0.8530000E+00, 0.8590000E+00, 0.8090000E+00, 0.8100000E+00, 0.8190000E+00, 0.8230000E+00, 0.8230000E+00, 0.8490000E+00, 0.8330000E+00, 0.8430000E+00, 0.7740000E+00, 0.7870000E+00, 0.7810000E+00, 0.7920000E+00, 0.8120000E+00, 0.8360000E+00, 0.8150000E+00, 0.8330000E+00, 0.8010000E+00, 0.8020000E+00, 0.7930000E+00, 0.7930000E+00, 0.8140000E+00, 0.8290000E+00, 0.8180000E+00, 0.8320000E+00, 0.8770000E+00, 0.8730000E+00, 0.8790000E+00, 0.8800000E+00, 0.8850000E+00, 0.8990000E+00, 0.8910000E+00, 0.9080000E+00, 0.7830000E+00, 0.7690000E+00, 0.7770000E+00, 0.7560000E+00, 0.7640000E+00, 0.7760000E+00, 0.7700000E+00, 0.7970000E+00, 0.8180000E+00, 0.8050000E+00, 0.8240000E+00, 0.8300000E+00, 0.8150000E+00, 0.8010000E+00, 0.8200000E+00, 0.8450000E+00, 0.8100000E+00, 0.8020000E+00, 0.8260000E+00, 0.8400000E+00, 0.8290000E+00, 0.8530000E+00, 0.8400000E+00, 0.8680000E+00, 0.7740000E+00, 0.7660000E+00, 0.7990000E+00, 0.8180000E+00, 0.8150000E+00, 0.8690000E+00, 0.8340000E+00, 0.8690000E+00, 0.7340000E+00, 0.7280000E+00, 0.7670000E+00, 0.7970000E+00, 0.7960000E+00, 0.8710000E+00, 0.8180000E+00, 0.8540000E+00, 0.6930000E+00, 0.6880000E+00, 0.7360000E+00, 0.7720000E+00, 0.7800000E+00, 0.8800000E+00, 0.8080000E+00, 0.8460000E+00, 0.6430000E+00, 0.6460000E+00, 0.6980000E+00, 0.7410000E+00, 0.7590000E+00, 0.8820000E+00, 0.7930000E+00, 0.8390000E+00, 0.5640000E+00, 0.5820000E+00, 0.6370000E+00, 0.6900000E+00, 0.7190000E+00, 0.8710000E+00, 0.7640000E+00, 0.8190000E+00, 0.4660000E+00, 0.4940000E+00, 0.5460000E+00, 0.6090000E+00, 0.6510000E+00, 0.8230000E+00, 0.7010000E+00, 0.7660000E+00, 0.3750000E+00, 0.4100000E+00, 0.4550000E+00, 0.5250000E+00, 0.5830000E+00, 0.7730000E+00, 0.6370000E+00, 0.7100000E+00, 0.2620000E+00, 0.3010000E+00, 0.3340000E+00, 0.4060000E+00, 0.4850000E+00, 0.6950000E+00, 0.5450000E+00, 0.6310000E+00, 0.1440000E+00, 0.1810000E+00, 0.2000000E+00, 0.2560000E+00, 0.3520000E+00, 0.5620000E+00, 0.4130000E+00, 0.5170000E+00, 0.6000000E-01, 0.7700000E-01, 0.8800000E-01, 0.1120000E+00, 0.1810000E+00, 0.3100000E+00, 0.2220000E+00, 0.3270000E+00/),(/nsizes,mb/))

  end subroutine

end module cldwtr
