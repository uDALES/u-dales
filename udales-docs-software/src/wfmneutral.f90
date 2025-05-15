!  This file is part of uDALES.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 2006-2021 the uDALES Team.
!
SUBROUTINE wfmneutral(hi,hj,hk,iout1,iout2,iomomflux,utang1,utang2,z0,n,ind,wforient)
   !wfmneutral
   !wf for momentum under neutral conditions
   !calculating wall function for momentum assuming neutral conditions
   !follow approach in wfuno
   !fluxes in m2/s2
   USE modglobal, ONLY : dzf,dzfi,dzh2i,dzhi,dzhiq,dy,dyi,dy2i,dyi5,dxf,dxh,dxfi,dxhi,dxh2i,ib,ie,jb,je,kb,ke,fkar,jmax,rk3step,kmax,jge,jgb
   USE modsubgriddata, ONLY:ekh, ekm
   USE modmpi, ONLY:myid
   USE initfac, ONLY:block
   USE modibmdata
   INTEGER i, j, k, jl, ju, kl, ku, il, iu, km, im, jm, ip, jp, kp

   REAL :: bcmomflux = 0. !temp storage for momentum flux
   REAL :: ctm = 0. !momentum transfer coefficient
   REAL :: dummy = 0. !for debugging
   REAL :: delta = 0. !distance from wall
   REAL :: logdz2 = 0. !log(delta/z0)**2
   REAL :: utang1Int !Interpolated 1st tangential velocity component needed for stability calculation (to T location)
   REAL :: utang2Int !Interpolated 2nd tangential velocity component needed for stability calculation (to T location)
   REAL :: fkar2 !fkar^2, von Karman constant squared
   REAL :: emmo = 0., epmo = 0., epom = 0., emom = 0., eopm = 0., eomm = 0., empo = 0.
   REAL :: umin = 0.0001 !m^2/s^2

   INTEGER, INTENT(in) :: hi !<size of halo in i
   INTEGER, INTENT(in) :: hj !<size of halo in j
   INTEGER, INTENT(in) :: hk !<size of halo in k
   REAL, INTENT(inout) :: iout1(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk) !updated prognostic tangential velocity (component1)
   REAL, INTENT(inout) :: iout2(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk) !updated prognostic tangential velocity (component2)
   REAL, INTENT(inout) :: iomomflux(ib - hi:ie + hi, jb - hj:je + hj, kb-hk:ke + hk) !a field to save the momentum flux
   REAL, INTENT(in)    :: z0
   REAL, INTENT(in)    :: utang1(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk) !tangential velocity field
   REAL, INTENT(in)    :: utang2(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk) !second tangential velocity field
   INTEGER, INTENT(in) :: n ! number of the block, used to get i,j,k-indeces
   INTEGER, INTENT(in) :: ind ! in case of y-wall (case 3x & 4x) "ind" is used for j-index, otherwise this is irrelevant
   INTEGER, INTENT(in) :: wforient !orientation of the facet see below:
   !frist digit, orientation of wall, determines iteration indices
   !second digit, if for momentum or for scalar (necessary because of staggered grid -> which variable to interpolate)
   !xlow=1,xup=2,yup=3,ylow=4,z=5
   !momentum=1
   fkar2 = fkar**2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CASES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CASES FOR MOMENTUM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SELECT CASE (wforient)
   !!!!!!!!!!!!!!!SPECIAL CASES FOR THE SURFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !can actually be made redundant and just be replaced by standard horizontal case (doesn't really matter though)
   ! SO: This code is essentially unchanged from uDALES v1, and should probably be moved out of this file in a later release.
   CASE (91) !surface momentum flux

      k = kb !
      km = k - 1 !
      il = ib
      iu = ie
      jl = jb
      ju = je

      delta = 0.5*dzf(k) !might need attention on streched grids! as well as the dzfi when updating up
      logdz2 = LOG(delta/z0)**2

      DO j = jl, ju !u component
         DO i = il, iu
            utang1Int = utang1(i, j, k)
            utang2Int = (utang2(i, j, k) + utang2(i - 1, j, k) + utang2(i, j + 1, k) + utang2(i - 1, j + 1, k))*0.25
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            ctm = fkar2/(logdz2)
            !dummy = (utang1Int**2)*ctm
            dummy = abs(utang1Int)*sqrt(utangInt)*ctm
            bcmomflux = SIGN(dummy, utang1Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
            emom = (dzf(km)*(ekm(i, j, k)*dxf(i - 1) + ekm(i - 1, j, k)*dxf(i)) + & ! dx is non-equidistant
                    dzf(k)*(ekm(i, j, km)*dxf(i - 1) + ekm(i - 1, j, km)*dxf(i)))*dxhi(i)*dzhiq(k)
            iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(i, j, km))*emom*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
         END DO
      END DO

      DO j = jl, ju !v component
         DO i = il, iu
            utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
            utang2Int = utang2(i, j, k)
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            ctm = fkar2/(logdz2)
            !dummy = (utang2Int**2)*ctm
            dummy = abs(utang2Int)*sqrt(utangInt)*ctm
            bcmomflux = SIGN(dummy, utang2Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
            eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
            iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
         END DO
      END DO

   END SELECT

END SUBROUTINE wfmneutral
