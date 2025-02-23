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
MODULE wallfunction_neutral
   IMPLICIT NONE
   SAVE

   CONTAINS

#if defined(_GPU)
      ATTRIBUTES(GLOBAL) SUBROUTINE wfmneutral_cuda(hi,hj,hk,iout1,iout2,iomomflux,utang1,utang2,z0,fkar2,n,ind,wforient)
         USE modcuda, ONLY: ib_d, ie_d, jb_d, je_d, kb_d, ke_d, dzf_d, dzfi_d, dzhi_d, dzhiq_d, dxf_d, dxhi_d, &
                            ekm_d, &
                            tidandstride
         IMPLICIT NONE

         INTEGER :: i, j, k, im, ip, jm, jp, km
         INTEGER :: tidx, tidy, tidz, stridex, stridey, stridez

         REAL :: bcmomflux !temp storage for momentum flux
         REAL :: ctm       !momentum transfer coefficient
         REAL :: dummy     !for debugging
         REAL :: delta     !distance from wall
         REAL :: logdz2    !log(delta/z0)**2
         REAL :: utang1Int !Interpolated 1st tangential velocity component needed for stability calculation (to T location)
         REAL :: utang2Int !Interpolated 2nd tangential velocity component needed for stability calculation (to T location)
         REAL :: utangInt  !Interpolated absolute tangential velocity
         REAL :: emom, eomm
         REAL :: umin = 0.0001 !m^2/s^2

         INTEGER, VALUE, INTENT(in)    :: hi
         INTEGER, VALUE, INTENT(in)    :: hj
         INTEGER, VALUE, INTENT(in)    :: hk
         REAL,           INTENT(inout) :: iout1(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d:ke_d + hk)
         REAL,           INTENT(inout) :: iout2(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d:ke_d + hk)
         REAL,           INTENT(inout) :: iomomflux(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d-hk:ke_d + hk)
         REAL,           INTENT(in)    :: utang1(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d - hk:ke_d + hk)
         REAL,           INTENT(in)    :: utang2(ib_d - hi:ie_d + hi, jb_d - hj:je_d + hj, kb_d - hk:ke_d + hk)
         REAL,    VALUE, INTENT(in)    :: z0
         REAL,    VALUE, INTENT(in)    :: fkar2
         INTEGER, VALUE, INTENT(in)    :: n
         INTEGER, VALUE, INTENT(in)    :: ind
         INTEGER, VALUE, INTENT(in)    :: wforient

         CALL tidandstride(tidx, tidy, tidz, stridex, stridey, stridez)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CASES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CASES FOR MOMENTUM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SELECT CASE (wforient)
         !!!!!!!!!!!!!!!SPECIAL CASES FOR THE SURFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         CASE (91) !surface momentum flux

            k = kb_d
            km = k - 1

            delta = 0.5*dzf_d(k)
            logdz2 = LOG(delta/z0)**2

            IF (tidz == k) THEN
               DO j = tidy, je_d, stridey
                  jp = j + 1
                  DO i = tidx, ie_d, stridex
                     im = i - 1
                     utang1Int = utang1(i, j, k)
                     utang2Int = (utang2(i, j, k) + utang2(im, j, k) + utang2(i, jp, k) + utang2(im, jp, k))*0.25
                     utangInt = max(umin, (utang1Int**2 + utang2Int**2))
                     ctm = fkar2/(logdz2)
                     dummy = abs(utang1Int)*sqrt(utangInt)*ctm
                     bcmomflux = SIGN(dummy, utang1Int)
                     iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi_d(k)
                     emom = (dzf_d(km)*(ekm_d(i, j, k)*dxf_d(im) + ekm_d(im, j, k)*dxf_d(i)) + &
                             dzf_d(k)*(ekm_d(i, j, km)*dxf_d(im) + ekm_d(im, j, km)*dxf_d(i)))*dxhi_d(i)*dzhiq_d(k)
                     iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(i, j, km))*emom*dzhi_d(k)*dzfi_d(k) - bcmomflux*dzfi_d(k)
                  END DO
               END DO
            END IF

            IF (tidz == k) THEN
               DO j = tidy, je_d, stridey
                  jm = j - 1
                  DO i = tidx, ie_d, stridex
                     ip = i + 1
                     utang1Int = (utang1(i, j, k) + utang1(i, jm, k) + utang1(ip, jm, k) + utang1(ip, j, k))*0.25
                     utang2Int = utang2(i, j, k)
                     utangInt = max(umin, (utang1Int**2 + utang2Int**2))
                     ctm = fkar2/(logdz2)
                     dummy = abs(utang2Int)*sqrt(utangInt)*ctm
                     bcmomflux = SIGN(dummy, utang2Int)
                     iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi_d(k)
                     eomm = (dzf_d(km)*(ekm_d(i, j, k) + ekm_d(i, jm, k)) + dzf_d(k)*(ekm_d(i, j, km) + ekm_d(i, jm, km)))*dzhiq_d(k)
                     iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi_d(k)*dzfi_d(k) - bcmomflux*dzfi_d(k)
                  END DO
               END DO
            END IF

         END SELECT

      END SUBROUTINE wfmneutral_cuda
#else
      SUBROUTINE wfmneutral(hi,hj,hk,iout1,iout2,iomomflux,utang1,utang2,z0,fkar2,n,ind,wforient)
         !wfmneutral
         !wf for momentum under neutral conditions 
         !calculating wall function for momentum assuming neutral conditions
         !follow approach in wfuno
         !fluxes in m2/s2
         USE modglobal,      ONLY: dzf, dzfi, dzhi, dzhiq, dxf, dxhi, ib, ie, jb, je, kb, ke
         USE modsubgriddata, ONLY: ekm
         IMPLICIT NONE

         INTEGER :: i, j, k, km, im, jm, ip, jp

         REAL :: bcmomflux !temp storage for momentum flux
         REAL :: ctm       !momentum transfer coefficient
         REAL :: dummy     !for debugging
         REAL :: delta     !distance from wall
         REAL :: logdz2    !log(delta/z0)**2
         REAL :: utang1Int !Interpolated 1st tangential velocity component needed for stability calculation (to T location)
         REAL :: utang2Int !Interpolated 2nd tangential velocity component needed for stability calculation (to T location)
         REAL :: utangInt  !Interpolated absolute tangential velocity
         REAL :: emom, eomm
         REAL :: umin = 0.0001 !m^2/s^2

         INTEGER, INTENT(in) :: hi !<size of halo in i
         INTEGER, INTENT(in) :: hj !<size of halo in j
         INTEGER, INTENT(in) :: hk !<size of halo in k
         REAL, INTENT(inout) :: iout1(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk) !updated prognostic tangential velocity (component1)
         REAL, INTENT(inout) :: iout2(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk) !updated prognostic tangential velocity (component2)
         REAL, INTENT(inout) :: iomomflux(ib - hi:ie + hi, jb - hj:je + hj, kb-hk:ke + hk) !a field to save the momentum flux
         REAL, INTENT(in)    :: utang1(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk)  !tangential velocity field
         REAL, INTENT(in)    :: utang2(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk)  !second tangential velocity field
         REAL, INTENT(in)    :: z0
         REAL, INTENT(in)    :: fkar2    !fkar^2, von Karman constant squared
         INTEGER, INTENT(in) :: n        ! number of the block, used to get i,j,k-indeces
         INTEGER, INTENT(in) :: ind      ! in case of y-wall (case 3x & 4x) "ind" is used for j-index, otherwise this is irrelevant
         INTEGER, INTENT(in) :: wforient !orientation of the facet see below:
         !frist digit, orientation of wall, determines iteration indices
         !second digit, if for momentum or for scalar (necessary because of staggered grid -> which variable to interpolate)
         !xlow=1,xup=2,yup=3,ylow=4,z=5
         !momentum=1

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CASES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CASES FOR MOMENTUM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SELECT CASE (wforient)
         !!!!!!!!!!!!!!!SPECIAL CASES FOR THE SURFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !can actually be made redundant and just be replaced by standard horizontal case (doesn't really matter though)
         ! SO: This code is essentially unchanged from uDALES v1, and should probably be moved out of this file in a later release.
         CASE (91) !surface momentum flux

            k = kb
            km = k - 1

            delta  = 0.5*dzf(k) !might need attention on streched grids! as well as the dzfi when updating up
            logdz2 = LOG(delta/z0)**2
            ctm    = fkar2/(logdz2)

            DO j = jb, je !u component
               jp = j + 1
               DO i = ib, ie
                  im = i - 1
                  utang1Int = utang1(i, j, k)
                  utang2Int = (utang2(i, j, k) + utang2(im, j, k) + utang2(i, jp, k) + utang2(im, jp, k))*0.25
                  utangInt = max(umin, (utang1Int**2 + utang2Int**2))
                  !dummy = (utang1Int**2)*ctm
                  dummy = abs(utang1Int)*sqrt(utangInt)*ctm
                  bcmomflux = SIGN(dummy, utang1Int)
                  iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
                  emom = (dzf(km)*(ekm(i, j, k)*dxf(im) + ekm(im, j, k)*dxf(i)) + & ! dx is non-equidistant
                          dzf(k)*(ekm(i, j, km)*dxf(im) + ekm(im, j, km)*dxf(i)))*dxhi(i)*dzhiq(k)
                  iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(i, j, km))*emom*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)
               END DO
            END DO

            DO j = jb, je !v component
               jm = j - 1
               DO i = ib, ie
                  ip = i + 1
                  utang1Int = (utang1(i, j, k) + utang1(i, jm, k) + utang1(ip, jm, k) + utang1(ip, j, k))*0.25
                  utang2Int = utang2(i, j, k)
                  utangInt = max(umin, (utang1Int**2 + utang2Int**2))
                  !dummy = (utang2Int**2)*ctm
                  dummy = abs(utang2Int)*sqrt(utangInt)*ctm
                  bcmomflux = SIGN(dummy, utang2Int)
                  iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
                  eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, jm, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, jm, km)))*dzhiq(k)
                  iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)
               END DO
            END DO

         END SELECT

   END SUBROUTINE wfmneutral
#endif

END MODULE wallfunction_neutral
