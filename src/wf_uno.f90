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
SUBROUTINE wfuno(hi,hj,hk,iout1,iout2,iot,iomomflux,iotflux,iocth,obcTfluxA,utang1,utang2,Tcell,Twall,z0,z0h,n,ind,wforient)
   !wfuno
   !calculating wall function for momentum and scalars following Cai2012&Uno1995, extension of Louis 1979 method to rough walls
   !fluxes in m2/s2 and Km/s
   USE modglobal, ONLY : dzf,dzfi,dzh2i,dzhi,dzhiq,dy,dyi,dy2i,dyi5,dxf,dxh,dxfi,dxhi,dxh2i,ib,ie,jb,je,kb,ke,fkar,grav,jmax,rk3step,kmax,jge,jgb
   USE modsubgriddata, ONLY:ekh, ekm
   USE modmpi, ONLY:myid
   USE initfac, ONLY:block
   USE modibmdata
   REAL, EXTERNAL :: unom
   INTEGER i, j, k, jl, ju, kl, ku, il, iu, km, im, jm, ip, jp, kp
   REAL :: Ribl0 = 0. !initial guess of Ribl based on Ts

   REAL :: bcTflux = 0. !temp storage for temperature flux
   REAL :: bcmomflux = 0. !temp storage for momentum flux
   REAL :: ctm = 0. !momentum transfer coefficient
   REAL :: cth = 0. !heat transfer coefficient
   REAL :: dummy = 0. !for debugging
   REAL :: delta = 0. !distance from wall
   REAL :: logdz = 0. !log(delta/z0)
   REAL :: logdzh = 0. !log(delta/z0h)
   REAL :: logzh = 0. !log(z0/z0h)
   REAL :: sqdz = 0. !sqrt(delta/z0)
   REAL :: utang1Int !Interpolated 1st tangential velocity component needed for stability calculation (to T location)
   REAL :: utang2Int !Interpolated 2nd tangential velocity component needed for stability calculation (to T location)
   REAL :: utangInt !Interpolated absolute tangential velocity
   REAL :: dT !Temperature difference between wall and cell
   REAL :: fkar2 !fkar^2, von Karman constant squared
   REAL :: emmo = 0., epmo = 0., epom = 0., emom = 0., eopm = 0., eomm = 0., empo = 0.
   REAL :: umin = 0.0001 !m^2/s^2

   INTEGER, INTENT(in) :: hi !<size of halo in i
   INTEGER, INTENT(in) :: hj !<size of halo in j
   INTEGER, INTENT(in) :: hk !<size of halo in k
   REAL, INTENT(out)   :: obcTfluxA !temperature flux of entire wall facet (double sum over indeces) [Km/s]
   REAL, INTENT(inout) :: iout1(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk) !updated prognostic tangential velocity (component1)
   REAL, INTENT(inout) :: iout2(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk) !updated prognostic tangential velocity (component2)
   REAL, INTENT(inout) :: iot(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk) !updated prognostic temperature
   REAL, INTENT(inout) :: iomomflux(ib - hi:ie + hi, jb - hj:je + hj, kb-hk:ke + hk) !a field to save the momentum flux
   REAL, INTENT(inout) :: iotflux(ib - hi:ie + hi, jb - hj:je + hj, kb-hk:ke + hk) !a field to save the heat flux
   REAL, INTENT(inout) :: iocth(ib - hi:ie + hi, jb - hj:je + hj, kb-hk:ke + hk) !heat transfer coefficient, used to calculate moisture flux
   REAL, INTENT(in)    :: Tcell(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk) !Temperature of fluid cell
   REAL, INTENT(in)    :: Twall !Temperature of surfaces !SINCE EVERY WALL HAS PRECISELY ONE TEMPERATURE (at the outside). CAREFUL IF THIS EVER CHANGES (i.e. multiple EB facets per wall)
   REAL, INTENT(in)    :: z0
   REAL, INTENT(in)    :: z0h
   REAL, INTENT(in)    :: utang1(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk) !tangential velocity field
   REAL, INTENT(in)    :: utang2(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk) !second tangential velocity field
   INTEGER, INTENT(in) :: n ! number of the block, used to get i,j,k-indeces
   INTEGER, INTENT(in) :: ind ! in case of y-wall (case 3x & 4x) "ind" is used for j-index, otherwise this is irrelevant
   INTEGER, INTENT(in) :: wforient !orientation of the facet see below:
   !frist digit, orientation of wall, determines iteration indices
   !second digit, if for momentum or for scalar (necessary because of staggered grid -> which variable to interpolate)
   !xlow=1,xup=2,yup=3,ylow=4,z=5
   !momentum=1,scalar=2
   fkar2 = fkar**2
   obcTfluxA = 0.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CASES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CASES FOR SCALARS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      logdz = LOG(delta/z0)
      logdzh = LOG(delta/z0h)
      logzh = LOG(z0/z0h)
      sqdz = SQRT(delta/z0)

      DO j = jl, ju !u component
         DO i = il, iu
            utang1Int = utang1(i, j, k)
            utang2Int = (utang2(i, j, k) + utang2(i - 1, j, k) + utang2(i, j + 1, k) + utang2(i - 1, j + 1, k))*0.25
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            dT = ((Tcell(i, j, k) + Tcell(i - 1, j, k)) - (Twall + Twall))*0.5
            Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
            ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
            !dummy = (utang1Int**2)*ctm
            dummy = abs(utang1Int)*sqrt(utangInt)*ctm
            bcmomflux = SIGN(dummy, utang1Int) !bcmomflux=u_star^2
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
            emom = (dzf(km)*(ekm(i, j, k)*dxf(i - 1) + ekm(i - 1, j, k)*dxf(i)) + & ! dx is non-equidistant
                    dzf(k)*(ekm(i, j, km)*dxf(i - 1) + ekm(i - 1, j, km)*dxf(i)))*dxhi(i)*dzhiq(k)
            iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(i, j, km))*emom*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
         END DO
      END DO

      DO j = jl, ju
         DO i = il, iu
            utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
            utang2Int = utang2(i, j, k)
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
            Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
            ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
            !dummy = (utang2Int**2)*ctm !save result and update field
            dummy = abs(utang2Int)*sqrt(utangInt)*ctm
            bcmomflux = SIGN(dummy, utang2Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
            eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
            iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
         END DO
      END DO

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   CASE (92) !surface temperature flux
      k = kb !block location
      ku = k !shear velocity location
      il = ib
      iu = ie
      jl = jb
      ju = je

      delta = dzf(k)*0.5
      logdz = LOG(delta/z0)
      logdzh = LOG(delta/z0h)
      logzh = LOG(z0/z0h)
      sqdz = SQRT(delta/z0)

      DO j = jl, ju
         DO i = il, iu
            utang1Int = (utang1(i, j, ku) + utang1(i + 1, j, ku))*0.5
            utang2Int = (utang2(i, j, ku) + utang2(i, j + 1, ku))*0.5
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            dT = (Tcell(i, j, ku) - Twall)
            Ribl0 = grav*delta*dT/(Twall*utangInt) !
            call unoh(bcTflux, cth, logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
            obcTfluxA = obcTfluxA + bcTflux
            iotflux(i, j, k) = iotflux(i, j, k) + bcTflux*dzfi(k)
            iot(i, j, ku) = iot(i, j, ku) + &
                            0.5*(dzf(k - 1)*ekh(i, j, k) + dzf(k)*ekh(i, j, k - 1))* & ! zero flux
                            (Tcell(i, j, k) - Tcell(i, j, k - 1))*dzh2i(ku)*dzfi(ku) &
                            - bcTflux*dzfi(k)
         END DO
      END DO

   END SELECT

END SUBROUTINE wfuno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!for scalar
!FUNCTION unoh(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !for heat, the bit that does not change no matter what wall
SUBROUTINE unoh(otf, octh, logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !for heat, the bit that does not change no matter what wall
   use modglobal, only : prandtlturb
!flux in Km/s
   IMPLICIT NONE
   REAL, INTENT(in) :: logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2
   REAL, INTENT(out) :: octh, otf
   REAL :: Ribl1, Fm, Fh, cm, ch, M, dTrough, cth
   REAL, PARAMETER :: b1 = 9.4 !parameters from Uno1995
   REAL, PARAMETER :: b2 = 4.7
   REAL, PARAMETER :: dm = 7.4
   REAL, PARAMETER :: dh = 5.3

   octh = 0.
   otf = 0.
   IF (Ribl0 > 0.) THEN !0.25 approx critical for bulk Richardson number  => stable
      Fm = 1./(1. + b2*Ribl0)**2 !Eq. 4
      Fh = Fm !Eq. 4
   ELSE ! => unstable
      cm = (dm*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
      ch = (dh*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
      Fm = 1. - (b1*Ribl0)/(1. + cm*SQRT(ABS(Ribl0))) !Eq. 3
      Fh = 1. - (b1*Ribl0)/(1. + ch*SQRT(ABS(Ribl0))) !Eq. 3
   END IF

   M = prandtlturb*logdz*SQRT(Fm)/Fh !Eq. 14

   Ribl1 = Ribl0 - Ribl0*prandtlturb*logzh/(prandtlturb*logzh + M) !Eq. 17

   !interate to get new Richardson number
   IF (Ribl1 > 0.) THEN !0.25 approx critical for bulk Richardson number  => stable
      Fm = 1./(1. + b2*Ribl1)**2 !Eq. 4
      Fh = Fm !Eq. 4
   ELSE ! => unstable
      cm = (dm*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
      ch = (dh*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
      Fm = 1. - (b1*Ribl1)/(1. + cm*SQRT(ABS(Ribl1))) !Eq. 3
      Fh = 1. - (b1*Ribl1)/(1. + ch*SQRT(ABS(Ribl1))) !Eq. 3
   END IF

   ! Uno (2)
   M = prandtlturb*logdz*SQRT(Fm)/Fh !Eq. 14
   dTrough = dT*1./(prandtlturb*logzh/M + 1.) !Eq. 13a
   octh = SQRT(utangInt)*fkar2/(logdz*logdz)*Fh/prandtlturb !Eq. 8
   otf = octh*dTrough !Eq. 2, Eq. 8

   ! ! Uno (8)
   ! octh = SQRT(utangInt)*fkar2/(logdz*logdzh)*Fh/prandtlturb !Eq. 8
   ! otf = octh*dT !Eq. 2, Eq. 8

END SUBROUTINE unoh

!!!!!!!!!!!!!
!for momentum
REAL FUNCTION unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !for momentum, this bit is not depended on orientation etc
   use modglobal, only : prandtlturb
!momentum flux in m2/s2
!dT,utang and logdzh are unused and could be removed
   IMPLICIT NONE
   REAL, INTENT(in) :: logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2
   REAL :: Ribl1, Fm, Fh, cm, ch, Ctm, M
   REAL, PARAMETER :: b1 = 9.4 !parameters from Uno1995
   REAL, PARAMETER :: b2 = 4.7
   REAL, PARAMETER :: dm = 7.4
   REAL, PARAMETER :: dh = 5.3
   !REAL, PARAMETER :: prandtlmol = 0.71

   IF (Ribl0 > 0.) THEN !0.25 approx critical for bulk Richardson number  => stable
      Fm = 1./(1. + b2*Ribl0)**2 !Eq. 4
      Fh = Fm !Eq. 4
   ELSE ! => unstable
      cm = (dm*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
      ch = (dh*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
      Fm = 1. - (b1*Ribl0)/(1. + cm*SQRT(ABS(Ribl0))) !Eq. 3
      Fh = 1. - (b1*Ribl0)/(1. + ch*SQRT(ABS(Ribl0))) !Eq. 3
   END IF

   M = prandtlturb*logdz*SQRT(Fm)/Fh !Eq. 14

   Ribl1 = Ribl0 - Ribl0*prandtlturb*logzh/(prandtlturb*logzh + M) !Eq. 17

   !interate to get new Richardson number
   IF (Ribl1 > 0.) THEN !0.25 approx critical for bulk Richardson number  => stable
      Fm = 1./(1. + b2*Ribl1)**2 !Eq. 4
   ELSE ! => unstable
      cm = (dm*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
      Fm = 1. - (b1*Ribl1)/(1. + cm*SQRT(ABS(Ribl1))) !Eq. 3
   END IF

   Ctm = fkar2/(logdz**2)*Fm !Eq. 7
   unom = Ctm !Eq. 2, Eq. 8
END FUNCTION unom
