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
SUBROUTINE wfGR(hi,hj,hk,ioq,ioqflux,icth,obcqfluxA,qcell,qwall,hurel,resc,ress,n,ind,wforient)
   !wfGR
   USE modglobal, ONLY : dzf,dzfi,dzh2i,dzhi,dzhiq,dzfi5,dy,dyi,dy2i,dyi5,dxf,dxh,dxfi,dxhi,dxh2i,dxfi5,ib,ie,jb,je,kb,ke,fkar,grav,jmax,rk3step
   USE modsubgriddata, ONLY:ekh
   USE modmpi, ONLY:myid
   USE initfac, ONLY:block, faclGR
   USE modibmdata
   use modfields, only : u0, v0, w0
   INTEGER i, j, k, jl, ju, kl, ku, il, iu, km, im, jm, ip, jp, kp

   REAL :: bcqflux = 0. !temp storage for temperature flux
   REAL :: bcmomflux = 0. !temp storage for momentum flux
   REAL :: dummy = 0. !for debugging
   REAL :: delta = 0. !distance from wall
   REAL :: fkar2 !fkar^2, von Karman constant squared
   REAL :: emmo = 0., epmo = 0., epom = 0., emom = 0., eopm = 0., eomm = 0., empo = 0.
   REAL :: umin = 0.0001 !m^2/s^2
   REAL :: cveg=0.8 !hardcoded for now, !fraction of GR covered in vegetation, should be made into a proper model parameter (-> modglobal)

   INTEGER, INTENT(in) :: hi !<size of halo in i
   INTEGER, INTENT(in) :: hj !<size of halo in j
   INTEGER, INTENT(in) :: hk !<size of halo in k
   REAL, INTENT(out)   :: obcqfluxA; !temperature flux of entire wall facet (double sum over indeces) [Km/s]
   REAL, INTENT(inout) :: ioq(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk) !updated prognostic temperature
   REAL, INTENT(inout) :: ioqflux(ib - hi:ie + hi, jb - hj:je + hj, kb-hk:ke + hk) !a field to save the momentum flux
   REAL, INTENT(in)   :: icth(ib - hi:ie + hi, jb - hj:je + hj, kb-hk:ke + hk) ; !heat transfer coefficient, used to calculate moisture flux
   REAL, INTENT(in)    :: qcell(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk) !Temperature of fluid cell
   !       real, intent(in)    :: Troof(ib:ie,jb:je,kb:ke)  !Temperature of horizontal surfaces (also includes roads)
   REAL, INTENT(in)    :: qwall
   REAL, INTENT(in)    :: hurel
   REAL, INTENT(in)    :: resc
   REAL, INTENT(in)    :: ress
   INTEGER, INTENT(in) :: n ! number of the block, used to get i,j,k-indeces
   INTEGER, INTENT(in) :: ind ! in case of y-wall (case 3x & 4x) "ind" is used for j-index, otherwise this is irrelevant
   INTEGER, INTENT(in) :: wforient !frist digit, orientation of wall, determines iteration idices and if Twall or Troof is used
   fkar2 = fkar**2
   obcqfluxA = 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CASES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CASES FOR SCALARS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SELECT CASE (wforient)
   CASE (12) !wall in yz -> wf in x (=vertical), lower wall, west wall
      i = block(n, 1) - 1 !wall property and fluid index
      ip = i + 1 !index to remove subgrid flux
      jl = MAX(block(n, 3) - myid*jmax, 1) ! starting j-index
      ju = MIN(block(n, 4) - myid*jmax, jmax) ! ending j-index
      kl = block(n, 5) ! starting k-index
      ku = block(n, 6) ! ending k-index

      if (faclGR(block(n, 8))) then
         DO k = kl, ku
            DO j = jl, ju
               bcqflux = min(0.,cveg*(qcell(i,j,k) - qwall)*1/(1/icth(i,j,k)+resc)+(1-cveg)*(qcell(i,j,k) - qwall*hurel)*1/(1/icth(i,j,k)+ress))
               obcqfluxA = obcqfluxA + bcqflux
               ioqflux(i, j, k) = ioqflux(i, j, k) + bcqflux*dxfi(i)
               ioq(i,j,k) = ioq(i,j,k) - bcqflux*dxfi(i) &
                          - 0.5*(ekh(ip,j,k)*dxf(i) + ekh(i,j,k)*dxf(ip))*(qcell(ip,j,k) - qcell(i,j,k))*dxh2i(ip)*dxfi(i) &
                          + (u0(ip, j, k)*(qcell(ip, j, k)*dxf(i) + qcell(i, j, k)*dxf(ip))*dxhi(ip))*dxfi5(i) &
                          - (u0(ip, j, k)*(qcell(i , j, k)*dxf(i) + qcell(i, j, k)*dxf(ip))*dxhi(ip))*dxfi5(i)
            END DO
         END DO

      else
         DO k = kl, ku
            DO j = jl, ju
               ioq(i,j,k) = ioq(i,j,k) &
                          - 0.5*(ekh(ip,j,k)*dxf(i) + ekh(i,j,k)*dxf(ip))*(qcell(ip,j,k) - qcell(i,j,k))*dxh2i(ip)*dxfi(i) &
                          + (u0(ip, j, k)*(qcell(ip, j, k)*dxf(i) + qcell(i, j, k)*dxf(ip))*dxhi(ip))*dxfi5(i) &
                          - (u0(ip, j, k)*(qcell(i , j, k)*dxf(i) + qcell(i, j, k)*dxf(ip))*dxhi(ip))*dxfi5(i)
            END DO
         END DO
      end if

!!! case 22 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   CASE (22)
      i = block(n, 2) + 1 !
      im = i - 1 !
      jl = MAX(block(n, 3) - myid*jmax, 1) ! starting j-index
      ju = MIN(block(n, 4) - myid*jmax, jmax) ! ending j-index
      kl = block(n, 5) ! starting k-index
      ku = block(n, 6) ! ending k-index

      if (faclGR(block(n, 9))) then
         DO k = kl, ku
            DO j = jl, ju
               bcqflux = min(0.,cveg*(qcell(i,j,k) - qwall)*1/(1/icth(i,j,k)+resc)+(1-cveg)*(qcell(i,j,k) - qwall*hurel)*1/(1/icth(i,j,k)+ress))
               obcqfluxA = obcqfluxA + bcqflux
               ioqflux(i, j, k) = ioqflux(i, j, k) + bcqflux*dxfi(i)
               ioq(i,j,k) = ioq(i,j,k) - bcqflux*dxfi(i) &
                           + 0.5*(ekh(i,j,k)*dxf(im) + ekh(im,j,k)*dxf(i))*(qcell(i,j,k) - qcell(im,j,k))*dxh2i(i)*dxfi(i) &
                           + (- u0(i, j, k)*(qcell(im, j, k)*dxf(i) + qcell(i, j, k)*dxf(im))*dxhi(i))*dxfi5(i) &
                           - (- u0(i, j, k)*(qcell(i , j, k)*dxf(i) + qcell(i, j, k)*dxf(im))*dxhi(i))*dxfi5(i)
            END DO
         END DO

      else
         DO k = kl, ku
            DO j = jl, ju
              ioq(i,j,k) = ioq(i,j,k) &
                         + 0.5*(ekh(i,j,k)*dxf(im) + ekh(im,j,k)*dxf(i))*(qcell(i,j,k) - qcell(im,j,k))*dxh2i(i)*dxfi(i) &
                         + (- u0(i, j, k)*(qcell(im, j, k)*dxf(i) + qcell(i, j, k)*dxf(im))*dxhi(i))*dxfi5(i) &
                         - (- u0(i, j, k)*(qcell(i , j, k)*dxf(i) + qcell(i, j, k)*dxf(im))*dxhi(i))*dxfi5(i)
            END DO
         END DO
      end if

!!! case 32 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   CASE (32) !wall in xz -> wf in y (=vertical) upper, north wall
      j = ind
      jm = j - 1
      il = block(n, 1)
      iu = block(n, 2)
      kl = block(n, 5)
      ku = block(n, 6)

      if (faclGR(block(n, 10))) then
         DO k = kl, ku
            DO i = il, iu
               bcqflux = min(0., cveg*(qcell(i,j,k) - qwall)*1/(1/icth(i,j,k)+resc)+(1-cveg)*(qcell(i,j,k) - qwall*hurel)*1/(1/icth(i,j,k)+ress))
               obcqfluxA = obcqfluxA + bcqflux
               ioqflux(i, j, k) = ioqflux(i, j, k) + bcqflux*dyi
               ioq(i, j, k) = ioq(i, j, k) -bcqflux*dyi &
                            + 0.5*(ekh(i, j, k) + ekh(i, jm, k))*(qcell(i, j, k) - qcell(i, jm, k))*dy2i &
                            + (- v0(i, j, k)*(qcell(i, jm, k) + qcell(i, j, k)))*dyi5 &
                            - (- v0(i, j, k)*(qcell(i, j , k) + qcell(i, j, k)))*dyi5

            END DO
         END DO

      else
         DO k = kl, ku
            DO i = il, iu
               ioq(i, j, k) = ioq(i, j, k) &
                            + 0.5*(ekh(i, j, k) + ekh(i, jm, k))*(qcell(i, j, k) - qcell(i, jm, k))*dy2i &
                            + (- v0(i, j, k)*(qcell(i, jm, k) + qcell(i, j, k)))*dyi5 &
                            - (- v0(i, j, k)*(qcell(i, j , k) + qcell(i, j, k)))*dyi5
            END DO
         END DO
      end if

!!! case 42 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   CASE (42) !wall in xz -> wf in y (=vertical) lower, south wall
      j = ind
      jp = j + 1
      il = block(n, 1)
      iu = block(n, 2)
      kl = block(n, 5)
      ku = block(n, 6)

      if (faclGR(block(n, 11))) then
         DO k = kl, ku
            DO i = il, iu
               bcqflux = min(0., cveg*(qcell(i,j,k) - qwall)*1/(1/icth(i,j,k)+resc)+(1-cveg)*(qcell(i,j,k) - qwall*hurel)*1/(1/icth(i,j,k)+ress))
               obcqfluxA = obcqfluxA + bcqflux
               ioqflux(i, j, k) = ioqflux(i, j, k) + bcqflux*dyi
               ioq(i, j, k) = ioq(i, j, k) - bcqflux*dyi &
                            - 0.5*(ekh(i, jp, k) + ekh(i, j, k))*(qcell(i, jp, k) - qcell(i, j, k))*dy2i &
                            + v0(i, jp, k)*(qcell(i, jp, k) + qcell(i, j, k))*dyi5 &
                            - v0(i, jp, k)*(qcell(i, j , k) + qcell(i, j, k))*dyi5

            END DO
         END DO

      else
         DO k = kl, ku
            DO i = il, iu
               ioq(i, j, k) = ioq(i, j, k) &
                            - 0.5*(ekh(i, jp, k) + ekh(i, j, k))*(qcell(i, jp, k) - qcell(i, j, k))*dy2i &
                            + v0(i, jp, k)*(qcell(i, jp, k) + qcell(i, j, k))*dyi5 &
                            - v0(i, jp, k)*(qcell(i, j , k) + qcell(i, j, k))*dyi5
            END DO
         END DO
      end if

!!! case 52 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   CASE (52) !wall in xy -> wf in z (=horizontal), top wall
      k = block(n, 6) + 1 !block location
      !if (.not.(k.gt.kmax)) then
      if (k > ke) return
      km = k - 1 !
      il = block(n, 1)
      iu = block(n, 2)
      jl = MAX(block(n, 3) - myid*jmax, 1)
      ju = MIN(block(n, 4) - myid*jmax, jmax)

      if (faclGR(block(n, 7))) then
         DO j = jl, ju
            DO i = il, iu
               bcqflux = min(0., cveg*(qcell(i,j,k) - qwall)*1/(1/icth(i,j,k)+resc)+(1-cveg)*(qcell(i,j,k) - qwall*hurel)*1/(1/icth(i,j,k)+ress))
               obcqfluxA = obcqfluxA + bcqflux
               ioqflux(i, j, k) = ioqflux(i, j, k) + bcqflux*dzfi(k)
               ioq(i, j, k) = ioq(i, j, k) - bcqflux*dzfi(k) &
                            + 0.5*(dzf(km)*ekh(i, j, k) + dzf(k)*ekh(i, j, km))*(qcell(i, j, k) - qcell(i, j, km))*dzh2i(k)*dzfi(k) &
                            + (- w0(i, j, k)*(qcell(i, j, km)*dzf(k) + qcell(i, j, k)*dzf(km))*dzhi(k))*dzfi5(k) &
                            - (- w0(i, j, k)*(qcell(i, j, k )*dzf(k) + qcell(i, j, k)*dzf(km))*dzhi(k))*dzfi5(k)

            END DO
         END DO

      else
         DO j = jl, ju
            DO i = il, iu
               ioq(i, j, k) = ioq(i, j, k) &
                              + 0.5*(dzf(km)*ekh(i, j, k) + dzf(k)*ekh(i, j, km))*(qcell(i, j, k) - qcell(i, j, km))*dzh2i(k)*dzfi(k) &
                              + (- w0(i, j, k)*(qcell(i, j, km)*dzf(k) + qcell(i, j, k)*dzf(km))*dzhi(k))*dzfi5(k) &
                              - (- w0(i, j, k)*(qcell(i, j, k )*dzf(k) + qcell(i, j, k)*dzf(km))*dzhi(k))*dzfi5(k)
            END DO
         END DO
      end if
      !end if

      ! CASE (62) !wall in xy -> wf in z (=horizontal), bottom wall
      !    k = block(n, 5) - 1 !block location
      !    !if (.not.(k.lt.0)) then
      !    if (k < kb) return
      !    kp = k + 1 !
      !    il = block(n, 1)
      !    iu = block(n, 2)
      !    jl = MAX(block(n, 3) - myid*jmax, 1)
      !    ju = MIN(block(n, 4) - myid*jmax, jmax)
      !
      !    if (faclGR(block(n, 12))) then
      !       DO j = jl, ju
      !          DO i = il, iu
      !             bcqflux = min(0., cveg*(qcell(i,j,k) - qwall)*1/(1/icth(i,j,k)+resc)+(1-cveg)*(qcell(i,j,k) - qwall*hurel)*1/(1/icth(i,j,k)+ress))
      !             obcqfluxA = obcqfluxA + bcqflux
      !             ioqflux(i, j, k) = ioqflux(i, j, k) + bcqflux*dzfi(k)
      !             ioq(i, j, k) = ioq(i, j, k) - bcqflux*dzfi(k) &
      !                          - 0.5*(dzf(kp)*ekh(i, j, k) + dzf(k)*ekh(i, j, kp))*(qcell(i, j, kp) - qcell(i, j, k))*dzh2i(k)*dzfi(k) &
      !                          + w0(i, j, kp)*(qcell(i, j, kp)*dzf(k) + qcell(i, j, k)*dzf(kp))*dzhi(kp)*dzfi5(k) &
      !                          - w0(i, j, kp)*(qcell(i, j, k )*dzf(k) + qcell(i, j, k)*dzf(kp))*dzhi(kp)*dzfi5(k)
      !
      !          END DO
      !       END DO
      !
      !    else
      !       DO j = jl, ju
      !          DO i = il, iu
      !             ioq(i, j, k) = ioq(i, j, k) &
      !                          - 0.5*(dzf(kp)*ekh(i, j, k) + dzf(k)*ekh(i, j, kp))*(qcell(i, j, kp) - qcell(i, j, k))*dzh2i(k)*dzfi(k) &
      !                          + w0(i, j, kp)*(qcell(i, j, kp)*dzf(k) + qcell(i, j, k)*dzf(kp))*dzhi(kp)*dzfi5(k) &
      !                          - w0(i, j, kp)*(qcell(i, j, k )*dzf(k) + qcell(i, j, k)*dzf(kp))*dzhi(kp)*dzfi5(k)
      !          END DO
      !       END DO
      !    end if
      !    !end if

END SELECT

END SUBROUTINE wfGR
