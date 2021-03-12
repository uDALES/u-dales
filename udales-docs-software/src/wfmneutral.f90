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
   REAL :: fkar2 = fkar**2 !fkar^2, von Karman constant squared
   REAL :: emmo = 0., epmo = 0., epom = 0., emom = 0., eopm = 0., eomm = 0., empo = 0.

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

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CASES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CASES FOR MOMENTUM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SELECT CASE (wforient)
   CASE (11) !wfmneutral11, case 11 , west wall
      i = block(n, 1) - 1 !fluid location (also where wall variables are stored)
      ip = i + 1 !inside wall, used for subtracting original diffusion term
      jl = MAX(block(n, 3) - myid*jmax, 1) + 1 ! starting j-index      !might cause problem when jl=1
      ju = MIN(block(n, 4) - myid*jmax, jmax) ! ending j-index     !might cause problem when ju=jmax
      kl = block(n, 5) ! starting k-index
      ku = block(n, 6) ! ending k-index

      delta = dxf(i)*0.5
      logdz2 = LOG(delta/z0)**2
            
      !v west
      DO k = kl, ku
         DO j = jl, ju

            utang1Int = utang1(i, j, k)

            ctm = fkar2/(logdz2)
            dummy = (utang1Int**2)*ctm
            bcmomflux = SIGN(dummy, utang1Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)

            epmo = 0.25*((ekm(i, j, k) + ekm(i, j - 1, k))*dxf(ip) + &
                         (ekm(ip, j, k) + ekm(ip, j - 1, k))*dxf(i))*dxhi(ip)

            iout1(i, j, k) = iout1(i, j, k) - (utang1(ip, j, k) - utang1(i, j, k))*epmo*dxhi(ip)*dxfi(i) - bcmomflux*dxfi(i) !

         END DO
      END DO

      !v west edge south
      j = MAX(block(n, 3) - myid*jmax, 1)
      DO k = kl, ku
         utang1Int = utang1(i, j, k)
         dummy = (utang1Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)*0.5 !only half since on edge of block (another half might come from another processor?)

         !   epmo = 0.5 * (ekm(ip,j,k)*dxf(i) + ekm(i,j,k)*dxf(ip)) * dxhi(ip)
         epmo = 0.25*((ekm(i, j, k) + ekm(i, j - 1, k))*dxf(ip) + &
                      (ekm(ip, j, k) + ekm(ip, j - 1, k))*dxf(i))*dxhi(ip)

         iout1(i, j, k) = iout1(i, j, k) - ((utang1(ip, j, k) - utang1(i, j, k))*epmo*dxhi(ip)*dxfi(i) - bcmomflux*dxfi(i))*0.5 ! remove standard diffusion apply only half of wall-flux since it's an edge
         !only half of the flux, since only half of the control-volume around v is touching this facet (other half is either in free air or touching another facet)
      END DO

      !v west edge north
      j = MIN(block(n, 4) - myid*jmax, jmax) + 1
      DO k = kl, ku

         utang1Int = utang1(i, j, k)
         dummy = (utang1Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)*0.5

         epmo = 0.25*((ekm(i, j, k) + ekm(i, j - 1, k))*dxf(ip) + &
                      (ekm(ip, j, k) + ekm(ip, j - 1, k))*dxf(i))*dxhi(ip)

         iout1(i, j, k) = iout1(i, j, k) - ((utang1(ip, j, k) - utang1(i, j, k))*epmo*dxhi(ip)*dxfi(i) - bcmomflux*dxfi(i))*0.5 ! %remove standard diffusion apply only half of wall-flux since it's an edge

      END DO

      !w west

      jl = MAX(block(n, 3) - myid*jmax, 1) !
      kl = block(n, 5) + 1 !
      DO k = kl, ku
         DO j = jl, ju

            utang2Int = utang2(i, j, k)
            ctm = fkar2/(logdz2)
            dummy = (utang2Int**2)*ctm
            bcmomflux = SIGN(dummy, utang2Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)

            epom = (dzf(k - 1)*(ekm(i, j, k)*dxf(ip) + ekm(ip, j, k)*dxf(i))*dxhi(ip) + &
                    dzf(k)*(ekm(i, j, k - 1)*dxf(ip) + ekm(ip, j, k - 1)*dxf(i))*dxhi(ip))*dzhiq(k)

            iout2(i, j, k) = iout2(i, j, k) - (utang2(ip, j, k) - utang2(i, j, k))*epom*dxhi(ip)*dxfi(i) - bcmomflux*dxfi(i) !
         END DO
      END DO

      !w west top edge
      k = block(n, 6) + 1 ! ending k-index
      km = k - 1
      DO j = jl, ju
         utang2Int = utang2(i, j, k)
         dummy = (utang2Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)*0.5

         epom = (dzf(km)*(ekm(i, j, k)*dxf(ip) + ekm(ip, j, k)*dxf(i))*dxhi(ip) + &
                 dzf(k)*(ekm(i, j, km)*dxf(ip) + ekm(ip, j, km)*dxf(i))*dxhi(ip))*dzhiq(k)

         iout2(i, j, k) = iout2(i, j, k) - ((utang2(ip, j, k) - utang2(i, j, k))*epom*dxhi(ip)*dxfi(i) - bcmomflux*dxfi(i))*0.5 !
      END DO

      !w west bottom edge
      k = block(n, 6)  ! ending k-index
      if (k.gt.0) then
      km = k - 1
      DO j = jl, ju
         utang2Int = utang2(i, j, k)
         dummy = (utang2Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)*0.5

         epom = (dzf(km)*(ekm(i, j, k)*dxf(ip) + ekm(ip, j, k)*dxf(i))*dxhi(ip) + &
                 dzf(k)*(ekm(i, j, km)*dxf(ip) + ekm(ip, j, km)*dxf(i))*dxhi(ip))*dzhiq(k)

         iout2(i, j, k) = iout2(i, j, k) - ((utang2(ip, j, k) - utang2(i, j, k))*epom*dxhi(ip)*dxfi(i) - bcmomflux*dxfi(i))*0.5 !
      END DO
      end if

      

!!! case 21 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !wfmneutral 21 !wall in yz -> wf in x (=vertical), upper wall, east wall
   CASE (21)
      !v east
      i = block(n, 2) + 1 !fluid
      im = i - 1 !inside block
      jl = MAX(block(n, 3) - myid*jmax, 1) + 1 ! starting j-index      !might cause problem when jl=1
      ju = MIN(block(n, 4) - myid*jmax, jmax) ! ending j-index     !might cause problem when ju=jmax
      kl = block(n, 5) ! starting k-index
      ku = block(n, 6) ! ending k-index

      delta = dxh(i)*0.5
      logdz2 = LOG(delta/z0)**2

      DO k = kl, ku
         DO j = jl, ju

            utang1Int = utang1(i, j, k)
            ctm = fkar2/(logdz2)
            dummy = (utang1Int**2)*ctm
            bcmomflux = SIGN(dummy, utang1Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)

            emmo = 0.25*((ekm(i, j, k) + ekm(i, j - 1, k))*dxf(im) + (ekm(im, j - 1, k) + ekm(im, j, k))*dxf(i))*dxhi(i) ! dx is non-equidistant

            iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(im, j, k))*emmo*dxhi(i)*dxfi(i) - bcmomflux*dxfi(i) !
         END DO
      END DO

      !v east edge south
      j = MAX(block(n, 3) - myid*jmax, 1)
      DO k = kl, ku
         utang1Int = utang1(i, j, k)
         dummy = (utang1Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)*0.5

         emmo = 0.25*((ekm(i, j, k) + ekm(i, j - 1, k))*dxf(im) + (ekm(im, j - 1, k) + ekm(im, j, k))*dxf(i))*dxhi(i) ! dx is non-equidistant

         iout1(i, j, k) = iout1(i, j, k) + ((utang1(i, j, k) - utang1(im, j, k))*emmo*dxhi(i)*dxfi(i) - bcmomflux*dxfi(i))*0.5 !

      END DO

      !v east edge north
      j = MIN(block(n, 4) - myid*jmax, jmax) + 1 !
      DO k = kl, ku
         utang1Int = utang1(i, j, k)
         dummy = (utang1Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)*0.5

         emmo = 0.25*((ekm(i, j, k) + ekm(i, j - 1, k))*dxf(im) + (ekm(im, j - 1, k) + ekm(im, j, k))*dxf(i))*dxhi(i) ! dx is non-equidistant

         iout1(i, j, k) = iout1(i, j, k) + ((utang1(i, j, k) - utang1(im, j, k))*emmo*dxhi(i)*dxfi(i) - bcmomflux*dxfi(i))*0.5 !

      END DO

      !w east
      jl = MAX(block(n, 3) - myid*jmax, 1) !
      kl = block(n, 5) + 1 !
      DO k = kl, ku
         DO j = jl, ju

            utang2Int = utang2(i, j, k)
            ctm = fkar2/(logdz2)
            dummy = (utang2Int**2)*ctm
            bcmomflux = SIGN(dummy, utang2Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)

            emom = (dzf(k - 1)*(ekm(i, j, k)*dxf(im) + ekm(im, j, k)*dxf(i))*dxhi(i) + &
                    dzf(k)*(ekm(i, j, k - 1)*dxf(im) + ekm(im, j, k - 1)*dxf(i))*dxhi(i))*dzhiq(k)

            iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(im, j, k))*emom*dxhi(i)*dxfi(i) - bcmomflux*dxfi(i) !
         END DO
      END DO
      !w east edge top
      k = block(n, 6) + 1 ! ending k-index
      DO j = jl, ju
         utang2Int = utang2(i, j, k)
         dummy = (utang2Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)*0.5

         emom = (dzf(k - 1)*(ekm(i, j, k)*dxf(im) + ekm(im, j, k)*dxf(i))*dxhi(i) + &
                 dzf(k)*(ekm(i, j, k - 1)*dxf(im) + ekm(im, j, k - 1)*dxf(i))*dxhi(i))*dzhiq(k)

         iout2(i, j, k) = iout2(i, j, k) + ((utang2(i, j, k) - utang2(im, j, k))*emom*dxhi(i)*dxfi(i) - bcmomflux*dxfi(i))*0.5 !

      END DO

     !w east edge bot
      k = block(n, 6)  ! 
      if (k.gt.0) then
      DO j = jl, ju
         utang2Int = utang2(i, j, k)
         dummy = (utang2Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)*0.5

         emom = (dzf(k - 1)*(ekm(i, j, k)*dxf(im) + ekm(im, j, k)*dxf(i))*dxhi(i) + &
                 dzf(k)*(ekm(i, j, k - 1)*dxf(im) + ekm(im, j, k - 1)*dxf(i))*dxhi(i))*dzhiq(k)

         iout2(i, j, k) = iout2(i, j, k) + ((utang2(i, j, k) - utang2(im, j, k))*emom*dxhi(i)*dxfi(i) - bcmomflux*dxfi(i))*0.5 !

      END DO
      end if

!!! case 31 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !wfuno 31
   CASE (31) !wall in xz -> wf in y (=vertical) upper, north wall

      j = ind
      jm = j - 1

      il = block(n, 1) + 1
      iu = block(n, 2)
      kl = block(n, 5)
      ku = block(n, 6)

      delta = 0.5*dy
      logdz2 = LOG(delta/z0)**2

      !u north
      DO k = kl, ku
         DO i = il, iu
            utang1Int = utang1(i, j, k)
            ctm = fkar2/(logdz2)
            dummy = (utang1Int**2)*ctm
            bcmomflux = SIGN(dummy, utang1Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi

            emmo = 0.25*((ekm(i, j, k) + ekm(i, jm, k))*dxf(i - 1) + (ekm(i - 1, jm, k) + ekm(i - 1, j, k))*dxf(i))*dxhi(i)

            iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(i, jm, k))*emmo*dy2i-bcmomflux*dyi !
         END DO
      END DO
      !u north east edge
      i = block(n, 2) + 1
      DO k = kl, ku
         utang1Int = utang1(i, j, k)
         dummy = (utang1Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi5
         emmo = 0.25*((ekm(i, j, k) + ekm(i, jm, k))*dxf(i - 1) + (ekm(i - 1, jm, k) + ekm(i - 1, j, k))*dxf(i))*dxhi(i)

         iout1(i, j, k) = iout1(i, j, k) + ((utang1(i, j, k) - utang1(i, jm, k))*emmo*dy2i-bcmomflux*dyi)*0.5 !

      END DO

      !u north west edge
      i = block(n, 1)
      DO k = kl, ku
         utang1Int = utang1(i, j, k)
         dummy = (utang1Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi5
         emmo = 0.25*((ekm(i, j, k) + ekm(i, jm, k))*dxf(i - 1) + (ekm(i - 1, jm, k) + ekm(i - 1, j, k))*dxf(i))*dxhi(i)
         iout1(i, j, k) = iout1(i, j, k) + ((utang1(i, j, k) - utang1(i, jm, k))*emmo*dy2i-bcmomflux*dyi)*0.5 !
      END DO

      !w north

      il = block(n, 1) !
      kl = block(n, 5) + 1 !
      DO k = kl, ku
         DO i = il, iu

            utang2Int = utang2(i, j, k)
            ctm = fkar2/(logdz2)
            dummy = (utang2Int**2)*ctm
            bcmomflux = SIGN(dummy, utang2Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi
            eomm = (dzf(k - 1)*(ekm(i, j, k) + ekm(i, jm, k)) + dzf(k)*(ekm(i, j, k - 1) + ekm(i, jm, k - 1)))*dzhiq(k) ! dz is non-eqidistant
            iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, jm, k))*eomm*dy2i-bcmomflux*dyi !
         END DO
      END DO

!w north edge top
      k = block(n, 6) + 1
      DO i = il, iu
         utang2Int = utang2(i, j, k)
         dummy = (utang2Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi5

         eomm = (dzf(k - 1)*(ekm(i, j, k) + ekm(i, jm, k)) + dzf(k)*(ekm(i, j, k - 1) + ekm(i, jm, k - 1)))*dzhiq(k) ! dz is non-eqidistant

         iout2(i, j, k) = iout2(i, j, k) + ((utang2(i, j, k) - utang2(i, jm, k))*eomm*dy2i-bcmomflux*dyi)*0.5 !

      END DO

!w north edge bot
      k = block(n, 6) 
     if (k.gt.0) then
      DO i = il, iu
         utang2Int = utang2(i, j, k)
         dummy = (utang2Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi5

         eomm = (dzf(k - 1)*(ekm(i, j, k) + ekm(i, jm, k)) + dzf(k)*(ekm(i, j, k - 1) + ekm(i, jm, k - 1)))*dzhiq(k) ! dz is non-eqidistant

         iout2(i, j, k) = iout2(i, j, k) + ((utang2(i, j, k) - utang2(i, jm, k))*eomm*dy2i-bcmomflux*dyi)*0.5 !

      END DO
end if


!!! case 41 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!`
      !wfuno41
   CASE (41) !wall in xz -> wf in y (=vertical) lower, south wall

      j = ind
      jp = j + 1
      il = block(n, 1) + 1
      iu = block(n, 2)
      kl = block(n, 5)
      ku = block(n, 6)

      delta = 0.5*dy
      logdz2 = LOG(delta/z0)**2
      
DO k = kl, ku
         DO i = il, iu

            utang1Int = utang1(i, j, k)
            ctm = fkar2/(logdz2)
            dummy = (utang1Int**2)*ctm
            bcmomflux = SIGN(dummy, utang1Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi

            empo = 0.25*((ekm(i, j, k) + ekm(i, jp, k))*dxf(i - 1) + &
                         (ekm(i - 1, j, k) + ekm(i - 1, jp, k))*dxf(i))*dxhi(i) ! dx is non-equidistant

            iout1(i, j, k) = iout1(i, j, k) - (utang1(i, jp, k) - utang1(i, j, k))*empo*dy2i-bcmomflux*dyi !

         END DO
      END DO

!u south edge west
      i = block(n, 1)
      DO k = kl, ku
         utang1Int = utang1(i, j, k)
         dummy = (utang1Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi5

         empo = 0.25*((ekm(i, j, k) + ekm(i, jp, k))*dxf(i - 1) + &
                      (ekm(i - 1, j, k) + ekm(i - 1, jp, k))*dxf(i))*dxhi(i) ! dx is non-equidistant

         iout1(i, j, k) = iout1(i, j, k) - ((utang1(i, jp, k) - utang1(i, j, k))*empo*dy2i-bcmomflux*dyi)*0.5 !

      END DO
!u south edge east
      i = block(n, 2) + 1
      DO k = kl, ku
         utang1Int = utang1(i, j, k)
         dummy = (utang1Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi5

         empo = 0.25*((ekm(i, j, k) + ekm(i, jp, k))*dxf(i - 1) + &
                      (ekm(i - 1, j, k) + ekm(i - 1, jp, k))*dxf(i))*dxhi(i) ! dx is non-equidistant

         iout1(i, j, k) = iout1(i, j, k) - ((utang1(i, jp, k) - utang1(i, j, k))*empo*dy2i-bcmomflux*dyi)*0.5 !

      END DO

!w south
      il = block(n, 1) !
      kl = block(n, 5) + 1 !
      DO k = kl, ku
         DO i = il, iu
            utang2Int = utang2(i, j, k)
            ctm = fkar2/(logdz2)
            dummy = (utang2Int**2)*ctm
            bcmomflux = SIGN(dummy, utang2Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi

            eopm = (dzf(k - 1)*(ekm(i, j, k) + ekm(i, jp, k)) + &
                    dzf(k)*(ekm(i, j, k - 1) + ekm(i, jp, k - 1)))*dzhiq(k)

            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, jp, k) - utang2(i, j, k))*eopm*dy2i-bcmomflux*dyi !
         END DO
      END DO
!w south edge top
      k = block(n, 6) + 1
      DO i = il, iu
         utang2Int = utang2(i, j, k)
         dummy = (utang2Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi5

         eopm = (dzf(k - 1)*(ekm(i, j, k) + ekm(i, jp, k)) + &
                 dzf(k)*(ekm(i, j, k - 1) + ekm(i, jp, k - 1)))*dzhiq(k)

         iout2(i, j, k) = iout2(i, j, k) - ((utang2(i, jp, k) - utang2(i, j, k))*eopm*dy2i-bcmomflux*dyi)*0.5 !

      END DO

!w south edge bot
      k = block(n, 6)
      if (k.gt.0) then 
      DO i = il, iu
         utang2Int = utang2(i, j, k)
         !call function repeatedly
         dummy = (utang2Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi5

         eopm = (dzf(k - 1)*(ekm(i, j, k) + ekm(i, jp, k)) + &
                 dzf(k)*(ekm(i, j, k - 1) + ekm(i, jp, k - 1)))*dzhiq(k)

         iout2(i, j, k) = iout2(i, j, k) - ((utang2(i, jp, k) - utang2(i, j, k))*eopm*dy2i-bcmomflux*dyi)*0.5 !

      END DO
      end if

!!!! case 51 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !wfuno51
   CASE (51) !wall in xy -> wf in z (=horizontal), top wall
      k = block(n, 6) + 1 !block location

      if (.not.(k.gt.kmax)) then
      km = k - 1 !shear velocity location
      il = block(n, 1) + 1
      iu = block(n, 2)
      jl = MAX(block(n, 3) - myid*jmax, 1)
      ju = MIN(block(n, 4) - myid*jmax, jmax)
 
      delta = 0.5*dzf(k)
      logdz2 = LOG(delta/z0)**2

      DO j = jl, ju
         DO i = il, iu
            utang1Int = utang1(i, j, k)
            ctm = fkar2/(logdz2)
            dummy = (utang1Int**2)*ctm
            bcmomflux = SIGN(dummy, utang1Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
            emom = (dzf(km)*(ekm(i, j, k)*dxf(i - 1) + ekm(i - 1, j, k)*dxf(i)) + &
                    dzf(k)*(ekm(i, j, km)*dxf(i - 1) + ekm(i - 1, j, km)*dxf(i)))*dxhi(i)*dzhiq(k)
            iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(i, j, km))*emom*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
         END DO
      END DO

!u top edge west
      i = block(n, 1)
      DO j = jl, ju
         utang1Int = utang1(i, j, k)
         dummy = (utang1Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)*0.5
         emom = (dzf(km)*(ekm(i, j, k)*dxf(i - 1) + ekm(i - 1, j, k)*dxf(i)) + &
                 dzf(k)*(ekm(i, j, km)*dxf(i - 1) + ekm(i - 1, j, km)*dxf(i)))*dxhi(i)*dzhiq(k)
         iout1(i, j, k) = iout1(i, j, k) + ((utang1(i, j, k) - utang1(i, j, km))*emom*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k))*0.5 !
      END DO

!u top edge east
  DO j = jl, ju
         utang1Int = utang1(i, j, k)
         dummy = (utang1Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)*0.5
         emom = (dzf(km)*(ekm(i, j, k)*dxf(i - 1) + ekm(i - 1, j, k)*dxf(i)) + &
                 dzf(k)*(ekm(i, j, km)*dxf(i - 1) + ekm(i - 1, j, km)*dxf(i)))*dxhi(i)*dzhiq(k)
         iout1(i, j, k) = iout1(i, j, k) + ((utang1(i, j, k) - utang1(i, j, km))*emom*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k))*0.5 !
      END DO

!v
      il = block(n, 1)
      jl = MAX(block(n, 3) - myid*jmax, 1) + 1
      DO j = jl, ju
         DO i = il, iu

            utang2Int = utang2(i, j, k)
            ctm = fkar2/(logdz2)
            dummy = (utang2Int**2)*ctm
            bcmomflux = SIGN(dummy, utang2Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
            eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
            iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
         END DO
      END DO

!v top edge south
      j = MAX(block(n, 3) - myid*jmax, 1)
      DO i=il,iu
         utang2Int = utang2(i, j, k)
         dummy = (utang2Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)*0.5
         eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
         iout2(i, j, k) = iout2(i, j, k) + ((utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k))*0.5 !
      END DO

!v top edge north
      j = MIN(block(n, 4) - myid*jmax, jmax) + 1

      DO i = il, iu
         utang2Int = utang2(i, j, k)
         dummy = (utang2Int**2)*fkar2/(logdz2)
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)*0.5
         eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
         iout2(i, j, k) = iout2(i, j, k) + ((utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k))*0.5 !
      END DO
end if
!!!!!!!!!!!!!!!SPECIAL CASES FOR THE SURFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!can actually be made redundant and just be replaced by standard horizontal case (doesn't really matter though)
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
            ctm = fkar2/(logdz2)
            dummy = (utang1Int**2)*ctm
            bcmomflux = SIGN(dummy, utang1Int) 

            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)

            emom = (dzf(km)*(ekm(i, j, k)*dxf(i - 1) + ekm(i - 1, j, k)*dxf(i)) + & ! dx is non-equidistant
                    dzf(k)*(ekm(i, j, km)*dxf(i - 1) + ekm(i - 1, j, km)*dxf(i)))*dxhi(i)*dzhiq(k)

            iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(i, j, km))*emom*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !

         END DO
      END DO

      DO j = jl, ju !v component
         DO i = il, iu

            utang2Int = utang2(i, j, k)
            ctm = fkar2/(logdz2)
            dummy = (utang2Int**2)*ctm 
            bcmomflux = SIGN(dummy, utang2Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
            eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)

            iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
         END DO
      END DO

      END SELECT

END SUBROUTINE wfmneutral
