SUBROUTINE wfGR(hi,hj,hk,ioq,ioqflux,icth,obcqfluxA,qcell,qwall,hurel,resc,ress,n,ind,wforient)
   !wfGR
   USE modglobal, ONLY : dzf,dzfi,dzh2i,dzhi,dzhiq,dy,dyi,dy2i,dyi5,dxf,dxh,dxfi,dxhi,dxh2i,ib,ie,jb,je,kb,ke,fkar,grav,jmax,rk3step
   USE modsubgriddata, ONLY:ekh
   USE modmpi, ONLY:myid
   USE initfac, ONLY:block
   USE modibmdata
   INTEGER i, j, k, jl, ju, kl, ku, il, iu, km, im, jm, ip, jp, kp

   REAL :: bcqflux = 0. !temp storage for temperature flux
   REAL :: bcmomflux = 0. !temp storage for momentum flux
   REAL :: dummy = 0. !for debugging
   REAL :: delta = 0. !distance from wall
   REAL :: fkar2 = fkar**2 !fkar^2, von Karman constant squared
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

      delta = dxf(i)*0.5 !
      DO k = kl, ku
         DO j = jl, ju
            bcqflux=min(0.,cveg*(qcell(i,j,k)-qwall)*1/(1/icth(i,j,k)+resc)+(1-cveg)*(qcell(i,j,k)-qwall*hurel)*1/(1/icth(i,j,k)+ress))

            obcqfluxA = obcqfluxA + bcqflux
            ioqflux(i, j, k) = ioqflux(i, j, k) + bcqflux*dxfi(i)  

            ioq(i,j,k)=ioq(i,j,k)-0.5*(ekh(ip,j,k)*dxf(i)+ekh(i,j,k)*dxf(ip))*(qcell(ip,j,k)-qcell(i,j,k))*dxh2i(ip)*dxfi(i)-bcqflux*dxfi(i) !

         END DO
      END DO

!!! case 22 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   CASE (22)
      i = block(n, 2) + 1 !
      im = i - 1 !
      jl = MAX(block(n, 3) - myid*jmax, 1) ! starting j-index
      ju = MIN(block(n, 4) - myid*jmax, jmax) ! ending j-index
      kl = block(n, 5) ! starting k-index
      ku = block(n, 6) ! ending k-index

      delta = dxh(i)*0.5
      DO k = kl, ku
         DO j = jl, ju

            !dq * 1/res, where res is in [s/m]
            bcqflux=min(0.,cveg*(qcell(i,j,k)-qwall)*1/(1/icth(i,j,k)+resc)+(1-cveg)*(qcell(i,j,k)-qwall*hurel)*1/(1/icth(i,j,k)+ress))

            obcqfluxA = obcqfluxA + bcqflux
            ioqflux(i, j, k) = ioqflux(i, j, k) + bcqflux*dxfi(i) 

           ioq(i,j,k) = ioq(i,j,k) +0.5*(ekh(i,j,k)*dxf(im)+ekh(im,j,k)*dxf(i))*(qcell(i,j,k)-qcell(im,j,k))*dxh2i(i) * dxfi(i) - bcqflux*dxfi(i)

         END DO
      END DO
!!! case 32 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   CASE (32) !wall in xz -> wf in y (=vertical) upper, north wall

      j = ind
      jm = j - 1
      il = block(n, 1)
      iu = block(n, 2)
      kl = block(n, 5)
      ku = block(n, 6)

      DO k = kl, ku
         DO i = il, iu
            
            bcqflux=min(0.,cveg*(qcell(i,j,k)-qwall)*1/(1/icth(i,j,k)+resc)+(1-cveg)*(qcell(i,j,k)-qwall*hurel)*1/(1/icth(i,j,k)+ress))

            obcqfluxA = obcqfluxA + bcqflux
            ioqflux(i, j, k) = ioqflux(i, j, k) + bcqflux*dyi  

            ioq(i, j, k) = ioq(i, j, k) + ( &
                           0.5*(ekh(i, j, k) + ekh(i, jm, k))*(qcell(i, j, k) - qcell(i, jm, k)))*dy2i &
                           -bcqflux*dyi
         END DO
      END DO

!!! case 42 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   CASE (42) !wall in xz -> wf in y (=vertical) lower, south wall

      j = ind
      jp = j + 1
      il = block(n, 1)
      iu = block(n, 2)
      kl = block(n, 5)
      ku = block(n, 6)

      DO k = kl, ku
         DO i = il, iu

            bcqflux=min(0.,cveg*(qcell(i,j,k)-qwall)*1/(1/icth(i,j,k)+resc)+(1-cveg)*(qcell(i,j,k)-qwall*hurel)*1/(1/icth(i,j,k)+ress))

            obcqfluxA = obcqfluxA + bcqflux
            ioqflux(i, j, k) = ioqflux(i, j, k) + bcqflux*dyi 

            ioq(i, j, k) = ioq(i, j, k) - &
                           0.5*(ekh(i, jp, k) + ekh(i, j, k))*(qcell(i, jp, k) - qcell(i, j, k))*dy2i &
                           -bcqflux*dyi
         END DO
      END DO

!!! case 52 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   CASE (52) !wall in xy -> wf in z (=horizontal), top wall

      k = block(n, 6) + 1 !block location
      km = k - 1 !
      il = block(n, 1)
      iu = block(n, 2)
      jl = MAX(block(n, 3) - myid*jmax, 1)
      ju = MIN(block(n, 4) - myid*jmax, jmax)

      delta = dzf(k)*0.5

      DO j = jl, ju
         DO i = il, iu
            bcqflux=min(0.,cveg*(qcell(i,j,k)-qwall)*1/(1/icth(i,j,k)+resc)+(1-cveg)*(qcell(i,j,k)-qwall*hurel)*1/(1/icth(i,j,k)+ress))

            obcqfluxA = obcqfluxA + bcqflux
            ioqflux(i, j, k) = ioqflux(i, j, k) + bcqflux*dzfi(k)

            ioq(i, j, k) = ioq(i, j, k) &
                           + 0.5*(dzf(km)*ekh(i, j, k) + dzf(k)*ekh(i, j, km))*(qcell(i, j, k) - qcell(i, j, km))*dzh2i(k)*dzfi(k) &
                           - bcqflux*dzfi(k)

         END DO
      END DO
END SELECT

END SUBROUTINE wfGR
