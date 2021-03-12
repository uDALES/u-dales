!! modinlet.f90 contains the method of Lund (1998) to a generate a turbulent inlet profile.
!! The velocity is extracted from a recycle plane, rescaled, and used as inlet condition
!! Note that due to the staggered grid arrangement the u-components are recycled from 
!! cell(irecycl,:,:), while the v- and w- components are read from cell (irecycle-1,:,:). 
!! This is because the inlet condition is located at i=ib for u, and at i=ib-1 for v and w.
!!
!! Also the method of Kong (2000) is added to generate a turbulent temperature inlet profile
!!
!!  \author Jasper Tomas,TU Delft, June 4th 2015
!!  \par Revision list
!!  \todo Documentation
!!
module modinlet
use modinletdata
implicit none
save
  public :: initinlet,exitinlet,momentumthickness,blthickness,dispthickness,writeinletfile,readinletfile,enthalpythickness,inletgen,inletgennotemp,zinterpolate1d,zinterpolatew1d,zinterpolatet1d,zinterpolate2d,blthicknesst,momentumthicknessexp,dispthicknessexp

contains
  subroutine initinlet
    use modglobal, only : ih,ib,ie,jh,jb,je,kb,ke,kh,iinletgen,iplane,xf,lstoreplane,nstore,Uinf,ltempeq,pi,zf,zh
    use modfields, only : um
    use modmpi, only : myid,nprocs
   
    implicit none
    real    :: pfi, epsi
    integer :: k


  if (iinletgen==1) then

  allocate(Utav(ib:ie,kb:ke))
  allocate(Uinl(kb:ke))
  allocate(Winl(kb:ke+1))
  allocate(Urec(kb:ke))
  allocate(Wrec(kb:ke+1))
  allocate(u0inletbc(jb:je,kb:ke))
  allocate(v0inletbc(jb:je,kb:ke))
  allocate(w0inletbc(jb:je,kb:ke+1))
  allocate(u0inletbcold(jb:je,kb:ke))
  allocate(v0inletbcold(jb:je,kb:ke))
  allocate(w0inletbcold(jb:je,kb:ke+1))
  allocate(uminletbc(jb:je,kb:ke))
  allocate(vminletbc(jb:je,kb:ke))
  allocate(wminletbc(jb:je,kb:ke+1))
  allocate(uaver(ib:ie,kb:ke))
  allocate(zirf(kb:ke))
  allocate(ziif(kb:ke))
  allocate(zirh(kb:ke+1))
  allocate(ziih(kb:ke+1))
  allocate(zorf(kb:ke))
  allocate(zoif(kb:ke))
  allocate(zorh(kb:ke+1))
  allocate(zoih(kb:ke+1))
  allocate(loclowif(kb:ke))
  allocate(locupif(kb:ke))
  allocate(loclowih(kb:ke+1))
  allocate(locupih(kb:ke+1))
  allocate(loclowof(kb:ke))
  allocate(locupof(kb:ke))
  allocate(loclowoh(kb:ke+1))
  allocate(locupoh(kb:ke+1))
  allocate(displ(ib:ie))
  allocate(displold(ib:ie))
  allocate(upupavinl(kb:ke))
  allocate(vpvpavinl(kb:ke))
  allocate(wpwpavinl(kb:ke))
  allocate(upwpavinl(kb:ke))
  allocate(thlpthlpavinl(kb:ke))
  allocate(thlpupavinl(kb:ke))
  allocate(thlpwpavinl(kb:ke))
  allocate(heavif(kb:ke))
  allocate(heavih(kb:ke+1))

  if (lstoreplane ) then
    allocate(storeu0inletbc(jb:je,kb:ke,1:nstore))
    allocate(storev0inletbc(jb:je,kb:ke,1:nstore))
    allocate(storew0inletbc(jb:je,kb:ke+1,1:nstore))
  end if 

  epsi = 0.25*di
  do k=kb,ke
    pfi = zf(k) -1.2*di -epsi
    if (pfi < -epsi) then
      heavif(k) = 1.
    elseif(pfi <= epsi) then
      heavif(k) = 0.5* ( 1. - (pfi / epsi) - (1./pi)*sin(pi*pfi/epsi))
    elseif (pfi > epsi) then
      heavif(k) = 0.
    end if  
  end do

  do k=kb,ke+1
    pfi = zh(k) -1.2*di -epsi
    if (pfi < -epsi) then
      heavih(k) = 1.
    elseif(pfi <= epsi) then
      heavih(k) = 0.5* ( 1. - (pfi / epsi) - (1./pi)*sin(pi*pfi/epsi))
    elseif (pfi > epsi) then
      heavih(k) = 0.
    end if
  end do



  if (ltempeq  ) then
    allocate(Ttav(ib:ie,kb:ke))
    allocate(taver(ib:ie,kb:ke))
    allocate(Tinl(kb:ke))
    allocate(Trec(kb:ke))
    allocate(t0inletbc(jb:je,kb:ke))
    allocate(t0inletbcold(jb:je,kb:ke))
    allocate(tminletbc(jb:je,kb:ke))
    allocate(zotr(kb:ke))
    allocate(zoti(kb:ke))
    allocate(loclowot(kb:ke))
    allocate(locupot(kb:ke))
    allocate(heavit(kb:ke))
  
    if (lstoreplane ) then
      allocate(storet0inletbc(jb:je,kb:ke,1:nstore))
    end if
! Heaviside function for temperature
    epsi = 0.25*dti
    do k=kb,ke
      pfi = zf(k) -1.2*dti -epsi
      if (pfi < -epsi) then
        heavit(k) = 1.
      elseif(pfi <= epsi) then
        heavit(k) = 0.5* ( 1. - (pfi / epsi) - (1./pi)*sin(pi*pfi/epsi))
      elseif (pfi > epsi) then
        heavit(k) = 0.
      end if
    end do

 
  end if

  displ=0.
  displold =0.

  irecy = ib+iplane        ! index of recycle plane equals iplane (read from namoptions)


  xfm  = sum(xf(ib:ie))/(ie-ib+1)            ! mean(xf) 
  xf2m = sum(xf(ib:ie)**2.)/(ie-ib+1)        ! mean(xf^2)
!  btime = timee                              ! this is done to make sure btime is set when avint is computed correctly at startup (only for RA) 

  
  else if (iinletgen == 2) then
    allocate(storeu0inletbc(jb:je,kb:ke,1:nstore))
    allocate(storev0inletbc(jb:je,kb:ke,1:nstore))
    allocate(storew0inletbc(jb:je,kb:ke+1,1:nstore))
    allocate(u0rot(1:nstore,jb-jh:je+jh,kb:ke))
    allocate(v0rot(1:nstore,jb-jh:je+jh,kb:ke))
    allocate(u0inletbc(jb:je,kb:ke))
    allocate(v0inletbc(jb:je,kb:ke))
    allocate(w0inletbc(jb:je,kb:ke+1))
    allocate(u0inletbcold(jb:je,kb:ke))
    allocate(v0inletbcold(jb:je,kb:ke))
    allocate(w0inletbcold(jb:je,kb:ke+1))
    allocate(uminletbc(jb:je,kb:ke))
    allocate(vminletbc(jb:je,kb:ke))
    allocate(wminletbc(jb:je,kb:ke+1))
    if (ltempeq ) then
      allocate(storet0inletbc(jb:je,kb:ke,1:nstore))
      allocate(t0inletbc(jb:je,kb:ke))
      allocate(t0inletbcold(jb:je,kb:ke))
      allocate(tminletbc(jb:je,kb:ke))
    end if
    !iangle = iangledeg * pi / 180.  ! convert degrees to radians
    irecy = ib+iplane
  ! read coordinates of inletprofile  
    call readzincoord
!    ddispdx      = 0.00038/Uinf        ! this value should becomputed from the w0 computed in the inletgenerator
    ddispdx      = wtop/Uinf            ! wtop is read from zgrid.inf 
    ddispdxold   = ddispdx              ! this value should becomputed from the w0 computed in the inletgenerator
!    inlfactor    = nprocs/nprocsinl     ! nprocs should be larger or equal to nprocsin!
!    write(6,*) 'inlfactor= ',inlfactor
  else
   return
  end if
 
  end subroutine initinlet

  subroutine inletgen
    use modglobal,   only : ib,ie,jb,je,jgb,jge,kb,ke,zf,zh,dzf,dzhi,timee,btime,totavtime,rk3step,dt,numol,iplane,lles,iinletgen,inletav,runavtime,Uinf,lwallfunc,linletRA,totinletav,lstoreplane,nstore,prandtlmoli,numol,grav,lbuoyancy,lfixinlet,luvolflowr,lfixutauin
    use modfields,   only : u0,v0,w0,thl0,wm,uprof
    use modsurfdata, only : thls,thl_top
    use modsave,     only : writerestartfiles
    use modmpi,      only : slabsum,myid
    implicit none
   
    real,dimension(ib:ib,jb:je,kb:ke)   :: uinletbc2   ! dummy variable
    real,dimension(ib:ib,jb:je,kb:ke)   :: tinletbc2   ! dummy variable
    real,dimension(jb:je,kb:ke)   :: uprec             ! velocity fluctuation (up_rec = u0 - Urec)
    real,dimension(jb:je,kb:ke)   :: vprec             ! velocity fluctuation (vp_rec = v0 - 0)
    real,dimension(jb:je,kb:ke+1) :: wprec             ! velocity fluctuation (wp_rec = w0 - Wrec)
    real,dimension(jb:je,kb:ke)   :: tprec             ! temperature fluctuation (tp_rec = t0 - Trec)
    real,dimension(jb:je,kb:ke)   :: upinli,vpinli     ! = gamma * (uprec,v interpolated to zii grid)
    real,dimension(jb:je,kb:ke)   :: tpinli            ! = lambda  * (tprec   interpolated to zii grid)
    real,dimension(jb:je,kb:ke)   :: upinlo,vpinlo     ! = gamma * (uprec,v interpolated to zoi grid)
    real,dimension(jb:je,kb:ke)   :: tpinlo            ! = lambda  * (tprec   interpolated to zoti grid)
    real,dimension(jb:je,kb:ke+1) :: wpinli            ! = gamma * (wprec   interpolated to zii grid)
    real,dimension(jb:je,kb:ke+1) :: wpinlo            ! = gamma * (wprec   interpolated to zoi grid)
    real,dimension(kb:ke)   :: udiff                   ! difference between Uinl and Urec
!    real,dimension(kb:ke)   :: Urecdiff                ! difference between Urec new and old
    real,dimension(kb:ke)   :: urav                    ! j-averaged u-velocity (not time-averaged) 
    real,dimension(kb:ke)   :: trav                    ! j-averaged temperature (not time-averaged) 
    real,dimension(kb:ke)   :: uravdzf                 ! j-averaged u-velocity (not time-averaged) times dzf
    real,dimension(kb:ke)   :: uinldzf                 ! j-averaged u-velocity (not time-averaged) times dzf
    real,dimension(kb:ke)   :: Urecdzf                 ! Urec times dzf
    real,dimension(kb:ke+1) :: wrav                    ! j-averaged w-velocity (not time-averaged) 
    real,dimension(kb:ke)   :: Uinli                   ! = gamma * (Urec interpolated to ziif grid points)
    real,dimension(kb:ke+1) :: Winli                   ! = gamma * (Wrec interpolated to ziih grid points)
    real,dimension(kb:ke)   :: Tinli                   ! = lambda  * (Trec interpolated to ziif grid points)
    real,dimension(kb:ke)   :: Uinlo                   ! = gamma * (Urec interpolated to zioif grid points)
    real,dimension(kb:ke+1) :: Winlo                   ! = gamma * (Wrec interpolated to zoih grid points)
    real,dimension(kb:ke)   :: Tinlo                   ! = lambda  * (Trec interpolated to zoti grid points)
    real,dimension(kb:ke)   :: wfuncf                  ! weight function at full level
    real,dimension(kb:ke+1) :: wfunch                  ! weight function at half level
    real,dimension(kb:ke)   :: wfunct                  ! weight function at full level
    real                    :: utaur2,utaui2           ! (utau)^2 at recycle station and inlet
    real                    :: gamm                    ! utaui / utaur
    real                    :: lamb                    ! ttaui / ttaur
    real                    :: avint,avinti            ! avering interval
    real                    :: alpha,beta              ! factors used in the Weight function
!    real                    :: totalu                  ! total u-velocity at outlet
    real                    :: Urectot                  ! total u-velocity at recycle plane
    real                    :: rk3coef
!    real                    :: di_test                 ! BL thickness as measured from Uinl
    real                    :: utop                    ! j-averaged top velocity
    real                    :: interval
    real                    :: dtinrk                  ! RK time step in inlet data
    real                    :: rk3coefin               ! Cumulative RK time step in inlet data
    real                    :: dr_old
    real                    :: scalef                      ! scale factor to scale instantaneous velocity profile with to get constant mass flux
    real                    :: totaluinl                   ! bulk velocity at the inlet
!    real                    :: q0                      ! heat flux
    integer i,j,k,kk,kdamp


   if (iinletgen == 1) then

   u0inletbcold = u0inletbc
   v0inletbcold = v0inletbc
   w0inletbcold = w0inletbc
   t0inletbcold = t0inletbc    ! temperature
   totaluold    = totalu
   displold     = displ
   ddispdxold   = ddispdx

   ! compute time-average velocities
    rk3coef = dt / (4. - dble(rk3step))
    if (rk3step==1) then
      deltat = rk3coef
    elseif (rk3step==2) then
      deltat = rk3coef - (dt/3.)
    elseif (rk3step==3) then
      deltat = rk3coef - (dt/2.)
    end if

    if (linletRA) then  ! this is a switch to use 'running average'
      avint = totinletav + timee-btime ! runav interval = averaging interval previuous sim  + current elapsed sim time
    else
      avint  = inletav
    end if
    avinti = 1./avint
    uaver=0. 
    taver=0. 
    do i=ib,ie
      call slabsum(uaver(i,:),kb,ke,  u0(i:i,jb:je,kb:ke),i,i,jb,je,kb,ke,i,i,jb,je,kb,ke)
      call slabsum(taver(i,:),kb,ke,thl0(i:i,jb:je,kb:ke),i,i,jb,je,kb,ke,i,i,jb,je,kb,ke)
    end do

    wrav=0.
    call slabsum(wrav(kb:ke+1),kb,ke,w0(irecy-1:irecy-1,jb:je,kb:ke+1),irecy-1,irecy-1,jb,je,kb,ke+1,irecy-1,irecy-1,jb,je,kb,ke+1)
    trav=0.
    call slabsum(trav(kb:ke),  kb,ke,thl0(irecy-1:irecy-1,jb:je,kb:ke), irecy-1,irecy-1,jb,je,kb,ke,  irecy-1,irecy-1,jb,je,kb,ke)

    uaver = uaver / (jge-jgb +1)                    ! average over j-direction
    taver = taver / (jge-jgb +1)                    ! average over j-direction
    urav = uaver(irecy,:)
    wrav = wrav / (jge-jgb +1)                    ! average over j-direction
    trav = trav / (jge-jgb +1)                    ! average over j-direction

    do k=kb,ke
      Urec(k) =  urav(k)*deltat*avinti + (1.-deltat*avinti)*Urec(k)
      Trec(k) =  trav(k)*deltat*avinti + (1.-deltat*avinti)*Trec(k)
    end do
    do k=kb,ke+1
      Wrec(k) =  wrav(k)*deltat*avinti + (1.-deltat*avinti)*Wrec(k)
    end do
    do k=kb,ke
    do i=ib,ie
      Utav(i,k) =  uaver(i,k)*deltat*avinti + (1.-deltat*avinti)*Utav(i,k)
      Ttav(i,k) =  taver(i,k)*deltat*avinti + (1.-deltat*avinti)*Ttav(i,k)
    end do
    end do
   
!    Urec = Urec +(Uinf-Urec(ke))     ! make sure at the recycle plane the top velocity equals Uinf
      

!    Urecdiff = Urecdiff - Urec 
!    if (myid==0) then
!      write(6,*) 'Urec_old - Urec_new (kb+40)=',Urecdiff(kb+40) 
!    end if

!! check if Urec contains NaN
!    if (myid==0) then
!      write(6,*) 'Checking Urec for NaN' 
!      do k=kb,ke
!        if (ISNAN(Urec(k))) then
!          write(6,*) 'Urec(k)=NaN at k=kb+', k-kb
!        end if
!      end do
!      write(6,*) 'Finished checking Urec for NaN'  
!    end if


!    if (myid==0) then
!      write(6,*) 'myid, Urec(ke)=',myid, Urec(ke)
!      write(6,*) 'wrav(ke), Wrec(ke)=',wrav(ke), Wrec(ke)
!      write(6,*) 'wrav(ke-1), Wrec(ke-1)=',wrav(ke-1), Wrec(ke-1)
!      write(6,*) 'wrav(ke-10), Wrec(ke-10)=',wrav(ke-10), Wrec(ke-10)
!      write(6,*) 'wrav(ke-30), Wrec(ke-30)=',wrav(ke-30), Wrec(ke-30)
!      write(6,*) 'wrav(kb+10), Wrec(kb+10)=',wrav(kb+10), Wrec(kb+10)
!      write(6,*) 'wrav(kb+11), Wrec(kb+11)=',wrav(kb+11), Wrec(kb+11)
!    end if


   ! compute velocity fluctuation at recycle station
    do k=kb,ke
      do j=jb,je
        uprec(j,k) = u0(irecy,j,k) - Urec(k)
        vprec(j,k) = v0(irecy-1,j,k)             ! mean v is zero
        tprec(j,k) = thl0(irecy-1,j,k) - Trec(k)
      end do
    end do

    do k=kb,ke+1
      do j=jb,je
        wprec(j,k) = w0(irecy-1,j,k) - Wrec(k)   ! note that w-velocity is taken at i=irecy-1 !!
      end do
    end do

 
    if (lwallfunc) then
      call wallawinlet(Urec(kb),dzf(kb),numol,utaur2)    ! compute wall shear stress at recycle station
    else
      utaur2 = 2.*numol*Urec(kb)/dzf(kb)
    end if
    utaur = sqrt(abs(utaur2))                          ! compute utau at recycle station
    ! heat flux at recycle station (isothermal wall) q = alpha * dT/dz = (nu/prandtl) * dT/dz 
!    q0 = numol*prandtlmoli*(Trec(kb) - Trec(kb-1)) * dzhi(kb) 
    q0 = numol*prandtlmoli*2*(Trec(kb) - thls) / dzf(kb) 
    ttaur = q0/utaur                ! ttau = q/(rho*cp*utau) =  (alpha *dT/dz) / utau
   ! compute momentum thickness at inlet and recycle plane
  
   if(lbuoyancy) then
     lmor = (thls* utaur**2 )/ (0.41 * grav * ttaur)      ! L = -T0*utau^3 / kappa*g*<w'T'> =  
!     write(6,*) 'Initial dr,myid, utaur, ttaur, Lmor =', dr,myid,utaur,ttaur,lmor
!     lmor = 0.3;
     lmoi = (thls* utaui**2 )/ (0.41 * grav * ttaui)      ! L = -T0*utau^3 / kappa*g*<w'T'> = 
!     lmoi = 0.3;
!     write(6,*) 'Initial di_test,myid, utaui, ttaui, Lmoi =', di_test,myid,utaui,ttaui,lmoi
     dr_old = dr

!     call blthicknessmo(dr,utaur,lmor)                    ! Also needed for momentumthickness
     call blthicknesst(dr,Urec,0.99)            ! changed back to this one (instead of the above)

!     call momentumthicknessmo(thetai,utaui,di,lmoi)
!     call momentumthicknessmo(thetar,utaur,dr,lmor)
     call momentumthicknessexp(thetai,Uinl)
     call momentumthicknessexp(thetar,Urec)
   else
!     call blthickness(dr,utaur)                           ! Also needed for momentumthickness
     call blthicknesst(dr,Urec,0.99)                        
!     call momentumthickness(thetai,utaui,di)
!     call momentumthickness(thetar,utaur,dr)
     call momentumthicknessexp(thetai,Uinl)
     call momentumthicknessexp(thetar,Urec)
   end if
   call enthalpythickness(thetati,Tinl,Uinl)
   call enthalpythickness(thetatr,Trec,Urec)
!   call blthickness(dr,utaur)
   call blthicknesst(dtr,Trec-thls,0.99)
   ! compute utau at inlet from interior field 
!    if (thetai == 0.) then
!      write(6,*) '!!! thetai = 0, myid=',myid  
!    else if (thetar == 0.) then
!      write(6,*) '!!! thetar = 0, myid=',myid      
!      thetar=0.00001 
!    else
!      utaui = utaur* (thetar/thetai)**(1./8.)    ! See Lund (1998): 'Similar to Ludwig-Tillmann correlation'
    if (.not.lfixutauin) then 
      utaui  = utaur* abs(thetar/thetai)**(1./8.)   ! See Lund (1998): 'Similar to Ludwig-Tillmann correlation'
    end if
    if (thetati == 0.) then
      thetati = 0.0000001
    end if  
    ttaui = ttaur*abs(thetatr/thetati)**(1./8.)   ! See Kong (2000):  
!    end if
    gamm = utaui / utaur                          ! Gamma in Lund (1998)
    if (ttaur == 0.) then
      ttaur = 0.0000001
    end if
    lamb = ttaui / ttaur                          ! Lambda in Kong (2000)

   ! compute inner scaling coordinates
    zirf = utaur*zf / numol                       ! inner scaling zf-coordinate at recycle station 
    zirh = utaur*zh / numol                       ! inner scaling zh-coordinate at recycle station
    ziif = utaui*zf / numol                       ! inner scaling zf-coordinate at inlet station 
    ziih = utaui*zh / numol                       ! inner scaling zh-coordinate at inlet station 
  
   ! compute outer scaling coordinates   
    zorf = zf / dr                                ! outer scaling zf-coor as measured from Uinldinate at recycle station 
    zorh = zh / dr                                ! outer scaling zh-coordinate at recycle station 
    zoif = zf / di                                ! outer scaling zf-coordinate at inlet station  (could be done once, actually..) 
    zoih = zh / di                                ! outer scaling zf-coordinate at inlet station  (could be done once, actually..)
    zotr = zf / dtr                               ! temperature outer scaling zf-coordinate at recycle station
    zoti = zf / dti                               ! temperature outer scaling zf-coordinate at inlet station

!!!!! Interpolation starts here

 !!! First inner coordinates
   ! determine which elements are needed when recycle velocity profile is interpolated on inlet plane
   ! for u(,v)-components (zf)
    do k=kb,ke
      do kk=kb,ke
        if (zirf(kk) >= ziif(k)) then
          locupif(k)  = kk
          loclowif(k) = kk-1
          exit
        elseif (kk==ke) then
          locupif(k)  = ke+1            ! this means extrapolation!
          loclowif(k) = ke-1            ! waarom niet ke? of wordt dit niet gebruikt?
        end if
      end do
    end do
   ! for w-components (zh)
    do k=kb,ke+1
      do kk=kb,ke+1
        if (zirh(kk) >= ziih(k)) then
          locupih(k)  = kk
          loclowih(k) = kk-1
          exit
        elseif (kk==ke+1) then
          locupih(k)  = ke+2            ! this means extrapolation!
          loclowih(k) = ke
        end if
      end do
    end do
 !!! Finished with inner coordinates

 !!! Do the same trick for outer coordinates
   ! determine which elements are needed when recycle velocity profile is interpolated on inlet plane
   ! for u(,v)-components (zf)
    do k=kb,ke
      do kk=kb,ke
        if (zorf(kk) >= zoif(k)) then
          locupof(k)  = kk
          loclowof(k) = kk-1
          exit
        elseif (kk==ke) then
          locupof(k)  = ke+1            ! this means extrapolation!
          loclowof(k) = ke-1
        end if
      end do
    end do
   ! for w-components (zh)
    do k=kb,ke+1
      do kk=kb,ke+1
        if (zorh(kk) >= zoih(k)) then
          locupoh(k)  = kk
          loclowoh(k) = kk-1
          exit
        elseif (kk==ke+1) then
          locupoh(k)  = ke+2            ! this means extrapolation!
          loclowoh(k) = ke
        end if
      end do
    end do

  !!! Finished with outer coordinates
  !!! Outer coordinates for temperature
    do k=kb,ke
      do kk=kb,ke
        if (zotr(kk) >= zoti(k)) then
          locupot(k)  = kk
          loclowot(k) = kk-1
          exit
        elseif (kk==ke) then
          locupot(k)  = ke+1            ! this means extrapolation!
          loclowot(k) = ke-1
        end if
      end do
    end do
  !!! Finished with outer coordinates temperature

!!! Now really interpolate
  !!! First inner coordinates
   ! compute Urec on zii grid
    do k=kb,ke
      if (locupif(k) == ke+1) then      ! indicator for extrapolation!
!        Uinli(k) = Urec(ke) + (Urec(ke) - Urec(ke-1)) / (zirf(ke)-zirf(ke-1)) * (ziif(k)-zirf(ke))
        Uinli(k) = Urec(ke)
        Tinli(k) = Trec(ke) 
      elseif (loclowif(k) == kb-1) then ! interprets this as extrapolation to bottom (use u=0 at z+=0)
        Uinli(k) = Urec(kb)/zirf(kb) * ziif(k)
!        Tinli(k) = thls + Trec(kb)/zirf(kb)*ziif(k)  
!        Tinli(k) = (Trec(kb)-thls)/zirf(kb)*ziif(k)  
        Tinli(k) = thls + (Trec(kb)-thls)/zirf(kb)*ziif(k)  
      else                            ! normal interpolation
        Uinli(k) = Urec(loclowif(k)) + (Urec(locupif(k)) - Urec(loclowif(k))) / (zirf(locupif(k)) - zirf(loclowif(k))) * (ziif(k) - zirf(loclowif(k)))
        Tinli(k) = Trec(loclowif(k)) + (Trec(locupif(k)) - Trec(loclowif(k))) / (zirf(locupif(k)) - zirf(loclowif(k))) * (ziif(k) - zirf(loclowif(k)))
        if ((ziif(k) .gt. zirf(locupif(k))) .or. (ziif(k) .lt. zirf(loclowif(k)))) then
          write(6,*) '!!!Mistake in Interpolation !!!!'
        end if
      end if
    end do

   ! compute Wrec on zii grid
    Winli(kb) = 0.0                      ! corresponds to ground level
    do k=kb+1,ke+1
      if (locupih(k) == ke+2) then     ! indicator for extrapolation!
!        Winli(k) = Wrec(ke+1) + (Wrec(ke+1) - Wrec(ke)) / (zirh(ke+1)-zirh(ke)) * (ziih(k)-zirh(ke+1))
        Winli(k) = Wrec(ke+1)
      else                            ! normal interpolation
        Winli(k) = Wrec(loclowih(k)) + (Wrec(locupih(k)) - Wrec(loclowih(k))) / (zirh(locupih(k)) - zirh(loclowih(k))) * (ziih(k) - zirh(loclowih(k)))
      end if
    end do

   ! compute u- and v- and t-fluctuation on zii grid
    do k=kb,ke
      if (locupif(k) == ke+1) then      ! indicator for extrapolation!
!        upinli(:,k) = uprec(:,ke) + (uprec(:,ke) - uprec(:,ke-1)) / (zirf(ke)-zirf(ke-1)) * (ziif(k)-zirf(ke))
        upinli(:,k) = 0.
        vpinli(:,k) = 0.
        tpinli(:,k) = 0.
      elseif (loclowif(k) == kb-1) then ! interprets this as extrapolation to bottom (use u=0 at z+=0)
        upinli(:,k) = uprec(:,kb)/zirf(kb) * ziif(k)
        vpinli(:,k) = vprec(:,kb)/zirf(kb) * ziif(k)
        tpinli(:,k) = tprec(:,kb)/zirf(kb) * ziif(k)
      else                            ! normal interpolation
        upinli(:,k) = uprec(:,loclowif(k)) + (uprec(:,locupif(k)) - uprec(:,loclowif(k))) / (zirf(locupif(k)) - zirf(loclowif(k))) * (ziif(k) - zirf(loclowif(k)))
        vpinli(:,k) = vprec(:,loclowif(k)) + (vprec(:,locupif(k)) - vprec(:,loclowif(k))) / (zirf(locupif(k)) - zirf(loclowif(k))) * (ziif(k) - zirf(loclowif(k)))
        tpinli(:,k) = tprec(:,loclowif(k)) + (tprec(:,locupif(k)) - tprec(:,loclowif(k))) / (zirf(locupif(k)) - zirf(loclowif(k))) * (ziif(k) - zirf(loclowif(k)))
      end if
    end do
   
   ! compute w-fluctuation on zii grid
    do k=kb+1,ke+1
!      if (locupih(k) == ke+1) then      ! indicator for extrapolation!
      if (locupih(k) == ke+2) then      ! indicator for extrapolation!
!        wpinli(:,k) = wprec(:,ke+1) + (wprec(:,ke+1) - wprec(:,ke)) / (zirh(ke+1)-zirh(ke)) * (ziih(k)-zirh(ke+1))
        wpinli(:,k) = 0.
      else                            ! normal interpolation
        wpinli(:,k) = wprec(:,loclowih(k)) + (wprec(:,locupih(k)) - wprec(:,loclowih(k))) / (zirh(locupih(k)) - zirh(loclowih(k))) * (ziih(k) - zirh(loclowih(k)))
      end if
    end do

 !! Finished with interpolating inner variables 
 !! Continue with interpolating outer variables
   ! compute Urec on zoi grid
    do k=kb,ke
      if (locupof(k) == ke+1) then      ! indicator for extrapolation!
!        Uinlo(k) = Urec(ke) + (Urec(ke) - Urec(ke-1)) / (zorf(ke)-zorf(ke-1)) * (zoif(k)-zorf(ke))
!        Uinlo(k) = Urec(ke)
        Uinlo(k) = Uinf
      elseif (loclowof(k) == kb-1) then ! interprets this as extrapolation to bottom (use u=0 at z+=0)
        Uinlo(k) = Urec(kb)/zorf(kb) * zoif(k)
      else                            ! normal interpolation
        Uinlo(k) = Urec(loclowof(k)) + (Urec(locupof(k)) - Urec(loclowof(k))) / (zorf(locupof(k)) - zorf(loclowof(k))) * (zoif(k) - zorf(loclowof(k)))
      end if
    end do
  
   ! compute Wrec on zii grid
    Winlo(kb) = 0.0                      ! corresponds to ground level
    do k=kb+1,ke+1
      if (locupoh(k) == ke+2) then     ! indicator for extrapolation!
!        Winlo(k) = Wrec(ke+1) + (Wrec(ke+1) - Wrec(ke)) / (zorh(ke+1)-zorh(ke)) * (zoih(k)-zorh(ke+1))
        Winlo(k) = Wrec(ke+1) 
      else                            ! normal interpolation
        Winlo(k) = Wrec(loclowoh(k)) + (Wrec(locupoh(k)) - Wrec(loclowoh(k))) / (zorh(locupoh(k)) - zorh(loclowoh(k))) * (zoih(k) - zorh(loclowoh(k)))
      end if
    end do

   ! compute u- and v-fluctuation on zoi grid
    do k=kb,ke
      if (locupof(k) == ke+1) then      ! indicator for extrapolation!
!        upinlo(:,k) = uprec(:,ke) + (uprec(:,ke) - uprec(:,ke-1)) / (zorf(ke)-zorf(ke-1)) * (zoif(k)-zorf(ke))
        upinlo(:,k) = 0.
        vpinlo(:,k) = 0.
      elseif (loclowof(k) == kb-1) then ! interprets this as extrapolation to bottom (use u=0 at z+=0)
        upinlo(:,k) = uprec(:,kb)/zorf(kb) * zoif(k)
        vpinlo(:,k) = vprec(:,kb)/zorf(kb) * zoif(k)
      else                            ! normal interpolation
        upinlo(:,k) = uprec(:,loclowof(k)) + (uprec(:,locupof(k)) - uprec(:,loclowof(k))) / (zorf(locupof(k)) - zorf(loclowof(k))) * (zoif(k) - zorf(loclowof(k)))
        vpinlo(:,k) = vprec(:,loclowof(k)) + (vprec(:,locupof(k)) - vprec(:,loclowof(k))) / (zorf(locupof(k)) - zorf(loclowof(k))) * (zoif(k) - zorf(loclowof(k)))
      end if
    end do

   ! compute w-fluctuation on zoi grid
    do k=kb+1,ke+1
      if (locupoh(k) == ke+2) then      ! indicator for extrapolation!
!        wpinlo(:,k) = wprec(:,ke+1) + (wprec(:,ke+1) - wprec(:,ke)) / (zorh(ke+1)-zorh(ke)) * (zoih(k)-zorh(ke+1))
        wpinlo(:,k) = 0.
      else                            ! normal interpolation
        wpinlo(:,k) = wprec(:,loclowoh(k)) + (wprec(:,locupoh(k)) - wprec(:,loclowoh(k))) / (zorh(locupoh(k)) - zorh(loclowoh(k))) * (zoih(k) - zorh(loclowoh(k)))
      end if
    end do
 !! Finished interpolating outer velocity variables
 !! Interpolating outer temperature
   ! mean temperature
    do k=kb,ke
      if (locupot(k) == ke+1) then      ! indicator for extrapolation!
        Tinlo(k) = thl_top
      elseif (loclowot(k) == kb-1) then ! interprets this as extrapolation to bottom (use Tinlo=thls at z+=0)
!        Tinlo(k) = Trec(kb)/zotr(kb) * zoti(k)
!        Tinlo(k) = (Trec(kb)-thls)/zotr(kb) * zoti(k)
        Tinlo(k) = thls + (Trec(kb)-thls)/zotr(kb) * zoti(k)
      else                            ! normal interpolation
        Tinlo(k) = Trec(loclowot(k)) + (Trec(locupot(k)) - Trec(loclowot(k))) / (zotr(locupot(k)) - zotr(loclowot(k))) * (zoti(k) - zotr(loclowot(k)))
      end if
    end do
   
   ! fluctuating temperature
    do k=kb,ke
      if (locupot(k) == ke+1) then      ! indicator for extrapolation!
!        upinlo(:,k) = uprec(:,ke) + (uprec(:,ke) - uprec(:,ke-1)) / (zorf(ke)-zorf(ke-1)) * (zoif(k)-zorf(ke))
        tpinlo(:,k) = 0.
      elseif (loclowot(k) == kb-1) then ! interprets this as extrapolation to bottom (use t=0 at z+=0)
        tpinlo(:,k) = tprec(:,kb)/zotr(kb) * zoti(k)
      else                            ! normal interpolation
        tpinlo(:,k) = tprec(:,loclowot(k)) + (tprec(:,locupot(k)) - tprec(:,loclowot(k))) / (zotr(locupot(k)) - zotr(loclowot(k))) * (zoti(k) - zotr(loclowot(k)))
      end if
    end do
     
 !! Finished interpolating out temperature
!!!!! Finished Interpolation! !!!!!


   ! compute rescaled inner variables ! Winli = Winli (interpolation is enough)
    Uinli = gamm* Uinli                   
    Tinli = lamb* Tinli + (1.-lamb)*thls     ! this is different for isoflux wall!                   
    upinli = gamm* upinli
    vpinli = gamm* vpinli
    wpinli = gamm* wpinli 
    tpinli = lamb* tpinli         ! See Kong (2000)
   ! compute rescaled outer variables ! Winlo = Winlo (interpolation is enough)
    Uinlo = gamm* Uinlo  + (1.- gamm)*Uinf
    Tinlo = lamb* Tinlo  + (1.- lamb)*thl_top 
!    Uinlo = gamm* Uinlo  + (1.- gamm)*Urec(ke) 
    upinlo = gamm* upinlo
    vpinlo = gamm* vpinlo
    wpinlo = gamm* wpinlo
    tpinlo = lamb* tpinlo         ! See Kong (2000)

!    utop = Uinlo(ke)
!    Uinlo = Uinlo +(Uinf-utop)     ! make sure at the inlet the mean top velocity equals Uinf
!! add defect velocity to make sure the j-averaged velocity at the top equals Uinf
!    utop = Uinlo(ke)
!    do k=kb,ke
!        Uinlo(k) = Uinlo(k)*Uinf/utop
!    end do


   ! Compute weight function (alpha=4, b=0.2)
    alpha = 4.
    beta = 0.2
    wfuncf = 0.5*(1. + tanh( alpha*(zoif-beta)/((1.-2.*beta)*zoif+beta) )/tanh(alpha) ) ! for full level height
    wfunch = 0.5*(1. + tanh( alpha*(zoih-beta)/((1.-2.*beta)*zoih+beta) )/tanh(alpha) ) ! for half level height
    wfunct = 0.5*(1. + tanh( alpha*(zoti-beta)/((1.-2.*beta)*zoti+beta) )/tanh(alpha) ) ! for temperature (full level height)
    do k=kb,ke
      if (wfuncf(k) .gt. 1.) then  
        wfuncf(k) = 1.
      end if
      if (wfunct(k) .gt. 1.) then  
        wfunct(k) = 1.
      end if
    end do
    do k=kb,ke+1
      if (wfunch(k) .gt. 1.) then  
        wfunch(k) = 1.
      end if
    end do
    

 
!    write(6,*) 'maxval(wfuncf)=', maxval(wfuncf)
!    write(6,*) 'maxval(wfunch)=', maxval(wfunch)


   ! Compute the velocity components for the inlet BC
    do k=kb,ke  
    do j=jb,je
!      u0inletbc(j,k) = (Uinli(k)+ upinli(j,k))*(1.-wfuncf(k)) +  (Uinlo(k) + upinlo(j,k))* wfuncf(k)
!      v0inletbc(j,k) =            vpinli(j,k) *(1.-wfuncf(k)) +              vpinlo(j,k) * wfuncf(k) 
!      t0inletbc(j,k) = (Tinli(k)+ tpinli(j,k))*(1.-wfunct(k)) +  (Tinlo(k) + tpinlo(j,k))* wfunct(k)      
      u0inletbc(j,k) = (Uinli(k)+ upinli(j,k)*heavif(k))*(1.-wfuncf(k)) +  (Uinlo(k) + upinlo(j,k)*heavif(k))* wfuncf(k)
      v0inletbc(j,k) =            vpinli(j,k)*heavif(k) *(1.-wfuncf(k)) +              vpinlo(j,k)*heavif(k) * wfuncf(k) 
      t0inletbc(j,k) = (Tinli(k)+ tpinli(j,k)*heavit(k))*(1.-wfunct(k)) +  (Tinlo(k) + tpinlo(j,k)*heavit(k))* wfunct(k)      
    end do
    end do
   
    do k=kb,ke+1 
    do j=jb,je
      w0inletbc(j,k) = (Winli(k)+ wpinli(j,k)*heavih(k))*(1-wfunch(k)) +  (Winlo(k) + wpinlo(j,k)*heavih(k))* wfunch(k)
    end do
    end do
    w0inletbc(:,kb) = 0.
    w0inletbc(:,ke+1) = 0.

!!    kdamp = kb + floor(0.75*(ke-kb+1))
!    kdamp = kb + 144  ! => zf = 2.24
!    do k=kdamp,ke  
!    do j=jb,je
!      if (u0inletbc(j,k) > Uinf) then
!        u0inletbc(j,k) = Uinf
!      end if
!    end do
!    end do



   ! Compute j-averaged inlet U  (used for compute thetai)
    uinletbc2(ib,jb:je,kb:ke) = u0inletbc(jb:je,kb:ke)  ! this is just a dummy variable to give uninletbc the right dimension in slabsum
    tinletbc2(ib,jb:je,kb:ke) = t0inletbc(jb:je,kb:ke)  ! this is just a dummy variable to give tninletbc the right dimension in slabsum
    urav = 0.
    trav = 0.
    call slabsum(urav  ,kb,ke,uinletbc2 ,ib,ib,jb,je,kb,ke,ib,ib,jb,je,kb,ke)
    call slabsum(trav  ,kb,ke,tinletbc2 ,ib,ib,jb,je,kb,ke,ib,ib,jb,je,kb,ke)
!    call slabsum(urav  ,kb,ke,u0 ,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,ib,ib,jb,je,kb,ke)
    urav = urav / (jge-jgb +1)                    ! average over j-direction
    trav = trav / (jge-jgb +1)                    ! average over j-direction

! determine bulk velocity of new profile
    do k=kb,ke
      uravdzf(k) = urav(k)*dzf(k)
    end do
    totalu    = sum(uravdzf(kb:ke))/(zh(ke+1)-zh(kb))      ! Area-averaged outflow velocity

! rescale the instantaneous profile to keep mass flux constant (tot avoid pressure fluctuations)
    if (luvolflowr ) then 
      do k=kb,ke
        uinldzf(k) = Uinl(k)*dzf(k)
      end do
      totaluinl = sum(uinldzf(kb:ke))/(zh(ke+1)-zh(kb))      ! Area-averaged inflow velocity
      scalef = totaluinl/totalu ! compute factor to scale the velocity profile with
      u0inletbc(:,:) = u0inletbc(:,:)*scalef  ! rescale the velocity profile to have constant mass-flux 
      urav(:) = urav(:)*scalef                ! also rescale the part that is added to the mean
    end if

!! add defect velocity to make sure the mass flow is the same as the initial mass flow
 !   u0inletbc = u0inletbc + (ubulk-totalu)
 !   urav      = urav      + (ubulk-totalu)

!! add defect velocity to make sure the j-averaged velocity at the top equals Uinf
!    utop = urav(ke)
!    do k=kb,ke
!      do j=jb,je
!        u0inletbc(j,k) = u0inletbc(j,k)*Uinf/utop
!      end do
!      urav(k) = urav(k)*Uinf/utop
!    end do
   
!    u0inletbc = u0inletbc + (Uinf-utop)
!    urav      = urav      + (Uinf-utop)


!    if (myid==0) then
!    write(6,*) 'u0inletbc(jb+2,ke)', u0inletbc
!    end if

   ! Compute j- and time-averaged  inlet U  (used for compute thetai)
    if (.not.lfixinlet) then  ! only update the average inlet profiles when lfixinlet .eqv..false.
      do k=kb,ke
        Uinl(k) =  urav(k)*deltat*avinti + (1.-deltat*avinti)*Uinl(k)
      end do
    end if
    do k=kb,ke
      Tinl(k) =  trav(k)*deltat*avinti + (1.-deltat*avinti)*Tinl(k)
    end do

!    utop = Uinl(ke)
!    Uinl = Uinl +(Uinf-utop)     ! make sure at the inlet the mean top velocity equals Uinf
!    uminletbc = uminletbc + (Uinf-utop)

! write inletplane to array (and to file after 1000 time steps)   
    if (lstoreplane ) then     
      storeu0inletbc(:,:,nstepread) = u0inletbc(:,:) 
      storev0inletbc(:,:,nstepread) = v0inletbc(:,:)
      storew0inletbc(:,:,nstepread) = w0inletbc(:,:)
      storet0inletbc(:,:,nstepread) = t0inletbc(:,:)
      nstepread = nstepread +1
      if (nstepread == nstore+1) then
        nfile = nfile +1       ! next file number
        call writeinletfile    ! write 1000 time steps to file
        call writerestartfiles 
        nstepread = 1          ! reset counter 
      end if   ! nstepread == 1001      
    end if ! lstoreplane 
    

    if (rk3step==1) then
      uminletbc = u0inletbc
      vminletbc = v0inletbc
      wminletbc = w0inletbc
      tminletbc = t0inletbc
    end if

   if (lbuoyancy ) then
!     call blthicknessmo(di_test,utaui,lmoi)
     call blthicknesst(di_test,Uinl,0.99)
!     call dispthicknessmo(displ)  ! needed in top BC
     call dispthicknessexp(displ)
   else
!     call blthickness(di_test,utaui)
     call blthicknesst(di_test,Uinl,0.99)
!     call dispthickness(displ)  ! needed in top BC
     call dispthicknessexp(displ)
   end if
   call blthicknesst(dti_test,Tinl-thls,0.99)
   if ((myid==0) .and. (rk3step==3)) then
     write(6,*) 'Inlet Gen: gamma,lambda=',gamm,lamb
     write(6,*) 'Inlet Gen: Uinl(ke),Tinl(ke)=',Uinl(ke),Tinl(ke)
     write(6,*) 'Inlet Gen: utaui,utaur =',utaui,utaur
     write(6,*) 'Inlet Gen: ttaui,ttaur =',ttaui,ttaur
     write(6,*) 'Inlet Gen: Lmoi,Lmor =',lmoi,lmor
     write(6,*) 'Inlet Gen: deltar, deltai_test', dr, di_test
     write(6,*) 'Inlet Gen: deltatr, deltati_test', dtr, dti_test
     write(6,*) 'Inlet Gen: d*i, d*r=',displ(ib),displ(irecy)
     write(6,*) 'Inlet Gen: thetai,thetar',thetai,thetar
     write(6,*) 'Inlet Gen: thetati,thetatr',thetati,thetatr
     if (luvolflowr ) then
       write(6,*) 'Inlet Gen: mass flux correction factor = ',scalef
!       write(6,*) 'Inlet Gen: mass flux                   = ',totalreadu
       write(6,*) 'Inlet Gen: mass flux                   = ',totaluinl
     end if
   end if
    
    elseif (iinletgen == 2) then
      if (myid==0) then
        write(6,*) 'nstepread=',nstepread
      end if
      u0inletbcold = u0inletbc
      v0inletbcold = v0inletbc
      w0inletbcold = w0inletbc
      t0inletbcold = t0inletbc

  ! determine time step interval in simulation
      rk3coef   = dt / (4. - dble(rk3step))
      if (rk3step==1) then
        deltat = rk3coef
      elseif (rk3step==2) then
        deltat = rk3coef   - (dt/3.)
      elseif (rk3step==3) then
        deltat = rk3coef - (dt/2.)
      end if
  ! determine time step interval in inlet data
      rk3coefin = dtin / (4. - dble(rk3stepin))
      if (rk3stepin==1) then
        dtinrk = rk3coefin
      elseif (rk3stepin==2) then
        dtinrk = rk3coefin - (dtin/3.)
      elseif (rk3stepin==3) then
        dtinrk = rk3coefin - (dtin/2.)
      end if
 
      interval = dtinrk - elapstep
      elapstep = elapstep + deltat
      if (elapstep > dtinrk) then      ! use new value at next time step
        nstepread = nstepread +1
        elapstep = mod(elapstep,dtinrk)
        rk3stepin = mod(rk3stepin,3) + 1
        rk3coefin = dtin / (4. - dble(rk3stepin))
        if (rk3stepin==1) then
          dtinrk = rk3coefin
        elseif (rk3stepin==2) then
          dtinrk = rk3coefin - (dtin/3.)
        elseif (rk3stepin==3) then
          dtinrk = rk3coefin - (dtin/2.)
        end if
        u0inletbc(:,:) = storeu0inletbc(:,:,nstepread)
        v0inletbc(:,:) = storev0inletbc(:,:,nstepread)
        w0inletbc(:,:) = storew0inletbc(:,:,nstepread)
        t0inletbc(:,:) = storet0inletbc(:,:,nstepread)
        if (nstepread == nstore) then
          nfile = nfile + 1
          call readinletfile
          call writerestartfiles
          nstepread = 0
        end if
        interval = dtinrk
        deltat = elapstep
!        write(6,*) 'dtinrk,deltat=', dtinrk,deltat
      end if
      u0inletbc(:,:) = (1. - deltat/interval)*u0inletbc(:,:) + (deltat/interval)*storeu0inletbc(:,:,nstepread+1)
      v0inletbc(:,:) = (1. - deltat/interval)*v0inletbc(:,:) + (deltat/interval)*storev0inletbc(:,:,nstepread+1)
      w0inletbc(:,:) = (1. - deltat/interval)*w0inletbc(:,:) + (deltat/interval)*storew0inletbc(:,:,nstepread+1)
      t0inletbc(:,:) = (1. - deltat/interval)*t0inletbc(:,:) + (deltat/interval)*storet0inletbc(:,:,nstepread+1)




!! massflow correction
      uinletbc2(ib,jb:je,kb:ke) = u0inletbc(jb:je,kb:ke)  ! this is just a dummy variable to give uninletbc the right dimension in slabsum
      urav = 0.
      call slabsum(urav  ,kb,ke,uinletbc2 ,ib,ib,jb,je,kb,ke,ib,ib,jb,je,kb,ke)
      urav = urav / (jge-jgb +1)                    ! average over j-direction

   ! determine bulk velocity of new (interpolated) profile
      do k=kb,ke
        uravdzf(k) = urav(k)*dzf(k)
      end do
      totalu    = sum(uravdzf(kb:ke))/(zh(ke+1)-zh(kb))      ! Area-averaged outflow velocity

   ! rescale the instantaneous profile to keep mass flux constant (tot avoid pressure fluctuations)
      scalef = totalreadu/totalu ! compute factor to scale the velocity profile with
      u0inletbc(:,:) = u0inletbc(:,:)*scalef  ! rescale the velocity profile to have constant mass-flux
!! end of massflow correction of interpolated streamwise velocity 

      if (rk3step==1) then
        uminletbc = u0inletbc
        vminletbc = v0inletbc
        wminletbc = w0inletbc
        tminletbc = t0inletbc
      end if
    end if  ! iinletgen

  end subroutine inletgen
 
  subroutine inletgennotemp
    use modglobal,   only : ib,ie,jb,je,jgb,jge,kb,ke,zf,zh,dzf,dzhi,timee,btime,totavtime,rk3step,dt,numol,iplane,lles,iinletgen,inletav,runavtime,Uinf,lwallfunc,linletRA,totinletav,lstoreplane,nstore,lfixinlet,lfixutauin,luvolflowr
    use modfields,   only : u0,v0,w0,wm,uprof
    use modsave,     only : writerestartfiles
    use modmpi,      only : slabsum,myid
    implicit none
   
    real,dimension(ib:ib,jb:je,kb:ke)   :: uinletbc2   ! dummy variable
    real,dimension(jb:je,kb:ke)   :: uprec             ! velocity fluctuation (up_rec = u0 - Urec)
    real,dimension(jb:je,kb:ke)   :: vprec             ! velocity fluctuation (vp_rec = v0 - 0)
    real,dimension(jb:je,kb:ke+1) :: wprec             ! velocity fluctuation (wp_rec = w0 - Wrec)
    real,dimension(jb:je,kb:ke)   :: upinli,vpinli     ! = gamma * (uprec,v interpolated to zii grid)
    real,dimension(jb:je,kb:ke)   :: upinlo,vpinlo     ! = gamma * (uprec,v interpolated to zoi grid)
    real,dimension(jb:je,kb:ke+1) :: wpinli            ! = gamma * (wprec   interpolated to zii grid)
    real,dimension(jb:je,kb:ke+1) :: wpinlo            ! = gamma * (wprec   interpolated to zoi grid)
    real,dimension(kb:ke)   :: udiff                   ! difference between Uinl and Urec
!    real,dimension(kb:ke)   :: Urecdiff                ! difference between Urec new and old
    real,dimension(kb:ke)   :: urav                    ! j-averaged u-velocity (not time-averaged) 
    real,dimension(kb:ke)   :: uravdzf                 ! j-averaged u-velocity (not time-averaged) times dzf
    real,dimension(kb:ke)   :: uinldzf                 ! j-averaged u-velocity (not time-averaged) times dzf
    real,dimension(kb:ke)   :: Urecdzf                 ! Urec times dzf
    real,dimension(kb:ke+1) :: wrav                    ! j-averaged w-velocity (not time-averaged) 
    real,dimension(kb:ke)   :: Uinli                   ! = gamma * (Urec interpolated to ziif grid points)
    real,dimension(kb:ke+1) :: Winli                   ! = gamma * (Wrec interpolated to ziih grid points)
    real,dimension(kb:ke)   :: Uinlo                   ! = gamma * (Urec interpolated to zioif grid points)
    real,dimension(kb:ke+1) :: Winlo                   ! = gamma * (Wrec interpolated to zoih grid points)
    real,dimension(kb:ke)   :: wfuncf                  ! weight function at full level
    real,dimension(kb:ke+1) :: wfunch                  ! weight function at half level
    real                    :: utaur2,utaui2           ! (utau)^2 at recycle station and inlet
    real                    :: gamm                    ! utaui / utaur
    real                    :: avint,avinti            ! avering interval
    real                    :: alpha,beta              ! factors used in the Weight function
!    real                    :: totalu                  ! total u-velocity at outlet
    real                    :: Urectot                  ! total u-velocity at recycle plane
    real                    :: rk3coef
!    real                    :: di_test                 ! BL thickness as measured from Uinl
    real                    :: utop                    ! j-averaged top velocity
    real                    :: interval
    real                    :: dtinrk                  ! RK time step in inlet data
    real                    :: rk3coefin               ! Cumulative RK time step in inlet data
    real                    :: dr_old
    real                    :: scalef                      ! scale factor to scale instantaneous velocity profile with to get constant mass flux
    real                    :: totaluinl                   ! bulk velocity at the inlet
    integer i,j,k,kk

   if (iinletgen == 1) then

   u0inletbcold = u0inletbc
   v0inletbcold = v0inletbc
   w0inletbcold = w0inletbc
   totaluold    = totalu
   displold     = displ
   ddispdxold   = ddispdx

   ! compute time-average velocities
    rk3coef = dt / (4. - dble(rk3step))
    if (rk3step==1) then
      deltat = rk3coef
    elseif (rk3step==2) then
      deltat = rk3coef - (dt/3.)
    elseif (rk3step==3) then
      deltat = rk3coef - (dt/2.)
    end if

    if (linletRA) then  ! this is a switch to use 'running average'
      avint = totinletav + timee-btime ! runav interval = averaging interval previuous sim  + current elapsed sim time
    else
      avint  = inletav
    end if
    avinti = 1./avint
    uaver=0. 
    do i=ib,ie
      call slabsum(uaver(i,:),kb,ke,u0(i:i,jb:je,kb:ke),i,i,jb,je,kb,ke,i,i,jb,je,kb,ke)
    end do

    wrav=0.
    call slabsum(wrav(kb:ke+1),kb,ke,w0(irecy-1:irecy-1,jb:je,kb:ke+1),irecy-1,irecy-1,jb,je,kb,ke+1,irecy-1,irecy-1,jb,je,kb,ke+1)

    uaver = uaver / (jge-jgb +1)                    ! average over j-direction
    urav = uaver(irecy,:)
    wrav = wrav / (jge-jgb +1)                    ! average over j-direction

    do k=kb,ke
      Urec(k) =  urav(k)*deltat*avinti + (1.-deltat*avinti)*Urec(k)
    end do
    do k=kb,ke+1
      Wrec(k) =  wrav(k)*deltat*avinti + (1.-deltat*avinti)*Wrec(k)
    end do
    do k=kb,ke
    do i=ib,ie
      Utav(i,k) =  uaver(i,k)*deltat*avinti + (1.-deltat*avinti)*Utav(i,k)
    end do
    end do
   


   ! compute velocity fluctuation at recycle station
    do k=kb,ke
      do j=jb,je
        uprec(j,k) = u0(irecy,j,k) - Urec(k)
        vprec(j,k) = v0(irecy-1,j,k)             ! mean v is zero
      end do
    end do

    do k=kb,ke+1
      do j=jb,je
        wprec(j,k) = w0(irecy-1,j,k) - Wrec(k)   ! note that w-velocity is taken at i=irecy-1 !!
      end do
    end do

 
    if (lwallfunc) then
      call wallawinlet(Urec(kb),dzf(kb),numol,utaur2)    ! compute wall shear stress at recycle station
    else
      utaur2 = 2.*numol*Urec(kb)/dzf(kb)
    end if
    utaur = sqrt(abs(utaur2))                          ! compute utau at recycle station
   ! compute momentum thickness at inlet and recycle plane
   
   dr_old = dr
!   call blthickness(dr,utaur)                     ! also needed for thetar
   call blthicknesst(dr,Urec,0.99)                 
!   call momentumthickness(thetai,utaui,di)        ! di is kept fixed
   call momentumthicknessexp(thetai,Uinl)       
!   call momentumthickness(thetar,utaur,dr)
   call momentumthicknessexp(thetar,Urec)       
!   call blthickness(dr,utaur)
  
    if (.not.lfixutauin) then 
      utaui  = utaur* abs(thetar/thetai)**(1./8.)   ! See Lund (1998): 'Similar to Ludwig-Tillmann correlation'
    end if
    gamm = utaui / utaur                          ! Gamma in Lund (1998)

   ! compute inner scaling coordinates
    zirf = utaur*zf / numol                       ! inner scaling zf-coordinate at recycle station 
    zirh = utaur*zh / numol                       ! inner scaling zh-coordinate at recycle station
    ziif = utaui*zf / numol                       ! inner scaling zf-coordinate at inlet station 
    ziih = utaui*zh / numol                       ! inner scaling zh-coordinate at inlet station 
  
   ! compute outer scaling coordinates   
    zorf = zf / dr                                ! outer scaling zf-coor as measured from Uinldinate at recycle station 
    zorh = zh / dr                                ! outer scaling zh-coordinate at recycle station 
    zoif = zf / di                                ! outer scaling zf-coordinate at inlet station  (could be done once, actually..) 
    zoih = zh / di                                ! outer scaling zf-coordinate at inlet station  (could be done once, actually..)

!!!!! Interpolation starts here

 !!! First inner coordinates
   ! determine which elements are needed when recycle velocity profile is interpolated on inlet plane
   ! for u(,v)-components (zf)
    do k=kb,ke
      do kk=kb,ke
        if (zirf(kk) >= ziif(k)) then
          locupif(k)  = kk
          loclowif(k) = kk-1
          exit
        elseif (kk==ke) then
          locupif(k)  = ke+1            ! this means extrapolation!
          loclowif(k) = ke-1            ! waarom niet ke? of wordt dit niet gebruikt?
        end if
      end do
    end do
   ! for w-components (zh)
    do k=kb,ke+1
      do kk=kb,ke+1
        if (zirh(kk) >= ziih(k)) then
          locupih(k)  = kk
          loclowih(k) = kk-1
          exit
        elseif (kk==ke+1) then
          locupih(k)  = ke+2            ! this means extrapolation!
          loclowih(k) = ke
        end if
      end do
    end do
 !!! Finished with inner coordinates

 !!! Do the same trick for outer coordinates
   ! determine which elements are needed when recycle velocity profile is interpolated on inlet plane
   ! for u(,v)-components (zf)
    do k=kb,ke
      do kk=kb,ke
        if (zorf(kk) >= zoif(k)) then
          locupof(k)  = kk
          loclowof(k) = kk-1
          exit
        elseif (kk==ke) then
          locupof(k)  = ke+1            ! this means extrapolation!
          loclowof(k) = ke-1
        end if
      end do
    end do
   ! for w-components (zh)
    do k=kb,ke+1
      do kk=kb,ke+1
        if (zorh(kk) >= zoih(k)) then
          locupoh(k)  = kk
          loclowoh(k) = kk-1
          exit
        elseif (kk==ke+1) then
          locupoh(k)  = ke+2            ! this means extrapolation!
          loclowoh(k) = ke
        end if
      end do
    end do

  !!! Finished with outer coordinates

!!! Now really interpolate
  !!! First inner coordinates
   ! compute Urec on zii grid
    do k=kb,ke
      if (locupif(k) == ke+1) then      ! indicator for extrapolation!
!        Uinli(k) = Urec(ke) + (Urec(ke) - Urec(ke-1)) / (zirf(ke)-zirf(ke-1)) * (ziif(k)-zirf(ke))
        Uinli(k) = Urec(ke)
      elseif (loclowif(k) == kb-1) then ! interprets this as extrapolation to bottom (use u=0 at z+=0)
        Uinli(k) = Urec(kb)/zirf(kb) * ziif(k)
      else                            ! normal interpolation
        Uinli(k) = Urec(loclowif(k)) + (Urec(locupif(k)) - Urec(loclowif(k))) / (zirf(locupif(k)) - zirf(loclowif(k)))* (ziif(k) - zirf(loclowif(k)))
        if ((ziif(k) .gt. zirf(locupif(k))) .or. (ziif(k) .lt. zirf(loclowif(k)))) then
          write(6,*) '!!!Mistake in Interpolation !!!!'
        end if
      end if
    end do

   ! compute Wrec on zii grid
    Winli(kb) = 0.0                      ! corresponds to ground level
    do k=kb+1,ke+1
      if (locupih(k) == ke+2) then     ! indicator for extrapolation!
!        Winli(k) = Wrec(ke+1) + (Wrec(ke+1) - Wrec(ke)) / (zirh(ke+1)-zirh(ke)) * (ziih(k)-zirh(ke+1))
        Winli(k) = Wrec(ke+1)
      else                            ! normal interpolation
        Winli(k) = Wrec(loclowih(k)) + (Wrec(locupih(k)) - Wrec(loclowih(k))) / (zirh(locupih(k)) - zirh(loclowih(k))) * (ziih(k) - zirh(loclowih(k)))
      end if
    end do

   ! compute u- and v- and t-fluctuation on zii grid
    do k=kb,ke
      if (locupif(k) == ke+1) then      ! indicator for extrapolation!
!        upinli(:,k) = uprec(:,ke) + (uprec(:,ke) - uprec(:,ke-1)) / (zirf(ke)-zirf(ke-1)) * (ziif(k)-zirf(ke))
        upinli(:,k) = 0.
        vpinli(:,k) = 0.
      elseif (loclowif(k) == kb-1) then ! interprets this as extrapolation to bottom (use u=0 at z+=0)
        upinli(:,k) = uprec(:,kb)/zirf(kb) * ziif(k)
        vpinli(:,k) = vprec(:,kb)/zirf(kb) * ziif(k)
      else                            ! normal interpolation
        upinli(:,k) = uprec(:,loclowif(k)) + (uprec(:,locupif(k)) - uprec(:,loclowif(k))) / (zirf(locupif(k)) - zirf(loclowif(k))) * (ziif(k) - zirf(loclowif(k)))
        vpinli(:,k) = vprec(:,loclowif(k)) + (vprec(:,locupif(k)) - vprec(:,loclowif(k))) / (zirf(locupif(k)) - zirf(loclowif(k))) * (ziif(k) - zirf(loclowif(k)))
      end if
    end do
   
   ! compute w-fluctuation on zii grid
    do k=kb+1,ke+1
!      if (locupih(k) == ke+1) then      ! indicator for extrapolation!
      if (locupih(k) == ke+2) then      ! indicator for extrapolation!
!        wpinli(:,k) = wprec(:,ke+1) + (wprec(:,ke+1) - wprec(:,ke)) / (zirh(ke+1)-zirh(ke)) * (ziih(k)-zirh(ke+1))
        wpinli(:,k) = 0.
      else                            ! normal interpolation
        wpinli(:,k) = wprec(:,loclowih(k)) + (wprec(:,locupih(k)) - wprec(:,loclowih(k))) / (zirh(locupih(k)) - zirh(loclowih(k))) * (ziih(k) - zirh(loclowih(k)))
      end if
    end do

 !! Finished with interpolating inner variables 
 !! Continue with interpolating outer variables
   ! compute Urec on zoi grid
    do k=kb,ke
      if (locupof(k) == ke+1) then      ! indicator for extrapolation!
!        Uinlo(k) = Urec(ke) + (Urec(ke) - Urec(ke-1)) / (zorf(ke)-zorf(ke-1)) * (zoif(k)-zorf(ke))
!        Uinlo(k) = Urec(ke)
        Uinlo(k) = Uinf
      elseif (loclowof(k) == kb-1) then ! interprets this as extrapolation to bottom (use u=0 at z+=0)
        Uinlo(k) = Urec(kb)/zorf(kb) * zoif(k)
      else                            ! normal interpolation
        Uinlo(k) = Urec(loclowof(k)) + (Urec(locupof(k)) - Urec(loclowof(k))) / (zorf(locupof(k)) - zorf(loclowof(k))) * (zoif(k) - zorf(loclowof(k)))
      end if
    end do
  
   ! compute Wrec on zii grid
    Winlo(kb) = 0.0                      ! corresponds to ground level
    do k=kb+1,ke+1
      if (locupoh(k) == ke+2) then     ! indicator for extrapolation!
!        Winlo(k) = Wrec(ke+1) + (Wrec(ke+1) - Wrec(ke)) / (zorh(ke+1)-zorh(ke)) * (zoih(k)-zorh(ke+1))
        Winlo(k) = Wrec(ke+1) 
      else                            ! normal interpolation
        Winlo(k) = Wrec(loclowoh(k)) + (Wrec(locupoh(k)) - Wrec(loclowoh(k))) / (zorh(locupoh(k)) - zorh(loclowoh(k))) * (zoih(k) - zorh(loclowoh(k)))
      end if
    end do

   ! compute u- and v-fluctuation on zoi grid
    do k=kb,ke
      if (locupof(k) == ke+1) then      ! indicator for extrapolation!
!        upinlo(:,k) = uprec(:,ke) + (uprec(:,ke) - uprec(:,ke-1)) / (zorf(ke)-zorf(ke-1)) * (zoif(k)-zorf(ke))
        upinlo(:,k) = 0.
        vpinlo(:,k) = 0.
      elseif (loclowof(k) == kb-1) then ! interprets this as extrapolation to bottom (use u=0 at z+=0)
        upinlo(:,k) = uprec(:,kb)/zorf(kb) * zoif(k)
        vpinlo(:,k) = vprec(:,kb)/zorf(kb) * zoif(k)
      else                            ! normal interpolation
        upinlo(:,k) = uprec(:,loclowof(k)) + (uprec(:,locupof(k)) - uprec(:,loclowof(k))) / (zorf(locupof(k)) - zorf(loclowof(k))) * (zoif(k) - zorf(loclowof(k)))
        vpinlo(:,k) = vprec(:,loclowof(k)) + (vprec(:,locupof(k)) - vprec(:,loclowof(k))) / (zorf(locupof(k)) - zorf(loclowof(k))) * (zoif(k) - zorf(loclowof(k)))
      end if
    end do

   ! compute w-fluctuation on zoi grid
    do k=kb+1,ke+1
      if (locupoh(k) == ke+2) then      ! indicator for extrapolation!
!        wpinlo(:,k) = wprec(:,ke+1) + (wprec(:,ke+1) - wprec(:,ke)) / (zorh(ke+1)-zorh(ke)) * (zoih(k)-zorh(ke+1))
        wpinlo(:,k) = 0.
      else                            ! normal interpolation
        wpinlo(:,k) = wprec(:,loclowoh(k)) + (wprec(:,locupoh(k)) - wprec(:,loclowoh(k))) / (zorh(locupoh(k)) - zorh(loclowoh(k))) * (zoih(k) - zorh(loclowoh(k)))
      end if
    end do
 !! Finished interpolating outer velocity variables
!!!!! Finished Interpolation! !!!!!


   ! compute rescaled inner variables ! Winli = Winli (interpolation is enough)
    Uinli = gamm* Uinli                   
    upinli = gamm* upinli
    vpinli = gamm* vpinli
    wpinli = gamm* wpinli 
   ! compute rescaled outer variables ! Winlo = Winlo (interpolation is enough)
    Uinlo = gamm* Uinlo  + (1.- gamm)*Uinf
!    Uinlo = gamm* Uinlo  + (1.- gamm)*Urec(ke) 
    upinlo = gamm* upinlo
    vpinlo = gamm* vpinlo
    wpinlo = gamm* wpinlo


   ! Compute weight function (alpha=4, b=0.2)
    alpha = 4.
    beta = 0.2
    wfuncf = 0.5*(1. + tanh( alpha*(zoif-beta)/((1.-2.*beta)*zoif+beta) )/tanh(alpha) ) ! for full level height
    wfunch = 0.5*(1. + tanh( alpha*(zoih-beta)/((1.-2.*beta)*zoih+beta) )/tanh(alpha) ) ! for half level height
    do k=kb,ke
      if (wfuncf(k) .gt. 1.) then  
        wfuncf(k) = 1.
      end if
    end do
    do k=kb,ke+1
      if (wfunch(k) .gt. 1.) then  
        wfunch(k) = 1.
      end if
    end do
 

   ! Compute the velocity components for the inlet BC
    do k=kb,ke  
    do j=jb,je
      u0inletbc(j,k) = (Uinli(k)+ upinli(j,k))*(1.-wfuncf(k)) +  (Uinlo(k) + upinlo(j,k))* wfuncf(k)
      v0inletbc(j,k) =            vpinli(j,k) *(1.-wfuncf(k)) +              vpinlo(j,k) * wfuncf(k) 
    end do
    end do


    do k=kb,ke+1 
    do j=jb,je
      w0inletbc(j,k) = (Winli(k)+ wpinli(j,k))*(1-wfunch(k)) +  (Winlo(k) + wpinlo(j,k))* wfunch(k)
    end do
    end do
    w0inletbc(:,kb) = 0.
    w0inletbc(:,ke+1) = 0.


   ! Compute j-averaged inlet U  (used for compute thetai)
    uinletbc2(ib,jb:je,kb:ke) = u0inletbc(jb:je,kb:ke)  ! this is just a dummy variable to give uninletbc the right dimension in slabsum
    urav = 0.
    call slabsum(urav  ,kb,ke,uinletbc2 ,ib,ib,jb,je,kb,ke,ib,ib,jb,je,kb,ke)
    urav = urav / (jge-jgb +1)                    ! average over j-direction

! determine bulk velocity of new profile
      do k=kb,ke
        uravdzf(k) = urav(k)*dzf(k)
      end do
      totalu = sum(uravdzf(kb:ke))/(zh(ke+1)-zh(kb))      ! Area-averaged outflow velocity

! correct instantaneous inflow velocity for constant mass-flux
      if (luvolflowr ) then
        do k=kb,ke
          uinldzf(k) = Uinl(k)*dzf(k)     
        end do
        totaluinl = sum(uinldzf(kb:ke))/(zh(ke+1)-zh(kb))      ! Area-averaged inflow velocity that should be kept
        scalef = totaluinl/totalu ! compute factor to scale the velocity profile with
        u0inletbc(:,:) = u0inletbc(:,:)*scalef  ! rescale the velocity profile to have constant mass-flux
        urav(:) = urav(:)*scalef                ! also rescale the part that is added to the mean
      end if

   ! Compute j- and time-averaged  inlet U  (used for compute thetai)
    if (.not.lfixinlet) then  ! only update the average inlet profiles when lfixinlet .eqv..false.
      do k=kb,ke
        Uinl(k) =  urav(k)*deltat*avinti + (1.-deltat*avinti)*Uinl(k)
      end do
    end if

! write inletplane to array (and to file after 1000 time steps)   
    if (lstoreplane ) then     
      storeu0inletbc(:,:,nstepread) = u0inletbc(:,:) 
      storev0inletbc(:,:,nstepread) = v0inletbc(:,:)
      storew0inletbc(:,:,nstepread) = w0inletbc(:,:)
      nstepread = nstepread +1
      if (nstepread == nstore+1) then
        nfile = nfile +1       ! next file number
        call writeinletfile    ! write 1000 time steps to file
        call writerestartfiles
        nstepread = 1          ! reset counter 
      end if   ! nstepread == 1001      
    end if ! lstoreplane 
    

    if (rk3step==1) then
      uminletbc = u0inletbc
      vminletbc = v0inletbc
      wminletbc = w0inletbc
    end if

!   call blthickness(di_test,utaui)
   call blthicknesst(di_test,Uinl,0.99)

!   call dispthickness(displ)  ! needed in top BC
   call dispthicknessexp(displ)  ! needed in top BC

   if ((myid==0) .and. (rk3step==3)) then
     write(6,*) 'Inlet Gen: gamma=',gamm
     write(6,*) 'Inlet Gen: Uinl(ke)=',Uinl(ke)
     write(6,*) 'Inlet Gen: utaui,utaur =',utaui,utaur
     write(6,*) 'Inlet Gen: deltar, deltai_test', dr, di_test
     write(6,*) 'Inlet Gen: d*i, d*r=',displ(ib),displ(irecy)
     write(6,*) 'Inlet Gen: thetai,thetar',thetai,thetar
     if (luvolflowr ) then
       write(6,*) 'Inlet Gen: mass flux correction factor = ',scalef
!       write(6,*) 'Inlet Gen: mass flux                   = ',totalreadu
       write(6,*) 'Inlet Gen: mass flux                   = ',totaluinl
     end if
   end if
    
    elseif (iinletgen == 2) then
      if (myid==0) then
        write(6,*) 'nstepread=',nstepread
      end if
      u0inletbcold = u0inletbc
      v0inletbcold = v0inletbc
      w0inletbcold = w0inletbc

  ! determine time step interval in simulation
      rk3coef   = dt / (4. - dble(rk3step))
      if (rk3step==1) then
        deltat = rk3coef
      elseif (rk3step==2) then
        deltat = rk3coef   - (dt/3.)
      elseif (rk3step==3) then
        deltat = rk3coef - (dt/2.)
      end if
  ! determine time step interval in inlet data
      rk3coefin = dtin / (4. - dble(rk3stepin))
      if (rk3stepin==1) then
        dtinrk = rk3coefin
      elseif (rk3stepin==2) then
        dtinrk = rk3coefin - (dtin/3.)
      elseif (rk3stepin==3) then
        dtinrk = rk3coefin - (dtin/2.)
      end if
 
      interval = dtinrk - elapstep
      elapstep = elapstep + deltat
      if (elapstep > dtinrk) then      ! use new value at next time step
        nstepread = nstepread +1
        elapstep = mod(elapstep,dtinrk)
        rk3stepin = mod(rk3stepin,3) + 1
        rk3coefin = dtin / (4. - dble(rk3stepin))
        if (rk3stepin==1) then
          dtinrk = rk3coefin
        elseif (rk3stepin==2) then
          dtinrk = rk3coefin - (dtin/3.)
        elseif (rk3stepin==3) then
          dtinrk = rk3coefin - (dtin/2.)
        end if
        u0inletbc(:,:) = storeu0inletbc(:,:,nstepread)
        v0inletbc(:,:) = storev0inletbc(:,:,nstepread)
        w0inletbc(:,:) = storew0inletbc(:,:,nstepread)
        if (nstepread == nstore) then
          nfile = nfile + 1
          call readinletfile
          call writerestartfiles
          nstepread = 0
        end if
        interval = dtinrk
        deltat = elapstep
!        write(6,*) 'dtinrk,deltat=', dtinrk,deltat
      end if
      u0inletbc(:,:) = (1. - deltat/interval)*u0inletbc(:,:) + (deltat/interval)*storeu0inletbc(:,:,nstepread+1)
      v0inletbc(:,:) = (1. - deltat/interval)*v0inletbc(:,:) + (deltat/interval)*storev0inletbc(:,:,nstepread+1)
      w0inletbc(:,:) = (1. - deltat/interval)*w0inletbc(:,:) + (deltat/interval)*storew0inletbc(:,:,nstepread+1)


!! massflow correction
      uinletbc2(ib,jb:je,kb:ke) = u0inletbc(jb:je,kb:ke)  ! this is just a dummy variable to give uninletbc the right dimension in slabsum
      urav = 0.
      call slabsum(urav  ,kb,ke,uinletbc2 ,ib,ib,jb,je,kb,ke,ib,ib,jb,je,kb,ke)
      urav = urav / (jge-jgb +1)                    ! average over j-direction

   ! determine bulk velocity of new (interpolated) profile
      do k=kb,ke
        uravdzf(k) = urav(k)*dzf(k)
      end do
      totalu    = sum(uravdzf(kb:ke))/(zh(ke+1)-zh(kb))      ! Area-averaged outflow velocity

   ! rescale the instantaneous profile to keep mass flux constant (tot avoid pressure fluctuations)
      scalef = totalreadu/totalu ! compute factor to scale the velocity profile with
      u0inletbc(:,:) = u0inletbc(:,:)*scalef  ! rescale the velocity profile to have constant mass-flux
!! end of massflow correction of interpolated streamwise velocity



      if (rk3step==1) then
        uminletbc = u0inletbc
        vminletbc = v0inletbc
        wminletbc = w0inletbc
      end if
    end if  ! iinletgen

  end subroutine inletgennotemp

  subroutine momentumthicknessexp(output,uinput)

    use modglobal, only : jb,kb,ke,dzf !,Uinf
    use modinletdata, only : ubulk
    use modmpi, only    : myid
    implicit none

       real, dimension(kb:ke), intent(in) :: uinput  !< input velocity
       real, intent(out)                  :: output  !< momentum thickness
       real, dimension(kb:ke)             :: mthick
!       real    :: umax
       integer :: k 
     
!      write(6,*) 'uinletbc(jb,ke),Uinl(ke)=', uinletbc(jb,ke),uinput(ke)
!       umax = maxval(uinput)
       do k=kb,ke
         mthick(k) = ((uinput(k)/uinput(ke)) - (uinput(k)**2. / uinput(ke)**2.) )*dzf(k)
!         mthick(k) = ((uinput(k)/umax) - (uinput(k)**2. / umax**2.) )*dzf(k)

       end do
       output   = sum(mthick)  ! momentum thickness

  end subroutine momentumthicknessexp

  subroutine momentumthickness(output,ustar,blth)

    use modglobal, only : pi,Uinf
!    use modinletdata, only : ubulk
 !   use modmpi, only    : myid
    implicit none

    real, intent(in)                   :: ustar      ! friction velocity
    real, intent(in)                   :: blth       ! boundary layer thickness
    real, intent(out)                  :: output     ! momentum thickness
    real :: B     = 5.0       ! Wake parameter
    real :: C     = 0.5       ! Coles parameter
    real :: kappa = 0.41      ! Von k�r�n constant
    real :: lam               ! = Uinf/ustar
   
    lam = Uinf / ustar
    output = ((1. + C)/(kappa*lam) - (1./(( kappa**2)*(lam**2)))*(2. + 2.*C*(1.852/pi +1.) + (3./2.)*(C**2)))* blth
    
  end subroutine momentumthickness

  subroutine momentumthicknessmo(output,ustar,blth,lmo)

    use modglobal, only : pi,Uinf
!    use modinletdata, only : ubulk
 !   use modmpi, only    : myid
    implicit none

    real, intent(in)                   :: ustar      ! friction velocity
    real, intent(in)                   :: lmo        ! Obukhov length
    real, intent(in)                   :: blth       ! boundary layer thickness
    real, intent(out)                  :: output     ! momentum thickness
    real :: B     = 5.0       ! Wake parameter
    real :: C     = 0.5       ! Coles parameter
    real :: kappa = 0.41      ! Von k�r�n constant
    real :: cmo   = 0.702     ! constant in MO theory (0.135*5.2)
    real :: lam               ! = Uinf/ustar

    lam = Uinf / ustar
    output = (1. + C + 0.5*cmo*blth/lmo)/(kappa*lam) - (1./(( kappa**2)*(lam**2)))*(2. + 2.*C*(1.852/pi +1.) + (3./2.)*(C**2) + (blth-0.25)*2.*cmo/lmo + (1. + 4./pi)*blth*C*cmo/lmo + (1./6.)*((cmo/lmo)**2)*(blth**2) )* blth

  end subroutine momentumthicknessmo

 
  subroutine enthalpythickness(output,tinput,uinput)

    use modglobal, only : jb,kb,ke,dzf !,Uinf
    use modinletdata, only : ubulk
    use modsurfdata, only : thls
    use modmpi, only    : myid
    implicit none

       real, dimension(kb:ke), intent(in) :: tinput  !< input temperature
       real, dimension(kb:ke), intent(in) :: uinput  !< input velocity
       real, intent(out)                  :: output  !< momentum thickness
       real, dimension(kb:ke)             :: ethick 
       real thlsdummy
       integer :: k 
       
       thlsdummy = thls
       if (tinput(ke) == thls) then
         thlsdummy = thls -0.000001
       end if 
       do k=kb,ke
!         ethick(k) = (uinput(k)/uinput(ke)) * ((tinput(k) - tinput(ke)) /(thls - tinput(ke)) )*dzf(k)
         ethick(k) = (uinput(k)/uinput(ke)) * ((tinput(k) - tinput(ke)) /(thlsdummy - tinput(ke)) )*dzf(k)

       end do
       output   = sum(ethick)  ! enthalpy thickness
       if (output==0.) then
         output= 0.000001
       end if

  end subroutine enthalpythickness


 
  subroutine dispthicknessexp(output)
! output is an array of length (ib:ie)) containing displacement thickness values
    use modglobal, only : ib,ie,kb,ke,dzf,xf!,Uinf
    implicit none

       real, dimension(ib:ie),       intent(out) :: output  !< dispacement thickness
       real, dimension(kb:ke)                    :: dthick
!       real    :: umax
       real    :: dispm
       real    :: disp2m
       real    :: xfdispm
       integer :: i,k 
     
     ! write(6,*) 'Uinl(ke)=', uinput(ke)
       do i=ib,ie
!       umax = maxval(Utav(i,:))
       do k=kb,ke
         dthick(k) = (1.- Utav(i,k)/ Utav(i,ke)) *dzf(k)       ! time-averaged, j-averaged velocity
!         dthick(k) = (1.- Utav(i,k)/ umax) *dzf(k)       ! time-averaged, j-averaged velocity

       end do
       output(i)   = sum(dthick)  ! displacement thickness
       end do

       dispm  = sum(output(ib:ie))/(ie-ib+1)      ! mean(displacement)
       disp2m = sum(output(ib:ie)**2.)/(ie-ib+1)   ! mean(displacement^2)
       xfdispm  = sum(xf(ib:ie)*output(ib:ie))/(ie-ib+1)     ! mean(xf*displ) 

       ddispdx = (xfdispm - (xfm*dispm)) / (xf2m -xfm**2.)    ! this is d/dx(delta*)
!       ddispdx = 0.    ! for the test
!       dinl    = dispm - ddispdx*xfm                         ! this is the starting value (delta* = dinl + d/dx(delta)*x)

  end subroutine dispthicknessexp

  subroutine dispthickness(output)
! output is an array of length (ib:ie)) containing displacement thickness values
    use modglobal, only : ib,ie,kb,ke,dzf,xf,Uinf,numol
    implicit none

       real, dimension(ib:ie),       intent(out) :: output  !< dispacement thickness
!       real, dimension(kb:ke)                    :: dthick
       real    :: dispm
       real    :: disp2m
       real    :: xfdispm
       real    :: ustar,blth
       real    :: B     = 5.0       ! Wake parameter
       real    :: C     = 0.5       ! Coles parameter
       real    :: kappa = 0.41      ! Von k�r�n constant
       real    :: lam               ! = Uinf/ustar
       integer :: i
   

       do i=ib,ie
         ustar     =  sqrt(abs(2*numol* Utav(i,kb)/dzf(kb)))              ! average streamwise friction
         lam       = Uinf / ustar
         blth      = (lam*numol/Uinf) * exp( kappa * (lam - B) - 2.*C)      ! See App. Lund et al.
         output(i) = ((1. + C) / (kappa*lam) ) * blth                 
       end do
       dispm  = sum(output(ib:ie))/(ie-ib+1)      ! mean(displacement)
       disp2m = sum(output(ib:ie)**2.)/(ie-ib+1)   ! mean(displacement^2)
       xfdispm  = sum(xf(ib:ie)*output(ib:ie))/(ie-ib+1)     ! mean(xf*displ)

       ddispdx = (xfdispm - (xfm*dispm)) / (xf2m -xfm**2.)    ! this is d/dx(delta*)


  end subroutine dispthickness

  subroutine dispthicknessmo(output)
! output is an array of length (ib:ie)) containing displacement thickness values
    use modglobal, only : ib,ie,kb,ke,dzf,xf,Uinf,numol,grav,prandtlmoli
    use modsurfdata, only : thls
    implicit none

       real, dimension(ib:ie),       intent(out) :: output  !< dispacement thickness
!       real, dimension(kb:ke)                    :: dthick
       real    :: dispm
       real    :: disp2m
       real    :: xfdispm
       real    :: ustar,tstar,blth
       real    :: B     = 5.0       ! Wake parameter
       real    :: C     = 0.5       ! Coles parameter
       real    :: kappa = 0.41      ! Von k�rm�n constant
       real    :: cmo   = 0.702     ! constant in MO theory (0.135*5.2)
       real    :: lam               ! = Uinf/ustar
       real    :: func,dfunc,utaunu,lmo
       integer :: i,n

       blth      = di      ! initial value 
       do i=ib,ie      
         ustar     = sqrt(abs(2.*numol* Utav(i,kb)/dzf(kb)))              ! average streamwise friction at x-location
         tstar     = numol*prandtlmoli* 2.*(Ttav(i,kb)-thls)/(dzf(kb)*ustar)    ! average shear temp. at x-location
         lmo       = (thls*ustar**2)/(kappa*grav*tstar)                  ! obukhov length at this x-location
         if ((lmo >= 10000.) .or. (lmo <= 0.01)) then
           lmo = 1000.
         end if
!         lmo = 0.3  !! TEMPORARY
         utaunu    = ustar / numol
         lam       = Uinf / ustar
         do n=1,10   ! Newton Raphson method to find BL height
!           write(6,*) 'blth,ustar,tstar,Lmo =',blth,ustar,tstar,lmo
           func   = log(blth) + (cmo*blth/lmo) + log(utaunu) - kappa*(lam-B) +2.*C
!           func   = log(blth) + log(utaunu) - kappa*(lam-B) +2.*C
           dfunc  = 1./blth + cmo/lmo
           blth   = blth - (func / dfunc)
           if (blth <= 0.) then
             blth = di
           end if
         end do
         output(i) = ((1. + C + 0.5*cmo*blth/lmo) / (kappa*lam) ) * blth
       end do
       dispm  = sum(output(ib:ie))/(ie-ib+1)      ! mean(displacement)
       disp2m = sum(output(ib:ie)**2.)/(ie-ib+1)   ! mean(displacement^2)
       xfdispm  = sum(xf(ib:ie)*output(ib:ie))/(ie-ib+1)     ! mean(xf*displ)

       ddispdx = (xfdispm - (xfm*dispm)) / (xf2m -xfm**2.)    ! this is d/dx(delta*)


  end subroutine dispthicknessmo



  ! thermal boundary layer thickness
  subroutine blthicknesst(output,uinput,criterion)

    use modglobal, only : kb,ke,zh,zf
    implicit none

       real, dimension(kb:ke), intent(in) :: uinput     !< input velocity
       real, intent(in)                   :: criterion  !< criterion for BL thickness computation (e.g. 0.95 or 0.99)
       real, intent(out)                  :: output     !< BL thickness based on input criterion
!       real, dimension(kb:ke)             :: mthick
       real                               :: ucrit
!       real                               :: umax
       integer :: k

!     umax = maxval(uinput)
     ucrit = uinput(ke)*criterion  ! Velocity at which BL-thickness is reached
!     ucrit = umax*criterion  ! Velocity at which BL-thickness is reached
     do k=kb,ke
       if (uinput(k) .GT. criterion*uinput(ke)) then
         if (k==kb) then
           output = zh(kb)+ (zf(k)-zh(k))/uinput(k)*ucrit ! interpolate z to BL-height
           exit
         else
           output = zf(k-1) + (zf(k)-zf(k-1))/(uinput(k)-uinput(k-1))*(ucrit-uinput(k-1)) !  interpolate z to BL-height
           exit
         end if
       else if (k==ke) then
         output = zf(ke)      ! maximum BL thickness
       end if
     end do
   end subroutine blthicknesst



!  subroutine blthickness(output,uinput,criterion)
!
!    use modglobal, only : kb,ke,zh,zf
!    implicit none
!
!       real, dimension(kb:ke), intent(in) :: uinput     !< input velocity
!       real, intent(in)                   :: criterion  !< criterion for BL thickness computation (e.g. 0.95 or 0.99)
!       real, intent(out)                  :: output     !< BL thickness based on input criterion
!!       real, dimension(kb:ke)             :: mthick
!       real                               :: ucrit
!!       real                               :: umax
!       integer :: k 
!     
!!     umax = maxval(uinput)
!     ucrit = uinput(ke)*criterion  ! Velocity at which BL-thickness is reached
!!     ucrit = umax*criterion  ! Velocity at which BL-thickness is reached
!     do k=kb,ke
!       if (uinput(k) .GT. criterion*uinput(ke)) then
!         if (k==kb) then
!           output = zh(kb)+ (zf(k)-zh(k))/uinput(k)*ucrit ! interpolate z to BL-height
!           exit
!         else  
!           output = zf(k-1) + (zf(k)-zf(k-1))/(uinput(k)-uinput(k-1))*(ucrit-uinput(k-1)) ! interpolate z to BL-height
!           exit
!         end if
!       else if (k==ke) then
!         output = zf(ke)      ! maximum BL thickness
!       end if
!     end do
!   end subroutine blthickness 

  subroutine blthickness(output,ustar)

    use modglobal, only : numol,Uinf
    implicit none

       real, intent(in)                   :: ustar      ! friction velocity
       real, intent(out)                  :: output     !< BL thickness based on law of the wake
!       real, dimension(kb:ke)             :: mthick
!       real                               :: ucrit
!       real                               :: umax
!       integer :: k
       real :: B     = 5.0       ! Wake parameter
       real :: C     = 0.5       ! Coles parameter
       real :: kappa = 0.41      ! Von k�r�n constant
       real :: lam               ! = Uinf/ustar

       lam = Uinf / ustar
       output = (lam*numol/Uinf) * exp( kappa * (lam - B) - 2.*C)      ! See App. Lund et al.

   end subroutine blthickness

  subroutine blthicknessmo(output,ustar,lmo)
! This routine compute the BL thicknes for a buoyancy affected boundary layer:
! Newton-Raphson method is used 
    use modglobal, only : numol,Uinf
    implicit none

       real, intent(in)                   :: ustar      ! friction velocity
       real, intent(in)                   :: lmo        ! Obukhov length
       real, intent(inout)                :: output     !< BL thickness based on law of the wake
!       real, dimension(kb:ke)             :: mthick
!       real                               :: ucrit
!       real                               :: umax
!       integer :: k
       real :: B     = 5.0       ! Wake parameter
       real :: C     = 0.5       ! Coles parameter
       real :: kappa = 0.41      ! Von k�r�n constant
       real :: cmo   = 0.702     ! Constant in MO theory (0.135*5.2)
       real :: lam               ! = Uinf/ustar
       real :: func,dfunc,utaunu
       integer :: n
       
       utaunu = ustar / numol 
       lam = Uinf / ustar
!       write(6,*) 'Initial delta, Lmo =', output,lmo
       do n=1,10
         func   = log(output) + (cmo*output/lmo) + log(utaunu) - kappa*(lam-B) +2.*C
!         func   = log(output) + log(utaunu) - kappa*(lam-B) +2.*C
         dfunc  = 1./output + cmo/lmo
         output = output - (func / dfunc)
         if (output <= 0.) then
           output = di
         end if
       end do
!       write(6,*) 'Computed delta, Lmo =', output,lmo

   end subroutine blthicknessmo



  subroutine wallawinlet(utan,dx,visc,tau)
! this should be the same as wallaw in modboundary!!! This routine is just
! copied to avoid circular dependencies
    implicit none

      real, intent(in)  :: utan,dx,visc
      real, intent(out) :: tau

      real    const1, const2, const3, const4
      real    tausub, taupow
      real    sub, dutan, utankr,utanabs
      real    aaa,bbb
      real    dxi

      parameter(aaa = 8.3)
      parameter(bbb = 0.1428571429)

      dxi = 1./dx
      const1 = 0.5 * (1. - bbb) * aaa ** ((1. + bbb) / (1. - bbb))
      const2 = (1. + bbb) / aaa
      const3 = aaa ** (2. / (1. - bbb))
      const4 = 2. / (1. + bbb)

      utanabs=abs(utan)
      utankr = 0.5 * visc * dxi * const3
      dutan  = utankr - utanabs
      sub    = max (sign(1.,dutan),0.)

      tausub    = 2. * visc * utanabs * dxi
!      taupow3   =   const1 * (visc * dxi)**(1.+bbb) + (const2 * (visc *
!      dxi)**bbb) * utanabs
      taupow    = ( const1 * (visc * dxi)**(1.+bbb) + (const2 * (visc * dxi)**bbb) * utanabs)** const4

!      if (taupow3<=0) then
!        write(6,*) 'taupow3 <=0!!!'
!      end if
      tau = sub * tausub + (1.- sub) * taupow
      tau = sign(tau,utan)  ! give tau the same sign as utan
      return
    end subroutine wallawinlet

  subroutine writeinletfile
    use modglobal, only : jb,je,kb,ke,cexpnr,ifoutput,nstore,ltempeq
    use modmpi,    only : cmyid,myid
!    use modinletdata, only : storeu0inletbc,storev0inletbc,storew0inletbc,nfile

    implicit none
    integer fileid
    integer j,k,n
    character(24) name

      name = 'inlet/inlet_    k   .'
      write (name(13:16)  ,'(i4.4)') nfile
      name(18:20)= cmyid
      name(22:24)= cexpnr

      write(6,*) 'Writing Inlet velocity: ', name
      open  (ifoutput,file=name,form='unformatted',position='append')

      write(ifoutput)  (((storeu0inletbc (j,k,n),j=jb,je),k=kb,ke),  n=1,nstore)
      write(ifoutput)  (((storev0inletbc (j,k,n),j=jb,je),k=kb,ke),  n=1,nstore)
      write(ifoutput)  (((storew0inletbc (j,k,n),j=jb,je),k=kb,ke+1),n=1,nstore)
      close (ifoutput)

      if (ltempeq ) then
        name = 'inlet/itemp_    k   .'
        write (name(13:16)  ,'(i4.4)') nfile
        name(18:20)= cmyid
        name(22:24)= cexpnr

        write(6,*) 'Writing Inlet temperature: ', name
        open  (ifoutput,file=name,form='unformatted',position='append')

        write(ifoutput)  (((storet0inletbc (j,k,n),j=jb,je),k=kb,ke),  n=1,nstore)
        close (ifoutput)
      end if

  end subroutine writeinletfile

  subroutine readinletfile
    use modglobal, only : ib,jb,je,jmax,kb,ke,cexpnr,ifinput,nstore,ltempeq,ntrun,zh,jgb,jge,jh 
    use modmpi,    only : cmyid,myid,nprocs,slabsum,excjs
!    use modinletdata, only : storeu0inletbc,storev0inletbc,storew0inletbc,nfile

    implicit none
    real, dimension(ib:ib,jb:jb+inlfactor*jmax-1,kbin:kein)        :: udummy
    real, dimension(kbin:kein)                                     :: uread
    real, dimension(kbin:kein)                                     :: ureaddzfin
    real, dimension(jb:jb+jtotin-1,kbin:kein,1:nstore)     :: storeu0inold
    real, dimension(jb:jb+jtotin-1,kbin:kein,1:nstore)     :: storev0inold
    real, dimension(jb:jb+jtotin-1,kbin:kein+1,1:nstore)   :: storew0inold
    real, dimension(jb:jb+jtotin-1,kbin:kein,1:nstore)     :: storet0inold
    real, dimension(jb:jb+jtotdum-1,kbin:kein,1:nstore)     :: storeu0indum
    real, dimension(jb:jb+jtotdum-1,kbin:kein,1:nstore)     :: storev0indum
    real, dimension(jb:jb+jtotdum-1,kbin:kein+1,1:nstore)   :: storew0indum
    real, dimension(jb:jb+jtotdum-1,kbin:kein,1:nstore)     :: storet0indum
    real, dimension(jb:je,kbin:kein,1:nstore)     :: storeu0innew
    real, dimension(jb:je,kbin:kein,1:nstore)     :: storev0innew
    real, dimension(jb:je,kbin:kein+1,1:nstore)   :: storew0innew
    real, dimension(jb:je,kbin:kein,1:nstore)     :: storet0innew
    integer filen,filee
    integer fileid
    integer j,k,n,js,jf,jfdum,jsdum
    character(24) name
      jfdum = jbdum-1  ! initial value
      do fileid = filenumstart, filenumstart+(filestoread-1)
        if (filen == -1) then 
          filen = nprocsinl-1                      ! -1 means the last proc (periodic)
        else
          filen = fileid - floor(real(fileid)/real(nprocsinl))*nprocsinl  ! loop over proc's
        end if
!        write(6,*) '!!!!! filen = ', filen
        name = 'inlet/inlet_    k   .'
        write (name(13:16)  ,'(i4.4)') nfile
        write (name(18:20)  ,'(i3.3)') filen
        name(22:24)= cexpnr
        write(6,*) 'Reading Inlet velocity: ', name
        open(unit=ifinput,file=name,form='unformatted')
        read(ifinput)  (((storeu0inold (j,k,n),j=jbin,jein),k=kbin,kein),  n=1,nstore)
        read(ifinput)  (((storev0inold (j,k,n),j=jbin,jein),k=kbin,kein),  n=1,nstore)
        read(ifinput)  (((storew0inold (j,k,n),j=jbin,jein),k=kbin,kein+1),n=1,nstore)
        close (ifinput)

        if (ltempeq ) then
          name = 'inlet/itemp_    k   .'
          write (name(13:16)  ,'(i4.4)') nfile
          write (name(18:20)  ,'(i3.3)') filen
          name(22:24)= cexpnr
          write(6,*) 'Reading Inlet temperature: ', name
          open(unit=ifinput,file=name,form='unformatted')
          read(ifinput)  (((storet0inold (j,k,n),j=jbin,jein),k=kbin,kein),  n=1,nstore)
          close (ifinput)
        end if


        ! determine start and end indices
        if (filen == procinlo) then
          js = jbeg
        else
          js = jbin
        end if
        if (filen == procinup) then
          jf = jend
        else
          jf = jein
        end if
        jsdum = jfdum + 1
        jfdum = jsdum + (jf-js)
!        if (jsdum >= 3) write(6,*) 'myid, jsdum = ',myid, jsdum 
!        if (jfdum >= 3) write(6,*) 'myid, jfdum = ',myid, jfdum 

 !!! put values from original in dummy variable
        storeu0indum(jsdum:jfdum,:,:)    = storeu0inold(js:jf,:,:)      ! s: start  f: final
        storev0indum(jsdum:jfdum,:,:)    = storev0inold(js:jf,:,:)      ! s: start  f: final
        storew0indum(jsdum:jfdum,:,:)    = storew0inold(js:jf,:,:)      ! s: start  f: final
        if (ltempeq ) then
          storet0indum(jsdum:jfdum,:,:)    = storet0inold(js:jf,:,:)      ! s: start  f: final
        end if
      end do  ! loop over original inlet files
  ! now interpolate in y
      call yinterpolate (storeu0indum,storeu0innew,kbin,kein)
      call yinterpolate (storev0indum,storev0innew,kbin,kein)
      call yinterpolate (storew0indum,storew0innew,kbin,kein+1)
      call yinterpolate (storet0indum,storet0innew,kbin,kein)

      if (.not.lzinzsim) then ! interpolate when zin =/ zsim
        call zinterpolate (storeu0innew(:,:,:),storeu0inletbc)   ! interpolate inlet profile to zgrid
        call zinterpolate (storev0innew(:,:,:),storev0inletbc)   ! interpolate inlet profile to zgrid
        call zinterpolatew(storew0innew(:,:,:),storew0inletbc)   ! interpolate inlet profile to zgrid
        if (ltempeq) then
          call zinterpolatet(storet0innew(:,:,:),storet0inletbc)   ! interpolate inlet profile to zgrid
        end if
      else
        storeu0inletbc(:,:,:) = storeu0inold(:,:,:)
        storev0inletbc(:,:,:) = storev0inold(:,:,:)
        storew0inletbc(:,:,:) = storew0inold(:,:,:)
        if (ltempeq) then
          storet0inletbc(:,:,:) = storet0inold(:,:,:)
        end if
      end if

      if (iangle/=0.0) then   ! modify for inflow angle
         do n=1,nstore
         do k=kb,ke
         do j=jb,je
           u0rot(n,j,k) = storeu0inletbc(j,k,n)         ! swap indices in order
           v0rot(n,j,k) = storev0inletbc(j,k,n)         ! to use excjs
         end do
         end do
         end do
         call excjs(u0rot, 1,nstore,jb,je,kb,ke,0,jh)
         call excjs(v0rot, 1,nstore,jb,je,kb,ke,0,jh)
!         write(6,*) 'v0rot(1,je+1,30) = ',v0rot(1,je+1,30)
         do n=1,nstore
         do k=kb,ke
         do j=jb,je
         ! apply horizontal rotation (neglecting the delta_x difference)
           storeu0inletbc(j,k,n) = u0rot(n,j,k)*cos(iangle) - 0.5*sin(iangle)*(v0rot(n,j,k)+v0rot(n,j+1,k))
           storev0inletbc(j,k,n) = v0rot(n,j,k)*cos(iangle) + 0.5*sin(iangle)*(u0rot(n,j,k)+u0rot(n,j-1,k))
         end do
         end do
         end do
      end if ! iangle =/0.0


  end subroutine readinletfile



  subroutine zinterpolate(input,output)
  use modglobal,      only : jb,je,kb,ke,zf,nstore
  implicit none
  real, dimension(jb:je,kbin:kein,1:nstore), intent(in)  :: input
  real, dimension(jb:je,kb:ke,1:nstore), intent(inout) :: output
  integer k

     do k=kb,ke
      if (linuf(k) == kein+1) then      ! indicator for extrapolation!
        output(:,k,:) = input(:,kein,:)
      elseif (linlf(k) == kbin-1) then ! interprets this as extrapolation to bottom (use u=0 at z+=0)
        output(:,k,:) = input(:,kbin,:)/zfin(kbin) * zf(k)
!        output(:,k,:) = input(:,kbin,:)   ! temeperature is not zero at the wall, so line above is wrong
      else                            ! normal interpolation
        output(:,k,:) = input(:,linlf(k),:) + (input(:,linuf(k),:) - input(:,linlf(k),:)) / (zfin(linuf(k)) - zfin(linlf(k))) * (zf(k) - zfin(linlf(k)))
        if ((zf(k) .gt. zfin(linuf(k))) .or. (zf(k) .lt. zfin(linlf(k)))) then
          write(6,*) '!!!Mistake in zinterpolate !!!!'
        end if
      end if
    end do


  end subroutine zinterpolate

  subroutine zinterpolate1d(input,output)
  use modglobal,      only : kb,ke,zf
  implicit none
  real, dimension(kbin:kein), intent(in)  :: input
  real, dimension(kb:ke), intent(inout) :: output
  integer k

     do k=kb,ke
      if (linuf(k) == kein+1) then      ! indicator for extrapolation!
        output(k) = input(kein)
      elseif (linlf(k) == kbin-1) then ! interprets this as extrapolation to bottom (use u=0 at z+=0)
        output(k) = input(kbin)/zfin(kbin) * zf(k)
!        output(:,k,:) = input(:,kbin,:)   ! temeperature is not zero at the wall, so line above is wrong
      else                            ! normal interpolation
        output(k) = input(linlf(k)) + (input(linuf(k)) - input(linlf(k))) / (zfin(linuf(k)) - zfin(linlf(k))) * (zf(k) - zf(linlf(k)))
        if ((zf(k) .gt. zfin(linuf(k))) .or. (zf(k) .lt. zfin(linlf(k)))) then
          write(6,*) '!!!Mistake in zinterpolate1d !!!!'
        end if
      end if
    end do

  end subroutine zinterpolate1d

 subroutine zinterpolate2d(input,output)
  use modglobal,      only : ib,ie,kb,ke,zf,nstore
  implicit none
  real, dimension(ib:ie,kbin:kein), intent(in)  :: input
  real, dimension(ib:ie,kb:ke), intent(inout) :: output
  integer k

     do k=kb,ke
      if (linuf(k) == kein+1) then      ! indicator for extrapolation!
        output(:,k) = input(:,kein)
      elseif (linlf(k) == kbin-1) then ! interprets this as extrapolation to bottom (use u=0 at z+=0)
        output(:,k) = input(:,kbin)/zfin(kbin) * zf(k)
!        output(:,k,:) = input(:,kbin,:)   ! temeperature is not zero at the wall, so line above is wrong
      else                            ! normal interpolation
        output(:,k) = input(:,linlf(k)) + (input(:,linuf(k)) - input(:,linlf(k))) / (zfin(linuf(k)) - zfin(linlf(k))) * (zf(k) - zf(linlf(k)))
        if ((zf(k) .gt. zfin(linuf(k))) .or. (zf(k) .lt. zfin(linlf(k)))) then
          write(6,*) '!!!Mistake in zinterpolate2d !!!!'
        end if
      end if
    end do


  end subroutine zinterpolate2d


  subroutine zinterpolatew(input,output)
  use modglobal,      only : jb,je,kb,ke,zh,nstore
  implicit none
  real, dimension(jb:je,kbin:kein+1,1:nstore), intent(in)  :: input
  real, dimension(jb:je,kb:ke+1,1:nstore), intent(inout) :: output
  integer k


     do k=kb,ke+1
      if (linuh(k) == kein+2) then      ! indicator for extrapolation!
        output(:,k,:) = input(:,kein+1,:)
      elseif (linlh(k) == kbin-1) then ! interprets this as extrapolation to bottom (use u=0 at z+=0)
        output(:,k,:) = input(:,kbin,:)  ! =0
!        output(:,k,:) = input(:,kbin,:)/zhin(kbin) * zh(k)
!        output(:,k,:) = input(:,kbin,:)   ! temeperature is not zero at the wall, so line above is wrong
      else                            ! normal interpolation
        output(:,k,:) = input(:,linlh(k),:) + (input(:,linuh(k),:) - input(:,linlh(k),:)) / (zhin(linuh(k)) - zhin(linlh(k))) * (zh(k) - zhin(linlh(k)))
        if ((zh(k) .gt. zhin(linuh(k))) .or. (zh(k) .lt. zhin(linlh(k)))) then
          write(6,*) '!!!Mistake in zinterpolatew !!!!'
        end if
      end if
    end do


  end subroutine zinterpolatew

  subroutine zinterpolatew1d(input,output)
  use modglobal,      only : kb,ke,zh
  implicit none
  real, dimension(kbin:kein+1), intent(in)  :: input
  real, dimension(kb:ke+1), intent(inout) :: output
  integer k


     do k=kb,ke+1
      if (linuh(k) == kein+2) then      ! indicator for extrapolation!
        output(k) = input(kein+1)
      elseif (linlh(k) == kbin-1) then ! interprets this as extrapolation to bottom (use u=0 at z+=0)
        output(k) = input(kbin) !=0
!        output(k) = input(kbin)/zhin(kbin) * zh(k)
!        output(:,k,:) = input(:,kbin,:)   ! temeperature is not zero at the wall, so line above is wrong
      else                            ! normal interpolation
        output(k) = input(linlh(k)) + (input(linuh(k)) - input(linlh(k))) / (zhin(linuh(k)) - zhin(linlh(k))) * (zh(k) - zh(linlh(k)))
        if ((zh(k) .gt. zhin(linuh(k))) .or. (zh(k) .lt. zhin(linlh(k)))) then
          write(6,*) '!!!Mistake in zinterpolatew1d !!!!'
        end if
      end if
    end do


  end subroutine zinterpolatew1d



  subroutine zinterpolatet(input,output)
  use modglobal,      only : jb,je,kb,ke,zf,nstore
  use modsurfdata,     only : thls
  implicit none
  real, dimension(jb:je,kbin:kein,1:nstore), intent(in)  :: input
  real, dimension(jb:je,kb:ke,1:nstore), intent(inout) :: output
  integer k

     do k=kb,ke
      if (linuf(k) == kein+1) then      ! indicator for extrapolation!
        output(:,k,:) = input(:,kein,:)
      elseif (linlf(k) == kbin-1) then ! interprets this as extrapolation to bottom (use u=0 at z+=0)
         output(:,k,:) = thls + (input(:,kb,:)-thls)/zfin(kbin)*zf(k)
!         output(:,k,:) = (input(:,kb,:)-thls)/zfin(kbin)*zf(k)
!        output(:,k,:) = input(:,kbin,:)/zfin(kbin) * zf(k)
!        output(:,k,:) = input(:,kbin,:)   ! temeperature is not zero at the wall, so line above is wrong
      else                            ! normal interpolation
        output(:,k,:) = input(:,linlf(k),:) + (input(:,linuf(k),:) - input(:,linlf(k),:)) / (zfin(linuf(k)) - zfin(linlf(k))) * (zf(k) - zfin(linlf(k)))
        if ((zf(k) .gt. zfin(linuf(k))) .or. (zf(k) .lt. zfin(linlf(k)))) then
          write(6,*) '!!!Mistake in zinterpolatet !!!!'
        end if
      end if
    end do


  end subroutine zinterpolatet
 
  subroutine zinterpolatet1d(input,output)
  use modglobal,      only : jb,je,kb,ke,zf,nstore
  use modsurfdata,     only : thls
  implicit none
  real, dimension(kbin:kein), intent(in)  :: input
  real, dimension(kb:ke), intent(inout) :: output
  integer k


     do k=kb,ke
      if (linuf(k) == kein+1) then      ! indicator for extrapolation!
        output(k) = input(kein)
      elseif (linlf(k) == kbin-1) then ! interprets this as extrapolation to bottom (use u=0 at z+=0)
         output(k) = (input(kb)-thls)/zfin(kbin)*zf(k)
!        output(:,k,:) = input(:,kbin,:)/zfin(kbin) * zf(k)
!        output(:,k,:) = input(:,kbin,:)   ! temeperature is not zero at the wall, so line above is wrong
      else                            ! normal interpolation
        output(k) = input(linlf(k)) + (input(linuf(k)) - input(linlf(k))) / (zfin(linuf(k)) - zfin(linlf(k))) * (zf(k) - zf(linlf(k)))
        if ((zf(k) .gt. zfin(linuf(k))) .or. (zf(k) .lt. zfin(linlf(k)))) then
          write(6,*) '!!!Mistake in zinterpolatet1d !!!!'
        end if
      end if
    end do


  end subroutine zinterpolatet1d

  subroutine yinterpolate(input,output,ks,kf)
  use modglobal, only : jb,je, nstore
  integer, intent(in) :: ks
  integer, intent(in) :: kf
  real, dimension(jbdum:jedum,ks:kf,1:nstore), intent(in)  :: input
  real, dimension(jb   :je   ,ks:kf,1:nstore), intent(inout) :: output
  integer j

    do j=jb,je
!      if (np==0 .and. yloclowf(j)==)
      output(j,:,:) = input(yloclowf(j),:,:) + (input(ylocupf(j),:,:) - input(yloclowf(j),:,:)) / (yfdum(ylocupf(j)) - yfdum(yloclowf(j))) * (yf(j) - yfdum(yloclowf(j)))
    end do
  end subroutine yinterpolate

  subroutine yinterpolateh(input,output,ks,kf)
  use modglobal, only : jb,je, nstore
  integer, intent(in) :: ks
  integer, intent(in) :: kf
  real, dimension(jbdum:jedum,ks:kf,1:nstore), intent(in)  :: input
  real, dimension(jb   :je   ,ks:kf,1:nstore), intent(inout) :: output
  integer j, jj

    do j=jb,je
!      if (np==0 .and. yloclowf(j)==)
      output(j,:,:) = input(yloclowh(j),:,:) + (input(ylocuph(j),:,:) - input(yloclowh(j),:,:)) / (yhdum(ylocuph(j)) - yhdum(yloclowh(j))) * (yh(j) - yhdum(yloclowh(j)))
    end do
  end subroutine yinterpolateh



  subroutine readzincoord
  use modglobal,   only :  kb,ke,kh,ifinput,zf,zh,ysize,jb,je,dy
  use modmpi,      only :  myid, mpi_integer,comm3d,mpierr,my_real,nprocs  
  implicit none     
  character(72) chmess
  character(20) namezinlet
  character(20) namezinfo
  integer ierr,k,kk,kmaxin,j,jj
  real ysizeproc

      namelist/INFO/nprocsinl,jgtotinl,kmaxin,dtin,wtop,totalreadu

  namezinlet = 'zgrid.inl'
  namezinfo  = 'zgrid.inf'
 
    if (myid==0) then
  

      open(ifinput,file=namezinfo,status='old',iostat=ierr)
        if (ierr /= 0) then
          write(0, *) 'ERROR: zgrid.inf does not exist'
          stop 1
        end if
        read (ifinput,INFO,iostat=ierr)
        if (ierr > 0) then
          write(0, *) 'Problem in zgrid.inf INFO'
          write(0, *) 'iostat error: ', ierr
          stop 1
        endif
        write(6,INFO)
      close(ifinput)
    end if
    kbin = 0
    kein = kmaxin-1
    call MPI_BCAST(nprocsinl,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(jgtotinl ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(kbin     ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(kein     ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(dtin     ,1,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(wtop     ,1,MY_REAL,0,comm3d,mpierr)
    call MPI_BCAST(totalreadu ,1,MY_REAL,0,comm3d,mpierr)
    allocate(zhin(kbin   :kein+1))
    allocate(zfin(kbin   :kein+1))
    allocate(dzfin(kbin-1:kein+1))
    allocate(dzhin(kbin  :kein+1))

    if (myid==0) then
      write(6,*) 'loading ',namezinlet
      open (ifinput,file=namezinlet)
      read(ifinput,'(a72)') chmess
      read(ifinput,'(a72)') chmess

      do k=kbin,kein
        read(ifinput,*) zfin(k)
      end do
      close(ifinput)

      zhin(kbin) = 0.0
      do k=kbin,kein
        zhin(k+1) = zhin(k) + 2.0*(zfin(k)-zhin(k))
      end do
      zfin(kein+kh)  = zfin(kein)+ 2.0*(zhin(kein+kh)-zfin(kein))

      do  k=kbin,kein
        dzfin(k) = zhin(k+1) - zhin(k)
      end do
      dzfin(kein+1) = dzfin(kein)
      dzfin(kbin-1) = dzfin(kbin)

      dzhin(kbin) = 2*zfin(kbin)
      do k=kbin+1,kein+kh
        dzhin(k) = zfin(k) - zfin(k-1)
      end do

    ! check if the inlet mesh and the simulation mesh differ
    do k=kb,min(ke,kein)
      if(abs(zfin(k)-zf(k)) > 1e-7) then
        lzinzsim = .false.
      end if
    enddo

    end if ! myid==0

  ! MPI broadcast kmax elements from zf
    call MPI_BCAST( zfin,kein-kbin+2,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST( zhin,kein-kbin+2,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dzfin,kein-kbin+3,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(dzhin,kein-kbin+2,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(lzinzsim,1,MPI_INTEGER      ,0,comm3d,mpierr)

    
    if (.not.lzinzsim) then 
      if (myid==0) then
        write(6,*) 'zgrid.inl does not equal zgrid.inp: Inlet will be interpolated in z'
      end if
!      allocate(linlf(kbin:kein)) 
!      allocate(linuf(kbin:kein)) 
!      allocate(linlh(kbin:kein+1)) 
!      allocate(linuh(kbin:kein+1)) 
      allocate(linlf(kb:ke)) 
      allocate(linuf(kb:ke)) 
      allocate(linlh(kb:ke+1)) 
      allocate(linuh(kb:ke+1)) 
   ! zf
      do k=kb,ke
        do kk=kbin,kein
          if (zfin(kk) >= zf(k)) then
            linuf(k)  = kk
            linlf(k)  = kk-1
            exit
          elseif (kk==kein) then
            linuf(k)  = kein+1            ! this means extrapolation!
            linlf(k)  = kein-1            ! waarom niet ke? of wordt dit niet gebruikt?
          end if
        end do
      end do
   ! for w-components (zh)
      do k=kb,ke+1
        do kk=kbin,kein+1
          if (zhin(kk) >= zh(k)) then
            linuh(k)  = kk
            linlh(k) = kk-1
            exit
          elseif (kk==kein+1) then
            linuh(k)  = kein+2            ! this means extrapolation!
            linlh(k)  = kein
          end if
        end do
      end do

     else  ! lzinzsim  -> grids are equal
       if (myid==0) then 
         write(6,*) 'zgrid.inl equals zgrid.inp: Inlet will not be interpolated in z'
       end if
     end if
    
! Now prepare everything for interpolation in y-direction
     jgbin = 1
     jgein = jgbin + jgtotinl -1
     jtotin = jgtotinl / nprocsinl
     jbin = 1
     jein = 1 + jtotin -1
     ysizeproc = ysize/nprocs
     dyin = ysize / jgtotinl
     jbdum = 1
     jtotdum = ceiling(ysizeproc/real(dyin))+1                    ! dummy indices
     jedum = jbdum + jtotdum-1
     
!     allocate(yf   (jb    :je)) 
     allocate(yf   (jb    :je+1)) 
     allocate(yh   (jb    :je+1)) 
     allocate(yfin (jgbin :jgein+1)) 
     allocate(yhin (jgbin-1 :jgein+1)) 
!     allocate(yfdum(jbdum :jedum)) 
     allocate(yfdum(jbdum :jedum+1)) 
     allocate(yhdum(jbdum :jedum+1)) 
     allocate(yloclowf(jb:je+1))
     allocate(ylocupf (jb:je+1))
     allocate(yloclowh(jb:je+1))
     allocate(ylocuph (jb:je+1))

    ! make global y-grid (equidistant) for inlet data
      
     do j = jgbin-1,jgein+1
       yhin(j) = (j-jgbin)*dyin
     end do
     do j = jgbin,jgein+1
       yfin(j) = yhin(j) + 0.5*dyin
     end do
    ! make new y-grid (equidistant)
     do j = jb,je+1
       yh(j) = myid*(ysize/nprocs) + (j-jb)*dy
       yf(j) = yh(j) + 0.5*dy
     end do

   ! check which original cells are needed for interpolation
      do j=jgein+1,jgbin,-1
        if (yhin(j)<= yh(jb)) then
          if (yfin(j)<= yf(jb)) then  
            procinlo = floor(real(j-jgbin)/real(jtotin))      ! this is the first cell to consider
            filenumstart = procinlo
            jgbeg = j
            jbeg = j-(procinlo*jtotin)
!            jend = jbeg+jtotdum-1
            jj = j+jtotdum-1
!            procinup = floor((j-jgbin)/real(jtotin))     
!            procinup = floor((j-jgbin+1)/real(jtotin))
            procinup = floor(real(jj-jgbin)/real(jtotin))
            filenumend = procinup
            jend = jj-(procinup*jtotin)
            procinup = procinup-floor(real(procinup)/real(nprocsinl))*nprocsinl  ! continue on first procinl again
          else
            if (j == jgbin) then
              jgbeg = j-1
              jbeg = jein
              procinlo = nprocsinl-1
              filenumstart = -1
              jj = j+jtotdum-2
              procinup = floor(real(jj-jgbin)/real(jtotin))
              filenumend = procinup
              jend = jj-(procinup*jtotin)
              procinup = procinup-floor(real(procinup)/real(nprocsinl))*nprocsinl  !continue on first procinl again
            else  
              procinlo = floor(real(j-jgbin-1)/real(jtotin))    ! One cell lower is needed
              filenumstart = procinlo
              jgbeg = j-1
              jbeg = j-(procinlo*jtotin)-1         
              jj = j+jtotdum-2
              procinup = floor(real(jj-jgbin)/real(jtotin))
              filenumend = procinup
              jend = jj-(procinup*jtotin)
              procinup = procinup-floor(real(procinup)/real(nprocsinl))*nprocsinl  ! continue on first procinl again
            end if ! j=jgbin
          end if
          exit
        end if
      end do

      write(6,*) '!! myid,procinlo,jbeg,procinup,jend,jgbeg = ',myid,procinlo,jbeg,procinup,jend,jgbeg


   ! make dummy y-grid (equidistant)
      do j = jbdum,jedum+1
        yhdum(j) = yhin(jgbeg) + (j-jbdum+1)*dyin
        yfdum(j) = yhdum(j) + 0.5*dyin
      end do
!      if (procoldup /= procoldlo) then
!        write(6,*) '!!! Start-cell and end-cell are not in the same file!!!'
!      end if
      filestoread = filenumend - filenumstart + 1 ! no. of files to be read
!      write(6,*) '!! procinlo,procinup = ',procinlo,procinup
!      write(6,*) '!! jbin,jein,jbeg,jend= ',jbin,jein,jbeg,jend

   ! for components defined on yf
    do j=jb,je
      do jj=jbdum+1,jedum+1
        if (yfdum(jj) >= yf(j)) then
          ylocupf(j)  = jj
          yloclowf(j) = jj-1
          exit
        end if
      end do
    end do
   ! for components defined on yh
    do j=jb,je+1
      do jj=jbdum+1,jedum+1
        if (yhdum(jj) >= yh(j)) then
          ylocuph(j)  = jj
          yloclowh(j) = jj-1
          exit
        end if
      end do
    end do

  end subroutine readzincoord


  subroutine exitinlet
  use modglobal,      only : iinletgen,lstoreplane,ltempeq

  if (iinletgen==1) then
    deallocate(Uinl,Winl,Urec,Wrec,u0inletbc,v0inletbc,w0inletbc,zirf,ziif,ziih,zirh,zorf,zoif,zorh,zoih,loclowif,locupif,loclowih,locupih,loclowof,locupof,loclowoh,locupoh,uminletbc,vminletbc,wminletbc,u0inletbcold,v0inletbcold,w0inletbcold,Utav,upupavinl,vpvpavinl,wpwpavinl,upwpavinl,thlpthlpavinl,thlpupavinl,thlpwpavinl)
    if (ltempeq ) then
       deallocate(t0inletbc,tminletbc,t0inletbcold,loclowot,locupot,zotr,zoti,Tinl,Trec)
    end if
    if (lstoreplane ) then 
      deallocate(storeu0inletbc,storev0inletbc,storew0inletbc)
       if (ltempeq ) then
         deallocate(storet0inletbc)
       end if
    end if 
  else if (iinletgen == 2) then
    deallocate(storeu0inletbc,storev0inletbc,storew0inletbc,u0inletbc,v0inletbc,w0inletbc,uminletbc,vminletbc,wminletbc,u0inletbcold,v0inletbcold,w0inletbcold)
    if (ltempeq ) then
       deallocate(t0inletbc,tminletbc,t0inletbcold,storet0inletbc)
    end if
  end if

  end subroutine exitinlet

end module
