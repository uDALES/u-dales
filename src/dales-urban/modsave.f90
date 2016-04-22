!> \file modsave.f90
!! Writes restart and data files.
!> 
!! modsave.f90 writes the restart and data files
!!  \author Jasper Tomas, June 4th 2015
!!  \todo documentation
!!  \par Revision list
!

module modsave


implicit none
! private
! public :: writerestartfiles, writedatafiles
save

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writerestartfiles
    use modsurfdata,only: ustar,thlflux,qtflux,svflux,dudz,dvdz,dthldz,dqtdz,ps,thls,qts,thvs,oblav,&
                          isurf
    use modfields, only : u0,v0,w0,thl0,qt0,ql0,ql0h,e120,dthvdz,presf,presh,sv0,mindist,wall,&
                          uav,vav,wav,uuav,vvav,wwav,uvav,uwav,vwav,thlav,thl2av,qtav,qlav,ql2av,qt2av,svav,sv2av,momthick,&
                          friction,displthick,pres0,viscratioav,thluav,thlvav,thlwav,qtuav,qtvav,qtwav,qluav,qlvav,qlwav,svuav,svvav,svwav,&
                          upupav,vpvpav,wpwpav,thlpthlpav,qlpqlpav,qtpqtpav,svpsvpav,upvpav,upwpav,vpwpav,thlpupav,thlpvpav,&
                          thlpwpav,qlpupav,qlpvpav,qlpwpav,qtpupav,qtpvpav,qtpwpav,svpupav,svpvpav,svpwpav,presav,&
                          uusgsav,vvsgsav,wwsgsav,uwsgsav,thlusgsav,thlwsgsav,qlusgsav,qlwsgsav,qtusgsav,qtwsgsav,svusgsav,svwsgsav,tkesgsav,&
                          strain2av,disssgsav,t_vav,tvmx,tvmy,tvmz,tsgsmx1,tsgsmx2,tsgsmy1,tsgsmy2,tsgsmz1,&
                          tsgsmz2,t_sgsav,nusgsav,tpm,t_pav,ttmx,ttmy,ttmz,t_tav,p_bav,d_sgsav,p_tav,tkeadv
    use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh,dsv,trestart,tnextrestart,dt_lim,timee,btime,xh,&
                          cexpnr,ntimee,rk3step,ifoutput,nsv,timeleft,dtheta,dqt,dt,cu,cv,ntrun,totavtime,&
                          linletgen,timee,runavtime,inletav,totinletav,linletRA,ltec3d,ltempeq,lmoist,jgb,jge,&
                          dzf,dzfi,dzhi,dxf,dxfi,dyi,dxhi,nstore,ltec2d,numol,dy2i,grav,libm,jmax,nblocks
    use modgenstat,only : lcanopy,ltkebudget,timecompl
    use modmpi,    only : cmyid,myid,slabsum,excjs
    use modsubgriddata, only : ekm
    use modibmdata,   only  : ibmxforcevol,block
    use modfielddump, only   : tec3d,tec2d
    use modinletdata, only   : Urec,Wrec,Uinl,Utav,QLinl,QTinl,QLrec,QTrec,QTtav,QLtav,Ttav,upupavinl,vpvpavinl,wpwpavinl,upwpavinl,&
                               thlpthlpavinl,thlpupavinl,thlpwpavinl,qlpqlpavinl,qlpupavinl,qlpwpavinl,qtpqtpavinl,qtpupavinl,qtpwpavinl,Tinl,Trec,nstepread

    implicit none
    real, dimension(kb:ke+1)            :: wpwpavinlh    ! dummy variable
    real, dimension(kb:ke+1)            :: upwpavinlh    ! dummy variable
    real, dimension(kb:ke+1)            :: thlpwpavinlh  ! dummy variable
    real, dimension(kb:ke+1)            :: qtpwpavinlh  ! dummy variable
    real, dimension(kb:ke+1)            :: qlpwpavinlh  ! dummy variable
    real, dimension(ib:ie,jb:je,kb:ke)  :: dissresav     ! average resolved dissipation: 2*nu*<Sij'*Sij'> = 2*nu*( <Sij*Sij> - <Sij>*<Sij> ) 
    real, dimension(ib:ie,jb:je,kb:ke)  :: disssgsfl     ! average subgrid visc. * average rate of strain squared : 2*<nu_t>*<Sij>*<Sij> 
    real, dimension(ib:ie,jb:je,kb:ke)  :: tke           ! tke = 0.5*<ui'ui'>
    real, dimension(ib:ie,jb:je,kb:ke)  :: mke           ! = <ui>d/dxj(<ui><uj>) + <ui>d/dxj(<ui'uj'>) = <ui>d/dxj(<ui*uj>)
    real, dimension(ib:ie+1,jb  :je,  kb:ke)    :: dummyx           
    real, dimension(ib:ie,  jb-1:je+1,kb:ke)    :: dummyy           
    real, dimension(ib:ie,  jb  :je,  kb:ke+1)  :: dummyz           
    real newtotavtime,timeinterval,strainav2
    logical :: lexitnow = .false.
    integer imin,ihour
    integer i,j,k,n,im,ip,jm,jp,jpp,km,kp,kpp,il,iu,jl,ju,kl,ku
    character(21) name,name2,name3,name4,linkname

    if (timee == 0) return
!    if (rk3step /=3) return
    if (linletgen==2 .and. nstepread==nstore) then                ! This overrules the need for rk3step to be 3 in case of reading inletfiles
      write(6,*) 'Writing restartfiles after reading in new inletfiles'
    else
      if (rk3step /=3) return   ! Normal check
    end if

    name = 'exit_now.'//cexpnr
    inquire(file=trim(name), EXIST=lexitnow)

    if (timee>=tnextrestart .or. lexitnow .or. nstepread == nstore+1) then
      tnextrestart = tnextrestart+trestart
! compute <u'u'> = <uu> - <u><u>
! first all 'squares': <x'x'> =(ib-1:ie+1,jb-1,je+1,kb-1,ke+1)
      upupav     = uuav - uav**2
      vpvpav     = vvav - vav**2
      wpwpav     = wwav - wav**2
      thlpthlpav = thl2av - thlav**2
      svpsvpav   = sv2av - svav**2
      qlpqlpav   = ql2av - qlav**2
      qtpqtpav   = qt2av - qtav**2
! Compute stresses and fluxes
! uw: ib:ie+1 jb:je kb:ke+1
    do k=kb,ke+1
      km = k-1
      do j=jb,je
        do i=ib,ie+1
          im = i-1
          upwpav(i,j,k) = uwav(i,j,k) - & 
                       0.25*(uav(i,j,km)*dzf(k) + uav(i,j,k)*dzf(km))*dzhi(k) * &     ! interpolate uav to edges
                            (wav(im,j,k)*dxf(i) + wav(i,j,k)*dxf(im))*dxhi(i)         ! interpolate wav to edges
        end do
      end do
    end do

! uv: ib:ie+1 jb:je kb:ke+1
    do k=kb,ke
      do j=jb,je+1
        jm = j-1
        do i=ib,ie+1
          im = i-1
          upvpav(i,j,k) = uvav(i,j,k) - &
                       0.25*(uav(i,jm,k)  + uav(i,j,k)) * &                           ! interpolate uav to edges
                            (vav(im,j,k)*dxf(i) + vav(i,j,k)*dxf(im))*dxhi(i)         ! interpolate vav to edges
        end do
      end do
    end do


! vw: ib:ie jb:je+1 kb:ke+1
    do k=kb,ke+1
      km = k-1
      do j=jb,je+1
        jm = j-1
        do i=ib,ie
          vpwpav(i,j,k) = vwav(i,j,k) - &
                       0.25*(wav(i,jm,k)  + wav(i,j,k)) * &                           ! interpolate wav to edges
                            (vav(i,j,km)*dzf(k) + vav(i,j,k)*dzf(km))*dzhi(k)         ! interpolate vav to edges
        end do
      end do
    end do

! thlu and svu: ib:ie+1 jb:je kb:ke  (located on u-faces)
    do k=kb,ke
      do j=jb,je
        do i=ib,ie+1
          im = i-1
          thlpupav(i,j,k) = thluav(i,j,k) - &
                        0.5* uav(i,j,k)  * &                                              ! no interpolation
                            (thlav(im,j,k)*dxf(i) + thlav(i,j,k)*dxf(im))*dxhi(i)         ! interpolate thlav to u-faces
          qlpupav(i,j,k) = qluav(i,j,k) - & 
                        0.5* qlav(i,j,k)  * & ! no interpolation
                            ( qlav(im,j,k)*dxf(i) + qlav(i,j,k)*dxf(im))*dxhi(i)         ! interpolate thlav to u-faces
          qtpupav(i,j,k) = qtuav(i,j,k) - & 
                        0.5* qtav(i,j,k)  * & ! no interpolation
                            ( qtav(im,j,k)*dxf(i) + qtav(i,j,k)*dxf(im))*dxhi(i)         ! interpolate thlav to u-faces
          do n=1,nsv
            svpupav(i,j,k,n) = svuav(i,j,k,n) - &
                          0.5* uav(i,j,k)  * &                                              ! no interpolation
                              (svav(im,j,k,n)*dxf(i) + svav(i,j,k,n)*dxf(im))*dxhi(i)        ! interpolate svav to u-faces
          end do
        end do
      end do
    end do
! thlv and svv: ib:ie jb:je+1 kb:ke  (located on v-faces)
    do k=kb,ke
      do j=jb,je+1
        jm = j-1
        do i=ib,ie
          thlpvpav(i,j,k) = thlvav(i,j,k) - &
                       0.5* vav(i,j,k) * &                                                ! no interpolation
                           (thlav(i,jm,k) + thlav(i,j,k))                                 ! interpolate thlav to v-faces
          
          qlpvpav(i,j,k) = qlvav(i,j,k) - &
                       0.5* vav(i,j,k) * &                                                ! no interpolation
                           (qlav(i,jm,k) + qlav(i,j,k))                                 ! interpolate thlav to v-faces
          
          qtpvpav(i,j,k) = qtvav(i,j,k) - &
                       0.5* vav(i,j,k) * &                                                ! no interpolation
                           (qtav(i,jm,k) + qtav(i,j,k))                                 ! interpolate thlav to v-faces
          

          do n=1,nsv
            svpvpav(i,j,k,n) = svvav(i,j,k,n) - &
                         0.5* vav(i,j,k) * &                                              ! no interpolation
                             (svav(i,jm,k,n) + svav(i,j,k,n))                             ! interpolate svav to v-faces
          end do
        end do
      end do
    end do
! thlw and svw: ib:ie jb:je kb:ke+1  (located on w-faces)
    do k=kb,ke+1
      km = k-1
      do j=jb,je
        do i=ib,ie
          thlpwpav(i,j,k) = thlwav(i,j,k) - &
                       0.5 * wav(i,j,k) * &                                               ! no interpolation
                            (thlav(i,j,km)*dzf(k) + thlav(i,j,k)*dzf(km))*dzhi(k)         ! interpolate thl to w-faces
          
          qlpwpav(i,j,k) = thlwav(i,j,k) - &
                       0.5 * wav(i,j,k) * &                                               ! no interpolation
                            (qlav(i,j,km)*dzf(k) + qlav(i,j,k)*dzf(km))*dzhi(k)         ! interpolate thl to w-faces
          
          qtpwpav(i,j,k) = qtwav(i,j,k) - &
                       0.5 * wav(i,j,k) * &                                               ! no interpolation
                            (qtav(i,j,km)*dzf(k) + qtav(i,j,k)*dzf(km))*dzhi(k)         ! interpolate thl to w-faces
          
          do n=1,nsv
            svpwpav(i,j,k,n) = svwav(i,j,k,n) - &
                         0.5 * wav(i,j,k) * &                                              ! no interpolation
                              (svav(i,j,km,n)*dzf(k) + svav(i,j,k,n)*dzf(km))*dzhi(k)      ! interpolate svav to w-faces
          end do
        end do
      end do
    end do

 
    do k=kb,ke
      km = k-1
      kp = k+1
      do j=jb,je
        jm = j-1
        jp = j+1
        do i=ib,ie
          im = i-1
          ip = i+1

          strainav2 =  ( &
            ((uav(ip,j,k)-uav(i,j,k))    *dxfi(i)     )**2    + &
            ((vav(i,jp,k)-vav(i,j,k))    *dyi         )**2    + &
            ((wav(i,j,kp)-wav(i,j,k))    *dzfi(k)     )**2    )

          strainav2 = strainav2 + 0.125 * ( &
            ((wav(i,j,kp)-wav(im,j,kp))   *dxhi(i)     + &
            (uav(i,j,kp)-uav(i,j,k))      *dzhi(kp)  )**2    + &
            ((wav(i,j,k)-wav(im,j,k))     *dxhi(i)     + &
            (uav(i,j,k)-uav(i,j,km))      *dzhi(k)   )**2    + &
            ((wav(ip,j,k)-wav(i,j,k))     *dxhi(ip)     + &
            (uav(ip,j,k)-uav(ip,j,km))    *dzhi(k)   )**2    + &
            ((wav(ip,j,kp)-wav(i,j,kp))   *dxhi(ip)     + &
            (uav(ip,j,kp)-uav(ip,j,k))    *dzhi(kp)  )**2    )

          strainav2 = strainav2 + 0.125 * ( &
            ((uav(i,jp,k)-uav(i,j,k))     *dyi     + &
            (vav(i,jp,k)-vav(im,jp,k))    *dxhi(i)        )**2    + &
            ((uav(i,j,k)-uav(i,jm,k))     *dyi     + &
            (vav(i,j,k)-vav(im,j,k))      *dxhi(i)        )**2    + &
            ((uav(ip,j,k)-uav(ip,jm,k))   *dyi     + &
            (vav(ip,j,k)-vav(i,j,k))      *dxhi(ip)       )**2    + &
            ((uav(ip,jp,k)-uav(ip,j,k))   *dyi     + &
            (vav(ip,jp,k)-vav(i,jp,k))    *dxhi(ip)       )**2    )

          strainav2 = strainav2 + 0.125 * ( &
            ((vav(i,j,kp)-vav(i,j,k))    *dzhi(kp) + &
            (wav(i,j,kp)-wav(i,jm,kp))   *dyi        )**2    + &
            ((vav(i,j,k)-vav(i,j,km))    *dzhi(k)+ &
            (wav(i,j,k)-wav(i,jm,k))     *dyi        )**2    + &
            ((vav(i,jp,k)-vav(i,jp,km))  *dzhi(k)+ &
            (wav(i,jp,k)-wav(i,j,k))     *dyi        )**2    + &
            ((vav(i,jp,kp)-vav(i,jp,k))  *dzhi(kp) + &
            (wav(i,jp,kp)-wav(i,j,kp))   *dyi        )**2    )

          dissresav(i,j,k) = 2.*numol  *(strain2av(i,j,k) - strainav2)  ! resolved dissipation 
          disssgsfl(i,j,k) = 2.*nusgsav(i,j,k)*strainav2                ! = 2*<nu_sgs>*<sij>*<sij>
          tke(i,j,k)       = 0.5*(0.5*(upupav(ip,j,k)+upupav(i,j,k)) + &
                                  0.5*(vpvpav(i,jp,k)+vpvpav(i,j,k)) + &
                                  0.5*(wpwpav(i,j,kp)+wpwpav(i,j,k)))
        end do
      end do
    end do

    if (ltkebudget ==.true.) then


      ! Tvav = (Tvm - <ui>*d/dxj(<Sij>)  ) + 2*nu*<Sij'Sij'>
      ! Tvm = Tvmx + Tvmy + Tvmz -> therefore: subtraction, then interpolation, then addition of 2*nu*<Sij'Sij'>
      do k=kb,ke
      km = k-1
      kp = k+1
      do j=jb,je
      jp = j+1
      jm = j-1
      do i=ib,ie
        im = i-1
        ip = i+1
!        t_vav(i,j,k) =  0.5*( (tvmx(i,j,k) - (                      &
        dummyx(i,j,k) =  (                      &
                  ( numol  * (uav(i+1,j,k)-uav(i,j,k))*dxfi(i) &
                   -numol * (uav(i,j,k)-uav(i-1,j,k))*dxfi(i-1) ) * 2. * dxhi(i) &
                  + &
                  ( numol * ( (uav(i,jp,k)-uav(i,j,k))   *dyi &
                            +(vav(i,jp,k)-vav(i-1,jp,k))*dxhi(i)) &
                   -numol * ( (uav(i,j,k)-uav(i,jm,k))   *dyi &
                            +(vav(i,j,k)-vav(i-1,j,k))  *dxhi(i)) &
                                       ) * dyi &
                  + &
                  ( numol * ( (uav(i,j,kp)-uav(i,j,k))   *dzhi(kp) &
                            +(wav(i,j,kp)-wav(i-1,j,kp))*dxhi(i)) &
                   -numol * ( (uav(i,j,k)-uav(i,j,km))   *dzhi(k) &
                            +(wav(i,j,k)-wav(i-1,j,k))  *dxhi(i)) &
                                       ) *dzfi(k) )             
      end do
      end do
      end do


 ! y-direction
      do k=kb,ke
      km = k-1
      kp = k+1
      do j=jb,je
      jp = j+1
      jm = j-1
      do i=ib,ie
        im = i-1
        ip = i+1
        dummyy(i,j,k) =  (                        &
                  ( numol * ( (vav(i+1,j,k)-vav(i,j,k))   *dxhi(i+1) &
                        +(uav(i+1,j,k)-uav(i+1,jm,k))*dyi) &
                -numol * ( (vav(i,j,k)-vav(i-1,j,k))   *dxhi(i) &
                        +(uav(i,j,k)-uav(i,jm,k))    *dyi) &
                           ) * dxfi(i) &        ! = d/dx( Km*(dv/dx + du/dy) )
                + &
              (numol * (vav(i,jp,k)-vav(i,j,k)) &
              -numol * (vav(i,j,k)-vav(i,jm,k))  ) * 2. * dy2i &        ! = d/dy( 2*Km*(dv/dy) )
                + &
              ( numol * ( (vav(i,j,kp)-vav(i,j,k))    *dzhi(kp) &
                        +(wav(i,j,kp)-wav(i,jm,kp))  *dyi) &
                -numol * ( (vav(i,j,k)-vav(i,j,km))    *dzhi(k) &
                        +(wav(i,j,k)-wav(i,jm,k))    *dyi)   &
                           ) * dzfi(k) )                    ! = d/dz( Km*(dv/dz + dw/dy) )
      end do
      end do
      end do
      
 ! z-direction
      do k=kb,ke
      km = k-1
      kp = k+1
      do j=jb,je
      jp = j+1
      jm = j-1
      do i=ib,ie
        im = i-1
        ip = i+1
             dummyz(i,j,k) = (                        &
                   ( numol * ( (wav(i+1,j,k)-wav(i,j,k))    *dxhi(i+1) &
                            +(uav(i+1,j,k)-uav(i+1,j,km)) *dzhi(k) ) &
                    -numol * ( (wav(i,j,k)-wav(i-1,j,k))    *dxhi(i) &
                            +(uav(i,j,k)-uav(i,j,km))     *dzhi(k) ) &
                             )*dxfi(i) &
                + &
                  ( numol * ( (wav(i,jp,k)-wav(i,j,k))     *dyi &
                            +(vav(i,jp,k)-vav(i,jp,km))   *dzhi(k) ) &
                    -numol * ( (wav(i,j,k)-wav(i,jm,k))     *dyi &
                            +(vav(i,j,k)-vav(i,j,km))     *dzhi(k) ) &
                             )*dyi &
                + &
                  ( numol * (wav(i,j,kp)-wav(i,j,k)) *dzfi(k) &
                   -numol * (wav(i,j,k)-wav(i,j,km)) *dzfi(km) ) * 2. &
                                               * dzhi(k))
      end do
      end do
      end do          

  ! BC's 
   tvmx   (ie+1,:,:) =  tvmx   (ie,:,:) 
   tsgsmx1(ie+1,:,:) =  tsgsmx1(ie,:,:)
   tsgsmx2(ie+1,:,:) =  tsgsmx2(ie,:,:)
   dummyx (ie+1,:,:) =  dummyx (ie,:,:)
   ttmx   (ie+1,:,:) =  ttmx   (ie,:,:) 
!   call excjs( tvmy   , ib,ie,jb-1,je+1,kb,ke,0,1)   ! jb-1 is not used
!   call excjs( tsgsmy1, ib,ie,jb-1,je+1,kb,ke,0,1)   ! jb-1 is not used
!   call excjs( tsgsmy2, ib,ie,jb-1,je+1,kb,ke,0,1)   ! jb-1 is not used
!   call excjs( dummyy,  ib,ie,jb-1,je+1,kb,ke,0,1)   ! jb-1 is not used
!   call excjs( ttmy   , ib,ie,jb-1,je+1,kb,ke,0,1)   ! jb-1 is not used
   call excjs( tvmy   , ib,ie,jb,je,kb,ke,0,1)   ! jb-1 is not used
   call excjs( tsgsmy1, ib,ie,jb,je,kb,ke,0,1)   ! jb-1 is not used
   call excjs( tsgsmy2, ib,ie,jb,je,kb,ke,0,1)   ! jb-1 is not used
   call excjs( dummyy,  ib,ie,jb,je,kb,ke,0,1)   ! jb-1 is not used
   call excjs( ttmy   , ib,ie,jb,je,kb,ke,0,1)   ! jb-1 is not used
   tvmz   (:,:,ke+1) =  tvmz   (:,:,ke) 
   tsgsmz1(:,:,ke+1) =  tsgsmz1(:,:,ke)
   tsgsmz2(:,:,ke+1) =  tsgsmz2(:,:,ke)
   dummyz (:,:,ke+1) =  dummyz(:,:,ke)
   ttmz   (:,:,ke+1) =  ttmz   (:,:,ke) 



      do k=kb,ke
      km = k-1
      kp = k+1
      do j=jb,je
      jp = j+1
      jm = j-1
      do i=ib,ie
        im = i-1
        ip = i+1
 
 ! Total viscous transport
        t_vav(i,j,k) =  0.5*( (tvmx(i, j,k) - dummyx(i,j,k) *uav(i, j,k))  +        &
                              (tvmx(ip,j,k) - dummyx(ip,j,k)*uav(ip,j,k)))         &
                      + 0.5*( (tvmy(i,j, k) - dummyy(i,j,k) *vav(i,j, k))  +        &
                              (tvmy(i,jp,k) - dummyy(i,jp,k)*vav(i,jp,k)))         &
                      + 0.5*( (tvmz(i,j,k ) - dummyz(i,j,k) *wav(i,j,k ))  +        &
                              (tvmz(i,j,kp) - dummyz(i,j,kp)*wav(i,j,kp)))         &
                      + dissresav(i,j,k)         ! d/dxj(2*nu*<ui'Sij'>) = <u_i*d/dxj(2*nu*Sij')> +2*nu*<Sij'Sij'>

! Now the same for subgrid stress        
! <d/dxj(2*u_i'*nu_t*Sij)'> = <u_i'*d/dxj(2*nu_t*Sij)'> + <(2*nu_t*Sij)'*Sij'>
!                           = <u_i*d/dxj(2*nu_t*Sij)> - <u_i>*<d/dxj(2*nu_t*Sij)>   + <2*nu_t*Sij*Sij> - <(2*nu_t*Sij)>*<Sij>
!                           = <u_i*d/dxj(2*nu_t*Sij)> - <u_i>*<d/dxj(2*nu_t*Sij)>   + <2*nu_t*Sij*Sij> - 2*<nu_t>*<Sij>*<Sij> - 2*<nu_t'*Sij'>*<Sij>


        t_sgsav(i,j,k) =  0.5*( (tsgsmx1(i,j,k) -  uav(i,j,k) *tsgsmx2(i,j,k))   +   &
                                (tsgsmx1(ip,j,k) - uav(ip,j,k)*tsgsmx2(ip,j,k)))     & 
                        +                                                            & ! = <2*nu_t*SijSij> - <2*nu_t*Sij>*<Sij>
                          0.5*( (tsgsmy1(i,j,k) -  vav(i,j,k) *tsgsmy2(i,j,k))   +   &
                                (tsgsmy1(i,jp,k) - vav(i,jp,k)*tsgsmy2(i,jp,k)))     & 
                        +                                                            & ! = <2*nu_t*SijSij> - <2*nu_t*Sij>*<Sij>        
                          0.5*( (tsgsmz1(i,j,k) -  vav(i,j,k) *tsgsmz2(i,j,k))   +   &
                                (tsgsmz1(i,j,kp) - vav(i,j,kp)*tsgsmz2(i,j,kp)))     &  
                        + disssgsav(i,j,k) - disssgsfl(i,j,k)                          ! -2*<nu_t'Sij'>*<Sij>  should still be added!

        d_sgsav(i,j,k)= - disssgsav(i,j,k) + disssgsfl(i,j,k)                          ! +2*<nu_t'Sij'>*<Sij>  should still be added! (is compensated with above)

! - <uj'*dp'/dxj> = - <uj*dp/dxj> + <uj>*d<p>/dxj
        t_pav(i,j,k)   = tpm(i,j,k) +                                                    &
                            0.5*(uav(i,j,k)  *(presav(i,j,k)-presav(i-1,j,k))*dxhi(i) + &
                                 uav(i+1,j,k)*(presav(i+1,j,k)-presav(i,j,k))*dxhi(i+1))&
                                   +                                                    &
                            0.5*(vav(i,j,k)  *(presav(i,j,k)-presav(i,j-1,k))*dyi +     &
                                 vav(i,j+1,k)*(presav(i,j+1,k)-presav(i,j,k))*dyi)      &
                                   +                                                    &
                            0.5*(wav(i,j,k)  *(presav(i,j,k)-presav(i,j,k-1))*dzhi(k) + &
                                 wav(i,j,k+1)*(presav(i,j,k+1)-presav(i,j,k))*dzhi(k+1))
! - d/dxj(<0.5*ui'ui'uj'>) = -<uj'd/dxj(<0.5*ui'ui'>) + <ui'uj'><Sij>
!                          = -<uj*d/dxj(0.5*ui'ui')> + <uj>*d/dxj(<0.5*ui'ui'> + <ui'uj'><Sij>) 


!        ttav(i,j,k)   = ttm(i,j,k) - 
 
! <advection term N.S. times ui> = MKE + A - Pshear - Tt
! Tt = -<ui'd/dxj(ui'uj')> = -<d/dxj(0.5*ui'ui'uj')> = A + MKE - Pshear - Total

        !Pshear =Ptav = -<ui'uj'>d/dxj(<Sij>) = -<ui'uj'>d<ui>/dxj
  
        p_tav(i,j,k)   =  - (                                                                      &
                         0.5 *(upupav(i,j,k)+upupav(ip,j,k))* (uav(ip,j,k)-uav(i,j,k))*dxfi(i)  + & ! <u'u'>*d<u>/dx     
                         0.25*(upvpav(i,j,k)  *(uav(i, j, k)-uav(i, jm,k) )*dyi +                 &
                               upvpav(i,jp,k) *(uav(i, jp,k)-uav(i, j, k) )*dyi +                 &
                               upvpav(ip,j,k) *(uav(ip,j, k)-uav(ip,jm,k) )*dyi +                 &
                               upvpav(ip,jp,k)*(uav(ip,jp,k)-uav(ip,j, k) )*dyi)                + & ! <u'v'>*d<u>/dy
                         0.25*(upwpav(i, j,k ) *(uav(i, j,k )-uav(i, j,km))*dzhi(k) +             &
                               upwpav(i, j,kp) *(uav(i, j,kp)-uav(i, j,k ))*dzhi(kp) +            &
                               upwpav(ip,j,k ) *(uav(ip,j,k )-uav(ip,j,km))*dzhi(k) +             & 
                               upwpav(ip,j,kp) *(uav(ip,j,kp)-uav(ip,j,k ))*dzhi(kp))           + & ! <u'w'>*d<u>/dz
                         0.25*(upvpav(i, j, k) *(vav(i, j, k)-vav(im,j ,k))*dxhi(i) +             & 
                               upvpav(ip,j, k) *(vav(ip,j, k)-vav(i, j ,k))*dxhi(ip) +            &
                               upvpav(i, jp,k) *(vav(i ,jp,k)-vav(im,jp,k))*dxhi(i) +             &
                               upvpav(ip,jp,k) *(vav(ip,jp,k)-vav(i, jp,k))*dxhi(ip))           + & ! <u'v'>*d<v>/dx
                         0.5 *(vpvpav(i,j,k)+vpvpav(i,jp,k))* (vav(i,jp,k)-vav(i,j,k))*dyi      + & ! <v'v'>*d<v>/dy
                         0.25*(vpwpav(i,j ,k ) *(vav(i,j ,k )-vav(i,j ,km))*dzhi(k) +             &
                               vpwpav(i,j ,kp) *(vav(i,j ,kp)-vav(i,j ,k ))*dzhi(kp) +            &
                               vpwpav(i,jp,k ) *(vav(i,jp,k )-vav(i,jp,km))*dzhi(k) +             & 
                               vpwpav(i,jp,kp) *(vav(i,jp,kp)-vav(i,jp,k ))*dzhi(kp))           + & ! <v'w'>*d<v>/dz
                         0.25*(upwpav(i, j, k) *(wav(i, j,k )-wav(im,j,k ))*dxhi(i) +             &
                               upwpav(ip,j, k) *(wav(ip,j,k )-wav(i, j,k ))*dxhi(ip) +            &
                               upwpav(i, j,kp) *(wav(i ,j,kp)-wav(im,j,kp))*dxhi(i) +             &
                               upwpav(ip,j,kp) *(wav(ip,j,kp)-wav(i, j,kp))*dxhi(ip))           + & ! <u'w'>*d<w>/dx
                         0.25*(vpwpav(i,j,k)  *(wav(i,j, k )-wav(i,jm,k ) )*dyi +                 &
                               vpwpav(i,jp,k) *(wav(i,jp,k )-wav(i,j, k ) )*dyi +                 &
                               vpwpav(ip,j,k) *(wav(i,j, kp)-wav(i,jm,kp) )*dyi +                 &
                               vpwpav(ip,jp,k)*(wav(i,jp,kp)-wav(i,j, kp) )*dyi)                + & ! <v'w'>*d<w>/dy
                         0.5 *(wpwpav(i,j,k)+wpwpav(i,j,kp))* (wav(i,j,kp)-wav(i,j,k))*dzfi(k) )    ! <w'w'>*d<w>/dz 
                           
        mke(i,j,k)    = 0.5*(uav(ip,j,k)+uav(i,j,k))*(uuav(ip,j,k)-uuav(i,j,k))*dxfi(i)  +        & ! <u>*d<uu>/dx                      
                        0.5*(uav(i, j,k)*(uvav(i ,jp,k)-uvav(i ,j,k))*dyi  +                      & ! <u>*d<uv>/dy                      
                             uav(ip,j,k)*(uvav(ip,jp,k)-uvav(ip,j,k))*dyi)               +        &
                        0.5*(uav(i, j,k)*(uwav(i ,j,kp)-uwav(i ,j,k))*dzfi(k)  +                  & ! <u>*d<uw>/dz                      
                             uav(ip,j,k)*(uwav(ip,j,kp)-uwav(ip,j,k))*dzfi(k))           +        &
                        0.5*(vav(i,j, k)*(uvav(ip,j ,k)-uvav(i,j ,k))*dxfi(i) +                   & ! <v>*d<uv>/dx
                             vav(i,jp,k)*(uvav(ip,jp,k)-uvav(i,jp,k))*dxfi(i))           +        & 
                        0.5*(vav(i,jp,k)+vav(i,j,k))*(vvav(i,jp,k)-vvav(i,j,k))*dyi      +        & ! <v>*d<vv>/dy
                        0.5*(vav(i,j ,k)*(vwav(i,j ,kp)-vwav(i,j ,k))*dzfi(k)  +                  & ! <v>*d<vw>/dz                      
                             vav(i,jp,k)*(vwav(i,jp,kp)-vwav(i,jp,k))*dzfi(k))           +        &
                        0.5*(wav(i,j,k )*(uwav(ip,j,k )-uwav(i,j,k ))*dxfi(i) +                   & ! <w>*d<uw>/dx
                             wav(i,j,kp)*(uwav(ip,j,kp)-uwav(i,j,kp))*dxfi(i))           +        &
                        0.5*(wav(i,j,k )*(vwav(i,jp,k )-vwav(i,j,k ))*dyi  +                      & ! <w>*d<vw>/dy                      
                             wav(i,j,kp)*(vwav(i,jp,kp)-vwav(i,j,kp))*dyi)               +        &
                        0.5*(wav(i,j,kp)+wav(i,j,k))*(wwav(i,j,kp)-wwav(i,j,k))*dzfi(k)             ! <w>*d<ww>/dz  

        tkeadv(i,j,k) = 0.5*(uav(i, j,k)*(tke(i, j,k)-tke(im,j,k))*dxhi(i) + &                      ! <u>*de/dx
                             uav(ip,j,k)*(tke(ip,j,k)-tke(i ,j,k))*dxhi(ip))         +            & !
                        0.5*(vav(i, j,k)*(tke(i,j ,k)-tke(i,jm,k))*dyi     + &                      ! <v>*de/dy
                             vav(i,jp,k)*(tke(i,jp,k)-tke(i,j ,k))*dyi)              +            & 
                        0.5*(wav(i,j,k )*(tke(i,j,k )-tke(i,j,km))*dzhi(k) + &                      ! <w>*de/dz
                             wav(i,j,kp)*(tke(i,j,kp)-tke(i,j,k ))*dzhi(kp))
! <advection term N.S. times ui> = MKE + A - Pshear - Tt
! Tt = -<ui'd/dxj(ui'uj')> = -<d/dxj(0.5*ui'ui'uj')> = A      +    MKE   -   Pshear  -   Total
!                                                    = tkeadv +    mke   -   p_tav   -   ttm
!        t_tav(i,j,k)   = tkeadv(i,j,k) + mke(i,j,k) - p_tav(i,j,k) - ttm(i,j,k)                 
        t_tav(i,j,k)   = tkeadv(i,j,k) + mke(i,j,k) - p_tav(i,j,k)  &
                          - 0.5*(ttmx(i,j,k) + ttmx(ip,j,k))        &
                          - 0.5*(ttmy(i,j,k) + ttmy(i,jp,k))        &
                          - 0.5*(ttmz(i,j,k) + ttmz(i,j,kp))                
         
        p_bav(i,j,k)   = (grav/thls)*0.5*(thlpwpav(i,j,k)+thlpwpav(i,j,kp))       
                      
      end do
      end do
      end do
  
    end if ! ltkebudget

    ! set all TKE transport terms to zero inside blocks
    if (libm == .true.) then
      do n=1,nblocks
        il = block(n,1)
        iu = block(n,2)
        kl = block(n,5)
        ku = block(n,6)
        jl = block(n,3)-myid*jmax
        ju = block(n,4)-myid*jmax
        if (ju < jb .or. jl > je) then
          cycle
        else
          if (ju > je) ju=je
          if (jl < jb) jl=jb
          tkeadv   (il:iu,jl:ju,kl:ku) = 0.
          t_tav    (il:iu,jl:ju,kl:ku) = 0.
          t_pav    (il:iu,jl:ju,kl:ku) = 0.
          t_vav    (il:iu,jl:ju,kl:ku) = 0.
          t_sgsav  (il:iu,jl:ju,kl:ku) = 0.
          p_bav    (il:iu,jl:ju,kl:ku) = 0.
          p_tav    (il:iu,jl:ju,kl:ku) = 0.
          dissresav(il:iu,jl:ju,kl:ku) = 0.
          d_sgsav  (il:iu,jl:ju,kl:ku) = 0.
          tkesgsav (il:iu,jl:ju,kl:ku) = 0.
          disssgsav(il:iu,jl:ju,kl:ku) = 0.
        end if
      end do
    end if ! libm


    
      name = 'initd        _   .'
      write (name(6:13)  ,'(i8.8)') ntrun
      name(15:17)= cmyid
      name(19:21)= cexpnr
      open  (ifoutput,file=name,form='unformatted',status='replace')

      write(ifoutput)  (((mindist(i,j,k),i=ib,ie  ),j=jb,je      ),k=kb,ke   )
      write(ifoutput)  ((((wall(i,j,k,n),i=ib,ie  ),j=jb,je      ),k=kb,ke   ),n=1,5)
      write(ifoutput)  (((u0 (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((v0 (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((w0    (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((pres0 (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((thl0  (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((e120  (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((ekm   (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((qt0   (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((ql0   (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((ql0h  (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  timee,  dt
      
      
      close (ifoutput)

      if (nsv>0) then
        name  = 'inits        _   .'
        write (name(6:13) ,'(i8.8)') ntrun
        name(15:17) = cmyid
        name(19:21) = cexpnr
        open  (ifoutput,file=name,form='unformatted')
        write(ifoutput) ((((sv0(i,j,k,n),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh),n=1,nsv)
        write(ifoutput)  timee

        close (ifoutput)
      end if

! Write mean values
!      newtotavtime  = totavtime + timee - btime   ! completed time in this simulation + totavtime (from previous sim)
      newtotavtime  = timecompl                    ! completed time in this simulation + totavtime (from previous sim)
      name2 = 'means   .'
      name2( 6:8 ) = cmyid
      name2(10:12) = cexpnr
      open  (ifoutput,file=name2,form='unformatted',status='replace')
      write(ifoutput)  newtotavtime, nsv 
      write(ifoutput)  (((uav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      write(ifoutput)  (((vav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      write(ifoutput)  (((wav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      write(ifoutput)  (((thlav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      write(ifoutput)  (((qlav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      write(ifoutput)  (((qtav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      write(ifoutput)  (((presav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      write(ifoutput) ((((svav(i,j,k,n),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   ),n=1,nsv)
      write(ifoutput)  (((viscratioav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      write(ifoutput)  (((uuav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      write(ifoutput)  (((vvav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      write(ifoutput)  (((wwav(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      write(ifoutput)  (((thl2av(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      write(ifoutput)  (((ql2av(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      write(ifoutput)  (((qt2av(i,j,k),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   )
      write(ifoutput) ((((sv2av(i,j,k,n),i=ib-ih,ie+ih  ),j=jb-jh,je+jh      ),k=kb-kh,ke+kh   ),n=1,nsv) 
      write(ifoutput)  (((uvav(i,j,k),i=ib,ie+ih  ),j=jb,je+jh      ),k=kb,ke   )
      write(ifoutput)  (((uwav(i,j,k),i=ib,ie+ih  ),j=jb,je         ),k=kb,ke+kh   )
      write(ifoutput)  (((vwav(i,j,k),i=ib,ie     ),j=jb,je+jh      ),k=kb,ke+kh   )
      write(ifoutput)  (((thluav(i,j,k),i=ib,ie     ),j=jb,je         ),k=kb,ke      )
      write(ifoutput)  (((thlvav(i,j,k),i=ib,ie     ),j=jb,je+jh      ),k=kb,ke      )
      write(ifoutput)  (((thlwav(i,j,k),i=ib,ie     ),j=jb,je         ),k=kb,ke+kh   )
      write(ifoutput)  (((qluav(i,j,k),i=ib,ie     ),j=jb,je         ),k=kb,ke      )
      write(ifoutput)  (((qlvav(i,j,k),i=ib,ie     ),j=jb,je+jh      ),k=kb,ke      )
      write(ifoutput)  (((qlwav(i,j,k),i=ib,ie     ),j=jb,je         ),k=kb,ke+kh   )
      write(ifoutput)  (((qtuav(i,j,k),i=ib,ie     ),j=jb,je         ),k=kb,ke      )
      write(ifoutput)  (((qtvav(i,j,k),i=ib,ie     ),j=jb,je+jh      ),k=kb,ke      )
      write(ifoutput)  (((qtwav(i,j,k),i=ib,ie     ),j=jb,je         ),k=kb,ke+kh   )
      write(ifoutput) ((((svuav(i,j,k,n),i=ib,ie     ),j=jb,je         ),k=kb,ke      ),n=1,nsv)
      write(ifoutput) ((((svvav(i,j,k,n),i=ib,ie     ),j=jb,je+jh      ),k=kb,ke      ),n=1,nsv)
      write(ifoutput) ((((svwav(i,j,k,n),i=ib,ie     ),j=jb,je         ),k=kb,ke+kh   ),n=1,nsv)

      close (ifoutput)

!     Write <x'y'> to file  (Don't read them for restart. This can all be derived from means.xxx).
      name2 = 'reyns   .'
      name2( 6:8 ) = cmyid
      name2(10:12) = cexpnr
      open  (ifoutput,file=name2,form='unformatted',status='replace')
      write(ifoutput)  newtotavtime, nsv
      write(ifoutput)  (((upupav(i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
      write(ifoutput)  (((vpvpav(i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
      write(ifoutput)  (((wpwpav(i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
      write(ifoutput)  (((thlpthlpav(i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
      write(ifoutput)  (((qlpqlpav(i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
      write(ifoutput)  (((qtpqtpav(i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
      write(ifoutput) ((((svpsvpav(i,j,k,n),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh),n=1,nsv)
      write(ifoutput) ((((svpsvpav(i,j,k,n),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh),n=1,nsv)
      write(ifoutput) ((((svpsvpav(i,j,k,n),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh),n=1,nsv)
      write(ifoutput)  (((upvpav(i,j,k),i=ib,ie+ih  ),j=jb,je+jh      ),k=kb,ke   )
      write(ifoutput)  (((upwpav(i,j,k),i=ib,ie+ih  ),j=jb,je         ),k=kb,ke+kh   )
      write(ifoutput)  (((vpwpav(i,j,k),i=ib,ie     ),j=jb,je+jh      ),k=kb,ke+kh   )
      write(ifoutput)  (((thlpupav(i,j,k),i=ib,ie+ih  ),j=jb,je         ),k=kb,ke      )
      write(ifoutput)  (((thlpvpav(i,j,k),i=ib,ie     ),j=jb,je+jh      ),k=kb,ke      )
      write(ifoutput)  (((thlpwpav(i,j,k),i=ib,ie     ),j=jb,je         ),k=kb,ke+kh   )
      write(ifoutput)  (((qlpupav(i,j,k),i=ib,ie+ih  ),j=jb,je         ),k=kb,ke      )
      write(ifoutput)  (((qlpvpav(i,j,k),i=ib,ie     ),j=jb,je+jh      ),k=kb,ke      )
      write(ifoutput)  (((qlpwpav(i,j,k),i=ib,ie     ),j=jb,je         ),k=kb,ke+kh   )
      write(ifoutput)  (((qtpupav(i,j,k),i=ib,ie+ih  ),j=jb,je         ),k=kb,ke      )
      write(ifoutput)  (((qtpvpav(i,j,k),i=ib,ie     ),j=jb,je+jh      ),k=kb,ke      )
      write(ifoutput)  (((qtpwpav(i,j,k),i=ib,ie     ),j=jb,je         ),k=kb,ke+kh   )
      write(ifoutput) ((((svpupav(i,j,k,n),i=ib,ie+ih  ),j=jb,je         ),k=kb,ke      ),n=1,nsv)
      write(ifoutput) ((((svpvpav(i,j,k,n),i=ib,ie     ),j=jb,je+jh      ),k=kb,ke      ),n=1,nsv)
      write(ifoutput) ((((svpwpav(i,j,k,n),i=ib,ie     ),j=jb,je         ),k=kb,ke+kh   ),n=1,nsv)
      close (ifoutput)


!     Write <x'y'>_SGS to file.
      name2 = 'SGS__   .'
      name2( 6:8 ) = cmyid
      name2(10:12) = cexpnr
      open  (ifoutput,file=name2,form='unformatted',status='replace')
      write(ifoutput)  newtotavtime, nsv
      write(ifoutput)  (((uusgsav(i,j,k)   ,i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
      write(ifoutput)  (((vvsgsav(i,j,k)   ,i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
      write(ifoutput)  (((wwsgsav(i,j,k)   ,i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
      write(ifoutput)  (((uwsgsav(i,j,k)   ,i=ib   ,ie+ih),j=jb   ,je   ),k=kb   ,ke+kh)
      write(ifoutput)  (((dissresav(i,j,k) ,i=ib   ,ie   ),j=jb   ,je   ),k=kb   ,ke   )   ! mean resolved dissipation
      write(ifoutput)  (((tkesgsav(i,j,k)  ,i=ib   ,ie   ),j=jb   ,je   ),k=kb   ,ke   )   ! mean estimated subgrid energy
      write(ifoutput)  (((disssgsav(i,j,k) ,i=ib   ,ie   ),j=jb   ,je   ),k=kb   ,ke   )   ! mean subgrid dissipation
      write(ifoutput)  (((strain2av(i,j,k) ,i=ib   ,ie   ),j=jb   ,je   ),k=kb   ,ke   )   ! <SijSij> (NOT <Sij><Sij> !!) (average over time)
      write(ifoutput)  (((nusgsav(i,j,k)   ,i=ib   ,ie   ),j=jb   ,je   ),k=kb   ,ke   )   ! <nu_sgs> (average over time)
      write(ifoutput)  (((thlusgsav(i,j,k) ,i=ib,ie+ih  ),j=jb,je         ),k=kb,ke      )
      write(ifoutput)  (((thlwsgsav(i,j,k) ,i=ib,ie     ),j=jb,je         ),k=kb,ke+kh   )
      write(ifoutput)  (((qlusgsav(i,j,k) ,i=ib,ie+ih  ),j=jb,je         ),k=kb,ke      )
      write(ifoutput)  (((qlwsgsav(i,j,k) ,i=ib,ie     ),j=jb,je         ),k=kb,ke+kh   )
      write(ifoutput)  (((qtusgsav(i,j,k) ,i=ib,ie+ih  ),j=jb,je         ),k=kb,ke      )
      write(ifoutput)  (((qtwsgsav(i,j,k) ,i=ib,ie     ),j=jb,je         ),k=kb,ke+kh   )
      write(ifoutput) ((((svusgsav(i,j,k,n),i=ib,ie+ih  ),j=jb,je         ),k=kb,ke      ),n=1,nsv)
      write(ifoutput) ((((svwsgsav(i,j,k,n),i=ib,ie     ),j=jb,je         ),k=kb,ke+kh   ),n=1,nsv)
      close (ifoutput)

      if (ltkebudget==.true.) then
        name2 = 'TKE__   .'
        name2( 6:8 ) = cmyid
        name2(10:12) = cexpnr
        open  (ifoutput,file=name2,form='unformatted',status='replace')
        write(ifoutput)  newtotavtime, nsv
        write(ifoutput)  (((tkeadv    (i,j,k)   ,i=ib,ie),j=jb,je),k=kb,ke)
        write(ifoutput)  (((t_pav     (i,j,k)   ,i=ib,ie),j=jb,je),k=kb,ke)
        write(ifoutput)  (((t_tav     (i,j,k)   ,i=ib,ie),j=jb,je),k=kb,ke)
        write(ifoutput)  (((t_vav     (i,j,k)   ,i=ib,ie),j=jb,je),k=kb,ke)
        write(ifoutput)  (((t_sgsav   (i,j,k)   ,i=ib,ie),j=jb,je),k=kb,ke)
        write(ifoutput)  (((p_bav     (i,j,k)   ,i=ib,ie),j=jb,je),k=kb,ke)
        write(ifoutput)  (((p_tav     (i,j,k)   ,i=ib,ie),j=jb,je),k=kb,ke)
        write(ifoutput)  (((-dissresav(i,j,k)   ,i=ib,ie),j=jb,je),k=kb,ke)
        write(ifoutput)  (((d_sgsav   (i,j,k)   ,i=ib,ie),j=jb,je),k=kb,ke)
        close (ifoutput)
      end if


! write average ibm force in x-direction to file
      if (lcanopy == .true.) then
        name2 = 'ibmx_   .'
        name2( 6:8 ) = cmyid
        name2(10:12) = cexpnr
        open  (ifoutput,file=name2,form='unformatted',status='replace')
        write(ifoutput) ((ibmxforcevol(i,k),i=ib-ih,ie+ih),k=kb-kh,ke+kh) 
        close (ifoutput)
      end if


      if (linletgen==1) then
        upupavinl=0.
        call slabsum(upupavinl  ,kb,ke,upupav(ib:ib,jb:je,kb:ke),ib,ib,jb,je,kb,ke,ib,ib,jb,je,kb,ke)
        vpvpavinl=0.
        call slabsum(vpvpavinl  ,kb,ke,vpvpav(ib:ib,jb:je,kb:ke),ib,ib,jb,je,kb,ke,ib,ib,jb,je,kb,ke)
        wpwpavinl=0.
        call slabsum(wpwpavinlh  ,kb,ke+1,wpwpav(ib:ib,jb:je,kb:ke+1),ib,ib,jb,je,kb,ke+1,ib,ib,jb,je,kb,ke+1)
        upwpavinl=0.
        call slabsum(upwpavinlh  ,kb,ke+1,upwpav(ib:ib,jb:je,kb:ke+1),ib,ib,jb,je,kb,ke+1,ib,ib,jb,je,kb,ke+1)
        thlpthlpavinl=0.
        call slabsum(thlpthlpavinl  ,kb,ke,thlpthlpav(ib:ib,jb:je,kb:ke),ib,ib,jb,je,kb,ke,ib,ib,jb,je,kb,ke)
        thlpupavinl=0.
        call slabsum(thlpupavinl  ,kb,ke,thlpupav(ib:ib,jb:je,kb:ke),ib,ib,jb,je,kb,ke,ib,ib,jb,je,kb,ke)
        thlpwpavinl=0.
        call slabsum(thlpwpavinlh ,kb,ke+1,thlpwpav(ib:ib,jb:je,kb:ke+1),ib,ib,jb,je,kb,ke+1,ib,ib,jb,je,kb,ke+1)
        qlpqlpavinl=0.
        call slabsum(qlpqlpavinl  ,kb,ke,qlpqlpav(ib:ib,jb:je,kb:ke),ib,ib,jb,je,kb,ke,ib,ib,jb,je,kb,ke)
        qlpupavinl=0.
        call slabsum(qlpupavinl  ,kb,ke,qlpupav(ib:ib,jb:je,kb:ke),ib,ib,jb,je,kb,ke,ib,ib,jb,je,kb,ke)
        qlpwpavinl=0.
        call slabsum(qlpwpavinlh ,kb,ke+1,qlpwpav(ib:ib,jb:je,kb:ke+1),ib,ib,jb,je,kb,ke+1,ib,ib,jb,je,kb,ke+1)
        thlpthlpavinl=0.
        call slabsum(qtpqtpavinl  ,kb,ke,qtpqtpav(ib:ib,jb:je,kb:ke),ib,ib,jb,je,kb,ke,ib,ib,jb,je,kb,ke)
        thlpupavinl=0.
        call slabsum(qtpupavinl  ,kb,ke,qtpupav(ib:ib,jb:je,kb:ke),ib,ib,jb,je,kb,ke,ib,ib,jb,je,kb,ke)
        thlpwpavinl=0.
        call slabsum(qtpwpavinlh ,kb,ke+1,qtpwpav(ib:ib,jb:je,kb:ke+1),ib,ib,jb,je,kb,ke+1,ib,ib,jb,je,kb,ke+1)
        do k = kb,ke 
          wpwpavinl(k)   = 0.5*(wpwpavinlh(k+1) + wpwpavinlh(k)) ! interpolate to zf
          upwpavinl(k)   = 0.5*(upwpavinlh(k+1) + upwpavinlh(k)) ! interpolate to zf
          thlpwpavinl(k) = 0.5*(thlpwpavinlh(k+1) + thlpwpavinlh(k)) ! interpolate to zf
          qlpwpavinl(k) = 0.5*(qlpwpavinlh(k+1) + qlpwpavinlh(k)) ! interpolate to zf
          qtpwpavinl(k) = 0.5*(qtpwpavinlh(k+1) + qtpwpavinlh(k)) ! interpolate to zf
        end do
        upupavinl = upupavinl /(jge-jgb+1)
        vpvpavinl = vpvpavinl /(jge-jgb+1)
        wpwpavinl = wpwpavinl /(jge-jgb+1)
        upwpavinl = wpwpavinl /(jge-jgb+1)
        thlpthlpavinl = thlpthlpavinl /(jge-jgb+1)
        thlpupavinl = thlpupavinl /(jge-jgb+1)
        thlpwpavinl = thlpwpavinl /(jge-jgb+1)
        qlpqlpavinl = qlpqlpavinl /(jge-jgb+1)
        qlpupavinl = qlpupavinl /(jge-jgb+1)
        qlpwpavinl = qlpwpavinl /(jge-jgb+1)
        qtpqtpavinl = qtpqtpavinl /(jge-jgb+1)
        qtpupavinl = qtpupavinl /(jge-jgb+1)
        qtpwpavinl = qtpwpavinl /(jge-jgb+1)

        if (myid==0) then
          if (linletRA ==.true.) then
            timeinterval = totinletav + timee -btime
          else
            timeinterval = inletav
          end if
          name4 = 'meaninlet        .   '
          write (name4(10:17)  ,'(i8.8)') ntrun
          name4(19:21) = cexpnr
          open (ifoutput,file=name4,form='unformatted',status='replace')
          write(ifoutput) timeinterval      ! length of averaging interval
          write(ifoutput) (Uinl(k),k=kb,ke   )
          write(ifoutput) (Urec(k),k=kb,ke   )
          write(ifoutput) (Wrec(k),k=kb,ke+1   ) 
          write(ifoutput)  ((Utav(i,k),i=ib,ie ),k=kb,ke   )
          close(ifoutput)

          name4 = 'varuinlet        .   '
          write (name4(10:17)  ,'(i8.8)') ntrun
          name4(19:21) = cexpnr
          open (ifoutput,file=name4,form='unformatted',status='replace')
          write(ifoutput) (upupavinl(k),k=kb,ke   )
          write(ifoutput) (vpvpavinl(k),k=kb,ke   )
          write(ifoutput) (wpwpavinl(k),k=kb,ke   )
          write(ifoutput) (upwpavinl(k),k=kb,ke   )
          close(ifoutput)

          if (ltempeq == .true.) then
            name4 = 'tempinlet        .   '
          write (name4(10:17)  ,'(i8.8)') ntrun
            name4(19:21) = cexpnr
            open (ifoutput,file=name4,form='unformatted',status='replace')
            write(ifoutput) timeinterval      ! length of averaging interval
            write(ifoutput) (Tinl(k),k=kb,ke   )
            write(ifoutput) (Trec(k),k=kb,ke   )
            write(ifoutput) ((Ttav(i,k),i=ib,ie ),k=kb,ke   )
            close(ifoutput)

            name4 = 'vartinlet        .   '
          write (name4(10:17)  ,'(i8.8)') ntrun
            name4(19:21) = cexpnr
            open (ifoutput,file=name4,form='unformatted',status='replace')
            write(ifoutput) (thlpthlpavinl(k),k=kb,ke   )
            write(ifoutput) (thlpupavinl(k),k=kb,ke   )
            write(ifoutput) (thlpwpavinl(k),k=kb,ke   )
            close(ifoutput)
  

          end if ! (ltempeq==true)
       
          if (lmoist == .true.) then
            name4 = 'moistinlet        .   '
          write (name4(10:17)  ,'(i8.8)') ntrun
            name4(19:21) = cexpnr
            open (ifoutput,file=name4,form='unformatted',status='replace')
            write(ifoutput) timeinterval      ! length of averaging interval
            write(ifoutput) (QLinl(k),k=kb,ke   )
            write(ifoutput) (QLrec(k),k=kb,ke   )
            write(ifoutput) ((QLtav(i,k),i=ib,ie ),k=kb,ke   )
            write(ifoutput) (QTinl(k),k=kb,ke   )
            write(ifoutput) (QTrec(k),k=kb,ke   )
            write(ifoutput) ((QTtav(i,k),i=ib,ie ),k=kb,ke   )
            close(ifoutput)

            name4 = 'varqinlet        .   '
          write (name4(10:17)  ,'(i8.8)') ntrun
            name4(19:21) = cexpnr
            open (ifoutput,file=name4,form='unformatted',status='replace')
            write(ifoutput) (qlpqlpavinl(k),k=kb,ke   )
            write(ifoutput) (qlpupavinl(k),k=kb,ke   )
            write(ifoutput) (qlpwpavinl(k),k=kb,ke   )
            write(ifoutput) (qtpqtpavinl(k),k=kb,ke   )
            write(ifoutput) (qtpupavinl(k),k=kb,ke   )
            write(ifoutput) (qtpwpavinl(k),k=kb,ke   )
            close(ifoutput)
           end if !lmoist
 
       end if !(myid==0)
      end if !(linletgen==1)

      if (lexitnow) then
        timeleft = 0  !jump out of the time loop
      end if
      if (lexitnow .and. myid == 0 ) then
        open(1, file=trim(name), status='old')
        close(1,status='delete')
        write(*,*) 'Stopped at t=',timee
      end if

      if (myid==0) then
        write(*,'(A,F15.7,A,I4)') 'dump at time = ',timee,' unit = ',ifoutput
      end if
      if (ltec3d == .true.) then
        call tec3d
      end if
      if (ltec2d == .true.) then
        call tec2d
      end if

      if (lcanopy==.true.) then
        call canopy
      end if

    end if

  end subroutine writerestartfiles

! canopy computes the form drag and skin-frictional drag components as well 
! as the moments due to these forces. According to Jackson (1981) the zero-plane 
! displacement d can be comuted by: d = all moments / all forces
! Moments due to: 1) form drag 2) skin friction on cube
! Forces:         1) form drag 2) skin friction on cube 3) skin friction on ground
  subroutine canopy
    use modmpi, only    : myid,comm3d,mpierr,mpi_sum,my_real,nprocs
    use modfields, only : uav,presav,svav,svpwpav,svwav,svpupav,svuav,&
                          svwsgsav,svusgsav
    use modglobal, only : nblocks,zf,dzf,dzfi,dxf,xh,dy,dyi,ie,ib,kb,ke,& 
                          jmax,jgb,jge,jb,je,numol,Uinf,ifoutput,linoutflow,&
                          ntrun, timee,zh,dzhi,dxh,nsv
    use modgenstat, only: lcanopy,istart,ni,timecompl   
    use modibmdata, only: block,nywallsp,nywallsm,ywallsp,ywallsm,&
                          ibmxforcevol
    
    real, allocatable :: f_ground(:)    ! force skin friction on ground
    real, allocatable :: f_p(:)         ! force from pressure difference on cube
    real, allocatable :: f_wallyp(:)    ! force from skin friction on y-plus wall
    real, allocatable :: f_wallym(:)    ! force from skin friction on y-min wall
    real, allocatable :: f_wallzp(:)    ! force from skin friction on z-plus wall
    real, allocatable :: f_ibm(:)       ! force from ibm
    real, allocatable :: f_ground_t(:)    ! sum over all procs
    real, allocatable :: f_p_t(:)         
    real, allocatable :: f_wallyp_t(:)    
    real, allocatable :: f_wallym_t(:)    
    real, allocatable :: f_wallzp_t(:)    
    real, allocatable :: f_ibm_t(:)    

    real, allocatable :: m_p(:)       ! moment from pressure difference on cube
    real, allocatable :: m_wallyp(:)  ! moment from skin friction on y-plus wall
    real, allocatable :: m_wallym(:)  ! moment from skin friction on y-min wall
    real, allocatable :: m_wallzp(:)  ! moment from skin friction on z-plus wall
    real, allocatable :: m_ibm(:)     ! moment from ibm force
    real, allocatable :: m_tot(:)     ! total of moments on this proc
    real, allocatable :: m_tot_t(:)   ! total of moments on all procs

    real, allocatable :: zero_plane(:)  ! moment from skin friction on z-plus wall
    real, allocatable :: c_f(:)         ! skin friction coefficient cf = tau_wall / (0.4 rho *Uinf^2)
    real, allocatable :: p_d(:)         ! form drag coefficient pf = f_p / (0.4 rho *Uinf^2)
    real, allocatable :: u_t(:)         ! friction velocity utau = sqrt((tau_wall + f_p)/rho)
    real, allocatable :: u_2h(:)        ! velocity at 2h on this proc
    real, allocatable :: u_2h_t(:)      ! velocity at 2h
    real, allocatable :: u_1p5h(:)      ! velocity at 1.5h on this proc
    real, allocatable :: u_1p5h_t(:)    ! velocity at 1.5h
    real, allocatable :: utop(:)        ! free stream velocity on this proc
    real, allocatable :: freestream(:)  ! free stream velocity

    real, allocatable :: cav(:,:)       ! volume-averaged concentration in canopy row
    real, allocatable :: adv_wc(:,:)    ! advective flux at top of canopy row
    real, allocatable :: tur_wc(:,:)    ! turbulent flux at top of canopy row
    real, allocatable :: sgs_wc(:,:)    ! sgs flux at top of canopy row
    real, allocatable :: adv_uc(:,:)    ! advective flux at street of canopy row
    real, allocatable :: tur_uc(:,:)    ! turbulent flux at street of canopy row
    real, allocatable :: sgs_uc(:,:)    ! sgs flux at street of canopy row

    ! totals
    real, allocatable :: cav_t(:,:)       ! volume-averaged concentration in canopy row
    real, allocatable :: adv_wc_t(:,:)    ! advective flux at top of canopy row
    real, allocatable :: tur_wc_t(:,:)    ! turbulent flux at top of canopy row
    real, allocatable :: sgs_wc_t(:,:)    ! sgs flux at top of canopy row
    real, allocatable :: adv_uc_t(:,:)    ! advective flux at street of canopy row
    real, allocatable :: tur_uc_t(:,:)    ! turbulent flux at street of canopy row
    real, allocatable :: sgs_uc_t(:,:)    ! sgs flux at street of canopy row

    real dummy,blockvol, freevol
    integer i,j,k,ii, xstat,il,iu,jl,ju,n,i1,i2,i3,i4,ktop,nb
    character(len=25) :: name1

 ! Compute forces and moments per repeating element in the canopy in x-direction
 ! Do loop over repeating elements
 !   - Check if a block wall falls within this repeating element
 !   - Compute force and moment for that wall

    xstat  = floor(real((ie-istart+1)/ni))  ! number of repeating elements in domain (verplaatsen naar modglobal?)  
    if (myid == 0) write(6,*) 'xstat = ', xstat

    allocate(f_ground(1:xstat))    ! force skin friction on ground
    allocate(f_p(1:xstat))         ! force from pressure difference on cube
    allocate(f_wallyp(1:xstat))    ! force from skin friction on y-plus wall
    allocate(f_wallym(1:xstat))    ! force from skin friction on y-min wall
    allocate(f_wallzp(1:xstat))    ! force from skin friction on z-plus wall
    allocate(f_ibm(1:xstat))       ! force from ibm force
    allocate(u_2h(1:xstat))        ! mean velocity at 2h on this proc
    allocate(u_1p5h(1:xstat))      ! mean velocity at 1.5h on this proc



    allocate(f_ground_t(1:xstat))    ! same but for all procs together
    allocate(f_p_t(1:xstat))         
    allocate(f_wallyp_t(1:xstat))    
    allocate(f_wallym_t(1:xstat))    
    allocate(f_wallzp_t(1:xstat))    
    allocate(f_ibm_t(1:xstat))    

    allocate(m_p(1:xstat))       ! moment from pressure difference on cube
    allocate(m_wallyp(1:xstat))  ! moment from skin friction on y-plus wall
    allocate(m_wallym(1:xstat))  ! moment from skin friction on y-min wall
    allocate(m_wallzp(1:xstat))  ! moment from skin friction on z-plus wall
    allocate(m_ibm(1:xstat))     ! moment from ibm
    allocate(m_tot(1:xstat))     ! total of moments on this proc
    allocate(m_tot_t(1:xstat))     ! total of moments on all procs together
    allocate(u_2h_t(1:xstat))    ! mean velocity at 2h
    allocate(u_1p5h_t(1:xstat))  ! mean velocity at 1.5h
    
    if (nsv>0) then
      allocate(cav   (1:xstat,1:nsv))  
      allocate(adv_wc(1:xstat,1:nsv))  
      allocate(tur_wc(1:xstat,1:nsv))  
      allocate(sgs_wc(1:xstat,1:nsv))  
      allocate(adv_uc(1:xstat,1:nsv))  
      allocate(tur_uc(1:xstat,1:nsv))  
      allocate(sgs_uc(1:xstat,1:nsv))  

      allocate(cav_t   (1:xstat,1:nsv))
      allocate(adv_wc_t(1:xstat,1:nsv))
      allocate(tur_wc_t(1:xstat,1:nsv))
      allocate(sgs_wc_t(1:xstat,1:nsv))
      allocate(adv_uc_t(1:xstat,1:nsv))
      allocate(tur_uc_t(1:xstat,1:nsv))
      allocate(sgs_uc_t(1:xstat,1:nsv))
      cav(:,:) =0.
      adv_wc(:,:) =0.
      tur_wc(:,:) =0.
      sgs_wc(:,:) =0.
      adv_uc(:,:) =0.
      tur_uc(:,:) =0.
      sgs_uc(:,:) =0.
   end if 


    allocate(utop(1:xstat))     
    allocate(freestream(1:xstat))     

    f_ground(:) = 0.
    f_p(:)      = 0.
    f_wallzp(:) = 0.
    f_wallyp(:) = 0.
    f_wallym(:) = 0.
    f_ibm(:)    = 0.
    m_p(:)      = 0.
    m_wallzp(:) = 0.
    m_wallyp(:) = 0.
    m_wallym(:) = 0.
    m_ibm(:)    = 0.
    m_tot(:)    = 0.
    utop(:)     = 0.
    freestream(:) = 0.
    do ii=1,xstat   ! loop over repeating elements
      il = istart + (ii-1)*ni
      iu = il + ni-1
!      write(6,*) '(jge-jgb+1)*dy = ', (jge-jgb+1)*dy
!      write(6,*) 'il, iu = ',il,iu
!      write(6,*)'xh(istart+ni)-xh(istart) = ',xh(istart+ni)-xh(istart)

      ! Skin frictional force on ground
      do j=jb,je
      do i=il,iu
!        f_ground(ii) = f_ground(ii) + numol*0.5*(sum(uav(i,jb:je,kb))+sum(uav(i+1,jb:je,kb)))/zf(kb)*dxf(i)*dy 
        f_ground(ii) = f_ground(ii) + numol*0.5*(uav(i,j,kb)+uav(i+1,j,kb))/zf(kb)*dxf(i)*dy   
      end do
      end do
    ! skin frictional force on top walls (using global block definition)
      do n=1,nblocks
!        if (block(n,1)>=il .and. block(n,1) <=iu) then ! check if the wall falls within the i-range (assumg blocks is not longer than iu)
        if (block(n,1)>=il .and. block(n,1) <iu) then ! check if the wall falls within the i-range (assumg blocks is not longer than iu)
          if (block(n,4) < jb+myid*jmax) cycle ! block ends before this proc
          if (block(n,3) > je+myid*jmax) cycle ! block begins after this proc
          jl = block(n,3) -myid*jmax       
          ju = block(n,4) -myid*jmax      
          if (jl<jb) jl=jb  ! keep j within this proc
! BUG!     if (ju>je) jl=je
          if (ju>je) ju=je
          k  = block(n,6)+1
          do i=block(n,1),block(n,2)
            do j=jl,ju
!              dummy           = numol*0.5*(uav(i,j,k)+uav(i+1,j,k))*dzfi(k)*dxf(i)*dy ! force due to skin friction at top of obstacle
              dummy           = numol*(uav(i,j,k)+uav(i+1,j,k))*dzfi(k)*dxf(i)*dy  ! force due to skin friction at top of obstacle 
              f_wallzp(ii) = f_wallzp(ii) + dummy
            end do 
          end do
          m_wallzp(ii) = f_wallzp(ii)*zh(k)          ! moment w.r.t. ground 
        ! pressure force
          i1 = block(n,1)-1
          i2 = block(n,1)
          i3 = block(n,2)
          i4 = block(n,2)+1
          do j=jl,ju
            do k=block(n,5),block(n,6)
!              dummy      = (0.5*(presav(i1,j,k)+presav(i2,j,k)) -&          ! force from pressure difference
!                            0.5*(presav(i3,j,k)+presav(i4,j,k)))*dy*dzf(k)
              dummy      = (presav(i1,j,k)-presav(i4,j,k))*dy*dzf(k)       !force from pressure difference
              f_p(ii) = f_p(ii) + dummy
              m_p(ii) = m_p(ii) + dummy*zf(k)   ! moment w.r.t. ground
            end do
          end do
        end if
      end do

!     ! skin frictional force on top walls (using local block definition)
!      do n=1,nzwallsshear
!        if (zwallsshear(n,2) >= il .and. zwallsshear(n,2) <= iu) then ! check if the wall falls within the i-range
!          k = zwallsshear(n,1)
!          do i=zwallsshear(n,2),zwallsshear(n,3)-1                      ! assume a block is not longer than the repeating element
!            jl = zwallsshear(n,4)
!            ju = zwallsshear(n,5)
!            if (jl < jb) then
!              jl = jb
!            end if
!            if (ju > je) then
!              ju = je
!            end if
!            do j=jl,ju
!              f_wallzp(xstat) = f_wallzp(xstat) + numol*0.5*(uav(i,j,k)+uav(i+1,j,k))*dzfi(k)*dxf(i)*dy 
!            end do
!          end do
!        end if
!      end do


    ! skin frictional force on side wall (yplus) (using local walls!!!)
      do n=1,nywallsp
        if (ywallsp(n,2) >= il .and. ywallsp(n,2) <= iu) then ! check if the wall falls within the i-range
          j = ywallsp(n,1) 
          do i=ywallsp(n,2),ywallsp(n,3)-1 
            do k=ywallsp(n,4),ywallsp(n,5)-1
!              dummy        = numol*0.5*(uav(i,j,k)+uav(i+1,j,k))*dyi*dxf(i)*dzf(k)
              dummy        = numol*(uav(i,j,k)+uav(i+1,j,k))*dyi*dxf(i)*dzf(k)
              f_wallyp(ii) = f_wallyp(ii) + dummy
              m_wallyp(ii) = m_wallyp(ii) + dummy*zf(k)
            end do
          end do
        end if
      end do

    ! skin frictional force on side wall (ymin) (using local walls!!!)
      do n=1,nywallsm
        if (ywallsm(n,2) >= il .and. ywallsm(n,2) <= iu) then ! check if the wall falls within the i-range
          j = ywallsm(n,1) 
          do i=ywallsm(n,2),ywallsm(n,3)-1 
            do k=ywallsm(n,4),ywallsm(n,5)-1
!              dummy        = numol*0.5*(uav(i,j,k)+uav(i+1,j,k))*dyi*dxf(i)*dzf(k)
              dummy        = numol*(uav(i,j,k)+uav(i+1,j,k))*dyi*dxf(i)*dzf(k)
              f_wallym(ii) = f_wallym(ii) + dummy
              m_wallym(ii) = m_wallym(ii) + dummy*zf(k)
            end do
          end do
        end if
      end do

! result in this loop needs to be divided by length and height of repeating
! element    
      do k=kb,ke
        do i=il,iu
          f_ibm(ii) = f_ibm(ii) - ibmxforcevol(i,k)        ! notice the minus sign
!          f_ibm(ii) = f_ibm(ii) - ibmxforce(i,k)*dxh(i)*dzf(k)        ! notice the minus sign  
          m_ibm(ii) = m_ibm(ii) - ibmxforcevol(i,k)*zf(k)  ! to get force ON obstacle 
!          m_ibm(ii) = m_ibm(ii) - ibmxforce(i,k)*dxh(i)*dzf(k) *zf(k)  ! to get force ON obstacle 
        end do
      end do


! assume equidistant mesh in x- and y-directions
      u_1p5h(ii) = 0.
      do k = kb,ke
        if (1.5 < zf(k)) then
          u_1p5h(ii)  =  sum(sum(uav(il:iu,jb:je,k-1),1),1) + ( sum(sum(uav(il:iu,jb:je,k),1),1) -  sum(sum(uav(il:iu,jb:je,k-1),1),1)) * dzhi(k) *&
                                         (1.5 - zf(k-1))  
          exit      ! if the location is found, exit
        end if
      end do

      u_2h(ii) = 0.
      do k = kb,ke
        if (2. < zf(k)) then
          u_2h(ii)  =  sum(sum(uav(il:iu,jb:je,k-1),1),1) + (sum(sum(uav(il:iu,jb:je,k),1),1) -  sum(sum(uav(il:iu,jb:je,k-1),1),1)) *dzhi(k) *&
                                         (2. - zf(k-1))
          exit      ! if the location is found, exit
        end if
      end do



!      if (linoutflow==.true.) then
!        utop(ii) = Uinf
!      else
        utop(ii) = 0.
        do j =jb,je
          do i =il,iu
            utop(ii) = utop(ii) + uav(i,j,ke)*dxf(i)
          end do
        end do
        utop(ii) = utop(ii) / ( (je-jb+1)*(xh(istart+ni)-xh(istart) ) )
!      end if !linoutflow

      if (nsv>0) then
!        cav(:,n) = 0.
        blockvol = 0.
        ktop = block(1,6) + 1
        do nb=1,nblocks
          if (block(nb,1)>=il .and. block(nb,1) <iu) then ! check if the wall falls within the i-range (assuming blocks is not longer than iu)       
            blockvol = blockvol + (xh(block(nb,2)+1) - xh(block(nb,1)))*(zh(block(nb,6)+1) - 0.)*(block(nb,4)+1 - block(nb,3))*dy
          end if          
        end do
        freevol = (jge-jgb+1)*dy *(xh(istart+ni)-xh(istart) )*(zh(block(1,6)+1) - 0.) & ! total volume - blocks volume: assume all blocks have same height
                    - blockvol
!        write(6,*) 'freevol, blockvol = ', freevol, blockvol

        do n=1,nsv
!          cav(ii,n) = sum(sum(sum(svav(il:iu,jb:je,0:block(1,6),n),1),1))/freevol
          do k=kb,block(1,6)
          do j=jb,je
          do i=il,iu
            cav(ii,n) = cav(ii,n) + svav(i,j,k,n)*dxf(i)*dy*dzf(k)/freevol
          end do
          end do
          end do

          do j=jb,je
          do i=il,iu
            tur_wc(ii,n) = tur_wc(ii,n) + svpwpav(i,j,ktop,n)*dxf(i)*dy/ &
                                          ((jge-jgb+1)*dy *(xh(istart+ni)-xh(istart) ))
            sgs_wc(ii,n) = sgs_wc(ii,n) + svwsgsav(i,j,ktop,n)*dxf(i)*dy/ &
                                          ((jge-jgb+1)*dy*(xh(istart+ni)-xh(istart) ))
            adv_wc(ii,n) = adv_wc(ii,n) + (svwav(i,j,ktop,n)-svpwpav(i,j,ktop,n))*dxf(i)*dy/ &
                                          ((jge-jgb+1)*dy*(xh(istart+ni)-xh(istart) ))
          end do
          end do
!          do j=kb,ktop-1
!          do i=il,iu
          do k=kb,ktop-1
          do j=jb,je
            tur_uc(ii,n) = tur_uc(ii,n) + svpupav(iu+1,j,k,n)*dzf(k)*dy/ &
                                          ((jge-jgb+1)*dy*(xh(istart+ni)-xh(istart) ))
            sgs_uc(ii,n) = sgs_uc(ii,n) + svusgsav(iu+1,j,k,n)*dzf(k)*dy/ &
                                          ((jge-jgb+1)*dy*(xh(istart+ni)-xh(istart) ))
            adv_uc(ii,n) = adv_uc(ii,n) + (svuav(iu+1,j,k,n)-svpupav(iu+1,j,k,n))*dzf(k)*dy/ &
                                          ((jge-jgb+1)*dy*(xh(istart+ni)-xh(istart) ))
          end do
          end do
        end do
      end if


    end do ! loop over repeating elementes (ii)

   
   ! take average   
    f_ground(:) = f_ground(:) / &
                  ((jge-jgb+1)*dy *(xh(istart+ni)-xh(istart) ))
    f_p(:)      = f_p(:)      / &
                  ((jge-jgb+1)*dy *(xh(istart+ni)-xh(istart) ))
    f_wallzp(:) = f_wallzp(:) / &
                  ((jge-jgb+1)*dy *(xh(istart+ni)-xh(istart) ))
    f_wallyp(:) = f_wallyp(:) / &
                  ((jge-jgb+1)*dy *(xh(istart+ni)-xh(istart) ))
    f_wallym(:) = f_wallym(:) / &
                  ((jge-jgb+1)*dy *(xh(istart+ni)-xh(istart) ))
    f_ibm(:)    = f_ibm(:)    / &
                  (xh(istart+ni)-xh(istart) )
    u_1p5h(:)   = u_1p5h(:)    / &
                  ((jge-jgb+1) *ni)
    u_2h(:)     = u_2h(:)      / &
                  ((jge-jgb+1) *ni)
    m_p(:)      = m_p(:)      / &
                  ((jge-jgb+1)*dy *(xh(istart+ni)-xh(istart) ))
    m_wallzp(:) = m_wallzp(:) / &
                  ((jge-jgb+1)*dy *(xh(istart+ni)-xh(istart) ))
    m_wallyp(:) = m_wallyp(:) / &
                  ((jge-jgb+1)*dy *(xh(istart+ni)-xh(istart) ))
    m_wallym(:) = m_wallym(:) / &
                  ((jge-jgb+1)*dy *(xh(istart+ni)-xh(istart) ))
    m_ibm(:)    = m_ibm(:)    / &
                  (xh(istart+ni)-xh(istart))


    m_tot(:) = m_p(:) + m_wallzp(:) + m_wallym(:) + m_wallyp(:) + m_ibm(:)

    
 ! sum the results of all processors
    f_ground_t(:) = 0.
    f_p_t(:)      = 0.
    f_wallzp_t(:) = 0.
    f_wallyp_t(:) = 0.
    f_wallym_t(:) = 0.
    f_ibm_t(:)    = 0.
    m_tot_t(:)    = 0.
    u_1p5h_t(:)   = 0.
    u_2h_t(:)     = 0.
    call MPI_ALLREDUCE(f_ground,f_ground_t,xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(f_p,     f_p_t     ,xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(f_wallzp,f_wallzp_t,xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(f_wallym,f_wallym_t,xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(f_wallyp,f_wallyp_t,xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(f_ibm,   f_ibm_t   ,xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(m_tot,   m_tot_t   ,xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(u_1p5h,  u_1p5h_t  ,xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(u_2h,    u_2h_t    ,xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
    call MPI_ALLREDUCE(utop,    freestream,xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
!    call MPI_ALLREDUCE(m_p,     m_p_t     ,xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
!    call MPI_ALLREDUCE(m_wallzp,m_wallzp_t,xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
!    call MPI_ALLREDUCE(m_wallym,m_wallym_t,xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
!    call MPI_ALLREDUCE(m_wallyp,m_wallyp_t,xstat,MY_REAL,MPI_SUM,comm3d,mpierr)

    if (nsv > 0) then
      cav_t(:,:)     = 0.
      adv_wc_t(:,:)  = 0.
      tur_wc_t(:,:)  = 0.
      sgs_wc_t(:,:)  = 0.
      adv_uc_t(:,:)  = 0.
      tur_uc_t(:,:)  = 0.
      sgs_uc_t(:,:)  = 0.
      do n=1,nsv
        call MPI_ALLREDUCE(cav   (:,n),cav_t   (:,n),xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
        call MPI_ALLREDUCE(adv_wc(:,n),adv_wc_t(:,n),xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
        call MPI_ALLREDUCE(tur_wc(:,n),tur_wc_t(:,n),xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
        call MPI_ALLREDUCE(sgs_wc(:,n),sgs_wc_t(:,n),xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
        call MPI_ALLREDUCE(adv_uc(:,n),adv_uc_t(:,n),xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
        call MPI_ALLREDUCE(tur_uc(:,n),tur_uc_t(:,n),xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
        call MPI_ALLREDUCE(sgs_uc(:,n),sgs_uc_t(:,n),xstat,MY_REAL,MPI_SUM,comm3d,mpierr)
      end do
    end if  ! nsv > 0
 
    freestream(:) = freestream(:)/nprocs
 
!    write(6,*) 'proc, sum(sum(ibmxforce)) = ', myid, sum(sum(ibmxforce(ib:ie,kb:ke),1),1) 
!    write(6,*) 'proc, sum(sum(ibmxforcevol)) = ', myid, sum(sum(ibmxforcevol(ib:ie,kb:ke),1),1) 
    if (myid ==0) then
  ! compute displacement thickness d = mom_tot / force_tot
      allocate(zero_plane(1:xstat))     ! zero plance displacement
      allocate(c_f(1:xstat))            ! friction coefficient
      allocate(p_d(1:xstat))            ! form drag coefficient
      allocate(u_t(1:xstat))            ! friction velocity
      zero_plane(:) = (m_tot_t(:) ) / &
                      (f_p_t(:) + f_wallzp_t(:) + f_wallym_t(:) + f_wallyp_t(:) + f_ground_t(:) + f_ibm_t(:)) 

      do i=1,xstat 
      ! skin-frictional drag (C_f = tau_w / (0.5*rho*Uinf^2) )
        c_f(i) = (f_ground_t(i) + f_wallzp_t(i) + f_wallyp_t(i) + f_wallym_t(i)) / & 
                 (0.5*freestream(i)**2)
      ! form drag: P_D = f_p / (0.5*rho*Uinf^2)
        p_d(i) = f_p_t(i) / (0.5*freestream(i)**2)
      ! friction velocity using pressure difference (and not ibm force)
        u_t(i) = sqrt(abs(f_ground_t(i) + f_wallzp_t(i) + f_wallyp_t(i) + & 
                      f_wallym_t(i)+f_p_t(i)))              

!      ! write to screen
!        write(6,*) 'i          = ', i 
!        write(6,*) 'f_ground_t = ', f_ground_t(i) 
!        write(6,*) 'f_p_t      = ', f_p_t(i) 
!        write(6,*) 'f_wallzp_t = ', f_wallzp_t(i) 
!        write(6,*) 'f_wallyp_t = ', f_wallyp_t(i) 
!        write(6,*) 'f_wallym_t = ', f_wallym_t(i) 
!        write(6,*) 'f_ibm_t    = ', f_ibm_t(i) 
!!        write(6,*) 'tau_tot    = ', f_ground_t(i) +  f_p_t(i) + f_wallzp_t(i) + f_wallyp_t(i) + f_wallym_t(i) + f_ibm_t(i) 
!!        write(6,*) 'tau_skin   = ', f_ground_t(i) +  f_wallzp_t(i) + f_wallyp_t(i) + f_wallym_t(i)
!!        write(6,*) 'ibmxforce(:,kb+10) = ', ibmxforce(:,kb+10) 
      end do
   
             

 
      ! write the streamwise data to file
      write(6,*) 'Writing average diagnostics to canopy.txt'
      open(ifoutput,file='canopy.txt',position='append')
      write(ifoutput,'(A)') &
      '#-----------------------------------------------#'
      write(ifoutput,'(A)') &
      '# 1_Elapsed time 2_Time step#'
      write(ifoutput,'(1e16.8,I8)'), timecompl, ntrun
      write(ifoutput,'(A)') &
      '# 1_row 2_<xf>     3_<U_inf>       4_<U_1p5h>      5_<U_2h>      6_F_ground/rho  7_F_top/rho     8_F_sides/rho   9_F_p/rho       10_F_ibm/rho    11_u_tau_p      12_u_tau_ibm    13_d'

      do i=1,xstat
        write(ifoutput,'(I3,12e16.8)') &
        i, &
        xh(istart) + (0.5+(i-1))*(xh(istart+ni)-xh(istart)), &
        freestream(i),&
        u_1p5h_t(i),&
        u_2h_t(i),&
        f_ground_t(i), & 
        f_wallzp_t(i), &
        f_wallyp_t(i) +  f_wallym_t(i), &
        f_p_t(i), &
        f_ibm_t(i) , &
        u_t(i), &
        sqrt(abs(f_ground_t(i) + f_wallzp_t(i) + f_wallyp_t(i) +  f_wallym_t(i) + f_ibm_t(i))),&
        zero_plane(i)
      end do
      close(ifoutput)

      if (nsv>0) then
        do n=1,nsv
          name1(1:13)= 'canopy_s_.txt' 
          write (name1(9:9)  ,'(i1.1)') n
          write(6,*) 'Writing average diagnostics to canopy_sx.txt'
          open(ifoutput,file=name1,position='append')
          write(ifoutput,'(A)') &
          '#-----------------------------------------------#'
          write(ifoutput,'(A)') &
          '# 1_Elapsed time 2_Time step#'
          write(ifoutput,'(1e16.8,I8)'), timecompl, ntrun
          write(ifoutput,'(A)') &
          '# 1_row 2_<xf>      3_<U_inf>       4_<U_1p5h>      5_<U_2h>        6_<c>_volume    7_<w><c>_top    8_<w''c''>_top   9_<w''c''>_sgs_top 10_<u><c>_str  11_<u''c''>_str   12_<u''c''>_sgs_str'

          do i=1,xstat
            write(ifoutput,'(I3,12e16.8)') &
            i, &
            xh(istart) + (0.5+(i-1))*(xh(istart+ni)-xh(istart)), &
            freestream(i),&
            u_1p5h_t(i),&
            u_2h_t(i),&
            cav_t(i,n),&
            adv_wc_t(i,n),&
            tur_wc_t(i,n),&
            sgs_wc_t(i,n),&
            adv_uc_t(i,n),&
            tur_uc_t(i,n),&
            sgs_uc_t(i,n) 
          end do !xstat
          close(ifoutput)
        end do ! nsv
      end if ! nsv > 0
 
    end if ! myid==0

 
  end subroutine canopy

  subroutine writedatafiles
    use modfields, only : u0,v0,w0,sv0,pres0,thl0,ql0,qt0,sv0
    use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh,dsv,trestart2,tnextrestart2,dt_lim,timee,btime,xh,&
                          cexpnr,ntimee,rk3step,ifoutput,nsv,timeleft,dtheta,dqt,dt,cu,cv,ntrun,&
                          totavtime,numol,lstore3d,lstorexz,nkplane,kplane,ltempeq,lmoist,dtmax,dzf,dzhi,&
                          jgb,jtot,jmax
    use modmpi,    only : cmyid,myid,nprocs,mpierr,my_real,comm3d,MPI_OFFSET_KIND,MPI_ORDER_FORTRAN,& 
                          MPI_REAL8,MPI_MODE_CREATE,MPI_MODE_WRONLY,MPI_INFO_NULL,MPI_STATUS_IGNORE,&
                          MPI_COMM_WORLD,MPI_DOUBLE_PRECISION
    use modsubgriddata, only : ekm
!     use mpi

    implicit none

    real, dimension(ib:ie,kb:ke)  :: uinter  ! velocity components interpolated to middle of xz-plane
!    real, dimension(ib:ie,kb:ke)  :: vinter
    real, dimension(ib:ie,kb:ke)  :: winter
    real, dimension(ib:ie,jb:je)  :: uinterxy  ! velocity components interpolated to middle of xz-plane
    real, dimension(ib:ie,jb:je)  :: vinterxy
!    real, dimension(ib:ie,jb:je)  :: winterxy
    real, dimension(ib:ie,jb:je)  :: thlinterxy
    real, dimension(ib:ie,jb:je)  :: qlinterxy
    real, dimension(ib:ie,jb:je)  :: qtinterxy
    real, dimension(ib:ie,jb:je,1:nsv)  :: svinterxy

    integer i,j,k,n,kl
    character(25) name



! MPI_IO
! only native data representation works!
    character(len=*), parameter :: datarep = 'native'
    integer, parameter :: array_ndim = 3
    integer, dimension(1:array_ndim) :: array_of_sizes, array_of_subsizes,array_of_starts
    integer :: subarray_handle1,subarray_handle2,filehandle,etype
    integer(kind=MPI_OFFSET_KIND) :: filesize,displacement
    real, dimension(0:(ie-ib),0:(je-jb),0:(ke-kb)) :: var

    if (timee == 0) return
    if (rk3step /=3) return

    if (mod(ntrun, (ceiling(trestart2/dtmax))) ==0.0) then ! if mod(ntrun/ ntrestart)==0, ntrestart = floor (trestart2/dtmax)
      if (lstore3d == .true.) then


! create subarray
       array_of_sizes(1) = ie-ib+1 +2
       array_of_sizes(2) = jmax + 2
       array_of_sizes(3) = ke-kb+1 +2
       array_of_subsizes(1) = ie-ib+1
       array_of_subsizes(2) = jmax
       array_of_subsizes(3) = ke-kb+1
       array_of_starts(1) = ib             
       array_of_starts(2) = jb              
       array_of_starts(3) = kb              


! first create subarray of local velocity field (without halo's)
       call MPI_TYPE_CREATE_SUBARRAY(array_ndim, array_of_sizes,array_of_subsizes, array_of_starts, &
!                                     MPI_ORDER_FORTRAN, MPI_REAL8, subarray_handle1, mpierr)
                                     MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, subarray_handle1, mpierr)
       call MPI_TYPE_COMMIT(subarray_handle1, mpierr)

! then create subarray for how this data fits in the file to be written (global)
       array_of_sizes(1) = ie-ib+1
       array_of_sizes(2) = jtot
       array_of_sizes(3) = ke-kb+1
       array_of_subsizes(1) = ie-ib+1
       array_of_subsizes(2) = jmax
       array_of_subsizes(3) = ke-kb+1
       array_of_starts(1) = 0             
       array_of_starts(2) = myid*jmax    
       array_of_starts(3) = 0         
       call MPI_TYPE_CREATE_SUBARRAY(array_ndim, array_of_sizes,array_of_subsizes, array_of_starts, &
!                                     MPI_ORDER_FORTRAN, MPI_REAL8, subarray_handle2, mpierr)
                                     MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, subarray_handle2, mpierr)
       call MPI_TYPE_COMMIT(subarray_handle2, mpierr)

!       etype = MY_REAL
!       etype = MPI_REAL8
       etype = MPI_DOUBLE_PRECISION

       name = 'uvel/uvel        .'
       write (name(10:17)  ,'(i8.8)') ntrun
       name(19:21)= cexpnr
      ! open file
       call MPI_FILE_OPEN(MPI_COMM_WORLD, name,MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, filehandle, mpierr)
!     ! guarantee overwriting
!       filesize = 0_MPI_OFFSET_KIND
!       call MPI_FILE_SET_SIZE(filehandle, filesize, mpierr)
!       displacement = 0_MPI_OFFSET_KIND
       displacement = 0
!       CALL MPI_FILE_SET_VIEW( filehandle, displacement, MPI_REAL8, subarray_handle2,  "native", MPI_INFO_NULL, mpierr )
       CALL MPI_FILE_SET_VIEW( filehandle, displacement, etype, subarray_handle2,  "native", MPI_INFO_NULL, mpierr )

!       call MPI_FILE_SET_VIEW(filehandle, displacement, etype, subarray_handle, "native", MPI_INFO_NULL, mpierr)
       write(6,*) 'myid, mpierr = ', myid, mpierr
       write(6,*) 'Before File_Write_all, myid = ', myid
       CALL MPI_FILE_WRITE_ALL( filehandle, u0, 1, subarray_handle1, MPI_STATUS_IGNORE, mpierr ) 

!       call MPI_FILE_WRITE_ALL(filehandle,var(:,:,:), product(array_of_subsizes), etype, MPI_STATUS_IGNORE, mpierr)
      ! close file
       call MPI_FILE_CLOSE(filehandle, mpierr)
       call MPI_TYPE_FREE(subarray_handle1, mpierr)
       call MPI_TYPE_FREE(subarray_handle2, mpierr)
       write(6,*) 'FINISHED MPI_IO, myid = ', myid


!!        open  (ifoutput,file=name,form='binary')
!!        write(ifoutput)  (((u0 (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
!!        open  (unit=15,file=name,form='unformatted')
!        open  (unit=15,file=name,form='unformatted',access='stream')
!        write(15) u0
!        close (ifoutput)
!
!        name(1:9)= 'vvel/vvel' 
!!        open  (ifoutput,file=name,form='binary')
!!        write(ifoutput)  (((v0 (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
!        open  (unit=15,file=name,form='unformatted',access='stream')
!        write(15) v0
!        close (15)
!
!        name(1:9)= 'wvel/wvel'
!!        open  (ifoutput,file=name,form='binary')
!!        write(ifoutput)  (((w0 (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
!        open  (unit=15,file=name,form='unformatted',access='stream')
!        write(15) w0
!        close (15)
! 
!        name(1:9)= 'pres/pres'
!!        open  (ifoutput,file=name,form='binary')
!!        write(ifoutput)  (((pres0 (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
!        open  (unit=15,file=name,form='unformatted',access='stream')
!        write(15) pres0  
!        close (15)
!      
!        name(1:9)= 'visc/visc'
!!        open  (ifoutput,file=name,form='binary')
!!        write(ifoutput)  ((((ekm(i,j,k)-numol)/numol,i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
!        open  (unit=15,file=name,form='unformatted',access='stream')
!        write(15) (ekm-numol)/numol
!        close (15)
!        
!        if (ltempeq == .true.) then
!          name(1:9)= 'thl_/thl_'
!!          open  (ifoutput,file=name,form='binary')
!!          write(ifoutput)  (((thl0(i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb-kh,ke+kh)
!          open  (unit=15,file=name,form='unformatted',access='stream')
!          write(15) thl0
!          close (15)
!        end if ! ltempeq = .true.
!
!        if (nsv>0) then
!         
!          do n=1,nsv 
!            name(1:9)= 'sca1/sca1'
!            write (name(4:4)  ,'(i1.1)') n
!            write (name(9:9)  ,'(i1.1)') n
!            open(unit=15,file=name,form='unformatted',access='stream')
!            write(15) sv0(:,:,:,n)
!            close(15)
!          end do



!          open  (ifoutput,file=name,form='binary')
!          write(ifoutput)  (((sv0(i,j,k,1),i=ib,ie),j=jb,je),k=kb,ke)
!          close (ifoutput)
!
!          name(1:9)= 'sca2/sca2'
!          open  (ifoutput,file=name,form='binary')
!          write(ifoutput)  (((sv0(i,j,k,2),i=ib,ie),j=jb,je),k=kb,ke)
!          close (ifoutput)
!
!          name(1:9)= 'sca3/sca3'
!          open  (ifoutput,file=name,form='binary')
!          write(ifoutput)  (((sv0(i,j,k,3),i=ib,ie),j=jb,je),k=kb,ke)
!          close (ifoutput)
! 
!          name(1:9)= 'sca4/sca4'
!          open  (ifoutput,file=name,form='binary')
!          write(ifoutput)  (((sv0(i,j,k,4),i=ib,ie),j=jb,je),k=kb,ke)
!          close (ifoutput)
!
!          name(1:9)= 'sca5/sca5'
!          open  (ifoutput,file=name,form='binary')
!          write(ifoutput)  (((sv0(i,j,k,5),i=ib,ie),j=jb,je),k=kb,ke)
!          close (ifoutput)


!        end if ! nsv > 0
      end if ! lstore3d == .true.

      if (nkplane> 0) then   ! horizontal planes (if nkplane>0)
      do k=1,nkplane
        do j=jb,je
        do i=ib,ie
          kl = kplane(k)
          uinterxy(i,j) = 0.25*((u0(i+1,j,kl)   + u0(i,j,kl))*dzf(kl-1) + &
                                (u0(i+1,j,kl-1) + u0(i,j,kl-1))*dzf(kl))*dzhi(kl)
          vinterxy(i,j) = 0.25*((v0(i,j+1,kl)   + v0(i,j,kl))*dzf(kl-1) + &
                                (v0(i,j+1,kl-1) + v0(i,j,kl-1))*dzf(k))*dzhi(kl)
!          winterxy(i,j) = 0.5*(w0(i,j,kplane(k)+1) + w0(i,j,kplane(k)))
          thlinterxy(i,j) = 0.5*(thl0(i,j,kl)*dzf(kl-1) + thl0(i,j,kl-1)*dzf(kl))*dzhi(kl)
          qlinterxy(i,j) = 0.5*(ql0(i,j,kl)*dzf(kl-1) + ql0(i,j,kl-1)*dzf(kl))*dzhi(kl)
          qtinterxy(i,j) = 0.5*(qt0(i,j,kl)*dzf(kl-1) + qt0(i,j,kl-1)*dzf(kl))*dzhi(kl)
          do n=1,nsv
            svinterxy(i,j,n) = 0.5*(sv0(i,j,kl,n)*dzf(kl-1) + sv0(i,j,kl-1,n)*dzf(kl))*dzhi(kl)
          end do
        end do
        end do
        
        name = 'xy_0/u_xy'
        write (name(3:4)  ,'(i2.2)') k
        write (name(10:17)  ,'(i8.8)') ntrun
        name(18:18) = '_   .'
        name(19:21) = cmyid
        name(22:22) = '.'
        write (name(23:25)  ,'(i3.3)') kplane(k)
        open  (ifoutput,file=name,form='binary')
        write(ifoutput)  ((uinterxy (i,j),i=ib,ie),j=jb,je)
        close (ifoutput)

        name(1:9) = 'xy_0/v_xy'
        write (name(3:4)  ,'(i2.2)') k
        open  (ifoutput,file=name,form='binary')
        write(ifoutput)  ((vinterxy (i,j),i=ib,ie),j=jb,je)
        close (ifoutput)

        name(1:9) = 'xy_0/w_xy'
        write (name(3:4)  ,'(i2.2)') k
        open  (ifoutput,file=name,form='binary')
!        write(ifoutput)  ((winterxy (i,j),i=ib,ie),j=jb,je)
        write(ifoutput)  ((w0 (i,j,kl),i=ib,ie),j=jb,je)
        close (ifoutput)

!        name(1:9) = 'xy_0/vrxy'
!        write (name(3:4)  ,'(i2.2)') k
!        open  (ifoutput,file=name,form='binary')
!        write(ifoutput)  (((ekm (i,j,kplane(k))-numol)/numol,i=ib,ie),j=jb,je)
!        close (ifoutput)

        if (ltempeq ==.true.) then
          name(1:9) = 'xy_0/T_xy'
          write (name(3:4)  ,'(i2.2)') k
          open  (ifoutput,file=name,form='binary')
          write(ifoutput)  ((thlinterxy (i,j),i=ib,ie),j=jb,je)
          close (ifoutput)
        end if ! ltempeq == .true.
        if (lmoist ==.true.) then
          name(1:9) = 'xy_0/QL_xy'
          write (name(3:4)  ,'(i2.2)') k
          open  (ifoutput,file=name,form='binary')
          write(ifoutput)  ((qlinterxy (i,j),i=ib,ie),j=jb,je)
          close (ifoutput)
          name(1:9) = 'xy_0/QT_xy'
          write (name(3:4)  ,'(i2.2)') k
          open  (ifoutput,file=name,form='binary')
          write(ifoutput)  ((qtinterxy (i,j),i=ib,ie),j=jb,je)
          close (ifoutput)

        end if ! lmoist == .true.

        if (nsv>0) then
          do n=1,nsv
            name(1:9) = 'xy_0/sca1'
            write (name(3:4)  ,'(i2.2)') k
            write (name(9:9)  ,'(i1.1)') n
            open  (ifoutput,file=name,form='binary')
!            write(ifoutput)  ((sv0 (i,j,kplane(k),n) ,i=ib,ie),j=jb,je)
            write(ifoutput)  ((svinterxy (i,j,n) ,i=ib,ie),j=jb,je)
            close (ifoutput)
          end do
        end if
      end do ! loop k=1:nkplane
      end if

!      if (lstorexz ==.true. .and. myid ==0) then
      if (lstorexz ==.true. .and. myid == int(nprocs/2.)) then  ! middle plane
        do k=kb,ke
        do i=ib,ie
          uinter(i,k) = 0.25*(u0(i+1,jb,k) + u0(i,jb,k)+u0(i+1,jb-1,k)+u0(i,jb-1,k));
!          vinter(i,k) = 0.5*(v0(i,jb,k) + v0(i,jb+1,k));
          winter(i,k) = 0.25*(w0(i,jb,k+1) + w0(i,jb,k)+w0(i,jb-1,k+1) + w0(i,jb-1,k));
        end do
        end do
        
        name = 'xz_0/u_xz'
        write (name(10:17)  ,'(i8.8)') ntrun
        open  (ifoutput,file=name,form='binary')
        write(ifoutput)  ((uinter (i,k),i=ib,ie),k=kb,ke)
        close (ifoutput)

        name = 'xz_0/v_xz'
        write (name(10:17)  ,'(i8.8)') ntrun
        open  (ifoutput,file=name,form='binary')
        write(ifoutput)  ((v0 (i,jb,k),i=ib,ie),k=kb,ke)
        close (ifoutput)

        name = 'xz_0/w_xz'
        write (name(10:17)  ,'(i8.8)') ntrun
        open  (ifoutput,file=name,form='binary')
        write(ifoutput)  ((winter (i,k),i=ib,ie),k=kb,ke)
        close (ifoutput)

        name = 'xz_0/p_xz'
        write (name(10:17)  ,'(i8.8)') ntrun
        open  (ifoutput,file=name,form='binary')
        write(ifoutput)  (( 0.5*(pres0 (i,jb,k)+pres0(i,jb-1,k)),i=ib,ie),k=kb,ke)
        close (ifoutput)

        name = 'xz_0/vrxz'
        write (name(10:17)  ,'(i8.8)') ntrun
        open  (ifoutput,file=name,form='binary')
        write(ifoutput)  (( (0.5*(ekm (i,jb,k)-numol)+(ekm(i,jb-1,k)-numol))  /numol,i=ib,ie),k=kb,ke)
        close (ifoutput)

        if (ltempeq ==.true.) then     
          name = 'xz_0/T_xz'
          write (name(10:17)  ,'(i8.8)') ntrun
          open  (ifoutput,file=name,form='binary')
          write(ifoutput)  (( 0.5*(thl0 (i,jb,k) + thl0 (i,jb-1,k)),i=ib,ie),k=kb,ke)
          close (ifoutput)
        end if ! ltempeq == .true.
        if (lmoist ==.true.) then     
          name = 'xz_0/QL_xz'
          write (name(10:17)  ,'(i8.8)') ntrun
          open  (ifoutput,file=name,form='binary')
          write(ifoutput)  (( 0.5*(ql0 (i,jb,k) + ql0 (i,jb-1,k)),i=ib,ie),k=kb,ke)
          close (ifoutput)
          name = 'xz_0/QT_xz'
          write (name(10:17)  ,'(i8.8)') ntrun
          open  (ifoutput,file=name,form='binary')
          write(ifoutput)  (( 0.5*(qt0 (i,jb,k) + qt0 (i,jb-1,k)),i=ib,ie),k=kb,ke)
          close (ifoutput)
        end if ! lmoist == .true.

        if (nsv>0) then
          do n=1,nsv
            name(1:9)= 'xz_0/sca1'
            write (name(9:9)  ,'(i1.1)') n
            write (name(10:17)  ,'(i8.8)') ntrun
            open  (ifoutput,file=name,form='binary')
            write(ifoutput)  (( 0.5*(sv0 (i,jb,k,n) + sv0 (i,jb-1,k,n)),i=ib,ie),k=kb,ke)
            close (ifoutput)
          end do
        end if

      end if ! lstorexz == .true.
 
    end if ! timee >= tnexrestart2


  end subroutine writedatafiles

end module modsave
