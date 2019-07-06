!> \file modetrees.f90
!!tg3315, ns4513, 20 Mar 2017 

!> Input trees into DALES model.

module modtrees
implicit none
save

contains
    subroutine createtrees
    use modglobal,  only : ltrees,ntrees,tree,cexpnr,ifinput,zh,zf,dzh,dzfi,dzhi,dzf,Qstar,&
                           dec,lad,kb,ke,cp,rhoa,ntree_max,dQdt,imax,jtot,tr_A,dy,xh
    use modfields,  only : um,vm,wm,thlm,qt0,svp,up,vp,wp,thlp,qtp,Rn,clai,qc,qa,ladzh,ladzf
    use modmpi,     only : myid,MPI_INTEGER,comm3d,mpierr,MY_REAL
    use modsurfdata,only : wtsurf
    use modibmdata ,only : bctfz
    implicit none
    integer :: n,k
    real :: Rq
    character(80) chmess

    if ((ltrees .eqv. .false.) .or. (ntrees==0)) return

      allocate(tree(ntrees,6))

      ! read global trees
      if(myid==0) then
        write(*,*), '1, myid, ntrees, ltrees, cexpnr', myid, ntrees, ltrees, cexpnr
        if (ntrees>0) then
          open (ifinput,file='trees.inp.'//cexpnr)
          read (ifinput,'(a80)') chmess
          read (ifinput,'(a80)') chmess
          do n=1,ntrees
            read (ifinput,*) &
                  tree(n,1), &
                  tree(n,2), &
                  tree(n,3), &
                  tree(n,4), &
                  tree(n,5), &
                  tree(n,6) 
          end do

          write (6,*) 'Tree number,   il, iu, jl, ju, kl, ku '
          do n=1,ntrees
            write (1,*) &
                  n , &
                  tree(n,1), &
                  tree(n,2), &
                  tree(n,3), &
                  tree(n,4), &
                  tree(n,5), &
                  tree(n,6)                 
          end do
        end if

      end if ! end if myid==0

      call MPI_BCAST(ntrees ,1,MPI_INTEGER ,0,comm3d,mpierr)
      call MPI_BCAST(tree ,6*ntrees,MPI_INTEGER ,0,comm3d,mpierr)

!      nlad = tree(1,6)+2-tree(1,5) !this could vary by tree...

      ! read leaf area density of trees ! tg3315 commented - for now assume constant lad
      !if (myid==0) then

      !  allocate(ladt(nlad)) 

      !  if (nlad>0) then
      !    open (ifinput,file='lad.inp.'//cexpnr)
      !    read (ifinput,'(a80)') chmess
      !    do n=1,nlad
      !      read (ifinput,*) &
      !            ladt(n)
      !      end do

      !    write (*,*) 'lad'
      !    do n=1,nlad
      !      write (1,*) &
      !            n , &
      !            ladt(n)
      !    end do
      !  end if

      ! input into leaf area density array (cell faces)
      !lad(tree(1,5):tree(1,6)+1) = ladt !defined at cell faces

      ! calculate cumulative leaf area index (cell faces)
      !do k = tree(1,5),tree(1,6)
      !  clai(k) = sum(0.5*(lad(k:tree(1,6))+lad(k+1:tree(1,6)+1))*(dzh(k:tree(1,6))))
      !end do

    ! calculate clai for tallest possible tree
    ! assume no grid stretching in lowest layer with trees
    if (myid==0) then

      ntree_max = maxval(tree(:,6))-minval(tree(:,5))+1

      ! calculate cumulative leaf area index (cell faces)
      !do k = 1,ntree_max+1
      !  clai(k) = (ntree_max-k+1)*lad*dzh(kb)
      !end do

      if (.false.) then
        !ladz(1:ntree_max) = (/  0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0, 1.3, 1.2, 0.7 /)
!        ladz(1:ntree_max) = (/  0.85/20, 0.885/20, 0.92/20, 1.02/20, 1.12/20,  1.31/20, 1.5/20, 1.85/20, 2.2/20, 2.45/20, 2.7/20, 2.85/20, 3./20, 2.975/20, 2.95/20, 2.85/20, 2.75/20, 2.325/20, 2./20, 1./20 /) !LAI =2 from Shaw, 1992
        ladzh(1:ntree_max) = (/  0.1, 0.11, 0.12, 0.14, 0.16, 0.19, 0.22, 0.265, 0.31, 0.335, 0.36, 0.37, 0.38, 0.3725, 0.365, 0.3375, 0.31, 0.2375, 0.165, 0.0875 /) !LAI =2 from Shaw, 1992
! assumes lad at +1 is 0.!
      else
        ladzh(1:ntree_max+1) = lad !+1 cos how clai is calculated using lad at cell faces so we need lad at top face
      end if

      do k=1,ntree_max
        ladzf(k) = 0.5*(ladzh(k)+ladzh(k+1))
      end do

      ! clai at cell faces using lad at cell centre
      do k = 1,ntree_max
        clai(k) = sum(0.5*(ladzh(k:ntree_max)+ladzh(k+1:ntree_max+1))*(dzf(tree(1,6)-ntree_max+k:tree(1,6))))
      end do

      ! Net radiation at cell faces! !W/m^2
      do k = 1,ntree_max+1
        Rn(k) = Qstar*exp(-dec*clai(k))
      end do

      ! Change in radiation over each layer (at cell centres) W/m^2
      do k = 1,ntree_max
        qc(k) = Rn(k+1) - Rn(k)
      end do

      ! incorporate the distributed storage term
      do k = 1,ntree_max
        Rq = qc(k)/Qstar
        qa(k) = (1 - Rq*0.11) * qc(k) - 0.11*Rq*dQdt + Rq*12.3
      end do

      write(*,*) 'ntree_max', ntree_max
      write(*,*) 'qa', qa(1:ntree_max)
      write(*,*) 'qc', qc(1:ntree_max)
      write(*,*) 'clai', clai(1:ntree_max+1)
      write(*,*) 'Rn', Rn(1:ntree_max+1)
      write(*,*) 'ladzf', ladzf
      write(*,*) 'ladzh', ladzh

    end if

    !apply storage to reset wtsurf and bctfz
    wtsurf = -((1-0.7)*Qstar-0.33*dQdt+38)/(rhoa*cp)
    bctfz  = -((1-0.7)*Qstar-0.33*dQdt+38)/(rhoa*cp)

    ! calc tree cover area
    if (myid==0) then
    do n = 1,ntrees
      tr_A = tr_A + (xh(tree(n,2)+1) - xh(tree(n,1)))*dy*(tree(n,4) - tree(n,3) + 1)
    end do
    end if

    write(*,*) 'wtsurf', wtsurf
    write(*,*) 'bctfz', bctfz

    ! broadcast variables
    call MPI_BCAST(ntree_max , 1,MY_REAL ,0,comm3d,mpierr)
    call MPI_BCAST(tr_A , 1,MY_REAL ,0,comm3d,mpierr)
    call MPI_BCAST(Rn , ke+1,MY_REAL ,0,comm3d,mpierr)
    call MPI_BCAST(qc , ke+1,MY_REAL ,0,comm3d,mpierr)
    call MPI_BCAST(clai  , ke+1,MY_REAL ,0,comm3d,mpierr)  
    call MPI_BCAST(qa  , ke+1,MY_REAL ,0,comm3d,mpierr)
    call MPI_BCAST(ladzf  , ke+1,MY_REAL ,0,comm3d,mpierr)
    call MPI_BCAST(ladzh , ke+1,MY_REAL ,0,comm3d,mpierr)

    write(*,*) 'tr_A', tr_A

    end subroutine createtrees

    subroutine trees
    use modglobal,  only : ib,ie,jb,je,kb,ke,dzf,xh,dxf,numol,prandtlmol,rlv,cp,jtot,tree,&
                           ntrees,ltrees,jtot,cd,ud,lmoist,nsv,dxf,dy,dzf,dzfi,zf,dy,zh,&
                           ltempeq,pref0,r_s,lad,rhoa,ntree_max,lsize,Qstar,dQdt,BCxs,tr_A,&
                           imax,jtot
    use modfields,  only : um,vm,wm,thlm,qtm,svp,up,vp,wp,thlp,qtp,svm,Rn,qc,qa,thlm,clai,ladzf,ladzh,&
                           tr_u,tr_v,tr_w,tr_qt,tr_qtR,tr_qtA,tr_thl,tr_sv,thlpcar
    use modmpi,     only : myid,nprocs,mpi_sum,mpierr,comm3d,mpierr,my_real,nprocs
    use modsurfdata,only : wtsurf,wttop,wqtop
    use modibmdata, only : bctfz
    use modsubgriddata, only : ekh
    implicit none
    integer :: i,j,k,n,m,il,iu,jl,ju,kl,ku
    real :: e_sat,e_vap,qh,qe,shade,r_a,s,D,numoli,gam,Rq,qhmin,qhmax
    real :: Vq_dum,Vq_dum2,Vq_dum3,VT_dum,VT_dum2,VT_dum3,wTbot_dum,wTbot_dum2,wTbot,wTbot_dum3
!wttop_dum,wttop_dum2,wqtop_dum,wqtop_dum2,wttop_dum3,wqtop_dum3


    numoli = 1/numol
   
    qhmin=0.;qhmax=0.

    ! Tree drag
    if (ltrees .eqv. .false.) return

    VT_dum2 = 0.
    Vq_dum2 = 0.
    wTbot_dum2 = 0.
    gam = 72.1886

    !write(*,*), 'myid, ntrees', myid, ntrees

    do n = 1,ntrees
     
    !  write(*,*), 'myid2, tree(n,1), tree(n,2), tree(n,3), tree(n,4), tree(n,5), tree(n,6)', myid, tree(n,1), tree(n,2), tree(n,3), tree(n,4), tree(n,5), tree(n,6)
 
      ! w drag
      il = tree(n,1)
      iu = tree(n,2)
      kl = tree(n,5)
      ku = tree(n,6) + 1
      jl = tree(n,3)-myid*jtot/nprocs
      ju = tree(n,4)-myid*jtot/nprocs
      if (ju < jb .or. jl > je) then
        cycle
      else
        if (ju > je) ju=je
        if (jl < jb) jl=jb  

        do k = kl,ku 
          do j = jl,ju
            do i = il,iu
              ! could do (il:iu,jl:ju,kl:ku)...?
              tr_w(i,j,k) = - cd * ladzh(ntree_max-(ku-k)+1) * wm(i,j,k) * &
                         sqrt( wm(i,j,k)**2 &
                        + (0.25*(um(i,j,k) + um(i+1,j,k) + um(i,j,k-1) + um(i+1,j,k-1)))**2 &
                        + (0.25*(vm(i,j,k) + vm(i,j+1,k) + vm(i,j,k-1) + vm(i,j+1,k-1)))**2 )
              wp(i,j,k) = wp(i,j,k) + tr_w(i,j,k) !cd * ladz(ntree_max-(ku-k)) * wm(i,j,k) * &
                        ! sqrt( wm(i,j,k)**2 &
                        !+ (0.25*(um(i,j,k) + um(i+1,j,k) + um(i,j,k-1) + um(i+1,j,k-1)))**2 &
                        !+ (0.25*(vm(i,j,k) + vm(i,j+1,k) + vm(i,j,k-1) + vm(i,j+1,k-1)))**2 )
            end do
          end do
        end do
      end if

      ! v drag
      il = tree(n,1)
      iu = tree(n,2)
      kl = tree(n,5)
      ku = tree(n,6)
      jl = tree(n,3) - myid*jtot/nprocs
      ju = tree(n,4) + 1 - myid*jtot/nprocs
      if (ju < jb .or. jl > je) then
        cycle
      else
        if (ju > je) ju=je
        if (jl < jb) jl=jb 

        do k = kl,ku
          do j = jl,ju
            do i = il,iu
              tr_v(i,j,k) = - cd * ladzf(ntree_max-(ku-k)) * vm(i,j,k) * &
                           sqrt( vm(i,j,k)**2 &
                          + (0.25*(um(i,j,k) + um(i+1,j,k) + um(i,j-1,k) + um(i+1,j-1,k)))**2 &
                          + (0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i,j-1,k) + wm(i,j-1,k+1)))**2 )
              vp(i,j,k) = vp(i,j,k) + tr_v(i,j,k) !- cd * ladz(ntree_max-(ku-k)) * vm(i,j,k) * &
                          ! sqrt( vm(i,j,k)**2 &
                          !+ (0.25*(um(i,j,k) + um(i+1,j,k) + um(i,j-1,k) + um(i+1,j-1,k)))**2 &
                          !+ (0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i,j-1,k) + wm(i,j-1,k+1)))**2 )
            end do
          end do
        end do
      end if

      ! u drag
      il = tree(n,1)
      iu = tree(n,2) + 1
      kl = tree(n,5)
      ku = tree(n,6)
      jl = tree(n,3) - myid*jtot/nprocs
      ju = tree(n,4) - myid*jtot/nprocs
      if (ju < jb .or. jl > je) then
        cycle
      else
        if (ju > je) ju=je
        if (jl < jb) jl=jb 
 
        do k = kl,ku
          do j = jl,ju
            do i = il,iu
              tr_u(i,j,k) = - cd * ladzf(ntree_max-(ku-k)) * um(i,j,k) * &
                           sqrt( um(i,j,k)**2 &
                          + (0.25*(vm(i,j,k) + vm(i,j+1,k) + vm(i-1,j,k) + vm(i-1,j+1,k)))**2 &
                          + (0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i-1,j,k) + wm(i-1,j,k+1)))**2 )
              up(i,j,k) = up(i,j,k) + tr_u(i,j,k) ! - cd * ladz(ntree_max-(ku-k)) * um(i,j,k) * &
                          ! sqrt( um(i,j,k)**2 &
                          !+ (0.25*(vm(i,j,k) + vm(i,j+1,k) + vm(i-1,j,k) + vm(i-1,j+1,k)))**2 &
                          !+ (0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i-1,j,k) + wm(i-1,j,k+1)))**2 )
            end do
          end do
        end do
      end if

      ! SCALARS

      ! Moisture and temp
      if (lmoist .and. ltempeq) then

        il = tree(n,1)
        iu = tree(n,2)
        kl = tree(n,5)
        ku = tree(n,6)
        jl = tree(n,3) - myid*jtot/nprocs
        ju = tree(n,4) - myid*jtot/nprocs
        if (ju < jb .or. jl > je) then
          cycle
        else
          if (ju > je) ju=je
          if (jl < jb) jl=jb 
       
          do k = kl,ku
            do j = jl,ju
              do i = il,iu

                 ! BOWEN RATIO
                 !qtp(i,j,k) = qtp(i,j,k) + qc(k) / ( rlv * 1.225 * (1. + Bowen) )
                 !thlp(i,j,k) = thlp(i,j,k) + qc(k) / ( (1. + 1. / Bowen) * 1.225 * cp )
                 
                 !esat_air = 100*6.112*exp((17.67*(thlm(i,j,k))-273.15)/(thlm(i,j,k)-29.65)) !Clausias-Clapeyron
                 !wsat = 0.622*psat_air/pref0 ! w_sat = 0.622*p_sat/p0
                 !D = (1-(qtm(i,j,k)/wsat_air))*psat_air ! (1-RH)*e_s
                 
                 e_sat = 100*6.112*exp((17.67*(thlm(i,j,k)-273.15))/(thlm(i,j,k)-29.65)) !Clausias-Clapeyron

                 e_vap = (qtm(i,j,k) * pref0) / (0.378 * qtm(i,j,k) + 0.622)
                 ! qtcheck.m and https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity

                 D = max(e_sat - e_vap,0.)

                 if (e_sat<e_vap) then
                   write(*,*) 'D, e_sat, e_vap, thlm, qtm', D, e_sat, e_vap, thlm(i,j,k), qtm(i,j,k)
                 end if

                 s = 1000*(4098*(0.6108*exp((17.27*(thlm(i,j,k)-273.15))/((thlm(i,j,k)-273.15)+273.3))))/(((thlm(i,j,k)-273.15)+273.3)**2)

                 r_a = 130*sqrt((lsize/(sqrt((0.5*(um(i,j,k)+um(i+1,j,k)))**2+(0.5*(vm(i,j,k)+vm(i,j+1,k)))**2+(0.5*(wm(i,j,k)+wm(i,j,k+1)))**2))))

!                 write(*,*) 'thlm, s, D, e_vap, e_sat, lad, ntree_max, r_a, r_s, gam'
!                 write(*,*) thlm(i,j,k), s, D, e_vap, e_sat, lad, ntree_max, r_a, r_s, gam

                 qe = (s*(qa(ntree_max-(ku-k))/(dzf(k)*ladzf(ntree_max-(ku-k)))) + 2*rhoa*cp*(1/r_a)*D) / (s + (gam*2*(r_a+r_s)/r_a))

                 !write(*,*) 'D', D
                 !write(*,*) 'qe', qe
                 !write(*,*) 'QA/ADV', s*(qa(ntree_max-(ku-k))/(dzf(k)*ladz(ntree_max-(ku-k))))/ (2*rhoa*cp*(1/r_a)*D)                

                 qh = qa(ntree_max-(ku-k))/(dzf(k)*ladzf(ntree_max-(ku-k))) - qe

                 tr_qt(i,j,k) = ladzf(ntree_max-(ku-k))*qe/(rhoa*rlv)
                 tr_qtR(i,j,k) = ladzf(ntree_max-(ku-k))*( (s*(qa(ntree_max-(ku-k))/(dzf(k)*ladzf(ntree_max-(ku-k))))) / (s + (gam*2*(r_a+r_s)/r_a)) ) /(rhoa*rlv)
                 tr_qtA(i,j,k) = ladzf(ntree_max-(ku-k))*( (s*(2*rhoa*cp*(1/r_a)*D)) / (s + (gam*2*(r_a+r_s)/r_a)) ) / (rhoa*rlv)
                 tr_thl(i,j,k) = ladzf(ntree_max-(ku-k))*qh/(rhoa*cp)

                 qtp(i,j,k) = qtp(i,j,k) + tr_qt(i,j,k) !ladz(ntree_max-(ku-k))*qe/(rhoa*rlv)
                 thlp(i,j,k) = thlp(i,j,k) + tr_thl(i,j,k) !ladz(ntree_max-(ku-k))*qh/(rhoa*cp)

                 VT_dum2 = VT_dum2 - tr_thl(i,j,k)*dzf(k)*dxf(i)*dy
                 Vq_dum2 = Vq_dum2 - tr_qt(i,j,k)*dzf(k)*dxf(i)*dy

              end do
            end do
          end do
        !end if
      
      Rq = Rn(ntree_max+1-(ku+1-kl)) / Qstar

      shade = ( (1 - Rq*0.11 ) * Rn(ntree_max+1-(ku+1-kl) ) -0.11*Rq*dQdt + Rq*12.3 )/ (rhoa*cp)

      wTbot_dum2 = wTbot_dum2 - shade*(xh(iu+1)-xh(il))*dy*(ju-jl+1)

      !write(*,*) 'Rn(ntree_max+1-(ku-kl) )', Rn(ntree_max+1-(ku+1-kl) )
      !write(*,*) 'shade', shade

      thlp(il:iu,jl:ju,kb+1) = thlp(il:iu,jl:ju,kb+1) + (bctfz + shade) * dzfi(kb+1)
 
      end if


    !   write(*,*) 'ntree_max+1-(ku-kl), Rn(ntree_max+1-(ku-kl))', ntree_max-(ku-kl), Rn(ntree_max-(ku-kl))
 
    !  write(*,*) 'Rn(1:ntree_max+1)', Rn(1:ntree_max+1)
    !  write(*,*) 'clai(1:ntree_max+1)', clai(1:ntree_max+1)
    !  write(*,*) 'qc(1:ntree_max+1)', qc(1:ntree_max+1)
      !write(*,*) 'bctfz', bctfz
      !write(*,*) 'shade', shade

      end if !lmoist

      ! Temp only
      !if (ltempeq .and. (.not. lmoist)) then
    
      !  il = tree(n,1)
      !  iu = tree(n,2)
      !  kl = tree(n,5)
      !  ku = tree(n,6)
      !  jl = tree(n,3) - myid*jtot/nprocs
      !  ju = tree(n,4) - myid*jtot/nprocs
      !  treec = (kl+ku) / 2
      !  if (ju < jb .or. jl > je) then
      !    cycle
      !  else
      !    if (ju > je) ju=je
      !    if (jl < jb) jl=jb 
      ! 
      !    do k = kl,ku
      !      do j = jl,ju
      !        do i = il,iu
      !            thlp(i,j,k) = thlp(i,j,k) + qc(k) / ( 1.225 * cp )
      !        end do
      !      end do
      !    end do
      !  end if
      !end if

      !if (ltempeq) then
      !  thlp(il:iu,jl:ju,kb) = thlp(il:iu,jl:ju,kb) + (wtsurf + shade)*dzfi(kb)
      !end if

      ! scalar deposition
      if (nsv>0) then

      if (BCxs==2) then ! no deposition effect from first or last tree to avoid interaction with BCs
        if (n==1) cycle
        if (n==ntrees) cycle
      end if
      
        il = tree(n,1)
        iu = tree(n,2)
        kl = tree(n,5)
        ku = tree(n,6)
        jl = tree(n,3) - myid*jtot/nprocs
        ju = tree(n,4) - myid*jtot/nprocs
        if (ju < jb .or. jl > je) then
          cycle
        else
          if (ju > je) ju=je
          if (jl < jb) jl=jb 

          do m = 1,nsv
            do k = kl,ku
              do j = jl,ju
                do i = il,iu
                  tr_sv(i,j,k,m) = - svm(i,j,k,m) * ladzf(ntree_max-(ku-k)) * ud
                  svp(i,j,k,m) = svp(i,j,k,m) + tr_sv(i,j,k,m)
                  !svp(i,j,k,m) = svp(i,j,k,m) - svm(i,j,k,m) * ladz(ntree_max-(ku-k)) * ud
                end do
              end do
            end do
          end do
        end if
      end if !nsv

    end do ! ntrees

    if (lmoist .and. ltempeq) then

      wTbot_dum3 = 0.
      wTbot_dum = 0.
      wTbot = 0.
      VT_dum3 = 0.
      Vq_dum3 = 0.
      VT_dum = 0.
      Vq_dum = 0.

      call MPI_ALLREDUCE(wTbot_dum2,   wTbot_dum3,1,MY_REAL,MPI_SUM,comm3d,mpierr)    
      call MPI_ALLREDUCE(VT_dum2, VT_dum3  ,1,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(Vq_dum2, Vq_dum3  ,1,MY_REAL,MPI_SUM,comm3d,mpierr)

      wTbot_dum = wTbot_dum3
      Vq_dum = Vq_dum3
      VT_dum = VT_dum3

      !write(*,*) 'wTbot_dum', wTbot_dum
      !write(*,*) 'Vq_dum', Vq_dum
      !write(*,*) 'VT_dum', VT_dum

      wTbot = wTbot_dum + ((xh(ie+1)-xh(ib))*jtot*dy - tr_A)*bctfz

      !write(*,*) 'wTbot', wTbot

      wttop = 0.5*(wTbot + VT_dum)/((xh(ie+1)-xh(ib))*jtot*dy)

      !write(*,*) 'wttop', wttop

!      wttop = (wttop_dum3 + ((xh(ie+1)-xh(ib))*jtot*dy - tr_A)*bctfz) / (2*(xh(ie+1)-xh(ib))*jtot*dy)

      thlpcar = wttop/(zh(ke+1)-zh(kb+1))
!      wqtop = (wqtop_dum3 ) / ((xh(ie+1)-xh(ib))*jtot*dy)    !2*
      wqtop = (Vq_dum ) / ((xh(ie+1)-xh(ib))*jtot*dy)

      !write(*,*) 'wqtop', wqtop

    end if

    end subroutine trees

  end module modtrees
