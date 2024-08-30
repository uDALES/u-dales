!> \file modetrees.f90
!!tg3315, ns4513, 20 Mar 2017 

!> Input trees into DALES model.

module modtrees
use mpi
implicit none
save

contains
    subroutine createtrees
    use modglobal,  only : ltrees,ntrees,tree,cexpnr,ifinput,zh,zf,dzh,dzfi,dzhi,dzf,Qstar,&
                           dec,lad,kb,ke,cp,rhoa,ntree_max,dQdt,tr_A,dy,xh
    use modfields,  only : um,vm,wm,thlm,qt0,svp,up,vp,wp,thlp,qtp,Rn,clai,qc,qa,ladzh,ladzf
    use modmpi,     only : myid,comm3d,mpierr,MY_REAL
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
        ! write(*,*) '1, myid, ntrees, ltrees, cexpnr', myid, ntrees, ltrees, cexpnr
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

          ! write (6,*) 'Tree number,   il, iu, jl, ju, kl, ku '
          ! do n=1,ntrees
          !   write (1,*) &
          !         n , &
          !         tree(n,1), &
          !         tree(n,2), &
          !         tree(n,3), &
          !         tree(n,4), &
          !         tree(n,5), &
          !         tree(n,6)                 
          ! end do
          close (ifinput)
        end if

      end if ! end if myid==0

      ! call MPI_BCAST(ntrees ,1,MPI_INTEGER ,0,comm3d,mpierr)
      call MPI_BCAST(tree ,6*ntrees,MPI_INTEGER ,0,comm3d,mpierr)

      !! capability to read lad of trees from a lad.inp.xxx file
      !! tg3315 commented - for now assume constant lad or hardcode as seen below
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

      ! initialise tree variables that are constant in time
      ! calculate clai for tallest possible tree
      ! assume no grid stretching in lowest layer with trees
      if (myid==0) then

        ntree_max = maxval(tree(:,6))-minval(tree(:,5))+1

        ! hard code non-uniform lad profiles - temporary
        if (.false.) then
          ladzh(1:ntree_max) = (/  0.1, 0.11, 0.12, 0.14, 0.16, 0.19, 0.22, 0.265, 0.31, 0.335, 0.36, 0.37, 0.38, 0.3725, 0.365, 0.3375, 0.31, 0.2375, 0.165, 0.0875 /) !LAI =2 from Shaw, 1992
          ! assumes lad at ntree_max+1 is 0.!
        else
          ladzh(1:ntree_max+1) = lad !up to ntree_max+1 cos ladzh is at cell faces
        end if

        ! interpolate to find lad at cell centres
        do k=1,ntree_max
          ladzf(k) = 0.5*(ladzh(k)+ladzh(k+1))
        end do

        ! clai at cell faces using lad at cell centre
        do k = 1,ntree_max
          clai(k) = sum(0.5*(ladzh(k:ntree_max)+ladzh(k+1:ntree_max+1))*(dzf(maxval(tree(:,6))-ntree_max+k:maxval(tree(:,6)))))
        end do

        ! Net radiation at cell faces! !W/m^2
        do k = 1,ntree_max+1
          Rn(k) = Qstar*exp(-dec*clai(k))
        end do

        ! Change in radiation over each layer (at cell centres) W/m^2
        do k = 1,ntree_max
          qc(k) = Rn(k+1) - Rn(k)
        end do

        ! incorporate the distributed storage term - updated for zero heat storage in tree
        do k = 1,ntree_max
          qa(k) = qc(k) !(1 - Rq*0.11) * qc(k) - 0.11*Rq*dQdt + Rq*12.3
        end do

        ! write relevant fields
        ! write(*,*) 'ntree_max', ntree_max
        ! write(*,*) 'qa', qa(1:ntree_max)
        ! write(*,*) 'qc', qc(1:ntree_max)
        ! write(*,*) 'clai', clai(1:ntree_max+1)
        ! write(*,*) 'Rn', Rn(1:ntree_max+1)
        ! write(*,*) 'ladzf', ladzf
        ! write(*,*) 'ladzh', ladzh

      end if

      !apply storage to set wtsurf and bctfz as a function of Qstar
      wtsurf = -((1-0.7)*Qstar-0.33*dQdt+38)/(rhoa*cp)
      bctfz  = -((1-0.7)*Qstar-0.33*dQdt+38)/(rhoa*cp)

      ! calc tree cover area - used to apply steady-state boundary condition
      if (myid==0) then
        do n = 1,ntrees
          tr_A = tr_A + (xh(tree(n,2)+1) - xh(tree(n,1)))*dy*(tree(n,4) - tree(n,3) + 1)
        end do
      end if

      ! write updated ground and roof heat fluxes
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

    end subroutine createtrees

    subroutine trees
    use modglobal,  only : ib,ie,jb,je,kb,ke,ih,jh,dzf,xh,dxf,numol,prandtlmol,rlv,cp,tree,&
                           ntrees,ltrees,itot,jtot,cd,ud,lmoist,nsv,dxf,dy,dzf,dzfi,zf,dy,zh,&
                           ltempeq,pref0,r_s,lad,rhoa,ntree_max,lsize,Qstar,dQdt,BCxs,tr_A,&
                           rslabs,kmax,rv,rd
    use modfields,  only : um,vm,wm,thlm,qtm,svp,up,vp,wp,thlp,qtp,svm,Rn,qc,qa,thlm,clai,&
                           ladzf,ladzh,IIcs,tr_omega,&
                           tr_u,tr_v,tr_w,tr_qt,tr_qtR,tr_qtA,tr_thl,tr_sv,thlpcar
    use modmpi,     only : myidx,myidy,mpi_sum,mpierr,comm3d,mpierr,my_real,nprocx,nprocy
    use modsurfdata,only : wtsurf,wttop,wqtop
    use modibmdata, only : bctfz
    use modsubgriddata, only : ekh
    implicit none
    integer :: i,j,k,n,m,il,iu,jl,ju,kl,ku
    real :: e_sat,e_vap,qh,qe,shade,r_a,s,D,numoli,gam,Rq,qhmin,qhmax,omega
    real :: Vq_dum,Vq_dum2,Vq_dum3,VT_dum,VT_dum2,VT_dum3,wTbot_dum,wTbot_dum2,wTbot,wTbot_dum3

    numoli = 1/numol
    gam = (cp*pref0*rv)/(rlv*rd)
    qhmin=0.;qhmax=0.

    if (ltrees .eqv. .false.) return

    ! dummy variables used for steady-state BC calculation
    VT_dum2 = 0.
    Vq_dum2 = 0.
    wTbot_dum2 = 0.

    ! loop over all tree canopies
    do n = 1,ntrees
 
      ! drag in z-direction
      il = tree(n,1) - myidx*itot/nprocx
      iu = tree(n,2) - myidx*itot/nprocx
      kl = tree(n,5)
      ku = tree(n,6) + 1
      jl = tree(n,3) - myidy*jtot/nprocy
      ju = tree(n,4) - myidy*jtot/nprocy
      if (iu < ib .or. il > ie .or. ju < jb .or. jl > je) then
        cycle
      else
        if (iu > ie) iu=ie
        if (il < ib) il=ib 
        if (ju > je) ju=je
        if (jl < jb) jl=jb  

        do k = kl,ku 
          do j = jl,ju
            do i = il,iu
              tr_w(i,j,k) = - cd * ladzh(ntree_max-(ku-k)+1) * wm(i,j,k) * &
                            sqrt( wm(i,j,k)**2 &
                            + (0.25*(um(i,j,k) + um(i+1,j,k) + um(i,j,k-1) + um(i+1,j,k-1)))**2 &
                            + (0.25*(vm(i,j,k) + vm(i,j+1,k) + vm(i,j,k-1) + vm(i,j+1,k-1)))**2 )
            end do
          end do
        end do
      end if

      ! drag in y-direction
      il = tree(n,1) - myidx*itot/nprocx
      iu = tree(n,2) - myidx*itot/nprocx
      kl = tree(n,5)
      ku = tree(n,6)
      jl = tree(n,3) - myidy*jtot/nprocy
      ju = tree(n,4) + 1 - myidy*jtot/nprocy
      if (iu < ib .or. il > ie .or. ju < jb .or. jl > je) then
        cycle
      else
        if (iu > ie) iu=ie
        if (il < ib) il=ib 
        if (ju > je) ju=je
        if (jl < jb) jl=jb 

        do k = kl,ku
          do j = jl,ju
            do i = il,iu
              tr_v(i,j,k) = - cd * ladzf(ntree_max-(ku-k)) * vm(i,j,k) * &
                           sqrt( vm(i,j,k)**2 &
                          + (0.25*(um(i,j,k) + um(i+1,j,k) + um(i,j-1,k) + um(i+1,j-1,k)))**2 &
                          + (0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i,j-1,k) + wm(i,j-1,k+1)))**2 )
            end do
          end do
        end do
      end if

      ! drag in x-direction
      il = tree(n,1) - myidx*itot/nprocx
      iu = tree(n,2) + 1 - myidx*itot/nprocx
      kl = tree(n,5)
      ku = tree(n,6)
      jl = tree(n,3) - myidy*jtot/nprocy
      ju = tree(n,4) - myidy*jtot/nprocy
      if (iu < ib .or. il > ie .or. ju < jb .or. jl > je) then
        cycle
      else
        if (iu > ie) iu=ie
        if (il < ib) il=ib 
        if (ju > je) ju=je
        if (jl < jb) jl=jb 
 
        do k = kl,ku
          do j = jl,ju
            do i = il,iu
              tr_u(i,j,k) = - cd * ladzf(ntree_max-(ku-k)) * um(i,j,k) * &
                           sqrt( um(i,j,k)**2 &
                          + (0.25*(vm(i,j,k) + vm(i,j+1,k) + vm(i-1,j,k) + vm(i-1,j+1,k)))**2 &
                          + (0.25*(wm(i,j,k) + wm(i,j,k+1) + wm(i-1,j,k) + wm(i-1,j,k+1)))**2 )
            end do
          end do
        end do
      end if

      ! scalar volumetric sources/sinks

      ! Canopy EB
      if (lmoist .and. ltempeq) then

        il = tree(n,1) - myidx*itot/nprocx
        iu = tree(n,2) - myidx*itot/nprocx
        kl = tree(n,5)
        ku = tree(n,6)
        jl = tree(n,3) - myidy*jtot/nprocy
        ju = tree(n,4) - myidy*jtot/nprocy
        if (iu < ib .or. il > ie .or. ju < jb .or. jl > je) then
          cycle
        else
          if (iu > ie) iu=ie
          if (il < ib) il=ib 
          if (ju > je) ju=je
          if (jl < jb) jl=jb 
       
          do k = kl,ku
            do j = jl,ju
              do i = il,iu

                 ! psychometrics
                 ! saturation vapour pressure pressure
                 e_sat = 610.8*exp((17.27*(thlm(i,j,k)-273.15))/(thlm(i,j,k)-35.85))

                 ! water vapour partial pressure
                 e_vap = (qtm(i,j,k) * pref0) / (0.378 * qtm(i,j,k) + 0.622)
                 ! qtcheck.m and https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity

                 ! vapour pressure deficit
                 D = max(e_sat - e_vap,0.)

                 ! output warning - cannot physically occur
                 if (e_sat<e_vap) then
                   write(*,*) 'D, e_sat, e_vap, thlm, qtm', D, e_sat, e_vap, thlm(i,j,k), qtm(i,j,k)
                 end if

                 ! slope of the curve relating saturation vapour pressure to temperature
                 s = (4098*e_sat)/((thlm(i,j,k)-35.85)**2)

                 ! aerodynamic resistance
                 r_a = 130*sqrt(lsize/(sqrt((0.5*(um(i,j,k)+um(i+1,j,k)))**2+(0.5*(vm(i,j,k)+vm(i,j+1,k)))**2+(0.5*(wm(i,j,k)+wm(i,j,k+1)))**2)))

                 ! decoupling factor
                 omega = 1/(1 + 2*(gam/(s+2*gam)) * ((r_s)/r_a) )

                 ! latent heat
                 qe = omega*(s/(s+2*gam))*(qa(ntree_max-(ku-k))/(dzf(k)*ladzf(ntree_max-(ku-k)))) + (1-omega)*(1/(gam*(r_s)))*rhoa*cp*D

                 ! sensible heat
                 qh = qa(ntree_max-(ku-k))/(dzf(k)*ladzf(ntree_max-(ku-k))) - qe

                 ! volumetric sinks/source of specific humidity and temp
                 tr_qt(i,j,k) = ladzf(ntree_max-(ku-k))*qe/(rhoa*rlv)
                 tr_qtR(i,j,k) = ladzf(ntree_max-(ku-k))*( omega*(s/(s+2*gam))*(qa(ntree_max-(ku-k))/(dzf(k)*ladzf(ntree_max-(ku-k)))))/(rhoa*rlv)
                 tr_qtA(i,j,k) = ladzf(ntree_max-(ku-k))*((1-omega)*(1/(gam*(r_s)))*rhoa*cp*D)/(rhoa*rlv)
                 tr_thl(i,j,k) = ladzf(ntree_max-(ku-k))*qh/(rhoa*cp)

                 tr_omega(i,j,k) = omega

                 ! add to rhs
                 qtp(i,j,k) = qtp(i,j,k) + tr_qt(i,j,k) !ladz(ntree_max-(ku-k))*qe/(rhoa*rlv)
                 thlp(i,j,k) = thlp(i,j,k) + tr_thl(i,j,k) !ladz(ntree_max-(ku-k))*qh/(rhoa*cp)

                 ! fill dummy variables for steady-state BCs
                 VT_dum2 = VT_dum2 - tr_thl(i,j,k)*dzf(k)*dxf(i)*dy
                 Vq_dum2 = Vq_dum2 - tr_qt(i,j,k)*dzf(k)*dxf(i)*dy

              end do
            end do
          end do

          ! fraction of total net radiation reaching the bottom of the tree
          Rq = Rn(ntree_max+1-(ku+1-kl)) / Qstar

          ! sensible heat flux from shaded surfaces
          shade = ( (1 - 0.7 ) * Rn(ntree_max+1-(ku+1-kl) ) - 0.33*Rq*dQdt + Rq*38 )/ (rhoa*cp)

          ! dummy var for steady-state BC
          wTbot_dum2 = wTbot_dum2 - shade*(xh(iu+1)-xh(il))*dy*(ju-jl+1)

          ! overwrite standard heat flux defined in modibm in shaded positions
          thlp(il:iu,jl:ju,kb+1) = thlp(il:iu,jl:ju,kb+1) + (bctfz + shade) * dzfi(kb+1)
 
         end if

      end if !lmoist and ltempeq

      ! scalar deposition
      if (nsv>0) then

        ! specific to infinite canyon study (tg3315 thesis - Chapter 5)
        if (BCxs==2) then ! no deposition effect from first or last tree to avoid interaction with BCs
          if (n==1) cycle
          if (n==ntrees) cycle
        end if
      
        il = tree(n,1) - myidx*itot/nprocx
        iu = tree(n,2) - myidx*itot/nprocx
        kl = tree(n,5)
        ku = tree(n,6)
        jl = tree(n,3) - myidy*jtot/nprocy
        ju = tree(n,4) - myidy*jtot/nprocy
        if (iu < ib .or. il > ie .or. ju < jb .or. jl > je) then
          cycle
        else
          if (iu > ie) iu=ie
          if (il < ib) il=ib 
          if (ju > je) ju=je
          if (jl < jb) jl=jb 

          do m = 1,nsv ! define this for scalar variables that are deposited
            do k = kl,ku
              do j = jl,ju
                do i = il,iu
                  tr_sv(i,j,k,m) = - svm(i,j,k,m) * ladzf(ntree_max-(ku-k)) * ud
                  svp(i,j,k,m) = svp(i,j,k,m) + tr_sv(i,j,k,m)
                end do
              end do
            end do
          end do
        end if
      end if !nsv

    end do ! ntrees

    ! define these outside loop to avoid double counting cell faces
    wp(ib:ie,jb:je,kb:ke) = wp(ib:ie,jb:je,kb:ke) + tr_w(ib:ie,jb:je,kb:ke)
    vp(ib:ie,jb:je,kb:ke) = vp(ib:ie,jb:je,kb:ke) + tr_v(ib:ie,jb:je,kb:ke)
    up(ib:ie,jb:je,kb:ke) = up(ib:ie,jb:je,kb:ke) + tr_u(ib:ie,jb:je,kb:ke)

    ! steady-state BC calc.
    if (lmoist .and. ltempeq) then

      ! dummy vars - not all used?
      wTbot_dum3 = 0.
      wTbot_dum = 0.
      wTbot = 0.
      VT_dum3 = 0.
      Vq_dum3 = 0.
      VT_dum = 0.
      Vq_dum = 0.

      ! sum dummy variables assigned in main loop
      call MPI_ALLREDUCE(wTbot_dum2,   wTbot_dum3,1,MY_REAL,MPI_SUM,comm3d,mpierr)    
      call MPI_ALLREDUCE(VT_dum2, VT_dum3  ,1,MY_REAL,MPI_SUM,comm3d,mpierr)
      call MPI_ALLREDUCE(Vq_dum2, Vq_dum3  ,1,MY_REAL,MPI_SUM,comm3d,mpierr)

      ! not necessary
      wTbot_dum = wTbot_dum3
      Vq_dum = Vq_dum3
      VT_dum = VT_dum3

      ! total flux of temp from ground and roofs accounting for shading
      wTbot = wTbot_dum + ((xh(ie+1)-xh(ib))*jtot*dy - tr_A)*bctfz

      ! temp flux at top equal to half of the total flux of heat (ground, roofs and canopies) - consistent with tg3315 thesis Chapter 4 - method iii
      wttop = 0.5*(wTbot + VT_dum)/((xh(ie+1)-xh(ib))*jtot*dy)

      ! volumetric sink term throughout domain to account for the remaining heating - ""
      thlpcar = wttop/((sum(real(IIcs(kb+1:ke))/rslabs)/real(kmax-1.))*(zh(ke+1)-zh(kb+1)))

      ! top flux of qt equal to total from canopies - therefore uniform flux profile
      wqtop = (Vq_dum ) / ((xh(ie+1)-xh(ib))*jtot*dy)

    end if

    end subroutine trees

  end module modtrees
