!> \file modetrees.f90
!!tg3315, ns4513, 20 Mar 2017 

!> Input trees into DALES model.

module modtrees
implicit none
save

contains
    subroutine createtrees
    use modglobal,  only : ltrees,ntrees,tree,cexpnr,ifinput,zh,zf,dzh,dzfi,dzhi,dzf,Rshade,sun,&
                           decay,ladt,kb,ke,cp
    use modfields,  only : u0,v0,w0,thl0,qt0,svp,up,vp,wp,thlp,qtp,Rn,lad,clai,qc
    use modmpi,     only : myid,MPI_INTEGER,comm3d,mpierr,MY_REAL
    use modsurfdata,only : wtsurf
    use modibmdata ,only : bctfz

    implicit none
    integer :: n,k,nlad
    real :: Rc,Rc_s
    character(80) chmess

    if (ltrees .eqv. .false.) return

      allocate(tree(ntrees,6)) 

      ! read global trees
      if(myid==0) then
        write(*,*) '1, myid, ntrees, ltrees, cexpnr', myid, ntrees, ltrees, cexpnr
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

      nlad = tree(1,6)+2-tree(1,5)

      ! read leaf area density of trees
      if (myid==0) then

        allocate(ladt(nlad)) 

        if (nlad>0) then
          open (ifinput,file='lad.inp.'//cexpnr)
          read (ifinput,'(a80)') chmess
          do n=1,nlad
            read (ifinput,*) &
                  ladt(n)
            end do

          write (*,*) 'lad'
          do n=1,nlad
            write (1,*) &
                  n , &
                  ladt(n)
          end do
        end if

      ! input into leaf area density array (cell faces)
      lad(tree(1,5):tree(1,6)+1) = ladt !defined at cell faces

      ! calculate cumulative leaf area intensity (cell faces)
      do k = tree(1,5),tree(1,6)
        clai(k) = sum(0.5*(lad(k:tree(1,6))+lad(k+1:tree(1,6)+1))*(dzh(k:tree(1,6))))
      end do

      ! Net radiation at cell faces! !W/m^2
      do k = tree(1,5),tree(1,6)+1
         Rn(k) = 0.3394 * sun * exp(- decay * clai(k))
      enddo
     
      ! divergence of radiative flux at cell centres W/m^3
      do k = tree(1,5),tree(1,6)
        Rc = Rn(k+1)-Rn(k)
        Rc_s = 0.11*Rc-12.3*Rc/Rn(tree(1,6)+1)
        qc(k) = (Rc-Rc_s) * dzhi(k)
      end do

      ! resultant surface flux underneath tree 
      Rshade = ((1-0.11)*Rn(tree(1,5)) + 12.3*Rn(tree(1,5))/Rn(tree(1,6)+1) ) / (1.225*cp)

    end if ! end if myid==0  

    wtsurf = -((1-0.7)*0.3988*sun+38)/(1.225*cp)
    bctfz  = -((1-0.7)*0.3988*sun+38)/(1.225*cp)

    ! broadcast variables
    call MPI_BCAST(lad   ,ke-kb+1,MY_REAL ,0,comm3d,mpierr)
    call MPI_BCAST(Rn ,ke-kb+1,MY_REAL ,0,comm3d,mpierr)
    call MPI_BCAST(qc ,ke-kb+1,MY_REAL ,0,comm3d,mpierr)
    call MPI_BCAST(clai  ,ke-kb+1,MY_REAL ,0,comm3d,mpierr)  
    call MPI_BCAST(Rshade ,1,MY_REAL ,0,comm3d,mpierr)  

    end subroutine createtrees

    subroutine trees
    use modglobal,  only : ib,ie,jb,je,kb,ke,dzf,numol,prandtlmol,rv,fkar,rlv,cp,jtot,tree,&
                           ntrees,ltrees,jtot,cd,ud,lmoist,nsv,dxf,dy,dzf,dzfi,zf,Rshade,&
                           ltempeq,Bowen,lnudge,lscasrcl,tnudge,nnudge
    use modfields,  only : u0,v0,w0,thl0,qt0,svp,up,vp,wp,thlp,qtp,sv0,Rn,qc,lad,thlm,&
                           thl0av,sv0av,qt0av
    use modmpi,     only : myid,nprocs
    use modsurfdata,only : ps,wtsurf
    use modsubgriddata, only : ekh
    implicit none
    integer :: i,j,k,n,m,il,iu,jl,ju,kl,ku,treec
    real :: ares,psat,pact,pdef,etdelta,evtrans,numoli

    numoli = 1/numol
   
    ! Tree drag
    if (ltrees .eqv. .false.) return

    !write(*,*) 'myid, ntrees', myid, ntrees
 
    do n = 1,ntrees
     
    !  write(*,*) 'myid2, tree(n,1), tree(n,2), tree(n,3), tree(n,4), tree(n,5), tree(n,6)', myid, tree(n,1), tree(n,2), tree(n,3), tree(n,4), tree(n,5), tree(n,6)
 
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
              wp(i,j,k) = wp(i,j,k) - cd * lad(k) * w0(i,j,k) * &
                         sqrt( w0(i,j,k)**2 &
                        + (0.25*(u0(i,j,k) + u0(i+1,j,k) + u0(i,j,k-1) + u0(i+1,j,k-1)))**2 &
                        + (0.25*(v0(i,j,k) + v0(i,j+1,k) + v0(i,j,k-1) + v0(i,j+1,k-1)))**2 )
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
              vp(i,j,k) = vp(i,j,k) - cd * lad(k) * v0(i,j,k) * &
                           sqrt( v0(i,j,k)**2 &
                          + (0.25*(u0(i,j,k) + u0(i+1,j,k) + u0(i,j-1,k) + u0(i+1,j-1,k)))**2 &
                          + (0.25*(w0(i,j,k) + w0(i,j,k+1) + w0(i,j-1,k) + w0(i,j-1,k+1)))**2 )
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
              up(i,j,k) = up(i,j,k) - cd * lad(k) * u0(i,j,k) * &
                           sqrt( u0(i,j,k)**2 &
                          + (0.25*(v0(i,j,k) + v0(i,j+1,k) + v0(i-1,j,k) + v0(i-1,j+1,k)))**2 &
                          + (0.25*(w0(i,j,k) + w0(i,j,k+1) + w0(i-1,j,k) + w0(i-1,j,k+1)))**2 )
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
        treec = (kl+ku) / 2
        if (ju < jb .or. jl > je) then
          cycle
        else
          if (ju > je) ju=je
          if (jl < jb) jl=jb 
       
          do k = kl,ku
            do j = jl,ju
              do i = il,iu
                 ! BOWEN RATIO
                 qtp(i,j,k) = qtp(i,j,k) + qc(k) / ( rlv * 1.225 * (1. + Bowen) )
                 thlp(i,j,k) = thlp(i,j,k) + qc(k) / ( (1. + 1. / Bowen) * 1.225 * cp )
              end do
            end do
          end do
        end if
      end if !lmoist

      ! Temp only
      if (ltempeq .and. (.not. lmoist)) then
    
        il = tree(n,1)
        iu = tree(n,2)
        kl = tree(n,5)
        ku = tree(n,6)
        jl = tree(n,3) - myid*jtot/nprocs
        ju = tree(n,4) - myid*jtot/nprocs
        treec = (kl+ku) / 2
        if (ju < jb .or. jl > je) then
          cycle
        else
          if (ju > je) ju=je
          if (jl < jb) jl=jb 
       
          do k = kl,ku
            do j = jl,ju
              do i = il,iu
                  thlp(i,j,k) = thlp(i,j,k) + qc(k) / ( 1.225 * cp )
              end do
            end do
          end do
        end if
      end if

      if (ltempeq) then
        thlp(il:iu,jl:ju,kb) = thlp(il:iu,jl:ju,kb) + (wtsurf + Rshade)*dzfi(kb)
      end if

      ! scalar deposition
      if (nsv>0) then
        il = tree(n,1)
        iu = tree(n,2)
        kl = tree(n,5)
        ku = tree(n,6)
        jl = tree(n,3) - myid*jtot/nprocs
        ju = tree(n,4) - myid*jtot/nprocs
        treec = (kl+ku) / 2
        if (ju < jb .or. jl > je) then
          cycle
        else
          if (ju > je) ju=je
          if (jl < jb) jl=jb 

          do m = 1,nsv
            do k = kl,ku
              do j = jl,ju
                do i = il,iu
                  svp(i,j,k,m) = svp(i,j,k,m) - sv0(i,j,k,m) * lad(k) * ud
                end do
              end do
            end do
          end do
        end if
      end if !nsv

    end do ! ntrees

    end subroutine trees

  end module modtrees
