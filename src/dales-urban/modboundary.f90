!> \file modboundary.f90
!>
!>
!! \par Revision list
!! All boundary conditions are in this file, except for immersed boundaries.
!! \par Authors
!! Jasper Tomas, March 31 2014
!!
module modboundary


  implicit none
  save
  private
  public :: initboundary, boundary, bottom, exitboundary,grwdamp, ksp,tqaver,cyclich,&
       bcp,bcpup,wallaw,wallawt,masscorr,&
       closurebc
  integer :: ksp = -1                 !<    lowest level of sponge layer
  real,allocatable :: tsc(:)          !<   damping coefficients to be used in grwdamp.
  real :: rnu0 = 2.75e-3
contains
  !>
  !! Initializing Boundary; specifically the sponge layer
  !>
  subroutine initboundary
    use modglobal, only : ib,kb,ke,kh,kmax,pi,zf,iplane
    use modinletdata, only : irecy
    implicit none

    real    :: zspb, zspt
    integer :: k
    allocate(tsc(kb:ke+kh))
    ! Sponge layer
    if (ksp==-1) then
       !      ksp  = min(3*kmax/4,kmax - 15)
       ksp  = (kb-1)+min(3*kmax/4,kmax - 15)
    end if

    zspb    = zf(ksp)
    zspt    = zf(ke)

    tsc(kb:ksp-1) = 0.0
    do k=ksp,ke
       tsc(k) = rnu0*sin(0.5*pi*(zf(k)-zspb)/(zspt-zspb))**2
    end do
    tsc(ke+1)=tsc(ke)
    irecy = ib + iplane
  end subroutine initboundary

  !>
  !! Execute boundary conditions
  subroutine boundary

    use modglobal, only : ib,ie,ih,jb,je,jgb,jge,jh,kb,ke,kh,linoutflow,dzf,zh,dy,&
         timee,ltempeq,lmoist
    use modfields, only : u0,w0,um,uout,uouttot
    use modmpi,    only : myid,slabsum
    use modinlet,  only : inletgen,inletgennotemp
    use modinletdata, only : irecy,ubulk,iangle
    use modsurface, only : getobl
    implicit none
    real massin, massout,massout2, masstopke,masstopke1
    integer i,k  

    call getobl

    if (linoutflow == .true.) then
!       uouttot = ubulk
uouttot = cos(iangle)*ubulk     
  if (ltempeq ==.true.) then
          call inletgen
       else 
          call inletgennotemp
       end if
       call topm
       call iolet
       call cyclicj 
    else
       call cyclicm
       call cyclich
       call topm
    endif
    call toph
    call bottom

  end subroutine boundary
  !> Cleans up after the run
  subroutine exitboundary
    implicit none
    deallocate(tsc)
  end subroutine exitboundary


  subroutine closurebc
    use modsubgriddata, only : ekm,ekh
    use modglobal,      only : ib,ie,jb,je,kb,ke,ih,jh,kh,numol,prandtlmoli,linoutflow,&
         linletgen, lzerogradtop
    use modmpi,         only : excjs
    integer i,j


    ! Top and bottom
    if (lzerogradtop==.true. .or. linletgen==1 .or. linletgen==2) then
       do j=jb-1,je+1
          do i=ib-1,ie+1
             ekm(i,j,ke+1)  = ekm(i,j,ke)                          ! zero-gradient top wall
             ekh(i,j,ke+1)  = ekh(i,j,ke)                          ! zero-gradient top wall
             ekm(i,j,kb-1)  = 2.*numol -ekm(i,j,kb)                ! no-slip lower wall
             ekh(i,j,kb-1) = (2.*numol*prandtlmoli) -ekh(i,j,kb)   ! no-slip lower wall
          end do
       end do
    else ! no slip top wall
       do j=jb-1,je+1
          do i=ib-1,ie+1
             ekm(i,j,ke+1)  = 2.*numol - ekm(i,j,ke)               ! no-slip top wall
             ekh(i,j,ke+1)  = (2.*numol*prandtlmoli) -ekh(i,j,ke)  ! no-slip top wall
             ekm(i,j,kb-1)  = 2.*numol -ekm(i,j,kb)                ! no-slip lower wall
             ekh(i,j,kb-1) = (2.*numol*prandtlmoli) -ekh(i,j,kb)   ! no-slip lower wall
          end do
       end do
    end if

    ! horizontal BC's
    if (linoutflow == .true.) then      ! inflow/outflow
       ekm(ib-1, :,:) = ekm(ib,:,:)
       ekm(ie+1,:,:) = ekm(ie, :,:)
       ekh(ib-1, :,:) = ekh(ib,:,:)
       ekh(ie+1,:,:) = ekh(ie, :,:)
    else
       ekm(ib-1, :,:) = ekm(ie,:,:)      ! periodic
       ekm(ie+1,:,:) = ekm(ib, :,:)
       ekh(ib-1, :,:) = ekh(ie,:,:) 
       ekh(ie+1,:,:) = ekh(ib, :,:)
    end if

    call excjs( ekm           , ib,ie,jb,je,kb-kh,ke+kh,ih,jh)
    call excjs( ekh           , ib,ie,jb,je,kb-kh,ke+kh,ih,jh)

  end subroutine closurebc



  !> Sets lateral periodic boundary conditions for the scalars
  subroutine cyclich
    use modglobal, only : ib,ie,jb,je,ih,jh,kb,ke,kh,nsv,dt,lscalinout,lscalrec,lmoistinout,ltempinout,rk3step,dxhi,ltempeq,&
                          ihc,jhc,khc
    use modfields, only : thl0,thlm,qt0,qtm,sv0,svm,svprof,thlprof,qtprof,uouttot
    use modmpi,    only : excjs
    real rk3coef
    integer k,n,m


    if (ltempinout == .true.) then
       ! New code for nonperiodic moisture and temperature
         do m = 1,ih
           thl0(ib-m,:,:)  = thl0(ie+1-m,:,:)
           thl0(ie+m,:,:)  = thl0(ib-1+m,:,:)
           thlm(ib-m,:,:)  = thlm(ie+1-m,:,:)
           thlm(ie+m,:,:)  = thlm(ib-1+m,:,:)
         end do
       ! j-direction
       call excjs( thl0           , ib,ie,jb,je,kb,ke+kh,ih,jh)
       call excjs( thlm           , ib,ie,jb,je,kb,ke+kh,ih,jh)
    end if

    if (lmoistinout == .true.) then
           do m = 1,ih
           qt0(ib-m,:,:)   = qt0(ie+1-m,:,:)
           qt0(ie+m,:,:)   = qt0(ib-1+m,:,:)
           qtm(ib-m,:,:)   = qtm(ie+1-m,:,:)
           qtm(ie+m,:,:)   = qtm(ib-1+m,:,:)
          end do
       call excjs( qt0            , ib,ie,jb,je,kb,ke+kh,ih,jh)
       call excjs( qtm            , ib,ie,jb,je,kb,ke+kh,ih,jh)
    end if



    if (lscalinout == .true.) then  ! convective outflow
       rk3coef = dt / (4. - dble(rk3step))
       do n=1,nsv
          do k=kb,ke+1
             sv0(ib-1,:,k,n) = 2*svprof(k,n) - sv0(ib,:,k,n)
             svm(ib-1,:,k,n) = 2*svprof(k,n) - svm(ib,:,k,n)
          end do
          sv0(ie+1,:,:,n)= sv0(ie,:,:,n) - (sv0(ie+1,:,:,n)-sv0(ie,:,:,n))*dxhi(ie+1)*rk3coef*uouttot
          svm(ie+1,:,:,n)= svm(ie,:,:,n) - (svm(ie+1,:,:,n)-svm(ie,:,:,n))*dxhi(ie+1)*rk3coef*uouttot
       enddo
    
    elseif (lscalrec == .true.) then 
      if (nsv>0) then
        rk3coef = dt / (4. - dble(rk3step))
        do m = 1,ihc ! loop over virtual cells
          do n=1,nsv-1
            sv0(ib-m,:,:,n+1) = sv0(ie+1-m,:,:,n)
            sv0(ie+m,:,:,n) = sv0(ib-1+m,:,:,n+1)
            svm(ib-m,:,:,n+1) = svm(ie+1-m,:,:,n)
            svm(ie+m,:,:,n) = svm(ib-1+m,:,:,n+1)
          end do

          ! zero conc. on scalar 1 !tg3315 should be changed to as above in
          ! lscalinout???
          sv0(ib-m,:,:,1) = 0
          svm(ib-m,:,:,1) = 0

          ! DIY outflow BC (advection step as linout) tg3315

          sv0(ie+m,:,:,nsv)=sv0(ie+1-m,:,:,nsv)-(sv0(ie+m,:,:,nsv)-sv0(ie+1-m,:,:,nsv))*dxhi(ie+m)*rk3coef*uouttot
          svm(ie+m,:,:,nsv)=svm(ie+1-m,:,:,nsv)-(svm(ie+m,:,:,nsv)-svm(ie+1-m,:,:,nsv))*dxhi(ie+m)*rk3coef*uouttot
        end do 
      end if

      else  ! 'normal' periodic boundary condition
      ! do m = 1,ih
        do m = 1,ihc
          sv0(ib-m,:,:,:) = sv0(ie+1-m,:,:,:)
          sv0(ie+m,:,:,:) = sv0(ib-1+m,:,:,:)
          svm(ib-m,:,:,:) = svm(ie+1-m,:,:,:)
          svm(ie+m,:,:,:) = svm(ib-1+m,:,:,:)
        end do
       
    end if

    do n=1,nsv
       call excjs( sv0(:,:,:,n) , ib,ie,jb,je,kb-khc,ke+khc,ihc,jhc)
       call excjs( svm(:,:,:,n) , ib,ie,jb,je,kb-khc,ke+khc,ihc,jhc)
    enddo


    return
  end subroutine cyclich

  !>set lateral periodic boundary conditions for momentum
  subroutine cyclicm

    use modglobal, only : ib,ie,jb,je,ih,jh,kb,ke,kh,jmax
    use modfields, only : u0,um,v0,vm,w0,wm,e120,e12m,shear
    use modsubgriddata, only : loneeqn,lsmagorinsky
    use modmpi,    only : excjs

    integer n,m

    do m = 1,ih

       u0(ib-m,:,:)   = u0(ie+1-m,:,:)
       u0(ie+m,:,:)   = u0(ib-1+m,:,:)
       v0(ib-m,:,:)   = v0(ie+1-m,:,:)
       v0(ie+m,:,:)   = v0(ib-1+m,:,:)
       w0(ib-m,:,:)   = w0(ie+1-m,:,:)
       w0(ie+m,:,:)   = w0(ib-1+m,:,:)
       um(ib-m,:,:)   = um(ie+1-m,:,:)
       um(ie+m,:,:)   = um(ib-1+m,:,:)
       vm(ib-m,:,:)   = vm(ie+1-m,:,:)
       vm(ie+m,:,:)   = vm(ib-1+m,:,:)
       wm(ib-m,:,:)   = wm(ie+1-m,:,:)
       wm(ie+m,:,:)   = wm(ib-1+m,:,:)

       e120(ib-m,:,:) = e120(ie+1-m,:,:)
       e120(ie+m,:,:) = e120(ib-1+m,:,:)
       e12m(ib-m,:,:) = e12m(ie+1-m,:,:)
       e12m(ie+m,:,:) = e12m(ib-1+m,:,:)

    end do

    call excjs( u0  , ib,ie,jb,je,kb,ke+kh,ih,jh)
    call excjs( v0  , ib,ie,jb,je,kb,ke+kh,ih,jh)
    call excjs( w0  , ib,ie,jb,je,kb,ke+kh,ih,jh)
    call excjs( um  , ib,ie,jb,je,kb,ke+kh,ih,jh)
    call excjs( vm  , ib,ie,jb,je,kb,ke+kh,ih,jh)
    call excjs( wm  , ib,ie,jb,je,kb,ke+kh,ih,jh)

    if (loneeqn==.true.) then
         e120(ib-m,:,:) = e120(ie+1-m,:,:)
         e120(ie+m,:,:) = e120(ib-1+m,:,:)
         e12m(ib-m,:,:) = e12m(ie+1-m,:,:)
         e12m(ie+m,:,:) = e12m(ib-1+m,:,:)
       call excjs( e120, ib,ie,jb,je,kb,ke+kh,ih,jh)
       call excjs( e12m, ib,ie,jb,je,kb,ke+kh,ih,jh)
       ! exchange shear components between processors
       do n=1,12 ! for all 12 components
          call excjs( shear(:,:,:,n) , ib,ie,jb,je,kb,ke,0,1)
       end do
    end if

    if (lsmagorinsky==.true.) then
       ! exchange shear components between processors
       do n=1,12 ! for all 12 components
          call excjs( shear(:,:,:,n) , ib,ie,jb,je,kb,ke,0,1)
       end do
    end if

    return
  end subroutine cyclicm

  !>set lateral periodic boundary conditions in j-direction (in case of inflow/outflow BC's in i-direction)
  subroutine cyclicj

    use modglobal, only : ib,ie,jb,je,ih,jh,kb,ke,kh,nsv,jmax,ltempeq,lmoist,ihc,jhc,khc
    use modfields, only : u0,um,v0,vm,w0,wm,e120,e12m,thl0,thlm,qt0,qtm,sv0,svm,shear
    use modsubgriddata, only : loneeqn,lsmagorinsky
    use modmpi,    only : excjs

    integer n

    ! Momentum
    call excjs( u0  , ib,ie,jb,je,kb,ke+kh,ih,jh)
    call excjs( v0  , ib,ie,jb,je,kb,ke+kh,ih,jh)
    call excjs( w0  , ib,ie,jb,je,kb,ke+kh,ih,jh)
    call excjs( um  , ib,ie,jb,je,kb,ke+kh,ih,jh)
    call excjs( vm  , ib,ie,jb,je,kb,ke+kh,ih,jh)
    call excjs( wm  , ib,ie,jb,je,kb,ke+kh,ih,jh)

    ! Heat
    if (ltempeq == .true.) then
       call excjs( thl0           , ib,ie,jb,je,kb,ke+kh,ih,jh)
       call excjs( thlm           , ib,ie,jb,je,kb,ke+kh,ih,jh)
    end if
    if (lmoist == .true.) then
       call excjs( qt0            , ib,ie,jb,je,kb,ke+kh,ih,jh)
       call excjs( qtm            , ib,ie,jb,je,kb,ke+kh,ih,jh)
    end if

    do n=1,nsv
       call excjs( sv0(:,:,:,n) , ib,ie,jb,je,kb-khc,ke+khc,ihc,jhc)
       call excjs( svm(:,:,:,n) , ib,ie,jb,je,kb-khc,ke+khc,ihc,jhc)
    enddo

    if (loneeqn==.true.) then
       call excjs( e120, ib,ie,jb,je,kb,ke+kh,ih,jh)
       call excjs( e12m, ib,ie,jb,je,kb,ke+kh,ih,jh)
       ! exchange shear components between processors
       do n=1,12 ! for all 12 components
          !!  call excjs( shear(:,:,:,n) , ib,ie,jb,je,kb,ke,0,jmax)
          call excjs( shear(:,:,:,n) , ib,ie,jb,je,kb,ke,0,1)
       end do
    end if

    if (lsmagorinsky==.true.) then  ! needed for wall-damping (lsmagorinsky and loneeqn)
       ! exchange shear components between processors
       do n=1,12 ! for all 12 components
          call excjs( shear(:,:,:,n) , ib,ie,jb,je,kb,ke,0,1)
       end do
    end if

    return
  end subroutine cyclicj

  !>set inlet and outlet boundary conditions in i-direction
  subroutine iolet

    use modglobal, only : dxhi,dxhci,xh,zh,ib,ie,jb,je,ih,jh,kb,ke,kh,nsv,rk3step,dt,linletgen,ltempeq,lmoist,ihc
    use modfields, only : u0,um,v0,vm,w0,wm,e120,e12m,thl0,thlm,qt0,qtm,sv0,svm,uprof,vprof,e12prof,thlprof,&
         qtprof,svprof,uouttot,wouttot
    use modmpi,    only : excjs,myid
    use modinletdata, only : u0inletbcold,v0inletbcold,w0inletbcold,uminletbc,vminletbc,wminletbc,totaluold,&
         t0inletbcold,tminletbc

    real rk3coef,uouttotold

    integer n,i,j,k,m


    rk3coef = dt / (4. - dble(rk3step))


    ! Inlet boundary is located at ib (not ib-1)!

    ! Inlet
    if (linletgen == 1 .or. linletgen ==2) then
       do j=jb,je
          do k=kb,ke
             u0(ib,j,k)=u0inletbcold(j,k)
             um(ib,j,k)=uminletbc(j,k)
             u0(ib-1,j,k)=2*u0(ib,j,k)-u0(ib+1,j,k)
             um(ib-1,j,k)=2*um(ib,j,k)-um(ib+1,j,k)

             v0(ib-1,j,k)   = v0inletbcold(j,k) 
             vm(ib-1,j,k)   = vminletbc(j,k)      

             ! to be changed in the future: e12 should be taken from recycle plane!
             e120(ib-1,j,k) = e120(ib,j,k)      ! extrapolate e12 from interior
             e12m(ib-1,j,k) = e12m(ib,j,k)      ! extrapolate e12 from interior
             
             do n=1,nsv
                do m = 1,ihc
                   sv0(ib-m,j,k,n) = 2*svprof(k,n) - sv0(ib+(m-1),j,k,n)
                   svm(ib-m,j,k,n) = 2*svprof(k,n) - svm(ib+(m-1),j,k,n)
                enddo
             enddo
          end do
          do k=kb,ke+1
             w0(ib-1,j,k)   = w0inletbcold(j,k) 
             wm(ib-1,j,k)   = wminletbc(j,k)
          end do
       end do

       ! Heat
       if (ltempeq ==.true.) then
          do k=kb,ke
             do j=jb,je
                thl0(ib-1,j,k) = t0inletbcold(j,k)
                thlm(ib-1,j,k) = tminletbc(j,k)
             end do
          end do
       end if


       if (lmoist ==.true.) then
          do k=kb,ke
             do j=jb,je
                qt0(ib-1,j,k) = 2*qtprof(k) - qt0(ib,j,k)  !watch!
                qtm(ib-1,j,k) = 2*qtprof(k) - qtm(ib,j,k)
             end do
          end do
       end if


    else  ! (if linetgen==.false.)

       do j=jb-1,je+1
          do k=kb,ke+1
             ! Momentum
             u0(ib,j,k)=uprof(k)
             um(ib,j,k)=uprof(k)

             v0(ib-1,j,k)   = 2*vprof(k) - v0(ib,j,k)      ! (v(ib)+v(ib-1))/2 = vprof
             vm(ib-1,j,k)   = 2*vprof(k) - vm(ib,j,k)      ! (v(ib)+v(ib-1))/2 = vprof

             e120(ib-1,j,k) = 2*e12prof(k) - e120(ib,j,k)      ! (e12(ib)+e12(ib-1))/2=e12prof
             e12m(ib-1,j,k) = 2*e12prof(k) - e12m(ib,j,k)      ! (e12(ib)+e12(ib-1))/2=e12prof
             
             ! commented by tg3315 !undone !i think should only be for nsv =1.
             !do n=1,nsv
             !   do m = 1,ihc
             !      sv0(ib-m,j,k,n) = 2*svprof(k,n) - sv0(ib+(m-1),j,k,n)
             !      svm(ib-m,j,k,n) = 2*svprof(k,n) - svm(ib+(m-1),j,k,n)
             !   enddo
             !enddo

             ! added tg3315
             do m = 1,ihc
                   sv0(ib-m,j,k,1) = 2*svprof(k,1) - sv0(ib+(m-1),j,k,1)
                   svm(ib-m,j,k,1) = 2*svprof(k,1) - svm(ib+(m-1),j,k,1)
             enddo



          enddo
       enddo

       ! Heat
       if (ltempeq == .true.) then
          do j=jb-1,je+1
             do k=kb,ke+1
                thl0(ib-1,j,k) = 2*thlprof(k) - thl0(ib,j,k) 
                thlm(ib-1,j,k) = 2*thlprof(k) - thlm(ib,j,k) 
             end do
          end do
       end if

       if (lmoist == .true.) then
          do j=jb-1,je+1
             do k=kb,ke+1
                qt0(ib-1,j,k) = 2*qtprof(k) - qt0(ib,j,k) 
                qtm(ib-1,j,k) = 2*qtprof(k) - qtm(ib,j,k) 
             end do
          end do
       end if


       u0(ib-1,:,:)   = 2*u0(ib,:,:)-u0(ib+1,:,:)  ! (u(ib+1)+u(ib-1))/2 = u(ib)   
       um(ib-1,:,:)   = 2*um(ib,:,:)-um(ib+1,:,:)    ! (u(ib+1)+u(ib-1))/2 = u(ib)
       w0(ib-1,:,:)   = -w0(ib,:,:)                  ! (w(ib)+w(ib-1))/2 = 0
       wm(ib-1,:,:)   = -wm(ib,:,:)
    end if ! linletgen==1 .or. linletgen==2

    ! Outlet
    ! Momentum
    v0(ie+1,:,:)   = v0(ie,:,:) - (v0(ie+1,:,:)-v0(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
    w0(ie+1,:,:)   = w0(ie,:,:) - (w0(ie+1,:,:)-w0(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
    vm(ie+1,:,:)   = vm(ie,:,:) - (vm(ie+1,:,:)-vm(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
    wm(ie+1,:,:)   = wm(ie,:,:) - (wm(ie+1,:,:)-wm(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
    e120(ie+1,:,:) = e120(ie,:,:) - (e120(ie+1,:,:)-e120(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
    e12m(ie+1,:,:) = e12m(ie,:,:) - (e12m(ie+1,:,:)-e12m(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot

    ! Heat
    if (ltempeq == .true.) then
       thl0(ie+1,:,:) = thl0(ie,:,:) - (thl0(ie+1,:,:)-thl0(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
       thlm(ie+1,:,:) = thlm(ie,:,:) - (thlm(ie+1,:,:)-thlm(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
    end if


    if (lmoist == .true.) then
       qt0(ie+1,:,:)  = qt0(ie,:,:) - (qt0(ie+1,:,:)-qt0(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
       qtm(ie+1,:,:)  = qtm(ie,:,:) - (qtm(ie+1,:,:)-qtm(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
    end if


    ! commented by tg3315 !undone !changed dxhi to dxhci!?
    !do n=1,nsv
    !   sv0(ie+1,:,:,n)= sv0(ie,:,:,n) - (sv0(ie+1,:,:,n)-sv0(ie,:,:,n))*dxhi(ie+1)*rk3coef*uouttot
    !   svm(ie+1,:,:,n)= svm(ie,:,:,n) - (svm(ie+1,:,:,n)-svm(ie,:,:,n))*dxhi(ie+1)*rk3coef*uouttot
    !end do

    !added tg3315
    sv0(ie+1,:,:,nsv)= sv0(ie,:,:,nsv) - (sv0(ie+1,:,:,nsv)-sv0(ie,:,:,nsv))*dxhci(ie+1)*rk3coef*uouttot
    svm(ie+1,:,:,nsv)= svm(ie,:,:,nsv) - (svm(ie+1,:,:,nsv)-svm(ie,:,:,nsv))*dxhci(ie+1)*rk3coef*uouttot

    return
  end subroutine iolet

  !> correct the u-velocity to get correct mass flow rate
  subroutine masscorr

    use modglobal, only : ib,ie,jb,je,ih,jh,kb,ke,kh,jgb,jge,dzf,zh,dy,dt,rk3step,massflowrate,&
         jmax,libm,dt,rk3step,linoutflow,lmassflowr
    use modfields, only : um,up,uout,uouttot
    use modmpi,    only : slabsum,myid

    real, dimension(kb:ke)             :: uoutold
    real udef,rk3coef,rk3coefi,massflowrateold
    integer i,j,k

    ! need to add .or. lscalrec??? tg3315
    if (linoutflow == .false. .and. lmassflowr == .true.) then
       rk3coef = dt / (4. - dble(rk3step))
       rk3coefi = 1 / rk3coef

       uout = 0.
       uoutold = 0.
       call slabsum(uout   ,kb,ke, up  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ie,ie,jb,je,kb,ke) ! determine horizontal (j) average outflow velocity diff
       call slabsum(uoutold,kb,ke, um  ,ib-ih,ie+ih,jb-jh,je+jh,kb-kh,ke+kh,ie,ie,jb,je,kb,ke) ! determine horizontal (j) average outflow velocity old

       do k=kb,ke
          uout(k)    = rk3coef*uout(k)   *dzf(k)*dy*1000.  ! mass flow rate through each slab (density = 1000 kg/m3)
          uoutold(k) =         uoutold(k)*dzf(k)*dy*1000.  ! mass flow rate through each slab (density = 1000 kg/m3) (previous time step)
       end do
       uouttot         = sum(uout(kb:ke))                 ! mass flow rate (at outlet)
       massflowrateold = sum(uoutold(kb:ke))              ! mass flow rate (at outlet) (previous time step)
       udef =  (massflowrate - (uouttot + massflowrateold))/(((jge-jgb+1)*dy)*(zh(ke+1)-zh(kb))*1000.)   !udef=massdef/(Area*density)
       do k = kb,ke
          do j = jb,je
             do i = ib,ie
                up(i,j,k)  = up(i,j,k)  + udef*rk3coefi
             end do
          end do
       end do
    end if

  end subroutine masscorr

  !>set boundary conditions pup,pvp,pwp in subroutine fillps in modpois.f90
  subroutine bcpup(pup,pvp,pwp,rk3coef)

    use modglobal, only : ib,ie,jb,je,ih,jh,kb,ke,kh,linoutflow,dxfi,linletgen,&
         Uinf,libm,jmax
    use modfields, only : pres0,up,vp,wp,um,w0,u0,uouttot
    use modmpi, only : excjs,myid
    use modinletdata,only : irecy,u0inletbc,ddispdx

    real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh), intent(inout) :: pup 
    real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh), intent(inout) :: pvp 
    real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh), intent(inout) :: pwp 

    real, intent(in) :: rk3coef
    real rk3coefi

    integer i,j,k

    rk3coefi = 1. / rk3coef
    ! again or scalrec??? tg3315
    if (linoutflow == .true.) then
       if (linletgen == 1 .or. linletgen == 2) then
          do j=jb,je
             do i=ib,ie
                pwp(i,j,kb)  = 0.
                pwp(i,j,ke+kh)= (Uinf*ddispdx ) *rk3coefi
             end do
          end do
          do k=kb,ke
             do j=jb,je
                pup(ie+1,j,k) = - (u0(ie+1,j,k)-u0(ie,j,k))*dxfi(ie)*uouttot + um(ie+1,j,k)*rk3coefi   ! du/dt +u*du/dx=0 -> pup(i)=um(i)/rk3coef -um(i)*(um(i)-um(i-1))/dxf(i-1)
                pup(ib,j,k) = u0inletbc(j,k)*rk3coefi
             end do
          end do
       else ! if not linletgen
          do j=jb,je
             do i=ib,ie
                pwp(i,j,kb)  = 0.
                pwp(i,j,ke+kh) = 0.
             end do
          end do
          do k=kb,ke
             do j=jb,je
                pup(ie+1,j,k) = - (u0(ie+1,j,k)-u0(ie,j,k))*dxfi(ie)*uouttot + um(ie+1,j,k)*rk3coefi   ! du/dt +u*du/dx=0 -> pup(i)=um(i)/rk3coef -um(i)*(um(i)-um(i-1))/dxf(i-1)
                pup(ib,j,k)   = pup(ib,j,k)-up(ib,j,k)  ! pup(ib)= up(ib) + um(ib)/rk3coef, where up should be zero!
             end do
          end do
       end if ! inletgen
    else ! if not linoutflow
       do j=jb,je
          do i=ib,ie
             pwp(i,j,kb)  = 0.
             pwp(i,j,ke+kh) = 0.
          end do
       end do
       do k=kb,ke
          do j=jb,je
             pup(ie+1,j,k) = pup(ib,j,k)                         ! cyclic 
          end do
       end do
    endif

    call excjs( pup          , ib,ie,jb,je,kb,ke+kh,ih,jh)  ! cyclic
    call excjs( pvp          , ib,ie,jb,je,kb,ke+kh,ih,jh)  ! cyclic
    call excjs( pwp          , ib,ie,jb,je,kb,ke+kh,ih,jh)  ! cyclic

  end subroutine bcpup


  !>set pressure boundary conditions
  subroutine bcp(p)

    use modglobal, only : ib,ie,jb,je,ih,jh,kb,ke,kh,linoutflow,dxfi,linletgen
    use modfields, only : pres0,up,u0,um,uouttot
    use modmpi, only : excj

    real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh), intent(inout) :: p !< pressure
    integer i,j,k
    ! again if linoutflow tg3315, should we remove sv from this and set to true?
    if (linoutflow == .true.) then
       do k=kb,ke
          do j=jb,je
             p(ib-1,j,k)    = p(ib,j,k)                       ! inflow:  dp/dn=0
             pres0(ib-1,j,k)= pres0(ib,j,k)                   ! inflow:  dp/dn=0
             p(ie+1,j,k)    = -p(ie,j,k)                      ! outflow: p=0
             pres0(ie+1,j,k)= -pres0(ie,j,k)                  ! outflow: p=0
             up(ie+1,j,k)  = - (u0(ie+1,j,k)-u0(ie,j,k))*dxfi(ie)*uouttot
          enddo
       enddo
    else
       do k=kb,ke
          do j=jb,je
             p(ib-1,j,k)=p(ie,j,k)
             p(ie+1,j,k)=p(ib,j,k)
          enddo
       enddo
    endif

    call excj( p ,ib-1,ie+1,jb-1,je+1,kb-1,ke+1)    ! cyclic
    call excj( pres0 ,ib-1,ie+1,jb-1,je+1,kb-1,ke+1)    ! cyclic

  end subroutine bcp

  !>
  !! grwdamp damps gravity waves in the upper part of the domain.
  !>
  !! The lower limit of the damping region is set by ksp
  !! Horizontal fluctuations at the top of the domain (for instance gravity waves)
  !! are damped out by a sponge layer through an additional forcing/source term.
  !! \latexonly
  !! \begin{eqnarray}
  !! \force{i}{sp}(z) &=& -\frac{1}{t^{\mr{sp}}}\left(\xav{\fav{u_i}}-\fav{u_i}\right), \\\\
  !!  \source{\varphi}{sp}(z) &=& -\frac{1}{t^{\mr{sp}}}\left(\xav{\varphi}-\varphi\right),
  !! \end{eqnarray}
  !! with $t^{\mr{sp}}$ a relaxation time scale that goes from
  !! $t^{\mr{sp}}_0=1/(2.75\times10^{-3})\mr{s}\approx 6$min at the top of the domain
  !! to infinity at the bottom of the sponge layer.
  !! \endlatexonly
  subroutine grwdamp
    use modglobal, only : ke,kmax,cu,cv,lcoriol,igrw_damp,geodamptime
    use modfields, only : up,vp,wp,thlp,qtp,u0,v0,w0,thl0,qt0, ug,vg,thl0av,qt0av,u0av,v0av
    implicit none

    integer k

    select case(igrw_damp)
    case(0) !do nothing
    case(1)
       do k=ksp,ke
          up(:,:,k)  = up(:,:,k)-(u0(:,:,k)-(u0av(k)-cu))*tsc(k)
          vp(:,:,k)  = vp(:,:,k)-(v0(:,:,k)-(v0av(k)-cv))*tsc(k)
          wp(:,:,k)  = wp(:,:,k)-w0(:,:,k)*tsc(k)
          thlp(:,:,k)= thlp(:,:,k)-(thl0(:,:,k)-thl0av(k))*tsc(k)
          qtp(:,:,k) = qtp(:,:,k)-(qt0(:,:,k)-qt0av(k))*tsc(k)
       end do
       if(lcoriol) then
          do k=ksp,ke
             up(:,:,k)  = up(:,:,k)-(u0(:,:,k)-(ug(k)-cu))*((1./(geodamptime*rnu0))*tsc(k))
             vp(:,:,k)  = vp(:,:,k)-(v0(:,:,k)-(vg(k)-cv))*((1./(geodamptime*rnu0))*tsc(k))
          end do
       end if
    case(2)
       do k=ksp,ke
          up(:,:,k)  = up(:,:,k)-(u0(:,:,k)-(ug(k)-cu))*tsc(k)
          vp(:,:,k)  = vp(:,:,k)-(v0(:,:,k)-(vg(k)-cv))*tsc(k)
          wp(:,:,k)  = wp(:,:,k)-w0(:,:,k)*tsc(k)
          thlp(:,:,k)= thlp(:,:,k)-(thl0(:,:,k)-thl0av(k))*tsc(k)
          qtp(:,:,k) = qtp(:,:,k)-(qt0(:,:,k)-qt0av(k))*tsc(k)
       end do
    case(3)
       do k=ksp,ke
          up(:,:,k)  = up(:,:,k)-(u0(:,:,k)-(u0av(k)-cu))*tsc(k)
          vp(:,:,k)  = vp(:,:,k)-(v0(:,:,k)-(v0av(k)-cv))*tsc(k)
          wp(:,:,k)  = wp(:,:,k)-w0(:,:,k)*tsc(k)
          thlp(:,:,k)= thlp(:,:,k)-(thl0(:,:,k)-thl0av(k))*tsc(k)
          qtp(:,:,k) = qtp(:,:,k)-(qt0(:,:,k)-qt0av(k))*tsc(k)
       end do
    case default
       stop "no gravity wave damping option selected"
    end select

    return
  end subroutine grwdamp

  !> Sets top boundary conditions for scalars
  subroutine toph

    use modglobal, only : ib,ie,jb,je,kb,ke,cu,cv,numol,prandtlmol,dzh,cu,cv,&
         dtheta,dqt,dsv,nsv,lzerogradtop,lzerogradtopscal,ltempeq,lmoist,khc,ltfluxtop,numoli,lwallfunc,prandtlmol,dzf
    use modfields, only : u0,v0,e120,um,vm,e12m,thl0,qt0,sv0,thlm,qtm,svm
    use modsurfdata,only: thlflux,qtflux,svflux,ustar,thl_top,qt_top,wttop
    implicit none

    real    :: flux
    integer :: i,j,jp,jm,n,m


    !! *** No-slip wall conditions ****
    !    do j=jb,je
    !       jp=j+1
    !       jm=j-1
    !
    !       do i=ib,ie
    !        thl0(i,j,ke+1) = thl0(i,j,ke) + (thlflux(i,j)*dzh(ke)*prandtlmol )/numol
    !        qt0(i,j,ke+1) = qt0(i,j,ke) + (qtflux(i,j)*dzh(ke)*prandtlmol )/numol
    !        sv0(i,j,ke+1,:) = sv0(i,j,ke,:) + (svflux(i,j,:)*dzh(ke)*prandtlmol )/numol
    !
    !        thlm(i,j,ke+1) = thlm(i,j,ke) + (thlflux(i,j)*dzh(ke)*prandtlmol )/numol
    !        qtm(i,j,ke+1) = qtm(i,j,ke) + (qtflux(i,j)*dzh(ke)*prandtlmol )/numol
    !        svm(i,j,ke+1,:) = svm(i,j,ke,:) + (svflux(i,j,:)*dzh(ke)*prandtlmol )/numol
    !       end do
    !     end do
    !! *** End of conditions for no-slip wall ****

    ! **  Top conditions :

    if (lzerogradtopscal) then
       if (ltempeq ==.true.) then
          thl0(:,:,ke+1) = thl0(:,:,ke)
          thlm(:,:,ke+1) = thlm(:,:,ke)
       end if
       if (lmoist == .true.) then
          qt0(:,:,ke+1) = qt0(:,:,ke)
          qtm(:,:,ke+1) = qtm(:,:,ke)
       end if
       do n=1,nsv 
          sv0(:,:,ke+1,n) = sv0(:,:,ke,n)
          svm(:,:,ke+1,n) = svm(:,:,ke,n)      
       end do
    else if (ltfluxtop == .true.) then
      thl0(:,:,ke+1) = wttop*dzf(ke)+thl0(:,:,ke)
      thlm(:,:,ke+1) = wttop*dzf(ke)+thlm(:,:,ke)   
    else
       if (ltempeq ==.true.) then
          thl0(:,:,ke+1) = 2.*thl_top - thl0(:,:,ke)                       ! T=thl_top
          thlm(:,:,ke+1) = 2.*thl_top - thlm(:,:,ke)                       ! T=thl_top
       end if
       if (lmoist ==.true.) then
          qt0(:,:,ke+1) = 2.*qt_top - qt0(:,:,ke)                       ! qt=qt_top
          qtm(:,:,ke+1) = 2.*qt_top - qtm(:,:,ke)                       ! qt=qt_top
       end if
       do m = 1,khc
          do n=1,nsv
             sv0(:,:,ke+m,n) = sv0(:,:,ke,n)
             svm(:,:,ke+m,n) = svm(:,:,ke,n)
          enddo
       enddo
    endif


    return
  end subroutine toph
  !> Sets top boundary conditions for momentum
  subroutine topm

    use modglobal,      only : ib,ie,jb,je,kb,ke,ih,jh,kh,cu,cv,numol,prandtlmol,dzh,cu,cv,lwallfunc,dzf,&
         e12min,dxfi,dxf,dxhi,xh,linletgen,jgb,jge,Uinf,numoli,dzfi,lzerogradtop
    use modfields,      only : u0,v0,w0,um,vm,wm,thl0,e120,e12m,shear,wout,wouttot
    use modsubgriddata, only : ekm
    use modinlet,       only : dispthickness
    use modinletdata,   only : Uinl,ddispdxold
    use modmpi,         only : slabsumi,myid
    use modsurfdata,    only : wtsurf
    implicit none

    !real, dimension(ib:ie) :: displ
    real    :: ucu,vcv,upcu,vpcv,fu,fv,nji
    integer :: i,j,jp,jm

    if (linletgen==1 .or. linletgen==2) then
       u0(:,:,ke+1)   = u0(:,:,ke)
       v0(:,:,ke+1)   = v0(:,:,ke)
       e120(:,:,ke+1) = e12min

       um(:,:,ke+1)   = um(:,:,ke)
       vm(:,:,ke+1)   = vm(:,:,ke)
       e12m(:,:,ke+1) = e12min
       do i=ib,ie
          w0(i,:,ke+1)   = Uinf*ddispdxold
          wm(i,:,ke+1)   = Uinf*ddispdxold
       end do
       call slabsumi(wout  ,ib,ie,w0  ,ib-ih,ie+ih,jb-jh,je+jh,kb-kh,ke+kh,ib,ie,jb,je,ke+1,ke+1) ! determine vertical (j) average outflow velocity
       nji = 1./ (jge-jgb+1)
       do i=ib,ie
          wout(i) = wout(i)*dxf(i)*nji
       end do
       wouttot = sum(wout(ib:ie))/(xh(ie+1)-xh(ib))      ! Area-averaged outflow velocity
    elseif(lzerogradtop) then
       u0(:,:,ke+1)   = u0(:,:,ke)
       um(:,:,ke+1)   = um(:,:,ke)

       v0(:,:,ke+1)   = v0(:,:,ke)
       vm(:,:,ke+1)   = vm(:,:,ke)

       w0(:,:,ke+1)   = 0.0
       wm(:,:,ke+1)   = 0.0
       !e120(:,:,ke+1) = e12min           ! free slip top wall
       !e12m(:,:,ke+1) = e12min
       e120(:,:,ke+1) = e120(:,:,ke)      ! zero grad.
       e12m(:,:,ke+1) = e120(:,:,ke)      ! zero grad.
    else ! no slip top
       if (lwallfunc == .true.) then
          do j=jb,je
             do i=ib,ie
                call wallaw(Uinf-u0(i,j,ke),0.5*(thl0(i,j,ke)*dxf(i-1)+thl0(i-1,j,ke)*dxf(i))*dxhi(i),wtsurf,dzf(ke),numol,shear(i,j,ke,4))    ! compute shear stress at upper wall (also used in v. Driest damping!)   
                u0(i,j,ke+1)   = u0(i,j,ke) + shear(i,j,ke,4)*dzh(ke+1)*numoli ! compute ghost value
                um(i,j,ke+1)   = um(i,j,ke) + shear(i,j,ke,4)*dzh(ke+1)*numoli ! compute ghost value

                call wallaw(-v0(i,j,ke),0.5*(thl0(i,j,ke)+thl0(i,j-1,ke)),wtsurf,dzf(ke),numol,shear(i,j,ke,8))    ! compute shear stress at upper wall (also used in v. Driest damping!)   
                v0(i,j,ke+1)   = v0(i,j,ke) + shear(i,j,ke,8)*dzh(ke+1)*numoli 
                vm(i,j,ke+1)   = vm(i,j,ke) + shear(i,j,ke,8)*dzh(ke+1)*numoli 

                w0(i,j,ke+1)   = 0.0
                wm(i,j,ke+1)   = 0.0
             end do
          end do
       else
          shear(:,:,ke,4)= 2.*numol*(Uinf-u0(:,:,ke))*dzfi(ke)     ! x-shear at upper wall (used in wall-damping function) 
          um(:,:,ke+1)   = 2.*Uinf - um(:,:,ke)
          u0(:,:,ke+1)   = 2.*Uinf - u0(:,:,ke)

          shear(:,:,ke,8)= -2.*numol*(v0(:,:,ke))*dzfi(ke)         ! y-shear at upper wall (used in wall-damping function) 
          v0(:,:,ke+1)   = -v0(:,:,ke)
          vm(:,:,ke+1)   = -vm(:,:,ke)

          w0(:,:,ke+1)   = 0.0
          wm(:,:,ke+1)   = 0.0
       end if
       !    e120(:,:,ke+1) = e12min           ! free slip top wall
       e120(:,:,ke+1) = e120(:,:,ke)      ! zero grad.
       e12m(:,:,ke+1) = e120(:,:,ke)      ! zero grad.
    end if


    return
  end subroutine topm

  !> Sets bottom boundary conditions for bottom wall
  subroutine bottom

    use modglobal, only : ib,ie,jb,je,kb,cu,cv,numol,prandtlmol,dzh,cu,cv,lwallfunc,lwallfuncs,&
                          dxf,dxhi,dzf,dzfi,numoli,ltempeq,khc,lmoist,ltfluxbot
    use modfields, only : u0,v0,e120,um,vm,w0,wm,e12m,thl0,qt0,sv0,thlm,qtm,svm,shear
    use modsurfdata,only: thlflux,qtflux,svflux,ustar,thvs,qts,wtsurf,wqsurf
    use modsubgriddata, only: ekm
    use modsurface, only: ustar
    implicit none

    real    :: ucu,vcv,upcu,vpcv,fu,fvi,flux
    integer :: i,j,jp,jm,m

    if (lwallfunc == .true.) then
       do j=jb,je
          do i=ib,ie

             call wallaw(u0(i,j,kb),0.5*(thl0(i,j,kb)*dxf(i-1)+thl0(i-1,j,kb)*dxf(i))*dxhi(i),wtsurf,dzf(kb),numol,shear(i,j,kb,3))

             !write(*,'(3A,2I5,F6.2)') 'x', 'y', 'Bottom shear x', i, j, shear(i,j,kb,3)
             !write(*,'(3A,2I5,F6.2)') 'x', 'y', 'Bottom shear x2', i, j, shear(i,j,kb,3)

             u0(i,j,kb-1)   = u0(i,j,kb) - shear(i,j,kb,3)*dzh(kb)*numoli ! compute ghost value
             !write(*,'(A,2(1pE9.2))') 'u0(kb), u0(kb-1)', u0(i,j,kb), u0(i,j,kb-1)
             um(i,j,kb-1)   = um(i,j,kb) - shear(i,j,kb,3)*dzh(kb)*numoli ! compute ghost value

      call wallaw(v0(i,j,kb),0.5*(thl0(i,j,kb)+thl0(i,j-1,kb)),wtsurf,dzf(kb),numol,shear(i,j,kb,7))   ! compute shear stress at lower wall (later on also used in v. Driest damping!)         
             v0(i,j,kb-1)   = v0(i,j,kb) - shear(i,j,kb,7)*dzh(kb)*numoli ! compute ghost value
             vm(i,j,kb-1)   = vm(i,j,kb) - shear(i,j,kb,7)*dzh(kb)*numoli ! compute ghost value

             w0(i,j,kb)     = 0.0     
             wm(i,j,kb)     = 0.0     

             do m=1,khc
                sv0(i,j,kb-m,:) = sv0(i,j,kb,:) + (svflux(i,j,:)*dzh(kb)*prandtlmol )*numoli
                svm(i,j,kb-m,:) = svm(i,j,kb,:) + (svflux(i,j,:)*dzh(kb)*prandtlmol )*numoli
             end do
          end do
       end do
    else  ! no wall function, just no-slip
       do j=jb,je
          do i=ib,ie
             !shear(i,j,kb,3)= 2.*numol*u0(i,j,kb)*dzfi(kb)     ! x-shear at lower wall (used in wall-damping function) 
              shear(i,j,kb,3)= 2.*numol*u0(i,j,kb)*dzfi(kb)     ! x-shear at lower wall (used in wall-damping function) 
             u0(i,j,kb-1)   = -u0(i,j,kb) !- fu*dzh(kb)/numol        ! BC for u
             um(i,j,kb-1)   = -um(i,j,kb) !- fu*dzh(kb)/numol        ! BC for u

             !shear(i,j,kb,7)= 2.*numol*v0(i,j,kb)*dzfi(kb)     ! y-shear at lower wall (used in wall-damping function) 
             shear(i,j,kb,7)= 2.*numol*v0(i,j,kb)*dzfi(kb)     ! y-shear at lower wall (used in wall-damping function) 
             v0(i,j,kb-1)   = -v0(i,j,kb)! - fv*dzh(kb)/numol        ! BC for v
             vm(i,j,kb-1)   = -vm(i,j,kb)! - fv*dzh(kb)/numol        ! BC for v

             w0(i,j,kb)     = 0.0     
             wm(i,j,kb)     = 0.0    

             do m=1,khc
                sv0(i,j,kb-m,:) = sv0(i,j,kb,:) + (svflux(i,j,:)*dzh(kb)*prandtlmol )*numoli
                svm(i,j,kb-m,:) = svm(i,j,kb,:) + (svflux(i,j,:)*dzh(kb)*prandtlmol )*numoli
             end do
          end do
       end do
    end if

    ! heat
    if (ltempeq == .true.) then
     if (ltfluxbot == .true.) then
       do j=jb,je
          do i=ib,ie
             thl0(i,j,kb-1) = thl0(i,j,kb) -wtsurf*dzf(kb)
             thlm(i,j,kb-1) = thlm(i,j,kb) -wtsurf*dzf(kb)
          end do
       end do
    else
      if (lwallfuncs == .true.) then
        do j=jb,je
        do i=ib,ie
          call wallawt(thl0(i,j,kb),thvs,dzf(kb),flux)
          thl0(i,j,kb-1)   = thl0(i,j,kb) - flux*dzh(kb)*prandtlmol*numoli !compute ghost value
          thlm(i,j,kb-1)   = thlm(i,j,kb) - flux*dzh(kb)*prandtlmol*numoli !compute ghost value
        end do
        end do
      else
        do j=jb,je
        do i=ib,ie    ! again if linoutflow tg3315, should we remove sv from this and set to truee
          thl0(i,j,kb-1) = 2.*thvs - thl0(i,j,kb) ! fixed temperature                                                                                                                                                                        
!          qt0(i,j,kb-1) = qt0(i,j,kb) + (qtflux(i,j)*dzh(kb)*prandtlmol
!          )*numoli
          thlm(i,j,kb-1) = 2.*thvs - thlm(i,j,kb) ! fixed temperature
!          qtm(i,j,kb-1) = qtm(i,j,kb) + (qtflux(i,j)*dzh(kb)*prandtlmol
!          )*numoli
        end do
        end do
      end if ! lwallfunc
    end if  ! ltfluxtop
  end if


 if (lmoist == .true.) then
     if (ltfluxbot == .true.) then
       do j=jb,je
          do i=ib,ie
             qt0(i,j,kb-1) = qt0(i,j,kb) -wqsurf*dzf(kb)
             qtm(i,j,kb-1) = qtm(i,j,kb) -wqsurf*dzf(kb)
          end do
       end do
    else
      if (lwallfuncs == .true.) then
        do j=jb,je
        do i=ib,ie
          call wallawt(qt0(i,j,kb),qts,dzf(kb),flux)
          qt0(i,j,kb-1)   = qt0(i,j,kb) - flux*dzh(kb)*prandtlmol*numoli
!compute ghost value
          qtm(i,j,kb-1)   = qtm(i,j,kb) - flux*dzh(kb)*prandtlmol*numoli
!compute ghost value
        end do
        end do
      else
        do j=jb,je
        do i=ib,ie
          qt0(i,j,kb-1) = 2.*qts - qt0(i,j,kb) ! fixed temperature                                                                                                                                                                        
!          qt0(i,j,kb-1) = qt0(i,j,kb) + (qtflux(i,j)*dzh(kb)*prandtlmol
!          )*numoli
          qtm(i,j,kb-1) = 2.*qts - qtm(i,j,kb) ! fixed temperature
!          qtm(i,j,kb-1) = qtm(i,j,kb) + (qtflux(i,j)*dzh(kb)*prandtlmol
!          )*numoli
        end do
        end do
      end if ! lwallfunc
    end if  ! ltfluxtop
  end if


    e120(:,:,kb-1) = e120(:,:,kb)
    e12m(:,:,kb-1) = e12m(:,:,kb)

    return
  end subroutine bottom


  subroutine wallawt(theta,thetasurf,dx,flux)
    use modglobal,       only : lMOST,prandtlmoli,numol
    use modsurfdata,     only : Csav,horvel
    implicit none

      real, intent(in)  :: theta,thetasurf,dx
      real, intent(out) :: flux


      if (lMOST == .true.) then      ! use Monin-Obukhov similarity
!          flux = -Csav*horvel*(theta-thetasurf)
          flux = Csav*horvel*(theta-thetasurf)
      else                           ! use NO wall function
!          flux = -numol*prandtlmoli*(theta-thetasurf)/(0.5*dx)
          flux = (theta-thetasurf)*prandtlmoli*numol/(0.5*dx)
      end if
    
  end subroutine wallawt



  subroutine wallaw(utan,thadj,flux,dx,visc,tau)
  


        use modglobal,       only : fkar,grav,ltempeq,lMOST
    use modsurfdata,     only : thls,z0,z0hav,Cmav,Csav,horvel

    implicit none

    real, intent(in)  :: utan,thadj,flux,dx,visc
    real, intent(out) :: tau

    real    const1, const2, const3, const4
    real    tausub, taupow
    real    sub, dutan, utankr,utanabs
    real    aaa,bbb
    real    dxi,dx5
    real    utan2, lmo, ust

    parameter(aaa = 8.3)
    parameter(bbb = 0.1428571429)

    dxi = 1./dx
    dx5  = 0.5*dx
         if (lMOST == .true.) then      ! use Monin-Obukhov similarity
!        tau = -Cmav*horvel*utan
        tau = Cmav*horvel*utan

!! multiply with numol, because it is later divided by it
!        tau = Cmav*horvel*utan*visc  
      else      

const1 = 0.5 * (1. - bbb) * aaa ** ((1. + bbb) / (1. - bbb))
    const2 = (1. + bbb) / aaa
    const3 = aaa ** (2. / (1. - bbb))
    const4 = 2. / (1. + bbb)

    utanabs=abs(utan) 
    utankr = 0.5 * visc * dxi * const3
    dutan  = utankr - utanabs
    sub    = max (sign(1.,dutan),0.)

    tausub    = 2. * visc * utanabs * dxi
    !      taupow3   =   const1 * (visc * dxi)**(1.+bbb) + (const2 * (visc * dxi)**bbb) * utanabs
    taupow    = ( const1 * (visc * dxi)**(1.+bbb) + (const2 * (visc * dxi)**bbb) * utanabs)** const4

    !      if (taupow3<=0) then
    !        write(6,*) 'taupow3 <=0!!!'
    !      end if
    tau = sub * tausub + (1.- sub) * taupow
    tau = sign(tau,utan)  ! give tau the same sign as utan
   end if 
   return
  end subroutine wallaw

!>Set thl, qt and sv(n) equal to slab average at level kmax
!>Set thl, qt and sv(n) equal to slab average at level kmax
  subroutine tqaver

  use modmpi,    only : comm3d,mpierr,my_real, mpi_sum
  use modglobal, only : ib,ie,jb,je,ih,jh,kb,ke,nsv,rslabs
  use modfields, only : thl0,qt0,sv0
  implicit none

  real thl0a, qt0a
  real thl0al, qt0al
  integer n
  real,allocatable, dimension(:) :: sv0al, sv0a
  allocate (sv0al(nsv),sv0a(nsv))

  thl0al=sum(thl0(ib:ie,jb:je,ke))
  qt0al =sum(qt0(ib:ie,jb:je,ke))

  do n=1,nsv
    sv0al(n) = sum(sv0(ib:ie,jb:je,ke,n))
  enddo

  call MPI_ALLREDUCE(thl0al, thl0a, 1,    MY_REAL, &
                         MPI_SUM, comm3d,mpierr)
  call MPI_ALLREDUCE(qt0al, qt0a , 1,     MY_REAL, &
                         MPI_SUM, comm3d,mpierr)
  if(nsv > 0) then
    call MPI_ALLREDUCE(sv0al, sv0a , nsv,   MY_REAL, &
                           MPI_SUM, comm3d,mpierr)
  end if


  thl0a=thl0a/rslabs
  qt0a =qt0a/rslabs
  sv0a = sv0a/rslabs

  thl0(ib:ie,jb:je,ke)=thl0a
  qt0(ib:ie,jb:je,ke) =qt0a
  do n=1,nsv
    sv0(ib:ie,jb:je,ke,n) = sv0a(n)
  enddo
  deallocate (sv0al,sv0a)

  return
  end subroutine tqaver



end module modboundary
