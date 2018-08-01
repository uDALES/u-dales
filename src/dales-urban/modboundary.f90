!> \file modboundary.f90
!>
!>
!! \par Revision list
!! All boundary conditions are in this file, except for immersed boundaries.
!! \par Authors
!!
module modboundary


  implicit none
  save
  private
  public :: initboundary, boundary, exitboundary,grwdamp, ksp,tqaver,cyclich,&
       bcp,bcpup,masscorr,&
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
!    use modsurface, only : getobl
    implicit none
    integer i,k  

!    call getobl

    if (linoutflow) then
!       uouttot = ubulk
uouttot = cos(iangle)*ubulk     
       if (ltempeq) then
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

  end subroutine boundary
  !> Cleans up after the run
  subroutine exitboundary
    implicit none
  end subroutine exitboundary


  subroutine closurebc
    use modsubgriddata, only : ekm,ekh
    use modglobal,      only : ib,ie,jb,je,kb,ke,ih,jh,kh,numol,prandtlmoli,linoutflow,&
         iinletgen, lzerogradtop
    use modmpi,         only : excjs
    integer i,j


    ! Top and bottom
    if ((lzerogradtop) .or. (iinletgen==1) .or. (iinletgen==2)) then
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
    if (linoutflow ) then      ! inflow/outflow
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
                          ihc,jhc,khc,lSIRANEinout,dy
    use modfields, only : thl0,thlm,qt0,qtm,sv0,svm,svprof,thlprof,qtprof,uouttot,um,u0,vm,v0
    use modmpi,    only : excjs,myid,nprocs
    use modinletdata,only : ubulk
    real rk3coef
    integer k,n,m

       ! New code for nonperiodic moisture and temperature !switch here is defined in opposite way to linoutflow...
         do m = 1,ih
           thl0(ib-m,:,:)  = thl0(ie+1-m,:,:)  
           thl0(ie+m,:,:)  = thl0(ib-1+m,:,:)
           thlm(ib-m,:,:)  = thlm(ie+1-m,:,:)
           thlm(ie+m,:,:)  = thlm(ib-1+m,:,:)
         end do
       ! j-direction
       call excjs( thl0           , ib,ie,jb,je,kb,ke+kh,ih,jh)
       call excjs( thlm           , ib,ie,jb,je,kb,ke+kh,ih,jh)

           do m = 1,ih
           qt0(ib-m,:,:)   = qt0(ie+1-m,:,:)
           qt0(ie+m,:,:)   = qt0(ib-1+m,:,:)
           qtm(ib-m,:,:)   = qtm(ie+1-m,:,:)
           qtm(ie+m,:,:)   = qtm(ib-1+m,:,:)
          end do
       call excjs( qt0            , ib,ie,jb,je,kb,ke+kh,ih,jh)
       call excjs( qtm            , ib,ie,jb,je,kb,ke+kh,ih,jh)

    if (lscalinout ) then  ! convective outflow
       rk3coef = dt / (4. - dble(rk3step))
       do n=1,nsv
          do k=kb,ke+1
             sv0(ib-1,:,k,n) = 2*svprof(k,n) - sv0(ib,:,k,n)
             svm(ib-1,:,k,n) = 2*svprof(k,n) - svm(ib,:,k,n)
          end do
          sv0(ie+1,:,:,n)= sv0(ie,:,:,n) - (sv0(ie+1,:,:,n)-sv0(ie,:,:,n))*dxhi(ie+1)*rk3coef*ubulk
          svm(ie+1,:,:,n)= svm(ie,:,:,n) - (svm(ie+1,:,:,n)-svm(ie,:,:,n))*dxhi(ie+1)*rk3coef*ubulk !changed from uouttot to ubulk here !tg3315 08/11/2017
       enddo
 
    elseif (lscalrec ) then ! recycling method for scalar fields following Matheou and Bowman (2015)
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

      elseif (lSIRANEinout) then ! tg3315 boundary condition to replicate the simulation set up in SIRANE
        if (nsv>0) then
          rk3coef = dt / (4. - dble(rk3step))
          do n=1,nsv   
            do m=1,ihc
              do k=kb,ke+1
                sv0(ib-m,:,k,n) = 2*svprof(k,n) - sv0(ib-m+1,:,k,n)  !scalars have two ghost cells...???
                svm(ib-m,:,k,n) = 2*svprof(k,n) - svm(ib-m+1,:,k,n)
              end do
!              sv0(ie+m,:,:,n)= sv0(ie+m-1,:,:,n) - (sv0(ie+m,:,:,n)-sv0(ie+m-1,:,:,n))*dxhi(ie+m)*rk3coef*ubulk
!              svm(ie+m,:,:,n)= svm(ie+m-1,:,:,n) - (svm(ie+m,:,:,n)-svm(ie+m-1,:,:,n))*dxhi(ie+m)*rk3coef*ubulk !changed from uouttot to ubulk here !tg3315 08/11/2017
!              sv0(ie+m,:,:,n)= sv0(ie+m-1,:,:,n) - (sv0(ie+m,:,:,n)-sv0(ie+m-1,:,:,n))*dxhi(ie+m)*rk3coef*u0(ie+m,:,:)
!              svm(ie+m,:,:,n)= svm(ie+m-1,:,:,n) - (svm(ie+m,:,:,n)-svm(ie+m-1,:,:,n))*dxhi(ie+m)*rk3coef*um(ie+m,:,:) !changed from uouttot to ubulk here !tg3315 08/11/2017
              svm(ie+m,:,:,n)= svm(ie+m-1,:,:,n)
              sv0(ie+m,:,:,n)= sv0(ie+m-1,:,:,n)
            end do !m, ihc
          end do !nsv 
        end if !nsv>0

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

    ! tg3315 commented this out for westerly wind direction... 17.05.18
    if (lSIRANEinout) then ! do this after excjs as these overwrite changes to ghost cells in y-direction
      do m=1,jhc
        do n=1,nsv
        if (myid==0) then
          do k=kb,ke+1
            sv0(:,jb-m,k,n) = 2*svprof(k,n) - sv0(:,jb-m+1,k,n)
            svm(:,jb-m,k,n) = 2*svprof(k,n) - svm(:,jb-m+1,k,n)
          end do
        end if
        if (myid==nprocs-1) then
          !sv0(:,je+m,:,n)= sv0(:,je+m-1,:,n) - (sv0(:,je+m,:,n)-sv0(:,je+m-1,:,n))*dy*rk3coef*ubulk
          !svm(:,je+m,:,n)= svm(:,je+m-1,:,n) - (svm(:,je+m,:,n)-svm(:,je+m-1,:,n))*dy*rk3coef*ubulk !changed from uouttot to ubulk here !tg3315 08/11/2017
    !      sv0(:,je+m,:,n)= sv0(:,je+m-1,:,n) - (sv0(:,je+m,:,n)-sv0(:,je+m-1,:,n))*dy*rk3coef*v0(:,je+m,:)
    !      svm(:,je+m,:,n)= svm(:,je+m-1,:,n) - (svm(:,je+m,:,n)-svm(:,je+m-1,:,n))*dy*rk3coef*vm(:,je+m,:) !changed from uouttot to ubulk here !tg3315 08/11/2017
          svm(:,je+m,:,n)= svm(:,je+m-1,:,n)
          sv0(:,je+m,:,n)= sv0(:,je+m-1,:,n)
          end if                                                                                
        end do !m ,jhc
      end do !nsv 
    end if ! lSIRANEinout

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

    if (loneeqn) then
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

    if (lsmagorinsky) then
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
    if (ltempeq ) then
       call excjs( thl0           , ib,ie,jb,je,kb,ke+kh,ih,jh)
       call excjs( thlm           , ib,ie,jb,je,kb,ke+kh,ih,jh)
    end if
    if (lmoist ) then
       call excjs( qt0            , ib,ie,jb,je,kb,ke+kh,ih,jh)
       call excjs( qtm            , ib,ie,jb,je,kb,ke+kh,ih,jh)
    end if

    do n=1,nsv
       call excjs( sv0(:,:,:,n) , ib,ie,jb,je,kb-khc,ke+khc,ihc,jhc)
       call excjs( svm(:,:,:,n) , ib,ie,jb,je,kb-khc,ke+khc,ihc,jhc)
    enddo

    if (loneeqn) then
       call excjs( e120, ib,ie,jb,je,kb,ke+kh,ih,jh)
       call excjs( e12m, ib,ie,jb,je,kb,ke+kh,ih,jh)
       ! exchange shear components between processors
       do n=1,12 ! for all 12 components
          !!  call excjs( shear(:,:,:,n) , ib,ie,jb,je,kb,ke,0,jmax)
          call excjs( shear(:,:,:,n) , ib,ie,jb,je,kb,ke,0,1)
       end do
    end if

    if (lsmagorinsky) then  ! needed for wall-damping (lsmagorinsky and loneeqn)
       ! exchange shear components between processors
       do n=1,12 ! for all 12 components
          call excjs( shear(:,:,:,n) , ib,ie,jb,je,kb,ke,0,1)
       end do
    end if

    return
  end subroutine cyclicj

  !>set inlet and outlet boundary conditions in i-direction
  subroutine iolet

    use modglobal, only : dxhi,dxhci,xh,zh,ib,ie,jb,je,ih,jh,kb,ke,kh,nsv,rk3step,dt,iinletgen,ltempeq,lmoist,ihc
    use modfields, only : u0,um,v0,vm,w0,wm,e120,e12m,thl0,thlm,qt0,qtm,sv0,svm,uprof,vprof,e12prof,thlprof,&
         qtprof,svprof,uouttot,wouttot
    use modmpi,    only : excjs,myid
    use modinletdata, only : u0inletbcold,v0inletbcold,w0inletbcold,uminletbc,vminletbc,wminletbc,totaluold,&
         t0inletbcold,tminletbc

    real rk3coef

    integer n,i,j,k,m


    rk3coef = dt / (4. - dble(rk3step))


    ! Inlet boundary is located at ib (not ib-1)!

    ! Inlet
    if ((iinletgen == 1) .or. (iinletgen == 2)) then
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
       if (ltempeq ) then
          do k=kb,ke
             do j=jb,je
                thl0(ib-1,j,k) = t0inletbcold(j,k)
                thlm(ib-1,j,k) = tminletbc(j,k)
             end do
          end do
       end if


       if (lmoist ) then
          do k=kb,ke
             do j=jb,je
                qt0(ib-1,j,k) = 2*qtprof(k) - qt0(ib,j,k)  !watch!
                qtm(ib-1,j,k) = 2*qtprof(k) - qtm(ib,j,k)
             end do
          end do
       end if


    else  ! (if linetgen.eqv..false.)

       do j=jb-1,je+1
          do k=kb,ke+1
             ! Momentum
             u0(ib,j,k)=uprof(k)
             um(ib,j,k)=uprof(k)

             v0(ib-1,j,k)   = 2*vprof(k) - v0(ib,j,k)      ! (v(ib)+v(ib-1))/2 = vprof
             vm(ib-1,j,k)   = 2*vprof(k) - vm(ib,j,k)      ! (v(ib)+v(ib-1))/2 = vprof

             e120(ib-1,j,k) = 2*e12prof(k) - e120(ib,j,k)      ! (e12(ib)+e12(ib-1))/2=e12prof
             e12m(ib-1,j,k) = 2*e12prof(k) - e12m(ib,j,k)      ! (e12(ib)+e12(ib-1))/2=e12prof
            
             !if (lscarec) then
             !  do m = 1,ihc
             !    sv0(ib-m,j,k,1) = 2*svprof(k,1) - sv0(ib+(m-1),j,k,1)
             !    svm(ib-m,j,k,1) = 2*svprof(k,1) - svm(ib+(m-1),j,k,1)
             !  enddo
             !else
               do n=1,nsv
                 do m = 1,ihc
                   sv0(ib-m,j,k,n) = 2*svprof(k,n) - sv0(ib+(m-1),j,k,n)
                   svm(ib-m,j,k,n) = 2*svprof(k,n) - svm(ib+(m-1),j,k,n)
                 end do
               end do
             !end if

          enddo
       enddo

       ! Heat
       if (ltempeq ) then
          do j=jb-1,je+1
             do k=kb,ke+1
                thl0(ib-1,j,k) = 2*thlprof(k) - thl0(ib,j,k) 
                thlm(ib-1,j,k) = 2*thlprof(k) - thlm(ib,j,k) 
             end do
          end do
       end if

       if (lmoist ) then
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
    end if ! iinletgen==1 .or. iinletgen==2

    ! Outlet
    ! Momentum
    v0(ie+1,:,:)   = v0(ie,:,:) - (v0(ie+1,:,:)-v0(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
    w0(ie+1,:,:)   = w0(ie,:,:) - (w0(ie+1,:,:)-w0(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
    vm(ie+1,:,:)   = vm(ie,:,:) - (vm(ie+1,:,:)-vm(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
    wm(ie+1,:,:)   = wm(ie,:,:) - (wm(ie+1,:,:)-wm(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
    e120(ie+1,:,:) = e120(ie,:,:) - (e120(ie+1,:,:)-e120(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
    e12m(ie+1,:,:) = e12m(ie,:,:) - (e12m(ie+1,:,:)-e12m(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot

    ! Heat
    if (ltempeq ) then
       thl0(ie+1,:,:) = thl0(ie,:,:) - (thl0(ie+1,:,:)-thl0(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
       thlm(ie+1,:,:) = thlm(ie,:,:) - (thlm(ie+1,:,:)-thlm(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
    end if


    if (lmoist ) then
       qt0(ie+1,:,:)  = qt0(ie,:,:) - (qt0(ie+1,:,:)-qt0(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
       qtm(ie+1,:,:)  = qtm(ie,:,:) - (qtm(ie+1,:,:)-qtm(ie,:,:))*dxhi(ie+1)*rk3coef*uouttot
    end if

    ! tg3315 !changed dxhi to dxhci!?
    do n=1,nsv
       sv0(ie+1,:,:,n)= sv0(ie,:,:,n) - (sv0(ie+1,:,:,n)-sv0(ie,:,:,n))*dxhci(ie+1)*rk3coef*uouttot
       svm(ie+1,:,:,n)= svm(ie,:,:,n) - (svm(ie+1,:,:,n)-svm(ie,:,:,n))*dxhci(ie+1)*rk3coef*uouttot
    end do

    return
  end subroutine iolet

  !> correct the u-velocity to get correct mass flow rate
  subroutine masscorr

    use modglobal, only : ib,ie,jb,je,ih,jh,kb,ke,kh,jgb,jge,dzf,zh,dy,dt,rk3step,massflowrate,&
         jmax,libm,dt,rk3step,linoutflow,lmassflowr
    use modfields, only : um,up,uout,uouttot,udef
    use modmpi,    only : slabsum,myid

    real, dimension(kb:ke)             :: uoutold
    real rk3coef,rk3coefi,massflowrateold
    integer i,j,k

    if ((.not.linoutflow).and. (lmassflowr)) then
       rk3coef = dt / (4. - dble(rk3step))
       rk3coefi = 1 / rk3coef

       udef = 0.
       uout = 0.
       uoutold = 0.
       call slabsum(uout   ,kb,ke, up  ,ib-ih,ie+ih,jb-jh,je+jh,kb,ke+kh,ie,ie,jb,je,kb,ke) ! determine horizontal (j) average outflow velocity diff
       call slabsum(uoutold,kb,ke, um  ,ib-ih,ie+ih,jb-jh,je+jh,kb-kh,ke+kh,ie,ie,jb,je,kb,ke) ! determine horizontal (j) average outflow velocity old

       do k=kb,ke
          uout(k)    = rk3coef*uout(k)   *dzf(k)*dy  ! mass flow rate through each slab (density = 1000 kg/m3)
          uoutold(k) =         uoutold(k)*dzf(k)*dy  ! mass flow rate through each slab (density = 1000 kg/m3) (previous time step)
       end do
       uouttot         = sum(uout(kb:ke))                 ! mass flow rate (at outlet)
       massflowrateold = sum(uoutold(kb:ke))              ! mass flow rate (at outlet) (previous time step)
       udef =  (massflowrate - (uouttot + massflowrateold))/(((jge-jgb+1)*dy)*(zh(ke+1)-zh(kb)))   !udef=massdef/(Area*density)
       do k = kb,ke
          do j = jb,je
             do i = ib,ie
                up(i,j,k) = up(i,j,k)  + udef*rk3coefi
             end do
          end do
       end do

    end if

  end subroutine masscorr

  !>set boundary conditions pup,pvp,pwp in subroutine fillps in modpois.f90
  subroutine bcpup(pup,pvp,pwp,rk3coef)

    use modglobal, only : ib,ie,jb,je,ih,jh,kb,ke,kh,linoutflow,dxfi,iinletgen,&
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
    if (linoutflow ) then
       if ((iinletgen == 1) .or. (iinletgen == 2)) then
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
       else ! if not iinletgen
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

    use modglobal, only : ib,ie,jb,je,ih,jh,kb,ke,kh,linoutflow,dxfi
    use modfields, only : pres0,up,u0,um,uouttot
    use modmpi, only : excj

    real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh), intent(inout) :: p !< pressure
    integer i,j,k

    if (linoutflow ) then
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
    use modglobal, only : ke,kmax,lcoriol,igrw_damp,geodamptime
    use modfields, only : up,vp,wp,thlp,qtp,u0,v0,w0,thl0,qt0, ug,vg,thl0av,qt0av,u0av,v0av
    use modmpi, only : myid    
    implicit none

    integer k
    !if (myid==0) then
    !write(*,*) "up before grwdamp",up(3,3,ke)
    !end if
    select case(igrw_damp)
    case(0) !do nothing
    case(1)
       do k=ksp,ke
          up(:,:,k)  = up(:,:,k)-(u0(:,:,k)-u0av(k))*tsc(k)
          vp(:,:,k)  = vp(:,:,k)-(v0(:,:,k)-v0av(k))*tsc(k)
          wp(:,:,k)  = wp(:,:,k)-w0(:,:,k)*tsc(k)
          thlp(:,:,k)= thlp(:,:,k)-(thl0(:,:,k)-thl0av(k))*tsc(k)
          qtp(:,:,k) = qtp(:,:,k)-(qt0(:,:,k)-qt0av(k))*tsc(k)
       end do
       if(lcoriol) then
          do k=ksp,ke
             up(:,:,k)  = up(:,:,k)-(u0(:,:,k)-ug(k))*((1./(geodamptime*rnu0))*tsc(k))
             vp(:,:,k)  = vp(:,:,k)-(v0(:,:,k)-vg(k))*((1./(geodamptime*rnu0))*tsc(k))
          end do
       end if
    case(2)
       do k=ksp,ke
          up(:,:,k)  = up(:,:,k)-(u0(:,:,k)-ug(k))*tsc(k)
          vp(:,:,k)  = vp(:,:,k)-(v0(:,:,k)-vg(k))*tsc(k)
          wp(:,:,k)  = wp(:,:,k)-w0(:,:,k)*tsc(k)
          thlp(:,:,k)= thlp(:,:,k)-(thl0(:,:,k)-thl0av(k))*tsc(k)
          qtp(:,:,k) = qtp(:,:,k)-(qt0(:,:,k)-qt0av(k))*tsc(k)
       end do
    case(3)
       do k=ksp,ke
          up(:,:,k)  = up(:,:,k)-(u0(:,:,k)-u0av(k))*tsc(k)
          vp(:,:,k)  = vp(:,:,k)-(v0(:,:,k)-v0av(k))*tsc(k)
          wp(:,:,k)  = wp(:,:,k)-w0(:,:,k)*tsc(k)
          thlp(:,:,k)= thlp(:,:,k)-(thl0(:,:,k)-thl0av(k))*tsc(k)
          qtp(:,:,k) = qtp(:,:,k)-(qt0(:,:,k)-qt0av(k))*tsc(k)
       end do
    case default
       stop "no gravity wave damping option selected"
    end select
!if (myid==0) then
!write(*,*) "up after grwdamp",up(3,3,ke)
!end if

    return
  end subroutine grwdamp

  !> Sets top boundary conditions for scalars
  subroutine toph

    use modglobal, only : ib,ie,jb,je,kb,ke,numol,prandtlmol,dzh,&
         dtheta,dqt,dsv,nsv,lzerogradtop,lzerogradtopscal,ltempeq,lmoist,khc,ltfluxtop,dzf,dzhi,&
         rk3step,dt,dzfi,dzh,dzh2i
    use modfields, only : u0,v0,e120,um,vm,e12m,thl0,qt0,sv0,thlm,qtm,svm
    use modsubgriddata, only : ekh
    use modsurfdata,only: ustar,thl_top,qt_top,wttop,wqtop
    implicit none

    real    :: flux,rk3coef
    integer :: i,j,jp,jm,n,m

    ! **  Top conditions :
    !1) zero gradient = zero flux (no perpendicular mean velocity)
    !2) given flux
    !3) given temperature at boundary
    if (lzerogradtopscal) then
       if (ltempeq ) then
          thl0(:,:,ke+1) = thl0(:,:,ke)
          thlm(:,:,ke+1) = thlm(:,:,ke)
       end if
       if (lmoist) then
          qt0(:,:,ke+1) = qt0(:,:,ke)
          qtm(:,:,ke+1) = qtm(:,:,ke)
       end if
       do n=1,nsv 
          sv0(:,:,ke+1,n) = sv0(:,:,ke,n)
          svm(:,:,ke+1,n) = svm(:,:,ke,n)      
       end do
    else if (ltfluxtop) then
       if (ltempeq) then
         rk3coef = dt / (4. - dble(rk3step))
        ! tg3315 edited 18/01/18 - also requires seeting value of ekh in startup
        thl0(:,:,ke+1) = thl0(:,:,ke) + dzh(ke+1) * wttop / ( dzhi(ke+1) * (0.5*(dzf(ke)*ekh(:,:,ke+1)+dzf(ke+1)*ekh(:,:,ke)) ) )
        thlm(:,:,ke+1) = thlm(:,:,ke) + dzh(ke+1) * wttop / ( dzhi(ke+1) * (0.5*(dzf(ke)*ekh(:,:,ke+1)+dzf(ke+1)*ekh(:,:,ke)) ) )
      end if
      if (lmoist) then
        qt0(:,:,ke+1) = qt0(:,:,ke) + dzh(ke+1) * wqtop / ( dzhi(ke+1) * (0.5*(dzf(ke)*ekh(:,:,ke+1)+dzf(ke+1)*ekh(:,:,ke)) ) )
        qtm(:,:,ke+1) = qtm(:,:,ke) + dzh(ke+1) * wqtop / ( dzhi(ke+1) * (0.5*(dzf(ke)*ekh(:,:,ke+1)+dzf(ke+1)*ekh(:,:,ke)) ) )
      end if
    else
       if (ltempeq ) then
          thl0(:,:,ke+1) = 2.*thl_top - thl0(:,:,ke)                       ! T=thl_top
          thlm(:,:,ke+1) = 2.*thl_top - thlm(:,:,ke)                       ! T=thl_top
       end if
       if (lmoist) then
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

    use modglobal,      only : ib,ie,jb,je,kb,ke,ih,jh,kh,dzh,dzf,&
         e12min,dxfi,dxf,dxhi,xh,iinletgen,jgb,jge,Uinf,numoli,dzfi,lzerogradtop
    use modfields,      only : u0,v0,w0,um,vm,wm,thl0,e120,e12m,wout,wouttot
    use modsubgriddata, only : ekm
    use modinlet,       only : dispthickness
    use modinletdata,   only : Uinl,ddispdxold
    use modmpi,         only : slabsumi,myid
    implicit none
    integer :: i
    real    :: nji
    if ((iinletgen==1) .or. (iinletgen==2)) then
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
       e120(:,:,ke+1) = e12min           ! free slip top wall
       e12m(:,:,ke+1) = e12min
    end if

    return
  end subroutine topm

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
