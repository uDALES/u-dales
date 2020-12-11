!!> \file modsubgrid.f90
!!!  Calculates and applies the Sub Filter Scale diffusion
!
!>
!!  \author Jasper Tomas, TU Delft
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!! Jasper Tomas: 
!!   -Implemented Vreman model (2004)
!!   -wall-damping applied for 1-equation and Smagorinsky models
!!   -factor 2 smaller constant in 1-equation and Smagorinsky models
!!  \todo Documentation
!!

module modsubgrid
  use modsubgriddata
  implicit none
  save
  public :: subgrid, initsubgrid, exitsubgrid, subgridnamelist

contains
  subroutine initsubgrid
    use modglobal, only : ih,ib,ie,jh,jb,je,kb,ke,kh,delta,zf,fkar, &
         pi,ifnamopt,fname_options
    use modmpi, only : myid

    implicit none

    integer   :: i, k

    real :: ceps, ch
    real :: mlen

    call subgridnamelist

    allocate(ekm(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(ekh(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))
    allocate(zlt(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(sbdiss(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(sbshr(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(sbbuo(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
    allocate(csz(ib-ih:ie+ih,kb:ke+kh))
    allocate(damp(ib:ie,jb:je,kb:ke))


    damp = 1. 
    cm = cf / (2. * pi) * (1.5*alpha_kolm)**(-1.5)

    !     ch   = 2. * alpha_kolm / beta_kolm
    ch   = prandtl
    ch2  = ch-ch1

    ceps = 2. * pi / cf * (1.5*alpha_kolm)**(-1.5)
    ce1  = (cn**2)* (cm/Rigc - ch1*cm )
    ce2  = ceps - ce1

    if(cs == -1.) then
       csz(:,:)  = (cm**3/ceps)**0.25   !< Smagorinsky constant
    else
       csz(:,:)  = cs
    end if

    !    if(lmason) then
    !      do k = kb,ke+kh
    !        do i=ib-ih,ie+ih
    !          mlen   = (1. / (csz(i,k) * delta(i,k))**nmason + 1. / (fkar * zf(k))**nmason)**(-1./nmason)
    !          csz(i,k) = mlen / delta(i,k)
    !        end do
    !      end do
    !    end if

    if (myid==0) then
       write (6,*) 'cf    = ',cf
       write (6,*) 'cm    = ',cm
       write (6,*) 'ch    = ',ch
       write (6,*) 'ch1   = ',ch1
       write (6,*) 'ch2   = ',ch2
       write (6,*) 'ceps  = ',ceps
       write (6,*) 'ceps1 = ',ce1
       write (6,*) 'ceps2 = ',ce2
       write (6,*) 'cs    = ',cs
       write (6,*) 'Rigc  = ',Rigc
    endif

  end subroutine initsubgrid

  subroutine subgridnamelist
    use modglobal, only : pi,ifnamopt,fname_options,lles
    use modmpi,    only : myid, nprocs, comm3d, mpierr, my_real, mpi_logical, mpi_integer

    implicit none

    integer :: ierr

    namelist/NAMSUBGRID/ &
         ldelta,lmason, cf,cn,Rigc,Prandtl,lsmagorinsky,lvreman,loneeqn,c_vreman,cs,nmason,lbuoycorr

    if(myid==0)then
       open(ifnamopt,file=fname_options,status='old',iostat=ierr)
       read (ifnamopt,NAMSUBGRID,iostat=ierr)
       if (ierr > 0) then
          write(0, *) 'ERROR: Problem in namoptions NAMSUBGRID'
          write(0, *) 'iostat error: ', ierr
          stop 1
       endif
       write(6 ,NAMSUBGRID)
       close(ifnamopt)
    end if

    call MPI_BCAST(ldelta     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lmason     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(nmason     ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(lsmagorinsky,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lvreman    ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(loneeqn    ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(c_vreman   ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(cs         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(cf         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(cn         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(Rigc       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(Prandtl    ,1,MY_REAL   ,0,comm3d,mpierr)
    prandtli = 1./Prandtl
    write(*,*) '1/prandtl',prandtli
    if ((lsmagorinsky) .or. (lvreman) .or. (loneeqn)) then
       lles =.true.
    endif


  end subroutine subgridnamelist

  subroutine subgrid

    ! Diffusion subroutines
    ! Thijs Heus, Chiel van Heerwaarden, 15 June 2007

    use modglobal, only : ih,jh,kh,nsv, lmoist,lles, ib,ie,jb,je,kb,ke,imax,jmax,kmax,&
         ihc,jhc,khc,ltempeq
    use modfields, only : up,vp,wp,e12p,thl0,thlp,qt0,qtp,sv0,svp,shear
    use modsurfdata,only : ustar,thlflux,qtflux,svflux
    use modmpi, only : myid, comm3d, mpierr, my_real,nprocs
    implicit none
    integer n, proces

    call closure
    call diffu(up)
    call diffv(vp)
    call diffw(wp)

    if (loneeqn) call diffe(e12p)   ! only in case of LES with 1-eq subgrid model

    if (ltempeq) call diffc(ih,jh,kh,thl0,thlp)
    if (lmoist)  call diffc(ih,jh,kh,qt0,qtp)
    do n=1,nsv
       call diffc(ihc,jhc,khc,sv0(:,:,:,n),svp(:,:,:,n))
    end do
    if (loneeqn) call sources       ! only in case of LES with 1-eq subgrid model
  end subroutine

  subroutine exitsubgrid
    implicit none
    deallocate(ekm,ekh,zlt,sbdiss,sbbuo,sbshr,csz)
  end subroutine exitsubgrid

  subroutine closure

    !-----------------------------------------------------------------|
    !                                                                 |
    !*** *closure*  calculates K-coefficients                         |
    !                                                                 |
    !      Hans Cuijpers   I.M.A.U.   06/01/1995                      |
    !                                                                 |
    !     purpose.                                                    |
    !     --------                                                    |
    !                                                                 |
    !     All the K-closure factors are calculated.                   |
    !                                                                 |
    !     ekm(i,j,k) = k sub m : for velocity-closure                 |
    !     ekh(i,j,k) = k sub h : for temperture-closure               |
    !     ekh(i,j,k) = k sub h = k sub c : for concentration-closure  |
    !                                                                 |
    !     We will use the next model for these factors:               |
    !                                                                 |
    !     k sub m = 0.12 * l * sqrt(E)                                |
    !                                                                 |
    !     k sub h = k sub c = ( 1 + (2*l)/D ) * k sub m               |
    !                                                                 |
    !           where : l = mixing length  ( in model = z2 )          |
    !                   E = subgrid energy                            |
    !                   D = grid-size distance                        |
    !                                                                 |
    !**   interface.                                                  |
    !     ----------                                                  |
    !                                                                 |
    !             *closure* is called from *program*.                 |
    !                                                                 |
    !-----------------------------------------------------------------|

    use modglobal,   only : ib,ie,jb,je,kb,ke,kh,ih,jh,jmax,delta,ekmin,grav, zf, fkar,jgb,jge,&
         dxf,dxf2,dxhi,dxfi,dy2,dyi,dyiq,dzf,dzf2,dzfi,dzhi,rk3step,rslabs, &
         numol,numoli,prandtlmoli,lles, rk3step,dxfiq,dzfiq,lbuoyancy,dzh
    use modfields,   only : dthvdz,e120,u0,v0,w0,thl0,mindist,wall,shear
    use modsurfdata, only : dudz,dvdz,thvs,ustar
    use modmpi,    only : excjs, myid, nprocs, comm3d, mpierr, my_real,mpi_sum,slabsumi
    use modboundary, only : closurebc 
    use modinletdata, only : utaui
    implicit none

    real, dimension(ib:ie) :: shearbot
    real    :: strain2,mlen,uhor,distplus,utaubot,a11,a12,a13, &
         a21,a22,a23,a31,a32,a33,aa,b11,b12,b13,b21,b22, &
         b23,b33,bb,const,const2
    integer :: i,j,k,kp,km,jp,jm,im,ip,iw,jw,kw,c1,c2

    !  if (lles  .and. rk3step == 1) then        ! compute ekm and ekh only once in complete RK-cycle
    if(lsmagorinsky) then
       do k = kb,ke
          kp=k+1
          km=k-1
          do i = ib,ie
             ip=i+1
             im=i-1
             mlen        = csz(i,k) * delta(i,k)
             do j = jb,je
                jp=j+1
                jm=j-1

                iw = wall(i,j,k,1)   ! indices of closest wall
                jw = wall(i,j,k,2)-myid*jmax   ! indices of closest wall in local j-index
                kw = wall(i,j,k,3)
                c1 = wall(i,j,k,4)   ! shear stress component
                c2 = wall(i,j,k,5)   ! shear stress component
                if ((jw >= jb-1) .and. (jw <= je+1)) then      ! check if jw is within the halo of this proc
                   !write(*,'(A,E9.2,A,E9.2,A,E9.2,A,E9.2)') 'component1:', c1, 'component2:', c2, 'shear c1:', shear(iw,jw,kw,c1), 'shear c2:', shear(iw,jw,kw,c2)
                   distplus = mindist(i,j,k)*sqrt(abs(shear(iw,jw,kw,c1))+abs(shear(iw,jw,kw,c2)))*numoli
                   damp(i,j,k) = sqrt(1. - exp((-distplus*0.04)**3.))            ! Wall-damping according to Piomelli
                   !    write(*,'(A,2(1pE9.2))') 'damp, distplus', damp(i,j,k), distplus
                else
                   damp(i,j,k) = 1.
                end if


                   strain2 =  ( &
                        ((u0(ip,j,k)-u0(i,j,k))    *dxfi(i)        )**2    + &
                        ((v0(i,jp,k)-v0(i,j,k))    *dyi        )**2    + &
                        ((w0(i,j,kp)-w0(i,j,k))    *dzfi(k)     )**2    )

                   strain2 = strain2 + 0.125 * ( &
                        ((w0(i,j,kp)-w0(im,j,kp))   *dxhi(i)     + &
                        (u0(i,j,kp)-u0(i,j,k))      *dzhi(kp)  )**2    + &
                        ((w0(i,j,k)-w0(im,j,k))     *dxhi(i)     + &
                        (u0(i,j,k)-u0(i,j,km))      *dzhi(k)   )**2    + &
                        ((w0(ip,j,k)-w0(i,j,k))     *dxhi(ip)     + &
                        (u0(ip,j,k)-u0(ip,j,km))    *dzhi(k)   )**2    + &
                        ((w0(ip,j,kp)-w0(i,j,kp))   *dxhi(ip)     + &
                        (u0(ip,j,kp)-u0(ip,j,k))    *dzhi(kp)  )**2    )

                   strain2 = strain2 + 0.125 * ( &
                        ((u0(i,jp,k)-u0(i,j,k))     *dyi     + &
                        (v0(i,jp,k)-v0(im,jp,k))    *dxhi(i)        )**2    + &
                        ((u0(i,j,k)-u0(i,jm,k))     *dyi     + &
                        (v0(i,j,k)-v0(im,j,k))      *dxhi(i)        )**2    + &
                        ((u0(ip,j,k)-u0(ip,jm,k))   *dyi     + &
                        (v0(ip,j,k)-v0(i,j,k))      *dxhi(ip)       )**2    + &
                        ((u0(ip,jp,k)-u0(ip,j,k))   *dyi     + &
                        (v0(ip,jp,k)-v0(i,jp,k))    *dxhi(ip)       )**2    )

                   strain2 = strain2 + 0.125 * ( &
                        ((v0(i,j,kp)-v0(i,j,k))    *dzhi(kp) + &
                        (w0(i,j,kp)-w0(i,jm,kp))   *dyi        )**2    + &
                        ((v0(i,j,k)-v0(i,j,km))    *dzhi(k)+ &
                        (w0(i,j,k)-w0(i,jm,k))     *dyi        )**2    + &
                        ((v0(i,jp,k)-v0(i,jp,km))  *dzhi(k)+ &
                        (w0(i,jp,k)-w0(i,j,k))     *dyi        )**2    + &
                        ((v0(i,jp,kp)-v0(i,jp,k))  *dzhi(kp) + &
                        (w0(i,jp,kp)-w0(i,j,kp))   *dyi        )**2    )

                ekm(i,j,k)  = (mlen*damp(i,j,k)) ** 2. * sqrt(2. * strain2)
                ekh(i,j,k)  = ekm(i,j,k) * prandtli
             end do
          end do
       end do
       damp(:,:,:) = max(damp(:,:,:),dampmin)
       ekm(:,:,:) = ekm(:,:,:) + numol                             ! add molecular viscosity
       ekh(:,:,:) = ekh(:,:,:) + numol*prandtlmoli                 ! add molecular diffusivity  
       !write(*,'(A,3(1pE9.2))') 'strain2, ekm, ekh', strain2, ekm(10,10,10), ekh(10,10,10)     

       !    ekm(:,:,:) = max(ekm(:,:,:),ekmin)
       !    ekh(:,:,:) = max(ekh(:,:,:),ekmin)
    elseif(lvreman) then

       if ((lbuoyancy) .and. (lbuoycorr)) then
       const = prandtli*grav/(thvs*sqrt(2.*3.))
         do k = kb,ke
          kp = k+1
          km = k-1
          do j = jb,je
             jp = j+1
             jm = j-1
             do i = ib,ie        ! aij = du_i / dx_i
                ip = i+1
                im = i-1
                a11 = (u0(ip,j,k) - u0(i,j,k)) * dxfi(i)

                a12 = (((v0(ip,jp,k)+ v0(ip,j,k))*dxf(i) + &
                     (v0(i,jp,k) + v0(i,j,k)) *dxf(ip)) * dxhi(ip) &  
                     - &
                     ((v0(i,jp,k) + v0(i,j,k)) *dxf(im)  + &
                     (v0(im,jp,k)+ v0(im,j,k))*dxf(i))  * dxhi(i)) * dxfiq(i)  

                a13 = (((w0(ip,j,kp)+ w0(ip,j,k))*dxf(i) + &
                     (w0(i,j,kp) + w0(i,j,k)) *dxf(ip)) * dxhi(ip) &
                     - &
                     ((w0(i,j,kp) + w0(i,j,k)) *dxf(im)  + &
                     (w0(im,j,kp)+ w0(im,j,k))*dxf(i))  * dxhi(i)) * dxfiq(i)

                a21 = (u0(ip,jp,k) + u0(i,jp,k) - u0(ip,jm,k) - u0(i,jm,k) )*dyiq   ! simplified after writing out interpolation.

                a22 = (v0(i,jp,k) - v0(i,j,k)) * dyi 

                a23 = (w0(i,jp,kp) + w0(i,jp,k) - w0(i,jm,kp) - w0(i,jm,k) )*dyiq   ! simplified after writing out interpolation.

                a31 = (((u0(ip,j,kp) + u0(i,j,kp))*dzf(k)  + &
                     (u0(ip,j,k)  + u0(i,j,k)) *dzf(kp)) * dzhi(kp) &
                     - &
                     ((u0(ip,j,k)  + u0(i,j,k)) *dzf(km) + &
                     (u0(ip,j,km) + u0(i,j,km))*dzf(k))  * dzhi(k)   )*dzfiq(k)

                a32 = (((v0(i,jp,kp) + v0(i,j,kp))*dzf(k)  + &
                     (v0(i,jp,k)  + v0(i,j,k)) *dzf(kp)) * dzhi(kp) &
                     - &
                     ((v0(i,jp,k)  + v0(i,j,k)) *dzf(km) + &
                     (v0(i,jp,km) + v0(i,j,km))*dzf(k))  * dzhi(k)   )*dzfiq(k)

                a33 = (w0(i,j,kp) - w0(i,j,k)) * dzfi(k)

                aa  = a11*a11 + a21*a21 + a31*a31 + &
                     a12*a12 + a22*a22 + a32*a32 + &
                     a13*a13 + a23*a23 + a33*a33

                b11 = dxf2(i)*a11*a11 + dy2*a21*a21 + dzf2(k)*a31*a31
                b22 = dxf2(i)*a12*a12 + dy2*a22*a22 + dzf2(k)*a32*a32 
                b12 = dxf2(i)*a11*a12 + dy2*a21*a22 + dzf2(k)*a31*a32 
                b33 = dxf2(i)*a13*a13 + dy2*a23*a23 + dzf2(k)*a33*a33 
                b13 = dxf2(i)*a11*a13 + dy2*a21*a23 + dzf2(k)*a31*a33
                b23 = dxf2(i)*a12*a13 + dy2*a22*a23 + dzf2(k)*a32*a33
                bb = b11*b22 - b12*b12 + b11*b33 - b13*b13 + b22*b33 - b23*b23

                dthvdz(i,j,k) = (thl0(i,j,k+1)-thl0(i,j,k-1))/(dzh(k+1)+dzh(k))
                if (dthvdz(i,j,k) <= 0) then
                  const2=(bb/aa)
                else
                 ! write(*,*) "const",const
                 ! write(*,*) "delta",delta
                  const2=(bb/aa)-(delta(i,k)**4)*dthvdz(i,j,k)*const
                  if (const2 <0.0) const2 = 0.0
                end if
                ekm(i,j,k)=c_vreman*sqrt(const2)
                ekh(i,j,k)=ekm(i,j,k)*prandtli 
             end do
          end do
       end do
       ekm(:,:,:) = ekm(:,:,:) + numol                             ! add molecular viscosity
       ekh(:,:,:) = ekh(:,:,:) + numol*prandtlmoli                 ! add molecular diffusivity
       
else  ! neutral case

    do k = kb,ke
      kp = k+1
      km = k-1
      do j = jb,je
        jp = j+1
        jm = j-1
        do i = ib,ie        ! aij = du_i / dx_i
          ip = i+1
          im = i-1
          a11 = (u0(ip,j,k) - u0(i,j,k)) * dxfi(i)

          a12 = (((v0(ip,jp,k)+ v0(ip,j,k))*dxf(i) + &
                  (v0(i,jp,k) + v0(i,j,k)) *dxf(ip)) * dxhi(ip) &  
                         - &
                  ((v0(i,jp,k) + v0(i,j,k)) *dxf(im)  + &
                   (v0(im,jp,k)+ v0(im,j,k))*dxf(i))  * dxhi(i)) * dxfiq(i)  

          a13 = (((w0(ip,j,kp)+ w0(ip,j,k))*dxf(i) + &
                  (w0(i,j,kp) + w0(i,j,k)) *dxf(ip)) * dxhi(ip) &
                         - &
                 ((w0(i,j,kp) + w0(i,j,k)) *dxf(im)  + &
                  (w0(im,j,kp)+ w0(im,j,k))*dxf(i))  * dxhi(i)) * dxfiq(i)

          a21 = (u0(ip,jp,k) + u0(i,jp,k) - u0(ip,jm,k) - u0(i,jm,k) )*dyiq   !simplified after writing out interpolation.

          a22 = (v0(i,jp,k) - v0(i,j,k)) * dyi 

          a23 = (w0(i,jp,kp) + w0(i,jp,k) - w0(i,jm,kp) - w0(i,jm,k) )*dyiq   !simplified after writing out interpolation.

          a31 = (((u0(ip,j,kp) + u0(i,j,kp))*dzf(k)  + &
                  (u0(ip,j,k)  + u0(i,j,k)) *dzf(kp)) * dzhi(kp) &
                         - &
                 ((u0(ip,j,k)  + u0(i,j,k)) *dzf(km) + &
                  (u0(ip,j,km) + u0(i,j,km))*dzf(k))  * dzhi(k)   )*dzfiq(k)

          a32 = (((v0(i,jp,kp) + v0(i,j,kp))*dzf(k)  + &
                  (v0(i,jp,k)  + v0(i,j,k)) *dzf(kp)) * dzhi(kp) &
                         - &
                 ((v0(i,jp,k)  + v0(i,j,k)) *dzf(km) + &
                  (v0(i,jp,km) + v0(i,j,km))*dzf(k))  * dzhi(k)   )*dzfiq(k)

          a33 = (w0(i,j,kp) - w0(i,j,k)) * dzfi(k)

          aa  = a11*a11 + a21*a21 + a31*a31 + &
                a12*a12 + a22*a22 + a32*a32 + &
                a13*a13 + a23*a23 + a33*a33
        
          b11 = dxf2(i)*a11*a11 + dy2*a21*a21 + dzf2(k)*a31*a31
          b22 = dxf2(i)*a12*a12 + dy2*a22*a22 + dzf2(k)*a32*a32 
          b12 = dxf2(i)*a11*a12 + dy2*a21*a22 + dzf2(k)*a31*a32 
          b33 = dxf2(i)*a13*a13 + dy2*a23*a23 + dzf2(k)*a33*a33 
          b13 = dxf2(i)*a11*a13 + dy2*a21*a23 + dzf2(k)*a31*a33
          b23 = dxf2(i)*a12*a13 + dy2*a22*a23 + dzf2(k)*a32*a33
          bb = b11*b22 - b12*b12 + b11*b33 - b13*b13 + b22*b33 - b23*b23
          if (bb < 0.00000001) then
            ekm(i,j,k) = 0.
            ekh(i,j,k) = 0.
          else
            ekm(i,j,k) = c_vreman*sqrt(bb / aa) 
            ekh(i,j,k) = ekm(i,j,k)*prandtli
          end if
        end do
      end do
    end do
!    ekm(:,:,:) = max(ekm(:,:,:),ekmin)
!    ekh(:,:,:) = max(ekh(:,:,:),ekmin)
    end if ! lbuoyancy
    ekm(:,:,:) = ekm(:,:,:) + numol                             ! add molecular viscosity
    ekh(:,:,:) = ekh(:,:,:) + numol*prandtlmoli                 ! add molecular diffusivity


        !do TKE scheme
    elseif (loneeqn ) then 
       do k=kb,ke
          do j=jb,je
             do i=ib,ie
                iw = wall(i,j,k,1)   ! indices of closest wall
                jw = wall(i,j,k,2)-myid*jmax   ! indices of closest wall in local j-index
                kw = wall(i,j,k,3)
                c1 = wall(i,j,k,4)   ! shear stress component
                c2 = wall(i,j,k,5)   ! shear stress component

                !ILS13 removed near-wall damping 25.06.2014
                !if (jw >= jb-1 .and. jw <= je+1) then      ! check if jw is within the halo of this proc
                !  distplus = mindist(i,j,k)*sqrt(abs(shear(iw,jw,kw,c1))+abs(shear(iw,jw,kw,c2)))*numoli
                !  damp(i,j,k) = sqrt(1. - exp((-distplus*0.04)**3.))            ! Wall-damping according to Piomelli
                !else 
                damp(i,j,k) = 1.
                !end if
                if ((ldelta) .or. (dthvdz(i,j,k)<=0)) then
                   zlt(i,j,k) = delta(i,k)
                   ekm(i,j,k) = cm * zlt(i,j,k) *damp(i,j,k)* e120(i,j,k) !* 0.5! LES with near-wall damping !!! added factor 0.5 for shear-driven flow
                   ekh(i,j,k) = (ch1 + ch2) * ekm(i,j,k)               ! maybe ekh should be calculated from (molecular) Prandtl number
                   ekm(i,j,k) = ekm(i,j,k) + numol                     ! add molecular viscosity
                   ekh(i,j,k) = ekh(i,j,k) + numol*prandtlmoli         ! add molecular diffusivity
                else
                   !            zlt(i,j,k) = min(delta(i,k),cn*e120(i,j,k)/sqrt(grav/thvs*abs(dthvdz(i,j,k))))
                   zlt(i,j,k) = min(delta(i,k),cn*e120(i,j,k)/sqrt(grav/thvs*abs(dthvdz(i,j,k))))   !thls is used
                   ekm(i,j,k) = cm * zlt(i,j,k) *damp(i,j,k)* e120(i,j,k) !* 0.5     ! LES with near-wall damping !!! added factor 0.5 for shear-driven flow
                   ekh(i,j,k) = (ch1 + ch2 * zlt(i,j,k)/delta(i,k)) * ekm(i,j,k) !  needed in LES!
                   ekm(i,j,k) = ekm(i,j,k) + numol                     ! add molecular viscosity
                   ekh(i,j,k) = ekh(i,j,k) + numol*prandtlmoli          ! add molecular diffusivity
               endif
             end do
          end do
       end do
  
               !   write(*,'(1A,F8.4)') 'ekh', ekh(2,3,3)

      damp(:,:,:) = max(damp(:,:,:),dampmin)
       !    ekm(:,:,:) = max(ekm(:,:,:),ekmin)
       !    ekh(:,:,:) = max(ekh(:,:,:),ekmin)
    else   ! no subgrid model (DNS!)
       ekm = numol
       ekh = numol*prandtlmoli
    end if

    !  write(*,'(A,3(1pE9.2))') 'strain2, ekm, ekh', strain2, ekm(10,10,10), ekh(10,10,10)     
    !do i=ib,ie  
    !write(*,'(A,1(1pE9.2))') 'damp, ekm', damp(i,32,kb), ekm(i,32,kb)   
    !end do

    !*************************************************************
    !     Set boundary condition for K-closure factors.          ! Also other BC's!!
    !*************************************************************
    call closurebc

    return
  end subroutine closure
  subroutine sources     ! only in case of LES computation


    !-----------------------------------------------------------------|
    !                                                                 |
    !*** *sources*                                                    |
    !      calculates various terms from the subgrid TKE equation     |
    !                                                                 |
    !     Hans Cuijpers   I.M.A.U.     06/01/1995                     |
    !                                                                 |
    !     purpose.                                                    |
    !     --------                                                    |
    !                                                                 |
    !      Subroutine sources calculates all other terms in the       |
    !      subgrid energy equation, except for the diffusion terms.   |
    !      These terms are calculated in subroutine diff.             |
    !                                                                 |
    !**   interface.                                                  |
    !     ----------                                                  |
    !                                                                 |
    !     *sources* is called from *program*.                         |
    !                                                                 |
    !-----------------------------------------------------------------|

    use modglobal,   only : ib,ie,jb,je,kb,ke,delta,dxhi,dxfi,dy,dyi,dzfi,dzhi,grav,numol,prandtlmol,&
         dzh, delta
    use modfields,   only : u0,v0,w0,e120,e12p,dthvdz,thl0,thvf
    use modsurfdata,  only : dudz,dvdz,thvs
!    use modmpi,       only : myid
    implicit none

    real    tdef2,prandtlmoli
    integer i,j,k,im,ip,jm,jp,km,kp

    prandtlmoli = 1./prandtlmol 

    ! Added by J. Tomas (thermodynamics routine is bypassed)
    ! thv does not exist (equals thl), therefore dthvdz is computed here
    ! IS+HJ: thermodynamics is back in now (which calculates dthvdz)
    !    do k=kb,ke
    !      do j=jb,je
    !        do i=ib,ie
    !          dthvdz(i,j,k) = (thl0(i,j,k+1)-thl0(i,j,k-1))/(dzh(k+1)+dzh(k))
    ! !         if (dthvdz(i,j,k) == 0.0) then
    ! !           write(6,*) 'dthvdz(i,j,k)=0,i,j,k',i,j,k
    ! !         end if
    !        end do
    !      end do
    !    end do
    ! End of addition by J. Tomas (thermodynamics routine is bypassed)

    do k=kb+1,ke
       do j=jb,je
          do i=ib,ie
             kp=k+1
             km=k-1
             jp=j+1
             jm=j-1
             ip=i+1
             im=i-1

             tdef2 = 2. * ( &
                  ((u0(ip,j,k)-u0(i,j,k))      * dxfi(i)  )**2 + &
                  ((v0(i,jp,k)-v0(i,j,k))      * dyi      )**2 + &
                  ((w0(i,j,kp)-w0(i,j,k))      * dzfi(k)  )**2   )

             tdef2 = tdef2 + 0.25 * ( &
                  ((w0(i,j,kp)-w0(im,j,kp))   * dxhi(i)       + &
                  (u0(i,j,kp)-u0(i,j,k))     * dzhi(kp) )**2 + &
                  !
                  ((w0(i,j,k)-w0(im,j,k))     * dxhi(i)       + &
                  (u0(i,j,k)-u0(i,j,km))     * dzhi(k)  )**2 + &
                  !
                  ((w0(ip,j,k)-w0(i,j,k))     * dxhi(ip)      + &
                  (u0(ip,j,k)-u0(ip,j,km))   * dzhi(k)  )**2 + &
                  !
                  ((w0(ip,j,kp)-w0(i,j,kp))   * dxhi(ip)      + &
                  (u0(ip,j,kp)-u0(ip,j,k))   * dzhi(kp) )**2   )

             tdef2 = tdef2 + 0.25 * ( &
                  ((u0(i,jp,k)-u0(i,j,k))     * dyi           + &
                  (v0(i,jp,k)-v0(im,jp,k))   * dxhi(i)  )**2 + &
                  !
                  ((u0(i,j,k)-u0(i,jm,k))     * dyi           + &
                  (v0(i,j,k)-v0(im,j,k))     * dxhi(i)  )**2 + &       
                  !
                  ((u0(ip,j,k)-u0(ip,jm,k))   * dyi           + &
                  (v0(ip,j,k)-v0(i,j,k))     * dxhi(ip))**2  + &
                  !
                  ((u0(ip,jp,k)-u0(ip,j,k))   * dyi           + &
                  (v0(ip,jp,k)-v0(i,jp,k))   * dxhi(ip))**2    )

             tdef2 = tdef2 + 0.25 * ( &
                  ((v0(i,j,kp)-v0(i,j,k))     * dzhi(kp)      + &
                  (w0(i,j,kp)-w0(i,jm,kp))   * dyi      )**2 + &
                  !             
                  ((v0(i,j,k)-v0(i,j,km))     * dzhi(k)       + &
                  (w0(i,j,k)-w0(i,jm,k))     * dyi      )**2 + &
                  !
                  ((v0(i,jp,k)-v0(i,jp,km))   * dzhi(k)       + &
                  (w0(i,jp,k)-w0(i,j,k))     * dyi      )**2 + &
                  !
                  ((v0(i,jp,kp)-v0(i,jp,k))   * dzhi(kp)      + &
                  (w0(i,jp,kp)-w0(i,j,kp))   * dyi      )**2   )


             !    sbshr(i,j,k)  = ekm(i,j,k)*tdef2/ ( 2*e120(i,j,k))
             !    sbbuo(i,j,k)  = -ekh(i,j,k)*grav/thvs*dthvdz(i,j,k)/ ( 2*e120(i,j,k))
             !    sbdiss(i,j,k) = - (ce1 + ce2*zlt(i,j,k)/delta(i,k)) * e120(i,j,k)**2 /(2.*zlt(i,j,k))
             sbshr(i,j,k)  = (ekm(i,j,k)-numol)*tdef2/ ( 2*e120(i,j,k))                                   ! subtract molecular viscosity
             !    sbbuo(i,j,k)  = -(ekh(i,j,k)-numol*prandtlmoli)*grav/thvs*dthvdz(i,j,k)/ ( 2*e120(i,j,k))     ! subtract molecular diffusivity
             sbbuo(i,j,k)  = -(ekh(i,j,k)-numol*prandtlmoli)*grav/thvs*dthvdz(i,j,k)/ ( 2*e120(i,j,k))     ! subtract molecular diffusivity and use thls instead of thvs (not defined)
             !    sbdiss(i,j,k) = - (ce1 + ce2*zlt(i,j,k)/delta(i,k)) * e120(i,j,k)**2 /(2.*damp*zlt(i,j,k))   ! add near-wall damping function
             ! added factor 2. for shear-driven flow
             sbdiss(i,j,k) = - 2. * (ce1 + ce2*zlt(i,j,k)/delta(i,k)) * e120(i,j,k)**2 /(2.*damp(i,j,k)*zlt(i,j,k))   ! add near-wall damping function !! added f
          end do
       end do
    end do
    !     ----------------------------------------------end i,j,k-loop
    !    special treatment for lowest level

    do j=jb,je
       do i=ib,ie
          jp=j+1
          jm=j-1


          tdef2= 2. * ( &
               ((u0(i+1,j,kb)-u0(i,j,kb))*dxfi(i))**2 &
               + ((v0(i,jp,kb)-v0(i,j,kb))*dyi)**2 &
               + ((w0(i,j,kb+1)-w0(i,j,kb))*dzfi(kb))**2   )

          tdef2 = tdef2 + ( 0.25*(w0(i+1,j,kb+1)-w0(i-1,j,kb+1))*dxfi(i) + &
               dudz(i,j)   )**2

          tdef2 = tdef2 +   0.25 *( &
               ((u0(i,jp,kb)-u0(i,j,kb))*dyi+(v0(i,jp,kb)-v0(i-1,jp,kb))*dxfi(i))**2 &
               +((u0(i,j,kb)-u0(i,jm,kb))*dyi+(v0(i,j,kb)-v0(i-1,j,kb))*dxfi(i))**2 &
               +((u0(i+1,j,kb)-u0(i+1,jm,kb))*dyi+(v0(i+1,j,kb)-v0(i,j,kb))*dxfi(i))**2 &
               +((u0(i+1,jp,kb)-u0(i+1,j,kb))*dyi+(v0(i+1,jp,kb)-v0(i,jp,kb))*dxfi(i))**2   )

          tdef2 = tdef2 + ( 0.25*(w0(i,jp,kb+1)-w0(i,jm,kb+1))*dyi + &
               dvdz(i,j)   )**2

          ! **  Include shear and buoyancy production terms and dissipation **

          sbshr(i,j,kb)  = ekm(i,j,kb)*tdef2/ ( 2*e120(i,j,kb))
          sbbuo(i,j,kb)  = -ekh(i,j,kb)*grav/thvf(kb)*dthvdz(i,j,kb)/ ( 2*e120(i,j,kb))
          sbdiss(i,j,kb) = - (ce1 + ce2*zlt(i,j,kb)/delta(i,kb)) * e120(i,j,kb)**2 /(2.*zlt(i,j,kb))
       end do
    end do

    !    ------------------------------------------------

    e12p(ib:ie,jb:je,kb:ke) = e12p(ib:ie,jb:je,kb:ke)+ &
         sbshr(ib:ie,jb:je,kb:ke)+sbbuo(ib:ie,jb:je,kb:ke)+sbdiss(ib:ie,jb:je,kb:ke)

    !write(*,'(A,3(1pE9.2))') 'for check if called and e12p =', e12p(32,32,kb)
    return

  end subroutine sources

  subroutine diffc (hi,hj,hk,putin,putout)

    use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh,dxh,dxhi,dxh2i,dxf,dxfi,dzf,dzfi,dyi,dy2i,&
         dzhi,dzh2i,jmax, numol, prandtlmoli,lles
    use modmpi, only : myid
    implicit none

    integer, intent(in) :: hi                                                  !<size of halo in i
    integer, intent(in) :: hj                                                  !<size of halo in j
    integer, intent(in) :: hk                                                  !<size of halo in k
    real, intent(in)    :: putin (ib-hi:ie+hi,jb-hj:je+hj,kb-hk:ke+hk)
    real, intent(inout) :: putout(ib-hi:ie+hi,jb-hj:je+hj,kb   :ke+hk)

    real    cekh
    integer i,j,k,jm,jp,km,kp

    if (lles) then

       do k=kb,ke
          kp=k+1
          km=k-1

          do j=jb,je
             jp=j+1
             jm=j-1

             do i=ib,ie
                putout(i,j,k) = putout(i,j,k) &
                     +  0.5 *  ( &
                     ( (ekh(i+1,j,k)*dxf(i)+ekh(i,j,k)*dxf(i+1)) *(putin(i+1,j,k)-putin(i,j,k))*dxh2i(i+1) &
                     -(ekh(i,j,k)*dxf(i-1)+ekh(i-1,j,k)*dxf(i))*(putin(i,j,k)-putin(i-1,j,k))*dxh2i(i)) * dxfi(i) &
                     + &
                     ( (ekh(i,jp,k)+ekh(i,j,k)) *(putin(i,jp,k)-putin(i,j,k)) &
                     -(ekh(i,j,k)+ekh(i,jm,k)) *(putin(i,j,k)-putin(i,jm,k)) )*dy2i &
                     + &
                     ( (dzf(kp)*ekh(i,j,k) + dzf(k)*ekh(i,j,kp)) &
                     *  (putin(i,j,kp)-putin(i,j,k)) * dzh2i(kp) &
                     - &
                     (dzf(km)*ekh(i,j,k) + dzf(k)*ekh(i,j,km)) &
                     *  (putin(i,j,k)-putin(i,j,km)) * dzh2i(k)           )*dzfi(k) &
                     ) 
             end do
          end do
       end do

    else ! DNS
       cekh = numol* prandtlmoli
       do k=kb,ke
          kp=k+1
          km=k-1

          do j=jb,je
             jp=j+1
             jm=j-1

             do i=ib,ie
                putout(i,j,k) = putout(i,j,k) &
                     +( &
                     ( cekh *(putin(i+1,j,k)-putin(i,j,k))*dxhi(i+1) &
                     -cekh*(putin(i,j,k)-putin(i-1,j,k))*dxhi(i)) * dxfi(i) &
                     + &
                     ( cekh *(putin(i,jp,k)-putin(i,j,k)) &
                     -cekh *(putin(i,j,k)-putin(i,jm,k)) )*dy2i &
                     + &
                     ( cekh * (putin(i,j,kp)-putin(i,j,k)) * dzhi(kp) &
                     -cekh * (putin(i,j,k)-putin(i,j,km)) * dzhi(k)           )*dzfi(k) &
                     ) 
             end do
          end do
       end do

    end if ! lles=.true.

  end subroutine diffc



  subroutine diffe(putout)

    use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh,dxf,dxfi,dxh2i,dzf,dzfi,&
         dy2i,dzhi,dzh2i,jmax
    use modfields, only : e120
    use modmpi,    only : myid
    implicit none

    real, intent(inout) :: putout(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
    integer             :: i,j,k,jm,jp,km,kp

    do k=kb,ke
       kp=k+1
       km=k-1

       do j=jb,je
          jp=j+1
          jm=j-1

          do i=ib,ie

             putout(i,j,k) = putout(i,j,k) &
                  +  1.0 *  ( &
                  ( (ekm(i+1,j,k)*dxf(i)+ekm(i,j,k)*dxf(i+1)) * (e120(i+1,j,k)-e120(i,j,k))*dxh2i(i+1) &
                  -(ekm(i,j,k)*dxf(i-1)+ekm(i-1,j,k)*dxf(i))*(e120(i,j,k)-e120(i-1,j,k))*dxh2i(i)) * dxfi(i) &
                  + &
                  ((ekm(i,jp,k)+ekm(i,j,k)) *(e120(i,jp,k)-e120(i,j,k)) &
                  -(ekm(i,j,k)+ekm(i,jm,k)) *(e120(i,j,k)-e120(i,jm,k)) )*dy2i &
                  + &
                  ((dzf(kp)*ekm(i,j,k) + dzf(k)*ekm(i,j,kp)) &
                  *(e120(i,j,kp)-e120(i,j,k)) *dzh2i(kp) &
                  -(dzf(km)*ekm(i,j,k) + dzf(k)*ekm(i,j,km)) &
                  *(e120(i,j,k)-e120(i,j,km)) *dzh2i(k)        )*dzfi(k) &
                  )      
          end do
       end do
    end do


  end subroutine diffe


  subroutine diffu (putout)

    use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh,kmax,dxhi,dxf,dxfi,lles,&          
         dzf,dzfi,dy,dyi,dy2i,dzhi,dzhiq,jmax,numol
    use modfields, only : u0,v0,w0
    use modsurfdata,only : ustar
    use modmpi, only    : myid
    implicit none

    real, intent(inout) :: putout(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
    real                :: emmo,emom,emop,empo
    real                :: fu,dummy
    real                :: ucu, upcu
    integer             :: i,j,k,jm,jp,km,kp

    if (lles) then

       do k=kb,ke
          kp=k+1
          km=k-1

          do j=jb,je
             jp=j+1
             jm=j-1

             do i=ib,ie

                emom = ( dzf(km) * ( ekm(i,j,k)*dxf(i-1)  + ekm(i-1,j,k)*dxf(i) )  + &             ! dx is non-equidistant
                     dzf(k)  * ( ekm(i,j,km)*dxf(i-1) + ekm(i-1,j,km)*dxf(i) ) )*dxhi(i) * dzhiq(k) 

                emop = ( dzf(kp) * ( ekm(i,j,k)*dxf(i-1)  + ekm(i-1,j,k)*dxf(i) )  + &              ! dx is non-equidistant
                     dzf(k)  * ( ekm(i,j,kp)*dxf(i-1) + ekm(i-1,j,kp)*dxf(i) ) )*dxhi(i) * dzhiq(kp)

                empo = 0.25 * ( ( ekm(i,j,k)+ekm(i,jp,k))*dxf(i-1) + (ekm(i-1,j,k)+ekm(i-1,jp,k))*dxf(i) ) * dxhi(i)  ! dx is non-equidistant

                emmo = 0.25 * ( ( ekm(i,j,k)+ekm(i,jm,k))*dxf(i-1)  +(ekm(i-1,jm,k)+ekm(i-1,j,k))*dxf(i) ) * dxhi(i)  ! dx is non-equidistant


                ! Discretized diffusion term
                putout(i,j,k) = putout(i,j,k) &
                     + &
                     ( ekm(i,j,k)  * (u0(i+1,j,k)-u0(i,j,k))*dxfi(i) &
                     -ekm(i-1,j,k)* (u0(i,j,k)-u0(i-1,j,k))*dxfi(i-1) ) * 2. * dxhi(i) &
                     + &
                     ( empo * ( (u0(i,jp,k)-u0(i,j,k))   *dyi &
                     +(v0(i,jp,k)-v0(i-1,jp,k))*dxhi(i)) &
                     -emmo * ( (u0(i,j,k)-u0(i,jm,k))   *dyi &
                     +(v0(i,j,k)-v0(i-1,j,k))  *dxhi(i)) &  
                     ) * dyi &
                     + &
                     ( emop * ( (u0(i,j,kp)-u0(i,j,k))   *dzhi(kp) &
                     +(w0(i,j,kp)-w0(i-1,j,kp))*dxhi(i)) &
                     -emom * ( (u0(i,j,k)-u0(i,j,km))   *dzhi(k) &
                     +(w0(i,j,k)-w0(i-1,j,k))  *dxhi(i)) &
                     ) *dzfi(k) 
             end do
          end do
       end do
    else ! DNS
       do k=kb,ke
          kp=k+1
          km=k-1

          do j=jb,je
             jp=j+1
             jm=j-1

             do i=ib,ie

                ! Discretized diffusion term
                putout(i,j,k) = putout(i,j,k) &
                     + &
                     ( numol  * (u0(i+1,j,k)-u0(i,j,k))*dxfi(i) &
                     -numol * (u0(i,j,k)-u0(i-1,j,k))*dxfi(i-1) ) * 2. * dxhi(i) &
                     + &
                     ( numol * ( (u0(i,jp,k)-u0(i,j,k))   *dyi &
                     +(v0(i,jp,k)-v0(i-1,jp,k))*dxhi(i)) &
                     -numol * ( (u0(i,j,k)-u0(i,jm,k))   *dyi &
                     +(v0(i,j,k)-v0(i-1,j,k))  *dxhi(i)) &  
                     ) * dyi &
                     + &
                     ( numol * ( (u0(i,j,kp)-u0(i,j,k))   *dzhi(kp) &
                     +(w0(i,j,kp)-w0(i-1,j,kp))*dxhi(i)) &
                     -numol * ( (u0(i,j,k)-u0(i,j,km))   *dzhi(k) &
                     +(w0(i,j,k)-w0(i-1,j,k))  *dxhi(i)) &
                     ) *dzfi(k) 
             end do
          end do
       end do

    end if   ! lles

  end subroutine diffu


  subroutine diffv (putout)

    use modglobal, only   : ib,ie,ih,jb,je,jh,kb,ke,kh,dxf,dxhi,dxfi,dzf,dzfi,dyi,&
         dy2i,dzhi,dzhiq,jmax,numol,lles
    use modfields, only   : u0,v0,w0
    use modsurfdata,only  : ustar
    use modmpi, only      : myid

    implicit none

    real, intent(inout) :: putout(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
    real                :: emmo, eomm,eomp,epmo
    real                :: fv, vcv,vpcv
    integer             :: i,j,k,jm,jp,km,kp

    if (lles) then

       do k=kb,ke
          kp=k+1
          km=k-1

          do j=jb,je
             jp=j+1
             jm=j-1

             do i=ib,ie

                eomm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                     dzf(k)  * ( ekm(i,j,km) + ekm(i,jm,km) ) ) * dzhiq(k)

                eomp = ( dzf(kp) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                     dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)

                emmo = 0.25 * ( ( ekm(i,j,k)+ekm(i,jm,k))*dxf(i-1)  +(ekm(i-1,jm,k)+ekm(i-1,j,k))*dxf(i) ) * dxhi(i)  ! dx is non-equidistant

                epmo = 0.25 * ( ( ekm(i,j,k)+ekm(i,jm,k))*dxf(i+1) + (ekm(i+1,jm,k)+ekm(i+1,j,k))*dxf(i) ) * dxhi(i+1)  ! dx is non-equidistant


                ! discretized diffusion term
                putout(i,j,k) = putout(i,j,k) &
                     + &
                     ( epmo * ( (v0(i+1,j,k)-v0(i,j,k))   *dxhi(i+1) &
                     +(u0(i+1,j,k)-u0(i+1,jm,k))*dyi) &
                     -emmo * ( (v0(i,j,k)-v0(i-1,j,k))   *dxhi(i) &
                     +(u0(i,j,k)-u0(i,jm,k))    *dyi) &
                     ) * dxfi(i) &        ! = d/dx( Km*(dv/dx + du/dy) )
                     + &
                     (ekm(i,j,k) * (v0(i,jp,k)-v0(i,j,k)) &
                     -ekm(i,jm,k)* (v0(i,j,k)-v0(i,jm,k))  ) * 2. * dy2i &        ! = d/dy( 2*Km*(dv/dy) )
                     + &
                     ( eomp * ( (v0(i,j,kp)-v0(i,j,k))    *dzhi(kp) &
                     +(w0(i,j,kp)-w0(i,jm,kp))  *dyi) &
                     -eomm * ( (v0(i,j,k)-v0(i,j,km))    *dzhi(k) &
                     +(w0(i,j,k)-w0(i,jm,k))    *dyi)   &
                     ) * dzfi(k)       ! = d/dz( Km*(dv/dz + dw/dy) )
             end do
          end do
       end do

    else  ! DNS

       do k=kb,ke
          kp=k+1
          km=k-1

          do j=jb,je
             jp=j+1
             jm=j-1

             do i=ib,ie

                putout(i,j,k) = putout(i,j,k) &
                     + &
                     ( numol * ( (v0(i+1,j,k)-v0(i,j,k))   *dxhi(i+1) &
                     +(u0(i+1,j,k)-u0(i+1,jm,k))*dyi) &
                     -numol * ( (v0(i,j,k)-v0(i-1,j,k))   *dxhi(i) &
                     +(u0(i,j,k)-u0(i,jm,k))    *dyi) &
                     ) * dxfi(i) &        ! = d/dx( Km*(dv/dx + du/dy) )
                     + &
                     (numol * (v0(i,jp,k)-v0(i,j,k)) &
                     -numol * (v0(i,j,k)-v0(i,jm,k))  ) * 2. * dy2i &        ! = d/dy( 2*Km*(dv/dy) )
                     + &
                     ( numol * ( (v0(i,j,kp)-v0(i,j,k))    *dzhi(kp) &
                     +(w0(i,j,kp)-w0(i,jm,kp))  *dyi) &
                     -numol * ( (v0(i,j,k)-v0(i,j,km))    *dzhi(k) &
                     +(w0(i,j,k)-w0(i,jm,k))    *dyi)   &
                     ) * dzfi(k)       ! = d/dz( Km*(dv/dz + dw/dy) )
             end do
          end do
       end do

    end if


  end subroutine diffv



  subroutine diffw(putout)

    use modglobal, only   : ib,ie,ih,jb,je,jh,kb,ke,kh,kmax,dxf,dxhi,dxfi,dy,&
         dyi,dy2i,dzf,dzfi,dzhi,dzhiq,jmax,numol,lles
    use modfields, only   : u0,v0,w0
    use modmpi, only      : myid
    implicit none

    !*****************************************************************

    real, intent(inout) :: putout(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)
    real                :: emom, eomm, eopm, epom
    integer             :: i,j,k,jm,jp,km,kp

    if (lles) then

       do k=kb+1,ke
          kp=k+1
          km=k-1
          do j=jb,je
             jp=j+1
             jm=j-1
             do i=ib,ie

                emom = ( dzf(km) * ( ekm(i,j,k)*dxf(i-1)  + ekm(i-1,j,k)*dxf(i) )*dxhi(i)  + &
                     dzf(k)  * ( ekm(i,j,km)*dxf(i-1) + ekm(i-1,j,km)*dxf(i) )*dxhi(i) ) * dzhiq(k)

                eomm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                     dzf(k)  * ( ekm(i,j,km) + ekm(i,jm,km) ) ) * dzhiq(k)

                eopm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jp,k)  )  + &
                     dzf(k)  * ( ekm(i,j,km) + ekm(i,jp,km) ) ) * dzhiq(k)

                epom = ( dzf(km) * ( ekm(i,j,k)*dxf(i+1)  + ekm(i+1,j,k)*dxf(i) )*dxhi(i+1)  + &
                     dzf(k)  * ( ekm(i,j,km)*dxf(i+1) + ekm(i+1,j,km)*dxf(i) )*dxhi(i+1) ) * dzhiq(k)


                ! discretized diffusion term
                putout(i,j,k) = putout(i,j,k) &
                     + &
                     ( epom * ( (w0(i+1,j,k)-w0(i,j,k))    *dxhi(i+1) &
                     +(u0(i+1,j,k)-u0(i+1,j,km)) *dzhi(k) ) &
                     -emom * ( (w0(i,j,k)-w0(i-1,j,k))    *dxhi(i) &
                     +(u0(i,j,k)-u0(i,j,km))     *dzhi(k) ) &
                     )*dxfi(i) &
                     + &
                     ( eopm * ( (w0(i,jp,k)-w0(i,j,k))     *dyi &
                     +(v0(i,jp,k)-v0(i,jp,km))   *dzhi(k) ) &
                     -eomm * ( (w0(i,j,k)-w0(i,jm,k))     *dyi &
                     +(v0(i,j,k)-v0(i,j,km))     *dzhi(k) ) &
                     )*dyi &
                     + &
                     ( ekm(i,j,k) * (w0(i,j,kp)-w0(i,j,k)) *dzfi(k) &
                     -ekm(i,j,km)* (w0(i,j,k)-w0(i,j,km)) *dzfi(km) ) * 2. &
                     * dzhi(k) 
             end do
          end do
       end do

    else ! DNS

       do k=kb+1,ke
          kp=k+1
          km=k-1
          do j=jb,je
             jp=j+1
             jm=j-1
             do i=ib,ie

                ! discretized diffusion term
                putout(i,j,k) = putout(i,j,k) &
                     + &
                     ( numol * ( (w0(i+1,j,k)-w0(i,j,k))    *dxhi(i+1) &
                     +(u0(i+1,j,k)-u0(i+1,j,km)) *dzhi(k) ) &
                     -numol * ( (w0(i,j,k)-w0(i-1,j,k))    *dxhi(i) &
                     +(u0(i,j,k)-u0(i,j,km))     *dzhi(k) ) &
                     )*dxfi(i) &
                     + &
                     ( numol * ( (w0(i,jp,k)-w0(i,j,k))     *dyi &
                     +(v0(i,jp,k)-v0(i,jp,km))   *dzhi(k) ) &
                     -numol * ( (w0(i,j,k)-w0(i,jm,k))     *dyi &
                     +(v0(i,j,k)-v0(i,j,km))     *dzhi(k) ) &
                     )*dyi &
                     + &
                     ( numol * (w0(i,j,kp)-w0(i,j,k)) *dzfi(k) &
                     -numol * (w0(i,j,k)-w0(i,j,km)) *dzfi(km) ) * 2. &
                     * dzhi(k) 
             end do
          end do
       end do

    end if

  end subroutine diffw

end module modsubgrid
