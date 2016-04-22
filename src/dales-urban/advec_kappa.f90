!> \file advec_kappa.f90
!!  Does advection with a kappa limiter scheme.
!! \par Revision list
!! variable x-grid now possible
!! \par Authors
!! \see Hundsdorfer et al 1995
!! Jasper Tomas June 4th 2015
!!
!! For advection of scalars that need to be strictly monotone (for example chemically reacting species) the kappa scheme has been implemented:
!! \latexonly
!! \begin{eqnarray}
!!  F_{i-\frac{1}{2}}^{\kappa} &=& \fav{u}_{i-\frac{1}{2}}
!!  \left[\phi_{i-1}+\frac{1}{2}\kappa_{i-\frac{1}{2}}\left(\phi_{i-1}-\phi_{i-2}\right)\right],
!! \end{eqnarray}
!! in case $\fav{u}>0$. $\kappa_{i-\smfrac{1}{2}}$ serves as a switch between higher order advection and first order upwind in case of strong upwind gradients of $\phi$.
!! \endlatexonly
!! This makes the scheme monotone, but also rather dissipative.
!!

!> Advection at cell center
  subroutine advecc_kappa(hi,hj,hk,var,varp)

!  use modglobal, only : i1,i2,ih,j1,j2,jh,k1,kmax,dxi,dyi,dzi
  use modglobal, only :  ib,ie,ihc,jb,je,jhc,kb,ke,khc,dxhci,dyi,dzhci,dxfc,dzfc,dxfci,dzfci,libm
  use modibmdata, only : nxwallsnorm,nywallsnorm,nzwallsnorm,xwallsnorm,& 
                         ywallsnorm,zwallsnorm,nywallsp,nywallsm,ywallsp,ywallsm
  use modfields, only :  u0, v0, w0
  implicit none
  real,external :: rlim
  integer, intent(in) :: hi                                                  !< size of halo in i
  integer, intent(in) :: hj                                                  !< size of halo in j
  integer, intent(in) :: hk                                                  !< size of halo in k
  real, dimension(ib-hi:ie+hi,jb-hj:je+hj,kb-hk:ke+hk), intent(in)  :: var !< Input: the cell centered field
  real, dimension(ib-hi:ie+hi,jb-hj:je+hj,kb   :ke+hk), intent(inout) :: varp !< Output: the tendency
!  real, dimension(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb-khc:ke+khc), intent(in)  :: var !< Input: the cell centered field
!  real, dimension(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc), intent(inout) :: varp !< Output: the tendency
  real, dimension(ib-hi:ie+hi,jb-hj:je+hj,kb   :ke+hk)      ::  duml ! 3d dummy variable: lower cell side
  real, dimension(ib-hi:ie+hi,jb-hj:je+hj,kb   :ke+hk)      ::  dumu ! 3d dummy variable: upper cell side  

  integer   i,j,k,il,iu,jl,ju,kl,ku,n
  real :: cf,d1,d2

  dumu(:,:,:) = 0.
  duml(:,:,:) = 0.
! -d(uc)/dx (stretched grid)
  do k=kb,ke
    do j=jb,je
      do i=ib,ie+1
        if (u0(i,j,k)>0) then
!          d1 = var(i-1,j,k)-var(i-2,j,k)
          d1 = (var(i-1,j,k) - var(i-2,j,k)) * dxhci(i-1)
!          d2 = var(i  ,j,k)-var(i-1,j,k)
          d2 = (var(i  ,j,k) - var(i-1,j,k)) * dxhci(i)  
          cf = var(i-1,j,k) 
        else
!          d1 = var(i  ,j,k)-var(i+1,j,k)
          d1 = (var(i  ,j,k) - var(i+1,j,k)) * dxhci(i+1)
!          d2 = var(i-1,j,k)-var(i  ,j,k)
          d2 = (var(i-1,j,k) - var(i  ,j,k)) * dxhci(i)
          cf = var(i  ,j,k)  
        end if
        cf = cf + dxfc(i)*rlim(d1,d2)
        dumu(i,j,k)   =   cf * u0(i,j,k) * dxfci(i)
        duml(i-1,j,k) = - cf * u0(i,j,k) * dxfci(i-1)
      end do
    end do
  end do

! Set advective fluxes to zero for IBM walls
  if (libm == .true.) then
    do n=1,nxwallsnorm
      i  = xwallsnorm(n,1)
      jl = xwallsnorm(n,2)
      ju = xwallsnorm(n,3)
      kl = xwallsnorm(n,4)
      ku = xwallsnorm(n,5)
      dumu(i-1,jl:ju,kl:ku) = 0. ! 'upper' flux in previous cell (i-1)
      duml(i  ,jl:ju,kl:ku) = 0. ! 'lower' flux in this cell (i)
    end do
  end if 
! Add (corrected fluxes to varp(i,j,k)
  varp(:,:,:) = varp(:,:,:) + dumu(:,:,:)+duml(:,:,:)



  dumu(:,:,:) = 0.
  duml(:,:,:) = 0.
! -d(vc)/dy (no stretched grid)
  do k=kb,ke
    do j=jb,je+1
      do i=ib,ie
        if (v0(i,j,k)>0) then
          d1 = var(i,j-1,k)-var(i,j-2,k)
          d2 = var(i,j  ,k)-var(i,j-1,k)
          cf = var(i,j-1,k)
        else
          d1 = var(i,j  ,k)-var(i,j+1,k)
          d2 = var(i,j-1,k)-var(i,j  ,k)
          cf = var(i,j  ,k)
        end if
        cf = cf + rlim(d1,d2)
        duml(i,j-1,k) = - cf * v0(i,j,k) * dyi
        dumu(i,j,k)   =   cf * v0(i,j,k) * dyi
      end do
    end do
  end do

! Set advective fluxes to zero for IBM walls
  if (libm == .true.) then
!    do n=1,nywallsnorm
!      j  = ywallsnorm(n,1)
!      il = ywallsnorm(n,2)
!      iu = ywallsnorm(n,3)
!      kl = ywallsnorm(n,4)
!      ku = ywallsnorm(n,5)
    do n=1,nywallsm
      j  = ywallsm(n,1)
      il = ywallsm(n,2)
      iu = ywallsm(n,3)-1  ! one less than shear component
      kl = ywallsm(n,4)
      ku = ywallsm(n,5)-1  ! one less than shear component
!      dumu(il:iu,j,kl:ku) = 0. ! 'upper' flux in previous cell (j-1)
      duml(il:iu,j  ,kl:ku) = 0. ! 'lower' flux in this cell (j)   
    end do
    do n=1,nywallsp
      j  = ywallsp(n,1)
      il = ywallsp(n,2)
      iu = ywallsp(n,3)-1  ! one less than shear component
      kl = ywallsp(n,4)
      ku = ywallsp(n,5)-1  ! one less than shear component
      dumu(il:iu,j,kl:ku) = 0. ! 'upper' flux in previous cell (j-1)
!      duml(il:iu,j  ,kl:ku) = 0. ! 'lower' flux in this cell (j)   
    end do
  end if
! Add (corrected fluxes to varp(i,j,k)
  varp(:,:,:) = varp(:,:,:) + dumu(:,:,:)+duml(:,:,:)



  dumu(:,:,:) = 0.
  duml(:,:,:) = 0.
! -d(wc)/dz (stretched grid)
!  do k=kb,ke+1
  do k=kb+1,ke+1
    do j=jb,je
      do i=ib,ie
        if (w0(i,j,k)>0) then
!          d1 = var(i,j,k-1)-var(i,j,k-2)
          d1 = (var(i,j,k-1) - var(i,j,k-2)) * dzhci(k-1)
!          d2 = var(i,j,k  )-var(i,j,k-1)
          d2 = (var(i,j,k  ) - var(i,j,k-1)) * dzhci(k)
          cf = var(i,j,k-1)
        else
!          d1 = var(i,j,k  )-var(i,j,k+1)
          d1 = (var(i,j,k  ) - var(i,j,k+1)) * dzhci(k+1)
!          d2 = var(i,j,k-1)-var(i,j,k  )
          d2 = (var(i,j,k-1) - var(i,j,k  )) * dzhci(k)
          cf = var(i,j,k  ) 
        end if
        cf = cf + dzfc(k)*rlim(d1,d2)
        duml(i,j,k-1) = - cf * w0(i,j,k) * dzfci(k-1)
        dumu(i,j,k)   =   cf * w0(i,j,k) * dzfci(k)
      end do
    end do
  end do


! Set advective fluxes to zero for IBM walls
  if (libm == .true.) then
    do n=1,nzwallsnorm
      k  = zwallsnorm(n,1)
      il = zwallsnorm(n,2)
      iu = zwallsnorm(n,3)
      jl = zwallsnorm(n,4)
      ju = zwallsnorm(n,5)
      if (k > kb) then
        dumu(il:iu,jl:ju,k-1) = 0.  ! 'upper' flux in previous cell (k-1)
      end if
      duml(il:iu,jl:ju,k  ) = 0.  ! 'lower' flux in this cell (k)
    end do
  end if
! Add (corrected fluxes to varp(i,j,k)
  varp(:,:,:) = varp(:,:,:) + dumu(:,:,:)+duml(:,:,:)


!  do j=2,j1
!    do i=2,i1
!      if (w0(i,j,2)>0) then
!        d1 = 0
!        d2 = var(i,j,1)-var(i,j,2)
!        cf = var(i,j,1)
!      else
!        d1 = var(i,j,2)-var(i,j,3)
!        d2 = var(i-1,j,1)-var(i  ,j,2)
!        cf = var(i  ,j,2)
!      end if
!      cf = cf + rlim(d1,d2)
!      varp(i,j,1) = varp(i,j,1) - cf * w0(i,j,2) * dzi
!      varp(i,j,2) = varp(i,j,2) + cf * w0(i,j,2) * dzi
!    end do
!  end do


  return
  end subroutine advecc_kappa


!> Determination of the limiter function
  real function rlim(d1,d2)
    use modglobal, only : eps1
    implicit none
    real, intent(in) :: d1 !< Scalar flux at 1.5 cells upwind
    real, intent(in) :: d2 !< Scalar flux at 0.5 cells upwind

    real ri,phir

    ri    = (d2+eps1)/(d1+eps1)
    phir  = max(0.,min(2.*ri,min(1./3.+2./3.*ri,2.)))
    rlim  = 0.5*phir*d1
    end function rlim



