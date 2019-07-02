!> \file scalsource.f90
!!tg3315, 8 Apr 2016 

!> Input a passive scalar into the first scalar field.

subroutine createscals

  use modglobal,  only : nsv,ib,ie,jb,je,kb,ke,ih,jh,kh,lscasrcr,cexpnr,ifinput,imax,jmax,jtot,&
                         libm,nblocks
  use modfields,  only : scar,scarl
  use modmpi,     only : myid,MY_REAL,comm3d,mpierr
  use initfac,    only : block
  implicit none
  integer :: i,j,k,n,m,p,il,iu,jl,ju

  if (lscasrcr .AND. nsv.gt.0) then

    ! read 2-D field of point sources for vehicular emissions network      
    open (ifinput, file='scals.inp.'//cexpnr)

    do j=jb,jtot
      read(ifinput, *), scar(:,j)
    end do

    scarl(ib:ie,jb:je) = scar(ib:ie,jb+myid*jmax:je+myid*jmax)

    ! Set scarl to 0 if set inside an obstacle
    if (libm) then
      do n=1,nblocks
        il = block(n,1)
        iu = block(n,2)
        jl = block(n,3)-myid*jmax
        ju = block(n,4)-myid*jmax
        if (ju < jb .or. jl > je) then
          cycle
        else
          if (ju > je) ju=je
          if (jl < jb) jl=jb
          scarl(il:iu,jl:ju) = 0.
        end if
      end do
    end if

  end if !lscasrcr

end subroutine

subroutine scalsource

  use modglobal,  only : nsv,ib,ie,jb,je,kb,ke,ih,jh,kh,ihc,jhc,khc,xf,zf,xh,zh,dy,jmax,lchem,&
                         xS,yS,zS,SS,sigS,lscasrc,lscasrcl,lscasrcr,libm,dxfi,dzfi,nblocks
  use modfields,  only : svp,svpp,scar,scarl
  use modmpi,     only : myid,mpierr,MY_REAL,comm3d,MPI_SUM
  use initfac, only : block

  implicit none
  integer :: i,j,k,n,il,iu,jl,ju,kl,ku
  real :: xL
  real :: dyi
  real :: ra2 = 0.
  real :: scalsum = 0.
  real :: scalsumt = 0.
  real :: Pi = 3.1415927

  dyi = 1./dy

  ! 2-D network of point sources at lowest level 
  if (lscasrcr .AND. nsv.gt.0) then
      if (lchem) then
        svp(ib:ie,jb:je,kb,1) = svp(ib:ie,jb:je,kb,1) + 0.9*scarl 
        svp(ib:ie,jb:je,kb,2) = svp(ib:ie,jb:je,kb,2) + 0.1*scarl
        svp(ib:ie,jb:je,kb,4) = svp(ib:ie,jb:je,kb,4) + scarl
      else
        svp(ib:ie,jb:je,kb,1) = svp(ib:ie,jb:je,kb,1) + scarl
      end if
  end if !lscasrcr

  scalsum = 0.
 
  !  Input passive scalar point sources
  if (lscasrc .AND. nsv.gt.0) then 
    do k=kb,ke
      do j=jb,je
        do i=ib,ie
               
          ra2 = (xf(i)-xS)**2 + ((j+myid*jmax-0.5)*dy-yS)**2 + (zf(k)-zS)**2
       
          if (ra2 .LE. 9*sigS**2) then
              
            scalsum = scalsum + dxfi(i) * dyi * dzfi(k) * SS*exp(-ra2/(2*sigS**2))
            svpp(i,j,k,1) = dxfi(i) * dyi * dzfi(k) * SS*exp(-ra2/(2*sigS**2))

          end if
        end do
      end do
    end do

    ! Set svpp to 0 when set inside an obstacle
    if (libm) then
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
          svpp(il:iu,jl:ju,kl:ku,:) = 0.
        end if
      end do
    end if

    ! Normalise scalar field to 1/s
    call MPI_ALLREDUCE(scalsum,scalsumt,1,MY_REAL,MPI_SUM,comm3d,mpierr)
    svpp(:,:,:,1) = svpp(:,:,:,1) / scalsumt
    svp(:,:,:,1) = svp(:,:,:,1) + svpp(:,:,:,1)

    scalsum = 0.

  end if !lscasrc

  ! Input passive scalar line sources
  if (lscasrcl .AND. nsv.gt.0) then

    if (nblocks>0) then

      do n = 1,nblocks+1
        if (n == nblocks+1) then   ! Added to run for pollutant in first canyon

          xL = xh(block(1,1) - (block(2,1) - block(1,2)+1)/2)
!          ra2 = (i - (block(1,1) - (block(2,1) - block(1,2))/2.0))**2 + (k)**2 !tg3315 commented for chem validation

        else !cycle through all other canyons

          xL = xh(block(n,2) + (block(2,1) - block(1,2)+1)/2)
!          ra2 =(xf(i) - xL)**2 + zf(k)**2

        end if


        do i=ib,ie
          do k=kb,ke

!           if (ra2 .LE. 12*sigS**2) then

              !scalsum = scalsum + dxf(i) * jmax * dy * dzf(k) * (SS/2*Pi*sigS**2) * exp(-ra2/(2*sigS**2))

              !tg3315 use this if we want to normalise th scalar conc.
!              scalsum = scalsum + ( (SS/4.) * &
!                        (erf((xh(i+1)-xL)/(sqrt(2.)*sigS)) - erf((xh(i)-xL)/(sqrt(2.)*sigS))) * &
!                        (erf((zh(k+1)-zh(kb))/(sqrt(2.)*sigS)) - erf((zh(k)-zh(kb))/(sqrt(2.)*sigS))) + &
! (SS/4.) * &
!                        (erf((xh(i+1)-xL)/(sqrt(2.)*sigS)) - erf((xh(i)-xL)/(sqrt(2.)*sigS))) * &
!                        (erf((zh(k+1)+zh(kb))/(sqrt(2.)*sigS)) - erf((zh(k)+zh(kb))/(sqrt(2.)*sigS))) ) &
!                        * dxfi(i) * dzfi(k)

              svpp(i,jb-jh:je+jh,k,1) = svpp(i,jb-jh:je+jh,k,1) + ( (SS/4.) * &
                        (erf((xh(i+1)-xL)/(sqrt(2.)*sigS)) - erf((xh(i)-xL)/(sqrt(2.)*sigS))) * &
                        (erf((zh(k+1)-zh(kb))/(sqrt(2.)*sigS)) - erf((zh(k)-zh(kb))/(sqrt(2.)*sigS))) + &
 (SS/4.) * &
                        (erf((xh(i+1)-xL)/(sqrt(2.)*sigS)) - erf((xh(i)-xL)/(sqrt(2.)*sigS))) * &
                        (erf((zh(k+1)+zh(kb))/(sqrt(2.)*sigS)) - erf((zh(k)+zh(kb))/(sqrt(2.)*sigS))) ) &
                        * dxfi(i) * dzfi(k)

!            end if
          end do
        end do
      end do
    end if !nblocks

    ! Set svpp to 0 when set inside an obstacle
    if (libm) then
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
          svpp(il:iu,jl:ju,kl:ku,:) = 0.
        end if
      end do
    end if !libm

    ! Normalise scalar field to 1/s
    call MPI_ALLREDUCE(scalsum,scalsumt,1,MY_REAL,MPI_SUM,comm3d,mpierr)

!    write(*,*), 'scalsum', scalsum

    if (lchem) then
      !svpp(:,:,:,1) = svpp(:,:,:,1)
      svp(:,:,:,1) = svp(:,:,:,1) + svpp(:,:,:,1)
      svp(:,:,:,2) = svp(:,:,:,2) + 0.1518 * svpp(:,:,:,1)
    else
!      svpp(:,:,:,1) = svpp(:,:,:,1)/ scalsumt !tg3315 not normalised 07/11/2017
      svp(:,:,:,1) = svp(:,:,:,1) + svpp(:,:,:,1)
      svp(:,:,:,2) = svp(:,:,:,2) + 0.1518 * svpp(:,:,:,1)
    end if

    svpp = 0.
    scalsum = 0.

  end if !lscasrcl

end subroutine scalsource
