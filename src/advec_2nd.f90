!> \file advec_2nd.f90
!!  Does advection with a 2nd order central differencing scheme.
!! \par Revision list
!! \par Authors
!! Second order central differencing can be used for variables where neither very
!! high accuracy nor strict monotonicity is necessary.
!! \latexonly
!!\begin{eqnarray}
!! F_{i-\frac{1}{2}}^{2nd} &=&
!!\fav{u}_{i-\frac{1}{2}}\frac{\phi_{i}+\phi_{i-1}}{2},
!!\end{eqnarray}
!! \endlatexonly
!!
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!

!> Advection at cell center
subroutine advecc_2nd(putin,putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzi5,dzf,dzh,leq
  use modfields, only : u0, v0, w0
  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the cell centered field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency

  integer :: i,j,k,ip,im,jp,jm,kp,km

  do k=1,kmax
  km=k-1
  kp=k+1
    do j=2,j1
    jm=j-1
    jp=j+1
      do i=2,i1
      im=i-1
      ip=i+1
        putout(i,j,k)  = putout(i,j,k)- (  &
              ( &
              u0(ip,j,k) * ( putin(ip,j,k) + putin(i,j,k) ) &
              -u0(i  ,j,k) * ( putin(im,j,k) + putin(i,j,k) ) &
              )* dxi5 &
            +( &
              v0(i,jp,k) * ( putin(i,jp,k) + putin(i,j,k) ) &
              -v0(i,j  ,k) * ( putin(i,jm,k) + putin(i,j,k) ) &
              )* dyi5 )

      end do
    end do
  end do

  if (leq) then

    do j=2,j1
    jm=j-1
    jp=j+1
      do i=2,i1
      im=i-1
      ip=i+1
        putout(i,j,1)  = putout(i,j,1)- ( &
                w0(i,j,2) * (putin(i,j,2) + putin(i,j,1) ) &
                ) * dzi5
        do k=2,kmax
        km=k-1
        kp=k+1
          putout(i,j,k)  = putout(i,j,k)- ( &
                w0(i,j,kp) * (putin(i,j,kp)+putin(i,j,k)) &
                -w0(i,j,k)   * (putin(i,j,km)+putin(i,j,k)) &
                )*dzi5
        end do
      end do
    end do

  else

    do j=2,j1
    jm=j-1
    jp=j+1
      do i=2,i1
      im=i-1
      ip=i+1
        putout(i,j,1)  = putout(i,j,1)- ( &
                w0(i,j,2) * (putin(i,j,2)*dzf(1) + putin(i,j,1)*dzf(2) ) / (2.*dzh(2)) &
                ) / dzf(1)

        do k=2,kmax
        km=k-1
        kp=k+1
          putout(i,j,k)  = putout(i,j,k)- ( &
                w0(i,j,kp) * (putin(i,j,kp)*dzf(k) + putin(i,j,k)*dzf(kp) ) / dzh(kp) &
               -w0(i,j,k ) * (putin(i,j,km)*dzf(k) + putin(i,j,k)*dzf(km) ) / dzh(k) &
                )/ (2. * dzf(k))
        end do
      end do
    end do

  end if

end subroutine advecc_2nd


!> Advection at the u point.
subroutine advecu_2nd(putin, putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxiq,dyiq,dziq,dzf,dzh,leq
  use modfields, only : u0, v0, w0
  implicit none


  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the u-field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency

  integer :: i,j,k,ip,im,jp,jm,kp,km

  do k=1,kmax
  km=k-1
  kp=k+1
    do j=2,j1
    jm=j-1
    jp=j+1
      do i=2,i1
      im=i-1
      ip=i+1
        putout(i,j,k)  = putout(i,j,k)- ( &
                ( &
                (putin(i,j,k)+putin(ip,j,k))*(u0(i,j,k)+u0(ip,j,k)) &
              -(putin(i,j,k)+putin(im,j,k))*(u0(i,j,k)+u0(im,j,k)) &
                )*dxiq &
                  +(  &
                (putin(i,j,k)+putin(i,jp,k))*(v0(i,jp,k)+v0(im,jp ,k)) &
              -(putin(i,j,k)+putin(i,jm,k))*(v0(i,j  ,k)+v0(im,j  ,k)) &
                )*dyiq )

      end do
    end do
  end do

  if (leq) then

    do j=2,j1
    jm=j-1
    jp=j+1
      do i=2,i1
      im=i-1
      ip=i+1
        putout(i,j,1)  = putout(i,j,1)-( &
            ( putin(i,j,2) + putin(i,j,1))*( w0(i,j,2)+ w0(im,j,2) ) &
            ) *dziq
        do k=2,kmax
        km=k-1
        kp=k+1

          putout(i,j,k)  = putout(i,j,k)- ( &
              (putin(i,j,k)+putin(i,j,kp) )*(w0(i,j,kp)+w0(im,j,kp)) &
              -(putin(i,j,k)+putin(i,j,km) )*(w0(i,j,k )+w0(im,j,k )) &
                  )*dziq
        end do
      end do
    end do

  else
    do j=2,j1
    jm=j-1
    jp=j+1
      do i=2,i1
      im=i-1
      ip=i+1
        putout(i,j,1)  = putout(i,j,1)-( &
              ( putin(i,j,2)*dzf(1) + putin(i,j,1)*dzf(2) ) / dzh(2) &
                *( w0(i,j,2)+ w0(im,j,2) ))/ (4.*dzf(1))
        do k=2,kmax
        km=k-1
        kp=k+1
          putout(i,j,k)  = putout(i,j,k)- ( &
                ( putin(i,j,kp)*dzf(k) + putin(i,j,k)*dzf(kp) ) / dzh(kp) &
                  *( w0(i,j,kp)+ w0(im,j,kp) ) &
               -( putin(i,j,k)*dzf(km) + putin(i,j,km)*dzf(k) ) / dzh(k) &
                  *( w0(i,j,k)  + w0(im,j,k)   ) &
                )/ (4.*dzf(k))
        end do
      end do
    end do




  end if

end subroutine advecu_2nd


!> Advection at the v point.
subroutine advecv_2nd(putin, putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxiq,dyiq,dziq,dzf,dzh,leq
  use modfields, only : u0, v0, w0
  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the v-field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency

  integer :: i,j,k,ip,im,jp,jm,kp,km
  do k=1,kmax
  km=k-1
  kp=k+1
    do j=2,j1
    jm=j-1
    jp=j+1
      do i=2,i1
      im=i-1
      ip=i+1

        putout(i,j,k)  = putout(i,j,k)- ( &
              ( &
              ( u0(ip,j,k)+u0(ip,jm,k))*(putin(i,j,k)+putin(ip,j,k)) &
              -(u0(i ,j,k)+u0(i ,jm,k))*(putin(i,j,k)+putin(im,j,k)) &
              )*dxiq &
              +( &
              ( v0(i,jp,k)+v0(i,j,k))*(putin(i,j,k)+putin(i,jp,k)) &
              -(v0(i,jm,k)+v0(i,j,k))*(putin(i,j,k)+putin(i,jm,k)) &
              )*dyiq )
      end do
    end do
  end do

  if (leq) then

    do j=2,j1
    jm=j-1
    jp=j+1
      do i=2,i1
      im=i-1
      ip=i+1
        putout(i,j,1)  = putout(i,j,1)- ( &
           (w0(i,j,2)+w0(i,jm,2))*(putin(i,j,2)+putin(i,j,1)) &
            )*dziq

        do k=2,kmax
        km=k-1
        kp=k+1
          putout(i,j,k)  = putout(i,j,k)- ( &
                ( w0(i,j,kp)+w0(i,jm,kp))*(putin(i,j,kp)+putin(i,j,k)) &
                -(w0(i,j,k) +w0(i,jm,k)) *(putin(i,j,km)+putin(i,j,k)) &
                )*dziq
        end do
      end do
    end do

  else
    do j=2,j1
    jm=j-1
    jp=j+1
      do i=2,i1
      im=i-1
      ip=i+1
        putout(i,j,1)  = putout(i,j,1)- ( &
           (w0(i,j,k+1)+w0(i,j-1,k+1)) &
          *(putin(i,j,k+1)*dzf(k)+putin(i,j,k)*dzf(k+1) )/ dzh(k+1) &
          -(w0(i,j,k)+w0(i,j-1,k)) &
          *(putin(i,j,k-1)*dzf(k)+putin(i,j,k)*dzf(k-1)) / dzh(k) &
           ) / (4. * dzf(k))
        do k=2,kmax
        km=k-1
        kp=k+1
          putout(i,j,k)  = putout(i,j,k)- ( &
                ( w0(i,j,kp)+w0(i,jm,kp))*(putin(i,j,kp)+putin(i,j,k)) &
                -(w0(i,j,k) +w0(i,jm,k)) *(putin(i,j,km)+putin(i,j,k)) &
                )*dziq
        end do
      end do
    end do

  end if

end subroutine advecv_2nd



!> Advection at the w point.
subroutine advecw_2nd(putin,putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxiq,dyiq,dziq,dzf,dzh,leq
  use modfields, only : u0, v0, w0
  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the w-field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency

  integer :: i,j,k,ip,im,jp,jm,kp,km

  if (leq) then


    do k=2,kmax
    km=k-1
    kp=k+1
      do j=2,j1
      jm=j-1
      jp=j+1
        do i=2,i1
        im=i-1
        ip=i+1

          putout(i,j,k)  = putout(i,j,k)- ( &
                ( &
                (putin(ip,j,k)+putin(i,j,k))*(u0(ip,j,k)+u0(ip,j,km)) &
              -(putin(im,j,k)+putin(i,j,k))*(u0(i  ,j,k)+u0(i  ,j,km)) &
                )*dxiq &
              + &
                ( &
                (putin(i,jp,k)+putin(i,j,k))*(v0(i,jp,k)+v0(i,jp,km)) &
              -(putin(i,jm,k)+putin(i,j,k))*(v0(i,j  ,k)+v0(i,j  ,km)) &
                )*dyiq &
              + &
                ( &
                (putin(i,j,k)+putin(i,j,kp) )*(w0(i,j,k) + w0(i,j,kp)) &
              -(putin(i,j,k)+putin(i,j,km) )*(w0(i,j,k) + w0(i,j,km)) &
                )*dziq &
                )

        end do
      end do
    end do
  else
    do k=2,kmax
    km=k-1
    kp=k+1
      do j=2,j1
        jm=j-1
        jp=j+1
        do i=2,i1
        im=i-1
        ip=i+1

          putout(i,j,k)  = - ( &
                ( &
                ( putin(ip,j,k) + putin(i,j,k) ) &
              *( dzf(km)*u0(ip,j,k) + dzf(k)*u0(ip,j,km) ) &
              -( putin(i,j,k) + putin(im,j,k) ) &
              *( dzf(km)*u0(i,j,k)+dzf(k)*u0(i ,j,km) ) &
                )*dxiq / dzh(k) &
              + &
                ( &
                ( putin(i,jp,k) + putin(i,j,k) ) &
              *( dzf(km)*v0(i,jp,k) + dzf(k)*v0(i,jp,km) ) &
              -( putin(i,j,k) + putin(i,j-1,k) ) &
              *( dzf(km)*v0(i,j,k) + dzf(k)*v0(i,j,km) ) &
                ) *dyiq / dzh(k) &
              + &
                ( &
                ( putin(i,j,k)+putin(i,j,kp) ) * (w0(i,j,k) + w0(i,j,kp) ) &
              -( putin(i,j,k)+putin(i,j,km) ) * (w0(i,j,k) + w0(i,j,km) ) &
                ) / (4. *dzh(k) ) &
                )

        end do
      end do
    end do

  end if

end subroutine advecw_2nd
