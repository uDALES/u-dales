!> \file advec_6th.f90
!!  Does advection with a 6th order central differencing scheme.
!! \par Revision list
!! \par Authors
!! \see Wicker and Scamarock 2002
!!
!! A higher-order accuracy in the calculation of the advection is reached with a
!! sixth order central differencing scheme.
!! \latexonly
!!!! \endlatexonly
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
subroutine advecc_6th(putin, putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi,dyi,dzf
  use modfields, only : u0, v0, w0

  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(out) :: putout

  integer :: i,j,k

!   if (leq) then

  do k=1,kmax
    do j=2,j1
      do i=2,i1

        if(k==1) then

          putout(i,j,k)  = putout(i,j,k)- ( &
                ( &
                    u0(i+1,j,k)/60. &
                    *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                    -u0(i,j,k)/60. &
                    *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                  )*dxi&
                +(&
                    v0(i,j+1,k)/60. &
                    *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                    -v0(i,j,k)/60. &
                    *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  )* dyi &
                +( &
                w0(i,j,k+1) * (putin(i,j,k+1) + putin(i,j,k)) &
                ) /(2*dzf(k)) &
                )

        elseif(k==2 .or. k==3 .or. k==kmax-1 .or. k==kmax) then
        !CvH do 2nd order for bottom and top

            putout(i,j,k)  = putout(i,j,k)- (  &
                  ( &
                      u0(i+1,j,k)/60. &
                      *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                      -u0(i,j,k)/60. &
                      *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                  )*dxi&
                +(&
                      v0(i,j+1,k)/60. &
                      *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                      -v0(i,j,k)/60. &
                      *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  )* dyi &
                +( &
                  w0(i,j,k+1) * (putin(i,j,k+1)+putin(i,j,k)) &
                  -w0(i,j,k)   * (putin(i,j,k-1)+putin(i,j,k)) &
                  )/(2*dzf(k)) &
                  )

        else

            putout(i,j,k)  = putout(i,j,k)- (  &
                  ( &
                      u0(i+1,j,k)/60. &
                      *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                      -u0(i,j,k)/60. &
                      *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                  )* dxi &
                +( &
                      v0(i,j+1,k)/60. &
                      *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                      -v0(i,j,k)/60. &
                      *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  )*dyi&
                +(&
                      w0(i,j,k+1)/60. &
                      *(37.*(putin(i,j,k+1)+putin(i,j,k))-8.*(putin(i,j,k+2)+putin(i,j,k-1))+(putin(i,j,k+3)+putin(i,j,k-2)))&
                      -w0(i,j,k)/60. &
                      *(37.*(putin(i,j,k)+putin(i,j,k-1))-8.*(putin(i,j,k+1)+putin(i,j,k-2))+(putin(i,j,k+2)+putin(i,j,k-3)))&
                  )/(dzf(k)) &
                  )
        end if

      end do
    end do
  end do

!   end if

end subroutine advecc_6th


!> Advection at the u point.
subroutine advecu_6th(putin,putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzf
  use modfields, only : u0, v0, w0

  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the u field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency

  integer :: i,j,k

!   if (leq) then

    do k=1,kmax
      do j=2,j1
        do i=2,i1

        if(k==1) then

          putout(i,j,k)  = putout(i,j,k)- ( &
                ( &
                    (u0(i+1,j,k)+u0(i,j,k))/60. &
                    *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i-1,j,k))/60. &
                    *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                )*dxi5&
              +(&
                    (v0(i,j+1,k)+v0(i-1,j+1,k))/60. &
                    *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i-1,j,k))/60. &
                    *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                )* dyi5 &
              +( &
                ( putin(i,j,k+1) + putin(i,j,k)) *(w0(i,j,k+1)+ w0(i-1,j,k+1)) &
                )/(4*dzf(k)) &
                )

        elseif(k==2 .or. k==3 .or. k==kmax-1 .or. k==kmax) then

          putout(i,j,k)  = putout(i,j,k)- ( &
                (&
                    (u0(i+1,j,k)+u0(i,j,k))/60. &
                    *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i-1,j,k))/60. &
                    *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                )*dxi5&
              +(&
                    (v0(i,j+1,k)+v0(i-1,j+1,k))/60. &
                    *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i-1,j,k))/60. &
                    *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                )* dyi5 &
              +( &
                (putin(i,j,k)+putin(i,j,k+1) )*(w0(i,j,k+1)+w0(i-1,j,k+1)) &
              -(putin(i,j,k)+putin(i,j,k-1) )*(w0(i,j,k  )+w0(i-1,j,k  )) &
                )/(4*dzf(k)) &
                )

        else

          putout(i,j,k)  = putout(i,j,k)- ( &
                  ( &
                      (u0(i+1,j,k)+u0(i,j,k))/60. &
                      *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                      -(u0(i,j,k)+u0(i-1,j,k))/60. &
                      *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                  )*dxi5&
                +(&
                      (v0(i,j+1,k)+v0(i-1,j+1,k))/60. &
                      *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                      -(v0(i,j,k)+v0(i-1,j,k))/60. &
                      *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  )*dyi5&
                + (&
                      (w0(i,j,k+1)+w0(i-1,j,k+1))/60. &
                      *(37.*(putin(i,j,k+1)+putin(i,j,k))-8.*(putin(i,j,k+2)+putin(i,j,k-1))+(putin(i,j,k+3)+putin(i,j,k-2)))&
                      -(w0(i,j,k)+w0(i-1,j,k))/60. &
                      *(37.*(putin(i,j,k)+putin(i,j,k-1))-8.*(putin(i,j,k+1)+putin(i,j,k-2))+(putin(i,j,k+2)+putin(i,j,k-3)))&
                  )/(2*dzf(k)) &
                  )

        end if

        end do
      end do
    end do

!   end if

end subroutine advecu_6th



!> Advection at the v point.
subroutine advecv_6th(putin, putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzf
  use modfields, only : u0, v0, w0
  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the v field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency

  integer :: i,j,k

!   if (leq) then

    do k=1,kmax
      do j=2,j1
        do i=2,i1

        if(k==1) then

          putout(i,j,k)  = putout(i,j,k)- ( &
                ( &
                    (u0(i+1,j,k)+u0(i+1,j-1,k))/60. &
                    *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i,j-1,k))/60. &
                    *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                 )*dxi5&
                +(&
                    (v0(i,j+1,k)+v0(i,j,k))/60. &
                    *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i,j-1,k))/60. &
                    *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  )* dyi5 &
                +( &
                  (w0(i,j,k+1)+w0(i,j-1,k+1)) *(putin(i,j,k+1)+putin(i,j,k)) &
                  ) /(4*dzf(k)) &
                  )

        elseif(k==2 .or. k==3 .or. k==kmax-1 .or. k==kmax) then

          putout(i,j,k)  = putout(i,j,k)- ( &
                ( &
                    (u0(i+1,j,k)+u0(i+1,j-1,k))/60. &
                    *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                    -(u0(i,j,k)+u0(i,j-1,k))/60. &
                    *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                 )*dxi5&
                +(&
                    (v0(i,j+1,k)+v0(i,j,k))/60. &
                    *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                    -(v0(i,j,k)+v0(i,j-1,k))/60. &
                    *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  )* dyi5 &
                +( &
                  (w0(i,j,k+1)+w0(i,j-1,k+1))*(putin(i,j,k+1)+putin(i,j,k)) &
                -(w0(i,j,k)  +w0(i,j-1,k))  *(putin(i,j,k-1)+putin(i,j,k)) &
                  )/(4*dzf(k)) &
                  )

        else

          putout(i,j,k)  = putout(i,j,k)- ( &
                  ( &
                      (u0(i+1,j,k)+u0(i+1,j-1,k))/60. &
                      *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                      -(u0(i,j,k)+u0(i,j-1,k))/60. &
                      *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                  )*dxi5&
                +(&
                      (v0(i,j+1,k)+v0(i,j,k))/60. &
                      *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                      -(v0(i,j,k)+v0(i,j-1,k))/60. &
                      *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  )*dyi5&
                +(&
                      (w0(i,j,k+1)+w0(i,j-1,k+1))/60. &
                      *(37.*(putin(i,j,k+1)+putin(i,j,k))-8.*(putin(i,j,k+2)+putin(i,j,k-1))+(putin(i,j,k+3)+putin(i,j,k-2)))&
                      -(w0(i,j,k)+w0(i,j-1,k))/60. &
                      *(37.*(putin(i,j,k)+putin(i,j,k-1))-8.*(putin(i,j,k+1)+putin(i,j,k-2))+(putin(i,j,k+2)+putin(i,j,k-3)))&
                  )/(2*dzf(k)) &
                  )

        end if

        end do
      end do
    end do

!   end if

end subroutine advecv_6th



!> Advection at the w point.
subroutine advecw_6th(putin, putout)

  use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi5,dyi5,dzh
  use modfields, only : u0, v0, w0
  implicit none

  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: putin !< Input: the w field
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout !< Output: the tendency

  integer :: i,j,k

!   if (leq) then

    do k=2,kmax
      do j=2,j1
        do i=2,i1

          if(k==2 .or. k==3 .or. k==kmax-1 .or. k==kmax) then
            putout(i,j,k)  = putout(i,j,k)- ( &
                 (&
                     (u0(i+1,j,k)+u0(i+1,j,k-1))/60. &
                     *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                     -(u0(i,j,k)+u0(i,j,k-1))/60. &
                     *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                 )*dxi5&
                +(&
                     (v0(i,j+1,k)+v0(i,j+1,k-1))/60. &
                     *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                     -(v0(i,j,k)+v0(i,j,k-1))/60. &
                     *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  )* dyi5 &
                + ( &
                  (putin(i,j,k)+putin(i,j,k+1) )*(w0(i,j,k) + w0(i,j,k+1)) &
                -(putin(i,j,k)+putin(i,j,k-1) )*(w0(i,j,k) + w0(i,j,k-1)) &
                  )/(4*dzh(k)) &
                  )
          else

            putout(i,j,k)  = putout(i,j,k)- ( &
                  (&
                      (u0(i+1,j,k)+u0(i+1,j,k-1))/60. &
                      *(37.*(putin(i+1,j,k)+putin(i,j,k))-8.*(putin(i+2,j,k)+putin(i-1,j,k))+(putin(i+3,j,k)+putin(i-2,j,k)))&
                      -(u0(i,j,k)+u0(i,j,k-1))/60. &
                      *(37.*(putin(i,j,k)+putin(i-1,j,k))-8.*(putin(i+1,j,k)+putin(i-2,j,k))+(putin(i+2,j,k)+putin(i-3,j,k)))&
                  )*dxi5&
                +(&
                      (v0(i,j+1,k)+v0(i,j+1,k-1))/60. &
                      *(37.*(putin(i,j+1,k)+putin(i,j,k))-8.*(putin(i,j+2,k)+putin(i,j-1,k))+(putin(i,j+3,k)+putin(i,j-2,k)))&
                      -(v0(i,j,k)+v0(i,j,k-1))/60. &
                      *(37.*(putin(i,j,k)+putin(i,j-1,k))-8.*(putin(i,j+1,k)+putin(i,j-2,k))+(putin(i,j+2,k)+putin(i,j-3,k)))&
                  )*dyi5&
                + (&
                      (w0(i,j,k)+w0(i,j,k+1))/60. &
                      *(37.*(putin(i,j,k+1)+putin(i,j,k))-8.*(putin(i,j,k+2)+putin(i,j,k-1))+(putin(i,j,k+3)+putin(i,j,k-2)))&
                      -(w0(i,j,k)+w0(i,j,k-1))/60. &
                      *(37.*(putin(i,j,k)+putin(i,j,k-1))-8.*(putin(i,j,k+1)+putin(i,j,k-2))+(putin(i,j,k+2)+putin(i,j,k-3)))&
                  )/(2*dzh(k)) &
                  )
          end if
        end do
      end do
     end do

!   end if

end subroutine advecw_6th