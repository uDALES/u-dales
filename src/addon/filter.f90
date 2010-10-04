subroutine filter(v2f, ndx)
  use modglobal, only : i1,ih,i2,j1,jh,j2
  ! top hat filter
  ! for now, only filter in horizontal
  implicit none

  integer, intent(in)    :: ndx
  real,    intent(inout) :: v2f(2-ih:i1+ih,2-jh:j1+jh)
  real                   :: v2fin(2-ih:i1+ih,2-jh:j1+jh)
  integer                :: i,j,m,n
  real                   :: weight

  v2fin(:,:) = v2f(:,:)
  v2f(:,:) = 0.
  
  weight   = ndx**2.

  ! Filter width is even
  if(mod(ndx,2) == 0) then
    do j = 2,j1
      do i = 2,i1
        do m = - ndx / 2, ndx / 2
          do n = - ndx / 2, ndx / 2
            ! Check if we are on the left edge or right edge
            if(m == - ndx / 2 .or. m == ndx / 2) then
              ! Check if we are on the top edge or bottom edge
              if(n == - ndx / 2 .or. n == ndx / 2) then
                v2f(i,j) = v2f(i,j) + 0.25 * v2fin(i+m,j+n)
              else
                v2f(i,j) = v2f(i,j) + 0.5 * v2fin(i+m,j+n)
              end if
            ! We are not on left or right edge
            else
              ! Check if we are on top edge or bottom edge
              if(n == - ndx / 2 .or. n == ndx / 2) then
                v2f(i,j) = v2f(i,j) + 0.5 * v2fin(i+m,j+n)
              ! We are not on any edge
              else
                v2f(i,j) = v2f(i,j) + 1.0 * v2fin(i+m,j+n)
              end if
            end if
          end do
        end do
        v2f(i,j) = v2f(i,j) / weight
      end do
    end do
  else
    ! Filter width is odd
    do j = 2,j1
      do i = 2,i1
        do m = - (ndx - 1) / 2, (ndx - 1) / 2
          do n = - (ndx - 1) / 2, (ndx - 1) / 2
             v2f(i,j) = v2f(i,j) + v2fin(i+m,j+n)
          end do
        end do
        v2f(i,j) = v2f(i,j) / weight
      end do
    end do
  end if

  return
end subroutine filter
