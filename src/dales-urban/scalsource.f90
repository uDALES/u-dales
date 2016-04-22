!> \file scalsource.f90
!!ls13, 6 Nov 2015 

!> some nice descr
!> Note: This does not (yet) check for sources inside blocks..

subroutine scalsource

  use modglobal,  only : nsv,ib,ie,jb,je,kb,ke,ih,jh,kh,ihc,jhc,khc,xf,zf,dy,jmax,xS,yS,zS,SS,sigS,lscasrc
  use modfields,  only : svp
  use modmpi,     only : myid
  implicit none
  integer :: i,j,k
  real :: ra2 = 0.
 
 if (lscasrc==.true. .AND. nsv.gt.0) then 
  do k=kb,ke
    do j=jb,je
      do i=ib,ie
               
         ra2 = (xf(i)-xS)**2 + ((j+myid*jmax-0.5)*dy-yS)**2 + (zf(k)-zS)**2
       
           if (ra2 .LE. 9*sigS**2) then


              svp(i,j,k,1) = svp(i,j,k,1) + SS*exp(-ra2/(2*sigS**2))
           
           ! tg3315
           !else if (ra2 .GT. 9*sigS**2) then

           !  svp(i,j,k,:) = svp(i,j,k,:)
           
           ! Write in svp = 0 if within block boundries.
           
           end if
      end do
    end do
  end do
end if
end subroutine scalsource
