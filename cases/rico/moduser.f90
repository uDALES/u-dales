!> \file moduser.f90
!! A dummy file for cases where one wants additional forcings
!----------------------------------------------------------------------------
! This file is part of DALES.
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
! Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!----------------------------------------------------------------------------
!
!
module moduser
implicit none
contains
subroutine force_user
  implicit none
  end subroutine force_user

subroutine rad_user
  implicit none
end subroutine rad_user

subroutine micro_user
  implicit none
end subroutine micro_user

subroutine initsurf_user
    use modglobal,   only : tmelt,bt,at,rd,rv,cp,es0,pref0
    use modsurfdata, only : ps,qts,thvs,thls
    real       :: exner, tsurf, es

    exner      = (ps / pref0)**(rd/cp)
    tsurf      = thls * exner
    es         = es0 * exp(at*(tsurf-tmelt) / (tsurf-bt))
    qts        = rd / rv * es / ps
    thvs = thls * (1. + (rv/rd - 1.) * qts)
end subroutine initsurf_user

subroutine surf_user
 use modglobal,  only : zf,i1,j1,i2,j2,grav,nsv,fkar,cv,cu
 use modsurfdata,only : ustar,dudz,dvdz,dqtdz,dthldz,&
                          svs,z0,qts,thls,thvs,thlflux,qtflux,svflux
  use modfields, only : u0,v0,thl0,qt0,sv0,u0av,v0av,qt0av,thl0av
  use modmpi,    only :  excj
  implicit none
  integer i, j, n
  real phihzf, phimzf
  real upcu, vpcv
  real horv2, horv, stab, obl
  real dthz1, dqz1, dsvz1
  real, parameter :: C_m = 0.001229, C_h = 0.001094, C_q = 0.001133

!***********************************************************************
!***  Calculate ust, tst, qst and obukhov-length iteratively   *********
!***********************************************************************

  dthz1 = thl0av(1)-thls
  dqz1  = qt0av(1)-qts
  horv = sqrt(u0av(1)**2 + v0av(1)**2)
  horv2 = u0av(1)**2 + v0av(1)**2
  stab  = dthz1+0.61*thvs*dqz1


  do j=2,j1
  do i=2,i1

    dthz1 = thl0(i,j,1) - thls
    dqz1  = qt0(i,j,1)  - qts
    upcu  = 0.5*(u0(i,j,1)+u0(i+1,j,1))+cu
    vpcv  = 0.5*(v0(i,j,1)+v0(i,j+1,1))+cv
    horv  = sqrt(upcu**2 + vpcv**2)
    horv2 = (upcu**2 + vpcv**2)
    stab  = dthz1+0.61*thvs*dqz1

    ustar(i,j) = sqrt(C_m*horv2)
    thlflux(i,j) = - C_h*horv*dthz1
    qtflux(i,j) = -C_q*horv*dqz1
    

    obl   = -ustar(i,j)**3/(fkar*(grav/thvs)*(thlflux(i,j)+0.61*thvs*qtflux(i,j)))

    if (stab < 0.) then
       phimzf = (1.-16.*zf(1)/obl)**(-0.25)
       phihzf = (1.-16.*zf(1)/obl)**(-0.50)
    endif

    if (stab == 0.) then
       phimzf = 1.
       phihzf = 1.
    endif

    if (stab > 0.) then
       phimzf = (1.+5.*zf(1)/obl)
       phihzf = (1.+8.*zf(1)/obl)
    endif
    


    dudz(i,j)   = ustar(i,j)*(phimzf/(fkar*zf(1)))*(upcu/horv)
    dvdz(i,j)   = ustar(i,j)*(phimzf/(fkar*zf(1)))*(vpcv/horv)
    dthldz(i,j) = - thlflux(i,j) / ustar(i,j) * phihzf / (fkar*zf(1))
    dqtdz (i,j) = - qtflux(i,j)  / ustar(i,j) * phihzf / (fkar*zf(1))

  end do
  end do

  do j=2,j1
    ustar(1,j)=ustar(i1,j)
  end do

  call excj( ustar  , 1, i2, 1, j2, 1,1)

  do n=1,nsv
    do j=2,j1
    do i=2,i1
      dsvz1 = sv0(i,j,1,n) - svs(n)
      svflux(i,j,n) = -C_q*horv*dsvz1
    enddo
    enddo
  enddo
end subroutine surf_user

end module moduser
