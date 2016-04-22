!> \file modsurface.f90
!!  Surface parameterization
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

!>
!! Surface routine including a land-surface scheme
!>
!! This module provides an interactive surface parameterization
!!
!! \par Revision list
!! \par Chiel van Heerwaarden
!! \todo documentation
!! \todo implement water reservoir at land surface for dew and interception water
!! \todo add moisture transport between soil layers
!!  \deprecated Modsurface replaces the old modsurf.f90

module modsurface
  use modsurfdata
  implicit none
  !public  :: initsurface, surface, exitsurface

  save

contains
  !> Reads the namelists and initialises the soil.
  subroutine initsurface

    use modglobal,  only : ib,ie,jb,je,jmax, ih, jh, cp, rlv, zf, nsv, ifnamopt, fname_options
    use modfields,  only : thl0, qt0
    use modmpi,     only : myid, nprocs, comm3d, mpierr, my_real, mpi_logical, mpi_integer

    implicit none

    integer   ::ierr
    namelist/NAMSURFACE/ & !< Soil related variables
         isurf,&
                                ! Prescribed values for isurf 2, 3, 4
         z0, thvs, thl_top, ps, ustin, wtsurf, wqsurf, wsvsurf

    ! 1    -   Initialize soil

    !if (isurf == 1) then

    ! 1.0  -   Read LSM-specific namelist

    if(myid==0)then
       open(ifnamopt,file=fname_options,status='old',iostat=ierr)
       read (ifnamopt,NAMSURFACE,iostat=ierr)
       if (ierr > 0) then
          print *, 'Problem in namoptions NAMSURFACE'
          print *, 'iostat error: ', ierr
          stop 'ERROR: Problem in namoptions NAMSURFACE'
       endif
       write(6 ,NAMSURFACE)
       close(ifnamopt)
    end if

    call MPI_BCAST(isurf        , 1       , MPI_INTEGER, 0, comm3d, mpierr)
    call MPI_BCAST(z0         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ustin      ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wtsurf     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wqsurf     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wsvsurf(1:nsv),nsv,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ps         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(thvs       ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(thl_top    ,1,MY_REAL   ,0,comm3d,mpierr)
    z0mav = z0
    z0hav = z0

    ! 3. Initialize surface layer
    allocate(ustar   (ib-1:ie+1,jb-1:je+1))

    allocate(z0m(ib-1:ie+1,jb-1:je+1))
    allocate(z0h(ib-1:ie+1,jb-1:je+1))
    allocate(obl(ib-1:ie+1,jb-1:je+1))
    allocate(tskin(ib-1:ie+1,jb-1:je+1))
    allocate(qskin(ib-1:ie+1,jb-1:je+1))
    allocate(Cm(ib-1:ie+1,jb-1:je+1))
    allocate(Cs(ib-1:ie+1,jb-1:je+1))

    allocate(dudz    (ib-1:ie+1,jb-1:je+1))
    allocate(dvdz    (ib-1:ie+1,jb-1:je+1))

    allocate(thlflux (ib-1:ie+1,jb-1:je+1))
    allocate(qtflux  (ib-1:ie+1,jb-1:je+1))
    allocate(dqtdz   (ib-1:ie+1,jb-1:je+1))
    allocate(dthldz  (ib-1:ie+1,jb-1:je+1))
    allocate(svflux  (ib-1:ie+1,jb-1:je+1,nsv))
    allocate(svs(nsv))

    z0m        = z0mav
    z0h        = z0hav
    ustar(:,:)=0   

    return
  end subroutine initsurface

  !> Calculates the interaction with the soil, the surface temperature and humidity, and finally the surface fluxes.
  subroutine surface
    use modglobal,  only : ib,ie,jb,je,ih, jh,kb, cp, rlv, fkar, zf, cu, cv, nsv, rk3step, timee, rslabs, pi, pref0, rd, rv, eps1, lneutral!, boltz, rhow
    use modfields,  only : thl0, qt0, u0, v0, rhof, ql0, exnf, presf, u0av, v0av
    use modmpi,     only : my_real, mpierr, comm3d, mpi_sum, myid, excj, excjs
    implicit none

    integer  :: i, j, n
    real     :: upcu, vpcv, horv, horvav
    real     :: phimzf, phihzf
    real     :: thlsl, qtsl

    real     :: ust,ustl
    real     :: wtsurfl, wqsurfl

    ! 2     -   Calculate the surface fluxes
   ! write(*,*) 'lneutral',lneutral
    if(lneutral) then
       obl(:,:) = -1.e10
       oblav    = -1.e10
    else
       call getobl
    end if

    call MPI_BCAST(oblav ,1,MY_REAL ,0,comm3d,mpierr)

    thlsl = 0.
    qtsl  = 0.

    do j = jb,je
       do i = ib,ie

          upcu   = 0.5 * (u0(i,j,kb) + u0(i+1,j,kb)) + cu
          vpcv   = 0.5 * (v0(i,j,kb) + v0(i,j+1,kb)) + cv
          horv   = sqrt(upcu ** 2. + vpcv ** 2.)
          horv   = max(horv, 0.1)
          horvav = sqrt(u0av(kb) ** 2. + v0av(kb) ** 2.)
          horvav = max(horvav, 0.1)

          if( isurf == 4) then
             if(lmostlocal) then
                ustar (i,j) = fkar * horv  / (log(zf(kb) / z0m(i,j)) - psim(zf(kb) / obl(i,j)) + psim(z0m(i,j) / obl(i,j)))
             else
                ustar (i,j) = fkar * horvav / (log(zf(kb) / z0m(i,j)) - psim(zf(kb) / obl(i,j)) + psim(z0m(i,j) / obl(i,j)))
             end if
          else
             ustar (i,j) = ustin
          end if

          ustar  (i,j) = max(ustar(i,j), 1.e-2)
          thlflux(i,j) = wtsurf
          qtflux (i,j) = wqsurf

          do n=1,nsv
             svflux(i,j,n) = wsvsurf(n)
          enddo
          if (obl(i,j) < 0.) then
             phimzf = (1.-16.*zf(kb)/obl(i,j))**(-0.25)
             !phimzf = (1. + 3.6 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
             phihzf = (1.-16.*zf(kb)/obl(i,j))**(-0.50)
             !phihzf = (1. + 7.9 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
          elseif (obl(i,j) > 0.) then
             phimzf = (1.+5.*zf(kb)/obl(i,j))
             phihzf = (1.+5.*zf(kb)/obl(i,j))
          else
             phimzf = 1.
             phihzf = 1.
          endif

          dudz  (i,j) = ustar(i,j) * phimzf / (fkar*zf(kb))*(upcu/horv)
          dvdz  (i,j) = ustar(i,j) * phimzf / (fkar*zf(kb))*(vpcv/horv)
          dthldz(i,j) = - thlflux(i,j) / ustar(i,j) * phihzf / (fkar*zf(kb))
          dqtdz (i,j) = - qtflux(i,j)  / ustar(i,j) * phihzf / (fkar*zf(kb))
          Cs(i,j) = fkar ** 2. / (log(zf(kb) / z0m(i,j)) - psim(zf(kb) / obl(i,j)) + psim(z0m(i,j) / obl(i,j))) / (log(zf(kb) / z0h(i,j)) - psih(zf(kb) / obl(i,j)) + psih(z0h(i,j) / obl(i,j)))

          tskin(i,j) = wtsurf / (Cs(i,j) * horv) + thl0(i,j,kb)
          qskin(i,j) = wqsurf / (Cs(i,j) * horv) + qt0(i,j,kb)
          thlsl      = thlsl + tskin(i,j)
          qtsl       = qtsl  + qskin(i,j)
       end do
    end do
   

  !     write(*,'(1A,F8.4)') 'upcu', upcu
  !     write(*,'(1A,F8.4)') 'vpcv', vpcv
!       write(*,'(1A,F8.4)') 'phimzf', phimzf
  !     write(*,'(1A,F8.4)') 'horv', horv
  !     write(*,'(1A,F8.4)') 'fkar', fkar
  !     write(*,'(1A,F8.4)') 'zf', zf(kb)
  !     write(*,'(1A,F8.4)') 'obl', obl(2,3)
  !     write(*,'(1A,F8.4)') 'ustar', ustar(2,3)
 !      write(*,'(1A,F8.4)') 'dudz', dudz(2,3)
 
    call MPI_ALLREDUCE(thlsl, thls, 1,  MY_REAL, MPI_SUM, comm3d,mpierr)
    call MPI_ALLREDUCE(qtsl, qts, 1,  MY_REAL, MPI_SUM, comm3d,mpierr)

    thls = thls / rslabs
    qts  = qts  / rslabs
    thvs = thls * (1. + (rv/rd - 1.) * qts)

    !call qtsurf

    ! Transfer ustar to neighbouring cells
    do j=jb-1,je+1
       ustar(ib-1,j)=ustar(ie,j)
       ustar(ie+1,j)=ustar(ib,j)
    end do

    call excj( ustar,ib-1,ie+1,jb-1,je+1,kb,kb)

    return

  end subroutine surface

subroutine getobl
    use modglobal,      only : ib,ie,jb,je,kb,ke,ih,jh,kh,zf,dzh,dzf,dxfi,dxf,xh,dxhi,lMOST,fkar,grav,jgb,jge,lbuoyancy
    use modfields,      only : thl0, u0, v0
    use modsurfdata,    only : thls, oblav, z0mav,z0hav,z0,Cmav,Csav,horvel
    use modmpi,         only : slabsumi,myid
    implicit none

    real, dimension(ib:ke) :: u1avdum
    real, dimension(ib:ie) :: v1avdum
    real, dimension(ib:ie) :: t1avdum
    real :: u1av, v1av,t1av,horv2,Rib,L,Lold,Lstart,Lend,fx,fxdif
    integer :: iter

    if (.not. (lMOST)) return


    ! determine average velocity at temperature at first level
    u1avdum = 0.
    v1avdum = 0.
    t1avdum = 0.
    call slabsumi(u1avdum,ib,ie,u0,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,ib,ie,jb,je,kb,kb)
    call slabsumi(v1avdum,ib,ie,v0,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,ib,ie,jb,je,kb,kb)
    call slabsumi(t1avdum,ib,ie,thl0,ib-1,ie+1,jb-1,je+1,kb-1,ke+1,ib,ie,jb,je,kb,kb)
    u1av = sum(u1avdum)/ ((ie-ib+1)*(jge-jgb+1))
    v1av = sum(v1avdum)/ ((ie-ib+1)*(jge-jgb+1))
    t1av = sum(t1avdum)/ ((ie-ib+1)*(jge-jgb+1))


    z0mav = z0
    z0hav = z0
     
    horv2 = u1av**2. + v1av**2.
    horv2 = max(horv2, 0.01)
    horvel = sqrt(horv2)         ! Average horizontal velocity at first level

    if (lbuoyancy == .false.) then
      oblav = -1.e10
    else   ! iterate to find L
      Rib   = grav / thls * zf(kb) * (t1av - thls) / horv2
!      write(6,*) 'horvel, u1av, t1av, Rib = ', horvel, u1av, t1av, Rib

      iter = 0
      L = oblav

      if(Rib * L < 0. .or. abs(L) == 1e5) then
        if(Rib > 0) L = 0.01
!        if(Rib < 0) L = -0.01
        if(Rib <= 0) L = -0.01
      end if

      do while (.true.)
        iter    = iter + 1
        Lold    = L
        fx      = Rib - zf(kb) / L * (log(zf(kb) / z0hav) - psih(zf(kb) / L) +psih(z0hav / L)) / (log(zf(kb) / z0mav) - psim(zf(kb) / L) + psim(z0mav / L)) **2.
        Lstart  = L - 0.001*L
        Lend    = L + 0.001*L
        fxdif   = ( (- zf(kb) / Lstart * (log(zf(kb) / z0hav) - psih(zf(kb) /Lstart) + psih(z0hav / Lstart)) / (log(zf(kb) / z0mav) - psim(zf(kb) / Lstart) +psim(z0mav / Lstart)) ** 2.) - (-zf(kb) / Lend * (log(zf(kb) / z0hav) -psih(zf(kb) / Lend) + psih(z0hav / Lend)) / (log(zf(kb) / z0mav) - psim(zf(kb) /Lend) + psim(z0mav / Lend)) ** 2.) ) / (Lstart - Lend)
        if (fxdif /= 0.) then
          L       = L - fx / fxdif
        end if
        if(Rib * L < 0. .or. abs(L) == 1e5) then
          if(Rib > 0) L = 0.01
!          if(Rib < 0) L = -0.01
          if(Rib <= 0) L = -0.01
        end if
!        if(abs((L - Lold) / Lold) < 1e-4) exit
!        if(abs((L - Lold) / Lold) < 1e-4) then
        if(abs((L - Lold) / Lold) < 1e-5) then
!          write(6,*) 'iter = ', iter
          exit
        end if
        if(iter > 1000) stop 'Obukhov length calculation does not converge!'
      end do

      if(L > 1e6)  L = 1e6
      if(L < -1e6) L = -1e6

      oblav = L
    end if ! lbuoyancy == .true.   

!   compute drag coefficients for momentum and heat
    Cmav = fkar ** 2. / (log(zf(kb) / z0mav) - psim(zf(kb) / oblav) + psim(z0mav/ oblav)) ** 2.
    Csav = fkar ** 2. / (log(zf(kb) / z0mav) - psim(zf(kb) / oblav) + psim(z0mav/ oblav)) / (log(zf(kb) / z0hav) - psih(zf(kb) / oblav) + psih(z0hav / oblav))

    if (myid ==0) then
      write(6,*) 'Obukhov length, Cmav, Csav = ', oblav, Cmav, Csav
    end if   

  end subroutine getobl




  !> Calculates the Obukhov length iteratively.
!  subroutine getobl
!    use modglobal, only : kb,zf, rv, rd, grav, rslabs, timee, cu, cv
!    use modfields, only : thl0av, qt0av, u0, v0, thl0, qt0, u0av, v0av
!    use modmpi,    only : my_real,mpierr,comm3d,mpi_sum,myid,excj
!    implicit none
!
 !   integer             :: i,j,iter
!    real                :: thv, thvsl, L, horv2, oblavl
!    real                :: Rib, Lstart, Lend, fx, fxdif, Lold
!    real                :: upcu, vpcv
!
!    thv    = thl0av(kb) * (1. + (rv/rd - 1.) * qt0av(kb))
!
!    horv2 = u0av(kb)**2. + v0av(kb)**2.
!    horv2 = max(horv2, 0.01)
!
!    Rib   = grav / thvs * zf(kb) * (thv - thvs) / horv2


   !    write(*,'(1A,F8.4)') 'Rib', Rib
  !     write(*,'(1A,F8.4)') 'grav', grav
 !      write(*,'(1A,F8.4)') 'thvs', thvs
    !   write(*,'(1A,F8.4)') 'thv', thv
       !write(*,'(1A,F8.4)') 'horv2', horv2



 !   iter =  0
 !   L = oblav
    !write(*,*) 'test why divided by zero'
    !write(*,*) 'oblav',oblav
    !write(*,'(1A,F8.4)') 'L', L
       !write(*,'(1A,F8.4)') 'z0hav', z0hav
       !write(*,'(1A,F8.4)') 'z0mav', z0mav


 !   if(Rib * L < 0. .or. abs(L) == 1e5) then
 !      if(Rib > 0) L = 0.01
 !      if(Rib < 0) L = -0.01
 !   end if

 !   do while (.true.)
 !      iter    = iter + 1
 !      Lold    = L
 !      fx      = Rib - zf(kb) / L * (log(zf(kb) / z0hav) - psih(zf(kb) / L) + psih(z0hav / L)) / (log(zf(kb) / z0mav) - psim(zf(kb) / L) + psim(z0mav / L)) ** 2.
 !      Lstart  = L - 0.001*L
 !      Lend    = L + 0.001*L
 !     !write(*,*) 'L',L
       !write(*,*) 'z0hav',z0hav
       !write(*,*) 'z0mav',z0mav
       !write(*,*) 'zf(kb)',zf(kb)
  !    fxdif   = ( (- zf(kb) / Lstart * (log(zf(kb) / z0hav) - psih(zf(kb) / Lstart) + psih(z0hav / Lstart)) / (log(zf(kb) / z0mav) - psim(zf(kb) / Lstart) + psim(z0mav / Lstart)) ** 2.) - (-zf(kb) / Lend * (log(zf(kb) / z0hav) - psih(zf(kb) / Lend) + psih(z0hav / Lend)) / (log(zf(kb) / z0mav) - psim(zf(kb) / Lend) + psim(z0mav / Lend)) ** 2.) ) / (Lstart - Lend)
 !      L       = L - fx / fxdif
 !      if(Rib * L < 0. .or. abs(L) == 1e5) then
 !         if(Rib > 0) L = 0.01
 !         if(Rib < 0) L = -0.01
!       end if
!       if(abs((L - Lold) / Lold) < 1e-4) exit
!       if(iter > 1000) stop 'Obukhov length calculation does not converge!'
!    end do
!
!    if(L > 1e6)  L = 1e6
!    if(L < -1e6) L = -1e6
!
!    obl(:,:) = L
!    oblav = L
!
!    return
!
!  end subroutine getobl

  function psim(zeta)
    implicit none

    real             :: psim
    real, intent(in) :: zeta
    real             :: x

    if(zeta <= 0) then
       x     = (1. - 16. * zeta) ** (0.25)
       psim  = 3.14159265 / 2. - 2. * atan(x) + log( (1.+x) ** 2. * (1. + x ** 2.) / 8.)
       ! CvH use Wilson, 2001 rather than Businger-Dyer for correct free convection limit
       !x     = (1. + 3.6 * abs(zeta) ** (2./3.)) ** (-0.5)
       !psim = 3. * log( (1. + 1. / x) / 2.)
    else
       psim  = -2./3. * (zeta - 5./0.35)*exp(-0.35 * zeta) - zeta - (10./3.) / 0.35
    end if

    return
  end function psim

  function psih(zeta)

    implicit none

    real             :: psih
    real, intent(in) :: zeta
    real             :: x

    if(zeta <= 0) then
       x     = (1. - 16. * zeta) ** (0.25)
       psih  = 2. * log( (1. + x ** 2.) / 2. )
       ! CvH use Wilson, 2001
       !x     = (1. + 7.9 * abs(zeta) ** (2./3.)) ** (-0.5)
       !psih  = 3. * log( (1. + 1. / x) / 2.)
    else
       psih  = -2./3. * (zeta - 5./0.35)*exp(-0.35 * zeta) - (1. + (2./3.) * zeta) ** (1.5) - (10./3.) / 0.35 + 1.
    end if

    return
  end function psih

  subroutine exitsurface
    implicit none
    return
  end subroutine exitsurface

end module modsurface
