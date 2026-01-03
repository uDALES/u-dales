!> \file tstep.f90
!!  Performs the time integration

!>
!!  Performs the time integration
!>
!! Tstep uses adaptive timestepping and 3rd order Runge Kutta time integration.
!! The adaptive timestepping chooses it's delta_t according to the courant number
!! and the diffusion number, depending on the advection scheme in use.
!!
!!  \author Jasper Tomas, TU Delft
!!  \author Chiel van Heerwaarden, Wageningen University
!!  \author Thijs Heus,MPI-M
!! \see Wicker and Skamarock 2002
!!  \par Revision list
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

!> Determine time step size dt in initialization and update time variables
!!
!! The size of the timestep Delta t is determined adaptively, and is limited by both the Courant-Friedrichs-Lewy criterion CFL
!! \latexonly
!! \begin{equation}
!! \CFL = \mr{max}\left(\left|\frac{u_i \Delta t}{\Delta x_i}\right|\right),
!! \end{equation}
!! and the diffusion number $d$. The timestep is further limited by the needs of other modules, e.g. the statistics.
!! \endlatexonly


module modtstep
  implicit none
  save
  contains

subroutine tstep_update


  use modglobal, only : ib,ie,jb,je,rk3step,timee,runtime,dtmax,dt,ntimee,ntrun,courant,diffnr,&
                        kb,ke,dx,dxi,dx2i,dyi,dy2i,dzh,dt_lim,ladaptive,timeleft,lwarmstart,&
                        dzh2i,rk3coef,rk3coefi
  use modfields, only : um,vm,wm
  use modsubgriddata, only : ekm,ekh
  use modmpi,    only : myid,comm3d,mpierr,mpi_max,my_real
  implicit none

  integer       :: i, j, k,imin,kmin
  real,save     :: courtot=-1.,diffnrtot=-1.
  real          :: courtotl,courold,diffnrtotl,diffnrold
!  logical,save  :: spinup=.true.
  logical,save  :: spinup=.false.


  if(lwarmstart) spinup = .false.
  rk3step = mod(rk3step,3) + 1
  if(rk3step == 1) then

    ! Initialization
    if (spinup) then
      write(6,*) '!spinup!'
      if (ladaptive) then
        courold = courtot
        diffnrold = diffnrtot
        courtotl=0.
        diffnrtotl = 0.
        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          courtotl = max(courtotl,(abs(um(i,j,k))*dxi + abs(vm(i,j,k))*dyi + abs(wm(i,j,k))/dzh(k))*dt)
!          diffnrtotl = max(diffnrtotl,  ekm(i,j,k)*(1/dzh(k)**2 + dxh2i(i) + dy2i)*dt )
          diffnrtotl = max(diffnrtotl,  ekm(i,j,k)*(dzh2i(k) + dx2i + dy2i)*dt, &
                                        ekh(i,j,k)*(dzh2i(k) + dx2i + dy2i)*dt )
        end do
        end do
        end do
        call MPI_ALLREDUCE(courtotl,courtot,1,MY_REAL,MPI_MAX,comm3d,mpierr)
        call MPI_ALLREDUCE(diffnrtotl,diffnrtot,1,MY_REAL,MPI_MAX,comm3d,mpierr)
        if ( diffnrold>0) then
          dt = min(dtmax,dt*courant/courtot,dt*diffnr/diffnrtot)
          if ((abs(courtot-courold)/courold<0.1) .and. (abs(diffnrtot-diffnrold)/diffnrold<0.1)) then
            spinup = .false.
          end if
        end if
        dt = dt
        dt_lim = timeleft
        timee   = timee  + dt
        timeleft = timeleft- dt
        ntimee  = ntimee + 1
        ntrun   = ntrun  + 1
      else
        dt = 2 * dt
        if (dt >= dtmax) then
          dt = dtmax
          spinup = .false.
        end if
      end if
    ! Normal time loop
    else  !spinup = .false.

      if (ladaptive) then
        courtotl=0.
        diffnrtotl = 1e-5
        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          courtotl = max(courtotl,(abs(um(i,j,k))*dxi + abs(vm(i,j,k))*dyi + abs(wm(i,j,k))/dzh(k))*dt)
          diffnrtotl = max(diffnrtotl,  ekm(i,j,k)*(dzh2i(k) + dx2i + dy2i)*dt,&
                                        ekh(i,j,k)*(dzh2i(k) + dx2i + dy2i)*dt )
!          if (diffnrtotl ==  ekh(i,j,k)*(dzh2i(k) + dxh2i(i) + dy2i)*dt) then
!           imin = i
!           kmin = k
!          end if
        end do
        end do
        end do
!     write(6,*) 'Peclet criterion at proc,i,k = ', myid,imin,kmin

        call MPI_ALLREDUCE(courtotl,courtot,1,MY_REAL,MPI_MAX,comm3d,mpierr)
        call MPI_ALLREDUCE(diffnrtotl,diffnrtot,1,MY_REAL,MPI_MAX,comm3d,mpierr)
        if (courtot <= 0) then
          write(6,*) 'courtot=0!'
        end if
        if (diffnrtot <= 0) then
          write(6,*) 'diffnrtot=0!'
        end if
        dt = min(dtmax,dt*courant/courtot,dt*diffnr/diffnrtot)
        timeleft=timeleft-dt
        dt_lim = timeleft
        timee   = timee  + dt
        ntimee  = ntimee + 1
        ntrun   = ntrun  + 1
      else
        dt = dtmax
        ntimee  = ntimee + 1
        ntrun   = ntrun  + 1
        timee   = timee  + dt
        timeleft=timeleft-dt
      end if
    end if
  end if

  if (rk3step == 0) then ! dt not defined yet
    rk3coef = 1.
  else
    rk3coef = dt / (4. - dble(rk3step))
  end if
  rk3coefi = 1. / rk3coef

end subroutine tstep_update


    !> Time integration is done by a third order Runge-Kutta scheme.
    !!
    !! \latexonly
    !! With $f^n(\phi^n)$ the right-hand side of the appropriate equation for variable
    !! $\phi=\{\fav{u},\fav{v},\fav{w},e^{\smfrac{1}{2}},\fav{\varphi}\}$, $\phi^{n+1}$
    !! at $t+\Delta t$ is calculated in three steps:
    !! \begin{eqnarray}
    !! \phi^{*} &=&\phi^n + \frac{\Delta t}{3}f^n(\phi^n)\nonumber\\\\
    !! \phi^{**} &=&\phi^{n} + \frac{\Delta t}{2}f^{*}(\phi^{*})\nonumber\\\\
    !! \phi^{n+1} &=&\phi^{n} + \Delta t f^{**}(\phi^{**}),
    !! \end{eqnarray}
    !! with the asterisks denoting intermediate time steps.
    !! \endlatexonly
    !! \see Wicker and Skamarock, 2002
    subroutine tstep_integrate
      use modglobal,      only : ib, ie, jb, je, kb, ke, rk3step, rk3coef, e12min, ltempeq, lmoist, nsv, &
                                 BCtopm, BCtopm_pressure, BCxm, BCxm_periodic, BCym, BCym_periodic, ierank, jerank, &
                                 iadv_thl, iadv_kappa, &
                                 dzhi, dyi, dxfi, rk3coefi
      use modsubgriddata, only : loneeqn
      use modchem,        only : chem
      ! use modpois,        only : pij
      ! use modmpi,         only : myid
      use modfields,      only : sv0
#if defined(_GPU)
      use modcuda,        only : up_d, vp_d, wp_d, um_d, vm_d, wm_d, u0_d, v0_d, w0_d, &
                                 e120_d, e12m_d, e12p_d, &
                                 thlp_d, thlm_d, thl0_d, thl0c_d, &
                                 qtp_d, qtm_d, qt0_d, &
                                 svp_d, svm_d, sv0_d
#else
     use modfields,       only : up, vp, wp, um, vm, wm, u0, v0, w0, &
                                 e120, e12m, e12p, &
                                 thlp, thlm, thl0, thl0c, &
                                 qtp, qtm, qt0, &
                                 svp, svm
#endif
      implicit none
      integer :: i, j, k, n

      !$acc kernels default(present)
      do k=kb,ke
        do j=jb,je
          do i=ib,ie
#if defined(_GPU)
            u0_d(i,j,k) = um_d(i,j,k) + rk3coef * up_d(i,j,k)
            v0_d(i,j,k) = vm_d(i,j,k) + rk3coef * vp_d(i,j,k)
            w0_d(i,j,k) = wm_d(i,j,k) + rk3coef * wp_d(i,j,k)
#else
            u0(i,j,k) = um(i,j,k) + rk3coef * up(i,j,k)
            v0(i,j,k) = vm(i,j,k) + rk3coef * vp(i,j,k)
            w0(i,j,k) = wm(i,j,k) + rk3coef * wp(i,j,k)
#endif
          end do
        end do
      end do
      !$acc end kernels

      if ((BCxm .ne. BCxm_periodic) .and. ierank) then
#if defined(_GPU)
        !$acc kernels default(present)
        u0_d(ie+1,jb:je,kb:ke) = um_d(ie+1,jb:je,kb:ke) + rk3coef * up_d(ie+1,jb:je,kb:ke)
        !$acc end kernels
#else
        u0(ie+1,jb:je,kb:ke) = um(ie+1,jb:je,kb:ke) + rk3coef * up(ie+1,jb:je,kb:ke)
#endif
      end if

      if ((BCym .ne. BCym_periodic) .and. jerank) then
#if defined(_GPU)
        !$acc kernels default(present)
        v0_d(ib:ie,je+1,kb:ke) = vm_d(ib:ie,je+1,kb:ke) + rk3coef * vp_d(ib:ie,je+1,kb:ke)
        !$acc end kernels
#else
        v0(ib:ie,je+1,kb:ke) = vm(ib:ie,je+1,kb:ke) + rk3coef * vp(ib:ie,je+1,kb:ke)
#endif
      end if

      if (BCtopm .eq. BCtopm_pressure) then
        ! do i=ib,ie
        !   do j=jb,je
        !     ! w0(i,j,ke+1) = w0(i,j,ke) - dzhi(ke)*((u0(i+1,j,ke)-u0(i,j,ke))*dxfi(i) + &
        !     !                                       (v0(i,j+1,ke)-v0(i,j,ke))*dyi)
        !     ! if (myid ==0 .and. (i==32 .and. j==1)) write(*,*) rk3coefi*(w0(i,j,ke) - dzhi(ke)*((u0(i+1,j,ke)-u0(i,j,ke))*dxfi(i) + &
        !     ! (v0(i,j+1,ke)-v0(i,j,ke))*dyi) - wm(i,j,ke+1)), &
        !     ! dpdztop(i,j), &
        !     ! 2*pij(ke)*dzhi(ke+1)
        !     !
        !     ! wp(i,j,ke+1) = rk3coefi*(w0(i,j,ke) - dzhi(ke)*((u0(i+1,j,ke)-u0(i,j,ke))*dxfi(i) + &
        !     !            (v0(i,j+1,ke)-v0(i,j,ke))*dyi) - wm(i,j,ke+1))
        !     ! wp(i,j,ke+1) = 2*pij(ke)*dzhi(ke+1)
        !   end do
        ! end do
#if defined(_GPU)
        !$acc kernels default(present)
        w0_d(ib:ie,jb:je,ke+1) = wm_d(ib:ie,jb:je,ke+1) + rk3coef * wp_d(ib:ie,jb:je,ke+1)
        !$acc end kernels
#else
        w0(ib:ie,jb:je,ke+1) = wm(ib:ie,jb:je,ke+1) + rk3coef * wp(ib:ie,jb:je,ke+1)
#endif
      end if


      if (loneeqn) then
        !$acc kernels default(present)
        do k=kb,ke
          do j=jb,je
            do i=ib,ie
#if defined(_GPU)
              e120_d(i,j,k) = e12m_d(i,j,k) + rk3coef * e12p_d(i,j,k)
              e120_d(i,j,k) = max(e12min,e120_d(i,j,k))
              e12m_d(i,j,k) = max(e12min,e12m_d(i,j,k))
#else
              e120(i,j,k) = e12m(i,j,k) + rk3coef * e12p(i,j,k)
              e120(i,j,k) = max(e12min,e120(i,j,k))
              e12m(i,j,k) = max(e12min,e12m(i,j,k))
#endif
            end do
          end do
        end do
        !$acc end kernels
      end if


      if (ltempeq) then
        !$acc kernels default(present)
        do k=kb,ke
          do j=jb,je
            do i=ib,ie
#if defined(_GPU)
              thl0_d(i,j,k) = thlm_d(i,j,k) + rk3coef * thlp_d(i,j,k)
#else
              thl0(i,j,k) = thlm(i,j,k) + rk3coef * thlp(i,j,k)
#endif
            end do
          end do
        end do
        !$acc end kernels

        if (iadv_thl == iadv_kappa) then
#if defined(_GPU)
          !$acc kernels default(present)
          thl0c_d(ib:ie,jb:je,kb:ke) = thl0_d(ib:ie,jb:je,kb:ke)
          !$acc end kernels
#else
          thl0c(ib:ie,jb:je,kb:ke) = thl0c(ib:ie,jb:je,kb:ke)
#endif
        end if
      end if


      if (lmoist) then
        !$acc kernels default(present)
        do k=kb,ke
          do j=jb,je
            do i=ib,ie
#if defined(_GPU)
              qt0_d(i,j,k) = qtm_d(i,j,k) + rk3coef * qtp_d(i,j,k)
#else
              qt0(i,j,k) = qtm(i,j,k) + rk3coef * qtp(i,j,k)
#endif
            end do
          end do
        end do
        !$acc end kernels
      end if


      if (nsv>0) then
        !$acc kernels default(present)
        do k=kb,ke
          do j=jb,je
            do i=ib,ie
              do n=1,nsv
#if defined(_GPU)
                sv0_d(i,j,k,n) = svm_d(i,j,k,n) + rk3coef * svp_d(i,j,k,n)
#else
                sv0(i,j,k,n) = svm(i,j,k,n) + rk3coef * svp(i,j,k,n)
#endif
              end do
            end do
          end do
        end do
        !$acc end kernels

        if (nsv==3) then
          call chem
        end if
      end if


#if !defined(_GPU)
      up=0.
      vp=0.
      wp=0.
      thlp=0.
      svp=0.
      e12p=0.
      qtp=0.
#endif

      if(rk3step == 3) then
#if defined(_GPU)
        um_d = u0_d
        vm_d = v0_d
        wm_d = w0_d
        thlm_d = thl0_d
        e12m_d = e120_d
        svm_d = sv0_d
        qtm_d = qt0_d
#else
        um = u0
        vm = v0
        wm = w0
        thlm = thl0
        e12m = e120
        svm = sv0
        qtm = qt0
#endif
      end if

    end subroutine tstep_integrate

end module modtstep