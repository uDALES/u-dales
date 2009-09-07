!> \file modfields.f90
!!  Declares, allocates and initializes the 3D fields

!>
!!  Declares, allocates and initializes the 3D fields
!>

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

module modfields

implicit none
save

  ! Prognostic variables

  real, allocatable :: um(:,:,:)        !<   x-component of velocity at time step t-1
  real, allocatable :: vm(:,:,:)        !<   y-component of velocity at time step t-1
  real, allocatable :: wm(:,:,:)        !<   z-component of velocity at time step t-1
  real, allocatable :: thlm(:,:,:)      !<   liq. water pot. temperature at time step t-1
  real, allocatable :: e12m(:,:,:)      !<   turb. kin. energy at time step t-1
  real, allocatable :: qtm(:,:,:)       !<   total specific humidity at time step t
  real, allocatable :: u0(:,:,:)        !<   x-component of velocity at time step t
  real, allocatable :: v0(:,:,:)        !<   y-component of velocity at time step t
  real, allocatable :: w0(:,:,:)        !<   z-component of velocity at time step t
  real, allocatable :: thl0(:,:,:)      !<   liq. water pot. temperature at time step t
  real, allocatable :: thl0h(:,:,:)     !<  3d-field of theta_l at half levels for kappa scheme
  real, allocatable :: qt0h(:,:,:)      !<  3d-field of q_tot   at half levels for kappa scheme
  real, allocatable :: e120(:,:,:)      !<   turb. kin. energy at time step t
  real, allocatable :: qt0(:,:,:)       !<   total specific humidity at time step t

  real, allocatable :: up(:,:,:)        !<   tendency of um
  real, allocatable :: vp(:,:,:)        !<   tendency of vm
  real, allocatable :: wp(:,:,:)        !<   tendency of wm
  real, allocatable :: thlp(:,:,:)      !<   tendency of thlm
  real, allocatable :: e12p(:,:,:)      !<   tendency of e12m
  real, allocatable :: qtp(:,:,:)       !<   tendency of qtm

  real, allocatable :: svm(:,:,:,:)   !<  scalar sv(n) at time step t-1
  real, allocatable :: sv0(:,:,:,:)   !<  scalar sv(n) at time step t
  real, allocatable :: svp(:,:,:,:)   !<  tendency of sv(n)



  ! Diagnostic variables

  real, allocatable :: ql0(:,:,:)  !<   liquid water content

  real, allocatable :: thv0h(:,:,:)!<   theta_v at half level

  real, allocatable :: whls(:)                       !<   large scale vert velocity at half levels

  real, allocatable :: presf(:)                      !<   hydrostatic pressure at full level
  real, allocatable :: presh(:)                      !<   hydrostatic pressure at half level
  real, allocatable :: exnf(:)                       !<   hydrostatic exner function at full level
  real, allocatable :: exnh(:)                       !<   hydrostatic exner function at half level

  real, allocatable :: rhof(:)                       !<   slab averaged density at full level
  real, allocatable :: qt0av(:)                      !<   slab averaged q_tot
  real, allocatable :: ql0av(:)                      !<   slab averaged q_liq

  real, allocatable :: thl0av(:)                     !<   slab averaged th_liq
  real, allocatable :: u0av(:)                       !<   slab averaged u
  real, allocatable :: v0av(:)                       !<   slab averaged v
  real, allocatable :: ug(:)                       !<   geostrophic u-wind
  real, allocatable :: vg(:)                       !<   geostrophic v-wind

  real, allocatable :: dpdxl(:)                      !<   large scale pressure x-gradient
  real, allocatable :: dpdyl(:)                      !<   large scale pressure y-gradient

  real, allocatable :: dthldxls(:)                   !<   large scale x-gradient of th_liq
  real, allocatable :: dthldyls(:)                   !<   large scale y-gradient of th_liq
  real, allocatable :: dqtdxls(:)                    !<   large scale x-gradient of q_tot
  real, allocatable :: dqtdyls(:)                    !<   large scale y-gradient of q_tot
  real, allocatable :: dqtdtls(:)                    !<   large scale y-gradient of q_tot
  real, allocatable :: dudxls(:)                     !<   large scale x-gradient of u

  real, allocatable :: dudyls(:)                     !<   large scale y-gradient of u
  real, allocatable :: dvdxls(:)                     !<   large scale x-gradient of v
  real, allocatable :: dvdyls(:)                     !<   large scale y-gradient of v
  real, allocatable :: wfls  (:)                     !<   large scale y-gradient of v
  real, allocatable :: ql0h(:,:,:)
  real, allocatable :: dthvdz(:,:,:)!<   theta_v at half level

  real, allocatable :: thlprof(:)                    !<   initial thl-profile
  real, allocatable :: qtprof(:)                     !<   initial qt-profile
  real, allocatable :: uprof(:)                      !<   initial u-profile
  real, allocatable :: vprof(:)                      !<   initial v-profile
  real, allocatable :: e12prof(:)                    !<   initial subgrid TKE profile
  real, allocatable :: sv0av(:,:)                  !<   slab average of sv(n)
  real, allocatable :: svprof(:,:)                 !<   initial sv(n)-profile

  real, allocatable :: thlpcar(:)                    !< prescribed radiatively forced thl tendency

contains
!> Allocate and initialize the prognostic variables
  subroutine initfields

    use modglobal, only : i1,ih,j1,jh,kmax,k1,nsv
    ! Allocation of prognostic variables
    implicit none

    allocate(um(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(vm(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(wm(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thlm(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(e12m(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(qtm(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(u0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(v0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(w0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thl0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thl0h(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(qt0h(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(e120(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(qt0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(up(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(vp(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(wp(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thlp(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(e12p(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(qtp(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(svm(2-ih:i1+ih,2-jh:j1+jh,k1,nsv))
    allocate(sv0(2-ih:i1+ih,2-jh:j1+jh,k1,nsv))
    allocate(svp(2-ih:i1+ih,2-jh:j1+jh,k1,nsv))

    ! Allocation of diagnostic variables
    allocate(ql0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thv0h(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(whls(k1))
    allocate(presf(k1))
    allocate(presh(k1))
    allocate(exnf(k1))
    allocate(exnh(k1))
    allocate(rhof(k1))
    allocate(qt0av(k1))
    allocate(ql0av(k1))
    allocate(thl0av(k1))
    allocate(u0av(k1))
    allocate(v0av(k1))
    allocate(ug(k1))
    allocate(vg(k1))
    allocate(dpdxl(k1))
    allocate(dpdyl(k1))
    allocate(dthldxls(k1))
    allocate(dthldyls(k1))
    allocate(dqtdxls(k1))
    allocate(dqtdyls(k1))
    allocate(dqtdtls(k1))
    allocate(dudxls(k1))
    allocate(dudyls(k1))
    allocate(dvdxls(k1))
    allocate(dvdyls(k1))
    allocate(wfls  (k1))
    allocate(ql0h(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(dthvdz(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thlprof(k1))
    allocate(qtprof(k1))
    allocate(uprof(k1))
    allocate(vprof(k1))
    allocate(e12prof(k1))
    allocate(sv0av(k1,nsv))
    allocate(svprof(k1,nsv))

    allocate(thlpcar(k1))

    um=0;u0=0;up=0
    vm=0;v0=0;vp=0
    wm=0;w0=0;wp=0
    thlm=0.;thl0=0.;thlp=0
    qtm=0;qt0=0;qtp=0
    e12m=0;e120=0;e12p=0
    svm=0;sv0=0;svp=0

    ql0=0;thv0h=0;thl0h=0;qt0h=0
    presf=0;presh=0;exnf=0;exnh=0;rhof=0    ! OG
    qt0av=0;ql0av=0;thl0av=0;u0av=0;v0av=0;sv0av=0
    thlprof=0;qtprof=0;uprof=0;vprof=0;e12prof=0;svprof=0
    ug=0;vg=0;dpdxl=0;dpdyl=0;wfls=0;whls=0;thlpcar = 0
    dthldxls=0;dthldyls=0;dqtdxls=0;dqtdyls=0;dudxls=0;dudyls=0;dvdxls=0;dvdyls=0
    dthvdz=0

  end subroutine initfields

!> Deallocate the fields
  subroutine exitfields
  implicit none

    deallocate(um,vm,wm,thlm,e12m,qtm,u0,v0,w0,thl0,thl0h,qt0h,e120,qt0)
    deallocate(up,vp,wp,thlp,e12p,qtp)
    deallocate(svm,sv0,svp)
    deallocate(ql0,ql0h,thv0h,dthvdz,whls,presf,presh,exnf,exnh,rhof,qt0av,ql0av,thl0av,u0av,v0av)
    deallocate(ug,vg,dpdxl,dpdyl,dthldxls,dthldyls,dqtdxls,dqtdyls,dqtdtls,dudxls,dudyls,dvdxls,dvdyls,wfls)
    deallocate(thlprof,qtprof,uprof,vprof,e12prof,sv0av,svprof)
    deallocate(thlpcar)

   end subroutine exitfields

end module modfields
