!> \file modstatsdump.f90
!! Dumps statistics of various fields
!>
!
!! \author Tom Grylls, ICL, May 2016
!! \par Revision list
!!   Dipanjan Majumdar, ICL (2024-2025)
!! \todo Documentation
!
! This file is part of uDALES (https://github.com/uDALES/u-dales).
!
! uDALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! uDALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) 2016- the uDALES Team, Imperial College London.
!
module modstatsdump

  use mpi
  use modglobal, only : dt,ltkedump,lmintdump
  implicit none
  private
  PUBLIC :: initstatsdump,statsdump,exitstatsdump
  save

  !NetCDF variables
  integer :: ncidtke,nrectke=0,nstattke=8,&
             ncidmint,nrecmint=0,nstatmint=6
  character(80) :: tkename = 'tkedump.xxx.nc'
  character(80) :: mintname = 'mintdump.xxx.xxx.xxx.nc'
  character(80),dimension(1,4) :: tncstattke
  character(80),dimension(1,4) :: tncstatmint

  integer :: klow,khigh
  real    :: tsamplep,tstatsdumpp

  contains

  subroutine initstatsdump
    use modmpi,     only : myid,cmyidx,cmyidy
    use modglobal,  only : imax,jmax,cexpnr,kb,ke
    use modstat_nc, only : ncinfo,open_nc,define_nc,writestat_dims_nc
    use modfields,  only : ncstattke,ncstatmint
    implicit none

    allocate(ncstattke(nstattke,4))
    allocate(ncstatmint(nstatmint,4))

    klow=kb
    khigh=ke

    if (lmintdump) then
      mintname(10:12) = cmyidx
      mintname(14:16) = cmyidy
      mintname(18:20) = cexpnr
      call ncinfo(tncstatmint(1,:),'time'      ,'Time'                        ,'s'      ,'time')
      call ncinfo(ncstatmint( 1,:),'ut'        ,'Streamwise velocity'         ,'m/s'    ,'mttt'  )
      call ncinfo(ncstatmint( 2,:),'vt'        ,'Spanwise velocity'           ,'m/s'    ,'tmtt'  )
      call ncinfo(ncstatmint( 3,:),'wt'        ,'Vertical velocity'           ,'m/s'    ,'ttmt'  )
      call ncinfo(ncstatmint( 4,:),'thlt'      ,'Temperature'                 ,'K'      ,'tttt'  )
      call ncinfo(ncstatmint( 5,:),'qtt'       ,'Moisture'                    ,'kg/kg'  ,'tttt'  )
      call ncinfo(ncstatmint( 6,:),'pt'        ,'Pressure'                    ,'m^2/s^2','tttt'  )

      call open_nc(mintname, ncidmint, nrecmint, n1=imax, n2=jmax, n3=khigh-klow+1)
      if (nrecmint==0) then
        call define_nc( ncidmint, 1, tncstatmint)
        call writestat_dims_nc(ncidmint)
      end if
      call define_nc( ncidmint, nstatmint, ncstatmint)
    end if

    !> Generate time, y and x averaged NetCDF for tke budget: tkedump.xxx.nc
    if (ltkedump) then
      tkename(9:11) = cexpnr
      call ncinfo(tncstattke(1,:),'time' ,'Time'                                 ,'s'       ,'time')
      call ncinfo(ncstattke( 1,:),'p_b'  ,'p_bant production or consumption term', 'm^2/s^3','tt'  )
      call ncinfo(ncstattke( 2,:),'t_p'  ,'total viscous transport (?)'          , 'm^2/s^3','tt'  )
      call ncinfo(ncstattke( 3,:),'adv'  ,'Advection by mean wind'               , 'm^2/s^3','tt'  )
      call ncinfo(ncstattke( 4,:),'t_t'  ,'Total turb???'                        , 'm^2/s^3','tt'  )
      call ncinfo(ncstattke( 5,:),'t_sgs','total SGS  term'                      , 'm^2/s^3','tt'  )
      call ncinfo(ncstattke( 6,:),'p_t'  ,'Shear production term'                , 'm^2/s^3','tt'  )
      call ncinfo(ncstattke( 7,:),'t_v'  ,'Resolved viscous dissipation term'    , 'm^2/s^3','tt'  )
      call ncinfo(ncstattke( 8,:),'d_sgs','SGS dissipation term'                 , 'm^2/s^3','tt'  )

      if (myid==0) then
        call open_nc(tkename, ncidtke, nrectke, n3=khigh-klow+1)
        if (nrectke==0) then
          call define_nc( ncidtke, 1, tncstattke)
          call writestat_dims_nc(ncidtke)
        end if
        call define_nc( ncidtke, nstattke, ncstattke)
      endif !myid==0
    endif

    !> Set times to zero so works for warm starts... could have issues with warmstarts here...
    tsamplep = 0.
    tstatsdumpp = 0.

  end subroutine initstatsdump

  !-------------------------
  !> Generate and write statistics into NetCDF file format
  !-------------------------

  subroutine statsdump
    use modfields,  only : um,vm,wm,qtm,thlm,pres0,ncstattke,ncstatmint,&
                           t_t,t_v,t_p,t_sgs,d_sgs,p_b,p_t,adv,&
                           umt,vmt,wmt,thlt,qtt,pt
    use modglobal,  only : ib,ie,jb,je,kb,ke,kh,imax,jmax,rk3step,&
                           timee,tsample,tstatsdump,tstatstart,&
                           ltempeq,lmoist
    use modstat_nc, only : writestat_nc,writestat_1D_nc
    use modmpi,     only : myid
    implicit none

    real, allocatable :: varstke(:,:), varsmint(:,:,:,:)
    real              :: tstatsdumppi

    if (timee < tstatstart) return

    if (.not.(lmintdump .or. ltkedump)) return

    if (.not. rk3step==3)  return

    if (tsamplep > tsample) then

      if (lmintdump) then

        tstatsdumppi = 1./tstatsdumpp

        !!>> CALCS FOR TIME DEPENDANT (AVERAGED) STATS
        ! Average 3-D fields in time
        umt(ib:ie,jb:je,kb:ke+kh) = (umt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + um(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        vmt(ib:ie,jb:je,kb:ke+kh) = (vmt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + vm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        wmt(ib:ie,jb:je,kb:ke+kh) = (wmt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        pt(ib:ie,jb:je,kb:ke+kh) = (pt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + pres0(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi

        if (ltempeq) then
          thlt(ib:ie,jb:je,kb:ke+kh) = (thlt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + thlm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        end if

        if (lmoist) then
          qtt(ib:ie,jb:je,kb:ke+kh) = (qtt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + qtm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        end if

      end if ! lmintdump

      tsamplep = dt
    else !timestatsdumpp < tsample

      tsamplep = tsamplep + dt

    endif

    if (tstatsdumpp > tstatsdump) then

      ! Write t-averaged statistics every tsample
      if (lmintdump) then
        allocate(varsmint(imax,jmax,khigh-klow+1,nstatmint))
        call writestat_nc(ncidmint,1,tncstatmint,(/timee/),nrecmint,.true.)
        varsmint(:,:,:,1)  = umt(ib:ie,jb:je,kb:ke)
        varsmint(:,:,:,2)  = vmt(ib:ie,jb:je,kb:ke)
        varsmint(:,:,:,3)  = wmt(ib:ie,jb:je,kb:ke)
        varsmint(:,:,:,4)  = thlt(ib:ie,jb:je,kb:ke)
        varsmint(:,:,:,5)  = qtt(ib:ie,jb:je,kb:ke)
        varsmint(:,:,:,6)  = pt(ib:ie,jb:je,kb:ke)
        call writestat_nc(ncidmint,nstatmint,ncstatmint,varsmint,nrecmint,imax,jmax,khigh-klow+1)
        deallocate(varsmint)
      end if !lmintdump

      if (ltkedump) then
        call tkestatsdump
        if (myid == 0) then
          call writestat_nc(ncidtke,1,tncstattke,(/timee/),nrectke,.true.)
          allocate(varstke(khigh-klow+1,nstattke))
          varstke(:,1) = p_b(kb:ke+kh)
          varstke(:,2) = t_p(kb:ke+kh)
          varstke(:,3) = adv(kb:ke+kh)
          varstke(:,4) = t_t(kb:ke+kh)
          varstke(:,5) = t_sgs(kb:ke+kh)
          varstke(:,6) = p_t(kb:ke+kh)
          varstke(:,7) = t_v(kb:ke+kh)
          varstke(:,8) = d_sgs(kb:ke+kh)
          call writestat_1D_nc(ncidtke,nstattke,ncstattke,varstke,nrectke,khigh-klow+1)
        end if !myid
      endif !ltkedump

      tstatsdumpp = dt

    else !tstatsdumpp < tstatsdump

      tstatsdumpp = tstatsdumpp + dt

    endif

  end subroutine statsdump

  !> tg3315 still under going work to be completed
  subroutine tkestatsdump

  use modfields,        only : u0,v0,w0,thl0,uav,vav,wav,uuav,vvav,wwav,uvav,uwav,vwav,thlav,thlthlav,pres0,thluav,thlvav,thlwav,&
                               upupav,vpvpav,wpwpav,thlpthlpav,upvpav,upwpav,vpwpav,thlpupav,thlpvpav,thlpwpav,presav,&
                               strain2av,disssgsav,t_vav,tvmx,tvmy,tvmz,tsgsmx1,tsgsmx2,tsgsmy1,tsgsmy2,tsgsmz1,t_sgsav,nusgsav,&
                               tpm,t_pav,ttmx,ttmy,ttmz,t_tav,p_bav,d_sgsav,p_tav,tkeadv,tsgsmz1,tsgsmz2,t_t,t_v,t_p,t_sgs,d_sgs,&
                               p_b,p_t,adv,IIc,IIcs
  use modglobal,        only : ib,ie,ih,jb,je,jgb,jge,dy,jh,ke,kb,kh,rk3step,timee,cexpnr,tsample,tstatsdump,jtot,imax,dzf,&
                               dzf,dzfi,dzhi,dxf,dxfi,dyi,dxhi,dy2i,grav,numol,ierank,jerank
  use modmpi,           only : myid,cmyid,my_real,mpi_sum,avey_ibm,mpierr,comm3d,excjs,avexy_ibm
  use modsurfdata,      only : thls
  use decomp_2d,        only : exchange_halo_z
  implicit none

  real, dimension(ib:ie,jb:je,kb:ke)  :: disssgsfl     ! average subgrid visc. * average rate of strain squared : 2*<nu_t>*<Sij>*<Sij>
  real, dimension(ib:ie,jb:je,kb:ke)  :: dissresav     ! average resolved dissipation: 2*nu*<Sij'*Sij'> = 2*nu*( <Sij*Sij> - <Sij>*<Sij> )
    real, dimension(ib:ie,jb:je,kb:ke)  :: tke           ! tke = 0.5*<ui'ui'>
    real, dimension(ib:ie,jb:je,kb:ke)  :: mke           ! = <ui>d/dxj(<ui><uj>) + <ui>d/dxj(<ui'uj'>) = <ui>d/dxj(<ui*uj>)
    real, dimension(ib-1:ie+1,jb-1:je+1,kb:ke)    :: dummyx
    real, dimension(ib-1:ie+1,jb-1:je+1,kb:ke)    :: dummyy
    real, dimension(ib:ie,  jb  :je,  kb:ke+1)  :: dummyz

  integer i,j,k,ip,im,jp,jm,kp,km
  real strainav2
  real dummy

    ! Tvav = (Tvm - <ui>*d/dxj(<Sij>)  ) + 2*nu*<Sij'Sij'>
    ! Tvm = Tvmx + Tvmy + Tvmz -> therefore: subtraction, then interpolation,
    ! then addition of 2*nu*<Sij'Sij'>
    do k=kb,ke
      km = k-1
      kp = k+1
      do j=jb,je
        jp = j+1
        jm = j-1
        do i=ib,ie
          im = i-1
          ip = i+1

!            t_vav(i,j,k) =  0.5*( (tvmx(i,j,k) - (                      &
             dummyx(i,j,k) =  (                      &
                              ( numol  * (uav(i+1,j,k)-uav(i,j,k))*dxfi(i) &
                              -numol * (uav(i,j,k)-uav(i-1,j,k))*dxfi(i-1) ) * 2. * dxhi(i) &
                              + &
                              ( numol * ( (uav(i,jp,k)-uav(i,j,k))   *dyi &
                              +(vav(i,jp,k)-vav(i-1,jp,k))*dxhi(i)) &
                              - numol * ( (uav(i,j,k)-uav(i,jm,k))   *dyi &
                              +(vav(i,j,k)-vav(i-1,j,k))  *dxhi(i)) &
                              ) * dyi &
                              + &
                              ( numol * ( (uav(i,j,kp)-uav(i,j,k))   *dzhi(kp) &
                              +(wav(i,j,kp)-wav(i-1,j,kp))*dxhi(i)) &
                              - numol * ( (uav(i,j,k)-uav(i,j,km))   *dzhi(k) &
                              +(wav(i,j,k)-wav(i-1,j,k))  *dxhi(i)) &
                              ) *dzfi(k) )
               ! y-direction
               dummyy(i,j,k) =  (                        &
                                ( numol * ( (vav(i+1,j,k)-vav(i,j,k))   *dxhi(i+1) &
                                +(uav(i+1,j,k)-uav(i+1,jm,k))*dyi) &
                                -numol * ( (vav(i,j,k)-vav(i-1,j,k))   *dxhi(i) &
                                +(uav(i,j,k)-uav(i,jm,k))    *dyi) &
                                ) * dxfi(i) &        ! = d/dx( Km*(dv/dx + du/dy) )
                                + &
                                (numol * (vav(i,jp,k)-vav(i,j,k)) &
                                -numol * (vav(i,j,k)-vav(i,jm,k))  ) * 2. * dy2i &        ! =d/dy( 2*Km*(dv/dy) )
                                + &
                                ( numol * ( (vav(i,j,kp)-vav(i,j,k))    *dzhi(kp) &
                                +(wav(i,j,kp)-wav(i,jm,kp))  *dyi) &
                                -numol * ( (vav(i,j,k)-vav(i,j,km))    *dzhi(k) &
                                +(wav(i,j,k)-wav(i,jm,k))    *dyi)   &
                                ) * dzfi(k) )                    ! = d/dz( Km*(dv/dz + dw/dy) )
               ! z-direction
               dummyz(i,j,k) = (                        &
                              ( numol * ( (wav(i+1,j,k)-wav(i,j,k))    *dxhi(i+1) &
                              +(uav(i+1,j,k)-uav(i+1,j,km)) *dzhi(k) ) &
                              -numol * ( (wav(i,j,k)-wav(i-1,j,k))    *dxhi(i) &
                              +(uav(i,j,k)-uav(i,j,km))     *dzhi(k) ) &
                             )*dxfi(i) &
                             + &
                             ( numol * ( (wav(i,jp,k)-wav(i,j,k))     *dyi &
                             +(vav(i,jp,k)-vav(i,jp,km))   *dzhi(k) ) &
                             -numol * ( (wav(i,j,k)-wav(i,jm,k))     *dyi &
                             +(vav(i,j,k)-vav(i,j,km))     *dzhi(k) ) &
                             )*dyi &
                             + &
                             ( numol * (wav(i,j,kp)-wav(i,j,k)) *dzfi(k) &
                             -numol * (wav(i,j,k)-wav(i,j,km)) *dzfi(km) ) * 2. &
                             * dzhi(k))

               strainav2 =  ( &
                            ((uav(ip,j,k)-uav(i,j,k))    *dxfi(i)     )**2    + &
                            ((vav(i,jp,k)-vav(i,j,k))    *dyi         )**2    + &
                            ((wav(i,j,kp)-wav(i,j,k))    *dzfi(k)     )**2    )

               strainav2 = strainav2 + 0.125 * ( &
                            ((wav(i,j,kp)-wav(im,j,kp))   *dxhi(i)     + &
                            (uav(i,j,kp)-uav(i,j,k))      *dzhi(kp)  )**2    + &
                            ((wav(i,j,k)-wav(im,j,k))     *dxhi(i)     + &
                            (uav(i,j,k)-uav(i,j,km))      *dzhi(k)   )**2    + &
                            ((wav(ip,j,k)-wav(i,j,k))     *dxhi(ip)     + &
                            (uav(ip,j,k)-uav(ip,j,km))    *dzhi(k)   )**2    + &
                            ((wav(ip,j,kp)-wav(i,j,kp))   *dxhi(ip)     + &
                            (uav(ip,j,kp)-uav(ip,j,k))    *dzhi(kp)  )**2    )

               strainav2 = strainav2 + 0.125 * ( &
                            ((uav(i,jp,k)-uav(i,j,k))     *dyi     + &
                            (vav(i,jp,k)-vav(im,jp,k))    *dxhi(i)        )**2    + &
                            ((uav(i,j,k)-uav(i,jm,k))     *dyi     + &
                            (vav(i,j,k)-vav(im,j,k))      *dxhi(i)        )**2    + &
                            ((uav(ip,j,k)-uav(ip,jm,k))   *dyi     + &
                            (vav(ip,j,k)-vav(i,j,k))      *dxhi(ip)       )**2    + &
                            ((uav(ip,jp,k)-uav(ip,j,k))   *dyi     + &
                            (vav(ip,jp,k)-vav(i,jp,k))    *dxhi(ip)       )**2    )

               strainav2 = strainav2 + 0.125 * ( &
                           ((vav(i,j,kp)-vav(i,j,k))    *dzhi(kp) + &
                           (wav(i,j,kp)-wav(i,jm,kp))   *dyi        )**2    + &
                           ((vav(i,j,k)-vav(i,j,km))    *dzhi(k)+ &
                           (wav(i,j,k)-wav(i,jm,k))     *dyi        )**2    + &
                           ((vav(i,jp,k)-vav(i,jp,km))  *dzhi(k)+ &
                           (wav(i,jp,k)-wav(i,j,k))     *dyi        )**2    + &
                           ((vav(i,jp,kp)-vav(i,jp,k))  *dzhi(kp) + &
                           (wav(i,jp,kp)-wav(i,j,kp))   *dyi        )**2    )

               dissresav(i,j,k) = 2.*numol  *(strain2av(i,j,k) - strainav2)  !resolved dissipation

          end do
        end do
      end do

      ! call excjs( tvmy   , ib,ie,jb,je,kb,ke,0,1)   ! jb-1 is not used
      ! call excjs( tsgsmy1, ib,ie,jb,je,kb,ke,0,1)   ! jb-1 is not used
      ! call excjs( tsgsmy2, ib,ie,jb,je,kb,ke,0,1)   ! jb-1 is not used
      ! call excjs( dummyy,  ib,ie,jb,je,kb,ke,0,1)   ! jb-1 is not used
      ! call excjs( ttmy   , ib,ie,jb,je,kb,ke,0,1)   ! jb-1 is not used

      call exchange_halo_z(tvmx, opt_zlevel=(/ih,jh,0/))
      call exchange_halo_z(tsgsmx1, opt_zlevel=(/ih,jh,0/))
      call exchange_halo_z(tsgsmx2, opt_zlevel=(/ih,jh,0/))
      call exchange_halo_z(dummyx, opt_zlevel=(/ih,jh,0/))
      call exchange_halo_z(ttmx, opt_zlevel=(/ih,jh,0/))

      call exchange_halo_z(tvmy, opt_zlevel=(/ih,jh,0/))
      call exchange_halo_z(tsgsmy1, opt_zlevel=(/ih,jh,0/))
      call exchange_halo_z(tsgsmy2, opt_zlevel=(/ih,jh,0/))
      call exchange_halo_z(dummyy, opt_zlevel=(/ih,jh,0/))
      call exchange_halo_z(ttmy, opt_zlevel=(/ih,jh,0/))

      ! BC's
      if (ierank) then
         tvmx   (ie+1,:,:) =  tvmx   (ie,:,:)
         tsgsmx1(ie+1,:,:) =  tsgsmx1(ie,:,:)
         tsgsmx2(ie+1,:,:) =  tsgsmx2(ie,:,:)
         dummyx (ie+1,:,:) =  dummyx (ie,:,:)
         ttmx   (ie+1,:,:) =  ttmx   (ie,:,:)
      end if

      if (jerank) then
         tvmy(:,je+1,:) = tvmy(:,je,:)
         tsgsmy1(:,je+1,:) = tsgsmy1(:,je,:)
         tsgsmy2(:,je+1,:) = tsgsmy2(:,je,:)
         dummyy(:,je+1,:) = dummyy(:,je,:)
         ttmy(:,je+1,:) = ttmy(:,je,:)
      end if

      tvmz   (:,:,ke+1) =  tvmz   (:,:,ke)
      tsgsmz1(:,:,ke+1) =  tsgsmz1(:,:,ke)
      tsgsmz2(:,:,ke+1) =  tsgsmz2(:,:,ke)
      dummyz (:,:,ke+1) =  dummyz(:,:,ke)
      ttmz   (:,:,ke+1) =  ttmz   (:,:,ke)

      do k=kb,ke
        km = k-1
        kp = k+1
        do j=jb,je
          jp = j+1
          jm = j-1
            do i=ib,ie
             im = i-1
             ip = i+1

             ! Total viscous dissipation
             t_vav(i,j,k) =  0.5*( (tvmx(i, j,k) - dummyx(i,j,k) *uav(i, j,k))  + &
                             (tvmx(ip,j,k) - dummyx(ip,j,k)*uav(ip,j,k))) &
                             + 0.5*( (tvmy(i,j, k) - dummyy(i,j,k) *vav(i,j, k))  + &
                             (tvmy(i,jp,k) - dummyy(i,jp,k)*vav(i,jp,k))) &
                             + 0.5*( (tvmz(i,j,k ) - dummyz(i,j,k) *wav(i,j,k ))  + &
                             (tvmz(i,j,kp) - dummyz(i,j,kp)*wav(i,j,kp))) &
                             + dissresav(i,j,k)         ! d/dxj(2*nu*<ui'Sij'>) = <u_i*d/dxj(2*nu*Sij')> +2*nu*<Sij'Sij'>

!      Now the same for subgrid stress
!      <d/dxj(2*u_i'*nu_t*Sij)'> = <u_i'*d/dxj(2*nu_t*Sij)'> + <(2*nu_t*Sij)'*Sij'>
!                                = <u_i*d/dxj(2*nu_t*Sij)> -
!                                  <u_i>*<d/dxj(2*nu_t*Sij)>   + <2*nu_t*Sij*Sij> -
!                                  <(2*nu_t*Sij)>*<Sij>
!                                = <u_i*d/dxj(2*nu_t*Sij)> -
!                                  <u_i>*<d/dxj(2*nu_t*Sij)>   + <2*nu_t*Sij*Sij> -
!                                  2*<nu_t>*<Sij>*<Sij> - 2*<nu_t'*Sij'>*<Sij>


     !---------------------------------------
     !Total subgrid TKE
     !---------------------------------------

             ! Mean SGS dissipation
             disssgsfl(i,j,k) = 2.*nusgsav(i,j,k)*strainav2 ! = 2*<nu_sgs>*<sij>*<sij>


             ! TKE
             tke(i,j,k)       = 0.5*(0.5*(upupav(ip,j,k)+upupav(i,j,k)) + &
                                  0.5*(vpvpav(i,jp,k)+vpvpav(i,j,k)) + &
                                  0.5*(wpwpav(i,j,kp)+wpwpav(i,j,k)))

             ! total SGS
             t_sgsav(i,j,k) =  0.5*( (tsgsmx1(i,j,k) -  uav(i,j,k) *tsgsmx2(i,j,k)) + &
                              (tsgsmx1(ip,j,k) - uav(ip,j,k)*tsgsmx2(ip,j,k))) &
                               + & ! = <2*nu_t*SijSij> - <2*nu_t*Sij>*<Sij>
                              0.5*( (tsgsmy1(i,j,k) -  vav(i,j,k) *tsgsmy2(i,j,k)) + &
                              (tsgsmy1(i,jp,k) - vav(i,jp,k)*tsgsmy2(i,jp,k))) &
                              + & ! = <2*nu_t*SijSij> - <2*nu_t*Sij>*<Sij>
                              0.5*( (tsgsmz1(i,j,k) -  vav(i,j,k) *tsgsmz2(i,j,k)) + &
                              (tsgsmz1(i,j,kp) - vav(i,j,kp)*tsgsmz2(i,j,kp))) &
                              + disssgsav(i,j,k) - disssgsfl(i,j,k)
             ! -2*<nu_t'Sij'>*<Sij>  should still be added!

             ! SGS dissipation
             d_sgsav(i,j,k)= - disssgsav(i,j,k) + disssgsfl(i,j,k)
             ! +2*<nu_t'Sij'>*<Sij>  should still be added! (is compensated with above)


     !---------------------------------------
     !Total pressure TKE
     !---------------------------------------

             ! Pressure correlation term
             ! - <uj'*dp'/dxj> = - <uj*dp/dxj> + <uj>*d<p>/dxj
             t_pav(i,j,k)   = tpm(i,j,k) + &
                              0.5*(uav(i,j,k)*(presav(i,j,k)-presav(i-1,j,k))*dxhi(i) + &
                              uav(i+1,j,k)*(presav(i+1,j,k)-presav(i,j,k))*dxhi(i+1)) &
                              + &
                              0.5*(vav(i,j,k)*(presav(i,j,k)-presav(i,j-1,k))*dyi + &
                              vav(i,j+1,k)*(presav(i,j+1,k)-presav(i,j,k))*dyi) &
                              + &
                              0.5*(wav(i,j,k)*(presav(i,j,k)-presav(i,j,k-1))*dzhi(k) + &
                              wav(i,j,k+1)*(presav(i,j,k+1)-presav(i,j,k))*dzhi(k+1))
                              ! - d/dxj(<0.5*ui'ui'uj'>) = -<uj'd/dxj(<0.5*ui'ui'>) + <ui'uj'><Sij>
!                             = -<uj*d/dxj(0.5*ui'ui')> + <uj>*d/dxj(<0.5*ui'ui'> +
!                             <ui'uj'><Sij>)


!            ttav(i,j,k)   = ttm(i,j,k) -


     !---------------------------------------
     !Total advection TKE
     !---------------------------------------

!            <advection term N.S. times ui> = MKE + A - Pshear - Tt
!            Tt = -<ui'd/dxj(ui'uj')> = -<d/dxj(0.5*ui'ui'uj')> = A + MKE - Pshear - Total

             !Pshear =Ptav = -<ui'uj'>d/dxj(<Sij>) = -<ui'uj'>d<ui>/dxj

             ! mechanical or shear production
             p_tav(i,j,k)    = - ( &
                               0.5 *(upupav(i,j,k)+upupav(ip,j,k))* (uav(ip,j,k)-uav(i,j,k))*dxfi(i)  + & ! <u'u'>*d<u>/dx
                               0.25*(upvpav(i,j,k)  *(uav(i, j, k)-uav(i, jm,k) )*dyi + &
                               upvpav(i,jp,k) *(uav(i, jp,k)-uav(i, j, k) )*dyi + &
                               upvpav(ip,j,k) *(uav(ip,j, k)-uav(ip,jm,k) )*dyi + &
                               upvpav(ip,jp,k)*(uav(ip,jp,k)-uav(ip,j, k) )*dyi) + & ! <u'v'>*d<u>/dy
                               0.25*(upwpav(i, j,k ) *(uav(i, j,k )-uav(i,j,km))*dzhi(k) + &
                               upwpav(i, j,kp) *(uav(i, j,kp)-uav(i, j,k))*dzhi(kp) + &
                               upwpav(ip,j,k ) *(uav(ip,j,k )-uav(ip,j,km))*dzhi(k) + &
                               upwpav(ip,j,kp) *(uav(ip,j,kp)-uav(ip,j,k))*dzhi(kp)) + & ! <u'w'>*d<u>/dz
                               0.25*(upvpav(i, j, k) *(vav(i, j, k)-vav(im,j,k))*dxhi(i) + &
                               upvpav(ip,j, k) *(vav(ip,j, k)-vav(i, j,k))*dxhi(ip) + &
                               upvpav(i, jp,k) *(vav(i,jp,k)-vav(im,jp,k))*dxhi(i) + &
                               upvpav(ip,jp,k) *(vav(ip,jp,k)-vav(i,jp,k))*dxhi(ip)) + & ! <u'v'>*d<v>/dx
                               0.5 *(vpvpav(i,j,k)+vpvpav(i,jp,k))*(vav(i,jp,k)-vav(i,j,k))*dyi + & ! <v'v'>*d<v>/dy
                               0.5 *(vpvpav(i,j,k)+vpvpav(i,jp,k))*(vav(i,jp,k)-vav(i,j,k))*dyi + & ! <v'v'>*d<v>/dy
                               0.25*(vpwpav(i,j ,k ) *(vav(i,j ,k )-vav(i,j,km))*dzhi(k) + &
                               vpwpav(i,j ,kp) *(vav(i,j ,kp)-vav(i,j ,k))*dzhi(kp) + &
                               vpwpav(i,jp,k ) *(vav(i,jp,k)-vav(i,jp,km))*dzhi(k) + &
                               vpwpav(i,jp,kp) *(vav(i,jp,kp)-vav(i,jp,k))*dzhi(kp)) + & ! <v'w'>*d<v>/dz
                               0.25*(upwpav(i, j, k) *(wav(i, j,k )-wav(im,j,k))*dxhi(i) + &
                               upwpav(ip,j, k) *(wav(ip,j,k )-wav(i, j,k))*dxhi(ip) + &
                               upwpav(i, j,kp) *(wav(i,j,kp)-wav(im,j,kp))*dxhi(i) + &
                               upwpav(ip,j,kp) *(wav(ip,j,kp)-wav(i,j,kp))*dxhi(ip)) + & ! <u'w'>*d<w>/dx
                               0.25*(vpwpav(i,j,k)  *(wav(i,j, k )-wav(i,jm,k ) )*dyi + &
                               vpwpav(i,jp,k) *(wav(i,jp,k )-wav(i,j, k ) )*dyi + &
                               vpwpav(ip,j,k) *(wav(i,j, kp)-wav(i,jm,kp) )*dyi + &
                               vpwpav(ip,jp,k)*(wav(i,jp,kp)-wav(i,j, kp) )*dyi) + & ! <v'w'>*d<w>/dy
                               0.5 *(wpwpav(i,j,k)+wpwpav(i,j,kp))*(wav(i,j,kp)-wav(i,j,k))*dzfi(k) ) ! <w'w'>*d<w>/dz

             ! Mean kinetic energy term (expected to be small).
             mke(i,j,k)      = 0.5*(uav(ip,j,k)+uav(i,j,k))*(uuav(ip,j,k)-uuav(i,j,k))*dxfi(i)  +        & !<u>*d<uu>/dx
                               0.5*(uav(i, j,k)*(uvav(i ,jp,k)-uvav(i ,j,k))*dyi  + & ! <u>*d<uv>/dy
                               uav(ip,j,k)*(uvav(ip,jp,k)-uvav(ip,j,k))*dyi) +        &
                               0.5*(uav(i, j,k)*(uwav(i ,j,kp)-uwav(i ,j,k))*dzfi(k)  + & ! <u>*d<uw>/dz
                               uav(ip,j,k)*(uwav(ip,j,kp)-uwav(ip,j,k))*dzfi(k)) +        &
                               0.5*(vav(i,j, k)*(uvav(ip,j ,k)-uvav(i,j ,k))*dxfi(i) + & ! <v>*d<uv>/dx
                               vav(i,jp,k)*(uvav(ip,jp,k)-uvav(i,jp,k))*dxfi(i)) + &
                               0.5*(vav(i,jp,k)+vav(i,j,k))*(vvav(i,jp,k)-vvav(i,j,k))*dyi +        & ! <v>*d<vv>/dy
                               0.5*(vav(i,j ,k)*(vwav(i,j ,kp)-vwav(i,j ,k))*dzfi(k)  + & ! <v>*d<vw>/dz
                               vav(i,jp,k)*(vwav(i,jp,kp)-vwav(i,jp,k))*dzfi(k)) + &
                               0.5*(wav(i,j,k )*(uwav(ip,j,k )-uwav(i,j,k ))*dxfi(i) + & ! <w>*d<uw>/dx
                               wav(i,j,kp)*(uwav(ip,j,kp)-uwav(i,j,kp))*dxfi(i)) +        &
                               0.5*(wav(i,j,k )*(vwav(i,jp,k )-vwav(i,j,k ))*dyi  + & ! <w>*d<vw>/dy
                               wav(i,j,kp)*(vwav(i,jp,kp)-vwav(i,j,kp))*dyi) +        &
                               0.5*(wav(i,j,kp)+wav(i,j,k))*(wwav(i,j,kp)-wwav(i,j,k))*dzfi(k) ! <w>*d<ww>/dz

             ! Advection of TKE
             tkeadv(i,j,k)   = 0.5*(uav(i, j,k)*(tke(i, j,k)-tke(im,j,k))*dxhi(i) + & ! <u>*de/dx
                               uav(ip,j,k)*(tke(ip,j,k)-tke(i ,j,k))*dxhi(ip)) +            & !
                               0.5*(vav(i, j,k)*(tke(i,j ,k)-tke(i,jm,k))*dyi     + & ! <v>*de/dy
                               vav(i,jp,k)*(tke(i,jp,k)-tke(i,j ,k))*dyi) +            &
                               0.5*(wav(i,j,k )*(tke(i,j,k )-tke(i,j,km))*dzhi(k) + & ! <w>*de/dz
                               wav(i,j,kp)*(tke(i,j,kp)-tke(i,j,k ))*dzhi(kp))

             ! <advection term N.S. times ui> = MKE + A - Pshear - Tt
             ! Tt = -<ui'd/dxj(ui'uj')> = -<d/dxj(0.5*ui'ui'uj')> = A      +    MKE   -
             ! Pshear  -   Total
             !                                                    = tkeadv +    mke   -
             !                                                    p_tav   -   ttm
             !        t_tav(i,j,k)   = tkeadv(i,j,k) + mke(i,j,k) - p_tav(i,j,k) - ttm(i,j,k)

             t_tav(i,j,k)   = tkeadv(i,j,k) + mke(i,j,k) - p_tav(i,j,k)  &
                              - 0.5*(ttmx(i,j,k) + ttmx(ip,j,k))        &
                              - 0.5*(ttmy(i,j,k) + ttmy(i,jp,k))        &
                              - 0.5*(ttmz(i,j,k) + ttmz(i,j,kp))

             p_bav(i,j,k)   = (grav/thls)*0.5*(thlpwpav(i,j,k)+thlpwpav(i,j,kp)) !use of thls here...????

          end do
        end do
      end do

    ! need updating tg3315
    call avexy_ibm(p_b(kb:ke+kh),p_bav(:,:,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc,IIcs,.true.)
    call avexy_ibm(t_p(kb:ke+kh),t_pav(:,:,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc,IIcs,.true.)
    call avexy_ibm(adv(kb:ke+kh),tkeadv(:,:,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc,IIcs,.true.)
    call avexy_ibm(t_t(kb:ke+kh),t_tav(:,:,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc,IIcs,.true.)
    call avexy_ibm(t_sgs(kb:ke+kh),t_sgsav(:,:,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc,IIcs,.true.)
    call avexy_ibm(p_t(kb:ke+kh),p_tav(:,:,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc,IIcs,.true.)
    call avexy_ibm(d_sgs(kb:ke+kh),d_sgsav(:,:,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc,IIcs,.true.)
     call avexy_ibm(t_v(kb:ke+kh),t_vav(:,:,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc,IIcs,.true.)

   end subroutine tkestatsdump

  !-------------------------
  !> Clean up when leaving the run
  !------------------------

  subroutine exitstatsdump
    use modstat_nc, only : exitstat_nc
    implicit none
    ! if (ltkedump) then
    !   call exitstat_nc(ncidtke)
    ! endif
  end subroutine exitstatsdump


end module modstatsdump
