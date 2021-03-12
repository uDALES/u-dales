!> \file modstatistics.f90
!!  Calculates field statistics to be written in modstatsdump.f90
!>
!!  \author Tom Grylls, ICL Dec 16 2016
!
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
module modstatistics

  use modglobal, only : dt,ltkedump
  use modmpi, only : myid
  implicit none
  private
  PUBLIC :: genstats,tkestats
  save

  !NetCDF variables
  integer :: klow,khigh,i,j,k
!  real    :: tsamplep,tstatsdumpp,tsample,tstatsdump

contains

  !-------------------------
  !> Calculate general stats
  !-------------------------

  subroutine genstats(tsamplep,tstatsdumpp,umint,vmint,wmint)

  use modfields,        only : um,up,vm,wm,thlm,uav,vav,wav,uuav,vvav,wwav,uvav,vwav,uwav,thlav,thlwav,thlthlav, &
                               upupav,vpvpav,wpwpav,upvpav,upwpav,vpwpav,thlpwpav
  use modglobal,        only : ib,ie,ih,jb,je,dy,jh,ke,kb,kh,rk3step,timee,cexpnr,tsample,tstatsdump,&
                               ltempeq,dxf,dzf,dzhi
  use modmpi,           only : myid,cmyid,my_real,mpi_sum,mpierr,comm3d
  implicit none

  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)     :: umint
  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)     :: vmint
  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)     :: wmint
  real :: tstatsdumppi,tsamplep,tstatsdumpp
  integer :: km

  tstatsdumppi = 1./tstatsdumpp

!  if (lydump) then
    if (.not. rk3step==3)  return
!      if (tsamplep > tsample) then

        !> Interpolate velocity fields to cell centers
!        do k=kb-kh,ke
!          do j=jb-jh,je+jh
!            do i=ib-ih,ie+ih
!              umint(i,j,k) = 0.5*(um(i,j,k)+um(i+1,j,k))
!              vmint(i,j,k) = 0.5*(vm(i,j,k)+vm(i,j+1,k))
!              wmint(i,j,k) = 0.5*(wm(i,j,k)+wm(i,j,k+1))
!            enddo
!          enddo
!        enddo

        do k=kb,ke
          do j=jb,je
            do i=ib,ie
              uav(i,j,k) = (uav(i,j,k)*(tstatsdumpp-tsamplep) + umint(i,j,k)*tsamplep)*tstatsdumppi
              vav(i,j,k) = (vav(i,j,k)*(tstatsdumpp-tsamplep) + vmint(i,j,k)*tsamplep)*tstatsdumppi
              wav(i,j,k) = (wav(i,j,k)*(tstatsdumpp-tsamplep) + wmint(i,j,k)*tsamplep)*tstatsdumppi
              uuav(i,j,k)  = (uuav(i,j,k)*(tstatsdumpp-tsamplep) + (umint(i,j,k)**2)*tsamplep)*tstatsdumppi
              vvav(i,j,k)  = (vvav(i,j,k)*(tstatsdumpp-tsamplep) + (vmint(i,j,k)**2)*tsamplep)*tstatsdumppi
              wwav(i,j,k)  = (wwav(i,j,k)*(tstatsdumpp-tsamplep) + (wmint(i,j,k)**2)*tsamplep)*tstatsdumppi    
              uvav(i,j,k)  = (uvav(i,j,k)*(tstatsdumpp-tsamplep) + umint(i,j,k)*vmint(i,j,k)*tsamplep)*tstatsdumppi
              vwav(i,j,k)  = (vwav(i,j,k)*(tstatsdumpp-tsamplep) + vmint(i,j,k)*wmint(i,j,k)*tsamplep)*tstatsdumppi
              uwav(i,j,k)  = (uwav(i,j,k)*(tstatsdumpp-tsamplep) + umint(i,j,k)*wmint(i,j,k)*tsamplep)*tstatsdumppi
              if (ltempeq) then
                thlav(i,j,k) = (thlav(i,j,k)*(tstatsdumpp-tsamplep) + thlm(i,j,k)*tsamplep)*tstatsdumppi
                thlwav(i,j,k) = (thlwav(i,j,k)*(tstatsdumpp-tsamplep) + thlm(i,j,k)*wmint(i,j,k)*tsamplep)*tstatsdumppi
                thlthlav(i,j,k) = (thlthlav(i,j,k)*(tstatsdumpp-tsamplep) + (thlm(i,j,k)**2)*tsamplep)*tstatsdumppi
              end if
            end do
          end do
        end do

        upupav = uuav - uav**2             ! overline(u'u') = overline(uu) - U^2
        vpvpav = vvav - vav**2             ! overline(v'v') = overline(vv) - V^2
        wpwpav = wwav - wav**2             ! overline(w'w') = overline(ww) - W^2
        upvpav = uvav - uav*vav            ! overline(u'v') = overline(uv) - U*V
        upwpav = uwav - uav*wav            ! overline(u'w') = overline(uw) - U*W
        vpwpav = vwav - vav*wav            ! overline(v'w') = overline(vw) - V*W

        ! thlw and svw: ib:ie jb:je kb:ke+1  (located on w-faces) !tg3315 BUT thlwav is on cell centre...
        do k=kb,ke+1
          km = k-1
          do j=jb,je
            do i=ib,ie
              thlpwpav(i,j,k) = thlwav(i,j,k) - &
                                0.5 * wav(i,j,k) * & ! no interpolation
                                (thlav(i,j,km)*dzf(k) + thlav(i,j,k)*dzf(km))*dzhi(k) ! interpolate thl to w-faces

!              qlpwpav(i,j,k) = thlwav(i,j,k) - &
!                                0.5 * wav(i,j,k) * & ! no interpolation
!                               (qlav(i,j,km)*dzf(k) + qlav(i,j,k)*dzf(km))*dzhi(k) ! interpolate thl to w-faces

!              qtpwpav(i,j,k) = qtwav(i,j,k) - &
!                                0.5 * wav(i,j,k) * & ! no interpolation
!                                (qtav(i,j,km)*dzf(k) + qtav(i,j,k)*dzf(km))*dzhi(k) ! interpolate thl to w-faces
!
!              do n=1,nsv
!                svpwpav(i,j,k,n) = svwav(i,j,k,n) - &
!                                   0.5 * wav(i,j,k) * & ! no interpolation
!                                   (svav(i,j,km,n)*dzf(k) + svav(i,j,k,n)*dzf(km))*dzhi(k) ! interpolate svav to w-faces
!              end do
            end do
          end do
        end do

        !> generate time averaged stats for TKE budget and call subroutine final field values
        if (ltkedump) then
          call tkestats(tsamplep,tstatsdumpp)
        end if

!        tsample = dt

!     else !timestatsdumpp < tsample

!       tsamplep = tsamplep + dt

!      end if
!    end if
!  end if

  end subroutine genstats

  !-------------------------
  !> Calculate TKE budget terms
  !-------------------------
 
  subroutine tkestats(tsamplep,tstatsdumpp) ! change of variable names not yet translated across to here ! tg3315 30/11/17

  use modfields,        only : u0,v0,w0,thlm,uyt,vyt,wyt,thlyt,pres0,&
                               tvmx,tvmy,tvmz,strain2av,tsgsmx1,tsgsmx2,tsgsmy1,tsgsmy2,&
                               tsgsmz1,tsgsmz2,pres0
  use modglobal,        only : ib,ie,ih,jb,je,jgb,jge,dy,jh,ke,kb,kh,rk3step,cexpnr,tsample,tstatsdump,dzf,zh,dxf,dzf,numol,&
                               dzfi,dxfi,dyi,dy2i,dxfiq,dxhiq,dyiq,dzfi5,dzh,dzf,dzhi,dzhiq,dxf,dxhi
  use modstat_nc,       only : writestat_nc
  use modsurfdata,      only : thls
  use modsubgriddata,   only : ekm
  implicit none

  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)  :: tekm  ! turbulent viscosity 
!  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)  :: emom   
!  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)  :: eomm 
!  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)  :: eopm 
!  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)  :: epom 
!  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)  :: emmo 
!  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)  :: eomp 
!  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)  :: epmo 
!  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)  :: emop 
!  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)  :: empo 
!  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)  :: tkesgs
!  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)  :: nusgs
!  real :: dummy                            
 
  integer i,j,k,im,ip,jm,jp,km,kp
  real tstatsdumppi,tsamplep,tstatsdumpp,strain2,tkesgs,nusgs,&
        emom,eomm,eopm,epom,emmo,eomp,epmo,emop,empo,dummy

  tekm(:,:,:) = ekm(:,:,:) - numol
  
  tstatsdumppi = 1./tstatsdumpp  

  !---------------------------------------
  ! Viscous transport TKE
  !---------------------------------------
 
  !> Time averaged viscous transport in x,y and z to be used to calculate the total viscous transport for TKE
  ! Tvmx at u-locations (ib:ih+ih:jb:je,kb:ke)
  ! This is similar to routine diffu time u_i

  do k = kb,ke
    kp=k+1
    km=k-1
    do j = jb,je
      jp=j+1
      jm=j-1
      do i = ib,ie
        ip=i+1
        im=i-1

        dummy =  u0(i,j,k)*(                           &
                  ( numol  * (u0(i+1,j,k)-u0(i,j,k))*dxfi(i) &
                   -numol * (u0(i,j,k)-u0(i-1,j,k))*dxfi(i-1) ) * 2. * dxhi(i) &
                  + &
                  ( numol * ( (u0(i,jp,k)-u0(i,j,k))   *dyi &
                            +(v0(i,jp,k)-v0(i-1,jp,k))*dxhi(i)) &
                   -numol * ( (u0(i,j,k)-u0(i,jm,k))   *dyi &
                            +(v0(i,j,k)-v0(i-1,j,k))  *dxhi(i)) &
                                        ) * dyi &
                  + &
                  ( numol * ( (u0(i,j,kp)-u0(i,j,k))   *dzhi(kp) &
                            +(w0(i,j,kp)-w0(i-1,j,kp))*dxhi(i)) &
                   -numol * ( (u0(i,j,k)-u0(i,j,km))   *dzhi(k) &
                            +(w0(i,j,k)-w0(i-1,j,k))  *dxhi(i)) & 
                                        ) *dzfi(k)  )
       
         tvmx(i,j,k)=(tvmx(i,j,k)*(tstatsdumpp-tsamplep) + dummy*tsamplep)*tstatsdumppi   ! update average tvmx

         ! Tvmv at v-locations (ib:ih:jb:je+1,kb:ke)
         ! This is similar to routine diffv time v
      
         dummy = v0(i,j,k) * (                            &
                  ( numol * ( (v0(i+1,j,k)-v0(i,j,k))   *dxhi(i+1) &
                           +(u0(i+1,j,k)-u0(i+1,jm,k))*dyi) &
                  -numol * ( (v0(i,j,k)-v0(i-1,j,k))   *dxhi(i) &
                           +(u0(i,j,k)-u0(i,jm,k))    *dyi) &
                                        ) * dxfi(i) &  ! = d/dx( Km*(dv/dx + du/dy) )
                + &
                  (numol * (v0(i,jp,k)-v0(i,j,k)) &
                  -numol * (v0(i,j,k)-v0(i,jm,k))  ) * 2. * dy2i &        ! = d/dy( 2*Km*(dv/dy) )
                + &
                  (numol * ( (v0(i,j,kp)-v0(i,j,k))    *dzhi(kp) &
                           +(w0(i,j,kp)-w0(i,jm,kp))  *dyi) &
                  -numol * ( (v0(i,j,k)-v0(i,j,km))    *dzhi(k) &
                           +(w0(i,j,k)-w0(i,jm,k))    *dyi) & 
                                        ) * dzfi(k) )      ! = d/dz( Km*(dv/dz + dw/dy) )
   
         tvmy(i,j,k)=(tvmy(i,j,k)*(tstatsdumpp-tsamplep) + dummy*tsamplep)*tstatsdumppi   ! update average tvmy

         ! Tvmz at w-locations (ib:ih:jb:je,kb:ke+kh)
         ! This is similar to routine diffw time w
        
         dummy = w0(i,j,k) * (                                         &
                   ( numol * ( (w0(i+1,j,k)-w0(i,j,k))    *dxhi(i+1) &
                             +(u0(i+1,j,k)-u0(i+1,j,km)) *dzhi(k) ) &
                    -numol * ( (w0(i,j,k)-w0(i-1,j,k))    *dxhi(i) &
                             +(u0(i,j,k)-u0(i,j,km))     *dzhi(k) ) &
                                       )*dxfi(i) &
                + &
                  ( numol * ( (w0(i,jp,k)-w0(i,j,k))     *dyi &
                             +(v0(i,jp,k)-v0(i,jp,km))   *dzhi(k) ) &
                   -numol * ( (w0(i,j,k)-w0(i,jm,k))     *dyi &
                             +(v0(i,j,k)-v0(i,j,km))     *dzhi(k) ) &
                                       )*dyi &
                + &
                  ( numol * (w0(i,j,kp)-w0(i,j,k)) *dzfi(k) &
                   -numol * (w0(i,j,k)-w0(i,j,km)) *dzfi(km) ) * 2. &
                                        * dzhi(k) )

         tvmz(i,j,k)=(tvmz(i,j,k)*(tstatsdumpp-tsamplep) + dummy*tsamplep)*tstatsdumppi   ! update average uwsgsav


         ! Compute stresses and fluxes at c.c. , also used in total viscous transport
         strain2 =  ( &
            ((u0(ip,j,k)-u0(i,j,k))    *dxfi(i)     )**2    + &
            ((v0(i,jp,k)-v0(i,j,k))    *dyi         )**2    + &
            ((w0(i,j,kp)-w0(i,j,k))    *dzfi(k)     )**2    )

          strain2 = strain2 + 0.125 * ( &
            ((w0(i,j,kp)-w0(im,j,kp))   *dxhi(i)     + &
            (u0(i,j,kp)-u0(i,j,k))      *dzhi(kp)  )**2    + &
            ((w0(i,j,k)-w0(im,j,k))     *dxhi(i)     + &
            (u0(i,j,k)-u0(i,j,km))      *dzhi(k)   )**2    + &
            ((w0(ip,j,k)-w0(i,j,k))     *dxhi(ip)     + &
            (u0(ip,j,k)-u0(ip,j,km))    *dzhi(k)   )**2    + &
            ((w0(ip,j,kp)-w0(i,j,kp))   *dxhi(ip)     + &
            (u0(ip,j,kp)-u0(ip,j,k))    *dzhi(kp)  )**2    )

          strain2 = strain2 + 0.125 * ( &
            ((u0(i,jp,k)-u0(i,j,k))     *dyi     + &
            (v0(i,jp,k)-v0(im,jp,k))    *dxhi(i)        )**2    + &
            ((u0(i,j,k)-u0(i,jm,k))     *dyi     + &
            (v0(i,j,k)-v0(im,j,k))      *dxhi(i)        )**2    + &
            ((u0(ip,j,k)-u0(ip,jm,k))   *dyi     + &
            (v0(ip,j,k)-v0(i,j,k))      *dxhi(ip)       )**2    + &
            ((u0(ip,jp,k)-u0(ip,j,k))   *dyi     + &
            (v0(ip,jp,k)-v0(i,jp,k))    *dxhi(ip)       )**2    )

          strain2 = strain2 + 0.125 * ( &
            ((v0(i,j,kp)-v0(i,j,k))    *dzhi(kp) + &
            (w0(i,j,kp)-w0(i,jm,kp))   *dyi        )**2    + &
            ((v0(i,j,k)-v0(i,j,km))    *dzhi(k)+ &
            (w0(i,j,k)-w0(i,jm,k))     *dyi        )**2    + &
            ((v0(i,jp,k)-v0(i,jp,km))  *dzhi(k)+ &
            (w0(i,jp,k)-w0(i,j,k))     *dyi        )**2    + &
            ((v0(i,jp,kp)-v0(i,jp,k))  *dzhi(kp) + &
            (w0(i,jp,kp)-w0(i,j,kp))   *dyi        )**2    )

          strain2av(i,j,k)=(strain2av(i,j,k)*(tstatsdumpp-tsamplep) + strain2*tsamplep)*tstatsdumppi   ! update average strain2av



        !--------------------------------------------------
        !> SGS TKE
        !--------------------------------------------------

        ! x-direction
        emom = ( dzf(km) * ( tekm(i,j,k)*dxf(i-1)  + tekm(i-1,j,k)*dxf(i) )  + &             ! dx is non-equidistant
                 dzf(k)  * ( tekm(i,j,km)*dxf(i-1) + tekm(i-1,j,km)*dxf(i) ) )*dxhi(i) * dzhiq(k)
        emop = ( dzf(kp) * ( tekm(i,j,k)*dxf(i-1)  + tekm(i-1,j,k)*dxf(i) )  + &              ! dx is non-equidistant
                 dzf(k)  * ( tekm(i,j,kp)*dxf(i-1) + tekm(i-1,j,kp)*dxf(i) ) )*dxhi(i) * dzhiq(kp)
        empo = 0.25 * ( ( tekm(i,j,k)+tekm(i,jp,k))*dxf(i-1) + (tekm(i-1,j,k)+tekm(i-1,jp,k))*dxf(i) ) * dxhi(i)  ! dx is non-equidistant
        emmo = 0.25 * ( ( tekm(i,j,k)+tekm(i,jm,k))*dxf(i-1)  +(tekm(i-1,jm,k)+tekm(i-1,j,k))*dxf(i) ) * dxhi(i)  ! dx is non-equidistant

!        dummy =  u0(i,j,k)*(                           &
        dummy =  (                           &
                   ( tekm(i,j,k)  * (u0(i+1,j,k)-u0(i,j,k))*dxfi(i) &
                    -tekm(i-1,j,k)* (u0(i,j,k)-u0(i-1,j,k))*dxfi(i-1) ) * 2. * dxhi(i) &
                  + &
                  ( empo * ( (u0(i,jp,k)-u0(i,j,k))   *dyi &
                            +(v0(i,jp,k)-v0(i-1,jp,k))*dxhi(i)) &
                    -emmo * ( (u0(i,j,k)-u0(i,jm,k))   *dyi &
                            +(v0(i,j,k)-v0(i-1,j,k))  *dxhi(i)) &
                                       ) * dyi &
                  + &
                  ( emop * ( (u0(i,j,kp)-u0(i,j,k))   *dzhi(kp) &
                            +(w0(i,j,kp)-w0(i-1,j,kp))*dxhi(i)) &
                    -emom * ( (u0(i,j,k)-u0(i,j,km))   *dzhi(k) &
                            +(w0(i,j,k)-w0(i-1,j,k))  *dxhi(i)) &
                                       ) *dzfi(k) )

        tsgsmx1(i,j,k)=(tsgsmx1(i,j,k)*(tstatsdumpp-tsamplep) + dummy*u0(i,j,k)*tsamplep)*tstatsdumppi   ! update average tsgsmx1
        tsgsmx2(i,j,k)=(tsgsmx2(i,j,k)*(tstatsdumpp-tsamplep) + dummy*tsamplep)*tstatsdumppi             ! update average tsgsmx2
        
        ! y-direction
          eomm = ( dzf(km) * ( tekm(i,j,k)  + tekm(i,jm,k)  )  + &
                      dzf(k)  * ( tekm(i,j,km) + tekm(i,jm,km) ) ) * dzhiq(k)
          eomp = ( dzf(kp) * ( tekm(i,j,k)  + tekm(i,jm,k)  )  + &
                      dzf(k)  * ( tekm(i,j,kp) + tekm(i,jm,kp) ) ) * dzhiq(kp)
          emmo = 0.25 * ( ( tekm(i,j,k)+tekm(i,jm,k))*dxf(i-1)  +(tekm(i-1,jm,k)+tekm(i-1,j,k))*dxf(i) ) * dxhi(i)  ! dx is non-equidistant
          epmo = 0.25 * ( ( tekm(i,j,k)+tekm(i,jm,k))*dxf(i+1) + (tekm(i+1,jm,k)+tekm(i+1,j,k))*dxf(i) ) * dxhi(i+1)  ! dx is non-equidistant

!       dummy = v0(i,j,k) * (                            &
       dummy = (                            &
               ( epmo * ( (v0(i+1,j,k)-v0(i,j,k))   *dxhi(i+1) &
                        +(u0(i+1,j,k)-u0(i+1,jm,k))*dyi) &
                -emmo * ( (v0(i,j,k)-v0(i-1,j,k))   *dxhi(i) &
                        +(u0(i,j,k)-u0(i,jm,k))    *dyi) &
                           ) * dxfi(i) &        ! = d/dx( Km*(dv/dx + du/dy) )
                + &
              (tekm(i,j,k) * (v0(i,jp,k)-v0(i,j,k)) &
              -tekm(i,jm,k)* (v0(i,j,k)-v0(i,jm,k))  ) * 2. * dy2i &        ! = d/dy( 2*Km*(dv/dy) )
                + &
              ( eomp * ( (v0(i,j,kp)-v0(i,j,k))    *dzhi(kp) &
                        +(w0(i,j,kp)-w0(i,jm,kp))  *dyi) &
                -eomm * ( (v0(i,j,k)-v0(i,j,km))    *dzhi(k) &
                        +(w0(i,j,k)-w0(i,jm,k))    *dyi)   &
                           ) * dzfi(k)  )     ! = d/dz( Km*(dv/dz + dw/dy) )

        tsgsmy1(i,j,k)=(tsgsmy1(i,j,k)*(tstatsdumpp-tsamplep) + dummy*v0(i,j,k)*tsamplep)*tstatsdumppi   ! update average tsgsmy1  = <v*d/dxj(2*nu*S2j)>
        tsgsmy2(i,j,k)=(tsgsmy2(i,j,k)*(tstatsdumpp-tsamplep) + dummy*tsamplep)*tstatsdumppi             ! update average tsgsmy2  = <d/dxj(2*nu*S2j)>

        ! z-direction
          emom = ( dzf(km) * ( tekm(i,j,k)*dxf(i-1)  + tekm(i-1,j,k)*dxf(i) )*dxhi(i)  + &
                      dzf(k)  * ( tekm(i,j,km)*dxf(i-1) + tekm(i-1,j,km)*dxf(i) )*dxhi(i) ) * dzhiq(k)
          eomm = ( dzf(km) * ( tekm(i,j,k)  + tekm(i,jm,k)  )  + &
                      dzf(k)  * ( tekm(i,j,km) + tekm(i,jm,km) ) ) * dzhiq(k)
          eopm = ( dzf(km) * ( tekm(i,j,k)  + tekm(i,jp,k)  )  + &
                      dzf(k)  * ( tekm(i,j,km) + tekm(i,jp,km) ) ) * dzhiq(k)
          epom = ( dzf(km) * ( tekm(i,j,k)*dxf(i+1)  + tekm(i+1,j,k)*dxf(i) )*dxhi(i+1)  + &
                      dzf(k)  * ( tekm(i,j,km)*dxf(i+1) + tekm(i+1,j,km)*dxf(i) )*dxhi(i+1) ) * dzhiq(k)

!        dummy = w0(i,j,k) * (                                         &
        dummy =   (                                         &
                  ( epom * ( (w0(i+1,j,k)-w0(i,j,k))    *dxhi(i+1) &
                            +(u0(i+1,j,k)-u0(i+1,j,km)) *dzhi(k) ) &
                    -emom * ( (w0(i,j,k)-w0(i-1,j,k))    *dxhi(i) &
                            +(u0(i,j,k)-u0(i,j,km))     *dzhi(k) ) &
                             )*dxfi(i) &
                + &
                  ( eopm * ( (w0(i,jp,k)-w0(i,j,k))     *dyi &
                            +(v0(i,jp,k)-v0(i,jp,km))   *dzhi(k) ) &
                    -eomm * ( (w0(i,j,k)-w0(i,jm,k))     *dyi &
                            +(v0(i,j,k)-v0(i,j,km))     *dzhi(k) ) &
                             )*dyi &
                + &
                  ( tekm(i,j,k) * (w0(i,j,kp)-w0(i,j,k)) *dzfi(k) &
                  -tekm(i,j,km)* (w0(i,j,k)-w0(i,j,km)) *dzfi(km) ) * 2. &
                                                              * dzhi(k))

        tsgsmz1(i,j,k)=(tsgsmz1(i,j,k)*(tstatsdumpp-tsamplep) + dummy*w0(i,j,k)*tsamplep)*tstatsdumppi   ! update average tsgsmz1 = <w*d/dxj(2*nu*S3j)>
        tsgsmz2(i,j,k)=(tsgsmz2(i,j,k)*(tstatsdumpp-tsamplep) + dummy*tsamplep)*tstatsdumppi            ! update average tsgsmz2 = <d/dxj(2*nu*S3j)>

        end do
      end do
    end do

  end subroutine tkestats

end module modstatistics
