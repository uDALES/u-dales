!> \file scalsource.f90

!> Input point, line or planar scalars
!  This file is part of uDALES.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 2006-2021 the uDALES Team.
!

!! ############  The below part is the scalsource script used by Tom Grylls. Now commented by DMajumdar.
!! If wish to use the below, do not forget to also uncomment Namelist reading of variables xS, yS, zS, SS and sigS
!! and broadcasting them in MPI in the modstartup.f90 file.  ############

! subroutine createscals

!   use modglobal,  only : nsv,ib,ie,jb,je,kb,ke,ih,jh,kh,lscasrcr,cexpnr,ifinput,imax,jmax,jtot,&
!                          libm,nblocks,xSa,ySa,zSa,lscasrc,xS,yS,zS
!   use modfields,  only : scar,scarl
!   use modmpi,     only : myid,MY_REAL,comm3d,mpierr
!   use initfac,    only : block
!   implicit none
!   integer :: i,j,k,n,m,p,il,iu,jl,ju,ku

!   ! if (lscasrc) then

!     ! allocate(xSa(1:nsv))
!     ! allocate(ySa(1:nsv))
!     ! allocate(zSa(1:nsv))

!     ! hard code point source locations
!     ! xSa(1:nsv) = (/ 453., 463., 523., 539., 407. /)
!     ! ySa(1:nsv) = (/ 375., 371., 327., 341., 325. /)
!     ! zSa(1:nsv) = (/ 4.   , 4.   , 4. , 4.  , 4.  /)

!     ! all point sources from position defined in namoptions
!     ! xSa = xS; ySa = yS; zSa = zS

!   ! end if

!   if (lscasrcr .AND. nsv.gt.0) then

!     ! read 2-D field of point sources for vehicular emissions network      
!     open (ifinput, file='scals.inp.'//cexpnr)

!     do j=jb,jtot
!       read(ifinput, *) scar(:,j)
!     end do

!     scarl(ib:ie,jb:je) = scar(ib:ie,jb+myid*jmax:je+myid*jmax)

!     ! Set scarl to 0 if set inside an obstacle
!     if (libm) then
!       do n=1,nblocks
!         il = block(n,1)
!         iu = block(n,2)
!         jl = block(n,3)-myid*jmax
!         ju = block(n,4)-myid*jmax
!         ku = block(n,6)
!         if (ju < jb .or. jl > je) then
!           cycle
!         else
!           if (ju > je) ju=je
!           if (jl < jb) jl=jb
!           ! addition due to all blocks at lowest layer
!           if (ku>0) then
!             scarl(il:iu,jl:ju) = 0.
!           end if
!         end if
!       end do
!     end if

!   end if !lscasrcr

! end subroutine

! subroutine scalsource

!   use modglobal,  only : nsv,ib,ie,jb,je,kb,ke,ih,jh,kh,ihc,jhc,khc,xf,zf,xh,zh,dy,jmax,lchem,&
!                          xS,yS,zS,SS,sigS,lscasrc,lscasrcl,lscasrcr,libm,dxfi,dzfi,nblocks,xSa,ySa,zSa
!   use modfields,  only : svp,svpp,scar,scarl
!   use modmpi,     only : myid,mpierr,MY_REAL,comm3d,MPI_SUM
!   use initfac, only : block

!   implicit none
!   integer :: i,j,k,n,il,iu,jl,ju,kl,ku,ncan
!   real :: xL,zL
!   real :: dyi
!   real :: ra2 = 0.
!   real :: scalsum = 0.
!   real :: scalsumt = 0.
!   real :: Pi = 3.1415927

!   dyi = 1./dy

!   ! 2-D network of point sources at lowest level 
!   if (lscasrcr .AND. nsv.gt.0) then
!       if (lchem) then
!         svp(ib:ie,jb:je,kb+1,1) = svp(ib:ie,jb:je,kb+1,1) + 0.522*scarl 
!         svp(ib:ie,jb:je,kb+1,2) = svp(ib:ie,jb:je,kb+1,2) + 0.2*scarl
!         svp(ib:ie,jb:je,kb+1,4) = svp(ib:ie,jb:je,kb+1,4) + scarl
!       else
!         svp(ib:ie,jb:je,kb+1,1) = svp(ib:ie,jb:je,kb+1,1) + scarl
!       end if
!   end if !lscasrcr

!   scalsum = 0.
 
!   !  Input passive scalar point sources
!   if (lscasrc .AND. nsv.gt.0) then 

!     do n=1,nsv 
!     do k=kb,ke
!       do j=jb,je
!         do i=ib,ie
               
!             ra2 = (xf(i)-xS)**2 + ((j+myid*jmax-0.5)*dy-yS)**2 + (zf(k)-zS)**2

!           if (ra2 .LE. 9*sigS**2) then
              
!             scalsum = scalsum + dxfi(i) * dyi * dzfi(k) * SS*exp(-ra2/(2*sigS**2))
!             svp(i,j,k,n) = svp(i,j,k,n) + dxfi(i) * dyi * dzfi(k) * SS*exp(-ra2/(2*sigS**2))

!           end if
!         end do
!       end do
!     end do
!     end do

!     ! Set svpp to 0 when set inside an obstacle
!     if (libm) then
!       do n=1,nblocks
!         il = block(n,1)
!         iu = block(n,2)
!         kl = block(n,5)
!         ku = block(n,6)
!         jl = block(n,3)-myid*jmax
!         ju = block(n,4)-myid*jmax
!         if (ju < jb .or. jl > je) then
!           cycle
!         else
!           if (ju > je) ju=je
!           if (jl < jb) jl=jb
!           svp(il:iu,jl:ju,kl:ku,:) = 0.
!         end if
!       end do
!     end if

!     ! Normalise scalar field to 1/s
!     !call MPI_ALLREDUCE(scalsum,scalsumt,1,MY_REAL,MPI_SUM,comm3d,mpierr)
!     !svpp(:,:,:,1) = svpp(:,:,:,1) / scalsumt
!     !svp(:,:,:,1) = svp(:,:,:,1) + svpp(:,:,:,1)

!     scalsum = 0.

!   end if !lscasrc

!   ! Input passive scalar line sources
!   if (lscasrcl .AND. nsv.gt.0) then

!   ncan = count(block(:,6)>0) !tg3315 update due to block at lowest level

!   zL = zf(kb+1)

!     if (nblocks>0) then

!       do n = 1,ncan-1 ! not ncan and ncan+1 because we do not want release in first and last canyon for BCxs==2
!         if (n == ncan+1) then   ! Added to run for pollutant in first canyon

!           !xL = xh(block(1,1) - (block(2,1) - block(1,2)+1)/2)
!           xL = xh(block(1,1)) - 0.5*(xh(block(2,1)) - xh(block(1,2)+1))
! !          ra2 = (i - (block(1,1) - (block(2,1) - block(1,2))/2.0))**2 + (k)**2 !tg3315 commented for chem validation

!         else !cycle through all other canyons

!           !xL = xh(block(n,2) + (block(2,1) - block(1,2)+1)/2)
!           xL = xh(block(n,2)+1) + 0.5*(xh(block(2,1)) - xh(block(1,2)+1))
! !          ra2 =(xf(i) - xL)**2 + zf(k)**2

!        end if

!         do i=ib,ie
!           do k=kb,ke

! !           if (ra2 .LE. 12*sigS**2) then

!               !scalsum = scalsum + dxf(i) * jmax * dy * dzf(k) * (SS/2*Pi*sigS**2) * exp(-ra2/(2*sigS**2))

!               !tg3315 use this if we want to normalise th scalar conc. !sums values in building too...
! !              scalsum = scalsum + dy * (je - jb +1) * ( (SS/4.) * &
! !                        (erf((xh(i+1)-xL)/(sqrt(2.)*sigS)) - erf((xh(i)-xL)/(sqrt(2.)*sigS))) * &
! !                        (erf((zh(k+1)-zh(kb+1))/(sqrt(2.)*sigS)) - erf((zh(k)-zh(kb+1))/(sqrt(2.)*sigS))) + &
! !                        (SS/4.) * &
! !                        (erf((xh(i+1)-xL)/(sqrt(2.)*sigS)) - erf((xh(i)-xL)/(sqrt(2.)*sigS))) * &
! !                        (erf((zh(k+1)+zh(kb+1))/(sqrt(2.)*sigS)) - erf((zh(k)+zh(kb+1))/(sqrt(2.)*sigS))) )
! !                        * dxfi(i) * dzfi(k)

!               svp(i,jb:je,k,1) = svp(i,jb:je,k,1) + ( (SS/4.) * & ! SS in g/ms... no normalisation
!                         (erf((xh(i+1)-xL)/(sqrt(2.)*sigS)) - erf((xh(i)-xL)/(sqrt(2.)*sigS))) * &
!                         (erf((zh(k+1)-zL)/(sqrt(2.)*sigS)) - erf((zh(k)-zL)/(sqrt(2.)*sigS))) + &
!                         (SS/4.) * & ! reflection from ground...
!                         (erf((xh(i+1)-xL)/(sqrt(2.)*sigS)) - erf((xh(i)-xL)/(sqrt(2.)*sigS))) * &
!                         (erf((zh(k)-2*(zh(k)-zh(kb+1))-zL)/(sqrt(2.)*sigS)) - erf((zh(k+1)-2*(zh(k+1)-zh(kb+1))-zL)/(sqrt(2.)*sigS))) ) &
!                         * dxfi(i) * dzfi(k)

!               !svp(i,jb:je,k,3) = svp(i,jb:je,k,3) + ( (SS/4.) * & ! SS in g/ms... no normalisation
!               !          (erf((xh(i+1)-xL)/(sqrt(2.)*sigS)) - erf((xh(i)-xL)/(sqrt(2.)*sigS))) * &
!               !          (erf((zh(k+1)-zL)/(sqrt(2.)*sigS)) - erf((zh(k)-zL)/(sqrt(2.)*sigS))) + &
!               !          (SS/4.) * & ! reflection from ground...
!               !          (erf((xh(i+1)-xL)/(sqrt(2.)*sigS)) - erf((xh(i)-xL)/(sqrt(2.)*sigS))) * &
!               !          (erf((zh(k)-2*(zh(k)-zh(kb+1))-zL)/(sqrt(2.)*sigS)) - erf((zh(k+1)-2*(zh(k+1)-zh(kb+1))-zL)/(sqrt(2.)*sigS))) ) &
!               !          * dxfi(i) * dzfi(k)

! !            end if
!           end do
!         end do
!       end do
!     end if !nblocks

!     ! Set svpp to 0 when set inside an obstacle
!     if (libm) then
!       do n=1,nblocks
!         il = block(n,1)
!         iu = block(n,2)
!         kl = kb
!         ku = block(n,6)
!         jl = block(n,3)-myid*jmax
!         ju = block(n,4)-myid*jmax
!         if (ju < jb .or. jl > je) then
!           cycle
!         else
!           if (ju > je) ju=je
!           if (jl < jb) jl=jb
!           svp(il:iu,jl:ju,kl:ku,:) = 0.
!         end if
!       end do
!     end if !libm

!     ! Normalise scalar field to 1/s
! !    call MPI_ALLREDUCE(scalsum,scalsumt,1,MY_REAL,MPI_SUM,comm3d,mpierr)

! !    write(*,*), 'scalsum', scalsum

! !    if (lchem) then
!       !svpp(:,:,:,1) = svpp(:,:,:,1)
! !      svp(:,:,:,1) = svp(:,:,:,1) + svpp(:,:,:,1)
!       !svp(:,:,:,2) = svp(:,:,:,2) + 0.1518 * svpp(:,:,:,1)
! !    else
! !      svpp(:,:,:,1) = svpp(:,:,:,1)/ scalsumt !tg3315 not normalised 07/11/2017
! !      svp(:,:,:,1) = svp(:,:,:,1) + svpp(:,:,:,1)
!       !svp(:,:,:,2) = svp(:,:,:,2) + 0.1518 * svpp(:,:,:,1)
! !    end if

!     svpp = 0.
!     scalsum = 0.

!   end if !lscasrcl

! end subroutine scalsource

!! ############  ****  ############


!! ############  Updated scalar source modelling. DMajumdar  ############

module modscalsource
  implicit none

  contains

    subroutine createscals
      use modglobal, only : nsv,cexpnr,ifinput,lscasrc,nscasrc,scasrcp,lscasrcl,nscasrcl,scasrcl, &
                            dx,dy,kb,ke,dzf,zh
      use modmpi,    only : myid,MY_REAL,comm3d,mpierr
      implicit none
      integer :: n,m,k
      real :: line_length2, minimum_sigma, relevant_dz, z_min, z_max
      character(len=2) :: cnsv
      character(80) :: chmess

      if (lscasrc .and. nsv > 0 .and. nscasrc <= 0) then
        if (myid == 0) write(6,*) 'Scalar point sources require nscasrc > 0.'
        stop 1
      end if

      if (lscasrcl .and. nsv > 0 .and. nscasrcl <= 0) then
        if (myid == 0) write(6,*) 'Scalar line sources require nscasrcl > 0.'
        stop 1
      end if

      if (lscasrc .AND. nsv.gt.0 .AND. nscasrc.gt.0) then

        allocate(scasrcp(nscasrc,5,nsv))

        ! read global scalar source locations
        if(myid==0) then
          do m=1,nsv
            write (cnsv, '(i2.1)') m
            open (ifinput,file='scalarsourcep.inp.'//trim(adjustl(cnsv))//'.'//cexpnr)
            read (ifinput,'(a80)') chmess
            read (ifinput,'(a80)') chmess
            do n=1,nscasrc
              read (ifinput,*) &
                    scasrcp(n,1,m), &
                    scasrcp(n,2,m), &
                    scasrcp(n,3,m), &
                    scasrcp(n,4,m), &
                    scasrcp(n,5,m)
            end do
            close (ifinput)
          end do

          ! write (6,*) 'nsv, source_n, xS, yS, zS, SS, sigS '
          ! do m=1,nsv
          !   do n=1,nscasrc
          !     write (6,*) &
          !           m , &
          !           n , &
          !           scasrcp(n,1,m), &
          !           scasrcp(n,2,m), &
          !           scasrcp(n,3,m), &
          !           scasrcp(n,4,m), &
          !           scasrcp(n,5,m)
          !   end do
          ! end do

        end if

        call MPI_BCAST(scasrcp ,5*nscasrc*nsv,MY_REAL ,0,comm3d,mpierr)

        do m = 1,nsv
          do n = 1,nscasrc
            ! On a stretched grid, use the thickness of the cell containing the
            ! source.  At an exact cell boundary, retain the larger neighbouring
            ! thickness.  A source outside the vertical domain uses the nearest
            ! boundary-cell thickness.
            relevant_dz = 0.
            do k = kb,ke
              if (scasrcp(n,3,m) >= zh(k) .and. scasrcp(n,3,m) <= zh(k + 1)) then
                relevant_dz = max(relevant_dz,dzf(k))
              end if
            end do
            if (relevant_dz <= 0.) then
              if (scasrcp(n,3,m) < zh(kb)) then
                relevant_dz = dzf(kb)
              else
                relevant_dz = dzf(ke)
              end if
            end if

            minimum_sigma = max(dx,dy,relevant_dz)
            if (scasrcp(n,5,m) < minimum_sigma) then
              if (myid == 0) write(6,*) 'Under-resolved scalar point source: scalar, source, sig, minimum sig = ', &
                                        m, n, scasrcp(n,5,m), minimum_sigma
              stop 1
            end if
          end do
        end do

      end if

      if (lscasrcl .AND. nsv.gt.0 .AND. nscasrcl.gt.0) then

        allocate(scasrcl(nscasrcl,8,nsv))

        ! read global scalar line source locations
        if(myid==0) then
          do m=1,nsv
            write (cnsv, '(i2.1)') m
            open (ifinput,file='scalarsourcel.inp.'//trim(adjustl(cnsv))//'.'//cexpnr)
            read (ifinput,'(a80)') chmess
            read (ifinput,'(a80)') chmess
            do n=1,nscasrcl
              read (ifinput,*) &
                    scasrcl(n,1,m), &
                    scasrcl(n,2,m), &
                    scasrcl(n,3,m), &
                    scasrcl(n,4,m), &
                    scasrcl(n,5,m), &
                    scasrcl(n,6,m), &
                    scasrcl(n,7,m), &
                    scasrcl(n,8,m)
            end do
            close (ifinput)
          end do

          ! write (6,*) 'nsv, source_n, xSb, ySb, zSb, xSe, ySe, zSe, SS, sigS '
          ! do m=1,nsv
          !   do n=1,nscasrcl
          !     write (6,*) &
          !           m , &
          !           n , &
          !           scasrcl(n,1,m), &
          !           scasrcl(n,2,m), &
          !           scasrcl(n,3,m), &
          !           scasrcl(n,4,m), &
          !           scasrcl(n,5,m), &
          !           scasrcl(n,6,m), &
          !           scasrcl(n,7,m), &
          !           scasrcl(n,8,m)
          !   end do
          ! end do
        end if

        call MPI_BCAST(scasrcl ,8*nscasrcl*nsv,MY_REAL ,0,comm3d,mpierr)

        do m = 1,nsv
          do n = 1,nscasrcl
            line_length2 = (scasrcl(n,4,m) - scasrcl(n,1,m))**2 + &
                          (scasrcl(n,5,m) - scasrcl(n,2,m))**2 + &
                          (scasrcl(n,6,m) - scasrcl(n,3,m))**2
            if (line_length2 <= 0.) then
              if (myid == 0) write(6,*) 'Scalar line-source endpoints must differ: scalar, source = ', m, n
              stop 1
            end if

            ! Use the largest vertical cell thickness crossed by the line.  For a
            ! line wholly outside the vertical domain, use the nearest boundary
            ! cell.  This prevents an unresolved Gaussian source cross-section.
            z_min = min(scasrcl(n,3,m),scasrcl(n,6,m))
            z_max = max(scasrcl(n,3,m),scasrcl(n,6,m))
            relevant_dz = 0.
            do k = kb,ke
              if (z_max >= zh(k) .and. z_min <= zh(k + 1)) then
                relevant_dz = max(relevant_dz,dzf(k))
              end if
            end do
            if (relevant_dz <= 0.) then
              if (z_max < zh(kb)) then
                relevant_dz = dzf(kb)
              else
                relevant_dz = dzf(ke)
              end if
            end if

            minimum_sigma = max(dx,dy,relevant_dz)
            if (scasrcl(n,8,m) < minimum_sigma) then
              if (myid == 0) write(6,*) 'Under-resolved scalar line source: scalar, source, sig, minimum sig = ', &
                                        m, n, scasrcl(n,8,m), minimum_sigma
              stop 1
            end if
          end do
        end do

      end if

      if ((lscasrc .or. lscasrcl) .and. nsv > 0) then
        call initialize_scalar_source_tendency
      end if

    end subroutine


    subroutine initialize_scalar_source_tendency
      use modglobal, only : pi, nsv, ib, ie, jb, je, kb, ke, zf, dx, dy, dxi, dyi, dzfi, &
                            lscasrc, lscasrcl, nscasrc, nscasrcl, scasrcp, scasrcl
      use modfields, only : scalar_source_tendency
      use decomp_2d, only : zstart
      implicit none

      integer :: i, j, k, n, ns, ig, jg
      real :: px, py, pz, vx, vy, vz, lsx, lsy, lsz
      real :: x_start, y_start, z_start, x_end, y_end, z_end
      real :: source_strength, sigma, sigma2, cutoff2, line_length2
      real :: dot_projection, ra2
      real, parameter :: sqrt_2pi = sqrt(2.*pi)

      allocate(scalar_source_tendency(ib:ie,jb:je,kb:ke,nsv))
      scalar_source_tendency = 0.

      ! Precompute passive-scalar point-source tendencies once.  Point sources
      ! are accumulated first to preserve the previous source ordering.
      if (lscasrc .and. nscasrc > 0) then
        do n = 1,nsv
          do ns = 1,nscasrc
            x_start = scasrcp(ns,1,n)
            y_start = scasrcp(ns,2,n)
            z_start = scasrcp(ns,3,n)
            source_strength = scasrcp(ns,4,n)
            sigma2 = scasrcp(ns,5,n)**2
            cutoff2 = 9.*sigma2

            do k = kb,ke
              pz = zf(k)
              do j = jb,je
                jg = zstart(2) + j - jb
                py = (real(jg) - 0.5)*dy
                do i = ib,ie
                  ig = zstart(1) + i - ib
                  px = (real(ig) - 0.5)*dx
                  ra2 = (px - x_start)**2 + (py - y_start)**2 + (pz - z_start)**2

                  if (ra2 <= cutoff2) then
                    scalar_source_tendency(i,j,k,n) = scalar_source_tendency(i,j,k,n) + &
                                                      dxi*dyi*dzfi(k)*source_strength*exp(-ra2/(2.*sigma2))
                  end if
                end do
              end do
            end do
          end do
        end do
      end if

      ! Precompute passive-scalar line-source tendencies once.
      if (lscasrcl .and. nscasrcl > 0) then
        do n = 1,nsv
          do ns = 1,nscasrcl
            x_start = scasrcl(ns,1,n)
            y_start = scasrcl(ns,2,n)
            z_start = scasrcl(ns,3,n)
            x_end = scasrcl(ns,4,n)
            y_end = scasrcl(ns,5,n)
            z_end = scasrcl(ns,6,n)
            source_strength = scasrcl(ns,7,n)
            sigma = scasrcl(ns,8,n)
            sigma2 = sigma**2
            cutoff2 = 9.*sigma2

            lsx = x_end - x_start
            lsy = y_end - y_start
            lsz = z_end - z_start
            line_length2 = lsx*lsx + lsy*lsy + lsz*lsz

            do k = kb,ke
              pz = zf(k)
              do j = jb,je
                jg = zstart(2) + j - jb
                py = (real(jg) - 0.5)*dy
                do i = ib,ie
                  ig = zstart(1) + i - ib
                  px = (real(ig) - 0.5)*dx

                  vx = px - x_start
                  vy = py - y_start
                  vz = pz - z_start
                  dot_projection = (vx*lsx + vy*lsy + vz*lsz)/line_length2

                  if (dot_projection < 0.) then
                    ra2 = (px - x_start)**2 + (py - y_start)**2 + (pz - z_start)**2
                  elseif (dot_projection > 1.) then
                    ra2 = (px - x_end)**2 + (py - y_end)**2 + (pz - z_end)**2
                  else
                    ra2 = (px - (x_start + dot_projection*lsx))**2 + &
                          (py - (y_start + dot_projection*lsy))**2 + &
                          (pz - (z_start + dot_projection*lsz))**2
                  end if

                  if (ra2 <= cutoff2) then
                    scalar_source_tendency(i,j,k,n) = scalar_source_tendency(i,j,k,n) + &
                      dxi*dyi*dzfi(k)*sqrt_2pi*source_strength*sigma* &
                      exp(-ra2/(2.*sigma2))*erf(sqrt((cutoff2 - ra2)/(2.*sigma2)))
                  end if
                end do
              end do
            end do
          end do
        end do
      end if
    end subroutine initialize_scalar_source_tendency


    subroutine scalsource
      use modglobal, only : nsv, ib, ie, jb, je, kb, ke
      use modfields, only : scalar_source_tendency
#if defined(_GPU)
      use modcuda,   only : svp_d, scalar_source_tendency_d
#else
      use modfields, only : svp
#endif
      implicit none

      integer :: i, j, k, n

      if (.not. allocated(scalar_source_tendency)) return

      !$acc parallel loop collapse(4) default(present)
      do n = 1,nsv
        do k = kb,ke
          do j = jb,je
            do i = ib,ie
#if defined(_GPU)
              svp_d(i,j,k,n) = svp_d(i,j,k,n) + scalar_source_tendency_d(i,j,k,n)
#else
              svp(i,j,k,n) = svp(i,j,k,n) + scalar_source_tendency(i,j,k,n)
#endif
            end do
          end do
        end do
      end do
      !$acc end parallel loop
    end subroutine scalsource


    subroutine exitscals
      use modfields, only : scalar_source_tendency
      implicit none
      if (allocated(scalar_source_tendency)) deallocate(scalar_source_tendency)
    end subroutine exitscals

end module modscalsource
