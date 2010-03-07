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

    use modglobal,  only : jmax, i1, i2, j1, j2, cp, rlv, zf, nsv, ifnamopt, fname_options
    use modraddata, only : iradiation
    use modfields,  only : thl0, qt0
    use modmpi,     only : myid, nprocs, comm3d, mpierr, my_real, mpi_logical, mpi_integer

    implicit none

    integer   :: i,j,k, ierr
    namelist/NAMSURFACE/ & !< Soil related variables
      isurf,tsoilav, tsoildeepav, phiwav, rootfav, &
      ! Land surface related variables
      lmostlocal, lsmoothflux, lneutral, z0mav, z0hav, rsisurf2, Cskinav, lambdaskinav, albedoav, Qnetav, cvegav, &
      ! Jarvis-Steward related variables
      rsminav, rssoilminav, LAIav, gDav, &
      ! Prescribed values for isurf 2, 3, 4
      z0, thls, ps, ustin, wtsurf, wqsurf, wsvsurf

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
    call MPI_BCAST(tsoilav      , ksoilmax, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(tsoildeepav  , 1       , MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(phiwav       , ksoilmax, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(rootfav      , ksoilmax, MY_REAL, 0, comm3d, mpierr)

    call MPI_BCAST(lmostlocal   , 1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(lsmoothflux  , 1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(lneutral     , 1, MPI_LOGICAL, 0, comm3d, mpierr)
    call MPI_BCAST(z0mav        , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(z0hav        , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(rsisurf2     , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(Cskinav      , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(lambdaskinav , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(albedoav     , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(Qnetav       , 1, MY_REAL, 0, comm3d, mpierr)

    call MPI_BCAST(rsminav      , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(rssoilminav  , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(cvegav       , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(LAIav        , 1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(gDav         , 1, MY_REAL, 0, comm3d, mpierr)

    call MPI_BCAST(z0         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ustin      ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wtsurf     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wqsurf     ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(wsvsurf(1:nsv),nsv,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(ps         ,1,MY_REAL   ,0,comm3d,mpierr)
    call MPI_BCAST(thls       ,1,MY_REAL   ,0,comm3d,mpierr)

    if((z0mav == -1 .and. z0hav == -1) .and. (z0 .ne. -1)) then
      z0mav = z0
      z0hav = z0
      write(6,*) "WARNING: z0m and z0h not defined, set equal to z0"
    end if

    if(isurf == 1) then
      if(tsoilav(1) == -1 .or. tsoilav(2) == -1 .or. tsoilav(3) == -1 .or. tsoilav(4) == -1) then
        stop "NAMSURFACE: tsoil is not set"
      end if
      if(tsoildeepav == -1) then
        stop "NAMSURFACE: tsoildeep is not set"
      end if
      if(phiwav(1) == -1 .or. phiwav(2) == -1 .or. phiwav(3) == -1 .or. phiwav(4) == -1) then
        stop "NAMSURFACE: phiw is not set"
      end if
      if(rootfav(1) == -1 .or. rootfav(2) == -1 .or. rootfav(3) == -1 .or. rootfav(4) == -1) then
        stop "NAMSURFACE: rootf is not set"
      end if
      if(Cskinav == -1) then
        stop "NAMSURFACE: Cskinav is not set"
      end if
      if(lambdaskinav == -1) then
        stop "NAMSURFACE: lambdaskinav is not set"
      end if
      if(albedoav == -1) then
        stop "NAMSURFACE: albedoav is not set"
      end if
      if(Qnetav == -1) then
        stop "NAMSURFACE: Qnetav is not set"
      end if
      if(rsminav == -1) then
        stop "NAMSURFACE: rsminav is not set"
      end if
      if(LAIav == -1) then
        stop "NAMSURFACE: LAIav is not set"
      end if
      if(gDav == -1) then
        stop "NAMSURFACE: gDav is not set"
      end if

    end if

    !if(isurf == 2) then
    !  if(rsisurf2 == -1 .and. (lsea .eqv. .false.)) then
    !    stop "NAMSURFACE: Set rsisurf2 if you use isurf = 2 over land "
    !  end if
    !end if


    if(isurf .ne. 3) then
      if(z0mav == -1) then
        stop "NAMSURFACE: z0mav is not set"
      end if
      if(z0hav == -1) then
        stop "NAMSURFACE: z0hav is not set"
      end if
    end if

    if(isurf == 1) then
    end if

    ! 1.1  -   Allocate arrays
    if(isurf == 1) then
      allocate(zsoil(ksoilmax))
      allocate(dzsoil(ksoilmax))
      allocate(dzsoilh(ksoilmax))

      allocate(lambda(i2,j2,ksoilmax))
      allocate(lambdah(i2,j2,ksoilmax))
      allocate(Dh(i2,j2,ksoilmax))
      allocate(phiw(i2,j2,ksoilmax))
      allocate(pCs(i2,j2,ksoilmax))
      allocate(rootf(i2,j2,ksoilmax))
      allocate(tsoil(i2,j2,ksoilmax))
      allocate(tsoildeep(i2,j2))
      allocate(phitot(i2,j2))

      ! 1.2   -  Initialize arrays
      ! First test, pick ECMWF config
      dzsoil(1) = 0.07
      dzsoil(2) = 0.21
      dzsoil(3) = 0.72
      dzsoil(4) = 1.89

      !! 1.3   -  Calculate vertical layer properties
      zsoil(1)  = dzsoil(1)
      do k = 2, ksoilmax
        zsoil(k) = zsoil(k-1) + dzsoil(k)
      end do

      do k = 1, ksoilmax-1
        dzsoilh(k) = 0.5 * (dzsoil(k+1) + dzsoil(k))
      end do
      dzsoilh(ksoilmax) = 0.5 * dzsoil(ksoilmax)

      ! 1.4   -   Set evaporation related properties
      ! Set water content of soil - constant in this scheme
      phiw(:,:,1) = phiwav(1)
      phiw(:,:,2) = phiwav(2)
      phiw(:,:,3) = phiwav(3)
      phiw(:,:,4) = phiwav(4)

      phitot = 0.0

      do k = 1, ksoilmax
        phitot(:,:) = phitot(:,:) + phiw(:,:,k) * dzsoil(k)
      end do

      phitot(:,:) = phitot(:,:) / zsoil(ksoilmax)

      ! Set root fraction per layer for short grass
      rootf(:,:,1) = rootfav(1)
      rootf(:,:,2) = rootfav(2)
      rootf(:,:,3) = rootfav(3)
      rootf(:,:,4) = rootfav(4)

      ! Calculate conductivity saturated soil
      lambdasat = lambdasm ** (1. - phi) * lambdaw ** (phi)

      tsoil(:,:,1)   = tsoilav(1)
      tsoil(:,:,2)   = tsoilav(2)
      tsoil(:,:,3)   = tsoilav(3)
      tsoil(:,:,4)   = tsoilav(4)
      tsoildeep(:,:) = tsoildeepav

      ! Calculate the soil heat capacity and conductivity based on water content
      ! CvH - If we need prognostic soil moisture at one point, these 2 loops should move to the surface function
      do k = 1, ksoilmax
        do j = 2, j1
          do i = 2, i1
            pCs(i,j,k)    = (1. - phi) * pCm + phiw(i,j,k) * pCw
            Ke            = log10(phiw(i,j,k) / phi) + 1.
            lambda(i,j,k) = Ke * (lambdasat - lambdadry) + lambdadry
          end do
        end do
      end do

      do k = 1, ksoilmax-1
        do j = 2, j1
          do i = 2, i1
            lambdah(i,j,k) = (lambda(i,j,k) * dzsoil(k+1) + lambda(i,j,k+1) * dzsoil(k)) / dzsoilh(k)
          end do
        end do
      end do

      do j = 2, j1
        do i = 2, i1
          lambdah(i,j,ksoilmax) = lambda(i,j,ksoilmax)
        end do
      end do

    end if

    ! 2    -   Initialize land surface
    ! 2.1  -   Allocate arrays
    if(isurf == 1) then
      allocate(Qnet(i2,j2))
      allocate(LE(i2,j2))
      allocate(H(i2,j2))
      allocate(G0(i2,j2))

      Qnet = Qnetav
    end if

    if(isurf <= 2) then
      allocate(rs(i2,j2))
      allocate(rsveg(i2,j2))
      allocate(rsmin(i2,j2))
      allocate(rssoil(i2,j2))
      allocate(rssoilmin(i2,j2))
      allocate(cveg(i2,j2))
      allocate(cliq(i2,j2))
      allocate(ra(i2,j2))
      allocate(tendskin(i2,j2))
      allocate(tskinm(i2,j2))
      allocate(Cskin(i2,j2))
      allocate(lambdaskin(i2,j2))
      allocate(LAI(i2,j2))
      allocate(gD(i2,j2))

      Cskin      = Cskinav
      lambdaskin = lambdaskinav
      rsmin      = rsminav
      rssoilmin  = rsminav
      LAI        = LAIav
      gD         = gDav

      cveg       = cvegav
      cliq       = 0.
    end if

    allocate(albedo(i2,j2))
    allocate(z0m(i2,j2))
    allocate(z0h(i2,j2))
    allocate(obl(i2,j2))
    allocate(tskin(i2,j2))
    allocate(qskin(i2,j2))
    allocate(Cm(i2,j2))
    allocate(Cs(i2,j2))

    if(iradiation == 1) then
      allocate(swdavn(i2,j2,nradtime))
      allocate(swuavn(i2,j2,nradtime))
      allocate(lwdavn(i2,j2,nradtime))
      allocate(lwuavn(i2,j2,nradtime))
    end if

    albedo     = albedoav


    z0m        = z0mav
    z0h        = z0hav

! CvH heterogeneous roughness test
!    do j=1,j2
!      do i=1,i2
!        if( mod(j-2, (jmax*nprocs) / 2) >= 0 .and. mod(j-2, (jmax*nprocs) / 2) < nprocs*jmax / 4 ) then
!          z0m(i,j) = z0m(i,j) * 10.
!          z0h(i,j) = z0h(i,j) * 10.
!        end if
!        !if(i == 2) write(6,*) "CvH", myid, j, z0m(i,j), nprocs, mod(j-2, (jmax*nprocs) / 2), nprocs*jmax / 4
!      end do
!    end do
!


    ! 3. Initialize surface layer
    allocate(ustar   (i2,j2))
    allocate(dudz    (i2,j2))
    allocate(dvdz    (i2,j2))
    allocate(thlflux (i2,j2))
    allocate(qtflux  (i2,j2))
    allocate(dqtdz   (i2,j2))
    allocate(dthldz  (i2,j2))
    allocate(svflux  (i2,j2,nsv))
    allocate(svs(nsv))

    ! CvH set initial values for rs and ra to be able to compute qskin
    if(isurf <= 2) then
      ra = 50.
      if(isurf == 1) then
        rs = 100.
      end if
    end if

    return
  end subroutine initsurface

!> Calculates the interaction with the soil, the surface temperature and humidity, and finally the surface fluxes.
  subroutine surface
    use modglobal,  only : rdt, i1, i2, j1, j2, cp, rlv, fkar, zf, cu, cv, nsv, rk3step, timee, rslabs, pi, pref0, rd, eps1, boltz
    use modraddata, only : iradiation, swu, swd, lwu, lwd, useMcICA
    use modfields,  only : thl0, qt0, u0, v0, rhof, ql0, exnf, presf
    use modmpi,     only : my_real, mpierr, comm3d, mpi_sum, myid, excj
    use moduser,   only : surf_user
    implicit none

    real     :: f1, f2, f3, f4 ! Correction functions for Jarvis-Stewart
    integer  :: i, j, k, n
    real     :: upcu, vpcv, horv, horvav
    real     :: phimzf, phihzf
    real     :: rk3coef, thlsl

    real     :: ust
    real     :: ustl, wtsurfl, wqsurfl

    real     :: swdav, swuav, lwdav, lwuav
    real     :: exner, exnera, tsurfm, Tatm, e, esat, qsat, desatdT, dqsatdT, Acoef, Bcoef
    real     :: fH, fLE, fLEveg, fLEsoil, fLEpot

    if (isurf==10) then
      call surf_user
      return
    end if
    ! 1     -   Calculate the surface temperature
    
    ! CvH start with computation of drag coefficients to allow for implicit solver
    if(isurf <= 2) then

      if(lneutral) then
        obl(:,:) = -1.e10
        oblav    = -1.e10
      else
        call getobl
      end if

      call MPI_BCAST(oblav ,1,MY_REAL ,0,comm3d,mpierr)

      do j = 2, j1
        do i = 2, i1
          if(isurf == 2) then
            rs(i,j) = rsisurf2
          else
            ! 2.1   -   Calculate the surface resistance 
            ! Stomatal opening as a function of incoming short wave radiation
            if (iradiation > 0) then
              f1  = 1. / min(1., (0.004 * max(0.,-swd(i,j,1)) + 0.05) / (0.81 * (0.004 * max(0.,-swd(i,j,1)) + 1.)))
            else
              f1  = 1.
            end if

            ! Soil moisture availability
            f2  = (phifc - phiwp) / (phitot(i,j) - phiwp)

            ! Response of stomata to vapor deficit of atmosphere
            esat = 0.611e3 * exp(17.2694 * (thl0(i,j,1) - 273.16) / (thl0(i,j,1) - 35.86))
            e    = qt0(i,j,1) * ps / 0.622
            f3   = 1. / exp(-gD(i,j) * (esat - e) / 100.)

            ! Response to temperature
            exnera  = (presf(1) / pref0) ** (rd/cp)
            Tatm    = exnera * thl0(i,j,1) + (rlv / cp) * ql0(i,j,1)
            f4      = 1./ (1. - 0.0016 * (298.0 - Tatm) ** 2.)

            rsveg(i,j)  = rsmin(i,j) / LAI(i,j) * f1 * f2 * f3 * f4

            ! 2.2   - Calculate soil resistance based on ECMWF method

            f2  = (phifc - phiwp) / (phiw(i,j,1) - phiwp)
            rssoil(i,j) = rssoilmin(i,j) * f2
            rssoil(i,j) = rssoilmin(i,j) * f2
          end if

          ! 3     -   Calculate the drag coefficient and aerodynamic resistance
          Cm(i,j) = fkar ** 2. / (log(zf(1) / z0m(i,j)) - psim(zf(1) / obl(i,j)) + psim(z0m(i,j) / obl(i,j))) ** 2.
          Cs(i,j) = fkar ** 2. / (log(zf(1) / z0m(i,j)) - psim(zf(1) / obl(i,j)) + psim(z0m(i,j) / obl(i,j))) / (log(zf(1) / z0h(i,j)) - psih(zf(1) / obl(i,j)) + psih(z0h(i,j) / obl(i,j)))

          !upcu  = 0.5 * (u0(i,j,1) + u0(i+1,j,1)) + cu
          !vpcv  = 0.5 * (v0(i,j,1) + v0(i,j+1,1)) + cv
          !horv  = sqrt(upcu ** 2. + vpcv ** 2.)
          !horv  = max(horv, 1.e-2)

          !CvH test smoothz0
          !ra(i,j) = 1. / ( Cs(i,j) * horv )
          horvav  = sqrt(u0av(1) ** 2. + v0av(1) ** 2.)
          ra(i,j) = 1. / ( Cs(i,j) * horvav )

        end do
      end do
    end if

    ! Solve the surface energy balance and the heat and moisture transport in the soil
    if(isurf == 1) then
      thlsl = 0.0
      do j = 2, j1
        do i = 2, i1
          ! 1.1   -   Calculate the heat transport properties of the soil
          ! CvH I put it in the init function, as we don't have prognostic soil moisture at this stage

          ! 1.2   -   Calculate the skin temperature as the top boundary conditions for heat transport
          if(iradiation > 0) then
            if(iradiation == 1 .and. useMcICA) then
              if(rk3step == 1) then
                swdavn(i,j,2:nradtime) = swdavn(i,j,1:nradtime-1)  
                swuavn(i,j,2:nradtime) = swuavn(i,j,1:nradtime-1)  
                lwdavn(i,j,2:nradtime) = lwdavn(i,j,1:nradtime-1)  
                lwuavn(i,j,2:nradtime) = lwuavn(i,j,1:nradtime-1)  

                swdavn(i,j,1) = swd(i,j,1) 
                swuavn(i,j,1) = swu(i,j,1) 
                lwdavn(i,j,1) = lwd(i,j,1) 
                lwuavn(i,j,1) = lwu(i,j,1) 

              end if
              
              swdav = sum(swdavn(i,j,:)) / nradtime
              swuav = sum(swuavn(i,j,:)) / nradtime
              lwdav = sum(lwdavn(i,j,:)) / nradtime
              lwuav = sum(lwuavn(i,j,:)) / nradtime

              Qnet(i,j) = -(swdav + swuav + lwdav + lwuav)
            else
              Qnet(i,j) = -(swd(i,j,1) + swu(i,j,1) + lwd(i,j,1) + lwu(i,j,1))
            end if
          end if

          ! CvH solve the surface temperature implicitly. Do not take into account variations in LWout.
          if(rk3step == 1 .and. timee > 0.) then
            tskinm(i,j) = tskin(i,j)
          end if

          exner   = (ps / pref0) ** (rd/cp)
          tsurfm  = tskinm(i,j) * exner
          
          esat    = 0.611e3 * exp(17.2694 * (tsurfm - 273.16) / (tsurfm - 35.86))
          qsat    = 0.622 * esat / ps
          desatdT = esat * (17.2694 / (tsurfm - 35.86) - 17.2694 * (tsurfm - 273.16) / (tsurfm - 35.86)**2.)
          dqsatdT = 0.622 * desatdT / ps

          ! First, remove LWup from Qnet calculation
          Qnet(i,j) = Qnet(i,j) + boltz * tsurfm ** 4.

          ! Calculate coefficients for surface fluxes
          fH      = rhof(1) * cp / ra(i,j)

          ! Allow for dew fall
          if(qsat - qt0(i,j,1) < 0.) then
            rsveg(i,j)  = 0.
            rssoil(i,j) = 0.
          end if

          fLEveg  = (1. - cliq(i,j)) * cveg(i,j) * rhof(1) * rlv / (ra(i,j) + rsveg(i,j))
          fLEsoil = (1. - cveg(i,j))             * rhof(1) * rlv / (ra(i,j) + rssoil(i,j))
          fLEpot  = cliq(i,j) * cveg(i,j)        * rhof(1) * rlv /  ra(i,j)
          
          fLE     = fLEveg + fLEsoil + fLEpot

          exnera  = (presf(1) / pref0) ** (rd/cp)
          Tatm    = exnera * thl0(i,j,1) + (rlv / cp) * ql0(i,j,1)
          
          ! Calculate coefficients for surface fluxes
          fH      = rhof(1) * cp / ra(i,j)

          ! Allow for dew fall
          if(qsat - qt0(i,j,1) < 0.) then
            rsveg(i,j)  = 0.
            rssoil(i,j) = 0.
          end if

          fLEveg  = (1. - cliq(i,j)) * cveg(i,j) * rhof(1) * rlv / (ra(i,j) + rsveg(i,j))
          fLEsoil = (1. - cveg(i,j))             * rhof(1) * rlv / (ra(i,j) + rssoil(i,j))
          fLEpot  = cliq(i,j) * cveg(i,j)        * rhof(1) * rlv /  ra(i,j)
          
          fLE     = fLEveg + fLEsoil + fLEpot

          exnera  = (presf(1) / pref0) ** (rd/cp)
          Tatm    = exnera * thl0(i,j,1) + (rlv / cp) * ql0(i,j,1)
          
          rk3coef = rdt / (4. - dble(rk3step))
          
          !Acoef   = Qnet(i,j) + fH * Tatm + fLE * (dqsatdT * tsurfm - qsat + qt0(i,j,1)) + lambdaskin(i,j) * tsoil(i,j,1)
          Acoef   = Qnet(i,j) - boltz * tsurfm ** 4. + 4. * boltz * tsurfm ** 4. / rk3coef + fH * Tatm + fLE * (dqsatdT * tsurfm - qsat + qt0(i,j,1)) + lambdaskin(i,j) * tsoil(i,j,1)
          !Bcoef   = fH + fLE * dqsatdT + lambdaskin(i,j)
          Bcoef   = 4. * boltz * tsurfm ** 3. / rk3coef + fH + fLE * dqsatdT + lambdaskin(i,j)
          
          !Acoef   = Qnet(i,j) + fH * Tatm + fLE * (dqsatdT * tsurfm - qsat + qt0(i,j,1)) + lambdaskin(i,j) * tsoil(i,j,1)
          Acoef   = Qnet(i,j) - boltz * tsurfm ** 4. + 4. * boltz * tsurfm ** 4. / rk3coef + fH * Tatm + fLE * (dqsatdT * tsurfm - qsat + qt0(i,j,1)) + lambdaskin(i,j) * tsoil(i,j,1)
          !Bcoef   = fH + fLE * dqsatdT + lambdaskin(i,j)
          Bcoef   = 4. * boltz * tsurfm ** 3. / rk3coef + fH + fLE * dqsatdT + lambdaskin(i,j)

          if (Cskin(i,j) == 0.) then
            tskin(i,j) = Acoef * Bcoef ** (-1.) / exner
          else
            tskin(i,j) = (1. + rk3coef / Cskin(i,j) * Bcoef) ** (-1.) * (tsurfm + rk3coef / Cskin(i,j) * Acoef) / exner
          end if

          Qnet(i,j)     = Qnet(i,j) - (boltz * tsurfm ** 4. + 4. * boltz * tsurfm ** 3. * (tskin(i,j) * exner - tsurfm) / rk3coef)
          !Qnet(i,j)     = Qnet(i,j) - boltz * (tskin(i,j) * exner) ** 4.
          G0(i,j)       = lambdaskin(i,j) * ( tskin(i,j) * exner - tsoil(i,j,1) )
          LE(i,j)       = - fLE * ( qt0(i,j,1) - (dqsatdT * (tskin(i,j) * exner - tsurfm) + qsat))
          if(LE(i,j) == 0.) then
            rs(i,j)     = 1.e8
          else
            rs(i,j)     = - rhof(1) * rlv * (qt0(i,j,1) - (dqsatdT * (tskin(i,j) * exner - tsurfm) + qsat)) / LE(i,j) - ra(i,j)
          end if

          H(i,j)        = - fH  * ( Tatm - tskin(i,j) * exner ) 
          tendskin(i,j) = Cskin(i,j) * (tskin(i,j) - tskinm(i,j)) * exner / rk3coef


          !write(6,*) "SEB: ", Qnet(i,j), H(i,j), LE(i,j), G0(i,j), tendskin(i,j), H(i,j)+LE(i,j)+G0(i,j)+tendskin(i,j)
          !write(6,*) "LEv, LEs ", -fLEveg * (qt0(i,j,1) - (dqsatdT * (tskin(i,j) * exner - tsurfm) + qsat)), -fLEsoil * ( qt0(i,j,1) - (dqsatdT * (tskin(i,j) * exner - tsurfm) + qsat))

          ! 1.4   -   Solve the diffusion equation for the heat transport
          tsoil(i,j,1) = tsoil(i,j,1) + rdt / pCs(i,j,1) * ( lambdah(i,j,ksoilmax) * (tsoil(i,j,2) - tsoil(i,j,1)) / dzsoilh(1) + G0(i,j) ) / dzsoil(1)
          do k = 2, ksoilmax-1
            tsoil(i,j,k) = tsoil(i,j,k) + rdt / pCs(i,j,k) * ( lambdah(i,j,k) * (tsoil(i,j,k+1) - tsoil(i,j,k)) / dzsoilh(k) - lambdah(i,j,k-1) * (tsoil(i,j,k) - tsoil(i,j,k-1)) / dzsoilh(k-1) ) / dzsoil(k)
          end do
          tsoil(i,j,ksoilmax) = tsoil(i,j,ksoilmax) + rdt / pCs(i,j,ksoilmax) * ( lambda(i,j,ksoilmax) * (tsoildeep(i,j) - tsoil(i,j,ksoilmax)) / dzsoil(ksoilmax) - lambdah(i,j,ksoilmax-1) * (tsoil(i,j,ksoilmax) - tsoil(i,j,ksoilmax-1)) / dzsoil(ksoilmax-1) ) / dzsoil(ksoilmax)

          thlsl = thlsl + tskin(i,j)
        end do
      end do

      call MPI_ALLREDUCE(thlsl, thls, 1,  MY_REAL, MPI_SUM, comm3d,mpierr)
      thls = thls / rslabs

      call qtsurf

    elseif(isurf == 2) then
      do j = 2, j1
        do i = 2, i1
          tskin(i,j) = thls
        end do
      end do

      call qtsurf

    end if

    ! 2     -   Calculate the surface fluxes
    if(isurf <= 2) then
      do j = 2, j1
        do i = 2, i1
          upcu   = 0.5 * (u0(i,j,1) + u0(i+1,j,1)) + cu
          vpcv   = 0.5 * (v0(i,j,1) + v0(i,j+1,1)) + cv
          horv   = sqrt(upcu ** 2. + vpcv ** 2.)
          horv   = max(horv, 1.e-2)
          horvav = sqrt(u0av(1) ** 2. + v0av(1) ** 2.)

          ustar  (i,j) = sqrt(Cm(i,j) * horv * horvav)
          thlflux(i,j) = - ( thl0(i,j,1) - tskin(i,j) ) / ra(i,j) 

          !CvH allow for dewfall at night, bypass stomatal resistance
          !if(qt0(i,j,1) - qskin(i,j) > 0.) then
          !  qtflux(i,j) = - (qt0(i,j,1)  - qskin(i,j)) / ra(i,j) 
          !else
          !  qtflux(i,j) = - (qt0(i,j,1)  - qskin(i,j)) / (ra(i,j) + rs(i,j))
          !end if
          qtflux(i,j) = - (qt0(i,j,1)  - qskin(i,j)) / ra(i,j)
          
          do n=1,nsv
            svflux(i,j,n) = wsvsurf(n) 
          enddo

          if (obl(i,j) < 0.) then
            !phimzf = (1.-16.*zf(1)/obl)**(-0.25)
            phimzf = (1. + 3.6 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
            !phihzf = (1.-16.*zf(1)/obl)**(-0.50)
            phihzf = (1. + 7.9 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
          elseif (obl(i,j) > 0.) then
            phimzf = (1.+5.*zf(1)/obl(i,j))
            phihzf = (1.+5.*zf(1)/obl(i,j))
          else
            phimzf = 1.
            phihzf = 1.
          endif

          dudz  (i,j) = ustar(i,j) * phimzf / (fkar*zf(1))*(upcu/horv)
          dvdz  (i,j) = ustar(i,j) * phimzf / (fkar*zf(1))*(vpcv/horv)
          dthldz(i,j) = - thlflux(i,j) / ustar(i,j) * phihzf / (fkar*zf(1))
          dqtdz (i,j) = - qtflux(i,j)  / ustar(i,j) * phihzf / (fkar*zf(1))

        end do
      end do

      if(lsmoothflux) then

        ustl    = sum(ustar  (2:i1,2:j1))
        wtsurfl = sum(thlflux(2:i1,2:j1))
        wqsurfl = sum(qtflux (2:i1,2:j1))

        call MPI_ALLREDUCE(ustl, ust, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)
        call MPI_ALLREDUCE(wtsurfl, wtsurf, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)
        call MPI_ALLREDUCE(wqsurfl, wqsurf, 1,  MY_REAL,MPI_SUM, comm3d,mpierr)

        wtsurf = wtsurf / rslabs
        wqsurf = wqsurf / rslabs

        do j = 2, j1
          do i = 2, i1

            thlflux(i,j) = wtsurf 
            qtflux (i,j) = wqsurf 

            do n=1,nsv
              svflux(i,j,n) = wsvsurf(n)
            enddo

            if (obl(i,j) < 0.) then
              !phimzf = (1.-16.*zf(1)/obl)**(-0.25)
              phimzf = (1. + 3.6 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
              !phihzf = (1.-16.*zf(1)/obl)**(-0.50)
              phihzf = (1. + 7.9 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
            elseif (obl(i,j) > 0.) then
              phimzf = (1.+5.*zf(1)/obl(i,j))
              phihzf = (1.+5.*zf(1)/obl(i,j))
            else
              phimzf = 1.
              phihzf = 1.
            endif

            upcu  = 0.5 * (u0(i,j,1) + u0(i+1,j,1)) + cu
            vpcv  = 0.5 * (v0(i,j,1) + v0(i,j+1,1)) + cv
            horv  = sqrt(upcu ** 2. + vpcv ** 2.)
            horv  = max(horv, 1.e-2)

            dudz  (i,j) = ustar(i,j) * phimzf / (fkar*zf(1))*(upcu/horv)
            dvdz  (i,j) = ustar(i,j) * phimzf / (fkar*zf(1))*(vpcv/horv)
            dthldz(i,j) = - thlflux(i,j) / ustar(i,j) * phihzf / (fkar*zf(1))
            dqtdz (i,j) = - qtflux(i,j)  / ustar(i,j) * phihzf / (fkar*zf(1))
          end do
        end do

      end if

    else

      if(lneutral) then
        obl(:,:) = -1.e10
        oblav    = -1.e10
      else
        call getobl
      end if

      thlsl = 0.
      do j = 2, j1
        do i = 2, i1

          upcu  = 0.5 * (u0(i,j,1) + u0(i+1,j,1)) + cu
          vpcv  = 0.5 * (v0(i,j,1) + v0(i,j+1,1)) + cv
          horv  = sqrt(upcu ** 2. + vpcv ** 2.)
          horv  = max(horv, 1.e-2)

          if( isurf == 4) then
            ustar (i,j) = fkar * horv / (log(zf(1) / z0m(i,j)) - psim(zf(1) / obl(i,j)) + psim(z0m(i,j) / obl(i,j)))
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
            !phimzf = (1.-16.*zf(1)/obl)**(-0.25)
            phimzf = (1. + 3.6 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
            !phihzf = (1.-16.*zf(1)/obl)**(-0.50)
            phihzf = (1. + 7.9 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
          elseif (obl(i,j) > 0.) then
            phimzf = (1.+5.*zf(1)/obl(i,j))
            phihzf = (1.+5.*zf(1)/obl(i,j))
          else
            phimzf = 1.
            phihzf = 1.
          endif

          dudz  (i,j) = ustar(i,j) * phimzf / (fkar*zf(1))*(upcu/horv)
          dvdz  (i,j) = ustar(i,j) * phimzf / (fkar*zf(1))*(vpcv/horv)
          dthldz(i,j) = - thlflux(i,j) / ustar(i,j) * phihzf / (fkar*zf(1))
          dqtdz (i,j) = - qtflux(i,j)  / ustar(i,j) * phihzf / (fkar*zf(1))

          Cs(i,j) = fkar ** 2. / (log(zf(1) / z0m(i,j)) - psim(zf(1) / obl(i,j)) + psim(z0m(i,j) / obl(i,j))) / (log(zf(1) / z0h(i,j)) - psih(zf(1) / obl(i,j)) + psih(z0h(i,j) / obl(i,j)))

          tskin(i,j) = wtsurf / (Cs(i,j) * horv) + thl0(i,j,1)
          thlsl      = thlsl + tskin(i,j)
        end do
      end do

      call MPI_ALLREDUCE(thlsl, thls, 1,  MY_REAL, MPI_SUM, comm3d,mpierr)
      thls = thls / rslabs

      call qtsurf

    end if

    ! Transfer ustar to neighbouring cells
    do j=1,j2
      ustar(1,j)=ustar(i1,j)
      ustar(i2,j)=ustar(2,j)
    end do

    call excj( ustar  , 1, i2, 1, j2, 1,1)

    return

  end subroutine surface

!> Calculate the surface humidity assuming saturation.
  subroutine qtsurf
    use modglobal,   only : tmelt,bt,at,rd,rv,cp,es0,pref0,rslabs,i1,j1
    use modfields,   only : qt0
    use modsurfdata, only : rs, ra
    use modmpi,      only : my_real,mpierr,comm3d,mpi_sum,myid

    implicit none
    real       :: exner, tsurf, qsatsurf, surfwet, es, qtsl
    integer    :: i,j

    if(isurf <= 2) then
      qtsl = 0.
      do j = 2, j1
        do i = 2, i1
          exner      = (ps / pref0)**(rd/cp)
          tsurf      = tskin(i,j) * exner
          es         = es0 * exp(at*(tsurf-tmelt) / (tsurf-bt))
          !qskin(i,j) = rd / rv * es / (ps-(1-rd/rv)*es)
          qsatsurf   = rd / rv * es / (ps-(1-rd/rv)*es)
          surfwet    = ra(i,j) / (ra(i,j) + rs(i,j))
          qskin(i,j) = surfwet * qsatsurf + (1. - surfwet) * qt0(i,j,1)
          qtsl       = qtsl + qskin(i,j)
          !write(6,*), ra(i,j), rs(i,j), surfwet
        end do
      end do

      call MPI_ALLREDUCE(qtsl, qts, 1,  MY_REAL, &
                         MPI_SUM, comm3d,mpierr)
      qts = qts / rslabs
    else
      exner = (ps/pref0)**(rd/cp)
      tsurf = thls*exner
      es    = es0*exp(at*(tsurf-tmelt)/(tsurf-bt))
      qts   = rd/rv*es/(ps-(1-rd/rv)*es)
      ! CvH Check to prevent collapse of spinup u* = 0 and  u = 0 free convection case
      qts   = max(qts, 0.)
    end if

    thvs = thls * (1. + (rv/rd - 1.) * qts)

    return

  end subroutine qtsurf

!> Calculates the Obuhkov length iteratively.
  subroutine getobl
    use modglobal, only : zf, rv, rd, grav, rslabs, i1, j1, i2, j2, timee, cu, cv
    use modfields, only : thl0av, qt0av, u0, v0, thl0, qt0
    use modmpi,    only : my_real,mpierr,comm3d,mpi_sum,myid,excj
    implicit none

    integer             :: i,j,iter
    real                :: thv, L, horv2, horv2l, oblavl
    real                :: Rib, Lstart, Lend, fx, fxdif, Lold

    if(lmostlocal) then

      oblavl = 0.

      do i=1,i2
        do j=1,j2
          thv    = thl0(i,j,1) * (1. + (rv/rd - 1.) * qt0(i,j,1))
          horv2 = u0(i,j,1)*u0(i,j,1) + v0(i,j,1)*v0(i,j,1)
          horv2 = max(horv2, 1.e-2)

          Rib   = grav / thvs * zf(1) * (thv - thvs) / horv2

          iter = 0
          L = obl(i,j)

          if(Rib * L < 0. .or. abs(L) == 1e5) then
            if(Rib > 0) L = 0.01
            if(Rib < 0) L = -0.01
          end if

          do while (.true.)
            iter    = iter + 1
            Lold    = L
            fx      = Rib - zf(1) / L * (log(zf(1) / z0h(i,j)) - psih(zf(1) / L) + psih(z0h(i,j) / L)) / (log(zf(1) / z0m(i,j)) - psim(zf(1) / L) + psim(z0m(i,j) / L)) ** 2.
            Lstart  = L - 0.001*L
            Lend    = L + 0.001*L
            fxdif   = ( (- zf(1) / Lstart * (log(zf(1) / z0h(i,j)) - psih(zf(1) / Lstart) + psih(z0h(i,j) / Lstart)) / (log(zf(1) / z0m(i,j)) - psim(zf(1) / Lstart) + psim(z0m(i,j) / Lstart)) ** 2.) - (-zf(1) / Lend * (log(zf(1) / z0h(i,j)) - psih(zf(1) / Lend) + psih(z0h(i,j) / Lend)) / (log(zf(1) / z0m(i,j)) - psim(zf(1) / Lend) + psim(z0m(i,j) / Lend)) ** 2.) ) / (Lstart - Lend)
            L       = L - fx / fxdif
            if(Rib * L < 0. .or. abs(L) == 1e5) then
              if(Rib > 0) L = 0.01
              if(Rib < 0) L = -0.01
            end if
            if(abs(L - Lold) < 0.0001) exit
          end do

          obl(i,j) = L

          oblavl   = oblavl + obl(i,j)

        end do
      end do

      do i=2,i1
        do j=2,j1
          oblavl = oblavl + obl(i,j)
        end do
      end do

      call MPI_ALLREDUCE(oblavl, oblav, 1,  MY_REAL, MPI_SUM, comm3d,mpierr)
      oblav = oblav / rslabs

    else

      thv    = thl0av(1) * (1. + (rv/rd - 1.) * qt0av(1))
      horv2l = sum( (u0(2:i1,2:j1,1) + cu ) ** 2.)  +  sum( (v0(2:i1,2:j1,1) + cv ) ** 2.)
      horv2l = max(horv2l, 1.e-2)

      call MPI_ALLREDUCE(horv2l, horv2, 1,  MY_REAL, MPI_SUM, comm3d,mpierr)
      horv2 = horv2 / rslabs

      Rib   = grav / thvs * zf(1) * (thv - thvs) / horv2

      iter = 0
      L = oblav

      if(Rib * L < 0. .or. abs(L) == 1e5) then
        if(Rib > 0) L = 0.01
        if(Rib < 0) L = -0.01
      end if

      do while (.true.)
        iter    = iter + 1
        Lold    = L
        fx      = Rib - zf(1) / L * (log(zf(1) / z0hav) - psih(zf(1) / L) + psih(z0hav / L)) / (log(zf(1) / z0mav) - psim(zf(1) / L) + psim(z0mav / L)) ** 2.
        Lstart  = L - 0.001*L
        Lend    = L + 0.001*L
        fxdif   = ( (- zf(1) / Lstart * (log(zf(1) / z0hav) - psih(zf(1) / Lstart) + psih(z0hav / Lstart)) / (log(zf(1) / z0mav) - psim(zf(1) / Lstart) + psim(z0mav / Lstart)) ** 2.) - (-zf(1) / Lend * (log(zf(1) / z0hav) - psih(zf(1) / Lend) + psih(z0hav / Lend)) / (log(zf(1) / z0mav) - psim(zf(1) / Lend) + psim(z0mav / Lend)) ** 2.) ) / (Lstart - Lend)
        L       = L - fx / fxdif
        if(Rib * L < 0. .or. abs(L) == 1e5) then
          if(Rib > 0) L = 0.01
          if(Rib < 0) L = -0.01
        end if
        if(abs(L - Lold) < 0.0001) exit
      end do

      obl   = L
      oblav = L

      !write(6,*) "CvH iter: ", iter, Rib

    end if

    return

  end subroutine getobl

  function psim(zeta)
    implicit none

    real             :: psim
    real, intent(in) :: zeta
    real             :: x

    if(zeta <= 0) then
      !x     = (1. - 16. * zeta) ** (0.25)
      !psim  = 3.14159265 / 2. - 2. * atan(x) + log( (1.+x) ** 2. * (1. + x ** 2.) / 8.)
      ! CvH use Wilson, 2001 rather than Businger-Dyer for correct free convection limit
      x     = (1. + 3.6 * abs(zeta) ** (2./3.)) ** (-0.5)
      psim = 3. * log( (1. + 1. / x) / 2.)
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
      !x     = (1. - 16. * zeta) ** (0.25)
      !psih  = 2. * log( (1. + x ** 2.) / 2. )
      ! CvH use Wilson, 2001
      x     = (1. + 7.9 * abs(zeta) ** (2./3.)) ** (-0.5)
      psih  = 3. * log( (1. + 1. / x) / 2.)
    else
      psih  = -2./3. * (zeta - 5./0.35)*exp(-0.35 * zeta) - (1. + (2./3.) * zeta) ** (1.5) - (10./3.) / 0.35 + 1.
    end if

    return

  end function psih
!>
  subroutine exitsurface
    implicit none
    return
  end subroutine exitsurface

end module modsurface
