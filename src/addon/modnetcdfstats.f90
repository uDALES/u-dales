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
module modnetcdfstats

use modglobal

implicit none
private
public :: initnetcdfstats, netcdfstats, exitnetcdfstats

save
  !namelist variables
  real                 :: dtav,tnext
  logical              :: lnetcdf  = .false.
  integer              :: ncklimit = 60

  !VARIABLES FOR STATISTICS
  !id of netcdf file
  integer :: ncid

  !id of dimensions
  integer :: xid, yid, zid, tid

  !id of variables (means)
  integer :: uavgid, vavgid, wavgid, thlavgid, qtavgid, eavgid

  !id of variables (variances)
  integer :: uvarid, vvarid, wvarid, thlvarid, qtvarid

  !id of variables (covariances)
  integer :: uvcovid, uwcovid, vwcovid
  integer :: vwcovsid
  integer :: uthlcovid, vthlcovid, wthlcovid
  integer :: wthlcovsid
  integer :: uqtcovid, vqtcovid, wqtcovid
  integer :: wqtcovsid
  integer :: thlqcovid

  !VARIABLES FOR MOVIE
  !id of netcdf
  integer :: ncidmovie

  !id of dimensions
  integer :: xidmovie, yidmovie, zidmovie, tidmovie

  !id of variables
  integer :: uidmovie, vidmovie, widmovie, thlidmovie, qtidmovie


  !VARIABLES FOR FIELD DUMP
  !id of netcdf
  integer :: ncidfieldu, ncidfieldv, ncidfieldw, ncidfieldthl, ncidfieldqt

  !id of dimensions
  integer :: xidfieldu,   yidfieldu,   zidfieldu,   tidfieldu
  integer :: xidfieldv,   yidfieldv,   zidfieldv,   tidfieldv
  integer :: xidfieldw,   yidfieldw,   zidfieldw,   tidfieldw
  integer :: xidfieldthl, yidfieldthl, zidfieldthl, tidfieldthl
  integer :: xidfieldqt,  yidfieldqt,  zidfieldqt,  tidfieldqt

  !id of variables
  integer :: uidfield, vidfield, widfield, thlidfield, qtidfield

  !COMMON VARIABLES
  !constants
  !integer ncklimit
  !parameter (ncklimit = 50)
  integer :: nccall = 1

  integer :: ncfieldflag
  parameter (ncfieldflag = 31)
  ! 1 = u
  ! 2 = v
  ! 4 = w
  ! 8 = thl
  ! 16 = qt

contains
  subroutine initnetcdfstats

    use typeSizes
    use netcdf
    use modmpi
    use modglobal

    implicit none

    integer             :: status
    character(len = 20) :: ncfile
    integer             :: n, ierr
    character(20)       :: name

    namelist/NAMNETCDFSTATS/ &
    dtav,lnetcdf,ncklimit

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMNETCDFSTATS,iostat=ierr)
      write(6, NAMNETCDFSTATS)
      close(ifnamopt)
    end if

    call MPI_BCAST(dtav       ,1,MY_REAL,     0,comm3d,mpierr)
    call MPI_BCAST(lnetcdf    ,1,MPI_LOGICAL, 0,comm3d,mpierr)
    call MPI_BCAST(ncklimit   ,1,MPI_INTEGER, 0,comm3d,mpierr)

    if(.not.(lnetcdf)) return
    tnext = dtav-1e-6
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'NETCDFSTATS: dtav should be a integer multiple of dtmax'
    end if

    ncfile = 'data/stats123.nc'
    write(ncfile(11:13),'(i3.3)') myid
    write(6,*) "NETCDFSTATS: Creating: ", ncfile

    !create file
    status = nf90_create(ncfile, nf90_clobber, ncid)
    if (status /= nf90_noerr) call nchandle_error(status)

    !create dimensions
    status = nf90_def_dim(ncid, "y", jmax, yid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_dim(ncid, "z", ncklimit, zid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_dim(ncid, "t", nf90_unlimited, tid)
    if (status /= nf90_noerr) call nchandle_error(status)

    !create variables
    !means
    status = nf90_def_var(ncid, "uavg", nf90_float, (/yid, zid, tid/), uavgid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "vavg", nf90_float, (/yid, zid, tid/), vavgid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "wavg", nf90_float, (/yid, zid, tid/), wavgid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "thlavg", nf90_float, (/yid, zid, tid/), thlavgid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "qtavg", nf90_float, (/yid, zid, tid/), qtavgid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "eavg", nf90_float, (/yid, zid, tid/), eavgid)
    if (status /= nf90_noerr) call nchandle_error(status)

    !variances
    status = nf90_def_var(ncid, "uvar", nf90_float, (/yid, zid, tid/), uvarid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "vvar", nf90_float, (/yid, zid, tid/), vvarid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "wvar", nf90_float, (/yid, zid, tid/), wvarid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "thlvar", nf90_float, (/yid, zid, tid/), thlvarid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "qtvar", nf90_float, (/yid, zid, tid/), qtvarid)
    if (status /= nf90_noerr) call nchandle_error(status)

    !covariances
    status = nf90_def_var(ncid, "vwcov", nf90_float, (/yid, zid, tid/), vwcovid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "vwcovs", nf90_float, (/yid, zid, tid/), vwcovsid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "wthlcov", nf90_float, (/yid, zid, tid/), wthlcovid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "wthlcovs", nf90_float, (/yid, zid, tid/), wthlcovsid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "wqtcov", nf90_float, (/yid, zid, tid/), wqtcovid)
    if (status /= nf90_noerr) call nchandle_error(status)
    status = nf90_def_var(ncid, "wqtcovs", nf90_float, (/yid, zid, tid/), wqtcovsid)
    if (status /= nf90_noerr) call nchandle_error(status)

    !turn off define mode
    status = nf90_enddef(ncid)
    if (status /= nf90_noerr) call nchandle_error(status)
  end subroutine initnetcdfstats


  subroutine netcdfstats

    use modglobal, only : rk3step,ntimee
    implicit none

    if (.not. lnetcdf) return
    if (rk3step/=3) return
    if (rk3step/=3) return
    if (timee==0) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+dtav
    dt_lim = min(dt_lim,tnext)

    call do_netcdfstats
    nccall = nccall + 1

  end subroutine netcdfstats


  subroutine do_netcdfstats

    use typeSizes
    use netcdf
    use modfields
    use modmpi, only : myid
    use modsurfdata
    use modsubgrid, only : ekm, ekh

    implicit none

    integer n,i,j,k
    integer status

    real, dimension(jmax,ncklimit) :: uavg, vavg, wavg, thlavg, qtavg, eavg, thlhavg, qthavg, vonwavg
    real, dimension(jmax,ncklimit) :: uvar, vvar, wvar, thlvar, qtvar
    real, dimension(jmax,ncklimit) :: vwcov, vwcovs
    real, dimension(jmax,ncklimit) :: wthlcov, wthlcovs, wqtcov, wqtcovs

    real  vonw(2-ih:i1+ih,2-jh:j1+jh,k1),putout(2-ih:i1+ih,2-jh:j1+jh,k1)

    uavg(:,:) = 0.0
    vavg(:,:) = 0.0
    wavg(:,:) = 0.0
    thlavg(:,:) = 0.0
    eavg(:,:) = 0.0
    qtavg(:,:) = 0.0

    thlhavg(:,:) = 0.0
    qthavg(:,:) = 0.0
    vonwavg(:,:) = 0.0

    uvar(:,:) = 0.0
    vvar(:,:) = 0.0
    wvar(:,:) = 0.0
    thlvar(:,:) = 0.0
    qtvar(:,:) = 0.0

    vwcov(:,:) = 0.0
    vwcovs(:,:) = 0.0
    wthlcov(:,:) = 0.0
    wthlcovs(:,:) = 0.0
    wqtcov(:,:) = 0.0
    wqtcovs(:,:) = 0.0

    !calculate averages and store them

    !Prepare data
    do k = 2,k1
      do j = 2,j1
        do i = 2,i1
          thl0h(i,j,k) = 0.5*(thl0(i,j,k) + thl0(i,j,k-1))
          qt0h(i,j,k) = 0.5*(qt0(i,j,k) + qt0(i,j,k-1))
          vonw(i,j,k) = 0.25*(v0(i,j,k) + v0(i,j+1,k) + v0(i,j,k-1) + v0(i,j+1,k-1))
        end do
      end do
    end do


    !LOOPS ARE NOT PUT IN FUNCTION BECAUSE OF ARRAY DEFINITIONS WHICH DIFFER AMONG VARIABLES!
    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          uavg(j,k) = uavg(j,k) + u0(i+1,j+1,k)
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          vavg(j,k) = vavg(j,k) + v0(i+1,j+1,k)
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          wavg(j,k) = wavg(j,k) + w0(i+1,j+1,k)
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          thlavg(j,k) = thlavg(j,k) + thl0(i+1,j+1,k)
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          qtavg(j,k) = qtavg(j,k) + qt0(i+1,j+1,k)
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          eavg(j,k) = eavg(j,k) + (e120(i+1,j+1,k))**2.0
        end do
      end do
    end do

    !create average for thl0h - only used for cov
    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          thlhavg(j,k) = thlhavg(j,k) + thl0h(i+1,j+1,k)
        end do
      end do
    end do

    !create average for qt0h - only used for cov
    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          qthavg(j,k) = qthavg(j,k) + qt0h(i+1,j+1,k)
        end do
      end do
    end do

    !create average for v projected on w - only used for cov, start at level 2
    do k = 2,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          vonwavg(j,k) = vonwavg(j,k) + vonw(i+1,j+1,k)
        end do
      end do
    end do

    uavg   = uavg / imax
    vavg   = vavg / imax
    wavg   = wavg / imax
    thlavg = thlavg / imax
    qtavg  = qtavg / imax
    eavg   = eavg / imax

    thlhavg = thlhavg / imax
    qthavg  = qthavg / imax
    vonwavg = vonwavg / imax

    status = nf90_put_var(ncid, uavgid, uavg, (/1,1,nccall/), (/jmax, ncklimit , 1/))
    if(status /= nf90_noerr) call nchandle_error(status)
    status = nf90_put_var(ncid, vavgid, vavg, (/1,1,nccall/), (/jmax, ncklimit , 1/))
    if(status /= nf90_noerr) call nchandle_error(status)
    status = nf90_put_var(ncid, wavgid, wavg, (/1,1,nccall/), (/jmax, ncklimit , 1/))
    if(status /= nf90_noerr) call nchandle_error(status)
    status = nf90_put_var(ncid, thlavgid, thlavg, (/1,1,nccall/), (/jmax, ncklimit , 1/))
    if(status /= nf90_noerr) call nchandle_error(status)
    status = nf90_put_var(ncid, qtavgid, qtavg, (/1,1,nccall/), (/jmax, ncklimit , 1/))
    if(status /= nf90_noerr) call nchandle_error(status)
    status = nf90_put_var(ncid, eavgid, eavg, (/1,1,nccall/), (/jmax, ncklimit , 1/))
    if(status /= nf90_noerr) call nchandle_error(status)


    !calculate variances and store them
    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          uvar(j,k) = uvar(j,k) + (u0(i+1,j+1,k)-uavg(j,k))**2.
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          vvar(j,k) = vvar(j,k) + (v0(i+1,j+1,k)-vavg(j,k))**2.
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          wvar(j,k) = wvar(j,k) + (w0(i+1,j+1,k)-wavg(j,k))**2.
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          thlvar(j,k) = thlvar(j,k) + (thl0(i+1,j+1,k)-thlavg(j,k))**2.
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          qtvar(j,k) = qtvar(j,k) + (qt0(i+1,j+1,k)-qtavg(j,k))**2.
        end do
      end do
    end do

    uvar = uvar / imax
    vvar = vvar / imax
    wvar = wvar / imax
    thlvar = thlvar / imax
    qtvar = qtvar / imax

    status = nf90_put_var(ncid, uvarid, uvar, (/1,1,nccall/), (/jmax, ncklimit , 1/))
    if(status /= nf90_noerr) call nchandle_error(status)
    status = nf90_put_var(ncid, vvarid, vvar, (/1,1,nccall/), (/jmax, ncklimit , 1/))
    if(status /= nf90_noerr) call nchandle_error(status)
    status = nf90_put_var(ncid, wvarid, wvar, (/1,1,nccall/), (/jmax, ncklimit , 1/))
    if(status /= nf90_noerr) call nchandle_error(status)
    status = nf90_put_var(ncid, thlvarid, thlvar, (/1,1,nccall/), (/jmax, ncklimit, 1/))
    if(status /= nf90_noerr) call nchandle_error(status)
    status = nf90_put_var(ncid, qtvarid, qtvar, (/1,1,nccall/), (/jmax, ncklimit, 1/))
    if(status /= nf90_noerr) call nchandle_error(status)

    !calculate covariances and store them
    do k = 2,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          vwcov(j,k) = vwcov(j,k) + (vonw(i+1,j+1,k)-vonwavg(j,k))*(w0(i+1,j+1,k)-wavg(j,k))
        end do
      end do
    end do


    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          if(k==1) then
            vwcovs(j,k) = vwcovs(j,k) - max(ustar(i,j)**2.0, 1.e-10)
          else
            vwcovs(j,k) = vwcovs(j,k) - 0.5*(ekm(i+1,j+1,k)+ekm(i+1,j+1,k-1)) * (0.5*(v0(i+1,j+1,k)+v0(i+1,j+2,k)) - 0.5*(v0(i+1,j+1,k-1)+v0(i+1,j+2,k-1))) / dz
          endif
        end do
      end do
    end do

    do k = 2,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          wthlcov(j,k) = wthlcov(j,k) + (w0(i+1,j+1,k)-wavg(j,k)) * (thl0h(i+1,j+1,k)-thlhavg(j,k))
        end do
      end do
    end do


    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          if(k==1) then
            wthlcovs(j,k) = wthlcovs(j,k) - tstar(i+1,j+1) * max(ustar(i,j),1.e-10)
          else
            wthlcovs(j,k) = wthlcovs(j,k) - 0.5*(ekh(i+1,j+1,k)+ekh(i+1,j+1,k-1)) * (thl0(i+1,j+1,k) - thl0(i+1,j+1,k-1)) / dz
          endif
        end do
      end do
    end do

    do k = 2,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          wqtcov(j,k) = wqtcov(j,k) + (w0(i+1,j+1,k)-wavg(j,k)) * (qt0h(i+1,j+1,k)-qthavg(j,k))
        end do
      end do
    end do

    do k = 1,ncklimit
      do j = 1,jmax
        do i = 1,imax
          !shift prognostic fields one step as 1st column
          !is dummy column because of MPI and periodicity
          if(k==1) then
            wqtcovs(j,k) = wqtcovs(j,k) - qstar(i+1,j+1) * max(ustar(i,j),1.e-10)
          else
            wqtcovs(j,k) = wqtcovs(j,k) - 0.5*(ekh(i+1,j+1,k)+ekh(i+1,j+1,k-1)) * (qt0(i+1,j+1,k) - qt0(i+1,j+1,k-1)) / dz
          endif
        end do
      end do
    end do

    vwcov = vwcov / imax
    vwcovs = vwcovs / imax
    wthlcov = wthlcov / imax
    wthlcovs = wthlcovs / imax
    wqtcov = wqtcov / imax
    wqtcovs = wqtcovs / imax

    status = nf90_put_var(ncid, vwcovid, vwcov, (/1,1,nccall/), (/jmax, ncklimit, 1/))
    if(status /= nf90_noerr) call nchandle_error(status)
    status = nf90_put_var(ncid, vwcovsid, vwcovs, (/1,1,nccall/), (/jmax, ncklimit, 1/))
    if(status /= nf90_noerr) call nchandle_error(status)
    status = nf90_put_var(ncid, wthlcovid, wthlcov, (/1,1,nccall/), (/jmax, ncklimit, 1/))
    if(status /= nf90_noerr) call nchandle_error(status)
    status = nf90_put_var(ncid, wthlcovsid, wthlcovs, (/1,1,nccall/), (/jmax, ncklimit, 1/))
    if(status /= nf90_noerr) call nchandle_error(status)
    status = nf90_put_var(ncid, wqtcovid, wqtcov, (/1,1,nccall/), (/jmax, ncklimit, 1/))
    if(status /= nf90_noerr) call nchandle_error(status)
    status = nf90_put_var(ncid, wqtcovsid, wqtcovs, (/1,1,nccall/), (/jmax, ncklimit, 1/))
    if(status /= nf90_noerr) call nchandle_error(status)

  end subroutine do_netcdfstats

  subroutine exitnetcdfstats

    use typeSizes
    use netcdf
    use modmpi, only : myid

    implicit none

    integer status

    if(.not.(lnetcdf)) return

    write(6,*) "NETCDFSTATS: Closing: ", myid

    status = nf90_close(ncid)
    if (status /= nf90_noerr) call nchandle_error(status)

  end subroutine exitnetcdfstats

  subroutine nchandle_error(status)

    use typeSizes
    use netcdf
    implicit none

    integer, intent(in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if

  end subroutine nchandle_error

end module modnetcdfstats
