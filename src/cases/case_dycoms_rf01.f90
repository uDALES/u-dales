!> \file case_dycoms_rf01.f90
!! DYCOMS-II RF01 case hook: the free-troposphere term of the Stevens et al.
!! (2005) radiation, added to the generic radpar tendency. Switched on with
!! `case_hooks = 'dycoms_rf01'` in &NAMCASES; parameters dycoms_D and
!! dycoms_qt_thres in &NAMCASE_DYCOMS.
!>
!! Ported from the rad_user of DALES dev 075ac4e cases/dycoms_rf02/moduser.f90,
!! with the changes listed at dycoms_rf01_radiation.
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
module case_dycoms_rf01
  use mpi   ! explicit MPI interfaces: gfortran rejects MPI_BCAST calls of different types in one file otherwise
  implicit none
  private
  public :: dycoms_rf01_config, dycoms_rf01_init, dycoms_rf01_radiation, dycoms_rf01_group_present
  save

  ! Case-prefixed names: tools/python/namelists.json maps variable names to groups globally.
  real :: dycoms_D        = 0.   !< large-scale divergence [1/s] (RF01: 3.75e-6); required
  real :: dycoms_qt_thres = 0.   !< qt defining the local inversion height [kg/kg] (RF01: 8e-3); required

contains

  !> Read &NAMCASE_DYCOMS on rank 0; ierr as of the namelist read (< 0: group absent).
  subroutine read_group(ierr)
    use modglobal, only : ifnamopt, fname_options
    implicit none
    integer, intent(out) :: ierr
    namelist/NAMCASE_DYCOMS/ dycoms_D, dycoms_qt_thres

    open(ifnamopt, file=fname_options, status='old', iostat=ierr)
    read (ifnamopt, NAMCASE_DYCOMS, iostat=ierr)
    close(ifnamopt)
  end subroutine read_group

  !> True on rank 0 if &NAMCASE_DYCOMS is in namoptions (for modcases:warn_unlisted_groups).
  logical function dycoms_rf01_group_present()
    implicit none
    integer :: ierr
    call read_group(ierr)
    dycoms_rf01_group_present = (ierr == 0)
  end function dycoms_rf01_group_present

  !> Read &NAMCASE_DYCOMS and check that the case can run (at registration, all ranks).
  subroutine dycoms_rf01_config
    use modglobal,  only : lmoist
    use modraddata, only : iradiation, irad_par
    use modmpi,     only : myid, comm3d, mpierr, my_real, mpi_integer
    implicit none
    integer :: ierr
    character(80) :: msg

    if (myid == 0) call read_group(ierr)
    call MPI_BCAST(ierr, 1, MPI_INTEGER, 0, comm3d, mpierr)   ! all ranks abort together
    if (ierr > 0) then
      write(msg, '(a,i0)') 'Problem in namoptions NAMCASE_DYCOMS, iostat = ', ierr
      call fail(trim(msg))
    end if
    call MPI_BCAST(dycoms_D,        1, MY_REAL, 0, comm3d, mpierr)
    call MPI_BCAST(dycoms_qt_thres, 1, MY_REAL, 0, comm3d, mpierr)

    if (dycoms_D <= 0. .or. dycoms_qt_thres <= 0.) &
      call fail('case hook dycoms_rf01 requires dycoms_D > 0 and dycoms_qt_thres > 0 in &NAMCASE_DYCOMS')
    if (.not. lmoist) call fail('case hook dycoms_rf01 requires lmoist = .true.')
    if (iradiation /= irad_par) call fail('case hook dycoms_rf01 requires iradiation = 2')
  end subroutine dycoms_rf01_config

  !> Warn (rank 0) if no column of the initial state has an inversion, i.e. the term is off.
  subroutine dycoms_rf01_init
    use modglobal, only : ib, ie, jb, je, kb, ke
    use modfields, only : qt0
    use modmpi,    only : myid, comm3d, mpierr, mpi_integer, mpi_max
    implicit none
    integer :: i, j, nloc, nglob

    nloc = 0
    do j = jb, je
      do i = ib, ie
        if (local_zi(qt0(i,j,kb:ke)) > 0.) nloc = 1
      end do
    end do
    call MPI_ALLREDUCE(nloc, nglob, 1, MPI_INTEGER, MPI_MAX, comm3d, mpierr)
    if (myid == 0 .and. nglob == 0) write(0, '(2a)') 'WARNING: case hook dycoms_rf01: no column of the ', &
      'initial qt crosses dycoms_qt_thres below the top level; the free-troposphere term is off'
  end subroutine dycoms_rf01_init

  !> Abort on all ranks with one message from rank 0.
  subroutine fail(msg)
    use modmpi, only : myid, abort_all_ranks
    implicit none
    character(*), intent(in) :: msg
    if (myid == 0) write(0, '(2a)') 'ERROR: ', msg
    call abort_all_ranks(1)
  end subroutine fail

  !> Add the free-troposphere term to thlprad, column by column.
  !!
  !! Flux form (D = dycoms_D): F(z) = rho_i cp D alpha_z [(z-zi)^(4/3)/4 + zi (z-zi)^(1/3)] for
  !! z > zi, 0 below, and thlprad(k) -= (F(zh(k+1)) - F(zh(k))) / (rho_i cp dzf(k)),
  !! i.e. the divergence of the analytic flux between the cell faces. The pointwise
  !! derivative used before, -(D/3)[(z-zi)^(1/3) + zi (z-zi)^(-2/3)], is singular
  !! at z -> zi+; the flux form is finite, conserves the column integral and equals
  !! the pointwise form to O(dz^2) away from zi. rho_i/rho ~ 1 and no 1/exner, as
  !! in DALES (alpha_z = 1 K m^-1/3, so ft_flux below is F / (rho_i cp)).
  !!
  !! Sign: COOLING. DALES adds +D/3[...], which doubles subsidence warming instead
  !! of balancing it. With thl = 297.5 + (z-zi)^(1/3), -dF/dz/(rho_i cp)
  !! = -(D/3) z (z-zi)^(-2/3) = w dthl/dz for w = -D z, i.e. exactly minus the
  !! subsidence warming.
  !!
  !! Two further fixes relative to DALES rad_user: zi is reset per column (DALES
  !! reuses the previous column's value) and the i loop covers the full range
  !! (DALES ran i = 2, j1).
  subroutine dycoms_rf01_radiation
    use modglobal,  only : ib, ie, jb, je, kb, ke, zh, dzf
    use modfields,  only : qt0, IIc
    use modraddata, only : thlprad
    implicit none
    real :: zi
    integer :: i, j, k

    do j = jb, je
      do i = ib, ie
        zi = local_zi(qt0(i,j,kb:ke))
        if (zi <= 0.) cycle
        do k = ke, kb, -1
          if (zh(k+1) <= zi) exit
          thlprad(i,j,k) = thlprad(i,j,k) - IIc(i,j,k) * &
                           (ft_flux(zh(k+1), zi) - ft_flux(zh(k), zi)) / dzf(k)
        end do
      end do
    end do
  end subroutine dycoms_rf01_radiation

  !> Local inversion height: highest downward crossing of dycoms_qt_thres, linear
  !> between levels; -1 if none. A column with qt > dycoms_qt_thres at the top
  !> level ke has no inversion inside the domain and is skipped (-1): no
  !> extrapolation above ke. Otherwise qt(k+1) <= dycoms_qt_thres < qt(k), so the
  !> interpolation never divides by zero.
  real function local_zi(qt)
    use modglobal, only : kb, ke, zf
    implicit none
    real, intent(in) :: qt(kb:ke)
    integer :: k

    local_zi = -1.
    if (qt(ke) > dycoms_qt_thres) return
    do k = ke - 1, kb, -1
      if (qt(k) > dycoms_qt_thres) then
        local_zi = zf(k) + (dycoms_qt_thres - qt(k)) * (zf(k+1) - zf(k)) / (qt(k+1) - qt(k))
        return
      end if
    end do
  end function local_zi

  !> Free-troposphere flux F(z) / (rho_i cp) [K m/s]; 0 at and below zi.
  real function ft_flux(z, zi)
    implicit none
    real, intent(in) :: z, zi

    if (z <= zi) then
      ft_flux = 0.
    else
      ft_flux = dycoms_D * ((z - zi)**(4./3.) / 4. + zi * (z - zi)**(1./3.))
    end if
  end function ft_flux

end module case_dycoms_rf01
