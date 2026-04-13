!> \file out_fields.f90
!! Writes instantaneous 3D fields to unified NetCDF files (one per x-processor column).
!>
!
!! Replaces modfielddump.f90 by gathering y-direction data to the column-head processor
!! (myidy==0) using writeoffset, producing ins_fields.xxx.expnr.nc files.
!!
!! \author Jingzi Huang, ICL (2024-2026)
!! \todo documentation
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
module out_fields

  use modglobal,  only : cexpnr, rk3step, ltempeq, lmoist, nsv, &
                         ib, ie, jb, je, kb, ke, ih, jh, kh, &
                         jtot, lfielddump, tfielddump, &
                         timee, fieldvars, dyi, dxfi, dzhi
  use modfields,  only : u0, v0, w0, thl0, qt0, ql0, sv0, pres0, div, dudx, dvdy, dwdz, &
                         tau_x, tau_y, tau_z, thl_flux
  use modibm,     only : mask_u, mask_v, mask_w, mask_c
  use modmpi,     only : myid, myidy, cmyidx, cmyidy, comm3d, mpierr, my_real, NPROCY
  use modstat_nc, only : ncinfo, open_nc, define_nc, writestat_dims_nc, writestat_nc, writeoffset

  implicit none
  private
  public :: out_initfields, out_fields_main, out_exitfields
  save

  integer :: xdim, ydim, zdim
  integer :: nvar
  integer :: ncid_fields, nrec_fields

  integer :: ilow, ihigh, jlow, jhigh, klow, khigh
  logical :: lhalos
  real    :: out_tnextfielddump

  character(80)              :: filename_fields
  character(80), dimension(1,4) :: tncname
  character(80), allocatable :: fldVars(:,:)

  ! Domain pointer type for field variables
  type domainptr
    real, pointer :: point(:,:,:)
  end type domainptr
  type(domainptr), dimension(30) :: pfields

contains

  subroutine out_initfields
    use modmpi, only : mpi_logical, mpi_integer, mpi_character
    implicit none
    integer :: ierr, n, jtot_out

    nvar = (LEN(trim(fieldvars))+1)/3

    if (nvar == 0) then
      lfielddump = .false.
      if (myid == 0) write(*,*) 'empty fieldvars therefore lfielddump = .false. and no instantaneous fields outputted'
      return
    end if

    allocate(fldVars(nvar,4))

    lhalos = .false.

    if (lhalos) then
      ilow  = ib - ih
      ihigh = ie + ih
      jlow  = jb - jh
      jhigh = je + jh
      klow  = kb - kh
      khigh = ke + kh
    else
      ilow  = ib
      ihigh = ie
      jlow  = jb
      jhigh = je
      klow  = kb
      khigh = ke
    end if

    call MPI_BCAST(klow,        1, MPI_INTEGER, 0, comm3d, ierr)
    call MPI_BCAST(khigh,       1, MPI_INTEGER, 0, comm3d, ierr)
    call MPI_BCAST(lfielddump,  1, MPI_LOGICAL, 0, comm3d, ierr)
    call MPI_BCAST(nvar,        1, MPI_INTEGER, 0, comm3d, mpierr)

    if (.not. lfielddump) return

    out_tnextfielddump = tfielddump

    xdim = ihigh - ilow + 1
    ydim = jhigh - jlow + 1
    zdim = khigh - klow + 1

    if (lhalos) then
      jtot_out = NPROCY * ydim
    else
      jtot_out = jtot
    end if

    call ncinfo(tncname(1,:), 'time', 'Time', 's', 'time')

    ! Set up variable descriptions and field pointers
    if (lhalos) then

      do n = 1, nvar
        select case(fieldvars(3*n-2:3*n-1))
        case('u0')
          call ncinfo(fldVars(n,:), 'u', 'West-East velocity', 'm/s', 'mttt')
          pfields(n)%point => u0(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
        case('v0')
          call ncinfo(fldVars(n,:), 'v', 'South-North velocity', 'm/s', 'tmtt')
          pfields(n)%point => v0(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
        case('w0')
          call ncinfo(fldVars(n,:), 'w', 'Vertical velocity', 'm/s', 'ttmt')
          pfields(n)%point => w0(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
        case('th')
          call ncinfo(fldVars(n,:), 'thl', 'Liquid water potential temperature', 'K', 'tttt')
          pfields(n)%point => thl0(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
        case('ql')
          call ncinfo(fldVars(n,:), 'ql', 'Liquid water mixing ratio', '1e-5kg/kg', 'tttt')
          pfields(n)%point => ql0(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
        case('qt')
          call ncinfo(fldVars(n,:), 'qt', 'Total water mixing ratio', '1e-5kg/kg', 'tttt')
          pfields(n)%point => qt0(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
        case('s1')
          call ncinfo(fldVars(n,:), 'sca1', 'scalar 1', 'M', 'tttt')
          pfields(n)%point => sv0(ib-ih:ie+ih, jb-jh:je+jh, kb:ke, 1)
        case('s2')
          call ncinfo(fldVars(n,:), 'sca2', 'scalar 2', 'M', 'tttt')
          pfields(n)%point => sv0(ib-ih:ie+ih, jb-jh:je+jh, kb:ke, 2)
        case('s3')
          call ncinfo(fldVars(n,:), 'sca3', 'scalar 3', 'M', 'tttt')
          pfields(n)%point => sv0(ib-ih:ie+ih, jb-jh:je+jh, kb:ke, 3)
        case('s4')
          call ncinfo(fldVars(n,:), 'sca4', 'scalar 4', 'M', 'tttt')
          pfields(n)%point => sv0(ib-ih:ie+ih, jb-jh:je+jh, kb:ke, 4)
        case('s5')
          call ncinfo(fldVars(n,:), 'sca5', 'scalar 5', 'M', 'tttt')
          pfields(n)%point => sv0(ib-ih:ie+ih, jb-jh:je+jh, kb:ke, 5)
        case('p0')
          call ncinfo(fldVars(n,:), 'pres', 'pressure field', 'M', 'tttt')
          pfields(n)%point => pres0(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
        case default
          call ncinfo(fldVars(n,:), 'u', 'West-East velocity', 'm/s', 'mttt')
          pfields(n)%point => u0(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
        end select
      end do

    else

      do n = 1, nvar
        select case(fieldvars(3*n-2:3*n-1))
        case('u0')
          call ncinfo(fldVars(n,:), 'u', 'West-East velocity', 'm/s', 'mttt')
          pfields(n)%point => u0(ib:ie, jb:je, kb:ke)
        case('v0')
          call ncinfo(fldVars(n,:), 'v', 'South-North velocity', 'm/s', 'tmtt')
          pfields(n)%point => v0(ib:ie, jb:je, kb:ke)
        case('w0')
          call ncinfo(fldVars(n,:), 'w', 'Vertical velocity', 'm/s', 'ttmt')
          pfields(n)%point => w0(ib:ie, jb:je, kb:ke)
        case('th')
          call ncinfo(fldVars(n,:), 'thl', 'Liquid water potential temperature', 'K', 'tttt')
          pfields(n)%point => thl0(ib:ie, jb:je, kb:ke)
        case('ql')
          call ncinfo(fldVars(n,:), 'ql', 'Liquid water mixing ratio', '1e-5kg/kg', 'tttt')
          pfields(n)%point => ql0(ib:ie, jb:je, kb:ke)
        case('qt')
          call ncinfo(fldVars(n,:), 'qt', 'Total water mixing ratio', '1e-5kg/kg', 'tttt')
          pfields(n)%point => qt0(ib:ie, jb:je, kb:ke)
        case('s1')
          call ncinfo(fldVars(n,:), 'sca1', 'scalar 1', 'M', 'tttt')
          pfields(n)%point => sv0(ib:ie, jb:je, kb:ke, 1)
        case('s2')
          call ncinfo(fldVars(n,:), 'sca2', 'scalar 2', 'M', 'tttt')
          pfields(n)%point => sv0(ib:ie, jb:je, kb:ke, 2)
        case('s3')
          call ncinfo(fldVars(n,:), 'sca3', 'scalar 3', 'M', 'tttt')
          pfields(n)%point => sv0(ib:ie, jb:je, kb:ke, 3)
        case('s4')
          call ncinfo(fldVars(n,:), 'sca4', 'scalar 4', 'M', 'tttt')
          pfields(n)%point => sv0(ib:ie, jb:je, kb:ke, 4)
        case('s5')
          call ncinfo(fldVars(n,:), 'sca5', 'scalar 5', 'M', 'tttt')
          pfields(n)%point => sv0(ib:ie, jb:je, kb:ke, 5)
        case('p0')
          call ncinfo(fldVars(n,:), 'pres', 'pressure field', 'M', 'tttt')
          pfields(n)%point => pres0(ib:ie, jb:je, kb:ke)
        case('tx')
          call ncinfo(fldVars(n,:), 'tau_x', 'stress x', 'M', 'mttt')
          pfields(n)%point => tau_x(ib:ie, jb:je, kb:ke)
        case('ty')
          call ncinfo(fldVars(n,:), 'tau_y', 'stress y', 'M', 'tmtt')
          pfields(n)%point => tau_y(ib:ie, jb:je, kb:ke)
        case('tz')
          call ncinfo(fldVars(n,:), 'tau_z', 'stress z', 'M', 'ttmt')
          pfields(n)%point => tau_z(ib:ie, jb:je, kb:ke)
        case('hf')
          call ncinfo(fldVars(n,:), 'thl_flux', 'heat flux', 'M', 'tttt')
          pfields(n)%point => thl_flux(ib:ie, jb:je, kb:ke)
        case('mu')
          call ncinfo(fldVars(n,:), 'mask_u', 'mask u', 'M', 'mttt')
          pfields(n)%point => mask_u(ib:ie, jb:je, kb:ke)
        case('mv')
          call ncinfo(fldVars(n,:), 'mask_v', 'mask v', 'M', 'tmtt')
          pfields(n)%point => mask_v(ib:ie, jb:je, kb:ke)
        case('mw')
          call ncinfo(fldVars(n,:), 'mask_w', 'mask w', 'M', 'ttmt')
          pfields(n)%point => mask_w(ib:ie, jb:je, kb:ke)
        case('mc')
          call ncinfo(fldVars(n,:), 'mask_c', 'mask c', 'M', 'tttt')
          pfields(n)%point => mask_c(ib:ie, jb:je, kb:ke)
        case('di')
          call ncinfo(fldVars(n,:), 'div', 'Divergence after pressure correction', 'M', 'tttt')
          pfields(n)%point => div(ib:ie, jb:je, kb:ke)
        case default
          call ncinfo(fldVars(n,:), 'u', 'West-East velocity', 'm/s', 'mttt')
          pfields(n)%point => u0(ib:ie, jb:je, kb:ke)
        end select
      end do

    end if

    ! Create NetCDF file: ins_fields.xxx.xxx.nc (one per x-processor column)
    filename_fields = 'ins_fields.xxx.xxx.nc'
    filename_fields(12:14) = cmyidx
    filename_fields(16:18) = cexpnr

    nrec_fields = 0
    if (myidy == 0) then
      call open_nc(filename_fields, ncid_fields, nrec_fields, n1=xdim, n2=jtot_out, n3=zdim)
      if (nrec_fields == 0) then
        call define_nc(ncid_fields, 1, tncname)
        call writestat_dims_nc(ncid_fields)
      end if
      call define_nc(ncid_fields, nvar, fldVars)
    end if

  end subroutine out_initfields


  subroutine out_fields_main
    implicit none
    integer :: i, j, k, n

    if (.not. lfielddump) return
    if (.not. ((timee >= out_tnextfielddump) .or. (rk3step == 0))) return
    if (rk3step /= 3 .and. rk3step /= 0) return

    ! Compute divergence (same as modfielddump)
    do k = kb, ke
      do j = jb, je
        do i = ib, ie
          dudx(i,j,k) = (u0(i+1,j,k) - u0(i,j,k)) * dxfi(i)
          dvdy(i,j,k) = (v0(i,j+1,k) - v0(i,j,k)) * dyi
          dwdz(i,j,k) = (w0(i,j,k+1) - w0(i,j,k)) * dzhi(k)
          div(i,j,k)  = (u0(i+1,j,k) - u0(i,j,k)) * dxfi(i) + &
                         (v0(i,j+1,k) - v0(i,j,k)) * dyi + &
                         (w0(i,j,k+1) - w0(i,j,k)) * dzhi(k)
        end do
      end do
    end do

    if (rk3step == 3) out_tnextfielddump = out_tnextfielddump + tfielddump

    ! Write time
    if (myidy == 0) call writestat_nc(ncid_fields, 'time', timee, nrec_fields, .true.)

    ! Write each field variable using writeoffset
    do n = 1, nvar
      call writeoffset(ncid_fields, trim(fldVars(n,1)), pfields(n)%point, nrec_fields, xdim, ydim, zdim)
    end do

  end subroutine out_fields_main


  subroutine out_exitfields
    use modstat_nc, only : exitstat_nc
    implicit none

    if (.not. lfielddump) return
    if (myidy == 0) call exitstat_nc(ncid_fields)
    if (allocated(fldVars)) deallocate(fldVars)
  end subroutine out_exitfields

end module out_fields
