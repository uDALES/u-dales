!> \file out_instant_fields.f90
!! Writes instantaneous 3D fields to unified NetCDF files (one per x-processor column).
!>
!
!! Replaces modfielddump.f90 by gathering y-direction data to the column-head processor
!! (myidy==0) using writeoffset, producing ins_fields.xxx.expnr.nc files.
!!
!! \author Jingzi Huang, ICL (2024-2026)
!! \author Dipanjan Majumdar, ICL (2023-2026)
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
module instant_fields_out

  use modglobal,  only : cexpnr, rk3step, &
                         ib, ie, jb, je, kb, ke, ih, jh, kh, &
                         jtot, lfielddump, tfieldstart, tfielddump, &
                         timee, btime, fieldvars, dyi, dxfi, dzhi
  use modfields,  only : u0, v0, w0, thl0, qt0, ql0, sv0, pres0, div, dudx, dvdy, dwdz, &
                         tau_x, tau_y, tau_z, thl_flux
  ! use modpois,    only : p, pup, pvp, pwp, rhs, dpupdx, dpvpdy, dpwpdz
  use modibm,     only : mask_u, mask_v, mask_w, mask_c
  use modmpi,     only : myidy, cmyidx
  use modstat_nc, only : ncinfo, open_nc, define_nc, writestat_dims_nc, writestat_nc, writeoffset

  implicit none
  private
  public :: ins_field_init, ins_field_main, ins_field_exit
  save

  integer :: xdim, ydim, zdim
  integer :: nfieldvars
  integer :: ncid_fields, nrec_fields

  integer :: ilow, ihigh, jlow, jhigh, klow, khigh
  logical :: lhalos_out
  real    :: out_tnextfielddump

  character(80)                 :: filename_fields
  character(80), dimension(1,4) :: tncname
  character(80), allocatable    :: fldVars(:,:)

  ! Domain pointer type for field variables
  type domainptr
    real, pointer :: point(:,:,:)
  end type domainptr
  type(domainptr), dimension(30) :: pfields

contains

  subroutine ins_field_init
    use instant, only: instant_validate_output_vars
    implicit none
    integer :: n, ydimtot

    if (.not. lfielddump) return

    call instant_validate_output_vars(nfieldvars, fieldvars, 'fieldvars', 'lfielddump')

    allocate(fldVars(nfieldvars,4))

    lhalos_out = .false.

    if (lhalos_out) then
      ilow  = ib - ih
      ihigh = ie + ih
      jlow  = jb - jh
      jhigh = je + jh
      klow  = kb - kh
      khigh = ke + kh
      ydimtot = jtot + 2*jh
    else
      ilow  = ib
      ihigh = ie
      jlow  = jb
      jhigh = je
      klow  = kb
      khigh = ke
      ydimtot = jtot
    end if

    if (tfieldstart .le. btime) then
      out_tnextfielddump = btime
    else
      out_tnextfielddump = tfieldstart
    end if

    xdim = ihigh - ilow + 1
    ydim = jhigh - jlow + 1
    zdim = khigh - klow + 1

    call ncinfo(tncname(1,:), 'time', 'Time', 's', 'time')

    ! Set up variable descriptions and field pointers
    if (lhalos_out) then

      do n = 1, nfieldvars
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
          call ncinfo(fldVars(n,:), 'ql', 'Liquid water mixing ratio', 'kg/kg', 'tttt')
          pfields(n)%point => ql0(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
        case('qt')
          call ncinfo(fldVars(n,:), 'qt', 'Total water mixing ratio', 'kg/kg', 'tttt')
          pfields(n)%point => qt0(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
        case('s1')
          call ncinfo(fldVars(n,:), 's1', 'scalar concentration field 1', 'g/m^3', 'tttt')
          pfields(n)%point => sv0(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh, 1)
        case('s2')
          call ncinfo(fldVars(n,:), 's2', 'scalar concentration field 2', 'g/m^3', 'tttt')
          pfields(n)%point => sv0(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh, 2)
        case('s3')
          call ncinfo(fldVars(n,:), 's3', 'scalar concentration field 3', 'g/m^3', 'tttt')
          pfields(n)%point => sv0(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh, 3)
        case('s4')
          call ncinfo(fldVars(n,:), 's4', 'scalar concentration field 4', 'g/m^3', 'tttt')
          pfields(n)%point => sv0(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh, 4)
        case('s5')
          call ncinfo(fldVars(n,:), 's5', 'scalar concentration field 5', 'g/m^3', 'tttt')
          pfields(n)%point => sv0(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh, 5)
        case('p0')
          call ncinfo(fldVars(n,:), 'p', 'Kinematic pressure field', 'm^2/s^2', 'tttt')
          pfields(n)%point => pres0(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
        ! case('pc')
        !   call ncinfo(fldVars( n,:),'pc','pressure correction','M','tttt')
        !   pfields(n)%point => p(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
        ! case('du')
        !   call ncinfo(fldVars( n,:),'dpupdx','','M','tttt')
        !   pfields(n)%point => dpupdx(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
        ! case('dv')
        !   call ncinfo(fldVars( n,:),'dpvpdy','','M','tttt')
        !   pfields(n)%point => dpvpdy(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
        ! case('dw')
        !   call ncinfo(fldVars( n,:),'dpwpdz','','M','tttt')
        !   pfields(n)%point => dpwpdz(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
        case default
          print *, "Invalid field variables name. Check namoptions setting for 'fieldvars'. &
                    &The variables can be only from the list: u0,v0,w0,p0,th,qt,ql,s1,s2,s3,s4,s5. &
                    &There should not be any space in the string."
          STOP 1
        end select
      end do

    else

      do n = 1, nfieldvars
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
          call ncinfo(fldVars(n,:), 'ql', 'Liquid water mixing ratio', 'kg/kg', 'tttt')
          pfields(n)%point => ql0(ib:ie, jb:je, kb:ke)
        case('qt')
          call ncinfo(fldVars(n,:), 'qt', 'Total water mixing ratio', 'kg/kg', 'tttt')
          pfields(n)%point => qt0(ib:ie, jb:je, kb:ke)
        case('s1')
          call ncinfo(fldVars(n,:), 's1', 'scalar concentration field 1', 'g/m^3', 'tttt')
          pfields(n)%point => sv0(ib:ie, jb:je, kb:ke, 1)
        case('s2')
          call ncinfo(fldVars(n,:), 's2', 'scalar concentration field 2', 'g/m^3', 'tttt')
          pfields(n)%point => sv0(ib:ie, jb:je, kb:ke, 2)
        case('s3')
          call ncinfo(fldVars(n,:), 's3', 'scalar concentration field 3', 'g/m^3', 'tttt')
          pfields(n)%point => sv0(ib:ie, jb:je, kb:ke, 3)
        case('s4')
          call ncinfo(fldVars(n,:), 's4', 'scalar concentration field 4', 'g/m^3', 'tttt')
          pfields(n)%point => sv0(ib:ie, jb:je, kb:ke, 4)
        case('s5')
          call ncinfo(fldVars(n,:), 's5', 'scalar concentration field 5', 'g/m^3', 'tttt')
          pfields(n)%point => sv0(ib:ie, jb:je, kb:ke, 5)
        case('p0')
          call ncinfo(fldVars(n,:), 'p', 'Kinematic pressure field', 'm^2/s^2', 'tttt')
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
        ! case('pc')
        !   call ncinfo(fldVars( n,:),'pc','pressure correction','M','tttt')
        !   pfields(n)%point => p(ib:ie,jb:je,kb:ke)
        ! case('pu')
        !   call ncinfo(fldVars( n,:),'pup','predicted u','M','mttt')
        !   pfields(n)%point => pup(ib:ie,jb:je,kb:ke)
        ! case('pv')
        !   call ncinfo(fldVars( n,:),'pvp','predicted v','M','tmtt')
        !   pfields(n)%point => pvp(ib:ie,jb:je,kb:ke)
        ! case('pw')
        !   call ncinfo(fldVars( n,:),'pwp','predicted w','M','ttmt')
        !   pfields(n)%point => pwp(ib:ie,jb:je,kb:ke)
        ! case('rs')
        !   call ncinfo(fldVars( n,:),'rhs','rhs of poisson equation','M','tttt')
        !   pfields(n)%point => rhs(ib:ie,jb:je,kb:ke)
        ! case('du')
        !   call ncinfo(fldVars( n,:),'dpupdx','','M','tttt')
        !   pfields(n)%point => dpupdx(ib:ie,jb:je,kb:ke)
        ! case('dv')
        !   call ncinfo(fldVars( n,:),'dpvpdy','','M','tttt')
        !   pfields(n)%point => dpvpdy(ib:ie,jb:je,kb:ke)
        ! case('dw')
        !   call ncinfo(fldVars( n,:),'dpwpdz','','M','tttt')
        !   pfields(n)%point => dpwpdz(ib:ie,jb:je,kb:ke)
        case('di')
          call ncinfo(fldVars(n,:), 'div', 'Divergence after pressure correction', 's^-1', 'tttt')
          pfields(n)%point => div(ib:ie, jb:je, kb:ke)
        ! case('ux')
        !   call ncinfo(fldVars( n,:),'dudx','','s^-1','tttt')
        !   pfields(n)%point => dudx(ib:ie,jb:je,kb:ke)
        ! case('vy')
        !   call ncinfo(fldVars( n,:),'dvdy','','s^-1','tttt')
        !   pfields(n)%point => dvdy(ib:ie,jb:je,kb:ke)
        ! case('wz')
        !   call ncinfo(fldVars( n,:),'dwdz','','s^-1','tttt')
        !   pfields(n)%point => dwdz(ib:ie,jb:je,kb:ke)
        case default
          print *, "Invalid field variables name. Check namoptions setting for 'fieldvars'. &
                    &The variables can be only from the list: u0,v0,w0,p0,th,qt,ql,s1,s2,s3,s4,s5,tx,ty,tz,hf,mu,mv,mw,mc,di. &
                    &There should not be any space in the string."
          STOP 1
        end select
      end do

    end if

    ! Create NetCDF file: ins_fields.xxx.xxx.nc (one per x-processor column)
    filename_fields = 'ins_fields.xxx.xxx.nc'
    filename_fields(12:14) = cmyidx
    filename_fields(16:18) = cexpnr

    nrec_fields = 0
    if (myidy == 0) then
      call open_nc(filename_fields, ncid_fields, nrec_fields, n1=xdim, n2=ydimtot, n3=zdim)
      if (nrec_fields == 0) then
        call define_nc(ncid_fields, 1, tncname)
        call writestat_dims_nc(ncid_fields)
      end if
      call define_nc(ncid_fields, nfieldvars, fldVars)
    end if

  end subroutine ins_field_init


  subroutine ins_field_main
    implicit none
    integer :: i, j, k, n

    if (.not. lfielddump) return
    if (.not. rk3step==3)  return
    if (.not. (timee >= out_tnextfielddump)) return

    ! To be uncommented it fo be output via ins_field_out
    ! do k = kb, ke
    !   do j = jb, je
    !     do i = ib, ie
    !       dudx(i,j,k) = (u0(i+1,j,k) - u0(i,j,k)) * dxfi(i)
    !       dvdy(i,j,k) = (v0(i,j+1,k) - v0(i,j,k)) * dyi
    !       dwdz(i,j,k) = (w0(i,j,k+1) - w0(i,j,k)) * dzhi(k)
    !       div(i,j,k)  = (u0(i+1,j,k) - u0(i,j,k)) * dxfi(i) + &
    !                      (v0(i,j+1,k) - v0(i,j,k)) * dyi + &
    !                      (w0(i,j,k+1) - w0(i,j,k)) * dzhi(k)
    !     end do
    !   end do
    ! end do

    out_tnextfielddump = out_tnextfielddump + tfielddump

    ! Write time
    if (myidy == 0) call writestat_nc(ncid_fields, 'time', timee, nrec_fields, .true.)

    ! Write each field variable using writeoffset
    do n = 1, nfieldvars
      call writeoffset(ncid_fields, trim(fldVars(n,1)), pfields(n)%point, nrec_fields, xdim, ydim, zdim)
    end do

  end subroutine ins_field_main


  subroutine ins_field_exit
    use modstat_nc, only : exitstat_nc
    implicit none
    if (.not. lfielddump) return
    if (myidy == 0) call exitstat_nc(ncid_fields)
    if (allocated(fldVars)) deallocate(fldVars)
  end subroutine ins_field_exit

end module instant_fields_out
