!> \file out_instants.f90
!! Writes instantaneous planes of various fields to different NetCDF files.
!>
!
!! Inspired from the slice routines in uDALES v2.2.0 modstatsdump.f90 written by Sam O. Owens, ICL (2024).
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
module instant
  use modglobal,  only : cexpnr, rk3step, ltempeq, lmoist, nsv,  &
                         fieldvars, lfielddump, &
                         slicevars, lislicedump, islice, nislice, ljslicedump, jslice, njslice, lkslicedump, kslice, nkslice, &
                         probevars, lprobedump, iprobe, jprobe, kprobe, nprobe, &
                         ib, ie, jb, je, kb, ke, ih, jh, kh, itot, jtot, ktot, &
                         tinstantstart, tinstantdump, tstatstart, tsample, dt, timee, btime, runtime, &
                         xf, yf, zf, xh, yh, zh
  use modfields,  only : um, vm, wm, thlm, qtm, svm, ql0, pres0, &
                         div, dudx, dvdy, dwdz, &
                         tau_x, tau_y, tau_z, thl_flux
  ! use modpois,    only : p, pup, pvp, pwp, rhs, dpupdx, dpvpdy, dpwpdz
  use modibm,     only : mask_u, mask_v, mask_w, mask_c
  use modmpi,     only : myid, myidx, myidy, cmyidx, cmyidy, comm3d, mpierr, my_real, nprocx, nprocy
  use decomp_2d,  only : zstart, zend
  use modstat_nc, only : ncinfo, open_nc, define_nc, writestat_dims_nc, writestat_nc, writeoffset, writeoffset_1dx
  use mpi
  use netcdf
  implicit none
  private
  public :: instant_init, instant_main, instant_exit
  save

  real    :: tnextinstantdump

  !! Field writing variables
  integer :: xdimfield, ydimfield, zdimfield
  integer                       :: nfieldvars
  character(80)                 :: filenamefield
  character(80), dimension(1,4) :: fieldtimeVar
  character(80), allocatable    :: fldVars(:,:)
  integer                       :: ncidfield, nrecfield

  ! Domain pointer type for field variables
  type domainptr
    real, pointer :: point(:,:,:)
  end type domainptr
  type(domainptr), dimension(30) :: pfields

  !! Slice writing variables
  integer :: xdim, ydim, zdim, kdim  ! Added kdim for multiple slices
  integer :: local_nislice, local_njslice  ! Number of islices on this X-processor (saved from create phase)
  integer                    :: nslicevars
  character(80)              :: filenameislice
  character(80)              :: filenamejslice
  character(80)              :: filenamekslice
  character(80)              :: slicetimeVar(1,4)
  character(80), allocatable :: isliceVars(:,:)
  character(80), allocatable :: jsliceVars(:,:)
  character(80), allocatable :: ksliceVars(:,:)
  integer                    :: ncidislice, nrecislice, ncidjslice, nrecjslice, ncidkslice, nreckslice

  !! Probe writing variables
  integer :: nprobevars
  character(80) :: filenameprobe
  integer :: ncidprobe
  integer :: nrecprobe = 0
  integer :: varid_time
  integer, allocatable :: varid_vals(:)

  contains

    subroutine instant_init
      implicit none
      if(.not.(lfielddump .or. lislicedump .or. ljslicedump .or. lkslicedump .or. lprobedump)) return

      if (tinstantstart .le. btime) then
        tnextinstantdump = btime
      else
        tnextinstantdump = tinstantstart
      end if

      call instant_field_init
      call instant_slice_init
      call instant_probe_init
    end subroutine instant_init

    subroutine instant_main
      implicit none
      if(.not.(lfielddump .or. lislicedump .or. ljslicedump .or. lkslicedump .or. lprobedump)) return
      
      if (.not. (timee >= tnextinstantdump)) return
      if (.not. rk3step==3)  return

      call instant_field_main
      call instant_slice_main
      call instant_probe_main

      tnextinstantdump = tnextinstantdump + tinstantdump
    end subroutine instant_main

    subroutine instant_exit
      implicit none
      call instant_field_exit
    end subroutine instant_exit

    subroutine instant_field_init
      implicit none
      integer :: ilow, ihigh, jlow, jhigh, klow, khigh
      integer :: n, ydimtot
      logical :: lhalos_out

      if (.not. lfielddump) return

      if(runtime <= tinstantstart) then
        if(myid==0) then
          write(*,*) "ERROR: no instantaneous 3D fields will be written as runtime <= tinstantstart. Note that runtime &
                      &must be greater than tinstantstart for wiriting ins_field.* files."
          write(*,*) "You have used runtime = ", runtime, ", tinstantstart = ", tinstantstart
          write(*,*) "Either correct the time settings or change all the lfielddump flag to false."
          stop 1
        end if
      end if

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

      xdimfield = ihigh - ilow + 1
      ydimfield = jhigh - jlow + 1
      zdimfield = khigh - klow + 1

      call ncinfo(fieldtimeVar(1,:), 'time', 'Time', 's', 'time')

      ! Set up variable descriptions and field pointers
      if (lhalos_out) then

        do n = 1, nfieldvars
          select case(fieldvars(3*n-2:3*n-1))
          case('u0')
            call ncinfo(fldVars(n,:), 'u', 'West-East velocity', 'm/s', 'mttt')
            pfields(n)%point => um(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
          case('v0')
            call ncinfo(fldVars(n,:), 'v', 'South-North velocity', 'm/s', 'tmtt')
            pfields(n)%point => vm(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
          case('w0')
            call ncinfo(fldVars(n,:), 'w', 'Vertical velocity', 'm/s', 'ttmt')
            pfields(n)%point => wm(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
          case('th')
            call ncinfo(fldVars(n,:), 'thl', 'Liquid water potential temperature', 'K', 'tttt')
            pfields(n)%point => thlm(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
          case('ql')
            call ncinfo(fldVars(n,:), 'ql', 'Liquid water mixing ratio', 'kg/kg', 'tttt')
            pfields(n)%point => ql0(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
          case('qt')
            call ncinfo(fldVars(n,:), 'qt', 'Total water mixing ratio', 'kg/kg', 'tttt')
            pfields(n)%point => qtm(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh)
          case('s1')
            call ncinfo(fldVars(n,:), 's1', 'scalar concentration field 1', 'g/m^3', 'tttt')
            pfields(n)%point => svm(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh, 1)
          case('s2')
            call ncinfo(fldVars(n,:), 's2', 'scalar concentration field 2', 'g/m^3', 'tttt')
            pfields(n)%point => svm(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh, 2)
          case('s3')
            call ncinfo(fldVars(n,:), 's3', 'scalar concentration field 3', 'g/m^3', 'tttt')
            pfields(n)%point => svm(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh, 3)
          case('s4')
            call ncinfo(fldVars(n,:), 's4', 'scalar concentration field 4', 'g/m^3', 'tttt')
            pfields(n)%point => svm(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh, 4)
          case('s5')
            call ncinfo(fldVars(n,:), 's5', 'scalar concentration field 5', 'g/m^3', 'tttt')
            pfields(n)%point => svm(ib-ih:ie+ih, jb-jh:je+jh, kb-kh:ke+kh, 5)
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
            pfields(n)%point => um(ib:ie, jb:je, kb:ke)
          case('v0')
            call ncinfo(fldVars(n,:), 'v', 'South-North velocity', 'm/s', 'tmtt')
            pfields(n)%point => vm(ib:ie, jb:je, kb:ke)
          case('w0')
            call ncinfo(fldVars(n,:), 'w', 'Vertical velocity', 'm/s', 'ttmt')
            pfields(n)%point => wm(ib:ie, jb:je, kb:ke)
          case('th')
            call ncinfo(fldVars(n,:), 'thl', 'Liquid water potential temperature', 'K', 'tttt')
            pfields(n)%point => thlm(ib:ie, jb:je, kb:ke)
          case('ql')
            call ncinfo(fldVars(n,:), 'ql', 'Liquid water mixing ratio', 'kg/kg', 'tttt')
            pfields(n)%point => ql0(ib:ie, jb:je, kb:ke)
          case('qt')
            call ncinfo(fldVars(n,:), 'qt', 'Total water mixing ratio', 'kg/kg', 'tttt')
            pfields(n)%point => qtm(ib:ie, jb:je, kb:ke)
          case('s1')
            call ncinfo(fldVars(n,:), 's1', 'scalar concentration field 1', 'g/m^3', 'tttt')
            pfields(n)%point => svm(ib:ie, jb:je, kb:ke, 1)
          case('s2')
            call ncinfo(fldVars(n,:), 's2', 'scalar concentration field 2', 'g/m^3', 'tttt')
            pfields(n)%point => svm(ib:ie, jb:je, kb:ke, 2)
          case('s3')
            call ncinfo(fldVars(n,:), 's3', 'scalar concentration field 3', 'g/m^3', 'tttt')
            pfields(n)%point => svm(ib:ie, jb:je, kb:ke, 3)
          case('s4')
            call ncinfo(fldVars(n,:), 's4', 'scalar concentration field 4', 'g/m^3', 'tttt')
            pfields(n)%point => svm(ib:ie, jb:je, kb:ke, 4)
          case('s5')
            call ncinfo(fldVars(n,:), 's5', 'scalar concentration field 5', 'g/m^3', 'tttt')
            pfields(n)%point => svm(ib:ie, jb:je, kb:ke, 5)
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
          ! case('di')
          !   call ncinfo(fldVars(n,:), 'div', 'Divergence after pressure correction', 's^-1', 'tttt')
          !   pfields(n)%point => div(ib:ie, jb:je, kb:ke)
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

      ! Create NetCDF file: ins_field.xxx.xxx.nc (one per x-processor column)
      filenamefield = 'ins_field.xxx.xxx.nc'
      filenamefield(11:13) = cmyidx
      filenamefield(15:17) = cexpnr

      nrecfield = 0
      if (myidy == 0) then
        call open_nc(filenamefield, ncidfield, nrecfield, n1=xdimfield, n2=ydimtot, n3=zdimfield)
        if (nrecfield == 0) then
          call define_nc(ncidfield, 1, fieldtimeVar)
          call writestat_dims_nc(ncidfield)
        end if
        call define_nc(ncidfield, nfieldvars, fldVars)
      end if

    end subroutine instant_field_init

    subroutine instant_field_main
      ! use modglobal, only : dxfi, dyi, dzhi
      implicit none
      integer :: i, j, k, n

      if (.not. lfielddump) return

      ! To be uncommented it fo be output via ins_field_out
      ! do k = kb, ke
      !   do j = jb, je
      !     do i = ib, ie
      !       dudx(i,j,k) = (um(i+1,j,k) - um(i,j,k)) * dxfi(i)
      !       dvdy(i,j,k) = (vm(i,j+1,k) - vm(i,j,k)) * dyi
      !       dwdz(i,j,k) = (wm(i,j,k+1) - wm(i,j,k)) * dzhi(k)
      !       div(i,j,k)  = (um(i+1,j,k) - um(i,j,k)) * dxfi(i) + &
      !                      (vm(i,j+1,k) - vm(i,j,k)) * dyi + &
      !                      (wm(i,j,k+1) - wm(i,j,k)) * dzhi(k)
      !     end do
      !   end do
      ! end do

      ! Write time
      if (myidy == 0) call writestat_nc(ncidfield, 'time', timee, nrecfield, .true.)

      ! Write each field variable using writeoffset
      do n = 1, nfieldvars
        call writeoffset(ncidfield, trim(fldVars(n,1)), pfields(n)%point, nrecfield, xdimfield, ydimfield, zdimfield)
      end do
    end subroutine instant_field_main

    subroutine instant_field_exit
      use modstat_nc, only : exitstat_nc
      implicit none
      if (.not. lfielddump) return
      if (myidy == 0) call exitstat_nc(ncidfield)
      if (allocated(fldVars)) deallocate(fldVars)
    end subroutine instant_field_exit


    subroutine instant_slice_init
      implicit none

      if(.not.(lislicedump .or. ljslicedump .or. lkslicedump)) return

      ! Validate slice inputs before proceeding
      call instant_validate_out_time('slices')
      call instant_validate_output_vars(nslicevars, slicevars, 'slicevars', 'one or more of lislicedump, ljslicedump and lkslicedump')  ! Validate slicevars only if at least one slice dump is true
      call instant_validate_slice_inputs(lkslicedump, nkslice, kslice, 'kslice')
      call instant_validate_slice_inputs(lislicedump, nislice, islice, 'islice') 
      call instant_validate_slice_inputs(ljslicedump, njslice, jslice, 'jslice')
      
      xdim = ie-ib+1
      ydim = je-jb+1
      zdim = ke-kb+1
      kdim = nkslice  ! Set kdim to number of kslices

      call ncinfo(slicetimeVar( 1,:), 'time', 'Time', 's', 'time')

      if (lislicedump) then
        allocate(isliceVars(nslicevars,4))
        call instant_ncdescription_islice
        call instant_create_ncislice    !> Generate sliced NetCDF: ins_islice.xxx.xxx.nc
        deallocate(isliceVars)
      end if
      
      if (ljslicedump) then
        allocate(jsliceVars(nslicevars,4))
        call instant_ncdescription_jslice
        call instant_create_ncjslice    ! Unified jslice creation (per myidy file, myidx==0)
        deallocate(jsliceVars)
      end if

      if (lkslicedump) then
        allocate(ksliceVars(nslicevars,4))
        call instant_ncdescription_kslice
        call instant_create_nckslice    !> Generate sliced NetCDF: ins_kslice.xxx.xxx.nc
        deallocate(ksliceVars)
      end if
    end subroutine instant_slice_init

    subroutine instant_slice_main
      implicit none
      if (lislicedump) call instant_write_islice
      if (ljslicedump) call instant_write_jslice
      if (lkslicedump) call instant_write_kslice
    end subroutine instant_slice_main

    subroutine instant_ncdescription_islice
      implicit none
      integer :: n
      do n=1,nslicevars
        select case(slicevars(3*n-2:3*n-1))
          case('u0')
            call ncinfo( isliceVars(n,:), 'u' , 'Streamwise velocity' , 'm/s' , 'tttt' )
          case('v0')
            call ncinfo( isliceVars(n,:), 'v' , 'Spanwise velocity'   , 'm/s' , 'tmtt' )
          case('w0')
            call ncinfo( isliceVars(n,:), 'w' , 'Vertical velocity'   , 'm/s' , 'ttmt' )
          case('p0')
            call ncinfo( isliceVars(n,:), 'p' , 'Kinematic Pressure'  , 'm^2/s^2' , 'tttt' )
          case('th')
            if (ltempeq) call ncinfo( isliceVars(n,:), 'thl' , 'Potential temperature' , 'K'     , 'tttt' )
          case('qt')
            if (lmoist)  call ncinfo( isliceVars(n,:), 'qt'  , 'Specific humidity'     , 'kg/kg' , 'tttt' )
          case('s1')
            if (nsv>0)   call ncinfo( isliceVars(n,:), 's1'  , 'Concentration field 1' , 'g/m^3' , 'tttt' )
          case('s2')
            if (nsv>1)   call ncinfo( isliceVars(n,:), 's2'  , 'Concentration field 2' , 'g/m^3' , 'tttt' )
          case('s3')
            if (nsv>2)   call ncinfo( isliceVars(n,:), 's3'  , 'Concentration field 3' , 'g/m^3' , 'tttt' )
          case('s4')
            if (nsv>3)   call ncinfo( isliceVars(n,:), 's4'  , 'Concentration field 4' , 'g/m^3' , 'tttt' )
          case default
            print *, "Invalid slice variables name. Check namoptions setting for 'slicevars'. &
                      &The variables can be only from the list: u0,v0,w0,p0,th,qt,s1,s2,s3,s4. &
                      &There should not be any space in the string. "
            STOP 1
        end select
      end do
    end subroutine instant_ncdescription_islice

    subroutine instant_ncdescription_jslice
      implicit none
      integer :: n
      do n=1,nslicevars
        select case(slicevars(3*n-2:3*n-1))
          case('u0')
            call ncinfo( jsliceVars(n,:), 'u' , 'Streamwise velocity' , 'm/s' , 'mttt' )
          case('v0')
            call ncinfo( jsliceVars(n,:), 'v' , 'Spanwise velocity '  , 'm/s' , 'tttt' )
          case('w0')
            call ncinfo( jsliceVars(n,:), 'w' , 'Vertical velocity'   , 'm/s' , 'ttmt' )
          case('p0')
            call ncinfo( jsliceVars(n,:), 'p' , 'Kinematic Pressure'  , 'm^2/s^2' , 'tttt' )
          case('th')
            if (ltempeq) call ncinfo( jsliceVars(n,:), 'thl' , 'Potential temperature' , 'K'     , 'tttt' )
          case('qt')
            if (lmoist)  call ncinfo( jsliceVars(n,:), 'qt'  , 'Specific humidity'     , 'kg/kg' , 'tttt' )
          case('s1')
            if (nsv>0)   call ncinfo( jsliceVars(n,:), 's1'  , 'Concentration field 1' , 'g/m^3' , 'tttt' )
          case('s2')
            if (nsv>1)   call ncinfo( jsliceVars(n,:), 's2'  , 'Concentration field 2' , 'g/m^3' , 'tttt' )
          case('s3')
            if (nsv>2)   call ncinfo( jsliceVars(n,:), 's3'  , 'Concentration field 3' , 'g/m^3' , 'tttt' )
          case('s4')
            if (nsv>3)   call ncinfo( jsliceVars(n,:), 's4'  , 'Concentration field 4' , 'g/m^3' , 'tttt' )
          case default
            print *, "Invalid slice variables name. Check namoptions setting for 'slicevars'. &
                      &The variables can be only from the list: u0,v0,w0,p0,th,qt,s1,s2,s3,s4. &
                      &There should not be any space in the string."
            STOP 1
        end select
      end do
    end subroutine instant_ncdescription_jslice

    subroutine instant_ncdescription_kslice
      implicit none
      integer :: n
      do n=1,nslicevars
        select case(slicevars(3*n-2:3*n-1))
          case('u0')
            call ncinfo( ksliceVars(n,:), 'u' , 'Streamwise velocity' , 'm/s' , 'mttt' )
          case('v0')
            call ncinfo( ksliceVars(n,:), 'v' , 'Spanwise velocity'   , 'm/s' , 'tmtt' )
          case('w0')
            call ncinfo( ksliceVars(n,:), 'w' , 'Vertical velocity'   , 'm/s' , 'tttt' )
          case('p0')
            call ncinfo( ksliceVars(n,:), 'p' , 'Kinematic Pressure'  , 'm^2/s^2' , 'tttt' )
          case('th')
            if (ltempeq) call ncinfo( ksliceVars( n,:), 'thl' , 'Potential temperature' , 'K'     , 'tttt' )
          case('qt')
            if (lmoist)  call ncinfo( ksliceVars( n,:), 'qt'  , 'Specific humidity'     , 'kg/kg' , 'tttt' )
          case('s1')
            if (nsv>0)   call ncinfo( ksliceVars( n,:), 's1'  , 'Concentration field 1' , 'g/m^3' , 'tttt' )
          case('s2')
            if (nsv>1)   call ncinfo( ksliceVars( n,:), 's2'  , 'Concentration field 2' , 'g/m^3' , 'tttt' )
          case('s3')
            if (nsv>2)   call ncinfo( ksliceVars( n,:), 's3'  , 'Concentration field 3' , 'g/m^3' , 'tttt' )
          case('s4')
            if (nsv>3)   call ncinfo( ksliceVars( n,:), 's4'  , 'Concentration field 4' , 'g/m^3' , 'tttt' )
          case default
            print *, "Invalid slice variables name. Check namoptions setting for 'slicevars'. &
                      &The variables can be only from the list: u0,v0,w0,p0,th,qt,s1,s2,s3,s4. &
                      &There should not be any space in the string."
            STOP 1
        end select
      end do
    end subroutine instant_ncdescription_kslice

    !! ## %% 1D parallel output creation for islice (y-direction processes only)
    subroutine instant_create_ncislice
      implicit none
      integer :: i
      logical :: has_islice
      
      ! Count how many islices are in this X-processor's range (save to module variable)
      local_nislice = 0
      has_islice = .false.
      do i = 1, nislice
        ! if ( (islice(i)-1)/(itot/nprocx) == myidx) then
        if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
          local_nislice = local_nislice + 1
          has_islice = .true.
        end if
      end do
      
      ! Only processors with islices AND myidy==0 create files
      if (has_islice .and. myidy == 0) then
        filenameislice = 'ins_islice.xxx.xxx.nc'
        filenameislice(12:14) = cmyidx   ! X-processor ID
        filenameislice(16:18) = cexpnr  ! Experiment number

        nrecislice = 0
        ! Use local_nislice (module variable) as first dimension, jtot as global y-dimension
        call open_nc(filenameislice, ncidislice, nrecislice, n1=local_nislice, n2=jtot, n3=zdim)
        if (nrecislice == 0) then
          call define_nc(ncidislice, 1, slicetimeVar)
          call writestat_dims_nc(ncidislice)
          ! Add islice-specific x coordinate information
          call instant_write_islice_xcoord_local(ncidislice, local_nislice, nislice, islice, myidx, nprocx, itot)
        end if
        call define_nc(ncidislice, nslicevars, isliceVars)
        
        ! write(*,'(A,A,A,I2,A)') '  Processor (myidx=', cmyidx, ') created islice file with ', local_nislice, ' slices'
      end if
    end subroutine instant_create_ncislice

    subroutine instant_create_ncjslice
      implicit none
      integer :: j
      logical :: has_jslice

      ! Count how many jslices are assigned to this Y-processor
      local_njslice = 0
      has_jslice = .false.
      do j = 1, njslice
        ! if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
        if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
          local_njslice = local_njslice + 1
          has_jslice = .true.
        end if
      end do
      
      ! Only processors with myidx==0 create files for their myidy column
      if (has_jslice .and. myidx == 0) then
        filenamejslice = 'ins_jslice.xxx.xxx.nc'
        filenamejslice(12:14) = cmyidy   ! Y-processor ID
        filenamejslice(16:18) = cexpnr  ! Experiment number
        
        nrecjslice = 0
        ! Use itot as x-dimension (global), local_n as jslice-dimension, zdim as z
        call open_nc(filenamejslice, ncidjslice, nrecjslice, n1=itot, n2=local_njslice, n3=zdim)
        if (nrecjslice==0) then
          call define_nc(ncidjslice, 1, slicetimeVar)
          call writestat_dims_nc(ncidjslice)
          call instant_write_jslice_ycoord_local(ncidjslice, local_njslice, njslice, jslice, myidy, nprocy, jtot)
        end if
        call define_nc(ncidjslice, nslicevars, jsliceVars)
        ! write(*,'(A,A,A,I2,A)') '  Processor (myidy=', cmyidy, ') created islice file with ', local_njslice, ' slices'
      end if
    end subroutine instant_create_ncjslice

    !! ## %% 1D parallel output creation (y-direction processes only)
    subroutine instant_create_nckslice
      implicit none

      if (myidy == 0) then
        filenamekslice = 'ins_kslice.xxx.xxx.nc'
        filenamekslice(12:14) = cmyidx   ! Only x-processor ID
        filenamekslice(16:18) = cexpnr  ! Experiment number

        nreckslice = 0
        ! Use kdim (nkslice) as the z-dimension, jtot as global y-dimension
        call open_nc(filenamekslice, ncidkslice, nreckslice, n1=xdim, n2=jtot, n3=kdim)
        if (nreckslice == 0) then
          call define_nc(ncidkslice, 1, slicetimeVar)
          call writestat_dims_nc(ncidkslice)
          ! Add kslice-specific z coordinate information
          call instant_write_kslice_zcoord(ncidkslice, nkslice, kslice)
        end if
        call define_nc(ncidkslice, nslicevars, ksliceVars)
      end if
    end subroutine instant_create_nckslice

    subroutine instant_write_islice
      implicit none
      real, allocatable :: tmp_slice(:,:,:)
      integer :: i, ii, ii_local, local_idx
      
      if (nislice == 0) return
      if (local_nislice == 0) return  ! No islices on this processor
      
      ! Allocate temporary array with local_nislice (module variable)
      allocate(tmp_slice(local_nislice, jb:je, kb:ke))
      tmp_slice = 0.0
      
      ! Write time variable (only myidy==0 writes)
      if (myidy == 0) then
        call writestat_nc(ncidislice, 'time', timee, nrecislice, .true.)
      end if
      
      ! u velocity (interpolated to cell centers in x-direction)
      if (present('u0')) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          ! if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)  ! Global X-position
            ! ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = 0.5*(um(ii_local,jb:je,kb:ke) + um(ii_local+1,jb:je,kb:ke))
          end if
        end do
        call writeoffset(ncidislice, 'u', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if

      ! v velocity
      if (present('v0')) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          ! if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
            ! ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = vm(ii_local,jb:je,kb:ke)
          end if
        end do
        call writeoffset(ncidislice, 'v', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      ! w velocity
      if (present('w0')) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          ! if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
            ! ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = wm(ii_local,jb:je,kb:ke)
          end if
        end do
        call writeoffset(ncidislice, 'w', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      ! p pressure
      if (present('p0')) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
          ! if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
            ! ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = pres0(ii_local,jb:je,kb:ke)
          end if
        end do
        call writeoffset(ncidislice, 'p', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if

      ! temperature
      if (present('th') .and. ltempeq) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
!          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
!            ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = thlm(ii_local,jb:je,kb:ke)
          end if
        end do
        call writeoffset(ncidislice, 'thl', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      ! moisture
      if (present('qt') .and. lmoist) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
!          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
!            ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = qtm(ii_local,jb:je,kb:ke)
          end if
        end do
        call writeoffset(ncidislice, 'qt', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      ! scalars s1-s4
      if (present('s1') .and. nsv>0) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
!          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
!            ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = svm(ii_local,jb:je,kb:ke,1)
          end if
        end do
        call writeoffset(ncidislice, 's1', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      if (present('s2') .and. nsv>1) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
!          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
!            ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = svm(ii_local,jb:je,kb:ke,2)
          end if
        end do
        call writeoffset(ncidislice, 's2', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      if (present('s3') .and. nsv>2) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
!          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
!            ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = svm(ii_local,jb:je,kb:ke,3)
          end if
        end do
        call writeoffset(ncidislice, 's3', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      if (present('s4') .and. nsv>3) then
        tmp_slice = 0.0
        local_idx = 0
        do i = 1, nislice
!          if ( (islice(i)-1)/(itot/nprocx) == myidx) then
          if (islice(i) >= zstart(1) .and. islice(i) <= zend(1)) then
            local_idx = local_idx + 1
            ii = islice(i)
!            ii_local = ii - (myidx * (itot/nprocx))
            ii_local = ii - zstart(1) + 1
              tmp_slice(local_idx,:,:) = svm(ii_local,jb:je,kb:ke,4)
          end if
        end do
        call writeoffset(ncidislice, 's4', tmp_slice, nrecislice, local_nislice, ydim, zdim)
      end if
      
      deallocate(tmp_slice)
    end subroutine instant_write_islice

    subroutine instant_write_jslice
      implicit none
      real, allocatable :: tmp_slice(:,:,:)
      integer :: j, jj, jj_local, local_idy
      
      if (njslice == 0) return
      if (local_njslice == 0) return

      allocate(tmp_slice(ib:ie, local_njslice, kb:ke))
      tmp_slice = 0.0

      ! Write time variable (only myidx==0 writes)
      if (myidx == 0) then
        call writestat_nc(ncidjslice, 'time', timee, nrecjslice, .true.)
      end if

      ! u velocity
      if (present('u0')) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = um(ib:ie, jj_local, kb:ke)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'u', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      ! v velocity (interpolated to cell centers in y-direction)
      if (present('v0')) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = 0.5*(vm(ib:ie, jj_local, kb:ke) + vm(ib:ie, jj_local+1, kb:ke))
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'v', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      ! w velocity
      if (present('w0')) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = wm(ib:ie, jj_local, kb:ke)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'w', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      ! p pressure
      if (present('p0')) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = pres0(ib:ie, jj_local, kb:ke)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'p', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      ! temperature
      if (present('th') .and. ltempeq) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = thlm(ib:ie, jj_local, kb:ke)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'thl', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      ! moisture
      if (present('qt') .and. lmoist) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = qtm(ib:ie, jj_local, kb:ke)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 'qt', tmp_slice, nrecjslice, xdim, local_njslice, zdim)
      end if

      ! scalars s1-s4
      if (present('s1') .and. nsv>0) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = svm(ib:ie, jj_local, kb:ke, 1)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 's1', tmp_slice, nrecjslice, xdim,  local_njslice, zdim)
      end if

      if (present('s2') .and. nsv>1) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = svm(ib:ie, jj_local, kb:ke, 2)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 's2', tmp_slice, nrecjslice, xdim,  local_njslice, zdim)
      end if

      if (present('s3') .and. nsv>2) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = svm(ib:ie, jj_local, kb:ke, 3)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 's3', tmp_slice, nrecjslice, xdim,  local_njslice, zdim)
      end if

      if (present('s4') .and. nsv>3) then
        tmp_slice = 0.0
        local_idy = 0
        do j = 1, njslice
!          if (jslice(j) >= ystart(2) .and. jslice(j) <= yend(2)) then
          if (jslice(j) >= zstart(2) .and. jslice(j) <= zend(2)) then
            local_idy = local_idy + 1
            jj = jslice(j)
!            jj_local = jj - ystart(2) + 1
            jj_local = jj - zstart(2) + 1
              tmp_slice(:, local_idy, :) = svm(ib:ie, jj_local, kb:ke, 4)
          end if
        end do
        call writeoffset_1dx(ncidjslice, 's4', tmp_slice, nrecjslice, xdim,  local_njslice, zdim)
      end if

      deallocate(tmp_slice)
    end subroutine instant_write_jslice

    subroutine instant_write_islice_xcoord_local(ncid, local_n, nislice_total, islice_positions, myidx_in, nprocx_in, itot_in)
      ! Write x-coordinates for LOCAL islice positions only
      implicit none
      integer, intent(in) :: ncid, local_n, nislice_total
      integer, dimension(nislice_total), intent(in) :: islice_positions
      integer, intent(in) :: myidx_in, nprocx_in, itot_in
      integer :: varid, ierr, i, local_idx
      real, allocatable :: x_islice_f(:), x_islice_h(:)
      integer, allocatable :: islice_indices(:)
      
      ! Allocate and fill x coordinates for LOCAL islices only
      allocate(x_islice_f(local_n))
      allocate(x_islice_h(local_n))
      allocate(islice_indices(local_n))
      
      local_idx = 0
      do i = 1, nislice_total
!        if ( (islice_positions(i)-1)/(itot_in/nprocx_in) == myidx_in) then
        if (islice_positions(i) >= zstart(1) .and. islice_positions(i) <= zend(1)) then
          local_idx = local_idx + 1
          x_islice_f(local_idx) = xf(islice_positions(i))  ! full level (cell center)
          x_islice_h(local_idx) = xh(islice_positions(i))  ! half level (cell edge)
          islice_indices(local_idx) = islice_positions(i)
        end if
      end do
      
      ! Write xt (full level)
      ierr = nf90_inq_varid(ncid, 'xt', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'islice_indices', islice_indices)
        ierr = nf90_put_att(ncid, varid, 'long_name', 'x-coordinate of i-slices (cell center)')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, x_islice_f)
      end if
      
      ! Write xm (half level)
      ierr = nf90_inq_varid(ncid, 'xm', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'islice_indices', islice_indices)
        ierr = nf90_put_att(ncid, varid, 'long_name', 'x-coordinate of i-slices (cell edge)')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, x_islice_h)
      end if
      
      deallocate(x_islice_f)
      deallocate(x_islice_h)
      deallocate(islice_indices)
    end subroutine instant_write_islice_xcoord_local

    subroutine instant_write_jslice_ycoord_local(ncid, local_n, njslice_total, jslice_positions, myidy_in, nprocy_in, jtot_in)
      ! Write y-coordinates for LOCAL jslice positions (round-robin assignment over myidy)
      implicit none
      integer, intent(in) :: ncid, local_n, njslice_total
      integer, dimension(njslice_total), intent(in) :: jslice_positions
      integer, intent(in) :: myidy_in, nprocy_in, jtot_in
      integer :: varid, ierr, j, local_idy
      real, allocatable :: y_jslice_f(:), y_jslice_h(:)
      integer, allocatable :: jslice_indices(:)

      allocate(y_jslice_f(local_n))
      allocate(y_jslice_h(local_n))
      allocate(jslice_indices(local_n))
      
      local_idy = 0
      do j = 1, njslice_total
!        if (jslice_positions(j) >= ystart(2) .and. jslice_positions(j) <= yend(2)) then
        if (jslice_positions(j) >= zstart(2) .and. jslice_positions(j) <= zend(2)) then
          local_idy = local_idy + 1
          y_jslice_f(local_idy) = yf(jslice_positions(j))
          y_jslice_h(local_idy) = yh(jslice_positions(j))
          jslice_indices(local_idy) = jslice_positions(j)
        end if
      end do

      ierr = nf90_inq_varid(ncid, 'yt', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'jslice_indices', jslice_indices)
        ierr = nf90_put_att(ncid, varid, 'long_name', 'y-coordinate of j-slices (cell center)')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, y_jslice_f)
      end if

      ierr = nf90_inq_varid(ncid, 'ym', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'jslice_indices', jslice_indices)
        ierr = nf90_put_att(ncid, varid, 'long_name', 'y-coordinate of j-slices (cell edge)')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, y_jslice_h)
      end if

      deallocate(y_jslice_f)
      deallocate(y_jslice_h)
      deallocate(jslice_indices)
    end subroutine instant_write_jslice_ycoord_local

    subroutine instant_write_kslice_zcoord(ncid, nkslice_count, kslice_positions)
      ! Update z-coordinate to reflect kslice levels instead of full z levels
      implicit none
      integer, intent(in) :: ncid, nkslice_count
      integer, dimension(nkslice_count), intent(in) :: kslice_positions
      integer :: varid, ierr, k
      real, allocatable :: z_kslice_f(:), z_kslice_h(:)
      
      ! Allocate and fill z coordinates for kslices
      allocate(z_kslice_f(nkslice_count))
      allocate(z_kslice_h(nkslice_count))
      do k = 1, nkslice_count
        z_kslice_f(k) = zf(kslice_positions(k))  ! full level (cell center)
        z_kslice_h(k) = zh(kslice_positions(k))  ! half level (cell edge)
      end do
      
      ! Write zt (full level)
      ierr = nf90_inq_varid(ncid, 'zt', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'kslice_indices', kslice_positions)
        ierr = nf90_put_att(ncid, varid, 'long_name', 'z-coordinate of k-slices (cell center)')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, z_kslice_f)
      end if
      
      ! Write zm (half level)
      ierr = nf90_inq_varid(ncid, 'zm', varid)
      if (ierr == nf90_noerr) then
        ierr = nf90_redef(ncid)
        ierr = nf90_put_att(ncid, varid, 'kslice_indices', kslice_positions)
        ierr = nf90_put_att(ncid, varid, 'long_name', 'z-coordinate of k-slices (cell edge)')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, z_kslice_h)
      end if
      
      deallocate(z_kslice_f)
      deallocate(z_kslice_h)
    end subroutine instant_write_kslice_zcoord

    subroutine instant_write_kslice
      implicit none
      real, allocatable :: tmp_slice(:,:,:)
      integer :: k, kk
      
      if (nkslice == 0) return
      
      ! Allocate temporary array to hold all kslice levels: (x, y, nkslice)
      allocate(tmp_slice(ib:ie, jb:je, nkslice))
      
      ! Write time variable (only myidy==0 writes)
      if (myidy == 0) then
        call writestat_nc(ncidkslice, 'time', timee, nreckslice, .true.)
      end if
      
      ! u velocity
      if (present('u0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = um(ib:ie, jb:je, kk)
        end do
        call writeoffset(ncidkslice, 'u', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      ! v velocity
      if (present('v0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = vm(ib:ie, jb:je, kk)
        end do
        call writeoffset(ncidkslice, 'v', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      ! w velocity (interpolated to cell centers)
      if (present('w0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = 0.5 * (wm(ib:ie, jb:je, kk) + wm(ib:ie, jb:je, kk+1))
        end do
        call writeoffset(ncidkslice, 'w', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if

      ! pressure
      if (present('p0')) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = pres0(ib:ie, jb:je, kk)
        end do
        call writeoffset(ncidkslice, 'p', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if


      ! temperature
      if (present('th') .and. ltempeq) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = thlm(ib:ie, jb:je, kk)
        end do
        call writeoffset(ncidkslice, 'thl', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      ! moisture
      if (present('qt') .and. lmoist) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = qtm(ib:ie, jb:je, kk)
        end do
        call writeoffset(ncidkslice, 'qt', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      ! scalars s1-s4
      if (present('s1') .and. nsv>0) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie, jb:je, kk, 1)
        end do
        call writeoffset(ncidkslice, 's1', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      if (present('s2') .and. nsv>1) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie, jb:je, kk, 2)
        end do
        call writeoffset(ncidkslice, 's2', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      if (present('s3') .and. nsv>2) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie, jb:je, kk, 3)
        end do
        call writeoffset(ncidkslice, 's3', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      if (present('s4') .and. nsv>3) then
        do k = 1, nkslice
          kk = kslice(k)
          tmp_slice(:,:,k) = svm(ib:ie, jb:je, kk, 4)
        end do
        call writeoffset(ncidkslice, 's4', tmp_slice, nreckslice, xdim, ydim, kdim)
      end if
      
      deallocate(tmp_slice)
    end subroutine instant_write_kslice

    logical function present(str)
      implicit none
      character(len=*), intent(in) :: str
      integer :: pos

      pos = index(slicevars, str)   ! returns the position of the substring str in slicevars. If not found, it returns 0.
      present = (pos > 0)
    end function present


    subroutine instant_probe_init
      use modglobal, only : ifinput
      implicit none
      integer :: n, ierr, vn
      integer :: point_dimid, time_dimid
      integer :: varid_xt, varid_xm, varid_yt, varid_ym, varid_zt, varid_zm
      real,    allocatable :: xt(:), xm(:), yt(:), ym(:), zt(:), zm(:)
      character(2)  :: varname
      character(80) :: chmess
      
      if (.not. lprobedump) return

      ! Validate probe inputs before proceeding
      call instant_validate_out_time('probes')
      call instant_validate_output_vars(nprobevars, probevars, 'probevars', 'lprobedump')

      if (nprobe <= 0) then
        write(0, *) 'ERROR: lprobedump=.true. but nprobe=', nprobe, ' (must be > 0)'
        stop 1
      end if

      ! read global probe point indices
      if(myid==0) then
        open (ifinput,file='probe.inp.'//cexpnr)
          read (ifinput,'(a80)') chmess
          do n = 1, nprobe
            read (ifinput,*) iprobe(n), jprobe(n), kprobe(n)
          end do
        close (ifinput)
      end if
      ! Broadcast the probe point indices to all processes
      call MPI_BCAST(iprobe, nprobe, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(jprobe, nprobe, MPI_INTEGER, 0, comm3d, mpierr)
      call MPI_BCAST(kprobe, nprobe, MPI_INTEGER, 0, comm3d, mpierr)
      
      do n = 1, nprobe
        if (iprobe(n) < 1 .or. iprobe(n) > itot) then
          write(0, *) 'ERROR: iprobe(', n, ')=', iprobe(n), ' is out of bounds (1 to ', itot, ')'
          stop 1
        end if
        if (jprobe(n) < 1 .or. jprobe(n) > jtot) then
          write(0, *) 'ERROR: jprobe(', n, ')=', jprobe(n), ' is out of bounds (1 to ', jtot, ')'
          stop 1
        end if
        if (kprobe(n) < 1 .or. kprobe(n) > ktot) then
          write(0, *) 'ERROR: kprobe(', n, ')=', kprobe(n), ' is out of bounds (1 to ', ktot, ')'
          stop 1
        end if
      end do
      
      if (myid == 0) then
          allocate(varid_vals(nprobevars))
          
          ! write(*,*) "=== Probe Output Init (Gather to Rank 0) ==="
          ! write(*,*) "Number of probes:", nprobe
          ! write(*,*) "Variables:", trim(probevars)
          
          filenameprobe = 'ins_probe.xxx.nc'
          filenameprobe(11:13) = cexpnr
          
          ! Create NetCDF file
          call check( nf90_create(trim(filenameprobe), NF90_CLOBBER, ncidprobe) )
          
          ! Define Dimensions
          call check( nf90_def_dim(ncidprobe, 'point', nprobe, point_dimid) )
          call check( nf90_def_dim(ncidprobe, 'time', NF90_UNLIMITED, time_dimid) )
          
          ! Define Variables Dimensions
          ! Coordinates (cell center and cell edge)
          call check( nf90_def_var(ncidprobe, 'xt', NF90_FLOAT, (/point_dimid/), varid_xt) )
          call check( nf90_put_att(ncidprobe, varid_xt, 'long_name', 'x-coordinate (cell center)') )
          call check( nf90_put_att(ncidprobe, varid_xt, 'units', 'm') )
          
          ! call check( nf90_def_var(ncidprobe, 'xm', NF90_FLOAT, (/point_dimid/), varid_xm) )
          ! call check( nf90_put_att(ncidprobe, varid_xm, 'long_name', 'x-coordinate (cell edge)') )
          ! call check( nf90_put_att(ncidprobe, varid_xm, 'units', 'm') )
          
          call check( nf90_def_var(ncidprobe, 'yt', NF90_FLOAT, (/point_dimid/), varid_yt) )
          call check( nf90_put_att(ncidprobe, varid_yt, 'long_name', 'y-coordinate (cell center)') )
          call check( nf90_put_att(ncidprobe, varid_yt, 'units', 'm') )
          
          ! call check( nf90_def_var(ncidprobe, 'ym', NF90_FLOAT, (/point_dimid/), varid_ym) )
          ! call check( nf90_put_att(ncidprobe, varid_ym, 'long_name', 'y-coordinate (cell edge)') )
          ! call check( nf90_put_att(ncidprobe, varid_ym, 'units', 'm') )

          call check( nf90_def_var(ncidprobe, 'zt', NF90_FLOAT, (/point_dimid/), varid_zt) )
          call check( nf90_put_att(ncidprobe, varid_zt, 'long_name', 'z-coordinate (cell center)') )
          call check( nf90_put_att(ncidprobe, varid_zt, 'units', 'm') )
          
          ! call check( nf90_def_var(ncidprobe, 'zm', NF90_FLOAT, (/point_dimid/), varid_zm) )
          ! call check( nf90_put_att(ncidprobe, varid_zm, 'long_name', 'z-coordinate (cell edge)') )
          ! call check( nf90_put_att(ncidprobe, varid_zm, 'units', 'm') )

          ! Time
          call check( nf90_def_var(ncidprobe, 'time', NF90_FLOAT, (/time_dimid/), varid_time) )
          call check( nf90_put_att(ncidprobe, varid_time, 'long_name', 'Time') )
          call check( nf90_put_att(ncidprobe, varid_time, 'units', 's') )

          ! Data Variables
          do vn = 1, nprobevars
             varname = probevars(3*vn-2 : 3*vn-1)
             call check( nf90_def_var(ncidprobe, trim(instant_get_nc_varname(varname)), NF90_FLOAT, &
                                      (/point_dimid, time_dimid/), varid_vals(vn)) )
             call instant_add_var_atts(ncidprobe, varid_vals(vn), varname)
          end do
          
          call check( nf90_enddef(ncidprobe) )
          
          ! Write Static Coordinates immediately
          allocate(xt(nprobe), xm(nprobe))
          allocate(yt(nprobe), ym(nprobe))
          allocate(zt(nprobe), zm(nprobe))
          do n=1, nprobe
             xt(n) = xf(iprobe(n))  ! cell center
            !  xm(n) = xh(iprobe(n))  ! cell edge
             yt(n) = yf(jprobe(n))  ! cell center
            !  ym(n) = yh(jprobe(n))  ! cell edge
             zt(n) = zf(kprobe(n))  ! cell center
            !  zm(n) = zh(kprobe(n))  ! cell edge
          end do
          
          call check( nf90_put_var(ncidprobe, varid_xt, xt, start=(/1/), count=(/nprobe/)) )
          ! call check( nf90_put_var(ncidprobe, varid_xm, xm, start=(/1/), count=(/nprobe/)) )
          call check( nf90_put_var(ncidprobe, varid_yt, yt, start=(/1/), count=(/nprobe/)) )
          ! call check( nf90_put_var(ncidprobe, varid_ym, ym, start=(/1/), count=(/nprobe/)) )
          call check( nf90_put_var(ncidprobe, varid_zt, zt, start=(/1/), count=(/nprobe/)) )
          ! call check( nf90_put_var(ncidprobe, varid_zm, zm, start=(/1/), count=(/nprobe/)) )
          call check( nf90_sync(ncidprobe) )
          
          deallocate(xt, xm, yt, ym, zt, zm)
      end if
      
    end subroutine instant_probe_init


    subroutine instant_probe_main
      implicit none
      real, allocatable :: send_buf(:), recv_buf(:)
      integer :: vn
      character(80) :: varname
      
      if (.not. lprobedump) return
         
      ! Only rank 0 tracks record number
      if (myid == 0) nrecprobe = nrecprobe + 1
      
      allocate(send_buf(nprobe))
      allocate(recv_buf(nprobe))
      
      ! Always write time first (Rank 0)
      if (myid == 0) then
          call check( nf90_put_var(ncidprobe, varid_time, timee, start=(/nrecprobe/)) )
      end if

      ! Loop over all requested variables
      do vn = 1, nprobevars
        varname = probevars(3*vn-2 : 3*vn-1)
        send_buf = 0.0
        recv_buf = 0.0
        
        call instant_gather_probe_var(varname, send_buf)
        
        ! Reduce to Rank 0
        call MPI_REDUCE(send_buf, recv_buf, nprobe, my_real, MPI_SUM, 0, comm3d, mpierr)
        
        if (myid == 0) then
            ! Use start/count to append to time dimension
            call check( nf90_put_var(ncidprobe, varid_vals(vn), recv_buf, &
                                    start=(/1, nrecprobe/), count=(/nprobe, 1/)) )
        end if
      end do
      
      if (myid == 0) then
        call check( nf90_sync(ncidprobe) )
      end if
      
      deallocate(send_buf, recv_buf)
    end subroutine instant_probe_main


    subroutine instant_gather_probe_var(vname, buf)
      implicit none
      character(len=*), intent(in) :: vname
      real, intent(out) :: buf(:)
      integer :: n, gi, gj, gk, li, lj, lk
      
      do n = 1, nprobe
         gi = iprobe(n)
         gj = jprobe(n)
         gk = kprobe(n)
         
         ! Check if local using decomp_2d zstart/zend
         if (gi >= zstart(1) .and. gi <= zend(1) .and. &
             gj >= zstart(2) .and. gj <= zend(2) .and. &
             gk >= zstart(3) .and. gk <= zend(3)) then
             
             li = gi - zstart(1) + 1
             lj = gj - zstart(2) + 1
             lk = gk - zstart(3) + 1
             
             select case (vname)
               case ('u0')
                  buf(n) = 0.5 * (um(li,lj,lk) + um(li+1,lj,lk))
               case ('v0')
                  buf(n) = 0.5 * (vm(li,lj,lk) + vm(li,lj+1,lk))
               case ('w0')
                  buf(n) = 0.5 * (wm(li,lj,lk) + wm(li,lj,lk+1))
               case ('p0')
                  buf(n) = pres0(li,lj,lk)
               case ('th')
                  if (ltempeq) buf(n) = thlm(li,lj,lk)
               case ('qt')
                  if (lmoist) buf(n) = qtm(li,lj,lk)
               case ('s1')
                  if (nsv>0) buf(n) = svm(li,lj,lk,1)
               case ('s2')
                  if (nsv>1) buf(n) = svm(li,lj,lk,2)
               case ('s3')
                  if (nsv>2) buf(n) = svm(li,lj,lk,3)
               case ('s4')
                  if (nsv>3) buf(n) = svm(li,lj,lk,4)
             end select
         endif
      end do
    end subroutine instant_gather_probe_var

    function instant_get_nc_varname(vname) result(ncname)
       implicit none
       character(len=*), intent(in) :: vname
       character(len=10) :: ncname
       select case(vname)
         case('u0')
            ncname = 'u'
         case('v0')
            ncname = 'v'
         case('w0')
            ncname = 'w'
         case('p0')
            ncname = 'p'
         case('th')
            ncname = 'thl'
         case('qt')
            ncname = 'qt'
         case('s1','s2','s3','s4')
            ncname = vname
         case default
            ncname = 'unknown'
       end select
    end function
    
    subroutine instant_add_var_atts(ncid, vid, vname)
       implicit none
       integer, intent(in) :: ncid, vid
       character(len=*), intent(in) :: vname
       integer :: ierr
       select case(vname)
         case('u0')
            ierr = nf90_put_att(ncid, vid, 'long_name', 'Streamwise velocity')
            ierr = nf90_put_att(ncid, vid, 'units', 'm/s')
         case('v0')
            ierr = nf90_put_att(ncid, vid, 'long_name', 'Spanwise velocity')
            ierr = nf90_put_att(ncid, vid, 'units', 'm/s')
         case('w0')
             ierr = nf90_put_att(ncid, vid, 'long_name', 'Vertical velocity')
             ierr = nf90_put_att(ncid, vid, 'units', 'm/s')
         case('p0')
             ierr = nf90_put_att(ncid, vid, 'long_name', 'Kinematic pressure')
             ierr = nf90_put_att(ncid, vid, 'units', 'm^2/s^2')
         case('th')
             ierr = nf90_put_att(ncid, vid, 'long_name', 'Potential temperature')
             ierr = nf90_put_att(ncid, vid, 'units', 'K')
         case('qt')
             ierr = nf90_put_att(ncid, vid, 'long_name', 'Specific humidity')
             ierr = nf90_put_att(ncid, vid, 'units', 'kg/kg')
         case('s1','s2','s3','s4')
             ierr = nf90_put_att(ncid, vid, 'long_name', 'Scalar concentration')
             ierr = nf90_put_att(ncid, vid, 'units', 'g/m^3')
       end select
    end subroutine instant_add_var_atts

    subroutine check(status)
      implicit none
      integer, intent(in) :: status
      if(status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        stop 2
      end if
    end subroutine check


    subroutine instant_validate_out_time(var_name)
      implicit none
      character(len=*), intent(in)  :: var_name
      if(runtime <= tstatstart) then
        if(myid==0) then
          write(*,*) "ERROR: no instantaneous ", trim(var_name), " will be written as runtime <= tstatstart. Note that runtime &
                      &must be greater than tstatstart for wiriting ", trim(var_name), " files."
          write(*,*) "You have used runtime = ", runtime, ", tstatstart = ", tstatstart
          write(*,*) "Either correct the time settings or change all the ", trim(var_name), " writing flags to false."
          stop 1
        end if
      end if
    end subroutine instant_validate_out_time

    subroutine instant_validate_output_vars(nvars, vars, var_name, var_flag)
      ! Count variables from vars string; stop with error if none found
      implicit none
      integer,          intent(out) :: nvars
      character(len=*), intent(in)  :: vars
      character(len=*), intent(in)  :: var_name
      character(len=*), intent(in)  :: var_flag

      nvars = (LEN(trim(vars))+1)/3
      if (nvars == 0) then
        if(myid==0) then
          print *, "ERROR: no '", trim(var_name), "' entries found, although ", trim(var_flag), " are set to true !!"
          stop 1
        end if
      end if
    end subroutine instant_validate_output_vars

    subroutine instant_validate_slice_inputs(slice_enabled, n_slice, slice_array, slice_name)
      ! Validate slice arrays: check count and non-zero entries
      implicit none
      logical, intent(in) :: slice_enabled
      integer, intent(in) :: n_slice
      integer, intent(in) :: slice_array(:)
      character(len=*), intent(in) :: slice_name
      
      if (slice_enabled) then
         if (n_slice <= 0) then
            write(0, *) 'ERROR: l', trim(slice_name), 'dump=.true. but n', trim(slice_name), '=', n_slice, ' (must be > 0)'
            stop 1
         end if
         if (count(slice_array(1:n_slice) > 0) /= n_slice) then
            write(0, *) 'ERROR: n', trim(slice_name), '=', n_slice, ' but only', count(slice_array(1:n_slice) > 0), &
                         ' valid (>0) entries found in ', trim(slice_name), ' array'
            write(0, *) 'Check that ', trim(slice_name), ' has exactly n', trim(slice_name), ' positive values'
            stop 1
         end if
      end if
    end subroutine instant_validate_slice_inputs

end module instant