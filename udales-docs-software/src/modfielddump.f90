!> \file modfielddump.f90
!!  Dumps 3D fields of several variables
!>
!!  Dumps 3D fields of several variables Written to wb*.myid.expnr
!! If netcdf is true, this module leads the fielddump.myid.expnr.nc output
!!  \author Jasper Tomas, TU Delft Match 31 2014
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!! Possibility to write tecplot (formatted!!) file
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
module modfielddump

  use mpi
  use modglobal, only : longint
  use modfields, only : ncname, ncname1, ncname2
  implicit none
  private
  PUBLIC :: initfielddump, fielddump,exitfielddump
  save
  !NetCDF variables
  integer :: ncid,nrec = 0
!  real, pointer :: point
  type domainptr
    real, pointer :: point(:,:,:)
  end type domainptr
  type(domainptr), dimension(30) :: pfields

  ! Fields that no model array holds. They are computed into dbuf at dump time,
  ! only for the codes listed in fieldvars; nothing is stored between dumps (#330).
  integer, parameter :: DERIV_DIV  = 1, DERIV_DUDX = 2, DERIV_DVDY = 3, DERIV_DWDZ = 4, &
                        DERIV_PRHS = 5, DERIV_PDX  = 6, DERIV_PDY  = 7, DERIV_PDZ  = 8
  integer :: derived_kind(30) = 0   !< per output field: 0 = points at a live model array
  integer :: derived_slot(30) = 0   !< per output field: its slice of dbuf
  integer :: nderived = 0
  real, allocatable, target :: dbuf(:,:,:,:)

  character(80) :: fname = 'fielddump.xxx.xxx.xxx.nc'
  !dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname

  integer :: ilow,ihigh,jlow,jhigh,klow,khigh,nvar
  logical :: ldiracc   = .false. !< switch for doing direct access writing (on/off)
  logical :: lbinary   = .false. !
  logical :: lhalos

contains
 !> Initializing fielddump. Read out the namelist, initializing the variables
  subroutine initfielddump
    use modmpi,   only   :mpierr,comm3d,mpi_logical,mpi_integer,cmyidx,cmyidy,mpi_character
    use modglobal,only   :cexpnr,kb,ke,fieldvars,ib,ie,jb,je,kb,ke, ih,jh,lfielddump,kh
    use modstat_nc,only  : open_nc, define_nc,ncinfo,writestat_dims_nc
    use modfields, only  : u0,v0,w0,thl0,sv0,ql0,qt0,pres0,tau_x, tau_y, tau_z, thl_flux
    use modpois,   only  : p, pup, pvp, pwp
    use modibm, only : mask_u, mask_v, mask_w, mask_c
    implicit none
    integer :: ierr, n

    !type(domainptr), dimension(nvar) :: pfields

    nvar = (LEN(trim(fieldvars))+1)/3

    if (nvar == 0) then
      lfielddump = .false.
      print *, 'empty fieldvars therefore lfielddump = .false. and no instantaneous fields outputted'
      return
    else
      allocate(ncname(nvar,4))
    end if




    lhalos = .false.

    if (lhalos) then
      ilow = ib-ih
      ihigh = ie+ih
      jlow = jb-jh
      jhigh = je+jh
      klow=kb-kh
      khigh=ke+kh
    else
      ilow = ib
      ihigh = ie
      jlow = jb
      jhigh = je
      klow=kb
      khigh=ke
    end if


!ils13 13.08.18: why is this broadcast, doesn't every processor do it anyway?
    call MPI_BCAST(klow        ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(khigh       ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(lfielddump  ,1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(ldiracc     ,1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(lbinary     ,1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(ncname     ,80,MPI_CHARACTER,0,comm3d,mpierr)
    call MPI_BCAST(nvar       ,1,MPI_INTEGER,0,comm3d,mpierr)

    !    dt_lim = min(dt_lim,tnext)

    if(.not.(lfielddump)) return

    fname(11:13) = cmyidx
    fname(15:17) = cmyidy
    fname(19:21) = cexpnr
    call ncinfo(tncname(1,:),'time','Time','s','time')

    ! tg3315 reads in fields specified by fieldvars
    if (lhalos) then

    do n=1,nvar
      select case(fieldvars(3*n-2:3*n-1))
      case('u0')
        call ncinfo(ncname( n,:),'u','West-East velocity','m/s','mttt')
        pfields(n)%point => u0(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      case('v0')
        call ncinfo(ncname( n,:),'v','South-North velocity','m/s','tmtt')
        pfields(n)%point => v0(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      case('w0')
        call ncinfo(ncname( n,:),'w','Vertical velocity','m/s','ttmt')
        pfields(n)%point => w0(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      case('th')
        call ncinfo(ncname( n,:),'thl','Liquid water potential temperature','K','tttt')
        pfields(n)%point => thl0(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      case('ql')
        call ncinfo(ncname( n,:),'ql','Liquid water mixing ratio','1e-5kg/kg','tttt')
        pfields(n)%point => ql0(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      case('qt')
        call ncinfo(ncname( n,:),'qt','Total water mixing ratio','1e-5kg/kg','tttt')
        pfields(n)%point => qt0(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      case('s1')
        call ncinfo(ncname( n,:),'sca1','scalar 1','M','tttt')
        pfields(n)%point => sv0(ib-ih:ie+ih,jb-jh:je+jh,kb:ke,1)
      case('s2')
        call ncinfo(ncname( n,:),'sca2','scalar 2','M','tttt')
        pfields(n)%point => sv0(ib-ih:ie+ih,jb-jh:je+jh,kb:ke,2)
      case('s3')
        call ncinfo(ncname( n,:),'sca3','scalar 3','M','tttt')
        pfields(n)%point => sv0(ib-ih:ie+ih,jb-jh:je+jh,kb:ke,3)
      case('s4')
        call ncinfo(ncname( n,:),'sca4','scalar 4','M','tttt')
        pfields(n)%point => sv0(ib-ih:ie+ih,jb-jh:je+jh,kb:ke,4)
      case('s5')
        call ncinfo(ncname( n,:),'sca5','scalar 5','M','tttt')
        pfields(n)%point => sv0(ib-ih:ie+ih,jb-jh:je+jh,kb:ke,5)
      case('p0')
        call ncinfo(ncname( n,:),'pres','pressure field','M','tttt')
        pfields(n)%point => pres0(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      ! case('pd')
      !   call ncinfo(ncname( n,:),'p','pressure correction','M','tttt')
      !   pfields(n)%point => p(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      ! case('pu')
      !   call ncinfo(ncname( n,:),'pup','predicted u','M','mttt')
      !   pfields(n)%point => pup(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      ! case('pv')
      !   call ncinfo(ncname( n,:),'pvp','predicted v','M','tmtt')
      !   pfields(n)%point => pvp(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      ! case('pw')
      !   call ncinfo(ncname( n,:),'pwp','predicted w','M','ttmt')
      !   pfields(n)%point => pwp(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      ! case('du')
      !   call ncinfo(ncname( n,:),'dpupdx','','M','tttt')
      !   pfields(n)%point => dpupdx(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      ! case('dv')
      !   call ncinfo(ncname( n,:),'dpvpdy','','M','tttt')
      !   pfields(n)%point => dpvpdy(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      ! case('dw')
      !   call ncinfo(ncname( n,:),'dpwpdz','','M','tttt')
      !   pfields(n)%point => dpwpdz(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      ! case default
      !   call ncinfo(ncname( n,:),'u','West-East velocity','m/s','mttt')
      !   pfields(n)%point => u0(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh)
      end select
    end do

    else

      do n=1,nvar
        select case(fieldvars(3*n-2:3*n-1))
        case('u0')
          call ncinfo(ncname( n,:),'u','West-East velocity','m/s','mttt')
          pfields(n)%point => u0(ib:ie,jb:je,kb:ke)
        case('v0')
          call ncinfo(ncname( n,:),'v','South-North velocity','m/s','tmtt')
          pfields(n)%point => v0(ib:ie,jb:je,kb:ke)
        case('w0')
          call ncinfo(ncname( n,:),'w','Vertical velocity','m/s','ttmt')
          pfields(n)%point => w0(ib:ie,jb:je,kb:ke)
        case('th')
          call ncinfo(ncname( n,:),'thl','Liquid water potential temperature','K','tttt')
          pfields(n)%point => thl0(ib:ie,jb:je,kb:ke)
        case('ql')
          call ncinfo(ncname( n,:),'ql','Liquid water mixing ratio','1e-5kg/kg','tttt')
          pfields(n)%point => ql0(ib:ie,jb:je,kb:ke)
        case('qt')
          call ncinfo(ncname( n,:),'qt','Total water mixing ratio','1e-5kg/kg','tttt')
          pfields(n)%point => qt0(ib:ie,jb:je,kb:ke)
        case('s1')
          call ncinfo(ncname( n,:),'sca1','scalar 1','M','tttt')
          pfields(n)%point => sv0(ib:ie,jb:je,kb:ke,1)
        case('s2')
          call ncinfo(ncname( n,:),'sca2','scalar 2','M','tttt')
          pfields(n)%point => sv0(ib:ie,jb:je,kb:ke,2)
        case('s3')
          call ncinfo(ncname( n,:),'sca3','scalar 3','M','tttt')
          pfields(n)%point => sv0(ib:ie,jb:je,kb:ke,3)
        case('s4')
          call ncinfo(ncname( n,:),'sca4','scalar 4','M','tttt')
          pfields(n)%point => sv0(ib:ie,jb:je,kb:ke,4)
        case('s5')
          call ncinfo(ncname( n,:),'sca5','scalar 5','M','tttt')
          pfields(n)%point => sv0(ib:ie,jb:je,kb:ke,5)
        case('p0')
          call ncinfo(ncname( n,:),'pres','pressure field','M','tttt')
          pfields(n)%point => pres0(ib:ie,jb:je,kb:ke)
        case('tx')
          call ncinfo(ncname( n,:),'tau_x','stress x','M','mttt')
          pfields(n)%point => tau_x(ib:ie,jb:je,kb:ke)
        case('ty')
          call ncinfo(ncname( n,:),'tau_y','stress y','M','tmtt')
          pfields(n)%point => tau_y(ib:ie,jb:je,kb:ke)
        case('tz')
          call ncinfo(ncname( n,:),'tau_z','stress z','M','ttmt')
          pfields(n)%point => tau_z(ib:ie,jb:je,kb:ke)
        case('hf')
          call ncinfo(ncname( n,:),'thl_flux','heat flux','M','tttt')
          pfields(n)%point => thl_flux(ib:ie,jb:je,kb:ke)
        case('mu')
          call ncinfo(ncname( n,:),'mask_u','mask u','M','mttt')
          pfields(n)%point => mask_u(ib:ie,jb:je,kb:ke)
        case('mv')
          call ncinfo(ncname( n,:),'mask_v','mask v','M','tmtt')
          pfields(n)%point => mask_v(ib:ie,jb:je,kb:ke)
        case('mw')
          call ncinfo(ncname( n,:),'mask_w','mask w','M','ttmt')
          pfields(n)%point => mask_w(ib:ie,jb:je,kb:ke)
        case('mc')
          call ncinfo(ncname( n,:),'mask_c','mask c','M','tttt')
          pfields(n)%point => mask_c(ib:ie,jb:je,kb:ke)
        case('pd')
          call ncinfo(ncname( n,:),'p','pressure correction of the last substep','m2/s2','tttt')
          pfields(n)%point => p(ib:ie,jb:je,kb:ke)
        case('pu')
          call ncinfo(ncname( n,:),'pup','predicted u of the last substep','m/s2','mttt')
          pfields(n)%point => pup(ib:ie,jb:je,kb:ke)
        case('pv')
          call ncinfo(ncname( n,:),'pvp','predicted v of the last substep','m/s2','tmtt')
          pfields(n)%point => pvp(ib:ie,jb:je,kb:ke)
        case('pw')
          call ncinfo(ncname( n,:),'pwp','predicted w of the last substep','m/s2','ttmt')
          pfields(n)%point => pwp(ib:ie,jb:je,kb:ke)
        ! Computed at dump time (see compute_derived); the model stores none of these.
        case('rs')
          call ncinfo(ncname( n,:),'rhs','rhs of the Poisson equation (divergence of the predicted velocity)','1/s2','tttt')
          call register_derived(n, DERIV_PRHS)
        case('du')
          call ncinfo(ncname( n,:),'dpupdx','x-term of the Poisson rhs','1/s2','tttt')
          call register_derived(n, DERIV_PDX)
        case('dv')
          call ncinfo(ncname( n,:),'dpvpdy','y-term of the Poisson rhs','1/s2','tttt')
          call register_derived(n, DERIV_PDY)
        case('dw')
          call ncinfo(ncname( n,:),'dpwpdz','z-term of the Poisson rhs','1/s2','tttt')
          call register_derived(n, DERIV_PDZ)
        case('di')
          call ncinfo(ncname( n,:),'div','Divergence after pressure correction','1/s','tttt')
          call register_derived(n, DERIV_DIV)
        case('ux')
          call ncinfo(ncname( n,:),'dudx','du/dx after pressure correction','1/s','tttt')
          call register_derived(n, DERIV_DUDX)
        case('vy')
          call ncinfo(ncname( n,:),'dvdy','dv/dy after pressure correction','1/s','tttt')
          call register_derived(n, DERIV_DVDY)
        case('wz')
          call ncinfo(ncname( n,:),'dwdz','dw/dz after pressure correction','1/s','tttt')
          call register_derived(n, DERIV_DWDZ)
        case default
          call ncinfo(ncname( n,:),'u','West-East velocity','m/s','mttt')
          pfields(n)%point => u0(ib:ie,jb:je,kb:ke)
        end select
      end do

    end if

    if (nderived > 0) then
      allocate(dbuf(ib:ie,jb:je,kb:ke,nderived)); dbuf = 0.
      do n=1,nvar
        if (derived_kind(n) /= 0) pfields(n)%point => dbuf(:,:,:,derived_slot(n))
      end do
    end if

    !call ncinfo(ncname( n,:),'u','West-East velocity','m/s','mttt')
    !call open_nc( fname, ncid, nrec, n1=imax+2, n2=jmax+2, n3=khigh-klow+1)
    call open_nc( fname, ncid, nrec, n1=ihigh-ilow+1, n2=jhigh-jlow+1, n3=khigh-klow+1)
    if (nrec==0) then
      call define_nc( ncid, 1, tncname)
      call writestat_dims_nc(ncid)
    end if

    call define_nc( ncid, nvar, ncname)
    !call open_nc( fname, ncid, nrec, n1=imax+2, n2=jmax+2, n3=khigh-klow+1)
    call open_nc( fname, ncid, nrec, n1=ihigh-ilow+1, n2=jhigh-jlow+1, n3=khigh-klow+1)

    if (nrec==0) then
      call define_nc( ncid, 1, tncname)
      call writestat_dims_nc(ncid)
    end if
    call define_nc( ncid, nvar, ncname)

     !   ! X-pencil
     !
     !   fname1(11:13) = cmyidx
     !   fname1(15:17) = cmyidy
     !   fname1(19:21) = cexpnr
     !   call ncinfo(tncname1(1,:),'time','Time','s','time')
     !   call ncinfo(ncname1(1,:),'u','West-East velocity','m/s','mttt')
     !
     !   !write(*,*) "done defining pfields"
     !
     !   call open_nc( fname1, ncid1, nrec, n1=imax1+2, n2=jmax1+2, n3=kmax1+2)
     !
     !   if (nrec==0) then
     !     call define_nc( ncid1, 1, tncname1)
     !     call writestat_dims_nc(ncid1)
     !   end if
     !   call define_nc( ncid1, nvar, ncname1)
     !   call open_nc( fname1, ncid1, nrec, n1=imax1+2, n2=jmax1+2, n3=kmax1+2)
     !   ! call open_nc( fname, ncid, nrec, n1=imax+2, n2=jmax+2, n3=khigh-klow+1)  !if want to print ghostcells
     !   if (nrec==0) then
     !     call define_nc( ncid1, 1, tncname1)
     !     call writestat_dims_nc(ncid1)
     !  end if
     !  call define_nc( ncid1, nvar, ncname1)
     !
     !
     !  ! Y-pencil
     !
     !  fname2(11:13) = cmyidx
     !  fname2(15:17) = cmyidy
     !  fname2(19:21) = cexpnr
     !  call ncinfo(tncname2(1,:),'time','Time','s','time')
     !  call ncinfo(ncname2( 1,:),'u','West-East velocity','m/s','mttt')
     !
     !  call open_nc( fname2, ncid2, nrec, n1=imax2+2, n2=jmax2+2, n3=kmax2+2)
     !  if (nrec==0) then
     !    call define_nc( ncid2, 1, tncname2)
     !    call writestat_dims_nc(ncid2)
     !  end if
     !  call define_nc( ncid2, nvar, ncname2)
     !  call open_nc( fname2, ncid2, nrec, n1=imax2+2, n2=jmax2+2, n3=kmax2+2)
     !  ! call open_nc( fname, ncid, nrec, n1=imax+2, n2=jmax+2, n3=khigh-klow+1)  !if want to print ghostcells
     !  if (nrec==0) then
     !    call define_nc( ncid2, 1, tncname2)
     !    call writestat_dims_nc(ncid2)
     ! end if
     ! call define_nc( ncid2, nvar, ncname2)

  end subroutine initfielddump

  !> Do fielddump. Collect data to truncated (2 byte) integers, and write them to file
  subroutine fielddump
    use modglobal, only : ib,ie,ih,jb,je,jh,ke,kb,kh,rk3step,timee,&
                          tfielddump, tnextfielddump, lfielddump
    !use modmpi,    only : myid,cmyid
    !use modsubgriddata, only : ekm,sbshr
    use modstat_nc, only : writestat_nc
    implicit none

    real, allocatable :: vars(:,:,:,:)
    integer n

    if (.not. ((timee>=tnextfielddump) .or. (rk3step==0))) return

    if (.not. lfielddump) return

    if (rk3step/=3 .and. rk3step/=0) return

    do n=1,nvar
      if (derived_kind(n) /= 0) call compute_derived(derived_kind(n), dbuf(:,:,:,derived_slot(n)))
    end do

    if (rk3step == 3) tnextfielddump=tnextfielddump+tfielddump

    if (lhalos) then
      allocate(vars(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh,nvar)); vars = 0;
      do n=1,nvar
        vars(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh,n) = pfields(n)%point
      end do

    else
      allocate(vars(ib:ie,jb:je,kb:ke,nvar)); vars = 0;
      do n=1,nvar
        vars(ib:ie,jb:je,kb:ke,n) = pfields(n)%point
      end do
    end if

   call writestat_nc(ncid,1,tncname,(/timee/),nrec,.true.)
   !call writestat_nc(ncid,nvar,ncname,vars,nrec,imax+2,jmax+2,khigh-klow+1)
   call writestat_nc(ncid,nvar,ncname,vars,nrec,ihigh-ilow+1,jhigh-jlow+1,khigh-klow+1)

   deallocate(vars)

  end subroutine fielddump

  !> Clean up when leaving the run
  subroutine exitfielddump
      use modglobal, only : lfielddump
      use modstat_nc, only : exitstat_nc
    implicit none

       if (lfielddump) call exitstat_nc(ncid)
       if (allocated(dbuf)) deallocate(dbuf)
  end subroutine exitfielddump

  !> Mark output field n as computed at dump time (one of the DERIV_* codes).
  subroutine register_derived(n, kind)
    implicit none
    integer, intent(in) :: n, kind
    nderived = nderived + 1
    derived_kind(n) = kind
    derived_slot(n) = nderived
  end subroutine register_derived

  !> Fill out (interior cells) with the diagnostic kind, one of the DERIV_* codes.
  !! The velocity derivatives span the cell (face to face), so they sit at the
  !! cell centre; the Poisson terms come from modpois.
  subroutine compute_derived(kind, out)
    use modfields, only : u0, v0, w0
    use modglobal, only : ib,ie,jb,je,kb,ke,dxfi,dyi,dzfi
    use modpois,   only : poisson_rhs_diagnostic
    implicit none
    integer, intent(in)  :: kind
    real,    intent(out) :: out(ib:ie,jb:je,kb:ke)
    integer i,j,k

    select case (kind)
    case (DERIV_DIV)
      do k=kb,ke
        do j=jb,je
          do i=ib,ie
            out(i,j,k) = (u0(i+1,j,k) - u0(i,j,k))*dxfi(i) &
                       + (v0(i,j+1,k) - v0(i,j,k))*dyi &
                       + (w0(i,j,k+1) - w0(i,j,k))*dzfi(k)
          end do
        end do
      end do
    case (DERIV_DUDX)
      do k=kb,ke
        do j=jb,je
          do i=ib,ie
            out(i,j,k) = (u0(i+1,j,k) - u0(i,j,k))*dxfi(i)
          end do
        end do
      end do
    case (DERIV_DVDY)
      do k=kb,ke
        do j=jb,je
          do i=ib,ie
            out(i,j,k) = (v0(i,j+1,k) - v0(i,j,k))*dyi
          end do
        end do
      end do
    case (DERIV_DWDZ)
      do k=kb,ke
        do j=jb,je
          do i=ib,ie
            out(i,j,k) = (w0(i,j,k+1) - w0(i,j,k))*dzfi(k)
          end do
        end do
      end do
    case (DERIV_PRHS)
      call poisson_rhs_diagnostic(0, out)
    case (DERIV_PDX)
      call poisson_rhs_diagnostic(1, out)
    case (DERIV_PDY)
      call poisson_rhs_diagnostic(2, out)
    case (DERIV_PDZ)
      call poisson_rhs_diagnostic(3, out)
    case default
      write(0, *) 'ERROR: fielddump: unknown derived field code ', kind
      stop 1
    end select
  end subroutine compute_derived

end module modfielddump
