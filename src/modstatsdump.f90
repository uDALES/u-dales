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
  use modglobal, only : dt,ltkedump,ltdump,lmintdump,ifoutput !,nstat
  use modmpi, only : myid
  implicit none
  private
  PUBLIC :: initstatsdump,statsdump,exitstatsdump
  save

  !NetCDF variables
  integer :: ncidtke,ncidkslice,ncidislice,ncidjslice,nrectke=0,&
             nreckslice=0,nrecislice=0,nrecjslice=0,nstattke=8,nstatkslice=5,nstatislice=5,nstatjslice=5,&
             ncidt,nrect=0,nstatt=32,ncidmint,nrecmint=0,nstatmint=6
  character(80) :: tkename = 'tkedump.xxx.nc'
  character(80) :: tname = 'tdump.xxx.xxx.xxx.nc'
  character(80) :: mintname = 'mintdump.xxx.xxx.xxx.nc'
  character(80) :: kslicename = 'kslicedump.xxx.xxx.xxx.nc'
  character(80) :: islicename = 'islicedump.xxx.xxx.xxx.nc'
  character(80) :: jslicename = 'jslicedump.xxx.xxx.xxx.nc'
  character(80),dimension(1,4) :: tncstattke
  character(80),dimension(1,4) :: tncstatkslice
  character(80),dimension(1,4) :: tncstatislice
  character(80),dimension(1,4) :: tncstatjslice
  character(80),dimension(1,4) :: tncstatt
  character(80),dimension(1,4) :: tncstatmint

  integer :: klow,khigh,i,j,k
  real    :: tsamplep,tstatsdumpp,tsample,tstatsdump

  integer :: isliceloc    ! local islice on core
  logical :: islicerank    ! cpu that islice is on
  integer :: jsliceloc    ! local jslice on core
  logical :: jslicerank    ! cpu that jslice is on


contains

  !--------------------------
  !> Initializing statsdump. Read out the namelist, initializing the variables
  !-------------------------

  subroutine initstatsdump
    use modmpi,   only : my_real,mpierr,comm3d,mpi_logical,mpi_integer,mpi_character,cmyid,cmyidx,cmyidy
    use modglobal,only : imax,jmax,kmax,cexpnr,ifnamopt,fname_options,ib,ie,jb,je,kb,ke,ladaptive,btime,&
                         nsv,lkslicedump,lislicedump,ljslicedump,ib,ie,islice,jslice
    use modstat_nc,only: open_nc, define_nc,ncinfo,writestat_dims_nc
    use modfields, only : ncstattke,ncstatkslice,ncstatislice,ncstatjslice,ncstatt,ncstatmint
    use decomp_2d, only : zstart, zend
    implicit none
    integer :: ierr

    namelist/NAMSTATSDUMP/ &
         tsample,klow,khigh,tstatsdump,ltkedump,ltdump,lmintdump    ! maybe removed; NAMSTATSDUMP is not in use anymore

    allocate(ncstattke(nstattke,4))
    allocate(ncstatkslice(nstatkslice,4))
    allocate(ncstatislice(nstatislice,4))
    allocate(ncstatjslice(nstatjslice,4))
    allocate(ncstatt(nstatt,4))
    allocate(ncstatmint(nstatmint,4))

    klow=kb
    khigh=ke

    if(myid==0)then
       open(ifnamopt,file=fname_options,status='old',iostat=ierr)
       read (ifnamopt,NAMSTATSDUMP,iostat=ierr)
       if (ierr > 0) then
          write(0, *) 'ERROR: Problem in namoptions NAMSTATSDUMP'
          write(0, *) 'iostat error: ', ierr
          stop 1
       endif
       !write(6 ,NAMSTATSDUMP)
       close(ifnamopt)
    end if

    call MPI_BCAST(klow        ,1,MPI_INTEGER,0,comm3d,ierr) !have to do this? just want nc for first CPU
    call MPI_BCAST(khigh       ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(nstatt      ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(nstatmint      ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(ncstattke   ,80,MPI_CHARACTER,0,comm3d,mpierr)
    call MPI_BCAST(ncstatt     ,80,MPI_CHARACTER,0,comm3d,mpierr)
    call MPI_BCAST(ncstatmint     ,80,MPI_CHARACTER,0,comm3d,mpierr)
    !call MPI_BCAST(ltdump      ,1,MPI_LOGICAL,0,comm3d,ierr)      ! maybe removed; unnecessary broadcast; this variable already broadcasted in modstartup
    call MPI_BCAST(lmintdump      ,1,MPI_LOGICAL,0,comm3d,ierr)

    !> Generate time averaged NetCDF: tdump.xxx.nc
    if (ltdump) then

      tname(7:9) = cmyidx
      tname(11:13) = cmyidy
      tname(15:17) = cexpnr
      call ncinfo(tncstatt(1,:),'time'      ,'Time'                        ,'s'      ,'time')
      call ncinfo(ncstatt( 1,:),'ut'        ,'Streamwise velocity'         ,'m/s'    ,'mttt'  )
      call ncinfo(ncstatt( 2,:),'vt'        ,'Spanwise velocity'           ,'m/s'    ,'tmtt'  )
      call ncinfo(ncstatt( 3,:),'wt'        ,'Vertical velocity'           ,'m/s'    ,'ttmt'  )
      call ncinfo(ncstatt( 4,:),'thlt'      ,'Temperature'                 ,'K'      ,'tttt'  )
      call ncinfo(ncstatt( 5,:),'qtt'       ,'Moisture'                    ,'kg/kg'  ,'tttt'  )
      call ncinfo(ncstatt( 6,:),'pt'        ,'Pressure'                    ,'m^2/s^2','tttt'  )
      call ncinfo(ncstatt( 7,:),'sca1t'     ,'Concentration field 1'       ,'g/m^3'  ,'tttt'  )
      call ncinfo(ncstatt( 8,:),'sca2t'     ,'Concentration field 2'       ,'g/m^3'  ,'tttt'  )
      call ncinfo(ncstatt( 9,:),'sca3t'     ,'Concentration field 3'       ,'g/m^3'  ,'tttt'  )
      call ncinfo(ncstatt(10,:),'sca4t'     ,'Concentration field 4'       ,'g/m^3'  ,'tttt'  )
      call ncinfo(ncstatt(11,:),'PSS'       ,'PSS defect'                  ,'gm/s'   ,'tttt'  )

      call ncinfo(ncstatt(12,:),'upwpt'     ,'Turbulent momentum flux'     ,'m^2/s^2','mtmt'  )
      call ncinfo(ncstatt(13,:),'vpwpt'     ,'Turbulent momentum flux'     ,'m^2/s^2','tmmt'  )
      call ncinfo(ncstatt(14,:),'upvpt'     ,'Turbulent momentum flux'     ,'m^2/s^2','mmtt'  )
      call ncinfo(ncstatt(15,:),'wpthlpt'   ,'Turbulent heat flux'         ,'K m/s'  ,'ttmt'  )
      call ncinfo(ncstatt(16,:),'wpsca1pt'  ,'Turbulent flux 1'            ,'gm/s'   ,'ttmt'  )
      call ncinfo(ncstatt(17,:),'wpsca2pt'  ,'Turbulent flux 2'            ,'gm/s'   ,'ttmt'  )
      call ncinfo(ncstatt(18,:),'wpsca3pt' ,'Turbulent flux 3'            ,'gm/s'   ,'ttmt'  )
      call ncinfo(ncstatt(19,:),'wpsca4pt'  ,'Turbulent flux 4'            ,'gm/s'   ,'ttmt'  )

      call ncinfo(ncstatt(20,:),'thlpthlpt','Temperature variance'        ,'K^2'    ,'tttt'  )
      call ncinfo(ncstatt(21,:),'upuptc'   ,'u variance'                  ,'m^2/s^2','tttt'  )
      call ncinfo(ncstatt(22,:),'vpvptc'   ,'v variance'                  ,'m^2/s^2','tttt'  )
      call ncinfo(ncstatt(23,:),'wpwptc'   ,'w variance'                  ,'m^2/s^2','tttt'  )
      call ncinfo(ncstatt(24,:),'tketc'    ,'TKE'                         ,'m^2/s^2','tttt'  )
      call ncinfo(ncstatt(25,:),'sca1psca1pt','Concentration variance 1'  ,'g^2/m^6','tttt'  )
      call ncinfo(ncstatt(26,:),'sca2psca2pt','Concentration variance 2'  ,'g^2/m^6','tttt'  )
      call ncinfo(ncstatt(27,:),'sca3psca3pt','Concentration variance 3'  ,'g^2/m^6','tttt'  )
      call ncinfo(ncstatt(28,:),'sca4psca4pt','Concentration variance 4'  ,'g^2/m^6','tttt'  )

      call ncinfo(ncstatt(29,:),'sv1sgs'   ,'SGS flux 1'                  ,'gm/s'   ,'ttmt'  )
      call ncinfo(ncstatt(30,:),'sv2sgs'   ,'SGS flux 2'                  ,'gm/s'   ,'ttmt'  )
      call ncinfo(ncstatt(31,:),'sv3sgs'   ,'SGS flux 3'                  ,'gm/s'   ,'ttmt'  )
      call ncinfo(ncstatt(32,:),'sv4sgs'   ,'SGS flux 4'                  ,'gm/s'   ,'ttmt'  )

      ! call ncinfo(ncstatt(33,:),'sca1t_max','Max concentration field 1'   ,'g/m^3'  ,'tttt'  )
      ! call ncinfo(ncstatt(34,:),'sca2t_max','Max concentration field 2'   ,'g/m^3'  ,'tttt'  )
      ! call ncinfo(ncstatt(35,:),'sca3t_max','Max concentration field 3'   ,'g/m^3'  ,'tttt'  )
      ! call ncinfo(ncstatt(36,:),'sca4t_max','Max concentration field 4'   ,'g/m^3'  ,'tttt'  )

!      if (myid==0) then
        call open_nc(tname, ncidt, nrect, n1=imax, n2=jmax, n3=khigh-klow+1)
        if (nrect==0) then
          call define_nc( ncidt, 1, tncstatt)
          call writestat_dims_nc(ncidt)
        end if
        call define_nc( ncidt, nstatt, ncstatt)
!      end if
    end if

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

    !      if (myid==0) then
        call open_nc(mintname, ncidmint, nrecmint, n1=imax, n2=jmax, n3=khigh-klow+1)
        if (nrecmint==0) then
          call define_nc( ncidmint, 1, tncstatmint)
          call writestat_dims_nc(ncidmint)
        end if
        call define_nc( ncidmint, nstatmint, ncstatmint)
    !      end if
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

    !> Generate sliced NetCDF: slicedump.xxx.xxx.nc
    if (lkslicedump) then

      kslicename(12:14) = cmyidx
      kslicename(16:18) = cmyidy
      kslicename(20:22) = cexpnr

      call ncinfo(tncstatkslice(1,:),'time'     ,'Time'   ,'s'   ,'time')
      call ncinfo(ncstatkslice( 1,:),'u_kslice'     ,'Streamwise velocity at kslice', '-', 'mt0t')
      call ncinfo(ncstatkslice( 2,:),'v_kslice'     ,'Spanwise velocity at kslice', '-', 'tm0t')
      call ncinfo(ncstatkslice( 3,:),'w_kslice'     ,'Vertical velocity at kslice', '-', 'tt0t')
      call ncinfo(ncstatkslice( 4,:),'thl_kslice'   ,'Potential temperature at kslice', '-', 'tt0t')
      call ncinfo(ncstatkslice( 5,:),'qt_kslice'    ,'Specific humidity at kslice', '-', 'tt0t')

      call open_nc(kslicename, ncidkslice, nreckslice, n1=imax, n2=jmax)

      if (nreckslice==0) then
        call define_nc( ncidkslice, 1, tncstatkslice)
        call writestat_dims_nc(ncidkslice)
      end if

      call define_nc( ncidkslice, nstatkslice, ncstatkslice)

    end if

    if (lislicedump) then

      islicename(12:14) = cmyidx
      islicename(16:18) = cmyidy
      islicename(20:22) = cexpnr

      call ncinfo(tncstatislice(1,:),'time'     ,'Time'   ,'s'   ,'time')
      call ncinfo(ncstatislice( 1,:),'u_islice'     ,'Streamwise velocity at islice', '-', '0ttt')
      call ncinfo(ncstatislice( 2,:),'v_islice'     ,'Spanwise velocity at islice', '-', '0mtt')
      call ncinfo(ncstatislice( 3,:),'w_islice'     ,'Vertical velocity at islice', '-', '0tmt')
      call ncinfo(ncstatislice( 4,:),'thl_islice'   ,'Potential temperature at islice', '-', '0ttt')
      call ncinfo(ncstatislice( 5,:),'qt_islice'    ,'Specific humidity at islice', '-', '0ttt')

      if ((islice >= zstart(1)) .and. (islice <= zend(1))) then
        islicerank = .true.
        isliceloc = islice - zstart(1) + 1
      else
        islicerank = .false.
      end if

      if (islicerank) then
        call open_nc(islicename, ncidislice, nrecislice, n2=jmax, n3=khigh-klow+1)
        if (nrecislice==0) then
          call define_nc( ncidislice, 1, tncstatislice)
          call writestat_dims_nc(ncidislice)
        end if
        call define_nc( ncidislice, nstatislice, ncstatislice)
      end if

    end if

    if (ljslicedump) then

      jslicename(12:14) = cmyidx
      jslicename(16:18) = cmyidy
      jslicename(20:22) = cexpnr

      call ncinfo(tncstatjslice(1,:),'time'     ,'Time'   ,'s'   ,'time')
      call ncinfo(ncstatjslice( 1,:),'u_jslice'     ,'Streamwise velocity at jslice', '-', 'm0tt')
      call ncinfo(ncstatjslice( 2,:),'v_jslice'     ,'Spanwise velocity at jslice', '-', 't0tt')
      call ncinfo(ncstatjslice( 3,:),'w_jslice'     ,'Vertical velocity at jslice', '-', 't0mt')
      call ncinfo(ncstatjslice( 4,:),'thl_jslice'   ,'Potential temperature at jslice', '-', 't0tt')
      call ncinfo(ncstatjslice( 5,:),'qt_jslice'    ,'Specific humidity at jslice', '-', 't0tt')

      if ((jslice >= zstart(2)) .and. (jslice <= zend(2))) then
        jslicerank = .true.
        jsliceloc = jslice - zstart(2) + 1
      else
        jslicerank = .false.
      end if

      if (jslicerank) then
         call open_nc(jslicename, ncidjslice, nrecjslice, n1=imax, n3=khigh-klow+1)
         if (nrecjslice==0) then
            call define_nc( ncidjslice, 1, tncstatjslice)
            call writestat_dims_nc(ncidjslice)
         end if
         call define_nc( ncidjslice, nstatjslice, ncstatjslice)
      end if

    end if

    !> Set times to zero so works for warm starts... could have issues with warmstarts here...
    tsamplep = 0.
    tstatsdumpp = 0.

  end subroutine initstatsdump

  !-------------------------
  !> Generate and write statistics into NetCDF file format
  !-------------------------

  subroutine statsdump

  use modfields,        only : um,up,vm,wm,svm,qtm,thlm,pres0,ncstattke,ncstatmint,&
                               ncstatkslice,ncstatislice,ncstatjslice,t_t,t_v,t_p,t_sgs,d_sgs,p_b,p_t,adv,&
                               IIc,IIu,IIv,IIw,IIuw,IIvw,IIct,IIwt,IIut,IIvt,IIuwt,IIuv,&
                               IIcs,IIws,IIus,IIvs,IIuws,IIvws,IIuvs,&
                               uwtik,vwtjk,uvtij,utik,wtik,wtjk,vtjk,utij,vtij,&
                               wthltk,wqttk,thlthlt,qtqtt,sv1sv1t,sv2sv2t,sv3sv3t,sv4sv4t,wmt,thltk,qttk,thlt,&
                               ncstatt,uutc,vvtc,wwtc,utc,vtc,wtc,&
                               umt,vmt,sv1t,sv2t,sv3t,sv4t,sv1tk,sv2tk,sv3tk,sv4tk,wsv1tk,wsv2tk,wsv3tk,wsv4tk,&
                               sv1sgst,sv2sgst,sv3sgst,sv4sgst,qtt,pt,PSSt !,sv1max,sv2max,sv3max,sv4max
  use modglobal,        only : ib,ie,ih,ihc,xf,xh,jb,je,jhc,jgb,jge,dy,dyi,jh,ke,kb,kh,khc,rk3step,&
                               timee,cexpnr,tsample,tstatsdump,tstatstart,jtot,imax,jmax,dzf,&
                               ltempeq,zh,dxf,dzf,dzh2i,lprofforc,lscasrcl,&
                               lkslicedump,lislicedump,ljslicedump,lchem,dzhi,dzfi,dzhiq,dxhi,lmoist,nsv,&
                               k1,JNO2,lchem,kslice,islice,jslice
!  use modsubgriddata,   only : ekm,sbshr
  use modstat_nc,       only : writestat_nc,writestat_1D_nc
  use modmpi,           only : myid,cmyid,my_real,mpi_sum,avey_ibm,mpierr,&
                               comm3d,nprocs,myidx,myidy
  use modsurfdata,      only : thls
  use modsubgrid,       only : ekh,ekm
  use modstatistics,    only : genstats,tkestats
  ! use, intrinsic :: ieee_arithmetic
  implicit none

  !> Create fields to be used in statistics

  ! interpolated fields
!  real, dimension(ib:ie,jb:je,kb:ke)     :: umc
!  real, dimension(ib:ie,jb:je,kb:ke)     :: vmc
!  real, dimension(ib:ie,jb:je,kb:ke)     :: wmc

 ! real, dimension(ib:ie,jb:je,kb:ke+kh)     :: thlk
 ! real, dimension(ib:ie,jb:je,kb:ke+kh)     :: qtk
 ! real, dimension(ib:ie,jb:je,kb:ke+kh)     :: uik
 ! real, dimension(ib:ie,jb:je,kb:ke+kh)     :: wik
 ! real, dimension(ib:ie,jb:je,kb:ke+kh)     :: vjk
 ! real, dimension(ib:ie,jb:je,kb:ke+kh)     :: wjk
 ! real, dimension(ib:ie,jb:je,kb:ke+kh)     :: uij
 ! real, dimension(ib:ie,jb:je,kb:ke+kh)     :: vij
 ! real, dimension(ib:ie,jb:je,kb:ke+kh)     :: uc
 ! real, dimension(ib:ie,jb:je,kb:ke+kh)     :: vc
 ! real, dimension(ib:ie,jb:je,kb:ke+kh)     :: wc

  real, allocatable     :: thlk(:,:,:)
  real, allocatable     :: qtk(:,:,:)
  real, allocatable     :: uik(:,:,:)
  real, allocatable     :: wik(:,:,:)
  real, allocatable     :: vjk(:,:,:)
  real, allocatable     :: wjk(:,:,:)
  real, allocatable     :: uij(:,:,:)
  real, allocatable     :: vij(:,:,:)
  real, allocatable     :: uc(:,:,:)
  real, allocatable     :: vc(:,:,:)
  real, allocatable     :: wc(:,:,:)

  ! SGS fluxes
  ! real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)     :: thlsgs
  ! real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)     :: qtsgs
  ! real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)     :: usgs
  ! real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)     :: vsgs
  ! real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)     :: wsgs

  real, allocatable     :: thlsgs(:,:,:)
  real, allocatable     :: qtsgs(:,:,:)
  real, allocatable     :: usgs(:,:,:)
  real, allocatable     :: vsgs(:,:,:)
  real, allocatable     :: wsgs(:,:,:)

  ! t-averaged fields
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: sv1k
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: sv2k
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: sv3k
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: sv4k
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: wpsv1p
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: wpsv2p
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: wpsv3p
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: wpsv4p
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: sv1sgs
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: sv2sgs
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: sv3sgs
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: sv4sgs

  real, allocatable     :: sv1k(:,:,:)
  real, allocatable     :: sv2k(:,:,:)
  real, allocatable     :: sv3k(:,:,:)
  real, allocatable     :: sv4k(:,:,:)
  real, allocatable     :: sv1sgs(:,:,:)
  real, allocatable     :: sv2sgs(:,:,:)
  real, allocatable     :: sv3sgs(:,:,:)
  real, allocatable     :: sv4sgs(:,:,:)
  real, allocatable     :: wpsv1p(:,:,:)
  real, allocatable     :: wpsv2p(:,:,:)
  real, allocatable     :: wpsv3p(:,:,:)
  real, allocatable     :: wpsv4p(:,:,:)
  real, allocatable     :: sv1psv1pt(:,:,:)
  real, allocatable     :: sv2psv2pt(:,:,:)
  real, allocatable     :: sv3psv3pt(:,:,:)
  real, allocatable     :: sv4psv4pt(:,:,:)
  real, allocatable     :: PSS(:,:,:)

  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: upwptik
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: vpwptjk
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: upvptij
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: wpthlptk
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: thlpthlpt
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: upuptc
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: vpvptc
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: wpwptc
  ! real, dimension(ib:ie,jb:je,kb:ke+kh)        :: tketc

  real, allocatable     :: upwptik(:,:,:)
  real, allocatable     :: vpwptjk(:,:,:)
  real, allocatable     :: upvptij(:,:,:)
  real, allocatable     :: wpthlptk(:,:,:)
  real, allocatable     :: thlpthlpt(:,:,:)
  real, allocatable     :: upuptc(:,:,:)
  real, allocatable     :: vpvptc(:,:,:)
  real, allocatable     :: wpwptc(:,:,:)
  real, allocatable     :: tketc(:,:,:)

  real, allocatable :: varstke(:,:),&
                       varkslice(:,:,:),varislice(:,:,:),varjslice(:,:,:),varst(:,:,:,:),varsmint(:,:,:,:)
  real    :: tstatsdumppi,emom
  integer :: i,j,k,ip,im,jp,jm,kp,km,n
  integer :: writecounter = 1
  integer :: reclength

  if (timee < tstatstart) return

  if (.not.(ltdump .or. lmintdump &
    .or. lkslicedump.or. lislicedump.or. ljslicedump)) return

  allocate(thlk(ib:ie,jb:je,kb:ke+kh))
  allocate(qtk(ib:ie,jb:je,kb:ke+kh))
  allocate(uik(ib:ie,jb:je,kb:ke+kh))
  allocate(wik(ib:ie,jb:je,kb:ke+kh))
  allocate(vjk(ib:ie,jb:je,kb:ke+kh))
  allocate(wjk(ib:ie,jb:je,kb:ke+kh))
  allocate(uij(ib:ie,jb:je,kb:ke+kh))
  allocate(vij(ib:ie,jb:je,kb:ke+kh))
  allocate(uc(ib:ie,jb:je,kb:ke+kh))
  allocate(vc(ib:ie,jb:je,kb:ke+kh))
  allocate(wc(ib:ie,jb:je,kb:ke+kh))

  allocate(thlsgs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
  allocate(qtsgs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
  allocate(usgs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
  allocate(vsgs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))
  allocate(wsgs(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh))

  allocate(sv1k(ib:ie,jb:je,kb:ke+kh))
  allocate(sv2k(ib:ie,jb:je,kb:ke+kh))
  allocate(sv3k(ib:ie,jb:je,kb:ke+kh))
  allocate(sv4k(ib:ie,jb:je,kb:ke+kh))
  allocate(sv1sgs(ib:ie,jb:je,kb:ke+kh))
  allocate(sv2sgs(ib:ie,jb:je,kb:ke+kh))
  allocate(sv3sgs(ib:ie,jb:je,kb:ke+kh))
  allocate(sv4sgs(ib:ie,jb:je,kb:ke+kh))
  allocate(wpsv1p(ib:ie,jb:je,kb:ke+kh))
  allocate(wpsv2p(ib:ie,jb:je,kb:ke+kh))
  allocate(wpsv3p(ib:ie,jb:je,kb:ke+kh))
  allocate(wpsv4p(ib:ie,jb:je,kb:ke+kh))
  allocate(sv1psv1pt(ib:ie,jb:je,kb:ke+kh))
  allocate(sv2psv2pt(ib:ie,jb:je,kb:ke+kh))
  allocate(sv3psv3pt(ib:ie,jb:je,kb:ke+kh))
  allocate(sv4psv4pt(ib:ie,jb:je,kb:ke+kh))
  allocate(PSS(ib:ie,jb:je,kb:ke+kh))

  allocate(upwptik(ib:ie,jb:je,kb:ke+kh))
  allocate(vpwptjk(ib:ie,jb:je,kb:ke+kh))
  allocate(upvptij(ib:ie,jb:je,kb:ke+kh))
  allocate(wpthlptk(ib:ie,jb:je,kb:ke+kh))
  allocate(thlpthlpt(ib:ie,jb:je,kb:ke+kh))
  allocate(upuptc(ib:ie,jb:je,kb:ke+kh))
  allocate(vpvptc(ib:ie,jb:je,kb:ke+kh))
  allocate(wpwptc(ib:ie,jb:je,kb:ke+kh))
  allocate(tketc(ib:ie,jb:je,kb:ke+kh))

  ! thlk=0.;qtk=0.;uik=0.;wik=0.;vjk=0.;wjk=0.;uij=0.;vij=0.;uc=0.;vc=0.;wc=0.;thlsgs=0.;qtsgs=0.;usgs=0.;vsgs=0.;wsgs=0.;sv1k=0.;sv2k=0.;sv3k=0.;sv4k=0.;sv1sgs=0.;sv2sgs=0.;sv3sgs=0.;sv4sgs=0.;wpsv1p=0.;wpsv2p=0.
  ! wpsv3p=0.;wpsv4p=0.;sv1psv1pt=0.;sv2psv2pt=0.;sv3psv3pt=0.;sv4psv4pt=0.;PSS=0.;upwptik=0.;vpwptjk=0.;upvptij=0.;wpthlptk=0.;thlpthlpt=0.;upuptc=0.;vpvptc=0.;wpwptc=0.;tketc=0.

  wpsv1p=0.;wpsv2p=0.;wpsv3p=0.;wpsv4p=0.
  sv1psv1pt=0.;sv2psv2pt=0.;sv3psv3pt=0.;sv4psv4pt=0.

  if (.not. rk3step==3)  return

  if (tsamplep > tsample) then

    if (ltdump .or. lmintdump) then

      tstatsdumppi = 1./tstatsdumpp

      !> Perform required interpolations for flux calculations
      !  tg3315 for non-equidistant x and z-grids this needs to change
      do k=kb,ke+kh
        do j=jb,je
          do i=ib,ie

            uik(i,j,k) = 0.5*dzhi(k)*(um(i,j,k)*dzf(k-1) + um(i,j,k-1)*dzf(k))
            wik(i,j,k) = 0.5*dxhi(i)*(wm(i,j,k)*dxf(i-1) + wm(i-1,j,k)*dxf(i))
            vjk(i,j,k) = 0.5*dzhi(k)*(vm(i,j,k)*dzf(k-1) + vm(i,j,k-1)*dzf(k))
            wjk(i,j,k) = 0.5*        (wm(i,j,k)          + wm(i,j-1,k))
            uij(i,j,k) = 0.5*        (um(i,j,k)          + um(i,j-1,k))
            vij(i,j,k) = 0.5*dxhi(i)*(vm(i,j,k)*dxf(i-1) + vm(i-1,j,k)*dxf(i))
            uc (i,j,k) = 0.5*        (um(i+1,j,k)        + um(i,j,k))
            vc (i,j,k) = 0.5*        (vm(i,j+1,k)        + vm(i,j,k))
            if (k==ke+kh) then
              wc(i,j,k) = wc(i,j,k-1)
            else
              wc(i,j,k) = 0.5*        (wm(i,j,k+1)        + wm(i,j,k))
            end if

            ! SGS fluxes
            ! interps ekm to cell corner (uw)
            emom = ( dzf(k-1) * ( ekm(i,j,k)*dxf(i-1)  + ekm(i-1,j,k)*dxf(i) )  + &
                     dzf(k)   * ( ekm(i,j,k-1)*dxf(i-1) + ekm(i-1,j,k-1)*dxf(i) ) )*dxhi(i) * dzhiq(k)
            usgs(i,j,k)  = emom * ( (um(i,j,k)-um(i,j,k-1)) *dzhi(k) &
                        +(wm(i,j,k)-wm(i-1,j,k))  *dxhi(i))

            ! interps ekm to cell corner (vw)
            emom = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i,j-1,k) )  + &
                     dzf(k)   * ( ekm(i,j,k-1) + ekm(i,j-1,k-1) ) ) * dzhiq(k)

            vsgs(i,j,k)  = emom * ( (vm(i,j,k)-vm(i,j,k-1)) *dzhi(k) &
                        +(wm(i,j,k)-wm(i,j-1,k))  *dyi)

         end do
       end do
     end do

     do k=kb,ke
       do j=jb,je
         do i=ib,ie
           wsgs(i,j,k) = ( ekm(i,j,k) * (wm(i,j,k+1)-wm(i,j,k)) *dzfi(k) &
                     -ekm(i,j,k-1)* (wm(i,j,k)-wm(i,j,k-1)) *dzfi(k-1) ) * 2. &
                     * dzhi(k) ! tg3315 check this
         end do
       end do
     end do

    if (ltempeq) then
      do k=kb,ke+kh
        do j=jb,je
          do i=ib,ie
              thlk(i,j,k) = 0.5*dzhi(k)*(thlm(i,j,k)*dzf(k-1) + thlm(i,j,k-1)*dzf(k))
          end do
        end do
      end do
      do k=kb,ke
        !> SGS fluxes
        thlsgs(:,:,k) = 0.5 * (dzf(k-1)*ekh(:,:,k) + dzf(k)*ekh(:,:,k-1)) &
                        * (thlm(:,:,k)-thlm(:,:,k-1)) * dzh2i(k)
      end do
    end if

    if (lmoist) then
      do k=kb,ke+kh
        do j=jb,je
          do i=ib,ie
              qtk(i,j,k) = 0.5*dzhi(k)*(qtm(i,j,k)*dzf(k-1) + qtm(i,j,k-1)*dzf(k))
          end do
        end do
      end do
      do k=kb,ke
        !> SGS fluxes
        qtsgs(:,:,k) = 0.5 * (dzf(k-1)*ekh(:,:,k) + dzf(k)*ekh(:,:,k-1)) &
                        * (qtm(:,:,k)-qtm(:,:,k-1)) * dzh2i(k)
      end do
    end if

    if (nsv>0) then
      do k=kb,ke
        do j=jb,je
          do i=ib,ie
              sv1k(i,j,k) = 0.5*dzhi(k)*(svm(i,j,k,1)*dzf(k-1) + svm(i,j,k-1,1)*dzf(k))
          end do
        end do
      end do
      do k=kb,ke
        sv1sgs(ib:ie,jb:je,k) = 0.5 * (dzf(k-1)*ekh(ib:ie,jb:je,k) + dzf(k)*ekh(ib:ie,jb:je,k-1)) &
                        * (svm(ib:ie,jb:je,k,1)-svm(ib:ie,jb:je,k-1,1)) * dzh2i(k)
      end do
    end if

    if (nsv>1) then
      do k=kb,ke+kh
        do j=jb,je
          do i=ib,ie
              sv2k(i,j,k) = 0.5*dzhi(k)*(svm(i,j,k,2)*dzf(k-1) + svm(i,j,k-1,2)*dzf(k))
          end do
        end do
      end do
      do k=kb,ke
        sv2sgs(ib:ie,jb:je,k) = 0.5 * (dzf(k-1)*ekh(ib:ie,jb:je,k) + dzf(k)*ekh(ib:ie,jb:je,k-1)) &
                        * (svm(ib:ie,jb:je,k,2)-svm(ib:ie,jb:je,k-1,2)) * dzh2i(k)
      end do
    end if

    if (nsv>2) then
      do k=kb,ke+kh
        do j=jb,je
          do i=ib,ie
              sv3k(i,j,k) = 0.5*dzhi(k)*(svm(i,j,k,3)*dzf(k-1) + svm(i,j,k-1,3)*dzf(k))
          end do
        end do
      end do
      do k=kb,ke
        sv3sgs(ib:ie,jb:je,k) = 0.5 * (dzf(k-1)*ekh(ib:ie,jb:je,k) + dzf(k)*ekh(ib:ie,jb:je,k-1)) &
                        * (svm(ib:ie,jb:je,k,3)-svm(ib:ie,jb:je,k-1,3)) * dzh2i(k)
      end do
    end if

    if (nsv>3) then
      do k=kb,ke+kh
        do j=jb,je
          do i=ib,ie
              sv4k(i,j,k) = 0.5*dzhi(k)*(svm(i,j,k,4)*dzf(k-1) + svm(i,j,k-1,4)*dzf(k))
          end do
        end do
      end do
      do k=kb,ke
        sv4sgs(ib:ie,jb:je,k) = 0.5 * (dzf(k-1)*ekh(ib:ie,jb:je,k) + dzf(k)*ekh(ib:ie,jb:je,k-1)) &
                        * (svm(ib:ie,jb:je,k,4)-svm(ib:ie,jb:je,k-1,4)) * dzh2i(k)
      end do
    end if

    if ((nsv>2) .and. (lchem .eqv. .true.)) then
      do k=kb,ke
        do j=jb,je
          do i=ib,ie
            if ((ABS(svm(i,j,k,2)) .gt. 1.e-40) .and. (IIc(i,j,k)==1)) then
              PSS(i,j,k) = ( ( (k1*(svm(i,j,k,1)/30.)*(svm(i,j,k,3)/48.))/(JNO2*(svm(i,j,k,2)/46.)) ) - 1 ) * 100
            end if
          end do
        end do
      end do
    end if

      !!>> CALCS FOR TIME DEPENDANT (AVERAGED) STATS

      ! Average 3-D fields in time
      ! tg3315 may be possible to do less calculations by splitting up
      if (ltdump .or. lmintdump) then
        uwtik(:,:,kb:ke+kh) = (uwtik(:,:,kb:ke+kh)*(tstatsdumpp-tsamplep) + wik(:,:,kb:ke+kh)*uik(:,:,kb:ke+kh)*tsamplep)*tstatsdumppi
        vwtjk(:,:,kb:ke+kh) = (vwtjk(:,:,kb:ke+kh)*(tstatsdumpp-tsamplep) + wjk(:,:,kb:ke+kh)*vjk(:,:,kb:ke+kh)*tsamplep)*tstatsdumppi
        uvtij(:,:,kb:ke+kh) = (uvtij(:,:,kb:ke+kh)*(tstatsdumpp-tsamplep) + uij(:,:,kb:ke+kh)*vij(:,:,kb:ke+kh)*tsamplep)*tstatsdumppi
        uutc(ib:ie,jb:je,kb:ke+kh) = (uutc(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + uc(ib:ie,jb:je,kb:ke+kh)*uc(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        vvtc(ib:ie,jb:je,kb:ke+kh) = (vvtc(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + vc(ib:ie,jb:je,kb:ke+kh)*vc(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        wwtc(ib:ie,jb:je,kb:ke+kh) = (wwtc(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wc(ib:ie,jb:je,kb:ke+kh)*wc(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        utik(:,:,kb:ke+kh) = (utik(:,:,kb:ke+kh)*(tstatsdumpp-tsamplep) + uik(:,:,kb:ke+kh)*tsamplep)*tstatsdumppi
        wtik(:,:,kb:ke+kh) = (wtik(:,:,kb:ke+kh)*(tstatsdumpp-tsamplep) + wik(:,:,kb:ke+kh)*tsamplep)*tstatsdumppi
        vtjk(:,:,kb:ke+kh) = (vtjk(:,:,kb:ke+kh)*(tstatsdumpp-tsamplep) + vjk(:,:,kb:ke+kh)*tsamplep)*tstatsdumppi
        wtjk(:,:,kb:ke+kh) = (wtjk(:,:,kb:ke+kh)*(tstatsdumpp-tsamplep) + wjk(:,:,kb:ke+kh)*tsamplep)*tstatsdumppi
        utij(:,:,kb:ke+kh) = (utij(:,:,kb:ke+kh)*(tstatsdumpp-tsamplep) + uij(:,:,kb:ke+kh)*tsamplep)*tstatsdumppi
        vtij(:,:,kb:ke+kh) = (vtij(:,:,kb:ke+kh)*(tstatsdumpp-tsamplep) + vij(:,:,kb:ke+kh)*tsamplep)*tstatsdumppi
        umt(ib:ie,jb:je,kb:ke+kh) = (umt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + um(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        vmt(ib:ie,jb:je,kb:ke+kh) = (vmt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + vm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        wmt(ib:ie,jb:je,kb:ke+kh) = (wmt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        utc(ib:ie,jb:je,kb:ke+kh) = (utc(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + uc(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        vtc(ib:ie,jb:je,kb:ke+kh) = (vtc(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + vc(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        wtc(ib:ie,jb:je,kb:ke+kh) = (wtc(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wc(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi

        pt(ib:ie,jb:je,kb:ke+kh) = (pt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + pres0(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi

        if (ltempeq) then
          wthltk(ib:ie,jb:je,kb:ke+kh) = (wthltk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + thlk(ib:ie,jb:je,kb:ke+kh)*wm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          thlthlt(ib:ie,jb:je,kb:ke+kh) = (thlthlt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + thlm(ib:ie,jb:je,kb:ke+kh)*thlm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          thltk(ib:ie,jb:je,kb:ke+kh) = (thltk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + thlk(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          thlt(ib:ie,jb:je,kb:ke+kh) = (thlt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + thlm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        end if

        if (lmoist) then
          wqttk(ib:ie,jb:je,kb:ke+kh) = (wqttk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + qtk(ib:ie,jb:je,kb:ke+kh)*wm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          qtqtt(ib:ie,jb:je,kb:ke+kh) = (qtqtt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + qtm(ib:ie,jb:je,kb:ke+kh)*qtm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          qttk(ib:ie,jb:je,kb:ke+kh) = (qttk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + qtk(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          qtt(ib:ie,jb:je,kb:ke+kh) = (qtt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + qtm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        end if

        if (nsv>0) then
          sv1t(ib:ie,jb:je,kb:ke+kh) = (sv1t(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + svm(ib:ie,jb:je,kb:ke+kh,1)*tsamplep)*tstatsdumppi
          sv1tk(ib:ie,jb:je,kb:ke+kh) = (sv1tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv1k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          wsv1tk(ib:ie,jb:je,kb:ke+kh) = (wsv1tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wm(ib:ie,jb:je,kb:ke+kh)*sv1k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          sv1sgst(ib:ie,jb:je,kb:ke+kh) = (sv1sgst(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv1sgs(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          sv1sv1t(ib:ie,jb:je,kb:ke+kh) = (sv1sv1t(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + svm(ib:ie,jb:je,kb:ke+kh,1)*svm(ib:ie,jb:je,kb:ke+kh,1)*tsamplep)*tstatsdumppi
          ! sv1max(ib:ie,jb:je,kb:ke) = max(sv1max(ib:ie,jb:je,kb:ke),svm(ib:ie,jb:je,kb:ke,1))
        end if

        if ((lchem .eqv. .true.) .and. (nsv>2)) then
          PSSt(ib:ie,jb:je,kb:ke+kh) = (PSSt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + PSS(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        end if

        if (nsv>1) then
          sv2t(ib:ie,jb:je,kb:ke+kh) = (sv2t(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + svm(ib:ie,jb:je,kb:ke+kh,2)*tsamplep)*tstatsdumppi
          sv2tk(ib:ie,jb:je,kb:ke+kh) = (sv2tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv2k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          wsv2tk(ib:ie,jb:je,kb:ke+kh) = (wsv2tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wm(ib:ie,jb:je,kb:ke+kh)*sv2k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          sv2sgst(ib:ie,jb:je,kb:ke+kh) = (sv2sgst(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv2sgs(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          sv2sv2t(ib:ie,jb:je,kb:ke+kh) = (sv2sv2t(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + svm(ib:ie,jb:je,kb:ke+kh,2)*svm(ib:ie,jb:je,kb:ke+kh,2)*tsamplep)*tstatsdumppi
          ! sv2max(ib:ie,jb:je,kb:ke) = max(sv2max(ib:ie,jb:je,kb:ke),svm(ib:ie,jb:je,kb:ke,2))
        end if

        if (nsv>2) then
          sv3t(ib:ie,jb:je,kb:ke+kh) = (sv3t(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + svm(ib:ie,jb:je,kb:ke+kh,3)*tsamplep)*tstatsdumppi
          sv3tk(ib:ie,jb:je,kb:ke+kh) = (sv3tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv3k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          wsv3tk(ib:ie,jb:je,kb:ke+kh) = (wsv3tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wm(ib:ie,jb:je,kb:ke+kh)*sv3k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          sv3sgst(ib:ie,jb:je,kb:ke+kh) = (sv3sgst(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv3sgs(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          sv3sv3t(ib:ie,jb:je,kb:ke+kh) = (sv3sv3t(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + svm(ib:ie,jb:je,kb:ke+kh,3)*svm(ib:ie,jb:je,kb:ke+kh,3)*tsamplep)*tstatsdumppi
          ! sv3max(ib:ie,jb:je,kb:ke) = max(sv3max(ib:ie,jb:je,kb:ke),svm(ib:ie,jb:je,kb:ke,3))
        end if

        if (nsv>3) then
          sv4t(ib:ie,jb:je,kb:ke+kh) = (sv4t(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + svm(ib:ie,jb:je,kb:ke+kh,4)*tsamplep)*tstatsdumppi
          sv4tk(ib:ie,jb:je,kb:ke+kh) = (sv4tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv4k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          wsv4tk(ib:ie,jb:je,kb:ke+kh) = (wsv4tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wm(ib:ie,jb:je,kb:ke+kh)*sv4k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          sv4sgst(ib:ie,jb:je,kb:ke+kh) = (sv4sgst(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv4sgs(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
          sv4sv4t(ib:ie,jb:je,kb:ke+kh) = (sv4sv4t(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + svm(ib:ie,jb:je,kb:ke+kh,4)*svm(ib:ie,jb:je,kb:ke+kh,4)*tsamplep)*tstatsdumppi
          ! sv4max(ib:ie,jb:je,kb:ke) = max(sv4max(ib:ie,jb:je,kb:ke),svm(ib:ie,jb:je,kb:ke,4))
        end if

      end if !ltdump

      ! Other 3-D fields specifically for tdump
      !if (ltdump) then
        ! bss116 already calculated above
        ! wmt(ib:ie,jb:je,kb:ke+kh) = (wmt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        ! sv1t(ib:ie,jb:je,kb:ke+kh) = (sv1t(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + svm(ib:ie,jb:je,kb:ke+kh,1)*tsamplep)*tstatsdumppi
        ! sv2t(ib:ie,jb:je,kb:ke+kh) = (sv2t(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + svm(ib:ie,jb:je,kb:ke+kh,2)*tsamplep)*tstatsdumppi
        ! sv3t(ib:ie,jb:je,kb:ke+kh) = (sv3t(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + svm(ib:ie,jb:je,kb:ke+kh,3)*tsamplep)*tstatsdumppi
        ! ! sv4t(ib:ie,jb:je,kb:ke+kh) = (sv4t(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + svm(ib:ie,jb:je,kb:ke+kh,4)*tsamplep)*tstatsdumppi
        ! sv1tk(ib:ie,jb:je,kb:ke+kh) = (sv1tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv1k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        ! sv2tk(ib:ie,jb:je,kb:ke+kh) = (sv2tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv2k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        ! sv3tk(ib:ie,jb:je,kb:ke+kh) = (sv3tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv3k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        ! sv4tk(ib:ie,jb:je,kb:ke+kh) = (sv4tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv4k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        ! wsv1tk(ib:ie,jb:je,kb:ke+kh) = (wsv1tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wm(ib:ie,jb:je,kb:ke+kh)*sv1k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        ! wsv2tk(ib:ie,jb:je,kb:ke+kh) = (wsv2tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wm(ib:ie,jb:je,kb:ke+kh)*sv2k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        ! wsv3tk(ib:ie,jb:je,kb:ke+kh) = (wsv3tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wm(ib:ie,jb:je,kb:ke+kh)*sv3k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        ! wsv4tk(ib:ie,jb:je,kb:ke+kh) = (wsv4tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wm(ib:ie,jb:je,kb:ke+kh)*sv4k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        ! sv1sgst(ib:ie,jb:je,kb:ke+kh) = (sv1sgst(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv1sgs(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        ! sv2sgst(ib:ie,jb:je,kb:ke+kh) = (sv2sgst(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv2sgs(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        ! sv3sgst(ib:ie,jb:je,kb:ke+kh) = (sv3sgst(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv3sgs(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        ! sv4sgst(ib:ie,jb:je,kb:ke+kh) = (sv4sgst(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv4sgs(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
      !end if ! ltdump

    end if ! ltdump .or. lmintdump

    if (lkslicedump) then
     allocate(varkslice(imax,jmax,nstatkslice))
     call writestat_nc(ncidkslice,1,tncstatkslice,(/timee/),nreckslice,.true.)
     varkslice(:,:,1) = um(ib:ie,jb:je,kslice)
     varkslice(:,:,2) = vm(ib:ie,jb:je,kslice)
     varkslice(:,:,3) = 0.5*(wm(ib:ie,jb:je,kslice)+wm(ib:ie,jb:je,kslice+1)) ! assumes equidistant
     varkslice(:,:,4) = thlm(ib:ie,jb:je,kslice)
     varkslice(:,:,5) = qtm(ib:ie,jb:je,kslice)
     call writestat_nc(ncidkslice,nstatkslice,ncstatkslice,varkslice,nreckslice,imax,jmax)

    endif

    if (lislicedump) then
      if (islicerank) then
        allocate(varislice(jmax,khigh-klow+1,nstatislice))
        call writestat_nc(ncidislice,1,tncstatislice,(/timee/),nrecislice,.true.)
        varislice(:,:,1) = 0.5*(um(isliceloc,jb:je,kb:ke)+um(isliceloc+1,jb:je,kb:ke))
        varislice(:,:,2) = vm(isliceloc,jb:je,kb:ke)
        varislice(:,:,3) = wm(isliceloc,jb:je,kb:ke)
        varislice(:,:,4) = thlm(isliceloc,jb:je,kb:ke)
        varislice(:,:,5) = qtm(isliceloc,jb:je,kb:ke)
        call writestat_nc(ncidislice,nstatislice,ncstatislice,varislice,nrecislice,jmax,khigh-klow+1)

      endif
    endif

    if (ljslicedump) then
       if (jslicerank) then
         allocate(varjslice(imax,khigh-klow+1,nstatjslice))
         call writestat_nc(ncidjslice,1,tncstatjslice,(/timee/),nrecjslice,.true.)
         varjslice(:,:,1) = um(ib:ie,jsliceloc,kb:ke)
         varjslice(:,:,2) = 0.5*(vm(ib:ie,jsliceloc,kb:ke)+vm(ib:ie,jsliceloc+1,kb:ke))
         varjslice(:,:,3) = wm(ib:ie,jsliceloc,kb:ke)
         varjslice(:,:,4) = thlm(ib:ie,jsliceloc,kb:ke)
         varjslice(:,:,5) = qtm(ib:ie,jsliceloc,kb:ke)
         call writestat_nc(ncidjslice,nstatjslice,ncstatjslice,varjslice,nrecjslice,imax,khigh-klow+1)

      endif
    endif

    if (ltkedump) then
      !call genstats(tsamplep,tstatsdumpp,umc,vmc,wmc)
    endif
    tsamplep = dt
  else !timestatsdumpp < tsample

    tsamplep = tsamplep + dt

  endif

  if (tstatsdumpp > tstatsdump) then

    ! Final calculations and write t-averaged statistics every tsample
    if (ltdump) then

    ! wpsv1p = wsv1tk - wmt*sv1tk
    ! wpsv2p = wsv2tk - wmt*sv2tk
    ! wpsv3p = wsv3tk - wmt*sv3tk
    ! wpsv4p = wsv4tk - wmt*sv4tk

    wpthlptk=0.;thlpthlpt=0.

    !> Turbulent fluxes
    upwptik = uwtik - utik*wtik
    vpwptjk = vwtjk - vtjk*wtjk
    upvptij = uvtij - utij*vtij
    if (ltempeq) then
      wpthlptk = wthltk - wmt*thltk
    end if
    if (nsv>0) then
      wpsv1p = wsv1tk - wmt*sv1tk
    end if
    if (nsv>1) then
      wpsv2p = wsv2tk - wmt*sv2tk
    end if
    if (nsv>2) then
      wpsv3p = wsv3tk - wmt*sv3tk
    end if
    if (nsv>3) then
      wpsv4p = wsv4tk - wmt*sv4tk
    end if

    !> Variances and TKE
    if (ltempeq) then
      thlpthlpt = thlthlt - thlt*thlt
    end if
    upuptc = uutc - utc*utc
    vpvptc = vvtc - vtc*vtc
    wpwptc = wwtc - wtc*wtc
    tketc = 0.5*(upuptc + vpvptc + wpwptc)
    if (nsv>0) then
      sv1psv1pt = sv1sv1t - sv1t*sv1t
    end if
    if (nsv>1) then
      sv2psv2pt = sv2sv2t - sv2t*sv2t
    end if
    if (nsv>2) then
      sv3psv3pt = sv3sv3t - sv3t*sv3t
    end if
    if (nsv>3) then
      sv4psv4pt = sv4sv4t - sv4t*sv4t
    end if

!      if (myid == 0) then
          allocate(varst(imax,jmax,khigh-klow+1,nstatt))
          call writestat_nc(ncidt,1,tncstatt,(/timee/),nrect,.true.)
          varst(:,:,:,1)  = umt(ib:ie,jb:je,kb:ke)
          varst(:,:,:,2)  = vmt(ib:ie,jb:je,kb:ke)
          varst(:,:,:,3)  = wmt(ib:ie,jb:je,kb:ke)
          varst(:,:,:,4)  = thlt(ib:ie,jb:je,kb:ke)
          varst(:,:,:,5)  = qtt(ib:ie,jb:je,kb:ke)
          varst(:,:,:,6)  = pt(ib:ie,jb:je,kb:ke)
          varst(:,:,:,7)  = sv1t(ib:ie,jb:je,kb:ke)
          varst(:,:,:,8)  = sv2t(ib:ie,jb:je,kb:ke)
          varst(:,:,:,9)  = sv3t(ib:ie,jb:je,kb:ke)
          varst(:,:,:,10) = sv4t(ib:ie,jb:je,kb:ke)
          varst(:,:,:,11) = PSSt(ib:ie,jb:je,kb:ke)

          varst(:,:,:,12) = upwptik(ib:ie,jb:je,kb:ke)
          varst(:,:,:,13) = vpwptjk(ib:ie,jb:je,kb:ke)
          varst(:,:,:,14) = upvptij(ib:ie,jb:je,kb:ke)
          varst(:,:,:,15) = wpthlptk(ib:ie,jb:je,kb:ke)
          varst(:,:,:,16) = wpsv1p(ib:ie,jb:je,kb:ke)
          varst(:,:,:,17) = wpsv2p(ib:ie,jb:je,kb:ke)
          varst(:,:,:,18) = wpsv3p(ib:ie,jb:je,kb:ke)
          varst(:,:,:,19) = wpsv4p(ib:ie,jb:je,kb:ke)

          varst(:,:,:,20) = thlpthlpt(ib:ie,jb:je,kb:ke)
          varst(:,:,:,21) = upuptc(ib:ie,jb:je,kb:ke)
          varst(:,:,:,22) = vpvptc(ib:ie,jb:je,kb:ke)
          varst(:,:,:,23) = wpwptc(ib:ie,jb:je,kb:ke)
          varst(:,:,:,24) = tketc(ib:ie,jb:je,kb:ke)
          varst(:,:,:,25) = sv1psv1pt(ib:ie,jb:je,kb:ke)
          varst(:,:,:,26) = sv2psv2pt(ib:ie,jb:je,kb:ke)
          varst(:,:,:,27) = sv3psv3pt(ib:ie,jb:je,kb:ke)
          varst(:,:,:,28) = sv4psv4pt(ib:ie,jb:je,kb:ke)

          varst(:,:,:,29) = sv1sgst(ib:ie,jb:je,kb:ke)
          varst(:,:,:,30) = sv2sgst(ib:ie,jb:je,kb:ke)
          varst(:,:,:,31) = sv3sgst(ib:ie,jb:je,kb:ke)
          varst(:,:,:,32) = sv4sgst(ib:ie,jb:je,kb:ke)

          ! varst(:,:,:,33) = sv1max(ib:ie,jb:je,kb:ke)
          ! varst(:,:,:,34) = sv2max(ib:ie,jb:je,kb:ke)
          ! varst(:,:,:,35) = sv3max(ib:ie,jb:je,kb:ke)
          ! varst(:,:,:,36) = sv4max(ib:ie,jb:je,kb:ke)

          ! do n=1,nstatt
          !   do k=kb,ke
          !     do j=jb,je
          !       do i=ib,ie
          !         if (.not.( varst(i,j,k,n) >= -huge(1.0) .and. varst(i,j,k,n) <= huge(1.0) )) then
          !           write(*,*) '1. Invalid value at myid=', myid, ', myidx=', myidx, ', myidy=', myidy, ', i=', i, ', j=', j, ', k=', k, ', n=', n, ', val=', varst(i,j,k,n)
          !         elseif (.not. ieee_is_finite( varst(i,j,k,n) )) then
          !           write(*,*) '2. Invalid value at myid=', myid, ', myidx=', myidx, ', myidy=', myidy, ', i=', i, ', j=', j, ', k=', k, ', n=', n, ', val=', varst(i,j,k,n)
          !         end if
          !       end do
          !     end do
          !   end do
          ! end do

          call writestat_nc(ncidt,nstatt,ncstatt,varst,nrect,imax,jmax,khigh-klow+1)
!        end if !myid
         deallocate(varst)
      end if !ltdump

      if (lmintdump) then
     !      if (myid == 0) then
            allocate(varsmint(imax,jmax,khigh-klow+1,nstatmint))
            call writestat_nc(ncidmint,1,tncstatmint,(/timee/),nrecmint,.true.)
            varsmint(:,:,:,1)  = umt(ib:ie,jb:je,kb:ke)
            varsmint(:,:,:,2)  = vmt(ib:ie,jb:je,kb:ke)
            varsmint(:,:,:,3)  = wmt(ib:ie,jb:je,kb:ke)
            varsmint(:,:,:,4)  = thlt(ib:ie,jb:je,kb:ke)
            varsmint(:,:,:,5)  = qtt(ib:ie,jb:je,kb:ke)
            varsmint(:,:,:,6)  = pt(ib:ie,jb:je,kb:ke)

            call writestat_nc(ncidmint,nstatmint,ncstatmint,varsmint,nrecmint,imax,jmax,khigh-klow+1)
     !        end if !myid
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

  deallocate(thlk,qtk,uik,wik,vjk,wjk,uij,vij,uc,vc,wc)
  deallocate(thlsgs,qtsgs,usgs,vsgs,wsgs)
  deallocate(sv1k,sv2k,sv3k,sv4k,sv1sgs,sv2sgs,sv3sgs,sv4sgs,PSS,wpsv1p,wpsv2p,wpsv3p,wpsv4p,sv1psv1pt,sv2psv2pt,sv3psv3pt,sv4psv4pt)
  deallocate(upwptik,vpwptjk,upvptij,wpthlptk,thlpthlpt,upuptc,vpvptc,wpwptc,tketc)

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
  use modsubgrid,       only : ekh
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
      use modglobal, only  : ltdump
    implicit none

!      if (ltkedump) then
!        call exitstat_nc(ncidtke)
!      endif

!       if (ltdump) then
!         call exitstat_nc(ncidt)
!       endif

  end subroutine exitstatsdump


end module modstatsdump
