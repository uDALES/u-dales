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

  use modglobal, only : longint,nvar
  use modfields, only : ncname
  implicit none
  private
  PUBLIC :: initfielddump, fielddump,exitfielddump,tec3d,tec2d
  save
  !NetCDF variables
  integer :: ncid,nrec = 0
  character(80) :: fname = 'fielddump.xxx.xxx.nc'
  !dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav,tnext
  integer :: klow,khigh
  logical :: lfielddump= .true. !< switch to enable the fielddump (on/off)
  logical :: ldiracc   = .false. !< switch for doing direct access writing (on/off)
  logical :: lbinary   = .false. !
 

contains
 !> Initializing fielddump. Read out the namelist, initializing the variables
  subroutine initfielddump
    use modmpi,   only :myid,my_real,mpierr,comm3d,mpi_logical,mpi_integer,cmyid
    use modglobal,only :imax,jmax,kmax,cexpnr,ifnamopt,fname_options,dtmax,dtav_glob,kb,ke, ladaptive,dt_lim,btime,lnetcdf,nsv
    use modstat_nc,only : open_nc, define_nc,ncinfo,writestat_dims_nc
    implicit none
    integer :: ierr


    namelist/NAMFIELDDUMP/ &
         dtav,lfielddump,ldiracc,lbinary,klow,khigh

    dtav=dtav_glob
    klow=kb
    khigh=ke
    if(myid==0)then
       open(ifnamopt,file=fname_options,status='old',iostat=ierr)
       read (ifnamopt,NAMFIELDDUMP,iostat=ierr)
       if (ierr > 0) then
          print *, 'Problem in namoptions NAMFIELDDUMP'
          print *, 'iostat error: ', ierr
          stop 'ERROR: Problem in namoptions NAMFIELDDUMP'
       endif
       write(6 ,NAMFIELDDUMP)
       close(ifnamopt)
    end if
    call MPI_BCAST(klow        ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(khigh       ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(dtav        ,1,MY_REAL    ,0,comm3d,ierr)
    call MPI_BCAST(lfielddump  ,1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(ldiracc     ,1,MPI_LOGICAL,0,comm3d,ierr)
    call MPI_BCAST(lbinary     ,1,MPI_LOGICAL,0,comm3d,ierr)
    
    tnext      = dtav   +btime
    if(.not.(lfielddump)) return
    !    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
       stop 'dtav should be a integer multiple of dtmax'
    end if
      if (lnetcdf) then
        fname(11:13) = cmyid
        fname(15:17) = cexpnr
        call ncinfo(tncname(1,:),'time','Time','s','time')
        call ncinfo(ncname( 1,:),'u','West-East velocity','m/s','mttt')
        call ncinfo(ncname( 2,:),'v','South-North velocity','m/s','tmtt')
        call ncinfo(ncname( 3,:),'w','Vertical velocity','m/s','ttmt')
        call ncinfo(ncname( 4,:),'qt','Total water mixing ratio','1e-5kg/kg','tttt')
        call ncinfo(ncname( 5,:),'ql','Liquid water mixing ratio','1e-5kg/kg','tttt')
        call ncinfo(ncname( 6,:),'thl','Liquid water potential temperature','K','tttt')
!        call ncinfo(ncname( 7,:),'ekm','subgrid eddy visucosity (for momentum)','m2/s','tttt')
!        call ncinfo(ncname( 8,:),'sbshr','sub grid shear','m2/s2','tttt')
       
        if (nsv>0) call ncinfo(ncname( 7,:),'sca1','scalar 1','M','tttt') !unit?
        if (nsv>1) call ncinfo(ncname( 8,:),'sca2','scalar 2','M','tttt')
        if (nsv>2) call ncinfo(ncname( 9,:),'sca3','scalar 3','M','tttt')
        if (nsv>3) call ncinfo(ncname(10,:),'sca4','scalar 4','M','tttt')
      
        call open_nc( fname, ncid, nrec, n1=imax, n2=jmax, n3=khigh-klow+1)
        if (nrec==0) then
          call define_nc( ncid, 1, tncname)
          call writestat_dims_nc(ncid)  
        end if
       call define_nc( ncid, nvar, ncname)
      end if

  end subroutine initfielddump

  !> Do fielddump. Collect data to truncated (2 byte) integers, and write them to file
  subroutine fielddump
    use modfields, only : u0,v0,w0,thl0,qt0,ql0,sv0  !ILS13 21.04.2015 changed to u0 from um  etc
    use modsurfdata,only : thls,qts,thvs
    use modglobal, only : ib,ie,ih,jb,je,jh,ke,kb,kh,rk3step,timee,dt_lim,cexpnr,ifoutput,imax,jmax, lnetcdf, tfielddump, tnextfielddump,nsv
    !use modmpi,    only : myid,cmyid
    !use modsubgriddata, only : ekm,sbshr
    use modstat_nc, only : writestat_nc
    use modmpi, only : myid,cmyid
    implicit none

   ! integer(KIND=selected_int_kind(6)), allocatable :: field(:,:,:),vars(:,:,:,:)
    real, allocatable :: field(:,:,:), vars(:,:,:,:)
    integer i,j,k
    integer :: writecounter = 1
    integer :: reclength
 !write(*,*) "timee",timee
 !write(*,*) "tnextfielddump",tnextfielddump
 !write(*,*) "lfielddump",lfielddump
 
 if (.not. (timee>=tnextfielddump)) return 

    if (.not. lfielddump) return
    if (rk3step/=3) return
    !if(timee<tnext) then
       !dt_lim = min(dt_lim,tnext-timee)
    !   return
    !end if
    tnextfielddump=tnextfielddump+tfielddump
    !tnext = tnext+dtav
    !dt_lim = minval((/dt_lim,tnext-timee/))

   ! allocate(field(ib-ih:ie+ih,jb-jh:je+jh,ke+1))
   ! I changed this, ILS13 
     allocate(field(ib-ih:ie+ih,jb-jh:je+jh,kb-kh:ke+kh))

   allocate(vars(imax,jmax,khigh-klow+1,nvar))

    reclength = imax*jmax*(khigh-klow+1)*2
    !write(*,*) 'um',1.0E2*um
    !field = NINT(1.0E2*um,selected_int_kind(4))
    field = u0
     !field(ib:ie,jb:je,klow)=-999
       if (lnetcdf) vars(:,:,:,1) = field(ib:ie,jb:je,klow:khigh) !need to change this as above?
     if (lbinary) then
     if (ldiracc) then
       open (ifoutput,file='wbuu.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
       write (ifoutput, rec=writecounter) field(ib:ie,jb:je,klow:khigh)
    else
       open  (ifoutput,file='wbuu.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
       write (ifoutput) (((field(i,j,k),i=ib,ie),j=jb,je),k=klow,khigh)
    end if
    close (ifoutput)
end if

    !field = NINT(1.0E2*vm,selected_int_kind(4))
    field = v0
    ! field(ib:ie,jb:je,klow)=-999
   if (lnetcdf) vars(:,:,:,2) = field(ib:ie,jb:je,klow:khigh)
if (lbinary) then    
if (ldiracc) then
       open (ifoutput,file='wbvv.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
       write (ifoutput, rec=writecounter) field(ib:ie,jb:je,klow:khigh)
    else
       open  (ifoutput,file='wbvv.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
       write (ifoutput) (((field(i,j,k),i=ib,ie),j=jb,je),k=klow,khigh)
    end if
    close (ifoutput)
end if

    !field = NINT(1.0E3*wm,selected_int_kind(4))
    field = w0
     !  field(ib:ie,jb:je,klow)=-999
  if (lnetcdf) vars(:,:,:,3) = field(ib:ie,jb:je,klow:khigh)
if (lbinary) then
    if (ldiracc) then
       open (ifoutput,file='wbww.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
       write (ifoutput, rec=writecounter) field(ib:ie,jb:je,klow:khigh)
    else
       open  (ifoutput,file='wbww.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
       write (ifoutput) (((field(i,j,k),i=ib,ie),j=jb,je),k=klow,khigh)
    end if
    close (ifoutput)
   end if

     !  field(ib:ie,jb:je,klow+1:khigh) = NINT(1.0E5*qtm(ib:ie,jb:je,klow+1:khigh),2)
     field(ib:ie,jb:je,klow+1:khigh) = qt0(ib:ie,jb:je,klow+1:khigh)
     !  field(ib:ie,jb:je,klow)=-999
  if (lnetcdf) vars(:,:,:,4) = field(ib:ie,jb:je,klow:khigh)
if (lbinary) then
    if (ldiracc) then
       open (ifoutput,file='wbqt.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
       write (ifoutput, rec=writecounter) field(ib:ie,jb:je,klow:khigh)
    else
       open  (ifoutput,file='wbqt.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
       write (ifoutput) (((field(i,j,k),i=ib,ie),j=jb,je),k=klow,khigh)
    end if
    close (ifoutput)
end if

    !field = NINT(ql0,2)
    field(ib:ie,jb:je,klow+1:khigh) = ql0(ib:ie,jb:je,klow+1:khigh)
    !  field(ib:ie,jb:je,klow)=-999
  if (lnetcdf) vars(:,:,:,5) = field(ib:ie,jb:je,klow:khigh)
if (lbinary) then    
if (ldiracc) then
       open (ifoutput,file='wbql.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
       write (ifoutput, rec=writecounter) field(ib:ie,jb:je,klow:khigh)
    else
       open  (ifoutput,file='wbql.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
       write (ifoutput) (((field(i,j,k),i=ib,ie),j=jb,je),k=klow,khigh)
    end if
    close (ifoutput)
end if

    !field = NINT(1.0E1*thlm,2)
    !field(ib:ie,jb:je,klow+1:khigh) = thlm(ib:ie,jb:je,klow+1:khigh)                    !thl0(ib:ie,jb:je,klow+1:khigh) !NINT(1.0E1*thlm(ib:ie,jb:je,klow+1:khigh),2)
     field = thl0-273.15
     write(*,*) 'fieldmax',MAXVAL(field)
    ! field(ib:ie,jb:je,klow)=-999
  if (lnetcdf) vars(:,:,:,6) = field(ib:ie,jb:je,klow:khigh)

if (lbinary) then
    if (ldiracc) then
       open (ifoutput,file='wbtl.'//cmyid//'.'//cexpnr,access='direct', form='unformatted', recl=reclength)
       write (ifoutput, rec=writecounter) field(ib:ie,jb:je,klow:khigh)
    else
       open  (ifoutput,file='wbtl.'//cmyid//'.'//cexpnr,form='unformatted',position='append')
       write (ifoutput) (((field(i,j,k),i=ib,ie),j=jb,je),k=klow,khigh)
    end if
    close (ifoutput)
end if

    if (lnetcdf) then        
    
     if (nsv>0) then 
     field = sv0(:,:,:,1)    !maybe sv0 ?
     vars(:,:,:,7) = field(ib:ie,jb:je,klow:khigh)
     end if
     if (nsv>1) then 
     field = sv0(:,:,:,2)    !maybe sv0 ?
     vars(:,:,:,8) = field(ib:ie,jb:je,klow:khigh)
     end if

     if (nsv>2) then 
     field = sv0(:,:,:,3)    !maybe sv0 ?
     vars(:,:,:,9) = field(ib:ie,jb:je,klow:khigh)
     end if

     if (nsv>3) then 
     field = sv0(:,:,:,4)    !maybe sv0 ?
     vars(:,:,:,10) = field(ib:ie,jb:je,klow:khigh)
     end if
end if


    !field = ekm
    !if (lnetcdf) vars(:,:,:,7) = field(ib:ie,jb:je,klow:khigh)

    !field = sbshr
    !if (lnetcdf) vars(:,:,:,8) = field(ib:ie,jb:je,klow:khigh)


     if(lnetcdf) then
           write(*,*) "write tnc"
           call writestat_nc(ncid,1,tncname,(/timee/),nrec,.true.)
           write(*,*) "write nc"
           write(*,*) "nvar", nvar
           call writestat_nc(ncid,nvar,ncname,vars,nrec,imax,jmax,khigh-klow+1)
         end if
        if(lbinary) then
         if (myid==0) then
           open(ifoutput, file='wbthls.'//cexpnr,form='formatted',position='append')
           write(ifoutput,'(F12.1 3F12.5)') timee,thls,qts,thvs
           close(ifoutput)
         end if
    writecounter=writecounter+1
end if
    deallocate(field,vars)

  end subroutine fielddump
  !> Clean up when leaving the run
  subroutine exitfielddump
      use modstat_nc, only : exitstat_nc
      use modglobal, only  : lnetcdf
    implicit none

       if(lfielddump .and. lnetcdf) call exitstat_nc(ncid)
  end subroutine exitfielddump

  subroutine tec2d
    use modfields, only    : um,vm,wm,thlm,qtm,ql0,pres0,&
                             thl0,thl0h,svm,e120,uav,vav,wav,thlav,svav,uuav,vvav,wwav,uvav,uwav,vwav,&
                             sv2av,thl2av,yplus,umint,vmint,wmint,upupav,vpvpav,wpwpav,upvpav,upwpav,vpwpav
    use modglobal, only    : ib,ie,jb,je,kb,ke,xf,dy,zf,jtot,jmax,dzh,zh,numol,nsv
    use modmpi, only       : cmyid,myid
    use modsubgriddata,only: ekm
    implicit none

    integer  ::  i,j,k,ksave
!    real     ::  upupav,vpvpav,wpwpav,upvpav,upwpav,vpwpav,thlpthlpav,svpsvpav
!    ! this is overline(u'u'), overline(u'v'),etc.
    real     ::  thlpthlpav  ! this is overline(u'u'), overline(u'v'),etc.

        open(unit=10,file='tec2D'//cmyid//'.tec')

        do j=jb,je
        do i=ib,ie
        do k=kb,ke
          umint(i,j,k) = 0.5 * (um(i+1,j,k) + um(i,j,k))
          vmint(i,j,k) = 0.5 * (vm(i,j+1,k) + vm(i,j,k))
          wmint(i,j,k) = 0.5 * (wm(i,j,k+1) + wm(i,j,k))
        end do
        end do
        end do


    if(myid==0)then
        write (10,*) 'title= "domain 2D"'
        write (10,*) 'variables = "z","x","y",'
        write (10,*) '"u","v","w","viscratio","thl"'!,"sca1","sca2","sca3","sca4","sca5"'
        if (nsv>0)  write (10,*) '"sca1"'
        if (nsv>1)  write (10,*) '"sca2"'
        if (nsv>2)  write (10,*) '"sca3"'
        if (nsv>3)  write (10,*) '"sca4"'
        if (nsv>4)  write (10,*) '"sca5"'
        write (10,*) 'zone t= "zone-title", f=point'

        write(10,*)'i = ', 1, ' , j = ', ie-ib+1 , ' , k = ', jtot
    endif

!        ksave = kb+19
        do j=jb,je
        do i=ib,ie
        do k=kb+19,kb+19
           write(10,3001) zf(k),xf(i),(j-jb+myid*jmax)*dy+0.5*dy,&
                          umint(i,j,k),vmint(i,j,k),wmint(i,j,k),&
                          (ekm(i,j,k)-numol)/numol,thlm(i,j,k)!,svm(i,j,k,1),svm(i,j,k,2),svm(i,j,k,3),svm(i,j,k,4),svm(i,j,k,5)!,svm(i,j,k,6)!,svm(i,j,k,7),&
           if (nsv>0) write(10,3002) svm(i,j,k,1)
           if (nsv>1) write(10,3003) svm(i,j,k,2)
           if (nsv>2) write(10,3004) svm(i,j,k,3)
           if (nsv>3) write(10,3005) svm(i,j,k,4)
           if (nsv>4) write(10,3006) svm(i,j,k,5)
        enddo
        enddo
        enddo
        close(10)
3001    format (8e20.12)
3002    format (1e20.12)
3003    format (1e20.12)
3004    format (1e20.12)
3005    format (1e20.12)
3006    format (1e20.12)


  end subroutine tec2d


  subroutine tec3d
    use modfields, only    : um,vm,wm,thlm,qtm,ql0,pres0,&
                             thl0,thl0h,svm,e120,uav,vav,wav,thlav,svav,uuav,vvav,wwav,uvav,uwav,vwav,&
                             sv2av,thl2av,yplus,umint,vmint,wmint,upupav,vpvpav,wpwpav,upvpav,upwpav,vpwpav
    use modglobal, only    : ib,ie,jb,je,kb,ke,xf,dy,zf,jtot,jmax,dzh,zh,numol,nsv
    use modmpi, only       : cmyid,myid
    use modsubgriddata,only: ekm
    implicit none

    integer  ::  i,j,k
    !    real     ::  upupav,vpvpav,wpwpav,upvpav,upwpav,vpwpav,thlpthlpav,svpsvpav  ! this is overline(u'u'), overline(u'v'),etc.
    real     ::  thlpthlpav  ! this is overline(u'u'), overline(u'v'),etc.

    open(unit=10,file='tec3D'//cmyid//'.tec')

    do j=jb,je
       do i=ib,ie
          do k=kb,ke
             umint(i,j,k) = 0.5 * (um(i+1,j,k) + um(i,j,k))
             vmint(i,j,k) = 0.5 * (vm(i,j+1,k) + vm(i,j,k))
             wmint(i,j,k) = 0.5 * (wm(i,j,k+1) + wm(i,j,k))
          end do
       end do
    end do


    if(myid==0)then
!        write (10,*) 'title= "domain 3"'
!        write (10,*) 'variables = "z","x","y",'
!        write (10,*) '"u","v","w","p","viscratio","nu_tot","nu_t","theta_l","scalar","e120",'
!        write (10,*) '"U","V","W","theta_l_av","scalar_av","u''u''","v''v''",'
!        write (10,*) '"w''w''","u''v''","u''w''","v''w''","theta_l''_sqr","scalar''_sqr","yplus"'
!        write (10,*) 'zone t= "zone-title", f=point' 
        write (10,*) 'title= "domain 3"'
        write (10,*) 'variables = "z","x","y",'
!        write (10,*) '"u","v","w","scalar1","scalar2","scalar3","scalar4","scalar5","viscratio"'
!        write (10,*) '"u","v","w","viscratio","thl","sca1","sca2","sca3","sca4"'
!        write (10,*) '"u","v","w","viscratio","thl"'!,"sca1","sca2","sca3","sca4","sca5","sca6"'!,"sca7","sca8","sca9"'
        write (10,*) '"u","v","w","viscratio","thl","qtm","ql0"'!,"sca1","sca2"'!,"sca3","sca4","sca5"' 
        if (nsv>0)  write (10,*) '"sca1"'
        if (nsv>1)  write (10,*) '"sca2"'
        if (nsv>2)  write (10,*) '"sca3"'
        if (nsv>3)  write (10,*) '"sca4"'
        if (nsv>4)  write (10,*) '"sca5"'
        write (10,*) 'zone t= "zone-title", f=point' 

       write(10,*)'i = ', ke-kb+1, ' , j = ', ie-ib+1 , ' , k = ', jtot

    endif

        do j=jb,je
        do i=ib,ie
        do k=kb,ke
!           upupav = uuav(i,j,k) - uav(i,j,k)**2             ! overline(u'u') = overline(uu) - U^2
!           vpvpav = vvav(i,j,k) - vav(i,j,k)**2             ! overline(v'v') = overline(vv) - V^2
!           wpwpav = wwav(i,j,k) - wav(i,j,k)**2             ! overline(w'w') = overline(ww) - W^2
!           upvpav = uvav(i,j,k) - uav(i,j,k)*vav(i,j,k)     ! overline(u'v') = overline(uv) - U*V
!           upwpav = uwav(i,j,k) - uav(i,j,k)*wav(i,j,k)     ! overline(u'w') = overline(uw) - U*W
!           vpwpav = vwav(i,j,k) - vav(i,j,k)*wav(i,j,k)     ! overline(v'w') = overline(vw) - V*W

!           thlpthlpav = thl2av(i,j,k) - thlav(i,j,k)**2     ! overline(thl'thl') = overline(thl^2) - THL^2
!           svpsvpav = sv2av(i,j,k,1) - sv2av(i,j,k,1)**2    ! overline(sv'sv') = overline(sv^2) - SV^2
           
!           write(10,3000) zf(k),xf(i),(j-jb+myid*jmax)*dy+0.5*dy,umint(i,j,k),vmint(i,j,k),wmint(i,j,k),pres0(i,j,k),&
!           (ekm(i,j,k)-numol)/numol,ekm(i,j,k),ekm(i,j,k)-numol,thl0(i,j,k),svm(i,j,k,1),e120(i,j,k),uav(i,j,k),&
!           vav(i,j,k),wav(i,j,k),thlav(i,j,k),svav(i,j,k,1),upupav,vpvpav,wpwpav,upvpav,&
!           upwpav,vpwpav,thlpthlpav,svpsvpav,yplus(i,j,k)
           write(10,3000) zf(k),xf(i),(j-jb+myid*jmax)*dy+0.5*dy,umint(i,j,k),vmint(i,j,k),wmint(i,j,k),(ekm(i,j,k)-numol)/numol,&
                          thlm(i,j,k),qtm(i,j,k),ql0(i,j,k)!,svm(i,j,k,1),svm(i,j,k,2)!,svm(i,j,k,3),svm(i,j,k,4),svm(i,j,k,5)!,svm(i,j,k,6)!,svm(i,j,k,7),&
           if (nsv>0) write(10,3002) svm(i,j,k,1)
           if (nsv>1) write(10,3003) svm(i,j,k,2)
           if (nsv>2) write(10,3004) svm(i,j,k,3)
           if (nsv>3) write(10,3005) svm(i,j,k,4)
           if (nsv>4) write(10,3006) svm(i,j,k,5)
!                          svm(i,j,k,8),svm(i,j,k,9)
!           write(10,3000) zf(k),xf(i),(j-jb+myid*jmax)*dy+0.5*dy,umint(i,j,k),vmint(i,j,k),wmint(i,j,k),&
!                          svm(i,j,k,1),svm(i,j,k,2),svm(i,j,k,3),svm(i,j,k,4),svm(i,j,k,5),(ekm(i,j,k)-numol)/numol 
!                          (ekm(i,j,k)-numol)/numol 
        enddo
        enddo
        enddo
        close(10)
3000    format (8e20.12)
3002    format (1e20.12)
3003    format (1e20.12)
3004    format (1e20.12)
3005    format (1e20.12)
3006    format (1e20.12)

  end subroutine tec3d

end module modfielddump
