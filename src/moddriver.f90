!! moddriverinlet.f90 contains the set up for creating a driver simulation to use an inlet conditions for a following simulation/
!! This creates a realistic turbulent inlet profiles as necessary.
!! The fields are extracted from a specified plane and stored in
!! respective files. These files are then used to create inlet
!! conditions for a following simualtion. Fields are linearly
!! interpolated in time where necessary.
!! Code set up is adapted from modinlet.f90 since similar processes are involved

!! \author Anton Esmail-Yakas, Imperial College London, August 5th 2017.
!! Edited by tg3315, ICL, May 2019.
!! \todo Documentation
!!       Remove unecessary "use" variables
!!       Remove unecessary commented lines
module moddriver
use modinletdata
implicit none
save
  public :: initdriver,exitdriver,readdriverfile,drivergen

contains
  subroutine initdriver
    use modglobal, only : ih,ib,ie,jh,jb,je,kb,ke,kh,jhc,khc,idriver,iplane,xf,lstoreplane,nstore,Uinf,ltempeq,lmoist,pi,zf,zh,driverstore,tdriverstart,tdriverdump,timeleft,dtdriver,nsv,timee
    use modfields, only : um
    use modmpi, only : myid,nprocs

    implicit none
    real    :: pfi, epsi
    integer :: k

    if (idriver==1) then
      if (tdriverstart < timee) then
        write(0, *) 'ERROR: tdriverstart must be greater than the elapsed time at the start of the simulation'
        stop 1
      end if
      tdriverdump = tdriverstart
   endif

    if (idriver==1) then
      allocate(storetdriver(1:driverstore))
      allocate(storeu0driver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      allocate(storeumdriver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      allocate(storev0driver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      allocate(storevmdriver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      allocate(storew0driver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      allocate(storewmdriver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      !allocate(storee120driver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      !allocate(storee12mdriver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      if (ltempeq  ) then
        allocate(storethl0driver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
        allocate(storethlmdriver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      end if
      if (lmoist  ) then
        allocate(storeqt0driver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
        allocate(storeqtmdriver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      end if
      if (nsv>0  ) then
        allocate(storesv0driver(jb-jhc:je+jhc,kb-khc:ke+khc,1:nsv,1:driverstore))
        allocate(storesvmdriver(jb-jhc:je+jhc,kb-khc:ke+khc,1:nsv,1:driverstore))
      end if

      irecydriver = iplane!+ib        ! index of recycle plane from driver simulation equals iplane (read from namoptions)

    else if (idriver == 2) then

      allocate(storetdriver(1:driverstore))
      allocate(storeu0driver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      allocate(storeumdriver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      allocate(storev0driver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      allocate(storevmdriver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      allocate(storew0driver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      allocate(storewmdriver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      !allocate(storee120driver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      !allocate(storee12mdriver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
      allocate(u0driver(jb-jh:je+jh,kb-kh:ke+kh))
      allocate(v0driver(jb-jh:je+jh,kb-kh:ke+kh))
      allocate(w0driver(jb-jh:je+jh,kb-kh:ke+kh))
      !allocate(e120driver(jb-jh:je+jh,kb-kh:ke+kh))
      allocate(umdriver(jb-jh:je+jh,kb-kh:ke+kh))
      allocate(vmdriver(jb-jh:je+jh,kb-kh:ke+kh))
      allocate(wmdriver(jb-jh:je+jh,kb-kh:ke+kh))
      !allocate(e12mdriver(jb-jh:je+jh,kb-kh:ke+kh))

      if (ltempeq ) then
        allocate(storethl0driver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
        allocate(storethlmdriver(jb-jh:je+kh,kb-kh:ke+kh,1:driverstore))
        allocate(thl0driver(jb-jh:je+jh,kb-kh:ke+kh))
        allocate(thlmdriver(jb-jh:je+jh,kb-kh:ke+kh))
      end if
      if (lmoist ) then
        allocate(storeqt0driver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
        allocate(storeqtmdriver(jb-jh:je+jh,kb-kh:ke+kh,1:driverstore))
        allocate(qt0driver(jb-jh:je+jh,kb-kh:ke+kh))
        allocate(qtmdriver(jb-jh:je+jh,kb-kh:ke+kh))
      end if
      if (nsv>0 ) then
        allocate(storesv0driver(jb-jhc:je+jhc,kb-khc:ke+khc,1:nsv,1:driverstore))
        allocate(storesvmdriver(jb-jhc:je+jhc,kb-khc:ke+khc,1:nsv,1:driverstore))
        allocate(sv0driver(jb-jhc:je+jhc,kb-khc:ke+khc,1:nsv))
        allocate(svmdriver(jb-jhc:je+jhc,kb-khc:ke+khc,1:nsv))
      end if
      irecydriver = iplane!+ib
    else
      return
    end if

  end subroutine initdriver

  subroutine drivergen
    use modglobal,   only : ib,ie,jb,je,jgb,jge,kb,ke,zf,zh,dzf,dzhi,timee,btime,totavtime,rk3step,&
                            dt,numol,iplane,lles,idriver,inletav,runavtime,Uinf,lwallfunc,linletRA,&
                            totinletav,lstoreplane,nstore,driverstore,prandtlmoli,numol,grav,lbuoyancy,&
                            lfixinlet,lfixutauin,tdriverstart,dtdriver,tdriverdump,ltempeq,lmoist,nsv
    use modfields,   only : u0,v0,w0,e120,thl0,qt0,wm,uprof
    use modsave,     only : writerestartfiles
    use modmpi,      only : slabsum,myid
    implicit none

    real :: inlrec                              ! time of last inlet record
    real :: elapsrec                            ! time elapsed in this inlet record  
    real :: dtint                               ! dt for linear interpolation
    REAL*8, PARAMETER :: eps = 1d-3
    integer i,j,k,kk,kdamp,x

    if (idriver == 1) then

      ! if (.not. (rk3step==3)) return
      if (.not. (timee>=tdriverstart)) return
      if (.not. (timee>=tdriverdump)) return
      if (nstepreaddriver>=driverstore) return
      if (nstepreaddriver==0) then
        ! tdriverdump = timee
        tdriverdump = tdriverstart
        ! tdriverstart = timee   !Update tdriverstart to the actual recorded value
        if ((myid==0) .and. (rk3step==3)) then
          write(6,*) '=================================================================='
          write(6,*) '*** Starting to write data for driver simulation ***'
          write(6,*) 'Driver recording variables:'
          write(6,'(A,F9.2,A,I4,A,F12.9)') ' Starting time: ',tdriverdump,' Stored time steps: ',driverstore,'     Inlet record intervals: ',dtdriver
          write(6,*) '=================================================================='
        end if
      end if
      if(rk3step==3) then
        nstepreaddriver = nstepreaddriver + 1
        tdriverdump = tdriverdump + dtdriver
        ! storetinlet(nstepreaddriver) = timee - tdriverstart
        call writedriverfile
      end if
   
    elseif (idriver == 2) then

      ! if (.not. rk3step==1) return
      if (timee>maxval(storetdriver)) then
        if(myid==0) then
          write(0,'(A,F9.2,A,F9.2)') 'timee: ',timee,'     Final inlet driver time:',maxval(storetdriver)
          write(0,'(A,I4,A,I4)') 'Inlet driver step: ',nstepreaddriver,'     Total inlet driver steps:',driverstore
        end if
        stop 'Time in simulation has exceeded the inlet information - no more inlet data available!'
      end if
      x = minloc(abs(storetdriver-timee),1)
      elapsrec = storetdriver(x) - timee
      if(myid==0) then
        ! if(rk3step==1) then
        ! write(6,*) '============ Inlet interpolating ============='
        ! write(6,*) 'Inlet interpolation time = ', elapsrec
        ! write(6,'(A,F9.4)') 'Inlet driver time stamp (x)  = ', storetdriver(x)
        ! write(6,'(A,F9.4)') 'Inlet driver time stamp (x+1) = ', storetdriver(x+1)
        ! write(6,'(A,F9.4)') 'Inlet driver time stamp (x-1) = ', storetdriver(x-1)
        ! write(6,'(A,E20.12)') 'Reading driver velocity: storeu0driver(je,ke,x) = ', storeu0driver(je,ke,x)
        ! write(6,*) 'Inlet step = ',nstepreaddriver
        ! end if
      end if

      if (abs(elapsrec) < eps) then

        if ((myid==0) .and. (rk3step==1)) then
          write(*,'(A,I5,A,F10.3,A)') '======= Inputs loaded from driver tstep ',x,' (at ',storetdriver(x),'s) ======='
        end if

        u0driver(:,:) = storeu0driver(:,:,x)
        v0driver(:,:) = storev0driver(:,:,x)
        w0driver(:,:) = storew0driver(:,:,x)
        !e120driver(:,:) = storee120driver(:,:,x)
        if (ltempeq) then
          thl0driver(:,:) = storethl0driver(:,:,x)
        end if
        if (lmoist) then
          qt0driver(:,:) = storeqt0driver(:,:,x)
        end if
        if (nsv>0) then
          sv0driver(:,:,:) = storesv0driver(:,:,:,x)
        end if
        nstepreaddriver = x
      elseif ((elapsrec > 0.) .and. (x == 1)) then
     
        if ((myid==0) .and. (rk3step==1)) then
          write(*,'(A,F10.3,A)') '======= Inputs loaded from the proceeding driver tstep 1 (at ',storetdriver(x),'s) ======='
        end if

        u0driver(:,:) = storeu0driver(:,:,x)
        v0driver(:,:) = storev0driver(:,:,x)
        w0driver(:,:) = storew0driver(:,:,x)
        ! e120driver(:,:) = storee120driver(:,:,x)
        if (ltempeq) then
          thl0driver(:,:) = storethl0driver(:,:,x)
        end if
        if (lmoist) then
          qt0driver(:,:) = storeqt0driver(:,:,x)
        end if
        if (nsv>0) then
          sv0driver(:,:,:) = storesv0driver(:,:,:,x)
        end if
        nstepreaddriver = x
      elseif (elapsrec < 0.) then

        if ((myid==0) .and. (rk3step==1)) then
          write(*,'(A,I5,A,F10.3,A,I5,A,F10.3,A)') '======= Inputs interpolated from driver tsteps ',x,' (',storetdriver(x),' s) and ',x+1,' (',storetdriver(x+1),' s) ======='
        end if

        dtint = (timee-storetdriver(x))/(storetdriver(x+1)-storetdriver(x))
        ! if(myid==0) then
        ! write(6,'(A,I4)') 'x: ', x
        ! write(6,'(A,F9.4)') 'dtint: ', dtint
        ! write(6,'(A,E20.12)') 'storeu0driver(1,32,x): ', storeu0driver(1,32,x)
        ! write(6,'(A,E20.12)') 'storeu0driver(1,32,x+1): ', storeu0driver(1,32,x+1)
        ! write(6,'(A,E20.12)') 'u0driver(1,32): ',  storeu0driver(1,32,x) + (storeu0driver(1,32,x+1)-storeu0driver(1,32,x))*dtint
        ! end if
        u0driver(:,:) = storeu0driver(:,:,x) + (storeu0driver(:,:,x+1)-storeu0driver(:,:,x))*dtint
        v0driver(:,:) = storev0driver(:,:,x) + (storev0driver(:,:,x+1)-storev0driver(:,:,x))*dtint
        w0driver(:,:) = storew0driver(:,:,x) + (storew0driver(:,:,x+1)-storew0driver(:,:,x))*dtint
        ! e120driver(:,:) = storee120driver(:,:,x) + (storee120driver(:,:,x+1)-storee120driver(:,:,x))*dtint
        if (ltempeq) then
          thl0driver(:,:) = storethl0driver(:,:,x) + (storethl0driver(:,:,x+1)-storethl0driver(:,:,x))*dtint
        end if
        if (lmoist) then
          qt0driver(:,:) = storeqt0driver(:,:,x) + (storeqt0driver(:,:,x+1)-storeqt0driver(:,:,x))*dtint
        end if
        if (nsv>0) then
          sv0driver(:,:,:) = storesv0driver(:,:,:,x) + (storesv0driver(:,:,:,x+1)-storesv0driver(:,:,:,x))*dtint
        end if
        nstepreaddriver = x
      elseif (elapsrec > 0.) then

        if ((myid==0) .and. (rk3step==1)) then
          write(*,'(A,I5,A,F10.3,A,I5,A,F10.3,A)') '======= Inputs interpolated from driver tsteps ',x,' (', storetdriver(x),' s) and ',x-1,' (',storetdriver(x-1),' s) ======='
        end if

        dtint = (timee-storetdriver(x-1))/(storetdriver(x)-storetdriver(x-1))
        u0driver(:,:) = storeu0driver(:,:,x-1) + (storeu0driver(:,:,x)-storeu0driver(:,:,x-1))*dtint
        v0driver(:,:) = storev0driver(:,:,x-1) + (storev0driver(:,:,x)-storev0driver(:,:,x-1))*dtint
        w0driver(:,:) = storew0driver(:,:,x-1) + (storew0driver(:,:,x)-storew0driver(:,:,x-1))*dtint
        ! e120driver(:,:) = storee120driver(:,:,x-1) + (storee120driver(:,:,x)-storee120driver(:,:,x-1))*dtint
        if (ltempeq) then
          thl0driver(:,:) = storethl0driver(:,:,x-1) + (storethl0driver(:,:,x)-storethl0driver(:,:,x-1))*dtint
        end if
        if (lmoist) then
          qt0driver(:,:) = storeqt0driver(:,:,x-1) + (storeqt0driver(:,:,x)-storeqt0driver(:,:,x-1))*dtint
        end if
        if (nsv>0) then
          sv0driver(:,:,:) = storesv0driver(:,:,:,x-1) + (storesv0driver(:,:,:,x)-storesv0driver(:,:,:,x-1))*dtint
        end if
        nstepreaddriver = x
      end if

      ! if(myid==0) then
        ! write(6,'(A,F9.4)') 'Simulation time: ', timee
        ! write(6,'(A,F9.4)') 'dtint: ', dtint
        ! write(6,*) 'Velocities interpolated:'
        ! write(6,'(A,e20.12)') 'storeu0driver(je,ke,x-1): ', storeu0driver(je,ke,x-1)
        ! write(6,'(A,e20.12)') 'storeu0driver(je,ke,x): ', storeu0driver(je,ke,x)
        ! write(6,'(A,e20.12)') 'storeu0driver(je,ke,x+1): ', storeu0driver(je,ke,x+1)
        ! write(6,'(A,e20.12)') 'Interpolated inlet velocity (jb,20): ', u0driver(jb,20)
        ! write(6,*) 'Temperatures interpolated:'
        ! write(6,'(A,e20.12)') 'storethl0driver(je,20,x-1): ', storethl0driver(jb,20,x-1)
        ! write(6,'(A,e20.12)') 'storethl0driver(je,20,x): ', storethl0driver(jb,20,x)
        ! write(6,'(A,e20.12)') 'storethl0driver(je,20,x+1): ', storethl0driver(jb,20,x+1)
        ! write(6,'(A,e20.12)') 'Interpolated inlet temperature (jb,20): ', thl0driver(jb,20)
      ! end if 

      ! umdriver = u0driver   ! MAYBE ITS BETTER TO WRITE THE M VARIABLES TO FILE TOO AND JUST READ THEM - THOUGH CURRENTLY THIS IS NOT DONE FOR RESTART FILES?? ae1212
      ! vmdriver = v0driver   ! EDIT READ AND WRITE INLET FILES (AND CHECK MODBOUNDARY & MODSURFACE) TO INCLUDE M VARIABLES
      ! wmdriver = w0driver
      ! thlmdriver = thl0driver
      ! qtmdriver = qt0driver

      if (rk3step==1) then
        umdriver = u0driver
        vmdriver = v0driver
        wmdriver = w0driver
        !e12mdriver = e120driver
        if (ltempeq) then
          thlmdriver = thl0driver
        end if
        if (lmoist) then
          qtmdriver = qt0driver
        end if
        if (nsv>0) then
          svmdriver = sv0driver
        end if
      end if

    else

      return

    end if  ! idrivergen

  end subroutine drivergen

  subroutine writedriverfile
    use modglobal, only : timee,tdriverstart,ib,ie,ih,jb,je,jh,kb,ke,kh,cexpnr,ifoutput,nstore,ltempeq,lmoist,driverstore,nsv
    use modfields, only : u0, v0, w0, e120, thl0, qt0, um, sv0
    use modmpi,    only : cmyid,myid
    use modinletdata, only : storetdriver,storeu0driver,storev0driver,storew0driver,storethl0driver,storeqt0driver,&
                             storesv0driver,nfile,nstepreaddriver
    implicit none
    integer :: fileid, IOS
    integer :: i,j,k,n
    integer :: filesizet, filesizev, filesizetest1, filesizetest2, filesizes
    character(15) :: name
    logical :: lexist
    real, allocatable :: arraysizetest(:,:)
    
    allocate(arraysizetest(jb-jh:je+jh,kb-kh:ke+kh))

    inquire(iolength=filesizet)(timee-tdriverstart)
    ! inquire(iolength=filesizetest1)(timee)
    ! inquire(iolength=filesizetest2)u0(1,1,1)
    inquire(iolength=filesizev)u0(irecydriver,:,:)
    inquire(iolength=filesizes)sv0(irecydriver,:,:,:)
    !!  
    ! if((myid==0) .and. (nstepreaddriver==1)) then
      ! write(6,*) 'inquire iolength ', filesizet
      ! write(6,*) 'inquire iolength test', filesizetest1
      ! write(6,*) 'inquire iolength test u', filesizetest2 
    ! end if

    ! inquire(iolength=filesizetest1)arraysizetest(:,:)
    ! filesizetest2 = (je-jb+2*jh)*(ke-kb+2*kh)

    ! if((myid==0) .and. (nstepreaddriver==1)) then
      ! write(6,*) 'je,jb,jh,ke,kb,kh', je,jb,jh,ke,kb,kh
      ! write(6,*) 'inquire iolength test 1 ', filesizetest1
      ! write(6,*) 'inquire iolength test 2', filesizetest2
      ! write(6,*) 'inquire iolength', filesizev
    ! end if

    if(myid==0) then
      write(6,*) '============ Writing driver files ============'
      write(*,*) 'Driver timestep: ', nstepreaddriver
    end if

    if(myid==0) then
      name = 'tdriver_   .'
      name(9:11)= cmyid
      name(13:15)= cexpnr
      ! name(15:18)= '.txt'
      inquire(file=name,exist=lexist)
      if (lexist) then
      ! write(6,*) 'Writing Time stamp to file: ', name
        open  (unit=11,file=name,form='unformatted',status='old',access='direct',recl=filesizet,action='write')
      else
        ! write(6,*) 'Creating Time stamp driver file: ', name
        open  (unit=11,file=name,form='unformatted',status='replace',access='direct',recl=filesizet,action='write',IOSTAT=IOS)
        ! if (IOS > 0) then
          ! write(6,*) 'IOS = ',IOS
        ! endif
      end if

      ! write(*,*) 'filesizet', filesizet
      ! write(ifoutput)  ( storetdriver (n),  n=1,nstore)
      ! write(6,'(A,F9.2)') 'Writing time stamp to file: ', timee-tdriverstart

      write(11,rec=nstepreaddriver)  ( timee-tdriverstart)
      close (unit=11)
      write(*,*) 'Driver time:' , timee-tdriverstart
    end if
    
    name = 'udriver_   .'
    ! write (name(13:16)  ,'(i4.4)') nfile
    name(9:11)= cmyid
    name(13:15)= cexpnr
    ! name(15:18)= '.txt'
    inquire(file=name,exist=lexist)
    if (lexist) then
      ! write(6,*) 'Writing Inlet u-velocity to file: ', name
      open  (unit=11,file=name,form='unformatted',status='old',access='direct',recl=filesizev,action='write')
    else
      ! write(6,*) 'Creating Inlet u-velocity inlet file: ', name
      open  (unit=11,file=name,form='unformatted',status='replace',access='direct',recl=filesizev,action='write')
    end if

    ! write(ifoutput)  (((storeu0driver (j,k,n),j=jb,je),k=kb,ke),  n=1,nstore) ! Nested implied do-loop
    !if(myid==0) then
      !write(6,'(A,e20.12)') 'Writing u0 to file. u0(irecydriver,je,ke)', u0(irecydriver,je,ke)
      !write(6,'(A,e20.12)') 'u0(irecydriver,jb,kb)', u0(irecydriver,jb,kb)
      !write(6,'(A,e20.12)') 'Writing thl0 to file. thl0(irecydriver-1,je,ke)', thl0(irecydriver-1,je,ke)
      !write(6,'(A,e20.12)') 'thl0(irecydriver-1,jb,kb)', thl0(irecydriver-1,jb,kb)
      ! write(6,*) 'irecydriver, je, ke, ib, jb, kb', irecydriver, je, ke, ib, jb, kb 
    !end if
    write(11,rec=nstepreaddriver)  (u0(irecydriver,:,:))      
    close (unit=11)
  
    name = 'vdriver_   .'
    ! write (name(13:16)  ,'(i4.4)') nfile
    name(9:11)= cmyid
    name(13:15)= cexpnr
    ! name(15:18)= '.txt'
    inquire(file=name,exist=lexist)
    if (lexist) then
      ! write(6,*) 'Writing Inlet v-velocity to file: ', name
      open  (unit=11,file=name,form='unformatted',status='old',action='write',access='direct',recl=filesizev)
    else
      ! write(6,*) 'Creating v-velocity inlet file: ', name
      open  (unit=11,file=name,form='unformatted',status='replace',access='direct',recl=filesizev,action='write')
    end if
    ! write(ifoutput)  (((storev0driver (j,k,n),j=jb,je),k=kb,ke),  n=1,nstore)
    ! '(F8.4)'
    write(11,rec=nstepreaddriver)  (v0(irecydriver,:,:)) !tg3315 removed irecydriver-1
    close (unit=11)
 
    name = 'wdriver_   .'
    ! write (name(13:16)  ,'(i4.4)') nfile
    name(9:11)= cmyid
    name(13:15)= cexpnr
    ! name(15:18)= '.txt'
    inquire(file=name,exist=lexist)
    if (lexist) then
      ! write(6,*) 'Writing Inlet w-velocity to file: ', name
      open  (unit=11,file=name,form='unformatted',status='old',action='write',access='direct',recl=filesizev)
    else
      ! write(6,*) 'Creating w-velocity inlet file: ', name
      open  (unit=11,file=name,form='unformatted',status='replace',access='direct',recl=filesizev,action='write')
    end if
    ! write(ifoutput)  (((storew0driver (j,k,n),j=jb,je),k=kb,ke+1),n=1,nstore)
    write(11,rec=nstepreaddriver)  (w0(irecydriver,:,:)) !tg3315 removed irecydriver-1
    close (unit=11)

    ! name = 'edriver_   .'
    ! write (name(13:16)  ,'(i4.4)') nfile
    ! name(9:11)= cmyid
    ! name(13:15)= cexpnr
    ! name(15:18)= '.txt'
    ! inquire(file=name,exist=lexist)
    ! if (lexist) then
      ! write(6,*) 'Writing Inlet w-velocity to file: ', name
      ! open(unit=11,file=name,form='unformatted',status='old',action='write',access='direct',recl=filesizev)
    ! else
      ! write(6,*) 'Creating w-velocity inlet file: ', name
      ! open(unit=11,file=name,form='unformatted',status='replace',access='direct',recl=filesizev,action='write')
    ! end if
    ! write(ifoutput)  (((storew0driver (j,k,n),j=jb,je),k=kb,ke+1),n=1,nstore)
    ! write(11,rec=nstepreaddriver)  (e120(irecydriver,:,:)) !tg3315 removed irecydriver-1
    ! close (unit=11)

    if (ltempeq ) then
      name = 'hdriver_   .'
      ! write (name(13:16)  ,'(i4.4)') nfile
      name(9:11)= cmyid
      name(13:15)= cexpnr
      ! name(15:18)= '.txt'
      inquire(file=name,exist=lexist)
      if (lexist) then
        ! write(6,*) 'Writing Inlet temperature to file: ', name
        open  (unit=11,file=name,form='unformatted',status='old',action='write',access='direct',recl=filesizev)
      else
        ! write(6,*) 'Creating temperature inlet file: ', name
        ! write(6,*) 'Creating w-velocity inlet file: ', name
        open  (unit=11,file=name,form='unformatted',status='replace',access='direct',recl=filesizev,action='write')
      end if
      ! write(ifoutput)  (((storew0driver (j,k,n),j=jb,je),k=kb,ke+1),n=1,nstore)
      write(11,rec=nstepreaddriver)  (thl0(irecydriver,:,:)) !tg3315 removed irecydriver-1
      close (unit=11)
    end if

    if (lmoist ) then
      name = 'qdriver_   .'
      ! write (name(13:16)  ,'(i4.4)') nfile
      name(9:11)= cmyid
      name(13:15)= cexpnr
      ! name(15:18)= '.txt'
      inquire(file=name,exist=lexist)
      if (lexist) then
        ! write(6,*) 'Writing Inlet temperature to file: ', name
        open  (unit=11,file=name,form='unformatted',status='old',action='write',access='direct',recl=filesizev)
      else
        ! write(6,*) 'Creating temperature inlet file: ', name
        ! write(6,*) 'Creating w-velocity inlet file: ', name
        open  (unit=11,file=name,form='unformatted',status='replace',access='direct',recl=filesizev,action='write')
      end if
      ! write(ifoutput)  (((storew0driver (j,k,n),j=jb,je),k=kb,ke+1),n=1,nstore)
      write(11,rec=nstepreaddriver)  (qt0(irecydriver,:,:)) !tg3315 removed irecydriver-1
      close (unit=11)
    end if

    if (nsv>0 ) then
      name = 'sdriver_   .'
      ! write (name(13:16)  ,'(i4.4)') nfile
      name(9:11)= cmyid
      name(13:15)= cexpnr
      ! name(15:18)= '.txt'
      inquire(file=name,exist=lexist)
      if (lexist) then
        ! write(6,*) 'Writing Inlet temperature to file: ', name
        open  (unit=11,file=name,form='unformatted',status='old',action='write',access='direct',recl=filesizes)
      else
        ! write(6,*) 'Creating temperature inlet file: ', name
        ! write(6,*) 'Creating w-velocity inlet file: ', name
        open  (unit=11,file=name,form='unformatted',status='replace',access='direct',recl=filesizes,action='write')
      end if
      ! write(ifoutput)  (((storew0driver (j,k,n),j=jb,je),k=kb,ke+1),n=1,nstore)
      write(11,rec=nstepreaddriver)  (sv0(irecydriver,:,:,:)) !tg3315 removed irecydriver-1
      close (unit=11)
    end if

  end subroutine writedriverfile

  subroutine readdriverfile
    use modfields, only : u0,sv0
    use modglobal, only : ib,jb,je,jmax,kb,ke,kh,jhc,khc,cexpnr,ifinput,driverstore,ltempeq,lmoist,zh,jgb,jge,jh,driverjobnr,nsv,timee,tdriverstart
    use modmpi,    only : cmyid,myid,nprocs,slabsum,excjs
    use modinletdata, only : storetdriver,storeu0driver,storev0driver,storew0driver,storethl0driver,storeqt0driver,storesv0driver,nfile
    implicit none
    integer :: filen,filee
    integer :: fileid, IOS, filesize, filesizes
    integer :: j,k,m,n,js,jf,jfdum,jsdum
    character(24) :: name

    if (myid==0) then
      write(*,*) '========================================================================'
      write(*,*) '*** Reading precursor driver simulation ***'
    end if

    name = 'tdriver_   .'
    ! write (name(13:16)  ,'(i4.4)') nfile
    name(9:11)= '000'
    ! write (name(18:20)  ,'(i3.3)') filen
    write (name(13:15)   ,'(i3.3)') driverjobnr

    inquire(file=name,size=filesize)

    if(myid==0) then
      write(6,*) 'Reading time stamps: ', name
      write(6,*) 'driverstore: ', driverstore
      write(6,*) 'File size of time in bytes (/8) = ', filesize
    endif
    ! driverstore = driverstore/4.
    ! write(6,*) 'driverstore: ', driverstore
    inquire(iolength=filesize)(timee-tdriverstart)
    open(unit=11,file=name,form='unformatted',status='old',action='read',access='direct',recl=filesize,IOSTAT=IOS)
    if(myid==0) then
      if (IOS > 0) then
        write(6,*) 'IOS = ',IOS
      endif
    endif
    do n =  1, driverstore
      read(11, rec=n, IOSTAT=IOS) storetdriver(n)
      if(myid==0) then
        if(IOS > 0) then
          write(6,*) 'IOS = ',IOS
        elseif (IOS<0) then
          write(6,*) 'n =', n
        end if
        write(6,'(A,e20.12)') ' Reading t:', storetdriver(n)
      end if
    end do
    storetdriver = storetdriver + timee !tg3315 added in case using a warmstart...
    close (unit=11)
    ! write(*,*) 'storetdriver', storetdriver
    ! end if
    name = 'udriver_   .'
    ! write (name(13:16)  ,'(i4.4)') nfile
    name(9:11)= cmyid
    ! write (name(18:20)  ,'(i3.3)') filen
    write (name(13:15)   ,'(i3.3)') driverjobnr
    write(6,*) 'Reading Driver u-velocity: ', name
    ! inquire(file=name,recl=filesize)
    inquire(iolength=filesize)u0(ib,:,:)
    write(6,*) 'record length ',filesize        
    open(unit=11,file=name,form='unformatted',status='old',action='read',access='direct',recl=filesize)
    do n = 1,driverstore
      read(11,rec=n)  ((storeu0driver (j,k,n),j=jb-jh,je+jh),k=kb-kh,ke+kh)
      ! if(myid==0) then
        ! write(6, '(A,e20.12)') 'Reading u(irecydriver, jb, kb)', storeu0driver(jb,kb,n)
      ! endif
    end do
    ! if(myid==0) then
      ! do k=ke,kb,-1
      !   write(6, '(A,e20.12)') 'Reading u(ib,1,:)', storeu0driver(jb,k,1)
      ! end do
    ! end if
    close (unit=11)
      
    name = 'vdriver_   .'
    ! write (name(13:16)  ,'(i4.4)') nfile
    name(9:11)= cmyid
    ! write (name(18:20)  ,'(i3.3)') filen
    write (name(13:15)   ,'(i3.3)') driverjobnr
    write(6,*) 'Reading Driver v-velocity: ', name
    ! inquire(file=name,recl=filesize)
    ! inquire(iolength=filesize)u0(ib,:,:)
    open(unit=11,file=name,form='unformatted',status='old',action='read',access='direct',recl=filesize)
    do n = 1,driverstore
      read(11,rec=n)  ((storev0driver (j,k,n),j=jb-jh,je+jh),k=kb-kh,ke+kh)
    end do
    close (unit=11)

    name = 'wdriver_   .'
    ! write (name(13:16)  ,'(i4.4)') nfile
    name(9:11)= cmyid
    ! write (name(18:20)  ,'(i3.3)') filen
    write (name(13:15)   ,'(i3.3)') driverjobnr
    write(6,*) 'Reading Driver w-velocity: ', name
    ! inquire(file=name,recl=filesize)
    ! inquire(iolength=filesize)u0(ib,:,:)
    open(unit=11,file=name,form='unformatted',status='old',action='read',access='direct',recl=filesize)
    do n = 1,driverstore
      read(11,rec=n)  ((storew0driver (j,k,n),j=jb-jh,je+jh),k=kb-kh,ke+kh)
    end do
    close (unit=11)

    !name = 'edriver_   .'
    ! write (name(13:16)  ,'(i4.4)') nfile
    ! name(9:11)= cmyid
    ! write (name(18:20)  ,'(i3.3)') filen
    ! write (name(13:15)   ,'(i3.3)') driverjobnr
    ! write(6,*) 'Reading Driver turbulent kinetic energy: ', name
    ! inquire(file=name,recl=filesize)
    ! inquire(iolength=filesize)u0(ib,:,:)
    ! open(unit=11,file=name,form='unformatted',status='old',action='read',access='direct',recl=filesize)
    ! do n = 1,driverstore
    ! read(11,rec=n)  ((storee120driver (j,k,n),j=jb-jh,je+jh),k=kb-kh,ke+kh)
    ! enddo
    ! close (unit=11)

    if (ltempeq ) then
      name = 'hdriver_   .'
      ! write (name(13:16)  ,'(i4.4)') nfile
      name(9:11)= cmyid
      ! write (name(18:20)  ,'(i3.3)') filen
      write (name(13:15)   ,'(i3.3)') driverjobnr
      write(6,*) 'Reading Driver temperature: ', name
      ! inquire(file=name,recl=filesize)
      open(unit=11,file=name,form='unformatted',status='old',action='read',access='direct',recl=filesize)
      do n = 1,driverstore
        read(11,rec=n)  ((storethl0driver (j,k,n),j=jb-jh,je+jh),k=kb-kh,ke+kh)
      end do
      !if(myid==0) then
      !  do k=ke,kb,-1
      !    write(6, '(A,e20.12)') 'Reading thl0(ib,1,:)', storethl0driver(jb,k,1)
      !  end do
      !end if

      close (unit=11)
    end if

    if (lmoist ) then
      name = 'qdriver_   .'
      ! write (name(13:16)  ,'(i4.4)') nfile
      name(9:11)= cmyid
      ! write (name(18:20)  ,'(i3.3)') filen
      write (name(13:15)   ,'(i3.3)') driverjobnr
      write(6,*) 'Reading Driver moisture: ', name
      ! inquire(file=name,recl=filesize)
      open(unit=11,file=name,form='unformatted',status='old',action='read',access='direct',recl=filesize)
      do n = 1,driverstore
        read(11,rec=n)  ((storeqt0driver (j,k,n),j=jb-jh,je+jh),k=kb-kh,ke+kh)
      end do
      close (unit=11)
    end if

    if (nsv>0 ) then
      name = 'sdriver_   .'
      ! write (name(13:16)  ,'(i4.4)') nfile
      name(9:11)= cmyid
      ! write (name(18:20)  ,'(i3.3)') filen
      write (name(13:15)   ,'(i3.3)') driverjobnr
      write(6,*) 'Reading Driver scalar: ', name
      ! inquire(file=name,recl=filesize)
      inquire(iolength=filesizes)sv0(ib,:,:,:)
      open(unit=12,file=name,form='unformatted',status='old',action='read',access='direct',recl=filesizes)
      do n = 1,driverstore
        read(12,rec=n)  (((storesv0driver (j,k,m,n),j=jb-jhc,je+jhc),k=kb-khc,ke+khc),m=1,nsv)
      end do
      close (unit=12)
    end if

  end subroutine readdriverfile

  subroutine exitdriver
    use modglobal,      only : idriver,lstoreplane,ltempeq,lmoist,nsv

    if (idriver==1) then
      if (lstoreplane ) then 
        deallocate(storetdriver,storeu0driver,storev0driver,storew0driver)!,storee120driver)
        if (ltempeq ) then
          deallocate(storethl0driver)
        end if
        if (lmoist ) then
          deallocate(storeqt0driver)
        end if
        if (nsv>0 ) then
          deallocate(storesv0driver)
        end if
      end if 
    else if (idriver == 2) then
      deallocate(storetdriver, storeu0driver,storev0driver,storew0driver,u0driver,v0driver,w0driver) !,e120driver,storee120driver)
      if (ltempeq ) then
        deallocate(storethl0driver,thl0driver)
      end if
      if (lmoist ) then
        deallocate(storeqt0driver,qt0driver)
      end if
      if (nsv>0 ) then
        deallocate(storesv0driver,sv0driver)
      end if
    end if

  end subroutine exitdriver

end module
