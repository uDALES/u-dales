!> \file modsave.f90
!! Writes restart and data files.
!> 
!! modsave.f90 writes the restart and data files
!!  \author Jasper Tomas, June 4th 2015
!!  \todo documentation
!!  \par Revision list
!

module modsave


implicit none
! private
! public :: writerestartfiles, writedatafiles
save

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writerestartfiles
    use modsurfdata,only: ustar,thlflux,qtflux,svflux,dudz,dvdz,dthldz,dqtdz,ps,thls,qts,thvs,oblav

    use modfields, only : u0,v0,w0,thl0,qt0,ql0,ql0h,e120,dthvdz,presf,presh,sv0,mindist,wall,&
                          uav,vav,wav,uuav,vvav,wwav,uvav,uwav,vwav,thlav,thl2av,qtav,qlav,ql2av,qt2av,svav,sv2av,momthick,&
                          friction,displthick,pres0,viscratioav,thluav,thlvav,thlwav,qtuav,qtvav,qtwav,qluav,qlvav,qlwav,svuav,svvav,svwav,&
                          upupav,vpvpav,wpwpav,thlpthlpav,qlpqlpav,qtpqtpav,svpsvpav,upvpav,upwpav,vpwpav,thlpupav,thlpvpav,&
                          thlpwpav,qlpupav,qlpvpav,qlpwpav,qtpupav,qtpvpav,qtpwpav,svpupav,svpvpav,svpwpav,presav,&
                          uusgsav,vvsgsav,wwsgsav,uwsgsav,thlusgsav,thlwsgsav,qlusgsav,qlwsgsav,qtusgsav,qtwsgsav,svusgsav,svwsgsav,tkesgsav,&
                          strain2av,disssgsav,t_vav,tvmx,tvmy,tvmz,tsgsmx1,tsgsmx2,tsgsmy1,tsgsmy2,tsgsmz1,&
                          tsgsmz2,t_sgsav,nusgsav,tpm,t_pav,ttmx,ttmy,ttmz,t_tav,p_bav,d_sgsav,p_tav,tkeadv
    use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh,trestart,tnextrestart,dt_lim,timee,btime,xh,&
                          cexpnr,ntimee,rk3step,ifoutput,nsv,timeleft,dt,ntrun,totavtime,&
                          iinletgen,timee,runavtime,inletav,totinletav,linletRA,ltempeq,lmoist,jgb,jge,&
                          dzf,dzfi,dzhi,dxf,dxfi,dyi,dxhi,nstore,numol,dy2i,grav,libm,jmax,nblocks
    use modmpi,    only : cmyid,myid,slabsum,excjs
    use modsubgriddata, only : ekm
    use modibmdata,   only  : ibmxforcevol
    use initfac , only : block
    use modinletdata, only   : Urec,Wrec,Uinl,Utav,QLinl,QTinl,QLrec,QTrec,QTtav,QLtav,Ttav,upupavinl,vpvpavinl,wpwpavinl,upwpavinl,&
                               thlpthlpavinl,thlpupavinl,thlpwpavinl,qlpqlpavinl,qlpupavinl,qlpwpavinl,qtpqtpavinl,qtpupavinl,qtpwpavinl,Tinl,Trec,nstepread

    implicit none
    logical :: lexitnow = .false.
    integer imin,ihour
    integer i,j,k,n,im,ip,jm,jp,jpp,km,kp,kpp,il,iu,jl,ju,kl,ku
    character(21) name,name2,name3,name4,linkname

    if (timee == 0) return
!    if (rk3step /=3) return
    if ((iinletgen==2) .and. (nstepread==nstore)) then                ! This overrules the need for rk3step to be 3 in case of reading inletfiles
      write(6,*) 'Writing restartfiles after reading in new inletfiles'
    else
      if (rk3step /=3) return   ! Normal check
    end if

    name = 'exit_now.'//cexpnr
    inquire(file=trim(name), EXIST=lexitnow)

    if (((timee>=tnextrestart)) .or. ((lexitnow) .or. (nstepread == nstore+1))) then
      tnextrestart = tnextrestart+trestart
    
      name = 'initd        _   .'
      write (name(6:13)  ,'(i8.8)') ntrun
      name(15:17)= cmyid
      name(19:21)= cexpnr
      open  (ifoutput,file=name,form='unformatted',status='replace')

      write(ifoutput)  (((mindist(i,j,k),i=ib,ie  ),j=jb,je      ),k=kb,ke   )
      write(ifoutput)  ((((wall(i,j,k,n),i=ib,ie  ),j=jb,je      ),k=kb,ke   ),n=1,5)
      write(ifoutput)  (((u0 (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((v0 (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((w0    (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((pres0 (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((thl0  (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((e120  (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((ekm   (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((qt0   (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((ql0   (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  (((ql0h  (i,j,k),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh)
      write(ifoutput)  timee,  dt
      
      if (myid==0) then
        write(*,*) '-------------------------'
        write(*,*) 'Saving initd restart file'
        write(*,*) 'ntrun ::: ', ntrun
        write(*,*) 'timee ::: ', timee
        write(*,*) '-------------------------'
      endif
      
      close (ifoutput)

      if (nsv>0) then
        name  = 'inits        _   .'
        write (name(6:13) ,'(i8.8)') ntrun
        name(15:17) = cmyid
        name(19:21) = cexpnr
        open  (ifoutput,file=name,form='unformatted')
        write(ifoutput) ((((sv0(i,j,k,n),i=ib-ih,ie+ih),j=jb-jh,je+jh),k=kb,ke+kh),n=1,nsv)
        write(ifoutput)  timee

        close (ifoutput)
      end if

      if (myid==0) then
        write(*,'(A,F15.7,A,I4)') 'dump at time = ',timee,' unit = ',ifoutput
      end if

    end if

  end subroutine writerestartfiles

end module modsave
