!> \file modparticles.f90
!!  Tracks the flow through Lagrangian particles

!>
!!  Tracks the flow through Lagrangian particles
!>
!!  \see Heus et al (JAS 2008), Verzijlbergh et al (ACP 2009)
!!  \author Gert Jan van Dijk, TU Delft
!!  \author Remco Verzijlbergh,TU Delft
!!  \author Harm Jonker, TU Delft
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  \todo Documentation
!!  \todo NetCDF implementation
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

module modparticles
  use modglobal, only : longint
implicit none
PRIVATE
PUBLIC :: initparticles, particles, exitparticles
SAVE

  logical :: lpartic    = .false.
  integer :: intmeth    =    3          !<     * Choose the integrationscheme you would like to use.
  logical :: lpartsgs   = .true.        !<  switch for particle subgrid velocity (on/off)
  logical :: lstat      = .false.
  real    :: timeav     = 3600.
  real    :: dtav       =   60.
  logical :: ldump      = .false.
  real    :: timedump   = 3600.         !<     * write the solution every dtwrite seconds
  integer :: npartdump  = 10
  character(30) :: startfilepart

  integer :: nsamples
  integer(kind=longint)    :: itimeav,idtav,itimedump,tnext,tnextwrite,tnextdump

  integer,parameter  :: inomove=0                       !< the options for the integrationscheme
  integer,parameter  :: irk3=3

  integer (KIND=selected_int_kind(10)):: idum = -12345  !< seed for gaussian random variable (must be negative integer)


  real :: ysizelocal
  real :: dsigma2dx_sgs=0, dsigma2dy_sgs=0, dsigma2dz_sgs=0, dsigma2dt_sgs=0,sigma2_new = 0
  real :: fce
  real,allocatable :: fsm(:)
  integer :: np

  TYPE :: particle_record
    real :: unique, tstart
    integer :: partstep
    real :: x,x_prev, xstart, ures, usgs, usgs_prev
    real :: y,y_prev, ystart, vres, vsgs, vsgs_prev
    real :: z,z_prev, zstart, wres, wsgs, wsgs_prev
    real ::sigma2_sgs

    TYPE (particle_record), POINTER:: next,prev
  end TYPE

  integer :: nrpartvar,ipunique,ipx,ipy,ipz,ipxstart,ipystart,ipzstart,iptsart,ipartstep,ipures,ipvres,ipwres
  integer :: ipuresprev,ipvresprev,ipwresprev,ipxprev,ipyprev,ipzprev,i,j,k,nsteps,nplisted
  integer :: ipusgs_prev, ipvsgs_prev, ipwsgs_prev, ipusgs, ipvsgs, ipwsgs
  integer :: ipsigma2_sgs

  TYPE (particle_record), POINTER:: head, tail

contains

!*****************************************************************************
  subroutine initparticles
    use modmpi,   only : myid,my_real,mpierr,comm3d,mpi_integer,mpi_logical,nprocs
    use modglobal,only : ifnamopt,fname_options,ifinput,dtmax,cexpnr,&
                         dx,dy,dzf,zh,kmax,k1,iexpnr,runtime,timee,ysize,dt_lim,btime,rtimee,tres

    implicit none

  ! LOCAL
    integer :: n, ierr
    real :: tstart, xstart, ystart, zstart
    type (particle_record), pointer:: particle
    namelist/NAMPARTICLES/ &
          lpartic, &
          intmeth, &
          lpartsgs, &
          lstat,&
          dtav, &
          timeav, &
          ldump,&
          timedump, &
          npartdump


    np = 0
    startfilepart = 'partstartpos.'//cexpnr
  !read namelist and broadcast it to all processors

      if(myid==0)then
        open(ifnamopt,file=fname_options,status='old',iostat=ierr)
        read (ifnamopt,NAMPARTICLES,iostat=ierr)
        if (ierr > 0) then
          print *, 'Problem in namoptions NAMPARTICLES'
          print *, 'iostat error: ', ierr
          stop 'ERROR: Problem in namoptions NAMPARTICLES'
        endif
        write(6 ,NAMPARTICLES)
        close(ifnamopt)
      end if


    call MPI_BCAST(lpartic     ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(timeav      ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(dtav        ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(timedump    ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(npartdump   ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(intmeth     ,1,MPI_INTEGER,0,comm3d,mpierr)
    call MPI_BCAST(lpartsgs    ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(lstat       ,1,MPI_LOGICAL,0,comm3d,mpierr)
    call MPI_BCAST(ldump       ,1,MPI_LOGICAL,0,comm3d,mpierr)

    if (.not.(lpartic)) return
      if(lstat) then
        idtav = dtav/tres
        itimeav = timeav/tres
        itimedump = timedump/tres

        tnext      = idtav   +btime
        tnextwrite = itimeav +btime
        nsamples = itimeav/idtav
        tnextdump  = itimedump+btime
        dt_lim = min(dt_lim,tnext)
      end if

    allocate(fsm(k1))

      !this is the file which contains the startpositions AND the number of particles             !#########

    open(ifinput,file='partstartpos.'//cexpnr,status='old',position = 'rewind',action='read')
    read(ifinput,'(I10.1)') np

    if ( np < 1 ) return

    ysizelocal = ysize/nprocs

    ipunique =  1
    ipx = 2
    ipy = 3
    ipz = 4
    ipxstart = 5
    ipystart = 6
    ipzstart =  7
    iptsart = 8
    ipures = 9
    ipvres = 10
    ipwres = 11
    ipusgs = 12
    ipvsgs = 13
    ipwsgs = 14
    ipusgs_prev = 15
    ipvsgs_prev = 16
    ipwsgs_prev = 17
    ipartstep = 18
    nrpartvar  = 18
    if (intmeth == irk3 .or.lpartsgs) then
      ipuresprev = nrpartvar+1
      ipvresprev = nrpartvar+2
      ipwresprev = nrpartvar+3
      ipxprev = nrpartvar+4
      ipyprev = nrpartvar+5
      ipzprev = nrpartvar+6
      nrpartvar=nrpartvar+6
    end if

    if (lpartsgs) then
      ipsigma2_sgs = nrpartvar + 1
      nrpartvar = nrpartvar + 1
    end if

    call particle_initlist()
    nplisted = 0


  ! initialize particle list & positions
    do n = 1, np
      read(ifinput,*) tstart, xstart, ystart, zstart

      if (floor(ystart/ysizelocal) == nprocs) ystart = 0.
      if (floor(ystart/ysizelocal) == myid) then
        call particle_add(particle)

        particle%unique = n + myid/1000.0
        particle%x = xstart/dx +2
        particle%y = (ystart-ysizelocal*myid)/dy+2
        do k=kmax,1,-1
          if (zh(k)<zstart) exit
        end do
        particle%z = k + (zstart-zh(k))/dzf(k)
        particle%xstart = xstart
        particle%ystart = ystart
        particle%zstart = zstart
        particle%tstart = tstart
        particle%ures = 0.
        particle%vres = 0.
        particle%wres = 0.
        particle%usgs = 0.
        particle%vsgs = 0.
        particle%wsgs = 0.
        particle%usgs_prev = 0.
        particle%vsgs_prev = 0.
        particle%wsgs_prev = 0.
        particle%x_prev = particle%x
        particle%y_prev = particle%y
        particle%z_prev = particle%z
        particle%partstep = 0
        particle%sigma2_sgs = epsilon(particle%sigma2_sgs)
      end if
    end do

    close(ifinput)

    write(6,*) 'timee :',rtimee,': proc ',myid,': #particles: ',nplisted

  end subroutine initparticles

!*****************************************************************************

  subroutine exitparticles

    use modmpi, only : myid
    use modglobal,only:rtimee
    implicit none


    if (.not.(lpartic)) return

    write(6,*) rtimee,": proc ",myid,": removing all particles from list..."
    call particle_quitlist()
    deallocate(fsm)
  end subroutine exitparticles

!*****************************************************************************

  subroutine particles
    use modglobal, only :dx,dy,dzf,k1,rk3step,rtimee,j2,i2
    implicit none
    type (particle_record), pointer:: particle

    if (.not.(lpartic)) return
    if ( np < 1 ) return
    if (rtimee == 0) return

    if (lpartsgs) call sgsinit

    particle => head
    do while( associated(particle) )
      if (  rtimee - particle%tstart >= 0 ) then
        particle%partstep = particle%partstep + 1

!interpolation of the velocity field
        if (  rtimee - particle%tstart >= 0 ) then
          particle%ures = velocity_ures(particle%x,particle%y,particle%z) / dx
          particle%vres = velocity_vres(particle%x,particle%y,particle%z) / dy
          particle%wres = velocity_wres(particle%x,particle%y,particle%z) / dzf(floor(particle%z))

          if (lpartsgs) then
            if (rk3step==1) then
              particle%usgs_prev = particle%usgs
              particle%vsgs_prev = particle%vsgs
              particle%wsgs_prev = particle%wsgs
            end if
            call sgshelpvar(particle)
            particle%usgs = velocity_usgs(particle) / dx
            particle%vsgs = velocity_vsgs(particle) / dy
            particle%wsgs = velocity_wsgs(particle) / dzf(floor(particle%z))
          end if
        end if
      end if


    particle => particle%next
    end do
!****Statistics******
    if (rk3step==3) then
      call statistics
      call writeparticles
    end if

!Time integration
    particle => head
    do while( associated(particle) )
      if (  rtimee - particle%tstart >= 0 ) then
        select case(intmeth)
        case(inomove)
                ! no movement
        case(irk3)
                call rk3(particle)
        case default
                stop 'PARTICLES ERROR: incorrect integration scheme'
        end select
      end if
    particle => particle%next
    end do

!Exchange particle to other processors
    call partcommunicate


  end subroutine particles

!######################################################################################################

!*****************************************************************************

  subroutine rk3(particle)
    use modglobal, only : rk3step,rdt,j1,dzf
    implicit none
    real :: rk3coef

    TYPE (particle_record), POINTER:: particle
    rk3coef = rdt / (4. - dble(rk3step))
    particle%x =  particle%x_prev + (particle%ures +particle%usgs) * rk3coef
    particle%y =  particle%y_prev + (particle%vres +particle%vsgs) * rk3coef
    particle%z =  particle%z_prev + (particle%wres +particle%wsgs) * rk3coef
    call checkbound(particle)
    if (floor(particle%z)/=floor(particle%z_prev)) then
      particle%z = floor(particle%z) + (particle%z -floor(particle%z))*dzf(floor(particle%z_prev))/dzf(floor(particle%z))
    end if

    if (rk3step==3) then
      particle%x_prev = particle%x
      particle%y_prev = particle%y
      particle%z_prev = particle%z
    end if

  end subroutine rk3

!*****************************************************************************
  function velocity_ures(x, y, z)

    use modglobal, only : cu,zh,zf,dzh,dzf
    use modfields, only : u0
    implicit none
    real, intent(in) ::  x, y, z
    integer :: xbottom, ybottom, zbottom
    real :: velocity_ures, deltax, deltay, deltaz

    xbottom = floor(x)
    ybottom = floor(y - 0.5)
    zbottom = floor(z - 0.5)
    deltax = x - xbottom
    deltay = y - 0.5 - ybottom

    if (zbottom==0)  then
      deltaz = 2.*z - 2.
      velocity_ures = (1-deltaz) * (1-deltay) * (1-deltax) * (-cu) + &
      &               (1-deltaz) * (1-deltay) * (  deltax) * (-cu) + &
      &               (1-deltaz) * (  deltay) * (1-deltax) * (-cu) + &
      &               (1-deltaz) * (  deltay) * (  deltax) * (-cu) + &
      &               (  deltaz) * (1-deltay) * (1-deltax) * u0(xbottom    , ybottom    , 1)  + &
      &               (  deltaz) * (1-deltay) * (  deltax) * u0(xbottom + 1, ybottom    , 1)  + &
      &               (  deltaz) * (  deltay) * (1-deltax) * u0(xbottom    , ybottom + 1, 1)  + &
      &               (  deltaz) * (  deltay) * (  deltax) * u0(xbottom + 1, ybottom + 1, 1)
    else
      deltaz = (zh(floor(z)) + dzf(floor(z))*(z-floor(z)) - zf(zbottom))/dzh(zbottom)
      velocity_ures =  (1-deltaz) * (1-deltay) * (1-deltax) * u0(xbottom    , ybottom    , zbottom    ) + &
      &                (1-deltaz) * (1-deltay) * (  deltax) * u0(xbottom + 1, ybottom    , zbottom    ) + &
      &                (1-deltaz) * (  deltay) * (1-deltax) * u0(xbottom    , ybottom + 1, zbottom    ) + &
      &                (1-deltaz) * (  deltay) * (  deltax) * u0(xbottom + 1, ybottom + 1, zbottom    ) + &
      &                (  deltaz) * (1-deltay) * (1-deltax) * u0(xbottom    , ybottom    , zbottom + 1) + &
      &                (  deltaz) * (1-deltay) * (  deltax) * u0(xbottom + 1, ybottom    , zbottom + 1) + &
      &                (  deltaz) * (  deltay) * (1-deltax) * u0(xbottom    , ybottom + 1, zbottom + 1) + &
      &                (  deltaz) * (  deltay) * (  deltax) * u0(xbottom + 1, ybottom + 1, zbottom + 1)
    end if

  end function velocity_ures

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function velocity_vres(x, y, z)

    use modglobal, only : cv,zh,zf,dzf,dzh
    use modfields, only : vm,v0

    real, intent(in) :: x, y, z

    !local::
    integer :: xbottom, ybottom, zbottom
    real ::  velocity_vres, deltax, deltay, deltaz

    xbottom = floor(x - 0.5)
    ybottom = floor(y)
    zbottom = floor(z - 0.5)

    deltax = x - 0.5 - xbottom
    deltay = y -  ybottom


    if (zbottom==0)  then
      deltaz = 2.*z - 2.
      velocity_vres = (1-deltaz) * (1-deltay) * (1-deltax) * (-cv) + &
      &               (1-deltaz) * (1-deltay) * (  deltax) * (-cv) + &
      &               (1-deltaz) * (  deltay) * (1-deltax) * (-cv) + &
      &               (1-deltaz) * (  deltay) * (  deltax) * (-cv) + &
      &               (  deltaz) * (1-deltay) * (1-deltax) * v0(xbottom    , ybottom    , 1)  + &
      &               (  deltaz) * (1-deltay) * (  deltax) * v0(xbottom + 1, ybottom    , 1)  + &
      &               (  deltaz) * (  deltay) * (1-deltax) * v0(xbottom    , ybottom + 1, 1)  + &
      &               (  deltaz) * (  deltay) * (  deltax) * v0(xbottom + 1, ybottom + 1, 1)
    else
      deltaz = (zh(floor(z)) + dzf(floor(z))*(z-floor(z)) - zf(zbottom))/dzh(zbottom)
      velocity_vres =  (1-deltaz) * (1-deltay) * (1-deltax) * v0(xbottom    , ybottom    , zbottom    ) + &
      &                (1-deltaz) * (1-deltay) * (  deltax) * v0(xbottom + 1, ybottom    , zbottom    ) + &
      &                (1-deltaz) * (  deltay) * (1-deltax) * v0(xbottom    , ybottom + 1, zbottom    ) + &
      &                (1-deltaz) * (  deltay) * (  deltax) * v0(xbottom + 1, ybottom + 1, zbottom    ) + &
      &                (  deltaz) * (1-deltay) * (1-deltax) * v0(xbottom    , ybottom    , zbottom + 1) + &
      &                (  deltaz) * (1-deltay) * (  deltax) * v0(xbottom + 1, ybottom    , zbottom + 1) + &
      &                (  deltaz) * (  deltay) * (1-deltax) * v0(xbottom    , ybottom + 1, zbottom + 1) + &
      &                (  deltaz) * (  deltay) * (  deltax) * v0(xbottom + 1, ybottom + 1, zbottom + 1)
    end if

  end function velocity_vres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function velocity_wres(x, y, z)

    use modfields, only : wm,w0

    real, intent(in) :: x, y, z

    !local::
    integer :: xbottom, ybottom, zbottom
    real :: velocity_wres, deltax, deltay, deltaz

    xbottom = floor(x - 0.5)
    ybottom = floor(y - 0.5)
    zbottom = floor(z )

    deltax = x - 0.5 - xbottom
    deltay = y - 0.5 - ybottom
    deltaz = z - zbottom

    !interpolation of 'w'.

    velocity_wres =  (1-deltaz) * (1-deltay) * (1-deltax) * w0(xbottom    , ybottom    , zbottom    ) + &
    &                (1-deltaz) * (1-deltay) * (  deltax) * w0(xbottom + 1, ybottom    , zbottom    ) + &
    &                (1-deltaz) * (  deltay) * (1-deltax) * w0(xbottom    , ybottom + 1, zbottom    ) + &
    &                (1-deltaz) * (  deltay) * (  deltax) * w0(xbottom + 1, ybottom + 1, zbottom    ) + &
    &                (  deltaz) * (1-deltay) * (1-deltax) * w0(xbottom    , ybottom    , zbottom + 1) + &
    &                (  deltaz) * (1-deltay) * (  deltax) * w0(xbottom + 1, ybottom    , zbottom + 1) + &
    &                (  deltaz) * (  deltay) * (1-deltax) * w0(xbottom    , ybottom + 1, zbottom + 1) + &
    &                (  deltaz) * (  deltay) * (  deltax) * w0(xbottom + 1, ybottom + 1, zbottom + 1)


  end function velocity_wres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function velocity_usgs(particle)

    use modglobal, only : dx,rdt
    implicit none

    !local::
    real :: velocity_usgs
    TYPE (particle_record), POINTER:: particle

      velocity_usgs = -0.5*fce*particle%usgs_prev*dx*rdt + &
                      0.5*(dsigma2dt_sgs*particle%usgs_prev*dx+dsigma2dx_sgs)*rdt + &
                      sqrt(fce*sigma2_new)*gauss1(idum) + &
                      particle%usgs_prev*dx


  end function velocity_usgs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function velocity_vsgs(particle)

    use modglobal, only : dy,rdt
    implicit none

    !local::
    real :: velocity_vsgs
    TYPE (particle_record), POINTER:: particle

      velocity_vsgs = -0.5*fce*particle%vsgs_prev*dy*rdt + &
                      0.5*(dsigma2dt_sgs*particle%vsgs_prev*dy+dsigma2dy_sgs)*rdt + &
                      sqrt(fce*sigma2_new)*gauss1(idum) + &
                      particle%vsgs_prev*dy


  end function velocity_vsgs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function velocity_wsgs(particle)

    use modglobal, only : dzf,rdt
    implicit none

    !local::
    real :: velocity_wsgs
    TYPE (particle_record), POINTER:: particle

      velocity_wsgs = -0.5*fce*particle%wsgs_prev*dzf(floor(particle%z))*rdt + &
                      0.5*(dsigma2dt_sgs*particle%wsgs_prev*dzf(floor(particle%z))+dsigma2dz_sgs)*rdt + &
                      sqrt(fce*sigma2_new)*gauss1(idum) + &
                      particle%wsgs_prev*dzf(floor(particle%z))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Remco Verzijlbergh: This is the original form of the sub-grid w-velocity equation. Tests have  !!!!!
!revealed that the form as given above (without the dsigma2dz_sgs-term) performs slightly better!!!!!
! with respect to accumulation/depletion of particles at certain regions                        !!!!!
!                        !!!!!
!      -0.5*fce*wsgs_previous*dz*dtpart + &                                    !!!!!
!                       0.5*(dsigma2dt_sgs*wsgs_previous*dz+dsigma2dz_sgs)*dtpart + &           !!!!!
!                       sqrt(fce*sigma2_new)*gauss1(idum) + &                                   !!!!!
!                       wsgs_previous*dz                                                        !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end function velocity_wsgs


  subroutine statistics
    use modglobal, only : kmax,k1,zf,dx,dy,dzf,ifoutput,cexpnr,timee,dt_lim,cu,cv
    use modmpi,    only : comm3d,mpierr,mpi_sum,my_real,mpi_integer,myid
    implicit none


    !definitions & initializations
    type (particle_record), pointer:: particle



    integer,allocatable, dimension(:) :: nrparticlesavl , nrparticlesav
     real,allocatable, dimension(:) :: uresavl  ,uresav,&
                                       vresavl  ,vresav,&
                                       wresavl  ,wresav,&
                                       u2resavl ,u2resav,&
                                       v2resavl ,v2resav,&
                                       w2resavl ,w2resav,&
                                       usgsavl  ,usgsav,&
                                       vsgsavl  ,vsgsav,&
                                       wsgsavl  ,wsgsav,&
                                       u2sgsavl ,u2sgsav,&
                                       v2sgsavl ,v2sgsav,&
                                       w2sgsavl ,w2sgsav,&
                                       tkesgsavl,tkesgsav


    if (.not. lstat) return
    if(timee<tnext .and. timee<tnextwrite) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if

     allocate(nrparticlesavl(k1) , nrparticlesav(k1))
     allocate( uresavl(k1)  ,uresav(k1),&
               vresavl(k1)  ,vresav(k1),&
               wresavl(k1)  ,wresav(k1),&
               u2resavl(k1) ,u2resav(k1),&
               v2resavl(k1) ,v2resav(k1),&
               w2resavl(k1) ,w2resav(k1),&
               usgsavl(k1)  ,usgsav(k1),&
               vsgsavl(k1)  ,vsgsav(k1),&
               wsgsavl(k1)  ,wsgsav(k1),&
               u2sgsavl(k1) ,u2sgsav(k1),&
               v2sgsavl(k1) ,v2sgsav(k1),&
               w2sgsavl(k1) ,w2sgsav(k1),&
               tkesgsavl(k1),tkesgsav(k1))


    nrparticlesavl=0
    uresavl = 0
    vresavl = 0
    wresavl = 0
    u2resavl = 0
    v2resavl = 0
    w2resavl = 0
    usgsavl = 0
    vsgsavl = 0
    wsgsavl = 0
    u2sgsavl = 0
    v2sgsavl = 0
    w2sgsavl = 0
    tkesgsavl = 0

    if (timee>=tnext) then
      tnext = tnext+idtav
      dt_lim = minval((/dt_lim,tnext-timee/))
      particle =>head
      do while (associated(particle))

        nrparticlesavl(floor(particle%z)) =  nrparticlesavl(floor(particle%z)) + 1
        uresavl(floor(particle%z)) =uresavl(floor(particle%z)) + particle%ures
        vresavl(floor(particle%z)) =vresavl(floor(particle%z)) + particle%vres
        wresavl(floor(particle%z)) =wresavl(floor(particle%z)) + particle%wres
        u2resavl(floor(particle%z)) =u2resavl(floor(particle%z)) + particle%ures**2
        v2resavl(floor(particle%z)) =v2resavl(floor(particle%z)) + particle%vres**2
        w2resavl(floor(particle%z)) =w2resavl(floor(particle%z)) + particle%wres**2

        usgsavl(floor(particle%z)) =usgsavl(floor(particle%z)) + particle%usgs
        vsgsavl(floor(particle%z)) =vsgsavl(floor(particle%z)) + particle%vsgs
        wsgsavl(floor(particle%z)) =wsgsavl(floor(particle%z)) + particle%wsgs
        u2sgsavl(floor(particle%z)) =u2sgsavl(floor(particle%z)) + particle%usgs**2
        v2sgsavl(floor(particle%z)) =v2sgsavl(floor(particle%z)) + particle%vsgs**2
        w2sgsavl(floor(particle%z)) =w2sgsavl(floor(particle%z)) + particle%wsgs**2


        if (lpartsgs) tkesgsavl(floor(particle%z)) =tkesgsavl(floor(particle%z)) +1.5*particle%sigma2_sgs


        particle => particle%next
      end do
      !MPI communication

      call MPI_ALLREDUCE(nrparticlesavl,nrparticlesav, k1,    MPI_INTEGER,  MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(uresavl,uresav, k1,    MY_REAL,  MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(vresavl,vresav, k1,    MY_REAL,  MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(wresavl,wresav, k1,    MY_REAL,  MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(u2resavl,u2resav, k1,    MY_REAL,  MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(v2resavl,v2resav, k1,    MY_REAL,  MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(w2resavl,w2resav, k1,    MY_REAL,  MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(usgsavl,usgsav, k1,    MY_REAL,  MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(vsgsavl,vsgsav, k1,    MY_REAL,  MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(wsgsavl,wsgsav, k1,    MY_REAL,  MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(u2sgsavl,u2sgsav, k1,    MY_REAL,  MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(v2sgsavl,v2sgsav, k1,    MY_REAL,  MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(w2sgsavl,w2sgsav, k1,    MY_REAL,  MPI_SUM, comm3d,mpierr)
      call MPI_ALLREDUCE(tkesgsavl,tkesgsav, k1,    MY_REAL,  MPI_SUM, comm3d,mpierr)

    end if

    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+itimeav

      !normalization
      where (nrparticlesav>0)
        uresav = uresav/nrparticlesav*dx + cu
        vresav = vresav/nrparticlesav*dy + cv
        wresav = wresav/nrparticlesav*dzf
        u2resav = u2resav/nrparticlesav*dx**2-uresav**2
        v2resav = v2resav/nrparticlesav*dy**2-vresav**2
        w2resav = w2resav/nrparticlesav*dzf**2-wresav**2
        usgsav = usgsav/nrparticlesav*dx
        vsgsav = vsgsav/nrparticlesav*dy
        wsgsav = wsgsav/nrparticlesav*dzf
        u2sgsav = u2sgsav/nrparticlesav*dx**2-usgsav**2
        v2sgsav = v2sgsav/nrparticlesav*dy**2-vsgsav**2
        w2sgsav = w2sgsav/nrparticlesav*dzf**2-wsgsav**2
        tkesgsav = tkesgsav/nrparticlesav


        nrparticlesav = nrparticlesav/nsamples
      end where

      !output
        if (myid==0) then
          open(ifoutput,file='partstat.'//cexpnr,position='append',action='write')
          write(ifoutput,'(//A,/A,F5.0,A,F7.0,A)') &
          '#--------------------------------------------------------'      &
          ,'#',(timeav),'--- AVERAGING TIMESTEP --- '      &
          ,timee      &
          ,'   SEC AFTER INITIALIZATION '

          write (ifoutput,'(A/2A/2A)') &
              '#--------------------------------------------------------' &
              ,'#LEV  HGHT   NRPART     URES         VRES         WRES        U2RES        V2RES        W2RES' &
              ,'           USGS         VSGS         WSGS        U2SGS        V2SGS        W2SGS        TKESGS' &
              ,'#      (M)    (-)    (----------- (M/S) ---------------)    (----------- (M/S)^2 ---------------)' &
              ,'    (----------- (M/S) ---------------)    (-------------- (M/S)^2 ------------------)'

          write(ifoutput,'(F8.2,I10,13E13.5)') &
                (zf       (k), nrparticlesav(k),uresav(k),vresav(k),wresav(k),u2resav(k),v2resav(k),w2resav(k),&
                usgsav(k),vsgsav(k),wsgsav(k),u2sgsav(k),v2sgsav(k),w2sgsav(k),tkesgsav(k),k=1,kmax)
          close(ifoutput)
        end if



      !reset averages
        uresavl = 0
        vresavl = 0
        wresavl = 0
        u2resavl = 0
        v2resavl = 0
        w2resavl = 0
        usgsavl = 0
        vsgsavl = 0
        wsgsavl = 0
        u2sgsavl = 0
        v2sgsavl = 0
        w2sgsavl = 0
        tkesgsavl = 0
        nrparticlesavl = 0


        uresav = 0
        vresav = 0
        wresav = 0
        u2resav = 0
        v2resav = 0
        w2resav = 0
        usgsav = 0
        vsgsav = 0
        wsgsav = 0
        u2sgsav = 0
        v2sgsav = 0
        w2sgsav = 0
        tkesgsav = 0
        nrparticlesav = 0


    end if
    deallocate(nrparticlesavl , nrparticlesav)
    deallocate( uresavl  ,uresav,&
              vresavl  ,vresav,&
              wresavl  ,wresav,&
              u2resavl ,u2resav,&
              v2resavl ,v2resav,&
              w2resavl ,w2resav,&
              usgsavl  ,usgsav,&
              vsgsavl  ,vsgsav,&
              wsgsavl  ,wsgsav,&
              u2sgsavl ,u2sgsav,&
              v2sgsavl ,v2sgsav,&
              w2sgsavl ,w2sgsav,&
              tkesgsavl,tkesgsav)

  end subroutine statistics

 !*****************************************************************************
  subroutine writeparticles
    use modmpi,    only : myid,cmyid
    use modglobal, only :  i2,jmax,j2,k1,dx,dy,dzf,dzh,zf,zh,es0,tmelt,rlv,rd,rv,cp,bt,at,cexpnr,ifoutput,timee,rtimee,dt_lim
    use modfields, only : qtm,thlm,presf, exnf
    use modsurfdata,only: thvs,thls,qts
    implicit none

    ! LOCAL
    real, allocatable,dimension (:,:) :: partdata
    integer(KIND=selected_int_kind(11)), allocatable, dimension (:) :: partids
    integer :: n,m,ndata
    type (particle_record), pointer:: particle
    real :: thlpart, qtpart, qlpart, thvpart, exnpart,prespart,tl,es,qs,qsl,b1

    ndata=npartdump

    if (.not. ldump) return
    if(timee<tnextdump) then
      dt_lim = minval((/dt_lim,tnextdump-timee/))
      return
    end if

    if (timee>=tnextdump) then
      tnextdump =tnextdump + itimedump
      allocate (partdata(ndata,nplisted))
      allocate (partids(nplisted))

      n = 0

      particle => head
      do while( associated(particle) )
        n = n + 1

        thlpart = scalintp(particle%x,particle%y,particle%z,thlm,thls)
        qtpart  = scalintp(particle%x,particle%y,particle%z, qtm,qts )

        exnpart = exnf(floor(particle%z))
        prespart= presf(floor(particle%z))

        tl = thlpart*exnpart
        es  = es0*exp(at*(tl-tmelt)/(tl-bt))
        qsl = rd/rv*es/(prespart-(1-rd/rv)*es)
        b1  = rlv/(rv*tl)*(rlv/(cp*tl))
        qs  = qsl*(1.+b1*qtpart)/(1.+b1*qsl)

        qlpart  = dim(qtpart-qs,0.0)
        thvpart = (thlpart+(rlv/cp)*qlpart/exnpart) * (1.+(rv/rd-1)*qtpart-rv/rd*qlpart)

         partids(n) = particle%unique
        if (ndata > 0) then
           partdata(1,n) = (particle%x-2)*dx
        endif
        if (ndata > 1) then
           partdata(2,n) = (jmax*myid+particle%y-2)*dy
        endif
        if (ndata > 2) then
           partdata(3,n) = zh(floor(particle%z)) + dzf(floor(particle%z))*(particle%z-floor(particle%z))
        endif
        if (ndata > 3) then
           partdata(4,n) = (particle%ures+particle%usgs)*dx
        endif
        if (ndata > 4) then
           partdata(5,n) = (particle%vres+particle%vsgs)*dy
        endif
        if (ndata > 5) then
           partdata(6,n) = (particle%wres+particle%wsgs)*dzf(floor(particle%z))
        endif
        if (ndata > 6) then
           partdata(7,n) = thlpart
        endif
        if (ndata > 7) then
           partdata(8,n) = thvpart
        endif
        if (ndata > 8) then
           partdata(9,n) = qtpart
        endif
        if (ndata > 9) then
           partdata(10,n)= qlpart
        endif
        particle => particle%next
      end do

      open(ifoutput,file='particles.'//cmyid//'.'//cexpnr,form='unformatted',position = 'append',action='write')
      write(ifoutput) rtimee,nplisted
      write(ifoutput) (partids(n),(partdata(m,n),m=1,ndata),n=1,nplisted)
      close(ifoutput)
      deallocate (partdata)
      deallocate (partids)

    end if

  end subroutine writeparticles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sgsinit
!calculation of several 3D fields necesarry for sgs model
!fce

  use modglobal, only : i1,j1,kmax,k1,rslabs
  use modfields, only : u0,v0,w0,e120
  use modmpi,    only : comm3d,my_real,mpi_sum,mpierr
  implicit none


  real,allocatable, dimension(:) :: &
        umslabl,vmslabl,um2slabl,vm2slabl,wm2slabl,tkesgsmslabl, &
        umslab ,vmslab ,um2slab ,vm2slab ,wm2slab ,tkesgsmslab,&
        tkeresmslab
  real :: fsmin=1

  allocate(umslabl(k1),vmslabl(k1),um2slabl(k1),vm2slabl(k1),wm2slabl(k1)&
        ,tkesgsmslabl(k1),umslab (k1),vmslab (k1),um2slab (k1),vm2slab (k1) &
        ,wm2slab (k1),tkesgsmslab(k1),tkeresmslab(k1))

  if (fsmin==1) then
    fsm = 1

    deallocate( &
        umslabl,vmslabl,um2slabl,vm2slabl,wm2slabl,tkesgsmslabl, &
        umslab ,vmslab ,um2slab ,vm2slab ,wm2slab ,tkesgsmslab,&
        tkeresmslab)
    return
  end if
  !calc <sigma2_sfs>
  do k=1,k1
    umslabl(k) = sum(u0(2:i1,2:j1,k))
    vmslabl(k) = sum(v0(2:i1,2:j1,k))
    um2slabl(k) = sum(u0(2:i1,2:j1,k)**2)
    vm2slabl(k) = sum(v0(2:i1,2:j1,k)**2)
    wm2slabl(k) = sum(w0(2:i1,2:j1,k)**2)
    tkesgsmslabl(k)=sum(e120(2:i1,2:j1,k)**2)
  end do
  call MPI_ALLREDUCE(umslabl,  umslab, k1,    MY_REAL, MPI_SUM, comm3d,mpierr)
  call MPI_ALLREDUCE(vmslabl,  vmslab, k1,    MY_REAL, MPI_SUM, comm3d,mpierr)
  call MPI_ALLREDUCE(um2slabl,um2slab, k1,    MY_REAL, MPI_SUM, comm3d,mpierr)
  call MPI_ALLREDUCE(vm2slabl,vm2slab, k1,    MY_REAL, MPI_SUM, comm3d,mpierr)
  call MPI_ALLREDUCE(wm2slabl,wm2slab, k1,    MY_REAL, MPI_SUM, comm3d,mpierr)

  call MPI_ALLREDUCE(tkesgsmslabl,tkesgsmslab,k1,MY_REAL,MPI_SUM, comm3d,mpierr)


  umslab=umslab/rslabs
  vmslab=vmslab/rslabs
  um2slab=um2slab/rslabs-umslab**2
  vm2slab=vm2slab/rslabs-vmslab**2
  wm2slab=wm2slab/rslabs
  tkesgsmslab=tkesgsmslab/rslabs

  do k=1,kmax
    tkeresmslab(k)=0.5*(um2slab(k)+vm2slab(k)+0.5*(wm2slab(k)+wm2slab(k+1)))

  end do
    tkeresmslab(k1)=0.5*(um2slab(k1)+vm2slab(k1)+0.5*wm2slab(k1))

  do k=1,k1
    fsm(k) = max(fsmin,tkesgsmslab(k)/(tkeresmslab(k)+tkesgsmslab(k)))

  end do


  deallocate( &
        umslabl,vmslabl,um2slabl,vm2slabl,wm2slabl,tkesgsmslabl, &
        umslab ,vmslab ,um2slab ,vm2slab ,wm2slab ,tkesgsmslab,&
        tkeresmslab)
  end subroutine sgsinit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sgshelpvar(particle)


    use modglobal,  only : dx,dy,dzf,zf,zh,dzh,grav,delta,timee,rdt
    use modfields,  only : dthvdz
    use modsubgriddata, only : ce1, ce2, cn
    use modsurfdata,only : thvs
    implicit none


    real    :: c0 = 6
      real    :: fs, ce, sgsdthvdz
      integer :: xbottom, ybottom, zbottom
      real    :: deltax,deltay,deltaz
    TYPE (particle_record), POINTER:: particle

    !This subroutine calculates the variables needed for the sgs-velocity equation.

    !
      xbottom = floor(particle%x-0.5)
      deltax  = particle%x-xbottom - 0.5
      ybottom = floor(particle%y-0.5)
      deltay  = particle%y-ybottom - 0.5
      zbottom = floor(particle%z-0.5)

    if (zbottom >= 1) then
      deltaz = (zh(floor(particle%z)) + dzf(floor(particle%z))*(particle%z-floor(particle%z)) - zf(zbottom))/dzh(zbottom)
      fs = (1-deltaz)*fsm(zbottom) + deltaz * fsm(zbottom+1)
      sgsdthvdz=grav/thvs*((1-deltaz) * (1-deltay) * (1-deltax) *   dthvdz(xbottom  , ybottom  , zbottom  ) + &
                          (1-deltaz) * (1-deltay) * (  deltax) *   dthvdz(xbottom+1, ybottom  , zbottom  ) + &
                          (1-deltaz) * (  deltay) * (1-deltax) *   dthvdz(xbottom  , ybottom+1, zbottom  ) + &
                          (1-deltaz) * (  deltay) * (  deltax) *   dthvdz(xbottom+1, ybottom+1, zbottom  ) + &
                          (  deltaz) * (1-deltay) * (1-deltax) *   dthvdz(xbottom  , ybottom  , zbottom+1) + &
                          (  deltaz) * (1-deltay) * (  deltax) *   dthvdz(xbottom+1, ybottom  , zbottom+1) + &
                          (  deltaz) * (  deltay) * (1-deltax) *   dthvdz(xbottom  , ybottom+1, zbottom+1) + &
                          (  deltaz) * (  deltay) * (  deltax) *   dthvdz(xbottom+1, ybottom+1, zbottom+1))

    else
    !BC for fs: fs=1 at z=0
     fs = fsm(zbottom+1)

    !fs = (1-deltat)* fsm(zbottom+1) + deltat * fs0(zbottom+1)
    !BC for dthvdz: dthvdz = cst at z=0
      sgsdthvdz =grav/thvs*((1-deltay) * (1-deltax) *    dthvdz(xbottom  , ybottom  , 1)  + &
                            (1-deltay) * (  deltax) *    dthvdz(xbottom+1, ybottom  , 1)  + &
                            (  deltay) * (1-deltax) *    dthvdz(xbottom  , ybottom+1, 1)  + &
                            (  deltay) * (  deltax) *    dthvdz(xbottom+1, ybottom+1, 1))

    end if
    sigma2_new = (2./3.)* sgstke(particle%x, particle%y, particle%z)


    if (1.5*cn**2*sigma2_new > delta(zbottom+1)**2*sgsdthvdz) then
    !unstable
      ce  = ce1 + ce2
      fce = 1.5*fs*c0*ce*sqrt(1.5*sigma2_new)/delta(zbottom+1)

    else
    !stable
      ce  = ce1 + ce2/delta(zbottom+1)*cn*sqrt(sigma2_new/sgsdthvdz)
      fce = 1.5*fs*c0*ce/cn*sqrt(sgsdthvdz)

    end if


    dsigma2dx_sgs = 2./(3*dx)*(sgstke(real(xbottom+1.5),particle%y,particle%z   ) - &
                              sgstke(real(xbottom+0.5)  ,particle%y,particle%z))
    dsigma2dy_sgs = 2./(3*dy)*(sgstke(particle%x, real(ybottom+1.5),particle%z) - &
                              sgstke(particle%x, real(ybottom+0.5)  ,particle%z))
    dsigma2dz_sgs = 2./(3*dzh(zbottom+1))* &
                             (sgstke(particle%x, particle%y, real(zbottom+1.5) ) - &
                              sgstke(particle%x, particle%y, real(zbottom+0.5)))

    dsigma2dt_sgs = (log(sigma2_new)-log(particle%sigma2_sgs))/rdt

    dsigma2dt_sgs = sign(1.,dsigma2dt_sgs)*min(abs(dsigma2dt_sgs),1./rdt)

    particle%sigma2_sgs = sigma2_new

  end subroutine sgshelpvar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function sgstke(x, y, z)

    use modglobal, only : e12min,i1,j1
    use modfields, only : e120
    implicit none

    real, intent(in) :: x, y, z

    !local::
    integer :: xbottom, ybottom, zbottom

    real :: sgstke, deltax, deltay, deltaz

    xbottom = floor(x - 0.5)
    ybottom = floor(y - 0.5)
    zbottom = floor(z - 0.5)
    deltax = x - 0.5 - xbottom
    deltay = y - 0.5 - ybottom
    deltaz = z - 0.5 - zbottom

    !interpolation of 'tke'.  For 't' = 'timee'
    xbottom = max(min(xbottom,i1),1)
    ybottom = max(min(xbottom,j1),1)
    if (zbottom == 0)  then
      sgstke =     (2-deltaz) * (1-deltay) * (1-deltax) * (  e120(xbottom    , ybottom    , 1)**2) + &
          &        (2-deltaz) * (1-deltay) * (  deltax) * (  e120(xbottom + 1, ybottom    , 1)**2) + &
          &        (2-deltaz) * (  deltay) * (1-deltax) * (  e120(xbottom    , ybottom + 1, 1)**2) + &
          &        (2-deltaz) * (  deltay) * (  deltax) * (  e120(xbottom + 1, ybottom + 1, 1)**2) + &
          &        (deltaz-1) * (1-deltay) * (1-deltax) *    e120(xbottom    , ybottom    , 2)**2  + &
          &        (deltaz-1) * (1-deltay) * (  deltax) *    e120(xbottom + 1, ybottom    , 2)**2  + &
          &        (deltaz-1) * (  deltay) * (1-deltax) *    e120(xbottom    , ybottom + 1, 2)**2  + &
          &        (deltaz-1) * (  deltay) * (  deltax) *    e120(xbottom + 1, ybottom + 1, 2)**2
    else
      sgstke =     (1-deltaz) * (1-deltay) * (1-deltax) *   e120(xbottom    , ybottom    , zbottom    )**2 + &
          &        (1-deltaz) * (1-deltay) * (  deltax) *   e120(xbottom + 1, ybottom    , zbottom    )**2 + &
          &        (1-deltaz) * (  deltay) * (1-deltax) *   e120(xbottom    , ybottom + 1, zbottom    )**2 + &
          &        (1-deltaz) * (  deltay) * (  deltax) *   e120(xbottom + 1, ybottom + 1, zbottom    )**2
      if(deltaz>0) then
        sgstke = sgstke + &
          &        (  deltaz) * (1-deltay) * (1-deltax) *   e120(xbottom    , ybottom    , zbottom + 1)**2 + &
          &        (  deltaz) * (1-deltay) * (  deltax) *   e120(xbottom + 1, ybottom    , zbottom + 1)**2 + &
          &        (  deltaz) * (  deltay) * (1-deltax) *   e120(xbottom    , ybottom + 1, zbottom + 1)**2 + &
          &        (  deltaz) * (  deltay) * (  deltax) *   e120(xbottom + 1, ybottom + 1, zbottom + 1)**2
      end if
    end if
    sgstke = max(sgstke,e12min**2)

  end function sgstke

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function gauss1(idum)



    !this function gives a random number from a gaussian distribution
    integer (KIND=selected_int_kind(10)):: idum, iset
    real :: gauss1, fac, gset, rsq, v1, v2
    save iset, gset
    data iset /0/
    rsq = 0.
    v1 = 0.
    v2 = 0.
    if (iset == 0) then
      do while (rsq >= 1 .or. rsq == 0)
        v1 = 2. * ran1(idum)-1.
        v2 = 2. * ran1(idum)-1.
        if (v1 == v2) write(*,*) 'v1 = v2'
        rsq = v1 * v1 + v2 * v2
      end do
      fac = sqrt(-2. * log(rsq) / rsq)
      gset = v1 * fac
      gauss1 = v2 * fac
      iset = 1
    else
      gauss1 = gset
      iset = 0
    end if
    return
  end function gauss1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function ran1(idum)



    !this function gives a random number from a uniform distribution (0..1)
    integer, parameter :: ntab = 32
    integer (KIND=selected_int_kind(10)):: idum, ia, im, iq, ir, iv(ntab), iy, ndiv,threshold=1
    real :: ran1, am, eps1, rnmx

    integer :: j, k
    save iv, iy

    ia = 16807
    im = 2147483647
    am = 1. / real(im)
    iq = 127773
    ir = 2836
    ndiv = 1 +  (im-1)/real(ntab)
    eps1 = 1.2E-7
    rnmx = 1. - eps1

    data iv /ntab*0/ , iy /0/

    if (idum <= 0 .OR. iy == 0 ) then
      idum = max(idum,threshold)
      do j = ntab + 8, 1, -1
        k = idum / real(iq)
        idum = ia * (idum - k * iq) - ir * k
        if (idum < 0) idum = idum + im
        if (j <= ntab ) iv(j) = idum
      end do
      iy = iv(1)

    end if
    k = idum / real(iq)
    idum = ia * (idum - k * iq) - ir * k
    if (idum <= 0) idum = idum + im
    j = 1 + iy / real(ndiv)
    iy = iv(j)
    iv(j) = idum
    ran1 = min(am*iy,rnmx)

      return
  end function ran1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine checkbound(particle)

    use modglobal, only : imax,k1
    implicit none

    type (particle_record), pointer:: particle

    !this subroutine also reflects the velocities at the bottom.

    particle%x      = modulo(particle%x-2,real(imax))+2
    particle%x_prev = modulo(particle%x_prev-2,real(imax))+2

    if (lpartsgs) then
      if (particle%z >= k1) then
          particle%wsgs = - abs(particle%wsgs)
      elseif (particle%z < 1) then
          particle%wsgs =  abs(particle%wsgs)
      end if
    end if

    if (particle%z >= k1) then
      particle%z = k1-0.0001
      particle%wres = - abs(particle%wres)
    elseif (particle%z < 1.01) then
      particle%z = abs(particle%z-1.01)+1.01
      particle%wres =  abs(particle%wres)
    end if

    if (particle%z_prev >= k1) then
      particle%z_prev = k1-0.0001
    elseif (particle%z_prev < 1.01) then
      particle%z_prev = abs(particle%z_prev-1.01)+1.01
    end if

  end subroutine checkbound

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function scalintp(x, y, z,field,field0)
    !interpolation of 'field' for 't' = 'timee'

    use modglobal, only :
    implicit none


    real, intent(in) :: x, y, z
    real, dimension(:,:,:), intent(in) :: field
    real, intent(in) :: field0

    !local::
    integer :: xbottom, ybottom, zbottom

    real :: scalintp, deltax, deltay, deltaz

    xbottom = floor(x - 0.5)
    ybottom = floor(y - 0.5)
    zbottom = floor(z - 0.5)

    deltax = x - 0.5 - xbottom
    deltay = y - 0.5 - ybottom
    deltaz = z - 0.5 - zbottom

    if (zbottom == 0)  then
      scalintp =     (1-deltaz) * (1-deltay) * (1-deltax) * (2*field0-field(xbottom    , ybottom    , 1)) + &
          &        (1-deltaz) * (1-deltay) * (  deltax) *   (2*field0-field(xbottom + 1, ybottom    , 1)) + &
          &        (1-deltaz) * (  deltay) * (1-deltax) *   (2*field0-field(xbottom    , ybottom + 1, 1)) + &
          &        (1-deltaz) * (  deltay) * (  deltax) *   (2*field0-field(xbottom + 1, ybottom + 1, 1)) + &
          &        (  deltaz) * (1-deltay) * (1-deltax) *   field(xbottom    , ybottom    , 1) + &
          &        (  deltaz) * (1-deltay) * (  deltax) *   field(xbottom + 1, ybottom    , 1) + &
          &        (  deltaz) * (  deltay) * (1-deltax) *   field(xbottom    , ybottom + 1, 1) + &
          &        (  deltaz) * (  deltay) * (  deltax) *   field(xbottom + 1, ybottom + 1, 1)
    else
      scalintp =     (1-deltaz) * (1-deltay) * (1-deltax) *   field(xbottom    , ybottom    , zbottom  ) + &
          &        (1-deltaz) * (1-deltay) * (  deltax) *   field(xbottom + 1, ybottom    , zbottom    ) + &
          &        (1-deltaz) * (  deltay) * (1-deltax) *   field(xbottom    , ybottom + 1, zbottom    ) + &
          &        (1-deltaz) * (  deltay) * (  deltax) *   field(xbottom + 1, ybottom + 1, zbottom    ) + &
          &        (  deltaz) * (1-deltay) * (1-deltax) *   field(xbottom    , ybottom    , zbottom + 1) + &
          &        (  deltaz) * (1-deltay) * (  deltax) *   field(xbottom + 1, ybottom    , zbottom + 1) + &
          &        (  deltaz) * (  deltay) * (1-deltax) *   field(xbottom    , ybottom + 1, zbottom + 1) + &
          &        (  deltaz) * (  deltay) * (  deltax) *   field(xbottom + 1, ybottom + 1, zbottom + 1)
    end if

  end function scalintp

!*****************************************************************************


  subroutine partcommunicate

    use modglobal, only : jmax,j2
    use modmpi,    only : comm3d,mpi_integer,mpi_status_size,mpierr,my_real,nbrtop,nbrbottom
    implicit none

    integer:: ii, n
    integer:: nrtosouth,nrtonorth
    integer:: nrfrsouth,nrfrnorth
    integer:: status(MPI_STATUS_SIZE)
    real, allocatable, dimension(:) :: buffsend, buffrecv
    TYPE (particle_record), POINTER:: particle,ptr

    nrtonorth = 0
    nrtosouth = 0

    particle => head
    do while( associated(particle) )
      if( particle%y >= j2 ) nrtonorth = nrtonorth + 1
      if( particle%y < 2 ) nrtosouth = nrtosouth + 1
      particle => particle%next
    end do

    call MPI_SENDRECV(nrtonorth,1,MPI_INTEGER,nbrtop,4, &
                      nrfrsouth,1,MPI_INTEGER,nbrbottom,4, &
                      comm3d, status, mpierr)

    call MPI_SENDRECV(nrtosouth,1,MPI_INTEGER,nbrbottom,5, &
                      nrfrnorth,1,MPI_INTEGER,nbrtop,5, &
                      comm3d, status, mpierr)

    if( nrtonorth > 0 ) allocate(buffsend(nrpartvar*nrtonorth))
    if( nrfrsouth > 0 ) allocate(buffrecv(nrpartvar*nrfrsouth))

    if( nrtonorth > 0 ) then
      particle => head
      ii = 0
      do while( associated(particle) )
        if( particle%y >= j2 ) then
          particle%y      = particle%y - jmax
          particle%y_prev = particle%y_prev - jmax
          call partfillbuffer(particle, buffsend(ii+1:ii+nrpartvar),ii,.true.)
          ptr => particle
          particle => particle%next
          call particle_delete(ptr)
          ii=ii+nrpartvar
        else
          particle => particle%next
        end if
      end do
    end if

    call MPI_SENDRECV(buffsend,nrpartvar*nrtonorth,MY_REAL,nbrtop,6, &
                      buffrecv,nrpartvar*nrfrsouth,MY_REAL,nbrbottom,6, &
                      comm3d, status, mpierr)

    ii = 0
    do n = 1,nrfrsouth
      call particle_add(particle)
      call partfillbuffer(particle, buffrecv(ii+1:ii+nrpartvar),ii,.false.)
      ii=ii+nrpartvar
    end do

    if( nrtonorth > 0 ) deallocate(buffsend)
    if( nrfrsouth > 0 ) deallocate(buffrecv)

    if( nrtosouth > 0 ) allocate(buffsend(nrpartvar*nrtosouth))
    if( nrfrnorth > 0 ) allocate(buffrecv(nrpartvar*nrfrnorth))

    if( nrtosouth > 0 ) then
      particle => head
      ii = 0
      do while( associated(particle) )
        if( particle%y < 2 ) then

          particle%y = particle%y + jmax
          particle%y_prev = particle%y_prev + jmax

          call partfillbuffer(particle, buffsend(ii+1:ii+nrpartvar),ii,.true.)

          ptr => particle
          particle => particle%next
          call particle_delete(ptr)
          ii=ii+nrpartvar
        else
          particle => particle%next
        end if
      end do

    end if

    call MPI_SENDRECV(buffsend,nrpartvar*nrtosouth,MY_REAL,nbrbottom,7, &
                      buffrecv,nrpartvar*nrfrnorth,MY_REAL,nbrtop,7, &
                      comm3d, status, mpierr)
    ii = 0
    do n = 1,nrfrnorth
      particle => head

      call particle_add(particle)
      call partfillbuffer(particle, buffrecv(ii+1:ii+nrpartvar),ii,.false.)
      ii=ii+nrpartvar
    end do

    if( nrtosouth > 0 ) deallocate(buffsend)
    if( nrfrnorth > 0 ) deallocate(buffrecv)

  end subroutine partcommunicate
!*****************************************************************************
  subroutine partfillbuffer(particle, buffer, n,send)
    implicit none


    logical,intent(in)::send
    integer,intent(in) :: n
    real,dimension(n+1:n+nrpartvar) :: buffer
    TYPE (particle_record), POINTER:: particle

    if (send) then
        buffer(n+ipunique)    = particle%unique
        buffer(n+ipx)         = particle%x
        buffer(n+ipy)         = particle%y
        buffer(n+ipz)         = particle%z
        buffer(n+ipures)      = particle%ures
        buffer(n+ipvres)      = particle%vres
        buffer(n+ipwres)      = particle%wres
        buffer(n+ipusgs)      = particle%usgs
        buffer(n+ipvsgs)      = particle%vsgs
        buffer(n+ipwsgs)      = particle%wsgs
        buffer(n+ipusgs_prev) = particle%usgs_prev
        buffer(n+ipvsgs_prev) = particle%vsgs_prev
        buffer(n+ipwsgs_prev) = particle%wsgs_prev
        buffer(n+ipxstart)    = particle%xstart
        buffer(n+ipystart)    = particle%ystart
        buffer(n+ipzstart)    = particle%zstart
        buffer(n+iptsart)     = particle%tstart
        buffer(n+ipartstep)   = particle%partstep
       if (intmeth == irk3 .or. lpartsgs) then
          buffer(n+ipxprev) = particle%x_prev
          buffer(n+ipyprev) = particle%y_prev
          buffer(n+ipzprev) = particle%z_prev
        end if
        if (lpartsgs) then
          buffer(n+ipsigma2_sgs)=particle%sigma2_sgs
        end if
    else
        particle%unique       = buffer(n+ipunique)
        particle%x            = buffer(n+ipx)
        particle%y            = buffer(n+ipy)
        particle%z            = buffer(n+ipz)
        particle%ures         = buffer(n+ipures)
        particle%vres         = buffer(n+ipvres)
        particle%wres         = buffer(n+ipwres)
        particle%usgs         = buffer(n+ipusgs)
        particle%vsgs         = buffer(n+ipvsgs)
        particle%wsgs         = buffer(n+ipwsgs)
        particle%usgs_prev    = buffer(n+ipusgs_prev)
        particle%vsgs_prev    = buffer(n+ipvsgs_prev)
        particle%wsgs_prev    = buffer(n+ipwsgs_prev)
        particle%xstart       = buffer(n+ipxstart)
        particle%ystart       = buffer(n+ipystart)
        particle%zstart       = buffer(n+ipzstart)
        particle%tstart       = buffer(n+iptsart)
        particle%partstep     = buffer(n+ipartstep)
        if (intmeth == irk3 .or. lpartsgs) then
          particle%x_prev =  buffer(n+ipxprev)
          particle%y_prev =  buffer(n+ipyprev)
          particle%z_prev =  buffer(n+ipzprev)
        end if
        if (lpartsgs) then
          particle%sigma2_sgs=buffer(n+ipsigma2_sgs)
        end if
    end if



    end subroutine partfillbuffer
!*****************************************************************************
  subroutine particle_initlist()
    implicit none
  ! LOCAL

    nullify(head) ! no particles in list
    nullify(tail) ! no particles in list
    nplisted = 0

  end subroutine particle_initlist
  !*****************************************************************************
  subroutine particle_quitlist()
     implicit none
 ! LOCAL

    do while (associated(tail))
        call particle_delete(tail)
    end do

  end subroutine particle_quitlist
  !*****************************************************************************
  subroutine particle_add(ptr)
    implicit none
  TYPE (particle_record), POINTER:: ptr
  ! LOCAL
  TYPE (particle_record), POINTER:: new_p

    if( .not. associated(head) ) then !empty list
      allocate(head)
      tail => head
      nullify(head%prev)
    else
      allocate(tail%next)
      new_p => tail%next
      new_p%prev => tail
      tail => new_p
    end if

    nplisted = nplisted + 1
    nullify(tail%next)
    ptr => tail

  end SUBROUTINE particle_add
  !*****************************************************************************
  subroutine particle_delete(ptr)
    implicit none
  TYPE (particle_record), POINTER:: ptr
  ! LOCAL
  TYPE (particle_record), POINTER:: next_p,prev_p,cur_p

    cur_p => ptr
    if( .not. associated(cur_p) ) then         !error in calling ptr
      write(6,*) 'WARNING: cannot delete empty pointer'
      return
    end if

    if( .not. associated(head) ) then         !empty list
      write(6,*) 'WARNING: cannot delete elements in an empty list'
    else
      if( .not. associated(cur_p%next) ) then   ! last in list
        if( .not. associated(cur_p%prev) ) then ! last element
          nullify(head)
          nullify(tail)
        else
          tail => cur_p%prev
          nullify(tail%next)
        end if
      else
        if( .not. associated(cur_p%prev) ) then ! first in list
          if( .not. associated(cur_p%next) ) then !last element
            nullify(head)
            nullify(tail)
          else
            head => cur_p%next
            nullify(head%prev)
          end if
        else
          next_p => cur_p%next
          prev_p => cur_p%prev
          next_p%prev => prev_p
          prev_p%next => next_p
        end if
      end if

      nplisted = nplisted - 1
      deallocate(cur_p)
    end if

  end subroutine particle_delete
  !*****************************************************************************
  subroutine particle_list()
    implicit none
  ! LOCAL
  TYPE (particle_record), POINTER:: ptr
  integer :: n

    write(6,*) 'particle list..'
    if( .not. associated(head) ) then
      write(6,*) 'empty list'
    else
      ptr => head
      n = 0
      do while (associated(ptr))
        n = n + 1
        write(6,*) n,ptr%unique
        ptr => ptr%next
      end do
    end if

  end subroutine particle_list
  !*****************************************************************************

end module modparticles

