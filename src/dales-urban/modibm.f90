!!> \file modibm.f90
!!!  adds forcing terms for immersed boundaries
!
!>
!!  \author Jasper Tomas, TU Delft, June 4th 2015
!
!
module modibm
use modibmdata
implicit none
save
  public :: createwalls,ibmshear,xwallshear,ywallshearplus,ywallshearmin,&
            zwallshear,ibmnorm,nearwall,nearwall2d,ibmforce

contains
    subroutine createwalls
   ! #wallsglobal is in global indices!
   ! xwallsglobal: 2nd index: iloc,jstart,jend,kstart,kend)
   ! ywallsglobal: 2nd index: jloc,istart,iend,kstart,kend)
   ! zwallsglobal: 2nd index: kloc,istart,iend,jstart,jend)
      use modglobal, only : ib,ie,jb,je,jgb,jge,kb,ke,jmax,nxwalls,nywalls,nzwalls,nblocks,&
                            nsv,cexpnr,ifinput,libm,ih,kh,lreadmean
      use modsurfdata, only : thls,qts
      use modfields, only : sv0,svm,thl0,thlm,qtp,qt0
      use modgenstat,only : lcanopy
      use modmpi,    only : myid,comm3d,mpierr,mpi_integer,nprocs,cmyid
      integer n,nn, pn, mn,jbeg,jend,nxn,nxs,nyn,nzn,nzs,iu,il,ju,jl,ku,kl,sc,&
              i,k
      character(80) chmess,name2

      if (libm==.false.) return      

      allocate(xwallsglobal(nxwalls,5))
      allocate(ywallsglobal(nywalls,5))
      allocate(zwallsglobal(nzwalls,5))
!      allocate(block(1,6))
      allocate(block(nblocks,6))
      allocate(ibmxforce(ib-ih:ie+ih,kb-kh:ke+kh))
      allocate(ibmxforcevol(ib-ih:ie+ih,kb-kh:ke+kh))
      allocate(ibmxforcevolp(ib-ih:ie+ih,kb-kh:ke+kh))
      ibmxforce(:,:) = 0.
      ibmxforcevol(:,:) = 0.
      ibmxforcevolp(:,:) = 0.

      ! Read average ibm force in x-direction to file
      if (lreadmean ==.true. .and. lcanopy == .true.) then
        write(6,*) 'Reading ibmx_XXX.XXX, proc = ',myid
        name2 = 'ibmx_   .'
        name2( 6:8 ) = cmyid
        name2(10:12) = cexpnr
        open(unit=ifinput,file=name2,form='unformatted')
        write(ifinput) ((ibmxforcevol(i,k),i=ib-ih,ie+ih),k=kb-kh,ke+kh)
        close (ifinput)
      end if


      ! read global blocks
!      if(myid==0 .and. libm==.true.)then
      if(myid==0) then
        if (nblocks>0) then
          open (ifinput,file='blocks.inp.'//cexpnr)
          read (ifinput,'(a80)') chmess
          read (ifinput,'(a80)') chmess
          do n=1,nblocks
            read (ifinput,*) &
                  block(n,1), &
                  block(n,2), &
                  block(n,3), &
                  block(n,4), &
                  block(n,5), &
                  block(n,6) 
          end do
          write (6,*) 'Block number,   il, iu, jl, ju, kl, ku '
          do n=1,nblocks
            write (6,*) &
                  n , &
                  block(n,1), &
                  block(n,2), &
                  block(n,3), &
                  block(n,4), &
                  block(n,5), &
                  block(n,6)                 
          end do

        end if
      end if ! end if myid==0

      call MPI_BCAST(block ,6*nblocks,MPI_INTEGER ,0,comm3d,mpierr)


  ! Create local walls
 

!needs attention  ILS13, 31.07.15


! For all blocks set the internal concentrations to zero and internal
! temperature to thls
  
      do n=1,nblocks 
        il = block(n,1)
        iu = block(n,2)
        kl = block(n,5)
        ku = block(n,6)
        jl = block(n,3)-myid*jmax
        ju = block(n,4)-myid*jmax
        if (ju < jb .or. jl > je) then 
          cycle
        else
          if (ju > je) ju=je
          if (jl < jb) jl=jb
          do sc=1,nsv
            sv0(il:iu,jl:ju,kl:ku,sc) = 0.
            svm(il:iu,jl:ju,kl:ku,sc) = 0.
          end do
          thl0(il:iu,jl:ju,kl:ku) = thls
          thlm(il:iu,jl:ju,kl:ku) = thls
        end if 
      end do

! First determine no. of local wall 

  ! First no. of local xwalls
      nxwallsnorm = 0   ! no. of local x-walls
      do n=1,nxwalls         
        if (xwallsglobal(n,3) < jb+myid*jmax) then
          ! no x-wall in the range of this proc
          cycle
        elseif (xwallsglobal(n,2) > je+myid*jmax ) then
          ! no x-wall in the range of this proc
          cycle
        else 
          nxwallsnorm = nxwallsnorm + 1
        end if 
      end do

    ! Then determine no. of local x-wall shear plane walls (staggered grid, so shear effect can be on next proc)
      nxwallsshear = 0   ! no. of x-walls with shear planes in this proc
      do n=1,nxwalls
        if (xwallsglobal(n,3) < jb-1+myid*jmax) then
          ! no x-wall shear plane in the range of this proc
          cycle
        elseif (xwallsglobal(n,2) > je+myid*jmax ) then
          ! no x-wall shear plane in the range of this proc
          cycle
        else
          nxwallsshear = nxwallsshear + 1
        end if
      end do

    ! Determine no. of local ywalls
      nywallsp = 0
      nywallsm = 0
      do n=1,nywalls
        if (ywallsglobal(n,1)-myid*jmax < je+1 .and. ywallsglobal(n,1)-myid*jmax > jb) then
          nywallsp = nywallsp + 1
          nywallsm = nywallsm + 1
        elseif (ywallsglobal(n,1)-myid*jmax == je+1) then
          nywallsm = nywallsm + 1
        elseif (ywallsglobal(n,1)-myid*jmax == jb) then
          nywallsp = nywallsp + 1
        end if
      end do
      nywallsnorm = nywallsp

    ! Determine no. of local z-walls. (Here z-velocity is zero)
      nzwallsnorm = 0   ! no. of local z-walls
      do n=1,nzwalls
        if (zwallsglobal(n,5) < jb+myid*jmax) then
          ! no z-wall in the range of this proc
          cycle
        elseif (zwallsglobal(n,4) > je+myid*jmax ) then
          ! no z-wall in the range of this proc
          cycle
        else
          nzwallsnorm = nzwallsnorm + 1
        end if
      end do

    ! Then determine no. of local z-wall shear plane walls (staggered grid, so shear effect can be on next proc)
      nzwallsshear = 0   ! no. of z-walls with shear planes in this proc
      do n=1,nzwalls
        if (zwallsglobal(n,5) < jb-1+myid*jmax) then
          ! no z-wall shear plane in the range of this proc
          cycle
        elseif (zwallsglobal(n,4) > je+myid*jmax ) then
          ! no z-wall shear plane in the range of this proc
          cycle
        else
          nzwallsshear = nzwallsshear + 1
        end if
      end do


    
  ! Now do the same from blocks
    ! first zero normal velocity components
      do n=1,nblocks                          ! first x and z walls
        if (block(n,4) < jb+myid*jmax) then
          ! no x-wall in the range of this proc
          cycle
        elseif (block(n,3) > je+myid*jmax ) then
          ! no x-wall in the range of this proc
          cycle
        else
          nxwallsnorm = nxwallsnorm + block(n,2)-block(n,1)+2  ! add appropriate no. of x-wall normal
          nzwallsnorm = nzwallsnorm + block(n,6)-block(n,5)+2  ! add appropriate no. of z-wall normal
        end if
      end do

      do n=1,nblocks                          ! then y-walls
        if (block(n,4) < jb-1+myid*jmax) then
          ! no x-wall in the range of this proc
          cycle
        elseif (block(n,3) > je+myid*jmax ) then
          ! no x-wall in the range of this proc
          cycle
        else
          jend = block(n,4)-myid*jmax+1
          jbeg = block(n,3)-myid*jmax
!          if (jend > je) jend = je-1
          if (jend > je) jend = je
          if (jbeg < jb) jbeg = jb
          nywallsnorm = nywallsnorm + jend - jbeg +1
        end if
      end do


    ! Then shear components (first shear due to y-walls)
      do n=1,nblocks

        if (myid==0 .and. block(n,4) == jge) then             ! periodicity!
            nywallsp = nywallsp + 1                             
        elseif (block(n,3) == jgb .and. myid==nprocs-1) then  ! periodicity!
            nywallsm = nywallsm + 1
        end if

        if (block(n,4) < jb-1+myid*jmax) then
          ! no y-wall shear plane in the range of this proc
          cycle
        elseif (block(n,3) > je+1+myid*jmax ) then
          ! no y-wall shear plane in the range of this proc
          cycle
        else      ! there is/are wall(s)
!          if (block(n,4)-myid*jmax == jb-1) then
!            nywallsp = nywallsp + 1
!          elseif (block(n,3)-myid*jmax == je+1) then
!            nywallsm = nywallsm + 1
!          else
            if (block(n,3)-myid*jmax <= je+1 .and. block(n,3)-myid*jmax > jb) then
              nywallsm = nywallsm + 1
            end if
            if (block(n,4)-myid*jmax <= je .and. block(n,4)-myid*jmax >=jb) then
              nywallsm = nywallsm + 1
            end if
            if (block(n,4)-myid*jmax >= jb-1 .and. block(n,4)-myid*jmax < je) then
              nywallsp = nywallsp + 1
            end if
            if (block(n,3)-myid*jmax >= jb .and. block(n,3)-myid*jmax <= je) then
              nywallsp = nywallsp + 1
            end if

!! following block was commented, uncommented it. ILS13  03.08.2015
!            ! special case: block ends at jge
!            if (myid==0) then
!              if (block(n,4) == jge) then
!                nywallsp = nywallsp + 1
!              end if
!            end if
!            ! special case: block begins at jgb
!            if (myid==nprocs-1) then
!             if (block(n,3) == jgb) then
!               nywallsm = nywallsm + 1
!              end if
!            end if  
!


! FOUT!
!            if (block(n,3)-myid*jmax <= je+1 .and. block(n,3)-myid*jmax > jb) then
!              nywallsm = nywallsm + 1
!            end if
!            if (block(n,3)-myid*jmax <= je .and. block(n,3)-myid*jmax >= jb) then
!              nywallsp = nywallsp + 1
!            end if
!            if (block(n,4)-myid*jmax >= jb-1 .and. block(n,4)-myid*jmax < je) then
!              nywallsm = nywallsm + 1
!            end if
!            if (block(n,4)-myid*jmax >= jb .and. block(n,4)-myid*jmax <= je) then
!              nywallsp = nywallsp + 1
!            end if

!



!          end if
        end if
      end do

   ! now do x/z-walls
      do n=1,nblocks
        if (block(n,4) < jb-1+myid*jmax) then
          ! no x/z-wall shear plane in the range of this proc
          cycle
        elseif (block(n,3) > je+myid*jmax ) then
          ! no x/z-wall shear plane in the range of this proc
          cycle
        else      ! there is/are wall(s)
          nxwallsshear = nxwallsshear + 2                       ! two shear x-walls per block
          if (block(n,5) < kb+1) then
            nzwallsshear = nzwallsshear +1  ! one less because lower wall already in global BC's
          else
            nzwallsshear = nzwallsshear +2  ! add appropriate no. of z-wall normal
          end if
        end if
      end do



      allocate(xwallsnorm(nxwallsnorm,5))          ! number of 'local' x-walls (zero x-velocity) is now known
      allocate(xwallsshear(nxwallsshear,5))        ! number of 'local' x-wall shear planes is now known
      allocate(ywallsp(nywallsp,5))                ! number of 'plus' walls is now known
      allocate(ywallsnorm(nywallsnorm,5))          ! number of zero normal velocity walls now known
      allocate(ywallsm(nywallsm,5))                ! number of 'min'  walls is now known
      allocate(zwallsnorm(nzwallsnorm,5))          ! number of 'local' z-walls (zero z-velocity) is now known
      allocate(zwallsshear(nzwallsshear,5))        ! number of 'local' z-wall shear planes is now known 
     
      allocate(bcthltypex(nxwallsshear))
      allocate(bcthltypeyp(nywallsp))
      allocate(bcthltypeym(nywallsm))
      allocate(bcthltypez(nzwallsshear))
      allocate(bcthlvaluex(nxwallsshear))
      allocate(bcthlvalueyp(nywallsp))
      allocate(bcthlvalueym(nywallsm))
      allocate(bcthlvaluez(nzwallsshear))

      allocate(bcqttypex(nxwallsshear))
      allocate(bcqttypeyp(nywallsp))
      allocate(bcqttypeym(nywallsm))
      allocate(bcqttypez(nzwallsshear))
      allocate(bcqtvaluex(nxwallsshear))
      allocate(bcqtvalueyp(nywallsp))
      allocate(bcqtvalueym(nywallsm))
      allocate(bcqtvaluez(nzwallsshear))



!  ! fixed temperature boundary conditions for temperature
!      bcthltypex = 1
!      bcthltypeyp = 1
!      bcthltypeym = 1
!      bcthltypez = 1
!      bcthlvaluex = thls
!      bcthlvalueyp = thls
!      bcthlvalueym = thls
!      bcthlvaluez = thls

 ! zero flux boundary conditions for temperature
      bcthltypex = 0
      bcthltypeyp = 0
      bcthltypeym = 0
      bcthltypez = 0     !fixed temperature at top (and bottom!) wall 
     bcthlvaluex = 0.0
      bcthlvalueyp = 0.0
      bcthlvalueym = 0.0
!      bcthlvaluez = 0.01
      bcthlvaluez = 0.0   ! surface value is used at top wall


 ! zero flux boundary conditions for moisture
      bcqttypex = 0
      bcqttypeyp = 0
      bcqttypeym = 0
      bcqttypez = 0     !fixed temperature at top (and bottom!) wall 
      bcqtvaluex = 0.0
      bcqtvalueyp = 0.0
      bcqtvalueym = 0.0
!     bcqtvaluez = 0.0
      bcqtvaluez = 0.0   ! surface value is used at top wall











    ! Do the loop again and determine again if wall exists on this proc
    ! and what the j-range is
      nxn = 0    ! index for 'local' x-walls
      do n=1,nxwalls
        if (xwallsglobal(n,3) < jb+myid*jmax) then
          ! no x-wall in the range of this proc
          cycle
        elseif (xwallsglobal(n,2) > je+myid*jmax ) then
          ! no x-wall in the range of this proc
          cycle
        else  ! There is a wall in the range of this proc
            nxn = nxn + 1
            xwallsnorm(nxn,:) = xwallsglobal(n,:)
            xwallsnorm(nxn,2) = xwallsnorm(nxn,2) - myid*jmax
            xwallsnorm(nxn,3) = xwallsnorm(nxn,3) - myid*jmax
          if (xwallsnorm(nxn,2) < jb) then   ! correct the starting j-index of the wall to minimum local value
            xwallsnorm(nxn,2) = jb
          end if
          if (xwallsnorm(nxn,3) > je) then  ! correct the starting j-index of the wall to maximum local value
            xwallsnorm(nxn,3) = je 
          end if
        end if
      end do
 
    ! Do the 'shear loop again and determine again if wall exists on this proc
    ! and what the j-range is (minor difference with above loop!!)
      nxs = 0    ! index for 'local' x-walls
      do n=1,nxwalls
        if (xwallsglobal(n,3) < jb-1+myid*jmax) then
          ! no x-wall in the range of this proc
          cycle
        elseif (xwallsglobal(n,2) > je+myid*jmax ) then
          ! no x-wall in the range of this proc
          cycle
        else  ! There is a wall in the range of this proc
            nxs = nxs + 1
            xwallsshear(nxs,1:5) = xwallsglobal(n,1:5)
            xwallsshear(nxs,2) = xwallsshear(nxs,2) - myid*jmax
            xwallsshear(nxs,3) = xwallsshear(nxs,3) - myid*jmax
            xwallsshear(nxs,5) = xwallsshear(nxs,5) +1     ! 'shear' components is always one more than 'normal' components
          if (xwallsshear(nxs,2) < jb) then   ! correct the starting j-index of the wall to minimum local value
            xwallsshear(nxs,2) = jb
          end if
          if (xwallsshear(nxs,3) >= je) then   ! correct the starting j-index of the wall to maximum local value
            xwallsshear(nxs,3) = je
          else
            xwallsshear(nxs,3) = xwallsshear(nxs,3) +1     ! 'shear' components is always one more than 'normal' components
          end if
        end if
      end do


    ! Do the loop again and determine 'plus' and 'min' y-walls 
      pn = 0  ! index for 'plus' array
      mn = 0  ! index for 'min'  array
      do n=1,nywalls
        if (ywallsglobal(n,1)-myid*jmax < je+1 .and. ywallsglobal(n,1)-myid*jmax > jb) then
          pn = pn + 1 
          mn = mn + 1
          ywallsp(pn,:) = ywallsglobal(n,:) 
          ywallsp(pn,1) = ywallsp(pn,1) -myid*jmax
          ywallsm(mn,:) = ywallsglobal(n,:)
          ywallsm(mn,1) = ywallsm(mn,1) -myid*jmax
          ywallsnorm(pn,:)= ywallsglobal(n,:)
          ywallsnorm(pn,1) = ywallsnorm(pn,1) -myid*jmax
        elseif (ywallsglobal(n,1)-myid*jmax == je+1) then
          mn = mn + 1
          ywallsm(mn,:) = ywallsglobal(n,:)
          ywallsm(mn,1) = ywallsm(mn,1) -myid*jmax
        elseif (ywallsglobal(n,1)-myid*jmax == jb) then
          pn = pn + 1 
          ywallsp(pn,:) = ywallsglobal(n,:)
          ywallsp(pn,1) = ywallsp(pn,1) -myid*jmax
          ywallsnorm(pn,:)= ywallsglobal(n,:)
          ywallsnorm(pn,1) = ywallsnorm(pn,1) -myid*jmax
        end if
      end do
      ywallsp(1:pn,3) = ywallsp(1:pn,3) +1    ! 'shear' components is always one more than 'normal' components
      ywallsp(1:mn,5) = ywallsp(1:mn,5) +1    ! 'shear' components is always one more than 'normal' components
 
    
 ! z-walls

    ! Do the loop again and determine again if wall exists on this proc
    ! and what the j-range is
      nzn = 0    ! index for 'local' xz-walls
      do n=1,nzwalls
        if (zwallsglobal(n,5) < jb+myid*jmax) then
          ! no z-wall in the range of this proc
          cycle
        elseif (zwallsglobal(n,4) > je+myid*jmax) then
          ! no z-wall in the range of this proc
          cycle
        else  ! There is a wall in the range of this proc
            nzn = nzn + 1
            zwallsnorm(nzn,:) =  zwallsglobal(n,:)
            zwallsnorm(nzn,4) = zwallsnorm(nzn,4) - myid*jmax
            zwallsnorm(nzn,5) = zwallsnorm(nzn,5) - myid*jmax
          if (zwallsnorm(nzn,4) < jb) then   ! correct the starting j-index of the wall to minimum local value
            zwallsnorm(nzn,4) = jb
          end if
          if (zwallsnorm(nzn,5) > je) then   ! correct the starting j-index of the wall to maximum local value
            zwallsnorm(nzn,5) = je
          end if
        end if
      end do

    ! Do the 'shear loop again and determine again if wall exists on this proc
    ! and what the j-range is (minor difference with above loop!!)
      nzs = 0    ! index for 'local' z-walls
      do n=1,nzwalls
        if (zwallsglobal(n,5) < jb-1+myid*jmax) then
          ! no z-wall in the range of this proc
          cycle
        elseif (zwallsglobal(n,4) > je+myid*jmax) then
          ! no x-wall in the range of this proc
          cycle
        else  ! There is a wall in the range of this proc
            nzs = nzs + 1
            zwallsshear(nzs,1:5) = zwallsglobal(n,1:5)
            zwallsshear(nzs,4) = zwallsshear(nzs,4) - myid*jmax
            zwallsshear(nzs,5) = zwallsshear(nzs,5) - myid*jmax
            zwallsshear(nzs,3) = zwallsshear(nzs,3) +1     ! 'shear' components is always one more than 'normal' components
          if (zwallsshear(nzs,4) < jb) then   ! correct the starting j-index of the wall to minimum local value
            zwallsshear(nzs,4) = jb
          end if
          if (zwallsshear(nzs,5) >= je) then   ! correct the starting j-index of the wall to maximum local value
            zwallsshear(nzs,5) = je
          else
            zwallsshear(nzs,5) = zwallsshear(nzs,5) +1     ! 'shear' components is always one more than 'normal' components
          end if
        end if
      end do

! Now determine all additional walls from blocks
      
    ! first zero normal velocity components
      do n=1,nblocks                                             ! for x- and z-walls
        if (block(n,4) < jb+myid*jmax) then
          ! no x/z-wall in the range of this proc
          cycle
        elseif (block(n,3) > je+myid*jmax ) then
          ! no x/z-wall in the range of this proc
          cycle
        else
          do nn=block(n,1),block(n,2)+1    ! loop over x-walls
            nxn=nxn+1                                            ! continuation of nxn! 
            xwallsnorm(nxn,1) = block(n,1)+(nn-block(n,1))       ! i-location
            xwallsnorm(nxn,4) = block(n,5)                       ! k_start
            xwallsnorm(nxn,5) = block(n,6)                       ! k_end
            xwallsnorm(nxn,2) = block(n,3)-myid*jmax             ! jstart
            xwallsnorm(nxn,3) = block(n,4)-myid*jmax             ! jend
            if (xwallsnorm(nxn,2) < jb)  xwallsnorm(nxn,2) = jb
            if (xwallsnorm(nxn,3) > je)  xwallsnorm(nxn,3) = je
          end do
          do nn=block(n,5),block(n,6)+1    ! loop over z-walls
            nzn=nzn+1
            zwallsnorm(nzn,1) = block(n,5)+(nn-block(n,5))       ! k-location
            zwallsnorm(nzn,2) = block(n,1)                       ! i_start
            zwallsnorm(nzn,3) = block(n,2)                       ! i_end
!            zwallsnorm(nzn,4) = block(n,3)                       ! j_start
!            zwallsnorm(nzn,5) = block(n,4)                       ! j_end
            zwallsnorm(nzn,4) = block(n,3)-myid*jmax             ! j_start
            zwallsnorm(nzn,5) = block(n,4)-myid*jmax             ! j_end
            if (zwallsnorm(nzn,4) < jb)  zwallsnorm(nzn,4) = jb
            if (zwallsnorm(nzn,5) > je)  zwallsnorm(nzn,5) = je
           end do
        end if
      end do

      nyn=pn
      do n=1,nblocks                                             ! Now for y-walls
        if (block(n,4) < jb-1+myid*jmax) then
          ! no y-wall in the range of this proc
          cycle
        elseif (block(n,3) > je+myid*jmax ) then
          ! no y-wall in the range of this proc
          cycle
        else
          do nn=block(n,3),block(n,4)+1    ! loop over y-walls
            if (nn-myid*jmax < jb .or. nn-myid*jmax>je) cycle
            nyn=nyn+1
            ywallsnorm(nyn,1) = block(n,3)+(nn-block(n,3))-myid*jmax  ! j-location
            ywallsnorm(nyn,2) = block(n,1)                      ! i-start
            ywallsnorm(nyn,3) = block(n,2)                      ! i-end
            ywallsnorm(nyn,4) = block(n,5)                      ! k-start
            ywallsnorm(nyn,5) = block(n,6)                      ! k-end
          end do
        end if
      end do
!      write(6,*) 'myid,nywallsnorm,nyn', myid,nywallsnorm,nyn



    ! Then shear components (first shear due to y-walls)
      do n=1,nblocks
        if (myid==0 .and. block(n,4) == jge) then
          pn = pn + 1
          ywallsp(pn,1) = jb
          ywallsp(pn,2) = block(n,1)
          ywallsp(pn,3) = block(n,2)+1
          ywallsp(pn,4) = block(n,5)
          ywallsp(pn,5) = block(n,6)+1
        elseif (block(n,3) == jgb .and. myid==nprocs-1) then
          mn = mn + 1
          ywallsm(mn,1) = je
          ywallsm(mn,2) = block(n,1)
          ywallsm(mn,3) = block(n,2)+1
          ywallsm(mn,4) = block(n,5)
          ywallsm(mn,5) = block(n,6)+1
        end if

        if (block(n,4) < jb-1+myid*jmax) then
          ! no y-wall shear plane in the range of this proc
          cycle
        elseif (block(n,3) > je+1+myid*jmax ) then
          ! no y-wall shear plane in the range of this proc
          cycle
        else      ! there is/are wall(s)
! First lower wall contributions

!            if (block(n,3)-myid*jmax < je+1 .and. block(n,3)-myid*jmax > jb) then 
            if (block(n,3)-myid*jmax <= je+1 .and. block(n,3)-myid*jmax > jb) then 
              mn = mn + 1
              if (mn > nywallsm) then
                write(6,*) 'WARNING!! mn > nywallsm'
              end if   
!              ywallsm(mn,1) = block(n,3)-myid*jmax
              ywallsm(mn,1) = block(n,3)-myid*jmax-1
              ywallsm(mn,2) = block(n,1)
              ywallsm(mn,3) = block(n,2)+1
              ywallsm(mn,4) = block(n,5)
              ywallsm(mn,5) = block(n,6)+1
            end if
!            if (block(n,3)-myid*jmax <= je .and. block(n,3)-myid*jmax >= jb) then
            if (block(n,3)-myid*jmax >= jb .and. block(n,3)-myid*jmax <= je) then
              pn = pn + 1
              if (pn > nywallsp) then
                write(6,*) 'WARNING!! pn > nywallsp'
              end if   
!              ywallsm(mn,1) = block(n,3)-myid*jmax
              ywallsp(pn,1) = block(n,3)-myid*jmax
              ywallsp(pn,2) = block(n,1)
              ywallsp(pn,3) = block(n,2)+1
              ywallsp(pn,4) = block(n,5)
              ywallsp(pn,5) = block(n,6)+1
            end if
! Then upper wall contributions

!            if (block(n,4)-myid*jmax > jb-1 .and. block(n,4)-myid*jmax < je) then
            if (block(n,4)-myid*jmax >= jb-1 .and. block(n,4)-myid*jmax < je) then
              pn = pn + 1
              if (pn > nywallsp) then
                write(6,*) 'WARNING!! pn > nywallsp'
              end if   
              ywallsp(pn,1) = block(n,4)-myid*jmax+1
              ywallsp(pn,2) = block(n,1)
              ywallsp(pn,3) = block(n,2)+1
              ywallsp(pn,4) = block(n,5)
              ywallsp(pn,5) = block(n,6)+1
            end if
!            if (block(n,4)-myid*jmax >= jb .and. block(n,4)-myid*jmax <= je) then
            if (block(n,4)-myid*jmax <= je .and. block(n,4)-myid*jmax >=jb) then
              mn = mn + 1
              if (mn > nywallsm) then
                write(6,*) 'WARNING!! mn > nywallsm'
              end if   
              ywallsm(mn,1) = block(n,4)-myid*jmax
              ywallsm(mn,2) = block(n,1)
              ywallsm(mn,3) = block(n,2)+1
              ywallsm(mn,4) = block(n,5)
              ywallsm(mn,5) = block(n,6)+1
            end if
!          end if
        end if
      end do
!      write(6,*) 'myid,nywallsp,pn',myid, nywallsp, pn
   
   ! now do x/z-walls
      do n=1,nblocks
        if (block(n,4) < jb-1+myid*jmax) then   ! if block ends at jb-1, there is only a line of v-components on this proc!
          ! no x/z-wall shear plane in the range of this proc
          cycle
        elseif (block(n,3) > je+myid*jmax ) then
          ! no x/z-wall shear plane in the range of this proc
          cycle
        else      ! there is/are wall(s)
          nxs = nxs + 1                       ! first x-wall of the block
          xwallsshear(nxs,1) = block(n,1)  
          xwallsshear(nxs,2) = block(n,3)-myid*jmax  
          xwallsshear(nxs,3) = block(n,4)-myid*jmax   
          xwallsshear(nxs,4) = block(n,5)  
          xwallsshear(nxs,5) = block(n,6) + 1                ! shear components is always one more  
!          if (xwallsshear(nxs,2) < jb) xwallsshear(nxs,2) = jb
          if (xwallsshear(nxs,2) < jb-1) xwallsshear(nxs,2) = jb-1
          if (xwallsshear(nxs,3) >= je) then
            xwallsshear(nxs,3) = je
          end if
          xwallsshear(nxs,3) = xwallsshear(nxs,3) + 1        ! shear components is always one more
          nxs = nxs + 1                       ! second x-wall of the block
          xwallsshear(nxs,1) = block(n,2)+1
          xwallsshear(nxs,2:5) = xwallsshear(nxs-1,2:5)
        ! now the z-walls (same criteria)
          if (block(n,5) > kb) then
            nzs = nzs + 1                       ! first z-wall is not on kb
            zwallsshear(nzs,1) = block(n,5)
            zwallsshear(nzs,2) = block(n,1)
            zwallsshear(nzs,3) = block(n,2) + 1                ! shear components is always one more
            zwallsshear(nzs,4) = block(n,3) -myid*jmax
            zwallsshear(nzs,5) = block(n,4) -myid*jmax 
!            if (zwallsshear(nzs,4) < jb) zwallsshear(nzs,4) = jb
            if (zwallsshear(nzs,4) < jb-1) zwallsshear(nzs,4) = jb-1
            if (zwallsshear(nzs,5) >= je) then
              zwallsshear(nzs,5) = je
            end if
            zwallsshear(nzs,5) = zwallsshear(nzs,5) + 1        ! shear components is always one more
            nzs = nzs + 1                       ! second z-wall of the block
            zwallsshear(nzs,1) = block(n,6)+1
            zwallsshear(nzs,2:5) = zwallsshear(nzs-1,2:5)
          else
            nzs = nzs + 1           ! only upper zwall, because lower wall is already in global BC
            zwallsshear(nzs,1) = block(n,6)+1
            zwallsshear(nzs,2) = block(n,1)
            zwallsshear(nzs,3) = block(n,2)+1                  ! shear components is always one more
            zwallsshear(nzs,4) = block(n,3)-myid*jmax
            zwallsshear(nzs,5) = block(n,4)-myid*jmax 
            if (zwallsshear(nzs,4) < jb) zwallsshear(nzs,4) = jb-1
            if (zwallsshear(nzs,5) >= je) then
              zwallsshear(nzs,5) = je
            end if
            zwallsshear(nzs,5) = zwallsshear(nzs,5) + 1        ! shear components is always one more
          end if
        end if
      end do 


!    write(6,*) 'Proc, nywallsp, nywallsm', myid, nywallsp,nywallsm 
!    if (myid==0) then
!      write(6,*) 'proc, ywallsm(1), j, il, iu', myid, ywallsm(1,1),ywallsm(1,2),ywallsm(1,3)
!      write(6,*) 'proc, ywallsm(2), j, il, iu', myid, ywallsm(2,1),ywallsm(2,2),ywallsm(2,3)
!      write(6,*) 'proc, ywallsp(1), j, il, iu', myid, ywallsp(1,1),ywallsp(1,2),ywallsp(1,3)
!    end if
!    if (myid==1) then
!      write(6,*) 'proc, ywallsp(1), j, il, iu', myid, ywallsp(1,1),ywallsp(1,2),ywallsp(1,3)
!      write(6,*) 'proc, ywallsp(2), j, il, iu', myid, ywallsp(2,1),ywallsp(2,2),ywallsp(2,3)
!    end if
!    if (myid==2) then
!      write(6,*) 'proc, ywallsm(1), j, il, iu', myid, ywallsm(1,1),ywallsm(1,2),ywallsm(1,3)
!      write(6,*) 'proc, ywallsm(2), j, il, iu', myid, ywallsm(2,1),ywallsm(2,2),ywallsm(2,3)
!    end if



    open(unit=11,file='wallsxn'//cmyid//'.txt',position='rewind')
    write(11,'(3A)') 'i, jl, ju, kl, ku'
    do n=1,nxwallsnorm
      write(11,'(5I)') xwallsnorm(n,1),xwallsnorm(n,2),xwallsnorm(n,3),xwallsnorm(n,4),xwallsnorm(n,5)
    end do  
 
    open(unit=11,file='wallsyn'//cmyid//'.txt',position='rewind')
    write(11,'(3A)') 'j, il, iu, kl, ku'
    do n=1,nywallsnorm
      write(11,'(5I)') ywallsnorm(n,1),ywallsnorm(n,2),ywallsnorm(n,3),ywallsnorm(n,4),ywallsnorm(n,5)
    end do  
 
    open(unit=11,file='wallszn'//cmyid//'.txt',position='rewind')
    write(11,'(3A)') 'k, il, iu, jl, ju'
    do n=1,nzwallsnorm
      write(11,'(5I)') zwallsnorm(n,1),zwallsnorm(n,2),zwallsnorm(n,3),zwallsnorm(n,4),zwallsnorm(n,5)
    end do   

    open(unit=11,file='wallsx_'//cmyid//'.txt',position='rewind')
    write(11,'(3A)') 'i, jl, ju+1, kl, ku+1'
    do n=1,nxwallsshear
      write(11,'(5I)') xwallsshear(n,1),xwallsshear(n,2),xwallsshear(n,3),xwallsshear(n,4),xwallsshear(n,5)
    end do   

    open(unit=11,file='wallsym'//cmyid//'.txt',position='rewind')
    write(11,'(3A)') 'j, il, iu+1, kl, ku+1'
    do n=1,nywallsm
      write(11,'(5I)') ywallsm(n,1),ywallsm(n,2),ywallsm(n,3),ywallsm(n,4),ywallsm(n,5)
    end do  
 
    open(unit=11,file='wallsyp'//cmyid//'.txt',position='rewind')
    write(11,'(3A)') 'j, il, iu+1, kl, ku+1'
    do n=1,nywallsp
      write(11,'(5I)') ywallsp(n,1),ywallsp(n,2),ywallsp(n,3),ywallsp(n,4),ywallsp(n,5)
    end do   

    open(unit=11,file='wallsz_'//cmyid//'.txt',position='rewind')
    write(11,'(3A)') 'k, il, iu+1, jl, ju+1'
    do n=1,nzwallsshear
      write(11,'(5I)') zwallsshear(n,1),zwallsshear(n,2),zwallsshear(n,3),zwallsshear(n,4),zwallsshear(n,5)
    end do  

 
    end subroutine createwalls

    subroutine ibmshear
    use modglobal, only : libm
    use modfields, only : up
    
    if (libm == .true.) then
 ! compute shear forces due to IBM
        call xwallshear
        call ywallshearplus   ! due to parallellisation differentiation between + and - side
        call ywallshearmin    ! due to parallellisation differentiation between + and - side
        call zwallshear

    end if
 
    end subroutine ibmshear

    subroutine xwallshear
       use modglobal, only : dzf,dzhiq,dzhi,dxf,dxfi,dxhi,dyi,lles,nsv,lwallfunc,numol,ltempeq,lmoist,&
                             ih,jh,kh,ihc,jhc,khc,dxh,dy,dt,totavtime,rk3step,ib,ie,kb,ke
       use modfields, only : um,up,v0,w0,vp,wp,shear,thl0,thlp,qt0,qtp,sv0,svp
       use modboundary, only    : wallaw
       use modsubgriddata, only : ekm,loneeqn
       use modsurfdata,    only : wtsurf

       real emmo, epmo, epom, emom
       integer i,j,k,n,nc,jl,ju,kl,ku,im,jm,jp,km



     if (lwallfunc==.true.) then
       do n=1,nxwallsshear   ! loop over all shear x-walls
         i  = xwallsshear(n,1)
         im = i-1
         jl = xwallsshear(n,2)    ! starting j-index
         ju = xwallsshear(n,3)    ! ending j-index
         kl = xwallsshear(n,4)    ! starting k-index
         ku = xwallsshear(n,5)    ! ending k-index

! first vp...
         do k=kl,ku-1
           km = k-1
           do j=jl,ju
             jm = j-1
!             write(*,*) "jp",jp
!             write(*,*) "jm",jm
!             write(*,*) "km",km

             epmo = 0.25 * ( ( ekm(im,j,k)+ekm(im,jm,k))*dxf(i) + &
                              (ekm(i,j,k)+ekm(i,jm,k))*dxf(im) ) * dxhi(i)  ! dx is non-equidistant
       
         ! subtract 'standard' diffusion term and add IBM diffusion term
             vp(im,j,k) = vp(im,j,k) +  (-(v0(i,j,k)-v0(im,j,k))*  epmo*dxhi(i) - shear(im,j,k,6)) * dxfi(im)
             vp(i,j,k)  = vp(i,j,k)  +  (-(v0(i,j,k)-v0(im,j,k))* -epmo*dxhi(i) - shear(i,j,k,5))  * dxfi(i)         
         ! Determining wall-shear stresses at immersed boundaries (also used in near-wall damping function in closure routine)
             call wallaw(v0(im,j,k),0.5*(thl0(im,j,k)+thl0(im,jm,k)),wtsurf,dxf(im),numol,shear(im,j,k, 6))
             call wallaw(v0(i, j,k),0.5*(thl0(i ,j,k)+thl0(i ,jm,k)),wtsurf,dxf(i ),numol,shear(i, j,k, 5))
           end do ! jl,ju
         end do ! kl,ku-1


  ! ...then wp
         do k=kl,ku
           km = k-1
           do j=jl,ju-1
             jm = j-1
             
             emom = ( dzf(km) * ( ekm(i,j,k)*dxf(im)  + ekm(im,j,k)*dxf(i) )*dxhi(i)  + &
                      dzf(k)  * ( ekm(i,j,km)*dxf(im) + ekm(im,j,km)*dxf(i) )*dxhi(i) ) * dzhiq(k)
 
         ! Determining wall-shear stresses at immersed boundaries (also used in near-wall damping function in closure routine)
             call wallaw(w0(im,j,k),0.5*(thl0(im,j,k)*dzf(km)+thl0(im,j,km)*dzf(k))*dzhi(k),wtsurf,dxf(im),numol,shear(im,j,k,10))
             call wallaw(w0(i ,j,k),0.5*(thl0(i ,j,k)*dzf(km)+thl0(i ,j,km)*dzf(k))*dzhi(k),wtsurf,dxf(i ),numol,shear(i ,j,k, 9))

         ! subtract 'standard' diffusion term and add IBM diffusion term
             wp(im,j,k) = wp(im,j,k) +  (-(w0(i,j,k)-w0(im,j,k))*  emom*dxhi(i) - shear(im,j,k,10))* dxfi(im)
             wp(i,j,k)  = wp(i,j,k)  +  (-(w0(i,j,k)-w0(im,j,k))* -emom*dxhi(i) - shear(i,j,k,9))  * dxfi(i)
  
           end do ! jl,ju
         end do ! kl,ku-1
         if (ltempeq) call xwallscalar(ih,jh,kh,thl0,thlp,bcthltypex(n),bcthlvaluex(n),n)
         if (lmoist) call xwallscalar(ih,jh,kh,qt0,qtp,bcqttypex(n),bcqtvaluex(n),n)
         do nc=1,nsv
           call xwallscalar(ihc,jhc,khc,sv0(:,:,:,nc),svp(:,:,:,nc),0,0.,n)         ! always use zero flux BC (for now)
         end do
         if (lles == .true. .and. loneeqn == .true.) then
           call xwalle12(n)
         end if
       end do ! n=1,nxwallsshear

     else ! no wall function
       do n=1,nxwallsshear   ! loop over all shear x-walls
         i  = xwallsshear(n,1)
         im = i-1
         jl = xwallsshear(n,2)    ! starting j-index
         ju = xwallsshear(n,3)    ! ending j-index
         kl = xwallsshear(n,4)    ! starting k-index
         ku = xwallsshear(n,5)    ! ending k-index
      
! first vp...
         do k=kl,ku-1
           km = k-1
           do j=jl,ju
             jm = j-1

             epmo = 0.25 * ( ( ekm(im,j,k)+ekm(im,jm,k))*dxf(i) + &
                              (ekm(i,j,k)+ekm(i,jm,k))*dxf(im) ) * dxhi(i)  ! dx is non-equidistant

         ! Determining wall-shear stresses at immersed boundaries (also used in near-wall damping function in closure routine)
!             shear(im,j,k,6)=  epmo*2.*v0(im,j,k)*dxfi(im) ! = -epmo*( v0(i+1,j,k)-v0(i,j,k))/dx^2  + epmo*( 2*v0(i,j,k))/dx^2
!             shear(i,j,k,5) =  emmo*2.*v0(i,j,k)*dxfi(i)   ! = -emmo*(-(v0(i,j,k)-v0(i-1,j,k))/dx^2 - emmo*( 2*v0(i,j,k))/dx^2
             shear(im,j,k,6)=  numol*2.*v0(im,j,k)*dxfi(im) ! = -epmo*( v0(i+1,j,k)-v0(i,j,k))/dx^2  + nu*( 2*v0(i,j,k))/dx^2
             shear(i,j,k,5) =  numol*2.*v0(i,j,k)*dxfi(i)   ! = -emmo*(-(v0(i,j,k)-v0(i-1,j,k))/dx^2 - nu*( 2*v0(i,j,k))/dx^2
         ! subtract 'standard' diffusion term and add IBM diffusion term
             vp(im,j,k) = vp(im,j,k) +  (-(v0(i,j,k)-v0(im,j,k))*  epmo*dxhi(i) - shear(im,j,k,6)) * dxfi(im)
             vp(i,j,k)  = vp(i,j,k)  +  (-(v0(i,j,k)-v0(im,j,k))* -epmo*dxhi(i) - shear(i,j,k,5))  * dxfi(i)
           end do
         end do
         
! ...then wp
         do k=kl,ku
           km = k-1
           do j=jl,ju-1
             jm = j-1

             emom = ( dzf(km) * ( ekm(i,j,k)*dxf(im)  + ekm(im,j,k)*dxf(i) )*dxhi(i)  + &
                      dzf(k)  * ( ekm(i,j,km)*dxf(im) + ekm(im,j,km)*dxf(i) )*dxhi(i) ) * dzhiq(k)

         ! Determining wall-shear stresses at immersed boundaries (also used in near-wall damping function in closure routine)
!             shear(im,j,k,10)= epom*2.*w0(im,j,k)*dxfi(im) ! = -epom*( w0(i+1,j,k)-w0(i,j,k))/dx^2  + epom*( 2*w0(i,j,k))/dx^2
!             shear(i,j,k,9)=   emom*2.*w0(i,j,k)*dxfi(i)   ! = -emom*(-(w0(i,j,k)-w0(i-1,j,k))/dx^2 - emom*( 2*w0(i,j,k))/dx^2
             shear(im,j,k,10)= numol*2.*w0(im,j,k)*dxfi(im) ! = -epom*( w0(i+1,j,k)-w0(i,j,k))/dx^2  + nu*( 2*w0(i,j,k))/dx^2
             shear(i,j,k,9)=   numol*2.*w0(i,j,k)*dxfi(i)   ! = -emom*(-(w0(i,j,k)-w0(i-1,j,k))/dx^2 - nu*( 2*w0(i,j,k))/dx^2

         ! subtract 'standard' diffusion term and add IBM diffusion term
             wp(im,j,k) = wp(im,j,k) +  (-(w0(i,j,k)-w0(im,j,k))*  emom*dxhi(i) - shear(im,j,k,10))* dxfi(im)
             wp(i,j,k)  = wp(i,j,k)  +  (-(w0(i,j,k)-w0(im,j,k))* -emom*dxhi(i) - shear(i,j,k,9))  * dxfi(i)

           end do ! jl,ju-1
         end do ! kl,ku
         if (ltempeq) call xwallscalar(ih,jh,kh,thl0,thlp,bcthltypex(n),bcthlvaluex(n),n)
         if (lmoist)  call xwallscalar(ih,jh,kh,qt0,qtp,bcqttypex(n),bcqtvaluex(n),n)
         do nc=1,nsv
           call xwallscalar(ihc,jhc,khc,sv0(:,:,:,nc),svp(:,:,:,nc),0,0.,n)
         end do
         if (lles == .true. .and. loneeqn == .true.) then
           call xwalle12(n)
         end if
       end do ! n=1,nxwallsshear
     end if

 

    end subroutine xwallshear
    
    subroutine xwallscalar(hi,hj,hk,putin,putout,bctype,bcvalue,n)
      ! nodig: nxwallsshear, xwallsshear(:,:)
       use modglobal, only : dxf,dxfi,dxhi,dxh2i,nsv,ib,ie,jb,je,kb,ke,prandtlmoli,numol
       use modsubgriddata, only : ekh
       
       integer i,j,k,jl,ju,kl,ku

       integer, intent(in) :: hi                                                  !<size of halo in i
       integer, intent(in) :: hj                                                  !<size of halo in j
       integer, intent(in) :: hk                                                  !<size of halo in k
       real, intent(in)    :: putin(ib-hi:ie+hi,jb-hj:je+hj,kb-hk:ke+hk)
       real, intent(inout) :: putout(ib-hi:ie+hi,jb-hj:je+hj,kb:ke+hk)
       real, intent(in)    :: bcvalue
       integer, intent(in)    :: bctype
       integer, intent(in)    :: n

       i  = xwallsshear(n,1)
       jl = xwallsshear(n,2)      ! starting j-index
       ju = xwallsshear(n,3)-1    ! ending j-index one less than shear components
       kl = xwallsshear(n,4)      ! starting k-index one less than shear components
       ku = xwallsshear(n,5)-1    ! ending k-index   

       if (bctype == 1) then    ! fixed value
         do k=kl,ku
           do j=jl,ju
              putout(i,j,k) = putout(i,j,k) + ( &
                          0.5*(ekh(i,j,k)*dxf(i-1)+ekh(i-1,j,k)*dxf(i))* &
                              (putin(i,j,k)-putin(i-1,j,k))*dxh2i(i) -   &
!                              (ekh(i,j,k)*dxf(i-1)+ekh(i-1,j,k)*dxf(i))* &
!                              (putin(i,j,k)-bcvalue)*dxfi(i)*dxhi(i)     &
                              2.*prandtlmoli*numol* &
                              (putin(i,j,k)-bcvalue)*dxfi(i)     &
                                    )*dxfi(i)               ! fixed value

              putout(i-1,j,k) = putout(i-1,j,k) + ( &
                         - 0.5*(ekh(i,j,k)*dxf(i-1)+ekh(i-1,j,k)*dxf(i)) *  &
                               (putin(i,j,k)-putin(i-1,j,k))*dxh2i(i) +     &
!                               (ekh(i,j,k)*dxf(i-1)+ekh(i-1,j,k)*dxf(i)) *  &
!                               (bcvalue-putin(i-1,j,k))*dxfi(i-1)*dxhi(i)   &
                               2.*prandtlmoli*numol*  &
                               (bcvalue-putin(i-1,j,k))*dxfi(i-1)   &
                                    )*dxfi(i-1)            ! fixed value
           end do
         end do
       else    ! bctype == 0 -> fixed flux
         do k=kl,ku
           do j=jl,ju
              putout(i,j,k) = putout(i,j,k) + ( &
                          0.5*(ekh(i,j,k)*dxf(i-1)+ekh(i-1,j,k)*dxf(i))* &
                              (putin(i,j,k)-putin(i-1,j,k))*dxh2i(i) -   &
                               bcvalue )*dxfi(i)               ! fixed flux

              putout(i-1,j,k) = putout(i-1,j,k) + ( &
                         - 0.5*(ekh(i,j,k)*dxf(i-1)+ekh(i-1,j,k)*dxf(i)) * &
                               (putin(i,j,k)-putin(i-1,j,k))*dxh2i(i) +      &
                                bcvalue )*dxfi(i-1)            ! fixed flux
           end do
         end do
       end if 



    end subroutine xwallscalar


    subroutine xwalle12(n)
      ! nodig: nxwallsshear, xwallsshear(:,:)
       use modglobal, only : dxf,dxfi,dxhi,dxh2i,dyi,dzhi,nsv,ib,ie,ih,jb,je,jh,kb,&
                             ke,kh,numol
       use modfields, only : e12p,e120,u0,v0,w0
       use modsubgriddata, only : ekm
       integer i,j,k,im,jl,ju,kl,ku,jm,jp,km,kp

       integer, intent(in)    :: n

       i  = xwallsshear(n,1)
       im = i-1
       jl = xwallsshear(n,2)      ! starting j-index
       ju = xwallsshear(n,3)-1    ! ending j-index one less than shear components
       kl = xwallsshear(n,4)      ! starting k-index one less than shear components
       ku = xwallsshear(n,5)-1    ! ending k-index

       do k=kl,ku
         kp = k+1
         km = k-1
         do j=jl,ju
            jp = j+1
            jm = j-1
            ! walls
            e12p(i,j,k) = e12p(i,j,k) + &
                            (ekm(i,j,k)*dxf(im)+ekm(im,j,k)*dxf(i))* &
                            (e120(i,j,k)-e120(im,j,k))*dxh2i(i)  *0.5*dxfi(i)             ! zero flux

            e12p(im,j,k) = e12p(im,j,k) + &
                             - (ekm(i,j,k)*dxf(im)+ekm(im,j,k)*dxf(i)) * &
                             (e120(i,j,k)-e120(im,j,k))*dxh2i(i)*0.5*dxfi(im)            ! zero flux
            ! sources
            
            e12p(i,j,k) = e12p(i,j,k) + 0.25*(ekm(i,j,k)-numol)/(2*e120(i,j,k))*(& 
             - ((w0(i,j,kp)-w0(im,j,kp))  * dxhi(i)      + &     ! subtract old term
               (u0(i,j,kp)-u0(i,j,k))     * dzhi(kp))**2 + &               
               (( 2.*w0(i,j,kp)           * dxfi(i))     + &     ! IBM for ux0 wall 
               (u0(i,j,kp)-u0(i,j,k))     * dzhi(kp))**2 + &
             
             - ((w0(i,j,k)-w0(im,j,k))    * dxhi(i)      + &     ! subtract old term
               (u0(i,j,k)-u0(i,j,km))     * dzhi(k)) **2 + &
               (( 2.*w0(i,j,k)            * dxfi(i))     + &     ! IBM for ux0 wall
               (u0(i,j,k)-u0(i,j,km))     * dzhi(k)) **2 + &

             - ((u0(i,jp,k)-u0(i,j,k))    * dyi          + &     ! subtract old term
               (v0(i,jp,k)-v0(im,jp,k))   * dxhi(i))**2  + &
               ((u0(i,jp,k)-u0(i,j,k))    * dyi          + &
               (2.*v0(i,jp,k)             * dxfi(i)))**2 + &     ! IBM for ux0 wall

             - ((u0(i,j,k)-u0(i,jm,k))    * dyi          + &     ! subtract old term
               (v0(i,j,k)-v0(im,j,k))     * dxhi(i))**2  + &
               ((u0(i,j,k)-u0(i,jm,k))    * dyi          + &
               (2.*v0(i,j,k)              * dxfi(i)))**2   )     ! IBM for ux0 wall

   
            e12p(im,j,k) = e12p(im,j,k) +0.25*(ekm(im,j,k)-numol)/(2*e120(im,j,k))*(&
             - ((w0(i,j,k)-w0(im,j,k))    * dxhi(i)      + &     ! subtract old term
               (u0(i,j,k)-u0(i,j,km))     * dzhi(k))**2  + &
               ((w0(i,j,k)-w0(im,j,k))    * dxhi(i)      + &     
               (-2.*w0(im,j,k)            * dxfi(im)))**2+ &     ! IBM for ux0 wall

             - ((w0(i,j,kp)-w0(im,j,kp))  * dxhi(i)      + &     ! subtract old term
               (u0(i,j,kp)-u0(i,j,k))     * dzhi(kp))**2 + &
               ((w0(i,j,kp)-w0(im,j,kp))  * dxhi(i)      + &
               (-2.*w0(im,j,kp)           * dxfi(im)))**2+ &     ! IBM for ux0 wall

             - ((u0(i,j,k)-u0(i,jm,k))    * dyi          + &     ! subtract old term
               (v0(i,j,k)-v0(im,j,k))     * dxhi(i) )**2 + &
               ((u0(i,j,k)-u0(i,jm,k))    * dyi          + &
               (-2.*v0(im,j,k)            * dxfi(im)))**2+ &     ! IBM for ux0 wall

             - ((u0(i,jp,k)-u0(i,j,k))    * dyi          + &     ! subtract old term
               (v0(i,jp,k)-v0(im,jp,k))   * dxhi(i))**2  + &
               ((u0(i,jp,k)-u0(i,j,k))    * dyi          + & 
               (-2.*v0(im,jp,k)           * dxfi(im)))**2  )     ! IBM for ux0 wall    
  
         end do
       end do
    end subroutine xwalle12


    subroutine ywallshearplus
      ! nodig: nywallsp, ywallsp(:,:)
       use modglobal, only : dzf,dzhiq,dzhi,dxf,dxhi,dy,dyi,nsv,lles,lwallfunc,numol,ltempeq,lmoist,&
                             ih,jh,kh,ihc,jhc,khc
       use modfields, only : u0,w0,up,wp,shear,thlp,thl0,qtp,qt0,sv0,svp
       use modboundary, only : wallaw
       use modsubgriddata, only : ekm,loneeqn
       use modsurfdata,     only : wtsurf
       real emmo, eomm 
       integer i,j,k,n,nc,il,iu,kl,ku,im,jm,km

     if (lwallfunc==.true.) then
       do n=1,nywallsp      ! loop over all shearplus y-walls
         j  = ywallsp(n,1)
         jm = j-1
         il = ywallsp(n,2)    ! starting i-index
         iu = ywallsp(n,3)    ! ending i-index
         kl = ywallsp(n,4)    ! starting k-index
         ku = ywallsp(n,5)    ! ending k-index
  ! first up ...
         do k=kl,ku-1
           km = k-1
           do i=il,iu
           im = i-1
           emmo = 0.25 * ( ( ekm(i,j,k)+ekm(i,jm,k))*dxf(im)  +&
                  (ekm(im,jm,k)+ekm(im,j,k))*dxf(i) ) * dxhi(i)  ! dx is non-equidistant
  ! Determining wall-shear stresses at immersed boundaries (also used in near-wall damping function in closure routine)
           call wallaw(u0(i,j,k),0.5*(thl0(i ,j,k)*dxf(i) + thl0(im,j,k)*dxf(im))*dxhi(i),wtsurf,dy,numol,shear(i,j,k,1))
  ! subtract 'standard' diffusion term and add IBM diffusion term
           up(i,j,k) = up(i,j,k) + (-(u0(i,j,k)-u0(i,jm,k))*-emmo*dyi - shear(i,j,k,1)) *dyi 
           end do ! il,iu
         end do ! kl,ku-1

  ! ...then wp
         do k=kl,ku
           km = k-1
           do i=il,iu-1
           im = i-1
           eomm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i,jm,km) ) ) * dzhiq(k)  ! dz is non-eqidistant

  ! Determining wall-shear stresses at immersed boundaries (also used in near-wall damping function in closure routine)
           call wallaw(w0(i,j,k),0.5*(thl0(i ,j,k)*dzf(k) + thl0(i,j,km)*dzf(km))*dzhi(k),wtsurf,dy,numol,shear(i,j,k,11))
  ! subtract 'standard' diffusion term and add IBM diffusion term
           wp(i,j,k) = wp(i,j,k) + (-(w0(i,j,k)-w0(i,jm,k))*-eomm*dyi - shear(i,j,k,11))*dyi 
           end do ! il,iu-1
         end do ! kl,ku
         if (ltempeq) call ywallscalarplus(ih,jh,kh,thl0,thlp,bcthltypeyp(n),bcthlvalueyp(n),n)
         if (lmoist) call ywallscalarplus(ih,jh,kh,qt0,qtp,bcqttypeyp(n),bcqtvalueyp(n),n)
         do nc=1,nsv
           call ywallscalarplus(ihc,jhc,khc,sv0(:,:,:,nc),svp(:,:,:,nc),0,0.,n)
         end do
         if (lles == .true. .and. loneeqn == .true.) then
           call ywalle12plus(n)
         end if
       end do ! 1,nywallsp

     else  ! no wall function

       do n=1,nywallsp      ! loop over all shearplus y-walls
         j  = ywallsp(n,1)
         jm = j-1
         il = ywallsp(n,2)    ! starting i-index
         iu = ywallsp(n,3)    ! ending i-index
         kl = ywallsp(n,4)    ! starting k-index
         ku = ywallsp(n,5)    ! ending k-index
  ! first up...

         do k=kl,ku-1
          km = k-1
         do i=il,iu
           im = i-1
           emmo = 0.25 * ( ( ekm(i,j,k)+ekm(i,jm,k))*dxf(im)  +&
                  (ekm(im,jm,k)+ekm(im,j,k))*dxf(i) ) * dxhi(i)  ! dx is non-equidistant
  ! Determining wall-shear stresses at immersed boundaries (also used in near-wall damping function in closure routine)
!           shear(i,j,k,1)= emmo*dyi*2.*u0(i,j,k)      ! new term: -(emmo/dy)*2*u0(i,j,k)
           shear(i,j,k,1)= numol*dyi*2.*u0(i,j,k)      ! new term: -(nu/dy)*2*u0(i,j,k)
  ! subtract 'standard' diffusion term and add IBM diffusion term
           up(i,j,k) = up(i,j,k) + (-(u0(i,j,k)-u0(i,jm,k))*-emmo*dyi - shear(i,j,k,1)) *dyi
           end do ! il,iu
         end do ! kl,ku-1
  ! ...then wp
         do k=kl,ku
          km = k-1
         do i=il,iu-1
           im = i-1
           eomm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                       dzf(k)  * ( ekm(i,j,km) + ekm(i,jm,km) ) ) * dzhiq(k)
  ! Determining wall-shear stresses at immersed boundaries (also used in near-wall damping function in closure routine)
!           shear(i,j,k,11)=  eomm*2.*w0(i,j,k)*dyi
           shear(i,j,k,11)=  numol*2.*w0(i,j,k)*dyi
  ! subtract 'standard' diffusion term and add IBM diffusion term
           wp(i,j,k) = wp(i,j,k) + (-(w0(i,j,k)-w0(i,jm,k))*-eomm*dyi - shear(i,j,k,11))*dyi
           end do ! il,iu-1
         end do ! kl,ku
         if (ltempeq) call ywallscalarplus(ih,jh,kh,thl0,thlp,bcthltypeyp(n),bcthlvalueyp(n),n)
         if (lmoist) call ywallscalarplus(ih,jh,kh,qt0,qtp,bcqttypeyp(n),bcqtvalueyp(n),n)
         do nc=1,nsv
           call ywallscalarplus(ihc,jhc,khc,sv0(:,:,:,nc),svp(:,:,:,nc),0,0.,n)
         end do
         if (lles == .true. .and. loneeqn == .true.) then
           call ywalle12plus(n)
         end if
       end do ! 1,nywallsp
    
     end if ! lwallfunc 
    
    end subroutine ywallshearplus
  
    subroutine ywallscalarplus(hi,hj,hk,putin,putout,bctype,bcvalue,n)
      ! nodig: nxwallsshear, xwallsshear(:,:)
       use modglobal, only : dyi,ib,ie,jb,je,kb,ke,numol,prandtlmoli
       use modsubgriddata, only : ekh
       integer i,j,k,il,iu,kl,ku

       integer, intent(in) :: hi                                                  !<size of halo in i
       integer, intent(in) :: hj                                                  !<size of halo in j
       integer, intent(in) :: hk                                                  !<size of halo in k
       real, intent(in)    :: putin (ib-hi:ie+hi,jb-hj:je+hj,kb-hk:ke+hk)
       real, intent(inout) :: putout(ib-hi:ie+hi,jb-hj:je+hj,kb   :ke+hk)
       real, intent(in)    :: bcvalue
       integer, intent(in)    :: bctype
       integer, intent(in)    :: n

       j  = ywallsp(n,1)
       il = ywallsp(n,2)      ! starting i-index
       iu = ywallsp(n,3)-1    ! ending i-index one less than shear components
       kl = ywallsp(n,4)      ! starting k-index one less than shear components
       ku = ywallsp(n,5)-1    ! ending k-index

      
       if (bctype == 1) then 
         do k=kl,ku
           do i=il,iu
!              putout(i,j-1,k) = putout(i,j-1,k) + ( &
!                          - 0.5*(ekh(i,j,k)+ekh(i,j-1,k))*(putin(i,j,k)-putin(i,j-1,k))*dyi  &    
!                              + (ekh(i,j,k)+ekh(i,j-1,k))*(bcvalue - putin(i,j-1,k))*dyi     &      ! fixed value
!                                                   )*dyi
              putout(i,j,k) = putout(i,j,k) + ( &
                            0.5*(ekh(i,j,k)+ekh(i,j-1,k))*(putin(i,j,k)-putin(i,j-1,k))*dyi  &    
!                              - (ekh(i,j,k)+ekh(i,j-1,k))*(putin(i,j,k)-bcvalue)*dyi     &      ! fixed value
                              - 2.*prandtlmoli*numol*(putin(i,j,k)-bcvalue)*dyi     &      ! fixed value
                                                   )*dyi
           end do 
         end do
       else  ! bctype == 0
         do k=kl,ku
           do i=il,iu
!              putout(i,j-1,k) = putout(i,j-1,k) + ( &
!                          - 0.5*(ekh(i,j,k)+ekh(i,j-1,k))*(putin(i,j,k)-putin(i,j-1,k))*dyi  &    
!                              +  bcvalue     &                                                      ! fixed flux
!                                                   )*dyi
              putout(i,j,k) = putout(i,j,k) + ( &
                           0.5*(ekh(i,j,k)+ekh(i,j-1,k))*(putin(i,j,k)-putin(i,j-1,k))*dyi  &    
                              -  bcvalue     &                                                      ! fixed flux
                                                   )*dyi
           end do
         end do
       end if
    end subroutine ywallscalarplus

    subroutine ywalle12plus(n)
      ! nodig: nxwallsshear, xwallsshear(:,:)
       use modglobal, only : dxhi,dyi,dy2i,dzhi,ib,ie,ih,jb,je,jh,kb,ke,kh,numol
       use modfields, only : e12p,e120,u0,v0,w0
       use modsubgriddata, only : ekm
       integer i,j,k,il,iu,kl,ku,im,ip,jm,km,kp

       integer, intent(in)    :: n

       j  = ywallsp(n,1)
       jm = j-1
       il = ywallsp(n,2)      ! starting i-index
       iu = ywallsp(n,3)-1    ! ending i-index one less than shear components
       kl = ywallsp(n,4)      ! starting k-index one less than shear components
       ku = ywallsp(n,5)-1    ! ending k-index


       do k=kl,ku
         kp = k+1
         km = k-1
         do i=il,iu
            im = i-1
            ip = i+1
            ! walls
            e12p(i,j,k) = e12p(i,j,k) + &
                            (ekm(i,j,k)+ekm(i,jm,k))*(e120(i,j,k)-e120(i,jm,k))*0.5*dy2i      ! zero flux
 
            e12p(i,j,k) = e12p(i,j,k) + 0.25*(ekm(i,j,k)-numol)/(2*e120(i,j,k))*(&
             - ((u0(i,j,k)-u0(i,jm,k))   * dyi          + &
               (v0(i,j,k)-v0(im,j,k))    * dxhi(i))**2  + &
               (2.*u0(i,j,k)             * dyi          + &
               (v0(i,j,k)-v0(im,j,k))    * dxhi(i))**2  + &

             - ((u0(ip,j,k)-u0(ip,jm,k)) * dyi          + &
               (v0(ip,j,k)-v0(i,j,k))    * dxhi(ip))**2 + &
               (2.*u0(ip,j,k)            * dyi          + &       ! IBM for vy0 wall (uyp-case)
               (v0(ip,j,k)-v0(i,j,k))    * dxhi(ip))**2 + &
             
             - ((v0(i,j,kp)-v0(i,j,k))   * dzhi(kp)     + &
               (w0(i,j,kp)-w0(i,jm,kp))  * dyi     )**2 + &
               ((v0(i,j,kp)-v0(i,j,k))   * dzhi(kp)     + &
               (2.*w0(i,j,kp)            * dyi)    )**2 + &

             - ((v0(i,j,k)-v0(i,j,km))   * dzhi(k)      + &
               (w0(i,j,k)-w0(i,jm,k))    * dyi     )**2 + &
               ((v0(i,j,k)-v0(i,j,km))   * dzhi(k)      + &
               (2.*w0(i,j,k)             * dyi)    )**2   )

         end do
       end do
    end subroutine ywalle12plus




    subroutine ywallshearmin
       use modglobal, only : dxf,dxhi,dy,dyi,dzhiq,dzf,dzhi,lles,nsv,lwallfunc,numol,ltempeq,lmoist,&
                             ih,jh,kh,ihc,jhc,khc
       use modfields, only : u0,w0,up,wp,shear,thl0,thlp,qt0,qtp,sv0,svp
       use modboundary, only : wallaw
       use modsubgriddata, only : ekm,loneeqn
       use modsurfdata,     only : wtsurf
       real empo,eopm
       integer i,j,k,n,nc,il,iu,kl,ku,im,jp,km

     if (lwallfunc==.true.) then
       do n=1,nywallsm      ! loop over all shearplus y-walls
!         j  = ywallsm(n,1)-1  ! this is 'min' side
         j  = ywallsm(n,1)  ! this is 'min' side
         jp = j+1
         il = ywallsm(n,2)    ! starting i-index
         iu = ywallsm(n,3)    ! ending i-index
         kl = ywallsm(n,4)    ! starting k-index
         ku = ywallsm(n,5)    ! ending k-index

  ! first up...
         do k=kl,ku-1
         km = k-1
         do i=il,iu
           im = i-1
           empo = 0.25 * ( ( ekm(i,j,k)+ekm(i,jp,k))*dxf(im) + &
                       (ekm(im,j,k)+ekm(im,jp,k))*dxf(i) ) * dxhi(i)  ! dx is non-equidistant
  ! Determining wall-shear stresses at immersed boundaries (also used in near-wall damping function in closure routine)
          call wallaw(u0(i,j,k),0.5*(thl0(i ,j,k)*dxf(i) + thl0(im,j,k)*dxf(im))*dxhi(i),wtsurf,dy,numol,shear(i,j,k,2))
  ! subtract 'standard' diffusion term and add IBM diffusion term
!           up(i,j,k) = up(i,j,k) + (-(u0(i,jp,k)-u0(i,j,k))*empo*dyi + shear(i,j,k,2)) *dyi
           up(i,j,k) = up(i,j,k) + (-(u0(i,jp,k)-u0(i,j,k))*empo*dyi - shear(i,j,k,2)) *dyi
           end do ! il,iu
         end do ! kl,ku-1

  ! ...then wp
         do k=kl,ku
         km = k-1
         do i=il,iu-1
           im = i-1
           eopm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jp,k)  )  + &
                       dzf(k)  * ( ekm(i,j,km) + ekm(i,jp,km) ) ) * dzhiq(k)

  ! Determining wall-shear stresses at immersed boundaries (also used in near-wall damping function in closure routine)
          call wallaw(w0(i,j,k),0.5*(thl0(i ,j,k)*dzf(k) + thl0(i,j,km)*dzf(km))*dzhi(k),wtsurf,dy,numol,shear(i,j,k,12))
  ! subtract 'standard' diffusion term and add IBM diffusion term
           wp(i,j,k) = wp(i,j,k) + (-(w0(i,jp,k)-w0(i,j,k))*eopm*dyi - shear(i,j,k,12)) *dyi
!           wp(i,j,k) = wp(i,j,k) + (-(w0(i,jp,k)-w0(i,j,k))*eopm*dyi + shear(i,j,k,12)) *dyi
           end do ! il,iu-1
         end do ! kl,ku

         if (ltempeq) call ywallscalarmin(ih,jh,kh,thl0,thlp,bcthltypeym(n),bcthlvalueym(n),n)
         if (lmoist)  call ywallscalarmin(ih,jh,kh,qt0,qtp,bcqttypeym(n),bcqtvalueym(n),n) 
         do nc=1,nsv
           call ywallscalarmin(ihc,jhc,khc,sv0(:,:,:,nc),svp(:,:,:,nc),0,0.,n)
         end do
         if (lles == .true. .and. loneeqn == .true.) then
           call ywalle12min(n)
         end if
       end do ! 1,nywallsp

     else ! no wall function

       do n=1,nywallsm      ! loop over all shearplus y-walls
!         j  = ywallsm(n,1)-1  ! this is 'min' side
         j  = ywallsm(n,1)  ! this is 'min' side
         jp = j+1
         il = ywallsm(n,2)    ! starting i-index
         iu = ywallsm(n,3)    ! ending i-index
         kl = ywallsm(n,4)    ! starting k-index
         ku = ywallsm(n,5)    ! ending k-index
  ! first up...
         do k=kl,ku-1
         km = k-1
         do i=il,iu
           im = i-1
           empo = 0.25 * ( ( ekm(i,j,k)+ekm(i,jp,k))*dxf(im) + &
                       (ekm(im,j,k)+ekm(im,jp,k))*dxf(i) ) * dxhi(i)  ! dx is non-equidistant
  ! Determining wall-shear stresses at immersed boundaries (also used in near-wall damping function in closure routine)
!           shear(i,j,k,2)  =  empo*dyi*2.*u0(i,j,k)      ! new term: -(empo/dy)*2*u0(i,j,k)
           shear(i,j,k,2)  =  numol*dyi*2.*u0(i,j,k)      ! new term: -(nu/dy)*2*u0(i,j,k)
  ! subtract 'standard' diffusion term and add IBM diffusion term
!           up(i,j,k) = up(i,j,k) + (-(u0(i,jp,k)-u0(i,j,k))*empo*dyi + shear(i,j,k,2)) *dyi
           up(i,j,k) = up(i,j,k) + (-(u0(i,jp,k)-u0(i,j,k))*empo*dyi - shear(i,j,k,2)) *dyi
           end do ! kl,ku-1
         end do ! il,iu
  ! ...then wp
         do k=kl,ku
         km = k-1
         do i=il,iu-1
           im = i-1
           eopm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jp,k)  )  + &
                       dzf(k)  * ( ekm(i,j,km) + ekm(i,jp,km) ) ) * dzhiq(k)
  ! Determining wall-shear stresses at immersed boundaries (also used in near-wall damping function in closure routine)
!           shear(i,j,k,12) =  eopm*2.*w0(i,j,k)*dyi
           shear(i,j,k,12) =  numol*2.*w0(i,j,k)*dyi
  ! subtract 'standard' diffusion term and add IBM diffusion term
           wp(i,j,k) = wp(i,j,k) + (-(w0(i,jp,k)-w0(i,j,k))*eopm*dyi - shear(i,j,k,12)) *dyi
!           wp(i,j,k) = wp(i,j,k) + (-(w0(i,jp,k)-w0(i,j,k))*eopm*dyi + shear(i,j,k,12)) *dyi
           end do ! il,iu-1
         end do ! kl,ku

         if (ltempeq)  call ywallscalarmin(ih,jh,kh,thl0,thlp,bcthltypeym(n),bcthlvalueym(n),n)
         if (lmoist)   call ywallscalarmin(ih,jh,kh,qt0,qtp,bcqttypeym(n),bcqtvalueym(n),n)
         
         do nc=1,nsv
           call ywallscalarmin(ihc,jhc,khc,sv0(:,:,:,nc),svp(:,:,:,nc),0,0.,n)
         end do
         if (lles == .true. .and. loneeqn == .true.) then
           call ywalle12min(n)
         end if
       end do ! 1,nywallsp

     end if ! lwallfunc

    end subroutine ywallshearmin

    subroutine ywallscalarmin(hi,hj,hk,putin,putout,bctype,bcvalue,n)
      ! nodig: nxwallsshear, xwallsshear(:,:)
       use modglobal, only : dyi,ib,ie,jb,je,kb,ke,prandtlmoli,numol
       use modsubgriddata, only : ekh
       integer i,j,k,il,iu,kl,ku

       integer, intent(in) :: hi                                                  !<size of halo in i
       integer, intent(in) :: hj                                                  !<size of halo in j
       integer, intent(in) :: hk                                                  !<size of halo in k
       real, intent(in)    :: putin (ib-hi:ie+hi,jb-hj:je+hj,kb-hk:ke+hk)
       real, intent(inout) :: putout(ib-hi:ie+hi,jb-hj:je+hj,kb   :ke+hk)
       real, intent(in)    :: bcvalue
       integer, intent(in)    :: bctype
       integer, intent(in)    :: n

       j  = ywallsm(n,1)
       il = ywallsm(n,2)      ! starting i-index
       iu = ywallsm(n,3)-1    ! ending i-index one less than shear components
       kl = ywallsm(n,4)      ! starting k-index one less than shear components
       ku = ywallsm(n,5)-1    ! ending k-index

       if (bctype == 1) then
         do k=kl,ku
           do i=il,iu
!             putout(i,j,k) = putout(i,j,k) + ( &
!                        0.5*(ekh(i,j,k)+ekh(i,j-1,k))*(putin(i,j,k)-putin(i,j-1,k))*dyi &     
!                           -(ekh(i,j,k)+ekh(i,j-1,k))*(putin(i,j,k)-bcvalue)*dyi        &      ! fixed value
!                                             )*dyi        
             putout(i,j,k) = putout(i,j,k) + ( &
                           - 0.5*(ekh(i,j,k)+ekh(i,j+1,k))*(putin(i,j+1,k)-putin(i,j,k))*dyi &     
!                           +(ekh(i,j,k)+ekh(i,j+1,k))*(bcvalue-putin(i,j,k))*dyi        &      ! fixed value
                           +2.*prandtlmoli*numol*(bcvalue-putin(i,j,k))*dyi        &      ! fixed value
                                             )*dyi        
           end do
         end do
       else ! bctype == 0
         do k=kl,ku
           do i=il,iu
!             putout(i,j,k) = putout(i,j,k) + ( &
!                        0.5*(ekh(i,j,k)+ekh(i,j-1,k))*(putin(i,j,k)-putin(i,j-1,k))*dyi &     
!                           - bcvalue           &                                               ! fixed flux
!                                             )*dyi        
             putout(i,j,k) = putout(i,j,k) + ( &
                           - 0.5*(ekh(i,j,k)+ekh(i,j+1,k))*(putin(i,j+1,k)-putin(i,j,k))*dyi &     
                           + bcvalue           &                                               ! fixed flux
                                             )*dyi        
           end do
         end do 
       end if
 
    end subroutine ywallscalarmin

!    subroutine ywalle12min(n)
!      ! nodig: nxwallsshear, xwallsshear(:,:)
!       use modglobal, only : dxhi,dyi,dy2i,dzhi,ib,ie,ih,jb,je,jh,kb,ke,kh,numol
!       use modfields, only : e12p,e120,u0,v0,w0
!       use modsubgriddata, only : ekm
!       integer i,j,k,il,iu,kl,ku,im,ip,jm,km,kp
!
!       integer, intent(in)    :: n
!
!       j  = ywallsm(n,1)
!       jm = j-1
!       il = ywallsm(n,2)      ! starting i-index
!       iu = ywallsm(n,3)-1    ! ending i-index one less than shear components
!       kl = ywallsm(n,4)      ! starting k-index one less than shear components
!       ku = ywallsm(n,5)-1    ! ending k-index
!
!
!       do k=kl,ku
!         kp = k+1
!         km = k-1
!         do i=il,iu
!            im = i-1
!            ip = i+1
!            ! walls
!            e12p(i,j-1,k) = e12p(i,jm,k) + &
!                            (ekm(i,j,k)+ekm(i,jm,k))*(e120(i,j,k)-e120(i,jm,k))*0.5*dy2i      ! zero flux
!
!            e12p(i,jm,k) = e12p(i,jm,k) + 0.25*(ekm(i,jm,k)-numol)/(2*e120(i,jm,k))*(&
!             - ((u0(i,j,k)-u0(i,jm,k))   * dyi          + &
!               (v0(i,j,k)-v0(im,j,k))    * dxhi(i) )**2 + &
!               (-2.*u0(i,jm,k)           * dyi          + &
!               (v0(i,j,k)-v0(im,j,k))    * dxhi(i) )**2 + &
!
!             - ((u0(ip,j,k)-u0(ip,jm,k)) * dyi          + &
!           end do
!         end do 
!       end if
! 
!    end subroutine ywalle12min(n)

    subroutine ywalle12min(n)
      ! nodig: nxwallsshear, xwallsshear(:,:)
       use modglobal, only : dxhi,dyi,dy2i,dzhi,ib,ie,ih,jb,je,jh,kb,ke,kh,numol
       use modfields, only : e12p,e120,u0,v0,w0
       use modsubgriddata, only : ekm
       integer i,j,k,il,iu,kl,ku,im,ip,jm,km,kp

       integer, intent(in)    :: n

       j  = ywallsm(n,1)
       jm = j-1
       il = ywallsm(n,2)      ! starting i-index
       iu = ywallsm(n,3)-1    ! ending i-index one less than shear components
       kl = ywallsm(n,4)      ! starting k-index one less than shear components
       ku = ywallsm(n,5)-1    ! ending k-index


       do k=kl,ku
         kp = k+1
         km = k-1
         do i=il,iu
            im = i-1
            ip = i+1
            ! walls
            e12p(i,j-1,k) = e12p(i,jm,k) + &
                            (ekm(i,j,k)+ekm(i,jm,k))*(e120(i,j,k)-e120(i,jm,k))*0.5*dy2i      ! zero flux

            e12p(i,jm,k) = e12p(i,jm,k) + 0.25*(ekm(i,jm,k)-numol)/(2*e120(i,jm,k))*(&
             - ((u0(i,j,k)-u0(i,jm,k))   * dyi          + &
               (v0(i,j,k)-v0(im,j,k))    * dxhi(i) )**2 + &
               (-2.*u0(i,jm,k)           * dyi          + &
               (v0(i,j,k)-v0(im,j,k))    * dxhi(i) )**2 + &

             - ((u0(ip,j,k)-u0(ip,jm,k)) * dyi          + &
               (v0(ip,j,k)-v0(i,j,k))    * dxhi(ip))**2 + &
               (-2.*u0(ip,jm,k)          * dyi          + &
               (v0(ip,j,k)-v0(i,j,k))    * dxhi(ip))**2 + &

             - ((v0(i,j,k)-v0(i,j,km))   * dzhi(k)      + &
               (w0(i,j,k)-w0(i,jm,k))    * dyi     )**2 + &
               ((v0(i,j,k)-v0(i,j,km))   * dzhi(k)      + &
               (-2.*w0(i,j,k)            * dyi)    )**2 + &

             - ((v0(i,j,kp)-v0(i,j,k))   * dzhi(kp)     + &
               (w0(i,j,kp)-w0(i,jm,kp))  * dyi     )**2 + &
               ((v0(i,j,kp)-v0(i,j,k))   * dzhi(kp)     + &
               (-2.*w0(i,j,kp)           * dyi)    )**2   )

         end do
       end do
    end subroutine ywalle12min


 
    subroutine zwallshear
       use modglobal, only : dzf,dzfi,dzhi,dzhiq,dxf,dxfi,dxhi,dyi,nsv,lles,lwallfunc,numol,ltempeq,lmoist,&
                             ih,jh,kh,ihc,jhc,khc
       use modfields, only : u0,v0,up,vp,shear,thl0,thlp,qt0,qtp,sv0,svp
       use modboundary, only : wallaw
       use modsubgriddata, only : ekm,loneeqn
       use modsurfdata,     only : wtsurf
       real emop,emom,eomp,eomm
       integer i,j,k,n,nc,il,iu,jl,ju,im,jm,km

     if (lwallfunc==.true.) then
       do n=1,nzwallsshear   ! loop over all shear z-walls
         k  = zwallsshear(n,1)
         km = k-1
         il = zwallsshear(n,2)    ! starting i-index
         iu = zwallsshear(n,3)    ! ending i-index
         jl = zwallsshear(n,4)    ! starting j-index
         ju = zwallsshear(n,5)    ! ending j-index

  ! first up...
         do j=jl,ju-1
           jm = j-1
           do i=il,iu
             im = i-1
             emop = ( dzf(k)  * ( ekm(i,j,km)*dxf(im)  + ekm(im,j,km)*dxf(i) )  + &              ! dx is non-equidistant
                      dzf(km) * ( ekm(i,j,k)*dxf(im) + ekm(im,j,k)*dxf(i) ) )*dxhi(i) * dzhiq(k)
              
             emom = ( dzf(km) * ( ekm(i,j,k)*dxf(im)  + ekm(im,j,k)*dxf(i) )  + &             ! dx is non-equidistant
                      dzf(k)  * ( ekm(i,j,km)*dxf(im) + ekm(im,j,km)*dxf(i) ) )*dxhi(i) * dzhiq(k)
   ! Determining wall-shear stresses at immersed boundaries (also used in near-wall damping function in closure routine)
             call wallaw(u0(i,j,km),0.5*(thl0(i,j,km)*dxf(im)+thl0(im,j,km)*dxf(i))*dxhi(i),wtsurf,dzf(km),numol,shear(i,j,km,4))
             call wallaw(u0(i,j,k ),0.5*(thl0(i,j,k )*dxf(im)+thl0(im,j,k )*dxf(i))*dxhi(i),wtsurf,dzf(k ),numol,shear(i,j,k, 3))
   ! subtract 'standard' diffusion term and add IBM diffusion term
             up(i,j,km) = up(i,j,km) + (-(u0(i,j,k)-u0(i,j,km))*  emop*dzhi(k) - shear(i,j,km,4))*dzfi(km)
             up(i,j,k)  = up(i,j,k)  + (-(u0(i,j,k)-u0(i,j,km))* -emop*dzhi(k) - shear(i,j,k,3)) *dzfi(k)
           end do ! i=il,iu
         end do ! j=jl,ju-1


  ! then vp...
         do j=jl,ju
           jm = j-1
           do i=il,iu-1
             im = i-1

             eomm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i,jm,km) ) ) * dzhiq(k)

   ! Determining wall-shear stresses at immersed boundaries (also used in near-wall damping function in closure routine)
             call wallaw(v0(i,j,km),0.5*(thl0(i,j,km)+thl0(i,jm,km)),wtsurf,dzf(km),numol,shear(i,j,km,8))
             call wallaw(v0(i,j,k ),0.5*(thl0(i,j,k )+thl0(i,jm,k )),wtsurf,dzf(k ),numol,shear(i,j,k, 7))
   ! subtract 'standard' diffusion term and add IBM diffusion term
             vp(i,j,km) = vp(i,j,km) + (-(v0(i,j,k)-v0(i,j,km))*  eomm*dzhi(k) - shear(i,j,km,8)) *dzfi(km)
             vp(i,j,k)  = vp(i,j,k)  + (-(v0(i,j,k)-v0(i,j,km))* -eomm*dzhi(k) - shear(i,j,k,7)) *dxfi(k)
           end do ! i=il,iu-1
         end do ! j=jl,ju
         if (ltempeq) call zwallscalar(ih,jh,kh,thl0,thlp,bcthltypez(n),bcthlvaluez(n),n)
         if (lmoist)  call zwallscalar(ih,jh,kh,qt0,qtp,bcqttypez(n),bcqtvaluez(n),n)
          do nc=1,nsv
           call zwallscalar(ihc,jhc,khc,sv0(:,:,:,nc),svp(:,:,:,nc),0,0.,n)
         end do
         if (lles == .true. .and. loneeqn == .true.) then
           call zwalle12(n)
         end if
       end do ! n=1,nzwallsshear

     else ! no wall function
  
       do n=1,nzwallsshear   ! loop over all shear z-walls
         k  = zwallsshear(n,1)
         km = k-1
         il = zwallsshear(n,2)    ! starting i-index
         iu = zwallsshear(n,3)    ! ending i-index
         jl = zwallsshear(n,4)    ! starting j-index
         ju = zwallsshear(n,5)    ! ending j-index

! first up...
         do j=jl,ju-1
           jm = j-1
           do i=il,iu
             im = i-1
             emop = ( dzf(k)  * ( ekm(i,j,km)*dxf(im)  + ekm(im,j,km)*dxf(i) )  + &              ! dx is non-equidistant
                      dzf(km) * ( ekm(i,j,k)*dxf(im) + ekm(im,j,k)*dxf(i) ) )*dxhi(i) * dzhiq(k)

   ! Determining wall-shear stresses at immersed boundaries (also used in near-wall damping function in closure routine)
!             shear(i,j,km,4) = emop*2.*u0(i,j,km)*dzfi(km)   ! new term: -(eomp/dzf(k)) *2*u0(i,j,k)
!             shear(i,j,k,3)  = emom*2.*u0(i,j,k)*dzfi(k)
             shear(i,j,k,3)  = numol*2.*u0(i,j,k)*dzfi(k)
!             shear(i,j,km,8) = numol*2.*u0(i,j,km)*dzfi(km)      ! new term: -(nu/dzf(k))*2*v0(i,j,k)
             shear(i,j,km,4) = numol*2.*u0(i,j,km)*dzfi(km)      ! new term: -(nu/dzf(k))*2*v0(i,j,k)
   ! subtract 'standard' diffusion term and add IBM diffusion term
             up(i,j,km) = up(i,j,km) + (-(u0(i,j,k)-u0(i,j,km))*  emop*dzhi(k) - shear(i,j,km,4))*dzfi(km)
             up(i,j,k)  = up(i,j,k)  + (-(u0(i,j,k)-u0(i,j,km))* -emop*dzhi(k) - shear(i,j,k,3)) *dzfi(k)
           end do ! i=il,iu
         end do ! j=jl,ju-1


         do j=jl,ju
           jm = j-1
           do i=il,iu-1
             im = i-1

             eomm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i,jm,km) ) ) * dzhiq(k)

   ! Determining wall-shear stresses at immersed boundaries (also used in near-wall damping function in closure routine)
!             shear(i,j,km,8) = eomp*2.*v0(i,j,km)*dzfi(km)      ! new term: -(eomp/dzf(k))*2*v0(i,j,k)
!             shear(i,j,km,4) = emop*2.*u0(i,j,km)*dzfi(km)   ! new term: -(eomp/dzf(k)) *2*u0(i,j,k)
             shear(i,j,k,7)  = numol*2.*v0(i,j,k)*dzfi(k)      ! new term: -(nu/dzf(k))*2*v0(i,j,k)
             shear(i,j,km,8) = numol*2.*v0(i,j,km)*dzfi(km)      ! new term: -(nu/dzf(k))*2*v0(i,j,k)
   ! subtract 'standard' diffusion term and add IBM diffusion term
             vp(i,j,km) = vp(i,j,km) + (-(v0(i,j,k)-v0(i,j,km))*  eomm*dzhi(k) - shear(i,j,km,8)) *dzfi(km)
             vp(i,j,k)  = vp(i,j,k)  + (-(v0(i,j,k)-v0(i,j,km))* -eomm*dzhi(k) - shear(i,j,k,7)) *dxfi(k)
           end do ! i=il,iu-1
         end do ! j=jl,ju
         if (ltempeq) call zwallscalar(ih,jh,kh,thl0,thlp,bcthltypez(n),bcthlvaluez(n),n)
         if (lmoist)  call zwallscalar(ih,jh,kh,qt0,qtp,bcqttypez(n),bcqtvaluez(n),n)
do nc=1,nsv
           call zwallscalar(ihc,jhc,khc,sv0(:,:,:,nc),svp(:,:,:,nc),0,0.,n)
         end do
         if (lles == .true. .and. loneeqn == .true.) then
           call zwalle12(n)
         end if
       end do ! n=1,nzwallsshear
     end if ! wallfunc

    end subroutine zwallshear

    subroutine zwallscalar(hi,hj,hk,putin,putout,bctype,bcvalue,n)
       use modglobal, only : dzf,dzfi,dzhi,dzh2i,ib,ie,jb,je,kb,ke,prandtlmoli,numol
       use modsubgriddata, only : ekh
       integer i,j,k,il,iu,jl,ju,km

       integer, intent(in) :: hi                                                  !<size of halo in i
       integer, intent(in) :: hj                                                  !<size of halo in j
       integer, intent(in) :: hk                                                  !<size of halo in k
       real, intent(in)    :: putin (ib-hi:ie+hi,jb-hj:je+hj,kb-hk:ke+hk)
       real, intent(inout) :: putout(ib-hi:ie+hi,jb-hj:je+hj,kb   :ke+hk)

       real, intent(in)    :: bcvalue
       integer, intent(in)    :: bctype
       integer, intent(in)    :: n

       k  = zwallsshear(n,1)
       km = k-1
       il = zwallsshear(n,2)       ! starting i-index
       iu = zwallsshear(n,3)-1     ! ending i-index one less than shear components
       jl = zwallsshear(n,4)       ! starting j-index
       ju = zwallsshear(n,5)-1     ! ending j-index one less than shear components

       if (bctype == 1) then
         do j=jl,ju
           do i=il,iu

              putout(i,j,k) = putout(i,j,k) + ( &
                          0.5*(dzf(km)*ekh(i,j,k) + dzf(k)*ekh(i,j,km))* &                
                              (putin(i,j,k)-putin(i,j,km))*dzh2i(k) -    &          
                              (dzf(km)*ekh(i,j,k) + dzf(k)*ekh(i,j,km))* &                 
                              (putin(i,j,k)-bcvalue)*dzfi(k)*dzhi(k)     &                ! fixed value
!                              2.*prandtlmoli*numol* &                 
!                              (putin(i,j,k)-bcvalue)*dzfi(k)     &                ! fixed value
                                  )*dzfi(k)  
 
              putout(i,j,k-1) = putout(i,j,k-1) + ( &
                          - 0.5*(dzf(km)*ekh(i,j,k) + dzf(k)*ekh(i,j,km))* &       
                                (putin(i,j,k)-putin(i,j,km))*dzh2i(k) +    &
                                (dzf(km)*ekh(i,j,k) + dzf(k)*ekh(i,j,km))* &
                                (bcvalue-putin(i,j,km))*dzfi(km)*dzhi(k)   &              ! fixed value  
!                                 2.*prandtlmoli*numol* &
!                                (bcvalue-putin(i,j,km))*dzfi(km)   &              ! fixed value  
                                   )*dzfi(km) 

           end do
         end do
       else  ! bctype == 0 
         do j=jl,ju
           do i=il,iu
             putout(i,j,k) = putout(i,j,k) + ( &
                         0.5*(dzf(km)*ekh(i,j,k) + dzf(k)*ekh(i,j,km))* &                 ! zero flux
                             (putin(i,j,k)-putin(i,j,km))*dzh2i(k) -     &
                              bcvalue)*dzfi(k) 

     ! this is correct, but removing a scalar from inside the block (without
     ! replenishing it) leads to problems, e.g. negative numbers. Which is
     ! especially bad for temperature etc..
     !        putout(i,j,k-1) = putout(i,j,k-1) + ( &
     !                    - 0.5*(dzf(km)*ekh(i,j,k) + dzf(k)*ekh(i,j,km))* &               ! zero flux
     !                          (putin(i,j,k)-putin(i,j,km))*dzh2i(k) +    &
     !                           bcvalue)*dzfi(km) 
           end do
         end do
       end if

     end subroutine zwallscalar



     subroutine zwalle12(n)
       use modglobal, only : dxhi,dzhi,dyi,dzf,dzfi,dzh2i,ib,ie,ih,jb,je,jh,kb,ke,kh,numol
       use modfields, only : e120, e12p,u0,v0,w0
       use modsubgriddata, only : ekm
       integer i,j,k,il,iu,jl,ju,im,ip,jm,jp,km,kp

       integer, intent(in)    :: n

       k  = zwallsshear(n,1)
       km = k-1
       kp = k+1
       il = zwallsshear(n,2)       ! starting i-index
       iu = zwallsshear(n,3)-1     ! ending i-index one less than shear components
       jl = zwallsshear(n,4)       ! starting j-index
       ju = zwallsshear(n,5)-1     ! ending j-index one less than shear components


       do j=jl,ju
            jm = j-1
            jp = j+1
         do i=il,iu
            im = i-1
            ip = i+1
            ! walls
            e12p(i,j,k) = e12p(i,j,k) + &
                           ( 0.5*(dzf(km)*ekm(i,j,k) + dzf(k)*ekm(i,j,km))* &                 ! zero flux
                           (e120(i,j,k)-e120(i,j,km))*dzh2i(k)*dzfi(k) )

            e12p(i,j,km) = e12p(i,j,km) + &
                           - ( 0.5*(dzf(km)*ekm(i,j,k) + dzf(k)*ekm(i,j,km))* &               ! zero flux
                           (e120(i,j,k)-e120(i,j,km))*dzh2i(k)*dzfi(km) )
            ! sources
            e12p(i,j,k) = e12p(i,j,k) + 0.25*(ekm(i,j,k)-numol)/(2*e120(i,j,k))*(&
             - ((w0(i,j,k)-w0(im,j,k))     * dxhi(i)       + &          ! subtract old term
               (u0(i,j,k)-u0(i,j,km))      * dzhi(k) )**2  + &
               ((w0(i,j,k)-w0(im,j,k))     * dxhi(i)       + &
               (2.*u0(i,j,k)               * dzfi(k)))**2  + &          ! add IBM term
              
             - ((w0(ip,j,k)-w0(i,j,k))     * dxhi(ip)      + &
               (u0(ip,j,k)-u0(ip,j,km))    * dzhi(k) )**2  + &
               ((w0(ip,j,k)-w0(i,j,k))     * dxhi(ip)      + &
               (2.*u0(ip,j,k)              * dzfi(k)))**2  + &

             - ((v0(i,j,k)-v0(i,j,km))     * dzhi(k)       + &
               (w0(i,j,k)-w0(i,jm,k))      * dyi     )**2  + &
               (2.*v0(i,j,k)               * dzfi(k)       + &
               (w0(i,j,k)-w0(i,jm,k))      * dyi     )**2  + &

             - ((v0(i,jp,k)-v0(i,jp,km))   * dzhi(k)       + &
               (w0(i,jp,k)-w0(i,j,k))      * dyi     )**2  + &
               (2.*v0(i,jp,k)              * dzfi(k)       + &
               (v0(i,jp,k)-v0(i,jp,km))    * dzhi(k) )**2    )

            e12p(i,j,km) = e12p(i,j,km) + 0.25*(ekm(i,j,km)-numol)/(2*e120(i,j,km))*(&
             - ((w0(i,j,k)-w0(im,j,k))     * dxhi(i)       + &
               (u0(i,j,k)-u0(i,j,km))      * dzhi(k)  )**2 + &
               ((w0(i,j,k)-w0(im,j,k))     * dxhi(i)       + &
               (-2.*u0(i,j,km)             * dzfi(km)))**2 + &
             
             - ((w0(ip,j,k)-w0(i,j,k))     * dxhi(ip)      + &
               (u0(ip,j,k)-u0(ip,j,km))    * dzhi(k)  )**2 + &
               ((w0(ip,j,k)-w0(i,j,k))     * dxhi(ip)      + &
               (-2.*u0(ip,j,km)            * dzfi(km)))**2 + &

             - ((v0(i,j,k)-v0(i,j,km))     * dzhi(k)       + &
               (w0(i,j,k)-w0(i,jm,k))      * dyi      )**2 + &
               (-2.*v0(i,j,km)             * dzfi(km)      + &
               (w0(i,j,k)-w0(i,jm,k))      * dyi      )**2 + &

             - ((v0(i,jp,k)-v0(i,jp,km))   * dzhi(k)       + &
               (w0(i,jp,k)-w0(i,j,k))      * dyi      )**2 + &
               (-2.*v0(i,jp,km)            * dzfi(km)      + &
               (w0(i,jp,k)-w0(i,j,k))      * dyi      )**2   )

         end do
       end do

     end subroutine zwalle12


     subroutine ibmforce
      use modglobal, only :ib,ie,ih,jb,je,jh,kb,ke,kh,rk3step,dt,libm,jmax,&
                            nblocks,nsv,ltempeq,lmoist,nsv,rk3step,ih,kh,dt,totavtime,&
                            dxh,dzf,dy,jgb,jge,timee,btime,startmean
      use modfields, only : up,vp,wp,um,vm,wm,svp,svm,thlp,qtp,dpdxl
      use modgenstat, only: timecompl
      use modmpi,     only: cmyid,myid
      real, dimension(ib-ih:ie+ih,kb-kh:ke+kh)          ::  dummy
      real rk3coef,rk3coefi,timecomplibm, timecomplibmplusdti,ylength_i,&
           rk3coefold,rk3dt,rk3dt_i,ibmxforcevoltot
      integer n,i,j,k,jl,ju,kl,ku,rk3stepold

!      if (libm == .true. .and. rk3step==3) then
      if (libm == .true. .and. timee-btime >= startmean) then
        rk3coef    = dt / (4. - dble(rk3step))
        rk3stepold = rk3step -1
        rk3coefold = dt / (4. - dble(rk3stepold))
        if (rk3stepold == 0) rk3coefold = 0.
        rk3dt = rk3coef - rk3coefold 
        rk3dt_i = 1. / rk3dt
        rk3coefi = 1. / rk3coef
!        timecomplibm     = timecompl -totavtime+rk3coefold
        timecomplibm     = timecompl +rk3coefold
        timecomplibmplusdti = 1./(timecomplibm+rk3dt)
!        if (myid==0) write(6,*) 'rk3coef, rk3coefold, rk3dt = ', rk3coef, rk3coefold, rk3dt
!        if (myid==0) write(6,*) 'timecompl = ', timecompl
!        write(6,*) 'rk3coef, rk3coefold, rk3dt = ', rk3coef, rk3coefold, rk3dt
 
        dummy(:,:) = 0.        ! set to zero initially

        do n=1,nxwallsnorm
          i  = xwallsnorm(n,1)
          jl = xwallsnorm(n,2)
          ju = xwallsnorm(n,3)
          kl = xwallsnorm(n,4)
          ku = xwallsnorm(n,5)
  ! first diagnostics: ibm force (used in canopy (modsave.f90))
          ylength_i = 1./((jge-jgb+1)*dy)
          do k=kl,ku
            do j=jl,ju
              dummy(i,k) =  dummy(i,k) + (-up(i,j,k)-um(i,j,k)*rk3coefi)*dxh(i)*dzf(k)*dy*ylength_i   ! This is a spanwise average
            end do
          end do
!          pup(i,jl:ju,kl:ku) = 0.
        end do  ! 1,nxwallsnorm
!        ibmxforce(:,:)=(ibmxforce(:,:)*(timecomplibm)+dummy(:,:)*dt)*timecomplibmplusdti
        ibmxforcevol(:,:)=(ibmxforcevol(:,:)*(timecomplibm)+dummy(:,:)*rk3dt)*timecomplibmplusdti
!        ibmxforce(:,:) = dummy(:,:)
!        ibmxforcevoltot = sum(sum(ibmxforcevol(ib:ie,kb:ke),1),1)
!        open(unit=11,file='ibmxforcevoltot'//cmyid//'.txt',position='append')
!        write(11,'(2e14.6)') timecomplibm+rk3dt,ibmxforcevoltot
!        close(11)

      end if ! libm
     end subroutine ibmforce




 
    subroutine ibmnorm
      use modglobal, only : ib,ie,ih,jb,je,jh,kb,ke,kh,rk3step,dt,libm,jmax,&
                            nblocks,nsv,ltempeq,lmoist,nsv,rk3step,ih,kh,dt,totavtime,&
                            dxh,dzf,dy
      use modfields, only : up,vp,wp,um,vm,wm,svp,svm,thlp,qtp
      use modmpi, only    : myid
      use modgenstat, only: timecompl
      real, dimension(ib-ih:ie+ih,kb-kh:ke+kh)          ::  dummy
      real rk3coef,rk3coefi,timecomplibm, timecomplibmplusdti
      integer n,i,j,k,il,iu,jl,ju,kl,ku,sc
      
      if (libm == .true.) then
        rk3coef = dt / (4. - dble(rk3step))
        rk3coefi = 1. / rk3coef
!        timecomplibm     = timecompl -totavtime
!        timecomplibmplusdti = 1./(timecomplibm+dt)
!        dummy(:,:) = 0.        ! set to zero initially
 

        do n=1,nxwallsnorm
          i  = xwallsnorm(n,1)
          jl = xwallsnorm(n,2)
          ju = xwallsnorm(n,3)
          kl = xwallsnorm(n,4)
          ku = xwallsnorm(n,5)
!  ! first diagnostics: ibm force (used in canopy (modsave.f90))
!          do k=kl,ku
!            do j=jl,ju
!              dummy(i,k) =  dummy(i,k) + (-up(i,j,k)-um(i,j,k)*rk3coefi)*dxh(i)*dzf(k)*dy   ! spanwise average should be taken in canopy (modsave.f90)
!            end do
!          end do
          up(i,jl:ju,kl:ku) = -um(i,jl:ju,kl:ku)*rk3coefi
!          pup(i,jl:ju,kl:ku) = 0.
        end do  ! 1,nxwallsnorm

!        if (rk3step==3) then
!          ibmxforce(:,:)=(ibmxforce(:,:)*(timecomplibm)+dummy(:,:)*dt)*timecomplibmplusdti
!        end if


        do n=1,nywallsnorm
          j  = ywallsnorm(n,1)
          il = ywallsnorm(n,2)
          iu = ywallsnorm(n,3)
          kl = ywallsnorm(n,4)
          ku = ywallsnorm(n,5)
!          write(6,*) 'myid,j=', myid,j 
          vp(il:iu,j,kl:ku) = -vm(il:iu,j,kl:ku)*rk3coefi
!          pvp(il:iu,j,kl:ku) = 0.
        end do
 
        do n=1,nzwallsnorm
          k  = zwallsnorm(n,1)
          il = zwallsnorm(n,2)
          iu = zwallsnorm(n,3)
          jl = zwallsnorm(n,4)
          ju = zwallsnorm(n,5)
          
          wp(il:iu,jl:ju,k) = -wm(il:iu,jl:ju,k)*rk3coefi
!          pwp(il:iu,jl:ju,k) = 0.
        end do

       
        if (ltempeq==.true. .and. lmoist==.true. .and. nsv>0) then
          do n=1,nblocks
            il = block(n,1)
            iu = block(n,2)
            kl = block(n,5)
            ku = block(n,6)
            jl = block(n,3)-myid*jmax
            ju = block(n,4)-myid*jmax
            if (ju < jb .or. jl > je) then
              cycle
            else
              if (ju > je) ju=je
              if (jl < jb) jl=jb
              thlp(il:iu,jl:ju,kl:ku)   = 0.
              svp (il:iu,jl:ju,kl:ku,:) = 0.
              qtp(il:iu,jl:ju,kl:ku)     = 0.
            end if
          end do
        elseif  (nsv>0 .and. lmoist==.true.) then
          do n=1,nblocks
            il = block(n,1)
            iu = block(n,2)
            kl = block(n,5)
            ku = block(n,6)
            jl = block(n,3)-myid*jmax
            ju = block(n,4)-myid*jmax
            if (ju < jb .or. jl > je) then
              cycle
            else
              if (ju > je) ju=je
              if (jl < jb) jl=jb
              svp (il:iu,jl:ju,kl:ku,:) = 0.
              qtp(il:iu,jl:ju,kl:ku)     = 0.
            end if
          end do
        elseif  (ltempeq==.true. .and. lmoist==.true.) then
          do n=1,nblocks
            il = block(n,1)
            iu = block(n,2)
            kl = block(n,5)
            ku = block(n,6)
            jl = block(n,3)-myid*jmax
            ju = block(n,4)-myid*jmax
            if (ju < jb .or. jl > je) then
              cycle
            else
              if (ju > je) ju=je
              if (jl < jb) jl=jb
              thlp(il:iu,jl:ju,kl:ku)   = 0.
              qtp(il:iu,jl:ju,kl:ku)     = 0. 
            end if
              
          end do
         elseif  (nsv>0 .and. ltempeq==.true.) then
          do n=1,nblocks
            il = block(n,1)
            iu = block(n,2)
            kl = block(n,5)
            ku = block(n,6)
            jl = block(n,3)-myid*jmax
            ju = block(n,4)-myid*jmax
            if (ju < jb .or. jl > je) then
              cycle
            else
              if (ju > je) ju=je
              if (jl < jb) jl=jb
              thlp(il:iu,jl:ju,kl:ku)   = 0.
              svp (il:iu,jl:ju,kl:ku,:) = 0.
            end if
          end do

          else
           if (nsv>0) then
           do n=1,nblocks
            il = block(n,1)
            iu = block(n,2)
            kl = block(n,5)
            ku = block(n,6)
            jl = block(n,3)-myid*jmax
            ju = block(n,4)-myid*jmax
            if (ju < jb .or. jl > je) then
              cycle
            else
              if (ju > je) ju=je
              if (jl < jb) jl=jb
              svp (il:iu,jl:ju,kl:ku,:) = 0.
            end if
           end do

           elseif (ltempeq) then
           do n=1,nblocks
            il = block(n,1)
            iu = block(n,2)
            kl = block(n,5)
            ku = block(n,6)
            jl = block(n,3)-myid*jmax
            ju = block(n,4)-myid*jmax
            if (ju < jb .or. jl > je) then
              cycle
            else
              if (ju > je) ju=je
              if (jl < jb) jl=jb
              thlp(il:iu,jl:ju,kl:ku)   = 0.
            end if
           end do

           elseif (lmoist) then
           do n=1,nblocks
            il = block(n,1)
            iu = block(n,2)
            kl = block(n,5)
            ku = block(n,6)
            jl = block(n,3)-myid*jmax
            ju = block(n,4)-myid*jmax
            if (ju < jb .or. jl > je) then
              cycle
            else
              if (ju > je) ju=je
              if (jl < jb) jl=jb
              qtp(il:iu,jl:ju,kl:ku)   = 0.
            end if
           end do
           end if

        end if 

      end if   ! libm == .true.

    end subroutine ibmnorm

!> Determines the distance to the nearest wall for each cell center (used in v. Driest damping function)
!> Output is a field with minimal distance to wall for each cell center
  subroutine nearwall

    use modglobal, only : ib,ie,jb,je, jgb,jge,jmax,kb,ke,xh,xf,dy,zh,zf,lwarmstart,&
                          nxwalls,nywalls,nzwalls,nblocks,libm,lzerogradtop,lwalldist
    use modsubgriddata, only : lsmagorinsky,loneeqn
    use modfields, only : mindist,wall
    use modibmdata, only: xwallsglobal,ywallsglobal,zwallsglobal,block
    use modmpi, only    : myid
    implicit none

    integer, allocatable :: ux0all(:,:,:), vy0all(:,:,:), wz0all(:,:,:)
    real, allocatable :: distxf(:,:), distxh(:,:), distyf(:,:), distyh(:,:), distzf(:,:), distzh(:,:),&
                         distxf2(:,:), distxh2(:,:), distyf2(:,:), distyh2(:,:), distzf2(:,:), distzh2(:,:), distance(:)
    real distx,disty,distz           ! distx is the distance to nearest x-wall, etc.
!    integer, allocatable :: optie(:)
    integer ic,jc,kc,i,j,k,optie,il,iu,jl,ju,kl,ku,n
    
!  if (lwarmstart==.true. .or. lles==.false. .or. lvreman==.true.) then
  if ((lsmagorinsky==.true. .or. loneeqn==.true.) .and. lwalldist==.true.) then
  if (myid==0) then
    write(6,*) 'Computing wall distances' 
  end if
!  if (lles==.false. .or. lvreman==.true.) then
  mindist=1.0e10
  
  allocate(ux0all(ib-1:ie+1,jgb-1:jge+1,kb-1:ke+1))   ! This contains ux0 + the lower and (possibly) the upper wall
  allocate(vy0all(ib-1:ie+1,jgb-1:jge+1,kb-1:ke+1))   ! This contains ux0 + the lower and (possibly) the upper wall
  allocate(wz0all(ib-1:ie+1,jgb-1:jge+1,kb-1:ke+1))   ! This contains ux0 + the lower and (possibly) the upper wall
  allocate(distxh(ib:ie,ib:ie+1))
  allocate(distxf(ib:ie,ib:ie+1))
  allocate(distyh(jb:je,jgb:jge+1))
  allocate(distyf(jb:je,jgb:jge+1))
  allocate(distzh(kb:ke,kb:ke+1))
  allocate(distzf(kb:ke,kb:ke+1))
  allocate(distxh2(ib:ie,ib:ie+1))
  allocate(distxf2(ib:ie,ib:ie+1))
  allocate(distyh2(jb:je,jgb:jge+1))
  allocate(distyf2(jb:je,jgb:jge+1))
  allocate(distzh2(kb:ke,kb:ke+1))
  allocate(distzf2(kb:ke,kb:ke+1))
  allocate(distance(4))
  
  ! initialize wall indicators
  ux0all = 0
  vy0all = 0
  wz0all = 0


! Determine for each cell face if an x/y/z-wall is present
! from immersed boundaries

  if (libm == .true.) then

! loop over global walls
!  write(6,*) 'xwallsglobal(1,1)=',xwallsglobal(1,1)
  do n=1,nxwalls
    i=xwallsglobal(n,1)
    do k=xwallsglobal(n,4),xwallsglobal(n,5)
      do j=xwallsglobal(n,2),xwallsglobal(n,3)
        ux0all(i,j,k) = 1 
      end do
    end do   
  end do
  
  do n=1,nywalls
    j=ywallsglobal(n,1)
    do k=ywallsglobal(n,4),ywallsglobal(n,5)
      do i=ywallsglobal(n,2),ywallsglobal(n,3)
        vy0all(i,j,k) = 1
      end do
    end do
  end do

  do n=1,nzwalls
    k=zwallsglobal(n,1)
    do j=zwallsglobal(n,4),zwallsglobal(n,5)
      do i=zwallsglobal(n,2),zwallsglobal(n,3)
        wz0all(i,j,k) = 1
      end do
    end do
  end do  

  ! do loop over blocks
  do n=1,nblocks
    il = block(n,1)
    iu = block(n,2)
    jl = block(n,3)
    ju = block(n,4)
    kl = block(n,5)
    ku = block(n,6)
    do k=kl,ku
      do j=jl,ju
        ux0all(il,j,k)=1      ! lower x-wall
        ux0all(iu+1,j,k)=1    ! upper x-wall
      end do
    end do
    do k=kl,ku
      do i=il,iu
        vy0all(i,jl,k)=1      ! lower y-wall
        vy0all(i,ju+1,k)=1    ! upper y-wall
      end do
    end do
    do j=jl,ju
      do i=il,iu
        wz0all(i,j,kl)=1      ! lower z-wall
        wz0all(i,j,ku+1)=1    ! upper z-wall
      end do
    end do
  end do   ! loop over nblocks

  end if ! libm = .true.

! add the global walls (probably upper and lower wall, or only lower wall)
  if (lzerogradtop==.true.) then
    do i=ib,ie    
    do j=jgb,jge
      wz0all(i,j,kb)=1           ! ground wall
    end do
    end do
  else 
    do i=ib,ie    
    do j=jgb,jge
      wz0all(i,j,kb)=1           ! ground wall
      wz0all(i,j,ke+1)=1;        ! top wall
    end do
    end do
  end if
 

  write(6,*) 'Determing distance matrices, proc=', myid
! Determine x-distance matrices:
  do ic=ib,ie    ! cell-center index
  do i=ib,ie+1   ! vertex-index (1 more than cell centers)
    distxh(ic,i) = xf(ic)-xh(i)
  end do
  end do 
 
  do ic=ib,ie    ! cell-center index
  do i=ib,ie+1   ! center -index 
    distxf(ic,i) = xf(ic)-xf(i)
  end do
  end do  

! Determine y-distance matrices:
  do jc=jb,je    ! cell-center index
  do j=jgb,jge+1 ! vertex-index (1 more than cell centers) (global index to make sure distance to all cells is determined)
    distyh(jc,j) = (jc+myid*jmax-j)*dy+0.5*dy
  end do
  end do 
  
  do jc=jb,je    ! cell-center index
  do j=jgb,jge+1 ! center-index  (global index to make sure distance to all cells is determined)
    distyf(jc,j) = (jc+myid*jmax-j)*dy
  end do
  end do 

! Determine z-distance matrices:
  do kc=kb,ke    ! cell-center index
  do k=kb,ke+1   ! vertex-index (1 more than cell centers)
    distzh(kc,k) = zf(kc)-zh(k)
  end do
  end do  
  
  do kc=kb,ke    ! cell-center index
  do k=kb,ke+1   ! vertex-index (1 more than cell centers)
    distzf(kc,k) = zf(kc)-zf(k)
  end do
  end do  

  distxh2=distxh**2
  distyh2=distyh**2
  distzh2=distzh**2
  distxf2=distxf**2
  distyf2=distyf**2
  distzf2=distzf**2


  write(6,*) 'Finished determing distance matrices, proc=', myid
  write(6,*) 'determing distance to nearest wall for each cell center, proc=', myid


! Loop over cells (ic,jc,kc) for which minimal wall-to-cell-center-distance needs to be determined
!  do jc=jgb,jge
  do kc=kb,ke
  do jc=jb,je
  do ic=ib,ie
! Determine distance between cc of cell (ic,jc,kc) and faces of all cells (i,j,k)
    do k=kb,ke+1            ! Level ke+1 is computed in a separate loop (only necessary with upper wall-> global approach=faster)
    do j=jgb,jge+1        ! loop goes up to jge+1 because jge+1 contains the last vy0-wall
    do i=ib,ie+1          ! loop goes up to ie+1 because ie+1 contains the last ux0-wall
      if (ux0all(i,j,k)==1 .OR. vy0all(i,j,k)==1 .OR. wz0all(i,j,k)==1) then
        distx=1.0e10   ! make sure distx is very large when no x-wall is present
        disty=1.0e10   ! make sure disty is very large when no y-wall is present
        distz=1.0e10   ! make sure distz is very large when no z-wall is present
        if (ux0all(i,j,k)==1) then
          distx = sqrt(distxh2(ic,i)+distyf2(jc,j)+distzf2(kc,k))
        end if
        if (vy0all(i,j,k)==1) then
          disty = sqrt(distxf2(ic,i)+distyh2(jc,j)+distzf2(kc,k))
        end if
        if (wz0all(i,j,k)==1) then
          distz = sqrt(distxf2(ic,i)+distyf2(jc,j)+distzh2(kc,k))
        end if
      else          ! no walls are present in cell (i,j,k) -> distance does not need to be determined for this cell
        cycle       ! go to next cell (i,j,k)
      end if 
! determine minimal wall distance between cc of (ic,jc,kc) and faces of cell (i,j,k)
      distance = (/ mindist(ic,jc,kc),distx,disty,distz /)
      optie = minloc(distance,1)
!        write(6,*) 'optie=', optie 

      if (optie==1) then
        cycle
      else if (optie==2) then
        mindist(ic,jc,kc) = distx
        wall(ic,jc,kc,2) = j
        wall(ic,jc,kc,3) = k
!        wall(ic,jc,kc,4) = 1     ! This means the wall closest to the cc of (ic,jc,kc) is at an x-wall at (i,j,k)        
        if (ic>=i) then
          wall(ic,jc,kc,1) = i
          wall(ic,jc,kc,4) = 5   ! shear component index (variable: shear)
          wall(ic,jc,kc,5) = 9   ! shear component index (variable: shear)
        else 
          wall(ic,jc,kc,1) = i-1 ! in the subgrid this stress is computed in the cell i-1
          wall(ic,jc,kc,4) = 6   ! shear component index (variable: shear)
          wall(ic,jc,kc,5) = 10   ! shear component index (variable: shear)
        end if
      else if (optie==3) then
        mindist(ic,jc,kc) = disty
        wall(ic,jc,kc,1) = i
        wall(ic,jc,kc,3) = k
!        wall(ic,jc,kc,4) = 2     ! This means the wall closest to the cc of (ic,jc,kc) is at a y-wall at (i,j,k)
        if (jc+myid*jmax>=j) then
          wall(ic,jc,kc,2) = j
          wall(ic,jc,kc,4) = 1   ! shear component index (variable: shear)
          wall(ic,jc,kc,5) = 11  ! shear component index (variable: shear)
        else 
          wall(ic,jc,kc,2) = j-1  ! in the subgrid this stress is computed in the cell j-1
          wall(ic,jc,kc,4) = 2   ! shear component index (variable: shear)
          wall(ic,jc,kc,5) = 12  ! shear component index (variable: shear)
        end if
      else if (optie==4) then
        mindist(ic,jc,kc) = distz
        wall(ic,jc,kc,1) = i
        wall(ic,jc,kc,2) = j
!        wall(ic,jc,kc,4) = 3     ! This means the wall closest to the cc of (ic,jc,kc) is at a z-wall at (i,j,k)
        if (kc>=k) then
          wall(ic,jc,kc,3) = k
          wall(ic,jc,kc,4) = 3   ! shear component index (variable: shear)
          wall(ic,jc,kc,5) = 7   ! shear component index (variable: shear)
        else 
          wall(ic,jc,kc,3) = k-1  ! in the subgrid this stress is computed in the cel k-1
          wall(ic,jc,kc,4) = 4   ! shear component index (variable: shear)
          wall(ic,jc,kc,5) = 8   ! shear component index (variable: shear)
        end if
      end if
!      mindist(ic,jc+myid*jmax,kc)=min(mindist(ic,jc+myid*jmax,kc),distx,disty,distz)   ! global j index
    end do
    end do
    end do
!    if (myid==0) write(6,*) 'finished for cell (ic,jc,kc)=',ic,jc,kc 
  end do
  end do
  end do

  write(6,*) 'Finished determing distance to nearest wall for each cell center, proc=', myid
!  write(6,*) 'mindist(ib,jb+myid*jmax,kb),mindist(ib,je+myid*jmax,kb)',mindist(ib,jb+myid*jmax,kb),mindist(ib,je+myid*jmax,kb)

  else
    return
  end if !(lwarmstart)

  deallocate(ux0all,vy0all,wz0all)
  deallocate(xwallsglobal,ywallsglobal,zwallsglobal,block)   ! used for determining boundaries 
  return
  end subroutine nearwall

!> !! Only for 2D simulations!!!!
!> Determines the distance to the nearest wall for each cell center (used in Piomelli damping function)
!> Output is a field with minimal distance to wall for each cell center
  subroutine nearwall2d

    use modglobal, only : ib,ie,jb,je, jgb,jge,jmax,kb,ke,xh,xf,dy,zh,zf,lwarmstart,lzerogradtop
    use modfields, only : mindist,wall
    use modsubgriddata, only : lsmagorinsky,loneeqn
    use modmpi, only    : myid
    implicit none

    real, allocatable :: ux0all(:,:,:), vy0all(:,:,:), wz0all(:,:,:)
    real, allocatable :: distxf(:,:), distxh(:,:), distyf(:,:), distyh(:,:), distzf(:,:), distzh(:,:),&
                         distxf2(:,:), distxh2(:,:), distyf2(:,:), distyh2(:,:), distzf2(:,:), distzh2(:,:), distance(:)
    real distx,disty,distz           ! distx is the distance to nearest x-wall, etc.
!    integer, allocatable :: optie(:)
    integer ic,jc,kc,i,j,k,optie
    
!  if (lwarmstart==.true. .or. lles==.false.) then
!  if (lwarmstart==.false. .or. lles==.false.) then
  if (lsmagorinsky==.true. .or. loneeqn==.true.) then
!  mindist=1.0e10    ! this should be removed next run and 3 lines above should be uncommented

  allocate(ux0all(ib-1:ie+1,jgb-1:jge+1,kb-1:ke+1))   ! This contains ux0 + the lower and (possibly) the upper wall
  allocate(wz0all(ib-1:ie+1,jgb-1:jge+1,kb-1:ke+1))   ! This contains ux0 + the lower and (possibly) the upper wall
  allocate(distxh(ib:ie,ib:ie+1))
  allocate(distxf(ib:ie,ib:ie+1))
  allocate(distzh(kb:ke,kb:ke+1))
  allocate(distzf(kb:ke,kb:ke+1))
  allocate(distxh2(ib:ie,ib:ie+1))
  allocate(distxf2(ib:ie,ib:ie+1))
  allocate(distzh2(kb:ke,kb:ke+1))
  allocate(distzf2(kb:ke,kb:ke+1))
  allocate(distance(3))
!  allocate(optie(4))

  
!  ux0all = ux0
!  wz0all = wz0

! add the global walls (probably upper and lower wall)
  if (lzerogradtop ==.true.) then
    do i=ib,ie    
    do j=jgb,jge
      wz0all(i,j,kb)=1;           ! ground wall
    end do
    end do
  else
    do i=ib,ie    
    do j=jgb,jge
      wz0all(i,j,kb)=1;           ! ground wall
      wz0all(i,j,ke+1)=1;         ! top wall
    end do
    end do
  end if
 

  write(6,*) 'Determing distance matrices, proc=', myid
! Determine x-distance matrices:
  do ic=ib,ie    ! cell-center index
  do i=ib,ie+1   ! vertex-index (1 more than cell centers)
    distxh(ic,i) = xf(ic)-xh(i)
  end do
  end do 
 
  do ic=ib,ie    ! cell-center index
  do i=ib,ie+1   ! center -index 
    distxf(ic,i) = xf(ic)-xf(i)
  end do
  end do  

  
! Determine z-distance matrices:
  do kc=kb,ke    ! cell-center index
  do k=kb,ke+1   ! vertex-index (1 more than cell centers)
    distzh(kc,k) = zf(kc)-zh(k)
  end do
  end do  
  
  do kc=kb,ke    ! cell-center index
  do k=kb,ke+1   ! vertex-index (1 more than cell centers)
    distzf(kc,k) = zf(kc)-zf(k)
  end do
  end do  

  distxh2=distxh**2
  distzh2=distzh**2
  distxf2=distxf**2
  distzf2=distzf**2


  write(6,*) 'Finished determing distance matrices, proc=', myid
  write(6,*) 'determing distance to nearest wall for each cell center, proc=', myid


! Loop over cells (ic,jc,kc) for which minimal wall-to-cell-center-distance needs to be determined
!  do jc=jgb,jge
  do kc=kb,ke
  do ic=ib,ie
! Determine distance between cc of cell (ic,jc,kc) and faces of all cells (i,j,k)
    do k=kb,ke+1            ! Level ke+1 is computed in a separate loop (only necessary with upper wall-> global approach=faster)
    do i=ib,ie+1          ! loop goes up to ie+1 because ie+1 contains the last ux0-wall
      if (ux0all(i,jb,k)==1 .OR. wz0all(i,jb,k)==1) then    ! assume this holds for all j
        distx=1.0e10   ! make sure distx is very large when no x-wall is present
        distz=1.0e10   ! make sure distz is very large when no z-wall is present
        if (ux0all(i,jb,k)==1) then
          distx = sqrt(distxh2(ic,i)+distzf2(kc,k))
        end if
        if (wz0all(i,jb,k)==1) then
          distz = sqrt(distxf2(ic,i)+distzh2(kc,k))
        end if
      else          ! no walls are present in cell (i,j,k) -> distance does not need to be determined for this cell
        cycle       ! go to next cell (i,j,k)
      end if 
! determine minimal wall distance between cc of (ic,jc,kc) and faces of cell (i,j,k)
      distance = (/ mindist(ic,jb,kc),distx,distz /)
      optie = minloc(distance,1)
      
      if (optie==1) then
        cycle
      else if (optie==2) then
        mindist(ic,:,kc) = distx
        wall(ic,:,kc,2) = jb
        wall(ic,:,kc,3) = k
!        wall(ic,jc,kc,4) = 1     ! This means the wall closest to the cc of (ic,jc,kc) is at an x-wall at (i,j,k)        
        if (ic>=i) then
          wall(ic,:,kc,1) = i
          wall(ic,:,kc,4) = 5   ! shear component index (variable: shear)
          wall(ic,:,kc,5) = 9   ! shear component index (variable: shear)
        else 
          wall(ic,:,kc,1) = i-1 ! in the subgrid this stress is computed in the cell i-1
          wall(ic,:,kc,4) = 6   ! shear component index (variable: shear)
          wall(ic,:,kc,5) = 10   ! shear component index (variable: shear)
        end if
      else if (optie==3) then
        mindist(ic,:,kc) = distz
        wall(ic,:,kc,1) = i
        wall(ic,:,kc,2) = jb
!        wall(ic,jc,kc,4) = 3     ! This means the wall closest to the cc of (ic,jc,kc) is at a z-wall at (i,j,k)
        if (kc>=k) then
          wall(ic,:,kc,3) = k
          wall(ic,:,kc,4) = 3   ! shear component index (variable: shear)
          wall(ic,:,kc,5) = 7   ! shear component index (variable: shear)
        else 
          wall(ic,:,kc,3) = k-1  ! in the subgrid this stress is computed in the cel k-1
          wall(ic,:,kc,4) = 4   ! shear component index (variable: shear)
          wall(ic,:,kc,5) = 8   ! shear component index (variable: shear)
        end if
      end if
!      mindist(ic,jc+myid*jmax,kc)=min(mindist(ic,jc+myid*jmax,kc),distx,disty,distz)   ! global j index
    end do
    end do
  end do
  end do

  write(6,*) 'Finished determing distance to nearest wall for each cell center, proc=', myid
  else
    return
  end if !(lwarmstart)

  deallocate(ux0all,wz0all)
  return
  end subroutine nearwall2d



end module modibm
