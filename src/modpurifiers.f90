!> \file modpurifs.f90                                                                              
!!tg3315, 3 Jul 2017 

!> Input air purifiers into DALES model.

module modpurifiers
use mpi
implicit none
save

contains
    subroutine createpurifiers
    use modglobal,  only : lpurif,npurif,purif,cexpnr,ifinput
    use modmpi,     only : myid,comm3d,mpierr

    implicit none
    integer :: n
    character(80) chmess

    if (lpurif .eqv. .false.) return
    
      allocate(purif(npurif,7))

      ! read global purifiers
      if(myid==0) then
        ! write(*,*) '1, myid, npurif, lpurif, cexpnr', myid, npurif, lpurif, cexpnr
        if (npurif>0) then
          open (ifinput,file='purifs.inp.'//cexpnr)
          read (ifinput,'(a80)') chmess
          read (ifinput,'(a80)') chmess
          do n=1,npurif
            read (ifinput,*) &
                  purif(n,1), &
                  purif(n,2), &
                  purif(n,3), &
                  purif(n,4), &
                  purif(n,5), &
                  purif(n,6), & 
                  purif(n,7)
          end do

          ! write (6,*) 'Purifier number,  il, iu, jl, ju, kl, ku, ipu '
          ! do n=1,npurif
          !   write (6,*) &
          !         n , &
          !         purif(n,1), &
          !         purif(n,2), &
          !         purif(n,3), &
          !         purif(n,4), &
          !         purif(n,5), &
          !         purif(n,6), &
          !         purif(n,7)                 
          ! end do
        end if
      ! write(*,*) 'Finished determining purifiers on myid == 0'
      end if ! end if myid==0

      call MPI_BCAST(purif  ,7*npurif,MPI_INTEGER ,0,comm3d,mpierr)
      call MPI_BCAST(npurif ,1,MPI_INTEGER ,0,comm3d,mpierr)

      end subroutine createpurifiers

      subroutine purifiers
      use modglobal,  only : ib,ie,jb,je,kb,ke,ih,xf,xh,zh,purif,npurif,lpurif,itot,jtot,Qpu,epu,dy,nsv
      use modfields,  only : um,vm,wm,u0,v0,w0,up,vp,wp,svp,svm,sv0
      use modmpi,     only : myidx,myidy,nprocx,nprocy

      implicit none
      integer :: i,j,k,n,il,iu,jl,ju,kl,ku
      real :: Apu,dpu,upu,vpu,wpu,udpu
      real, allocatable :: inpu  (:,:,:,:)

      if (lpurif .eqv. .false.) return

      do n = 1, npurif

        upu=0.
        vpu=0.
        wpu=0.

        ! Determine flow direction through purifier=
        select case(purif(n,7))

        ! calculate cross-sectional area, length of purifier and flow speed
        case(1)
        Apu = (purif(n,4)-purif(n,3)+1)*dy*(zh(purif(n,6)+1)-zh(purif(n,5)))
        dpu = xh(purif(n,2)+1) - xh(purif(n,1))
        upu = Qpu/Apu
        case(2)
        Apu = (purif(n,4)-purif(n,3)+1)*dy*(zh(purif(n,6)+1)-zh(purif(n,5)))
        dpu = xh(purif(n,2)+1) - xh(purif(n,1))
        upu = -Qpu/Apu
        case(3)
        Apu = (xh(purif(n,2)+1)-xh(purif(n,1)))*(zh(purif(n,6)+1)-zh(purif(n,5)))
        dpu = (purif(n,4)-purif(n,3)+1)*dy
        vpu = Qpu/Apu
        case(4)
        Apu = (xh(purif(n,2)+1)-xh(purif(n,1)))*(zh(purif(n,6)+1)-zh(purif(n,5)))
        dpu = (purif(n,4)-purif(n,3)+1)*dy
        vpu = -Qpu/Apu
        case(5)
        Apu = (purif(n,4)-purif(n,3)+1)*dy*(xh(purif(n,2)+1)-xh(purif(n,1)))
        dpu = zh(purif(n,6)+1) - zh(purif(n,5))
        wpu = Qpu/Apu
        case(6)
        Apu = (purif(n,4)-purif(n,3)+1)*dy*(xh(purif(n,2)+1)-xh(purif(n,1)))
        dpu = zh(purif(n,6)+1) - zh(purif(n,5))
        wpu = -Qpu/Apu
        case(7) ! purifier that takes from above and below and feeds outwards
        Apu = (purif(n,4)-purif(n,3)+1)*dy*(xh(purif(n,2)+1)-xh(purif(n,1)))*2 ! *2 due to geometry of purifier
        dpu = zh(purif(n,6)+1) - zh(purif(n,5))/2 ! half distance to travel before purification?
        udpu = Qpu/Apu ! new term to be used at purifier faces
        case(8) ! one cell purifier that intakes from all side and outputs upwards
        Apu = 2*(purif(n,4)-purif(n,3)+1)*dy*(xh(purif(n,2)+1)-xh(purif(n,1))) + &
              2*(xh(purif(n,2)+1)-xh(purif(n,1)))*(zh(purif(n,6)+1)-zh(purif(n,5))) ! sum up 4 horizontal cell faces
        dpu = 0.5*( (purif(n,4)-purif(n,3)+1)*dy + xh(purif(n,2)+1) - xh(purif(n,1))  ) ! average distance from outside cell to centre cell in horizontal directions
        udpu = Qpu/Apu ! new term to be used at purifier faces
        end select

        ! enforce flowrate

        ! u flowrate
        il = purif(n,1) - myidx*itot/nprocx
        iu = purif(n,2) + 1 - myidx*itot/nprocx
        kl = purif(n,5)
        ku = purif(n,6)
        jl = purif(n,3) - myidy*jtot/nprocy
        ju = purif(n,4) - myidy*jtot/nprocy

        if (iu < ib .or. il > ie .or. ju < jb .or. jl > je) then
          cycle
        else
          if (iu > ie) iu=ie
          if (il < ib) il=ib 
          if (ju > je) ju=je
          if (jl < jb) jl=jb 
 
          up(il:iu,jl:ju,kl:ku) = 0.
          um(il:iu,jl:ju,kl:ku) = upu
          u0(il:iu,jl:ju,kl:ku) = upu ! tg3315 !WARNING changed these to u0 as realised RK3... change back to um if not working... !tg3315 14.05.18 have made it do both due to role in RK3, ..m is more crucial to change and may need ..0         
        end if

        ! v flowrate
        il = purif(n,1) - myidx*itot/nprocx
        iu = purif(n,2) - myidx*itot/nprocx
        kl = purif(n,5)
        ku = purif(n,6)
        jl = purif(n,3) - myidy*jtot/nprocy
        ju = purif(n,4) + 1 - myidy*jtot/nprocy
        if (iu < ib .or. il > ie .or. ju < jb .or. jl > je) then
          cycle
        else
          if (iu > ie) iu=ie
          if (il < ib) il=ib 
          if (ju > je) ju=je
          if (jl < jb) jl=jb 

          vp(il:iu,jl:ju,kl:ku) = 0.
          vm(il:iu,jl:ju,kl:ku) = vpu
          v0(il:iu,jl:ju,kl:ku) = vpu
        end if

        ! w flowrate
        il = purif(n,1) - myidx*itot/nprocx
        iu = purif(n,2) - myidx*itot/nprocx
        kl = purif(n,5)
        ku = purif(n,6) + 1
        jl = purif(n,3) - myidy*jtot/nprocy
        ju = purif(n,4) - myidy*jtot/nprocy
        if (iu < ib .or. il > ie .or. ju < jb .or. jl > je) then
          cycle
        else
          if (iu > ie) iu=ie
          if (il < ib) il=ib 
          if (ju > je) ju=je
          if (jl < jb) jl=jb 

          wp(il:iu,jl:ju,kl:ku)  = 0.
          wm(il:iu,jl:ju,kl:ku)  = wpu
          w0(il:iu,jl:ju,kl:ku)  = wpu
        end if

        ! Scalars
        il = purif(n,1) - myidx*itot/nprocx
        iu = purif(n,2) - myidx*itot/nprocx
        kl = purif(n,5)
        ku = purif(n,6)
        jl = purif(n,3) - myidy*jtot/nprocy
        ju = purif(n,4) - myidy*jtot/nprocy
        if (iu < ib .or. il > ie .or. ju < jb .or. jl > je) then
          cycle
        else
          if (iu > ie) iu=ie
          if (il < ib) il=ib 
          if (ju > je) ju=je
          if (jl < jb) jl=jb

          where (sv0(:,:,:,1) < 0.) sv0(:,:,:,1)=0.     !must do this in tstep after svo = ... 

          ! calculate concentration at purifier inlet
          allocate(inpu(il:iu,jl:ju,kl:ku,1:nsv))
          inpu=0.
          select case(purif(n,7))
            case(1)
            do i=il,iu
              inpu(i,:,:,:) = svm(il-1,jl:ju,kl:ku,:) !tg3315 also changed svm -> sv0 change back if stops working ! tg3315 14.05.18 I have undone this now as believe that ..m represents the latest confirmed value whilst ..0 is updated through RK3 timestep. If causes issues then stop!
            end do
            case(2)
            do i=il,iu
              inpu(i,:,:,:) = svm(iu+1,jl:ju,kl:ku,:)
            end do
            case(3)
            do j=jl,ju
              inpu(:,j,:,:) = svm(il:iu,jl-1,kl:ku,:)
            end do
            case(4)
            do j=jl,ju
             inpu(:,j,:,:) = svm(il:iu,ju+1,kl:ku,:)
            end do
            case(5)
            do k=kl,ku
              inpu(:,:,k,:) = svm(il:iu,jl:ju,kl-1,:)
            end do
            case(6)
            do k=kl,ku
              inpu(:,:,k,:) = svm(il:iu,jl:ju,ku+1,:)
            end do
            case(7)
              inpu(:,:,ku,:) = svm(il:iu,jl:ju,ku+1,:)
              inpu(:,:,kl,:) = svm(il:iu,jl:ju,kl-1,:)
              wm(il:iu,jl:ju,ku+1) = - udpu
              w0(il:iu,jl:ju,ku+1) = - udpu
              wm(il:iu,jl:ju,kl)   = udpu
              w0(il:iu,jl:ju,kl)   = udpu
            case(8)
              !> initially coded for purifiers of cell size = 1 (il=iu,jl=ju,kl=ku)
              inpu(il,jl,:,:) = 0.25 * ( svm(il-1,jl,kl:ku,:) + svm(iu+1,jl,kl:ku,:) + svm(il,jl-1,kl:ku,:) + svm(il,ju+1,kl:ku,:) )
              um(il,jl:ju,kl:ku) = udpu
              u0(il,jl:ju,kl:ku) = udpu
              vm(il:iu,jl,kl:ku) = udpu
              v0(il:iu,jl,kl:ku) = udpu
              um(iu+1,jl:ju,kl:ku) = - udpu
              u0(iu+1,jl:ju,kl:ku) = - udpu
              vm(il:iu,ju+1,kl:ku) = - udpu
              v0(il:iu,ju+1,kl:ku) = - udpu
          end select

          ! apply sink term to purify at given efficiency
          svp(il:iu,jl:ju,kl:ku,1) = svp(il:iu,jl:ju,kl:ku,1) - &
                                       (Qpu/Apu) * epu * inpu(:,:,:,1) / dpu

          if (nsv>1) then
            svp(il:iu,jl:ju,kl:ku,2) = svp(il:iu,jl:ju,kl:ku,2) - &
                                         (Qpu/Apu) * 0.7 * inpu(:,:,:,2) / dpu
          end if

          if (nsv>3) then
            svp(il:iu,jl:ju,kl:ku,4) = svp(il:iu,jl:ju,kl:ku,4) - &
                                         (Qpu/Apu) * 0.65 * inpu(:,:,:,4) / dpu
          end if

          deallocate(inpu)

       end if

  end do ! npurif

  end subroutine purifiers

end module modpurifiers
