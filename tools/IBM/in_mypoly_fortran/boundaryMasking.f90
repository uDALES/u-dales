!> \file boundaryMasking.f90
!! Module to generate boundary masks for immersed boundary method
!>
!
!! \author Dipanjan Majumdar, ICL (2023-2026)
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

module boundaryMasking
  implicit none

  contains

    subroutine boundaryMasks(fluid_IB, solid_IB, fluid_IB_xyz, nfluid_IB, &
                             uvwc, itot, jtot, ktot, xgrid, ygrid, zgrid, &
                             solid_x, include_diagonals, stl_ground)
      implicit none
      character,                           intent(in)  :: uvwc
      integer,                             intent(in)  :: itot, jtot, ktot
      real,                                intent(in)  :: xgrid(itot), ygrid(jtot), zgrid(ktot)
      logical,  dimension(itot,jtot,ktot), intent(in)  :: solid_x
      logical,                             intent(in)  :: include_diagonals, stl_ground
      
      logical,  dimension(itot,jtot,ktot), intent(out) :: fluid_IB, solid_IB
      real,     allocatable,               intent(out) :: fluid_IB_xyz(:,:)
      integer,                             intent(out) :: nfluid_IB

      logical,  allocatable, dimension(:,:,:)          :: fluid_x, solid_x_w
      integer                                          :: ix, iy, iz, counter
      ! character*14                                     :: filename

      allocate(fluid_x(itot,jtot,ktot))

      do iy = 1,jtot
        do iz = 1,ktot
          do ix = 1,itot
            fluid_x(ix,iy,iz) = .not.(solid_x(ix,iy,iz))
          end do
        end do
      end do

      if (uvwc=='w') then

        allocate(solid_x_w(itot,jtot,ktot))
        solid_x_w = solid_x
        if (stl_ground) then
          do iy = 1,jtot
            do ix = 1,itot
              solid_x_w(ix,iy,1) = .true. !! % Bottom is always solid for w
            end do
          end do
        end if

        call getBoundaryCells(fluid_IB, solid_IB, itot, jtot, ktot, fluid_x, solid_x_w, include_diagonals)

        do iy = 1,jtot
          do ix = 1,itot
            fluid_IB(ix,iy,1) = .false. !! % Bottom is always solid for w
          end do
        end do
        deallocate(solid_x_w)

      else        !! For u, v and c

        call getBoundaryCells(fluid_IB, solid_IB, itot, jtot, ktot, fluid_x, solid_x, include_diagonals)

        if (stl_ground) then
          do iy = 1,jtot
            do ix = 1,itot
              if (.not.(solid_x(ix,iy,1))) then
                fluid_IB(ix,iy,1) = .true.
              end if
            end do
          end do
        end if

      end if

      nfluid_IB = count(fluid_IB)

      allocate(fluid_IB_xyz(nfluid_IB,3))
      counter = 1
      do iy = 1,jtot
        do iz = 1,ktot
          do ix = 1,itot
            if (fluid_IB(ix,iy,iz)) then
              fluid_IB_xyz(counter,:) = (/xgrid(ix), ygrid(iy), zgrid(iz)/)
              counter = counter + 1
            end if
          end do
        end do
      end do
        
      ! !$ call OMP_SET_NUM_THREADS(2)
      ! !$OMP parallel
      !     !$OMP sections private(ix,iy,iz,filename)
          
      !     !$OMP section
      !     filename = 'fluid_IB_' // uvwc // '.txt'
      !     open(unit=3,file=filename)
      !     do iy = 1,jtot
      !         do iz = 1,ktot
      !             do ix = 1,itot
      !                 if (fluid_IB(ix,iy,iz)) then
      !                         write(unit=3,fmt='(i1,A)',advance='no') 1, NEW_LINE('a')
      !                 else
      !                         write(unit=3,fmt='(i1,A)',advance='no') 0, NEW_LINE('a')
      !                 end if
      !             end do
      !         end do
      !     end do
      !     close(unit=3)

      !     !$OMP section
      !     filename = 'solid_IB_' // uvwc // '.txt'
      !     open(unit=4,file=filename)
      !     do iy = 1,jtot
      !         do iz = 1,ktot
      !             do ix = 1,itot
      !                 if (solid_IB(ix,iy,iz)) then
      !                         write(unit=4,fmt='(i1,A)',advance='no') 1, NEW_LINE('a')
      !                 else
      !                         write(unit=4,fmt='(i1,A)',advance='no') 0, NEW_LINE('a')
      !                 end if
      !             end do
      !         end do
      !     end do
      !     close(unit=4)

      !     !$OMP end sections
      ! !$OMP end parallel

      deallocate(fluid_x)
    end subroutine boundaryMasks


    subroutine getBoundaryCells(fluid_IB, solid_IB, itot, jtot, ktot, fluid, solid, include_diagonals)
      implicit none

      integer,                            intent(in)  :: itot, jtot, ktot
      logical, dimension(itot,jtot,ktot), intent(in)  :: fluid, solid
      logical,                            intent(in)  :: include_diagonals

      logical, dimension(itot,jtot,ktot), intent(out) :: fluid_IB
      logical, dimension(itot,jtot,ktot), intent(out) :: solid_IB

      integer :: ix, iy, iz

      ! Initializing with false
      do iy = 1,jtot
        do iz = 1,ktot
          do ix = 1,itot
            fluid_IB(ix,iy,iz) = .false.
            solid_IB(ix,iy,iz) = .false.
          end do
        end do
      end do

      do iy = 1,jtot
        do iz = 1,ktot
          do ix = 1,itot
            if (fluid(ix,iy,iz)) then

              ! Identify fluid IB points
              if (ix/=1) then
                if (solid(ix-1,iy,iz)) then
                  fluid_IB(ix,iy,iz) = .true.
                end if
              end if
              if (ix/=itot) then
                if (solid(ix+1,iy,iz)) then
                  fluid_IB(ix,iy,iz) = .true.
                end if
              end if
              if (iy/=1) then
                if (solid(ix,iy-1,iz)) then
                  fluid_IB(ix,iy,iz) = .true.
                end if
              end if
              if (iy/=jtot) then
                if (solid(ix,iy+1,iz)) then
                  fluid_IB(ix,iy,iz) = .true.
                end if
              end if
              if (iz/=1) then
                if (solid(ix,iy,iz-1)) then
                  fluid_IB(ix,iy,iz) = .true.
                end if
              end if
              if (iz/=ktot) then
                if (solid(ix,iy,iz+1)) then
                  fluid_IB(ix,iy,iz) = .true.
                end if
              end if

              if (include_diagonals) then

                if (ix/=1 .and. iy/=1) then
                  if (solid(ix-1,iy-1,iz)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=1 .and. iy/=jtot) then
                  if (solid(ix-1,iy+1,iz)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=itot .and. iy/=1) then
                  if (solid(ix+1,iy-1,iz)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=itot .and. iy/=jtot) then
                  if (solid(ix+1,iy+1,iz)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if

                if (ix/=1 .and. iz/=1) then
                  if (solid(ix-1,iy,iz-1)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=1 .and. iz/=ktot) then
                  if (solid(ix-1,iy,iz+1)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=itot .and. iz/=1) then
                  if (solid(ix+1,iy,iz-1)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=itot .and. iz/=ktot) then
                  if (solid(ix+1,iy,iz+1)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if

                if (iy/=1 .and. iz/=1) then
                  if (solid(ix,iy-1,iz-1)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (iy/=1 .and. iz/=ktot) then
                  if (solid(ix,iy-1,iz+1)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (iy/=jtot .and. iz/=1) then
                  if (solid(ix,iy+1,iz-1)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (iy/=jtot .and. iz/=ktot) then
                  if (solid(ix,iy+1,iz+1)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if

                if (ix/=1 .and. iy/=1 .and. iz/=1) then
                  if (solid(ix-1,iy-1,iz-1)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=itot .and. iy/=1 .and. iz/=1) then
                  if (solid(ix+1,iy-1,iz-1)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=1 .and. iy/=jtot .and. iz/=1) then
                  if (solid(ix-1,iy+1,iz-1)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=itot .and. iy/=jtot .and. iz/=1) then
                  if (solid(ix+1,iy+1,iz-1)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=1 .and. iy/=1 .and. iz/=ktot) then
                  if (solid(ix-1,iy-1,iz+1)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=itot .and. iy/=1 .and. iz/=ktot) then
                  if (solid(ix+1,iy-1,iz+1)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=1 .and. iy/=jtot .and. iz/=ktot) then
                  if (solid(ix-1,iy+1,iz+1)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=itot .and. iy/=jtot .and. iz/=ktot) then
                  if (solid(ix+1,iy+1,iz+1)) then
                    fluid_IB(ix,iy,iz) = .true.
                  end if
                end if

              end if

            else
              ! Identify solid IB points
              if (ix/=1) then
                if (fluid(ix-1,iy,iz)) then
                  solid_IB(ix,iy,iz) = .true.
                end if
              end if
              if (ix/=itot) then
                if (fluid(ix+1,iy,iz)) then
                  solid_IB(ix,iy,iz) = .true.
                end if
              end if
              if (iy/=1) then
                if (fluid(ix,iy-1,iz)) then
                  solid_IB(ix,iy,iz) = .true.
                end if
              end if
              if (iy/=jtot) then
                if (fluid(ix,iy+1,iz)) then
                  solid_IB(ix,iy,iz) = .true.
                end if
              end if
              if (iz/=1) then
                if (fluid(ix,iy,iz-1)) then
                  solid_IB(ix,iy,iz) = .true.
                end if
              end if
              if (iz/=ktot) then
                if (fluid(ix,iy,iz+1)) then
                  solid_IB(ix,iy,iz) = .true.
                end if
              end if

              if (include_diagonals) then

                if (ix/=1 .and. iy/=1) then
                  if (fluid(ix-1,iy-1,iz)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=1 .and. iy/=jtot) then
                  if (fluid(ix-1,iy+1,iz)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=itot .and. iy/=1) then
                  if (fluid(ix+1,iy-1,iz)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=itot .and. iy/=jtot) then
                  if (fluid(ix+1,iy+1,iz)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if

                if (ix/=1 .and. iz/=1) then
                  if (fluid(ix-1,iy,iz-1)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=1 .and. iz/=ktot) then
                  if (fluid(ix-1,iy,iz+1)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=itot .and. iz/=1) then
                  if (fluid(ix+1,iy,iz-1)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=itot .and. iz/=ktot) then
                  if (fluid(ix+1,iy,iz+1)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if

                if (iy/=1 .and. iz/=1) then
                  if (fluid(ix,iy-1,iz-1)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (iy/=1 .and. iz/=ktot) then
                  if (fluid(ix,iy-1,iz+1)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (iy/=jtot .and. iz/=1) then
                  if (fluid(ix,iy+1,iz-1)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (iy/=jtot .and. iz/=ktot) then
                  if (fluid(ix,iy+1,iz+1)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if

                if (ix/=1 .and. iy/=1 .and. iz/=1) then
                  if (fluid(ix-1,iy-1,iz-1)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=itot .and. iy/=1 .and. iz/=1) then
                  if (fluid(ix+1,iy-1,iz-1)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=1 .and. iy/=jtot .and. iz/=1) then
                  if (fluid(ix-1,iy+1,iz-1)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=itot .and. iy/=jtot .and. iz/=1) then
                  if (fluid(ix+1,iy+1,iz-1)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=1 .and. iy/=1 .and. iz/=ktot) then
                  if (fluid(ix-1,iy-1,iz+1)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=itot .and. iy/=1 .and. iz/=ktot) then
                  if (fluid(ix+1,iy-1,iz+1)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=1 .and. iy/=jtot .and. iz/=ktot) then
                  if (fluid(ix-1,iy+1,iz+1)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if
                if (ix/=itot .and. iy/=jtot .and. iz/=ktot) then
                  if (fluid(ix+1,iy+1,iz+1)) then
                    solid_IB(ix,iy,iz) = .true.
                  end if
                end if

              end if

            end if
          end do
        end do
      end do

    end subroutine getBoundaryCells

end module boundaryMasking