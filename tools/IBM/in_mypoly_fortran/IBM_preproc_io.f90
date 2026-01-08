!> \file IBM_preproc_io.f90
!! Module for input/output routines for IBM preprocessor
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

module IBM_preproc_io
    use omp_lib
    implicit none

    contains

    subroutine num2logical(lvar, ivar, varname)
        implicit none
        integer,          intent(in)  :: ivar
        character(len=*), intent(in)  :: varname
        logical,          intent(out) :: lvar
        if (ivar==1) then
            lvar = .true.
        elseif (ivar==0) then
            lvar = .false.
        else
            print *, "Incorrect input for '" // varname // "'."
        end if
    end subroutine num2logical


    !! Read STL data from text files vertices_file and facets_file, and z grid
    subroutine read_data(vertices_file,n_vert,facets_file,n_fcts,zf_file,zh_file,ktot, &
        vertices,facets,incenters,faceNormals,zf,zh)
        implicit none

        integer, intent(in) :: n_vert, n_fcts, ktot
        character*9 , intent(in) :: facets_file
        character*12, intent(in) :: vertices_file 
        character*10, intent(in) :: zf_file, zh_file
        real, dimension(n_vert*3), intent(out) :: vertices
        integer, dimension(n_fcts*3), intent(out) :: facets
        real, dimension(n_fcts*3), intent(out) :: incenters, faceNormals
        real, intent(out) :: zf(ktot), zh(ktot)

        integer :: i, kk

        !$ call OMP_SET_NUM_THREADS(4)
        !$OMP parallel
        !$OMP sections private(i,kk)

        !$OMP section
        open(unit=1,file=vertices_file)
        do i=1,n_vert
        read(unit=1,fmt='(f15.10,x,f15.10,x,f15.10)') vertices(3*i-2), vertices(3*i-1), vertices(3*i)
        end do
        close(unit=1)

        !$OMP section
        15 format(i8,x,i8,x,i8,x,f15.10,x,f15.10,x,f15.10,x,f15.10,x,f15.10,x,f15.10)
        open(unit=2,file=facets_file)
        kk=1
        do i=1,n_fcts
        read(unit=2,fmt=15) facets(kk), facets(kk+1), facets(kk+2), &
                    incenters(kk), incenters(kk+1), incenters(kk+2), &
                    faceNormals(kk), faceNormals(kk+1), faceNormals(kk+2)
        kk = kk+3
        end do
        close(unit=2)

        !$OMP section
        open(unit=3,file=zf_file)
        do i = 1,ktot    
        read(unit=3,fmt='(f15.10)') zf(i)
        end do
        close(unit=3)

        !$OMP section
        open(unit=4,file=zh_file)
        do i = 1,ktot   
        read(unit=4,fmt='(f15.10)') zh(i)
        end do
        close(unit=4)

        !$OMP end sections
        !$OMP end parallel

    end subroutine read_data



    !! Subroutine to print STL data
    subroutine print_stl_data(n_vert,vertices,n_fcts,facets,incenters,faceNormals)
        implicit none

        integer, intent(in) :: n_vert, n_fcts
        real, dimension(n_vert*3), intent(in) :: vertices
        integer, dimension(n_fcts*3), intent(in) :: facets
        real, dimension(n_fcts*3), intent(in) :: incenters, faceNormals

        integer :: i, kk

        10 format(f15.10,x,f15.10,x,f15.10)

        kk=1
        print *,"VERTICES"
        do i=1,n_vert
            write(*,fmt=10,advance='no') vertices(kk), vertices(kk+1), vertices(kk+2)
            write(*,*)
            kk=kk+3
        end do

        kk=1
        print *,"FACETS"
        do i=1,n_fcts
            write(*,fmt='(i5,x,i5,x,i5)',advance='no') facets(kk), facets(kk+1), facets(kk+2)
            write(*,*)
            kk=kk+3
        end do

        kk=1
        print *,"IN-CENTERS"
        do i=1,n_fcts
            write(*,fmt=10,advance='no') incenters(kk), incenters(kk+1), incenters(kk+2)
            write(*,*)
            kk=kk+3
        end do

        kk=1
        print *,"FACE NORMALS"
        do i=1,n_fcts
            write(*,fmt=10,advance='no') faceNormals(kk), faceNormals(kk+1), faceNormals(kk+2)
            write(*,*)
            kk=kk+3
        end do

    end subroutine print_stl_data



    subroutine print_solid_flags(itot,jtot,ktot,solid_u,solid_v,solid_w,solid_c)
        implicit none

        integer, intent(in) :: itot, jtot, ktot
        logical, dimension(itot,jtot,ktot), intent(in) :: solid_u, solid_v, solid_w, solid_c

        integer :: ix, iy, iz

        !$ call OMP_SET_NUM_THREADS(4)
        !$OMP parallel
            !$OMP sections private(ix,iy,iz)
                
                !$OMP section
                open(unit=1,file='flag_u.txt')
                do iy = 1,jtot
                    do iz = 1,ktot
                        do ix = 1,itot
                            if (solid_u(ix,iy,iz)) then
                                    write(unit=1,fmt='(i1,A)',advance='no') 1, NEW_LINE('a')
                            else
                                    write(unit=1,fmt='(i1,A)',advance='no') 0, NEW_LINE('a')
                            end if
                        end do
                    end do
                end do
                close(unit=1)
                
                !$OMP section
                open(unit=2,file='flag_v.txt')
                do iy = 1,jtot
                    do iz = 1,ktot
                        do ix = 1,itot
                            if (solid_v(ix,iy,iz)) then
                                    write(unit=2,fmt='(i1,A)',advance='no') 1, NEW_LINE('a')
                            else
                                    write(unit=2,fmt='(i1,A)',advance='no') 0, NEW_LINE('a')
                            end if
                        end do
                    end do
                end do
                close(unit=2)
                
                !$OMP section
                open(unit=3,file='flag_w.txt')
                do iy = 1,jtot
                    do iz = 1,ktot
                        do ix = 1,itot
                            if (solid_w(ix,iy,iz)) then
                                    write(unit=3,fmt='(i1,A)',advance='no') 1, NEW_LINE('a')
                            else
                                    write(unit=3,fmt='(i1,A)',advance='no') 0, NEW_LINE('a')
                            end if
                        end do
                    end do
                end do
                close(unit=3)
                
                !$OMP section
                open(unit=4,file='flag_c.txt')
                do iy = 1,jtot
                    do iz = 1,ktot
                        do ix = 1,itot
                            if (solid_c(ix,iy,iz)) then
                                    write(unit=4,fmt='(i1,A)',advance='no') 1, NEW_LINE('a')
                            else
                                    write(unit=4,fmt='(i1,A)',advance='no') 0, NEW_LINE('a')
                            end if
                        end do
                    end do
                end do
                close(unit=4)

            !$OMP end sections
        !$OMP end parallel

    end subroutine print_solid_flags


    subroutine print_solid_points_index(itot, jtot, ktot, solid_u, solid_v, solid_w, solid_c, &
                                        filename_u, filename_v, filename_w, filename_c)
        implicit none
        integer,                            intent(in) :: itot, jtot, ktot
        logical, dimension(itot,jtot,ktot), intent(in) :: solid_u, solid_v, solid_w, solid_c
        character(len=*),                   intent(in) :: filename_u, filename_v, filename_w, filename_c

        !$ call OMP_SET_NUM_THREADS(4)
        !$OMP parallel
            !$OMP sections

                !$OMP section
                call write_position_index(itot, jtot, ktot, solid_u, filename_u, 1)

                !$OMP section
                call write_position_index(itot, jtot, ktot, solid_v, filename_v, 2)

                !$OMP section
                call write_position_index(itot, jtot, ktot, solid_w, filename_w, 3)

                !$OMP section
                call write_position_index(itot, jtot, ktot, solid_c, filename_c, 4)

            !$OMP end sections
        !$OMP end parallel

    end subroutine print_solid_points_index


    subroutine print_IB_index(itot, jtot, ktot, fluid_IB_u, fluid_IB_v, fluid_IB_w, fluid_IB_c, &
                                   filename_u, filename_v, filename_w, filename_c)
        implicit none
        integer,                            intent(in) :: itot, jtot, ktot
        logical, dimension(itot,jtot,ktot), intent(in) :: fluid_IB_u, fluid_IB_v, fluid_IB_w, fluid_IB_c
        character(len=*),                   intent(in) :: filename_u, filename_v, filename_w, filename_c

        !$ call OMP_SET_NUM_THREADS(4)
        !$OMP parallel
            !$OMP sections

                !$OMP section
                call write_position_index(itot, jtot, ktot, fluid_IB_u, filename_u, 1)

                !$OMP section
                call write_position_index(itot, jtot, ktot, fluid_IB_v, filename_v, 2)

                !$OMP section
                call write_position_index(itot, jtot, ktot, fluid_IB_w, filename_w, 3)

                !$OMP section
                call write_position_index(itot, jtot, ktot, fluid_IB_c, filename_c, 4)

            !$OMP end sections
        !$OMP end parallel

    end subroutine print_IB_index


    subroutine write_position_index(itot, jtot, ktot, field, filename, fid)
        implicit none
        integer,                            intent(in) :: itot, jtot, ktot
        logical, dimension(itot,jtot,ktot), intent(in) :: field
        character(len=*),                   intent(in) :: filename
        integer,                            intent(in) :: fid

        integer :: ix, iy, iz

        open(unit=fid,file=filename)
        write(unit=fid,fmt='(a18)') '# position (i,j,k)'
        do iy = 1,jtot
            do iz = 1,ktot
                do ix = 1,itot
                    if (field(ix,iy,iz)) then
                            write(unit=fid,fmt='(i4,x,i4,x,i4,A)',advance='no') ix, iy, iz, NEW_LINE('a')
                    end if
                end do
            end do
        end do
        close(unit=fid)
        ! write(*,*) 'Written ' // filename
    end subroutine write_position_index

    
end module IBM_preproc_io
