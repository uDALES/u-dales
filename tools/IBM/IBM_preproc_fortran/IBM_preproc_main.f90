!> \file IBM_preproc_main.f90
!! Main program for IBM preprocessor in Fortran
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

program IBM_preproc
    use omp_lib
    use ibm_preproc_mod
    implicit none

    integer                              :: n_vert, n_fcts
    real   , allocatable, dimension(:)   :: vertices_flat, faceNormals_flat, incenters_flat
    integer, allocatable, dimension(:)   :: facets_flat
    real   , allocatable, dimension(:,:) :: vertices, faceNormals, incenters
    integer, allocatable, dimension(:,:) :: facets
    real                                 :: dx, dy
    integer                              :: itot, jtot, ktot
    real   , allocatable, dimension(:)   :: zt, zm
    real   ,              dimension(3)   :: Ray_dir_u, Ray_dir_v, Ray_dir_w, Ray_dir_c
    integer                              :: stl_ground_i, diag_neighbs_i, periodic_x_i, periodic_y_i
    logical                              :: stl_ground, diag_neighbs, periodic_x, periodic_y
    integer                              :: n_threads
    real                                 :: tol

    integer,              dimension(13)  :: counts

    integer                              :: i

    open(unit=50,file='inmypoly_inp_info.txt')
        read(unit=50,fmt='(f15.10,x,f15.10)') dx, dy
        read(unit=50,fmt='(i5,x,i5,x,i5)') itot, jtot, ktot
        read(unit=50,fmt='(f15.10)') tol
        read(unit=50,fmt='(f15.10,x,f15.10,x,f15.10)') Ray_dir_u(1), Ray_dir_u(2), Ray_dir_u(3)
        read(unit=50,fmt='(f15.10,x,f15.10,x,f15.10)') Ray_dir_v(1), Ray_dir_v(2), Ray_dir_v(3)
        read(unit=50,fmt='(f15.10,x,f15.10,x,f15.10)') Ray_dir_w(1), Ray_dir_w(2), Ray_dir_w(3)
        read(unit=50,fmt='(f15.10,x,f15.10,x,f15.10)') Ray_dir_c(1), Ray_dir_c(2), Ray_dir_c(3)
        read(unit=50,fmt='(i8,x,i8)') n_vert, n_fcts
        read(unit=50,fmt='(i4)') n_threads
        read(unit=50,fmt='(i1,x,i1,x,i1,x,i1)') stl_ground_i, diag_neighbs_i, periodic_x_i, periodic_y_i
    close(unit=50)

    call num2logical(stl_ground  , stl_ground_i  , 'stl_ground'  )
    call num2logical(diag_neighbs, diag_neighbs_i, 'diag_neighbs')
    call num2logical(periodic_x  , periodic_x_i  , 'periodic_x'  )
    call num2logical(periodic_y  , periodic_y_i  , 'periodic_y'  )

    allocate(zt(ktot), zm(ktot))
    allocate(vertices_flat(n_vert*3), facets_flat(n_fcts*3), &
             incenters_flat(n_fcts*3), faceNormals_flat(n_fcts*3))
    allocate(vertices(n_vert,3), facets(n_fcts,3), &
             incenters(n_fcts,3), faceNormals(n_fcts,3))

    call read_data('vertices.txt', n_vert, 'faces.txt', n_fcts, 'zfgrid.txt', 'zhgrid.txt', ktot, &
                    vertices_flat, facets_flat, incenters_flat, faceNormals_flat, zt, zm)

    !$ call OMP_SET_NUM_THREADS(n_threads)
    !$OMP parallel do default(shared) private(i) schedule(static)
    do i = 1, n_vert
        vertices(i,1) = vertices_flat(3*i-2)
        vertices(i,2) = vertices_flat(3*i-1)
        vertices(i,3) = vertices_flat(3*i)
    end do
    !$OMP end parallel do

    !$OMP parallel do default(shared) private(i) schedule(static)
    do i = 1, n_fcts
        facets(i,1) = facets_flat(3*i-2)
        facets(i,2) = facets_flat(3*i-1)
        facets(i,3) = facets_flat(3*i)
        incenters(i,1) = incenters_flat(3*i-2)
        incenters(i,2) = incenters_flat(3*i-1)
        incenters(i,3) = incenters_flat(3*i)
        faceNormals(i,1) = faceNormals_flat(3*i-2)
        faceNormals(i,2) = faceNormals_flat(3*i-1)
        faceNormals(i,3) = faceNormals_flat(3*i)
    end do
    !$OMP end parallel do

    call run_ibm_preproc_core(n_vert, n_fcts, vertices, facets, incenters, faceNormals, &
                                  dx, dy, itot, jtot, ktot, zt, zm, &
                                  Ray_dir_u, Ray_dir_v, Ray_dir_w, Ray_dir_c, &
                                  stl_ground, diag_neighbs, periodic_x, periodic_y, &
                                  n_threads, tol, counts)

    open(unit=100,file='info_fort.txt')
        write(unit=100,fmt='(A)') 'nfcts  nsolpts_u  nsolpts_v  nsolpts_w  nsolpts_c  ' // &
                                'nbndpts_u  nbndpts_v  nbndpts_w  nbndpts_c  ' // &
                                'nfctsecs_u  nfctsecs_v  nfctsecs_w  nfctsecs_c'
        write(unit=100,fmt='(i12,x,i12,x,i12,x,i12,x,i12,x,i12,x,i12,x,i12,x,i12,x,i12,x,i12,x,i12,x,i12)') counts
    close(unit=100)

    deallocate(zt, zm)
    deallocate(vertices_flat, facets_flat, incenters_flat, faceNormals_flat)
    deallocate(vertices, facets, incenters, faceNormals)
end program IBM_preproc
