!> \file matchFacetsCells.f90
!! Module to match immersed boundary facets to fluid cells in a
!! staggered grid arrangement (u, v, w, c).
!>
!
!! \author Sam O. Owens, ICL (2020-2024)
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

module matchFacets2Cells
   use omp_lib
   implicit none

   ! Derived type to hold thread-local results
   type :: thread_results_type
      real, dimension(:), allocatable :: areas
      integer, dimension(:), allocatable :: facet_ids
      integer, dimension(:), allocatable :: bnd_pts
      real, dimension(:), allocatable :: distances
   end type thread_results_type

   contains


   subroutine readGeometry(fname_faces, nFaces, fname_vertices, nVertices, &
                           connectivityList, faceNormal, vertices)
     character(9), intent(in) :: fname_faces
     character(12), intent(in) :: fname_vertices
     integer, intent(in) :: nFaces, nVertices
     integer, intent(out) :: connectivityList(nFaces,3)
     real   , intent(out) :: faceNormal(nFaces,3), vertices(nVertices,3)
     real :: incenter(nFaces,3) ! not currently used by this program
     integer, parameter :: ifinput = 1
     integer :: n

     open (ifinput, file=fname_faces)
     do n = 1,nFaces
       read (ifinput, *) connectivityList(n,1), connectivityList(n,2), connectivityList(n,3), &
                         incenter(n,1), incenter(n,2), incenter(n,3), &
                         faceNormal(n,1), faceNormal(n,2), faceNormal(n,3)
     end do
     close (ifinput)

     open (ifinput, file=fname_vertices)
     do n = 1,nVertices
       read (ifinput, *) vertices(n,1), vertices(n,2), vertices(n,3)
     end do
     close (ifinput)

   end subroutine readGeometry


   subroutine readBoundaryPoints(xgrid, ygrid, zgrid, itot, jtot, ktot, &
      fname_fluid_boundary, nfluid_IB, fname_solid_boundary, nsolid_IB, &
          fluid_IB,         solid_IB,        fluid_IB_xyz)

   character(20), intent(in) :: fname_fluid_boundary, fname_solid_boundary
   integer, intent(in) :: itot, jtot, ktot
   real, intent(in) :: xgrid(itot), ygrid(jtot), zgrid(ktot)
   integer, intent(in) :: nfluid_IB, nsolid_IB
   logical, dimension(itot,jtot,ktot), intent(out) :: fluid_IB, solid_IB
   integer, allocatable :: fluid_IB_ijk(:,:), solid_IB_ijk(:,:)
   real, allocatable, intent(out) :: fluid_IB_xyz(:,:)
   integer :: fluid_IB_read, solid_IB_read, i, j, k, n
   character(80) :: chmess
   integer, parameter :: ifinput = 1

   allocate(fluid_IB_ijk(nfluid_IB,3), solid_IB_ijk(nsolid_IB,3), fluid_IB_xyz(nfluid_IB,3))

   open (ifinput, file=fname_fluid_boundary)
   read (ifinput, '(a80)') chmess
   do n = 1,nfluid_IB
     read (ifinput, *) fluid_IB_ijk(n,1), fluid_IB_ijk(n,2), fluid_IB_ijk(n,3)
     fluid_IB_xyz(n,:) = (/xgrid(fluid_IB_ijk(n,1)), ygrid(fluid_IB_ijk(n,2)), zgrid(fluid_IB_ijk(n,3))/)
     fluid_IB(fluid_IB_ijk(n,1), fluid_IB_ijk(n,2), fluid_IB_ijk(n,3)) = .true.
   end do
   close (ifinput)

   open (ifinput, file=fname_solid_boundary)
   read (ifinput, '(a80)') chmess
   do n = 1,nsolid_IB
     read (ifinput, *) solid_IB_ijk(n,1), solid_IB_ijk(n,2), solid_IB_ijk(n,3)
     solid_IB(solid_IB_ijk(n,1), solid_IB_ijk(n,2), solid_IB_ijk(n,3)) = .true.
   end do
   close (ifinput)

   deallocate(fluid_IB_ijk, solid_IB_ijk)


   end subroutine readBoundaryPoints


   subroutine matchFacetsToCells(connectivityList, faceNormal, nFaces, vertices, nVertices, &
      fluid_IB, solid_IB, fluid_IB_xyz, nfluid_IB, &
      xgrid, ygrid, zgrid, itot, jtot, ktot, &
      diag_neighbs, periodic_x, periodic_y, n_threads, &
      secfacids, secbndptids, secareas, bnddst, nfacsecs)

    use, intrinsic :: ieee_arithmetic
    integer, intent(in) :: nFaces, nVertices, nfluid_IB, itot, jtot, ktot, n_threads
    integer, intent(in) :: connectivityList(nFaces,3)
    real   , intent(in) :: faceNormal(nFaces,3), vertices(nVertices,3), fluid_IB_xyz(nfluid_IB,3)
    real   , intent(in) :: xgrid(itot), ygrid(jtot), zgrid(ktot)
    logical, intent(in), dimension(itot,jtot,ktot) :: fluid_IB, solid_IB
    logical, intent(in) :: diag_neighbs, periodic_x, periodic_y
    integer, dimension(:), allocatable, intent(out) :: secfacids, secbndptids
    real, dimension(:), allocatable, intent(out) :: secareas, bnddst
    integer, intent(out) :: nfacsecs
    !$ type(thread_results_type), dimension(:), allocatable :: thread_data
    !$ integer :: actual_threads
    integer, dimension(:,:), allocatable :: clipFaces
    real, dimension(:,:), allocatable :: clipVertices, projVert, clipFaceNormal
    logical :: il_comp(itot+1), iu_comp(itot+1), jl_comp(jtot+1), ju_comp(jtot+1), kl_comp(ktot+1), ku_comp(ktot+1)
    logical :: search_adj
    integer :: il, iu, jl, ju, kl, ku, nClipFaces, nClipVertices, id, loc
    real :: dx, dy, dz, xmin, xmax, ymin, ymax, zmin, zmax, xl, xu, yl, yu, zl, zu, tol, &
            planes(6,4), area, area_miss, planeNormal(3), xproj, yproj, zproj, proj, projVec(3), projArea, &
            dist, dists(27), angle, angles(27), BI(3), BIs(27,3)
    real, dimension(3) :: xyz, xyz1, xyz2, xyz3, xyz4, xyz5, xyz6, xyz7, xyz8, xyz9, &
                          xyz10, xyz11, xyz12, xyz13, xyz14, xyz15, xyz16, xyz17, xyz18, &
                          xyz19, xyz20, xyz21, xyz22, xyz23, xyz24, xyz25, xyz26, xyz27
    integer :: n, m, i, j, k, p, dir, ids(2), thread_id

    allocate(secareas(0))
    allocate(secfacids(0))
    allocate(secbndptids(0))
    allocate(bnddst(0))

    !$ actual_threads = 1
    ! Get actual number of threads that will be used
    !$ call OMP_SET_NUM_THREADS(n_threads)
    !$OMP parallel
    !$ actual_threads = omp_get_num_threads()
    !$OMP end parallel

    ! Initialize thread-local data arrays
    !$ allocate(thread_data(actual_threads))
    !$ do i = 1, actual_threads
    !$     allocate(thread_data(i)%areas(0))
    !$     allocate(thread_data(i)%facet_ids(0))
    !$     allocate(thread_data(i)%bnd_pts(0))
    !$     allocate(thread_data(i)%distances(0))
    !$ end do

    dx = xgrid(2)-xgrid(1)
    dy = ygrid(2)-ygrid(1)
    dz = zgrid(2)-zgrid(1)

    tol = 1e-8 ! machine precision errors
    area_miss = 0.0 ! Initialize for OpenMP reduction

    !$OMP parallel do default(none) &
    !$OMP shared(nFaces, connectivityList, faceNormal, vertices, xgrid, ygrid, zgrid, tol, &
    !$OMP        fluid_IB, solid_IB, fluid_IB_xyz, nfluid_IB, itot, jtot, ktot, &
    !$OMP        diag_neighbs, periodic_x, periodic_y, dx, dy, dz, &
    !$OMP        thread_data, actual_threads) &
    !$OMP private(xmin, xmax, ymin, ymax, zmin, zmax, &
    !$OMP         il_comp, iu_comp, jl_comp, ju_comp, kl_comp, ku_comp, &
    !$OMP         il, iu, jl, ju, kl, ku, xl, xu, yl, yu, zl, zu, planes, &
    !$OMP         clipVertices, clipFaces, clipFaceNormal, nClipFaces, nClipVertices, &
    !$OMP         projVert, planeNormal, xproj, yproj, zproj, proj, projVec, projArea, &
    !$OMP         area, dist, dists, angle, angles, BI, BIs, &
    !$OMP         xyz, xyz1, xyz2, xyz3, xyz4, xyz5, xyz6, xyz7, xyz8, xyz9, &
    !$OMP         xyz10, xyz11, xyz12, xyz13, xyz14, xyz15, xyz16, xyz17, xyz18, &
    !$OMP         xyz19, xyz20, xyz21, xyz22, xyz23, xyz24, xyz25, xyz26, xyz27, &
    !$OMP         search_adj, id, loc, i, j, k, m, p, dir, ids, thread_id) &
    !$OMP reduction(+:area_miss) schedule(dynamic)
    do n=1,nFaces
      ! Get thread ID (1-based)
      !$ thread_id = omp_get_thread_num() + 1
      
      ! no shear stress in normal direction
      if ((xgrid(1) == 0. .and. all(abs(abs(faceNormal(n, :)) - (/1.,0.,0./)) < tol)) .or. &
          (ygrid(1) == 0. .and. all(abs(abs(faceNormal(n, :)) - (/0.,1.,0./)) < tol)) .or. &
          (zgrid(1) == 0. .and. all(abs(abs(faceNormal(n, :)) - (/0.,0.,1./)) < tol))) cycle

      xmin = minval(vertices(connectivityList(n,:),1))
      ymin = minval(vertices(connectivityList(n,:),2))
      zmin = minval(vertices(connectivityList(n,:),3))
      xmax = maxval(vertices(connectivityList(n,:),1))
      ymax = maxval(vertices(connectivityList(n,:),2))
      zmax = maxval(vertices(connectivityList(n,:),3))

      ! ignore facets on the ground and facing down
      if ((abs(zmin) < epsilon(zmin) .and. abs(zmax) < epsilon(zmax)) .and. &
       all(abs(faceNormal(n, :) - (/0.,0.,-1./)) < epsilon(faceNormal(n, 1)))) cycle

      where (xmin >= (/xgrid   -dx/2., xgrid(itot)+dx/2./)-tol)
         il_comp = .true.
      elsewhere
         il_comp = .false.
      end where
      where (xmax <= (/xgrid(1)-dx/2., xgrid      +dx/2./)+tol)
         iu_comp = .true.
      elsewhere
         iu_comp = .false.
      end where
      where (ymin >= (/ygrid   -dy/2., ygrid(jtot)+dy/2./)-tol)
         jl_comp = .true.
      elsewhere
         jl_comp = .false.
      end where
      where (ymax <= (/ygrid(1)-dy/2., ygrid      +dy/2./)+tol)
         ju_comp = .true.
      elsewhere
         ju_comp = .false.
      end where
      where (zmin >= (/zgrid   -dz/2., zgrid(ktot)+dz/2./)-tol)
         kl_comp = .true.
      elsewhere
         kl_comp = .false.
      end where
      where (zmax <= (/zgrid(1)-dz/2., zgrid      +dz/2./)+tol)
         ku_comp = .true.
      elsewhere
         ku_comp = .false.
      end where

      il = findloc(il_comp, .true., 1, back=.true.)
      iu = findloc(iu_comp, .true., 1)
      jl = findloc(jl_comp, .true., 1, back=.true.)
      ju = findloc(ju_comp, .true., 1)
      kl = findloc(kl_comp, .true., 1, back=.true.)
      ku = findloc(ku_comp, .true., 1)

      ! Not sure these are needed?
      if (xmax > xgrid(itot) + dx/2) iu = itot
      if (ymax > ygrid(jtot) + dy/2) ju = jtot

      if (il==0) then
         !write(*,*) "Warning: skipping facet ", n, " as it is out of bounds in lower x direction."
         cycle
      end if
      if (iu==0) then
         !write(*,*) "Warning: skipping facet ", n, " as it is out of bounds in upper x direction."
         cycle
      end if
      if (jl==0) then
         !write(*,*) "warning: skipping facet ", n, " as it is out of bounds in lower y direction."
         cycle
      end if
      if (ju==0) then
         !write(*,*) "warning: skipping facet ", n, " as it is out of bounds in upper y direction."
         cycle
      end if
      if (kl==0) then
         !write(*,*) "warning: skipping facet ", n, " as it is out of bounds in lower z direction."
         cycle
      end if
      if (ku==0) then
         !write(*,*) "warning: skipping facet ", n, " as it is out of bounds in upper z direction."
         cycle
      end if

      ! Facet exists in cell N+1.
      ! Freqently occurs for u and v grids, currently the solution is to
      ! double cell 1 areas in periodic cases (because cell N+1 = cell 1).
      ! This causes small errors due to the shape of facets.
      if (iu > itot) iu = itot
      if (ju > jtot) ju = jtot

      do i=il,iu
         do j=jl,ju
            do k=kl,ku
               if (.not.(fluid_IB(i,j,k) .or. solid_IB(i,j,k))) cycle
               ! Define corners of cube
               xl = xgrid(i) - dx/2. - tol
               xu = xgrid(i) + dx/2. + tol
               yl = ygrid(j) - dy/2. - tol
               yu = ygrid(j) + dy/2. + tol
               zl = zgrid(k) - dz/2. - tol
               zu = zgrid(k) + dz/2. + tol

               planes(1,:) = (/ 1., 0., 0., xu/)
               planes(2,:) = (/-1., 0., 0.,-xl/)
               planes(3,:) = (/ 0., 1., 0., yu/)
               planes(4,:) = (/ 0.,-1., 0.,-yl/)
               planes(5,:) = (/ 0., 0., 1., zu/)
               planes(6,:) = (/ 0., 0.,-1.,-zl/)

               call sutherlandHodgman3D(vertices(connectivityList(n,:),:), 3, planes, 6, clipVertices)

               nClipVertices = size(clipVertices, 1)

               if (nClipVertices < 3) then
                  !nClipFaces = 0
                  !allocate(clipFaces(nClipFaces,3))
                  !deallocate(clipVertices)
                  cycle
               elseif (nClipVertices == 3) then
                  area = 0.5*norm2(cross_product(clipVertices(2,:) - clipVertices(1,:), &
                     clipVertices(3,:) - clipVertices(1,:)))
                  ! remove anything below 1 square centimetre, as we only write to this precision.
                  if (area < 1e-5) cycle

                  nClipFaces = 1
                  allocate(clipFaces(nClipFaces,3))
                  clipFaces(1,:) = (/1, 2, 3/)

               elseif (nClipVertices > 3) then
                  xproj = dot_product(faceNormal(n,:), (/1., 0., 0./))
                  yproj = dot_product(faceNormal(n,:), (/0., 1., 0./))
                  zproj = dot_product(faceNormal(n,:), (/0., 0., 1./))
                  projvec = abs((/xproj, yproj, zproj/))
                  dir = maxloc(projvec, 1)
                  proj = projvec(dir)

                  if (dir==0) then
                     write(*,*) "something wrong with finding direction to project in"
                  elseif (dir==1) then
                     planeNormal = (/1., 0., 0./)
                     ids = (/2, 3/)
                  elseif (dir==2) then
                     planeNormal = (/0., 1., 0./)
                     ids = (/1, 3/)
                  elseif (dir==3) then
                     planeNormal = (/0., 0., 1./)
                     ids = (/1, 2/)
                  end if

                  allocate(projVert(nClipVertices,3))
                  do m = 1,nClipVertices
                     projVert(m, :) = clipVertices(m, :) - &
                        dot_product(clipVertices(m, :), planeNormal) * planeNormal
                  end do

                  projArea = abs(polyarea(projVert(:,ids(1)), projVert(:,ids(2)), nClipVertices))

                  deallocate(projVert)

                  area = projArea / proj

                  if (area < 1e-5) cycle

                  nClipFaces = nClipVertices - 2
                  allocate(clipFaces(nClipFaces,3))
                  do m=1,nClipFaces
                     clipFaces(m,:) = (/1, m+1, m+2/)
                  end do

               else
                 write(*,*) "something wrong with clipped polygon"
               end if

               if (((xgrid(i) == 0.) .and. periodic_x) .or. ((ygrid(j) == 0.) .and. periodic_y)) then ! Account for periodicity - flux at point N+1
                  area = area * 2.
               end if

               dists = ieee_value(dists, ieee_quiet_nan)
               angles = ieee_value(angles, ieee_quiet_nan)
               BIs = ieee_value(BIs, ieee_quiet_nan)
               search_adj = .false.

               allocate(clipFaceNormal(nClipFaces,3))
               do m=1,nClipFaces
                  clipFaceNormal(m, :) = faceNormal(n, :)
               end do

               if (fluid_IB(i,j,k)) then
                  xyz1 = (/xgrid(i), ygrid(j), zgrid(k)/)
                  loc = findloc(ismember_rows(fluid_IB_xyz, xyz1), .true., 1)
                  !facet_section(3) = loc;
                  call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                     xyz1, dist, BI)
                  angle = dot_product(faceNormal(n,:), (xyz1 - BI)/norm2(xyz1 - BI))

                  if (abs(angle - 1.) < epsilon(angle)) then ! Wall-normal defined, use this cell
                  !if (abs(angle - 1.) < epsilon(angle) .and. dist > 0.05*exp(1.)) then ! Wall-normal defined, use this cell
                     id = 1 ! not necessary?
                     xyz = xyz1
                     !write(*,*) "normal found"
                  else
                     !write(*,*) "normal not found, searching adjacent cells"
                     ! Not normal, search adjacent fluid IB cells
                     search_adj = .true.
                     if (dist > 0) then
                     !if (dist > 0.05*exp(1.)) then
                        ! Include in comparison
                        dists(1) = dist
                        angles(1) = angle
                        BIs(1,:) = BI
                     end if
                  end if
               end if

               if (solid_IB(i,j,k) .or. search_adj) then

                   if (i /= 1) then
                      if (fluid_IB(i-1,j,k)) then
                        xyz2 = (/xgrid(i-1), ygrid(j), zgrid(k)/)
                        call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                           xyz2, dist, BI)
                        dists(2) = dist
                        angles(2) = dot_product(faceNormal(n,:), (xyz2 - BI)/norm2(xyz2 - BI))
                        BIs(2,:) = BI
                     end if !(fluid_IB(i-1,j,k))
                  end if !(i /= 1)

                  if (i /= itot) then
                     if (fluid_IB(i+1,j,k)) then
                        xyz3 = (/xgrid(i+1), ygrid(j), zgrid(k)/)
                        call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                           xyz3, dist, BI)
                        dists(3) = dist
                        angles(3) = dot_product(faceNormal(n,:), (xyz3 - BI)/norm2(xyz3 - BI))
                        BIs(3,:) = BI
                     end if
                  end if

                  if (j /= 1) then
                     if (fluid_IB(i,j-1,k)) then
                        xyz4 = (/xgrid(i), ygrid(j-1), zgrid(k)/)
                        call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                           xyz4, dist, BI)
                        dists(4) = dist
                        angles(4) = dot_product(faceNormal(n,:), (xyz4 - BI)/norm2(xyz4 - BI))
                        BIs(4,:) = BI
                     end if
                  end if

                  if (j /= jtot) then
                     if (fluid_IB(i,j+1,k)) then
                        xyz5 = (/xgrid(i), ygrid(j+1), zgrid(k)/)
                        call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                           xyz5, dist, BI)
                        dists(5) = dist
                        angles(5) = dot_product(faceNormal(n,:), (xyz5 - BI)/norm2(xyz5 - BI))
                        BIs(5,:) = BI
                     end if
                  end if

                  if (k /= 1) then
                     if (fluid_IB(i,j,k-1)) then
                        xyz6 = (/xgrid(i), ygrid(j), zgrid(k-1)/)
                        call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                           xyz6, dist, BI)
                        dists(6) = dist
                        angles(6) = dot_product(faceNormal(n,:), (xyz6 - BI)/norm2(xyz6 - BI))
                        BIs(6,:) = BI
                     end if
                  end if

                  if (k /= ktot) then
                     if (fluid_IB(i,j,k+1)) then
                        xyz7 = (/xgrid(i), ygrid(j), zgrid(k+1)/)
                        call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                           xyz7, dist, BI)
                        dists(7) = dist
                        angles(7) = dot_product(faceNormal(n,:), (xyz7 - BI)/norm2(xyz7 - BI))
                        BIs(7,:) = BI
                     end if
                  end if

                  if (diag_neighbs) then
                     if (i/=1 .and. j/=1) then
                        if (fluid_IB(i-1,j-1,k)) then
                           xyz8 = (/xgrid(i-1), ygrid(j-1), zgrid(k)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz8, dist, BI)
                           dists(8) = dist
                           angles(8) = dot_product(faceNormal(n,:), (xyz8 - BI)/norm2(xyz8 - BI))
                           BIs(8,:) = BI
                        end if
                     end if

                     if (i/=1 .and. j/=jtot) then
                        if (fluid_IB(i-1,j+1,k)) then
                           xyz9 = (/xgrid(i-1), ygrid(j+1), zgrid(k)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz9, dist, BI)
                           dists(9) = dist
                           angles(9) = dot_product(faceNormal(n,:), (xyz9 - BI)/norm2(xyz9 - BI))
                           BIs(9,:) = BI
                        end if
                     end if

                     if (i/=itot .and. j/=1) then
                        if (fluid_IB(i+1,j-1,k)) then
                           xyz10 = (/xgrid(i+1), ygrid(j-1), zgrid(k)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz10, dist, BI)
                           dists(10) = dist
                           angles(10) = dot_product(faceNormal(n,:), (xyz10 - BI)/norm2(xyz10 - BI))
                           BIs(10,:) = BI
                        end if
                     end if

                     if (i/=itot .and. j/=jtot) then
                        if (fluid_IB(i+1,j+1,k)) then
                           xyz11 = (/xgrid(i+1), ygrid(j+1), zgrid(k)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz11, dist, BI)
                           dists(11) = dist
                           angles(11) = dot_product(faceNormal(n,:), (xyz11 - BI)/norm2(xyz11 - BI))
                           BIs(11,:) = BI
                        end if
                     end if

                     if (i/=1 .and. k/=1) then
                        if (fluid_IB(i-1,j,k-1)) then
                           xyz12 = (/xgrid(i-1), ygrid(j), zgrid(k-1)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz12, dist, BI)
                           dists(12) = dist
                           angles(12) = dot_product(faceNormal(n,:), (xyz12 - BI)/norm2(xyz12 - BI))
                           BIs(12,:) = BI
                        end if
                     end if

                     if (i/=1 .and. k/=ktot) then
                        if (fluid_IB(i-1,j,k+1)) then
                           xyz13 = (/xgrid(i-1), ygrid(j), zgrid(k+1)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz13, dist, BI)
                           dists(13) = dist
                           angles(13) = dot_product(faceNormal(n,:), (xyz13 - BI)/norm2(xyz13 - BI))
                           BIs(13,:) = BI
                        end if
                     end if

                     if (i/=itot .and. k/=1) then
                        if (fluid_IB(i+1,j,k-1)) then
                           xyz14 = (/xgrid(i+1), ygrid(j), zgrid(k-1)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz14, dist, BI)
                           dists(14) = dist
                           angles(14) = dot_product(faceNormal(n,:), (xyz14 - BI)/norm2(xyz14 - BI))
                           BIs(14,:) = BI
                        end if
                     end if

                     if (i/=itot .and. k/=ktot) then
                        if (fluid_IB(i+1,j,k+1)) then
                           xyz15 = (/xgrid(i+1), ygrid(j), zgrid(k+1)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz15, dist, BI)
                           dists(15) = dist
                           angles(15) = dot_product(faceNormal(n,:), (xyz15 - BI)/norm2(xyz15 - BI))
                           BIs(15,:) = BI
                        end if
                     end if

                     if (j/=1 .and. k/=1) then
                        if (fluid_IB(i,j-1,k-1)) then
                           xyz16 = (/xgrid(i), ygrid(j-1), zgrid(k-1)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz16, dist, BI)
                           dists(16) = dist
                           angles(16) = dot_product(faceNormal(n,:), (xyz16 - BI)/norm2(xyz16 - BI))
                           BIs(16,:) = BI
                        end if
                     end if

                     if (j/=1 .and. k/=ktot) then
                        if (fluid_IB(i,j-1,k+1)) then
                           xyz17 = (/xgrid(i), ygrid(j-1), zgrid(k+1)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz17, dist, BI)
                           dists(17) = dist
                           angles(17) = dot_product(faceNormal(n,:), (xyz17 - BI)/norm2(xyz17 - BI))
                           BIs(17,:) = BI
                        end if
                     end if

                     if (j/=jtot .and. k/=1) then
                        if (fluid_IB(i,j+1,k-1)) then
                           xyz18 = (/xgrid(i), ygrid(j+1), zgrid(k-1)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz18, dist, BI)
                           dists(18) = dist
                           angles(18) = dot_product(faceNormal(n,:), (xyz18 - BI)/norm2(xyz18 - BI))
                           BIs(18,:) = BI
                        end if
                     end if

                     if (j/=jtot .and. k/=ktot) then
                        if (fluid_IB(i,j+1,k+1)) then
                           xyz19 = (/xgrid(i), ygrid(j+1), zgrid(k+1)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz19, dist, BI)
                           dists(19) = dist
                           angles(19) = dot_product(faceNormal(n,:), (xyz19 - BI)/norm2(xyz19 - BI))
                           BIs(19,:) = BI
                        end if
                     end if

                     if (i/=1 .and. j/=1 .and. k/=1) then
                        if (fluid_IB(i-1,j-1,k-1)) then
                           xyz20 = (/xgrid(i-1), ygrid(j-1), zgrid(k-1)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz20, dist, BI)
                           dists(20) = dist
                           angles(20) = dot_product(faceNormal(n,:), (xyz20 - BI)/norm2(xyz20 - BI))
                           BIs(20,:) = BI
                        end if
                     end if

                     if (i/=itot .and. j/=1 .and. k/=1) then
                        if (fluid_IB(i+1,j-1,k-1)) then
                           xyz21 = (/xgrid(i+1), ygrid(j-1), zgrid(k-1)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz21, dist, BI)
                           dists(21) = dist
                           angles(21) = dot_product(faceNormal(n,:), (xyz21 - BI)/norm2(xyz21 - BI))
                           BIs(21,:) = BI
                        end if
                     end if

                     if (i/=1 .and. j/=jtot .and. k/=1) then
                        if (fluid_IB(i-1,j+1,k-1)) then
                           xyz22 = (/xgrid(i-1), ygrid(j+1), zgrid(k-1)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz22, dist, BI)
                           dists(22) = dist
                           angles(22) = dot_product(faceNormal(n,:), (xyz22 - BI)/norm2(xyz22 - BI))
                           BIs(22,:) = BI
                        end if
                     end if

                     if (i/=itot .and. j/=jtot .and. k/=1) then
                        if (fluid_IB(i+1,j+1,k-1)) then
                           xyz23 = (/xgrid(i+1), ygrid(j+1), zgrid(k-1)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz23, dist, BI)
                           dists(23) = dist
                           angles(23) = dot_product(faceNormal(n,:), (xyz23 - BI)/norm2(xyz23 - BI))
                           BIs(23,:) = BI
                        end if
                     end if

                     if (i/=1 .and. j/=1 .and. k/=ktot) then
                        if (fluid_IB(i-1,j-1,k+1)) then
                           xyz24 = (/xgrid(i-1), ygrid(j-1), zgrid(k+1)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz24, dist, BI)
                           dists(24) = dist
                           angles(24) = dot_product(faceNormal(n,:), (xyz24 - BI)/norm2(xyz24 - BI))
                           BIs(24,:) = BI
                        end if
                     end if

                     if (i/=itot .and. j/=1 .and. k/=ktot) then
                        if (fluid_IB(i+1,j-1,k+1)) then
                           xyz25 = (/xgrid(i+1), ygrid(j-1), zgrid(k+1)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz25, dist, BI)
                           dists(25) = dist
                           angles(25) = dot_product(faceNormal(n,:), (xyz25 - BI)/norm2(xyz25 - BI))
                           BIs(25,:) = BI
                        end if
                     end if

                     if (i/=1 .and. j/=jtot .and. k/=ktot) then
                        if (fluid_IB(i-1,j+1,k+1)) then
                           xyz26 = (/xgrid(i-1), ygrid(j+1), zgrid(k+1)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz26, dist, BI)
                           dists(26) = dist
                           angles(26) = dot_product(faceNormal(n,:), (xyz26 - BI)/norm2(xyz26 - BI))
                           BIs(26,:) = BI
                        end if
                     end if

                     if (i/=itot .and. j/=jtot .and. k/=ktot) then
                        if (fluid_IB(i+1,j+1,k+1)) then
                           xyz27 = (/xgrid(i+1), ygrid(j+1), zgrid(k+1)/)
                           call fastPoint2TriMesh(clipFaces, clipFaceNormal, nClipFaces, clipVertices, nClipVertices, &
                              xyz27, dist, BI)
                           dists(27) = dist
                           angles(27) = dot_product(faceNormal(n,:), (xyz27 - BI)/norm2(xyz27 - BI))
                           BIs(27,:) = BI
                        end if
                     end if

                  end if ! diag_neighbs

                  ! do p = 1,27
                  !   if (dists(p) < 0.05*exp(1.)) then
                  !       dists(p) = ieee_value(dists(p), ieee_quiet_nan)
                  !    end if
                  ! end do

                  id = maxloc(abs(angles) / (dists / (dx*dy*dz)**(1./3.)), 1)
                  dist = dists(id)
                  BI = BIs(id,:)

                  if (isnan(dist)) then
                     !write(*,*) "facet ", n, " in cell ", i, j, k, " could not find a cell to give flux to"
                     area_miss = area_miss + area
                     deallocate(clipFaces)
                     deallocate(clipFaceNormal)
                     cycle
                  end if

                  select case (id)
                  case(1)
                     xyz = xyz1
                  case(2)
                     xyz = xyz2
                  case(3)
                     xyz = xyz3
                  case(4)
                     xyz = xyz4
                  case(5)
                     xyz = xyz5
                  case(6)
                     xyz = xyz6
                  case(7)
                     xyz = xyz7
                  case(8)
                     xyz = xyz8
                  case(9)
                     xyz = xyz9
                  case(10)
                     xyz = xyz10
                  case(11)
                     xyz = xyz11
                  case(12)
                     xyz = xyz12
                  case(13)
                     xyz = xyz13
                  case(14)
                     xyz = xyz14
                  case(15)
                     xyz = xyz15
                  case(16)
                     xyz = xyz16
                  case(17)
                     xyz = xyz17
                  case(18)
                     xyz = xyz18
                  case(19)
                     xyz = xyz19
                  case(20)
                     xyz = xyz20
                  case(21)
                     xyz = xyz21
                  case(22)
                     xyz = xyz22
                  case(23)
                     xyz = xyz23
                  case(24)
                     xyz = xyz24
                  case(25)
                     xyz = xyz25
                  case(26)
                     xyz = xyz26
                  case(27)
                     xyz = xyz27

                  end select

               end if !(solid_IB(i,j,k) .or. search_adj)

               if (isnan(dist) .or. abs(dist)<tol) dist = 0.1 ! 0.1 m minimum distance to avoid singularities
               if (isnan(area)) area = 0.0

               loc = findloc(ismember_rows(fluid_IB_xyz, xyz), .true., 1)

               !! For serial run version: Append to global arrays
               ! call appendToArray1D_real(secareas, area)
               ! call appendToArray1D_integer(secfacids, n)
               ! call appendToArray1D_integer(secbndptids, loc)
               ! call appendToArray1D_real(bnddst, abs(dist))
               
               ! Alternatively for OpenMP run: Append to thread-specific arrays
               !$ call appendToThreadResults(thread_data(thread_id), area, n, loc, abs(dist))

               deallocate(clipFaces)
               deallocate(clipFaceNormal)
            end do
         end do
      end do
    end do
    !$OMP end parallel do

    ! Merge and sort results from all threads
    !$ call mergeAndSortThreadResults(thread_data, actual_threads, secfacids, secareas, secbndptids, bnddst)

    nfacsecs = size(secfacids,1)

   write(*,*) "Total area missing flux: ", area_miss, " m^2"

end subroutine matchFacetsToCells


subroutine appendToThreadResults(thread_data, area, facet_id, bnd_pt, distance)
   ! Append results to thread's local arrays efficiently
   implicit none
   type(thread_results_type), intent(inout) :: thread_data
   real, intent(in) :: area, distance
   integer, intent(in) :: facet_id, bnd_pt
   
   real, dimension(:), allocatable :: temp_real
   integer, dimension(:), allocatable :: temp_int
   integer :: n
   
   ! Get current size
   n = size(thread_data%areas)
   
   ! Expand areas array
   allocate(temp_real(n+1))
   if (n > 0) temp_real(1:n) = thread_data%areas
   temp_real(n+1) = area
   call move_alloc(temp_real, thread_data%areas)
   
   ! Expand facet_ids array
   allocate(temp_int(n+1))
   if (n > 0) temp_int(1:n) = thread_data%facet_ids
   temp_int(n+1) = facet_id
   call move_alloc(temp_int, thread_data%facet_ids)
   
   ! Expand bnd_pts array
   allocate(temp_int(n+1))
   if (n > 0) temp_int(1:n) = thread_data%bnd_pts
   temp_int(n+1) = bnd_pt
   call move_alloc(temp_int, thread_data%bnd_pts)
   
   ! Expand distances array
   allocate(temp_real(n+1))
   if (n > 0) temp_real(1:n) = thread_data%distances
   temp_real(n+1) = distance
   call move_alloc(temp_real, thread_data%distances)
   
end subroutine appendToThreadResults


subroutine mergeAndSortThreadResults(thread_data, num_threads, merged_facet_ids, merged_areas, merged_bnd_pts, merged_distances)
   ! Merge thread-local arrays efficiently and sort by facet ID
   implicit none
   type(thread_results_type), dimension(:), intent(in) :: thread_data
   integer, intent(in) :: num_threads
   integer, dimension(:), allocatable, intent(out) :: merged_facet_ids, merged_bnd_pts
   real, dimension(:), allocatable, intent(out) :: merged_areas, merged_distances
   integer :: total_results, current_pos, thread_id, thread_size
   
   ! Count total elements efficiently
   total_results = 0
   do thread_id = 1, num_threads
      total_results = total_results + size(thread_data(thread_id)%areas)
   end do
   
   if (total_results == 0) then
      allocate(merged_facet_ids(0), merged_areas(0), merged_bnd_pts(0), merged_distances(0))
      return
   end if
   
   ! Allocate merged arrays
   allocate(merged_facet_ids(total_results))
   allocate(merged_areas(total_results))
   allocate(merged_bnd_pts(total_results))
   allocate(merged_distances(total_results))
   
   ! Efficiently concatenate thread arrays
   current_pos = 1
   do thread_id = 1, num_threads
      thread_size = size(thread_data(thread_id)%areas)
      if (thread_size > 0) then
         merged_facet_ids(current_pos:current_pos+thread_size-1) = thread_data(thread_id)%facet_ids
         merged_areas(current_pos:current_pos+thread_size-1) = thread_data(thread_id)%areas
         merged_bnd_pts(current_pos:current_pos+thread_size-1) = thread_data(thread_id)%bnd_pts
         merged_distances(current_pos:current_pos+thread_size-1) = thread_data(thread_id)%distances
         current_pos = current_pos + thread_size
      end if
   end do
   
   ! Sort by facet ID for deterministic output
   call sortArraysByFacetId(merged_facet_ids, merged_areas, merged_bnd_pts, merged_distances)
   
end subroutine mergeAndSortThreadResults


subroutine sortArraysByFacetId(facet_ids, areas, bnd_pts, distances)
   ! Sort four arrays based on facet IDs
   implicit none
   integer, dimension(:), allocatable, intent(inout) :: facet_ids, bnd_pts
   real, dimension(:), allocatable, intent(inout) :: areas, distances
   integer :: n, i, j, temp_facet_id, temp_bnd_pt
   real :: temp_area, temp_dist
   
   n = size(facet_ids)
   
   ! Insertion sort based on facet_ids
   do i = 2, n
      ! Store current element values
      temp_facet_id = facet_ids(i)
      temp_area = areas(i)
      temp_bnd_pt = bnd_pts(i)
      temp_dist = distances(i)
      
      j = i - 1
      
      ! Shift elements that are greater than temp_facet_id
      do while (j >= 1 .and. facet_ids(j) > temp_facet_id)
         facet_ids(j + 1) = facet_ids(j)
         areas(j + 1) = areas(j)
         bnd_pts(j + 1) = bnd_pts(j)
         distances(j + 1) = distances(j)
         j = j - 1
      end do
      
      ! Insert the stored elements at correct position
      facet_ids(j + 1) = temp_facet_id
      areas(j + 1) = temp_area
      bnd_pts(j + 1) = temp_bnd_pt
      distances(j + 1) = temp_dist
   end do
   
end subroutine sortArraysByFacetId


subroutine writeFacetSections(secfacids, secareas, secbndptids, bnddst, nfacsecs, fname_facet_sections, fid)
   integer, intent(in) :: nfacsecs
   integer, intent(in), dimension(nfacsecs) :: secfacids, secbndptids
   real   , intent(in), dimension(nfacsecs) :: secareas, bnddst
   character(20), intent(in) :: fname_facet_sections
   integer, intent(in) :: fid
   integer :: n

   open (unit=fid,file=fname_facet_sections,action="write")
   write(fid,*) "# facet      area flux point distance"
   do n=1,nfacsecs
      ! Formatting assumes: #facets < 10 million, #fluid boundary points < 1 billion,
      ! section area < 1000 m^2 (rounded to cm^2), and distance < 1000m
      ! if (bnddst(n) < 0.05*exp(1.)) then
      !    write(*,*) bnddst(n)
      !  end if
      write(unit=fid,fmt='(i8,f10.4,i11,f9.4)') secfacids(n), secareas(n), secbndptids(n), bnddst(n)
   end do
   close (fid)
end subroutine writeFacetSections


subroutine print_facet_sections(nfacsecs_u, nfacsecs_v, nfacsecs_w, nfacsecs_c, &
                                secfacids_u, secareas_u, secbndptids_u, bnddst_u, &
                                secfacids_v, secareas_v, secbndptids_v, bnddst_v, &
                                secfacids_w, secareas_w, secbndptids_w, bnddst_w, &
                                secfacids_c, secareas_c, secbndptids_c, bnddst_c, &
                                filename_u, filename_v, filename_w, filename_c)
   implicit none
   integer,                        intent(in) :: nfacsecs_u, nfacsecs_v, nfacsecs_w, nfacsecs_c
   integer, dimension(nfacsecs_u), intent(in) :: secfacids_u, secbndptids_u
   integer, dimension(nfacsecs_v), intent(in) :: secfacids_v, secbndptids_v
   integer, dimension(nfacsecs_w), intent(in) :: secfacids_w, secbndptids_w
   integer, dimension(nfacsecs_c), intent(in) :: secfacids_c, secbndptids_c
   real   , dimension(nfacsecs_u), intent(in) :: secareas_u, bnddst_u
   real   , dimension(nfacsecs_v), intent(in) :: secareas_v, bnddst_v
   real   , dimension(nfacsecs_w), intent(in) :: secareas_w, bnddst_w
   real   , dimension(nfacsecs_c), intent(in) :: secareas_c, bnddst_c
   character(len=*),               intent(in) :: filename_u, filename_v, filename_w, filename_c

   !$ call OMP_SET_NUM_THREADS(4)
   !$OMP parallel default(shared)
      !$OMP sections

         !$OMP section
         call writeFacetSections(secfacids_u, secareas_u, secbndptids_u, bnddst_u, nfacsecs_u, filename_u, 10)

         !$OMP section
         call writeFacetSections(secfacids_v, secareas_v, secbndptids_v, bnddst_v, nfacsecs_v, filename_v, 20)

         !$OMP section
         call writeFacetSections(secfacids_w, secareas_w, secbndptids_w, bnddst_w, nfacsecs_w, filename_w , 30)

         !$OMP section
         call writeFacetSections(secfacids_c, secareas_c, secbndptids_c, bnddst_c, nfacsecs_c, filename_c, 40)

      !$OMP end sections
   !$OMP end parallel

end subroutine print_facet_sections


subroutine fastPoint2TriMesh(connectivityList, faceNormal, nFaces, vertices, nVertices, &
    queryPoint, minDistance, interceptPoint)
   ! Calculates the minimum distance between a query point and a triangulation (using the incenters)
   ! designed for a single point and a small triangulation (3 <= size <= 6), so probably won't scale well
   ! connectivityList: nFaces x 3(vertices of triangle)
   ! vertices: nVertices x 3(x,y,z)
   ! faceNormal: nFaces x 3(x,y,z)
   ! points: 1 x 3
   use, intrinsic :: ieee_arithmetic
   integer, intent(in) :: nFaces, nVertices, connectivityList(nFaces, 3)
   real, intent(in) :: faceNormal(nFaces, 3), vertices(nVertices, 3), queryPoint(3)!, incenter(nFaces, 3)
   real, intent(out) :: minDistance, interceptPoint(3)
   real :: perpDist, distances(nFaces), projectedPoints(nFaces, 3), &
   sgn, cor1(3), cor2(3), cor3(3), cor4(3), distCornersToEdges(3), edgePoints(3, 3)
   integer :: n, loc
   logical :: check, any_valid

   do n = 1,nFaces
      !perpDist = dot_product(queryPoint - incenter(n, :), faceNormal(n, :))
      perpDist = dot_product(queryPoint - vertices(connectivityList(n,1), :), faceNormal(n, :))
      !sgn = 0
      if (perpDist > 0) then ! the point is 'outside' the triangle
         cor1 = vertices(connectivityList(n, 1), :)
         cor2 = vertices(connectivityList(n, 2), :)
         cor3 = vertices(connectivityList(n, 3), :)
         cor4 = queryPoint - (perpDist * faceNormal(n, :))
         projectedPoints(n, :) = cor4

         ! check if the point is within the triangle using the barycenter
         check = isInsideTriangle(cor1, cor2, cor3, cor4)

         if (check .eqv. .false.) then
            ! The projected point is not within the triangle so the nearest point will lie on the nearest line segment.
            ! Check nearest point to each line segment/edge of face and use minimum distance as the correct projection point
            call projectPointToLineSegment(cor1, cor2, cor4, distCornersToEdges(1), edgePoints(1, :))
            call projectPointToLineSegment(cor2, cor3, cor4, distCornersToEdges(2), edgePoints(2, :))
            call projectPointToLineSegment(cor1, cor3, cor4, distCornersToEdges(3), edgePoints(3, :))
            loc = minloc(distCornersToEdges, 1)
            projectedPoints(n, :) = edgePoints(loc, :)
         end if

         distances(n) = norm2(queryPoint - projectedPoints(n, :))

      else
         distances(n) = ieee_value(distances(n), ieee_quiet_nan)
      end if

   end do

   any_valid = any(.not. ieee_is_nan(distances))
   if (any_valid) then
      loc = minloc(distances, dim=1, mask=.not.ieee_is_nan(distances))
   else
      loc = -1
   end if
   
   if (loc > 0) then
      minDistance = distances(loc)
      interceptPoint = projectedPoints(loc, :)
   else
      minDistance = ieee_value(minDistance, ieee_quiet_nan)
      interceptPoint = 0.0
   end if

end subroutine fastPoint2TriMesh


logical function isInsideTriangle(A, B, C, P)
   real, intent(in), dimension(3) :: A, B, C, P
   real, dimension(3) :: v0, v1, v2
   real :: u, v, dot00, dot01, dot02, dot11, dot12, invDenom

   v0 = C - A
   v1 = B - A
   v2 = P - A

   dot00 = dot_product(v0, v0)
   dot01 = dot_product(v0, v1)
   dot02 = dot_product(v0, v2)
   dot11 = dot_product(v1, v1)
   dot12 = dot_product(v1, v2)

   invDenom = 1. / (dot00 * dot11 - dot01 * dot01)
   u = (dot11 * dot02 - dot01 * dot12) * invDenom
   v = (dot00 * dot12 - dot01 * dot02) * invDenom
   isInsideTriangle = .false.
   if ((u >= 0) .and. (v >= 0) .and. (u + v < 1)) isInsideTriangle = .true.

end function isInsideTriangle


subroutine projectPointToLineSegment(A,B,p,dist,q)
   real, intent(in), dimension(3) :: A, B, p
   real, intent(out) :: dist, q(3)
   real :: AB(3), Ap(3), AB_squared, t

    AB = B - A
    AB_squared = dot_product(AB, AB)

    if (AB_squared == 0) then
        t = 0
    else
        Ap = p - A
        t = dot_product(Ap,AB) / AB_squared
    end if

    if (t < 0.0) then
        q = A
    elseif (t > 1.0) then
        q = B
    else
        q = A + t * AB;
    end if

    dist = norm2(p - q)

end subroutine projectPointToLineSegment



subroutine sutherlandHodgman3D(subjectPolygon, nVertices, clipPlanes, nPlanes, clippedPolygon)
    ! subjectPolygon: nVertices x 3
    ! clipPlanes: nPlanes x 4, where columns are n1, n2, n3, d (where (n1,n2,n3) dot (x,y,z) = d is the plane equation)
   integer, intent(in) :: nVertices, nPlanes
   real, intent(in) :: subjectPolygon(nVertices, 3), clipPlanes(nPlanes, 4)
   real, allocatable, intent(out) :: clippedPolygon(:,:)
   real, allocatable :: inputList(:,:), tempPolygon(:,:), catPolygon(:,:)
   real :: previousVertex(3), intersection(3)
   integer :: i, j, sizeCatPolygon

   ! nVertices = size(subjectPolygon, 1)
   ! nPlanes = size(clipPlanes, 1)

   allocate(tempPolygon(size(subjectPolygon, 1), 3))
   tempPolygon = subjectPolygon

    do i = 1, nPlanes
      !write(*,*) "i", i
        allocate(inputList(size(tempPolygon, 1), 3))
        inputList = tempPolygon

        deallocate(tempPolygon)
        allocate(tempPolygon(0, 3))
        if (size(inputList, 1) > 0) then
            previousVertex = inputList(size(inputList, 1), :)
        end if

        !write(*,*) "inputList", inputList
        !write(*,*) "previousVertex", previousVertex
        do j = 1, size(inputList, 1)
            if (inside(inputList(j, :), clipPlanes(i, :))) then
               !write(*,*) "a"
                if (.not. inside(previousVertex, clipPlanes(i, :))) then
                   !write(*,*) "b"
                    intersection = computeIntersection(clipPlanes(i, :), previousVertex, inputList(j, :))
                    !write(*,*) "intersection", intersection
                    if (.not.(any(ismember_rows(tempPolygon, intersection)))) then
                       call appendToArray2D(tempPolygon, intersection)
                    end if
                    ! sizeCatPolygon = size(tempPolygon, 1)
                    ! allocate(catPolygon(sizeCatPolygon, 3))
                    ! catPolygon = tempPolygon
                    ! deallocate(tempPolygon)
                    ! allocate(tempPolygon(sizeCatPolygon+1, 3))
                    ! tempPolygon(1:sizeCatPolygon, :) = catPolygon
                    ! tempPolygon(sizeCatPolygon+1, :) = intersection
                    ! write(*,*) tempPolygon
                    ! deallocate(catPolygon)
                end if

                if (.not.(any(ismember_rows(tempPolygon, inputList(j, :))))) then
                   call appendToArray2D(tempPolygon, inputList(j, :))
                end if
                ! sizeCatPolygon = size(tempPolygon, 1)
                ! allocate(catPolygon(sizeCatPolygon, 3))
                ! catPolygon = tempPolygon
                ! deallocate(tempPolygon)
                ! allocate(tempPolygon(sizeCatPolygon+1, 3))
                ! tempPolygon(1:sizeCatPolygon, :) = catPolygon
                ! tempPolygon(sizeCatPolygon+1, :) = inputList(j, :)
                ! write(*,*) tempPolygon
                ! deallocate(catPolygon)

            elseif (inside(previousVertex, clipPlanes(i, :))) then
               !write(*,*) "c"
                intersection = computeIntersection(clipPlanes(i, :), previousVertex, inputList(j, :))
                if (.not.(any(ismember_rows(tempPolygon, intersection)))) then
                   !write(*,*) "intersection", intersection
                   call appendToArray2D(tempPolygon, intersection)
                end if
                ! sizeCatPolygon = size(tempPolygon, 1)
                ! allocate(catPolygon(sizeCatPolygon, 3))
                ! catPolygon = tempPolygon
                ! deallocate(tempPolygon)
                ! allocate(tempPolygon(sizeCatPolygon+1, 3))
                ! tempPolygon(1:sizeCatPolygon, :) = catPolygon
                ! tempPolygon(sizeCatPolygon+1, :) = intersection
                ! write(*,*) tempPolygon
                ! deallocate(catPolygon)

            end if
            !rite(*,*) "tempPolygon", tempPolygon
            previousVertex = inputList(j, :)
        end do
        deallocate(inputList)
    end do

    allocate(clippedPolygon(size(tempPolygon, 1), 3))
    clippedPolygon = tempPolygon

end subroutine sutherlandHodgman3D


subroutine appendToArray1D_real(array, item)
   !integer, intent(in) :: nVertices
   real, allocatable, intent(inout) :: array(:)
   real, intent(in) :: item
   real, dimension(:), allocatable :: tempArray
   integer sizeArray

   sizeArray = size(array, 1)
   ! check that size(item,1) = 1
   allocate(tempArray(sizeArray))
   tempArray = array
   deallocate(array)
   allocate(array(sizeArray+1))
   array(1:sizeArray) = tempArray
   array(sizeArray+1) = item
   deallocate(tempArray)

end subroutine appendToArray1D_real


subroutine appendToArray1D_integer(array, item)
   !integer, intent(in) :: nVertices
   integer, allocatable, intent(inout) :: array(:)
   integer, intent(in) :: item
   integer, dimension(:), allocatable :: tempArray
   integer sizeArray

   sizeArray = size(array, 1)
   ! check that size(item,1) = 1
   allocate(tempArray(sizeArray))
   tempArray = array
   deallocate(array)
   allocate(array(sizeArray+1))
   array(1:sizeArray) = tempArray
   array(sizeArray+1) = item
   deallocate(tempArray)

end subroutine appendToArray1D_integer


subroutine appendToArray2D(array, item)
   !integer, intent(in) :: nVertices
   real, allocatable, intent(inout) :: array(:,:)
   real, dimension(:), intent(in) :: item
   real, dimension(:,:), allocatable :: tempArray
   integer sizeArray1, sizeArray2

   sizeArray1 = size(array, 1)
   sizeArray2 = size(array, 2)
   ! check that size(item,1) = sizeArray2
   allocate(tempArray(sizeArray1, sizeArray2))
   tempArray = array
   deallocate(array)
   allocate(array(sizeArray1+1, sizeArray2))
   array(1:sizeArray1, :) = tempArray
   array(sizeArray1+1, :) = item
   deallocate(tempArray)

end subroutine appendToArray2D


function computeIntersection(plane, point1, point2) ! point1 and point2 are points on a line
   real :: computeIntersection(3), point3(3), dist
   real, intent(in) :: point1(3), point2(3), plane(4)

   point3 = plane(1:3) * plane(4)
   dist = dot_product(plane(1:3), point3 - point1) / dot_product(plane(1:3), point2 - point1)
   computeIntersection = point1 + dist * (point2 - point1)

end function computeIntersection


logical function inside(point, plane)
   real :: point(3), plane(4)

   if ((dot_product(point, plane(1:3)) - plane(4)) <= 0.) then
      inside = .true.
   else
      inside = .false.
   end if

end function

function cross_product(a,b)
  ! Calculate the cross product (a x b)
  implicit none
  real, dimension(3) :: cross_product
  real, dimension(3), intent(in) :: a, b

  cross_product(1) = a(2)*b(3) - a(3)*b(2)
  cross_product(2) = a(3)*b(1) - a(1)*b(3)
  cross_product(3) = a(1)*b(2) - a(2)*b(1)

end function cross_product


FUNCTION polyarea(x, y, nb) RESULT(fn_val)

! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-04  Time: 12:24:06

IMPLICIT NONE

REAL, INTENT(IN)     :: x(:)
REAL, INTENT(IN)     :: y(:)
INTEGER, INTENT(IN)  :: nb
REAL                 :: fn_val

!*****************************************************************

!   GIVEN A SEQUENCE OF NB POINTS (X(I),Y(I)),  polyarea COMPUTES THE AREA
! BOUNDED BY THE CLOSED POLYGONAL CURVE WHICH PASSES THROUGH THE POINTS IN
! THE ORDER THAT THEY ARE INDEXED.  THE FINAL POINT OF THE CURVE IS ASSUMED
! TO BE THE FIRST POINT GIVEN.  THEREFORE, IT NEED NOT BE LISTED AT THE END
! OF X AND Y.  THE CURVE IS NOT REQUIRED TO BE SIMPLE.  e.g. It may cross over
! itself.

!*****************************************************************

INTEGER  :: i, n, nm1
REAL     :: a

n = nb
IF (x(1) == x(n) .AND. y(1) == y(n)) n = n - 1

SELECT CASE (n)
  CASE (:2)
     fn_val = 0.0

  CASE (3)
     fn_val= 0.5*((x(2) - x(1))*(y(3) - y(1)) - (x(3) - x(1))*(y(2) - y(1)))

  CASE DEFAULT
    nm1 = n - 1
    a = x(1)*(y(2) - y(n)) + x(n)*(y(1) - y(nm1))
    DO  i = 2, nm1
      a = a + x(i)*(y(i+1) - y(i-1))
    END DO
    fn_val = 0.5*a
END SELECT

RETURN
END FUNCTION polyarea


function ismember_rows(a,b) result(c)
implicit none
real, intent(in) :: a(:,:), b(:)
logical, allocatable :: c(:)

integer :: n, m, i

n = size(a,1)
m = size(a,2)
if (size(b) /= m) stop "b size is wrong"

allocate( c(n) )
do i = 1, n
    c(i) = all( b(:) == a(i,:) )
end do

end function ismember_rows


end module matchFacets2Cells
