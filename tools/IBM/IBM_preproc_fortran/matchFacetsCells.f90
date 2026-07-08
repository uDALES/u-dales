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
   use iso_fortran_env, only: int64
   use omp_lib
   implicit none

   ! Derived type to hold thread-local results
   type :: thread_results_type
      real, dimension(:), allocatable :: areas
      integer, dimension(:), allocatable :: facet_ids
      integer, dimension(:), allocatable :: bnd_pts
      real, dimension(:), allocatable :: distances
      integer :: nused = 0
      integer :: capacity = 0
   end type thread_results_type

   ! Sparse boundary-cell index. Records are ordered as the original
   ! matchFacetsToCells cell loop: i outer, j middle, k inner.
   type :: boundary_index_type
      integer, dimension(:,:), allocatable :: ij_start
      integer, dimension(:), allocatable :: k
      logical, dimension(:), allocatable :: is_fluid
      logical, dimension(:), allocatable :: is_solid
      integer, dimension(:), allocatable :: fluid_id
      integer, dimension(:), allocatable :: cand_start
      integer, dimension(:), allocatable :: cand_code
      integer, dimension(:), allocatable :: cand_i
      integer, dimension(:), allocatable :: cand_j
      integer, dimension(:), allocatable :: cand_k
      integer, dimension(:), allocatable :: cand_fluid_id
      integer(int64), dimension(:), allocatable :: fluid_keys
   end type boundary_index_type

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
    type(boundary_index_type) :: bidx
    !$ type(thread_results_type), dimension(:), allocatable :: thread_data
    !$ integer :: actual_threads
    real, dimension(:,:), allocatable :: clipVertices
    real, dimension(:), allocatable :: x_edges, y_edges, z_edges
    logical :: search_adj
    integer :: il, iu, jl, ju, kl, ku, nClipVertices, id, loc
    integer :: bnd_first, bnd_last, bnd_idx, cand_pos, cand_code, candidate_ids(27)
    real :: dx, dy, dz, xmin, xmax, ymin, ymax, zmin, zmax, xl, xu, yl, yu, zl, zu, tol, &
            planes(6,4), area, area_miss, xproj, yproj, zproj, proj, projVec(3), projArea, &
            dist, dists(27), angle, angles(27), BI(3), BIs(27,3)
    real, dimension(3) :: xyz, xyz1
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
    !$     thread_data(i)%nused = 0
    !$     thread_data(i)%capacity = 0
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

    allocate(x_edges(itot+1), y_edges(jtot+1), z_edges(ktot+1))
    x_edges(1:itot) = xgrid - dx/2.
    x_edges(itot+1) = xgrid(itot) + dx/2.
    y_edges(1:jtot) = ygrid - dy/2.
    y_edges(jtot+1) = ygrid(jtot) + dy/2.
    z_edges(1:ktot) = zgrid - dz/2.
    z_edges(ktot+1) = zgrid(ktot) + dz/2.

    call buildBoundaryIndex(fluid_IB, solid_IB, itot, jtot, ktot, diag_neighbs, bidx)

    !$OMP parallel do default(none) &
    !$OMP shared(nFaces, connectivityList, faceNormal, vertices, xgrid, ygrid, zgrid, tol, &
    !$OMP        fluid_IB, solid_IB, itot, jtot, ktot, &
    !$OMP        diag_neighbs, periodic_x, periodic_y, dx, dy, dz, &
    !$OMP        x_edges, y_edges, z_edges, bidx, thread_data, actual_threads) &
    !$OMP private(xmin, xmax, ymin, ymax, zmin, zmax, &
    !$OMP         il, iu, jl, ju, kl, ku, xl, xu, yl, yu, zl, zu, planes, &
    !$OMP         clipVertices, nClipVertices, &
    !$OMP         xproj, yproj, zproj, proj, projVec, projArea, &
    !$OMP         area, dist, dists, angle, angles, BI, BIs, &
    !$OMP         xyz, xyz1, search_adj, id, loc, candidate_ids, &
    !$OMP         bnd_first, bnd_last, bnd_idx, cand_pos, cand_code, &
    !$OMP         i, j, k, m, p, dir, ids, thread_id) &
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

      il = lastEdgeLeq(x_edges, xmin + tol)
      iu = firstEdgeGeq(x_edges, xmax - tol)
      jl = lastEdgeLeq(y_edges, ymin + tol)
      ju = firstEdgeGeq(y_edges, ymax - tol)
      kl = lastEdgeLeq(z_edges, zmin + tol)
      ku = firstEdgeGeq(z_edges, zmax - tol)

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
            bnd_first = lowerBoundInt(bidx%k, bidx%ij_start(i,j), &
                                      bidx%ij_start(i,j+1)-1, kl)
            bnd_last = upperBoundInt(bidx%k, bidx%ij_start(i,j), &
                                     bidx%ij_start(i,j+1)-1, ku) - 1

            do bnd_idx=bnd_first,bnd_last
               k = bidx%k(bnd_idx)
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
                  deallocate(clipVertices)
                  cycle
               elseif (nClipVertices == 3) then
                  area = 0.5*norm2(cross_product(clipVertices(2,:) - clipVertices(1,:), &
                     clipVertices(3,:) - clipVertices(1,:)))
                  ! remove anything below 1 square centimetre, as we only write to this precision.
                  if (area < 1e-5) then
                     deallocate(clipVertices)
                     cycle
                  end if

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
                     ids = (/2, 3/)
                  elseif (dir==2) then
                     ids = (/1, 3/)
                  elseif (dir==3) then
                     ids = (/1, 2/)
                  end if

                  projArea = abs(polyarea(clipVertices(:,ids(1)), clipVertices(:,ids(2)), nClipVertices))

                  area = projArea / proj

                  if (area < 1e-5) then
                     deallocate(clipVertices)
                     cycle
                  end if

               else
                 write(*,*) "something wrong with clipped polygon"
               end if

               if (((xgrid(i) == 0.) .and. periodic_x) .or. ((ygrid(j) == 0.) .and. periodic_y)) then ! Account for periodicity - flux at point N+1
                  area = area * 2.
               end if

               dists = ieee_value(dists, ieee_quiet_nan)
               angles = ieee_value(angles, ieee_quiet_nan)
               BIs = ieee_value(BIs, ieee_quiet_nan)
               candidate_ids = 0
               search_adj = .false.

               if (bidx%is_fluid(bnd_idx)) then
                  xyz1 = (/xgrid(i), ygrid(j), zgrid(k)/)
                  loc = bidx%fluid_id(bnd_idx)
                  call fastPoint2ClippedFacet(faceNormal(n,:), clipVertices, nClipVertices, &
                     xyz1, dist, BI)
                  angle = dot_product(faceNormal(n,:), (xyz1 - BI)/norm2(xyz1 - BI))

                  if (abs(angle - 1.) < epsilon(angle)) then ! Wall-normal defined, use this cell
                  !if (abs(angle - 1.) < epsilon(angle) .and. dist > 0.05*exp(1.)) then ! Wall-normal defined, use this cell
                     id = 1 ! not necessary?
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
                        candidate_ids(1) = loc
                     end if
                  end if
               end if

               if (bidx%is_solid(bnd_idx) .or. search_adj) then

                  do cand_pos=bidx%cand_start(bnd_idx), bidx%cand_start(bnd_idx+1)-1
                     cand_code = bidx%cand_code(cand_pos)
                     xyz = (/xgrid(bidx%cand_i(cand_pos)), &
                             ygrid(bidx%cand_j(cand_pos)), &
                             zgrid(bidx%cand_k(cand_pos))/)
                     call fastPoint2ClippedFacet(faceNormal(n,:), clipVertices, nClipVertices, &
                        xyz, dist, BI)
                     dists(cand_code) = dist
                     angles(cand_code) = dot_product(faceNormal(n,:), (xyz - BI)/norm2(xyz - BI))
                     BIs(cand_code,:) = BI
                     candidate_ids(cand_code) = bidx%cand_fluid_id(cand_pos)
                  end do
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
                     deallocate(clipVertices)
                     cycle
                  end if

                  loc = candidate_ids(id)

               end if !(solid_IB(i,j,k) .or. search_adj)

               if (isnan(dist) .or. abs(dist)<1e-4) dist = 0.1 ! 0.1 m minimum distance to avoid singularities
               if (isnan(area)) area = 0.0

               !! For serial run version: Append to global arrays
               ! call appendToArray1D_real(secareas, area)
               ! call appendToArray1D_integer(secfacids, n)
               ! call appendToArray1D_integer(secbndptids, loc)
               ! call appendToArray1D_real(bnddst, abs(dist))
               
               ! Alternatively for OpenMP run: Append to thread-specific arrays
               !$ call appendToThreadResults(thread_data(thread_id), area, n, loc, abs(dist))

               deallocate(clipVertices)
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


subroutine buildBoundaryIndex(fluid_IB, solid_IB, itot, jtot, ktot, diag_neighbs, bidx)
   implicit none
   integer, intent(in) :: itot, jtot, ktot
   logical, intent(in), dimension(itot,jtot,ktot) :: fluid_IB, solid_IB
   logical, intent(in) :: diag_neighbs
   type(boundary_index_type), intent(out) :: bidx

   integer, allocatable :: ij_count(:,:), ij_next(:,:)
   integer :: i, j, k, idx, pos, nBoundary, nFluid, fid
   integer :: totalCand, cand_count

   allocate(ij_count(itot,jtot))
   ij_count = 0
   nBoundary = 0
   nFluid = 0

   do j=1,jtot
      do k=1,ktot
         do i=1,itot
            if (fluid_IB(i,j,k)) nFluid = nFluid + 1
            if (fluid_IB(i,j,k) .or. solid_IB(i,j,k)) then
               ij_count(i,j) = ij_count(i,j) + 1
               nBoundary = nBoundary + 1
            end if
         end do
      end do
   end do

   allocate(bidx%fluid_keys(nFluid))
   fid = 0
   do j=1,jtot
      do k=1,ktot
         do i=1,itot
            if (fluid_IB(i,j,k)) then
               fid = fid + 1
               bidx%fluid_keys(fid) = cellKey(i, j, k, itot, ktot)
            end if
         end do
      end do
   end do

   allocate(bidx%ij_start(itot,jtot+1))
   pos = 1
   do i=1,itot
      do j=1,jtot
         bidx%ij_start(i,j) = pos
         pos = pos + ij_count(i,j)
      end do
      bidx%ij_start(i,jtot+1) = pos
   end do

   allocate(bidx%k(nBoundary), bidx%is_fluid(nBoundary), bidx%is_solid(nBoundary))
   allocate(bidx%fluid_id(nBoundary))
   allocate(ij_next(itot,jtot))
   ij_next = bidx%ij_start(:,1:jtot)

   do i=1,itot
      do j=1,jtot
         do k=1,ktot
            if (fluid_IB(i,j,k) .or. solid_IB(i,j,k)) then
               idx = ij_next(i,j)
               bidx%k(idx) = k
               bidx%is_fluid(idx) = fluid_IB(i,j,k)
               bidx%is_solid(idx) = solid_IB(i,j,k)
               if (fluid_IB(i,j,k)) then
                  bidx%fluid_id(idx) = fluidIdFromKey(bidx%fluid_keys, &
                     cellKey(i, j, k, itot, ktot))
               else
                  bidx%fluid_id(idx) = 0
               end if
               ij_next(i,j) = idx + 1
            end if
         end do
      end do
   end do

   allocate(bidx%cand_start(nBoundary+1))
   totalCand = 0
   do i=1,itot
      do j=1,jtot
         do idx=bidx%ij_start(i,j), bidx%ij_start(i,j+1)-1
            bidx%cand_start(idx) = totalCand + 1
            cand_count = 0
            call addReceiverCandidatesForCell(fluid_IB, itot, jtot, ktot, i, j, bidx%k(idx), &
               diag_neighbs, bidx%fluid_keys, cand_count)
            totalCand = totalCand + cand_count
         end do
      end do
   end do
   bidx%cand_start(nBoundary+1) = totalCand + 1

   allocate(bidx%cand_code(totalCand), bidx%cand_i(totalCand), bidx%cand_j(totalCand))
   allocate(bidx%cand_k(totalCand), bidx%cand_fluid_id(totalCand))
   totalCand = 0
   do i=1,itot
      do j=1,jtot
         do idx=bidx%ij_start(i,j), bidx%ij_start(i,j+1)-1
            call addReceiverCandidatesForCell(fluid_IB, itot, jtot, ktot, i, j, bidx%k(idx), &
               diag_neighbs, bidx%fluid_keys, totalCand, bidx%cand_code, bidx%cand_i, &
               bidx%cand_j, bidx%cand_k, bidx%cand_fluid_id)
         end do
      end do
   end do

   deallocate(ij_count, ij_next)
end subroutine buildBoundaryIndex


subroutine addReceiverCandidatesForCell(fluid_IB, itot, jtot, ktot, i, j, k, diag_neighbs, &
   fluid_keys, cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
   implicit none
   integer, intent(in) :: itot, jtot, ktot, i, j, k
   logical, intent(in), dimension(itot,jtot,ktot) :: fluid_IB
   logical, intent(in) :: diag_neighbs
   integer(int64), intent(in), dimension(:) :: fluid_keys
   integer, intent(inout) :: cand_count
   integer, intent(inout), optional, dimension(:) :: cand_code, cand_i, cand_j, cand_k, cand_fluid_id

   call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i-1, j, k, 2, fluid_keys, &
      cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
   call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i+1, j, k, 3, fluid_keys, &
      cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
   call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i, j-1, k, 4, fluid_keys, &
      cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
   call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i, j+1, k, 5, fluid_keys, &
      cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
   call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i, j, k-1, 6, fluid_keys, &
      cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
   call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i, j, k+1, 7, fluid_keys, &
      cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)

   if (diag_neighbs) then
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i-1, j-1, k, 8, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i-1, j+1, k, 9, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i+1, j-1, k, 10, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i+1, j+1, k, 11, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i-1, j, k-1, 12, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i-1, j, k+1, 13, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i+1, j, k-1, 14, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i+1, j, k+1, 15, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i, j-1, k-1, 16, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i, j-1, k+1, 17, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i, j+1, k-1, 18, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i, j+1, k+1, 19, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i-1, j-1, k-1, 20, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i+1, j-1, k-1, 21, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i-1, j+1, k-1, 22, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i+1, j+1, k-1, 23, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i-1, j-1, k+1, 24, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i+1, j-1, k+1, 25, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i-1, j+1, k+1, 26, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
      call addCandidateIfFluid(fluid_IB, itot, jtot, ktot, i+1, j+1, k+1, 27, fluid_keys, &
         cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
   end if
end subroutine addReceiverCandidatesForCell


subroutine addCandidateIfFluid(fluid_IB, itot, jtot, ktot, ci, cj, ck, code, fluid_keys, &
   cand_count, cand_code, cand_i, cand_j, cand_k, cand_fluid_id)
   implicit none
   integer, intent(in) :: itot, jtot, ktot, ci, cj, ck, code
   logical, intent(in), dimension(itot,jtot,ktot) :: fluid_IB
   integer(int64), intent(in), dimension(:) :: fluid_keys
   integer, intent(inout) :: cand_count
   integer, intent(inout), optional, dimension(:) :: cand_code, cand_i, cand_j, cand_k, cand_fluid_id

   if (ci < 1 .or. ci > itot) return
   if (cj < 1 .or. cj > jtot) return
   if (ck < 1 .or. ck > ktot) return
   if (.not. fluid_IB(ci,cj,ck)) return

   cand_count = cand_count + 1
   if (present(cand_code)) then
      cand_code(cand_count) = code
      cand_i(cand_count) = ci
      cand_j(cand_count) = cj
      cand_k(cand_count) = ck
      cand_fluid_id(cand_count) = fluidIdFromKey(fluid_keys, cellKey(ci, cj, ck, itot, ktot))
   end if
end subroutine addCandidateIfFluid


integer(int64) function cellKey(i, j, k, itot, ktot)
   implicit none
   integer, intent(in) :: i, j, k, itot, ktot

   cellKey = int(i, int64) + int(k-1, int64) * int(itot, int64) + &
             int(j-1, int64) * int(itot, int64) * int(ktot, int64)
end function cellKey


integer function fluidIdFromKey(fluid_keys, key)
   implicit none
   integer(int64), intent(in), dimension(:) :: fluid_keys
   integer(int64), intent(in) :: key
   integer :: lo, hi, mid

   lo = 1
   hi = size(fluid_keys)
   fluidIdFromKey = 0
   do while (lo <= hi)
      mid = (lo + hi) / 2
      if (fluid_keys(mid) == key) then
         fluidIdFromKey = mid
         return
      elseif (fluid_keys(mid) < key) then
         lo = mid + 1
      else
         hi = mid - 1
      end if
   end do
end function fluidIdFromKey


integer function lowerBoundInt(values, first, last, target)
   implicit none
   integer, intent(in), dimension(:) :: values
   integer, intent(in) :: first, last, target
   integer :: lo, hi, mid

   lo = first
   hi = last + 1
   do while (lo < hi)
      mid = (lo + hi) / 2
      if (values(mid) < target) then
         lo = mid + 1
      else
         hi = mid
      end if
   end do
   lowerBoundInt = lo
end function lowerBoundInt


integer function upperBoundInt(values, first, last, target)
   implicit none
   integer, intent(in), dimension(:) :: values
   integer, intent(in) :: first, last, target
   integer :: lo, hi, mid

   lo = first
   hi = last + 1
   do while (lo < hi)
      mid = (lo + hi) / 2
      if (values(mid) <= target) then
         lo = mid + 1
      else
         hi = mid
      end if
   end do
   upperBoundInt = lo
end function upperBoundInt


integer function lastEdgeLeq(edges, value)
   implicit none
   real, intent(in), dimension(:) :: edges
   real, intent(in) :: value
   integer :: lo, hi, mid

   lo = 1
   hi = size(edges) + 1
   do while (lo < hi)
      mid = (lo + hi) / 2
      if (edges(mid) <= value) then
         lo = mid + 1
      else
         hi = mid
      end if
   end do
   lastEdgeLeq = lo - 1
end function lastEdgeLeq


integer function firstEdgeGeq(edges, value)
   implicit none
   real, intent(in), dimension(:) :: edges
   real, intent(in) :: value
   integer :: lo, hi, mid

   lo = 1
   hi = size(edges) + 1
   do while (lo < hi)
      mid = (lo + hi) / 2
      if (edges(mid) < value) then
         lo = mid + 1
      else
         hi = mid
      end if
   end do
   if (lo > size(edges)) then
      firstEdgeGeq = 0
   else
      firstEdgeGeq = lo
   end if
end function firstEdgeGeq


subroutine appendToThreadResults(thread_data, area, facet_id, bnd_pt, distance)
   implicit none
   type(thread_results_type), intent(inout) :: thread_data
   real, intent(in) :: area, distance
   integer, intent(in) :: facet_id, bnd_pt

   integer :: n

   if (thread_data%nused >= thread_data%capacity) then
      call reserveThreadResults(thread_data, max(1024, max(1, 2*thread_data%capacity)))
   end if

   n = thread_data%nused + 1
   thread_data%areas(n) = area
   thread_data%facet_ids(n) = facet_id
   thread_data%bnd_pts(n) = bnd_pt
   thread_data%distances(n) = distance
   thread_data%nused = n
end subroutine appendToThreadResults


subroutine reserveThreadResults(thread_data, new_capacity)
   implicit none
   type(thread_results_type), intent(inout) :: thread_data
   integer, intent(in) :: new_capacity

   real, dimension(:), allocatable :: temp_real
   integer, dimension(:), allocatable :: temp_int
   integer :: n

   if (new_capacity <= thread_data%capacity) return

   n = thread_data%nused

   allocate(temp_real(new_capacity))
   if (n > 0) temp_real(1:n) = thread_data%areas(1:n)
   call move_alloc(temp_real, thread_data%areas)

   allocate(temp_int(new_capacity))
   if (n > 0) temp_int(1:n) = thread_data%facet_ids(1:n)
   call move_alloc(temp_int, thread_data%facet_ids)

   allocate(temp_int(new_capacity))
   if (n > 0) temp_int(1:n) = thread_data%bnd_pts(1:n)
   call move_alloc(temp_int, thread_data%bnd_pts)

   allocate(temp_real(new_capacity))
   if (n > 0) temp_real(1:n) = thread_data%distances(1:n)
   call move_alloc(temp_real, thread_data%distances)

   thread_data%capacity = new_capacity
end subroutine reserveThreadResults


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
      total_results = total_results + thread_data(thread_id)%nused
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
      thread_size = thread_data(thread_id)%nused
      if (thread_size > 0) then
         merged_facet_ids(current_pos:current_pos+thread_size-1) = thread_data(thread_id)%facet_ids(1:thread_size)
         merged_areas(current_pos:current_pos+thread_size-1) = thread_data(thread_id)%areas(1:thread_size)
         merged_bnd_pts(current_pos:current_pos+thread_size-1) = thread_data(thread_id)%bnd_pts(1:thread_size)
         merged_distances(current_pos:current_pos+thread_size-1) = thread_data(thread_id)%distances(1:thread_size)
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
   integer, dimension(:), allocatable :: counts, next_pos, sorted_facet_ids, sorted_bnd_pts
   real, dimension(:), allocatable :: sorted_areas, sorted_distances
   integer :: n, i, facet_id, pos, max_facet_id
   
   n = size(facet_ids)
   if (n <= 1) return

   max_facet_id = maxval(facet_ids)
   allocate(counts(max_facet_id), next_pos(max_facet_id))
   counts = 0
   do i = 1, n
      counts(facet_ids(i)) = counts(facet_ids(i)) + 1
   end do

   pos = 1
   do i = 1, max_facet_id
      next_pos(i) = pos
      pos = pos + counts(i)
   end do

   allocate(sorted_facet_ids(n), sorted_bnd_pts(n), sorted_areas(n), sorted_distances(n))
   do i = 1, n
      facet_id = facet_ids(i)
      pos = next_pos(facet_id)
      sorted_facet_ids(pos) = facet_ids(i)
      sorted_areas(pos) = areas(i)
      sorted_bnd_pts(pos) = bnd_pts(i)
      sorted_distances(pos) = distances(i)
      next_pos(facet_id) = pos + 1
   end do

   call move_alloc(sorted_facet_ids, facet_ids)
   call move_alloc(sorted_areas, areas)
   call move_alloc(sorted_bnd_pts, bnd_pts)
   call move_alloc(sorted_distances, distances)
   
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


subroutine fastPoint2ClippedFacet(facetNormal, vertices, nVertices, queryPoint, minDistance, interceptPoint)
   ! Same triangle-fan search as fastPoint2TriMesh for a clipped facet whose
   ! fan connectivity is (1,2,3), (1,3,4), ...
   use, intrinsic :: ieee_arithmetic
   integer, intent(in) :: nVertices
   real, intent(in) :: facetNormal(3), vertices(nVertices, 3), queryPoint(3)
   real, intent(out) :: minDistance, interceptPoint(3)
   real :: perpDist, distances(max(1, nVertices-2)), projectedPoints(max(1, nVertices-2), 3), &
      cor1(3), cor2(3), cor3(3), cor4(3), distCornersToEdges(3), edgePoints(3, 3)
   integer :: n, nFaces, loc
   logical :: check, any_valid

   nFaces = nVertices - 2
   if (nFaces < 1) then
      minDistance = ieee_value(minDistance, ieee_quiet_nan)
      interceptPoint = 0.0
      return
   end if

   do n = 1,nFaces
      perpDist = dot_product(queryPoint - vertices(1, :), facetNormal)
      if (perpDist > 0) then ! the point is 'outside' the triangle
         cor1 = vertices(1, :)
         cor2 = vertices(n+1, :)
         cor3 = vertices(n+2, :)
         cor4 = queryPoint - (perpDist * facetNormal)
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

   any_valid = any(.not. ieee_is_nan(distances(1:nFaces)))
   if (any_valid) then
      loc = minloc(distances(1:nFaces), dim=1, mask=.not.ieee_is_nan(distances(1:nFaces)))
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

end subroutine fastPoint2ClippedFacet


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
   implicit none
   integer, intent(in) :: nVertices, nPlanes
   real, intent(in) :: subjectPolygon(nVertices, 3), clipPlanes(nPlanes, 4)
   real, allocatable, intent(out) :: clippedPolygon(:,:)
   real :: inputList(max(1, nVertices + 2*nPlanes + 2), 3)
   real :: tempPolygon(max(1, nVertices + 2*nPlanes + 2), 3)
   real :: previousVertex(3), intersection(3)
   integer :: i, j, inputCount, tempCount

   tempPolygon(1:nVertices, :) = subjectPolygon
   tempCount = nVertices

    do i = 1, nPlanes
      !write(*,*) "i", i
        inputCount = tempCount
        if (inputCount > 0) then
            inputList(1:inputCount, :) = tempPolygon(1:inputCount, :)
            previousVertex = inputList(inputCount, :)
        end if
        tempCount = 0

        !write(*,*) "inputList", inputList
        !write(*,*) "previousVertex", previousVertex
        do j = 1, inputCount
            if (inside(inputList(j, :), clipPlanes(i, :))) then
               !write(*,*) "a"
                if (.not. inside(previousVertex, clipPlanes(i, :))) then
                   !write(*,*) "b"
                    intersection = computeIntersection(clipPlanes(i, :), previousVertex, inputList(j, :))
                    !write(*,*) "intersection", intersection
                    call appendUniqueRow3(tempPolygon, tempCount, intersection)
                end if

                call appendUniqueRow3(tempPolygon, tempCount, inputList(j, :))

            elseif (inside(previousVertex, clipPlanes(i, :))) then
               !write(*,*) "c"
                intersection = computeIntersection(clipPlanes(i, :), previousVertex, inputList(j, :))
                !write(*,*) "intersection", intersection
                call appendUniqueRow3(tempPolygon, tempCount, intersection)

            end if
            !rite(*,*) "tempPolygon", tempPolygon
            previousVertex = inputList(j, :)
        end do
    end do

    allocate(clippedPolygon(tempCount, 3))
    if (tempCount > 0) clippedPolygon = tempPolygon(1:tempCount, :)

end subroutine sutherlandHodgman3D


subroutine appendUniqueRow3(array, nRows, item)
   implicit none
   real, intent(inout) :: array(:,:)
   integer, intent(inout) :: nRows
   real, intent(in) :: item(3)
   integer :: row

   do row = 1,nRows
      if (array(row,1) == item(1) .and. array(row,2) == item(2) .and. array(row,3) == item(3)) return
   end do

   if (nRows >= size(array, 1)) then
      write(*,*) "sutherlandHodgman3D output exceeds local buffer"
      stop
   end if

   nRows = nRows + 1
   array(nRows, :) = item

end subroutine appendUniqueRow3


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
