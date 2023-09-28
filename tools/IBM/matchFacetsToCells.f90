program run
   implicit none

   real ::tol, Ray_dir_u(3), Ray_dir_v(3), Ray_dir_w(3), Ray_dir_c(3)
   integer :: nVertices, nFaces, stl_ground_i, diag_neighbs_i, n_threads
   logical :: stl_ground, diag_neighbs, periodic_x, periodic_y
   real :: dx, dy  !, dz
   integer :: itot, jtot, ktot, i, j, k, nfluid_IB_c
   real, allocatable, dimension(:) :: xf, xh, yf, yh, zf, zh
   real, allocatable, dimension(:,:) :: vertices, faceNormals
   integer, allocatable, dimension(:,:) :: faces
   logical, dimension(:,:,:), allocatable :: fluid_IB_c, solid_IB_c
   real, dimension(:,:), allocatable :: fluid_IB_xyz_c
   integer, dimension(:), allocatable :: secfacids_c, secbndptids_c
   real, dimension(:), allocatable :: secareas_c, bnddst_c

   open(unit=50,file='inmypoly_inp_info.txt')
   read(unit=50,fmt='(f15.10,x,f15.10)') dx, dy  !, dz
   read(unit=50,fmt='(i5,x,i5,x,i5)') itot, jtot, ktot
   read(unit=50,fmt='(f15.10)') tol
   read(unit=50,fmt='(f15.10,x,f15.10,x,f15.10)') Ray_dir_u(1), Ray_dir_u(2), Ray_dir_u(3)
   read(unit=50,fmt='(f15.10,x,f15.10,x,f15.10)') Ray_dir_v(1), Ray_dir_v(2), Ray_dir_v(3)
   read(unit=50,fmt='(f15.10,x,f15.10,x,f15.10)') Ray_dir_w(1), Ray_dir_w(2), Ray_dir_w(3)
   read(unit=50,fmt='(f15.10,x,f15.10,x,f15.10)') Ray_dir_c(1), Ray_dir_c(2), Ray_dir_c(3)
   read(unit=50,fmt='(i8,x,i8)') nVertices, nFaces
   read(unit=50,fmt='(i4)') n_threads
   read(unit=50,fmt='(i1,x,i1)') stl_ground_i, diag_neighbs_i
   close(unit=50)

   allocate(xf(itot))
   allocate(xh(itot))
   allocate(yf(jtot))
   allocate(yh(jtot))
   allocate(zf(ktot))
   allocate(zh(ktot))
   allocate(faces(nFaces, 3))
   allocate(faceNormals(nFaces, 3))
   allocate(vertices(nVertices, 3))
   allocate(fluid_IB_c(itot,jtot,ktot), solid_IB_c(itot,jtot,ktot))
   allocate(fluid_IB_xyz_c(nfluid_IB_c,3))

   xf(1:itot) = [((i-1)*dx+(dx/2.0), i=1,itot)]
   xh(1:itot) = [((i-1)*dx, i=1,itot)]

   yf(1:jtot) = [((j-1)*dy+(dy/2.0), j=1,jtot)]
   yh(1:jtot) = [((j-1)*dy, j=1,jtot)]

   open(unit=3,file='zfgrid.txt')
   do i = 1,ktot
      read(unit=3,fmt='(f15.10)') zf(i)
   end do
   close(unit=3)

   open(unit=4,file='zhgrid.txt')
   do k = 1,ktot
      read(unit=4,fmt='(f15.10)') zh(i)
   end do
   close(unit=4)

   ! dx = 1
   ! dy = 1
   ! dz = 1
   ! do i=1,itot
   !   xgrid(i) = (i-1)*dx+dx/2
   ! end do
   ! do j=1,jtot
   !   ygrid(j) = (j-1)*dy+dy/2
   ! end do
   ! do k=1,ktot
   !   zgrid(k) = (k-1)*dz+dz/2
   ! end do

   ! write(*,*) xgrid
   ! write(*,*) ygrid
   ! write(*,*) zgrid

   ! periodic_x = .false.
   ! periodic_y = .false.
   ! diag_neighbs = .false.

   if (stl_ground_i==1) then
        stl_ground = .true.
    elseif (stl_ground_i==0) then
        stl_ground = .false.
    else
        print *, "Incorrect input for 'stl_ground'."
    end if

    if (diag_neighbs_i==1) then
        diag_neighbs = .true.
    elseif (diag_neighbs_i==0) then
        diag_neighbs = .false.
    else
        print *, "Incorrect input for 'diag_neighbs'."
    end if

   ! ! grid
   ! itot = 4
   ! jtot = 4
   ! ktot = 4

   nfluid_IB_c = 24
   periodic_x = .false.
   periodic_y = .false.

   call readGeometry('faces.txt', 'faceNormals.txt', nFaces, 'vertices.txt', nVertices, &
     faces, faceNormals, vertices)

   call readBoundaryPoints(xf, yf, zf, itot, jtot, ktot, &
     'fluid_IB_c.txt', 'solid_IB_c.txt', 'fluid_boundary_c.txt' , nfluid_IB_c, &
      fluid_IB_c,       solid_IB_c,       fluid_IB_xyz_c)

   call matchFacetsToCells(faces, faceNormals, nFaces, vertices, nVertices, &
      fluid_IB_c, solid_IB_c, fluid_IB_xyz_c, nfluid_IB_c, xf, yf, zf, itot, jtot, ktot, &
      diag_neighbs, periodic_x, periodic_y, secfacids_c, secbndptids_c, secareas_c, bnddst_c)

   contains


   subroutine readGeometry(fname_faces, fname_faceNormals, nFaces, fname_vertices, nVertices, &
                                 faces,       faceNormals,               vertices)
     character(9), intent(in) :: fname_faces
     character(15), intent(in) :: fname_faceNormals
     character(12), intent(in) :: fname_vertices
     integer, intent(in) :: nFaces, nVertices
     integer, intent(out) :: faces(nFaces,3)
     real   , intent(out) :: faceNormals(nFaces,3), vertices(nVertices,3)
     integer, parameter :: ifinput = 1
     integer :: n

     open (ifinput, file=fname_faces)
     do n = 1,nFaces
       read (ifinput, *) faces(n,1), faces(n,2), faces(n,3)
     end do
     close (ifinput)

     write(*,*) "faces", faces


     open (ifinput, file=fname_faceNormals)
     do n = 1,nFaces
       read (ifinput, *) faceNormals(n,1), faceNormals(n,2), faceNormals(n,3)
     end do
     close (ifinput)

     write(*,*) "faceNormals", faceNormals

     open (ifinput, file=fname_vertices)
     do n = 1,nVertices
       read (ifinput, *) vertices(n,1), vertices(n,2), vertices(n,3)
     end do
     close (ifinput)

     write(*,*) "vertices", vertices

   end subroutine readGeometry


   subroutine readBoundaryPoints(xgrid, ygrid, zgrid, itot, jtot, ktot, &
     fname_fluid_IB, fname_solid_IB,  fname_fluid_boundary, nfluid_IB, &
         fluid_IB,         solid_IB,        fluid_IB_xyz)

   character(14), intent(in) :: fname_fluid_IB, fname_solid_IB
   character(20), intent(in) :: fname_fluid_boundary
   integer, intent(in) :: itot, jtot, ktot, nfluid_IB
   real, intent(in) :: xgrid(itot), ygrid(jtot), zgrid(ktot)
   logical, dimension(itot,jtot,ktot), intent(out) :: fluid_IB, solid_IB
   integer, dimension(nfluid_IB,3) :: fluid_IB_ijk
   real, dimension(nfluid_IB,3), intent(out) :: fluid_IB_xyz
   integer :: fluid_IB_read, solid_IB_read, i, j, k, n
   character(80) :: chmess
   integer, parameter :: ifinput = 1

   open (ifinput, file=fname_fluid_boundary)
   read (ifinput, '(a80)') chmess
   do n = 1,nfluid_IB
     read (ifinput, *) fluid_IB_ijk(n,1), fluid_IB_ijk(n,2), fluid_IB_ijk(n,3)
     fluid_IB_xyz(n,:) = (/xgrid(fluid_IB_ijk(n,1)), ygrid(fluid_IB_ijk(n,2)), zgrid(fluid_IB_ijk(n,3))/)
     write(*,*) fluid_IB_ijk(n,:)
     write(*,*) fluid_IB_xyz(n,:)
   end do
   close (ifinput)

   write(*,*) "fluid_IB_xyz", fluid_IB_xyz

   open (ifinput, file='fluid_IB_c.txt')
   do j = 1,jtot
      do k = 1,ktot
           do i = 1,itot
             read (ifinput, *) fluid_IB_read
             if (fluid_IB_read == 1) then
                 fluid_IB(i,j,k) = .true.
             else
                 fluid_IB(i,j,k) = .false.
             end if
           end do
      end do
   end do
   close (ifinput)

   open (ifinput, file='solid_IB_c.txt')
   do j = 1,jtot
      do k = 1,ktot
           do i = 1,itot
             read (ifinput, *) solid_IB_read
             if (solid_IB_read == 1) then
                 solid_IB(i,j,k) = .true.
             else
                 solid_IB(i,j,k) = .false.
             end if
           end do
      end do
   end do
   close (ifinput)


   end subroutine readBoundaryPoints


   subroutine matchFacetsToCells(faces, faceNormal, nFaces, vertices, nVertices, &
      fluid_IB, solid_IB, fluid_IB_xyz, nfluid_IB, &
      xgrid, ygrid, zgrid, itot, jtot, ktot, &
      diag_neighbs, periodic_x, periodic_y, &
      secfacids, secbndptids, secareas, bnddst)

    use, intrinsic :: ieee_arithmetic
    integer, intent(in) :: nFaces, nVertices, nfluid_IB, itot, jtot, ktot
    integer, intent(in) :: faces(nFaces,3)
    real   , intent(in) :: faceNormal(nFaces,3), vertices(nVertices,3), fluid_IB_xyz(nfluid_IB,3)
    real   , intent(in) :: xgrid(itot), ygrid(jtot), zgrid(ktot)
    logical, intent(in), dimension(itot,jtot,ktot) :: fluid_IB, solid_IB
    logical, intent(in) :: diag_neighbs, periodic_x, periodic_y
    integer, dimension(:), allocatable, intent(out) :: secfacids, secbndptids
    real, dimension(:), allocatable, intent(out) :: secareas, bnddst
    integer, dimension(:,:), allocatable :: clipFaces!, fluid_IB_ijk
    !logical, dimension(:,:,:), allocatable :: fluid_IB, solid_IB
    !real, dimension(:,:), allocatable :: incenter, faceNormal, vertices, fluid_IB_xyz, solid_IB_xyz
    real, dimension(:,:), allocatable :: clipVertices, projVert, clipFaceNormal
   ! real, dimension(:), allocatable :: xgrid, ygrid, zgrid
    logical :: il_comp(itot+1), iu_comp(itot+1), jl_comp(jtot+1), ju_comp(jtot+1), kl_comp(ktot+1), ku_comp(ktot+1)
    logical :: search_adj
    integer :: il, iu, jl, ju, kl, ku, nClipFaces, nClipVertices, id, loc
    real :: dx, dy, dz, xmin, xmax, ymin, ymax, zmin, zmax, xl, xu, yl, yu, zl, zu, tol, &
            planes(6,4), area, planeNormal(3), xproj, yproj, zproj, proj, projVec(3), projArea, &
            dist, dists(27), angle, angles(27), BI(3), BIs(27,3)
    real, dimension(3) :: xyz, xyz1, xyz2, xyz3, xyz4, xyz5, xyz6, xyz7, xyz8, xyz9, &
                          xyz10, xyz11, xyz12, xyz13, xyz14, xyz15, xyz16, xyz17, xyz18, &
                          xyz19, xyz20, xyz21, xyz22, xyz23, xyz24, xyz25, xyz26, xyz27
    integer :: n, m, i, j, k, dir, ids(2)


    !allocate(il_comp(itot+1), iu_comp(jtot+1), jl_comp(jtot+1), ju_comp(jtot+1), kl_comp(ktot+1), ku_comp(ktot+1))

    !
    ! allocate(fluid_IB_ijk(nfluid_IB,3), fluid_IB_xyz(nfluid_IB,3), solid_IB_xyz(nsolid_IB,3))
    ! open (ifinput, file='fluid_boundary_c.txt')
    ! read (ifinput, '(a80)') chmess
    ! do n = 1,nfluid_IB
    !   read (ifinput, *) fluid_IB_ijk(n,1), fluid_IB_ijk(n,2), fluid_IB_ijk(n,3)
    !   fluid_IB_xyz(n,:) = (/xgrid(fluid_IB_ijk(n,1)), ygrid(fluid_IB_ijk(n,2)), zgrid(fluid_IB_ijk(n,3))/)
    ! end do
    ! close (ifinput)
    !
    ! write(*,*) "fluid_IB_xyz", fluid_IB_xyz
    !
    !
    ! allocate(fluid_IB(itot,jtot,ktot), solid_IB(itot,jtot,ktot))
    ! open (ifinput, file='fluid_IB_c.txt')
    ! count = 1
    ! do j = 1,jtot
    !     do k = 1,ktot
    !         do i = 1,itot
    !            read (ifinput, *) fluid_IB_read
    !            if (fluid_IB_read == 1) then
    !               fluid_IB(i,j,k) = .true.
    !            else
    !               fluid_IB(i,j,k) = .false.
    !            end if
    !            count = count+1
    !         end do
    !     end do
    ! end do
    ! close (ifinput)
    !
    ! open (ifinput, file='solid_IB_c.txt')
    ! count = 1
    ! do j = 1,jtot
    !     do k = 1,ktot
    !         do i = 1,itot
    !            read (ifinput, *) solid_IB_read
    !            if (solid_IB_read == 1) then
    !               solid_IB(i,j,k) = .true.
    !            else
    !               solid_IB(i,j,k) = .false.
    !            end if
    !            count = count+1
    !         end do
    !     end do
    ! end do
    ! close (ifinput)

    ! nFaces = 10
    ! nVertices = 12
    ! allocate(faces(nFaces, 3))
    ! allocate(vertices(nVertices, 3))
    ! allocate(faceNormal(nFaces, 3))
    ! allocate(incenter(nFaces, 3))
    !
    ! open (ifinput, file='faces.txt')
    ! do n = 1,nFaces
    !   read (ifinput, *) faces(n,1), faces(n,2), faces(n,3)
    ! end do
    ! close (ifinput)
    !
    !
    ! open (ifinput, file='faceNormals.txt')
    ! do n = 1,nFaces
    !   read (ifinput, *) faceNormal(n,1), faceNormal(n,2), faceNormal(n,3)
    ! end do
    ! close (ifinput)
    !
    ! open (ifinput, file='vertices.txt')
    ! do n = 1,nVertices
    !   read (ifinput, *) vertices(n,1), vertices(n,2), vertices(n,3)
    ! end do
    ! close (ifinput)

    allocate(secareas(0))
    allocate(secfacids(0))
    allocate(secbndptids(0))
    allocate(bnddst(0))

    dx = xgrid(2)-xgrid(1)
    dy = ygrid(2)-ygrid(1)
    dz = zgrid(2)-zgrid(1)

    do n=1,nFaces
      write(*,*) "facet", n
      ! no shear stress in normal direction
      if ((xgrid(1) == 0. .and. all(abs(abs(faceNormal(n, :)) - (/1.,0.,0./)) < 1e-8)) .or. &
          (ygrid(1) == 0. .and. all(abs(abs(faceNormal(n, :)) - (/0.,1.,0./)) < 1e-8)) .or. &
          (zgrid(1) == 0. .and. all(abs(abs(faceNormal(n, :)) - (/0.,0.,1./)) < 1e-8))) cycle

      xmin = minval(vertices(faces(n,:),1))
      ymin = minval(vertices(faces(n,:),2))
      zmin = minval(vertices(faces(n,:),3))
      xmax = maxval(vertices(faces(n,:),1))
      ymax = maxval(vertices(faces(n,:),2))
      zmax = maxval(vertices(faces(n,:),3))

      write(*,*) "xmin, xmax, ymin, ymax, zmin, zmax", xmin, xmax, ymin, ymax, zmin, zmax

      ! ignore facets on the ground and facing down
      if ((abs(zmin) < epsilon(zmin) .and. abs(zmax) < epsilon(zmax)) .and. &
       all(abs(faceNormal(n, :) - (/0.,0.,-1./)) < epsilon(faceNormal(n, 1)))) cycle

       write(*,*) "here"

      tol = 1e-10
      where (xmin >= (/xgrid   -dx/2, xgrid(itot)+dx/2/)-tol)
         il_comp = .true.
      elsewhere
         il_comp = .false.
      end where
      where (xmax <= (/xgrid(1)-dx/2, xgrid      +dx/2/)+tol)
         iu_comp = .true.
      elsewhere
         iu_comp = .false.
      end where
      where (ymin >= (/ygrid   -dy/2, ygrid(jtot)+dy/2/)-tol)
         jl_comp = .true.
      elsewhere
         jl_comp = .false.
      end where
      where (ymax <= (/ygrid(1)-dy/2, ygrid      +dy/2/)+tol)
         ju_comp = .true.
      elsewhere
         ju_comp = .false.
      end where
      where (zmin >= (/zgrid   -dz/2, zgrid(ktot)+dz/2/)-tol)
         kl_comp = .true.
      elsewhere
         kl_comp = .false.
      end where
      where (zmax <= (/zgrid(1)-dz/2, zgrid      +dz/2/)+tol)
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

      if ((il==0) .or. (iu==0)) write(*,*) "problem with x coord"
      if ((jl==0) .or. (ju==0)) write(*,*) "problem with y coord"
      if ((kl==0) .or. (ku==0)) write(*,*) "problem with z coord"

      ! write(*,*) il_comp, il
      ! write(*,*) iu_comp, iu
      ! write(*,*) jl_comp, jl
      ! write(*,*) ju_comp, ju
      ! write(*,*) kl_comp, kl
      ! write(*,*) ku_comp, ku

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
               write(*,*) "i,j,k", i, j, k
               ! Define corners of cube
               xl = xgrid(i) - dx/2
               xu = xgrid(i) + dx/2
               yl = ygrid(j) - dy/2
               yu = ygrid(j) + dy/2
               zl = zgrid(k) - dz/2
               zu = zgrid(k) + dz/2

               planes(1,:) = (/ 1., 0., 0., xu/)
               planes(2,:) = (/-1., 0., 0.,-xl/)
               planes(3,:) = (/ 0., 1., 0., yu/)
               planes(4,:) = (/ 0.,-1., 0.,-yl/)
               planes(5,:) = (/ 0., 0., 1., zu/)
               planes(6,:) = (/ 0., 0.,-1.,-zl/)

               call sutherlandHodgman3D(vertices(faces(n,:),:), 3, planes, 6, clipVertices)
               nClipVertices = size(clipVertices, 1)
               if (nClipVertices == 0) then
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
                  projvec = (/xproj, yproj, zproj/)
                  dir = maxloc(abs(projvec), 1)
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

                  projArea = polyarea(projVert(:,ids(1)), projVert(:,ids(2)), nClipVertices)
                  deallocate(projVert)

                  area = projArea / proj

                  if (area < 1e-5) cycle

                  if (nClipVertices==4) then
                     nClipFaces = 2
                     allocate(clipFaces(nClipFaces,3))
                     clipFaces(1,:) = (/1, 2, 3/)
                     clipFaces(2,:) = (/1, 3, 4/)

                  elseif (nClipVertices==5) then
                     nClipFaces = 3
                     allocate(clipFaces(nClipFaces,3))
                     clipFaces(1,:) = (/1, 2, 3/)
                     clipFaces(2,:) = (/1, 3, 4/)
                     clipFaces(3,:) = (/1, 4, 5/)

                  elseif (nClipVertices==6) then
                     nClipFaces = 4
                     allocate(clipFaces(nClipFaces,3))
                     clipFaces(1,:) = (/1, 2, 3/)
                     clipFaces(2,:) = (/1, 3, 4/)
                     clipFaces(3,:) = (/1, 4, 5/)
                     clipFaces(4,:) = (/1, 5, 6/)

                  else
                     write(*,*) "size of clipped polygon greater than 6"
                  end if

               else
                 write(*,*) "something wrong with clipped polygon"
               end if

               !write(*,*) "clipFaces", clipFaces
               !write(*,*) "area", area

               if (((xgrid(i) == 0.) .and. periodic_x) .or. ((ygrid(j) == 0.) .and. periodic_y)) then ! Account for periodicity - flux at point N+1
                  area = area * 2.
               end if

               call appendToArray1D_real(secareas, area)
               call appendToArray1D_integer(secfacids, n)

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
                     id = 1 ! not necessary?
                     xyz = xyz1
                     write(*,*) "normal found"
                  else
                     write(*,*) "normal not found, searching adjacent cells"
                     ! Not normal, search adjacent fluid IB cells
                     search_adj = .true.
                     if (dist > 0) then
                        ! Include in comparison
                        dists(1) = dist
                        angles(1) = angle
                        BIs(1,:) = BI
                     end if
                  end if
                  !write(*,*) "dist", dist
               end if

               if (solid_IB(i,j,k) .or. search_adj) then
                   if (solid_IB(i,j,k)) then
                      write(*,*) "solid_IB"
                  !    xyz1 = (/xgrid(i), ygrid(j), zgrid(k)/)
                  !    loc = findloc(ismember_rows(fluid_IB_xyz, xyz1), .true., 1)
                   end if
                   if (search_adj) then
                      write(*,*) "search_adj"
                   end if


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
                           dists(17) = dist;
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
                           dists(22) = dist;
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
                           dists(25) = dist;
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
                           dists(27) = dist;
                           angles(27) = dot_product(faceNormal(n,:), (xyz27 - BI)/norm2(xyz27 - BI))
                           BIs(27,:) = BI
                        end if
                     end if

                  end if ! diag_neighbs

                  id = maxloc(abs(angles) / (dists / (dx*dy*dz)**(1./3.)), 1)
                  dist = dists(id)
                  BI = BIs(id,:)

                  if (isnan(dist)) then
                     write(*,*) "facet ", n, " in cell ", i, j, k, " could not find a cell to give flux to"
                     write(*,*) "ensure diag_neighbs = true"
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

               loc = findloc(ismember_rows(fluid_IB_xyz, xyz), .true., 1)
               call appendToArray1D_integer(secbndptids, loc)
               call appendToArray1D_real(bnddst, abs(dist))

               deallocate(clipFaces)
               deallocate(clipFaceNormal)
            end do
         end do
      end do

    end do

    write(*,*) "secfacids ", "secareas ", "secbndptids ", "bnddst"
    do n=1,size(secfacids,1)
      write(*,*) secfacids(n), secareas(n), secbndptids(n), bnddst(n)
   end do

end subroutine matchFacetsToCells

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
   logical :: check

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

      ! if all(isnan(distances))
      loc = minloc(distances, 1)
      minDistance = distances(loc)
      interceptPoint = projectedPoints(loc, :)

   end do

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
        allocate(inputList(size(tempPolygon, 1), 3))
        inputList = tempPolygon

        deallocate(tempPolygon)
        allocate(tempPolygon(0, 3))
        if (size(inputList, 1) > 0) then
            previousVertex = inputList(size(inputList, 1), :)
        end if

        do j = 1, size(inputList, 1)
            if (inside(inputList(j, :), clipPlanes(i, :))) then
                if (.not. inside(previousVertex, clipPlanes(i, :))) then
                    intersection = computeIntersection(clipPlanes(i, :), previousVertex, inputList(j, :))
                    call appendToArray2D(tempPolygon, intersection)
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
                call appendToArray2D(tempPolygon, inputList(j, :))
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
                intersection = computeIntersection(clipPlanes(i, :), previousVertex, inputList(j, :))
                call appendToArray2D(tempPolygon, intersection)
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


end program run
