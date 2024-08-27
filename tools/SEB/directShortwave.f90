program run
   implicit none

   integer :: nFaces, nVertices
   real :: irradiance, resolution
   real, dimension(3) :: nsun
   real   , allocatable, dimension(:,:)   :: vertices, faceNormal, incenter
   integer, allocatable, dimension(:,:)   :: connectivityList
   real   , allocatable, dimension(:)     :: Sdir
   integer :: n

   open(unit=50,file='info_directShortwave.txt')
   read(unit=50,fmt='(i8,x,i8)') nFaces, nVertices
   read(unit=50,fmt='(f15.10,x,f15.10,f15.10)') nsun(1), nsun(2), nsun(3)
   read(unit=50,fmt='(f15.10)') irradiance
   read(unit=50,fmt='(f15.10)') resolution
   close(unit=50)

   ! write(*,*) "nFaces, nVertices", nFaces, nVertices
   ! write(*,*) "nsun", nsun
   ! write(*,*) "irradiance", irradiance
   ! write(*,*) "resolution", resolution

   !resolution = 1

   allocate(connectivityList(nFaces, 3))
   allocate(incenter(nFaces, 3))
   allocate(faceNormal(nFaces, 3))
   allocate(vertices(nVertices, 3))
   allocate(Sdir(nFaces))

   call readGeometry('faces.txt', nFaces, 'vertices.txt', nVertices, &
     connectivityList, incenter, faceNormal, vertices)

   ! write(*,*) "faces"
   ! do n=1,nFaces
   !    write(*,*) n, connectivityList(n,:), incenter(n,:), faceNormal(n,:)
   ! end do
   !
   ! write(*,*) "vertices"
   ! do n=1,nVertices
   !    write(*,*) n, vertices(n,:)
   ! end do

   nsun = nsun / norm2(nsun) ! make sure the vector to the sun is normalized
   call calculateDirectShortwave(connectivityList, incenter, faceNormal, nFaces, vertices, nVertices, &
      nsun, irradiance, resolution, Sdir)

   call writeDirectShortwave(Sdir, nFaces)


contains

   subroutine readGeometry(fname_faces, nFaces, fname_vertices, nVertices, &
                           connectivityList, incenter, faceNormal, vertices)
     character(9), intent(in) :: fname_faces
     character(12), intent(in) :: fname_vertices
     integer, intent(in) :: nFaces, nVertices
     integer, intent(out) :: connectivityList(nFaces,3)
     real   , intent(out) :: faceNormal(nFaces,3), incenter(nFaces,3), vertices(nVertices,3)
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


   subroutine calculateDirectShortwave(connectivityList, incenter, faceNormal, nFaces, vertices, nVertices, &
      nsun, irradiance, resolution, Sdir)
      integer, intent(in) :: nFaces, nVertices
      integer, intent(in), dimension(nFaces,3)    :: connectivityList
      real   , intent(in), dimension(nFaces,3)    :: incenter, faceNormal
      real   , intent(in), dimension(nVertices,3) :: vertices
      real   , intent(in) :: nsun(3), irradiance, resolution ! nsun must be normalized
      real   , intent(out), dimension(nFaces) :: Sdir
      !real   , dimension(nFaces) :: Sdir
      real   , dimension(nFaces,3) :: planeIncenter
      real   , dimension(nFaces) :: areas, distIncenter, projAreas
      integer, dimension(nFaces) :: visibility, sortedFaces ! if visiblity is false then face is self-shaded, i.e. not visible to the sun
      real   , dimension(nVertices,3) :: planeVertices, projVertices
      real   , dimension(nVertices,2) :: locVertices
      real   , dimension(nVertices) :: distVertices
      logical, dimension(:,:), allocatable :: mask
      integer, dimension(:,:), allocatable :: maskIDs
      real   , dimension(:), allocatable :: locCoord1, locCoord2
      integer, dimension(:), allocatable :: counts
      real :: xmin, xmax, xrange, ymin, ymax, yrange, temp
      real, dimension(3) :: p0, u1, u2
      real, dimension(2) :: cor1, cor2, cor3, cor4
      real, dimension(3,3) :: matrix, invMatrix
      integer :: i, j, n, m, id_temp, size_xi, size_eta
      logical :: flag
      real :: start, finish

      ! projection
      xmin = minval(vertices(:,1))
      xmax = maxval(vertices(:,1))
      xrange = xmax - xmin
      ymin = minval(vertices(:,2))
      ymax = maxval(vertices(:,2))
      yrange = ymax - ymin
      p0 = (/xmin+xrange/2., ymin+yrange/2., 0./) + 3.*max(xrange,yrange) * nsun
      ! write(*,*) "p0", p0

      ! define orthonormal basis on the plane
      u1 = -(/nsun(2), -nsun(1), 0./)
      u1 = u1 / norm2(u1)
      u2 = cross_product(u1, nsun)

      ! write(*,*) "u1", u1
      ! write(*,*) "u2", u2

      ! M(1,:) = u1
      ! M(2,:) = u2
      ! M(3,:) = -nsun
      ! Note difference to MATLAB here - matrix has been transposed because need to left-multiply inverse
      matrix(:,1) = u1
      matrix(:,2) = u2
      matrix(:,3) = -nsun

      ! write(*,*) "M", matrix
      ! write(*,*) matrix(1,:)
      ! write(*,*) matrix(2,:)
      ! write(*,*) matrix(3,:)

      invMatrix = matinv3(matrix)

      ! write(*,*) "matmul(M, Minv)", matmul(matrix, invMatrix)
      !
      ! write(*,*) "Minv"
      ! write(*,*) invMatrix(1,:)
      ! write(*,*) invMatrix(2,:)
      ! write(*,*) invMatrix(3,:)
            ! solve matrix system M * projIncenter = (incenter - p0)

      ! calculate areas and determine visiblity

      sortedFaces(1) = 1

      do n=1,nFaces
         areas(n) = 0.5*norm2(cross_product(vertices(connectivityList(n,2),:) - vertices(connectivityList(n,1),:), &
                                            vertices(connectivityList(n,3),:) - vertices(connectivityList(n,1),:)))
         if (dot_product(faceNormal(n,:), nsun) > 0)  visibility(n) = 1
         planeIncenter(n,:) = matmul(invMatrix, incenter(n,:) - p0)
         distIncenter(n) = planeIncenter(n,3)
         !write(*,*) n, distIncenter(n)
         ! sort
         if (n > 1) then
            m = n
            do while (m > 1 .and. distIncenter(m) > distIncenter(m - 1)) ! note in descending order (>)
            ! Swap array(m) and array(m - 1)
            temp = distIncenter(m)
            distIncenter(m) = distIncenter(m - 1)
            distIncenter(m - 1) = temp
            ! Update the corresponding index
            id_temp = sortedFaces(m)
            sortedFaces(m) = sortedFaces(m - 1)
            sortedFaces(m - 1) = id_temp
            m = m - 1
            end do
            ! Update the index for the newly inserted element
            sortedFaces(m) = n
         end if
      end do

      ! open (unit=11,file='sortedFaces_fort.txt',action="write")
      ! do n=1,nFaces
      !    !write (10,*) secfacids(n), secareas(n), secbndptids(n), bnddst(n)
      !    ! Formatting assumes: #facets < 10 million, #fluid boundary points < 1 billion,
      !    ! section area < 1000 m^2 (rounded to cm^2), and distance < 1000m
      !    write(unit=11,fmt='(i4)') sortedFaces(n)
      ! end do
      ! close (11)

      ! do n=1,nFaces
      !    write(*,*) sortedFaces(n), distIncenter(n)
      ! end do

      do n=1,nVertices
         planeVertices(n,:) = matmul(invMatrix, vertices(n,:) - p0)
         distVertices(n) = planeVertices(n,3)
         projVertices(n,:) = vertices(n,:) + distVertices(n) * nsun
         !write(*,*) n, projVertices(n,:)
      end do

      locVertices(:,1) = planeVertices(:,1) + abs(minval(planeVertices(:,1)))
      locVertices(:,2) = planeVertices(:,2) + abs(minval(planeVertices(:,2)))

      ! do n=1,nVertices
      !    write(*,*) n, locVertices(n,:)
      ! end do

      ! construct mask
      size_xi = ceiling(maxval(locVertices(:,1) / resolution))
      size_eta = ceiling(maxval(locVertices(:,2) / resolution))
      !write(*,*) "sizeMask1, sizeMask2", sizeMask1, sizeMask2
      !allocate(mask(size_eta, size_xi))
      allocate(maskIDs(size_eta, size_xi))

      ! allocate(locCoord1(sizeMask1))
      ! allocate(locCoord2(sizeMask2))
      !
      ! ! cell centres
      ! do i=1,sizeMask1
      !    locCoord1(i) = resolution/2. + resolution*i
      !    !write(*,*) locCoord1(i)
      ! end do
      !
      ! do j=1,sizeMask2
      !    locCoord2(j) = resolution/2. + resolution*j
      !    !write(*,*) locCoord2(j)
      ! end do

      ! ! polygon scan conversion
      ! do n=1,nFaces
      !   write(*,*) "n", n
      !    m = sortedFaces(n)
      !    cor1 = locVertices(connectivityList(m, 1), :)
      !    cor2 = locVertices(connectivityList(m, 2), :)
      !    cor3 = locVertices(connectivityList(m, 3), :)
      !    mask = .false.
      !
      !    ! poly2mask
      !    do i=1,sizeMask1
      !       do j=1,sizeMask2
      !          cor4 = (/locCoord1(i), locCoord2(j)/)
      !          if (isInsideTriangle(cor1, cor2, cor3, cor4)) then
      !              !mask(i,j) = .true.
      !              maskIDs(i,j) = m * visibility(m)
      !           end if
      !       end do
      !    end do
      !
      !    !where (mask) maskIDs = m * visibility(m)
      ! end do

         ! do i=1,sizeMask1
         !    write(*,*) "i", i
         !    do j=1,sizeMask2
         !       n = 0
         !       flag = .true.
         !       do while(flag .and. n < nFaces) ! looping closest facets first
         !          n = n+1
         !          m = sortedFaces(n)
         !          cor1 = locVertices(connectivityList(m, 1), :)
         !          cor2 = locVertices(connectivityList(m, 2), :)
         !          cor3 = locVertices(connectivityList(m, 3), :)
         !          cor4 = (/locCoord1(i), locCoord2(j)/)
         !          if (isInsideTriangle(cor1, cor2, cor3, cor4)) then
         !             maskIDs(i,j) = m * visibility(m)
         !             flag = .false. ! once facet found, no need to loop any more
         !          end if
         !       end do
         !    end do
         ! end do

         call cpu_time(start)

         ! polygon scan conversion
         do n=1,nFaces
           !write(*,*) "n", n
            m = sortedFaces(n)
            !mask = .false.
            call poly2maskIDs(locVertices(connectivityList(m, :), 1) / resolution, &
                           locVertices(connectivityList(m, :), 2) / resolution, size_eta, &
                           size_xi, maskIDs, m*visibility(m))

            ! call poly2maskIDs([locVertices(connectivityList(m, :), 1) / resolution, &
            !                    locVertices(connectivityList(m, 1), 1) / resolution], &
            !                   [locVertices(connectivityList(m, :), 2) / resolution, &
            !                    locVertices(connectivityList(m, 1), 2) / resolution], &
            !                   size_eta, size_xi, mask, maskIDs, m*visibility(m));

            !write(*,*) mask
            !maskIDs(i,j) = m * visibility(m)

            !where (mask) maskIDs = m * visibility(m)
            !write(*,*) maskIDs
         end do

         !write(*,*) shape(mask)

      allocate(counts(nFaces))
      counts = 0
      do i=1,size(maskIDs,1)
         do j=1,size(maskIDs,2)
            if (maskIDs(i,j) > 0) counts(maskIDs(i,j)) = counts(maskIDs(i,j)) + 1
         end do
      end do

      ! open (unit=11,file='counts_fort.txt',action="write")
      ! do n=1,nFaces
      !    !write (10,*) secfacids(n), secareas(n), secbndptids(n), bnddst(n)
      !    ! Formatting assumes: #facets < 10 million, #fluid boundary points < 1 billion,
      !    ! section area < 1000 m^2 (rounded to cm^2), and distance < 1000m
      !    write(unit=11,fmt='(i8)') counts(n)
      ! end do
      ! close (11)

            call cpu_time(finish)

      projAreas = counts * resolution**2

      Sdir = irradiance * projAreas / areas

      ! write(*,*) "n", "count", "projArea", "Sdir"
      ! do n=1,1000
      !    write(*,*) n, counts(n), projAreas(n), Sdir(n)
      ! end do

      print '("Time = ",f10.3," seconds.")',finish-start

   end subroutine calculateDirectShortwave


   subroutine writeDirectShortwave(Sdir, nFaces)
      integer, intent(in) :: nFaces
      real   , intent(in), dimension(nFaces) :: Sdir
      !character(25), intent(in) :: fname_Sdir
      integer :: n

      open (unit=10,file='Sdir.txt',action="write")
      do n=1,nFaces
         write(unit=10,fmt='(f8.2)') Sdir(n)
      end do
      close (10)

      write(*,*) "written Sdir.txt"

   end subroutine writeDirectShortwave


   ! logical function isInsideTriangle(A, B, C, P)
   !    real, intent(in), dimension(2) :: A, B, C, P
   !    real, dimension(2) :: v0, v1, v2
   !    real :: u, v, dot00, dot01, dot02, dot11, dot12, invDenom
   !
   !    v0 = C - A
   !    v1 = B - A
   !    v2 = P - A
   !
   !    dot00 = dot_product(v0, v0)
   !    dot01 = dot_product(v0, v1)
   !    dot02 = dot_product(v0, v2)
   !    dot11 = dot_product(v1, v1)
   !    dot12 = dot_product(v1, v2)
   !
   !    invDenom = 1. / (dot00 * dot11 - dot01 * dot01)
   !    u = (dot11 * dot02 - dot01 * dot12) * invDenom
   !    v = (dot00 * dot12 - dot01 * dot02) * invDenom
   !    isInsideTriangle = .false.
   !    if ((u >= 0) .and. (v >= 0) .and. (u + v < 1)) isInsideTriangle = .true.
   !
   ! end function isInsideTriangle


   function cross_product(a,b)
     ! Calculate the cross product (a x b)
     implicit none
     real, dimension(3) :: cross_product
     real, dimension(3), intent(in) :: a, b

     cross_product(1) = a(2)*b(3) - a(3)*b(2)
     cross_product(2) = a(3)*b(1) - a(1)*b(3)
     cross_product(3) = a(1)*b(2) - a(2)*b(1)

   end function cross_product


   function matinv3(A) result(B)
     !! calculates the inverse of a 3×3 matrix.
     real, intent(in) :: A(3, 3) !! matrix
     real             :: B(3, 3) !! inverse matrix
     real             :: detinv

     !inverse determinant of the matrix
     detinv = 1/(A(1, 1)*A(2, 2)*A(3, 3) - A(1, 1)*A(2, 3)*A(3, 2) &
     - A(1, 2)*A(2, 1)*A(3, 3) + A(1, 2)*A(2, 3)*A(3, 1) &
     + A(1, 3)*A(2, 1)*A(3, 2) - A(1, 3)*A(2, 2)*A(3, 1))

     !inverse of the matrix
     B(1, 1) = +detinv*(A(2, 2)*A(3, 3) - A(2, 3)*A(3, 2))
     B(2, 1) = -detinv*(A(2, 1)*A(3, 3) - A(2, 3)*A(3, 1))
     B(3, 1) = +detinv*(A(2, 1)*A(3, 2) - A(2, 2)*A(3, 1))
     B(1, 2) = -detinv*(A(1, 2)*A(3, 3) - A(1, 3)*A(3, 2))
     B(2, 2) = +detinv*(A(1, 1)*A(3, 3) - A(1, 3)*A(3, 1))
     B(3, 2) = -detinv*(A(1, 1)*A(3, 2) - A(1, 2)*A(3, 1))
     B(1, 3) = +detinv*(A(1, 2)*A(2, 3) - A(1, 3)*A(2, 2))
     B(2, 3) = -detinv*(A(1, 1)*A(2, 3) - A(1, 3)*A(2, 1))
     B(3, 3) = +detinv*(A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1))
   end function


  !  subroutine poly2mask(x, y, M, N, BW)
  ! implicit none
  ! integer, intent(in) :: M, N
  ! real(8), intent(in) :: x(:), y(:)
  ! logical, intent(out) :: BW(M, N)
  ! integer :: edgeRows, edgeColumns, newLength, i
  ! real(8) :: xNew(:), yNew(:)
  ! real(8) :: xpt(:), ypt(:)
  ! integer :: scale
  ! integer :: sizeX
  ! real(8), dimension(:), allocatable :: xLinePts, yLinePts
  ! integer :: borderSize
  ! integer :: xLength
  ! integer, dimension(N) :: minY, maxY
  ! integer :: xUse, yUse
  ! integer :: pt, diff
  !
  ! ! Check that the number of inputs is 4
  ! if (size(x) /= size(y)) then
  !   write(*,*) "Error: x and y vectors must have the same size."
  !   return
  ! end if
  !
  ! ! Initialize xNew and yNew
  ! edgeRows = size(x)
  ! edgeColumns = size(y)
  ! newLength = edgeRows + 1
  ! allocate(xNew(newLength))
  ! allocate(yNew(newLength))
  !
  ! ! Check if the last point is different from the first
  ! if (x(edgeRows) /= x(1) .or. y(edgeRows) /= y(1)) then
  !   xNew(1:edgeRows) = x(1:edgeRows)
  !   xNew(edgeRows + 1) = x(1)
  !   yNew(1:edgeRows) = y(1:edgeRows)
  !   yNew(edgeRows + 1) = y(1)
  !   call poly2MaskFromIntegralImagePortable(xNew, yNew, M, N, BW)
  ! else
  !   call poly2MaskFromIntegralImagePortable(x, y, M, N, BW)
  ! end if
  !
  ! ! Deallocate arrays
  ! deallocate(xNew)
  ! deallocate(yNew)
  !
  ! end subroutine poly2mask

  ! subroutine poly2maskIDs(xpt, ypt, M, N, out, outIDs, id)
  !     real  , intent(in) :: xpt(:), ypt(:) ! assumes x(end) == x(1)
  !     real, allocatable, dimension(:) :: x, y
  !     integer, intent(in) :: M, N, id
  !     logical, intent(out) :: out(M, N)
  !     integer, intent(out) :: outIDs(M, N)
  !     integer :: sizeX
  !     real    :: scale = 5.0
  !     integer :: minY(N), maxY(N)
  !
  !     ! Initialize minY and maxY arrays
  !     minY = 0
  !     maxY = 0
  !     sizeX = size(xpt,1) ! number of vertices
  !     allocate(x(sizeX), y(sizeX))
  !     x = xpt
  !     y = ypt
  !     ! x(1:sizeX-1) = xpt
  !     ! y(1:sizeX-1) = ypt
  !     ! x(sizeX) = xpt(1)
  !     ! y(sizeX) = ypt(1)
  !
  !     ! Create edges in mask to be used during parity scan
  !     call poly2edgelist(x, y, sizeX, scale, M, N, out, minY, maxY)
  !     ! Perform the parity scan over the output mask 'out'
  !     call parityScan(out, N, minY, maxY, outIDs, id)
  !
  !     deallocate(x, y)
  !
  ! end subroutine poly2maskIDs

    !subroutine poly2maskIDs(xpt, ypt, M, N, out, outIDs, id)
    subroutine poly2maskIDs(xpt, ypt, M, N, outIDs, id)
        real  , intent(in) :: xpt(:), ypt(:) ! assumes x(end) != x(1)
        real, allocatable, dimension(:) :: x, y
        integer, intent(in) :: M, N, id
        !logical, intent(out) :: out(M, N)
        logical :: out(M, N)
        integer, intent(out) :: outIDs(M, N)
        integer :: sizeX
        real    :: scale = 5.0
        integer :: minY(N), maxY(N)

        ! Initialize minY and maxY arrays
        minY = 0
        maxY = 0
        sizeX = size(xpt,1)+1 ! number of vertices
        allocate(x(sizeX), y(sizeX))
        x(1:sizeX-1) = xpt
        y(1:sizeX-1) = ypt
        x(sizeX) = xpt(1)
        y(sizeX) = ypt(1)

        ! Create edges in mask to be used during parity scan
        call poly2edgelist(x, y, sizeX, scale, M, N, out, minY, maxY)
        ! Perform the parity scan over the output mask 'out'
        call parityScan(out, N, minY, maxY, outIDs, id)

        deallocate(x, y)

    end subroutine poly2maskIDs

    subroutine poly2edgelist(x, y, xLength, scale, M, N, out, minY, maxY)
        real  , intent(inout) :: x(:), y(:)
        integer, intent(in) :: xLength, M, N
        real   , intent(in) :: scale
        logical, intent(out) :: out(M, N)
        integer, intent(out) :: minY(N), maxY(N)
        real, allocatable, dimension(:) :: xLinePts, yLinePts
        integer :: borderSize
        integer :: i, pt, xUse, yUse
        integer :: tempDX, tempDY
        real :: diff, scaledDown

        borderSize = 0
        do i = 1, xLength
           x(i) = floor(scale * (x(i) - 0.5) + 0.5) + 1.0
           y(i) = floor(scale * (y(i) - 0.5) + 0.5) + 1.0
           if (i > 1) then
             tempDX = nint(abs(x(i) - x(i - 1)))
             tempDY = nint(abs(y(i) - y(i - 1)))
             if (tempDX >= tempDY) then
                borderSize = borderSize + tempDX + 1
             else
                borderSize = borderSize + tempDY + 1
             end if
           end if
        end do

        ! Fill xLinePts and yLinePts with border information
        allocate(xLinePts(borderSize), yLinePts(borderSize))
        call makeBorder(x, y, xLength, borderSize, xLinePts, yLinePts)

         do pt = 1, borderSize - 1
           ! Check if x changes
           diff = xLinePts(pt + 1) - xLinePts(pt)
           if (abs(diff) >= 1.0) then
             yLinePts(pt) = min(yLinePts(pt), yLinePts(pt + 1))
             if (diff < 0.0) then
                xLinePts(pt) = xLinePts(pt) - 1.0
             end if
             scaledDown = (xLinePts(pt) + (scale - 1.0) / 2.0) / scale
             if (abs(scaledDown - dble(floor(scaledDown))) < 1.0 / (scale * 50.0)) then
                xLinePts(pt) = scaledDown
                yLinePts(pt) = ceiling((yLinePts(pt) + (scale - 1.0) / 2.0) / scale)
                ! Set xUse and yUse to these scaled down x and y values for indexing into the boolean array.
                xUse = nint(xLinePts(pt)) + 1
                yUse = nint(yLinePts(pt)) + 1
               if (.not. (xUse < 2 .or. xUse > N + 1)) then
                 if (yUse > M + 1) then
                   if (minY(xUse - 1) == 0) then
                     minY(xUse - 1) = M + 1
                   else
                     minY(xUse - 1) = min(minY(xUse - 1), M + 1)
                   end if
                   if (maxY(xUse - 1) == 0) then
                     maxY(xUse - 1) = M + 1
                   else
                     maxY(xUse - 1) = max(maxY(xUse - 1), M + 1)
                    end if
               else
                  yUse = max(2, yUse)
                  out(yUse - 1, xUse - 1) = .not. out(yUse - 1, xUse - 1)
                  if (minY(xUse - 1) == 0) then
                     minY(xUse - 1) = yUse - 1
                  else
                     minY(xUse - 1) = min(minY(xUse - 1), yUse - 1)
                  end if
                  if (maxY(xUse - 1) == 0) then
                     maxY(xUse - 1) = yUse
                  else
                     maxY(xUse - 1) = max(maxY(xUse - 1), yUse)
                  end if
               end if
              end if

           end if
         end if
       end do

       deallocate(xLinePts, yLinePts)

    end subroutine poly2edgelist

    subroutine makeBorder(x, y, xLength, borderSize, xLinePts, yLinePts)
        real   , intent(in) :: x(:), y(:)
        integer, intent(in) :: xLength, borderSize
        real, intent(out) :: xLinePts(borderSize), yLinePts(borderSize)
        integer :: borderPosition
        integer :: i
        integer :: dx, dy
        real    :: m
        integer :: xit, yit
        real :: xVal, yVal
        real :: x1, y1, x2, y2

        borderPosition = 1

        do i = 2, xLength
           call intLine(x(i-1),y(i-1),x(i),y(i), xLinePts, yLinePts, borderPosition)
        end do
        !
        !    x1 = x(i-1)
        !    y1 = y(i-1)
        !    x2 = x(i)
        !    y2 = y(i)
        !     !call intLine(x(i-1), y(i-1), x(i), y(i), xLinePts, yLinePts, borderPosition)
        !     dx = nint(abs(x2 - x1))
        !     dy = nint(abs(y2 - y1))
        !
        !     if (dx == 0 .and. dy == 0) then
        !         xLinePts(borderPosition) = x1
        !         yLinePts(borderPosition) = y1
        !         borderPosition = borderPosition + 1
        !         cycle
        !     end if
        !
        !     if (dx >= dy) then
        !         m = (y2 - y1) / (x2 - x1)
        !         if (x2 > x1) then
        !             !do xVal = x1, x2
        !             do xit = int(x1), int(x2)
        !                 xVal = xit + x1-int(x1)
        !                 yVal = nint(y1 + m * (xVal - x1))
        !                 xLinePts(borderPosition) = xVal
        !                 yLinePts(borderPosition) = yVal
        !                 borderPosition = borderPosition + 1
        !             end do
        !         else
        !             !do xVal = x1, x2, -1
        !             do xit = int(x1), int(x2), -1
        !                 xVal = xit + x1-int(x1)
        !                 yVal = nint(y1 + m * (xVal - x1))
        !                 xLinePts(borderPosition) = xVal
        !                 yLinePts(borderPosition) = yVal
        !                 borderPosition = borderPosition + 1
        !             end do
        !         end if
        !     else
        !         m = (x2 - x1) / (y2 - y1)
        !         if (y2 > y1) then
        !             !do yVal = y1, y2
        !             do yit = int(y1), int(y2)
        !                 yVal = yit + y1-int(y1)
        !                 xVal = nint(x1 + m * (yVal - y1))
        !                 xLinePts(borderPosition) = xVal
        !                 yLinePts(borderPosition) = yVal
        !                 borderPosition = borderPosition + 1
        !             end do
        !         else
        !            !do yVal = y1, y2, -1
        !            do yit = int(y1), int(y2), -1
        !                 yVal = yit + y1-int(y1)
        !                 xVal = nint(x1 + m * (yVal - y1))
        !                 xLinePts(borderPosition) = xVal
        !                 yLinePts(borderPosition) = yVal
        !                 borderPosition = borderPosition + 1
        !             end do
        !         end if
        !     end if
        ! end do

    end subroutine makeBorder

    subroutine intLine(x1, y1, x2, y2, xLinePts, yLinePts, borderPosition)
        real   , intent(in) :: x1, y1, x2, y2
        real   , intent(out) :: xLinePts(:), yLinePts(:)
        integer, intent(inout) :: borderPosition
        integer :: dx, dy
        real    :: m
        integer :: xit, yit
        real :: xVal, yVal

        dx = nint(abs(x2 - x1))
        dy = nint(abs(y2 - y1))

        if (dx == 0 .and. dy == 0) then
            xLinePts(borderPosition) = x1
            yLinePts(borderPosition) = y1
            borderPosition = borderPosition + 1
            return
        end if

        if (dx >= dy) then
            m = (y2 - y1) / (x2 - x1)
            if (x2 > x1) then
                !do xVal = x1, x2
                do xit = int(x1), int(x2)
                    xVal = xit + x1-int(x1)
                    yVal = nint(y1 + m * (xVal - x1))
                    xLinePts(borderPosition) = xVal
                    yLinePts(borderPosition) = yVal
                    borderPosition = borderPosition + 1
                end do
            else
                !do xVal = x1, x2, -1
                do xit = int(x1), int(x2), -1
                    xVal = xit + x1-int(x1)
                    yVal = nint(y1 + m * (xVal - x1))
                    xLinePts(borderPosition) = xVal
                    yLinePts(borderPosition) = yVal
                    borderPosition = borderPosition + 1
                end do
            end if
        else
            m = (x2 - x1) / (y2 - y1)
            if (y2 > y1) then
                !do yVal = y1, y2
                do yit = int(y1), int(y2)
                    yVal = yit + y1-int(y1)
                    xVal = nint(x1 + m * (yVal - y1))
                    xLinePts(borderPosition) = xVal
                    yLinePts(borderPosition) = yVal
                    borderPosition = borderPosition + 1
                end do
            else
               !do yVal = y1, y2, -1
               do yit = int(y1), int(y2), -1
                    yVal = yit + y1-int(y1)
                    xVal = nint(x1 + m * (yVal - y1))
                    xLinePts(borderPosition) = xVal
                    yLinePts(borderPosition) = yVal
                    borderPosition = borderPosition + 1
                end do
            end if
        end if
    end subroutine intLine

    subroutine parityScan(out, N, minY, maxY, outIDs, id)
        logical, intent(inout) :: out(:,:)
        integer, intent(out) :: outIDs(:,:)
        integer, intent(in) :: N, id
        integer, intent(in) :: minY(N), maxY(N)
        logical :: pixel
        integer :: c, r!, id_first, id_last

        do c = 1, N
            pixel = .false.
            do r = minY(c), maxY(c) - 1
                if (out(r, c)) pixel = .not. pixel
                out(r, c) = pixel
                if (out(r,c)) outIDs(r, c) = id
            end do
            ! if (maxY(c) > 0) then
            !    id_first = findloc(out(minY(c):maxY(c)-1,c), .true., 1) + minY(c) - 1
            !    id_last  = findloc(out(minY(c):maxY(c)-1,c), .true., 1, back=.true.) + minY(c) - 1
            !    out(id_last, c) = .false.
            !    if (id_first+1 <= id_last-1) then
            !       out(id_first+1:id_last-1, c) = .true.
            !    end if
            ! end if
        end do
    end subroutine parityScan



end program run
