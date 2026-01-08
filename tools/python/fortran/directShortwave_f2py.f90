module directshortwave_f2py_mod
   implicit none
contains

   subroutine calculate_direct_shortwave(connectivityList, incenter, faceNormal, vertices, nFaces, nVertices, nsun, irradiance, resolution, Sdir)
      integer, intent(in) :: nFaces, nVertices
      integer, intent(in), dimension(nFaces, 3) :: connectivityList
      real   , intent(in), dimension(nFaces, 3) :: incenter, faceNormal
      real   , intent(in), dimension(nVertices, 3) :: vertices
      real   , intent(in) :: nsun(3), irradiance, resolution
      real   , intent(out), dimension(nFaces) :: Sdir

      real   , dimension(:, :), allocatable :: planeIncenter
      real   , dimension(:), allocatable :: areas, distIncenter, projAreas
      integer, dimension(:), allocatable :: visibility, sortedFaces
      real   , dimension(:, :), allocatable :: planeVertices, projVertices
      real   , dimension(:, :), allocatable :: locVertices
      real   , dimension(:), allocatable :: distVertices
      integer, dimension(:, :), allocatable :: maskIDs
      real   , dimension(:), allocatable :: locCoord1, locCoord2
      integer, dimension(:), allocatable :: counts
      real :: xmin, xmax, xrange, ymin, ymax, yrange, temp
      real, dimension(3) :: p0, u1, u2, nsun_unit
      real, dimension(3, 3) :: matrix, invMatrix
      integer :: i, j, n, m, id_temp, size_xi, size_eta
      real :: start, finish

      allocate(planeIncenter(nFaces, 3))
      allocate(areas(nFaces))
      allocate(distIncenter(nFaces))
      allocate(projAreas(nFaces))
      allocate(visibility(nFaces))
      allocate(sortedFaces(nFaces))
      allocate(planeVertices(nVertices, 3))
      allocate(projVertices(nVertices, 3))
      allocate(locVertices(nVertices, 2))
      allocate(distVertices(nVertices))

      xmin = minval(vertices(:,1))
      xmax = maxval(vertices(:,1))
      xrange = xmax - xmin
      ymin = minval(vertices(:,2))
      ymax = maxval(vertices(:,2))
      yrange = ymax - ymin
      nsun_unit = nsun / norm2(nsun)
      p0 = (/xmin+xrange/2., ymin+yrange/2., 0./) + 3.*max(xrange,yrange) * nsun_unit

      u1 = -(/nsun_unit(2), -nsun_unit(1), 0./)
      u1 = u1 / norm2(u1)
      u2 = cross_product(u1, nsun_unit)

      matrix(:,1) = u1
      matrix(:,2) = u2
      matrix(:,3) = -nsun_unit

      invMatrix = matinv3(matrix)

      sortedFaces(1) = 1
      do n=1,nFaces
         areas(n) = 0.5*norm2(cross_product(vertices(connectivityList(n,2),:) - vertices(connectivityList(n,1),:), &
                                            vertices(connectivityList(n,3),:) - vertices(connectivityList(n,1),:)))
         visibility(n) = 0
         if (dot_product(faceNormal(n,:), nsun_unit) > 0)  visibility(n) = 1
         planeIncenter(n,:) = matmul(invMatrix, incenter(n,:) - p0)
         distIncenter(n) = planeIncenter(n,3)
         if (n > 1) then
            m = n
            do while (m > 1 .and. distIncenter(m) > distIncenter(m - 1))
               temp = distIncenter(m)
               distIncenter(m) = distIncenter(m - 1)
               distIncenter(m - 1) = temp
               id_temp = sortedFaces(m)
               sortedFaces(m) = sortedFaces(m - 1)
               sortedFaces(m - 1) = id_temp
               m = m - 1
            end do
            sortedFaces(m) = n
         end if
      end do

      do n=1,nVertices
         planeVertices(n,:) = matmul(invMatrix, vertices(n,:) - p0)
         distVertices(n) = planeVertices(n,3)
         projVertices(n,:) = vertices(n,:) + distVertices(n) * nsun_unit
      end do

      locVertices(:,1) = planeVertices(:,1) + abs(minval(planeVertices(:,1)))
      locVertices(:,2) = planeVertices(:,2) + abs(minval(planeVertices(:,2)))

      size_xi = ceiling(maxval(locVertices(:,1)) / resolution)
      size_eta = ceiling(maxval(locVertices(:,2)) / resolution)
      allocate(maskIDs(size_eta, size_xi))
      maskIDs = 0

      allocate(locCoord1(size_xi))
      allocate(locCoord2(size_eta))
      locCoord1 = (/ ( (i - 0.5) * resolution, i = 1, size_xi ) /)
      locCoord2 = (/ ( (i - 0.5) * resolution, i = 1, size_eta ) /)

      call cpu_time(start)
      do n=1,nFaces
         m = sortedFaces(n)
         call poly2maskIDs(locVertices(connectivityList(m, :), 1) / resolution, &
                           locVertices(connectivityList(m, :), 2) / resolution, size_eta, &
                           size_xi, maskIDs, m*visibility(m))
      end do

      allocate(counts(nFaces))
      counts = 0
      do i=1,size(maskIDs,1)
         do j=1,size(maskIDs,2)
            if (maskIDs(i,j) > 0) counts(maskIDs(i,j)) = counts(maskIDs(i,j)) + 1
         end do
      end do

      call cpu_time(finish)
      projAreas = counts * resolution**2
      Sdir = irradiance * projAreas / areas

      deallocate(planeIncenter, areas, distIncenter, projAreas, visibility, sortedFaces)
      deallocate(planeVertices, projVertices, locVertices, distVertices, maskIDs)
      deallocate(locCoord1, locCoord2, counts)

   end subroutine calculate_direct_shortwave


   function cross_product(a,b)
     implicit none
     real, dimension(3) :: cross_product
     real, dimension(3), intent(in) :: a, b

     cross_product(1) = a(2)*b(3) - a(3)*b(2)
     cross_product(2) = a(3)*b(1) - a(1)*b(3)
     cross_product(3) = a(1)*b(2) - a(2)*b(1)
   end function cross_product


   function matinv3(A) result(B)
     real, intent(in) :: A(3, 3)
     real             :: B(3, 3)
     real             :: detinv

     detinv = 1/(A(1, 1)*A(2, 2)*A(3, 3) - A(1, 1)*A(2, 3)*A(3, 2) &
     - A(1, 2)*A(2, 1)*A(3, 3) + A(1, 2)*A(2, 3)*A(3, 1) &
     + A(1, 3)*A(2, 1)*A(3, 2) - A(1, 3)*A(2, 2)*A(3, 1))

     B(1, 1) = +detinv*(A(2, 2)*A(3, 3) - A(2, 3)*A(3, 2))
     B(2, 1) = -detinv*(A(2, 1)*A(3, 3) - A(2, 3)*A(3, 1))
     B(3, 1) = +detinv*(A(2, 1)*A(3, 2) - A(2, 2)*A(3, 1))
     B(1, 2) = -detinv*(A(1, 2)*A(3, 3) - A(1, 3)*A(3, 2))
     B(2, 2) = +detinv*(A(1, 1)*A(3, 3) - A(1, 3)*A(3, 1))
     B(3, 2) = -detinv*(A(1, 1)*A(3, 2) - A(1, 2)*A(3, 1))
     B(1, 3) = +detinv*(A(1, 2)*A(2, 3) - A(1, 3)*A(2, 2))
     B(2, 3) = -detinv*(A(1, 1)*A(2, 3) - A(1, 3)*A(2, 1))
     B(3, 3) = +detinv*(A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1))
   end function matinv3


    subroutine poly2maskIDs(xpt, ypt, M, N, outIDs, id)
        real  , intent(in) :: xpt(:), ypt(:)
        real, allocatable, dimension(:) :: x, y
        integer, intent(in) :: M, N, id
        logical :: out(M, N)
        integer, intent(out) :: outIDs(M, N)
        integer :: sizeX
        real    :: scale = 5.0
        integer :: minY(N), maxY(N)

        minY = 0
        maxY = 0
        sizeX = size(xpt,1)+1
        allocate(x(sizeX), y(sizeX))
        x(1:sizeX-1) = xpt
        y(1:sizeX-1) = ypt
        x(sizeX) = xpt(1)
        y(sizeX) = ypt(1)

        call poly2edgelist(x, y, sizeX, scale, M, N, out, minY, maxY)
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

        allocate(xLinePts(borderSize), yLinePts(borderSize))
        call makeBorder(x, y, xLength, borderSize, xLinePts, yLinePts)

        do pt = 1, borderSize - 1
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

        borderPosition = 1
        do i = 2, xLength
           call intLine(x(i-1),y(i-1),x(i),y(i), xLinePts, yLinePts, borderPosition)
        end do

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
                do xit = int(x1), int(x2)
                    xVal = xit + x1-int(x1)
                    yVal = nint(y1 + m * (xVal - x1))
                    xLinePts(borderPosition) = xVal
                    yLinePts(borderPosition) = yVal
                    borderPosition = borderPosition + 1
                end do
            else
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
                do yit = int(y1), int(y2)
                    yVal = yit + y1-int(y1)
                    xVal = nint(x1 + m * (yVal - y1))
                    xLinePts(borderPosition) = xVal
                    yLinePts(borderPosition) = yVal
                    borderPosition = borderPosition + 1
                end do
            else
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
        integer :: c, r

        do c = 1, N
            pixel = .false.
            do r = minY(c), maxY(c) - 1
                if (out(r, c)) pixel = .not. pixel
                out(r, c) = pixel
                if (out(r,c)) outIDs(r, c) = id
            end do
        end do
    end subroutine parityScan

end module directshortwave_f2py_mod
