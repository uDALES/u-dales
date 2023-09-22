program matchFacetsToCells
   implicit none

    ! Declarations
    real, dimension(:, :), allocatable :: subjectPolygon, clippedPolygon
    real, dimension(:, :), allocatable :: clipPlanes
    integer :: nVertices, nPlanes, i, j

    ! Define the subject polygon (a 3D triangle)
    nVertices = 3
    allocate(subjectPolygon(nVertices, 3))
    subjectPolygon(1,:) = (/0.0, 0.0, 0.0/)
    subjectPolygon(2,:) = (/4.0, 0.0, 0.0/)
    subjectPolygon(3,:) = (/4.0, 24.0, 0.0/)

    ! Define clipping planes (to clip the triangle)
    nPlanes = 6
    allocate(clipPlanes(nPlanes, 4))
    clipPlanes(1,:) = (/1.0, 0.0, 0.0, 2.0/)
    clipPlanes(2,:) = (/-1.0, 0.0, 0.0, -1.0/)
    clipPlanes(3,:) = (/0.0, 1.0, 0.0, 6.0/)
    clipPlanes(4,:) = (/0.0, -1.0, 0.0, -5.0/)
    clipPlanes(5,:) = (/0.0, 0.0, 1.0, 1.0/)
    clipPlanes(6,:) = (/0.0, 0.0, -1.0, 0.0/)

    ! Call the polygon clipping subroutine
    call sutherlandHodgman3D(subjectPolygon, nVertices, clipPlanes, nPlanes, clippedPolygon)

    ! Print the clipped polygon
    print *, "Clipped Polygon:"
    do i = 1, size(clippedPolygon, 1)
        print *, clippedPolygon(i, :)
    end do

    ! Deallocate memory
    deallocate(subjectPolygon, clipPlanes, clippedPolygon)

contains

subroutine sutherlandHodgman3D(subjectPolygon, nVertices, clipPlanes, nPlanes, clippedPolygon)
    ! subjectPolygon: nVertices x 3
    ! clipPlanes: nPlanes x 4, where columns are n1, n2, n3, d (where (n1,n2,n3) dot (x,y,z) = d is the plane equation)
   integer, intent(in) :: nVertices, nPlanes
   real, intent(in) :: subjectPolygon(nVertices, 3), clipPlanes(nPlanes, 4)
   real, allocatable, intent(out) :: clippedPolygon(:,:)
   real, allocatable :: inputList(:,:), tempPolygon(:,:), catPolygon(:,:)
   real :: previousVertex(3), intersection(3)
   integer :: sizeCatPolygon

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
                    call concatenateToPolygon(tempPolygon, intersection)
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
                call concatenateToPolygon(tempPolygon, inputList(j, :))
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
                call concatenateToPolygon(tempPolygon, intersection)
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


subroutine concatenateToPolygon(polygon, point)
   !integer, intent(in) :: nVertices
   real, allocatable, intent(inout) :: polygon(:,:)
   real, dimension(3), intent(in) :: point
   real, dimension(:,:), allocatable :: catPolygon
   integer sizePolygon

   sizePolygon = size(polygon, 1)
   allocate(catPolygon(sizePolygon, 3))
   catPolygon = polygon
   deallocate(polygon)
   allocate(polygon(sizePolygon+1, 3))
   polygon(1:sizePolygon, :) = catPolygon
   polygon(sizePolygon+1, :) = point
   deallocate(catPolygon)

end subroutine ConcatenateToPolygon


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

end program matchFacetsToCells
