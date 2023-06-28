module in_mypoly_functions

    use omp_lib

    contains

    real function max_facet_side(n_vert,vertices,n_fcts,facets)
        implicit none

        integer, intent(in) :: n_vert, n_fcts
        integer, dimension(n_fcts*3), intent(in) :: facets
        real, dimension(n_vert*3), intent(in) :: vertices

        integer :: i, j, p, q, r
        real :: distance
        real, dimension(3) :: side

        max_facet_side = 0;

        do i = 1,n_fcts

            p = facets(3*(i-1)+1)
            q = facets(3*(i-1)+2)
            r = facets(3*(i-1)+3)
            
            distance = 0.0
            do j = 1,3
            distance = distance + ( vertices(3*(p-1)+j) - vertices(3*(q-1)+j) )**2
            end do
            side(1) = SQRT(distance)

            distance = 0.0
            do j = 1,3
            distance = distance + ( vertices(3*(q-1)+j) - vertices(3*(r-1)+j) )**2
            end do
            side(2) = SQRT(distance)

            distance = 0.0
            do j = 1,3
            distance = distance + ( vertices(3*(r-1)+j) - vertices(3*(p-1)+j) )**2
            end do
            side(3) = SQRT(distance)
    
            if (max_facet_side<MAXVAL(side)) then
                max_facet_side = MAXVAL(side)
            end if

        end do

    end function max_facet_side

    !! Subroutine to check if grid points of a given mesh defined by xgrid, ygrid and zgrid 
    !! are inside or outside a 3D STL defined by the vertices, facets and faceNormals.
    function is_grid_in_mypoly_func(n_vert,vertices,n_fcts,facets,incenters,faceNormals, &
        nx,xgrid,ny,ygrid,nz,zgrid,Ray_dir,L_char,max_height,tol)
        implicit none

        integer, intent(in) :: n_vert, n_fcts, nx, ny, nz
        real, intent(in) :: L_char, max_height, tol
        real, dimension(3), intent(in) :: Ray_dir
        real, dimension(n_vert*3), intent(in) :: vertices
        integer, dimension(n_fcts*3), intent(in) :: facets
        real, dimension(n_fcts*3), intent(in) :: incenters, faceNormals
        real, dimension(nx), intent(in) :: xgrid
        real, dimension(ny), intent(in) :: ygrid
        real, dimension(nz), intent(in) :: zgrid
        
        logical, dimension(nx,ny,nz) :: is_grid_in_mypoly_func

        ! logical, external :: is_in_mypoly

        integer :: ix, iy, iz
        real, dimension(3) :: Origin
        
        !$ call OMP_SET_NUM_THREADS(8)
        !$OMP parallel do default(shared) private(ix,iy,iz,Origin) schedule(dynamic)
        do iy = 1,ny
            do iz = 1,nz
                if (zgrid(iz)>max_height) then
                    is_grid_in_mypoly_func(:,iy,iz) = .false.
                else
                    do ix = 1,nx
                        Origin(1) = xgrid(ix)
                        Origin(2) = ygrid(iy)
                        Origin(3) = zgrid(iz)
                        is_grid_in_mypoly_func(ix,iy,iz) = is_in_mypoly(n_vert,vertices,n_fcts,facets, &
                        incenters,faceNormals,Origin,Ray_dir,L_char,tol);
                    end do
                end if
            end do
        end do
        !$OMP end parallel do
    end function is_grid_in_mypoly_func



    !! Subroutine to check if grid points of a given mesh defined by xgrid, ygrid and zgrid 
    !! are inside or outside a 3D STL defined by the vertices, facets and faceNormals.
    subroutine is_grid_in_mypoly_sub(solid,n_vert,vertices,n_fcts,facets,incenters,faceNormals, &
        nx,xgrid,ny,ygrid,nz,zgrid,Ray_dir,L_char,max_height,tol)
        implicit none

        integer, intent(in) :: n_vert, n_fcts, nx, ny, nz
        real, intent(in) :: L_char, max_height, tol
        real, dimension(3), intent(in) :: Ray_dir
        real, dimension(n_vert*3), intent(in) :: vertices
        integer, dimension(n_fcts*3), intent(in) :: facets
        real, dimension(n_fcts*3), intent(in) :: incenters, faceNormals
        real, dimension(nx), intent(in) :: xgrid
        real, dimension(ny), intent(in) :: ygrid
        real, dimension(nz), intent(in) :: zgrid
        logical, dimension(nx,ny,nz), intent(out) :: solid

        ! logical, external :: is_in_mypoly

        integer :: ix, iy, iz
        real, dimension(3) :: Origin
        
        do iy = 1,ny
            do iz = 1,nz
                if (zgrid(iz)>max_height) then
                    solid(:,iy,iz) = .false.
                else 
                    do ix = 1,nx
                        Origin(1) = xgrid(ix)
                        Origin(2) = ygrid(iy)
                        Origin(3) = zgrid(iz)
                        solid(ix,iy,iz) = is_in_mypoly(n_vert,vertices,n_fcts,facets,incenters,faceNormals, &
                        Origin,Ray_dir,L_char,tol);
                    end do
                end if
            end do
        end do
    end subroutine is_grid_in_mypoly_sub



    !! Function to check if a given point (Origin) is inside or outside a 3D STL defined by the 
    !! vertices, facets and faceNormals.
    logical function is_in_mypoly(n_vert,vertices,n_fcts,facets,incenters,faceNormals,Origin,Ray_dir,L_char,tol)
        implicit none

        integer, intent(in) :: n_vert, n_fcts
        real, intent(in) ::  L_char, tol
        real, dimension(n_vert*3), intent(in) :: vertices
        integer, dimension(n_fcts*3), intent(in) :: facets
        real, dimension(n_fcts*3), intent(in) :: incenters, faceNormals
        real, dimension(3), intent(in) :: Origin, Ray_dir

        ! logical, external :: point_triangle_intersect, point_on_triangle

        real, dimension(3) :: vert1, vert2, vert3
        real, dimension(n_fcts*9) :: facet_intersect_old
        integer :: i, i_facet, kk, j, counter, count2, max_i(1)
        logical :: withinBox, intersect_on(2)

        !! Take x axis, or y axis or positive z axis even with carefull cohice of this ray direction there 
        !! might be wrong prediction. To resolve this one must shoot at the least three rays along the three 
        !! cartresian direction and decide based on the majority outcome.
        !! Dir = (/0,0,1/)       !! along positive z-axis. should be decided carefully. If changed the first if clause
        !! inside 'iloop' below also needs to be modified accordingly or commented
        
        is_in_mypoly = .false.
        max_i = MAXLOC(ABS(Ray_dir))
        counter = 0

        kk = 1
        iloop: do i_facet = 1,n_fcts

            select case(max_i(1))
                case(1) 
                    if (incenters(kk+1) > Origin(2)-L_char .and. incenters(kk+1) < Origin(2)+L_char .and. &
                        incenters(kk+2) > Origin(3)-L_char .and. incenters(kk+2) < Origin(3)+L_char) then
                        withinBox = .true.
                    else
                        withinBox = .false.
                    end if
                case(2)
                    if (incenters(kk) > Origin(1)-L_char .and. incenters(kk) < Origin(1)+L_char .and. &
                        incenters(kk+2) > Origin(3)-L_char .and. incenters(kk+2) < Origin(3)+L_char) then
                        withinBox = .true.
                    else
                        withinBox = .false.
                    end if
                case(3)
                    if (incenters(kk) > Origin(1)-L_char .and. incenters(kk) < Origin(1)+L_char .and. &
                        incenters(kk+1) > Origin(2)-L_char .and. incenters(kk+1) < Origin(2)+L_char) then
                        withinBox = .true.
                    else
                        withinBox = .false.
                    end if
                case default
                    print *,'Error!! The pointLineSegmentIntersect function is valid only for 3D inputs.'
            end select

            !! not in the user choice. The below if clause should be decided carefully and must be modified it the 
            !! direction of the ray from origin point is changed
            if (withinBox) then
                
                do i = 1,3
                    vert1(i) = vertices(3*(facets(kk)-1) + i)
                    vert2(i) = vertices(3*(facets(kk+1)-1) + i)
                    vert3(i) = vertices(3*(facets(kk+2)-1) + i)
                end do
                
                intersect_on = point_triangle_intersect(Origin,Ray_dir,vert1,vert2,vert3, &
                                                        incenters(kk:kk+2),faceNormals(kk:kk+2),tol)

                if (intersect_on(1)) then

                    !! origin points lying on the surface of spherical objects  
                    ! if (point_on_triangle(Origin,vert1, vert2, vert3, incenters(kk:kk+2), faceNormals(kk:kk+2), tol)) then
                    if (intersect_on(2)) then                        
                    
                        is_in_mypoly = .true. !! inside
                        return

                    else
                            
                        if (counter == 0) then
                            counter=counter+1;
                            facet_intersect_old(9*(counter-1)+1:9*(counter-1)+3) = vert1
                            facet_intersect_old(9*(counter-1)+4:9*(counter-1)+6) = vert2
                            facet_intersect_old(9*(counter-1)+7:9*counter) = vert3
                        else

                            count2 = 0

                            do_j: do j = 1,3*counter
                                
                                if (vert1(1) == facet_intersect_old(3*(j-1)+1) .and. vert1(2) == facet_intersect_old(3*(j-1)+2) &
                                .and. vert1(3) == facet_intersect_old(3*j)) then
                                    count2 = 1
                                    exit do_j
                                elseif (vert2(1) == facet_intersect_old(3*(j-1)+1) &
                                    .and. vert2(2) == facet_intersect_old(3*(j-1)+2) &
                                    .and. vert2(3) == facet_intersect_old(3*j)) then
                                    count2 = 1
                                    exit do_j
                                elseif (vert3(1) == facet_intersect_old(3*(j-1)+1) &
                                    .and. vert3(2) == facet_intersect_old(3*(j-1)+2) &
                                    .and. vert3(3) == facet_intersect_old(3*j)) then
                                    count2 = 1
                                    exit do_j
                                end if

                            end do do_j

                            if (count2 == 0) then
                                counter=counter+1
                                facet_intersect_old(9*(counter-1)+1:9*(counter-1)+3) = vert1
                                facet_intersect_old(9*(counter-1)+4:9*(counter-1)+6) = vert2
                                facet_intersect_old(9*(counter-1)+7:9*counter) = vert3
                            end if

                        end if

                    end if
                end if
            end if

            kk = kk + 3

        end do iloop

        if(is_in_mypoly .neqv. .true.) then
            if (MOD(counter,2) == 0) then
                is_in_mypoly = .false. !! outside
            else
                is_in_mypoly = .true. !! inside
            end if
        end if

    end function is_in_mypoly



    !! Function to check if a given point (Origin) lie on a triangular plane element formed by the
    !! three points A, B and C (given by vertA, vertB and vertC) in 3D 
    !! logical function point_on_triangle(Origin, vertA, vertB, vertC, incenter, faceNormal, tol)
    logical function point_on_triangle(Origin, edge1, edge2, ray2, incenter, faceNormal, tol)
        implicit none

        !! real, dimension(3), intent(in) :: Origin, vertA, vertB, vertC, incenter, faceNormal
        real, dimension(3), intent(in) :: Origin, edge1, edge2, ray2, incenter, faceNormal
        real, intent(in) :: tol

        !! real, dimension(3) :: Vec, edge1, edge2, ray2, deno
        real, dimension(3) :: Vec, deno
        real :: max_deno, neu_1, neu_2, a, b
        integer :: i, max_i(1)

        Vec = Origin-incenter
        do i=1,3
            if (ABS(Vec(i)) < tol) then
                Vec(i) = 0
            end if
        end do

        if (ABS(DOT_PRODUCT(Vec,faceNormal))<tol) then

            !! edge1 = vertB - vertA
            !! edge2 = vertC - vertA
            !! ray2 = Origin - vertA

            deno(1) = edge1(1)*edge2(2) - edge1(2)*edge2(1)
            deno(2) = edge1(1)*edge2(3) - edge1(3)*edge2(1)
            deno(3) = edge1(2)*edge2(3) - edge1(3)*edge2(2)

            max_i = MAXLOC(ABS(deno))
            max_deno = deno(max_i(1))

            select case(max_i(1))
                case(1)
                    neu_1 = -edge2(1)*ray2(2)+edge2(2)*ray2(1)
                    neu_2 = -ray2(1)*edge1(2) + ray2(2)*edge1(1)
                case(2)
                    neu_1 = -edge2(1)*ray2(3)+edge2(3)*ray2(1)
                    neu_2 = -ray2(1)*edge1(3)+ray2(3)*edge1(1)
                case(3)
                    neu_1 = -edge2(2)*ray2(3)+edge2(3)*ray2(2)
                    neu_2 = -ray2(2)*edge1(3)+ray2(3)*edge1(2)
                case default
                    print *, "Error! 'point_on_triangle' function is valid only for 3D inputs."
            end select

            if (max_deno/=0) then
                a = neu_1/max_deno
                b = neu_2/max_deno
                if (a>=0 .and. b>=0 .and. a+b<=1) then
                    point_on_triangle = .true.
                else
                    point_on_triangle = .false.
                end if
            else
                point_on_triangle = .false.
            end if

        else

            point_on_triangle = .false.

        end if

    end function point_on_triangle



    !! Function to check if a ray (given by Dir) from a origin point (given by Origin) intersects the 
    !! triangular element formed by the three points A, B and C (given by vertA, vertB and vertC) in 3D 
    function point_triangle_intersect(Origin, Dir, vertA, vertB, vertC, incenter, faceNormal, tol)
        implicit none

        real, dimension(3), intent(in) :: Origin, Dir, vertA, vertB, vertC, incenter, faceNormal
        real, intent(in) :: tol

        logical, dimension(2) :: point_triangle_intersect

        ! real, external :: determinant
        ! logical, external :: point_line_segment_intersect

        real, dimension(3) :: edge1, edge2, ray2
        real :: det1, det_t, det_a, det_b, t, a, b
        integer :: i

        edge1 = vertB - vertA
        edge2 = vertC - vertA
        
        ray2 = Origin - vertA
        do i=1,3
            if (ABS(ray2(i)) < tol) then
                ray2(i) = 0
            end if
        end do

        det1 = determinant(-Dir,edge1,edge2)
        det_t = determinant(ray2,edge1,edge2)
        det_a = determinant(-Dir,ray2,edge2)
        det_b = determinant(-Dir,edge1,ray2)

        if (ABS(det1) > tol) then   ! the ray intersects the plane of the triangle

            t = det_t/det1
            a = det_a/det1
            b = det_b/det1

            ! ray intersects the plane within the domain bounded by the triangle (inclusive of edge and vertices)
            if (a>=0 .and. b>=0 .and. a+b<=1 .and. t>=0) then 
                point_triangle_intersect(1) = .true.
                if (point_on_triangle(Origin, edge1, edge2, ray2, incenter, faceNormal, tol)) then ! Origin point lies on the triangle
                    point_triangle_intersect(2) = .true.
                else        ! %Origin point dos not lie on the triangle
                    point_triangle_intersect(2) = .false.
                end if
            else            ! ray DOES NOT intersect the plane within the domain bounded by the triangle
                point_triangle_intersect(1) = .false.
                point_triangle_intersect(2) = .false.
            end if

        else if (ABS(det_a)<tol .and. ABS(det_b)<tol .and. ABS(det_t)<tol) then ! ray is parallel to the plane of triangle and lie on the plane

            if (point_on_triangle(Origin, edge1, edge2, ray2, incenter, faceNormal, tol)) then ! Origin point lies on the triangle
                point_triangle_intersect(1) = .true.
                point_triangle_intersect(2) = .true.
            else        ! % Origin point is outside the triangle
                point_triangle_intersect(1) = .false.
                point_triangle_intersect(2) = .false.
            end if
            
            ! if (point_line_segment_intersect(Origin,Dir,vertA,vertB)) then
            !     point_triangle_intersect = .true.
            ! elseif (point_line_segment_intersect(Origin,Dir,vertB,vertC)) then
            !     point_triangle_intersect = .true.
            ! elseif (point_line_segment_intersect(Origin,Dir,vertC,vertA)) then
            !     point_triangle_intersect = .true.
            ! else                    ! ray DOES NOT intersect any of the triangle edges
            !     point_triangle_intersect = .false.
            ! end if

        else        ! The ray is parallel to the plane of triangle with a offset and DOES NOT intersect it

            point_triangle_intersect(1) = .false.
            point_triangle_intersect(2) = .false.

        end if

    end function point_triangle_intersect



    !! Function to check if a ray (given by Dir) from a origin point (given by Origin) intersects the 
    !! line segment between the two points A and B(given by vert1 and vert2) in 3D 
    logical function point_line_segment_intersect(Origin, Dir, vert1, vert2)
        implicit none

        real, dimension(3), intent(in) :: Origin, Dir, vert1, vert2

        ! real, external :: norm_cross

        integer :: max_i(1)
        real :: max_deno, neu_1, neu_2, t_1, t_2
        real, dimension(3) :: line, ray2, deno

        line = vert2 - vert1        !! get the vector from points A to B which forms the line segment
        ray2 = Origin - vert1       !! get the vector from points A to Origin

        deno(1) = -Dir(1)*line(2) + Dir(2)*line(1)
        deno(2) = -Dir(1)*line(3) + Dir(3)*line(1)
        deno(3) = -Dir(2)*line(3) + Dir(3)*line(2)

        max_i = MAXLOC(ABS(deno))
        max_deno = deno(max_i(1))

        select case(max_i(1))
            case(1) 
                neu_1 = -line(1)*ray2(2) + line(2)*ray2(1)
                neu_2 = ray2(1)*Dir(2) - ray2(2)*Dir(1)
            case(2)
                neu_1 = -line(1)*ray2(3) + line(3)*ray2(1)
                neu_2 = ray2(1)*Dir(3) - ray2(3)*Dir(1)
            case(3)
                neu_1 = -line(2)*ray2(3) + line(3)*ray2(2)
                neu_2 = ray2(2)*Dir(3) - ray2(3)*Dir(2)
            case default
                print *,'Error!! The pointLineSegmentIntersect function is valid only for 3D inputs.'
        end select

        if (max_deno/=0) then

            t_1 = neu_1/max_deno
            t_2 = neu_2/max_deno
            if(t_1>=0 .and. t_2>=0 .and. t_2<=1) then
                point_line_segment_intersect = .true.
            else
                point_line_segment_intersect = .false.
            end if

        else if (max_deno==0 .and. neu_1==0 .and. neu_2==0) then

            if (line(1)==0 .and. line(2)==0 .and. line(3)==0) then
                print *, 'The verices forming the line segment must have different coordinates'
            else if (Dir(1)==0 .and. Dir(2)==0 .and. Dir(3)==0) then
                print *, 'Magnitude of the direction vector cannot be zero'
            else if(norm_cross(ray2,line) == 0) then

                max_i = MAXLOC(line)
                if(DOT_PRODUCT(Dir,line)>0 .and. ray2(max_i(1))/line(max_i(1)) <= 1 ) then
                    point_line_segment_intersect = .true.
                elseif (DOT_PRODUCT(Dir,line)<0 .and. ray2(max_i(1))/line(max_i(1)) >= 0 ) then
                    point_line_segment_intersect = .true.
                else
                    point_line_segment_intersect = .false.
                end if

            else
                point_line_segment_intersect = .false.
            end if

        else

            point_line_segment_intersect = .false.

        end if

    end function point_line_segment_intersect



    !! Function to calculate the norm of the cross product of two vectors (must be 3x1 vectors)
    real function norm_cross(A,B)
        implicit none
        
        real, dimension(3), intent(in) :: A, B

        norm_cross = SQRT((A(2)*B(3) - A(3)*B(2))**2 + (-A(1)*B(3) + A(3)*B(1))**2 + (A(1)*B(2) - A(2)*B(1))**2)

    end function norm_cross



    !! Function to calculate determinant of a 3x3 matrix (M) with column vectors A, B and C (each 3x1) such that M = [A B C]
    real function determinant(A, B, C)
        implicit none

        real, dimension(3), intent(in) :: A, B, C

        determinant = (A(2)*B(3) - A(3)*B(2))*C(1) + (A(3)*B(1) - A(1)*B(3))*C(2) + (A(1)*B(2) - A(2)*B(1))*C(3)

    end function determinant

end module in_mypoly_functions
