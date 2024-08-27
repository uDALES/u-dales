function array_shifted = generateStaggeredCubesVariableNew(Hx, Hy, Hz, Cx, Cy, Nx, Ny, Lx, Ly, divisions, tol)
    scale = [Hx Hy Hz];
    shift = [];

    for i=1:Nx
        if ~mod(i,2)
            shift = [shift; [i*(Cx+Hx)-Hx/2-Cx/2, 0, Hz/2]];
        end

        for j=1:Ny
            if mod(i,2)
                shift = [shift; [i*(Cx+Hx)-Hx/2-Cx/2, j*(Cy+Hy)-Hy/2-Hy/4, Hz/2]];
            else
                shift = [shift; [i*(Cx+Hx)-Hx/2-Cx/2, j*(Cy+Hy), Hz/2]];
            end

        end
    end

    array = operateUnitCube(scale, shift, divisions);

%     figure; 
%     trisurf(array);
%     axis equal

    faceNormal = array.faceNormal;

    points = array.Points;
    % shift a tiny amount in each direction apart from -z to make sure
    % solid points on boundary get labelled as such by inpolyhedron
    for n=1:size(array,1)
        if all(abs(faceNormal(n,:) - [1 0 0]) <= tol) % x
            points(array.ConnectivityList(n,:),:) = ...
                array.Points(array.ConnectivityList(n,:),:) + tol*[1 0 0];

        elseif all(abs(faceNormal(n,:) - [-1 0 0]) <= tol) % x
            points(array.ConnectivityList(n,:),:) = ...
                array.Points(array.ConnectivityList(n,:),:) + tol*[-1 0 0];

        elseif all(abs(faceNormal(n,:) - [0 1 0]) <= tol) % x
            points(array.ConnectivityList(n,:),:) ...
                = array.Points(array.ConnectivityList(n,:),:) + tol*[0 1 0];

        elseif all(abs(faceNormal(n,:) - [0 -1 0]) <= tol) % x
            points(array.ConnectivityList(n,:),:) ...
                = array.Points(array.ConnectivityList(n,:),:) + tol*[0 -1 0];

        elseif all(abs(faceNormal(n,:) - [0 0 1]) <= tol) % x
            points(array.ConnectivityList(n,:),:) ...
                = array.Points(array.ConnectivityList(n,:),:) + tol*[0 0 1];

        elseif all(abs(faceNormal(n,:) - [0 0 -1]) <= tol) % -z

        else
            disp("Normal not in x, y, or z direction")
        end
    end

    % Move points outside domain to the edge - this creates 1D facets!
%     points_og = points;
%     points_remove = [];
%     for n=1:size(points_og,1)
%         if points(n,2) < 0
%             point = points(n,:);
%             point(2) = 0;
%             if ~ismember(point, points_og, 'rows')
%                 points(n,:) = point;
%             else
%                 points_remove = [points_remove; n];
%             end
%         elseif points(n,2) > Ly
%             point = points(n,:);
%             point(2) = Ly;
%             if ~ismember(point, points_og, 'rows')
%                 points(n,:) = point;
%             else
%                 points_remove = [points_remove; n];
%             end
%         end
% 
%     end

%     figure; 
%     trisurf(triangulation(array.ConnectivityList, points));
%     axis equal

    % remove duplicate points
    %     points_remove = [];
%     for n=1:size(points,1)
%         if points(n,2) < 0
%             points_remove = [points_remove; n];
%         end
%         if points(n,2) > Ly
%             points_remove = [points_remove; n];
%         end
% 
%     end

%     % see if any facets are 1D - if so, remove
%     faces_remove = [];
%     for n=1:size(array)
%         verts = points(array.ConnectivityList(n,:), :);
%         area = 1/2*norm(cross(verts(2,:)-verts(1,:),verts(3,:)-verts(1,:)))
%         if (area < 1e-5)
%             faces_remove = [faces_remove; n];
%         end
%     end

    conn_list = array.ConnectivityList;

    if (divisions<2)
        for n=1:size(points,1)
            if points(n,2) < 0
                points(n,2) = 0;
            elseif points(n,2) > Ly
                points(n,2) = Ly;
            end
        end
    else
        points_remove = [];
        for n=1:size(points,1)
            if points(n,2) <= 0
                points_remove = [points_remove; n];
            end
            if points(n,2) >= Ly
                points_remove = [points_remove; n];
            end

        end

        faces_remove = [];
        for n=1:size(array)
            if all(ismember(array.ConnectivityList(n,:), points_remove))
                faces_remove = [faces_remove; n];
            end
        end
        conn_list(faces_remove,:) = [];
    end

    %points(points_remove,:) = [];
    %[conn_list,points] = patchCleanUnused(conn_list, points);

%     figure; 
%     trisurf(triangulation(conn_list, points));
%     axis equal

    array_shifted = triangulation(conn_list, points);


%      figure; 
%      trisurf(array_shifted);
%      
%      axis equal
% 
%     xlim([0 Lx])
    %ylim([0 Ly])
    %zlim([0 Lz])
