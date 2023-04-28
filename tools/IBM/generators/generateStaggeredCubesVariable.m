function array_shifted = generateStaggeredCubesVariable(Hx, Hy, Hz, Cx, Cy, Nx, Ny, Lx, Ly, divisions, tol)
    scale = [Hx Hy Hz];
    shift = [];

    for i=1:Nx
        if ~mod(i,2)
            shift = [shift; [(2*i-1)*Cx, 0, Hz/2]];
        end

        for j=1:Ny
            if mod(i,2)
                shift = [shift; [(2*i-1)*Cx, (2*j-1)*Cy, Hz/2]];
            else
                shift = [shift; [(2*i-1)*Cx, (2*j)*Cy, Hz/2]];
            end

        end
    end

    array = operateUnitCube(scale, shift, divisions);
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

    % Move points outside domain to the edge
    for n=1:size(points,1)
        if points(n,2) <= 0
            points(n,2) = 0;
        end
        if points(n,2) >= Ly
            points(n,2) = Ly;
        end

    end

    array_shifted = triangulation(array.ConnectivityList, points);


    figure; 
    trisurf(array_shifted);
    
    axis equal
    xlim([0 Lx])
    ylim([0 Ly])
    %zlim([0 Lz])
