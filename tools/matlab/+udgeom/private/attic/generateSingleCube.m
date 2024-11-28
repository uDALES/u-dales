function array_shifted = generateSingleCube(Hx, Hy, Hz, Lx, Ly, divisions, shift, tol)
    scale = [Hx Hy Hz];
%     shift = [];
%     for i=1:Nx
%         for j=1:Ny
%             shift = [shift; [(2*i-1)*Hx, (2*j-1)*Hy, Hz/2]];
%         end
%     end
    %shift = [Lx/2 Ly/2 Hz/2];

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

    array_shifted = triangulation(array.ConnectivityList, points);

    figure; 
    trisurf(array);
    
    axis equal
    xmax = Lx;
    ymax = Ly;
    zmax = 1.5*Hz;
    xlim([0 xmax])
    ylim([0 ymax])
    zlim([0 zmax])
end