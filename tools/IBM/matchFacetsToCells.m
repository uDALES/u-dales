function [facet_sections] = matchFacetsToCells(TR, fluid_IB, solid_IB, ...
    fluid_IB_xyz, solid_IB_xyz, xgrid, ygrid, zgrid, include_diagonals)

% Assumes fluid_IB etc have been created already

% fluid_IB_xyz_cut = zeros(0,3);
% fluid_IB_facids = cell(0,0);
% fluid_IB_facareas = cell(0,0);
% fluid_IB_dists = cell(0,0);

dx = xgrid(2) - xgrid(1);
dy = ygrid(2) - ygrid(1);
dz = zgrid(2) - zgrid(1);

itot = length(xgrid);
jtot = length(ygrid);
ktot = length(zgrid);

facet_sections = [];

figure
for facet=1:size(TR.ConnectivityList,1)
%     disp(num2str(facet))
%     abs(TR.faceNormal(facet))
%     xgrid(1)
%     abs(abs(TR.faceNormal(facet)) - [1 0 0])
    if (xgrid(1) == 0 && all(abs(abs(TR.faceNormal(facet)) - [1 0 0]) <= 1e-8)) ...
    || (ygrid(1) == 0 && all(abs(abs(TR.faceNormal(facet)) - [0 1 0]) <= 1e-8)) ...
    || (zgrid(1) == 0 && all(abs(abs(TR.faceNormal(facet)) - [0 0 1]) <= 1e-8))
        continue
    end

    % Only search grid cells that can possibly contain it
    xmin = min(TR.Points(TR.ConnectivityList(facet,:),1));
    ymin = min(TR.Points(TR.ConnectivityList(facet,:),2));
    zmin = min(TR.Points(TR.ConnectivityList(facet,:),3));
    xmax = max(TR.Points(TR.ConnectivityList(facet,:),1));
    ymax = max(TR.Points(TR.ConnectivityList(facet,:),2));
    zmax = max(TR.Points(TR.ConnectivityList(facet,:),3));

    if ((abs(zmin) < eps) && (abs(zmax) < eps) && ... % at ground level
            all(abs(TR.faceNormal(facet) - [0 0 -1]) < eps)) % pointing down
        disp(['Facet ' num2str(facet) ' is on the ground and pointing down - ignore'])
        continue
    end

    il = find(xmin >= [xgrid-dx/2, xgrid(end)+dx/2], 1, 'last');
    iu = find(xmax <= [xgrid(1)-dx/2, xgrid + dx/2], 1, 'first');
    jl = find(ymin >= [ygrid-dy/2, ygrid(end)+dy/2], 1, 'last');
    ju = find(ymax <= [ygrid(1)-dy/2, ygrid + dy/2], 1, 'first');
    kl = find(zmin >= [zgrid-dz/2, zgrid(end)+dz/2], 1, 'last'); % Think about non-equidistance
    ku = find(zmax <= [zgrid(1)-dz/2, zgrid + dz/2], 1, 'first');

    %tol = 1e-6;
%     il = find(xmin - [xgrid-dx/2, xgrid(end)+dx/2] > tol, 1, 'last');
%     iu = find(xmax - [xgrid(1)-dx/2, xgrid + dx/2] < tol, 1, 'first')
%     jl = find(ymin - [ygrid-dy/2, ygrid(end)+dy/2] > tol, 1, 'last')
%     ju = find(ymax - [ygrid(1)-dy/2, ygrid + dy/2] < tol, 1, 'first')
%     kl = find(zmin - [zgrid-dz/2, zgrid(end)+dz/2] > tol, 1, 'last') % Think about non-equidistance
%     ku = find(zmax - [zgrid(1)-dz/2, zgrid + dz/2] < tol, 1, 'first')

    % facets on domain edge
    if il < 1; il = 1; end
    if jl < 1; jl = 1; end
    if kl < 1; kl = 1; end
    if iu > length(xgrid); iu = length(xgrid); end
    if ju > length(ygrid); ju = length(ygrid); end
    

    stopflag = false;    
    for i=il:iu
        for j=jl:ju
            for k=kl:ku
                if ~(fluid_IB(i,j,k) || solid_IB(i,j,k))
                    continue
                end


                % Define corners of cube
                xl = xgrid(i) - dx/2;
                xu = xgrid(i) + dx/2;
                yl = ygrid(j) - dy/2;
                yu = ygrid(j) + dy/2;
                zl = zgrid(k) - dz/2;
                zu = zgrid(k) + dz/2;

                % Define planes corresponding to cell faces
                xrange = [xl, xu];
                yrange = [yl, yu];
                zrange = [zl, zu];

                planes = [1  0  0  xrange(2);
                    -1  0  0 -xrange(1);
                    0  1  0  yrange(2);
                    0 -1  0 -yrange(1);
                    0  0  1  zrange(2);
                    0  0 -1 -zrange(1)];

                verts = TR.Points(TR.ConnectivityList(facet,:),:); %nverts x npoints (3x3)
                clip = sutherlandHodgman3D(verts, planes);

                tol = 1e-10;
                if ~isempty(clip)
                    if (size(clip,1) == 3) % triangle
                        area = 1/2*norm(cross(clip(2,:)-clip(1,:),clip(3,:)-clip(1,:)));
                        
                        if (area < tol)
                            disp(['Facet ' num2str(facet) ' in cell ' num2str(i) ',' num2str(j) ',' num2str(k) ' has zero area'])
                            continue
                        end

                        tri = triangulation([1 2 3], clip);

                    elseif(size(clip,1) > 3) % other polygon
                        xangle = dot(TR.faceNormal(facet), [1 0 0]);
                        yangle = dot(TR.faceNormal(facet), [0 1 0]);
                        zangle = dot(TR.faceNormal(facet), [0 0 1]);
                        [angle, dir] = max(abs([xangle, yangle, zangle]));

                        if dir == 1 % Drop x coordinate, project to yz plane
                            planeNormal = [1 0 0];
                            ids = [2 3];
                        elseif dir == 2
                            planeNormal = [0 1 0]; % Drop y coordinate, project to xz plane
                            ids = [1 3];
                        elseif dir == 3
                            planeNormal = [0 0 1]; % Drop z coordinate, project to xy plane
                            ids = [1 2];
                        end

                        projVert = zeros(size(clip,1),3);
                        for m = 1:size(clip,1)
                            projVert(m,:) = clip(m,:) - dot(clip(m,:), planeNormal) * planeNormal;
                        end

                        projArea = polyarea(projVert(:,ids(1)), projVert(:,ids(2)));
                        area = projArea / angle;

                        if (area < tol)
                            disp(['Facet ' num2str(facet) ' in cell ' num2str(i) ',' num2str(j) ',' num2str(k) ' has zero area'])
                            continue
                        end

                        switch size(clip,1)
                            case 4
                                tri = triangulation([1 2 3;1 3 4], clip);
                            case 5
                                tri = triangulation([1 2 3;1 3 4;1 4 5], clip);
                            case 6
                                tri = triangulation([1 2 3;1 3 4;1 4 5;1 5 6], clip);
                        end

                    else
                        error('something wrong with clipped polygon')
                    end
                    
                    % Formalise this
                    if ((xgrid(i) == 0. || ygrid(j) == 0.))% || zgrid(k) == 0.))% && k==1)
                        area = area * 2;
                        % Account for periodicity - flux at point N+1
                        % Or due to facet at ground level
                    end

                    facet_section = zeros(1,9);
                    facet_section(1) = facet;
                    facet_section(2) = area;

                    dists = nan(27,1);
                    angles = nan(27,1);
                    typeids = nan(27,1);
                    BIs = nan(27,3);
                    search_adj = false;

                    if fluid_IB(i,j,k)
                        [~, loc] = ismember([xgrid(i), ygrid(j), zgrid(k)], fluid_IB_xyz, 'rows');
                        facet_section(3) = loc; 

                        xyz27 = [xgrid(i), ygrid(j), zgrid(k)];
                        [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz27, 'UseSubSurface', false);
                        angle = dot(TR.faceNormal(facet), (xyz27 - BI)/vecnorm((xyz27 - BI)));
                        %if (typeid == 3) 
                        % Wall-normal defined, use this cell
                        if (abs(angle - 1) < eps)
                            id = 27;
                            xyz = xyz27;
                        else
                            % Not normal, search adjacent fluid IB cells
                            search_adj = true;
                            if (dist > 0)
                                % Include in comparison
                                dists(27) = dist;
                                typeids(27) = typeid;
                                angles(27) = dot(TR.faceNormal(facet), (xyz27 - BI)/vecnorm((xyz27 - BI)));
                                BIs(27,:) = BI;
                            end
                        end
                    end

                    if (solid_IB(i,j,k) || search_adj)
                        if (solid_IB(i,j,k))
                            [~, loc] = ismember([xgrid(i), ygrid(j), zgrid(k)], solid_IB_xyz, 'rows');
                            facet_section(4) = loc;
                        end

                        if i~=1
                            if fluid_IB(i-1,j,k)
                                xyz1 = [xgrid(i-1), ygrid(j), zgrid(k)];
                                [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz1, 'UseSubSurface', false);
                                dists(1) = dist;
                                typeids(1) = typeid;
                                angles(1) = dot(TR.faceNormal(facet), (xyz1- BI)/vecnorm((xyz1 - BI)));
                                BIs(1,:) = BI;
                            end
                        end
                        if i~=itot
                            if fluid_IB(i+1,j,k)
                                xyz2 = [xgrid(i+1), ygrid(j), zgrid(k)];
                                [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz2, 'UseSubSurface', false);
                                dists(2) = dist;
                                typeids(2) = typeid;
                                angles(2) = dot(TR.faceNormal(facet), (xyz2- BI)/vecnorm((xyz2 - BI)));
                                BIs(2,:) = BI;
                            end
                        end
                        if j~=1
                            if fluid_IB(i,j-1,k)
                                xyz3 = [xgrid(i), ygrid(j-1), zgrid(k)];
                                [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz3, 'UseSubSurface', false);
                                dists(3) = dist;
                                typeids(3) = typeid;
                                angles(3) = dot(TR.faceNormal(facet), (xyz3- BI)/vecnorm((xyz3 - BI)));
                                BIs(3,:) = BI;
                            end
                        end
                        if j~=jtot
                            if fluid_IB(i,j+1,k)
                                xyz4 = [xgrid(i), ygrid(j+1), zgrid(k)];
                                [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz4, 'UseSubSurface', false);
                                dists(4) = dist;
                                typeids(4) = typeid;
                                angles(4) = dot(TR.faceNormal(facet), (xyz4- BI)/vecnorm((xyz4 - BI)));
                                BIs(4,:) = BI;
                            end
                        end
                        if k~=1
                            if fluid_IB(i,j,k-1)
                                xyz5 = [xgrid(i), ygrid(j), zgrid(k-1)];
                                [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz5, 'UseSubSurface', false);
                                dists(5) = dist;
                                typeids(5) = typeid;
                                angles(5) = dot(TR.faceNormal(facet), (xyz5 - BI)/vecnorm((xyz5 - BI)));
                                BIs(5,:) = BI;
                            end
                        end
                        if k~=ktot
                            if fluid_IB(i,j,k+1)
                                xyz6 = [xgrid(i), ygrid(j), zgrid(k+1)];
                                [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz6, 'UseSubSurface', false);
                                dists(6) = dist;
                                typeids(6) = typeid;
                                angles(6) = dot(TR.faceNormal(facet), (xyz6- BI)/vecnorm((xyz6 - BI)));
                                BIs(6,:) = BI;
                            end
                        end

                        if (include_diagonals)
                            if (i~=1 && j~=1)
                                if fluid_IB(i-1,j-1,k)
                                    xyz7 = [xgrid(i-1), ygrid(j-1), zgrid(k)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz7, 'UseSubSurface', false);
                                    dists(7) = dist;
                                    typeids(7) = typeid;
                                    angles(7) = dot(TR.faceNormal(facet), (xyz7 - BI)/vecnorm((xyz7 - BI)));
                                end
                            end
                            if (i~=1 && j~=jtot)
                                if fluid_IB(i-1,j+1,k)
                                    xyz8 = [xgrid(i-1), ygrid(j+1), zgrid(k)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz8, 'UseSubSurface', false);
                                    dists(8) = dist;
                                    typeids(8) = typeid;
                                    angles(8) = dot(TR.faceNormal(facet), (xyz8 - BI)/vecnorm((xyz8 - BI)));
                                end
                            end
                            if (i~=itot && j~=1)
                                if fluid_IB(i+1,j-1,k)
                                    xyz9 = [xgrid(i+1), ygrid(j-1), zgrid(k)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz9, 'UseSubSurface', false);
                                    dists(9) = dist;
                                    typeids(9) = typeid;
                                    angles(9) = dot(TR.faceNormal(facet), (xyz9 - BI)/vecnorm((xyz9 - BI)));
                                end
                            end
                            if (i~=itot && j~=jtot)
                                if fluid_IB(i+1,j+1,k)
                                    xyz10 = [xgrid(i+1), ygrid(j+1), zgrid(k)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz10, 'UseSubSurface', false);
                                    dists(10) = dist;
                                    typeids(10) = typeid;
                                    angles(10) = dot(TR.faceNormal(facet), (xyz10 - BI)/vecnorm((xyz10 - BI)));
                                end
                            end

                            if (i~=1 && k~=1)
                                if fluid_IB(i-1,j,k-1)
                                    xyz11 = [xgrid(i-1), ygrid(j), zgrid(k-1)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz11, 'UseSubSurface', false);
                                    dists(11) = dist;
                                    typeids(11) = typeid;
                                    angles(11) = dot(TR.faceNormal(facet), (xyz11- BI)/vecnorm((xyz11 - BI)));
                                end
                            end
                            if (i~=1 && k~=ktot)
                                if fluid_IB(i-1,j,k+1)
                                    xyz12 = [xgrid(i-1), ygrid(j), zgrid(k+1)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz12, 'UseSubSurface', false);
                                    dists(12) = dist;
                                    typeids(12) = typeid;
                                    angles(12) = dot(TR.faceNormal(facet), (xyz12 - BI)/vecnorm((xyz12 - BI)));
                                end
                            end
                            if (i~=itot && k~=1)
                                if fluid_IB(i+1,j,k-1)
                                    xyz13 = [xgrid(i+1), ygrid(j), zgrid(k-1)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz13, 'UseSubSurface', false);
                                    dists(13) = dist;
                                    typeids(13) = typeid;
                                    angles(13) = dot(TR.faceNormal(facet), (xyz13 - BI)/vecnorm((xyz13 - BI)));
                                end
                            end
                            if (i~=itot && k~=ktot)
                                if fluid_IB(i+1,j,k+1)
                                    xyz14 = [xgrid(i+1), ygrid(j), zgrid(k+1)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz14, 'UseSubSurface', false);
                                    dists(14) = dist;
                                    typeids(14) = typeid;
                                    angles(14) = dot(TR.faceNormal(facet), (xyz14 - BI)/vecnorm((xyz14 - BI)));
                                end
                            end

                            if (j~=1 && k~=1)
                                if fluid_IB(i,j-1,k-1)
                                    xyz15 = [xgrid(i), ygrid(j-1), zgrid(k-1)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz15, 'UseSubSurface', false);
                                    dists(15) = dist;
                                    typeids(15) = typeid;
                                    angles(15) = dot(TR.faceNormal(facet), (xyz15 - BI)/vecnorm((xyz15 - BI)));
                                end
                            end
                            if (j~=1 && k~=ktot)
                                if fluid_IB(i,j-1,k+1)
                                    xyz16 = [xgrid(i), ygrid(j-1), zgrid(k+1)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz16, 'UseSubSurface', false);
                                    dists(16) = dist;
                                    typeids(16) = typeid;
                                    angles(16) = dot(TR.faceNormal(facet), (xyz16 - BI)/vecnorm((xyz16 - BI)));
                                end
                            end
                            if (j~=jtot && k~=1)
                                if fluid_IB(i,j+1,k-1)
                                    xyz17 = [xgrid(i), ygrid(j+1), zgrid(k-1)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz17, 'UseSubSurface', false);
                                    dists(17) = dist;
                                    typeids(17) = typeid;
                                    angles(17) = dot(TR.faceNormal(facet), (xyz17 - BI)/vecnorm((xyz17 - BI)));
                                end
                            end
                            if (j~=jtot && k~=ktot)
                                if fluid_IB(i,j+1,k+1)
                                    xyz18 = [xgrid(i), ygrid(j+1), zgrid(k+1)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz18, 'UseSubSurface', false);
                                    dists(18) = dist;
                                    typeids(18) = typeid;
                                    angles(18) = dot(TR.faceNormal(facet), (xyz18 - BI)/vecnorm((xyz18 - BI)));
                                end
                            end

                            if (i~=1 && j~=1 && k~=1)
                                if fluid_IB(i-1,j-1,k-1)
                                    xyz19 = [xgrid(i-1), ygrid(j-1), zgrid(k-1)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz19, 'UseSubSurface', false);
                                    dists(19) = dist;
                                    typeids(19) = typeid;
                                    angles(19) = dot(TR.faceNormal(facet), (xyz19 - BI)/vecnorm((xyz19 - BI)));
                                end
                            end
                            if (i~=itot && j~=1 && k~=1)
                                if fluid_IB(i+1,j-1,k-1)
                                    xyz20 = [xgrid(i+1), ygrid(j-1), zgrid(k-1)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz20, 'UseSubSurface', false);
                                    dists(20) = dist;
                                    typeids(20) = typeid;
                                    angles(20) = dot(TR.faceNormal(facet), (xyz20 - BI)/vecnorm((xyz20 - BI)));
                                end
                            end
                            if (i~=1 && j~=jtot && k~=1)
                                if fluid_IB(i-1,j+1,k-1)
                                    xyz21 = [xgrid(i-1), ygrid(j+1), zgrid(k-1)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz21, 'UseSubSurface', false);
                                    dists(21) = dist;
                                    typeids(21) = typeid;
                                    angles(21) = dot(TR.faceNormal(facet), (xyz21 - BI)/vecnorm((xyz21 - BI)));
                                end
                            end
                            if (i~=itot && j~=jtot && k~=1)
                                if fluid_IB(i+1,j+1,k-1)
                                    xyz22 = [xgrid(i+1), ygrid(j+1), zgrid(k-1)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz22, 'UseSubSurface', false);
                                    dists(22) = dist;
                                    typeids(22) = typeid;
                                    angles(22) = dot(TR.faceNormal(facet), (xyz22 - BI)/vecnorm((xyz22 - BI)));
                                end
                            end

                            if (i~=1 && j~=1 && k~=ktot)
                                if fluid_IB(i-1,j-1,k+1)
                                    xyz23 = [xgrid(i-1), ygrid(j-1), zgrid(k+1)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz23, 'UseSubSurface', false);
                                    dists(23) = dist;
                                    typeids(23) = typeid;
                                    angles(23) = dot(TR.faceNormal(facet), (xyz23 - BI)/vecnorm((xyz23 - BI)));
                                end
                            end
                            if (i~=itot && j~=1 && k~=ktot)
                                if fluid_IB(i+1,j-1,k+1)
                                    xyz24 = [xgrid(i+1), ygrid(j-1), zgrid(k+1)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz24, 'UseSubSurface', false);
                                    dists(24) = dist;
                                    typeids(24) = typeid;
                                    angles(24) = dot(TR.faceNormal(facet), (xyz24 - BI)/vecnorm((xyz24 - BI)));
                                end
                            end
                            if (i~=1 && j~=jtot && k~=ktot)
                                if fluid_IB(i-1,j+1,k+1)
                                    xyz25 = [xgrid(i-1), ygrid(j+1), zgrid(k+1)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz25, 'UseSubSurface', false);
                                    dists(25) = dist;
                                    typeids(25) = typeid;
                                    angles(25) = dot(TR.faceNormal(facet), (xyz25 - BI)/vecnorm((xyz25 - BI)));
                                end
                            end
                            if (i~=itot && j~=jtot && k~=ktot)
                                if fluid_IB(i+1,j+1,k+1)
                                    xyz26 = [xgrid(i+1), ygrid(j+1), zgrid(k+1)];
                                    [dist, BI, ~, typeid] = point2trimesh('Faces', tri.ConnectivityList, 'Vertices', tri.Points, 'QueryPoints', xyz26, 'UseSubSurface', false);
                                    dists(26) = dist;
                                    typeids(26) = typeid;
                                    angles(26) = dot(TR.faceNormal(facet), (xyz26 - BI)/vecnorm((xyz26 - BI)));
                                end
                            end
                        end

                        % Use maximum dot product
                        [~, id] = max(abs(angles) ./ (dists / (dx*dy*dz)^(1/3))); % minimise both angles and dists
                        dist = dists(id);
                        BI = BIs(id,:);
                        if isnan(dist)
                            % point2trimesh found BI but not dist?
                            % Ignore it and take maximum again
                            angles(id) = nan;
                            [~, id] = max(abs(angles));
                            dist = dists(id);
                            BI = BIs(id,:);
                        end

                        if (isnan(dist))
                            disp(['Facet ' num2str(facet) ' in cell ' num2str(i) ',' num2str(j) ',' num2str(k) ' has no neighbouring fluid cells'])
                            % probably internal, or if w grid, facets at
                            % the bottom of the domain, but check.
                            % Need to remove these sections from
                            % the total area calculation.
                            %clf
                            
                            if (xgrid(1) == 0)
                                title('u')
                            elseif (ygrid(1) == 0)
                                title('v')
                            elseif (zgrid(1) == 0)
                                title('w')
                            else
                                title('c')
                            end

                            view(3)
                            patch('Faces', TR.ConnectivityList, 'Vertices', TR.Points, 'FaceColor', ones(3,1)*69/100, 'FaceAlpha', 0.5)
                            %patch('Faces', [1 2 3], 'Vertices', TR.Points(TR.ConnectivityList(facet,:),:), 'FaceColor', [1,0,0], 'FaceAlpha', 0.5)
                            hold on
                            incenter = TR.incenter(1); faceNormal = TR.faceNormal;
                            quiver3(incenter(1), incenter(2), incenter(3), faceNormal(1), faceNormal(2), faceNormal(3))
                            V = [xl yl zl; xu yl zl; xu yu zl; xl yu zl; ...
                                xl yl zu; xu yl zu; xu yu zu; xl yu zu];
                            F = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
                            patch('Faces', F, 'Vertices', V, 'FaceColor', [1,1,1], 'FaceAlpha', 0.5)
                            trisurf(tri)
                            axis equal
                            xlabel('x')
                            ylabel('y')
                            zlabel('z')
                            xlim([0 xgrid(end)])
                            ylim([0 ygrid(end)])
                            zlim([0 zgrid(end)])
                            %drawnow
                            %pause(5)
                            continue
                        end

                    end

                    switch id
                        case 1
                            xyz = xyz1;
                        case 2
                            xyz = xyz2;
                        case 3
                            xyz = xyz3;
                        case 4
                            xyz = xyz4;
                        case 5
                            xyz = xyz5;
                        case 6
                            xyz = xyz6;
                        case 7
                            xyz = xyz7;
                        case 8
                            xyz = xyz8;
                        case 9
                            xyz = xyz9;
                        case 10
                            xyz = xyz10;
                        case 11
                            xyz = xyz11;
                        case 12
                            xyz = xyz12;
                        case 13
                            xyz = xyz13;
                        case 14
                            xyz = xyz14;
                        case 15
                            xyz = xyz15;
                        case 16
                            xyz = xyz16;
                        case 17
                            xyz = xyz17;
                        case 18
                            xyz = xyz18;
                        case 19
                            xyz = xyz19;
                        case 20
                            xyz = xyz20;
                        case 21
                            xyz = xyz21;
                        case 22
                            xyz = xyz22;
                        case 23
                            xyz = xyz23;
                        case 24
                            xyz = xyz24;
                        case 25
                            xyz = xyz25;
                        case 26
                            xyz = xyz26;
                        case 27
                            xyz = xyz27;
                    end

                    [~, loc] = ismember(xyz, fluid_IB_xyz, 'rows');
                    facet_section(5) = loc;
                    facet_section(6) = abs(dist);
                    facet_section(7:9) = BI;
                    facet_sections = [facet_sections; facet_section];
                     
%                     [id, loc] = ismember(xyz, fluid_IB_xyz_cut, 'rows');
%                     if id % the cell has already been logged as IB
%                         fluid_IB_facids{loc}(end+1) = facet;
%                         fluid_IB_facareas{loc}(end+1) = area;
%                         fluid_IB_dists{loc}(end+1) = dist;
%                     else % log the cell
%                         fluid_IB_xyz_cut = [fluid_IB_xyz_cut; xyz];
%                         fluid_IB_facids{end+1} = facet;
%                         fluid_IB_facareas{end+1} = area;
%                         fluid_IB_dists{end+1} = dist;
%                     end

%                     if (angles(id) ~= 1)
%                         stopflag = true;
%                     end

                    if stopflag
                        %figure
                        clf
                        view(3)
                        %patch('Faces', TR.ConnectivityList, 'Vertices', TR.Points, 'FaceColor', ones(3,1)*69/100, 'FaceAlpha', 0.5)
                        patch('Faces', [1 2 3], 'Vertices', TR.Points(TR.ConnectivityList(facet,:),:), 'FaceColor', [1,0,0], 'FaceAlpha', 0.5)
                        hold on
                        incenter = TR.incenter(1); faceNormal = TR.faceNormal;
                        quiver3(incenter(1), incenter(2), incenter(3), faceNormal(1), faceNormal(2), faceNormal(3))
                        V = [xl yl zl; xu yl zl; xu yu zl; xl yu zl; ...
                            xl yl zu; xu yl zu; xu yu zu; xl yu zu];
                        F = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
                        patch('Faces', F, 'Vertices', V, 'FaceColor', [1,1,1], 'FaceAlpha', 0.5)

                        scatter3(xgrid(i),ygrid(j),zgrid(k),20,[0,0,0],'filled')
                        scatter3(xyz(1),xyz(2),xyz(3),20,[0,0,0],'filled')
                        vec = xyz - BI;
                        quiver3(BI(:,1), BI(:,2), BI(:,3), vec(:,1), vec(:,2), vec(:,3),'off')
                        if ~fluid_IB(i,j,k)
                            xl = xyz(1)-dx/2; xu = xyz(1)+dx/2;
                            yl = xyz(2)-dy/2; yu = xyz(2)+dy/2;
                            zl = xyz(3)-dz/2; zu = xyz(3)+dz/2;
                            V = [xl yl zl; xu yl zl; xu yu zl; xl yu zl; ...
                                xl yl zu; xu yl zu; xu yu zu; xl yu zu];
                            patch('Faces', F, 'Vertices', V, 'FaceColor', [1,1,1], 'FaceAlpha', 0.5)
                        end
                        %patch('Faces', 1:size(clip,1), 'Vertices', clip, 'FaceColor', [0,0,1], 'FaceAlpha', 0.5)
                        trisurf(tri)
                        axis equal
                        %                     xlim([xgrid(il-2) xgrid(iu+2)])
                        %                     ylim([ygrid(jl-2) ygrid(ju+2)])
                        %                     zlim([zgrid(kl-2) zgrid(ku+2)])

                        xlabel('x')
                        ylabel('y')
                        zlabel('z')
                        %xlim([0 xgrid(end)])
                        %ylim([0 ygrid(end)])
                        %zlim([0 zgrid(end)])
                        drawnow
                        pause(1)
                        %return
                    end
                end
            end
        end
    end
end

% %%
% figure
% 
% %i=13; j=57; k=66;
% %facet = 895;
% 
% view(3)
% patch('Faces', [1 2 3], 'Vertices', TR.Points(TR.ConnectivityList(facet,:),:), 'FaceColor', [1,0,0], 'FaceAlpha', 0.5)
% %patch('Faces', TR.ConnectivityList, 'Vertices', TR.Points, 'FaceColor', ones(3,1)*69/100, 'FaceAlpha', 1)
% hold on
% 
% xl = xgrid(i) - dx/2;
% xu = xgrid(i) + dx/2;
% yl = ygrid(j) - dy/2;
% yu = ygrid(j) + dy/2;
% zl = zgrid(k) - dz/2;
% zu = zgrid(k) + dz/2;
% 
% V = [xl yl zl; xu yl zl; xu yu zl; xl yu zl; ...
%     xl yl zu; xu yl zu; xu yu zu; xl yu zu];
% F = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
% patch('Faces', F, 'Vertices', V, 'FaceColor', [1,0,1], 'FaceAlpha', 0.5)
% 
% if (solid(i,j,k))
% if i~=1
%     if fluid(i-1,j,k)
%         xls = xgrid(i-1)-dx/2; xus = xgrid(i-1)+dx/2;
%         yls = ygrid(j)-dy/2; yus = ygrid(j)+dy/2;
%         zls = zgrid(k)-dz/2; zus = zgrid(k)+dz/2;
%         V = [xls yls zls; xus yls zls; xus yus zls; xls yus zls; ...
%             xls yls zus; xus yls zus; xus yus zus; xls yus zus];
%         patch('Faces', F, 'Vertices', V, 'FaceColor', [1,1,1], 'FaceAlpha', 0.5)
%     end
% end
% if i~=itot
%     if fluid(i+1,j,k)
%         xls = xgrid(i+1)-dx/2; xus = xgrid(i+1)+dx/2;
%         yls = ygrid(j)-dy/2; yus = ygrid(j)+dy/2;
%         zls = zgrid(k)-dz/2; zus = zgrid(k)+dz/2;
%         V = [xls yls zls; xus yls zls; xus yus zls; xls yus zls; ...
%             xls yls zus; xus yls zus; xus yus zus; xls yus zus];
%         patch('Faces', F, 'Vertices', V, 'FaceColor', [1,1,1], 'FaceAlpha', 0.5)
%     end
% end
% if j~=1
%     if fluid(i,j-1,k)
%         xls = xgrid(i)-dx/2; xus = xgrid(i)+dx/2;
%         yls = ygrid(j-1)-dy/2; yus = ygrid(j-1)+dy/2;
%         zls = zgrid(k)-dz/2; zus = zgrid(k)+dz/2;
%         V = [xls yls zls; xus yls zls; xus yus zls; xls yus zls; ...
%             xls yls zus; xus yls zus; xus yus zus; xls yus zus];
%         patch('Faces', F, 'Vertices', V, 'FaceColor', [1,1,1], 'FaceAlpha', 0.5)
%     end
% end
% if j~=jtot
%     if fluid(i,j+1,k)
%         xls = xgrid(i)-dx/2; xus = xgrid(i)+dx/2;
%         yls = ygrid(j+1)-dy/2; yus = ygrid(j+1)+dy/2;
%         zls = zgrid(k)-dz/2; zus = zgrid(k)+dz/2;
%         V = [xls yls zls; xus yls zls; xus yus zls; xls yus zls; ...
%             xls yls zus; xus yls zus; xus yus zus; xls yus zus];
%         patch('Faces', F, 'Vertices', V, 'FaceColor', [1,1,1], 'FaceAlpha', 0.5)
%     end
% end
% if k~=1
%     if fluid(i,j,k-1)
%         xls = xgrid(i)-dx/2; xus = xgrid(i)+dx/2;
%         yls = ygrid(j)-dy/2; yus = ygrid(j)+dy/2;
%         zls = zgrid(k-1)-dz/2; zus = zgrid(k-1)+dz/2;
%         V = [xls yls zls; xus yls zls; xus yus zls; xls yus zls; ...
%             xls yls zus; xus yls zus; xus yus zus; xls yus zus];
%         patch('Faces', F, 'Vertices', V, 'FaceColor', [1,1,1], 'FaceAlpha', 0.5)
%     end
% end
% if k~=ktot
%     if fluid(i,j,k+1)
%         xls = xgrid(i)-dx/2; xus = xgrid(i)+dx/2;
%         yls = ygrid(j)-dy/2; yus = ygrid(j)+dy/2;
%         zls = zgrid(k+1)-dz/2; zus = zgrid(k+1)+dz/2;
%         V = [xls yls zls; xus yls zls; xus yus zls; xls yus zls; ...
%             xls yls zus; xus yls zus; xus yus zus; xls yus zus];
%         patch('Faces', F, 'Vertices', V, 'FaceColor', [1,1,1], 'FaceAlpha', 0.5)
%     end
% end
% end
% 
% %patch('Faces', 1:size(clip,1), 'Vertices', clip, 'FaceColor', [0,0,1], 'FaceAlpha', 0.5)
% axis equal
% % xlim([xgrid(il-2) xgrid(iu+2)])
% % ylim([ygrid(jl-2) ygrid(ju+2)])
% % zlim([zgrid(kl-2) zgrid(ku+2)])
% xlabel('x')
% ylabel('y')
% zlabel('z')
% pause(0.25)
% drawnow