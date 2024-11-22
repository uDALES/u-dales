%% MAKE SURE TO CHANGE PATHS

function TR = generateGround(TR_noground, xsize, ysize, edgelength)
    % Add facets at ground level
    if size(TR_noground,1) > 1
        ground_ogptids = find(TR_noground.Points(:,3) == 0);
        ground_points = TR_noground.Points(ground_ogptids, 1:2);

        TR_noground_edges = edges(TR_noground);
        is_ground_edge = ismember(TR_noground_edges(:, 1), ground_ogptids) & ...
                        ismember(TR_noground_edges(:, 2), ground_ogptids);

        filtered_edges = TR_noground_edges(is_ground_edge, :);
        [~, ground_constraints] = ismember(filtered_edges, ground_ogptids);
    else
        ground_points = [];
        ground_constraints = [];
    end

    % add regular points on the edge
    % always start at the origin - can be changed for special simulations
    x0 = 0;
    y0 = 0;

    nedges_x = xsize / edgelength;
    nedges_y = ysize / edgelength;

    % Edge of domain - contrained (perhaps not necessary in some cases but makes sure domain edge is defined)
    ground_points_edges = [];
    ground_constraints_edges = [];
    % left to right
    for i=0:nedges_x
        ground_points_edges = [ground_points_edges; x0+i*edgelength y0];
        nground_points_edges = size(ground_points_edges,1);
        ground_constraints_edges = [ground_constraints_edges; [nground_points_edges nground_points_edges+1]];
    end
    for i=1:nedges_y
        ground_points_edges = [ground_points_edges; x0+xsize y0+i*edgelength];
        nground_points_edges = size(ground_points_edges,1);
        ground_constraints_edges = [ground_constraints_edges; [nground_points_edges nground_points_edges+1]];
    end
    for i=1:nedges_x
        ground_points_edges = [ground_points_edges; x0+xsize-i*edgelength y0+ysize];
        nground_points_edges = size(ground_points_edges,1);
        ground_constraints_edges = [ground_constraints_edges; [nground_points_edges nground_points_edges+1]];
    end
    for i=1:nedges_y-1
        ground_points_edges = [ground_points_edges; x0 y0+ysize-i*edgelength];
        if i<nedges_y-1
            nground_points_edges = size(ground_points_edges,1);
            ground_constraints_edges = [ground_constraints_edges; [nground_points_edges nground_points_edges+1]];
        end
    end
    ground_constraints_edges = [ground_constraints_edges; [nground_points_edges+1 1]];

    % interior points - not contrained
    ground_points_inside = [];
    for i=1:nedges_x
        for j=1:nedges_y
            point_inside = [x0+i*edgelength y0+j*edgelength];
            ground_points_inside = [ground_points_inside; point_inside];
        end
    end

    points_DT = [ground_points; ground_points_edges; ground_points_inside];
    constraints_DT = [ground_constraints; ground_constraints_edges+size(ground_points,1)];
    warning('off')
    DT = delaunayTriangulation(points_DT, constraints_DT);
    warning('on')

    % figure
    % scatter(points_DT(:,1), points_DT(:,2))
    % hold on
    % triplot(DT)
    % axis equal tight
    
    % Make 3D triangulation
    TR_ground_points = [DT.Points, zeros(size(DT.Points, 1), 1)];
    TR_ground_3D = triangulation(DT.ConnectivityList, TR_ground_points);
    % figure
    % trisurf(TR_ground_3D)

    if size(TR_noground,1) > 1
        % Merge triangulations
        [F,V] = joinElementSets({TR_noground.ConnectivityList, TR_ground_3D.ConnectivityList}, ...
            {TR_noground.Points, TR_ground_3D.Points});
        [F,V] = mergeVertices(F,V);
        TR_combined = triangulation(F,V);

        % figure
        % trisurf(TR_combined)

        % Remove facets below buildings
        ground_facids = all(ismember(TR_combined.ConnectivityList, find(TR_combined.Points(:,3) == 0)), 2);
        in = false(size(F, 1),1);
        in(ground_facids) = inpolyhedron(F, V, TR_combined.incenter(find(ground_facids)));
        F(in,:) = [];

        % Remove facets with repeated vertex ids
        n_remove = [];
        for n=1:size(F,1)
            if size(unique(V(F(n,:),:), 'rows'), 1) < 3
                n_remove = [n_remove, n];
            end
        end

        F(n_remove,:) = [];

        % Remove unused vertices
        [F,V] = patchCleanUnused(F,V);
        TR = triangulation(F,V);
    else
        TR = TR_ground_3D;
    end

end