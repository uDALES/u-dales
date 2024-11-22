function TR = generateCanyons(B, W, H, L, Nx, Ny)
    % Generate base canyon
    points_floor1 = [0 0 0; W/2 0 0; W/2 L 0; 0 L 0];
    points_floor2 = points_floor1 + [W/2+B 0 0];
    points_wall1 = [0 0 0; 0 0 H; 0 L H; 0 L 0] + [W/2 0 0];
    points_wall2 = flip(points_wall1,1) + [B 0 0];
    points_roof = [0 0 0; B 0 0; B L 0; 0 L 0] + [0 0 H] + [W/2 0 0];

%     figure
%     view(3)
%     axis equal
%     patch('XData', points_floor1(:,1), 'YData', points_floor1(:,2), 'ZData', points_floor1(:,3))
%     patch('XData', points_floor2(:,1), 'YData', points_floor2(:,2), 'ZData', points_floor2(:,3))
%     patch('XData', points_wall1(:,1), 'YData', points_wall1(:,2), 'ZData', points_wall1(:,3))
%     patch('XData', points_wall2(:,1), 'YData', points_wall2(:,2), 'ZData', points_wall2(:,3))
%     patch('XData', points_roof(:,1), 'YData', points_roof(:,2), 'ZData', points_roof(:,3))

    %conn_unit = [1 2 3; 1 3 4];
    points_unit = [points_floor1; points_wall1; points_roof; points_wall2; points_floor2];
    conn = [];
    points = [];

    for i=1:Nx
        for j=1:Ny
            for n=1:5
                locs = zeros(1,4);
                for m=1:4
                    %loc = 4*(n-1) + m;
                    point = points_unit(4*(n-1) + m, :) + [(i-1)*(B+W) (j-1)*L 0];
                    if isempty(points)
                        points = [points; point];
                        loc = 1;
                    else
                        [lia, locb] = ismember(point, points, 'rows'); % find first match
                        if ~lia % it is not already in the list
                            points = [points; point];
                            loc = size(points, 1);
                        else % it is, so use the index
                            loc = locb;
                        end
                    end
                    locs(m) = loc;
                end
                conn = [conn; [locs(1) locs(2) locs(3); locs(1) locs(3) locs(4)]];
            end
        end
    end

    conn
    points
    TR = triangulation(conn, points);
    figure
    trisurf(TR);
    axis equal

    % conn = conn_unit;
    % points = points_unit;
    % for i=1:Nx
    %     for j=1:Ny
    %         if (i==1 && j==1) 
    %             continue
    %         end
    %         points = [points; points_unit + [(i-1)*(B+W) (j-1)*L]];
    %         conn = [conn; 1 2 3; 1 3 4] + max(conn, [], 'all')]

    %end
end
