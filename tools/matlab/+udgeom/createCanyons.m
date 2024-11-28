function geom = createCanyons(xsize, ysize, B, W, H, shift, edgelength, rotate90)
% createCanyons    creates one-dimensional street canyons
%
%    geom = createCanyons(xsize, ysize, B, W, H, shift, edgelength) returns
%    a geom instance that can be saved to an stl file.
%        xsize:      length of the domain in x-direction
%        ysize:      length of the domain in y-direction
%        B:          building width
%        W:          street width
%        H:          building height
%        shift:      shifts the canyons to the right in the x-direction 
%        edgelength: the length of individual facets. Best taken as xsize
%                    (or ysize) divided by an integer number.
%        rotate90:   boolean that allows one to rotate the domain once it
%                    has been generated. Set this parameter to false under
%                    normal conditions.

L = (B+W); % canyon length

divisions = round(B/edgelength);
TR = generateCanyons(B, W, H, L, xsize, ysize, divisions, shift);

if rotate90
    points = TR.Points;
    angle = 90;
    % Rotate about the z axis
    Rz = [ cosd(angle), -sind(angle), 0 ;
        sind(angle), cosd(angle),  0 ;
        0,           0,            1 ];

    points = points * Rz' + [xsize 0 0];
    TR = triangulation(TR.ConnectivityList, points);
end

geom = udgeom.udgeom(TR);

% ----------------------------------------------------------------------%

    function TR = generateCanyons(B, W, H, L, Lx, Ly, divisions, shift)

        if divisions > 2
            error('divisions must be <= 2')
        end

        if divisions==2
            L = L/2;
        end

        Nx = Lx/(B+W);
        Ny = Ly/L;

        if mod(Nx,1) ~= 0 || mod(Ny,1) ~= 0
            error('The domain size should be a multiple of canyon width')
        end

        % Generate base canyon
        points_floor1 = [0 0 0; W/2 0 0; W/2 L 0; 0 L 0];
        points_floor2 = points_floor1 + [W/2+B 0 0];
        points_wall1 = [0 0 0; 0 0 H; 0 L H; 0 L 0] + [W/2 0 0];
        points_wall2 = flip(points_wall1,1) + [B 0 0];
        points_roof = [0 0 0; B 0 0; B L 0; 0 L 0] + [0 0 H] + [W/2 0 0];


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

        xsize = max(points(:,1));
        if shift > 0
            for n=1:size(points,1)
                if (points(n,1) > 0) && (points(n,1) < xsize)
                    points(n,1) = points(n,1) + shift;
                end
            end
        end

        warning('off')
        TR = triangulation(conn, points);
        warning('on')
        if divisions==1
            TR = divideFaces(TR);
        elseif divisions==2
            TR = divideFaces(TR);
        end

    end
end
