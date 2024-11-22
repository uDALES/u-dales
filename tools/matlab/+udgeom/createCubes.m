function geom = createCubes(xsize, ysize, Hx, Hy, Hz, Cx, Cy, geom_option, edgelength)
% createCubes    creates cubes, either a single one or an array of cubes.
%
%    geom = createCubes(xsize, ysize, Hx, Hy, Hz, Cx, Cy, geom_option, edgelength)
%    returns a geom instance that can be saved to an stl file.
%        xsize:       length of the domain in x-direction
%        ysize:       length of the domain in y-direction
%        Hx:          cube x-length
%        Hy:          cube y-length
%        Hz:          cube height
%        Cx:          cube spacing in x-direction
%        Cy:          cube spacing in y-direction
%        geom_option: the type of geometry to create:
%                                'S': single cube
%                               'AC': aligned cubes
%                               'SC': staggered cubes
%        edgelength: the length of individual facets. Best taken as xsize
%                    (or ysize) divided by an integer number.

divisions = round(Hx/edgelength); % number of additional times to divide the cube (starting from once on each face)
tol_TR = 0; % this will shift the cubes slightly in the face normal direction - can potentially be removed if not needed

if strcmp(geom_option, 'SC')
    TR_og = generateStaggeredCubes(Hx, Hy, Hz, Cx, Cy, xsize, ysize, divisions, tol_TR);
elseif strcmp(geom_option, 'AC')
    TR_og = generateAlignedCubes(Hx, Hy, Hz, Cx, Cy, xsize, ysize, divisions, tol_TR);
elseif strcmp(geom_option, 'S')
    shift = [xsize/2 ysize/2 0]; % user-defined, here shift to centre of domain
    TR_og = generateSingleCube(Hx, Hy, Hz, shift, divisions, tol_TR);
end

stl_ground = true; % switch for adding ground facets

if stl_ground
    TR = generateGround(TR_og, xsize, ysize, edgelength);
else
    % remove any existing floor facets
    floor_ids = all(abs(TR_og.faceNormal - [0 0 -1]) <= 1e-8, 2); % facet normal in -z direction
    conn = TR_og.ConnectivityList;
    points = TR_og.Points;
    conn(floor_ids,:) = [];
    TR = triangulation(conn, points);
end

geom = udgeom.udgeom(TR);

% ---------------------------------------------------------------------- %

    function array_shifted = generateStaggeredCubes(Hx, Hy, Hz, Cx, Cy, Lx, Ly, divisions, tol)
        Nx = Lx / (Hx + Cx); Ny = Ly / (Hy + Cy);

        if mod(Nx,1) ~= 0 || mod(Ny,1) ~= 0
            error('The domain size should be a multiple of cube width + canyon width')
        end

        scale = [Hx Hy Hz];
        shift = [];

        for i=1:Nx
            if ~mod(i,2)
                shift = [shift; [i*(Cx+Hx)-Hx/2-Cx/2, 0, Hz/2]];
            end

            for j=1:Ny
                if mod(i,2)
                    shift = [shift; [i*(Cx+Hx)-Hx/2-Cx/2, j*(Cy+Hy)-Hy/2-Cy/2, Hz/2]];
                else
                    shift = [shift; [i*(Cx+Hx)-Hx/2-Cx/2, j*(Cy+Hy), Hz/2]];
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

        conn = array.ConnectivityList;

        % remove points and facets outside and on the edge of the domain
        if (divisions<2)
            for n=1:size(points,1)
                if points(n,2) < 0
                    points(n,2) = 0;
                elseif points(n,2) > Ly
                    points(n,2) = Ly;
                end
            end
        end

        points_remove = [];
        for n=1:size(points,1)
            if points(n,2) <= 0
                points_remove = [points_remove; n];
            elseif points(n,2) >= Ly
                points_remove = [points_remove; n];
            end
        end

        faces_remove = [];
        for n=1:size(array, 1)
            if all(ismember(array.ConnectivityList(n,:), points_remove))
                faces_remove = [faces_remove; n];
            end
        end
        conn(faces_remove,:) = [];


        [conn,points] = patchCleanUnused(conn,points);

        array_shifted = triangulation(conn, points);

    end

    % ----------------------------------------------------------------- %

    function array_shifted = generateAlignedCubes(Hx, Hy, Hz, Cx, Cy, Lx, Ly, divisions, tol)
        Nx = Lx / (Hx + Cx); Ny = Ly / (Hy + Cy);

        if mod(Nx,1) ~= 0 || mod(Ny,1) ~= 0
            error('The domain size should be a multiple of cube width + canyon width')
        end

        scale = [Hx Hy Hz];
        shift = [];

        for i=1:Nx
            for j=1:Ny
                shift = [shift; [i*(Cx+Hx)-Hx/2-Cx/2, j*(Cy+Hy)-Hy/2-Cy/2, Hz/2]];
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

        array_shifted = triangulation(array.ConnectivityList, points);

    end

    % ----------------------------------------------------------------- %

    function array_shifted = generateSingleCube(Hx, Hy, Hz, shift, divisions, tol)
        shift = shift + [0 0 Hz/2];
        scale = [Hx Hy Hz];

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
    end

    % ----------------------------------------------------------------- %

    function TR_out = operateUnitCube(scale, shift, divisions)
        % shift defines where the centre of the cubes are
        % scale defines how large they are

        TR_in = unitCube();

        % this cube is centred on the origin and has side length 1
        % scaling is done first, then shifting
        % for a regular array of cube size(s) size Hx, Hy, Hz
        % scale = [Hx Hy Hz], shift(n,:) = [2*n*Hx 2*n*Hy Hz/2];

        TR = TR_in;
        for n=1:divisions
            TR = divideFaces(TR);
        end

        cube1_vert = TR.Points;
        cube1_conn = TR.ConnectivityList;

        ncubes = size(shift,1);
        cubes_vert = [];
        cubes_conn = [];
        for n=1:ncubes
            cubes_vert = [cubes_vert; cube1_vert.*scale + shift(n,:)];
            cubes_conn = [cubes_conn; cube1_conn + size(TR.Points,1) * (n-1)];
        end

        TR_out = triangulation(cubes_conn, cubes_vert);

    end
end


