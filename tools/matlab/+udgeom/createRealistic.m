function geom = createRealistic(stlfile, xsize, ysize, shift, edgelength)
% createRealistic    creates a realistic urban surface based on an stl file
%                    that contains the buildings. The function adds the
%                    ground surface.  
%
%    geom = createRealistic(stlfile, xsize, ysize, shift, edgelength)
%    returns a geom instance that can be saved to an stl file. 
%        stlfile:    the STL file that contains the buildings (NOT the
%                    ground)
%        xsize:      length of the new domain in x-direction
%        ysize:      length of the new domain in y-direction
%        shift:      array that shifts the geometry. shift[1], shift[2] and
%                    shift[3] represent, respectively, the shift in x-, y-
%                    and z-direction.
%        edgelength: the length of individual facets. Best taken as xsize
%                    (or ysize) divided by an integer number.

TR_og = stlread(stlfile);

% translation, e.g. to make domain larger in x direction
points = TR_og.Points + shift;
TR_shifted = triangulation(TR_og.ConnectivityList, points);

stl_ground = true;

if stl_ground
    % Add ground facets
    TR = generateGround(TR_shifted, xsize, ysize, edgelength);
else
    % remove any existing floor facets
    floor_ids = all(abs(TR_shifted.faceNormal - [0 0 -1]) <= 0, 2); % facet normal in -z direction
    conn = TR_shifted.ConnectivityList;
    conn(floor_ids,:) = [];
    warning('off')
    TR = triangulation(conn, points);
    warning('on')
end

geom = udgeom.udgeom(TR);