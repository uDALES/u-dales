% udgeom Geometry class for uDALES
%    The udgeom class contains the triangulated surface.

% uDALES (https://github.com/uDALES/u-dales).

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Copyright (C) 2016-2021 the uDALES Team.

classdef udgeom < handle
   % Geometry class for uDALES

   properties (Hidden = true, SetAccess = protected)
      path;                    % Path to geometry.
      cpath;                   % Current path.
   end

   properties (SetAccess = protected)
      stl;                     % stl of the geometry
      outline;                 % outline edges calculated using calculateOutline
      outline2d;               % cached 2D polygon outlines and centroids (cell array of structs with 'polygon' and 'centroid' fields)
      building_outlines;       % cached 3D outline edges for individual buildings (cell array of edge arrays, lazy loaded)
      buildings;               % cell array of individual building triangulations (lazy loaded)
      face_to_building_map;    % array mapping each face index to its building ID (lazy loaded)
   end   
   
   methods
      % geometry constructor
      
      function obj = udgeom(varargin)
         % Class constructor.
         %
         % varargin : Path for geometry
         %            Or a triangulation object
         
         % store current and simulation directory
         obj.cpath = pwd;
         
         dapath = obj.cpath;
         if (nargin == 1)
             % Check if the input argument is a string
                if ischar(varargin{1}) || isstring(varargin{1})
                    dapath = varargin{1};
                % Check if the input argument is a triangulation object
                elseif isa(varargin{1}, 'triangulation')
                    obj.stl = varargin{1};
                    % Calculate and store outline edges with default angle threshold
                    [obj.outline, ~] = udgeom.calculateOutline(obj.stl, 45);
                else
                    error('Input must be either a string or a triangulation object.');
                end
         end
         cd([dapath]);
         obj.path = pwd;
      end
      
      % ------------------------------------------------------------- %
      
      function gopath(obj)
         % Go to simulation path.
         % gopath(obj)
         
         cd(obj.path);
      end
      
      % ------------------------------------------------------------- %
      
      function gohome(obj)
         % Go to work path.
         % gohome(obj)
         
         cd(obj.cpath);
      end
      
      % ------------------------------------------------------------- %
      
      function chcpath(obj, newpath)
         
         % Change work path.
         % chcpath(obj, newpath)
         
         here = pwd;
         cd(newpath)
         obj.cpath = pwd;
         cd(here);
      end
      
      % -------------------------------------------------------------- %
      
      function load(obj, filename)
         % load an STL file.
         % 
         % example:
         %   obj.load(obj, filename)
          
         obj.gopath()
         obj.stl = stlread(filename);
         obj.gohome()
         
         % Calculate and store outline edges with default angle threshold
         if ~isempty(obj.stl)
             [obj.outline, ~] = udgeom.calculateOutline(obj.stl, 45);
         end
        % invalidate cached 2D outlines and buildings mapping
        obj.outline2d = [];
        obj.building_outlines = [];
        obj.buildings = [];
        obj.face_to_building_map = [];
      end
      
      % -------------------------------------------------------------- %
      
      function save(obj, filename)
         % save an STL file.
         %
         % example:
         %   obj.save(obj, filename)

         obj.gopath()
         stlwrite(obj.stl, filename)
         obj.gohome()
      end
      
      % -------------------------------------------------------------- %
      
      function building_components = get_buildings(obj)
         % Get individual building components with lazy loading.
         %
         % building_components = get_buildings(obj)
         %   Returns cell array of triangulation objects, one for each building
         %
         % Example:
         %   buildings = geom.get_buildings();
         %   num_buildings = length(buildings);
         
         if isempty(obj.stl)
             warning('No STL geometry loaded. Load geometry first.');
             building_components = {};
             return;
         end
         
         % Lazy load buildings if not already computed
            if isempty(obj.buildings)
                [obj.buildings, obj.face_to_building_map] = udgeom.splitBuildings(obj.stl);
                % invalidate cached 2D and 3D outlines when buildings are (re)computed
                obj.outline2d = [];
                obj.building_outlines = [];
            end
         
         building_components = obj.buildings;
      end
      
      % -------------------------------------------------------------- %
      
      function face_map = get_face_to_building_map(obj)
         % Get mapping from face indices to building IDs.
         %
         % face_map = get_face_to_building_map(obj)
         %   Returns array where face_map(i) is the building ID for face i
         %
         % Example:
         %   face_map = geom.get_face_to_building_map();
         %   building_1_faces = find(face_map == 1);
         
         if isempty(obj.stl)
             warning('No STL geometry loaded. Load geometry first.');
             face_map = [];
             return;
         end
         
         % Ensure buildings and mapping are computed
         if isempty(obj.face_to_building_map)
             obj.get_buildings(); % This will compute both buildings and mapping
         end
         
         face_map = obj.face_to_building_map;
      end
      
      % -------------------------------------------------------------- %
      
      function outlines = get_building_outlines(obj, angle_threshold)
         % Get 3D outline edges for individual buildings with lazy loading.
         %
         % outlines = get_building_outlines(obj)
         %   Returns cell array of outline edge arrays, one per building
         %   Uses default angle threshold of 45 degrees
         %
         % outlines = get_building_outlines(obj, angle_threshold)
         %   Uses custom angle threshold for edge detection
         %
         % Each element of outlines{i} is an Nx2 array where each row contains
         % the indices of two vertices that form an outline edge for building i.
         %
         % Example:
         %   outlines = geom.get_building_outlines();
         %   edges_building_1 = outlines{1};
         %
         % SEE ALSO: calculateOutline, get_buildings, calculate_outline2d
         
         if isempty(obj.stl)
             warning('No STL geometry loaded. Load geometry first.');
             outlines = {};
             return;
         end
         
         if nargin < 2
             angle_threshold = 45; % Default angle threshold
         end
         
         % Lazy load building outlines if not already computed
         if isempty(obj.building_outlines)
             % Get individual buildings
             blds = obj.get_buildings();
             if isempty(blds)
                 outlines = {};
                 return;
             end
             
             % Calculate outline for each building
             num_buildings = length(blds);
             obj.building_outlines = cell(num_buildings, 1);
             
             for i = 1:num_buildings
                 b = blds{i};
                 if isempty(b)
                     obj.building_outlines{i} = [];
                     continue;
                 end
                 
                 % Accept both new struct-wrapped triangulation and raw triangulation
                 if isstruct(b) && isfield(b, 'triangulation')
                     TR = b.triangulation;
                 elseif isa(b, 'triangulation')
                     TR = b;
                 else
                     obj.building_outlines{i} = [];
                     continue;
                 end
                 
                 % Calculate outline edges for this building using calculateOutline
                 try
                     [edges, ~] = udgeom.calculateOutline(TR, angle_threshold);
                     obj.building_outlines{i} = edges;
                 catch
                     obj.building_outlines{i} = [];
                 end
             end
         end
         
         outlines = obj.building_outlines;
      end
      
      % -------------------------------------------------------------- %
      
      function calculate_outline_edges(obj, angle_threshold)
         % Recalculate outline edges with a custom angle threshold.
         %
         % calculate_outline_edges(obj, angle_threshold)
         %   angle_threshold: Angle threshold in degrees for detecting sharp edges
         %
         % Example:
         %   obj.calculate_outline_edges(30);  % More aggressive edge detection
         
         if isempty(obj.stl)
             warning('No STL geometry loaded. Load geometry first.');
             return;
         end
         
         if nargin < 2
             angle_threshold = 45; % Default angle threshold
         end
         
         [obj.outline, ~] = udgeom.calculateOutline(obj.stl, angle_threshold);
      end

      % -------------------------------------------------------------- %

      function outline2d = calculate_outline2d(obj)
         % Calculate 2D building outlines and centroids for all buildings.
         % outlines = calculate_outline2d(obj)
         %   Returns a cell array where each entry is a struct with fields:
         %     polygon : Nx2 array of boundary points defining the 2D polygon (may be empty)
         %     centroid : 1x2 array of the polygon centroid

         % If we've already computed outlines, return cached value
         if ~isempty(obj.outline2d)
             outline2d = obj.outline2d;
             return;
         end

         % Initialize
         outline2d = {};
         blds = obj.get_buildings();
         if isempty(blds)
             return;
         end

         num_buildings = length(blds);
         for i = 1:num_buildings
             b = blds{i};
             if isempty(b)
                 outline2d{i} = struct('polygon', [], 'centroid', [NaN NaN]);
                 continue;
             end

             % Accept both new struct-wrapped triangulation and raw triangulation
             if isstruct(b) && isfield(b, 'triangulation')
                 TR = b.triangulation;
             elseif isa(b, 'triangulation')
                 TR = b;
             else
                 outline2d{i} = struct('polygon', [], 'centroid', [NaN NaN]);
                 continue;
             end

             pts = TR.Points;
             faces = TR.ConnectivityList;

             % Basic validation
             if isempty(pts) || size(pts,1) < 3 || isempty(faces)
                 % Ensure centroid is a 1x2 row vector
                 if isempty(pts)
                     centroid = [NaN NaN];
                 else
                     centroid = mean(pts(:,1:2), 1);
                     centroid = centroid(:)';  % Force row vector
                 end
                 outline2d{i} = struct('polygon', [], 'centroid', centroid);
                 continue;
             end

             % Project to 2D and attempt to get the external boundary using freeBoundary
             pts2 = pts(:,1:2);
             try
                 tri2 = triangulation(faces, pts2);
                 boundary_edges = freeBoundary(tri2);

                 if ~isempty(boundary_edges)
                     % Use indices returned by freeBoundary to select boundary points
                     % (matches behaviour in the last committed plot_building_ids)
                     polygon = tri2.Points(boundary_edges(:,1), :);
                     centroid = mean(polygon, 1);
                 else
                     % Fallback: no free boundary found -> use centroid of all points
                     polygon = [];
                     centroid = mean(pts2, 1);
                 end
             catch
                 % On error (degenerate geometry etc.) fall back to centroid
                 polygon = [];
                 centroid = mean(pts2, 1);
             end

             % Ensure centroid is always a 1x2 row vector
             centroid = centroid(:)';
             outline2d{i} = struct('polygon', polygon, 'centroid', centroid);
         end

         % Cache result on the object for faster repeated access
         obj.outline2d = outline2d;
      end
      
      % -------------------------------------------------------------- %
      
      function show(obj, varargin)
         % plot the geometry
         %
         % show(obj, colorbuildings)
         %       colorbuildings (optional): boolean parameter on whether
         %                                  to colour buildings. This
         %                                  parameter is true by default.
         %                                  Needs to be set to false for 
         %                                  large geometries.
         %
         % examples:
         %   obj.show();
         %   obj.show(false);

         color_buildings = true;
         if ~isempty(varargin)
             color_buildings = varargin{1};
         end

         faceNormals = faceNormal(obj.stl);
         incenters = incenter(obj.stl);
         
         figure
         patch('Faces', obj.stl.ConnectivityList, 'Vertices', obj.stl.Points, ...
               'FaceColor', ones(3,1)*0.85, 'FaceAlpha', 1) 
         hold on
         
         if (color_buildings)
             for i=1:length(incenters(:,1))
                 if incenters(i,3)>0
                     A = obj.stl.Points(obj.stl.ConnectivityList(i,1),:);
                     B = obj.stl.Points(obj.stl.ConnectivityList(i,2),:);
                     C = obj.stl.Points(obj.stl.ConnectivityList(i,3),:);
                     fill3([A(1) B(1) C(1)], [A(2) B(2) C(2)], [A(3) B(3) C(3)], ...
                         [0.73 0.83 0.96])
                 end
             end
         end
         
         quiver3(incenters(:,1), incenters(:,2), incenters(:,3), ...
                 faceNormals(:,1), faceNormals(:,2), faceNormals(:,3), 0.2)
         view(3)
         axis equal tight
         set(gca,"Box","on")
         set(gca,"BoxStyle","full")
         hold off
      end
      
      % -------------------------------------------------------------- %
      
      function show_outline(obj, angle_threshold)
         % Plot the geometry outline edges
         %
         % show_outline(obj) plots the precomputed outline edges of the geometry
         % show_outline(obj, angle_threshold) recalculates outline with custom threshold
         %
         % Parameters:
         %   angle_threshold (optional): Angle threshold in degrees for edge detection
         %                              If not provided, uses precomputed outline
         %
         % Examples:
         %   obj.show_outline();     % Use precomputed outline
         %   obj.show_outline(30);   % Recalculate with 30° threshold

         if isempty(obj.stl)
             warning('No STL geometry loaded. Load geometry first.');
             return;
         end

         % Handle optional angle_threshold argument
         if nargin < 2
             angle_threshold = [];
         end

         % Determine which outline to use
         if isempty(angle_threshold)
             % Use precomputed outline if available
             if isempty(obj.outline)
                 % Calculate with default threshold if not available
                 obj.calculate_outline_edges(45);
             end
             outline_edges = obj.outline;
         else
             % Recalculate with custom threshold
             [outline_edges, ~] = udgeom.calculateOutline(obj.stl, angle_threshold);
         end

         if isempty(outline_edges)
             warning('No outline edges found.');
             return;
         end

         % Create figure and plot mesh and outline
         figure
         hold on;

         % Plot the mesh without edge lines
         patch('Faces', obj.stl.ConnectivityList, 'Vertices', obj.stl.Points, ...
               'FaceColor', [0.85 0.85 0.85], ...
               'EdgeColor', 'none', ...
               'FaceAlpha', 0.8);

         % Get points
         points = obj.stl.Points;

         % Prepare line segment coordinates for efficient plotting
         n_edges = size(outline_edges, 1);
         coords = nan(3*n_edges, 3);

         % Fill coordinate array (point1, point2, NaN for each edge)
         for i = 1:n_edges
             idx = 3*(i-1) + 1;
             coords(idx,:)   = points(outline_edges(i,1),:);
             coords(idx+1,:) = points(outline_edges(i,2),:);
             % coords(idx+2,:) remains NaN for line separation
         end

         % Plot outline edges with thinner lines
         line(coords(:,1), coords(:,2), coords(:,3), ...
              'Color', 'k', ...
              'LineStyle', '-');

         % Set up the plot
         view(3);
         axis equal tight;
         grid on;
         xlabel('x [m]');
         ylabel('y [m]');
         zlabel('z [m]');
         set(gca, 'Box', 'on');
         set(gca, 'BoxStyle', 'full');
         
         hold off;
      end
   end
end
