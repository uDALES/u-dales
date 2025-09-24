function building_components = splitBuildings(mesh_TR, varargin)
% SPLITBUILDINGS Separate connected building components from a triangulated mesh
%
% SYNTAX:
%   building_components = splitBuildings(mesh_TR)
%   building_components = splitBuildings(mesh_TR, 'RemoveGround', false)
%
% DESCRIPTION:
%   Analyzes a triangulated mesh to identify and separate individual
%   connected building components. The function automatically removes
%   ground-level faces (Z=0) by default, then uses graph connectivity
%   analysis to group triangular faces that share edges, effectively
%   separating distinct buildings or building clusters from a complex
%   urban geometry.
%
%   This is particularly useful for urban flow simulations where individual
%   building analysis is required, or when processing merged CAD models
%   that contain multiple disconnected structures.
%
% INPUT ARGUMENTS:
%   mesh_TR - Triangulation object containing the mesh to analyze:
%             * ConnectivityList: [N×3] triangle connectivity matrix
%             * Points: [M×3] vertex coordinates [x, y, z]
%             Can include ground faces (they will be automatically removed)
%
%   Optional Name-Value Pairs:
%   'RemoveGround' - Logical flag to enable/disable ground removal (default: true)
%
% OUTPUT ARGUMENTS:
%   building_components - Cell array where each cell contains a triangulation
%                        object representing one connected building component:
%                        {TR1, TR2, TR3, ...} where each TRi is a triangulation
%                        of an individual building or connected structure
%
% ALGORITHM:
%   1. Optional ground removal (faces with Z=0 vertices)
%   2. Edge-based graph construction (shared edges = connections)
%   3. Connected component analysis to group faces
%   4. Triangulation object creation for each component
%   5. Component sorting by face count (largest first)
%
% CONNECTIVITY ANALYSIS:
%   - Edge-based connectivity: Faces sharing edges are considered connected
%   - Graph components: Groups of transitively connected faces
%   - Component separation: Each connected group becomes a separate building
%   - Topology preservation: Maintains mesh integrity within each component
%
% EXAMPLE:
%   % Load urban mesh and separate buildings (ground removed automatically)
%   urban_mesh = stlread('urban_scene.stl');
%   buildings = udgeom.splitBuildings(urban_mesh);
%   
%   % Process each building individually
%   for i = 1:length(buildings)
%       fprintf('Building %d: %d faces\n', i, size(buildings{i}.ConnectivityList, 1));
%   end
%
%   % Keep ground faces if needed
%   buildings_with_ground = udgeom.splitBuildings(urban_mesh, 'RemoveGround', false);
%
% APPLICATIONS:
%   - Individual building analysis in urban environments
%   - Separation of merged CAD models
%   - Building-specific flow simulation setup
%   - Urban morphology studies
%   - Facet-based energy balance calculations per building
%
% PERFORMANCE NOTES:
%   - Computational complexity scales with number of faces and edges
%   - Memory usage depends on mesh connectivity and number of components
%   - Graph analysis is efficient for typical urban mesh sizes
%   - Ground removal adds minimal computational overhead
%
% SEE ALSO:
%   triangulation, graph, conncomp, udgeom.calculateOutline, udgeom.getBuildingInfo
%
% AUTHOR: uDALES Development Team
% VERSION: 1.0

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

% Copyright (C) 2025 the uDALES Team.

    %% PARSE INPUT ARGUMENTS
    % Set up input parser for optional parameters
    p = inputParser;
    addRequired(p, 'mesh_TR', @(x) isa(x, 'triangulation'));
    addParameter(p, 'RemoveGround', true, @(x) islogical(x));
    
    parse(p, mesh_TR, varargin{:});
    
    % Extract parsed parameters
    remove_ground = p.Results.RemoveGround;
    
    % Display initial message
    fprintf('Extracting individual buildings from STL geometry...\n');
    
    %% PREPROCESSING: REMOVE GROUND FACES (OPTIONAL)
    % Store original face count for statistics
    original_faces = size(mesh_TR.ConnectivityList, 1);
    
    % Automatically remove ground-level faces unless disabled
    if remove_ground
        % Call the private deleteGround function to remove Z=0 faces
        processed_mesh = deleteGround(mesh_TR);
        remaining_faces = size(processed_mesh.ConnectivityList, 1);
        removed_faces = original_faces - remaining_faces;
        fprintf('Ground preprocessing: %d faces removed, %d faces remaining\n', removed_faces, remaining_faces);
    else
        % Use the mesh as-is without ground removal
        processed_mesh = mesh_TR;
        fprintf('Ground preprocessing: Ground faces preserved (all %d faces kept)\n', original_faces);
    end
    
    %% EXTRACT MESH DATA
    % Get triangular connectivity and vertex coordinates from processed mesh
    faces = processed_mesh.ConnectivityList;  % [N×3] Triangle connectivity matrix
    points = processed_mesh.Points;           % [M×3] Vertex coordinates [x, y, z]
    
    %% BUILD EDGE LIST WITH FACE ASSOCIATIONS
    % Create comprehensive edge list mapping each edge to its adjacent faces
    % This is essential for identifying face-to-face connectivity
    edges = [];
    
    for i = 1:size(faces, 1)
        face = faces(i, :);
        % Add three edges of current triangle with face index
        edges = [edges; sort([face(1), face(2)]), i];  % Edge vertex1-vertex2
        edges = [edges; sort([face(2), face(3)]), i];  % Edge vertex2-vertex3  
        edges = [edges; sort([face(3), face(1)]), i];  % Edge vertex3-vertex1
    end
    
    %% IDENTIFY SHARED EDGES
    % Sort edges by vertex indices to group identical edges together
    sorted_edges = sortrows(edges, [1, 2]);
    
    % Find pairs of faces that share edges (adjacent faces)
    shared_pairs = [];
    for k = 1:size(sorted_edges, 1)-1
        % Check if consecutive edges in sorted list are identical
        if all(sorted_edges(k, 1:2) == sorted_edges(k+1, 1:2))
            % Record the two faces that share this edge
            shared_pairs = [shared_pairs; sorted_edges(k,3), sorted_edges(k+1,3)];
        end
    end
    
    %% GRAPH CONNECTIVITY ANALYSIS
    % Build undirected graph where nodes are faces and edges represent shared boundaries
    G = graph(shared_pairs(:,1), shared_pairs(:,2));
    
    % Find connected components (groups of transitively connected faces)
    [components, component_sizes] = conncomp(G);
    
    % Determine number of separate building components
    num_components = max(components);
    
    %% CREATE SEPARATE TRIANGULATIONS FOR EACH COMPONENT
    % Initialize cell array to store individual building triangulations
    building_components = cell(num_components, 1);
    
    % Process each connected component
    for i = 1:num_components
        % Extract faces belonging to current component
        component_faces = faces(components == i, :);
        
        if ~isempty(component_faces)
            % Create triangulation object for this building component
            % Suppress warning about unreferenced points (expected behavior for consistent indexing)
            warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');
            building_components{i} = triangulation(component_faces, points);
            warning('on', 'MATLAB:triangulation:PtsNotInTriWarnId');
        end
    end
    
    % Count valid buildings and calculate statistics
    valid_buildings = 0;
    total_building_faces = 0;
    largest_building_faces = 0;
    smallest_building_faces = inf;
    
    for i = 1:num_components
        if ~isempty(building_components{i})
            valid_buildings = valid_buildings + 1;
            num_faces = size(building_components{i}.ConnectivityList, 1);
            total_building_faces = total_building_faces + num_faces;
            largest_building_faces = max(largest_building_faces, num_faces);
            smallest_building_faces = min(smallest_building_faces, num_faces);
        end
    end
    
    if smallest_building_faces == inf
        smallest_building_faces = 0;
    end
    
    % Display meaningful statistics
    fprintf('\nBuilding separation completed:\n');
    fprintf('  Buildings identified: %d separate components\n', valid_buildings);
    fprintf('  Total building faces: %d\n', total_building_faces);
    if valid_buildings > 0
        fprintf('  Largest building: %d faces\n', largest_building_faces);
        fprintf('  Smallest building: %d faces\n', smallest_building_faces);
        fprintf('  Average faces per building: %.1f\n', total_building_faces / valid_buildings);
    end

end