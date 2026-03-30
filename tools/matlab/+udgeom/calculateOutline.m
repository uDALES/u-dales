function [boundary_edges, mesh_data] = calculateOutline(TR, angle_threshold)
% CALCULATEOUTLINE Calculate boundary and sharp edges for mesh outlining
%
% SYNTAX:
%   [boundary_edges, mesh_data] = udgeom.calculateOutline(TR)
%   [boundary_edges, mesh_data] = udgeom.calculateOutline(TR, angle_threshold)
%
% DESCRIPTION:
%   Analyzes a triangulated mesh to identify edges that should be outlined,
%   including both true boundary edges and sharp geometric features based
%   on face normal angles.
%
% INPUT ARGUMENTS:
%   TR              - Triangulation object containing:
%                     * ConnectivityList: [N×3] triangle connectivity matrix
%                     * Points: [M×3] vertex coordinates [x, y, z]
%   angle_threshold - (Optional) Angle threshold in degrees for detecting
%                     sharp edges. Default: 45 degrees
%
% OUTPUT ARGUMENTS:
%   boundary_edges  - [K×2] Array of vertex pairs defining edges to outline
%   mesh_data       - Structure containing mesh information:
%                     * faces: Triangle connectivity matrix
%                     * vertices: Vertex coordinates
%                     * normals: Face normal vectors
%                     * edge_map: Map from edges to adjacent faces
%
% ALGORITHM:
%   1. Computes face normal vectors for each triangle
%   2. Maps each edge to its adjacent faces
%   3. Identifies boundary edges (1 adjacent face)
%   4. Identifies sharp edges (angle between normals > threshold)
%
% EXAMPLE:
%   % Basic usage
%   [edges, data] = udgeom.calculateOutline(triangulation_obj);
%   
%   % With custom threshold
%   [edges, data] = udgeom.calculateOutline(triangulation_obj, 30);
%
% SEE ALSO:
%   udbase.plot_outline, plot_mesh_with_outline
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

    % Set default angle threshold if not provided
    if nargin < 2 || isempty(angle_threshold)
        angle_threshold = 45;  % Degrees
    end
    
    % Extract mesh data from triangulation object
    faces = TR.ConnectivityList;  % Triangle connectivity [N×3]
    vertices = TR.Points;         % Vertex coordinates [M×3]
    
    %% STEP 1: COMPUTE FACE NORMAL VECTORS
    % Calculate unit normal vectors for each triangular face
    % Normal direction determines face orientation for edge angle calculations
    normals = zeros(size(faces,1), 3);
    for i = 1:size(faces,1)
        % Get three vertices of the triangle
        v1 = vertices(faces(i,1),:);
        v2 = vertices(faces(i,2),:);
        v3 = vertices(faces(i,3),:);
        
        % Compute normal using cross product of two edge vectors
        normal_vec = cross(v2-v1, v3-v1);
        normals(i,:) = normal_vec / norm(normal_vec);  % Normalize to unit vector
    end
    
    %% STEP 2: BUILD EDGE-TO-FACE MAPPING
    % Create a map from each edge to the faces that contain it
    % This allows identification of boundary edges (1 face) and internal edges (2+ faces)
    edge_map = containers.Map;
    
    for i = 1:size(faces,1)
        % Extract the three edges of the current triangle
        % Sort vertex indices to ensure consistent edge representation
        edges = [sort([faces(i,1), faces(i,2)]);    % Edge 1-2
                 sort([faces(i,2), faces(i,3)]);    % Edge 2-3
                 sort([faces(i,3), faces(i,1)])];   % Edge 3-1
        
        % Add current face to each edge's face list
        for j = 1:3
            key = sprintf('%d-%d', edges(j,1), edges(j,2));
            if isKey(edge_map, key)
                edge_map(key) = [edge_map(key), i];  % Append face ID
            else
                edge_map(key) = [i];                 % Initialize with first face ID
            end
        end
    end
    
    %% STEP 3: IDENTIFY BOUNDARY AND SHARP EDGES
    % Determine which edges should be outlined based on two criteria:
    % 1. Boundary edges: belong to only one face (mesh boundary)
    % 2. Sharp edges: angle between adjacent face normals > threshold
    
    boundary_edges = [];
    keys = edge_map.keys();
    
    for k = 1:length(keys)
        face_ids = edge_map(keys{k});
        
        if length(face_ids) == 1
            % TRUE BOUNDARY EDGE: belongs to only one face
            boundary_edges = [boundary_edges; sscanf(keys{k}, '%d-%d')'];
            
        else
            % INTERNAL EDGE: check if it's a sharp geometric feature
            % Calculate maximum angle between any pair of adjacent faces
            max_angle = 0;
            for i = 1:length(face_ids)-1
                for j = i+1:length(face_ids)
                    n1 = normals(face_ids(i),:);
                    n2 = normals(face_ids(j),:);
                    
                    % Calculate angle between normal vectors
                    cos_theta = dot(n1, n2);
                    cos_theta = min(max(cos_theta, -1), 1);  % Clamp to [-1,1] for numerical stability
                    angle = acosd(cos_theta);
                    max_angle = max(max_angle, angle);
                end
            end
            
            % Mark as boundary if angle exceeds threshold (sharp feature)
            if max_angle > angle_threshold
                boundary_edges = [boundary_edges; sscanf(keys{k}, '%d-%d')'];
            end
        end
    end
    
    % Package mesh data for output
    mesh_data.faces = faces;
    mesh_data.vertices = vertices;
    mesh_data.normals = normals;
    mesh_data.edge_map = edge_map;
    mesh_data.angle_threshold = angle_threshold;

end