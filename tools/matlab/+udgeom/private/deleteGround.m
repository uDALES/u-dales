function [filtered_TR, face_mapping] = deleteGround(TR)
% DELETEGROUND Remove ground-level faces from a triangulated mesh
%
% SYNTAX:
%   filtered_TR = deleteGround(TR)
%
% DESCRIPTION:
%   Removes triangular faces that lie entirely at ground level (Z=0) from
%   a triangulated mesh. This is commonly used in urban flow simulations
%   to eliminate ground plane elements and focus analysis on building
%   structures and elevated geometry.
%
%   The function performs intelligent mesh cleaning by:
%   1. Identifying faces where ALL vertices have Z=0 coordinates
%   2. Removing these ground-level faces from the mesh
%   3. Compacting the vertex list to remove unused vertices
%   4. Reindexing face connectivity to maintain mesh integrity
%   5. Providing statistics on the filtering results
%
% INPUT ARGUMENTS:
%   TR - Triangulation object containing the original mesh with:
%        * ConnectivityList: [N×3] triangle connectivity matrix
%        * Points: [M×3] vertex coordinates [x, y, z]
%
% OUTPUT ARGUMENTS:
%   filtered_TR - New triangulation object with ground faces removed:
%                 * Reduced connectivity list (ground faces eliminated)
%                 * Compacted vertex list (unused vertices removed)
%                 * Proper reindexing to maintain mesh integrity
%   face_mapping - Array mapping filtered face indices to original face indices:
%                  * Length equals number of faces in filtered_TR
%                  * face_mapping(i) gives original face index for filtered face i
%                  * Enables reconstruction of full geometry mapping
%
% ALGORITHM:
%   1. Extract mesh connectivity and vertex coordinates
%   2. Identify ground faces (all vertices at Z=0)
%   3. Filter out ground faces from connectivity matrix
%   4. Determine which vertices are still referenced
%   5. Create compact vertex array with only used vertices
%   6. Remap face indices to new compact vertex numbering
%   7. Construct new triangulation with filtered data
%
% MESH FILTERING CRITERIA:
%   - Ground faces: Triangles where ALL three vertices have Z = 0
%   - Preserved faces: Any triangle with at least one vertex above Z = 0
%   - Vertex compaction: Removes vertices not referenced by any face
%
% OPERATION:
%   Function executes silently and returns the filtered triangulation
%   without generating any plots or console output.
%
% EXAMPLE:
%   % Load a mesh containing buildings and ground
%   TR_original = stlread('urban_scene.stl');
%   
%   % Remove ground plane for building analysis
%   [TR_buildings, face_mapping] = deleteGround(TR_original);
%   
%   % Result: Clean mesh with mapping to original face indices
%   original_face_1 = face_mapping(1);  % Original index of first filtered face
%
% APPLICATIONS:
%   - Urban flow simulation preprocessing
%   - Building geometry isolation from terrain
%   - CAD model cleaning for CFD analysis
%   - Mesh preparation for facet-based simulations
%
% NOTES:
%   - Function assumes ground level is exactly at Z=0
%   - Preserves mesh topology and connectivity integrity
%   - Automatically handles vertex reindexing
%   - Creates visualization for quality verification
%
% SEE ALSO:
%   triangulation, stlread, trisurf, udgeom.createFlatSurface
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

    %% EXTRACT MESH DATA
    % Get triangular connectivity and vertex coordinates from input triangulation
    faces = TR.ConnectivityList;  % [N×3] Triangle connectivity matrix
    points = TR.Points;           % [M×3] Vertex coordinates [x, y, z]
    z_coords = points(:, 3);      % Extract Z-coordinates for ground detection
    
    %% IDENTIFY GROUND FACES
    % Find faces where ALL vertices lie exactly at ground level (Z=0)
    % This ensures only true ground plane triangles are removed
    faces_z0 = all(z_coords(faces) == 0, 2);
    
    % Create logical mask for faces to preserve (non-ground faces)
    faces_to_keep = faces(~faces_z0, :);
    
    %% COMPACT VERTEX LIST
    % Determine which vertices are still referenced by remaining faces
    % This step is crucial for memory efficiency and proper mesh structure
    used_vertex_indices = unique(faces_to_keep(:));
    
    % Extract only the vertices that are actually used in the filtered mesh
    new_points = points(used_vertex_indices, :);
    
    %% REINDEX FACE CONNECTIVITY
    % Map old vertex indices to new compact indices to maintain mesh integrity
    % This is necessary because vertex removal changes the numbering scheme
    [~, new_indices] = ismember(faces_to_keep, used_vertex_indices);
    new_faces = reshape(new_indices, size(faces_to_keep));
    
    %% CREATE FILTERED TRIANGULATION
    % Construct new triangulation object with ground faces removed
    filtered_TR = triangulation(new_faces, new_points);
    
    %% CREATE FACE MAPPING
    % Create mapping from filtered face indices to original face indices
    original_face_indices = find(~faces_z0);  % Indices of non-ground faces in original mesh
    face_mapping = original_face_indices;     % Direct mapping array
    
    %% FILTERING COMPLETE
    % Ground filtering completed with mapping information

end