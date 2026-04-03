function harvest_udgeom_reference(stl_path, output_path)
%HARVEST_UDGEOM_REFERENCE Serialize udgeom outputs for offline Python tests.
%
% This helper is only used to generate committed MATLAB reference data for the
% Python udgeom unit tests. Normal test runs consume the saved JSON output and
% do not require MATLAB.

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'tools', 'matlab'));

[geom_dir, geom_name, geom_ext] = fileparts(stl_path);
if isempty(geom_ext)
    geom_file = geom_name;
else
    geom_file = [geom_name geom_ext];
end

geom = udgeom.udgeom(geom_dir);
geom.load(geom_file);

faces = geom.stl.ConnectivityList;
points = geom.stl.Points;

p1 = points(faces(:,1), :);
p2 = points(faces(:,2), :);
p3 = points(faces(:,3), :);

face_centers = (p1 + p2 + p3) ./ 3.0;
face_incenters = incenter(geom.stl);
face_normals = faceNormal(geom.stl);
cross_vec = cross(p2 - p1, p3 - p1, 2);
face_areas = 0.5 .* sqrt(sum(cross_vec.^2, 2));
total_area = sum(face_areas);
signed_volume = sum(dot(p1, cross(p2, p3, 2), 2)) / 6.0;
volume = abs(signed_volume);

% Watertight iff every unique edge appears exactly twice.
edge_list = [
    sort(faces(:, [1 2]), 2);
    sort(faces(:, [2 3]), 2);
    sort(faces(:, [3 1]), 2)
];
[~, ~, edge_ids] = unique(edge_list, 'rows');
edge_counts = accumarray(edge_ids, 1);
is_watertight = all(edge_counts == 2);

bounds = [min(points, [], 1); max(points, [], 1)];
try
    face_to_building_map = geom.get_face_to_building_map();
    buildings = geom.get_buildings();
    outline2d = geom.calculate_outline2d();
    building_outlines = geom.get_building_outlines();
catch
    face_to_building_map = zeros(1, size(faces, 1));
    buildings = {};
    outline2d = {};
    building_outlines = {};
end
outline_edges = geom.get_outline();

building_summaries = cell(length(buildings), 1);
for i = 1:length(buildings)
    b = buildings{i};
    if isstruct(b) && isfield(b, 'triangulation')
        tr = b.triangulation;
        original_face_indices = b.original_face_indices;
    else
        tr = b;
        original_face_indices = [];
    end

    if isempty(tr)
        building_summaries{i} = struct( ...
            'n_faces', 0, ...
            'n_vertices', 0, ...
            'bounds', [], ...
            'original_face_indices', []);
        continue;
    end

    building_summaries{i} = struct( ...
        'n_faces', size(tr.ConnectivityList, 1), ...
        'n_vertices', size(tr.Points, 1), ...
        'bounds', [min(tr.Points, [], 1); max(tr.Points, [], 1)], ...
        'original_face_indices', original_face_indices(:)');
end

outline2d_summaries = cell(length(outline2d), 1);
for i = 1:length(outline2d)
    poly = outline2d{i}.polygon;
    centroid = outline2d{i}.centroid;
    outline2d_summaries{i} = struct( ...
        'polygon', poly, ...
        'centroid', centroid(:)');
end

reference = struct();
reference.n_faces = size(faces, 1);
reference.n_vertices = size(points, 1);
reference.bounds = bounds;
reference.face_centers = face_centers;
reference.face_incenters = face_incenters;
reference.face_normals = face_normals;
reference.face_areas = face_areas;
reference.total_area = total_area;
reference.volume = volume;
reference.is_watertight = is_watertight;
reference.face_to_building_map = face_to_building_map(:)';
reference.buildings = {building_summaries{:}};
reference.outline2d = {outline2d_summaries{:}};
reference.building_outlines = {building_outlines{:}};
reference.outline_edges = outline_edges;

json_text = jsonencode(reference, PrettyPrint=true);
fid = fopen(output_path, 'w');
if fid < 0
    error('Could not open output file: %s', output_path);
end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, json_text, 'char');
