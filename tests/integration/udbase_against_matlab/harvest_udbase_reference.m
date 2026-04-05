function harvest_udbase_reference(case_dir, output_path)
%HARVEST_UDBASE_REFERENCE Serialize selected udbase outputs for Python tests.
%
% This helper is only used to generate committed MATLAB reference data for the
% Python udbase parity tests. Normal unit tests consume the saved JSON output
% and do not require MATLAB.

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'tools', 'matlab'));

case_dir = char(case_dir);
[~, case_name] = fileparts(case_dir);
expnr = str2double(case_name);

sim = udbase(expnr, case_dir);
facsec_c = sim.facsec.c;
frontal = sim.calculate_frontal_properties();

reference = struct();
reference.case = case_name;
reference.facsec_c = struct( ...
    'facid', facsec_c.facid(:)', ...
    'area', facsec_c.area(:)', ...
    'locs', facsec_c.locs, ...
    'distance', facsec_c.distance(:)');
reference.frontal = struct( ...
    'skylinex', frontal.skylinex, ...
    'skyliney', frontal.skyliney, ...
    'Afx', frontal.Afx, ...
    'Afy', frontal.Afy, ...
    'brx', frontal.brx, ...
    'bry', frontal.bry);

json_text = jsonencode(reference, PrettyPrint=true);
fid = fopen(output_path, 'w');
if fid < 0
    error('Could not open output file: %s', output_path);
end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, json_text, 'char');
