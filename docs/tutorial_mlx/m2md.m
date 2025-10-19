%% m2md.m - Convert MATLAB Live Scripts to Markdown for uDALES Documentation
%
% This script automates the conversion of uDALES tutorial live scripts (.m files) 
% to markdown (.md files) for documentation purposes. It performs two main tasks:
%
% 1. EXPORT: Converts MATLAB live scripts to markdown using MATLAB's export function
% 2. POST-PROCESS: Cleans up the generated markdown to match documentation standards
%
% WORKFLOW:
% - Exports live scripts to markdown files
% - Anonymizes identifying paths for generic documentation  
% - Standardizes code block formatting for consistency with manual docs
% - Fixes markdown formatting issues for MkDocs compatibility
% - Removes unnecessary markup and escaping for cleaner rendering
%
% INPUT FILES (must be in current directory):
% - facets_tutorial.m     : Working with uDALES facet data
% - fields_tutorial.m     : Working with uDALES field data  
% - udbase_tutorial.m     : Introduction to udbase class
% - geometry_tutorial.m   : Working with uDALES geometry
% - utility_tutorial.m    : Tutorial for utility functions
%
% OUTPUT FILES:
% - udales-facets-tutorial.md    : Markdown version (local, anonymized)
% - udales-fields-tutorial.md    : Markdown version (local, anonymized)
% - udales-udbase-tutorial.md    : Markdown version (local, anonymized) 
% - udales-geometry-tutorial.md  : Markdown version (local, anonymized)
% - udales-utility-tutorial.md  : Markdown version (local, anonymized)
%
% POST-PROCESSING OPERATIONS:
% 1. Path anonymization: Replaces specific file paths with generic placeholders
%
% USAGE:
%   Run this script from the /docs/tutorial_mlx/ directory
%
% AUTHORS: uDALES Development Team
% DATE: 2024-2025

%% MATLAB Version Check
% This script requires MATLAB R2025b or later for:
% - export() function with markdown format support for live scripts

fprintf('Checking MATLAB version compatibility...\n');
matlabVersion = version('-release');
matlabYear = str2double(matlabVersion(1:4));
matlabRelease = matlabVersion(5);

% Define minimum requirements: R2025b
minYear = 2025;
minRelease = 'b';

% Check version compatibility
isCompatible = false;
if matlabYear > minYear
    isCompatible = true;
elseif matlabYear == minYear && strcmp(matlabRelease, 'b')
    isCompatible = true;
end

if ~isCompatible
    fprintf('\n❌ ERROR: MATLAB Version Incompatibility\n');
    fprintf('Current version: MATLAB %s\n', matlabVersion);
    fprintf('Required version: MATLAB R2025b or later\n\n');
    fprintf('Please upgrade MATLAB to continue using this script.\n');
    fprintf('The export() function with markdown format was introduced in R2025b.\n');
    error('Incompatible MATLAB version detected. Script terminated.');
end

fprintf('✅ MATLAB %s detected - Version compatible\n\n', matlabVersion);

%% Configuration: Define live scripts to convert and processing parameters
liveScripts = {'facets_tutorial.m', 'fields_tutorial.m', 'udbase_tutorial.m', 'geometry_tutorial.m', 'utility_tutorial.m'};
% Maximum number of references to process (e.g., [1] through [maxrefs])
maxrefs = 10;

% Generate output filenames (with udales- prefix and underscore->hyphen conversion)
outputFileNames = cell(size(liveScripts));
for i = 1:length(liveScripts)
    outputFileNames{i} = ['udales-', strrep(strrep(liveScripts{i}, '_', '-'), '.m', '.md')];
end

%% Step 1: Export live scripts to markdown
fprintf('=== MATLAB LIVE SCRIPT TO MARKDOWN CONVERTER ===\n');
fprintf('Converting %d live scripts to markdown...\n\n', length(liveScripts));

% Loop through each live script and export it to Markdown
for i = 1:length(liveScripts)
    % Check if file exists
    if ~exist(liveScripts{i}, 'file')
        fprintf('Warning: File %s not found, skipping...\n', liveScripts{i});
        continue;
    end
    
    fprintf('Processing: %s\n', liveScripts{i});
    
    try
        % Export the live script to Markdown
        export(liveScripts{i}, outputFileNames{i}, 'Format', 'markdown');
        fprintf('Successfully exported: %s -> %s\n', liveScripts{i}, outputFileNames{i});
    catch ME
        fprintf('Error exporting %s: %s\n', liveScripts{i}, ME.message);
    end
end

%% Step 2: Post-process markdown files for documentation standards
fprintf('\n=== POST-PROCESSING MARKDOWN FILES ===\n');
fprintf('Cleaning up markdown files for documentation consistency...\n\n');

for i = 1:length(liveScripts)
    % Use the pre-generated output filename
    outputFileName = outputFileNames{i};
    
    % Check if the markdown file exists
    if ~exist(outputFileName, 'file')
        fprintf('Warning: Markdown file %s not found, skipping post-processing...\n', outputFileName);
        continue;
    end
    
    fprintf('Post-processing: %s\n', outputFileName);
    
    try
        % Read the file content
        fileID = fopen(outputFileName, 'r', 'n', 'UTF-8');
        if fileID == -1
            error('Could not open file %s', outputFileName);
        end
        content = fread(fileID, '*char')';
        fclose(fileID);
        
        % Perform search and replace operations for documentation standards
        
        % 1. ANONYMIZE PATHS: Replace specific addpath lines with generic template
        %    Converts: addpath('C:\Users\...\u-dales\tools\matlab')
        %    To:       addpath('path_to_udales\tools\matlab')
        content = regexprep(content, 'addpath\(''[^'']*''\)', 'addpath(''path_to_udales\\tools\\matlab'')');
        
        % 2. ANONYMIZE EXPERIMENT DIRECTORIES: Replace specific paths, preserve experiment numbers  
        %    Converts: expdir = 'C:\Users\...\any_folder\065';
        %    To:       expdir = 'path_to_experiments\065';  
        %    Matches any path ending with slash/backslash + number
        content = regexprep(content, 'expdir = ''[^'']*[\\\/](\d+)'';', 'expdir = ''path_to_experiments\\$1'';');
        
        % output the modified content back to the file

        fileID = fopen(outputFileName, 'w', 'n', 'UTF-8');
        if fileID == -1
            error('Could not open file %s for writing', outputFileName);
        end

        fwrite(fileID, content);
        fclose(fileID);

        fprintf('Successfully post-processed: %s\n', outputFileName);
        
    catch ME
        fprintf('Error post-processing %s: %s\n', outputFileName, ME.message);
    end
end

fprintf('\n=== CONVERSION COMPLETE ===\n');
fprintf('All markdown files have been generated and post-processed.\n');
fprintf('Files are ready for documentation deployment.\n\n');
fprintf('NEXT STEPS:\n');
fprintf('1. Review generated .md files for content accuracy\n');
fprintf('2. Copy files to appropriate documentation directories if needed\n');
fprintf('3. Test markdown rendering in your documentation system\n');
fprintf('4. Run the Python post-processor to apply final linting and anchor fixes:\n');
fprintf('   From this directory run: python convert_all_md.py\n');

% Attempt to run the Python post-processor automatically. This requires
% that a suitable Python is available on the PATH. The call is optional and
% will not abort the MATLAB script if it fails; errors are printed for
% diagnostics.
try
    fprintf('\nNote: the Python post-processor will OVERWRITE the .md files in this directory\n');
    fprintf('Running Python post-processor: convert_all_md.py\n');
    [status, cmdout] = system('python convert_all_md.py');
    if status == 0
        fprintf('convert_all_md.py completed successfully.\n');
        if ~isempty(cmdout)
            fprintf('%s\n', cmdout);
        end
    else
        fprintf('convert_all_md.py returned non-zero exit status %d\n', status);
        if ~isempty(cmdout)
            fprintf('Output:\n%s\n', cmdout);
        end
        fprintf('You can re-run the script manually: python convert_all_md.py\n');
    end
catch ME
    fprintf('Failed to run convert_all_md.py from MATLAB: %s\n', ME.message);
    fprintf('Please run it manually from the shell: python convert_all_md.py\n');
end