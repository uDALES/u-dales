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
% - Fixes LaTeX formatting issues for proper rendering
%
% INPUT FILES (must be in current directory):
% - facets_tutorial.m     : Working with uDALES facet data
% - fields_tutorial.m     : Working with uDALES field data  
% - udbase_tutorial.m     : Introduction to udbase class
% - geometry_tutorial.m   : Working with uDALES geometry
%
% OUTPUT FILES:
% - udales-facets-tutorial.md    : Markdown version (local, anonymized)
% - udales-fields-tutorial.md    : Markdown version (local, anonymized)
% - udales-udbase-tutorial.md    : Markdown version (local, anonymized) 
% - udales-geometry-tutorial.md  : Markdown version (local, anonymized)
%
% POST-PROCESSING OPERATIONS:
% 1. Path anonymization: Replaces specific file paths with generic placeholders
% 2. Code block standardization: Converts matlabTextOutput -> text blocks
% 3. LaTeX cleanup: Removes trailing spaces that break KaTeX compilation
%
% USAGE:
%   Run this script from the /docs/tutorial_mlx/ directory
%
% AUTHORS: uDALES Development Team
% DATE: 2024-2025

%% Configuration: Define live scripts to convert
liveScripts = {'facets_tutorial.m', 'fields_tutorial.m', 'udbase_tutorial.m', 'geometry_tutorial.m'};

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
        content = regexprep(content, "addpath\('[^']*'\)", "addpath('path_to_udales\\tools\\matlab')");
        
        % 2. ANONYMIZE EXPERIMENT DIRECTORIES: Replace specific paths, preserve experiment numbers  
        %    Converts: expdir = 'C:\Users\...\experiments\065';
        %    To:       %expdir = 'path_to_experiments\065';
        content = regexprep(content, "expdir = '[^']*[\\/](\\d+)';", "%expdir = 'path_to_experiments\\$1';");
        
        % 3. STANDARDIZE CODE BLOCKS: Convert MATLAB-specific output blocks to generic text
        %    Converts: ```matlabTextOutput  ->  ```text
        %    This ensures consistency with manually maintained documentation
        content = regexprep(content, '```matlabTextOutput', '```text');
        
        % 4. FIX LATEX FORMATTING: Remove trailing spaces after $$ that break KaTeX compilation
        %    Converts: $$<space><space>  ->  $$
        content = regexprep(content, '\$\$\s+$', '$$', 'lineanchors');
        
        % Write the modified content back to file
        fileID = fopen(outputFileName, 'w', 'n', 'UTF-8');
        if fileID == -1
            error('Could not write to file %s', outputFileName);
        end
        fwrite(fileID, content, 'char');
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