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
% 2. Code block standardization: Converts matlabTextOutput -> text blocks
% 3. LaTeX cleanup: Removes trailing spaces that break KaTeX compilation
% 4. 4-backtick fix: Converts ````markdown to ```markdown and removes document wrapper backticks
% 5. Backslash cleanup: Removes unnecessary escaping in regular text (post\-processing -> post-processing, function\_name -> function_name)
% 6. Anchor management: Adds readable anchor IDs to headings and converts table of contents links to match
% 7. Reference formatting: Fixes LaTeX references (\[1\] -> [1])
%
% USAGE:
%   Run this script from the /docs/tutorial_mlx/ directory
%
% AUTHORS: uDALES Development Team
% DATE: 2024-2025

%% MATLAB Version Check
% This script requires MATLAB R2025b or later for:
% - export() function with markdown format support for live scripts
% - Enhanced containers.Map functionality
% - Robust string and regexp handling for post-processing

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

%% Helper function to escape special regex characters
function escaped = regexpescape(str)
    % Escape special regex characters for literal matching
    escaped = regexprep(str, '([\[\](){}.*+?^$|\\])', '\\$1');
end

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
        
        % 3. STANDARDIZE CODE BLOCKS: Convert MATLAB-specific output blocks to generic text
        %    Converts: ```matlabTextOutput  ->  ```text
        %    This ensures consistency with manually maintained documentation
        content = regexprep(content, '```matlabTextOutput', '```text');
        
        % 4. FIX LATEX FORMATTING: Remove trailing spaces after $$ that break KaTeX compilation
        %    Converts: $$<space><space>  ->  $$
        content = regexprep(content, '\$\$\s+$', '$$', 'lineanchors');
        
        % 5. FIX 4-BACKTICK CODE BLOCKS: Convert 4 backticks to 3 for MkDocs compatibility
        %    Converts: ````markdown  ->  ```markdown
        %    Converts: ````text     ->  ```text
        %    Also removes document wrapper backticks at start/end
        content = regexprep(content, '````(\w*)', '```$1');
        % Remove wrapping backticks at document start and end
        content = regexprep(content, '^```markdown\s*\n', '', 'once');
        content = regexprep(content, '\n```\s*$', '', 'once');
        
        % 6. REMOVE UNNECESSARY BACKSLASH ESCAPING: Clean up escaped hyphens and underscores
        %    Converts: post\-processing  ->  post-processing
        %    Converts: time\-averaged    ->  time-averaged
        %    Converts: function\_name    ->  function_name
        content = regexprep(content, '(\w)\\-(\w)', '$1-$2');
        content = regexprep(content, '(\w)\\- ', '$1- ');
        content = regexprep(content, '(\w)\\_(\w)', '$1_$2');
        
        % 7. ADD LIST INDENTATION TO MATCH MANUAL DOCS FORMAT
        %    Converts: -<space> -> <space><space><space>-<space> (3-space indentation)
        fprintf('  - Adding list indentation for manual docs compatibility...\n');
        content = regexprep(content, '^- ', '   - ', 'lineanchors');
        
        % Save original markdown for debugging (before anchor processing)
        % 8. MULTI-PASS ANCHOR AND HEADING PROCESSING SYSTEM
        %    This system extracts anchors, maps them to headings, and creates descriptive anchors
        
        % PASS 1: Extract ToC function references and anchor tag mappings
        
        lines = strsplit(content, '\n');
        tocFunctionRefs = {};  % Store {functionName, anchorId, linkText, lineNum}
        anchorTagMap = containers.Map();  % Map anchor tag ID -> line number of tag
        
        for lineNum = 1:length(lines)
            line = lines{lineNum};
            
            % Find ToC function links: [**function_name**](#H_xxxx)
            funcLinkMatches = regexp(line, '\[\*\*([a-zA-Z][a-zA-Z0-9_]*)\*\*\]\(#([^)]+)\)', 'tokens');
            for matchIdx = 1:length(funcLinkMatches)
                functionName = funcLinkMatches{matchIdx}{1};
                anchorId = funcLinkMatches{matchIdx}{2};
                linkText = funcLinkMatches{matchIdx}{1}; % Use function name as link text
                tocFunctionRefs{end+1} = {functionName, anchorId, linkText, lineNum};
            end
            
            % Find anchor tags: <a id="H_xxxx"></a>
            anchorTagMatches = regexp(line, '<a id="([^"]*)">', 'tokens');
            for matchIdx = 1:length(anchorTagMatches)
                tagId = anchorTagMatches{matchIdx}{1};
                anchorTagMap(tagId) = lineNum;
            end
        end
        
        % PASS 2: Extract all headings and their line numbers
        
        headings = {};  % Store {headingText, lineNum, headingLevel, functionName}
        
        for lineNum = 1:length(lines)
            line = lines{lineNum};
            % Match headings: # heading text
            headingMatch = regexp(line, '^(#+)\s*(.+)$', 'tokens', 'once');
            if ~isempty(headingMatch)
                headingLevel = length(headingMatch{1});
                headingText = strtrim(headingMatch{2});
                % Remove any existing anchor tags from heading text
                headingText = regexprep(headingText, '\s*\{#[^}]+\}$', '');
                
                % Extract function name from heading (first word before colon)
                funcMatch = regexp(headingText, '^([a-zA-Z][a-zA-Z0-9_]*)', 'tokens', 'once');
                functionName = '';
                if ~isempty(funcMatch)
                    functionName = funcMatch{1};
                end
                
                headings{end+1} = {headingText, lineNum, headingLevel, functionName};
            end
        end
        
        % PASS 3: Create function-to-descriptive-anchor mapping
        
        functionToAnchor = containers.Map();  % Dictionary: functionName -> descriptiveAnchor
        
        % For each ToC function reference, find the corresponding heading
        for tocIdx = 1:length(tocFunctionRefs)
            functionName = tocFunctionRefs{tocIdx}{1};
            anchorId = tocFunctionRefs{tocIdx}{2};
            
            % Find the anchor tag line number
            if isKey(anchorTagMap, anchorId)
                anchorTagLineNum = anchorTagMap(anchorId);
                
                % Find the next heading after this anchor tag
                bestHeading = [];
                minDistance = inf;
                
                for headingIdx = 1:length(headings)
                    headingLineNum = headings{headingIdx}{2};
                    headingFuncName = headings{headingIdx}{4};
                    
                    % Must be after the anchor tag and match function name
                    if headingLineNum > anchorTagLineNum && strcmp(headingFuncName, functionName)
                        distance = headingLineNum - anchorTagLineNum;
                        if distance < minDistance
                            minDistance = distance;
                            bestHeading = headings{headingIdx};
                        end
                    end
                end
                
                % If we found the matching heading, create descriptive anchor
                if ~isempty(bestHeading)
                    headingText = bestHeading{1};
                    
                    % Create descriptive anchor (like manual docs format)
                    % Convert "method_name: description text" -> "method_name-description-text"
                    descriptiveAnchor = lower(headingText);
                    descriptiveAnchor = regexprep(descriptiveAnchor, '[ :,()]+', '-');
                    descriptiveAnchor = regexprep(descriptiveAnchor, '[^a-z0-9\-_]', '');
                    descriptiveAnchor = regexprep(descriptiveAnchor, '-+', '-');
                    descriptiveAnchor = regexprep(descriptiveAnchor, '^-+|-+$', '');
                    
                    functionToAnchor(anchorId) = descriptiveAnchor;
                end
            end
        end
        
        % PASS 4: Update all anchor references to use descriptive anchors
        
        updatedLines = lines;  % Copy of lines to modify
        updatedCount = 0;
        
        for lineNum = 1:length(updatedLines)
            line = updatedLines{lineNum};
            originalLine = line;
            
            % Find and replace all anchor references in this line
            linkMatches = regexp(line, '\[([^\]]+)\]\(#([^)]+)\)', 'tokens');
            for matchIdx = length(linkMatches):-1:1  % Process in reverse to maintain positions
                linkText = linkMatches{matchIdx}{1};
                anchorId = linkMatches{matchIdx}{2};
                
                if isKey(functionToAnchor, anchorId)
                    descriptiveAnchor = functionToAnchor(anchorId);
                    
                    % Replace the anchor reference
                    oldPattern = sprintf('\\[%s\\]\\(#%s\\)', regexpescape(linkText), regexpescape(anchorId));
                    newPattern = sprintf('[%s](#%s)', linkText, descriptiveAnchor);
                    line = regexprep(line, oldPattern, newPattern);
                    updatedCount = updatedCount + 1;
                end
            end
            
            if ~strcmp(line, originalLine)
                updatedLines{lineNum} = line;
            end
        end
        
        if updatedCount > 0
            fprintf('    Updated %d anchor references to descriptive format\n', updatedCount);
        end
        
        % Reconstruct content from updated lines
        content = strjoin(updatedLines, '\n');
        
        % 8. FIX LATEX REFERENCE FORMATTING: Clean up escaped citation references only
        %    Converts: \[1\]  ->  [1] (in citations, not units)
        %    Converts: \[2\]  ->  [2] (in citations, not units)
        %    Ultra-simple approach: literal string replacement for citations up to maxrefs
        for refnum = 1:maxrefs
            oldref = sprintf('\\[%d\\]', refnum);
            newref = sprintf('[%d]', refnum);
            content = strrep(content, oldref, newref);
        end
        
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