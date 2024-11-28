function str2md = processEquations(str2md, format)
% Copyright 2020 The MathWorks, Inc.
% 
% For Github users: 
% format = 'github_math'
% Use https://latex.codecogs.com
% see. http://idken.net/posts/2017-02-28-math_github/ (Japanese)
%
% format = 'github'
% Use GitHub capability to display equations (first version became
% available in May 2022)
% Leave inline equation as it is (文中の数式は latex で $equation$ なのでそのまま)
% and $$equation$$ will be changed to
% $$
% equation
% $$
%
% For Qiita users: (Qiita platform renders equations via mathML)
% Leave inline equation as it is (文中の数式は latex で $equation$ なのでそのまま)
% and $$equation$$ will be changed to
% ```math
% equation
% ```

switch format
    case 'qiita'
        str2md = regexprep(str2md,"[^`]?\$\$([^$]+)\$\$[^`]?",newline+"```math" + newline + "$1" + newline + "```");
    case 'github'
        str2md = regexprep(str2md,"[^`]?\$\$([^$]+)\$\$[^`]?",newline+"$$" + newline + "$1" + newline + "$$");
    case 'github_math'
        tt = regexp(str2md,"[^`]?\$\$([^$]+)\$\$[^`]?", 'tokens');
        idx = cellfun(@iscell,tt); 
        % if tt contains 0x0 string, horzcat(tt{:}) generates string vector
        % whereas if tt with cell only, horzcat(tt{:}) generates cell
        % vector... so.
        parts = horzcat(tt{idx});
        for ii=1:length(parts)
            eqncode = replace(parts{ii},string(newline)," ");
            eqncode = replace(eqncode," ", "&space;");
            partsMD = "<img src=""https://latex.codecogs.com/gif.latex?" ...
                + eqncode + """/>";
            str2md = replace(str2md, "$$"+parts(ii)+"$$", partsMD);
        end
        
        % Inline
        tt = regexp(str2md,"[^`$]\$([^$]+)\$[^`$]", 'tokens');
        parts = horzcat(tt{:});
        for ii=1:length(parts)
            eqncode = replace(parts(ii)," ", "&space;");
            partsMD = "<img src=""https://latex.codecogs.com/gif.latex?\inline&space;" ...
                + eqncode + """/>";
            str2md = replace(str2md, "$"+parts(ii)+"$", partsMD);
        end
end
