function mdfile = latex2markdown(filename,options)
%  Copyright 2020-2022 The MathWorks, Inc.

% What is arguments?
% see: https://jp.mathworks.com/help/matlab/matlab_prog/argument-validation-functions.html
arguments
    filename (1,1) string
    options.outputfilename char = filename
    options.format char {mustBeMember(options.format,{'qiita','github','github_math'})} = 'github'
    options.png2jpeg logical = false
    options.tableMaxWidth (1,1) double = 20
    options.ToC logical = false
end

% Latex filename
[filepath,name,ext] = fileparts(filename);

if filepath == ""
    filepath = pwd;
end
    
if ext == "" % if without extention, add .tex
    latexfile = fullfile(filepath, name + ".tex");
else %
    latexfile = filename;
end

% check if latexfile exists
assert(exist(latexfile,'file')>0, ...
    latexfile + " does not exist. " ...
    + "If you haven't generate latex file from a live script please do so.");

% Import latex file and have it as a string variable
% If the mlx contains double byte charactors, it needs to use fgets.
% Otherwise the following does the work.
% str = string(fileread(latexfile));
fid = fopen(latexfile,'r','n','UTF-8');
str = string;
tmp = fgets(fid);
while tmp > 0
    str = append(str,string(tmp));
    tmp = fgets(fid);
end
fclose(fid);

% Extract body from latex
str = extractBetween(str,"\begin{document}","\end{document}");

% Delete table of contents: 目次は現時点で削除（TODO）
% ex: \label{H_D152BAC0}
str = regexprep(str,"\\matlabtableofcontents{([^{}]+)}", "");
str = regexprep(str,"\\label{[a-zA-Z_0-9]+}","");


%% ToC
% \matlabtitle{タイトル}
% \matlabheading{セクション１}
% \matlabheadingtwo{サブセクション１}
% \matlabheadingthree{サブサブセクション}
% 
% # Table of contents
% - [セクション](#セクション)
%   - [サブセクション](#サブセクション)
%      - [サブサブセクション](#サブサブセクション)

% Delete strings for literal display from generating ToC
% namely..
% \begin{matlabcode}(code)\end{matlabcode}
% \begin{verbatim}(code)\end{verbatim}
% \begin{matlaboutput}(code)\end{matlaboutput}
tmp = regexprep(str,"\\begin{verbatim}.*?\\end{verbatim}","");
tmp = regexprep(tmp,"\\begin{matlabcode}.*?\\end{matlabcode}","");
tmp = regexprep(tmp,"\\begin{matlaboutput}.*?\\end{matlaboutput}","");
toc_str = regexp(tmp,"\\matlabheading(|two|three){([^{}]+)}","match");
toc_id = regexp(tmp,"\\matlabheading(?:|two|three){([^{}]+)}","tokens");

% check if any duplicate id for toc
toc_id = string(toc_id);
if length(toc_id) ~= length(unique(toc_id))
    warning("latex2markdown:ToCdupID","Duplication in section title is found. Some hyperlinks in ToC may not work properly.")
end

% generate ToC with hyperlink for markdown
toc_md = regexprep(toc_str,"\\matlabheading{([^{}]+)}","- [$1](#$1)");
toc_md = regexprep(toc_md,"\\matlabheadingtwo{([^{}]+)}","  - [$1](#$1)");
toc_md = regexprep(toc_md,"\\matlabheadingthree{([^{}]+)}","    - [$1](#$1)");

% IDs need to be all lowercase
% spcaes in IDs need to be replaced with - (the code here is not clean...)
ids = regexp(toc_md,"\(#.*\)","match"); % 
for ii=1:length(ids) % for each IDs
    tmp1 = ids{ii}; 
    tmp2 = replace(tmp1," ","-"); % space is replased by -.
    tmp2 = lower(tmp2); % lowercase
    % replace ID string with a new string
    toc_md = replace(toc_md,tmp1,tmp2);
end
% Note: the code below replaces id itself with -.
% toc_md = regexprep(toc_md,"\(#\S*(\s+)\S*\)","-");
% it would be nice if we can perform any operation on the token directly,
% eg. toc_md = regexprep(toc_str,"\\matlabheading{([^{}]+)}","- [$1](#fun($1))");

toc_md = ["# Table of contents", toc_md]; % add tile
% join the strings
toc_md = join(toc_md,newline);

%% Devide the body into each environment
% Preprocess 1:
% Add 'newline' to the end of the following.
% \end{lstlisting}, \end{verbatim}, \end{matlabcode}, \end{matlaboutput},\end{center}
% \end{matlabtableoutput}, \end{matlabsymbolicoutput}  \vspace{1em}
% 要素ごとに分割しやすいように \end の後に改行が無い場合は１つ追加
str = replace(str,"\end{lstlisting}"+newline,"\end{lstlisting}"+newline+newline);
str = replace(str,"\end{verbatim}"+newline,"\end{verbatim}"+newline+newline);
str = replace(str,"\end{matlabcode}"+newline,"\end{matlabcode}"+newline+newline);
str = replace(str,"\end{matlaboutput}"+newline,"\end{matlaboutput}"+newline+newline);
str = replace(str,"\end{matlabtableoutput}"+newline,"\end{matlabtableoutput}"+newline+newline);
str = replace(str,"\end{matlabsymbolicoutput}"+newline,"\end{matlabsymbolicoutput}"+newline+newline);
str = replace(str,"\end{center}"+newline,"\end{center}"+newline+newline);
str = replace(str,"\vspace{1em}"+newline,"\vspace{1em}"+newline+newline);

% Preprocess 2:
% Replace more than three \n to \n\n.
% 3行以上の空白は2行にしておく
str = regexprep(str,'\n{3,}','\n\n');
% Devide them into parts by '\n\n'
% 空白行で要素に分割
str = strsplit(str,'\n\n')';

% Preprocess 3:
% The following environments will be merge into one string
% for the ease of process.
% MATLABコードが複数行にわたるとうまく処理できないので 対応する \end が見つかるまで連結処理
% \begin{verbatim}
% \begin{lstlisting}
% \begin{matlabcode}
% \begin{matlaboutput}
% \begin{matlabtableoutput}
% \begin{matlabsymbolicoutput}
str = mergeSameEnvironments(str,"lstlisting");
str = mergeSameEnvironments(str,"verbatim");
str = mergeSameEnvironments(str,"matlabcode");
str = mergeSameEnvironments(str,"matlaboutput");
str = mergeSameEnvironments(str,"matlabtableoutput");
str = mergeSameEnvironments(str,"matlabsymbolicoutput");

%% Let's convert latex to markdown
% 1: Process parts that require literal output.
[str, idxLiteral] = processLiteralOutput(str);

% 2: Process that other parts
str2md = str(~idxLiteral);
str2md = processDocumentOutput(str2md,options.tableMaxWidth);

% Equations (数式部分)
str2md = processEquations(str2md, options.format);

% includegraphics (画像部分)
str2md = processincludegraphics(str2md, options.format, options.png2jpeg, name, filepath);

% Apply vertical/horizontal space
% markdown: two spaces for linebreak
% latex: \vspace{1em}
% latex: \hskip1em
str2md = regexprep(str2md,"\\vspace{1em}","  ");
str2md = regexprep(str2md,"\\hskip1em","  ");
str(~idxLiteral) = str2md;

%% Done! Merge them together
strmarkdown = join(str,newline);

% Add toc at the top 
if options.ToC
    strmarkdown = [toc_md; newline; strmarkdown];
end

%% File outputファイル出力
mdfile = options.outputfilename + ".md";
fileID = fopen(mdfile,'w');
fprintf(fileID,'%s\n',strmarkdown);
fclose(fileID);

disp("Coverting latex to markdown is complete");
disp(mdfile);
disp("Note: Related images are saved in " + name + "_images");
