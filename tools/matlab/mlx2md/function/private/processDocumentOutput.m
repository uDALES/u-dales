function str2md = processDocumentOutput(str2md,tableMaxWidth)
% Copyright 2020 The MathWorks, Inc.
    
%% 2-1: Fix latex conventions for non-literal parts
% ^ (live script) -> \textasciicircum{} (latex)
str2md = replace(str2md,"\textasciicircum{}","^");
% _ (live script) -> \_ (latex) example: test_case -> test\_ case
str2md = replace(str2md,"\_","_");
% / backslash (live script) -> \textbackslash{} (latex)
str2md = replace(str2md,"\textbackslash{}","\");
% > (live script) -> \textgreater{} (latex)
str2md = replace(str2md,"\textgreater{}",">");
% < (live script) -> \textless{} (latex)
str2md = replace(str2md,"\textless{}","<");
% $ (live script) -> \$ (latex)
str2md = replace(str2md,"\$","$");
% % (live script) -> \% (latex)
str2md = replace(str2md,"\%","%");

% These will be left as they are till the end of this function
% since these affect the markdown format 
% { (live script) -> \} (latex) (leave it till end)
% } (live script) -> \{ (latex) (leave it till end)

% To deal with \{ and \} inside other commands, it's easier to
% make regular expression if we change these to letters. (will change it back later)
str2md = replace(str2md,"\{", "BackslashCurlyBlacketOpen");
str2md = replace(str2md,"\}", "BackslashCurlyBlacketClose");

%% 2-2: Text decoration (文字装飾部分)
% \textbf{bold}
% \textit{italic}
% \underline{underline}
% \texttt{equispace}
% and all the possible conbinations of these four.

% Need to keep this execution sequence
str2md = regexprep(str2md,"\\textbf{([^{}]+)}","**$1**");
str2md = regexprep(str2md,"\\textit{([^{}]+)}","*$1*");
str2md = regexprep(str2md,"\\underline{([^{}]+)}","$1"); % Ignore underline (下線は無視）
str2md = regexprep(str2md,"\\texttt{(\*{0,3})([^*{}]+)(\*{0,3})}","$1`$2`$3");

% Note on the processing \texttt
% Example: 
% str = "\\texttt{\\textbf{EquispaceBold}}";
% str = regexprep(str,"\\textbf{([^{}]+)}","**$1**");
% str = regexprep(str,"\\texttt{([^{}]+)}","`$1`");
% gives
% `**EquispaceBold**`
% which does not work. ` ` needs to be most inside.
% `` が最も外側にくるが一番内側にある必要がある。

%% 2-3: Hyperlinks (ハイパーリンク)
% Markdown: [string](http://xxx.com)
% latex: \href{http://xxx.com}{string}
str2md = regexprep(str2md,"\\href{([^{}]+)}{([^{}]+)}","[$2]($1)");
str2md = regexprep(str2md,"\\hyperref\[(.*?)\]\{(.*?)\}","[$2]($1)");

%% 2-4: Titile and headings (見出し部分)
str2md = regexprep(str2md,"\\matlabtitle{([^{}]+)}","# $1");
str2md = regexprep(str2md,"\\matlabheading{([^{}]+)}","# $1");
str2md = regexprep(str2md,"\\matlabheadingtwo{([^{}]+)}","## $1");
str2md = regexprep(str2md,"\\matlabheadingthree{([^{}]+)}","### $1");

% Put \{ and \{ back.
str2md = replace(str2md,"BackslashCurlyBlacketOpen","\{");
str2md = replace(str2md, "BackslashCurlyBlacketClose","\}");

%% 2-5: Quotation (引用パラグラフ)
% Markdown: >
% Latex:
% \begin{par}
% \begin{center}
% xxxx
% \end{center}
% \end{par}
% Note: \includegraphics is an exception
idxNonGraphics = ~contains(str2md,"\includegraphics");
str2md(idxNonGraphics) = replace(str2md(idxNonGraphics),...
    "\begin{par}"+newline+"\begin{center}"+newline,"> ");

%% 2-6: Delete unnecessary commands (不要コマンドを削除)
% Commands to specify the text position
str2md = erase(str2md,"\begin{par}");
str2md = erase(str2md,"\end{par}");
str2md = erase(str2md,"\begin{flushleft}");
str2md = erase(str2md,"\end{flushleft}");
str2md = erase(str2md,"\begin{flushright}");
str2md = erase(str2md,"\end{flushright}");
str2md = erase(str2md,"\begin{center}");
str2md = erase(str2md,"\end{center}");

%% 2-7: Unordered list (リスト)
% markdown: add - to each item
% latex:
%      \begin{itemize}
%      \setlength{\itemsep}{-1ex}
%         \item{\begin{flushleft} リスト１ \end{flushleft}}
%         \item{\begin{flushleft} リスト２ \end{flushleft}}
%         \item{\begin{flushleft} リスト３ \end{flushleft}}
%      \end{itemize}
str2md = erase(str2md,"\setlength{\itemsep}{-1ex}"+newline);
itemizeIdx = contains(str2md,["\begin{itemize}","\end{itemize}"]);
itemsParts = str2md(itemizeIdx);
partsMarkdown = regexprep(itemsParts,"\\item{([^{}]+)}","- $1");
partsMarkdown = erase(partsMarkdown,["\begin{itemize}","\end{itemize}"]);
str2md(itemizeIdx) = partsMarkdown;

%% 2-8: Ordered list (数付きリスト)
% markdown: 1. itemname
% latex:
%      \begin{enumerate}
%      \setlength{\itemsep}{-1ex}
%         \item{\begin{flushleft} リスト１ \end{flushleft}}
%         \item{\begin{flushleft} リスト２ \end{flushleft}}
%         \item{\begin{flushleft} リスト３ \end{flushleft}}
%      \end{enumerate}
str2md = erase(str2md,"\setlength{\itemsep}{-1ex}"+newline);
itemizeIdx = contains(str2md,["\begin{enumerate}","\end{enumerate}"]);
itemsParts = str2md(itemizeIdx);
partsMarkdown = regexprep(itemsParts,"\\item{([^{}]+)}","1. $1");% Any numder works
partsMarkdown = erase(partsMarkdown,["\begin{enumerate}","\end{enumerate}"]);
str2md(itemizeIdx) = partsMarkdown;

%% 2-9: Symbolic output (シンボリック出力)
% markdown: inline equation
% latex:
% \begin{matlabsymbolicoutput}
% ans =
%     $\displaystyle -\cos \left(x\right)$
% \end{matlabsymbolicoutput}
%
% and
%
% \begin{matlabsymbolicoutput}
% a = 
%     $\displaystyle \left(\begin{array}{cccc}
% \cos \left(\theta \right) & -\sin \left(\theta \right) & 0 & 0\\
% \sin \left(\theta \right) & \cos \left(\theta \right) & 0 & 0\\
% 0 & 0 & 1 & 0\\
% 0 & 0 & 0 & 1
% \end{array}\right)$
% \end{matlabsymbolicoutput}

symoutIdx = contains(str2md,["\begin{matlabsymbolicoutput}","\end{matlabsymbolicoutput}"]);
symoutParts = str2md(symoutIdx);
tmp = erase(symoutParts,"\begin{matlabsymbolicoutput}"+newline);
tmp = replace(tmp,"$\displaystyle","$$");
partsMarkdown = replace(tmp,"$"+newline+"\end{matlabsymbolicoutput}","$$");
str2md(symoutIdx) = partsMarkdown;
% NOTE: This part will be processed by processEquations.m

%% 2-10: table output (table 型データの出力)
% markdown:
% | TH left | TH center | TH right |
% | :--- | :---: | ---: |
% | TD | TD | TD |
% | TD | TD | TD |
% latex:
% \begin{matlabtableoutput}
% {
% \begin{tabular} {|l|c|r|}\hline
% \mlcell{TD} & \mlcell{TD} & \mlcell{TD} \\ \hline
% \mlcell{TD} & \mlcell{TD} & \mlcell{TD} \\
% \hline
% \end{tabular}
% }
% \end{matlabtableoutput}
idxTblOutput = startsWith(str2md,"\begin{matlabtableoutput}"+newline);
tableLatex = extractBetween(str2md(idxTblOutput),...
    "\begin{tabular}","\end{tabular}");

tableMD = strings(sum(idxTblOutput),1);
for ii=1:sum(idxTblOutput)
    tablecontents = split(tableLatex(ii),"\hline");
    formatLatex = tablecontents(1); % {|l|c|r|}
    headerLatex = tablecontents(2); % \mlcell{TD} & \mlcell{TD} & \mlcell{TD} \\ \hline
    bodyLatex = tablecontents(3:end); % and the rest.
    
    format = regexp(formatLatex,"\{([^{}]+)}",'tokens');
    format = format{:};
    format = replace(format, "c",":--:");
    format = replace(format, "l",":--");
    format = replace(format, "r","--:");
    
    % MultiColumn is not standard in markdown.
    % It only happens as a variable name in MATLAB
    % so adding special case for processing headerLatex
    multicol = regexp(headerLatex,"\\multicolumn{(\d)+}",'tokens');
    tmp = regexp(headerLatex,"\\mlcell{(.|\s)*?} (?:&|\\\\)",'tokens');
    if isempty(multicol)
        header = "|" + join([tmp{:}],"|") + "|";
    else
        nRepeat = double(multicol{:});
        header = "|" + join([tmp{:}],"|") + join(repmat("|",1,nRepeat));
    end
    
    body = string;
    for jj=1:length(bodyLatex)
        tmp = regexp(bodyLatex(jj),"\\mlcell{(.|\s)*?} (?:&|\\\\)",'tokens');
        if isempty(tmp)
            break;
        end
        tmp = cellfun(@(str1) cutStringLength(str1, tableMaxWidth), tmp, 'UniformOutput', false);
        
        % Adding escape to text that affects markdown table (\n and |)
        tmp = cellfun(@(str1) replace(str1,"|","\|"), tmp, 'UniformOutput', false);
        tmp = cellfun(@(str1) replace(str1,newline,"<br>"), tmp, 'UniformOutput', false);
        body = body + "|" + join([tmp{:}],"|") + "|" + newline;
    end
    
    tableMD(ii) = strjoin([header,format,body],newline);
end
str2md(idxTblOutput) = tableMD;

%% 2-11: Hyperlinks 2 (ハイパーリンク)
% Markdown: [string](http://xxx.com)
% latex: \href{http://xxx.com}{string}
str2md = regexprep(str2md,"\\href{([^{}]+)}{([^{}]+)}","[$2]($1)");
% str2md = regexprep(str2md,"\\hyperref\[(.*?)\]\{(.*?)\}","[$2]($1)");

%% finish up
str2md = replace(str2md,"\{","{");
str2md = replace(str2md,"\}","}");

end

function str2 = cutStringLength(str1, N)
    
    str2 = str1;
    if strlength(str1) > N
       tmp = char(str1);
       str2 = string(tmp(1:N)) + "...";
    end
end