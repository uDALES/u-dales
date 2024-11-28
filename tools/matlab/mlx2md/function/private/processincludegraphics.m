function str = processincludegraphics(str,format,png2jpeg,filename,filepath)
% Copyright 2020 The MathWorks, Inc.

% Note: There are two cases in the tex
% 1: inserted image: \includegraphics[width=\maxwidth{64.52584044154541em}]{image_0}
% 2: generated figure: \includegraphics[width=\maxwidth{52.78474661314601em}]{figure_0.png}
%
% Inserted images needs to 

% markdown (GitHub): ![string]('path to a image')
% latex では \includegraphics[width=\maxwidth{56.196688409433015em}]{filename}
imageIdx = contains(str,"\includegraphics");
imageParts = str(imageIdx);

% When exported latex from live script, figures and inserted images
% are saved in 'imagedir' as image files.
% latex を生成した時点で Figure 等は画像としてimagedir に保存されている
if isMATLABReleaseOlderThan("R2023b")
    imagedir = filename + "_images/";
else
    imagedir = filename + "_media/";
end
% imagedir = filename + "_images/"; (pre R2023b)
imagedir = strrep(imagedir, '\', '/');

% for each images
for ii=1:length(imageParts)
    fileid = regexp(imageParts(ii),"\\includegraphics\[[^\]]+\]{([^{}]+)}", "tokens");
    [~,fileid_wo_ext,~] = fileparts(fileid{:});
    imagefilename = ls(fullfile(filepath,imagedir,fileid_wo_ext + ".*")); % get the actual filename with extention
    
    % Compress PNG images as JPEG
    if png2jpeg
        [~,imagefilename_wo_ext,ext] = fileparts(imagefilename);
        if strcmp(ext,'.png')
            I = imread(fullfile(imagedir,imagefilename));
            imagefilename = [imagefilename_wo_ext,'_png.jpg'];
            imwrite(I,fullfile(imagedir,imagefilename),'Quality',85);
        end
    end
    
    switch format
        case 'qiita'
            % Qiita に移行する際は、画像ファイルを該当箇所に drag & drop する必要
            % 幅指定する場合には
            % <img src="" alt="attach:cat" title="attach:cat" width=500px>
            imageParts(ii) = regexprep(imageParts(ii),"\\includegraphics\[[^\]]+\]{"+fileid{:}+"}",...
                "<--" + newline ...
                + "**Please drag & drop an image file here**" + newline ...
                + "Filename: **"+imagedir+imagefilename + "**" + newline ...
                + "If you want to set the image size use the following command" + newline ...
                + "<img src="" alt=""attach:cat"" title=""attach:cat"" width=500px>" + newline ...
                + "-->");
            
        case 'github'
            %  ![string]('path to a image')
            imageParts(ii) = regexprep(imageParts(ii),"\\includegraphics\[[^\]]+\]{"+fileid{:}+"}",...
                "!["+imagefilename+"]("+imagedir+imagefilename+")");
    end
end

str(imageIdx) = imageParts;