clear all
close all

addpath('../matlab/mlx2md/function');

latex2markdown('udbase_tutorial.tex');
latex2markdown('geometry_tutorial.tex');
latex2markdown('fields_tutorial.tex');
latex2markdown('facets_tutorial.tex');

%%
% after running this script, you will have to manually correct the
% hyperlinks in the markdown files, which show as [text][ref]. This can be
% easily done in VS code using the autocomplete functionality. Simply
% select ref and then type # and select the correct section. In VScode you
% can check the display of the markdown file using CTRL-SHIFT-V (in Windows).
