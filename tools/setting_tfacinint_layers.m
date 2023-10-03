%% Making Tfacinnit_layers files.

expnrbase = '209';
expnrnew = '995';

facT_filepath  = ['/media/chris/Project4/outputs/' expnrbase '/facT.' expnrbase '.nc'];
Tfac = ncread(facT_filepath,'T');
Tfacinit_layers = Tfac(:, :, end);
fname = ['/media/chris/Project3/uDALES2.0/experiments/' expnrnew '/Tfacinit_layers.inp.' expnrnew];
fileID = fopen(fname, 'W');
fprintf(fileID, '# Initial facet tempereatures in radiative equilibrium\n');
fclose(fileID);
dlmwrite(fname, Tfacinit_layers, '-append','delimiter',' ','precision','%4f');