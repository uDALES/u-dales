%% Setting green facets
fpath = ['/media/chris/Project3/uDALES2.0/experiments/' expnr '/']; 
green_facets = dlmread([fpath '/green_facets.' expnr]);
facet_types = ones(nfcts,1);
facet_types(green_facets) = 12;

filename_facets = [fpath 'facets.inp.' expnr];
fileID_facets = fopen(filename_facets,'W');
fprintf(fileID_facets, '# type, normal\n');
fclose(fileID_facets);
dlmwrite(filename_facets, [facet_types, TR.faceNormal], '-append','delimiter',' ','precision','%1.0f');

disp('Done setting types.');