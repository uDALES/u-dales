%% Setting green facets

facet_types = ones(nfcts,1);

filename_facets = [fpath 'facets.inp.' expnr];
fileID_facets = fopen(filename_facets,'W');
fprintf(fileID_facets, '# type, normal\n');
fclose(fileID_facets);
dlmwrite(filename_facets, [facet_types, TR.faceNormal], '-append','delimiter',' ','precision','%1.0f');

albedos = [];
typids = factypes(:,1);
for i = 1:length(facet_types)
    my_typid = facet_types(i);
    alb = factypes(find(typids==my_typid),5);
    albedos = [albedos;alb];
end


