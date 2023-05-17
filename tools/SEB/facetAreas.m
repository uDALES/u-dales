V = TR.Points;
F = TR.ConnectivityList;

nV = size(V,1);
nF = size(F,1);

for i = 1:nF
    verts = V(F(i,:), :);
    nml = cross(verts(2,:)-verts(1,:),verts(3,:)-verts(1,:));
    A(i) = 0.5 * norm(nml);
end


filename_facetarea = [fpath 'facetarea.inp'];
fileID_facetarea = fopen(filename_facetarea,'W');
fprintf(fileID_facetarea, '# facet area\n');
fclose(fileID_facetarea);
dlmwrite(filename_facetarea, A', '-append','delimiter',' ','precision',8);