function STLtoView3D(infile, outfile, outformat, maxD, row, col)
% Reads an STL (.stl) file and writes a file in View3D (.vs3) format.
% infile: path to STL file.
% outfile: path to write View3D file to.
% outformat: view3d output format. 0: text, 1: binary, 2: sparse.
% row should default to 0, if not then this is the facet for which we
% calculate view factors.

TR = stlread(infile);
F = TR.ConnectivityList;
V = TR.Points;
nV = size(V,1);
nF = size(F,1);
fID = fopen(outfile,'w');

if maxD < Inf
    fprintf(fID,['T\r\nC out=' num2str(outformat) ' maxD=' num2str(maxD) ' row=' num2str(row) ' col=' num2str(col) '\r\nF 3\r\n']);
else
    fprintf(fID,['T\r\nC out=' num2str(outformat) ' row=' num2str(row) ' col=' num2str(col) '\r\nF 3\r\n']);
end

fprintf(fID,'! %4s %6s %6s %6s\r\n','#', 'x', 'y', 'z');
fprintf(fID,'V %4d %6.6f %6.6f %6.6f\r\n',[(1:nV)', V]');
fprintf(fID,'! %4s %6s %6s %6s %6s %6s %6s %6s %6s\r\n','#', 'v1', 'v2', 'v3', 'v4', 'base', 'cmb', 'emit', 'name' );
fprintf(fID,'S %4d %6d %6d %6d %6d %6d %6d %6d %6df\r\n',[(1:nF)', [F zeros(nF,1)], zeros(nF,3), (1:nF)']');
fprintf(fID,'End of Data\r\n');
fclose(fID);
end