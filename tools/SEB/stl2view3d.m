V = TR.Points;
F = TR.ConnectivityList;

nV = size(V,1);
nF = size(F,1);


%%
fID = fopen([fpath '/facets.vs3'],'w');
fprintf(fID,'T \r\nF  3\r\n');
fprintf(fID,'! %4s %6s %6s %6s\r\n','#', 'x', 'y', 'z');
fprintf(fID,'V %4d %6d %6d %6d\r\n',[(1:nV)', V]');
fprintf(fID,'! %4s %6s %6s %6s %6s %6s %6s %6s %6s\r\n','#', 'v1', 'v2', 'v3', 'v4', 'base', 'cmb', 'emit', 'name' );
fprintf(fID,'S %4d %6d %6d %6d %6d %6d %6d %6d %6d\r\n',[(1:nF)', [F zeros(nF,1)], zeros(nF,3), (1:nF)']');
fprintf(fID,'End of Data\r\n');
fclose(fID);
