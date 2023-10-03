% fprintf('Writing inmypoly_inp_info.txt ...\n')
disp(fpath)
fileID = fopen([fpath 'inmypoly_inp_info.txt'],'w');
fprintf(fileID,'%15.10f %15.10f\n',[dx dy]');
fprintf(fileID,'%5d %5d %5d\n',[itot jtot ktot]');
fprintf(fileID,'%15.10f\n',tol_mypoly);
fprintf(fileID,'%15.10f %15.10f %15.10f\n',Dir_ray_u);
fprintf(fileID,'%15.10f %15.10f %15.10f\n',Dir_ray_v);
fprintf(fileID,'%15.10f %15.10f %15.10f\n',Dir_ray_w);
fprintf(fileID,'%15.10f %15.10f %15.10f\n',Dir_ray_c);
fprintf(fileID,'%8d %8d\n',[size(TR.Points, 1) size(TR.ConnectivityList, 1)]);
fprintf(fileID,'%4d\n',n_threads);
fprintf(fileID,'%d %d\n',[stl_ground diag_neighbs]);
fclose(fileID);

%
fileID = fopen([fpath 'info_matchFacetsToCells.txt'],'w');
fprintf(fileID,'%15.10f %15.10f\n',[dx dy]');
fprintf(fileID,'%5d %5d %5d\n',[itot jtot ktot]');
fprintf(fileID,'%8d %8d\n',[size(TR.ConnectivityList, 1), size(TR.Points, 1)]);
fprintf(fileID,'%d %d %d\n',[periodic_x, periodic_y, diag_neighbs]);
fclose(fileID);

% fprintf('Writing zhgrid.txt ...\n')
fileID = fopen([fpath 'zhgrid.txt'],'w');
fprintf(fileID,'%15.10f\n',zgrid_w');
fclose(fileID);

% fprintf('Writing zfgrid.txt ...\n')
fileID = fopen([fpath 'zfgrid.txt'],'w');
fprintf(fileID,'%15.10f\n',zgrid_c');
fclose(fileID);

% fprintf('Writing vertices.txt ...\n')
fileID = fopen([fpath 'vertices.txt'],'w');
fprintf(fileID,'%15.10f %15.10f %15.10f\n',TR.Points');
fclose(fileID);

% fprintf('Writing faces.txt ...\n')
fileID = fopen([fpath 'faces.txt'],'w');
fprintf(fileID,'%8d %8d %8d %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n',[TR.ConnectivityList TR.incenter TR.faceNormal]');
fclose(fileID);
