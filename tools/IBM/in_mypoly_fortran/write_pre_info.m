% fprintf('Writing inmypoly_inp_info.txt ...\n')
fileID = fopen([fpath 'inmypoly_inp_info.txt'],'w');
fprintf(fileID,'%15.10f %15.10f %15.10f\n',[dx dy dz]');
fprintf(fileID,'%5d %5d %5d\n',[itot jtot ktot]');
fprintf(fileID,'%15.10f\n',tol_mypoly);
fprintf(fileID,'%15.10f %15.10f %15.10f\n',Dir_ray_u);
fprintf(fileID,'%15.10f %15.10f %15.10f\n',Dir_ray_v);
fprintf(fileID,'%15.10f %15.10f %15.10f\n',Dir_ray_w);
fprintf(fileID,'%15.10f %15.10f %15.10f\n',Dir_ray_c);
fprintf(fileID,'%8d %8d\n',[size(TR.Points, 1) size(TR.ConnectivityList, 1)]);
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

% fprintf('Writing Stl_data.txt ...\n')
fileID = fopen([fpath 'Stl_data.txt'],'w');
fprintf(fileID,'%8d %8d %8d %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n',[TR.ConnectivityList TR.incenter TR.faceNormal]');  
fclose(fileID);