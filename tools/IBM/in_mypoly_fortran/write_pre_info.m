% fprintf('Writing inmypoly_inp_info.txt ...\n')
fileID = fopen([fpath 'inmypoly_inp_info.txt'],'w');
fprintf(fileID,'%3s\n',expnr);
fprintf(fileID,'%15.10f %15.10f %15.10f\n',[r.dx r.dy r.dz]');
fprintf(fileID,'%5d %5d %5d\n',[r.itot r.jtot r.ktot]');
fprintf(fileID,'%15.10f\n',tol_mypoly);
fprintf(fileID,'%15.10f %15.10f %15.10f\n',Dir_ray_u);
fprintf(fileID,'%15.10f %15.10f %15.10f\n',Dir_ray_v);
fprintf(fileID,'%15.10f %15.10f %15.10f\n',Dir_ray_w);
fprintf(fileID,'%15.10f %15.10f %15.10f\n',Dir_ray_c);
fprintf(fileID,'%5d %5d\n',[length(V) length(F)]);
fclose(fileID);

% fprintf('Writing zhgrid.txt ...\n')
fileID = fopen([fpath 'zhgrid.txt'],'w');
fprintf(fileID,'%15.10f\n',r.zh');
fclose(fileID);

% fprintf('Writing zfgrid.txt ...\n')
fileID = fopen([fpath 'zfgrid.txt'],'w');
fprintf(fileID,'%15.10f\n',r.zf');
fclose(fileID);

% fprintf('Writing vertices.txt ...\n')
fileID = fopen([fpath 'vertices.txt'],'w');
fprintf(fileID,'%15.10f %15.10f %15.10f\n',V');
fclose(fileID);

% fprintf('Writing Stl_data.txt ...\n')
fileID = fopen([fpath 'Stl_data.txt'],'w');
fprintf(fileID,'%5d %5d %5d %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n',[F TR.incenter TR.faceNormal]');  
fclose(fileID);