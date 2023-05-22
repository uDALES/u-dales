%% Setting the material types of each facet 
% as well as setting the initial temps


if ~lfacTlyrs
    Tfacinit = ones(nfcts,1).*facT;
    fname = [fpath, 'Tfacinit.inp.', expnr];
    fileID = fopen(fname, 'W');
    fprintf(fileID, '# Initial facet tempereatures in radiative equilibrium\n');
    fclose(fileID);
    dlmwrite(fname, Tfacinit, '-append','delimiter',' ','precision','%4f')
else
    Tfac = ncread(facT_file, 'T');
    Tfacinit_layers = Tfac(:, :, end); 
    fname = [fpath, 'Tfacinit_layers.inp.', expnr];
    fileID = fopen(fname, 'W');
    fprintf(fileID, '# Initial facet tempereatures in radiative equilibrium\n');
    fclose(fileID);
    dlmwrite(fname, Tfacinit_layers, '-append','delimiter',' ','precision','%4f')
end 
