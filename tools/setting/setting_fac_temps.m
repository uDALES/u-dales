%% Setting the material types of each facet 
% as well as setting the initial temps

nwalllayers = 3; % Number of layers in the facet
Tinit = 301; % Initial facet temp [k]

Tfacinit = ones(nfcts,1).*Tinit;
Tfacinit_layers = ones(nfcts,nwalllayers).*Tinit;

fname = [fpath, 'Tfacinit.inp.', expnr];
fileID = fopen(fname, 'W');
fprintf(fileID, '# Initial facet tempereatures in radiative equilibrium\n');
fclose(fileID);
dlmwrite(fname, Tfacinit, '-append','delimiter',' ','precision','%4f')


%% Tfacinint_layers is only needed for a warm start and will need to look
% at this more carefully when we do this. 
fname = [fpath, 'Tfacinit_layers.inp.', expnr];
fileID = fopen(fname, 'W');
fprintf(fileID, '# Initial facet tempereatures in radiative equilibrium\n');
fclose(fileID);
dlmwrite(fname, Tfacinit_layers, '-append','delimiter',' ','precision','%4f')
