%% Setting the initial profile and the x and z grids
z    = zf;  % Real space locations of computational z-grid [m]
thl0 = 301; % Liquid potential temperature [k]
qt0  = 0.01; % Humidity [-]
v0   = 0;   % [m/s]
u0   = 1;   % [m/s]
tke0 = 0;   % Turbulence kinetic energy[m^2/s^2]

pr = zeros(length(zf),6);
pr(:,1) = zf;
pr(:,2) = thl0;
pr(:,3) = qt0;
pr(:,4) = u0;
pr(:,5) = v0;
pr(:,6) = tke0;

prof = fopen([fpath 'prof.inp.' expnr], 'w');
fprintf(prof, '%-12s\n', '# SDBL flow');
fprintf(prof, '%-60s\n', '# z thl qt u v tke');
fprintf(prof, '%-20.15f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f\n', pr');
fclose(prof);
 
xgrid = fopen([fpath 'xgrid.inp.' expnr], 'w');
fprintf(xgrid, '%12s\n', '#     x-grid');
fprintf(xgrid, '%12s\n', '#           ');
fprintf(xgrid, '%-20.15f\n', xf);
fclose(xgrid);

zgrid = fopen([fpath 'zgrid.inp.' expnr], 'w');
fprintf(zgrid, '%12s\n', '#     z-grid');
fprintf(zgrid, '%12s\n', '#           ');
fprintf(zgrid, '%-20.15f\n', zf);
fclose(zgrid);