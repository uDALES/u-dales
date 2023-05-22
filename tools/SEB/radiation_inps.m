%% Setting constants


az = solaz - xaz;
nsun = [sind(zen)*cosd(az), -sind(zen)*sind(az), cosd(zen)];
F = TR.ConnectivityList;
V = TR.Points;


%% Assumes view factors have been calculated using View3D

vf = dlmread([fpath 'vf.txt'], ' ', 2, 0);

vf(end,:) = [];

%%

if ~lvfsparse
    ncid = netcdf.create([fpath 'vf.nc.inp.' expnr], 'NC_WRITE');
    dimidrow = netcdf.defDim(ncid,'rows', nfcts);
    dimidcol = netcdf.defDim(ncid,'columns', nfcts);
    varid = netcdf.defVar(ncid,'view factor','NC_FLOAT',[dimidrow dimidcol]);
    netcdf.endDef(ncid);
    netcdf.putVar(ncid,varid,vf);
    netcdf.close(ncid);
else 
    vftol = zeros(size(vf));
    vftol(vf > min_vf) = vf(vf > min_vf);
    vfsparse = sparse(double(vftol));
    [i,j,s] = find(vfsparse);
    fID = fopen([fpath 'vfsparse.inp.' expnr], 'w');
    fprintf(fID, '%d %d %.6f \n', [i, j, s]');
end 


%%

svf = max(1 - sum(vf, 2), 0);
fname = [fpath 'svf.inp.' expnr];
fileID = fopen(fname,'W');
fprintf(fileID, '# sky view factors\n');
fclose(fileID);
dlmwrite(fname, svf, '-append','delimiter',' ','precision','%4f')

%albedo = ones(size(svf))*alb;
S = directShortwave(F,V,nsun,I,res);
K = reflectedShortwave(S,Dsky,vf,svf,albedos);

fname = [fpath 'netsw.inp.' expnr];
fileID = fopen(fname, 'w');
fprintf(fileID,'# %4s\n','net shortwave on facets [W/m2] (including reflections and diffusive)');
fprintf(fileID,'%6d\n', K);
fclose(fileID);