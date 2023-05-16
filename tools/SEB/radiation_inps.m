%% Setting constants

zen =45 ;
solaz = 235;
xaz = 90;
az = solaz - xaz;
nsun = [sin(zen)*cos(az), -sin(zen)*sin(az), cos(zen)];
F = TR.ConnectivityList;
V = TR.Points;
I = 800;
Dsky = 100;
res = 0.01;
alb = 0.5;

%% Assumes view factors have been calculated using View3D

vf = dlmread([fpath 'vf.txt'], ' ', 2, 0);

vf(end,:) = [];

%%

lwritenetcdf = 1;

nfcts = size(vf,2);

if lwritenetcdf
    ncid = netcdf.create([fpath 'vf.nc.inp.' expnr], 'NC_WRITE');
    dimidrow = netcdf.defDim(ncid,'rows', nfcts);
    dimidcol = netcdf.defDim(ncid,'columns', nfcts);
    varid = netcdf.defVar(ncid,'view factor','NC_FLOAT',[dimidrow dimidcol]);
    netcdf.endDef(ncid);
    netcdf.putVar(ncid,varid,vf);
    netcdf.close(ncid);
end

%% Remove view factors below a certain tolerance

tol = 0.001;
vftol = zeros(size(vf));
vftol(vf > tol) = vf(vf > tol);
vfsparse = sparse(double(vftol));
[i,j,s] = find(vfsparse);
%prec = '%.6f';
fID = fopen([fpath 'vfsparse.inp.' expnr], 'w');
fprintf(fID, '%d %d %.6f \n', [i, j, s]');


%%

svf = max(1 - sum(vfsparse, 2), 0);
fname = [fpath 'svf.inp.' expnr];
fileID = fopen(fname,'W');
fprintf(fileID, '# sky view factors\n');
fclose(fileID);
dlmwrite(fname, svf, '-append','delimiter',' ','precision','%4f')

albedo = ones(size(svf))*alb;
S = directShortwave(F,V,nsun,I,res);
K = reflectedShortwave(S,Dsky,vf,svf,albedo);

fname = [fpath 'netsw.inp.' expnr];
fileID = fopen(fname, 'w');
fprintf(fileID,'# %4s\n','net shortwave on facets [W/m2] (including reflections and diffusive)');
fprintf(fileID,'%6d\n', K);
fclose(fileID);