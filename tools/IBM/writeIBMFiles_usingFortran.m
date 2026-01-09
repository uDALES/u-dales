% Assumes the following variables exist in the workspace:
% TR: triangulation describing the geometry.
% itot, jtot, ktot, dx, dy: number of grid points along x, y, and
%   z-directions; and grid spacing along x and y. Assumes uniform meshing
%   along x and y.
% fpath: directory to write files to
% stl_ground: true if there are facets on ground level (recommended).
% diag_neighbs: determines whether non-directly-adjacent cells count
%   as neighbours. In situations with facets aligning nicely with the grid,
%   this probably is not needed. However, if not used then it is possible
%   that some facet sections will not find an adjacent fluid boundary point
%   to give the flux to, thus their effect will not be felt by the flow.
%   suggestion: leave off, and see if area_facets_c == area_facets.
% periodic_x, periodic_y: boolean for periodic boundary conditions in x/y.

% It will write out the following files:
% solid_u/v/w: list of indices of solid points (inside the geometry).
% fluid_boundary_u/v/w: list of indices of fluid points that have solid neighbours.
% facet_sections_u/v/w: list of facet section information: facet id,

nfcts = size(TR,1);
Dir_ray_u = [0 0 1];
Dir_ray_v = [0 0 1];
Dir_ray_w = [0 0 1];
Dir_ray_c = [0 0 1];
tol_mypoly = 5e-4;
ldebugplot = true;

currentPath = pwd;
stk = dbstack; activeFilename = which(stk(1).file);
[folder, ~, ~] = fileparts(activeFilename);


disp('Determining both fluid/solid points and facet sections using Fortran.')

%%Compile Fortran script

n_threads = 8;
in_mypoly_fortran_path = [folder '/in_mypoly_fortran/'];
addpath(in_mypoly_fortran_path)
cd(in_mypoly_fortran_path);
system('gfortran -O3 -fopenmp in_mypoly_functions.f90 boundaryMasking.f90 matchFacetsCells.f90 IBM_preproc_io.f90 IBM_preproc_main.f90 -o pre.exe');
copyfile('pre.exe', fpath)
delete pre.exe in_mypoly_functions.mod boundaryMasking.mod matchFacets2Cells.mod IBM_preproc_io.mod


%%Write input files to run Fortran script

% Wriring parameters for Fortran run
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
fprintf(fileID,'%d %d %d %d\n',[stl_ground diag_neighbs periodic_x periodic_y]);
fclose(fileID);

% Wriring z-grid information for Fortran run
fileID = fopen([fpath 'zhgrid.txt'],'w');
fprintf(fileID,'%15.10f\n',zgrid_w');
fclose(fileID);

fileID = fopen([fpath 'zfgrid.txt'],'w');
fprintf(fileID,'%15.10f\n',zgrid_c');
fclose(fileID);

% Wriring STL geometry information for Fortran run
fileID = fopen([fpath 'vertices.txt'],'w');
fprintf(fileID,'%15.10f %15.10f %15.10f\n',TR.Points');
fclose(fileID);

fileID = fopen([fpath 'faces.txt'],'w');
fprintf(fileID,'%8d %8d %8d %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n',[TR.ConnectivityList TR.incenter TR.faceNormal]');
fclose(fileID);


%%Run Fortran script
cd(fpath)
if lwindows
    system('pre.exe');
else
    system('./pre.exe');
end
delete pre.exe inmypoly_inp_info.txt faces.txt vertices.txt zfgrid.txt zhgrid.txt;
cd(currentPath)


%% Creating ncounts variable necessary for updating namoptions
ncounts = readmatrix('info_fort.txt', 'Range', [2,1]);

%% Plot

if ldebugplot

    % u-grid
    solid_ijk_u = readmatrix([fpath 'solid_u.txt'],'Range', 2);
    fluid_IB_ijk_u = readmatrix([fpath 'fluid_boundary_u.txt'],'Range', 2);
    fluid_IB_xyz_u = [xgrid_u(fluid_IB_ijk_u(:,1))', ygrid_u(fluid_IB_ijk_u(:,2))', zgrid_u(fluid_IB_ijk_u(:,3))'];

    % v-grid
    solid_ijk_v = readmatrix([fpath 'solid_v.txt'],'Range', 2);
    fluid_IB_ijk_v = readmatrix([fpath 'fluid_boundary_v.txt'],'Range', 2);
    fluid_IB_xyz_v = [xgrid_v(fluid_IB_ijk_v(:,1))', ygrid_v(fluid_IB_ijk_v(:,2))', zgrid_v(fluid_IB_ijk_v(:,3))'];

    % w-grid
    solid_ijk_w = readmatrix([fpath 'solid_w.txt'],'Range', 2);
    fluid_IB_ijk_w = readmatrix([fpath 'fluid_boundary_w.txt'],'Range', 2);
    fluid_IB_xyz_w = [xgrid_w(fluid_IB_ijk_w(:,1))', ygrid_w(fluid_IB_ijk_w(:,2))', zgrid_w(fluid_IB_ijk_w(:,3))'];

    % c-grid
    solid_ijk_c = readmatrix([fpath 'solid_c.txt'],'Range', 2);
    fluid_IB_ijk_c = readmatrix([fpath 'fluid_boundary_c.txt'],'Range', 2);
    fluid_IB_xyz_c = [xgrid_c(fluid_IB_ijk_c(:,1))', ygrid_c(fluid_IB_ijk_c(:,2))', zgrid_c(fluid_IB_ijk_c(:,3))'];


    %%Plot solid points
    fig = figure('Visible', 'off');
    
    patch('Faces', TR.ConnectivityList, 'Vertices', TR.Points, 'FaceColor', ones(3,1)*0.85, 'FaceAlpha', 1)
    hold on
    incenters = TR.incenter;
    faceNormals = TR.faceNormal;
    quiver3(incenters(:,1), incenters(:,2), incenters(:,3), faceNormals(:,1), faceNormals(:,2), faceNormals(:,3), 0)
    view(3)
    axis equal tight
    
    %%solid points
    if ~isempty(solid_ijk_u)
       scatter3(xgrid_u(solid_ijk_u(:,1)),ygrid_u(solid_ijk_u(:,2)),zgrid_u(solid_ijk_u(:,3)),10,[0,0,1],'filled')
    end
    if ~isempty(solid_ijk_v)
       scatter3(xgrid_v(solid_ijk_v(:,1)),ygrid_v(solid_ijk_v(:,2)),zgrid_v(solid_ijk_v(:,3)),10,[0,0,1],'filled')
    end
    if ~isempty(solid_ijk_w)
       scatter3(xgrid_w(solid_ijk_w(:,1)),ygrid_w(solid_ijk_w(:,2)),zgrid_w(solid_ijk_w(:,3)),10,[0,0,1],'filled')
    end
    if ~isempty(solid_ijk_c)
       scatter3(xgrid_c(solid_ijk_c(:,1)),ygrid_c(solid_ijk_c(:,2)),zgrid_c(solid_ijk_c(:,3)),10,[0,0,1],'filled')
    end
    
    set(fig, 'Visible', 'on');
    savefig(fig, 'solid_points.fig');  % Save figure
    close(fig);  % Close the invisible figure
    
    
    %%Plot fluid boundary points
    fig = figure('Visible', 'off');
    
    patch('Faces', TR.ConnectivityList, 'Vertices', TR.Points, 'FaceColor', ones(3,1)*0.85, 'FaceAlpha', 1)
    hold on
    view(3)
    axis equal tight
    
    %%fluid boundary points
    scatter3(fluid_IB_xyz_u(:,1),fluid_IB_xyz_u(:,2),fluid_IB_xyz_u(:,3),10,[0,0,1],'filled')
    scatter3(fluid_IB_xyz_v(:,1),fluid_IB_xyz_v(:,2),fluid_IB_xyz_v(:,3),10,[0,0,1],'filled')
    scatter3(fluid_IB_xyz_w(:,1),fluid_IB_xyz_w(:,2),fluid_IB_xyz_w(:,3),10,[0,0,1],'filled')
    scatter3(fluid_IB_xyz_c(:,1),fluid_IB_xyz_c(:,2),fluid_IB_xyz_c(:,3),10,[0,0,1],'filled')
    
    set(fig, 'Visible', 'on');
    savefig(fig, 'fluid_boundary_points.fig');  % Save figure
    close(fig);  % Close the invisible figure
end
