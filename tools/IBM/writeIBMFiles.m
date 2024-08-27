% addpath('./inpolyhedron/')
% addpath('./point2trimesh/')
% addpath('./in_mypoly/')

% Assumes the following variables exist in the workspace:
% TR: triangulation describing the geometry.
% xgrid_u, xgrid_v, xgrid_w, xgrid_c: the x coordinates on which u,v,w,
%   and scalars are defined (and similar arrays for y and z coords).
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
%   Calculated using inpolyhedron if lmypoly=false, or using Dipanjan's routine
%   if lmypoly=true.
% fluid_boundary_u/v/w: list of indices of fluid points that have solid neighbours.
% facet_sections_u/v/w: list of facet section information: facet id,
%   section area, fluid boundary point id, distance to section/surface as a whole (_2)

nfcts = size(TR,1);
Dir_ray_u = [0 0 1];
Dir_ray_v = [0 0 1];
Dir_ray_w = [0 0 1];
Dir_ray_c = [0 0 1];
tol_mypoly = 5e-4;
tol_inpolyhedron = 5e-2; %z0=5e-2; tol_inpolyhedron = ceil(z0*exp(1)*1000)/1000;
lcheckdist = false;

currentPath = pwd;
stk = dbstack; activeFilename = which(stk(1).file);
[folder, ~, ~] = fileparts(activeFilename);

% Formatting assumes: #facets < 10 million, #boundary points < 1 billion, 
% section area < 1000 m^2 (rounded to cm^2), and distance < 1000m
facet_sections_header = ' # facet      area flux point distance\n';
facet_sections_format = '%8d %9.4f %10d %8.4f\n';

%% Fluid/solid points
if lmypolyfortran
    disp('Determining fluid/solid points using Fortran.')
    % Compile
    n_threads = 8;
    in_mypoly_fortran_path = [folder '/in_mypoly_fortran/'];
    addpath(in_mypoly_fortran_path)
    cd(in_mypoly_fortran_path);
    system('gfortran -O2 -fopenmp in_mypoly_functions.f90 ibm_necessary_functions.f90 IBM_flagging.f90 -o pre.exe');
    copyfile('pre.exe', fpath)
    delete pre.exe in_mypoly_functions.mod ibm_necessary_functions.mod

    % Write input files
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

    fileID = fopen([fpath 'zhgrid.txt'],'w');
    fprintf(fileID,'%15.10f\n',zgrid_w');
    fclose(fileID);

    fileID = fopen([fpath 'zfgrid.txt'],'w');
    fprintf(fileID,'%15.10f\n',zgrid_c');
    fclose(fileID);

    fileID = fopen([fpath 'vertices.txt'],'w');
    fprintf(fileID,'%15.10f %15.10f %15.10f\n',TR.Points');
    fclose(fileID);

    fileID = fopen([fpath 'faces.txt'],'w');
    fprintf(fileID,'%8d %8d %8d %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n',[TR.ConnectivityList TR.incenter TR.faceNormal]');
    fclose(fileID);

    % Run
    cd(fpath)
    if lwindows
        system('pre.exe');
    else
        system('./pre.exe');
    end
    delete pre.exe inmypoly_inp_info.txt faces.txt vertices.txt zfgrid.txt zhgrid.txt;
    cd(currentPath)

    %% u-grid
    solid_ijk_u = readmatrix([fpath 'solid_u.txt'],'Range', 2);
    fluid_IB_ijk_u = readmatrix([fpath 'fluid_boundary_u.txt'],'Range', 2);
    %fluid_IB_ijk_u = sortrows(fluid_IB_ijk_u,3);
    fluid_boundary_u = fluid_IB_ijk_u;
    fluid_IB_xyz_u = [xgrid_u(fluid_IB_ijk_u(:,1))', ygrid_u(fluid_IB_ijk_u(:,2))', zgrid_u(fluid_IB_ijk_u(:,3))'];
    solid_IB_ijk_u = readmatrix([fpath 'solid_boundary_u.txt'],'Range', 2);
    %solid_IB_ijk_u = sortrows(solid_IB_ijk_u,3);
    if isempty(solid_IB_ijk_u)
        solid_IB_xyz_u = [];
    else
        solid_IB_xyz_u = [xgrid_u(solid_IB_ijk_u(:,1))', ygrid_u(solid_IB_ijk_u(:,2))', zgrid_u(solid_IB_ijk_u(:,3))'];
    end

    solid_u = false(itot,jtot,ktot);
    for n=1:size(solid_ijk_u,1)
        solid_u(solid_ijk_u(n,1), solid_ijk_u(n,2), solid_ijk_u(n,3)) = true;
    end
    fluid_u = ~solid_u;

    fluid_IB_u = false(itot,jtot,ktot);
    for n=1:size(fluid_IB_ijk_u,1)
        fluid_IB_u(fluid_IB_ijk_u(n,1), fluid_IB_ijk_u(n,2), fluid_IB_ijk_u(n,3)) = true;
    end

    solid_IB_u = false(itot,jtot,ktot);
    for n=1:size(solid_IB_ijk_u,1)
        solid_IB_u(solid_IB_ijk_u(n,1), solid_IB_ijk_u(n,2), solid_IB_ijk_u(n,3)) = true;
    end

%     formatSpec = '%d';
%     fileID = fopen([fpath 'flag_u.txt'],'r');
%     flag_u = fscanf(fileID,formatSpec);
%     fclose(fileID);
%     count = 1;
%     solid_u = false(itot,jtot,ktot);
%     for iy = 1:jtot
%         for iz = 1:ktot
%             for ix = 1:itot
%                 if (flag_u(count) == 1)
%                     solid_u(ix,iy,iz) = true;
%                 else
%                     solid_u(ix,iy,iz) = false;
%                 end
%                 count = count+1;
%             end
%         end
%     end
% 
%     fileID = fopen([fpath 'fluid_IB_u.txt'],'r');
%     fluid_IB_u_fort = fscanf(fileID,formatSpec);
%     fclose(fileID);
%     count = 1;
%     for iy = 1:jtot
%         for iz = 1:ktot
%             for ix = 1:itot
%                 if (fluid_IB_u_fort(count) == 1)
%                     fluid_IB_u(ix,iy,iz) = true;
%                 else
%                     fluid_IB_u(ix,iy,iz) = false;
%                 end
%                 count = count+1;
%             end
%         end
%     end
% 
%     fileID = fopen([fpath 'solid_IB_u.txt'],'r');
%     solid_IB_u_fort = fscanf(fileID,formatSpec);
%     fclose(fileID);
%     count = 1;
%     for iy = 1:jtot
%         for iz = 1:ktot
%             for ix = 1:itot
%                 if (solid_IB_u_fort(count) == 1)
%                     solid_IB_u(ix,iy,iz) = true;
%                 else
%                     solid_IB_u(ix,iy,iz) = false;
%                 end
%                 count = count+1;
%             end
%         end
%     end



    %% v-grid
    solid_ijk_v = readmatrix([fpath 'solid_v.txt'],'Range', 2);
    fluid_IB_ijk_v = readmatrix([fpath 'fluid_boundary_v.txt'],'Range', 2);
    %fluid_IB_ijk_v = sortrows(fluid_IB_ijk_v,3);
    fluid_boundary_v = fluid_IB_ijk_v;
    fluid_IB_xyz_v = [xgrid_v(fluid_IB_ijk_v(:,1))', ygrid_v(fluid_IB_ijk_v(:,2))', zgrid_v(fluid_IB_ijk_v(:,3))'];
    solid_IB_ijk_v = readmatrix([fpath 'solid_boundary_v.txt'],'Range', 2);
    %solid_IB_ijk_v = sortrows(solid_IB_ijk_v,3);
    if isempty(solid_IB_ijk_v)
        solid_IB_xyz_v = [];
    else
        solid_IB_xyz_v = [xgrid_v(solid_IB_ijk_v(:,1))', ygrid_v(solid_IB_ijk_v(:,2))', zgrid_v(solid_IB_ijk_v(:,3))'];
    end
    
    solid_v = false(itot,jtot,ktot);
    for n=1:size(solid_ijk_v,1)
        solid_v(solid_ijk_v(n,1), solid_ijk_v(n,2), solid_ijk_v(n,3)) = true;
    end
    fluid_v = ~solid_v;

    fluid_IB_v = false(itot,jtot,ktot);
    for n=1:size(fluid_IB_ijk_v,1)
        fluid_IB_v(fluid_IB_ijk_v(n,1), fluid_IB_ijk_v(n,2), fluid_IB_ijk_v(n,3)) = true;
    end

    solid_IB_v = false(itot,jtot,ktot);
    for n=1:size(solid_IB_ijk_v,1)
        solid_IB_v(solid_IB_ijk_v(n,1), solid_IB_ijk_v(n,2), solid_IB_ijk_v(n,3)) = true;
    end

%     formatSpec = '%d';
%     fileID = fopen([fpath 'flag_v.txt'],'r');
%     flag_v = fscanf(fileID,formatSpec);
%     fclose(fileID);
%     count = 1;
%     for iy = 1:jtot
%         for iz = 1:ktot
%             for ix = 1:itot
%                 if (flag_v(count) == 1)
%                     solid_v(ix,iy,iz) = true;
%                 else
%                     solid_v(ix,iy,iz) = false;
%                 end
%                 count = count+1;
%             end
%         end
%     end
% 
%     fileID = fopen([fpath 'fluid_IB_v.txt'],'r');
%     fluid_IB_v_fort = fscanf(fileID,formatSpec);
%     fclose(fileID);
%     count = 1;
%     for iy = 1:jtot
%         for iz = 1:ktot
%             for ix = 1:itot
%                 if (fluid_IB_v_fort(count) == 1)
%                     fluid_IB_v(ix,iy,iz) = true;
%                 else
%                     fluid_IB_v(ix,iy,iz) = false;
%                 end
%                 count = count+1;
%             end
%         end
%     end
% 
%     fileID = fopen([fpath 'solid_IB_v.txt'],'r');
%     solid_IB_v_fort = fscanf(fileID,formatSpec);
%     fclose(fileID);
%     count = 1;
%     for iy = 1:jtot
%         for iz = 1:ktot
%             for ix = 1:itot
%                 if (solid_IB_v_fort(count) == 1)
%                     solid_IB_v(ix,iy,iz) = true;
%                 else
%                     solid_IB_v(ix,iy,iz) = false;
%                 end
%                 count = count+1;
%             end
%         end
%     end

    %% w-grid
    % w
    solid_ijk_w = readmatrix([fpath 'solid_w.txt'],'Range', 2);
    fluid_IB_ijk_w = readmatrix([fpath 'fluid_boundary_w.txt'],'Range', 2);
    %fluid_IB_ijk_w = sortrows(fluid_IB_ijk_w,3);
    fluid_boundary_w = fluid_IB_ijk_w;
    fluid_IB_xyz_w = [xgrid_w(fluid_IB_ijk_w(:,1))', ygrid_w(fluid_IB_ijk_w(:,2))', zgrid_w(fluid_IB_ijk_w(:,3))'];
    solid_IB_ijk_w = readmatrix([fpath 'solid_boundary_w.txt'],'Range', 2);
    %solid_IB_ijk_w = sortrows(solid_IB_ijk_w,3);
    if isempty(solid_IB_ijk_w)
        solid_IB_xyz_w = [];
    else
        solid_IB_xyz_w = [xgrid_w(solid_IB_ijk_w(:,1))', ygrid_w(solid_IB_ijk_w(:,2))', zgrid_w(solid_IB_ijk_w(:,3))'];
    end
    
    % Convert from sparse (list) format to mask (3D array)
    solid_w = false(itot,jtot,ktot);
    for n=1:size(solid_ijk_w,1)
        solid_w(solid_ijk_w(n,1), solid_ijk_w(n,2), solid_ijk_w(n,3)) = true;
    end
    fluid_w = ~solid_w;

    fluid_IB_w = false(itot,jtot,ktot);
    for n=1:size(fluid_IB_ijk_w,1)
        fluid_IB_w(fluid_IB_ijk_w(n,1), fluid_IB_ijk_w(n,2), fluid_IB_ijk_w(n,3)) = true;
    end

    solid_IB_w = false(itot,jtot,ktot);
    for n=1:size(solid_IB_ijk_w,1)
        solid_IB_w(solid_IB_ijk_w(n,1), solid_IB_ijk_w(n,2), solid_IB_ijk_w(n,3)) = true;
    end
    
%     formatSpec = '%d';
%     fileID = fopen([fpath 'flag_w.txt'],'r');
%     flag_w = fscanf(fileID,formatSpec);
%     fclose(fileID);
%     count = 1;
%     for iy = 1:jtot
%         for iz = 1:ktot
%             for ix = 1:itot
%                 if (flag_w(count) == 1)
%                             solid_w(ix,iy,iz) = true;
%                 else
%                             solid_w(ix,iy,iz) = false;
%                 end
%                 count = count+1;
%             end
%         end
%     end
% 
%     fileID = fopen([fpath 'fluid_IB_w.txt'],'r');
%     fluid_IB_w_fort = fscanf(fileID,formatSpec);
%     fclose(fileID);
%     count = 1;
%     for iy = 1:jtot
%         for iz = 1:ktot
%             for ix = 1:itot
%                 if (fluid_IB_w_fort(count) == 1)
%                             fluid_IB_w(ix,iy,iz) = true;
%                 else
%                             fluid_IB_w(ix,iy,iz) = false;
%                 end
%                 count = count+1;
%             end
%         end
%     end
% 
%     fileID = fopen([fpath 'solid_IB_w.txt'],'r');
%     solid_IB_w_fort = fscanf(fileID,formatSpec);
%     fclose(fileID);
%     count = 1;
%     for iy = 1:jtot
%         for iz = 1:ktot
%             for ix = 1:itot
%                 if (solid_IB_w_fort(count) == 1)
%                             solid_IB_w(ix,iy,iz) = true;
%                 else
%                             solid_IB_w(ix,iy,iz) = false;
%                 end
%                 count = count+1;
%             end
%         end
%     end

    %% c-grid
    solid_ijk_c = readmatrix([fpath 'solid_c.txt'],'Range', 2);
    fluid_IB_ijk_c = readmatrix([fpath 'fluid_boundary_c.txt'],'Range', 2);
    %fluid_IB_ijk_c = sortrows(fluid_IB_ijk_c,3);
    fluid_boundary_c = fluid_IB_ijk_c;
    fluid_IB_xyz_c = [xgrid_c(fluid_IB_ijk_c(:,1))', ygrid_c(fluid_IB_ijk_c(:,2))', zgrid_c(fluid_IB_ijk_c(:,3))'];
    solid_IB_ijk_c = readmatrix([fpath 'solid_boundary_c.txt'],'Range', 2);
    %solid_IB_ijk_c = sortrows(solid_IB_ijk_c,3);
    if isempty(solid_IB_ijk_c)
        solid_IB_xyz_c = [];
    else
        solid_IB_xyz_c = [xgrid_c(solid_IB_ijk_c(:,1))', ygrid_c(solid_IB_ijk_c(:,2))', zgrid_c(solid_IB_ijk_c(:,3))'];
    end
    
    % Convert from sparse (list) format to mask (3D array)
    solid_c = false(itot,jtot,ktot);
    for n=1:size(solid_ijk_c,1)
        solid_c(solid_ijk_c(n,1), solid_ijk_c(n,2), solid_ijk_c(n,3)) = true;
    end
    fluid_c = ~solid_c;

    fluid_IB_c = false(itot,jtot,ktot);
    for n=1:size(fluid_IB_ijk_c,1)
        fluid_IB_c(fluid_IB_ijk_c(n,1), fluid_IB_ijk_c(n,2), fluid_IB_ijk_c(n,3)) = true;
    end

    solid_IB_c = false(itot,jtot,ktot);
    for n=1:size(solid_IB_ijk_c,1)
        solid_IB_c(solid_IB_ijk_c(n,1), solid_IB_ijk_c(n,2), solid_IB_ijk_c(n,3)) = true;
    end

%     formatSpec = '%d';
%     fileID = fopen([fpath 'flag_c.txt'],'r');
%     flag_c = fscanf(fileID,formatSpec);
%     fclose(fileID);
%     count = 1;
%     for iy = 1:jtot
%         for iz = 1:ktot
%             for ix = 1:itot
%                 if (flag_c(count) == 1)
%                             solid_c(ix,iy,iz) = true;
%                 else
%                             solid_c(ix,iy,iz) = false;
%                 end
%                 count = count+1;
%             end
%         end
%     end

%     fileID = fopen([fpath 'fluid_IB_c.txt'],'r');
%     fluid_IB_c_fort = fscanf(fileID,formatSpec);
%     fclose(fileID);
%     count = 1;
%     for iy = 1:jtot
%         for iz = 1:ktot
%             for ix = 1:itot
%                 if (fluid_IB_c_fort(count) == 1)
%                             fluid_IB_c(ix,iy,iz) = true;
%                 else
%                             fluid_IB_c(ix,iy,iz) = false;
%                 end
%                 count = count+1;
%             end
%         end
%     end
% 
%     fileID = fopen([fpath 'solid_IB_c.txt'],'r');
%     solid_IB_c_fort = fscanf(fileID,formatSpec);
%     fclose(fileID);
%     count = 1;
%     for iy = 1:jtot
%         for iz = 1:ktot
%             for ix = 1:itot
%                 if (solid_IB_c_fort(count) == 1)
%                             solid_IB_c(ix,iy,iz) = true;
%                 else
%                             solid_IB_c(ix,iy,iz) = false;
%                 end
%                 count = count+1;
%             end
%         end
%     end

else
    %% Solid/fluid detection
    if lmypoly
        max_height = max(TR.Points(:,3)) + tol_mypoly;
        L_char = max_facet_size(TR.Points,TR.ConnectivityList) + tol_mypoly;
    end

    if lmypoly
        disp('Determining fluid/solid points using MATLAB.')
        solid_u = in_grid_mypoly(TR.Points,TR.ConnectivityList, ...
            TR.incenter,TR.faceNormal,xgrid_u,ygrid_u,zgrid_u,Dir_ray_u,L_char,max_height,tol_mypoly);
    else
        disp('Determining fluid/solid points using MATLAB (inpolyhedron).')
        solid_u = inpolyhedron(TR.ConnectivityList, TR.Points, xgrid_u, ygrid_u, zgrid_u, ...
           'FACENORMALS', TR.faceNormal, 'TOL', tol_inpolyhedron);
        solid_u = permute(solid_u, [2 1 3]);
        %solid_u(1,:,:) = 0;

        % check for fluid points very close to surface
        if lcheckdist
            fluid_u = ~solid_u;
            [fluid_IB_u, ~] = getBoundaryCells(xgrid_u, ygrid_u, zgrid_u, fluid_u, solid_u, diag_neighbs);
            if (stl_ground)
                fluid_u_1 = fluid_u(:,:,1);
                fluid_IB_u_1 = fluid_IB_u(:,:,1);
                fluid_IB_u_1(fluid_u_1) = true;
                fluid_IB_u(:,:,1) = fluid_IB_u_1;
            end

            [fluid_IB_i_u, fluid_IB_j_u, fluid_IB_k_u] = ind2sub(size(fluid_IB_u), find(fluid_IB_u));
            fluid_IB_xyz_u = [xgrid_u(fluid_IB_i_u)', ygrid_u(fluid_IB_j_u)', zgrid_u(fluid_IB_k_u)'];

            [dists, ~, ~, ~] = point2trimesh('Faces', TR.ConnectivityList, 'Vertices', TR.Points, ...
                'QueryPoints', fluid_IB_xyz_u, 'UseSubSurface', false);

            for n=1:size(fluid_IB_xyz_u,1)
                if abs(dists(n)) < tol_inpolyhedron
                    [~, i] = ismember(fluid_IB_xyz_u(n,1), xgrid_u);
                    [~, j] = ismember(fluid_IB_xyz_u(n,2), ygrid_u);
                    [~, k] = ismember(fluid_IB_xyz_u(n,3), zgrid_u);
                    solid_u(i,j,k) = true;
                end
            end
        end
    end

    % check for acute corners - single column of solid cells
    solid_u_copy = solid_u;
    for k=1:ktot
        for j=1:jtot
            for i=1:itot
                if solid_u(i,j,k)
                    if i ~= 1 && i ~= itot
                        if ~solid_u(i+1,j,k) && ~solid_u(i-1,j,k)
                            solid_u_copy(i,j,k) = false;
                        end
                    end
                    if j ~= 1 && j ~= jtot
                        if ~solid_u(i,j+1,k) && ~solid_u(i,j-1,k)
                            solid_u_copy(i,j,k) = false;
                        end
                    end
                end
            end
        end
    end
    solid_u = solid_u_copy;

    % Convert from mask (3D array) to sparse (list) format
    [solid_i_u, solid_j_u, solid_k_u] = ind2sub(size(solid_u), find(solid_u));
    solid_ijk_u = [solid_i_u, solid_j_u, solid_k_u];
    filename_u = [fpath 'solid_u.txt'];
    fileID_u = fopen(filename_u,'W');
    fprintf(fileID_u, '# position (i,j,k)\n');
    fclose(fileID_u);
    dlmwrite(filename_u, solid_ijk_u, '-append','delimiter',' ');
    disp('Written solid_u.txt')
    fluid_u = ~solid_u;

    if lmypoly
        solid_v = in_grid_mypoly(TR.Points,TR.ConnectivityList,TR.incenter,TR.faceNormal,xgrid_v,ygrid_v,zgrid_v,Dir_ray_v,L_char,max_height,tol_mypoly);
    else
        solid_v = inpolyhedron(TR.ConnectivityList, TR.Points, xgrid_v, ygrid_v, zgrid_v, ...
            'FACENORMALS', TR.faceNormal, 'TOL', tol_inpolyhedron);
        solid_v = permute(solid_v, [2 1 3]);
        %solid_v(:,1,:) = 0;

        if lcheckdist
            fluid_v = ~solid_v;
            [fluid_IB_v, ~] = getBoundaryCells(xgrid_v, ygrid_v, zgrid_v, fluid_v, solid_v, diag_neighbs);
            if (stl_ground)
                fluid_v_1 = fluid_v(:,:,1);
                fluid_IB_v_1 = fluid_IB_v(:,:,1);
                fluid_IB_v_1(fluid_v_1) = true;
                fluid_IB_v(:,:,1) = fluid_IB_v_1;
            end

            [fluid_IB_i_v, fluid_IB_j_v, fluid_IB_k_v] = ind2sub(size(fluid_IB_v), find(fluid_IB_v));
            fluid_IB_xyz_v = [xgrid_v(fluid_IB_i_v)', ygrid_v(fluid_IB_j_v)', zgrid_v(fluid_IB_k_v)'];

            [dists, ~, ~, ~] = point2trimesh('Faces', TR.ConnectivityList, 'Vertices', TR.Points, ...
                'QueryPoints', fluid_IB_xyz_v, 'UseSubSurface', false);

            for n=1:size(fluid_IB_xyz_v,1)
                if abs(dists(n)) < tol_inpolyhedron
                    [~, i] = ismember(fluid_IB_xyz_v(n,1), xgrid_v);
                    [~, j] = ismember(fluid_IB_xyz_v(n,2), ygrid_v);
                    [~, k] = ismember(fluid_IB_xyz_v(n,3), zgrid_v);
                    solid_v(i,j,k) = true;
                end
            end
        end
    end

    % check for acute corners - single column of solid cells
    solid_v_copy = solid_v;
    for k=1:ktot
        for j=1:jtot
            for i=1:itot
                if solid_v(i,j,k)
                    if i ~= 1 && i ~= itot
                        if ~solid_v(i+1,j,k) && ~solid_v(i-1,j,k)
                            solid_v_copy(i,j,k) = false;
                        end
                    end
                    if j ~= 1 && j ~= jtot
                        if ~solid_v(i,j+1,k) && ~solid_v(i,j-1,k)
                            solid_v_copy(i,j,k) = false;
                        end
                    end
                end
            end
        end
    end
    solid_v = solid_v_copy;

    [solid_i_v, solid_j_v, solid_k_v] = ind2sub(size(solid_v), find(solid_v));
    solid_ijk_v = [solid_i_v, solid_j_v, solid_k_v];
    filename_v = [fpath 'solid_v.txt'];
    fileID_v = fopen(filename_v,'W');
    fprintf(fileID_v, '# position (i,j,k)\n');
    fclose(fileID_v);
    dlmwrite(filename_v, solid_ijk_v, '-append','delimiter',' ');
    disp('Written solid_v.txt')
    fluid_v = ~solid_v;

    if lmypoly
        solid_w = in_grid_mypoly(TR.Points,TR.ConnectivityList, ...
            TR.incenter,TR.faceNormal,xgrid_w,ygrid_w,zgrid_w,Dir_ray_w,L_char,max_height,tol_mypoly);
    else
        solid_w = inpolyhedron(TR.ConnectivityList, TR.Points, xgrid_w, ygrid_w, zgrid_w, ...
            'FACENORMALS', TR.faceNormal, 'TOL', tol_inpolyhedron);
        solid_w = permute(solid_w, [2 1 3]);

        if lcheckdist
            if (stl_ground)
                solid_w_b = solid_w;
                solid_w_b(:,:,1) = true;
            end
            fluid_w = ~solid_w;

            [fluid_IB_w, ~] = getBoundaryCells(xgrid_w, ygrid_w, zgrid_w, fluid_w, solid_w_b, diag_neighbs);
            fluid_IB_w(:,:,1) = 0; % Bottom is always solid

            [fluid_IB_i_w, fluid_IB_j_w, fluid_IB_k_w] = ind2sub(size(fluid_IB_w), find(fluid_IB_w));
            fluid_IB_xyz_w = [xgrid_w(fluid_IB_i_w)', ygrid_w(fluid_IB_j_w)', zgrid_w(fluid_IB_k_w)'];

            [dists, ~, ~, ~] = point2trimesh('Faces', TR.ConnectivityList, 'Vertices', TR.Points, ...
                'QueryPoints', fluid_IB_xyz_w, 'UseSubSurface', false);

            for n=1:size(fluid_IB_xyz_w,1)
                if abs(dists(n)) < tol_inpolyhedron
                    [~, i] = ismember(fluid_IB_xyz_w(n,1), xgrid_w);
                    [~, j] = ismember(fluid_IB_xyz_w(n,2), ygrid_w);
                    [~, k] = ismember(fluid_IB_xyz_w(n,3), zgrid_w);
                    solid_w(i,j,k) = true;
                end
            end
        end
    end

    % check for acute corners - single column of solid cells
    solid_w_copy = solid_w;
    for k=1:ktot
        for j=1:jtot
            for i=1:itot
                if solid_w(i,j,k)
                    if i ~= 1 && i ~= itot
                        if ~solid_w(i+1,j,k) && ~solid_w(i-1,j,k)
                            solid_w_copy(i,j,k) = false;
                        end
                    end
                    if j ~= 1 && j ~= jtot
                        if ~solid_w(i,j+1,k) && ~solid_w(i,j-1,k)
                            solid_w_copy(i,j,k) = false;
                        end
                    end
                end
            end
        end
    end
    solid_w = solid_w_copy;

    [solid_i_w, solid_j_w, solid_k_w] = ind2sub(size(solid_w), find(solid_w));
    solid_ijk_w = [solid_i_w, solid_j_w, solid_k_w];
    filename_w = [fpath 'solid_w.txt'];
    fileID_w = fopen(filename_w,'W');
    fprintf(fileID_w, '# position (i,j,k)\n');
    fclose(fileID_w);
    dlmwrite(filename_w, solid_ijk_w, '-append','delimiter',' ');
    disp('Written solid_w.txt')

    if (stl_ground)
        solid_w_b = solid_w;
        solid_w_b(:,:,1) = true;
    end
    fluid_w = ~solid_w;

    %% Boundary points
    % u
    [fluid_IB_u, solid_IB_u] = getBoundaryCells(xgrid_u, ygrid_u, zgrid_u, fluid_u, solid_u, diag_neighbs);
    if (stl_ground)
        fluid_u_1 = fluid_u(:,:,1);
        fluid_IB_u_1 = fluid_IB_u(:,:,1);
        fluid_IB_u_1(fluid_u_1) = true;
        fluid_IB_u(:,:,1) = fluid_IB_u_1;
    end

    [fluid_IB_i_u, fluid_IB_j_u, fluid_IB_k_u] = ind2sub(size(fluid_IB_u), find(fluid_IB_u));
    fluid_IB_xyz_u = [xgrid_u(fluid_IB_i_u)', ygrid_u(fluid_IB_j_u)', zgrid_u(fluid_IB_k_u)'];

    [solid_IB_i_u, solid_IB_j_u, solid_IB_k_u] = ind2sub(size(solid_IB_u), find(solid_IB_u));
    solid_IB_xyz_u = [xgrid_u(solid_IB_i_u)', ygrid_u(solid_IB_j_u)', zgrid_u(solid_IB_k_u)'];

    fluid_IB_ijk_u = [fluid_IB_i_u, fluid_IB_j_u, fluid_IB_k_u];
    solid_IB_ijk_u = [solid_IB_i_u, solid_IB_j_u, solid_IB_k_u];

    fluid_boundary_u = fluid_IB_ijk_u;
    filename_u = [fpath 'fluid_boundary_u.txt'];
    fileID_u = fopen(filename_u,'W');
    fprintf(fileID_u, '# position (i,j,k)\n');
    fclose(fileID_u);
    dlmwrite(filename_u, fluid_boundary_u, '-append','delimiter',' ');
    disp('Written fluid_boundary_u.txt')

    solid_boundary_u = solid_IB_ijk_u;
    filename_u = [fpath 'solid_boundary_u.txt'];
    fileID_u = fopen(filename_u,'W');
    fprintf(fileID_u, '# position (i,j,k)\n');
    fclose(fileID_u);
    dlmwrite(filename_u, solid_boundary_u, '-append','delimiter',' ');
    disp('Written solid_boundary_u.txt')

    % v
    [fluid_IB_v, solid_IB_v] = getBoundaryCells(xgrid_v, ygrid_v, zgrid_v, fluid_v, solid_v, diag_neighbs);
    if (stl_ground)
        fluid_v_1 = fluid_v(:,:,1);
        fluid_IB_v_1 = fluid_IB_v(:,:,1);
        fluid_IB_v_1(fluid_v_1) = true;
        fluid_IB_v(:,:,1) = fluid_IB_v_1;
    end

    [fluid_IB_i_v, fluid_IB_j_v, fluid_IB_k_v] = ind2sub(size(fluid_IB_v), find(fluid_IB_v));
    fluid_IB_xyz_v = [xgrid_v(fluid_IB_i_v)', ygrid_v(fluid_IB_j_v)', zgrid_v(fluid_IB_k_v)'];

    [solid_IB_i_v, solid_IB_j_v, solid_IB_k_v] = ind2sub(size(solid_IB_v), find(solid_IB_v));
    solid_IB_xyz_v = [xgrid_v(solid_IB_i_v)', ygrid_v(solid_IB_j_v)', zgrid_v(solid_IB_k_v)'];

    fluid_IB_ijk_v = [fluid_IB_i_v, fluid_IB_j_v, fluid_IB_k_v];
    solid_IB_ijk_v = [solid_IB_i_v, solid_IB_j_v, solid_IB_k_v];

    fluid_boundary_v = fluid_IB_ijk_v;
    filename_v = [fpath 'fluid_boundary_v.txt'];
    fileID_v = fopen(filename_v,'W');
    fprintf(fileID_v, '# position (i,j,k)\n');
    fclose(fileID_v);
    dlmwrite(filename_v, fluid_boundary_v, '-append','delimiter',' ');
    disp('Written fluid_boundary_v.txt')

    solid_boundary_v = solid_IB_ijk_v;
    filename_v = [fpath 'solid_boundary_v.txt'];
    fileID_v = fopen(filename_v,'W');
    fprintf(fileID_v, '# position (i,j,k)\n');
    fclose(fileID_v);
    dlmwrite(filename_v, solid_boundary_v, '-append','delimiter',' ');
    disp('Written solid_boundary_v.txt')

    % w
    [fluid_IB_w, solid_IB_w] = getBoundaryCells(xgrid_w, ygrid_w, zgrid_w, fluid_w, solid_w_b, diag_neighbs);
    fluid_IB_w(:,:,1) = 0; % Bottom is always solid

    [fluid_IB_i_w, fluid_IB_j_w, fluid_IB_k_w] = ind2sub(size(fluid_IB_w), find(fluid_IB_w));
    fluid_IB_xyz_w = [xgrid_w(fluid_IB_i_w)', ygrid_w(fluid_IB_j_w)', zgrid_w(fluid_IB_k_w)'];

    [solid_IB_i_w, solid_IB_j_w, solid_IB_k_w] = ind2sub(size(solid_IB_w), find(solid_IB_w));
    solid_IB_xyz_w = [xgrid_w(solid_IB_i_w)', ygrid_w(solid_IB_j_w)', zgrid_w(solid_IB_k_w)'];

    fluid_IB_ijk_w = [fluid_IB_i_w, fluid_IB_j_w, fluid_IB_k_w];
    solid_IB_ijk_w = [solid_IB_i_w, solid_IB_j_w, solid_IB_k_w];

    fluid_boundary_w = fluid_IB_ijk_w;
    filename_w = [fpath 'fluid_boundary_w.txt'];
    fileID_w = fopen(filename_w,'W');
    fprintf(fileID_w, '# position (i,j,k)\n');
    fclose(fileID_w);
    dlmwrite(filename_w, fluid_boundary_w, '-append','delimiter',' ');
    disp('Written fluid_boundary_w.txt')

    solid_boundary_w = solid_IB_ijk_w;
    filename_w = [fpath 'solid_boundary_w.txt'];
    fileID_w = fopen(filename_w,'W');
    fprintf(fileID_w, '# position (i,j,k)\n');
    fclose(fileID_w);
    dlmwrite(filename_w, solid_boundary_w, '-append','delimiter',' ');
    disp('Written solid_boundary_w.txt')

    % c
    if lmypoly
        solid_c = in_grid_mypoly(TR.Points,TR.ConnectivityList,TR.incenter,TR.faceNormal,xgrid_c,ygrid_c,zgrid_c,Dir_ray_c,L_char,max_height,tol_mypoly);
    else
        solid_c = inpolyhedron(TR.ConnectivityList, TR.Points, xgrid_c, ygrid_c, zgrid_c, ...
            'FACENORMALS', TR.faceNormal, 'TOL', tol_inpolyhedron);
        solid_c = permute(solid_c, [2 1 3]);

        if lcheckdist
            fluid_c = ~solid_c;
            [fluid_IB_c, ~] = getBoundaryCells(xgrid_c, ygrid_c, zgrid_c, fluid_c, solid_c, diag_neighbs);
            if (stl_ground)
                fluid_c_1 = fluid_c(:,:,1);
                fluid_IB_c_1 = fluid_IB_c(:,:,1);
                fluid_IB_c_1(fluid_c_1) = true;
                fluid_IB_c(:,:,1) = fluid_IB_c_1;
            end

            [fluid_IB_i_c, fluid_IB_j_c, fluid_IB_k_c] = ind2sub(size(fluid_IB_c), find(fluid_IB_c));
            fluid_IB_xyz_c = [xgrid_c(fluid_IB_i_c)', ygrid_c(fluid_IB_j_c)', zgrid_c(fluid_IB_k_c)'];

            [dists, ~, ~, ~] = point2trimesh('Faces', TR.ConnectivityList, 'Vertices', TR.Points, ...
                'QueryPoints', fluid_IB_xyz_c, 'UseSubSurface', false);

            for n=1:size(fluid_IB_xyz_c,1)
                if abs(dists(n)) < tol_inpolyhedron
                    [~, i] = ismember(fluid_IB_xyz_c(n,1), xgrid_c);
                    [~, j] = ismember(fluid_IB_xyz_c(n,2), ygrid_c);
                    [~, k] = ismember(fluid_IB_xyz_c(n,3), zgrid_c);
                    solid_c(i,j,k) = true;
                end
            end
        end
    end

    % check for acute corners - single column of solid cells
    solid_c_copy = solid_c;
    %solid_c_corners = false(itot,jtot,ktot);

    for k=1:ktot
        for j=1:jtot
            for i=1:itot
                if solid_c(i,j,k)
                    if i ~= 1 && i ~= itot
                        if ~solid_c(i+1,j,k) && ~solid_c(i-1,j,k)
                            solid_c_copy(i,j,k) = false;
                            %solid_c_corners(i,j,k) = true;
                        end
                    end
                    if j ~= 1 && j ~= jtot
                        if ~solid_c(i,j+1,k) && ~solid_c(i,j-1,k)
                            solid_c_copy(i,j,k) = false;
                            %solid_c_corners(i,j,k) = true;
                        end
                    end
                end
            end
        end
    end
    solid_c = solid_c_copy;

    [solid_i_c, solid_j_c, solid_k_c] = ind2sub(size(solid_c), find(solid_c));
    solid_ijk_c = [solid_i_c, solid_j_c, solid_k_c];
    filename_c = [fpath 'solid_c.txt'];
    fileID_c = fopen(filename_c,'W');
    fprintf(fileID_c, '# position (i,j,k)\n');
    fclose(fileID_c);
    dlmwrite(filename_c, solid_ijk_c, '-append','delimiter',' ');
    disp('Written solid_c.txt')
    fluid_c = ~solid_c;

    [fluid_IB_c, solid_IB_c] = getBoundaryCells(xgrid_c, ygrid_c, zgrid_c, fluid_c, solid_c, diag_neighbs);
    if (stl_ground)
        fluid_c_1 = fluid_c(:,:,1);
        fluid_IB_c_1 = fluid_IB_c(:,:,1);
        fluid_IB_c_1(fluid_c_1) = true;
        fluid_IB_c(:,:,1) = fluid_IB_c_1;
    end

    [fluid_IB_i_c, fluid_IB_j_c, fluid_IB_k_c] = ind2sub(size(fluid_IB_c), find(fluid_IB_c));
    fluid_IB_xyz_c = [xgrid_c(fluid_IB_i_c)', ygrid_c(fluid_IB_j_c)', zgrid_c(fluid_IB_k_c)'];

    [solid_IB_i_c, solid_IB_j_c, solid_IB_k_c] = ind2sub(size(solid_IB_c), find(solid_IB_c));
    solid_IB_xyz_c = [xgrid_c(solid_IB_i_c)', ygrid_c(solid_IB_j_c)', zgrid_c(solid_IB_k_c)'];

    fluid_IB_ijk_c = [fluid_IB_i_c, fluid_IB_j_c, fluid_IB_k_c];
    solid_IB_ijk_c = [solid_IB_i_c, solid_IB_j_c, solid_IB_k_c];

    fluid_boundary_c = fluid_IB_ijk_c;
    filename_c = [fpath 'fluid_boundary_c.txt'];
    fileID_c = fopen(filename_c,'W');
    fprintf(fileID_c, '# position (i,j,k)\n');
    fclose(fileID_c);
    dlmwrite(filename_c, fluid_boundary_c, '-append','delimiter',' ');
    disp('Written fluid_boundary_c.txt')

    solid_boundary_c = solid_IB_ijk_c;
    filename_c = [fpath 'solid_boundary_c.txt'];
    fileID_c = fopen(filename_c,'W');
    fprintf(fileID_c, '# position (i,j,k)\n');
    fclose(fileID_c);
    dlmwrite(filename_c, solid_boundary_c, '-append','delimiter',' ');
    disp('Written solid_boundary_c.txt')

end

if lmatchFacetsToCellsFortran
    disp('Determining facet sections using Fortran.')
%     in_mypoly_fortran_path = [folder '/in_mypoly_fortran/'];
%     addpath(in_mypoly_fortran_path)
    cd(folder);
    % Needs fluid_IB_u.txt and fluid_boundary_u.txt to be defined.
    system('gfortran -O2 matchFacetsToCells.f90 -o MFTC.exe');
    copyfile('MFTC.exe', fpath)
    delete MFTC.exe

    fileID = fopen([fpath 'info_matchFacetsToCells.txt'],'w');
    fprintf(fileID,'%15.10f %15.10f\n',[dx dy]');
    fprintf(fileID,'%5d %5d %5d\n',[itot jtot ktot]');
    fprintf(fileID,'%8d %8d\n',[size(TR.ConnectivityList, 1), size(TR.Points, 1)]);
    fprintf(fileID,'%8d %8d %8d %8d\n',[size(fluid_IB_ijk_u,1) size(fluid_IB_ijk_v,1) size(fluid_IB_ijk_w,1) size(fluid_IB_ijk_c,1)]);
    fprintf(fileID,'%8d %8d %8d %8d\n',[size(solid_IB_ijk_u,1) size(solid_IB_ijk_v,1) size(solid_IB_ijk_w,1) size(solid_IB_ijk_c,1)]);
    fprintf(fileID,'%d %d %d\n',[periodic_x, periodic_y, diag_neighbs]);
    fclose(fileID);

    fileID = fopen([fpath 'zhgrid.txt'],'w');
    fprintf(fileID,'%15.10f\n',zgrid_w');
    fclose(fileID);

    fileID = fopen([fpath 'zfgrid.txt'],'w');
    fprintf(fileID,'%15.10f\n',zgrid_c');
    fclose(fileID);

    fileID = fopen([fpath 'vertices.txt'],'w');
    fprintf(fileID,'%15.10f %15.10f %15.10f\n',TR.Points');
    fclose(fileID);

    fileID = fopen([fpath 'faces.txt'],'w');
    fprintf(fileID,'%8d %8d %8d %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n',[TR.ConnectivityList TR.incenter TR.faceNormal]');
    fclose(fileID);

    cd(fpath)
    if lwindows
        system('MFTC.exe');
    else
        system('./MFTC.exe');
    end
    delete MFTC.exe info_matchFacetsToCells.txt faces.txt vertices.txt zfgrid.txt zhgrid.txt;
    cd(currentPath)

    facet_sections_u_fromfile = readmatrix([fpath 'facet_sections_u.txt'],'Range', 2);
    facet_sections_u = NaN(size(facet_sections_u_fromfile,1), 9);
    facet_sections_u(:,[1,2,5,6]) = facet_sections_u_fromfile;

    facet_sections_v_fromfile = readmatrix([fpath 'facet_sections_v.txt'],'Range', 2);
    facet_sections_v = NaN(size(facet_sections_v_fromfile,1), 9);
    facet_sections_v(:,[1,2,5,6]) = facet_sections_v_fromfile;

    facet_sections_w_fromfile = readmatrix([fpath 'facet_sections_w.txt'],'Range', 2);
    facet_sections_w = NaN(size(facet_sections_w_fromfile,1), 9);
    facet_sections_w(:,[1,2,5,6]) = facet_sections_w_fromfile;

    facet_sections_c_fromfile = readmatrix([fpath 'facet_sections_c.txt'],'Range', 2);
    facet_sections_c = NaN(size(facet_sections_c_fromfile,1), 9);
    facet_sections_c(:,[1,2,5,6]) = facet_sections_c_fromfile;

else
    %% Facet sections
    disp('Determining facet sections using MATLAB.')
    if (calculate_facet_sections_uvw)
        % u
        facet_sections_u = matchFacetsToCells(...
            TR, fluid_IB_u, solid_IB_u, fluid_IB_xyz_u, solid_IB_xyz_u, xgrid_u, ygrid_u, zgrid_u, diag_neighbs, periodic_x, periodic_y);

        % For debugging - area 'visible' to fluid for each facet.
        area_facets_u = zeros(nfcts,1);
        for n=1:nfcts
            area_facets_u(n) = sum(facet_sections_u(facet_sections_u(:,1)==n, 2));
        end
        area_fluid_IB_u = sum(area_facets_u);

        filename_u = [fpath 'facet_sections_u.txt'];
        fileID_u = fopen(filename_u,'W');
        fprintf(fileID_u, facet_sections_header);
        fprintf(fileID_u, facet_sections_format, facet_sections_u(:,[1,2,5,6])');
        fclose(fileID_u);
        disp('Written facet_sections_u.txt')

        % v
        facet_sections_v = matchFacetsToCells(...
            TR, fluid_IB_v, solid_IB_v, fluid_IB_xyz_v, solid_IB_xyz_v, xgrid_v, ygrid_v, zgrid_v, diag_neighbs, periodic_x, periodic_y);
        area_facets_v = zeros(nfcts,1);
        for n=1:nfcts
            area_facets_v(n) = sum(facet_sections_v(facet_sections_v(:,1)==n, 2));
        end
        area_fluid_IB_v = sum(area_facets_v);

        filename_v = [fpath 'facet_sections_v.txt'];
        fileID_v = fopen(filename_v,'W');
        fprintf(fileID_v, facet_sections_header);
        fprintf(fileID_v, facet_sections_format, facet_sections_v(:,[1,2,5,6])');
        fclose(fileID_v);
        disp('Written facet_sections_v.txt')

        % w
        facet_sections_w = matchFacetsToCells(...
            TR, fluid_IB_w, solid_IB_w, fluid_IB_xyz_w, solid_IB_xyz_w, xgrid_w, ygrid_w, zgrid_w, diag_neighbs, periodic_x, periodic_y);
        area_facets_w = zeros(nfcts,1);
        for n=1:nfcts
            area_facets_w(n) = sum(facet_sections_w(facet_sections_w(:,1)==n, 2));
        end
        area_fluid_IB_w = sum(area_facets_w);

        filename_w = [fpath 'facet_sections_w.txt'];
        fileID_w = fopen(filename_w,'W');
        fprintf(fileID_w, facet_sections_header);
        fprintf(fileID_w, facet_sections_format, facet_sections_w(:,[1,2,5,6])');
        fclose(fileID_w);
        disp('Written facet_sections_w.txt')
    else
        facet_sections_u = [];
        facet_sections_v = [];
        facet_sections_w = [];
    end

    % c
    if (calculate_facet_sections_c)
        facet_sections_c = matchFacetsToCells(...
            TR, fluid_IB_c, solid_IB_c, fluid_IB_xyz_c, solid_IB_xyz_c, xgrid_c, ygrid_c, zgrid_c, diag_neighbs, periodic_x, periodic_y);
        area_facets_c = zeros(nfcts,1);
        for n=1:nfcts
            area_facets_c(n) = sum(facet_sections_c(facet_sections_c(:,1)==n, 2));
        end
        area_fluid_IB_c = sum(area_facets_c);

        filename_c = [fpath 'facet_sections_c.txt'];
        fileID_c = fopen(filename_c,'W');
        fprintf(fileID_c, facet_sections_header);
        fprintf(fileID_c, facet_sections_format, facet_sections_c(:,[1,2,5,6])');
        fclose(fileID_c);
        disp('Written facet_sections_c.txt')
    else
        facet_sections_c = [];
    end
end

lBImin = false;
% Instead of using distance of boundary point to facet section, use distance of
% boundary point to ALL facets
if lBImin
    facet_sections_u_2 = facet_sections_u;
    [alldists, allBIs, ~, ~] = point2trimesh('Faces', TR.ConnectivityList, 'Vertices', TR.Points, 'QueryPoints', fluid_IB_xyz_u, 'UseSubSurface', false);
    for n=1:size(facet_sections_u,1)
        m = facet_sections_u(n,5);
        facet_sections_u_2(n,6) = alldists(m);
        facet_sections_u_2(n,7:9) = allBIs(m,:);
    end
    filename_u = [fpath 'facet_sections_u_2.txt'];
    fileID_u = fopen(filename_u,'W');
    fprintf(fileID_u, facet_sections_header);
    fprintf(fileID_u, facet_sections_format, facet_sections_u_2(:,[1,2,5,6])');
    fclose(fileID_u);
    disp('Written facet_sections_u_2.txt')

    facet_sections_v_2 = facet_sections_v;
    [alldists, allBIs, ~, ~] = point2trimesh('Faces', TR.ConnectivityList, 'Vertices', TR.Points, 'QueryPoints', fluid_IB_xyz_v, 'UseSubSurface', false);
    for n=1:size(facet_sections_v,1)
        m = facet_sections_v(n,5);
        facet_sections_v_2(n,6) = alldists(m);
        facet_sections_v_2(n,7:9) = allBIs(m,:);
    end
    filename_v = [fpath 'facet_sections_v_2.txt'];
    fileID_v = fopen(filename_v,'W');
    fprintf(fileID_v, facet_sections_header);
    fprintf(fileID_v, facet_sections_format, facet_sections_v_2(:,[1,2,5,6])');
    fclose(fileID_v);
    disp('Written facet_sections_v_2.txt')

    facet_sections_w_2 = facet_sections_w;
    [alldists, allBIs, ~, ~] = point2trimesh('Faces', TR.ConnectivityList, 'Vertices', TR.Points, 'QueryPoints', fluid_IB_xyz_w, 'UseSubSurface', false);
    for n=1:size(facet_sections_w,1)
        m = facet_sections_w(n,5);
        facet_sections_w_2(n,6) = alldists(m);
        facet_sections_w_2(n,7:9) = allBIs(m,:);
    end
    filename_w = [fpath 'facet_sections_w_2.txt'];
    fileID_w = fopen(filename_w,'W');
    fprintf(fileID_w, facet_sections_header);
    fprintf(fileID_w, facet_sections_format, facet_sections_w_2(:,[1,2,5,6])');
    fclose(fileID_w);
    disp('Written facet_sections_w_2.txt')

    facet_sections_c_2 = facet_sections_c;
    [alldists, allBIs, ~, ~] = point2trimesh('Faces', TR.ConnectivityList, 'Vertices', TR.Points, 'QueryPoints', fluid_IB_xyz_c, 'UseSubSurface', false);
    for n=1:size(facet_sections_c,1)
        m = facet_sections_c(n,5);
        facet_sections_c_2(n,6) = alldists(m);
        facet_sections_c_2(n,7:9) = allBIs(m,:);
    end
    filename_c = [fpath 'facet_sections_c_2.txt'];
    fileID_c = fopen(filename_c,'W');
    fprintf(fileID_c, facet_sections_header);
    fprintf(fileID_c, facet_sections_format, facet_sections_c_2(:,[1,2,5,6])');
    fclose(fileID_c);
    disp('Written facet_sections_c_2.txt')
end

%%
disp(['nfcts = ', num2str(nfcts)])
disp(['nsolpts_u = ', num2str(size(solid_ijk_u,1))])
disp(['nsolpts_v = ', num2str(size(solid_ijk_v,1))])
disp(['nsolpts_w = ', num2str(size(solid_ijk_w,1))])
disp(['nsolpts_c = ', num2str(size(solid_ijk_c,1))])
disp(['nbndpts_u = ', num2str(size(fluid_IB_xyz_u,1))])
disp(['nbndpts_v = ', num2str(size(fluid_IB_xyz_v,1))])
disp(['nbndpts_w = ', num2str(size(fluid_IB_xyz_w,1))])
disp(['nbndpts_c = ', num2str(size(fluid_IB_xyz_c,1))])
disp(['nfctsecs_u = ', num2str(size(facet_sections_u,1))])
disp(['nfctsecs_v = ', num2str(size(facet_sections_v,1))])
disp(['nfctsecs_w = ', num2str(size(facet_sections_w,1))])
disp(['nfctsecs_c = ', num2str(size(facet_sections_c,1))])

filename_info = [fpath 'info.txt'];
fileID_info = fopen(filename_info,'W');
fprintf(fileID_info, ['nfcts = ', num2str(nfcts), '\n']);
fprintf(fileID_info, ['nsolpts_u = ', num2str(size(solid_ijk_u,1)), '\n']);
fprintf(fileID_info, ['nsolpts_v = ', num2str(size(solid_ijk_v,1)), '\n']);
fprintf(fileID_info, ['nsolpts_w = ', num2str(size(solid_ijk_w,1)), '\n']);
fprintf(fileID_info, ['nsolpts_c = ', num2str(size(solid_ijk_c,1)), '\n']);
fprintf(fileID_info, ['nbndpts_u = ', num2str(size(fluid_IB_xyz_u,1)), '\n']);
fprintf(fileID_info, ['nbndpts_v = ', num2str(size(fluid_IB_xyz_v,1)), '\n']);
fprintf(fileID_info, ['nbndpts_w = ', num2str(size(fluid_IB_xyz_w,1)), '\n']);
fprintf(fileID_info, ['nbndpts_c = ', num2str(size(fluid_IB_xyz_c,1)), '\n']);
fprintf(fileID_info, ['nfctsecs_u = ', num2str(size(facet_sections_u,1)), '\n']);
fprintf(fileID_info, ['nfctsecs_v = ', num2str(size(facet_sections_v,1)), '\n']);
fprintf(fileID_info, ['nfctsecs_w = ', num2str(size(facet_sections_w,1)), '\n']);
fprintf(fileID_info, ['nfctsecs_c = ', num2str(size(facet_sections_c,1)), '\n']);
fclose(fileID_info);

%% Clean up exp directory
if lmypolyfortran
    cd(fpath)
    delete flag_u.txt flag_v.txt flag_w.txt flag_c.txt;
    delete fluid_IB_u.txt fluid_IB_v.txt fluid_IB_w.txt fluid_IB_c.txt;
    delete solid_IB_u.txt solid_IB_v.txt solid_IB_w.txt solid_IB_c.txt;
    delete solid_boundary_u.txt solid_boundary_v.txt solid_boundary_w.txt solid_boundary_c.txt;
    cd(currentPath)
end

%% Plot
figure

patch('Faces', TR.ConnectivityList, 'Vertices', TR.Points, 'FaceColor', ones(3,1)*0.85, 'FaceAlpha', 1)
hold on
incenters = TR.incenter;
faceNormals = TR.faceNormal;
%quiver3(incenters(:,1), incenters(:,2), incenters(:,3), faceNormals(:,1), faceNormals(:,2), faceNormals(:,3), 0)
view(3)

axis equal tight

% xlim([0 xsize])
% ylim([0 ysize])
% zlim([0 zsize])

% Uncomment to view solid points
%if isempty(solid_IB_ijk_u)
%    scatter3(xgrid_u([]),ygrid_u([]),zgrid_u([]),10,[0,0,1],'filled')
%else
%    scatter3(xgrid_u(solid_ijk_u(:,1)),ygrid_u(solid_ijk_u(:,2)),zgrid_u(solid_ijk_u(:,3)),10,[0,0,1],'filled')
%end
%if isempty(solid_IB_ijk_v)
%    scatter3(xgrid_v([]),ygrid_v([]),zgrid_v([]),10,[0,0,1],'filled')
%else
%    scatter3(xgrid_v(solid_ijk_v(:,1)),ygrid_v(solid_ijk_v(:,2)),zgrid_v(solid_ijk_v(:,3)),10,[0,0,1],'filled')
%end
%if isempty(solid_IB_ijk_w)
%    scatter3(xgrid_w([]),ygrid_w([]),zgrid_w([]),10,[0,0,1],'filled')
%else
%    scatter3(xgrid_w(solid_ijk_w(:,1)),ygrid_w(solid_ijk_w(:,2)),zgrid_w(solid_ijk_w(:,3)),10,[0,0,1],'filled')
%end
%if isempty(solid_IB_ijk_c)
%    scatter3(xgrid_c([]),ygrid_c([]),zgrid_c([]),10,[0,0,1],'filled')
%else
%    scatter3(xgrid_c(solid_ijk_c(:,1)),ygrid_c(solid_ijk_c(:,2)),zgrid_c(solid_ijk_c(:,3)),10,[0,0,1],'filled')
%end

% Uncomment to view fluid boundary points
%scatter3(fluid_IB_xyz_u(:,1),fluid_IB_xyz_u(:,2),fluid_IB_xyz_u(:,3),10,[0,0,1],'filled')
%scatter3(fluid_IB_xyz_v(:,1),fluid_IB_xyz_v(:,2),fluid_IB_xyz_v(:,3),10,[0,0,1],'filled')
%scatter3(fluid_IB_xyz_w(:,1),fluid_IB_xyz_w(:,2),fluid_IB_xyz_w(:,3),10,[0,0,1],'filled')
%scatter3(fluid_IB_xyz_c(:,1),fluid_IB_xyz_c(:,2),fluid_IB_xyz_c(:,3),10,[0,0,1],'filled')
