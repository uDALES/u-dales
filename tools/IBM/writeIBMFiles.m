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
%% Determine total area of surface

if lmypolyfortran
    write_pre_info;
    if lwindows
        in_mypoly_fortran_path = [DA_TOOLSDIR '/IBM/in_mypoly_fortran/'];
        cd(in_mypoly_fortran_path)
        system('gfortran -O2 -fopenmp in_mypoly_functions.f90 IBM_flagging.f90 -o pre.exe');
        copyfile('pre.exe',fpath); % remember to build pre.exe in local system. gfortran -O2 -fopenmp in_mypoly_functions.f90 IBM_flagging.f90 -o pre.exe
        delete pre.exe in_mypoly_functions.mod;
        cd(fpath)
        system('pre.exe'); 
        delete pre.exe inmypoly_inp_info.txt Stl_data.txt vertices.txt zfgrid.txt zhgrid.txt;
    else
        cd(DA_EXPDIR)
        cd ..
        in_mypoly_command = ['u-dales/tools/IBM/in_mypoly_fortran/pre_run.sh experiments/' expnr];
        system(in_mypoly_command);
        cd(fpath)
    end
else
    if lmypoly
        max_height = max(TR.Points(:,3)) + tol_mypoly;
        L_char = max_facet_size(TR.Points,TR.ConnectivityList) + tol_mypoly;
    end
end

%% Calculate u
disp('Determining solid points for u-grid.')
if lmypolyfortran
    formatSpec = '%d';
    fileID = fopen([fpath 'flag_u.txt'],'r');
    flag_u = fscanf(fileID,formatSpec);
    fclose(fileID);
    count = 1;
    for iy = 1:r.jtot
        for iz = 1:r.ktot
            for ix = 1:r.itot
                if (flag_u(count) == 1)
                            solid_u(ix,iy,iz) = true;
                else
                            solid_u(ix,iy,iz) = false;
                end
                count = count+1;
            end
        end
    end
    solid_ijk_u = readmatrix('solid_u.txt');
    disp('Written solid_u.txt')
else
    if lmypoly
        solid_u = in_grid_mypoly(TR.Points,TR.ConnectivityList, ...
            TR.incenter,TR.faceNormal,xgrid_u,ygrid_u,zgrid_u,Dir_ray_u,L_char,max_height,tol_mypoly);
    else
        solid_u = inpolyhedron(TR.ConnectivityList, TR.Points, ...
            xgrid_u, ygrid_u, zgrid_u, 'FACENORMALS', TR.faceNormal);
        solid_u = permute(solid_u, [2 1 3]);
        %solid_u(1,:,:) = 0;
    end
    
    [solid_i_u, solid_j_u, solid_k_u] = ind2sub(size(solid_u), find(solid_u));
    solid_ijk_u = [solid_i_u, solid_j_u, solid_k_u];
    filename_u = [fpath 'solid_u.txt'];
    fileID_u = fopen(filename_u,'W');
    fprintf(fileID_u, '# position (i,j,k)\n');
    fclose(fileID_u);
    dlmwrite(filename_u, solid_ijk_u, '-append','delimiter',' ');
    disp('Written solid_u.txt')
    
end

fluid_u = ~solid_u;

%% Boundary masks
disp('Determining fluid boundary points for u-grid.')
[fluid_IB_u, solid_IB_u] = getBoundaryCells(xgrid_u, ygrid_u, zgrid_u, fluid_u, solid_u, diag_neighbs);
if (stl_ground)
    fluid_u_1 = fluid_u(:,:,1);
    fluid_IB_u_1 = fluid_IB_u(:,:,1);
    fluid_IB_u_1(fluid_u_1) = true;
    fluid_IB_u(:,:,1) = fluid_IB_u_1;
end

% Boundary coordinates
[fluid_IB_i_u, fluid_IB_j_u, fluid_IB_k_u] = ind2sub(size(fluid_IB_u), find(fluid_IB_u));
fluid_IB_xyz_u = [xgrid_u(fluid_IB_i_u)', ygrid_u(fluid_IB_j_u)', zgrid_u(fluid_IB_k_u)'];

[solid_IB_i_u, solid_IB_j_u, solid_IB_k_u] = ind2sub(size(solid_IB_u), find(solid_IB_u));
solid_IB_xyz_u = [xgrid_u(solid_IB_i_u)', ygrid_u(solid_IB_j_u)', zgrid_u(solid_IB_k_u)'];

% Fluid boundary points
fluid_IB_ijk_u = [fluid_IB_i_u, fluid_IB_j_u, fluid_IB_k_u];
fluid_boundary_u = fluid_IB_ijk_u;
filename_u = [fpath 'fluid_boundary_u.txt'];
fileID_u = fopen(filename_u,'W');
fprintf(fileID_u, '# position (i,j,k)\n');
fclose(fileID_u);
dlmwrite(filename_u, fluid_boundary_u, '-append','delimiter',' ');
disp('Written fluid_boundary_u.txt')

%% Facet sections
if (calculate_facet_sections_uvw)
    disp('Determining facet sections for u-grid.')
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
    fprintf(fileID_u, '# facet, area, fluid boundary point, distance\n');
    fprintf(fileID_u, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_u(:,[1,2,5,6])');
    fclose(fileID_u);
    disp('Written facet_sections_u.txt')

    lBImin_u = true;
    % Instead of using distance of boundary point to facet section, use distance of
    % boundary point to ALL facets
    if lBImin_u
        facet_sections_u_2 = facet_sections_u;
        [alldists, allBIs, ~, ~] = point2trimesh('Faces', TR.ConnectivityList, 'Vertices', TR.Points, 'QueryPoints', fluid_IB_xyz_u, 'UseSubSurface', false);
        for n=1:size(facet_sections_u,1)
            m = facet_sections_u(n,5);
            facet_sections_u_2(n,6) = alldists(m);
            facet_sections_u_2(n,7:9) = allBIs(m,:);
        end
        filename_u = [fpath 'facet_sections_u_2.txt'];
        fileID_u = fopen(filename_u,'W');
        fprintf(fileID_u, '# facet, area, fluid boundary point, distance\n');
        fprintf(fileID_u, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_u_2(:,[1,2,5,6])');
        fclose(fileID_u);
        disp('Written facet_sections_u_2.txt')

        facet_sections_u_3 = facet_sections_u;
        for n=1:size(facet_sections_u,1)
            facet = facet_sections_u(n,1);
            m = facet_sections_u(n,5); % fluid boundary point
            inputs.faces = TR.ConnectivityList(facet,:);
            inputs.nodes = TR.Points;
            inputs.face_mean_nodes = TR.incenter(facet);
            inputs.face_normals = TR.faceNormal(facet);
            [dist, BI, ~] = fastPoint2TriMesh(inputs, fluid_boundary_u(m,:), 0, 0);
            facet_sections_u_3(n,6) = dist;
            facet_sections_u_3(n,7:9) = BI;
        end

        filename_u = [fpath 'facet_sections_u_3.txt'];
        fileID_u = fopen(filename_u,'W');
        fprintf(fileID_u, '# facet, area, fluid boundary point, distance\n');
        fprintf(fileID_u, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_u_3(:,[1,2,5,6])');
        fclose(fileID_u);
        disp('Written facet_sections_u_3.txt')

    end
else
    facet_sections_u = [];
end

%% Calculate v
disp('Determining solid points for v-grid.')
if lmypolyfortran
    formatSpec = '%d';
    fileID = fopen([fpath 'flag_v.txt'],'r');
    flag_v = fscanf(fileID,formatSpec);
    fclose(fileID);
    count = 1;
    for iy = 1:r.jtot
        for iz = 1:r.ktot
            for ix = 1:r.itot
                if (flag_v(count) == 1)
                            solid_v(ix,iy,iz) = true;
                else
                            solid_v(ix,iy,iz) = false;
                end
                count = count+1;
            end
        end
    end
    solid_ijk_v = readmatrix('solid_v.txt');
    disp('Written solid_v.txt')
else
    if lmypoly
        solid_v = in_grid_mypoly(TR.Points,TR.ConnectivityList,TR.incenter,TR.faceNormal,xgrid_v,ygrid_v,zgrid_v,Dir_ray_v,L_char,max_height,tol_mypoly);
    else
        solid_v = inpolyhedron(TR.ConnectivityList, TR.Points, ...
            xgrid_v, ygrid_v, zgrid_v, 'FACENORMALS', TR.faceNormal);
        solid_v = permute(solid_v, [2 1 3]);
        %solid_v(:,1,:) = 0;
    end
    
    [solid_i_v, solid_j_v, solid_k_v] = ind2sub(size(solid_v), find(solid_v));
    solid_ijk_v = [solid_i_v, solid_j_v, solid_k_v];
    filename_v = [fpath 'solid_v.txt'];
    fileID_v = fopen(filename_v,'W');
    fprintf(fileID_v, '# position (i,j,k)\n');
    fclose(fileID_v);
    dlmwrite(filename_v, solid_ijk_v, '-append','delimiter',' ');
    disp('Written solid_v.txt')
    
end

fluid_v = ~solid_v;

%% Boundary masks
disp('Determining fluid boundary points for v-grid.')
% Boundary masks
[fluid_IB_v, solid_IB_v] = getBoundaryCells(xgrid_v, ygrid_v, zgrid_v, fluid_v, solid_v, diag_neighbs);
if (stl_ground)
    fluid_v_1 = fluid_v(:,:,1);
    fluid_IB_v_1 = fluid_IB_v(:,:,1);
    fluid_IB_v_1(fluid_v_1) = true;
    fluid_IB_v(:,:,1) = fluid_IB_v_1;
end

% Boundary coordinates
[fluid_IB_i_v, fluid_IB_j_v, fluid_IB_k_v] = ind2sub(size(fluid_IB_v), find(fluid_IB_v));
fluid_IB_xyz_v = [xgrid_v(fluid_IB_i_v)', ygrid_v(fluid_IB_j_v)', zgrid_v(fluid_IB_k_v)'];

[solid_IB_i_v, solid_IB_j_v, solid_IB_k_v] = ind2sub(size(solid_IB_v), find(solid_IB_v));
solid_IB_xyz_v = [xgrid_v(solid_IB_i_v)', ygrid_v(solid_IB_j_v)', zgrid_v(solid_IB_k_v)'];

fluid_IB_ijk_v = [fluid_IB_i_v, fluid_IB_j_v, fluid_IB_k_v];
fluid_boundary_v = fluid_IB_ijk_v;
filename_v = [fpath 'fluid_boundary_v.txt'];
fileID_v = fopen(filename_v,'W');
fprintf(fileID_v, '# position (i,j,k)\n');
fclose(fileID_v);
dlmwrite(filename_v, fluid_boundary_v, '-append','delimiter',' ');
disp('Written fluid_boundary_v.txt')

%% Facet sections
if (calculate_facet_sections_uvw)
    disp('Determining facet sections for v-grid.')
    facet_sections_v = matchFacetsToCells(...
        TR, fluid_IB_v, solid_IB_v, fluid_IB_xyz_v, solid_IB_xyz_v, xgrid_v, ygrid_v, zgrid_v, diag_neighbs, periodic_x, periodic_y);

    area_facets_v = zeros(nfcts,1);
    for n=1:nfcts
        area_facets_v(n) = sum(facet_sections_v(facet_sections_v(:,1)==n, 2));
    end
    area_fluid_IB_v = sum(area_facets_v);

    filename_v = [fpath 'facet_sections_v.txt'];
    fileID_v = fopen(filename_v,'W');
    fprintf(fileID_v, '# facet, area, fluid boundary point, distance\n');
    fprintf(fileID_v, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_v(:,[1,2,5,6])');
    fclose(fileID_v);
    disp('Written facet_sections_v.txt')

    lBImin_v = true;
    % Instead of using distance of boundary point to facet section, use distance of
    % boundary point to ALL facets
    if lBImin_v
        facet_sections_v_2 = facet_sections_v;
        [alldists, allBIs, ~, ~] = point2trimesh('Faces', TR.ConnectivityList, 'Vertices', TR.Points, 'QueryPoints', fluid_IB_xyz_v, 'UseSubSurface', false);
        for n=1:size(facet_sections_v,1)
            m = facet_sections_v(n,5);
            facet_sections_v_2(n,6) = alldists(m);
            facet_sections_v_2(n,7:9) = allBIs(m,:);
        end
        filename_v = [fpath 'facet_sections_v_2.txt'];
        fileID_v = fopen(filename_v,'W');
        fprintf(fileID_v, '# facet, area, fluid boundary point, distance\n');
        fprintf(fileID_v, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_v_2(:,[1,2,5,6])');
        fclose(fileID_v);
        disp('Written facet_sections_v_2.txt')

        facet_sections_v_3 = facet_sections_v;
        for n=1:size(facet_sections_v,1)
            facet = facet_sections_v(n,1);
            m = facet_sections_v(n,5); % fluid boundary point
            inputs.faces = TR.ConnectivityList(facet,:);
            inputs.nodes = TR.Points;
            inputs.face_mean_nodes = TR.incenter(facet);
            inputs.face_normals = TR.faceNormal(facet);
            [dist, BI, ~] = fastPoint2TriMesh(inputs, fluid_boundary_v(m,:), 0, 0);
            facet_sections_v_3(n,6) = dist;
            facet_sections_v_3(n,7:9) = BI;
        end

        filename_v = [fpath 'facet_sections_v_3.txt'];
        fileID_v = fopen(filename_v,'W');
        fprintf(fileID_v, '# facet, area, fluid boundary point, distance\n');
        fprintf(fileID_v, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_v_3(:,[1,2,5,6])');
        fclose(fileID_v);
        disp('Written facet_sections_v_3.txt')
    end
else
    facet_sections_v = [];
end

%% Calculate w
disp('Determining solid points for w-grid.')
if lmypolyfortran
    formatSpec = '%d';
    fileID = fopen([fpath 'flag_w.txt'],'r');
    flag_w = fscanf(fileID,formatSpec);
    fclose(fileID);
    count = 1;
    for iy = 1:r.jtot
        for iz = 1:r.ktot
            for ix = 1:r.itot
                if (flag_w(count) == 1)
                            solid_w(ix,iy,iz) = true;
                else
                            solid_w(ix,iy,iz) = false;
                end
                count = count+1;
            end
        end
    end
    solid_ijk_w = readmatrix('solid_w.txt');
    disp('Written solid_w.txt')
else
    if lmypoly
        solid_w = in_grid_mypoly(TR.Points,TR.ConnectivityList, ...
            TR.incenter,TR.faceNormal,xgrid_w,ygrid_w,zgrid_w,Dir_ray_w,L_char,max_height,tol_mypoly);
    else
        solid_w = inpolyhedron(TR.ConnectivityList, TR.Points, ...
           xgrid_w, ygrid_w, zgrid_w, 'FACENORMALS', TR.faceNormal);
        solid_w = permute(solid_w, [2 1 3]);
    end
    
    [solid_i_w, solid_j_w, solid_k_w] = ind2sub(size(solid_w), find(solid_w));
    solid_ijk_w = [solid_i_w, solid_j_w, solid_k_w];
    filename_w = [fpath 'solid_w.txt'];
    fileID_w = fopen(filename_w,'W');
    fprintf(fileID_w, '# position (i,j,k)\n');
    fclose(fileID_w);
    dlmwrite(filename_w, solid_ijk_w, '-append','delimiter',' ');
    disp('Written solid_w.txt')
    
end

if (stl_ground)
    solid_w_b = solid_w;
    solid_w_b(:,:,1) = true;
end

fluid_w = ~solid_w;

%% Boundary points
disp('Determining fluid boundary points for w-grid.')
% Boundary masks
[fluid_IB_w, solid_IB_w] = getBoundaryCells(xgrid_w, ygrid_w, zgrid_w, fluid_w, solid_w_b, diag_neighbs);
fluid_IB_w(:,:,1) = 0; % Bottom is always solid

% Boundary coordinates
[fluid_IB_i_w, fluid_IB_j_w, fluid_IB_k_w] = ind2sub(size(fluid_IB_w), find(fluid_IB_w));
fluid_IB_xyz_w = [xgrid_w(fluid_IB_i_w)', ygrid_w(fluid_IB_j_w)', zgrid_w(fluid_IB_k_w)'];

[solid_IB_i_w, solid_IB_j_w, solid_IB_k_w] = ind2sub(size(solid_IB_w), find(solid_IB_w));
solid_IB_xyz_w = [xgrid_w(solid_IB_i_w)', ygrid_w(solid_IB_j_w)', zgrid_w(solid_IB_k_w)'];

fluid_IB_ijk_w = [fluid_IB_i_w, fluid_IB_j_w, fluid_IB_k_w];
fluid_boundary_w = fluid_IB_ijk_w;
filename_w = [fpath 'fluid_boundary_w.txt'];
fileID_w = fopen(filename_w,'W');
fprintf(fileID_w, '# position (i,j,k)\n');
fclose(fileID_w);
dlmwrite(filename_w, fluid_boundary_w, '-append','delimiter',' ');
disp('Written fluid_boundary_w.txt')

%% Facet sections
if (calculate_facet_sections_uvw)
    disp('Determining facet sections for w-grid.')
    facet_sections_w = matchFacetsToCells(...
        TR, fluid_IB_w, solid_IB_w, fluid_IB_xyz_w, solid_IB_xyz_w, xgrid_w, ygrid_w, zgrid_w, diag_neighbs, periodic_x, periodic_y);

    area_facets_w = zeros(nfcts,1);
    for n=1:nfcts
        area_facets_w(n) = sum(facet_sections_w(facet_sections_w(:,1)==n, 2));
    end
    area_fluid_IB_w = sum(area_facets_w);

    filename_w = [fpath 'facet_sections_w.txt'];
    fileID_w = fopen(filename_w,'W');
    fprintf(fileID_w, '# facet, area, fluid boundary point, distance\n');

    fprintf(fileID_w, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_w(:,[1,2,5,6])');
    fclose(fileID_w);
    disp('Written facet_sections_w.txt')

    lBImin_w = true;
    % Instead of using distance of boundary point to facet section, use distance of
    % boundary point to ALL facets
    if lBImin_w
        facet_sections_w_2 = facet_sections_w;
        [alldists, allBIs, ~, ~] = point2trimesh('Faces', TR.ConnectivityList, 'Vertices', TR.Points, 'QueryPoints', fluid_IB_xyz_w, 'UseSubSurface', false);
        for n=1:size(facet_sections_w,1)
            m = facet_sections_w(n,5);
            facet_sections_w_2(n,6) = alldists(m);
            facet_sections_w_2(n,7:9) = allBIs(m,:);
        end
        filename_w = [fpath 'facet_sections_w_2.txt'];
        fileID_w = fopen(filename_w,'W');
        fprintf(fileID_w, '# facet, area, fluid boundary point, distance\n');
        fprintf(fileID_w, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_w_2(:,[1,2,5,6])');
        fclose(fileID_w);
        disp('Written facet_sections_w_2.txt')

        facet_sections_w_3 = facet_sections_w;
        for n=1:size(facet_sections_w,1)
            facet = facet_sections_w(n,1);
            m = facet_sections_w(n,5); % fluid boundary point
            inputs.faces = TR.ConnectivityList(facet,:);
            inputs.nodes = TR.Points;
            inputs.face_mean_nodes = TR.incenter(facet);
            inputs.face_normals = TR.faceNormal(facet);
            [dist, BI, ~] = fastPoint2TriMesh(inputs, fluid_boundary_w(m,:), 0, 0);
            facet_sections_w_3(n,6) = dist;
            facet_sections_w_3(n,7:9) = BI;
        end

        filename_w = [fpath 'facet_sections_w_3.txt'];
        fileID_w = fopen(filename_w,'W');
        fprintf(fileID_w, '# facet, area, fluid boundary point, distance\n');
        fprintf(fileID_w, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_w_3(:,[1,2,5,6])');
        fclose(fileID_w);
        disp('Written facet_sections_w_3.txt')
    end
else
    facet_sections_w = [];
end

%% Calculate c
disp('Determining solid points for c-grid.')
if lmypolyfortran
    formatSpec = '%d';
    fileID = fopen([fpath 'flag_c.txt'],'r');
    flag_c = fscanf(fileID,formatSpec);
    fclose(fileID);
    count = 1;
    for iy = 1:r.jtot
        for iz = 1:r.ktot
            for ix = 1:r.itot
                if (flag_c(count) == 1)
                            solid_c(ix,iy,iz) = true;
                else
                            solid_c(ix,iy,iz) = false;
                end
                count = count+1;
            end
        end
    end
    solid_ijk_c = readmatrix('solid_c.txt');
    disp('Written solid_c.txt')
else
    if lmypoly
        solid_c = in_grid_mypoly(TR.Points,TR.ConnectivityList,TR.incenter,TR.faceNormal,xgrid_c,ygrid_c,zgrid_c,Dir_ray_c,L_char,max_height,tol_mypoly);
    else
        solid_c = inpolyhedron(TR.ConnectivityList, TR.Points, ...
            xgrid_c, ygrid_c, zgrid_c, 'FACENORMALS', TR.faceNormal);
        solid_c = permute(solid_c, [2 1 3]);
    end
    
    [solid_i_c, solid_j_c, solid_k_c] = ind2sub(size(solid_c), find(solid_c));
    solid_ijk_c = [solid_i_c, solid_j_c, solid_k_c];
    filename_c = [fpath 'solid_c.txt'];
    fileID_c = fopen(filename_c,'W');
    fprintf(fileID_c, '# position (i,j,k)\n');
    fclose(fileID_c);
    dlmwrite(filename_c, solid_ijk_c, '-append','delimiter',' ');
    disp('Written solid_c.txt')
    
end

fluid_c = ~solid_c;

%% Boundary points
disp('Determining fluid boundary points for c-grid.')
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
fluid_boundary_c = fluid_IB_ijk_c;
filename_c = [fpath 'fluid_boundary_c.txt'];
fileID_c = fopen(filename_c,'W');
fprintf(fileID_c, '# position (i,j,k), distance to surface, reconstruction point location\n');
fclose(fileID_c);
dlmwrite(filename_c, fluid_boundary_c, '-append','delimiter',' ');
disp('Written fluid_boundary_c.txt')

%% Facet sections
if (calculate_facet_sections_c)
    disp('Determining facet sections for c-grid.')
    facet_sections_c = matchFacetsToCells(...
        TR, fluid_IB_c, solid_IB_c, fluid_IB_xyz_c, solid_IB_xyz_c, xgrid_c, ygrid_c, zgrid_c, diag_neighbs, periodic_x, periodic_y);

    area_facets_c = zeros(nfcts,1);
    for n=1:nfcts
        area_facets_c(n) = sum(facet_sections_c(facet_sections_c(:,1)==n, 2));
    end
    area_fluid_IB_c = sum(area_facets_c);

    filename_c = [fpath 'facet_sections_c.txt'];
    fileID_c = fopen(filename_c,'W');
    fprintf(fileID_c, '# facet, area, fluid boundary point, distance\n');
    fprintf(fileID_c, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_c(:,[1,2,5,6])');
    fclose(fileID_c);
    disp('Written facet_sections_c.txt')

    lBImin_c = true;
    % Instead of using distance of boundary point to facet section, use distance of
    % boundary point to ALL facets
    if lBImin_c
        facet_sections_c_2 = facet_sections_c;
        [alldists, allBIs, ~, ~] = point2trimesh('Faces', TR.ConnectivityList, 'Vertices', TR.Points, 'QueryPoints', fluid_IB_xyz_c, 'UseSubSurface', false);
        for n=1:size(facet_sections_c,1)
            m = facet_sections_c(n,5);
            facet_sections_c_2(n,6) = alldists(m);
            facet_sections_c_2(n,7:9) = allBIs(m,:);
        end
        filename_c = [fpath 'facet_sections_c_2.txt'];
        fileID_c = fopen(filename_c,'W');
        fprintf(fileID_c, '# facet, area, fluid boundary point, distance\n');
        fprintf(fileID_c, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_c_2(:,[1,2,5,6])');
        fclose(fileID_c);
        disp('Written facet_sections_c_2.txt')

        facet_sections_c_3 = facet_sections_c;
        for n=1:size(facet_sections_c,1)
            facet = facet_sections_c(n,1);
            m = facet_sections_c(n,5); % fluid boundary point
            inputs.faces = TR.ConnectivityList(facet,:);
            inputs.nodes = TR.Points;
            inputs.face_mean_nodes = TR.incenter(facet);
            inputs.face_normals = TR.faceNormal(facet);
            [dist, BI, ~] = fastPoint2TriMesh(inputs, fluid_boundary_c(m,:), 0, 0);
            facet_sections_c_3(n,6) = dist;
            facet_sections_c_3(n,7:9) = BI;
        end

        filename_c = [fpath 'facet_sections_c_3.txt'];
        fileID_c = fopen(filename_c,'W');
        fprintf(fileID_c, '# facet, area, fluid boundary point, distance\n');
        fprintf(fileID_c, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_c_3(:,[1,2,5,6])');
        fclose(fileID_c);
        disp('Written facet_sections_c_3.txt')
    end
else
    facet_sections_c = [];
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
    delete flag_u.txt flag_v.txt flag_w.txt flag_c.txt;
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

xlim([0 xsize])
ylim([0 ysize])
zlim([0 zsize])

scatter3(X_u(solid_u), Y_u(solid_u), Z_u(solid_u), 10,[0,0,1],'filled')
scatter3(X_v(solid_v), Y_v(solid_v), Z_v(solid_v), 10,[0,0,1],'filled')
scatter3(X_w(solid_w), Y_w(solid_w), Z_w(solid_w), 10,[0,0,1],'filled')
scatter3(X_c(solid_c), Y_c(solid_c), Z_c(solid_c), 10,[0,0,1],'filled')

%scatter3(X_u(fluid_IB_u), Y_u(fluid_IB_u), Z_u(fluid_IB_u), 10,[0,0,1],'filled')
%scatter3(X_v(fluid_IB_v), Y_v(fluid_IB_v), Z_v(fluid_IB_v), 10,[0,0,1],'filled')
%scatter3(X_w(fluid_IB_w), Y_u(fluid_IB_w), Z_u(fluid_IB_w), 10,[0,0,1],'filled')
%scatter3(X_c(fluid_IB_c), Y_u(fluid_IB_c), Z_u(fluid_IB_c), 10,[0,0,1],'filled')

%scatter3(fluid_IB_xyz_u(:,1),fluid_IB_xyz_u(:,2),fluid_IB_xyz_u(:,3),10,[0,0,1],'filled')
%scatter3(fluid_IB_xyz_v(:,1),fluid_IB_xyz_v(:,2),fluid_IB_xyz_v(:,3),10,[0,0,1],'filled')
%scatter3(fluid_IB_xyz_w(:,1),fluid_IB_xyz_w(:,2),fluid_IB_xyz_w(:,3),10,[0,0,1],'filled')
%scatter3(fluid_IB_xyz_c(:,1),fluid_IB_xyz_c(:,2),fluid_IB_xyz_c(:,3),10,[0,0,1],'filled')

% %% u
% scatter3(fluid_IB_xyz_u(:,1),fluid_IB_xyz_u(:,2),fluid_IB_xyz_u(:,3),10,[0,0,1],'filled')
% %scatter3(solid_IB_xyz_u(:,1),solid_IB_xyz_u(:,2),solid_IB_xyz_u(:,3),10,[0,0,1],'filled')
% %quiver3(fluid_IB_BI_u(:,1), fluid_IB_BI_u(:,2), fluid_IB_BI_u(:,3), fluid_IB_vec_u(:,1), fluid_IB_vec_u(:,2), fluid_IB_vec_u(:,3),'off')
%scatter3(fluid_IB_rec_u(:,1),fluid_IB_rec_u(:,2),fluid_IB_rec_u(:,3),10,[0,0,1],'filled')
%quiver3(fluid_IB_xyz_u(:,1), fluid_IB_xyz_u(:,2), fluid_IB_xyz_u(:,3), fluid_IB_rec_vec_u(:,1), fluid_IB_rec_vec_u(:,2), fluid_IB_rec_vec_u(:,3),'off')
%quiver3(fluid_IB_xyz_u(:,1), fluid_IB_xyz_u(:,2), fluid_IB_xyz_u(:,3), fluid_IB_rec_vec_u(:,1), fluid_IB_rec_vec_u(:,2), fluid_IB_rec_vec_u(:,3),'off')
%
% %% v
% scatter3(fluid_IB_xyz_v(:,1),fluid_IB_xyz_v(:,2),fluid_IB_xyz_v(:,3),10,[0,0,1],'filled')
% %scatter3(solid_IB_xyz_v(:,1),solid_IB_xyz_v(:,2),solid_IB_xyz_v(:,3),10,[0,0,1],'filled')
% quiver3(fluid_IB_BI_v(:,1), fluid_IB_BI_v(:,2), fluid_IB_BI_v(:,3), fluid_IB_vec_v(:,1), fluid_IB_vec_v(:,2), fluid_IB_vec_v(:,3),'off')
% %scatter3(fluid_IB_rec_v(:,1),fluid_IB_rec_v(:,2),fluid_IB_rec_v(:,3),10,[0,0,1],'filled')
% %quiver3(fluid_IB_xyz_v(:,1), fluid_IB_xyz_v(:,2), fluid_IB_xyz_v(:,3), fluid_IB_rec_vec_v(:,1), fluid_IB_rec_vec_v(:,2), fluid_IB_rec_vec_v(:,3),'off')
%
% %% w
% scatter3(fluid_IB_xyz_w(:,1),fluid_IB_xyz_w(:,2),fluid_IB_xyz_w(:,3),10,[0,0,1],'filled')
% %scatter3(solid_IB_xyz_w(:,1),solid_IB_xyz_w(:,2),solid_IB_xyz_w(:,3),10,[0,0,1],'filled')
% quiver3(fluid_IB_BI_w(:,1), fluid_IB_BI_w(:,2), fluid_IB_BI_w(:,3), fluid_IB_vec_w(:,1), fluid_IB_vec_w(:,2), fluid_IB_vec_w(:,3),'off')
% %scatter3(fluid_IB_rec_w(:,1),fluid_IB_rec_w(:,2),fluid_IB_rec_w(:,3),10,[0,0,1],'filled')
% %quiver3(fluid_IB_xyz_w(:,1), fluid_IB_xyz_w(:,2), fluid_IB_xyz_w(:,3), fluid_IB_rec_vec_w(:,1), fluid_IB_rec_vec_w(:,2), fluid_IB_rec_vec_w(:,3),'off')
%
% %% c
% scatter3(fluid_IB_xyz_c(:,1),fluid_IB_xyz_c(:,2),fluid_IB_xyz_c(:,3),10,[0,0,1],'filled')
% %scatter3(solid_IB_xyz_c(:,1),solid_IB_xyz_c(:,2),solid_IB_xyz_c(:,3),10,[0,0,1],'filled')
% quiver3(fluid_IB_BI_c(:,1), fluid_IB_BI_c(:,2), fluid_IB_BI_c(:,3), fluid_IB_vec_c(:,1), fluid_IB_vec_c(:,2), fluid_IB_vec_c(:,3),'off')
% %scatter3(fluid_IB_rec_c(:,1),fluid_IB_rec_c(:,2),fluid_IB_rec_c(:,3),10,[0,0,1],'filled')
% %quiver3(fluid_IB_xyz_c(:,1), fluid_IB_xyz_c(:,2), fluid_IB_xyz_c(:,3), fluid_IB_rec_vec_c(:,1), fluid_IB_rec_vec_c(:,2), fluid_IB_rec_vec_c(:,3),'off')
