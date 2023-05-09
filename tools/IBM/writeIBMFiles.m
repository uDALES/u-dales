addpath('./inpolyhedron/')
addpath('./point2trimesh/')
addpath('./in_mypoly/')

% Assumes the following variables have been already defined:
% TR: triangulation describing the entire geometry, including ground facets
%   note this shouldn't be closed - the bottom of buildings are not present.
% TR_noground: A CLOSED version of the geometry TR, with the ground facets
%   not present and the bottom of buildings represented by facets.
%   Eventually it should be possible to go between these two.
% xgrid_u, xgrid_v, xgrid_w, xgrid_c: the x coordinates on which u,v,w, 
%   and scalars are defined (and similar arrays for y and z coords).
% fpath: directory to write files to
% lgroundfacets: true if there are facets on ground level (recommended).
 
% It will write out the following files:
% solid_u/v/w: list of indices of solid points (inside the geometry).
%   Calculated using inpolyhedron if lmypoly=false, or using Dipanjan's routine
%   if lmypoly=true, which requires a closed triangulation (DON'T THINK
%   THIS IS TRUE ANY MORE).
% fluid_boundary_u/v/w: list of indices of fluid points that have solid neighbours.
%   `include_diagonals` determines whether non-directly-adjacent cells are 'neighbours'.
% facet_sections_u/v/w: list of facet section information: facet id,
%   section area, fluid boundary point id, distance to section/surface as a whole (_2)

% Note for now a tolerance of 1e-10 is used here and there to avoid precision errors.
nfcts = size(TR,1);
%% Determine total area of surface
tol = 1e-10;
%tol = 0;
area_facets = zeros(size(TR,1),1);
for facet=1:size(TR.ConnectivityList,1)
    verts = TR.Points(TR.ConnectivityList(facet,:), :);
    nml = cross(verts(2,:)-verts(1,:),verts(3,:)-verts(1,:));
    mag = norm(nml);
    % Only include facets that are not at z=0 AND whose normal is in -z dir
%     if ~all(abs(nml/mag - [0 0 -1]) < tol) && ~all(verts(:,3) < tol)
         area_facets(facet) = 1/2*mag;
%     end
end

%% Calculate u
if lmypoly
    max_height = max(TR.Points(:,3));
    solid_u = in_grid_mypoly(TR_noground.Points,TR_noground.ConnectivityList, ...
        TR_noground.incenter,TR_noground.faceNormal,xgrid_u,ygrid_u,zgrid_u,L_char,max_height);
else
    solid_u = inpolyhedron(TR_noground.ConnectivityList, TR_noground.Points, ...
        xgrid_u, ygrid_u, zgrid_u, 'FACENORMALS', TR_noground.faceNormal);
    solid_u = permute(solid_u, [2 1 3]);
end

fluid_u = ~solid_u;

%% Boundary masks
[fluid_IB_u, solid_IB_u] = getBoundaryCells(xgrid_u, ygrid_u, zgrid_u, fluid_u, solid_u, include_diagonals);
if (lgroundfacets)
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

% Facet sections
facet_sections_u = matchFacetsToCells(...
    TR, fluid_IB_u, solid_IB_u, fluid_IB_xyz_u, solid_IB_xyz_u, xgrid_u, ygrid_u, zgrid_u, include_diagonals);
%area_fluid_IB_u = sum(facet_sections_u(:,2)); % The total area for each IB cell should equal the sum of the facet

area_facets_u = zeros(nfcts,1);
for n=1:nfcts
    area_facets_u(n) = sum(facet_sections_u(facet_sections_u(:,1)==n, 2));
end
area_fluid_IB_u = sum(area_facets_u);

%% Write u
% Solid points
[solid_i_u, solid_j_u, solid_k_u] = ind2sub(size(solid_u), find(solid_u));
solid_ijk_u = [solid_i_u, solid_j_u, solid_k_u];
filename_u = [fpath 'solid_u.txt'];
fileID_u = fopen(filename_u,'W');
fprintf(fileID_u, '# position (i,j,k)\n');
fclose(fileID_u);
dlmwrite(filename_u, solid_ijk_u, '-append','delimiter',' ');

% Fluid boundary points
fluid_IB_ijk_u = [fluid_IB_i_u, fluid_IB_j_u, fluid_IB_k_u];
fluid_boundary_u = [fluid_IB_ijk_u];%, fluid_IB_rec_u];
filename_u = [fpath 'fluid_boundary_u.txt'];
fileID_u = fopen(filename_u,'W');
fprintf(fileID_u, '# position (i,j,k)\n');
fclose(fileID_u);
dlmwrite(filename_u, fluid_boundary_u, '-append','delimiter',' ');

% Facet sections
filename_u = [fpath 'facet_sections_u.txt'];
fileID_u = fopen(filename_u,'W');
fprintf(fileID_u, '# facet, area, fluid boundary point, boundary intercept\n');
%fprintf(fileID_u, '%-2d %-4.4f %-4d %-4.4f %-4.4f %-4.4f\n', facet_sections_u(:,[1,2,5,7:9])');
fprintf(fileID_u, '%-2d %-4.4f %-4d %-4.4f\n', facet_sections_u(:,[1,2,5,6])');
fclose(fileID_u);

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
    fprintf(fileID_u, '# facet, area, fluid boundary point, boundary intercept\n');
    %fprintf(fileID_u, '%-2d %-4.4f %-4d %-4.4f %-4.4f %-4.4f\n', facet_sections_u_2(:,[1,2,5,7:9])');
    fprintf(fileID_u, '%-2d %-4.4f %-4d %-4.4f\n', facet_sections_u_2(:,[1,2,5,6])');
    fclose(fileID_u);
end

%% Calculate v

if lmypoly
    solid_v = in_grid_mypoly(TR_noground.Points,TR_noground.ConnectivityList,TR_noground.incenter,TR_noground.faceNormal,xgrid_v,ygrid_v,zgrid_v,L_char,max_height);
else
    solid_v = inpolyhedron(TR_noground.ConnectivityList, TR_noground.Points, ...
        xgrid_v, ygrid_v, zgrid_v, 'FACENORMALS', TR_noground.faceNormal);
    solid_v = permute(solid_v, [2 1 3]);
end

fluid_v = ~solid_v;
%%
% Boundary masks
[fluid_IB_v, solid_IB_v] = getBoundaryCells(xgrid_v, ygrid_v, zgrid_v, fluid_v, solid_v, include_diagonals);
if (lgroundfacets)
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

% Facet sections
facet_sections_v = matchFacetsToCells(...
    TR, fluid_IB_v, solid_IB_v, fluid_IB_xyz_v, solid_IB_xyz_v, xgrid_v, ygrid_v, zgrid_v, include_diagonals);
%area_fluid_IB_v = sum(facet_sections_v(:,2)); % The total area for each IB cell should equal the sum of the facet

area_facets_v = zeros(nfcts,1);
for n=1:nfcts
    area_facets_v(n) = sum(facet_sections_v(facet_sections_v(:,1)==n, 2));
end
area_fluid_IB_v = sum(area_facets_v);


%% Write v
% Solid points
[solid_i_v, solid_j_v, solid_k_v] = ind2sub(size(solid_v), find(solid_v));
solid_ijk_v = [solid_i_v, solid_j_v, solid_k_v];
filename_v = [fpath 'solid_v.txt'];
fileID_v = fopen(filename_v,'W');
fprintf(fileID_v, '# position (i,j,k)\n');
fclose(fileID_v);
dlmwrite(filename_v, solid_ijk_v, '-append','delimiter',' ');

% Fluid boundary points
fluid_IB_ijk_v = [fluid_IB_i_v, fluid_IB_j_v, fluid_IB_k_v];
fluid_boundary_v = [fluid_IB_ijk_v];%, fluid_IB_rec_v];
filename_v = [fpath 'fluid_boundary_v.txt'];
fileID_v = fopen(filename_v,'W');
fprintf(fileID_v, '# position (i,j,k)\n');
fclose(fileID_v);
dlmwrite(filename_v, fluid_boundary_v, '-append','delimiter',' ');

% Facet sections
filename_v = [fpath 'facet_sections_v.txt'];
fileID_v = fopen(filename_v,'W');
fprintf(fileID_v, '# facet, area, fluid boundary point, boundary intercept\n');
%fprintf(fileID_v, '%-2d %-4.4f %-4d %-4.4f %-4.4f %-4.4f\n', facet_sections_v(:,[1,2,5,7:9])');
fprintf(fileID_v, '%-2d %-4.4f %-4d %-4.4f\n', facet_sections_v(:,[1,2,5,6])');
fclose(fileID_v);

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
    fprintf(fileID_v, '# facet, area, fluid boundary point, boundary intercept\n');
    %fprintf(fileID_v, '%-2d %-4.4f %-4d %-4.4f %-4.4f %-4.4f\n', facet_sections_v_2(:,[1,2,5,7:9])');
    fprintf(fileID_v, '%-2d %-4.4f %-4d %-4.4f\n', facet_sections_v_2(:,[1,2,5,6])');
    fclose(fileID_v);
end

%% Calculate w
if lmypoly
    solid_w = in_grid_mypoly(TR_noground.Points,TR_noground.ConnectivityList, ...
        TR_noground.incenter,TR_noground.faceNormal,xgrid_w,ygrid_w,zgrid_w,L_char,max_height);
else
    solid_w = inpolyhedron(TR_noground.ConnectivityList, TR_noground.Points, ...
        xgrid_w, ygrid_w, zgrid_w, 'FACENORMALS', TR_noground.faceNormal);
    solid_w = permute(solid_w, [2 1 3]);
end

if (lgroundfacets)
    solid_w_b = solid_w;
    solid_w_b(:,:,1) = true;
end

fluid_w = ~solid_w;

%%
% Boundary masks
[fluid_IB_w, solid_IB_w] = getBoundaryCells(xgrid_w, ygrid_w, zgrid_w, fluid_w, solid_w_b, include_diagonals);
fluid_IB_w(:,:,1) = 0; % Bottom layer is solid

% Boundary coordinates
[fluid_IB_i_w, fluid_IB_j_w, fluid_IB_k_w] = ind2sub(size(fluid_IB_w), find(fluid_IB_w));
fluid_IB_xyz_w = [xgrid_w(fluid_IB_i_w)', ygrid_w(fluid_IB_j_w)', zgrid_w(fluid_IB_k_w)'];

[solid_IB_i_w, solid_IB_j_w, solid_IB_k_w] = ind2sub(size(solid_IB_w), find(solid_IB_w));
solid_IB_xyz_w = [xgrid_w(solid_IB_i_w)', ygrid_w(solid_IB_j_w)', zgrid_w(solid_IB_k_w)'];

%%
% Facet sections
facet_sections_w = matchFacetsToCells(...
    TR, fluid_IB_w, solid_IB_w, fluid_IB_xyz_w, solid_IB_xyz_w, xgrid_w, ygrid_w, zgrid_w, include_diagonals);
area_fluid_IB_w = sum(facet_sections_w(:,2)); % The total area for each IB cell should equal the sum of the facet

area_facets_w = zeros(nfcts,1);
for n=1:nfcts
    area_facets_w(n) = sum(facet_sections_w(facet_sections_w(:,1)==n, 2));
end
area_fluid_IB_w = sum(area_facets_w);

%% Write w
% Solid points
[solid_i_w, solid_j_w, solid_k_w] = ind2sub(size(solid_w), find(solid_w));
solid_ijk_w = [solid_i_w, solid_j_w, solid_k_w];
filename_w = [fpath 'solid_w.txt'];
fileID_w = fopen(filename_w,'W');
fprintf(fileID_w, '# position (i,j,k)\n');
fclose(fileID_w);
dlmwrite(filename_w, solid_ijk_w, '-append','delimiter',' ');

% Fluid boundary points
fluid_IB_ijk_w = [fluid_IB_i_w, fluid_IB_j_w, fluid_IB_k_w];
fluid_boundary_w = [fluid_IB_ijk_w];%, fluid_IB_rec_w];
filename_w = [fpath 'fluid_boundary_w.txt'];
fileID_w = fopen(filename_w,'W');
fprintf(fileID_w, '# position (i,j,k)\n');
fclose(fileID_w);
dlmwrite(filename_w, fluid_boundary_w, '-append','delimiter',' ');

% Facet sections
filename_w = [fpath 'facet_sections_w.txt'];
fileID_w = fopen(filename_w,'W');
fprintf(fileID_w, '# facet, area, fluid boundary point, boundary intercept\n');
%fprintf(fileID_w, '%-2d %-4.4f %-4d %-4.4f %-4.4f %-4.4f\n', facet_sections_w(:,[1,2,5,7:9])');
fprintf(fileID_w, '%-2d %-4.4f %-4d %-4.4f\n', facet_sections_w(:,[1,2,5,6])');
fclose(fileID_w);
%dlmwrite(filename_w, facet_sections_w(:,[1,2,5,7:9]), '-append','delimiter',' ','precision',8);

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
    fprintf(fileID_w, '# facet, area, fluid boundary point, boundary intercept\n');
    %fprintf(fileID_w, '%-2d %-4.4f %-4d %-4.4f %-4.4f %-4.4f\n', facet_sections_w_2(:,[1,2,5,7:9])');
    fprintf(fileID_w, '%-2d %-4.4f %-4d %-4.4f\n', facet_sections_w_2(:,[1,2,5,6])');
    fclose(fileID_w);
end


%% Calculate c
if lmypoly
    solid_c = in_grid_mypoly(TR_noground.Points,TR_noground.ConnectivityList,TR_noground.incenter,TR_noground.faceNormal,xgrid_c,ygrid_c,zgrid_c,L_char,max_height);
else
    solid_c = inpolyhedron(TR_noground.ConnectivityList, TR_noground.Points, ...
        xgrid_c, ygrid_c, zgrid_c, 'FACENORMALS', TR_noground.faceNormal);
    solid_c = permute(solid_c, [2 1 3]);
end

fluid_c = ~solid_c;

%% Boundary masks
[fluid_IB_c, solid_IB_c] = getBoundaryCells(xgrid_c, ygrid_c, zgrid_c, fluid_c, solid_c, include_diagonals);
if (lgroundfacets)
    fluid_c_1 = fluid_c(:,:,1);
    fluid_IB_c_1 = fluid_IB_c(:,:,1);
    fluid_IB_c_1(fluid_c_1) = true;
    fluid_IB_c(:,:,1) = fluid_IB_c_1;
end

% Boundary coordinates
[fluid_IB_i_c, fluid_IB_j_c, fluid_IB_k_c] = ind2sub(size(fluid_IB_c), find(fluid_IB_c));
fluid_IB_xyz_c = [xgrid_c(fluid_IB_i_c)', ygrid_c(fluid_IB_j_c)', zgrid_c(fluid_IB_k_c)'];

[solid_IB_i_c, solid_IB_j_c, solid_IB_k_c] = ind2sub(size(solid_IB_c), find(solid_IB_c));
solid_IB_xyz_c = [xgrid_c(solid_IB_i_c)', ygrid_c(solid_IB_j_c)', zgrid_c(solid_IB_k_c)'];

% Facet sections
facet_sections_c = matchFacetsToCells(...
    TR, fluid_IB_c, solid_IB_c, fluid_IB_xyz_c, solid_IB_xyz_c, xgrid_c, ygrid_c, zgrid_c, include_diagonals);
%area_fluid_IB_c = sum(facet_sections_c(:,2)); % The total area for each IB cell should equal the sum of the facet

area_facets_c = zeros(nfcts,1);
for n=1:nfcts
    area_facets_c(n) = sum(facet_sections_c(facet_sections_c(:,1)==n, 2));
end
area_fluid_IB_c = sum(area_facets_c);

%% Write c
% Solid points
[solid_i_c, solid_j_c, solid_k_c] = ind2sub(size(solid_c), find(solid_c));
solid_ijk_c = [solid_i_c, solid_j_c, solid_k_c];
filename_c = [fpath 'solid_c.txt'];
fileID_c = fopen(filename_c,'W');
fprintf(fileID_c, '# position (i,j,k)\n');
fclose(fileID_c);
dlmwrite(filename_c, solid_ijk_c, '-append','delimiter',' ');

% Fluid boundary points
fluid_IB_ijk_c = [fluid_IB_i_c, fluid_IB_j_c, fluid_IB_k_c];
fluid_boundary_c = [fluid_IB_ijk_c];%, fluid_IB_rec_c];
filename_c = [fpath 'fluid_boundary_c.txt'];
fileID_c = fopen(filename_c,'W');
fprintf(fileID_c, '# position (i,j,k), distance to surface, reconstruction point location\n');
fclose(fileID_c);
dlmwrite(filename_c, fluid_boundary_c, '-append','delimiter',' ');

% Facet sections
filename_c = [fpath 'facet_sections_c.txt'];
fileID_c = fopen(filename_c,'W');
fprintf(fileID_c, '# facet, area, fluid boundary point, boundary intercept\n');
%fprintf(fileID_c, '%-2d %-4.4f %-5d %-4.4f %-4.4f %-4.4f\n', facet_sections_c(:,[1,2,5,7:9])');
fprintf(fileID_c, '%-2d %-4.4f %-4d %-4.4f\n', facet_sections_c(:,[1,2,5,6])');
fclose(fileID_c);
%dlmwrite(filename_c, facet_sections_c(:,[1,2,5,7:9]), '-append','delimiter',' ','precision',4);

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
    fprintf(fileID_c, '# facet, area, fluid boundary point, boundary intercept\n');
    %fprintf(fileID_c, '%-2d %-4.4f %-4d %-4.4f %-4.4f %-4.4f\n', facet_sections_c_2(:,[1,2,5,7:9])');
    fprintf(fileID_c, '%-2d %-4.4f %-4d %-4.4f\n', facet_sections_c_2(:,[1,2,5,6])');
    fclose(fileID_c);
    %dlmwrite(filename_c, facet_sections_c_2(:,[1,2,5,7:9]), '-append','delimiter',' ','precision',8);
end

%% Reconstruction points - this could be done in fortran
%fluid_IB_rec_u = generateReconstructionPoints(TR, facet_sections_c, TR.faceNormal, fluid_IB_ijk_c, xgrid_c, ygrid_c, zgrid_c);


%% For View3D
fID = fopen([fpath '/facets.vs3'],'w');
fprintf(fID,'T \r\nF  3\r\n');
fprintf(fID,'! %4s %6s %6s %6s\r\n','#', 'x', 'y', 'z');
fprintf(fID,'V %4d %6d %6d %6d\r\n',[(1:size(TR.Points))', TR.Points]');
fprintf(fID,'! %4s %6s %6s %6s %6s %6s %6s %6s %6s\r\n','#', 'v1', 'v2', 'v3', 'v4', 'base', 'cmb', 'emit', 'name' );
fprintf(fID,'S %4d %6d %6d %6d %6d %6d %6d %6d %6d\r\n',[(1:nfcts)', [TR.ConnectivityList            zeros(nfcts,1)], zeros(nfcts,3), (1:nfcts)']');
fprintf(fID,'End of Data\r\n');
fclose(fID);

%%
filename_facets = [fpath 'facets.inp'];
fileID_facets = fopen(filename_facets,'W');
fprintf(fileID_facets, '# type, normal\n');
fclose(fileID_facets);
dlmwrite(filename_facets, [ones(nfcts,1), TR.faceNormal], '-append','delimiter',' ','precision','%1.0f');

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

%% Plot
figure
%trisurf(TR)

patch('Faces', TR.ConnectivityList, 'Vertices', TR.Points, 'FaceColor', ones(3,1)*0.85, 'FaceAlpha', 1)
hold on
% incenters = TR.incenter;
% faceNormals = TR.faceNormal;
% quiver3(incenters(:,1), incenters(:,2), incenters(:,3), faceNormals(:,1), faceNormals(:,2), faceNormals(:,3), 0)
view(3)

axis equal tight

%xlim([0 Lx])
%ylim([0 Ly])
%zlim([0 Lz])

scatter3(X_u(solid_u), Y_u(solid_u), Z_u(solid_u), 10,[0,0,1],'filled') 
%scatter3(X_v(solid_v), Y_v(solid_v), Z_v(solid_v), 10,[0,0,1],'filled')
%scatter3(X_w(solid_w), Y_w(solid_w), Z_w(solid_w), 10,[0,0,1],'filled')
%scatter3(X_c(solid_c), Y_c(solid_c), Z_c(solid_c), 10,[0,0,1],'filled')

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
