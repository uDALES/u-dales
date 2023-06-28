%% writeIBMfilestiledfast

% addpath('./inpolyhedron/')
% addpath('./point2trimesh/')
% addpath('./in_mypoly/')

% Assumes the following variables exist in the workspace:
% TRtile: triangulation describing the geometry.
% xgrid_ut, xgrid_vt, xgrid_wtt, xgrid_ct: the x coordinates on which u,v,w,
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

%TRtile = stlread([fpath2 r2.stl_file]);
TRtile = TR2;
%nfcts = size(TRtile,1);
nfcts = nfcts2;

% c-grid (scalars/pressure)
xgrid_ct = r2.xf;
ygrid_ct = r2.yf;
zgrid_ct = r2.zf;
[X_ct,Y_ct,Z_ct] = ndgrid(xgrid_ct,ygrid_ct,zgrid_ct);

% u-grid
xgrid_ut = r2.xh(1:end-1); %Work smarter not harder
ygrid_ut = r2.yf;
zgrid_ut = r2.zf;
[X_ut,Y_ut,Z_ut] = ndgrid(xgrid_ut,ygrid_ut,zgrid_ut);

% v-grid
xgrid_vt = r2.xf;
ygrid_vt = r2.yh(1:end-1);
zgrid_vt = r2.zf;
[X_vt,Y_vt,Z_vt] = ndgrid(xgrid_vt,ygrid_vt,zgrid_vt);

% w-grid
xgrid_wt = r2.xf;
ygrid_wt = r2.yf;
zgrid_wt = r2.zh;
[X_wt,Y_wt,Z_wt] = ndgrid(xgrid_wt,ygrid_wt,zgrid_wt);


%% Determine total area of surface

if lmypoly
    Dir_ray_u = [0 0 1];
    Dir_ray_v = [0 0 1];
    Dir_ray_w = [0 0 1];
    Dir_ray_c = [0 0 1];
    tol_mypoly = 1e-7;
    max_height = max(TRtile.Points(:,3)) + (zgrid_ut(2)-zgrid_ut(1));
    L_char = 2 * max_facet_size(TRtile.Points,TRtile.ConnectivityList);
end

% if lmypolyfortran
%     write_pre_info;
%     if lwindows
%         in_mypoly_fortran_path = [DA_TOOLSDIR '/IBM/in_mypoly_fortran/'];
%         cd(in_mypoly_fortran_path)
%         system('gfortran -O2 -fopenmp in_mypoly_functions.f90 IBM_flagging.f90 -o pre.exe');
%         copyfile('pre.exe',fpath); % remember to build pre.exe in local system. gfortran -O2 -fopenmp in_mypoly_functions.f90 IBM_flagging.f90 -o pre.exe
%         delete pre.exe in_mypoly_functions.mod;
%         cd(fpath)
%         system('pre.exe'); 
%         delete pre.exe inmypoly_inp_info.txt Stl_data.txt vertices.txt zfgrid.txt zhgrid.txt;
%     else
%         cd(DA_EXPDIR)
%         cd ..
%         in_mypoly_command = ['u-dales/tools/IBM/in_mypoly_fortran/pre_run.sh experiments/' expnr2];
%         system(in_mypoly_command);
%         cd(fpath)
%     end
% else
%     if lmypoly
%         max_height = max(TR.Points(:,3)) + tol_mypoly;
%         L_char = max_facet_size(TR.Points,TR.ConnectivityList) + tol_mypoly;
%     end
% end

%% Calculate u
disp('Determing solid points for u-grid.')
if lmypoly
    solid_u = in_grid_mypoly(TRtile.Points,TRtile.ConnectivityList, ...
        TRtile.incenter,TRtile.faceNormal,xgrid_ut,ygrid_ut,zgrid_ut,Dir_ray_u,L_char,max_height,tol_mypoly);
else
    solid_u = inpolyhedron(TRtile.ConnectivityList, TRtile.Points, ...
        xgrid_ut, ygrid_ut, zgrid_utt, 'FACENORMALS', TRtile.faceNormal);
    solid_u = permute(solid_u, [2 1 3]);
end

fluid_u = ~solid_u;

%% Boundary masks
disp('Determing fluid boundary points for u-grid.')
[fluid_IB_u, solid_IB_u] = getBoundaryCells(xgrid_ut, ygrid_ut, zgrid_ut, fluid_u, solid_u, diag_neighbs);
if (stl_ground)
    fluid_u_1 = fluid_u(:,:,1);
    fluid_IB_u_1 = fluid_IB_u(:,:,1);
    fluid_IB_u_1(fluid_u_1) = true;
    fluid_IB_u(:,:,1) = fluid_IB_u_1;
end

% Boundary coordinates
[fluid_IB_i_u, fluid_IB_j_u, fluid_IB_k_u] = ind2sub(size(fluid_IB_u), find(fluid_IB_u));
fluid_IB_xyz_u = [xgrid_ut(fluid_IB_i_u)', ygrid_ut(fluid_IB_j_u)', zgrid_ut(fluid_IB_k_u)'];

[solid_IB_i_u, solid_IB_j_u, solid_IB_k_u] = ind2sub(size(solid_IB_u), find(solid_IB_u));
solid_IB_xyz_u = [xgrid_ut(solid_IB_i_u)', ygrid_ut(solid_IB_j_u)', zgrid_ut(solid_IB_k_u)'];

% Facet sections
disp('Determining facet sections for u-grid.')
facet_sections_u = matchFacetsToCells(...
    TRtile, fluid_IB_u, solid_IB_u, fluid_IB_xyz_u, solid_IB_xyz_u, xgrid_ut, ygrid_ut, zgrid_ut, diag_neighbs, periodic_x, periodic_y);

%% For debugging - area 'visible' to fluid for each facet.
area_facets_u = zeros(nfcts,1);
for n=1:nfcts
    area_facets_u(n) = sum(facet_sections_u(facet_sections_u(:,1)==n, 2));
end
area_fluid_IB_u = sum(area_facets_u);

%% Write u
% Solid points
[solid_i_u, solid_j_u, solid_k_u] = ind2sub(size(solid_u), find(solid_u));
solid_ijk_u = [solid_i_u, solid_j_u, solid_k_u];


nfcts = nfcts*xtiles*ytiles;
solid_ijk_u_original = solid_ijk_u;
solid_ijk_u_final = [];
for i = 0:xtiles-1
    for j = 0:ytiles-1
        solid_ijk_u_new = solid_ijk_u_original;
        solid_ijk_u_new(:,1) = solid_ijk_u_original(:,1)+i*(length(xgrid_ut));
        solid_ijk_u_new(:,2) = solid_ijk_u_original(:,2)+j*(length(ygrid_ut));
        solid_ijk_u_final = [solid_ijk_u_final; solid_ijk_u_new];
    end 
end
solid_ijk_u = solid_ijk_u_final;
filename_u = [fpath 'solid_u.txt'];
fileID_u = fopen(filename_u,'W');
fprintf(fileID_u, '# position (i,j,k)\n');
fclose(fileID_u);
dlmwrite(filename_u, solid_ijk_u, '-append','delimiter',' ');
disp('Written solid_u.txt')

% Fluid boundary points
fluid_IB_ijk_u = [fluid_IB_i_u, fluid_IB_j_u, fluid_IB_k_u];
fluid_boundary_u = fluid_IB_ijk_u;

fluid_IB_ijk_u_original = fluid_IB_ijk_u;
%reduced_fluid_IB_u = fluid_IB_ijk_u_original;
%reduced_fluid_IB_u(reduced_fluid_IB_u(:,1) == min(reduced_fluid_IB_u(:,1)),:) = [];
fluid_IB_ijk_u_final = [];
ubps_tot = [];
% for i = 0:xtiles-1
%     for j = 0:ytiles-1
%         if i == 0    
%             fluid_IB_ijk_u_new = fluid_IB_ijk_u_original;
%             fluid_IB_ijk_u_new(:,1) = fluid_IB_ijk_u_original(:,1)+i*(length(xgrid_ut));
%             fluid_IB_ijk_u_new(:,2) = fluid_IB_ijk_u_original(:,2)+j*(length(ygrid_ut));
%         else 
%             fluid_IB_ijk_u_new = reduced_fluid_IB_u;
%             fluid_IB_ijk_u_new(:,1) = reduced_fluid_IB_u(:,1)+i*(length(xgrid_ut));
%             fluid_IB_ijk_u_new(:,2) = reduced_fluid_IB_u(:,2)+j*(length(ygrid_ut));
%         end
%         if i+j==0
%             ubps_added = 0;
%         else 
%             ubps_added = length(fluid_IB_ijk_u_new(:,1));
%         end    
%         ubps_tot = [ubps_tot, ubps_added];
%         fluid_IB_ijk_u_final = [fluid_IB_ijk_u_final; fluid_IB_ijk_u_new];
%     end 
% end

for i = 0:xtiles-1
    for j = 0:ytiles-1   
        fluid_IB_ijk_u_new = fluid_IB_ijk_u_original;
        fluid_IB_ijk_u_new(:,1) = fluid_IB_ijk_u_original(:,1)+i*(length(xgrid_ut));
        fluid_IB_ijk_u_new(:,2) = fluid_IB_ijk_u_original(:,2)+j*(length(ygrid_ut));
        if i+j==0
            ubps_added = 0;
        else 
            ubps_added = length(fluid_IB_ijk_u_new(:,1));
        end    
        ubps_tot = [ubps_tot, ubps_added];
        fluid_IB_ijk_u_final = [fluid_IB_ijk_u_final; fluid_IB_ijk_u_new];
    end 
end

ubps_tot = cumsum(ubps_tot);
fluid_boundary_u = fluid_IB_ijk_u_final;
fluid_IB_ijk_u = fluid_boundary_u;
filename_u = [fpath 'fluid_boundary_u.txt'];
fileID_u = fopen(filename_u,'W');
fprintf(fileID_u, '# position (i,j,k)\n');
fclose(fileID_u);
dlmwrite(filename_u, fluid_boundary_u, '-append','delimiter',' ');
disp('Written fluid_boundary_u.txt')

% Facet sections

facet_sections_u_original = facet_sections_u;
facet_sections_u_final = [];
count = 1;
nfcts_u_original = max(facet_sections_u(:,1));
nfbp_u_original = max(facet_sections_u(:,5));
while count <= xtiles*ytiles
    facet_sections_u_new = facet_sections_u_original;
    facet_sections_u_new(:,1) = facet_sections_u_original(:,1)+(count-1)*nfcts_u_original;
    facet_sections_u_new(:,5) = facet_sections_u_original(:,5)+ubps_tot(count);
    facet_sections_u_final = [facet_sections_u_final; facet_sections_u_new];
    count = count + 1; 
end 



% for i = 0:xtiles-1
%     for j = 0:ytiles-1
%         facet_sections_u_new = facet_sections_u_original;
%         facet_sections_u_new(:,1) = facet_sections_u_original(:,1)+count*nfcts_u_original;
%         facet_sections_u_new(:,5) = facet_sections_u_original(:,5)+count*nfbp_u_original;
%         count = count + 1;
%         facet_sections_u_final = [facet_sections_u_final; facet_sections_u_new];
%     end 
% end
facet_sections_u = facet_sections_u_final;
filename_u = [fpath 'facet_sections_u.txt'];
fileID_u = fopen(filename_u,'W');
fprintf(fileID_u, '# facet, area, fluid boundary point, distance\n');
fprintf(fileID_u, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_u(:,[1,2,5,6])');
fclose(fileID_u);
disp('Written facet_sections_u.txt')

lBImin_u = false;
% Instead of using distance of boundary point to facet section, use distance of
% boundary point to ALL facets
if lBImin_u
    facet_sections_u_2 = facet_sections_u;
    [alldists, allBIs, ~, ~] = point2trimesh('Faces', TRtile.ConnectivityList, 'Vertices', TRtile.Points, 'QueryPoints', fluid_IB_xyz_u, 'UseSubSurface', false);
    for n=1:size(facet_sections_u,1)
        m = facet_sections_u(n,5);
        facet_sections_u_2(n,6) = alldists(m);
        facet_sections_u_2(n,7:9) = allBIs(m,:);
    end

%     facet_sections_u_original = facet_sections_u;
%     facet_sections_u_final = [];
%     count = 0;
%     nfcts_u_original = max(facet_sections_u(:,1));
%     nfbp_u_original = max(facet_sections_u(:,5));
%     for i = 0:xtiles-1
%         for j = 0:ytiles-1
%             facet_sections_u_new = facet_sections_u_original;
%             facet_sections_u_new(:,1) = facet_sections_u_original(:,1)+count*nfcts_u_original;
%             facet_sections_u_new(:,5) = facet_sections_u_original(:,5)+count*nfbp_u_original;
%             count = count + 1;
%             facet_sections_u_final = [facet_sections_u_final; facet_sections_u_new];
%         end 
%     end
% 
%     facet_sections_u_2 = facet_sections_u_2_final
    filename_u = [fpath 'facet_sections_u_2.txt'];
    fileID_u = fopen(filename_u,'W');
    fprintf(fileID_u, '# facet, area, fluid boundary point, distance\n');
    fprintf(fileID_u, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_u_2(:,[1,2,5,6])');
    fclose(fileID_u);
    disp('Written facet_sections_u_2.txt')
end

%% Calculate v
disp('Determing solid points for v-grid.')
if lmypoly
    solid_v = in_grid_mypoly(TRtile.Points,TRtile.ConnectivityList,TRtile.incenter,TRtile.faceNormal,xgrid_vt,ygrid_vt,zgrid_vt,Dir_ray_v,L_char,max_height,tol_mypoly);
else
    solid_v = inpolyhedron(TRtile.ConnectivityList, TRtile.Points, ...
        xgrid_vt, ygrid_vt, zgrid_vt, 'FACENORMALS', TRtile.faceNormal);
    solid_v = permute(solid_v, [2 1 3]);
end

fluid_v = ~solid_v;
%%
disp('Determing fluid boundary points for v-grid.')
% Boundary masks
[fluid_IB_v, solid_IB_v] = getBoundaryCells(xgrid_vt, ygrid_vt, zgrid_vt, fluid_v, solid_v, diag_neighbs);
if (stl_ground)
    fluid_v_1 = fluid_v(:,:,1);
    fluid_IB_v_1 = fluid_IB_v(:,:,1);
    fluid_IB_v_1(fluid_v_1) = true;
    fluid_IB_v(:,:,1) = fluid_IB_v_1;
end

% Boundary coordinates
[fluid_IB_i_v, fluid_IB_j_v, fluid_IB_k_v] = ind2sub(size(fluid_IB_v), find(fluid_IB_v));
fluid_IB_xyz_v = [xgrid_vt(fluid_IB_i_v)', ygrid_vt(fluid_IB_j_v)', zgrid_vt(fluid_IB_k_v)'];

[solid_IB_i_v, solid_IB_j_v, solid_IB_k_v] = ind2sub(size(solid_IB_v), find(solid_IB_v));
solid_IB_xyz_v = [xgrid_vt(solid_IB_i_v)', ygrid_vt(solid_IB_j_v)', zgrid_vt(solid_IB_k_v)'];

%% Facet sections
disp('Determing facet sections for v-grid.')
facet_sections_v = matchFacetsToCells(...
    TRtile, fluid_IB_v, solid_IB_v, fluid_IB_xyz_v, solid_IB_xyz_v, xgrid_vt, ygrid_vt, zgrid_vt, diag_neighbs, periodic_x, periodic_y);

area_facets_v = zeros(nfcts,1);
for n=1:nfcts
    area_facets_v(n) = sum(facet_sections_v(facet_sections_v(:,1)==n, 2));
end
area_fluid_IB_v = sum(area_facets_v);


%% Write v
% Solid points
[solid_i_v, solid_j_v, solid_k_v] = ind2sub(size(solid_v), find(solid_v));
solid_ijk_v = [solid_i_v, solid_j_v, solid_k_v];

solid_ijk_v_original = solid_ijk_v;
solid_ijk_v_final = [];
for i = 0:xtiles-1
    for j = 0:ytiles-1
        solid_ijk_v_new = solid_ijk_v_original;
        solid_ijk_v_new(:,1) = solid_ijk_v_original(:,1)+i*(length(xgrid_vt));
        solid_ijk_v_new(:,2) = solid_ijk_v_original(:,2)+j*(length(ygrid_vt)-1);
        solid_ijk_v_final = [solid_ijk_v_final; solid_ijk_v_new];
    end 
end
solid_ijk_v = solid_ijk_v_final;
filename_v = [fpath 'solid_v.txt'];
fileID_v = fopen(filename_v,'W');
fprintf(fileID_v, '# position (i,j,k)\n');
fclose(fileID_v);
dlmwrite(filename_v, solid_ijk_v, '-append','delimiter',' ');
disp('Written solid_v.txt')

% Fluid boundary points
fluid_IB_ijk_v = [fluid_IB_i_v, fluid_IB_j_v, fluid_IB_k_v];
fluid_boundary_v = fluid_IB_ijk_v;

fluid_IB_ijk_v_original = fluid_IB_ijk_v;
% reduced_fluid_IB_v = fluid_IB_ijk_v_original;
% reduced_fluid_IB_v(reduced_fluid_IB_v(:,2) == min(reduced_fluid_IB_v(:,2)),:) = [];
fluid_IB_ijk_v_final = [];
vbps_tot = [];
for i = 0:xtiles-1
    for j = 0:ytiles-1   
        fluid_IB_ijk_v_new = fluid_IB_ijk_v_original;
        fluid_IB_ijk_v_new(:,1) = fluid_IB_ijk_v_original(:,1)+i*(length(xgrid_vt));
        fluid_IB_ijk_v_new(:,2) = fluid_IB_ijk_v_original(:,2)+j*(length(ygrid_vt));
        if i+j==0
            vbps_added = 0;
        else 
            vbps_added = length(fluid_IB_ijk_v_new(:,2));
        end    
        vbps_tot = [vbps_tot, vbps_added];
        fluid_IB_ijk_v_final = [fluid_IB_ijk_v_final; fluid_IB_ijk_v_new];
    end 
end
% for i = 0:xtiles-1
%     for j = 0:ytiles-1
%         if j == 0    
%             fluid_IB_ijk_v_new = fluid_IB_ijk_v_original;
%             fluid_IB_ijk_v_new(:,1) = fluid_IB_ijk_v_original(:,1)+i*(length(xgrid_vt));
%             fluid_IB_ijk_v_new(:,2) = fluid_IB_ijk_v_original(:,2)+j*(length(ygrid_vt)-1);
%         else 
%             fluid_IB_ijk_v_new = reduced_fluid_IB_v;
%             fluid_IB_ijk_v_new(:,1) = reduced_fluid_IB_v(:,1)+i*(length(xgrid_vt));
%             fluid_IB_ijk_v_new(:,2) = reduced_fluid_IB_v(:,2)+j*(length(ygrid_vt)-1);
%         end
%         if i+j==0
%             vbps_added = 0;
%         else 
%             vbps_added = length(fluid_IB_ijk_v_new(:,2));
%         end    
%         vbps_tot = [vbps_tot, vbps_added];
%         fluid_IB_ijk_v_final = [fluid_IB_ijk_v_final; fluid_IB_ijk_v_new];
%     end 
% end
vbps_tot = cumsum(vbps_tot);
fluid_boundary_v = fluid_IB_ijk_v_final;
fluid_IB_ijk_v = fluid_boundary_v;
filename_v = [fpath 'fluid_boundary_v.txt'];
fileID_v = fopen(filename_v,'W');
fprintf(fileID_v, '# position (i,j,k)\n');
fclose(fileID_v);
dlmwrite(filename_v, fluid_boundary_v, '-append','delimiter',' ');
disp('Written fluid_boundary_v.txt')

% Facet sections

facet_sections_v_original = facet_sections_v;
facet_sections_v_final = [];
nfcts_v_original = max(facet_sections_v(:,1));
nfbp_v_original = max(facet_sections_v(:,5));
count = 1;
while count <= xtiles*ytiles
    facet_sections_v_new = facet_sections_v_original;
    facet_sections_v_new(:,1) = facet_sections_v_original(:,1)+(count-1)*nfcts_v_original;
    facet_sections_v_new(:,5) = facet_sections_v_original(:,5)+vbps_tot(count);
    facet_sections_v_final = [facet_sections_v_final; facet_sections_v_new];
    count = count + 1; 
end 

% for i = 0:xtiles-1
%     for j = 0:ytiles-1
%         facet_sections_v_new = facet_sections_v_original;
%         facet_sections_v_new(:,1) = facet_sections_v_original(:,1)+count*nfcts_v_original;
%         facet_sections_v_new(:,5) = facet_sections_v_original(:,5)+count*nfbp_v_original;
%         count = count + 1;
%         facet_sections_v_final = [facet_sections_v_final; facet_sections_v_new];
%     end 
% end

facet_sections_v = facet_sections_v_final;
filename_v = [fpath 'facet_sections_v.txt'];
fileID_v = fopen(filename_v,'W');
fprintf(fileID_v, '# facet, area, fluid boundary point, distance\n');
fprintf(fileID_v, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_v(:,[1,2,5,6])');
fclose(fileID_v);
disp('Written facet_sections_v.txt')

lBImin_v = false;
% Instead of using distance of boundary point to facet section, use distance of
% boundary point to ALL facets
if lBImin_v
    facet_sections_v_2 = facet_sections_v;
    [alldists, allBIs, ~, ~] = point2trimesh('Faces', TRtile.ConnectivityList, 'Vertices', TRtile.Points, 'QueryPoints', fluid_IB_xyz_v, 'UseSubSurface', false);
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
end

%% Calculate w
disp('Determing solid points for w-grid.')
if lmypoly
    solid_w = in_grid_mypoly(TRtile.Points,TRtile.ConnectivityList, ...
        TRtile.incenter,TRtile.faceNormal,xgrid_wt,ygrid_wt,zgrid_wt,Dir_ray_w,L_char,max_height,tol_mypoly);
else
    solid_w = inpolyhedron(TRtile.ConnectivityList, TRtile.Points, ...
        xgrid_wt, ygrid_wt, zgrid_wt, 'FACENORMALS', TRtile.faceNormal);
    solid_w = permute(solid_w, [2 1 3]);
end

if (stl_ground)
    solid_w_b = solid_w;
    solid_w_b(:,:,1) = true;
end

fluid_w = ~solid_w;

%% Boundary points
disp('Determing fluid boundary points for w-grid.')
% Boundary masks
[fluid_IB_w, solid_IB_w] = getBoundaryCells(xgrid_wt, ygrid_wt, zgrid_wt, fluid_w, solid_w_b, diag_neighbs);
fluid_IB_w(:,:,1) = 0; % Bottom is always solid

% Boundary coordinates
[fluid_IB_i_w, fluid_IB_j_w, fluid_IB_k_w] = ind2sub(size(fluid_IB_w), find(fluid_IB_w));
fluid_IB_xyz_w = [xgrid_wt(fluid_IB_i_w)', ygrid_wt(fluid_IB_j_w)', zgrid_wt(fluid_IB_k_w)'];

[solid_IB_i_w, solid_IB_j_w, solid_IB_k_w] = ind2sub(size(solid_IB_w), find(solid_IB_w));
solid_IB_xyz_w = [xgrid_wt(solid_IB_i_w)', ygrid_wt(solid_IB_j_w)', zgrid_wt(solid_IB_k_w)'];

%% Facet sections
disp('Determing facet sections for w-grid.')
facet_sections_w = matchFacetsToCells(...
    TRtile, fluid_IB_w, solid_IB_w, fluid_IB_xyz_w, solid_IB_xyz_w, xgrid_wt, ygrid_wt, zgrid_wt, diag_neighbs, periodic_x, periodic_y);

area_facets_w = zeros(nfcts,1);
for n=1:nfcts
    area_facets_w(n) = sum(facet_sections_w(facet_sections_w(:,1)==n, 2));
end
area_fluid_IB_w = sum(area_facets_w);

%% Write w
% Solid points
[solid_i_w, solid_j_w, solid_k_w] = ind2sub(size(solid_w), find(solid_w));
solid_ijk_w = [solid_i_w, solid_j_w, solid_k_w];

solid_ijk_w_original = solid_ijk_w;
solid_ijk_w_final = [];
for i = 0:xtiles-1
    for j = 0:ytiles-1
        solid_ijk_w_new = solid_ijk_w_original;
        solid_ijk_w_new(:,1) = solid_ijk_w_original(:,1)+i*(length(xgrid_wt));
        solid_ijk_w_new(:,2) = solid_ijk_w_original(:,2)+j*(length(ygrid_wt));
        solid_ijk_w_final = [solid_ijk_w_final; solid_ijk_w_new];
    end 
end
solid_ijk_w = solid_ijk_w_final;
filename_w = [fpath 'solid_w.txt'];
fileID_w = fopen(filename_w,'W');
fprintf(fileID_w, '# position (i,j,k)\n');
fclose(fileID_w);
dlmwrite(filename_w, solid_ijk_w, '-append','delimiter',' ');
disp('Written solid_w.txt')

% Fluid boundary points
fluid_IB_ijk_w = [fluid_IB_i_w, fluid_IB_j_w, fluid_IB_k_w];
fluid_boundary_w = fluid_IB_ijk_w;

fluid_IB_ijk_w_original = fluid_IB_ijk_w;
fluid_IB_ijk_w_final = [];
for i = 0:xtiles-1
    for j = 0:ytiles-1
        fluid_IB_ijk_w_new = fluid_IB_ijk_w_original;
        fluid_IB_ijk_w_new(:,1) = fluid_IB_ijk_w_original(:,1)+i*(length(xgrid_wt));
        fluid_IB_ijk_w_new(:,2) = fluid_IB_ijk_w_original(:,2)+j*(length(ygrid_wt));
        fluid_IB_ijk_w_final = [fluid_IB_ijk_w_final; fluid_IB_ijk_w_new];
    end 
end
fluid_boundary_w = fluid_IB_ijk_w_final;
fluid_IB_ijk_w = fluid_boundary_w;
filename_w = [fpath 'fluid_boundary_w.txt'];
fileID_w = fopen(filename_w,'W');
fprintf(fileID_w, '# position (i,j,k)\n');
fclose(fileID_w);
dlmwrite(filename_w, fluid_boundary_w, '-append','delimiter',' ');
disp('Written fluid_boundary_w.txt')

% Facet sections

facet_sections_w_original = facet_sections_w;
facet_sections_w_final = [];
count = 0;
nfcts_w_original = max(facet_sections_w(:,1));
nfbp_w_original = max(facet_sections_w(:,5));
for i = 0:xtiles-1
    for j = 0:ytiles-1
        facet_sections_w_new = facet_sections_w_original;
        facet_sections_w_new(:,1) = facet_sections_w_original(:,1)+count*nfcts_w_original;
        facet_sections_w_new(:,5) = facet_sections_w_original(:,5)+count*nfbp_w_original;
        count = count + 1;
        facet_sections_w_final = [facet_sections_w_final; facet_sections_w_new];
    end 
end

facet_sections_w = facet_sections_w_final;
filename_w = [fpath 'facet_sections_w.txt'];
fileID_w = fopen(filename_w,'W');
fprintf(fileID_w, '# facet, area, fluid boundary point, distance\n');

fprintf(fileID_w, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_w(:,[1,2,5,6])');
fclose(fileID_w);
disp('Written facet_sections_w.txt')

lBImin_w = false;
% Instead of using distance of boundary point to facet section, use distance of
% boundary point to ALL facets
if lBImin_w
    facet_sections_w_2 = facet_sections_w;
    [alldists, allBIs, ~, ~] = point2trimesh('Faces', TRtile.ConnectivityList, 'Vertices', TRtile.Points, 'QueryPoints', fluid_IB_xyz_w, 'UseSubSurface', false);
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
end


%% Calculate c
disp('Determing solid points for c-grid.')
if lmypoly
    solid_c = in_grid_mypoly(TRtile.Points,TRtile.ConnectivityList,TRtile.incenter,TRtile.faceNormal,xgrid_ct,ygrid_ct,zgrid_ct,Dir_ray_c,L_char,max_height,tol_mypoly);
else
    solid_c = inpolyhedron(TRtile.ConnectivityList, TRtile.Points, ...
        xgrid_ct, ygrid_ct, zgrid_ct, 'FACENORMALS', TRtile.faceNormal);
    solid_c = permute(solid_c, [2 1 3]);
end

fluid_c = ~solid_c;

%% Boundary masks
disp('Determing fluid boundary points for c-grid.')
[fluid_IB_c, solid_IB_c] = getBoundaryCells(xgrid_ct, ygrid_ct, zgrid_ct, fluid_c, solid_c, diag_neighbs);
if (stl_ground)
    fluid_c_1 = fluid_c(:,:,1);
    fluid_IB_c_1 = fluid_IB_c(:,:,1);
    fluid_IB_c_1(fluid_c_1) = true;
    fluid_IB_c(:,:,1) = fluid_IB_c_1;
end

% Boundary coordinates
[fluid_IB_i_c, fluid_IB_j_c, fluid_IB_k_c] = ind2sub(size(fluid_IB_c), find(fluid_IB_c));
fluid_IB_xyz_c = [xgrid_ct(fluid_IB_i_c)', ygrid_ct(fluid_IB_j_c)', zgrid_ct(fluid_IB_k_c)'];

[solid_IB_i_c, solid_IB_j_c, solid_IB_k_c] = ind2sub(size(solid_IB_c), find(solid_IB_c));
solid_IB_xyz_c = [xgrid_ct(solid_IB_i_c)', ygrid_ct(solid_IB_j_c)', zgrid_ct(solid_IB_k_c)'];

%% Facet sections
disp('Determing facet sections for c-grid.')
facet_sections_c = matchFacetsToCells(...
    TRtile, fluid_IB_c, solid_IB_c, fluid_IB_xyz_c, solid_IB_xyz_c, xgrid_ct, ygrid_ct, zgrid_ct, diag_neighbs, periodic_x, periodic_y);

area_facets_c = zeros(nfcts,1);
for n=1:nfcts
    area_facets_c(n) = sum(facet_sections_c(facet_sections_c(:,1)==n, 2));
end
area_fluid_IB_c = sum(area_facets_c);

%% Write c
% Solid points
[solid_i_c, solid_j_c, solid_k_c] = ind2sub(size(solid_c), find(solid_c));
solid_ijk_c = [solid_i_c, solid_j_c, solid_k_c];

solid_ijk_c_original = solid_ijk_c;
solid_ijk_c_final = [];
for i = 0:xtiles-1
    for j = 0:ytiles-1
        solid_ijk_c_new = solid_ijk_c_original;
        solid_ijk_c_new(:,1) = solid_ijk_c_original(:,1)+i*(length(xgrid_ct));
        solid_ijk_c_new(:,2) = solid_ijk_c_original(:,2)+j*(length(ygrid_ct));
        solid_ijk_c_final = [solid_ijk_c_final; solid_ijk_c_new];
    end 
end
solid_ijk_c = solid_ijk_c_final;

filename_c = [fpath 'solid_c.txt'];
fileID_c = fopen(filename_c,'W');
fprintf(fileID_c, '# position (i,j,k)\n');
fclose(fileID_c);
dlmwrite(filename_c, solid_ijk_c, '-append','delimiter',' ');
disp('Written solid_c.txt')

% Fluid boundary points
fluid_IB_ijk_c = [fluid_IB_i_c, fluid_IB_j_c, fluid_IB_k_c];
fluid_boundary_c = fluid_IB_ijk_c;

fluid_IB_ijk_c_original = fluid_IB_ijk_c;
fluid_IB_ijk_c_final = [];
for i = 0:xtiles-1
    for j = 0:ytiles-1
        fluid_IB_ijk_c_new = fluid_IB_ijk_c_original;
        fluid_IB_ijk_c_new(:,1) = fluid_IB_ijk_c_original(:,1)+i*(length(xgrid_ct));
        fluid_IB_ijk_c_new(:,2) = fluid_IB_ijk_c_original(:,2)+j*(length(ygrid_ct));
        fluid_IB_ijk_c_final = [fluid_IB_ijk_c_final; fluid_IB_ijk_c_new];
    end 
end
fluid_boundary_c = fluid_IB_ijk_c_final;
fluid_IB_ijk_c = fluid_boundary_c;
filename_c = [fpath 'fluid_boundary_c.txt'];
fileID_c = fopen(filename_c,'W');
fprintf(fileID_c, '# position (i,j,k), distance to surface, reconstruction point location\n');
fclose(fileID_c);
dlmwrite(filename_c, fluid_boundary_c, '-append','delimiter',' ');
disp('Written fluid_boundary_c.txt')

% Facet sections

facet_sections_c_original = facet_sections_c;
facet_sections_c_final = [];
count = 0;
nfcts_c_original = max(facet_sections_c(:,1));
nfbp_c_original = max(facet_sections_c(:,5));
for i = 0:xtiles-1
    for j = 0:ytiles-1
        facet_sections_c_new = facet_sections_c_original;
        facet_sections_c_new(:,1) = facet_sections_c_original(:,1)+count*nfcts_c_original;
        facet_sections_c_new(:,5) = facet_sections_c_original(:,5)+count*nfbp_c_original;
        count = count + 1;
        facet_sections_c_final = [facet_sections_c_final; facet_sections_c_new];
    end 
end

facet_sections_c = facet_sections_c_final;
filename_c = [fpath 'facet_sections_c.txt'];
fileID_c = fopen(filename_c,'W');
fprintf(fileID_c, '# facet, area, fluid boundary point, distance\n');
fprintf(fileID_c, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_c(:,[1,2,5,6])');
fclose(fileID_c);
disp('Written facet_sections_c.txt')

lBImin_c = false;
% Instead of using distance of boundary point to facet section, use distance of
% boundary point to ALL facets
if lBImin_c
    facet_sections_c_2 = facet_sections_c;
    [alldists, allBIs, ~, ~] = point2trimesh('Faces', TRtile.ConnectivityList, 'Vertices', TRtile.Points, 'QueryPoints', fluid_IB_xyz_c, 'UseSubSurface', false);
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
%fprintf(fileID_info, ['nbndpts_u = ', num2str(size(fluid_IB_xyz_u,1)), '\n']);
fprintf(fileID_info, ['nbndpts_u = ', num2str(size(fluid_IB_ijk_u,1)), '\n']);
% fprintf(fileID_info, ['nbndpts_v = ', num2str(size(fluid_IB_xyz_v,1)), '\n']);
% fprintf(fileID_info, ['nbndpts_w = ', num2str(size(fluid_IB_xyz_w,1)), '\n']);
% fprintf(fileID_info, ['nbndpts_c = ', num2str(size(fluid_IB_xyz_c,1)), '\n']);
fprintf(fileID_info, ['nbndpts_v = ', num2str(size(fluid_IB_ijk_v,1)), '\n']);
fprintf(fileID_info, ['nbndpts_w = ', num2str(size(fluid_IB_ijk_w,1)), '\n']);
fprintf(fileID_info, ['nbndpts_c = ', num2str(size(fluid_IB_ijk_c,1)), '\n']);
fprintf(fileID_info, ['nfctsecs_u = ', num2str(size(facet_sections_u,1)), '\n']);
fprintf(fileID_info, ['nfctsecs_v = ', num2str(size(facet_sections_v,1)), '\n']);
fprintf(fileID_info, ['nfctsecs_w = ', num2str(size(facet_sections_w,1)), '\n']);
fprintf(fileID_info, ['nfctsecs_c = ', num2str(size(facet_sections_c,1)), '\n']);
fclose(fileID_info);

%% Plot
figure

%patch('Faces', TRtile.ConnectivityList, 'Vertices', TRtile.Points, 'FaceColor', ones(3,1)*0.85, 'FaceAlpha', 1)
patch('Faces', TRtile.ConnectivityList, 'Vertices', TRtile.Points, 'FaceColor', ones(3,1)*0.85, 'FaceAlpha', 1)
hold on
incenters = TRtile.incenter;
faceNormals = TRtile.faceNormal;
%quiver3(incenters(:,1), incenters(:,2), incenters(:,3), faceNormals(:,1), faceNormals(:,2), faceNormals(:,3), 0)
view(3)

axis equal tight

% xlim([0 Lx])
% ylim([0 Ly])
% zlim([0 Lz])

% solid_ut = solid_u;
% solid_vt = solid_v;
% solid_wt = solid_w;
% solid_ct = solid_c;
% 
% fluid_IB_ut = fluid_IB_u;
% fluid_IB_vt = fluid_IB_v;
% fluid_IB_wt = fluid_IB_w;
% fluid_IB_ct = fluid_IB_c;

scatter3(X_ut(solid_u), Y_ut(solid_u), Z_ut(solid_u), 10,[0,0,1],'filled')
%scatter3(X_vtt(solid_v), Y_vt(solid_v), Z_vt(solid_v), 10,[0,0,1],'filled')
%scatter3(X_wtt(solid_w), Y_wt(solid_w), Z_wt(solid_w), 10,[0,0,1],'filled')
%scatter3(X_ct(solid_c), Y_ct(solid_c), Z_ct(solid_c), 10,[0,0,1],'filled')

%scatter3(X_ut(fluid_IB_u), Y_ut(fluid_IB_u), Z_ut(fluid_IB_u), 10,[0,0,1],'filled')
%scatter3(X_vt(fluid_IB_v), Y_vt(fluid_IB_v), Z_vt(fluid_IB_v), 10,[0,0,1],'filled')
%scatter3(X_wt(fluid_IB_w), Y_ut(fluid_IB_w), Z_ut(fluid_IB_w), 10,[0,0,1],'filled')
%scatter3(X_ct(fluid_IB_c), Y_ut(fluid_IB_c), Z_ut(fluid_IB_c), 10,[0,0,1],'filled')

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
%scatter3(fluid_IB_xyz_c(:,1),fluid_IB_xyz_c(:,2),fluid_IB_xyz_c(:,3),10,[0,0,1],'filled')
%scatter3(solid_IB_xyz_c(:,1),solid_IB_xyz_c(:,2),solid_IB_xyz_c(:,3),10,[0,0,1],'filled')
 %quiver3(fluid_IB_BI_c(:,1), fluid_IB_BI_c(:,2), fluid_IB_BI_c(:,3), fluid_IB_vec_c(:,1), fluid_IB_vec_c(:,2), fluid_IB_vec_c(:,3),'off')
% %scatter3(fluid_IB_rec_c(:,1),fluid_IB_rec_c(:,2),fluid_IB_rec_c(:,3),10,[0,0,1],'filled')
% quiver3(fluid_IB_xyz_c(:,1), fluid_IB_xyz_c(:,2), fluid_IB_xyz_c(:,3), fluid_IB_rec_vec_c(:,1), fluid_IB_rec_vec_c(:,2), fluid_IB_rec_vec_c(:,3),'off')