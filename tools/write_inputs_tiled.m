%% write_inputs_tiled 

% uDALES (https://github.com/uDALES/u-dales).

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Copyright (C) 2016-2023 the uDALES Team.

% This script is run by the bash script da_inp.sh.
% It used to generate the necessary input files for uDALES.
tic
expnr = '048';
expnr2 = '048';
tiled =false;
xtiles = 1;
ytiles = 1;
% DA_EXPDIR = getenv('DA_EXPDIR');
% DA_TOOLSDIR = getenv('DA_TOOLSDIR');
DA_EXPDIR = '/media/chris/Project3/uDALES2.0/experiments'
DA_TOOLSDIR = '/media/chris/Project3/uDALES2.0/u-dales/tools'
addpath(genpath([DA_TOOLSDIR '/']));
addpath([DA_TOOLSDIR '/IBM/'])
addpath([DA_TOOLSDIR '/SEB/'])
addpath([DA_TOOLSDIR '/setting/'])
exppath = [DA_EXPDIR '/'];
fpath = [DA_EXPDIR '/' expnr '/'];
fpath2 = [DA_EXPDIR '/' expnr2 '/'];
addpath(fpath2);
cd(fpath)


r = preprocessing(expnr, exppath); % reads namoptions file and creates the object r
preprocessing.set_defaults(r);
preprocessing.generate_xygrid(r);
preprocessing.generate_zgrid(r);
preprocessing.generate_lscale(r)
preprocessing.write_lscale(r)
disp(['Written lscal.inp.', r.expnr])
preprocessing.generate_prof(r);
preprocessing.write_prof(r);
disp(['Written prof.inp.', r.expnr])

r2 = preprocessing(expnr2,exppath); % for the tiled one
preprocessing.set_defaults(r2);
preprocessing.generate_xygrid(r2);
preprocessing.generate_zgrid(r2);

% if r.nsv>0
%     preprocessing.generate_scalar(r);
%     preprocessing.write_scalar(r);
%     disp(['Written scalar.inp.', r.expnr])
% end

if isfile(['factypes.inp.', expnr])
    r.factypes = dlmread(['factypes.inp.', r.expnr],'',3,0);
else
    preprocessing.write_factypes(r)
    disp(['Written factypes.inp.', r.expnr])
end


if r.libm
    %% Read the .stl file and write necessary ibm files
    TR = stlread(r.stl_file);
    TR2 = stlread(r2.stl_file);
    F = TR.ConnectivityList;
    V = TR.Points;
    F2 = TR2.ConnectivityList;
    V2 = TR2.Points;
    %%
    area_facets = facetAreas(F, V); % Useful for checking if area_fluid_IB_c == sum(area_facets)
    area_facets2 = facetAreas(F2,V2);
    %%

    % Set facet types
    nfcts = size(TR.ConnectivityList,1);
    nfcts2 = size(TR2.ConnectivityList,1);
    preprocessing.set_nfcts(r, nfcts);
    preprocessing.set_nfcts(r2, nfcts2);
    facet_types = ones(nfcts,1); % facet_types are to be user-defined - defaults to type 1 (concrete)
    preprocessing.write_facets(r, facet_types, TR.faceNormal);
    disp(['Written facets.inp.', r.expnr])

    calculate_facet_sections_uvw = r.iwallmom > 1;
    calculate_facet_sections_c = r.ltempeq || r.lmoist
    % c-grid (scalars/pressure)
    xgrid_c = r.xf;
    ygrid_c = r.yf;
    zgrid_c = r.zf;

    [X_c,Y_c,Z_c] = ndgrid(xgrid_c,ygrid_c,zgrid_c);

    % u-grid
    xgrid_u = r.xh(1:end-1);
    ygrid_u = r.yf;
    zgrid_u = r.zf;
    [X_u,Y_u,Z_u] = ndgrid(xgrid_u,ygrid_u,zgrid_u);

    % v-grid
    xgrid_v = r.xf;
    ygrid_v = r.yh(1:end-1);
    zgrid_v = r.zf;
    [X_v,Y_v,Z_v] = ndgrid(xgrid_v,ygrid_v,zgrid_v);

    % w-grid
    xgrid_w = r.xf;
    ygrid_w = r.yf;
    zgrid_w = r.zh(1:end-1);
    [X_w,Y_w,Z_w] = ndgrid(xgrid_w,ygrid_w,zgrid_w);

    diag_neighbs = r.diag_neighbs;
    stl_ground = r.stl_ground;
    periodic_x = r.BCxm == 1;
    periodic_y = r.BCym == 1;
    xsize = r.xlen;
    ysize = r.ylen;
    zsize = r.zsize;
    %lmypoly = 1; % remove eventually

%         Dir_ray_u = [0 0 1];
%         Dir_ray_v = [0 0 1];
%         Dir_ray_w = [0 0 1];
%         Dir_ray_c = [0 0 1];
%         tol_mypoly = 1e-4;
% 
%         write_pre_info;
%     lmypolyfortran = 1; lmypoly = 0;		% remove eventually
%     lwindows = false;
%     Dir_ray_u = [0 0 1];
%     Dir_ray_v = [0 0 1];
%     Dir_ray_w = [0 0 1];
%     Dir_ray_c = [0 0 1];
%     tol_mypoly = 5e-4;
    writetiledIBMFiles
    copy_command = ['cp ' r.geom_path 'solid_* ' r.geom_path 'fluid_boundary_* ' fpath];
    system(copy_command);
    copy_command = ['cp ' r.geom_path 'fluid_boundary_* ' fpath];
    system(copy_command);
    if calculate_facet_sections_uvw
        copy_command = ['cp ' r.geom_path 'facet_sections_u* ' fpath];
        system(copy_command);
        copy_command = ['cp ' r.geom_path 'facet_sections_v* ' fpath];
        system(copy_command);
        copy_command = ['cp ' r.geom_path 'facet_sections_w* ' fpath];
        system(copy_command);
    end
    if calculate_facet_sections_c
        copy_command = ['cp ' r.geom_path 'facet_sections_c* ' fpath];
        system(copy_command);
    end
end

%% Set facet types
nfcts = size(TR.ConnectivityList,1);
%preprocessing.addvar(r, 'nfcts', nfcts);
preprocessing.set_nfcts(r, nfcts);
facet_types = ones(nfcts,1); % facet_types are to be user-defined - defaults to type 1 (concrete)
preprocessing.write_facets(r, facet_types, TR.faceNormal);
preprocessing.write_facetarea(r, area_facets);
%%
toc 
if r.lEB
    %% Write STL in View3D input format
    fpath_facets_view3d = [fpath 'facets_tile.vs3']; %Write the tile vs3 to the main exp
    STLtoView3D(r2.stl_file, fpath_facets_view3d); %Make a vs3 file of the tile

    %% Calculate view factors
    % Add check to see if View3D exists in the tools directory.
    view3d_exe = [DA_TOOLSDIR '/View3D/build/src/view3d'];
    fpath_vf = [fpath 'vf_tile.txt']; %Write to main geom folder
    vftile = view3d(view3d_exe, fpath_facets_view3d, fpath_vf); %Now we have a reduced size vf
    %%
    vf = diagonal_matrix_maker(vftile,xtiles,ytiles);
    toc
    svftile = max(1-sum(vftile,2),0);
    svf = max(1 - sum(vf, 2), 0);
    preprocessing.write_svf(r, svf);
    preprocessing.write_svf(r2,svftile);
    %%
    if ~r.lvfsparse
        preprocessing.write_vf(r, vf)
        disp(['Written vf.nc.inp.', r.expnr])
    else
        vfsparse = sparse(double(vf));
        preprocessing.write_vfsparse(r, vfsparse);
        disp(['Written vfsparse.inp.', r.expnr])
    end
    %% Set facet types
    nfcts = size(TR.ConnectivityList,1);
    preprocessing.set_nfcts(r, nfcts);
    facet_types = ones(nfcts,1); % facet_types are to be user-defined - defaults to type 1 (concrete)
    preprocessing.write_facets(r, facet_types, TR.faceNormal);

        %% Calculate direct solar radiation (Sdir)
        disp('Calculating direct solar radiation.')
        azimuth = r.solarazimuth - r.xazimuth;
        nsun = [sind(r.solarzenith)*cosd(azimuth), -sind(r.solarzenith)*sind(azimuth), cosd(r.solarzenith)];
        show_plot_2d = false; % User-defined
        show_plot_3d = true;  % User-defined
        %Sdir = directShortwave(F, V, nsun, r.I, r.psc_res, show_plot_2d, show_plot_3d);
        Sdirtile = directShortwave(F2,V2,nsun, r.I, r.psc_res, show_plot_2d, show_plot_3d);
        toc 
        %% Calculate net shortwave radiation (Knet)
        disp('Calculating net shortwave radiation.')
        albedos2 = preprocessing.generate_albedos(r2, facet_types);
        %Knet = netShortwave(Sdir, r.Dsky, vf, svf, albedos);
        Knet_tile = netShortwave(Sdirtile, r2.Dsky, vftile, svftile, albedos2);
        tiles = xtiles*ytiles;
        Knet = zeros(length(Knet_tile)*tiles,1);
        for i = 0:tiles-1
            xl = 1+i*length(Knet_tile);
            xu = length(Knet_tile) + i*length(Knet_tile);
            Knet(xl:xu) = Knet_tile;
        end 
        toc
        preprocessing.write_netsw(r2, Knet_tile);
        preprocessing.write_netsw(r, Knet);
        disp(['Written netsw.inp.', r.expnr])
    end

    %% Write initial facet temperatures
    if (r.lEB || r.iwallmom == 2 || r.iwalltemp == 2)
        disp('Setting initial facet temperatures.')
        facT = r.facT;
        nfaclyrs = r.nfaclyrs;
        facT_file = r.facT_file;
        lfacTlyrs = r.lfacTlyrs;
        if ~r.lfacTlyrs
            Tfacinit = ones(nfcts,1) .* r.facT;;
            preprocessing.write_Tfacinit(r, Tfacinit);
            disp(['Written Tfacinit.inp.', r.expnr]);
            % Could always read in facet temperature as layers, defaulting to linear?
        else
            Tfac = ncread(r.facT_file, 'T');
            Tfacinit_layers = Tfac(:, :, end);
            preprocessing.write_Tfacinit_layers(r, Tfacinit_layers);
            disp(['Written Tfacinit_layers.inp.', r.expnr]);
        end
    end
%% Setting vars
lamdba_calculation
setting_types

%% Determine effective albedo
%efctvalb = 1-sum(Knet)/sum(Sdir+r.Dsky*svf)
efctvalb = 1-sum(Knet_tile)/sum(Sdirtile+r2.Dsky*svftile);
preprocessing.write_efalb(r,efctvalb);
toc
%%
function newmatrix = diagonal_matrix_maker(m, nxtiles, nytiles)
    % Determine the size of the matrix m
    [rowm, colm] = size(m);
    % Specify the desired size of matrix B (m*m)
    tiles = nxtiles*nytiles;
    totx = rowm*tiles;
    toty = colm*tiles;
    B = zeros(totx,toty);
    for i = 0:tiles-1
            il = 1+i*rowm;
            iu = rowm+i*rowm;
            jl = 1+i*colm;
            ju = colm+i*colm;
            B(il:iu,jl:ju) = m;
    end
    % Display the resulting matrix B
    newmatrix = B;
end 