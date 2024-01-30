%% write_inputs

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

expnr = '901';
%
DA_EXPDIR = getenv('DA_EXPDIR');
DA_TOOLSDIR = getenv('DA_TOOLSDIR');
addpath(genpath([DA_TOOLSDIR '/']));
addpath([DA_TOOLSDIR '/IBM/'])
addpath([DA_TOOLSDIR '/SEB/'])
exppath = [DA_EXPDIR '/'];
fpath = [DA_EXPDIR '/' expnr '/'];
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

if r.nsv>0
    preprocessing.generate_scalar(r);
    preprocessing.write_scalar(r);
    disp(['Written scalar.inp.', r.expnr])
end

if r.ltrees || r.ltreesfile
    disp('Generating trees')
    preprocessing.generate_trees_from_namoptions(r);
    preprocessing.write_trees(r);
    disp(['Written trees.inp.', r.expnr])
end

if isfile(['factypes.inp.', expnr])
    r.factypes = dlmread(['factypes.inp.', r.expnr],'',3,0);
else
    preprocessing.write_factypes(r)
    disp(['Written factypes.inp.', r.expnr])
end

if r.libm
    %% Read the .stl file and write necessary ibm files
    TR = stlread(r.stl_file);
    F = TR.ConnectivityList;
    V = TR.Points;
    area_facets = facetAreas(F, V); % Useful for checking if area_fluid_IB_c == sum(area_facets)

    % Set facet types
    nfcts = size(TR.ConnectivityList,1);
    preprocessing.set_nfcts(r, nfcts);
    if r.read_types
        facet_types = dlmread(r.types_path, '', 1, 0);
    else
        facet_types = ones(nfcts,1); % defaults to type 1
    end
    preprocessing.write_facets(r, facet_types, TR.faceNormal);
    disp(['Written facets.inp.', r.expnr])

    calculate_facet_sections_uvw = r.iwallmom > 1;
    calculate_facet_sections_c = r.ltempeq || r.lmoist;
    lwindows = false;
    if r.gen_geom
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
        itot = r.itot;
        jtot = r.jtot;
        ktot = r.ktot;
        dx = r.dx;
        dy = r.dy;

        lmypolyfortran = 1; lmypoly = 0;		% remove eventually
        lmatchFacetsToCellsFortran = 1;

        writeIBMFiles; % Could turn into a function and move writing to this script
    else
        if isempty(r.geom_path)
            error('Need to specify the path to geometry files')
        end
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

    preprocessing.write_facetarea(r, area_facets); % always write facet area
    if r.lEB
        %% Write STL in View3D input format
        fpath_facets_view3d = [fpath 'facets.vs3'];
        STLtoView3D(r.stl_file, fpath_facets_view3d, r.view3d_out, r.maxD);

        %% Calculate view factors
        % remember to build View3D in local system windows/linux
        % Add check to see if View3D exists in the tools directory.
        if lwindows
            view3d_exe = [DA_TOOLSDIR '/View3D/src/View3D.exe'];
        else
            view3d_exe = [DA_TOOLSDIR '/View3D/build/src/view3d'];
        end

        %vf = view3d(view3d_exe, fpath_facets_view3d, fpath_vf);
        if r.calc_vf % run view3d
            if r.view3d_out == 0 % text
                fpath_vf = [fpath 'vf.txt'];
            elseif r.view3d_out == 1 % binary
                fpath_vf = [fpath 'vf.bin'];
            elseif r.view3d_out == 2 % sparse
                fpath_vf = [fpath 'vfsparse.inp.' r.expnr];
            end

            view3d_execution_command = [view3d_exe ' ' fpath_facets_view3d ' ' fpath_vf];
            system(view3d_execution_command);
        else % view3d has already been run - read output file
            fpath_vf = r.vf_path;
        end

        if r.view3d_out == 0
            vf = dlmread(fpath_vf, ' ', 2, 0);
            vf(end,:) = [];
            vf = sparse(vf);
        elseif r.view3d_out == 1
            fID = fopen(fpath_vf);
            area = fread(fID, 8+r.nfcts, 'single'); % first 8 bytes are header
            vf = fread(fID, r.nfcts^2, 'single');
            fclose(fID);
            %vf = reshape(vf, [nfcts, nfcts])'; % transpose to get in row-major order
            %vf = sparse(vf);
            %svf = max(1 - sum(vf, 2), 0);
            % since doing above sometimes exceeds array size...
            vf = reshape(vf, [r.nfcts, r.nfcts]);
            vf = sparse(vf);
            vf = vf';
        elseif r.view3d_out == 2
            ijs = dlmread(fpath_vf, ' ', 0, 0);
            vf = sparse(ijs(:,1), ijs(:,2), ijs(:,3));
        end

        svf = max(1 - sum(vf, 2), 0);
        preprocessing.write_svf(r, svf);

        % write uDALES view factor file
        if (r.view3d_out == 0 || r.view3d_out == 1)  % view3d_out==2 is already sparse
            if ~r.lvfsparse
                preprocessing.write_vf(r, vf)
                disp(['Written vf.nc.inp.', r.expnr])
            else
                %vfsparse = sparse(double(vf));
                preprocessing.write_vfsparse(r, vf);
                disp(['Written vfsparse.inp.', r.expnr])
            end
            if r.calc_vf
                delete(fpath_vf) % remove view3d output file
            end
        elseif (r.view3d_out == 2 && ~r.calc_vf)
            copyfile(fpath_vf, [fpath 'vfsparse.inp.' r.expnr]);
        end

        %% Calculate shortwave radiation
        albedos = preprocessing.generate_albedos(r, facet_types);
        resolution   = r.psc_res;
        xazimuth     = r.xazimuth;
        ltimedepsw   = r.ltimedepsw;
        ldirectShortwaveFortran = r.ldirectShortwaveFortran;
        lscatter = true;

        if ltimedepsw
            runtime = r.runtime;
            dtSP    = r.dtSP;
            start = datetime(r.year, r.month, r.day, r.hour, r.minute, r.second);
            longitude = r.longitude;
            latitude  = r.latitude;
            timezone  = r.timezone;
            elevation = r.elevation;
        else
            lcustomsw = r.lcustomsw;
            if lcustomsw
                solarazimuth = r.solarazimuth;
                solarzenith  = r.solarzenith;
                irradiance   = r.I;
                Dsky         = r.Dsky;
            else
                start = datetime(r.year, r.month, r.day, r.hour, r.minute, r.second);
                longitude = r.longitude;
                latitude  = r.latitude;
                timezone  = r.timezone;
                elevation = r.elevation;
            end
        end

        shortwave;
        preprocessing.write_netsw(r, Knet(:,1));
        disp(['Written netsw.inp.', r.expnr])

        if r.ltimedepsw
            preprocessing.write_timedepsw(r, tSP, Knet)
        end
    end

    %% Write initial facet temperatures
    if (r.lEB || r.iwallmom == 2 || r.iwalltemp == 2)
        disp('Setting initial facet temperatures.')
        facT = r.facT;
        nfaclyrs = r.nfaclyrs;
        facT_file = r.facT_file;
        lfacTlyrs = r.lfacTlyrs;
        if ~r.lfacTlyrs
            Tfacinit = ones(nfcts,1) .* r.facT;
            preprocessing.write_Tfacinit(r, Tfacinit)
            disp(['Written Tfacinit.inp.', r.expnr])
            % Could always read in facet temperature as layers, defaulting to linear?
        else
            Tfac = ncread(r.facT_file, 'T');
            Tfacinit_layers = Tfac(:, :, end);
            preprocessing.write_Tfacinit_layers(r, Tfacinit_layers)
            disp(['Written Tfacinit_layers.inp.', r.expnr])
        end

    end
end
