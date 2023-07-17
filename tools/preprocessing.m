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

% Copyright (C) 2016-2021 the uDALES Team.

classdef preprocessing < dynamicprops
    % Class for pre-processing in uDALES
    %
    properties (Hidden = true, SetAccess = protected)
        path;                    % Path to simulations.
        cpath;                   % Current path.
        buf;                     % All purpose buffer.
        expnr;
    end

    methods (Static)

        function obj = preprocessing(expnr, varargin)
            % Class constructor.
            %
            % expnr       : Integer equal to simulation number.
            % varargin{2} : Path to simulations.
            % varargin includes the explicitly declared input, thus
            % varargin{1}=expnr

            %read namoptions

            TOKENS = '(.*)\=(.*)';
            DROP = '(\s*|\[|\]|''|''|;)';
            WHITE = '\s*';
            VOID = '';

            obj.cpath = pwd;
            if (nargin==0)
                disp('expnr input not found. . Exiting ...')
                return
            elseif nargin > 1
                dapath = [varargin{1},num2str(expnr)];
            else
                dapath = num2str(varargin{1});
            end
            cd(dapath);
            obj.path = pwd;

            expstr = num2str(expnr, '%3.3d');
            fid = fopen(['namoptions.', expstr]);

            if (fid == -1)
                disp(['namoptions.', expstr, ' not found. Exiting ...'])
                return
            end

            obj.expnr = expstr;

            line = fgets(fid);

            % load each line of simparams_* into object
            while min(line) > 0
                toks = regexp(line, TOKENS, 'tokens');
                if ~isempty(toks)
                    lhs = regexprep(toks{1}{1}, WHITE, VOID);
                    h = addprop(obj, lhs);
                    rhs=regexprep(toks{1}{2}, WHITE, VOID);
                    if strcmp(rhs,'.false.')
                        obj.(lhs) = 0;
                    elseif strcmp(rhs,'.true.')
                        obj.(lhs) = 1;
                    elseif all(isstrprop(rhs,'digit'))  %check if it is a number
                        obj.(lhs) = str2double(rhs);
                    else
                        if ~isnan(str2double(rhs))
                            % number or logical
                            obj.(lhs) = str2double(rhs);
                        else
                            % filename
                            obj.(lhs) = rhs;
                        end
                    end
                    h.SetAccess = 'protected';
                end
                line = fgets(fid);
            end

            cd(obj.cpath); % tg3315 23.04.18 - so does not change directory everytime

        end

        function addvar(obj, lhs, var)
            % Install a variable in the object.
            % addvar(obj, lhs, var)
            %
            % lhs : Identifier for variable.
            % var : R-value (value of variable).

            if (~isprop(obj, lhs))
                obj.addprop(lhs);
                %                 h.SetAccess = 'public';
                obj.(lhs) = var;
            end

        end

        function gopath(obj)
            % Go to simulation path.
            % gopath(obj)

            cd(obj.path);
        end

        function gohome(obj)
            % Go to work path.
            % gohome(obj)

            cd(obj.cpath);
        end

        function chcpath(obj, newpath)

            % Change work path.
            % chcpath(obj, newpath)

            here = pwd;
            cd(newpath)
            obj.cpath = pwd;
            cd(here);
        end

        function set_defaults(obj)
            %% &RUN
            preprocessing.addvar(obj, 'ltrees', 0) % switch for trees (not implemented)
            if obj.ltrees
                error('Trees not currently implemented')
%                 preprocessing.addvar(obj, 'tree_dz',0)
%                 preprocessing.addvar(obj, 'tree_dx',0)
%                 preprocessing.addvar(obj, 'tree_h',0)
%                 preprocessing.addvar(obj, 'tree_w',0)
%                 preprocessing.addvar(obj, 'tree_b',0)
%
%                 preprocessing.addvar(obj, 'nt1',0)
%                 preprocessing.addvar(obj, 'md',0)
%                 preprocessing.addvar(obj, 'ww',0)
%                 preprocessing.addvar(obj, 'lw',0)
%                 preprocessing.addvar(obj, 'nt2',0)
            end

            preprocessing.addvar(obj, 'lpurif', 0) % switch for purifiers (not implemented)
            if obj.lpurif
                error('Purifiers not currently implemented')
%                 if obj.lcanyons
%                     preprocessing.addvar(obj, 'purif_dz', 1)  % purifier starting point from bottom
%                     preprocessing.addvar(obj, 'purif_dx', 3)  % distance from block
%                     preprocessing.addvar(obj, 'purif_h', 3)   % purifier height
%                     preprocessing.addvar(obj, 'purif_w', 0)   % purifier width
%                     preprocessing.addvar(obj, 'purif_dy', 1)  % depth of purifier (in y)
%                     preprocessing.addvar(obj, 'purif_sp', 31) % spacing between purifiers
%                     preprocessing.addvar(obj, 'purif_i', 1)   % case for purifier (1 = +ve x, 2 = -ve x, 3 = +ve y etc.)
%                     preprocessing.addvar(obj, 'npurif', obj.jtot / (obj.npurif_dy + obj.purif_sp))
%
%                     if ceil(npurif) ~= floor(npurif)
%                         lp = 0:obj.tot / 2;
%                         indp = rem(obj.jtot / 2, lp) == 0;
%                         errp = ([lp(indp); (obj.jtot / 2) ./ lp(indp)]);
%                         disp('Purifier layout does not fit grid')
%                         disp(['sum widths to: ' num2str(errp(1,:))])
%                         disp(['Current width: ' num2str(obj.purif_dy + obj.purif_sp)])
%                         error('Incorrect purifier layout')
%                     end
%                 else
%                     error('Must use lcanyons configuration to use purifiers')
%                 end
            end

            preprocessing.addvar(obj, 'luoutflowr', 0) % switch that determines whether u-velocity is corrected to get a fixed outflow rate
            preprocessing.addvar(obj, 'lvoutflowr', 0) % switch that determines whether v-velocity is corrected to get a fixed outflow rate.
            preprocessing.addvar(obj, 'luvolflowr', 0) % switch that determines whether u-velocity is corrected to get a fixed volume flow rate.
            preprocessing.addvar(obj, 'lvvolflowr', 0) % switch that determines whether v-velocity is corrected to get a fixed volume flow rate.

            %% &DOMAIN
            preprocessing.addvar(obj, 'itot', 64)  % # cells in x-direction
            preprocessing.addvar(obj, 'xlen', 64) % domain size in x-direction
            preprocessing.addvar(obj, 'jtot', 64)  % # cells in y-direction
            preprocessing.addvar(obj, 'ylen', 64) % domain size in y-direction
            preprocessing.addvar(obj, 'ktot', 96)  % # cells in z-direction

            preprocessing.addvar(obj, 'dx', obj.xlen / obj.itot)
            preprocessing.addvar(obj, 'dy', obj.ylen / obj.jtot)

            %% BCs

            preprocessing.addvar(obj, 'BCxm', 1);
            preprocessing.addvar(obj, 'BCym', 1);


            %% &ENERGYBALANCE
            preprocessing.addvar(obj, 'lEB', 0)
            preprocessing.addvar(obj, 'lfacTlyrs', 0)
            

            %% &WALLS
            preprocessing.addvar(obj, 'iwallmom', 3)
            preprocessing.addvar(obj, 'iwalltemp', 1)

            %% &PHYSICS
            preprocessing.addvar(obj, 'ltempeq', 0)
            preprocessing.addvar(obj, 'lmoist', 0)
            preprocessing.addvar(obj, 'lchem' , 0) % switch for chemistry (not implemented)
            preprocessing.addvar(obj, 'lprofforc', 0)  % switch for 1D geostrophic forcing
            preprocessing.addvar(obj, 'lcoriol', 0)    % switch for coriolis forcing
            preprocessing.addvar(obj, 'idriver', 0)    % case for driver simulations | 1 - writes driver files | 2 - reads driver files

            if (not(obj.luoutflowr) && not(obj.lvoutflowr) && not(obj.luvolflowr) && not(obj.lvvolflowr) && not(obj.lprofforc) && not(obj.lcoriol) && (obj.idriver~=2))
                preprocessing.addvar(obj, 'ldp', 1)
                disp('No forcing switch config. setup and not a driven simulation so initial velocities and pressure gradients applied.')
            else
                preprocessing.addvar(obj, 'ldp', 0)
            end

            if (obj.ltempeq == 0 || obj.iwalltemp==1 && obj.iwallmom==2)
                obj.iwallmom = 3;
            end

            %% &INPS
            preprocessing.addvar(obj, 'zsize', 96) % domain size in z-direction
            preprocessing.addvar(obj, 'lzstretch', 0) % switch for stretching z grid
            preprocessing.addvar(obj, 'stl_file', '')
            preprocessing.addvar(obj, 'gen_geom', true) % generate the geometry from scratch
            preprocessing.addvar(obj, 'geom_path', '') % if not generating the geometry, the path to the geometry files
            preprocessing.addvar(obj, 'diag_neighbs', true)
            preprocessing.addvar(obj, 'stl_ground', true) % Does STL include facets at ground 
            
            if obj.lzstretch
                preprocessing.addvar(obj, 'stretchconst', 0.01)
                preprocessing.addvar(obj, 'lstretchexp', 0)
                preprocessing.addvar(obj, 'lstretchtanh', 0)
                preprocessing.addvar(obj, 'lstretch2tanh', 0)
                preprocessing.addvar(obj, 'hlin', 0)
                preprocessing.addvar(obj, 'dzlin', 0)
                preprocessing.addvar(obj, 'dz', obj.dzlin)
            else
                preprocessing.addvar(obj, 'dz', obj.zsize / obj.ktot)
            end

            if obj.lEB
                preprocessing.addvar(obj, 'maxlen', 10); % maximum size of facets
            else
                preprocessing.addvar(obj, 'maxlen', inf);
            end

            preprocessing.addvar(obj, 'u0', 0) % initial u-velocity - also applied as geostrophic term where applicable
            preprocessing.addvar(obj, 'v0', 0) % initial v-velocity - also applied as geostrophic term where applicable
            preprocessing.addvar(obj, 'tke', 0)
            preprocessing.addvar(obj, 'dpdx', 0) % dp/dx [Pa/m]
            preprocessing.addvar(obj, 'dpdy', 0) % dp/dy [Pa/m]
            preprocessing.addvar(obj, 'thl0', 288) % temperature at lowest level
            preprocessing.addvar(obj, 'qt0', 0)    % specific humidity

            preprocessing.addvar(obj, 'nsv', 0)         % number of scalar variables (not implemented)
            if obj.nsv>0
                preprocessing.addvar(obj, 'sv10', 0)        % first scalar variable initial/ background conc.
                preprocessing.addvar(obj, 'sv20', 0)        % second scalar variable initial/ background conc.
                preprocessing.addvar(obj, 'sv30', 0)        % third scalar variable initial/ background conc.
                preprocessing.addvar(obj, 'sv40', 0)        % fourth scalar variable initial/ background conc.
            	preprocessing.addvar(obj, 'sv50', 0)        % fifth scalar variable initial/ background conc.
            	preprocessing.addvar(obj, 'lscasrc', 0)     % switch for scalar point source
            	preprocessing.addvar(obj, 'lscasrcl', 0)    % switch for scalar line source
            	preprocessing.addvar(obj, 'lscasrcr', 0)    % switch for network of scalar point source
            	preprocessing.addvar(obj, 'xS', -1)         % x-position of scalar point source [m]
            	preprocessing.addvar(obj, 'yS', -1)         % y-position of scalar point source [m]
            	preprocessing.addvar(obj, 'xS', -1)         % z-position of scalar point source [m]
           	 preprocessing.addvar(obj, 'SS', -1)         % source strength of scalar line/ point source
            	preprocessing.addvar(obj, 'sigS', -1)       % standard deviation/ spread of scalar line/ point source
            	if ((obj.lscasrc) && any([obj.xS==-1 obj.yS==-1 obj.zS==-1 obj.SS==-1 obj.sigS==-1]))
                    error('Must set non-zero xS, yS, zS, SS and sigS for scalar point source')
            	end
            	if ((obj.lscasrcl) && any([obj.SS==-1 obj.sigS==-1]))
                    error('Must set non-zero SS and sigS for scalar line source')
            	end
            	if obj.lscasrcr
                    error('Network of point sources not currently implemented')
            	end
            end

            preprocessing.addvar(obj, 'lapse', 0)  % lapse rate [K/s]
            preprocessing.addvar(obj, 'w_s',0) % subsidence [*units?*]
            preprocessing.addvar(obj, 'R',0)   % radiative forcing [*units?*]

            preprocessing.addvar(obj, 'libm', 1)

            if obj.lEB
                preprocessing.addvar(obj, 'solarazimuth', 135); % azimuth angle
                preprocessing.addvar(obj, 'xazimuth', 90);   % azimuth of x-direction wrt N. Default: x = East
                                                             % north -> xazimuth = 0;
                                                             % east  ->            90;
                                                             % south ->            180;
                                                             % west  ->            270;
                preprocessing.addvar(obj, 'solarzenith', 28.4066); % zenith angle
                preprocessing.addvar(obj, 'centerweight', 12 / 32);
                preprocessing.addvar(obj, 'cornerweight', (1 - obj.centerweight) / 4);
                preprocessing.addvar(obj, 'I', 800); % Direct normal irradiance [W/m2]
                preprocessing.addvar(obj, 'Dsky', 418.8041); % Diffuse incoming radiation [W/m2]
                preprocessing.addvar(obj, 'lvfsparse', false) % Switch for turning on lvfsparse
                preprocessing.addvar(obj, 'psc_res', 0.01); % Poly scan conversion resolution,lower number gives better results for solar radiation calculation 
                %preprocessing.addvar(obj, 'min_vf', 0.01); % Any vf below this is ignored in sparse format
            end
            preprocessing.addvar(obj, 'facT', 288.) % Initial facet temperatures.
            preprocessing.addvar(obj, 'nfaclyrs', 3) % Number of facet layers
            
            preprocessing.addvar(obj, 'nfcts', 0)

            preprocessing.generate_factypes(obj)
            preprocessing.addvar(obj, 'facT_file', '')

        end

        function generate_factypes(obj)
            K = obj.nfaclyrs;
            factypes = [];

            % Bounding walls (bw)
            id_bw  = -101;
            lGR_bw = 0;
            z0_bw  = 0;
            z0h_bw = 0;
            al_bw  = 0.5;
            em_bw  = 0.85;
            D_bw   = 0.;
            d_bw   = D_bw / K;
            C_bw   = 0.;
            l_bw   = 0.;
            k_bw   = 0.;
            bw = [id_bw, lGR_bw, z0_bw, z0h_bw, al_bw, em_bw, d_bw * ones(1,K), C_bw * ones(1,K), l_bw * ones(1,K), k_bw * ones(1,K+1)];
            factypes = [factypes; bw];

            % Floors (f)
            id_f  = -1;
            lGR_f = 0;
            z0_f  = 0.05;
            z0h_f = 0.00035;
            al_f  = 0.5;
            em_f  = 0.85;
            D_f   = 0.5;
            d_f   = D_f / K;
            C_f   = 1.875e6;
            l_f   = 0.75;
            k_f   = 0.4e-6;
            if (K == 3)
                % Reproduce the original factypes.inp (d_f not constant for each layer)
                f =  [id_f, lGR_f, z0_f, z0h_f, al_f, em_f, 0.1, 0.2, 0.2, C_f * ones(1,K), l_f * ones(1,K), k_f * ones(1,K+1)];
            else
                f =  [id_f, lGR_f, z0_f, z0h_f, al_f, em_f, d_f * ones(1,K), C_f * ones(1,K), l_f * ones(1,K), k_f * ones(1,K+1)];
            end
            factypes = [factypes; f];

            % Dummy (dm)
            id_dm  = 0;
            lGR_dm = 0;
            z0_dm  = 0;
            z0h_dm = 0;
            al_dm  = 0;
            em_dm  = 0;
            D_dm   = 0.3;
            d_dm = D_dm / K;
            C_dm = 1.875e6;
            l_dm = 0.75;
            k_dm = 0.4e-6;
            dm = [id_dm, lGR_dm, z0_dm, z0h_dm, al_dm, em_dm, d_dm * ones(1,K), C_dm * ones(1,K), l_dm * ones(1,K), k_dm * ones(1,K+1)];
            factypes = [factypes; dm];

            % Concrete (c)
            id_c  = 1;
            lGR_c = 0;
            z0_c  = 0.05;
            z0h_c = 0.00035;
            al_c = 0.5;
            em_c = 0.85;
            D_c = 0.36;
            d_c = D_c / K;
            C_c = 2.5e6;
            l_c = 1;
            k_c = 0.4e-6;
            c = [id_c, lGR_c, z0_c, z0h_c, al_c, em_c, d_c * ones(1,K), C_c * ones(1,K), l_c * ones(1,K), k_c * ones(1,K+1)];
            factypes = [factypes; c];

            % Brick (b)
            id_b  = 2;
            lGR_b = 0;
            z0_b  = 0.05;
            z0h_b = 0.00035;
            al_b = 0.5;
            em_b = 0.85;
            D_b = 0.36;
            d_b = D_b / K;
            C_b = 2.766667e6;
            l_b = 0.83;
            k_b = 0.3e-6;
            b = [id_b, lGR_b, z0_b, z0h_b, al_b, em_b, d_b * ones(1,K), C_b * ones(1,K), l_b * ones(1,K), k_b * ones(1,K+1)];
            factypes = [factypes; b];

            % Stone (s)
            id_s  = 3;
            lGR_s = 0;
            z0_s  = 0.05;
            z0h_s = 0.00035;
            al_s = 0.5;
            em_s = 0.85;
            D_s = 0.36;
            d_s = D_s / K;
            C_s = 2.19e6;
            l_s = 2.19;
            k_s = 1e-6;
            s = [id_s, lGR_s, z0_s, z0h_s, al_s, em_s, d_s * ones(1,K), C_s * ones(1,K), l_s * ones(1,K), k_s * ones(1,K+1)];
            factypes = [factypes; s];

            % Wood (w)
            id_w  = 4;
            lGR_w = 0;
            z0_w  = 0.05;
            z0h_w = 0.00035;
            al_w = 0.5;
            em_w = 0.85;
            D_w = 0.36;
            d_w = D_w / K;
            C_w = 1e6;
            l_w = 0.1;
            k_w = 0.1e-6;
            w = [id_w, lGR_w, z0_w, z0h_w, al_w, em_w, d_w * ones(1,K), C_w * ones(1,K), l_w * ones(1,K), k_w * ones(1,K+1)];
            factypes = [factypes; w];

            % GR1
            id_GR1 = 11;
            lGR_GR1 = 1;
            z0_GR1 = 0.05;
            z0h_GR1 = 0.00035;
            al_GR1 = 0.25;
            em_GR1 = 0.95;
            D_GR1 = 0.6;
            d_GR1 = D_GR1 / K;
            C_GR1 = 5e6;
            l_GR1 = 2;
            k_GR1 = 0.4e-6;
            GR1 = [id_GR1, lGR_GR1, z0_GR1, z0h_GR1, al_GR1, em_GR1, d_GR1 * ones(1,K), C_GR1 * ones(1,K), l_GR1 * ones(1,K), k_GR1 * ones(1,K+1)];
            factypes = [factypes; GR1];

            % GR2
            id_GR2 = 12;
            lGR_GR2 = 1;
            z0_GR2 = 0.05;
            z0h_GR2 = 0.00035;
            al_GR2 = 0.35;
            em_GR2 = 0.90;
            D_GR2 = 0.6;
            d_GR2 = D_GR2 / K;
            C_GR2 = 2e6;
            l_GR2 = 0.8;
            k_GR2 = 0.4e-6;
            GR2 = [id_GR2, lGR_GR2, z0_GR2, z0h_GR2, al_GR2, em_GR2, d_GR2 * ones(1,K), C_GR2 * ones(1,K), l_GR2 * ones(1,K), k_GR2 * ones(1,K+1)];
            factypes = [factypes; GR2];

            preprocessing.addvar(obj, 'factypes', factypes);
        end

        function write_facets(obj, types, normals)
            fname = ['facets.inp.' obj.expnr];
            fileID = fopen(fname,'W');
            fprintf(fileID, '# type, normal\n');
            fprintf(fileID, '%-4d %-4.4f %-4.4f %-4.4f\n', [types normals]');
            fclose(fileID);
        end

        function write_factypes(obj)
            K = obj.nfaclyrs;

            fname = ['factypes.inp.', obj.expnr];

            dheaderstring = '';
            for k=1:K
                dheaderstring = [dheaderstring, sprintf('  d%d [m]',k)];
            end

            Cheaderstring = '';
            for k=1:K
                Cheaderstring = [Cheaderstring, sprintf('  C%d [J/(K m^3)]', k)];
            end

            lheaderstring = '';
            for k=1:K
                lheaderstring = [lheaderstring, sprintf('  l%d [W/(m K)]', k)];
            end

            kheaderstring = '';
            for k=1:K+1
                kheaderstring = [kheaderstring, sprintf('  k%d [W/(m K)]', k)];
            end

            fileID = fopen(fname,'W');
            fprintf(fileID, ['# walltype, ', num2str(K), ' layers per type where layer 1 is the outdoor side and layer ', num2str(K), ' is indoor side\n']);
            fprintf(fileID, '# 0=default dummy, -1=asphalt floors;-101=concrete bounding walls;1=concrete;2=bricks;3=stone;4=painted wood;11=GR1; 12=GR2\n');
            fprintf(fileID, ['# wallid  lGR  z0 [m]  z0h [m]  al [-]  em [-]', dheaderstring, Cheaderstring, lheaderstring, kheaderstring, '\n']);

            valstring1 = '%8d  %3d  %6.2f  %7.5f  %6.2f  %6.2f';

            valstring2 = '';
            for k=1:K
                valstring2 = [valstring2, '  %6.2f'];
            end
            for k=1:K
                valstring2 = [valstring2, '  %14.0f'];
            end
            for k=1:K
                valstring2 = [valstring2, ' %13.4f'];
            end
            for k=1:K+1
                valstring2 = [valstring2, ' %13.8f'];
            end

            valstring = [valstring1, valstring2];

            [nfactypes, ~] = size(obj.factypes);

            for i = 1:nfactypes
                fprintf(fileID, sprintf(valstring, obj.factypes(i,:)));
                fprintf(fileID, '\n');
            end

            fclose(fileID);
        end

        function albedos = generate_albedos(obj, facet_types)
            albedos = [];
            typeids = obj.factypes(:,1);
            for i = 1:obj.nfcts
                my_typid = facet_types(i);
                albedo = obj.factypes(find(typeids == my_typid),5);
                albedos = [albedos; albedo];
            end
        end

        function plot_profiles(obj)
            figure
            subplot(141)
            plot(obj.pr(:, 2), 1:obj.ktot)
            title('Temperature')

            subplot(142)
            plot(obj.ls(:, 10), 1:obj.ktot)
            title('Radiative forcing')

            subplot(143)
            plot(obj.ls(:, 6), 1:obj.ktot)
            title('Subsidence')

            subplot(144)
            plot(obj.ls(:,2), 1:obj.ktot)
            hold on
            plot(obj.ls(:,3), 1:obj.ktot, 'r--')
            title('Velocity')
            legend('u', 'v')
        end

        function generate_xygrid(obj)
             preprocessing.addvar(obj, 'xf', 0.5 * obj.dx : obj.dx : obj.xlen - 0.5 * obj.dx);
             preprocessing.addvar(obj, 'yf', 0.5 * obj.dy : obj.dy : obj.ylen - 0.5 * obj.dy);
             preprocessing.addvar(obj, 'xh', 0 : obj.dx : obj.xlen);
             preprocessing.addvar(obj, 'yh', 0 : obj.dy : obj.ylen);
        end

        function write_xgrid(obj)
            xgrid = fopen(['xgrid.inp.' obj.expnr], 'w');
            fprintf(xgrid, '%12s\n', '#     x-grid');
            fprintf(xgrid, '%12s\n', '#           ');
            fprintf(xgrid, '%-20.15f\n', obj.xf);
            fclose(xgrid);
        end

        function generate_zgrid(obj)
            if ~obj.lzstretch
                preprocessing.addvar(obj, 'zf', 0.5 * obj.dz:obj.dz:obj.zsize - 0.5 * obj.dz);
                preprocessing.addvar(obj, 'zh', 0:obj.dz:obj.zsize);
                preprocessing.addvar(obj, 'dzf', obj.zh(2:end) - obj.zh(1:end - 1));
            else
                if obj.lstretchexp
                   preprocessing.stretch_exp(obj)
                elseif obj.lstretchtanh
                    preprocessing.stretch_tanh(obj)
                elseif obj.lstretch2tanh
                    preprocessing.stretch_2tanh(obj)
                else
                    error('Invalid stretch');
                end
            end
        end

        function stretch_exp(obj)
            il = round(obj.hlin / obj.dzlin);
            ir  = obj.ktot - il;

            preprocessing.addvar(obj, 'zf', zeros(1, obj.ktot));
            preprocessing.addvar(obj, 'dzf', zeros(1, obj.ktot));
            preprocessing.addvar(obj, 'zh', zeros(1, obj.ktot+1));

            obj.zf(1:il) = 0.5 * obj.dzlin : obj.dzlin : obj.hlin;
            obj.zh(1:il+1) = 0 : obj.dzlin : obj.hlin;

            gf = obj.stretchconst;

            while true
                obj.zh(il + 1:end) = obj.zh(il + 1) + (obj.zsize - obj.zh(il+1)) * (exp(gf * (0:1:ir) / (ir)) - 1)/(exp(gf) - 1); %dh has been replaced by zsize
                if (obj.zh(il+2) - obj.zh(il + 1)) < obj.dzlin
                    gf = gf - 0.01; %make sufficiently small steps to avoid an initial bump in dz
                else
                    if (obj.zh(end) - obj.zh(end - 1)) > 3 * obj.dzlin
                        disp('Warnning: final grid spacing large - consider reducing domain height')
                    end
                    break
                end
            end

            for i = 1:obj.ktot
                obj.zf(i) = (obj.zh(i) + obj.zh(i+1)) / 2 ;
                obj.dzf(i) = obj.zh(i+1) - obj.zh(i);
            end
        end

        function stretch_tanh(obj)
            il = round(obj.hlin / obj.dzlin);
            ir  = obj.ktot - il;

            preprocessing.addvar(obj, 'zf', zeros(obj.ktot, 1));
            preprocessing.addvar(obj, 'dzf', zeros(obj.ktot, 1));
            preprocessing.addvar(obj, 'zh', zeros(obj.ktot + 1, 1));

            obj.zf(1:il) = 0.5 * obj.dzlin : obj.dzlin : obj.hlin;
            obj.zh(1:il+1) = 0 : obj.dzlin : obj.hlin;

            gf = obj.stretchconst;

            while true
                obj.zh(il + 1:end) = obj.zh(il + 1) + (obj.zsize - obj.zh(il + 1)) * (1 - tanh(gf * (1 - 2 * (0:1:ir)' / (2*ir))) / tanh(gf));

            if (obj.zh(il + 2) - obj.zh(il + 1)) < obj.dzlin
                gf = gf - 0.01; % make sufficiently small steps to avoid an initial bump in dz
            else
                if (obj.zh(end) - obj.zh(end - 1)) > 3 * obj.dzlin
                disp('Warning: final grid spacing large - consider reducing domain height')
                end
                break
            end

            end

            for i = 1:obj.ktot
                obj.zf(i) = 0.5 * (obj.zh(i) + obj.zh(i+1));
                obj.dzf(i) = obj.zh(i+1) - obj.zh(i);
            end
        end

        function stretch_2tanh(obj)
            il = round(obj.hlin / obj.dzlin);
            ir  = obj.ktot - il;

            preprocessing.addvar(obj, 'zf', zeros(obj.ktot, 1));
            preprocessing.addvar(obj, 'dzf', zeros(obj.ktot, 1));
            preprocessing.addvar(obj, 'zh', zeros(obj.ktot+1, 1));

            obj.zf(1:il) = 0.5 * obj.dzlin:obj.dzlin:obj.hlin;
            obj.zh(1:il+1) = 0:obj.dzlin:obj.hlin;
            gf = obj.stretchconst;

            while true
                obj.zh(il+1:end) = obj.zh(il+1) + (obj.zsize - obj.zh(il+1)) / 2 * (1 - tanh(gf * (1 - 2 * (0:1:ir)'/(ir))) / tanh(gf));
                if (obj.zh(il + 2) - obj.zh(il + 1)) < obj.dzlin
                    gf = gf - 0.01; % make sufficiently small steps to avoid an initial bump in dz
                else
                    if (max(diff(obj.zh))) > 3 * obj.dzlin
                        disp('Warning: final grid spacing large - consider reducing domain height')
                    end
                    break
                end
            end

            for i = 1:obj.ktot
                obj.zf(i) = (obj.zh(i) + obj.zh(i+1)) / 2 ;
                obj.dzf(i) = obj.zh(i + 1) - obj.zh(i);
            end
        end

        function write_zgrid(obj)
            zgrid = fopen(['zgrid.inp.' obj.expnr], 'w');
            fprintf(zgrid, '%12s\n', '#     z-grid');
            fprintf(zgrid, '%12s\n', '#           ');
            fprintf(zgrid, '%-20.15f\n', obj.zf);
            fclose(zgrid);
        end

        function generate_lscale(obj)
            if ((obj.luoutflowr || obj.lvoutflowr) + (obj.luvolflowr || obj.lvvolflowr) + obj.lprofforc + obj.lcoriol + obj.ldp) > 1
                error('More than one forcing type specified')
            end

            preprocessing.addvar(obj, 'ls', zeros(length(obj.zf), 10));
            obj.ls(:,1) = obj.zf;
            obj.ls(:,6) = obj.w_s;
            obj.ls(:,10) = obj.R;
            if obj.lprofforc || obj.lcoriol
                obj.ls(:,2) = obj.u0;
                obj.ls(:,3) = obj.v0;
            elseif obj.ldp
                obj.ls(:,4) = obj.dpdx;
                obj.ls(:,5) = obj.dpdy;
            end
        end

        function write_lscale(obj)
            lscale = fopen(['lscale.inp.' obj.expnr], 'w');
            fprintf(lscale, '%-12s\n', '# SDBL flow');
            fprintf(lscale, '%-60s\n', '# z uq vq pqx pqy wfls dqtdxls dqtdyls dqtdtls dthlrad');
            fprintf(lscale, '%-20.15f %-12.6f %-12.6f %-12.9f %-12.6f %-15.9f %-12.6f %-12.6f %-12.6f %-17.12f\n', obj.ls');
            fclose(lscale);
        end

        function generate_prof(obj)
            preprocessing.addvar(obj, 'pr', zeros(length(obj.zf), 6));
            obj.pr(:,1) = obj.zf;

            if obj.lapse
                thl = zeros(obj.ktot, 1);
                thl(1) = obj.thl0;
                for k = 1:obj.ktot - 1
                    thl(k+1) = thl(k) + obj.lapse * obj.zsize / obj.ktot;
                end
                obj.pr(:,2) = thl;
            else
                obj.pr(:,2) = obj.thl0;
            end

            obj.pr(:,3) = obj.qt0;
            obj.pr(:,4) = obj.u0;
            obj.pr(:,5) = obj.v0;
            obj.pr(:,6) = obj.tke;
        end

        function write_prof(obj)
            prof = fopen(['prof.inp.' obj.expnr], 'w');
            fprintf(prof, '%-12s\n', '# SDBL flow');
            fprintf(prof, '%-60s\n', '# z thl qt u v tke');
            fprintf(prof, '%-20.15f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f\n', obj.pr');
            fclose(prof);
        end

        function generate_scalar(obj)
            preprocessing.addvar(obj, 'sc', zeros(length(obj.zf), obj.nsv+1));
            obj.sc(:,1) = obj.zf;
            if obj.nsv>0
                obj.sc(:,2) = obj.sv10;
            end
            if obj.nsv>1
                obj.sc(:,3) = obj.sv20;
            end
            if obj.nsv>2
                obj.sc(:,4) = obj.sv30;
            end
            if obj.nsv>3
                obj.sc(:,5) = obj.sv40;
            end
            if obj.nsv>4
                obj.sc(:,6) = obj.sv50;
            end
        end

        function write_scalar(obj)
            scalar = fopen(['scalar.inp.' obj.expnr], 'w');
            fprintf(scalar, '%-12s\n', '# SDBL flow');
            fprintf(scalar, '%-60s\n', '# z scaN,  N=1,2...nsv');
            fprintf(scalar, ['%-20.15f' repmat(' %-14.10f',[1,obj.nsv])  '\n'], obj.sc');
            fclose(scalar);
        end

        function set_nfcts(obj, nfcts)
            obj.nfcts = nfcts;
        end

        function write_vf(obj, vf)
            ncid = netcdf.create(['vf.nc.inp.' num2str(obj.expnr)], 'NC_WRITE');
            dimidrow = netcdf.defDim(ncid,'rows', obj.nfcts);
            dimidcol = netcdf.defDim(ncid,'columns', obj.nfcts);
            varid = netcdf.defVar(ncid,'view factor','NC_FLOAT',[dimidrow dimidcol]);
            netcdf.endDef(ncid);
            netcdf.putVar(ncid,varid,vf);
            netcdf.close(ncid);
        end

        function write_vfsparse(obj, vfsparse)
            [i,j,s] = find(vfsparse);
            fID = fopen([fpath 'vfsparse.inp.' num2str(obj.expnr)], 'w');
            fprintf(fID, '%d %d %.6f \n', [i, j, s]');
        end

        function write_svf(obj, svf)
            fname = ['svf.inp.' num2str(obj.expnr)];
            fileID = fopen(fname,'W');
            fprintf(fileID, '# sky view factors\n');
            fclose(fileID);
            dlmwrite(fname, svf, '-append','delimiter',' ','precision','%4f')
        end

        function write_facetarea(obj, facetarea)
            fname = ['facetarea.inp.' num2str(obj.expnr)];
            fileID = fopen(fname,'W');
            fprintf(fileID, '# area of facets\n');
            fclose(fileID);
            dlmwrite(fname, facetarea,'-append','delimiter',' ','precision','%4f')
        end

        function write_netsw(obj, Knet)
            fname = ['netsw.inp.' obj.expnr];
            fileID = fopen(fname, 'w');
            fprintf(fileID,'# %4s\n','net shortwave on facets [W/m2] (including reflections and diffusive)');
            fprintf(fileID,'%6.4f\n', Knet);
            fclose(fileID);
        end

        function write_Tfacinit(obj, Tfacinit)
            fname = ['Tfacinit.inp.' obj.expnr];
            fileID = fopen(fname, 'W');
            fprintf(fileID, '# Initial facet tempereatures in radiative equilibrium\n');
            fclose(fileID);
            dlmwrite(fname, Tfacinit, '-append','delimiter',' ','precision','%4f')
        end

        function write_Tfacinit_layers(obj, Tfacinit_layers)
           fname = ['Tfacinit_layers.inp.' obj.expnr];
           fileID = fopen(fname, 'W');
            fprintf(fileID, '# Initial facet tempereatures in radiative equilibrium\n');
            fclose(fileID);
            dlmwrite(fname, Tfacinit_layers, '-append','delimiter',' ','precision','%4f')
        end
    end

    methods (Static, Access = protected)

        function data = loadvar(filename, svar)

            % Load netcdf data.
            % loadvar(filename, svar)
            %
            % filename : Name of file (string).
            % svar     : Variable name (string).
            %
            % data     : Array of doubles.

            data = [];
            ncid = netcdf.open(filename,'NC_NOWRITE');
            [~, nvars] = netcdf.inq(ncid);
            found = false;
            for n=0:nvars-1
                if (strcmp(netcdf.inqVar(ncid,n),svar))
                    data = netcdf.getVar(ncid,n);
                    found = true;
                end
            end
            netcdf.close(ncid);

            if (~found)
                disp(['variable ', svar, ' not found!']);
            end
        end
    end
end
