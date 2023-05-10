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
                    %h.SetAccess = 'protected';
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

        function set_defaults(obj, ncpus)
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
            preprocessing.addvar(obj, 'imax', 64)  % # cells in x-direction
            preprocessing.addvar(obj, 'xsize', 64) % domain size in x-direction
            preprocessing.addvar(obj, 'jtot', 64)  % # cells in y-direction
            preprocessing.addvar(obj, 'ysize', 64) % domain size in y-direction
            preprocessing.addvar(obj, 'kmax', 96)  % # cells in z-direction

            if ceil(obj.jtot / ncpus) ~= floor (obj.jtot / ncpus)
                disp(['Possible jtot: ' num2str([2 3 4 5 6 7 8] * ncpus)])
                error('No. CPUs does not fit j grid size')
            end

            preprocessing.addvar(obj, 'dx', obj.xsize / obj.imax)
            preprocessing.addvar(obj, 'dy', obj.ysize / obj.jtot)

            %% &ENERGYBALANCE
            preprocessing.addvar(obj, 'lEB', 0)

            %% &WALLS
            preprocessing.addvar(obj, 'iwalltemp', 1)

            %% &PHYSICS
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

            %% &INPS
            preprocessing.addvar(obj, 'zsize', 96) % domain size in z-direction
            preprocessing.addvar(obj, 'lzstretch', 0) % switch for stretching z grid
            if obj.lzstretch
                preprocessing.addvar(obj, 'stretchconst', 0.01)
                preprocessing.addvar(obj, 'lstretchexp', 0)
                preprocessing.addvar(obj, 'lstretchtanh', 0)
                preprocessing.addvar(obj, 'lstretch2tanh', 0)
                preprocessing.addvar(obj, 'hlin', 0)
                preprocessing.addvar(obj, 'dzlin', 0)
                preprocessing.addvar(obj, 'dz', obj.dzlin)
            else
                preprocessing.addvar(obj, 'dz', obj.zsize / obj.kmax)
            end

            %if obj.lEB
                %preprocessing.addvar(obj, 'maxsize', 10); % maximum size of facets
            %else
                preprocessing.addvar(obj, 'maxsize', inf);
            %end
            
            preprocessing.addvar(obj, 'lmaxsizewalls', 0) % if maxsize affects buildings walls in addition to floors

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

            % Blocks
            preprocessing.addvar(obj, 'lblocksfile', 0) % switch for using blocks from a file
            if obj.lblocksfile
                preprocessing.addvar(obj, 'blocksfile', '') % name of blocks file
            end

            % Blocks
            preprocessing.addvar(obj, 'ltxtblocks', 0) % switch for using blocks from a file
            preprocessing.addvar(obj, 'ltxtmat', 0) % switch for using blocks from a file
            preprocessing.addvar(obj, 'loverhangs', 0) % switch for overhangs
            preprocessing.addvar(obj, 'loverlayers', 0) % switch for overlayers
            

            if obj.ltxtblocks
                preprocessing.addvar(obj, 'txtblocksfile', '') % name of blocks file
                if obj.ltxtmat
                    preprocessing.addvar(obj, 'txtmatfile', '') % name of materials file
                end                
                if obj.loverhangs
                    preprocessing.addvar(obj, 'overhangsfile', '') % name of overhangs file
                end
                if obj.loverlayers
                    preprocessing.addvar(obj, 'noverlayers', 1) % switch for overlayers
                    preprocessing.addvar(obj, 'overlayerfile', '') % name of overlayer file 
                end
            end
            preprocessing.addvar(obj, 'lpadoutlet', 0) % Remove blocks at outlet
            preprocessing.addvar(obj, 'ltxtTG', 0) % switch for using blocks from a file
            
            
            preprocessing.addvar(obj, 'lflat', 0) % switch for no blocks or facets (as in DALES)
            preprocessing.addvar(obj, 'lfloors', 0) % switch for just floors, no buildings.

            if (obj.lEB && (obj.lflat || obj.lfloors))
                error('Energy balance currently not implemented for flat domain')
                % Eventually lfloors should work with lEB.
            end

            preprocessing.addvar(obj, 'lcube', 0)   % switch for linear cubes
            preprocessing.addvar(obj, 'lstaggered', 0) % switch for staggered cubes
            preprocessing.addvar(obj, 'lstaggeredv', 0) % switch for staggered cubes with variable parameters
            preprocessing.addvar(obj, 'lcanyons', 0) % switch for infinite canyons
            if obj.nsv>0
                if (obj.lscasrcl && not(obj.lcanyons))
                    error('Scalar line sources only implemented for lcanyons')
                end
            end

            if (obj.lcube || obj.lstaggered || obj.lstaggeredv || obj.lcanyons)
                preprocessing.addvar(obj, 'blockheight', 16) % block height
                preprocessing.addvar(obj, 'blockwidth', 16)  % block width
                preprocessing.addvar(obj, 'canyonwidth', 16) % canyonwidth
            end

            preprocessing.addvar(obj, 'llidar', 0)
            if obj.llidar
                preprocessing.addvar(obj, 'sourcename', '')
                preprocessing.addvar(obj, 'dxinp', 1) % resolution of image [m/pixel]
                preprocessing.addvar(obj, 'dyinp', 1)
                preprocessing.addvar(obj, 'dzinp', 1)
                preprocessing.addvar(obj, 'centeri', 0) % center of area of interest in original image [pixel]
                preprocessing.addvar(obj, 'centerj', 0)
                preprocessing.addvar(obj, 'maxh', 0) % maximum height of buildings in image [m]
                preprocessing.addvar(obj, 'pad', 5) % padding. A padding of 0 makes only sense for idealised cases. There should be no building at domain edge
                preprocessing.addvar(obj, 'smallarea', round(150 / (obj.dx * obj.dy))) % objects smaller than this will be deleted
            end

            if obj.lEB
                preprocessing.addvar(obj, 'solarazimuth', 135); % solar azimuth angle
                preprocessing.addvar(obj, 'solarzenith', 28.4066); % solar zenith angle
                preprocessing.addvar(obj, 'centerweight', 12 / 32);
                preprocessing.addvar(obj, 'cornerweight', (1 - obj.centerweight) / 4);
                preprocessing.addvar(obj, 'I', 184.8775); % Direct solar irradiation [W/m2]
                preprocessing.addvar(obj, 'Dsk', 418.8041); % Diffuse incoming radiation [W/m2]
                preprocessing.addvar(obj, 'xazimuth', 90); % azimuth of x-direction wrt N. Default: x = East
                % north -> xazimuth = 0;
                % east  ->            90;
                % south ->            180;
                % west  ->            270;
                preprocessing.addvar(obj, 'solaz', obj.solarazimuth - obj.xazimuth + 90); % solar azimuth angle
                preprocessing.addvar(obj,  'lboundingwalls', 1)
            end
            preprocessing.addvar(obj, 'facT', 288.) % Initial facet temperatures.
            preprocessing.addvar(obj, 'nwalllayers', 3) % Number of facet layers
            
            
            preprocessing.addvar(obj, 'ltimedepsw', 0)
            preprocessing.addvar(obj, 'ltimedeplw', 0)
            preprocessing.addvar(obj, 'ltimedepnudge', 0)
            preprocessing.addvar(obj, 'starttime', 0)
            
            preprocessing.addvar(obj, 'ntimedepsw', 24)
            preprocessing.addvar(obj, 'ntimedeplw', 24)
            preprocessing.addvar(obj, 'nweatherdata', 24)
            
            preprocessing.addvar(obj, 'nblocks', 0)
            preprocessing.addvar(obj, 'nfcts', 0)

            %preprocessing.generate_walltypes(obj)
            preprocessing.addvar(obj, 'walltypes', [])

        end

        function generate_walltypes(obj)
            K = obj.nwalllayers;
            walltypes = [];

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
            walltypes = [walltypes; bw];
            size(bw)

%             % Floors (f)
%             id_f  = -1;
%             lGR_f = 0;
%             z0_f  = 0.05;
%             z0h_f = 0.00035;
%             al_f  = 0.5;
%             em_f  = 0.85;
%             D_f   = 0.5;
%             d_f   = D_f / K;
%             C_f   = 1.875e6;
%             l_f   = 0.75;
%             k_f   = 0.4e-6;
%             if (K == 3)
%                 % Reproduce the original walltypes.inp (d_f not constant for each layer)
%                 f =  [id_f, lGR_f, z0_f, z0h_f, al_f, em_f, 0.1, 0.2, 0.2, C_f * ones(1,K), l_f * ones(1,K), k_f * ones(1,K+1)];
%             else
%                 f =  [id_f, lGR_f, z0_f, z0h_f, al_f, em_f, d_f * ones(1,K), C_f * ones(1,K), l_f * ones(1,K), k_f * ones(1,K+1)];
%             end
%             walltypes = [walltypes; f];

%             id_f1  = -1;
%             lGR_f1 = 0;
%             z0_f1  = 0.05;
%             z0h_f1 = 0.00035;
%             al_f1  = 0.29;
%             em_f1  = 0.9;
%             %D_f1   = 0.15;
%             D_f1   = 1;
%             d_f1   = D_f1 / K;
%             C_f1   = 2120 * 920;
%             l_f1   = 0.7;
%             k_f1   = l_f1 / C_f1;
%             f1 =  [id_f1, lGR_f1, z0_f1, z0h_f1, al_f1, em_f1, d_f1 * ones(1,K), C_f1 * ones(1,K), l_f1 * ones(1,K), k_f1 * ones(1,K+1)];
%             walltypes = [walltypes; f1];
            
            id_f2  = -1;
            lGR_f2 = 0;
            z0_f2 = 0.05;
            z0h_f2 = 0.00035;
            al_f2  = 0.25;
            em_f2  = 0.9;
            D_f2   = 0.5;
            %d_f2   = D_f2 / (K-1);
            d_f2   = D_f2 / K;
            C_f2   = 1840 * 1500;
            l_f2   = 1.0;
            k_f2   = l_f2 / C_f2;
            %HTC_f2 = 5;
            f2 =  [id_f2, lGR_f2, z0_f2, z0h_f2, al_f2, em_f2, d_f2 * ones(1,K), C_f2 * ones(1,K), l_f2 * ones(1,K), k_f2 * ones(1,K+1)];
            walltypes = [walltypes; f2];
            %size(f2)

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
            walltypes = [walltypes; dm];
            size(dm)

%             % Concrete (c)
%             id_c  = 1;
%             lGR_c = 0;
%             z0_c  = 0.05;
%             z0h_c = 0.00035;
%             al_c = 0.5;
%             em_c = 0.85;
%             D_c = 0.36;
%             d_c = D_c / K;
%             C_c = 2.5e6;
%             l_c = 1;
%             k_c = 0.4e-6;
%             c = [id_c, lGR_c, z0_c, z0h_c, al_c, em_c, d_c * ones(1,K), C_c * ones(1,K), l_c * ones(1,K), k_c * ones(1,K+1)];
%             walltypes = [walltypes; c];
            
%             % Concrete walls (c1)
%             id_c1  = 1;
%             lGR_c1 = 0;
%             z0_c1  = 0.05;
%             z0h_c1 = 0.00035;
%             al_c1 = 0.33;
%             em_c1 = 0.9;
%             D_c1 = 0.15;
%             d_c1 = D_c1 / K;
%             C_c1 = 2200 * 880;
%             l_c1 = 1.28;
%             k_c1 = l_c1 / C_c1;
%             c1 = [id_c1, lGR_c1, z0_c1, z0h_c1, al_c1, em_c1, d_c1 * ones(1,K), C_c1 * ones(1,K), l_c1 * ones(1,K), k_c1 * ones(1,K+1)];
%             walltypes = [walltypes; c1];

            % Concrete walls (c1)
            id_c2  = 1;
            lGR_c2 = 0;
            z0_c2 = 0.05;
            z0h_c2 = 0.00035;
            al_c2 = 0.35;
            em_c2 = 0.9;
            
            D_c2 = 0.25;
            D_ins = 0.25;
            d_tot = (D_c2 + D_ins)  / K;
            
            K_c2 = D_c2 / d_tot;
            K_ins = D_ins / d_tot;
            
%             K_c2 = floor(D_c2 / d_tot)
%             K_ins = floor(D_ins / d_tot)
%             
%             remd_c2 = D_c2 - K_c2 * d_tot
%             remd_ins = D_ins - K_ins * d_tot
            
              
%             K_c2 = D_c2 / d_tot % assumes integer
%             K_ins = D_ins / d_tot % assumes integer
%             d_c2 = D_c2 / K_c2
%             d_ins = D_ins / K_ins
            
            C_c2 = 2200 * 880;
            l_c2 = 1.28;
            k_c2 = l_c2 / C_c2;
            
            C_ins = 12 * 1030;
            l_ins = 0.042;
            k_ins = l_ins / C_ins;
            
            %HTC_c2 = 4;
            c2 = [id_c2, lGR_c2, z0_c2, z0h_c2, al_c2, em_c2, [d_tot * ones(1,K_c2), d_tot * ones(1,K_ins)], [C_c2 * ones(1,K_c2), C_ins * ones(1,K_ins)], [l_c2 * ones(1,K_c2), l_ins * ones(1,K_ins)], [k_c2 * ones(1,K_c2), k_ins * ones(1,K_ins+1)]];
            walltypes = [walltypes; c2];
            size(c2)
            
            % Brick (b)
%             id_b  = 2;
%             lGR_b = 0;
%             z0_b  = 0.05;
%             z0h_b = 0.00035;
%             al_b = 0.5;
%             em_b = 0.85;
%             D_b = 0.36;
%             d_b = D_b / K;
%             C_b = 2.766667e6;
%             l_b = 0.83;
%             k_b = 0.3e-6;
%             b = [id_b, lGR_b, z0_b, z0h_b, al_b, em_b, d_b * ones(1,K), C_b * ones(1,K), l_b * ones(1,K), k_b * ones(1,K+1)];
%             walltypes = [walltypes; b];

            % Stone (s)
%             id_s  = 3;
%             lGR_s = 0;
%             z0_s  = 0.05;
%             z0h_s = 0.00035;
%             al_s = 0.5;
%             em_s = 0.85;
%             D_s = 0.36;
%             d_s = D_s / K;
%             C_s = 2.19e6;
%             l_s = 2.19;
%             k_s = 1e-6;
%             s = [id_s, lGR_s, z0_s, z0h_s, al_s, em_s, d_s * ones(1,K), C_s * ones(1,K), l_s * ones(1,K), k_s * ones(1,K+1)];
%             walltypes = [walltypes; s];

            % Wood (w)
%             id_w  = 4;
%             lGR_w = 0;
%             z0_w  = 0.05;
%             z0h_w = 0.00035;
%             al_w = 0.5;
%             em_w = 0.85;
%             D_w = 0.36;
%             d_w = D_w / K;
%             C_w = 1e6;
%             l_w = 0.1;
%             k_w = 0.1e-6;
%             w = [id_w, lGR_w, z0_w, z0h_w, al_w, em_w, d_w * ones(1,K), C_w * ones(1,K), l_w * ones(1,K), k_w * ones(1,K+1)];
%             walltypes = [walltypes; w];
            
            % Sand
%             id_sand = 6;
%             lGR_sand = 0;
%             z0_sand = 0.05;
%             z0h_sand = 0.00035;
%             al_sand = 0.30;
%             em_sand = 0.90;
%             D_sand = 1;
%             d_sand = D_sand / K;
%             C_sand = 1.28e6;
%             l_sand = 0.3;
%             k_sand = l_sand / C_sand;
%             sand = [id_sand, lGR_sand, z0_sand, z0h_sand, al_sand, em_sand, d_sand * ones(1,K), C_sand * ones(1,K), l_sand * ones(1,K), k_sand * ones(1,K+1)];
%             walltypes = [walltypes; sand];

            % GR1
%             id_GR1 = 11;
%             lGR_GR1 = 1;
%             z0_GR1 = 0.05;
%             z0h_GR1 = 0.00035;
%             al_GR1 = 0.25;
%             em_GR1 = 0.95;
%             D_GR1 = 0.6;
%             d_GR1 = D_GR1 / K;
%             C_GR1 = 5e6;
%             l_GR1 = 2;
%             k_GR1 = 0.4e-6;
%             GR1 = [id_GR1, lGR_GR1, z0_GR1, z0h_GR1, al_GR1, em_GR1, d_GR1 * ones(1,K), C_GR1 * ones(1,K), l_GR1 * ones(1,K), k_GR1 * ones(1,K+1)];
%             walltypes = [walltypes; GR1];

            % GR2
%             id_GR2 = 12;
%             lGR_GR2 = 1;
%             z0_GR2 = 0.05;
%             z0h_GR2 = 0.00035;
%             al_GR2 = 0.35;
%             em_GR2 = 0.90;
%             D_GR2 = 0.6;
%             d_GR2 = D_GR2 / K;
%             C_GR2 = 2e6;
%             l_GR2 = 0.8;
%             k_GR2 = 0.4e-6;
%             GR2 = [id_GR2, lGR_GR2, z0_GR2, z0h_GR2, al_GR2, em_GR2, d_GR2 * ones(1,K), C_GR2 * ones(1,K), l_GR2 * ones(1,K), k_GR2 * ones(1,K+1)];
%             walltypes = [walltypes; GR2];
            
            % Water
%             id_wat = 20;
%             lGR_wat = 1;
%             z0_wat = 0.05;
%             z0h_wat = 0.00035;
%             al_wat = 0.06;
%             em_wat = 0.98;
%             D_wat = 1;
%             d_wat = D_wat / K;
%             C_wat = 4.2e6;
%             l_wat = 0.6;
%             k_wat = l_wat / C_wat;
%             wat = [id_wat, lGR_wat, z0_wat, z0h_wat, al_wat, em_wat, d_wat * ones(1,K), C_wat * ones(1,K), l_wat * ones(1,K), k_wat * ones(1,K+1)];
%             walltypes = [walltypes; wat];
            
%             % Concrete walls (c1)
%             id_c1  = 7;
%             lGR_c1 = 0;
%             z0_c1  = 0.05;
%             z0h_c1 = 0.00035;
%             al_c1 = 0.33;
%             em_c1 = 0.9;
%             D_c1 = 0.15;
%             d_c1 = D_c1 / K;
%             C_c1 = 2200 * 880;
%             l_c1 = 1.28;
%             k_c1 = l_c1 / C_c1;
%             c1 = [id_c1, lGR_c1, z0_c1, z0h_c1, al_c1, em_c1, d_c1 * ones(1,K), C_c1 * ones(1,K), l_c1 * ones(1,K), k_c1 * ones(1,K+1)];
%             walltypes = [walltypes; c1];
            
%             %Stone floors (f1)
%             id_f1  = 8;
%             lGR_f1 = 0;
%             z0_f1  = 0.05;
%             z0h_f1 = 0.00035;
%             al_f1  = 0.29;
%             em_f1  = 0.9;
%             D_f1   = 0.15;
%             D_f1   = 1;
%             d_f1   = D_f1 / K;
%             C_f1   = 2120 * 920;
%             l_f1   = 0.7;
%             k_f1   = l_f1 / C_f1;
%             f1 =  [id_f1, lGR_f1, z0_f1, z0h_f1, al_f1, em_f1, d_f1 * ones(1,K), C_f1 * ones(1,K), l_f1 * ones(1,K), k_f1 * ones(1,K+1)];
%             walltypes = [walltypes; f1];

            obj.walltypes = walltypes;
            %preprocessing.addvar(obj, 'walltypes', walltypes);
        end

        function write_walltypes(obj)
            K = obj.nwalllayers;

            fname = ['walltypes.inp.', obj.expnr];

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
                fw_d = 6 + floor(log10(k));
                valstring2 = [valstring2, ['  %' num2str(fw_d) '.3f']];
            end
            for k=1:K
                fw_C = 14 + floor(log10(k));
                valstring2 = [valstring2, ['  %' num2str(fw_C) '.0f']];
            end
            for k=1:K
                fw_l = 13 + floor(log10(k));
                valstring2 = [valstring2, [' %' num2str(fw_l) '.4f']];
            end
            for k=1:K+1
                fw_k = 13 + floor(log10(k));
                valstring2 = [valstring2, [' %' num2str(fw_k) '.8f']];
            end

            valstring = [valstring1, valstring2];

            [nwalltypes, ~] = size(obj.walltypes);

            for i = 1:nwalltypes
                fprintf(fileID, sprintf(valstring, obj.walltypes(i,:)));
                fprintf(fileID, '\n');
            end

            fclose(fileID);
        end

        function plot_profiles(obj)
            figure
            subplot(141)
            plot(obj.pr(:, 2), 1:obj.kmax)
            title('Temperature')

            subplot(142)
            plot(obj.ls(:, 10), 1:obj.kmax)
            title('Radiative forcing')

            subplot(143)
            plot(obj.ls(:, 6), 1:obj.kmax)
            title('Subsidence')

            subplot(144)
            plot(obj.ls(:,2), 1:obj.kmax)
            hold on
            plot(obj.ls(:,3), 1:obj.kmax, 'r--')
            title('Velocity')
            legend('u', 'v')
        end

        function generate_xygrid(obj)
             preprocessing.addvar(obj, 'xf', 0.5 * obj.dx : obj.dx : obj.xsize - 0.5 * obj.dx);
             preprocessing.addvar(obj, 'yf', 0.5 * obj.dy : obj.dy : obj.ysize - 0.5 * obj.dy);
             preprocessing.addvar(obj, 'xh', 0 : obj.dx : obj.xsize);
             preprocessing.addvar(obj, 'yh', 0 : obj.dy : obj.ysize);
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
            ir  = obj.kmax - il;

            preprocessing.addvar(obj, 'zf', zeros(obj.kmax, 1));
            preprocessing.addvar(obj, 'dzf', zeros(obj.kmax, 1));
            preprocessing.addvar(obj, 'zh', zeros(obj.kmax+1, 1));

            obj.zf(1:il) = 0.5 * obj.dzlin : obj.dzlin : obj.hlin;
            obj.zh(1:il+1) = 0 : obj.dzlin : obj.hlin;

            gf = obj.stretchconst;

            while true
                obj.zh(il + 1:end) = obj.zh(il + 1) + (obj.zsize - obj.zh(il+1)) * (exp(gf * (0:1:ir) / (ir)) - 1)/(exp(gf) - 1); %dh has been replaced by zsize
                if (obj.zh(il+2) - obj.zh(il + 1)) < obj.dzlin
                    gf = gf - 0.01; %make sufficiently small steps to avoid an initial bump in dz
                else
                    if (obj.zh(end) - obj.zh(end - 1)) > 3 * obj.dzlin
                        disp('Warning: final grid spacing large - consider reducing domain height')
                    end
                    break
                end
            end

            for i = 1:obj.kmax
                obj.zf(i) = (obj.zh(i) + obj.zh(i+1)) / 2 ;
                obj.dzf(i) = obj.zh(i+1) - obj.zh(i);
            end
        end

        function stretch_tanh(obj)
            il = round(obj.hlin / obj.dzlin);
            ir  = obj.kmax - il;

            preprocessing.addvar(obj, 'zf', zeros(obj.kmax, 1));
            preprocessing.addvar(obj, 'dzf', zeros(obj.kmax, 1));
            preprocessing.addvar(obj, 'zh', zeros(obj.kmax + 1, 1));

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

            for i = 1:obj.kmax
                obj.zf(i) = 0.5 * (obj.zh(i) + obj.zh(i+1));
                obj.dzf(i) = obj.zh(i+1) - obj.zh(i);
            end
        end

        function stretch_2tanh(obj)
            il = round(obj.hlin / obj.dzlin);
            ir  = obj.kmax - il;

            preprocessing.addvar(obj, 'zf', zeros(obj.kmax, 1));
            preprocessing.addvar(obj, 'dzf', zeros(obj.kmax, 1));
            preprocessing.addvar(obj, 'zh', zeros(obj.kmax+1, 1));

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

            for i = 1:obj.kmax
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
            fprintf(lscale, '%-20.15f %-12.6f %-12.6f %-12.6f %-12.6f %-15.9f %-12.6f %-12.6f %-12.6f %-17.12f\n', obj.ls');
            fclose(lscale);
        end

        function generate_prof(obj)
            preprocessing.addvar(obj, 'pr', zeros(length(obj.zf), 6));
            obj.pr(:,1) = obj.zf;

            if obj.lapse
                thl = zeros(obj.kmax, 1);
                thl(1) = obj.thl0;
                for k = 1:obj.kmax - 1
                    thl(k+1) = thl(k) + obj.lapse * obj.zsize / obj.kmax;
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
        
        function write_timedepnudge(obj)
            fname = ['timedepnudge3.inp.', obj.expnr];
            fileID = fopen(fname, 'W');
            
            thl_nudge = circshift(obj.timedepTair, -obj.starttime);
            qt_nudge = ones(obj.ntimedepnudge,1) * 0.01; % Do this properly
            u_nudge = ones(obj.ntimedepnudge,1) * obj.u0;
            v_nudge = ones(obj.ntimedepnudge,1) * obj.v0;
            
            timedepnudge = zeros(obj.kmax, 5);
            for n=1:obj.ntimedepnudge
                n
                t = obj.timedeptime(n)+172800;
                fileID = fopen(fname, 'a');
                fprintf(fileID, 'height    thl        qt       u        v\n');
                fprintf(fileID, ['#    ', num2str(t), '\n']);
                timedepnudge(:,1) = obj.zf;
                timedepnudge(:,2) = thl_nudge(n); % Theta
                timedepnudge(:,3) = qt_nudge(n); % Moisture
                timedepnudge(:,4) = u_nudge(n); % u
                timedepnudge(:,5) = v_nudge(n); % v
                %timedepnudge(1,2:5) = 0.;
                fclose(fileID);
                dlmwrite(fname, timedepnudge, '-append','delimiter',' ','precision','%11.6f')
                fileID = fopen(fname, 'a');
                fprintf(fileID, '--------------------------------------------\n');
            end
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

        function generate_topo_from_bl(obj)
            preprocessing.addvar(obj, 'topomask', zeros(obj.jtot, obj.imax));
            preprocessing.addvar(obj, 'topo', zeros(obj.jtot, obj.imax));
            for n = 1:size(obj.bl, 1)
                obj.topo(obj.bl(n,3):obj.bl(n,4),obj.bl(n,1):obj.bl(n,2)) = obj.zh(obj.bl(n,6) + 1);
                obj.topomask(obj.bl(n,3):obj.bl(n,4),obj.bl(n,1):obj.bl(n,2)) = 1;
            end
        end

        function generate_topo_from_txt(obj)
            preprocessing.addvar(obj, 'topomask', zeros(obj.jtot, obj.imax));
            preprocessing.addvar(obj, 'topo', zeros(obj.jtot, obj.imax));
            topo_dum = load([obj.txtblocksfile],'-ascii');
            %obj.topo = round(topo_dum);
            %obj.topomask(topo_dum>0) = 1;
            
            %if obj.lEB
               %obj.topo(2:end-1,2:end-1) = obj.dz * floor(topo_dum(2:end-1,2:end-1) / obj.dz);
               %obj.topomask(obj.topo > 0) = 1;
            %else
                %obj.topo = round(topo_dum);
                obj.topo = obj.dz * floor(topo_dum / obj.dz); % Round to grid.
                obj.topomask(obj.topo > 0) = 1;
%                 obj.topo = obj.topo';
%                 obj.topomask = obj.topomask';
            %end
            
            if (obj.lpadoutlet)
                padding = 8;
               obj.topo(:,end-padding:end) = 0;
               obj.topomask(:,end-padding:end) = 0;
            end
        end
        
        function generate_topomat_from_txt(obj)
            preprocessing.addvar(obj, 'topomat', -ones(obj.jtot, obj.imax));             
            topomat = load([obj.txtmatfile],'-ascii');
            % -1 = floors
            % 0 = terrain
            % 1 = buildings
            %obj.topomat(obj.topo > 0) = topomat(obj.topo > 0); % Need to think about this when we have floors of differing material
            obj.topomat = topomat;
        end
        
        function generate_overhangs(obj)
            preprocessing.addvar(obj, 'overhangs', zeros(obj.jtot, obj.imax));
            overhangs = load([obj.overhangsfile],'-ascii'); 
            obj.overhangs = obj.dz * floor(overhangs / obj.dz); % round to the grid
            preprocessing.addvar(obj, 'overhangs_og', obj.overhangs);
            % Following line causes issues
            obj.overhangs(obj.overhangs > 0 & obj.overhangs == obj.topo) = obj.overhangs(obj.overhangs > 0 & obj.overhangs == obj.topo) - obj.dz;
            %obj.overhangs(obj.overhangs > 0 & obj.overhangs == obj.topo) = 0; obj.topo(obj.overhangs > 0 & obj.overhangs == obj.topo) = 0; % Remove very thin overhangs % Remove very thin overhangs
            %obj.overhangs(obj.topomask > 0) = obj.overhangs(obj.topomask > 0); % Make sure all overhangs are below a block
        end
        
         function generate_overlayer(obj)
%              if obj.noverlayers == 1
%                 preprocessing.addvar(obj, 'overlayer', zeros(obj.jtot, obj.imax));
%                 overlayer = load([obj.overlayerfile],'-ascii');
%                 obj.overlayer = overlayer;
%             
%                 preprocessing.addvar(obj, 'overlayermask', zeros(obj.jtot, obj.imax));
%                 obj.overlayermask(obj.overlayer > 0) = 1;
%                 %preprocessing.addvar(obj, 'layerheight', max(obj.topo, [], 'all')); % Check this is always true
%                 %preprocessing.addvar(obj, 'layermax', max(obj.topo, [], 'all')/obj.dz); % Check this is always true
%              else % Can probably be combined eventually
                 preprocessing.addvar(obj, 'overlayers', zeros(obj.jtot, obj.imax, obj.noverlayers));
                 load([obj.overlayerfile], 'overlayers');
                 obj.overlayers = overlayers;
                 
                 preprocessing.addvar(obj, 'overlayermask', zeros(obj.jtot, obj.imax, obj.noverlayers));
                 for n = 1:obj.noverlayers
                    overlayermask = zeros(obj.jtot, obj.imax);
                    size(obj.overlayers)
                    overlayer = squeeze(obj.overlayers(:,:,n));
                    size(overlayer)
                    overlayermask(overlayer > 0) = 1;
                    obj.overlayermask(:,:,n) = overlayer;
                 end
             %end
        end
        
        function generate_topo_from_txt_TG(obj)
            % Assumes terrain and geometry in two seperate files
%             % Combined
%             preprocessing.addvar(obj, 'topomask', zeros(obj.jtot, obj.imax));
%             preprocessing.addvar(obj, 'topo', zeros(obj.jtot, obj.imax));
%             topo_dum = load([obj.txtblocksfile],'-ascii');
%             obj.topo = topo_dum;
%             obj.topomask(topo_dum>0) = 1;
            
            % Terrain
            preprocessing.addvar(obj, 'topomask_terrain', zeros(obj.jtot, obj.imax));
            preprocessing.addvar(obj, 'topo_terrain', zeros(obj.jtot, obj.imax));
            topo_dum = load([obj.txtterrainfile],'-ascii');
            obj.topo_terrain = topo_dum;
            obj.topomask_terrain(topo_dum>0) = 1;
            
            % Geometry
            preprocessing.addvar(obj, 'topomask_geometry', zeros(obj.jtot, obj.imax));
            preprocessing.addvar(obj, 'topo_geometry', zeros(obj.jtot, obj.imax));
            topo_dum = load([obj.txtgeometryfile],'-ascii');
            obj.topo_geometry = topo_dum;
            obj.topomask_geometry(topo_dum>0) = 1;
            
            % Combined
            preprocessing.addvar(obj, 'topomask', obj.topo_terrain); % Assumes all geometry is on terrain
            preprocessing.addvar(obj, 'topo', obj.topo_terrain + obj.topo_geometry);            
        end
        
        function set_walltypes_from_topomat(obj)
            
           for n = 1:obj.nblockstotal
                il = obj.blocks(n,1);
                iu = obj.blocks(n,2) + 1;                
                jl = obj.blocks(n,3);
                ju = obj.blocks(n,4) + 1;
                mat = obj.topomat(jl, il);
                if (mat == -1) % No terrain or geometry: ocean
                    for m = 1:5
                        facid = obj.blocks(n,6+m);
                        if facid > 0
                            obj.facets(facid,2) = 20;
                        end
                    end
                elseif (mat == 0) % Terrain: sand
                    for m = 1:5
                        facid = obj.blocks(n,6+m);
                        if facid > 0
                            obj.facets(facid,2) = 6;
                        end
                    end
               elseif (mat == 1) % Streets
                    for m = 1:5
                        facid = obj.blocks(n,6+m);
                        if facid > 0
                            obj.facets(facid,2) = 8;
                        end
                    end
                
                elseif (mat == 2) % Bulding walls and roofs
                    for m = 1:5
                        facid = obj.blocks(n,6+m);
                        if facid > 0
                            obj.facets(facid,2) = 7;
                        end
                    end
                end
           end
        end
        
        function set_walltypes_TG(obj)
            % Set floors to water
            obj.floorfacets(:,2) = 20;
            obj.facets(end-obj.nfloorfacets+1:end,:) = obj.floorfacets;
                        
            for n = 1:obj.nblocks
                il = obj.blocks(n,1);
                iu = obj.blocks(n,2) + 1;                
                jl = obj.blocks(n,3);
                ju = obj.blocks(n,4) + 1;
                if (obj.topomask_terrain(il,jl) && ~obj.topomask_geometry(il,jl))
                    for m = 1:5
                        facid = obj.blocks(n,6+m);
                        if facid > 0
                           obj.facets(facid,2) = 3;
                        end
                    end
                end
            end
                        
        end

        function generate_bl_from_namoptions(obj)
            imax = obj.imax;
            jtot = obj.jtot;
            blockwidth = obj.blockwidth;
            blockheight = obj.blockheight;
            canyonwidth = obj.canyonwidth;
            nrows =  imax / (blockwidth + canyonwidth);

            if ceil(nrows) ~= floor(nrows)
                l = 0:0.5 * imax;
                ind = rem(0.5 * imax, l) == 0;
                err = ([l(ind); (0.5 * imax) ./ l(ind)]);
                disp('Block system does not fit grid')
                disp(['sum widths to: ' num2str(err(1,:))])
                disp(['Current width: ' num2str(blockwidth + canyonwidth)])
                error('Incorrect block system')
            end

            if obj.lstaggered
                nrows = imax / (blockwidth * 2);
                ncolumns = jtot / (blockwidth * 2);
                bl = zeros(nrows * ncolumns + nrows / 2, 13);
                bl(:, 5) = 0;
                bl(:, 6) = blockheight - 1;

                for n = 1:nrows
                    for nn = 0:ncolumns
                        if mod(n, 2) == 0
                            if nn == 0
                                bl(length(nonzeros(bl(:,3))) + 1, 3) = 1;
                                bl(length(nonzeros(bl(:,4))) + 1, 4) = 0.5 * blockwidth;
                            elseif nn == ncolumns
                                bl(length(nonzeros(bl(:,3))) + 1, 3) = jtot - 0.5 * blockwidth + 1;
                                bl(length(nonzeros(bl(:,4))) + 1, 4) = jtot;
                            else
                                bl(length(nonzeros(bl(:,3))) + 1, 3) = 1 + nn * blockwidth * 2 - 0.5 * blockwidth;
                                bl(length(nonzeros(bl(:,4))) + 1, 4) = bl(length(nonzeros(bl(:,4))) + 1, 3) + blockwidth - 1;
                            end
                            bl(length(nonzeros(bl(:,1))) + 1, 1) = - 0.5 * blockwidth + (2 * n - 1) * blockwidth + 1;
                            bl(length(nonzeros(bl(:,2))) + 1, 2) = bl(length(nonzeros(bl(:,2)))+1,1) + blockwidth - 1;
                        end
                    end
                    for nn = 0:ncolumns - 1
                        if mod(n,2) ~= 0
                            bl(length(nonzeros(bl(:,3))) + 1, 3) = 1 + blockwidth + nn * blockwidth * 2 - 0.5 * blockwidth;
                            bl(length(nonzeros(bl(:,4))) + 1, 4) = bl(length(nonzeros(bl(:,4))) + 1, 3) + blockwidth - 1;
                            bl(length(nonzeros(bl(:,1))) + 1, 1) = - 0.5 * blockwidth + (2 * n - 1) * blockwidth + 1;
                            bl(length(nonzeros(bl(:,2))) + 1, 2) = bl(length(nonzeros(bl(:,2))) +1 , 1) + blockwidth - 1;
                        end

                    end
                end
                
            elseif obj.lstaggeredv
                Nx = imax / (blockwidth + canyonwidth);
                Ny = jtot / (blockwidth + canyonwidth);
                bl = zeros(Nx * Ny + Nx / 2, 13);
                bl(:, 5) = 0;
                bl(:, 6) = blockheight - 1;
                
                even_xl = 0.5 * canyonwidth;
                even_xu = even_xl + blockwidth;
                odd_xl = even_xl + blockwidth + canyonwidth;
                odd_xu = odd_xl + blockwidth;
                
                for i = 0:(Nx-1)
                    %for j = 0:Ny-1
                    if mod(i, 2) == 0
                        i = i / 2;
                        for j = 0:(Ny-1)
                            bl(length(nonzeros(bl(:,1))) + 1, 1) = even_xl + i * 2 * (blockwidth + canyonwidth) + 1;
                            bl(length(nonzeros(bl(:,2))) + 1, 2) = even_xu + i * 2 * (blockwidth + canyonwidth);
                            bl(length(nonzeros(bl(:,3))) + 1, 3) = 0.5 * canyonwidth + j * (blockwidth + canyonwidth) + 1;
                            bl(length(nonzeros(bl(:,4))) + 1, 4) = 0.5 * canyonwidth + j * (blockwidth + canyonwidth) + blockwidth;
                        end
                    else
                        i = (i-1)/2;
                        bl(length(nonzeros(bl(:,1))) + 1, 1) = odd_xl + i * 2 * (blockwidth + canyonwidth) + 1;
                        bl(length(nonzeros(bl(:,2))) + 1, 2) = odd_xu + i * 2 * (blockwidth + canyonwidth);
                        bl(length(nonzeros(bl(:,3))) + 1, 3) = 1;
                        bl(length(nonzeros(bl(:,4))) + 1, 4) = 0.5 * blockwidth;
                        for j = 0:(Ny-2)
                            bl(length(nonzeros(bl(:,1))) + 1, 1) = odd_xl + i * 2 * (blockwidth + canyonwidth) + 1;
                            bl(length(nonzeros(bl(:,2))) + 1, 2) = odd_xu + i * 2 * (blockwidth + canyonwidth);
                            bl(length(nonzeros(bl(:,3))) + 1, 3) = 0.5 * blockwidth + canyonwidth + j * (blockwidth + canyonwidth) + 1;
                            bl(length(nonzeros(bl(:,4))) + 1, 4) = 0.5 * blockwidth + canyonwidth + j * (blockwidth + canyonwidth) + blockwidth;
                        end
                        bl(length(nonzeros(bl(:,1))) + 1, 1) = odd_xl + i * 2 * (blockwidth + canyonwidth) + 1;
                        bl(length(nonzeros(bl(:,2))) + 1, 2) = odd_xu + i * 2 * (blockwidth + canyonwidth);
                        bl(length(nonzeros(bl(:,3))) + 1, 3) = jtot - 0.5 * blockwidth + 1;
                        bl(length(nonzeros(bl(:,4))) + 1, 4) = jtot;
                    end
                end

            elseif obj.lcube
                xsize = obj.xsize; zsize = obj.zsize;
                kmax = obj.kmax;
                nrows = imax / (blockwidth * 2);
                ncolumns = jtot / (blockwidth * 2);
                bl = zeros(nrows * ncolumns, 13);
                for n = 1 : nrows
                    for nn = 0 : ncolumns - 1
                        bl((n - 1) * ncolumns + nn + 1, 1) = -0.5 * blockwidth + (2 * n - 1) * blockwidth + 1;
                        bl((n - 1) * ncolumns + nn + 1, 2) = bl((n - 1) * ncolumns + nn + 1, 1) + blockwidth - 1;
                        bl((n - 1) * ncolumns + nn + 1, 5) = 0;
                        bl((n - 1) * ncolumns + nn + 1, 6) = ceil(blockwidth * (xsize / imax) / (zsize / kmax));
                        bl((n - 1) * ncolumns + nn + 1, 3) = 1 + blockwidth / 2 + nn * blockwidth * 2;
                        bl((n - 1) * ncolumns + nn + 1, 4) = bl((n - 1) * ncolumns + nn + 1, 3) + blockwidth - 1;
                    end
                end

            elseif obj.lcanyons
                bl = zeros(nrows, 13);
                bl(1:nrows, 1) = (0.5 * canyonwidth + 1 : canyonwidth + blockwidth : imax - 0.5 * canyonwidth)';
                bl(1:nrows, 2) = bl(1:nrows, 1) + blockwidth - 1;
                bl(:, 3) = 1;
                bl(:, 4) = jtot;
                bl(:, 5) = 0;
                bl(1:nrows, 6) = blockheight - 1;
            end

            preprocessing.addvar(obj, 'bl', bl)
        end

        function generate_bl_from_file(obj)
            bl = dlmread(obj.blocksfile, '', 2, 0);
            preprocessing.addvar(obj, 'bl', bl)
        end

        function generate_topo_from_LIDAR(obj)
            maxh = obj.maxh;
            dx = obj.dx;
            dy = obj.dy;
            dz = obj.dz;
            ni = obj.imax;
            nj = obj.jtot;
            dxinp =obj.dxinp;
            dyinp = obj.dyinp;
            centeri = obj.centeri;
            centerj = obj.centerj;
            pad = obj.pad;
            smallarea = obj.smallarea;
            sourcename = obj.sourcename;
            A = imread(sourcename);  %read topo image
            [njorig, niorig, ~]=size(A);


            % since its a greyscale image all 3 rgb channels have the same value
            % substract value from 255 and scale by max value to get topography
            topot=(255-double(A(:,:,1)))/255*maxh;

            % if ltestplot
            % figure; imagesc(topot)
            % end

            % interpolate to a coarser grid if necessary
            if dx~=dxinp || dy~=dyinp %can't use imresize, have to do it manually since x and y scaling might be different
                xxxo=dxinp/2:dxinp:(niorig*dxinp); % original pixel centres
                yyyo=dyinp/2:dyinp:(njorig*dyinp);
                xxx=dx/2:dx:(niorig*dxinp); %desired pixel centres
                yyy=dy/2:dy:(njorig*dxinp);
                [X,Y]=meshgrid(xxx,yyy);
                topot=interp2(yyyo,xxxo,topot',Y,X,'nearest');
                clear xxxo yyyo X Y
            end

            % if ltestplot
            % figure; imagesc(topot);
            % end


            % upper and lower coordinates to select
            ip=round(centeri/dx)+ni/2-1-pad;
            im=round(centeri/dx)-ni/2+pad;
            jp=round(centerj/dy)+nj/2-1-pad;
            jm=round(centerj/dy)-nj/2+pad;

            topo=zeros(nj,ni);
            topo(pad+(1:(nj-2*pad)),pad+(1:(ni-2*pad)))=topot(jm:jp,im:ip);

            % if ltestplot
            % figure; imagesc(xf,yf,topo)
            % end

            % round height to grid
            topo=round(topo/dz)*dz;
            topoinit=topo;


            %% manual mainpulation of the topograpy come here, if necessary
            %test for varying building height
            % topo(170:180,173:209)=topo(170:180,173:209)+10;
            % topo(180:184,291:331)=topo(180:184,291:331)+12.5;
            % topo(143:148,173:214)=topo(143:148,173:214)+7.5;


            %% processing
            %% fill big holes
            topoh=imfill(topo,'holes');
            to=topoh-topo;
            toi=imbinarize(to);

            toi=bwareaopen(toi,smallarea);
            topo=topoh-toi.*to;
            clear to toi topoh topot

            % if ltestplot
            % figure; imagesc(xf,yf,topo); title('after filling holes')
            % end

            % create building mask
            topomask=imbinarize(topo);

            %% remove small objects
            topomask=bwareaopen(topomask,smallarea);

            topo=topo.*topomask;

            % if ltestplot
            %     figure
            %     imagesc(xf,yf,topo)
            %     title('after removing small objects')
            %     xlim([xh(1) xh(end)])
            %     ylim([yh(1) yh(end)])
            % end

            %% fill 1 cell gaps  %potentially problematic if the gap is between blocks of different height
            [topomask, topo] = fillgaps(topomask,topo,nj,ni,pad);




            % if ltestplot
            %     figure
            %     imagesc(xf,yf,topo)
            %     title('after filling gaps (size 1)')
            %     xlim([xh(1) xh(end)])
            %     ylim([yh(1) yh(end)])
            % end
            % if pad==0  %buildings to the edge, make sure there is no weird gaps at domain edge
            % [topomask, topo] = smoothborders(topomask,topo,nj,ni,pad);
            % end
            % if ltestplot
            %     figure
            %     imagesc(xf,yf,topo)
            %     title('after smoothing borders')
            %     xlim([xh(1) xh(end)])
            %     ylim([yh(1) yh(end)])
            % end

            %% remove cells with only 1 neighbour
            data=topo;
            datamask=topomask;

            c2=1;
            while c2~=0  %repeat until there is no cell with only one neighbour left
                c2=0;
                for i=2:ni-1
                    for j=2:nj-1
                        count=0;
                        if datamask(j,i)>0
                            c=datamask(j,i-1)+datamask(j,i+1)+datamask(j+1,i)+datamask(j-1,i);
                            if c>1
                                continue
                            else
                                c2=c2+1;
                                datamask(j,i)=0;
                                data(j,i)=0;
                            end
                        end
                    end
                end
            end

            topo=data;
            topomask=datamask;
            preprocessing.addvar(obj, 'topo', flipud(topo))
            preprocessing.addvar(obj, 'topomask', flipud(topomask))

            function [mask, image ] = fillgaps( mask,image,nj,ni,pad)
                %fills horizontal and vertical 1D gaps of width 1
                while true
                    change=false;
                    for j=pad+1:nj-pad-1
                        for i=pad+1:ni-pad-1
                            if j==1 || i==1 %do nothing at domain edge
                                continue

                            elseif sum(sum(mask(j-1:j,i-1:i+1)))==5 && mask(j,i)==0 && mask(j+1,i)==0 %towards south. If it has 5 block neighbours, is not a block and j+1 is also not a block then fill the gap towards the south. remember that images are upside down, so j+1 is going down
                                jj=j;
                                while true
                                    mask(jj,i)=1;
                                    image(jj,i)=image(jj-1,i);
                                    jj=jj+1;
                                    if jj==nj %reached end of domain
                                        break
                                    end
                                    if ~(sum(sum(mask(jj-1:jj,i-1:i+1)))==5 && mask(jj,i)==0 && mask(jj+1,i)==0)
                                        break
                                    end
                                end
                                change=true;
                            elseif sum(sum(mask(j:j+1,i-1:i+1)))==5 && mask(j,i)==0 && mask(j-1,i)==0 %towards north
                                jj=j;
                                while true
                                    mask(jj,i)=1;
                                    image(jj,i)=image(jj+1,i);
                                    jj=jj-1;
                                    if jj==1 %reached end of domain
                                        break
                                    end
                                    if ~(sum(sum(mask(jj:jj+1,i-1:i+1)))==5 && mask(jj,i)==0 && mask(jj-1,i)==0)
                                        break
                                    end
                                end
                                change=true;
                            elseif sum(sum(mask(j-1:j+1,i-1:i)))==5 && mask(j,i)==0 && mask(j,i+1)==0 %towards east
                                ii=i;
                                while true
                                    mask(j,ii)=1;
                                    image(j,ii)=image(j,ii-1);
                                    ii=ii+1;
                                    if ii==ni %reached end of domain
                                        break
                                    end
                                    if ~(sum(sum(mask(j-1:j+1,ii-1:ii)))==5 && mask(j,ii)==0 && mask(j,ii+1)==0)
                                        break
                                    end
                                end
                                change=true;
                            elseif sum(sum(mask(j-1:j+1,i:i+1)))==5 && mask(j,i)==0 && mask(j,i-1)==0 %towards west
                                ii=i;
                                while true
                                    mask(j,ii)=1;
                                    image(j,ii)=image(j,ii+1);
                                    ii=ii-1;
                                    if ii==1 %reached end of domain
                                        break
                                    end
                                    if ~(sum(sum(mask(j-1:j+1,ii:ii+1)))==5 && mask(j,ii)==0 && mask(j,ii-1)==0)
                                        break
                                    end
                                end
                            end

                        end


                    end
                    if ~change
                        break
                    end
                end
            end
        end

        function plot_bl(obj)
            figure
            %title('Blocks', 'interpreter','latex')
            view(3)
            if (obj.lstaggered || obj.lcube || obj.lcanyons)
                clr = [0.85,0.85,0.85];
                for i = 1:size(obj.bl, 1)
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1)  obj.xh(obj.bl(i,1))], [obj.yh(obj.bl(i,3))  obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.zh(obj.bl(i,6)+1)  obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1)], clr)
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,3))  obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3))], [obj.bl(i,5)  obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], clr)
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.bl(i,5)  obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], clr)
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) ], [obj.yh(obj.bl(i,4)+1)  obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3))], [obj.bl(i,5)  obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], clr)
                    patch([obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.bl(i,5)  obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], clr)
                    % patch([xh(1) xh(end) xh(end)  xh(1)], [yh(1)  yh(1) yh(end) yh(end)], [zh(1)  zh(1) zh(1) zh(1)], [245 245 245] ./ 255)
                end

            else
%                 for i = 1:size(obj.bl, 1)
%                     patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,1))], [obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1)], 'w')
%                     patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3))], [obj.bl(i,5) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], 'w')
%                     patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.bl(i,5) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], 'w')
%                     patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1))], [obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3))], [obj.bl(i,5) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], 'w')
%                     patch([obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.bl(i,5) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], 'w')
%                 end
            end

            zlim([0 obj.zh(end)]); %/(r.blockheight-1))
            xlim([0 obj.xh(end)]); %/(r.blockheight-1))
            ylim([0 obj.yh(end)]); %/(r.blockheight-1))

            set(gca,'ticklabelinterpreter','latex')
            xlabel('$x(\mathrm{m})$','interpreter','latex')
            ylabel('$y(\mathrm{m})$','interpreter','latex')
            zlabel('$z(\mathrm{m})$','interpreter','latex')
            set(gca,'BoxStyle','full','Box','on')
            daspect([1 1 1])
            grid on
            set(gca, 'FontSize', 12)

        end

        function plot_blocks(obj)
            figure
            %title('Blocks', 'interpreter', 'latex')
            view(3)

            if ~obj.lflat
                clr = [0.85, 0.85, 0.85];
                for i = 1:obj.nblockstotal
                    il = obj.blocks(i,1);
                    iu = obj.blocks(i,2);
                    jl = obj.blocks(i,3);
                    ju = obj.blocks(i,4);
                    kl = obj.blocks(i,5);
                    ku = obj.blocks(i,6);

                    if i <= obj.nblocks
                        patch([obj.xh(il)   obj.xh(iu+1) obj.xh(iu+1) obj.xh(il)]  , [obj.yh(jl)   obj.yh(jl)   obj.yh(ju+1) obj.yh(ju+1)], [obj.zh(ku+1) obj.zh(ku+1) obj.zh(ku+1) obj.zh(ku+1)], clr)
                        patch([obj.xh(il)   obj.xh(il)   obj.xh(iu+1) obj.xh(iu+1)], [obj.yh(jl)   obj.yh(jl)   obj.yh(jl)   obj.yh(jl)],   [obj.zh(kl)   obj.zh(ku+1) obj.zh(ku+1) obj.zh(kl)], clr)
                        patch([obj.xh(il)   obj.xh(il)   obj.xh(iu+1) obj.xh(iu+1)], [obj.yh(ju+1) obj.yh(ju+1) obj.yh(ju+1) obj.yh(ju+1)], [obj.zh(kl)   obj.zh(ku+1) obj.zh(ku+1) obj.zh(kl)], clr)
                        patch([obj.xh(il)   obj.xh(il)   obj.xh(il)   obj.xh(il)]  , [obj.yh(ju+1) obj.yh(ju+1) obj.yh(jl)   obj.yh(jl)],   [obj.zh(kl)   obj.zh(ku+1) obj.zh(ku+1) obj.zh(kl)], clr)
                        patch([obj.xh(iu+1) obj.xh(iu+1) obj.xh(iu+1) obj.xh(iu+1)], [obj.yh(jl)   obj.yh(jl)   obj.yh(ju+1) obj.yh(ju+1)], [obj.zh(kl)   obj.zh(ku+1) obj.zh(ku+1) obj.zh(kl)], clr)
                    else
                        patch([obj.xh(il)   obj.xh(iu+1) obj.xh(iu+1) obj.xh(il)]  , [obj.yh(jl)   obj.yh(jl)   obj.yh(ju+1) obj.yh(ju+1)], [obj.zh(ku+1) obj.zh(ku+1) obj.zh(ku+1) obj.zh(ku+1)], clr)
                        patch([obj.xh(il)   obj.xh(il)   obj.xh(iu+1) obj.xh(iu+1)], [obj.yh(jl)   obj.yh(jl)   obj.yh(jl)   obj.yh(jl)],   [0            obj.zh(ku+1) obj.zh(ku+1)            0], clr)
                        patch([obj.xh(il)   obj.xh(il)   obj.xh(iu+1) obj.xh(iu+1)], [obj.yh(ju+1) obj.yh(ju+1) obj.yh(ju+1) obj.yh(ju+1)], [0            obj.zh(ku+1) obj.zh(ku+1)            0], clr)
                        patch([obj.xh(il)   obj.xh(il)   obj.xh(il)   obj.xh(il)]  , [obj.yh(ju+1) obj.yh(ju+1) obj.yh(jl)   obj.yh(jl)],   [0            obj.zh(ku+1) obj.zh(ku+1)            0], clr)
                        patch([obj.xh(iu+1) obj.xh(iu+1) obj.xh(iu+1) obj.xh(iu+1)], [obj.yh(jl)   obj.yh(jl)   obj.yh(ju+1) obj.yh(ju+1)], [0            obj.zh(ku+1) obj.zh(ku+1)            0], clr)
                    end
                end
            end

            zlim([0 obj.zh(end)]); %/(r.blockheight-1))
            xlim([0 obj.xh(end)]); %/(r.blockheight-1))
            ylim([0 obj.yh(end)]); %/(r.blockheight-1))

            set(gca,'ticklabelinterpreter','latex')
            xlabel('$x(\mathrm{m})$','interpreter','latex')
            ylabel('$y(\mathrm{m})$','interpreter','latex')
            zlabel('$z(\mathrm{m})$','interpreter','latex')
            set(gca,'BoxStyle','full','Box','on')
            daspect([1 1 1])
            set(gca, 'FontSize', 12)
            grid on
        end
        
        function plot_blocks_new(obj)
            f = figure;
%             if obj.lEB % light comes from sun's position
%                 phi = (90 - obj.solaz) * pi / 180;
%                 theta = obj.solarzenith * pi / 180;
%                 x = cos(phi) * sin(theta);
%                 y = sin(phi) * sin(theta);
%                 z = cos(theta);
%                 light('Position', [x,y,z], 'Style', 'infinite');
%             else
%                 light('Position', [-0.75 -0.35 1], 'Style', 'infinite');
%             end
            
            %clr = [0.85, 0.85, 0.85];
            nblockstotal = size(obj.blocks,1);
            for i = 1:nblockstotal
            %for i = 1:obj.nblocks
%                 if i <= obj.nblocks
%                     clr = [0.85, 0.85, 0.85];
%                 else
%                     clr = [0, 0.4470, 0.7410];
%                 end
                il = obj.blocks(i,1);
                iu = obj.blocks(i,2);
                jl = obj.blocks(i,3);
                ju = obj.blocks(i,4);
                kl = obj.blocks(i,5) + 1;
                ku = obj.blocks(i,6) + 1;  % k-index starts at 0 in uDALES
                
                [~,~,facid] = find(obj.blocks(i,7:end), 1, 'first');
                %if ~isempty(facid)
%                     switch obj.facets(facid, 2) % wall type of facet
%                         case 20 % Water
%                             clr = [0, 0.4470, 0.7410];
%                         case 6 % Sand
%                             clr = [0.9290, 0.6940, 0.1250];
%                         case 1 % Concrete
%                             clr = [0.85, 0.85, 0.85];                            
%                         case 7 % Concrete
%                             clr = [0.85, 0.85, 0.85];                           
%                         case 8 % Floors
%                             clr = [0.5, 0.5, 0.5];
                            
                    %end
                %else % Internal block
                    %clr = [1, 1, 1];
                %end
                clr = [.5 .5 .5];
                V = [obj.xh(il) obj.yh(jl) obj.zh(kl); obj.xh(iu+1) obj.yh(jl)  obj.zh(kl); obj.xh(iu+1) obj.yh(ju+1)   obj.zh(kl); obj.xh(il) obj.yh(ju+1)   obj.zh(kl); ... % Lower 4 vertices
                    obj.xh(il) obj.yh(jl) obj.zh(ku+1); obj.xh(iu+1) obj.yh(jl) obj.zh(ku+1); obj.xh(iu+1) obj.yh(ju+1) obj.zh(ku+1); obj.xh(il) obj.yh(ju+1) obj.zh(ku+1)];    % Upper 4 vertices
                
                %    Bottom,  East,    North,   West,    South,   Top
                F = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
                patch('Faces',F,'Vertices',V,'EdgeColor','k','FaceColor',clr,'FaceLighting','flat','AmbientStrength', 0.3, 'SpecularStrength', 0.6, 'DiffuseStrength', 0.7);
            end
            
            axis equal
            grid on
            view(3)
            xlim([obj.xh(1), obj.xh(end)])
            ylim([obj.yh(1), obj.yh(end)])
            zlim([obj.zh(1), obj.zh(end)])
            xlabel('$x \mathrm{\ [m]}$', 'interpreter', 'latex')
            hXLabel = get(gca,'XLabel');
            set(hXLabel,'rotation', 25)
            ylabel('$y \mathrm{\ [m]}$', 'interpreter', 'latex')
            hYLabel = get(gca,'YLabel');
            set(hYLabel,'rotation', -35)
            zlabel('$z \mathrm{\ [m]}$', 'interpreter', 'latex')
            set(gca,'BoxStyle','full','Box','on')
            set(gca,'ticklabelinterpreter','latex')
            set(gca, 'FontSize', 12)

            
        end

        function plot_facets(obj, facids)
            figure
            %cmap = colormap('parula');

            if ~obj.lflat
                top = 1; west = 2; east = 3; north = 4; south = 5; bot = 6;
                
                %for i = 1:obj.nfcts
                for i = facids
                    
                    if (obj.facets(i,2) == 1) % concrete
                        clr = [0.75,0.75,0.75];
                    elseif (obj.facets(i,2) == -1) % roads
                        clr = [0.75,0.75,0.75];
                    elseif (obj.facets(i,2) == 11) % green
                        clr = [0.4660, 0.6740, 0.1880];
                    elseif (obj.facets(i,2) == 20) % blue
                        clr = [0, 0.4470, 0.7410];
                    elseif(obj.facets(i,2) == 3) % yellow
                        clr = [0.9290, 0.6940, 0.1250];    
                    elseif(obj.facets(i,2) == -101)
                        clr = [1, 1, 1];
                    end
                    
                    il = obj.facets(i, 6); iu = obj.facets(i, 7);
                    jl = obj.facets(i, 8); ju = obj.facets(i, 9);
                    kl = obj.facets(i, 10); ku = obj.facets(i, 11);
                    if i <= obj.nblockfcts
                        switch obj.facets(i, 1)
                            case {east, west}
                                x = [obj.xh(il), obj.xh(il), obj.xh(il), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(ju), obj.yh(ju)];
                                if kl == 0 && ku == 0
                                    z = [0, 0, 0, 0];
                                elseif kl == 0 && ku ~= 0
                                    z = [0, obj.zh(ku), obj.zh(ku), 0];
                                else
                                    z = [obj.zh(kl), obj.zh(ku), obj.zh(ku), obj.zh(kl)];
                                end
                            case {north, south}
                                x = [obj.xh(il), obj.xh(iu), obj.xh(iu), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(jl), obj.yh(jl)];
                                if kl == 0 && ku == 0
                                    z = [0, 0, 0, 0];
                                elseif kl == 0 && ku ~= 0
                                    z = [0, 0, obj.xh(ku), obj.xh(ku)];
                                else
                                    z = [obj.zh(kl), obj.zh(kl), obj.zh(ku), obj.zh(ku)];
                                end
                            case {top, bot}
                                x = [obj.xh(il), obj.xh(iu), obj.xh(iu), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(ju), obj.yh(ju)];

                                if kl == 0
                                    z = [0, 0, 0, 0];
                                else
                                    z = [obj.zh(kl), obj.zh(kl), obj.zh(kl), obj.zh(kl)];
                                end
                        end

                    elseif i <= obj.nblockfcts + obj.nboundingwallfacets
                        switch obj.facets(i, 1)
                            case east
                                ju = ju + 1;
                                kl = kl + 1;
                                ku = ku + 2;
                                x = [obj.xh(il), obj.xh(il), obj.xh(il), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(ju), obj.yh(ju)];
                                if kl == 0 && ku == 0
                                    z = [0, 0, 0, 0];
                                elseif kl == 0 && ku ~= 0
                                    z = [0, obj.zh(ku), obj.zh(ku), 0];
                                else
                                    z = [obj.zh(kl), obj.zh(ku), obj.zh(ku), obj.zh(kl)];
                                end
                            case west
                                il = il + 1;
                                iu = iu + 1;
                                ju = ju + 1;
                                kl = kl + 1;
                                ku = ku + 2;
                                x = [obj.xh(il), obj.xh(il), obj.xh(il), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(ju), obj.yh(ju)];
                                if kl == 0 && ku == 0
                                    z = [0, 0, 0, 0];
                                elseif kl == 0 && ku ~= 0
                                    z = [0, obj.zh(ku), obj.zh(ku), 0];
                                else
                                    z = [obj.zh(kl), obj.zh(ku), obj.zh(ku), obj.zh(kl)];
                                end


                            case north
                                iu = iu + 1;
                                kl = kl + 1;
                                ku = ku + 2;
                                x = [obj.xh(il), obj.xh(iu), obj.xh(iu), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(jl), obj.yh(jl)];
                                if kl == 0 && ku == 0
                                    z = [0, 0, 0, 0];
                                elseif kl == 0 && ku ~= 0
                                    z = [0, 0, obj.xh(ku), obj.xh(ku)];
                                else
                                    z = [obj.zh(kl), obj.zh(kl), obj.zh(ku), obj.zh(ku)];
                                end

                            case south
                                iu = iu + 1;
                                jl = jl + 1;
                                ju = ju + 1;
                                kl = kl + 1;
                                ku = ku + 2;
                                x = [obj.xh(il), obj.xh(iu), obj.xh(iu), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(jl), obj.yh(jl)];
                                if kl == 0 && ku == 0
                                    z = [0, 0, 0, 0];
                                elseif kl == 0 && ku ~= 0
                                    z = [0, 0, obj.xh(ku), obj.xh(ku)];
                                else
                                    z = [obj.zh(kl), obj.zh(kl), obj.zh(ku), obj.zh(ku)];
                                end


                            case {top, bot}
                                x = [obj.xh(il), obj.xh(iu), obj.xh(iu), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(ju), obj.yh(ju)];

                                if kl == 0
                                    z = [0, 0, 0, 0];
                                else
                                    z = [obj.zh(kl), obj.zh(kl), obj.zh(kl), obj.zh(kl)];
                                end
                        end
                    else
                        iu = iu + 1;
                        ju = ju + 1;
                        switch obj.facets(i, 1)
                            case {east, west}
                                x = [obj.xh(il), obj.xh(il), obj.xh(il), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(ju), obj.yh(ju)];
                                if kl == 0 && ku == 0
                                    z = [0, 0, 0, 0];
                                elseif kl == 0 && ku ~= 0
                                    z = [0, obj.zh(ku), obj.zh(ku), 0];
                                else
                                    z = [obj.zh(kl), obj.zh(ku), obj.zh(ku), obj.zh(kl)];
                                end
                            case {north, south}
                                x = [obj.xh(il), obj.xh(iu), obj.xh(iu), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(jl), obj.yh(jl)];
                                if kl == 0 && ku == 0
                                    z = [0, 0, 0, 0];
                                elseif kl == 0 && ku ~= 0
                                    z = [0, 0, obj.xh(ku), obj.xh(ku)];
                                else
                                    z = [obj.zh(kl), obj.zh(kl), obj.zh(ku), obj.zh(ku)];
                                end
                            case {top, bot}
                                x = [obj.xh(il), obj.xh(iu), obj.xh(iu), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(ju), obj.yh(ju)];

                                if kl == 0
                                    z = [0, 0, 0, 0];
                                else
                                    z = [obj.zh(kl), obj.zh(kl), obj.zh(kl), obj.zh(kl)];
                                end
                        end
                    end

                    %patch(x,y,z,clr,'FaceLighting','none','FaceAlpha',1);
                    patch('XData',x,'YData',y,'ZData',z,'EdgeColor','k','FaceColor',clr,'FaceLighting','flat','AmbientStrength', 0.3, 'SpecularStrength', 0.6, 'DiffuseStrength', 0.7);
                    %d = [0, 0, 0]; a = 0.25;
%                     switch(obj.facets(i, 1))
%                         case top
%                             d(3) = a * obj.dz(1);
%                         case west
%                             d(1) = -a * obj.dx(1);
%                         case east
%                             d(1) = a * obj.dx(1);
%                         case south
%                             d(2) = -a * obj.dy(1);
%                         case north
%                             d(2) = a * obj.dy(1);
%                     end

                    %text(mean(x) + d(1), mean(y) + d(2), mean(z) + d(3), num2str(i), 'horizontalalignment', 'center')
                    %hold on
                    %title('Facets')
                end

            end
            view(3)
            set(gca,'ticklabelinterpreter','latex')
            xlabel('$x(\mathrm{m})$','interpreter','latex')
            ylabel('$y(\mathrm{m})$','interpreter','latex')
            zlabel('$z(\mathrm{m})$','interpreter','latex')
            set(gca,'BoxStyle','full','Box','on')
            daspect([1 1 1])
            xlim([0 obj.xh(end)])
            ylim([0 obj.yh(end)])
            zlim([0 obj.zh(end)])
            set(gca, 'FontSize', 12)
            grid on
        end      

        function makeblocks(obj)
            topomask = obj.topomask;
            topo = obj.topo;
            over = obj.overhangs;
            imax = obj.imax;
            jtot = obj.jtot;
            dz = obj.dz;

            if (obj.lzstretch && (max(max(topo))>obj.hlin))
            	error('ERROR: hlin must be greater than the maximum vertical building index')
            end

            maxnrblocks = sum(topomask(:)); %allocate arrays with maximum size they can possibly have, reduce size later
            xmin = zeros(maxnrblocks,1); %store lower x bound of blocks
            xmax = zeros(maxnrblocks,1); %store upper x bound of blocks
            ymin = zeros(maxnrblocks,1); %store lower y bound of blocks
            ymax = zeros(maxnrblocks,1); %store upper y bound of blocks
            zmin = zeros(maxnrblocks,1); %store lower z bound of blocks
            zmax = zeros(maxnrblocks,1); %store upper z bound of blocks

            blchecked=zeros(size(topomask)); %mask of blocks which have already been checked
            indexmask=zeros(size(topomask)); %mask of the indeces of the (new) blocks
                
            topo_over = topo - over;
            count = 0;
            for i = 1:imax %loop over all x
                xminl = i;
                xmaxl = i;
                j = 1;
                while j <= jtot %loop along j
                    if topomask(j, i) == 0 %not a building
                        j = j + 1; %check next j
                        continue
                    else
                        yminl = j;
                        %heightgradient = diff(topo(j:end, i));
                        heightgradient = diff(topo(j:end, i) - over(j:end, i)); % If overhangs != 0 then we need to divide the block
                        if isempty(heightgradient) % if at the end of the y-direction (je)
                            heightchangey = j;
                        elseif all(heightgradient == 0) % same block/ floor until the end of the domain
                            if ~obj.lmaxsizewalls
                                heightchangey = jtot;
                            else
                                heightchangey = min(obj.maxsize + j, jtot);
                            end
                        else
                            if ~obj.lmaxsizewalls
                                heightchangey = find(heightgradient ~=0 , 1) + j - 1;  %last cell with same height as j (i.e. there is a height change betwen this and the next cell
                            else
                                heightchangey = min(obj.maxsize, find(heightgradient ~=0 , 1)) + j - 1;
                            end
                        end
                        
                         % If at any point from ymin to ymax there is a difference in height with the adjacent cells in x direction,
                         % or there is an overhang adjacent 
                        for jj = yminl+1:heightchangey
                            if (i == 1)
                                if ((any(topo(jj-1:jj,i+1) > 0)) && (topo(jj, i+1) ~= topo(jj-1,i+1))) || ((any(over(jj-1:jj,i+1) > 0)) && (over(jj, i+1) ~= over(jj-1,i+1)))
                                    heightchangey = jj - 1;
                                    break
                                end
                            elseif (i == imax)
                                if (any(topo(jj-1:jj,i-1) > 0) && (topo(jj,i-1)~= topo(jj-1,i-1))) || (any(over(jj-1:jj,i-1) > 0) && (over(jj,i-1)~= over(jj-1,i-1)))
                                    heightchangey = jj-1;
                                    break
                                end
                            else
                                if ((any(topo(jj-1:jj,i-1) > 0) && (topo(jj,i-1) ~= topo(jj-1,i-1))) || ((any(topo(jj-1:jj,i+1)>0)) && (topo(jj,i+1)~= topo(jj-1,i+1)))) ...
                                || ((any(over(jj-1:jj,i-1) > 0) && (over(jj,i-1) ~= over(jj-1,i-1))) || ((any(over(jj-1:jj,i+1)>0)) && (over(jj,i+1)~= over(jj-1,i+1))))
                                    heightchangey = jj-1; % overwrite heightchangey as we need to truncate block earlier!
                                    break
                                end
                            end
                        end

                        ymaxl = heightchangey;  % end of the current block (either floor or block with different height comes next)
                        ztemp = zeros(6,1);
                        ztempoverhangs = zeros(6,1);

                        % tg3315 changed from commented below as can have yminl==1 and
                        % ymaxl==nj
                        if yminl == 1
                            ztemp(2,1) = NaN;
                            ztempoverhangs(2,1) = NaN;
                        else
                            ztemp(2,1) = topo(yminl-1,xminl); % Topo height in -y direction 
                            ztempoverhangs(2,1) = over(yminl-1,xminl); % Overhang height in -y direction
                        end

                        if xminl == 1
                            ztemp(3,1) = NaN;
                            ztempoverhangs(3,1) = NaN;
                        else
                            ztemp(3,1) = topo(yminl, xminl-1); % Topo height in -x direction
                            ztempoverhangs(3,1) = over(yminl, xminl-1); % Overhang height in -x direction
                        end

                        if xmaxl == imax
                            ztemp(4,1) = NaN;
                            ztempoverhangs(4,1) = NaN;
                        else
                            ztemp(4,1) = topo(yminl, xminl+1); % Topo height in +x direction 
                            ztempoverhangs(4,1) = over(yminl, xminl+1); % Overhang height in +x direction
                        end

                        if ymaxl == jtot
                            ztemp(5,1) = NaN;
                            ztempoverhangs(5,1) = NaN;
                        else
                            %ztemp(5,1) = topo(yminl+1, xminl); % Topo height in +y direction
                            ztemp(5,1) = topo(ymaxl+1, xminl); % Height of next block. Not sure about this but if not done get facets that are too large
                            %ztempoverhangs(5,1) = over(yminl+1, xminl); % Overhang height in +y direction 
                            ztempoverhangs(5,1) = over(ymaxl+1, xminl); % Height of next overhang. Not sure about this but if not done see some missing blocks
                        end
                        
%                         if ymaxl ~= jtot
%                             zcutmax = topo(ymaxl+1, xminl); % Height of next block
%                         else
%                             zcutmax = 0;
%                         end

                         %zcuttemp = sort(ztemp(ztemp < topo(yminl, xminl) & ztemp > 0));
                         %zcuttemp = sort([ztemp(ztemp < topo(yminl, xminl) & ztemp > 0); zcutmax]);
                        zcuttemp = sort([ztemp(ztemp < topo(yminl, xminl) & ztemp > 0); ...
                                         ztempoverhangs(ztempoverhangs < topo(yminl, xminl) & ztempoverhangs > 0); ...
                                         over(yminl, xminl)]);
                        
                        zcuttemp = unique(zcuttemp);                        
                        zcut = [0; zcuttemp; topo(yminl,xminl)];
                        zcut = unique(zcut);
                        
                        % Also max vertical facet size < maxsize
%                         'before'
%                         zcut
                        if obj.lmaxsizewalls
                            p=1;
                            while p < length(zcut)
                                gap = zcut(p+1) - zcut(p);
                                if gap > obj.maxsize %&& mod(gap - obj.maxsize, dz)==0
                                    zcut = [zcut(1:p); zcut(p) + obj.maxsize; zcut(p+1:end)];
                                end
                                p=p+1;
                            end
                        end
%                         'after'
%                         zcut
                        
                        
                        for kc = 1:size(zcut, 1)-1
                            count = count+1;
                            blchecked(yminl:ymaxl,xminl:xmaxl) = blchecked(yminl:ymaxl,xminl:xmaxl)+1;
                            indexmask(yminl:ymaxl,xminl:xmaxl) = count;
                            xmin(count) = xminl;
                            xmax(count) = xmaxl;
                            ymin(count) = yminl;
                            ymax(count) = ymaxl;
                            %if overhangs(yminl,xminl) > 0
                                %zmin(count) = 0;
                                %zmax(count) = 0;
                            %else
                            if zcut(kc,1) == 0 % kc==1
                                zmin(count) = 0;
                            else
                                %disp('zcut(kc,1) ~= 0')
                                zmin(count) = zcut(kc,1) / dz;
                            end
                            zmax(count) = zcut(kc+1,1) / dz - 1;
                            %zmin(count)
                            if zmin(count) > zmax(count)
                                disp('zmin > zmax')
                            end
                            %end
                        end
                        j = heightchangey + 1;  % move to index after current block
                    end
                end
            end

            %shorten arrays
            ymin((count+1):end)=[];
            ymax((count+1):end)=[];
            xmin((count+1):end)=[];
            xmax((count+1):end)=[];
            zmin((count+1):end)=[];
            zmax((count+1):end)=[];

            %% add y slices of same dimension, aggregate along x
            %keep old xmax etc.. to keep better track whats happening
            xmax2 = xmax; ymax2 = ymax; zmax2 = zmax;
            xmin2 = xmin; ymin2 = ymin; zmin2 = zmin;
            dsize = 1;
            %dsize = 0;
            sizeold = count;
            while dsize > 0  %do as long as there are unmerged blocks % tg3315 need to also check here that you do not merge if they have a block above...!
                i = 1;
                while 1  %do along x
                    a = xmax2(i);  %upper x index of this block
                    bv = find(xmin2 == (a+1));  %all blocks with a lower x bound 1 bigger than this blocks upper x bound
                    b2 = bv(ymin2(bv) == ymin2(i)); %all of those blocks with also the same lower y bound
                    b3 = b2(zmin2(b2) == zmin2(i)); %all of those with same lower z bound - should return just 1 block (the adjacent block in the x direction)
                    if ~isempty(b3)
                        %[xmin2(i), xmax2(i), ymin2(i), ymax2(i), zmin2(i), zmax2(i)]
                        %[xmin2(b3), xmax2(b3), ymin2(b3), ymax2(b3), zmin2(b3), zmax2(b3)] 
                        %if all(ymax2(b3) == ymax2(i)) && all(zmax2(b3) == zmax2(i)) && topo_over(ymin2(b3),xmin2(b3)) == topo_over(ymin2(i), xmin2(i)) %&& all(zmin2(b2)==zmin2(i)) %if they also have the same upper y bound and the same height
                        if all(ymax2(b3) == ymax2(i)) && all(zmax2(b3) == zmax2(i)) && topo(ymin2(b3),xmin2(b3)) == topo(ymin2(i), xmin2(i)) ...%&& all(zmin2(b2)==zmin2(i)) ... %if they also have the same upper y bound and the same height
                                                                                    && over(ymin2(b3),xmin2(b3)) == over(ymin2(i), xmin2(i))
                            % additional check to make sure we do not merge blocks in a way that causes internal-external facets at this point
                            %if (topo_over(max(1,ymin2(b3)-1),xmax2(b3)) == topo_over(max(1,ymin2(i)-1),xmax2(i))) && (topo_over(min(ymax2(b3)+1, jtot),xmax2(b3)) == topo_over(min(ymax2(i)+1, jtot), xmax2(i)))
                            if (topo(max(1,ymin2(b3)-1),xmax2(b3)) == topo(max(1,ymin2(i)-1),xmax2(i))) && (topo(min(ymax2(b3)+1, jtot),xmax2(b3)) == topo(min(ymax2(i)+1, jtot), xmax2(i))) ...
                            && (over(max(1,ymin2(b3)-1),xmax2(b3)) == over(max(1,ymin2(i)-1),xmax2(i))) && (over(min(ymax2(b3)+1, jtot),xmax2(b3)) == over(min(ymax2(i)+1, jtot), xmax2(i)))    
                                
                                if ~obj.lmaxsizewalls
                                    xmax2(i) = xmax2(b3); %merge
                                    xmax2(b3) = []; ymax2(b3) = []; zmax2(b3)  =[]; xmin2(b3) = []; ymin2(b3) = []; zmin2(b3) = []; %remove the just merged block from list of blocks
                                else
                                    if xmax2(b3) - xmin2(i) < obj.maxsize 
                                        xmax2(i) = xmax2(b3); %merge
                                        xmax2(b3) = []; ymax2(b3) = []; zmax2(b3)  =[]; xmin2(b3) = []; ymin2(b3) = []; zmin2(b3) = []; %remove the just merged block from list of blocks
                                    end
                                end
                                
                            end
                        end
                        %end
                    end
                    i = i+1; %check furhter along x
                    if i >= length(xmin2)
                        break
                    end
                end
                dsize = sizeold-length(xmax2);
                sizeold = length(xmax2);
            end
            count2 = sizeold;            

            %make fields again
            datamean1 = zeros(size(topomask));
            datamean2 = zeros(size(topomask));
            %datamean3 = zeros(size(topomask));
            indexmask2 = zeros(size(topomask));

            %make new matrices
            for i=1:count
                datamean1(ymin(i):ymax(i),xmin(i):xmax(i)) = zmax(i);
            end
            for i=1:count2
                indexmask2(ymin2(i):ymax2(i),xmin2(i):xmax2(i)) = i; % What happens with overhangs...
                datamean2(ymin2(i):ymax2(i),xmin2(i):xmax2(i)) = zmax2(i);
            end

            %% split slices up according to rules in x and y, do z later
            xmin3 = zeros(maxnrblocks,1);
            xmax3 = zeros(maxnrblocks,1);
            ymin3 = zeros(maxnrblocks,1);
            ymax3 = zeros(maxnrblocks,1);
            zmin3 = zeros(maxnrblocks,1);
            zmax3 = zeros(maxnrblocks,1);

            xmin3(1:count2) = xmin2;
            xmax3(1:count2) = xmax2;
            ymin3(1:count2) = ymin2;
            ymax3(1:count2) = ymax2;
            zmin3(1:count2) = zmin2;
            zmax3(1:count2) = zmax2;


            % SO: I think the following does need to be done any more: 
            % The thing we want to avoid (for EB) is having overlapping facets
            % of different sizes that are not entirely internal.
            % I have hopefully implemented this earlier on in makeblocks,
            % so it will be true even if lEB = false, which results in more facets
            % but the geometry is now unaffected by the lEB switch, which has
            % some advantages.
            % The reason for removing this (other than it possibly being unnecessary)
            % is that it causes problems sometimes, i.e. gaps in the geometry.                     
            %if obj.lEB
            if false
                indexmask3 = indexmask2;
                change = true;
                count3 = count2;

                maxx = 0;
                cmax = 0;

                while change  %do until there is no changes anymore
                    change = false;
                    for c = 1:count3 %check if all blocks have same dimension as neighbour
                        %only have to check x-1 and x+1, since there can't be a
                        %neightbour on y-1 or y+1 (due to the inital slicing along y)
                        %index of blocks to the left

                        if ~(xmin3(c) == 1) %on left domain boundary, don't need to check for neighbouring blocks
                            illb = indexmask3(ymin3(c):ymax3(c),xmin3(c)-1);
                            iillb = find(illb>0,1);
                            if ~isempty(illb)
                                ilb = illb(iillb);
                            else
                                ilb = 0;
                            end
                        else
                            ilb = 0;
                        end

                        if ~(xmax3(c) == imax) %on right domain boundary, don't check for neighbouring blocsk
                            irrb = indexmask3(ymin3(c):ymax3(c),xmax3(c)+1);
                            iirrb = find(irrb>0,1);
                            if ~isempty(irrb)
                                irb = irrb(iirrb);
                            else
                                irb = 0;
                            end
                        else
                            irb = 0;
                        end
                        %   irt=indexmask(ymin(c),xmax(c)+1);
                        if (irb == 0) %no block to the right
                            yminr = ymin3(c); %use y coordinates of this block
                            ymaxr = ymax3(c);
                        else
                            yminr = ymin3(irb);  %use y coordinates of block to the right
                            ymaxr = ymax3(irb);
                        end

                        if (ilb == 0) %no block to the left
                            ymaxl = ymax3(c); %use y coordinates of this block
                            yminl = ymin3(c);
                        else
                            ymaxl = ymax3(ilb); %use y coordinates of block to the left
                            yminl = ymin3(ilb);
                        end

                        if xmin(3) > maxx
                            maxx = xmin3(c);
                            %maxx/ni;
                        end

                        if c/count3 > cmax
                            cmax=c/count3;
                        end

                        %
                        %

                        if (ymaxl < ymax3(c))% & zmax3(ilb)==zmax3(c)  %this block extends further (up) than left neighbour, this block has to be split
                            %x  x  x  x      %x  x  x  x
                            %x    [c] x      %x    [n] x
                            %x [l][c] x      %x [l][c] x
                            %x [l][c] x ===> %x [l][c] x
                            %x [l][c] x      %x [l][c] x
                            %x  x  x  x      %x  x  x  x
                            change = true;            %keep checking
                            count3 = count3+1;        %add new block
                            ymin3(count3) = ymaxl+1;  %new block starts 1 above xmax of left neighbour
                            ymax3(count3) = ymax3(c); %new block ends at xmax of this block
                            xmin3(count3) = xmin3(c); %new block has same x coordinates as this block
                            xmax3(count3) = xmax3(c);
                            zmin3(count3) = zmin3(c); %keep the same heights
                            zmax3(count3) = zmax3(c);
                            indexmask3(ymin3(count3):ymax3(count3),xmin3(count3):xmax3(count3)) = count3;
                            ymax3(c) = ymaxl;         %shorten this block
                            break                     %continue forloop from the start

                        elseif (yminl > ymin3(c))% & zmin3(ilb)==zmin3(c)  %this block extends further (down) than left neighbour, this block has to be split
                            %x  x  x  x      %x  x  x  x
                            %x [l][c] x      %x [l][c] x
                            %x [l][c] x      %x [l][c] x
                            %x [l][c] x ===> %x [l][c] x
                            %x    [c] x      %x    [n] x
                            %x  x  x  x      %x  x  x  x
                            change = true; %keep checking
                            count3 = count3+1; %add new block
                            ymin3(count3) = ymin3(c); %new block starts at same ymin as this
                            ymax3(count3) = yminl-1;  %new block ends 1 below ymin of left neighbour
                            xmin3(count3) = xmin3(c); %new block has same x coordinates as this block
                            xmax3(count3) = xmax3(c);
                            zmin3(count3) = zmin3(c); %keep the same heights
                            zmax3(count3) = zmax3(c);
                            indexmask3(ymin3(count3):ymax3(count3),xmin3(count3):xmax3(count3)) = count3;  %overwrite the index in the indexmask with the new one
                            ymin3(c)=yminl;         %shorten this block
                            break                   %continue forloop from the start

                        elseif (yminr > ymin3(c))% & zmin3(irb)==zmin3(c) %this block extends further (down) than reight neighbour, this block has to be split
                            %x  x  x  x      %x  x  x  x
                            %x [c][r] x      %x [c][r] x
                            %x [c][r] x      %x [c][r] x
                            %x [c][r] x ===> %x [c][r] x
                            %x [c]    x      %x [n]    x
                            %x  x  x  x      %x  x  x  x
                            change = true;            %keep checking
                            count3 = count3+1;        %add new block
                            ymin3(count3) = ymin3(c); %new block starts at same ymin as this
                            ymax3(count3) = yminr-1;  %new block ends 1 below ymin of reight neighbour
                            xmin3(count3) = xmin3(c); %new block has same x coordinates as this block
                            xmax3(count3) = xmax3(c);
                            zmin3(count3) = zmin3(c); %keep the same heights
                            zmax3(count3) = zmax3(c);
                            indexmask3(ymin3(count3):ymax3(count3),xmin3(count3):xmax3(count3)) = count3; %overwrite the index in the indexmask with the new one
                            ymin3(c) = yminr;         %shorten this block
                            break                     %continue for-loop from the start

                        elseif (ymaxr < ymax3(c))% & zmax3(irb)==zmax3(c) %this block extends further (up) than reight neighbour, this block has to be split
                            %x  x  x  x      %x  x  x  x
                            %x [c]    x      %x [n]    x
                            %x [c][r] x      %x [c][r] x
                            %x [c][r] x ===> %x [c][r] x
                            %x [c][r] x      %x [c][r] x
                            %x  x  x  x      %x  x  x  x
                            change = true;             %keep checking
                            count3 = count3+1;         %add new block
                            ymin3(count3) = ymaxr+1 ;  %new block starts 1 above xmax of reight neighbour
                            ymax3(count3) = ymax3(c);  %new block ends at xmax of this block
                            xmin3(count3) = xmin3(c);  %new block has same x coordinates as this block
                            xmax3(count3) = xmax3(c);
                            zmin3(count3) = zmin3(c); %keep the same heights
                            zmax3(count3) = zmax3(c);
                            indexmask3(ymin3(count3):ymax3(count3),xmin3(count3):xmax3(count3)) = count3;
                            ymax3(c) = ymaxr; %shorten this block
                            break %continue forloop from the start

                        end

                    end


                end
                ymin3((count3+1):end) = [];
                ymax3((count3+1):end) = [];
                xmin3((count3+1):end) = [];
                xmax3((count3+1):end) = [];
                zmin3((count3+1):end) = [];
                zmax3((count3+1):end) = [];


            else
                count3 = count2;
                xmin3 = xmin2;
                xmax3 = xmax2;
                ymin3 = ymin2;
                ymax3 = ymax2;
                zmin3 = zmin2;
                zmax3 = zmax2;

            end

            %% aggregate building internal blocks  (only rectangles)
            %  1 2 3 4     1 2 3 4
            %  5 6 7 8     5i 6 6 7
            %  9 a b c  -> 8 6 6 9
            %  d e f g     a b c d

            %use that all neighbours have the same size in 1 dimension % tg3315 mot
            %true now...

            datamean4 = zeros(size(topomask));
            indexmask4 = zeros(size(topomask));
            internalmask = zeros(size(topomask));

            for i = 1:count3
                indexmask4(ymin3(i):ymax3(i), xmin3(i):xmax3(i)) = i;
                datamean4(ymin3(i):ymax3(i), xmin3(i):xmax3(i)) = zmax3(i);
            end

            for k = 1:count3
                i = indexmask4(ymin3(k), xmin3(k)); %index of this block (k = i ?might be unnecessary)
                %check if it is internal
                xind = xmin3(i):xmax3(i);
                yind = ymin3(i):ymax3(i);
                temp1=zeros(length(yind),2);
                temp2=zeros(2,length(xind));
                if xmin3(i) == 1 %at edge
                    temp1(:,1) = 999999;    %dummy value
                else
                    temp1(:,1) = topo(yind,xmin3(i)) - topo(yind,xmin3(i)-1); % tg3315 changed so do not get internal blocks with adjacent different sizes %indexmask4(yind, xmin3(i)-1);
                end
                if xmax3(i) == imax %at edge
                    temp1(:,2) = 999999;
                else
                    temp1(:,2) = topo(yind,xmax3(i)) - topo(yind,xmax3(i)+1); %indexmask4(yind, xmax3(i)+1);
                end
                if ymin3(i) == 1 %at edge
                    temp2(1,:) = 999999;    %dummy value
                else
                    temp2(1,:) = topo(ymin3(i),xind) - topo(ymin3(i)-1,xind); %indexmask4(ymin3(i)-1,xind)    ;
                end
                if ymax3(i) == jtot %at edge
                    temp2(2,:) = 999999;    %dummy value
                else
                    temp2(2,:) = topo(ymax3(i),xind) - topo(ymax3(i)+1,xind); %indexmask4(ymax3(i)+1,xind)  ;
                end

                temp=[temp1(:)' temp2(:)'];
                if all(temp == 0) % tg3315 switched this indexing around... %temp>0) %it's internal
                    internalmask(yind,xind) = 1;
                end
            end
            externalmask = topomask - internalmask;

            % merge internal blocks
            count5 = count3;
            xmin5 = xmin3;
            xmax5 = xmax3;
            ymin5 = ymin3;
            ymax5 = ymax3;
            zmax5 = zmax3;
            zmin5 = zmin3;
            
            indexmask5 = indexmask4;
            
                       
            % SO: Not sure indexmask should be used when we have blocks on top of each other.
            % I have removed the following section for now and it seems ok, 
            % but if problems appear then this needs to be looked at again.
            %change = true;
            change = false;
            %in y first
            while change
                change=false;
                for k=1:count5
                    i=indexmask5(ymin5(k),xmin5(k));
                    if internalmask(ymin5(i),xmin5(i))
                        if internalmask(ymin5(i)-1,xmin5(i))  %can only have 1 neighbour in y direction at this stage
                            i2=indexmask5(ymin5(i)-1,xmin5(i));
                            ymin5(i)=ymin5(i2);
                            indexmask5(ymin5(i):ymax5(i),xmin5(i):xmax5(i))=i;
                            indexmask5(indexmask5>i2)=indexmask5(indexmask5>i2)-1;
                            count5=count5-1;
                            xmin5(i2)=[];
                            xmax5(i2)=[];
                            ymin5(i2)=[];
                            ymax5(i2)=[];
                            zmax5(i2)=[];
                            zmin5(i2)=[];
                            change=true;
                            %count5
                            break
                        end
                    end
                end
            end

            %make blocks
            datamean5 = zeros(size(topomask));
            dataind = zeros(size(topomask));
            for i = 1:count5
                datamean5(ymin5(i):ymax5(i),xmin5(i):xmax5(i)) = zmax5(i);
                dataind(ymin5(i):ymax5(i),xmin5(i):xmax5(i)) = i;
            end

            ymin6 = ymin5;
            ymax6 = ymax5;
            xmin6 = xmin5;
            xmax6 = xmax5;
            zmin6 = zmin5;
            zmax6 = zmax5;

            %% flip the whole y-dimension!!
            ymin6f = jtot - ymax6 + 1;
            ymax6f = jtot - ymin6 + 1;
            ymin2f = jtot - ymax2 + 1;
            ymax2f = jtot - ymin2 + 1;

            dummy = ones(count5, 5);
            blocks = [xmin5 xmax5 ymin6 ymax6 zmin5 zmax5 dummy];


            %write block file to use for test of ray intersection (without roads that will be added to blocks.inp)
            %already raise by 1 (because of roads, is done to blocks.inp.xxx later)
            %give each building a unique number, might differ from results in block2fac
            %but it doesn't matter

            blockindexmask = dataind;
            nbl = size(xmin5, 1);
            buildingindexmask = bwlabel(blockindexmask);
            buildingindexlist = zeros(1, nbl);

            xuu = xmax5;
            yuu = ymax6;
            zuu = zmax5;

            for i = 1:nbl
                buildingindexlist(i) = buildingindexmask(yuu(i), xuu(i));
            end            

            
            buildings = [xmin5  xmax5  ymin6 ymax6 zmin5+1 zmax5+1 buildingindexlist'];
            
            % Remove blocks
            overhangblockids_remove = [];
            for i=1:size(blocks,1)
                il = blocks(i,1); iu = blocks(i,2);
                jl = blocks(i,3); ju = blocks(i,4);
                kl = blocks(i,5); ku = blocks(i,6);
                
                if all(over(jl:ju, il:iu), 'all')
                    if kl > 0 && obj.zh(ku+2) <= over(jl,il) % Remove blocks less than overhang height but not at ground level
                        overhangblockids_remove = [overhangblockids_remove; i];
                    end
                end          
            end
            
            blocks(overhangblockids_remove,:) = [];
            buildings(overhangblockids_remove,:) = [];
            
            %% Shrink overhang blocks
            overhangblockids = [];           

            OH = zeros(imax, jtot);
            BH = zeros(imax, jtot);
            M = zeros(imax, jtot);
            IM = zeros(imax, jtot);
            for i=1:size(blocks,1)
                il = blocks(i,1); iu = blocks(i,2);
                jl = blocks(i,3); ju = blocks(i,4);
                kl = blocks(i,5); ku = blocks(i,6);
                M(il:iu, jl:ju) = 1;
                IM(il:iu, jl:ju) = i; % What about when blocks are on top of each other?               
                BH(il:iu, jl:ju) = max(BH(il:iu, jl:ju), ku+1);
                
                if all(over(jl:ju, il:iu), 'all')
                    OH(il:iu, jl:ju) = over(jl:ju, il:iu)' / dz;
                    if kl == 0 % Shrink lowest block in stack
                        %OH(il:iu, jl:ju) = ku+1; % Maybe simpler to just to overhang topo                       
                        blocks(i,5) = -1;
                        blocks(i,6) = -1;
                        overhangblockids = [overhangblockids; i];
%                     elseif kl > 0 && ku+1 <= over(jl,il) % Remove blocks less than overhang height but not at ground level
%                         %OH(il:iu, jl:ju) = ku+1; % Maybe simpler to just to overhang topo
%                         overhangblockids_remove = [overhangblockids_remove; i];
                    end
                end               
            end


            preprocessing.addvar(obj, 'BH', BH);
            preprocessing.addvar(obj, 'M', M);
            preprocessing.addvar(obj, 'IM', IM);
            preprocessing.addvar(obj, 'OH', OH);
            
            obj.nblocks = size(blocks, 1);
            preprocessing.addvar(obj, 'blocks', blocks)
            preprocessing.addvar(obj, 'buildings', buildings)
            preprocessing.addvar(obj, 'overhangblockids', overhangblockids)
            %preprocessing.addvar(obj, 'overhangblockids_above', overhangblockids_above)
        end
        
        function makeblocks_overlayer(obj, n)
            overlayermask = squeeze(obj.overlayermask(:,:,n));
            overlayer = squeeze(obj.overlayers(:,:,n));
            %over = obj.overhangs;
            
            layerheight = max(obj.blocks(:,6))+1; % Height of existing blocks in cells - include previous overblocks
            if n==1
                prevlayer = obj.topo;
            else
                %prevlayer = obj.topo + squeeze(obj.overlayermask(:,:,n-1)); % Actual height of the previous layer
                prevlayer = obj.layerheights(n-1) + squeeze(obj.overlayermask(:,:,n-1));
            end
            
            imax = obj.imax;
            jtot = obj.jtot;
            dz = obj.dz;
            
            
            obj.layerheights(n) = layerheight;
            % Define using obj.BH + layerheight(n-1);

            if (obj.lzstretch && (max(max(overlayer))>obj.hlin))
            	error('ERROR: hlin must be greater than the maximum vertical building index')
            end

            maxnrblocks = sum(overlayermask(:)); %allocate arrays with maximum size they can possibly have, reduce size later
            xmin = zeros(maxnrblocks,1); %store lower x bound of blocks
            xmax = zeros(maxnrblocks,1); %store upper x bound of blocks
            ymin = zeros(maxnrblocks,1); %store lower y bound of blocks
            ymax = zeros(maxnrblocks,1); %store upper y bound of blocks
            zmin = zeros(maxnrblocks,1); %store lower z bound of blocks
            zmax = zeros(maxnrblocks,1); %store upper z bound of blocks

            blchecked=zeros(size(overlayermask)); %mask of blocks which have already been checked
            indexmask=zeros(size(overlayermask)); %mask of the indeces of the (new) blocks
                
            count = 0;
            for i = 1:imax %loop over all x
                xminl = i;
                xmaxl = i;
                j = 1;
                while j <= jtot %loop along j
                    if overlayermask(j, i) == 0 %not a building
                        j = j + 1; %check next j
                        continue
                    else
                        yminl = j;
                        %heightgradient = diff(topo(j:end, i));
                         % If there are blocks below then we need to divide the block. i.e. when overlayer = topo.
                             heightgradient = diff(overlayer(j:end, i));
                             if isempty(heightgradient) % if at the end of the y-direction (je)
                                 heightchangey = j;
                             elseif all(heightgradient == 0) % same block/ floor until the end of the domain
                                 heightchangey = jtot;
                                 %heightchangey = min(j+obj.maxsize, jtot);
                             else
                                 heightchangey = find(heightgradient ~=0 , 1) + j - 1;  %last cell with same height as j (i.e. there is a height change betwen this and the next cell
                                 %heightchangey = min(obj.maxsize, find(heightgradient ~=0 , 1)) + j - 1;  %last cell with same height as j (i.e. there is a height change betwen this and the next cell                               
                             end
                             
                             %If there is a block below somewhere in the y direction but not at the current location
                             if any(prevlayer(j:heightchangey,i) == layerheight)
                                 if all(prevlayer(j:heightchangey,i) == layerheight) || all(prevlayer(j:heightchangey,i) ~= layerheight) 
                                     % No change of situation
                                 
                                 else % Need to split at first change
                                    if prevlayer(j,i) == layerheight % If current position is over a block    
                                        pos = find(prevlayer(j:heightchangey,i) ~= layerheight, 1)-1; % last cell that doesn't have a block below - subtract one here and not above because diff is a forward difference.
                                        heightchangey = min(heightchangey, pos + j - 1);
                                    else
                                       pos = find(prevlayer(j:heightchangey,i) == layerheight, 1)-1; % last cell that doesn't have a block below - subtract one here and not above because diff is a forward difference.
                                       heightchangey = min(heightchangey, pos + j - 1);
                                    end
                                 end
                             end

                         % If at any point from ymin to ymax there is a difference in height with the adjacent cells in x direction,
                         % or there is an overhang adjacent 
                        for jj = yminl+1:heightchangey
                            if (i == 1)
                                if ((any(overlayer(jj-1:jj,i+1) > 0)) && (overlayer(jj, i+1) ~= overlayer(jj-1,i+1))) %Need to add condition
                                    heightchangey = jj - 1;
                                    break
                                end
                            elseif (i == imax)
                                if (any(overlayer(jj-1:jj,i-1) > 0) && (overlayer(jj,i-1)~= overlayer(jj-1,i-1))) %Need to add condition
                                    heightchangey = jj-1;
                                    break
                                end
                            else
                                if ((any(overlayer(jj-1:jj,i-1) > 0) && (overlayer(jj,i-1) ~= overlayer(jj-1,i-1))) || ((any(overlayer(jj-1:jj,i+1) > 0)) && (overlayer(jj,i+1)~= overlayer(jj-1,i+1)))) ...
                                || ((prevlayer(jj,i-1) ~= layerheight && prevlayer(jj-1,i-1) == layerheight) || (prevlayer(jj,i-1) == layerheight && prevlayer(jj-1,i-1) ~= layerheight) || (prevlayer(jj,i+1) ~= layerheight && prevlayer(jj-1,i+1) == layerheight) || (prevlayer(jj,i+1) == layerheight && prevlayer(jj-1,i+1) ~= layerheight))
                                    heightchangey = jj-1; % overwrite heightchangey as we need to truncate block earlier!
                                    break
                                end
                            end
                        end

                        ymaxl = heightchangey;  % end of the current block (either floor or block with different height comes next)
                        ztemp = zeros(6,1);
                        %ztempoverhangs = zeros(6,1);

                        % tg3315 changed from commented below as can have yminl==1 and
                        % ymaxl==nj
                        if yminl == 1
                            ztemp(2,1) = NaN;
                            %ztempoverhangs(2,1) = NaN;
                        else
                            ztemp(2,1) = overlayer(yminl-1,xminl); % Topo height in -y direction 
                            %ztempoverhangs(2,1) = over(yminl-1,xminl); % Overhang height in -y direction
                        end

                        if xminl == 1
                            ztemp(3,1) = NaN;
                            %ztempoverhangs(3,1) = NaN;
                        else
                            ztemp(3,1) = overlayer(yminl, xminl-1); % Topo height in -x direction
                            %ztempoverhangs(3,1) = over(yminl, xminl-1); % Overhang height in -x direction
                        end

                        if xmaxl == imax
                            ztemp(4,1) = NaN;
                            %ztempoverhangs(4,1) = NaN;
                        else
                            ztemp(4,1) = overlayer(yminl, xminl+1); % Topo height in +x direction 
                            %ztempoverhangs(4,1) = over(yminl, xminl+1); % Overhang height in +x direction 
                        end

                        if ymaxl == jtot
                            ztemp(5,1) = NaN;
                            %ztempoverhangs(5,1) = NaN;
                        else
                            %ztemp(5,1) = overlayer(yminl+1, xminl); % Topo height in +y direction
                            ztemp(5,1) = overlayer(ymaxl+1, xminl); % Height of next block. Not sure about this...
                            %ztempoverhangs(5,1) = over(yminl+1, xminl); % Overhang height in +y direction 
                        end
                        
%                         if ymaxl ~= jtot
%                             zcutmax = topo(ymaxl+1, xminl); % Height of next block
%                         else
%                             zcutmax = 0;
%                         end

                         zcuttemp = sort(ztemp(ztemp < overlayer(yminl, xminl) & ztemp > 0));
                         %zcuttemp = sort([ztemp(ztemp < topo(yminl, xminl) & ztemp > 0); zcutmax]);
                        %zcuttemp = sort([ztemp(ztemp < overlayer(yminl, xminl) & ztemp > 0); ...
                                         %ztempoverhangs(ztempoverhangs < overlayer(yminl, xminl) & ztempoverhangs > 0); ...
                                         %over(yminl, xminl)]);
                                     
                        zcuttemp = unique(zcuttemp);                        
                        zcut = [0; zcuttemp; overlayer(yminl,xminl)];
                        zcut = unique(zcut);
                        for kc = 1:size(zcut, 1)-1
                            count = count+1;
                            blchecked(yminl:ymaxl,xminl:xmaxl) = blchecked(yminl:ymaxl,xminl:xmaxl)+1;
                            indexmask(yminl:ymaxl,xminl:xmaxl) = count;
                            xmin(count) = xminl;
                            xmax(count) = xmaxl;
                            ymin(count) = yminl;
                            ymax(count) = ymaxl;
                            %if overhangs(yminl,xminl) > 0
                                %zmin(count) = 0;
                                %zmax(count) = 0;
                            %else
                            if zcut(kc,1) == 0 % kc==1
                                zmin(count) = 0;
                            else
                                %disp('zcut(kc,1) ~= 0')
                                zmin(count) = zcut(kc,1) / dz;
                            end
                            zmax(count) = zcut(kc+1,1) / dz - 1;
                            %zmin(count)
                            if zmin(count) > zmax(count)
                                disp('zmin > zmax')
                            end
                            %end
                        end
                        j = heightchangey + 1;  % move to index after current block
                    end
                end
            end

            %shorten arrays
            ymin((count+1):end)=[];
            ymax((count+1):end)=[];
            xmin((count+1):end)=[];
            xmax((count+1):end)=[];
            zmin((count+1):end)=[];
            zmax((count+1):end)=[];

            %% add y slices of same dimension, aggregate along x
            %keep old xmax etc.. to keep better track whats happening
            xmax2 = xmax; ymax2 = ymax; zmax2 = zmax;
            xmin2 = xmin; ymin2 = ymin; zmin2 = zmin;
            dsize = 1;
            sizeold = count;
            while dsize > 0  %do as long as there are unmerged blocks % tg3315 need to also check here that you do not merge if they have a block above...!
                i = 1;
                while 1  %do along x
                    a = xmax2(i);  %upper x index of this block
                    bv = find(xmin2 == (a+1));  %all blocks with a lower x bound 1 bigger than this blocks upper x bound
                    b2 = bv(ymin2(bv) == ymin2(i)); %all of those blocks with also the same lower y bound
                    b3 = b2(zmin2(b2) == zmin2(i));
                    if ~isempty(b3) 
                        %if length(b3) == 1
                        if (all(ymax2(b3) == ymax2(i)) && all(zmax2(b3) == zmax2(i)) && overlayer(ymin2(b3),xmin2(b3)) == overlayer(ymin2(i), xmin2(i))) ...
                            && ~((prevlayer(ymin2(b3), xmin2(b3)) ~= layerheight && prevlayer(ymin2(i), xmin2(i)) == layerheight) || (prevlayer(ymin2(b3), xmin2(b3)) == layerheight && prevlayer(ymin2(i), xmin2(i)) ~= layerheight))%&& all(zmin2(b2)==zmin2(i)) %if they also have the same upper y bound and the same height
                            % additional check to make sure we do not merge blocks in a way that causes internal-external facets at this point
                            %'Merging'
                            if (overlayer(max(1,ymin2(b3)-1),xmax2(b3)) == overlayer(max(1,ymin2(i)-1),xmax2(i))) && (overlayer(min(ymax2(b3)+1, jtot),xmax2(b3)) == overlayer(min(ymax2(i)+1, jtot), xmax2(i)))
                                xmax2(i) = xmax2(b3); %merge
                                xmax2(b3) = []; ymax2(b3) = []; zmax2(b3)  =[]; xmin2(b3) = []; ymin2(b3) = []; zmin2(b3) = []; %remove the just merged block from list of blocks
                            end
                        end
                        %end
                    end
                    i = i+1; %check furhter along x
                    if i >= length(xmin2)
                        break
                    end
                end
                dsize = sizeold-length(xmax2);
                sizeold = length(xmax2);
            end
            count2 = sizeold;            

            %make fields again
            datamean1 = zeros(size(overlayermask));
            datamean2 = zeros(size(overlayermask));
            %datamean3 = zeros(size(topomask));
            indexmask2 = zeros(size(overlayermask));

            %make new matrices
            for i=1:count
                datamean1(ymin(i):ymax(i),xmin(i):xmax(i)) = zmax(i);
            end
            for i=1:count2
                indexmask2(ymin2(i):ymax2(i),xmin2(i):xmax2(i)) = i;
                datamean2(ymin2(i):ymax2(i),xmin2(i):xmax2(i)) = zmax2(i);
            end

            %% split slices up according to rules in x and y, do z later
            xmin3 = zeros(maxnrblocks,1);
            xmax3 = zeros(maxnrblocks,1);
            ymin3 = zeros(maxnrblocks,1);
            ymax3 = zeros(maxnrblocks,1);
            zmin3 = zeros(maxnrblocks,1);
            zmax3 = zeros(maxnrblocks,1);

            xmin3(1:count2) = xmin2;
            xmax3(1:count2) = xmax2;
            ymin3(1:count2) = ymin2;
            ymax3(1:count2) = ymax2;
            zmin3(1:count2) = zmin2;
            zmax3(1:count2) = zmax2;

            
            if obj.lEB
            %if false
                indexmask3 = indexmask2;
                change = true;
                count3 = count2;

                maxx = 0;
                cmax = 0;

                while change  %do until there is no changes anymore
                    change = false;
                    for c = 1:count3 %check if all blocks have same dimension as neighbour
                        %only have to check x-1 and x+1, since there can't be a
                        %neightbour on y-1 or y+1 (due to the inital slicing along y)
                        %index of blocks to the left

                        if ~(xmin3(c) == 1) %on left domain boundary, don't need to check for neighbouring blocks
                            illb = indexmask3(ymin3(c):ymax3(c),xmin3(c)-1);
                            iillb = find(illb>0,1);
                            if ~isempty(illb)
                                ilb = illb(iillb);
                            else
                                ilb = 0;
                            end
                        else
                            ilb = 0;
                        end

                        if ~(xmax3(c) == imax) %on right domain boundary, don't check for neighbouring blocsk
                            irrb = indexmask3(ymin3(c):ymax3(c),xmax3(c)+1);
                            iirrb = find(irrb>0,1);
                            if ~isempty(irrb)
                                irb = irrb(iirrb);
                            else
                                irb = 0;
                            end
                        else
                            irb = 0;
                        end
                        %   irt=indexmask(ymin(c),xmax(c)+1);
                        if (irb == 0) %no block to the right
                            yminr = ymin3(c); %use y coordinates of this block
                            ymaxr = ymax3(c);
                        else
                            yminr = ymin3(irb);  %use y coordinates of block to the right
                            ymaxr = ymax3(irb);
                        end

                        if (ilb == 0) %no block to the left
                            ymaxl = ymax3(c); %use y coordinates of this block
                            yminl = ymin3(c);
                        else
                            ymaxl = ymax3(ilb); %use y coordinates of block to the left
                            yminl = ymin3(ilb);
                        end

                        if xmin(3) > maxx
                            maxx = xmin3(c);
                            %maxx/ni;
                        end

                        if c/count3 > cmax
                            cmax=c/count3;
                        end

                        %
                        %

                        if (ymaxl < ymax3(c))% & zmax3(ilb)==zmax3(c)  %this block extends further (up) than left neighbour, this block has to be split
                            %x  x  x  x      %x  x  x  x
                            %x    [c] x      %x    [n] x
                            %x [l][c] x      %x [l][c] x
                            %x [l][c] x ===> %x [l][c] x
                            %x [l][c] x      %x [l][c] x
                            %x  x  x  x      %x  x  x  x
                            change = true;            %keep checking
                            count3 = count3+1;        %add new block
                            ymin3(count3) = ymaxl+1;  %new block starts 1 above xmax of left neighbour
                            ymax3(count3) = ymax3(c); %new block ends at xmax of this block
                            xmin3(count3) = xmin3(c); %new block has same x coordinates as this block
                            xmax3(count3) = xmax3(c);
                            zmin3(count3) = zmin3(c); %keep the same heights
                            zmax3(count3) = zmax3(c);
                            indexmask3(ymin3(count3):ymax3(count3),xmin3(count3):xmax3(count3)) = count3;
                            ymax3(c) = ymaxl;         %shorten this block
                            break                     %continue forloop from the start

                        elseif (yminl > ymin3(c))% & zmin3(ilb)==zmin3(c)  %this block extends further (down) than left neighbour, this block has to be split
                            %x  x  x  x      %x  x  x  x
                            %x [l][c] x      %x [l][c] x
                            %x [l][c] x      %x [l][c] x
                            %x [l][c] x ===> %x [l][c] x
                            %x    [c] x      %x    [n] x
                            %x  x  x  x      %x  x  x  x
                            change = true; %keep checking
                            count3 = count3+1; %add new block
                            ymin3(count3) = ymin3(c); %new block starts at same ymin as this
                            ymax3(count3) = yminl-1;  %new block ends 1 below ymin of left neighbour
                            xmin3(count3) = xmin3(c); %new block has same x coordinates as this block
                            xmax3(count3) = xmax3(c);
                            zmin3(count3) = zmin3(c); %keep the same heights
                            zmax3(count3) = zmax3(c);
                            indexmask3(ymin3(count3):ymax3(count3),xmin3(count3):xmax3(count3)) = count3;  %overwrite the index in the indexmask with the new one
                            ymin3(c)=yminl;         %shorten this block
                            break                   %continue forloop from the start

                        elseif (yminr > ymin3(c))% & zmin3(irb)==zmin3(c) %this block extends further (down) than reight neighbour, this block has to be split
                            %x  x  x  x      %x  x  x  x
                            %x [c][r] x      %x [c][r] x
                            %x [c][r] x      %x [c][r] x
                            %x [c][r] x ===> %x [c][r] x
                            %x [c]    x      %x [n]    x
                            %x  x  x  x      %x  x  x  x
                            change = true;            %keep checking
                            count3 = count3+1;        %add new block
                            ymin3(count3) = ymin3(c); %new block starts at same ymin as this
                            ymax3(count3) = yminr-1;  %new block ends 1 below ymin of reight neighbour
                            xmin3(count3) = xmin3(c); %new block has same x coordinates as this block
                            xmax3(count3) = xmax3(c);
                            zmin3(count3) = zmin3(c); %keep the same heights
                            zmax3(count3) = zmax3(c);
                            indexmask3(ymin3(count3):ymax3(count3),xmin3(count3):xmax3(count3)) = count3; %overwrite the index in the indexmask with the new one
                            ymin3(c) = yminr;         %shorten this block
                            break                     %continue for-loop from the start

                        elseif (ymaxr < ymax3(c))% & zmax3(irb)==zmax3(c) %this block extends further (up) than reight neighbour, this block has to be split
                            %x  x  x  x      %x  x  x  x
                            %x [c]    x      %x [n]    x
                            %x [c][r] x      %x [c][r] x
                            %x [c][r] x ===> %x [c][r] x
                            %x [c][r] x      %x [c][r] x
                            %x  x  x  x      %x  x  x  x
                            change = true;             %keep checking
                            count3 = count3+1;         %add new block
                            ymin3(count3) = ymaxr+1 ;  %new block starts 1 above xmax of reight neighbour
                            ymax3(count3) = ymax3(c);  %new block ends at xmax of this block
                            xmin3(count3) = xmin3(c);  %new block has same x coordinates as this block
                            xmax3(count3) = xmax3(c);
                            zmin3(count3) = zmin3(c); %keep the same heights
                            zmax3(count3) = zmax3(c);
                            indexmask3(ymin3(count3):ymax3(count3),xmin3(count3):xmax3(count3)) = count3;
                            ymax3(c) = ymaxr; %shorten this block
                            break %continue forloop from the start

                        end

                    end


                end
                ymin3((count3+1):end) = [];
                ymax3((count3+1):end) = [];
                xmin3((count3+1):end) = [];
                xmax3((count3+1):end) = [];
                zmin3((count3+1):end) = [];
                zmax3((count3+1):end) = [];


            else
                count3 = count2;
                xmin3 = xmin2;
                xmax3 = xmax2;
                ymin3 = ymin2;
                ymax3 = ymax2;
                zmin3 = zmin2;
                zmax3 = zmax2;

            end

            %% aggregate building internal blocks  (only rectangles)
            %  1 2 3 4     1 2 3 4
            %  5 6 7 8     5i 6 6 7
            %  9 a b c  -> 8 6 6 9
            %  d e f g     a b c d

            %use that all neighbours have the same size in 1 dimension % tg3315 mot
            %true now...

            datamean4 = zeros(size(overlayermask));
            indexmask4 = zeros(size(overlayermask));
            internalmask = zeros(size(overlayermask));

            for i = 1:count3
                indexmask4(ymin3(i):ymax3(i), xmin3(i):xmax3(i)) = i;
                datamean4(ymin3(i):ymax3(i), xmin3(i):xmax3(i)) = zmax3(i);
            end

            for k = 1:count3
                i = indexmask4(ymin3(k), xmin3(k)); %index of this block (k = i ?might be unnecessary)
                %check if it is internal
                xind = xmin3(i):xmax3(i);
                yind = ymin3(i):ymax3(i);
                temp1=zeros(length(yind),2);
                temp2=zeros(2,length(xind));
                if xmin3(i) == 1 %at edge
                    temp1(:,1) = 999999;    %dummy value
                else
                    temp1(:,1) = overlayer(yind,xmin3(i)) - overlayer(yind,xmin3(i)-1); % tg3315 changed so do not get internal blocks with adjacent different sizes %indexmask4(yind, xmin3(i)-1);
                end
                if xmax3(i) == imax %at edge
                    temp1(:,2) = 999999;
                else
                    temp1(:,2) = overlayer(yind,xmax3(i)) - overlayer(yind,xmax3(i)+1); %indexmask4(yind, xmax3(i)+1);
                end
                if ymin3(i) == 1 %at edge
                    temp2(1,:) = 999999;    %dummy value
                else
                    temp2(1,:) = overlayer(ymin3(i),xind) - overlayer(ymin3(i)-1,xind); %indexmask4(ymin3(i)-1,xind)    ;
                end
                if ymax3(i) == jtot %at edge
                    temp2(2,:) = 999999;    %dummy value
                else
                    temp2(2,:) = overlayer(ymax3(i),xind) - overlayer(ymax3(i)+1,xind); %indexmask4(ymax3(i)+1,xind)  ;
                end

                temp=[temp1(:)' temp2(:)'];
                if all(temp == 0) % tg3315 switched this indexing around... %temp>0) %it's internal
                    internalmask(yind,xind) = 1;
                end
            end
            externalmask = overlayermask - internalmask;

            % merge internal blocks
            count5 = count3;
            xmin5 = xmin3;
            xmax5 = xmax3;
            ymin5 = ymin3;
            ymax5 = ymax3;
            zmax5 = zmax3;
            zmin5 = zmin3;
            
            indexmask5 = indexmask4;
            change = true;
            %change = false; % SO changed to false so that 'internal' blocks are not removed when setting maxsize
            %in y first
            while change
                change=false;
                for k=1:count5
                    i=indexmask5(ymin5(k),xmin5(k));
                    if internalmask(ymin5(i),xmin5(i))
                        if internalmask(ymin5(i)-1,xmin5(i))  %can only have 1 neighbour in y direction at this stage
                            i2=indexmask5(ymin5(i)-1,xmin5(i));
                            ymin5(i)=ymin5(i2);
                            indexmask5(ymin5(i):ymax5(i),xmin5(i):xmax5(i))=i;
                            indexmask5(indexmask5>i2)=indexmask5(indexmask5>i2)-1;
                            count5=count5-1;
                            xmin5(i2)=[];
                            xmax5(i2)=[];
                            ymin5(i2)=[];
                            ymax5(i2)=[];
                            zmax5(i2)=[];
                            zmin5(i2)=[];
                            change=true;
                            %count5
                            break
                        end
                    end
                end
            end

            %make blocks
            datamean5 = zeros(size(overlayermask));
            dataind = zeros(size(overlayermask));
            for i = 1:count5
                datamean5(ymin5(i):ymax5(i),xmin5(i):xmax5(i)) = zmax5(i);
                dataind(ymin5(i):ymax5(i),xmin5(i):xmax5(i)) = i;
            end

            ymin6 = ymin5;
            ymax6 = ymax5;
            xmin6 = xmin5;
            xmax6 = xmax5;
            zmin6 = zmin5;
            zmax6 = zmax5;

            %% flip the whole y-dimension!!
            ymin6f = jtot - ymax6 + 1;
            ymax6f = jtot - ymin6 + 1;
            ymin2f = jtot - ymax2 + 1;
            ymax2f = jtot - ymin2 + 1;

            dummy = ones(count5, 5);
            overblocks = [xmin5 xmax5 ymin6 ymax6 zmin5+layerheight zmax5+layerheight dummy];


            %write block file to use for test of ray intersection (without roads that will be added to blocks.inp)
            %already raise by 1 (because of roads, is done to blocks.inp.xxx later)
            %give each building a unique number, might differ from results in block2fac
            %but it doesn't matter

            blockindexmask = dataind;
            nbl = size(xmin5, 1);
            buildingindexmask = bwlabel(blockindexmask);
            buildingindexlist = zeros(1, nbl);

            xuu = xmax5;
            yuu = ymax6;
            zuu = zmax5;

            for i = 1:nbl
                buildingindexlist(i) = buildingindexmask(yuu(i), xuu(i));
            end            

            
            
            obj.blocks = [obj.blocks; overblocks];
            obj.overblocks = [obj.overblocks; overblocks];
            obj.nblockslayers(n) = size(overblocks,1);
            obj.noverblocks = obj.noverblocks + size(overblocks,1);
            
            obj.nblocks = obj.nblocks + size(overblocks,1); % Don't update til the end     
            size(obj.blocks)
            obj.nblocks
            obj.overbuildings = [obj.overbuildings; [xmin5  xmax5  ymin6 ymax6 zmin5+1+layerheight zmax5+1+layerheight buildingindexlist']];
            %obj.buildings = [obj.buildings; overbuildings];
            %obj.nbuildings = obj.nbuildings+size(overbuildings,1);
            %preprocessing.addvar(obj, 'noverblocks', size(overblocks,1));
            %preprocessing.addvar(obj, 'overblocks', overblocks)
            %preprocessing.addvar(obj, 'buildings', [xmin5  xmax5  ymin6 ymax6 zmin5+1 zmax5+1 buildingindexlist'])
            %preprocessing.addvar(obj, 'overhangblockids', overhangblockids)
            %preprocessing.addvar(obj, 'overhangblockids_above', overhangblockids_above)
            
        end

        function block2fac(obj)
            imax = obj.imax;
            jtot = obj.jtot;
            blocks = obj.blocks;
            nblocks = obj.nblocks;
            blocks(:, [5,6]) = blocks(:,[5,6]) + 1; % indices for k start at zero?
            %obj.blocks = blocks;
            %% create Mask-matrix
            % this new mask is in x,y coordinates, not y,x coordinates as before
            M = obj.M;%zeros(imax, jtot);
            IM = obj.IM;%zeros(imax , jtot);
            BH = obj.BH;%zeros(imax, jtot);
            OH = obj.OH;
%             for i = 1:size(blocks,1)-obj.noverblocks
%                 xl = blocks(i,1);
%                 xu = blocks(i,2);
%                 yl = blocks(i,3);
%                 yu = blocks(i,4);
%                 M(xl:xu, yl:yu) = 1;
%                 IM(xl:xu, yl:yu) = i; % What about when blocks are on top of each other?
%                 BH(xl:xu, yl:yu) = max(BH(xl:xu, yl:yu), blocks(i,6));
%             end
%             preprocessing.addvar(obj, 'BH', BH)
            
            top = 1; west = 2; east = 3; north = 4; south=5; bot = 6;

            il = 1; iu = 2;
            jl = 3; ju = 4;
            kl = 5; ku = 6;

            nfcts = (nblocks-obj.noverblocks) * 6;
            facets = int32(zeros(nfcts, 10));
            %obj.overhangblockids
            overhangfacetids = [];
            %overhangfacetids_above = [];
            for j = 1:nblocks-obj.noverblocks
                for k = top:bot
                    i = (j - 1) * 6 + k;
                    if any(ismember(obj.overhangblockids, j))
                        overhangfacetids = [overhangfacetids; i];
                        %continue
                    end
                    facets(i, 1) = k; % orientation
                    % for all orientations apart from bottom, wall type =
                    % blk(j, 6+5) = blk(j,11) = 1. For bottom, wall type =
                    % blk(j, 6+6) = blk(j,12) = 0. This depends on what
                    % bllk was defined as before, and needs to be changed
                    % so that it is easier to change wall type.
                    facets(i,2) = blocks(j, 6 + min(k,5)); % wall type
                    facets(i,3) = j; % blockid
                    switch(k)
                        case top
                            facets(i, 5 + il) = blocks(j, il);
                            facets(i, 5 + jl) = blocks(j, jl);
                            facets(i, 5 + kl) = blocks(j, ku) + 1;
                            facets(i, 5 + iu) = blocks(j, iu) + 1;
                            facets(i, 5 + ju) = blocks(j, ju) + 1;
                            facets(i, 5 + ku) = blocks(j, ku) + 1;
                        case west
                            facets(i, 5 + il) = blocks(j, il);
                            facets(i, 5 + jl) = blocks(j, jl);
                            facets(i, 5 + kl) = blocks(j, kl);
                            facets(i, 5 + iu) = blocks(j, il);
                            facets(i, 5 + ju) = blocks(j, ju) + 1;
                            facets(i, 5 + ku) = blocks(j, ku) + 1;
                        case east
                            facets(i, 5 + il) = blocks(j, iu) + 1;
                            facets(i, 5 + jl) = blocks(j, jl);
                            facets(i, 5 + kl) = blocks(j, kl);
                            facets(i, 5 + iu) = blocks(j, iu) + 1;
                            facets(i, 5 + ju) = blocks(j, ju) + 1;
                            facets(i, 5 + ku) = blocks(j, ku) + 1;
                        case north
                            facets(i, 5 + il) = blocks(j, il);
                            facets(i, 5 + jl) = blocks(j, ju) + 1;
                            facets(i, 5 + kl) = blocks(j, kl);
                            facets(i, 5 + iu) = blocks(j, iu) + 1;
                            facets(i, 5 + ju) = blocks(j, ju) + 1;
                            facets(i, 5 + ku) = blocks(j, ku) + 1;
                        case south
                            facets(i, 5 + il) = blocks(j, il);
                            facets(i, 5 + jl) = blocks(j, jl);
                            facets(i, 5 + kl) = blocks(j, kl);
                            facets(i, 5 + iu) = blocks(j, iu) + 1;
                            facets(i, 5 + ju) = blocks(j, jl);
                            facets(i, 5 + ku) = blocks(j, ku) + 1;
                        case bot
                            facets(i, 5 + il) = blocks(j, il);
                            facets(i, 5 + jl) = blocks(j, jl);
                            facets(i, 5 + kl) = blocks(j, kl);
                            facets(i, 5 + iu) = blocks(j, iu) + 1;
                            facets(i, 5 + ju) = blocks(j, ju) + 1;
                            facets(i, 5 + ku) = blocks(j, kl);
                    end
                end
            end          
            preprocessing.addvar(obj, 'overhangfacetids', overhangfacetids)
            
            intern = zeros(nfcts, 2);
            j = 1;
            
%             for i = 1:nfcts
%                 switch(facets(i, 1))
%                     %only one test necessary, since the whole edge is internal, or not
%                     %remember that facets stores facet coordinates, not block coordinates
%                     case 1
%                         ids = find(ismember(facets(:,[1,6:9]), facets(i,[1,6:9]),'rows')); 
%                         if length(ids) ~= 1 % there are top facets immediately below or above.
%                               kcurrent = facets(i,5+kl);
%                               ks = facets(ids, 5+kl);
%                               if kcurrent == min(ks) % If it is the lowest in the stack
%                                   facets(i,5) = 1;
%                                   intern(j,1) = facets(i, 3);
%                                   intern(j,2) = IM(facets(i, 5 + il), facets(i, 5 + jl));
%                                   j = j + 1;
%                               end
%                         end
%                     case 2 % west
%                         if facets(i, 5+il)-1 >= 1 %not at the domain edge
%                             %if M(facets(i, 5 +il)-1 , facets(i, 5+jl)) % if there is a block in the -y direction
%                                 
%                                ids = find(ismember(facets(:,6:11), facets(i,6:11),'rows'))'; % Find facets with exactly the same coordinates
%                                if length(ids) > 1
%                                    for id = ids%(2:end) % remove duplicates
%                                        facets(id,5) = 1;
%                                        intern(j,1) = facets(id, 3);
%                                        %intern(j,2) = IM(facets(id, 5+il)-1 ,facets(id, 5+jl));
%                                        intern(j,2) = IM(facets(id, 5+il) ,facets(id, 5+jl));
%                                        j = j + 1;
%                                    end
%                                end
%                                
% %                                 if ismember(i, overhangfacetids_above)
% %                                     if M(facets(id, 5+il)-1-1, facets(id, 5+jl)) % If there is a block adjacent in -x direction (need to be careful if there is terrain)                                  
% %                                         facets(i,5) = 1;
% %                                         intern(j,1) = facets(i, 3);
% %                                         intern(j,2) = IM(facets(i, 5+il)-1, facets(i, 5+jl));
% %                                         j = j + 1;
% %                                     end
% %                                 end
%                                
%                             %end
%                         end
%                     case 3 % east
%                         if facets(i, 5 + iu) <= imax %not at the domain edge
%                             %if M(facets(i, 5+iu), facets(i, 5+jl)) % if building
% 
%                                 ids = find(ismember(facets(:,6:11), facets(i,6:11),'rows'))';
% 
%                                 if length(ids) > 1
%                                     %i
%                                     %disp('Duplicate facet')
%                                     for id = ids%(2:end)
%                                         facets(id,5) = 1;
%                                         intern(j,1) = facets(id, 3);
%                                         %intern(j,2) = IM(facets(id, 5+iu), facets(id, 5+jl));
%                                         intern(j,2) = IM(facets(id, 5+iu)-1, facets(id, 5+jl));
%                                         j = j + 1;
%                                     end
%                                 end
%                                 
% %                                 if ismember(i, overhangfacetids_above)
% %                                     if M(facets(id, 5+iu)+1, facets(id, 5+jl)) % If there is a block adjacent in +x direction (need to be careful if there is terrain)                                  
% %                                         facets(i,5) = 1;
% %                                         intern(j,1) = facets(i, 3);
% %                                         intern(j,2) = IM(facets(i, 5+iu), facets(i, 5+jl));
% %                                         j = j + 1;
% %                                     end
% %                                 end
%                                 
%                             %end
%                         end
%                     case 4 % north
%                         if facets(i, 5 + ju) <= jtot %not at the domain edge
%                             %if M(facets(i, 5+il), facets(i, 5+ju))
%                                 ids = find(ismember(facets(:,6:11), facets(i,6:11),'rows'))';
%                                 if length(ids) > 1
%                                     %i
%                                     %disp('Duplicate facet')
%                                     for id = ids%(2:end)
%                                         facets(id,5) = 1;
%                                         intern(j,1) = facets(id, 3);
%                                         %intern(j,2) = IM(facets(id, 5+il), facets(id, 5+ju));
%                                         intern(j,2) = IM(facets(id, 5+il), facets(id, 5+ju)-1);
%                                         j = j + 1;
%                                     end
%                                 end
%                                  
% %                                 if ismember(i, overhangfacetids_above)
% %                                     if M(facets(i, 5+il), facets(i, 5+ju)+1) % If there is a block adjacent in +y direction (need to be careful if there is terrain)                                  
% %                                         facets(i,5) = 1;
% %                                         intern(j,1) = facets(i, 3);
% %                                         intern(j,2) = IM(facets(i, 5+il), facets(i, 5+ju));
% %                                         j = j + 1;
% %                                     end
% %                                 end
%                                 
%                             %end
%                         end
%                     case 5 % south
%                         if facets(i, 5 + jl) - 1 >= 1 %not at the domain edge
%                             %if M(facets(i, 5+il), facets(i, 5+jl)-1)
%                                 ids = find(ismember(facets(:,6:11), facets(i,6:11),'rows'))';
%                                 if length(ids) > 1
%                                     %i
%                                     %disp('Duplicate facet')
%                                     for id = ids%(2:end)
%                                         facets(id, 5) = 1;
%                                         intern(j, 1) = facets(id, 3);
%                                         %intern(j, 2) = IM(facets(id, 5+il), facets(id, 5+jl)-1);
%                                         intern(j, 2) = IM(facets(id, 5+il), facets(id, 5+jl));
%                                         j = j + 1;
%                                     end
%                                 end
%                                 
% %                                 if ismember(i, overhangfacetids_above)
% %                                     if M(facets(id, 5+il), facets(id, 5+jl)-1-1) % If there is a block adjacent in -y direction (need to be careful if there is terrain)                                  
% %                                         facets(i,5) = 1;
% %                                         intern(j,1) = facets(i, 3);
% %                                         intern(j,2) = IM(facets(i, 5+il), facets(i, 5+jl)-1);
% %                                         j = j + 1;
% %                                     end
% %                                 end
%                                 
%                             %end
%                         end
% %                      case 6
% %                          if all(obj.overhangs(facets(i,5+jl):facets(i,5+ju)-1, facets(i,5+il):facets(i,5+iu)-1))
% %                              ids = find(ismember(facets(:,[1,6:9]), facets(i,[1,6:9]),'rows'));
% %                              if length(ids) ~= 1 % there are bottom facets immediately below or above
% %                                  kcurrent = facets(i,5+kl);
% %                                  ks = facets(ids, 5+kl);
% %                                  if kcurrent == max(ks) % Highest bottom facet in the stack: remove
% %                                      facets(i,5) = 1;
% %                                      intern(j,1) = facets(i, 3);
% %                                      intern(j,2) = IM(facets(i, 5 + il), facets(i, 5 + jl));
% %                                      j = j + 1;
% %                                  end
% %                              end
% %                          else % Remove if no overhang
% %                              facets(i, 5) = 1;
% %                              intern(j, 1) = facets(i, 3);
% %                              intern(j,2) = IM(facets(i, 5 + il), facets(i, 5 + jl));
% %                              j = j + 1;
% %                          end
%                 end
%             end
                       
            for i = 1:nfcts
                if any(ismember(overhangfacetids, i))
                    if (facets(i, 1) ~= 1)
                        facets(i,5) = 1;
                        intern(j,1) = facets(i, 3);
                        intern(j,2) = IM(facets(i, 5 + il), facets(i, 5 + jl));
                        j = j + 1;
                    else
                        facets(i,5+kl) = 0;
                        facets(i,5+ku) = 0;
                    end
                    continue
                end
                
                switch(facets(i, 1))
                    %only one test necessary, since the whole edge is internal, or not
                    %remember that facets stores facet coordinates, not block coordinates
                    case 1
                        if facets(i,5+ku)-1 ~= BH(facets(i,5+il),facets(i,5+jl)) % If not at the top of the geometry
%                          ids = find(ismember(facets(:,[1,6:9]), facets(i,[1,6:9]),'rows'));
%                          if length(ids) ~= 1
%                              kcurrent = facets(i,5+kl);
%                              ks = facets(ids, 5+kl);
%                              if kcurrent ~= max(ks) % There are lower top facets in the stack, so remove
                                 facets(i,5) = 1;
                                 intern(j,1) = facets(i, 3);
                                 intern(j,2) = IM(facets(i, 5 + il), facets(i, 5 + jl));
                                 j = j + 1;
%                              end
%                          end
                         end
                    case 2 % west
                        if facets(i, 5 + il) - 1 >= 1 %not at the domain edge
                            % If there is a block in the -x direction
                            if M(facets(i, 5 + il) - 1 , facets(i, 5 + jl)) % There is a block in the -x direction
                                % If block in -x direction is higher than current block or the current facet is lower than the next block, and the current facet is higher than any overhang.
                                if (BH(facets(i, 5+il)-1, facets(i, 5+jl)) >= BH(facets(i, 5+il), facets(i, 5+jl)) || facets(i,5+ku)-1 <= BH(facets(i, 5+il)-1, facets(i, 5+jl))) && facets(i,5+kl) >= OH(facets(i, 5+il)-1, facets(i, 5+jl))+1 % If facet is lower than block in -x direction
                                    facets(i,5) = 1;
                                    intern(j,1) = facets(i, 3);
                                    intern(j,2) = IM(facets(i, 5 + il),facets(i, 5 + jl));
                                    j = j + 1;
                                    
                                end
                            end
                        end

                    case 3 % east
                        if facets(i, 5 + iu) <= imax %not at the domain edge
                            if M(facets(i, 5 +  iu), facets(i, 5 + jl)) % There is a block in the +x direction
                                % If block in +x direction is higher
                                if (BH(facets(i, 5+iu), facets(i, 5+jl)) >= BH(facets(i, 5+iu)-1, facets(i, 5+jl)) || facets(i,5+ku)-1 <= BH(facets(i, 5+iu), facets(i, 5+jl))) && facets(i,5+kl) >= OH(facets(i, 5+iu), facets(i, 5+jl))+1
                                %if obj.topo(facets(i, 5 + jl), facets(i, 5 +  iu)) >= obj.topo(facets(i, 5 + jl), facets(i, 5 + iu)-1)
                                    facets(i,5) = 1;
                                    intern(j,1) = facets(i, 3);
                                    intern(j,2) = IM(facets(i, 5 + iu)-1, facets(i, 5 + jl));
                                    j = j + 1;
                                    
                                end
                            end
                        end
                    case 4 % north
                        if facets(i, 5 + ju) <= jtot %not at the domain edge
                            if M(facets(i, 5 + il), facets(i, 5 + ju)) % There is a block in the +y direction
                                % If block in +y direction is higher than current block
                                if (BH(facets(i, 5+il), facets(i, 5+ju)) >= BH(facets(i, 5+il), facets(i, 5+ju)-1) || facets(i,5+ku)-1 <= BH(facets(i, 5+il), facets(i, 5+ju))) && facets(i,5+kl) >= OH(facets(i, 5+il), facets(i, 5+ju))+1
                                %if obj.topo(facets(i, 5 + ju), facets(i, 5 + il)) >= obj.topo(facets(i, 5 + ju)-1, facets(i, 5 + il))
                                    facets(i,5) = 1;
                                    intern(j,1) = facets(i, 3);
                                    intern(j,2) = IM(facets(i, 5+il), facets(i, 5+ju)-1);
                                    j = j + 1;

                                end
                            end
                        end
                    case 5 % south
                        if facets(i, 5 + jl) - 1 >= 1 %not at the domain edge
                            if M(facets(i, 5 + il), facets(i, 5 + jl) - 1) % There is a block in the -y direction
                                % If block in -y direction is higher
                                if (BH(facets(i, 5+il), facets(i, 5+jl)-1) >= BH(facets(i, 5+il), facets(i, 5+jl)) || facets(i,5+ku)-1 <= BH(facets(i, 5+il), facets(i, 5+jl)-1)) && facets(i,5+kl) >= OH(facets(i, 5+il), facets(i, 5+jl)-1)+1
                                %if obj.topo( facets(i, 5 + jl) - 1, facets(i, 5 + il)) >= obj.topo(facets(i, 5 + jl), facets(i, 5 + il))
                                    facets(i, 5) = 1;
                                    intern(j, 1) = facets(i, 3);
                                    intern(j, 2) = IM(facets(i, 5+il), facets(i, 5+jl));
                                    j = j + 1;
                                end                           
                            end
                        end
                    case 6 % bottom
                         if all(OH(facets(i,5+il):facets(i,5+iu)-1, facets(i,5+jl):facets(i,5+ju)-1))
                             if facets(i,5+kl)-1 ~= OH(facets(i,5+il),facets(i,5+jl)) % If not at the overhang
%                              ids = find(ismember(facets(:,[1,6:9]), facets(i,[1,6:9]),'rows'));
%                              if length(ids) ~= 1 % there are bottom facets immediately below or above
%                                  kcurrent = facets(i,5+kl);
%                                  ks = facets(ids, 5+kl);
%                                  if kcurrent ~= min(ks) % There are lower bottom facets in the stack, so remove
                                      facets(i,5) = 1;
                                      intern(j,1) = facets(i, 3);
                                      intern(j,2) = IM(facets(i, 5 + il), facets(i, 5 + jl));
                                      j = j + 1;
%                                  end
                              end
                         else % Remove if no overhang
                             %if (facets(i,5+ku) ~= 1) % Don't remove ground level facets (for generating closed stl)
                             facets(i,5) = 1;
                             intern(j,1) = facets(i, 3);
                             intern(j,2) = IM(facets(i, 5 + il), facets(i, 5 + jl));
                             j = j + 1;
                             %end                               
                         end
                        
                end
            end

            %disp('No. facets removed:')
            nintern = sum(facets(:,5));
            intern = intern(1:nintern,:);
               
            bblk0 = facets(:, 3);
            bblk1 = bblk0;

            for n=1:nintern
                %list=find(bblk0==intern(n,2)); % SO: Think IM causes problems with blocks on top of each other so just use block id
                list=find(bblk0==intern(n,1));
                list2=find(bblk1==bblk1(list(1))); % Find where facet block ids match the internal facet block id
                bblk1(list2)=intern(n,1);
            end

            % assign building id to facets
            bblku = unique(bblk1);
            nbld = length(bblku);
            nbuildings = length(bblku);
            for n = 1:nbld
                facets(bblk1 == bblku(n), 4) = n;
            end
            

            %% remove downward facets and internal facets
            % facets format: orientation, walltype, blockid, buildingid, isinternal
            %              il, iu, jl, ju, kl, ku
            
            %for n=1:length(obj.overhangfacetsids_above)
                
            
            %if obj.loverhangs
                %notoverhangfacetids = setdiff(1:nfcts, overhangfacetids);
                %sel = intersect(find(facets(:,5) ~= 1), notoverhangfacetids);
            %else
                %sel = find(facets(:,1) ~= bot & facets(:,5) ~= 1)
                sel = find(facets(:,5) ~= 1); 
            %end
            
            % Also remove non-top facets of overhang blocks, and keep bottom facets of overlying block.
            
            facets = facets(sel, :);

            %% remove facets at domain edge

            sel = find(~((facets(:,6) == 1 & facets(:,7) == 1) | (facets(:,6)== imax + 1 & facets(:,7) == imax + 1) | (facets(:,8) == 1 & facets(:,9) == 1) | (facets(:,8) == jtot + 1 & facets(:,9) == jtot + 1)));
            facets = facets(sel, :);
            %if obj.lEB && (length(sel) ~= nfcts)
                  %error("Can't have blocks on edge of domain when using energy balance")
            %end

%             overhangfacets = [];
%             for i=1:length(obj.overhangblockids)
%                 nbuildings = nbuildings+1;
%                 id = obj.overhangblockids(i);
%                 %overhangfacets(i,:) = [1,1,id,0, 0, blocks(id, 1), blocks(id, 2)+1, blocks(id,3),  blocks(id,4)+1, blocks(id,5), blocks(id,5)];
%                 overhangfacets(i,:) = [1,1,id,nbuildings, 0, blocks(id, 1), blocks(id, 2)+1, blocks(id,3),  blocks(id,4)+1, blocks(id,5), blocks(id,5)];
%             end
%             facets = [facets; overhangfacets];
            
            obj.nfcts = size(facets, 1);
            preprocessing.addvar(obj, 'facets', facets)
            preprocessing.addvar(obj, 'nblockfcts', size(facets, 1))
            preprocessing.addvar(obj, 'nbuildings', nbuildings);
        end
        
        function block2fac_overlayer(obj, layer)
            layerheight = obj.layerheights(layer);
            overlayer = squeeze(obj.overlayers(:,:,layer));
            if layer==1 
                prevlayer = obj.topo;
            else
                prevlayer = obj.topo + squeeze(obj.overlayers(:,:,layer-1));
            end          
            
            imax = obj.imax;
            jtot = obj.jtot;
            
            if layer==1
                blockslayer = obj.overblocks(1:obj.nblockslayers(1),:);
            else
                blockslayer = obj.overblocks(sum(obj.nblockslayers(1:layer-1))+1:sum(obj.nblockslayers(1:layer)),:);
            end
            
            nblockslayer = obj.nblockslayers(layer);
            blockslayer(:, [5,6]) = blockslayer(:,[5,6]) + 1; % indices for k start at zero?
            %obj.blocks = blocks;
            %% create Mask-matrix
            % this new mask is in x,y coordinates, not y,x coordinates as before
            M = zeros(imax, jtot);
            IM = zeros(imax , jtot);
            BH = zeros(imax, jtot);
            for i = 1:nblockslayer
                xl = blockslayer(i,1);
                xu = blockslayer(i,2);
                yl = blockslayer(i,3);
                yu = blockslayer(i,4);
                M(xl:xu, yl:yu) = 1;
                IM(xl:xu, yl:yu) = i; % What about when blocks are on top of each other?
                BH(xl:xu, yl:yu) = max(BH(xl:xu, yl:yu), blockslayer(i,6));
            end
            %preprocessing.addvar(obj, 'BH', BH)
            
            top = 1; west = 2; east = 3; north = 4; south=5; bot = 6;

            il = 1; iu = 2;
            jl = 3; ju = 4;
            kl = 5; ku = 6;

            nfcts = nblockslayer * 6;
            facets = int32(zeros(nfcts, 10));

            for j = 1:nblockslayer
                for k = top:bot
                    i = (j - 1) * 6 + k;
                    facets(i, 1) = k; % orientation
                    % for all orientations apart from bottom, wall type =
                    % blk(j, 6+5) = blk(j,11) = 1. For bottom, wall type =
                    % blk(j, 6+6) = blk(j,12) = 0. This depends on what
                    % bllk was defined as before, and needs to be changed
                    % so that it is easier to change wall type.
                    facets(i,2) = blockslayer(j, 6 + min(k,5)); % wall type
                    facets(i,3) = j; % blockid
                    switch(k)
                        case top
                            facets(i, 5 + il) = blockslayer(j, il);
                            facets(i, 5 + jl) = blockslayer(j, jl);
                            facets(i, 5 + kl) = blockslayer(j, ku) + 1;
                            facets(i, 5 + iu) = blockslayer(j, iu) + 1;
                            facets(i, 5 + ju) = blockslayer(j, ju) + 1;
                            facets(i, 5 + ku) = blockslayer(j, ku) + 1;
                        case west
                            facets(i, 5 + il) = blockslayer(j, il);
                            facets(i, 5 + jl) = blockslayer(j, jl);
                            facets(i, 5 + kl) = blockslayer(j, kl);
                            facets(i, 5 + iu) = blockslayer(j, il);
                            facets(i, 5 + ju) = blockslayer(j, ju) + 1;
                            facets(i, 5 + ku) = blockslayer(j, ku) + 1;
                        case east
                            facets(i, 5 + il) = blockslayer(j, iu) + 1;
                            facets(i, 5 + jl) = blockslayer(j, jl);
                            facets(i, 5 + kl) = blockslayer(j, kl);
                            facets(i, 5 + iu) = blockslayer(j, iu) + 1;
                            facets(i, 5 + ju) = blockslayer(j, ju) + 1;
                            facets(i, 5 + ku) = blockslayer(j, ku) + 1;
                        case north
                            facets(i, 5 + il) = blockslayer(j, il);
                            facets(i, 5 + jl) = blockslayer(j, ju) + 1;
                            facets(i, 5 + kl) = blockslayer(j, kl);
                            facets(i, 5 + iu) = blockslayer(j, iu) + 1;
                            facets(i, 5 + ju) = blockslayer(j, ju) + 1;
                            facets(i, 5 + ku) = blockslayer(j, ku) + 1;
                        case south
                            facets(i, 5 + il) = blockslayer(j, il);
                            facets(i, 5 + jl) = blockslayer(j, jl);
                            facets(i, 5 + kl) = blockslayer(j, kl);
                            facets(i, 5 + iu) = blockslayer(j, iu) + 1;
                            facets(i, 5 + ju) = blockslayer(j, jl);
                            facets(i, 5 + ku) = blockslayer(j, ku) + 1;
                        case bot
                            facets(i, 5 + il) = blockslayer(j, il);
                            facets(i, 5 + jl) = blockslayer(j, jl);
                            facets(i, 5 + kl) = blockslayer(j, kl);
                            facets(i, 5 + iu) = blockslayer(j, iu) + 1;
                            facets(i, 5 + ju) = blockslayer(j, ju) + 1;
                            facets(i, 5 + ku) = blockslayer(j, kl);
                    end
                end
            end          
            
            % First remove the intersection between the layer below and above.
            rem =[];
            for i=1:obj.nfcts
                if obj.facets(i,1) == 1
                    %find(ismember(facets(:,6:11), obj.facets(i,6:11),'rows'))
%                     if ~isempty(find(ismember(facets(:,6:11), obj.facets(i,6:11),'rows')))
%                         rem = [rem; i];
%                     end
                     % If facet height is equal to layer height and there are overlying blocks at the facet location
                     if (obj.facets(i,5+kl)-1 == layerheight) && all(overlayer(obj.facets(i,5+jl):obj.facets(i,5+ju)-1, obj.facets(i,5+il):obj.facets(i,5+iu)-1), 'all')
                        rem = [rem; i];
                     end
                end
            end
            
            obj.facets(rem,:) = [];
            obj.nfcts = obj.nfcts - length(rem);
            obj.nblockfcts = obj.nblockfcts - length(rem);
            
            intern = zeros(nfcts, 2);
            j = 1;
            
%             for i = 1:nfcts
%                 switch(facets(i, 1))
%                     %only one test necessary, since the whole edge is internal, or not
%                     %remember that facets stores facet coordinates, not block coordinates
%                     case 1
%                         ids = find(ismember(facets(:,[1,6:9]), facets(i,[1,6:9]),'rows')); 
%                         if length(ids) ~= 1 % there are top facets immediately below or above.
%                               kcurrent = facets(i,5+kl);
%                               ks = facets(ids, 5+kl);
%                               if kcurrent == min(ks) % If it is the lowest in the stack
%                                   facets(i,5) = 1;
%                                   intern(j,1) = facets(i, 3);
%                                   intern(j,2) = IM(facets(i, 5 + il), facets(i, 5 + jl));
%                                   j = j + 1;
%                               end
%                         end
%                     case 2 % west
%                         if facets(i, 5+il)-1 >= 1 %not at the domain edge
%                             %if M(facets(i, 5 +il)-1 , facets(i, 5+jl)) % if there is a block in the -y direction
%                                 
%                                ids = find(ismember(facets(:,6:11), facets(i,6:11),'rows'))'; % Find facets with exactly the same coordinates
%                                if length(ids) > 1
%                                    for id = ids%(2:end) % remove duplicates
%                                        facets(id,5) = 1;
%                                        intern(j,1) = facets(id, 3);
%                                        %intern(j,2) = IM(facets(id, 5+il)-1 ,facets(id, 5+jl));
%                                        intern(j,2) = IM(facets(id, 5+il) ,facets(id, 5+jl));
%                                        j = j + 1;
%                                    end
%                                end
%                                
% %                                 if ismember(i, overhangfacetids_above)
% %                                     if M(facets(id, 5+il)-1-1, facets(id, 5+jl)) % If there is a block adjacent in -x direction (need to be careful if there is terrain)                                  
% %                                         facets(i,5) = 1;
% %                                         intern(j,1) = facets(i, 3);
% %                                         intern(j,2) = IM(facets(i, 5+il)-1, facets(i, 5+jl));
% %                                         j = j + 1;
% %                                     end
% %                                 end
%                                
%                             %end
%                         end
%                     case 3 % east
%                         if facets(i, 5 + iu) <= imax %not at the domain edge
%                             %if M(facets(i, 5+iu), facets(i, 5+jl)) % if building
% 
%                                 ids = find(ismember(facets(:,6:11), facets(i,6:11),'rows'))';
% 
%                                 if length(ids) > 1
%                                     %i
%                                     %disp('Duplicate facet')
%                                     for id = ids%(2:end)
%                                         facets(id,5) = 1;
%                                         intern(j,1) = facets(id, 3);
%                                         %intern(j,2) = IM(facets(id, 5+iu), facets(id, 5+jl));
%                                         intern(j,2) = IM(facets(id, 5+iu)-1, facets(id, 5+jl));
%                                         j = j + 1;
%                                     end
%                                 end
%                                 
% %                                 if ismember(i, overhangfacetids_above)
% %                                     if M(facets(id, 5+iu)+1, facets(id, 5+jl)) % If there is a block adjacent in +x direction (need to be careful if there is terrain)                                  
% %                                         facets(i,5) = 1;
% %                                         intern(j,1) = facets(i, 3);
% %                                         intern(j,2) = IM(facets(i, 5+iu), facets(i, 5+jl));
% %                                         j = j + 1;
% %                                     end
% %                                 end
%                                 
%                             %end
%                         end
%                     case 4 % north
%                         if facets(i, 5 + ju) <= jtot %not at the domain edge
%                             %if M(facets(i, 5+il), facets(i, 5+ju))
%                                 ids = find(ismember(facets(:,6:11), facets(i,6:11),'rows'))';
%                                 if length(ids) > 1
%                                     %i
%                                     %disp('Duplicate facet')
%                                     for id = ids%(2:end)
%                                         facets(id,5) = 1;
%                                         intern(j,1) = facets(id, 3);
%                                         %intern(j,2) = IM(facets(id, 5+il), facets(id, 5+ju));
%                                         intern(j,2) = IM(facets(id, 5+il), facets(id, 5+ju)-1);
%                                         j = j + 1;
%                                     end
%                                 end
%                                  
% %                                 if ismember(i, overhangfacetids_above)
% %                                     if M(facets(i, 5+il), facets(i, 5+ju)+1) % If there is a block adjacent in +y direction (need to be careful if there is terrain)                                  
% %                                         facets(i,5) = 1;
% %                                         intern(j,1) = facets(i, 3);
% %                                         intern(j,2) = IM(facets(i, 5+il), facets(i, 5+ju));
% %                                         j = j + 1;
% %                                     end
% %                                 end
%                                 
%                             %end
%                         end
%                     case 5 % south
%                         if facets(i, 5 + jl) - 1 >= 1 %not at the domain edge
%                             %if M(facets(i, 5+il), facets(i, 5+jl)-1)
%                                 ids = find(ismember(facets(:,6:11), facets(i,6:11),'rows'))';
%                                 if length(ids) > 1
%                                     %i
%                                     %disp('Duplicate facet')
%                                     for id = ids%(2:end)
%                                         facets(id, 5) = 1;
%                                         intern(j, 1) = facets(id, 3);
%                                         %intern(j, 2) = IM(facets(id, 5+il), facets(id, 5+jl)-1);
%                                         intern(j, 2) = IM(facets(id, 5+il), facets(id, 5+jl));
%                                         j = j + 1;
%                                     end
%                                 end
%                                 
% %                                 if ismember(i, overhangfacetids_above)
% %                                     if M(facets(id, 5+il), facets(id, 5+jl)-1-1) % If there is a block adjacent in -y direction (need to be careful if there is terrain)                                  
% %                                         facets(i,5) = 1;
% %                                         intern(j,1) = facets(i, 3);
% %                                         intern(j,2) = IM(facets(i, 5+il), facets(i, 5+jl)-1);
% %                                         j = j + 1;
% %                                     end
% %                                 end
%                                 
%                             %end
%                         end
%                      case 6
%                          % If not at layer height, remove
%                          if facets(i,5+kl)-1 ~= obj.layermax % -1 because blocks have been raised here
%                              facets(i, 5) = 1;
%                              intern(j, 1) = facets(i, 3);
%                              intern(j,2) = IM(facets(i, 5 + il), facets(i, 5 + jl));
%                              j = j + 1;
%                          else % If at layer height
%                              % If the facet aligns with the topomax
%                              if facets(i,5+kl)-1 == obj.topo(facets(i,5+jl), facets(i,5+il))
%                                  facets(i, 5) = 1;
%                                  intern(j, 1) = facets(i, 3);
%                                  intern(j,2) = IM(facets(i, 5 + il), facets(i, 5 + jl));
%                                  j = j + 1;
%                              end
%                          end
%                 end
%             end
                       
            for i = 1:nfcts
                switch(facets(i, 1))
                    %only one test necessary, since the whole edge is internal, or not
                    %remember that facets stores facet coordinates, not block coordinates
                    case 1
                        ids = find(ismember(facets(:,[1,6:9]), facets(i,[1,6:9]),'rows'));
                        if length(ids) ~= 1
                            kcurrent = facets(i,5+kl);
                            ks = facets(ids, 5+kl);
                            if kcurrent ~= max(ks) % There are lower top facets in the stack, so remove
                                facets(i,5) = 1;
                                intern(j,1) = facets(i, 3);
                                intern(j,2) = IM(facets(i, 5 + il), facets(i, 5 + jl));
                                j = j + 1;
                            end
                        end
                    case 2 % west
                        if facets(i, 5 + il) - 1 >= 1 %not at the domain edge
                            % If there is a block in the -x direction
                            if M(facets(i, 5 + il) - 1 , facets(i, 5 + jl)) % There is a block in the -x direction
                                % If block in -x direction is higher than current block or the current facet is lower than the next block, and the current facet is higher than any overhang.
                                if (BH(facets(i, 5+il)-1, facets(i, 5+jl)) >= BH(facets(i, 5+il), facets(i, 5+jl)) || facets(i,5+ku)-1 <= BH(facets(i, 5+il)-1, facets(i, 5+jl)))% && facets(i,5+kl) >= OH(facets(i, 5+il)-1, facets(i, 5+jl))+1 % If facet is lower than block in -x direction
                                    facets(i,5) = 1;
                                    intern(j,1) = facets(i, 3);
                                    intern(j,2) = IM(facets(i, 5 + il),facets(i, 5 + jl));
                                    j = j + 1;
                                    
                                end
                            end
                        end

                    case 3 % east
                        if facets(i, 5 + iu) <= imax %not at the domain edge
                            if M(facets(i, 5 +  iu), facets(i, 5 + jl)) % There is a block in the +x direction
                                % If block in +x direction is higher
                                if (BH(facets(i, 5+iu), facets(i, 5+jl)) >= BH(facets(i, 5+iu)-1, facets(i, 5+jl)) || facets(i,5+ku)-1 <= BH(facets(i, 5+iu), facets(i, 5+jl)))% && facets(i,5+kl) >= OH(facets(i, 5+iu), facets(i, 5+jl))+1
                                %if obj.topo(facets(i, 5 + jl), facets(i, 5 +  iu)) >= obj.topo(facets(i, 5 + jl), facets(i, 5 + iu)-1)
                                    facets(i,5) = 1;
                                    intern(j,1) = facets(i, 3);
                                    intern(j,2) = IM(facets(i, 5 + iu)-1, facets(i, 5 + jl));
                                    j = j + 1;
                                    
                                end
                            end
                        end
                    case 4 % north
                        if facets(i, 5 + ju) <= jtot %not at the domain edge
                            if M(facets(i, 5 + il), facets(i, 5 + ju)) % There is a block in the +y direction
                                % If block in +y direction is higher than current block
                                if (BH(facets(i, 5+il), facets(i, 5+ju)) >= BH(facets(i, 5+il), facets(i, 5+ju)-1) || facets(i,5+ku)-1 <= BH(facets(i, 5+il), facets(i, 5+ju)))% && facets(i,5+kl) >= OH(facets(i, 5+il), facets(i, 5+ju))+1
                                %if obj.topo(facets(i, 5 + ju), facets(i, 5 + il)) >= obj.topo(facets(i, 5 + ju)-1, facets(i, 5 + il))
                                    facets(i,5) = 1;
                                    intern(j,1) = facets(i, 3);
                                    intern(j,2) = IM(facets(i, 5+il), facets(i, 5+ju)-1);
                                    j = j + 1;

                                end
                            end
                        end
                    case 5 % south
                        if facets(i, 5 + jl) - 1 >= 1 %not at the domain edge
                            if M(facets(i, 5 + il), facets(i, 5 + jl) - 1) % There is a block in the -y direction
                                % If block in -y direction is higher
                                if (BH(facets(i, 5+il), facets(i, 5+jl)-1) >= BH(facets(i, 5+il), facets(i, 5+jl)) || facets(i,5+ku)-1 <= BH(facets(i, 5+il), facets(i, 5+jl)-1))% && facets(i,5+kl) >= OH(facets(i, 5+il), facets(i, 5+jl)-1)+1
                                %if obj.topo( facets(i, 5 + jl) - 1, facets(i, 5 + il)) >= obj.topo(facets(i, 5 + jl), facets(i, 5 + il))
                                    facets(i, 5) = 1;
                                    intern(j, 1) = facets(i, 3);
                                    intern(j, 2) = IM(facets(i, 5+il), facets(i, 5+jl));
                                    j = j + 1;
                                end                           
                            end
                        end
                    case 6 % bottom
                         % If not at layer height, remove
                         if facets(i,5+kl)-1 ~= layerheight % -1 because blocks have been raised here
                             facets(i, 5) = 1;
                             intern(j, 1) = facets(i, 3);
                             intern(j,2) = IM(facets(i, 5 + il), facets(i, 5 + jl));
                             j = j + 1;
                         else % If at layer height
                             % If the facet aligns perfectly with the topomax
                             %if facets(i,5+kl)-1 == obj.topo(facets(i,5+jl), facets(i,5+il))
                             if facets(i,5+kl)-1 == prevlayer(facets(i,5+jl), facets(i,5+il))
                                 facets(i, 5) = 1;
                                 intern(j, 1) = facets(i, 3);
                                 intern(j,2) = IM(facets(i, 5 + il), facets(i, 5 + jl));
                                 j = j + 1;
                             end
                         end
                        
                end
            end
            
            disp('No. facets removed:')
            j
            nintern = sum(facets(:,5));
            intern = intern(1:nintern,:);

            bblk0 = facets(:, 3);
            bblk1 = bblk0;

            for n=1:nintern
                list=find(bblk0==intern(n,2));
                list2=find(bblk1==bblk1(list(1)));
                bblk1(list2)=intern(n,1);
            end

            % assign building id to facets
            bblku = unique(bblk1);
            nbld = length(bblku);
            nbuildings = length(bblku);
            for n = 1:nbld
                facets(bblk1 == bblku(n), 4) = n;
            end
            

            %% remove downward facets and internal facets
            % facets format: orientation, walltype, blockid, buildingid, isinternal
            %              il, iu, jl, ju, kl, ku
            
            %for n=1:length(obj.overhangfacetsids_above)
                
            
            %if obj.loverhangs
                %notoverhangfacetids = setdiff(1:nfcts, overhangfacetids);
                %sel = intersect(find(facets(:,5) ~= 1), notoverhangfacetids);
            %else
                %sel = find(facets(:,1) ~= bot & facets(:,5) ~= 1) 
            %end
            sel = find(facets(:,5) ~= 1);
            
            % Also remove non-top facets of overhang blocks, and keep bottom facets of overlying block.
            
            facets = facets(sel, :);

            %% remove facets at domain edge

            sel = find(~((facets(:,6) == 1 & facets(:,7) == 1) | (facets(:,6)== imax + 1 & facets(:,7) == imax + 1) | (facets(:,8) == 1 & facets(:,9) == 1) | (facets(:,8) == jtot + 1 & facets(:,9) == jtot + 1)));
            facets = facets(sel, :);
            %if obj.lEB && (length(sel) ~= nfcts)
                  %error("Can't have blocks on edge of domain when using energy balance")
            %end
            

                    
            
             
            if layer == 1
                facets(:,3) = facets(:,3) + obj.nunderblocks;
                facets(:,4) = facets(:,4) + obj.nbuildings;
            else
                facets(:,3) = facets(:,3) + obj.nunderblocks+sum(obj.nblockslayers(1:layer-1));
                facets(:,4) = facets(:,4) + obj.nbuildings;
            end
            
            obj.overlayerfacets = [obj.overlayerfacets; facets];
            obj.noverlayerfacets = obj.noverlayerfacets + size(facets,1);
            %preprocessing.addvar(obj, 'overlayerfacets', facets);
            %preprocessing.addvar(obj, 'noverlayerfacets', size(facets,1));
            obj.facets = [obj.facets; facets];
            obj.nfcts = obj.nfcts + size(facets,1);
            obj.nblockfcts = obj.nblockfcts + size(facets,1);
            %preprocessing.addvar(obj, 'nbuildings', nbuildings);
            obj.nbuildings = obj.nbuildings+nbuildings;
        end

        function addboundingwalls(obj)
            % Add a bounding wall around the domain used in radiation calculations only
            % THIS DOES NOT YET DEAL WITH PERIODIC GEOMETRY (I.E. CANYONS)
            imax = obj.imax;
            jtot = obj.jtot;
            maxsize = obj.maxsize;
            blocks = obj.blocks;
            height = floor(median(blocks(:, 6)));

            nxwalls = ceil(jtot / maxsize);
            remx = rem(jtot, maxsize);
            nywalls = ceil(imax / maxsize);
            remy = rem(imax, maxsize);
            nzw = ceil((height + 1) / maxsize);
            remz = rem((height + 1), maxsize);

            nboundingwallfacets = 2 * nzw * (nxwalls + nywalls);
            boundingwallfacets = zeros(nboundingwallfacets, 11);

            if remx > 0
                for j = 1:nzw
                    if ((j == nzw) && (remz > 0))
                        lh = height - remz + 1;
                        uh = height;
                    else
                        lh = (j - 1) * maxsize;
                        uh = j * maxsize - 1;
                    end
                    for i = 1:(nxwalls-1)
                        boundingwallfacets((i-1) * nzw + j, 6:11) = [1, 1, (i - 1) * maxsize + 1, i * maxsize, lh, uh];
                        boundingwallfacets((i-1) * nzw + j + nxwalls * nzw, 6:11) =  [imax, imax, (i - 1) * maxsize + 1, i * maxsize, lh, uh];
                    end
                    boundingwallfacets((nxwalls - 1) * nzw + j, 6:11) = [1, 1, jtot - remx + 1, jtot, lh, uh];
                    boundingwallfacets((nxwalls - 1) * nzw + j + nxwalls * nzw, 6:11) = [imax, imax, jtot - remx + 1, jtot, lh, uh];
                end
            else
                disp('remx <= 0')
                for j = 1:nzw
                    if ((j == nzw) && (remz > 0))
                        lh = height - remz + 1;
                        uh = height;
                    else
                        %lh = (j - 1) * nzw;
                        lh = (j - 1) * maxsize;
                        %uh = j * nzw - 1;
                        uh = j * maxsize - 1;
                    end
                    for i = 1:nxwalls
                        boundingwallfacets((i - 1) * nzw + j, 6:11) = [1, 1, (i - 1) * maxsize + 1, i * maxsize, lh, uh];
                        boundingwallfacets((i - 1) * nzw + j + nxwalls * nzw, 6:11) = [imax, imax, (i - 1) * maxsize + 1, i * maxsize, lh, uh];
                    end
                end
            end

            if remy > 0
                for j = 1:nzw
                    if ((j == nzw) && (remz > 0))
                        lh = height - remz + 1;
                        uh = height;
                    else
                        lh = (j - 1) * maxsize;
                        uh = j * maxsize - 1;
                    end
                    for i = 1:(nywalls - 1)
                        boundingwallfacets(2 * nzw * nxwalls + (i - 1) * nzw + j, 6:11) = [(i - 1) * maxsize + 1, i * maxsize, jtot, jtot, lh, uh];
                        boundingwallfacets(2 * nzw * nxwalls + (i - 1) * nzw + j + nywalls * nzw, 6:11) = [(i - 1) * maxsize + 1, i * maxsize, 1, 1, lh, uh];
                    end
                    boundingwallfacets(2 * nzw * nxwalls + (nywalls - 1) * nzw + j, 6:11) = [imax - remy + 1, imax, jtot, jtot, lh, uh];
                    boundingwallfacets(2 * nzw * nxwalls + (nywalls - 1) * nzw + j + nywalls * nzw, 6:11) = [imax - remy + 1, imax, 1, 1, lh, uh];
                end
            else
                disp('remy <= 0')
                for j = 1:nzw
                    if ((j == nzw) && (remz > 0))
                        lh = height - remz + 1;
                        uh = height;
                    else
                        %lh = (j - 1) * nzw;
                        lh = (j - 1) * maxsize;
                        %uh = j * nzw - 1;
                        uh = j * maxsize - 1;
                    end
                    for i = 1:nywalls
                        boundingwallfacets(2 * nzw * nxwalls + (i - 1) * nzw + j, 6:11) = [(i - 1) * maxsize + 1, i * maxsize, jtot, jtot, lh, uh];
                        boundingwallfacets(2 * nzw * nxwalls + (i - 1) * nzw + j + nywalls * nzw, 6:11) = [(i - 1) * maxsize + 1, i * maxsize, 1, 1, lh, uh];
                    end
                end
            end


            for i = 1:(nxwalls * nzw)
                % East
                boundingwallfacets(i, 1:4) =  [3, -101, i, -101]; % west, facing east.
                % Add +1 and +1+1 - why is this done?
%                 boundingwallfacets(i, 9) = boundingwallfacets(i, 9) + 1;
%                 boundingwallfacets(i, 10) = boundingwallfacets(i, 10) + 1;
%                 boundingwallfacets(i, 11) = boundingwallfacets(i, 11) + 2;

                % West
                boundingwallfacets(nxwalls * nzw + i, 1:4) = [2, -101, i + nxwalls * nzw, -101]; % east, facing west.
                %Add +1 and +1+1 - why is this done?
%                 boundingwallfacets(nxwalls * nzw + i, 6) = boundingwallfacets(nxwalls * nzw + i, 6) + 1;
%                 boundingwallfacets(nxwalls * nzw + i, 7) = boundingwallfacets(nxwalls * nzw + i, 7) + 1;
%                 boundingwallfacets(nxwalls * nzw + i, 9) = boundingwallfacets(nxwalls * nzw + i, 9) + 1;
%                 boundingwallfacets(nxwalls * nzw + i, 10) = boundingwallfacets(nxwalls * nzw + i, 10) + 1;
%                 boundingwallfacets(nxwalls * nzw + i, 11) = boundingwallfacets(nxwalls * nzw + i, 11) + 2;
            end
            for i = 1:(nywalls * nzw)
                % South
                boundingwallfacets(2 * (nxwalls * nzw) + i, 1:4) = [5, -101, 2 * nxwalls * nzw + i, -101]; %north, facing south.
                %Add +1 and +1+1 - why is this done?
%                 boundingwallfacets(2 * (nxwalls * nzw) + i, 7) = boundingwallfacets(2 * (nxwalls * nzw) + i, 7) + 1;
%                 boundingwallfacets(2 * (nxwalls * nzw) + i, 8) = boundingwallfacets(2 * (nxwalls * nzw) + i, 8) + 1;
%                 boundingwallfacets(2 * (nxwalls * nzw) + i, 9) = boundingwallfacets(2 * (nxwalls * nzw) + i, 9) + 1;
%                 boundingwallfacets(2 * (nxwalls * nzw) + i, 10) = boundingwallfacets(2 * (nxwalls * nzw) + i, 10) + 1;
%                 boundingwallfacets(2 * (nxwalls * nzw) + i, 11) = boundingwallfacets(2 * (nxwalls * nzw) + i, 11) + 2;

                % North
                boundingwallfacets(2 * (nxwalls * nzw) + (nywalls * nzw) + i, 1:4) = [4, -101, 2 * nxwalls * nzw + i + nywalls * nzw, -101]; %south, facing north.
                %Add +1 and +1+1 - why is this done?
%                 boundingwallfacets(2 * (nxwalls * nzw) + (nywalls * nzw) + i, 7) = boundingwallfacets(2 * (nxwalls * nzw) + (nywalls * nzw) + i, 7) + 1;
%                 boundingwallfacets(2 * (nxwalls * nzw) + (nywalls * nzw) + i, 10) = boundingwallfacets(2 * (nxwalls * nzw) + (nywalls * nzw) + i, 10) + 1;
%                 boundingwallfacets(2 * (nxwalls * nzw) + (nywalls * nzw) + i, 11) = boundingwallfacets(2 * (nxwalls * nzw) + (nywalls * nzw) + i, 11) + 2;
            end
            obj.facets(end + 1:end + nboundingwallfacets, :) = boundingwallfacets;
            preprocessing.addvar(obj, 'boundingwallfacets', boundingwallfacets);
            preprocessing.addvar(obj, 'nboundingwallfacets', nboundingwallfacets);
        end

        function createfloors(obj)
            % Create floors
            % fill space between blocks with floors (streets etc.)
            % floors only have x and y coordinates, with building ID = -101
            % obj.floorfacets holds floor facets information.
            % Format: orientation, wall type, block id, building id,
            % isinternal, il, iu, jl, ju, kl, ku
            imax = obj.imax;
            jtot = obj.jtot;
            maxsize = obj.maxsize;
            blocks = obj.blocks;
            nblocks = obj.nblocks;
            M = ones(imax, jtot);
            BI = zeros(imax, jtot); %block index mask
            corm = zeros(imax, jtot); %mask with all wall-floor corners
            cornm = zeros(imax, jtot);
            
            if obj.lfloors
                if maxsize ~= inf
                    nfloorfacets_x = ceil(imax / maxsize);
                    nfloorfacets_y = ceil(jtot / maxsize);
                    nfloorfacets = nfloorfacets_x * nfloorfacets_y;
                    floorfacets = zeros(nfloorfacets, 11);
                    for i = 1:nfloorfacets_x - 1
                        xl = (i - 1) * maxsize + 1;
                        xu = i * maxsize;
                        for j = 1:nfloorfacets_y - 1
                            yl = (j - 1) * maxsize + 1;
                            yu = j * maxsize;
                            floorfacets((i - 1) * nfloorfacets_y + j, :) = [1, -1, i, -1, 0, xl, xu, yl, yu, 0, 0];
                            blocks((i - 1) * nfloorfacets_y + j, :) = [xl, xu, yl, yu, 0, 0, (i - 1) * nfloorfacets_y + j, 0, 0, 0];
                        end
                    end
                    i = nfloorfacets_x;
                    rem_x = rem(imax, maxsize);
                    xl = imax - rem_x + 1;
                    xu = imax;
                    for j = 1:nfloorfacets_y - 1
                        yl = (j - 1) * maxsize + 1;
                        yu = j * maxsize;
                        floorfacets((i - 1) * nfloorfacets_y + j, :) = [1, -1, i, -1, 0, xl, xu, yl, yu, 0, 0];
                        blocks((i - 1) * nfloorfacets_y + j, :) = [xl, xu, yl, yu, 0, 0, (i - 1) * nfloorfacets_y + j, 0, 0, 0];
                    end
                    j = nfloorfacets_y;
                    rem_y = rem(jtot, maxsize);
                    yl = jtot - rem_y + 1;
                    yu = jtot;
                    for i = 1:nfloorfacets_x - 1
                        xl = (i - 1) * maxsize + 1;
                        xu = i * maxsize;
                        floorfacets((i - 1) * nfloorfacets_y + j, :) = [1, -1, i, -1, 0, xl, xu, yl, yu, 0, 0];
                        blocks((i - 1) * nfloorfacets_y + j, :) = [xl, xu, yl, yu, 0, 0, (i - 1) * nfloorfacets_y + j, 0, 0, 0];
                    end
                    
                    i = nfloorfacets_x;
                    j = nfloorfacets_y;
                    xl = imax - rem_x + 1;
                    xu = imax;
                    yl = jtot - rem_y + 1;
                    yu = jtot;
                    floorfacets((i - 1) * nfloorfacets_y + j, :) = [1, -1, i, -1, 0, xl, xu, yl, yu, 0, 0];
                    blocks((i - 1) * nfloorfacets_y + j, :) = [xl, xu, yl, yu, 0, 0, (i - 1) * nfloorfacets_y + j, 0, 0, 0];
                    
                else
                    nfloorfacets = 1;
                    floorfacets = [1, -1, 1, -1, 0, 1, imax, 1, jtot, 0, 0];
                    blocks = [1, imax, 1, jtot, 0, 0, 1, 0, 0, 0, 0];
                end
            else
                BI = zeros(obj.imax, obj.jtot);
                BH = zeros(obj.imax, obj.jtot);
                for i = 1:nblocks-obj.noverblocks
                    xl = blocks(i,1);
                    xu = blocks(i,2);
                    yl = blocks(i,3);
                    yu = blocks(i,4);
                    zl = blocks(i,5)+1;
                    zu = blocks(i,6)+1;
                    M(xl:xu, yl:yu) = 0;
                    for x=xl:xu
                        for y=yl:yu
                            if BI(x,y) == 0
                                BI(x,y) = i;
                                BH(x,y) = blocks(BI(x,y),5)+1;
                            elseif zl < blocks(BI(x,y), 5)
                                BI(x,y) = i;
                                BH(x,y) = blocks(BI(x,y),5)+1;
                            end
                        end
                    end
                end
                NM = 1 - M;
                
                % make them around blocks first, with 1 blocksize wide
                %numbers indicate floors, "--" and "|" walls  (only 1 blocksize wide, just two
                %number used to make diagonal clear)
                %
                %   ----------------------
                % |
                % |    -------------------
                % |   |3311111111111111111
                % |   |3311111111111111111
                % |   |22
                % |   |22
                % |   |22
                %after detsub
                %   ----------------------
                % |
                % |    -------------------
                % |   |2111111111111111111
                % |   |2211111111111111111
                % |   |22
                % |   |22
                % |   |22
                
                maxblocks = sum(M(:)); %there cant be more blocks then number of grid cells
                %floors = NaN(maxblocks,4); %allocate a maximum size for floors, reduce size later. xy coordinates of corners, counterclockwise, starting nortwest
                floorfacets = [zeros(maxblocks, 5), NaN(maxblocks, 4), zeros(maxblocks, 2)];
                c = 0;
                M2 = M;
                iM = zeros(size(M)); %save indeces
                %for i = 1:nblocks-obj.noverblocks
                BIlist = unique(BI(BI ~= 0));
                preprocessing.addvar(obj, 'BIlist', BIlist)
                for i=BIlist'
                    xl = blocks(i, 1);
                    xu = blocks(i, 2);
                    yl = blocks(i, 3);
                    yu = blocks(i, 4);
                    
                    %west
                    if xl - 1 >= 1 %not at domain edge
                        if BI(xl - 1, yl) == 0 % left neighbour is a floor
                            c = c + 1;
                            M2(xl - 1, yl:yu) = 0; %set to 2 for later check if it is a corner
                            iM(xl - 1, yl:yu) = c;
                            floorfacets(c, 6:9) = [xl - 1, xl - 1, yl, yu];
                            if (yl - 1 >= 1) && (xu + 1 <= imax) && (yu + 1 <= jtot) && (xl - 1 >= 1)
                                if BI(xl - 1, yl - 1) > 0 %corner with a north wall
                                    cornm(xl - 1, yl) = 8;
                                    floorfacets(c, 6:9) = [xl-1, xl-1, yl+1, yu];
                                    floorfacets(c+1, 6:9) = [xl-1, xl-1, yl, yl];
                                    c = c+1;
                                elseif BI(xl - 1, yu + 1) > 0 %corner with a south wall
                                    cornm(xl - 1, yu) = 10;
                                    floorfacets(c, 6:9) = [xl-1, xl-1, yl, yu-1];
                                    floorfacets(c+1, 6:9) = [xl-1, xl-1, yu, yu];
                                    c = c+1;
                                end
                            end
                        end
                    end
                    %east
                    if xu + 1 <= imax %not at domain edge
                        if BI(xu + 1, yl) == 0
                            c = c + 1;
                            M2(xu + 1, yl:yu) = 0;
                            iM(xu + 1, yl:yu) = c;
                            floorfacets(c, 6:9) = [xu + 1, xu + 1, yl, yu];
                            if (yl - 1 >= 1) && (xu + 1 <= imax) && (yu + 1 <= jtot) && (xl - 1 >= 1)
                                if BI(xu + 1, yl - 1) > 0 %corner with a north wall
                                    cornm(xu + 1, yl) = 12;
                                    floorfacets(c, 6:9) = [xu+1, xu+1, yl+1, yu];
                                    floorfacets(c+1, 6:9) = [xu+1, xu+1, yl, yl];
                                    c = c+1;
                                elseif BI(xu+1, yu + 1) > 0  %corner with a south wall
                                    cornm(xu+1, yu) = 15;
                                    floorfacets(c, 6:9) = [xu+1, xu+1, yl, yu-1];
                                    floorfacets(c+1, 6:9) = [xu+1, xu+1, yu, yu];
                                    c = c+1;
                                end
                            end
                        end
                    end
                    %north
                    if yu + 1 <= jtot %not aat domain edge
                        if BI(xu, yu + 1) == 0
                            c = c + 1;
                            M2(xl:xu, yu + 1) = 0;
                            iM(xl:xu, yu + 1) = c;
                            floorfacets(c, 6:9) = [xl, xu, yu + 1, yu + 1];
                            if (yl - 1 >= 1) && (xu + 1 <= imax) && (yu + 1 <= jtot) && (xl - 1 >= 1)
                                if BI(xl - 1, yu + 1) > 0 %corner with an east wall
                                    cornm(xl, yu + 1) = 12;
                                    floorfacets(c, 6:9) = [xl+1, xu, yu+1, yu+1];
                                    floorfacets(c+1, 6:9) = [xl, xl, yu+1, yu+1];
                                    c = c + 1;
                                elseif BI(xu + 1, yu + 1) > 0 %corner with a west wall
                                    cornm(xu, yu + 1) = 8;
                                    floorfacets(c, 6:9) = [xl, xu-1, yu+1, yu+1];
                                    floorfacets(c+1, 6:9) = [xu, xu, yu+1, yu+1];
                                    c = c + 1;
                                end
                            end
                        end
                    end
                    %south
                    if yl - 1 >= 1 %not at domain edge
                        if BI(xu, yl - 1) == 0
                            c = c + 1;
                            M2(xl:xu, yl - 1) = 0;
                            iM(xl:xu, yl - 1) = c;
                            %floors(c, :) = [xl, xu, yl - 1, yl - 1];
                            floorfacets(c, 6:9) = [xl, xu, yl - 1, yl - 1];
                            if (yl - 1 >= 1) && (xu + 1 <= imax) && (yu + 1 <= jtot) && (xl - 1 >= 1)
                                if BI(xl - 1, yl - 1) > 0 %corner with an east wall
                                    cornm(xl , yl - 1) = 15;
                                    floorfacets(c, 6:9) = [xl+1, xu, yl-1, yl-1];
                                    floorfacets(c+1, 6:9) = [xl, xl, yl-1, yl-1];
                                    c = c + 1;
                                elseif BI(xu + 1, yl - 1) > 0 %corner with a west wall
                                    cornm(xu, yl - 1) = 10;
                                    floorfacets(c, 6:9) = [xl, xu-1, yl-1, yl-1];
                                    floorfacets(c+1, 6:9) = [xu, xu, yl-1, yl-1];
                                    c = c + 1;
                                end
                            end
                        end
                    end
                end
                
                corm(iM > 0) = 1;
                
                
                %% remove identical facets in corners (if it's a 1x1 facet in both cases)
                %truncate matrix
                lnan2 = find(isnan(floorfacets(:,6)));
                if ~isempty(lnan2)
                    floorfacets(lnan2(1):lnan2(end), :) = [];
                end
                
                floorfacets = unique(floorfacets,'rows','stable');
                nfloorfacets = size(floorfacets, 1);
                count = 1;
                while count <= nfloorfacets
                    i = count;
                    if sum(floorfacets(:, 6) <= floorfacets(i, 6) & floorfacets(:, 7) >= floorfacets(i, 7) & floorfacets(:, 8) <= floorfacets(i, 8) & floorfacets(:, 9) >= floorfacets(i, 9)) > 1
                        floorfacets(i, :) = []; %this floor is contained within another and can be removed
                        nfloorfacets = nfloorfacets - 1;
                    %elseif floorfacets(i,6) == floorfacets(i,6) 
                    else
                        count = count + 1;
                    end
                end
                
                c = size(floorfacets,1);
                %% Make floors
                % make them in 1D first (fixed x, along y)
                while any(M2(:) > 0)
                    for i = 1:imax
                        ls = find(M2(i, :) == 1);
                        if ~isempty(ls)
                            first = ls(1);
                            if length(ls) > 1
                                last = ls(find(diff(ls)~=1, 1));
                                if isempty(last) % there is a floor all the way to the domain edge
                                    % there could already be a floor that borders a block, so check whether there already is a floor
                                    if sum(M2(i,ls(end)+1:jtot))==0
                                        last = min(ls(end),first + maxsize - 1);
                                    else
                                        last = min(jtot,first + maxsize - 1); %create floor to the domain end or as long as the limit allows
                                    end
                                else
                                    last = min(last,first + maxsize - 1);
                                end
                            else
                                last = first;
                            end
                            c = c + 1;
                            floorfacets(c, :) = [NaN(1, 5), i, i, first, last, NaN(1, 2)];
                            M2(i, first:last) = 0;
                            iM(i, first:last) = c;
                        end
                    end
                end
                
                lnan2 = find(isnan(floorfacets(:,6)));
                if ~isempty(lnan2)
                    floorfacets(lnan2(1):lnan2(end), :) = [];
                end
                
                nslice = size(floorfacets, 1);
                dsize = 1;
                sizeold = nslice;
                while dsize>0
                    i = 1;
                    while 1
                        a = floorfacets(i, 7);
                        if corm(a, floorfacets(i, 8)) %his is a floor that belongs to a corner with a wall, don't merge with others
                            i = i + 1;
                            continue
                        end
                        bv = find(floorfacets(:, 6) == (a + 1)); %all floors with xl == this floors xu+1
                        b2 = bv(find((floorfacets(bv, 8) == floorfacets(i, 8)) & (floorfacets(bv, 9) == floorfacets(i, 9)))); %all of these floors witch also have the same y dimensions (should only be one)
                        if ~isempty(b2) && ~corm(floorfacets(b2, 6), floorfacets(b2, 8))
                            floorfacets(i, 7) = floorfacets(b2, 7);
                            floorfacets(b2, :) = [];
                        end
                        i = i + 1;
                        if i >= size(floorfacets, 1)
                            break
                        end
                    end
                    dsize = sizeold - size(floorfacets, 1);
                    sizeold = size(floorfacets, 1);
                end
                
                ls = 999;
                while ~isempty(ls)
                    nfloorfacets = size(floorfacets, 1);
                    ls = find(floorfacets(:, 7) - floorfacets(:, 6) > maxsize);
                    floorfacets = [floorfacets; NaN(length(ls), 11)];
                    for i = 1:length(ls)
                        ind = ls(i);
                        floorfacets(nfloorfacets + i, :) = floorfacets(ind, :);
                        floorfacets(ind, 7) = floorfacets(ind, 6) + maxsize - 1;
                        floorfacets(nfloorfacets + i, 6) = floorfacets(ind, 7) + 1;
                    end
                end
                
                % merge floors in y, where possible and as long as smaller than maxsize, don't merge triple corners
                change = true;
                while change
                    change = false;
                    for j = 1:size(floorfacets, 1)
                        il = floorfacets(j, 6);
                        iu = floorfacets(j, 7);
                        jl = floorfacets(j, 8);
                        ju = floorfacets(j, 9);
                        if sum(sum(cornm(il:iu, jl:ju))) == 0 %no triple corner somewhere on this floor facet, try to merge along y
                            flu = find(floorfacets(:, 6) == il & floorfacets(:, 7) == iu & floorfacets(:, 8) == ju + 1); %floor with same x dimension on ju+1
                            fll = find(floorfacets(:, 6) == il & floorfacets(:, 7) == iu & floorfacets(:, 9) == jl-1); %floor with same x dimension on jl-1
                            if ~isempty(flu)
                                ilu = floorfacets(flu, 6);
                                iuu = floorfacets(flu, 7);
                                jlu = floorfacets(flu, 8);
                                juu = floorfacets(flu, 9);
                                if sum(sum(cornm(ilu:iuu, jlu:juu))) == 0 && (floorfacets(flu, 9) - floorfacets(j, 8) + 1 < maxsize)
                                    floorfacets(j, 9) = floorfacets(flu, 9);
                                    floorfacets(flu, :) = [];
                                    change = true;
                                end
                            elseif ~isempty(fll)
                                ill = floorfacets(fll, 6);
                                iul = floorfacets(fll, 7);
                                jll = floorfacets(fll, 8);
                                jul = floorfacets(fll, 9);
                                if sum(sum(cornm(ill:iul, jll:jul))) == 0 && (floorfacets(j, 9) - floorfacets(fll, 8) + 1 < maxsize)
                                    floorfacets(j, 8) = floorfacets(fll, 8);
                                    floorfacets(fll, :) = [];
                                    change = true;
                                end
                            end
                            
                        end
                        if change
                            break
                        end
                    end
                end
                               
                preprocessing.addvar(obj, 'cornm', cornm);
                nfloorfacets = size(floorfacets, 1);
                
                for i = 1:nfloorfacets
                    %floorfacets(i, 1:5) = [1, -1, i, -1, 0]; % SO: block ID should be nblocks + i
                    floorfacets(i, 1:5) = [1, -1, nblocks + i, -1, 0];
                    floorfacets(i, 10:11) = [0, 0];
                end
                               
                %disp([num2str(obj.nfcts) ' facets, of which: ' num2str(obj.nblockfcts) ' from buildings, ' num2str(obj.nboundingwallfacets) ' from walls, ' num2str(obj.nfloorfacets) ' from floors.'])
                
                nblocks = obj.nblocks;
                blocks(:, 5:6) = blocks(:, 5:6) + 1;
                
                blocks(nblocks + 1 : nblocks + nfloorfacets, :) = zeros(nfloorfacets, 11);
                blocks(nblocks + 1 : nblocks + nfloorfacets, 1:6) = floorfacets(:, 6:11);
                blocks(:, 7:12) = zeros(nblocks + nfloorfacets, 6);
                
                
                facets = obj.facets;
                for i = 1:obj.nblockfcts
                    blocks(facets(i, 3), facets(i, 1) + 6) = i;
                end
                
                % add floors below buildings
                % do we need/want this? DALES will loop over more blocks but not really
                % do anything, on the other hand the statistics might look better?
                % If we use this, then the number of blocks in modibm should
                % be different from the number of blocks in masking matrices to avoid
                % looping. Also these blocks don't have corresponding facets and thus
                % access element 0 of any facet array in DALES
                
                j = nblocks + 1;
                for i = (obj.nblockfcts + obj.nboundingwallfacets + 1):(obj.nblockfcts + obj.nboundingwallfacets + nfloorfacets) %for floors
                    blocks(j, 7) = i;
                    j = j + 1;
                end
            
            end

            obj.nfcts = obj.nblockfcts + obj.nboundingwallfacets + nfloorfacets;
            obj.facets = [obj.facets; floorfacets];
            obj.blocks = blocks;
            preprocessing.addvar(obj, 'floorfacets', floorfacets)
            preprocessing.addvar(obj, 'nfloorfacets', nfloorfacets)
            preprocessing.addvar(obj, 'nblockstotal', obj.nblocks + nfloorfacets);
        end

        function vsolc(obj)
            xh = obj.xh;
            yh = obj.yh;
            zh = obj.zh;
            blocks = obj.blocks; nblocks = obj.nblocks;
            facets = obj.facets; nfcts  = obj.nfcts;
            floorfacets = obj.floorfacets; nfloorfacets = obj.nfloorfacets;
            boundingwallfacets = obj.boundingwallfacets; nboundingwallfacets = obj.nboundingwallfacets; nbw = obj.nboundingwallfacets;
            buildings = obj.buildings;
            centerweight = obj.centerweight; cornerweight = obj.cornerweight;

%             buildings(:, 7) = 0;
%             for i = 1:nblocks
%                 index = find(buildings(:, 1) <= buildings(i, 1) & buildings(:, 2) >= buildings(i, 2) & buildings(:, 3) <= buildings(i, 3) & buildings(:, 4) >= buildings(i, 4));
%                 if buildings(index, 7) ~=0 && buildings(index, 7) ~= facets(blocks(i, 7), 4)
%                     disp('sth went probably wrong, neighbouring blocks appear to belong to different buildings')
%                 end
%                 buildings(index, 7) = facets(blocks(i, 7), 4);
%             end
%             obj.buildings = buildings;

            %% assign azimuthal angle based on orientation
            wallaz = zeros(nfcts, 1); %wall azimuthal angles

            for i = 1:nfcts
                if facets(i, 1) == 1 %horizontal surface
                    wallaz(i) = obj.solaz; %azimuthal angle irrelevant for horizontal surface. Just set it to "solaz" so it passes self-shading test.
                elseif facets(i, 1) == 2 % west
                    wallaz(i) = 270;
                elseif facets(i, 1) == 3 % east
                    wallaz(i) = 90;
                elseif facets(i, 1) == 4 % north
                    wallaz(i) = 0;
                elseif facets(i,1) == 5
                    wallaz(i) = 180;  %=5, south
                end
            end

            solarzenith = obj.solarzenith; solaz = obj.solaz;
            %% vector to sun
            x = sin(solarzenith / 360 * 2 * pi) * cos((solaz - 90) / 360 * 2 * pi);  %-90 since it's from north, not from east (= our x coordinate)
            y = sin(solarzenith / 360 * 2 * pi) * sin((solaz + 90) / 360 * 2 * pi);
            z = cos(solarzenith / 360 * 2 * pi);
            v1 = [x, y, z];

            % facets is the first 4 columns of obj.facets so can just replace it
            % with obj.facets where it appears.



            %% create blocks to test for intersection
            % coordinates in physical space (not indeces)
            % blocks is blocks so can just replace it with obj.blocks
            %preprocessing.addvar(obj, 'blocks_phys', zeros(obj.nblocks + obj.nboundingwallfacets, 6));
            blocks_phys = zeros(nblocks + nboundingwallfacets, 6);
            for k = 1:nblocks
                xl = xh(blocks(k, 1));
                xu = xh(blocks(k,2) + 1);
                yl = yh(blocks(k,3));
                yu = yh(blocks(k,4) + 1);
                zl = zh(blocks(k,5) + 1);
                zu = zh(blocks(k,6) + 2);
                blocks_phys(k, :) = [xl, xu, yl, yu, zl, zu];
            end

            for k = 1:nboundingwallfacets
                fi3 = facets(k + (nfcts - nfloorfacets - nboundingwallfacets), 1);
                il = boundingwallfacets(k, 6);
                iu = boundingwallfacets(k, 7);
                jl = boundingwallfacets(k, 8);
                ju = boundingwallfacets(k, 9);
                kl = boundingwallfacets(k, 10) + 1;
                ku = boundingwallfacets(k, 11) + 1;

                if (fi3 == 2)
                    xl = xh(end);
                    xu = xh(end) + 0.1; %to make the bounding wall 3D, give them a thickness
                    yl = yh(jl);
                    yu = yh(ju + 1);
                    zl = zh(kl + 1);
                    zu = zh(ku + 2);
                elseif (fi3 == 3)
                    xl = xh(1) - 0.1; %to make the bounding wall 3D
                    xu = xh(1);
                    yl = yh(jl);
                    yu = yh(ju + 1);
                    zl = zh(kl + 1);
                    zu = zh(ku + 2);
                elseif (fi3 == 4)
                    xl = xh(il);
                    xu = xh(iu + 1);
                    yl = yh(1) - 0.1; %to make the bounding wall 3D
                    yu = yh(1);
                    zl = zh(kl + 1);
                    zu = zh(ku + 2);
                else %if (fi==5)
                    xl = xh(il);
                    xu = xh(iu + 1);
                    yl = yh(end);
                    yu = yh(end) + 0.1; %to make the bounding wall 3D
                    zl = zh(kl + 1);
                    zu = zh(ku + 2);
                end
                blocks_phys(k + nblocks, :) = [xl, xu, yl, yu, zl, zu];
            end

            preprocessing.addvar(obj, 'blocks_phys', blocks_phys);
            boundingwalls = boundingwallfacets(:, 6:end);

            gl = zeros(nfloorfacets, 6); %[xl xu yl yu zl zu] %space coordinates of floors
            for j = 1:nfloorfacets
                xl = xh(floorfacets(j, 6));
                xu = xh(floorfacets(j, 7) + 1);
                yl = yh(floorfacets(j, 8));
                yu = yh(floorfacets(j, 9) + 1);
                zl = zh(0 + 1); %to make the floors 3D (not really necessary)
                zu = zh(0 + 2);
                gl(j, :) = [xl, xu, yl, yu, zl, zu];
            end
            floors = floorfacets(:, 6:end);
            %bl = blocks_phys;


            %% Calculate Ray-Block-Intersection
            wsl = ones(nfcts, 1); %potentially sunlit surface (not self shaded)
            asl = zeros(nfcts, 1);
            walltheta = zeros(nfcts, 1); %theta

            for i = 1:nfcts
                if abs(wallaz(i) - solaz) >= 90  %angle between sun and facet-normal larger 90
                    wsl(i) = 0; %self shading
                    continue
                end

                % get facet center and corners
                cornm = obj.cornm;
                delta = 0.01;
                [ndim, ~, co] = preprocessing.detsub(i,facets,blocks,floors,boundingwalls,cornm,nblocks,xh,yh,zh,delta);
                %ndim=(dim1*dim2), number of cells of that facet
                %co = returns xyz-coordinates of center and 4 corners clockwise from
                %bottom left slightly shifter and
                %returns xyz-coordinates of center and 4 corners clockwise from
                %bottom left

                %count how much area cannot see sun (i.e. view is blocked)
                [ as ] = preprocessing.prblckd(i,-999,co,ndim,true,v1,-999,-999,facets,centerweight,cornerweight,nblocks,nboundingwallfacets,blocks_phys);
                asl(i) = 1 - as; %fraction of facet in sunlight
            end

            %walltheta is called Xi in the thesis
            %wallaz is called Omega_i
            for i = 1:nfcts
                if wsl(i) == 0 %shaded
                    walltheta(i) = 90;
                elseif facets(i,1) == 1 %horizontal
                    walltheta(i) = 0;
                else
                    walltheta(i) = abs(solaz - wallaz(i));
                end
            end

            I = obj.I;
            %% Direct shortwave
            Sdir = zeros(nfcts, 1);
            for i = 1:nfcts  %direct from sky
                if facets(i, 1) == 1 %horizontal
                    phi = 0;
                else
                    phi = 90; %vertical
                end
                Sdir(i) = I * cos((solarzenith - phi) / 360 * 2 * pi) * cos(walltheta(i) / 360 * 2 * pi) * asl(i);
            end

            preprocessing.addvar(obj, 'asl', asl);
            preprocessing.addvar(obj, 'Sdir', Sdir);
        end
        
        function generateBlocksPhys(obj)
            %% create blocks to test for intersection (not floors)
            % coordinates in physical space (not indeces)
            % blocks is blocks so can just replace it with obj.blocks
            %preprocessing.addvar(obj, 'blocks_phys', zeros(obj.nblocks + obj.nboundingwallfacets, 6));
            
            blocks = obj.blocks;
            nblocks = obj.nblocks;
            nfcts = obj.nfcts;
            nboundingwallfacets = obj.nboundingwallfacets;
            boundingwallfacets = obj.boundingwallfacets;
            floorfacets = obj.floorfacets;
            nfloorfacets = obj.nfloorfacets;
            xh = obj.xh; yh = obj.yh; zh = obj.zh;
            
            blocks_phys = zeros(nblocks + nboundingwallfacets, 6);
                        
            for k = 1:nblocks
                xl = xh(blocks(k, 1));
                xu = xh(blocks(k,2) + 1);
                yl = yh(blocks(k,3));
                yu = yh(blocks(k,4) + 1);
                zl = zh(blocks(k,5) + 1);
                zu = zh(blocks(k,6) + 2);
                blocks_phys(k, :) = [xl, xu, yl, yu, zl, zu];
            end

            for k = 1:nboundingwallfacets
                fi3 = facets(k + (nfcts - nfloorfacets - nboundingwallfacets), 1);
                il = boundingwallfacets(k, 6);
                iu = boundingwallfacets(k, 7);
                jl = boundingwallfacets(k, 8);
                ju = boundingwallfacets(k, 9);
                kl = boundingwallfacets(k, 10) + 1;
                ku = boundingwallfacets(k, 11) + 1;

                if (fi3 == 2)
                    xl = xh(end);
                    xu = xh(end) + 0.1; %to make the bounding wall 3D, give them a thickness
                    yl = yh(jl);
                    yu = yh(ju + 1);
                    zl = zh(kl + 1);
                    zu = zh(ku + 2);
                elseif (fi3 == 3)
                    xl = xh(1) - 0.1; %to make the bounding wall 3D
                    xu = xh(1);
                    yl = yh(jl);
                    yu = yh(ju + 1);
                    zl = zh(kl + 1);
                    zu = zh(ku + 2);
                elseif (fi3 == 4)
                    xl = xh(il);
                    xu = xh(iu + 1);
                    yl = yh(1) - 0.1; %to make the bounding wall 3D
                    yu = yh(1);
                    zl = zh(kl + 1);
                    zu = zh(ku + 2);
                else %if (fi==5)
                    xl = xh(il);
                    xu = xh(iu + 1);
                    yl = yh(end);
                    yu = yh(end) + 0.1; %to make the bounding wall 3D
                    zl = zh(kl + 1);
                    zu = zh(ku + 2);
                end
                blocks_phys(k + nblocks, :) = [xl, xu, yl, yu, zl, zu];
            end
            
            floors_phys = zeros(nfloorfacets, 6); %[xl xu yl yu zl zu] %space coordinates of floors
            for j = 1:nfloorfacets
                xl = xh(floorfacets(j, 6));
                xu = xh(floorfacets(j, 7) + 1);
                yl = yh(floorfacets(j, 8));
                yu = yh(floorfacets(j, 9) + 1);
                zl = zh(0 + 1); %to make the floors 3D (not really necessary)
                zu = zh(0 + 2);
                floors_phys(j, :) = [xl, xu, yl, yu, zl, zu];
            end

            preprocessing.addvar(obj, 'blocks_phys', blocks_phys);
            preprocessing.addvar(obj, 'floors_phys', floors_phys);
        end
        
        function readWeatherFile(obj)
            %range = ['A3629:Q3652', num2str(1+obj.nweatherdata)];
            timedepdata = readtable(obj.weatherfname);
            ids = find(table2array(timedepdata(:,1)) == 10600); % June - change later
            preprocessing.addvar(obj, 'timedeptime', table2array(timedepdata(ids, 'TIME')))
            preprocessing.addvar(obj, 'timedepTair', table2array(timedepdata(ids, 'TAIR')) + 273.15)
            preprocessing.addvar(obj, 'timedepI', table2array(timedepdata(ids, 'HELIOM')))
            preprocessing.addvar(obj, 'timedepDsky', table2array(timedepdata(ids, 'DIFSOLAR')))
            preprocessing.addvar(obj, 'timedepLsky', table2array(timedepdata(ids, 'IRSKYT')))
            preprocessing.addvar(obj, 'timedepzenith', table2array(timedepdata(ids, 'SOLAR')))
            preprocessing.addvar(obj, 'timedepazimuth', table2array(timedepdata(ids, 'SOLAR_1'))+90)
        end
        
        function vsolcAccelerated_td(obj)
            nblocks = obj.nblocks;
            facets = obj.facets; nfcts  = obj.nfcts;
            nboundingwallfacets = obj.nboundingwallfacets; %nbw = obj.nboundingwallfacets;
            centerweight = obj.centerweight; cornerweight = obj.cornerweight;
            blocks_phys = obj.blocks_phys;
            facetdim = obj.facetdim;
            facetcorner = obj.facetcorner;
              
            %% Loop over solar positions
            %wsl_array = ones(nfcts, npos); %potentially sunlit surface (not self shaded)
            asl_td = zeros(nfcts, obj.nweatherdata);
            %walltheta_array = zeros(nfcts, npos); %theta
            Sdir_td = zeros(nfcts, obj.nweatherdata);
            %wallaz_array = zeros(nfcts, npos); %wall azimuthal angles
            In_td = zeros(1, obj.nweatherdata);

            timedepazimuth_shift = circshift(obj.timedepazimuth, -obj.starttime);
            timedepzenith_shift = circshift(obj.timedepzenith, -obj.starttime);
            timedepI_shift = circshift(obj.timedepI, -obj.starttime);
            
            min_blocks_x = blocks_phys(:,1); min_blocks_y = blocks_phys(:,3); min_blocks_z = blocks_phys(:,5);
            max_blocks_x = blocks_phys(:,2); max_blocks_y = blocks_phys(:,4); max_blocks_z = blocks_phys(:,6);  
            
            for p = 1:obj.nweatherdata
                p
                solaz = timedepazimuth_shift(p) - obj.xazimuth + 90; % Correct for rotation of buildings
                solarzenith = timedepzenith_shift(p);
                I = timedepI_shift(p);
                
                if solarzenith > 89
                    continue;
                end
                
                %% assign azimuthal angle based on orientation
                wallaz = zeros(nfcts, 1); %wall azimuthal angles
                
%                 for i = 1:nfcts
%                     if facets(i, 1) == 1 %horizontal surface
%                         wallaz(i) = solaz; %azimuthal angle irrelevant for horizontal surface. Just set it to "solaz" so it passes self-shading test.
%                     elseif facets(i, 1) == 2 % west
%                         wallaz(i) = 270;
%                     elseif facets(i, 1) == 3 % east
%                         wallaz(i) = 90;
%                     elseif facets(i, 1) == 4 % north
%                         wallaz(i) = 0;
%                     elseif facets(i,1) == 5
%                         wallaz(i) = 180;  %=5, south
%                     end
%                 end
                
                for i = 1:nfcts
                    if facets(i, 1) == 1 %horizontal surface
                        wallaz(i) = solaz; %azimuthal angle irrelevant for horizontal surface. Just set it to "solaz" so it passes self-shading test.
                    elseif facets(i, 1) == 2 % west
                        wallaz(i) = 270 - obj.xazimuth + 90;
                    elseif facets(i, 1) == 3 % east
                        wallaz(i) = 90 - obj.xazimuth + 90;
                    elseif facets(i, 1) == 4 % north
                        wallaz(i) = 0 - obj.xazimuth + 90;
                    elseif facets(i,1) == 5 % south
                        wallaz(i) = 180 - obj.xazimuth + 90;
                    elseif facets(i,1) == 6 % bot
                        wallaz(i) = solaz;
                    end
                end
                
                %% vector to sun
                x = sin(solarzenith / 360 * 2 * pi) * cos((solaz - 90) / 360 * 2 * pi);  %-90 since it's from north, not from east (= our x coordinate)
                y = sin(solarzenith / 360 * 2 * pi) * sin((solaz + 90) / 360 * 2 * pi);
                z = cos(solarzenith / 360 * 2 * pi);
                v1 = [x, y, z];
                                           
                %% Calculate Ray-Block-Intersection
                wsl = ones(nfcts, 1); %potentially sunlit surface (not self shaded)
                asl = zeros(nfcts, 1);
                walltheta = zeros(nfcts, 1); %theta
                
                parfor i = 1:nfcts
                    if abs(wallaz(i) - solaz) >= 90  %angle between sun and facet-normal larger 90
                        wsl(i) = 0; %self shading
                        continue
                    end
                    
                    % get facet center and corners
                    %cornm = obj.cornm;
                    %delta = 0.01;
                    ndim = facetdim(i);
                    co = squeeze(facetcorner(i,:,:));
                    %[ndim, ~, co] = preprocessing.detsub(i,facets,blocks,floors,boundingwalls,cornm,nblocks,xh,yh,zh,delta);
                    %ndim=(dim1*dim2), number of cells of that facet
                    %co = returns xyz-coordinates of center and 4 corners clockwise from
                    %bottom left slightly shifter and
                    %returns xyz-coordinates of center and 4 corners clockwise from
                    %bottom left
                    
                    %count how much area cannot see sun (i.e. view is blocked)
                    %[ as ] = preprocessing.prblckd(i,-999,co,ndim,true,v1,-999,-999,facets,centerweight,cornerweight,nblocks,nboundingwallfacets,blocks_phys);
                    
                    bi=facets(i,3); %block index
                    ci=facets(i,4); %building index (-1 for roads, -101 for bounding wall)
                    flag = 0;
                    testcrit = 0;
                    as = 0;
                    if ndim>testcrit  %also test corner, otherwise only test center
                        for k=1:5
                            facetpoint=co(k,:);
                            
                            or = facetpoint;
                            d = v1;
                            [tmin_vec,flag_vec] = rayBoxVec(min_blocks_x, min_blocks_y, min_blocks_z, max_blocks_x, max_blocks_y, max_blocks_z, or(1), or(2), or(3), d(1), d(2), d(3));                    

                            flag_vec(bi) = 0;
                            intersection = repmat(facetpoint, [length(tmin_vec),1]) + tmin_vec * v1;
                            flag_vec(intersection(:,3) < facetpoint(3)) = 0;
                            
                            if any(flag_vec)
                              if k==1
                                  as=centerweight;
                              else
                                  as=as+cornerweight;
                              end
                            end
                            
%                             for n=1:nblocks %check if any block intersects
%                                 if ((n==bi) && (ci>0))  %don't intersect with itself
%                                     continue
%                                 end
%                                 
%                                 %[flag,dint] = preprocessing.rbi(facetpoint, v1, blocks_phys(n,:));
%                                 [flag,dint] = rayBoxIntersection(facetpoint, v1, blocks_phys(n,:));
%                                 
%                                 intersection = facetpoint + dint*v1;
%                                 if intersection(3)<facetpoint(3) %downstream direction of sun, thus not blocking the sun
%                                     flag=0;
%                                 end
%                                 
%                                 if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
%                                     if k==1
%                                         as=centerweight;
%                                     else
%                                         as=as+cornerweight;
%                                     end
%                                     %disp(['Loop: blocked by ' num2str(n)])
%                                     break %get out of inner for loop
%                                 end
%                             end
%                             
%                             if ~flag
%                                 for m=1:nboundingwallfacets %check if any bounding wall intersects
%                                     
%                                     %[flag,dint] = preprocessing.rbi(facetpoint, v1, blocks_phys(m+nblocks,:));
%                                     [flag,dint] = rayBoxIntersection(facetpoint, v1, blocks_phys(m+nblocks,:));
%                                     
%                                     intersection = facetpoint + dint*v1;
%                                     
%                                     if intersection(3)<facetpoint(3) %downstream direction of sun
%                                         flag=0;
%                                     end
%                                     if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
%                                         if k==1 %it's the center of the facet
%                                             as=centerweight;
%                                         else
%                                             as=as+cornerweight;
%                                         end
%                                         break %get out of inner for loop
%                                     end
%                                 end
%                             end
                         end
                    else
                        % %disp('small area, checking only center')
                        % %disp(['facet: ' num2str(i)])
                        facetpoint=co(1,:);  %center of facet
                        for n=1:nblocks %check if any block intersects
                            if ((n==bi) && (ci>0))  %don't intersect with itself
                                continue
                            end
                            
                            %[flag,dint] = preprocessing.rbi(facetpoint, v1, blocks_phys(n,:));
                            [flag,dint] = rayBoxIntersection(facetpoint, v1, blocks_phys(n,:));
                            
                            intersection = facetpoint + dint*v1;
                            if intersection(3)<facetpoint(3) %downstream direction of sun, thus not blocking the sun
                                % intersectionss(:)=0;
                                flag=0;
                            end
                            if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                                as=1;
                                break %get out of inner for loop
                            end
                        end
                        
                        if ~flag
                            for k=1:nboundingwallfacets %check if any bounding wall intersects
                                
                                %[flag,dint] = preprocessing.rbi(facetpoint, v1, blocks_phys(k+nblocks,:));
                                [flag,dint] = rayBoxIntersection(facetpoint, v1, blocks_phys(k+nblocks,:));
                                
                                intersection = facetpoint + dint*v1;
                                
                                if intersection(3)<facetpoint(3) %downstream direction of sun
                                    flag=0;
                                end
                                if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                                    as=1;
                                    break %get out of inner for loop
                                end
                            end
                        end
                    end
                    
                    asl(i) = 1 - as; %fraction of facet in sunlight
                end
                
                %walltheta is called Xi in the thesis
                %wallaz is called Omega_i
                for i = 1:nfcts
                    if wsl(i) == 0 %shaded
                        walltheta(i) = 90;
                    elseif facets(i,1) == 1 %horizontal
                        walltheta(i) = 0;
                    else
                        walltheta(i) = abs(solaz - wallaz(i));
                    end
                end
                               
                %% Direct shortwave
                Sdir = zeros(nfcts, 1);
                for i = 1:nfcts  %direct from sky
                    if facets(i, 1) == 1 %horizontal
                        phi = 0;
                    else
                        phi = 90; %vertical
                    end
                    Sdir(i) = I * cos((solarzenith - phi) / 360 * 2 * pi) * cos(walltheta(i) / 360 * 2 * pi) * asl(i);
                end
                
                Sdir_td(:,p) = Sdir;
                asl_td(:,p) = asl;
                In_td(p) = I;
                
            end
            
            preprocessing.addvar(obj, 'asl_td', asl_td);
            preprocessing.addvar(obj, 'Sdir_td', Sdir_td);
            preprocessing.addvar(obj, 'In_td', In_td);
        end
        
        function vsolcAccelerated(obj)
            xh = obj.xh;
            yh = obj.yh;
            zh = obj.zh;
            blocks = obj.blocks; nblocks = obj.nblocks;
            facets = obj.facets; nfcts  = obj.nfcts;
            floorfacets = obj.floorfacets; nfloorfacets = obj.nfloorfacets;
            boundingwallfacets = obj.boundingwallfacets; nboundingwallfacets = obj.nboundingwallfacets; %nbw = obj.nboundingwallfacets;
            buildings = obj.buildings;
            centerweight = obj.centerweight; cornerweight = obj.cornerweight;
               
%           Ignore reference to buildings as they are unused in preprocessing.
%             buildings(:, 7) = 0;
%             for i = 1:nblocks
%                 index = find(buildings(:, 1) <= buildings(i, 1) & buildings(:, 2) >= buildings(i, 2) & buildings(:, 3) <= buildings(i, 3) & buildings(:, 4) >= buildings(i, 4));
%                 %if buildings(index, 7) ~=0 && buildings(index, 7) ~= facets(blocks(i, 7), 4)
%                     %disp('sth went probably wrong, neighbouring blocks appear to belong to different buildings')
%                 %end
%                 buildings(index, 7) = facets(blocks(i, 7), 4);
%             end
%             obj.buildings = buildings;

            %% assign azimuthal angle based on orientation
            wallaz = zeros(nfcts, 1); %wall azimuthal angles
            solarzenith = obj.solarzenith; solaz = obj.solaz;
            
            for i = 1:nfcts
                if facets(i, 1) == 1 %horizontal surface
                    wallaz(i) = solaz; %azimuthal angle irrelevant for horizontal surface. Just set it to "solaz" so it passes self-shading test.
                elseif facets(i, 1) == 2 % west
                    wallaz(i) = 270 - obj.xazimuth + 90;
                elseif facets(i, 1) == 3 % east
                    wallaz(i) = 90 - obj.xazimuth + 90;
                elseif facets(i, 1) == 4 % north
                    wallaz(i) = 0 - obj.xazimuth + 90;
                elseif facets(i,1) == 5 % south
                    wallaz(i) = 180 - obj.xazimuth + 90;
                elseif facets(i,1) == 6 % bot
                    wallaz(i) = inf;
                end
            end

            
            %% vector to sun
            x = sin(solarzenith / 360 * 2 * pi) * cos((solaz - 90) / 360 * 2 * pi);  %-90 since it's from north, not from east (= our x coordinate)
            y = sin(solarzenith / 360 * 2 * pi) * sin((solaz + 90) / 360 * 2 * pi);
            z = cos(solarzenith / 360 * 2 * pi);
            v1 = [x, y, z];

            % facets is the first 4 columns of obj.facets so can just replace it
            % with obj.facets where it appears.

            %boundingwalls = boundingwallfacets(:, 6:end);


            %floors = floorfacets(:, 6:end);
            
            
            blocks_phys = obj.blocks_phys;
            %floors_phys = obj.floors_phys;

            %% Calculate Ray-Block-Intersection
            wsl = ones(nfcts, 1); %potentially sunlit surface (not self shaded)
            asl = zeros(nfcts, 1);
            walltheta = zeros(nfcts, 1); %theta

            for i = 1:nfcts
                if abs(wallaz(i) - solaz) >= 90  %angle between sun and facet-normal larger 90
                    wsl(i) = 0; %self shading
                    continue
                end

                % get facet center and corners
                %cornm = obj.cornm;
                %delta = 0.01;
                ndim = obj.facetdim(i);
                co = squeeze(obj.facetcorner(i,:,:));
                %[ndim, ~, co] = preprocessing.detsub(i,facets,blocks,floors,boundingwalls,cornm,nblocks,xh,yh,zh,delta);
                %ndim=(dim1*dim2), number of cells of that facet
                %co = returns xyz-coordinates of center and 4 corners clockwise from
                %bottom left slightly shifter and
                %returns xyz-coordinates of center and 4 corners clockwise from
                %bottom left

                %count how much area cannot see sun (i.e. view is blocked)
                %[ as ] = preprocessing.prblckd(i,-999,co,ndim,true,v1,-999,-999,facets,centerweight,cornerweight,nblocks,nboundingwallfacets,blocks_phys);
                
                bi=facets(i,3); %block index
                ci=facets(i,4); %building index (-1 for roads, -101 for bounding wall)
                flag = 0;
                testcrit = 0;
                as = 0;
                if ndim>testcrit  %also test corner, otherwise only test center
                    for k=1:5
                        facetpoint=co(k,:);
                        
%                         min = [blocks_phys(:,1), blocks_phys(:,3), blocks_phys(:,5)];
%                         max = [blocks_phys(:,2), blocks_phys(:,4), blocks_phys(:,6)];
%                         
%                         or = facetpoint;
%                         d = v1;
%                         
%                         [tmin_vec,flag_vec] = rayBoxVec(min(:,1)', min(:,2)', min(:,3)', ...
%                             max(:,1)', max(:,2)', max(:,3)', ...
%                             or(1), or(2), or(3)  , ...
%                             d(1), d(2), d(3));
% 
%                         flag_vec(bi) = 0;
%                         
%                         %size(tmin_vec)
%                         %size(v1)
%                         intersection = facetpoint + tmin_vec' * v1;
%                         %size(intersection)
%                         flag_vec(intersection(:,3) < facetpoint(3)) = 0;
% 
%                         if any(flag_vec)
%                             if k==1
%                                 as=centerweight;
%                             else
%                                 %disp(['Vector: blocked by ' num2str(find(flag_vec,1,'first'))])
%                                 as=as+cornerweight;
%                             end
%                         end
                        
                        for n=1:nblocks %check if any block intersects
                            if ((n==bi) && (ci>0))  %don't intersect with itself
                                continue
                            end
                            
                            %[flag,dint] = preprocessing.rbi(facetpoint, v1, blocks_phys(n,:));
                            [flag,dint] = rayBoxIntersection(facetpoint, v1, blocks_phys(n,:));
                            
                            intersection = facetpoint + dint*v1;
                            if intersection(3)<facetpoint(3) %downstream direction of sun, thus not blocking the sun
                                flag=0;
                            end
                            
                            if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                                if k==1
                                    as=centerweight;
                                else
                                    as=as+cornerweight;
                                end
                                %disp(['Loop: blocked by ' num2str(n)])
                                break %get out of inner for loop
                            end
                        end
                        
                        if ~flag
                            for m=1:nboundingwallfacets %check if any bounding wall intersects
                                
                                %[flag,dint] = preprocessing.rbi(facetpoint, v1, blocks_phys(m+nblocks,:));
                                [flag,dint] = rayBoxIntersection(facetpoint, v1, blocks_phys(m+nblocks,:));
                                
                                intersection = facetpoint + dint*v1;
                                
                                if intersection(3)<facetpoint(3) %downstream direction of sun
                                    flag=0;
                                end
                                if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                                    if k==1 %it's the center of the facet
                                        as=centerweight;
                                    else
                                        as=as+cornerweight;
                                    end
                                    break %get out of inner for loop
                                end
                            end
                        end
                    end
                else
                    % %disp('small area, checking only center')
                    % %disp(['facet: ' num2str(i)])
                    facetpoint=co(1,:);  %center of facet
                    for n=1:nblocks %check if any block intersects
                        if ((n==bi) && (ci>0))  %don't intersect with itself
                            continue
                        end
                        
                        %[flag,dint] = preprocessing.rbi(facetpoint, v1, blocks_phys(n,:));
                        [flag,dint] = rayBoxIntersection(facetpoint, v1, blocks_phys(n,:));
                        
                        intersection = facetpoint + dint*v1;
                        if intersection(3)<facetpoint(3) %downstream direction of sun, thus not blocking the sun
                            % intersectionss(:)=0;
                            flag=0;
                        end
                        if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                            as=1;
                            break %get out of inner for loop
                        end
                    end
                    
                    if ~flag
                        for k=1:nboundingwallfacets %check if any bounding wall intersects
                            
                            %[flag,dint] = preprocessing.rbi(facetpoint, v1, blocks_phys(k+nblocks,:));
                            [flag,dint] = rayBoxIntersection(facetpoint, v1, blocks_phys(k+nblocks,:));
                            
                            intersection = facetpoint + dint*v1;
                            
                            if intersection(3)<facetpoint(3) %downstream direction of sun
                                flag=0;
                            end
                            if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                                as=1;
                                break %get out of inner for loop
                            end
                        end
                    end
                end
                
                asl(i) = 1 - as; %fraction of facet in sunlight
            end

            %walltheta is called Xi in the thesis
            %wallaz is called Omega_i
            for i = 1:nfcts
                if wsl(i) == 0 %shaded
                    walltheta(i) = 90;
                elseif facets(i,1) == 1 %horizontal
                    walltheta(i) = 0;
                else
                    walltheta(i) = abs(solaz - wallaz(i));
                end
            end

            I = obj.I;
            %% Direct shortwave
            Sdir = zeros(nfcts, 1);
            for i = 1:nfcts  %direct from sky
                if facets(i, 1) == 1 %horizontal
                    phi = 0;
                else
                    phi = 90; %vertical
                end
                Sdir(i) = I * cos((solarzenith - phi) / 360 * 2 * pi) * cos(walltheta(i) / 360 * 2 * pi) * asl(i);
            end

            preprocessing.addvar(obj, 'asl', asl);
            preprocessing.addvar(obj, 'Sdir', Sdir);
        end
        
        function plot_Sdir(obj)
            figure
            %cmap = colormap('parula');

            if ~obj.lflat
                top = 1; west = 2; east = 3; north = 4; south = 5; bot = 6;
                
                for i = 1:obj.nfcts
                    clr = [obj.Sdir(i), obj.Sdir(i), obj.Sdir(i)] / obj.I;
                    
                    il = obj.facets(i, 6); iu = obj.facets(i, 7);
                    jl = obj.facets(i, 8); ju = obj.facets(i, 9);
                    kl = obj.facets(i, 10); ku = obj.facets(i, 11);
                    if i <= obj.nblockfcts
                        switch obj.facets(i, 1)
                            case {east, west}
                                x = [obj.xh(il), obj.xh(il), obj.xh(il), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(ju), obj.yh(ju)];
                                if kl == 0 && ku == 0
                                    z = [0, 0, 0, 0];
                                elseif kl == 0 && ku ~= 0
                                    z = [0, obj.zh(ku), obj.zh(ku), 0];
                                else
                                    z = [obj.zh(kl), obj.zh(ku), obj.zh(ku), obj.zh(kl)];
                                end
                            case {north, south}
                                x = [obj.xh(il), obj.xh(iu), obj.xh(iu), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(jl), obj.yh(jl)];
                                if kl == 0 && ku == 0
                                    z = [0, 0, 0, 0];
                                elseif kl == 0 && ku ~= 0
                                    z = [0, 0, obj.xh(ku), obj.xh(ku)];
                                else
                                    z = [obj.zh(kl), obj.zh(kl), obj.zh(ku), obj.zh(ku)];
                                end
                            case {top, bot}
                                x = [obj.xh(il), obj.xh(iu), obj.xh(iu), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(ju), obj.yh(ju)];

                                if kl == 0
                                    z = [0, 0, 0, 0];
                                else
                                    z = [obj.zh(kl), obj.zh(kl), obj.zh(kl), obj.zh(kl)];
                                end
                        end

                    elseif i <= obj.nblockfcts + obj.nboundingwallfacets
                        switch obj.facets(i, 1)
                            case east
                                ju = ju + 1;
                                kl = kl + 1;
                                ku = ku + 2;
                                x = [obj.xh(il), obj.xh(il), obj.xh(il), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(ju), obj.yh(ju)];
                                if kl == 0 && ku == 0
                                    z = [0, 0, 0, 0];
                                elseif kl == 0 && ku ~= 0
                                    z = [0, obj.zh(ku), obj.zh(ku), 0];
                                else
                                    z = [obj.zh(kl), obj.zh(ku), obj.zh(ku), obj.zh(kl)];
                                end
                            case west
                                il = il + 1;
                                iu = iu + 1;
                                ju = ju + 1;
                                kl = kl + 1;
                                ku = ku + 2;
                                x = [obj.xh(il), obj.xh(il), obj.xh(il), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(ju), obj.yh(ju)];
                                if kl == 0 && ku == 0
                                    z = [0, 0, 0, 0];
                                elseif kl == 0 && ku ~= 0
                                    z = [0, obj.zh(ku), obj.zh(ku), 0];
                                else
                                    z = [obj.zh(kl), obj.zh(ku), obj.zh(ku), obj.zh(kl)];
                                end


                            case north
                                iu = iu + 1;
                                kl = kl + 1;
                                ku = ku + 2;
                                x = [obj.xh(il), obj.xh(iu), obj.xh(iu), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(jl), obj.yh(jl)];
                                if kl == 0 && ku == 0
                                    z = [0, 0, 0, 0];
                                elseif kl == 0 && ku ~= 0
                                    z = [0, 0, obj.xh(ku), obj.xh(ku)];
                                else
                                    z = [obj.zh(kl), obj.zh(kl), obj.zh(ku), obj.zh(ku)];
                                end

                            case south
                                iu = iu + 1;
                                jl = jl + 1;
                                ju = ju + 1;
                                kl = kl + 1;
                                ku = ku + 2;
                                x = [obj.xh(il), obj.xh(iu), obj.xh(iu), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(jl), obj.yh(jl)];
                                if kl == 0 && ku == 0
                                    z = [0, 0, 0, 0];
                                elseif kl == 0 && ku ~= 0
                                    z = [0, 0, obj.xh(ku), obj.xh(ku)];
                                else
                                    z = [obj.zh(kl), obj.zh(kl), obj.zh(ku), obj.zh(ku)];
                                end


                            case {top, bot}
                                x = [obj.xh(il), obj.xh(iu), obj.xh(iu), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(ju), obj.yh(ju)];

                                if kl == 0
                                    z = [0, 0, 0, 0];
                                else
                                    z = [obj.zh(kl), obj.zh(kl), obj.zh(kl), obj.zh(kl)];
                                end
                        end
                    else
                        iu = iu + 1;
                        ju = ju + 1;
                        switch obj.facets(i, 1)
                            case {east, west}
                                x = [obj.xh(il), obj.xh(il), obj.xh(il), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(ju), obj.yh(ju)];
                                if kl == 0 && ku == 0
                                    z = [0, 0, 0, 0];
                                elseif kl == 0 && ku ~= 0
                                    z = [0, obj.zh(ku), obj.zh(ku), 0];
                                else
                                    z = [obj.zh(kl), obj.zh(ku), obj.zh(ku), obj.zh(kl)];
                                end
                            case {north, south}
                                x = [obj.xh(il), obj.xh(iu), obj.xh(iu), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(jl), obj.yh(jl)];
                                if kl == 0 && ku == 0
                                    z = [0, 0, 0, 0];
                                elseif kl == 0 && ku ~= 0
                                    z = [0, 0, obj.xh(ku), obj.xh(ku)];
                                else
                                    z = [obj.zh(kl), obj.zh(kl), obj.zh(ku), obj.zh(ku)];
                                end
                            case {top, bot}
                                x = [obj.xh(il), obj.xh(iu), obj.xh(iu), obj.xh(il)];
                                y = [obj.yh(jl), obj.yh(jl), obj.yh(ju), obj.yh(ju)];

                                if kl == 0
                                    z = [0, 0, 0, 0];
                                else
                                    z = [obj.zh(kl), obj.zh(kl), obj.zh(kl), obj.zh(kl)];
                                end
                        end
                    end

                    patch(x,y,z, clr,'FaceLighting','none','FaceAlpha',1);
                    d = [0, 0, 0]; a = 0.25;
                    switch(obj.facets(i, 1))
                        case top
                            d(3) = a * obj.dz(1);
                        case west
                            d(1) = -a * obj.dx(1);
                        case east
                            d(1) = a * obj.dx(1);
                        case south
                            d(2) = -a * obj.dy(1);
                        case north
                            d(2) = a * obj.dy(1);
                    end

                    %text(mean(x) + d(1), mean(y) + d(2), mean(z) + d(3), num2str(i), 'horizontalalignment', 'center')
                    %hold on
                    %title('Facets')
                end

            end
            view(180 - obj.solaz, 90-obj.solarzenith)
            set(gca,'ticklabelinterpreter','latex')
            xlabel('$x(\mathrm{m})$','interpreter','latex')
            ylabel('$y(\mathrm{m})$','interpreter','latex')
            zlabel('$z(\mathrm{m})$','interpreter','latex')
            set(gca,'BoxStyle','full','Box','on')
            daspect([1 1 1])
            xlim([0 obj.xh(end)])
            ylim([0 obj.yh(end)])
            zlim([0 obj.zh(end)])
            set(gca, 'FontSize', 12)
            grid on
        end
        
        function plot_shading(obj)
            scale = 2;
            scalef = 1.5;
            h = figure;
            view(3);
            xlabel('x')
            ylabel('y')
            zlabel('z')
            axis equal
            set(gcf,'units','centimeters','position', [2 2 2+14.5*scale 2+14.5*scale/2]);
            set(h,'PaperPosition',[2 2 2+14.5*scale 2+14.5*scale/2]);
            set(h,'PaperUnits','centimeters');
            set(h,'renderer','painters');

            %plotshadowline

            dx = ones(obj.imax, 1) * obj.dx;
            dy = ones(obj.jtot, 1) * obj.dy;
            dz = obj.dzf;

            for i=1:obj.nfcts
                bi = obj.facets(i,3); %block index

                fi = obj.facets(i,1); %facet index
                % facets is the first 4 columns of obj.facets
                % floors is floors, so floors(:, 1) = obj.floorfacets(:, 6)  - add 5
                % blocks is blocks so can just replace it with obj.blocks
                % boundingwalls is boundingwalls, so boundingwalls(:, 1) = obj.boundingwallfacets(:, 6) - add 5

                if (obj.facets(i,4) <= -100)  %it is a bounding wall
                    il = obj.boundingwallfacets(bi, 6);
                    iu = obj.boundingwallfacets(bi, 7);
                    jl = obj.boundingwallfacets(bi, 8);
                    ju = obj.boundingwallfacets(bi, 9);
                    kl = obj.boundingwallfacets(bi, 10) + 1;
                    ku = obj.boundingwallfacets(bi, 11) + 1;

                    if (fi == 2)
                        xl = obj.xf(iu) + 0.5 * dx(iu);
                        xu = obj.xf(iu) + 0.5 * dx(iu);
                        yl = obj.yf(jl) - 0.5 * dy(jl);
                        yu = obj.yf(ju) + 0.5 * dy(ju);
                        zl = obj.zf(kl) - 0.5 * dz(kl);
                        zu = obj.zf(ku) + 0.5 * dz(ku);
                    elseif (fi == 3)
                        xl = obj.xf(il) - dx(il);
                        xu = obj.xf(il) - 0.5 * dx(il);
                        yl = obj.yf(jl) - 0.5 * dy(jl);
                        yu = obj.yf(ju) + 0.5 * dy(ju);
                        zl = obj.zf(kl) - 0.5 * dz(kl);
                        zu = obj.zf(ku) + 0.5 * dz(ku);
                    elseif (fi==4)
                        xl = obj.xf(il) - 0.5 * dx(il);
                        xu = obj.xf(iu) + 0.5 * dx(iu);
                        yl = obj.yf(jl) - dy(jl);
                        yu = obj.yf(jl) - 0.5 * dy(jl);
                        zl = obj.zf(kl) - 0.5 * dz(kl);
                        zu = obj.zf(ku) + 0.5 * dz(ku);
                    else %if (fi==5)
                        xl = obj.xf(il) - 0.5 * dx(il);
                        xu = obj.xf(iu) + 0.5 * dx(iu);
                        yl = obj.yf(ju) + 0.5 * dy(ju);
                        yu = obj.yf(ju) + dy(ju);
                        zl = obj.zf(kl) - 0.5 * dz(kl);
                        zu = obj.zf(ku) + 0.5 * dz(ku);
                    end
                    disp([il, iu, jl, ju, kl, ku])
                    disp([xl, xu, yl, yu, zl, zu])

                elseif (obj.facets(i,4) < 0 && obj.facets(i,4) > -100) %it is a floor, not a building
                    bif = bi - obj.nblocks;
                    il = obj.floorfacets(bif, 6);
                    iu = obj.floorfacets(bif, 7);
                    jl = obj.floorfacets(bif, 8);
                    ju = obj.floorfacets(bif, 9);
                    %xl = obj.xf(il) - 0.5 * dx(il);
                    xu = obj.xf(iu) + 0.5 * dx(iu);
                    yl = obj.yf(jl) - 0.5 * dy(jl);
                    yu = obj.yf(ju) + 0.5 * dy(ju);
                    zl = 0;
                    zu = 0;

                else
                    il = obj.blocks(bi, 1);
                    iu = obj.blocks(bi, 2);
                    jl = obj.blocks(bi, 3);
                    ju = obj.blocks(bi, 4);
                    kl = obj.blocks(bi, 5) + 1;
                    ku = obj.blocks(bi, 6) + 1;

                    xl = obj.xf(il) - 0.5 * dx(il);
                    xu = obj.xf(iu) + 0.5 * dx(iu);
                    yl = obj.yf(jl) - 0.5 * dy(jl);
                    yu = obj.yf(ju) + 0.5 * dy(ju);
                    zl = obj.zf(kl) - 0.5 * dz(kl);
                    zu = obj.zf(ku) + 0.5 * dz(ku);
                end

                switch obj.facets(i, 1)
                    case {1} %top
                        x = [xl, xu, xu, xl];
                        y = [yl, yl, yu, yu];
                        z = [zu, zu, zu, zu];
                    case {2} %west
                        x = [xl, xl, xl, xl];
                        y = [yl, yl, yu, yu];
                        z = [zl, zu, zu, zl];
                    case {3} %east
                        x = [xu, xu, xu, xu];
                        y = [yl, yl, yu, yu];
                        z = [zl, zu, zu, zl];
                    case {4} %north
                        x = [xl, xl, xu, xu];
                        y = [yu, yu, yu, yu];
                        z = [zl, zu, zu, zl];
                    case {5} %south
                        x = [xl, xl, xu, xu];
                        y = [yl, yl, yl, yl];
                        z = [zl, zu, zu, zl];
                end

                %shaded vs non-shaded
                greymax = 0.15; %maximum percent of blackness, [1 1 1] is white;
                cs = obj.asl(i) * [1-greymax 1-greymax 1-greymax] + [greymax greymax greymax];
                patch(x, y, z, cs,'FaceLighting','none');
                hold on

                %quiver3(face(92,1,2,1),face(92,1,2,2),face(92,1,2,3), v1(1),v1(2),v1(3), 8);
                %quiver3(face(92,1,3,1),face(92,1,3,2),face(92,1,3,3), v1(1),v1(2),v1(3), 8);
                % quiver3(face(92,2,2,1),face(92,2,2,2),face(92,2,2,3), v1(1),v1(2),v1(3), 8);
                %  quiver3(face(92,2,3,1),face(92,2,3,2),face(92,2,3,3), v1(1),v1(2),v1(3), 8);


                %   title('sunlit fraction')
            end

            %HS=scatter3(shadowscatter(:,1),shadowscatter(:,2),shadowscatter(:,3),'filled','dk');

            %                 ffv=134;
            %                 quiver3((xb(fctl(ffv, 6))+xb(fctl(ffv, 7)))/2,(yb(fctl(ffv, 8))+yb(fctl(ffv, 9)))/2,zb(fctl(ffv, 10))+0.05,v1(1),v1(2),v1(3),3,'r','LineWidth',2,'MaxHeadSize',1)
            %                 quiver3(xb(fctl(ffv, 6)),yb(fctl(ffv, 8)),zb(fctl(ffv, 10))+0.05,v1(1),v1(2),v1(3),1.5,'r','LineWidth',2,'MaxHeadSize',1)
            %                 quiver3(xb(fctl(ffv, 7)),yb(fctl(ffv, 8)),zb(fctl(ffv, 10))+0.05,v1(1),v1(2),v1(3),0.5,'r','LineWidth',2,'MaxHeadSize',1)
            %                 quiver3(xb(fctl(ffv, 6)),yb(fctl(ffv, 9)),zb(fctl(ffv, 10))+0.05,v1(1),v1(2),v1(3),3,'r','LineWidth',2,'MaxHeadSize',1)
            %                 quiver3(xb(fctl(ffv, 7)),yb(fctl(ffv, 9)),zb(fctl(ffv, 10))+0.05,v1(1),v1(2),v1(3),3,'r','LineWidth',2,'MaxHeadSize',1)
            %
            %                 greycm=zeros(9,3);
            %                 cticks=zeros(9,1);
            %                 ctickl=zeros(9,1);
            %                 for l=0:8
            %                     cs=l*0.125*[1-greymax 1-greymax 1-greymax] + [greymax greymax greymax];
            %                     cticks(l+1)=0.15+(l)*0.85/9+0.85/18;
            %                     ctickl(l+1)=l*0.125;
            %                     greycm(l+1,:)=cs;
            %                 end
            %                 cm=colormap(greycm);
            %                 c=colorbar;
            %                 c.Title.Interpreter='latex';
            %                 set(c, 'TickLabelInterpreter', 'latex')
            %                 set(c,'FontSize',12*scalef)
            %                 caxis([0.15 1])
            %                 c.Ticks=cticks
            %                 c.TickDirection='out'
            %                 c.TickLabels=num2cell(ctickl)
            %                 c.Title.String='$f_{\mathrm{e}}$ $[-]$'
            %                 %c.Position(1)=c.Position(1)+0.00
            %                 c.Position(2)=c.Position(2)+0.02
            %                 c.Position(4)=c.Position(4)-0.2
            %
            %                 view(3)
            %                 axis equal
            %                 xlim([0 xb(end)])
            %                 ylim([0 yb(end)])
            %                 % zlim([0 72])
            %                 zlim([0 max(zmax5)])
            %                 set(gca,'TickLabelInterpreter','latex')
            %                 set(gca,'FontSize',12*scalef)
            %                 %h1=gca;
            %                 % h1.Position=[0.08 0.1100 0.78 0.8150];
            %                 %caxis([0 71])
            %                 xlabel('x [-]','Interpreter','latex','FontSize',12*scalef)
            %                 ylabel('y [-]','Interpreter','latex','FontSize',12*scalef)
            %                 zlabel('z [-]','Interpreter','latex','FontSize',12*scalef)
            %                 set(gcf, 'Color', 'w');
            %export_fig facetshadowalt.eps


            %%
            scale=2;
            scalef=1.5;
            h=figure;
            %cm=colormap(MPL_gist_heat(31:end,:));
            %cm=colormap(MPL_viridis(30:end,:));
            %cm=colormap(flipud(perscolormaps('MPL_BuPu')));
            %cm=colormap(cm(1:end-2,:));

            set(gcf,'units','centimeters','position', [2 2 2+14.5*scale 2+14.5*scale/2]);
            set(h,'PaperPosition',[2 2 2+14.5*scale 2+14.5*scale/2]);
            set(h,'PaperUnits','centimeters');
            set(h,'renderer','painters');

            for i=1:obj.nfcts
                bi = obj.facets(i,3); %block index
                fi = obj.facets(i,1); %facet index
                % facets is the first 4 columns of obj.facets
                % floors is floors, so floors(:, 1) = obj.floorfacets(:, 6)  - add 5
                % blocks is blocks so can just replace it with obj.blocks
                % boundingwalls is boundingwalls, so boundingwalls(:, 1) = obj.boundingwallfacets(:, 6) - add 5
                if (obj.facets(i,4) < 0 && obj.facets(i,4) > -100) %it is a floor, not a building
                    bif = bi - obj.nblocks;
                    il = obj.floorfacets(bif, 6);
                    iu = obj.floorfacets(bif, 7);
                    jl = obj.floorfacets(bif, 8);
                    ju = obj.floorfacets(bif, 9);
                    xl = obj.xf(il) - 0.5 * dx(il);
                    xu = obj.xf(iu) + 0.5 * dx(iu);
                    yl = obj.yf(jl) - 0.5 * dy(jl);
                    yu = obj.yf(ju) + 0.5 * dy(ju);
                    zl = 0;
                    zu = 0;
                elseif (obj.facets(i,4) <= -100)  %it is a bounding wall
                    il = obj.boundingwallfacets(bi, 6);
                    iu = obj.boundingwallfacets(bi, 7);
                    jl = obj.boundingwallfacets(bi, 8);
                    ju = obj.boundingwallfacets(bi, 9);
                    kl = obj.boundingwallfacets(bi, 10) + 1;
                    ku = obj.boundingwallfacets(bi, 11) + 1;

                    if (fi == 2)
                        xl = obj.xf(iu) + 0.5 * dx(iu);
                        xu = obj.xf(iu) + 0.5 * dx(iu);
                        yl = obj.yf(jl) - 0.5 * dy(jl);
                        yu = obj.yf(ju) + 0.5 * dy(ju);
                        zl = obj.zf(kl) - 0.5 * dz(kl);
                        zu = obj.zf(ku) + 0.5 * dz(ku);
                    elseif (fi == 3)
                        xl = obj.xf(il) - dx(il);
                        xu = obj.xf(il) - 0.5 * dx(il);
                        yl = obj.yf(jl) - 0.5 * dy(jl);
                        yu = obj.yf(ju) + 0.5 * dy(ju);
                        zl = obj.zf(kl) - 0.5 * dz(kl);
                        zu = obj.zf(ku) + 0.5 * dz(ku);
                    elseif (fi==4)
                        xl = obj.xf(il) - 0.5 * dx(il);
                        xu = obj.xf(iu) + 0.5 * dx(iu);
                        yl = obj.yf(jl) - dy(jl);
                        yu = obj.yf(jl) - 0.5 * dy(jl);
                        zl = obj.zf(kl) - 0.5 * dz(kl);
                        zu = obj.zf(ku) + 0.5 * dz(ku);
                    else %if (fi==5)
                        xl = obj.xf(il) - 0.5 * dx(il);
                        xu = obj.xf(iu) + 0.5 * dx(iu);
                        yl = obj.yf(ju) + 0.5 * dy(ju);
                        yu = obj.yf(ju) + dy(ju);
                        zl = obj.zf(kl) - 0.5 * dz(kl);
                        zu = obj.zf(ku) + 0.5 * dz(ku);
                    end
                else
                    il = obj.blocks(bi, 1);
                    iu = obj.blocks(bi, 2);
                    jl = obj.blocks(bi, 3);
                    ju = obj.blocks(bi, 4);
                    kl = obj.blocks(bi, 5) + 1;
                    ku = obj.blocks(bi, 6) + 1;

                    xl = obj.xf(il) - 0.5 * dx(il);
                    xu = obj.xf(iu) + 0.5 * dx(iu);
                    yl = obj.yf(jl) - 0.5 * dy(jl);
                    yu = obj.yf(ju) + 0.5 * dy(ju);
                    zl = obj.zf(kl) - 0.5 * dz(kl);
                    zu = obj.zf(ku) + 0.5 * dz(ku);
                end

                switch obj.facets(i, 1)
                    case {1} %top
                        x = [xl, xu, xu, xl];
                        y = [yl, yl, yu, yu];
                        z = [zu, zu, zu, zu];
                    case {2} %west
                        x = [xl, xl, xl, xl];
                        y = [yl, yl, yu, yu];
                        z = [zl, zu, zu, zl];
                    case {3} %east
                        x = [xu, xu, xu, xu];
                        y = [yl, yl, yu, yu];
                        z = [zl, zu, zu, zl];
                    case {4} %north
                        x = [xl, xl, xu, xu];
                        y = [yu, yu, yu, yu];
                        z = [zl, zu, zu, zl];
                    case {5} %south
                        x = [xl, xl, xu, xu];
                        y = [yl, yl, yl, yl];
                        z = [zl, zu, zu, zl];
                end
                %shaded vs non-shaded
                %cs=Sdir(i)/I*[0.85 0.85 0.85] + [0.15 0.15 0.15];
                cm = colormap('parula');
                cs = cm(max(1, min(size(cm,1), round(obj.Sdir(i) / obj.I * size(cm, 1)))),:);

                patch(x,y,z, cs,'FaceLighting','none');
                hold on
                %        quiver3(face(92,1,2,1),face(92,1,2,2),face(92,1,2,3), v1(1),v1(2),v1(3), 1);
                %        quiver3(face(92,1,3,1),face(92,1,3,2),face(92,1,3,3), v1(1),v1(2),v1(3), 1);
                %        quiver3(face(92,2,2,1),face(92,2,2,2),face(92,2,2,3), v1(1),v1(2),v1(3), 1);
                %        quiver3(face(92,2,3,1),face(92,2,3,2),face(92,2,3,3), v1(1),v1(2),v1(3), 1);
                %   title('Sdir')

            end



            %                 for i=1:nfcts
            %                     bi=facets(i,3); %block index
            %                     fi=facets(i,1); %facet index
            %
            %                     if (facets(i,4)<0 && facets(i,4)>-100) %it is a floor, not a building
            %                         il=floors(bi,1);
            %                         iu=floors(bi,2);
            %                         jl=floors(bi,3);
            %                         ju=floors(bi,4);
            %                         xl=xc(il)-dx(il)/2;
            %                         xu=xc(iu)+dx(iu)/2;
            %                         yl=yc(jl)-dy(jl)/2;
            %                         yu=yc(ju)+dy(ju)/2;
            %                         zl=0;
            %                         zu=0;
            %                     elseif (facets(i,4)<=-100)  %it is a bounding wall
            %                         il=boundingwalls(bi,1);
            %                         iu=boundingwalls(bi,2);
            %                         jl=boundingwalls(bi,3);
            %                         ju=boundingwalls(bi,4);
            %                         kl=boundingwalls(bi,5)+1;
            %                         ku=boundingwalls(bi,6)+1;
            %
            %                         if (fi==2)
            %                             xl=xc(iu)+dx(iu)/2;
            %                             xu=xc(iu)+dx(iu);
            %                             yl=yc(jl)-dy(jl)/2;
            %                             yu=yc(ju)+dy(ju)/2;
            %                             zl=zc(kl)-dz(kl)/2;
            %                             zu=zc(ku)+dz(ku)/2;
            %                         elseif (fi==3)
            %                             xl=xc(il)-dx(il);
            %                             xu=xc(il)-dx(il)/2;
            %                             yl=yc(jl)-dy(jl)/2;
            %                             yu=yc(ju)+dy(ju)/2;
            %                             zl=zc(kl)-dz(kl)/2;
            %                             zu=zc(ku)+dz(ku)/2;
            %                         elseif (fi==4)
            %                             xl=xc(il)-dx(il)/2;
            %                             xu=xc(iu)+dx(iu)/2;
            %                             yl=yc(jl)-dy(jl);
            %                             yu=yc(jl)-dy(jl)/2;
            %                             zl=zc(kl)-dz(kl)/2;
            %                             zu=zc(ku)+dz(ku)/2;
            %                         else %if (fi==5)
            %                             xl=xc(il)-dx(il)/2;
            %                             xu=xc(iu)+dx(iu)/2;
            %                             yl=yc(ju)+dy(ju)/2;
            %                             yu=yc(ju)+dy(ju);
            %                             zl=zc(kl)-dz(kl)/2;
            %                             zu=zc(ku)+dz(ku)/2;
            %                         end
            %                     else
            %                         il=B(bi,1);
            %                         iu=B(bi,2);
            %                         jl=B(bi,3);
            %                         ju=B(bi,4);
            %                         kl=B(bi,5)+1;
            %                         ku=B(bi,6)+1;
            %
            %                         xl=xc(il)-dx(il)/2;
            %                         xu=xc(iu)+dx(iu)/2;
            %                         yl=yc(jl)-dy(jl)/2;
            %                         yu=yc(ju)+dy(ju)/2;
            %                         zl=zc(kl)-dz(kl)/2;
            %                         zu=zc(ku)+dz(ku)/2;
            %                     end
            %
            %                     switch facets(i, 1)
            %                         case {1} %top
            %                             x = [xl, xu, xu, xl];
            %                             y = [yl, yl, yu, yu];
            %                             z = [zu, zu, zu, zu];
            %                         case {2} %west
            %                             x = [xl, xl, xl, xl];
            %                             y = [yl, yl, yu, yu];
            %                             z = [zl, zu, zu, zl];
            %                         case {3} %east
            %                             x = [xu, xu, xu, xu];
            %                             y = [yl, yl, yu, yu];
            %                             z = [zl, zu, zu, zl];
            %                         case {4} %north
            %                             x = [xl, xl, xu, xu];
            %                             y = [yu, yu, yu, yu];
            %                             z = [zl, zu, zu, zl];
            %                         case {5} %south
            %                             x = [xl, xl, xu, xu];
            %                             y = [yl, yl, yl, yl];
            %                             z = [zl, zu, zu, zl];
            %                     end
            %                     if (i==99)
            %                         x=[9, 10, 11 , 9];
            %                         y=[8,  8, 9 , 9];
            %                     elseif (i==103)
            %                         x=[10, 11, 11, 10];
            %                         y=[7, 7, 9, 8];
            %                     end
            %
            %                 end




            c=colorbar;
            c.Title.Interpreter='latex';
            set(c, 'TickLabelInterpreter', 'latex')
            set(c,'FontSize',12*scalef)
            c.Title.String = '$S$ $[\mathrm{W}\mathrm{m}^{-2}]$';
            %c.Position(1)=c.Position(1)+0.00
            c.Position(2)=c.Position(2)+0.02;
            c.Position(4)=c.Position(4)-0.2;


            caxis([0 obj.I])
            view(3)
            axis equal
            xlim([0 obj.xh(end)])
            ylim([0 obj.yh(end)])
            % zlim([0 72])
            zlim([0 obj.zh(end)])
            set(gca,'TickLabelInterpreter','latex')
            set(gca,'FontSize',12*scalef)
            %h1=gca;
            % h1.Position=[0.08 0.1100 0.78 0.8150];
            %caxis([0 71])
            xlabel('x [-]','Interpreter','latex','FontSize',12*scalef)
            ylabel('y [-]','Interpreter','latex','FontSize',12*scalef)
            zlabel('z [-]','Interpreter','latex','FontSize',12*scalef)
            set(gcf, 'Color', 'w');
            %export_fig facetSdir.eps
        end

        
        function generateFacetAreas(obj, lcutcorners)
            % returns ndim=(dim1*dim2), area, center & corner of every facet (could be preprocessed and
            % saved)
            blocks = obj.blocks;
            facets = obj.facets;
            floors = obj.floorfacets(:, 6:end);
            boundingwalls = obj.boundingwallfacets(:, 6:end);
            cornm = obj.cornm;
            nblocks = obj.nblocks;
            nfcts = obj.nfcts;
            xb = obj.xh;
            yb = obj.yh;
            zb = obj.zh;
            delta = 0.01;
                 
            facetarea = zeros(nfcts,1);
            facetdim = zeros(nfcts,1);
            facetcorner = zeros(nfcts, 9, 3);
            facetcornervis = facetcorner;
            
            for i = 1:nfcts
                
                bi=facets(i,3); %block index
                fi=facets(i,1); %facet index
                ci=facets(i,4); %building index (-1 for roads, -101 for bounding wall)
                
                il=-1; iu=-1;jl=-1;ju=-1;kl=-1;ku=-1;it=-1;jt=-1;kt=-1;
                co=-1;
                ndim=-1;
                area=-1;
                
                if (ci<=-100) %it is a bounding wall
                    kl=boundingwalls(bi,5)+1+1;
                    ku=boundingwalls(bi,6)+1+1;
                    if (fi==2)  %east, facing west
                        jl=boundingwalls(bi,3);     %lower y index of floor facet 1
                        ju=boundingwalls(bi,4);     %upper y index of floor facet 1
                        it=boundingwalls(bi,1);
                        
                        x=xb(it+1);
                        
                        ndim=((ju-jl)+1)*((ku-kl)+1);
                        
                    elseif (fi==3)  %west, facing east
                        jl=boundingwalls(bi,3);     %lower y index of floor facet 1
                        ju=boundingwalls(bi,4);     %upper y index of floor facet 1
                        it=boundingwalls(bi,2);
                        
                        x=xb(it);
                        
                        ndim=((ju-jl)+1)*((ku-kl)+1);
                    elseif (fi==4)  %south, facing north
                        il=boundingwalls(bi,1);
                        iu=boundingwalls(bi,2);
                        jt=boundingwalls(bi,4);
                        
                        y=yb(jt);
                        ndim=((iu-il)+1)*((ku-kl)+1);
                        
                    else %if (fi==5)  %north, facing south
                        jt=boundingwalls(bi,3);
                        il=boundingwalls(bi,1);
                        iu=boundingwalls(bi,2);
                        
                        y=yb(jt+1);
                        
                        ndim=((iu-il)+1)*((ku-kl)+1);
                    end
                elseif (ci>=0) %deal with floors separately
                    if ((fi==1) || (fi==6))  %top / bottom
                        il=blocks(bi,1);     %lower x index of facet 1
                        iu=blocks(bi,2);     %upper x index of facet 1
                        jl=blocks(bi,3);     %lower y index of facet 1
                        ju=blocks(bi,4);     %upper y index of facet 1
                        
                        if fi==1
                            kt=blocks(bi,6)+1;
                            z=zb(kt+1);
                        else
                            kt=blocks(bi,5)+1;
                            z=zb(kt);
                        end
                        
                        ndim=((iu-il)+1)*((ju-jl)+1);
                        
                    elseif ((fi==2) || (fi==3))  %west / east
                        jl=blocks(bi,3);     %lower y index of facet 1
                        ju=blocks(bi,4);     %upper y index of facet 1
                        kl=blocks(bi,5)+1;     %lower z index of facet 1
                        ku=blocks(bi,6)+1;     %upper z index of facet 1
                        
                        if fi==2  %west, facing west
                            it=blocks(bi,1);
                            x=xb(it);
                        else %east, facing east
                            it=blocks(bi,2);
                            x=xb(it+1);
                        end
                        
                        ndim=((ju-jl)+1)*((ku-kl)+1);
                        
                        
                    elseif ((fi==4) || (fi==5)) %north / south
                        il=blocks(bi,1) ;    %lower x index of facet 1
                        iu=blocks(bi,2) ;    %upper x index of facet 1
                        kl=blocks(bi,5)+1 ;    %lower y index of facet 1
                        ku=blocks(bi,6)+1 ;    %upper y index of facet 1
                        
                        if fi==4 %north, facing north
                            jt=blocks(bi,4);
                            y=yb(jt+1);
                        else  %south, facing south
                            jt=blocks(bi,3);
                            y=yb(jt);
                        end
                        
                        ndim=((iu-il)+1)*((ku-kl)+1);
                    end
                end
                
                %return xyz-coordinates of center and 4 corners clockwise from bottom left
                %move coordinates of corners very slightly to the interior of the facet (by
                %delta)
                if (ci>=0 || ci<=-100) %it's not a floor -> it's a building or bounding wall
                    
                    if ((fi==1) || (fi==6))
                        co=[(xb(iu+1)+xb(il))*.5, (yb(ju+1)+yb(jl))*.5  , z; ... %center
                            xb(il)+delta        , yb(jl)+delta          , z; ... %4 corners
                            xb(il)+delta        , yb(ju+1)-delta        , z; ...
                            xb(iu+1)-delta      , yb(ju+1)-delta        , z; ...
                            xb(iu+1)-delta      , yb(jl)+delta          , z; ...
                            xb(il)        , yb(jl)          , z; ... %4 corners
                            xb(il)        , yb(ju+1)        , z; ...
                            xb(iu+1)      , yb(ju+1)        , z; ...
                            xb(iu+1)      , yb(jl)          , z; ...
                            ];
                        area=(xb(iu+1)-xb(il))*(yb(ju+1)-yb(jl));
                        
                    elseif ((fi==2) || (fi==3))
                        
                        co=[ x  , (yb(ju+1)+yb(jl))*.5  , (zb(ku+1)+zb(kl))*.5; ... %center
                            x  , yb(jl)+delta          , zb(kl)+delta        ; ... %4 corners
                            x  , yb(jl)+delta          , zb(ku+1)-delta      ; ...
                            x  , yb(ju+1)-delta        , zb(ku+1)-delta      ; ...
                            x  , yb(ju+1)-delta        , zb(kl)+delta        ; ...
                            x  , yb(jl)          , zb(kl)        ; ... %4 corners
                            x  , yb(jl)          , zb(ku+1)      ; ...
                            x  , yb(ju+1)        , zb(ku+1)      ; ...
                            x  , yb(ju+1)        , zb(kl)        ; ...
                            ];
                        area=(zb(ku+1)-zb(kl))*(yb(ju+1)-yb(jl));
                        
                    elseif ((fi==4) || (fi==5))
                        
                        
                        co=[(xb(iu+1)+xb(il))*.5     , y  , (zb(ku+1)+zb(kl))*.5  ; ... %center
                            xb(il)+delta             , y  , zb(kl)+delta          ; ... %4 corners
                            xb(il)+delta             , y  , zb(ku+1)-delta        ; ...
                            xb(iu+1)-delta           , y  , zb(ku+1)-delta        ; ...
                            xb(iu+1)-delta           , y  , zb(kl)+delta          ; ...
                            xb(il)             , y  , zb(kl)          ; ... %4 corners
                            xb(il)             , y  , zb(ku+1)        ; ...
                            xb(iu+1)           , y  , zb(ku+1)        ; ...
                            xb(iu+1)           , y  , zb(kl)          ; ...
                            ];
                        area=(xb(iu+1)-xb(il))*(zb(ku+1)-zb(kl));
                    end
                
                    covis = co;
                else  % it is a floor, not a building
                    bif = bi - nblocks; % index in floor list
                    il=floors(bif,1);     %lower x index of floor facet 1
                    iu=floors(bif,2);     %upper x index of floor facet 1
                    jl=floors(bif,3);     %lower y index of floor facet 1
                    ju=floors(bif,4);     %upper y index of floor facet 1
                    %z=1;
                    z=zb(2);
                    ndim=((iu-il)+1)*((ju-jl)+1);
                    co=[(xb(iu+1)+xb(il))*.5, (yb(ju+1)+yb(jl))*.5  , z; ... %center
                        xb(il)+delta        , yb(jl)+delta          , z; ... %4 corners
                        xb(il)+delta        , yb(ju+1)-delta        , z; ...
                        xb(iu+1)-delta      , yb(ju+1)-delta        , z; ...
                        xb(iu+1)-delta      , yb(jl)+delta          , z; ...
                        xb(il)        , yb(jl)          , z; ... %4 corners
                        xb(il)        , yb(ju+1)        , z; ...
                        xb(iu+1)      , yb(ju+1)        , z; ...
                        xb(iu+1)      , yb(jl)          , z; ...
                        ];
                    area=(xb(iu+1)-xb(il))*(yb(ju+1)-yb(jl));
                    %if it is 1*1 nothing has to be done since it matches both walls
                    if lcutcorners
                        covis = co;
                        if ju-jl>0  %floor facet along a west or east wall
                            %test if it is a wall-wall-floor corner, if so, the floor facet needs
                            %to be cut diagonally, see "createfloors.m"
                            if cornm(iu,ju)==10  %south/west corner
                                %JUST REMOVE IT ONE dx OR dy FROM THE OTHER BLOCK
                                co(2+1,2)=yb(ju)+delta;
                                co(6+1,2)=yb(ju);
                                area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
                            elseif cornm(iu,ju)==15  %south/east corner
                                co(3+1,2)=yb(ju)+delta;
                                co(7+1,2)=yb(ju);
                                area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
                            end
                            if cornm(iu,jl)==8  %north/west corner
                                co(1+1,2)=yb(jl+1)-delta;
                                co(5+1,2)=yb(jl+1);
                                area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
                            elseif cornm(iu,jl)==12  %north/east corner
                                co(4+1,2)=yb(jl+1)-delta;
                                co(8+1,2)=yb(jl+1);
                                area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
                            end
                            
                        elseif iu-il>0 %floor facet along along a south or north wall
                            if cornm(iu,ju)==10  %south/west corner
                                co(4+1,1)=xb(iu)+delta;
                                co(8+1,1)=xb(iu);
                                area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
                            elseif cornm(iu,ju)==8  %north/west corner
                                co(3+1,1)=xb(iu)+delta;
                                co(7+1,1)=xb(iu);
                                area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
                            end
                            if cornm(il,ju)==12  %north/east corner
                                co(2+1,1)=xb(il+1)-delta ;
                                co(6+1,1)=xb(il+1) ;
                                area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
                            elseif cornm(il,ju)==15  %south/east corner
                                co(1+1,1)=xb(il+1)-delta;
                                co(5+1,1)=xb(il+1);
                                area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
                            end
                        end
                    else
                        covis = co;
                    end
                end
                
                facetdim(i) = ndim;
                facetarea(i) = area;
                facetcorner(i,:,:) = co;
                facetcornervis(i,:,:) = covis;
                
            end
            preprocessing.addvar(obj, 'facetdim', facetdim)
            preprocessing.addvar(obj, 'facetarea', facetarea)
            preprocessing.addvar(obj, 'facetcorner', facetcorner)
            preprocessing.addvar(obj, 'facetcornervis', facetcornervis)
        end
        
        function generateUniqueFacetCorners(obj)
            
            facetcorner = obj.facetcornervis;
            nfcts = obj.nfcts;
            
            surfaces = zeros(nfcts, 4);
            
            corners = squeeze(facetcorner(1,6:9,:));
            corners_unique = corners([1,4,3,2],:);
            surfaces(1,1:4) = [1,4,3,2];
            
            for i=2:nfcts
                %i
                corners = squeeze(facetcorner(i,6:9,:));
                
                for n = [1,4,3,2] % Ivo's corners defined clockwise but anticlockwise in view3d
                    pos = find(ismember(corners_unique, corners(n,:),'rows'));
                    if isempty(pos)
                        corners_unique = [corners_unique; corners(n,:)];
                        surfaces(i,n) = size(corners_unique,1);
                    else
                        surfaces(i,n) = pos;
                    end
                end
            end
            
            preprocessing.addvar(obj, 'facetcornerunique', corners_unique)
            preprocessing.addvar(obj, 'surfaces', surfaces)
            %ncorners = size(corners_unique,1);
        end
        
        function plot_facets_new(obj)
            surfaces = obj.surfaces;
            corners_unique = obj.facetcornerunique;
            
            clr = [.75 .75 .75];
            %clr_triplet = 1 ./ (1:obj.nfcts)';
            %clr = [clr_triplet; clr_triplet; clr_triplet];
            %clr_index = (1:obj.nfcts)';
            %clr = obj.Sdir;
            %clr_k = [1, 1, 1]
            
            figure
            light('Position', [-0.75 -0.35 1], 'Style', 'infinite');
            p = patch('Faces',surfaces,'Vertices',corners_unique,'FaceVertexCData', clr, 'FaceLighting','flat','AmbientStrength', 0.3, 'SpecularStrength', 0.6, 'DiffuseStrength', 0.7);%, 'FaceAlpha', 0.5);
%             patch('Faces',surfaces,'Vertices',corners_unique, 'FaceVertexCData',zeros(size(corners_unique,1),3),'EdgeColor','flat','FaceColor', 'flat');% 'FaceLighting','flat','AmbientStrength', 0.3, 'SpecularStrength', 0.6, 'DiffuseStrength', 0.7);
%             clear S
%             S.Vertices = corners_unique;
%             S.Faces = surfaces;
%             %S.FaceVertexCData = 0.5 * ones(size(surfaces,1),1);
%             nsurfaces = size(surfaces,1);
%             S.FaceVertexCData = randi(nsurfaces,nsurfaces,1);%(1:nsurfaces)';
%             S.FaceColor = 'flat';
%             S.EdgeColor = 'red';
%             S.LineWidth = 10;
%             patch(S)
            
            shading flat
            view(3)
            %view(180 - obj.solaz, 90-obj.solarzenith)
            axis equal
            grid on
            xlabel('$x(\mathrm{m})$','interpreter','latex')
            ylabel('$y(\mathrm{m})$','interpreter','latex')
            zlabel('$z(\mathrm{m})$','interpreter','latex')
            %colorbar
        end
        
        function vfc(obj)
            xh = obj.xh;
            yh = obj.yh;
            zh = obj.zh;
            blocks = obj.blocks; nblocks = obj.nblocks;
            facets = obj.facets; nfcts  = obj.nfcts;
            floorfacets = obj.floorfacets;  floors = floorfacets(:, 6:end); %nfloorfacets = obj.nfloorfacets;
            boundingwallfacets = obj.boundingwallfacets; boundingwalls = boundingwallfacets(:, 6:end); nboundingwallfacets = obj.nboundingwallfacets; %nbw = nboundingwallfacets;
            blocks_phys = obj.blocks_phys;
            cornm = obj.cornm;
            centerweight = obj.centerweight; cornerweight = obj.cornerweight; delta = 0.01;
            
            %% View factors
            % Calculate fraction of facets seeing each other (not blocked)
            % Calculate orientation, no overlap on orthogonal planes
            % calculate viewfactors
            vf = zeros(nfcts, nfcts);
            pf1sf2 = zeros(nfcts, nfcts);
            A = zeros(nfcts, 1);
            %tim = zeros(obj.nfcts - 1, 1);
            %%
            %do dummy parallelisation by starting multiple matlab and only iterate over
            %a fraction and then write to file and merge files
            for i = 1:(nfcts - 1)
                i
                bi = facets(i, 3); %block index
                fi = facets(i, 1); %facet index
                ci = facets(i, 4); %building index (-1 for roads, -99 for bounding wall)
                %[ndima, areaa, coa] = preprocessing.detsub(obj, i);
                
                [ ndima, areaa, coa] = preprocessing.detsub(i,facets,blocks,floors,boundingwalls,cornm,nblocks,xh,yh,zh,delta);
%                 ndima = obj.facetdim(i);
%                 areaa = obj.facetarea(i);
%                 coa   = obj.facetcorners(i,:);
                for j = (i + 1):nfcts
                    bi2 = facets(j, 3); %block index
                    fi2 = facets(j, 1); %facet index
                    ci2 = facets(j, 4); %building index (-1 for roads, -99 for bounding wall)
                    if ((fi2 == fi) || ((ci2 * bi2) == (ci * bi)))  %facets looking in same direction or on same block&building (ci*bi is unique), bi alone could also be a floor or a bounding wall
                        continue
                    end
                    
                    %[ndimb, areab, cob] = preprocessing.detsub(obj, j);
                    [ ndimb, areab, cob] = preprocessing.detsub(j,facets,blocks,floors,boundingwalls,cornm,nblocks,xh,yh,zh,delta);
                    %[prblckdij] = preprocessing.prblckd(obj, i, j, coa, ndima, false, -999, cob, ndimb);
                    [prblckdij] = preprocessing.prblckd(i,j,coa,ndima,false,-999,cob,ndimb,facets,centerweight,cornerweight,nblocks,nboundingwallfacets,blocks_phys);
                    c = 1 - prblckdij;
                    pf1sf2(i, j) = c;
                    pf1sf2(j, i) = c;
                    A(i) = areaa;
                    A(j) = areab;
                    
                    if prblckdij < 0.98 %view is not blocked, vf calclation
                        if (prblckdij >= (2 * cornerweight - eps)) %at least two corners have to be blocked, otherwise it cannot possibly be blocked by itself
                            %slice up
                            %pf1sf2u also has to include check for center, in case center is
                            %blocked
                            [coaa, cobb, pf1sf2u] = preprocessing.slice(fi, fi2, coa, cob, cornerweight, centerweight);
                            c = c + pf1sf2u;
                            pf1sf2(i, j) = c;
                            pf1sf2(j, i) = c;
                        else
                            coaa = coa(6:9, :);
                            cobb = cob(6:9, :);
                        end
                        
                        %determine if they form a corner/have common edge
                        %determine the vertices pair that forms that corner
                        vcorner = [0, 0, 1]; %corner vertices for i and j & sign of correction (as always, clockwise from bottom left)
                        glpo = 6; %order of gauss-legendre polynomial
                        glpop = 15; %more accurate, order of gauss-legendre polynomial
                        if (fi == 1) %fi is horizontal, test if corner with a vertical wall
                            if  (fi2 == 2)
                                if all(coaa(4, :) == cobb(1, :)) && all(coaa(3, :) == cobb(4, :)) %they share one edge completely
                                    vcorner = [3, 4, -1];
                                    glpo = glpop;
                                    % elseif (coaa(4,1)==cobb(1,1)) && (coaa(4,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %     glpo=18;  %NEEDS TO BE AN EVEN NUMBER!!!
                                end
                            elseif (fi2 == 3)
                                if all(coaa(1, :) == cobb(1, :)) && all(coaa(2, :) == cobb(4, :)) %they share one edge completely
                                    vcorner = [1, 4, 1];
                                    glpo = glpop;
                                    % elseif (coaa(1,1)==cobb(1,1)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %     glpo=18;
                                end
                            elseif (fi2==4)
                                if all(coaa(1,:) == cobb(1,:)) && all(coaa(4,:)==cobb(4,:)) %they share one edge completely
                                    corner = [4, 4, -1];
                                    glpo = glpop;
                                    %  elseif (coaa(1,2)==cobb(1,2)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %      glpo=18;
                                end
                            elseif (fi2 == 5)
                                if all(coaa(2, :) == cobb(1,:)) && all(coaa(3, :) == cobb(4, :)) %they share one edge completely
                                    vcorner = [2, 4, 1];
                                    glpo = glpop;
                                    %  elseif (coaa(2,2)==cobb(1,2)) && (coaa(2,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %      glpo=18;
                                end
                            end
                        elseif (fi == 2)
                            if  (fi2 == 4)
                                if all(coaa(1, :) == cobb(4, :)) && all(coaa(2, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [1, 3, 1];
                                    glpo = glpop;
                                elseif (coaa(1, 1) == cobb(4, 1)) && (coaa(1, 2) ==cobb(4, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 5)
                                if all(coaa(4, :) == cobb(4, :)) && all(coaa(3, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [3, 3, -1];
                                    glpo = glpop;
                                elseif (coaa(4,1)==cobb(4,1)) && (coaa(4,2)==cobb(4,2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 1)
                                if all(coaa(1, :) == cobb(4, :)) && all(coaa(4, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [4, 3, -1];
                                    glpo = glpop;
                                    %                    elseif (coaa(1,1)==cobb(4,1)) && (coaa(1,3)==cobb(4,3)) %they potentially share parts of an edge
                                    %                        glpo=18;
                                end
                            end
                        elseif (fi == 3)
                            if  (fi2 == 4)
                                if all(coaa(1, :) == cobb(1, :)) && all(coaa(2, :) == cobb(2,:)) %they share one edge completely
                                    vcorner = [1, 1, -1];
                                    glpo = glpop;
                                elseif (coaa(1,1)==cobb(1,1)) && (coaa(1,2)==cobb(1,2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 5)
                                if all(coaa(4, :) == cobb(1, :)) && all(coaa(3, :) == cobb(2, :)) %they share one edge completely
                                    vcorner = [3, 1, 1];
                                    glpo = glpop;
                                elseif (coaa(4, 1) == cobb(1, 1)) && (coaa(4, 2) == cobb(1, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 1)
                                if all(coaa(1, :) == cobb(1, :)) && all(coaa(4, :) == cobb(2, :)) %they share one edge completely
                                    vcorner = [4, 1, 1];
                                    glpo = glpop;
                                    %  elseif (coaa(1,1)==cobb(1,1)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %      glpo=18;
                                end
                            end
                        elseif (fi == 4)
                            if  (fi2 == 2)
                                if all(coaa(4, :) == cobb(1, :)) && all(coaa(3, :) == cobb(2, :)) %they share one edge completely
                                    vcorner = [3, 1, 1];
                                    glpo = glpop;
                                elseif (coaa(4, 1) == cobb(1, 1)) && (coaa(4, 2) == cobb(1, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 3)
                                if all(coaa(1, :) == cobb(1, :)) && all(coaa(2, :) == cobb(2, :)) %they share one edge completely
                                    vcorner = [1, 1, -1];
                                    glpo = glpop;
                                elseif (coaa(1, 1) == cobb(1, 1)) && (coaa(1, 2) == cobb(1, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 1)
                                if all(coaa(1, :) == cobb(1, :)) && all(coaa(4, :) == cobb(4, :)) %they share one edge completely
                                    vcorner = [4, 4, -1];
                                    glpo = glpop;
                                    %    elseif (coaa(1,2)==cobb(1,2)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %       glpo=18;
                                end
                            end
                        elseif (fi == 5)
                            if  (fi2 == 2)
                                if all(coaa(4, :) == cobb(4, :)) && all(coaa(3, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [3, 3, -1];
                                    glpo = glpop;
                                elseif (coaa(4, 1) == cobb(4, 1)) && (coaa(4, 2) == cobb(4, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 3)
                                if all(coaa(1, :) == cobb(4, :)) && all(coaa(2, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [1, 3, 1];
                                    glpo = glpop;
                                elseif (coaa(1, 1) == cobb(4, 1)) && (coaa(1, 2) == cobb(4, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 1)
                                if all(coaa(1, :) == cobb(2, :)) && all(coaa(4, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [4, 2, 1];
                                    glpo = glpop;
                                    % elseif (coaa(1,2)==cobb(2,2)) && (coaa(1,3)==cobb(2,3)) %they potentially share parts of an edge
                                    %     glpo=18
                                end
                            end
                        end
                        
                        %calculate viewfactor
                        [F12, F21] = preprocessing.ViewFactor(coaa, cobb, areaa, areab, glpo, vcorner);
                        %increase glpo if view factors are big
                        if F12 > 0.5 || F21 > 0.5
                            [F12, F21] = preprocessing.ViewFactor(coaa, cobb, areaa, areab, glpo + 10, vcorner);
                        end
                        
                        
                        % round to 1promille accuracy? i.e. cut smaller 1promille completely
                        % this violates reciprocity though
                        
                        vf(i,j) = F12 * c;
                        vf(j,i) = F21 * c;
                        %   vf(i,j)=floor(F12*c*100)/100;
                        %  vf(j,i)=floor(F21*c*100)/100;
                    end
                end
            end
            
            vf = single(vf);
            %pf1sf2o=single(pf1sf2);
            blub = sum(vf, 2);
            lblub = find(blub > 1);
            if ~isempty(lblub)
                %[maxvf, indexmaxvf] = max(blub);
                for i = 1:length(lblub)
                    vf(lblub(i), :) = vf(lblub(i),:)/blub(lblub(i));
                end
            end
            preprocessing.addvar(obj, 'vf', vf)
            preprocessing.addvar(obj, 'svf', max(1 - sum(vf, 2), 0));
            preprocessing.addvar(obj, 'facetarea', A);
        end
        
        function vfcAccelerated(obj)
            %xh = obj.xh;
            %yh = obj.yh;
            %zh = obj.zh;
            %blocks = obj.blocks; 
            nblocks = obj.nblocks;
            facets = obj.facets; nfcts = obj.nfcts;
            %floorfacets = obj.floorfacets;  %floors = floorfacets(:, 6:end); %nfloorfacets = obj.nfloorfacets;
            %boundingwallfacets = obj.boundingwallfacets; %boundingwalls = boundingwallfacets(:, 6:end); nboundingwallfacets = obj.nboundingwallfacets; %nbw = nboundingwallfacets;
            blocks_phys = obj.blocks_phys;
            %cornm = obj.cornm;
            centerweight = obj.centerweight; cornerweight = obj.cornerweight; %delta = 0.01;
            facetarea = obj.facetarea;
            facetdim = obj.facetdim;
            facetcorner = obj.facetcorner;
            
            min_blocks_x = blocks_phys(:,1); min_blocks_y = blocks_phys(:,3); min_blocks_z = blocks_phys(:,5);
            max_blocks_x = blocks_phys(:,2); max_blocks_y = blocks_phys(:,4); max_blocks_z = blocks_phys(:,6);                  
            
            %% View factors
            % Calculate fraction of facets seeing each other (not blocked)
            % Calculate orientation, no overlap on orthogonal planes
            % calculate viewfactors
            vf = zeros(nfcts, nfcts);
            %pf1sf2 = zeros(nfcts, nfcts);
            %A = zeros(nfcts, 1);
            %tim = zeros(obj.nfcts - 1, 1);
            %%
            %do dummy parallelisation by starting multiple matlab and only iterate over
            %a fraction and then write to file and merge files
            
            % Only calculate view factors for building facets
            sel = find(facets(:,2) == 1);
            size(sel)
            
            
            for i = 1:(nfcts - 1)
            %for i = sel'
            %for i = 1;
            %for i = (nfcts+1)/2:(nfcts+1)/2+5
                i
                bi = facets(i, 3); %block
                fi = facets(i, 1); %orientation
                ci = facets(i, 4); %building index (-1 for roads, -99 for bounding wall)
                %[ndima, areaa, coa] = preprocessing.detsub(obj, i);

%                 xli = facets(i,6); xui = facets(i,7); 
%                 yli = facets(i,8); yui = facets(i,9); 
%                 zli = facets(i,10); zui = facets(i,11);
                
                %[ ndima, areaa, coa] = preprocessing.detsub(i,facets,blocks,floors,boundingwalls,cornm,nblocks,xh,yh,zh,delta);
                ndima = facetdim(i);
                areaa = facetarea(i);
                coa = squeeze(facetcorner(i,:,:));
                %% Geometry tests
                facids = (i+1):nfcts;
                facids_o = (i+1):nfcts;
                facids_i = (i+1):nfcts;
                facids_og = (i+1):nfcts;
                % Elimination due to orientation of i
                
                % Must be different orientation;
                facids_o = facids_o(facets(facids_o,1) ~= facets(i,1));
                facids_o = facids_o(facets(facids_o,3) .* facets(facids_o,4) ~= facets(i,3) * facets(i,4));
                
                facids = facids(facets(facids,1) ~= facets(i,1));
                facids = facids(facets(facids,3) .* facets(facids,4) ~= facets(i,3) * facets(i,4));
                
                switch fi
                    case 1 % Top: maximum zs must be greater than z.
                        facids_i = facids_i(facets(facids_i,11) > facets(i,10));
                        facids = facids(facets(facids,11) > facets(i,10));
                    case 4 % North: maximum ys must be greater than y
                        facids_i = facids_i(facets(facids_i,9) > facets(i,8));
                        facids = facids(facets(facids,9) > facets(i,8));
                    case 3 % East: maximum xs must be greater than x
                        facids_i = facids_i(facets(facids_i,7) > facets(i,6));
                        facids = facids(facets(facids,7) > facets(i,6));
                    case 5 % South: minimum ys must be less than y
                        facids_i = facids_i(facets(facids_i,8) < facets(i,9));
                        facids = facids(facets(facids,8) < facets(i,9));
                    case 2 % West: minimum xs must be less than x
                        facids_i =  facids_i(facets(facids_i,6) < facets(i,7));
                        facids =  facids(facets(facids,6) < facets(i,7));
                    case 6 % Bottom: minimum zs must be less than z
                        facids_i = facids_i(facets(facids_i,10) < facets(i,11));
                        facids = facids(facets(facids,10) < facets(i,11));
                end
                
                facids_j = [];
                facids_jog = [];
                % Elimination due to orientation of j
                if fi ~= 1 % Top: maximum z must be greater than zs;
                     %fT = facids(facets(facids,1) == 1);
                     %facids_j = [facids_j, facids(facets(i,11) > facets(fT, 10))]
                     facids_j = [facids_j, intersect(facids(facets(facids,1) == 1), facids(facets(i,11) > facets(facids, 10)))];
                     facids_jog = [facids_j, intersect(facids_og(facets(facids_og,1) == 1), facids_og(facets(i,11) > facets(facids_og, 10)))];
                end
                if fi ~= 4 % North:
                    %fN = facids(facets(facids,1) == 4);
                    %facids_j  = [facids_j, facids(facets(i,9) > facets(fN,8))]
                    facids_j  = [facids_j, intersect(facids(facets(facids,1) == 4), facids(facets(i,9) > facets(facids,8)))];                    
                    facids_jog  = [facids_jog, intersect(facids_og(facets(facids_og,1) == 4), facids_og(facets(i,9) > facets(facids_og,8)))];
                end
                if fi ~= 3 % East:
                    %fE = facids(facets(facids,1) == 3);
                    %facids_j = [facids_j, facids(facets(i,7) > facets(fE,6))]
                    facids_j = [facids_j, intersect(facids(facets(facids,1) == 3), facids(facets(i,7) > facets(facids,6)))];
                    facids_jog = [facids_jog, intersect(facids_og(facets(facids_og,1) == 3), facids_og(facets(i,7) > facets(facids_og,6)))];
                end
                if fi ~= 5 % South
                    %fS = facids(facets(facids,1) == 5);
                    %facids_j = [facids_j, facids(facets(i,8) < facets(fS,9))];
                    facids_j = [facids_j, intersect(facids(facets(facids,1) == 5), facids(facets(i,8) < facets(facids,9)))];
                    facids_jog = [facids_jog, intersect(facids_og(facets(facids_og,1) == 5), facids_og(facets(i,8) < facets(facids_og,9)))];
                end
                if fi ~= 2 % West
                    %fW = facids(facets(facids,1) == 2);
                    %facids_j = [facids_j, facids(facets(i,6) < facets(fW,7))];
                    facids_j = [facids_j, intersect(facids(facets(facids,1) == 2), facids(facets(i,6) < facets(facids,7)))];
                    facids_jog = [facids_jog, intersect(facids_og(facets(facids_og,1) == 2), facids_og(facets(i,6) < facets(facids_og,7)))];
                end
                if fi ~= 6 % Bottom
                     facids_j = [facids_j, intersect(facids(facets(facids,1) == 6), facids(facets(i,10) < facets(facids, 11)))];
                     facids_jog = [facids_j, intersect(facids_og(facets(facids_og,1) == 6), facids_og(facets(i,10) < facets(facids_og, 11)))];
                end
                
                facids = intersect(facids_j, facids);
                facids_compare = intersect(intersect(facids_o, facids_i), facids_jog);
                size(facids)
                size(facids_compare)
                
                if isempty(facids)
                     disp('No calculation')
                     continue
                end
                
                
                %count=0;
                for j = facids
                %for j = (i + 1):nfcts
                    %j
%                     if ~ismember(j,facids_o)
%                         disp('Facets cant see each other according to vectorial method because same orientation/block')
%                     end
%                     if ~ismember(j,facids_i)
%                         disp('Facets cant see each other according to vectorial method due to i')
%                     end
%                     if ~ismember(j,facids_j)
%                         disp('Facets cant see each other according to vectorial method due to j')
%                     end
%                    bi2 = facets(j, 3); %block
                    fi2 = facets(j, 1); %orientation
%                    ci2 = facets(j, 4); %building index (-1 for roads, -99 for bounding wall)
%                     if ((fi2 == fi) || ((ci2 * bi2) == (ci * bi)))  %facets looking in same direction or on same block&building (ci*bi is unique), bi alone could also be a floor or a bounding wall
%                         %disp('Facets cant see each other according to pairwise method because same orientation/block')
%                         continue
%                     end
%                     
%                     % More tests
%                     xlj = facets(j,6); xuj = facets(j,7); ylj = facets(j,8); yuj = facets(j,9); zlj = facets(j,10); zuj = facets(j,11);
%                     switch fi
%                         case 1 % Top
%                             if (zuj <= zli)
%                                 continue
%                                 %disp('Facets cant see each other according to pairwise method due to i')
%                             end
%                         case 4 % North
%                             if (yuj <= yli)
%                                 continue
%                                 %disp('Facets cant see each other according to pairwise method due to i')
%                             end
%                         case 3 % East
%                             if (xuj <= xli)
%                                 continue
%                                 %disp('Facets cant see each other according to pairwise method due to i')
%                             end
%                         case 5 % South
%                             if (ylj >= yui)
%                                 continue
%                                 %disp('Facets cant see each other according to pairwise method due to i')
%                             end
%                         case 2 % West
%                             if (xlj >= xui)
%                                 continue
%                                 %disp('Facets cant see each other according to pairwise method due to i')
%                             end
%                     end
%                     
%                     switch fi2
%                         case 1 % Top
%                             if (zui <= zlj)
%                                 continue
%                                 %disp('Facets cant see each other according to pairwise method due to j')
%                             end
%                         case 4 % North
%                             if (yui <= ylj)
%                                 continue
%                                 %disp('Facets cant see each other according to pairwise method due to j')
%                             end
%                         case 3 % East
%                             if (xui <= xlj)
%                                 continue
%                                 %disp('Facets cant see each other according to pairwise method due to j')
%                             end
%                         case 5 % South
%                             if (yli >= yuj)
%                                 continue
%                                 %disp('Facets cant see each other according to pairwise method due to j')
%                             end
%                         case 2 % West
%                             if (xli >= xuj)
%                                 continue
%                                 %disp('Facets cant see each other according to pairwise method due to j')
%                             end
%                     end
% 
%                     count = count+1;

                    %% Calculate view factor
                    %[ ndimb, areab, cob] = preprocessing.detsub(j,facets,blocks,floors,boundingwalls,cornm,nblocks,xh,yh,zh,delta);

                    ndimb = facetdim(j);
                    areab = facetarea(j);
                    cob = squeeze(facetcorner(j,:,:));
                    %[prblckdij] = preprocessing.prblckd(i,j,coa,ndima,false,-999,cob,ndimb,facets,centerweight,cornerweight,nblocks,nboundingwallfacets,blocks_phys);
                    %prblckdij = percentageBlocked(i,j,coa,ndima,false,-999,cob,ndimb,facets,centerweight,cornerweight,nblocks,nboundingwallfacets,blocks_phys);
                    
                    testcrit = 0;
                    prblckdij = 0;
                    flag = 0;
                    
                    if (ndima+ndimb)>(2*testcrit)  %also test corner, otherwise only test center
                        %disp('test corners')
                        %disp(['% blocked so far ' num2str(prcntgblckd)])
                        %find orientation to know which corners to compare
                        ordera=[1 2 3 4 5];
                        cas=fi*10+fi2;
                        %disp('case')
                        switch cas
                            case 12, orderb=[1 3 4 5 2];
                            case 13, orderb=[1 2 5 4 3];
                            case 15, orderb=[1 3 2 5 4];
                            case 21, orderb=[1 5 2 3 4];
                            case 24, orderb=[1 5 4 3 2];
                            case 26, orderb=[1 2 5 4 3];
                            case 31, orderb=[1 2 5 4 3];
                            case 35, orderb=[1 5 4 3 2];
                            case 36, orderb=[1 5 2 3 4];
                            case 42, orderb=[1 5 4 3 2];
                            case 46, orderb=[1 3 2 5 4];
                            case 51, orderb=[1 3 2 5 4];
                            case 53, orderb=[1 5 4 3 2];
                            case 62, orderb=[1 2 5 4 3];
                            case 63, orderb=[1 3 4 5 2];
                            case 64, orderb=[1 3 2 5 4];
                            otherwise %14, 41, 23, 32, 25, 52, 34, 43, 45, 54
                                orderb=[1 2 3 4 5];
                        end
                        
                        for k=1:5
                            test=1;
                            oa=ordera(k);
                            ob=orderb(k);
                            facetpoint=coa(oa,:);
                            facetpoint2=cob(ob,:);
                            v=[facetpoint2(1)-facetpoint(1),facetpoint2(2)-facetpoint(2),facetpoint2(3)-facetpoint(3)];
                            dint=norm(v);
                            v1=v/dint;
                            % Add test for bottom facets
                            if ((fi==1) && (v(3)<=0)) %horizontal facets can't see facet that is lower.
                                test=0;
                            elseif ((fi==2) && (v(1)>=0)) %west facets cant see anything to the east
                                test=0;
                            elseif ((fi==3) && (v(1)<=0)) %east facets cant see anything to the west
                                test=0;
                            elseif ((fi==4) && (v(2)<=0)) %north facets cant see anything to the south
                                test=0;
                            elseif ((fi==5) && (v(2)>=0)) %south facets cant see anything to the north
                                test=0;
                            elseif ((fi~=1) && (fi2==1) && (v(3)>0)) %vertical walls cant see any horizontal facet that is higher
                                test=0;
                                % the following cases should be blocked by a block anyway, but
                                % this is faster
                            elseif ((fi==4) || (fi==5)) && (fi2==3) && (v(1)>0) %north/south can't see an east facet if it is east of it
                                test=0;
                            elseif ((fi==4) || (fi==5)) && (fi2==2) && (v(1)<0) %north/south can't see a west facet if it is west of it
                                test=0;
                            elseif ((fi==2) || (fi==3)) && (fi2==4) && (v(2)>0) %west/east can't see a north facet if it is north of it
                                test=0;
                            elseif ((fi==2) || (fi==3)) && (fi2==5) && (v(2)<0) %west/east can't see a south facet if it is south of it
                                test=0;
                            end
                            
                            if test
%                                 min_blocks = [blocks_phys(:,1), blocks_phys(:,3), blocks_phys(:,5)];
%                                 max_blocks = [blocks_phys(:,2), blocks_phys(:,4), blocks_phys(:,6)];
                                
                                or = facetpoint;
                                d = v1;
                                
%                                 [tmin_vec,flag_vec] = rayBoxVec(min_blocks(:,1)', min_blocks(:,2)', min_blocks(:,3)', ...
%                                                                               max_blocks(:,1)', max_blocks(:,2)', max_blocks(:,3)', ...
%                                                                               or(1), or(2), or(3)  , ...
%                                                                               d(1), d(2), d(3));
                                                                          
                                [tmin_vec,flag_vec] = rayBoxVec(min_blocks_x, min_blocks_y, min_blocks_z, max_blocks_x, max_blocks_y, max_blocks_z, or(1), or(2), or(3), d(1), d(2), d(3));
                                
%                                 gmin_blocks = gpuArray(min_blocks); gmax_blocks = gpuArray(max_blocks);
%                                 [gtmin,gflag] = arrayfun(@rayBoxGPU, gmin_blocks(:,1)', gmin_blocks(:,2)', gmin_blocks(:,3)', ...
%                                                                      gmax_blocks(:,1)', gmax_blocks(:,2)', gmax_blocks(:,3)', ...
%                                                                      or(1), or(2), or(3)  , ...
%                                                                      d(1), d(2), d(3));
%                                 tmin_vec = gather(gtmin);
%                                 flag_vec = gather(gflag);
                                
                                flag_vec(bi) = 0;
                                
                                flag_vec(tmin_vec <= eps) = 0;
                                flag_vec(tmin_vec >= dint - eps) = 0;
                                                                
                                if any(flag_vec)
                                    if k==1
                                        prblckdij=centerweight;
                                    else
                                        %disp(['Vector: blocked by ' num2str(find(flag_vec,1,'first'))])
                                        prblckdij=prblckdij+cornerweight;
                                    end
                                end
                                
%                                 for n=1:nblocks %check if any block intersects
%                                     if ((n==bi) && (ci>0))  %don't intersect with itself
%                                         continue
%                                     end
%                                     
%                                     %[flag,dints] = preprocessing.rbi(facetpoint, v1, bl(n,:));
%                                     [flag,dints] = rayBoxIntersection(facetpoint, v1, blocks_phys(n,:));
%                                     
%                                     if dints<=(0+eps)  %intersection is downstream
%                                         flag=0;
%                                     elseif dints>=(dint-eps) %block is behind target subfacet
%                                         flag=0;
%                                     end
%                                     if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
%                                         if k==1
%                                             prblckdij=centerweight;
%                                         else
%                                             %disp(['Loop: blocked by ' num2str(n)])
%                                             prblckdij=prblckdij+cornerweight;
%                                         end
%                                         break %get out of inner for loop
%                                     end
%                                 end
% 
                            else %facets can't possibly see each other
                                if k==1
                                    prblckdij=centerweight;
                                else
                                    prblckdij=prblckdij+cornerweight;
                                end
                            end
                        end
                        
                    else %small area, checking only center
                        % %disp('small area, checkig only center')
                        test=1;
                        facetpoint=coa(1,:);
                        facetpoint2=cob(1,:);
                        v=[facetpoint2(1)-facetpoint(1),facetpoint2(2)-facetpoint(2),facetpoint2(3)-facetpoint(3)];
                        dint=norm(v);
                        v1=v/dint;
                        
                        %abs(v)
                        %%disp(['v1: ' num2str(v1)])
                        %test if facets can potentially see each other
                        if ((fi==1) && (v(3)<=0)) %horizontal facets can't see facets whos center is lower
                            test=0;
                        elseif ((fi==2) && (v(1)>=0)) %west facets cant see anything to the east
                            test=0;
                        elseif ((fi==3) && (v(1)<=0)) %east facets cant see anything to the west
                            test=0;
                        elseif ((fi==4) && (v(2)<=0)) %north facets cant see anything to the south
                            test=0;
                        elseif ((fi==5) && (v(2)>=0)) %south facets cant see anything to the north
                            test=0;
                        elseif ((fi~=1) && (fi2==1) && (v(3)>0)) %vertical walls cant see any horizontal facet that is higher
                            test=0;
                            % the following cases should be blocked by a block anyway, but
                            % this is faster
                        elseif ((fi==4) || (fi==5)) && (fi2==3) && (v(1)>0) %north/south can't see an east facet if it is east of it
                            test=0;
                        elseif ((fi==4) || (fi==5)) && (fi2==2) && (v(1)<0) %north/south can't see a west facet if it is west of it
                            test=0;
                        elseif ((fi==2) || (fi==3)) && (fi2==4) && (v(2)>0) %west/east can't see a north facet if it is north of it
                            test=0;
                        elseif ((fi==2) || (fi==3)) && (fi2==5) && (v(2)<0) %west/east can't see a south facet if it is south of it
                            test=0;
                        end
                        
                        %disp(['test: ' num2str(test)])
                        
                        if test %facets could see each other, check if blocked
                            %calculate actual intersection distance
                            
                            for n=1:nblocks %check if any block intersects
                                if ((n==bi) && (ci>0))  %don't intersect with itself
                                    continue
                                end
                                
                                [flag,dints] = preprocessing.rbi(facetpoint, v1, blocks_phys(n,:));
                                
                                if dints<=0  %intersection is downstream
                                    flag=0;
                                elseif dints>=dint %block is behind target subfacet
                                    flag=0;
                                end
                                if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                                    %  %disp(['intersection with block: ' num2str(n)])
                                    prblckdij=1;
                                    break %get out of inner for loop
                                end
                            end
                            
                        else %facets can't possibly see each other
                            prblckdij=1;
                        end
                    end
                    
                    %%
                    c = 1 - prblckdij;
                    %pf1sf2(i, j) = c;
                    %pf1sf2(j, i) = c;
                    %A(i) = areaa;
                    %A(j) = areab;

                    if prblckdij < 0.98 %view is not blocked, vf calclation
                        if (prblckdij >= (2 * cornerweight - eps)) %at least two corners have to be blocked, otherwise it cannot possibly be blocked by itself
                            %slice up
                            %pf1sf2u also has to include check for center, in case center is
                            %blocked
                            [coaa, cobb, pf1sf2u] = preprocessing.slice(fi, fi2, coa, cob, cornerweight, centerweight);
                            c = c + pf1sf2u;
                            %pf1sf2(i, j) = c;
                            %pf1sf2(j, i) = c;
                        else
                            coaa = coa(6:9, :);
                            cobb = cob(6:9, :);
                        end

                        %determine if they form a corner/have common edge
                        %determine the vertices pair that forms that corner
                        vcorner = [0, 0, 1]; %corner vertices for i and j & sign of correction (as always, clockwise from bottom left)
                        glpo = 6; %order of gauss-legendre polynomial
                        glpop = 15; %more accurate, order of gauss-legendre polynomial
                        if (fi == 1) %fi is horizontal, test if corner with a vertical wall
                            if  (fi2 == 2)
                                if all(coaa(4, :) == cobb(1, :)) && all(coaa(3, :) == cobb(4, :)) %they share one edge completely
                                    vcorner = [3, 4, -1];
                                    glpo = glpop;
                                    % elseif (coaa(4,1)==cobb(1,1)) && (coaa(4,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %     glpo=18;  %NEEDS TO BE AN EVEN NUMBER!!!
                                end
                            elseif (fi2 == 3)
                                if all(coaa(1, :) == cobb(1, :)) && all(coaa(2, :) == cobb(4, :)) %they share one edge completely
                                    vcorner = [1, 4, 1];
                                    glpo = glpop;
                                    % elseif (coaa(1,1)==cobb(1,1)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %     glpo=18;
                                end
                            elseif (fi2==4)
                                if all(coaa(1,:) == cobb(1,:)) && all(coaa(4,:)==cobb(4,:)) %they share one edge completely
                                    vcorner = [4, 4, -1];
                                    glpo = glpop;
                                    %  elseif (coaa(1,2)==cobb(1,2)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %      glpo=18;
                                end
                            elseif (fi2 == 5)
                                if all(coaa(2, :) == cobb(1,:)) && all(coaa(3, :) == cobb(4, :)) %they share one edge completely
                                    vcorner = [2, 4, 1];
                                    glpo = glpop;
                                    %  elseif (coaa(2,2)==cobb(1,2)) && (coaa(2,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %      glpo=18;
                                end
                            end
                        elseif (fi == 2)
                            if  (fi2 == 4)
                                if all(coaa(1, :) == cobb(4, :)) && all(coaa(2, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [1, 3, 1];
                                    glpo = glpop;
                                elseif (coaa(1, 1) == cobb(4, 1)) && (coaa(1, 2) ==cobb(4, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 5)
                                if all(coaa(4, :) == cobb(4, :)) && all(coaa(3, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [3, 3, -1];
                                    glpo = glpop;
                                elseif (coaa(4,1)==cobb(4,1)) && (coaa(4,2)==cobb(4,2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 1)
                                if all(coaa(1, :) == cobb(4, :)) && all(coaa(4, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [4, 3, -1];
                                    glpo = glpop;
                                    %                    elseif (coaa(1,1)==cobb(4,1)) && (coaa(1,3)==cobb(4,3)) %they potentially share parts of an edge
                                    %                        glpo=18;
                                end
                             elseif (fi2 == 6)
                                if all(coaa(2, :) == cobb(4, :)) && all(coaa(3, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [2, 3, 1];
                                    glpo = glpop;
                                    %                    elseif (coaa(1,1)==cobb(4,1)) && (coaa(1,3)==cobb(4,3)) %they potentially share parts of an edge
                                    %                        glpo=18;
                                end
                            end
                        elseif (fi == 3)
                            if  (fi2 == 4)
                                if all(coaa(1, :) == cobb(1, :)) && all(coaa(2, :) == cobb(2,:)) %they share one edge completely
                                    vcorner = [1, 1, -1];
                                    glpo = glpop;
                                elseif (coaa(1,1)==cobb(1,1)) && (coaa(1,2)==cobb(1,2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 5)
                                if all(coaa(4, :) == cobb(1, :)) && all(coaa(3, :) == cobb(2, :)) %they share one edge completely
                                    vcorner = [3, 1, 1];
                                    glpo = glpop;
                                elseif (coaa(4, 1) == cobb(1, 1)) && (coaa(4, 2) == cobb(1, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 1)
                                if all(coaa(1, :) == cobb(1, :)) && all(coaa(4, :) == cobb(2, :)) %they share one edge completely
                                    vcorner = [4, 1, 1];
                                    glpo = glpop;
                                    %  elseif (coaa(1,1)==cobb(1,1)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %      glpo=18;
                                end
                            elseif (fi2 == 6)
                                if all(coaa(2, :) == cobb(1, :)) && all(coaa(3, :) == cobb(2, :)) %they share one edge completely
                                    vcorner = [2, 1, -1];
                                    glpo = glpop;
                                    %  elseif (coaa(1,1)==cobb(1,1)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %      glpo=18;
                                end
                            end
                        elseif (fi == 4)
                            if  (fi2 == 2)
                                if all(coaa(4, :) == cobb(1, :)) && all(coaa(3, :) == cobb(2, :)) %they share one edge completely
                                    vcorner = [3, 1, 1];
                                    glpo = glpop;
                                elseif (coaa(4, 1) == cobb(1, 1)) && (coaa(4, 2) == cobb(1, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 3)
                                if all(coaa(1, :) == cobb(1, :)) && all(coaa(2, :) == cobb(2, :)) %they share one edge completely
                                    vcorner = [1, 1, -1];
                                    glpo = glpop;
                                elseif (coaa(1, 1) == cobb(1, 1)) && (coaa(1, 2) == cobb(1, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 1)
                                if all(coaa(1, :) == cobb(1, :)) && all(coaa(4, :) == cobb(4, :)) %they share one edge completely
                                    vcorner = [4, 4, -1];
                                    glpo = glpop;
                                    %    elseif (coaa(1,2)==cobb(1,2)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %       glpo=18;
                                end
                            elseif (fi2 == 6)
                                if all(coaa(2, :) == cobb(1, :)) && all(coaa(3, :) == cobb(4, :)) %they share one edge completely
                                    vcorner = [2, 4, 1];
                                    glpo = glpop;
                                    %    elseif (coaa(1,2)==cobb(1,2)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %       glpo=18;
                                end
                            end
                        elseif (fi == 5)
                            if  (fi2 == 2)
                                if all(coaa(4, :) == cobb(4, :)) && all(coaa(3, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [3, 3, -1];
                                    glpo = glpop;
                                elseif (coaa(4, 1) == cobb(4, 1)) && (coaa(4, 2) == cobb(4, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 3)
                                if all(coaa(1, :) == cobb(4, :)) && all(coaa(2, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [1, 3, 1];
                                    glpo = glpop;
                                elseif (coaa(1, 1) == cobb(4, 1)) && (coaa(1, 2) == cobb(4, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 1)
                                if all(coaa(1, :) == cobb(2, :)) && all(coaa(4, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [4, 2, 1];
                                    glpo = glpop;
                                    % elseif (coaa(1,2)==cobb(2,2)) && (coaa(1,3)==cobb(2,3)) %they potentially share parts of an edge
                                    %     glpo=18
                                end
                            elseif (fi2 == 6)
                                if all(coaa(2, :) == cobb(2, :)) && all(coaa(3, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [2, 2, -1];
                                    glpo = glpop;
                                    % elseif (coaa(1,2)==cobb(2,2)) && (coaa(1,3)==cobb(2,3)) %they potentially share parts of an edge
                                    %     glpo=18
                                end
                            end
                        elseif (fi == 6)
                            if  (fi2 == 2)
                                if all(coaa(4, :) == cobb(2, :)) && all(coaa(3, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [3, 2, 1];
                                    glpo = glpop;
                                    % elseif (coaa(4,1)==cobb(1,1)) && (coaa(4,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %     glpo=18;  %NEEDS TO BE AN EVEN NUMBER!!!
                                end
                            elseif (fi2 == 3)
                                if all(coaa(1, :) == cobb(2, :)) && all(coaa(2, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [1, 2, -1];
                                    glpo = glpop;
                                    % elseif (coaa(1,1)==cobb(1,1)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %     glpo=18;
                                end
                            elseif (fi2==4)
                                if all(coaa(1,:) == cobb(2,:)) && all(coaa(4,:)==cobb(3,:)) %they share one edge completely
                                    vcorner = [4, 2, 1];
                                    glpo = glpop;
                                    %  elseif (coaa(1,2)==cobb(1,2)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %      glpo=18;
                                end
                            elseif (fi2 == 5)
                                if all(coaa(2, :) == cobb(2,:)) && all(coaa(3, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [2, 2, -1];
                                    glpo = glpop;
                                    %  elseif (coaa(2,2)==cobb(1,2)) && (coaa(2,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %      glpo=18;
                                end
                            end                            
                        end

                        %calculate viewfactor
                        [F12, F21] = preprocessing.ViewFactor(coaa, cobb, areaa, areab, glpo, vcorner);
                        %increase glpo if view factors are big
                        if F12 > 0.5 || F21 > 0.5
                            [F12, F21] = preprocessing.ViewFactor(coaa, cobb, areaa, areab, glpo + 10, vcorner);
                        end


                        % round to 1promille accuracy? i.e. cut smaller 1promille completely
                        % this violates reciprocity though
                        vf(i,j) = F12 * c;
                        vf(j,i) = F21 * c;
                        %   vf(i,j)=floor(F12*c*100)/100;
                        %  vf(j,i)=floor(F21*c*100)/100; 
                    end
                end
                %count
            end
                       
            vf = single(vf);
            %pf1sf2o=single(pf1sf2);
            blub = sum(vf, 2); % Sum of view factors
            lblub = find(blub > 1); % Position of where view factor sum is greater than 1
            if ~isempty(lblub)
                %[maxvf, indexmaxvf] = max(blub);
                for i = 1:length(lblub) % For each facet with view factors summing to more than 1
                    vf(lblub(i), :) = vf(lblub(i),:)/blub(lblub(i)); % Divide by the sum
                end
            end          
            
            preprocessing.addvar(obj, 'vf', vf)
            preprocessing.addvar(obj, 'svf', max(1 - sum(vf, 2), 0));
            %preprocessing.addvar(obj, 'facetarea', A);
        end
        
        function vfcAcceleratedPar(obj)

            nblocks = obj.nblocks;
            facets = obj.facets; nfcts = obj.nfcts;
            blocks_phys = obj.blocks_phys;
            centerweight = obj.centerweight; cornerweight = obj.cornerweight; %delta = 0.01;
            facetarea = obj.facetarea;
            facetdim = obj.facetdim;
            facetcorner = obj.facetcorner;
            
            min_blocks_x = blocks_phys(:,1); min_blocks_y = blocks_phys(:,3); min_blocks_z = blocks_phys(:,5);
            max_blocks_x = blocks_phys(:,2); max_blocks_y = blocks_phys(:,4); max_blocks_z = blocks_phys(:,6);                  
            
            vf = zeros(nfcts, nfcts);
            vfp = zeros(nfcts, nfcts);
          
            parfor i = 1:(nfcts - 1)
                i
                bi = facets(i, 3); %block
                fi = facets(i, 1); %orientation
                ci = facets(i, 4); %building index (-1 for roads, -99 for bounding wall)

                ndima = facetdim(i);
                areaa = facetarea(i);
                coa = squeeze(facetcorner(i,:,:));
                %% Geometry tests
                facids = (i+1):nfcts;
                facids_o = (i+1):nfcts;
                facids_i = (i+1):nfcts;
                facids_og = (i+1):nfcts;
                % Elimination due to orientation of i
                
                % Must be different orientation;
                facids_o = facids_o(facets(facids_o,1) ~= facets(i,1));
                facids_o = facids_o(facets(facids_o,3) .* facets(facids_o,4) ~= facets(i,3) * facets(i,4));
                
                facids = facids(facets(facids,1) ~= facets(i,1));
                facids = facids(facets(facids,3) .* facets(facids,4) ~= facets(i,3) * facets(i,4));
                
                switch fi
                    case 1 % Top: maximum zs must be greater than z.
                        facids_i = facids_i(facets(facids_i,11) > facets(i,10));
                        facids = facids(facets(facids,11) > facets(i,10));
                    case 4 % North: maximum ys must be greater than y
                        facids_i = facids_i(facets(facids_i,9) > facets(i,8));
                        facids = facids(facets(facids,9) > facets(i,8));
                    case 3 % East: maximum xs must be greater than x
                        facids_i = facids_i(facets(facids_i,7) > facets(i,6));
                        facids = facids(facets(facids,7) > facets(i,6));
                    case 5 % South: minimum ys must be less than y
                        facids_i = facids_i(facets(facids_i,8) < facets(i,9));
                        facids = facids(facets(facids,8) < facets(i,9));
                    case 2 % West: minimum xs must be less than x
                        facids_i =  facids_i(facets(facids_i,6) < facets(i,7));
                        facids =  facids(facets(facids,6) < facets(i,7));
                    case 6 % Bottom: minimum zs must be less than z
                        facids_i = facids_i(facets(facids_i,10) < facets(i,11));
                        facids = facids(facets(facids,10) < facets(i,11));
                end
                
                facids_j = [];
                facids_jog = [];
                % Elimination due to orientation of j
                if fi ~= 1 % Top: maximum z must be greater than zs;
                     facids_j = [facids_j, intersect(facids(facets(facids,1) == 1), facids(facets(i,11) > facets(facids, 10)))];
                end
                if fi ~= 4 % North:
                    facids_j  = [facids_j, intersect(facids(facets(facids,1) == 4), facids(facets(i,9) > facets(facids,8)))];                    
                end
                if fi ~= 3 % East:
                    facids_j = [facids_j, intersect(facids(facets(facids,1) == 3), facids(facets(i,7) > facets(facids,6)))];
                end
                if fi ~= 5 % South
                    facids_j = [facids_j, intersect(facids(facets(facids,1) == 5), facids(facets(i,8) < facets(facids,9)))];
                end
                if fi ~= 2 % West
                    facids_j = [facids_j, intersect(facids(facets(facids,1) == 2), facids(facets(i,6) < facets(facids,7)))];
                end
                if fi ~= 6 % Bottom
                     facids_j = [facids_j, intersect(facids(facets(facids,1) == 6), facids(facets(i,10) < facets(facids, 11)))];
                end
                
                facids = intersect(facids_j, facids);              
                if isempty(facids)
                     disp('No calculation')
                     continue
                end
                
                vfj = zeros(1,nfcts);
                vfjp = zeros(1,nfcts);
                for dum = 1:length(facids)
                    j = facids(dum);
                    fi2 = facets(j, 1); %orientation

                    ndimb = facetdim(j);
                    areab = facetarea(j);
                    cob = squeeze(facetcorner(j,:,:));
                    
                    testcrit = 0;
                    prblckdij = 0;
                    
                    if (ndima+ndimb)>(2*testcrit)  %also test corner, otherwise only test center
                        ordera=[1 2 3 4 5];
                        cas=fi*10+fi2;
                        switch cas
                            case 12, orderb=[1 3 4 5 2];
                            case 13, orderb=[1 2 5 4 3];
                            case 15, orderb=[1 3 2 5 4];
                            case 21, orderb=[1 5 2 3 4];
                            case 24, orderb=[1 5 4 3 2];
                            case 26, orderb=[1 2 5 4 3];
                            case 31, orderb=[1 2 5 4 3];
                            case 35, orderb=[1 5 4 3 2];
                            case 36, orderb=[1 5 2 3 4];
                            case 42, orderb=[1 5 4 3 2];
                            case 46, orderb=[1 3 2 5 4];
                            case 51, orderb=[1 3 2 5 4];
                            case 53, orderb=[1 5 4 3 2];
                            case 62, orderb=[1 2 5 4 3];
                            case 63, orderb=[1 3 4 5 2];
                            case 64, orderb=[1 3 2 5 4];
                            otherwise %14, 41, 23, 32, 25, 52, 34, 43, 45, 54
                                orderb=[1 2 3 4 5];
                        end
                        
                        for k=1:5
                            test=1;
                            oa=ordera(k);
                            ob=orderb(k);
                            facetpoint=coa(oa,:);
                            facetpoint2=cob(ob,:);
                            v=[facetpoint2(1)-facetpoint(1),facetpoint2(2)-facetpoint(2),facetpoint2(3)-facetpoint(3)];
                            dint=norm(v);
                            v1=v/dint;
                            % Add test for bottom facets
                            if ((fi==1) && (v(3)<=0)) %horizontal facets can't see facet that is lower.
                                test=0;
                            elseif ((fi==2) && (v(1)>=0)) %west facets cant see anything to the east
                                test=0;
                            elseif ((fi==3) && (v(1)<=0)) %east facets cant see anything to the west
                                test=0;
                            elseif ((fi==4) && (v(2)<=0)) %north facets cant see anything to the south
                                test=0;
                            elseif ((fi==5) && (v(2)>=0)) %south facets cant see anything to the north
                                test=0;
                            elseif ((fi~=1) && (fi2==1) && (v(3)>0)) %vertical walls cant see any horizontal facet that is higher
                                test=0;
                                % the following cases should be blocked by a block anyway, but
                                % this is faster
                            elseif ((fi==4) || (fi==5)) && (fi2==3) && (v(1)>0) %north/south can't see an east facet if it is east of it
                                test=0;
                            elseif ((fi==4) || (fi==5)) && (fi2==2) && (v(1)<0) %north/south can't see a west facet if it is west of it
                                test=0;
                            elseif ((fi==2) || (fi==3)) && (fi2==4) && (v(2)>0) %west/east can't see a north facet if it is north of it
                                test=0;
                            elseif ((fi==2) || (fi==3)) && (fi2==5) && (v(2)<0) %west/east can't see a south facet if it is south of it
                                test=0;
                            end
                            
                            if test
                                or = facetpoint;
                                d = v1;
                                
                                [tmin_vec,flag_vec] = rayBoxVec(min_blocks_x, min_blocks_y, min_blocks_z, max_blocks_x, max_blocks_y, max_blocks_z, or(1), or(2), or(3), d(1), d(2), d(3));
                                                                
                                flag_vec(bi) = 0;
                                flag_vec(tmin_vec <= eps) = 0;
                                flag_vec(tmin_vec >= dint - eps) = 0;
                                                                
                                if any(flag_vec)
                                    if k==1
                                        prblckdij=centerweight;
                                    else
                                        prblckdij=prblckdij+cornerweight;
                                    end
                                end
                                
                            else %facets can't possibly see each other
                                if k==1
                                    prblckdij=centerweight;
                                else
                                    prblckdij=prblckdij+cornerweight;
                                end
                            end
                        end
                        
                    else %small area, checking only center. Currently hard-coded to never happen but should do vectorisation eventualy
                        % %disp('small area, checkig only center'
                        test=1;
                        facetpoint=coa(1,:);
                        facetpoint2=cob(1,:);
                        v=[facetpoint2(1)-facetpoint(1),facetpoint2(2)-facetpoint(2),facetpoint2(3)-facetpoint(3)];
                        dint=norm(v);
                        v1=v/dint;
                        
                        %abs(v)
                        %%disp(['v1: ' num2str(v1)])
                        %test if facets can potentially see each other
                        if ((fi==1) && (v(3)<=0)) %horizontal facets can't see facets whos center is lower
                            test=0;
                        elseif ((fi==2) && (v(1)>=0)) %west facets cant see anything to the east
                            test=0;
                        elseif ((fi==3) && (v(1)<=0)) %east facets cant see anything to the west
                            test=0;
                        elseif ((fi==4) && (v(2)<=0)) %north facets cant see anything to the south
                            test=0;
                        elseif ((fi==5) && (v(2)>=0)) %south facets cant see anything to the north
                            test=0;
                        elseif ((fi~=1) && (fi2==1) && (v(3)>0)) %vertical walls cant see any horizontal facet that is higher
                            test=0;
                            % the following cases should be blocked by a block anyway, but
                            % this is faster
                        elseif ((fi==4) || (fi==5)) && (fi2==3) && (v(1)>0) %north/south can't see an east facet if it is east of it
                            test=0;
                        elseif ((fi==4) || (fi==5)) && (fi2==2) && (v(1)<0) %north/south can't see a west facet if it is west of it
                            test=0;
                        elseif ((fi==2) || (fi==3)) && (fi2==4) && (v(2)>0) %west/east can't see a north facet if it is north of it
                            test=0;
                        elseif ((fi==2) || (fi==3)) && (fi2==5) && (v(2)<0) %west/east can't see a south facet if it is south of it
                            test=0;
                        end
                                                
                        if test %facets could see each other, check if blocked
                            %calculate actual intersection distance
                            
                            for n=1:nblocks %check if any block intersects
                                if ((n==bi) && (ci>0))  %don't intersect with itself
                                    continue
                                end
                                
                                [flag,dints] = preprocessing.rbi(facetpoint, v1, blocks_phys(n,:));
                                
                                if dints<=0  %intersection is downstream
                                    flag=0;
                                elseif dints>=dint %block is behind target subfacet
                                    flag=0;
                                end
                                if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                                    %  %disp(['intersection with block: ' num2str(n)])
                                    prblckdij=1;
                                    break %get out of inner for loop
                                end
                            end
                            
                        else %facets can't possibly see each other
                            prblckdij=1;
                        end
                    end
                    
                    c = 1 - prblckdij;
                    if prblckdij < 0.98 %view is not blocked, vf calclation
                        if (prblckdij >= (2 * cornerweight - eps)) %at least two corners have to be blocked, otherwise it cannot possibly be blocked by itself
                            %slice up
                            %pf1sf2u also has to include check for center, in case center is
                            %blocked
                            [coaa, cobb, pf1sf2u] = preprocessing.slice(fi, fi2, coa, cob, cornerweight, centerweight);
                            c = c + pf1sf2u;
                        else
                            coaa = coa(6:9, :);
                            cobb = cob(6:9, :);
                        end

                        %determine if they form a corner/have common edge
                        %determine the vertices pair that forms that corner
                        vcorner = [0, 0, 1]; %corner vertices for i and j & sign of correction (as always, clockwise from bottom left)
                        glpo = 6; %order of gauss-legendre polynomial
                        glpop = 15; %more accurate, order of gauss-legendre polynomial
                        if (fi == 1) %fi is horizontal, test if corner with a vertical wall
                            if  (fi2 == 2)
                                if all(coaa(4, :) == cobb(1, :)) && all(coaa(3, :) == cobb(4, :)) %they share one edge completely
                                    vcorner = [3, 4, -1];
                                    glpo = glpop;
                                    % elseif (coaa(4,1)==cobb(1,1)) && (coaa(4,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %     glpo=18;  %NEEDS TO BE AN EVEN NUMBER!!!
                                end
                            elseif (fi2 == 3)
                                if all(coaa(1, :) == cobb(1, :)) && all(coaa(2, :) == cobb(4, :)) %they share one edge completely
                                    vcorner = [1, 4, 1];
                                    glpo = glpop;
                                    % elseif (coaa(1,1)==cobb(1,1)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %     glpo=18;
                                end
                            elseif (fi2==4)
                                if all(coaa(1,:) == cobb(1,:)) && all(coaa(4,:)==cobb(4,:)) %they share one edge completely
                                    vcorner = [4, 4, -1];
                                    glpo = glpop;
                                    %  elseif (coaa(1,2)==cobb(1,2)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %      glpo=18;
                                end
                            elseif (fi2 == 5)
                                if all(coaa(2, :) == cobb(1,:)) && all(coaa(3, :) == cobb(4, :)) %they share one edge completely
                                    vcorner = [2, 4, 1];
                                    glpo = glpop;
                                    %  elseif (coaa(2,2)==cobb(1,2)) && (coaa(2,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %      glpo=18;
                                end
                            end
                        elseif (fi == 2)
                            if  (fi2 == 4)
                                if all(coaa(1, :) == cobb(4, :)) && all(coaa(2, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [1, 3, 1];
                                    glpo = glpop;
                                elseif (coaa(1, 1) == cobb(4, 1)) && (coaa(1, 2) ==cobb(4, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 5)
                                if all(coaa(4, :) == cobb(4, :)) && all(coaa(3, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [3, 3, -1];
                                    glpo = glpop;
                                elseif (coaa(4,1)==cobb(4,1)) && (coaa(4,2)==cobb(4,2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 1)
                                if all(coaa(1, :) == cobb(4, :)) && all(coaa(4, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [4, 3, -1];
                                    glpo = glpop;
                                    %                    elseif (coaa(1,1)==cobb(4,1)) && (coaa(1,3)==cobb(4,3)) %they potentially share parts of an edge
                                    %                        glpo=18;
                                end
                             elseif (fi2 == 6)
                                if all(coaa(2, :) == cobb(4, :)) && all(coaa(3, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [2, 3, 1];
                                    glpo = glpop;
                                    %                    elseif (coaa(1,1)==cobb(4,1)) && (coaa(1,3)==cobb(4,3)) %they potentially share parts of an edge
                                    %                        glpo=18;
                                end
                            end
                        elseif (fi == 3)
                            if  (fi2 == 4)
                                if all(coaa(1, :) == cobb(1, :)) && all(coaa(2, :) == cobb(2,:)) %they share one edge completely
                                    vcorner = [1, 1, -1];
                                    glpo = glpop;
                                elseif (coaa(1,1)==cobb(1,1)) && (coaa(1,2)==cobb(1,2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 5)
                                if all(coaa(4, :) == cobb(1, :)) && all(coaa(3, :) == cobb(2, :)) %they share one edge completely
                                    vcorner = [3, 1, 1];
                                    glpo = glpop;
                                elseif (coaa(4, 1) == cobb(1, 1)) && (coaa(4, 2) == cobb(1, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 1)
                                if all(coaa(1, :) == cobb(1, :)) && all(coaa(4, :) == cobb(2, :)) %they share one edge completely
                                    vcorner = [4, 1, 1];
                                    glpo = glpop;
                                    %  elseif (coaa(1,1)==cobb(1,1)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %      glpo=18;
                                end
                            elseif (fi2 == 6)
                                if all(coaa(2, :) == cobb(1, :)) && all(coaa(3, :) == cobb(2, :)) %they share one edge completely
                                    vcorner = [2, 1, -1];
                                    glpo = glpop;
                                    %  elseif (coaa(1,1)==cobb(1,1)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %      glpo=18;
                                end
                            end
                        elseif (fi == 4)
                            if  (fi2 == 2)
                                if all(coaa(4, :) == cobb(1, :)) && all(coaa(3, :) == cobb(2, :)) %they share one edge completely
                                    vcorner = [3, 1, 1];
                                    glpo = glpop;
                                elseif (coaa(4, 1) == cobb(1, 1)) && (coaa(4, 2) == cobb(1, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 3)
                                if all(coaa(1, :) == cobb(1, :)) && all(coaa(2, :) == cobb(2, :)) %they share one edge completely
                                    vcorner = [1, 1, -1];
                                    glpo = glpop;
                                elseif (coaa(1, 1) == cobb(1, 1)) && (coaa(1, 2) == cobb(1, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 1)
                                if all(coaa(1, :) == cobb(1, :)) && all(coaa(4, :) == cobb(4, :)) %they share one edge completely
                                    vcorner = [4, 4, -1];
                                    glpo = glpop;
                                    %    elseif (coaa(1,2)==cobb(1,2)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %       glpo=18;
                                end
                            elseif (fi2 == 6)
                                if all(coaa(2, :) == cobb(1, :)) && all(coaa(3, :) == cobb(4, :)) %they share one edge completely
                                    vcorner = [2, 4, 1];
                                    glpo = glpop;
                                    %    elseif (coaa(1,2)==cobb(1,2)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %       glpo=18;
                                end
                            end
                        elseif (fi == 5)
                            if  (fi2 == 2)
                                if all(coaa(4, :) == cobb(4, :)) && all(coaa(3, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [3, 3, -1];
                                    glpo = glpop;
                                elseif (coaa(4, 1) == cobb(4, 1)) && (coaa(4, 2) == cobb(4, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 3)
                                if all(coaa(1, :) == cobb(4, :)) && all(coaa(2, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [1, 3, 1];
                                    glpo = glpop;
                                elseif (coaa(1, 1) == cobb(4, 1)) && (coaa(1, 2) == cobb(4, 2)) %they potentially share parts of an edge
                                    glpo = 18;
                                end
                            elseif (fi2 == 1)
                                if all(coaa(1, :) == cobb(2, :)) && all(coaa(4, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [4, 2, 1];
                                    glpo = glpop;
                                    % elseif (coaa(1,2)==cobb(2,2)) && (coaa(1,3)==cobb(2,3)) %they potentially share parts of an edge
                                    %     glpo=18
                                end
                            elseif (fi2 == 6)
                                if all(coaa(2, :) == cobb(2, :)) && all(coaa(3, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [2, 2, -1];
                                    glpo = glpop;
                                    % elseif (coaa(1,2)==cobb(2,2)) && (coaa(1,3)==cobb(2,3)) %they potentially share parts of an edge
                                    %     glpo=18
                                end
                            end
                        elseif (fi == 6)
                            if  (fi2 == 2)
                                if all(coaa(4, :) == cobb(2, :)) && all(coaa(3, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [3, 2, 1];
                                    glpo = glpop;
                                    % elseif (coaa(4,1)==cobb(1,1)) && (coaa(4,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %     glpo=18;  %NEEDS TO BE AN EVEN NUMBER!!!
                                end
                            elseif (fi2 == 3)
                                if all(coaa(1, :) == cobb(2, :)) && all(coaa(2, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [1, 2, -1];
                                    glpo = glpop;
                                    % elseif (coaa(1,1)==cobb(1,1)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %     glpo=18;
                                end
                            elseif (fi2==4)
                                if all(coaa(1,:) == cobb(2,:)) && all(coaa(4,:)==cobb(3,:)) %they share one edge completely
                                    vcorner = [4, 2, 1];
                                    glpo = glpop;
                                    %  elseif (coaa(1,2)==cobb(1,2)) && (coaa(1,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %      glpo=18;
                                end
                            elseif (fi2 == 5)
                                if all(coaa(2, :) == cobb(2,:)) && all(coaa(3, :) == cobb(3, :)) %they share one edge completely
                                    vcorner = [2, 2, -1];
                                    glpo = glpop;
                                    %  elseif (coaa(2,2)==cobb(1,2)) && (coaa(2,3)==cobb(1,3)) %they potentially share parts of an edge
                                    %      glpo=18;
                                end
                            end                            
                        end

                        %calculate viewfactor
                        [F12, F21] = preprocessing.ViewFactor(coaa, cobb, areaa, areab, glpo, vcorner);
                        %increase glpo if view factors are big
                        if F12 > 0.5 || F21 > 0.5
                            [F12, F21] = preprocessing.ViewFactor(coaa, cobb, areaa, areab, glpo + 10, vcorner);
                        end

                        vfj(j) = F12 * c;
                        vfjp(j) = F21 * c;
                    end
                end

                vf(i,:) = vfj;
                vfp(i,:) = vfjp;
            end
                       
            vf = single(vf + vfp');
            blub = sum(vf, 2); % Sum of view factors
            lblub = find(blub > 1); % Position of where view factor sum is greater than 1
            if ~isempty(lblub)
                for i = 1:length(lblub) % For each facet with view factors summing to more than 1
                    vf(lblub(i), :) = vf(lblub(i),:)/blub(lblub(i)); % Divide by the sum
                end
            end          
            
            preprocessing.addvar(obj, 'vf', vf)
            preprocessing.addvar(obj, 'svf', max(1 - sum(vf, 2), 0));
        end

        function write_svf(obj)
            fname = ['svf.inp.' num2str(obj.expnr)];
            fileID = fopen(fname,'W');
            fprintf(fileID, '# sky view factors\n');
            fclose(fileID);
            dlmwrite(fname, obj.svf, '-append','delimiter',' ','precision','%4f')
        end

        function write_vf(obj)
            ncid = netcdf.create(['vf.nc.inp.' num2str(obj.expnr)], 'NC_WRITE');
            dimidrow = netcdf.defDim(ncid,'rows', obj.nfcts);
            dimidcol = netcdf.defDim(ncid,'columns', obj.nfcts);
            varid = netcdf.defVar(ncid,'view factor','NC_FLOAT',[dimidrow dimidcol]);
            netcdf.endDef(ncid);
            netcdf.putVar(ncid,varid,obj.vf);
            netcdf.close(ncid);
        end

        function write_facetarea(obj)
            fname = ['facetarea.inp.' num2str(obj.expnr)];
            fileID = fopen(fname,'W');
            fprintf(fileID, '# area of facets\n');
            fclose(fileID);
            dlmwrite(fname, obj.facetarea,'-append','delimiter',' ','precision','%4f')
        end

        function rayit(obj)
            facets = obj.facets;
            walltypes = obj.walltypes;
             function [fct, wall] = loadfacets()
                %M = dlmread(['facets.inp.' num2str(expnr)],'',1,0);
                vars = {'or', 'wlid', 'blk', 'bld'};
                for n = 1:length(vars)
                    fct.(vars{n}) = facets(:, n);
                end

                %disp(obj.walltypes)
                %M = dlmread(['walltypes.inp.', obj.expnr],'',3,0);
                %disp(M)
                % wallid  lGR   z0 [m]     z0h [m]     al [-]   em [-]   d1 [m]  d2 [m]    cp1 [J/kgK]  cp2 [J/kgK]   l1 [W/(m K)]  l2 [W/(m K)]    k1 [W/(m K)]    k1 [W/(m K)]
                vars = {'id','lGR','z0','z0h','al', 'em', 'd1', 'd2', 'cp1', 'cp2', 'l1', 'l2', 'k1', 'k2'};
                for n = 1:length(vars)
                    wall.(vars{n}) = walltypes(:, n);
                end

                % check whether all referenced wall ids have been defined
                idref = unique(fct.wlid);
                iddef = unique(wall.id');

                if length(wall.id) > length(iddef)
                    disp('ERROR: multiple definitions of a wall id')
                    return
                end
                if length(idref) > length(iddef)
                    disp('ERROR: more walltypes used than defined')
                    return
                end


                % assign the wall type to each of the facets -- this way it is easy to know
                % in which structure the proporties can be found.
                fct.wltp = zeros(size(fct.or));

                for m = 1:size(fct.wlid, 1)
                    if fct.wlid(m) >= 0 && fct.wlid(m) <= 10 && any(ismember(iddef, fct.wlid(m)))
                        fct.wltp(n) = 1; %normal
                    elseif fct.wlid(m) > 10 && any(ismember(iddef, fct.wlid(m)))
                        fct.wltp(n) = 2; %Green roof
                    elseif fct.wlid(m) < 0 && any(ismember(iddef, fct.wlid(m)))
                        if fct.wlid(m) <= -99
                            fct.wltp(n) = 3; %bounding wall
                        else
                            fct.wltp(n) = 4; %floor
                        end
                    else
                        disp(['ERROR: walltype of facet ',num2str(m),' not defined'])
                    end
                end


                for n = 1:length(idref)
                    if (ismember(idref(n), wall.id))
                        fct.wltp(fct.wlid == idref(n)) = 1;
                    elseif (ismember(idref(n), gr.id))
                        fct.wltp(fct.wlid == idref(n)) = 2;
                    else
                        sprintf('ERROR. wall id %8d not defined. Correct wall and greenroof property input files', idref(n))
                        return
                    end
                end
            end


        %facets = obj.facets;
        [fct, wall] = loadfacets();



        %[sortt, sorti]=sort(facets(:,1));  %sort by walltype

        %% shortwave
        %

        %Z=35;          %zenith angle of the sun (could be function of time, location etc)
        %Dsky=zeros(obj.nfcts,1);
        %Denv=zeros(nfcts,1);
        nfcts = obj.nfcts;
        albedo = zeros(nfcts,1);
        emissivity = zeros(nfcts,1);
        Sdir = obj.Sdir;
        Dsk = obj.Dsk;
        svf = obj.svf;

        for i=1:nfcts
            j=find(fct.wlid(i)==wall.id);
            albedo(i)=wall.al(j);
            emissivity(i)=wall.em(j);
        end
        preprocessing.addvar(obj, 'emissivity', emissivity);

        %isroof
        isnotroof=ones(nfcts,1);
        isnotroof(find(facets(:,1)==1 & facets(:,4)>0))=0;

        %diffuse flux from sky and other walls
        Kinnew=zeros(nfcts,1);

        Kininit=(1-albedo).*(Sdir+Dsk.*svf);
        Koutinit=albedo.*(Sdir+Dsk.*svf);
        %find(Kininit == 0)
        %total radiation absorbed and reflected
        sum(Kininit+Koutinit);
        %total radiation absorbed
        sum(Kininit);
        %total radiation reflected
        sum(Koutinit);



        totincrease=zeros(10,1);
        increase=zeros(10,1);
        Kout=zeros(10,1);
        Kout(2)=sum(Koutinit);


        Kin=Kininit; %Kin §adds up
        Kinold=Kininit;
        Koutold=Koutinit; %Kout is wiped clean with every reflection
        Koutnew=zeros(nfcts,1);
        Kintemp=zeros(nfcts,1);
        Kouttemp=zeros(nfcts,1);
        count=0;
        storerad=zeros(nfcts,nfcts);
        Kiterin=zeros(300,1);
        Kiterout=zeros(300,1);
        itermaxdiff=zeros(300,1);
        itermaxdiffloc=zeros(300,1);
        itermaxrelloc=zeros(300,1);
        moep=zeros(nfcts,20);
        blub=zeros(nfcts,20);
        moep(:,1)=Kininit;

        A = obj.facetarea;
        vf = obj.vf;
        while true
            count=count+1;
            for i=1:nfcts
                Kintemp(i)=0;
                Kouttemp(i)=0;
                for j=1:nfcts %sum all the incoming radiation on i, originally reflected from all the other j facets ("radiation reflected on j" x "what perecentage does i take of j's vision")
                    inc=Koutold(j)*A(j)/A(i)*vf(j,i);  %[W/m2]
                    storerad(i,j)=inc;
                    Kintemp(i)=Kintemp(i)+(1-albedo(i))*inc;
                    Kouttemp(i)=Kouttemp(i)+albedo(i)*inc;
                end
                Kin(i)=Kin(i)+Kintemp(i); %add newly absorbed radiation to already existing one
                Koutnew(i)=Kouttemp(i); %save newly reflected radiation for next iteration
            end
            % Kiterin(count)=Kin(1);
            % Kiterout(count)=Koutnew(1);
            % [itermaxdiff(count),itermaxdiffloc(count)]=max(Kintemp(i));
            % [~,itermaxrelloc(count)]=max(Koutnew./Koutold);
            % max(Koutnew./Koutold)
            % if all(Koutnew./Koutold<0.01)
            %     break
            % end
            moep(:, count+1) = Kin;
            size(Kinold==0)
            disp(max((Kin(Kinold>0) - Kinold(Kinold>0)) ./ Kinold(Kinold>0)))
            if (max((Kin(Kinold>0) - Kinold(Kinold>0)) ./ Kinold(Kinold>0)) < 0.01)
                break
            end
            %    if (max(Koutnew-Koutold)<0.01 )
            %        disp('reached absolute tolerance')
            %        break
            %    end
            Kinold = Kin;
            Koutold = Koutnew; %overwrite reflected radiation with new value
        end

        preprocessing.addvar(obj, 'Kin', Kin);
        end

        function write_netsw(obj)
            fname = ['netsw.inp.' obj.expnr];
            fileID = fopen(fname, 'w');
            fprintf(fileID,'# %4s\n','net shortwave on facets [W/m2] (including reflections and diffusive)');
            fprintf(fileID,'%6d\n', obj.Kin);
            fclose(fileID);
        end
        
        function rayit_td(obj)
            %starttime_past_midnight = obj.solposstart * 3600 / 100;
            %[~, startindex] = min(abs(obj.time_past_midnight  - obj.solposstart * 3600 / 100))

            facets = obj.facets;
            walltypes = obj.walltypes;
            function [fct, wall] = loadfacets()
                %M = dlmread(['facets.inp.' num2str(expnr)],'',1,0);
                vars = {'or', 'wlid', 'blk', 'bld'};
                for n = 1:length(vars)
                    fct.(vars{n}) = facets(:, n);
                end

                % wallid  lGR   z0 [m]     z0h [m]     al [-]   em [-]   d1 [m]  d2 [m]    cp1 [J/kgK]  cp2 [J/kgK]   l1 [W/(m K)]  l2 [W/(m K)]    k1 [W/(m K)]    k1 [W/(m K)]
                vars = {'id','lGR','z0','z0h','al', 'em', 'd1', 'd2', 'cp1', 'cp2', 'l1', 'l2', 'k1', 'k2'};
                for n = 1:length(vars)
                    wall.(vars{n}) = walltypes(:, n);
                end

                % check whether all referenced wall ids have been defined
                idref = unique(fct.wlid);
                iddef = unique(wall.id');

                if length(wall.id) > length(iddef)
                    disp('ERROR: multiple definitions of a wall id')
                    return
                end
                if length(idref) > length(iddef)
                    disp('ERROR: more walltypes used than defined')
                    return
                end


                % assign the wall type to each of the facets -- this way it is easy to know
                % in which structure the proporties can be found.
                fct.wltp = zeros(size(fct.or));

                for m = 1:size(fct.wlid, 1)
                    if fct.wlid(m) >= 0 && fct.wlid(m) <= 10 && any(ismember(iddef, fct.wlid(m)))
                        fct.wltp(n) = 1; %normal
                    elseif fct.wlid(m) > 10 && any(ismember(iddef, fct.wlid(m)))
                        fct.wltp(n) = 2; %Green roof
                    elseif fct.wlid(m) < 0 && any(ismember(iddef, fct.wlid(m)))
                        if fct.wlid(m) <= -99
                            fct.wltp(n) = 3; %bounding wall
                        else
                            fct.wltp(n) = 4; %floor
                        end
                    else
                        disp(['ERROR: walltype of facet ',num2str(m),' not defined'])
                    end
                end


                for n = 1:length(idref)
                    if (ismember(idref(n), wall.id))
                        fct.wltp(fct.wlid == idref(n)) = 1;
                    elseif (ismember(idref(n), gr.id))
                        fct.wltp(fct.wlid == idref(n)) = 2;
                    else
                        sprintf('ERROR. wall id %8d not defined. Correct wall and greenroof property input files', idref(n))
                        return
                    end
                end
            end

            [fct, wall] = loadfacets();


            nfcts = obj.nfcts;
            albedo = zeros(nfcts,1);
            emissivity = zeros(nfcts,1);


            for i=1:nfcts
                j=find(fct.wlid(i)==wall.id);
                albedo(i)=wall.al(j);
                emissivity(i)=wall.em(j);
            end


            Sdir_td = obj.Sdir_td;
            Dsk = obj.Dsk;
            svf = obj.svf;

            %D_array = obj.ashraeC * Sdir_array;

            Kin_td = zeros(nfcts, obj.nweatherdata);
            %Dsky_array = zeros(nfcts, npos);

%             for i = 1:nfcts
%                Dsky_array(i, :) = D_array(i, :) * svf(i);
%             end

            timedepzenith_shift = circshift(obj.timedepzenith, -obj.starttime);
            timedepDsky_shift = circshift(obj.timedepDsky, -obj.starttime);
            
            parfor p = 1:obj.nweatherdata
                p
                solarzenith  = timedepzenith_shift(p);

                if solarzenith > 89
                    continue;
                end

                Sdir = Sdir_td(:, p);
                %D = D_array(:, p);
                %Dsky = Dsky_array(:, p);
                %Dsky = obj.ashraeC * obj.In_array(p) * svf;
                Dsky = timedepDsky_shift(p) * svf;
                %Kininit = (1 - albedo) .* (Sdir + D .* svf); % thesis -> Knet = Kin - Kout =  (1 - al) * Kin
                Kininit = (1 - albedo) .* (Sdir + Dsky);
                %Koutinit = albedo .* (Sdir + D .* svf); % thesis -> Kout = al * Kin
                Koutinit = albedo .* (Sdir + Dsky);


                Kin = Kininit; %Kin §adds up
                Kinold = Kininit;
                Koutold = Koutinit; %Kout is wiped clean with every reflection
                Koutnew = zeros(nfcts,1);
                Kintemp = zeros(nfcts,1);
                Kouttemp = zeros(nfcts,1);
                count = 0;
                storerad = zeros(nfcts,nfcts);

                moep = zeros(nfcts,20);
                moep(:,1) = Kininit;

                %A = obj.facetarea;
                vf = obj.vf;
                while true
                    count = count+1;
                    for i = 1:nfcts
                        Kintemp(i) = 0;
                        Kouttemp(i) = 0;
                        for j = 1:nfcts %sum all the incoming radiation on i, originally reflected from all the other j facets ("radiation reflected on j" x "what perecentage does i take of j's vision")
                            %inc = Koutold(j) * A(j) / A(i) * vf(j,i);  %[W/m2]
                            inc = Koutold(j) * vf(i,j);  %[W/m2]
                            storerad(i,j) = inc;
                            Kintemp(i) = Kintemp(i) + (1 - albedo(i)) * inc;
                            Kouttemp(i) = Kouttemp(i) + albedo(i) * inc;
                        end
                        Kin(i) = Kin(i) + Kintemp(i); %add newly absorbed radiation to already existing one
                        Koutnew(i) = Kouttemp(i); %save newly reflected radiation for next iteration
                    end

                    moep(:, count+1) = Kin;
                    max((Kin - Kinold) ./ Kinold)

%                     if (Kinold == 0)
%                         find(Kinold == 0)
%                         return
%                     end

                    if (max((Kin - Kinold) ./ Kinold) < 0.01)
                        break
                    end
                    %    if (max(Koutnew-Koutold)<0.01 )
                    %        disp('reached absolute tolerance')
                    %        break
                    %    end
                    Kinold = Kin;
                    Koutold = Koutnew; %overwrite reflected radiation with new value
                    Kin_td(:, p) = Kin;
                end

            end


            preprocessing.addvar(obj, 'Kin_td', Kin_td)
            %preprocessing.addvar(obj, 'Dsky_td', Dsky_array)
            preprocessing.addvar(obj, 'emissivity', emissivity);
        end

        function write_timedepsw(obj)
            % Interpolate onto energy balance timestep
            %timeEB = 0:obj.dtEB:(86400 - obj.dtEB);
            timeEB = 0:obj.dtEB:86400;
            %Kin_td_interp_p = interp1(0:3600:(86400-3600), obj.Kin_td', timeEB);
            Kin_td_interp_p = interp1(0:3600:86400, [obj.Kin_td, obj.Kin_td(:,1)]', timeEB, 'makima');
            %obj.Kin_td_interp = Kin_td_interp_p(1:end-1,:)';
            %preprocessing.addvar(obj, 'Kin_td_interp', Kin_td_interp_p')
                      
            fname = ['timedepsw3.inp.' obj.expnr];
            fileID = fopen(fname, 'w');
            fprintf(fileID,'# %4s\n','net shortwave on facets [W/m2] (including reflections and diffusive)');
            for n = 1:obj.ntimedepsw
                n
                %t = (n-1) * 86400 / (obj.nweatherdata);
                t = timeEB(n) + 172800;
                fprintf(fileID, ['#    ', num2str(t), '\n']);
                fclose(fileID);
                %dlmwrite(fname, obj.Kin_td(:,n)', '-append','delimiter',' ','precision','%11.6f')
                dlmwrite(fname, Kin_td_interp_p(n,:), '-append','delimiter',' ','precision','%11.6f')
                fileID = fopen(fname, 'a');
            end
            %t = 86400;
            %fprintf(fileID, ['#    ', num2str(t), '\n']);
            fclose(fileID);
        end
        
        function write_timedeplw(obj)
            fname = ['timedeplw.inp.' obj.expnr];
            fileID = fopen(fname, 'w');
            fprintf(fileID,'# %4s\n','sky longwave on facets [W/m2]');
            fprintf(fileID,'# %4s\n','time skyLW');
            timedeplw_shift = circshift(obj.timedepLsky, -obj.starttime);
            fclose(fileID);
            dlmwrite(fname, [obj.timedeptime, timedeplw_shift], 'delimiter', ' ', '-append')
            %t = 86400;
            %fprintf(fileID, ['#    ', num2str(t), '\n']);
            
        end
        

        function generate_Tfacinit(obj, iss)
            %% Inital Temperature and longwave
            %if iss = true, solve energy budget equation in an initial steady state
            %assume 0 wall heatflux
            %assume 0 latent heatflux
            %assume constant air temperature and heat transfer coefficient
            %assume constant longwave
            %solve for Tinit and Linit
            %assume absorbtivity is equal to emissivity
            %=>
            %K+L=H
            %iss=0;
            if iss
                nfcts = obj.nfcts;

                Tair = 300;  %K, air temperature
                Tinitial = 300; %K, initial facet temperature

                % Tinc=1*ones(nfcts,1); %incremental temperature change of facets if not in equilibrium
                % tolerance=2.5;  %W/m2, if below this threshold, change will be made to facet temperature

                Tinc = 2 * ones(nfcts, 1); %incremental temperature change of facets if not in equilibrium
                tolerance = 2.5 * Tinc(1);  %W/m2, if below this threshold, change will be made to facet temperature
                %the two are somewhat related, an increase of 1K will result in a longwave
                %change of approximately 5W/m2, e.g.:
                %if we are 4.56K off we correct approximately  (1-5/22.8)*2 = 1.5616
                %now we are 3K off and correct approximately for (1-5/15)*2 = 1.3333
                %now we are 1.6666 off and we correct for 0.8 (had we corrected for 2 at
                %this step, we would have overshot and the solution would oscilate)
                %now we are satisfied at 0.866 off which is less than 1degree
                %if for any


                Tterminate = 0; %Terminate if the absolute temperature change between iterations is equal or below this value
                %somewhat redundant since this basically means all the facets are within the tolerance

                vf = obj.vf;
                svf = obj.svf;
                facetarea = obj.facetarea;
                Kin = obj.Kin;
                emissivity = obj.emissivity;
                absorptivity = obj.emissivity;
                sigma = 5.67e-8;
                Lsk = obj.skyLW;
                hc = 0;  %(rho*cp)/R    ~(1.2*1000)/100=12   R~100

                Tinit = ones(nfcts, 1) * Tinitial;
                Lsky = Lsk .* svf;

                %Loutinit=emissivity.*sigma.*Tinit.^4;
                Told = Tinit;
                Tnew = zeros(nfcts, 100);
                Lin = zeros(nfcts, 1);
                %b
                k = 0;
                while true  %what happens to the reflected longwave? i.e. (1-absorptivity)*Lin
                    %for k=1:10
                    %count=count+1;
                    k = k + 1;
                    %maxchange=0;
                    change = 0;
                    Lin(:) = 0;
                    for i = 1:nfcts
                        for j = 1:nfcts %sum all the incoming radiation on i, originally reflected from all the other j facets ("radiation reflected on j" x "what perecentage does i take of j's vision")
                            inc = Told(j)^4 * emissivity(j) * sigma * vf(j, i) * facetarea(j) / facetarea(i);
                            Lin(i) = Lin(i) + inc;
                        end

                        %calculate energy balance
                        eb = Kin(i) + absorptivity(i) * (Lin(i) + Lsky(i)) - hc * (Told(i) - Tair) - emissivity(i) * sigma * Told(i)^4;

                        if eb < -tolerance %if energy balance is negative, facet is too hot

                            Tnew(i,k) = Told(i) - Tinc(i) * (1 + tolerance / eb); %remove incremental temperature, scale by deviation from accepted tolerance in W/m2

                        elseif eb > tolerance
                            Tnew(i, k) = Told(i) + Tinc(i) * (1 - tolerance / eb);
                        else %do nothing, within tolerance
                            Tnew(i, k) = Told(i);
                        end

                        % maxchange=max(maxchange,abs(Tnew(i,k)-Told(i)));
                        change = change + abs(Tnew(i, k) - Told(i));

                    end

                    Told = Tnew(:,k);
                    %plot(1:nfcts,Tnew(:,k),'-x')

                    if change <= Tterminate %maxchange<=Tterminate
                        break
                    elseif k > 99
                        break
                    end
                end
                preprocessing.addvar(obj, 'Tfacinit', Tnew(:, k));
            else
                preprocessing.addvar(obj, 'Tfacinit', obj.facT * ones(obj.nfcts, 1));
            end
%             if ltestplot
%                 figure
%                 hold on
%                 plot(Tnew(sorti,k)-Tinit(sorti),'r-x')
%                 plot(blub(sorti,k)-Tinit(sorti),'b-x')
%                 xlabel('facet NR')
%                 ylabel('T_{end}-T_{init}')
%                 ax1 = gca; % current axes
%                 set(ax1,'Xtick',1:1:nfcts)
%                 xticklabels(sortt)
%                 ax1_pos = ax1.Position; % position of first axes
%                 xlabel('orientation')
%                 ylabel('T_{new}-T_{init}')
%                 ax2 = axes('Position',ax1_pos,...
%                     'YAxisLocation', 'right',...
%                     'XAxisLocation', 'top', ...
%                     'Color', 'None');
%                 set (ax2, 'XLim', get (ax1, 'XLim'), 'Layer', 'top');
%                 set (ax2, 'YLim', get (ax1, 'YLim'), 'Layer', 'top');
%                 set(ax2,'Xtick',1:1:nfcts)
%                 set(ax2,'Xticklabels',sorti)
%                 set(ax2,'FontSize',8)
%                 xlabel('facet index')
%
%
%
%                 figure
%                 hold on
%                 plot(Tnew(15,1:k),'-x')
%                 plot(Tnew(41,1:k),'-x')
%                 plot(Tnew(51,1:k),'-x')
%                 plot(Tnew(58,1:k),'-x')
%                 plot(Tnew(68,1:k),'-x')
%                 xlabel('iteration')
%                 ylabel('facet temperature')
%                 legend('15','41','51','58','68')
%             end
        end

        function write_Tfacinit(obj)
            fname = ['Tfacinit.inp.', obj.expnr];
            fileID = fopen(fname, 'W');
            fprintf(fileID, '# Initial facet tempereatures in radiative equilibrium\n');
            fclose(fileID);
            dlmwrite(fname, obj.Tfacinit, '-append','delimiter',' ','precision','%4f')

%             fname = ['Tfacinitnudged.inp.' obj.expnr];
%             fileID = fopen(fname, 'W');
%             fprintf(fileID, '# Initial facet tempereatures in radiative equilibrium nudged to 288\n');
%             fclose(fileID);
%             blub = Tnew(Tnew > 0);
%             blublub = mean(blub(:));
%             blub = Tnew - (Tnew - blublub) * 0.5;
%             blub(Tnew == 0) = 0;
%             dlmwrite(fname, blub, '-append', 'delimiter',' ','precision','%4f')
        end

        function write_blocks(obj)
            blocks = fopen( ['blocks.inp.' obj.expnr], 'w');
            fprintf(blocks, '# Block data\n');
            fprintf(blocks, '#  il\t   iu\t   jl\t   ju\t   kl\t   ku\t dtop\t dwest\t deast\t dnor\t dsou dbot\n');
            fprintf(blocks, '%4d\t%4d\t%4d\t%4d\t%4d\t%4d\t%4d\t%4d\t%4d\t%4d\t%4d\t%4d\n', obj.blocks');
            fclose(blocks);
        end

        function write_facets(obj)
            fname = ['facets.inp.' num2str(obj.expnr)];
            fileID = fopen(fname,'w');
            fprintf(fileID,'# %4s %6s %6s %6s\n','or', 'wl', 'blk', 'bld');
            fprintf(fileID,'%6d %6d %6d %6d\n', obj.facets(:, 1:4)');
            fclose(fileID);
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

        function [coaa, cobb, pf1sf2u] = slice(fi,fi2,coa,cob,cornerweight,centerweight)
            %in case one facet cuts through
            %the plane of the other facet:
            %                    |          .-----------.
            %   ---------------  |    .or.  |           |
            %                    |          | . . . . . |.------.
            %                    |          |___________||      |
            %                                            |______|
            %
            % if this is the case slice facet up (max 2 slices) -> coaa/cobb
            % count the number of slices -> m
            % increase the percentage they see of each other, since the two corners were blocked but now
            % they are not included in the vf-calculation anymore -> pf1sf2u


            %1=top,2=west face,3=east face,4=north face, 5=south face
            pf1sf2u = 0;
            coaa = coa(6:9, :);
            cobb = cob(6:9, :);
            nottested = 1; %necessary to check for i and j (alternatively the if's in each case could be combined)
            if ((fi*fi2 == 6) || (fi * fi2 == 20)) %opposite each other, do nothing;
                return
            elseif (fi == 1) %fi2=2,3,4,5

                if (cobb(1, 3) < coaa(1, 3)) %j starts lower than height of i
                    cobb(1, 3) = coaa(1, 3); %slice j on same height as i
                    cobb(4, 3) = coaa(1, 3); %slice j on same height as i
                    pf1sf2u = 2 * cornerweight; %add 2*cornerweight since lower corner now is visible
                    nottested = 0;
                    if  cob(1,3) <= coaa(1, 3) %if center of j is lower than height also add centerweight
                        pf1sf2u = 2 * cornerweight + centerweight;
                    end
                end

            elseif (fi == 2) %fi2=1,4,5
                if (((fi2 == 1) && (coaa(1, 1) < cobb(4, 1))) || (((fi2 == 4) || (fi2 == 5))&& (coaa(1, 1) < cobb(4, 1))))%my x<xmax
                    cobb(3, 1) = coaa(1, 1);
                    cobb(4, 1) = coaa(1, 1);
                    pf1sf2u = 2 * cornerweight;
                    nottested = 0;
                    if coaa(1, 1) <= cob(1, 1) %if j center is east of i wall
                        pf1sf2u = 2 * cornerweight + centerweight;
                    end
                    % elseif (((fi2==4) || (fi2==5))&& (coa(1,1)<cob(4,1)))
                    %     cobb=cob(2:5,:);
                    %     cobb(3,1)=coaa(1,1);
                    %     cobb(4,1)=coaa(1,1);
                    %     pf1sf2u=2*cornerweight;
                end
            elseif (fi == 3) %fi2=1,4,5
                if (((fi2 == 1) && (coaa(1, 1) > cobb(1, 1))) || (((fi2 == 4) || (fi2 == 5)) && (coaa(1, 1) > cobb(1, 1)))) %my x>xmin
                    cobb(1, 1) = coaa(1, 1);
                    cobb(2, 1) = coaa(1, 1);
                    pf1sf2u = 2 * cornerweight;
                    nottested = 0;
                    if coaa(1, 1) >= cob(1, 1) %if j center is west of i wall
                        pf1sf2u = 2 * cornerweight + centerweight  ;
                    end
                    % elseif (((fi2==4) || (fi2==5)) && (coa(1,1)>cob(1,1)))
                    %     cobb=cob(2:5,:);
                    %     cobb(1,1)=coaa(1,1);
                    %     cobb(2,1)=coaa(1,1);
                    %     pf1sf2u=2*cornerweight;
                end

            elseif (fi == 4) %fi2=1,2,3
                if ((fi2 == 1) && (coaa(1, 2) > cobb(1, 2))) %my y>ymin
                    cobb(1, 2) = coaa(1, 2);
                    cobb(4, 2) = coaa(1, 2);
                    pf1sf2u = 2 * cornerweight;
                    nottested = 0;
                    if coaa(1, 2) >= cob(1, 2) %if j center is south of i wall
                        pf1sf2u = 2 * cornerweight + centerweight;
                    end
                elseif (((fi2 == 2) || (fi2 == 3)) && (coaa(1, 2) > cobb(1, 2)))

                    cobb(1, 2) = coaa(1, 2);
                    cobb(2, 2) = coaa(1, 2);
                    pf1sf2u = 2 * cornerweight;
                    nottested = 0;
                    if coaa(1, 2) >= cob(1, 2) %if j center is south of i wall
                        pf1sf2u = 2 * cornerweight + centerweight;
                    end
                end

            elseif (fi == 5) %fi2=1,2,3
                if ((fi2 == 1) && (coaa(1, 2) < cobb(2, 2))) %my y<ymax
                    cobb(2, 2) = coaa(1, 2);
                    cobb(3, 2) = coaa(1, 2);
                    pf1sf2u = 2 * cornerweight;
                    nottested = 0;
                    if coaa(1, 2) <= cob(1, 2) %if j center is north of i wall
                        pf1sf2u = 2 * cornerweight + centerweight;
                    end
                elseif (((fi2 == 2) || (fi2 == 3)) && (coaa(1, 2) < cobb(3, 2))) %ILS13 11.11.17, was cobb(2,2)
                    cobb(3,2) = coaa(1,2);
                    cobb(4,2) = coaa(1,2);
                    pf1sf2u = 2 * cornerweight;
                    nottested = 0;
                    if coaa(1, 2) <= cob(1, 2) %if j center is north of i wall
                        pf1sf2u = 2 * cornerweight + centerweight;
                    end
                end
            end
            if nottested
                if (fi2 == 1)  %fi=2,3,4,5
                    if (coaa(1, 3) < cobb(1, 3)) %i starts lower than height of j
                        coaa(1, 3) = cobb(1, 3); %slice i on same height as j
                        coaa(4, 3) = cobb(1, 3); %slice i on same height as j
                        pf1sf2u = 2 * cornerweight; %add 2*cornerweight since lower corner now is visible
                        if  coa(1,3) <= cobb(1, 3) %if center of i is lower than height also add centerweight
                            pf1sf2u = 2 * cornerweight + centerweight;
                        end
                    end
                    %end

                elseif (fi2 == 2) %fi=1,4,5
                    if (((fi == 1) && (cobb(1, 1) < coaa(4, 1))) || (((fi == 4) || (fi == 5))&& (cobb(1, 1)<coaa(4, 1))))%my x<xmax
                        coaa(3, 1)=cobb(1, 1);
                        coaa(4, 1)=cobb(1, 1);
                        pf1sf2u = 2 * cornerweight;

                        if cobb(1, 1) <= coa(1, 1) %if i center is east of j wall
                            pf1sf2u = 2 * cornerweight + centerweight;
                        end

                    end

                elseif (fi2 == 3) %fi=1,4,5
                    if (((fi == 1) && (cobb(1, 1) > coaa(1, 1))) || (((fi == 4) || (fi == 5)) && (cobb(1, 1) > coaa(1, 1)))) %my x>xmin
                        coaa(1, 1) = cobb(1, 1);
                        coaa(2, 1) = cobb(1, 1);
                        pf1sf2u = 2 * cornerweight;
                        if cobb(1, 1) >= coa(1, 1) %if i center is west of j wall
                            pf1sf2u = 2 * cornerweight + centerweight;
                        end
                    end

                elseif (fi2 == 4) %fi=1,2,3
                    if ((fi == 1) && (cobb(1, 2) > coaa(1, 2))) %my y>ymin
                        coaa(1, 2) = cobb(1, 2);
                        coaa(4, 2) = cobb(1, 2);
                        pf1sf2u = 2 * cornerweight;
                        if cobb(1, 2) >= coa(1, 2) %if i center is south of j wall
                            pf1sf2u = 2 * cornerweight + centerweight;
                        end
                    elseif (((fi == 2) || (fi == 3)) && (cobb(1, 2) > coaa(1, 2)))
                        coaa(1, 2) = cobb(1, 2);
                        coaa(2, 2) = cobb(1, 2);
                        pf1sf2u = 2 * cornerweight;
                        if cobb(1,2) >= coa(1,2) %if i center is south of j wall
                            pf1sf2u = 2 * cornerweight + centerweight  ;
                        end
                    end
                elseif (fi2 == 5)  %fi=1,2,3
                    if ((fi == 1) && (cobb(1, 2) < coaa(2, 2))) %my y<ymax
                        coaa(2, 2) = cobb(1,2);
                        coaa(3,2) = cobb(1,2);
                        pf1sf2u = 2 * cornerweight;
                        if cobb(1, 2) <= coa(1, 2) %if j center is north of i wall
                            pf1sf2u = 2 * cornerweight + centerweight;
                        end
                    elseif (((fi == 2) || (fi == 3)) && (cobb(1, 2) < coaa(2, 2)))
                        coaa(3, 2) = cobb(1, 2);
                        coaa(4, 2) = cobb(1, 2);
                        pf1sf2u = 2 * cornerweight;
                        if cob(1, 2) <= cob(1, 2) %if j center is north of i wall
                            pf1sf2u = 2 * cornerweight + centerweight;
                        end
                    end
                end
            end
        end

        function [F12, F21, A1, A2] = ViewFactor(C1, C2, A1, A2, GP, vcorner)
            %ViewFactor calculates view factors between two planar surfaces in 3D.
            %   [F12,F21,A1,A2] = ViewFactor(C1,C2) returns the view factors between
            %   planer surfaces whose vertices are stored in coordinate arrays
            %   C1 and C2.  Cn are of the form [x1,y1,z1; x2,y2,z2; x3,y3,z3;....],
            %   where (xn,yn,zn) is the coordinate of the nth vertex.  The vertices
            %   should be in order as  encountered on a trip around the perimeter, and
            %   each vertex should be counted only once. The integration is carried
            %   out with 7-point Gauss-Legendre quadrature.

            %
            % adapted from  Matt Fig, Date: 09/20/2016

            if nargin < 3
                GP = 7;
            end
            [A, W] = Gauss(GP); % Weights and abscissa for Gauss-Legendre quadrature

            if vcorner(1) == 0
                L1 = 4; %size(C1,1);   % Number of vertices.
                L2=[4 4 4 4];
                %L2=[size(C2,1) size(C2,1) size(C2,1) size(C2,1)];
                % Close the boundary.
                C1(5,:) = C1(1,:);
                C2(5,:) = C2(1,:);
                cc = 0;
            else
                dorder=[1 2 3 4];
                L1 = 4; %size(C1,1);   % Number of vertices.
                L2=[4 4 4 3];
                order1=circshift(dorder,4-vcorner(1));
                order2=circshift(dorder,4-vcorner(2));
                C1=C1(order1,:);
                C2=C2(order2,:);
                C1(5,:) = C1(1,:);
                C2(5,:) = C2(1,:);
                Lc=norm(C1(4,:)-C1(1,:));
                cc = Lc^2*(1.5-log(Lc))*4*sign(vcorner(3));
            end

            % A1 = rectarea3d(C1);
            % A2 = rectarea3d(C2);
            % A1 = polyarea3d(C1);
            % A2 = polyarea3d(C2);
            S=0;

            for ii = 1:L1 % Loop over segment pairs.
                P1 = C1(ii,:);
                P2 = C1(ii+1,:);

                for jj = 1:L2(ii)
                    SM = 0;
                    P3 = C2(jj,:);
                    P4 = C2(jj+1,:);
                    % Next perform the Gauss-Legendre quadrature
                    for kk = 1:GP
                        SM = SM + sum(W(kk)*W.*F(A(kk),A,P1,P2,P3,P4));
                    end

                    S = S + SM;
                end
            end

            %Calculation of the view factors
            %cp = cross(C1(2,:) - C1(1,:),C1(3,:) - C1(1,:));
            %A1 = polyarea3d(C1,cp/(norm(cp)));
            %cp = cross(C2(2,:)-C2(1,:),C2(3,:)-C2(1,:));
            %A2 = polyarea3d(C2,cp/(norm(cp)));

            F12 = (abs(S+cc))/(8*pi*A1);
            F21 = (abs(S+cc))/(8*pi*A2);

            function INT  = F(s,t,P1,P2,P3,P4)
                % Integrand for contour integral.  We parameterize the linesegments on
                % s,t=[-1,1] for easy 2D Gaussian quadrature..
                % Parametric equations for the first line segment.
                t_ = (t+1)/2;
                x1 = P1(1)+(P2(1)-P1(1))*t_;
                y1 = P1(2)+(P2(2)-P1(2))*t_;
                z1 = P1(3)+(P2(3)-P1(3))*t_;
                % Parametric equations for the second line segment.
                s_ = (s+1)/2;
                x2 = P3(1)+(P4(1)-P3(1))*s_;
                y2 = P3(2)+(P4(2)-P3(2))*s_;
                z2 = P3(3)+(P4(3)-P3(3))*s_;
                % Distance between seg 1 and seg 2.
                R = sqrt((x2-x1).^2 + (y2-y1).^2 + (z2-z1).^2);
                % The integrand
                INT = log(R)*((P2(1)-P1(1))*(P4(1)-P3(1))+(P2(2)-P1(2))*...
                    (P4(2)-P3(2))+(P2(3)-P1(3))*(P4(3)-P3(3)));


                % function [area] = polyarea3d(V,N)
                % % Calculates the area of the polygon given by the matrix of vertices and
                % % the normal vector N.  This assumes the vertices are closed, i.e., the
                % % last row in V is equal to the first row.  Algorithm developed by
                % % Daniel Sunday and available in c++ here:
                % % http://geomalgorithms.com/a01-_area.html
                % n = size(V,1)-1;
                % % select largest abs coordinate to ignore for projection
                % ax = abs(N(1));
                % ay = abs(N(2));
                % az = abs(N(3));
                %
                % if (ax > ay)&&(ax > az)
                %     area = (sum(V(2:n,2) .* (V(3:n+1,3) - V(1:n-1,3))) + ...
                %            (V(n+1,2) .* (V(2,3) - V(n,3))))/ (2 * N(1));
                % elseif (ay > az)&&(ay > ax)
                %     area = (sum((V(2:n,3) .* (V(3:n+1,1) - V(1:n-1,1)))) + ...
                %            (V(n+1,3) * (V(2,1) - V(n,1))))/ (2 * N(2));
                % else
                %     area = (sum((V(2:n,1) .* (V(3:n+1,2) - V(1:n-1,2)))) + ...
                %            (V(n+1,1) * (V(2,2) - V(n,2))))/(2 * N(3));
                % end

                % function [area] = rectarea3d(V)
                % % if only rectangles are used
                % %area=(0.2+norm(V(1,:)-V(2,:)))*(0.2+norm(V(2,:)-V(3,:)));
                % l1=norm(V(1,:)-V(2,:));
                % l2=norm(V(2,:)-V(3,:));
                % area=l1*l2;
            end

            function [abscissa, weights] = Gauss(n)
                % Generates the abscissa and weights for a Gauss-Legendre quadrature.
                % Reference:  Numerical Recipes in Fortran 77, Cornell press.
                abscissa = zeros(n,1);  % Preallocations.
                weights = abscissa;
                m = (n+1)/2;
                for iii=1:m
                    z = cos(pi*(iii-.25)/(n+.5)); % Initial estimate.
                    z1 = z+1;

                    while abs(z-z1)>eps
                        p1 = 1;
                        p2 = 0;

                        for jjj = 1:n
                            p3 = p2;
                            p2 = p1;
                            p1 = ((2*jjj-1)*z*p2-(jjj-1)*p3)/jjj; % The Legendre polynomial.
                        end

                        pp = n*(z*p1-p2)/(z^2-1);   % The L.P. derivative.
                        z1 = z;
                        z = z1-p1/pp;
                    end

                    abscissa(iii) = -z;      % Build up the abscissas.
                    abscissa(n+1-iii) = z;
                    weights(iii) = 2/((1-z^2)*(pp^2));  % Build up the weights.
                    weights(n+1-iii) = weights(iii);
                end
            end
        end

        function [ prcntgblckd ] = prblckd(i,j,coa,ndima,sun,v1,cob,ndimb,facets,centerweight,cornerweight,nblocks,nbw,bl)
            %prcntgblckd     = percentage of view blocked
            %i               = index of facet 1
            %j               = index of facet 2      (only if sun==false)
            %coa             = corners of facet 1
            %ndima           = dimension of facet 1
            %sun             = calculation if facet can see sun, not if two facets can
            %                  see each other
            %v1              = vector to the sun
            %cob             = corners of facet 2    (only if sun==false)
            %ndimb           = dimension of facet 2  (only if sun==false)
            %facets               = facets
            %centerweight    = % of viewfield if centers see each other
            %cornerweight    = % per corner if they see each other
            %nblocks         = number of blocks
            %nbw             = number of bounding walls
            %bl              = list of blocks
            testcrit = 0; %=4;
            prcntgblckd = 0;
            flag = 0;
            if sun  %calculation between facet 1 and the sun
                bi=facets(i,3); %block index
                ci=facets(i,4); %building index (-1 for roads, -101 for bounding wall)
                if ndima>testcrit  %also test corner, otherwise only test center
                    for k=1:5
                        facetpoint=coa(k,:);
                        for n=1:nblocks %check if any block intersects
                            if ((n==bi) && (ci>0))  %don't intersect with itself
                                continue
                            end

                            [flag,dint] = preprocessing.rbi(facetpoint, v1, bl(n,:));

                            intersection = facetpoint + dint*v1;
                            if intersection(3)<facetpoint(3) %downstream direction of sun, thus not blocking the sun
                                flag=0;
                            end
                            if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                                if k==1
                                    prcntgblckd=centerweight;
                                else
                                    prcntgblckd=prcntgblckd+cornerweight;
                                end
                                break %get out of inner for loop
                            end
                        end

                        if ~flag
                            for m=1:nbw %check if any bounding wall intersects

                                [flag,dint] = preprocessing.rbi(facetpoint, v1, bl(m+nblocks,:));

                                intersection = facetpoint + dint*v1;

                                if intersection(3)<facetpoint(3) %downstream direction of sun
                                    flag=0;
                                end
                                if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                                    if k==1 %it's the center of the facet
                                        prcntgblckd=centerweight;
                                    else
                                        prcntgblckd=prcntgblckd+cornerweight;
                                    end
                                    break %get out of inner for loop
                                end
                            end
                        end
                    end
                else
                    % %disp('small area, checking only center')
                    % %disp(['facet: ' num2str(i)])
                    facetpoint=coa(1,:);  %center of facet
                    for n=1:nblocks %check if any block intersects
                        if ((n==bi) && (ci>0))  %don't intersect with itself
                            continue
                        end

                        [flag,dint] = preprocessing.rbi(facetpoint, v1, bl(n,:));

                        intersection = facetpoint + dint*v1;
                        if intersection(3)<facetpoint(3) %downstream direction of sun, thus not blocking the sun
                            % intersectionss(:)=0;
                            flag=0;
                        end
                        if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                            prcntgblckd=1;
                            break %get out of inner for loop
                        end
                    end

                    if ~flag
                        for k=1:nbw %check if any bounding wall intersects

                            [flag,dint] = preprocessing.rbi(facetpoint, v1, bl(k+nblocks,:));

                            intersection = facetpoint + dint*v1;

                            if intersection(3)<facetpoint(3) %downstream direction of sun
                                flag=0;
                            end
                            if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                                prcntgblckd=1;
                                break %get out of inner for loop
                            end
                        end
                    end

                end
            else %calculation between facet 1 and facet 2
                %disp('facet to facet')
                %disp(['% blocked so far ' num2str(prcntgblckd)])

                bi=facets(i,3); %block index
                ci=facets(i,4); %building index (-1 for roads, -99 for bounding wall)
                bi2=facets(j,3); %block index
                fi=facets(i,1); %facet index
                fi2=facets(j,1); %facet index

                if (ndima+ndimb)>(2*testcrit)  %also test corner, otherwise only test center
                    %disp('test corners')
                    %disp(['% blocked so far ' num2str(prcntgblckd)])
                    %find orientation to know which corners to compare
                    ordera=[1 2 3 4 5];
                    cas=fi*10+fi2;
                    %disp('case')
                    switch cas
                        case 12, orderb=[1 3 4 5 2];
                        case 13, orderb=[1 2 5 4 3];
                        case 15, orderb=[1 3 2 5 4];
                        case 21, orderb=[1 5 2 3 4];
                        case 24, orderb=[1 5 4 3 2];
                        case 31, orderb=[1 2 5 4 3];
                        case 35, orderb=[1 5 4 3 2];
                        case 42, orderb=[1 5 4 3 2];
                        case 51, orderb=[1 3 2 5 4];
                        case 53, orderb=[1 5 4 3 2];
                        otherwise %14, 41, 23, 32, 25, 52, 34, 43, 45, 54
                            orderb=[1 2 3 4 5];
                    end

                    for k=1:5
                        test=1;
                        oa=ordera(k);
                        ob=orderb(k);
                        facetpoint=coa(oa,:);
                        facetpoint2=cob(ob,:);
                        v=[facetpoint2(1)-facetpoint(1),facetpoint2(2)-facetpoint(2),facetpoint2(3)-facetpoint(3)];
                        dint=norm(v);
                        v1=v/dint;
                        if ((fi==1) && (v(3)<=0)) %horizontal facets can't see facets whos center is lower
                            test=0;
                        elseif ((fi==2) && (v(1)>=0)) %west facets cant see anything to the east
                            test=0;
                        elseif ((fi==3) && (v(1)<=0)) %east facets cant see anything to the west
                            test=0;
                        elseif ((fi==4) && (v(2)<=0)) %north facets cant see anything to the south
                            test=0;
                        elseif ((fi==5) && (v(2)>=0)) %south facets cant see anything to the north
                            test=0;
                        elseif ((fi~=1) && (fi2==1) && (v(3)>0)) %vertical walls cant see any horizontal facet that is higher
                            test=0;
                            % the following cases should be blocked by a block anyway, but
                            % this is faster
                        elseif ((fi==4) || (fi==5)) && (fi2==3) && (v(1)>0) %north/south can't see an east facet if it is east of it
                            test=0;
                        elseif ((fi==4) || (fi==5)) && (fi2==2) && (v(1)<0) %north/south can't see a west facet if it is west of it
                            test=0;
                        elseif ((fi==2) || (fi==3)) && (fi2==4) && (v(2)>0) %west/east can't see a north facet if it is north of it
                            test=0;
                        elseif ((fi==2) || (fi==3)) && (fi2==5) && (v(2)<0) %west/east can't see a south facet if it is south of it
                            test=0;
                        end
                        if test
                            %                 disp(['can potentially see each other' num2str(k)])
                            %                 disp(['% blocked so far ' num2str(prcntgblckd)])
                            for n=1:nblocks %check if any block intersects
                                if ((n==bi) && (ci>0))  %don't intersect with itself
                                    continue
                                end

                                [flag,dints] = preprocessing.rbi(facetpoint, v1, bl(n,:));

                                if dints<=(0+eps)  %intersection is downstream
                                    flag=0;
                                elseif dints>=(dint-eps) %block is behind target subfacet
                                    flag=0;
                                end
                                if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                                    if k==1
                                        prcntgblckd=centerweight;
                                    else
                                        %disp(['blocked by ' num2str(n)])
                                        prcntgblckd=prcntgblckd+cornerweight;
                                    end
                                    break %get out of inner for loop
                                end
                            end
                        else %facets can't possibly see each other
                            %disp('facets cant possibly see each other')
                            if k==1
                                prcntgblckd=centerweight;
                            else
                                prcntgblckd=prcntgblckd+cornerweight;
                            end
                        end
                    end

                else %small area, checking only center
                    % %disp('small area, checkig only center')
                    test=1;
                    facetpoint=coa(1,:);
                    facetpoint2=cob(1,:);
                    v=[facetpoint2(1)-facetpoint(1),facetpoint2(2)-facetpoint(2),facetpoint2(3)-facetpoint(3)];
                    dint=norm(v);
                    v1=v/dint;

                    %abs(v)
                    %%disp(['v1: ' num2str(v1)])
                    %test if facets can potentially see each other
                    if ((fi==1) && (v(3)<=0)) %horizontal facets can't see facets whos center is lower
                        test=0;
                    elseif ((fi==2) && (v(1)>=0)) %west facets cant see anything to the east
                        test=0;
                    elseif ((fi==3) && (v(1)<=0)) %east facets cant see anything to the west
                        test=0;
                    elseif ((fi==4) && (v(2)<=0)) %north facets cant see anything to the south
                        test=0;
                    elseif ((fi==5) && (v(2)>=0)) %south facets cant see anything to the north
                        test=0;
                    elseif ((fi~=1) && (fi2==1) && (v(3)>0)) %vertical walls cant see any horizontal facet that is higher
                        test=0;
                        % the following cases should be blocked by a block anyway, but
                        % this is faster
                    elseif ((fi==4) || (fi==5)) && (fi2==3) && (v(1)>0) %north/south can't see an east facet if it is east of it
                        test=0;
                    elseif ((fi==4) || (fi==5)) && (fi2==2) && (v(1)<0) %north/south can't see a west facet if it is west of it
                        test=0;
                    elseif ((fi==2) || (fi==3)) && (fi2==4) && (v(2)>0) %west/east can't see a north facet if it is north of it
                        test=0;
                    elseif ((fi==2) || (fi==3)) && (fi2==5) && (v(2)<0) %west/east can't see a south facet if it is south of it
                        test=0;
                    end

                    %%disp(['test: ' num2str(test)])

                    if test %facets could see each other, check if blocked
                        %calculate actual intersection distance
                        for n=1:nblocks %check if any block intersects
                            if ((n==bi) && (ci>0))  %don't intersect with itself
                                continue
                            end

                            [flag,dints] = preprocessing.rbi(facetpoint, v1, bl(n,:));

                            if dints<=0  %intersection is downstream
                                flag=0;
                            elseif dints>=dint %block is behind target subfacet
                                flag=0;
                            end
                            if flag==1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
                                %  %disp(['intersection with block: ' num2str(n)])
                                prcntgblckd=1;
                                break %get out of inner for loop
                            end
                        end

                    else %facets can't possibly see each other
                        prcntgblckd=1;
                    end
                end
            end
        end

        function [flag ,tmin] = rbi(origin, direction, vec)
            %  Ray/box intersection using the Smits' algorithm
            %
            % Input:
            %    origin.
            %    direction.
            %    vec = [xl,xu,yl,yu,zl,zu]  (=Block)
            % Output:
            %    flag: (0) Reject, (1) Intersect.
            %    tmin: distance from the ray origin.

            vmin=[vec(1) vec(3) vec(5)];
            vmax=[vec(2) vec(4) vec(6)];


            if ((direction(1) == 0) && (((vmin(1) - origin(1))==0) || ((vmax(1) - origin(1))==0))) %comparing floats!  %0/0 is undefined whereas 1/0 = inf, -1/0=-inf works
                %this will only be the case if the ray is tangential -> no intersection
                %(or infinite intersections...)
                tmin=0;
                tmax=inf;
            elseif (direction(1) >= 0)
                tmin = (vmin(1) - origin(1)) / direction(1);
                tmax = (vmax(1) - origin(1)) / direction(1);
            else
                tmin = (vmax(1) - origin(1)) / direction(1);
                tmax = (vmin(1) - origin(1)) / direction(1);
            end

            if ((direction(2) == 0) && (((vmin(2) - origin(2))==0) || ((vmax(2) - origin(2))==0))) %comparing floats!
                tymin=0;
                tymax=inf;
            elseif (direction(2) >= 0)
                tymin = (vmin(2) - origin(2)) / direction(2);
                tymax = (vmax(2) - origin(2)) / direction(2);
            else
                tymin = (vmax(2) - origin(2)) / direction(2);
                tymax = (vmin(2) - origin(2)) / direction(2);
            end

            if ( (tmin > tymax) || (tymin > tmax) )
                flag = 0;
                tmin = NaN;  %originally -1
                return;
            end

            if (tymin > tmin)
                tmin = tymin;
            end

            if (tymax < tmax)
                tmax = tymax;
            end

            if ((direction(3) == 0) && (((vmin(3) - origin(3))==0) || ((vmax(3) - origin(3))==0))) %comparing floats!
                tzmin=0;
                tzmax=inf;
            elseif (direction(3) >= 0)
                tzmin = (vmin(3) - origin(3)) / direction(3);
                tzmax = (vmax(3) - origin(3)) / direction(3);
            else
                tzmin = (vmax(3) - origin(3)) / direction(3);
                tzmax = (vmin(3) - origin(3)) / direction(3);
            end


            if ((tmin > tzmax) || (tzmin > tmax))
                flag = 0;
                tmin = -1;
                return;
            end

            if (tzmin > tmin)
                tmin = tzmin;
            end

            if (tzmax < tmax)
                tmax = tzmax;
            end

            % if( (tmin < t1) && (tmax > t0) )
            flag = 1;
            % else
            %    flag = 0;
            %    tmin = -1;
            % end;
        end

        function [ ndim, area, co] = detsub(i,facets,blocks,floors,boundingwalls,cornm,nblocks,xb,yb,zb,delta)
            % returns ndim=(dim1*dim2), area, center & corner of every facet (could be preprocessed and
            % saved)

            bi=facets(i,3); %block index
            fi=facets(i,1); %facet index
            ci=facets(i,4); %building index (-1 for roads, -101 for bounding wall)

            il=-1; iu=-1;jl=-1;ju=-1;kl=-1;ku=-1;it=-1;jt=-1;kt=-1;
            co=-1;
            ndim=-1;
            area=-1;




            if (ci<=-100) %it is a bounding wall
                kl=boundingwalls(bi,5)+1+1;
                ku=boundingwalls(bi,6)+1+1;
                if (fi==2)  %east, facing west
                    jl=boundingwalls(bi,3);     %lower y index of floor facet 1
                    ju=boundingwalls(bi,4);     %upper y index of floor facet 1
                    it=boundingwalls(bi,1);

                    x=xb(it+1);

                    ndim=((ju-jl)+1)*((ku-kl)+1);

                elseif (fi==3)  %west, facing east
                    jl=boundingwalls(bi,3);     %lower y index of floor facet 1
                    ju=boundingwalls(bi,4);     %upper y index of floor facet 1
                    it=boundingwalls(bi,2);

                    x=xb(it);

                    ndim=((ju-jl)+1)*((ku-kl)+1);
                elseif (fi==4)  %south, facing north
                    il=boundingwalls(bi,1);
                    iu=boundingwalls(bi,2);
                    jt=boundingwalls(bi,4);

                    y=yb(jt);
                    ndim=((iu-il)+1)*((ku-kl)+1);

                else %if (fi==5)  %north, facing south
                    jt=boundingwalls(bi,3);
                    il=boundingwalls(bi,1);
                    iu=boundingwalls(bi,2);

                    y=yb(jt+1);

                    ndim=((iu-il)+1)*((ku-kl)+1);
                end
            elseif (ci>=0) %deal with floors separately
                if (fi==1)  %top
                    il=blocks(bi,1);     %lower x index of facet 1
                    iu=blocks(bi,2);     %upper x index of facet 1
                    jl=blocks(bi,3);     %lower y index of facet 1
                    ju=blocks(bi,4);     %upper y index of facet 1
                    kt=blocks(bi,6)+1;
                    z=zb(kt+1);

                    ndim=((iu-il)+1)*((ju-jl)+1);

                elseif ((fi==2) || (fi==3))  %west / east
                    jl=blocks(bi,3);     %lower y index of facet 1
                    ju=blocks(bi,4);     %upper y index of facet 1
                    kl=blocks(bi,5)+1;     %lower z index of facet 1
                    ku=blocks(bi,6)+1;     %upper z index of facet 1

                    if fi==2  %west, facing west
                        it=blocks(bi,1);
                        x=xb(it);
                    else %east, facing east
                        it=blocks(bi,2);
                        x=xb(it+1);
                    end

                    ndim=((ju-jl)+1)*((ku-kl)+1);


                elseif ((fi==4) || (fi==5)) %north / south
                    il=blocks(bi,1) ;    %lower x index of facet 1
                    iu=blocks(bi,2) ;    %upper x index of facet 1
                    kl=blocks(bi,5)+1 ;    %lower y index of facet 1
                    ku=blocks(bi,6)+1 ;    %upper y index of facet 1

                    if fi==4 %north, facing north
                        jt=blocks(bi,4);
                        y=yb(jt+1);
                    else  %south, facing south
                        jt=blocks(bi,3);
                        y=yb(jt);
                    end

                    ndim=((iu-il)+1)*((ku-kl)+1);
                end
            end

            %return xyz-coordinates of center and 4 corners clockwise from bottom left
            %move coordinates of corners very slightly to the interior of the facet (by
            %delta)
            if (ci>=0 || ci<=-100) %it's not a floor -> it's a building or bounding wall

                if (fi==1)
                    co=[(xb(iu+1)+xb(il))*.5, (yb(ju+1)+yb(jl))*.5  , z; ... %center
                        xb(il)+delta        , yb(jl)+delta          , z; ... %4 corners
                        xb(il)+delta        , yb(ju+1)-delta        , z; ...
                        xb(iu+1)-delta      , yb(ju+1)-delta        , z; ...
                        xb(iu+1)-delta      , yb(jl)+delta          , z; ...
                        xb(il)        , yb(jl)          , z; ... %4 corners
                        xb(il)        , yb(ju+1)        , z; ...
                        xb(iu+1)      , yb(ju+1)        , z; ...
                        xb(iu+1)      , yb(jl)          , z; ...
                        ];
                    area=(xb(iu+1)-xb(il))*(yb(ju+1)-yb(jl));

                elseif ((fi==2) || (fi==3))

                    co=[ x  , (yb(ju+1)+yb(jl))*.5  , (zb(ku+1)+zb(kl))*.5; ... %center
                        x  , yb(jl)+delta          , zb(kl)+delta        ; ... %4 corners
                        x  , yb(jl)+delta          , zb(ku+1)-delta      ; ...
                        x  , yb(ju+1)-delta        , zb(ku+1)-delta      ; ...
                        x  , yb(ju+1)-delta        , zb(kl)+delta        ; ...
                        x  , yb(jl)          , zb(kl)        ; ... %4 corners
                        x  , yb(jl)          , zb(ku+1)      ; ...
                        x  , yb(ju+1)        , zb(ku+1)      ; ...
                        x  , yb(ju+1)        , zb(kl)        ; ...
                        ];
                    area=(zb(ku+1)-zb(kl))*(yb(ju+1)-yb(jl));

                elseif ((fi==4) || (fi==5))


                    co=[(xb(iu+1)+xb(il))*.5     , y  , (zb(ku+1)+zb(kl))*.5  ; ... %center
                        xb(il)+delta             , y  , zb(kl)+delta          ; ... %4 corners
                        xb(il)+delta             , y  , zb(ku+1)-delta        ; ...
                        xb(iu+1)-delta           , y  , zb(ku+1)-delta        ; ...
                        xb(iu+1)-delta           , y  , zb(kl)+delta          ; ...
                        xb(il)             , y  , zb(kl)          ; ... %4 corners
                        xb(il)             , y  , zb(ku+1)        ; ...
                        xb(iu+1)           , y  , zb(ku+1)        ; ...
                        xb(iu+1)           , y  , zb(kl)          ; ...
                        ];
                    area=(xb(iu+1)-xb(il))*(zb(ku+1)-zb(kl));
                end

            else  % it is a floor, not a building
                bif = bi - nblocks; % index in floor list
                il=floors(bif,1);     %lower x index of floor facet 1
                iu=floors(bif,2);     %upper x index of floor facet 1
                jl=floors(bif,3);     %lower y index of floor facet 1
                ju=floors(bif,4);     %upper y index of floor facet 1
                z=1;
                ndim=((iu-il)+1)*((ju-jl)+1);
                co=[(xb(iu+1)+xb(il))*.5, (yb(ju+1)+yb(jl))*.5  , z; ... %center
                    xb(il)+delta        , yb(jl)+delta          , z; ... %4 corners
                    xb(il)+delta        , yb(ju+1)-delta        , z; ...
                    xb(iu+1)-delta      , yb(ju+1)-delta        , z; ...
                    xb(iu+1)-delta      , yb(jl)+delta          , z; ...
                    xb(il)        , yb(jl)          , z; ... %4 corners
                    xb(il)        , yb(ju+1)        , z; ...
                    xb(iu+1)      , yb(ju+1)        , z; ...
                    xb(iu+1)      , yb(jl)          , z; ...
                    ];
                area=(xb(iu+1)-xb(il))*(yb(ju+1)-yb(jl));
                %if it is 1*1 nothing has to be done since it matches both walls
                if ju-jl>0  %floor facet along a west or east wall
                    %test if it is a wall-wall-floor corner, if so, the floor facet needs
                    %to be cut diagonally, see "createfloors.m"
                    if cornm(iu,ju)==10  %south/west corner
                        %JUST REMOVE IT ONE dx OR dy FROM THE OTHER BLOCK
                        co(2+1,2)=yb(ju)+delta;
                        co(6+1,2)=yb(ju);
                        area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
                    elseif cornm(iu,ju)==15  %south/east corner
                        co(3+1,2)=yb(ju)+delta;
                        co(7+1,2)=yb(ju);
                        area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
                    end
                    if cornm(iu,jl)==8  %north/west corner
                        co(1+1,2)=yb(jl+1)-delta;
                        co(5+1,2)=yb(jl+1);
                        area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
                    elseif cornm(iu,jl)==12  %north/east corner
                        co(4+1,2)=yb(jl+1)-delta;
                        co(8+1,2)=yb(jl+1);
                        area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
                    end

                elseif iu-il>0 %floor facet along along a south or north wall
                    if cornm(iu,ju)==10  %south/west corner
                        co(4+1,1)=xb(iu)+delta;
                        co(8+1,1)=xb(iu);
                        area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
                    elseif cornm(iu,ju)==8  %north/west corner
                        co(3+1,1)=xb(iu)+delta;
                        co(7+1,1)=xb(iu);
                        area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
                    end
                    if cornm(il,ju)==12  %north/east corner
                        co(2+1,1)=xb(il+1)-delta ;
                        co(6+1,1)=xb(il+1) ;
                        area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
                    elseif cornm(il,ju)==15  %south/east corner
                        co(1+1,1)=xb(il+1)-delta;
                        co(5+1,1)=xb(il+1);
                        area=area-(xb(2)-xb(1))*(yb(2)-yb(1))/2;
                    end
                end
            end
        end

    end
end
