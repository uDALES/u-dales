classdef da_pp < dynamicprops
    % Class for pre-processing in uDALES 
    %
    properties (Hidden = true, SetAccess = protected)
        path;                    % Path to simulations.
        cpath;                   % Current path.
        buf;                     % All purpose buffer.
        expnr;
    end
    
    properties (SetAccess = protected)
        g = 9.81;                % Gravitational acceleration.
        
    end
    methods (Static)
        
        function obj = da_pp(expnr, varargin)
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
        
        
        
        %make grid
        
        %do other things
        %  function f = loadfd(obj, field, vars, varargin)
        
        % Loads field files into struct obj.bd.(svar)
        % loadfd(obj, field, vars, varargin)
        %
        % field       : Name of field, e.g. 'F000100s'.
        % vars        : Cell arrary or requested variables.
        %               E.g. {'u', 'w'}.
        % varargin{1} : Custom name for property. (string).
        
        
        %  end
        
        % ------------------------------------------------------- %
        
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
        
        % ------------------------------------------------------- %
        
        function gopath(obj)
            % Go to simulation path.
            % gopath(obj)
            
            cd(obj.path);
        end
        
        % ------------------------------------------------------- %
        
        function gohome(obj)
            % Go to work path.
            % gohome(obj)
            
            cd(obj.cpath);
        end
        
        % ------------------------------------------------------- %
        
        function chcpath(obj, newpath)
            
            % Change work path.
            % chcpath(obj, newpath)
            
            here = pwd;
            cd(newpath)
            obj.cpath = pwd;
            cd(here);
        end
        
        function set_defaults(obj, ncpus)
            % Domain
            da_pp.addvar(obj, 'imax', 64)  % # cells in x-direction
            da_pp.addvar(obj, 'xsize', 64) % Domain size in x-direction
            da_pp.addvar(obj, 'jtot', 64)  % # cells in j-direction
            da_pp.addvar(obj, 'ysize', 64) % Domain size in y-direction
            da_pp.addvar(obj, 'kmax', 32)  % # cells in k-direction
            da_pp.addvar(obj, 'zsize', 32) % Domain size in z-direction
            da_pp.addvar(obj, 'imin', 1)
            da_pp.addvar(obj, 'jmin', 1)
            da_pp.addvar(obj, 'kmin', 1)
            da_pp.addvar(obj, 'dx', obj.xsize / obj.imax)
            da_pp.addvar(obj, 'dy', obj.ysize / obj.jtot)
            da_pp.addvar(obj, 'dz', obj.zsize / obj.kmax)
            da_pp.addvar(obj, 'lzstretch', 0) % switch for stretching z grid
            da_pp.addvar(obj, 'stretchconst', 0)
            
            % LIDAR parameters
            da_pp.addvar(obj, 'llidar', 0)
            da_pp.addvar(obj, 'sourcename', '')
            da_pp.addvar(obj, 'dxinp', 1) % resolution of image [m/pixel]
            da_pp.addvar(obj, 'dyinp', 1)
            da_pp.addvar(obj, 'dzinp', 1)
            da_pp.addvar(obj, 'centeri', 0) % center of area of interest in original image [pixel]
            da_pp.addvar(obj, 'centerj', 0)
            da_pp.addvar(obj, 'maxh', 0) % magimum height of buildings in image [m]
            da_pp.addvar(obj, 'pad', 5) % padding. A padding of 0 makes only sense for idealised cases. There should be no building at domain edge
            da_pp.addvar(obj, 'smallarea', round(150 / (obj.dx * obj.dy))) % objects smaller than this will be deleted
            
            if ceil(obj.jtot / ncpus) ~= floor (obj.jtot / ncpus)
                disp(['Possible jtot: ' num2str([2 3 4 5 6 7 8] * ncpus)])
                error('No. CPUs does not fit j grid size')
            else
                %disp('cpus and jtot successful')
            end
            
            % Blocks
            da_pp.addvar(obj, 'lblocksfile', 0) % switch for using blocks from a file
            da_pp.addvar(obj, 'blocksfile', '') % name of blocks file
            
            da_pp.addvar(obj, 'lcastro', 0) % switch for staggered cubes
            da_pp.addvar(obj, 'lcube', 0)   % switch for linear cubes
            da_pp.addvar(obj, 'lblocks', 0) % switch for infinite blocks
            da_pp.addvar(obj, 'lflat', 0) % switch for flat domain

            %da_pp.addvar(obj, 'bl', [])
            da_pp.addvar(obj, 'blocks', [])

            da_pp.addvar(obj, 'blockheight', 16) % block height
            da_pp.addvar(obj, 'blockwidth', 16)  % block width
            da_pp.addvar(obj, 'canyonwidth', 16) % canyonwidth
            
            % Profiles
            %da_pp.addvar(obj, 'lmassflowr', 0) % switch for constant mass
            %flow rate.
            da_pp.addvar(obj, 'luoutflowr', 0) % switch for constant volume flow rate in u.
            da_pp.addvar(obj, 'lvoutflowr', 0) % switch for constant volume flow rate in u.
            da_pp.addvar(obj, 'luvolflowr', 0) % switch for constant volume flow rate in u.
            da_pp.addvar(obj, 'lvvolflowr', 0) % switch for constant volume flow rate in v. 
            da_pp.addvar(obj, 'lprofforc', 0)  % switch for 1D geostrophic forcing
            da_pp.addvar(obj, 'lcoriol', 0)    % switch for coriolis forcing
            if (not(obj.luoutflowr) && not(obj.lvoutflowr) && not(obj.luvolflowr) && not(obj.lvvolflowr) && not(obj.lprofforc) && not(obj.lcoriol))
                da_pp.addvar(obj, 'ldp', 1)
                %disp('No forcing switch config. setup so initial velocities and pressure gradients applied.')
            else
                da_pp.addvar(obj, 'ldp', 0)
            end
            
            % For forcing, need to check the value given in physics namelist 

            da_pp.addvar(obj, 'u0', 0) % initial u-velocity - also applied as geostrophic term where applicable
            da_pp.addvar(obj, 'v0', 0) % initial v-velocity - also applied as geostrophic term where applicable
            
            da_pp.addvar(obj, 'dpdx', 0) % dp/dx [Pa/m]
            da_pp.addvar(obj, 'dpdy', 0) % dp/dy [Pa/m]

            % 
            da_pp.addvar(obj, 'thl0', 288) % temperature at lowest level
            da_pp.addvar(obj, 'qt0', 0)    % specific humidity
            
            
            da_pp.addvar(obj, 'sv10', 0)   % scalar
            da_pp.addvar(obj, 'sv20', 0)   % scalar
            da_pp.addvar(obj, 'lapse', 0)  % lapse rate [K/s]
            
            % Other
            da_pp.addvar(obj, 'w_s',0) % subsidence [*units?*]
            da_pp.addvar(obj, 'R',0)   % radiative forcing [*units?*]

            % Walls
            da_pp.addvar(obj, 'z0horiz', 0.01)
            da_pp.addvar(obj, 'z0hhoriz', 0.000067)
            da_pp.addvar(obj, 'Thoriz', 288)
            da_pp.addvar(obj, 'Twest', 288)
            da_pp.addvar(obj, 'Teast', 288)
            da_pp.addvar(obj, 'Tnorth', 288)
            da_pp.addvar(obj, 'Tsouth', 288)
                  
            % Pollutants/chemistry
            da_pp.addvar(obj, 'nsv', 0)
            da_pp.addvar(obj, 'lchem' , 0) % switch for chemistry
            da_pp.addvar(obj, 'NOb' , 0)
            da_pp.addvar(obj, 'NO2b', 0)
            da_pp.addvar(obj, 'O3b', 0)
            
            % Trees
            da_pp.addvar(obj, 'ltrees', 0)
            if obj.ltrees
                da_pp.addvar(obj, 'tree_dz',0)
                da_pp.addvar(obj, 'tree_dx',0)
                da_pp.addvar(obj, 'tree_h',0)
                da_pp.addvar(obj, 'tree_w',0)
                da_pp.addvar(obj, 'tree_b',0)

                da_pp.addvar(obj, 'nt1',0)
                da_pp.addvar(obj, 'md',0)
                da_pp.addvar(obj, 'ww',0)
                da_pp.addvar(obj, 'lw',0)
                da_pp.addvar(obj, 'nt2',0)

                % Different cases?
            end

            
            % Purifiers
            da_pp.addvar(obj, 'lpurif', 0)
            if obj.lpurif
                if obj.lblocks
                    da_pp.addvar(obj, 'purif_dz', 1)  % purifier starting point from bottom
                    da_pp.addvar(obj, 'purif_dx', 3)  % distance from block
                    da_pp.addvar(obj, 'purif_h', 3)   % purifier height
                    da_pp.addvar(obj, 'purif_w', 0)   % purifier width
                    da_pp.addvar(obj, 'purif_dy', 1)  % depth of purifier (in y)
                    da_pp.addvar(obj, 'purif_sp', 31) % spacing between purifiers
                    da_pp.addvar(obj, 'purif_i', 1)   % case for purifier (1 = +ve x, 2 = -ve x, 3 = +ve y etc.)
                    da_pp.addvar(obj, 'npurif', obj.jtot / (obj.npurif_dy + obj.purif_sp))

                    if ceil(npurif) ~= floor(npurif)
                        lp = 0:obj.tot / 2;
                        indp = rem(obj.jtot / 2, lp) == 0;
                        errp = ([lp(indp); (obj.jtot / 2) ./ lp(indp)]);
                        disp('Purifier layout does not fit grid')
                        disp(['sum widths to: ' num2str(errp(1,:))])
                        disp(['Current width: ' num2str(obj.purif_dy + obj.purif_sp)])
                        error('Incorrect purifier layout')
                    else
                        %disp('Successful purifer layout')
                        %disp(['Number of purifiers: ' num2str(npurif*2+1)])
                    end
                else
                    error('Must use lblocks configuration to use purifiers')
                end
            end
            
            % EB
            da_pp.addvar(obj, 'solaz', 135); % azimuth angle
            da_pp.addvar(obj, 'Z', 28.4066); % zenith angle
            %da_pp.addvar(obj, 'delta', 0.01); % small adjustment for position of rad corners
            da_pp.addvar(obj, 'centerweight', 12 / 32);
            da_pp.addvar(obj, 'cornerweight', (1 - obj.centerweight) / 4);
            da_pp.addvar(obj, 'I', 184.8775); % Direct solar irradiation [W/m2]
            da_pp.addvar(obj, 'Dsk', 418.8041); % Diffuse incoming radiation [W/m2]
            
            if obj.lEB
                da_pp.addvar(obj, 'maxsize', 10); % Add to namoptions
            else
                da_pp.addvar(obj, 'maxsize', inf);
            end
            
            %da_pp.addvar(obj, 'blocks', []);
            da_pp.addvar(obj, 'facets', []);           
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
             da_pp.addvar(obj, 'xf', 0.5 * obj.dx : obj.dx : obj.xsize - 0.5 * obj.dx); 
             da_pp.addvar(obj, 'yf', 0.5 * obj.dy : obj.dy : obj.ysize - 0.5 * obj.dy);
             da_pp.addvar(obj, 'xh', 0 : obj.dx : obj.xsize);
             da_pp.addvar(obj, 'yh', 0 : obj.dy : obj.ysize);
        end
        
        function write_xgrid(obj)
            xgrid = fopen(['xgrid.inp.' obj.expnr], 'w');
            fprintf(xgrid, '%12s\n', '#     x-grid');
            fprintf(xgrid, '%12s\n', '#           ');
            fprintf(xgrid, '%-20.15f\n', obj.xf);
            fclose(xgrid);
            %disp(['... written xgrid.inp.' obj.expnr]) 
        end
        
        function generate_zgrid(obj)
            if ~obj.lzstretch
                da_pp.addvar(obj, 'zf', 0.5 * obj.dz : obj.dz : obj.zsize - 0.5 * obj.dz);
                da_pp.addvar(obj, 'zh', 0 : obj.dz : obj.zsize);
                da_pp.addvar(obj, 'dzf', obj.zh(2:end) - obj.zh(1:end - 1));
            else
                if obj.lstretchexp
                   da_pp.stretch_exp(obj)
                elseif obj.lstretchtanh
                    da_pp.stretch_tanh(obj)
                elseif obj.lstretch2tanh
                    da_pp.stretch_2tanh(obj)                   
                else
                    error('Invalid stretch');
                end
            end
        end
        
        function stretch_exp(obj)
            il = round(obj.maxh / obj.dz); % so: what should maxh be replaced with?
            ir  = obj.kmax - il;
            
            da_pp.addvar(obj, 'zf', zeros(obj.kmax, 1));
            da_pp.addvar(obj, 'dzf', zeros(obj.kmax, 1));
            da_pp.addvar(obj, 'zh', zeros(obj.kmax+1, 1));
            
            obj.zf(1:il) = obj.dz / 2 : obj.dz : obj.maxh; % so: what should maxh be replaced with?
            obj.zh(1:il+1) = 0 : obj.dz : obj.maxh; % so: what should maxh be replaced with?
            
            gf = obj.stretchconst;
            
            while true
                obj.zh(il + 1:end) = obj.zh(il + 1) + (obj.zsize - obj.zh(il+1)) * (exp(gf * (0:1:ir) / (ir)) - 1)/(exp(gf) - 1); %dh has been replaced by zsize                
                if (obj.zh(il+2) - obj.zh(il + 1)) < obj.dz
                    gf = gf - 0.01; %make sufficiently small steps to avoid an initial bump in dz
                    disp(['Decreasing stretchconst to:' num2str(gf)])
                    
                else
                    if (obj.zh(end) - obj.zh(end - 1)) > 3 * obj.dz
                        disp('WARNING: final grid spacing large! Consider reducing domain height')
                    end
                    break
                end
                
            end
            
            for i = 1:obj.kmax
                obj.zf(i) = (obj.zh(i) + obj.zh(i+1)) / 2 ;
                obj.dzf(i) = obj.zh(i+1) - obj.zh(i);
            end
            
            %disp(['growth factor ~' num2str((obj.dzf(il+1)/obj.dzf(il)-1)*100) '-' num2str((obj.dzf(end)/obj.dzf(end-1)-1)*100) '%'])
        end
        
        function stretch_tanh(obj)
            il = round(obj.maxh / obj.dz); % so: what should maxh be replaced with?
            ir  = obj.kmax - il;
            
            da_pp.addvar(obj, 'zf', zeros(obj.kmax, 1));
            da_pp.addvar(obj, 'dzf', zeros(obj.kmax, 1));
            da_pp.addvar(obj, 'zh', zeros(obj.kmax+1, 1));
            
            obj.zf(1:il) = obj.dz / 2 : obj.dz : obj.maxh; % so: what should maxh be replaced with?
            obj.zh(1:il+1) = 0 : obj.dz : obj.maxh; % so: what should maxh be replaced with?
            
            gf = obj.stretchconst;

            while true
                obj.zh(il + 1:end) = obj.zh(il + 1) + (obj.zsize - obj.zh(il + 1)) * (1 - tanh(gf * (1 - 2 * (0:1:ir)' / (2*ir))) / tanh(gf));

            if (obj.zh(il + 2) - obj.zh(il + 1)) < obj.dz
                gf = gf - 0.01; %make sufficiently small steps to avoid an initial bump in dz
                %disp(['Decreasing stretchconst to:' num2str(gf)])

            else
                if (obj.zh(end) - obj.zh(end - 1)) > 3 * obj.dz
                disp('WARNING: final grid spacing large! Consider reducing domain height') 
                end
                break
            end

            end

            for i = 1:obj.kmax
                obj.zf(i) = (obj.zh(i) + obj.zh(i+1)) / 2 ;
                obj.dzf(i) = obj.zh(i+1) - obj.zh(i);
            end

            %disp(['growth factor ~' num2str((obj.dzf(il+1) / obj.dzf(il) - 1) * 100) '-' num2str((obj.dzf(end) / dzf(end-1)-1)*100) '%'])
        end
        
        function stretch_2tanh(obj)
            il = round(maxh / obj.dz); % so: what should maxh be replaced with here?
            ir  = obj.kmax - il;
            
            da_pp.addvar(obj, 'zf', zeros(obj.kmax, 1));
            da_pp.addvar(obj, 'dzf', zeros(obj.kmax, 1));
            da_pp.addvar(obj, 'zh', zeros(obj.kmax+1, 1));
            
            obj.zf(1:il) = obj.dz / 2 : obj.dz : obj.maxh; % so: what should maxh be replaced with?
            obj.zh(1:il+1) = 0 : obj.dz : obj.maxh; % so: what should maxh be replaced with?
            
            gf = obj.stretchconst;
            
            while true
                obj.zh(il+1:end) = obj.zh(il+1) + (obj.zsize - obj.zh(il+1)) / 2 * (1 - tanh(gf * (1 - 2 * (0:1:ir)'/(ir))) / tanh(gf));
                
                if (obj.zh(il + 2) - obj.zh(il + 1)) < obj.dz
                    gf = gf - 0.01; %make sufficiently small steps to avoid an initial bump in dz
                    %disp(['Decreasing stretchconst to:' num2str(gf)])
                    
                else
                    if (max(diff(obj.zh))) > 3 * dz
                        disp('WARNING: final grid spacing large! Consider reducing domain height')
                    end
                    break
                end               
            end
            
            for i = 1:obj.zmax
                obj.zf(i) = (obj.zh(i) + obj.zh(i+1)) / 2 ;
                obj.dzf(i) = obj.zh(i + 1) - obj.zh(i);
            end            
            %disp(['growth factor ~' num2str((obj.dzf(il+1)/obj.dzf(il)-1)*100) '-' num2str((obj.dzf(end)/obj.dzf(end-1)-1)*100) '%'])            
        end
        
        function write_zgrid(obj)
            zgrid = fopen(['zgrid.inp.' obj.expnr], 'w');
            fprintf(zgrid, '%12s\n', '#     z-grid');
            fprintf(zgrid, '%12s\n', '#           ');
            fprintf(zgrid, '%-20.15f\n', obj.zf);
            fclose(zgrid);
            %disp(['... written zgrid.inp.' obj.expnr]) 
        end
        
        function generate_lscale(obj)
            if ((obj.luoutflowr || obj.lvoutflowr) + (obj.luvolflowr || obj.lvvolflowr) + obj.lprofforc + obj.lcoriol + obj.ldp) > 1
                error('More than one forcing type specified')
            end
            
            da_pp.addvar(obj, 'ls', zeros(length(obj.zf), 10));
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
            %disp(['... written lscale.inp.' obj.expnr]) 
        end
        
        function generate_prof(obj)
            da_pp.addvar(obj, 'pr', zeros(length(obj.zf), 6));
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
            %pr(:,6) = tke;
        end
        
        function write_prof(obj)
            prof = fopen(['prof.inp.' obj.expnr], 'w');
            fprintf(prof, '%-12s\n', '# SDBL flow');
            fprintf(prof, '%-60s\n', '# z thl qt u v tke');
            fprintf(prof, '%-20.15f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f\n', obj.pr');
            fclose(prof);
            %disp(['... written prof.inp.' obj.expnr])
        end
        
        function generate_scalar(obj)
            if obj.lchem
                % Do something
            else
                % Do something else
            end
            % For the moment:
            
            if obj.nsv > 0 
                da_pp.addvar(obj, 'sc', zeros(length(obj.zf), 5));
                obj.sc(:,1) = obj.zf;
                obj.sc(:,2) = obj.sv10;
                obj.sc(:,3) = obj.sv20;
            end
        end
        
        function write_scalar(obj)
            if obj.nsv > 0
                scalar = fopen(['scalar.inp.' obj.expnr], 'w');
                fprintf(scalar, '%-12s\n', '# SDBL flow');
                fprintf(scalar, '%-60s\n', '# z sca1 sca2 sca3 sca4');
                fprintf(scalar, '%-20.15f %-14.10f %-14.10f %-14.10f %-14.10f\n', obj.sc');
                fclose(scalar);
                %disp(['... written scalar.inp.' obj.expnr])
            end
        end
        
        function generate_topo_from_bl(obj)
            da_pp.addvar(obj, 'topomask', zeros(obj.jtot, obj.imax));
            da_pp.addvar(obj, 'topo', zeros(obj.jtot, obj.imax));
            
            % if no blocks add lowest level
            if isnan(obj.bl)
                
            else
                for n = 1:size(obj.bl, 1)
                    obj.topo(obj.bl(n,3):obj.bl(n,4),obj.bl(n,1):obj.bl(n,2)) = obj.zh(obj.bl(n,6) + 1);
                    obj.topomask(obj.bl(n,3):obj.bl(n,4),obj.bl(n,1):obj.bl(n,2)) = 1;
                end
            end
        end
        
        function generate_bl_from_namoptions(obj)
            %aspectratio = r.zh(r.blockheight+1)/(r.xh(r.canyonwidth+1)-r.xh(1));
            %aspectratio = obj.zh(obj.blockheight + 1) / (obj.xh(obj.canyonwidth + 1) - obj.xh(1));
            %nrows = ie/(r.blockwidth+r.canyonwidth);
            
            imax = obj.imax;
            jtot = obj.jtot;
            blockwidth = obj.blockwidth;
            blockheight = obj.blockheight;
            canyonwidth = obj.canyonwidth;
            nrows =  imax / (blockwidth + canyonwidth);
            %da_pp.addvar(obj, 'nrows', obj.imax / (obj.blockwidth + obj.canyonwidth));
            
            %if obj.lflat

            if ceil(nrows) ~= floor(nrows)
                %l = 0:ie/2;
                l = 0:0.5 * imax;
                %ind = rem(ie/2,l)==0; %// logical index that tells if remainder is zero or not
                ind = rem(0.5 * imax, l) == 0;
                %err = ([l(ind); (ie/2)./l(ind)]);
                err = ([l(ind); (0.5 * imax) ./ l(ind)]);
                disp('Block system does not fit grid')
                disp(['sum widths to: ' num2str(err(1,:))])
                %disp(['Current width: ' num2str(r.blockwidth + r.canyonwidth)])
                disp(['Current width: ' num2str(blockwidth + canyonwidth)])
                error('Incorrect block system')
            else
                %disp('Successful block network')
                %disp(['aspect ratio: ' num2str(aspectratio)])
            end 
            
            jmin = obj.jmin;
            
            %if obj.lflat
                %bl = [];
            if obj.lcastro
                %nrows = ie/(r.blockwidth*2);
                nrows = imax / (blockwidth * 2);
                %ncolumns = je/(r.blockwidth*2);
                ncolumns = jtot / (blockwidth * 2);
                %da_pp.addvar(obj, 'ncolumns', obj.jtot / (obj.blockwidth * 2));
                bl = zeros(nrows * ncolumns + nrows / 2, 13);
                %obj.bl = zeros(obj.nrows * obj.ncolumns + obj.nrows / 2, 13);
                bl(:, 5) = 0;
                bl(:, 6) = blockheight - 1;
                

                for n = 1:nrows
                    for nn = 0:ncolumns
                        if mod(n, 2) == 0
                            if nn == 0
                                %bl(length(nonzeros(bl(:,3)))+1,3) = jb;
                                bl(length(nonzeros(bl(:,3))) + 1, 3) = jmin;
                                %bl(length(nonzeros(bl(:,4))) + 1, 4) = r.blockwidth/2;
                                bl(length(nonzeros(bl(:,4))) + 1, 4) = 0.5 * blockwidth;
                            elseif nn == ncolumns
                                %bl(length(nonzeros(bl(:,3)))+1,3) = je-r.blockwidth/2+1;
                                bl(length(nonzeros(bl(:,3))) + 1, 3) = jtot - 0.5 * blockwidth + 1;
                                %bl(length(nonzeros(bl(:,4)))+1,4) = je;
                                bl(length(nonzeros(bl(:,4))) + 1, 4) = jtot;
                            else
                                %bl(length(nonzeros(bl(:,3)))+1,3) = jb + nn * r.blockwidth*2 - r.blockwidth/2;
                                bl(length(nonzeros(bl(:,3))) + 1, 3) = jmin + nn * blockwidth * 2 - 0.5 * blockwidth;
                                %bl(length(nonzeros(bl(:,4)))+1,4) = bl(length(nonzeros(bl(:,4)))+1,3) + r.blockwidth - 1;
                                bl(length(nonzeros(bl(:,4))) + 1, 4) = bl(length(nonzeros(bl(:,4))) + 1, 3) + blockwidth - 1;
                            end
                            %bl(length(nonzeros(bl(:,1)))+1,1) = -0.5*r.blockwidth + (2*n-1) * r.blockwidth + 1;
                            bl(length(nonzeros(bl(:,1))) + 1, 1) = - 0.5 * blockwidth + (2 * n - 1) * blockwidth + 1;
                            %bl(length(nonzeros(bl(:,2)))+1,2) = bl(length(nonzeros(bl(:,2)))+1,1) + r.blockwidth - 1;
                            bl(length(nonzeros(bl(:,2))) + 1, 2) = bl(length(nonzeros(bl(:,2)))+1,1) + blockwidth - 1;
                        end
                    end
                    for nn = 0:ncolumns - 1
                        if mod(n,2) ~= 0
                            %bl(length(nonzeros(bl(:,3)))+1,3) = jb + r.blockwidth + nn * r.blockwidth*2 - r.blockwidth/2;
                            bl(length(nonzeros(bl(:,3))) + 1, 3) = jmin + blockwidth + nn * blockwidth * 2 - 0.5 * blockwidth;
                            %bl(length(nonzeros(bl(:,4)))+1,4) = bl(length(nonzeros(bl(:,4)))+1,3) + r.blockwidth - 1;
                            bl(length(nonzeros(bl(:,4))) + 1, 4) = bl(length(nonzeros(bl(:,4))) + 1, 3) + blockwidth - 1;
                            %bl(length(nonzeros(bl(:,1)))+1,1) = -0.5*r.blockwidth + (2*n-1) * r.blockwidth + 1;
                            bl(length(nonzeros(bl(:,1))) + 1, 1) = - 0.5 * blockwidth + (2 * n - 1) * blockwidth + 1;
                            %bl(length(nonzeros(bl(:,2)))+1,2) = bl(length(nonzeros(bl(:,2)))+1,1) + r.blockwidth - 1;
                            bl(length(nonzeros(bl(:,2))) + 1, 2) = bl(length(nonzeros(bl(:,2))) +1 , 1) + blockwidth - 1;
                        end

                    end
                end

            elseif obj.lcube
                xsize = obj.xsize; zsize = obj.zsize;
                kmax = obj.kmax;
                %nrows = ie/(r.blockwidth*2);
                nrows = imax / (blockwidth * 2);
                %ncolumns = je/(r.blockwidth*2);
                ncolumns = jtot / (blockwidth * 2);
                %da_pp.addvar(obj, 'ncolumns', obj.jtot / (obj.blockwidth * 2));
                bl = zeros(nrows * ncolumns, 13);
                %obj.bl = zeros(obj.nrows * obj.ncolumns, 13);
                for n = 1 : nrows
                    for nn = 0 : ncolumns - 1
                        %bl((n - 1) * ncolumns + nn + 1, 1) = -0.5 * r.blockwidth + (2 * n - 1) * r.blockwidth + 1;
                        bl((n - 1) * ncolumns + nn + 1, 1) = -0.5 * blockwidth + (2 * n - 1) * blockwidth + 1;
                        %bl((n-1)*ncolumns+nn+1,2) = bl((n-1)*ncolumns+nn+1,1) + r.blockwidth - 1;
                        bl((n - 1) * ncolumns + nn + 1, 2) = bl((n - 1) * ncolumns + nn + 1, 1) + blockwidth - 1;
                        %bl((n-1)*ncolumns+nn+1,5) = 0;
                        bl((n - 1) * ncolumns + nn + 1, 5) = 0;
                        %bl((n-1)*ncolumns+nn+1,6) = ceil(r.blockwidth*(h/ie)/(hz/ke));
                        bl((n - 1) * ncolumns + nn + 1, 6) = ceil(blockwidth * (xsize / imax) / (zsize / kmax));
                        %bl((n-1)*ncolumns+nn+1,3) = jb + r.blockwidth/2 + nn * r.blockwidth*2;
                        bl((n - 1) * ncolumns + nn + 1, 3) = jmin + blockwidth / 2 + nn * blockwidth * 2;
                        %bl((n-1)*ncolumns+nn+1,4) = bl((n-1)*ncolumns+nn+1,3) + r.blockwidth - 1;
                        bl((n - 1) * ncolumns + nn + 1, 4) = bl((n - 1) * ncolumns + nn + 1, 3) + blockwidth - 1;
                    end
                end

            elseif obj.lblocks % should be lcanyons;
                %bl = zeros(nrows, 13);
                bl = zeros(nrows, 13);
                %bl(1:nrows,1) = (r.canyonwidth/2+1:r.canyonwidth+r.blockwidth:ie-r.canyonwidth/2)';
                bl(1:nrows, 1) = (0.5 * canyonwidth + 1 : canyonwidth + blockwidth : imax - 0.5 * canyonwidth)';
                %bl(1:nrows,2) = bl(1:nrows,1) + r.blockwidth - 1;
                bl(1:nrows, 2) = bl(1:nrows, 1) + blockwidth - 1;
                %bl(:,3) = jb;
                bl(:, 3) = jmin;
                %bl(:,4) = je;
                bl(:, 4) = jtot;
                bl(:, 5) = 0;
                %bl(1:nrows,6) = r.blockheight-1;
                bl(1:nrows, 6) = blockheight - 1; 
            end 
            
            da_pp.addvar(obj, 'bl', bl)

            
            %maximum size of floors and bounding walls (cells in each dimension)
            
            %obj.blocks = obj.bl;
            
            % if no blocks add lowest level
%             if isnan(obj.blocks)
%             else
%                 for n = 1:size(obj.blocks,1)
%                     topo(obj.blocks(n,3):obj.blocks(n,4),obj.blocks(n,1):obj.blocks(n,2)) = zh(obj.blocks(n,6)+1);
%                     topomask(obj.blocks(n,3):obj.blocks(n,4),obj.blocks(n,1):obj.blocks(n,2)) = 1;
%                 end
%             end
            
            %bl2blocks_temp
            %makeblocks
            
            %da_pp.generate_topomask_from_blocks(obj)
            %da_pp.makeblocks(obj)  
        end
        
        function generate_bl_from_file(obj)
            obj.bl = dlmread(obj.blocksfile, '', 2, 0);
        end
        
%         function generate_topo_from_LIDAR(obj, sourcename)
%             %disp(obj.sourcename)
%             A = imread(sourcename);  %read topo image
%             [njorig, niorig, ~] = size(A);
%             
% 
%             % since its a greyscale image all 3 rgb channels have the same value
%             % substract value from 255 and scale by max value to get topography
%             topot = (255 - double(A(:,:,1))) / 255 * obj.maxh;
%             
%             %if ltestplot
%                 %figure; imagesc(topot)
%             %end
%             
%             % interpolate to a coarser grid if necessary
%             if obj.dx ~= obj.dxinp || obj.dy ~= obj.dyinp %can't use imresize, have to do it manually since x and y scaling might be different
%                 xxxo = 0.5 * obj.dxinp:obj.dxinp:(niorig * obj.dxinp); % original pixel centres
%                 yyyo = 0.5 * obj.dyinp:obj.dyinp:(njorig * obj.dyinp);
%                 xxx = 0.5 * obj.dx:obj.dx:(niorig * obj.dxinp); %desired pixel centres
%                 yyy = 0.5 * obj.dy:obj.dy:(njorig * obj.dxinp);
%                 [X,Y] = meshgrid(xxx, yyy);
%                 topot = interp2(yyyo, xxxo, topot', Y, X, 'nearest');
%                 clear xxxo yyyo X Y
%             end
%             
%             %if ltestplot
%                 %figure; imagesc(topot);
%             %end
%             
%             
%             % upper and lower coordinates to select
%             ip = round(obj.centeri / obj.dx) + obj.imax / 2 - 1 - obj.pad;
%             im = round(obj.centeri / obj.dx) - obj.imax / 2 + obj.pad;
%             jp = round(obj.centerj / obj.dy) + obj.jtot / 2 - 1 - obj.pad;
%             jm = round(obj.centerj / obj.dy) - obj.jtot / 2 + obj.pad;
%             
%             da_pp.addvar(obj, 'topo', zeros(obj.jtot, obj.imax));
%             obj.topo(obj.pad + (1:(obj.jtot - 2 * obj.pad)), obj.pad + (1:(obj.imax - 2 * obj.pad))) = topot(jm:jp, im:ip);
%             
% %             if ltestplot
% %                 figure; imagesc(xf,yf,topo)
% %             end
%             
%             % round height to grid
%             obj.topo = round(obj.topo / obj.dz) * obj.dz;
%             topoinit = obj.topo;
%             
%             
%             %% manual mainpulation of the topograpy come here, if necessary
%             %test for varying building height
%             % topo(170:180,173:209)=topo(170:180,173:209)+10;
%             % topo(180:184,291:331)=topo(180:184,291:331)+12.5;
%             % topo(143:148,173:214)=topo(143:148,173:214)+7.5;
%             
%             
%             %% processing
%             %% fill big holes
%             topoh = imfill(obj.topo, 'holes');
%             to = topoh - obj.topo;
%             toi = imbinarize(to);
%             
%             toi = bwareaopen(toi, obj.smallarea);
%             topo = topoh - toi.*to;
%             clear to toi topoh topot
%             
% %             if ltestplot
% %                 figure; imagesc(xf,yf,topo); title('after filling holes')
% %             end
%             
%             % create building mask
%             da_pp.addvar(obj, 'topomask', imbinarize(topo));
%             
%             %% remove small objects
%             obj.topomask = bwareaopen(obj.topomask, obj.smallarea);
%             
%             obj.topo = obj.topo .* obj.topomask;
%             
% %             if ltestplot
% %                 figure
% %                 imagesc(xf,yf,topo)
% %                 title('after removing small objects')
% %                 xlim([xh(1) xh(end)])
% %                 ylim([yh(1) yh(end)])
% %             end
%             
%             %% fill 1 cell gaps  %potentially problematic if the gap is between blocks of different height
%             
%             %function [mask, image ] = fillgaps( mask,image,nj,ni,pad)
%                 %fills horizontal and vertical 1D gaps of width 1
%                 while true
%                     change = false;
%                     for j = obj.pad + 1:obj.jtot - obj.pad - 1
%                         for i = obj.pad + 1:obj.imax - obj.pad - 1
%                             if j == 1 || i == 1 %do nothing at domain edge
%                                 continue
%                                 
%                             elseif sum(sum(obj.topomask(j-1:j,i-1:i+1))) == 5 && obj.topomask(j,i) == 0 && obj.topomask(j+1,i) == 0 %towards south. If it has 5 block neighbours, is not a block and j+1 is also not a block then fill the gap towards the south. remember that images are upside down, so j+1 is going down
%                                 jj = j;
%                                 while true
%                                     obj.topomask(jj, i)=1;
%                                     obj.topo(jj,i) = obj.topo(jj-1,i);
%                                     jj=jj+1;
%                                     if jj==obj.jtot %reached end of domain
%                                         break
%                                     end
%                                     if ~(sum(sum(obj.topomask(jj-1:jj,i-1:i+1)))==5 && obj.topomask(jj,i)==0 && obj.topomask(jj+1,i)==0)
%                                         break
%                                     end
%                                 end
%                                 change = true;
%                             elseif sum(sum(obj.topomask(j:j+1,i-1:i+1)))==5 && obj.topomask(j,i) == 0 && obj.topomask(j-1,i)==0 %towards north
%                                 jj=j;
%                                 while true
%                                     obj.topomask(jj,i)=1;
%                                     obj.topo(jj,i)=obj.topo(jj+1,i);
%                                     jj=jj-1;
%                                     if jj==1 %reached end of domain
%                                         break
%                                     end
%                                     if ~(sum(sum(obj.topomask(jj:jj+1,i-1:i+1)))==5 && obj.topomask(jj,i)==0 && obj.topomask(jj-1,i)==0)
%                                         break
%                                     end
%                                 end
%                                 change=true;
%                             elseif sum(sum(obj.topomask(j-1:j+1,i-1:i)))==5 && obj.topomask(j,i)==0 && obj.topomask(j,i+1)==0 %towards east
%                                 ii=i;
%                                 while true
%                                     obj.topomask(j,ii)=1;
%                                     obj.topo(j,ii)=obj.topo(j,ii-1);
%                                     ii=ii+1;
%                                     if ii==obj.imax %reached end of domain
%                                         break
%                                     end
%                                     if ~(sum(sum(obj.topomask(j-1:j+1,ii-1:ii)))==5 && obj.topomask(j,ii)==0 && obj.topomask(j,ii+1)==0)
%                                         break
%                                     end
%                                 end
%                                 change=true;
%                             elseif sum(sum(obj.topomask(j-1:j+1,i:i+1)))==5 && obj.topomask(j,i)==0 && obj.topomask(j,i-1)==0 %towards west
%                                 ii=i;
%                                 while true
%                                     obj.topomask(j,ii)=1;
%                                     obj.topo(j,ii)=obj.topo(j,ii+1);
%                                     ii=ii-1;
%                                     if ii==1 %reached end of domain
%                                         break
%                                     end
%                                     if ~(sum(sum(obj.topomask(j-1:j+1,ii:ii+1)))==5 && obj.topomask(j,ii)==0 && obj.topomask(j,ii-1)==0)
%                                         break
%                                     end
%                                 end
%                             end
%                             
%                         end
%                         
%                         
%                     end
%                     if ~change
%                         break
%                     end
%                 end
%             %end
%             
%             
%             %[topomask, topo] = fillgaps(topomask,topo,nj,ni,pad);
% %             if ltestplot
% %                 figure
% %                 imagesc(xf,yf,topo)
% %                 title('after filling gaps (size 1)')
% %                 xlim([xh(1) xh(end)])
% %                 ylim([yh(1) yh(end)])
% %             end
%             
%             % NEED TO IMPLEMENT SMOOTHBORDERS
%             %if pad == 0  %buildings to the edge, make sure there is no weird gaps at domain edge
%                 %[topomask, topo] = smoothborders(obj.topomask, obj.topo,nj,ni,pad);
%             %end
% %             if ltestplot
% %                 figure
% %                 imagesc(xf,yf,topo)
% %                 title('after smoothing borders')
% %                 xlim([xh(1) xh(end)])
% %                 ylim([yh(1) yh(end)])
% %             end
%             
%             %% remove cells with only 1 neighbour
%             data = obj.topo;
%             datamask = obj.topomask;
%             
%             c2 = 1;
%             while c2 ~= 0  %repeat until there is no cell with only one neighbour left
%                 c2 = 0;
%                 for i = 2:obj.imax - 1
%                     for j = 2:obj.jtot - 1
%                         count = 0;
%                         if datamask(j,i) > 0
%                             c = datamask(j,i-1) + datamask(j,i+1) + datamask(j+1,i) + datamask(j-1,i);
%                             if c > 1
%                                 continue
%                             else
%                                 c2 = c2 + 1;
%                                 datamask(j,i) = 0;
%                                 data(j,i) = 0;
%                             end
%                         end
%                     end
%                 end
%             end
%             
%             obj.topo = data;
%             obj.topomask = datamask;
%             
% %             if ltestplot
% %                 figure
% %                 title('after removing cells with only 1 neighbour')
% %                 imagesc(xf,yf,topo)
% %                 xlim([xh(1) xh(end)])
% %                 ylim([yh(1) yh(end)])
% %             end
%             %%
%             %if lhqplot
%                 %cd(outputdir)
%                 
% %                 h=figure;
% %                 hp1=subplot(1,2,1);
% %                 imagesc(obj.xf,obj.yf,flipud(topoinit))
% %                 set(gca,'YDir','normal','TickLabelInterpreter','latex')
% %                 xlim([obj.xh(1) obj.xh(end)])
% %                 ylim([obj.yh(1) obj.yh(end)])
% %                 caxis([0 45])
% %                 axis equal tight
% %                 xlabel('x [m]','Interpreter','latex','FontSize',12)
% %                 ylabel('y [m]','Interpreter','latex','FontSize',12)
% %                 hp1.XTick=[0 200 400 600 800 960];
% %                 hp1.YTick=[0 200 400 480];
% %                 set(gca,'TickLabelInterpreter','latex')
% %                 set(gca,'FontSize',12)
%                 %set(gca,XTick,[0 200 400 600 800 960])
%                 %set(gca,YTick,[0 200 400 480])
%                 %set(gca,'YTickLabel',[]);
%                 %set(gca,'XTickLabel',[]);
%                 
% %                 hp2=subplot(1,2,2);
% %                 imagesc(obj.xf,obj.yf,flipud(obj.topo))
% %                 set(gca,'YDir','normal','TickLabelInterpreter','latex')
% %                 xlim([obj.xh(1) obj.xh(end)])
% %                 ylim([obj.yh(1) obj.yh(end)])
% %                 caxis([0 45])
% %                 axis equal tight
% %                 xlabel('x [m]','Interpreter','latex','FontSize',12)
% %                 set(gca,'YTickLabel',[]);
% %                 hp2.XTick=[0 200 400 600 800 960];
% %                 hp1.Position=[0.09 0.1100 0.35 0.8150];
% %                 hp2.Position=[0.511 0.1100 0.35 0.8150];
% %                 hcb=colorbar('Position',[0.92 0.425 0.03 0.15]);
% %                 colormap(flipud(bone)) %colormap(flipud(gray))
% %                 hcb.Label.Interpreter='latex';
% %                 hcb.TickLabelInterpreter='latex';
% %                 title(hcb,'height [m]','Interpreter','latex','FontSize',12)
% %                 set(gca,'TickLabelInterpreter','latex')
% %                 set(gca,'FontSize',12)
% %                 set(gcf, 'Color', 'w');
% %                 
%                 %ylabel('y [m]','Interpreter','latex','FontSize',12)
%                 %
%                 %export_fig topocompbone -eps -png
%                 %print -depsc2 ICtopocompbone.eps
%                 %print -dpng ICtopocompbone.png
%                 %print -dpdf ICtopocompbone.pdf
%                 
%                 %%
% %                 h=figure;
% %                 set(gcf,'units','centimeters','position',[0 0 14.5 14.5]);
% %                 set(h,'PaperPosition',[0 0 14.5 14.5]);
% %                 set(h,'PaperUnits','centimeters');
% %                 
% %                 imagesc(obj.xf,obj.yf,flipud(topoinit))
% %                 set(gca,'YDir','normal','TickLabelInterpreter','latex')
% %                 xlim([obj.xh(1) obj.xh(end)])
% %                 ylim([obj.yh(1) obj.yh(end)])
% %                 caxis([0 45])
% %                 axis equal tight
% %                 xlabel('','Interpreter','latex','FontSize',10)
% %                 ylabel('y [m]','Interpreter','latex','FontSize',10)
% %                 h.Children.YTick=[0 200 400 480];
% %                 set(gca,'XTickLabel',[]);
% %                 set(gca,'TickLabelInterpreter','latex')
% %                 set(gca,'FontSize',12)
% %                 colormap(flipud(bone)) %colormap(flipud(gray))
% %                 hcb=colorbar('northoutside');
% %                 hcb.Label.Interpreter='latex';
% %                 hcb.TickLabelInterpreter='latex';
% %                 title(hcb,'height [m]','Interpreter','latex','FontSize',10)
% %                 set(gca,'TickLabelInterpreter','latex')
% %                 set(gca,'FontSize',12)
% %                 set(gcf, 'Color', 'w');
%                 %colormap(flipud(gray))
%                 
%                 %export_fig inittopobone -eps -png
%                 
%                 
% %                 h=figure;
% %                 set(gcf,'units','centimeters','position',[0 0 14.5 14.5]);
% %                 set(h,'PaperPosition',[0 0 14.5 14.5]);
% %                 set(h,'PaperUnits','centimeters');
% %                 imagesc(obj.xf,obj.yf,flipud(obj.topo))
% %                 set(gca,'YDir','normal','TickLabelInterpreter','latex')
% %                 xlim([obj.xh(1) obj.xh(end)])
% %                 ylim([obj.yh(1) obj.yh(end)])
% %                 caxis([0 45])
% %                 axis equal tight
% %                 h.Children.XTick=[0 200 400 600 800 960];
% %                 xlabel('x [m]','Interpreter','latex','FontSize',10)
% %                 ylabel('y [m]','Interpreter','latex','FontSize',10)
% %                 colormap(flipud(bone))
% %                 h.Children.YTick=[0 200 400 480];
% %                 %h.Position=[0.09 0.1100 0.35 0.8150];
% %                 %h.Position=[0.511 0.1100 0.35 0.8150];
% %                 %hcb=colorbar('Position',[0.92 0.425 0.03 0.15]);
% %                 set(gca,'TickLabelInterpreter','latex')
% %                 set(gca,'FontSize',10)
% %                 set(gcf, 'Color', 'w');
%                 
%                 %ylabel('y [m]','Interpreter','latex','FontSize',12)
%                 %
%                 %export_fig topobone -eps -png
%                 
%                 %cd(parentdir)
%             %end
%         end
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
da_pp.addvar(obj, 'topo', topo)
da_pp.addvar(obj, 'topomask', topomask)

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

% if ltestplot
%     figure
%     title('after removing cells with only 1 neighbour')
%     imagesc(xf,yf,topo)
%     xlim([xh(1) xh(end)])
%     ylim([yh(1) yh(end)])
% end
% %%
% if lhqplot
%     cd(outputdir)
%   
%     h=figure;
%     hp1=subplot(1,2,1);
%     imagesc(xf,yf,flipud(topoinit))
%     set(gca,'YDir','normal','TickLabelInterpreter','latex')
%     xlim([xh(1) xh(end)])
%     ylim([yh(1) yh(end)])
%     caxis([0 45])
%     axis equal tight
%     xlabel('x [m]','Interpreter','latex','FontSize',12)
%     ylabel('y [m]','Interpreter','latex','FontSize',12)
%     hp1.XTick=[0 200 400 600 800 960];
%     hp1.YTick=[0 200 400 480];
%     set(gca,'TickLabelInterpreter','latex')
%     set(gca,'FontSize',12)
%     %set(gca,XTick,[0 200 400 600 800 960])
%     %set(gca,YTick,[0 200 400 480])
%     %set(gca,'YTickLabel',[]);
%     %set(gca,'XTickLabel',[]);
%     
%     hp2=subplot(1,2,2);
%     imagesc(xf,yf,flipud(topo))
%     set(gca,'YDir','normal','TickLabelInterpreter','latex')
%     xlim([xh(1) xh(end)])
%     ylim([yh(1) yh(end)])
%     caxis([0 45])
%     axis equal tight
%     xlabel('x [m]','Interpreter','latex','FontSize',12)
%     set(gca,'YTickLabel',[]);
%     hp2.XTick=[0 200 400 600 800 960];
%     hp1.Position=[0.09 0.1100 0.35 0.8150];
%     hp2.Position=[0.511 0.1100 0.35 0.8150];
%     hcb=colorbar('Position',[0.92 0.425 0.03 0.15]);
%     colormap(flipud(bone)) %colormap(flipud(gray))
%     hcb.Label.Interpreter='latex';
%     hcb.TickLabelInterpreter='latex';
%     title(hcb,'height [m]','Interpreter','latex','FontSize',12)
%     set(gca,'TickLabelInterpreter','latex')
%     set(gca,'FontSize',12)
%     set(gcf, 'Color', 'w');
%     
%     %ylabel('y [m]','Interpreter','latex','FontSize',12)
%     %
%     export_fig topocompbone -eps -png
%     %print -depsc2 ICtopocompbone.eps
%     %print -dpng ICtopocompbone.png
%     %print -dpdf ICtopocompbone.pdf
%     
%     %%
%     h=figure;
%     set(gcf,'units','centimeters','position',[0 0 14.5 14.5]);
%     set(h,'PaperPosition',[0 0 14.5 14.5]);
%     set(h,'PaperUnits','centimeters');
% 
%     imagesc(xf,yf,flipud(topoinit))
%     set(gca,'YDir','normal','TickLabelInterpreter','latex')
%     xlim([xh(1) xh(end)])
%     ylim([yh(1) yh(end)])
%     caxis([0 45])
%     axis equal tight
%     xlabel('','Interpreter','latex','FontSize',10)
%     ylabel('y [m]','Interpreter','latex','FontSize',10)
%     h.Children.YTick=[0 200 400 480];
%     set(gca,'XTickLabel',[]);
%     set(gca,'TickLabelInterpreter','latex')
%     set(gca,'FontSize',12)
%     colormap(flipud(bone)) %colormap(flipud(gray))
%     hcb=colorbar('northoutside')
%     hcb.Label.Interpreter='latex';
%     hcb.TickLabelInterpreter='latex';
%     title(hcb,'height [m]','Interpreter','latex','FontSize',10)
%     set(gca,'TickLabelInterpreter','latex')
%     set(gca,'FontSize',12)
%     set(gcf, 'Color', 'w');
%      %colormap(flipud(gray))
% 
%   export_fig inittopobone -eps -png
% 
%     
%     h=figure;
%     set(gcf,'units','centimeters','position',[0 0 14.5 14.5]);
%     set(h,'PaperPosition',[0 0 14.5 14.5]);
%     set(h,'PaperUnits','centimeters');
%     imagesc(xf,yf,flipud(topo))
%     set(gca,'YDir','normal','TickLabelInterpreter','latex')
%     xlim([xh(1) xh(end)])
%     ylim([yh(1) yh(end)])
%     caxis([0 45])
%     axis equal tight
%     h.Children.XTick=[0 200 400 600 800 960];
%     xlabel('x [m]','Interpreter','latex','FontSize',10)
%     ylabel('y [m]','Interpreter','latex','FontSize',10)
%     colormap(flipud(bone))
%     h.Children.YTick=[0 200 400 480];
%     %h.Position=[0.09 0.1100 0.35 0.8150];
%     %h.Position=[0.511 0.1100 0.35 0.8150];
%     %hcb=colorbar('Position',[0.92 0.425 0.03 0.15]);
%     set(gca,'TickLabelInterpreter','latex')
%     set(gca,'FontSize',10)
%     set(gcf, 'Color', 'w');
%     
%     %ylabel('y [m]','Interpreter','latex','FontSize',12)
%     %
%    export_fig topobone -eps -png
% 
%  cd(parentdir)   
% end
end
        end
        
        
        function plot_bl(obj)
            figure
            title('Blocks (old)')
            view(52, 23)
            if (obj.lcastro || obj.lcube || obj.lblocks)
                for i = 1:size(obj.bl, 1)
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1)  obj.xh(obj.bl(i,1))], [obj.yh(obj.bl(i,3))  obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.zh(obj.bl(i,6)+1)  obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1)], [245 245 245] ./ 255)
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,3))  obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3))], [obj.bl(i,5)  obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], [245 245 245] ./ 255)
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.bl(i,5)  obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], [245 245 245] ./ 255)
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) ], [obj.yh(obj.bl(i,4)+1)  obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3))], [obj.bl(i,5)  obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], [245 245 245] ./ 255)
                    patch([obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.bl(i,5)  obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], [245 245 245] ./ 255)
                    % patch([xh(1) xh(end) xh(end)  xh(1)], [yh(1)  yh(1) yh(end) yh(end)], [zh(1)  zh(1) zh(1) zh(1)], [245 245 245] ./ 255)
                end
                
            elseif obj.lflat
                
            else
                for i = 1:size(obj.bl, 1)
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,1))], [obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1)], 'w')
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3))], [obj.bl(i,5) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], 'w')
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.bl(i,5) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], 'w')
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1))], [obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3))], [obj.bl(i,5) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], 'w')
                    patch([obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.bl(i,5) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], 'w')
                end
            end
            
            zlim([0 obj.zh(end)]); %/(r.blockheight-1))
            xlim([0 obj.xh(end)]); %/(r.blockheight-1))
            ylim([0 obj.yh(end)]); %/(r.blockheight-1))

            set(gca,'ticklabelinterpreter','latex')
            xlabel('x [m]','interpreter','latex')
            ylabel('y [m]','interpreter','latex')
            zlabel('z [m]','interpreter','latex')
            set(gca,'BoxStyle','full','Box','on')
            daspect([1 1 1])
            
        end
        
        
        function plot_blocks(obj)
            figure
            title('Blocks')
            view(52, 23)
        for i = 1:obj.nblockstotal
            il = obj.blocks(i,1);
            iu = obj.blocks(i,2);
            jl = obj.blocks(i,3);
            ju = obj.blocks(i,4);
            kl = obj.blocks(i,5);
            ku = obj.blocks(i,6);
            
            if i <= obj.nblocks
%             if il == 0
%                 il = il + 1;
%             end
%             if jl == 0
%                 jl = jl + 1;
%             end
%             if kl == 0
%                 kl = kl + 1;
%             end
%             if ku == 0
%                 ku = ku + 1;
%             end
                           
            patch([obj.xh(il)   obj.xh(iu+1) obj.xh(iu+1) obj.xh(il)]  , [obj.yh(jl)   obj.yh(jl)   obj.yh(ju+1) obj.yh(ju+1)], [obj.zh(ku+1) obj.zh(ku+1) obj.zh(ku+1) obj.zh(ku+1)], [245 245 245] ./ 255)
            patch([obj.xh(il)   obj.xh(il)   obj.xh(iu+1) obj.xh(iu+1)], [obj.yh(jl)   obj.yh(jl)   obj.yh(jl)   obj.yh(jl)],   [obj.zh(kl)   obj.zh(ku+1) obj.zh(ku+1) obj.zh(kl)], [245 245 245] ./ 255)
            patch([obj.xh(il)   obj.xh(il)   obj.xh(iu+1) obj.xh(iu+1)], [obj.yh(ju+1) obj.yh(ju+1) obj.yh(ju+1) obj.yh(ju+1)], [obj.zh(kl)   obj.zh(ku+1) obj.zh(ku+1) obj.zh(kl)], [245 245 245] ./ 255)
            patch([obj.xh(il)   obj.xh(il)   obj.xh(il)   obj.xh(il)]  , [obj.yh(ju+1) obj.yh(ju+1) obj.yh(jl)   obj.yh(jl)],   [obj.zh(kl)   obj.zh(ku+1) obj.zh(ku+1) obj.zh(kl)], [245 245 245] ./ 255)
            patch([obj.xh(iu+1) obj.xh(iu+1) obj.xh(iu+1) obj.xh(iu+1)], [obj.yh(jl)   obj.yh(jl)   obj.yh(ju+1) obj.yh(ju+1)], [obj.zh(kl)   obj.zh(ku+1) obj.zh(ku+1) obj.zh(kl)], [245 245 245] ./ 255)
            
            else
            patch([obj.xh(il)   obj.xh(iu+1) obj.xh(iu+1) obj.xh(il)]  , [obj.yh(jl)   obj.yh(jl)   obj.yh(ju+1) obj.yh(ju+1)], [obj.zh(ku+1) obj.zh(ku+1) obj.zh(ku+1) obj.zh(ku+1)], [245 245 245] ./ 255)
            patch([obj.xh(il)   obj.xh(il)   obj.xh(iu+1) obj.xh(iu+1)], [obj.yh(jl)   obj.yh(jl)   obj.yh(jl)   obj.yh(jl)],   [0            obj.zh(ku+1) obj.zh(ku+1)            0], [245 245 245] ./ 255)
            patch([obj.xh(il)   obj.xh(il)   obj.xh(iu+1) obj.xh(iu+1)], [obj.yh(ju+1) obj.yh(ju+1) obj.yh(ju+1) obj.yh(ju+1)], [0            obj.zh(ku+1) obj.zh(ku+1)            0], [245 245 245] ./ 255)
            patch([obj.xh(il)   obj.xh(il)   obj.xh(il)   obj.xh(il)]  , [obj.yh(ju+1) obj.yh(ju+1) obj.yh(jl)   obj.yh(jl)],   [0            obj.zh(ku+1) obj.zh(ku+1)            0], [245 245 245] ./ 255)
            patch([obj.xh(iu+1) obj.xh(iu+1) obj.xh(iu+1) obj.xh(iu+1)], [obj.yh(jl)   obj.yh(jl)   obj.yh(ju+1) obj.yh(ju+1)], [0            obj.zh(ku+1) obj.zh(ku+1)            0], [245 245 245] ./ 255)    
            end
                
        end
            
            zlim([0 obj.zh(end)]); %/(r.blockheight-1))
            xlim([0 obj.xh(end)]); %/(r.blockheight-1))
            ylim([0 obj.yh(end)]); %/(r.blockheight-1))

            set(gca,'ticklabelinterpreter','latex')
            xlabel('x [m]','interpreter','latex')
            ylabel('y [m]','interpreter','latex')
            zlabel('z [m]','interpreter','latex')
            set(gca,'BoxStyle','full','Box','on')
            daspect([1 1 1])
        end
        
        function generate_facets(obj)
            da_pp.block2fac(obj)
            da_pp.addvar(obj, 'nboundingwallfacets', 0)
            if obj.lEB
                da_pp.addboundingwalls(obj)
            end
            da_pp.createfloors(obj);
        end
              
        
        function plot_facets(obj)
            figure
            cmap = colormap('parula');
            top = 1; west = 2; east = 3; north = 4; south = 5; bot = 6;
            for i = 1:obj.nfcts
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
                %ci = min(floor(double(obj.facets(i, 2)) / double(max(obj.facets(:, 2))) * length(cmap)) + 1, length(cmap))
                ci = 64;
                patch(x,y,z, cmap(ci, :),'FaceLighting','none');
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
                
                text(mean(x) + d(1), mean(y) + d(2), mean(z) + d(3), num2str(i), 'horizontalalignment', 'center')
                hold on
                title('Facet number')
            end
            view(3)
            xlabel('x')
            ylabel('y')
            zlabel('z')
            axis equal
            xlim([0 obj.xh(end)])
            ylim([0 obj.yh(end)])
            zlim([0 obj.zh(end)])
        end
            

                       
            %disp(obj.facets)
%             figure
%             cmap = colormap('parula');
%             top = 1; west = 2; east = 3; north = 4; south = 5; bot = 6;
%             il = 1; iu = 2; jl = 3; ju = 4; kl = 5; ku = 6;
%             for i = 1:obj.nfcts
%                 switch obj.facets(i, 1)
%                     case {top, bot}
%                         x = [obj.xh(obj.facets(i, 5+il)), obj.xh(obj.facets(i, 5+iu)), obj.xh(obj.facets(i, 5+iu)), obj.xh(obj.facets(i, 5+il))];
%                         y = [obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+ju)), obj.yh(obj.facets(i, 5+ju))];
%                         z = [obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+kl))];
%                     case {west, east}
%                         x = [obj.xh(obj.facets(i, 5+il)), obj.xh(obj.facets(i, 5+il)), obj.xh(obj.facets(i, 5+il)), obj.xh(obj.facets(i, 5+il))];
%                         y = [obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+ju)), obj.yh(obj.facets(i, 5+ju))];
%                         z = [obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+ku)), obj.zh(obj.facets(i, 5+ku)), obj.zh(obj.facets(i, 5+kl))];
%                     case {north, south}
%                         x = [obj.xh(obj.facets(i, 5+il)), obj.xh(obj.facets(i, 5+iu)), obj.xh(obj.facets(i, 5+iu)), obj.xh(obj.facets(i, 5+il))];
%                         y = [obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+jl))];
%                         z = [obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+ku)), obj.zh(obj.facets(i, 5+ku))];
%                 end
%                 
%                 subplot(1, 3, 1)
%                 ci = min(floor(double(obj.facets(i, 1)) / 6 * length(cmap)) + 1, length(cmap));
%                 patch(x,y,z, cmap(ci, :),'FaceLighting','none');
%                 hold on
%                 title('Orientation')
%                 
%                 subplot(1, 3, 2)
%                 ci = min(floor(double(obj.facets(i, 3)) / obj.nblocks * length(cmap)) + 1, length(cmap));
%                 patch(x,y,z, cmap(ci, :),'FaceLighting','none');
%                 hold on
%                 title('Blocks')
%                 
%                 subplot(1, 3, 3)
%                 ci = min(floor(double(obj.facets(i, 4))/ obj.nbuildings * length(cmap)) + 1, length(cmap));
%                 patch(x,y,z, cmap(ci, :),'FaceLighting','none');
%                 hold on
%                 title('Buildings')
%             end
            
%             for n =1:3
%                 subplot(1,3,n)
%                 view(3)
%                 xlabel('x')
%                 ylabel('y')
%                 zlabel('z')
%                 axis equal
%                 xlim([0 obj.xh(end)])
%                 ylim([0 obj.yh(end)])
%                 zlim([0 obj.zh(end)])
%             end
%             
%             
%             % Plot building map
%             figure
%             subplot(1,2,1)
%             xlabel('x')
%             ylabel('y')
%             axis equal
%             title('Building id')
%             for i = find(obj.facets(:,1) == top)'
%                 x = [obj.xh(obj.facets(i, 5+il)), obj.xh(obj.facets(i, 5+iu)), obj.xh(obj.facets(i, 5+iu)), obj.xh(obj.facets(i, 5+il))];
%                 y = [obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+ju)), obj.yh(obj.facets(i, 5+ju))];
%                 
%                 ci = min(floor(double(obj.facets(i, 4)) / obj.nbuildings * length(cmap)) + 1, length(cmap));
%                 patch(x, y, cmap(ci, :),'FaceLighting','none');
%                 hold on
%                 text(mean(x), mean(y), num2str(obj.facets(i, 4), '%8d'), ...
%                     'horizontalalignment', 'center')
%             end
%            
%             subplot(1,2,2);
%             for i = 1:obj.nfcts
%                 switch obj.facets(i, 1)
%                     case {top, bot}
%                         x = [obj.xh(obj.facets(i, 5+il)), obj.xh(obj.facets(i, 5+iu)), obj.xh(obj.facets(i, 5+iu)), obj.xh(obj.facets(i, 5+il))];
%                         y = [obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+ju)), obj.yh(obj.facets(i, 5+ju))];
%                         z = [obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+kl))];
%                     case {west, east}
%                         x = [obj.xh(obj.facets(i, 5+il)), obj.xh(obj.facets(i, 5+il)), obj.xh(obj.facets(i, 5+il)), obj.xh(obj.facets(i, 5+il))];
%                         y = [obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+ju)), obj.yh(obj.facets(i, 5+ju))];
%                         z = [obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+ku)), obj.zh(obj.facets(i, 5+ku)), obj.zh(obj.facets(i, 5+kl))];
%                     case {north, south}
%                         x = [obj.xh(obj.facets(i, 5+il)), obj.xh(obj.facets(i, 5+iu)), obj.xh(obj.facets(i, 5+iu)), obj.xh(obj.facets(i, 5+il))];
%                         y = [obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+jl))];
%                         z = [obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+ku)), obj.zh(obj.facets(i, 5+ku))];
%                 end
%                 
%                 ci = min(floor(double(obj.facets(i, 2)) / double(max(obj.facets(:, 2))) * length(cmap)) + 1, length(cmap));
%                 patch(x,y,z, cmap(ci, :),'FaceLighting','none');
%                 d = [0, 0, 0]; a = 0.25;
%                 switch(obj.facets(i, 1))
%                     case top
%                         d(3) = a * obj.dz(1);
%                     case west
%                         d(1) = -a * obj.dx(1);
%                     case east
%                         d(1) = a * obj.dx(1);
%                     case south
%                         d(2) = -a * obj.dy(1);
%                     case north
%                         d(2) = a * obj.dy(1);
%                 end
%                 text(mean(x) + d(1), mean(y) + d(2), mean(z) + d(3), num2str(obj.facets(i, 2)), 'horizontalalignment', 'center')
%                 hold on
%                 title('Wall type')
%             end
%             
%             view(3)
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             axis equal
%             xlim([0 obj.xh(end)])
%             ylim([0 obj.yh(end)])
%             zlim([0 obj.zh(end)])
%             % colorbar
%             % return
%             
%             % Plot facets
%             figure
%             for i = 1:obj.nfcts
%                 switch obj.facets(i, 1)
%                     case {top, bot}
%                         x = [obj.xh(obj.facets(i, 5+il)), obj.xh(obj.facets(i, 5+iu)), obj.xh(obj.facets(i, 5+iu)), obj.xh(obj.facets(i, 5+il))];
%                         y = [obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+ju)), obj.yh(obj.facets(i, 5+ju))];
%                         z = [obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+kl))];
%                     case {west, east}
%                         x = [obj.xh(obj.facets(i, 5+il)), obj.xh(obj.facets(i, 5+il)), obj.xh(obj.facets(i, 5+il)), obj.xh(obj.facets(i, 5+il))];
%                         y = [obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+ju)), obj.yh(obj.facets(i, 5+ju))];
%                         z = [obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+ku)), obj.zh(obj.facets(i, 5+ku)), obj.zh(obj.facets(i, 5+kl))];
%                     case {north, south}
%                         x = [obj.xh(obj.facets(i, 5+il)), obj.xh(obj.facets(i, 5+iu)), obj.xh(obj.facets(i, 5+iu)), obj.xh(obj.facets(i, 5+il))];
%                         y = [obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+jl)), obj.yh(obj.facets(i, 5+jl))];
%                         z = [obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+kl)), obj.zh(obj.facets(i, 5+ku)), obj.zh(obj.facets(i, 5+ku))];
%                 end
%                 
%                 ci = min(floor(double(obj.facets(i, 2)) / double(max(obj.facets(:, 2))) * length(cmap)) + 1, length(cmap));
%                 patch(x,y,z, cmap(ci, :),'FaceLighting','none');
%                 d = [0, 0, 0]; a = 0.25;
%                 switch(fctl(i, 1))
%                     case top
%                         d(3) = a * obj.dz(1);
%                     case west
%                         d(1) = -a * obj.dx(1);
%                     case east
%                         d(1) = a * obj.dx(1);
%                     case south
%                         d(2) = -a * obj.dy(1);
%                     case north
%                         d(2) = a * obj.dy(1);
%                 end
%                 text(mean(x) + d(1), mean(y) + d(2), mean(z) + d(3), num2str(i), ...
%                     'horizontalalignment', 'center')
%                 hold on
%                 title('facet nr')
%             end
%             view(3)
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             axis equal
%             xlim([0 obj.xh(end)])
%             ylim([0 obj.yh(end)])
%             zlim([0 obj.zh(end)])
        
        
        function generate_EB(obj, ltestplot, lhqplot, lwritefile)
            %da_pp.addvar(obj, 'walltypes', dlmread(['walltypes.inp.', obj.expnr],'',3,0));
            da_pp.vsolc(obj, ltestplot, lhqplot)
            da_pp.vfc(obj, lwritefile)
            da_pp.rayit(obj, ltestplot, lhqplot, lwritefile)
            da_pp.generate_Tfacinit(obj, ltestplot, lwritefile)
        end
        
        function makeblocks(obj)
            topomask = obj.topomask;
            topo = obj.topo;
            imax = obj.imax;
            jtot = obj.jtot;
            dz = obj.dz;
            
            maxnrblocks = sum(topomask(:)); %allocate arrays with maximum size they can possibly have, reduce size later
            xmin = zeros(maxnrblocks,1); %store lower x bound of blocks
            xmax = zeros(maxnrblocks,1); %store upper x bound of blocks
            ymin = zeros(maxnrblocks,1); %store lower y bound of blocks
            ymax = zeros(maxnrblocks,1); %store upper y bound of blocks
            zmin = zeros(maxnrblocks,1); %store lower z bound of blocks
            zmax = zeros(maxnrblocks,1); %store upper z bound of blocks
            
            blchecked=zeros(size(topomask)); %mask of blocks which have already been checked
            indexmask=zeros(size(topomask)); %mask of the indeces of the (new) blocks
            
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
                        heightgradient = diff(topo(j:end, i));                       
                        if isempty(heightgradient) % if at the end of the y-direction (je)
                            heightchangey = j;
                        elseif all(heightgradient == 0) % same block/ floor until the end of the domain
                            heightchangey = jtot;
                        else
                            heightchangey = find(heightgradient ~=0 , 1) + j - 1;  %last cell with same height as j (i.e. there is a height change betwen this and the next cell
                        end
                        for jj = yminl+1:heightchangey                          
                            if (i == 1)
                                if ((any(topo(jj-1:jj,i+1) > 0)) && (topo(jj, i+1) ~= topo(jj-1,i+1)))
                                    heightchangey = jj - 1;
                                    break                                  
                                end                           
                            elseif (i == imax)                               
                                if (any(topo(jj-1:jj,i-1) > 0) && (topo(jj,i-1)~= topo(jj-1,i-1) ) )
                                    heightchangey = jj-1;
                                    break                                   
                                end                               
                            else                                
                                if (any(topo(jj-1:jj,i-1) > 0) && (topo(jj,i-1) ~= topo(jj-1,i-1))) || ((any(topo(jj-1:jj,i+1)>0)) && (topo(jj,i+1)~= topo(jj-1,i+1)))                               
                                    heightchangey = jj-1; % overwrite heightchangey as we need to truncate block earlier!                                   
                                    break                                    
                                end                                
                            end                           
                        end
                        
                        ymaxl = heightchangey;  % end of the current block (either floor or block with different height comes next)                        
                        ztemp = zeros(6,1);
                        
                        % tg3315 changed from commented below as can have yminl==1 and
                        % ymaxl==nj
                        if yminl == 1
                            ztemp(2,1) = NaN;
                        else
                            ztemp(2,1) = topo(yminl-1,xminl);
                        end
                        
                        if xminl == 1
                            ztemp(3,1) = NaN;
                        else
                            ztemp(3,1) = topo(yminl, xminl-1);
                        end
                        
                        if xmaxl == imax
                            ztemp(4,1) = NaN;
                        else
                            ztemp(4,1) = topo(yminl, xminl+1);
                        end
                        
                        if ymaxl == jtot
                            ztemp(5,1) = NaN;
                        else
                            ztemp(5,1) = topo(yminl+1, xminl);
                        end
                        zcuttemp = sort(ztemp(ztemp<topo(yminl, xminl) & ztemp>0));
                        zcut = [0; zcuttemp; topo(yminl,xminl)];
                        for kc = 1:size(zcut, 1)-1                           
                            count = count+1;                            
                            blchecked(yminl:ymaxl,xminl:xmaxl) = blchecked(yminl:ymaxl,xminl:xmaxl)+1;
                            indexmask(yminl:ymaxl,xminl:xmaxl) = count;
                            xmin(count) = xminl;
                            xmax(count) = xmaxl;
                            ymin(count) = yminl;
                            ymax(count) = ymaxl;
                            if zcut(kc,1) == 0
                                zmin(count) = 0;
                            else
                                zmin(count) = zcut(kc,1) / dz + 1;
                            end
                            zmax(count) = zcut(kc+1,1) / dz;                            
                        end                       
                        j = heightchangey + 1;  % move to index after current block                       
                    end
                end
            end
%             if ltestplot
%                 figure
%                 imagesc(xf,yf,indexmask)
%                 set(gca,'YDir','normal')
%                 xlim([xh(1) xh(end)])
%                 ylim([yh(1) yh(end)])
%                 title('after cutting into y slices')
%             end
            
            %shorten arrays
            ymin((count+1):end)=[];
            ymax((count+1):end)=[];
            xmin((count+1):end)=[];
            xmax((count+1):end)=[];
            zmin((count+1):end)=[];
            zmax((count+1):end)=[];
            
            %disp(['Number of blocks after cutting into y slices: ' num2str(count)])
            
            
            
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
                        if all(ymax2(b3) == ymax2(i)) && all(zmax2(b3) == zmax2(i)) && topo(ymin2(b3),xmin2(b3)) == topo(ymin2(i), xmin2(i)) %&& all(zmin2(b2)==zmin2(i)) %if they also have the same upper y bound and the same height
                            % additional check to make sure we do not merge blocks in a way that causes internal-external facets at this point
                            if (topo(max(1,ymin2(b3)-1),xmax2(b3)) == topo(max(1,ymin2(i)-1),xmax2(i))) && (topo(min(ymax2(b3)+1, jtot),xmax2(b3)) == topo(min(ymax2(i)+1, jtot), xmax2(i)))
                                xmax2(i) = xmax2(b3); %merge
                                xmax2(b3) = []; ymax2(b3) = []; zmax2(b3)  =[]; xmin2(b3) = []; ymin2(b3) = []; zmin2(b3) = []; %remove the just merged block from list of blocks
                            end
                        end
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
            %disp(['Number of blocks after merging y slices of same size along x: ' num2str(count2)])
            
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
                indexmask2(ymin2(i):ymax2(i),xmin2(i):xmax2(i)) = i;
                datamean2(ymin2(i):ymax2(i),xmin2(i):xmax2(i)) = zmax2(i);
            end
            
%             if ltestplot
%                 figure
%                 imagesc(xf,yf,indexmask2)
%                 set(gca,'YDir','normal')
%                 xlim([xh(1) xh(end)])
%                 ylim([yh(1) yh(end)])
%                 title('after merging y slices of same size along x')
%             end
            
            %% tg and bs 30.05.19 - try to split blocks vertically too
            
            % xmax25=xmax2; ymax25=ymax2; zmax25=zmax2;
            % xmin25=xmin2; ymin25=ymin2; zmin25=zmin2;
            %
            % for i=1:ni %loop over all x
            
            
            %% At this point can form block file that works for old cold
            
            % Edit for NY in 1930s topology
            % dummy=ones(5,count2);
            %
            % for k=1:length(zmax2)
            %     if zmax2(k)>12
            %         zmax2(k) = 12+round(12*rand(1));
            %     end
            % end
            %
            % fileID = fopen([outputdir '/blocksOLD.inp.' num2str(expnr)],'w');
            % fprintf(fileID,'# Block data\n');
            % fprintf(fileID,'#  il\t   iu\t   jl\t   ju\t   kl\t   ku\t dtop\t dwest\t deast\t dnor\t dsou\n');
            % fprintf(fileID,'%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\n',[xmin2';xmax2';ymin2';ymax2';zmin2';zmax2';dummy]);
            % fclose(fileID);
            
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
                %disp(['Number of blocks after applying rules for radiation in x and y : ' num2str(count3)])
                ymin3((count3+1):end) = [];
                ymax3((count3+1):end) = [];
                xmin3((count3+1):end) = [];
                xmax3((count3+1):end) = [];
                zmin3((count3+1):end) = [];
                zmax3((count3+1):end) = [];
                
                
%                 if ltestplot
%                     figure
%                     imagesc(indexmask3)
%                     hold on
%                     for i = 1:count3-1
%                         rectangle('Position',[xmin3(i)-0.5 ymin3(i)-0.5 xmax3(i)-xmin3(i)+1 ymax3(i)-ymin3(i)+1])
%                     end
%                     axis equal tight
%                     hold off
%                     title('after applying rules for radiation in x and y')
%                     
%                 end                               
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
%             if ltestplot
%                 figure
%                 imagesc(datamean4)
%             end
            
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
             
%             if ltestplot
%                 figure
%                 imagesc(indexmask4.*internalmask)
%                 hold on
%                 for i=1:count3-1
%                     rectangle('Position',[xmin3(i)-0.5 ymin3(i)-0.5 xmax3(i)-xmin3(i)+1 ymax3(i)-ymin3(i)+1])
%                 end
%                 axis equal tight
%                 hold off
%                 title('internal blocks after applying rules for radiation in x and y')
%             end
            
            % tg3315 may need to develop this further for geometries with VERY uneven
            % roof tops... Currently will not allow internal blocks below but this does
            % satisy the radiation laws in ther vertical direction...
            
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
            %disp(['Number of blocks after merging internal blocks along y: ' num2str(length(xmin5))])
            
%             if ltestplot
%                 figure
%                 imagesc(indexmask5.*internalmask)
%                 hold on
%                 for i=1:count5
%                     rectangle('Position',[xmin5(i)-0.5 ymin5(i)-0.5 xmax5(i)-xmin5(i)+1 ymax5(i)-ymin5(i)+1])
%                 end
%                 axis equal tight
%                 hold off
%                 title('Internal blocks after merging internal blocks')
%             end
            
            %make blocks
            datamean5 = zeros(size(topomask));
            dataind = zeros(size(topomask));
            for i = 1:count5
                datamean5(ymin5(i):ymax5(i),xmin5(i):xmax5(i)) = zmax5(i);
                dataind(ymin5(i):ymax5(i),xmin5(i):xmax5(i)) = i;
            end
                       
%             if ltestplot
%                 figure
%                 subplot(1,3,1)
%                 imagesc(indexmask5)
%                 hold on
%                 for i=1:count5
%                     rectangle('Position',[xmin5(i)-0.5 ymin5(i)-0.5 xmax5(i)-xmin5(i)+1 ymax5(i)-ymin5(i)+1])
%                 end
%                 axis equal tight
%                 hold off
%                 title('index: blocks after merging merging internal blocks')
%                 set(gca,'YDir','normal')
%                                
%                 subplot(1,3,2)
%                 imagesc(datamean5)
%                 axis equal tight
%                 hold off
%                 title('height: blocks after merging merging internal blocks')
%                 set(gca,'YDir','normal')
%                                 
%                 subplot(1,3,3)
%                 imagesc(dataind)
%                 axis equal tight
%                 hold off
%                 title('index: blocks after merging merging internal blocks')
%                 set(gca,'YDir','normal')
%             end
                        
            ymin6 = ymin5;
            ymax6 = ymax5;
            xmin6 = xmin5;
            xmax6 = xmax5;
            zmin6 = zmin5;
            zmax6 = zmax5;
            
%             if lhqplot
%                 cd(outputdir)
%                 h=figure;
%                 set(gcf,'units','centimeters','position',[0 0 14.5 14.5]);
%                 set(h,'PaperPosition',[0 0 14.5 14.5]);
%                 set(h,'PaperUnits','centimeters');
%                 imagesc(xf,yf,flipud(datamean2))
%                 set(gca,'YDir','normal','TickLabelInterpreter','latex')
%                 h1=gca;
%                 h1.Position=[0.08 0.1100 0.78 0.8150];
%                 xlim([xh(1) xh(end)])
%                 ylim([yh(1) yh(end)])
%                 caxis([0 1])
%                 xlabel('x [m]','Interpreter','latex','FontSize',12)
%                 ylabel('y [m]','Interpreter','latex','FontSize',12)
%                 
%                 colormap(flipud(bone)) %colormap(flipud(gray))
%                 hcb=colorbar('Position',[0.92 0.15 0.03 0.69]);
%                 hcb.Label.Interpreter='latex';
%                 hcb.TickLabelInterpreter='latex';
%                 title(hcb,'height [m]','Interpreter','latex','FontSize',12)
%                 
%                 hold on
%                 for i=1:count5-1
%                     rectangle('Position',[2*(xmin6(i))-4 960-2*(ymin6(i))+4-2*(ymax6(i)-ymin6(i)+1) 2*(xmax6(i)-xmin6(i)+1) 2*(ymax6(i)-ymin6(i)+1)])
%                 end
%                 axis equal tight
%                 
%                 hold off
%                 
%                 print -depsc2 topo.eps
%                 print -dpng topo.png
%                 cd(parentdir)
%             end
            
            %% flip the whole y-dimension!!
            
            ymin6f = jtot - ymax6 + 1;
            ymax6f = jtot - ymin6 + 1;
            ymin2f = jtot - ymax2 + 1;
            ymax2f = jtot - ymin2 + 1;
            
            %%
            % wtmean5=ones(5,length(xmin6));
            %
            % fileID = fopen([tempdir '/blocksmeannointernalflip.inp'],'w');
            % fprintf(fileID,'# Block data\n');
            % fprintf(fileID,'# block indices,                     wall type\n');
            % fprintf(fileID,'#  il\t   iu\t   jl\t   ju\t   kl\t   ku\t ttop\t twest\t teast\t tnor\t tsou\n');
            % fprintf(fileID,'%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\n',[xmin6';xmax6';ymin6f';ymax6f';zmin6';zmax6';wtmean5]);
            % fclose(fileID);
            %%
            % wtmean=ones(5,count2);
            %
            % fileID = fopen([tempdir '/bbri.inp'],'w');
            % fprintf(fileID,'# Block data\n');
            % fprintf(fileID,'# block indices,                     wall type\n');
            % fprintf(fileID,'#  il\t   iu\t   jl\t   ju\t   kl\t   ku\t ttop\t twest\t teast\t tnor\t tsou\n');
            % fprintf(fileID,'%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\n',[xmin2';xmax2';ymin2f';ymax2f';zmin2';zmax2';wtmean]);
            % fclose(fileID);
            
                      
            %write block file to extend with roads and bounding walls etc            
            dummy = ones(count5, 5);
            
            %fileID = fopen([tempdir '/blocks.inp.' num2str(expnr)],'w');
            %fprintf(fileID,'# Block data\n');
            %fprintf(fileID,'#  il\t   iu\t   jl\t   ju\t   kl\t   ku\t dtop\t dwest\t deast\t dnor\t dsou\n');
            %fprintf(fileID,'%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\n',[xmin5';xmax5';ymin6';ymax6';zmin5';zmax5';dummy]);
            %fclose(fileID);
            
            blocks = [xmin5 xmax5 ymin6 ymax6 zmin5 zmax5 dummy];
           
            
            %write block file to use for test of ray intersection (without roads that will be added to blocks.inp)
            %already raise by 1 (because of roads, is done to blocks.inp.xxx later)
            %give each building a unique number, might differ from results in block2fac
            %but it doesn't matter
            
            blockindexmask = dataind;
            nbl = size(xmin5, 1);
            buildingindexmask = bwlabel(blockindexmask);
            buildingindexlist = zeros(1, nbl);
            
%             if ltestplot
%                 figure
%                 imagesc(buildingindexmask)
%                 axis equal tight
%                 hold off
%                 title('building indeces')
%                 set(gca,'YDir','normal')
%             end
            
            xuu = xmax5;
            yuu = ymax6;
            zuu = zmax5;
            %end
            
            for i = 1:nbl
                buildingindexlist(i) = buildingindexmask(yuu(i), xuu(i));
            end
            
            obj.nblocks = size(blocks, 1);
            obj.blocks = blocks;
            da_pp.addvar(obj, 'buildings', [xmin5  xmax5  ymin6 ymax6 zmin5+1 zmax5+1 buildingindexlist']);            
        end
        
        function block2fac(obj)
            imax = obj.imax;
            jtot = obj.jtot;
            %% read blocks
            %nheader = 2;
%             try %in case file is empty -> no blocks
%                 blocks = obj.blocks;
%             catch
%                 blocks = [];
%             end
            blocks = obj.blocks;
            
            nblocks = obj.nblocks;
            %nblks = size(blocks, 1);
            %if nblocks > 0
            blocks(:, [5,6]) = blocks(:,[5,6]) + 1; % indices for k start at zero?
            %end
            %obj.blocks = blocks;
            disp(blocks)
                            
            %% create Mask-matrix
            % this new mask is in x,y coordinates, not y,x coordinates as before
            M = zeros(imax, jtot);
            IM = zeros(imax , jtot);
            for i = 1:size(blocks,1)
                xl = blocks(i,1);
                xu = blocks(i,2);
                yl = blocks(i,3);
                yu = blocks(i,4);
                M(xl:xu, yl:yu) = 1;
                IM(xl:xu, yl:yu) = i;
            end
            
            % figure
            % imagesc(M)
            % set(gca,'YDir','normal')
            
            %% define all facets
            
            % fctl format: orientation, walltype, blockid, buildingid, isinternal, il, iu, jl, ju, kl, ku
            
            top = 1; west = 2; east = 3; north = 4; south=5; bot = 6;
            
            il = 1; iu = 2;
            jl = 3; ju = 4;
            kl = 5; ku = 6;
            
            nfcts = nblocks * 6;
            facets = int32(zeros(nfcts, 10));
            for j = 1:nblocks
                for k = top:bot
                    i = (j - 1) * 6 + k;
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
            
            
            
            %% alternative after merging internal blocks
            %disp('determine internal facets')
            intern = zeros(nfcts, 2);
            j = 1;
            for i = 1:nfcts
                switch(facets(i, 1))
                    %only one test necessary, since the whole edge is internal, or not
                    %remember that facets stores facet coordinates, not block coordinates
                    case 2
                        if facets(i, 5 + il) - 1 >= 1 %not at the domain edge
                            if M(facets(i, 5 + il) - 1 , facets(i, 5 + jl))
                                facets(i,5) = 1;
                                intern(j,1) = facets(i, 3);
                                intern(j,2) = IM(facets(i, 5 + il) - 1 ,facets(i, 5 + jl));
                                j = j + 1;
                            end
                        end
                    case 3
                        if facets(i, 5 + iu) <= imax %not at the domain edge
                            if M(facets(i, 5 +  iu), facets(i, 5 + jl))
                                facets(i,5) = 1;
                                intern(j,1) = facets(i, 3);
                                intern(j,2) = IM(facets(i, 5 + iu), facets(i, 5 + jl));
                                j = j + 1;
                            end
                        end
                    case 4
                        if facets(i, 5 + ju) <= jtot %not at the domain edge
                            if M(facets(i, 5 + il), facets(i, 5 + ju))
                                facets(i,5) = 1;
                                intern(j,1) = facets(i, 3);
                                intern(j,2) = IM(facets(i, 5+il), facets(i, 5+ju));
                                j = j + 1;
                            end
                        end
                    case 5
                        if facets(i, 5 + jl) - 1 >= 1 %not at the domain edge
                            if M(facets(i, 5 + il), facets(i, 5 + jl) - 1)
                                facets(i, 5) = 1;
                                intern(j, 1) = facets(i, 3);
                                intern(j, 2) = IM(facets(i, 5+il), facets(i, 5+jl) - 1);
                                j = j + 1;
                            end
                        end
                end
            end
            nintern = sum(facets(:,5));
            intern = intern(1:nintern,:);
            
            %disp('remove downward and internal facets')
            bblk0 = facets(:, 3);
            bblk1 = bblk0;
            
            for n=1:nintern
                % find all facets belonging to block associated with the facet in
                % the second column of intern, and replace these with the block id
                % from the facet in the first column (contagion)
                list=find(bblk0==intern(n,2));
                list2=find(bblk1==bblk1(list(1)));
                bblk1(list2)=intern(n,1);
            end            
            
            % assign building id to facets
            bblku = unique(bblk1);
            nbld = length(bblku);
            da_pp.addvar(obj, 'nbuildings', length(bblku));
            for n = 1:nbld
                facets(bblk1 == bblku(n), 4) = n;
            end
            
            %% remove downward facets and internal facets
            % facets format: orientation, walltype, blockid, buildingid, isinternal
            %              il, iu, jl, ju, kl, ku
            sel = find(facets(:,1) ~= bot & facets(:,5) ~= 1);
            facets = facets(sel, :);
            nfcts=size(facets,1);
            %% remove facets at domain edge (not actually done)
            
            sel = find(~((facets(:,6) == 1 & facets(:,7) == 1) | (facets(:,6)== imax + 1 & facets(:,7) == imax + 1) | (facets(:,8) == 1 & facets(:,9) == 1) | (facets(:,8) == jtot + 1 & facets(:,9) == jtot + 1)));
            
            if obj.lEB && (length(sel) ~= nfcts)
                  error("Can't have blocks on edge of domain when using energy balance")
            end

            
            obj.nfcts = size(facets,1);
            da_pp.addvar(obj, 'nblockfcts', size(facets, 1))
             
            
%             for i = 1:obj.nblockfcts
%                 obj.blocks(facets(i, 3), facets(i, 1) + 6) = i;
%             end
            
            %obj.blocks = blocks;
            obj.facets = facets;
            %disp(obj.facets)
        end
        
        function addboundingwalls(obj)
            % Add a bounding wall around the domain used in radiation calculations only
            % THIS DOES NOT YET DEAL WITH PERIODIC GEOMETRY (I.E. CANYONS)
                                   
            %fctl = obj.facets; % Need to only use obj.facets eventually. 
            
            imax = obj.imax;
            jtot = obj.jtot;
            maxsize = obj.maxsize;
            blocks = obj.blocks;
            height = floor(median(blocks(:, 6)));
            %nboundingwallfacets = obj.nboundingwallfacets;
            
            
            %if obj.nblocks > 0
                if any(blocks(:,1) == 1)
                    %disp('Building(s) at lower x domain edge')
                end
                if any(blocks(:,2) == imax)
                    %disp('Building(s) at upper x domain edge')
                end
                if any(blocks(:,3) == 1)
                    %disp('Building(s) at lower y domain edge')
                end
                if any(blocks(:,4) == jtot)
                    %disp('Building(s) at upper y domain edge')
                end
            %end
            
            
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
                for j = 1:nzw
                    if ((j == nzw) && (remz > 0))
                        lh = height - remz + 1;
                        uh = height;
                    else
                        lh = (j - 1) * nzw;
                        uh = j * nzw - 1;
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
                for j = 1:nzw
                    if ((j == nzw) && (remz > 0))
                        lh = height - remz + 1;
                        uh = height;
                    else
                        lh = (j - 1) * nzw;
                        uh = j * nzw - 1;
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
            da_pp.addvar(obj, 'boundingwallfacets', boundingwallfacets);
            obj.nboundingwallfacets = nboundingwallfacets;
        end
        
        function createfloors(obj)
            imax = obj.imax;
            jtot = obj.jtot;
            maxsize = obj.maxsize;
            % Create floors
            % fill space between blocks with floors (streets etc.)
            % floors only have x and y coordinates, with building ID = -101
            
            % obj.floorfacets holds floor facets information.
            % Format: orientation, wall type, block id, building id,
            % isinternal, il, iu, jl, ju, kl, ku
            
            %try
            blocks = obj.blocks;
            %catch
            %blocks = [];
            %end
            
            
            if obj.lflat
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
                    rem_x = rem(imax, maxsize)
                    xl = imax - rem_x + 1
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
                    %nblockstotal = 1;
                    blocks = [1, imax, 1, jtot, 0, 0, 1, 0, 0, 0, 0];
                end

            else
                nblocks = obj.nblocks;
                M = ones(imax, jtot);
                BI = zeros(imax, jtot); %block index mask
                corm = zeros(imax, jtot); %mask with all wall-floor corners
                cornm = zeros(imax, jtot);
                
                for i = 1:nblocks
                    xl = blocks(i,1);
                    xu = blocks(i,2);
                    yl = blocks(i,3);
                    yu = blocks(i,4);
                    M(xl:xu, yl:yu) = 0;
                    BI(xl:xu, yl:yu) = i;
                end
                NM = 1 - M;
                
                %             if ltestplot
                %                 % figure
                %                 % imagesc(xb,yb,M')
                %                 % axis equal tight
                %                 % set(gca,'YDir','normal')
                %                 % title('building mask')
                %                 figure
                %                 imagesc(xb+1,yb-1,BI')
                %                 axis equal tight
                %                 set(gca,'YDir','normal')
                %                 title('building index mask')
                %             end
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
                for i = 1:nblocks
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
                                elseif BI(xl - 1, yu + 1) > 0 %corner with a south wall
                                    cornm(xl - 1, yu) = 10;
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
                                elseif BI(xu+1, yu + 1) > 0  %corner with a south wall
                                    cornm(xu+1, yu) = 15;
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
                                elseif BI(xu + 1, yu + 1) > 0 %corner with a west wall
                                    cornm(xu, yu + 1) = 8;
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
                                elseif BI(xu + 1, yl - 1) > 0 %corner with a west wall
                                    cornm(xu, yl - 1) = 10;
                                end
                            end
                        end
                    end
                end
                
                corm(iM > 0) = 1;
                
                
                %             if ltestplot
                %                 % figure
                %                 % imagesc(xb,yb,M2')
                %                 % axis equal tight
                %                 % set(gca,'YDir','normal')
                %                 % title('floor mask after adding floors around buildings')
                %
                %                 figure
                %                 imagesc(xb,yb,iM')
                %                 axis equal tight
                %                 set(gca,'YDir','normal')
                %                 title('floor indeces after adding floors around buildings')
                %
                %                 figure
                %                 blub=corm'+NM'.*0.5;
                %                 imagesc(xb,yb,blub)
                %                 axis equal tight
                %                 set(gca,'YDir','normal')
                %                 title('mask of all corners floor-wall')
                %
                %                 figure
                %                 imagesc(xb,yb,obj.cornm')
                %                 axis equal tight
                %                 set(gca,'YDir','normal')
                %                 title('mask of all corners floor-wall-wall')
                %             end
                
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
                                if isempty(last)
                                    last = min(jtot,first + maxsize - 1);
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
                % if ltestplot
                % figure
                % imagesc(xb,yb,M2')
                % axis equal tight
                % set(gca,'YDir','normal')
                % title('floormask after filling')
                % figure
                % imagesc(xb,yb,iM')
                % axis equal tight
                % set(gca,'YDir','normal')
                % title('floor indeces')
                % end
                
                lnan2 = find(isnan(floorfacets(:,6)));
                if ~isempty(lnan2)
                    floorfacets(lnan2(1):lnan2(end), :) = [];
                end
                
                % figure
                % subplot(1,2,1)
                % imagesc(M')
                % set(gca,'YDir','normal')
                % axis equal tight
                % title('before')
                % subplot(1,2,2)
                % imagesc(M2')
                % set(gca,'YDir','normal')
                % axis equal tight
                % title('after')
                
                
                % New format
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
                
                da_pp.addvar(obj, 'cornm', cornm);
                nfloorfacets = size(floorfacets, 1);
                
                for i = 1:nfloorfacets
                    floorfacets(i, 1:5) = [1, -1, i, -1, 0]; % SO: block ID should be nblocks + i
                    floorfacets(i, 10:11) = [0, 0];
                end
                

                %disp([num2str(obj.nfcts) ' facets, of which: ' num2str(obj.nblockfcts) ' from buildings, ' num2str(obj.nboundingwallfacets) ' from walls, ' num2str(obj.nfloorfacets) ' from floors.'])
                
                nblocks = obj.nblocks;
                %if nblocks > 0
                    blocks(:, 5:6) = blocks(:, 5:6) + 1;
                %end
                
                blocks(nblocks + 1 : nblocks + nfloorfacets, :) = zeros(nfloorfacets, 11);
                blocks(nblocks + 1 : nblocks + nfloorfacets, 1:6) = floorfacets(:, 6:11);
                blocks(:, 7:11) = zeros(nblocks + nfloorfacets, 5);
                
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
            floorfacets
            blocks
            da_pp.addvar(obj, 'floorfacets', floorfacets)
            da_pp.addvar(obj, 'nfloorfacets', nfloorfacets)
            da_pp.addvar(obj, 'nblockstotal', obj.nblocks + nfloorfacets);
        end
        
        function vsolc(obj)
            xh = obj.xh;
            yh = obj.yh;
            zh = obj.zh;
            blocks = obj.blocks; nblocks = obj.nblocks
            facets = obj.facets; nfcts  = obj.nfcts;
            floorfacets = obj.floorfacets; nfloorfacets = obj.nfloorfacets;
            boundingwallfacets = obj.boundingwallfacets; nboundingwallfacets = obj.nboundingwallfacets; nbw = obj.nboundingwallfacets;
            buildings = obj.buildings;
            centerweight = obj.centerweight; cornerweight = obj.cornerweight;
            
            buildings(:, 7) = 0;
            for i = 1:nblocks
                index = find(buildings(:, 1) <= buildings(i, 1) & buildings(:, 2) >= buildings(i, 2) & buildings(:, 3) <= buildings(i, 3) & buildings(:, 4) >= buildings(i, 4));
                if buildings(index, 7) ~=0 && buildings(index, 7) ~= facets(blocks(i, 7), 4)
                    disp('sth went probably wrong, neighbouring blocks appear to belong to different buildings')
                end
                buildings(index, 7) = facets(blocks(i, 7), 4);
            end
            obj.buildings = buildings;
            
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

            Z = obj.Z; solaz = obj.solaz;
            %% vector to sun
            x = sin(Z / 360 * 2 * pi) * cos((solaz - 90) / 360 * 2 * pi);  %-90 since it's from north, not from east (= our x coordinate)
            y = sin(Z / 360 * 2 * pi) * sin((solaz + 90) / 360 * 2 * pi);
            z = cos(Z / 360 * 2 * pi);
            v1 = [x, y, z];
            
            % facets is the first 4 columns of obj.facets so can just replace it
            % with obj.facets where it appears.
            
            
            
            %% create blocks to test for intersection
            % coordinates in physical space (not indeces)
            % blocks is blocks so can just replace it with obj.blocks 
            %da_pp.addvar(obj, 'blocks_phys', zeros(obj.nblocks + obj.nboundingwallfacets, 6));
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
            
            
            % boundingwalls is boundingwalls, and boundingwalls(:, 1) =
            % obj.boundingwallfacets(:, 6) - add 5 to the index
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
            
            da_pp.addvar(obj, 'blocks_phys', blocks_phys);
            boundingwalls = boundingwallfacets(:, 6:end);
            
            % floors is floors, so floors(:, 1) = obj.floorfacets(:, 6)  - add 5
            % to the index
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
            bl = blocks_phys;

            
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
                [ndim, ~, co] = da_pp.detsub(i,facets,blocks,floors,boundingwalls,cornm,xh,yh,zh,delta);
                %[ndim, ~, co] = da_pp.detsub(obj, i); %determine corners
                %ndim=(dim1*dim2), number of cells of that facet
                %co = returns xyz-coordinates of center and 4 corners clockwise from
                %bottom left slightly shifter and
                %returns xyz-coordinates of center and 4 corners clockwise from
                %bottom left
                
                %count how much area cannot see sun (i.e. view is blocked)
                [ as ] = da_pp.prblckd(i,-999,co,ndim,true,v1,-999,-999,facets,centerweight,cornerweight,nblocks,nboundingwallfacets,blocks_phys);
                %[ as ] = da_pp.prblckd(obj, i, -999, co, ndim, true, v1, -999, -999);
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
            %da_pp.addvar(obj, 'Sdir', zeros(obj.nfcts, 1)); %direct solar radiation on every facet            
            for i = 1:nfcts  %direct from sky
                if facets(i, 1) == 1 %horizontal
                    phi = 0;
                else
                    phi = 90; %vertical
                end
                Sdir(i) = I * cos((Z - phi) / 360 * 2 * pi) * cos(walltheta(i) / 360 * 2 * pi) * asl(i);
            end
            
            da_pp.addvar(obj, 'asl', asl);
            da_pp.addvar(obj, 'Sdir', Sdir);

%             if ltestplot
%                 % facets is the first 4 columns of obj.facets
%                 % floors is floors, so floors(:, 1) = obj.floorfacets(:, 6)  - add 5
%                 % blocks is blocks so can just replace it with obj.blocks 
%                 % boundingwalls is boundingwalls, so boundingwalls(:, 1) = obj.boundingwallfacets(:, 6) - add 5 
% %                 whichblock = 17;
% %                 wbl = find(obj.facets(:,3) == whichblock & obj.facets(:,4) > 0);
% %                 wfl = obj.facets(wbl, :);
%                 scale=2;
%                 scalef=1.5;
%                 h=figure;
%                 set(gcf,'units','centimeters','position',[2 2 2+14.5*scale 2+14.5*scale]);
%                 set(h,'PaperPosition',[2 2 2+14.5*scale 2+14.5*scale]);
%                 set(h,'PaperUnits','centimeters');
%                 set(h,'renderer','painters');
%                 
%                 dx = ones(obj.imax, 1) * obj.dx;
%                 dy = ones(obj.jtot, 1) * obj.dy;
%                 dz = obj.dzf;
%                 
%                 for i=1:obj.nfcts
%                     bi = obj.facets(i,3); %block index
%                     fi = obj.facets(i,1); %orientation
%                     
%                     if (obj.facets(i,4) == -1) %it is a floor, not a building
%                         il = obj.floorfacets(bi, 6);
%                         iu = obj.floorfacets(bi, 7);
%                         jl = obj.floorfacets(bi, 8);
%                         ju = obj.floorfacets(bi, 9);
%                         xl = xh(il);
%                         xu = xh(iu+1);
%                         yl = yh(jl);
%                         yu = yh(ju+1);
%                         zl = obj.zh(0+1+1); %no 0th index
%                         zu = obj.zh(0+1+1);
%                     elseif (obj.facets(i,4) == -101)  %it is a bounding wall
%                         il = obj.boundingwallfacets(bi, 6);
%                         iu = obj.boundingwallfacets(bi, 7);
%                         jl = obj.boundingwallfacets(bi, 8);
%                         ju = obj.boundingwallfacets(bi, 9);
%                         kl = obj.boundingwallfacets(bi, 10) + 1;
%                         ku = obj.boundingwallfacets(bi, 11) + 1 + 1;
%                         
%                         if (fi == 2)
%                             xl = obj.xf(iu) + dx(iu) / 2;
%                             xu = obj.xf(iu) + dx(iu);
%                             yl = obj.yf(jl) - dy(jl) / 2;
%                             yu = obj.yf(ju) + dy(ju) / 2;
%                             zl = obj.zf(kl) - dz(kl) / 2;
%                             zu = obj.zf(ku) + dz(ku) / 2;
%                         elseif (fi == 3)
%                             xl = obj.xf(il) - dx(il);
%                             xu = obj.xf(il) - dx(il) / 2;
%                             yl = obj.yf(jl) - dy(jl) / 2;
%                             yu = obj.yf(ju) + dy(ju) / 2;
%                             zl = obj.zf(kl) - dz(kl) / 2;
%                             zu = obj.zf(ku) + dz(ku) / 2;
%                         elseif (fi == 4)
%                             xl = obj.xf(il) - dx(il) / 2;
%                             xu = obj.xf(iu) + dx(iu) / 2;
%                             yl = obj.yf(jl) - dy(jl);
%                             yu = obj.yf(jl) - dy(jl) / 2;
%                             zl = obj.zf(kl) - dz(kl) / 2;
%                             zu = obj.zf(ku) + dz(ku) / 2;
%                         else %if (fi==5)
%                             xl = obj.xf(il) - dx(il) / 2;
%                             xu = obj.xf(iu) + dx(iu) / 2;
%                             yl = obj.yf(ju) + dy(ju) / 2;
%                             yu = obj.yf(ju) + dy(ju);
%                             zl = obj.zf(kl) - dz(kl) / 2;
%                             zu = obj.zf(ku) + dz(ku) / 2;
%                         end
%                     else
%                         il = obj.blocks(bi, 1);
%                         iu = obj.blocks(bi, 2);
%                         jl = obj.blocks(bi, 3);
%                         ju = obj.blocks(bi, 4);
%                         kl = obj.blocks(bi, 5) + 1; %no 0th index
%                         ku = obj.blocks(bi, 6) + 1; %no 0th index
%                         
%                         xl = xh(il);
%                         xu = xh(iu + 1);
%                         yl = obj.yf(jl);
%                         yu = obj.yf(ju + 1);
%                         zl = xh(kl) - 0.25; %move a bit down for better plots
%                         zu = xh(ku + 1);
%                     end
%                     
%                     switch obj.facets(i, 1)
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
%                     
%                     %shaded vs non-shaded
%                     cs = obj.asl(i) * [0.75 0.75 0.75] + [0.25 0.25 0.25];
%                     
%                     patch(x,y,z, cs,'FaceLighting','none');
%                     hold on
%                 end
%                 
%                 for i=1:size(wfl,1)
%                     xl=xb(wfl(i,6));
%                     xu=xb(wfl(i,7));
%                     yl=yb(wfl(i,8));
%                     yu=yb(wfl(i,9));
%                     zl=zb(wfl(i,10)+1);
%                     zu=zb(wfl(i,11)+1);
%                     
%                     if (wsl(wbl(i))) %if not self shaded
%                         switch wfl(i,1)
%                             case 1
%                                 quiver3(xl,yl,zu,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3(xl,yu,zu,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3(xu,yu,zu,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3(xu,yl,zu,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3((xu+xl)/2,(yl+yu)/2,zu,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                             case 2
%                                 quiver3(xl,yl,zu,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3(xl,yu,zu,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3(xl,yu,zl,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3(xl,yl,zl,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3(xl,(yl+yu)/2,(zl+zu)/2,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                             case 3
%                                 quiver3(xu,yl,zu,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3(xu,yu,zu,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3(xu,yu,zl,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3(xu,yl,zl,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3(xu,(yl+yu)/2,(zl+zu)/2,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                             case 4
%                                 quiver3(xl,yu,zu,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3(xl,yu,zl,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3(xu,yu,zl,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3(xu,yu,zu,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3((xu+xl)/2,yu,(zu+zl)/2,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                             case 5
%                                 quiver3(xl,yl,zu,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3(xl,yl,zl,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3(xu,yu,zl,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3(xu,yl,zu,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                                 quiver3((xu+xl)/2,yl,(zu+zl)/2,v1(1),v1(2),v1(3),15.123,'r','LineWidth',2,'MaxHeadSize',1)
%                         end
%                     end
%                     
%                 end
%                 %asl(wbl)
%                 
%                 
%                 view(v1)
%                 axis equal
%                 xlim([-1 xb(end)+1])
%                 ylim([-1 yb(end)+1])
%                 %IC25
%                 %zlim([0 72])
%                 zlim([0 max(zmax5+1)+5])
%                 set(gca,'TickLabelInterpreter','latex')
%                 set(gca,'FontSize',12*scalef)
%                 xlabel('x [m]','Interpreter','latex','FontSize',12*scalef)
%                 ylabel('y [m]','Interpreter','latex','FontSize',12*scalef)
%                 zlabel('z [m]','Interpreter','latex','FontSize',12*scalef)
%                 set(gcf, 'Color', 'c');
%                 % export_fig shading.eps
%             end
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
                if (obj.facets(i,4) < 0 && obj.facets(i,4) > -100) %it is a floor, not a building
                    il = obj.floorfacets(bi, 6);
                    iu = obj.floorfacets(bi, 7);
                    jl = obj.floorfacets(bi, 8);
                    ju = obj.floorfacets(bi, 9);
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
                    il = obj.floorfacets(bi, 6);
                    iu = obj.floorfacets(bi, 7);
                    jl = obj.floorfacets(bi, 8);
                    ju = obj.floorfacets(bi, 9);
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
        
%         function [ndim, area, co] = detsub(obj, i)
%             % returns ndim=(dim1*dim2), area, center & corner of every facet (could be preprocessed and
%             % saved)
%             
%             bi = obj.facets(i, 3); %block index
%             fi = obj.facets(i, 1); %facet index
%             ci = obj.facets(i, 4); %building index (-1 for roads, -101 for bounding wall)
%             
%             il=-1; iu=-1; jl=-1; ju=-1; kl=-1; ku=-1; it=-1; jt=-1; kt=-1;
%             co = -1;
%             ndim = -1;
%             area = -1;
%                                  
%             if (ci <= -100) %it is a bounding wall
%                 kl = obj.boundingwallfacets(bi, 10) + 1 + 1;
%                 ku = obj.boundingwallfacets(bi, 11) + 1 + 1;
%                 if (fi == 2)  %east, facing west
%                     jl = obj.boundingwallfacets(bi, 8);     %lower y index of floor facet 1
%                     ju = obj.boundingwallfacets(bi, 9);     %upper y index of floor facet 1
%                     it = obj.boundingwallfacets(bi, 6);                    
%                     x = obj.xh(it + 1);                    
%                     ndim = ((ju - jl) + 1) * ((ku - kl) + 1);                    
%                 elseif (fi == 3)  %west, facing east
%                     jl = obj.boundingwallfacets(bi, 8);     %lower y index of floor facet 1
%                     ju = obj.boundingwallfacets(bi, 9);     %upper y index of floor facet 1
%                     it = obj.boundingwallfacets(bi, 7);                    
%                     x = obj.xh(it);
%                     ndim = ((ju - jl) + 1) * ((ku - kl) + 1);
%                 elseif (fi == 4)  %south, facing north
%                     il = obj.boundingwallfacets(bi, 6);
%                     iu = obj.boundingwallfacets(bi, 7);
%                     jt = obj.boundingwallfacets(bi, 9);                   
%                     y = obj.yh(jt);
%                     ndim = ((iu - il) + 1) * ((ku - kl) + 1);                    
%                 else %if (fi==5)  %north, facing south
%                     jt = obj.boundingwallfacets(bi, 8);
%                     il = obj.boundingwallfacets(bi, 6);
%                     iu = obj.boundingwallfacets(bi, 7);                    
%                     y = obj.yh(jt + 1);
%                     ndim=((iu - il) + 1) * ((ku - kl) + 1);
%                 end
%             elseif (ci >= 0) %deal with floors separately
%                 if (fi == 1)  %top
%                     il = obj.blocks(bi, 1);     %lower x index of facet 1
%                     iu = obj.blocks(bi, 2);     %upper x index of facet 1
%                     jl = obj.blocks(bi, 3);     %lower y index of facet 1
%                     ju = obj.blocks(bi, 4);     %upper y index of facet 1
%                     kt = obj.blocks(bi, 6) + 1;
%                     z = obj.zh(kt + 1);                    
%                     ndim=((iu-il)+1)*((ju-jl)+1);                   
%                 elseif ((fi == 2) || (fi == 3))  %west / east
%                     jl = obj.blocks(bi, 3);     %lower y index of facet 1
%                     ju = obj.blocks(bi, 4);     %upper y index of facet 1
%                     kl = obj.blocks(bi, 5) + 1;     %lower z index of facet 1
%                     ku = obj.blocks(bi, 6) + 1;     %upper z index of facet 1
%                     
%                     if fi == 2  %west, facing west
%                         it = obj.blocks(bi, 1);
%                         x = obj.xh(it);
%                     else %east, facing east
%                         it = obj.blocks(bi, 2);
%                         x = obj.xh(it + 1);
%                     end
%                     
%                     ndim = ((ju - jl) + 1) *( (ku - kl) + 1);
%                                        
%                 elseif ((fi == 4) || (fi == 5)) %north / south
%                     il = obj.blocks(bi, 1) ;    %lower x index of facet 1
%                     iu = obj.blocks(bi, 2) ;    %upper x index of facet 1
%                     kl = obj.blocks(bi, 5) + 1 ;    %lower y index of facet 1
%                     ku = obj.blocks(bi, 6) + 1 ;    %upper y index of facet 1
%                     
%                     if fi == 4 %north, facing north
%                         jt = obj.blocks(bi, 4);
%                         y = obj.yh(jt + 1);
%                     else  %south, facing south
%                         jt = obj.blocks(bi, 3);
%                         y = obj.yh(jt);
%                     end
%                     
%                     ndim=((iu - il) + 1) * ((ku - kl) + 1);
%                 end
%             end
%             
%             %return xyz-coordinates of center and 4 corners clockwise from bottom left
%             %move coordinates of corners very slightly to the interior of the facet (by
%             %delta)
%             delta = 0.01;
%             if (ci >= 0 || ci <= -100) %it's not a floor -> it's a building or bounding wall              
%                 if (fi == 1)
%                     co = [(obj.xh(iu + 1) + obj.xh(il)) * 0.5, (obj.yh(ju + 1) + obj.yh(jl)) * 0.5, z; ... %center
%                            obj.xh(il) + delta, obj.yh(jl) + delta, z; ... %4 corners
%                            obj.xh(il) + delta, obj.yh(ju + 1) - delta, z; ...
%                            obj.xh(iu + 1) - delta, obj.yh(ju + 1) - delta, z; ...
%                            obj.xh(iu + 1) - delta, obj.yh(jl) + delta, z; ...
%                            obj.xh(il), obj.yh(jl), z; ... %4 corners
%                            obj.xh(il), obj.yh(ju + 1), z; ...
%                            obj.xh(iu + 1), obj.yh(ju + 1), z; ...
%                            obj.xh(iu + 1), obj.yh(jl), z; ...
%                            ];
%                     area = (obj.xh(iu + 1) - obj.xh(il)) * (obj.yh(ju + 1) - obj.yh(jl));
%                     
%                 elseif ((fi == 2) || (fi == 3))                    
%                     co = [x, (obj.yh(ju + 1) + obj.yh(jl))  * 0.5, (obj.zh(ku + 1) + obj.zh(kl)) * 0.5; ... %center
%                           x, obj.yh(jl) + delta, obj.zh(kl) + delta; ... %4 corners
%                           x, obj.yh(jl) + delta, obj.zh(ku + 1) - delta; ...
%                           x, obj.yh(ju + 1) - delta, obj.zh(ku + 1) - delta; ...
%                           x, obj.yh(ju + 1) - delta, obj.zh(kl) + delta; ...
%                           x, obj.yh(jl), obj.zh(kl); ... %4 corners
%                           x, obj.yh(jl), obj.zh(ku + 1); ...
%                           x, obj.yh(ju + 1), obj.zh(ku + 1); ...
%                           x, obj.yh(ju + 1), obj.zh(kl); ...
%                           ];
%                     area = (obj.zh(ku + 1) - obj.zh(kl)) * (obj.yh(ju + 1) - obj.yh(jl));
%                     
%                 elseif ((fi == 4) || (fi == 5))                
%                     co = [(obj.xh(iu + 1) + obj.xh(il)) * 0.5, y, (obj.zh(ku + 1) + obj.zh(kl)) * 0.5; ... %center
%                           obj.xh(il) + delta, y, obj.zh(kl) + delta; ... %4 corners
%                           obj.xh(il) + delta, y, obj.zh(ku + 1) - delta; ...
%                           obj.xh(iu + 1) - delta, y, obj.zh(ku + 1) - delta; ...
%                           obj.xh(iu + 1) - delta, y, obj.zh(kl) + delta; ...
%                           obj.xh(il), y, obj.zh(kl); ... %4 corners
%                           obj.xh(il), y, obj.zh(ku + 1); ...
%                           obj.xh(iu + 1), y, obj.zh(ku + 1); ...
%                           obj.xh(iu + 1), y, obj.zh(kl); ...
%                           ];
%                     area = (obj.xh(iu + 1) - obj.xh(il)) * (obj.zh(ku + 1) - obj.zh(kl));
%                 end
%                 
%             else  % it is a floor, not a building
%                 il = obj.floorfacets(bi, 6);     %lower x index of floor facet 1
%                 iu = obj.floorfacets(bi ,7);     %upper x index of floor facet 1
%                 jl = obj.floorfacets(bi, 8);     %lower y index of floor facet 1
%                 ju = obj.floorfacets(bi, 9);     %upper y index of floor facet 1
%                 z = 1;
%                 ndim = ((iu - il) + 1) * ((ju - jl) + 1);
%                 co = [(obj.xh(iu + 1) + obj.xh(il)) * 0.5, (obj.yh(ju + 1) + obj.yh(jl)) * 0.5, z; ... %center
%                       obj.xh(il) + delta, obj.yh(jl) + delta, z; ... %4 corners
%                       obj.xh(il) + delta, obj.yh(ju + 1) - delta, z; ...
%                       obj.xh(iu + 1) - delta, obj.yh(ju + 1) - delta, z; ...
%                       obj.xh(iu + 1) - delta, obj.yh(jl) + delta, z; ...
%                       obj.xh(il), obj.yh(jl), z; ... %4 corners
%                       obj.xh(il), obj.yh(ju + 1), z; ...
%                       obj.xh(iu + 1), obj.yh(ju+1), z; ...
%                       obj.xh(iu + 1), obj.yh(jl), z; ...
%                       ];
%                 area = (obj.xh(iu + 1) - obj.xh(il)) * (obj.yh(ju + 1) - obj.yh(jl));
%                 %if it is 1*1 nothing has to be done since it matches both walls
%                 if ju - jl > 0  %floor facet along a west or east wall
%                     %test if it is a wall-wall-floor corner, if so, the floor facet needs
%                     %to be cut diagonally, see "createfloors.m"
%                     if obj.cornm(iu, ju) == 10  %south/west corner
%                         %JUST REMOVE IT ONE dx OR dy FROM THE OTHER BLOCK
%                         co(2+1, 2) = obj.yh(ju) + delta;
%                         co(6+1, 2) = obj.yh(ju);
%                         area = area - (obj.xh(2) - obj.xh(1)) * (obj.yh(2) - obj.yh(1)) / 2;
%                     elseif obj.cornm(iu, ju) == 15  %south/east corner
%                         co(3+1, 2) = obj.yh(ju) + delta;
%                         co(7+1, 2) = obj.yh(ju);
%                         area = area - (obj.xh(2) - obj.xh(1)) * (obj.yh(2) - obj.yh(1)) / 2;
%                     end
%                     if obj.cornm(iu, jl) == 8  %north/west corner
%                         co(1+1, 2) = obj.yh(jl + 1) - delta;
%                         co(5+1, 2) = obj.yh(jl + 1);
%                         area = area - (obj.xh(2) - obj.xh(1)) * (obj.yh(2) - obj.yh(1)) / 2;
%                     elseif obj.cornm(iu, jl) == 12  %north/east corner
%                         co(4+1,2) = obj.yh(jl + 1) - delta;
%                         co(8+1,2) = obj.yh(jl + 1);
%                         area = area - (obj.xh(2) - obj.xh(1)) * (obj.yh(2) - obj.yh(1)) / 2;
%                     end
%                     
%                 elseif iu - il > 0 %floor facet along along a south or north wall
%                     if obj.cornm(iu, ju) == 10  %south/west corner
%                         co(4+1, 1) = obj.xh(iu) + delta;
%                         co(8+1, 1) = obj.xh(iu);
%                         area = area - (obj.xh(2) - obj.xh(1)) * (obj.yh(2) - obj.yh(1)) / 2;
%                     elseif obj.cornm(iu,ju)==8  %north/west corner
%                         co(3+1,1) = obj.xh(iu) + delta;
%                         co(7+1,1) = obj.xh(iu);
%                         area = area - (obj.xh(2) - obj.xh(1)) * (obj.yh(2) - obj.yh(1)) / 2;
%                     end
%                     if obj.cornm(il, ju) == 12  %north/east corner
%                         co(2+1, 1) = obj.xh(il + 1) - delta ;
%                         co(6+1, 1) = obj.xh(il + 1) ;
%                         area = area - (obj.xh(2) - obj.xh(1)) * (obj.yh(2) - obj.yh(1)) / 2;
%                     elseif obj.cornm(il, ju) == 15  %south/east corner
%                         co(1+1, 1) = obj.xh(il + 1) - delta;
%                         co(5+1, 1) = obj.xh(il + 1);
%                         area = area - (obj.xh(2) - obj.xh(1)) * (obj.yh(2) - obj.yh(1)) / 2;
%                     end
%                 end
%             end
%         end
        
%         function [flag ,tmin] = rbi(origin, direction, vec)
%             %  Ray/box intersection using the Smits' algorithm
%             %
%             % Input:
%             %    origin.
%             %    direction.
%             %    vec = [xl,xu,yl,yu,zl,zu]  (=Block)
%             % Output:
%             %    flag: (0) Reject, (1) Intersect.
%             %    tmin: distance from the ray origin.
%             
%             vmin=[vec(1) vec(3) vec(5)];
%             vmax=[vec(2) vec(4) vec(6)];
%             
%             
%             if ((direction(1) == 0) && (((vmin(1) - origin(1))==0) || ((vmax(1) - origin(1))==0))) %comparing floats!  %0/0 is undefined whereas 1/0 = inf, -1/0=-inf works
%                 %this will only be the case if the ray is tangential -> no intersection
%                 %(or infinite intersections...)
%                 tmin=0;
%                 tmax=inf;
%             elseif (direction(1) >= 0)
%                 tmin = (vmin(1) - origin(1)) / direction(1);
%                 tmax = (vmax(1) - origin(1)) / direction(1);
%             else
%                 tmin = (vmax(1) - origin(1)) / direction(1);
%                 tmax = (vmin(1) - origin(1)) / direction(1);
%             end
%             
%             if ((direction(2) == 0) && (((vmin(2) - origin(2))==0) || ((vmax(2) - origin(2))==0))) %comparing floats!
%                 tymin=0;
%                 tymax=inf;
%             elseif (direction(2) >= 0)
%                 tymin = (vmin(2) - origin(2)) / direction(2);
%                 tymax = (vmax(2) - origin(2)) / direction(2);
%             else
%                 tymin = (vmax(2) - origin(2)) / direction(2);
%                 tymax = (vmin(2) - origin(2)) / direction(2);
%             end
%             
%             if ( (tmin > tymax) || (tymin > tmax) )
%                 flag = 0;
%                 tmin = NaN;  %originally -1
%                 return;
%             end
%             
%             if (tymin > tmin)
%                 tmin = tymin;
%             end
%             
%             if (tymax < tmax)
%                 tmax = tymax;
%             end
%             
%             if ((direction(3) == 0) && (((vmin(3) - origin(3))==0) || ((vmax(3) - origin(3))==0))) %comparing floats!
%                 tzmin=0;
%                 tzmax=inf;
%             elseif (direction(3) >= 0)
%                 tzmin = (vmin(3) - origin(3)) / direction(3);
%                 tzmax = (vmax(3) - origin(3)) / direction(3);
%             else
%                 tzmin = (vmax(3) - origin(3)) / direction(3);
%                 tzmax = (vmin(3) - origin(3)) / direction(3);
%             end
%             
%             
%             if ((tmin > tzmax) || (tzmin > tmax))
%                 flag = 0;
%                 tmin = -1;
%                 return;
%             end
%             
%             if (tzmin > tmin)
%                 tmin = tzmin;
%             end
%             
%             if (tzmax < tmax)
%                 tmax = tzmax;
%             end
%             
%             % if( (tmin < t1) && (tmax > t0) )
%             flag = 1;
%             % else
%             %    flag = 0;
%             %    tmin = -1;
%             % end;
%         end
        
%         function [ prcntgblckd ] = prblckd(obj, i, j, coa, ndima, sun, v1, cob, ndimb)
%             %prcntgblckd     = percentage of view blocked
%             %i               = index of facet 1
%             %j               = index of facet 2      (only if sun==false)
%             %coa             = corners of facet 1
%             %ndima           = dimension of facet 1
%             %sun             = calculation if facet can see sun, not if two facets can
%             %                  see each other
%             %v1              = vector to the sun
%             %cob             = corners of facet 2    (only if sun==false)
%             %ndimb           = dimension of facet 2  (only if sun==false)
%             %facets               = facets
%             %centerweight    = % of viewfield if centers see each other
%             %cornerweight    = % per corner if they see each other
%             %nblocks         = number of blocks
%             %nbw             = number of bounding walls
%             %bl              = list of blocks
%             
%             testcrit = 0; %=4;
%             prcntgblckd = 0;
%             flag = 0;
%         
%             if sun  %calculation between facet 1 and the sun
%                 disp('running prblckd')
%                 bi = obj.facets(i, 3); %block index
%                 ci = obj.facets(i, 4); %building index (-1 for roads, -101 for bounding wall)
%                 if ndima > testcrit  %also test corner, otherwise only test center                    
%                     for k = 1:5
%                         facetpoint = coa(k,:);
%                         for n = 1:obj.nblocks %check if any block intersects
%                             if ((n == bi) && (ci > 0))  %don't intersect with itself
%                                 continue
%                             end
%                             
%                             [flag, dint] = da_pp.rbi(facetpoint, v1, obj.blocks_phys(n,:));
%                             
%                             intersection = facetpoint + dint * v1;
%                             if intersection(3) < facetpoint(3) %downstream direction of sun, thus not blocking the sun
%                                 flag = 0;
%                             end
%                             if flag == 1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
%                                 if k == 1
%                                     prcntgblckd = obj.centerweight;
%                                 else
%                                     prcntgblckd = prcntgblckd + obj.cornerweight;
%                                 end
%                                 break %get out of inner for loop
%                             end
%                         end
%                         
%                         if ~flag
%                             for m = 1:obj.nboundingwallfacets %check if any bounding wall intersects                     
%                                 [flag, dint] = da_pp.rbi(facetpoint, v1, obj.blocks_phys(m + obj.nblocks, :));
%                                 
%                                 intersection = facetpoint + dint * v1;
%                                 
%                                 if intersection(3) < facetpoint(3) %downstream direction of sun
%                                     flag = 0;
%                                 end
%                                 if flag == 1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
%                                     if k == 1 %it's the center of the facet
%                                         prcntgblckd = obj.centerweight;
%                                     else
%                                         prcntgblckd = prcntgblckd + obj.cornerweight;
%                                     end
%                                     break %get out of inner for loop
%                                 end
%                             end
%                         end
%                     end
%                 else
%                     % %disp('small area, checking only center')
%                     % %disp(['facet: ' num2str(i)])
%                     facetpoint = coa(1, :);  %center of facet
%                     for n = 1:obj.nblocks %check if any block intersects
%                         if ((n == bi) && (ci > 0))  %don't intersect with itself
%                             continue
%                         end
%                         
%                         [flag, dint] = da_pp.rbi(facetpoint, v1, obj.blocks_phys(n,:));
%                         
%                         intersection = facetpoint + dint*v1;
%                         if intersection(3)<facetpoint(3) %downstream direction of sun, thus not blocking the sun
%                             % intersectionss(:)=0;
%                             flag = 0;
%                         end
%                         if flag == 1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
%                             prcntgblckd = 1;
%                             break %get out of inner for loop
%                         end
%                     end
%                     
%                     if ~flag
%                         for k = 1:obj.nboundingwallfacets %check if any bounding wall intersects
%                             
%                             [flag, dint] = da_pp.rbi(facetpoint, v1, obj.blocks_phys(k+nblocks,:));
%                             
%                             intersection = facetpoint + dint * v1;
%                             
%                             if intersection(3) < facetpoint(3) %downstream direction of sun
%                                 flag = 0;
%                             end
%                             if flag == 1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
%                                 prcntgblckd = 1;
%                                 break %get out of inner for loop
%                             end
%                         end
%                     end                   
%                 end
%             else %calculation between facet 1 and facet 2
%                 %disp('facet to facet')
%                 %disp(['% blocked so far ' num2str(prcntgblckd)])
%                 
%                 bi = obj.facets(i,3); %block index
%                 ci = obj.facets(i,4); %building index (-1 for roads, -99 for bounding wall)
%                 bi2 = obj.facets(j,3); %block index
%                 fi = obj.facets(i,1); %facet index
%                 fi2 = obj.facets(j,1); %facet index
%                 
%                 if (ndima + ndimb) > (2 * testcrit)  %also test corner, otherwise only test center
%                     %disp('test corners')
%                     %disp(['% blocked so far ' num2str(prcntgblckd)])
%                     %find orientation to know which corners to compare
%                     ordera = [1 2 3 4 5];
%                     cas = fi * 10 + fi2;
%                     %disp('case')
%                     switch cas
%                         case 12, orderb = [1 3 4 5 2];
%                         case 13, orderb = [1 2 5 4 3];
%                         case 15, orderb = [1 3 2 5 4];
%                         case 21, orderb = [1 5 2 3 4];
%                         case 24, orderb = [1 5 4 3 2];
%                         case 31, orderb = [1 2 5 4 3];
%                         case 35, orderb = [1 5 4 3 2];
%                         case 42, orderb = [1 5 4 3 2];
%                         case 51, orderb = [1 3 2 5 4];
%                         case 53, orderb = [1 5 4 3 2];
%                         otherwise %14, 41, 23, 32, 25, 52, 34, 43, 45, 54
%                             orderb = [1 2 3 4 5];
%                     end
%                     
%                     for k = 1:5
%                         test = 1;
%                         oa = ordera(k);
%                         ob = orderb(k);
%                         facetpoint = coa(oa, :);
%                         facetpoint2 = cob(ob, :);
%                         v = [facetpoint2(1) - facetpoint(1), facetpoint2(2) - facetpoint(2), facetpoint2(3) - facetpoint(3)];
%                         dint = norm(v);
%                         v1 = v / dint;
%                         if ((fi == 1) && (v(3) <= 0)) %horizontal facets can't see facets whos center is lower
%                             test = 0;
%                         elseif ((fi == 2) && (v(1) >= 0)) %west facets cant see anything to the east
%                             test = 0;
%                         elseif ((fi == 3) && (v(1) <= 0)) %east facets cant see anything to the west
%                             test = 0;
%                         elseif ((fi == 4) && (v(2) <= 0)) %north facets cant see anything to the south
%                             test = 0;
%                         elseif ((fi == 5) && (v(2) >= 0)) %south facets cant see anything to the north
%                             test = 0;
%                         elseif ((fi ~= 1) && (fi2 == 1) && (v(3) > 0)) %vertical walls cant see any horizontal facet that is higher
%                             test = 0;
%                             % the following cases should be blocked by a block anyway, but
%                             % this is faster
%                         elseif ((fi == 4) || (fi == 5)) && (fi2 == 3) && (v(1) > 0) %north/south can't see an east facet if it is east of it
%                             test = 0;
%                         elseif ((fi == 4) || (fi == 5)) && (fi2 == 2) && (v(1) < 0) %north/south can't see a west facet if it is west of it
%                             test = 0;
%                         elseif ((fi == 2) || (fi == 3)) && (fi2 == 4) && (v(2) > 0) %west/east can't see a north facet if it is north of it
%                             test = 0;
%                         elseif ((fi == 2) || (fi == 3)) && (fi2 == 5) && (v(2) < 0) %west/east can't see a south facet if it is south of it
%                             test = 0;
%                         end
%                         if test
%                             %                 disp(['can potentially see each other' num2str(k)])
%                             %                 disp(['% blocked so far ' num2str(prcntgblckd)])
%                             for n = 1:obj.nblocks %check if any block intersects
%                                 if ((n == bi) && (ci > 0))  %don't intersect with itself
%                                     continue
%                                 end
%                                 
%                                 [flag, dints] = da_pp.rbi(facetpoint, v1, obj.blocks_phys(n,:));
%                                 
%                                 if dints <= (0 + eps)  %intersection is downstream
%                                     flag = 0;
%                                 elseif dints >= (dint - eps) %block is behind target subfacet
%                                     flag = 0;
%                                 end
%                                 if flag == 1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
%                                     if k == 1
%                                         prcntgblckd = obj.centerweight;
%                                     else
%                                         %disp(['blocked by ' num2str(n)])
%                                         prcntgblckd = prcntgblckd + obj.cornerweight;
%                                     end
%                                     break %get out of inner for loop
%                                 end
%                             end
%                         else %facets can't possibly see each other
%                             %disp('facets cant possibly see each other')
%                             if k == 1
%                                 prcntgblckd = obj.centerweight;
%                             else
%                                 prcntgblckd = prcntgblckd + obj.cornerweight;
%                             end
%                         end
%                     end
%                     
%                 else %small area, checking only center
%                     % %disp('small area, checkig only center')
%                     test = 1;
%                     facetpoint = coa(1, :);
%                     facetpoint2 = cob(1, :);
%                     v = [facetpoint2(1) - facetpoint(1), facetpoint2(2) - facetpoint(2), facetpoint2(3) - facetpoint(3)];
%                     dint = norm(v);
%                     v1 = v / dint;
%                     
%                     %abs(v)
%                     %%disp(['v1: ' num2str(v1)])
%                     %test if facets can potentially see each other
%                     if ((fi == 1) && (v(3) <= 0)) %horizontal facets can't see facets whos center is lower
%                         test = 0;
%                     elseif ((fi == 2) && (v(1) >= 0)) %west facets cant see anything to the east
%                         test = 0;
%                     elseif ((fi == 3) && (v(1) <= 0)) %east facets cant see anything to the west
%                         test = 0;
%                     elseif ((fi == 4) && (v(2) <= 0)) %north facets cant see anything to the south
%                         test = 0;
%                     elseif ((fi == 5) && (v(2) >= 0)) %south facets cant see anything to the north
%                         test = 0;
%                     elseif ((fi ~= 1) && (fi2 == 1) && (v(3) > 0)) %vertical walls cant see any horizontal facet that is higher
%                         test = 0;
%                         % the following cases should be blocked by a block anyway, but
%                         % this is faster
%                     elseif ((fi == 4) || (fi == 5)) && (fi2 == 3) && (v(1) > 0) %north/south can't see an east facet if it is east of it
%                         test = 0;
%                     elseif ((fi == 4) || (fi == 5)) && (fi2 == 2) && (v(1) < 0) %north/south can't see a west facet if it is west of it
%                         test = 0;
%                     elseif ((fi == 2) || (fi == 3)) && (fi2 == 4) && (v(2) > 0) %west/east can't see a north facet if it is north of it
%                         test = 0;
%                     elseif ((fi == 2) || (fi == 3)) && (fi2 == 5) && (v(2) < 0) %west/east can't see a south facet if it is south of it
%                         test = 0;
%                     end
%                     
%                     %%disp(['test: ' num2str(test)])
%                     
%                     if test %facets could see each other, check if blocked
%                         %calculate actual intersection distance
%                         for n = 1:obj.nblocks %check if any block intersects
%                             if ((n == bi) && (ci > 0))  %don't intersect with itself
%                                 continue
%                             end
%                             
%                             [flag, dints] = da_pp.rbi(facetpoint, v1, obj.blocks_phys(n,:));
%                             
%                             if dints <= 0  %intersection is downstream
%                                 flag = 0;
%                             elseif dints >= dint %block is behind target subfacet
%                                 flag = 0;
%                             end
%                             if flag == 1 %found an intersection, i.e. facet is shaded, don't test against any more blocks
%                                 %  %disp(['intersection with block: ' num2str(n)])
%                                 prcntgblckd = 1;
%                                 break %get out of inner for loop
%                             end
%                         end
%                         
%                     else %facets can't possibly see each other
%                         prcntgblckd = 1;
%                     end
%                 end
%             end          
%         end
        
      

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


            %xb = obj.xh; yb = obj.yh; zb = obj.zh; 
            %nfcts = obj.nfcts;
            %% View factors
            % Calculate fraction of facets seeing each other (not blocked)
            % Calculate orientation, no overlap on orthogonal planes
            % calculate viewfactors
            vf = zeros(nfcts, nfcts);            
            pf1sf2 = zeros(nfcts, nfcts);
            %da_pp.addvar(obj, 'facetarea', zeros(obj.nfcts, 1));
            A = zeros(nfcts, 1);
            %tim = zeros(obj.nfcts - 1, 1);
            %%
            %do dummy parallelisation by starting multiple matlab and only iterate over
            %a fraction and then write to file and merge files
            for i = 1:(nfcts - 1)
                bi = facets(i, 3); %block index
                fi = facets(i, 1); %facet index
                ci = facets(i, 4); %building index (-1 for roads, -99 for bounding wall)
                %[ndima, areaa, coa] = da_pp.detsub(obj, i);
                
                [ ndima, areaa, coa] = da_pp.detsub(i,facets,blocks,floors,boundingwalls,cornm,xh,yh,zh,delta);
                for j = (i + 1):nfcts
                    bi2 = facets(j, 3); %block index
                    fi2 = facets(j, 1); %facet index
                    ci2 = facets(j, 4); %building index (-1 for roads, -99 for bounding wall)
                    if ((fi2 == fi) || ((ci2 * bi2) == (ci * bi)))  %facets looking in same direction or on same block&building (ci*bi is unique), bi alone could also be a floor or a bounding wall
                        continue
                    end
                    
                    %[ndimb, areab, cob] = da_pp.detsub(obj, j);
                    [ ndimb, areab, cob] = da_pp.detsub(j,facets,blocks,floors,boundingwalls,cornm,xh,yh,zh,delta);
                    %[prblckdij] = da_pp.prblckd(obj, i, j, coa, ndima, false, -999, cob, ndimb);
                    [prblckdij] = da_pp.prblckd(i,j,coa,ndima,false,-999,cob,ndimb,facets,centerweight,cornerweight,nblocks,nboundingwallfacets,blocks_phys);
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
                            [coaa, cobb, pf1sf2u] = da_pp.slice(fi, fi2, coa, cob, cornerweight, centerweight);
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
                        [F12, F21] = da_pp.ViewFactor(coaa, cobb, areaa, areab, glpo, vcorner);
                        %increase glpo if view factors are big
                        if F12 > 0.5 || F21 > 0.5
                            [F12, F21] = da_pp.ViewFactor(coaa, cobb, areaa, areab, glpo + 10, vcorner);
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
            da_pp.addvar(obj, 'vf', vf)
            da_pp.addvar(obj, 'svf', max(1 - sum(vf, 2), 0));
            da_pp.addvar(obj, 'facetarea', A);
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
        
%         function rayit(obj)
%             function [fct, wall] = loadfacets()
%                 %M = dlmread(['facets.inp.' num2str(expnr)],'',1,0);
%                 vars = {'or', 'wlid', 'blk', 'bld'};
%                 for n = 1:length(vars)
%                     fct.(vars{n}) = obj.facets(:, n);
%                 end
%                 
%                 %disp(obj.walltypes)
%                 %M = dlmread(['walltypes.inp.', obj.expnr],'',3,0);
%                 %disp(M)
%                 % wallid  lGR   z0 [m]     z0h [m]     al [-]   em [-]   d1 [m]  d2 [m]    cp1 [J/kgK]  cp2 [J/kgK]   l1 [W/(m K)]  l2 [W/(m K)]    k1 [W/(m K)]    k1 [W/(m K)]
%                 vars = {'id','lGR','z0','z0h','al', 'em', 'd1', 'd2', 'cp1', 'cp2', 'l1', 'l2', 'k1', 'k2'};
%                 for n = 1:length(vars)
%                     wall.(vars{n}) = obj.walltypes(:, n);
%                 end
%                 
%                 % check whether all referenced wall ids have been defined
%                 idref = unique(fct.wlid);
%                 iddef = unique(wall.id');
%                 
%                 if length(wall.id) > length(iddef)
%                     disp('ERROR: multiple definitions of a wall id')
%                     return
%                 end
%                 if length(idref) > length(iddef)
%                     disp('ERROR: more walltypes used than defined')
%                     return
%                 end
%                 
%                 
%                 % assign the wall type to each of the facets -- this way it is easy to know
%                 % in which structure the proporties can be found.
%                 fct.wltp = zeros(size(fct.or));
%                 
%                 for m = 1:size(fct.wlid, 1)
%                     if fct.wlid(m) >= 0 && fct.wlid(m) <= 10 && any(ismember(iddef, fct.wlid(m)))
%                         fct.wltp(n) = 1; %normal
%                     elseif fct.wlid(m) > 10 && any(ismember(iddef, fct.wlid(m)))
%                         fct.wltp(n) = 2; %Green roof
%                     elseif fct.wlid(m) < 0 && any(ismember(iddef, fct.wlid(m)))
%                         if fct.wlid(m) <= -99
%                             fct.wltp(n) = 3; %bounding wall
%                         else
%                             fct.wltp(n) = 4; %floor
%                         end
%                     else
%                         disp(['ERROR: walltype of facet ',num2str(m),' not defined'])
%                     end
%                 end
%                 
%                 
%                 for n = 1:length(idref)
%                     if (ismember(idref(n), wall.id))
%                         fct.wltp(fct.wlid == idref(n)) = 1;
%                     elseif (ismember(idref(n), gr.id))
%                         fct.wltp(fct.wlid == idref(n)) = 2;
%                     else
%                         sprintf('ERROR. wall id %8d not defined. Correct wall and greenroof property input files', idref(n))
%                         return
%                     end
%                 end
%             end
%             
%             
%             
%             %% iterate shortwave and longwave radiation
%             %% iterate over every facet and over entire system until energy is balanced
%                         
%             %% derived quantities
%             
%             [fct, wall] = loadfacets();
%             [sortt, sorti] = sort(obj.facets(:, 1));  %sort by walltype
%             
%             %% shortwave
%             %
%             
%             %Z=35;          %zenith angle of the sun (could be function of time, location etc)
%             Dsky = zeros(obj.nfcts, 1);
%             Denv = zeros(obj.nfcts, 1);
%             
%             da_pp.addvar(obj, 'albedo', zeros(obj.nfcts, 1));
%             da_pp.addvar(obj, 'emissivity', zeros(obj.nfcts, 1));
%             
%             for i = 1:obj.nfcts
%                 j = find(fct.wlid(i)==wall.id);
%                 obj.albedo(i) = wall.al(j);
%                 obj.emissivity(i) = wall.em(j);
%             end
%             
%             %isroof
%             isnotroof = ones(obj.nfcts,1);
%             isnotroof(find(obj.facets(:, 1) == 1 & obj.facets(:,4) > 0)) = 0;
%             
%             %diffuse flux from sky and other walls
%             Kinnew = zeros(obj.nfcts, 1);
%             
%             da_pp.addvar(obj, 'Kininit', (1 - obj.albedo) .* (obj.Sdir + obj.Dsk .* obj.svf));
%             Koutinit = obj.albedo .* (obj.Sdir + obj.Dsk .* obj.svf);
%             
% %             if ltestplot
% %                 figure
% %                 a=(1-albedo).*Sdir;
% %                 b=(1-albedo).*Dsk.*svf;
% %                 plot(1:nfcts,a(sorti),1:nfcts,b(sorti))
% %                 ax1 = gca; % current axes
% %                 set(ax1,'Xtick',1:1:nfcts)
% %                 xticklabels(sortt)
% %                 ax1_pos = ax1.Position; % position of first axes
% %                 xlabel('orientation')
% %                 ylabel('S and D')
% %                 ax2 = axes('Position',ax1_pos,...
% %                     'YAxisLocation', 'right',...
% %                     'XAxisLocation', 'top', ...
% %                     'Color', 'None');
% %                 set (ax2, 'XLim', get (ax1, 'XLim'), 'Layer', 'top');
% %                 set (ax2, 'YLim', get (ax1, 'YLim'), 'Layer', 'top');
% %                 set(ax2,'Xtick',1:1:nfcts)
% %                 set(ax2,'Xticklabels',sorti)
% %                 set(ax2,'FontSize',8)
% %                 xlabel('facet index')
% %             end
%             
%             %total radiation absorbed and reflected
%             sum(obj.Kininit + Koutinit);
%             %total radiation absorbed
%             sum(obj.Kininit);
%             %total radiation reflected
%             sum(Koutinit);
%             
%             
%             
%             totincrease = zeros(10, 1);
%             increase = zeros(10, 1);
%             Kout = zeros(10, 1);
%             Kout(2) = sum(Koutinit);
%             
%             
%             da_pp.addvar(obj, 'Kin', obj.Kininit); %Kin §adds up
%             Kinold = obj.Kininit;
%             Koutold = Koutinit; %Kout is wiped clean with every reflection
%             Koutnew = zeros(obj.nfcts, 1);
%             Kintemp = zeros(obj.nfcts, 1);
%             Kouttemp = zeros(obj.nfcts, 1);
%             count = 0;
%             storerad = zeros(obj.nfcts, obj.nfcts);
%             Kiterin = zeros(300, 1);
%             Kiterout = zeros(300, 1);
%             itermaxdiff = zeros(300, 1);
%             itermaxdiffloc = zeros(300, 1);
%             itermaxrelloc = zeros(300, 1);
%             moep = zeros(obj.nfcts, 20);
%             blub = zeros(obj.nfcts, 20);
%             moep(:, 1) = obj.Kininit;
%             
%             A = obj.facetarea;
%             albedo = obj.albedo;
%             vfo = obj.vfo;
%             while true
%                 count = count + 1;
%                 for i = 1:obj.nfcts
%                     Kintemp(i) = 0;
%                     Kouttemp(i) = 0;
%                     for j = 1:obj.nfcts %sum all the incoming radiation on i, originally reflected from all the other j facets ("radiation reflected on j" x "what perecentage does i take of j's vision")
%                         inc = Koutold(j) * A(j) / A(i) * vfo(j,i);  %[W/m2]
%                         storerad(i, j) = inc;
%                         Kintemp(i) = Kintemp(i) + (1 - albedo(i)) * inc;
%                         Kouttemp(i) = Kouttemp(i) + albedo(i) * inc;
%                     end
%                     obj.Kin(i) = obj.Kin(i) + Kintemp(i); %add newly absorbed radiation to already existing one
%                     Koutnew(i) = Kouttemp(i); %save newly reflected radiation for next iteration
%                 end
%                 % Kiterin(count)=Kin(1);
%                 % Kiterout(count)=Koutnew(1);
%                 % [itermaxdiff(count),itermaxdiffloc(count)]=max(Kintemp(i));
%                 % [~,itermaxrelloc(count)]=max(Koutnew./Koutold);
%                 % max(Koutnew./Koutold)
%                 % if all(Koutnew./Koutold<0.01)
%                 %     break
%                 % end
%                 moep(:, count + 1) = obj.Kin;
%                 if (max((obj.Kin - Kinold) ./ Kinold) < 0.01)
%                     %disp(['reached relative tolerance after ' num2str(count) ' iterations'])
%                     break
%                 end
%                 %    if (max(Koutnew-Koutold)<0.01 )
%                 %        disp('reached absolute tolerance')
%                 %        break
%                 %    end
%                 Kinold = obj.Kin;
%                 Koutold = Koutnew; %overwrite reflected radiation with new value
%             end
%             %%
%             %if lhqplot
% %                 scale=2;
% %                 scalef=1.5;
% %                 h= figure;
% %                 set(gcf,'units','centimeters','position',[0 0 14.5*scale 14.5*scale]);
% %                 set(h,'PaperPosition',[0 0 14.5*scale 14.5*scale]);
% %                 set(h,'PaperUnits','centimeters');
% %                 set(h,'renderer','painters');
% %                 select=140:-1:6;
% %                 select=[1:5 select];
% %                 hs1=subplot(1,2,1)
% %                 hold on
% %                 stylez={':x';':o';':d';':s';':+'};
% %                 cm=hsv;
% %                 colint=floor(length(cm)/5);
% %                 for k=1:length(select)
% %                     h=select(k);
% %                     plot(0:1:count,moep(h,1:count+1),stylez{obj.facets(h,1)},'color',cm(1+obj.facets(h,1)*colint,:),'Linewidth',2)
% %                 end
% %                 xlim([0 count])
% %                 xlabel('iteration $n$','Interpreter','latex','FontSize',18)
% %                 ylabel('$K^{\downarrow}$ [Wm$^{-2}$]','Interpreter','latex','FontSize',18)
% %                 xticks(0:count)
% %                 xticklabels({'init' '1' '2' '3' '4' '5'})
% %                 set(hs1,'TickLabelInterpreter','latex')
% %                 set(hs1,'FontSize',12)
% %                 hs=subplot(1,2,2);
% %                 for r=2:count+1
% %                     blub(:,r)=((moep(:,r)-moep(:,r-1))./moep(:,r-1));
% %                 end
% %                 blub(:,1)=NaN;
% %                 hold on
% %                 for k=1:length(select)
% %                     h=select(k);
% %                     plot(0:1:count,blub(h,1:count+1),stylez{obj.facets(h,1)},'color',cm(1+obj.facets(h,1)*colint,:),'Linewidth',2)
% %                 end
% %                 plot([0 count],[0.01 0.01],'k:','Linewidth',2)
% %                 xlim([0 count])
% %                 xlabel('iteration $n$','Interpreter','latex','FontSize',18)
% %                 xticks(0:count)
% %                 xticklabels({'init' '1' '2' '3' '4' '5'})
% %                 ylabel('$(K^{\downarrow}_n-K^{\downarrow}_{n-1})/K^{\downarrow}_{n-1}$ [-]','Interpreter','latex','FontSize',18)
% %                 hs.YScale='log';
% %                 ylim([0.0001 10])
% %                 yticklabels({'0.0001' '0.001' '0.01' '0.1' '1' '10'})
% %                 set(hs,'TickLabelInterpreter','latex')
% %                 set(hs,'FontSize',12)
% %                 set(gcf, 'Color', 'w');
% %                 legend('roof/road','west','east','north','south','Location','NorthEast')
% %                 %export_fig rayswcov.eps
% %             %end
% %             
% %             %if lhqplot
% %                 scale=2;
% %                 scalef=1.5;
% %                 h= figure;
% %                 set(gcf,'units','centimeters','position',[0 0 14.5*scale 14.5*scale]);
% %                 set(h,'PaperPosition',[0 0 14.5*scale 14.5*scale]);
% %                 set(h,'PaperUnits','centimeters');
% %                 set(h,'renderer','painters');
% %                 a=(1-obj.albedo).*obj.Sdir;
% %                 b=(1-obj.albedo).*obj.Dsk.*obj.svf;
% %                 plot(1:obj.nfcts,obj.Kin(sorti),'x',1:obj.nfcts, obj.Kininit(sorti),'+',1:obj.nfcts,a(sorti),'o',1:obj.nfcts,b(sorti),'d')
% %                 legend('Total K','Initial K','Initial Sdir','Initial D')
% %                 ax1 = gca; % current axes
% %                 set(ax1,'Xtick',1:1:obj.nfcts)
% %                 set(ax1,'FontSize',9)
% %                 xticklabels(sortt)
% %                 ax1_pos = ax1.Position; % position of first axes
% %                 xlabel('orientation','Interpreter','latex','FontSize',12)
% %                 ylabel('Radiation [$W/m^2$]','Interpreter','latex','FontSize',12)
% %                 ax2 = axes('Position',ax1_pos,...
% %                     'YAxisLocation', 'right',...
% %                     'XAxisLocation', 'top', ...
% %                     'Color', 'None');
% %                 set (ax2, 'XLim', get (ax1, 'XLim'), 'Layer', 'top');
% %                 set (ax2, 'YLim', get (ax1, 'YLim'), 'Layer', 'top');
% %                 set(ax2,'YTickLabel','')
% %                 set(ax2,'Xtick',1:3:obj.nfcts)
% %                 set(ax2,'Xticklabels',sorti(1:3:obj.nfcts),'TickLabelInterpreter','latex')
% %                 set(ax2,'FontSize',9)
% %                 xl2=xlabel('facet index','FontSize',12);
% %                 xl2.Position(2)=xl2.Position(2)+15;
% %                 offset = repmat(ax2.YTick(end)+20,1,numel(ax2.XTick));
% %                 % create new lables:
% %                 % text(ax1.XTick(2:3:nfcts),offset,num2str(sorti(2:3:nfcts)),'HorizontalAlign','center','FontSize',9,'Interpreter','latex')
% %                 
% %                 offset = repmat(ax2.YTick(end)+30,1,numel(ax2.XTick)-1);
% %                 % create new lables:
% %                 %   text(ax1.XTick(3:3:nfcts),offset,num2str(sorti(3:3:nfcts)),'HorizontalAlign','center','FontSize',9,'Interpreter','latex')
% %                 set(gcf, 'Color', 'w');
% %                 %export_fig shortwaveit.eps
% %                 
% %                 
% %                 scale=2;
% %                 scalef=1.5;
% %                 h= figure;
% %                 set(gcf,'units','centimeters','position',[0 0 14.5*scale 14.5*scale]);
% %                 set(h,'PaperPosition',[0 0 14.5*scale 14.5*scale]);
% %                 set(h,'PaperUnits','centimeters');
% %                 set(h,'renderer','painters');
% %                 a=(1-obj.albedo).*obj.Sdir;
% %                 b=(1-obj.albedo).*obj.Dsk.*obj.svf;
% %                 plot(1:obj.nfcts,blub(sorti,2),'x',1:obj.nfcts,blub(sorti,3),'+',1:obj.nfcts,blub(sorti,4),'o',1:obj.nfcts,blub(sorti,5),'d',1:obj.nfcts,blub(sorti,6),'s')
% %                 legend('Total K','Initial K','Initial Sdir','Initial D')
% %                 ax1 = gca; % current axes
% %                 ax1.YScale='log';
% %                 set(ax1,'Xtick',1:1:obj.nfcts)
% %                 set(ax1,'FontSize',9)
% %                 xticklabels(sortt)
% %                 ax1_pos = ax1.Position; % position of first axes
% %                 xlabel('orientation','Interpreter','latex','FontSize',12)
% %                 ylabel('Radiation [$W/m^2$]','Interpreter','latex','FontSize',12)
% %                 ax2 = axes('Position',ax1_pos,...
% %                     'YAxisLocation', 'right',...
% %                     'XAxisLocation', 'top', ...
% %                     'Color', 'None');
% %                 set (ax2, 'XLim', get (ax1, 'XLim'), 'Layer', 'top');
% %                 set (ax2, 'YLim', get (ax1, 'YLim'), 'Layer', 'top');
% %                 set(ax2,'YTickLabel','')
% %                 %set(ax2,'Xtick',1:3:nfcts)
% %                 %set(ax2,'Xticklabels',sorti(1:3:nfcts),'TickLabelInterpreter','latex')
% %                 set(ax2,'Xticklabels','','TickLabelInterpreter','latex')
% %                 set(ax2,'FontSize',9)
% %                 xl2=xlabel('facet index','FontSize',12);
% %                 xl2.Position(2)=xl2.Position(2)+15;
% %                 offset = repmat(ax2.YTick(end)+20,1,numel(ax2.XTick));
% %                 % create new lables:
% %                 % text(ax1.XTick(2:3:nfcts),offset,num2str(sorti(2:3:nfcts)),'HorizontalAlign','center','FontSize',9,'Interpreter','latex')
% %                 
% %                 offset = repmat(ax2.YTick(end)+30,1,numel(ax2.XTick)-1);
% %                 % create new lables:
% %                 %   text(ax1.XTick(3:3:nfcts),offset,num2str(sorti(3:3:nfcts)),'HorizontalAlign','center','FontSize',9,'Interpreter','latex')
% %                 set(gcf, 'Color', 'w');
% %                 %export_fig shortwaveit2.eps
% %                 
% %                 % ax3 = axes('Position',ax1_pos,...
% %                 %     'YAxisLocation', 'right',...
% %                 %     'XAxisLocation', 'top', ...
% %                 %     'Color', 'None');
% %                 % set (ax3, 'XLim', get (ax1, 'XLim'), 'Layer', 'top');
% %                 % set (ax3, 'YLim', get (ax1, 'YLim'), 'Layer', 'top');
% %                 % ax3.YTickLabel='';
% %                 % ax3.YTick=''
% %                 % set(ax3,'Xtick',2:2:nfcts)
% %                 % set(ax3,'Xticklabels',sorti(2:2:nfcts))
% %                 % set(ax3,'FontSize',8)
% %                 % box(ax3(1),'off')
% %                 % ax3.TickLength=[0;0]
% %                 %ax3.Position(4)=ax3.Position(4)*1.02;
% %                 
% %                 
% %                 
% %                 
% %                 
% %                 figure
% %                 plot(obj.Kin(sorti)./obj.Kininit(sorti))
% %                 ax1 = gca; % current axes
% %                 set(ax1,'Xtick',1:1:obj.nfcts)
% %                 xticklabels(sortt)
% %                 ax1_pos = ax1.Position; % position of first axes
% %                 xlabel('orientation')
% %                 ylabel('K/Kinit')
% %                 ax2 = axes('Position',ax1_pos,...
% %                     'YAxisLocation', 'right',...
% %                     'XAxisLocation', 'top', ...
% %                     'Color', 'None');
% %                 set (ax2, 'XLim', get (ax1, 'XLim'), 'Layer', 'top');
% %                 set (ax2, 'YLim', get (ax1, 'YLim'), 'Layer', 'top');
% %                 set(ax2,'Xtick',1:1:obj.nfcts)
% %                 set(ax2,'Xticklabels',sorti)
% %                 set(ax2,'FontSize',8)
% %                 xlabel('facet index')
% %                 
% %                 figure
% %                 plot(obj.Kin(sorti)-obj.Kininit(sorti))
% %                 ax1 = gca; % current axes
% %                 set(ax1,'Xtick',1:1:obj.nfcts)
% %                 xticklabels(sortt)
% %                 ax1_pos = ax1.Position; % position of first axes
% %                 xlabel('orientation')
% %                 ylabel('K-Kinit')
% %                 ax2 = axes('Position',ax1_pos,...
% %                     'YAxisLocation', 'right',...
% %                     'XAxisLocation', 'top', ...
% %                     'Color', 'None');
% %                 set (ax2, 'XLim', get (ax1, 'XLim'), 'Layer', 'top');
% %                 set (ax2, 'YLim', get (ax1, 'YLim'), 'Layer', 'top');
% %                 set(ax2,'Xtick',1:1:obj.nfcts)
% %                 set(ax2,'Xticklabels',sorti)
% %                 set(ax2,'FontSize',8)
% %                 xlabel('facet index')
% %                 
%             %end
%         end

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
        albedo=zeros(nfcts,1);
        emissivity=zeros(nfcts,1);
        Sdir = obj.Sdir;
        Dsk = obj.Dsk;
        svf = obj.svf;
        
        for i=1:nfcts
            j=find(fct.wlid(i)==wall.id);
            albedo(i)=wall.al(j);
            emissivity(i)=wall.em(j);
        end
        da_pp.addvar(obj, 'emissivity', emissivity);
        
        %isroof
        isnotroof=ones(nfcts,1);
        isnotroof(find(facets(:,1)==1 & facets(:,4)>0))=0;
        
        %diffuse flux from sky and other walls
        Kinnew=zeros(nfcts,1);
        
        Kininit=(1-albedo).*(Sdir+Dsk.*svf);
        Koutinit=albedo.*(Sdir+Dsk.*svf);
        
        %     if ltestplot
        %         figure
        %         a=(1-albedo).*Sdir;
        %         b=(1-albedo).*Dsk.*svf;
        %         plot(1:nfcts,a(sorti),1:nfcts,b(sorti))
        %         ax1 = gca; % current axes
        %         set(ax1,'Xtick',1:1:nfcts)
        %         xticklabels(sortt)
        %         ax1_pos = ax1.Position; % position of first axes
        %         xlabel('orientation')
        %         ylabel('S and D')
        %         ax2 = axes('Position',ax1_pos,...
        %             'YAxisLocation', 'right',...
        %             'XAxisLocation', 'top', ...
        %             'Color', 'None');
        %         set (ax2, 'XLim', get (ax1, 'XLim'), 'Layer', 'top');
        %         set (ax2, 'YLim', get (ax1, 'YLim'), 'Layer', 'top');
        %         set(ax2,'Xtick',1:1:nfcts)
        %         set(ax2,'Xticklabels',sorti)
        %         set(ax2,'FontSize',8)
        %         xlabel('facet index')
        %     end
        
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
            moep(:,count+1)=Kin;
            if (max((Kin-Kinold)./Kinold)<0.01)
                disp(['reached relative tolerance after ' num2str(count) ' iterations'])
                break
            end
            %    if (max(Koutnew-Koutold)<0.01 )
            %        disp('reached absolute tolerance')
            %        break
            %    end
            Kinold=Kin;
            Koutold=Koutnew; %overwrite reflected radiation with new value
        end
    
        da_pp.addvar(obj, 'Kin', Kin);
        end
        
        
        function write_netsw(obj)
            fname = ['netsw.inp.' obj.expnr];
            fileID = fopen(fname, 'w');
            fprintf(fileID,'# %4s\n','net shortwave on facets [W/m2] (including reflections and diffusive)');
            fprintf(fileID,'%6d\n', obj.Kin);
            fclose(fileID);   
        end
        
        
        function generate_Tfacinit(obj, iss)
            %% Inital Temperature and longwave
            %solve energy budget equation in an initial steady state
            %assume 0 wall heatflux
            %assume 0 latent heatflux
            %assume constant air temperature and heat transfer coefficient
            %assume constant longwave
            %solve for Tinit and Linit
            %assume absorbtivity is equal to emissivity
            %=>
            %K+L=H
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
                Lsk = 350;
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
            
            
                da_pp.addvar(obj, 'Tfacinit', Tnew(:, k));
                        
            else 
                da_pp.addvar(obj, 'Tfacinit', 288 * ones(obj.nfcts, 1));
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
            fprintf(blocks, '#  il\t   iu\t   jl\t   ju\t   kl\t   ku\t dtop\t dwest\t deast\t dnor\t dsou\n');
            fprintf(blocks, '%4d\t%4d\t%4d\t%4d\t%4d\t%4d\t%4d\t%4d\t%4d\t%4d\t%4d\n', obj.blocks');
            fclose(blocks);
            %disp(['... written blocks.inp.' obj.expnr])
        end
        
        function write_facets(obj)
            fname = ['facets.inp.' num2str(obj.expnr)];
            fileID = fopen(fname,'w');
            fprintf(fileID,'# %4s %6s %6s %6s\n','or', 'wl', 'blk', 'bld');
            fprintf(fileID,'%6d %6d %6d %6d\n', obj.facets(:, 1:4)');
            fclose(fileID); 
        end
        
        function generate_trees(obj, lwritefile)
            da_pp.addvar(obj, 'nrows', obj.imax / (obj.blockwidth + obj.canyonwidth));
            if obj.lblocks
                if obj.ltrees
                   % Implement trees - haven't done so far because in
                   % da_inp, xh is being changed.
                   da_pp.addvar(obj, 'trees', zeros(obj.nrows, 6));
                end
            end
            
            if lwritefile
                if obj.ltrees
                    da_pp.write_trees(obj)
                end
            end
        end
        
        function write_trees(obj)
            tree_write = fopen( ['trees.inp.' obj.expnr], 'w');
            fprintf(tree_write, '%-12s\n', '# Tree location');
            fprintf(tree_write, '%-60s\n', '#  il  iu  jl  ju  kl  ku');
            fprintf(tree_write, '%-3.0f %-3.0f %-3.0f %-3.0f %-3.0f %-3.0f\n', obj.trees');
            fclose(tree_write);
            %disp(['... written trees.inp.' obj.expnr]) 
        end
        
        function generate_purifs(obj)
            % Uses obj.bl, so should be run after having run
            % generate_blocks (or reading in a blocks file)
            da_pp.addvar(obj, 'nrows', obj.imax / (obj.blockwidth + obj.canyonwidth));
            if obj.lblocks
                if obj.lpurif
                    %purifs = zeros(obj.nrows * 2 * obj.npurif, 7);
                    da_pp.addvar(obj, 'purifs', zeros(obj.nrows * 2 * obj.npurif, 7));
                    for i = 1:obj.nrows
                        for j = 1:obj.npurif
                            obj.purifs((i - 1) * obj.npurif + j,1) = obj.bl(obj.nrows + i, 2) - obj.purif_dx - obj.purif_w;
                            obj.purifs((i - 1) * obj.npurif + j,2) = obj.bl(nrows+i,2) - obj.purif_dx;
                            if j==1
                                %obj.purifs((i - 1) * obj.npurif + j, 3) = ((je/npurif)/2);
                                obj.purifs((i - 1) * obj.npurif + j, 3) = ((obj.jtot / obj. npurif) / 2);
                                obj.purifs((i - 1) * obj.npurif + j, 4) = obj.purifs((i - 1) * obj.npurif + j, 3) + obj.purif_dy;
                            else
                                obj.purifs((i - 1) * obj.npurif + j, 3) = obj.purifs((i - 1) * obj.npurif + j - 1, 3) + obj.purif_dy + obj.purif_sp;
                                obj.purifs((i - 1) * obj.npurif + j, 4) = obj.purifs((i - 1) * obj.npurif + j - 1, 4) + obj.purif_dy + obj.purif_sp;
                            end
                        end
                    end
                    for i = 1:obj.nrows
                        for j = 1:obj.npurif
                            obj.purifs(obj.nrows * obj.npurif + (i - 1) * obj.npurif + j,1) =  obj.bl(obj.nrows + 1 + i, 1) + obj.purif_dx;
                            obj.purifs(obj.nrows * obj.npurif + (i - 1) * obj.npurif + j,2) =  obj.bl(obj.nrows + 1 + i, 1) + obj.purif_dx + obj.purif_w;
                            if j == 1
                                obj.purifs(obj.nrows * obj.npurif + (i - 1) * obj.npurif + j, 3) = ((obj.tot / obj.npurif) / 2);
                                obj.purifs(obj.nrows * obj.npurif + (i - 1) * obj.npurif + j, 4) = obj.purifs(( i - 1) * obj.npurif + j, 3) + obj.purif_dy;
                            else
                                obj.purifs(obj.nrows * obj.npurif + (i - 1) * obj.npurif + j, 3) = obj.purifs((i - 1) * obj.npurif + j - 1, 3) + obj.purif_dy + obj.purif_sp;
                                obj.purifs(obj.nrows * obj.npurif + (i - 1) * obj.npurif + j, 4) = obj.purifs((i - 1) * obj.npurif + j - 1, 4) + obj.purif_dy + obj.purif_sp;
                            end
                        end
                    end

                    obj.purifs(:,5) = obj.purif_dz;
                    obj.purifs(:,6) = obj.purif_dz + obj.purif_h;
                    obj.purifs(:,7) = obj.purif_i;
                end
            end
            
            if lwritefile
                if obj.lpurif
                    da_pp.write_purifs(obj)
                end
            end
        end
        
        function write_purifs(obj)
            purif_write = fopen( ['purifs.inp.' obj.expnr], 'w');
            fprintf(purif_write, '%-12s\n', '# Purifier location');
            fprintf(purif_write, '%-60s\n', '#  il  iu  jl  ju  kl  ku  ipu');
            fprintf(purif_write, '%-3.0f %-3.0f %-3.0f %-3.0f %-3.0f %-3.0f %-3.0f\n', obj.purifs');
            fclose(purif_write);
            %disp(['... written purifs.inp.' obj.expnr]) 
        end
        
        function plot_domain(obj)
            figure
            view(52, 23)
            if (obj.lcastro || obj.lcube || obj.lblocks)
                for i = 1:size(obj.bl, 1)
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1)  obj.xh(obj.bl(i,1))], [obj.yh(obj.bl(i,3))  obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.zh(obj.bl(i,6)+1)  obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1)], [245 245 245] ./ 255)
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,3))  obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3))], [obj.bl(i,5)  obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], [245 245 245] ./ 255)
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.bl(i,5)  obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], [245 245 245] ./ 255)
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) ], [obj.yh(obj.bl(i,4)+1)  obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3))], [obj.bl(i,5)  obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], [245 245 245] ./ 255)
                    patch([obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.bl(i,5)  obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], [245 245 245] ./ 255)
                    % patch([xh(1) xh(end) xh(end)  xh(1)], [yh(1)  yh(1) yh(end) yh(end)], [zh(1)  zh(1) zh(1) zh(1)], [245 245 245] ./ 255)
                end

            elseif obj.lflat

            else
                for i = 1:size(obj.bl, 1)
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,1))], [obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1)], 'w')
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3))], [obj.bl(i,5) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], 'w')
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.bl(i,5) obj.zh(bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], 'w')
                    patch([obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1)) obj.xh(obj.bl(i,1))], [obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3))], [obj.bl(i,5) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], 'w')
                    patch([obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) obj.xh(obj.bl(i,2)+1) ], [obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,3)) obj.yh(obj.bl(i,4)+1) obj.yh(obj.bl(i,4)+1)], [obj.bl(i,5) obj.zh(obj.bl(i,6)+1) obj.zh(obj.bl(i,6)+1) obj.bl(i,5)], 'w')
                end
            end

            if obj.ltrees
                for i = 1:size(obj.trees,1)
                    patch([obj.xh(obj.trees(i,1)) obj.xh(obj.trees(i,2)+1) obj.xh(obj.trees(i,2)+1)  obj.xh(obj.trees(i,1))], [obj.yh(obj.trees(i,3)) obj.yh(obj.trees(i,3)) obj.yh(obj.trees(i,4)+1) obj.yh(obj.trees(i,4)+1)], [obj.zh(obj.trees(i,6)+1) obj.zh(obj.trees(i,6)+1) obj.zh(obj.trees(i,6)+1) obj.zh(obj.trees(i,6)+1)], [128 229 78] ./ 255)
                    patch([obj.xh(obj.trees(i,1)) obj.xh(obj.trees(i,1)) obj.xh(obj.trees(i,2)+1) obj.xh(obj.trees(i,2)+1) ], [obj.yh(obj.trees(i,3)) obj.yh(obj.trees(i,3)) obj.yh(obj.trees(i,3)) obj.yh(obj.trees(i,3))], [obj.zh(obj.trees(i,5)) obj.zh(obj.trees(i,6)+1) obj.zh(obj.trees(i,6)+1) obj.zh(obj.trees(i,5))], [128 229 78] ./ 255)
                    patch([obj.xh(obj.trees(i,1)) obj.xh(obj.trees(i,1)) obj.xh(obj.trees(i,2)+1) obj.xh(obj.trees(i,2)+1) ], [obj.yh(obj.trees(i,4)+1) obj.yh(obj.trees(i,4)+1) obj.yh(obj.trees(i,4)+1) obj.yh(obj.trees(i,4)+1)], [obj.zh(obj.trees(i,5)) obj.zh(obj.trees(i,6)+1) obj.zh(obj.trees(i,6)+1) obj.zh(obj.trees(i,5))],[128 229 78] ./ 255 )
                    patch([obj.xh(obj.trees(i,1)) obj.xh(obj.trees(i,1)) obj.xh(trees(i,1)) obj.xh(trees(i,1)) ], [obj.yh(obj.trees(i,4)+1) obj.yh(obj.trees(i,4)+1) obj.yh(obj.trees(i,3)) obj.yh(obj.trees(i,3))], [obj.zh(obj.trees(i,5)) obj.zh(obj.trees(i,6)+1) obj.zh(obj.trees(i,6)+1) obj.zh(obj.trees(i,5))],[128 229 78] ./ 255 )
                    patch([obj.xh(obj.trees(i,2)+1) obj.xh(obj.trees(i,2)+1) obj.xh(obj.trees(i,2)+1) obj.xh(obj.trees(i,2)+1)], [obj.yh(obj.trees(i,3)) obj.yh(obj.trees(i,3)) obj.yh(obj.trees(i,4)+1) obj.yh(obj.trees(i,4)+1)], [obj.zh(obj.trees(i,5)) obj.zh(obj.trees(i,6)+1) obj.zh(obj.trees(i,6)+1) obj.zh(obj.trees(i,5))], [128 229 78] ./ 255)

                end
            end

            if obj.lpurif
                for i = 1:size(obj.purifs,1)
                    patch([obj.xh(obj.purifs(i,1)) obj.xh(obj.purifs(i,2)+1) obj.xh(obj.purifs(i,2)+1) obj.xh(obj.purifs(i,1))], [obj.yh(obj.purifs(i,3)) obj.yh(obj.purifs(i,3)) obj.yh(obj.purifs(i,4)+1) obj.yh(obj.purifs(i,4)+1)], [obj.zh(obj.purifs(i,6)+1) obj.zh(obj.purifs(i,6)+1) obj.zh(obj.purifs(i,6)+1) obj.zh(obj.purifs(i,6)+1)], [245 245 245] ./ 255)
                    patch([obj.xh(obj.purifs(i,1)) obj.xh(obj.purifs(i,1)) obj.xh(obj.purifs(i,2)+1) obj.xh(obj.purifs(i,2)+1)], [obj.yh(obj.purifs(i,3)) obj.yh(obj.purifs(i,3)) obj.yh(purifs(i,3)) obj.yh(obj.purifs(i,3))], [obj.zh(obj.purifs(i,5)) obj.zh(obj.purifs(i,6)+1) obj.zh(obj.purifs(i,6)+1) obj.zh(obj.purifs(i,5))], [245 245 245] ./ 255)
                    patch([obj.xh(obj.purifs(i,1)) obj.xh(obj.purifs(i,1)) obj.xh(obj.purifs(i,2)+1) obj.xh(obj.purifs(i,2)+1)], [obj.yh(obj.purifs(i,4)+1) obj.yh(obj.purifs(i,4)+1) obj.yh(obj.purifs(i,4)+1) obj.yh(obj.purifs(i,4)+1)], [obj.zh(obj.purifs(i,5)) obj.zh(obj.purifs(i,6)+1) obj.zh(obj.purifs(i,6)+1) obj.zh(obj.purifs(i,5))], [245 245 245] ./ 255)
                    patch([obj.xh(obj.purifs(i,1)) obj.xh(obj.purifs(i,1)) obj.xh(obj.purifs(i,1)) obj.xh(obj.purifs(i,1))], [obj.yh(obj.purifs(i,4)+1) obj.yh(obj.purifs(i,4)+1) obj.yh(obj.purifs(i,3)) obj.yh(obj.purifs(i,3))], [obj.zh(obj.purifs(i,5)) obj.zh(obj.purifs(i,6)+1) obj.zh(obj.purifs(i,6)+1) obj.zh(obj.purifs(i,5))], [245 245 245] ./ 255)
                    patch([obj.xh(obj.purifs(i,2)+1) obj.xh(obj.purifs(i,2)+1) obj.xh(obj.purifs(i,2)+1) obj.xh(obj.purifs(i,2)+1)], [obj.yh(obj.purifs(i,3)) obj.yh(obj.purifs(i,3)) obj.yh(obj.purifs(i,4)+1) obj.yh(obj.purifs(i,4)+1)], [obj.zh(obj.purifs(i,5)) obj.zh(obj.purifs(i,6)+1) obj.zh(obj.purifs(i,6)+1) obj.zh(obj.purifs(i,5))], [245 245 245] ./ 255)

                end
            end

            zlim([0 obj.zh(end)]); %/(r.blockheight-1))
            xlim([0 obj.xh(end)]); %/(r.blockheight-1))
            ylim([0 obj.yh(end)]); %/(r.blockheight-1))

            set(gca,'ticklabelinterpreter','latex')
            xlabel('x [m]','interpreter','latex')
            ylabel('y [m]','interpreter','latex')
            zlabel('z [m]','interpreter','latex')
            set(gca,'BoxStyle','full','Box','on')
            daspect([1 1 1])
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
                            
                            [flag,dint] = da_pp.rbi(facetpoint, v1, bl(n,:));
                            
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
                                
                                [flag,dint] = da_pp.rbi(facetpoint, v1, bl(m+nblocks,:));
                                
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
                        
                        [flag,dint] = da_pp.rbi(facetpoint, v1, bl(n,:));
                        
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
                            
                            [flag,dint] = da_pp.rbi(facetpoint, v1, bl(k+nblocks,:));
                            
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
                                
                                [flag,dints] = da_pp.rbi(facetpoint, v1, bl(n,:));
                                
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
                            
                            [flag,dints] = da_pp.rbi(facetpoint, v1, bl(n,:));
                            
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
        
        function [ ndim, area, co] = detsub(i,facets,blocks,floors,boundingwalls,cornm,xb,yb,zb,delta)
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
                il=floors(bi,1);     %lower x index of floor facet 1
                iu=floors(bi,2);     %upper x index of floor facet 1
                jl=floors(bi,3);     %lower y index of floor facet 1
                ju=floors(bi,4);     %upper y index of floor facet 1
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