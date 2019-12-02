classdef da_pp < dynamicprops
    % Class for storing input parameters to uDALES.
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
                        % number or logical
                        obj.(lhs) = str2double(rhs);
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
            
            if ceil(obj.jtot / ncpus) ~= floor (obj.jtot / ncpus)
                disp(['Possible jtot: ' num2str([2 3 4 5 6 7 8] * ncpus)])
                error('No. CPUs does not fit j grid size')
            else
                disp('cpus and jtot successful')
            end
            
            % Blocks
            da_pp.addvar(obj, 'lcastro', 0) % switch for staggered cubes
            da_pp.addvar(obj, 'lcube', 0)   % switch for linear cubes
            da_pp.addvar(obj, 'lblocks', 0) % switch for infinite blocks
            if (not(obj.lcastro) && not(obj.lcube) && not(obj.lblocks))
                da_pp.addvar(obj, 'lflat',1)
                disp('No standard block config. setup so flat domain assumed.')
            else
                da_pp.addvar(obj, 'lflat',0)
            end

            da_pp.addvar(obj, 'blockheight', 16) % block height
            da_pp.addvar(obj, 'blockwidth', 16)  % block width
            da_pp.addvar(obj, 'canyonwidth', 16) % canyonwidth
            
            % Profiles
            da_pp.addvar(obj, 'lmassflowr', 0) % switch for constant mass flow rate
            da_pp.addvar(obj, 'lprofforc', 0)  % switch for 1D geostrophic forcing
            da_pp.addvar(obj, 'lcoriol', 0)    % switch for coriolis forcing
            if (not(obj.lmassflowr) && not(obj.lprofforc) && not(obj.lcoriol))
                da_pp.addvar(obj, 'ldp', 1)
                disp('No forcing switch config. setup so initial velocities and pressure gradients applied.')
            else
                da_pp.addvar(obj, 'ldp', 0)
            end

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
                        disp('Purifier layoud does not fit grid')
                        disp(['sum widths to: ' num2str(errp(1,:))])
                        disp(['Current width: ' num2str(obj.purif_dy + obj.purif_sp)])
                        error('Incorrect purifier layout')
                    else
                        disp('Successful purifer layout')
                        disp(['Number of purifiers: ' num2str(npurif*2+1)])
                    end
                else
                    error('Must use lblocks configuration to use purifiers')
                end
            end
        end
        
        function generate_xygrid(obj, lwritefile)
             da_pp.addvar(obj, 'xf', 0.5 * obj.dx : obj.dx : obj.xsize - 0.5 * obj.dx); 
             da_pp.addvar(obj, 'yf', 0.5 * obj.dy : obj.dy : obj.ysize - 0.5 * obj.dy);
             da_pp.addvar(obj, 'xh', 0 : obj.dx : obj.xsize);
             da_pp.addvar(obj, 'yh', 0 : obj.dy : obj.ysize);
             
             if lwritefile
                da_pp.write_xgrid(obj);
             end
        end
        
        function write_xgrid(obj)
            xgrid = fopen(['xgrid.inp.' obj.expnr], 'w');
            fprintf(xgrid, '%12s\n', '#     x-grid');
            fprintf(xgrid, '%12s\n', '#           ');
            fprintf(xgrid, '%-20.15f\n', obj.xf);
            fclose(xgrid);
            disp(['... written xgrid.inp.' obj.expnr]) 
        end
        
        function generate_zgrid(obj, lwritefile)
            if ~obj.lzstretch
                da_pp.addvar(obj, 'zf', 0.5 * obj.dz : obj.dz : obj.zsize - 0.5 * obj.dz);
                da_pp.addvar(obj, 'zh', 0 : obj.dz : obj.zsize);
            else
                % Implement z-stretching here
            end
            
            if lwritefile
                da_pp.write_zgrid(obj);
            end
        end
        
        function write_zgrid(obj)
            zgrid = fopen(['zgrid.inp.' obj.expnr], 'w');
            fprintf(zgrid, '%12s\n', '#     z-grid');
            fprintf(zgrid, '%12s\n', '#           ');
            fprintf(zgrid, '%-20.15f\n', obj.zf);
            fclose(zgrid);
            disp(['... written zgrid.inp.' obj.expnr]) 
        end
        
        function generate_lscale(obj, lwritefile)
            if (obj.lmassflowr + obj.lprofforc + obj.lcoriol + obj.ldp) > 1
                error('More than one forcing specified')
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
            
            if lwritefile
                da_pp.write_lscale(obj);
            end
        end
        
        function write_lscale(obj)
            lscale = fopen(['lscale.inp.' obj.expnr], 'w');
            fprintf(lscale, '%-12s\n', '# SDBL flow');
            fprintf(lscale, '%-60s\n', '# z uq vq pqx pqy wfls dqtdxls dqtdyls dqtdtls dthlrad');
            fprintf(lscale, '%-20.15f %-12.6f %-12.6f %-12.6f %-12.6f %-15.9f %-12.6f %-12.6f %-12.6f %-17.12f\n', obj.ls');
            fclose(lscale);
            disp(['... written lscale.inp.' obj.expnr]) 
        end
        
        function generate_prof(obj, lwritefile)
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
            
            if lwritefile
                da_pp.write_prof(obj)
            end
        end
        
        function write_prof(obj)
            prof = fopen(['prof.inp.' obj.expnr], 'w');
            fprintf(prof, '%-12s\n', '# SDBL flow');
            fprintf(prof, '%-60s\n', '# z thl qt u v tke');
            fprintf(prof, '%-20.15f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f\n', obj.pr');
            fclose(prof);
            disp(['... written prof.inp.' obj.expnr])
        end
        
        function generate_scalar(obj, lwritefile)
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
            
            if lwritefile
                if obj.nsv > 0
                    da_pp.write_scalar(obj)
                end
            end
        end
        
        function write_scalar(obj)
            scalar = fopen(['scalar.inp.' obj.expnr], 'w');
            fprintf(scalar, '%-12s\n', '# SDBL flow');
            fprintf(scalar, '%-60s\n', '# z sca1 sca2 sca3 sca4');
            fprintf(scalar, '%-20.15f %-14.10f %-14.10f %-14.10f %-14.10f\n', obj.sc');
            fclose(scalar);
            disp(['... written scalar.inp.' obj.expnr]) 
        end
        
        function generate_blocks(obj, lwritefile)
            %aspectratio = r.zh(r.blockheight+1)/(r.xh(r.canyonwidth+1)-r.xh(1));
            aspectratio = obj.zh(obj.blockheight + 1) / (obj.xh(obj.canyonwidth + 1) - obj.xh(1));
            %nrows = ie/(r.blockwidth+r.canyonwidth);
            da_pp.addvar(obj, 'nrows', obj.imax / (obj.blockwidth + obj.canyonwidth));
            
            if obj.lflat
    
            elseif ceil(obj.nrows) ~= floor(obj.nrows)
                %l = 0:ie/2;
                l = 0:obj.imax / 2;
                %ind = rem(ie/2,l)==0; %// logical index that tells if remainder is zero or not
                ind = rem(obj.imax / 2, l) == 0;
                %err = ([l(ind); (ie/2)./l(ind)]);
                err = ([l(ind); (obj.imax / 2) ./ l(ind)]);
                disp('Block system does not fit grid')
                disp(['sum widths to: ' num2str(err(1,:))])
                %disp(['Current width: ' num2str(r.blockwidth + r.canyonwidth)])
                disp(['Current width: ' num2str(obj.blockwidth + obj.canyonwidth)])
                error('Incorrect block system')
            else
                disp('Successful block network')
                disp(['aspect ratio: ' num2str(aspectratio)])
            end 
            
            if obj.lflat
    
            da_pp.addvar(obj, 'bl', []);

            elseif obj.lcastro
                %nrows = ie/(r.blockwidth*2);
                obj.nrows = obj.imax / (obj.blockwidth * 2);
                %ncolumns = je/(r.blockwidth*2);
                da_pp.addvar(obj, 'ncolumns', obj.jtot / (obj.blockwidth * 2));
                %bl = zeros(obj.nrows * obj.ncolumns + obj.nrows / 2, 13);
                da_pp.addvar(obj, 'bl', obj.nrows * obj.ncolumns + obj.nrows / 2, 13);
                obj.bl(:,5) = 0;
                obj.bl(:,6) = obj.blockheight - 1;
                
                for n = 1:obj.nrows
                    for nn = 0:obj.ncolumns
                        if mod(n,2) == 0
                            if nn == 0
                                %bl(length(nonzeros(bl(:,3)))+1,3) = jb;
                                obj.bl(length(nonzeros(obj.bl(:,3))) + 1, 3) = obj.jmin;
                                %bl(length(nonzeros(bl(:,4))) + 1, 4) = r.blockwidth/2;
                                obj.bl(length(nonzeros(obj.bl(:,4))) + 1, 4) = obj.blockwidth / 2;
                            elseif nn == obj.ncolumns
                                %bl(length(nonzeros(bl(:,3)))+1,3) = je-r.blockwidth/2+1;
                                obj.bl(length(nonzeros(obj.bl(:,3))) + 1, 3) = obj.jtot - obj.blockwidth / 2 + 1;
                                %bl(length(nonzeros(bl(:,4)))+1,4) = je;
                                obj.bl(length(nonzeros(obj.bl(:,4))) + 1, 4) = obj.jtot;
                            else
                                %bl(length(nonzeros(bl(:,3)))+1,3) = jb + nn * r.blockwidth*2 - r.blockwidth/2;
                                obj.bl(length(nonzeros(obj.bl(:,3))) + 1, 3) = obj.jmin + nn * obj.blockwidth * 2 - obj.blockwidth / 2;
                                %bl(length(nonzeros(bl(:,4)))+1,4) = bl(length(nonzeros(bl(:,4)))+1,3) + r.blockwidth - 1;
                                obj.bl(length(nonzeros(obj.bl(:,4))) + 1, 4) = bl(length(nonzeros(obj.bl(:,4))) + 1, 3) + obj.blockwidth - 1;
                            end
                            %bl(length(nonzeros(bl(:,1)))+1,1) = -0.5*r.blockwidth + (2*n-1) * r.blockwidth + 1;
                            obj.bl(length(nonzeros(obj.bl(:,1))) + 1, 1) = - 0.5 * obj.blockwidth + (2 * n - 1) * obj.blockwidth + 1;
                            %bl(length(nonzeros(bl(:,2)))+1,2) = bl(length(nonzeros(bl(:,2)))+1,1) + r.blockwidth - 1;
                            obj.bl(length(nonzeros(obj.bl(:,2))) + 1, 2) = bl(length(nonzeros(obj.bl(:,2)))+1,1) + obj.blockwidth - 1;
                        end
                    end
                    for nn = 0:obj.ncolumns - 1
                        if mod(n,2) ~= 0
                            %bl(length(nonzeros(bl(:,3)))+1,3) = jb + r.blockwidth + nn * r.blockwidth*2 - r.blockwidth/2;
                            obj.bl(length(nonzeros(obj.bl(:,3))) + 1, 3) = obj.jmin + obj.blockwidth + nn * obj.blockwidth * 2 - obj.blockwidth / 2;
                            %bl(length(nonzeros(bl(:,4)))+1,4) = bl(length(nonzeros(bl(:,4)))+1,3) + r.blockwidth - 1;
                            obj.bl(length(nonzeros(obj.bl(:,4))) + 1, 4) = bl(length(nonzeros(obj.bl(:,4))) + 1,3) + obj.blockwidth - 1;
                            %bl(length(nonzeros(bl(:,1)))+1,1) = -0.5*r.blockwidth + (2*n-1) * r.blockwidth + 1;
                            obj.bl(length(nonzeros(obj.bl(:,1))) + 1, 1) = - 0.5 * obj.blockwidth + (2 * n - 1) * obj.blockwidth + 1;
                            %bl(length(nonzeros(bl(:,2)))+1,2) = bl(length(nonzeros(bl(:,2)))+1,1) + r.blockwidth - 1;
                            obj.bl(length(nonzeros(obj.bl(:,2))) + 1, 2) = bl(length(nonzeros(bl(:,2))) +1 , 1) + obj.blockwidth - 1;
                        end

                    end
                end

            elseif obj.lcube
                %nrows = ie/(r.blockwidth*2);
                obj.nrows = obj.imax / (obj.blockwidth * 2);
                %ncolumns = je/(r.blockwidth*2);
                da_pp.addvar(obj, 'ncolumns', obj.jtot / (obj.blockwidth * 2));
                %bl = zeros(nrows*ncolumns,11);
                da_pp.addvar(obj, 'bl', zeros(obj.nrows * obj.ncolumns, 13));
                for n = 1 : obj.nrows
                    for nn = 0 : obj.ncolumns - 1
                        %bl((n - 1) * ncolumns + nn + 1, 1) = -0.5 * r.blockwidth + (2 * n - 1) * r.blockwidth + 1;
                        obj.bl((n - 1) * obj.ncolumns + nn + 1, 1) = -0.5 * obj.blockwidth + (2 * n - 1) * obj.blockwidth + 1;
                        %bl((n-1)*ncolumns+nn+1,2) = bl((n-1)*ncolumns+nn+1,1) + r.blockwidth - 1;
                        obj.bl((n-1) * obj.ncolumns + nn + 1, 2) = obj.bl((n - 1) * obj.ncolumns + nn + 1, 1) + obj.blockwidth - 1;
                        %bl((n-1)*ncolumns+nn+1,5) = 0;
                        obj.bl((n - 1) * obj.ncolumns + nn + 1, 5) = 0;
                        %bl((n-1)*ncolumns+nn+1,6) = ceil(r.blockwidth*(h/ie)/(hz/ke));
                        obj.bl((n - 1) * obj.ncolumns + nn + 1, 6) = ceil(obj.blockwidth * (obj.xsize / obj.imax) / (obj.zsize / obj.kmax));
                        %bl((n-1)*ncolumns+nn+1,3) = jb + r.blockwidth/2 + nn * r.blockwidth*2;
                        obj.bl((n - 1) * obj.ncolumns + nn + 1, 3) = obj.jmin + obj.blockwidth / 2 + nn * obj.blockwidth * 2;
                        %bl((n-1)*ncolumns+nn+1,4) = bl((n-1)*ncolumns+nn+1,3) + r.blockwidth - 1;
                        obj.bl((n - 1) * obj.ncolumns + nn + 1, 4) = obj.bl((n - 1) * obj.ncolumns + nn + 1, 3) + obj.blockwidth - 1;
                    end
                end

            elseif obj.lblocks
                %bl = zeros(nrows, 13);
                da_pp.addvar(obj, 'bl', zeros(nrows, 13));
                %bl(1:nrows,1) = (r.canyonwidth/2+1:r.canyonwidth+r.blockwidth:ie-r.canyonwidth/2)';
                obj.bl(1:obj.nrows, 1) = (obj.canyonwidth / 2 + 1 : obj.canyonwidth + obj.blockwidth : obj.imax - obj.canyonwidth / 2)';
                %bl(1:nrows,2) = bl(1:nrows,1) + r.blockwidth - 1;
                obj.bl(1:obj.nrows, 2) = obj.bl(1:obj.nrows,1) + obj.blockwidth - 1;
                %bl(:,3) = jb;
                obj.bl(:,3) = obj.jmin;
                %bl(:,4) = je;
                obj.bl(:,4) = obj.jmin;
                obj.bl(:,5) = 0;
                %bl(1:nrows,6) = r.blockheight-1;
                obj.bl(1:obj.nrows, 6) = obj.blockheight - 1; 
            end 
            
            da_pp.addvar(obj, 'blocks', zeros(1,11));
            da_pp.addvar(obj, 'facets', zeros(1,4));
            bl2blocks_temp
            makeblocks
            block2fac
            if obj.lEB
                addboundingwalls
            else
                nwallfcts = 0;
            end
            createfloors
            
            if lwritefile
                da_pp.write_blocks(obj);
            end
            
            if lwritefile
                da_pp.write_facets(obj);
            end
            
            if ~obj.lEB
                if lwritefile
                    da_pp.write_Tfac(obj); % switch in write Tfac
                end
            end
            
        end
        
        function write_blocks(obj)
            blocks = fopen( ['blocks.inp.' obj.expnr], 'w');
            fprintf(blocks, '# Block data\n');
            fprintf(blocks, '#  il\t   iu\t   jl\t   ju\t   kl\t   ku\t dtop\t dwest\t deast\t dnor\t dsou\n');
            fprintf(blocks, '%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\n', obj.blocks');
            fclose(blocks);
            disp(['... written blocks_old.inp.' obj.expnr])
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
            disp(['... written trees.inp.' obj.expnr]) 
        end
        
        function generate_purifs(obj, lwritefile)
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
            disp(['... written purifs.inp.' obj.expnr]) 
        end
        
        function plot_domain(obj)
            figure
            view(52,23)
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
        
        
    end
    
end