% udbase base class for uDALES simulations

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

% Aug 2024, Maarten van Reeuwijk, Jingzi Huang. first version
% Dec 2024, Maarten van Reeuwijk, Chris Wilson. added facet functionality.

% Copyright (C) 2016-2024 the uDALES Team.

classdef udbase < dynamicprops
    %  Base class for uDALES simulations

    properties (SetAccess = protected)
        expnr;                   % the simulation ID
        geom;                    % geometry instance of udgeom

        % grid
        xm;                      % x-location for u location
        ym;                      % y-location for v location
        zm;                      % z-location for w location
        xt;                      % x-location for v,w,c location
        yt;                      % y-location for u,w,c location
        zt;                      % z-location for u,v,c location
        dx;                      % x-increment
        dy;                      % y-increment
        dzm;                     % z-increment for w locations
        dzt;                     % z-increment for u,v,c locations
        
        % 3D logical arrays that identify the cells that are solid. Note
        % that there are four arrays since the underlying grid is staggered.
        Su;                      % 3d array for variables at u-locations
        Sv;                      % 3d array for variables at v-locations
        Sw;                      % 3d array for variables at w-locations
        Sc;                      % 3d array for variables at c-locations

        % facet information
        facs;                    % structure with facet info
        factypes;                % structure with facet type info
        facsec;                  % structure with fecet section info
    end
    properties (Hidden = true, SetAccess = protected)
        path;                    % Path to simulations.
        cpath;                   % Current path.
        buf;                     % All purpose buffer.

        % standard filenames
        fnamoptions      = 'namoptions'
        fprof            = 'prof.inp'
        fxytdump         = 'xytdump';
        ftdump           = 'tdump';
        ffielddump       = 'fielddump';
        fislicedump      = 'islicedump';
        fjslicedump      = 'jslicedump';
        fkslicedump      = 'kslicedump';
        fsolid           = 'solid';
        ffacEB           = 'facEB';
        ffacT            = 'facT';
        ffac             = 'fac';
        ffacets          = 'facets.inp';
        ffactypes        = 'factypes.inp';
        ffacetarea       = 'facetarea.inp';
        ffluid_boundary  = 'fluid_boundary';
        ffacet_sections  = 'facet_sections';

        % logicals for whether input files are present
        lfprof           = true;
        lfsolid          = true;
        lffacEB          = true;
        lffacT           = true;
        lffac            = true;
        lffacets         = true;
        lffactypes       = true;
        lffacetarea      = true;
        lffluid_boundary = true;
        lffacet_sections = true;
        lfgeom           = true; % stl geometry
    end

    methods
        function obj = udbase(expnr, varargin)
            % Class constructor.
            %
            %
            % udbase(expnr, dapath, load_preprocdata) 
            %    expnr:                       experiment number
            %    dapath (optional):           path to the experiment
            %
            % Example:
            %   obj = udbase(expnr, '/path/experiments/../100');

            obj.expnr = expnr;

            % store current and simulation directory
            obj.cpath = pwd;

            dapath = obj.cpath;
            if (nargin>=2)
               dapath = varargin{1};
            end

            % go to the simulation directory  
            cd(dapath);
            obj.path = pwd;

            % ------ %
            % open and read the namoptions file
            
            TOKENS = '(.*)\=(.*)';
            DROP = '(\s*|\[|\]|''|''|;)';
            WHITE = '\s*';
            VOID = '';

            expstr = num2str(expnr, '%3.3d');
            fid = fopen([obj.fnamoptions, '.', expstr]);

            if (fid == -1) % file not present
                disp([obj.fnamoptions, '.',  expstr, ' not found. Exiting ...'])
                return
            end

            obj.expnr = expstr;

            line = fgets(fid);

            % load each line of namoptions file into object
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
            
            % ------- %
            % load grid

            obj.dx = obj.xlen/obj.itot;
            obj.dy = obj.ylen/obj.itot;
            obj.xm = (0:obj.itot-1)' * obj.dx;
            obj.ym = (0:obj.jtot-1)' * obj.dy;
            obj.xt = obj.xm + 0.5 * obj.dx;
            obj.yt = obj.ym + 0.5 * obj.dy;

            % load prof.inp if present
            try
                data = obj.load_prof();
            catch
                obj.lfprof = false;
            end

            if (obj.lfprof) 
                obj.zt  = data(:, 1); % cell centred coordinates provided.
                obj.zm  = [0; 0.5 * (obj.zt(1:end-1)+obj.zt(2:end)); obj.zsize];
                obj.dzt = diff(obj.zm);
                obj.dzm = [2*obj.zt(1); diff(obj.zt)];

                % make sure zm arrays have the same length as zt arrays
                % assumes that the first ghost point is outputted
                obj.zm  = obj.zm(1:end-1); 
            else % file not present
                warning([obj.fprof, '.', expstr, ' not found. Assuming equidistant grid.'])
                dz = obj.zsize/obj.ktot;

                obj.zm = (0:obj.ktot-1)' * dz;
                obj.zt = obj.zm + 0.5 * dz;

                obj.dzm = dz * ones(obj.ktot, 1);
                obj.dzt = dz * ones(obj.ktot, 1);
            end

            % ------ %
            % read the stl file of the geometry if present

            if isprop(obj, 'stl_file')  && exist(obj.stl_file, 'file')
                geom = udgeom.udgeom(obj.path);
                geom.load(obj.stl_file);
                obj.geom = geom;
            else
                obj.lfgeom = false;
                warning('STL file not found, implying a flat surface. obj.geom is empty')
                geom = [];
            end
            obj.geom = geom;

            % ------ %
            % load 3D arrays that indicate the location of solid regions.
            % [note that there are four arrays since the domain is
            % staggered].

            try
                obj.Su = obj.load_solid([obj.fsolid, '_u.txt'], obj.itot, obj.jtot, obj.ktot);
                obj.Sv = obj.load_solid([obj.fsolid, '_v.txt'], obj.itot, obj.jtot, obj.ktot);
                obj.Sw = obj.load_solid([obj.fsolid, '_w.txt'], obj.itot, obj.jtot, obj.ktot);
                obj.Sc = obj.load_solid([obj.fsolid, '_c.txt'], obj.itot, obj.jtot, obj.ktot);
            catch
                warning(['One or more', obj.fsolid, '_(u,v,w,c).', expstr, ' files not found.']);
                obj.lfsolid = false;
            end

            % ------ %
            % Read facet data

            try % this file is only available when SEB is on.
                facareas = obj.load_facetarea();
                facs.area = facareas;
            catch
                obj.lffacetarea = false;
            end

            try
                facets = obj.load_facets();
                facettypeids = facets(:,1);
                facs.typeid = facettypeids;
            catch
                obj.lffacets = false;
            end

            % ------ %
            % Read facet type data
            try
                factypedata = obj.load_factypes();
                factypes.id = factypedata(:,1);

                % add names to facet type information
                key = [0, -1, -101, 1, 2, 3, 4, 11, 12];
                names = ["Default dummy" "Asphalt floor" ...
                         "Concrete bounding wall" "Concrete" ...
                         "Brick" "Stone" "Painted wood" ...
                         "Green 1" "Green 2"];
                % dict = dictionary(key, names); % dict only available >= Matlab2022b
                dict = containers.Map(key, names);
                newkey = setdiff(unique(factypes.id), key);
                for i = 1:length(newkey)
                    dict(newkey(i)) = "Custom walltype";
                end
                %factypes.name = dict(factypes.id);% dict only available >= Matlab2022b
                factypes.name = string(arrayfun(@(x) dict(x), factypes.id, 'UniformOutput', false));
                factypes.lGR = factypedata(:, 2);
                factypes.z0 = factypedata(:, 3);
                factypes.z0h = factypedata(:, 4);
                factypes.al = factypedata(:, 5);
                factypes.em = factypedata(:, 6);
                factypes.d = factypedata(:, 7:7+obj.nfaclyrs-1);
                factypes.C = factypedata(:, 7+obj.nfaclyrs:7+2*obj.nfaclyrs-1);
                factypes.lam = factypedata(:, 7+2*obj.nfaclyrs:7+3*obj.nfaclyrs);
            catch
                obj.lffactypes = false;
            end

            % ------ %
            % load information about the facet sections
            try
                facsec.u = obj.load_facsec('u');
                facsec.v = obj.load_facsec('v');
                facsec.w = obj.load_facsec('w');
                facsec.c = obj.load_facsec('c');

                % Add facet variables to the instance
                obj.facs = facs;
                obj.factypes = factypes;
                obj.facsec = facsec;
            catch
                obj.lffacet_sections = false;
            end

            % ------ %
            % store colors for use in plot_fac_types
            fig = figure;
            ax = gca;
            colors = ax.ColorOrder;
            close(fig);

            h = addprop(obj, 'colors');
            h.Hidden = true;
            h.SetAccess = 'protected';
            obj.colors = colors;
            
            obj.gohome()
        end
    
        % ------------------------------------------------------------- %

        function data = load_stat_xyt(obj, varargin)
            % A method to retrieve plane- and time-averaged 1D statistics
            % information from the xytdump file.  
            %
            % load_stat_xyt(OBJ) displays the variables in the xytdump file
            %
            % load_stat_xyt(OBJ, svar) retrieves a variable from the xytdump file
            %
            % Example (view contents of output):
            %   obj = udbase(expnr);
            %   obj.load_stat_xyt();

            obj.gopath();

            fname = [obj.fxytdump, '.', obj.expnr, '.nc'];
            ffullname = fullfile(obj.path, fname);
            data = udbase.load_ncdata(ffullname, varargin);

            obj.gohome();
        end

        % ------------------------------------------------------------- %

        function data = load_stat_t(obj, varargin)
            % A method to retrieve time-averaged statistics from the tdump file 
            %
            % load_stat_t(OBJ) displays the variables in the tdump file
            %
            % load_stat_t(OBJ, svar) retrieves a variable from the tdump file
            %
            % Example (view contents of output):
            %   obj = udbase(expnr);
            %   obj.load_stat_t();

            obj.gopath();

            fname = [obj.ftdump, '.', obj.expnr, '.nc'];
            ffullname = fullfile(obj.path, fname);
            data = udbase.load_ncdata(ffullname, varargin);

            obj.gohome();
        end

        % ------------------------------------------------------------- %

        function data = load_field(obj, varargin)
            %  A method to retrieve instantaneous 3D data from from the fielddump file.
            %
            %  load_field(OBJ) displays the variables in the fielddump file
            %
            %  load_field(OBJ, svar) retrieves variable svar from the fielddump file
            %
            %  Example (view contents of output):
            %    obj = udbase(expnr);
            %    obj.load_field();

            obj.gopath()

            fname = [obj.ffielddump, '.', obj.expnr, '.nc'];
            ffullname = fullfile(obj.path, fname);
            data = udbase.load_ncdata(ffullname, varargin);

            obj.gohome();
        end

        % ------------------------------------------------------------- %

        function data = load_slice(obj, option, varargin)
            % A method to retrieve instantaneous 2D slices from from the slicedump file.
            %
            % load_slice(obj) displays information, option chooses plane.
            %
            % load_slice(obj, ...) 
            %
            % Example (view contents of output of horizontal slice):
            %    obj = udbase(expnr);
            %    obj.load_slice('k');

            obj.gopath();

            switch option
                case 'i'
                    fname = [obj.fislicedump, '.', obj.expnr, '.nc'];
                case 'j'
                    fname = [obj.fjslicedump, '.', obj.expnr, '.nc'];
                case 'k'
                    fname = [obj.fkslicedump, '.', obj.expnr, '.nc'];
                otherwise
                    error('Invalid option. Use ''i'', ''j'', or ''k''.');
            end
            ffullname = fullfile(obj.path, fname);
            data = udbase.load_ncdata(ffullname, varargin);

            obj.gohome();
        end

             % ------------------------------------------------------------- %

        function data = load_fac_eb(obj, varargin)
            % A method to retrieve facet data of a surface energy balance
            % term from the facEB file.   
            %
            % load_fac_eb(OBJ) displays the variables in the file
            %
            % load_fac_eb(OBJ, svar) retrieves a variable from the file
            %
            % Example (view contents of output):
            %   obj = udbase(expnr);
            %   obj.load_fac_eb();

            obj.gopath();
            fname = [obj.ffacEB, '.', obj.expnr, '.nc'];
            ffullname = fullfile(obj.path, fname);
            data = udbase.load_ncdata(ffullname, varargin);
            obj.gohome();
        end

        % ------------------------------------------------------------- %

        function data = load_fac_temperature(obj, varargin)
            % A method to retrieve temperature and temperature gradient
            % data from the facT file.   
            %
            % load_fac_temperature(OBJ) displays the variables in the file
            %
            % load_fac_temperature(OBJ, svar) retrieves a variable from the file
            %
            % Example (view contents of output):
            %   obj = udbase(expnr);
            %   obj.load_fac_temperature();

            obj.gopath();
            fname = [obj.ffacT, '.', obj.expnr, '.nc'];
            ffullname = fullfile(obj.path, fname);
            data = udbase.load_ncdata(ffullname, varargin);
            obj.gohome();
        end

        % ------------------------------------------------------------- %
    
        function data = load_fac_momentum(obj, varargin)
            % A method to retrieve facet data for pressure and shear
            % information from the fac file.  
            %
            % load_fac_momentum(OBJ) displays the variables in the file
            %
            % load_fac_momentum(OBJ, svar) retrieves a variable from the file
            %
            % Example (view contents of output):
            %   obj = udbase(expnr);
            %   obj.load_fac_momentum();

            obj.gopath();
            fname = [obj.ffac, '.', obj.expnr, '.nc'];
            ffullname = fullfile(obj.path, fname);
            data = udbase.load_ncdata(ffullname, varargin);
            obj.gohome();
        end
        
        % ------------------------------------------------------------- %

        function data = assign_prop_to_fac(obj, strprop)
            % Method for assigning properties of a material (stored in
            % factypes) to each facet for visualisation or calculation..   
            %
            % assign_prop_to_fac(OBJ, strprop) assigns the property strprop
            %                                  to the appropriate facet.
            %
            % Example (assign albedo to each facet):
            %   obj = udbase(expnr);
            %   al = obj.assign_prop_to_fac('al');

            % Function only works when required data has been loaded.
            if (~obj.lffactypes)
                error(['This method requires the file ', obj.ffactypes, '.', obj.expnr])
            end

            % Determine the dimensions of obj.factypes.(strprop)
            values = obj.factypes.(strprop);
            dims = size(values);
            [found, idx] = ismember(obj.facs.typeid, obj.factypes.id);

            if isvector(values)
               data = NaN(size(obj.facs.typeid)); % Initialize with NaN for unmatched cases
               data(found) = values(idx(found));
            elseif ismatrix(values)
               % Case 2: values is a 2D array
               data = NaN(size(obj.facs.typeid, 1), dims(2)); % Initialize for 2D array
               data(found, :) = values(idx(found), :); % Assign corresponding rows
            else
                error('Unsupported dimensions for the property array.');
            end
        end

        % ------------------------------------------------------------- %

        function fld = convert_fac_to_field(obj, var, facsec, dz)
            % Method for converting a facet variable to a density in a 3D
            % field. 
            %
            % convert_fac_to_field(OBJ, var) converts the facet variable var  
            %                                to a 3D field density.   
            %
            % Example:
            %   fld = obj.convert_fac_to_field(var);

             % Function only works when required data has been loaded.
            if (~obj.lffacet_sections)
                error(['This method requires the files ', ...
                       obj.ffacet_sections, '_(u,v,w,c).', obj.expnr, ' and ', ...
                       obj.ffluid_boundary, '_(u,v,w,c).', obj.expnr, ...
                       ' to execute.']);
            end

            fld = zeros(obj.itot, obj.jtot, obj.ktot, 'single');
            facsec_ind = sub2ind(size(fld), ...
                                 facsec.locs(:, 1), ...
                                 facsec.locs(:, 2), ...
                                 facsec.locs(:, 3));

            % loop over all facet sections and create density field for var
            for m = 1:length(facsec.area) 
                facid = facsec.facid(m);
                fld(facsec_ind(m)) = fld(facsec_ind(m)) ...
                                   + var(facid) .* facsec.area(m) ...
                                   / (obj.dx * obj.dy * dz(facsec.locs(m, 3)));
            end
        end

        % ------------------------------------------------------------- %

        function seb = load_seb(obj)
            % A method to retrieve all surface energy balance terms on each
            % of the facets as a function of time. 
            %
            % load_seb(OBJ) loads the terms in the surface energy balance.
            %
            % Example:
            %   obj = udbase(expnr);
            %   obj.load_seb();

            t = load_fac_eb(obj, 't');
            K = load_fac_eb(obj, 'netsw');
            Lin = load_fac_eb(obj, 'LWin');
            Lout = load_fac_eb(obj, 'LWout');
            H = load_fac_eb(obj, 'hf');
            E = load_fac_eb(obj, 'ef');

            % load data required to calculate G
            T = load_fac_temperature(obj, 'T');
            dTdz = load_fac_temperature(obj, 'dTdz');
            lam = obj.assign_prop_to_fac('lam');

            % Manipulate data
            L = Lin-Lout;
            G = -lam(:, 1).*squeeze(dTdz(:,1,:)); % calculate ground heat flux
            Tsurf = squeeze(T(:,1,:));
            seb.Kstar = K;
            seb.Lstar = L;
            seb.Lin = Lin;
            seb.Lout = Lout;
            seb.H = -H;
            seb.E = -E;
            seb.G = G;
            seb.Tsurf = Tsurf;
            seb.t = t;
        end 
        % ------------------------------------------------------------- %

        function seb_av = area_average_seb(obj, seb)
            % A method for calculating the area-averaged surface energy
            % balance.
            %
            % area_average_seb(OBJ, seb) averages the variables in seb over the
            % area.
            %
            % Example:
            %   obj = udbase(expnr);
            %   seb = obj.load_seb()
            %   seb_av = obj.area_average_seb(seb);

            Kav = obj.area_average_fac(seb.Kstar);
            Lav = obj.area_average_fac(seb.Lstar);
            Linav = obj.area_average_fac(seb.Lin);
            Loutav = obj.area_average_fac(seb.Lout);
            Hav = obj.area_average_fac(seb.H);
            Eav = obj.area_average_fac(seb.E);
            Gav = obj.area_average_fac(seb.G);
    
            seb_av.Kstar = Kav;
            seb_av.Lstar = Lav;
            seb_av.Lin = Linav;
            seb_av.Lout = Loutav;
            seb_av.H = Hav;
            seb_av.E = Eav;
            seb_av.G = Gav;
            seb_av.t = seb.t;
        end

        % ------------------------------------------------------------- %

        function plot_fac(obj, var)
            % A method for plotting a facet variable var as a 3D surface
            %
            % plot_fac(OBJ, var) plots variable var
            %
            % Example (plot net shortwave radiation):
            %   obj.plot_fac(K); 

            % Function only works when required data has been loaded.
            if (~obj.lfgeom)
                error('This method requires a geometry (STL) file.');
            end 

            cons = obj.geom.stl.ConnectivityList;
            points = obj.geom.stl.Points;
            patch('Faces', cons, 'Vertices', points, 'FaceVertexCData',var,'FaceColor','flat', 'EdgeColor', 'None')
        end 

        % ------------------------------------------------------------- %

        function plot_fac_type(obj)
            % A method for plotting the different surface types used in a
            % geometry.
            %
            % plot_fac_types(OBJ) plots the surface types 
            %
            % Example:
            %   obj.plot_fac_type();

            % Function only works when required data has been loaded.
            if (~obj.lfgeom)
                error('This method requires a geometry (STL) file.');
            end  
            if (~obj.lffacets)
                error(['This method requires the file ', ...
                       obj.ffacets, '.', obj.expnr])
            end            
            if (~obj.lffactypes)
                error(['This method requires the file ', ...
                       obj.ffactypes, '.', obj.expnr])
            end
     
            facids = obj.facs.typeid;
            labs = obj.factypes.name;
            typeids = obj.factypes.id;
            unique_ids = unique(facids);
            coloursneeded = length(unique(facids));
            colors = obj.colors;
            
            if length(colors(:,1))<coloursneeded
                disp('Too many surfaces requested to colour uniquely.')
                return
            end 

            cs = [];
            for i = 1:length(facids)
                typeid = facids(i);
                colorind = typeid==unique_ids;
                cs = [cs;colors(colorind,:)];
            end 
            cons = obj.geom.stl.ConnectivityList;
            ps = obj.geom.stl.Points;
            hold on
            labels = [];
            for m = 1:length(unique_ids)
                id = unique_ids(m);
                msel = facids==id;
                n = find(id==typeids, 1);
                labels = [labels;{labs(n)}];
                tcons = cons(msel,:);
                c = colors(m,:);
                patch('Faces', tcons, 'Vertices', ps, 'FaceVertexCData',c, ...
                      'FaceColor','flat', 'EdgeColor', 'None')
            end
            legend(labels)
        end

        % ------------------------------------------------------------- %
        
        function res = calculate_frontal_properties(obj)
            % A method to calculate the skyline, frontal areas and blockage
            % ratios in the x- and y-direction. 
            %
            % res = calculate_frontal_properties(OBJ) executes the method
            %                               and returns a structure with
            %                               the skyline, frontal areas and
            %                               blockage ratios in x- and
            %                               y-direction. 
            % Example:
            %   res = obj.calculate_frontal_properties();

            % Function only works when required data has been loaded.
            if (~obj.lfgeom)
                error('This method requires a geometry (STL) file.');
            end 
            if (~obj.lffacet_sections)
                error(['This method requires the files ', ...
                       obj.ffacet_sections, '_(u,v,w,c).', obj.expnr, ' and ', ...
                       obj.ffluid_boundary, '_(u,v,w,c).', obj.expnr, ...
                       ' to execute.']);
            end

            norms = obj.geom.stl.faceNormal;

            % create surface quantity for projected area in x- and y-direction
            phix = -min(dot(norms, repmat([1, 0, 0], size(norms, 1), 1), 2),0);
            phiy = -min(dot(norms, repmat([0, 1, 0], size(norms, 1), 1), 2),0);

            % convert to a density field
            rhoLx = obj.convert_fac_to_field(phix,obj.facsec.c,obj.dzt);
            rhoLy = obj.convert_fac_to_field(phiy,obj.facsec.c,obj.dzt);

            % Calculate indicator functions for blockage and plot them
            Ibx = double(squeeze(sum(rhoLx,1))>0);
            Iby = double(squeeze(sum(rhoLy,2))>0);

            % Integrate to get the frontal areas and blockage ratio:
            Afx = 0; brx = 0;
            Afy = 0; bry = 0;
            for k = 1:length(obj.dzt)
                Afx = Afx + sum(sum(rhoLx(:,:,k)*obj.dx*obj.dy*obj.dzt(k)));
                Afy = Afy + sum(sum(rhoLy(:,:,k)*obj.dx*obj.dy*obj.dzt(k)));
                brx = brx + sum(Ibx(:,k))*obj.dy*obj.dzt(k);
                bry = bry + sum(Iby(:,k))*obj.dx*obj.dzt(k);
            end
            brx = brx / (obj.ylen * obj.zsize);
            bry = bry / (obj.xlen * obj.zsize);

            res.skylinex = Ibx;
            res.skyliney = Iby;
            res.Afx = Afx;
            res.Afy = Afy;
            res.brx = brx;
            res.bry = bry;

            disp(['x-direction: frontal area = ', num2str(Afx, '%8.1f'), ...
                  ' m2, blockage ratio = ', num2str(brx, '%8.3f')])
            disp(['y-direction: frontal area = ', num2str(Afy, '%8.1f'), ...
                  ' m2, blockage ratio = ', num2str(bry, '%8.3f')])
        end

        % ------------------------------------------------------------- %

        function av = area_average_fac(obj, var, varargin)
            % A method for area-averaging a facet quantity, either over all
            % facets or over a selection.  
            %
            % area_average_fac(OBJ, var) area-averages variable var over
            % all facets.
            %
            % area_average_fac(OBJ, var, sel) area-averages variable var
            % over the facet indices in sel
            %
            % Example:
            %   K_av = obj.area_average_fac(K);   

            % Function only works when required data has been loaded.
            if (~obj.lffacetarea)
                error(['This method requires the file ', obj.ffacetarea, '.', obj.expnr])
            end

            if isempty(varargin)   
                sel = 1:length(obj.facs.area);
            else
                sel = varargin{1};
            end 

            % the statement below works rather remarkably -- needs testing.
            tot_var = sum(var(sel, :, :, :).*obj.facs.area(sel));
            tot_area = sum(obj.facs.area(sel));
            av = squeeze(tot_var/tot_area);
        end 

        % ------------------------------------------------------------- %

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

        % ------------------------------------------------------------- %

        function gopath(obj)
            % Go to simulation path.
            % gopath(obj)

            cd(obj.path);
        end

        % ------------------------------------------------------------- %

        function gohome(obj)
            % Go to work path.
            % gohome(obj)

            cd(obj.cpath);
        end

        % ------------------------------------------------------------- %

        function chcpath(obj, newpath)

            % Change work path.
            % chcpath(obj, newpath)

            here = pwd;
            cd(newpath)
            obj.cpath = pwd;
            cd(here);
        end
    end

    % ----------------------------------------------------------------- %

    methods (Access = protected)
        function data = load_prof(obj)
            % A method to retrieve information from prof.inp file

            fname = [obj.fprof '.' obj.expnr];
            ffullname = fullfile(obj.path, fname);
            data = readmatrix(ffullname,'Range',[3,1], 'FileType', 'text');
        end 

        % ------------------------------------------------------------- %
        
        function data = load_factypes(obj)
            % A method to retrieve information from factypes.inp file

            fname = [obj.ffactypes '.' obj.expnr];
            ffullname = fullfile(obj.path, fname);
            data = readmatrix(ffullname,'Range',[4,1], 'FileType', 'text');
        end 

        % ------------------------------------------------------------- %

        function data = load_facets(obj)
             % A method to retrieve information from facets.inp file

            fname = [obj.ffacets '.' obj.expnr];
            ffullname = fullfile(obj.path, fname);
            data = readmatrix(ffullname,'Range',[2,1], 'FileType', 'text');
        end 

        % ------------------------------------------------------------- %

        function data = load_facetarea(obj)
            % A method to retrieve information from facetarea.inp file

            fname = [obj.ffacetarea '.' obj.expnr];
            ffullname = fullfile(obj.path, fname);
            data = readmatrix(ffullname,'Range',[2,1], 'FileType', 'text');
        end 

        % ------------------------------------------------------------- %

        function data = load_facsec(obj, strvar)
            % A method to retrieve information about facet section for
            % variable strvar
            
            % load facet section data for variable strvar
            fname = [obj.path '/' obj.ffacet_sections '_' strvar '.txt'];
            facsecs = readmatrix(fname,'Range',[2,1], 'FileType', 'text');

            % load associated facet section fluid boundary points
            ffullname = [obj.path '/' obj.ffluid_boundary '_', strvar, '.txt'];
            fluid_boundary = readmatrix(ffullname,'Range',[2,1], 'FileType', 'text');

            % assign data to structure
            data.facid = facsecs(:,1);
            data.area = facsecs(:,2);
            data.locs = fluid_boundary(facsecs(:,3),:);
            data.distance = facsecs(:,4);
        end
    end

    % ----------------------------------------------------------------- %

    methods (Static, Access = protected)

        function M = load_solid(filename, Nx, Ny, Nz)
            % load the cells containing buildings

            % read solid_X.txt file which contains the indices where
            % buildings are located.  
            opts = detectImportOptions(filename);
            opts.DataLines = [2 Inf];
            B = readmatrix(filename, opts);

            % create a 3D array logical array to indicate the solid cells
            M = false([Nx, Ny, Nz]);
            S = sub2ind(size(M), B(:, 1), B(:, 2), B(:, 3));
            M(S) = true;
        end

        % ------------------------------------------------------------- %

        function data = loadvar(filename, svar)

            % Load netcdf data.
            % loadvar(filename, svar)
            %
            % filename : Name of file (string).
            % svar     : Variable name (string).
            %
            % data     : Array of doubles.

            data = [];

            try
                ncid = netcdf.open(filename,'NC_NOWRITE');
            catch
                disp(['Not found in directory: ', filename])
                data = [];
                return
            end

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

        % ------------------------------------------------------------- %

        function data = load_ncdata(ffullname, svar)
            % Load data from a netCDF file
            %
            % load_ncdata(ffullname, svar)
            %    ffullname: filename, potentially including full path
            %    svar:      the variable to load

            if (isempty(svar))
                % show file contents
                udbase.ncinfo(ffullname)
                data = [];
                return
            else
                svar = svar{1};
            end

            data = udbase.loadvar(ffullname, svar);
        end

        % -------------------------------------------------------------- %

        function ncinfo(ffullname)
            % Displays the data from a netCDF file in a table
            %
            % ncinfo(ffullname)
            %    ffullname: filename, potentially including full path

            [~, fileName, fileExt] = fileparts(ffullname);
            fname = [fileName, fileExt];

            % Get information about the NetCDF file
            try
                info = ncinfo(ffullname);
            catch
                disp(['Not found in directory: ', ffullname])
                %data = [];
                return
            end       

            % Initialize cell arrays to store information
            varNames = cell(length(info.Variables), 1);
            varSizes = cell(length(info.Variables), 1);
            varDims = cell(length(info.Variables), 1);
            longNames = cell(length(info.Variables), 1);
            units = cell(length(info.Variables), 1);

            % Loop through each variable to extract information
            for i = 1:length(info.Variables)
                varNames{i} = info.Variables(i).Name;  % Variable name
                varSizes{i} = info.Variables(i).Size;  % Variable size

                % Extract the dimensions
                dims = {info.Variables(i).Dimensions.Name};  % Variable dimensions names
                varDims{i} = strjoin(dims, ', ');  % Combine dimensions into a string

                % Initialize long_name and units to empty in case they are not present
                longNameAttr = '';
                unitsAttr = '';

                % Loop through the attributes to find long_name and units
                for j = 1:length(info.Variables(i).Attributes)
                    if strcmp(info.Variables(i).Attributes(j).Name, 'longname')
                        longNameAttr = info.Variables(i).Attributes(j).Value;
                    elseif strcmp(info.Variables(i).Attributes(j).Name, 'units')
                        unitsAttr = info.Variables(i).Attributes(j).Value;
                    end
                end

                longNames{i} = longNameAttr;  % Store the long name
                units{i} = unitsAttr;  % Store the units
            end

            % Convert sizes to strings for display
            varSizesStr = cellfun(@(x) sprintf('%dx', x), varSizes, 'UniformOutput', false);
            varSizesStr = cellfun(@(x) x(1:end-1), varSizesStr, 'UniformOutput', false);  % Remove trailing 'x'

            % Create the table
            T = table(varNames, longNames, units, varSizesStr, varDims, ...
                'VariableNames', {'Name', 'Description', 'Units', 'Size', 'Dimensions'});

            % convert cell arrays to character arrays
            T.('Name') = char(T.('Name'));
            T.('Description') = char(T.('Description'));
            T.('Units') = char(T.('Units'));
            T.('Size') = char(T.('Size'));
            T.('Dimensions') = char(T.('Dimensions'));

            % sort table
            T = sortrows(T, 'Name');
            
            % Display the table
            disp(['Contents of ', fname, ':'])
            disp(T);
        end
    end

    % ----------------------------------------------------------------- %

    methods (Static)
        function av = time_average(var, time, varargin)
            % A method for averaging facet variables in time. The time
            % index is assumed to be the last index of the array.

            if isempty(varargin)            
                indstart = 1;
                indstop = length(time);
            else
                tstart = varargin{1};
                tstop = varargin{2};
                indstart = find(time>=tstart, 1, 'first');
                indstop = find(time>=tstop, 1, 'first');
            end 

            av = nan;
            if length(size(var))==2
                av = mean(var(:,indstart:indstop),2);
            elseif length(size(var)) == 3
                av = mean(var(:,:,indstart:indstop),3);
            elseif length(size(var)) == 4
                av = mean(var(:,:,:, indstart:indstop),4);                
            else
                disp('Dimensions of input variable not compatible with function.')
                return
            end     
        end
    end
end