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

% Copyright (C) 2016-2024 the uDALES Team.

classdef udbase < dynamicprops
    %  Base class for uDALES simulations

    properties (SetAccess = protected)
        expnr;                   % the simulation ID
        geom;                    % geometry instance of udgeom

        % 3D logical arrays that identify the cells that are solid. Note
        % that there are four arrays since the underlying grid is staggered.
        Su;                      % 3d array for variables at u-locations
        Sv;                      % 3d array for variables at v-locations
        Sw;                      % 3d array for variables at w-locations
        Sc;                      % 3d array for variables at c-locations
    end
    properties (Hidden = true, SetAccess = protected)
        path;                    % Path to simulations.
        cpath;                   % Current path.
        buf;                     % All purpose buffer.

        % standard filenames
        fxytdump   = 'xytdump';
        ftdump     = 'tdump';
        ffielddump = 'fielddump';
        fislicedump = 'islicedump';
        fjslicedump = 'jslicedump';
        fkslicedump = 'kslicedump';
        fsolid = 'solid';
    end

    methods
        function obj = udbase(expnr, varargin)
            % Class constructor.
            %
            % expnr    : Integer equal to simulation number.
            % varargin : Path to simulations.

            obj.expnr = expnr;

            % store current and simulation directory
            obj.cpath = pwd;

            dapath = obj.cpath;
            if (nargin>1)
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
            
            % ------ %
            % read the stl file of the geometry if present

            if isprop(obj, 'stl_file')  && exist(obj.stl_file, 'file')
                geom = udgeom.udgeom(obj.path);
                geom.load(obj.stl_file);
                obj.geom = geom;
            else
                disp('STL file not found, implying a flat surface. obj.geom is empty')
                 geom = [];
            end
            obj.geom = geom;

            % ------ %
            % load 3D arrays that indicate the location of solid regions.
            % [note that there are four arrays since the domain is
            % staggered]. 

            obj.Su = obj.load_solid([obj.fsolid, '_u.txt'], obj.itot, obj.jtot, obj.ktot);
            obj.Sv = obj.load_solid([obj.fsolid, '_v.txt'], obj.itot, obj.jtot, obj.ktot);
            obj.Sw = obj.load_solid([obj.fsolid, '_w.txt'], obj.itot, obj.jtot, obj.ktot);
            obj.Sc = obj.load_solid([obj.fsolid, '_c.txt'], obj.itot, obj.jtot, obj.ktot);

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

end