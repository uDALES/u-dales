classdef da_pp < dynamicprops
    %description goes here
    %
    properties (Hidden = true, SetAccess = protected)
        path;                    % Path to simulations.
        cpath;                   % Current path.
        buf;                     % All purpose buffer.
    end
    
    properties (SetAccess = protected)
        g = 9.81;                % Gravitational acceleration.
        
    end
    methods (Static)
        
        function obj = da_pp(expnr, varargin);
            % Class constructor.
            % Reads simparams_expnr.m and fills object parameters.
            % spbase(expnr, varargin)
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
            
            
            expstr = num2str( expnr, '%3.3d');
            fid = fopen(['namoptions.', expstr]);
            
            if (fid == -1)
                disp(['namoptions.', expstr, ' not found. Exiting ...'])
                return
            end
            
            line = fgets(fid);
            
            % load each line of simparams_* into object
            while min(line) > 0;
                toks = regexp(line, TOKENS, 'tokens');
                if ~isempty(toks)
                    lhs = regexprep(toks{1}{1}, WHITE, VOID);
                    h = addprop(obj, lhs);
                    rhs=regexprep(toks{1}{2}, WHITE, VOID);
                    if strcmp(rhs,'.false.'); % string
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