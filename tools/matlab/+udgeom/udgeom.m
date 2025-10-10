% udgeom Geometry class for uDALES
%    The udgeom class contains the triangulated surface.

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

classdef udgeom < handle
   % Geometry class for uDALES

   properties (Hidden = true, SetAccess = protected)
      path;                    % Path to geometry.
      cpath;                   % Current path.
   end

   properties (SetAccess = protected)
      stl;                     % stl of the geometry
   end   
   
   methods
      % geometry constructor
      
      function obj = udgeom(varargin)
         % Class constructor.
         %
         % varargin : Path for geometry
         %            Or a triangulation object
         
         % store current and simulation directory
         obj.cpath = pwd;
         
         dapath = obj.cpath;
         if (nargin == 1)
             % Check if the input argument is a string
                if ischar(varargin{1}) || isstring(varargin{1})
                    dapath = varargin{1};
                % Check if the input argument is a triangulation object
                elseif isa(varargin{1}, 'triangulation')
                    obj.stl = varargin{1};
                    % You can add more logic here if needed for triangulation object input
                else
                    error('Input must be either a string or a triangulation object.');
                end
         end
         cd([dapath]);
         obj.path = pwd;
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
      
      % -------------------------------------------------------------- %
      
      function load(obj, filename)
         % load an STL file.
         % 
         % example:
         %   obj.load(obj, filename)
          
         obj.gopath()
         obj.stl = stlread(filename);
         obj.gohome()
      end
      
      % -------------------------------------------------------------- %
      
      function save(obj, filename)
         % save an STL file.
         %
         % example:
         %   obj.save(obj, filename)

         obj.gopath()
         stlwrite(obj.stl, filename)
         obj.gohome()
      end
      
      % -------------------------------------------------------------- %
      
      function show(obj, varargin)
         % plot the geometry
         %
         % show(obj, color_buildings, plot_quiver)
         %
         %   color_buildings (optional): boolean parameter on whether
         %                               to color buildings. Default =
         %                               true. Needs to be set to false for 
         %                               large geometries.
         %
         %   plot_quiver (optional): boolean parameter on whether
         %                           to plot quiver arrows. Default = true.
         %
         % Examples:
         %   obj.show();                   % both true
         %   obj.show(false);              % color_buildings = false, plot_quiver = true
         %   obj.show(true, false);        % color_buildings = true, plot_quiver = false
    
         % --- Default values ---
         color_buildings = true;
         plot_quiver = true;

         % --- Handle variable arguments ---
         if ~isempty(varargin)
             color_buildings = varargin{1};
         end
         if numel(varargin) >= 2
             plot_quiver = varargin{2};
         end

         faceNormals = faceNormal(obj.stl);
         incenters = incenter(obj.stl);
         
         figure
         patch('Faces', obj.stl.ConnectivityList, 'Vertices', obj.stl.Points, ...
               'FaceColor', ones(3,1)*0.85, 'FaceAlpha', 1) 
         hold on
         
         if (color_buildings)
             for i=1:length(incenters(:,1))
                 if incenters(i,3)>0
                     A = obj.stl.Points(obj.stl.ConnectivityList(i,1),:);
                     B = obj.stl.Points(obj.stl.ConnectivityList(i,2),:);
                     C = obj.stl.Points(obj.stl.ConnectivityList(i,3),:);
                     fill3([A(1) B(1) C(1)], [A(2) B(2) C(2)], [A(3) B(3) C(3)], ...
                         [0.73 0.83 0.96])
                 end
             end
         end
         
         if (plot_quiver)
             quiver3(incenters(:,1), incenters(:,2), incenters(:,3), ...
                 faceNormals(:,1), faceNormals(:,2), faceNormals(:,3), 0.2)
         end

         view(3)
         axis equal tight
         set(gca,"Box","on")
         set(gca,"BoxStyle","full")
      end
   end
end
