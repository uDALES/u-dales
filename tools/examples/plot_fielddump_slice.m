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

function plot_fielddump_slice(expnr, field_var, slice_var, slice_id, time_id)

% Plots a slice of fielddump in 2D and 3D.
%% Usage: matlab -nosplash -nodesktop -r "cd('tools/examples'); plot_fielddump_slice('$expnr','$field_var','$slice_var',$slice_id,$time_id); quit"

this_dir = pwd;
exp_dir = [this_dir, '/../../outputs/', expnr];
filepath = [exp_dir, '/fielddump.', expnr, '.nc'];
field = ncread(filepath, field_var);
zt = ncread(filepath, 'zt');
xt = ncread(filepath, 'xt');
yt = ncread(filepath, 'yt');
zm = ncread(filepath, 'zm');
xm = ncread(filepath, 'xm');
ym = ncread(filepath, 'ym');
time = ncread(filepath, 'time');

% Assume equidistant grid.
dx = xm(2); dy = ym(2); dz = zm(2);
xh = [xm; xm(end) + dx]; yh = [ym; ym(end) + dy]; zh = [zm; zm(end) + dz];
xf = xt; yf = yt; zf = zt;

%% Load blocks
blocks = dlmread([exp_dir, '/blocks.inp.', expnr],'',2,0);

%% Plot slice
f_2D = figure('visible', 'off');

if (field_var == 'u')
    x = xh(1:end-1) - 0.5*dx; y = yf - 0.5*dy; z = zf - 0.5*dz;
elseif (field_var == 'v')
    x = xf - 0.5*dx; y = yh(1:end-1) - 0.5*dy; z = zf - 0.5*dz;
elseif (field_var == 'w')
    x = xf - 0.5*dx; y = yf - 0.5*dy; z = zh(1:end-1) - 0.5*dz;
end

if (slice_var == 'x')
    slice_val = x(slice_id)+0.5*dx;
    pcolor(y,z,squeeze(field(slice_id,:,:,time_id))')
    axis equal
    xlim([yh(1), y(end)])
    ylim([zh(1), z(end)])
    xlabel('$y$ [m]', 'interpreter', 'latex')
    ylabel('$z$ [m]', 'interpreter', 'latex')
elseif (slice_var == 'y')
    slice_val = y(slice_id)+0.5*dy;
    pcolor(x,z,squeeze(field(:,slice_id,:,time_id))')
    axis equal
    xlim([xh(1), x(end)])
    ylim([zh(1), z(end)])
    xlabel('$x$ [m]', 'interpreter', 'latex')
    ylabel('$z$ [m]', 'interpreter', 'latex')
elseif (slice_var == 'z')
    slice_val = z(slice_id)+0.5*dz;
    pcolor(x,y,squeeze(field(:,:,slice_id,time_id))')
    axis equal
    xlim([xh(1), x(end)])
    ylim([yh(1), y(end)])
    xlabel('$x$ [m]', 'interpreter', 'latex')
    ylabel('$y$ [m]', 'interpreter', 'latex')
end

shading interp
title(['$', field_var, '(', slice_var, '=', num2str(slice_val), '\mathrm{\ m}, t=', num2str(time(time_id)), '\mathrm{\ s})$'], 'interpreter', 'latex')
cmap = redblue();
colormap(cmap)
c = colorbar;
c.Title.String = '$\mathrm{\ m\ s^{-1}}$';
c.Title.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
c.FontSize = 12;
c.TickLabelInterpreter = 'latex';
field_max = max(abs(squeeze(field(:,:,:,time_id))), [], 'all');
caxis([-field_max, field_max]);
set(gca, 'FontSize', 12)
set(gca,'ticklabelinterpreter','latex')

% Plot blocks slice
% Assumes all buildings have kl = 1 (no stacks of blocks).
% Plot buildings starting from z = 0 rather than z = 1.

clr = [0.75, 0.75, 0.75];
if (slice_var == 'x')
    for i = 1:size(blocks,1)
        il = blocks(i,1);
        iu = blocks(i,2);
        jl = blocks(i,3);
        ju = blocks(i,4);
        %kl = blocks(i,5) + 1;
        ku = blocks(i,6) + 1; % k-index starts at 0 in uDALES
        
        if ((xh(il) <= slice_val) && (slice_val <= xh(iu+1)))
            X = [yh(ju+1) yh(ju+1) yh(jl)   yh(jl)];
            Y = [0        zh(ku+1) zh(ku+1) 0     ];
            patch('XData', X, 'YData', Y, 'EdgeColor','none','FaceColor',clr)
        end
    end
elseif (slice_var == 'y')
    for i = 1:size(blocks,1)
        il = blocks(i,1);
        iu = blocks(i,2);
        jl = blocks(i,3);
        ju = blocks(i,4);
        %kl = blocks(i,5) + 1;
        ku = blocks(i,6) + 1; % k-index starts at 0 in uDALES
        
        if ((yh(jl) <= slice_val) && (slice_val <= yh(ju+1)))
            X = [xh(il)   xh(il)   xh(iu+1) xh(iu+1)];
            Y = [0        zh(ku+1) zh(ku+1) 0       ];
            patch('XData', X, 'YData', Y, 'EdgeColor','none','FaceColor',clr)
        end
    end
elseif (slice_var == 'z')
    for i = 1:length(blocks)
        il = blocks(i,1);
        iu = blocks(i,2);
        jl = blocks(i,3);
        ju = blocks(i,4);
        %kl = blocks(i,5) + 1;
        ku = blocks(i,6) + 1; % k-index starts at 0 in uDALES
               
        if (slice_val <= zh(ku+1))
            X = [xh(il)   xh(iu+1) xh(iu+1) xh(il)  ];
            Y = [yh(jl)   yh(jl)   yh(ju+1) yh(ju+1)];
            patch('XData', X, 'YData', Y, 'EdgeColor','none','FaceColor',clr)
        end
    end
end

filename = ['fielddump_slice_2D_', expnr, '_', field_var, '(', slice_var, '=', num2str(slice_val), 'm,t=', num2str(time(time_id)), 's)'];
saveas(f_2D,[exp_dir, '/', filename, '.png'])

%% Plot slice in 3D space
f_3D = figure('visible', 'off');
light('Position', [-0.75 -0.35 1], 'Style', 'infinite');

% Plot blocks
% Assumes all buildings have kl = 1 (no stacks of blocks).
% Plot buildings starting from z = 0 rather than z = 1.
% Plot floors with side faces.

clr = [0.85, 0.85, 0.85];
for i = 1:length(blocks)
    il = blocks(i,1);
    iu = blocks(i,2);
    jl = blocks(i,3);
    ju = blocks(i,4);
    %kl = blocks(i,5) + 1;
    ku = blocks(i,6) + 1;  % k-index starts at 0 in uDALES
    
    V = [xh(il) yh(jl)        0; xh(iu+1) yh(jl)        0; xh(iu+1) yh(ju+1)        0; xh(il) yh(ju+1)        0; ... % Lower 4 vertices
         xh(il) yh(jl) zh(ku+1); xh(iu+1) yh(jl) zh(ku+1); xh(iu+1) yh(ju+1) zh(ku+1); xh(il) yh(ju+1) zh(ku+1)];    % Upper 4 vertices

    %    Bottom,  East,    North,   West,    South,   Top
    F = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];     
    patch('Faces',F,'Vertices',V,'EdgeColor','none','FaceColor',clr,'FaceLighting','flat','AmbientStrength', 0.3, 'SpecularStrength', 0.6, 'DiffuseStrength', 0.7);
end

axis equal
hold on
view(3)
xlim([xh(1), xh(end)])
ylim([yh(1), yh(end)])
zlim([zh(1), zh(end)])
xlabel('$x$ [m]', 'interpreter', 'latex')
hXLabel = get(gca,'XLabel');
set(hXLabel,'rotation', 25)
ylabel('$y$ [m]', 'interpreter', 'latex')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation', -35)
zlabel('$z$ [m]', 'interpreter', 'latex')

% Plot slice

if (field_var == 'u')
    x = xh(1:end-1) - 0.5*dx; y = yf - 0.5*dy; z = zf - 0.5*dz;
elseif (field_var == 'v')
    x = xf - 0.5*dx; y = yh(1:end-1) - 0.5*dy; z = zf - 0.5*dz;
elseif (field_var == 'w')
    x = xf - 0.5*dx; y = yf - 0.5*dy; z = zh(1:end-1) - 0.5*dz;
end

if (slice_var == 'x')
    slice_val = x(slice_id)+0.5*dx;
    [Y,Z] = meshgrid(y,z);
    X = slice_val + zeros(size(Y));
    s = surf(X,Y,Z,squeeze(field(slice_id,:,:,time_id))');
elseif (slice_var == 'y')
    slice_val = y(slice_id)+0.5*dy;
    [X,Z] = meshgrid(x,z);
    Y = slice_val + zeros(size(Z));
    s = surf(X,Y,Z,squeeze(field(:,slice_id,:,time_id))');
elseif (slice_var == 'z')
    slice_val = z(slice_id)+0.5*dz;
    [X,Y] = meshgrid(x,y);
    Z = slice_val + zeros(size(X));
    s = surf(X,Y,Z,squeeze(field(:,:,slice_id,time_id))');
end

set(s, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceLighting', 'flat', 'AmbientStrength', 1)
title(['$', field_var, '(', slice_var, '=', num2str(slice_val), '\mathrm{\ m}, t=', num2str(time(time_id)), '\mathrm{\ s})$'], 'interpreter', 'latex')
cmap = redblue();
colormap(cmap)
c = colorbar;
c.Title.String = '$\mathrm{\ m\ s^{-1}}$';
c.Title.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
c.FontSize = 12;
field_max = max(abs(squeeze(field(:,:,:,time_id))), [], 'all');
caxis([-field_max, field_max]);
set(gca, 'FontSize', 12)
set(gca,'ticklabelinterpreter','latex')

filename = ['fielddump_slice_3D_', expnr, '_', field_var, '(', slice_var, '=', num2str(slice_val), 'm,t=', num2str(time(time_id)), 's)'];
saveas(f_3D,[exp_dir, '/', filename, '.png'])

function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.
%   Adam Auton, 9th October 2009
if nargin < 1, m = size(get(gcf,'colormap'),1); end
if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end
c = [r g b];
end

end
