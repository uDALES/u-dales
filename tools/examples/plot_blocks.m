function plot_blocks(expnr)

%% Usage: matlab -nosplash -nodesktop -r "cd('tools/examples'); plot_blocks('$expnr'); quit"

this_dir = pwd;
exp_dir = [this_dir, '/../../examples/', expnr];

blocks = dlmread([exp_dir, '/blocks.inp.', expnr],'',2,0);
xf = dlmread([exp_dir, '/xgrid.inp.', expnr],'',2,0);
zf = dlmread([exp_dir, '/zgrid.inp.', expnr],'',2,0);
xh = zeros(length(xf)+1,1);
zh = zeros(length(zf)+1,1);

% Define xh and yh in the same way as initglobal.
xh(1) = 0;
for i = 1:length(xf)
    xh(i+1) = xh(i) + 2.0*(xf(i) - xh(i));
end

zh(1) = 0;
for k = 1:length(zf)
    zh(k+1) = zh(k) + 2.0*(zf(k) - zh(k));
end

ysize = 64;
jtot = 64;

% read namoptions to get ysize and jtot.
fid = fopen([exp_dir, '/namoptions.', expnr]);

if (fid == -1)
    disp(['namoptions.', expnr, ' not found. Exiting ...'])
    return
end

TOKENS = '(.*)\=(.*)';
WHITE = '\s*';
VOID = '';
line = fgets(fid);
while min(line) > 0
    toks = regexp(line, TOKENS, 'tokens');
    if ~isempty(toks)
        lhs = regexprep(toks{1}{1}, WHITE, VOID);
        rhs = regexprep(toks{1}{2}, WHITE, VOID);
        if strcmp(lhs,'ysize')
            ysize = str2double(rhs);
        elseif strcmp(lhs,'jtot')
            jtot = str2double(rhs);
        end
    end
    line = fgets(fid);
end

dy = ysize / jtot;
yh = 0:dy:ysize;

% Plot blocks
% Assumes all buildings have kl = 1 (no stacks of blocks).
% Plot buildings starting from z = 0 rather than z = 1.
% Plot floors with side faces.

f = figure('visible', 'off');
light('Position', [-0.75 -0.35 1], 'Style', 'infinite');

clr = [0.85, 0.85, 0.85];
for i = 1:size(blocks,1)
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
grid on
view(3)
xlim([xh(1), xh(end)])
ylim([yh(1), yh(end)])
zlim([zh(1), zh(end)])
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

filename = ['blocks.', expnr];
saveas(f,[exp_dir, '/', filename, '.png'])

end