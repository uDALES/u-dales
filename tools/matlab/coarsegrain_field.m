function var_filtered = coarsegrain_field(var, Lflt, dx, xm, ym)
% coarsegrain_field  Apply spatial filtering to 3D field data using FFT-based convolution
%
%   var_filtered = coarsegrain_field(var, Lflt, dx, xm, ym)
%
% This function applies spatial coarse-graining filters to 3D field data.
% Multiple filter sizes are applied simultaneously, creating a 4D output
% where the 4th dimension corresponds to different physical filter lengths.
%
% Inputs
%   var     - 3D field data with dimensions [itot, jtot, ktot] where
%             the first two dimensions are horizontal (x, y) and the
%             third is vertical (z) or time
%   Lflt    - Array of filter lengths in physical units (meters)
%            
%   dx      - Grid spacing in x-direction (meters)
%   xm      - x-coordinates of grid points (meters)
%   ym      - y-coordinates of grid points (meters)
%
% Outputs
%   var_filtered - 4D filtered data with dimensions [itot, jtot, ktot, length(Lflt)]
%                  where the 4th dimension corresponds to different filter sizes
%
% Algorithm
%   - Converts physical filter lengths to grid cell numbers (Ng = round(Lflt/dx))
%   - Creates square filters with half-width Lflt for each filter size
%   - Uses distance from domain boundaries to handle periodic conditions
%   - Applies FFT-based convolution for computational efficiency
%   - Normalizes each filter to preserve mean values
%   - Displays conversion information and actual filter sizes used
%
% Example:
%   % Apply multiple filter sizes to velocity field
%   filter_lengths = [10, 20, 40, 80, 160]; % Physical lengths in meters
%   u_filtered = coarsegrain_fields(u_data, filter_lengths, obj.dx, obj.xm, obj.ym);
%
% See also: apply_fft_filter

% Validate inputs
if ~isnumeric(var) || ndims(var) ~= 3
    error('Input var must be a 3D numeric array');
end

if ~isnumeric(Lflt) || ~isvector(Lflt) || any(Lflt <= 0)
    error('Lflt must be a vector of positive numbers (physical lengths in meters)');
end

% Get grid information
if isempty(xm) || isempty(ym)
    error('Grid information (xm, ym) not available.');
end

% Convert physical filter lengths to grid cell numbers
% Ng represents width of filter in grid cells
Ng = round(Lflt/2 / dx); % Round to nearest integer for grid cells
Ng = max(Ng, 1); % Ensure minimum filter size is 1 grid cell

% Display conversion information
% fprintf('Converting physical filter lengths to grid cells:\n');
for i = 1:length(Lflt)
    actual_Lflt = Ng(i) * dx; % Actual physical length after rounding
%    fprintf('  Lflt=%.1fm -> Ng=%d cells (actual=%.1fm, dx=%.2fm)\n', ...
%            Lflt(i), Ng(i), actual_Lflt, dx);
end

% Create coordinate grids
[Xm, Ym] = ndgrid(xm - xm(1), ym - ym(1));
Lx = xm(end);
Ly = ym(end);

% Calculate distance from boundaries (for periodic conditions)
Dx = min(abs(Xm), abs(Xm - Lx));
Dy = min(abs(Ym), abs(Ym - Ly));

% Initialize output array
var_filtered = zeros(size(var, 1), size(var, 2), size(var, 3), length(Lflt));

% Apply filtering for each filter size
fprintf('Applying coarse-graining filters...\n');
tic;

for n = 1:length(Lflt)
    % Use the calculated Ng value and actual physical length
    actual_Lflt = Ng(n) * dx;
    
    % Create 2D filter
    f2d = zeros(length(xm), length(ym));
    f2d(max(Dx, Dy) < actual_Lflt) = 1;
    f2d = f2d / sum(f2d(:)); % normalize the filter
    
    % Apply FFT-based filtering
    var_filtered(:, :, :, n) = apply_fft_filter(f2d, var);
    
    fprintf('  Filter %d/%d (Lflt=%.1fm) completed\n', ...
            n, length(Lflt), actual_Lflt*2);
end

elapsed_time = toc;
fprintf('Coarse-graining completed in %.2f seconds\n', elapsed_time);
end

function gf = apply_fft_filter(f2d, g)
% apply_fft_filter  Apply 2D filter to 3D data using FFT convolution
%
%   gf = apply_fft_filter(f2d, g)
%
% Helper function that applies a 2D filter to each horizontal slice
% of a 3D field using FFT-based convolution.
%
% Inputs
%   f2d  - 2D filter kernel [itot, jtot]
%   g    - 3D field data [itot, jtot, ktot]
%
% Outputs
%   gf   - Filtered 3D field data [itot, jtot, ktot]

% Compute 2D FFT of filter once
f2dhat = fft2(f2d);

% Initialize output
gf = zeros(size(g));

% Apply filter to each vertical level
for k = 1:size(g, 3)
    gf(:, :, k) = real(ifft2(f2dhat .* fft2(g(:, :, k))));
end
end