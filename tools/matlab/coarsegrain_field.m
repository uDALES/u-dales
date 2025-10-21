function var_filtered = coarsegrain_field(var, Lflt, xm, ym)
% coarsegrain_field  Apply 2D spatial filter to 3D field data.
%
%   var_filtered = coarsegrain_field(var, Lflt, xm, ym)
%
% This function applies 2D spatial (x-y) coarse-graining filters to 3D field data.
% Multiple filter sizes are applied simultaneously, creating a 4D output
% where the 4th dimension corresponds to different filter sizes.
%
% Inputs
%   var     - 3D field data with dimensions [itot, jtot, ktot] where
%             the first two dimensions are horizontal (x, y) and the
%             third is vertical (z) or time. Field is assumed periodic.
%   Lflt    - Array of filter lengths in physical units (meters)
%   xm      - x-coordinates of grid points (meters)
%   ym      - y-coordinates of grid points (meters)
%            
% Outputs
%   var_filtered - 4D filtered data with dimensions [itot, jtot, ktot, length(Lflt)]
%                  where the 4th dimension corresponds to different filter sizes
%
% Algorithm
%   - Converts physical filter lengths to grid cell numbers (Ng = round(Lflt/dx))
%   - Works at discrete level with normalized periodic filters  
%   - Uses FFT-based convolution for computational efficiency
%
% Example:
%   % Apply multiple filter sizes to velocity field
%   filter_lengths = [10, 20, 40, 80, 160]; % Physical lengths in meters
%   u_filtered = coarsegrain_field(u_data, filter_lengths, xm, ym);

% Validate inputs
if ~isnumeric(var) || ndims(var) ~= 3
    error('Input var must be a 3D numeric array');
end

if ~isnumeric(Lflt) || ~isvector(Lflt) || any(Lflt <= 0)
    error('Lflt must be a vector of positive numbers (physical lengths in meters)');
end

% Get grid information
if isempty(xm) || isempty(ym)
    error('Grid information (xm, ym) not provided.');
end

% Convert physical filter lengths to grid cell numbers
dx = xm(2) - xm(1); % x-direction is assumed to be equidistant
dy = ym(2) - ym(1); % y-direction is assumed to be equidistant

% Compute half-widths in i- and j-directions
Ngx = max(round(Lflt/2 / dx), 1); % half-width in i (x) grid cells
Ngy = max(round(Lflt/2 / dy), 1); % half-width in j (y) grid cells

% Get grid dimensions
[itot, jtot, ktot] = size(var);

% Initialize output array
var_filtered = zeros(itot, jtot, ktot, length(Lflt));

tic;
for n = 1:length(Lflt) % Loop over all filter sizes
    ngx = Ngx(n);
    ngy = Ngy(n);
    actual_Lflt_x = ngx * dx * 2; % Actual physical length in x (full width)
    actual_Lflt_y = ngy * dy * 2; % Actual physical length in y (full width)
    
    % Vectorized creation of 2D filter kernel for periodic domain
    [I, J] = ndgrid(1:itot, 1:jtot);
    di = min(I - 1, itot - (I - 1)); % periodic distance in i-direction
    dj = min(J - 1, jtot - (J - 1)); % periodic distance in j-direction
    f2d = (di <= ngx) & (dj <= ngy); % logical mask for filter    
    
    % Normalize the filter to preserve mean values
    f2d = double(f2d / sum(f2d(:)));
    
    % Apply FFT-based filtering using periodic convolution
    var_filtered(:, :, :, n) = apply_fft_filter_periodic(f2d, var);
    fprintf(' Filter %d/%d (Lflt_x=%.1fm, Lflt_y=%.1fm) completed\n', ...
            n, length(Lflt(n)), actual_Lflt_x, actual_Lflt_y );
end
elapsed_time = toc;
fprintf('Coarse-graining completed in %.2f seconds\n', elapsed_time);
end

function gf = apply_fft_filter_periodic(f2d, g)
% apply_fft_filter_periodic  Apply 2D filter to 3D data using periodic FFT convolution
%
%   gf = apply_fft_filter_periodic(f2d, g)
%
% Helper function that applies a 2D filter to each horizontal slice
% of a 3D field using FFT-based convolution with proper periodic handling.
%
% Inputs
%   f2d  - 2D filter kernel [itot, jtot]
%   g    - 3D field data [itot, jtot, ktot]
%
% Outputs
%   gf   - Filtered 3D field data [itot, jtot, ktot]

% The FFT convolution automatically handles periodic boundary conditions
% when the filter and data have the same size

% Compute 2D FFT of filter once
f2dhat = fft2(f2d);

% Initialize output
gf = zeros(size(g));

% Apply filter to each vertical level
for k = 1:size(g, 3)
    % FFT-based convolution with periodic boundaries
    gf(:, :, k) = real(ifft2(f2dhat .* fft2(g(:, :, k))));
end
end
