function [Xmean, Ymean, cov] = merge_stat_cov(X, Y, XpYp, n)
% merge_stat_cov  Merge short-time statistics into longer-time averages
%
%   [Xmean, Ymean, cov] = merge_stat_cov(X, Y, XpYp, n)
%
% This utility combines short-time (instantaneous or high-frequency)
% statistics into longer-time averaged statistics. The routine groups the
% input time series into non-overlapping windows of length n and computes
% the mean of each variable as well as the total covariance inside each
% window. It is intended for workflows that produce many short-time
% statistics which then need to be merged into coarser, long-time
% statistics without re-reading all original data.
%
% Inputs
%   X       - First variable time series. The last dimension is assumed
%             to be time (e.g. size(X) = [Nx, Ny, Nt] or [Nt]).
%   Y       - Second variable time series (same shape as X).
%   XpYp    - Pre-computed instantaneous covariance series of the same
%             size as X and Y (i.e. instantaneous contribution that
%             should be added to the product of fluctuations). For
%             variance use XpYp = XpXp where XpXp is instantaneous
%             variance contribution.
%   n       - Window length (number of consecutive time samples to
%             average together). If n >= Nt the entire record is
%             averaged and a single output time sample is returned.
%
% Outputs
%   Xmean   - Time-averaged X inside each non-overlapping window
%             (same dimensions as X but with the time dimension replaced
%             by the number of windows).
%   Ymean   - Time-averaged Y inside each window (same shape as Xmean).
%   cov     - Time-averaged total covariance inside each window. The
%             total covariance inside a window is computed as
%                 mean( (X - Xmean_window).*(Y - Ymean_window) + XpYp )
%             where the mean is taken over the time samples in the window.
%
% Example:
%   % Merge every 10 samples:
%   [Xm, Ym, C] = merge_stat_cov(X, Y, zeros(size(X)), 10);
%
% See also: merge_stat_var

% Get the number of time points
dims = size(X);
t = dims(end);   % time is always the last dimension

% If n is larger than available time points, use all time points
if n >= t
    n = t;
end

% Calculate number of averaged sets and points to discard
N = fix(t/n);       % number of sets after averaging
shift = rem(t,n);   % remaining time points to discard from start

% Create arrays for storing results with appropriate dimensions
output_dims = dims;
output_dims(end) = N;  % last dimension becomes number of averaging sets

% Initialize output arrays to store results
Xmean = zeros(output_dims);
Ymean = zeros(output_dims);
cov = zeros(output_dims);

% Process each time averaging window
for i = 1:N
    % Calculate start and end indices for current window
    st = 1 + shift + (i-1)*n;
    ed = n + shift + (i-1)*n;
    
    % Create indexing cell array for all dimensions
    indices = arrayfun(@(x) ':', 1:ndims(X)-1, 'UniformOutput', false);
    time_window = {st:ed};  % Time window for current averaging period
    
    % Extract data for current time window using dynamic indexing
    X_window = X(indices{:}, time_window{:});
    Y_window = Y(indices{:}, time_window{:});
    XpYp_window = XpYp(indices{:}, time_window{:});
    
    % Calculate means for current window along time dimension
    X_mean = mean(X_window, ndims(X));
    Y_mean = mean(Y_window, ndims(Y));
    
    % Store means in output arrays
    Xmean(indices{:}, i) = X_mean;
    Ymean(indices{:}, i) = Y_mean;
    
    % Initialize array for covariance calculation
    temp_cov = zeros(size(X_window));
    
    % Calculate covariance for each time point in window
    for j = 1:n
        % Get instantaneous values
        X_instant = X_window(indices{:}, j);
        Y_instant = Y_window(indices{:}, j);
        XpYp_instant = XpYp_window(indices{:}, j);
        
        % Calculate fluctuations from window mean
        X_fluct = X_instant - X_mean;
        Y_fluct = Y_instant - Y_mean;
        
        % Calculate total covariance: (X-Xmean)*(Y-Ymean) + XpYp
        temp_cov(indices{:}, j) = X_fluct .* Y_fluct + XpYp_instant;
    end
    
    % Average covariance over time window
    cov(indices{:}, i) = mean(temp_cov, ndims(X));
end
end