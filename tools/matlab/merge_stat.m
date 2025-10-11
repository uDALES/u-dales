function varargout = merge_stat(varargin)
% merge_stat  Merge short-time statistics into longer-time averages
%
%   Xmean = merge_stat(X, n)
%       Computes time-averaged mean for a single variable X.
%       Groups the time series into non-overlapping windows of length n
%       and computes statistics inside each window.
%
%   [Xmean, var] = merge_stat(X, XpXp, n)
%       Computes time-averaged mean and variance for a single variable X.
%       Groups the time series into non-overlapping windows of length n
%       and computes statistics inside each window.
%
%   [Xmean, Ymean, cov] = merge_stat(X, Y, XpYp, n)  
%       Computes time-averaged means and covariance for two variables X and Y.
%       Groups the time series into non-overlapping windows and computes
%       statistics inside each window.
%
% Inputs:
%   X    - First variable time series (time in final dimension)
%   Y    - Second variable time series (same shape as X) [two-variable case only]
%   XpXp - Instantaneous variance contribution (same shape as X) [single-variable]
%   XpYp - Instantaneous covariance contribution (same shape as X and Y) [two-variable]  
%   n    - Window length (number of time samples per averaged window)
%
% Outputs:
%   Xmean - Time-averaged X in each window
%   Ymean - Time-averaged Y in each window [two-variable case only]
%   var   - Time-averaged variance in each window [single-variable case]
%   cov   - Time-averaged covariance in each window [two-variable case]
%
% Examples:
%    X_avg = merge_stat(X, 20);
%   [X_avg, X_var] = merge_stat(X, XpXp, 20);
%   [X_avg, Y_avg, XY_cov] = merge_stat(X, Y, XpYp, 50);
if nargin == 2
    % Single-variable case: compute mean and variance
    [Xmean, ~] = merge_stat_var(varargin{1}, zeros(size(varargin{1})), varargin{2});
    varargout{1} = Xmean;

elseif nargin == 3
    % Single-variable case: compute mean and variance
    [Xmean, var] = merge_stat_var(varargin{1}, varargin{2}, varargin{3});
    varargout{1} = Xmean;
    varargout{2} = var;

elseif nargin == 4
    % Two-variable case: compute means and covariance
    [Xmean, Ymean, cov] = merge_stat_cov(varargin{1}, varargin{2}, varargin{3}, varargin{4});
    varargout{1} = Xmean;
    varargout{2} = Ymean;
    varargout{3} = cov;

else
    error('merge_stat requires either 3 arguments (X, XpXp, n) or 4 arguments (X, Y, XpYp, n)');
end
end

% -------------------------------------------------------------------------
function [Xmean, var] = merge_stat_var(X, XpXp, n)
% merge_stat_var  Merge short-time variance statistics into longer-time averages
%
%   [Xmean, var] = merge_stat_var(X, XpXp, n)
%
% Thin wrapper that computes time-averaged mean and variance for a
% single variable by delegating to merge_stat_cov with Y==X.

[Xmean, ~, var] = merge_stat_cov(X, X, XpXp, n);
end

% -------------------------------------------------------------------------
function [Xmean, Ymean, cov] = merge_stat_cov(X, Y, XpYp, n)
% merge_stat_cov  Merge short-time statistics into longer-time averages
%
%   [Xmean, Ymean, cov] = merge_stat_cov(X, Y, XpYp, n)
%
% Combine short-time statistics into longer-time averages by grouping the
% inputs into non-overlapping windows of length n and computing means and
% total covariance inside each window. Time is always the last dimension.

% Get the number of time points (last dimension)
dims = size(X);
t = dims(end);

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

    % Create indexing cell array for all non-time dimensions
    nd = ndims(X);
    indices = repmat({':'}, 1, nd-1);
    time_window = st:ed;

    % Extract data for current time window using dynamic indexing
    X_window = X(indices{:}, time_window);
    Y_window = Y(indices{:}, time_window);
    XpYp_window = XpYp(indices{:}, time_window);

    % Calculate means for current window along time dimension
    X_mean = mean(X_window, nd);
    Y_mean = mean(Y_window, nd);

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
    cov(indices{:}, i) = mean(temp_cov, nd);
end
end