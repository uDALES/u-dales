function [Xmean, var] = merge_stat_var(X, XpXp, n)
% merge_stat_var  Merge short-time variance statistics into longer-time averages
%
%   [Xmean, var] = merge_stat_var(X, XpXp, n)
%
% Helper wrapper that computes time-averaged mean and variance for a
% single variable. This is a thin wrapper around merge_stat_cov where
% X and Y are the same variable. The total variance inside each window
% is computed as mean( (X - Xmean).^2 + XpXp ) where XpXp is the
% instantaneous contribution (e.g. subgrid or measurement variance).
%
% Inputs
%   X     - Variable time series with time in the last dimension.
%   XpXp  - Instantaneous contribution to variance (same shape as X).
%   n     - Window length in number of time samples.
%
% Outputs
%   Xmean - Time-averaged X in each non-overlapping window.
%   var   - Time-averaged total variance in each window.
%
% Example:
%   % Average X every 20 samples with no instantaneous term:
%   [Xm, V] = merge_stat_var(X, zeros(size(X)), 20);
%
% See also: merge_stat_cov

% Delegate to merge_stat_cov with X==Y
[Xmean, ~, var] = merge_stat_cov(X, X, XpXp, n);
end