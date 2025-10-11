function varargout = time_average(varargin)
% time_average  Time-average variables over all available time intervals
%
%   [Xmean, var] = time_average(X)
%       Compute the time-averaged mean and variance of X over the full
%       time record (time is assumed to be the last dimension).
%
%   [Xmean, Ymean, cov] = time_average(X, Y)
%       Compute the time-averaged means of X and Y and their covariance
%       over the full time record.
%
% Inputs:
%   X   - First variable time series. The final array dimension is assumed
%         to be time (e.g. size(X) = [..., Nt]).
%   Y   - Second variable time series (same shape as X) [optional]
%
% Outputs:
%   For single variable (1 input):
%     Xmean - Time-averaged X over the full record
%     var   - Time-averaged variance over the full record
%
%   For two variables (2 inputs):
%     Xmean - Time-averaged X over the full record
%     Ymean - Time-averaged Y over the full record
%     cov   - Time-averaged covariance over the full record
%
% Examples:
%   [X_avg, X_var] = time_average(X);
%   [X_avg, Y_avg, XY_cov] = time_average(X, Y);

if nargin == 1
    % Single variable case: time_average(X)
    X = varargin{1};
    
    % Get the number of time points and use all of them
    dims = size(X);
    n = dims(end);   % time is the last dimension
    
    % Use merge_stat_var with zero instantaneous variance contribution
    XpXp = zeros(size(X));
    [Xmean, var] = merge_stat_var(X, XpXp, n);
    
    varargout{1} = Xmean;
    varargout{2} = var;
    
elseif nargin == 2
    % Two variable case: time_average(X, Y)
    X = varargin{1};
    Y = varargin{2};
    
    % Get the number of time points and use all of them
    dims = size(X);
    n = dims(end);   % time is the last dimension
    
    % Use merge_stat_cov with zero instantaneous covariance contribution
    XpYp = zeros(size(X));
    [Xmean, Ymean, cov] = merge_stat_cov(X, Y, XpYp, n);
    
    varargout{1} = Xmean;
    varargout{2} = Ymean;
    varargout{3} = cov;
    
else
    error('time_average requires 1 or 2 input arguments. Usage: time_average(X) or time_average(X, Y)');
end

end