function [stats, F]=dfa(data, boxes, opt)
% This function computes the scaling exponent by means of detrended fluctuation analysis.
% Input arguments: 
% 'data' - (required) data vector to compute scaling exponent of; 
% 'boxes' - (required) vector of scales for detrending; must be integers;  
% 'opt' - (optional) degree of polynomial for detrending; default - 1
% 
% Output variables: 
% 'stats' - [scaling_exponent intercept norm_of_residuals]
% 'F' - vestor of residual variance after detrending on scales specifed in
% 'boxes'; useful for graphical output

% version history
% 20200722
% code cleaning 
% 20200429
% updated help contents
% added copyright line
% 20200318 v3.0
% force input data to single vector
% removed 'total' from output
% reduced the number of intermediate variables for speed
% 20200306 v2.0
% changed use of 'opt' input argument;
% updated intial check for input arguments; 
% updated code for usage of 'detrend' to do detrending on column vectors; 
% optimized calculation of residual variance


narginchk(2, 3)
% make sense of input arguments
if nargin==2 % polynomial fit degree not specified
    opt=1; % default: linear
else % nargin==3, but varargin{3} is not 1 or 0
    if opt<1 
        opt=1; 
        warning 'Polynomial degree for detrending must be at least 1 (linear)';
    end
end

data=data(:); boxes=boxes(:); % force to column vector regardless initial format

% prealocate F for speed
F=zeros(length(boxes), 1);
% integrate the time series; use this as input
intser=cumsum(detrend(data, 0));

for idx=1:length(boxes) % process one scale at a time
    nmax=floor(length(intser)/boxes(idx)); % max number of boxes that fit in the input dataset
    
    % truncate data to integer number of boxes and reshape for detrending (columns of box length)
    idata=reshape(intser(1:nmax*boxes(idx)), boxes(idx), nmax); 
    
    % perform detrending by column
    idata=detrend(idata, opt);
    
    % store RMS fluctuation in F - second output
    F(idx)=sqrt(mean(idata(:).^2)); 
end % next box size

% compute scaling exponent (alpha), intercept (beta), and the norm of residuals (S.normr)
[p, S]=polyfit(log10(boxes(:)), log10(F(:)), 1);

% compile output
stats=[p S.normr]; % final output: [alpha beta normr]

end % end of function
