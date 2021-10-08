function op=speaktroughs(data, per, first_midnight, opt) 
% SPEAKTROUGHS locates the epochs with highest level of
% activity, and the epochs with lowest level of activity
% within each day recorded. Input data ('data' - required) must be column vector.
% 'per' (required) specifies the period in number of samples/cycle.
% 'first_midnight' (optional) specifes the beginning of the first cycle in case
% of unregistered data. If not specified, it is assumed that the data is
% registered and the first cycle start with the first sample. 'opt' (optional) specifies the width
% of peak and trough in hours, respectively (default [10 5]).
%
% The output is sorted as follows: op=[M10 M10L L5 L5L RA per];
% Each column contains the raw data for the most active, or least active epochs
% and the location of the peak/trough for each day. 
% The peaks and troughs are detected after smoothing with
% sliding Gauss windows, width for each specified by 'opt'. RA - relative amplitude calculated as: 
% RA=(M10-L5)/(M10+L5)
%
%

% Version history: 
% 20201103 - version 4
% removed graphical output otpion - all plotting included in the main app
% 20200906 - version 3.0
% included calculation of periodogram on 3-day running window in standard
% output; changed periodogram estimation method to Lomb-Scargle (NB: the
% periodogram focuses on circadian period for speed); changed 
% peaks&troughs plot (#2) for compatibility with b/w printing; changed output to include 
% all values for average circadian profile as last line in the matrix;
% updated plots #2-3 to display the values calculated on the average
% circadian profile. 
% 20200904 - version 2.8
% changes the display in last plot to include reference values from healthy
% controls (similar to middle plot)
% 20200827 - version 2.75
% changed labels on x axis for first 2 graphs
% 20200713 - version 2.7
% changed order of operations to first reshape then smooth for detection of
% peaks and troughs; general cleaning of code; *set size of figure
% porportional to the number of days recorded*
% 20200603 - version 2.6
% added option to request the handle to figure (optional) to use in the
% report; when requested, the figure is generated without displaying it 
% 20200601 - version 2.5
% corrected a bug that allowed L5 to occur later than M10 within the day
% updated reference values for HC; update help header
% 20200521 - version 2.4
% upon request from potential user: addedd reference values for healthy
% controls: shaded areas in subplot #2
% added line for M10L and L5L for the average circadian profile
% 20200504 - version 2.3
% renamed to use only for Ahead (a_*.m)
% simplified output from structure to matrix [M10 M10L L5 L5L RA] - updated
% code for 'feature_extraction.m'
% 20200423 - version 2.2
% set the limits for hourly averages for M10 and L5 to not exceeed size of
% input data
% 20200416 - version 2.1
% Changed subplot(3, 1, 3) to display chi-square spectrogram
% 20190207 - version 2.0
% Updated the calculation of M10 and L5 to project the indices back to
% original dataset to calculate hourly averages for peaks and troughs.
% Added the calculation of RA as default. Included option for graphic
% display. Restricted input data format to column vector.
% Inception not recorded - version 1.0
 
narginchk(2, 4);
nargoutchk(1, 1); % only one output argument allowed

if nargin==4 
    if length(opt)~=2
%          warning 'Peak/trough width incorrectly specified! Option set to default [10 5].'
         opt=[10 5];
    end
end
    
if nargin<4
    opt=[10 5]; % default for width of peaks and troughs 
    if nargin==2
        first_midnight=1; % set start of recording to reference point
    end
end

if length(data)<per
    error 'Input data shorter than 1 period. You may want to check the input arguments?'
end

if mod(per, 24)~=0
    error 'Specified period does not match 24h!'
end

if size(data, 2)>1
    error 'This version of the software can only handle data in vector format.'
end

% preparatory steps
np=floor(length(data(first_midnight:end, :))/per); % number of integer periods in trimmed data
opt=opt/24*per; % translate into number of samples

% trim data to integer number of periods [first_midnight:last_midnight]
sdm=data(first_midnight+[0:(per*np-1)]);  % trimmed for peaks: midgnight-midnight

sdl=[nan(max(0, per/2-first_midnight), 1);data(max(1, first_midnight-per/2):end)]; % pad with NaNs in the beginning if needed
sdl=sdl(1:per*np); % trimmed for troughs: midday-midday

% reshape to 1 col/day(cycle) and locate M10 and L5 for each day (cycle)
smdata=reshape(sdm, per, np); % M10
[~, M10L]=max(smoothdata(smdata, 'gaussian', 1+opt(1)), [], 'omitnan'); % detect peak location for each day

sldata=reshape(sdl, per, np); % L5
[~, L5L]=min(smoothdata(sldata, 'gaussian', 1+opt(2)), [], 'omitnan'); % detect trough location for each day
L5L=L5L-per/2; % correct for initial shift

% project peak and trough location onto *original* input data
M10lims=M10L+first_midnight-1+[0:np-1]*per+repmat([-opt(1); opt(1)]/2, 1, np);
L5lims=L5L+first_midnight-1 +[0:np-1]*per+repmat([-opt(2); opt(2)]/2, 1, np); 

% the indices may exceed the size of input data --> trim to fit 
M10lims(1)=max(1, M10lims(1)); M10lims(end)=min(length(data), M10lims(end));
L5lims(1)=max(1, L5lims(1)); L5lims(end)=min(length(data), L5lims(end)); 
    
% calculate hourly averages for M10 and L5
for idx=1:np
    M10s(idx)=60*mean(data(max(1, M10lims(1, idx)):min(length(data), M10lims(2, idx)))); % M10, hourly average
    L5s(idx)=60*mean(data(max(1, L5lims(1, idx)):min(length(data), L5lims(2, idx)))); % L5, hourly average
end

% round location to 0-24h system, 4 significant decimals
M10L=round(M10L(:)/per, 4)*24; L5L=round(L5L(:)/per, 4)*24; 

% run calculations for Lomb-Scargle periodogram
if np>1 % several days analyzed
    d=[smdata; [smdata(:, 2:end) smdata(:, 1)]; [smdata(:, 3:end) smdata(:, 1:2)]]; % wrap around to keep the number of days on display
else % average profile analyzed
    d=repmat(smdata, 3, 1);
end

for idx=1:np
    f(:, idx)=sfastlomb(d(:, idx), 1:per*3, 0.01, 10); % LS periodogram for 3 consecutive days
end
[~, imax]=max(f); % locate maxima
imax=imax/30*24; % translate to hours  

% calculation for average profile
% reshape to 1 col/day(cycle) and locate M10 and L5 for each day (cycle)
[pM10, pM10L]=max(smoothdata(mean(smdata, 2, 'omitnan'), 'gaussian', 1+opt(1)), [], 'omitnan'); % detect peak location for average circadian profile
[pL5, pL5L]=min(smoothdata(mean(sldata, 2, 'omitnan'), 'gaussian', 1+opt(2)), [], 'omitnan'); % detect trough location for average circadian profile
pL5L=pL5L-per/2; % correct for initial shift
% round location to 0-24h system, 4 significant decimals
pM10L=round(pM10L(:)/per, 4)*24; pL5L=round(pL5L(:)/per, 4)*24;
        
d=data(first_midnight-1+[1:per*np]); % trim input data for plotting purposes
f=sfastlomb(d, 1:length(d), 0.01, 10); % LS periodogram for entire recording, restricted to 1% of total length
[~, fmax]=max(f); % locate main circadian period
fmax=fmax/(np*10)*24; % translate to hours

% compile output: [M10 M10L L5 L5L RA per] for individual days
op=[M10s(:) M10L(:) L5s(:) L5L(:) (M10s(:)-L5s(:))./(M10s(:)+L5s(:)) imax']; 

% add parameters for average circadian profile
op=[op; [pM10 pM10L pL5 pL5L (pM10-pL5)/(pM10+pL5) fmax]];

end % end of function