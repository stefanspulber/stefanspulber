function [output_matrix, hdr]=a_feature_extraction(datain)
% Extract features from actigraphy recording. Input must be [timestamp(0-1) activity].
% Output order (table format): 
% [ID, age, sex,] ...
% ndays, period, ...
% M10, sdM10, cvM10, M10L, sdM10L, trendM10L, ...
% L5, sdL5, cvL5, L5L, sdL5L, trendL5L, ...
% RA, sdRA, trendRA, ...
% aM10, aM10L, aL5, aL5L, aRA, ...
% alpha_full (4min-24h), alpha_short (4min-4h), alpha_long(4-24h), ...
% IV5, IV30, IV60, IS5, IS30, IS60 

%   Copyright 2020 NorthernLight Diagnostics AB
% 20200907 - updated code to accommodate changes in a_peaktroughs: all
% calculations (including periodogram and all calculations on average
% circadian profile) are performed in a_peaktroughs v3.x; removed progress
% bar (too fast to need one)
% 20200620 - increased resolution for DFA; updated limit long/short
% 20200523 - added 'omitnan' in all mean/std calculations
% 20200504 - updated code to use Ahead-dedicated functions (naming
% convention: 'a_*.m')
% 20200311 - updated code to not display the actual operation performed;
% updated code for DFA calculations (faster processing)
% 20020123 - moved the trimming as first step; added number of days
% analyzed as first output argument (before per)
    
    op=a_speaktroughs(datain(:, 2), 1440, find(datain(:, 1)==0, 1, 'first'), [10 5], 'none'); 
    dint.speaktroughs=op(1:end-1, :);% M10, L5, RA, for each day separately
    dint.apeaktroughs=op(end, :); % M10, L5, RA on average 24h profile

    % trim data to first and last midnight for further calculations
    onoff=[find(datain(:, 1)==0, 1, 'first') find(datain(:, 1)==0, 1, 'last')-1];
    cdata=datain(onoff(1):onoff(2), 2);
    
    % calculate trend in L5L and M10L after removing outliers, i.e. points
    % diverging by more than 3xSD from the average of the rest
    dint.trendM10L=mean(abs(diff(rmoutliers(dint.speaktroughs(:, 2), 'mean'))));
    dint.trendL5L=mean(abs(diff(rmoutliers(dint.speaktroughs(:, 4), 'mean'))));
      
    
%     avg=mean(reshape(cdata, 1440, length(cdata)/1440), 2, 'omitnan');     % 24h-profile
%     dint.apeaktroughs=a_speaktroughs([avg(721:end); avg; avg(1:720)], 1440, 721, [10 5], 'none'); % wrap around 50% to supply enough datapoints
    
    % DFA
    alpha=a_dfa(cdata, round(2.^[2:.25:10.5])); % 4min-24h
    alpha1=a_dfa(cdata, round(2.^[2:.25:7.75])); % 4min-3.5h
    alpha2=a_dfa(cdata, round(2.^[8:.25:10.5])); % 4.25-24h
    dint.dfa=[alpha(1) alpha1(1) alpha2(1)]; % store scaling coeffs

    % IV, IS, etc.
    op=a_idvar(cdata, 1440, [5 30 60]);
    dint.idvar.iv=op.IV;
    dint.idvar.is=op.IS;
    
    % do calculations for output
    output_matrix=[length(cdata)/1440, dint.apeaktroughs(end) ...
        mean(dint.speaktroughs(:, 1), 'omitnan') std(dint.speaktroughs(:, 1), 'omitnan') std(dint.speaktroughs(:, 1), 'omitnan')/mean(dint.speaktroughs(:, 1), 'omitnan') ...
        mean(dint.speaktroughs(:, 2), 'omitnan') std(dint.speaktroughs(:, 2), 'omitnan') dint.trendM10L ...
        mean(dint.speaktroughs(:, 3), 'omitnan') std(dint.speaktroughs(:, 3), 'omitnan') std(dint.speaktroughs(:, 3), 'omitnan')/mean(dint.speaktroughs(:, 3), 'omitnan') ...
        mean(dint.speaktroughs(:, 4), 'omitnan') std(dint.speaktroughs(:, 4), 'omitnan') dint.trendL5L ...
        mean(dint.speaktroughs(:, 5), 'omitnan') std(dint.speaktroughs(:, 5), 'omitnan') ...
        dint.apeaktroughs(1) dint.apeaktroughs(2) dint.apeaktroughs(3) dint.apeaktroughs(4) dint.apeaktroughs(5) ...
        dint.dfa ...
        dint.idvar.iv(:, 2)' ...
        dint.idvar.is(:, 2)']; % mRR];
    
    % prepare header list; remember to update this as you add features to extract
    hdr={'ndays', 'period', ...
        'M10', 'sdM10', 'cvM10', 'M10L', 'sdM10L', 'trendM10L', ...
        'L5', 'sdL5', 'cvL5', 'L5L', 'sdL5L', 'trendL5L', ...
        'RA', 'sdRA', ...
        'aM10', 'aM10L', 'aL5', 'aL5L', 'aRA', ...
        'alpha_full', 'alpha_short', 'alpha_long', ...
        'IV5', 'IV30', 'IV60', 'IS5', 'IS30', 'IS60'};

end % end of function
