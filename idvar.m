function op=idvar(datain, origper, rsample)
% Intradaily variability (IV) and interdaily stability (IS)
% The function calculates indices of sample-to-sample and across-days variability over
% timeseries resampled as specified. 
% Input data must be vector, and 
% must be longer than 3 times the maximum resampling value. The period of 
% original timeseries (required) is specified as 'origper'. The 'rsample'
% vector (required) specifies the bin sizes for resampling. 
% The output is a structure with the following fields: .iv, .is storing, 
% respectively, intradaily variabilityand interdaily stability.  
% 
% Algorithms implemented as described in :
% Goncalves BSB, Adamowicz T, Mazzilli Luzada F, Moreno CR, and Fontenele
% Araujo F (2015) A fresh look at the use of nonparametric analysis in
% actimetry. Sleep Medicine Reviews 20:84-91
% 
% and: 
% Faedda GL, Ohashi K, Hernandez M, McGreenery CE, Grant MC, Baroni A,
% Polcari A, Teicher MH (2016). Actigraph measures discriminate pediatric bipolar
% disorder from attention-deficit/hyperactivity disorder and typically developing 
% controls. Journal of child psychology and psychiatry, and allied
% disciplines doi:10.1111/jcpp.12520


% version history
% 20210212 - otimized calculations -> 10x faster
% 20200620 - reduced output arguments
% 20200504 - general clean-up of the code, including removal of unnecesary
% calculation of TDCOV and MSSD
% inception not recorded

narginchk(3, 3); 

if min(size(rsample))~=1 % resample vector specification
    error 'Specify resample argument as vector!'
end

if 3*max(rsample)>length(datain(:)) % make sure resampling is reasonable
    error 'Resample vector incorrectly specified: allow a minimum of 3 bins after resampling!'
end 

% check whether all values in rsample vector are divisors of the original
% period and trim as apropriate
t=find(mod(repmat(origper, numel(rsample), 1), rsample(:))==0); 
if numel(t)<numel(rsample) % not all resampling bins yield integer periods
    rsample=rsample(t); % trim resample bin specification
end 

datain=datain(:); % make input data to vector 
datain=datain(1:origper*(floor(length(datain)/origper))); % trim data to fit an integer number of periods periods

% preallocate output structure
op=struct('iv', [], 'is', [], 'ivd', []); 

for ridx=1:numel(rsample) 
    
    % resample input data
    d=sum(reshape(datain, rsample(ridx), []))';
    
    % reconstitute days
    dd=reshape(d, origper/rsample(ridx), []);
    
    % compute average profile
    dp=mean(dd, 2);
    
    % calculate output
    op.iv(ridx, :)=[rsample(ridx) mean((diff(d)).^2)/mean((d-mean(d)).^2)];
    op.is(ridx, :)=[rsample(ridx) mean((dp-mean(dp)).^2)/mean((d-mean(d)).^2)];
    op.ivd(ridx, :)=mean((diff(dd)).^2)./mean((dd-mean(dd)).^2);
end % next resampling bin

end % end of function
